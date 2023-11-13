use bio;
use bio::alignment::{Alignment, AlignmentOperation};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use debruijn::dna_string::DnaString;
use bio::alignment::pairwise::banded;
use bio::alignment::sparse::HashMapFx;

pub type Scoring = bio::alignment::pairwise::Scoring<bio::alignment::pairwise::MatchParams>;

/// Alignment against a reference
/// This struct pack the index of the reference and the alignment together.
#[derive(Debug, Clone)]
pub struct AlignmentPacket {
    pub ref_idx: usize, // index of reference
    pub alignment: Alignment,
}

impl AlignmentPacket {
    pub fn new(ref_idx: usize, alignment: Alignment) -> Self {
        AlignmentPacket { ref_idx, alignment }
    }

    /// Compute the edit distance implied by the alignment.
    /// The clipped regions are not considered in the computation
    pub fn edit_distance(&self) -> usize {
        let mut d = 0;
        for &op in &self.alignment.operations {
            d += match op {
                AlignmentOperation::Match => 0,
                AlignmentOperation::Xclip(_) => 0,
                AlignmentOperation::Yclip(_) => 0,
                _ => 1,
            };
        }
        d
    }
}

pub struct Aligner<'a> {
    pub refs: Vec<DnaString>,
    pub kmers_hash: HashMapFx<&'a [u8], Vec<(usize, u32)>>,
    scoring: Scoring,
    pub k: usize,
    pub w: usize,
}
impl<'a> Aligner<'a> {
    /// Creates a new AlignHelper instance given the set of reference sequences,
    /// the scoring parameters, and the banding parameters (kmer length and window size)
    ///
    /// # Arguments
    ///
    /// * `refs` - vector of reference sequences (stored as String)
    /// * `scoring` - Scoring struct
    /// * `k` - kmer length for constructing the band (see bio::alignment::pairwise::banded)
    /// * `w` - window size for constructing the band (see bio::alignment::pairwise::banded)
    ///
    pub fn new(refs: &'a [&str], scoring: Scoring, k: usize, w: usize) -> Self {
        // Make sure that the scoring implies that the alignment is local in the reference y
        assert!(scoring.yclip_prefix == 0);
        assert!(scoring.yclip_suffix == 0);
        let mut kmers_hash: HashMapFx<&'a [u8], Vec<(usize, u32)>> = HashMapFx::default();
        let mut dna_strings = Vec::new();
        for (ref_idx, ref_seq) in refs.iter().enumerate() {
            let dna_string = DnaString::from_dna_string(ref_seq);
            dna_strings.push(dna_string);
            let seq = ref_seq.as_bytes();
            for i in 0..(seq.len() + 1).saturating_sub(k) {
                kmers_hash
                    .entry(&seq[i..i + k])
                    .or_insert_with(Vec::new)
                    .push((ref_idx, i as u32));
            }
        }
        Aligner {
            refs: dna_strings,
            kmers_hash,
            scoring,
            k,
            w,
        }
    }
    pub fn set_scoring(&mut self, scoring: Scoring) {
        // Make sure that the scoring implies that the alignment is local in the reference y
        assert!(scoring.yclip_prefix == 0);
        assert!(scoring.yclip_suffix == 0);
        self.scoring = scoring;
    }
    /// Align the read with at least a minimum SW score. The read is aligned with the
    /// sequences in refs and the first acceptable alignment is returned. We use the
    /// banded aligner from rust-bio (see bio::alignment::pairwise::banded). The function
    /// returns Some(AlignmentPacket) if we could find a good enough alignment, otherwise None.
    /// First we find all the kmer matches between the read and all the refs using the
    /// precomputed 'kmers_hash'. We sort the matches based on the number of kmer
    /// matches per reference and then call the aligner.
    ///
    /// # Arguments
    ///
    /// * `read` - Read object.
    /// * `min_align_score` - Minimum SW score for accepting an alignment
    ///
    pub fn find_read_matches(
        &self,
        read: &str,
        min_align_score: i32,
    ) -> Option<AlignmentPacket> {
        let mut all_kmer_matches: Vec<(i32, usize, usize, usize)> = Vec::new();
        let read_seq = read.to_string();
        let read_seq = read_seq.as_bytes();
        let mut ref_kmer_counts = vec![0i32; self.refs.len()];
        for i in 0..(read_seq.len() + 1).saturating_sub(self.k) {
            let slc = &read_seq[i..i + self.k];
            if let Some(matches) = self.kmers_hash.get(slc) {
                for &ref_pos in matches {
                    all_kmer_matches.push((0, ref_pos.0, i, ref_pos.1 as usize));
                    ref_kmer_counts[ref_pos.0] += 1;
                }
            }
        }
        for i in 0..all_kmer_matches.len() {
            all_kmer_matches[i].0 = -ref_kmer_counts[all_kmer_matches[i].1];
        }
        all_kmer_matches.sort_unstable();
        let mut aligner = banded::Aligner::with_scoring(self.scoring, self.k, self.w);
        for (ref_idx, match_group) in &all_kmer_matches.into_iter().group_by(|x| x.1) {
            let matches = match_group
                .into_iter()
                .map(|x| (x.2 as u32, x.3 as u32))
                .collect();
            let ref_seq = &self.refs[ref_idx].to_string();
            let alignment = aligner.custom_with_expanded_matches(
                read_seq,
                ref_seq.as_bytes(),
                matches,
                Some(1),
                true,
            );
            if alignment.score > min_align_score {
                let align_pack = AlignmentPacket::new(ref_idx, alignment);
                // assert!(read.validate_alignment(&align_pack));
                return Some(align_pack);
            }
        }
        None
    }
}

#[test]
fn test_find_read_matches(){
    let refs = ["CCCCACGTCCACGTGGGGGGA", "ACGTACGTACGTGGGGG"];

    let scoring = Scoring::from_scores(-4, -2, 2, -2).xclip(-5).yclip(0);
    let mut align_helper = Aligner::new(&refs, scoring, 8, 8);

    let quals = vec![30; 8];
    let reads = vec![
        "CCCACCTACGTGGG",
        "TTTTTTTT",
    ];
    let matches = align_helper.find_read_matches(reads[0], 15);
    // assert_eq!(matches.unwrap().alignment.score, 16);
    dbg!(matches);

}

#[test]
fn test_find_read_matches_1() {
    let refs = ["CCCCACGTACGTGGGGGGA", "ACGTACGTACGTGGGGG"];

    let scoring = Scoring::from_scores(-4, -2, 2, -2).xclip(-5).yclip(0);
    let mut align_helper = Aligner::new(&refs, scoring, 3, 8);

    let quals = vec![30; 8];
    let reads = vec![
        "ACGTACCC",
        "TTTTTTTT",
        "ACGTACGT",
        "ACGAACGA"
    ];

    for read in &reads {
        let matches = align_helper.find_read_matches(read, 100);
        assert!(matches.is_none());
    }

    let matches = align_helper.find_read_matches(&reads[2], 15);
    assert_eq!(matches.unwrap().alignment.score, 16);

    align_helper.set_scoring(Scoring::from_scores(-4, -1, 2, -2).xclip(-5).yclip(0));
    let cigars = vec!["6=2X", "", "8=", "3=1X3=1X", "6=2X", "2X6="];
    let dists = vec![2, 0, 0, 2, 2, 2, 0];
    for i in 0..reads.len() {
        let matches = align_helper.find_read_matches(&reads[i], 0);
        if i == 1 {
            assert!(matches.is_none());
        } else {
            let unwrapped = matches.unwrap();
            assert_eq!(unwrapped.alignment.cigar(false), cigars[i]);
            assert_eq!(unwrapped.edit_distance(), dists[i]);
        }
    }

    align_helper.set_scoring(Scoring::from_scores(-4, -2, 2, -2).xclip(-3).yclip(0));
    let cigars = vec!["6=2S", "", "8=", "3=1X3=1X", "6=2S", "2S6=", "1="];
    let dists = vec![0, 0, 0, 2, 0, 0, 0];
    for i in 0..reads.len() {
        let matches = align_helper.find_read_matches(&reads[i], 0);
        if i == 1 {
            assert!(matches.is_none());
        } else {
            let unwrapped = matches.unwrap();
            assert_eq!(unwrapped.alignment.cigar(false), cigars[i]);
            assert_eq!(unwrapped.edit_distance(), dists[i]);
        }
    }
}