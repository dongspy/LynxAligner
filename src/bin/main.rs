use bio::alignment::pairwise::Scoring;
use aligners::aligner::Aligner;

fn main() {
    let refs = ["CCCCACGTCCACGTGGGGGGA", "ACGTACGTACGTGGGGG"];

    let scoring = Scoring::from_scores(-4, -2, 2, -2).xclip(-5).yclip(0);
    let aligner = Aligner::new(&refs, scoring, 8, 8);
    let reads = vec![
        "CCCACCTACGTGGG",
    ];
    let matches = aligner.find_read_matches(reads[0], 15);
    dbg!(matches);
}
