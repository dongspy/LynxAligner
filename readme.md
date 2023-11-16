# LynxAligner

 > Notice: This project is still under active development and not guaranteed to have a stable API.

LynxAligner, inspired by the agility and precision of a lynx, is a sequence alignment software developed using Rust and the seed-extend algorithm, designed for efficient and accurate bioinformatics data analysis.

The seed-extend model is a popular approach in bioinformatics for sequence alignment, especially useful for handling large genomic datasets efficiently. This model typically works in two stages:

* Seed: Identifies short, exact matches (seeds) between two sequences or within a database. These seeds are potential starting points for longer alignments.

* Extend: Expands these seeds in both directions, aligning additional nucleotides until the alignment score falls below a certain threshold.


This approach is effective in balancing the speed (by quickly identifying potential alignment regions through seeding) and accuracy (by extending and scoring these regions) of the alignment process.

## Install

```toml
# Cargo.toml
[dependencies]
aligners = {git="https://github.com/dongspy/aligners.git"}
```

## Usage

```rust
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

```
