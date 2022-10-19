use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::*;

fn main() {

    let x = b"AAAAAAA";
    let y = b"AABBBAA";
    let z = b"AABCBAA";

    let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
    let mut aligner = Aligner::new(scoring, x);
    // z differs from x in 3 locations
    assert_eq!(aligner.global(z).alignment().score, 1);
    aligner.global(y).add_to_graph();
    // z differs from x and y's partial order alignment by 1 base
    assert_eq!(aligner.global(z).alignment().score, 5);
}