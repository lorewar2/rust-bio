use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::*;

fn main() {

    let x = b"ABCDEFG";
    let y = b"AABBBAA";
    let z = b"AABCBAA";

    let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
    let mut aligner = Aligner::new(scoring, x);
    
    //aligner.global(y).add_to_graph();
    aligner.global(z).alignment().score;

}