use bio::alignment::pairwise::Scoring;
use bio::alignment::{poa::*, TextSlice};

fn main() {

    let x = b"ABCDEFG";
    let y = b"AABBBAA";
    let z = b"AABCBAA";

    //turn the sequences into homopolymer compressed
    let mut x_hps = HomopolymerSequence::new();
    let mut y_hps = HomopolymerSequence::new();
    let mut z_hps = HomopolymerSequence::new();
    x_hps.read_query(&x.to_vec());
    y_hps.read_query(&y.to_vec());
    z_hps.read_query(&z.to_vec());

    //scoring definition
    let scoring = Scoring::new(-4, -2, |a: u8, b: u8| if a == b { 2i32 } else { -4i32 });
    let scoring_hps = Scoring::new(-4, -2, |a: u8, b: u8| if a == b { 2i32 } else { -4i32 });
    //normal aligner (not compressed)
    let mut aligner = Aligner::new(scoring, x);

    aligner.global(y).add_to_graph();
    aligner.global(z).alignment().score;

    //hps aligner (compressed)
    let mut aligner_hps = Aligner::new(scoring_hps, &x_hps.bases);
    
    aligner_hps.global(&y_hps.bases).add_to_graph();
    aligner_hps.global(&z_hps.bases).add_to_graph();
    let hps_graph = aligner_hps.graph();

    println!("{}",hps_graph.node_count());

}

pub struct HomopolymerSequence {
    pub bases: Vec<u8>,
    pub frequencies: Vec<u8>,
}

impl HomopolymerSequence {
    fn new() -> Self{
        HomopolymerSequence {
            bases: Vec::new(),
            frequencies: Vec::new(),
        }
    }
    fn read_query(&mut self, query: &Vec<u8>) -> &mut Self {
        let mut prev_base = 0;
        for &base in query{
            if prev_base == base{
                if let Some(last) = self.frequencies.last_mut() {
                    *last = *last + 1;
                }
            }
            else {
                self.bases.push(base);
                self.frequencies.push(1);
            }
            prev_base = base;
        }
        self
    }
    fn print_sequence(&self){
        println!("base sequence{:?}", self.bases);
        println!("frequency sequence{:?}", self.frequencies);
    }
}

