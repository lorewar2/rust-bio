use bio::alignment::pairwise::Scoring;
use bio::alignment::{poa::*, TextSlice};
use petgraph::dot::{Dot, Config};
fn main() {
    let test = vec![
        "GCC",
        "CTGCC",
    ];
    let example4 = vec![
        "TGTACNTGTTTGTGAGGCTA",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
        "AGTTCCTGCTGCGTTTGCTGGACTTATGTACTTGTTTGTGAGGCAA",
        "AAGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGGTTTGTGNAGGCAA",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTNAGGCAA",
        "AGTTCCTGCTGCGTTTGCT",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTT",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
        "AGTTNCTGNTGNGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
        "GTACNTGTTTGTGAGGCTA",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
        "AGTTCCTGCTGCTTTTGCTGGACTGATGTACTTGATTGTGAGGCAA",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
        "AGTTCCTGCTGCGCTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGCGGCAA",
        "AGTCCTGCGCGTTTGCGGACGGATGTACTTGTTGTGAGGCAA",
        "GCAA",
        "GGCAA",
        "CTGATGTACTTGTTGTGAGGGCAA",
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA",
        "GTTCTGCCTGCGTTTGCTGAACTGATGTACTTGTTAGTAAGCAA",
        "CGTTACTGCGGGGTTTGCTGGACTCATGACTTTGTTNGTAGGCAA",
    ];
    let seqvec = test;
    
    let scoring = Scoring::new(-2, -2, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
    let mut i = 0;
    let mut aligner = Aligner::new(scoring, seqvec[0].as_bytes());
    for seq in seqvec{
        if i != 0 {
            aligner.global(seq.as_bytes()).add_to_graph(i);
            println!("{}", seq);
        }
        i += 1;
        
    }
    aligner.poa.consensus();
    let mut index: u8 = 0;
    for node in aligner.poa.node_seq_tracker {
        println!("Node index: {} ", index);
        for j in node {
            println!("seq {} " , j);
        }
        index += 1;
    }
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

