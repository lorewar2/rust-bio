use bio::alignment::pairwise::Scoring;
use bio::alignment::{poa::*, TextSlice};
use bio::data_structures::qgram_index::Match;
use petgraph::dot::{Dot, Config};
fn main() {
    let test = vec![
        "GCCCTC",
        "AAAT"
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
    run(test);
}

fn run(seqvec : Vec<&str>) {
    ////////////////////////
    //normal poa alignment//
    ////////////////////////
    
    let scoring = Scoring::new(-2, -2, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
    let mut i = 0;
    let mut aligner = Aligner::new(scoring, seqvec[0].as_bytes());
    for seq in &seqvec{
        if i != 0 {
            aligner.global(seq.as_bytes()).add_to_graph(i);
        }
        i += 1;
    }
    let normal_consensus = aligner.poa.consensus();
    //get scores of sequences compared to normal consensus
    let mut normal_score = 0;
    for seq in &seqvec{
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(normal_consensus.len(), seq.len(), -2, -2, &score);
        let alignment = aligner.global(&normal_consensus, seq.as_bytes());
        normal_score += alignment.score;
    }
    ////////////////////////////
    //compressed poa alignment//
    ////////////////////////////
    
    //make homopolymer compression vector
    let mut homopolymer_vec: Vec<HomopolymerSequence> = vec![];
    for seq in &seqvec{
        homopolymer_vec.push(HomopolymerSequence::new(seq.as_bytes()));
    }

    let scoring = Scoring::new(-2, -2, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
    let mut i = 0;
    let mut aligner = Aligner::new(scoring, &homopolymer_vec[0].bases);
    for homopolymer_seq in &homopolymer_vec{
        if i != 0 {
            aligner.global(&homopolymer_seq.bases).add_to_graph(i);
        }
        i += 1;
    }
    let homopolymer_consensus = aligner.poa.consensus();

    //use homopolymer compressions sequences to make expanded consensus
    let mut homopolymer_score = 0;

    for homopolymer_seq in &homopolymer_vec{
        let mut consensus_index = 0;
        let mut homopolymer_index = 0;
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(homopolymer_consensus.len(), homopolymer_seq.bases.len(), -2, -2, &score);
        let alignment = aligner.global(&homopolymer_consensus, &homopolymer_seq.bases);
        consensus_index = alignment.xstart;
        homopolymer_index = alignment.ystart;
        for op in alignment.operations {
            match op {
                Match => {
                    if (homopolymer_consensus[consensus_index] == homopolymer_seq.bases[homopolymer_index]){
                        println!("ssss");
                    }
                    println!("{} Match {}", homopolymer_consensus[consensus_index], homopolymer_seq.bases[homopolymer_index]);
                    //consensus_index += 1;
                    //homopolymer_index += 1;
                },
                Subst => {
                    println!("{} UnMatch {}", homopolymer_consensus[consensus_index], homopolymer_seq.bases[homopolymer_index]);
                    //consensus_index += 1;
                    //homopolymer_index += 1;
                },
                Del => {
                    println!("Del {}", homopolymer_consensus[consensus_index]);
                    //consensus_index += 1;
                },
                Ins => {
                    println!("Ins {}", homopolymer_seq.bases[homopolymer_index]);
                    //homopolymer_index += 1;
                },
            }
        }
        homopolymer_score += alignment.score;
    }
    //print the results
    print!("Normal consensus:\t");
    for i in &normal_consensus{
        print!("{}", *i as char);
    }
    println!("");
    println!("Normal consensus score:\t{}", normal_score);
    print!("Homopolymer consensus:\t");
    for i in &homopolymer_consensus{
        print!("{}", *i as char);
    }
    println!("");
    println!("Homopolymer consensus score: {}", homopolymer_score);
}
pub struct HomopolymerSequence {
    pub bases: Vec<u8>,
    pub frequencies: Vec<u8>,
}

impl HomopolymerSequence {
    fn new(query: TextSlice) -> Self{
        let mut temp_bases = vec![];
        let mut temp_frequencies = vec![];
        let mut prev_base = 0;
        for &base in query{
            if prev_base == base{
                if let Some(last) = temp_frequencies.last_mut() {
                    *last = *last + 1;
                }
            }
            else {
                temp_bases.push(base);
                temp_frequencies.push(1);
            }
            prev_base = base;
        }
        HomopolymerSequence {
            bases: temp_bases,
            frequencies: temp_frequencies,
        }
    }
    fn read_query(&mut self, query: TextSlice) -> &mut Self {
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

