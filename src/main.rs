use bio::alignment::pairwise::Scoring;
use bio::alignment::{poa::*, TextSlice};
use petgraph::dot::{Dot, Config};
use std::{
    fs::File,
    io::{prelude::*, BufReader},
    path::Path,
};


const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = -2;
const MATCH: i32 = 1;
const MISMATCH: i32 = -1;

fn get_fasta_sequences(filename: impl AsRef<Path>) -> Vec<String> {
    let mut tempvec: Vec<String> = vec![];
    let mut seqvec: Vec<String> = vec![];
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    let lines: Vec<String> = buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();
    //get the indices of >
    let mut indices: Vec<usize> = vec![];
    let mut index = 0;
    for line in &lines{
        if line.as_bytes()[0] as char == '>'{
            //println!("{:?}", line);
            indices.push(index);
        }
        index += 1;
    }
    //join the lines between >s and remove > lines
    let mut prev_index: usize = 0;
    for index in indices {
        if prev_index != 0 && index != 0{
            tempvec.push(lines[(prev_index + 1)..index].join(""));
        }
        prev_index = index;
    }
    //reverse complement every other line
    let mut index = 0;
    for seq in &tempvec{
        if index % 2 != 0 {
            let mut tempseq: Vec<char> = vec![];
            let iterator = seq.chars().rev().into_iter();
            for char in iterator{
                tempseq.push(match char {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => ' ',
                });
            }
            seqvec.push(tempseq.iter().cloned().collect::<String>());
        }
        else {
            seqvec.push((*seq.clone()).to_string());
        }
        index += 1;
    }
    //remove the last line
    seqvec.pop();

    seqvec
}

fn main() {
    let seqvec = get_fasta_sequences("./data/65874.fasta");

    let test = vec![
        "AAAATGCC".to_string(),
        "AAAAT".to_string(),
        "AAAAT".to_string(),
    ];
    let example4 = vec![
        "TGTACNTGTTTGTGAGGCTA".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTTATGTACTTGTTTGTGAGGCAA".to_string(),
        "AAGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGGTTTGTGNAGGCAA".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTNAGGCAA".to_string(),
        "AGTTCCTGCTGCGTTTGCT".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTT".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA".to_string(),
        "AGTTNCTGNTGNGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA".to_string(),
        "GTACNTGTTTGTGAGGCTA".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA".to_string(),
        "AGTTCCTGCTGCTTTTGCTGGACTGATGTACTTGATTGTGAGGCAA".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA".to_string(),
        "AGTTCCTGCTGCGCTTGCTGGACTGATGTACTTGTTTGTGAGGCAA".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGCGGCAA".to_string(),
        "AGTCCTGCGCGTTTGCGGACGGATGTACTTGTTGTGAGGCAA".to_string(),
        "GCAA".to_string(),
        "GGCAA".to_string(),
        "CTGATGTACTTGTTGTGAGGGCAA".to_string(),
        "AGTTCCTGCTGCGTTTGCTGGACTGATGTACTTGTTTGTGAGGCAA".to_string(),
        "GTTCTGCCTGCGTTTGCTGAACTGATGTACTTGTTAGTAAGCAA".to_string(),
        "CGTTACTGCGGGGTTTGCTGGACTCATGACTTTGTTNGTAGGCAA".to_string(),
    ];
    run(seqvec);
}

fn get_consensus_score(seqvec : &Vec<String>, consensus: &Vec<u8>) -> i32{
    let mut consensus_score = 0;
    for seq in seqvec{
        let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
        let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(consensus.len(), seq.len(), GAP_OPEN, GAP_EXTEND, &score);
        let alignment = aligner.global(&consensus, seq.as_bytes());
        consensus_score += alignment.score;
    }
    consensus_score
}

fn get_expanded_consensus(homopolymer_vec: Vec<HomopolymerSequence>, homopolymer_consensus: &Vec<u8>) -> (Vec<u8>, Vec<u32>, i32)  {
    //use homopolymer compressions sequences to make expanded consensus //make function
    let mut homopolymer_score = 0;
    let mut homopolymer_consensus_freq: Vec<u32> = vec![0; homopolymer_consensus.len()];
    for homopolymer_seq in &homopolymer_vec{
        let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
        let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(homopolymer_consensus.len(), homopolymer_seq.bases.len(), GAP_OPEN, GAP_EXTEND, &score);
        let alignment = aligner.global(&homopolymer_consensus, &homopolymer_seq.bases);
        let mut consensus_index = alignment.xstart;
        let mut homopolymer_index = alignment.ystart;
        //println!("start index consensus {}", consensus_index);
        //println!("start index sequence {}", homopolymer_index);
        for op in alignment.operations {
            match op {
                bio::alignment::AlignmentOperation::Match => {
                    //println!("{} Match {}", homopolymer_consensus[consensus_index], homopolymer_seq.bases[homopolymer_index]);
                    homopolymer_consensus_freq[consensus_index] += homopolymer_seq.frequencies[homopolymer_index];
                    homopolymer_index += 1;
                    consensus_index += 1;
                },
                bio::alignment::AlignmentOperation::Subst => {
                    //println!("{} MisMatch {}", homopolymer_consensus[consensus_index], homopolymer_seq.bases[homopolymer_index]);
                    homopolymer_index += 1;
                    consensus_index += 1;
                },
                bio::alignment::AlignmentOperation::Del => {
                    //println!("Del {}", homopolymer_seq.bases[homopolymer_index]);
                    homopolymer_index += 1;
                },
                bio::alignment::AlignmentOperation::Ins => {
                    //println!("Ins {}", homopolymer_consensus[consensus_index]);
                    consensus_index += 1;
                    
                },
                _ => {},
            }
        }
        homopolymer_score += alignment.score;
    }

    //make the expanded consensus using the frequencies
    //++ average ++
    let mut expanded_consensus: Vec<u8> = vec![];
    for i in 0..homopolymer_consensus.len(){
        expanded_consensus.push(homopolymer_consensus[i]);
        for _ in 0..(homopolymer_consensus_freq[i] / homopolymer_vec.len() as u32) as i32 - 1 {
            expanded_consensus.push(homopolymer_consensus[i]);
        }
    }
    (expanded_consensus, homopolymer_consensus_freq, homopolymer_score)
}

fn run(seqvec: Vec<String>) {
    ////////////////////////
    //normal poa alignment//
    ////////////////////////
    
    let scoring = Scoring::new(GAP_OPEN, GAP_EXTEND, |a: u8, b: u8| if a == b { MATCH } else { MISMATCH });
    let mut seqnum: u8 = 0;
    let mut aligner = Aligner::new(scoring, seqvec[0].as_bytes());
    for seq in &seqvec{
        if seqnum != 0 {
            aligner.global(seq.as_bytes()).add_to_graph(seqnum);
        }
        seqnum += 1;
    }
    let normal_consensus = aligner.poa.consensus();

    //get scores of sequences compared to normal consensus 
    let normal_score = get_consensus_score(&seqvec, &normal_consensus);

    ////////////////////////////
    //compressed poa alignment//
    ////////////////////////////
    
    //make homopolymer compression vector
    let mut homopolymer_vec: Vec<HomopolymerSequence> = vec![];
    for seq in &seqvec{
        homopolymer_vec.push(HomopolymerSequence::new(seq.as_bytes()));
    }

    let scoring = Scoring::new(GAP_OPEN, GAP_EXTEND, |a: u8, b: u8| if a == b { MATCH } else { MISMATCH });
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
    let (expanded_consensus, homopolymer_consensus_freq, homopolymer_score) =  get_expanded_consensus(homopolymer_vec, &homopolymer_consensus);
    //get the scores of expanded consensus compared to sequences
    let expanded_score = get_consensus_score(&seqvec, &expanded_consensus);

    //print the results
    print!("Normal consensus:\t\t");
    for i in &normal_consensus{
        print!("{}", *i as char);
    }
    print!("\nNormal consensus score:\t\t{}", normal_score);
    print!("\nHomopolymer consensus:\t\t");
    for i in &homopolymer_consensus{
        print!("{}", *i as char);
    }
    print!("\nHomopolymer score: \t\t{}", homopolymer_score);
    print!("\nHomopolymer consensus freq:\t");
    print!("{:?}", homopolymer_consensus_freq);
    print!("\nExpanded consensus:\t\t");
    for i in &expanded_consensus{
        print!("{}", *i as char);
    }
    print!("\nExpanded consensus score:\t{}", expanded_score);
    println!("");

}
pub struct HomopolymerSequence {
    pub bases: Vec<u8>,
    pub frequencies: Vec<u32>,
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
}

