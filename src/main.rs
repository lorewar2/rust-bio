use bio::alignment::pairwise::Scoring;
use bio::alignment::{bandedpoa::*, TextSlice}; //use bio::alignment::{poa::*, TextSlice};
use std::{
    fs::File,
    fs::OpenOptions,
    io::{prelude::*, BufReader},
    path::Path,
};
use chrono;
use rand::{Rng,SeedableRng};
use rand::rngs::StdRng;

const GAP_OPEN: i32 = -4;
const GAP_EXTEND: i32 = -2;
const MATCH: i32 = 2;
const MISMATCH: i32 = -4;
const FILENAME: &str = "./data/65874.fasta";
const SEED: u64 = 1337;

fn main() {
    //let seqvec = get_fasta_sequences_from_file(FILENAME);
    let seqvec = get_random_sequences_from_generator(100, 10);
    println!("generated string: {}", seqvec[0]);
    run(seqvec);
    //to get consensus score from file (abPOA test)
    //let abpoa_consensus = get_consensus_from_file("./data/cons.fa");
    //let abpoa_consensus_score = get_consensus_score(&seqvec, &abpoa_consensus);
    //println!("abpoa score: {:?}", abpoa_consensus_score);
    //write_filtered_data_fasta_file("./results/filtered_data.fa", &seqvec);
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
        println!("Sequence {} processed", seqnum);
    }
    let normal_consensus;
    (normal_consensus, _) = aligner.bandedpoa.consensus(); //just poa

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
        println!("Sequence {} processed", seqnum);
    }
    let homopolymer_consensus;
    (homopolymer_consensus, _) = aligner.bandedpoa.consensus(); //poa

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

    //write results to file
    write_scores_result_file("./results/results.txt", normal_score, homopolymer_score, expanded_score);
    write_consensus_fasta_file("./results/consensus.fa", &normal_consensus, &homopolymer_consensus, &expanded_consensus);

}

fn get_random_sequences_from_generator(sequence_length: i32, num_of_sequences: i32) -> Vec<String> {
    let mut rng = StdRng::seed_from_u64(SEED);
    //vector to save all the sequences 
    let mut randomvec: Vec<String> = vec![];
    //generate the first sequence of random bases of length sequence_length
    let mut firstseq: Vec<char> = vec![];
    for _ in 0..sequence_length{
        firstseq.push(match rng.gen_range(0..4) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'X'
        });
    }
    randomvec.push(firstseq.iter().collect::<String>());
    //loop for 10 
    for _ in 0..num_of_sequences{
        //clone the sequence
        let mut mutseq = firstseq.clone();
        //mutate the all the bases with 0.05 chance
        for i in 0..mutseq.len() {
            match rng.gen_range(0..20){
                0 => {
                    mutseq[i] = match rng.gen_range(0..4){
                        0 => 'A',
                        1 => 'C',
                        2 => 'G',
                        3 => 'T',
                        _ => 'X'
                    }
                },
                _ => {}
            }
        }
        //put indels at location with chance 0.1 
        for i in 0..mutseq.len() {
            let mean_value: f64 = 1.5;
            //get length of the indel geometric distributed mean value 1.5
            let indel_length: usize  = ((1.0 - rng.gen::<f64>()).ln() / (1.00 - (1.00 / mean_value) as f64).ln()).ceil() as usize;
            match rng.gen_range(0..20){
                //insertion of elements
                0 => {
                    if i + indel_length < mutseq.len(){
                        for _ in 0..indel_length{
                            mutseq.insert(i + 1, mutseq[i]);
                        }
                    }
                },
                //deletion of elements
                1 => {
                    if i + indel_length < mutseq.len(){
                        for _ in 0..indel_length{
                            mutseq.remove(i);
                        }
                    }
                }
                _ => {}
            }
        }
        println!("{:?}", mutseq.iter().collect::<String>());
        //insert to vector
        randomvec.push(mutseq.iter().collect::<String>());
    }
    randomvec

}

fn write_scores_result_file(filename: impl AsRef<Path>, normal_score: i32, homopolymer_score: i32, expanded_score: i32) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
            "{:?} \nFILE: {}\nNormal score:\t\t\t{}\nHomopolymer score:\t\t{}\nExpanded score:\t\t\t{}",
            chrono::offset::Local::now(), FILENAME, normal_score, homopolymer_score, expanded_score)
            .expect("result file cannot be written");
}

fn write_consensus_fasta_file(filename: impl AsRef<Path>, normal_consensus: &Vec<u8>, homopolymer_consensus: &Vec<u8>, expanded_consensus: &Vec<u8>) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n>Normal consensus:\n{}\n>Homopolymer consensus:\n{}\n>Expanded consensus:\n{}",
        chrono::offset::Local::now(), FILENAME, std::str::from_utf8(normal_consensus).unwrap(), std::str::from_utf8(homopolymer_consensus).unwrap(), std::str::from_utf8(expanded_consensus).unwrap())
        .expect("result file cannot be written");
}

fn write_filtered_data_fasta_file(filename: impl AsRef<Path>, seqvec: &Vec<String>) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    let mut index = 1;
    for seq in seqvec {
        writeln!(file,
            ">seq {}\n{}",
            index, seq)
            .expect("result file cannot be written");
        index += 1;
    }
}

fn get_consensus_from_file(filename: impl AsRef<Path>) -> Vec<u8> {
    let mut consensus: Vec<u8> = vec![];
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    let lines: Vec<String> = buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();
    consensus = (*lines[1].as_bytes()).to_vec();
    consensus
}

fn get_fasta_sequences_from_file(filename: impl AsRef<Path>) -> Vec<String> {
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
    //get rid of the last incomplete reading
    seqvec.pop();
    //sort the vector by size
    seqvec.sort_by_key(|seq| seq.len());
    //drop the sequences which are > 1.8x median size
    let mut drop_index = seqvec.len();
    let median_size: f32 = seqvec[(seqvec.len() / 2) - 1].len() as f32;
    for index in (seqvec.len() / 2)..(seqvec.len() - 1) {
        if seqvec[index].len() as f32 > (median_size * 1.8) {
            drop_index = index;
            break;
        }
    }
    for _ in drop_index..seqvec.len() {
        seqvec.pop();
    }
    // rearrange the seq vector median first and rest according median size difference
    seqvec.sort_by(|a, b| ((a.len() as f32 - median_size).abs()).partial_cmp(&(b.len() as f32 - median_size).abs()).unwrap());
    seqvec
}

fn get_pair_score(seq1 : String, seq2: String){
    let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
    let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(seq1.len(), seq2.len(), GAP_OPEN, GAP_EXTEND, &score);
    let alignment = aligner.global(seq1.as_bytes(), seq2.as_bytes());
    println!("{}", alignment.score);
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
        for _ in 0..(homopolymer_consensus_freq[i] as f32 / homopolymer_vec.len() as f32).round() as i32 - 1 {
            expanded_consensus.push(homopolymer_consensus[i]);
        }
    }
    (expanded_consensus, homopolymer_consensus_freq, homopolymer_score)
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

