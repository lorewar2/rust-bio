use bio::alignment::pairwise::Scoring;
use bio::alignment::{poa::*, TextSlice}; //bandedpoa/ poa
use std::{
    fs::File,
    fs::OpenOptions,
    io::{prelude::*, BufReader},
    path::Path,
};
use chrono;
use rand::{Rng,SeedableRng};
use rand::rngs::StdRng;
use petgraph::dot::{Dot, Config};

const GAP_OPEN: i32 = -4;
const GAP_EXTEND: i32 = -2;
const MATCH: i32 = 2;
const MISMATCH: i32 = -4;
const FILENAME: &str = "./data/65874.fasta";
const SEED: u64 = 3;
const CONSENSUS_METHOD: u8 = 1; //0==average 1==median

fn main() {
    let seqvec = get_fasta_sequences_from_file(FILENAME);
    //let seqvec = get_random_sequences_from_generator(1000, 10);
    //println!("generated string: {}", seqvec[0]);
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
            aligner.global(seq.as_bytes()).add_to_graph();
        }
        seqnum += 1;
        println!("Sequence {} processed", seqnum);
    }
    let normal_consensus;
    (normal_consensus, _) = aligner.poa.consensus(); //just poa

    //get scores of sequences compared to normal consensus 
    let normal_score = get_consensus_score(&seqvec, &normal_consensus);
    //let graph = aligner.poa.graph;
    //println!("normal graph \n {:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]));
    ////////////////////////////
    //compressed poa alignment//
    ////////////////////////////
    
    //make homopolymer compression vector
    let mut homopolymer_vec: Vec<HomopolymerSequence> = vec![];
    for seq in &seqvec{
        homopolymer_vec.push(HomopolymerSequence::new(seq.as_bytes()));
    }

    let scoring = Scoring::new(GAP_OPEN, GAP_EXTEND, |a: u8, b: u8| if a == b { MATCH } else { MISMATCH });
    seqnum = 0;
    let mut aligner = Aligner::new(scoring, &homopolymer_vec[0].bases);
    for homopolymer_seq in &homopolymer_vec{
        if seqnum != 0 {
            aligner.global(&homopolymer_seq.bases).add_to_graph();
        }
        seqnum += 1;
        println!("Sequence {} processed", seqnum);
    }
    let homopolymer_consensus;
    (homopolymer_consensus, _) = aligner.poa.consensus(); //poa
    //let graph = aligner.poa.graph;
    //println!("homopolymer graph \n {:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]));
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

    //align normal consensus with expanded consensus for debugging
    println!("normal vs expanded consensus");
    let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
    let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(normal_consensus.len(), expanded_consensus.len(), GAP_OPEN, GAP_EXTEND, &score);
    let alignment = aligner.global(&normal_consensus, &expanded_consensus);
    //print_alignment_with_count(&normal_consensus,&expanded_consensus, &alignment, &homopolymer_consensus_freq);
    let (normal_alignment, expanded_alignment) = get_alignment_vectors(&normal_consensus,&expanded_consensus, &alignment);
    //write results to file
    write_scores_result_file("./results/results.txt", normal_score, homopolymer_score, expanded_score);
    write_consensus_fasta_file("./results/consensus.fa", &normal_consensus, &homopolymer_consensus, &expanded_consensus);
    write_alignment_data_fasta_file("./results/consensus.fa", &normal_alignment, &expanded_alignment);

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
    //randomvec.push(firstseq.iter().collect::<String>());
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
            let mean_value: f64 = 1.5; //2.0 before
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
fn write_alignment_data_fasta_file(filename: impl AsRef<Path>, normal_consensus: &Vec<u8>, expanded_consensus: &Vec<u8>) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n>Normal consensus vs Expanded consensus:",
        chrono::offset::Local::now(), FILENAME)
        .expect("result file cannot be written");
    let mut index = 0;
    while index < normal_consensus.len() {
        //get 50 characters
        let mut temp_vec1 = vec!();
        let mut temp_vec2 = vec!();
        for i in index..index + 150{
            if i < normal_consensus.len() {
                temp_vec1.push(normal_consensus[i]);
                temp_vec2.push(expanded_consensus[i]);
            }
            else {
                break;
            }
        }
        index = index + 150;
        writeln!(file,
            "{}\n{}\n",
            std::str::from_utf8(&temp_vec1).unwrap(), std::str::from_utf8(&temp_vec2).unwrap())
            .expect("result file cannot be written");
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


fn get_expanded_consensus(homopolymer_vec: Vec<HomopolymerSequence>, homopolymer_consensus: &Vec<u8>) -> (Vec<u8>, Vec<Vec<u32>>, i32)  {
    //use homopolymer compressions sequences to make expanded consensus //make function
    let mut homopolymer_score = 0;
    let mut homopolymer_consensus_freq: Vec<Vec<u32>> = vec![vec![0; homopolymer_vec.len()]; homopolymer_consensus.len()];
    let mut i = 0;
    for homopolymer_seq in &homopolymer_vec{
        let mut sequence_base_freq: Vec<u32> = vec![0; homopolymer_seq.bases.len()];
        let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
        let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(homopolymer_consensus.len(), homopolymer_seq.bases.len(), GAP_OPEN, GAP_EXTEND, &score);
        let alignment = aligner.global(&homopolymer_consensus, &homopolymer_seq.bases);
        let mut consensus_index = alignment.xstart;
        let mut homopolymer_index = alignment.ystart;
        //println!("start index consensus {}", consensus_index);
        //println!("start index sequence {}", homopolymer_index);
        for op in &alignment.operations {
            match op {
                bio::alignment::AlignmentOperation::Match => {
                    //println!("{} Match {}", homopolymer_consensus[consensus_index], homopolymer_seq.bases[homopolymer_index]);                    
                    homopolymer_consensus_freq[consensus_index][i] += homopolymer_seq.frequencies[homopolymer_index];
                    sequence_base_freq[homopolymer_index] += homopolymer_seq.frequencies[homopolymer_index];
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

            //print_u8_consensus(&homopolymer_consensus);
            //print_u8_consensus(&homopolymer_consensus);
        
        //println!("sequence {}", i);
        //println!("{:?}", sequence_base_freq);
        //print_alignment_with_count(homopolymer_consensus, &homopolymer_seq.bases, &alignment, &sequence_base_freq);
        homopolymer_score += alignment.score;
        i += 1;

    }

    //make the expanded consensus using the frequencies
    let mut expanded_consensus: Vec<u8> = vec![];
    let mut repetitions: Vec<f32> = vec![0.0; homopolymer_consensus.len()];
    //++ average ++
    if CONSENSUS_METHOD == 0 {
        for i in 0..homopolymer_vec.len() {
            for j in 0..homopolymer_consensus.len() {
                repetitions[j] += homopolymer_consensus_freq[j][i] as f32; 
            }
        }
        
        for j in 0..homopolymer_consensus.len() {
            for _ in 0..((repetitions[j] / homopolymer_vec.len() as f32).round() as usize) {
                expanded_consensus.push(homopolymer_consensus[j]);
            }
        }
    }
    //++ median ++ 
    if CONSENSUS_METHOD == 1 {
        //reorder the homopolymer_consensus_freq by ascending order
        for i in 0..homopolymer_consensus.len(){
            homopolymer_consensus_freq[i].sort();
            repetitions[i] = homopolymer_consensus_freq[i][(homopolymer_vec.len() / 2) as usize] as f32;
        }
        for j in 0..homopolymer_consensus.len() {
            for _ in 0..((repetitions[j]).round() as usize) {
                expanded_consensus.push(homopolymer_consensus[j]);
            }
        }
    }
    
    
    (expanded_consensus, homopolymer_consensus_freq, homopolymer_score)
}

fn get_alignment_vectors(vector1: &Vec<u8>, vector2: &Vec<u8>, alignment: &bio::alignment::Alignment) -> (Vec<u8>, Vec<u8>){
    let mut vec1_representation = vec![];
    let mut vec2_representation = vec![];
    let mut vec1_index: usize = alignment.xstart;
    let mut vec2_index: usize = alignment.ystart;
    for op in &alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                vec1_representation.push(vector1[vec1_index]);
                vec1_index += 1;
                vec2_representation.push(vector2[vec2_index]);
                vec2_index += 1;
                
            },
            bio::alignment::AlignmentOperation::Subst => {
                vec1_representation.push(vector1[vec1_index]);
                vec1_index += 1;
                vec2_representation.push(vector2[vec2_index]);
                vec2_index += 1;
                //println!("mismatch, {},{}",vec1_index, vec2_index);

            },
            bio::alignment::AlignmentOperation::Del => {
                vec1_representation.push(45);
                vec2_representation.push(vector2[vec2_index]);
                //println!("del, {},{}",vec1_index, vec2_index);
                vec2_index += 1;
               
            },
            bio::alignment::AlignmentOperation::Ins => {
                vec1_representation.push(vector1[vec1_index]);
                vec1_index += 1;
                vec2_representation.push(45);
            },
            _ => {},
        }
    }
    (vec1_representation, vec2_representation)
}
//structs here
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

//print stuff here
fn print_u8_consensus(vector: &Vec<u8>) {
    for base in vector {
        match base {
            55 => print!("_"),
            65 => print!("A"),
            67 => print!("C"),
            71 => print!("G"),
            84 => print!("T"),
            _ => {},
        }
    }
    println!("");
}

fn print_alignment_with_count(vector1: &Vec<u8>, vector2: &Vec<u8>, alignment: &bio::alignment::Alignment, count: &Vec<u32>){
    let mut vec1_representation = vec![];
    let mut vec2_representation = vec![];
    let mut vec1_index: usize = alignment.xstart;
    let mut vec2_index: usize = alignment.ystart;
    let mut prev_base: u8 = 0; //expanded consensus
    let mut same_base = false;
    let mut count_representation = vec![];
    let mut count_index = 0;
    for op in &alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                vec1_representation.push(vector1[vec1_index]);
                vec1_index += 1;
                vec2_representation.push(vector2[vec2_index]);
                if vector2[vec2_index] == prev_base {
                    same_base = true;
                }
                else {
                    same_base = false;
                }
                vec2_index += 1;
                
                
            },
            bio::alignment::AlignmentOperation::Subst => {
                vec1_representation.push(vector1[vec1_index]);
                vec1_index += 1;
                vec2_representation.push(vector2[vec2_index]);
                if vector2[vec2_index] == prev_base {
                    same_base = true;
                }
                else {
                    same_base = false;
                }
                vec2_index += 1;
                //println!("mismatch, {},{}",vec1_index, vec2_index);
                

            },
            bio::alignment::AlignmentOperation::Del => {
                vec1_representation.push(55);
                vec2_representation.push(vector2[vec2_index]);
                //println!("del, {},{}",vec1_index, vec2_index);
                if vector2[vec2_index] == prev_base {
                    same_base = true;
                }
                else {
                    same_base = false;
                }
                vec2_index += 1;
                
            },
            bio::alignment::AlignmentOperation::Ins => {
                vec1_representation.push(vector1[vec1_index]);
                vec1_index += 1;
                vec2_representation.push(55);
                same_base = true;
            },
            _ => {},
        }
        if same_base{
            count_representation.push(45);
        }
        else {
            count_representation.push(count[count_index]);
            count_index += 1;
        }
        if vec2_index != 0 {
            prev_base = vector2[vec2_index - 1];
        }
    }
    //print_u8_consensus(&vec1_representation);
    //print_u8_consensus(&vec2_representation);
    print_consensus_with_count(&vec1_representation, &vec2_representation, &count_representation);
}

fn print_consensus_with_count(normal: &Vec<u8>, expanded: &Vec<u8>, count_representation: &Vec<u32>) {
    for i in 0..normal.len(){
        match normal[i] {
            55 => print!("_"),
            65 => print!("A"),
            67 => print!("C"),
            71 => print!("G"),
            84 => print!("T"),
            _ => {},
        }
        print!(" ");
        match expanded[i] {
            55 => print!("_"),
            65 => print!("A"),
            67 => print!("C"),
            71 => print!("G"),
            84 => print!("T"),
            _ => {},
        }
        print!(" {}\n", count_representation[i]);
    }
}
fn print_count(count_representation: &Vec<u32>){
    for i in 0..count_representation.len() {
        if count_representation[i] != 45 {
            print!("{},", count_representation[i]);
        }
        else {
            print!("__,");
        }
    }
    println!("");
}