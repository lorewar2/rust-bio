use bio::alignment::{poa::*, TextSlice, pairwise::Scoring};
use std::{fs, fs::File, fs::OpenOptions, io::{prelude::*, BufReader}, path::Path, fmt, cmp, collections::HashMap};
use chrono;
use rand::{Rng, SeedableRng, rngs::StdRng, seq};
use petgraph::{Directed, Graph, Incoming, Outgoing, Direction, dot::Dot, graph::NodeIndex, visit::Topo};
use statrs::function::factorial::binomial;
use logaddexp::LogAddExp;
use libm::exp;
use queues::*;

const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = -3;
const MATCH: i32 = 2;
const MISMATCH: i32 = -4;
const SEED: u64 = 6;
const CONSENSUS_METHOD: u8 = 1; //0==average 1==median //2==mode
const ERROR_PROBABILITY: f64 = 0.85;
const HOMOPOLYMER_DEBUG: bool = false;
const HOMOPOLYMER: bool = false;
const QUALITY_SCORE: bool = false;
const NUM_OF_ITER_FOR_PARALLEL: usize = 10;
const NUM_OF_ITER_FOR_ZOOMED_GRAPHS: usize = 4;
const USEPACBIODATA: bool = false;
const PACBIOALLFILES: bool = false;
const USER_DEFINED: bool = false;
const ERROR_LINE_NUMBER: usize = 10; //default 10
const PRINT_ALL: bool = true;
const RANDOM_SEQUENCE_LENGTH: usize = 10;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 2;

// file names input
const INPUT_FILE_NAME: &str = "11928566";
const INPUT_READ_FOLDER_PATH: &str = "./data/PacBioReads/";
const INPUT_CONSENSUS_FOLDER_PATH: &str = "./data/PacBioConsensus/";
const INPUT_PACBIO_ERROR_FOLDER_PATH: &str = "./data/PacBioError/";


// file names output
const OUTPUT_RESULT_PATH: &str = "./results/";
fn main() {
    if PACBIOALLFILES {
        // get the file names of reads
        let mut file_name_vec: Vec<(String, usize)> = vec![];
        let paths = fs::read_dir(INPUT_PACBIO_ERROR_FOLDER_PATH).unwrap();
        let mut total_length = 0;
        for path in paths {
            let mut temp = path.unwrap().file_name().to_string_lossy().to_string();
            match temp.find("pdf") {
                Some(_) => {
                    temp.truncate(temp.len() - 8);
                    let split = temp.split("_").collect::<Vec<&str>>();
                    file_name_vec.push((split[0].to_string(), split[1].to_string().parse::<usize>().unwrap()));
                    total_length += 1;
                },
                None => {}
            };
        }
        println!("{:?}", file_name_vec);
        let mut index = 1;
        // run all reads
        for file_name in file_name_vec {
            let mut seqvec;
            println!("Running read {}, {}/{}", file_name.0, index, total_length);
            //create a folder 
            fs::create_dir([OUTPUT_RESULT_PATH, &file_name.0].concat()).ok();
            //create corrosponding result files
            let (output_debug_file_name, output_consensus_file_name, output_scores_file_name, output_normal_graph_file_name, output_homopolymer_graph_file_name, output_quality_graph_file_name, output_quality_file_name) = create_required_result_files(&[OUTPUT_RESULT_PATH, &file_name.0].concat());
            let pac_bio_consensus =  get_consensus_from_file([INPUT_CONSENSUS_FOLDER_PATH, &file_name.0, ".fastq"].concat());
            seqvec = get_fasta_sequences_from_file([INPUT_READ_FOLDER_PATH, &file_name.0, ".fasta"].concat());
            seqvec = check_the_scores_and_change_alignment(seqvec, &pac_bio_consensus);
            seqvec.insert(0, pac_bio_consensus);
            run(seqvec, [INPUT_CONSENSUS_FOLDER_PATH, &file_name.0, ".fastq"].concat().to_string(), output_debug_file_name, output_consensus_file_name, output_scores_file_name, output_normal_graph_file_name, output_homopolymer_graph_file_name, output_quality_graph_file_name, output_quality_file_name, file_name.1);
            index += 1;
        }
    }
    else {
        let mut seqvec = vec![];
        // read the pacbio consensus
        if USEPACBIODATA {
            //create a folder 
            fs::create_dir([OUTPUT_RESULT_PATH, INPUT_FILE_NAME].concat()).ok();
            //create corrosponding result files
            let (output_debug_file_name, output_consensus_file_name, output_scores_file_name, output_normal_graph_file_name, output_homopolymer_graph_file_name, output_quality_graph_file_name, output_quality_file_name) = create_required_result_files(&[OUTPUT_RESULT_PATH, INPUT_FILE_NAME].concat());
            let pac_bio_consensus =  get_consensus_from_file([INPUT_CONSENSUS_FOLDER_PATH, INPUT_FILE_NAME, ".fastq"].concat());
            seqvec = get_fasta_sequences_from_file([INPUT_READ_FOLDER_PATH, INPUT_FILE_NAME, ".fasta"].concat());
            seqvec = check_the_scores_and_change_alignment(seqvec, &pac_bio_consensus);
            seqvec.insert(0, pac_bio_consensus);
            run(seqvec, [INPUT_CONSENSUS_FOLDER_PATH, INPUT_FILE_NAME, ".fastq"].concat().to_string(), output_debug_file_name, output_consensus_file_name, output_scores_file_name, output_normal_graph_file_name, output_homopolymer_graph_file_name, output_quality_graph_file_name, output_quality_file_name, ERROR_LINE_NUMBER);
        }
        else if USER_DEFINED {
            //create a folder 
            fs::create_dir([OUTPUT_RESULT_PATH, "user"].concat()).ok();
            //create corrosponding result files
            let (output_debug_file_name, output_consensus_file_name, output_scores_file_name, output_normal_graph_file_name, output_homopolymer_graph_file_name, output_quality_graph_file_name, output_quality_file_name) = create_required_result_files(&[OUTPUT_RESULT_PATH, "user"].concat());
            //seqvec.push("SCATGGTGGCTTACTTCGTCAAAGCATGCAAGCCAAGAAGGCGACAGAGAGAGAGAGAGAGAGGGAGAE".to_string());
            seqvec.push("SCATGGTGGCTTACTTCGTACAAAGCATGCAAGCCAAGAAGGCGACAGAGCAGAAGAGAGTAGAGGGAGAE".to_string());
            seqvec.push("SCATGGTGGCTTACTTCGTCAAGCATGCAAGCAAGAAGGCGACAGAGAGAAAGGAAGAGAGAGGGAGAE".to_string());
            seqvec.push("SCATGGGTGGCTTACTTCGTCAAAGCATGCAAGCCAAGAAGCGACAGAGAGTAAAGAGAGAGAGAGGTGAE".to_string());
            seqvec.push("SCATGGTGGCTTACTTCGTCAAAAGATGCAAGCCAAGAAGGCGACAGAGAGAAAGAGAGAGGAGGGAGAE".to_string());
            seqvec.push("SCAAGTGGCTTACTTCGTCAAAGCATGCAAGCCAGAAGCGACGAGAGAAAGAGAGAGAGGGAGAE".to_string());
            seqvec.push("SCATGTGCGCTTACTTCGTCAAAGCATGCAAGGCCAGAAGGCGACAGACGAGAAAGAGAGAGAGGGAGAE".to_string());
            seqvec.push("SCATGGTGGCTTTACTTTCGTCAAAGCATTCAAGCCAAGAACGAACAGAGAGAAAGAGAGAGAGGGGAAE".to_string());
            seqvec.push("SACCATGGTGCTTACTGCAGTTCAAAGCATTGCAAGCCAAGAAGGCGACAGAGAGGAAAGAGATGAGAGGGAGGAE".to_string());
            seqvec.push("SCATGGTGCTTACTTCGTCAAAGCATGCAGCCCAAGAAGGAGCGACAGGAGAGAAAGAGAGGAGAGGAGAE".to_string());
            seqvec.push("SCATGGTGGCTTACTTCGGTCAAAGCATTGCAAGCCAAGAAGGCGACAGAGAGAAGAGAAAGAGAGAGAGGGAGAE".to_string());
            seqvec.push("SCTTGGTGCTTACTCGGCAAAGCATGATCAAAGCCAAGAAGGCGACAGAGAAAAGAGAGAGAGGGAGAE".to_string());
            seqvec.push("SCTGGTGGACTTACTTCGTCAAAGCATGCAAGCCAAGAAGGCCGACAGAGAGTAAAGAAGAGAGAGGGAGAE".to_string());
            seqvec.push("SCATGGTGCTTACTTCGTCAAAGTCATGCAAGCCAAGAAGGCGACAGAGAGAAAGAGAGAGAGGGAGAGE".to_string());
            seqvec.push("SCATGGTGGGCTTACTTCGTCAAGCATGCAAGCCACAGAAGGCGACAGAGAGAAGAGAGAGAGGGAGAE".to_string());
            seqvec.push("SCATGTGGCTTACTTCGTCCAAGCATGCAAGCCAAAAGGCGACAGAGAGAAAGAGAGAGAGGGGAGAE".to_string());
            run(seqvec, [INPUT_CONSENSUS_FOLDER_PATH, INPUT_FILE_NAME].concat().to_string(), output_debug_file_name, output_consensus_file_name, output_scores_file_name, output_normal_graph_file_name, output_homopolymer_graph_file_name, output_quality_graph_file_name, output_quality_file_name, ERROR_LINE_NUMBER);
        }
        else {
            //create a folder 
            fs::create_dir([OUTPUT_RESULT_PATH, "random"].concat()).ok();
            //create corrosponding result files
            let (output_debug_file_name, output_consensus_file_name, output_scores_file_name, output_normal_graph_file_name, output_homopolymer_graph_file_name, output_quality_graph_file_name, output_quality_file_name) = create_required_result_files(&[OUTPUT_RESULT_PATH, "random"].concat());
            seqvec = get_random_sequences_from_generator(RANDOM_SEQUENCE_LENGTH, NUMBER_OF_RANDOM_SEQUENCES);
            run(seqvec, [INPUT_CONSENSUS_FOLDER_PATH, INPUT_FILE_NAME].concat().to_string(), output_debug_file_name, output_consensus_file_name, output_scores_file_name, output_normal_graph_file_name, output_homopolymer_graph_file_name, output_quality_graph_file_name, output_quality_file_name, ERROR_LINE_NUMBER);
        }
    }
}

fn run (seqvec: Vec<String>, input_consensus_file_name: String, output_debug_file_name: String, 
    output_consensus_file_name: String, output_scores_file_name: String, output_normal_graph_file_name: String, 
    output_homopolymer_graph_file_name: String, output_quality_graph_file_name: String, output_quality_file_name: String, error_line_number: usize) {
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
    let normal_topology;
    (normal_consensus, normal_topology) = aligner.poa.consensus(); //just poa

    // get scores of sequences compared to normal consensus 
    let normal_score = get_consensus_score(&seqvec, &normal_consensus);
    // get the normal graph
    let normal_graph = aligner.graph();
    
    if HOMOPOLYMER {
    ////////////////////////////
    //compressed poa alignment//
    ////////////////////////////
        // make homopolymer compression vector
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
        let homopolymer_topology;
        (homopolymer_consensus, homopolymer_topology) = aligner.poa.consensus();
        // get graph
        let homopolymer_graph: &Graph<u8, i32, Directed, usize> = aligner.graph();
        // use homopolymer compressions sequences to make expanded consensus
        let (homopolymer_consensus_freq, homopolymer_score) = get_aligned_homopolymer_sequences_to_homopolymer_consensus(&homopolymer_vec, &homopolymer_consensus);
        let (expanded_consensus, homopolymer_expanded) =  calculate_and_get_expansion (&homopolymer_vec, &homopolymer_consensus, &homopolymer_consensus_freq);
        // get the scores of expanded consensus compared to sequences
        let expanded_score = get_consensus_score(&seqvec, &expanded_consensus);
        println!("expanded score:{}", expanded_score);

        if HOMOPOLYMER_DEBUG {
            let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
            let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(normal_consensus.len(), expanded_consensus.len(), GAP_OPEN, GAP_EXTEND, &score);
            let alignment = aligner.global(&normal_consensus, &expanded_consensus);
            let mut saved_indices: IndexStruct;
            saved_indices = get_indices_for_homopolymer_debug(&alignment, &homopolymer_expanded, &normal_topology, &homopolymer_topology);
            let (normal_rep, expanded_rep, count_rep) 
                = get_alignment_with_count_for_debug (&normal_consensus,&expanded_consensus, &alignment, &homopolymer_consensus_freq, seqnum as usize);
            //write results to file
            write_scores_result_file(output_scores_file_name, normal_score, homopolymer_score, expanded_score);
            //modify the graphs to indicate 
            saved_indices = modify_and_write_the_graphs_and_get_zoomed_graphs(output_normal_graph_file_name.clone(), output_homopolymer_graph_file_name, saved_indices, normal_graph, homopolymer_graph);
            write_alignment_and_zoomed_graphs_fasta_file(output_consensus_file_name, &normal_rep, &expanded_rep,
                &count_rep, seqnum as usize, &saved_indices);
        }
    }
    /////////////////////////////
    //quality score calculation//
    /////////////////////////////
    if QUALITY_SCORE {
        let mut invalid_info: Vec<(usize, usize, bool, bool, bool, bool)> = vec![]; //index, node_id, parallel invalid, match invalid, quality invalid, error pacbio
        // calculate and get the quality scores
        let (quality_scores, parallel_validity, base_count_vec, mut debug_strings) = get_consensus_quality_scores(seqnum as usize, &normal_consensus, &normal_topology, normal_graph);
        if USEPACBIODATA {
            let mut parallel_count = 0;
            let mut mismatch_count = 0;
            let mut quality_count = 0;
            let (pacbio_quality_scores, mismatch_indices, aligned_pacbio_bases, calc_error_line_number) = get_quality_score_aligned (get_consensus_from_file(input_consensus_file_name.clone()), &normal_consensus, get_quality_from_file(input_consensus_file_name.clone()), error_line_number);
            // get invalid indices is quality score is too low
            for index in 0..quality_scores.len() {
                if parallel_validity[index] == true {
                    invalid_info.push((index, normal_topology[index], true, false, false, false));
                    parallel_count += 1;
                }
                if mismatch_indices.contains(&index) {
                    match invalid_info.iter().position(|&r| r.0 == index) {
                        Some(position) => {invalid_info[position].3 = true;},
                        None => {invalid_info.push((index, normal_topology[index], false, true, false, false));},
                    }
                    mismatch_count += 1;
                }
                if quality_scores[index] < 10.0 {
                    match invalid_info.iter().position(|&r| r.0 == index) {
                        Some(position) => {invalid_info[position].4 = true;},
                        None => {invalid_info.push((index, normal_topology[index], false, false, true, false));},
                    }
                    quality_count += 1;
                }
                if index == calc_error_line_number {
                    match invalid_info.iter().position(|&r| r.0 == index) {
                        Some(position) => {invalid_info[position].5 = true;},
                        None => {invalid_info.push((index, normal_topology[index], false, false, false, true));},
                    }
                }
            }
            println!("INVALID COUNT: {} parallel err: {} mismatch err: {} quality err: {}", invalid_info.len(), parallel_count, mismatch_count, quality_count);
            debug_strings.push(format!("INVALID COUNT: {} parallel err: {} mismatch err: {} quality err: {}", invalid_info.len(), parallel_count, mismatch_count, quality_count));
            // write the zoomed in graphs for invalid and low quality entries.
            if invalid_info.len() < 500 {
                let temp_quality_error = write_quality_scores_to_file(output_quality_file_name, &quality_scores, &normal_consensus, &normal_topology, &invalid_info, &base_count_vec, &pacbio_quality_scores, &aligned_pacbio_bases);
                let temp_graph_error = write_zoomed_quality_score_graphs (output_quality_graph_file_name, &invalid_info, &quality_scores, &base_count_vec, normal_graph, &pacbio_quality_scores);
                debug_strings = [debug_strings, temp_quality_error, temp_graph_error].concat();
                write_debug_data_to_file(output_debug_file_name, debug_strings);
            }
        }
        else {
            for index in 0..quality_scores.len() {
                if parallel_validity[index] == true {
                    invalid_info.push((index, normal_topology[index], true, false, false, false));
                }
                if index == error_line_number {
                    match invalid_info.iter().position(|&r| r.0 == index) {
                        Some(position) => {invalid_info[position].5 = true;},
                        None => {invalid_info.push((index, normal_topology[index], false, false, false, true));},
                    }
                }
            }
            println!("INVALID COUNT: {}", invalid_info.len());
            if invalid_info.len() < 100 {
                let temp_quality_error = write_quality_scores_to_file(output_quality_file_name, &quality_scores, &normal_consensus, &normal_topology, &invalid_info, &base_count_vec, &vec![34, 34], &vec![65; quality_scores.len()]);
                let temp_graph_error = write_zoomed_quality_score_graphs (output_quality_graph_file_name, &invalid_info, &quality_scores, &base_count_vec, normal_graph, &vec![34, 34]);
                debug_strings = [debug_strings, temp_quality_error, temp_graph_error].concat();
                write_quality_score_graph(output_normal_graph_file_name, normal_graph);
                write_debug_data_to_file(output_debug_file_name, debug_strings);
            }
        }
    }
    /////////////////////////////
    //alternate aligners       //
    /////////////////////////////
    let (aligned_dp, dp_score) = normal_dp (&seqvec[0].as_bytes().to_vec(), &seqvec[1].as_bytes().to_vec());
    println!("consensus: {:?} score {}", aligned_dp, dp_score);

    let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
    let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(seqvec[0].len(), seqvec[1].len(), GAP_OPEN, GAP_EXTEND, &score);
    let alignment = aligner.global(&seqvec[0].as_bytes().to_vec(), &seqvec[1].as_bytes().to_vec());
    println!("score of original {}", alignment.score);

    //convert the two sequences to homopolymer vec
    let mut homopolymer_vec_x: Vec<HomopolymerCell> = convert_sequence_to_homopolymer (seqvec[0].clone());
    let mut homopolymer_vec_y: Vec<HomopolymerCell> = convert_sequence_to_homopolymer (seqvec[1].clone());

    for base in &homopolymer_vec_x {
        //println!("{} {}", base.base, base.frequency);
    }
    println!("");
    for base in &homopolymer_vec_y {
        //println!("{} {}", base.base, base.frequency);
    }
    //homopolymer_dp (&homopolymer_vec_x, &homopolymer_vec_y);
    // divide the sequence in to two 
    /* 
    divide_pipeline(&seqvec);
    let bfs_con = bfs_consensus(normal_graph, &seqvec);
    let bfs_con_score = get_consensus_score(&seqvec, &bfs_con);

    for base in &bfs_con {
        print!("{}", *base as char);
    }
    println!("");
    for base in &normal_consensus {
        print!("{}", *base as char);
    }
    println!("");
    println!("normal score: {}", normal_score);
    println!("bfs score: {}", bfs_con_score);
    //println!("{}", format!("{:?}", Dot::new(&normal_graph.map(|_, n| (*n) as char, |_, e| *e)))); 
    let (topology_consensus, _) = topology_cut_consensus(&seqvec);
    let topology_score = get_consensus_score(&seqvec, &topology_consensus);
    //let (mod_heavy_consensus, _) = heavy_bundle_modified_consensus(&seqvec);
    //let mod_heavy_score = get_consensus_score(&seqvec, &mod_heavy_consensus);
    println!("normal score: {}", normal_score);
    println!("topo score: {}", topology_score);
    //println!("mod heavy score: {}", mod_heavy_score);
    // score of calcualted 
    println!("normal score:\t{}", normal_score);
    // score of pacbio
    let pacbio_score = get_consensus_score(&seqvec, &seqvec[0].as_bytes().to_vec());
    println!("pacbio score:\t{}", pacbio_score);
    // score of pacbio error corrected
    let mut modified_pacbio = seqvec[0].as_bytes().to_vec();
    modified_pacbio[15157] = 65;
    let modified_pacbio_score = get_consensus_score(&seqvec, &modified_pacbio);
    println!("modified pacbio score:\t{}", modified_pacbio_score);
    */
}

pub fn convert_sequence_to_homopolymer (sequence: String) -> Vec<HomopolymerCell> {
    let mut homopolymer_vec: Vec<HomopolymerCell> = vec![];
    let mut prev_base: u8 = 0;
    let mut frequency: usize = 1;
    for base in sequence.as_bytes().to_vec() {
        if prev_base == base {
            frequency += 1;
        }
        else if prev_base != 0 {
            homopolymer_vec.push(HomopolymerCell::new(prev_base, frequency));
            frequency = 1;
        }
        prev_base = base;
    }
    homopolymer_vec.push(HomopolymerCell::new(prev_base, frequency));
    homopolymer_vec
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

pub struct HomopolymerCell {
    pub base: u8,
    pub frequency: usize,
}

impl HomopolymerCell {
    fn new(base: u8, frequency: usize) -> Self{
        HomopolymerCell {
            base: base,
            frequency: frequency,
        }
    }
}

fn homopolymer_dp (homo_x: &Vec<HomopolymerCell>, homo_y: &Vec<HomopolymerCell>) -> (Vec<u8>, isize) {
    let mut align_vec: Vec<u8> = Vec::new();
    // calculation variables
    // initialize the weight matrix with (0, and size)
    let mut match_matrix: Vec<Vec<(isize, usize)>> = vec![vec![(0, 0); homo_y.len() + 1]; homo_x.len() + 1]; // match or mismatch diagonal edges
    let mut del_matrix: Vec<Vec<(isize, usize)>> = vec![vec![(0, 0); homo_y.len() + 1]; homo_x.len() + 1];  // x deletions right direction edges
    let mut ins_matrix: Vec<Vec<(isize, usize)>> = vec![vec![(0, 0); homo_y.len() + 1]; homo_x.len() + 1]; // x insertion down direction edges
    // initialize the backtrace matrix with ms, ds, and is
    let mut back_matrix: Vec<Vec<char>> = vec![vec!['m'; homo_y.len() + 1]; homo_x.len() + 1];
    for i in 1..homo_y.len() + 1 {
        back_matrix[0][i] = 'd';
        let temp_value = del_matrix[0][i - 1].0 + (GAP_EXTEND as isize) * (homo_y[i - 1].frequency as isize);
        del_matrix[0][i].0 = temp_value; ins_matrix[0][i].0 = temp_value; match_matrix[0][i].0 = temp_value;
        del_matrix[0][i].1 = homo_y[i - 1].frequency; ins_matrix[0][i].1 = homo_y[i - 1].frequency; match_matrix[0][i].1 = homo_y[i - 1].frequency;
    }
    for i in 1..homo_x.len() + 1 {
        back_matrix[i][0] = 'i';
        let temp_value = ins_matrix[i - 1][0].0 + (GAP_EXTEND as isize) * (homo_x[i - 1].frequency as isize);
        ins_matrix[i][0].0 = temp_value; del_matrix[i][0].0 = temp_value; match_matrix[i][0].0 = temp_value;
        ins_matrix[i][0].1 = homo_x[i - 1].frequency; del_matrix[i][0].1 = homo_x[i - 1].frequency; match_matrix[i][0].1 = homo_x[i - 1].frequency;
    }
    // calculations
    // filling out score matrices and back matrix
    for i in 1..homo_x.len() + 1 {
        for j in 1..homo_y.len() + 1 {
            // fill del matrix 
            // get j - 1 score from same matrix with gap extend
            let temp_del_score = del_matrix[i][j - 1].0 + ((GAP_EXTEND as isize) * (homo_y[j - 1].frequency as isize));
            // get j - 1 score from match matrix with gap open penalty
            let temp_match_score = match_matrix[i][j - 1].0 + GAP_OPEN as isize + ((GAP_EXTEND as isize) * ((homo_y[j - 1].frequency - 1) as isize));
            // insert the max
            del_matrix[i][j].0 = cmp::max(temp_del_score, temp_match_score);
            del_matrix[i][j].1 = homo_y[j - 1].frequency;

            // fill ins matrix
            // get i - 1 score from the same matrix
            let temp_ins_score = ins_matrix[i - 1][j].0 + ((GAP_EXTEND as isize) * (homo_x[i - 1].frequency as isize));
            // get i - 1 score from the match matrix with gap open penalty
            let temp_match_score = match_matrix[i - 1][j].0 + GAP_OPEN as isize + ((GAP_EXTEND as isize) * ((homo_x[i - 1].frequency - 1) as isize));
            // insert the max
            ins_matrix[i][j].0 = cmp::max(temp_ins_score, temp_match_score);
            ins_matrix[i][j].1 = homo_x[i - 1].frequency;

            // fill match matrix
            // get the i,j from the insertion matrix
            let temp_ins_score = ins_matrix[i][j].0;
            // get the i,j from the deletion matrix
            let temp_del_score = del_matrix[i][j].0;
            // get the match from i-1,j-1 from match matrix with match score or mismatch score
            let temp_match_score;
            if homo_x[i - 1].base == homo_y[j - 1].base {
                temp_match_score = match_matrix[i - 1][j - 1].0 + ((MATCH as isize) * (cmp::min(homo_x[i - 1].frequency as isize, homo_y[j - 1].frequency as isize)));
            }
            else {
                temp_match_score = match_matrix[i - 1][j - 1].0 + ((MISMATCH as isize) * (cmp::min(homo_x[i - 1].frequency as isize, homo_y[j - 1].frequency as isize)));
            }
            // insert the max
            match_matrix[i][j].0 = cmp::max(temp_match_score, cmp::max(temp_ins_score, temp_del_score));
            if (temp_match_score >= temp_ins_score) && (temp_match_score >= temp_del_score) {
                back_matrix[i][j] = 'm';
                match_matrix[i][j].1 = cmp::min(homo_x[i - 1].frequency, homo_y[j - 1].frequency);
            }
            else if temp_ins_score > temp_del_score {
                back_matrix[i][j] = 'i';
                match_matrix[i][j].1 = homo_x[i - 1].frequency;
            }
            else {
                back_matrix[i][j] = 'd';
                match_matrix[i][j].1 = homo_y[j - 1].frequency;
            }
        }
    }
    println!("ins score");
    for i in 0..homo_x.len() + 1 {
        for j in 0..homo_y.len() + 1 {
            print!("{:>3}", ins_matrix[i][j].0);
        }
        println!("");
    }
    println!("ins freq");
    for i in 0..homo_x.len() + 1 {
        for j in 0..homo_y.len() + 1 {
            print!("{:>3}", ins_matrix[i][j].1);
        }
        println!("");
    }
    println!("mat score");
    for i in 0..homo_x.len() + 1 {
        for j in 0..homo_y.len() + 1 {
            print!("{:>3}", match_matrix[i][j].0);
        }
        println!("");
    }
    println!("mat freq");
    for i in 0..homo_x.len() + 1 {
        for j in 0..homo_y.len() + 1 {
            print!("{:>3}", match_matrix[i][j].1);
        }
        println!("");
    }
    println!("del score");
    for i in 0..homo_x.len() + 1 {
        for j in 0..homo_y.len() + 1 {
            print!("{:>3}", del_matrix[i][j].0);
        }
        println!("");
    }
    println!("del freq");
    for i in 0..homo_x.len() + 1 {
        for j in 0..homo_y.len() + 1 {
            print!("{:>3}", del_matrix[i][j].1);
        }
        println!("");
    }
    // back tracing using back matrix and filling out align_vec
    let mut i = homo_x.len();
    let mut j = homo_y.len();
    let score = match_matrix[i][j].0;
    loop {
        match back_matrix[i][j] {
            'i' => {
                i = i - 1;
                for _ in 0..homo_x[i].frequency {
                    align_vec.push(homo_x[i].base);
                }
            },
            'm' => {
                i = i - 1;
                j = j - 1;
                for _ in 0..cmp::min(homo_x[i].frequency, homo_y[i].frequency) {
                    align_vec.push(homo_x[i].base);
                }
            },
            'd' => {
                j = j - 1;
                for _ in 0..homo_y[j].frequency {
                    align_vec.push(homo_y[j].base);
                }
            },
            _ => (),
        }
        if i == 0 && j == 0 {
            break;
        }
    }
    (align_vec, score)
}

fn normal_dp (seq_x: &Vec<u8>, seq_y: &Vec<u8>) -> (Vec<u8>, isize) {
    // variables to save results
    let mut align_vec: Vec<u8> = Vec::new();

    // calculation variables
    // initialize the weight matrix with zeros
    let mut match_matrix: Vec<Vec<isize>> = vec![vec![0; seq_y.len() + 1]; seq_x.len() + 1]; // match or mismatch diagonal edges
    let mut del_matrix: Vec<Vec<isize>> = vec![vec![0; seq_y.len() + 1]; seq_x.len() + 1];  // x deletions right direction edges
    let mut ins_matrix: Vec<Vec<isize>> = vec![vec![0; seq_y.len() + 1]; seq_x.len() + 1]; // x insertion down direction edges
    // initialize the backtrace matrix with ms, ds, and is
    let mut back_matrix: Vec<Vec<char>> = vec![vec!['m'; seq_y.len() + 1]; seq_x.len() + 1];
    back_matrix[0][1] = 'd'; let temp_value = GAP_OPEN as isize; del_matrix[0][1] = temp_value; ins_matrix[0][1] = temp_value; match_matrix[0][1] = temp_value;
    back_matrix[1][0] = 'i'; let temp_value = GAP_OPEN as isize; ins_matrix[1][0] = temp_value; del_matrix[1][0] = temp_value; match_matrix[1][0] = temp_value;
    for i in 2..seq_y.len() + 1 {back_matrix[0][i] = 'd'; let temp_value = del_matrix[0][i - 1] + GAP_EXTEND as isize; del_matrix[0][i] = temp_value; ins_matrix[0][i] = temp_value; match_matrix[0][i] = temp_value;}
    for i in 2..seq_x.len() + 1 {back_matrix[i][0] = 'i'; let temp_value = ins_matrix[i - 1][0] + GAP_EXTEND as isize; ins_matrix[i][0] = temp_value; del_matrix[i][0] = temp_value; match_matrix[i][0] = temp_value;}

    // calculations
    // filling out score matrices and back matrix
    for i in 1..seq_x.len() + 1 {
        for j in 1..seq_y.len() + 1 {
            // fill del matrix 
            // get j - 1 score from same matrix with gap extend
            let temp_del_score = del_matrix[i][j - 1] + GAP_EXTEND as isize;
            // get j - 1 score from match matrix with gap open penalty
            let temp_match_score = match_matrix[i][j - 1] + GAP_OPEN as isize;
            // insert the max
            del_matrix[i][j] = cmp::max(temp_del_score, temp_match_score);

            // fill ins matrix
            // get i - 1 score from the same matrix
            let temp_ins_score = ins_matrix[i - 1][j] + GAP_EXTEND as isize;
            // get i - 1 score from the match matrix with gap open penalty
            let temp_match_score = match_matrix[i - 1][j] + GAP_OPEN as isize;
            // insert the max
            ins_matrix[i][j] = cmp::max(temp_ins_score, temp_match_score);

            // fill match matrix
            // get the i,j from the insertion matrix
            let temp_ins_score = ins_matrix[i][j];
            // get the i,j from the deletion matrix
            let temp_del_score = del_matrix[i][j];
            // get the match from i-1,j-1 from match matrix with match score or mismatch score
            let temp_match_score;
            if seq_x[i - 1] == seq_y[j - 1] {
                temp_match_score = match_matrix[i - 1][j - 1] + MATCH as isize;
            }
            else {
                temp_match_score = match_matrix[i - 1][j - 1] + MISMATCH as isize;
            }
            // insert the max
            match_matrix[i][j] = cmp::max(temp_match_score, cmp::max(temp_ins_score, temp_del_score));
            if (temp_match_score >= temp_ins_score) && (temp_match_score >= temp_del_score) {
                back_matrix[i][j] = 'm';
            }
            else if temp_ins_score > temp_del_score {
                back_matrix[i][j] = 'i';
            }
            else {
                back_matrix[i][j] = 'd';
            }
        }
    }
    println!("ins score");
    for i in 0..seq_x.len() + 1 {
        for j in 0..seq_y.len() + 1 {
            print!("{:>3}", ins_matrix[i][j]);
        }
        println!("");
    }
    println!("match score");
    for i in 0..seq_x.len() + 1 {
        for j in 0..seq_y.len() + 1 {
            print!("{:>3}", match_matrix[i][j]);
        }
        println!("");
    }
    println!("ins score");
    for i in 0..seq_x.len() + 1 {
        for j in 0..seq_y.len() + 1 {
            print!("{:>3}", del_matrix[i][j]);
        }
        println!("");
    }
    // back tracing using back matrix and filling out align_vec
    let mut i = seq_x.len();
    let mut j = seq_y.len();
    let score = match_matrix[i][j];
    loop {
        match back_matrix[i][j] {
            'i' => {
                i = i - 1;
                align_vec.push(seq_x[i]);
            },
            'm' => {
                i = i - 1;
                j = j - 1;
                align_vec.push(seq_x[i]);
            },
            'd' => {
                j = j - 1;
                align_vec.push(seq_y[j]);
            },
            _ => (),
        }
        if i == 0 && j == 0 {
            break;
        }
    }
    (align_vec.into_iter().rev().collect(), score)
}

pub fn divide_pipeline (seqvec: &Vec<String>) {
    let mut original_sequences = vec![]; //sequence, number sequence
    for seq in seqvec{
        original_sequences.push((*seq.clone().as_bytes()).to_vec());
    }
    // get middle k-mer if available search for a 5-mer in middle 20 while increasing the mismatch count to 2, 
    let mut test_sequence = vec![];
    for seq in &original_sequences {
        let seq_mid = seq.len() / 2;
        test_sequence = [test_sequence, seq[seq_mid - 10..seq_mid + 10].to_vec()].concat();
    }
    //println!("{:?}", test_sequence);
    // go through all the possibilities max k_mer with min required distance and max frequency in a loop trial and error
    // hamming distance one, maximize k, while only one frequency available
    let mut frequent_kmer;
    let mut frequent_indices;
    let mut k_mer: usize = 5;
    loop {
        (frequent_kmer, frequent_indices) = divide_find_frequent_kmers_with_mismatches(&test_sequence, k_mer, 1);
        if frequent_indices.len() == 1 && frequent_indices[0].len() <= seqvec.len() || k_mer >= 9 {
            break;
        }
        println!("Searching k_mer size... {}", k_mer);
        k_mer += 1;
    }
    //println!("{:?} {:?}", frequent_kmer, frequent_indices);
    // find the missing sections kmer by increasing the required_distance
    // what are the missing sequences
    let mut valid_indices = vec![false; seqvec.len()];
    let mut final_indices: Vec<usize> = vec![0; seqvec.len()];
    for index in 0..frequent_indices[0].len() {
        if frequent_indices[0][index] % 20 <= (20 - k_mer) {
            valid_indices[frequent_indices[0][index] / 20] = true;
            final_indices[frequent_indices[0][index] / 20] = frequent_indices[0][index] % 20;
        }
    }
    //println!("valid indices {:?} {:?}", valid_indices, final_indices);
    // find them make a function to search for a string by increasing the hamming distance
    for index in 0..valid_indices.len() {
        if !valid_indices[index] {
            // get the middle
            let seq_mid = original_sequences[index].len() / 2;
            let seq_part = original_sequences[index][seq_mid - 10..seq_mid + 10].to_vec();
            // find the k_mer index
            let obtained_index = divide_find_pattern_in_sequence_with_mismatches (&seq_part, &frequent_kmer[0]);
            // put it in the vector
            final_indices[index] = obtained_index;
        }
    }
    //println!("valid indices {:?} {:?}", valid_indices, final_indices);
    // from the kmers divide the sequences
    // put the divided slices into vector
    let mut first_half_slices = vec![];
    let mut second_half_slices = vec![];
    for index in 0..original_sequences.len() {
        let first_break_pos = (original_sequences[index].len() / 2) - 10;
        let first_half = [vec![83] , original_sequences[index][0..first_break_pos + final_indices[index]].to_vec(), vec![69]].concat();
        let second_half = [vec![83] ,original_sequences[index][first_break_pos + final_indices[index] + frequent_kmer[0].len()..original_sequences[index].len()].to_vec(), vec![69]].concat();;
        first_half_slices.push(first_half);
        second_half_slices.push(second_half);
    }
    // attach anchors S E

    // do poa and get consensus of the slices
    let scoring = Scoring::new(GAP_OPEN, GAP_EXTEND, |a: u8, b: u8| if a == b { MATCH } else { MISMATCH });
    let mut seqnum: u8 = 0;
    let mut aligner = Aligner::new(scoring, &first_half_slices[0]);
    for seq in &first_half_slices{
        if seqnum != 0 {
            aligner.global(seq).add_to_graph();
        }
        seqnum += 1;
        println!("Sequence {} processed", seqnum);
    }
    let firsthalf_consensus;
    (firsthalf_consensus, _) = aligner.poa.consensus(); //just poa

    let scoring = Scoring::new(GAP_OPEN, GAP_EXTEND, |a: u8, b: u8| if a == b { MATCH } else { MISMATCH });
    let mut seqnum: u8 = 0;
    let mut aligner = Aligner::new(scoring, &second_half_slices[0]);
    for seq in &second_half_slices{
        if seqnum != 0 {
            aligner.global(seq).add_to_graph();
        }
        seqnum += 1;
        println!("Sequence {} processed", seqnum);
    }
    let secondhalf_consensus;
    (secondhalf_consensus, _) = aligner.poa.consensus(); //just poa
    // get the full consensus by adding kmer and slices
    let divide_consensus = [firsthalf_consensus, frequent_kmer[0].clone(), secondhalf_consensus].concat();
    
    // remove S and E
    let mut modified_divide_consensus = vec![];
    for base in &divide_consensus {
        match base {
            83 => {},
            69 => {},
            x => {modified_divide_consensus.push(*x)},
        }
    }

    // check the scores of one division
    let divide_consensus_score = get_consensus_score(&seqvec, &modified_divide_consensus);
    for base in &modified_divide_consensus {
        print!("{}", *base as char);
    }
    println!("");
    println!("divide score {}", divide_consensus_score);
}

pub fn divide_find_pattern_in_sequence_with_mismatches (sequence: &Vec<u8>, pattern: &Vec<u8>) -> usize {
    let mut required_distance = 0;
    let location_index: usize;
    'bigloop: loop {
        for index in 0..(sequence.len() - pattern.len() + 1) {
            let current_pattern = sequence[index..(index + pattern.len())].to_vec();
            if divide_min_hamming_distance(&current_pattern, &pattern) <= required_distance {
                location_index = index;
                break 'bigloop;
            }
        }
        required_distance += 1;
    }
    location_index
}

pub fn divide_find_frequent_kmers_with_mismatches (sequence: &Vec<u8>, k_mer: usize, required_distance: usize) -> (Vec<Vec<u8>>, Vec<Vec<usize>>) {
    // initialize stuff
    let mut frequent_patterns: Vec<Vec<u8>> = vec![];
    let number_of_bases: i32 = 4;
    let number_of_total_kmers = (number_of_bases.pow(k_mer as u32)) as usize;
    let mut frequency_array: Vec<usize> = vec![0; number_of_total_kmers];
    let mut obtained_indices_array: Vec<Vec<usize>> = vec![vec![]; number_of_total_kmers];
    let mut max_obatained_indices: Vec<Vec<usize>> = vec![];
    let mut max_frequency: usize = 0;
    
    // go through all the kmers in the sequence
    for index in 0..number_of_total_kmers {
        // get the pattern
        let pattern = divide_number_to_pattern(index, k_mer);
        //println!("pattern {:?}", pattern);
        // count the frequency of the pattern
        let (obtained_frequency, obtained_indices) = divide_approximate_pattern_count(&sequence, &pattern, required_distance);
        frequency_array[index] = obtained_frequency;
        obtained_indices_array[index] = obtained_indices;
        if obtained_frequency >= max_frequency {
            max_frequency = obtained_frequency;
            //println!("max frequency pattern {:?} {}", pattern, max_frequency);
        }
    }
    // find the max count in frequency array and get those patterns, get the positions by iterating through it.
    let max_indices = frequency_array
                            .iter()
                            .enumerate()
                            .filter_map(|(i, &r)| if r == max_frequency { Some(i) } else { None })
                            .collect::<Vec<_>>();
    for index in max_indices {
        let pattern = divide_number_to_pattern(index, k_mer);
        frequent_patterns.push(pattern.clone());
        max_obatained_indices.push(obtained_indices_array[index].clone());
    }
    (frequent_patterns, max_obatained_indices)
}

pub fn divide_approximate_pattern_count (sequence: &Vec<u8>, pattern: &Vec<u8>, required_distance: usize) -> (usize, Vec<usize>) {
    let mut k_mer_count: usize = 0;
    let mut location_indices: Vec<usize> = vec![];
    for index in 0..(sequence.len() - pattern.len() + 1) {
        let current_pattern = sequence[index..(index + pattern.len())].to_vec();
        if divide_min_hamming_distance(&current_pattern, &pattern) <= required_distance {
            k_mer_count += 1;
            location_indices.push(index);
        }
    }
    (k_mer_count, location_indices)
}

pub fn divide_number_to_pattern (number: usize, k: usize) -> Vec<u8> {
    let mut processing_number = number;
    let mut pattern: Vec<u8> = vec![];
    for _ in 0..k {
        let remainder = processing_number % 4;
        match remainder {
            0 => {pattern.push(65);},
            1 => {pattern.push(67);},
            2 => {pattern.push(71);},
            3 => {pattern.push(84);},
            _ => {},
        }
        processing_number = processing_number / 4;
    }
    pattern.reverse();
    pattern
}

pub fn divide_pattern_to_number (pattern: Vec<u8>) -> usize {
    let mut prefix_value: usize = 0;
    let symbol_value: usize;
    if pattern.len() == 0 {
        return 0
    }
    symbol_value = divide_symbol_to_number(pattern[pattern.len() - 1]);
    if pattern.len() >= 1 {
        prefix_value = 4 * divide_pattern_to_number(pattern[0..pattern.len() - 1].to_vec());
    }
    //println!(" {} {} ", symbol_value, prefix_value);
    symbol_value + prefix_value
}

pub fn divide_symbol_to_number (symbol: u8) -> usize {
    let mut number: usize = 0;
    match symbol {
        65 => {number = 0;},
        67 => {number = 1;},
        71 => {number = 2;},
        84 => {number = 3;},
        _ => {}
    }
    number
}

pub fn divide_min_hamming_distance (seq1: &Vec<u8>, seq2: &Vec<u8>) -> usize {
    let mut hamming_distance: usize = 0;
    // go through them one by one and count the differences
    for index in 0..seq1.len() {
        if seq1[index] != seq2[index] {
            hamming_distance += 1;
        }
    }
    hamming_distance
}

pub fn bfs_consensus (graph: &Graph<u8, i32, Directed, usize>, seq_vec: &Vec<String>) -> Vec<u8> {
    let mut bfs_vector:Vec<(usize, usize, bool)> = vec![]; //(node_number, depth, visited)
    // get the initial node from topology sort
    let mut topologically_ordered = Topo::new(graph);
    let head_node_index = topologically_ordered.next(graph).unwrap().index();
    // process it
    let head_node_content = (head_node_index, 0, true);
    // make the vector
    bfs_vector.push(head_node_content);
    //println!("{}", format!("{:?}", Dot::new(&graph.map(|_, n| (*n) as char, |_, e| *e))));
    // run bfs on the inital node
    bfs_vector = bfs(graph, head_node_index, 0, bfs_vector);

    // sort the bfs vector by depth
    bfs_vector.sort_by_key(|r| r.1);

    // reposition the nodes with depth
    bfs_vector = bfs_reposition(graph, bfs_vector);

    // get a single consensus
    let (bfs_consensus_long_vec, _) = bfs_get_single_consensus (graph, &bfs_vector);

    // get a single consensus which has the highest score
    let bfs_consensus_short = bfs_filter(&bfs_consensus_long_vec, seq_vec);

    // print the vector
    let mut current_depth = 100;
    for entry in bfs_vector {
        if current_depth != entry.1 {
            println!("");
            print!("depth: {}", entry.1);
        }
        print!(" {}[{}] ", graph.raw_nodes()[entry.0].weight as char, entry.0);
        current_depth = entry.1;
    }
    println!("");
    
    //println!("{:?}", bfs_consensus_short);
    bfs_consensus_short
}

pub fn bfs_filter (bfs_consensus_long_vec: &Vec<Vec<u8>>, seq_vec: &Vec<String>) -> Vec<u8> {
    let mut bfs_consensus_short: Vec<u8>;
    let mut index = 5 as i32; // get the middle consensus
    let score_up = get_consensus_score(seq_vec, &bfs_consensus_long_vec[(index + 1) as usize]);
    let score_down = get_consensus_score(seq_vec, &bfs_consensus_long_vec[(index - 1) as usize]);
    //println!("score up {} score down {}", score_up, score_down);
    let next_step;
    let mut current_score;
    let mut prev_score;
    if score_up > score_down {
        next_step = 1;
        current_score = score_up;
        bfs_consensus_short = bfs_consensus_long_vec[(index + 1) as usize].clone();
    }
    else {
        next_step = -1;
        current_score = score_down;
        bfs_consensus_short = bfs_consensus_long_vec[(index - 1) as usize].clone();
    }
    index += next_step;

    loop {
        prev_score = current_score;
        //println!("index {}", index);
        current_score = get_consensus_score(seq_vec, &bfs_consensus_long_vec[index as usize]);
        if (prev_score > current_score) || (index + next_step >= bfs_consensus_long_vec.len() as i32) || (index + next_step < 0) {
            //println!("current score: {} prev score {}", current_score, prev_score);
            break;
        }
        else {
            bfs_consensus_short = bfs_consensus_long_vec[index as usize].clone();
        }
        index += next_step;
    }
    bfs_consensus_short
}

pub fn bfs_get_single_consensus (graph: &Graph<u8, i32, Directed, usize>, bfs_vector: &Vec<(usize, usize, bool)>) -> (Vec<Vec<u8>>, Vec<Vec<usize>>) {
    let mut bfs_consensus_long: Vec<Vec<u8>> = vec![vec![]; 10];
    let mut bfs_topology_long: Vec<Vec<usize>> = vec![vec![]; 10];
    let mut nothing_available: usize = 0;
    let mut depth: usize = 0;
    loop {
        // get all the nodes of that depth
        let current_depth_nodes = bfs_vector
                                            .iter()
                                            .enumerate()
                                            .filter_map(|(_, &r)| if r.1 == depth { Some(r.0) } else { None })
                                            .collect::<Vec<_>>();
        if current_depth_nodes.len() == 0 {
            nothing_available += 1; 
        }
        else {
            nothing_available = 0;
            let mut max_number_of_seq: usize = 0;
            let mut max_node = 0;
            // find the node with highest number of seq passing through
            for depth_node in current_depth_nodes {
                let number_of_seq = find_the_seq_passing_through(depth_node, graph);
                if number_of_seq > max_number_of_seq {
                    max_number_of_seq = number_of_seq;
                    max_node = depth_node;
                }
            }
            for filter_node in 0..10 {
                if max_number_of_seq >= filter_node {
                    bfs_consensus_long[filter_node].push(graph.raw_nodes()[max_node].weight);
                    bfs_topology_long[filter_node].push(max_node); 
                }
            } 
        }
        if nothing_available > 2000 {
            break;
        }
        depth += 1;
    }
    (bfs_consensus_long, bfs_topology_long)
}

pub fn bfs_reposition (graph: &Graph<u8, i32, Directed, usize>, mut bfs_vector: Vec<(usize, usize, bool)>) -> Vec<(usize, usize, bool)> {
    for i in 0..bfs_vector.len() {
        //println!("current node: {}, depth {}", bfs_vector[i].0, bfs_vector[i].1);
        // get the current depth
        let current_depth = bfs_vector[i].1;
        // get the current depth nodes
        let mut current_depth_nodes = bfs_vector
                                                .iter()
                                                .enumerate()
                                                .filter_map(|(_, &r)| if r.1 == current_depth { Some(r.0) } else { None })
                                                .collect::<Vec<_>>();
        //println!("current depth nodes: {:?}", current_depth_nodes);
        for j in i + 1..bfs_vector.len() {
            let mut mutually_exclusive = true;
            // get the node index
            let node_index = bfs_vector[j].0;
            // check if in current depth nodes
            if current_depth_nodes.contains(&node_index) {
                continue;
            }
            // check if mutually exclusive to each other //add if it is
            for index in 0..current_depth_nodes.len() {
                let compare_node = current_depth_nodes[index];
                // get the parents
                let compare_node_parents = get_direction_nodes(Incoming, 1, vec![], compare_node, graph);
                // get the childen
                let compare_node_children = get_direction_nodes(Outgoing, 1, vec![], compare_node, graph);
                // check if node index is parent or child of the parallel nodes
                if compare_node_parents.contains(&node_index) || compare_node_children.contains(&node_index) {
                    mutually_exclusive = false;
                    break;
                }
            }
            // if not mutually exclusive break
            if !mutually_exclusive {
                break;
            }
            else {
                //println!("added node: {}", node_index);
                current_depth_nodes.push(node_index);
            }
        }
        // modify the bfs vector with current_depth_nodes
        // get the indices of bfs vector
        let bfs_vector_modify_indices = bfs_vector
                                                .iter()
                                                .enumerate()
                                                .filter_map(|(i, &r)| if current_depth_nodes.contains(&r.0) { Some(i) } else { None })
                                                .collect::<Vec<_>>();
        // modify them
        for index in bfs_vector_modify_indices {
            //println!("processing node {}", bfs_vector[index].0);
            bfs_vector[index].1 = current_depth;
        }
    }
    bfs_vector
}

pub fn bfs (graph: &Graph<u8, i32, Directed, usize>, head_index: usize, head_depth: usize, mut bfs_vector: Vec<(usize, usize, bool)>) -> Vec<(usize, usize, bool)> {
    // queue for saving stuff
    let mut test_index = 0;
    let mut bfs_queue: Queue<(usize, usize)> = queue![]; //node id is saved in the queue
    bfs_queue.add((head_index, head_depth)).unwrap();
    'test: while bfs_queue.size() > 0 {
        let (node_index, node_depth) = bfs_queue.remove().unwrap();
        // println!("bfs on {}", node_index);
        // sort the bfs vector by depth

        // print the vector
        //let mut current_depth = 100;
        /*for entry in &bfs_vector {
            if current_depth != entry.1 {
                println!("");
                print!("depth: {}", entry.1);
            }
            print!(" {}[{}] ", graph.raw_nodes()[entry.0].weight as char, entry.0);
            current_depth = entry.1;
        }
        println!("");*/

         // get the children of the node
        let children = get_direction_nodes(Outgoing, 1, vec![], node_index, graph);
        //println!("bfs children {:?}", children);
        // process all childeren first
        for child in &children {
            match bfs_vector.iter().position(|r| r.0 == *child) {
                Some(x) => {
                    // this node is already visited
                    // increase depth of all the nodes in this depth except the parent of this node
                    if bfs_vector[x].1 < node_depth + 1 {
                        //println!("increase node{}", bfs_vector[x].0 );
                        let parent_bfs_index = bfs_vector.iter().position(|r| r.0 == node_index).unwrap();
                        let parent_bfs_depth = bfs_vector[parent_bfs_index].1;
                        bfs_vector = bfs_depth_increase(graph, node_index, bfs_vector, *child, parent_bfs_depth + 1);
                        if test_index > 100 {
                            break 'test;
                        }
                        test_index += 1;
                    }
                },
                None => {
                    let parent_bfs_index = bfs_vector.iter().position(|r| r.0 == node_index).unwrap();
                    let parent_bfs_depth = bfs_vector[parent_bfs_index].1;
                    bfs_queue.add((*child, parent_bfs_depth + 1)).unwrap();
                    bfs_vector.push((*child, parent_bfs_depth + 1, true));
                }
            }
            
        }
    }
    /* 
    // run bfs on children
    for child in &children {
        if !avoid_children.contains(child) {
            //get parents depth before doing this as it could have been changed
            let parent_bfs_index = bfs_vector.iter().position(|r| r.0 == node_index).unwrap();
            let parent_bfs_depth = bfs_vector[parent_bfs_index].1;
            bfs_vector = bfs(graph, *child, parent_bfs_depth + 1, bfs_vector);
        }
    }
    */
    bfs_vector
}

pub fn bfs_depth_increase (graph: &Graph<u8, i32, Directed, usize>, node_index: usize, mut bfs_vector: Vec<(usize, usize, bool)>, node_to_increase: usize, depth_to_increase: usize) -> Vec<(usize, usize, bool)> {
    // find the current depth of node_to_increase
    //println!("I am node {:?}", node_index);
    //println!("original node to increase {}", node_to_increase);
    let node_to_increase_index = bfs_vector.iter().position(|r| r.0 == node_to_increase).unwrap();
    let node_to_increase_depth = bfs_vector[node_to_increase_index].1;
    
    // find all the nodes with that depth
    let nodes_to_increase_indices = bfs_vector
                    .iter()
                    .enumerate()
                    .filter_map(|(_, &r)| if r.1 == node_to_increase_depth { Some(r.0) } else { None })
                    .collect::<Vec<_>>();
    
    //println!("nodes to increase {:?}", nodes_to_increase_indices);
    // find the node which is ancestor of node_index who has the node_to_increase_depth
    let mut search_radius = 1;
    let mut ancestor_index = 0;
    let mut ancestors = vec![node_index];
    loop {
        //println!("ancestor searching {:?}", ancestors);
        match bfs_vector.iter().position(|r| (ancestors.contains(&r.0) && r.1 == node_to_increase_depth)) {
            Some(x) => {
                ancestor_index = bfs_vector[x].0;
                break;
            },
            None => {
                search_radius += 1;
            }
        }
        if search_radius > 10 {
            break;
        }
        ancestors = get_xiterations_direction_nodes(Incoming, search_radius, vec![], node_index, graph);
    }
    //println!("ancestor node {}", ancestor_index);
    // update the depth of all nodes in that depth and their children except the ancestor nodes
    for iter_index in nodes_to_increase_indices {
        if iter_index != ancestor_index {
            // increase the depth of the index
            // find the position of the index in bfs thing
            let iter_index_pos = bfs_vector.iter().position(|r| (r.0 == iter_index)).unwrap();
            //println!("increasing parent node {} from {} to {}", bfs_vector[iter_index_pos].0, bfs_vector[iter_index_pos].1, depth_to_increase);
            bfs_vector[iter_index_pos].1 = depth_to_increase;

            let mut search_radius = 1;
            let mut none_in_this_iteration;
            loop {
                none_in_this_iteration = true;
                let decendants = get_xiterations_direction_nodes(Outgoing, search_radius, vec![], iter_index, graph);
                // find the decendants in the bfs vector and modify them
                for decendant in decendants {
                    match bfs_vector.iter().position(|r| (r.0 == decendant)) {
                        Some(x) => {
                            none_in_this_iteration = false;
                            //println!("increasing child node {} from {} to {}", bfs_vector[x].0, bfs_vector[x].1, depth_to_increase + search_radius);
                            bfs_vector[x].1 = depth_to_increase + search_radius;
                            
                        },
                        None => {

                        }
                    }
                }
                //if none in this iter break
                if none_in_this_iteration {
                    break;
                }
                if search_radius > 10 {
                    break;
                }
                search_radius += 1;
                
            }
        }
    }
    bfs_vector
}

pub fn topology_cut_consensus (seqvec: &Vec<String>) -> (Vec<u8>, Vec<usize>) {
    let mut output: Vec<u8> = vec![];
    let mut topopos: Vec<usize> = vec![];
    let mut topologically_ordered_indices = vec![];
    
    //get the poa graph
    let scoring = Scoring::new(GAP_OPEN, GAP_EXTEND, |a: u8, b: u8| if a == b { MATCH } else { MISMATCH });
    let mut seqnum: u8 = 0;
    let mut aligner = Aligner::new(scoring, seqvec[0].as_bytes());
    for seq in seqvec{
        if seqnum != 0 {
            aligner.global(seq.as_bytes()).add_to_graph();
        }
        seqnum += 1;
        println!("Sequence {} processed", seqnum);
    }
    let graph = aligner.graph();
    //get the topologically ordered vector
    let mut topologically_ordered = Topo::new(graph);

    while let Some(node) = topologically_ordered.next(graph) {
        topologically_ordered_indices.push(node.index());
    }
    //for each node check the parallel counts and save the parallel nodes
    let mut node_neighbour_counts_parallel: Vec<(usize, usize, (Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>), (Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>))> = vec![];
    for topologically_ordered_index in topologically_ordered_indices {
        // variables for fun
        let mut acgt_count = (vec![], vec![], vec![], vec![]);
        let mut acgt_nodes = (vec![], vec![], vec![], vec![]);

        // get largest incoming as parent
        let mut max_parent_weight = 0;
        let mut chosen_parent: Option<usize> = None;
        let mut parent_candidates = graph.neighbors_directed(NodeIndex::new(topologically_ordered_index), Incoming);
        while let Some(parent_candidate) = parent_candidates.next() {
            match graph.find_edge(parent_candidate, NodeIndex::new(topologically_ordered_index)) {
                Some(edge) => {
                    if max_parent_weight < *graph.edge_weight(edge).unwrap() {
                        max_parent_weight = *graph.edge_weight(edge).unwrap();
                        chosen_parent = Some(parent_candidate.index());
                    }
                }
                None => {},
            }
        }
        // get largest outgoing as child
        let mut max_child_weight = 0;
        let mut chosen_child: Option<usize> = None;
        let mut child_candidates = graph.neighbors_directed(NodeIndex::new(topologically_ordered_index), Outgoing);
        while let Some(child_candidate) = child_candidates.next() {
            match graph.find_edge(NodeIndex::new(topologically_ordered_index), child_candidate) {
                Some(edge) => {
                    if max_child_weight < *graph.edge_weight(edge).unwrap() {
                        max_child_weight = *graph.edge_weight(edge).unwrap();
                        chosen_child = Some(child_candidate.index());
                    }
                }
                None => {},
            }
        }
        // create skip nodes from parent and child
        let mut skip_nodes = vec![topologically_ordered_index];
        match chosen_child {
            Some(x) => {skip_nodes.push(x)},
            None => {},
        }
        match chosen_parent {
            Some(x) => {skip_nodes.push(x)},
            None => {},
        }
        
        // get the parallel nodes
        let (temp_parallel_nodes, temp_parallel_num_incoming_seq, _) = get_parallel_nodes_with_topology_cut(skip_nodes, seqnum as usize, topologically_ordered_index, chosen_parent, chosen_child, graph);
        let mut parallel_nodes = vec![];
        let mut parallel_num_incoming_seq = vec![];

        // remove recurring 
        for index in 0..temp_parallel_nodes.len() {
            if !parallel_nodes.contains(&temp_parallel_nodes[index]) {
                parallel_nodes.push(temp_parallel_nodes[index].clone());
                parallel_num_incoming_seq.push(temp_parallel_num_incoming_seq[index].clone());
            }
        }
        //println!("parallel nodes: {:?}", parallel_nodes);
        for index in 0..parallel_nodes.len() {
            match graph.raw_nodes()[parallel_nodes[index]].weight {
                65 => {
                    acgt_count.0.push(parallel_num_incoming_seq[index]);
                    acgt_nodes.0.push(parallel_nodes[index]);
                },
                67 => {
                    acgt_count.1.push(parallel_num_incoming_seq[index]);
                    acgt_nodes.1.push(parallel_nodes[index]);
                },
                71 => {
                    acgt_count.2.push(parallel_num_incoming_seq[index]);
                    acgt_nodes.2.push(parallel_nodes[index]);
                },
                84 => {
                    acgt_count.3.push(parallel_num_incoming_seq[index]);
                    acgt_nodes.3.push(parallel_nodes[index]);
                },
                _ => {}
            }
        }
        node_neighbour_counts_parallel.push((topologically_ordered_index, 12345678, acgt_count, acgt_nodes));
    }
    //for every node in topological list, get the outgoing nodes, filter ones with the corrosponding base, if multiple select highest edge one
    for entry in node_neighbour_counts_parallel.clone() {
        // get all the outgoing neighbours
        let mut neighbours = graph.neighbors_directed(NodeIndex::new(entry.0), Outgoing);
        // find which one is the max of summed up ACGT counts
        let mut acgt_count = [0, 0, 0, 0];
        let mut acgt_nodes = [vec![], vec![], vec![], vec![]];
        while let Some(neighbour) = neighbours.next() {
            // find the location of the neightbour node in the array
            let position_in_node_neighbour_counts_parallel_array = node_neighbour_counts_parallel.iter().position(|r| r.0 == neighbour.index()).unwrap();
            // update count
            for node_index in 0..node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.0.len() {
                //if !acgt_nodes[0].contains(&node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.0[node_index]) {
                    acgt_nodes[0].push(node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.0[node_index]);
                    acgt_count[0] += node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].2.0[node_index];
                //}
            }
            for node_index in 0..node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.1.len() {
                //if !acgt_nodes[1].contains(&node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.1[node_index]) {
                    acgt_nodes[1].push(node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.1[node_index]);
                    acgt_count[1] += node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].2.1[node_index];
                //}
            }
            for node_index in 0..node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.2.len() {
                //if !acgt_nodes[2].contains(&node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.2[node_index]) {
                    acgt_nodes[2].push(node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.2[node_index]);
                    acgt_count[2] += node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].2.2[node_index];
                //}
            }
            for node_index in 0..node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.3.len() {
                //if !acgt_nodes[3].contains(&node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.3[node_index]) {
                    acgt_nodes[3].push(node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].3.3[node_index]);
                    acgt_count[3] += node_neighbour_counts_parallel[position_in_node_neighbour_counts_parallel_array].2.3[node_index];
                //}   
            }
        }
        let mut max_base = 0;
        let mut max_node = 0;
        // find the max base (from ACGT)
        if acgt_count[0] >= cmp::max(cmp::max(acgt_count[1], acgt_count[2]), acgt_count[3]) {
            max_base = 0;
        }
        else if acgt_count[1] >= cmp::max(cmp::max(acgt_count[0], acgt_count[2]), acgt_count[3]) {
            max_base = 1;
        }
        else if acgt_count[2] >= cmp::max(cmp::max(acgt_count[1], acgt_count[0]), acgt_count[3]) {
            max_base = 2;
        }
        else if acgt_count[3] >= cmp::max(cmp::max(acgt_count[1], acgt_count[2]), acgt_count[0]) {
            max_base = 3;
        }
        
        // filter out the other nodes
        let mut important_nodes = [acgt_nodes[0].clone(), acgt_nodes[1].clone(), acgt_nodes[2].clone(), acgt_nodes[3].clone()].concat().to_vec();
        // get the max outgoing node 
        // get largest outgoing as child
        let mut max_weight = 0;
        let mut chosen_node = 0;
        for important_node in &important_nodes {
            match graph.find_edge(NodeIndex::new(entry.0), NodeIndex::new(*important_node)) {
                Some(edge) => {
                    if max_weight < *graph.edge_weight(edge).unwrap() {
                        max_weight = *graph.edge_weight(edge).unwrap();
                        chosen_node = *important_node;
                    }
                }
                None => {},
            }
        }
        let mut chosen_base = 0;
        if acgt_nodes[0].contains(&chosen_node) {
            chosen_base = 0;
        }
        if acgt_nodes[1].contains(&chosen_node) {
            chosen_base = 1;
        }
        if acgt_nodes[2].contains(&chosen_node) {
            chosen_base = 2;
        }
        if acgt_nodes[3].contains(&chosen_node) {
            chosen_base = 3;
        }
        match max_base {
            0 => {important_nodes = acgt_nodes[0].clone()},
            1 => {important_nodes = acgt_nodes[1].clone()},
            2 => {important_nodes = acgt_nodes[2].clone()},
            3 => {important_nodes = acgt_nodes[3].clone()},
            _ => {},
        }
        // among the filtered nodes select highest weighted edge
        let mut highest_weight_index = (0, 12345678);
        for important_node in &important_nodes {
            // get the weight from node to important node
            match graph.find_edge(NodeIndex::new(entry.0), NodeIndex::new(*important_node)) {
                Some(edge) => {
                    if highest_weight_index.0 < *graph.edge_weight(edge).unwrap() {
                        highest_weight_index.0 = *graph.edge_weight(edge).unwrap();
                        highest_weight_index.1 = *important_node;
                        max_node = *important_node;
                    }
                }
                None => {},
            }
        }
        
        if chosen_base != max_base {
            // check if the chosen_ node is a predecessor of max_node if so ignore
            let back_nodes = get_direction_nodes(Incoming, 3, vec![], max_node, graph);
            if back_nodes.contains(&chosen_node) {
                highest_weight_index.1 = chosen_node;
            }
        }
        //double checking 
        if PRINT_ALL {
            println!("I am {}", entry.0);
            println!("chose base (normal) {} weight {} node {}", chosen_base, max_weight, chosen_node);
            println!("max base (topo) {} {}", max_base, max_node);
            println!("count {:?}", acgt_count);
            println!("nodes {:?}", acgt_nodes);
            println!("");
        }
        // save the passing node
        let position = node_neighbour_counts_parallel.iter().position(|r| r.0 == entry.0).unwrap();
        node_neighbour_counts_parallel[position].1 = highest_weight_index.1;

        

    }
    for entry in &node_neighbour_counts_parallel {
        println!("{} -> {}", entry.0, entry.1);
    } 
    //start with topologically sorted top node
    let mut current_node = node_neighbour_counts_parallel[0].clone();
    loop {
        output.push(graph.raw_nodes()[current_node.0].weight);
        topopos.push(current_node.0);
        if current_node.1 != 12345678 {
            // iterate to the next node
            let position = node_neighbour_counts_parallel.iter().position(|r| r.0 == current_node.1).unwrap();
            current_node = node_neighbour_counts_parallel[position].clone();
        }
        else {
            //find a node to traverse, if none available break
            let mut highest_weight_index = (0, 12345678);
            let mut neighbours = graph.neighbors_directed(NodeIndex::new(current_node.0), Outgoing);
            while let Some(neighbour) = neighbours.next() {
                match graph.find_edge(NodeIndex::new(current_node.0), neighbour) {
                    Some(edge) => {
                        if highest_weight_index.0 < *graph.edge_weight(edge).unwrap() {
                            highest_weight_index.0 = *graph.edge_weight(edge).unwrap();
                            highest_weight_index.1 = neighbour.index();
                        }
                    }
                    None => {},
                }
            }
            match node_neighbour_counts_parallel.iter().position(|r| r.0 == highest_weight_index.1) {
                Some(x) => {current_node = node_neighbour_counts_parallel[x].clone();},
                None => {break;}
            };
        }
    }
    //println!("{}", format!("{:?}", Dot::new(&graph.map(|_, n| (*n) as char, |_, e| *e))));
    (output, topopos)
}

fn heavy_bundle_modified_consensus (seqvec: &Vec<String>) -> (Vec<u8>, Vec<usize>) {
    //get the normal consensus and graph
    let scoring = Scoring::new(GAP_OPEN, GAP_EXTEND, |a: u8, b: u8| if a == b { MATCH } else { MISMATCH });
    let mut seqnum: u8 = 0;
    let mut aligner = Aligner::new(scoring, seqvec[0].as_bytes());
    for seq in seqvec{
        if seqnum != 0 {
            aligner.global(seq.as_bytes()).add_to_graph();
        }
        seqnum += 1;
        println!("Sequence {} processed", seqnum);
    }
    let consensus;
    let topology;
    (consensus, topology) = aligner.poa.consensus(); //just poa
    let graph = aligner.graph();
    let mut nodes_to_change_and_by_what: Vec<(usize, usize)> = vec![];
    let mut changed_stuff = false;
    //run all the consensus through get indices
    for i in 0..consensus.len() {
        // skip the indices which are in the passed consensus
        let skip_nodes: Vec<usize> = topology[0 .. i + 1].to_vec();
        // new method using topology cut
        let mut target_node_parent = None;
        let mut target_node_child = None;
        if i != 0 {
            target_node_parent = Some(topology[i - 1]);
        }
        if i != consensus.len() - 1 {
            target_node_child = Some(topology[i + 1]);
        }
        let (parallel_nodes, parallel_num_incoming_seq, _) = get_parallel_nodes_with_topology_cut (skip_nodes, seqvec.len(),  topology[i], target_node_parent, target_node_child, graph);
        //println!("base: {} parallel nodes {:?} count {:?}", consensus[i] as char, parallel_nodes, parallel_num_incoming_seq);
        // check the parallel nodes bases
        // number of As number of Cs number of Gs number of Ts
        let mut acgt_count = [0, 0, 0, 0];
        let mut acgt_nodes = [vec![], vec![], vec![], vec![]];
        let mut target_base_index = 0;
        match consensus[i] {
            65 => {target_base_index = 0;},
            67 => {target_base_index = 1;},
            71 => {target_base_index = 2;},
            84 => {target_base_index = 3;},
            _ => {}
        }
        for index in 0..parallel_nodes.len() {
            match graph.raw_nodes()[parallel_nodes[index]].weight {
                65 => {
                    acgt_count[0] += parallel_num_incoming_seq[index];
                    acgt_nodes[0].push(parallel_nodes[index]);
                },
                67 => {
                    acgt_count[1] += parallel_num_incoming_seq[index];
                    acgt_nodes[1].push(parallel_nodes[index]);
                },
                71 => {
                    acgt_count[2] += parallel_num_incoming_seq[index];
                    acgt_nodes[2].push(parallel_nodes[index]);
                },
                84 => {
                    acgt_count[3] += parallel_num_incoming_seq[index];
                    acgt_nodes[3].push(parallel_nodes[index]);
                },
                _ => {}
            }
        }
        /*if acgt_count[target_base_index] < acgt_count[0] {
            consensus[i] = 65;
            changed_stuff = true;
        }
        else if acgt_count[target_base_index] < acgt_count[1] {
            consensus[i] = 67;
            changed_stuff = true;
        }
        else if acgt_count[target_base_index] < acgt_count[2] {
            consensus[i] = 71;
            changed_stuff = true;
        }
        else if acgt_count[target_base_index] < acgt_count[3] {
            consensus[i] = 84;
            changed_stuff = true;
        }*/
    
        // determine if change is required and to what nodes and to what number and save them
        if acgt_count[target_base_index] < acgt_count[0] {
            changed_stuff = true;
            for node in &acgt_nodes[0] {
                nodes_to_change_and_by_what.push((*node, acgt_count[0]));
            }
        }
        else if acgt_count[target_base_index] < acgt_count[1] {
            changed_stuff = true;
            for node in &acgt_nodes[1] {
                nodes_to_change_and_by_what.push((*node, acgt_count[1]));
            }
        }
        else if acgt_count[target_base_index] < acgt_count[2] {
            changed_stuff = true;
            for node in &acgt_nodes[2] {
                nodes_to_change_and_by_what.push((*node, acgt_count[2]));
            }
        }
        else if acgt_count[target_base_index] < acgt_count[3] {
            changed_stuff = true;
            for node in &acgt_nodes[3] {
                nodes_to_change_and_by_what.push((*node, acgt_count[3]));
            }
        }
    }
    // change the graph
    println!("CHANGED STUFF {} {:?}", changed_stuff, nodes_to_change_and_by_what);
    let mut node_neighbour_values = vec![];
    for (node, value) in nodes_to_change_and_by_what {
        // find the outgoing edges
        let mut neighbours = graph.neighbors_directed(NodeIndex::new(node), Outgoing);
        let mut max_weight = 0;
        // find the max weight of the outgoing edges
        while let Some(neighbour) = neighbours.next() {
            match graph.find_edge(NodeIndex::new(node), neighbour) {
                Some(edge) => {
                    if max_weight <= *graph.edge_weight(edge).unwrap() {
                        max_weight = *graph.edge_weight(edge).unwrap();
                        node_neighbour_values.push((node, neighbour.index(), value));
                    }
                }
                None => {},
            }
        }
    }
    // increase the weights
    for node_neighbour_value in node_neighbour_values {
        aligner.poa.change_edge_weight(node_neighbour_value.0, node_neighbour_value.1, node_neighbour_value.2 as i32);
    }
    // get the consensus again and return it
    let (consensus, topology) = aligner.poa.consensus();
    (consensus, topology)
}

fn create_required_result_files (path: &str) -> (String, String, String, String, String, String, String) {
    let output_debug_file_name: String = [path, "/debug.txt"].concat();
    let output_consensus_file_name: String = [path, "/consensus.fa"].concat();
    let output_scores_file_name: String = [path, "/results.txt"].concat();
    let output_normal_graph_file_name: String = [path, "/normal_graph.fa"].concat();
    let output_homopolymer_graph_file_name: String = [path, "/homopolymer_graph.fa"].concat();
    let output_quality_graph_file_name: String = [path, "/quality_graphs.fa"].concat();
    let output_quality_file_name: String = [path, "/quality_scores.fa"].concat();
    File::create([path, "/consensus.fa"].concat()).ok();
    File::create([path, "/filtered_data.fa"].concat()).ok();
    File::create([path, "/results.txt"].concat()).ok();
    File::create([path, "/normal_graph.fa"].concat()).ok();
    File::create([path, "/homopolymer_graph.fa"].concat()).ok();
    File::create([path, "/quality_scores.fa"].concat()).ok();
    File::create([path, "/quality_graphs.fa"].concat()).ok();
    File::create([path, "/debug.txt"].concat()).ok();
    (output_debug_file_name, output_consensus_file_name, output_scores_file_name, output_normal_graph_file_name, output_homopolymer_graph_file_name, output_quality_graph_file_name, output_quality_file_name)
}

fn check_the_scores_and_change_alignment (seqvec: Vec<String>, pacbio_consensus: &String) -> Vec<String> {
    let mut invert: bool = false;
    let mut seqvec2: Vec<String> = vec![];
    // check the scores for 3 sequences
    let mut index = 0;
    for seq in &seqvec{
        let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
        let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(pacbio_consensus.len(), seq.len(), GAP_OPEN, GAP_EXTEND, &score);
        let alignment = aligner.global(&pacbio_consensus.as_bytes(), &seq.as_bytes());
        if alignment.score < 1000 {
            invert = true;
            break;
        }
        else if index > 3 {
            break;
        }
        index += 1;
    }
    if invert {
        println!("Scores are too low, inverting sequences.");
        //reverse complement every line
        for seq in &seqvec {
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
            seqvec2.push(tempseq.iter().cloned().collect::<String>());
        }
    }
    else {
        seqvec2 = seqvec;
    }
    seqvec2
}

fn write_debug_data_to_file (filename: impl AsRef<Path>, debug_strings: Vec<String>) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n DEBUG",
        chrono::offset::Local::now(), INPUT_FILE_NAME)
        .expect("result file cannot be written");
    for line in debug_strings {
        writeln!(file,
            "{}",
            line)
            .expect("result file cannot be written");
    }
}

fn get_quality_score_aligned (pacbio_consensus: String, calculated_consensus: &Vec<u8>, pacbio_quality_scores: String, error_line_number: usize) -> (Vec<usize>, Vec<usize>, Vec<u8>, usize) {
    let mut consensus_match_invalid_indices: Vec<usize> = vec![];
    let pacbio_consensus_vec: Vec<u8> = pacbio_consensus.bytes().collect();
    let pacbio_quality_scores_vec: Vec<char> =  pacbio_quality_scores.chars().collect();
    let mut aligned_pacbio_scores_vec: Vec<usize> = vec![];
    let mut aligned_pacbio_bases:Vec<u8> = vec![];
    let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
    let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(calculated_consensus.len(), pacbio_consensus_vec.len(), GAP_OPEN, GAP_EXTEND, &score);
    let alignment = aligner.global(&calculated_consensus, &pacbio_consensus_vec);
    let mut calc_index = alignment.xstart;
    let mut pacbio_index = alignment.ystart;
    let mut calc_error_line_number: usize = 0;
    for op in &alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                aligned_pacbio_scores_vec.push(pacbio_quality_scores_vec[pacbio_index] as usize);
                aligned_pacbio_bases.push(pacbio_consensus_vec[pacbio_index]);
                pacbio_index += 1;
                calc_index += 1;
            },
            bio::alignment::AlignmentOperation::Subst => {
                aligned_pacbio_scores_vec.push(pacbio_quality_scores_vec[pacbio_index] as usize);
                aligned_pacbio_bases.push(pacbio_consensus_vec[pacbio_index]);
                consensus_match_invalid_indices.push(calc_index);
                pacbio_index += 1;
                calc_index += 1;
            },
            bio::alignment::AlignmentOperation::Del => {
                pacbio_index += 1;
            },
            bio::alignment::AlignmentOperation::Ins => {
                consensus_match_invalid_indices.push(calc_index);
                aligned_pacbio_bases.push(126);
                aligned_pacbio_scores_vec.push(33);
                calc_index += 1;
            },
            _ => {},
        }
        if error_line_number == pacbio_index {
            calc_error_line_number = calc_index
        }
    }
    println!("Error line number {} corrosponds to calculated number {}", error_line_number, calc_error_line_number);
    (aligned_pacbio_scores_vec, consensus_match_invalid_indices, aligned_pacbio_bases, calc_error_line_number)
}
fn get_consensus_quality_scores(seq_num: usize, consensus: &Vec<u8>, topology: &Vec<usize>, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<f64>, Vec<bool>, Vec<Vec<usize>>, Vec<String>) {
    let mut quality_scores: Vec<f64> = vec![];
    let mut validity: Vec<bool> = vec![];
    let mut base_count_vec: Vec<Vec<usize>> = vec![];
    let mut debug_strings: Vec<String> = vec![];
    let mut temp_string: String;
    //run all the consensus through get indices
    for i in 0..consensus.len() {
        // skip the indices which are in the passed consensus
        let skip_nodes: Vec<usize> = topology[0 .. i + 1].to_vec();
        // new method using topology cut
        let mut target_node_parent = None;
        let mut target_node_child = None;
        if i != 0{
            target_node_parent = Some(topology[i - 1]);
        }
        if i != consensus.len() - 1 {
            target_node_child = Some(topology[i + 1]);
        }
        temp_string = format!("BASE NODE: {} ({})", consensus[i] as char, topology[i]);
        println!("{}", temp_string);
        debug_strings.push(temp_string.clone());
        let (parallel_nodes, parallel_num_incoming_seq, temp_debug_strings) = get_parallel_nodes_with_topology_cut (skip_nodes, seq_num,  topology[i], target_node_parent, target_node_child, graph);
        debug_strings = [debug_strings, temp_debug_strings].concat();
        let (temp_quality_score, temp_count_mismatch, temp_base_counts, temp_debug_strings) = base_quality_score_calculation (seq_num, parallel_nodes, parallel_num_incoming_seq, consensus[i], graph);
        debug_strings = [debug_strings, temp_debug_strings].concat();
        quality_scores.push(temp_quality_score);
        validity.push(temp_count_mismatch);
        base_count_vec.push(temp_base_counts);
    }
    (quality_scores, validity, base_count_vec, debug_strings)
}

fn get_parallel_nodes_with_topology_cut (skip_nodes: Vec<usize>, total_seq: usize, target_node: usize, target_node_parent: Option<usize>, target_node_child: Option<usize>, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>, Vec<String>) {
    // vector initialization
    let mut debug_strings: Vec<String> = vec![];
    let mut topology = Topo::new(graph);
    let mut topologically_ordered_nodes = Vec::new();
    let mut parallel_nodes: Vec<usize> = vec![];
    let mut parallel_node_parents: Vec<usize> = vec![];
    let mut parallel_num_incoming_seq: Vec<usize> = vec![];
    let mut direction: Option<Direction> = None;
    let temp_string;
    //print stuff
    if PRINT_ALL {
        println!("NODE CHECKING FOR PARALLEL: {}, base {}", target_node, graph.raw_nodes()[target_node].weight as char);
    }

    // make a topologically ordered list
    while let Some(node) = topology.next(graph) {
        topologically_ordered_nodes.push(node.index());
    }
    // find the position of the target node, its child and parent in topology list
    let target_node_topological_position = topologically_ordered_nodes.iter().position(|&r| r == target_node).unwrap();
    let target_child_topological_position = match target_node_child { 
        Some(child_node) => {topologically_ordered_nodes.iter().position(|&r| r == child_node).unwrap()},
        None => {direction = Some(Incoming); target_node_topological_position}
    };
    let target_parent_topological_position = match target_node_parent { 
        Some(parent_node) => {topologically_ordered_nodes.iter().position(|&r| r == parent_node).unwrap()},
        None => {direction = Some(Outgoing); target_node_topological_position}
    };
    // choose a direction with the least amount of intermediate nodes
    if (direction == None) && (topologically_ordered_nodes[target_parent_topological_position..target_node_topological_position].len() > topologically_ordered_nodes[target_node_topological_position..target_child_topological_position].len()) {
        direction = Some(Outgoing);
    }
    else if direction == None {
        direction = Some(Incoming);
    }
    match direction {
        Some(x) => {
            if x == Incoming {
                temp_string = format!("Going backwards");
                if PRINT_ALL {
                    println!("{}", temp_string);
                }
                debug_strings.push(temp_string.clone());
            }
            else {
                temp_string = format!("Going forward");
                if PRINT_ALL {
                    println!("{}", temp_string);
                }
                debug_strings.push(temp_string.clone());
            }
        }
        None => {}
    }
    // check if the target node corrosponds with all the sequences
    let num_seq_through_target_base = find_the_seq_passing_through (target_node, graph);

    if num_seq_through_target_base == total_seq {
        parallel_nodes.push(target_node);
        if USEPACBIODATA {
            parallel_num_incoming_seq.push(num_seq_through_target_base - 1);
        }
        parallel_num_incoming_seq.push(num_seq_through_target_base);
        return (parallel_nodes, parallel_num_incoming_seq, debug_strings);
    }
    // go back skip_count and go forward skip_count + 3 and check if parent and child are before and after target_node_position,
    // iterate skip_count until all sequences are found, break on 5
    let mut seq_found_so_far = num_seq_through_target_base;
    let mut bubble_size = 1;
    while seq_found_so_far < total_seq  && bubble_size < NUM_OF_ITER_FOR_PARALLEL {
        let temp_debug_strings;
        (parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, temp_debug_strings) = move_in_direction_and_find_crossing_nodes (&skip_nodes, total_seq, direction.unwrap(), parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, target_node, bubble_size, &topologically_ordered_nodes, target_node_topological_position, graph);
        debug_strings = [debug_strings, temp_debug_strings].concat();
        bubble_size += 1;
    }
    if USEPACBIODATA {
        parallel_num_incoming_seq.push(num_seq_through_target_base - 1);
    }
    parallel_num_incoming_seq.push(num_seq_through_target_base);
    parallel_nodes.push(target_node);
    (parallel_nodes, parallel_num_incoming_seq, debug_strings)
}

fn move_in_direction_and_find_crossing_nodes (skip_nodes: &Vec<usize>, total_seq: usize, direction: Direction, mut parallel_nodes: Vec<usize>, mut parallel_node_parents: Vec<usize>, mut parallel_num_incoming_seq: Vec<usize>, mut seq_found_so_far: usize, focus_node: usize, bubble_size: usize, topologically_ordered_nodes: &Vec<usize>, target_node_position: usize, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>, Vec<usize>, usize, Vec<String>) {
    let mut debug_strings: Vec<String> = vec![];
    let mut temp_string: String;
    // get a list of x back_iterations back nodes
    let back_nodes_list = get_xiterations_direction_nodes(direction, bubble_size, vec![], focus_node, graph);
    // get a list of all forward nodes 0..(back_iterations + 3) for all the back_nodes
    let mut edge_nodes_list: Vec<usize> = vec![];
    for back_node in &back_nodes_list {
        let temp_forward_list = get_direction_nodes (direction.opposite(), bubble_size + 3, vec![], *back_node, graph);
        for temp_forward_node in &temp_forward_list {
            if !edge_nodes_list.contains(temp_forward_node) {
                edge_nodes_list.push(*temp_forward_node);
            }
        }
    }
    temp_string = format!("Iteration: {} BackNodes: {:?} CheckNodes: {:?}", bubble_size, back_nodes_list, edge_nodes_list);
    if PRINT_ALL {
        println!("{}", temp_string);
    }
    debug_strings.push(temp_string.clone());
    
    // get the two slices of topologically_ordered_list back front
    let mut slice: Vec<Vec<usize>> = [topologically_ordered_nodes[0..target_node_position].to_vec(), topologically_ordered_nodes[target_node_position + 1..topologically_ordered_nodes.len()].to_vec()].to_vec();
    // for debugging
    if slice[0].len() > 10 {
        temp_string = format!("Back slice {:?}", slice[0][(slice[0].len() - 10)..slice[0].len()].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }
    else {
        temp_string = format!("Back slice {:?}", slice[0][0..slice[0].len()].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }
    if slice[1].len() > 10 {
        temp_string = format!("Front slice {:?}", slice[1][0..10].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }
    else {
        temp_string = format!("Front slice {:?}", slice[1][0..slice[1].len()].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }

    if direction == Outgoing {
        slice.reverse();
    }
    //iterate through edge nodes obtained
    for edge_node in &edge_nodes_list {
        // get the parents of the edge node
        let edge_node_parents = get_direction_nodes (direction, 1, vec![], *edge_node, graph);
        'parent_loop: for edge_node_parent in &edge_node_parents {
            // if the parent is in back section and node is in front section add to parallel nodes or if both parent and target is in intermediate add to parallel loop
            if slice[0].contains(edge_node_parent) && slice[1].contains(edge_node) && (*edge_node_parent != focus_node) {
                // edge node parent check
                if parallel_nodes.contains(edge_node) && parallel_node_parents.contains(edge_node_parent) {
                    // go through the parallel nodes and if there is a match check if the same parent and continue if so
                    for index in 0..parallel_nodes.len() {
                        if (parallel_nodes[index] == *edge_node) && (parallel_node_parents[index] == *edge_node_parent) {
                            continue 'parent_loop;
                        }
                    }
                }
                // target node front of parallel node check
                if direction == Incoming {
                    if get_direction_nodes(Outgoing, 4, vec![],  *edge_node, graph).contains(&focus_node) {
                        continue;
                    }
                    if skip_nodes.contains(edge_node) {
                        continue;
                    }
                }
                else {
                    if get_direction_nodes(Outgoing, 4, vec![], focus_node, graph).contains(&edge_node) {
                        continue;
                    }
                    if skip_nodes.contains(edge_node_parent) {
                        continue;
                    }
                }
                // all found 
                if seq_found_so_far >= total_seq {
                    break;
                }
                parallel_nodes.push(*edge_node);
                parallel_node_parents.push(*edge_node_parent);
                temp_string = format!("success node {} parent/child {}\n", *edge_node, *edge_node_parent);
                if PRINT_ALL {
                    println!("{}", temp_string);
                }
                debug_strings.push(temp_string.clone());
                // get the edge weight and add to seq_found_so_far
                let mut incoming_weight = 0;
                if direction == Incoming {
                    let mut edges = graph.edges_connecting(NodeIndex::new(*edge_node_parent), NodeIndex::new(*edge_node));
                    while let Some(edge) = edges.next() {
                        incoming_weight += edge.weight().clone();
                    }
                }
                else {
                    let mut edges = graph.edges_connecting(NodeIndex::new(*edge_node), NodeIndex::new(*edge_node_parent));
                    while let Some(edge) = edges.next() {
                        incoming_weight += edge.weight().clone();
                    }
                }
                parallel_num_incoming_seq.push(incoming_weight as usize);
                seq_found_so_far += incoming_weight as usize;
            }
            
        }
    }
    (parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, debug_strings)
}

fn get_direction_nodes (direction: Direction, iteration: usize, mut direction_node_list: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    //forward outgoing
    //backward incoming
    if iteration <= 0 {
        return direction_node_list;
    }
    //get the back nodes of the target
    let mut direction_neighbours = graph.neighbors_directed(NodeIndex::new(focus_node), direction);
    //iterate through the neighbours
    while let Some(direction_neighbour) = direction_neighbours.next() {
        if !direction_node_list.contains(&direction_neighbour.index()){
            direction_node_list.push(direction_neighbour.index());
            direction_node_list = get_direction_nodes (direction, iteration - 1, direction_node_list, direction_neighbour.index(), graph);
        }
    }
    direction_node_list
}

fn get_xiterations_direction_nodes (direction: Direction ,iteration: usize, mut direction_node_list: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    if iteration <= 0 {
        return direction_node_list;
    }
    //get the back nodes of the target
    let mut direction_neighbours = graph.neighbors_directed(NodeIndex::new(focus_node), direction);
    //iterate through the neighbours
    while let Some(direction_neighbour) = direction_neighbours.next() {
        if iteration == 1 {
            if !direction_node_list.contains(&direction_neighbour.index()){
                direction_node_list.push(direction_neighbour.index());
            }
        }
        direction_node_list = get_xiterations_direction_nodes (direction, iteration - 1, direction_node_list, direction_neighbour.index(), graph);
    }
    direction_node_list
}

fn base_quality_score_calculation (mut total_seq: usize, indices_of_parallel_nodes: Vec<usize>, seq_through_parallel_nodes: Vec<usize>, base: u8, graph: &Graph<u8, i32, Directed, usize>) -> (f64, bool, Vec<usize>, Vec<String>) {
    if USEPACBIODATA {
        total_seq -= 1;
    }
    //variable initialization
    let mut debug_strings: Vec<String> = Vec::new();
    let mut count_mismatch: bool = false;
    let error_score: f64;
    let quality_score;
    let base_counts: Vec<usize>;

    let ln_prob_base_a = 0.25_f64.ln();
    let ln_prob_base_c = 0.25_f64.ln();
    let ln_prob_base_g = 0.25_f64.ln();
    let ln_prob_base_t = 0.25_f64.ln();
    
    let mut base_a_count = 0;
    let mut base_c_count = 0;
    let mut base_g_count = 0;
    let mut base_t_count = 0;
    //find out how many sequences run through each base
    //match the indices to the base and ++
    for index in 0..indices_of_parallel_nodes.len() {
        match graph.raw_nodes()[indices_of_parallel_nodes[index]].weight {
            65 => {
                base_a_count += seq_through_parallel_nodes[index];
            },
            67 => {
                base_c_count += seq_through_parallel_nodes[index];
            },
            71 => {
                base_g_count += seq_through_parallel_nodes[index];
            },
            84 => {
                base_t_count += seq_through_parallel_nodes[index];
            },
            _ => {
                //nothing
                },
        }
    }
    // save the base counts for debug
    base_counts = [base_a_count, base_c_count, base_g_count, base_t_count].to_vec();
    if (base_a_count + base_c_count + base_g_count + base_t_count) != (total_seq) {
        count_mismatch = true;
    }
    match count_mismatch {
        true => {
            let temp_string = format!("base counts A:{} C:{} G:{} T:{} MISMATCHHHH!!!!!!!!!!!!!!!!!!!!! \n", base_a_count, base_c_count, base_g_count, base_t_count);
            println!("{}", temp_string);
            debug_strings.push(temp_string.clone());
        },
        false => {
            let temp_string = format!("base counts A:{} C:{} G:{} T:{}\n", base_a_count, base_c_count, base_g_count, base_t_count);
            println!("{}", temp_string);
            debug_strings.push(temp_string.clone());
        }
    }
    
    
    //calculate all the probablilities
    let ln_prob_data_given_a = calculate_binomial(total_seq, base_a_count, ERROR_PROBABILITY).ln();
    let ln_prob_data_given_c = calculate_binomial(total_seq, base_c_count, ERROR_PROBABILITY).ln();
    let ln_prob_data_given_g = calculate_binomial(total_seq, base_g_count, ERROR_PROBABILITY).ln();
    let ln_prob_data_given_t = calculate_binomial(total_seq, base_t_count, ERROR_PROBABILITY).ln();
    //println!("D|A:{} D|C:{} D|G:{} D|T:{}", ln_prob_data_given_a, ln_prob_data_given_c, ln_prob_data_given_g, ln_prob_data_given_t);
    //get the error score *changed to log space*
    let mut ln_sum_of_probablities = ln_prob_data_given_a + ln_prob_base_a;
    ln_sum_of_probablities = ln_sum_of_probablities.ln_add_exp(ln_prob_data_given_c + ln_prob_base_c);
    ln_sum_of_probablities = ln_sum_of_probablities.ln_add_exp(ln_prob_data_given_g + ln_prob_base_g);
    ln_sum_of_probablities = ln_sum_of_probablities.ln_add_exp(ln_prob_data_given_t + ln_prob_base_t);

    error_score =  1.0 - match base {
        65 => {
            //println!("Focus base : A" );
            exp(ln_prob_data_given_a + ln_prob_base_a - ln_sum_of_probablities)
        },
        67 => {
            //println!("Focus base : C" );
            exp(ln_prob_data_given_c + ln_prob_base_c - ln_sum_of_probablities)
        },
        71 => {
            //println!("Focus base : G" );
            exp(ln_prob_data_given_g + ln_prob_base_g - ln_sum_of_probablities)
        },
        84 => {
            //println!("Focus base : T" );
            exp(ln_prob_data_given_t + ln_prob_base_t - ln_sum_of_probablities)
        },
        _ => {0.0},
    };
    quality_score = (-10.00) * error_score.log10();
    //println!("quality score: {}", quality_score);
    //println!("");
    (quality_score, count_mismatch, base_counts, debug_strings)
}

fn calculate_binomial (n: usize, k: usize, prob: f64) -> f64 {
    let binomial_coeff = binomial(n as u64, k as u64);
    let success: f64 = prob.powf(k as f64);
    let failure: f64 = (1.00 - prob).powf((n - k) as f64);
    binomial_coeff * success * failure
}

fn find_the_seq_passing_through (target_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> usize {
    let node_index = NodeIndex::new(target_node);
    //edges directed toward the base
    let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(node_index, Incoming).collect();
    let mut incoming_weight = 0;
    for incoming_node in incoming_nodes {
        let mut edges = graph.edges_connecting(incoming_node, node_index);
        while let Some(edge) = edges.next() {
            incoming_weight += edge.weight().clone();
        }
    }
    //edges directed from the base
    let outgoing_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(node_index, Outgoing).collect();
    let mut outgoing_weight = 0;
    for outgoing_node in outgoing_nodes {
        let mut edges = graph.edges_connecting(node_index, outgoing_node);
        while let Some(edge) = edges.next() {
            outgoing_weight += edge.weight().clone();
        }
    }
    cmp::max(outgoing_weight, incoming_weight) as usize
}

fn find_neighbouring_indices (num_of_iterations: usize, focus_node: usize, graph: &Graph<u8, i32, Directed, usize> ) -> Vec<usize> {
    let mut indices: Vec<usize> = vec![];
    if num_of_iterations <= 0 {
        return indices;
    }
    let mut immediate_neighbours = graph.neighbors_undirected(NodeIndex::new(focus_node));
    while let Some(neighbour_node) = immediate_neighbours.next() {
        if !indices.contains(&neighbour_node.index()){
            indices.push(neighbour_node.index());
        }
        let obtained_indices = find_neighbouring_indices(num_of_iterations - 1, neighbour_node.index(), graph);
        for obtained_index in obtained_indices {
            if !indices.contains(&obtained_index){
                indices.push(obtained_index);
            }
        }
    }
    indices 
}

fn modify_dot_graph_with_highlight (mut dot: String, focus_node: &usize, error_count: &usize, description_type: usize) -> String {
    match dot.find(&format!(" {} [", focus_node)) {
        Some(mut x) => {
            while dot.chars().nth(x).unwrap() != '[' {
                x += 1;
            }
            //normal_dot.replace_range(x + 17..x + 17,&format!("color = \"red\" style = \"filled\"").to_string());
            match description_type {
                0 => {dot.replace_range(x + 12..x + 12,&format!("Mismatch {}: ", error_count).to_string());},
                1 => {dot.replace_range(x + 12..x + 12,&format!("Insert {}: ", error_count).to_string());},
                2 => {dot.replace_range(x + 12..x + 12,&format!("Delete {}: ", error_count).to_string());},
                3 => {dot.replace_range(x + 12..x + 12,&format!("Target node{}: ", error_count).to_string());},
                _ => {}
            };
            
        },
        None => {}
    };
    dot
}

fn get_zoomed_graph_section (normal_graph: &Graph<u8, i32, Directed, usize>, focus_node: &usize, error_count: &usize, description_type: usize)-> String {
    let mut graph_section= "".to_string();
    let normal_dot = format!("{:?}", Dot::new(&normal_graph.map(|_, n| (*n) as char, |_, e| *e)));
    let displaying_nodes: Vec<usize> = find_neighbouring_indices (NUM_OF_ITER_FOR_ZOOMED_GRAPHS, *focus_node, normal_graph);
    let mut graph_section_nodes: String = "".to_string();
    let mut graph_section_edges: String = "".to_string();
    //find the position in the dot file and add to the graph section
    for node in &displaying_nodes {
        //get the nodes from the dot file
        match normal_dot.find(&format!(" {} [", node)) {
            Some(start) => {
                let mut char_seq = vec![];
                let mut end = start;
                loop {
                    char_seq.push(normal_dot.chars().nth(end).unwrap());
                    if normal_dot.chars().nth(end).unwrap() == '\n' {
                        break;
                    }
                    end += 1;
                }
                //add from start to end to the graph_section string
                graph_section_nodes = char_seq.iter().collect::<String>();
            },
            None => {}
        }
        graph_section = format!("{}{}", graph_section, graph_section_nodes).to_string();
        graph_section_nodes = "".to_string();
    }
    for node in &displaying_nodes {
        //get the edges from the dot file
        let edge_entries: Vec<usize> = normal_dot.match_indices(&format!(" {} ->", node)).map(|(i, _)|i).collect();
        for edge in edge_entries {
            let mut char_seq = vec![];
                let mut end = edge;
                loop {
                    char_seq.push(normal_dot.chars().nth(end).unwrap());
                    if normal_dot.chars().nth(end).unwrap() == '\n' {
                        break;
                    }
                    end += 1;
                }
                //add from start to end to the graph_section string
                graph_section_edges = format!("{}{}", graph_section_edges, char_seq.iter().collect::<String>()).to_string();
        }
        graph_section = format!("{}{}", graph_section, graph_section_edges).to_string();
        graph_section_edges = "".to_string();
    }
    //modifying the section graph with the highlight
    graph_section = modify_dot_graph_with_highlight (graph_section, focus_node, error_count, description_type);
    //make it a dot graph
    graph_section = format!("digraph {{\n{} }}", graph_section);
    graph_section
}

fn get_random_sequences_from_generator(sequence_length: usize, num_of_sequences: usize) -> Vec<String> {
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

fn calculate_and_get_expansion (homopolymer_vec: &Vec<HomopolymerSequence>, homopolymer_consensus: &Vec<u8>, homopolymer_consensus_freq: &Vec<Vec<u32>>) -> (Vec<u8>, Vec<u8>) {
    //make the expanded consensus using the frequencies
    let mut expanded_consensus: Vec<u8> = vec![];
    let mut repetitions: Vec<f32> = vec![0.0; homopolymer_consensus.len()];
    let mut homopolymervec_expanded: Vec<u8> = vec![];

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
            homopolymervec_expanded.push((repetitions[j] / homopolymer_vec.len() as f32).round() as u8);
        }
    }
    //++ median ++ 
    if CONSENSUS_METHOD == 1 {
        let mut cpy_homopolymner_consensus_freq = homopolymer_consensus_freq.clone();
        //reorder the homopolymer_consensus_freq by ascending order
        for i in 0..cpy_homopolymner_consensus_freq.len(){
            cpy_homopolymner_consensus_freq[i].sort();
            repetitions[i] = cpy_homopolymner_consensus_freq[i][(homopolymer_vec.len() / 2) as usize] as f32;
        }
        for j in 0..homopolymer_consensus.len() {
            for _ in 0..((repetitions[j]).round() as usize) {
                expanded_consensus.push(homopolymer_consensus[j]);
            }
            homopolymervec_expanded.push((repetitions[j]).round() as u8);
        }
    }
    //++ mode ++ if unwrap fails median used
    if CONSENSUS_METHOD == 2 {
        let mut cpy_homopolymner_consensus_freq = homopolymer_consensus_freq.clone();
        //get the mode of each base
        for i in 0..homopolymer_consensus.len(){
            cpy_homopolymner_consensus_freq[i].sort();
            let mut counts = HashMap::new();
            let mode = cpy_homopolymner_consensus_freq[i].iter().copied().max_by_key(|&n| {
                let count = counts.entry(n).or_insert(0);
                *count += 1;
                *count
            });
            match mode {
                Some(value) => {
                    repetitions[i] = value as f32;
                },
                None => {
                    repetitions[i] = cpy_homopolymner_consensus_freq[i][(homopolymer_vec.len() / 2) as usize] as f32;
                }
            }
        }
        for j in 0..homopolymer_consensus.len() {
            for _ in 0..((repetitions[j]).round() as usize) {
                expanded_consensus.push(homopolymer_consensus[j]);
            }
            homopolymervec_expanded.push((repetitions[j]).round() as u8);
        }
    }
    (expanded_consensus, homopolymervec_expanded)
}

fn get_aligned_homopolymer_sequences_to_homopolymer_consensus (homopolymer_vec: &Vec<HomopolymerSequence>, homopolymer_consensus: &Vec<u8>) -> (Vec<Vec<u32>>, i32)  {
    //use homopolymer compressions sequences to make expanded consensus //make function
    let mut homopolymer_score = 0;
    let mut homopolymer_consensus_freq: Vec<Vec<u32>> = vec![vec![0; homopolymer_vec.len()]; homopolymer_consensus.len()];
    let mut i = 0;
    for homopolymer_seq in homopolymer_vec{
        let mut sequence_base_freq: Vec<u32> = vec![0; homopolymer_seq.bases.len()];
        let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
        let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(homopolymer_consensus.len(), homopolymer_seq.bases.len(), GAP_OPEN, GAP_EXTEND, &score);
        let alignment = aligner.global(&homopolymer_consensus, &homopolymer_seq.bases);
        let mut consensus_index = alignment.xstart;
        let mut homopolymer_index = alignment.ystart;
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
        homopolymer_score += alignment.score;
        i += 1;

    }
    (homopolymer_consensus_freq, homopolymer_score)
}

fn get_indices_for_homopolymer_debug(alignment: &bio::alignment::Alignment, homopolymer_expand: &Vec<u8>, normal_topo: &Vec<usize>, homopolymer_topo: &Vec<usize>) 
                            -> IndexStruct {

    let mut normal_mismatches: Vec<usize> = vec![];
    let mut normal_insertions: Vec<usize> = vec![];
    let mut normal_deletions: Vec<usize> = vec![];

    let mut expanded_mismatches: Vec<usize> = vec![];
    let mut expanded_insertions: Vec<usize> = vec![];
    let mut expanded_deletions: Vec<usize> = vec![];

    let mut homopolymer_mismatches: Vec<usize> = vec![];
    let mut homopolymer_insertions: Vec<usize> = vec![];
    let mut homopolymer_deletions: Vec<usize> = vec![];

    let mut alignment_mismatches: Vec<usize> = vec![];
    let mut alignment_insertions: Vec<usize> = vec![];
    let mut alignment_deletions: Vec<usize> = vec![];
    let mut alignment_index: usize = 0;

    let mut normal_index: usize = alignment.xstart;
    let mut expanded_index: usize = alignment.ystart;
    for op in &alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                normal_index += 1;
                expanded_index += 1;
            },
            bio::alignment::AlignmentOperation::Subst => {
                normal_mismatches.push(normal_index);
                expanded_mismatches.push(expanded_index);
                alignment_mismatches.push(alignment_index);
                normal_index += 1;
                expanded_index += 1;
            },
            bio::alignment::AlignmentOperation::Del => {
                expanded_insertions.push(expanded_index);
                normal_deletions.push(normal_index);
                alignment_deletions.push(alignment_index);
                expanded_index += 1;
            },
            bio::alignment::AlignmentOperation::Ins => {
                normal_insertions.push(normal_index);
                expanded_deletions.push(expanded_index);
                alignment_insertions.push(alignment_index);
                normal_index += 1;
            },
            _ => {},
        }
        alignment_index += 1;
    }
    // calculate the homopolymer stuff from expanded stuff
    for expanded_entry in &expanded_mismatches {
        let mut index = 0;
        let mut homopolymer_index = 0;
        for homopolymer_freq in homopolymer_expand {
            if expanded_entry <= &index {
                homopolymer_mismatches.push(homopolymer_index);
                break;
            }
            index += *homopolymer_freq as usize;
            homopolymer_index += 1;
        }
    }
    for expanded_entry in &expanded_insertions {
        let mut index = 0;
        let mut homopolymer_index = 0;
        for homopolymer_freq in homopolymer_expand {
            if expanded_entry <= &index {
                homopolymer_insertions.push(homopolymer_index);
                break;
            }
            index += *homopolymer_freq as usize;
            homopolymer_index += 1;
        }
    }
    for expanded_entry in &expanded_deletions {
        let mut index = 0;
        let mut homopolymer_index = 0;
        for homopolymer_freq in homopolymer_expand {
            if expanded_entry <= &index {
                homopolymer_deletions.push(homopolymer_index);
                break;
            }
            index += *homopolymer_freq as usize;
            homopolymer_index += 1;
        }
    }
    //calculate normal and homopolymer positions in respective graphs using the topology indices
    for index in 0..normal_mismatches.len() {
        normal_mismatches[index] = normal_topo[normal_mismatches[index] as usize] as usize;
    }
    for index in 0..normal_insertions.len() {
        normal_insertions[index] = normal_topo[normal_insertions[index] as usize] as usize;
    }
    for index in 0..normal_deletions.len() {
        normal_deletions[index] = normal_topo[normal_deletions[index] as usize] as usize;
    }
    for index in 0..homopolymer_mismatches.len() {
        homopolymer_mismatches[index] = homopolymer_topo[homopolymer_mismatches[index] as usize] as usize;
    }
    for index in 0..homopolymer_insertions.len() {
        homopolymer_insertions[index] = homopolymer_topo[homopolymer_insertions[index] as usize] as usize;
    }
    for index in 0..homopolymer_deletions.len() {
        homopolymer_deletions[index] = homopolymer_topo[homopolymer_deletions[index] as usize] as usize;
    }
    IndexStruct {
        normal_mismatch_indices: normal_mismatches,
        normal_mismatch_graph_sections: vec![],
        normal_insert_indices: normal_insertions,
        normal_insert_graph_sections: vec![],
        normal_del_indices: normal_deletions,
        normal_del_graph_sections: vec![],
        homopolymer_mismatch_indices: homopolymer_mismatches,
        homopolymer_mismatch_graph_sections: vec![],
        homopolymer_insert_indices: homopolymer_insertions,
        homopolymer_insert_graph_sections: vec![],
        homopolymer_del_indices: homopolymer_deletions,
        homopolymer_del_graph_sections: vec![],
        aligned_mismatch_indices: alignment_mismatches,
        aligned_insert_indices: alignment_insertions,
        aligned_del_indices: alignment_deletions,
    }
}
fn get_alignment_with_count_for_debug(vector1: &Vec<u8>, vector2: &Vec<u8>, alignment: &bio::alignment::Alignment, count: &Vec<Vec<u32>>, sequence_num: usize) -> (Vec<u8>, Vec<u8>, Vec<Vec<u32>>){
    let mut vec1_representation = vec![];
    let mut vec2_representation = vec![];
    let mut vec1_index: usize = alignment.xstart;
    let mut vec2_index: usize = alignment.ystart;
    let mut prev_base: u8 = 0; //expanded consensus
    let mut same_base = false;
    //initializae count_representation
    let mut count_representation: Vec<Vec<u32>> = vec![vec![]; sequence_num];
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
            for i in 0..sequence_num {
                count_representation[i].push(0);
            }
            
        }
        else {
            for i in 0..sequence_num {
                count_representation[i].push(count[count_index][i]);
            }
            count_index += 1;
        }
        if vec2_index != 0 {
            prev_base = vector2[vec2_index - 1];
        }
    }
    //write
    (vec1_representation, vec2_representation, count_representation)
    //write_alignment_data_fasta_file("./results/consensus.fa", &vec1_representation, &vec2_representation, &count_representation, sequence_num);
}

//structs here
pub struct IndexStruct {
    pub normal_mismatch_indices: Vec<usize>,
    pub normal_mismatch_graph_sections: Vec<String>,
    pub normal_insert_indices: Vec<usize>,
    pub normal_insert_graph_sections: Vec<String>,
    pub normal_del_indices: Vec<usize>,
    pub normal_del_graph_sections: Vec<String>,

    pub homopolymer_mismatch_indices: Vec<usize>,
    pub homopolymer_mismatch_graph_sections: Vec<String>,
    pub homopolymer_insert_indices: Vec<usize>,
    pub homopolymer_insert_graph_sections: Vec<String>,
    pub homopolymer_del_indices: Vec<usize>,
    pub homopolymer_del_graph_sections: Vec<String>,

    pub aligned_mismatch_indices: Vec<usize>,
    pub aligned_insert_indices: Vec<usize>,
    pub aligned_del_indices: Vec<usize>,
}
impl fmt::Display for IndexStruct {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "nm{:?} ni{:?} nd{:?} hm{:?} hi{:?} hd{:?} am{:?} ai{:?} ad{:?}",
            self.normal_mismatch_indices, self.normal_insert_indices, self.normal_del_indices,
            self.homopolymer_mismatch_indices, self.homopolymer_insert_indices, self.homopolymer_del_indices,
            self.aligned_mismatch_indices, self.aligned_insert_indices, self.aligned_del_indices
            )
    }
}
//file stuff

fn get_quality_from_file (filename: impl AsRef<Path>) -> String {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    let lines: Vec<String> = buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();
    let quality: String = lines[3].clone().into();
    quality
}

fn get_consensus_from_file (filename: impl AsRef<Path>) -> String {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    let lines: Vec<String> = buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();
    let consensus: String = lines[1].clone().into();
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

//write stuff here 

fn write_alignment_and_zoomed_graphs_fasta_file (filename: impl AsRef<Path>, normal: &Vec<u8>, expanded: &Vec<u8>, count_representation: &Vec<Vec<u32>>, sequence_num: usize, saved_indices: &IndexStruct) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n>Normal consensus vs Expanded consensus with counts:",
        chrono::offset::Local::now(), INPUT_FILE_NAME)
        .expect("result file cannot be written");

    let mut index = 0;
    println!("{}", saved_indices);
    while index + 50 < normal.len() {
        let mut write_string: Vec<String> = vec![];
        //count the mismatches, inserts and del in that range
        let mut mismatch_count = 0;
        let mut mismatch_struct_indices: Vec<usize> = vec![];
        let mut insert_count = 0;
        let mut insert_struct_indices: Vec<usize> = vec![];
        let mut del_count = 0;
        let mut del_struct_indices: Vec<usize> = vec![];
        for mismatch_index in index..index + 50 {
            if saved_indices.aligned_mismatch_indices.contains(&mismatch_index) {
                mismatch_count += 1;
                mismatch_struct_indices.push(saved_indices.aligned_mismatch_indices.iter().position(|&r| r == mismatch_index).unwrap());
            }
            if saved_indices.aligned_insert_indices.contains(&mismatch_index) {
                insert_count += 1;
                insert_struct_indices.push(saved_indices.aligned_insert_indices.iter().position(|&r| r == mismatch_index).unwrap());
            }
            if saved_indices.aligned_del_indices.contains(&mismatch_index) {
                del_count += 1;
                del_struct_indices.push(saved_indices.aligned_del_indices.iter().position(|&r| r == mismatch_index).unwrap());
            }
        }
        write_string.push(format!("\n{}~{} out of {} (mismatches:{}, inserts:{}, deletions:{})\n", index, index + 50, normal.len(), mismatch_count, insert_count, del_count).to_string());
        write_string.push("normal:".to_string());
        for i in index..index + 50 {
            match normal[i] {
                55 => write_string.push("  _".to_string()),
                65 => write_string.push("  A".to_string()),
                67 => write_string.push("  C".to_string()),
                71 => write_string.push("  G".to_string()),
                84 => write_string.push("  T".to_string()),
                _ => {},
            }
            if saved_indices.aligned_mismatch_indices.contains(&i) {
                write_string.push("*".to_string());
            }
            else if saved_indices.aligned_insert_indices.contains(&i) {
                write_string.push("%".to_string());
            }
            else if saved_indices.aligned_del_indices.contains(&i) {
                write_string.push("?".to_string());
            }
            else {
                write_string.push(" ".to_string());
            }
        }
        write_string.push(format!("\nexpand:"));
        for i in index..index + 50 {
            match expanded[i] {
                55 => write_string.push("  _".to_string()),
                65 => write_string.push("  A".to_string()),
                67 => write_string.push("  C".to_string()),
                71 => write_string.push("  G".to_string()),
                84 => write_string.push("  T".to_string()),
                _ => {},
            }
            if saved_indices.aligned_mismatch_indices.contains(&i) {
                write_string.push("*".to_string());
            }
            else if saved_indices.aligned_insert_indices.contains(&i) {
                write_string.push("%".to_string());
            }
            else if saved_indices.aligned_del_indices.contains(&i) {
                write_string.push("?".to_string());
            }
            else {
                write_string.push(" ".to_string());
            }
        }
        write_string.push("\n".to_string());
        for j in 0..sequence_num {
            write_string.push(format!("seq{:>3}:", j).to_string());
            for i in index..index + 50 {
                write_string.push(format!("{:>3},", count_representation[j][i]).to_string()); 
            }
            write_string.push("\n".to_string());
        }
        write_string.push("\n".to_string());
        //push the graphs to the vector
        for entry in &mismatch_struct_indices {
            write_string.push("\nNormal Graph: mismatch\n".to_string());
            write_string.push(saved_indices.normal_mismatch_graph_sections[*entry].clone());
        }
        for entry in &insert_struct_indices {
            write_string.push("\nNormal Graph: insert\n".to_string());
            write_string.push(saved_indices.normal_insert_graph_sections[*entry].clone());
        }
        for entry in &del_struct_indices {
            write_string.push("\nNormal Graph: del\n".to_string());
            write_string.push(saved_indices.normal_del_graph_sections[*entry].clone());
        }
        for entry in &mismatch_struct_indices {
            write_string.push("\nHomopolymer Graph: mismatch\n".to_string());
            write_string.push(saved_indices.homopolymer_mismatch_graph_sections[*entry].clone());
        }
        for entry in &insert_struct_indices {
            write_string.push("\nHomopolymer Graph: del\n".to_string());
            write_string.push(saved_indices.homopolymer_del_graph_sections[*entry].clone());
        }
        for entry in &del_struct_indices {
            write_string.push("\nHomopolymer Graph: insert\n".to_string());
            write_string.push(saved_indices.homopolymer_insert_graph_sections[*entry].clone());
        }
        index = index + 50;
        for entry in write_string{
            //print!("{}", entry);
            write!(file,
                "{}",
                entry)
                .expect("result file cannot be written");
        }
    }
}

fn write_quality_score_graph (normal_filename: impl AsRef<Path>, normal_graph: &Graph<u8, i32, Directed, usize>) {
    let normal_dot = format!("{:?}", Dot::new(&normal_graph.map(|_, n| (*n) as char, |_, e| *e)));
    let mut normal_file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(normal_filename)
        .unwrap();
    writeln!(normal_file,
        "{:?} \nFILE: {}\n{}",
        chrono::offset::Local::now(), INPUT_FILE_NAME, normal_dot)
        .expect("result file cannot be written");
}

fn modify_and_write_the_graphs_and_get_zoomed_graphs (normal_filename: impl AsRef<Path>, homopolymer_filename: impl AsRef<Path>, mut saved_indices: IndexStruct,
    normal_graph: &Graph<u8, i32, Directed, usize>, homopolymer_graph: &Graph<u8, i32, Directed, usize>)
    -> IndexStruct {
    let mut count = 0;
    let mut normal_dot = format!("{:?}", Dot::new(&normal_graph.map(|_, n| (*n) as char, |_, e| *e)));
    let mut homopolymer_dot = format!("{:?}", Dot::new(&homopolymer_graph.map(|_, n| (*n) as char, |_, e| *e)));

    for index in &saved_indices.normal_mismatch_indices {
        let graph_section = get_zoomed_graph_section (&normal_graph, &index, &count, 0);
        saved_indices.normal_mismatch_graph_sections.push(graph_section);
        normal_dot = modify_dot_graph_with_highlight(normal_dot, &index, &count, 0);
        count += 1;
    }
    count = 0;
    for index in &saved_indices.normal_insert_indices {
        let graph_section = get_zoomed_graph_section (&normal_graph, &index, &count, 1);
        saved_indices.normal_insert_graph_sections.push(graph_section);
        normal_dot = modify_dot_graph_with_highlight(normal_dot, &index, &count, 1);
        count += 1;
    }
    count = 0;
    for index in &saved_indices.normal_del_indices {
        let graph_section = get_zoomed_graph_section (&normal_graph, &index, &count, 2);
        saved_indices.normal_del_graph_sections.push(graph_section);
        normal_dot = modify_dot_graph_with_highlight(normal_dot, &index, &count, 2);
        count += 1;
    }
    count = 0;
    for index in &saved_indices.homopolymer_mismatch_indices {
        let graph_section = get_zoomed_graph_section (&homopolymer_graph, &index, &count, 0);
        saved_indices.homopolymer_mismatch_graph_sections.push(graph_section);
        homopolymer_dot = modify_dot_graph_with_highlight(homopolymer_dot, &index, &count, 0);
        count += 1;
    }
    count = 0;
    for index in &saved_indices.homopolymer_insert_indices {
        let graph_section = get_zoomed_graph_section (&homopolymer_graph, &index, &count, 1);
        saved_indices.homopolymer_insert_graph_sections.push(graph_section);
        homopolymer_dot = modify_dot_graph_with_highlight(homopolymer_dot, &index, &count, 1);
        count += 1;
    }
    count = 0;
    for index in &saved_indices.homopolymer_del_indices {
        let graph_section = get_zoomed_graph_section (&homopolymer_graph, &index, &count, 2);
        saved_indices.homopolymer_del_graph_sections.push(graph_section);
        homopolymer_dot = modify_dot_graph_with_highlight(homopolymer_dot, &index, &count, 2);
        count += 1;
    } 
    let mut normal_file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(normal_filename)
        .unwrap();
    writeln!(normal_file,
        "{:?} \nFILE: {}\n{}",
        chrono::offset::Local::now(), INPUT_FILE_NAME, normal_dot)
        .expect("result file cannot be written");
    let mut homopolymer_file = OpenOptions::new()
    .write(true)
    .append(true)
    .open(homopolymer_filename)
    .unwrap();
    writeln!(homopolymer_file,
        "{:?} \nFILE: {}\n{}",
        chrono::offset::Local::now(), INPUT_FILE_NAME, homopolymer_dot)
        .expect("result file cannot be written");
    saved_indices
    //println!("{}", normal_dot);
    //println!("{}", homopolymer_dot);
}
fn write_scores_result_file(filename: impl AsRef<Path>, normal_score: i32, homopolymer_score: i32, expanded_score: i32) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
            "{:?} \nFILE: {}\nNormal score:\t\t\t{}\nHomopolymer score:\t\t{}\nExpanded score:\t\t\t{}",
            chrono::offset::Local::now(), INPUT_FILE_NAME, normal_score, homopolymer_score, expanded_score)
            .expect("result file cannot be written");
}

fn write_quality_scores_to_file (filename: impl AsRef<Path>, quality_scores: &Vec<f64>, consensus: &Vec<u8>, topology: &Vec<usize>, invalid_info: &Vec<(usize, usize, bool, bool, bool, bool)>, base_count_vec: &Vec<Vec<usize>>, pacbioquality: &Vec<usize>, aligned_pacbio_bases: &Vec<u8>) -> Vec<String> {
    let mut saved_preceeding: Vec<String> = ["".to_string(), "".to_string(), "".to_string()].to_vec();
    let mut saved_proceeding: Vec<String> = vec![];
    let mut error_found = false;
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n>Quality score data & graphs:",
        chrono::offset::Local::now(), INPUT_FILE_NAME)
        .expect("result file cannot be written");
    for index in 0..consensus.len() {
        let mut print_info = (false, false, false, false);
        match invalid_info.iter().position(|&r| r.0 == index) {
            Some(position) => {
                print_info = (invalid_info[position].2, invalid_info[position].3, invalid_info[position].4, invalid_info[position].5);
                if invalid_info[position].5 {
                    error_found = true;
                }
            },
            None => {},
        }
        let output_string = format!("{} {} [{:>6}]\t\t -> {:>8.3}[{}] \tinvalidity =[p,m,q,e] {:?}\t base_counts = ACGT{:?}", consensus[index] as char, aligned_pacbio_bases[index] as char, topology[index], quality_scores[index], (pacbioquality[index % pacbioquality.len()] - 33), print_info, base_count_vec[index]);
        if error_found == false {
            //rearrange preceding strings
            saved_preceeding[0] = saved_preceeding[1].clone();
            saved_preceeding[1] = saved_preceeding[2].clone();
            //save to preceding vector
            saved_preceeding[2] = output_string.clone();
        }
        else if saved_proceeding.len() < 4{
            //save the proceeding vector
            saved_proceeding.push(output_string.clone());
        }   
        writeln!(file, "{}", output_string).expect("result file cannot be written");
    }
    [saved_preceeding, saved_proceeding].concat()
}

fn write_zoomed_quality_score_graphs (filename: impl AsRef<Path>, invalid_info: &Vec<(usize, usize, bool, bool, bool, bool)>, quality_scores: &Vec<f64>, base_count_vec: &Vec<Vec<usize>>, graph: &Graph<u8, i32, Directed, usize>, pacbioquality: &Vec<usize>) -> Vec<String>{
    let mut return_string = vec![];
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n>Quality score data & graphs:",
        chrono::offset::Local::now(), INPUT_FILE_NAME)
        .expect("result file cannot be written");
    for entry in invalid_info {
        let graph_section = get_zoomed_graph_section (graph, &entry.1, &0, 3);
        let temp_output = format!("\nnode_index:{}\t\tnode_base:{}\t\tquality_score:{:.3}[{}] invalidity:[p,m,q,e]{:?}\tbase_count:ACGT{:?}\t\n{}\n",
        entry.1, graph.raw_nodes()[entry.1].weight as char, quality_scores[entry.0], (pacbioquality[entry.0 % pacbioquality.len()] - 33), (entry.2, entry.3, entry.4, entry.5), base_count_vec[entry.0], graph_section);
        writeln!(file, "{}", temp_output).expect("result file cannot be written");
        if entry.5 {
            return_string.push(temp_output.clone());
        }
    }
    return_string
}