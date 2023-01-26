use bio::alignment::pairwise::Scoring;
use bio::alignment::{poa::*, TextSlice}; //bandedpoa/ poa
use bio::io::fasta::Index;
use petgraph::visit::IntoNeighborsDirected;
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
use petgraph::{Directed, Graph, Incoming, Outgoing};
use std::collections::HashMap;
use petgraph::graph::NodeIndex;
use std::{fmt, cmp};
use statrs::function::factorial::binomial;
use num_traits::pow;

const GAP_OPEN: i32 = -4;
const GAP_EXTEND: i32 = -2;
const MATCH: i32 = 2;
const MISMATCH: i32 = -4;
const FILENAME: &str = "./data/161808690.fasta";
const SEED: u64 = 0;
const CONSENSUS_METHOD: u8 = 1; //0==average 1==median //2==mode
const ERROR_PROBABILITY: f64 = 0.90;
const DEBUG: bool = false;

fn main() {
    //let seqvec = get_fasta_sequences_from_file(FILENAME);
    let seqvec = get_random_sequences_from_generator(100, 10);
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
    let normal_topology;
    (normal_consensus, normal_topology) = aligner.poa.consensus(); //just poa

    //get scores of sequences compared to normal consensus 
    let normal_score = get_consensus_score(&seqvec, &normal_consensus);
    //get the normal graph
    let normal_graph = aligner.graph();
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
    let homopolymer_topology;
    (homopolymer_consensus, homopolymer_topology) = aligner.poa.consensus(); //poa
    //get graph
    let homopolymer_graph: &Graph<u8, i32, Directed, usize> = aligner.graph();
    //use homopolymer compressions sequences to make expanded consensus
    let (expanded_consensus, homopolymer_consensus_freq, homopolymer_score, homopolymer_expanded) =  get_expanded_consensus(homopolymer_vec, &homopolymer_consensus);
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
    print!("\nExpanded consensus:\t\t");
    for i in &expanded_consensus{
        print!("{}", *i as char);
    }
    print!("\nExpanded consensus score:\t{}", expanded_score);
    println!("");

    //DEBUGGING
    if DEBUG == true {
        let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
        let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(normal_consensus.len(), expanded_consensus.len(), GAP_OPEN, GAP_EXTEND, &score);
        let alignment = aligner.global(&normal_consensus, &expanded_consensus);
        let mut saved_indices: IndexStruct;
        saved_indices = get_indices_for_debug(&normal_consensus,&expanded_consensus, &alignment, &homopolymer_expanded, &normal_topology, &homopolymer_topology);
        let (normal_rep, expanded_rep, count_rep) 
            = get_alignment_with_count_for_debug(&normal_consensus,&expanded_consensus, &alignment, &homopolymer_consensus_freq, seqnum as usize);
        //write results to file
        write_scores_result_file("./results/results.txt", normal_score, homopolymer_score, expanded_score);
        //modify the graphs to indicate 
        saved_indices = modify_and_write_the_graphs_and_get_zoomed_graphs("./results/normal_graph.fa", "./results/homopolymer_graph.fa", saved_indices, normal_graph, homopolymer_graph);
        write_alignment_and_zoomed_graphs_fasta_file("./results/consensus.fa", &normal_rep, &expanded_rep,
            &count_rep, seqnum as usize, &saved_indices);
    }
}

fn base_quality_score_calculation (total_seq: usize, base: u8, graph_index: usize, graph: &Graph<u8, i32, Directed, usize>) -> f64 {
    //variable initialization
    let mut error_score: f64 = 0.0;
    let mut quality_score = 0.0;
    let mut num_seq_through_base = 0;
    let prob_base_a = 0.25;
    let prob_base_c = 0.25;
    let prob_base_g = 0.25;
    let prob_base_t = 0.25;
    let mut prob_data_given_a = 0.0;
    let mut prob_data_given_c = 0.0;
    let mut prob_data_given_g = 0.0;
    let mut prob_data_given_t = 0.0;
    let node_index = NodeIndex::new(graph_index);

    //find out how many sequences run through the base in graph
    //edges directed toward the base
    let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(node_index, Incoming).collect();
    let mut incoming_weight = 0;
    for incoming_node in incoming_nodes {
        let mut edges = graph.edges_connecting(node_index, incoming_node);
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
    //get the number of seq passing through this base
    num_seq_through_base = cmp::max(outgoing_weight, incoming_weight);
    //calculate all the probablilities

    //get the probability base given data
    let sum_of_probabilities = (prob_data_given_a * prob_base_a) + (prob_data_given_c * prob_base_c) + (prob_data_given_g * prob_base_g) + (prob_data_given_t * prob_base_t);
    error_score = 1.0 - match base {
        65 => {
            (prob_data_given_a * prob_base_a) / sum_of_probabilities
        },
        67 => {
            (prob_data_given_c * prob_base_c) / sum_of_probabilities
        },
        71 => {
            (prob_data_given_g * prob_base_g) / sum_of_probabilities
        },
        84 => {
            (prob_data_given_t * prob_base_t) / sum_of_probabilities
        },
        _ => {0.0},
    };
    quality_score = (-10.00) * error_score.log10(); 
    quality_score
}

fn calculate_binomial (n: usize, k: usize, prob: f64) -> f64 {
    let binomial_coeff = binomial(n as u64, k as u64);
    let success: f64 = prob.powf(n as f64);
    let failure: f64 = (1.00 - prob).powf(k as f64);
    binomial_coeff * success * failure
}

fn write_alignment_and_zoomed_graphs_fasta_file(filename: impl AsRef<Path>, normal: &Vec<u8>, expanded: &Vec<u8>, count_representation: &Vec<Vec<u32>>, sequence_num: usize, saved_indices: &IndexStruct){
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n>Normal consensus vs Expanded consensus with counts:",
        chrono::offset::Local::now(), FILENAME)
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
        chrono::offset::Local::now(), FILENAME, normal_dot)
        .expect("result file cannot be written");
    let mut homopolymer_file = OpenOptions::new()
    .write(true)
    .append(true)
    .open(homopolymer_filename)
    .unwrap();
    writeln!(homopolymer_file,
        "{:?} \nFILE: {}\n{}",
        chrono::offset::Local::now(), FILENAME, homopolymer_dot)
        .expect("result file cannot be written");
    saved_indices
    //println!("{}", normal_dot);
    //println!("{}", homopolymer_dot);
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
                _ => {}
            };
            
        },
        None => {}
    };
    dot
}

fn get_zoomed_graph_section (normal_graph: &Graph<u8, i32, Directed, usize>, focus_node: &usize, error_count: &usize, description_type: usize)-> String {
    let mut graph_section= "".to_string();
    let mut normal_dot = format!("{:?}", Dot::new(&normal_graph.map(|_, n| (*n) as char, |_, e| *e)));
    let displaying_nodes: Vec<usize> = find_neighbouring_indices (2, *focus_node, normal_graph);
    let mut graph_section_nodes: String = "".to_string();
    let mut graph_section_edges: String = "".to_string();
    //find the position in the dot file and add to the graph section
    for node in &displaying_nodes {
        //get the nodes from the dot file
        match normal_dot.find(&format!(" {} [", node)) {
            Some(start) => {
                let mut char_seq = vec![];
                let mut end = start;
                while true {
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
                while true {
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

fn find_neighbouring_indices (num_of_iterations: usize, focus_node: usize, graph: &Graph<u8, i32, Directed, usize> ) -> Vec<usize> {
    let mut indices: Vec<usize> = vec![];
    if num_of_iterations <= 0 {
        return indices;
    }
    let mut immediate_neighbours = graph.neighbors_undirected(NodeIndex::new(focus_node));
    while let Some(neighbour_node) = immediate_neighbours.next() {
        if indices.contains(&neighbour_node.index()) == false{
            indices.push(neighbour_node.index());
        }
        let obtained_indices = find_neighbouring_indices(num_of_iterations - 1, neighbour_node.index(), graph);
        for obtained_index in obtained_indices {
            if indices.contains(&obtained_index) == false {
                indices.push(obtained_index);
            }
        }
    }
    indices 
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


fn get_expanded_consensus(homopolymer_vec: Vec<HomopolymerSequence>, homopolymer_consensus: &Vec<u8>) -> (Vec<u8>, Vec<Vec<u32>>, i32, Vec<u8>)  {
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
        homopolymer_score += alignment.score;
        i += 1;

    }

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
            homopolymervec_expanded.push(((repetitions[j] / homopolymer_vec.len() as f32).round() as u8));
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
    
    (expanded_consensus, homopolymer_consensus_freq, homopolymer_score, homopolymervec_expanded)
}

fn get_indices_for_debug(normal: &Vec<u8>, expanded: &Vec<u8>, alignment: &bio::alignment::Alignment, homopolymer_expand: &Vec<u8>, normal_topo: &Vec<usize>, homopolymer_topo: &Vec<usize>) 
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
    //print
    //print_consensus_with_count(&vec1_representation, &vec2_representation, &count_representation, sequence_num);
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

//write stuff here 
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

//print stuff here
fn print_consensus_with_count(normal: &Vec<u8>, expanded: &Vec<u8>, count_representation: &Vec<Vec<u32>>, sequence_num: usize) {
    let mut index = 0;
    while index + 50 < normal.len() {
        println!("{}~{} out of {}", index, index + 50, normal.len());
        print!("normal:");
        for i in index..index + 50 {
            match normal[i] {
                55 => print!("  _ "),
                65 => print!("  A "),
                67 => print!("  C "),
                71 => print!("  G "),
                84 => print!("  T "),
                _ => {},
            }
        }
        print!("\nexpand:");
        for i in index..index + 50 {
            match expanded[i] {
                55 => print!("  _ "),
                65 => print!("  A "),
                67 => print!("  C "),
                71 => print!("  G "),
                84 => print!("  T "),
                _ => {},
            }
        }
        print!("\n");
        
        for j in 0..sequence_num {
            print!("seq{:>3}:", j);
            for i in index..index + 50 {
                print!("{:>3},", count_representation[j][i]); 
            }
            print!("\n");
        }
        print!("\n");
        index = index + 50;
    }
}

