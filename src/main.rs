use bio::alignment::{poa::*, TextSlice, pairwise::Scoring};
use std::{fs::File, fs::OpenOptions, io::{prelude::*, BufReader}, path::Path, fmt, cmp, collections::HashMap};
use chrono;
use rand::{Rng, SeedableRng, rngs::StdRng};
use petgraph::{Directed, Graph, Incoming, Outgoing, Direction, dot::Dot, graph::NodeIndex, visit::Topo};
use statrs::function::factorial::binomial;
use logaddexp::LogAddExp;
use libm::exp;

const GAP_OPEN: i32 = -4;
const GAP_EXTEND: i32 = -2;
const MATCH: i32 = 2;
const MISMATCH: i32 = -4;
const SEED: u64 = 4;
const CONSENSUS_METHOD: u8 = 1; //0==average 1==median //2==mode
const ERROR_PROBABILITY: f64 = 0.90;
const HOMOPOLYMER_DEBUG: bool = false;
const HOMOPOLYMER: bool = false;
const QUALITY_SCORE: bool = true;
const NUM_OF_ITER_FOR_PARALLEL: usize = 10;
const NUM_OF_ITER_FOR_ZOOMED_GRAPHS: usize = 4;
const USEPACBIODATA: bool = true;
const REVERSE_COMPLEMENT: bool = true;

// file names
const FILENAME: &str = "./data/PacBioReads/155060338.fasta";
const CONSENSUS_FILENAME: &str = "./data/PacBioConsensus/155060338.fastq";
const DEBUG_FILE: &str = "./results/debug.txt";

fn main() {
    let mut seqvec;
    //read the pacbio consensus
    if USEPACBIODATA {
        let pac_bio_consensus =  get_consensus_from_file(CONSENSUS_FILENAME);
        seqvec = vec![pac_bio_consensus];
        seqvec = [seqvec, get_fasta_sequences_from_file(FILENAME)].concat();
        //println!("{}", get_consensus_from_file(CONSENSUS_FILENAME));
    }
    else {
        seqvec = get_random_sequences_from_generator(2000, 10);
    }
    //check_the_alignment_pacbio(seqvec);
    run(seqvec);
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

    if HOMOPOLYMER {
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
        (homopolymer_consensus, homopolymer_topology) = aligner.poa.consensus();
        //get graph
        let homopolymer_graph: &Graph<u8, i32, Directed, usize> = aligner.graph();
        //use homopolymer compressions sequences to make expanded consensus
        let (homopolymer_consensus_freq, homopolymer_score) = get_aligned_homopolymer_sequences_to_homopolymer_consensus(&homopolymer_vec, &homopolymer_consensus);
        let (expanded_consensus, homopolymer_expanded) =  calculate_and_get_expansion (&homopolymer_vec, &homopolymer_consensus, &homopolymer_consensus_freq);
        //get the scores of expanded consensus compared to sequences
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
            write_scores_result_file("./results/results.txt", normal_score, homopolymer_score, expanded_score);
            //modify the graphs to indicate 
            saved_indices = modify_and_write_the_graphs_and_get_zoomed_graphs("./results/normal_graph.fa", "./results/homopolymer_graph.fa", saved_indices, normal_graph, homopolymer_graph);
            write_alignment_and_zoomed_graphs_fasta_file("./results/consensus.fa", &normal_rep, &expanded_rep,
                &count_rep, seqnum as usize, &saved_indices);
        }
    }
    /////////////////////////////
    //quality score calculation//
    /////////////////////////////
    if QUALITY_SCORE {
        let mut invalid_info: Vec<(usize, usize, bool, bool, bool)> = vec![]; //index, node_id, parallel invalid, match invalid, quality invalid
        // calculate and get the quality scores
        let (quality_scores, parallel_validity, base_count_vec, mut debug_strings) = get_consensus_quality_scores(seqnum as usize, &normal_consensus, &normal_topology, normal_graph);
        if USEPACBIODATA {
            let mut parallel_count = 0;
            let mut mismatch_count = 0;
            let mut quality_count = 0;
            let (pacbio_quality_scores, mismatch_indices, aligned_pacbio_bases) = get_quality_score_aligned (get_consensus_from_file(CONSENSUS_FILENAME), &normal_consensus, get_quality_from_file(CONSENSUS_FILENAME));
            // get invalid indices is quality score is too low
            for index in 0..quality_scores.len() {
                if parallel_validity[index] == true {
                    invalid_info.push((index, normal_topology[index], true, false, false));
                    parallel_count += 1;
                }
                if mismatch_indices.contains(&index) {
                    match invalid_info.iter().position(|&r| r.0 == index) {
                        Some(position) => {invalid_info[position].3 = true;},
                        None => {invalid_info.push((index, normal_topology[index], false, true, false));},
                    }
                    mismatch_count += 1;
                }
                if quality_scores[index] < 30.0 {
                    match invalid_info.iter().position(|&r| r.0 == index) {
                        Some(position) => {invalid_info[position].4 = true;},
                        None => {invalid_info.push((index, normal_topology[index], false, false, true));},
                    }
                    quality_count += 1;
                }
            }
            println!("INVALID COUNT: {} parallel err: {} mismatch err: {} quality err: {}", invalid_info.len(), parallel_count, mismatch_count, quality_count);
            debug_strings.push(format!("INVALID COUNT: {} parallel err: {} mismatch err: {} quality err: {}", invalid_info.len(), parallel_count, mismatch_count, quality_count));
            // write the zoomed in graphs for invalid and low quality entries.
            write_quality_scores_to_file("./results/quality_scores.fa", &quality_scores, &normal_consensus, &normal_topology, &invalid_info, &base_count_vec, &pacbio_quality_scores, &aligned_pacbio_bases);
            write_zoomed_quality_score_graphs ("./results/quality_graphs.fa", &invalid_info, &quality_scores, &base_count_vec, normal_graph, &pacbio_quality_scores);
            write_debug_data_to_file(DEBUG_FILE, debug_strings);
        }
        else {
            for index in 0..quality_scores.len() {
                if parallel_validity[index] == true {
                    invalid_info.push((index, normal_topology[index], true, false, false));
                }
            }
            println!("INVALID COUNT: {}", invalid_info.len());
            write_quality_scores_to_file("./results/quality_scores.fa", &quality_scores, &normal_consensus, &normal_topology, &invalid_info, &base_count_vec, &vec![34, 34], &vec![65; quality_scores.len()]);
            write_zoomed_quality_score_graphs ("./results/quality_graphs.fa", &invalid_info, &quality_scores, &base_count_vec, normal_graph, &vec![34, 34]);
            write_quality_score_graph("./results/normal_graph.fa", normal_graph);
            write_debug_data_to_file(DEBUG_FILE, debug_strings);
        }
    }
}

fn write_debug_data_to_file (filename: impl AsRef<Path>, debug_strings: Vec<String>) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n DEBUG",
        chrono::offset::Local::now(), FILENAME)
        .expect("result file cannot be written");
    for line in debug_strings {
        writeln!(file,
            "{}",
            line)
            .expect("result file cannot be written");
    }
}

fn get_quality_score_aligned (pacbio_consensus: String, calculated_consensus: &Vec<u8>, pacbio_quality_scores: String) -> (Vec<usize>, Vec<usize>, Vec<u8>) {
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
                aligned_pacbio_bases.push(126);
                aligned_pacbio_scores_vec.push(33);
                pacbio_index += 1;
            },
            bio::alignment::AlignmentOperation::Ins => {
                calc_index += 1;
            },
            _ => {},
        }
    }
    println!("pacbio index {} pacbio length {}", pacbio_index, pacbio_consensus_vec.len());
    (aligned_pacbio_scores_vec, consensus_match_invalid_indices, aligned_pacbio_bases)
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
                println!("{}", temp_string);
                debug_strings.push(temp_string.clone());
            }
            else {
                temp_string = format!("Going forward");
                println!("{}", temp_string);
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
    println!("{}", temp_string);
    debug_strings.push(temp_string.clone());
    
    // get the two slices of topologically_ordered_list back front
    let mut slice: Vec<Vec<usize>> = [topologically_ordered_nodes[0..target_node_position].to_vec(), topologically_ordered_nodes[target_node_position + 1..topologically_ordered_nodes.len()].to_vec()].to_vec();
    // for debugging
    if slice[0].len() > 10 {
        temp_string = format!("Back slice {:?}", slice[0][(slice[0].len() - 10)..slice[0].len()].to_vec());
        println!("{}", temp_string);
        debug_strings.push(temp_string.clone());
    }
    else {
        temp_string = format!("Back slice {:?}", slice[0][0..slice[0].len()].to_vec());
        println!("{}", temp_string);
        debug_strings.push(temp_string.clone());
    }
    if slice[1].len() > 10 {
        temp_string = format!("Front slice {:?}", slice[1][0..10].to_vec());
        println!("{}", temp_string);
        debug_strings.push(temp_string.clone());
    }
    else {
        temp_string = format!("Front slice {:?}", slice[1][0..slice[1].len()].to_vec());
        println!("{}", temp_string);
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
                println!("{}", temp_string);
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
    let mut seqvec2: Vec<String> = vec![];
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
    if REVERSE_COMPLEMENT {
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

//write stuff here 

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
        /* 
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
        */
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
        chrono::offset::Local::now(), FILENAME, normal_dot)
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

fn write_quality_scores_to_file (filename: impl AsRef<Path>, quality_scores: &Vec<f64>, consensus: &Vec<u8>, topology: &Vec<usize>, invalid_info: &Vec<(usize, usize, bool, bool, bool)>, base_count_vec: &Vec<Vec<usize>>, pacbioquality: &Vec<usize>, aligned_pacbio_bases: &Vec<u8>) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n>Quality score data & graphs:",
        chrono::offset::Local::now(), FILENAME)
        .expect("result file cannot be written");
    for index in 0..consensus.len() {
        let mut print_info = (false, false, false);
        match invalid_info.iter().position(|&r| r.0 == index) {
            Some(position) => {print_info = (invalid_info[position].2, invalid_info[position].3, invalid_info[position].4)},
            None => {},
        }
        if pacbioquality[index % pacbioquality.len()] != 33 {
            writeln!(file,
                "{} {} [{:>6}]\t\t -> {:>8.3}[{}] \tinvalidity =[p,m,q] {:?}\t base_counts = ACGT{:?}",
                consensus[index] as char, aligned_pacbio_bases[index] as char, topology[index], quality_scores[index], (pacbioquality[index % pacbioquality.len()] - 33), print_info, base_count_vec[index])
                .expect("result file cannot be written");
        }
    }
}

fn write_zoomed_quality_score_graphs (filename: impl AsRef<Path>, invalid_info: &Vec<(usize, usize, bool, bool, bool)>, quality_scores: &Vec<f64>, base_count_vec: &Vec<Vec<usize>>, graph: &Graph<u8, i32, Directed, usize>, pacbioquality: &Vec<usize>) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n>Quality score data & graphs:",
        chrono::offset::Local::now(), FILENAME)
        .expect("result file cannot be written");
    for entry in invalid_info {
        let graph_section = get_zoomed_graph_section (graph, &entry.1, &0, 3);
        writeln!(file,
            "\nnode_index:{}\t\tnode_base:{}\t\tquality_score:{:.3}[{}] invalidity:[p,m,q]{:?}\tbase_count:ACGT{:?}\t\n{}\n",
            entry.1, graph.raw_nodes()[entry.1].weight as char, quality_scores[entry.0], (pacbioquality[entry.0 % pacbioquality.len()] - 33), (entry.2, entry.3, entry.4), base_count_vec[entry.0], graph_section)
            .expect("result file cannot be written");
    }
}