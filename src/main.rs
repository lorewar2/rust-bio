use bio::alignment::pairwise::Scoring;
use bio::alignment::{poa::*, TextSlice}; //bandedpoa/ poa
use bio::io::fasta::Index;
use ndarray::indices;
use petgraph::visit::{IntoNeighborsDirected, IntoEdgesDirected};
use std::{
    fs::File,
    fs::OpenOptions,
    io::{prelude::*, BufReader},
    path::Path,
};
use chrono;
use rand::{Rng,SeedableRng, seq};
use rand::rngs::StdRng;
use petgraph::dot::{Dot, Config};
use petgraph::{Directed, Graph, Incoming, Outgoing};
use std::collections::HashMap;
use petgraph::graph::NodeIndex;
use std::{fmt, cmp};
use statrs::function::factorial::binomial;
use num_traits::pow;
use logaddexp::LogAddExp;
use libm::exp;
use petgraph::visit::Topo;

const GAP_OPEN: i32 = -4;
const GAP_EXTEND: i32 = -2;
const MATCH: i32 = 2;
const MISMATCH: i32 = -4;
const FILENAME: &str = "./data/PacBioReads/141232172.fasta";
const CONSENSUS_FILENAME: &str = "./data/PacBioConsensus/141232172.fastq";
const SEED: u64 = 1;
const CONSENSUS_METHOD: u8 = 1; //0==average 1==median //2==mode
const ERROR_PROBABILITY: f64 = 0.90;
const QUALITY_SCORE: bool = true;
const HOMOPOLYMER_DEBUG: bool = false;
const HOMOPOLYMER: bool = false;
const NUM_OF_ITER_FOR_ZOOMED_GRAPHS: usize = 4;
const USEPACBIODATA: bool = true;

fn main() {
    let mut seqvec;
    //read the pacbio consensus
    if USEPACBIODATA {
        let pac_bio_consensus =  get_consensus_from_file(CONSENSUS_FILENAME);
        seqvec = vec![pac_bio_consensus];
        seqvec = [seqvec, get_fasta_sequences_from_file(FILENAME)].concat();
    }
    else {
        seqvec = get_random_sequences_from_generator(1000, 10);
    }
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
        (homopolymer_consensus, homopolymer_topology) = aligner.poa.consensus(); //poa
        //get graph
        let homopolymer_graph: &Graph<u8, i32, Directed, usize> = aligner.graph();
        //use homopolymer compressions sequences to make expanded consensus
        let (expanded_consensus, homopolymer_consensus_freq, homopolymer_score, homopolymer_expanded) =  get_expanded_consensus(homopolymer_vec, &homopolymer_consensus);
        //get the scores of expanded consensus compared to sequences
        let expanded_score = get_consensus_score(&seqvec, &expanded_consensus);

        if HOMOPOLYMER_DEBUG {
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
    /////////////////////////////
    //quality score calculation//
    /////////////////////////////
    let mut invalid_indices: Vec<(usize, usize)> = vec![];
    // calculate and get the quality scores
    let (quality_scores, validity, base_count_vec) = get_consensus_quality_scores(seqnum as usize, &normal_consensus, &normal_topology, normal_graph);
    // get invalid indices is quality score is too low
    for index in 0..quality_scores.len() {
        if //validity[index] == true &&
            quality_scores[index] <= 30.00 {
            invalid_indices.push((index, normal_topology[index]));
        }
    }
    println!("ERROR COUNT: {}", invalid_indices.len());
    //get the pacbio quality
    if USEPACBIODATA {
        // write the zoomed in graphs for invalid and low quality entries.
        let pacbio_quality_scores = get_quality_from_file(CONSENSUS_FILENAME);
        write_quality_scores_to_file("./results/quality_scores.fa", &quality_scores, &normal_consensus, &normal_topology, &validity, &base_count_vec, &pacbio_quality_scores);
        write_zoomed_quality_score_graphs ("./results/quality_graphs.fa", &invalid_indices, &quality_scores, &base_count_vec, normal_graph, &pacbio_quality_scores);
    }
    else {
        write_quality_scores_to_file("./results/quality_scores.fa", &quality_scores, &normal_consensus, &normal_topology, &validity, &base_count_vec, &" ".to_string());
        write_zoomed_quality_score_graphs ("./results/quality_graphs.fa", &invalid_indices, &quality_scores, &base_count_vec, normal_graph, &" ".to_string());
        write_quality_score_graph("./results/normal_graph.fa", normal_graph);
    }
}



fn get_consensus_quality_scores(mut seq_num: usize, consensus: &Vec<u8>, topology: &Vec<usize>, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<f64>, Vec<bool>, Vec<Vec<usize>>) {
    let mut quality_scores: Vec<f64> = vec![];
    let mut validity: Vec<bool> = vec![];
    let mut base_count_vec: Vec<Vec<usize>> = vec![];
    //run all the consensus through get indices
    for i in 0..consensus.len() {
        // skip the indices which are in the passed consensus
        let mut skip_nodes: Vec<usize> = topology[0 .. i + 1].to_vec();
        // new method using topology cut
        let mut target_node_parent = None;
        if i != 0 {
            target_node_parent = Some(topology[i - 1]);
        } 
        let (parallel_nodes, parallel_num_incoming_seq) = get_parallel_nodes_with_topology_cut (seq_num, skip_nodes, topology[i], target_node_parent, graph);
        let (temp_quality_score, temp_count_mismatch, temp_base_counts) = base_quality_score_calculation(seq_num, parallel_nodes, parallel_num_incoming_seq, consensus[i], graph);
        quality_scores.push(temp_quality_score);
        validity.push(temp_count_mismatch);
        base_count_vec.push(temp_base_counts);
    }
    (quality_scores, validity, base_count_vec)
}

fn get_parallel_nodes_with_topology_cut (total_seq: usize, mut skip_nodes: Vec<usize>, target_node: usize, target_node_parent: Option<usize>, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>) {
    // vector initialization
    let mut topology = Topo::new(graph);
    let mut topologically_ordered_nodes = Vec::new();
    let mut parallel_nodes: Vec<usize> = vec![];
    let mut parallel_node_parents: Vec<usize> = vec![];
    let mut parallel_num_incoming_seq: Vec<usize> = vec![];

    // make a topologically ordered list
    while let Some(node) = topology.next(graph) {
        topologically_ordered_nodes.push(node.index());
    }
    // find the position of the target node in topology list
    let target_node_topological_position = topologically_ordered_nodes.iter().position(|&r| r == target_node).unwrap();
    // get a list of back nodes to skip
    skip_nodes = [skip_nodes, get_back_nodes(3, vec![], target_node, graph)].concat();
    

    // check if the target node corrosponds with all the sequences
    let num_seq_through_target_base = find_the_seq_passing_through (target_node, graph);

    if num_seq_through_target_base == total_seq {
        parallel_nodes.push(target_node);
        if USEPACBIODATA {
            parallel_num_incoming_seq.push(num_seq_through_target_base - 1);
        }
        else {
            parallel_num_incoming_seq.push(num_seq_through_target_base);
        }
        return (parallel_nodes, parallel_num_incoming_seq);
    }
    // go back skip_count and go forward skip_count + 3 and check if parent and child are before and after target_node_position,
    // iterate skip_count until all sequences are found, break on 5
    let mut seq_found_so_far = num_seq_through_target_base;
    let mut bubble_size = 1;
    while seq_found_so_far < total_seq  && bubble_size < 5 {
        let mut skip_forward: Vec<usize> = vec![];
        let mut skip_forward_parent: Vec<usize> = vec![];
        (parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far) = check_neighbours_and_find_crossing_nodes (parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, target_node, target_node_parent, bubble_size, &topologically_ordered_nodes, target_node_topological_position, graph);
        //check if the there are any nodes which are not parallel (in sequence) and remove them
        for parallel_node in &parallel_nodes {
            if !skip_nodes.contains(parallel_node) && !skip_forward.contains(parallel_node) {
                skip_forward = [skip_forward, get_forward_nodes(2, vec![], *parallel_node, graph)].concat();
                let temp_skip_forward_parent = vec![*parallel_node; skip_forward.len() - skip_forward_parent.len()];
                skip_forward_parent = [skip_forward_parent, temp_skip_forward_parent].concat();
                //println!("added to skip {:?} {} {}", get_forward_nodes(2, vec![], *parallel_node, graph), skip_forward_parent.len(), skip_forward.len());
            }
        }
        // remove forward skip nodes if parent is the parallel node ** obsolete **
        while skip_forward.iter().any(|&i| parallel_nodes.contains(&i)) {
            let position_parallel = parallel_nodes.iter().position(|&r| skip_forward.contains(&r)).unwrap();
            let position_forward = skip_forward.iter().position(|&i| parallel_nodes.contains(&i)).unwrap();
            if parallel_node_parents[position_parallel] == skip_forward_parent[position_forward] {
                //println!("removing this due to forward: {}", parallel_nodes[position_parallel]);
                parallel_nodes.remove(position_parallel);
                parallel_node_parents.remove(position_parallel);
                seq_found_so_far -= parallel_num_incoming_seq[position_parallel];
                parallel_num_incoming_seq.remove(position_parallel);
            }
            else {
                break;
            }
        }
        // remove the skip nodes if present in parallel nodes ** obsolete **
        while skip_nodes.iter().any(|&i| parallel_nodes.contains(&i)) {
            let position = parallel_nodes.iter().position(|&r| skip_nodes.contains(&r)).unwrap();
            //println!("removing this: {}", parallel_nodes[position]);
            seq_found_so_far -= parallel_num_incoming_seq[position];
            parallel_nodes.remove(position);
            parallel_node_parents.remove(position);
            parallel_num_incoming_seq.remove(position);
        }
        bubble_size += 1;
    }
    if USEPACBIODATA {
        parallel_num_incoming_seq.push(num_seq_through_target_base - 1);
    }
    else {
        parallel_num_incoming_seq.push(num_seq_through_target_base);
    }
    parallel_nodes.push(target_node);
    (parallel_nodes, parallel_num_incoming_seq)
}

fn check_neighbours_and_find_crossing_nodes (mut parallel_nodes: Vec<usize>, mut parallel_node_parents: Vec<usize>, mut parallel_num_incoming_seq: Vec<usize>, mut seq_found_so_far: usize, focus_node: usize, focus_node_parent: Option<usize>, bubble_size: usize, topologically_ordered_nodes: &Vec<usize>, target_node_position: usize, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>, Vec<usize>, usize) {
    // get a list of x bubble edge nodes
    //let edge_nodes_list = find_nth_iteration_neighbouring_indices(bubble_size, vec![], focus_node, graph);
    
    // get a list of x back_iterations back nodes
    let back_nodes_list = get_xiterations_back_nodes(bubble_size, vec![], focus_node, graph);
    // get a list of all forward nodes 0..(back_iterations + 3) for all the back_nodes
    let mut edge_nodes_list: Vec<usize> = vec![];
    for back_node in &back_nodes_list {
        let temp_forward_list = get_forward_nodes (bubble_size + 3, vec![], *back_node, graph);
        for temp_forward_node in &temp_forward_list {
            if !edge_nodes_list.contains(temp_forward_node) {
                edge_nodes_list.push(*temp_forward_node);
            }
        }
    }
    //println!("{} {:?} {:?}", bubble_size, back_nodes_list, edge_nodes_list);
    // get the intermidiate slice between node_position and its parent
    let mut target_node_parent_position = 0;
    let intermediate_slice = match focus_node_parent {
        Some(x) => {
            target_node_parent_position = topologically_ordered_nodes.iter().position(|&r| r == x).unwrap();
            topologically_ordered_nodes[target_node_parent_position..target_node_position].to_vec()
        },
        None => {vec![]},
    };
    // get the two slices of topologically_ordered_list
    let back_slice = topologically_ordered_nodes[0..target_node_parent_position + 1].to_vec();
    let front_slice = topologically_ordered_nodes[target_node_position + 1..topologically_ordered_nodes.len()].to_vec();
    if(back_slice.len() > 5 && front_slice.len() > 5){
        //println!("back slice {:?}\nfront slice {:?}", back_slice[(back_slice.len()-5)..back_slice.len()].to_vec(), front_slice[0..5].to_vec());
    }
    
    //println!("intermediate slice {:?}", intermediate_slice);
    //iterate through edge nodes obtained
    for edge_node in &edge_nodes_list {
        // get the parents of the edge node
        let edge_node_parents = get_back_nodes (1, vec![], *edge_node, graph);
        'parent_loop: for edge_node_parent in &edge_node_parents {
            // if the parent is in back section and node is in front section add to parallel nodes or if both parent and target is in intermediate add to parallel loop
            if (back_slice.contains(edge_node_parent) && front_slice.contains(edge_node) && (*edge_node_parent != focus_node)) ||
                (intermediate_slice.contains(edge_node_parent) && intermediate_slice.contains(edge_node)) ||
                (back_slice.contains(edge_node_parent) && intermediate_slice.contains(edge_node)) {
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
                if get_forward_nodes(4, vec![],  *edge_node, graph).contains(&focus_node) {
                    continue;
                }
                parallel_nodes.push(*edge_node);
                parallel_node_parents.push(*edge_node_parent);
                //print!("success node {} parent {}\n", *edge_node, *edge_node_parent);
                // get the edge weight and add to seq_found_so_far
                let mut incoming_weight = 0;
                let mut edges = graph.edges_connecting(NodeIndex::new(*edge_node_parent), NodeIndex::new(*edge_node));
                while let Some(edge) = edges.next() {
                    incoming_weight += edge.weight().clone();
                }
                parallel_num_incoming_seq.push(incoming_weight as usize);
                seq_found_so_far += incoming_weight as usize;
            }
            
        }
    }
    (parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far)
}

fn find_nth_iteration_neighbouring_indices (num_of_iterations: usize, mut neighbour_nodes: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize> ) -> Vec<usize> {
    if num_of_iterations <= 0 {
        return neighbour_nodes;
    }
    let mut immediate_neighbours = graph.neighbors_undirected(NodeIndex::new(focus_node));
    while let Some(neighbour_node) = immediate_neighbours.next() {
        if num_of_iterations == 1 {
            if !neighbour_nodes.contains(&neighbour_node.index()){
                neighbour_nodes.push(neighbour_node.index());
            }
        }
        neighbour_nodes = find_nth_iteration_neighbouring_indices(num_of_iterations - 1, neighbour_nodes, neighbour_node.index(), graph);   
    }
    neighbour_nodes
}

fn get_forward_nodes (iteration: usize, mut forward_node_list: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    if iteration <= 0 {
        return forward_node_list;
    }
    //get the back nodes of the target
    let mut forward_neighbours = graph.neighbors_directed(NodeIndex::new(focus_node), Outgoing);
    //iterate through the neighbours
    while let Some(forward_neighbour) = forward_neighbours.next() {
        if !forward_node_list.contains(&forward_neighbour.index()){
            forward_node_list.push(forward_neighbour.index());
            forward_node_list = get_forward_nodes (iteration - 1, forward_node_list, forward_neighbour.index(), graph);
        }
    }
    forward_node_list
}

fn get_back_nodes (iteration: usize, mut back_node_list: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    if iteration <= 0 {
        return back_node_list;
    }
    //get the back nodes of the target
    let mut back_neighbours = graph.neighbors_directed(NodeIndex::new(focus_node), Incoming);
    //iterate through the neighbours
    while let Some(back_neighbour) = back_neighbours.next() {
        if !back_node_list.contains(&back_neighbour.index()){
            back_node_list.push(back_neighbour.index());
            back_node_list = get_back_nodes (iteration - 1, back_node_list, back_neighbour.index(), graph);
        }
    }
    back_node_list
}

fn get_xiterations_back_nodes (iteration: usize, mut back_node_list: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    if iteration <= 0 {
        return back_node_list;
    }
    //get the back nodes of the target
    let mut back_neighbours = graph.neighbors_directed(NodeIndex::new(focus_node), Incoming);
    //iterate through the neighbours
    while let Some(back_neighbour) = back_neighbours.next() {
        if iteration == 1 {
            if !back_node_list.contains(&back_neighbour.index()){
                back_node_list.push(back_neighbour.index());
            }
        }
        back_node_list = get_xiterations_back_nodes (iteration - 1, back_node_list, back_neighbour.index(), graph);
    }
    back_node_list
}

fn base_quality_score_calculation (total_seq: usize, indices_of_parallel_nodes: Vec<usize>, seq_through_parallel_nodes: Vec<usize>, base: u8, graph: &Graph<u8, i32, Directed, usize>) -> (f64, bool, Vec<usize>) {
    //variable initialization
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
    if (base_a_count + base_c_count + base_g_count + base_t_count) != (total_seq - 1) {
        count_mismatch = true;
    }
    match count_mismatch {
        true => {
            //println!("base counts A:{} C:{} G:{} T:{} MISMATCHHHH!!!!!!!!!!!!!!!!!!!!!", base_a_count, base_c_count, base_g_count, base_t_count);
        },
        false => {
            //println!("base counts A:{} C:{} G:{} T:{}", base_a_count, base_c_count, base_g_count, base_t_count);
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
    (quality_score, count_mismatch, base_counts)
}

fn calculate_binomial (n: usize, k: usize, prob: f64) -> f64 {
    let binomial_coeff = binomial(n as u64, k as u64);
    let success: f64 = prob.powf(k as f64);
    let failure: f64 = (1.00 - prob).powf((n - k) as f64);
    binomial_coeff * success * failure
}

fn get_indices_of_parallel_nodes_of_target (mut skip_nodes: Vec<usize>, total_seq: usize, target_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>, usize) {
    // initialize the vectors
    let mut parallel_nodes = vec![target_node];
    let mut parallel_num_incoming_seq = vec![];
    // find out how many sequences run through the target_node in graph
    let num_seq_through_target_base = find_the_seq_passing_through (target_node, graph);
    parallel_num_incoming_seq.push(num_seq_through_target_base);
    if num_seq_through_target_base == total_seq {
        // nothing to do return the num_seq_corrosponding to other bases as 0
        (parallel_nodes, parallel_num_incoming_seq, num_seq_through_target_base)
    }
    else {
        // populate a list of back nodes and add them to the skip list
        // iterate until all the parallel nodes are found, while skipping the amount you go back 
        skip_nodes = [skip_nodes, get_back_nodes(3, vec![], target_node, graph)].concat();
        //println!("SKIP NODES {:?}", skip_nodes);
        let mut seq_found_so_far = num_seq_through_target_base;
        let mut prev_back_nodes: Vec<usize> = vec![target_node];
        let mut skip_nodes_internal: Vec<usize> = vec![];
        let mut skip_nodes_internal_parents: Vec<usize> = vec![];
        let mut skip_by_index = 1;
        while seq_found_so_far < total_seq  && skip_by_index < 5{
            (parallel_nodes, parallel_num_incoming_seq, prev_back_nodes, seq_found_so_far, skip_nodes, skip_nodes_internal, skip_nodes_internal_parents) = go_back_one_and_get_the_target_indices(parallel_nodes, parallel_num_incoming_seq, &prev_back_nodes, skip_nodes, graph, seq_found_so_far, skip_by_index, skip_nodes_internal, skip_nodes_internal_parents);
            skip_by_index += 1;
            // remove the skip nodes if present in parallel nodes ** obsolete **
            while skip_nodes[1..skip_nodes.len()].iter().any(|&i| parallel_nodes.contains(&i)) {
                let position = parallel_nodes.iter().position(|&r| skip_nodes.contains(&r)).unwrap();
                parallel_nodes.remove(position);
            }
        }
        (parallel_nodes, parallel_num_incoming_seq, seq_found_so_far)
    }
}

fn go_back_one_and_get_the_target_indices (mut parallel_nodes: Vec<usize>, mut parallel_num_incoming_seq: Vec<usize>, prev_back_nodes: &Vec<usize>, 
                                            skip_nodes: Vec<usize>, graph: &Graph<u8, i32, Directed, usize>, mut seq_found_so_far: usize, 
                                            skip_by_index: usize, mut skip_nodes_internal: Vec<usize>, mut skip_nodes_internal_parents: Vec<usize>) 
                                            -> (Vec<usize>, Vec<usize>, Vec<usize>, usize, Vec<usize>, Vec<usize>, Vec<usize>) {
    // initialize the vectors
    print!("iteration: {} prev back nodes: {:?}", skip_by_index, prev_back_nodes);
    let mut current_back_nodes: Vec<usize> = vec![];
    // populate a list of nodes which are one edge behind the prev_back_nodes
    for prev_back_node in prev_back_nodes {
        let node_index = NodeIndex::new(*prev_back_node);
        let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(node_index, Incoming).collect();
        for incoming_node in incoming_nodes {
            current_back_nodes.push(incoming_node.index());
        }
    }
    print!(" current back nodes: {:?}\n", current_back_nodes);
    // go through the list while finding parallel nodes break
    for current_back_node in &current_back_nodes {
        //get a list of nodes which are skip_by_index edges away from the back node and not in the nodes_travelled
        let (acquired_parallel_nodes, acquired_parallel_num_seq_incoming, acquired_parallel_parent_nodes) = find_most_front_neighbour_indices(skip_by_index, current_back_node, graph);
        println!("acquired nodes {:?}",  acquired_parallel_nodes);
        for index in 0..acquired_parallel_nodes.len() {
            //check if the skip nodes parents 
            if !skip_nodes.contains(&acquired_parallel_nodes[index]) {  
                if skip_nodes_internal.contains(&acquired_parallel_nodes[index]) {
                    //get the position of the skip node internal 
                    let position = skip_nodes_internal.iter().position(|&r| r == acquired_parallel_nodes[index]).unwrap();
                    // check if the parent of the acquired parallel node is the same as skip nodes parent, skip if yes, else put it in
                    if skip_nodes_internal_parents[position] == acquired_parallel_parent_nodes[index] {
                        continue;
                    }
                }
                parallel_nodes.push(acquired_parallel_nodes[index]);
                parallel_num_incoming_seq.push(acquired_parallel_num_seq_incoming[index]);
                skip_nodes_internal.push(acquired_parallel_nodes[index]);
                skip_nodes_internal_parents.push(acquired_parallel_parent_nodes[index]);
                seq_found_so_far += acquired_parallel_num_seq_incoming[index];
            }
        }
    }
    (parallel_nodes, parallel_num_incoming_seq, current_back_nodes, seq_found_so_far, skip_nodes, skip_nodes_internal, skip_nodes_internal_parents)
}

fn find_most_front_neighbour_indices (num_of_iterations: usize, focus_node: &usize, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    // initialize required vectors
    let mut parallel_indices: Vec<usize> = vec![];
    let mut parallel_num_seq_incoming: Vec<usize> = vec![];
    let mut parallel_parent_nodes: Vec<usize> = vec![];
    // break when the recursion reach target
    if num_of_iterations <= 0 {
        return (parallel_indices, parallel_num_seq_incoming, parallel_parent_nodes);
    }
    // get the nodes with a directed edge from the focus node
    let mut front_neighbours = graph.neighbors_directed(NodeIndex::new(*focus_node), Outgoing);
    //iterate through the neighbours
    while let Some(front_neighbour) = front_neighbours.next() {
        // save the unique final node index to indices vector and edge weight to seq_incoming 
        if num_of_iterations == 1 {
            parallel_parent_nodes.push(*focus_node);
            parallel_indices.push(front_neighbour.index());
            // get the weight
            let mut incoming_weight = 0;
            let mut edges = graph.edges_connecting(NodeIndex::new(*focus_node), front_neighbour);
            while let Some(edge) = edges.next() {
                incoming_weight += edge.weight().clone();
            }
            parallel_num_seq_incoming.push(incoming_weight as usize);
        }
        // iterate through neighbours of neighours
        let (obtained_parallel_indices, obtained_parallel_num_seq_incoming, obtained_parallel_parent_nodes)  = find_most_front_neighbour_indices(num_of_iterations - 1, &front_neighbour.index(), graph);
        // save the obtained values
        if num_of_iterations - 1 == 1 {
            parallel_indices = [parallel_indices, obtained_parallel_indices].concat();
            parallel_num_seq_incoming = [parallel_num_seq_incoming, obtained_parallel_num_seq_incoming].concat();
            parallel_parent_nodes = [parallel_parent_nodes, obtained_parallel_parent_nodes].concat();
        }
    }
    (parallel_indices, parallel_num_seq_incoming, parallel_parent_nodes)
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

fn write_quality_score_graph (normal_filename: impl AsRef<Path>, normal_graph: &Graph<u8, i32, Directed, usize>) {
    let mut count = 0;
    let mut normal_dot = format!("{:?}", Dot::new(&normal_graph.map(|_, n| (*n) as char, |_, e| *e)));
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

fn write_quality_scores_to_file (filename: impl AsRef<Path>, quality_scores: &Vec<f64>, consensus: &Vec<u8>, topology: &Vec<usize>, validity: &Vec<bool>, base_count_vec: &Vec<Vec<usize>>, pacbioquality: &String) {
    let pacbiochar: Vec<char> = pacbioquality.chars().collect();
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    for index in 0..consensus.len() {
        writeln!(file,
            "{}[{:>6}]\t\t -> {:>8.3}[{}] \tvalid = {}\t base_counts = ACGT{:?}",
            consensus[index] as char, topology[index], quality_scores[index], (pacbiochar[index % pacbiochar.len()] as u8 - 33), !validity[index], base_count_vec[index])
            .expect("result file cannot be written");
    }
}

fn write_zoomed_quality_score_graphs (filename: impl AsRef<Path>, write_indices: &Vec<(usize, usize)>, quality_scores: &Vec<f64>, base_count_vec: &Vec<Vec<usize>>, graph: &Graph<u8, i32, Directed, usize>, pacbioquality: &String) {
    let pacbiochar: Vec<char> = pacbioquality.chars().collect();
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(filename)
        .unwrap();
    writeln!(file,
        "{:?}\nFILE: {}\n>Quality score data & graphs:",
        chrono::offset::Local::now(), FILENAME)
        .expect("result file cannot be written");
    for index in write_indices {
        let graph_section = get_zoomed_graph_section (graph, &index.1, &0, 3);
        writeln!(file,
            "\nnode_index:{}\t\tnode_base:{}\t\tquality_score:{:.3}[{}]\tbase_count:ACGT{:?}\t\n{}\n",
            index.1, graph.raw_nodes()[index.1].weight as char, quality_scores[index.0], (pacbiochar[index.0 % pacbiochar.len()] as u8 - 33), base_count_vec[index.0], graph_section)
            .expect("result file cannot be written");
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

