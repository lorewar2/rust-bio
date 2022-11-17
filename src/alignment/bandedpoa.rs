
use std::cmp::{max, Ordering};

use crate::utils::TextSlice;

use crate::alignment::pairwise::{MatchFunc, Scoring};

use petgraph::Direction::Outgoing;
use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;

use petgraph::{Directed, Graph, Incoming};

pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs
pub type POAGraph = Graph<u8, i32, Directed, usize>;

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
#[derive(Debug, Clone)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
}


pub struct Alignment {
    pub score: i32,
    //    xstart: Edge,
    pub operations: Vec<AlignmentOperation>,
}

#[derive(Debug, Clone)]
pub struct TracebackCell {
    score: i32,
    op: AlignmentOperation,
}

impl Ord for TracebackCell {
    fn cmp(&self, other: &TracebackCell) -> Ordering {
        self.score.cmp(&other.score)
    }
}

impl PartialOrd for TracebackCell {
    fn partial_cmp(&self, other: &TracebackCell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TracebackCell {
    fn eq(&self, other: &TracebackCell) -> bool {
        self.score == other.score
    }
}

//impl Default for TracebackCell { }

impl Eq for TracebackCell {}

pub struct Traceback {
    rows: usize,
    cols: usize,

    // store the last visited node in topological order so that
    // we can index into the end of the alignment when we backtrack
    last: NodeIndex<usize>,
    matrix: Vec<Vec<TracebackCell>>,
}

impl Traceback {
    /// Create a Traceback matrix with given maximum sizes
    ///
    /// # Arguments
    ///
    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    fn with_capacity(m: usize, n: usize) -> Self {
        let matrix = vec![
            vec![
                TracebackCell {
                    score: 0,
                    op: AlignmentOperation::Match(None)
                };
                n + 1
            ];
            m + 1
        ];
        Traceback {
            rows: m,
            cols: n,
            last: NodeIndex::new(0),
            matrix,
        }
    }

    /// Populate the edges of the traceback matrix
    fn initialize_scores(&mut self, gap_open: i32) {
        for (i, row) in self
            .matrix
            .iter_mut()
            .enumerate()
            .take(self.rows + 1)
            .skip(1)
        {
            // TODO: these should be -1 * distance from head node
            row[0] = TracebackCell {
                score: (i as i32) * gap_open, // gap_open penalty
                op: AlignmentOperation::Del(None),
            };
        }
        for j in 1..=self.cols {
            self.matrix[0][j] = TracebackCell {
                score: (j as i32) * gap_open,
                op: AlignmentOperation::Ins(None),
            };
        }
    }

    fn new() -> Self {
        Traceback {
            rows: 0,
            cols: 0,
            last: NodeIndex::new(0),
            matrix: Vec::new(),
        }
    }

    fn set(&mut self, i: usize, j: usize, cell: TracebackCell) {
        self.matrix[i][j] = cell;
    }

    fn get(&self, i: usize, j: usize) -> &TracebackCell {
        &self.matrix[i][j]
    }

    pub fn print(&self, g: &Graph<u8, i32, Directed, usize>, query: TextSlice) {
        let (m, n) = (g.node_count(), query.len());
        print!(".\t");
        for base in query.iter().take(n) {
            print!("{:?}\t", *base);
        }
        for i in 0..m {
            print!("\n{:?}\t", g.raw_nodes()[i].weight);
            for j in 0..n {
                print!("{}.\t", self.get(i + 1, j + 1).score);
            }
        }
        println!();
    }

    pub fn print_operation(&self, g: &Graph<u8, i32, Directed, usize>, query: TextSlice) {
        let (m, n) = (g.node_count(), query.len());
        print!(".\t");
        for base in query.iter().take(n) {
            print!("{:?}\t", *base);
        }
        for i in 0..m {
            print!("\n{:?}\t", g.raw_nodes()[i].weight);
            for j in 0..n {
                let operation = &self.get(i + 1, j + 1).op;
                match operation {
                    AlignmentOperation::Match(x) => match x {
                        Some(x) => print!("m{},{}\t",x.0,x.1),
                        None => print!("mn\t"),
                    },
                    AlignmentOperation::Del(x) => match x {
                        Some(x) => print!("d{},{}\t",x.0,x.1),
                        None => print!("dn\t"),
                    },
                    AlignmentOperation::Ins(x) => match x {
                        Some(x) => print!("i{}\t",x),
                        None => print!("in\t"),
                    }
                }
            }
        }
        println!();
    }

    pub fn alignment(&self) -> Alignment {
        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];

        // Now backtrack through the matrix to construct an optimal path
        let mut i = self.last.index() + 1;
        let mut j = self.cols;

        while i > 0 || j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.matrix[i][j].op.clone());
            match self.matrix[i][j].op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Del(Some((p, _))) => {
                    i = p + 1;
                }
                AlignmentOperation::Ins(Some(p)) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Match(None) => {
                    i -= 1; // break;
                    j -= 1;
                }
                AlignmentOperation::Del(None) => {
                    i -= 1; // j -= 1;
                }
                AlignmentOperation::Ins(None) => {
                    j -= 1; // i -= 1;
                }
            }
        }

        ops.reverse();

        Alignment {
            score: self.matrix[self.last.index() + 1][self.cols].score,
            operations: ops,
        }
    }
}

/// A partially ordered aligner builder
///
/// Uses consuming builder pattern for constructing partial order alignments with method chaining
pub struct BandedAligner<F: MatchFunc> {
    traceback: Traceback,
    query: Vec<u8>,
    pub poa: Poa<F>,
    consensus: Vec<u8>,
    band: Band,
    k: usize,
    w: usize,
}

impl<F: MatchFunc> Aligner<F> {
    /// Create new instance.
    pub fn new(scoring: Scoring<F>, reference: TextSlice) -> Self {
        println!("POA Aligner constructor run"); //added
        //println!("Reference to vector: {:?}", reference.to_vec()); //added
        Aligner {
            traceback: Traceback::new(),
            query: reference.to_vec(),
            poa: Poa::from_string(scoring, reference),
        }
    }

    /// Add the alignment of the last query to the graph.
    pub fn add_to_graph(&mut self, seq_number: u8) -> &mut Self {
        let alignment = self.traceback.alignment();
        self.poa.add_alignment(&alignment, &self.query, seq_number);
        //println!("{:?}", Dot::with_config(&self.poa.graph, &[Config::EdgeNoLabel])); //added
        self
    }

    /// Return alignment of last added query against the graph.
    pub fn alignment(&self) -> Alignment {
        self.traceback.alignment()
    }

    /// Globally align a given query against the graph.
    pub fn global(&mut self, query: TextSlice) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.global(query);
        //self.traceback.print(&self.poa.graph, query);
        //self.traceback.print_operation(&self.poa.graph, query);
        //println!(" ");
        /* 
        for op in self.traceback.alignment().operations{
            match op {
                AlignmentOperation::Match(x) => match x {
                    Some(x) => print!("m{},{}\t",x.0,x.1),
                    None => print!("mn\t"),
                },
                AlignmentOperation::Del(x) => match x {
                    Some(x) => print!("d{},{}\t",x.0,x.1),
                    None => print!("dn\t"),
                },
                AlignmentOperation::Ins(x) => match x {
                    Some(x) => print!("i{}\t",x),
                    None => print!("in\t"),
                }
            }
        }
        */
        self
    }

    /// Return alignment graph.
    pub fn graph(&self) -> &POAGraph {
        //println!("Reference Graph: {:?}", Dot::with_config(&self.poa.graph, &[Config::EdgeIndexLabel])); //added
        &self.poa.graph
    }
}

/// A partially ordered alignment graph
///
/// A directed acyclic graph datastructure that represents the topology of a
/// traceback matrix.
pub struct Poa<F: MatchFunc> {
    scoring: Scoring<F>,
    pub graph: POAGraph,
    pub node_seq_tracker: Vec<Vec<u8>>,
    pub consensus: Vec<u8>,
}

impl<F: MatchFunc> Poa<F> {
    /// Create a new aligner instance from the directed acyclic graph of another.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `poa` - the partially ordered reference alignment
    pub fn new(scoring: Scoring<F>, graph: POAGraph) -> Self {
        Poa { scoring, graph, node_seq_tracker: Vec::new() }
    }

    /// Create a new POA graph from an initial reference sequence and alignment penalties.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `reference` - a reference TextSlice to populate the initial reference graph
    pub fn from_string(scoring: Scoring<F>, seq: TextSlice) -> Self {
        let mut node_seq_tracker = Vec::new();
        //println!("from_string function to make partial order graph-petgraph lib"); //added
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        node_seq_tracker.push(vec![0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node_seq_tracker.push(vec![0]);
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        //println!("Reference Graph: {:?}", Dot::with_config(&graph, &[Config::EdgeIndexLabel])); //added
        Poa { scoring, graph, node_seq_tracker }
    }

    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    pub fn global(&self, query: TextSlice) -> Traceback {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.initialize_scores(self.scoring.gap_open);
        //println!("Printing the  empty Initialized matrix");//added
        //traceback.print(&self.graph, query);//added

        traceback.set(
            0,
            0,
            TracebackCell {
                score: 0,
                op: AlignmentOperation::Match(None),
            },
        );
        // construct the score matrix (O(n^2) space)
        //println!("Topological sort of reference graph!!"); //added
        let mut topo = Topo::new(&self.graph);
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            //println!("Previous Index Reference Node being processed index:{} base:{}", node.index(), r); //added
            let i = node.index() + 1;
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            //println!("Nodes with directed edges to current node:{:?}", prevs); //added
            // query base and its index in the DAG (traceback matrix rows)
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    //println!("NO directed, matching score from i = 0 and j = {} : {} + {}", j - 1 ,traceback.get(0, j - 1).score , self.scoring.match_fn.score(r, *q));
                    TracebackCell {
                        score: traceback.get(0, j - 1).score + self.scoring.match_fn.score(r, *q),
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = TracebackCell {
                        score: MIN_SCORE,
                        op: AlignmentOperation::Match(None),
                    };
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        max_cell = max(
                            max_cell,
                            max(
                                TracebackCell {
                                    score: traceback.get(i_p, j - 1).score
                                        + self.scoring.match_fn.score(r, *q),
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: traceback.get(i_p, j).score + self.scoring.gap_open,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };

                let score = max(
                    max_cell,
                    TracebackCell {
                        score: traceback.get(i, j - 1).score + self.scoring.gap_open,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                //println!("{}", score.score);
                traceback.set(i, j, score);
            }
        }

        traceback
    }
    
    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    pub fn add_alignment(&mut self, aln: &Alignment, seq: TextSlice, seq_number: u8) {
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut i: usize = 0;
        let mut start_seq_unmatched: bool = false;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => { //previously i += 1;
                    //println!("");
                    let node: NodeIndex<usize> = NodeIndex::new(0);
                    if (seq[i] != self.graph.raw_nodes()[0].weight) && (seq[i] != b'X') {
                        //println!("Start node mismatch with sequence, do something");
                        let new_node = self.graph.add_node(seq[i]);
                        self.node_seq_tracker.push(vec![seq_number]);
                        if start_seq_unmatched == true {
                            //println!("matchn making edge from {}->{}", seq[i], self.graph.raw_nodes()[prev.index()].weight);
                            self.graph.add_edge(prev, new_node, 1);
                            start_seq_unmatched = false;
                        }
                        prev = new_node;
                    }
                    else {
                        //println!("Start node match with sequence, do nothing");
                        if start_seq_unmatched == true {
                            //println!("matchn making edge from {}->{}", seq[i], self.graph.raw_nodes()[prev.index()].weight);
                            self.graph.add_edge(prev, node, 1);
                            prev = node;
                            start_seq_unmatched = false;
                        }
                    }
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {   
                    let node = NodeIndex::new(*p);
                    if (seq[i] != self.graph.raw_nodes()[*p].weight) && (seq[i] != b'X') {
                        //println!("mpx making edge from {}->{}",self.graph.raw_nodes()[prev.index()].weight, seq[i]);
                        let node = self.graph.add_node(seq[i]);
                        self.node_seq_tracker.push(vec![seq_number]);
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                //println!("mpm incr edge from {}->{}",self.graph.raw_nodes()[prev.index()].weight, seq[i]);
                                *self.graph.edge_weight_mut(edge).unwrap() += 1;
                                self.node_seq_tracker[node.index()].push(seq_number);
                            }
                            None => {
                                // where the previous node was newly added
                                //println!("mpm making edge from {}->{}",self.graph.raw_nodes()[prev.index()].weight, seq[i]);
                                self.graph.add_edge(prev, node, 1);
                                self.node_seq_tracker[node.index()].push(seq_number);
                            }
                        }
                        prev = NodeIndex::new(*p);
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => { // previously just i += 1
                    //insertion at the start, take care if the sequence
                    let node = self.graph.add_node(seq[i]);
                    self.node_seq_tracker.push(vec![seq_number]);
                    if start_seq_unmatched == true {
                        //println!("insn making edge from {}->{}", seq[i], self.graph.raw_nodes()[prev.index()].weight);
                        self.graph.add_edge(prev, node, 1);
                    }
                    prev = node;
                    start_seq_unmatched = true;
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    if start_seq_unmatched == true {
                        let node = NodeIndex::new(0);
                        //println!("insn making edge from {}->{}", seq[i], self.graph.raw_nodes()[prev.index()].weight);
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                        start_seq_unmatched = false;
                    }
                    let node = self.graph.add_node(seq[i]);
                    self.node_seq_tracker.push(vec![seq_number]);
                    //println!("insp making edge from {}->{}", self.graph.raw_nodes()[prev.index()].weight, seq[i]);
                    self.graph.add_edge(prev, node, 1); 
                    prev = node;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes
            }
        }
    }
    pub fn consensus(&self) -> Vec<u8> {
        let mut output: Vec<u8> = vec![];
        let mut topo = Topo::new(&self.graph);
        let mut topo_indices = Vec::new();
        let mut max_index = 0;
        let mut max_score = 0.0;

        while let Some(node) = topo.next(&self.graph) {
            topo_indices.push(node);
            if max_index < node.index(){
                max_index = node.index()
            }
        }
        topo_indices.reverse();
        //define score and nextinpath vectors with capacity of num nodes.
        let mut scores: Vec<f64> = vec![0.0; max_index + 1];
        let mut next_in_path: Vec<usize> = vec![0; max_index + 1];
        //iterate thorugh the nodes in revere
        for node in topo_indices{
            //print!("\nstart node: {:?}", self.graph.raw_nodes()[node.index()].weight);
            let mut best_weight_score_edge: (i32, f64, usize) = (-1 , -1.0, 123456789);
            //let mut outEdges = self.graph.neighbors_directed(node, Outgoing).detach();
            let mut neighbour_nodes = self.graph.neighbors_directed(node, Outgoing);
            while let Some(neighbour_node) = neighbour_nodes.next() {
                //print!(" end node: {:?}", self.graph.raw_nodes()[neighbour_node.index()].weight);
                let mut edges = self.graph.edges_connecting(node, neighbour_node);
                let mut weight: i32 = 0;
                while let Some(edge) = edges.next() {
                    weight += edge.weight().clone();
                    //print!(" Edge of weight {}", weight);
                }
                let weight_score_edge = (weight, scores[neighbour_node.index()], neighbour_node.index());
                if weight_score_edge > best_weight_score_edge{
                    best_weight_score_edge = weight_score_edge;
                }
            }
            //save score and traceback
            if best_weight_score_edge.0 as f64 + best_weight_score_edge.1 > max_score{
                max_score = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
            }
            scores[node.index()] = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
            next_in_path[node.index()] = best_weight_score_edge.2;
        }
        //println!("{:?}", scores);
        //println!("{:?}", next_in_path);
        let mut pos = scores.iter().position(|&r| r == max_score).unwrap();
        //using traceback print out the max sequence
        //println!("Consensus");
        while pos != 123456789 {
            output.push(self.graph.raw_nodes()[pos].weight);
            pos = next_in_path[pos];
        }
        output
    }
}

impl Band {
    // Create new Band instance with given size
    //
    // # Arguments
    //
    // * `m` - the expected size of x
    // * `n` - the expected size of y
    //
    fn new(m: usize, n: usize) -> Self {
        Band {
            rows: m + 1,
            cols: n + 1,
            ranges: vec![m + 1..0; n + 1],
        }
    }

    // Add cells around a kmer of length 'k', starting at 'start', which are within a
    // distance of 'w' in x or y directions to the band.
    fn add_kmer(&mut self, start: (u32, u32), k: usize, w: usize) {
        let (r, c) = (start.0 as usize, start.1 as usize);
        // println!("{} {} {}", r, k, self.rows);
        debug_assert!(r + k <= self.rows);
        debug_assert!(c + k <= self.cols);

        if k == 0 {
            return;
        }

        let i = r.saturating_sub(w);
        for j in c.saturating_sub(w)..min(c + w + 1, self.cols) {
            self.ranges[j].start = min(self.ranges[j].start, i);
        }

        let mut i = r.saturating_sub(w);
        for j in min(c + w, self.cols)..min(c + k + w, self.cols) {
            self.ranges[j].start = min(self.ranges[j].start, i);
            i += 1;
        }

        let mut i = r + w + k;
        let mut j = (c + k - 1).saturating_sub(w);
        loop {
            if j <= c.saturating_sub(w) {
                break;
            }
            j -= 1;
            i -= 1;
            self.ranges[j].end = max(self.ranges[j].end, min(i, self.rows));
        }

        let i = min(r + w + k, self.rows);
        for j in (c + k - 1).saturating_sub(w)..min(c + k + w, self.cols) {
            self.ranges[j].end = max(self.ranges[j].end, i);
        }
    }

    // Add cells around a specific position to the band. An cell which is within 'w' distance
    // in x or y directions are added
    fn add_entry(&mut self, pos: (u32, u32), w: usize) {
        let (r, c) = (pos.0 as usize, pos.1 as usize);

        let istart = r.saturating_sub(w);
        let iend = min(r + w + 1, self.rows);
        for j in c.saturating_sub(w)..min(c + w + 1, self.cols) {
            self.ranges[j].start = min(self.ranges[j].start, istart);
            self.ranges[j].end = max(self.ranges[j].end, iend);
        }
    }

    // Each gap generates a line from the start to end.
    fn add_gap(&mut self, start: (u32, u32), end: (u32, u32), w: usize) {
        let nrows = end.0 - start.0;
        let ncols = end.1 - start.1;
        if nrows > ncols {
            for r in start.0..end.0 {
                let c = start.1 + (end.1 - start.1) * (r - start.0) / (end.0 - start.0);
                self.add_entry((r, c), w);
            }
        } else {
            for c in start.1..end.1 {
                let r = start.0 + (end.0 - start.0) * (c - start.1) / (end.1 - start.1);
                self.add_entry((r, c), w);
            }
        }
    }

    // The band needs to start either at (0,0) or at a point that is zero score from (0,0).
    // This naturally sets the start positions correctly for global, semiglobal and local
    // modes. Similarly the band has to either end at (m,n) or at a point from which there is
    // a zero score path to (m,n).
    //
    // At the minimum, irrespective of the score (0,0)->start or end->(m,n), we extend the band
    // diagonally for a length "lazy_extend"(2k) or when it hits the corner, whichever happens first
    //
    // start - the index of the first matching kmer in LCSk++
    // end - the index of the last matching kmer in LCSk++
    //
    fn set_boundaries<F: MatchFunc>(
        &mut self,
        start: (u32, u32),
        end: (u32, u32),
        k: usize,
        w: usize,
        scoring: &Scoring<F>,
    ) {
        let lazy_extend: usize = 2 * k;

        // -------------- START --------------
        // Nothing to do if the start is already at (0,0)
        let (r, c) = (start.0 as usize, start.1 as usize);
        if !(r == 0usize && c == 0usize) {
            let mut score_to_start = if r > 0 { scoring.xclip_prefix } else { 0i32 };
            score_to_start += if c > 0 { scoring.yclip_prefix } else { 0i32 };

            if score_to_start == 0 {
                // Just do a "lazy_extend"
                // First diagonally
                let d = min(lazy_extend, min(r, c));
                self.add_kmer(((r - d) as u32, (c - d) as u32), d, w);

                // If we hit one of the edges before completing lazy_extend
                self.add_gap(
                    (
                        r.saturating_sub(lazy_extend) as u32,
                        c.saturating_sub(lazy_extend) as u32,
                    ),
                    ((r - d) as u32, (c - d) as u32),
                    w,
                );
            } else {
                // we need to find a zero cost cell

                // First try the diagonal
                let diagonal_score = match r.cmp(&c) {
                    // We will hit (r-c, 0)
                    Ordering::Greater => scoring.xclip_prefix,
                    // We will hit (0, c-r)
                    Ordering::Less => scoring.yclip_prefix,
                    Ordering::Equal => 0,
                };

                if diagonal_score == 0 {
                    let d = min(r, c);
                    self.add_kmer(((r - d) as u32, (c - d) as u32), d, w);
                    // Make sure we do at least "lazy_extend" extension
                    let start = (
                        r.saturating_sub(lazy_extend) as u32,
                        c.saturating_sub(lazy_extend) as u32,
                    );
                    let end = ((r - d) as u32, (c - d) as u32);
                    if (start.0 <= end.0) && (start.1 <= end.1) {
                        self.add_gap(start, end, w);
                    }
                } else {
                    // Band to origin
                    self.add_gap((0u32, 0u32), start, w);
                }
            }
        }

        // -------------- END --------------
        // Nothing to do if the last kmer ends at (m, n)
        let (r, c) = (end.0 as usize + k, end.1 as usize + k);
        debug_assert!(r <= self.rows);
        debug_assert!(c <= self.cols);
        if !(r == self.rows && c == self.cols) {
            let mut score_from_end = if r == self.rows {
                0
            } else {
                scoring.xclip_suffix
            };
            score_from_end += if c == self.cols {
                0
            } else {
                scoring.yclip_suffix
            };

            if score_from_end == 0 {
                // Just a lazy_extend
                let d = min(lazy_extend, min(self.rows - r, self.cols - c));
                self.add_kmer((r as u32, c as u32), d, w);

                let r1 = min(self.rows, r + d) - 1;
                let c1 = min(self.cols, c + d) - 1;
                let r2 = min(self.rows, r + lazy_extend);
                let c2 = min(self.cols, c + lazy_extend);
                if (r1 <= r2) && (c1 <= c2) {
                    self.add_gap((r1 as u32, c1 as u32), (r2 as u32, c2 as u32), w);
                }
            } else {
                // we need to find a zero cost cell

                // First try the diagonal
                let dr = self.rows - r;
                let dc = self.cols - c;
                let diagonal_score = match dr.cmp(&dc) {
                    // We will hit (r+dc, self.cols)
                    Ordering::Greater => scoring.xclip_suffix,
                    // We will hit (self.rows, c+dr)
                    Ordering::Less => scoring.yclip_suffix,
                    // We will hit the corner
                    Ordering::Equal => 0,
                };

                if diagonal_score == 0 {
                    let d = min(dr, dc);
                    self.add_kmer((r as u32, c as u32), d, w);
                    // Make sure we do at least "lazy_extend" extension
                    let r1 = min(self.rows, r + d) - 1;
                    let c1 = min(self.cols, c + d) - 1;
                    let r2 = min(self.rows, r + lazy_extend);
                    let c2 = min(self.cols, c + lazy_extend);
                    if (r1 <= r2) && (c1 <= c2) {
                        self.add_gap((r1 as u32, c1 as u32), (r2 as u32, c2 as u32), w);
                    }
                } else {
                    // Band to lower right corner
                    let rows = self.rows as u32;
                    let cols = self.cols as u32;
                    self.add_gap((r as u32, c as u32), (rows as u32, cols as u32), w);
                }
            }
        }
    }

    fn create<F: MatchFunc>(
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        k: usize,
        w: usize,
        scoring: &Scoring<F>,
    ) -> Band {
        let matches = sparse::find_kmer_matches(x, y, k);
        Band::create_with_matches(x, y, k, w, scoring, &matches)
    }

    fn create_with_prehash<F: MatchFunc>(
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        k: usize,
        w: usize,
        scoring: &Scoring<F>,
        y_kmer_hash: &HashMapFx<&[u8], Vec<u32>>,
    ) -> Band {
        let matches = sparse::find_kmer_matches_seq2_hashed(x, y_kmer_hash, k);
        Band::create_with_matches(x, y, k, w, scoring, &matches)
    }

    fn create_with_matches<F: MatchFunc>(
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        k: usize,
        w: usize,
        scoring: &Scoring<F>,
        matches: &[(u32, u32)],
    ) -> Band {
        if matches.is_empty() {
            let mut band = Band::new(x.len(), y.len());
            band.full_matrix();
            return band;
        }

        let match_score = match scoring.match_scores {
            Some((m, _)) => m,
            None => DEFAULT_MATCH_SCORE,
        };

        let res = sparse::sdpkpp(
            matches,
            k,
            match_score as u32,
            scoring.gap_open,
            scoring.gap_extend,
        );
        Band::create_from_match_path(x, y, k, w, scoring, &res.path, matches)
    }

    fn create_from_match_path<F: MatchFunc>(
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        k: usize,
        w: usize,
        scoring: &Scoring<F>,
        path: &[usize],
        matches: &[(u32, u32)],
    ) -> Band {
        let mut band = Band::new(x.len(), y.len());

        if matches.is_empty() {
            band.full_matrix();
            return band;
        }

        let ps = path[0];
        let pe = path[path.len() - 1];

        // Set the boundaries
        band.set_boundaries(matches[ps], matches[pe], k, w, scoring);
        let mut prev: Option<(u32, u32)> = None;

        for &idx in path {
            let curr = matches[idx];
            if curr.continues(prev) {
                let p = prev.unwrap();
                band.add_entry((p.0 + k as u32, p.1 + k as u32), w);
            } else {
                if let Some(p) = prev {
                    band.add_gap((p.0 + (k - 1) as u32, p.1 + (k - 1) as u32), curr, w)
                }
                band.add_kmer(curr, k, w);
            }
            prev = Some(curr);
        }
        band
    }

    fn full_matrix(&mut self) {
        self.ranges.clear();
        self.ranges.resize(self.cols, 0..self.rows);
    }

    fn num_cells(&self) -> usize {
        let mut banded_cells = 0;
        for j in 0..self.ranges.len() {
            banded_cells += self.ranges[j].end.saturating_sub(self.ranges[j].start);
        }
        banded_cells
    }

    #[allow(dead_code)]
    fn visualize(&self) {
        let mut view = vec!['.'; self.rows * self.cols];
        let index = |i, j| i * self.cols + j;
        for j in 0..self.ranges.len() {
            let range = &self.ranges[j];
            for i in range.start..range.end {
                view[index(i, j)] = 'x';
            }
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                print!("{}", view[index(i, j)]);
            }
            println!();
        }
    }

    #[allow(dead_code)]
    fn stat(&self) {
        let total_cells = self.rows * self.cols;
        let banded_cells = self.num_cells();
        let percent_cells = (banded_cells as f64) / (total_cells as f64) * 100.0;
        println!(
            " {} of {} cells are in the band ({2:.2}%)",
            banded_cells, total_cells, percent_cells
        );
    }
}
