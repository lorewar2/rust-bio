
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
pub struct Aligner<F: MatchFunc> {
    traceback: Traceback,
    query: Vec<u8>,
    pub poa: Poa<F>,
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
    pub fn add_to_graph(&mut self) -> &mut Self {
        let alignment = self.traceback.alignment();
        self.poa.add_alignment(&alignment, &self.query);
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
    pub fn add_alignment(&mut self, aln: &Alignment, seq: TextSlice) {
        let head = Topo::new(&self.graph).next(&self.graph).unwrap();
        let mut prev: NodeIndex<usize> = NodeIndex::new(head.index());
        let mut i: usize = 0;
        let mut edge_not_connected: bool = false;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => {
                    let node: NodeIndex<usize> = NodeIndex::new(0);
                    if (seq[i] != self.graph.raw_nodes()[head.index()].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        prev = node;
                    }
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                        edge_not_connected = false;
                    }
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(*p);
                    if (seq[i] != self.graph.raw_nodes()[*p].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                *self.graph.edge_weight_mut(edge).unwrap() += 1;
                            }
                            None => {
                                if prev.index() != head.index() {
                                    self.graph.add_edge(prev, node, 1);
                                }
                            }
                        }
                        prev = NodeIndex::new(*p);
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    let node = self.graph.add_node(seq[i]);
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);    
                    }
                    prev = node;
                    edge_not_connected = true;
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes
            }
        }
    }
    pub fn consensus(&self) -> (Vec<u8>, Vec<usize>) {
        let mut output: Vec<u8> = vec![];
        let mut topopos: Vec<usize> = vec![];
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
        let mut weight_scores: Vec<i32> = vec![0; max_index + 1];
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
            weight_scores[node.index()] = best_weight_score_edge.0;
        }
        let mut pos = scores.iter().position(|&r| r == max_score).unwrap();
        //calculate the start weight score
        let mut consensus_started: bool = false;
        let weight_average = scores[pos] / scores.len() as f64;
        let weight_threshold = weight_average as i32 / 2; 
        while pos != 123456789 {
            //continue if starting weight score is too low
            if consensus_started == false && weight_scores[pos] < weight_threshold {
                pos = next_in_path[pos];
                continue;
            }
            consensus_started = true;
            topopos.push(pos as usize);
            output.push(self.graph.raw_nodes()[pos].weight);
            pos = next_in_path[pos];
        }
        (output, topopos)
    }
}
