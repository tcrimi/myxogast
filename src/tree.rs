extern crate serde_json;

//use std::fmt::Write;
use std::io::Write;
use std::collections::BTreeMap;
use self::serde_json::Value as JSON_Val;
use seq::*;
use align::*;
use matrix::*;
use std::iter::Iterator;
use std::rc::Rc;


#[derive(Debug, Clone)]
pub enum SeqNode {
    Nil,

    Frag {
        id: u32,
        val: Sequence,
        llocal: bool,
        rlocal: bool,

        next: Rc<SeqNode>
    },

    //Dist { id: u32, scores: ProbMatr },

    Branch {
        id: u32,
        llocal: bool,
        rlocal: bool,
        members: Vec<SeqNode>
    }
}

#[derive(Debug)]
pub struct SeqGraph {
    root: SeqNode,
    names: BTreeMap<u32, String>
}

/// GraphPath - used for Iterator trait on SeqGraph
#[derive(Debug)]
pub struct GraphPath<'a> {
    curr_node: &'a SeqNode,
    path: Vec<u32>,
    pos: usize
}

#[derive(Debug)]
pub enum SeqErr {
    BadJsonElement,
    Ambiguous,
    StringExpected,
    Unsupported,
    Mismatch
}


impl SeqNode {

    fn dispatch( idx: &mut u32, names : &mut BTreeMap<u32, String>, elem: &JSON_Val, next: Rc<SeqNode> )
                 -> Result<SeqNode, SeqErr> {
        *idx += 1;

        match elem {
            &JSON_Val::String(ref s) => SeqNode::read_str( idx, names, s, next ),
            &JSON_Val::Object(ref map) => SeqNode::read_obj( idx, names, map, next ),
            &JSON_Val::Array(ref l) => SeqNode::read_list( idx, names, l, 0 ),
            _ => Err(SeqErr::BadJsonElement)
        }
    }

    fn read_str( idx: &mut u32, names : &mut BTreeMap<u32, String>, s: &String, next: Rc<SeqNode> )
                 -> Result<SeqNode, SeqErr> {
        Ok(SeqNode::Frag{ id: *idx, val: Sequence::from_str(&s).unwrap(),
                          llocal: false,
                          rlocal: false,
                          next: next.clone() })
    }

    fn read_list( idx: &mut u32, names : &mut BTreeMap<u32, String>, l: &Vec<JSON_Val>, pos: usize )
                  -> Result<SeqNode, SeqErr> {
        *idx += 1;
        if pos < l.len() {
            let next = SeqNode::read_list( idx, names, l, pos+1 ).unwrap();
            SeqNode::dispatch( idx, names, &l[pos], Rc::new(next) )
        } else {
            Ok( SeqNode::Nil )
        }
    }

    fn read_obj( idx: &mut u32, names: &mut BTreeMap<u32, String>, map: &BTreeMap<String, JSON_Val>, next: Rc<SeqNode> )
                 -> Result<SeqNode, SeqErr> {

        let _ = match map.get("id") {
            Some(s) => {
                match s {
                    &JSON_Val::String(ref s2) => { names.insert(*idx, s2.clone()); () },
                    _ => return Err(SeqErr::StringExpected)
                }},
            None => ()
        };

        if map.contains_key("seq") {
            if map.contains_key("dist") || map.contains_key("branch") {
                Err(SeqErr::Ambiguous)
            } else {
                match map.get("seq").unwrap() {
                    &JSON_Val::String(ref s2) => {
                        Ok( SeqNode::Frag { id: *idx, val: Sequence::from_str(&s2).unwrap(),
                                            llocal: false, rlocal: false, next: next.clone() } )
                    },
                    _ => Err(SeqErr::BadJsonElement)
                }
            }
        } else if map.contains_key("branch") {
            if map.contains_key("dist") {
                Err(SeqErr::Ambiguous)
            } else {
                let members = match map.get("branch").unwrap() {
                    &JSON_Val::Array(ref _l) => {
                        let l : &Vec<JSON_Val> = _l;
                        let mut m = Vec::new();
                        for x in l {
                            *idx += 1;
                            m.push( SeqNode::dispatch( idx, names, x, next.clone() ).unwrap());
                        }
                        Ok(m)
                        //Ok( l.iter().map(|x| SeqNode::dispatch( idx, names, x, next ).unwrap() ).collect() )
                    },
                    _ => Err(SeqErr::BadJsonElement)
                }.unwrap();
                *idx += 1;
                Ok( SeqNode::Branch { id: *idx, members: members, llocal: false, rlocal: false } )
            }
        } else if map.contains_key("dist") {
            Err(SeqErr::Unsupported)
        } else {
            Err(SeqErr::BadJsonElement)
        }
    }

    pub fn iden(&self) -> Option<u32> {
        match self {
            &SeqNode::Nil => None,
            &SeqNode::Frag { id: ref id, ..} => Some(id.clone()),
            &SeqNode::Branch { id: ref id, ..} => Some(id.clone())
        }
    }
}

impl SeqGraph {
    pub fn from_json( serialized : &str ) -> Result<SeqGraph, SeqErr> {
        let mut names : BTreeMap<u32, String> = BTreeMap::new();
        let mut idx = 0;
        
        let value = serde_json::from_str(serialized).unwrap();
        let mut idx = 0u32;
        let tree = SeqNode::dispatch( &mut idx, &mut names, &value, Rc::new(SeqNode::Nil) ).unwrap();
        Ok( SeqGraph { root: tree, names: names } )
    }

    fn _max_len(n: &SeqNode) -> usize {
        match n {
            &SeqNode::Nil => 0,
            &SeqNode::Frag { val: ref val, next: ref next, ..} => { val.len() + SeqGraph::_max_len(next) },
            &SeqNode::Branch { members: ref members, ..} => { members.iter()
                                                              .map( |n| SeqGraph::_max_len(n) )
                                                              .max()
                                                              .unwrap() }
        }
    }
    /// max_len - the maximum sequence length encoded by this graph
    pub fn max_len(&self) -> usize {
        SeqGraph::_max_len( &self.root )
    }


    fn _align__global_max(node: &SeqNode, query: &Sequence, m: &mut Matrix<Cell>,
                          base_params: &AlnParams, start: i32, path: &mut Vec<u32>, pos: usize )
                          -> Option<(/*score*/ i32, /*path*/ Vec<u32>)> {
        match node {
            &SeqNode::Nil => Some((0, path.clone())),
            &SeqNode::Frag { id: ref id, val: ref val, next: ref next, ..} => {
                path.push( id.clone() );

                let t = align_matrix( val, query, &base_params, Some(start), m ).unwrap();
                let score = Cell::unpack(&m[t]).unwrap().1;
                let (next_score, _) = SeqGraph::_align__global_max( next, query, m, base_params,
                                                                    start + val.len() as i32, path,
                                                                    pos + 1 )
                    .unwrap();
                Some((score + next_score, path.to_vec()))
            },
            &SeqNode::Branch { members: ref members, ..} => {
                path.push( node.iden().unwrap() );

                let mut best_node = &SeqNode::Nil;
                let mut best_score = i32::min_value();
                let mut best_path = Vec::new();

                for n in members {
                    path.truncate( pos + 1 );
                    let (next_score, _) = SeqGraph::_align__global_max( n, query, m, base_params,
                                                                        start, path, pos + 1 )
                        .unwrap();
                    if next_score > best_score {
                        best_score = next_score;
                        best_node = &n;
                        best_path = path.to_vec();
                    }
                }
                best_path.push( best_node.iden().unwrap() );
                Some((best_score, best_path.to_vec()))
            }
        }
    }
    /// SeqGraph::align__global_max -- align query to graph, testing every possible branch to
    ///   find the global maximum.
    pub fn align__global_max(&self, query: &Sequence, base_params: &AlnParams )
                            -> Option<(/*path*/ Vec<u32>, /*padded_ref*/ Sequence, /*padded_query*/ Sequence)> {
        let mut _path = Vec::new();
        let ref_len = self.max_len();
        let mut m = Matrix::<Cell>::new( Cell(0), ref_len + 2, query.len() + 2 );
        match SeqGraph::_align__global_max( &self.root, query,  &mut m, base_params, 0,
                                             &mut _path, 0) {
            Some((score, path)) => {
                println!("align__global_max -> {:?}", (score,&path) );
                let mut full_ref_v = Vec::new();
                for s in GraphPath::from_graph( self, path.to_vec() ) {
                    full_ref_v.extend( s.0 );
                }
                let full_ref = Sequence(full_ref_v);
                let (_, t) = m.max();

                let (fwd_r, fwd_q) = aln_from_coord( &t.0, &t.1, &1, &full_ref, &query, &m );
                let (rev_r1, rev_q1) = aln_from_coord( &t.0, &t.1, &(-1), &full_ref, &query, &m );
                let mut rev_r2 = rev_r1.reverse().0;
                if t.0 <= full_ref.len() as i32 {
                    rev_r2.push( full_ref[t.0-1] );
                }

                let mut rev_q2 = rev_q1.reverse().0;
                if t.1 <= query.len() as i32 {
                    rev_q2.push( query[t.1-1] );
                }
                Some((path, Sequence(rev_r2) + fwd_r, Sequence(rev_q2) + fwd_q))
            },
            None => None
        }
    }
}

impl<'a> GraphPath<'a> {
    pub fn from_graph(graph: &'a SeqGraph, path: Vec<u32>) -> GraphPath {
        GraphPath {
            curr_node: &graph.root,
            path: path,
            pos: 0
        }
    }

    fn _next(&mut self) -> Option<Sequence> {
        if self.pos < self.path.len() {
            assert_eq!( self.path[self.pos], self.curr_node.iden().unwrap() );

            match self.curr_node {
                &SeqNode::Frag { val: ref val, next: ref next, ..} => {
                    println!("found frag: {:?}", val);
                    self.pos += 1;
                    self.curr_node = next;
                    Some(val.clone())
                },
                &SeqNode::Branch { members: ref members, ..} => {
                    assert_eq!( self.curr_node.iden(), Some(self.path[self.pos]) );
                    self.pos += 1;
                    for n in members {
                        if n.iden() == Some(self.path[self.pos]) {
                            self.curr_node = n;
                            return self._next()
                        }
                    }
                    panic!("couldn't find thing");
                },
                &SeqNode::Nil => None
            }
        } else {
            None
        }
    }
}

impl<'a> Iterator for GraphPath<'a> {
    type Item = Sequence;

    fn next(&mut self) -> Option<Sequence> {
        if self.pos < self.path.len() {
            self._next()
        } else {
            None
        }
    }
}
