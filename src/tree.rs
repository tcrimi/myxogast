extern crate serde_json;

use std::collections::BTreeMap;
use self::serde_json::Value as JSON_Val;
use seq::*;
use align::*;
use matrix::*;


#[derive(Debug, Clone)]
pub enum SeqNode {
    Frag { id: u32, val: Sequence },
    Splat,
    //Dist { id: u32, scores: ProbMatr },
    List { id: u32, members: Vec<SeqNode> },
    Branch { id: u32, members: Vec<SeqNode> }
}

#[derive(Debug)]
pub enum SeqErr {
    Ambiguous,
    BadJsonElement,
    Unsupported,
    StringExpected,
    Unknown
}

#[derive(Debug)]
pub struct SeqGraph {
    root: SeqNode,
    names: BTreeMap<u32, String>
}

impl SeqNode {
    fn read_item( idx: &mut u32, names : &mut BTreeMap<u32, String>, elem: &JSON_Val )
                   -> Result<SeqNode, SeqErr> {
        let ret = match elem {
            &JSON_Val::String(ref s) => Ok(SeqNode::Frag{id: *idx, val: Sequence::from_str(&s).unwrap() }),
            &JSON_Val::Object(ref obj) => SeqNode::read_btree( idx, names, &obj ),
            _ => Err(SeqErr::BadJsonElement)
        };
        *idx += 1;
        ret
    }

    fn read_list( idx: &mut u32, names : &mut BTreeMap<u32, String>, elem: &JSON_Val )
                   -> Result<Vec<SeqNode>, SeqErr> {
        match elem {
            &JSON_Val::Array(ref l) => {
                let mut ret = Vec::new();
                for (i, x) in l.iter().enumerate() {
                    ret.push( SeqNode::read_item( idx, names, x ).unwrap() );
                    *idx += 1;
                };
                Ok(ret)
            },
            _ => Err(SeqErr::BadJsonElement)
        }
    }

    fn read_btree( idx: &mut u32, names : &mut BTreeMap<u32, String>, map : &BTreeMap<String, JSON_Val>)
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
                        let ret = SeqNode::Frag { id: *idx, val: Sequence::from_str(&s2).unwrap()};
                        *idx += 1;
                        Ok(ret) },
                    _ => Err(SeqErr::StringExpected)
                }
            }
        } else if map.contains_key("branch") {
            if map.contains_key("dist") {
                Err(SeqErr::Ambiguous)
            } else {
                let v = try!{ SeqNode::read_list( idx, names, &map.get("branch").unwrap()) };
                let ret = SeqNode::Branch { id: *idx+1, members: v };
                *idx += 2;
                Ok(ret) }
        } else if map.contains_key("dist") {
            Err(SeqErr::Unsupported)
        } else {
            Err(SeqErr::Unknown)
        }
    }

    pub fn iden(&self) -> Option<u32> {
        match self {
            &SeqNode::Frag{ id: ref id, ..} => Some(id.clone()),
            &SeqNode::Branch { id: ref id, ..} => Some(id.clone()),
            &SeqNode::List { id: ref id, ..} => Some(id.clone()),
            &SeqNode::Splat => None
        }
    }
}

impl SeqGraph {
    pub fn from_json( serialized : &str ) -> Result<SeqGraph, SeqErr> {
        let mut names : BTreeMap<u32, String> = BTreeMap::new();
        let mut idx = 0;
        
        let value = serde_json::from_str(serialized).unwrap();
        let mut idx = 0u32;
        let tree = try!{
            match value {
                JSON_Val::Array(_) => {
                    let l = try!{ SeqNode::read_list(&mut idx, &mut names, &value) };
                    Ok( SeqNode::List{ id: idx, members: l })
                },
                _ => SeqNode::read_item(&mut idx, &mut names, &value)
            }
        };
        Ok( SeqGraph { root: tree, names: names } )
    }

    

    fn _longest_path(&self, node: &SeqNode, path: &mut Vec<u32>, path_len: usize, base_ct: usize )
                     -> (/*len_bases*/ usize, /*len_nodes*/ usize, /*path*/ Vec<u32>) {

        match node {
            &SeqNode::Frag{ id: ref id, val: ref val } => {
                if path_len < path.len() {
                    path[ path_len ] = *id;
                } else {
                    path.push( *id );
                };

                (base_ct + val.len(), path_len, path.clone()) },

            &SeqNode::Splat => (base_ct, path_len, path.clone()),

            &SeqNode::Branch { id: ref id, members: ref members } => {
                members.iter()
                    .map( |n| self._longest_path( n, path, path_len, base_ct ) )
                    .max_by_key( |r| r.0 )
                    .unwrap()
            },

            &SeqNode::List { id: ref id, members: ref members } => {
                if members.is_empty() {
                    (base_ct, path_len, path.clone())
                } else {
                    // this ain't pretty:
                    let mut new_members = members.clone();
                    let head = new_members.remove(0);
                    let (head_bases, head_nodes, mut head_path) = self._longest_path( &head, path, path_len, base_ct );
                    // FIXME: how do we assign an id here?
                    self._longest_path( &SeqNode::List{id: u32::max_value(), members: new_members}, &mut head_path,
                                         head_nodes, head_bases )
                }
            },

            //&SeqNode::Dist { id: id, scores: scores } => panic!("dist unsupported")

        }
    }

    pub fn longest_path(&self) -> (/*len_bases*/ usize, /*path*/ Vec<u32>) {
        let (len_bases, len_nodes, path) = self._longest_path( &self.root, &mut Vec::<u32>::new(), 0, 0 );
        (len_bases, path[..len_nodes].to_vec())
    }

    fn _align_by_path( node: &SeqNode, query: &Sequence, m: &mut Matrix<Cell>, base_params: &AlnParams,
                       start: i32, path: &Vec<u32>, pos: usize ) -> (/*offset*/ i32, /*score*/ i32) {

        println!("path: {:?}, pos: {}, expected id: {}, id: {:?}", path, pos, path[pos], node.iden());
        println!("matrix:\n{:?}", m);
     
        match node.iden() {
            Some(id) => assert_eq!( id, path[pos] ),
            None => ()
        }

        assert!( pos < path.len(), "passed end of path");
 
        match node {
            &SeqNode::Frag{ id: ref id, val: ref val } => {
                // ALIGNMENT HAPPENS HERE
                assert_eq!( *id, path[pos] );
                let t = align_matrix( val, query, &base_params, Some(start), m ).unwrap();
                (val.len() as i32, Cell::unpack(&m[t]).unwrap().1 )
            },

            &SeqNode::Splat => panic!("splat outside list"),

            &SeqNode::Branch { id: ref id, members: ref members } => {
                assert_eq!( *id, path[pos] );
                if pos+1 >= path.len() && members.len() > 0 {
                    panic!("bad path: truncated");
                }
                for n in members {
                    match n.iden() {
                        Some(id) => {
                            if id == path[pos+1] {
                                return SeqGraph::_align_by_path( n, query, m, &base_params, start, path, pos+1 )
                            }},
                        None => panic!("splat not allowed in branch")
                    }
                }
                panic!(format!("bad path: didn't find id {} in path", path[pos+1]));
            },

            &SeqNode::List { id: ref id, members: ref members } => {
                // workhorse loop -- resposible for handling splats, positional math, etc
                let mut offset = 0i32;
                let mut params : AlnParams = base_params.clone();
                let mut last_score = -1i32;
                for (i, n) in members.iter().enumerate() {
                    match n {
                        &SeqNode::Frag{..} => {
                            let (f_len, last_score) = SeqGraph::_align_by_path( n, query, m, &params, start,
                                                                                path, pos + 1 + i );
                            offset += f_len;
                            params = AlnParams::copy_but_llocal( &params, false );
                        },
                        &SeqNode::Splat => {
                            params = AlnParams::copy_but_llocal( &params, true );
                        },
                        &SeqNode::Branch{..} => {
                            let (f_len, last_score) = SeqGraph::_align_by_path( n, query, m, &params, start,
                                                                                path, pos + 1 + i );
                            offset += f_len;
                            params = AlnParams::copy_but_llocal( &params, false );
                        },
                        &SeqNode::List {..} => {
                            panic!("lists not permitted within lists");
                        }
                    }
                }
                (offset, last_score)
            }
        }
    }

    pub fn align(&self, query: &Sequence, path: &Vec<u32>, base_params: &AlnParams)
                 -> Option<(/*padded_ref*/ Sequence, /*padded_query*/ Sequence)> {
        let (ref_len, _) = self.longest_path();
        let mut m = Matrix::<Cell>::new( Cell(0), ref_len + 2, query.len() + 2 );

        println!("ref_len: {}", ref_len);
        let (len, score) = SeqGraph::_align_by_path( &self.root, query, &mut m, &base_params, 0, path, 0 );
        println!("len: {}, score: {}", len, score);
        None
    }

    /// find the "local maxima" path -- ie, traverse the graph, trying each option at each branch point
    ///  to find the best path, without worrying about combinitorics
    fn _local_max( node: &SeqNode, query: &Sequence, m: &mut Matrix<Cell>, base_params: &AlnParams,
                   start: i32, path: &mut Vec<u32> )
                   -> (/*offset*/ i32, /*score*/ i32) {
        println!("node: {:?}, query: {}, start: {}, path: {:?}\n", node, query, start, path);

        match node {
            &SeqNode::Frag{ id: ref id, val: ref val } => {
                path.push( *id );
                let t = align_matrix( val, query, &base_params, Some(start), m ).unwrap();
                (start + val.len() as i32, Cell::unpack(&m[t]).unwrap().1)
            }

            &SeqNode::Branch { id: ref id, members: ref members } => {
                path.push( *id );

                // our vec @path gets a little tricky here, as subsequent recursive calls
                //   will extend it, and we only want to choose one of the possible extensions.
                //   the solution is to keep a copy of the best path, and truncate the working
                //   vec to its original state each time
                //
                let len_before = path.len();
                let mut best_path : Vec<u32> = Vec::new();
                let mut best_x = -1;
                let mut best_id = 0;
                let mut best_score = i32::min_value();

                for n in members {
                    match n.iden() {
                        Some(nid) => {
                            println!("branch trying id={}", nid);

                            let (next_x, next_score) = SeqGraph::_local_max( n, query, m, base_params,
                                                                             start, path );
                            if next_score > best_score {
                                println!("id {} score={} > old score {}", nid, next_score, best_score);
                                best_path = path.to_owned();

                                best_id = nid;
                                best_score = next_score;
                            }

                            path.truncate(len_before);
                        },
                        None => {},
                    }
                }
                path.extend( best_path[len_before..].iter().cloned() );
                println!("after branch, path={:?}", path);
                (best_x, best_score)
            },

            &SeqNode::List { id: ref id, members: ref members } => {
                path.push( *id );

                let mut acc = 0i32;
                let mut score = -1;
                for n in members {
                    match n.iden() {
                        Some(nid) => {
                            println!("list iter on id {}", nid);

                            let (len, score) = SeqGraph::_local_max( n, query, m, base_params, start+acc, path );
                            acc += len;
                        },
                        None => {}
                    };
                }
                (acc, score)
            },

            &SeqNode::Splat => SeqGraph::_local_max( node, query, m, base_params, start, path )
        }
    }

    pub fn local_max(&self, query: &Sequence, base_params: &AlnParams ) -> Option<(Sequence, Sequence)> {
        let (ref_len, _) = self.longest_path();
        let mut m = Matrix::<Cell>::new( Cell(0), ref_len + 2, query.len() + 2 );
        let mut v = Vec::new();
        SeqGraph::_local_max( &self.root, query, &mut m, base_params, 0, &mut v );
        println!("-- finished - path={:?}", v);
        None
    }
}
