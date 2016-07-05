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

    fn dispatch( idx: &mut u32, names : &mut BTreeMap<u32, String>, elem: &JSON_Val, next: SeqNode )
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
                          next: Rc::new(next) })
    }

    fn read_list( idx: &mut u32, names : &mut BTreeMap<u32, String>, l: &Vec<JSON_Val>, pos: usize )
                  -> Result<SeqNode, SeqErr> {
        *idx += 1;
        if pos < l.len() {
            let next = SeqNode::read_list( idx, names, l, pos+1 ).unwrap();
            SeqNode::dispatch( idx, names, &l[pos], next )
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
                            m.push( SeqNode::dispatch( idx, names, x, next ).unwrap());
                        }
                        Ok(m)
                        //Ok( l.iter().map(|x| SeqNode::dispatch( idx, names, x, next ).unwrap() ).collect() )
                    },
                    _ => Err(SeqErr::BadJsonElement)
                }.unwrap();
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
        let tree = SeqNode::dispatch( &mut idx, &mut names, &value, SeqNode::Nil ).unwrap();
        Ok( SeqGraph { root: tree, names: names } )
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
            self.pos += 1;

            match self.curr_node {
                &SeqNode::Frag { val: ref val, next: ref next, ..} => {
                    self.curr_node = next;
                    Some(val.clone())
                },
                &SeqNode::Branch { members: ref members, ..} => {
                    for n in members {
                        if n.iden().unwrap() == self.path[self.pos] {
                            self.curr_node = n;
                            return self._next()
                        }
                    }
                    None
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
