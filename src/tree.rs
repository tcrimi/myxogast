extern crate serde_json;

use std::collections::BTreeMap;
use self::serde_json::Value as JSON_Val;
use seq::*;
use align::*;



#[derive(Debug)]
pub enum SeqNode {
    Frag { id: u32, val: Sequence },
    Splat,
    Dist { id: u32, scores: ProbMatr },
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
    tree: SeqNode,
    names: BTreeMap<u32, String>
}

impl SeqNode {
    fn read_item( idx: u32, names : &mut BTreeMap<u32, String>, elem: &JSON_Val ) -> Result<SeqNode, SeqErr> {
        match elem {
            &JSON_Val::String(ref s) => Ok(SeqNode::Frag{id: 0, val: Sequence::from_str(&s).unwrap() }),
            &JSON_Val::Object(ref obj) => SeqNode::read_btree( idx, names, &obj ),
            _ => Err(SeqErr::BadJsonElement)
        }
    }

    fn read_list( idx: u32, names : &mut BTreeMap<u32, String>, elem: &JSON_Val ) -> Result<Vec<SeqNode>, SeqErr> {
        match elem {
            &JSON_Val::Array(ref l) => 
                Ok( l.iter().map( |x| SeqNode::read_item(idx, names, x).unwrap() ).collect() ),
            _ => Err(SeqErr::BadJsonElement)
        }
    }

    fn read_btree( idx : u32, names : &mut BTreeMap<u32, String>, map : &BTreeMap<String, JSON_Val>)
                   -> Result<SeqNode, SeqErr> {

        let _ = match map.get("id") {
            Some(s) => {
                match s {
                    &JSON_Val::String(ref s2) => { names.insert(idx, s2.clone()); () },
                    _ => return Err(SeqErr::StringExpected)
                }},
            None => ()
        };

        if map.contains_key("seq") {
            if map.contains_key("dist") || map.contains_key("branch") {
                Err(SeqErr::Ambiguous)
            } else {
                match map.get("seq").unwrap() {
                    &JSON_Val::String(ref s2) => Ok( SeqNode::Frag { id: idx, val: Sequence::from_str(&s2).unwrap()}),
                    _ => Err(SeqErr::StringExpected)
                }
            }
        } else if map.contains_key("branch") {
            if map.contains_key("dist") {
                Err(SeqErr::Ambiguous)
            } else {
                let v = try!{ SeqNode::read_list(idx, names, &map.get("branch").unwrap()) };
                Ok( SeqNode::Branch { id: idx, members: v })
            }
        } else if map.contains_key("dist") {
            Err(SeqErr::Unsupported)
        } else {
            Err(SeqErr::Unknown)
        }
    }

}

impl SeqGraph {
    pub fn from_json( serialized : &str ) -> Result<SeqGraph, SeqErr> {
        let mut names : BTreeMap<u32, String> = BTreeMap::new();
        let mut idx = 0;
        
        let value = serde_json::from_str(serialized).unwrap();
        let tree = try!{
            match value {
                JSON_Val::Array(_) => {
                    let l = try!{ SeqNode::read_list(0, &mut names, &value) };
                    Ok( SeqNode::List{ id: idx, members: l })
                },
                _ => SeqNode::read_item(0, &mut names, &value)
            }
        };
        Ok( SeqGraph { tree: tree, names: names } )
    }
}
