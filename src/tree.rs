use seq::*;
use align::*;


#[derive(Debug,RustcDecodable,RustcEncodable)]
pub enum SeqNode {
    Frag { id: u32, val: Sequence },
    Splat,
    Dist { id: u32, scores: ProbMatr },
    List { id: u32, members: Vec<SeqNode> },
    Branch { id: u32, members: Vec<SeqNode> }
}
