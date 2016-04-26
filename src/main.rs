extern crate rustc_serialize;
use rustc_serialize::json;
use std::ops::{Index,IndexMut};
use std::clone::Clone;
//use core::cmp::Eq;

const MMERS : usize = 4;

pub type Mmer = u8;

#[derive(Debug,RustcDecodable,RustcEncodable)]
pub struct Sequence( Vec<Mmer> );

impl Sequence {
    pub fn from_str( seq : &str ) -> Result<Sequence, String> {
        let x = seq.to_uppercase()
            .chars()
            .map( |ch| match ch {
                'A' => Ok(0),
                'T' => Ok(1),
                'G' => Ok(2),
                'C' => Ok(3),
                _   => Err(format!("unrecognized base: {}", ch))
            } )
            .collect();
        match x {
            Ok(arr) => Ok( Sequence(arr) ),
            Err(msg) => Err(msg)
        }
    }
    fn len(&self) -> usize { self.0.len() }
}


#[derive(Debug,RustcDecodable,RustcEncodable)]
pub struct Matrix<T:Clone> {
    width:   usize,
    height:  usize,
    data :   Vec<T>
}

impl<T:Clone> Matrix<T> {
    pub fn new( init: T, w: usize, h: usize ) -> Matrix<T> {
        Matrix {
            width:  w,
            height: h,
            data : {
                // FIXME: why doesn't the vec! macro work here?  this seems inefficient
                let mut v : Vec<T> = Vec::new();
                for _ in 0..(w * h) {
                    v.push(init.clone());
                }
                v
            }
        }
    }
}

impl<T:Clone> Index<(i32, i32)> for Matrix<T> {
    type Output = T;
    fn index<'a>(&'a self, _index: (i32, i32)) -> &'a T {
        let (a, b) = _index;
        let a2 : usize = if a < 0 {self.width - (a.abs() as usize)} else {a as usize};
        let b2 : usize = if b < 0 {self.height - (b.abs() as usize)} else {b as usize};
        return &self.data[ a2 * self.width + b2 ];
    }
}

impl<T:Clone> IndexMut<(i32, i32)> for Matrix<T> {
    fn index_mut<'a>(&'a mut self, _index: (i32, i32)) -> &'a mut T {
        let (a, b) = _index;
        let a2 : usize = if a < 0 {self.width - (a.abs() as usize)} else {a as usize};
        let b2 : usize = if b < 0 {self.height - (b.abs() as usize)} else {b as usize};
        return &mut self.data[ a2 * self.width + b2 ];
    }
}


// used to store HMMer-style probabilities; dimensions: BASE-COUNT x REF-LENGTH
pub type ProbMatr = Matrix<f32>;


#[derive(Debug,RustcDecodable,RustcEncodable)]
enum SeqNode {
    Frag { id: u32, val: Sequence },
    Splat,
    Dist { id: u32, scores: ProbMatr },
    List { id: u32, members: Vec<SeqNode> },
    Branch { id: u32, members: Vec<SeqNode> }
}


pub struct AlnParams {
    llocal:     bool,           // global or local on the left?
    rlocal:     bool,           // global or local on the right?
    max_indel:  Option<u8>,  // maximum number of indels before short-circuiting
    gap_open:   i16,          // penalty for opening a gap
    gap_ext:    i16,          // gap extention penalty
}

impl AlnParams {
    pub fn new(  _llocal: Option<bool>, _rlocal: Option<bool>, _max_indel: Option<u8>,
                 _gap_open: Option<i16>, _gap_ext: Option<i16> ) -> AlnParams {
        AlnParams {
            llocal:     match _llocal { Some(b) => b, None => true },
            rlocal:     match _rlocal { Some(b) => b, None => true },
            max_indel:  _max_indel,
            gap_open:   match _gap_open { Some(g) => g, None => -1 },
            gap_ext:    match _gap_ext { Some(g) => g, None => -1 },

        }
    }
}

#[derive(Debug)]
pub enum AlnState {
    Match,
    Mismatch,
    Ins,
    Del
}

// NOTE: I'm a little shocked I can't just #derive this...
impl PartialEq for AlnState {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (&AlnState::Match, &AlnState::Match) => true,
            (&AlnState::Mismatch, &AlnState::Mismatch) => true,
            (&AlnState::Ins, &AlnState::Ins) => true,
            (&AlnState::Del, &AlnState::Del) => true,
            _ => false }
    }
}
impl Eq for AlnState {} 

// bit-packing operations
// FIXME: i16 might not be enough, so I should investigate a different division, eg
//        4 bits for AlnState and 28 bits get munged into an i32
pub fn pack_cell( state : AlnState, score : i16 ) -> i32 {
    let _state : i32 = match state {
        AlnState::Match => (0 << 16) | 0xffff,
        AlnState::Mismatch => (1 << 16) | 0xffff,
        AlnState::Ins => (2 << 16) | 0xffff,
        AlnState::Del => (3 << 16) | 0xffff };
    let _score : i32 = (score as i32) | 0xffff0000;
    (_score & _state)
}

pub fn unpack_cell( packed : i32 ) -> Result<(AlnState, i16), ()> {
    let state = try!(
        match packed >> 16 {
            0 => Ok(AlnState::Match),
            1 => Ok(AlnState::Mismatch),
            2 => Ok(AlnState::Ins),
            3 => Ok(AlnState::Del),
            _ => Err(()) }
    );
    Ok( (state, packed as i16) )
}

pub fn align( reference: Sequence, query: Sequence, params: AlnParams ) -> Option<()> {
    let mut m = Matrix::<i32>::new( 0, reference.len() + 1, query.len() + 1 );
    
    None
}

pub fn align_hmm( reference: ProbMatr, query: Sequence, params: AlnParams ) -> Option<()> {
    let mut m = Matrix::<f32>::new( 0., reference.width + 1, query.len() + 1 );
    None
}


fn main() {
    let x = SeqNode::Frag { id: 0, val: Sequence( vec![0,1,2,3] ) };
    println!("Hello World: {:?}", x );
    println!("{:?}", json::encode(&x));
    let mut y = Matrix::<u32>::new(0, 10, 10);
    y[(1,0)] = 5u32;
    println!("5 == {}?", y[(1,0)]);

    let q = Sequence::from_str("atgc");
    println!("{:?}", q);
    let r = Sequence::from_str("atgcn");
    println!("{:?}", r);

    let _ = pack_cell( AlnState::Match, -1 );
    let _ = pack_cell( AlnState::Mismatch, -1 );
    let _ = pack_cell( AlnState::Del, -1 );
    assert_eq!( unpack_cell( pack_cell( AlnState::Ins, -10 ) ).unwrap(), (AlnState::Ins, -10) );
    assert_eq!( unpack_cell( pack_cell( AlnState::Match, 222 ) ).unwrap(), (AlnState::Match, 222) );

}
