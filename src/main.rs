extern crate rustc_serialize;
use rustc_serialize::json;
use std::ops::{Index,IndexMut};
use std::clone::Clone;

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
}



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



pub enum Reference {
    ProbMatr( Matrix<f32> ),  // used to store HMMer-style probabilities
    Simple( Sequence )        // simple reference sequence
}


#[derive(Debug,RustcDecodable,RustcEncodable)]
enum SeqNode {
    Frag { id: u32, val: Sequence },
    Splat,
    Dist { id: u32, scores: [f32; MMERS] },
    List { id: u32, members: Vec<SeqNode> },
    Branch { id: u32, members: Vec<SeqNode> }
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
}
