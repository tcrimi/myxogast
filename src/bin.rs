extern crate rustc_serialize;
extern crate myxogast;

use myxogast::align::*;
use myxogast::tree::*;
use myxogast::seq::*;
use myxogast::matrix::*;


use rustc_serialize::json;
use std::ops::{Index,IndexMut};
use std::clone::Clone;
use std::fmt;
use std::fmt::Debug;
use std::cmp::{PartialOrd,Ordering,max};





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

    assert_eq!( Cell::unpack( &Cell::pack( &AlnState::Ins, &-10 ) ).unwrap(), (AlnState::Ins, -10) );
    assert_eq!( Cell::unpack( &Cell::pack( &AlnState::Match, &222 ) ).unwrap(), (AlnState::Match, 222) );

    let params = AlnParams {
        llocal:    false,
        rlocal:    false,
        max_indel: None,
        gap_open:  -1,
        gap_ext:   -1,
        mismatch:  -1,
        equal:     1 };

    let s = "ATGCATGC";
    assert_eq!( format!("{}", Sequence::from_str(s).unwrap()), s);

    let reference = Sequence::from_str("AAAAATGCTCGAAAAAAAA").unwrap();
    let query = Sequence::from_str("TGCTCG").unwrap();

    let (r, q) = align( &reference, &query, &params ).unwrap();
    println!("R: {}\nQ: {}\n", reference, query);


    let ref2 = Sequence::from_str("ATGCAT").unwrap();
    let q2 = Sequence::from_str("ATGCA").unwrap();

    let (r2, q2) = align( &ref2, &q2, &params ).unwrap();
    println!("R: {}\nQ: {}", r2, q2);

    /*
    let x_scores = Matrix { width: x.width, height: x.height,
                            data: x.data.iter().map( |c| Cell::unpack(c.clone()).unwrap().1 ).collect() };
    println!("\n{:?}", x_scores );
     */
}
