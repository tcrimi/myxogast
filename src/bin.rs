extern crate rustc_serialize;
extern crate myxogast;
extern crate bio;

use std::io;
use std::fs::File;
use rustc_serialize::json;
use std::ops::{Index,IndexMut};
use std::clone::Clone;
use std::fmt;
use std::str;
use std::fmt::Debug;
use std::cmp::{PartialOrd,Ordering,max};

use bio::io::fasta;

use myxogast::align::*;
use myxogast::tree::*;
use myxogast::seq::*;
use myxogast::matrix::*;



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

    println!("testing AlnState::Ins, -10");
    assert_eq!( Cell::unpack( &Cell::pack( &AlnState::Ins, &-10 ) ).unwrap(), (AlnState::Ins, -10) );
    println!("testing Match, 222");
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
    assert_eq!( r, reference );
    assert_eq!( q, Sequence::from_str("-----TGCTCG--------").unwrap() );

    let ref2 = Sequence::from_str("ATGCAT").unwrap();
    let query2 = Sequence::from_str("ATGCA").unwrap();
    let (r2, q2) = align( &ref2, &query2, &params ).unwrap();
    assert_eq!( r2, ref2 );
    assert_eq!( q2, Sequence::from_str("ATGCA-").unwrap());


    let params_loc = AlnParams {
        llocal:    true,
        rlocal:    true,
        max_indel: None,
        gap_open:  -1,
        gap_ext:   -1,
        mismatch:  -1,
        equal:     1 };

    // rust-bio
    let reader = fasta::Reader::from_file("e_coli.fasta").unwrap();
    for gene in reader.records() {
        let seq = gene.unwrap();

        // FIXME: ugh ..
        let seq_str : String = String::from_utf8_lossy( &seq.seq() ).into_owned();

        let e_coli = Sequence::from_str(seq_str.as_str()).unwrap();

        let sub_ec = Sequence::from_str("CGAAGTGTTTGTGATTGGCGTCGGTGGCGTTGGCGGTGCGCTGCTGGAGCAACTGAA").unwrap();

        let (ec_r, ec_q) = align( &e_coli, &sub_ec, &params_loc ).unwrap();
        println!("R: {}", ec_r);
        println!("Q: {}", ec_q);

        //let s2 = seq.to_vec();
        //println!("{}", Sequence(s2));
    }

    /*
    let x_scores = Matrix { width: x.width, height: x.height,
                            data: x.data.iter().map( |c| Cell::unpack(c.clone()).unwrap().1 ).collect() };
    println!("\n{:?}", x_scores );
     */
}
