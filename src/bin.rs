extern crate rustc_serialize;
extern crate myxogast;
extern crate bio;

extern crate argparse;

use std::io;
use std::fs::File;
use rustc_serialize::json;
use rustc_serialize::json::Json;
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
 
use argparse::{ArgumentParser, StoreTrue, Store};


fn main() {
    let mut ref_fname : String = String::new();
    let mut query_fname : String = String::new();
    { // scope block?
        let mut parser = ArgumentParser::new();
        parser.refer(&mut ref_fname)
            .add_argument("ref", Store, "refernce JSON file")
            .required();
        parser.refer(&mut query_fname)
            .add_argument("query", Store, "query file")
            .required();
        parser.parse_args_or_exit();
    }
    println!("fname: {}", ref_fname);
    println!("fname: {}", query_fname);

    //let graph = None;

    let reader = fasta::Reader::from_file(query_fname).unwrap();
    for gene in reader.records() {
        let seq = gene.unwrap();
 
        // FIXME: ugh ..
        let seq_str : String = String::from_utf8_lossy( &seq.seq() ).into_owned();
 
        let query = Sequence::from_str(seq_str.as_str()).unwrap();
    }
}
