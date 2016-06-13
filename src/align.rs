use std::fmt;
use std::fmt::Debug;
use std::cmp::min;
use matrix::Matrix;
use seq::*;

#[derive(Clone)]
pub struct Cell(pub i32);

/// Cell: bit-pack maximal alignment state and score into 32 bits
///
/// note that a simple Smith-Waterman doesn't require storing the "direction" of the
/// previous max in the cell, however, we need it to implement separate gap-opening
/// and gap-extension penalties
///
impl Cell {
    // bit-packing operations
    // FIXME: i16 might not be enough, so I should investigate a different division, eg
    //        4 bits for AlnState and 28 bits get munged into an i32
    pub fn pack( state : &AlnState, score : &i16 ) -> Cell {
        let _state : i32 = match *state {
            AlnState::Nil =>  0xffff,
            AlnState::Match => (1 << 16) | 0xffff,
            AlnState::Mismatch => (2 << 16) | 0xffff,
            AlnState::Ins => (3 << 16) | 0xffff,
            AlnState::Del => (4 << 16) | 0xffff };
        let _score : i32 = (*score as i32) | 0xffff0000;
        Cell(_score & _state)
    }

    pub fn unpack( c : &Cell ) -> Result<(AlnState, i16), ()> {
        let Cell(packed) = *c;
        let state = try!(
            match packed >> 16 {
                0 => Ok(AlnState::Nil),
                1 => Ok(AlnState::Match),
                2 => Ok(AlnState::Mismatch),
                3 => Ok(AlnState::Ins),
                4 => Ok(AlnState::Del),
                _ => Err(()) }
        );
        Ok( (state, packed as i16) )
    }
}

impl fmt::Debug for Cell {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (a, b) = Cell::unpack( &self ).unwrap();
        write!(f, "{:?}({})", a, b)
    }
}

// used to store HMMer-style probabilities; dimensions: BASE-COUNT x REF-LENGTH
pub type ProbMatr = Matrix<f32>;

pub struct AlnParams {
    pub llocal:     bool,        // global or local on the left?
    pub rlocal:     bool,        // global or local on the right?
    pub max_indel:  Option<u8>,  // maximum number of indels before short-circuiting
    pub gap_open:   i16,         // penalty for opening a gap
    pub gap_ext:    i16,         // gap extention penalty
    pub mismatch:   i16,         // mismatch penalty
    pub equal:      i16,         // match score, but "match" is a keyword
}

impl AlnParams {
    pub fn new(  _llocal: Option<bool>, _rlocal: Option<bool>, _max_indel: Option<u8>,
                 _gap_open: Option<i16>, _gap_ext: Option<i16>, _mismatch: i16, _equal: Option<i16> ) -> AlnParams {
        AlnParams {
            llocal:     match _llocal { Some(b) => b, None => true },
            rlocal:     match _rlocal { Some(b) => b, None => true },
            max_indel:  _max_indel,
            gap_open:   match _gap_open { Some(g) => g, None => -1 },
            gap_ext:    match _gap_ext { Some(g) => g, None => -1 },
            mismatch:   _mismatch,
            equal:      match _equal { Some(g) => g, None => 1 }, // "match" is a keyword
        }
    }
}

#[derive(Debug,Clone)]
pub enum AlnState {
    Nil,
    Match,
    Mismatch,
    Ins,
    Del
}

// NOTE: I'm a little shocked I can't just #derive this...
impl PartialEq for AlnState {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (&AlnState::Nil, &AlnState::Nil) => true,
            (&AlnState::Match, &AlnState::Match) => true,
            (&AlnState::Mismatch, &AlnState::Mismatch) => true,
            (&AlnState::Ins, &AlnState::Ins) => true,
            (&AlnState::Del, &AlnState::Del) => true,
            _ => false }
    }
}
impl Eq for AlnState {} 

/// Dynamic-programming alignment
/// params can specify a max_indel, in which case this can return None if a
///   solution can't be found with fewer indels
///
/// returns Option(m, i, j), where m is the alignment matrix, and i,j is the
///   the location of the cell with the highest value
///
fn align_matrix( reference: &Sequence, query: &Sequence, params: &AlnParams ) -> Option<(Matrix<Cell>, i32, i32)> {
    let mut m = Matrix::<Cell>::new( Cell(0), reference.len() + 2, query.len() + 2 );
    let ref_len : i32 = reference.len() as i32;
    let query_len : i32 = query.len() as i32;

    let mut best_val : i16 = i16::min_value();
    let mut best_loc : (i32, i32) = (0,0);

    // initialize edges
    if !params.llocal {
        for i in 0 .. ref_len + 1 { m[ (i as i32, 0) ] = Cell::pack( &AlnState::Nil, &(-i as i16)); }
    }
    for j in 0i16 .. (query_len + 1) as i16 { m[ (0, j as i32) ] = Cell::pack( &AlnState::Nil, &(-j as i16)); }

    for i in 1 .. ref_len + 1 {
        for j in 1 .. query_len + 1 {

            let (dstate, del) = Cell::unpack( &m[ (i-1, j) ] ).unwrap();
            let del_score = del + match (params.rlocal && i == ref_len, dstate) {
                (true, _) => 0,
                (_, AlnState::Del) => params.gap_ext,
                _ => params.gap_open
            };

            let (istate, ins) = Cell::unpack( &m[ (i, j-1) ] ).unwrap();
            let ins_score = ins + if istate == AlnState::Ins { params.gap_ext } else { params.gap_open };

            let (_, diag) = Cell::unpack( &m[ (i-1, j-1) ] ).unwrap();
            let diag_score = diag + if reference[i-1] == query[j-1] { params.equal } else { params.mismatch };

            let (a,b) = {
                if diag_score >= del_score && diag_score >= ins_score {
                    (if reference[i-1] == query[j-1] {AlnState::Match} else {AlnState::Mismatch}, diag_score)
                } else if del_score > diag_score && del_score >= ins_score {
                    (AlnState::Del, del_score)
                } else {
                    (AlnState::Ins, ins_score)
                }};

            if b > best_val {
                best_val = b;
                best_loc = (i, j);
            }

            m[ (i, j) ] = Cell::pack( &a, &b );
        }
    };

    println!("-- finished:\n");
    println!("{:?}\n", Matrix { width: m.width, height: m.height,
                                data: m.data.iter().map( |c| Cell::unpack(&c).unwrap().1 ).collect() });
    println!("best_loc: {},{}", best_loc.0, best_loc.1);
    Some((m, best_loc.0, best_loc.1))
}

pub fn align_hmm( reference: ProbMatr, query: Sequence, params: AlnParams ) -> Option<()> {
    //let mut m = Matrix::<f32>::new( 0., reference.width + 1, query.len() + 1 );
    None
}

fn aln_from_coord( st_i : &i32, st_j : &i32, _inc : &i32, reference : &Sequence,
                   query : &Sequence, alignment : &Matrix<Cell> )
                   -> (Sequence, Sequence) {

    let ref_len : i32 = reference.len() as i32;
    let query_len : i32 = query.len() as i32;

    let mut padded_ref : Vec<Mmer> = Vec::with_capacity(reference.len());
    let mut padded_query : Vec<Mmer> = Vec::with_capacity(query.len());
    let inc = *_inc;
    let mut i = *st_i;
    let mut j = *st_j;

    let max_i = if inc < 0 { ref_len + 1 } else { ref_len };
    let max_j = if inc < 0 { query_len + 1 } else { query_len };

    while i < max_i && i > 1 && j < max_j && j > 1 {

        let (_, m_sc) = Cell::unpack( &alignment[ (i+inc, j+inc) ] ).unwrap();
        let (_, r_sc) = Cell::unpack( &alignment[ (i+inc, j) ] ).unwrap();
        let (_, d_sc) = Cell::unpack( &alignment[ (i, j+inc) ] ).unwrap();

        if m_sc >= r_sc && m_sc >= d_sc {
            i += inc;
            j += inc;
            padded_ref.push( reference[i-1] );
            padded_query.push( query[j-1] );

        } else if r_sc >= m_sc && r_sc >= d_sc {
            i += inc;
            padded_ref.push( reference[i-1] );
            padded_query.push(HYPHEN);

        } else {
            j += inc;
            padded_ref.push(HYPHEN);
            padded_query.push( query[j-1] );
        }
    }

    while i < ref_len && i > 1 {
        i += inc;
        padded_ref.push( reference[i-1] );
        padded_query.push(HYPHEN);
    }

    while j < query_len && j > 1 {
        j += inc;
        padded_ref.push(HYPHEN);
        padded_query.push( query[j-1] );
    }

    (Sequence(padded_ref), Sequence(padded_query))
}

pub fn align( reference: &Sequence, query: &Sequence, params: &AlnParams )
              -> Option<(Sequence, Sequence)> {
    match align_matrix( reference, query, params ) {
        Some( (m, st_i, st_j) ) => {
            println!("val @ max: {}", Cell::unpack( &m[(st_i as i32, st_j as i32)] ).unwrap().1 );
            for j in 0 .. m.height {
                let mut disp : Vec<i16> = Vec::with_capacity(m.width);
                for i in 0 .. min( m.width, 45 ) {
                    disp.push( Cell::unpack( &m[(i as i32, j as i32)] ).unwrap().1 );
                }
                println!("X: {:?}", disp);
            }

            let (fwd_r, fwd_q) = aln_from_coord( &st_i, &st_j, &1, &reference, &query, &m );
            let (rev_r1, rev_q1) = aln_from_coord( &st_i, &st_j, &(-1), &reference, &query, &m );
            let mut rev_r2 = rev_r1.reverse().0;
            if st_i <= reference.len() as i32 {
                rev_r2.push( reference[st_i-1] );
            }

            let mut rev_q2 = rev_q1.reverse().0;
            if st_j <= query.len() as i32 {
                rev_q2.push( query[st_j-1] );
            }
            
            Some( (Sequence(rev_r2) + fwd_r, Sequence(rev_q2) + fwd_q) ) },

        None => None
    }
}
