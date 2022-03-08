//! C bindings for block aligner.
//!
//! Generics are monomorphised manually.
//!
//! Nucleotide and arbitrary byte alignment do not have bindings yet.

use std::ffi::c_void;

use crate::scan_block::*;
use crate::scores::*;
use crate::cigar::*;

// avoid generics by using void pointer and monomorphism
/// A handle for a block in block aligner.
pub type BlockHandle = *mut c_void;

/// Represents a range that has inclusive lower and upper bounds.
#[derive(Copy, Clone, PartialEq)]
#[repr(C)]
pub struct SizeRange {
    pub min: usize,
    pub max: usize
}


// AAMatrix

/// Create a new simple AAMatrix with custom match and mismatch scores.
///
/// Note that the match score must be positive and the mismatch score must be negative.
#[no_mangle]
pub unsafe extern fn block_new_simple_aamatrix(match_score: i8, mismatch_score: i8) -> *mut AAMatrix {
    let matrix = Box::new(AAMatrix::new_simple(match_score, mismatch_score));
    Box::into_raw(matrix)
}

/// Frees an AAMatrix.
#[no_mangle]
pub unsafe extern fn block_free_simple_aamatrix(matrix: *mut AAMatrix) {
    drop(Box::from_raw(matrix));
}


// CIGAR

/// Create a new empty CIGAR string.
#[no_mangle]
pub unsafe extern fn block_new_cigar(query_len: usize, reference_len: usize) -> *mut Cigar {
    let cigar = Box::new(Cigar::new(query_len, reference_len));
    Box::into_raw(cigar)
}

/// Get the operation at a certain index in a CIGAR string.
#[no_mangle]
pub unsafe extern fn block_get_cigar(cigar: *const Cigar, i: usize) -> OpLen {
    let cigar_str = &*cigar;
    cigar_str.get(i)
}

/// Get the length of a CIGAR string.
#[no_mangle]
pub unsafe extern fn block_len_cigar(cigar: *const Cigar) -> usize {
    let cigar_str = &*cigar;
    cigar_str.len()
}

/// Frees a CIGAR string.
#[no_mangle]
pub unsafe extern fn block_free_cigar(cigar: *mut Cigar) {
    drop(Box::from_raw(cigar));
}


// PaddedBytes

/// Create a new empty padded amino acid string.
#[no_mangle]
pub unsafe extern fn block_new_padded_aa(len: usize, max_size: usize) -> *mut PaddedBytes {
    let padded_bytes = Box::new(PaddedBytes::new::<AAMatrix>(len, max_size));
    Box::into_raw(padded_bytes)
}

/// Write to a padded amino acid string.
#[no_mangle]
pub unsafe extern fn block_set_bytes_padded_aa(padded: *mut PaddedBytes, s: *const u8, len: usize, max_size: usize) {
    let bytes = std::slice::from_raw_parts(s, len);
    let padded_bytes = &mut *padded;
    padded_bytes.set_bytes::<AAMatrix>(bytes, max_size);
}

/// Frees a padded amino acid string.
#[no_mangle]
pub unsafe extern fn block_free_padded_aa(padded: *mut PaddedBytes) {
    drop(Box::from_raw(padded));
}


// Block

macro_rules! gen_functions {
    ($new_name:ident, $new_doc:expr,
     $align_name:ident, $align_doc:expr,
     $res_name:ident, $res_doc:expr,
     $trace_name:ident, $trace_doc:expr,
     $free_name:ident, $free_doc:expr,
     $matrix:ty, $trace:literal, $x_drop:literal) => {
        #[doc = $new_doc]
        #[no_mangle]
        pub unsafe extern fn $new_name(query_len: usize,
                                       reference_len: usize,
                                       max_size: usize) -> BlockHandle {
            let aligner = Box::new(Block::<$trace, $x_drop>::new(query_len, reference_len, max_size));
            Box::into_raw(aligner) as BlockHandle
        }

        #[doc = $align_doc]
        #[no_mangle]
        pub unsafe extern fn $align_name(b: BlockHandle,
                                         q: *const PaddedBytes,
                                         r: *const PaddedBytes,
                                         m: *const $matrix,
                                         g: Gaps,
                                         s: SizeRange,
                                         x: i32) {
            let aligner = &mut *(b as *mut Block<$trace, $x_drop>);
            aligner.align(&*q, &*r, &*m, g, s.min..=s.max, x);
        }

        #[doc = $res_doc]
        #[no_mangle]
        pub unsafe extern fn $res_name(b: BlockHandle) -> AlignResult {
            let aligner = &*(b as *const Block<$trace, $x_drop>);
            aligner.res()
        }

        #[doc = $trace_doc]
        #[no_mangle]
        pub unsafe extern fn $trace_name(b: BlockHandle, query_idx: usize, reference_idx: usize, cigar: *mut Cigar) {
            let aligner = &*(b as *const Block<$trace, $x_drop>);
            aligner.trace().cigar(query_idx, reference_idx, &mut *cigar);
        }

        #[doc = $free_doc]
        #[no_mangle]
        pub unsafe extern fn $free_name(b: BlockHandle) {
            drop(Box::from_raw(b as *mut Block<$trace, $x_drop>));
        }
    };
}

gen_functions!(
    block_new_aa,
    "Create a new block aligner instance for global alignment of amino acid strings (no traceback).",
    block_align_aa,
    "Global alignment of two amino acid strings (no traceback).",
    block_res_aa,
    "Retrieves the result of global alignment of two amino acid strings (no traceback).",
    _block_cigar_aa,
    "Don't use.",
    block_free_aa,
    "Frees the block used for global alignment of two amino acid strings (no traceback).",
    AAMatrix, false, false
);

gen_functions!(
    block_new_aa_xdrop,
    "Create a new block aligner instance for X-drop alignment of amino acid strings (no traceback).",
    block_align_aa_xdrop,
    "X-drop alignment of two amino acid strings (no traceback).",
    block_res_aa_xdrop,
    "Retrieves the result of X-drop alignment of two amino acid strings (no traceback).",
    _block_cigar_aa_xdrop,
    "Don't use.",
    block_free_aa_xdrop,
    "Frees the block used for X-drop alignment of two amino acid strings (no traceback).",
    AAMatrix, false, true
);

gen_functions!(
    block_new_aa_trace,
    "Create a new block aligner instance for global alignment of amino acid strings, with traceback.",
    block_align_aa_trace,
    "Global alignment of two amino acid strings, with traceback.",
    block_res_aa_trace,
    "Retrieves the result of global alignment of two amino acid strings, with traceback.",
    block_cigar_aa_trace,
    "Retrieves the resulting CIGAR string from global alignment of two amino acid strings, with traceback.",
    block_free_aa_trace,
    "Frees the block used for global alignment of two amino acid strings, with traceback.",
    AAMatrix, true, false
);

gen_functions!(
    block_new_aa_trace_xdrop,
    "Create a new block aligner instance for X-drop alignment of amino acid strings, with traceback.",
    block_align_aa_trace_xdrop,
    "X-drop alignment of two amino acid strings, with traceback.",
    block_res_aa_trace_xdrop,
    "Retrieves the result of X-drop alignment of two amino acid strings, with traceback.",
    block_cigar_aa_trace_xdrop,
    "Retrieves the resulting CIGAR string from X-drop alignment of two amino acid strings, with traceback.",
    block_free_aa_trace_xdrop,
    "Frees the block used for X-drop alignment of two amino acid strings, with traceback.",
    AAMatrix, true, true
);
