/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Block-aligner library wrapper
 */
#include "benchmark_wfalm.h"
#include "external/block-aligner/c/block_aligner.h"

/*
 * Adapt CIGAR
 */
//void benchmark_blockaligner_adapt_cigar(
//    char* const pattern,
//    const int pattern_length,
//    char* const text,
//    const int text_length,
//    Cigar* const cigar,
//    const size_t cigar_len,
//    char* const cigar_operations,
//    int* const num_cigar_operations) {
//  // Note: 'M' signals either a match or mismatch
//  char ops_char[] = {' ', 'M', 'I', 'D'};
//  // Traverse all CIGAR ops
//  int pattern_pos = 0;
//  int text_pos = 0;
//  int cigar_idx = 0, j, i;
//  for (j=0;j<cigar_len;j++) {
//    OpLen o = block_get_cigar(cigar,j);
//    switch (ops_char[o.op]) {
//      case 'D':
//        text_pos += o.len;
//        for (i=0;i<o.len;++i) cigar_operations[cigar_idx++] = 'I';
//        break;
//      case 'I':
//        pattern_pos += o.len;
//        for (i=0;i<o.len;++i) cigar_operations[cigar_idx++] = 'D';
//        break;
//      case 'M':
//        for (i=0;i<o.len;++i,++text_pos,++pattern_pos) {
//          if (pattern[pattern_pos] != text[text_pos]) {
//            cigar_operations[cigar_idx++] = 'X';
//          } else {
//            cigar_operations[cigar_idx++] = 'M';
//          }
//        }
//        break;
//      default:
//        break;
//    }
//  }
//  // Return
//  *num_cigar_operations = cigar_idx;
//}
/*
 * Benchmark Block-aligner
 */
void benchmark_blockaligner_global_affine(
    align_input_t* const align_input,
    affine_penalties_t* const penalties,
    const int block_size) {
//  // Configure global alignment with traceback
//  SizeRange range = {.min = block_size, .max = block_size};
//  Gaps gaps = {.open = -penalties->gap_opening, .extend = -penalties->gap_extension};
//  AAMatrix* m_matrix = block_new_simple_aamatrix(0,-penalties->mismatch);
//  PaddedBytes* a = block_new_padded_aa(align_input->pattern_length, range.max);
//  PaddedBytes* b = block_new_padded_aa(align_input->text_length, range.max);
//  block_set_bytes_padded_aa(a,(const uint8_t*)align_input->pattern,align_input->pattern_length,range.max);
//  block_set_bytes_padded_aa(b,(const uint8_t*)align_input->text,align_input->text_length, range.max);
//  // Align
//  timer_start(&align_input->timer);
//  BlockHandle block = block_new_aa_trace(align_input->pattern_length,align_input->text_length,range.max);
//  block_align_aa_trace(block,a,b,m_matrix,gaps,range,0);
//  AlignResult res = block_res_aa_trace(block);
//  timer_stop(&align_input->timer);
//  // Allocate CIGAR
//  const int max_cigar_length = align_input->pattern_length + align_input->text_length;
//  cigar_t cigar = {
//      .operations = malloc(max_cigar_length),
//      .begin_offset = 0,
//      .end_offset = 0,
//  };
//  // Adapt CIGAR
//  Cigar* ba_cigar = block_new_cigar(res.query_idx, res.reference_idx);
//  block_cigar_aa_trace(block, res.query_idx,res.reference_idx,ba_cigar);
//  size_t ba_cigar_len = block_len_cigar(ba_cigar);
//  benchmark_blockaligner_adapt_cigar(
//      align_input->pattern,align_input->pattern_length,
//      align_input->text,align_input->text_length,
//      ba_cigar,ba_cigar_len,
//      cigar.operations,&cigar.end_offset);
//  // DEBUG
//  if (align_input->debug_flags) {
//    benchmark_check_alignment(align_input,&cigar);
//  }
//  // Output
//  if (align_input->output_file) {
//    benchmark_print_output(align_input,gap_affine,false,&cigar);
//  }
//  // Free
//  block_free_cigar(ba_cigar);
//  block_free_aa_trace(block);
//  block_free_padded_aa(a);
//  block_free_padded_aa(b);
//  free(cigar.operations);
}

