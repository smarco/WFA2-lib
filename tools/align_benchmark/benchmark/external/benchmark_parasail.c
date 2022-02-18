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
 * DESCRIPTION: Parasail library wrapper
 */

#include "benchmark_parasail.h"
#include "benchmark/benchmark_check.h"
#include "external/parasail/parasail.h"

/*
 * Parasail
 */
void benchmark_parasail_parse_cigar(
    align_input_t* const align_input,
    affine_penalties_t* const penalties,
    parasail_result_t* const parasail_result,
    parasail_cigar_t* const parasail_cigar,
    cigar_t* const cigar) {
  // Allocate
  cigar->operations = malloc(align_input->pattern_length+align_input->text_length);
  cigar->begin_offset = 0;
  cigar->end_offset = 0;
  // Decode all CIGAR operations
  int i;
  for (i=0;i<parasail_cigar->len;++i) {
    // Fetch operation and length
    char operation = parasail_cigar_decode_op(parasail_cigar->seq[i]);
    const uint32_t length = parasail_cigar_decode_len(parasail_cigar->seq[i]);
    // Adapt operation '='
    if (operation == '=') operation = 'M';
    else if (operation == 'D') operation = 'I';
    else if (operation == 'I') operation = 'D';
    // Dump into the cigar
    int j;
    for (j=0;j<length;++j) {
      cigar->operations[(cigar->end_offset)++] = operation;
    }
    cigar->operations[cigar->end_offset] = '\0';
  }
//  // DEBUG
//  const int score_parasail = parasail_result_get_score(parasail_result);
//  const int score_cigar = cigar_score_gap_affine(cigar,penalties);
//  if (score_parasail != score_cigar) {
//    char* const cigar_str = parasail_cigar_decode(parasail_cigar);
//    fprintf(stderr,"COMPUTED=%d\t",score_cigar);
//    cigar_print(stderr,cigar);
//    fprintf(stderr,"\nPARASAIL=%d\t%s\n",score_parasail,cigar_str);
//    free(cigar_str);
//  }
}
void benchmark_parasail_nw_stripped(
    align_input_t* const align_input,
    affine_penalties_t* const penalties) {
  // Create parasail-matrix
  parasail_matrix_t* const matrix =
      parasail_matrix_create("ACGT", -penalties->match, -penalties->mismatch);
  parasail_result_t* result = NULL;
  // Align
  timer_start(&align_input->timer);
  result = parasail_nw_trace_striped_avx2_256_16(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      penalties->gap_opening + penalties->gap_extension,
      penalties->gap_extension,matrix);
  // Traceback
  parasail_cigar_t* parasail_cigar = parasail_result_get_cigar(
      result,align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,matrix);
  timer_stop(&align_input->timer);
  // Adapt CIGAR
  cigar_t cigar;
  benchmark_parasail_parse_cigar(
      align_input,penalties,result,parasail_cigar,&cigar);
  // Debug alignment
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Free
  free(cigar.operations);
  parasail_result_free(result);
  parasail_matrix_free(matrix);
}
void benchmark_parasail_nw_scan(
    align_input_t* const align_input,
    affine_penalties_t* const penalties) {
  // Create parasail-matrix
  parasail_matrix_t* const matrix =
      parasail_matrix_create("ACGT", -penalties->match, -penalties->mismatch);
  parasail_result_t* result = NULL;
  // Align
  timer_start(&align_input->timer);
  result = parasail_nw_trace_scan_avx2_256_16(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      penalties->gap_opening + penalties->gap_extension,
      penalties->gap_extension,matrix);
  // Traceback
  parasail_cigar_t* parasail_cigar = parasail_result_get_cigar(
      result,align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,matrix);
  timer_stop(&align_input->timer);
  // Adapt CIGAR
  cigar_t cigar;
  benchmark_parasail_parse_cigar(
      align_input,penalties,result,
      parasail_cigar,&cigar);
  // Debug alignment
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Free
  free(cigar.operations);
  parasail_result_free(result);
  parasail_matrix_free(matrix);
}
void benchmark_parasail_nw_diag(
    align_input_t* const align_input,
    affine_penalties_t* const penalties) {
  // Create parasail-matrix
  parasail_matrix_t* const matrix =
      parasail_matrix_create("ACGT", -penalties->match, -penalties->mismatch);
  parasail_result_t* result = NULL;
  // Align
  timer_start(&align_input->timer);
  result = parasail_nw_trace_diag_avx2_256_16(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      penalties->gap_opening + penalties->gap_extension,
      penalties->gap_extension,matrix);
  // Traceback
  parasail_cigar_t* parasail_cigar = parasail_result_get_cigar(
      result,align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,matrix);
  timer_stop(&align_input->timer);
  // Adapt CIGAR
  cigar_t cigar;
  benchmark_parasail_parse_cigar(
      align_input,penalties,result,
      parasail_cigar,&cigar);
  // Debug alignment
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Free
  free(cigar.operations);
  parasail_result_free(result);
  parasail_matrix_free(matrix);
}
void benchmark_parasail_nw_banded(
    align_input_t* const align_input,
    affine_penalties_t* const penalties,
    const int bandwidth) {
  // Create parasail-matrix
  parasail_matrix_t* const matrix =
      parasail_matrix_create("ACGT", -penalties->match, -penalties->mismatch);
  parasail_result_t* result = NULL;
  // Align
  timer_start(&align_input->timer);
  result = parasail_nw_banded(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      penalties->gap_opening + penalties->gap_extension,
      penalties->gap_extension,bandwidth,matrix);
  timer_stop(&align_input->timer);
  // NOTE:
  //   No CIGAR is produced, just score
  // Free
  parasail_result_free(result);
  parasail_matrix_free(matrix);
}

