/*
 *                             The MIT License
 *
 * Wavefront Alignments Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignments Algorithms.
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
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Gap-affine alignment algorithms wrapper (including WFA)
 */

#include "benchmark/benchmark_check.h"
#include "benchmark/benchmark_gap_affine.h"
#include "gap_affine/affine_matrix.h"
#include "gap_affine/swg.h"
#include "wavefront/wavefront_align.h"

/*
 * Benchmark SWG
 */
void benchmark_gap_affine_swg(
    align_input_t* const align_input,
    affine_penalties_t* const penalties) {
  // Allocate
  affine_matrix_t affine_matrix;
  affine_matrix_allocate(
      &affine_matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t cigar;
  cigar_allocate(&cigar,
      align_input->pattern_length+align_input->text_length,
      align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  swg_compute(&affine_matrix,penalties,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,&cigar);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  if (align_input->output_file) {
    const int score = cigar_score_gap_affine(&cigar,penalties);
    benchmark_print_alignment_short(align_input->output_file,score,&cigar);
  }
  // Free
  affine_matrix_free(&affine_matrix,align_input->mm_allocator);
  cigar_free(&cigar);
}
void benchmark_gap_affine_swg_endsfree(
    align_input_t* const align_input,
    affine_penalties_t* const penalties,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free) {
  // Allocate
  affine_matrix_t affine_matrix;
  affine_matrix_allocate(
      &affine_matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t cigar;
  cigar_allocate(&cigar,
      align_input->pattern_length+align_input->text_length,
      align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  swg_compute_endsfree(&affine_matrix,penalties,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      pattern_begin_free,pattern_begin_free,
      text_begin_free,text_end_free,&cigar);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  if (align_input->output_file) {
    const int score = cigar_score_gap_affine(&cigar,penalties);
    benchmark_print_alignment_short(align_input->output_file,score,&cigar);
  }
  // Free
  affine_matrix_free(&affine_matrix,align_input->mm_allocator);
  cigar_free(&cigar);
}
void benchmark_gap_affine_swg_banded(
    align_input_t* const align_input,
    affine_penalties_t* const penalties,
    const int bandwidth) {
  // Allocate
  affine_matrix_t affine_matrix;
  affine_matrix_allocate(
      &affine_matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t cigar;
  cigar_allocate(&cigar,
      align_input->pattern_length+align_input->text_length,
      align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  swg_compute_banded(&affine_matrix,penalties,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,bandwidth,&cigar);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  if (align_input->output_file) {
    const int score = cigar_score_gap_affine(&cigar,penalties);
    benchmark_print_alignment_short(align_input->output_file,score,&cigar);
  }
  // Free
  affine_matrix_free(&affine_matrix,align_input->mm_allocator);
  cigar_free(&cigar);
}
void benchmark_gap_affine_wavefront(
    align_input_t* const align_input,
    affine_penalties_t* const penalties) {
  // Clear & resize
  wavefront_aligner_t* const wf_aligner = align_input->wf_aligner;
  wavefront_aligner_resize(wf_aligner,
      align_input->pattern_length,align_input->text_length);
  // Align
  timer_start(&align_input->timer);
  wavefront_align(wf_aligner,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&wf_aligner->cigar);
  }
  if (align_input->output_file) {
    const int score_only = (wf_aligner->alignment_scope == compute_score);
    const int score = (score_only) ?
        wf_aligner->cigar.score : cigar_score_gap_affine(&wf_aligner->cigar,penalties);
    benchmark_print_alignment_short(align_input->output_file,score,&wf_aligner->cigar);
  }
}
