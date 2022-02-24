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
 * DESCRIPTION: WFAlm library wrapper
 */
#include "benchmark_wfalm.h"

/*
 * Benchmark WFAlm Bindings
 */
void benchmark_wfalm_bridge_global_affine(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int mismatch,
    const int gap_opening,
    const int gap_extension,
    char* const cigar_operations,
    int* const num_cigar_operations);
void benchmark_wfalm_bridge_global_affine_lowmem(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int mismatch,
    const int gap_opening,
    const int gap_extension,
    char* const cigar_operations,
    int* const num_cigar_operations);

/*
 * Benchmark WFAlm
 */
void benchmark_wfalm_global_affine(
    align_input_t* const align_input,
    affine_penalties_t* const penalties) {
  // Parameters
  const int max_cigar_length = align_input->pattern_length + align_input->text_length;
  cigar_t cigar = {
      .operations = malloc(max_cigar_length),
      .begin_offset = 0,
      .end_offset = 0,
  };
  // Align
  timer_start(&align_input->timer);
  benchmark_wfalm_bridge_global_affine(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      penalties->mismatch,penalties->gap_opening,
      penalties->gap_extension,
      cigar.operations,&cigar.end_offset);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,gap_affine,false,&cigar);
  }
  // Free
  free(cigar.operations);
}
void benchmark_wfalm_global_affine_lowmem(
    align_input_t* const align_input,
    affine_penalties_t* const penalties) {
  // Parameters
  const int max_cigar_length = align_input->pattern_length + align_input->text_length;
  cigar_t cigar = {
      .operations = malloc(max_cigar_length),
      .begin_offset = 0,
      .end_offset = 0,
  };
  // Align
  timer_start(&align_input->timer);
  benchmark_wfalm_bridge_global_affine_lowmem(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      penalties->mismatch,penalties->gap_opening,
      penalties->gap_extension,
      cigar.operations,&cigar.end_offset);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,gap_affine,false,&cigar);
  }
  // Free
  free(cigar.operations);
}
