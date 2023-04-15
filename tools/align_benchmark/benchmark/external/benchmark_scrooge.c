/*
 *  Wavefront Alignment Algorithms
 *  Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of Wavefront Alignment Algorithms.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Wavefront Alignment Algorithms Benchmark
 */

#include "benchmark/benchmark_utils.h"
#include "benchmark/benchmark_edit.h"
#include "benchmark/benchmark_check.h"

/*
 * Brige
 */
void benchmark_scrooge_bridge(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    char* const edit_operations,
    int* const num_edit_operations,
    uint64_t* const time_ns);

/*
 * Benchmark Scrooge
 */
void benchmark_scrooge(
    align_input_t* const align_input) {
  // Parameters
  const int max_cigar_length = align_input->pattern_length + align_input->text_length;
  cigar_t* const cigar = cigar_new(2*max_cigar_length);
  // Align
  uint64_t time_ns = 0;
  benchmark_scrooge_bridge(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      cigar->operations,&cigar->end_offset,&time_ns);
  counter_add(&align_input->timer.time_ns,time_ns);
  cigar_score_edit(cigar);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,edit,false,cigar);
  }
  // Free
  cigar_free(cigar);
}
