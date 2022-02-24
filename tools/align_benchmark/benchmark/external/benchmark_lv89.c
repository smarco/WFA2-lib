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

#include "external/lv89/lv89.c"
#include "external/lv89/lv89.h"

/*
 * Benchmark LV89
 */
void benchmark_lv89(
    align_input_t* const align_input) {
  // Allocate
  uint8_t* mem = (uint8_t*)malloc((align_input->pattern_length+align_input->text_length)*16);
  // Align
  timer_start(&align_input->timer);
  const int score = lv_ed(
      align_input->pattern_length,align_input->pattern,
      align_input->text_length,align_input->text,
      false,mem);
  timer_stop(&align_input->timer);
  // Free
  free(mem);
  // Output. NOTE: No CIGAR is produced, just score
  cigar_t cigar = {.begin_offset=0,.end_offset=0,.score=score};
  if (align_input->output_file) {
    benchmark_print_output(align_input,edit,true,&cigar);
  }
}
