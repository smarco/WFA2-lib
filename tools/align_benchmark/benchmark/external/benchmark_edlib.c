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

#include "external/edlib/edlib/include/edlib.h"

/*
 * Benchmark EdLib
 */
void benchmark_edlib_adapt_cigar(
    align_input_t* const align_input,
    char* const edlib_cigar,
    cigar_t* const cigar) {
  // Decode all CIGAR operations
  const int cigar_length = strlen(edlib_cigar);
  int chars_read = 0;
  while (chars_read < cigar_length) {
    // Read operation
    int length, n;
    char operation;
    sscanf(edlib_cigar+chars_read,"%d%c%n",&length,&operation,&n);
    chars_read += n;
    // Adapt operation encoding
    if (operation=='=') operation = 'M';
    else if (operation=='D') operation = 'I';
    else if (operation=='I') operation = 'D';
    // Dump operation
    int i;
    for (i=0;i<length;++i) {
      cigar->operations[(cigar->end_offset)++] = operation;
    }
  }
}
void benchmark_edlib(align_input_t* const align_input) {
  // Parameters
  EdlibAlignResult result;
  char* edlib_cigar = NULL;
  // Align
  timer_start(&align_input->timer);
  result = edlibAlign(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      edlibNewAlignConfig(-1,EDLIB_MODE_NW,EDLIB_TASK_PATH,NULL,0));
  edlib_cigar = edlibAlignmentToCigar(
      result.alignment,result.alignmentLength,EDLIB_CIGAR_EXTENDED); // Traceback
  timer_stop(&align_input->timer);
  // Adapt CIGAR
  cigar_t cigar;
  cigar.operations = malloc(align_input->pattern_length+align_input->text_length);
  cigar.begin_offset = 0;
  cigar.end_offset = 0;
  benchmark_edlib_adapt_cigar(align_input,edlib_cigar,&cigar);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,edit,false,&cigar);
  }
  // Free
  free(edlib_cigar);
  free(cigar.operations);
  edlibFreeAlignResult(result);
}
