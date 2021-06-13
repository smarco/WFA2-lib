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
 * DESCRIPTION: Benchmark utils
 */

#include "benchmark/benchmark_utils.h"
#include "alignment/score_matrix.h"
#include "edit/edit_dp.h"
#include "gap_lineal/nw.h"
#include "gap_affine/affine_matrix.h"
#include "gap_affine/swg.h"

/*
 * Setup
 */
void benchmark_align_input_clear(
    align_input_t* const align_input) {
  // Accuracy Stats
  counter_reset(&(align_input->align));
  counter_reset(&(align_input->align_correct));
  counter_reset(&(align_input->align_score));
  counter_reset(&(align_input->align_score_total));
  counter_reset(&(align_input->align_score_diff));
  counter_reset(&(align_input->align_cigar));
  counter_reset(&(align_input->align_bases));
  counter_reset(&(align_input->align_matches));
  counter_reset(&(align_input->align_mismatches));
  counter_reset(&(align_input->align_del));
  counter_reset(&(align_input->align_ins));
}
/*
 * Display
 */
void benchmark_print_alignment(
    FILE* const stream,
    align_input_t* const align_input,
    const int score_computed,
    cigar_t* const cigar_computed,
    const int score_correct,
    cigar_t* const cigar_correct) {
  // Print Sequence
  fprintf(stream,"ALIGNMENT (#%d)\n",align_input->sequence_id);
  fprintf(stream,"  PATTERN  %s\n",align_input->pattern);
  fprintf(stream,"  TEXT     %s\n",align_input->text);
  // Print CIGARS
  if (cigar_computed != NULL && score_computed != -1) {
    fprintf(stream,"    COMPUTED\tscore=%d\t",score_computed);
    cigar_print(stream,cigar_computed,true);
    fprintf(stream,"\n");
  }
  if (cigar_computed != NULL) {
    cigar_print_pretty(stream,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,
        cigar_computed,align_input->mm_allocator);
  }
  if (cigar_correct != NULL && score_correct != -1) {
    fprintf(stream,"    CORRECT \tscore=%d\t",score_correct);
    cigar_print(stream,cigar_correct,true);
    fprintf(stream,"\n");
  }
  if (cigar_correct != NULL) {
    cigar_print_pretty(stream,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,
        cigar_correct,align_input->mm_allocator);
  }
}
void benchmark_print_alignment_short(
    FILE* const stream,
    const int score,
    cigar_t* const cigar) {
  fprintf(stream,"%d ",score);
  cigar_print(stream,cigar,true);
  fprintf(stream,"\n");
}
/*
 * Stats
 */
void benchmark_print_stats(
    FILE* const stream,
    align_input_t* const align_input,
    const bool print_wf_stats) {
  // General stats
  fprintf(stream,"[Accuracy]\n");
  fprintf(stream," => Alignments.Correct     ");
  counter_print(stream,&align_input->align_correct,&align_input->align,"alg       ",true);
  fprintf(stream," => Score.Correct          ");
  counter_print(stream,&align_input->align_score,&align_input->align,"alg       ",true);
  fprintf(stream,"   => Score.Total          ");
  counter_print(stream,&align_input->align_score_total,NULL,"score uds.",true);
  fprintf(stream,"     => Score.Diff         ");
  counter_print(stream,&align_input->align_score_diff,&align_input->align_score_total,"score uds.",true);
  fprintf(stream," => CIGAR.Correct          ");
  counter_print(stream,&align_input->align_cigar,&align_input->align,"alg       ",true);
  fprintf(stream,"   => CIGAR.Matches        ");
  counter_print(stream,&align_input->align_matches,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Mismatches     ");
  counter_print(stream,&align_input->align_mismatches,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Insertions     ");
  counter_print(stream,&align_input->align_ins,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Deletions      ");
  counter_print(stream,&align_input->align_del,&align_input->align_bases,"bases     ",true);
}






