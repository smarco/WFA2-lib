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
 * DESCRIPTION: Dynamic-programming alignment algorithm for computing
 *   gap-lineal pairwise alignment (Needleman-Wunsch - NW)
 */

#include "gap_lineal/nw.h"

/*
 * NW Traceback
 */
void nw_traceback(
    score_matrix_t* const score_matrix,
    lineal_penalties_t* const penalties,
    cigar_t* const cigar) {
  // Parameters
  int** const dp = score_matrix->columns;
  char* const operations = cigar->operations;
  cigar->end_offset = cigar->max_operations;
  int op_sentinel = cigar->end_offset - 1;
  // Compute traceback
  int h = score_matrix->num_columns-1;
  int v = score_matrix->num_rows-1;
  while (h>0 && v>0) {
    if (dp[h][v] == dp[h][v-1]+penalties->deletion) {
      operations[op_sentinel--] = 'D';
      --v;
    } else if (dp[h][v] == dp[h-1][v]+penalties->insertion) {
      operations[op_sentinel--] = 'I';
      --h;
    } else {
      operations[op_sentinel--] =
          (dp[h][v] == dp[h-1][v-1]+penalties->mismatch) ? 'X' : 'M';
      --h;
      --v;
    }
  }
  while (h>0) {operations[op_sentinel--] = 'I'; --h;}
  while (v>0) {operations[op_sentinel--] = 'D'; --v;}
  cigar->begin_offset = op_sentinel + 1;
}
void nw_compute(
    score_matrix_t* const score_matrix,
    lineal_penalties_t* const penalties,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar) {
  // Parameters
  int** dp = score_matrix->columns;
  int h, v;
  // Init DP (No ends-free)
  dp[0][0] = 0;
  for (v=1;v<=pattern_length;++v) {
    dp[0][v] = dp[0][v-1] + penalties->deletion;
  }
  for (h=1;h<=text_length;++h) {
    dp[h][0] = dp[h-1][0] + penalties->insertion;
  }
  // Compute DP
  for (h=1;h<=text_length;++h) {
    for (v=1;v<=pattern_length;++v) {
      int min = dp[h-1][v-1] + ((pattern[v-1]==text[h-1]) ? 0 : penalties->mismatch); // Misms
      min = MIN(min,dp[h-1][v]+penalties->insertion); // Ins
      min = MIN(min,dp[h][v-1]+penalties->deletion); // Del
      dp[h][v] = min;
    }
  }
  // Compute traceback
  nw_traceback(score_matrix,penalties,cigar);
}
