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
 * DESCRIPTION: BitPal wrapper
 */

#include "benchmark_bitpal.h"

/*
 * BitPal
 */
int bitwise_alignment_m0_x1_g1(char * s1,char* s2,int words);
void benchmark_bitpal_m0_x1_g1(
    align_input_t* const align_input) {
  // Align
  timer_start(&align_input->timer);
  bitwise_alignment_m0_x1_g1(
      align_input->pattern,align_input->text,
      MAX(align_input->pattern_length,align_input->text_length)/63 + 1);
  timer_stop(&align_input->timer);
  // NOTE:
  //   No CIGAR is produced, just score
}
int bitwise_alignment_m1_x4_g2(char * s1,char* s2,int words);
void benchmark_bitpal_m1_x4_g2(
    align_input_t* const align_input) {
  // Align
  timer_start(&align_input->timer);
  bitwise_alignment_m1_x4_g2(
      align_input->pattern,align_input->text,
      MAX(align_input->pattern_length,align_input->text_length)/63 + 1);
  timer_stop(&align_input->timer);
  // NOTE:
  //   No CIGAR is produced, just score
}
