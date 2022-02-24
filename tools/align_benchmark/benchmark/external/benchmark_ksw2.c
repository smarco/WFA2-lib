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
 * DESCRIPTION: KSW2 library wrapper
 */

#include "benchmark_ksw2.h"
#include "benchmark/benchmark_check.h"
#include "external/ksw2/ksw2.h"

/*
 * ksw2
 */
void benchmark_ksw2_adapt_cigar(
    align_input_t* const align_input,
    ksw_extz_t* const ez,
    cigar_t* const cigar) {
  // Adapt CIGAR
  cigar->operations = malloc(2*(align_input->text_length+align_input->pattern_length));
  cigar->begin_offset = 0;
  cigar->end_offset = 0;
  int i, j;
  for (i=0;i<ez->n_cigar;++i) {
    const int length = ez->cigar[i]>>4;
    const char operation = "MDI"[ez->cigar[i]&0xf];
    for (j=0;j<length;++j) {
      cigar->operations[(cigar->end_offset)++] = operation;
    }
  }
  cigar->operations[cigar->end_offset] = '\0';
  // Add mismatches
  cigar_add_mismatches(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,cigar);
}
/*
 * ksw2
 */
void benchmark_ksw2_extz2_sse(
    align_input_t* const align_input,
    affine_penalties_t* const penalties,
    const bool approximate_max__drop,
    const int band_width,
    const int zdrop) {
  // Parameters
  const char* tseq = align_input->text;
  const char* qseq = align_input->pattern;
  const int sc_mch = -penalties->match;
  const int sc_mis = -penalties->mismatch;
  const int gapo = penalties->gap_opening;
  const int gape = penalties->gap_extension;
  // Prepare data
  void *km = 0;
  int flag = 0;
  int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
  int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
  int tlen = strlen(tseq), qlen = strlen(qseq);
  uint8_t *ts, *qs, c[256];
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  memset(c, 4, 256);
  c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
  c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
  ts = (uint8_t*)malloc(tlen);
  qs = (uint8_t*)malloc(qlen);
  for (i = 0; i < tlen; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
  for (i = 0; i < qlen; ++i) qs[i] = c[(uint8_t)qseq[i]];
  if (approximate_max__drop) flag |= KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP;
  // Align
  timer_start(&align_input->timer);
  ksw_extz2_sse(
      km,qlen,(uint8_t*)qseq,tlen,(uint8_t*)tseq,5,
      mat,gapo,gape,band_width,zdrop,0,flag,&ez);
  timer_stop(&align_input->timer);
  // Adapt CIGAR
  cigar_t cigar;
  benchmark_ksw2_adapt_cigar(align_input,&ez,&cigar);
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
  free(ez.cigar);
  free(ts);
  free(qs);
}
void benchmark_ksw2_extd2_sse(
    align_input_t* const align_input,
    affine2p_penalties_t* const penalties,
    const bool approximate_max__drop,
    const int band_width,
    const int zdrop) {
  // Parameters
  const char *tseq = align_input->text;
  const char *qseq = align_input->pattern;
  const int sc_mch = -penalties->match;
  const int sc_mis = -penalties->mismatch;
  const int gapo1 = penalties->gap_opening1;
  const int gape1 = penalties->gap_extension1;
  const int gapo2 = penalties->gap_opening2;
  const int gape2 = penalties->gap_extension2;
  // Prepare data
  void *km = 0;
  int flag = 0;
  int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
  int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
  int tlen = strlen(tseq), qlen = strlen(qseq);
  uint8_t *ts, *qs, c[256];
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  memset(c, 4, 256);
  c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
  c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
  ts = (uint8_t*)malloc(tlen);
  qs = (uint8_t*)malloc(qlen);
  for (i = 0; i < tlen; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
  for (i = 0; i < qlen; ++i) qs[i] = c[(uint8_t)qseq[i]];
  if (approximate_max__drop) flag |= KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP;
  // Align
  timer_start(&align_input->timer);
  ksw_extd2_sse(
      km,qlen,(uint8_t*)qseq,tlen,(uint8_t*)tseq,5,
      mat,gapo1,gape1,gapo2,gape2,band_width,zdrop,0,flag,&ez);
  timer_stop(&align_input->timer);
  // Adapt CIGAR
  cigar_t cigar;
  benchmark_ksw2_adapt_cigar(align_input,&ez,&cigar);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,gap_affine_2p,false,&cigar);
  }
  // Free
  free(cigar.operations);
  free(ez.cigar);
  free(ts);
  free(qs);
}

