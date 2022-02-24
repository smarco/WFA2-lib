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
 * DESCRIPTION: GABA library wrapper
 */

#include "benchmark_gaba.h"
#include "benchmark/benchmark_check.h"
#include "external/libgaba/gaba_wrap.h"      /* multiple target configuration: gcc example.c libgaba.a */
#include "external/libgaba/gaba_parse.h"     /* parser: contains gaba_dp_print_cigar_forward and so on */

/*
 * Libgaba
 */
int printer(
    void* arg1,
    uint64_t len,
    char c) {
  cigar_t* const cigar = (cigar_t*)arg1;
  int i;
  for (i=0;i<len;++i) {
    cigar->operations[cigar->end_offset+i] = c;
  }
  cigar->end_offset += len;
  // Ok
  return 0;
}
void benchmark_gaba_aband(
    align_input_t* const align_input,
    affine_penalties_t* const penalties) {
  // Parameters
  const int padding = 32;
  // Initialize sequences
  uint8_t* const a = calloc(align_input->pattern_length + 2*padding,1);
  uint8_t* const b = calloc(align_input->text_length + 2*padding,1);
  uint8_t* const t = calloc(2*padding,1); // tail array
  int i;
  for (i=0;i<align_input->pattern_length;++i) {
    switch (align_input->pattern[i]) {
      case 'A': a[i] = 1; break;
      case 'C': a[i] = 2; break;
      case 'G': a[i] = 4; break;
      case 'T': a[i] = 8; break;
    }
  }
  for (;i<align_input->pattern_length+padding;++i) a[i] = 1|2|4|8;
  for (i=0;i<align_input->text_length;++i) {
    switch (align_input->text[i]) {
    case 'A': b[i] = 1; break;
    case 'C': b[i] = 2; break;
    case 'G': b[i] = 4; break;
    case 'T': b[i] = 8; break;
    }
  }
  for (;i<align_input->text_length+padding;++i) b[i] = 1|2|4|8;
  // Initialize a global context
  int match_penalty = (penalties->match==0) ? 1 : -penalties->match;
  gaba_t* ctx = gaba_init(
      GABA_PARAMS(.xdrop = 100,
          GABA_SCORE_SIMPLE(match_penalty,penalties->mismatch,
                            penalties->gap_opening,penalties->gap_extension)));
  struct gaba_section_s asec = gaba_build_section(0, a, align_input->pattern_length+padding);
  struct gaba_section_s bsec = gaba_build_section(2, b, align_input->text_length+padding);
  struct gaba_section_s tail = gaba_build_section(4, t+padding, padding);

  // Align
  cigar_t cigar;
  gaba_dp_t* dp;
  struct gaba_alignment_s* alignment = NULL;
  timer_start(&align_input->timer);
  // Fill root
  dp = gaba_dp_init(ctx);
  struct gaba_section_s *ap = &asec, *bp = &bsec;
  struct gaba_fill_s *f = gaba_dp_fill_root(dp, ap, 0, bp, 0, UINT32_MAX);
  struct gaba_fill_s const *m = f; // Track max
  //  // Until X-drop condition is detected
  //  while((f->status & GABA_TERM) == 0) {
  //    if(f->status & GABA_UPDATE_A) { ap = &tail; }
  //    if(f->status & GABA_UPDATE_B) { bp = &tail; }
  //    f = gaba_dp_fill(dp, f, ap, bp, UINT32_MAX);
  //    m = f->max > m->max ? f : m;
  //  }
  // Fill-in the tail of the banded matrix
  uint32_t flag = GABA_TERM;
  do {
    if (f->status & GABA_UPDATE_A) { ap = &tail; }
    if (f->status & GABA_UPDATE_B) { bp = &tail; }
    flag |= f->status & ( GABA_UPDATE_A | GABA_UPDATE_B );
    f = gaba_dp_fill(dp, f, ap, bp, UINT32_MAX);
    m = (f->max > m->max)? f : m;
  } while (!(flag & f->status));

  alignment = gaba_dp_trace(dp,m,NULL);
  // Allocate CIGAR
  cigar.operations = malloc(2*alignment->plen);
  cigar.begin_offset = 0;
  cigar.end_offset = 0;
  // Trace alignment
  gaba_print_cigar_forward(printer,(void *)(&cigar),alignment->path,0,alignment->plen);
  timer_stop(&align_input->timer);
  // Adapt CIGAR (Add mismatches)
  cigar_add_mismatches(
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,&cigar);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,gap_affine,false,&cigar);
  }
  // Free
  free(a);
  free(b);
  free(t);
  free(cigar.operations);
  gaba_dp_res_free(dp,alignment);
  gaba_dp_clean(dp);
  gaba_clean(ctx);
}

