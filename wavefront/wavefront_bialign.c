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
 */

#include "wavefront_bialign.h"

#include "wavefront_align.h"
#include "wavefront_extend.h"
#include "wavefront_compute.h"
#include "wavefront_compute_edit.h"
#include "wavefront_compute_linear.h"
#include "wavefront_compute_affine.h"
#include "wavefront_compute_affine2p.h"
#include "wavefront_backtrace.h"

/*
 * Config
 */
#define WF_BIALIGN_FALLBACK_MIN_SCORE 100

/*
 * Debug
 */
void wavefront_bialign_debug(
    wf_bialign_breakpoint_t* const breakpoint,
    const int rlevel) {
  // Parameters
  const int breakpoint_h = WAVEFRONT_H(breakpoint->k_forward,breakpoint->offset_forward);
  const int breakpoint_v = WAVEFRONT_V(breakpoint->k_forward,breakpoint->offset_forward);
  // Prinf debug info
  fprintf(stderr,"[WFA::BiAlign][Recursion=%d] ",rlevel);
  int i; for (i=0;i<rlevel;++i) fprintf(stderr,"   ");
  fprintf(stderr,"Breakpoint at (h,v,score,comp) = (%d,%d,%d,",
      breakpoint_h,breakpoint_v,breakpoint->score);
  switch (breakpoint->component) {
    case affine2p_matrix_M:  fprintf(stderr,"M");  break;
    case affine2p_matrix_I1: fprintf(stderr,"I1"); break;
    case affine2p_matrix_I2: fprintf(stderr,"I2"); break;
    case affine2p_matrix_D1: fprintf(stderr,"D1"); break;
    case affine2p_matrix_D2: fprintf(stderr,"D2"); break;
    default: fprintf(stderr,"?"); break;
  }
  fprintf(stderr,")\n");
}
/*
 * Utils
 */
int wavefront_bialign_gap_opening_adjustment(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric) {
  switch (distance_metric) {
    case gap_affine:
      return wf_aligner->penalties.gap_opening1;
    case gap_affine_2p:
      return MAX(wf_aligner->penalties.gap_opening1,wf_aligner->penalties.gap_opening2);
    case indel:
    case edit:
    case gap_linear:
    default:
      return 0;
  }
}
/*
 * Bialign check breakpoints
 */
void wavefront_bialign_breakpoint_i2i(
    wavefront_aligner_t* const wf_aligner,
    const bool breakpoint_forward,
    const int score_0,
    const int score_1,
    wavefront_t* const iwf_0,
    wavefront_t* const iwf_1,
    const int k_0,
    const int k_1,
    const bool gap_affine_2p,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const int text_length = wf_aligner->text_length;
  const int gap_open = (!gap_affine_2p) ?
      wf_aligner->penalties.gap_opening1 :
      wf_aligner->penalties.gap_opening2;
  // Check breakpoint i2i
  const wf_offset_t ioffset_0 = iwf_0->offsets[k_0];
  const wf_offset_t ioffset_1 = iwf_1->offsets[k_1];
  const int ih_0 = WAVEFRONT_H(k_0,ioffset_0);
  const int ih_1 = WAVEFRONT_H(k_1,ioffset_1);
  if (ih_0 + ih_1 >= text_length && score_0 + score_1 - gap_open < breakpoint->score) {
    if (breakpoint_forward) {
      breakpoint->score_forward = score_0;
      breakpoint->score_reverse = score_1;
      breakpoint->k_forward = k_0;
      breakpoint->k_reverse = k_1;
      breakpoint->offset_forward = ih_0;
      breakpoint->offset_reverse = ih_1;
    } else {
      breakpoint->score_forward = score_1;
      breakpoint->score_reverse = score_0;
      breakpoint->k_forward = k_1;
      breakpoint->k_reverse = k_0;
      breakpoint->offset_forward = ih_1;
      breakpoint->offset_reverse = ih_0;
    }
    breakpoint->score = score_0 + score_1 - gap_open;
    breakpoint->component = (!gap_affine_2p) ? affine2p_matrix_I1 : affine2p_matrix_I2;
  }
}
void wavefront_bialign_breakpoint_d2d(
    wavefront_aligner_t* const wf_aligner,
    const bool breakpoint_forward,
    const int score_0,
    const int score_1,
    wavefront_t* const dwf_0,
    wavefront_t* const dwf_1,
    const int k_0,
    const int k_1,
    const bool gap_affine_2p,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const int text_length = wf_aligner->text_length;
  const int gap_open = (!gap_affine_2p) ?
      wf_aligner->penalties.gap_opening1 :
      wf_aligner->penalties.gap_opening2;
  // Check breakpoint d2d
  const wf_offset_t doffset_0 = dwf_0->offsets[k_0];
  const wf_offset_t doffset_1 = dwf_1->offsets[k_1];
  const int dh_0 = WAVEFRONT_H(k_forward,doffset_0);
  const int dh_1 = WAVEFRONT_H(k_reverse,doffset_1);
  if (dh_0 + dh_1 >= text_length && score_0 + score_1 - gap_open < breakpoint->score) {
    if (breakpoint_forward) {
      breakpoint->score_forward = score_0;
      breakpoint->score_reverse = score_1;
      breakpoint->k_forward = k_0;
      breakpoint->k_reverse = k_1;
      breakpoint->offset_forward = dh_0;
      breakpoint->offset_reverse = dh_1;
    } else {
      breakpoint->score_forward = score_1;
      breakpoint->score_reverse = score_0;
      breakpoint->k_forward = k_1;
      breakpoint->k_reverse = k_0;
      breakpoint->offset_forward = dh_1;
      breakpoint->offset_reverse = dh_0;
    }
    breakpoint->score = score_0 + score_1 - gap_open;
    breakpoint->component = (!gap_affine_2p) ? affine2p_matrix_D1 : affine2p_matrix_D2;
  }
}
void wavefront_bialign_breakpoint_m2m(
    wavefront_aligner_t* const wf_aligner,
    const bool breakpoint_forward,
    const int score_0,
    const int score_1,
    wavefront_t* const mwf_0,
    wavefront_t* const mwf_1,
    const int k_0,
    const int k_1,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const int text_length = wf_aligner->text_length;
  // Check breakpoint m2m
  const wf_offset_t moffset_0 = mwf_0->offsets[k_0];
  const wf_offset_t moffset_1 = mwf_1->offsets[k_1];
  const int mh_0 = WAVEFRONT_H(k_forward,moffset_0);
  const int mh_1 = WAVEFRONT_H(k_reverse,moffset_1);
  if (mh_0 + mh_1 >= text_length && score_0 + score_1 < breakpoint->score) {
    if (breakpoint_forward) {
      breakpoint->score_forward = score_0;
      breakpoint->score_reverse = score_1;
      breakpoint->k_forward = k_0;
      breakpoint->k_reverse = k_1;
      breakpoint->offset_forward = moffset_0;
      breakpoint->offset_reverse = moffset_1;
    } else {
      breakpoint->score_forward = score_1;
      breakpoint->score_reverse = score_0;
      breakpoint->k_forward = k_1;
      breakpoint->k_reverse = k_0;
      breakpoint->offset_forward = moffset_1;
      breakpoint->offset_reverse = moffset_0;
    }
    breakpoint->score = score_0 + score_1;
    breakpoint->component = affine2p_matrix_M;
  }
}
/*
 * Bialign find overlaps
 */
void wavefront_bialign_overlap_diagonal(
    wavefront_aligner_t* const wf_aligner_0,
    wavefront_aligner_t* const wf_aligner_1,
    const int score_0,
    const int score_1,
    const int k_0,
    const bool breakpoint_forward,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const int pattern_length = wf_aligner_0->pattern_length;
  const int text_length = wf_aligner_0->text_length;
  const int max_score_scope = wf_aligner_0->wf_components.max_score_scope;
  const distance_metric_t distance_metric = wf_aligner_0->penalties.distance_metric;
  // Fetch wavefronts-0
  const int score_mod_0 = score_0 % max_score_scope;
  wavefront_t* const mwf_0 = wf_aligner_0->wf_components.mwavefronts[score_mod_0];
  wavefront_t* d1wf_0 = NULL, *i1wf_0 = NULL;
  if (distance_metric >= gap_affine) {
    d1wf_0 = wf_aligner_0->wf_components.d1wavefronts[score_mod_0];
    i1wf_0 = wf_aligner_0->wf_components.i1wavefronts[score_mod_0];
  }
  wavefront_t* d2wf_0 = NULL, *i2wf_0 = NULL;
  if (distance_metric == gap_affine_2p) {
    d2wf_0 = wf_aligner_0->wf_components.d2wavefronts[score_mod_0];
    i2wf_0 = wf_aligner_0->wf_components.i2wavefronts[score_mod_0];
  }
  // Traverse all scores-1
  const int k_1 = WAVEFRONT_K_INVERSE(k_0,pattern_length,text_length);
  int i;
  for (i=0;i<max_score_scope;++i) {
    // Compute score
    const int score_i = score_1 - i;
    if (score_i < 0) break;
    const int score_mod_i = score_i % max_score_scope;
    // Check breakpoint m2m
    wavefront_t* const mwf_1 = wf_aligner_1->wf_components.mwavefronts[score_mod_i];
    if (mwf_1 != NULL && mwf_1->lo <= k_1 && k_1 <= mwf_1->hi) {
      wavefront_bialign_breakpoint_m2m(wf_aligner_0,breakpoint_forward,
          score_0,score_i,mwf_0,mwf_1,k_0,k_1,breakpoint);
    }
    if (distance_metric <= gap_linear) continue;
    // Check breakpoint d2d
    wavefront_t* const d1wf_1 = wf_aligner_1->wf_components.d1wavefronts[score_mod_i];
    if (d1wf_0 != NULL && d1wf_1 != NULL && d1wf_1->lo <= k_1 && k_1 <= d1wf_1->hi) {
      wavefront_bialign_breakpoint_d2d(wf_aligner_0,breakpoint_forward,
          score_0,score_i,d1wf_0,d1wf_1,k_0,k_1,false,breakpoint);
    }
    // Check breakpoint i2i
    wavefront_t* const i1wf_1 = wf_aligner_1->wf_components.i1wavefronts[score_mod_i];
    if (i1wf_0 != NULL && i1wf_1 != NULL && i1wf_1->lo <= k_1 && k_1 <= i1wf_1->hi) {
      wavefront_bialign_breakpoint_i2i(wf_aligner_0,breakpoint_forward,
          score_0,score_i,i1wf_0,i1wf_1,k_0,k_1,false,breakpoint);
    }
    if (distance_metric == gap_affine) continue;
    // Check breakpoint d2d
    wavefront_t* const d2wf_1 = wf_aligner_1->wf_components.d2wavefronts[score_mod_i];
    if (d2wf_0 != NULL && d2wf_1 != NULL && d2wf_1->lo <= k_1 && k_1 <= d2wf_1->hi) {
      wavefront_bialign_breakpoint_d2d(wf_aligner_0,breakpoint_forward,
          score_0,score_i,d2wf_0,d2wf_1,k_0,k_1,true,breakpoint);
    }
    // Check breakpoint i2i
    wavefront_t* const i2wf_1 = wf_aligner_1->wf_components.i2wavefronts[score_mod_i];
    if (i2wf_0 != NULL && i2wf_1 != NULL && i2wf_1->lo <= k_1 && k_1 <= i2wf_1->hi) {
      wavefront_bialign_breakpoint_i2i(wf_aligner_0,breakpoint_forward,
          score_0,score_i,i2wf_0,i2wf_1,k_0,k_1,true,breakpoint);
    }
  }
}
void wavefront_bialign_overlap(
    wavefront_aligner_t* const wf_aligner_0,
    wavefront_aligner_t* const wf_aligner_1,
    const int score_0,
    const int score_1,
    const bool breakpoint_forward,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const int max_score_scope = wf_aligner_0->wf_components.max_score_scope;
  // Fetch m-wavefronts-0
  const int score_mod_0 = score_0 % max_score_scope;
  wavefront_t* const mwf_0 = wf_aligner_0->wf_components.mwavefronts[score_mod_0];
  if (mwf_0 == NULL) return;
  // Traverse all diagonals and look for overlaps
  int k_0;
  for (k_0=mwf_0->lo;k_0<=mwf_0->hi;k_0++) {
    // Find overlaps on diagonal
    wf_bialign_breakpoint_t diagonal_breakpoint = { .score = INT_MAX };
    wavefront_bialign_overlap_diagonal(
        wf_aligner_0,wf_aligner_1,score_0,score_1,k_0,
        breakpoint_forward,&diagonal_breakpoint);
    if (diagonal_breakpoint.score < breakpoint->score) {
      *breakpoint = diagonal_breakpoint;
    }
  }
}
/*
 * Bialign breakpoint detection
 */
void wavefront_bialign_find_breakpoint_init(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const distance_metric_t distance_metric,
    alignment_form_t* const form,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end) {
  // Resize wavefront aligner
  wavefront_aligner_resize(wf_aligner_forward,pattern,pattern_length,text,text_length,false);
  wavefront_aligner_resize(wf_aligner_reverse,pattern,pattern_length,text,text_length,true);
  // Configure form forward and reverse
  alignment_span_t span_forward =
      (form->pattern_begin_free > 0 || form->text_begin_free > 0) ? alignment_endsfree : alignment_end2end;
  alignment_form_t form_forward = {
      .span = span_forward,
      .pattern_begin_free = form->pattern_begin_free,
      .pattern_end_free = 0,
      .text_begin_free = form->text_begin_free,
      .text_end_free = 0,
  };
  alignment_span_t span_reverse =
      (form->pattern_end_free > 0 || form->text_end_free > 0) ? alignment_endsfree : alignment_end2end;
  alignment_form_t form_reverse = {
      .span = span_reverse,
      .pattern_begin_free = form->pattern_end_free,
      .pattern_end_free = 0,
      .text_begin_free = form->text_end_free,
      .text_end_free = 0,
  };
  // Configure WF-compute function (global)
  switch (distance_metric) {
    case indel:
    case edit:
      wf_aligner_forward->align_status.wf_align_compute = &wavefront_compute_edit;
      break;
    case gap_linear:
      wf_aligner_forward->align_status.wf_align_compute = &wavefront_compute_linear;
      break;
    case gap_affine:
      wf_aligner_forward->align_status.wf_align_compute = &wavefront_compute_affine;
      break;
    case gap_affine_2p:
      wf_aligner_forward->align_status.wf_align_compute = &wavefront_compute_affine2p;
      break;
    default:
      fprintf(stderr,"[WFA] Distance function not implemented\n");
      exit(1);
      break;
  }
  // Initialize wavefront (forward)
  wf_aligner_forward->alignment_form = form_forward;
  wf_aligner_forward->component_begin = component_begin;
  if (span_forward == alignment_end2end) {
    wavefront_align_end2end_initialize(wf_aligner_forward);
  } else {
    wavefront_align_endsfree_initialize(wf_aligner_forward,pattern_length,text_length);
  }
  // Initialize wavefront (reverse)
  wf_aligner_reverse->alignment_form = form_reverse;
  wf_aligner_reverse->component_begin = component_end;
  if (span_reverse == alignment_end2end) {
    wavefront_align_end2end_initialize(wf_aligner_reverse);
  } else {
    wavefront_align_endsfree_initialize(wf_aligner_reverse,pattern_length,text_length);
  }
}
void wavefront_bialign_find_breakpoint(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const distance_metric_t distance_metric,
    alignment_form_t* const form,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Init bialignment
  wavefront_bialign_find_breakpoint_init(
      wf_aligner_forward,wf_aligner_reverse,
      pattern,pattern_length,text,text_length,
      distance_metric,form,component_begin,component_end);
  // Compute wavefronts of increasing score until both wavefronts overlap
  const int max_antidiagonal = DPMATRIX_ANTIDIAGONAL(pattern_length,text_length) - 1;
  void (*wf_align_compute)(wavefront_aligner_t* const,const int) = wf_aligner_forward->align_status.wf_align_compute;
  breakpoint->score = INT_MAX;
  int score_forward = 0, score_reverse = 0;
  int forward_max_ak = wavefront_extend_end2end_max(wf_aligner_forward,score_forward);
  int reverse_max_ak = wavefront_extend_end2end_max(wf_aligner_reverse,score_reverse);
  int max_ak;
  bool last_wf_forward;
  while (true) {
    // Check if they are close to collision
    if (forward_max_ak + reverse_max_ak >= max_antidiagonal) break;
    // Compute-next & extend wavefront forward
    ++score_forward;
    (*wf_align_compute)(wf_aligner_forward,score_forward);
    max_ak = wavefront_extend_end2end_max(wf_aligner_forward,score_forward);
    if (forward_max_ak < max_ak) forward_max_ak = max_ak;
    last_wf_forward = true;
    // Check if they are close to collision
    if (forward_max_ak + reverse_max_ak >= max_antidiagonal) break;
    // Compute-next & extend wavefront reverse
    ++score_reverse;
    (*wf_align_compute)(wf_aligner_reverse,score_reverse);
    max_ak = wavefront_extend_end2end_max(wf_aligner_reverse,score_reverse);
    if (reverse_max_ak < max_ak) reverse_max_ak = max_ak;
    last_wf_forward = false;
  }
  // Advance until overlap is found
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  const int gap_opening = wavefront_bialign_gap_opening_adjustment(wf_aligner_forward,distance_metric);
  while (true) {
    if (last_wf_forward) {
      // Check overlapping wavefronts
      const int min_score_reverse = (score_reverse > max_score_scope-1) ? score_reverse - (max_score_scope-1) : 0;
      if (score_forward + min_score_reverse - gap_opening >= breakpoint->score) break; // Done!
      wavefront_bialign_overlap(wf_aligner_forward,wf_aligner_reverse,score_forward,score_reverse,true,breakpoint);
      // Compute-next and extend reverse-wavefront
      ++score_reverse;
      (*wf_align_compute)(wf_aligner_reverse,score_reverse);
      wavefront_extend_end2end(wf_aligner_reverse,score_reverse);
    }
    // Check overlapping wavefronts
    const int min_score_forward = (score_forward > max_score_scope-1) ? score_forward - (max_score_scope-1) : 0;
    if (min_score_forward + score_reverse - gap_opening >= breakpoint->score) break; // Done!
    wavefront_bialign_overlap(wf_aligner_reverse,wf_aligner_forward,score_reverse,score_forward,false,breakpoint);
    // Compute-next and extend forward-wavefront
    ++score_forward;
    (*wf_align_compute)(wf_aligner_forward,score_forward);
    wavefront_extend_end2end(wf_aligner_forward,score_forward);
    // Enable always
    last_wf_forward = true;
  }
}
/*
 * Bidirectional Alignment
 */
int wavefront_align_unidirectional(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);
void wavefront_bialign_init_half_0(
    alignment_form_t* const global_form,
    alignment_form_t* const half_form) {
  // Align half_0
  const alignment_span_t span_0 =
      (global_form->pattern_begin_free > 0 ||
       global_form->text_begin_free > 0) ?
           alignment_endsfree : alignment_end2end;
  half_form->span = span_0;
  half_form->pattern_begin_free = global_form->pattern_begin_free;
  half_form->pattern_end_free = 0;
  half_form->text_begin_free = global_form->text_begin_free;
  half_form->text_end_free = 0;
}
void wavefront_bialign_init_half_1(
    alignment_form_t* const global_form,
    alignment_form_t* const half_form) {
  // Align half_0
  const alignment_span_t span_1 =
      (global_form->pattern_begin_free > 0 ||
       global_form->text_begin_free > 0) ?
           alignment_endsfree : alignment_end2end;
  half_form->span = span_1;
  half_form->pattern_begin_free = 0;
  half_form->pattern_end_free = global_form->pattern_end_free;
  half_form->text_begin_free = 0;
  half_form->text_end_free = global_form->text_end_free;
}
void wavefront_bialign(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    alignment_form_t* const form,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int score_remaining,
    cigar_t* const cigar,
    const int rlevel) {
  // Trivial cases
  if (text_length == 0) {
    cigar_append_deletion(cigar,pattern_length);
    return;
  } else if (pattern_length == 0) {
    cigar_append_insertion(cigar,text_length);
    return;
  }
  // Fallback to regular WFA
  if (score_remaining <= WF_BIALIGN_FALLBACK_MIN_SCORE) {
    // Align the remaining
    wf_aligner->component_begin = component_begin;
    wf_aligner->component_end = component_end;
    wf_aligner->alignment_form = *form;
    wavefront_align_unidirectional(wf_aligner,
        pattern,pattern_length,
        text,text_length);
    cigar_append(cigar,&wf_aligner->cigar);
    return;
  }
  // Find breakpoint in the alignment
  wf_bialign_breakpoint_t breakpoint;
  wavefront_bialign_find_breakpoint(
      wf_aligner->aligner_forward,wf_aligner->aligner_reverse,
      pattern,pattern_length,text,text_length,
      wf_aligner->penalties.distance_metric,
      form,component_begin,component_end,&breakpoint);
  const int breakpoint_h = WAVEFRONT_H(breakpoint.k_forward,breakpoint.offset_forward);
  const int breakpoint_v = WAVEFRONT_V(breakpoint.k_forward,breakpoint.offset_forward);
  // DEBUG
  if (wf_aligner->system.verbose == 1) wavefront_bialign_debug(&breakpoint,rlevel);
  // Align half_0
  alignment_form_t form_0;
  wavefront_bialign_init_half_0(form,&form_0);
  wavefront_bialign(
      wf_aligner,pattern,breakpoint_v,text,breakpoint_h,
      &form_0,component_begin,breakpoint.component,
      breakpoint.score_forward,cigar,rlevel+1);
  // Align half_1
  alignment_form_t form_1;
  wavefront_bialign_init_half_1(form,&form_1);
  wavefront_bialign(wf_aligner,
      pattern+breakpoint_v,pattern_length-breakpoint_v,
      text+breakpoint_h,text_length-breakpoint_h,
      &form_1,breakpoint.component,component_end,
      breakpoint.score_reverse,cigar,rlevel+1);
  // Set score
  cigar->score = -breakpoint.score;
}


