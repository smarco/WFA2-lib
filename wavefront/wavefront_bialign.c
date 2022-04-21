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
  fprintf(stderr,"[WFA::BiAlign] [Recursion=%d] ",rlevel);
  int i; for (i=0;i<rlevel;++i) fprintf(stderr,"   ");
  fprintf(stderr,"Breakpoint at (h,v,score) = (%d,%d,%d)\n",
      breakpoint_h,breakpoint_v,breakpoint->score);
}
/*
 * Bialign check breakpoints
 */
void wavefront_bialign_breakpoint_i2i(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wavefront_t* const iwf_forward,
    wavefront_t* const iwf_reverse,
    const int k_forward,
    const int k_reverse,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const int text_length = wf_aligner_forward->text_length;
  const int gap_opening = wf_aligner_forward->penalties.gap_opening1;
  // Check breakpoint i2i
  const wf_offset_t ioffset_forward = iwf_forward->offsets[k_forward];
  const wf_offset_t ioffset_reverse = iwf_reverse->offsets[k_reverse];
  const int ih_forward = WAVEFRONT_H(k_forward,ioffset_forward);
  const int ih_reverse = WAVEFRONT_H(k_reverse,ioffset_reverse);
  if (ih_forward + ih_reverse >= text_length &&
      score_forward + score_reverse - gap_opening < breakpoint->score) {
    breakpoint->score_forward = score_forward;
    breakpoint->score_reverse = score_reverse;
    breakpoint->score = score_forward + score_reverse - gap_opening;
    breakpoint->k_forward = k_forward;
    breakpoint->k_reverse = k_reverse;
    breakpoint->offset_forward = ih_forward;
    breakpoint->offset_reverse = ih_reverse;
    breakpoint->component = affine_matrix_I;
  }
}
void wavefront_bialign_breakpoint_d2d(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wavefront_t* const dwf_forward,
    wavefront_t* const dwf_reverse,
    const int k_forward,
    const int k_reverse,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const int text_length = wf_aligner_forward->text_length;
  const int gap_opening = wf_aligner_forward->penalties.gap_opening1;
  // Check breakpoint d2d
  const wf_offset_t doffset_forward = dwf_forward->offsets[k_forward];
  const wf_offset_t doffset_reverse = dwf_reverse->offsets[k_reverse];
  const int dh_forward = WAVEFRONT_H(k_forward,doffset_forward);
  const int dh_reverse = WAVEFRONT_H(k_reverse,doffset_reverse);
  if (dh_forward + dh_reverse >= text_length &&
      score_forward + score_reverse - gap_opening < breakpoint->score) {
    breakpoint->score_forward = score_forward;
    breakpoint->score_reverse = score_reverse;
    breakpoint->score = score_forward + score_reverse - gap_opening;
    breakpoint->k_forward = k_forward;
    breakpoint->k_reverse = k_reverse;
    breakpoint->offset_forward = dh_forward;
    breakpoint->offset_reverse = dh_reverse;
    breakpoint->component = affine_matrix_D;
  }
}
void wavefront_bialign_breakpoint_m2m(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wavefront_t* const mwf_forward,
    wavefront_t* const mwf_reverse,
    const int k_forward,
    const int k_reverse,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const int text_length = wf_aligner_forward->text_length;
  // Check breakpoint m2m
  const wf_offset_t moffset_forward = mwf_forward->offsets[k_forward];
  const wf_offset_t moffset_reverse = mwf_reverse->offsets[k_reverse];
  const int mh_forward = WAVEFRONT_H(k_forward,moffset_forward);
  const int mh_reverse = WAVEFRONT_H(k_reverse,moffset_reverse);
  if (mh_forward + mh_reverse >= text_length &&
      score_forward + score_reverse < breakpoint->score) {
    breakpoint->score_forward = score_forward;
    breakpoint->score_reverse = score_reverse;
    breakpoint->score = score_forward + score_reverse;
    breakpoint->k_forward = k_forward;
    breakpoint->k_reverse = k_reverse;
    breakpoint->offset_forward = moffset_forward;
    breakpoint->offset_reverse = moffset_reverse;
    breakpoint->component = affine_matrix_M;
  }
}
/*
 * Bialign find overlaps
 */
void wavefront_bialign_forward_overlap_diagonal(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    const int k_forward,
    wf_bialign_breakpoint_t* const best_breakpoint) {
  // Parameters
  const int pattern_length = wf_aligner_forward->pattern_length;
  const int text_length = wf_aligner_forward->text_length;
  const bool memory_modular = wf_aligner_forward->wf_components.memory_modular;
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  // Fetch wavefronts forward
  const int score_forward_mod = (memory_modular) ? (score_forward % max_score_scope) : score_forward;
  wavefront_t* const mwf_forward = wf_aligner_forward->wf_components.mwavefronts[score_forward_mod];
  wavefront_t* const dwf_forward = wf_aligner_forward->wf_components.d1wavefronts[score_forward_mod];
  wavefront_t* const iwf_forward = wf_aligner_forward->wf_components.i1wavefronts[score_forward_mod];
  // Traverse all reverse scores
  const int k_reverse = WAVEFRONT_K_INVERSE(k_forward,pattern_length,text_length);
  int i;
  for (i=0;i<max_score_scope;++i) {
    // Compute score
    const int score_rev = score_reverse - i;
    if (score_rev < 0) break;
    const int score_reverse_mod = (memory_modular) ? (score_rev % max_score_scope) : score_rev;
    // Fetch wavefronts reverse
    wavefront_t* const mwf_reverse = wf_aligner_reverse->wf_components.mwavefronts[score_reverse_mod];
    wavefront_t* const dwf_reverse = wf_aligner_reverse->wf_components.d1wavefronts[score_reverse_mod];
    wavefront_t* const iwf_reverse = wf_aligner_reverse->wf_components.i1wavefronts[score_reverse_mod];
    // Check breakpoint m2m
    if (mwf_reverse != NULL && mwf_reverse->lo <= k_reverse && k_reverse <= mwf_reverse->hi) {
      wavefront_bialign_breakpoint_m2m(
          wf_aligner_forward,wf_aligner_reverse,
          score_forward,score_rev,
          mwf_forward,mwf_reverse,
          k_forward,k_reverse,best_breakpoint);
    }
    // Check breakpoint d2d
    if (dwf_forward != NULL && dwf_reverse != NULL && dwf_reverse->lo <= k_reverse && k_reverse <= dwf_reverse->hi) {
      wavefront_bialign_breakpoint_d2d(
          wf_aligner_forward,wf_aligner_reverse,
          score_forward,score_rev,
          dwf_forward,dwf_reverse,
          k_forward,k_reverse,best_breakpoint);
    }
    // Check breakpoint i2i
    if (iwf_forward != NULL && iwf_reverse != NULL && iwf_reverse->lo <= k_reverse && k_reverse <= iwf_reverse->hi) {
      wavefront_bialign_breakpoint_i2i(
          wf_aligner_forward,wf_aligner_reverse,
          score_forward,score_rev,
          iwf_forward,iwf_reverse,
          k_forward,k_reverse,best_breakpoint);
    }
  }
}
void wavefront_bialign_reverse_overlap_diagonal(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    const int k_reverse,
    wf_bialign_breakpoint_t* const best_breakpoint) {
  // Parameters
  const int pattern_length = wf_aligner_forward->pattern_length;
  const int text_length = wf_aligner_forward->text_length;
  const bool memory_modular = wf_aligner_forward->wf_components.memory_modular;
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  // Fetch wavefronts reverse
  const int score_reverse_mod = (memory_modular) ? (score_reverse % max_score_scope) : score_reverse;
  wavefront_t* const mwf_reverse = wf_aligner_reverse->wf_components.mwavefronts[score_reverse_mod];
  wavefront_t* const dwf_reverse = wf_aligner_reverse->wf_components.d1wavefronts[score_reverse_mod];
  wavefront_t* const iwf_reverse = wf_aligner_reverse->wf_components.i1wavefronts[score_reverse_mod];
  // Traverse all reverse scores
  const int k_forward = WAVEFRONT_K_INVERSE(k_reverse,pattern_length,text_length);
  int i;
  for (i=0;i<max_score_scope;++i) {
    // Compute score
    const int score_for = score_forward - i;
    if (score_for < 0) break;
    const int score_for_mod = (memory_modular) ? (score_for % max_score_scope) : score_for;
    // Fetch wavefronts forward
    wavefront_t* const mwf_forward = wf_aligner_forward->wf_components.mwavefronts[score_for_mod];
    wavefront_t* const dwf_forward = wf_aligner_forward->wf_components.d1wavefronts[score_for_mod];
    wavefront_t* const iwf_forward = wf_aligner_forward->wf_components.i1wavefronts[score_for_mod];
    // Check breakpoint m2m
    if (mwf_forward != NULL && mwf_forward->lo <= k_forward && k_forward <= mwf_forward->hi) {
      wavefront_bialign_breakpoint_m2m(
          wf_aligner_forward,wf_aligner_reverse,
          score_for,score_reverse,
          mwf_forward,mwf_reverse,
          k_forward,k_reverse,best_breakpoint);
    }
    // Check breakpoint d2d
    if (dwf_reverse != NULL && dwf_forward != NULL && dwf_forward->lo <= k_forward && k_forward <= dwf_forward->hi) {
      wavefront_bialign_breakpoint_d2d(
          wf_aligner_forward,wf_aligner_reverse,
          score_for,score_reverse,
          dwf_forward,dwf_reverse,
          k_forward,k_reverse,best_breakpoint);
    }
    // Check breakpoint i2i
    if (iwf_reverse != NULL && iwf_forward != NULL && iwf_forward->lo <= k_forward && k_forward <= iwf_forward->hi) {
      wavefront_bialign_breakpoint_i2i(
          wf_aligner_forward,wf_aligner_reverse,
          score_for,score_reverse,
          iwf_forward,iwf_reverse,
          k_forward,k_reverse,best_breakpoint);
    }
  }
}
void wavefront_bialign_forward_overlap(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const bool memory_modular = wf_aligner_forward->wf_components.memory_modular;
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  // DEBUG
  //wavefront_aligner_print(stderr,wf_aligner_forward,score_forward,score_forward,6,0); // DEBUG
  //wavefront_aligner_print(stderr,wf_aligner_reverse,score_reverse,score_reverse,6,0); // DEBUG
  // Fetch m-wavefronts forward
  const int score_forward_mod = (memory_modular) ? (score_forward % max_score_scope) : score_forward;
  wavefront_t* const mwf_forward = wf_aligner_forward->wf_components.mwavefronts[score_forward_mod];
  if (mwf_forward == NULL) return;
  // Traverse all diagonals and look for overlaps
  int k_forward;
  for (k_forward=mwf_forward->lo;k_forward<=mwf_forward->hi;k_forward++) {
    // Find overlaps on diagonal
    wf_bialign_breakpoint_t diagonal_breakpoint = { .score = INT_MAX };
    wavefront_bialign_forward_overlap_diagonal(
        wf_aligner_forward,wf_aligner_reverse,
        score_forward,score_reverse,k_forward,
        &diagonal_breakpoint);
    if (diagonal_breakpoint.score < breakpoint->score) {
      *breakpoint = diagonal_breakpoint;
    }
  }
}
void wavefront_bialign_reverse_overlap(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const bool memory_modular = wf_aligner_forward->wf_components.memory_modular;
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  // DEBUG
  //wavefront_aligner_print(stderr,wf_aligner_forward,score_forward,score_forward,6,0); // DEBUG
  //wavefront_aligner_print(stderr,wf_aligner_reverse,score_reverse,score_reverse,6,0); // DEBUG
  // Fetch m-wavefronts forward
  const int score_reverse_mod = (memory_modular) ? (score_reverse % max_score_scope) : score_reverse;
  wavefront_t* const mwf_reverse = wf_aligner_reverse->wf_components.mwavefronts[score_reverse_mod];
  if (mwf_reverse == NULL) return;
  // Traverse all diagonals and look for overlaps
  int k_reverse;
  for (k_reverse=mwf_reverse->lo;k_reverse<=mwf_reverse->hi;k_reverse++) {
    // Find overlaps on diagonal
    wf_bialign_breakpoint_t diagonal_breakpoint = { .score = INT_MAX };
    wavefront_bialign_reverse_overlap_diagonal(
        wf_aligner_forward,wf_aligner_reverse,
        score_forward,score_reverse,
        k_reverse,&diagonal_breakpoint);
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
    const affine_matrix_type component_begin,
    const affine_matrix_type component_end) {
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
    const affine_matrix_type component_begin,
    const affine_matrix_type component_end,
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
  const int gap_opening = wf_aligner_forward->penalties.gap_opening1;
  while (true) {
    if (last_wf_forward) {
      // Check overlapping wavefronts
      const int min_score_reverse = (score_reverse > max_score_scope-1) ? score_reverse - (max_score_scope-1) : 0;
      if (score_forward + min_score_reverse - gap_opening >= breakpoint->score) break; // Done!
      wavefront_bialign_forward_overlap(wf_aligner_forward,wf_aligner_reverse,score_forward,score_reverse,breakpoint);
      // Compute-next and extend reverse-wavefront
      ++score_reverse;
      (*wf_align_compute)(wf_aligner_reverse,score_reverse);
      wavefront_extend_end2end(wf_aligner_reverse,score_reverse);
    }
    // Check overlapping wavefronts
    const int min_score_forward = (score_forward > max_score_scope-1) ? score_forward - (max_score_scope-1) : 0;
    if (min_score_forward + score_reverse - gap_opening >= breakpoint->score) break; // Done!
    wavefront_bialign_reverse_overlap(wf_aligner_forward,wf_aligner_reverse,score_forward,score_reverse,breakpoint);
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
    const affine_matrix_type component_begin,
    const affine_matrix_type component_end,
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


