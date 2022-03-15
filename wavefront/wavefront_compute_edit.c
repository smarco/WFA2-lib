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
 * DESCRIPTION: WaveFront alignment module for computing wavefronts (edit/indel)
 */

#include "utils/string_padded.h"
#include "wavefront_compute.h"
#include "wavefront_backtrace_offload.h"

/*
 * Compute Kernels
 */
void wavefront_compute_indel_idm(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const wf_offset_t* const prev_offsets = wf_prev->offsets;
  wf_offset_t* const curr_offsets = wf_curr->offsets;
  // Loop peeling (k=lo)
  curr_offsets[lo] = prev_offsets[lo+1];
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE
  for (k=lo+1;k<=hi-1;++k) {
    // Compute maximum offset
    const wf_offset_t ins = prev_offsets[k-1] + 1;
    const wf_offset_t del = prev_offsets[k+1];
    wf_offset_t max = MAX(del,ins);
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    curr_offsets[k] = max;
  }
  // Loop peeling (k=hi)
  curr_offsets[hi] = prev_offsets[hi-1] + 1;
}
void wavefront_compute_edit_idm(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const wf_offset_t* const prev_offsets = wf_prev->offsets;
  wf_offset_t* const curr_offsets = wf_curr->offsets;
  // Loop peeling (k=lo)
  curr_offsets[lo] = prev_offsets[lo+1];
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE
  for (k=lo+1;k<=hi-1;++k) {
    // Compute maximum offset
    const wf_offset_t ins = prev_offsets[k-1]; // Lower
    const wf_offset_t del = prev_offsets[k+1]; // Upper
    const wf_offset_t misms = prev_offsets[k]; // Mid
    wf_offset_t max = MAX(del,MAX(ins,misms)+1);
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    curr_offsets[k] = max;
  }
  // Loop peeling (k=hi)
  curr_offsets[hi] = prev_offsets[hi-1] + 1;
}
/*
 * Compute Kernel (Piggyback)
 */
void wavefront_compute_indel_idm_piggyback(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi,
    const int score) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // Previous WF
  const wf_offset_t* const prev_offsets = wf_prev->offsets;
  const pcigar_t* const prev_pcigar = wf_prev->bt_pcigar;
  const bt_block_idx_t* const prev_bt_idx = wf_prev->bt_prev;
  // Current WF
  wf_offset_t* const curr_offsets = wf_curr->offsets;
  pcigar_t* const curr_pcigar = wf_curr->bt_pcigar;
  bt_block_idx_t* const curr_bt_idx = wf_curr->bt_prev;
  // Loop peeling (k=lo)
  curr_offsets[lo] = prev_offsets[lo+1];
  curr_pcigar[lo] = PCIGAR_PUSH_BACK_DEL(prev_pcigar[lo+1]);
  curr_bt_idx[lo] = prev_bt_idx[lo+1];
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE // Ifs predicated by the compiler
  for (k=lo+1;k<=hi-1;++k) {
    // Compute maximum offset
    const wf_offset_t ins = prev_offsets[k-1] + 1;
    const wf_offset_t del = prev_offsets[k+1];
    wf_offset_t max = MAX(del,ins);
    // Update pcigar & bt-block
    if (max == del) {
      curr_pcigar[k] = PCIGAR_PUSH_BACK_DEL(prev_pcigar[k+1]);
      curr_bt_idx[k] = prev_bt_idx[k+1];
    } else { // max == ins
      curr_pcigar[k] = PCIGAR_PUSH_BACK_INS(prev_pcigar[k-1]);
      curr_bt_idx[k] = prev_bt_idx[k-1];
    }
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    curr_offsets[k] = max;
  }
  // Loop peeling (k=hi)
  curr_offsets[hi] = prev_offsets[hi-1] + 1;
  curr_pcigar[hi] = PCIGAR_PUSH_BACK_INS(prev_pcigar[hi-1]);
  curr_bt_idx[hi] = prev_bt_idx[hi-1];
  // Offload backtrace (if necessary)
  if (score % PCIGAR_MAX_LENGTH == 0) {
    wavefront_backtrace_offload_blocks_linear(
        wf_aligner,curr_offsets,curr_pcigar,curr_bt_idx,lo,hi);
  }
}
void wavefront_compute_edit_idm_piggyback(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi,
    const int score) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // Previous WF
  const wf_offset_t* const prev_offsets = wf_prev->offsets;
  const pcigar_t* const prev_pcigar = wf_prev->bt_pcigar;
  const bt_block_idx_t* const prev_bt_idx = wf_prev->bt_prev;
  // Current WF
  wf_offset_t* const curr_offsets = wf_curr->offsets;
  pcigar_t* const curr_pcigar = wf_curr->bt_pcigar;
  bt_block_idx_t* const curr_bt_idx = wf_curr->bt_prev;
  // Loop peeling (k=lo)
  curr_offsets[lo] = prev_offsets[lo+1];
  curr_pcigar[lo] = PCIGAR_PUSH_BACK_DEL(prev_pcigar[lo+1]);
  curr_bt_idx[lo] = prev_bt_idx[lo+1];
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE // Ifs predicated by the compiler
  for (k=lo+1;k<=hi-1;++k) {
    // Compute maximum offset
    const wf_offset_t ins = prev_offsets[k-1] + 1; // Lower
    const wf_offset_t del = prev_offsets[k+1];     // Upper
    const wf_offset_t misms = prev_offsets[k] + 1; // Mid
    wf_offset_t max = MAX(del,MAX(ins,misms));
    // Update pcigar & bt-block
    if (max == ins) {
      curr_pcigar[k] = PCIGAR_PUSH_BACK_INS(prev_pcigar[k-1]);
      curr_bt_idx[k] = prev_bt_idx[k-1];
    }
    if (max == del) {
      curr_pcigar[k] = PCIGAR_PUSH_BACK_DEL(prev_pcigar[k+1]);
      curr_bt_idx[k] = prev_bt_idx[k+1];
    }
    if (max == misms) {
      curr_pcigar[k] = PCIGAR_PUSH_BACK_MISMS(prev_pcigar[k]);
      curr_bt_idx[k] = prev_bt_idx[k];
    }
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    curr_offsets[k] = max;
  }
  // Loop peeling (k=hi)
  curr_offsets[hi] = prev_offsets[hi-1] + 1;
  curr_pcigar[hi] = PCIGAR_PUSH_BACK_INS(prev_pcigar[hi-1]);
  curr_bt_idx[hi] = prev_bt_idx[hi-1];
  // Offload backtrace (if necessary)
  if (score % PCIGAR_MAX_LENGTH == 0) {
    wavefront_backtrace_offload_blocks_linear(
        wf_aligner,curr_offsets,curr_pcigar,curr_bt_idx,lo,hi);
  }
}
/*
 * Compute next wavefront
 */
void wavefront_compute_edit(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  // Compute scores
  int score_prev = score - 1;
  int score_curr = score;
  if (wf_components->memory_modular) { // Modular wavefront
    score_prev = score_prev % wf_components->max_score_scope;
    score_curr = score_curr % wf_components->max_score_scope;
    if (wf_components->mwavefronts[score_curr]) { // Free
      wavefront_slab_free(wf_aligner->wavefront_slab,wf_components->mwavefronts[score_curr]);
    }
  }
  // Fetch previous wavefront, compute limits & initialize
  wavefront_t* const wf_prev = wf_components->mwavefronts[score_prev];
  const int lo = wf_prev->lo - 1;
  const int hi = wf_prev->hi + 1;
  //  wf_components->historic_min_lo = min_lo;
  //  wf_components->historic_max_hi = max_hi;
  wf_prev->offsets[lo] = WAVEFRONT_OFFSET_NULL;
  wf_prev->offsets[hi] = WAVEFRONT_OFFSET_NULL;
  // Allocate output wavefront
  wavefront_t* const wf_curr = wavefront_slab_allocate(wf_aligner->wavefront_slab,lo-1,hi+1);
  wf_components->mwavefronts[score_curr] = wf_curr;
  wf_components->mwavefronts[score_curr]->lo = lo;
  wf_components->mwavefronts[score_curr]->hi = hi;
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) { 
    // For heuristics computations (as a buffer)
    wavefront_components_resize_null__victim(wf_components,lo-1,hi+1);
  }
  // Compute next wavefront
  if (wf_aligner->wf_components.bt_piggyback) {
    if (wf_aligner->penalties.distance_metric == indel) {
      wavefront_compute_indel_idm_piggyback(wf_aligner,wf_prev,wf_curr,lo,hi,score);
    } else {
      wavefront_compute_edit_idm_piggyback(wf_aligner,wf_prev,wf_curr,lo,hi,score);
    }
  } else {
    if (wf_aligner->penalties.distance_metric == indel) {
      wavefront_compute_indel_idm(wf_aligner,wf_prev,wf_curr,lo,hi);
    } else {
      wavefront_compute_edit_idm(wf_aligner,wf_prev,wf_curr,lo,hi);
    }
  }
  // Trim wavefront ends
  wavefront_compute_trim_ends(wf_aligner,wf_curr);
}


