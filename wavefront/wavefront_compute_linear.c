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
 * DESCRIPTION: WaveFront alignment module for computing wavefronts (gap-linear)
 */

#include "utils/string_padded.h"
#include "wavefront_compute.h"
#include "wavefront_backtrace_offload.h"

/*
 * Compute Kernels
 */
void wavefront_compute_linear_idm(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // In Offsets
  const wf_offset_t* const m_misms = wavefront_set->in_mwavefront_misms->offsets;
  const wf_offset_t* const m_open1 = wavefront_set->in_mwavefront_open1->offsets;
  // Out Offsets
  wf_offset_t* const out_m = wavefront_set->out_mwavefront->offsets;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE
  for (k=lo;k<=hi;++k) {
    // Compute maximum Offset
    const wf_offset_t ins1 = m_open1[k-1];
    const wf_offset_t del1 = m_open1[k+1];
    const wf_offset_t misms = m_misms[k];
    wf_offset_t max = MAX(del1,MAX(misms,ins1)+1);
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    out_m[k] = max;
  }
}
/*
 * Compute Kernel (Piggyback)
 */
void wavefront_compute_linear_idm_piggyback(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // In M
  const wf_offset_t* const m_misms = wavefront_set->in_mwavefront_misms->offsets;
  const pcigar_t* const m_misms_bt_pcigar = wavefront_set->in_mwavefront_misms->bt_pcigar;
  const bt_block_idx_t* const m_misms_bt_prev = wavefront_set->in_mwavefront_misms->bt_prev;
  // In I/D
  const wf_offset_t* const m_open1 = wavefront_set->in_mwavefront_open1->offsets;
  const pcigar_t* const m_open1_bt_pcigar = wavefront_set->in_mwavefront_open1->bt_pcigar;
  const bt_block_idx_t* const m_open1_bt_prev = wavefront_set->in_mwavefront_open1->bt_prev;
  // Out
  wf_offset_t* const out_m = wavefront_set->out_mwavefront->offsets;
  pcigar_t* const out_m_bt_pcigar = wavefront_set->out_mwavefront->bt_pcigar;
  bt_block_idx_t* const out_m_bt_prev = wavefront_set->out_mwavefront->bt_prev;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE // Ifs predicated by the compiler
  for (k=lo;k<=hi;++k) {
    // Compute maximum Offset
    const wf_offset_t ins1 = m_open1[k-1] + 1;
    const wf_offset_t del1 = m_open1[k+1];
    const wf_offset_t misms = m_misms[k] + 1;
    wf_offset_t max = MAX(del1,MAX(misms,ins1));
    // Update pcigar & bt-block
    if (max == ins1) {
      out_m_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(m_open1_bt_pcigar[k-1]);
      out_m_bt_prev[k] = m_open1_bt_prev[k-1];
    }
    if (max == del1) {
      out_m_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(m_open1_bt_pcigar[k+1]);
      out_m_bt_prev[k] = m_open1_bt_prev[k+1];
    }
    if (max == misms) {
      out_m_bt_pcigar[k] = PCIGAR_PUSH_BACK_MISMS(m_misms_bt_pcigar[k]);
      out_m_bt_prev[k] = m_misms_bt_prev[k];
    }
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    out_m[k] = max;
  }
  // Offload backtrace
  wavefront_backtrace_offload_linear(wf_aligner,wavefront_set,lo,hi);
}
/*
 * Compute next wavefront
 */
void wavefront_compute_linear(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Select wavefronts
  wavefront_set_t wavefront_set;
  wavefront_compute_fetch_input(wf_aligner,&wavefront_set,score);
  // Check null wavefronts
  if (wavefront_set.in_mwavefront_misms->null &&
      wavefront_set.in_mwavefront_open1->null) {
    wavefront_compute_allocate_output_null(wf_aligner,score); // Null s-wavefront
    return;
  }
  // Set limits
  int hi, lo;
  wavefront_compute_limits(wf_aligner,&wavefront_set,&lo,&hi);
  // Allocate wavefronts
  wavefront_compute_allocate_output(wf_aligner,&wavefront_set,score,lo,hi);
  // Init wavefront ends
  wavefront_compute_init_ends(wf_aligner,&wavefront_set,lo,hi);
  // Compute next wavefront
  if (wf_aligner->wf_components.bt_piggyback) {
    wavefront_compute_linear_idm_piggyback(wf_aligner,&wavefront_set,lo,hi);
  } else {
    wavefront_compute_linear_idm(wf_aligner,&wavefront_set,lo,hi);
  }
  // Trim wavefront ends
  wavefront_compute_trim_ends_set(wf_aligner,&wavefront_set);
}


