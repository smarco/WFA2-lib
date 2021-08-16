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
 * DESCRIPTION: WaveFront alignment module for computing wavefronts (gap-affine-2p)
 */

#include "utils/string_padded.h"
#include "wavefront_compute.h"

/*
 * Compute Kernels (Piggyback)
 */
void wavefront_compute_affine2p_idm_bounded(
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // In Offsets
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_mwavefront_sub,m_sub);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_mwavefront_gap1,m_open1);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_mwavefront_gap2,m_open2);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_i1wavefront_ext,i1_ext);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_i2wavefront_ext,i2_ext);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_d1wavefront_ext,d1_ext);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_d2wavefront_ext,d2_ext);
  // Out Offsets
  wf_offset_t* const out_m = wavefront_set->out_mwavefront->offsets;
  wf_offset_t* const out_i1 = wavefront_set->out_i1wavefront->offsets;
  wf_offset_t* const out_i2 = wavefront_set->out_i2wavefront->offsets;
  wf_offset_t* const out_d1 = wavefront_set->out_d1wavefront->offsets;
  wf_offset_t* const out_d2 = wavefront_set->out_d2wavefront->offsets;
  // Compute-Next kernel loop
  int k;
  for (k=lo;k<=hi;++k) {
    // Update I1
    const wf_offset_t ins1_o = WF_COND_FETCH(m_open1,k-1);
    const wf_offset_t ins1_e = WF_COND_FETCH(i1_ext,k-1);
    const wf_offset_t ins1 = MAX(ins1_o,ins1_e) + 1;
    out_i1[k] = ins1;
    // Update I2
    const wf_offset_t ins2_o = WF_COND_FETCH(m_open2,k-1);
    const wf_offset_t ins2_e = WF_COND_FETCH(i2_ext,k-1);
    const wf_offset_t ins2 = MAX(ins2_o,ins2_e) + 1;
    out_i2[k] = ins2;
    // Update I
    const wf_offset_t ins = MAX(ins1,ins2);
    // Update D1
    const wf_offset_t del1_o = WF_COND_FETCH(m_open1,k+1);
    const wf_offset_t del1_e = WF_COND_FETCH(d1_ext,k+1);
    const wf_offset_t del1 = MAX(del1_o,del1_e);
    out_d1[k] = del1;
    // Update D2
    const wf_offset_t del2_o = WF_COND_FETCH(m_open2,k+1);
    const wf_offset_t del2_e = WF_COND_FETCH(d2_ext,k+1);
    const wf_offset_t del2 = MAX(del2_o,del2_e);
    out_d2[k] = del2;
    // Update D
    const wf_offset_t del = MAX(del1,del2);
    // Update M
    const wf_offset_t sub = WF_COND_FETCH_INC(m_sub,k,1);
    out_m[k] = MAX(del,MAX(sub,ins));
  }
}
void wavefront_compute_affine2p_idm_unbounded(
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // In Offsets
  const wf_offset_t* const m_sub = wavefront_set->in_mwavefront_sub->offsets;
  const wf_offset_t* const m_open1 = wavefront_set->in_mwavefront_gap1->offsets;
  const wf_offset_t* const m_open2 = wavefront_set->in_mwavefront_gap2->offsets;
  const wf_offset_t* const i1_ext = wavefront_set->in_i1wavefront_ext->offsets;
  const wf_offset_t* const i2_ext = wavefront_set->in_i2wavefront_ext->offsets;
  const wf_offset_t* const d1_ext = wavefront_set->in_d1wavefront_ext->offsets;
  const wf_offset_t* const d2_ext = wavefront_set->in_d2wavefront_ext->offsets;
  // Out Offsets
  wf_offset_t* const out_m = wavefront_set->out_mwavefront->offsets;
  wf_offset_t* const out_i1 = wavefront_set->out_i1wavefront->offsets;
  wf_offset_t* const out_i2 = wavefront_set->out_i2wavefront->offsets;
  wf_offset_t* const out_d1 = wavefront_set->out_d1wavefront->offsets;
  wf_offset_t* const out_d2 = wavefront_set->out_d2wavefront->offsets;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE
  for (k=lo;k<=hi;++k) {
    // Update I1
    const wf_offset_t ins1_o = m_open1[k-1];
    const wf_offset_t ins1_e = i1_ext[k-1];
    const wf_offset_t ins1 = MAX(ins1_o,ins1_e) + 1;
    out_i1[k] = ins1;
    // Update I2
    const wf_offset_t ins2_o = m_open2[k-1];
    const wf_offset_t ins2_e = i2_ext[k-1];
    const wf_offset_t ins2 = MAX(ins2_o,ins2_e) + 1;
    out_i2[k] = ins2;
    // Update I
    const wf_offset_t ins = MAX(ins1,ins2);
    // Update D1
    const wf_offset_t del1_o = m_open1[k+1];
    const wf_offset_t del1_e = d1_ext[k+1];
    const wf_offset_t del1 = MAX(del1_o,del1_e);
    out_d1[k] = del1;
    // Update D2
    const wf_offset_t del2_o = m_open2[k+1];
    const wf_offset_t del2_e = d2_ext[k+1];
    const wf_offset_t del2 = MAX(del2_o,del2_e);
    out_d2[k] = del2;
    // Update D
    const wf_offset_t del = MAX(del1,del2);
    // Update M
    const wf_offset_t sub = m_sub[k] + 1;
    out_m[k] = MAX(del,MAX(sub,ins));
  }
}
/*
 * Compute Kernel (Piggyback)
 */
void wavefront_compute_affine2p_idm_piggyback_bounded(
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // In Offsets
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_mwavefront_sub,m_sub);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_mwavefront_gap1,m_open1);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_mwavefront_gap2,m_open2);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_i1wavefront_ext,i1_ext);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_i2wavefront_ext,i2_ext);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_d1wavefront_ext,d1_ext);
  WF_DECLARE_OFFSETS__LIMITS(wavefront_set->in_d2wavefront_ext,d2_ext);
  // Out Offsets
  wf_offset_t* const out_m = wavefront_set->out_mwavefront->offsets;
  wf_offset_t* const out_i1 = wavefront_set->out_i1wavefront->offsets;
  wf_offset_t* const out_i2 = wavefront_set->out_i2wavefront->offsets;
  wf_offset_t* const out_d1 = wavefront_set->out_d1wavefront->offsets;
  wf_offset_t* const out_d2 = wavefront_set->out_d2wavefront->offsets;
  // In BT-pcigar
  const pcigar_t* const m_sub_bt_pcigar   = wavefront_set->in_mwavefront_sub->bt_pcigar;
  const pcigar_t* const m_open1_bt_pcigar = wavefront_set->in_mwavefront_gap1->bt_pcigar;
  const pcigar_t* const m_open2_bt_pcigar = wavefront_set->in_mwavefront_gap2->bt_pcigar;
  const pcigar_t* const i1_ext_bt_pcigar  = wavefront_set->in_i1wavefront_ext->bt_pcigar;
  const pcigar_t* const i2_ext_bt_pcigar  = wavefront_set->in_i2wavefront_ext->bt_pcigar;
  const pcigar_t* const d1_ext_bt_pcigar  = wavefront_set->in_d1wavefront_ext->bt_pcigar;
  const pcigar_t* const d2_ext_bt_pcigar  = wavefront_set->in_d2wavefront_ext->bt_pcigar;
  // In BT-prev
  const block_idx_t* const m_sub_bt_prev   = wavefront_set->in_mwavefront_sub->bt_prev;
  const block_idx_t* const m_open1_bt_prev = wavefront_set->in_mwavefront_gap1->bt_prev;
  const block_idx_t* const m_open2_bt_prev = wavefront_set->in_mwavefront_gap2->bt_prev;
  const block_idx_t* const i1_ext_bt_prev  = wavefront_set->in_i1wavefront_ext->bt_prev;
  const block_idx_t* const i2_ext_bt_prev  = wavefront_set->in_i2wavefront_ext->bt_prev;
  const block_idx_t* const d1_ext_bt_prev  = wavefront_set->in_d1wavefront_ext->bt_prev;
  const block_idx_t* const d2_ext_bt_prev  = wavefront_set->in_d2wavefront_ext->bt_prev;
  // Out BT-pcigar
  pcigar_t* const out_m_bt_pcigar   = wavefront_set->out_mwavefront->bt_pcigar;
  pcigar_t* const out_i1_bt_pcigar  = wavefront_set->out_i1wavefront->bt_pcigar;
  pcigar_t* const out_i2_bt_pcigar  = wavefront_set->out_i2wavefront->bt_pcigar;
  pcigar_t* const out_d1_bt_pcigar  = wavefront_set->out_d1wavefront->bt_pcigar;
  pcigar_t* const out_d2_bt_pcigar  = wavefront_set->out_d2wavefront->bt_pcigar;
  // Out BT-prev
  block_idx_t* const out_m_bt_prev  = wavefront_set->out_mwavefront->bt_prev;
  block_idx_t* const out_i1_bt_prev = wavefront_set->out_i1wavefront->bt_prev;
  block_idx_t* const out_i2_bt_prev = wavefront_set->out_i2wavefront->bt_prev;
  block_idx_t* const out_d1_bt_prev = wavefront_set->out_d1wavefront->bt_prev;
  block_idx_t* const out_d2_bt_prev = wavefront_set->out_d2wavefront->bt_prev;
  // Compute-Next kernel loop
  int k;
  for (k=lo;k<=hi;++k) {
    /*
     * Insertion Block
     */
    // Update I1
    const wf_offset_t ins1_o = WF_COND_FETCH_INC(m_open1,k-1,1);
    const wf_offset_t ins1_e = WF_COND_FETCH_INC(i1_ext,k-1,1);
    const wf_offset_t ins1 = MAX(ins1_o,ins1_e);
    if (ins1 >= 0) {
      if (ins1 == ins1_e) {
        out_i1_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(i1_ext_bt_pcigar[k-1]);
        out_i1_bt_prev[k] = i1_ext_bt_prev[k-1];
      } else { // ins1 == ins1_o
        out_i1_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(m_open1_bt_pcigar[k-1]);
        out_i1_bt_prev[k] = m_open1_bt_prev[k-1];
      }
    } else {
      out_i1_bt_pcigar[k] = 0;
      out_i1_bt_prev[k] = 0;
    }
    out_i1[k] = ins1;
    // Update I2
    const wf_offset_t ins2_o = WF_COND_FETCH_INC(m_open2,k-1,1);
    const wf_offset_t ins2_e = WF_COND_FETCH_INC(i2_ext,k-1,1);
    const wf_offset_t ins2 = MAX(ins2_o,ins2_e);
    if (ins2 >= 0) {
      if (ins2 == ins2_e) {
        out_i2_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(i2_ext_bt_pcigar[k-1]);
        out_i2_bt_prev[k] = i2_ext_bt_prev[k-1];
      } else { // ins2 == ins2_o
        out_i2_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(m_open2_bt_pcigar[k-1]);
        out_i2_bt_prev[k] = m_open2_bt_prev[k-1];
      }
    } else {
      out_i2_bt_pcigar[k] = 0;
      out_i2_bt_prev[k] = 0;
    }
    out_i2[k] = ins2;
    // Update I
    const wf_offset_t ins = MAX(ins1,ins2);
    /*
     * Deletion Block
     */
    // Update D1
    const wf_offset_t del1_o = WF_COND_FETCH(m_open1,k+1);
    const wf_offset_t del1_e = WF_COND_FETCH(d1_ext,k+1);
    const wf_offset_t del1 = MAX(del1_o,del1_e);
    if (del1 >= 0) {
      if (del1 == del1_e) {
        out_d1_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(d1_ext_bt_pcigar[k+1]);
        out_d1_bt_prev[k] = d1_ext_bt_prev[k+1];
      } else { // del1 == del1_o
        out_d1_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(m_open1_bt_pcigar[k+1]);
        out_d1_bt_prev[k] = m_open1_bt_prev[k+1];
      }
    } else {
      out_d1_bt_pcigar[k] = 0;
      out_d1_bt_prev[k] = 0;
    }
    out_d1[k] = del1;
    // Update D2
    const wf_offset_t del2_o = WF_COND_FETCH(m_open2,k+1);
    const wf_offset_t del2_e = WF_COND_FETCH(d2_ext,k+1);
    const wf_offset_t del2 = MAX(del2_o,del2_e);
    if (del2 >= 0) {
      if (del2 == del2_e) {
        out_d2_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(d2_ext_bt_pcigar[k+1]);
        out_d2_bt_prev[k] = d2_ext_bt_prev[k+1];
      } else { // del2 == del2_o
        out_d2_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(m_open2_bt_pcigar[k+1]);
        out_d2_bt_prev[k] = m_open2_bt_prev[k+1];
      }
    } else {
      out_d2_bt_pcigar[k] = 0;
      out_d2_bt_prev[k] = 0;
    }
    out_d2[k] = del2;
    // Update D
    const wf_offset_t del = MAX(del1,del2);
    // Update M
    const wf_offset_t sub = WF_COND_FETCH_INC(m_sub,k,1);
    const wf_offset_t max = MAX(del,MAX(sub,ins));
    if (max >= 0) {
      if (max == sub) {
        out_m_bt_pcigar[k] = m_sub_bt_pcigar[k];
        out_m_bt_prev[k] = m_sub_bt_prev[k];
      } else if (max == del2) {
        out_m_bt_pcigar[k] = out_d2_bt_pcigar[k];
        out_m_bt_prev[k] = out_d2_bt_prev[k];
      } else if (max == del1) {
        out_m_bt_pcigar[k] = out_d1_bt_pcigar[k];
        out_m_bt_prev[k] = out_d1_bt_prev[k];
      } else if (max == ins2) {
        out_m_bt_pcigar[k] = out_i2_bt_pcigar[k];
        out_m_bt_prev[k] = out_i2_bt_prev[k];
      } else { // max == ins1
        out_m_bt_pcigar[k] = out_i1_bt_pcigar[k];
        out_m_bt_prev[k] = out_i1_bt_prev[k];
      }
      // Coming from I/D -> X is fake to represent gap-close
      // Coming from M -> X is real to represent mismatch
      out_m_bt_pcigar[k] = PCIGAR_PUSH_BACK_MISMS(out_m_bt_pcigar[k]);
    } else {
      out_m_bt_pcigar[k] = 0;
      out_m_bt_prev[k] = 0;
    }
    out_m[k] = max;
  }
}
void wavefront_compute_affine2p_idm_piggyback_unbounded(
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // In Offsets
  const wf_offset_t* const m_sub   = wavefront_set->in_mwavefront_sub->offsets;
  const wf_offset_t* const m_open1 = wavefront_set->in_mwavefront_gap1->offsets;
  const wf_offset_t* const m_open2 = wavefront_set->in_mwavefront_gap2->offsets;
  const wf_offset_t* const i1_ext  = wavefront_set->in_i1wavefront_ext->offsets;
  const wf_offset_t* const i2_ext  = wavefront_set->in_i2wavefront_ext->offsets;
  const wf_offset_t* const d1_ext  = wavefront_set->in_d1wavefront_ext->offsets;
  const wf_offset_t* const d2_ext  = wavefront_set->in_d2wavefront_ext->offsets;
  // Out Offsets
  wf_offset_t* const out_m  = wavefront_set->out_mwavefront->offsets;
  wf_offset_t* const out_i1 = wavefront_set->out_i1wavefront->offsets;
  wf_offset_t* const out_i2 = wavefront_set->out_i2wavefront->offsets;
  wf_offset_t* const out_d1 = wavefront_set->out_d1wavefront->offsets;
  wf_offset_t* const out_d2 = wavefront_set->out_d2wavefront->offsets;
  // In BT-pcigar
  const pcigar_t* const m_sub_bt_pcigar   = wavefront_set->in_mwavefront_sub->bt_pcigar;
  const pcigar_t* const m_open1_bt_pcigar = wavefront_set->in_mwavefront_gap1->bt_pcigar;
  const pcigar_t* const m_open2_bt_pcigar = wavefront_set->in_mwavefront_gap2->bt_pcigar;
  const pcigar_t* const i1_ext_bt_pcigar  = wavefront_set->in_i1wavefront_ext->bt_pcigar;
  const pcigar_t* const i2_ext_bt_pcigar  = wavefront_set->in_i2wavefront_ext->bt_pcigar;
  const pcigar_t* const d1_ext_bt_pcigar  = wavefront_set->in_d1wavefront_ext->bt_pcigar;
  const pcigar_t* const d2_ext_bt_pcigar  = wavefront_set->in_d2wavefront_ext->bt_pcigar;
  // In BT-prev
  const block_idx_t* const m_sub_bt_prev   = wavefront_set->in_mwavefront_sub->bt_prev;
  const block_idx_t* const m_open1_bt_prev = wavefront_set->in_mwavefront_gap1->bt_prev;
  const block_idx_t* const m_open2_bt_prev = wavefront_set->in_mwavefront_gap2->bt_prev;
  const block_idx_t* const i1_ext_bt_prev  = wavefront_set->in_i1wavefront_ext->bt_prev;
  const block_idx_t* const i2_ext_bt_prev  = wavefront_set->in_i2wavefront_ext->bt_prev;
  const block_idx_t* const d1_ext_bt_prev  = wavefront_set->in_d1wavefront_ext->bt_prev;
  const block_idx_t* const d2_ext_bt_prev  = wavefront_set->in_d2wavefront_ext->bt_prev;
  // Out BT-pcigar
  pcigar_t* const out_m_bt_pcigar   = wavefront_set->out_mwavefront->bt_pcigar;
  pcigar_t* const out_i1_bt_pcigar  = wavefront_set->out_i1wavefront->bt_pcigar;
  pcigar_t* const out_i2_bt_pcigar  = wavefront_set->out_i2wavefront->bt_pcigar;
  pcigar_t* const out_d1_bt_pcigar  = wavefront_set->out_d1wavefront->bt_pcigar;
  pcigar_t* const out_d2_bt_pcigar  = wavefront_set->out_d2wavefront->bt_pcigar;
  // Out BT-prev
  block_idx_t* const out_m_bt_prev  = wavefront_set->out_mwavefront->bt_prev;
  block_idx_t* const out_i1_bt_prev = wavefront_set->out_i1wavefront->bt_prev;
  block_idx_t* const out_i2_bt_prev = wavefront_set->out_i2wavefront->bt_prev;
  block_idx_t* const out_d1_bt_prev = wavefront_set->out_d1wavefront->bt_prev;
  block_idx_t* const out_d2_bt_prev = wavefront_set->out_d2wavefront->bt_prev;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE
  for (k=lo;k<=hi;++k) {
    // Update I1
    const wf_offset_t ins1_o = m_open1[k-1];
    const wf_offset_t ins1_e = i1_ext[k-1];
    wf_offset_t ins1;
    if (ins1_e >= ins1_o) { // MAX (predicated by the compiler)
      ins1 = ins1_e + 1;
      out_i1_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(i1_ext_bt_pcigar[k-1]);
      out_i1_bt_prev[k] = i1_ext_bt_prev[k-1];
    } else {
      ins1 = ins1_o + 1;
      out_i1_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(m_open1_bt_pcigar[k-1]);
      out_i1_bt_prev[k] = m_open1_bt_prev[k-1];
    }
    out_i1[k] = ins1;
    // Update I2
    const wf_offset_t ins2_o = m_open2[k-1];
    const wf_offset_t ins2_e = i2_ext[k-1];
    wf_offset_t ins2;
    if (ins2_e >= ins2_o) { // MAX (predicated by the compiler)
      ins2 = ins2_e + 1;
      out_i2_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(i2_ext_bt_pcigar[k-1]);
      out_i2_bt_prev[k] = i2_ext_bt_prev[k-1];
    } else {
      ins2 = ins2_o + 1;
      out_i2_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(m_open2_bt_pcigar[k-1]);
      out_i2_bt_prev[k] = m_open2_bt_prev[k-1];
    }
    out_i2[k] = ins2;
    // Update I
    const wf_offset_t ins = MAX(ins1,ins2);
    // Update D1
    const wf_offset_t del1_o = m_open1[k+1];
    const wf_offset_t del1_e = d1_ext[k+1];
    wf_offset_t del1;
    if (del1_e >= del1_o) { // MAX (predicated by the compiler)
      del1 = del1_e;
      out_d1_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(d1_ext_bt_pcigar[k+1]);
      out_d1_bt_prev[k] = d1_ext_bt_prev[k+1];
    } else {
      del1 = del1_o;
      out_d1_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(m_open1_bt_pcigar[k+1]);
      out_d1_bt_prev[k] = m_open1_bt_prev[k+1];
    }
    out_d1[k] = del1;
    // Update D2
    const wf_offset_t del2_o = m_open2[k+1];
    const wf_offset_t del2_e = d2_ext[k+1];
    wf_offset_t del2;
    if (del2_e >= del2_o) { // MAX (predicated by the compiler)
      del2 = del2_e;
      out_d2_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(d2_ext_bt_pcigar[k+1]);
      out_d2_bt_prev[k] = d2_ext_bt_prev[k+1];
    } else {
      del2 = del2_o;
      out_d2_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(m_open2_bt_pcigar[k+1]);
      out_d2_bt_prev[k] = m_open2_bt_prev[k+1];
    }
    out_d2[k] = del2;
    // Update D
    const wf_offset_t del = MAX(del1,del2);
    // Update M
    const wf_offset_t sub = m_sub[k] + 1;
    const wf_offset_t max = MAX(del,MAX(sub,ins));
    if (max == ins1) {
      out_m_bt_pcigar[k] = out_i1_bt_pcigar[k];
      out_m_bt_prev[k] = out_i1_bt_prev[k];
    }
    if (max == ins2) {
      out_m_bt_pcigar[k] = out_i2_bt_pcigar[k];
      out_m_bt_prev[k] = out_i2_bt_prev[k];
    }
    if (max == del1) {
      out_m_bt_pcigar[k] = out_d1_bt_pcigar[k];
      out_m_bt_prev[k] = out_d1_bt_prev[k];
    }
    if (max == del2) {
      out_m_bt_pcigar[k] = out_d2_bt_pcigar[k];
      out_m_bt_prev[k] = out_d2_bt_prev[k];
    }
    if (max == sub) {
      out_m_bt_pcigar[k] = m_sub_bt_pcigar[k];
      out_m_bt_prev[k] = m_sub_bt_prev[k];
    }
    // Coming from I/D -> X is fake to represent gap-close
    // Coming from M -> X is real to represent mismatch
    out_m_bt_pcigar[k] = PCIGAR_PUSH_BACK_MISMS(out_m_bt_pcigar[k]);
    out_m[k] = max;
  }
}
/*
 * Wavefront Propagate Backtrace (attending the Piggyback)
 */
void wavefront_compute_affine2p_idm_piggyback_offload(
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi,
    wf_backtrace_buffer_t* const bt_buffer) {
  // Parameters
  pcigar_t* const out_m_bt_pcigar   = wavefront_set->out_mwavefront->bt_pcigar;
  block_idx_t* const out_m_bt_prev  = wavefront_set->out_mwavefront->bt_prev;
  pcigar_t* const out_i1_bt_pcigar  = wavefront_set->out_i1wavefront->bt_pcigar;
  block_idx_t* const out_i1_bt_prev = wavefront_set->out_i1wavefront->bt_prev;
  pcigar_t* const out_i2_bt_pcigar  = wavefront_set->out_i2wavefront->bt_pcigar;
  block_idx_t* const out_i2_bt_prev = wavefront_set->out_i2wavefront->bt_prev;
  pcigar_t* const out_d1_bt_pcigar  = wavefront_set->out_d1wavefront->bt_pcigar;
  block_idx_t* const out_d1_bt_prev = wavefront_set->out_d1wavefront->bt_prev;
  pcigar_t* const out_d2_bt_pcigar  = wavefront_set->out_d2wavefront->bt_pcigar;
  block_idx_t* const out_d2_bt_prev = wavefront_set->out_d2wavefront->bt_prev;
  // Out Offsets
  wf_offset_t* const out_m  = wavefront_set->out_mwavefront->offsets;
  wf_offset_t* const out_i1 = wavefront_set->out_i1wavefront->offsets;
  wf_offset_t* const out_i2 = wavefront_set->out_i2wavefront->offsets;
  wf_offset_t* const out_d1 = wavefront_set->out_d1wavefront->offsets;
  wf_offset_t* const out_d2 = wavefront_set->out_d2wavefront->offsets;
  // Check PCIGAR buffers full and off-load if needed
  int k;
  for (k=lo;k<=hi;++k) {
    if (PCIGAR_IS_ALMOST_FULL(out_i1_bt_pcigar[k]) && out_i1[k]>=0) {
      wf_backtrace_buffer_store_block_bt(bt_buffer,out_i1_bt_pcigar+k,out_i1_bt_prev+k);
    }
    if (PCIGAR_IS_ALMOST_FULL(out_i2_bt_pcigar[k]) && out_i2[k]>=0) {
      wf_backtrace_buffer_store_block_bt(bt_buffer,out_i2_bt_pcigar+k,out_i2_bt_prev+k);
    }
    if (PCIGAR_IS_ALMOST_FULL(out_d1_bt_pcigar[k]) && out_d1[k]>=0) {
      wf_backtrace_buffer_store_block_bt(bt_buffer,out_d1_bt_pcigar+k,out_d1_bt_prev+k);
    }
    if (PCIGAR_IS_ALMOST_FULL(out_d2_bt_pcigar[k]) && out_d2[k]>=0) {
      wf_backtrace_buffer_store_block_bt(bt_buffer,out_d2_bt_pcigar+k,out_d2_bt_prev+k);
    }
    if (PCIGAR_IS_ALMOST_FULL(out_m_bt_pcigar[k]) && out_m[k]>=0) {
      wf_backtrace_buffer_store_block_bt(bt_buffer,out_m_bt_pcigar+k,out_m_bt_prev+k);
    }
  }
}
/*
 * Compute Wavefront (IDM)
 */
void wavefront_compute_affine2p_idm(
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Compute loop peeling limits [max_lo,min_hi] (dense region where all the offsets exists
  int min_hi, max_lo;
  wavefront_compute_limits_dense(wavefront_set,gap_affine_2p,&max_lo,&min_hi);
  // Compute wavefronts (prologue)
  wavefront_compute_affine2p_idm_bounded(wavefront_set,lo,max_lo-1);
  // Compute wavefronts (core)
  wavefront_compute_affine2p_idm_unbounded(wavefront_set,max_lo,min_hi);
  // Compute wavefronts (epilogue)
  wavefront_compute_affine2p_idm_bounded(wavefront_set,min_hi+1,hi);
}
void wavefront_compute_affine2p_idm_piggyback(
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi,
    wf_backtrace_buffer_t* const bt_buffer) {
  // Compute loop peeling limits [max_lo,min_hi] (dense region where all the offsets exists)
  int min_hi, max_lo;
  wavefront_compute_limits_dense(wavefront_set,gap_affine_2p,&max_lo,&min_hi);
  // Compute wavefronts (prologue)
  wavefront_compute_affine2p_idm_piggyback_bounded(wavefront_set,lo,max_lo-1);
  // Compute wavefronts (core)
  wavefront_compute_affine2p_idm_piggyback_unbounded(wavefront_set,max_lo,min_hi);
  // Compute wavefronts (epilogue)
  wavefront_compute_affine2p_idm_piggyback_bounded(wavefront_set,min_hi+1,hi);
  // Offload Backtrace
  wavefront_compute_affine2p_idm_piggyback_offload(wavefront_set,lo,hi,bt_buffer);
}
/*
 * Compute next wavefront
 */
void wavefront_compute_affine2p(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Select wavefronts
  wavefront_set_t wavefront_set;
  wavefront_aligner_fetch_input(wf_aligner,&wavefront_set,score);
  // Check null wavefronts
  if (wavefront_set.in_mwavefront_sub->null &&
      wavefront_set.in_mwavefront_gap1->null &&
      wavefront_set.in_mwavefront_gap2->null &&
      wavefront_set.in_i1wavefront_ext->null &&
      wavefront_set.in_i2wavefront_ext->null &&
      wavefront_set.in_d1wavefront_ext->null &&
      wavefront_set.in_d2wavefront_ext->null) {
    wavefront_aligner_allocate_output_null(wf_aligner,score); // Null s-wavefront
    return;
  }
  // Set limits
  int hi, lo;
  wavefront_compute_limits(&wavefront_set,gap_affine_2p,&lo,&hi);
  // Allocate wavefronts
  wavefront_aligner_allocate_output(wf_aligner,&wavefront_set,score,lo,hi);
  // Compute next wavefront
  if (wf_aligner->wf_components.bt_piggyback) {
    wavefront_compute_affine2p_idm_piggyback(&wavefront_set,lo,hi,wf_aligner->wf_components.bt_buffer);
  } else {
    wavefront_compute_affine2p_idm(&wavefront_set,lo,hi);
  }
  // Trim wavefront ends
  wavefront_aligner_trim_ends(wf_aligner,score);
}



