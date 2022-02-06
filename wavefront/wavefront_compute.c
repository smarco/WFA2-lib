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
 * DESCRIPTION: WaveFront alignment module for computing wavefronts
 */

#include "utils/string_padded.h"
#include "gap_affine2p/affine2p_penalties.h"
#include "wavefront_compute.h"

/*
 * Compute limits
 */
void wavefront_compute_limits(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    int* const lo,
    int* const hi) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  const wavefront_t* const m_misms = wavefront_set->in_mwavefront_misms;
  const wavefront_t* const m_open1 = wavefront_set->in_mwavefront_open1;
  // Gap-lineal
  int min_lo = m_misms->lo;
  int max_hi = m_misms->hi;
  if (min_lo > m_open1->lo-1) min_lo = m_open1->lo-1;
  if (max_hi < m_open1->hi+1) max_hi = m_open1->hi+1;
  if (distance_metric==gap_lineal) {
    *lo = min_lo;
    *hi = max_hi;
    return;
  }
  // Parameters
  const wavefront_t* const i1_ext = wavefront_set->in_i1wavefront_ext;
  const wavefront_t* const d1_ext = wavefront_set->in_d1wavefront_ext;
  // Gap-affine
  if (min_lo > i1_ext->lo+1) min_lo = i1_ext->lo+1;
  if (max_hi < i1_ext->hi+1) max_hi = i1_ext->hi+1;
  if (min_lo > d1_ext->lo-1) min_lo = d1_ext->lo-1;
  if (max_hi < d1_ext->hi-1) max_hi = d1_ext->hi-1;
  if (distance_metric==gap_affine) {
    *lo = min_lo;
    *hi = max_hi;
    return;
  }
  // Parameters
  const wavefront_t* const m_open2 = wavefront_set->in_mwavefront_open2;
  const wavefront_t* const i2_ext = wavefront_set->in_i2wavefront_ext;
  const wavefront_t* const d2_ext = wavefront_set->in_d2wavefront_ext;
  // Gap-affine-2p
  if (min_lo > m_open2->lo-1) min_lo = m_open2->lo-1;
  if (max_hi < m_open2->hi+1) max_hi = m_open2->hi+1;
  if (min_lo > i2_ext->lo+1) min_lo = i2_ext->lo+1;
  if (max_hi < i2_ext->hi+1) max_hi = i2_ext->hi+1;
  if (min_lo > d2_ext->lo-1) min_lo = d2_ext->lo-1;
  if (max_hi < d2_ext->hi-1) max_hi = d2_ext->hi-1;
  *lo = min_lo;
  *hi = max_hi;
}
/*
 * Input wavefronts (fetch)
 */
wavefront_t* wavefront_compute_get_mwavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->mwavefronts[score] == NULL ||
          wf_components->mwavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->mwavefronts[score];
}
wavefront_t* wavefront_compute_get_i1wavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->i1wavefronts[score] == NULL ||
          wf_components->i1wavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->i1wavefronts[score];
}
wavefront_t* wavefront_compute_get_i2wavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->i2wavefronts[score] == NULL ||
          wf_components->i2wavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->i2wavefronts[score];
}
wavefront_t* wavefront_compute_get_d1wavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->d1wavefronts[score] == NULL ||
          wf_components->d1wavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->d1wavefronts[score];
}
wavefront_t* wavefront_compute_get_d2wavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->d2wavefronts[score] == NULL ||
          wf_components->d2wavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->d2wavefronts[score];
}
void wavefront_compute_fetch_input(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    const int score) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Compute scores
  const wavefronts_penalties_t* const penalties = &(wf_aligner->penalties);
  int mismatch = score - penalties->mismatch;
  int gap_open1 = score - penalties->gap_opening1 - penalties->gap_extension1;
  int gap_extend1 = score - penalties->gap_extension1;
  int gap_open2 = score - penalties->gap_opening2 - penalties->gap_extension2;
  int gap_extend2 = score - penalties->gap_extension2;
  // Modular wavefront
  if (wf_components->memory_modular) {
    const int max_score_scope = wf_components->max_score_scope;
    if (mismatch > 0) mismatch =  mismatch % max_score_scope;
    if (gap_open1 > 0) gap_open1 = gap_open1 % max_score_scope;
    if (gap_extend1 > 0) gap_extend1 = gap_extend1 % max_score_scope;
    if (gap_open2 > 0) gap_open2 = gap_open2 % max_score_scope;
    if (gap_extend2 > 0) gap_extend2 = gap_extend2 % max_score_scope;
  }
  // Fetch wavefronts
  wavefront_set->in_mwavefront_misms = wavefront_compute_get_mwavefront(wf_components,mismatch);
  wavefront_set->in_mwavefront_open1 = wavefront_compute_get_mwavefront(wf_components,gap_open1);
  if (distance_metric==gap_lineal) return;
  wavefront_set->in_i1wavefront_ext = wavefront_compute_get_i1wavefront(wf_components,gap_extend1);
  wavefront_set->in_d1wavefront_ext = wavefront_compute_get_d1wavefront(wf_components,gap_extend1);
  if (distance_metric==gap_affine) return;
  wavefront_set->in_mwavefront_open2 = wavefront_compute_get_mwavefront(wf_components,gap_open2);
  wavefront_set->in_i2wavefront_ext = wavefront_compute_get_i2wavefront(wf_components,gap_extend2);
  wavefront_set->in_d2wavefront_ext = wavefront_compute_get_d2wavefront(wf_components,gap_extend2);
}
/*
 * Output wavefronts (allocate)
 */
void wavefront_compute_free_output(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_slab_t* const wavefront_slab = wf_aligner->wavefront_slab;
  // Free
  if (wf_components->mwavefronts[score]) wavefront_slab_free(wavefront_slab,wf_components->mwavefronts[score]);
  if (distance_metric==gap_lineal) return;
  if (wf_components->i1wavefronts[score]) wavefront_slab_free(wavefront_slab,wf_components->i1wavefronts[score]);
  if (wf_components->d1wavefronts[score]) wavefront_slab_free(wavefront_slab,wf_components->d1wavefronts[score]);
  if (distance_metric==gap_affine) return;
  if (wf_components->i2wavefronts[score]) wavefront_slab_free(wavefront_slab,wf_components->i2wavefronts[score]);
  if (wf_components->d2wavefronts[score]) wavefront_slab_free(wavefront_slab,wf_components->d2wavefronts[score]);
}
void wavefront_compute_allocate_output_null(
    wavefront_aligner_t* const wf_aligner,
    int score) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Modular wavefront
  if (wf_components->memory_modular) {
    score = score % wf_components->max_score_scope;
    wavefront_compute_free_output(wf_aligner,score);
  }
  // Nullify Wavefronts
  wf_components->mwavefronts[score] = NULL;
  if (distance_metric==gap_lineal) return;
  wf_components->i1wavefronts[score] = NULL;
  wf_components->d1wavefronts[score] = NULL;
  if (distance_metric==gap_affine) return;
  wf_components->i2wavefronts[score] = NULL;
  wf_components->d2wavefronts[score] = NULL;
}
void wavefront_compute_allocate_output(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    int score,
    const int lo,
    const int hi) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_slab_t* const wavefront_slab = wf_aligner->wavefront_slab;
  const int max_score_scope = wf_components->max_score_scope;
  // Compute effective hi/lo dimensions (after padding to avoid compute-kernel peeling)
  const int effective_lo = lo - (max_score_scope+1);
  const int effective_hi = hi + (max_score_scope+1);
  const int padded_lo = MIN(effective_lo,wf_components->historic_min_lo);
  const int padded_hi = MAX(effective_hi,wf_components->historic_max_hi);
  wf_components->historic_min_lo = padded_lo;
  wf_components->historic_max_hi = padded_hi;
  // Resize null/victim wavefronts
  wavefront_components_resize_null__victim(wf_components,padded_lo,padded_hi);
  // Modular wavefront
  if (wf_components->memory_modular) {
    score = score % wf_components->max_score_scope;
    wavefront_compute_free_output(wf_aligner,score);
  }
  // Allocate M-Wavefront
  wavefront_set->out_mwavefront = wavefront_slab_allocate(wavefront_slab,padded_lo,padded_hi);
  wf_components->mwavefronts[score] = wavefront_set->out_mwavefront;
  wf_components->mwavefronts[score]->lo = lo;
  wf_components->mwavefronts[score]->hi = hi;
  if (distance_metric==gap_lineal) return;
  // Allocate I1-Wavefront
  if (!wavefront_set->in_mwavefront_open1->null || !wavefront_set->in_i1wavefront_ext->null) {
    wavefront_set->out_i1wavefront = wavefront_slab_allocate(wavefront_slab,padded_lo,padded_hi);
    wf_components->i1wavefronts[score] = wavefront_set->out_i1wavefront;
    wf_components->i1wavefronts[score]->lo = lo;
    wf_components->i1wavefronts[score]->hi = hi;
  } else {
    wavefront_set->out_i1wavefront = wf_components->wavefront_victim;
    wf_components->i1wavefronts[score] = NULL;
  }
  // Allocate D1-Wavefront
  if (!wavefront_set->in_mwavefront_open1->null || !wavefront_set->in_d1wavefront_ext->null) {
    wavefront_set->out_d1wavefront = wavefront_slab_allocate(wavefront_slab,padded_lo,padded_hi);
    wf_components->d1wavefronts[score] = wavefront_set->out_d1wavefront;
    wf_components->d1wavefronts[score]->lo = lo;
    wf_components->d1wavefronts[score]->hi = hi;
  } else {
    wavefront_set->out_d1wavefront = wf_components->wavefront_victim;
    wf_components->d1wavefronts[score] = NULL;
  }
  if (distance_metric==gap_affine) return;
  // Allocate I2-Wavefront
  if (!wavefront_set->in_mwavefront_open2->null || !wavefront_set->in_i2wavefront_ext->null) {
    wavefront_set->out_i2wavefront = wavefront_slab_allocate(wavefront_slab,padded_lo,padded_hi);
    wf_components->i2wavefronts[score] = wavefront_set->out_i2wavefront;
    wf_components->i2wavefronts[score]->lo = lo;
    wf_components->i2wavefronts[score]->hi = hi;
  } else {
    wavefront_set->out_i2wavefront = wf_components->wavefront_victim;
    wf_components->i2wavefronts[score] = NULL;
  }
  // Allocate D2-Wavefront
  if (!wavefront_set->in_mwavefront_open2->null || !wavefront_set->in_d2wavefront_ext->null) {
    wavefront_set->out_d2wavefront = wavefront_slab_allocate(wavefront_slab,padded_lo,padded_hi);
    wf_components->d2wavefronts[score] = wavefront_set->out_d2wavefront;
    wf_components->d2wavefronts[score]->lo = lo;
    wf_components->d2wavefronts[score]->hi = hi;
  } else {
    wavefront_set->out_d2wavefront = wf_components->wavefront_victim;
    wf_components->d2wavefronts[score] = NULL;
  }
}
/*
 * Initialize wavefronts ends
 */
void wavefront_compute_init_ends_wf_lower(
    wavefront_t* const wavefront,
    const int min_lo) {
  // Check initialization (lowest element)
  if (wavefront->wf_elements_init_min <= min_lo) return;
  // Initialize lower elements
  wf_offset_t* const offsets = wavefront->offsets;
  const int min_init = MIN(wavefront->wf_elements_init_min,wavefront->lo);
  int k;
  for (k=min_lo;k<min_init;++k) {
    offsets[k] = WAVEFRONT_OFFSET_NULL;
  }
  // Set new minimum
  wavefront->wf_elements_init_min = min_lo;
}
void wavefront_compute_init_ends_wf_higher(
    wavefront_t* const wavefront,
    const int max_hi) {
  // Check initialization (highest element)
  if (wavefront->wf_elements_init_max >= max_hi) return;
  // Initialize lower elements
  wf_offset_t* const offsets = wavefront->offsets;
  const int max_init = MAX(wavefront->wf_elements_init_max,wavefront->hi);
  int k;
  for (k=max_init+1;k<=max_hi;++k) {
    offsets[k] = WAVEFRONT_OFFSET_NULL;
  }
  // Set new maximum
  wavefront->wf_elements_init_max = max_hi;
}
void wavefront_compute_init_ends(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Init missing elements, instead of loop peeling (M)
  const bool m_misms_null = wavefront_set->in_mwavefront_misms->null;
  const bool m_gap1_null = wavefront_set->in_mwavefront_open1->null;
  if (!m_misms_null) {
    wavefront_compute_init_ends_wf_higher(wavefront_set->in_mwavefront_misms,hi);
    wavefront_compute_init_ends_wf_lower(wavefront_set->in_mwavefront_misms,lo);
  }
  if (!m_gap1_null) {
    wavefront_compute_init_ends_wf_higher(wavefront_set->in_mwavefront_open1,hi+1);
    wavefront_compute_init_ends_wf_lower(wavefront_set->in_mwavefront_open1,lo-1);
  }
  if (distance_metric==gap_lineal) return;
  // Init missing elements, instead of loop peeling (Open1/I1/D1)
  const bool i1_ext_null = wavefront_set->in_i1wavefront_ext->null;
  const bool d1_ext_null = wavefront_set->in_d1wavefront_ext->null;
  if (!i1_ext_null) {
    wavefront_compute_init_ends_wf_higher(wavefront_set->in_i1wavefront_ext,hi);
    wavefront_compute_init_ends_wf_lower(wavefront_set->in_i1wavefront_ext,lo-1);
  }
  if (!d1_ext_null) {
    wavefront_compute_init_ends_wf_higher(wavefront_set->in_d1wavefront_ext,hi+1);
    wavefront_compute_init_ends_wf_lower(wavefront_set->in_d1wavefront_ext,lo);
  }
  if (distance_metric==gap_affine) return;
  // Init missing elements, instead of loop peeling (Open2/I2/D2)
  const bool m_gap2_null = wavefront_set->in_mwavefront_open2->null;
  const bool i2_ext_null = wavefront_set->in_i2wavefront_ext->null;
  const bool d2_ext_null = wavefront_set->in_d2wavefront_ext->null;
  if (!m_gap2_null) {
    wavefront_compute_init_ends_wf_higher(wavefront_set->in_mwavefront_open2,hi+1);
    wavefront_compute_init_ends_wf_lower(wavefront_set->in_mwavefront_open2,lo-1);
  }
  if (!i2_ext_null) {
    wavefront_compute_init_ends_wf_higher(wavefront_set->in_i2wavefront_ext,hi);
    wavefront_compute_init_ends_wf_lower(wavefront_set->in_i2wavefront_ext,lo-1);
  }
  if (!d2_ext_null) {
    wavefront_compute_init_ends_wf_higher(wavefront_set->in_d2wavefront_ext,hi+1);
    wavefront_compute_init_ends_wf_lower(wavefront_set->in_d2wavefront_ext,lo);
  }
}
/*
 * Trim wavefronts ends
 */
void wavefront_compute_trim_ends_wf(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  wf_offset_t* const offsets = wavefront->offsets;
  // Trim from hi
  int k;
  const int lo = wavefront->lo;
  for (k=wavefront->hi;k>=lo;--k) {
    // Fetch offset
    const wf_offset_t offset = offsets[k];
    // Check boundaries
    const uint32_t h = WAVEFRONT_H(k,offset); // Make unsigned to avoid checking negative
    const uint32_t v = WAVEFRONT_V(k,offset); // Make unsigned to avoid checking negative
    if (h <= text_length && v <= pattern_length) break;
  }
  wavefront->hi = k; // Set new hi
  wavefront->wf_elements_init_max = k;
  // Trim from lo
  const int hi = wavefront->hi;
  for (k=wavefront->lo;k<=hi;++k) {
    // Fetch offset
    const wf_offset_t offset = offsets[k];
    // Check boundaries
    const uint32_t h = WAVEFRONT_H(k,offset); // Make unsigned to avoid checking negative
    const uint32_t v = WAVEFRONT_V(k,offset); // Make unsigned to avoid checking negative
    if (h <= text_length && v <= pattern_length) break;
  }
  wavefront->lo = k; // Set new lo
  wavefront->wf_elements_init_min = k;
}
void wavefront_compute_trim_ends(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Trim ends from non-null WFs
  if (wavefront_set->out_mwavefront) wavefront_compute_trim_ends_wf(wf_aligner,wavefront_set->out_mwavefront);
  if (distance_metric==gap_lineal) return;
  if (wavefront_set->out_i1wavefront) wavefront_compute_trim_ends_wf(wf_aligner,wavefront_set->out_i1wavefront);
  if (wavefront_set->out_d1wavefront) wavefront_compute_trim_ends_wf(wf_aligner,wavefront_set->out_d1wavefront);
  if (distance_metric==gap_affine) return;
  if (wavefront_set->out_i2wavefront) wavefront_compute_trim_ends_wf(wf_aligner,wavefront_set->out_i2wavefront);
  if (wavefront_set->out_d2wavefront) wavefront_compute_trim_ends_wf(wf_aligner,wavefront_set->out_d2wavefront);
}
/*
 * Backtrace offloading
 */
void wavefront_compute_offload_backtrace_blocks_selective(
    wf_offset_t* const out_offsets,
    pcigar_t* const out_bt_pcigar,
    bt_block_idx_t* const out_bt_prev,
    const int lo,
    const int hi,
    const pcigar_t occupation_mask,
    wf_backtrace_buffer_t* const bt_buffer) {
  // Fetch BT-buffer free memory
  int bt_blocks_available;
  bt_block_t* bt_block_mem;
  bt_block_idx_t global_pos = wf_backtrace_buffer_get_mem(bt_buffer,&bt_block_mem,&bt_blocks_available);
  bt_block_idx_t current_pos = global_pos;
  const int max_pos = current_pos + bt_blocks_available;
  // Check PCIGAR buffers full and off-load if needed
  int k;
  for (k=lo;k<=hi;++k) {
    if (PCIGAR_IS_UTILISED(out_bt_pcigar[k],occupation_mask) && out_offsets[k]>=0) {
      // Store
      bt_block_mem->pcigar = out_bt_pcigar[k];
      bt_block_mem->prev_idx = out_bt_prev[k];
      bt_block_mem++;
      // Reset
      out_bt_pcigar[k] = 0;
      out_bt_prev[k] = current_pos;
      current_pos++;
      // Update pos
      if (current_pos >= max_pos) {
        wf_backtrace_buffer_add_used(bt_buffer,current_pos-global_pos);
        global_pos = wf_backtrace_buffer_get_mem(bt_buffer,&bt_block_mem,&bt_blocks_available);
      }
    }
  }
  wf_backtrace_buffer_add_used(bt_buffer,current_pos-global_pos);
}
void wavefront_compute_offload_backtrace_blocks_all(
    wf_offset_t* const out_offsets,
    pcigar_t* const out_bt_pcigar,
    bt_block_idx_t* const out_bt_prev,
    const int lo,
    const int hi,
    wf_backtrace_buffer_t* const bt_buffer) {
  // Offload all BT-blocks (no matter the occupancy)
  int k = lo;
  while (k <= hi) {
    // Fetch BT-buffer free memory
    int bt_blocks_available;
    bt_block_t* bt_block_mem;
    const bt_block_idx_t global_pos = wf_backtrace_buffer_get_mem(bt_buffer,&bt_block_mem,&bt_blocks_available);
    // Offload as many BT-blocks as possible
    const int offloaded_blocks = MIN(hi-k+1,bt_blocks_available);
    int offset;
    PRAGMA_LOOP_VECTORIZE
    for (offset=0;offset<offloaded_blocks;++offset) {
      // Store
      bt_block_mem[offset].pcigar = out_bt_pcigar[k+offset];
      bt_block_mem[offset].prev_idx = out_bt_prev[k+offset];
      // Reset
      out_bt_pcigar[k+offset] = 0;
      out_bt_prev[k+offset] = global_pos+offset;
    }
    // Update offloaded
    wf_backtrace_buffer_add_used(bt_buffer,offloaded_blocks);
    k += offloaded_blocks;
  }
}
int wavefront_compute_offload_backtrace_blocks(
    wavefront_aligner_t* const wf_aligner,
    wf_offset_t* const out_offsets,
    pcigar_t* const out_bt_pcigar,
    bt_block_idx_t* const out_bt_prev,
    const int lo,
    const int hi) {
  // Parameters
  const wavefront_memory_t wavefront_memory = wf_aligner->memory_mode;
  wf_backtrace_buffer_t* const bt_buffer = wf_aligner->wf_components.bt_buffer;
  // Select memory-mode
  switch (wavefront_memory) {
    case wavefront_memory_med:
      wavefront_compute_offload_backtrace_blocks_all(
          out_offsets,out_bt_pcigar,out_bt_prev,lo,hi,bt_buffer);
      return 0; // Maximum occupancy (all empty)
      break;
    case wavefront_memory_low:
      wavefront_compute_offload_backtrace_blocks_selective(
          out_offsets,out_bt_pcigar,out_bt_prev,
          lo,hi,PCIGAR_HALF_FULL_MASK,bt_buffer);
      return PCIGAR_MAX_LENGTH/2; // Maximum occupancy (at least 2 ops)
    case wavefront_memory_ultralow:
      wavefront_compute_offload_backtrace_blocks_selective(
          out_offsets,out_bt_pcigar,out_bt_prev,
          lo,hi,PCIGAR_FULL_MASK,bt_buffer);
      return PCIGAR_MAX_LENGTH-2; // Maximum occupancy (at least 2 ops)
    default:
      fprintf(stderr,"[WFA::compute] Wrong memory-mode\n");
      exit(1);
      return 0;
  }
}
void wavefront_compute_offload_backtrace_occupation(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Select distance metric
  int occ_max_m = 0;
  int occ_max_i1 = 0, occ_max_i2 = 0;
  int occ_max_d1 = 0, occ_max_d2 = 0;
  if (distance_metric==gap_lineal) {
    // Parameters
    const wavefront_t* const m_misms = wavefront_set->in_mwavefront_misms;
    const wavefront_t* const m_open1 = wavefront_set->in_mwavefront_open1;
    // Compute BT occupancy maximum
    if (!m_open1->null) occ_max_d1 = m_open1->bt_occupancy_max + 1;
    occ_max_m = m_misms->bt_occupancy_max;
    if (occ_max_m < occ_max_d1) occ_max_m = occ_max_d1;
    ++occ_max_m;
    // Set new occupancy
    wavefront_set->out_mwavefront->bt_occupancy_max = occ_max_m;
  } else if (distance_metric==gap_affine) {
    // Parameters
    const wavefront_t* const m_misms = wavefront_set->in_mwavefront_misms;
    const wavefront_t* const m_open1 = wavefront_set->in_mwavefront_open1;
    const wavefront_t* const i1_ext = wavefront_set->in_i1wavefront_ext;
    const wavefront_t* const d1_ext = wavefront_set->in_d1wavefront_ext;
    // Compute BT occupancy maximum
    if (!m_open1->null) {
      occ_max_i1 = m_open1->bt_occupancy_max + 1;
      occ_max_d1 = m_open1->bt_occupancy_max + 1;
    }
    if (!i1_ext->null) occ_max_i1 = MAX(occ_max_i1,i1_ext->bt_occupancy_max+1);
    if (!d1_ext->null) occ_max_d1 = MAX(occ_max_d1,d1_ext->bt_occupancy_max+1);
    occ_max_m = m_misms->bt_occupancy_max;
    if (occ_max_m < occ_max_i1) occ_max_m = occ_max_i1;
    if (occ_max_m < occ_max_d1) occ_max_m = occ_max_d1;
    ++occ_max_m;
    // Set new occupancy
    wavefront_set->out_i1wavefront->bt_occupancy_max = occ_max_i1;
    wavefront_set->out_d1wavefront->bt_occupancy_max = occ_max_d1;
    wavefront_set->out_mwavefront->bt_occupancy_max = occ_max_m;
  } else { // distance_metric == gap_affine_2p
    // Parameters
    const wavefront_t* const m_misms = wavefront_set->in_mwavefront_misms;
    const wavefront_t* const m_open1 = wavefront_set->in_mwavefront_open1;
    const wavefront_t* const m_open2 = wavefront_set->in_mwavefront_open2;
    const wavefront_t* const i1_ext = wavefront_set->in_i1wavefront_ext;
    const wavefront_t* const i2_ext = wavefront_set->in_i2wavefront_ext;
    const wavefront_t* const d1_ext = wavefront_set->in_d1wavefront_ext;
    const wavefront_t* const d2_ext = wavefront_set->in_d2wavefront_ext;
    // Compute BT occupancy maximum (I)
    if (!m_open1->null) {
      occ_max_i1 = m_open1->bt_occupancy_max + 1;
      occ_max_d1 = m_open1->bt_occupancy_max + 1;
    }
    // Compute BT occupancy maximum (D)
    if (!i1_ext->null) occ_max_i1 = MAX(occ_max_i1,i1_ext->bt_occupancy_max+1);
    if (!d1_ext->null) occ_max_d1 = MAX(occ_max_d1,d1_ext->bt_occupancy_max+1);
    if (!m_open2->null) {
      occ_max_i2 = m_open2->bt_occupancy_max + 1;
      occ_max_d2 = m_open2->bt_occupancy_max + 1;
    }
    if (!i2_ext->null) occ_max_i2 = MAX(occ_max_i2,i2_ext->bt_occupancy_max+1);
    if (!d2_ext->null) occ_max_d2 = MAX(occ_max_d2,d2_ext->bt_occupancy_max+1);
    // Compute BT occupancy maximum (M)
    occ_max_m = m_misms->bt_occupancy_max;
    if (occ_max_m < occ_max_i1) occ_max_m = occ_max_i1;
    if (occ_max_m < occ_max_i2) occ_max_m = occ_max_i2;
    if (occ_max_m < occ_max_d1) occ_max_m = occ_max_d1;
    if (occ_max_m < occ_max_d2) occ_max_m = occ_max_d2;
    ++occ_max_m;
    // Set new occupancy
    wavefront_set->out_i1wavefront->bt_occupancy_max = occ_max_i1;
    wavefront_set->out_i2wavefront->bt_occupancy_max = occ_max_i2;
    wavefront_set->out_d1wavefront->bt_occupancy_max = occ_max_d1;
    wavefront_set->out_d2wavefront->bt_occupancy_max = occ_max_d2;
    wavefront_set->out_mwavefront->bt_occupancy_max = occ_max_m;
  }
}
void wavefront_compute_offload_backtrace(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Compute maximum occupancy
  wavefront_compute_offload_backtrace_occupation(wf_aligner,wavefront_set);
  /*
   * Offload if necessary (Edit or Gap-Lineal)
   */
  const wavefront_t* const out_m = wavefront_set->out_mwavefront;
  if (!out_m->null && out_m->bt_occupancy_max >= PCIGAR_MAX_LENGTH) {
    wf_offset_t* const out_m  = wavefront_set->out_mwavefront->offsets;
    pcigar_t* const out_m_bt_pcigar = wavefront_set->out_mwavefront->bt_pcigar;
    bt_block_idx_t* const out_m_bt_prev = wavefront_set->out_mwavefront->bt_prev;
    wavefront_set->out_mwavefront->bt_occupancy_max =
        wavefront_compute_offload_backtrace_blocks(
            wf_aligner,out_m,out_m_bt_pcigar,out_m_bt_prev,lo,hi);
  }
  if (distance_metric==gap_lineal) return;
  /*
   * Offload if necessary (Gap-Affine)
   */
  const wavefront_t* const out_i1 = wavefront_set->out_i1wavefront;
  if (!out_i1->null && out_i1->bt_occupancy_max >= PCIGAR_MAX_LENGTH) {
    wf_offset_t* const out_i1 = wavefront_set->out_i1wavefront->offsets;
    pcigar_t* const out_i1_bt_pcigar = wavefront_set->out_i1wavefront->bt_pcigar;
    bt_block_idx_t* const out_i1_bt_prev = wavefront_set->out_i1wavefront->bt_prev;
    wavefront_set->out_i1wavefront->bt_occupancy_max =
        wavefront_compute_offload_backtrace_blocks(
            wf_aligner,out_i1,out_i1_bt_pcigar,out_i1_bt_prev,lo,hi);
  }
  const wavefront_t* const out_d1 = wavefront_set->out_d1wavefront;
  if (!out_d1->null && out_d1->bt_occupancy_max >= PCIGAR_MAX_LENGTH) {
    wf_offset_t* const out_d1 = wavefront_set->out_d1wavefront->offsets;
    pcigar_t* const out_d1_bt_pcigar  = wavefront_set->out_d1wavefront->bt_pcigar;
    bt_block_idx_t* const out_d1_bt_prev = wavefront_set->out_d1wavefront->bt_prev;
    wavefront_set->out_d1wavefront->bt_occupancy_max =
        wavefront_compute_offload_backtrace_blocks(
            wf_aligner,out_d1,out_d1_bt_pcigar,out_d1_bt_prev,lo,hi);
  }
  if (distance_metric==gap_affine) return;
  /*
   * Offload if necessary (Gap-Affine-2p)
   */
  const wavefront_t* const out_i2 = wavefront_set->out_i2wavefront;
  if (!out_i2->null && out_i2->bt_occupancy_max >= PCIGAR_MAX_LENGTH) {
    wf_offset_t* const out_i2 = wavefront_set->out_i2wavefront->offsets;
    pcigar_t* const out_i2_bt_pcigar = wavefront_set->out_i2wavefront->bt_pcigar;
    bt_block_idx_t* const out_i2_bt_prev = wavefront_set->out_i2wavefront->bt_prev;
    wavefront_set->out_i2wavefront->bt_occupancy_max =
        wavefront_compute_offload_backtrace_blocks(
            wf_aligner,out_i2,out_i2_bt_pcigar,out_i2_bt_prev,lo,hi);
  }
  const wavefront_t* const out_d2 = wavefront_set->out_d2wavefront;
  if (!out_d2->null && out_d2->bt_occupancy_max >= PCIGAR_MAX_LENGTH) {
    wf_offset_t* const out_d2 = wavefront_set->out_d2wavefront->offsets;
    pcigar_t* const out_d2_bt_pcigar = wavefront_set->out_d2wavefront->bt_pcigar;
    bt_block_idx_t* const out_d2_bt_prev = wavefront_set->out_d2wavefront->bt_prev;
    wavefront_set->out_d2wavefront->bt_occupancy_max =
        wavefront_compute_offload_backtrace_blocks(
            wf_aligner,out_d2,out_d2_bt_pcigar,out_d2_bt_prev,lo,hi);
  }
}

