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
    const wavefront_set_t* const wavefront_set,
    const distance_metric_t distance_metric,
    int* const lo,
    int* const hi) {
  // Gap-lineal
  int min_lo = wavefront_set->in_mwavefront_sub->lo;
  int max_hi = wavefront_set->in_mwavefront_sub->hi;
  if (min_lo > wavefront_set->in_mwavefront_gap1->lo) min_lo = wavefront_set->in_mwavefront_gap1->lo;
  if (max_hi < wavefront_set->in_mwavefront_gap1->hi) max_hi = wavefront_set->in_mwavefront_gap1->hi;
  if (distance_metric==gap_lineal) {
    *lo = min_lo-1; *hi = max_hi+1; return;
  }
  // Gap-affine
  if (min_lo > wavefront_set->in_i1wavefront_ext->lo) min_lo = wavefront_set->in_i1wavefront_ext->lo;
  if (min_lo > wavefront_set->in_d1wavefront_ext->lo) min_lo = wavefront_set->in_d1wavefront_ext->lo;
  if (max_hi < wavefront_set->in_i1wavefront_ext->hi) max_hi = wavefront_set->in_i1wavefront_ext->hi;
  if (max_hi < wavefront_set->in_d1wavefront_ext->hi) max_hi = wavefront_set->in_d1wavefront_ext->hi;
  if (distance_metric==gap_affine) {
    *lo = min_lo-1; *hi = max_hi+1; return;
  }
  // Gap-affine-2p
  if (min_lo > wavefront_set->in_mwavefront_gap2->lo) min_lo = wavefront_set->in_mwavefront_gap2->lo;
  if (min_lo > wavefront_set->in_i2wavefront_ext->lo) min_lo = wavefront_set->in_i2wavefront_ext->lo;
  if (min_lo > wavefront_set->in_d2wavefront_ext->lo) min_lo = wavefront_set->in_d2wavefront_ext->lo;
  if (max_hi < wavefront_set->in_mwavefront_gap2->hi) max_hi = wavefront_set->in_mwavefront_gap2->hi;
  if (max_hi < wavefront_set->in_i2wavefront_ext->hi) max_hi = wavefront_set->in_i2wavefront_ext->hi;
  if (max_hi < wavefront_set->in_d2wavefront_ext->hi) max_hi = wavefront_set->in_d2wavefront_ext->hi;
  *lo = min_lo-1;
  *hi = max_hi+1;
}
void wavefront_compute_limits_dense(
    const wavefront_set_t* const wavefront_set,
    const distance_metric_t distance_metric,
    int* const lo,
    int* const hi) {
  // Gap-lineal
  int max_lo = wavefront_set->in_mwavefront_sub->lo;
  int min_hi = wavefront_set->in_mwavefront_sub->hi;
  if (!wavefront_set->in_mwavefront_gap1->null && max_lo < wavefront_set->in_mwavefront_gap1->lo+1) max_lo = wavefront_set->in_mwavefront_gap1->lo+1;
  if (!wavefront_set->in_mwavefront_gap1->null && min_hi > wavefront_set->in_mwavefront_gap1->hi-1) min_hi = wavefront_set->in_mwavefront_gap1->hi-1;
  if (distance_metric==gap_lineal) {
    *lo = max_lo; *hi = min_hi; return;
  }
  // Gap-affine
  if (!wavefront_set->in_i1wavefront_ext->null && max_lo < wavefront_set->in_i1wavefront_ext->lo+1) max_lo = wavefront_set->in_i1wavefront_ext->lo+1;
  if (!wavefront_set->in_d1wavefront_ext->null && max_lo < wavefront_set->in_d1wavefront_ext->lo-1) max_lo = wavefront_set->in_d1wavefront_ext->lo-1;
  if (!wavefront_set->in_i1wavefront_ext->null && min_hi > wavefront_set->in_i1wavefront_ext->hi+1) min_hi = wavefront_set->in_i1wavefront_ext->hi+1;
  if (!wavefront_set->in_d1wavefront_ext->null && min_hi > wavefront_set->in_d1wavefront_ext->hi-1) min_hi = wavefront_set->in_d1wavefront_ext->hi-1;
  if (distance_metric==gap_affine) {
    *lo = max_lo; *hi = min_hi; return;
  }
  // Gap-affine-2p
  if (!wavefront_set->in_mwavefront_gap2->null && max_lo < wavefront_set->in_mwavefront_gap2->lo+1) max_lo = wavefront_set->in_mwavefront_gap2->lo+1;
  if (!wavefront_set->in_i2wavefront_ext->null && max_lo < wavefront_set->in_i2wavefront_ext->lo+1) max_lo = wavefront_set->in_i2wavefront_ext->lo+1;
  if (!wavefront_set->in_d2wavefront_ext->null && max_lo < wavefront_set->in_d2wavefront_ext->lo-1) max_lo = wavefront_set->in_d2wavefront_ext->lo-1;
  if (!wavefront_set->in_mwavefront_gap2->null && min_hi > wavefront_set->in_mwavefront_gap2->hi-1) min_hi = wavefront_set->in_mwavefront_gap2->hi-1;
  if (!wavefront_set->in_i2wavefront_ext->null && min_hi > wavefront_set->in_i2wavefront_ext->hi+1) min_hi = wavefront_set->in_i2wavefront_ext->hi+1;
  if (!wavefront_set->in_d2wavefront_ext->null && min_hi > wavefront_set->in_d2wavefront_ext->hi-1) min_hi = wavefront_set->in_d2wavefront_ext->hi-1;
  *lo = max_lo;
  *hi = min_hi;
}
/*
 * Input wavefronts (fetch)
 */
wavefront_t* wavefront_components_get_mwavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->mwavefronts[score] == NULL ||
          wf_components->mwavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->mwavefronts[score];
}
wavefront_t* wavefront_components_get_i1wavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->i1wavefronts[score] == NULL ||
          wf_components->i1wavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->i1wavefronts[score];
}
wavefront_t* wavefront_components_get_i2wavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->i2wavefronts[score] == NULL ||
          wf_components->i2wavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->i2wavefronts[score];
}
wavefront_t* wavefront_components_get_d1wavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->d1wavefronts[score] == NULL ||
          wf_components->d1wavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->d1wavefronts[score];
}
wavefront_t* wavefront_components_get_d2wavefront(
    wavefront_components_t* const wf_components,
    const int score) {
  return (score < 0 ||
          wf_components->d2wavefronts[score] == NULL ||
          wf_components->d2wavefronts[score]->null) ?
      wf_components->wavefront_null : wf_components->d2wavefronts[score];
}
void wavefront_aligner_fetch_input(
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
  wavefront_set->in_mwavefront_sub = wavefront_components_get_mwavefront(wf_components,mismatch);
  wavefront_set->in_mwavefront_gap1 = wavefront_components_get_mwavefront(wf_components,gap_open1);
  if (distance_metric==gap_lineal) return;
  wavefront_set->in_i1wavefront_ext = wavefront_components_get_i1wavefront(wf_components,gap_extend1);
  wavefront_set->in_d1wavefront_ext = wavefront_components_get_d1wavefront(wf_components,gap_extend1);
  if (distance_metric==gap_affine) return;
  wavefront_set->in_mwavefront_gap2 = wavefront_components_get_mwavefront(wf_components,gap_open2);
  wavefront_set->in_i2wavefront_ext = wavefront_components_get_i2wavefront(wf_components,gap_extend2);
  wavefront_set->in_d2wavefront_ext = wavefront_components_get_d2wavefront(wf_components,gap_extend2);
}
/*
 * Output wavefronts (allocate)
 */
void wavefront_aligner_free_output(
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
void wavefront_aligner_allocate_output_null(
    wavefront_aligner_t* const wf_aligner,
    int score) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Modular wavefront
  if (wf_components->memory_modular) {
    score = score % wf_components->max_score_scope;
    wavefront_aligner_free_output(wf_aligner,score);
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
void wavefront_aligner_allocate_output(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    int score,
    const int lo,
    const int hi) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_slab_t* const wavefront_slab = wf_aligner->wavefront_slab;
  // Resize null/victim wavefronts
  wavefront_components_resize_null__victim(wf_components,lo,hi);
  // Modular wavefront
  if (wf_components->memory_modular) {
    score = score % wf_components->max_score_scope;
    wavefront_aligner_free_output(wf_aligner,score);
  }
  // Allocate M-Wavefront
  wavefront_set->out_mwavefront = wavefront_slab_allocate(wavefront_slab,lo,hi);
  wf_components->mwavefronts[score] = wavefront_set->out_mwavefront;
  if (distance_metric==gap_lineal) return;
  // Allocate I1-Wavefront
  if (!wavefront_set->in_mwavefront_gap1->null || !wavefront_set->in_i1wavefront_ext->null) {
    wavefront_set->out_i1wavefront = wavefront_slab_allocate(wavefront_slab,lo,hi);
    wf_components->i1wavefronts[score] = wavefront_set->out_i1wavefront;
  } else {
    wavefront_set->out_i1wavefront = wf_components->wavefront_victim;
    wf_components->i1wavefronts[score] = NULL;
  }
  // Allocate D1-Wavefront
  if (!wavefront_set->in_mwavefront_gap1->null || !wavefront_set->in_d1wavefront_ext->null) {
    wavefront_set->out_d1wavefront = wavefront_slab_allocate(wavefront_slab,lo,hi);
    wf_components->d1wavefronts[score] = wavefront_set->out_d1wavefront;
  } else {
    wavefront_set->out_d1wavefront = wf_components->wavefront_victim;
    wf_components->d1wavefronts[score] = NULL;
  }
  if (distance_metric==gap_affine) return;
  // Allocate I2-Wavefront
  if (!wavefront_set->in_mwavefront_gap2->null || !wavefront_set->in_i2wavefront_ext->null) {
    wavefront_set->out_i2wavefront = wavefront_slab_allocate(wavefront_slab,lo,hi);
    wf_components->i2wavefronts[score] = wavefront_set->out_i2wavefront;
  } else {
    wavefront_set->out_i2wavefront = wf_components->wavefront_victim;
    wf_components->i2wavefronts[score] = NULL;
  }
  // Allocate D2-Wavefront
  if (!wavefront_set->in_mwavefront_gap2->null || !wavefront_set->in_d2wavefront_ext->null) {
    wavefront_set->out_d2wavefront = wavefront_slab_allocate(wavefront_slab,lo,hi);
    wf_components->d2wavefronts[score] = wavefront_set->out_d2wavefront;
  } else {
    wavefront_set->out_d2wavefront = wf_components->wavefront_victim;
    wf_components->d2wavefronts[score] = NULL;
  }
}
/*
 * Trim wavefronts ends
 */
void wavefront_aligner_trim_ends_wf(
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
}
void wavefront_aligner_trim_ends(
    wavefront_aligner_t* const wf_aligner,
    int score) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Modular wavefront
  if (wf_components->memory_modular) {
    score = score % wf_components->max_score_scope;
  }
  // Free
  if (wf_components->mwavefronts[score]) wavefront_aligner_trim_ends_wf(wf_aligner,wf_components->mwavefronts[score]);
  if (distance_metric==gap_lineal) return;
  if (wf_components->i1wavefronts[score]) wavefront_aligner_trim_ends_wf(wf_aligner,wf_components->i1wavefronts[score]);
  if (wf_components->d1wavefronts[score]) wavefront_aligner_trim_ends_wf(wf_aligner,wf_components->d1wavefronts[score]);
  if (distance_metric==gap_affine) return;
  if (wf_components->i2wavefronts[score]) wavefront_aligner_trim_ends_wf(wf_aligner,wf_components->i2wavefronts[score]);
  if (wf_components->d2wavefronts[score]) wavefront_aligner_trim_ends_wf(wf_aligner,wf_components->d2wavefronts[score]);
}



