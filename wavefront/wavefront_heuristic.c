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
 * DESCRIPTION: Support functions for wavefront heuristic strategies
 */

#include "wavefront_heuristic.h"
#include "wavefront_aligner.h"

/*
 * Setup
 */
void wavefront_heuristic_set_none(
    wavefront_heuristic_t* const wf_heuristic) {
  wf_heuristic->strategy = wf_heuristic_none;
}
void wavefront_heuristic_set_banded_static(
    wavefront_heuristic_t* const wf_heuristic,
    const int band_min_k,
    const int band_max_k) {
  wf_heuristic->strategy = wf_heuristic_banded_static;
  wf_heuristic->min_k = band_min_k;
  wf_heuristic->max_k = band_max_k;
}
void wavefront_heuristic_set_banded_adaptive(
    wavefront_heuristic_t* const wf_heuristic,
    const int band_min_k,
    const int band_max_k,
    const int steps_between_cutoffs) {
  wf_heuristic->strategy = wf_heuristic_banded_adaptive;
  wf_heuristic->min_k = band_min_k;
  wf_heuristic->max_k = band_max_k;
  wf_heuristic->steps_between_cutoffs = steps_between_cutoffs;
  // Internals
  wf_heuristic->steps_wait = steps_between_cutoffs;
}
void wavefront_heuristic_set_wfadaptive(
    wavefront_heuristic_t* const wf_heuristic,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int steps_between_cutoffs) {
  wf_heuristic->strategy = wf_heuristic_wfadaptive;
  wf_heuristic->min_wavefront_length = min_wavefront_length;
  wf_heuristic->max_distance_threshold = max_distance_threshold;
  wf_heuristic->steps_between_cutoffs = steps_between_cutoffs;
  // Internals
  wf_heuristic->steps_wait = steps_between_cutoffs;
}
void wavefront_heuristic_set_xdrop(
    wavefront_heuristic_t* const wf_heuristic,
    const int xdrop,
    const int steps_between_cutoffs) {
  wf_heuristic->strategy = wf_heuristic_xdrop;
  wf_heuristic->xdrop = xdrop;
  wf_heuristic->steps_between_cutoffs = steps_between_cutoffs;
  // Internals
  wf_heuristic->steps_wait = steps_between_cutoffs;
  wf_heuristic->max_sw_score = 0;
  wf_heuristic->max_sw_score_offset = WAVEFRONT_OFFSET_NULL;
  wf_heuristic->max_sw_score_k = DPMATRIX_DIAGONAL_NULL;
}
void wavefront_heuristic_set_zdrop(
    wavefront_heuristic_t* const wf_heuristic,
    const int zdrop,
    const int steps_between_cutoffs) {
  wf_heuristic->strategy = wf_heuristic_zdrop;
  wf_heuristic->zdrop = zdrop;
  wf_heuristic->steps_between_cutoffs = steps_between_cutoffs;
  // Internals
  wf_heuristic->steps_wait = steps_between_cutoffs;
  wf_heuristic->max_sw_score = 0;
  wf_heuristic->max_sw_score_offset = WAVEFRONT_OFFSET_NULL;
  wf_heuristic->max_sw_score_k = DPMATRIX_DIAGONAL_NULL;
}
void wavefront_heuristic_clear(
    wavefront_heuristic_t* const wf_heuristic) {
  // Internals
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
  wf_heuristic->max_sw_score = 0;
  wf_heuristic->max_sw_score_offset = WAVEFRONT_OFFSET_NULL;
  wf_heuristic->max_sw_score_k = DPMATRIX_DIAGONAL_NULL;
}
/*
 * Utils
 */
int wf_compute_distance_end2end(
    const wf_offset_t offset,
    const int k,
    const int pattern_length,
    const int text_length) {
  const int left_v = pattern_length - WAVEFRONT_V(k,offset);
  const int left_h = text_length - WAVEFRONT_H(k,offset);
  return (offset >= 0) ? MAX(left_v,left_h) : -WAVEFRONT_OFFSET_NULL;
}
int wf_compute_distance_endsfree(
    const wf_offset_t offset,
    const int k,
    const int pattern_length,
    const int text_length,
    const int pattern_end_free,
    const int text_end_free) {
  const int left_v = pattern_length - WAVEFRONT_V(k,offset);
  const int left_h = text_length - WAVEFRONT_H(k,offset);
  const int left_v_endsfree = left_v - pattern_end_free;
  const int left_h_endsfree = left_h - text_end_free;
  const int dist_up = MAX(left_h,left_v_endsfree);
  const int dist_down = MAX(left_v,left_h_endsfree);
  return (offset >= 0) ? MIN(dist_up,dist_down) : -WAVEFRONT_OFFSET_NULL;
}
int wf_compute_sw_score(
    const int wf_score,
    const wf_offset_t offset,
    const int k) {
  const int v = WAVEFRONT_V(k,offset);
  const int h = WAVEFRONT_H(k,offset);
  return (offset >= 0) ? MIN(v,h) - wf_score : WAVEFRONT_OFFSET_NULL;
}
int wf_compute_sw_score_single_gap(
    const int gap_extension,
    const wf_offset_t wf1_offset,
    const int wf1_k,
    const wf_offset_t wf2_offset,
    const int wf2_k) {
  const int diff_h = WAVEFRONT_H(wf2_k,wf2_offset) - WAVEFRONT_H(wf1_k,wf1_offset);
  const int diff_v = WAVEFRONT_V(wf2_k,wf2_offset) - WAVEFRONT_V(wf1_k,wf1_offset);
  const int abs_diff = (diff_h >= diff_v) ? diff_h-diff_v : diff_v-diff_h;
  return abs_diff * gap_extension;
}
void wavefront_heuristic_cutoff_equate(
    wavefront_t* const wavefront_dst,
    wavefront_t* const wavefront_src) {
  if (wavefront_dst != NULL) {
    if (wavefront_src->lo > wavefront_dst->lo) wavefront_dst->lo = wavefront_src->lo;
    if (wavefront_src->hi < wavefront_dst->hi) wavefront_dst->hi = wavefront_src->hi;
    if (wavefront_dst->lo > wavefront_dst->hi) wavefront_dst->null = true;
    // Save min/max WF initialized
    wavefront_dst->wf_elements_init_min = wavefront_dst->lo;
    wavefront_dst->wf_elements_init_max = wavefront_dst->hi;
  }
}
/*
 * Cut-offs Banded
 */
void wavefront_cufoff_banded_static(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront) {
  // Parameters
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  // Check wavefront limits
  if (wavefront->lo < wf_heuristic->min_k) wavefront->lo = wf_heuristic->min_k;
  if (wavefront->hi > wf_heuristic->max_k) wavefront->hi = wf_heuristic->max_k;
}
void wavefront_cufoff_banded_adaptive(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  // Check steps
  if (wf_heuristic->steps_wait > 0) return;
  // Check wavefront length
  const int lo = wavefront->lo;
  const int hi = wavefront->hi;
  const int wf_length = hi - lo + 1;
  if (wf_length < 4) return; // We cannot do anything here
  // Adjust the band
  const wf_offset_t* const offsets = wavefront->offsets;
  const int max_wf_length = wf_heuristic->max_k - wf_heuristic->min_k + 1;
  if (wf_length > max_wf_length) {
    // Sample wavefront
    const int leeway = (wf_length - max_wf_length) / 2;
    const int quarter = wf_length / 4;
    const int dist_p0 = wf_compute_distance_end2end(
        offsets[lo],lo,pattern_length,text_length);
    const int dist_p1 = wf_compute_distance_end2end(
        offsets[lo+quarter],lo+quarter,pattern_length,text_length);
    const int dist_p2 = wf_compute_distance_end2end(
        offsets[lo+2*quarter],lo+2*quarter,pattern_length,text_length);
    const int dist_p3 = wf_compute_distance_end2end(
        offsets[hi],hi,pattern_length,text_length);
    // Heuristically decide where to place the band
    int new_lo = lo;
    if (dist_p0 > dist_p3) new_lo += leeway;
    if (dist_p1 > dist_p2) new_lo += leeway;
    // Set wavefront limits
    wavefront->lo = new_lo;
    if (wavefront->lo < lo) wavefront->lo = lo;
    wavefront->hi = new_lo + max_wf_length - 1;
    if (wavefront->hi > hi) wavefront->hi = hi;
  }
  // Set wait steps (don't repeat this heuristic often)
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
}
/*
 * Cut-off Wavefront Adaptive
 */
int wavefront_compute_distance_end2end(
    wavefront_t* const wavefront,
    const int pattern_length,
    const int text_length,
    wf_offset_t* const distances) {
  // Compute min-distance
  const wf_offset_t* const offsets = wavefront->offsets;
  int k, min_distance = MAX(pattern_length,text_length);
  PRAGMA_LOOP_VECTORIZE
  for (k=wavefront->lo;k<=wavefront->hi;++k) {
    const int distance = wf_compute_distance_end2end(
        offsets[k],k,pattern_length,text_length);
    distances[k] = distance;
    min_distance = MIN(min_distance,distance);
  }
  return min_distance;
}
int wavefront_compute_distance_endsfree(
    wavefront_t* const wavefront,
    const int pattern_length,
    const int text_length,
    const int pattern_end_free,
    const int text_end_free,
    wf_offset_t* const distances) {
  // Compute min-distance
  const wf_offset_t* const offsets = wavefront->offsets;
  int k, min_distance = MAX(pattern_length,text_length);
  PRAGMA_LOOP_VECTORIZE
  for (k=wavefront->lo;k<=wavefront->hi;++k) {
    const int distance = wf_compute_distance_endsfree(
        offsets[k],k,pattern_length,text_length,
        pattern_end_free,text_end_free);
    distances[k] = distance;
    min_distance = MIN(min_distance,distance);
  }
  return min_distance;
}
void wavefront_cufoff_wfadaptive_reduce(
    wavefront_t* const wavefront,
    const wf_offset_t* const distances,
    const int min_distance,
    const int max_distance_threshold,
    const int min_k,
    const int max_k) {
  int k;
  // Reduce from bottom
  const int top_limit = MIN(max_k,wavefront->hi); // Preserve target-diagonals
  int lo_reduced = wavefront->lo;
  for (k=wavefront->lo;k<top_limit;++k) {
    if (distances[k] - min_distance  <= max_distance_threshold) break;
    ++lo_reduced;
  }
  wavefront->lo = lo_reduced;
  // Reduce from top
  const int botton_limit = MAX(min_k,wavefront->lo); // Preserve target-diagonals
  int hi_reduced = wavefront->hi;
  for (k=wavefront->hi;k>botton_limit;--k) {
    if (distances[k] - min_distance <= max_distance_threshold) break;
    --hi_reduced;
  }
  wavefront->hi = hi_reduced;
}
void wavefront_cufoff_wfadaptive(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const int min_wavefront_length = wf_aligner->heuristic.min_wavefront_length;
  const int max_distance_threshold = wf_aligner->heuristic.max_distance_threshold;
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  // Check steps
  if (wf_heuristic->steps_wait > 0) return;
  // Check minimum wavefront length
  const int base_hi = wavefront->hi;
  const int base_lo = wavefront->lo;
  if ((base_hi - base_lo + 1) < min_wavefront_length) return;
  // Use victim as temporal buffer
  wavefront_components_resize_null__victim(&wf_aligner->wf_components,base_lo-1,base_hi+1);
  wf_offset_t* const distances = wf_aligner->wf_components.wavefront_victim->offsets;
  // Compute distance & cut-off
//  const int pattern_end_free = wf_aligner->alignment_form.pattern_end_free;
//  const int text_end_free = wf_aligner->alignment_form.text_end_free;
//  if ((wf_aligner->alignment_form.span == alignment_end2end) ||
//      (pattern_end_free==0 && text_end_free==0)) {
    const int min_distance = wavefront_compute_distance_end2end(
        wavefront,pattern_length,text_length,distances);
    // Cut-off wavefront
    const int alignment_k = DPMATRIX_DIAGONAL(text_length,pattern_length);
    wavefront_cufoff_wfadaptive_reduce(
        wavefront,distances,min_distance,max_distance_threshold,
        alignment_k,alignment_k);
//  }
//  else {
//    const int min_distance = wavefront_compute_distance_endsfree(
//        wavefront,pattern_length,text_length,
//        pattern_end_free,text_end_free,distances);
//    // Cut-off wavefront
//    const int alignment_k = DPMATRIX_DIAGONAL(text_length,pattern_length);
//    wavefront_cufoff_wfadaptive_reduce(
//        wavefront,distances,min_distance,max_distance_threshold,
//        alignment_k-text_end_free,alignment_k+pattern_end_free);
//  }
  // Set wait steps (don't repeat this heuristic often)
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
}
/*
 * Drops
 */
void wavefront_compute_sw_scores(
    wavefront_t* const wavefront,
    const int score,
    wf_offset_t* const sw_scores,
    wf_offset_t* const max_sw_score,
    wf_offset_t* const max_sw_score_k) {
  // Compute min-distance
  const wf_offset_t* const offsets = wavefront->offsets;
  int k, score_max = -score, score_max_k = 0;
  PRAGMA_LOOP_VECTORIZE
  for (k=wavefront->lo;k<=wavefront->hi;++k) {
    const int sw_score = wf_compute_sw_score(score,offsets[k],k);
    sw_scores[k] = sw_score;
    if (score_max < sw_score) {
      score_max = sw_score;
      score_max_k = k;
    }
  }
  *max_sw_score = score_max;
  *max_sw_score_k = score_max_k;
}
void wavefront_cufoff_xdrop(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront,
    const int score) {
  // Parameters
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  const int base_hi = wavefront->hi;
  const int base_lo = wavefront->lo;
  // Check steps
  if (wf_heuristic->steps_wait > 0) return;
  // Use victim as temporal buffer
  wavefront_components_resize_null__victim(&wf_aligner->wf_components,base_lo-1,base_hi+1);
  wf_offset_t* const sw_scores = wf_aligner->wf_components.wavefront_victim->offsets;
  // Compute SW scores (classic scores)
  wf_offset_t current_max_sw_score;
  wf_offset_t current_max_sw_score_k;
  wavefront_compute_sw_scores(
      wavefront,score,sw_scores,
      &current_max_sw_score,&current_max_sw_score_k);
  // Apply X-Drop
  const int xdrop = wf_heuristic->xdrop;
  const int max_sw_score = wf_heuristic->max_sw_score;
  const int max_sw_score_k = wf_heuristic->max_sw_score_k;
  if (max_sw_score_k != DPMATRIX_DIAGONAL_NULL) {
    // Reduce from bottom
    int k, lo_reduced = wavefront->lo;
    for (k=wavefront->lo;k<=wavefront->hi;++k) {
      if ((int)sw_scores[k] >= max_sw_score - xdrop) break;
      ++lo_reduced;
    }
    wavefront->lo = lo_reduced;
    // Reduce from top
    int hi_reduced = wavefront->hi;
    for (k=wavefront->hi;k>=wavefront->lo;--k) {
      if ((int)sw_scores[k] >= max_sw_score - xdrop) break;
      --hi_reduced;
    }
    wavefront->hi = hi_reduced;
  }
  // Update maximum score observed
  wf_heuristic->max_sw_score = current_max_sw_score;
  wf_heuristic->max_sw_score_k = current_max_sw_score_k;
  // Set wait steps (don't repeat this heuristic often)
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
}
void wavefront_cufoff_zdrop(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront,
    const int score) {
  // Parameters
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  const int base_hi = wavefront->hi;
  const int base_lo = wavefront->lo;
  // Check steps
  if (wf_heuristic->steps_wait > 0) return;
  // Use victim as temporal buffer
  wavefront_components_resize_null__victim(&wf_aligner->wf_components,base_lo-1,base_hi+1);
  wf_offset_t* const sw_scores = wf_aligner->wf_components.wavefront_victim->offsets;
  // Compute SW scores (classic scores)
  wf_offset_t current_max_sw_score;
  wf_offset_t current_max_sw_score_k;
  wavefront_compute_sw_scores(
      wavefront,score,sw_scores,
      &current_max_sw_score,&current_max_sw_score_k);
  // Apply X-Drop
  const wf_offset_t* const offsets = wavefront->offsets;
  wavefronts_penalties_t* const penalties = &wf_aligner->penalties;
  const int gap_extension = (penalties->gap_extension1 > 0) ? penalties->gap_extension1 : 1;
  const int zdrop = wf_heuristic->zdrop;
  const int max_sw_score = wf_heuristic->max_sw_score;
  const int max_sw_score_k = wf_heuristic->max_sw_score_k;
  const int max_sw_score_offset = wf_heuristic->max_sw_score_offset;
  if (max_sw_score_k != DPMATRIX_DIAGONAL_NULL) {
    // Reduce from bottom
    int k, lo_reduced = wavefront->lo;
    for (k=wavefront->lo;k<=wavefront->hi;++k) {
      if (offsets[k] < 0) {
        ++lo_reduced;
        continue;
      }
      const int single_gap = wf_compute_sw_score_single_gap(
          gap_extension,max_sw_score_offset,max_sw_score_k,offsets[k],k);
      if ((int)sw_scores[k] > max_sw_score - (zdrop + single_gap)) break;
      ++lo_reduced;
    }
    wavefront->lo = lo_reduced;
    // Reduce from top
    int hi_reduced = wavefront->hi;
    for (k=wavefront->hi;k>=wavefront->lo;--k) {
      if (offsets[k] < 0) {
        --hi_reduced;
        continue;
      }
      const int single_gap = wf_compute_sw_score_single_gap(
          gap_extension,max_sw_score_offset,max_sw_score_k,offsets[k],k);
      if ((int)sw_scores[k] > max_sw_score - (zdrop + single_gap)) break;
      --hi_reduced;
    }
    wavefront->hi = hi_reduced;
    // Update maximum score observed
    if (current_max_sw_score > wf_heuristic->max_sw_score) {
      wf_heuristic->max_sw_score = current_max_sw_score;
      wf_heuristic->max_sw_score_k = current_max_sw_score_k;
      wf_heuristic->max_sw_score_offset = max_sw_score_offset;
    }
  } else {
    // Update maximum score observed
    wf_heuristic->max_sw_score = current_max_sw_score;
    wf_heuristic->max_sw_score_k = current_max_sw_score_k;
    wf_heuristic->max_sw_score_offset = max_sw_score_offset;
  }
  // Set wait steps (don't repeat this heuristic often)
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
}
/*
 * Cut-offs dispatcher
 */
bool wavefront_heuristic_cufoff(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int score_mod) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_components->mwavefronts[score_mod];
  if (mwavefront==NULL || mwavefront->lo > mwavefront->hi) return false; // Not dropped
  // Cut-off m-wavefront
  --(wf_heuristic->steps_wait);
  if (wf_heuristic->strategy == wf_heuristic_banded_static) {
    wavefront_cufoff_banded_static(wf_aligner,mwavefront);
  }
  if (wf_heuristic->strategy == wf_heuristic_banded_adaptive) {
    wavefront_cufoff_banded_adaptive(wf_aligner,mwavefront);
  }
  if (wf_heuristic->strategy == wf_heuristic_wfadaptive) {
    wavefront_cufoff_wfadaptive(wf_aligner,mwavefront);
  }
  if (wf_heuristic->strategy == wf_heuristic_xdrop) {
    wavefront_cufoff_xdrop(wf_aligner,mwavefront,score);
  }
  if (wf_heuristic->strategy == wf_heuristic_zdrop) {
    wavefront_cufoff_zdrop(wf_aligner,mwavefront,score);
  }
  // Check wavefront length
  if (mwavefront->lo > mwavefront->hi) return true; // Dropped alignment
  // Save min/max WF initialized
  mwavefront->wf_elements_init_min = mwavefront->lo;
  mwavefront->wf_elements_init_max = mwavefront->hi;
  // Equate other wavefronts
  if (distance_metric <= gap_linear) return false; // Not dropped
  // Cut-off the other wavefronts (same dimensions as M)
  wavefront_t* const i1wavefront = wf_components->i1wavefronts[score_mod];
  wavefront_t* const d1wavefront = wf_components->d1wavefronts[score_mod];
  wavefront_heuristic_cutoff_equate(i1wavefront,mwavefront);
  wavefront_heuristic_cutoff_equate(d1wavefront,mwavefront);
  if (distance_metric == gap_affine) return false; // Not dropped
  wavefront_t* const i2wavefront = wf_components->i2wavefronts[score_mod];
  wavefront_t* const d2wavefront = wf_components->d2wavefronts[score_mod];
  wavefront_heuristic_cutoff_equate(i2wavefront,mwavefront);
  wavefront_heuristic_cutoff_equate(d2wavefront,mwavefront);
  // Return
  return false; // Not dropped
  // DEBUG
  //  if (wf_aligner->system.verbose) {
  //    const int wf_length_base = hi_base-lo_base+1;
  //    const int wf_length_reduced = mwavefront->hi-mwavefront->lo+1;
  //    fprintf(stderr,"[WFA::Heuristic] Heuristic from %d to %d offsets (%2.2f%%)\n",
  //        wf_length_base,wf_length_reduced,100.0f*(float)wf_length_reduced/(float)wf_length_base);
  //  }
}
/*
 * Display
 */
void wavefront_heuristic_print(
    FILE* const stream,
    wavefront_heuristic_t* const wf_heuristic) {
  // Select heuristic strategy
  if (wf_heuristic->strategy == wf_heuristic_none) {
    fprintf(stream,"(none)");
  } else if (wf_heuristic->strategy == wf_heuristic_banded_static) {
    fprintf(stream,"(banded-static,%d,%d)",
        wf_heuristic->min_k,
        wf_heuristic->max_k);
  } else if (wf_heuristic->strategy == wf_heuristic_banded_adaptive) {
    fprintf(stream,"(banded-adapt,%d,%d,%d)",
        wf_heuristic->min_k,
        wf_heuristic->max_k,
        wf_heuristic->steps_between_cutoffs);
  } else if (wf_heuristic->strategy == wf_heuristic_wfadaptive) {
    fprintf(stream,"(wf-adapt,%d,%d,%d)",
        wf_heuristic->min_wavefront_length,
        wf_heuristic->max_distance_threshold,
        wf_heuristic->steps_between_cutoffs);
  } else if (wf_heuristic->strategy == wf_heuristic_xdrop) {
    fprintf(stream,"(xdrop,%d,%d)",
        wf_heuristic->xdrop,
        wf_heuristic->steps_between_cutoffs);
  } else if (wf_heuristic->strategy == wf_heuristic_zdrop) {
    fprintf(stream,"(zdrop,%d,%d)",
        wf_heuristic->zdrop,
        wf_heuristic->steps_between_cutoffs);
  }
}
