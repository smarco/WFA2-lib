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
 * DESCRIPTION: Support functions for wavefront reduction strategies
 */

#include "wavefront_reduction.h"
#include "wavefront_aligner.h"

/*
 * Setup
 */
void wavefront_reduction_set_none(
    wavefront_reduction_t* const wavefront_reduction) {
  wavefront_reduction->reduction_strategy = wavefront_reduction_none;
}
void wavefront_reduction_set_dynamic(
    wavefront_reduction_t* const wavefront_reduction,
    const int min_wavefront_length,
    const int max_distance_threshold) {
  wavefront_reduction->reduction_strategy = wavefront_reduction_dynamic;
  wavefront_reduction->min_wavefront_length = min_wavefront_length;
  wavefront_reduction->max_distance_threshold = max_distance_threshold;
}
/*
 * Distances
 */
int wavefront_distance_to_target_global(
    const int pattern_length,
    const int text_length,
    const wf_offset_t offset,
    const int k) {
  const int v = WAVEFRONT_V(k,offset);
  const int h = WAVEFRONT_H(k,offset);
  const int left_v = pattern_length - v;
  const int left_h = text_length - h;
  //return left_v + left_h;
  return MAX(left_v,left_h);
}
int wavefront_distance_to_target_semiglobal(
    const int pattern_length,
    const int text_length,
    const wf_offset_t offset,
    const int k) {
  return pattern_length - WAVEFRONT_V(k,offset);
}
/*
 * Reduce wavefront
 */
void wavefront_reduce_equate(
    wavefront_t* const wavefront_dst,
    wavefront_t* const wavefront_src) {
  if (wavefront_dst!=NULL) {
    if (wavefront_src->lo > wavefront_dst->lo) wavefront_dst->lo = wavefront_src->lo;
    if (wavefront_src->hi < wavefront_dst->hi) wavefront_dst->hi = wavefront_src->hi;
    if (wavefront_dst->lo > wavefront_dst->hi) wavefront_dst->null = true;
  }
}
void wavefront_reduce_mwavefront_global(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront,
    const int pattern_length,
    const int text_length,
    const int min_distance,
    const int max_distance_threshold,
    const int alignment_k) {
  // Parameters
  const wf_offset_t* const offsets = wavefront->offsets;
  int k;
  // Reduce from bottom
  const int top_limit = MIN(alignment_k-1,wavefront->hi); // Preserve target-diagonal
  int lo_reduced = wavefront->lo;
  for (k=wavefront->lo;k<top_limit;++k) {
    const int distance = wavefront_distance_to_target_global(pattern_length,text_length,offsets[k],k);
    if (distance - min_distance  <= max_distance_threshold) break;
    ++lo_reduced;
  }
  wavefront->lo = lo_reduced;
  // Reduce from top
  const int botton_limit = MAX(alignment_k+1,wavefront->lo); // Preserve target-diagonal
  int hi_reduced = wavefront->hi;
  for (k=wavefront->hi;k>botton_limit;--k) {
    const int distance = wavefront_distance_to_target_global(pattern_length,text_length,offsets[k],k);
    if (distance - min_distance <= max_distance_threshold) break;
    --hi_reduced;
  }
  wavefront->hi = hi_reduced;
  // Check hi/lo range
  if (wavefront->lo > wavefront->hi) {
    wavefront->null = true;
  }
}
void wavefront_reduce(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length,
    const int score) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->distance_metric;
  const int min_wavefront_length = wf_aligner->reduction.min_wavefront_length;
  const int max_distance_threshold = wf_aligner->reduction.max_distance_threshold;
  const int alignment_k = WAVEFRONT_DIAGONAL(text_length,pattern_length);
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->mwavefronts[score];
  if (mwavefront==NULL) return;
  if ((mwavefront->hi - mwavefront->lo + 1) < min_wavefront_length) return;
  // Compute min-distance
  const wf_offset_t* const offsets = mwavefront->offsets;
  int min_distance = MAX(pattern_length,text_length);
  int k;
  for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
    const int distance = wavefront_distance_to_target_global(
        pattern_length,text_length,offsets[k],k);
    min_distance = MIN(min_distance,distance);
  }
  // Reduce m-wavefront
  wavefront_reduce_mwavefront_global(
      wf_aligner,mwavefront,
      pattern_length,text_length,
      min_distance,max_distance_threshold,
      alignment_k);
  if (distance_metric <= gap_lineal) return;
  // Reduce the other wavefronts (same dimensions as M-reduced)
  wavefront_t* const i1wavefront = wf_aligner->i1wavefronts[score];
  wavefront_t* const d1wavefront = wf_aligner->d1wavefronts[score];
  wavefront_reduce_equate(i1wavefront,mwavefront);
  wavefront_reduce_equate(d1wavefront,mwavefront);
  if (distance_metric == gap_affine) return;
  wavefront_t* const i2wavefront = wf_aligner->i2wavefronts[score];
  wavefront_t* const d2wavefront = wf_aligner->d2wavefronts[score];
  wavefront_reduce_equate(i2wavefront,mwavefront);
  wavefront_reduce_equate(d2wavefront,mwavefront);
}
