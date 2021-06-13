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
 * DESCRIPTION: WaveFront aligner components
 */

#include "wavefront_components.h"

/*
 * Setup
 */
void wavefront_components_allocate(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric) {
  // Parameters
  mm_allocator_t* const mm_allocator = wf_aligner->mm_allocator;
  // Wavefronts components
  const int num_wavefronts = wf_aligner->num_wavefronts;
  const bool init_wf = wf_aligner->memory_modular;
  wf_aligner->mwavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
  if (distance_metric==edit || distance_metric==gap_lineal) {
    wf_aligner->i1wavefronts = NULL;
    wf_aligner->d1wavefronts = NULL;
    wf_aligner->i2wavefronts = NULL;
    wf_aligner->d2wavefronts = NULL;
    return;
  }
  wf_aligner->i1wavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
  wf_aligner->d1wavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
  if (distance_metric==gap_affine) {
    wf_aligner->i2wavefronts = NULL;
    wf_aligner->d2wavefronts = NULL;
    return;
  }
  wf_aligner->i2wavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
  wf_aligner->d2wavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
}
void wavefront_components_clear(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric) {
  // Parameters
  const int num_wavefronts = wf_aligner->num_wavefronts;
  // Wavefronts components
  memset(wf_aligner->mwavefronts,0,num_wavefronts*sizeof(wavefront_t*));
  if (distance_metric==edit || distance_metric==gap_lineal) return;
  memset(wf_aligner->i1wavefronts,0,num_wavefronts*sizeof(wavefront_t*));
  memset(wf_aligner->d1wavefronts,0,num_wavefronts*sizeof(wavefront_t*));
  if (distance_metric==gap_affine) return;
  memset(wf_aligner->i2wavefronts,0,num_wavefronts*sizeof(wavefront_t*));
  memset(wf_aligner->d2wavefronts,0,num_wavefronts*sizeof(wavefront_t*));
}
void wavefront_components_free(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric) {
  // Parameters
  mm_allocator_t* const mm_allocator = wf_aligner->mm_allocator;
  // Wavefronts components
  mm_allocator_free(mm_allocator,wf_aligner->mwavefronts);
  if (distance_metric==edit || distance_metric==gap_lineal) return;
  mm_allocator_free(mm_allocator,wf_aligner->i1wavefronts);
  mm_allocator_free(mm_allocator,wf_aligner->d1wavefronts);
  if (distance_metric==gap_affine) return;
  mm_allocator_free(mm_allocator,wf_aligner->i2wavefronts);
  mm_allocator_free(mm_allocator,wf_aligner->d2wavefronts);
}
/*
 * Compute dimensions
 */
void wavefront_components_dimensions_edit(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  // Dimensions
  if (wf_aligner->memory_modular) {
    *max_score_scope = 2;
    *num_wavefronts = 2;
  } else {
    *max_score_scope = -1;
    *num_wavefronts = MAX(pattern_length,text_length);
  }
}
void wavefront_components_dimensions_lineal(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  // Parameters
  wavefronts_penalties_t* const penalties = &wf_aligner->penalties;
  // Dimensions
  if (wf_aligner->memory_modular) {
    *max_score_scope = MAX(penalties->mismatch,penalties->gap_opening1) + 1;
    *num_wavefronts = *max_score_scope;
  } else {
    const int abs_seq_diff = ABS(pattern_length-text_length);
    const int max_score_misms = MIN(pattern_length,text_length) * penalties->mismatch;
    const int max_score_indel = penalties->gap_opening1 * abs_seq_diff;
    *max_score_scope = -1;
    *num_wavefronts = max_score_misms + max_score_indel;
  }
}
void wavefront_components_dimensions_affine(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  // Parameters
  wavefronts_penalties_t* const penalties = &wf_aligner->penalties;
  // Dimensions
  if (wf_aligner->memory_modular) {
    const int max_score_scope_indel = penalties->gap_opening1+penalties->gap_extension1;
    *max_score_scope = MAX(max_score_scope_indel,penalties->mismatch) + 1;
    *num_wavefronts = *max_score_scope;
  } else {
    const int abs_seq_diff = ABS(pattern_length-text_length);
    const int max_score_misms = MIN(pattern_length,text_length) * penalties->mismatch;
    const int max_score_indel = penalties->gap_opening1 + abs_seq_diff * penalties->gap_extension1;
    *max_score_scope = -1;
    *num_wavefronts = max_score_misms + max_score_indel;
  }
}
void wavefront_components_dimensions_affine2p(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  // Parameters
  wavefronts_penalties_t* const penalties = &wf_aligner->penalties;
  // Dimensions
  if (wf_aligner->memory_modular) {
    const int max_score_scope_indel =
        MAX(penalties->gap_opening1+penalties->gap_extension1,
            penalties->gap_opening2+penalties->gap_extension2);
    *max_score_scope = MAX(max_score_scope_indel,penalties->mismatch) + 1;
    *num_wavefronts = *max_score_scope;
  } else {
    const int abs_seq_diff = ABS(pattern_length-text_length);
    const int max_score_misms = MIN(pattern_length,text_length) * penalties->mismatch;
    const int max_score_indel1 = penalties->gap_opening1 + abs_seq_diff * penalties->gap_extension1;
    const int max_score_indel2 = penalties->gap_opening2 + abs_seq_diff * penalties->gap_extension2;
    const int max_score_indel = MIN(max_score_indel1,max_score_indel2);
    *max_score_scope = -1;
    *num_wavefronts = max_score_misms + max_score_indel;
  }
}
void wavefront_components_dimensions(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric,
    const int pattern_length,
    const int text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  switch (distance_metric) {
    case edit:
      wavefront_components_dimensions_edit(wf_aligner,
          pattern_length,text_length,
          max_score_scope,
          num_wavefronts);
      break;
    case gap_lineal:
      wavefront_components_dimensions_lineal(wf_aligner,
          pattern_length,text_length,
          max_score_scope,
          num_wavefronts);
      break;
    case gap_affine:
      wavefront_components_dimensions_affine(wf_aligner,
          pattern_length,text_length,
          max_score_scope,
          num_wavefronts);
      break;
    case gap_affine_2p:
      wavefront_components_dimensions_affine2p(wf_aligner,
          pattern_length,text_length,
          max_score_scope,
          num_wavefronts);
      break;
  }
}

