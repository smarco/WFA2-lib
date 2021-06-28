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
 * DESCRIPTION: WaveFront aligner data structure
 */

#include "wavefront_aligner.h"
#include "wavefront_reduction.h"
#include "wavefront_components.h"

/*
 * Configuration
 */
#define WF_NULL_INIT_LO     (-1024)
#define WF_NULL_INIT_HI     ( 1024)
#define WF_NULL_INIT_LENGTH WAVEFRONT_LENGTH(WF_NULL_INIT_LO,WF_NULL_INIT_HI)

/*
 * Penalties
 */
void wavefront_set_penalties(
    wavefront_aligner_t* const wf_aligner,
    wavefront_aligner_attr_t* const attributes) {
  switch (attributes->distance_metric) {
    case edit: break;
    case gap_lineal:
      wavefronts_penalties_set_lineal(
          &wf_aligner->penalties,
          &attributes->lineal_penalties,
          wavefronts_penalties_shifted_penalties);
      break;
    case gap_affine:
      wavefronts_penalties_set_affine(
          &wf_aligner->penalties,
          &attributes->affine_penalties,
          wavefronts_penalties_shifted_penalties);
      break;
    case gap_affine_2p:
      wavefronts_penalties_set_affine2p(
          &wf_aligner->penalties,
          &attributes->affine2p_penalties,
          wavefronts_penalties_shifted_penalties);
      break;
  }
}
/*
 * Setup
 */
wavefront_aligner_t* wavefront_aligner_new(
    const int pattern_length,
    const int text_length,
    wavefront_aligner_attr_t* attributes) {
  // Attributes
  if (attributes == NULL) attributes = &wavefront_aligner_attr_default;
  const bool score_only = (attributes->alignment_scope == compute_score);
  const bool memory_modular = attributes->low_memory || score_only;
  const bool bt_piggyback = attributes->low_memory && !score_only;
  // MM
  mm_allocator_t* mm_allocator = attributes->mm_allocator;
  bool mm_allocator_own = false;
  if (mm_allocator == NULL) {
    mm_allocator = mm_allocator_new(BUFFER_SIZE_4M);
    mm_allocator_own = true;
  }
  // Handler
  wavefront_aligner_t* const wf_aligner = mm_allocator_alloc(mm_allocator,wavefront_aligner_t);
  wf_aligner->mm_allocator = mm_allocator;
  wf_aligner->mm_allocator_own = mm_allocator_own;
  wf_aligner->wavefront_slab = wavefront_slab_new(1000,bt_piggyback,mm_allocator);
  wf_aligner->bt_buffer = (bt_piggyback) ? wf_backtrace_buffer_new(mm_allocator) : NULL;
  // Configuration
  wf_aligner->pattern_length = pattern_length;
  wf_aligner->text_length = text_length;
  wf_aligner->distance_metric = attributes->distance_metric;
  wf_aligner->alignment_scope = attributes->alignment_scope;
  wf_aligner->alignment_form = attributes->alignment_form;
  wf_aligner->memory_modular = memory_modular;
  wf_aligner->bt_piggyback = bt_piggyback;
  wavefront_set_penalties(wf_aligner,attributes); // Set penalties
  // Reduction strategy
  if (attributes->reduction.reduction_strategy == wavefront_reduction_adaptive) {
    wavefront_reduction_set_adaptive(
        &wf_aligner->reduction,
        attributes->reduction.min_wavefront_length,
        attributes->reduction.max_distance_threshold);
  } else { // wavefront_reduction_none
    wavefront_reduction_set_none(&wf_aligner->reduction);
  }
  // Allocate victim wavefront
  wavefront_t* const wavefront_victim = mm_allocator_alloc(mm_allocator,wavefront_t);
  wavefront_allocate(wavefront_victim,WF_NULL_INIT_LENGTH,bt_piggyback,mm_allocator);
  wavefront_init_victim(wavefront_victim,WF_NULL_INIT_LO,WF_NULL_INIT_HI);
  wf_aligner->wavefront_victim = wavefront_victim;
  // Allocate null wavefront
  wavefront_t* const wavefront_null = mm_allocator_alloc(mm_allocator,wavefront_t);
  wavefront_allocate(wavefront_null,WF_NULL_INIT_LENGTH,bt_piggyback,mm_allocator);
  wavefront_init_null(wavefront_null,WF_NULL_INIT_LO,WF_NULL_INIT_HI);
  wf_aligner->wavefront_null = wavefront_null;
  // Allocate wavefronts
  wavefront_components_dimensions(
      wf_aligner,attributes->distance_metric,
      pattern_length,text_length,
      &wf_aligner->max_score_scope,
      &wf_aligner->num_wavefronts); // Compute dimensions
  wavefront_components_allocate(wf_aligner,attributes->distance_metric);
  // CIGAR
  cigar_allocate(&wf_aligner->cigar,2*(pattern_length+text_length),mm_allocator);
  // Limits
  wf_aligner->limit_probe_interval = WF_LIMIT_PROBE_INTERVAL_DEFAULT;
  wf_aligner->max_memory_used = attributes->max_memory_used;
  wf_aligner->max_resident_memory = WF_MAX_MEMORY_RESIDENT_DEFAULT;
  // Return
  return wf_aligner;
}
void wavefront_aligner_reap(
    wavefront_aligner_t* const wf_aligner) {
  if (wf_aligner->bt_buffer) wf_backtrace_buffer_reap(wf_aligner->bt_buffer); // BT-Buffer
  wavefront_slab_reap(wf_aligner->wavefront_slab,wf_slab_reap_all); // Clear Slab
}
void wavefront_aligner_clear(
    wavefront_aligner_t* const wf_aligner) {
  // Clear wavefront components
  if (wf_aligner->memory_modular) wavefront_components_clear(wf_aligner,wf_aligner->distance_metric);
  // Clear CIGAR
  cigar_clear(&wf_aligner->cigar);
  // Clear BT-Buffer and Slab
  if (wf_aligner->bt_buffer) wf_backtrace_buffer_clear(wf_aligner->bt_buffer);
  wavefront_slab_clear(wf_aligner->wavefront_slab);
}
void wavefront_aligner_resize(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length) {
  // Resize wavefront components
  int num_wavefronts;
  wavefront_components_dimensions(
      wf_aligner,wf_aligner->distance_metric,
      pattern_length,text_length,
      &wf_aligner->max_score_scope,&num_wavefronts); // Compute dimensions
  if (num_wavefronts > wf_aligner->num_wavefronts) {
    wf_aligner->num_wavefronts = num_wavefronts;
    wavefront_components_free(wf_aligner,wf_aligner->distance_metric);
    wavefront_components_allocate(wf_aligner,wf_aligner->distance_metric);
  } else {
    if (wf_aligner->memory_modular) wavefront_components_clear(wf_aligner,wf_aligner->distance_metric);
  }
  // Resize CIGAR
  cigar_resize(&wf_aligner->cigar,2*(pattern_length+text_length));
  // Clear BT-Buffer and Slab
  if (wf_aligner->bt_buffer) wf_backtrace_buffer_clear(wf_aligner->bt_buffer);
  wavefront_slab_clear(wf_aligner->wavefront_slab);
  // Resize pattern & text
  wf_aligner->pattern_length = pattern_length;
  wf_aligner->text_length = text_length;
}
void wavefront_aligner_delete(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  mm_allocator_t* const mm_allocator = wf_aligner->mm_allocator;
  // Null wavefront
  wavefront_free(wf_aligner->wavefront_null,mm_allocator);
  mm_allocator_free(mm_allocator,wf_aligner->wavefront_null);
  // Victim wavefront
  wavefront_free(wf_aligner->wavefront_victim,mm_allocator);
  mm_allocator_free(mm_allocator,wf_aligner->wavefront_victim);
  // Wavefront components
  wavefront_components_free(wf_aligner,wf_aligner->distance_metric);
  // CIGAR
  cigar_free(&wf_aligner->cigar);
  // BT-Buffer and Slab
  if (wf_aligner->bt_buffer) wf_backtrace_buffer_delete(wf_aligner->bt_buffer);
  wavefront_slab_delete(wf_aligner->wavefront_slab);
  // MM
  const bool mm_allocator_own = wf_aligner->mm_allocator_own;
  mm_allocator_free(mm_allocator,wf_aligner); // Handler
  if (mm_allocator_own) {
    mm_allocator_delete(mm_allocator);
  }
}
/*
 * Configuration
 */
void wavefront_aligner_set_alignment_end_to_end(
    wavefront_aligner_t* const wf_aligner) {
  wf_aligner->alignment_form.span = alignment_end2end;
}
void wavefront_aligner_set_alignment_free_ends(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free) {
  wf_aligner->alignment_form.span = alignment_endsfree;
  wf_aligner->alignment_form.pattern_begin_free = pattern_begin_free;
  wf_aligner->alignment_form.pattern_end_free = pattern_end_free;
  wf_aligner->alignment_form.text_begin_free = text_begin_free;
  wf_aligner->alignment_form.text_end_free = text_end_free;
}
void wavefront_aligner_set_reduction_none(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_reduction_set_none(&wf_aligner->reduction);
}
void wavefront_aligner_set_reduction_adaptive(
    wavefront_aligner_t* const wf_aligner,
    const int min_wavefront_length,
    const int max_distance_threshold) {
  wavefront_reduction_set_adaptive(&wf_aligner->reduction,
      min_wavefront_length,max_distance_threshold);
}
void wavefront_aligner_set_max_alignment_score(
    wavefront_aligner_t* const wf_aligner,
    const int max_alignment_score) {
  wf_aligner->alignment_form.max_alignment_score = max_alignment_score;
}
void wavefront_aligner_set_max_memory_used(
    wavefront_aligner_t* const wf_aligner,
    const uint64_t max_memory_used) {
  wf_aligner->max_memory_used = max_memory_used;
}
/*
 * Utils
 */
uint64_t wavefront_aligner_get_size(
    wavefront_aligner_t* const wf_aligner) {
  const uint64_t bt_buffer_size = (wf_aligner->bt_buffer) ?
      wf_backtrace_buffer_get_size(wf_aligner->bt_buffer) : 0;
  const uint64_t slab_size = wavefront_slab_get_size(wf_aligner->wavefront_slab);
  return bt_buffer_size + slab_size;
}


