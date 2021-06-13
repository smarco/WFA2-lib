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

#include "gap_affine2p/affine2p_penalties.h"
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
 * Default parameters
 */
wavefront_aligner_attr_t wavefront_aligner_attr_default = {
    // Distance model & Penalties
    .distance_metric = gap_affine, // TODO: DU wants lineal
    .alignment_scope = alignment_scope_alignment, // TODO: I don't know what DU wants, honestly
    .lineal_penalties = {
        .match = 0,
        .mismatch = 4,
        .insertion = 2,
        .deletion  = 2,
    },
    .affine_penalties = {
        .match = 0,
        .mismatch = 4,
        .gap_opening = 6,
        .gap_extension = 2, // 1
    },
    .affine2p_penalties = {
        .match = 0,
        .mismatch = 4,
        .gap_opening1 = 6,
        .gap_extension1 = 2,
        .gap_opening2 = 24,
        .gap_extension2 = 1, // TODO: Should be 2 (???)
    },
    // Reduction
    .reduction = {
        .reduction_strategy = wavefront_reduction_none, // TODO: DU wants dynamic
        .min_wavefront_length = -1, // 10,
        .max_distance_threshold = -1, // 50,
    },
    // Memory model
    .low_memory = false,  // TODO: DU wants (low_memory if length >= 10000) // FIXME
    // MM
    .mm_allocator = NULL, // Use private MM
};

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
  const bool score_only = (attributes->alignment_scope == alignment_scope_score);
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
  wf_aligner->memory_modular = memory_modular;
  wf_aligner->bt_piggyback = bt_piggyback;
  wavefront_set_penalties(wf_aligner,attributes); // Set penalties
  // Reduction strategy
  if (attributes->reduction.reduction_strategy == wavefront_reduction_dynamic) {
    wavefront_reduction_set_dynamic(
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
  cigar_allocate(&wf_aligner->cigar,pattern_length+text_length,mm_allocator);
  // Return
  return wf_aligner;
}
void wavefront_aligner_clear(
    wavefront_aligner_t* const wf_aligner) {
  // Clear wavefront components
  if (wf_aligner->memory_modular) wavefront_components_clear(wf_aligner,wf_aligner->distance_metric);
  // Clear CIGAR
  cigar_clear(&wf_aligner->cigar);
  // BT-Buffer
  if (wf_aligner->bt_buffer) wf_backtrace_buffer_clear(wf_aligner->bt_buffer);
  // Clear Slab
  wavefront_slab_clear(wf_aligner->wavefront_slab,false);
}
void wavefront_aligner_clear__resize(
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
  if (pattern_length+text_length > wf_aligner->pattern_length+wf_aligner->text_length) {
    cigar_free(&wf_aligner->cigar); // FIXME
    cigar_allocate(&wf_aligner->cigar,pattern_length+text_length,wf_aligner->mm_allocator); // FIXME
  } else {
    cigar_clear(&wf_aligner->cigar);
  }
  // BT-Buffer
  if (wf_aligner->bt_buffer) wf_backtrace_buffer_clear(wf_aligner->bt_buffer);
  // Clear Slab
  wavefront_slab_clear(wf_aligner->wavefront_slab,false);
  // Resize
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
  // MM
  if (wf_aligner->bt_buffer) wf_backtrace_buffer_delete(wf_aligner->bt_buffer);
  wavefront_slab_delete(wf_aligner->wavefront_slab); // Slab
  if (wf_aligner->mm_allocator_own) {
    mm_allocator_delete(wf_aligner->mm_allocator);
  }
  mm_allocator_free(mm_allocator,wf_aligner); // Handler
}


