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

#ifndef WAVEFRONT_ALIGNER_H_
#define WAVEFRONT_ALIGNER_H_

#include "utils/commons.h"
#include "utils/heatmap.h"
#include "system/profiler_counter.h"
#include "system/profiler_timer.h"
#include "system/mm_allocator.h"
#include "system/mm_stack.h"
#include "alignment/cigar.h"
#include "wavefront_slab.h"
#include "wavefront_penalties.h"
#include "wavefront_attributes.h"
#include "wavefront_components.h"

/*
 * Wavefront Aligner
 */
typedef struct _wavefront_aligner_t {
  // Sequences
  char* pattern;                               // Pattern sequence (padded)
  int pattern_length;                          // Pattern length
  char* text;                                  // Text sequence (padded)
  int text_length;                             // Text length
  // Alignment Attributes
  alignment_scope_t alignment_scope;           // Alignment scope (score only or full-CIGAR)
  alignment_form_t alignment_form;             // Alignment form (end-to-end/ends-free)
  wavefronts_penalties_t penalties;            // Alignment penalties
  wavefront_reduction_t reduction;             // Reduction parameters
  // Wavefront components
  wavefront_components_t wf_components;        // Wavefront components
  // CIGAR
  cigar_t cigar;                               // Alignment CIGAR
  // MM
  bool mm_allocator_own;                       // Ownership of MM-Allocator
  mm_allocator_t* mm_allocator;                // MM-Allocator
  wavefront_slab_t* wavefront_slab;            // MM-Wavefront-Slab (Allocates/Reuses the individual wavefronts)
  // Display
  wavefront_plot_params_t plot_params;         // Wavefront plot parameters
  wavefront_plot_t wf_plot;                    // Wavefront plot
  // System
  alignment_system_t system;                   // System related parameters
} wavefront_aligner_t;

/*
 * Setup
 */
wavefront_aligner_t* wavefront_aligner_new(
    wavefront_aligner_attr_t* attributes);
void wavefront_aligner_resize(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length);
void wavefront_aligner_reap(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_delete(
    wavefront_aligner_t* const wf_aligner);

/*
 * Configuration
 */
void wavefront_aligner_set_alignment_end_to_end(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_set_alignment_free_ends(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free);

void wavefront_aligner_set_reduction_none(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_set_reduction_adaptive(
    wavefront_aligner_t* const wf_aligner,
    const int min_wavefront_length,
    const int max_distance_threshold);

void wavefront_aligner_set_max_alignment_score(
    wavefront_aligner_t* const wf_aligner,
    const int max_alignment_score);
void wavefront_aligner_set_max_memory_used(
    wavefront_aligner_t* const wf_aligner,
    const uint64_t max_memory_used);

/*
 * Utils
 */
uint64_t wavefront_aligner_get_size(
    wavefront_aligner_t* const wf_aligner);

/*
 * Display
 */
void wavefront_aligner_print_status(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner,
    const int current_score);

#endif /* WAVEFRONT_ALIGNER_H_ */
