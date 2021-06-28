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
#include "system/profiler_counter.h"
#include "system/profiler_timer.h"
#include "system/mm_allocator.h"
#include "system/mm_stack.h"
#include "alignment/cigar.h"
#include "wavefront_slab.h"
#include "wavefront_penalties.h"
#include "wavefront_attributes.h"

/*
 * Wavefront Aligner
 */
typedef struct {
  // Dimensions
  int pattern_length;                          // Pattern length
  int text_length;                             // Text length
  int max_score_scope;                         // Maximum score-difference between dependent wavefronts
  // Alignment Attributes
  distance_metric_t distance_metric;           // Alignment metric/distance used
  alignment_scope_t alignment_scope;           // Alignment scope (score only or full-CIGAR)
  alignment_form_t alignment_form;             // Alignment form (end-to-end/ends-free)
  wavefronts_penalties_t penalties;            // Alignment penalties
  // Wavefront Attributes
  wavefront_reduction_t reduction;             // Reduction parameters
  bool memory_modular;                         // Memory strategy (modular wavefronts)
  bool bt_piggyback;                           // Backtrace Piggyback
  // Wavefront components
  int num_wavefronts;                          // Total number of allocated wavefronts
  wavefront_t** mwavefronts;                   // M-wavefronts
  wavefront_t** i1wavefronts;                  // I1-wavefronts
  wavefront_t** i2wavefronts;                  // I2-wavefronts
  wavefront_t** d1wavefronts;                  // D1-wavefronts
  wavefront_t** d2wavefronts;                  // D2-wavefronts
  wavefront_t* wavefront_null;                 // Null wavefront (orthogonal reading)
  wavefront_t* wavefront_victim;               // Dummy wavefront (orthogonal writing)
  // CIGAR
  cigar_t cigar;                               // Alignment CIGAR
  wf_backtrace_buffer_t* bt_buffer;            // Backtrace Buffer
  // MM
  bool mm_allocator_own;                       // Ownership of MM-Allocator
  mm_allocator_t* mm_allocator;                // MM-Allocator
  wavefront_slab_t* wavefront_slab;            // MM-Wavefront-Slab (Allocates/Reuses the individual wavefronts)
  // Limits
  int limit_probe_interval;                    // Score-ticks to check limits
  uint64_t max_memory_used;                    // Maximum memory allowed to used before quit
  uint64_t max_resident_memory;                // Maximum memory allowed to be buffered before reap
} wavefront_aligner_t;

/*
 * Setup
 */
wavefront_aligner_t* wavefront_aligner_new(
    const int pattern_length,
    const int text_length,
    wavefront_aligner_attr_t* attributes);
void wavefront_aligner_reap(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_clear(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_resize(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length);
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

#endif /* WAVEFRONT_ALIGNER_H_ */
