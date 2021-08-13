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
 * DESCRIPTION: WaveFront aligner data structure attributes
 */

#ifndef WAVEFRONT_ATTRIBUTES_H_
#define WAVEFRONT_ATTRIBUTES_H_

#include "utils/commons.h"
#include "alignment/cigar.h"
#include "gap_affine/affine_penalties.h"
#include "gap_affine2p/affine2p_penalties.h"
#include "gap_lineal/lineal_penalties.h"
#include "system/mm_allocator.h"

#include "wavefront_penalties.h"
#include "wavefront_reduction.h"
#include "wavefront_plot.h"
#include "wavefront_display.h"

/*
 * Alignment scope
 */
typedef enum {
  compute_score,          // Only distance/score
  compute_alignment,      // Full alignment CIGAR
} alignment_scope_t;
typedef enum {
  alignment_end2end,       // End-to-end alignment (aka global)
  alignment_endsfree,      // Ends-free alignment  (semiglobal, glocal, etc)
} alignment_span_t;
typedef struct {
  // Mode
  alignment_span_t span;   // Alignment form (End-to-end/Ends-free)
  // Ends-free
  int pattern_begin_free;  // Allow free-gap at the beginning of the pattern
  int pattern_end_free;    // Allow free-gap at the end of the pattern
  int text_begin_free;     // Allow free-gap at the beginning of the text
  int text_end_free;       // Allow free-gap at the end of the text
  // Limits
  int max_alignment_score; // Maximum score allowed before quit
} alignment_form_t;

/*
 * Alignment system configuration
 */
typedef struct {
  // Global
  int global_probe_interval;          // Score-ticks interval to check any limits
  // BT-Buffer compacting
  int bt_compact_probe_interval;      // Score-ticks interval to check BT-buffer compacting
  uint64_t bt_compact_max_memory;     // Maximum BT-buffer memory (allowed before trying compacting)
  uint64_t bt_compact_max_memory_eff; // Effective maximum BT-buffer memory
  // Memory
  uint64_t max_memory_used;           // Maximum memory allowed to used before quit
  uint64_t max_memory_resident;       // Maximum memory allowed to be buffered before reap
  // Misc
  bool verbose;                       // Verbose (regulates messages during alignment)
} alignment_system_t;

/*
 * Wavefront Aligner Attributes
 */
typedef struct {
  // Distance model
  distance_metric_t distance_metric;         // Alignment metric/distance used
  alignment_scope_t alignment_scope;         // Alignment scope (score only or full-CIGAR)
  alignment_form_t alignment_form;           // Alignment mode (end-to-end/ends-free)
  // Penalties
  lineal_penalties_t lineal_penalties;       // Gap-lineal penalties (placeholder)
  affine_penalties_t affine_penalties;       // Gap-affine penalties (placeholder)
  affine2p_penalties_t affine2p_penalties;   // Gap-affine-2p penalties (placeholder)
  // Reduction strategy
  wavefront_reduction_t reduction;           // Wavefront reduction
  // Memory model
  bool low_memory;                           // Use low-memory strategy (modular wavefronts and piggyback)
  // External MM (instead of allocating one inside)
  mm_allocator_t* mm_allocator;              // MM-Allocator
  // Display
  wavefront_plot_params_t plot_params;       // Wavefront plot
  // System
  alignment_system_t system;                 // System related parameters
} wavefront_aligner_attr_t;

/*
 * Default parameters
 */
extern wavefront_aligner_attr_t wavefront_aligner_attr_default;

#endif /* WAVEFRONT_ATTRIBUTES_H_ */
