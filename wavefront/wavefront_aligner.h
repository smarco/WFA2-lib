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
 * DESCRIPTION: WaveFront aligner data structure
 */

#ifndef WAVEFRONT_ALIGNER_H_
#define WAVEFRONT_ALIGNER_H_

#include "utils/commons.h"
#include "utils/heatmap.h"
#include "utils/string_padded.h"
#include "system/profiler_counter.h"
#include "system/profiler_timer.h"
#include "system/mm_allocator.h"
#include "system/mm_stack.h"
#include "alignment/cigar.h"
#include "wavefront_slab.h"
#include "wavefront_penalties.h"
#include "wavefront_attributes.h"
#include "wavefront_components.h"
#include "wavefront_align.h"
#include "wavefront_bialign.h"

/*
 * Error codes & messages
 */
#define WF_STATUS_SUCCESSFUL            0
#define WF_STATUS_IN_PROGRESS           1
#define WF_STATUS_HEURISTICALY_DROPPED -1
#define WF_STATUS_MAX_SCORE_REACHED    -2
#define WF_STATUS_OOM                  -3
extern char* wf_error_msg[5];
char* wavefront_align_strerror(const int wf_error_code);

/*
 * Alignment status
 */
typedef struct _wavefront_aligner_t wavefront_aligner_t;
typedef struct {
  // Status
  int status;                                                     // Status code
  int score;                                                      // Current WF-alignment score
  uint64_t memory_used;                                           // Total memory used
  // Wavefront alignment functions
  void (*wf_align_compute)(wavefront_aligner_t* const,const int); // WF Compute function
  int (*wf_align_extend)(wavefront_aligner_t* const,const int);   // WF Extend function
} wavefront_align_status_t;

/*
 * Wavefront Aligner
 */
typedef struct _wavefront_aligner_t {
  // Status
  wavefront_align_status_t align_status;      // Current alignment status
  // Sequences
  strings_padded_t* sequences;                // Padded sequences
  char* pattern;                              // Pattern sequence (padded)
  int pattern_length;                         // Pattern length
  char* text;                                 // Text sequence (padded)
  int text_length;                            // Text length
  // Sequences (reversed)
  strings_padded_t* sequences_rev;            // Padded sequences reversed
  char* pattern_rev;                          // Pattern sequence (padded & reversed)
  int pattern_length_rev;                     // Pattern length reversed
  char* text_rev;                             // Text sequence (padded & reversed)
  int text_length_rev;                        // Text length reversed
  // Alignment Attributes
  alignment_scope_t alignment_scope;          // Alignment scope (score only or full-CIGAR)
  alignment_form_t alignment_form;            // Alignment form (end-to-end/ends-free)
  wavefronts_penalties_t penalties;           // Alignment penalties
  wavefront_heuristic_t heuristic;            // Heuristic's parameters
  wavefront_memory_t memory_mode;             // Wavefront memory strategy (modular wavefronts and piggyback)
  // Custom function to compare sequences
  alignment_match_funct_t match_funct;        // Custom matching function (match(v,h,args))
  void* match_funct_arguments;                // Generic arguments passed to matching function (args)
  // Wavefront components
  wavefront_components_t wf_components;       // Wavefront components
  affine2p_matrix_type component_begin;       // Alignment begin component
  affine2p_matrix_type component_end;         // Alignment end component
  wavefront_pos_t alignment_end_pos;          // Alignment end position
  // Bidirectional Alignment
  bool bidirectional_alignment;               // Enable bidirectional WFA alignment
  wavefront_aligner_t* aligner_forward;       // Forward aligner
  wavefront_aligner_t* aligner_reverse;       // Reverse aligner
  wf_bialign_breakpoint_t bialign_breakpoint; // Breakpoint of two wavefronts (bialigner)
  // CIGAR
  cigar_t cigar;                              // Alignment CIGAR
  // MM
  bool mm_allocator_own;                      // Ownership of MM-Allocator
  mm_allocator_t* mm_allocator;               // MM-Allocator
  wavefront_slab_t* wavefront_slab;           // MM-Wavefront-Slab (Allocates/Reuses the individual wavefronts)
  // Display
  wavefront_plot_params_t plot_params;        // Wavefront plot parameters
  wavefront_plot_t wf_plot;                   // Wavefront plot
  // System
  alignment_system_t system;                  // System related parameters
} wavefront_aligner_t;

/*
 * Setup
 */
wavefront_aligner_t* wavefront_aligner_new(
    wavefront_aligner_attr_t* attributes);
void wavefront_aligner_resize(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const bool reverse_sequences);
void wavefront_aligner_reap(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_delete(
    wavefront_aligner_t* const wf_aligner);

/*
 * Span configuration
 */
void wavefront_aligner_set_alignment_end_to_end(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_set_alignment_free_ends(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free);

/*
 * Heuristic configuration
 */
void wavefront_aligner_set_heuristic_none(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_set_heuristic_banded_static(
    wavefront_aligner_t* const wf_aligner,
    const int band_min_k,
    const int band_max_k);
void wavefront_aligner_set_heuristic_banded_adaptive(
    wavefront_aligner_t* const wf_aligner,
    const int band_min_k,
    const int band_max_k,
    const int score_steps);
void wavefront_aligner_set_heuristic_wfadaptive(
    wavefront_aligner_t* const wf_aligner,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int score_steps);
void wavefront_aligner_set_heuristic_xdrop(
    wavefront_aligner_t* const wf_aligner,
    const int xdrop,
    const int score_steps);
void wavefront_aligner_set_heuristic_zdrop(
    wavefront_aligner_t* const wf_aligner,
    const int ydrop,
    const int score_steps);

/*
 * Match-funct configuration
 */
void wavefront_aligner_set_match_funct(
    wavefront_aligner_t* const wf_aligner,
    int (*match_funct)(int,int,void*),
    void* const match_funct_arguments);

/*
 * System configuration
 */
void wavefront_aligner_set_max_alignment_score(
    wavefront_aligner_t* const wf_aligner,
    const int max_alignment_score);
void wavefront_aligner_set_max_memory(
    wavefront_aligner_t* const wf_aligner,
    const uint64_t max_memory_resident,
    const uint64_t max_memory_abort);

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
