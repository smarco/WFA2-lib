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

#include "wavefront_aligner.h"
#include "wavefront_components.h"
#include "wavefront_heuristic.h"
#include "wavefront_plot.h"

/*
 * Configuration
 */
#define PATTERN_LENGTH_INIT 1000
#define TEXT_LENGTH_INIT    1000
#define SEQUENCES_PADDING     10

/*
 * Error messages
 */
char* wf_error_msg[] =
{
  /* WF_STATUS_OOM                  == -3 */ "[WFA] Alignment failed. Maximum memory threshold reached",
  /* WF_STATUS_MAX_SCORE_REACHED    == -2 */ "[WFA] Alignment failed. Maximum score reached",
  /* WF_STATUS_HEURISTICALY_DROPPED == -1 */ "[WFA] Alignment dropped heuristically",
  /* WF_STATUS_SUCCESSFUL           ==  0 */ "[WFA] Alignment successful",
  /* WF_STATUS_IN_PROGRESS          ==  1 */ "[WFA] Alignment in progress",
};
char* wavefront_align_strerror(const int wf_error_code) {
  return wf_error_msg[wf_error_code+3];
}
/*
 * Alignment status
 */
void wavefront_align_status_clear(
    wavefront_align_status_t* const wf_align_status) {
  wf_align_status->status = WF_STATUS_IN_PROGRESS;
  wf_align_status->score = 0;
}
/*
 * Setup
 */
void wavefront_aligner_set_penalties(
    wavefront_aligner_t* const wf_aligner,
    wavefront_aligner_attr_t* const attributes) {
  switch (attributes->distance_metric) {
    case indel:
      wavefronts_penalties_set_indel(&wf_aligner->penalties);
      break;
    case edit:
      wavefronts_penalties_set_edit(&wf_aligner->penalties);
      break;
    case gap_linear:
      wavefronts_penalties_set_linear(
          &wf_aligner->penalties,
          &attributes->linear_penalties,
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
void wavefront_aligner_set_heuristic(
    wavefront_aligner_t* const wf_aligner,
    wavefront_aligner_attr_t* const attributes) {
  // Parameters
  wavefront_heuristic_t* const wf_heuristic = &attributes->heuristic;
  // Select and configure heuristics
  if (wf_heuristic->strategy == wf_heuristic_none) {
    wavefront_heuristic_set_none(&wf_aligner->heuristic);
  } else {
    if ((wf_heuristic->strategy & wf_heuristic_banded_static) != 0) {
      wavefront_heuristic_set_banded_static(&wf_aligner->heuristic,
          wf_heuristic->min_k,wf_heuristic->max_k);
    }
    if ((wf_heuristic->strategy & wf_heuristic_banded_adaptive) != 0) {
      wavefront_heuristic_set_banded_adaptive(&wf_aligner->heuristic,
          wf_heuristic->min_k,wf_heuristic->max_k,wf_heuristic->steps_between_cutoffs);
    }
    if ((wf_heuristic->strategy & wf_heuristic_wfadaptive) != 0) {
      wavefront_heuristic_set_wfadaptive(
          &wf_aligner->heuristic,wf_heuristic->min_wavefront_length,
          wf_heuristic->max_distance_threshold,wf_heuristic->steps_between_cutoffs);
    }
    if ((wf_heuristic->strategy & wf_heuristic_xdrop) != 0) {
      wavefront_heuristic_set_xdrop(&wf_aligner->heuristic,
          wf_heuristic->xdrop,wf_heuristic->steps_between_cutoffs);
    }
    if ((wf_heuristic->strategy & wf_heuristic_zdrop) != 0) {
      wavefront_heuristic_set_zdrop(&wf_aligner->heuristic,
          wf_heuristic->zdrop,wf_heuristic->steps_between_cutoffs);
    }
  }
}
void wavefront_aligner_set_system(
    wavefront_aligner_t* const wf_aligner,
    alignment_system_t* const system) {
  // Copy all parameters
  wf_aligner->system = *system;
  // Reset effective limits
  wf_aligner->system.max_memory_compact = BUFFER_SIZE_256M;
  wf_aligner->system.max_memory_resident = BUFFER_SIZE_256M + BUFFER_SIZE_256M;
  switch (wf_aligner->memory_mode) {
    case wavefront_memory_med:
      wf_aligner->system.max_partial_compacts = 4;
      break;
    case wavefront_memory_low:
      wf_aligner->system.max_partial_compacts = 1;
      break;
    default:
      break;
  }
  // Profile
  timer_reset(&wf_aligner->system.timer);
}
wavefront_aligner_t* wavefront_aligner_new(
    wavefront_aligner_attr_t* attributes) {
  // Attributes
  const int pattern_length = PATTERN_LENGTH_INIT;
  const int text_length = TEXT_LENGTH_INIT;
  if (attributes == NULL) attributes = &wavefront_aligner_attr_default;
  const bool score_only = (attributes->alignment_scope == compute_score);
  const bool memory_modular = (attributes->memory_mode > 0) || score_only;
  const bool bt_piggyback = (attributes->memory_mode > 0) && !score_only;
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
  // Configuration
  wf_aligner->pattern_length = pattern_length;
  wf_aligner->text_length = text_length;
  wf_aligner->sequences = NULL;
  wf_aligner->alignment_scope = attributes->alignment_scope;
  wf_aligner->alignment_form = attributes->alignment_form;
  wavefront_aligner_set_penalties(wf_aligner,attributes);
  // Memory mode
  wf_aligner->memory_mode = attributes->memory_mode;
  // Heuristic strategy
  wavefront_aligner_set_heuristic(wf_aligner,attributes);
  // Custom matching functions
  wf_aligner->match_funct = attributes->match_funct;
  wf_aligner->match_funct_arguments = attributes->match_funct_arguments;
  // Wavefront components
  wavefront_components_allocate(
      &wf_aligner->wf_components,pattern_length,text_length,
      &wf_aligner->penalties,memory_modular,bt_piggyback,mm_allocator);
  // CIGAR
  cigar_allocate(&wf_aligner->cigar,2*(pattern_length+text_length),mm_allocator);
  // Display
  wf_aligner->plot_params = attributes->plot_params;
  if (attributes->plot_params.plot_enabled) {
    wavefront_plot_allocate(&wf_aligner->wf_plot,
        wf_aligner->penalties.distance_metric,
        pattern_length,text_length,
        &wf_aligner->plot_params);
  }
  // System
  wavefront_aligner_set_system(wf_aligner,&attributes->system);
  // Return
  return wf_aligner;
}
void wavefront_aligner_resize(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  // Configure sequences and status
  wf_aligner->pattern_length = pattern_length;
  wf_aligner->text_length = text_length;
  if (wf_aligner->match_funct == NULL) {
    if (wf_aligner->sequences != NULL) strings_padded_delete(wf_aligner->sequences);
    wf_aligner->sequences = strings_padded_new_rhomb(
            pattern,pattern_length,text,text_length,
            SEQUENCES_PADDING,wf_aligner->mm_allocator);
    wf_aligner->pattern = wf_aligner->sequences->pattern_padded;
    wf_aligner->text = wf_aligner->sequences->text_padded;
  } else {
    wf_aligner->sequences = NULL;
    wf_aligner->pattern = NULL;
    wf_aligner->text = NULL;
  }
  wavefront_align_status_clear(&wf_aligner->align_status);
  // Heuristics clear
  wavefront_heuristic_clear(&wf_aligner->heuristic);
  // Wavefront components
  wavefront_components_resize(&wf_aligner->wf_components,
      pattern_length,text_length,&wf_aligner->penalties);
  // CIGAR
  cigar_resize(&wf_aligner->cigar,2*(pattern_length+text_length));
  // Slab
  wavefront_slab_clear(wf_aligner->wavefront_slab);
  // Display
  if (wf_aligner->plot_params.plot_enabled) {
    wavefront_plot_free(&wf_aligner->wf_plot);
    wavefront_plot_allocate(&wf_aligner->wf_plot,
        wf_aligner->penalties.distance_metric,
        pattern_length,text_length,
        &wf_aligner->plot_params);
  }
  // System
  wavefront_aligner_set_system(wf_aligner,&wf_aligner->system);
}
void wavefront_aligner_reap(
    wavefront_aligner_t* const wf_aligner) {
  // Padded sequences
  if (wf_aligner->sequences != NULL) {
    strings_padded_delete(wf_aligner->sequences);
  }
  // Wavefront components
  wavefront_components_reap(&wf_aligner->wf_components);
  // Slab
  wavefront_slab_reap(wf_aligner->wavefront_slab,wf_slab_reap_all);
}
void wavefront_aligner_delete(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  mm_allocator_t* const mm_allocator = wf_aligner->mm_allocator;
  // Padded sequences
  if (wf_aligner->sequences != NULL) {
    strings_padded_delete(wf_aligner->sequences);
  }
  // Wavefront components
  wavefront_components_free(&wf_aligner->wf_components);
  // CIGAR
  cigar_free(&wf_aligner->cigar);
  // Slab
  wavefront_slab_delete(wf_aligner->wavefront_slab);
  // Display
  if (wf_aligner->plot_params.plot_enabled) {
    wavefront_plot_free(&wf_aligner->wf_plot);
  }
  // MM
  const bool mm_allocator_own = wf_aligner->mm_allocator_own;
  mm_allocator_free(mm_allocator,wf_aligner); // Handler
  if (mm_allocator_own) {
    mm_allocator_delete(mm_allocator);
  }
}
/*
 * Span configuration
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
/*
 * Heuristic configuration
 */
void wavefront_aligner_set_heuristic_none(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_heuristic_set_none(&wf_aligner->heuristic);
}
void wavefront_aligner_set_heuristic_banded_static(
    wavefront_aligner_t* const wf_aligner,
    const int band_min_k,
    const int band_max_k) {
  wavefront_heuristic_set_banded_static(&wf_aligner->heuristic,band_min_k,band_max_k);
}
void wavefront_aligner_set_heuristic_banded_adaptive(
    wavefront_aligner_t* const wf_aligner,
    const int band_min_k,
    const int band_max_k,
    const int score_steps) {
  wavefront_heuristic_set_banded_adaptive(&wf_aligner->heuristic,band_min_k,band_max_k,score_steps);
}
void wavefront_aligner_set_heuristic_wfadaptive(
    wavefront_aligner_t* const wf_aligner,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int score_steps) {
  wavefront_heuristic_set_wfadaptive(&wf_aligner->heuristic,min_wavefront_length,max_distance_threshold,score_steps);
}
void wavefront_aligner_set_heuristic_xdrop(
    wavefront_aligner_t* const wf_aligner,
    const int xdrop,
    const int score_steps) {
  wavefront_heuristic_set_xdrop(&wf_aligner->heuristic,xdrop,score_steps);
}
void wavefront_aligner_set_heuristic_zdrop(
    wavefront_aligner_t* const wf_aligner,
    const int ydrop,
    const int score_steps) {
  wavefront_heuristic_set_zdrop(&wf_aligner->heuristic,ydrop,score_steps);
}
/*
 * Match-funct configuration
 */
void wavefront_aligner_set_match_funct(
    wavefront_aligner_t* const wf_aligner,
    int (*match_funct)(int,int,void*),
    void* const match_funct_arguments) {
  wf_aligner->match_funct = match_funct;
  wf_aligner->match_funct_arguments = match_funct_arguments;
}
/*
 * System configuration
 */
void wavefront_aligner_set_max_alignment_score(
    wavefront_aligner_t* const wf_aligner,
    const int max_alignment_score) {
  wf_aligner->alignment_form.max_alignment_score = max_alignment_score;
}
void wavefront_aligner_set_max_memory(
    wavefront_aligner_t* const wf_aligner,
    const uint64_t max_memory_compact,
    const uint64_t max_memory_resident,
    const uint64_t max_memory_abort) {
  wf_aligner->system.max_memory_compact = max_memory_compact;
  wf_aligner->system.max_memory_resident = max_memory_resident;
  wf_aligner->system.max_memory_abort = max_memory_abort;
}
/*
 * Utils
 */
uint64_t wavefront_aligner_get_size(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  // Compute size
  const uint64_t bt_buffer_size = (wf_components->bt_buffer) ?
      wf_backtrace_buffer_get_size_allocated(wf_components->bt_buffer) : 0;
  const uint64_t slab_size = wavefront_slab_get_size(wf_aligner->wavefront_slab);
  return bt_buffer_size + slab_size;
}
/*
 * Display
 */
void wavefront_aligner_print_status(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  // Approximate progress
  const int dist_total = MAX(wf_aligner->text_length,wf_aligner->pattern_length);
  int s = (wf_components->memory_modular) ? score%wf_components->max_score_scope : score;
  wavefront_t* wavefront = wf_components->mwavefronts[s];
  if (wavefront==NULL && s>0) {
    s = (wf_components->memory_modular) ? (score-1)%wf_components->max_score_scope : (score-1);
    wavefront = wf_components->mwavefronts[s];
  }
  int dist_max = -1, wf_len = -1, k;
  if (wavefront!=NULL) {
    wf_offset_t* const offsets = wavefront->offsets;
    for (k=wavefront->lo;k<=wavefront->hi;++k) {
      const int dist = MAX(WAVEFRONT_V(k,offsets[k]),WAVEFRONT_H(k,offsets[k]));
      dist_max = MAX(dist_max,dist);
    }
    wf_len = wavefront->hi-wavefront->lo+1;
  }
  // Memory used
  const uint64_t slab_size = wavefront_slab_get_size(wf_aligner->wavefront_slab);
  const uint64_t bt_buffer_used = (wf_components->bt_buffer) ?
      wf_backtrace_buffer_get_size_used(wf_components->bt_buffer) : 0;
  // Print one-line status
  fprintf(stream,
      "[WFA] SequenceLength=(%d,%d) Score %d (~ %2.3f%% aligned). "
      "MemoryUsed(WF-Slab,BT-buffer)=(%lu MB,%lu MB). "
      "Wavefronts ~ %2.3f Moffsets\n",
      wf_aligner->pattern_length,wf_aligner->text_length,
      score,(dist_max>=0) ? (100.0f*(float)dist_max/(float)dist_total) : -1.0f,
      CONVERT_B_TO_MB(slab_size),CONVERT_B_TO_MB(bt_buffer_used),
      (wf_len>=0) ? (float)wf_len/1000000.0f : -1.0f);
}


