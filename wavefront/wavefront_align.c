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
 * DESCRIPTION: WaveFront alignment module for sequence pairwise alignment
 */

#include "utils/string_padded.h"
#include "wavefront_align.h"
#include "wavefront_extend.h"
#include "wavefront_compute.h"
#include "wavefront_compute_affine.h"
#include "wavefront_compute_affine2p.h"
#include "wavefront_backtrace.h"

/*
 * Error messages
 */
char* wf_error_msg[] =
{
  /* WF_ALIGN_SUCCESSFUL */ "[WFA] Alignment successful",
  /* WF_ALIGN_MAX_SCORE */  "[WFA] Alignment failed. Maximum score reached",
  /* WF_ALIGN_OOM */        "[WFA] Alignment failed. Maximum memory threshold reached",
};
char* wavefront_align_strerror(
    const int wf_error_code) {
  return wf_error_msg[-wf_error_code];
}

/*
 * Checks
 */
void wavefront_form_endsfree_check(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length) {
  alignment_form_t* const form = &wf_aligner->alignment_form;
  if (form->pattern_begin_free > pattern_length ||
      form->pattern_end_free > pattern_length ||
      form->text_begin_free > text_length ||
      form->text_end_free > text_length) {
    fprintf(stderr,"[WFA] Ends-free parameters must be not larger than the sequences "
        "(P0=%d,Pf=%d,T0=%d,Tf=%d). Must be (P0<=|P|,Pf<=|P|,T0<=|T|,Tf<=|T|) where (|P|,|T|)=(%d,%d)\n",
        form->pattern_begin_free,form->pattern_end_free,
        form->text_begin_free,form->text_end_free,
        pattern_length,text_length);
    exit(1);
  }
}
/*
 * Limits
 */
int wavefront_align_reached_limits(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Check alignment-score limit
  if (score >= wf_aligner->alignment_form.max_alignment_score) {
    wf_aligner->cigar.score = wf_aligner->alignment_form.max_alignment_score;
    return WF_ALIGN_MAX_SCORE;
  }
  // Global probing interval
  alignment_system_t* const system = &wf_aligner->system;
  if ((score%system->probe_interval_global) != 0) return 0;
  if (system->verbose) {
    wavefront_aligner_print_status(stderr,wf_aligner,score); // DEBUG
  }
  // BT-Buffer
  wavefront_components_t*const wf_components = &wf_aligner->wf_components;
  if (wf_components->bt_buffer!=NULL && (score%system->probe_interval_compact)==0) {
    uint64_t bt_memory = wf_backtrace_buffer_get_size_used(wf_components->bt_buffer);
    // Check BT-buffer memory
    if (bt_memory > system->max_memory_compact) {
      // Compact BT-buffer
      wavefront_components_compact_bt_buffer(wf_components,score,wf_aligner->system.verbose);
      // Set new buffer limit
      bt_memory = wf_backtrace_buffer_get_size_used(wf_components->bt_buffer);
      uint64_t proposed_mem = (double)bt_memory * TELESCOPIC_FACTOR;
      if (system->max_memory_compact < proposed_mem && proposed_mem < system->max_memory_abort) {
        proposed_mem = system->max_memory_compact;
      }
    }
  }
  // Check overall memory used
  const uint64_t wf_memory_used = wavefront_aligner_get_size(wf_aligner);
  if (wf_memory_used > system->max_memory_abort) {
    return WF_ALIGN_OOM;
  }
  // Otherwise OK
  return 0;
}
/*
 * End-to-end alignment (Global)
 */
void wavefront_align_end2end_initialize(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  const int max_score_scope = wf_components->max_score_scope;
  // Init wavefronts
  const int effective_lo = -(max_score_scope+1);
  const int effective_hi = (max_score_scope+1);
  wf_components->mwavefronts[0] = wavefront_slab_allocate(
      wf_aligner->wavefront_slab,effective_lo,effective_hi);
  wf_components->mwavefronts[0]->offsets[0] = 0;
  wf_components->mwavefronts[0]->lo = 0;
  wf_components->mwavefronts[0]->hi = 0;
  // Store initial BT-piggypack element
  if (wf_components->bt_piggyback) {
    const bt_block_idx_t block_idx = wf_backtrace_buffer_init_block(wf_components->bt_buffer,0,0);
    wf_components->mwavefronts[0]->bt_pcigar[0] = 0;
    wf_components->mwavefronts[0]->bt_prev[0] = block_idx;
  }
  // Nullify unused WFs
  if (distance_metric==edit || distance_metric==gap_lineal) return;
  wf_components->d1wavefronts[0] = NULL;
  wf_components->i1wavefronts[0] = NULL;
  if (distance_metric==gap_affine) return;
  wf_components->d2wavefronts[0] = NULL;
  wf_components->i2wavefronts[0] = NULL;
}
bool wavefront_align_end2end_terminate(
    wavefront_aligner_t* const wf_aligner,
    const int score_final) {
  // Parameters
  const char* const pattern = wf_aligner->pattern;
  const int pattern_length = wf_aligner->pattern_length;
  const char* const text = wf_aligner->text;
  const int text_length = wf_aligner->text_length;
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const int alignment_k = WAVEFRONT_DIAGONAL(text_length,pattern_length);
  const wf_offset_t alignment_offset = WAVEFRONT_OFFSET(text_length,pattern_length);
  int score = score_final;
  // Check wavefront
  if (wf_components->memory_modular) score = score % wf_components->max_score_scope;
  wavefront_t* const mwavefront = wf_components->mwavefronts[score];
  if (mwavefront==NULL) return false;
  // Check limits
  wf_offset_t* const offsets = mwavefront->offsets;
  if (mwavefront->lo > alignment_k || alignment_k > mwavefront->hi) return false;
  // Check offset
  const wf_offset_t offset = offsets[alignment_k];
  if (offset < alignment_offset) return false; // end2end termination condition
  // DEBUG
  // wavefront_aligner_print(stderr,wf_aligner,0,score_final,6,0);
  // Retrieve alignment
  if (wf_aligner->alignment_scope == compute_score) {
    cigar_clear(&wf_aligner->cigar);
  } else {
    if (wf_components->bt_piggyback) {
      // Fetch backtrace from buffer and recover alignment
      wf_backtrace_buffer_recover_cigar(
          wf_components->bt_buffer,
          pattern,pattern_length,text,text_length,
          wf_aligner->match_funct,wf_aligner->match_funct_arguments,
          alignment_k,alignment_offset,
          mwavefront->bt_pcigar[alignment_k],
          mwavefront->bt_prev[alignment_k],
          &wf_aligner->cigar);
    } else {
      // Backtrace alignment
      wavefront_backtrace_affine(wf_aligner,score_final,alignment_k,alignment_offset);
    }
  }
  // Set score & finish
  wf_aligner->cigar.score = -score_final;
  return true;
}
/*
 * Ends-free alignment (Semiglobal/glocal/...)
 */
void wavefront_align_endsfree_initialize(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  const int text_begin_free = wf_aligner->alignment_form.text_begin_free;
  const int pattern_begin_free = wf_aligner->alignment_form.pattern_begin_free;
  const int max_score_scope = wf_components->max_score_scope;
  // Init wavefront zero
  const int effective_lo = -pattern_begin_free - (max_score_scope+1);
  const int effective_hi = text_begin_free + (max_score_scope+1);
  wf_components->mwavefronts[0] = wavefront_slab_allocate(
      wf_aligner->wavefront_slab,effective_lo,effective_hi);
  wf_components->mwavefronts[0]->offsets[0] = 0;
  wf_components->mwavefronts[0]->lo = -pattern_begin_free;
  wf_components->mwavefronts[0]->hi = text_begin_free;
  // Store initial BT-piggypack element
  if (wf_components->bt_piggyback) {
    const bt_block_idx_t block_idx = wf_backtrace_buffer_init_block(wf_components->bt_buffer,0,0);
    wf_components->mwavefronts[0]->bt_pcigar[0] = 0;
    wf_components->mwavefronts[0]->bt_prev[0] = block_idx;
  }
  // Init text begin-free
  int h;
  for (h=1;h<=text_begin_free;++h) {
    const int k = WAVEFRONT_DIAGONAL(h,0);
    wf_components->mwavefronts[0]->offsets[k] = WAVEFRONT_OFFSET(h,0);
    if (wf_components->bt_piggyback) {
      const bt_block_idx_t block_idx = wf_backtrace_buffer_init_block(wf_components->bt_buffer,0,h);
      wf_components->mwavefronts[0]->bt_pcigar[k] = 0;
      wf_components->mwavefronts[0]->bt_prev[k] = block_idx;
    }
  }
  // Init pattern begin-free
  int v;
  for (v=1;v<=pattern_begin_free;++v) {
    const int k = WAVEFRONT_DIAGONAL(0,v);
    wf_components->mwavefronts[0]->offsets[k] = WAVEFRONT_OFFSET(0,v);
    if (wf_components->bt_piggyback) {
      const bt_block_idx_t block_idx = wf_backtrace_buffer_init_block(wf_components->bt_buffer,v,0);
      wf_components->mwavefronts[0]->bt_pcigar[k] = 0;
      wf_components->mwavefronts[0]->bt_prev[k] = block_idx;
    }
  }
  // Nullify unused WFs
  if (distance_metric==edit || distance_metric==gap_lineal) return;
  wf_components->d1wavefronts[0] = NULL;
  wf_components->i1wavefronts[0] = NULL;
  if (distance_metric==gap_affine) return;
  wf_components->d2wavefronts[0] = NULL;
  wf_components->i2wavefronts[0] = NULL;
}
bool wavefront_align_endsfree_terminate(
    wavefront_aligner_t* const wf_aligner,
    const int score_final) {
  // Parameters
  const char* const pattern = wf_aligner->pattern;
  const int pattern_length = wf_aligner->pattern_length;
  const char* const text = wf_aligner->text;
  const int text_length = wf_aligner->text_length;
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  int score = score_final;
  // Check wavefront
  if (wf_components->memory_modular) score = score % wf_components->max_score_scope;
  wavefront_t* const mwavefront = wf_components->mwavefronts[score];
  if (mwavefront==NULL) return false;
  // Check end reached
  if (mwavefront->k_alignment_end==WAVEFRONT_DIAGONAL_NULL) return false;
  // DEBUG
  //wavefront_aligner_print(stderr,wf_aligner,0,score_final,6,0);
  // Retrieve alignment
  if (wf_aligner->alignment_scope == compute_score) {
    cigar_clear(&wf_aligner->cigar);
  } else {
    const int alignment_k = mwavefront->k_alignment_end;
    const wf_offset_t alignment_offset = mwavefront->offsets[alignment_k];
    if (wf_components->bt_piggyback) {
      // Fetch backtrace from buffer and recover alignment
      wf_backtrace_buffer_recover_cigar(
          wf_components->bt_buffer,
          pattern,pattern_length,text,text_length,
          wf_aligner->match_funct,wf_aligner->match_funct,
          alignment_k,alignment_offset,
          mwavefront->bt_pcigar[alignment_k],
          mwavefront->bt_prev[alignment_k],
          &wf_aligner->cigar);
    } else {
      // Backtrace alignment
      wavefront_backtrace_affine(wf_aligner,score,alignment_k,alignment_offset);
    }
  }
  // Set score & finish
  wf_aligner->cigar.score = -score_final;
  return true;
}
/*
 * General Alignment
 */
int wavefront_align_sequences(
    wavefront_aligner_t* const wf_aligner,
    void (*wavefront_align_initialize)(wavefront_aligner_t*),
    bool (*wavefront_align_terminate)(wavefront_aligner_t* const,const int),
    void (*wavefront_align_compute)(wavefront_aligner_t* const,const int),
    void (*wavefront_align_extend)(wavefront_aligner_t* const,const int)) {
  // Parameters
  char* const pattern = wf_aligner->pattern;
  char* const text = wf_aligner->text;
  const bool plot = wf_aligner->plot_params.plot_enabled;
  // Initialize wavefront
  (*wavefront_align_initialize)(wf_aligner);
  // PROFILE
  if (plot) wavefront_plot(wf_aligner,pattern,text,0);
  // Compute wavefronts of increasing score
  int score = 0;
  while (true) {
    // Exact extend s-wavefront
    (*wavefront_align_extend)(wf_aligner,score);
    // Check termination condition
    if ((*wavefront_align_terminate)(wf_aligner,score)) break;
    // Compute (s+1)-wavefront
    ++score;
    (*wavefront_align_compute)(wf_aligner,score);
    // Probe limits
    if (score >= wf_aligner->alignment_form.max_alignment_score) {
      wf_aligner->cigar.score = wf_aligner->alignment_form.max_alignment_score;
      return WF_ALIGN_MAX_SCORE;
    }
    if (wavefront_align_reached_limits(wf_aligner,score)) {
      return WF_ALIGN_OOM;
    }
    // PROFILE
    if (plot) wavefront_plot(wf_aligner,pattern,text,score);
    // DEBUG
    //wavefront_aligner_print(stderr,wf_aligner,0,score,6,0);
  }
  // Return OK
  return WF_ALIGN_SUCCESSFUL;
}
/*
 * Wavefront Alignment
 */
int wavefront_align(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  // DEBUG
  profiler_timer_t timer;
  if (wf_aligner->system.verbose >= 2) {
    timer_reset(&timer);
    timer_start(&timer);
    if (wf_aligner->system.verbose >= 3) {
      wavefront_report_verbose_begin(stderr,wf_aligner,pattern,pattern_length,text,text_length);
    }
  }
  // Resize wavefront aligner
  wavefront_aligner_resize(wf_aligner,pattern_length,text_length);
  // Init padded strings
  strings_padded_t* const sequences = NULL;
  if (wf_aligner->match_funct == NULL) {
    // Set sequences
    strings_padded_t* const sequences =
        strings_padded_new_rhomb(
            pattern,pattern_length,text,text_length,
            WAVEFRONT_PADDING,wf_aligner->mm_allocator);
    wf_aligner->pattern = sequences->pattern_padded;
    wf_aligner->text = sequences->text_padded;
  } else {
    wf_aligner->pattern = NULL;
    wf_aligner->text = NULL;
  }
  // Wavefront functions
  void (*wavefront_align_initialize)(wavefront_aligner_t*);
  bool (*wavefront_align_terminate)(wavefront_aligner_t* const,const int);
  void (*wavefront_align_compute)(wavefront_aligner_t* const,const int);
  void (*wavefront_align_extend)(wavefront_aligner_t* const,const int);
  // Select wavefront functions
  switch (wf_aligner->penalties.distance_metric) {
    case gap_affine: wavefront_align_compute = &wavefront_compute_affine; break;
    case gap_affine_2p: wavefront_align_compute = &wavefront_compute_affine2p; break;
    default:
      fprintf(stderr,"[WFA] Distance function not implemented yet\n"); exit(1);
      break;
  }
  const bool end2end = (wf_aligner->alignment_form.span == alignment_end2end);
  if (end2end) {
    wavefront_align_initialize = &wavefront_align_end2end_initialize;
    wavefront_align_terminate = &wavefront_align_end2end_terminate;
    wavefront_align_extend = &wavefront_extend_end2end;
  } else {
    wavefront_form_endsfree_check(wf_aligner,pattern_length,text_length);
    wavefront_align_initialize = &wavefront_align_endsfree_initialize;
    wavefront_align_terminate = &wavefront_align_endsfree_terminate;
    wavefront_align_extend = &wavefront_extend_endsfree;
  }
  if (wf_aligner->match_funct != NULL) {
    wavefront_align_extend = &wavefront_extend_custom;
  }
  // Wavefront align sequences
  const int wf_status = wavefront_align_sequences(
      wf_aligner,
      wavefront_align_initialize,
      wavefront_align_terminate,
      wavefront_align_compute,
      wavefront_align_extend);
  // Free padded strings
  if (sequences!=NULL) strings_padded_delete(sequences);
  // Reap if maximum resident memory is reached
  const uint64_t wf_memory_used = wavefront_aligner_get_size(wf_aligner);
  const bool reap_memory = wf_memory_used > wf_aligner->system.max_memory_resident;
  if (reap_memory) wavefront_aligner_reap(wf_aligner);
  // DEBUG
  if (wf_aligner->system.verbose >= 2) {
    timer_stop(&timer);
    if (wf_aligner->system.verbose == 2) {
      wavefront_report_lite(stderr,wf_aligner,
          pattern,pattern_length,text,text_length,
          wf_status,wf_memory_used,&timer);
    } else {
      wavefront_report_verbose_end(stderr,wf_aligner,wf_status,wf_memory_used,&timer);
    }
  }
  // Return
  return wf_status;
}

