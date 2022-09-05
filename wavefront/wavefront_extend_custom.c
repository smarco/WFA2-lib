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
 * DESCRIPTION: WFA module for the "extension" of exact matches
 */

#include "wavefront_extend_custom.h"
#include "wavefront_compute.h"
#include "wavefront_termination.h"

#ifdef WFA_PARALLEL
#include <omp.h>
#endif

/*
 * Extend Custom (kernel)
 */
bool wavefront_extend_matches_custom(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi,
    const bool endsfree) {
  // Parameters (custom matching function)
  alignment_match_funct_t match_funct = wf_aligner->match_funct;
  void* const func_arguments = wf_aligner->match_funct_arguments;
  // Extend diagonally each wavefront point
  wf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=lo;k<=hi;++k) {
    // Check offset
    wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Count equal characters
    int v = WAVEFRONT_V(k,offset);
    int h = WAVEFRONT_H(k,offset);
    while (match_funct(v,h,func_arguments)) {
      h++; v++; offset++;
    }
    // Update offset
    offsets[k] = offset;
    // Check ends-free reaching boundaries
    if (endsfree && wavefront_termination_endsfree(wf_aligner,mwavefront,score,k,offset)) {
      return true; // Quit (we are done)
    }
  }
  // Alignment not finished
  return false;
}
/*
 * Extend Custom
 */
bool wavefront_extend_custom_dispatcher(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score) {
  const bool endsfree = (wf_aligner->alignment_form.span == alignment_endsfree);
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  bool end_reached = false;
  if (num_threads == 1) {
    // Extend wavefront
    end_reached = wavefront_extend_matches_custom(wf_aligner,mwavefront,score,lo,hi,endsfree);
  } else {
#ifdef WFA_PARALLEL
    // Extend wavefront in parallel
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      wavefront_compute_thread_limits(
          omp_get_thread_num(),omp_get_num_threads(),lo,hi,&t_lo,&t_hi);
      if (wavefront_extend_matches_custom(wf_aligner,mwavefront,score,t_lo,t_hi,endsfree)) {
        end_reached = true;
      }
    }
#endif
  }
  // Return
  return end_reached;
}
int wavefront_extend_custom(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Compute score
  const bool memory_modular = wf_aligner->wf_components.memory_modular;
  const int max_score_scope = wf_aligner->wf_components.max_score_scope;
  const int score_mod = (memory_modular) ? score % max_score_scope : score;
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score_mod];
  if (mwavefront == NULL) {
    // Check alignment feasibility (for heuristic variants that can lead to no solution)
    if (wf_aligner->align_status.num_null_steps > wf_aligner->wf_components.max_score_scope) {
      wf_aligner->align_status.status = WF_STATUS_UNFEASIBLE;
      wf_aligner->align_status.score = score;
      return 1; // Done
    }
    return 0; // Not done
  }
  // Extend (dispatcher)
  bool end_reached = wavefront_extend_custom_dispatcher(wf_aligner,mwavefront,score);
  if (wf_aligner->alignment_form.span == alignment_end2end) {
    end_reached = wavefront_termination_end2end(wf_aligner,mwavefront,score,score_mod);
  }
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_END_REACHED;
    wf_aligner->align_status.score = score;
    return 1; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    wavefront_heuristic_cufoff(wf_aligner,score,score_mod);
  }
  return 0; // Not done
}


