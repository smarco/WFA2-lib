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
 * DESCRIPTION: WaveFront-Alignment module for the "extension" of exact matches
 */

#include "utils/string_padded.h"
#include "wavefront_extend.h"
#include "wavefront_align.h"
#include "wavefront_compute.h"
#include "wavefront_heuristic.h"

/*
 * Wavefront check termination (detect end of alignment)
 */
bool wavefront_extend_end2end_check_termination(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // Check wavefront limits
  wf_offset_t* const offsets = mwavefront->offsets;
  const int alignment_k = WAVEFRONT_DIAGONAL(text_length,pattern_length);
  if (mwavefront->lo > alignment_k || alignment_k > mwavefront->hi) return false; // Not done
  // Check offset
  const wf_offset_t offset = offsets[alignment_k];
  const wf_offset_t alignment_offset = WAVEFRONT_OFFSET(text_length,pattern_length);
  if (offset < alignment_offset) return false; // Not done
  // We are done
  mwavefront->k_alignment_end = alignment_k;
  return true;
}
bool wavefront_extend_endsfree_check_termination(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const wf_offset_t offset,
    const int k) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // Check ends-free reaching boundaries
  const int h_pos = WAVEFRONT_H(k,offset);
  const int v_pos = WAVEFRONT_V(k,offset);
  if (h_pos >= text_length) { // Text is aligned
    // Is Pattern end-free?
    const int pattern_left = pattern_length - v_pos;
    const int pattern_end_free = wf_aligner->alignment_form.pattern_end_free;
    if (pattern_left <= pattern_end_free) {
      mwavefront->k_alignment_end = k;
      return true; // Quit (we are done)
    }
  }
  if (v_pos >= pattern_length) { // Pattern is aligned
    // Is text end-free?
    const int text_left = text_length - h_pos;
    const int text_end_free = wf_aligner->alignment_form.text_end_free;
    if (text_left <= text_end_free) {
      mwavefront->k_alignment_end = k;
      return true; // Quit (we are done)
    }
  }
  // Not done
  return false;
}
/*
 * Wavefront offset extension comparing characters
 */
bool wavefront_extend_matches_packed(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const bool endsfree) {
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront==NULL) return false;
  // Extend diagonally each wavefront point
  wf_offset_t* const offsets = mwavefront->offsets;
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  int k;
  for (k=lo;k<=hi;++k) {
    // Check offset & positions
    //   - No offset should be out of boundaries !(h>tlen,v>plen)
    //   - if (h==tlen,v==plen) extension won't increment (sentinels)
    wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Fetch pattern/text blocks
    uint64_t* pattern_blocks = (uint64_t*)(wf_aligner->pattern+WAVEFRONT_V(k,offset));
    uint64_t* text_blocks = (uint64_t*)(wf_aligner->text+WAVEFRONT_H(k,offset));
    // Compare 64-bits blocks
    uint64_t cmp = *pattern_blocks ^ *text_blocks;
    while (__builtin_expect(cmp==0,0)) {
      // Increment offset (full block)
      offset += 8;
      // Next blocks
      ++pattern_blocks;
      ++text_blocks;
      // Compare
      cmp = *pattern_blocks ^ *text_blocks;
    }
    // Count equal characters
    const int equal_right_bits = __builtin_ctzl(cmp);
    const int equal_chars = DIV_FLOOR(equal_right_bits,8);
    offset += equal_chars;
    // Update offset
    offsets[k] = offset;
    // Check ends-free reaching boundaries
    if (endsfree && wavefront_extend_endsfree_check_termination(wf_aligner,mwavefront,offset,k)) {
      return true; // Quit (we are done)
    }
  }
  // Check end-to-end finished
  if (!endsfree) {
    return wavefront_extend_end2end_check_termination(wf_aligner,mwavefront);
  }
  // Alignment not finished
  return false;
}
bool wavefront_extend_matches_custom(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const bool endsfree) {
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront==NULL) return false;
  // Parameters (custom matching function)
  alignment_match_funct_t match_funct = wf_aligner->match_funct;
  void* const func_arguments = wf_aligner->match_funct_arguments;
  // Extend diagonally each wavefront point
  wf_offset_t* const offsets = mwavefront->offsets;
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
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
    if (endsfree && wavefront_extend_endsfree_check_termination(wf_aligner,mwavefront,offset,k)) {
      return true; // Quit (we are done)
    }
  }
  // Check end-to-end finished
  if (!endsfree) {
    return wavefront_extend_end2end_check_termination(wf_aligner,mwavefront);
  }
  // Alignment not finished
  return false;
}
/*
 * Wavefront exact "extension"
 */
bool wavefront_extend_end2end(
    wavefront_aligner_t* const wf_aligner,
    int score) {
  // Modular wavefront
  if (wf_aligner->wf_components.memory_modular) score = score % wf_aligner->wf_components.max_score_scope;
  // Extend wavefront
  const bool end_reached = wavefront_extend_matches_packed(wf_aligner,score,false);
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_SUCCESSFUL;
    return true; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    const bool alignment_dropped = wavefront_heuristic_cufoff(wf_aligner,score);
    if (alignment_dropped) {
      wf_aligner->align_status.status = WF_STATUS_HEURISTICALY_DROPPED;
      return true; // Done
    }
  }
  return false; // Not done
}
bool wavefront_extend_endsfree(
    wavefront_aligner_t* const wf_aligner,
    int score) {
  // Modular wavefront
  if (wf_aligner->wf_components.memory_modular) score = score % wf_aligner->wf_components.max_score_scope;
  // Extend wavefront
  const bool end_reached = wavefront_extend_matches_packed(wf_aligner,score,true);
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_SUCCESSFUL;
    return true; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    const bool alignment_dropped = wavefront_heuristic_cufoff(wf_aligner,score);
    if (alignment_dropped) {
      wf_aligner->align_status.status = WF_STATUS_HEURISTICALY_DROPPED;
      return true; // Done
    }
  }
  return false; // Not done
}
bool wavefront_extend_custom(
    wavefront_aligner_t* const wf_aligner,
    int score) {
  // Modular wavefront
  if (wf_aligner->wf_components.memory_modular) score = score % wf_aligner->wf_components.max_score_scope;
  // Extend wavefront
  const bool endsfree = (wf_aligner->alignment_form.span == alignment_endsfree);
  const bool end_reached = wavefront_extend_matches_custom(wf_aligner,score,endsfree);
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_SUCCESSFUL;
    return true; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    const bool alignment_dropped = wavefront_heuristic_cufoff(wf_aligner,score);
    if (alignment_dropped) {
      wf_aligner->align_status.status = WF_STATUS_HEURISTICALY_DROPPED;
      return true; // Done
    }
  }
  return false; // Not done
}


