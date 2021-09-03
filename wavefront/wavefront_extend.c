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
 * DESCRIPTION: WaveFront-Alignment module for the "extension" of exact matches
 */

#include "utils/string_padded.h"
#include "wavefront_extend.h"
#include "wavefront_compute.h"
#include "wavefront_reduction.h"

/*
 * Wavefront offset extension comparing characters
 *   // TODO Avoid register spilling in x86
 */
bool wavefront_extend_packed(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const bool endsfree) {
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront==NULL) return false;
  // Extend diagonally each wavefront point
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  wf_offset_t* const offsets = mwavefront->offsets;
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  int k;
  for (k=lo;k<=hi;++k) {
    // Fetch offset & positions
    //   - No offset should be out of boundaries !(h>tlen,v>plen)
    //   - if (h==tlen,v==plen) extension won't increment (sentinels)
    wf_offset_t offset = offsets[k];
    const uint32_t h = WAVEFRONT_H(k,offset); // Make unsigned to avoid checking negative
    if (h > text_length) {
      offsets[k] = WAVEFRONT_OFFSET_NULL;
      continue;
    }
    const uint32_t v = WAVEFRONT_V(k,offset); // Make unsigned to avoid checking negative
    if (v > pattern_length) {
      offsets[k] = WAVEFRONT_OFFSET_NULL;
      continue;
    }
    // Fetch pattern/text blocks
    uint64_t* pattern_blocks = (uint64_t*)(wf_aligner->pattern+v);
    uint64_t* text_blocks = (uint64_t*)(wf_aligner->text+h);
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
    if (endsfree) {
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
    }
  }
  // No end reached
  return false;
}
/*
 * Wavefront exact "extension"
 */
void wavefront_extend_end2end(
    wavefront_aligner_t* const wf_aligner,
    int score) {
  // Modular wavefront
  if (wf_aligner->wf_components.memory_modular) score = score % wf_aligner->wf_components.max_score_scope;
  // Extend wavefront
  wavefront_extend_packed(wf_aligner,score,false);
  // Reduce wavefront adaptively
  if (wf_aligner->reduction.reduction_strategy == wavefront_reduction_adaptive) {
    wavefront_reduce(wf_aligner,score);
  }
}
void wavefront_extend_endsfree(
    wavefront_aligner_t* const wf_aligner,
    int score) {
  // Modular wavefront
  if (wf_aligner->wf_components.memory_modular) score = score % wf_aligner->wf_components.max_score_scope;
  // Extend wavefront
  const bool end_reached = wavefront_extend_packed(wf_aligner,score,true);
  // Reduce wavefront adaptively
  if (!end_reached && wf_aligner->reduction.reduction_strategy == wavefront_reduction_adaptive) {
    wavefront_reduce(wf_aligner,score);
  }
}


