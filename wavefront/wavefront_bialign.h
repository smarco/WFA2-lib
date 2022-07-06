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
 */

#ifndef WAVEFRONT_WAVEFRONT_BIALIGN_H_
#define WAVEFRONT_WAVEFRONT_BIALIGN_H_

#include "utils/commons.h"
#include "alignment/affine_penalties.h"
#include "alignment/cigar.h"
#include "wavefront_offset.h"
#include "wavefront_attributes.h"

// Wavefront ahead definition
typedef struct _wavefront_aligner_t wavefront_aligner_t;

typedef struct {
  // Scores
  int score;                      // Score total
  int score_forward;              // Score (forward)
  int score_reverse;              // Score (reverse)
  // Location
  int k_forward;                  // Breakpoint diagonal (forward)
  int k_reverse;                  // Breakpoint diagonal (reverse)
  wf_offset_t offset_forward;     // Offset (forward)
  wf_offset_t offset_reverse;     // Offset (reverse)
  affine2p_matrix_type component; // Component (M/I/D)
} wf_bialign_breakpoint_t;

/*
 * Bidirectional WFA
 */
void wavefront_bialign(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    alignment_form_t* const form,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int score_remaining,
    cigar_t* const cigar,
    const int rlevel);

#endif /* WAVEFRONT_WAVEFRONT_BIALIGN_H_ */
