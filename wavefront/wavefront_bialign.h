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

typedef struct {
  wavefront_aligner_t* aligner_forward;       // Forward aligner
  wavefront_aligner_t* aligner_reverse;       // Reverse aligner
  wavefront_aligner_t* aligner_subsidiary;    // Subsidiary aligner
  wf_bialign_breakpoint_t bialign_breakpoint; // Breakpoint of two wavefronts (bialigner)
} wavefront_bialigner_t;

/*
 * Setup
 */
wavefront_bialigner_t* wavefront_bialigner_new(
    wavefront_aligner_attr_t* const attributes);
void wavefront_bialigner_reap(
    wavefront_bialigner_t* const wf_bialigner);
void wavefront_bialigner_delete(
    wavefront_bialigner_t* const wf_bialigner);

/*
 * Accessors
 */
uint64_t wavefront_bialigner_get_size(
    wavefront_bialigner_t* const wf_bialigner);
void wavefront_bialigner_heuristic_inherit(
    wavefront_bialigner_t* const wf_bialigner,
    wavefront_heuristic_t* const heuristic);

/*
 * Bidirectional WFA
 */
void wavefront_bialign(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);

#endif /* WAVEFRONT_WAVEFRONT_BIALIGN_H_ */
