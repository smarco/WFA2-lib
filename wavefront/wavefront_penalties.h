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
 * DESCRIPTION: WaveFront penalties handling module
 */

#ifndef WAVEFRONT_WAVEFRONT_PENALTIES_H_
#define WAVEFRONT_WAVEFRONT_PENALTIES_H_

#include "gap_lineal/lineal_penalties.h"
#include "gap_affine/affine_penalties.h"
#include "gap_affine2p/affine2p_penalties.h"

/*
 * Penalty adaptation strategy
 */
typedef enum {
  wavefronts_penalties_force_zero_match,
  wavefronts_penalties_shifted_penalties,
} wf_penalties_strategy_type;

/*
 * Wavefront Penalties
 */
typedef struct {
  // int match;          // (M = 0)
  int mismatch;          // (X > 0)
  int gap_opening1;      // (O1 > 0)
  int gap_extension1;    // (E1 > 0)
  int gap_opening2;      // (O2 > 0)
  int gap_extension2;    // (E2 > 0)
} wavefronts_penalties_t;

/*
 * Penalties adjustment
 */
void wavefronts_penalties_set_lineal(
    wavefronts_penalties_t* const wavefronts_penalties,
    lineal_penalties_t* const lineal_penalties,
    const wf_penalties_strategy_type penalties_strategy);
void wavefronts_penalties_set_affine(
    wavefronts_penalties_t* const wavefronts_penalties,
    affine_penalties_t* const affine_penalties,
    const wf_penalties_strategy_type penalties_strategy);
void wavefronts_penalties_set_affine2p(
    wavefronts_penalties_t* const wavefronts_penalties,
    affine2p_penalties_t* const affine2p_penalties,
    const wf_penalties_strategy_type penalties_strategy);

#endif /* WAVEFRONT_WAVEFRONT_PENALTIES_H_ */
