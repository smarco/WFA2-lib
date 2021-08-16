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
 * DESCRIPTION: Support functions for wavefront reduction strategies
 */

#ifndef WAVEFRONT_REDUCTION_H_
#define WAVEFRONT_REDUCTION_H_

#include "utils/commons.h"

// Wavefront ahead definition
typedef struct _wavefront_aligner_t wavefront_aligner_t;

/*
 * Wavefront Reduction
 */
typedef enum {
  wavefront_reduction_none,
  wavefront_reduction_adaptive,
} wavefront_reduction_type;
typedef struct {
  wavefront_reduction_type reduction_strategy; // Reduction strategy
  int min_wavefront_length;                    // Adaptive: Minimum wavefronts length to reduce
  int max_distance_threshold;                  // Adaptive: Maximum distance between offsets allowed
} wavefront_reduction_t;

/*
 * Setup
 */
void wavefront_reduction_set_none(
    wavefront_reduction_t* const wavefront_reduction);
void wavefront_reduction_set_adaptive(
    wavefront_reduction_t* const wavefront_reduction,
    const int min_wavefront_length,
    const int max_distance_threshold);

/*
 * Reduce wavefront
 */
void wavefront_reduce(
    wavefront_aligner_t* const wf_aligner,
    const int score);

#endif /* WAVEFRONT_REDUCTION_H_ */
