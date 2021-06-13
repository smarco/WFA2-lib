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
 * DESCRIPTION: WaveFront components module
 */

#ifndef WAVEFRONT_WAVEFRONT_COMPONENTS_H_
#define WAVEFRONT_WAVEFRONT_COMPONENTS_H_

#include "wavefront_aligner.h"

/*
 * Setup
 */
void wavefront_components_allocate(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric);
void wavefront_components_clear(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric);
void wavefront_components_free(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric);

/*
 * Compute dimensions
 */
void wavefront_components_dimensions(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric,
    const int pattern_length,
    const int text_length,
    int* const max_score_scope,
    int* const num_wavefronts);

#endif /* WAVEFRONT_WAVEFRONT_COMPONENTS_H_ */
