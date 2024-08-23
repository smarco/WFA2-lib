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

#ifndef WAVEFRONT_EXTEND_AVX_H_
#define WAVEFRONT_EXTEND_AVX_H_

#if __AVX2__

#include "wavefront_aligner.h"

void wavefront_extend_matches_packed_end2end_avx2(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int lo,
    const int hi);

wf_offset_t wavefront_extend_matches_packed_end2end_max_avx2(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int lo,
    const int hi);

bool wavefront_extend_matches_packed_endsfree_avx2(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi);

#if __AVX512CD__ && __AVX512VL__
void wavefront_extend_matches_packed_end2end_avx512(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int lo,
    const int hi);

wf_offset_t wavefront_extend_matches_packed_end2end_max_avx512(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int lo,
    const int hi);

bool wavefront_extend_matches_packed_endsfree_avx512(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi);
#endif
#endif // AVX2

#endif /* WAVEFRONT_EXTEND_AVX_H_ */
