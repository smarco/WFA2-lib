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

#ifndef CIGAR_UTILS_H_
#define CIGAR_UTILS_H_

#include "cigar.h"

/*
 * Compare & Copy
 */
int cigar_cmp(
    const cigar_t* const cigar_a,
    const cigar_t* const cigar_b);
void cigar_copy(
    cigar_t* const cigar_dst,
    const cigar_t* const cigar_src);

/*
 * Mismatch Discovery
 */
void cigar_discover_mismatches(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar);

/*
 * Maxtrim
 *   Reduce the CIGAR to the maximal scoring sequence, starting from
 *   the beginning, under a given distance function
 */
bool cigar_maxtrim_gap_linear(
    cigar_t* const cigar,
    const linear_penalties_t* const penalties);
bool cigar_maxtrim_gap_affine(
    cigar_t* const cigar,
    const affine_penalties_t* const penalties);
bool cigar_maxtrim_gap_affine2p(
    cigar_t* const cigar,
    const affine2p_penalties_t* const penalties);

/*
 * Local Alignment Extraction
 */
bool cigar_maxlocal_gap_affine2p(
    cigar_t* const cigar,
    const affine2p_penalties_t* const penalties,
    const int pattern_length,
    const int text_length);

#endif /* CIGAR_UTILS_H_ */
