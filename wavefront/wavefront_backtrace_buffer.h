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
 * DESCRIPTION: WaveFront backtrace buffer to store bactrace-blocks
 */

#ifndef WAVEFRONT_BACKTRACE_BUFFER_H_
#define WAVEFRONT_BACKTRACE_BUFFER_H_

#include "utils/commons.h"
#include "utils/vector.h"
#include "system/mm_allocator.h"
#include "alignment/cigar.h"
#include "wavefront_pcigar.h"

/*
 * Constants
 */
#define WF_BACKTRACE_PREV_NULL UINT32_MAX

/*
 * Backtrace Block
 */
typedef uint32_t block_idx_t;
typedef struct {
  pcigar_t pcigar;
  block_idx_t prev_idx;
} __attribute__((packed)) wf_backtrace_block_t;

/*
 * Backtrace Buffer
 */
typedef struct {
  // Current pointer
  int segment_idx; // Current segment idx
  int segment_pos; // Current free position within segment
  // Buffer
  vector_t* segments;   // Memory segments (wf_backtrace_block_t*)
  vector_t* palignment; // Temporal buffer to store final alignment (pcigar_t)
  // MM
  mm_allocator_t* mm_allocator;
} wf_backtrace_buffer_t;

/*
 * Setup
 */
wf_backtrace_buffer_t* wf_backtrace_buffer_new();
void wf_backtrace_buffer_clear(
    wf_backtrace_buffer_t* const bt_buffer);
void wf_backtrace_buffer_reap(
    wf_backtrace_buffer_t* const bt_buffer);
void wf_backtrace_buffer_delete(
    wf_backtrace_buffer_t* const bt_buffer);

/*
 * Store blocks
 */
void wf_backtrace_buffer_store_block(
    wf_backtrace_buffer_t* const bt_buffer,
    pcigar_t* const pcigar,
    block_idx_t* const prev_idx);
void wf_backtrace_buffer_store_starting_block(
    wf_backtrace_buffer_t* const bt_buffer,
    const int v,
    const int h,
    pcigar_t* const pcigar,
    block_idx_t* const prev_idx);

/*
 * Recover CIGAR
 */
void wf_backtrace_buffer_recover_cigar(
    wf_backtrace_buffer_t* const bt_buffer,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int alignment_k,
    const int alignment_offset,
    const pcigar_t pcigar_last,
    const block_idx_t prev_idx_last,
    cigar_t* const cigar);

/*
 * Compact
 */
void wf_backtrace_buffer_mark_backtrace(
    wf_backtrace_buffer_t* const bt_buffer,
    const block_idx_t idx_last);
void wf_backtrace_buffer_compact(
    wf_backtrace_buffer_t* const bt_buffer);

/*
 * Utils
 */
uint64_t wf_backtrace_buffer_get_size(
    wf_backtrace_buffer_t* const bt_buffer);

#endif /* WAVEFRONT_BACKTRACE_BUFFER_H_ */
