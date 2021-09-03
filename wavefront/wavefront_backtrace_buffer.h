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
 * Backtrace Block
 */
typedef uint32_t block_idx_t; // Up to 2^31-1 (1 bit for marking) references (~16GB of pCIGARs)
#define WF_BACKTRACE_BLOCK_MAX_IDX   ((1ul<<31)-1)
#define WF_BACKTRACE_BLOCK_IDX_NULL  UINT32_MAX

typedef struct {
  pcigar_t pcigar;            // Packed CIGAR
  block_idx_t prev_idx;       // Index of the previous CIGAR-block
} __attribute__((packed)) wf_backtrace_block_t;

// Marked blocks
#define WF_BACKTRACE_BLOCK_MARK_MASK            (1ul<<31)
#define WF_BACKTRACE_BLOCK_MARK(block_idx)      ((block_idx) | WF_BACKTRACE_BLOCK_MARK_MASK)
#define WF_BACKTRACE_BLOCK_UNMARK(block_idx)    ((block_idx) & ~(WF_BACKTRACE_BLOCK_MARK_MASK))
#define WF_BACKTRACE_BLOCK_IS_MARKED(block_idx) (((block_idx) & WF_BACKTRACE_BLOCK_MARK_MASK)!=0)

/*
 * Backtrace Buffer
 */
typedef struct {
  // Locator
  int segment_idx;                  // Current segment idx
  int segment_pos;                  // Current free position within segment
  wf_backtrace_block_t* block_next; // Next block free
  // Buffer
  vector_t* segments;   // Memory segments (wf_backtrace_block_t*)
  vector_t* palignment; // Temporal buffer to store final alignment (pcigar_t)
  // MM
  mm_allocator_t* mm_allocator;
} wf_backtrace_buffer_t;

/*
 * Setup
 */
wf_backtrace_buffer_t* wf_backtrace_buffer_new(
    mm_allocator_t* const mm_allocator);
void wf_backtrace_buffer_clear(
    wf_backtrace_buffer_t* const bt_buffer);
void wf_backtrace_buffer_reap(
    wf_backtrace_buffer_t* const bt_buffer);
void wf_backtrace_buffer_delete(
    wf_backtrace_buffer_t* const bt_buffer);

/*
 * Accessors
 */
void wf_backtrace_buffer_add_used(
    wf_backtrace_buffer_t* const bt_buffer,
    const int used);
block_idx_t wf_backtrace_buffer_get_mem(
    wf_backtrace_buffer_t* const bt_buffer,
    wf_backtrace_block_t** const bt_block_mem,
    int* const bt_blocks_available);

/*
 * Store blocks
 */
void wf_backtrace_buffer_store_block_bt(
    wf_backtrace_buffer_t* const bt_buffer,
    pcigar_t* const pcigar,
    block_idx_t* const prev_idx);
void wf_backtrace_buffer_store_block_init(
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
    const block_idx_t bt_block_idx);
block_idx_t wf_backtrace_buffer_translate_idx(
    wf_backtrace_buffer_t* const bt_buffer,
    const block_idx_t bt_block_idx);
void wf_backtrace_buffer_compact_marked(
    wf_backtrace_buffer_t* const bt_buffer_full,
    wf_backtrace_buffer_t* const bt_buffer_compacted,
    const bool verbose);

/*
 * Utils
 */
uint64_t wf_backtrace_buffer_get_size_allocated(
    wf_backtrace_buffer_t* const bt_buffer);
uint64_t wf_backtrace_buffer_get_size_used(
    wf_backtrace_buffer_t* const bt_buffer);

#endif /* WAVEFRONT_BACKTRACE_BUFFER_H_ */
