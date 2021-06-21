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

#include "wavefront_backtrace_buffer.h"

/*
 * Config
 */
#define BT_BUFFER_SEGMENT_LENGTH BUFFER_SIZE_8M

/*
 * BT-Block Segments
 */
void wf_backtrace_buffer_add_segment(
    wf_backtrace_buffer_t* const bt_buffer) {
  wf_backtrace_block_t* const bt_segment = mm_allocator_calloc(
      bt_buffer->mm_allocator,BT_BUFFER_SEGMENT_LENGTH,wf_backtrace_block_t,false);
  vector_insert(bt_buffer->segments,bt_segment,wf_backtrace_block_t*);
}
/*
 * Setup
 */
wf_backtrace_buffer_t* wf_backtrace_buffer_new(
    mm_allocator_t* const mm_allocator) {
  // Alloc
  wf_backtrace_buffer_t* const bt_buffer =
      mm_allocator_alloc(mm_allocator,wf_backtrace_buffer_t);
  bt_buffer->mm_allocator = mm_allocator;
  // Initialize
  bt_buffer->segment_idx = 0;
  bt_buffer->segment_pos = 1;  // Discard block-0
  bt_buffer->segments = vector_new(10,wf_backtrace_block_t*);
  wf_backtrace_buffer_add_segment(bt_buffer); // Add initial segment
  bt_buffer->palignment = vector_new(100,pcigar_t);
  // Return
  return bt_buffer;
}
void wf_backtrace_buffer_clear(
    wf_backtrace_buffer_t* const bt_buffer) {
  bt_buffer->segment_idx = 0;
  bt_buffer->segment_pos = 1; // Discard block-0
  vector_set_used(bt_buffer->segments,1);
}
void wf_backtrace_buffer_reap(
    wf_backtrace_buffer_t* const bt_buffer) {
  // Reap segments beyond the first
  const int num_segments = vector_get_used(bt_buffer->segments);
  wf_backtrace_block_t** const bt_blocks =
      vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  int i;
  for (i=1;i<num_segments;++i) {
    mm_allocator_free(bt_buffer->mm_allocator,bt_blocks[i]);
  }
  vector_set_used(bt_buffer->segments,1);
  // Clear
  bt_buffer->segment_idx = 0;
  bt_buffer->segment_pos = 1; // Discard block-0
}
void wf_backtrace_buffer_delete(
    wf_backtrace_buffer_t* const bt_buffer) {
  // Free segments
  const int num_segments = vector_get_used(bt_buffer->segments);
  wf_backtrace_block_t** const bt_blocks =
      vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  int i;
  for (i=0;i<num_segments;++i) {
    mm_allocator_free(bt_buffer->mm_allocator,bt_blocks[i]);
  }
  // Free handlers
  vector_delete(bt_buffer->segments);
  vector_delete(bt_buffer->palignment);
  mm_allocator_free(bt_buffer->mm_allocator,bt_buffer);
}
/*
 * Accessors
 */
wf_backtrace_block_t* wf_backtrace_buffer_get_block(
    wf_backtrace_buffer_t* const bt_buffer,
    const block_idx_t block_idx) {
  // Compute location
  const int segment_idx = block_idx/BT_BUFFER_SEGMENT_LENGTH;
  const int segment_pos = block_idx%BT_BUFFER_SEGMENT_LENGTH;
  wf_backtrace_block_t** const segments_mem = vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  return &(segments_mem[segment_idx][segment_pos]);
}
void wf_backtrace_buffer_store_block(
    wf_backtrace_buffer_t* const bt_buffer,
    pcigar_t* const pcigar,
    block_idx_t* const prev_idx) {
  // Check occupancy
  if (bt_buffer->segment_pos >= BT_BUFFER_SEGMENT_LENGTH) {
    bt_buffer->segment_pos = 0; // Reset
    ++(bt_buffer->segment_idx);
    // Check segments
    if (bt_buffer->segment_idx >= vector_get_used(bt_buffer->segments)) {
      wf_backtrace_buffer_add_segment(bt_buffer);
    }
  }
  // Parameters
  const int segment_idx = bt_buffer->segment_idx;
  const int segment_pos = bt_buffer->segment_pos;
  // Store BT-block
  wf_backtrace_block_t** const segments_mem = vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  segments_mem[segment_idx][segment_pos].pcigar = *pcigar;
  segments_mem[segment_idx][segment_pos].prev_idx = *prev_idx;
  // Reset pcigar & set current position
  *pcigar = 0;
  *prev_idx = (segment_idx*BT_BUFFER_SEGMENT_LENGTH) + segment_pos;
  // Next
  ++(bt_buffer->segment_pos);
}
/*
 * Recover CIGAR
 */
void wf_backtrace_buffer_recover_cigar(
    wf_backtrace_buffer_t* const bt_buffer,
    const pcigar_t pcigar_last,
    const block_idx_t prev_idx_last,
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    cigar_t* const cigar) {
  // Clear temporal buffer
  vector_t* const palignment = bt_buffer->palignment;
  vector_clear(palignment);
  // Traverse-back the BT-blocks and store all the pcigars
  wf_backtrace_block_t bt_block_last = {
      .pcigar = pcigar_last,
      .prev_idx = prev_idx_last
  };
  wf_backtrace_block_t* bt_block = &bt_block_last;
  vector_insert(palignment,bt_block->pcigar,pcigar_t);
  while (bt_block->prev_idx > 0) {
    bt_block = wf_backtrace_buffer_get_block(bt_buffer,bt_block->prev_idx);
    vector_insert(palignment,bt_block->pcigar,pcigar_t);
  }
  // Clear cigar
  cigar->begin_offset = 0;
  // Traverse-forward the pcigars and recover the cigar
  char* cigar_buffer = cigar->operations;
  const int num_palignment_blocks = vector_get_used(palignment);
  pcigar_t* const palignment_blocks = vector_get_mem(palignment,pcigar_t);
  affine_matrix_type current_matrix_type = affine_matrix_M;
  int h = 0, v = 0, cigar_block_length = 0, i;
  for (i=num_palignment_blocks-1;i>=0;--i) {
    // Recover block
    pcigar_recover(
        palignment_blocks[i],
        pattern,pattern_length,
        text,text_length,
        &v,&h,
        cigar_buffer,&cigar_block_length,
        &current_matrix_type);
    // Update CIGAR
    cigar_buffer += cigar_block_length;
  }
  // Account for last stroke of matches
  const int num_matches = pcigar_recover_extend(
      pattern,pattern_length,text,text_length,v,h,cigar_buffer);
  v += num_matches;
  h += num_matches;
  cigar_buffer += num_matches;
  // Account for last stroke of insertion/deletion
  while (v < pattern_length) {*cigar_buffer = 'D'; ++cigar_buffer; ++v;};
  while (h < text_length) {*cigar_buffer = 'I'; ++cigar_buffer; ++h;};
  // Close CIGAR
  cigar->end_offset = cigar_buffer - cigar->operations;
  *cigar_buffer = '\0';
}
/*
 * Utils
 */
uint64_t wf_backtrace_buffer_get_size(
    wf_backtrace_buffer_t* const bt_buffer) {
  const uint64_t segments_used = vector_get_used(bt_buffer->segments);
  return segments_used*BT_BUFFER_SEGMENT_LENGTH;
}



