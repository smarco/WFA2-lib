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
#define BT_BUFFER_SEGMENT_LENGTH BUFFER_SIZE_1M

#define BT_BUFFER_IDX(segment_idx,segment_offset) ((segment_idx)*BT_BUFFER_SEGMENT_LENGTH) + (segment_offset)

/*
 * BT-Block Segments
 */
void wf_backtrace_buffer_segment_add(
    wf_backtrace_buffer_t* const bt_buffer) {
  wf_backtrace_block_t* const bt_segment = mm_allocator_calloc(
      bt_buffer->mm_allocator,BT_BUFFER_SEGMENT_LENGTH,wf_backtrace_block_t,false);
  vector_insert(bt_buffer->segments,bt_segment,wf_backtrace_block_t*);
}
void wf_backtrace_buffer_segment_reserve(
    wf_backtrace_buffer_t* const bt_buffer) {
  // Reset position
  bt_buffer->segment_offset = 0;
  ++(bt_buffer->segment_idx);
  // Check segments
  if (bt_buffer->segment_idx >= vector_get_used(bt_buffer->segments)) {
    // Check segment position
    const uint64_t block_idx = ((uint64_t)bt_buffer->segment_idx+1) * BT_BUFFER_SEGMENT_LENGTH;
    if (block_idx >= WF_BTBLOCK_IDX_NULL) {
      fprintf(stderr,"[WFA::BacktraceBuffer] Reached maximum addressable index"); exit(-1);
    }
    // Add segment
    wf_backtrace_buffer_segment_add(bt_buffer);
  }
  // Set pointer to next block free
  wf_backtrace_block_t** const segments =
      vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  bt_buffer->block_next = segments[bt_buffer->segment_idx];
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
  bt_buffer->segment_offset = 0;
  bt_buffer->segments = vector_new(10,wf_backtrace_block_t*);
  wf_backtrace_buffer_segment_add(bt_buffer); // Add initial segment
  bt_buffer->alignment_init_pos = vector_new(100,wf_backtrace_init_pos_t);
  bt_buffer->alignment_packed = vector_new(100,pcigar_t);
  bt_buffer->block_next = vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*)[0];
  // Return
  return bt_buffer;
}
void wf_backtrace_buffer_clear(
    wf_backtrace_buffer_t* const bt_buffer) {
  bt_buffer->segment_idx = 0;
  bt_buffer->segment_offset = 0;
  bt_buffer->block_next = vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*)[0];
  vector_clear(bt_buffer->alignment_init_pos);
}
void wf_backtrace_buffer_reap(
    wf_backtrace_buffer_t* const bt_buffer) {
  // Reap segments beyond the first
  const int num_segments = vector_get_used(bt_buffer->segments);
  wf_backtrace_block_t** const segments =
      vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  int i;
  for (i=1;i<num_segments;++i) {
    mm_allocator_free(bt_buffer->mm_allocator,segments[i]);
  }
  vector_set_used(bt_buffer->segments,1);
  // Clear
  bt_buffer->segment_idx = 0;
  bt_buffer->segment_offset = 0;
  bt_buffer->block_next = vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*)[0];
}
void wf_backtrace_buffer_delete(
    wf_backtrace_buffer_t* const bt_buffer) {
  // Free segments
  const int num_segments = vector_get_used(bt_buffer->segments);
  wf_backtrace_block_t** const segments =
      vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  int i;
  for (i=0;i<num_segments;++i) {
    mm_allocator_free(bt_buffer->mm_allocator,segments[i]);
  }
  // Free handlers
  vector_delete(bt_buffer->segments);
  vector_delete(bt_buffer->alignment_init_pos);
  vector_delete(bt_buffer->alignment_packed);
  mm_allocator_free(bt_buffer->mm_allocator,bt_buffer);
}
/*
 * Accessors
 */
wf_backtrace_block_t* wf_backtrace_buffer_get_block(
    wf_backtrace_buffer_t* const bt_buffer,
    const block_idx_t block_idx) {
  // Compute location
  const int segment_idx = block_idx / BT_BUFFER_SEGMENT_LENGTH;
  const int segment_offset = block_idx % BT_BUFFER_SEGMENT_LENGTH;
  wf_backtrace_block_t** const segments =
      vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  return &(segments[segment_idx][segment_offset]);
}
void wf_backtrace_buffer_add_used(
    wf_backtrace_buffer_t* const bt_buffer,
    const int used) {
  // Next
  bt_buffer->segment_offset += used;
  bt_buffer->block_next += used;
  // Reserve
  if (bt_buffer->segment_offset >= BT_BUFFER_SEGMENT_LENGTH) {
    wf_backtrace_buffer_segment_reserve(bt_buffer);
  }
}
block_idx_t wf_backtrace_buffer_get_mem(
    wf_backtrace_buffer_t* const bt_buffer,
    wf_backtrace_block_t** const bt_block_mem,
    int* const bt_blocks_available) {
  // Parameters
  const int segment_idx = bt_buffer->segment_idx;
  const int segment_offset = bt_buffer->segment_offset;
  // Get total available blocks
  *bt_block_mem = bt_buffer->block_next;
  *bt_blocks_available = BT_BUFFER_SEGMENT_LENGTH - bt_buffer->segment_offset;
  // Return current global position
  return BT_BUFFER_IDX(segment_idx,segment_offset);
}
/*
 * Store blocks
 */
void wf_backtrace_buffer_store_block(
    wf_backtrace_buffer_t* const bt_buffer,
    const pcigar_t pcigar,
    const block_idx_t prev_idx) {
  // Store BT-block
  bt_buffer->block_next->pcigar = pcigar;
  bt_buffer->block_next->prev_idx = prev_idx;
  // Next
  ++(bt_buffer->block_next);
  ++(bt_buffer->segment_offset);
  // Reserve
  if (bt_buffer->segment_offset >= BT_BUFFER_SEGMENT_LENGTH) {
    wf_backtrace_buffer_segment_reserve(bt_buffer);
  }
}
void wf_backtrace_buffer_store_block_bt(
    wf_backtrace_buffer_t* const bt_buffer,
    pcigar_t* const pcigar,
    block_idx_t* const prev_idx) {
  // Parameters
  const int segment_idx = bt_buffer->segment_idx;
  const int segment_offset = bt_buffer->segment_offset;
  // Store BT-block
  wf_backtrace_buffer_store_block(bt_buffer,*pcigar,*prev_idx);
  // Reset pcigar & set current position
  *pcigar = 0;
  *prev_idx = BT_BUFFER_IDX(segment_idx,segment_offset);
}
void wf_backtrace_buffer_store_block_init(
    wf_backtrace_buffer_t* const bt_buffer,
    const int v,
    const int h,
    pcigar_t* const pcigar,
    block_idx_t* const prev_idx) {
  // Parameters
  const int segment_idx = bt_buffer->segment_idx;
  const int segment_offset = bt_buffer->segment_offset;
  // Store initial position (v,h)
  const int init_position_offset = vector_get_used(bt_buffer->alignment_init_pos);
  wf_backtrace_init_pos_t init_pos = { .v = v, .h = h };
  vector_insert(bt_buffer->alignment_init_pos,init_pos,wf_backtrace_init_pos_t);
  // Store BT-block (Index to initial position,NULL prev)
  wf_backtrace_buffer_store_block(bt_buffer,init_position_offset,WF_BTBLOCK_IDX_NULL);
  // Set current position
  *pcigar = 0;
  *prev_idx = BT_BUFFER_IDX(segment_idx,segment_offset);
}
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
    cigar_t* const cigar) {
  // Clear temporal buffer
  vector_t* const alignment_packed = bt_buffer->alignment_packed;
  vector_clear(alignment_packed);
  // Traverse-back the BT-blocks and store all the pcigars
  wf_backtrace_block_t bt_block_last = {
      .pcigar = pcigar_last,
      .prev_idx = prev_idx_last
  };
  wf_backtrace_block_t* bt_block = &bt_block_last;
  while (bt_block->prev_idx != WF_BTBLOCK_IDX_NULL) {
    vector_insert(alignment_packed,bt_block->pcigar,pcigar_t);
    bt_block = wf_backtrace_buffer_get_block(bt_buffer,bt_block->prev_idx);
  }
  // Clear cigar
  char* cigar_buffer = cigar->operations;
  cigar->begin_offset = 0;
  // Fetch initial coordinate
  const int init_position_offset = bt_block->pcigar;
  wf_backtrace_init_pos_t* const backtrace_init_pos =
      vector_get_elm(bt_buffer->alignment_init_pos,init_position_offset,wf_backtrace_init_pos_t);
  // Add init insertions/deletions
  int i;
  int v = backtrace_init_pos->v;
  int h = backtrace_init_pos->h;
  for (i=0;i<h;++i) {*cigar_buffer = 'I'; ++cigar_buffer;};
  for (i=0;i<v;++i) {*cigar_buffer = 'D'; ++cigar_buffer;};
  // Traverse-forward the pcigars and recover the cigar
  const int num_palignment_blocks = vector_get_used(alignment_packed);
  pcigar_t* const palignment_blocks = vector_get_mem(alignment_packed,pcigar_t);
  affine_matrix_type current_matrix_type = affine_matrix_M;
  for (i=num_palignment_blocks-1;i>=0;--i) {
    // Recover block
    int cigar_block_length = 0;
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
  while (h < text_length) {*cigar_buffer = 'I'; ++cigar_buffer; ++h;};
  while (v < pattern_length) {*cigar_buffer = 'D'; ++cigar_buffer; ++v;};
  // Close CIGAR
  *cigar_buffer = '\0';
  cigar->end_offset = cigar_buffer - cigar->operations;
}
/*
 * Compact
 */
void wf_backtrace_buffer_mark_backtrace(
    wf_backtrace_buffer_t* const bt_buffer,
    const block_idx_t bt_block_idx,
    bitmap_t* const bitmap) {
  // Traverse-back the BT-blocks while not marked
  wf_backtrace_block_t bt_block_last = { .prev_idx = bt_block_idx };
  wf_backtrace_block_t* bt_block = &bt_block_last;
  // Check marked and fetch previous (until already marked or NULL is found)
  while (bt_block->prev_idx!=WF_BTBLOCK_IDX_NULL &&
        !bitmap_check__set(bitmap,bt_block->prev_idx)) {
    // Fetch previous BT-block
    const block_idx_t prev_idx = bt_block->prev_idx;
    bt_block = wf_backtrace_buffer_get_block(bt_buffer,prev_idx);
  }
}
void wf_backtrace_buffer_compact_marked(
    wf_backtrace_buffer_t* const bt_buffer,
    bitmap_t* const bitmap,
    const bool verbose) {
  // Parameters
  const int num_segments = vector_get_used(bt_buffer->segments);
  wf_backtrace_block_t** const segments = vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  // Sentinels
  block_idx_t read_segidx = 0, read_offset = 0, read_global_pos = 0;
  block_idx_t write_segidx = 0, write_offset = 0, write_global_pos = 0;
  wf_backtrace_block_t* read_block = segments[0];
  wf_backtrace_block_t* write_block = segments[0];
  // Traverse all BT-blocks from the beginning (stored marked)
  const block_idx_t max_block_idx = BT_BUFFER_IDX(bt_buffer->segment_idx,bt_buffer->segment_offset);
  while (read_global_pos < max_block_idx) {
    // Check marked block
    if (bitmap_is_set(bitmap,read_global_pos)) {
      // Store pcigar in compacted BT-buffer
      write_block->pcigar = read_block->pcigar;
      // Translate and store index in compacted BT-buffer
      write_block->prev_idx = (read_block->prev_idx==WF_BTBLOCK_IDX_NULL) ?
          WF_BTBLOCK_IDX_NULL : bitmap_erank(bitmap,read_block->prev_idx);
      // Next write
      ++write_offset; ++write_block; ++write_global_pos;
      if (write_offset >= BT_BUFFER_SEGMENT_LENGTH) {
        // Next segment
        write_block = segments[++write_segidx];
        write_offset = 0;
      }
    }
    // Next read
    ++read_offset; ++read_block; ++read_global_pos;
    if (read_offset >= BT_BUFFER_SEGMENT_LENGTH) {
      // Next segment
      if (++read_segidx >= num_segments) break;
      read_block = segments[read_segidx];
      read_offset = 0;
    }
  }
  // Update next BT-buffer index
  bt_buffer->segment_offset = write_offset;
  bt_buffer->segment_idx = write_segidx;
  bt_buffer->block_next = write_block;
  // DEBUG
  if (verbose) {
    fprintf(stderr,"[WFA::BacktraceBuffer] Compacted from %lu MB to %lu MB (%2.2f%%)\n",
        CONVERT_B_TO_MB(read_global_pos*sizeof(wf_backtrace_block_t)),
        CONVERT_B_TO_MB(write_global_pos*sizeof(wf_backtrace_block_t)),
        100.0f*(float)write_global_pos/(float)read_global_pos);
  }
}
/*
 * Utils
 */
uint64_t wf_backtrace_buffer_get_used(
    wf_backtrace_buffer_t* const bt_buffer) {
  const block_idx_t max_block_idx = BT_BUFFER_IDX(bt_buffer->segment_idx,bt_buffer->segment_offset);
  return max_block_idx;
}
uint64_t wf_backtrace_buffer_get_size_allocated(
    wf_backtrace_buffer_t* const bt_buffer) {
  const uint64_t segments_used = vector_get_used(bt_buffer->segments);
  return segments_used*BT_BUFFER_SEGMENT_LENGTH*sizeof(wf_backtrace_block_t);
}
uint64_t wf_backtrace_buffer_get_size_used(
    wf_backtrace_buffer_t* const bt_buffer) {
  const block_idx_t max_block_idx = BT_BUFFER_IDX(bt_buffer->segment_idx,bt_buffer->segment_offset);
  return max_block_idx*sizeof(wf_backtrace_block_t);
}



