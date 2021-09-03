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

#define BT_BUFFER_IDX(segment_idx,segment_pos) ((segment_idx)*BT_BUFFER_SEGMENT_LENGTH) + (segment_pos)

#define BT_BUFFER_PCIGAR_PB(h,v)      (((uint64_t)h << 32) | (uint64_t)v)
#define BT_BUFFER_PCIGAR_PB_H(pcigar) UINT64_TO_UINT32_MSB(pcigar);
#define BT_BUFFER_PCIGAR_PB_V(pcigar) UINT64_TO_UINT32_LSB(pcigar);

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
  bt_buffer->segment_pos = 0;
  ++(bt_buffer->segment_idx);
  // Check segments
  if (bt_buffer->segment_idx >= vector_get_used(bt_buffer->segments)) {
    // Check segment position
    const block_idx_t block_idx = (bt_buffer->segment_idx+1) * BT_BUFFER_SEGMENT_LENGTH;
    if (block_idx >= WF_BACKTRACE_BLOCK_MAX_IDX) {
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
  bt_buffer->segment_pos = 0;
  bt_buffer->segments = vector_new(10,wf_backtrace_block_t*);
  wf_backtrace_buffer_segment_add(bt_buffer); // Add initial segment
  bt_buffer->palignment = vector_new(100,pcigar_t);
  bt_buffer->block_next = vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*)[0];
  // Return
  return bt_buffer;
}
void wf_backtrace_buffer_clear(
    wf_backtrace_buffer_t* const bt_buffer) {
  bt_buffer->segment_idx = 0;
  bt_buffer->segment_pos = 0;
  bt_buffer->block_next = vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*)[0];
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
  bt_buffer->segment_pos = 0;
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
  const int segment_idx = block_idx / BT_BUFFER_SEGMENT_LENGTH;
  const int segment_pos = block_idx % BT_BUFFER_SEGMENT_LENGTH;
  wf_backtrace_block_t** const segments =
      vector_get_mem(bt_buffer->segments,wf_backtrace_block_t*);
  return &(segments[segment_idx][segment_pos]);
}
void wf_backtrace_buffer_add_used(
    wf_backtrace_buffer_t* const bt_buffer,
    const int used) {
  // Next
  bt_buffer->segment_pos += used;
  bt_buffer->block_next += used;
  // Reserve
  if (bt_buffer->segment_pos >= BT_BUFFER_SEGMENT_LENGTH) {
    wf_backtrace_buffer_segment_reserve(bt_buffer);
  }
}
block_idx_t wf_backtrace_buffer_get_mem(
    wf_backtrace_buffer_t* const bt_buffer,
    wf_backtrace_block_t** const bt_block_mem,
    int* const bt_blocks_available) {
  // Parameters
  const int segment_idx = bt_buffer->segment_idx;
  const int segment_pos = bt_buffer->segment_pos;
  // Get total available blocks
  *bt_block_mem = bt_buffer->block_next;
  *bt_blocks_available = BT_BUFFER_SEGMENT_LENGTH - bt_buffer->segment_pos;
  // Return current global position
  return BT_BUFFER_IDX(segment_idx,segment_pos);
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
  ++(bt_buffer->segment_pos);
  // Reserve
  if (bt_buffer->segment_pos >= BT_BUFFER_SEGMENT_LENGTH) {
    wf_backtrace_buffer_segment_reserve(bt_buffer);
  }
}
void wf_backtrace_buffer_store_block_bt(
    wf_backtrace_buffer_t* const bt_buffer,
    pcigar_t* const pcigar,
    block_idx_t* const prev_idx) {
  // Parameters
  const int segment_idx = bt_buffer->segment_idx;
  const int segment_pos = bt_buffer->segment_pos;
  // Store BT-block
  wf_backtrace_buffer_store_block(bt_buffer,*pcigar,*prev_idx);
  // Reset pcigar & set current position
  *pcigar = 0;
  *prev_idx = BT_BUFFER_IDX(segment_idx,segment_pos);
}
void wf_backtrace_buffer_store_block_init(
    wf_backtrace_buffer_t* const bt_buffer,
    const int v,
    const int h,
    pcigar_t* const pcigar,
    block_idx_t* const prev_idx) {
  // Parameters
  const int segment_idx = bt_buffer->segment_idx;
  const int segment_pos = bt_buffer->segment_pos;
  // Store BT-block
  wf_backtrace_buffer_store_block(bt_buffer,BT_BUFFER_PCIGAR_PB(h,v),WF_BACKTRACE_BLOCK_IDX_NULL);
  // Set current position
  *pcigar = 0;
  *prev_idx = BT_BUFFER_IDX(segment_idx,segment_pos);
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
  vector_t* const palignment = bt_buffer->palignment;
  vector_clear(palignment);
  // Traverse-back the BT-blocks and store all the pcigars
  wf_backtrace_block_t bt_block_last = {
      .pcigar = pcigar_last,
      .prev_idx = prev_idx_last
  };
  wf_backtrace_block_t* bt_block = &bt_block_last;
  while (bt_block->prev_idx != WF_BACKTRACE_BLOCK_IDX_NULL) {
    vector_insert(palignment,bt_block->pcigar,pcigar_t);
    bt_block = wf_backtrace_buffer_get_block(bt_buffer,bt_block->prev_idx);
  }
  // Clear cigar
  char* cigar_buffer = cigar->operations;
  cigar->begin_offset = 0;
  // Set initial coordinate and add init insertions/deletions
  int i;
  int v = BT_BUFFER_PCIGAR_PB_V(bt_block->pcigar);
  int h = BT_BUFFER_PCIGAR_PB_H(bt_block->pcigar);
  for (i=0;i<h;++i) {*cigar_buffer = 'I'; ++cigar_buffer;};
  for (i=0;i<v;++i) {*cigar_buffer = 'D'; ++cigar_buffer;};
  // Traverse-forward the pcigars and recover the cigar
  const int num_palignment_blocks = vector_get_used(palignment);
  pcigar_t* const palignment_blocks = vector_get_mem(palignment,pcigar_t);
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
    const block_idx_t bt_block_idx) {
  // Traverse-back the BT-blocks while not marked
  wf_backtrace_block_t bt_block_last = { .prev_idx = bt_block_idx };
  wf_backtrace_block_t* bt_block = &bt_block_last;
  while (!WF_BACKTRACE_BLOCK_IS_MARKED(bt_block->prev_idx)) {
    const block_idx_t prev_idx = bt_block->prev_idx;
    bt_block->prev_idx = WF_BACKTRACE_BLOCK_MARK(prev_idx);
    bt_block = wf_backtrace_buffer_get_block(bt_buffer,prev_idx);
  }
}
block_idx_t wf_backtrace_buffer_translate_idx(
    wf_backtrace_buffer_t* const bt_buffer,
    const block_idx_t bt_block_idx) {
  // Fetch block
  wf_backtrace_block_t* const bt_block = wf_backtrace_buffer_get_block(bt_buffer,bt_block_idx);
  // Return translation
  return bt_block->prev_idx;
}
void wf_backtrace_buffer_compact_marked(
    wf_backtrace_buffer_t* const bt_buffer_full,
    wf_backtrace_buffer_t* const bt_buffer_compacted,
    const bool verbose) {
  // Parameters
  const int num_segments = vector_get_used(bt_buffer_full->segments);
  wf_backtrace_block_t** const segments =
     vector_get_mem(bt_buffer_full->segments,wf_backtrace_block_t*);
  // Sentinels
  block_idx_t read_seg = 0, read_pos = 0, read_pos_global = 0;
  block_idx_t write_pos_global = 0;
  wf_backtrace_block_t* read_block = segments[0];
  // Traverse all BT-blocks from the beginning (stored marked)
  const block_idx_t max_block_idx_full = BT_BUFFER_IDX(bt_buffer_full->segment_idx,bt_buffer_full->segment_pos);
  while (read_pos_global < max_block_idx_full) {
    // Check marked block
    if (WF_BACKTRACE_BLOCK_IS_MARKED(read_block->prev_idx)) {
      // Translate old index
      block_idx_t block_idx_tr;
      if (read_block->prev_idx==WF_BACKTRACE_BLOCK_IDX_NULL) {
        block_idx_tr = WF_BACKTRACE_BLOCK_IDX_NULL;
      } else {
        block_idx_tr = wf_backtrace_buffer_translate_idx(
            bt_buffer_full,WF_BACKTRACE_BLOCK_UNMARK(read_block->prev_idx));
      }
      // Store in compacted BT-buffer
      wf_backtrace_buffer_store_block(
          bt_buffer_compacted,read_block->pcigar,block_idx_tr);
      // Store translation in old block
      read_block->prev_idx = write_pos_global++;
    }
    // Next read
    ++read_pos; ++read_pos_global; ++read_block;
    if (read_pos >= BT_BUFFER_SEGMENT_LENGTH) {
      if (++read_seg >= num_segments) break;
      read_pos = 0; // Clear position
      read_block = segments[read_seg]; // Next segment
    }
  }
  // DEBUG
  if (verbose) {
    const block_idx_t max_block_idx_compacted =
        BT_BUFFER_IDX(bt_buffer_compacted->segment_idx,bt_buffer_compacted->segment_pos);
    fprintf(stderr,"[WFA::BacktraceBuffer] Compacted from %lu MB to %lu MB (%2.2f%%)\n",
        CONVERT_B_TO_MB(max_block_idx_full*sizeof(wf_backtrace_block_t)),
        CONVERT_B_TO_MB(max_block_idx_compacted*sizeof(wf_backtrace_block_t)),
        100.0f*(float)max_block_idx_compacted/(float)max_block_idx_full);
  }
}
/*
 * Utils
 */
uint64_t wf_backtrace_buffer_get_size_allocated(
    wf_backtrace_buffer_t* const bt_buffer) {
  const uint64_t segments_used = vector_get_used(bt_buffer->segments);
  return segments_used*BT_BUFFER_SEGMENT_LENGTH*sizeof(wf_backtrace_block_t);
}
uint64_t wf_backtrace_buffer_get_size_used(
    wf_backtrace_buffer_t* const bt_buffer) {
  const block_idx_t max_block_idx = BT_BUFFER_IDX(bt_buffer->segment_idx,bt_buffer->segment_pos);
  return max_block_idx*sizeof(wf_backtrace_block_t);
}



