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
 * DESCRIPTION: Basic bitmap datastructure (static)
 */

#ifndef BITMAP_H_
#define BITMAP_H_

/*
 * Includes
 */
#include "utils/commons.h"
#include "system/mm_allocator.h"

#define BITMAP_BLOCK_ELEMENTS 64
#define BITMAP_BLOCK_MASK     0x0000000000000001ul

/*
 * Bitmap
 */
typedef struct {
  uint64_t counter;
  uint64_t bitmap;
} bitmap_block_t;
typedef struct {
  // Bitmap
  uint64_t num_blocks;
  bitmap_block_t* bitmap_blocks;
  // MM
  mm_allocator_t* mm_allocator;
} bitmap_t;

/*
 * Setup
 */
bitmap_t* bitmap_new(
    const uint64_t length,
    mm_allocator_t* const mm_allocator);
void bitmap_delete(
    bitmap_t* const bitmap);

/*
 * Accessors
 */
void bitmap_set(
    bitmap_t* const bitmap,
    const uint64_t pos);
bool bitmap_is_set(
    bitmap_t* const bitmap,
    const uint64_t pos);
bool bitmap_check__set(
    bitmap_t* const bitmap,
    const uint64_t pos);

/*
 * Rank
 */
void bitmap_update_counters(
    bitmap_t* const bitmap);
uint64_t bitmap_erank(
    bitmap_t* const bitmap,
    const uint64_t pos);

#endif /* BITMAP_H_ */
