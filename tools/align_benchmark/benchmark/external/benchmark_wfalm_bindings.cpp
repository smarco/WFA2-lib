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
 * DESCRIPTION: WFAlm library C++/C bridge
 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <bits/stdint-intn.h>
#include <string>
#include <vector>
#include <set>
#include <tuple>

#include "external/wfalm/wfa_lm.hpp"

/*
 * Benchmark SeqAn Adapt CIGAR
 */
void benchmark_wfalm_adapt_cigar(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    std::vector<wfalm::CIGAROp> cigar,
    char* const cigar_operations,
    int* const num_cigar_operations) {
  int pattern_pos = 0;
  int text_pos = 0;
  int cigar_idx = 0, i;
  for (auto c : cigar) {
    switch (c.op) {
      case 'D':
        text_pos += c.len;
        for (i=0;i<c.len;++i) cigar_operations[cigar_idx++] = 'I';
        break;
      case 'I':
        pattern_pos += c.len;
        for (i=0;i<c.len;++i) cigar_operations[cigar_idx++] = 'D';
        break;
      case 'M':
        for (i=0;i<c.len;++i,++text_pos,++pattern_pos) {
          if (pattern[pattern_pos] != text[text_pos]) {
            cigar_operations[cigar_idx++] = 'X';
          } else {
            cigar_operations[cigar_idx++] = 'M';
          }
        }
        break;
      default:
        break;
    }
  }
  // Return
  *num_cigar_operations = cigar_idx;
}
/*
 * Benchmark SeqAn
 */
extern "C" void benchmark_wfalm_bridge_global_affine(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int mismatch,
    const int gap_opening,
    const int gap_extension,
    char* const cigar_operations,
    int* const num_cigar_operations) {
  // Configure
  std::string patternStr(pattern,pattern_length);
  std::string textStr(text,text_length);
  int prune = -1;
  wfalm::WFScores scores(mismatch,gap_opening,gap_extension);
  // Align
  std::vector<wfalm::CIGAROp> cigar;
  int32_t score;
  std::tie(cigar,score) = wfalm::wavefront_align(patternStr,textStr,scores,prune);
  // Adapt CIGAR
  benchmark_wfalm_adapt_cigar(
      pattern,pattern_length,text,text_length,
      cigar,cigar_operations,num_cigar_operations);
}
extern "C" void benchmark_wfalm_bridge_global_affine_lowmem(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int mismatch,
    const int gap_opening,
    const int gap_extension,
    char* const cigar_operations,
    int* const num_cigar_operations) {
  // Configure
  std::string patternStr(pattern,pattern_length);
  std::string textStr(text,text_length);
  int prune = -1;
  wfalm::WFScores scores(mismatch,gap_opening,gap_extension);
  // Align
  std::vector<wfalm::CIGAROp> cigar;
  int32_t score;
  std::tie(cigar,score) = wfalm::wavefront_align_low_mem(patternStr,textStr,scores,prune);
  // Adapt CIGAR
  benchmark_wfalm_adapt_cigar(
      pattern,pattern_length,text,text_length,
      cigar,cigar_operations,num_cigar_operations);
}


