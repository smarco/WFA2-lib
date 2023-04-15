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
 * DESCRIPTION: SeqAn library C++/C bridge
 */

#include <iostream>
#include <seqan/align.h>
#include <seqan/align/align_base.h>

using namespace seqan;
using namespace std;

typedef String<char> TSequence;
typedef Align<TSequence,ArrayGaps> TAlign;

typedef typename Row<TAlign>::Type TRow;
typedef typename Position<TAlign>::Type TPosition;
typedef typename Iterator<typename Row<TAlign>::Type const, Standard>::Type TIter;

/*
 * Benchmark SeqAn Adapt CIGAR
 */
void benchmark_seqan_parse_cigar(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    TAlign& align,
    char* const edit_operations,
    int* const num_edit_operations) {
  // Extract Alignment Rows
  TRow& row1 = row(align,0);
  TRow& row2 = row(align,1);
  TIter row1_iter = iter(row1,0);
  TIter row2_iter = iter(row2,0);
  // Adapt CIGAR
  int cigar_pos = 0, pattern_pos = 0, text_pos = 0;
  while (pattern_pos < pattern_length ||
         text_pos < text_length) {
    if (pattern_pos < pattern_length && isGap(row2_iter)) {
      edit_operations[cigar_pos++] = 'D';
      row1_iter++; row2_iter++; pattern_pos++;
    } else if (text_pos < text_length && isGap(row1_iter)) {
      edit_operations[cigar_pos++] = 'I';
      row1_iter++; row2_iter++; text_pos++;
    } else {
      edit_operations[cigar_pos++] =
          (pattern[pattern_pos]==text[text_pos]) ? 'M' : 'X';
      pattern_pos++; row1_iter++;
      text_pos++; row2_iter++;
    }
  }
  *num_edit_operations = cigar_pos;
}
/*
 * Benchmark SeqAn
 */
extern "C" void benchmark_seqan_bridge_global_edit(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    char* const edit_operations,
    int* const num_edit_operations) {
  TSequence seq1 = pattern;
  TSequence seq2 = text;
  // Align
  TAlign align;
  resize(rows(align),2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  globalAlignment(align,Score<int,Simple>(0,-1,-1));
  // Adapt alignment CIGAR
  benchmark_seqan_parse_cigar(
      pattern,pattern_length,text,text_length,
      align,edit_operations,num_edit_operations);
}
extern "C" int benchmark_seqan_bridge_global_edit_bpm(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    char* const edit_operations,
    int* const num_edit_operations) {
  TSequence seq1 = pattern;
  TSequence seq2 = text;
  // Align
  return globalAlignmentScore(seq1,seq2,MyersBitVector());
}
extern "C" void benchmark_seqan_bridge_global_lineal(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int match,
    const int mismatch,
    const int insertion,
    const int deletion,
    char* const edit_operations,
    int* const num_edit_operations) {
  TSequence seq1 = pattern;
  TSequence seq2 = text;
  // Align
  TAlign align;
  resize(rows(align),2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  globalAlignment(align,
      Score<int,Simple>(-match,-mismatch,-insertion,-deletion),LinearGaps());
  // Adapt alignment CIGAR
  benchmark_seqan_parse_cigar(
      pattern,pattern_length,text,text_length,
      align,edit_operations,num_edit_operations);
}
extern "C" void benchmark_seqan_bridge_global_affine(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int match,
    const int mismatch,
    const int gap_opening,
    const int gap_extension,
    char* const edit_operations,
    int* const num_edit_operations) {
  TSequence seq1 = pattern;
  TSequence seq2 = text;
  // Align
  TAlign align;
  resize(rows(align),2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  globalAlignment(align,
      Score<int,Simple>(-match,-mismatch,-gap_extension,-gap_opening-gap_extension),AffineGaps());
  // Adapt alignment CIGAR
  benchmark_seqan_parse_cigar(
      pattern,pattern_length,text,text_length,
      align,edit_operations,num_edit_operations);
}


