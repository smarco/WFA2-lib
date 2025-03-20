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

#include "cigar_utils.h"
#include "utils/commons.h"

/*
 * Compare & Copy
 */
int cigar_cmp(
    const cigar_t* const cigar_a,
    const cigar_t* const cigar_b) {
  // Compare lengths
  const int length_cigar_a = cigar_a->end_offset - cigar_a->begin_offset;
  const int length_cigar_b = cigar_b->end_offset - cigar_b->begin_offset;
  if (length_cigar_a != length_cigar_b) return length_cigar_a - length_cigar_b;
  // Compare operations
  char* const operations_a = cigar_a->operations + cigar_a->begin_offset;
  char* const operations_b = cigar_b->operations + cigar_b->begin_offset;
  int i;
  for (i=0;i<length_cigar_a;++i) {
    if (operations_a[i] != operations_b[i]) {
      return operations_a[i] - operations_b[i];
    }
  }
  // Equal
  return 0;
}
void cigar_copy(
    cigar_t* const cigar_dst,
    const cigar_t* const cigar_src) {
  cigar_dst->max_operations = cigar_src->max_operations;
  cigar_dst->begin_offset = cigar_src->begin_offset;
  cigar_dst->end_offset = cigar_src->end_offset;
  cigar_dst->score = cigar_src->score;
  memcpy(cigar_dst->operations+cigar_src->begin_offset,
         cigar_src->operations+cigar_src->begin_offset,
         cigar_src->end_offset-cigar_src->begin_offset);
}
/*
 * Mismatch Discovery
 */
void cigar_discover_mismatches(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar) {
  // Refine adding mismatches
  int i, p=0, t=0;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    // Check limits
    if (p >= pattern_length || t >= text_length) break;
    switch (cigar->operations[i]) {
      case 'M':
        cigar->operations[i] = (pattern[p]==text[t]) ? 'M' : 'X';
        ++p; ++t;
        break;
      case 'I':
        ++t;
        break;
      case 'D':
        ++p;
        break;
      default:
        fprintf(stderr,"[CIGAR] Wrong edit operation\n");
        exit(1);
        break;
    }
  }
  while (p < pattern_length) { cigar->operations[i++] = 'D'; ++p; };
  while (t < text_length) { cigar->operations[i++] = 'I'; ++t; };
  cigar->end_offset = i;
  cigar->operations[cigar->end_offset] = '\0';
  //  // DEBUG
  //  printf("Score=%ld\nPath-length=%" PRIu64 "\nCIGAR=%s\n",
  //      gaba_alignment->score,gaba_alignment->plen,
  //      cigar->operations);
}
/*
 * Maxtrim
 *   Reduce the CIGAR to the maximal scoring sequence, starting from
 *   the beginning, under a given distance function
 */
bool cigar_maxtrim_gap_linear(
    cigar_t* const cigar,
    const linear_penalties_t* const penalties) {
  // Parameters
  const char* const operations = cigar->operations;
  const int begin_offset = cigar->begin_offset;
  const int end_offset = cigar->end_offset;
  const int match_score = (penalties->match!=0) ? penalties->match : -1;
  // Max-score
  int max_score = 0, max_score_offset = begin_offset, max_end_v = 0, max_end_h = 0;
  // Traverse all cigar
  int score = 0, end_v = 0, end_h = 0, i;
  for (i=begin_offset;i<end_offset;++i) {
    // Update score
    switch (operations[i]) {
      case 'M':
        score -= match_score;
        ++end_v; ++end_h;
        break;
      case 'X':
        score -= penalties->mismatch;
        ++end_v; ++end_h;
        break;
      case 'I':
        score -= penalties->indel;
        ++end_h;
        break;
      case 'D':
        score -= penalties->indel;
        ++end_v;
        break;
    }
    // Compare max
    if (max_score < score) {
      max_score = score;
      max_score_offset = i;
      max_end_v = end_v;
      max_end_h = end_h;
    }
  }
  // Keep the max-scoring part of the cigar
  const bool cigar_trimmed = (max_score_offset != end_offset-1);
  if (max_score == 0) {
    cigar_clear(cigar);
  } else {
    cigar->operations[max_score_offset+1] = '\0';
    cigar->end_offset = max_score_offset + 1;
    cigar->score = max_score;
    cigar->end_v = max_end_v;
    cigar->end_h = max_end_h;
  }
  // Return
  return cigar_trimmed;
}
bool cigar_maxtrim_gap_affine(
    cigar_t* const cigar,
    const affine_penalties_t* const penalties) {
  // Parameters
  const char* const operations = cigar->operations;
  const int begin_offset = cigar->begin_offset;
  const int end_offset = cigar->end_offset;
  const int match_score = (penalties->match!=0) ? penalties->match : -1;
  // Max-score
  int max_score = 0, max_score_offset = begin_offset, max_end_v = 0, max_end_h = 0;
  // Traverse all cigar
  char last_op = '\0';
  int score = 0, end_v = 0, end_h = 0, i;
  for (i=begin_offset;i<end_offset;++i) {
    // Update score
    switch (operations[i]) {
      case 'M':
        score -= match_score;
        ++end_v; ++end_h;
        break;
      case 'X':
        score -= penalties->mismatch;
        ++end_v; ++end_h;
        break;
      case 'I':
        score -= penalties->gap_extension + ((last_op=='I') ? 0 : penalties->gap_opening);
        ++end_h;
        break;
      case 'D':
        score -= penalties->gap_extension + ((last_op=='D') ? 0 : penalties->gap_opening);
        ++end_v;
        break;
    }
    last_op = operations[i];
    // Compare max
    if (max_score < score) {
      max_score = score;
      max_score_offset = i;
      max_end_v = end_v;
      max_end_h = end_h;
    }
  }
  // Keep the max-scoring part of the cigar
  const bool cigar_trimmed = (max_score_offset != end_offset-1);
  if (max_score == 0) {
    cigar_clear(cigar);
  } else {
    cigar->operations[max_score_offset+1] = '\0';
    cigar->end_offset = max_score_offset + 1;
    cigar->score = max_score;
    cigar->end_v = max_end_v;
    cigar->end_h = max_end_h;
  }
  // Return
  return cigar_trimmed;
}
int cigar_maxtrim_gap_affine2p_score_op(
    const char operation,
    const int length,
    const affine2p_penalties_t* const penalties,
    int* const end_v,
    int* const end_h) {
  switch (operation) {
    case 'M': {
      *end_v += length; *end_h += length;
      const int match_score = (penalties->match!=0) ? penalties->match : -1;
      return match_score*length;
    }
    case 'X':
      *end_v += length; *end_h += length;
      return penalties->mismatch*length;
    case 'D': {
      *end_v += length;
      const int score1 = penalties->gap_opening1 + penalties->gap_extension1*length;
      const int score2 = penalties->gap_opening2 + penalties->gap_extension2*length;
      return MIN(score1,score2);
    }
    case 'I': {
      *end_h += length;
      const int score1 = penalties->gap_opening1 + penalties->gap_extension1*length;
      const int score2 = penalties->gap_opening2 + penalties->gap_extension2*length;
      return MIN(score1,score2);
    }
    default:
      fprintf(stderr,"[CIGAR] Computing CIGAR score: Unknown operation\n");
      exit(1);
  }
}
bool cigar_maxtrim_gap_affine2p(
    cigar_t* const cigar,
    const affine2p_penalties_t* const penalties) {
  // Parameters
  const char* const operations = cigar->operations;
  const int begin_offset = cigar->begin_offset;
  const int end_offset = cigar->end_offset;
  if (begin_offset >= end_offset) return false;
  // Max-score
  int max_score = 0, max_score_offset = begin_offset, max_end_v = 0, max_end_h = 0;
  // Traverse all cigar
  char last_op = '\0';
  int score = 0, end_v = 0, end_h = 0, op_length = 0;
  int i;
  for (i=begin_offset;i<end_offset;++i) {
    // Account for operation
    const char operation = operations[i];
    if (operation != last_op && last_op != '\0') {
      score -= cigar_maxtrim_gap_affine2p_score_op(last_op,op_length,penalties,&end_v,&end_h);
      op_length = 0;
      // Compare max
      if (max_score < score) {
        max_score = score;
        max_score_offset = i - 1;
        max_end_v = end_v;
        max_end_h = end_h;
      }
    }
    last_op = operation;
    ++op_length;
  }
  // Account for last operation
  score -= cigar_maxtrim_gap_affine2p_score_op(last_op,op_length,penalties,&end_v,&end_h);
  if (max_score < score) {
    max_score = score;
    max_score_offset = end_offset - 1;
    max_end_v = end_v;
    max_end_h = end_h;
  }
  // Keep the max-scoring part of the cigar
  const bool cigar_trimmed = (max_score_offset != end_offset-1);
  if (max_score == 0) {
    cigar_clear(cigar);
  } else {
    cigar->operations[max_score_offset+1] = '\0';
    cigar->end_offset = max_score_offset + 1;
    cigar->score = max_score;
    cigar->end_v = max_end_v;
    cigar->end_h = max_end_h;
  }
  // Return
  return cigar_trimmed;
}
/*
 * Local Alignment Extraction
 */
typedef struct {
  int score;
  int begin_offset;
  int end_offset;
  int begin_v;
  int begin_h;
  int end_v;
  int end_h;
} maxlocal_alg_t;
bool cigar_maxlocal_gap_affine2p(
    cigar_t* const cigar,
    const affine2p_penalties_t* const penalties,
    const int pattern_length,
    const int text_length) {
  // Parameters
  const char* const operations = cigar->operations;
  const int begin_offset = cigar->begin_offset;
  const int end_offset = cigar->end_offset;
  if (begin_offset >= end_offset) return false;
  // Current local alignment (la)
  maxlocal_alg_t la = {
      .begin_offset = begin_offset,
      .end_offset = begin_offset,
      .score = 0,
      .begin_v = 0,
      .begin_h = 0,
      .end_v = 0,
      .end_h = 0
  };
  // Max local alignment (maxla)
  maxlocal_alg_t maxla = la;
  // Traverse all cigar
  char last_op = '\0';
  int op_length = 0, i;
  for (i=begin_offset;i<end_offset;++i) {
    // Account for operation
    const char operation = operations[i];
    if (operation != last_op && last_op != '\0') {
      la.score -= cigar_maxtrim_gap_affine2p_score_op(last_op,op_length,penalties,&la.end_v,&la.end_h);
      op_length = 0;
      // Check maximum local alignment
      if (la.score > maxla.score) {
        la.end_offset = i;
        maxla = la;
      } else if (la.score < 0) { // Check negative score (reset)
        la.score = 0;
        la.begin_offset = i;
        la.end_offset = i;
        la.begin_v = la.end_v;
        la.begin_h = la.end_h;
      }
    }
    last_op = operation;
    ++op_length;
  }
  // Account for the last operation
  la.score -= cigar_maxtrim_gap_affine2p_score_op(last_op,op_length,penalties,&la.end_v,&la.end_h);
  if (la.score > maxla.score) { // Check maximum local alignment
    maxla = la;
    maxla.end_offset = end_offset;
  }
  // Keep the max-scoring part of the cigar
  const bool cigar_trimmed = (maxla.end_offset != end_offset-1);
  if (maxla.score == 0) {
    cigar_clear(cigar);
  } else {
    const int begin_ins = la.begin_h;
    const int begin_del = la.begin_v;
    // Relocate local alignment in buffer
    const int local_cigar_length = maxla.end_offset - maxla.begin_offset;
    memcpy(cigar->operations+begin_ins+begin_del,
           cigar->operations+maxla.begin_offset,
           local_cigar_length);
    // Add initial indels
    char* buffer = cigar->operations;
    memset(buffer,'I',begin_ins); buffer += begin_ins;
    memset(buffer,'D',begin_del); buffer += begin_del + local_cigar_length;
    // Add final indels
    const int end_ins = text_length - la.end_h + 1;
    const int end_del = pattern_length - la.end_v + 1;
    memset(buffer,'I',end_ins); buffer += end_ins;
    memset(buffer,'D',end_del); buffer += end_del;
    *buffer = '\0';
    // Set offsets
    cigar->begin_offset = 0;
    cigar->end_offset = buffer - cigar->operations;
    cigar->end_v = pattern_length;
    cigar->end_h = text_length;
    cigar->score = cigar_score_gap_affine2p(cigar,penalties);
  }
  // Return
  return cigar_trimmed;
}
