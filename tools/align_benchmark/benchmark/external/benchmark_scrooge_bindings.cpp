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

#include <string.h>
#include <scrooge/genasm_cpu.hpp>

using namespace std;

/*
 * Adapt CIGAR
 */
void benchmark_scrooge_adapt_cigar(
    const int pattern_length,
    const int text_length,
    char* const scooge_cigar,
    char* const edit_operations,
    int* const num_edit_operations) {
  // Init
  *num_edit_operations = 0;
  // Decode all CIGAR operations
  const int cigar_length = strlen(scooge_cigar);
  int chars_read = 0, plen = 0, tlen = 0;
  int length, n, i;
  char operation;
  while (chars_read < cigar_length) {
    // Read operation
    sscanf(scooge_cigar+chars_read,"%d%c%n",&length,&operation,&n);
    chars_read += n;
    // Adapt operation encoding
    if (operation=='=' || operation=='M') {
      operation = 'M'; plen+=length; tlen+=length;
    } else if (operation=='X') {
      operation = 'X'; plen+=length; tlen+=length;
    } else if (operation=='D') {
      operation = 'I'; tlen+=length;
    } else if (operation=='I') {
      operation = 'D'; plen+=length;
    } else {
      fprintf(stderr,"[Scrooge] Computing CIGAR score: Unknown operation '%c'\n",operation);
      exit(1);
    }
    // Dump operation
    for (i=0;i<length;++i) {
      edit_operations[(*num_edit_operations)++] = operation;
    }
  }
  // Add final indel
  if (plen < pattern_length) {
    for (i=0;i<pattern_length-plen;++i) {
      edit_operations[(*num_edit_operations)++] = 'D';
    }
  }
  if (tlen < text_length) {
    for (i=0;i<text_length-tlen;++i) {
      edit_operations[(*num_edit_operations)++] = 'I';
    }
  }
}

/*
 * Benchmark Scrooge
 */
extern "C" void benchmark_scrooge_bridge(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    char* const edit_operations,
    int* const num_edit_operations,
    uint64_t* const time_ns) {
  // Parameters
  std::vector<std::string> texts;
  std::vector<std::string> queries;

  // Prepare data
  texts.push_back(string(text,text_length));
  queries.push_back(string(pattern,pattern_length));

  // Align
  long long core_algorithm_ns;
  vector<Alignment_t> alignments = genasm_cpu::align_all(texts,queries,1,&core_algorithm_ns);
  //fprintf(stderr,"[Scrooge::debug]>>%s\n",alignments[0].cigar.c_str());
  *time_ns = core_algorithm_ns;

  // Process alignment
  benchmark_scrooge_adapt_cigar(
      pattern_length,text_length,(char*)alignments[0].cigar.c_str(),
      edit_operations,num_edit_operations);
}


