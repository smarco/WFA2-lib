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
 * DESCRIPTION: WFA Sample-Code
 */

#include <iostream>
#include <string>
#include "bindings/cpp/WFAligner.hpp"

using namespace std;

int main(int argc,char* argv[]) {
  // Patter & Text
  string pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
  string text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";

  // Create a WFAligner
  WFAlignerGapAffine aligner(4,6,2);
  // Align
  aligner.align(pattern,text);
  cout << "WFA-Alignemnt returns score " << aligner.get_alignment_score() << endl;

  // Print CIGAR
  string cigar = aligner.get_alignment_cigar();
  cout << "PATTERN: " << pattern  << endl;
  cout << "TEXT: " << text  << endl;
  cout << "CIGAR: " << cigar  << endl;
}
