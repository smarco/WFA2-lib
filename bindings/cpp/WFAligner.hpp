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
 * DESCRIPTION: C++ bindings for the WaveFront Alignment modules
 */

#ifndef BINDINGS_CPP_WFALIGNER_HPP_
#define BINDINGS_CPP_WFALIGNER_HPP_

#include <string>

extern "C" {
  #include "../../wavefront/wavefront_aligner.h"
}

/*
 * Namespace
 */
namespace wfa {

/*
 * General Wavefront Aligner
 */
class WFAligner {
public:
  // Configuration
  enum MemoryModel {
    WavefrontMemoryFull,
    WavefrontMemoryHigh,
    WavefrontMemoryMed,
    WavefrontMemoryLow
  };
  void setReductionNone();
  void setReductionAdaptive(
      const int minWavefrontLength,
      const int maxDistanceThreshold);
  void setMatchFunct(
      int (*matchFunct)(int,int,void*),
      void* matchFunctArguments);
  void setMaxAlignmentScore(
      const int maxAlignmentScore);
  void setMaxMemory(
      const uint64_t maxMemoryCompact,
      const uint64_t maxMemoryResident,
      const uint64_t maxMemoryAbort);
  // Accessors
  int getAlignmentScore();
  void getAlignmentCigar(
      char** const cigarOperations,
      int* cigarLength);
  std::string getAlignmentCigar();
  // Align End-to-end
  int alignEnd2End(
      const int patternLength,
      const int textLength);
  int alignEnd2End(
      const char* const pattern,
      const int patternLength,
      const char* const text,
      const int textLength);
  int alignEnd2End(
      std::string& pattern,
      std::string& text);
  // Align Ends-free
  int alignEndsFree(
      const int patternLength,
      const int patternBeginFree,
      const int patternEndFree,
      const int textLength,
      const int textBeginFree,
      const int textEndFree);
  int alignEndsFree(
      const char* const pattern,
      const int patternLength,
      const int patternBeginFree,
      const int patternEndFree,
      const char* const text,
      const int textLength,
      const int textBeginFree,
      const int textEndFree);
  int alignEndsFree(
      std::string& pattern,
      const int patternBeginFree,
      const int patternEndFree,
      std::string& text,
      const int textBeginFree,
      const int textEndFree);
  // Display & errors
  char* strError(
      const int wfErrorCode);
protected:
  wavefront_aligner_attr_t attributes;
  wavefront_aligner_t* wfAligner;
  // Setup
  WFAligner(
      const bool onlyScore = false,
      const MemoryModel memoryModel = WavefrontMemoryFull);
  ~WFAligner();
};
/*
 * Gap-Affine Aligner
 */
class WFAlignerGapAffine : public WFAligner {
public:
  WFAlignerGapAffine(
      const int mismatch,
      const int gapOpening,
      const int gapExtension,
      const bool onlyScore = false,
      const MemoryModel memoryModel = WFAligner::WavefrontMemoryFull);
};
/*
 * Gap-Affine Dual-Cost (2 pieces) Aligner
 */
class WFAlignerGapAffine2Pieces : public WFAligner {
public:
  WFAlignerGapAffine2Pieces(
      const int mismatch,
      const int gapOpening1,
      const int gapExtension1,
      const int gapOpening2,
      const int gapExtension2,
      const bool onlyScore = false,
      const MemoryModel memoryModel = WFAligner::WavefrontMemoryFull);
};

} /* namespace wfa */

#endif /* BINDINGS_CPP_WFALIGNER_HPP_ */
