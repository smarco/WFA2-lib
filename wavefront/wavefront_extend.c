/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 * Copyright (c) 2022 by Quim Aguado-Puig <quim.aguado.p@gmail.com>
 * Copyright (c) 2022 by Christos Matzoros <christos.matzoros@bsc.es>
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
 * DESCRIPTION: WaveFront-Alignment module for the "extension" of exact matches
 */

#include "utils/string_padded.h"
#include "wavefront_extend.h"
#include "wavefront_align.h"
#include "wavefront_compute.h"
#include "wavefront_heuristic.h"

#if __AVX2__
#include <immintrin.h>
#endif

#ifdef WFA_PARALLEL
#include <omp.h>
#endif

/*
 * Termination (detect end of alignment)
 */
bool wavefront_extend_end2end_check_termination(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int score_mod) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const affine2p_matrix_type component_end = wf_aligner->component_end;
  const int alignment_k = DPMATRIX_DIAGONAL(text_length,pattern_length);
  const wf_offset_t alignment_offset = DPMATRIX_OFFSET(text_length,pattern_length);
  // Select end component
  switch (component_end) {
    case affine2p_matrix_M: {
      // Check diagonal/offset
      if (mwavefront->lo > alignment_k || alignment_k > mwavefront->hi) return false; // Not done
      const wf_offset_t moffset = mwavefront->offsets[alignment_k];
      if (moffset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    case affine2p_matrix_I1: {
      // Fetch I1-wavefront & check diagonal/offset
      wavefront_t* const i1wavefront = wf_aligner->wf_components.i1wavefronts[score_mod];
      if (i1wavefront == NULL || i1wavefront->lo > alignment_k || alignment_k > i1wavefront->hi) return false; // Not done
      const wf_offset_t i1offset = i1wavefront->offsets[alignment_k];
      if (i1offset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    case affine2p_matrix_I2: {
      // Fetch I2-wavefront & check diagonal/offset
      wavefront_t* const i2wavefront = wf_aligner->wf_components.i2wavefronts[score_mod];
      if (i2wavefront == NULL || i2wavefront->lo > alignment_k || alignment_k > i2wavefront->hi) return false; // Not done
      const wf_offset_t i2offset = i2wavefront->offsets[alignment_k];
      if (i2offset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    case affine2p_matrix_D1: {
      // Fetch D1-wavefront & check diagonal/offset
      wavefront_t* const d1wavefront = wf_aligner->wf_components.d1wavefronts[score_mod];
      if (d1wavefront == NULL || d1wavefront->lo > alignment_k || alignment_k > d1wavefront->hi) return false; // Not done
      const wf_offset_t d1offset = d1wavefront->offsets[alignment_k];
      if (d1offset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    case affine2p_matrix_D2: {
      // Fetch D2-wavefront & check diagonal/offset
      wavefront_t* const d2wavefront = wf_aligner->wf_components.d2wavefronts[score_mod];
      if (d2wavefront == NULL || d2wavefront->lo > alignment_k || alignment_k > d2wavefront->hi) return false; // Not done
      const wf_offset_t d2offset = d2wavefront->offsets[alignment_k];
      if (d2offset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    default:
      break;
  }
  return false;
}
bool wavefront_extend_endsfree_check_termination(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int k,
    const wf_offset_t offset) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // Check ends-free reaching boundaries
  const int h_pos = WAVEFRONT_H(k,offset);
  const int v_pos = WAVEFRONT_V(k,offset);
  if (h_pos >= text_length) { // Text is aligned
    // Is Pattern end-free?
    const int pattern_left = pattern_length - v_pos;
    const int pattern_end_free = wf_aligner->alignment_form.pattern_end_free;
    if (pattern_left <= pattern_end_free) {
      #ifdef WFA_PARALLEL
      #pragma omp critical
      #endif
      {
        wf_aligner->alignment_end_pos.score = score;
        wf_aligner->alignment_end_pos.k = k;
        wf_aligner->alignment_end_pos.offset = offset;
      }
      return true; // Quit (we are done)
    }
  }
  if (v_pos >= pattern_length) { // Pattern is aligned
    // Is text end-free?
    const int text_left = text_length - h_pos;
    const int text_end_free = wf_aligner->alignment_form.text_end_free;
    if (text_left <= text_end_free) {
      #ifdef WFA_PARALLEL
      #pragma omp critical
      #endif
      {
        wf_aligner->alignment_end_pos.score = score;
        wf_aligner->alignment_end_pos.k = k;
        wf_aligner->alignment_end_pos.offset = offset;
      }
      return true; // Quit (we are done)
    }
  }
  // Not done
  return false;
}
/*
 * Extend kernel
 */
FORCE_INLINE wf_offset_t wavefront_extend_matches_packed_kernel(
    wavefront_aligner_t* const wf_aligner,
    const int k,
    wf_offset_t offset) {
  // Fetch pattern/text blocks
  uint64_t* pattern_blocks = (uint64_t*)(wf_aligner->pattern+WAVEFRONT_V(k,offset));
  uint64_t* text_blocks = (uint64_t*)(wf_aligner->text+WAVEFRONT_H(k,offset));
  // Compare 64-bits blocks
  uint64_t cmp = *pattern_blocks ^ *text_blocks;
  while (__builtin_expect(cmp==0,0)) {
    // Increment offset (full block)
    offset += 8;
    // Next blocks
    ++pattern_blocks;
    ++text_blocks;
    // Compare
    cmp = *pattern_blocks ^ *text_blocks;
  }
  // Count equal characters
  const int equal_right_bits = __builtin_ctzl(cmp);
  const int equal_chars = DIV_FLOOR(equal_right_bits,8);
  offset += equal_chars;
  // Return extended offset
  return offset;
}

#if __AVX2__
FORCE_INLINE  __m256i avx2_lzcnt_epi32 (__m256i v) {
#if __AVX512CD__ && __AVX512VL__
  return _mm256_lzcnt_epi32(v);
#else
  // Emulate clz for AVX2: https://stackoverflow.com/a/58827596
  v = _mm256_andnot_si256(_mm256_srli_epi32(v, 8), v); // keep 8 MSB
  v = _mm256_castps_si256(_mm256_cvtepi32_ps(v)); // convert an integer to float
  v = _mm256_srli_epi32(v, 23); // shift down the exponent
  v = _mm256_subs_epu16(_mm256_set1_epi32(158), v); // undo bias
  v = _mm256_min_epi16(v, _mm256_set1_epi32(32)); // clamp at 32
  return v;
#endif
}
#endif

/*
 * Wavefront offset extension comparing characters
 *   Remember:
 *   - No offset is out of boundaries !(h>tlen,v>plen)
 *   - if (h==tlen,v==plen) extension won't increment (sentinels)
 */
FORCE_NO_INLINE void wavefront_extend_matches_packed_end2end(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int lo,
    const int hi) {
#if __AVX2__
  // Parameters
  wf_offset_t* const offsets = mwavefront->offsets;
  int k_min = lo;
  int k_max = hi;
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const char* pattern = wf_aligner->pattern;
  const char* text = wf_aligner->text;

  // Extend diagonally each wavefront point
  int k;

  const __m256i vector_null = _mm256_set1_epi32(WAVEFRONT_OFFSET_NULL);
  const __m256i fours = _mm256_set1_epi32(4);
  const __m256i eights = _mm256_set1_epi32(8);
  const __m256i vecShuffle = _mm256_set_epi8(28,29,30,31,24,25,26,27,
                                           20,21,22,23,16,17,18,19,
                                           12,13,14,15, 8, 9,10,11,
                                           4 , 5, 6, 7, 0, 1, 2 ,3);

  const int elems_per_register = 8;

  int num_of_diagonals = k_max-k_min+1;
  int loop_peeling_iters = num_of_diagonals%elems_per_register;

  for(int i=k_min;i<k_min+loop_peeling_iters;i++){
    const uint32_t v = offsets[i] - i; // Make unsigned to avoid checking negative
    if (v >= pattern_length) continue;
    const uint32_t h = offsets[i]; // Make unsigned to avoid checking negative
    if (h >= text_length) continue;

    // Fetch pattern/text blocks
    uint64_t* pattern_blocks = (uint64_t*)(pattern+v);
    uint64_t* text_blocks = (uint64_t*)(text+h);
    uint64_t pattern_block = *pattern_blocks;
    uint64_t text_block = *text_blocks;
    // Compare 64-bits blocks
    uint64_t cmp = pattern_block ^ text_block;

    while (__builtin_expect(!cmp,0)) {
      // Increment offset (full block)
      offsets[i] += 8;
      // Next blocks
      ++pattern_blocks;
      ++text_blocks;
      // Fetch
      pattern_block = *pattern_blocks;
      text_block = *text_blocks;

      // Compare
      cmp = pattern_block ^ text_block;
    }

    // Count equal characters
    const int equal_right_bits = __builtin_ctzl(cmp);
    const int equal_chs = DIV_FLOOR(equal_right_bits,8);

    offsets[i] += equal_chs;
  }

  if (num_of_diagonals <= 8) return;

  k_min+=loop_peeling_iters;

  __m256i ks = _mm256_set_epi32 (k_min+7,k_min+6,k_min+5,k_min+4,k_min+3,
                                 k_min+2,k_min+1,k_min);

  // Main SIMD extension loop
  for (k=k_min;k<=k_max;k+=elems_per_register) {
    __m256i offsets_vector = _mm256_lddqu_si256 ((__m256i*)&offsets[k]);
    __m256i h_vector = offsets_vector;
    __m256i v_vector = _mm256_sub_epi32(offsets_vector, ks);
    ks =_mm256_add_epi32 (ks, eights);

    // NULL offsets will read at index 0 (avoid segfaults)
    __m256i null_mask = _mm256_cmpeq_epi32(offsets_vector, vector_null);
    v_vector = _mm256_andnot_si256(null_mask, v_vector);
    h_vector = _mm256_andnot_si256(null_mask, h_vector);

    __m256i pattern_vector = _mm256_i32gather_epi32((int const*)&pattern[0],v_vector, 1);
    __m256i text_vector = _mm256_i32gather_epi32((int const*)&text[0],h_vector, 1);

    // Change endianess to make the xor + clz character comparison
    pattern_vector = _mm256_shuffle_epi8(pattern_vector, vecShuffle);
    text_vector = _mm256_shuffle_epi8(text_vector, vecShuffle);
    __m256i xor_result_vector = _mm256_xor_si256(pattern_vector,text_vector);
    __m256i clz_vector = avx2_lzcnt_epi32 (xor_result_vector);

    // Divide clz by 8 to get the number of equal characters
    // Assume there are sentinels on sequences so we won't count characters
    // outside the sequences
    __m256i equal_chars =  _mm256_srli_epi32 (clz_vector,3);
    offsets_vector =  _mm256_add_epi32 (offsets_vector, equal_chars);

    v_vector = _mm256_add_epi32 (v_vector, fours);
    h_vector = _mm256_add_epi32 (h_vector, fours);

    // Lanes to continue == 0xffffffff, other lanes = 0
    __m256i vector_mask = _mm256_cmpeq_epi32(equal_chars, fours);

    _mm256_storeu_si256((__m256i*)&offsets[k], offsets_vector);

    int mask = _mm256_movemask_epi8(vector_mask);

    if(mask == 0) continue;

    // ctz(0) is undefined
    while (mask != 0) {
      int tz = __builtin_ctz(mask);
      int curr_k = k + (tz/4);

      const uint32_t v = offsets[curr_k] - (curr_k); // Make unsigned to avoid checking negative
      const uint32_t h = offsets[curr_k]; // Make unsigned to avoid checking negative

      if ((v < pattern_length) && (h < text_length)) {
        // Fetch pattern/text blocks
        uint64_t* pattern_blocks = (uint64_t*)(pattern+v);
        uint64_t* text_blocks = (uint64_t*)(text+h);
        uint64_t pattern_block = *pattern_blocks;
        uint64_t text_block = *text_blocks;
        // Compare 64-bits blocks
        uint64_t cmp = pattern_block ^ text_block;

        while (__builtin_expect(!cmp,0)) {
          // Increment offset (full block)
          offsets[curr_k] += 8;
          // Next blocks
          ++pattern_blocks;
          ++text_blocks;
          // Fetch
          pattern_block = *pattern_blocks;
          text_block = *text_blocks;

          // Compare
          cmp = pattern_block ^ text_block;
        }

        // Count equal characters
        const int equal_right_bits = __builtin_ctzl(cmp);
        const int equal_chs = DIV_FLOOR(equal_right_bits,8);

        offsets[curr_k] += equal_chs;
      }
      mask &= (0xfffffff0 << tz);
    }
  }
#else // No AVX2 available
  wf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=lo;k<=hi;++k) {
    // Fetch offset
    const wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Extend offset
    offsets[k] = wavefront_extend_matches_packed_kernel(wf_aligner,k,offset);
  }
#endif
}

FORCE_NO_INLINE wf_offset_t wavefront_extend_matches_packed_max(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int lo,
    const int hi) {
  wf_offset_t* const offsets = mwavefront->offsets;
  wf_offset_t max_antidiag = 0;
  int k;
  for (k=lo;k<=hi;++k) {
    // Fetch offset
    const wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Extend offset
    offsets[k] = wavefront_extend_matches_packed_kernel(wf_aligner,k,offset);
    // Compute max
    const wf_offset_t antidiag = WAVEFRONT_ANTIDIAGONAL(k,offsets[k]);
    if (max_antidiag < antidiag) max_antidiag = antidiag;
  }
  return max_antidiag;
}
FORCE_NO_INLINE bool wavefront_extend_matches_packed_endsfree(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi) {
  wf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=lo;k<=hi;++k) {
    // Fetch offset
    wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Extend offset
    offset = wavefront_extend_matches_packed_kernel(wf_aligner,k,offset);
    offsets[k] = offset;
    // Check ends-free reaching boundaries
    if (wavefront_extend_endsfree_check_termination(wf_aligner,mwavefront,score,k,offset)) {
      return true; // Quit (we are done)
    }
  }
  // Alignment not finished
  return false;
}
bool wavefront_extend_matches_custom(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi,
    const bool endsfree) {
  // Parameters (custom matching function)
  alignment_match_funct_t match_funct = wf_aligner->match_funct;
  void* const func_arguments = wf_aligner->match_funct_arguments;
  // Extend diagonally each wavefront point
  wf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=lo;k<=hi;++k) {
    // Check offset
    wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Count equal characters
    int v = WAVEFRONT_V(k,offset);
    int h = WAVEFRONT_H(k,offset);
    while (match_funct(v,h,func_arguments)) {
      h++; v++; offset++;
    }
    // Update offset
    offsets[k] = offset;
    // Check ends-free reaching boundaries
    if (endsfree && wavefront_extend_endsfree_check_termination(wf_aligner,mwavefront,score,k,offset)) {
      return true; // Quit (we are done)
    }
  }
  // Alignment not finished
  return false;
}
/*
 * Wavefront exact "extension"
 */
int wavefront_extend_end2end_max(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    int* const max_antidiagonal) {
  // Compute score
  const bool memory_modular = wf_aligner->wf_components.memory_modular;
  const int max_score_scope = wf_aligner->wf_components.max_score_scope;
  const int score_mod = (memory_modular) ? score % max_score_scope : score;
  *max_antidiagonal = 0; // Init
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score_mod];
  if (mwavefront == NULL) {
    // Check alignment feasibility (for heuristic variants that can lead to no solution)
    if (wf_aligner->align_status.num_null_steps > wf_aligner->wf_components.max_score_scope) {
      wf_aligner->align_status.status = WF_STATUS_UNFEASIBLE;
      wf_aligner->align_status.score = score;
      return 1; // Done
    }
    return 0; // Not done
  }
  // Multithreading dispatcher
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  wf_offset_t max_antidiag = 0;
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  if (num_threads == 1) {
    // Extend wavefront
    max_antidiag = wavefront_extend_matches_packed_max(wf_aligner,mwavefront,lo,hi);
  } else {
#ifdef WFA_PARALLEL
    // Extend wavefront in parallel
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      wavefront_compute_thread_limits(
          omp_get_thread_num(),omp_get_num_threads(),lo,hi,&t_lo,&t_hi);
      wf_offset_t t_max_antidiag = wavefront_extend_matches_packed_max(wf_aligner,mwavefront,t_lo,t_hi);
      #ifdef WFA_PARALLEL
      #pragma omp critical
      #endif
      {
        if (t_max_antidiag > max_antidiag) max_antidiag = t_max_antidiag;
      }
    }
#endif
  }
  // Check end-to-end finished
  const bool end_reached = wavefront_extend_end2end_check_termination(wf_aligner,mwavefront,score,score_mod);
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_END_REACHED;
    wf_aligner->align_status.score = score;
    return 1; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    wavefront_heuristic_cufoff(wf_aligner,score,score_mod);
  }
  *max_antidiagonal = max_antidiag;
  return 0; // Not done
}
int wavefront_extend_end2end(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Compute score
  const bool memory_modular = wf_aligner->wf_components.memory_modular;
  const int max_score_scope = wf_aligner->wf_components.max_score_scope;
  const int score_mod = (memory_modular) ? score % max_score_scope : score;
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score_mod];
  if (mwavefront == NULL) {
    // Check alignment feasibility (for heuristic variants that can lead to no solution)
    if (wf_aligner->align_status.num_null_steps > wf_aligner->wf_components.max_score_scope) {
      wf_aligner->align_status.status = WF_STATUS_UNFEASIBLE;
      wf_aligner->align_status.score = score;
      return 1; // Done
    }
    return 0; // Not done
  }
  // Multithreading dispatcher
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  bool end_reached = false;
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  if (num_threads == 1) {
    // Extend wavefront
    wavefront_extend_matches_packed_end2end(wf_aligner,mwavefront,lo,hi);
  } else {
#ifdef WFA_PARALLEL
    // Extend wavefront in parallel
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      wavefront_compute_thread_limits(
          omp_get_thread_num(),omp_get_num_threads(),lo,hi,&t_lo,&t_hi);
      wavefront_extend_matches_packed_end2end(wf_aligner,mwavefront,t_lo,t_hi);
    }
#endif
  }
  // Check end-to-end finished
  end_reached = wavefront_extend_end2end_check_termination(wf_aligner,mwavefront,score,score_mod);
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_END_REACHED;
    wf_aligner->align_status.score = score;
    return 1; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    wavefront_heuristic_cufoff(wf_aligner,score,score_mod);
  }
  return 0; // Not done
}
int wavefront_extend_endsfree(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Modular wavefront
  const bool memory_modular = wf_aligner->wf_components.memory_modular;
  const int max_score_scope = wf_aligner->wf_components.max_score_scope;
  const int score_mod = (memory_modular) ? score % max_score_scope : score;
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score_mod];
  if (mwavefront == NULL) {
    // Check alignment feasibility (for heuristic variants that can lead to no solution)
    if (wf_aligner->align_status.num_null_steps > wf_aligner->wf_components.max_score_scope) {
      wf_aligner->align_status.status = WF_STATUS_UNFEASIBLE;
      wf_aligner->align_status.score = score;
      return 1; // Done
    }
    return 0; // Not done
  }
  // Multithreading dispatcher
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  bool end_reached = false;
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  if (num_threads == 1) {
    // Extend wavefront
    end_reached = wavefront_extend_matches_packed_endsfree(wf_aligner,mwavefront,score,lo,hi);
  } else {
#ifdef WFA_PARALLEL
    // Extend wavefront in parallel
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      wavefront_compute_thread_limits(
          omp_get_thread_num(),omp_get_num_threads(),lo,hi,&t_lo,&t_hi);
      if (wavefront_extend_matches_packed_endsfree(wf_aligner,mwavefront,score,t_lo,t_hi)) {
        end_reached = true;
      }
    }
#endif
  }
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_END_REACHED;
    wf_aligner->align_status.score = score;
    return 1; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    wavefront_heuristic_cufoff(wf_aligner,score,score_mod);
  }
  return 0; // Not done
}
int wavefront_extend_custom(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Compute score
  const bool memory_modular = wf_aligner->wf_components.memory_modular;
  const int max_score_scope = wf_aligner->wf_components.max_score_scope;
  const int score_mod = (memory_modular) ? score % max_score_scope : score;
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score_mod];
  if (mwavefront == NULL) {
    // Check alignment feasibility (for heuristic variants that can lead to no solution)
    if (wf_aligner->align_status.num_null_steps > wf_aligner->wf_components.max_score_scope) {
      wf_aligner->align_status.status = WF_STATUS_UNFEASIBLE;
      wf_aligner->align_status.score = score;
      return 1; // Done
    }
    return 0; // Not done
  }
  // Multithreading dispatcher
  const bool endsfree = (wf_aligner->alignment_form.span == alignment_endsfree);
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  bool end_reached = false;
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  if (num_threads == 1) {
    // Extend wavefront
    end_reached = wavefront_extend_matches_custom(wf_aligner,mwavefront,score,lo,hi,endsfree);
  } else {
#ifdef WFA_PARALLEL
    // Extend wavefront in parallel
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      wavefront_compute_thread_limits(
          omp_get_thread_num(),omp_get_num_threads(),lo,hi,&t_lo,&t_hi);
      if (wavefront_extend_matches_custom(wf_aligner,mwavefront,score,t_lo,t_hi,endsfree)) {
        end_reached = true;
      }
    }
#endif
  }
  // Check end-to-end finished
  if (!endsfree) {
    end_reached = wavefront_extend_end2end_check_termination(wf_aligner,mwavefront,score,score_mod);
  }
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_END_REACHED;
    wf_aligner->align_status.score = score;
    return 1; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    wavefront_heuristic_cufoff(wf_aligner,score,score_mod);
  }
  return 0; // Not done
}


