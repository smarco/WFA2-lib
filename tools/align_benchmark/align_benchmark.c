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
 * DESCRIPTION: Wavefront Alignment Algorithms benchmarking tool
 */

#define EXTERNAL_BENCHMARKS

#include "utils/commons.h"
#include "system/profiler_timer.h"

#include "alignment/score_matrix.h"
#include "edit/edit_dp.h"
#include "gap_linear/nw.h"
#include "gap_affine/swg.h"
#include "gap_affine2p/affine2p_matrix.h"
#include "gap_affine2p/affine2p_dp.h"
#include "wavefront/wavefront_align.h"

#include "benchmark/benchmark_indel.h"
#include "benchmark/benchmark_edit.h"
#include "benchmark/benchmark_gap_linear.h"
#include "benchmark/benchmark_gap_affine.h"
#include "benchmark/benchmark_gap_affine2p.h"

#ifdef EXTERNAL_BENCHMARKS
#include "benchmark/external/benchmark_bitpal.h"
#include "benchmark/external/benchmark_blockaligner.h"
#include "benchmark/external/benchmark_daligner.h"
#include "benchmark/external/benchmark_diffutils.h"
#include "benchmark/external/benchmark_edlib.h"
#include "benchmark/external/benchmark_gaba.h"
#include "benchmark/external/benchmark_ksw2.h"
#include "benchmark/external/benchmark_lv89.h"
#include "benchmark/external/benchmark_parasail.h"
#include "benchmark/external/benchmark_seqan.h"
#include "benchmark/external/benchmark_wfalm.h"
#endif

/*
 * Algorithms
 */
typedef enum {
  // Test
  alignment_test,
  // Indel
  alignment_indel_wavefront,
  // Edit
  alignment_edit_bpm,
  alignment_edit_dp,
  alignment_edit_dp_banded,
  alignment_edit_wavefront,
  // Gap-linear
  alignment_gap_linear_nw,
  alignment_gap_linear_wavefront,
  // Gap-affine
  alignment_gap_affine_swg,
  alignment_gap_affine_swg_endsfree,
  alignment_gap_affine_swg_banded,
  alignment_gap_affine_wavefront,
  // Gap-affine dual-cost
  alignment_gap_affine2p_dp,
  alignment_gap_affine2p_wavefront,
#ifdef EXTERNAL_BENCHMARKS
  // External algorithms
  alignment_bitpal_edit,
  alignment_bitpal_scored,
  alignment_blockaligner,
  alignment_daligner,
  alignment_diffutils,
  alignment_edlib,
  alignment_gaba_aband,
  alignment_ksw2_extz2_sse,
  alignment_ksw2_extd2_sse,
  alignment_lv89,
  alignment_parasail_nw_stripped,
  alignment_parasail_nw_scan,
  alignment_parasail_nw_diag,
  alignment_parasail_nw_banded,
  alignment_seqan_edit,
  alignment_seqan_edit_bpm,
  alignment_seqan_lineal,
  alignment_seqan_affine,
  alignment_wfalm,
  alignment_wfalm_lowmem
#endif
} alignment_algorithm_type;
bool align_benchmark_is_wavefront(
    const alignment_algorithm_type algorithm) {
  return algorithm == alignment_indel_wavefront ||
         algorithm == alignment_edit_wavefront ||
         algorithm == alignment_gap_linear_wavefront ||
         algorithm == alignment_gap_affine_wavefront ||
         algorithm == alignment_gap_affine2p_wavefront;
}
/*
 * Generic parameters
 */
typedef struct {
  // Algorithm
  alignment_algorithm_type algorithm;
  // Input
  char *input_filename;
  char *output_filename;
  bool output_full;
  // Penalties
  linear_penalties_t linear_penalties;
  affine_penalties_t affine_penalties;
  affine2p_penalties_t affine2p_penalties;
  // Alignment form
  bool endsfree;
  double pattern_begin_free;
  double text_begin_free;
  double pattern_end_free;
  double text_end_free;
  // Wavefront parameters
  bool wfa_score_only;
  wavefront_reduction_type wfa_reduction_type;
  int wfa_min_wf_length;
  int wfa_max_dist_th;
  wavefront_memory_t wfa_memory_mode;
  alignment_match_funct_t wfa_match_funct;
  void* wfa_match_funct_arguments;
  uint64_t wfa_max_memory;
  // Other algorithms parameters
  int bandwidth;
#ifdef EXTERNAL_BENCHMARKS
  int ba_block_size;
  bool ksw2_approx_max__drop;
  int ksw2_bandwidth;
  int ksw2_zdrop;
#endif
  // Misc
  bool check_display;
  bool check_correct;
  bool check_score;
  bool check_alignments;
  int check_metric;
  int check_bandwidth;
  int plot;
  // Profile
  profiler_timer_t timer_global;
  // System
  int progress;
  int verbose;
} benchmark_args;
benchmark_args parameters = {
  // Algorithm
  .algorithm = alignment_test,
  // Input
  .input_filename = NULL,
  .output_filename = NULL,
  .output_full = false,
  // Penalties
  .linear_penalties = {
      .match = 0,
      .mismatch = 4,
      .indel = 2,
  },
  .affine_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening = 6,
      .gap_extension = 2,
  },
  .affine2p_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening1 = 6,
      .gap_extension1 = 2,
      .gap_opening2 = 24,
      .gap_extension2 = 1,
  },
  // Alignment form
  .endsfree = false,
  .pattern_begin_free = 0.0,
  .text_begin_free = 0.0,
  .pattern_end_free = 0.0,
  .text_end_free = 0.0,
  // Wavefront parameters
  .wfa_score_only = false,
  .wfa_reduction_type = wavefront_reduction_none,
  .wfa_min_wf_length = 10,
  .wfa_max_dist_th = 50,
  .wfa_memory_mode = wavefront_memory_high,
  .wfa_match_funct = NULL,
  .wfa_match_funct_arguments = NULL,
  .wfa_max_memory = UINT64_MAX,
  // Other algorithms parameters
  .bandwidth = -1,
#ifdef EXTERNAL_BENCHMARKS
  .ba_block_size = 256,
  .ksw2_approx_max__drop = false,
  .ksw2_bandwidth = -1,
  .ksw2_zdrop = -1,
#endif
  // Misc
  .check_bandwidth = -1,
  .check_display = false,
  .check_correct = false,
  .check_score = false,
  .check_alignments = false,
  .check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE,
  .plot = 0,
  // System
  .progress = 10000,
  .verbose = 0,
};

/*
 * Benchmark UTest
 */
void align_pairwise_test() {
  // Patters & Texts
//  char * pattern = "GATTACA";
//  char * text = "GATCACTA";
  char* pattern = "AAGGGTAATCTAAGTGTCTGGTCCTTTGTCATTCTGACTTTCTTCATAATGTGATCTCCTCACCTC";
  char* text =    "AAGGGCAATTTAAGTGTCTGGTCCTTTGTCATTCGACTTTCTTCATAATGTGATCTTCCTCACCTC";

  // MMAllocator
  mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  // Penalties
  linear_penalties_t linear_penalties = {
      .match = 0,
      .mismatch = 4,
      .indel = 2,
  };
  affine_penalties_t affine_penalties = {
      .match = 0,
      .mismatch = 4, //9,
      .gap_opening = 6, //13,
      .gap_extension = 2,
  };
  // Ends
  const int pattern_begin_free = 0;
  const int pattern_end_free = 0;
  const int text_begin_free = 0;
  const int text_end_free = 0;
  const bool endsfree =
      pattern_begin_free>0 || pattern_end_free>0 ||
      text_begin_free>0 || text_end_free>0;
  /*
   * Gap-Affine
   */
  // Allocate
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine;
  attributes.linear_penalties = linear_penalties;
  attributes.affine_penalties = affine_penalties;
  // attributes.affine2p_penalties = affine2p_penalties;
  attributes.reduction.reduction_strategy = wavefront_reduction_none;
  attributes.reduction.min_wavefront_length = 256;
  attributes.reduction.max_distance_threshold = 4096;
  attributes.alignment_scope = compute_alignment; // compute_score
  attributes.memory_mode = wavefront_memory_med;
  attributes.alignment_form.span = (endsfree) ? alignment_endsfree : alignment_end2end;
  attributes.alignment_form.pattern_begin_free = pattern_begin_free;
  attributes.alignment_form.pattern_end_free = pattern_end_free;
  attributes.alignment_form.text_begin_free = text_begin_free;
  attributes.alignment_form.text_end_free = text_end_free;
  attributes.plot_params.plot_enabled = false;
  attributes.mm_allocator = mm_allocator;
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  // Align
  wavefront_align(wf_aligner,
      pattern,strlen(pattern),text,strlen(text));
  // CIGAR
  fprintf(stderr,">> WFA2");
  cigar_print_pretty(stderr,
      pattern,strlen(pattern),text,strlen(text),
      &wf_aligner->cigar,mm_allocator);
  fprintf(stderr,"SCORE: %d \n",cigar_score_gap_affine(&wf_aligner->cigar,&affine_penalties));
  // Plot
  if (attributes.plot_params.plot_enabled) {
    FILE* const wf_plot = fopen("test.wfa","w");
    wavefront_plot_print(wf_plot,wf_aligner);
    fclose(wf_plot);
  }
  // Free
  wavefront_aligner_delete(wf_aligner);
  mm_allocator_delete(mm_allocator);
}
/*
 * Simplest Extend-matching function (for testing purposes)
 */
typedef struct {
  char* pattern;
  int pattern_length;
  char* text;
  int text_length;
} match_function_params_t;
match_function_params_t match_function_params;
int match_function(int v,int h,void* arguments) {
  // Extract parameters
  match_function_params_t* match_arguments = (match_function_params_t*)arguments;
  // Check match
  if (v > match_arguments->pattern_length || h > match_arguments->text_length) return 0;
  return (match_arguments->pattern[v] == match_arguments->text[h]);
}
/*
 * Configuration
 */
wavefront_aligner_t* align_benchmark_configure_wf(
    align_input_t* const align_input,
    mm_allocator_t* const mm_allocator) {
  // Set attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.memory_mode = parameters.wfa_memory_mode;
  attributes.mm_allocator = mm_allocator;
  if (parameters.wfa_score_only) {
    attributes.alignment_scope = compute_score;
  }
  // WF-Reduction
  if (parameters.wfa_reduction_type == wavefront_reduction_adaptive) {
    attributes.reduction.reduction_strategy = wavefront_reduction_adaptive;
    attributes.reduction.min_wavefront_length = parameters.wfa_min_wf_length;
    attributes.reduction.max_distance_threshold = parameters.wfa_max_dist_th;
  } else {
    attributes.reduction.reduction_strategy = wavefront_reduction_none;
    attributes.reduction.min_wavefront_length = -1;
    attributes.reduction.max_distance_threshold = -1;
  }
  // Select flavor
  switch (parameters.algorithm) {
    case alignment_indel_wavefront:
      attributes.distance_metric = indel;
      break;
    case alignment_edit_wavefront:
      attributes.distance_metric = edit;
      break;
    case alignment_gap_linear_wavefront:
      attributes.distance_metric = gap_linear;
      attributes.linear_penalties = parameters.linear_penalties;
      break;
    case alignment_gap_affine_wavefront:
      attributes.distance_metric = gap_affine;
      attributes.affine_penalties = parameters.affine_penalties;
      break;
    case alignment_gap_affine2p_wavefront:
      attributes.distance_metric = gap_affine_2p;
      attributes.affine2p_penalties = parameters.affine2p_penalties;
      break;
    default:
      return NULL; // No WF selected
      break;
  }
  // Select alignment form
  attributes.alignment_form.span = (parameters.endsfree) ? alignment_endsfree : alignment_end2end;
  // Misc
  if (parameters.wfa_match_funct_arguments != NULL) {
    attributes.match_funct = parameters.wfa_match_funct;
    attributes.match_funct_arguments = parameters.wfa_match_funct_arguments;
  }
  attributes.plot_params.plot_enabled = (parameters.plot > 0);
  attributes.plot_params.resolution_points = parameters.plot;
  attributes.system.verbose = parameters.verbose;
  attributes.system.max_memory_abort = parameters.wfa_max_memory;
  // Allocate
  return wavefront_aligner_new(&attributes);
}
void align_benchmark_configure_global(
    align_input_t* const align_input) {
  // Clear
  benchmark_align_input_clear(align_input);
  // Alignment form
  align_input->ends_free = parameters.endsfree;
  // Output
  if (parameters.output_filename == NULL) {
    align_input->output_file = NULL;
  } else {
    align_input->output_file = fopen(parameters.output_filename, "w");
  }
  align_input->output_full = parameters.output_full;
  // MM
  align_input->mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  // WFA
  if (align_benchmark_is_wavefront(parameters.algorithm)) {
    align_input->wf_aligner = align_benchmark_configure_wf(align_input,align_input->mm_allocator);
  } else {
    align_input->wf_aligner = NULL;
  }
  // PROFILE/STATS
  timer_reset(&align_input->timer);
  // DEBUG
  align_input->debug_flags = 0;
  align_input->debug_flags |= parameters.check_metric;
  if (parameters.check_display) align_input->debug_flags |= ALIGN_DEBUG_DISPLAY_INFO;
  if (parameters.check_correct) align_input->debug_flags |= ALIGN_DEBUG_CHECK_CORRECT;
  if (parameters.check_score) align_input->debug_flags |= ALIGN_DEBUG_CHECK_SCORE;
  if (parameters.check_alignments) align_input->debug_flags |= ALIGN_DEBUG_CHECK_ALIGNMENT;
  align_input->check_linear_penalties = &parameters.linear_penalties;
  align_input->check_affine_penalties = &parameters.affine_penalties;
  align_input->check_affine2p_penalties = &parameters.affine2p_penalties;
  align_input->check_bandwidth = parameters.check_bandwidth;
  align_input->verbose = parameters.verbose;
}
void align_benchmark_configure_local(
    align_input_t* const align_input) {
  // Ends-free configuration
  if (parameters.endsfree) {
    const int plen = align_input->pattern_length;
    const int tlen = align_input->text_length;
    align_input->pattern_begin_free = nominal_prop_u32(plen,parameters.pattern_begin_free);
    align_input->pattern_end_free = nominal_prop_u32(plen,parameters.pattern_end_free);
    align_input->text_begin_free = nominal_prop_u32(tlen,parameters.text_begin_free);
    align_input->text_end_free = nominal_prop_u32(tlen,parameters.text_end_free);
    if (align_benchmark_is_wavefront(parameters.algorithm)) {
      wavefront_aligner_set_alignment_free_ends(align_input->wf_aligner,
          align_input->pattern_begin_free,align_input->pattern_end_free,
          align_input->text_begin_free,align_input->text_end_free);
    }
  }
  // Custom extend-match function
  if (parameters.wfa_match_funct != NULL) {
    match_function_params.pattern = align_input->pattern;
    match_function_params.pattern_length = align_input->pattern_length;
    match_function_params.text = align_input->text;
    match_function_params.text_length = align_input->text_length;
  }
}
/*
 * I/O
 */
bool align_benchmark_read_input(
    FILE* input_file,
    char** line1,
    char** line2,
    size_t* line1_allocated,
    size_t* line2_allocated,
    const int seqs_processed,
    align_input_t* const align_input) {
  // Parameters
  int line1_length=0, line2_length=0;
  // Read queries
  line1_length = getline(line1,line1_allocated,input_file);
  if (line1_length==-1) return false;
  line2_length = getline(line2,line2_allocated,input_file);
  if (line1_length==-1) return false;
  // Configure input
  align_input->sequence_id = seqs_processed;
  align_input->pattern = *line1 + 1;
  align_input->pattern_length = line1_length - 2;
  align_input->pattern[align_input->pattern_length] = '\0';
  align_input->text = *line2 + 1;
  align_input->text_length = line2_length - 2;
  align_input->text[align_input->text_length] = '\0';
  return true;
}
/*
 * Display
 */
void align_benchmark_print_progress(
    align_input_t* const align_input,
    const int seqs_processed) {
  const uint64_t time_elapsed_alg = timer_get_total_ns(&(align_input->timer));
  const float rate_alg = (float)seqs_processed/(float)TIMER_CONVERT_NS_TO_S(time_elapsed_alg);
  fprintf(stderr,"...processed %d reads (alignment = %2.3f seq/s)\n",seqs_processed,rate_alg);
}
void align_benchmark_print_results(
    align_input_t* const align_input,
    const int seqs_processed,
    const bool print_stats) {
  // Print benchmark results
  fprintf(stderr,"[Benchmark]\n");
  fprintf(stderr,"=> Total.reads            %d\n",seqs_processed);
  fprintf(stderr,"=> Time.Benchmark      ");
  timer_print(stderr,&parameters.timer_global,NULL);
  fprintf(stderr,"  => Time.Alignment    ");
  timer_print(stderr,&align_input->timer,&parameters.timer_global);
  // Print Stats
  if (parameters.check_display || parameters.check_correct || parameters.check_score || parameters.check_alignments) {
    const bool print_wf_stats = (parameters.algorithm == alignment_gap_affine_wavefront);
    benchmark_print_stats(stderr,align_input,print_wf_stats);
  }
}
void align_benchmark_plot_wf(
    align_input_t* const align_input,
    const int seq_id) {
  // Setup filename
  char filename[100];
  sprintf(filename,"%s.seq%03d.wfa",parameters.input_filename,seq_id);
  // Open file
  FILE* const wf_plot = fopen(filename,"w");
  wavefront_plot_print(wf_plot,align_input->wf_aligner);
  fclose(wf_plot);
}
/*
 * Benchmark
 */
void align_benchmark() {
  // Parameters
  FILE *input_file = NULL;
  char *line1 = NULL, *line2 = NULL;
  size_t line1_allocated=0, line2_allocated=0;
  align_input_t align_input;
  // PROFILE
  timer_reset(&(parameters.timer_global));
  timer_start(&(parameters.timer_global));
  // Initialize files and configure align-benchmark
  input_file = fopen(parameters.input_filename, "r");
  if (input_file == NULL) {
    fprintf(stderr,"Input file '%s' couldn't be opened\n",parameters.input_filename);
    exit(1);
  }
  // Global configuration
  align_benchmark_configure_global(&align_input);
  // Read-align loop
  int seqs_processed = 0, progress = 0;
  while (true) {
    // Read input sequence-pair
    const bool input_read = align_benchmark_read_input(
        input_file,&line1,&line2,&line1_allocated,
        &line2_allocated,seqs_processed,&align_input);
    if (!input_read) break;
    // Sequence-dependent configuration
    align_benchmark_configure_local(&align_input);
    // Align queries using DP
    switch (parameters.algorithm) {
      /*
       * Algorithms
       */
      // Indel
      case alignment_indel_wavefront:
        benchmark_indel_wavefront(&align_input);
        break;
      // Edit
      case alignment_edit_bpm:
        benchmark_edit_bpm(&align_input);
        break;
      case alignment_edit_dp:
        benchmark_edit_dp(&align_input);
        break;
      case alignment_edit_dp_banded:
        benchmark_edit_dp_banded(&align_input,parameters.bandwidth);
        break;
      case alignment_edit_wavefront:
        benchmark_edit_wavefront(&align_input);
        break;
      // Gap-linear
      case alignment_gap_linear_nw:
        benchmark_gap_linear_nw(&align_input,&parameters.linear_penalties);
        break;
      case alignment_gap_linear_wavefront:
        benchmark_gap_linear_wavefront(&align_input,&parameters.linear_penalties);
        break;
      // Gap-affine
      case alignment_gap_affine_swg:
        benchmark_gap_affine_swg(&align_input,&parameters.affine_penalties);
        break;
      case alignment_gap_affine_swg_endsfree:
        benchmark_gap_affine_swg_endsfree(
            &align_input,&parameters.affine_penalties);
        break;
      case alignment_gap_affine_swg_banded:
        benchmark_gap_affine_swg_banded(&align_input,
            &parameters.affine_penalties,parameters.bandwidth);
        break;
      case alignment_gap_affine_wavefront:
        benchmark_gap_affine_wavefront(&align_input,&parameters.affine_penalties);
        break;
      // Gap-affine 2p
      case alignment_gap_affine2p_dp:
        benchmark_gap_affine2p_dp(&align_input,&parameters.affine2p_penalties);
        break;
      case alignment_gap_affine2p_wavefront:
        benchmark_gap_affine2p_wavefront(&align_input,&parameters.affine2p_penalties);
        break;
#ifdef EXTERNAL_BENCHMARKS
      /*
       * External Algorithms
       */
      case alignment_bitpal_edit:
        benchmark_bitpal_m0_x1_g1(&align_input);
        break;
      case alignment_bitpal_scored:
        benchmark_bitpal_m1_x4_g2(&align_input);
        break;
      case alignment_blockaligner:
        benchmark_blockaligner_global_affine(
            &align_input,&parameters.affine_penalties,
            parameters.ba_block_size);
        break;
      case alignment_daligner:
        benchmark_daligner(&align_input);
        break;
      case alignment_diffutils:
        benchmark_diffutils(&align_input,true);
        break;
      case alignment_edlib:
        benchmark_edlib(&align_input);
        break;
      case alignment_gaba_aband:
        benchmark_gaba_aband(&align_input,&parameters.affine_penalties);
        break;
      case alignment_ksw2_extz2_sse:
        benchmark_ksw2_extz2_sse(
            &align_input,&parameters.affine_penalties,
            parameters.ksw2_approx_max__drop,
            parameters.ksw2_bandwidth,parameters.ksw2_zdrop);
        break;
      case alignment_ksw2_extd2_sse:
        benchmark_ksw2_extd2_sse(
            &align_input,&parameters.affine2p_penalties,
            parameters.ksw2_approx_max__drop,
            parameters.ksw2_bandwidth,parameters.ksw2_zdrop);
        break;
      case alignment_lv89:
        benchmark_lv89(&align_input);
        break;
      case alignment_parasail_nw_stripped:
        benchmark_parasail_nw_stripped(&align_input,&parameters.affine_penalties);
        break;
      case alignment_parasail_nw_scan:
        benchmark_parasail_nw_scan(&align_input,&parameters.affine_penalties);
        break;
      case alignment_parasail_nw_diag:
        benchmark_parasail_nw_diag(&align_input,&parameters.affine_penalties);
        break;
      case alignment_parasail_nw_banded:
        benchmark_parasail_nw_banded(&align_input,&parameters.affine_penalties,parameters.bandwidth);
        break;
      case alignment_seqan_edit:
        benchmark_seqan_global_edit(&align_input);
        break;
      case alignment_seqan_edit_bpm:
        benchmark_seqan_global_edit_bpm(&align_input);
        break;
      case alignment_seqan_lineal:
        benchmark_seqan_global_lineal(&align_input,&parameters.linear_penalties);
        break;
      case alignment_seqan_affine:
        benchmark_seqan_global_affine(&align_input,&parameters.affine_penalties);
        break;
      case alignment_wfalm:
        benchmark_wfalm_global_affine(&align_input,&parameters.affine_penalties);
        break;
      case alignment_wfalm_lowmem:
        benchmark_wfalm_global_affine_lowmem(&align_input,&parameters.affine_penalties);
        break;
#endif
      default:
        fprintf(stderr,"Algorithm not implemented\n");
        exit(1);
        break;
    }
    // Update progress
    ++seqs_processed;
    if (++progress == parameters.progress) {
      progress = 0;
      align_benchmark_print_progress(&align_input,seqs_processed);
    }
    // DEBUG mm_allocator_print(stderr,align_input.mm_allocator,true);
    // Plot
    if (parameters.plot > 0) align_benchmark_plot_wf(&align_input,seqs_processed);
  }
  timer_stop(&(parameters.timer_global));
  // Print benchmark results
  align_benchmark_print_results(&align_input,seqs_processed,true);
  // Free
  fclose(input_file);
  if (align_input.output_file != NULL) fclose(align_input.output_file);
  if (align_input.wf_aligner) wavefront_aligner_delete(align_input.wf_aligner);
  mm_allocator_delete(align_input.mm_allocator);
  free(line1);
  free(line2);
}
/*
 * Generic Menu
 */
void usage() {
  fprintf(stderr,
      "USE: ./align_benchmark -a <algorithm> -i <input>                        \n"
      "      Options::                                                         \n"
      "        [Algorithm]                                                     \n"
      "          --algorithm|a <algorithm>                                     \n"
      "            [Indel (Longest Common Subsequence)]                        \n"
      "              indel-wfa                                                 \n"
      "            [Edit (Levenshtein)]                                        \n"
      "              edit-bpm                                                  \n"
      "              edit-dp                                                   \n"
      "              edit-dp-banded                                            \n"
      "              edit-wfa                                                  \n"
      "            [Gap-linear (Needleman-Wunsch)]                             \n"
      "              gap-linear-nw                                             \n"
      "              gap-linear-wfa                                            \n"
      "              gap-linear-wfa-adaptive                                   \n"
      "            [Gap-affine (Smith-Waterman-Gotoh)]                         \n"
      "              gap-affine-swg                                            \n"
      "              gap-affine-swg-banded                                     \n"
      "              gap-affine-wfa                                            \n"
      "            [Gap-affine-2pieces (Concave 2-pieces)]                     \n"
      "              gap-affine2p-dp                                           \n"
      "              gap-affine2p-wfa                                          \n"
#ifdef EXTERNAL_BENCHMARKS
      "            [External/BitPal]                                           \n"
      "              bitpal-edit          (Edit)[score-only]                   \n"
      "              bitpal-scored        (Gap-linear)[score-only]             \n"
      "            [External/BlockAligner]                                     \n"
      "              block-aligner        (Gap-affine)                         \n"
      "            [External/Daligner]                                         \n"
      "              daligner             (Edit)                               \n"
      "            [External/Diffutils]                                        \n"
      "              diffutils            (Edit)                               \n"
      "            [External/Edlib]                                            \n"
      "              edlib                (Edit)                               \n"
      "            [External/GABA]                                             \n"
      "              gaba-aband           (Gap-affine)                         \n"
      "            [External/KSW2]                                             \n"
      "              ksw2-extz2-sse       (Gap-affine)                         \n"
      "              ksw2-extd2-sse       (Gap-affine-2pieces)                 \n"
      "            [External/LV89]                                             \n"
      "              lv89                 (Edit)[score-only]                   \n"
      "            [External/Parasail]                                         \n"
      "              parasail-nw-stripped (Gap-affine)                         \n"
      "              parasail-nw-scan     (Gap-affine)                         \n"
      "              parasail-nw-diag     (Gap-affine)                         \n"
      "              parasail-nw-banded   (Gap-affine)[score-only]             \n"
      "            [External/SeqAn]                                            \n"
      "              seqan-edit           (Edit)                               \n"
      "              seqan-edit-bpm       (Edit)[score-only]                   \n"
      "              seqan-lineal         (Gap-linear)                         \n"
      "              seqan-affine         (Gap-affine)                         \n"
      "            [External/WFAlm]                                            \n"
      "              wfalm                (Gap-affine)                         \n"
      "              wfalm-lowmem         (Gap-affine)[low-mem]                \n"
#endif
      "        [Input & Output]                                                \n"
      "          --input|i <File>                                              \n"
      "          --output|o <File>                                             \n"
      "          --output-full <File>                                          \n"
      "        [Penalties & Span]                                              \n"
      "          --linear-penalties|p M,X,O                                    \n"
      "          --affine-penalties|g M,X,O,E                                  \n"
      "          --affine2p-penalties M,X,O1,E1,O2,E2                          \n"
      "          --ends-free P0,Pf,T0,Tf                                       \n"
      "        [Wavefront parameters]                                          \n"
      "          --wfa-score-only                                              \n"
      "          --wfa-memory-mode 'high'|'med'|'low'|'ultralow'               \n"
      "          --wfa-reduction 'none'|'adaptive'                             \n"
      "          --wfa-reduction-parameters  <P1>,<P2>                         \n"
      "            [Adaptive]                                                  \n"
      "              P1 = minimum-wavefront-length                             \n"
      "              P2 = maxumum-difference-distance                          \n"
      "        [Other Parameters]                                              \n"
      "          --bandwidth <INT>                                             \n"
#ifdef EXTERNAL_BENCHMARKS
      "          --ba-block-size <INT>                                         \n"
      "          --ksw2-approx-max-drop                                        \n"
      "          --ksw2-bandwidth <INT>                                        \n"
      "          --ksw2-zdrop <INT>                                            \n"
#endif
      "        [Misc]                                                          \n"
      "          --check|c 'correct'|'score'|'alignment'                       \n"
      "          --check-distance 'indel'|'edit'|'linear'|'affine'|'affine2p'  \n"
      "          --check-bandwidth <INT>                                       \n"
      "          --plot                                                        \n"
      "        [System]                                                        \n"
      "          --max-memory <bytes>                                          \n"
      "          --progress|P <integer>                                        \n"
      "          --help|h                                                      \n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* Input */
    { "algorithm", required_argument, 0, 'a' },
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "output-full", required_argument, 0, 800 },
    /* Penalties */
    { "linear-penalties", required_argument, 0, 'p' },
    { "affine-penalties", required_argument, 0, 'g' },
    { "affine2p-penalties", required_argument, 0, 900 },
    { "ends-free", required_argument, 0, 901 },
    /* Wavefront parameters */
    { "wfa-score-only", no_argument, 0, 1000 },
    { "wfa-memory-mode", required_argument, 0, 1001 },
    { "wfa-reduction", required_argument, 0, 1002 },
    { "wfa-reduction-parameters", required_argument, 0, 1003 },
    { "wfa-custom-match-funct", no_argument, 0, 1004 },
    { "wfa-max-memory", required_argument, 0, 1005 },
    /* Other alignment parameters */
    { "bandwidth", required_argument, 0, 2000 },
#ifdef EXTERNAL_BENCHMARKS
    { "ba-block-size", required_argument, 0, 2001 },
    { "ksw2-approx-max-drop", no_argument, 0, 2002 },
    { "ksw2-bandwidth", required_argument, 0, 2003 },
    { "ksw2-zdrop", required_argument, 0, 2004 },
#endif
    /* Misc */
    { "check", optional_argument, 0, 'c' },
    { "check-distance", required_argument, 0, 3001 },
    { "check-bandwidth", required_argument, 0, 3002 },
    { "plot", optional_argument, 0, 3003 },
    /* System */
    { "progress", required_argument, 0, 'P' },
    { "verbose", optional_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  if (argc <= 1) {
    usage();
    exit(0);
  }
  while (1) {
    c=getopt_long(argc,argv,"a:i:o:p:g:P:c:v::h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /*
     * Algorithm
     */
    case 'a': {
      /* Test bench */
      if (strcmp(optarg,"test")==0) {
        parameters.algorithm = alignment_test;
      // Indel
      } else if (strcmp(optarg,"indel-wfa")==0) {
        parameters.algorithm = alignment_indel_wavefront;
      // Edit
      } else if (strcmp(optarg,"edit-bpm")==0) {
        parameters.algorithm = alignment_edit_bpm;
      } else if (strcmp(optarg,"edit-dp")==0) {
        parameters.algorithm = alignment_edit_dp;
      } else if (strcmp(optarg,"edit-dp-banded")==0) {
        parameters.algorithm = alignment_edit_dp_banded;
      } else if (strcmp(optarg,"edit-wfa")==0) {
        parameters.algorithm = alignment_edit_wavefront;
      // Gap-Linear
      } else if (strcmp(optarg,"gap-linear-nw")==0 ||
                 strcmp(optarg,"gap-linear-dp")==0) {
        parameters.algorithm = alignment_gap_linear_nw;
      } else if (strcmp(optarg,"gap-linear-wfa")==0) {
        parameters.algorithm = alignment_gap_linear_wavefront;
      // Gap-Affine
      } else if (strcmp(optarg,"gap-affine-swg")==0 ||
                 strcmp(optarg,"gap-affine-dp")==0) {
        parameters.algorithm = alignment_gap_affine_swg;
      } else if (strcmp(optarg,"gap-affine-swg-banded")==0 ||
                 strcmp(optarg,"gap-affine-dp-banded")==0) {
        parameters.algorithm = alignment_gap_affine_swg_banded;
      } else if (strcmp(optarg,"gap-affine-wfa")==0) {
        parameters.algorithm = alignment_gap_affine_wavefront;
      // Gap-Affine 2-Pieces
      } else if (strcmp(optarg,"gap-affine2p-dp")==0) {
        parameters.algorithm = alignment_gap_affine2p_dp;
      } else if (strcmp(optarg,"gap-affine2p-wfa")==0) {
        parameters.algorithm = alignment_gap_affine2p_wavefront;
#ifdef EXTERNAL_BENCHMARKS
      /*
       * External Algorithm
       */
      // External (BitPal)
      } else if (strcmp(optarg,"bitpal-edit")==0) {
        parameters.algorithm = alignment_bitpal_edit;
      } else if (strcmp(optarg,"bitpal-scored")==0) {
        parameters.algorithm = alignment_bitpal_scored;
      // External (BlockAligner)
      } else if (strcmp(optarg,"block-aligner")==0) {
        parameters.algorithm = alignment_blockaligner;
      // External (Daligner)
      } else if (strcmp(optarg,"daligner")==0) {
        parameters.algorithm = alignment_daligner;
      // External (Diffutils)
      } else if (strcmp(optarg,"diffutils")==0) {
        parameters.algorithm = alignment_diffutils;
      // External (Edlib)
      } else if (strcmp(optarg,"edlib")==0) {
        parameters.algorithm = alignment_edlib;
      // External (Gaba)
      } else if (strcmp(optarg,"gaba-aband")==0) {
        parameters.algorithm = alignment_gaba_aband;
      // External (KSW2)
      } else if (strcmp(optarg,"ksw2-extz2-sse")==0) {
        parameters.algorithm = alignment_ksw2_extz2_sse;
      } else if (strcmp(optarg,"ksw2-extd2-sse")==0) {
        parameters.algorithm = alignment_ksw2_extd2_sse;
      // External (LV89)
      } else if (strcmp(optarg,"lv89")==0) {
        parameters.algorithm = alignment_lv89;
      // External (Parasail)
      } else if (strcmp(optarg,"parasail-nw-stripped")==0) {
        parameters.algorithm = alignment_parasail_nw_stripped;
      } else if (strcmp(optarg,"parasail-nw-scan")==0) {
        parameters.algorithm = alignment_parasail_nw_scan;
      } else if (strcmp(optarg,"parasail-nw-diag")==0) {
        parameters.algorithm = alignment_parasail_nw_diag;
      } else if (strcmp(optarg,"parasail-nw-banded")==0) {
        parameters.algorithm = alignment_parasail_nw_banded;
      // External (SeqAn)
      } else if (strcmp(optarg,"seqan-edit")==0) {
        parameters.algorithm = alignment_seqan_edit;
      } else if (strcmp(optarg,"seqan-edit-bpm")==0) {
        parameters.algorithm = alignment_seqan_edit_bpm;
      } else if (strcmp(optarg,"seqan-lineal")==0) {
        parameters.algorithm = alignment_seqan_lineal;
      } else if (strcmp(optarg,"seqan-affine")==0) {
        parameters.algorithm = alignment_seqan_affine;
      } else if (strcmp(optarg,"wfalm")==0) {
        parameters.algorithm = alignment_wfalm;
      } else if (strcmp(optarg,"wfalm-lowmem")==0) {
        parameters.algorithm = alignment_wfalm_lowmem;
#endif
      } else {
        fprintf(stderr,"Algorithm '%s' not recognized\n",optarg);
        exit(1);
      }
      break;
    }
    /*
     * Input/Output
     */
    case 'i':
      parameters.input_filename = optarg;
      break;
    case 'o':
      parameters.output_filename = optarg;
      break;
    case 800: // --output-full
      parameters.output_filename = optarg;
      parameters.output_full = true;
      break;
    /*
     * Penalties
     */
    case 'p': { // --linear-penalties M,X,O
      char* sentinel = strtok(optarg,",");
      parameters.linear_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.linear_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.linear_penalties.indel = atoi(sentinel);
      break;
    }
    case 'g': { // --affine-penalties M,X,O,E
      char* sentinel = strtok(optarg,",");
      parameters.affine_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_opening = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_extension = atoi(sentinel);
      break;
    }
    case 900: { // --affine2p-penalties M,X,O1,E1,O2,E2
      char* sentinel = strtok(optarg,",");
      parameters.affine2p_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_opening1 = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_extension1 = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_opening2 = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_extension2 = atoi(sentinel);
      break;
    }
    case 901: { // --ends-free P0,Pf,T0,Tf
      parameters.endsfree = true;
      char* sentinel = strtok(optarg,",");
      parameters.pattern_begin_free = atof(sentinel);
      sentinel = strtok(NULL,",");
      parameters.pattern_end_free = atof(sentinel);
      sentinel = strtok(NULL,",");
      parameters.text_begin_free = atof(sentinel);
      sentinel = strtok(NULL,",");
      parameters.text_end_free = atof(sentinel);
      break;
    }
    /*
     * Wavefront parameters
     */
    case 1000: // --wfa-score-only
      parameters.wfa_score_only = true;
      break;
    case 1001: // --wfa-memory-mode in {'high','med','low','ultralow'}
      if (strcmp(optarg,"high")==0) {
        parameters.wfa_memory_mode = wavefront_memory_high;
      } else if (strcmp(optarg,"med")==0) {
        parameters.wfa_memory_mode = wavefront_memory_med;
      } else if (strcmp(optarg,"low")==0) {
        parameters.wfa_memory_mode = wavefront_memory_low;
      } else if (strcmp(optarg,"ultralow")==0) {
        parameters.wfa_memory_mode = wavefront_memory_ultralow;
      } else {
        fprintf(stderr,"Option '--wfa-memory-mode' must be in {'high','med','low','ultralow'}\n");
        exit(1);
      }
      break;
    case 1002: // --wfa-reduction in {'none','adaptive'}
      if (strcmp(optarg,"none")==0) {
        parameters.wfa_reduction_type = wavefront_reduction_none;
      } else if (strcmp(optarg,"adaptive")==0) {
        parameters.wfa_reduction_type = wavefront_reduction_adaptive;
      } else {
        fprintf(stderr,"Option '--wf-reduction' must be in {'none','adaptive'}\n");
        exit(1);
      }
      break;
    case 1003: { // --wfa-reduction-parameters  <P1>,<P2>
      char* sentinel = strtok(optarg,",");
      const int p1 = atoi(sentinel);
      parameters.wfa_min_wf_length = p1;
      sentinel = strtok(NULL,",");
      const int p2 = atoi(sentinel);
      parameters.wfa_max_dist_th = p2;
      break;
    }
    case 1004: // --wfa-custom-match-funct
      parameters.wfa_match_funct = match_function;
      parameters.wfa_match_funct_arguments = &match_function_params;
      break;
    case 1005:
      parameters.wfa_max_memory = atol(optarg);
      break;
    /*
     * Other alignment parameters
     */
    case 2000: // --bandwidth
      parameters.bandwidth = atoi(optarg);
      break;
#ifdef EXTERNAL_BENCHMARKS
    case 2001: // --ba-block-size
      parameters.ba_block_size = atoi(optarg);
      break;
    case 2002: // --ksw2-approx-max-drop
      parameters.ksw2_approx_max__drop = true;
      break;
    case 2003: // --ksw2-bandwidth
      parameters.ksw2_bandwidth = atoi(optarg);
      break;
    case 2004: // --ksw2-zdrop
      parameters.ksw2_zdrop = atoi(optarg);
      break;
#endif
    /*
     * Misc
     */
    case 'c':
      if (optarg ==  NULL) { // default = score
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"display")==0) {
        parameters.check_display = true;
      } else if (strcasecmp(optarg,"correct")==0) {
        parameters.check_correct = true;
        parameters.check_score = false;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"score")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"alignment")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = true;
      } else {
        fprintf(stderr,"Option '--check' must be in {'correct','score','alignment'}\n");
        exit(1);
      }
      break;
    case 3001: // --check-distance in {'indel','edit','linear','affine','affine2p'}
      if (strcasecmp(optarg,"indel")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_INDEL;
      } else if (strcasecmp(optarg,"edit")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_EDIT;
      } else if (strcasecmp(optarg,"linear")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_LINEAR;
      } else if (strcasecmp(optarg,"affine")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE;
      } else if (strcasecmp(optarg,"affine2p")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE2P;
      } else {
        fprintf(stderr,"Option '--check-distance' must be in {'indel','edit','linear','affine','affine2p'}\n");
        exit(1);
      }
      break;
    case 3002: // --check-bandwidth
      parameters.check_bandwidth = atoi(optarg);
      break;
    case 3003: // --plot
      parameters.plot = (optarg==NULL) ? 1000 : atoi(optarg);
      break;
    /*
     * System
     */
    case 'P':
      parameters.progress = atoi(optarg);
      break;
    case 'v':
      if (optarg==NULL) {
        parameters.verbose = 1;
      } else {
        parameters.verbose = atoi(optarg);
        if (parameters.verbose < 0 || parameters.verbose > 3) {
          fprintf(stderr,"Option '--verbose' must be in {0,1,2,3}\n");
          exit(1);
        }
      }
      break;
    case 'h':
      usage();
      exit(1);
    // Other
    case '?': default:
      fprintf(stderr,"Option not recognized \n");
      exit(1);
    }
  }
  // General checks
  if (parameters.algorithm!=alignment_test && parameters.input_filename==NULL) {
    fprintf(stderr,"Option --input is required \n");
    exit(1);
  }
  // Check 'ends-free' parameter
  if (parameters.endsfree) {
    switch (parameters.algorithm) {
      case alignment_gap_affine_swg:
        parameters.algorithm = alignment_gap_affine_swg_endsfree;
        break;
      case alignment_indel_wavefront:
      case alignment_edit_wavefront:
      case alignment_gap_linear_wavefront:
      case alignment_gap_affine_wavefront:
      case alignment_gap_affine2p_wavefront:
        break;
      default:
        fprintf(stderr,"Ends-free variant not implemented for the selected algorithm\n");
        exit(1);
        break;
    }
  }
  // Check 'bandwidth' parameter
  switch (parameters.algorithm) {
    case alignment_edit_dp_banded:
    case alignment_gap_affine_swg_banded:
    case alignment_parasail_nw_banded:
      if (parameters.bandwidth == -1) {
        fprintf(stderr,"Parameter 'bandwidth' has to be provided for banded algorithms\n");
        exit(1);
      }
      break;
    default:
      if (parameters.bandwidth == -1) {
        fprintf(stderr,"Parameter 'bandwidth' has no effect with the selected algorithm\n");
        exit(1);
      }
      break;
  }
}
int main(int argc,char* argv[]) {
  // Parsing command-line options
  parse_arguments(argc,argv);
  // Select option
  if (parameters.algorithm == alignment_test) {
    align_pairwise_test();
  } else {
    align_benchmark();
  }
}
