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
#include <omp.h>

#include "utils/commons.h"
#include "utils/sequence_buffer.h"
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
  // I/O
  char *input_filename;
  char *output_filename;
  bool output_full;
  // I/O internals
  FILE* input_file;
  char* line1;
  char* line2;
  size_t line1_allocated;
  size_t line2_allocated;
  FILE* output_file;
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
  wf_heuristic_strategy wfa_heuristic;
  int wfa_heuristic_p1;
  int wfa_heuristic_p2;
  int wfa_heuristic_p3;
  wavefront_memory_t wfa_memory_mode;
  alignment_match_funct_t wfa_match_funct;
  void* wfa_match_funct_arguments;
  uint64_t wfa_max_memory;
  // Other algorithms parameters
  int bandwidth;
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
  int num_threads;
  int batch_size;
  int progress;
  int verbose;
} benchmark_args;
benchmark_args parameters = {
  // Algorithm
  .algorithm = alignment_test,
  // I/O
  .input_filename = NULL,
  .output_filename = NULL,
  .output_full = false,
  .output_file = NULL,
  // I/O internals
  .input_file = NULL,
  .line1 = NULL,
  .line2 = NULL,
  .line1_allocated = 0,
  .line2_allocated = 0,
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
  .wfa_heuristic = wf_heuristic_none,
  .wfa_heuristic_p1 = -1,
  .wfa_heuristic_p2 = -1,
  .wfa_heuristic_p3 = -1,
  .wfa_memory_mode = wavefront_memory_high,
  .wfa_match_funct = NULL,
  .wfa_match_funct_arguments = NULL,
  .wfa_max_memory = UINT64_MAX,
  // Other algorithms parameters
  .bandwidth = -1,
  // Misc
  .check_bandwidth = -1,
  .check_display = false,
  .check_correct = false,
  .check_score = false,
  .check_alignments = false,
  .check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE,
  .plot = 0,
  // System
  .num_threads = 1,
  .batch_size = 10000,
  .progress = 10000,
  .verbose = 0,
};

/*
 * Benchmark UTest
 */
void align_pairwise_test() {
  // Patters & Texts
  char * pattern = "GATTACA";
  char * text = "GATCACTA";

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
  attributes.heuristic.strategy = wf_heuristic_none;
  attributes.heuristic.min_wavefront_length = 256;
  attributes.heuristic.max_distance_threshold = 4096;
  attributes.heuristic.steps_between_cutoffs = 10;
  attributes.alignment_scope = compute_alignment; // compute_score
  attributes.memory_mode = wavefront_memory_med;
  attributes.alignment_form.span = (endsfree) ? alignment_endsfree : alignment_end2end;
  attributes.alignment_form.pattern_begin_free = pattern_begin_free;
  attributes.alignment_form.pattern_end_free = pattern_end_free;
  attributes.alignment_form.text_begin_free = text_begin_free;
  attributes.alignment_form.text_end_free = text_end_free;
  attributes.plot_params.plot_enabled = false;
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  // Align
  wavefront_align(wf_aligner,
      pattern,strlen(pattern),text,strlen(text));
  // CIGAR
  fprintf(stderr,">> WFA2");
  cigar_print_pretty(stderr,
      pattern,strlen(pattern),text,strlen(text),
      &wf_aligner->cigar,wf_aligner->mm_allocator);
  fprintf(stderr,"SCORE: %d \n",cigar_score_gap_affine(&wf_aligner->cigar,&affine_penalties));
  // Plot
  if (attributes.plot_params.plot_enabled) {
    FILE* const wf_plot = fopen("test.wfa","w");
    wavefront_plot_print(wf_plot,wf_aligner);
    fclose(wf_plot);
  }
  // Free
  wavefront_aligner_delete(wf_aligner);
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
wavefront_aligner_t* align_input_configure_wavefront(
    align_input_t* const align_input,
    mm_allocator_t* const mm_allocator) {
  // Set attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.memory_mode = parameters.wfa_memory_mode;
  attributes.mm_allocator = mm_allocator;
  if (parameters.wfa_score_only) {
    attributes.alignment_scope = compute_score;
  }
  // WF-Heuristic
  switch (parameters.wfa_heuristic) {
    case wf_heuristic_none:
      attributes.heuristic.strategy = wf_heuristic_none;
      break;
    case wf_heuristic_banded_static:
      attributes.heuristic.strategy = wf_heuristic_banded_static;
      attributes.heuristic.min_k = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_k = parameters.wfa_heuristic_p2;
      break;
    case wf_heuristic_banded_adaptive:
      attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
      attributes.heuristic.min_k = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_k = parameters.wfa_heuristic_p2;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p3;
      break;
    case wf_heuristic_wfadaptive:
      attributes.heuristic.strategy = wf_heuristic_wfadaptive;
      attributes.heuristic.min_wavefront_length = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_distance_threshold = parameters.wfa_heuristic_p2;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p3;
      break;
    case wf_heuristic_xdrop:
      attributes.heuristic.strategy = wf_heuristic_xdrop;
      attributes.heuristic.xdrop = parameters.wfa_heuristic_p1;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p2;
      break;
    case wf_heuristic_zdrop:
      attributes.heuristic.strategy = wf_heuristic_zdrop;
      attributes.heuristic.zdrop = parameters.wfa_heuristic_p1;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p2;
      break;
    default:
      break;
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
void align_input_configure_global(
    align_input_t* const align_input) {
  // Clear
  benchmark_align_input_clear(align_input);
  // Penalties
  align_input->linear_penalties = parameters.linear_penalties;
  align_input->affine_penalties = parameters.affine_penalties;
  align_input->affine2p_penalties = parameters.affine2p_penalties;
  // Alignment form
  align_input->ends_free = parameters.endsfree;
  // Output
  align_input->output_file = parameters.output_file;
  align_input->output_full = parameters.output_full;
  // MM
  align_input->mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  // WFA
  if (align_benchmark_is_wavefront(parameters.algorithm)) {
    align_input->wf_aligner = align_input_configure_wavefront(align_input,align_input->mm_allocator);
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
void align_input_configure_local(
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
void align_benchmark_free(
    align_input_t* const align_input) {
  if (align_input->wf_aligner) wavefront_aligner_delete(align_input->wf_aligner);
  mm_allocator_delete(align_input->mm_allocator);
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
    const int seqs_processed) {
  const uint64_t time_elapsed_alg = timer_get_current_total_ns(&parameters.timer_global);
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
  if (parameters.num_threads == 1) {
    fprintf(stderr,"  => Time.Alignment    ");
    timer_print(stderr,&align_input->timer,&parameters.timer_global);
  } else {
    for (int i=0;i<parameters.num_threads;++i) {
      fprintf(stderr,"  => Time.Alignment.Thread.%0d    ",i);
      timer_print(stderr,&align_input[i].timer,&parameters.timer_global);
    }
  }
  // Print Stats
  const bool checks_enabled =
      parameters.check_display || parameters.check_correct ||
      parameters.check_score || parameters.check_alignments;
  if (checks_enabled && parameters.num_threads==1) {
    const bool print_wf_stats = (parameters.algorithm == alignment_gap_affine_wavefront);
    benchmark_print_stats(stderr,align_input,print_wf_stats);
  }
}
void align_benchmark_plot_wf(
    align_input_t* const align_input,
    const int seq_id) {
  // Setup filename
  char filename[500];
  if (parameters.output_filename != NULL) {
    sprintf(filename,"%s.%03d.wfa",parameters.output_filename,seq_id);
  } else {
    sprintf(filename,"%s.%03d.wfa",parameters.input_filename,seq_id);
  }
  // Open file
  FILE* const wf_plot = fopen(filename,"w");
  wavefront_plot_print(wf_plot,align_input->wf_aligner);
  fclose(wf_plot);
}
/*
 * Benchmark
 */
void align_benchmark_run_algorithm(
    align_input_t* const align_input) {
  // Sequence-dependent configuration
  align_input_configure_local(align_input);
  // Select algorithm
  switch (parameters.algorithm) {
    // Indel
    case alignment_indel_wavefront:
      benchmark_indel_wavefront(align_input);
      break;
    // Edit
    case alignment_edit_bpm:
      benchmark_edit_bpm(align_input);
      break;
    case alignment_edit_dp:
      benchmark_edit_dp(align_input);
      break;
    case alignment_edit_dp_banded:
      benchmark_edit_dp_banded(align_input,parameters.bandwidth);
      break;
    case alignment_edit_wavefront:
      benchmark_edit_wavefront(align_input);
      break;
    // Gap-linear
    case alignment_gap_linear_nw:
      benchmark_gap_linear_nw(align_input,&parameters.linear_penalties);
      break;
    case alignment_gap_linear_wavefront:
      benchmark_gap_linear_wavefront(align_input,&parameters.linear_penalties);
      break;
    // Gap-affine
    case alignment_gap_affine_swg:
      benchmark_gap_affine_swg(align_input,&parameters.affine_penalties);
      break;
    case alignment_gap_affine_swg_endsfree:
      benchmark_gap_affine_swg_endsfree(
          align_input,&parameters.affine_penalties);
      break;
    case alignment_gap_affine_swg_banded:
      benchmark_gap_affine_swg_banded(align_input,
          &parameters.affine_penalties,parameters.bandwidth);
      break;
    case alignment_gap_affine_wavefront:
      benchmark_gap_affine_wavefront(align_input,&parameters.affine_penalties);
      break;
    // Gap-affine 2p
    case alignment_gap_affine2p_dp:
      benchmark_gap_affine2p_dp(align_input,&parameters.affine2p_penalties);
      break;
    case alignment_gap_affine2p_wavefront:
      benchmark_gap_affine2p_wavefront(align_input,&parameters.affine2p_penalties);
      break;
    default:
      fprintf(stderr,"Algorithm not implemented\n");
      exit(1);
      break;
  }
}
void align_benchmark_sequential() {
  // PROFILE
  timer_reset(&parameters.timer_global);
  timer_start(&parameters.timer_global);
  // I/O files
  parameters.input_file = fopen(parameters.input_filename, "r");
  if (parameters.input_file == NULL) {
    fprintf(stderr,"Input file '%s' couldn't be opened\n",parameters.input_filename);
    exit(1);
  }
  if (parameters.output_filename != NULL) {
    parameters.output_file = fopen(parameters.output_filename, "w");
  }
  // Global configuration
  align_input_t align_input;
  align_input_configure_global(&align_input);
  // Read-align loop
  int seqs_processed = 0, progress = 0;
  while (true) {
    // Read input sequence-pair
    const bool input_read = align_benchmark_read_input(
        parameters.input_file,&parameters.line1,&parameters.line2,
        &parameters.line1_allocated,&parameters.line2_allocated,
        seqs_processed,&align_input);
    if (!input_read) break;
    // Execute the selected algorithm
    align_benchmark_run_algorithm(&align_input);
    // Update progress
    ++seqs_processed;
    if (++progress == parameters.progress) {
      progress = 0;
      align_benchmark_print_progress(seqs_processed);
    }
    // DEBUG
    // mm_allocator_print(stderr,align_input.mm_allocator,true);
    // Plot
    if (parameters.plot > 0) align_benchmark_plot_wf(&align_input,seqs_processed);
  }
  // Print benchmark results
  timer_stop(&parameters.timer_global);
  align_benchmark_print_results(&align_input,seqs_processed,true);
  // Free
  align_benchmark_free(&align_input);
  fclose(parameters.input_file);
  if (parameters.output_file) fclose(parameters.output_file);
  free(parameters.line1);
  free(parameters.line2);
}
void align_benchmark_parallel() {
  // PROFILE
  timer_reset(&parameters.timer_global);
  timer_start(&parameters.timer_global);
  // Open input file
  parameters.input_file = fopen(parameters.input_filename, "r");
  if (parameters.input_file == NULL) {
    fprintf(stderr,"Input file '%s' couldn't be opened\n",parameters.input_filename);
    exit(1);
  }
  if (parameters.output_filename != NULL) {
    parameters.output_file = fopen(parameters.output_filename, "w");
  }
  // Global configuration
  align_input_t align_input[parameters.num_threads];
  for (int tid=0;tid<parameters.num_threads;++tid) {
    align_input_configure_global(align_input+tid);
  }
  // Read-align loop
  sequence_buffer_t* const sequence_buffer = sequence_buffer_new(2*parameters.batch_size,100);
  int seqs_processed = 0, progress = 0, seqs_batch = 0;
  while (true) {
    // Read batch-input sequence-pair
    sequence_buffer_clear(sequence_buffer);
    for (seqs_batch=0;seqs_batch<parameters.batch_size;++seqs_batch) {
      const bool seqs_pending = align_benchmark_read_input(
          parameters.input_file,&parameters.line1,&parameters.line2,
          &parameters.line1_allocated,&parameters.line2_allocated,
          seqs_processed,align_input);
      if (!seqs_pending) break;
      // Add pair pattern-text
      sequence_buffer_add_pair(sequence_buffer,
          align_input->pattern,align_input->pattern_length,
          align_input->text,align_input->text_length);
    }
    if (seqs_batch == 0) break;
    // Parallel processing of the sequences batch
    #pragma omp parallel num_threads(parameters.num_threads)
    {
      int tid = omp_get_thread_num();
      #pragma omp for
      for (int seq_idx=0;seq_idx<seqs_batch;++seq_idx) {
        // Configure sequence
        sequence_offset_t* const offset = sequence_buffer->offsets + seq_idx;
        align_input[tid].sequence_id = seqs_processed;
        align_input[tid].pattern = sequence_buffer->buffer + offset->pattern_offset;
        align_input[tid].pattern_length = offset->pattern_length;
        align_input[tid].text = sequence_buffer->buffer + offset->text_offset;
        align_input[tid].text_length = offset->text_length;
        // Execute the selected algorithm
        align_benchmark_run_algorithm(align_input+tid);
      }
    }
    // Update progress
    seqs_processed += seqs_batch;
    progress += seqs_batch;
    if (progress >= parameters.progress) {
      progress -= parameters.progress;
      align_benchmark_print_progress(seqs_processed);
    }
    // DEBUG
    // mm_allocator_print(stderr,align_input.mm_allocator,true);
  }
  // Print benchmark results
  timer_stop(&parameters.timer_global);
  align_benchmark_print_results(align_input,seqs_processed,true);
  // Free
  for (int tid=0;tid<parameters.num_threads;++tid) {
    align_benchmark_free(align_input+tid);
  }
  sequence_buffer_delete(sequence_buffer);
  fclose(parameters.input_file);
  if (parameters.output_file) fclose(parameters.output_file);
  free(parameters.line1);
  free(parameters.line2);
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
      "            [Gap-affine (Smith-Waterman-Gotoh)]                         \n"
      "              gap-affine-swg                                            \n"
      "              gap-affine-swg-banded                                     \n"
      "              gap-affine-wfa                                            \n"
      "            [Gap-affine-2pieces (Concave 2-pieces)]                     \n"
      "              gap-affine2p-dp                                           \n"
      "              gap-affine2p-wfa                                          \n"
      "        [Input & Output]                                                \n"
      "          --input|i <File>                                              \n"
      "          --output|o <File>                                             \n"
      "          --output-full <File>                                          \n"
      "        [Penalties & Span]                                              \n"
      "          --linear-penalties|p M,X,I                                    \n"
      "          --affine-penalties|g M,X,O,E                                  \n"
      "          --affine2p-penalties M,X,O1,E1,O2,E2                          \n"
      "          --ends-free P0,Pf,T0,Tf                                       \n"
      "        [Wavefront parameters]                                          \n"
      "          --wfa-score-only                                              \n"
      "          --wfa-memory-mode 'high'|'med'|'low'                          \n"
      "          --wfa-heuristic <Strategy>                                    \n"
      "          --wfa-heuristic-parameters  <P1>,<P2>[,<P3>]                  \n"
      "            [Strategy='banded-static']                                  \n"
      "              P1 = minimum-diagonal-band (e.g., -100)                   \n"
      "              P2 = maximum-diagonal-band (e.g., +100)                   \n"
      "            [Strategy='banded-adaptive']                                \n"
      "              P1 = minimum-diagonal-band (e.g., -100)                   \n"
      "              P2 = maximum-diagonal-band (e.g., +100)                   \n"
      "              P3 = steps-between-cutoffs                                \n"
      "            [Strategy='wfa-adaptive']                                   \n"
      "              P1 = minimum-wavefront-length                             \n"
      "              P2 = maximum-difference-distance                          \n"
      "              P3 = steps-between-cutoffs                                \n"
      "            [Strategy='xdrop']                                          \n"
      "              P1 = x-drop                                               \n"
      "              P2 = steps-between-cutoffs                                \n"
      "            [Strategy='zdrop']                                          \n"
      "              P1 = z-drop                                               \n"
      "              P2 = steps-between-cutoffs                                \n"
      "          --wfa-max-memory <Bytes>                                      \n"
      "        [Other Parameters]                                              \n"
      "          --bandwidth <INT>                                             \n"
      "        [Misc]                                                          \n"
      "          --check|c 'correct'|'score'|'alignment'                       \n"
      "          --check-distance 'indel'|'edit'|'linear'|'affine'|'affine2p'  \n"
      "          --check-bandwidth <INT>                                       \n"
      "          --plot                                                        \n"
      "        [System]                                                        \n"
      "          --num-threads|t <INT>                                         \n"
      "          --batch-size <INT>                                            \n"
    //"          --progress|P <INT>                                            \n"
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
    { "wfa-heuristic", required_argument, 0, 1002 },
    { "wfa-heuristic-parameters", required_argument, 0, 1003 },
    { "wfa-custom-match-funct", no_argument, 0, 1004 },
    { "wfa-max-memory", required_argument, 0, 1005 },
    /* Other alignment parameters */
    { "bandwidth", required_argument, 0, 2000 },
    /* Misc */
    { "check", optional_argument, 0, 'c' },
    { "check-distance", required_argument, 0, 3001 },
    { "check-bandwidth", required_argument, 0, 3002 },
    { "plot", optional_argument, 0, 3003 },
    /* System */
    { "num-threads", required_argument, 0, 't' },
    { "batch-size", required_argument, 0, 4000 },
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
    c=getopt_long(argc,argv,"a:i:o:p:g:P:c:v::t:h",long_options,&option_index);
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
    case 'p': { // --linear-penalties M,X,I
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
    case 1001: // --wfa-memory-mode in {'high','med','low'}
      if (strcmp(optarg,"high")==0) {
        parameters.wfa_memory_mode = wavefront_memory_high;
      } else if (strcmp(optarg,"med")==0) {
        parameters.wfa_memory_mode = wavefront_memory_med;
      } else if (strcmp(optarg,"low")==0) {
        parameters.wfa_memory_mode = wavefront_memory_low;
      } else {
        fprintf(stderr,"Option '--wfa-memory-mode' must be in {'high','med','low'}\n");
        exit(1);
      }
      break;
    case 1002: // --wfa-heuristic in {'none'|'banded-static'|'banded-adaptive'|'wfa-adaptive'|'xdrop'|'zdrop'}
      if (strcmp(optarg,"none")==0) {
        parameters.wfa_heuristic = wf_heuristic_none;
      } else if (strcmp(optarg,"banded-static")==0 || strcmp(optarg,"banded")==0) {
        parameters.wfa_heuristic = wf_heuristic_banded_static;
      } else if (strcmp(optarg,"banded-adaptive")==0) {
        parameters.wfa_heuristic = wf_heuristic_banded_adaptive;
      } else if (strcmp(optarg,"wfa-adaptive")==0) {
        parameters.wfa_heuristic = wf_heuristic_wfadaptive;
      } else if (strcmp(optarg,"xdrop")==0) {
        parameters.wfa_heuristic = wf_heuristic_xdrop;
      } else if (strcmp(optarg,"zdrop")==0) {
        parameters.wfa_heuristic = wf_heuristic_zdrop;
      } else {
        fprintf(stderr,"Option '--wf-heuristic' must be in {'none'|'banded-static'|'banded-adaptive'|'wfa-adaptive'|'xdrop'|'zdrop'}\n");
        exit(1);
      }
      break;
    case 1003: { // --wfa-heuristic-parameters  <P1>,<P2>[,<P3>]
      char* sentinel = strtok(optarg,",");
      const int p1 = atoi(sentinel);
      parameters.wfa_heuristic_p1 = p1;
      sentinel = strtok(NULL,",");
      const int p2 = atoi(sentinel);
      parameters.wfa_heuristic_p2 = p2;
      sentinel = strtok(NULL,",");
      if (sentinel != NULL) {
        const int p3 = atoi(sentinel);
        parameters.wfa_heuristic_p3 = p3;
      }
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
    case 't': // --num-threads
      parameters.num_threads = atoi(optarg);
      break;
    case 4000: // --batch-size
      parameters.batch_size = atoi(optarg);
      break;
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
  // Checks general
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
      if (parameters.bandwidth == -1) {
        fprintf(stderr,"Parameter 'bandwidth' has to be provided for banded algorithms\n");
        exit(1);
      }
      break;
    default:
      if (parameters.bandwidth != -1) {
        fprintf(stderr,"Parameter 'bandwidth' has no effect with the selected algorithm\n");
        exit(1);
      }
      break;
  }
  // Check 'wfa-heuristic'
  switch (parameters.wfa_heuristic) {
    case wf_heuristic_banded_static:
    case wf_heuristic_xdrop:
    case wf_heuristic_zdrop:
      if (parameters.wfa_heuristic_p1 == -1 ||
          parameters.wfa_heuristic_p2 == -1) {
        fprintf(stderr,"Heuristic requires parameters '--wfa-heuristic-parameters' <P1>,<P2>\n");
        exit(1);
      }
      break;
    case wf_heuristic_banded_adaptive:
    case wf_heuristic_wfadaptive:
      if (parameters.wfa_heuristic_p1 == -1 ||
          parameters.wfa_heuristic_p2 == -1 ||
          parameters.wfa_heuristic_p3 == -1) {
        fprintf(stderr,"Heuristic requires parameters '--wfa-heuristic-parameters' <P1>,<P2>,<P3>\n");
        exit(1);
      }
      break;
    default:
      break;
  }
  // Checks parallel
  if (parameters.num_threads > 1) {
    if (parameters.plot > 0) {
      fprintf(stderr,"Parameter 'plot' disabled for parallel executions\n");
      parameters.plot = 0;
    }
  }
}
int main(int argc,char* argv[]) {
  // Parsing command-line options
  parse_arguments(argc,argv);
  // Select option
  if (parameters.algorithm == alignment_test) {
    align_pairwise_test();
  } else {
    // Execute benchmark
    if (parameters.num_threads == 1) {
      align_benchmark_sequential();
    } else {
      align_benchmark_parallel();
    }
  }
}
