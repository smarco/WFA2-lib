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
 * DESCRIPTION: Wavefront Alignments Algorithms benchmarking tool
 */

#include "utils/commons.h"
#include "system/profiler_timer.h"

#include "alignment/score_matrix.h"
#include "edit/edit_dp.h"
#include "gap_lineal/nw.h"
#include "gap_affine/swg.h"

#include "benchmark/benchmark_edit.h"
#include "benchmark/benchmark_gap_lineal.h"
#include "benchmark/benchmark_gap_affine.h"
#include "benchmark/benchmark_gap_affine2p.h"

#include "gap_affine2p/affine2p_penalties.h"
#include "gap_affine2p/affine2p_matrix.h"
#include "gap_affine2p/affine2p_dp.h"

#include "wavefront/wavefront_align.h"

/*
 * Algorithms
 */
typedef enum {
  alignment_test,
  alignment_edit_dp,
  alignment_edit_dp_banded,
  alignment_edit_wavefront,
  alignment_gap_lineal_nw,
  alignment_gap_lineal_wavefront,
  alignment_gap_affine_swg,
  alignment_gap_affine_swg_endsfree,
  alignment_gap_affine_swg_banded,
  alignment_gap_affine_wavefront,
  alignment_gap_affine2p_dp,
  alignment_gap_affine2p_wavefront
} alignment_algorithm_type;

/*
 * Generic parameters
 */
typedef struct {
  // Algorithm
  alignment_algorithm_type algorithm;
  // Input
  char *input_filename;
  char *output_filename;
  // Penalties
  lineal_penalties_t lineal_penalties;
  affine_penalties_t affine_penalties;
  affine2p_penalties_t affine2p_penalties;
  // Wavefront parameters
  bool score_only;
  wavefront_reduction_type reduction_type;
  int min_wavefront_length;
  int max_distance_threshold;
  bool low_memory;
  bool endsfree;
  int pattern_begin_free;
  int text_begin_free;
  int pattern_end_free;
  int text_end_free;
  // Misc
  int bandwidth;
  bool check_correct;
  bool check_score;
  bool check_alignments;
  int check_metric;
  int check_bandwidth;
  int plot;
  // Profile
  profiler_timer_t timer_global;
  // System
  uint64_t max_memory;
  int progress;
  bool verbose;
} benchmark_args;
benchmark_args parameters = {
  // Algorithm
  .algorithm = alignment_test,
  // Input
  .input_filename=NULL,
  .output_filename=NULL,
  // Penalties
  .lineal_penalties = {
      .match = 0,
      .mismatch = 4,
      .insertion = 2,
      .deletion  = 2,
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
  // Wavefront parameters
  .score_only = false,
  .reduction_type = wavefront_reduction_none,
  .min_wavefront_length = 10,
  .max_distance_threshold = 50,
  .low_memory = false,
  .endsfree = false,
  .pattern_begin_free = 0,
  .text_begin_free = 0,
  .pattern_end_free = 0,
  .text_end_free = 0,
  // Misc
  .bandwidth = 10,
  .check_bandwidth = -1,
  .check_correct = false,
  .check_score = false,
  .check_alignments = false,
  .check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE,
  .plot = 0,
  // System
  .max_memory = UINT64_MAX,
  .progress = 10000,
  .verbose = false
};

/*
 * Benchmark UTest
 */
void align_pairwise_test() {
  // Patters & Texts
//  char * pattern = "GATTACA";
//  char * text = "GATCACTA";
    char * pattern = "GCAGAGAATTACGACCGGCTCGCTGAATTGCGAAG";
    char * text = "GCGAGAATTACGACCGGCTCGCTGAATTGCGCGAAG";

  // MMAllocator
  mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  // Penalties
  affine_penalties_t affine_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening = 6,
      .gap_extension = 2,
  };
  affine2p_penalties_t affine2p_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening1 = 6,
      .gap_extension1 = 2,
      .gap_opening2 = 7000,
      .gap_extension2 = 3000,
  };
  // Ends
  const int pattern_begin_free = 0;
  const int pattern_end_free = 0;
  const int text_begin_free = 0;
  const int text_end_free = 0;
//  const int pattern_begin_free = 20;
//  const int pattern_end_free = 20;
//  const int text_begin_free = 20;
//  const int text_end_free = 20;
  const bool endsfree =
      pattern_begin_free>0 ||
      pattern_end_free>0   ||
      text_begin_free>0    ||
      text_end_free>0;

//  /*
//   * SWG
//   */
//  // Allocate
//  affine_matrix_t matrix;
//  affine_matrix_allocate(&matrix,strlen(pattern)+1,strlen(text)+1,mm_allocator);
//  cigar_t cigar;
//  cigar_allocate(&cigar,strlen(pattern)+strlen(text),mm_allocator);
//  // Align
//  swg_compute_endsfree(
//      &matrix,&affine_penalties,
//      pattern,strlen(pattern),text,strlen(text),
//      pattern_begin_free,pattern_end_free,
//      text_begin_free,text_end_free,&cigar);
//  cigar_print_pretty(stderr,
//      pattern,strlen(pattern),text,strlen(text),
//      &cigar,mm_allocator);
//  fprintf(stderr,"SCORE: %d \n",cigar_score_gap_affine(&cigar,&affine_penalties));
//  // Free
//  affine_matrix_free(&matrix,mm_allocator);
//  cigar_free(&cigar);

//  /*
//   * SWG 2P
//   */
//  // Allocate
//  affine2p_matrix_t matrix;
//  affine2p_matrix_allocate(&matrix,strlen(pattern)+1,strlen(text)+1,mm_allocator);
//  cigar_t cigar;
//  cigar_allocate(&cigar,strlen(pattern)+strlen(text),mm_allocator);
//  // Align
//  affine2p_dp_compute(
//      &matrix,&affine2p_penalties,
//      pattern,strlen(pattern),
//      text,strlen(text),&cigar);
//  cigar_print_pretty(stderr,
//      pattern,strlen(pattern),text,strlen(text),
//      &cigar,mm_allocator);
//  fprintf(stderr,"SCORE: %d \n",cigar_score_gap_affine2p(&cigar,&affine2p_penalties));
//  // Free
//  affine2p_matrix_free(&matrix,mm_allocator);
//  cigar_free(&cigar);

  /*
   * Gap-Affine
   */
  // Allocate
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine;
  attributes.affine_penalties = affine_penalties;
  // attributes.affine2p_penalties = affine2p_penalties;
  attributes.reduction.reduction_strategy = wavefront_reduction_adaptive; //wavefront_reduction_none; // wavefront_reduction_adaptive
  attributes.reduction.min_wavefront_length = 10;
  attributes.reduction.max_distance_threshold = 50;
  attributes.alignment_scope = compute_alignment; // compute_score
  attributes.low_memory = false;
  attributes.alignment_form.span = (endsfree) ? alignment_endsfree : alignment_end2end;
  attributes.alignment_form.pattern_begin_free = pattern_begin_free;
  attributes.alignment_form.pattern_end_free = pattern_end_free;
  attributes.alignment_form.text_begin_free = text_begin_free;
  attributes.alignment_form.text_end_free = text_end_free;
  attributes.plot_params.plot_enabled = true;
  attributes.mm_allocator = mm_allocator;
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  // Align
  wavefront_align(wf_aligner,
      pattern,strlen(pattern),text,strlen(text));
  // CIGAR
  cigar_print_pretty(stderr,
      pattern,strlen(pattern),text,strlen(text),
      &wf_aligner->cigar,mm_allocator);
  fprintf(stderr,"SCORE: %d \n",cigar_score_gap_affine2p(&wf_aligner->cigar,&affine2p_penalties));
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
 * Configuration
 */
wavefront_aligner_t* align_benchmark_configure_wf(
    align_input_t* const align_input,
    mm_allocator_t* const mm_allocator) {
  // Set attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.low_memory = parameters.low_memory;
  attributes.mm_allocator = mm_allocator;
  if (parameters.score_only) {
    attributes.alignment_scope = compute_score;
  }
  // WF-Reduction
  if (parameters.reduction_type == wavefront_reduction_adaptive) {
    attributes.reduction.reduction_strategy = wavefront_reduction_adaptive;
    attributes.reduction.min_wavefront_length = parameters.min_wavefront_length;
    attributes.reduction.max_distance_threshold = parameters.max_distance_threshold;
  } else {
    attributes.reduction.reduction_strategy = wavefront_reduction_none;
    attributes.reduction.min_wavefront_length = -1;
    attributes.reduction.max_distance_threshold = -1;
  }
  // Select flavor
  switch (parameters.algorithm) {
    case alignment_edit_wavefront:
      attributes.distance_metric = edit;
      break;
    case alignment_gap_lineal_wavefront:
      attributes.distance_metric = gap_lineal;
      attributes.lineal_penalties = parameters.lineal_penalties;
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
  if (parameters.endsfree) {
    attributes.alignment_form.span = alignment_endsfree;
    attributes.alignment_form.pattern_begin_free = parameters.pattern_begin_free;
    attributes.alignment_form.text_begin_free = parameters.text_begin_free;
    attributes.alignment_form.pattern_end_free = parameters.pattern_end_free;
    attributes.alignment_form.text_end_free = parameters.text_end_free;
  } else {
    attributes.alignment_form.span = alignment_end2end;
  }
  // Misc
  attributes.plot_params.plot_enabled = (parameters.plot > 0);
  attributes.plot_params.resolution_points = parameters.plot;
  attributes.system.verbose = parameters.verbose;
  attributes.system.max_memory_used = parameters.max_memory;
  // Allocate
  return wavefront_aligner_new(&attributes);
}
void align_benchmark_configure(
    align_input_t* const align_input) {
  // Clear
  benchmark_align_input_clear(align_input);
  // Output
  if (parameters.output_filename == NULL) {
    align_input->output_file = NULL;
  } else {
    align_input->output_file = fopen(parameters.output_filename, "w");
  }
  // MM
  align_input->mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  align_input->wf_aligner = align_benchmark_configure_wf(align_input,align_input->mm_allocator);
  // PROFILE/STATS
  timer_reset(&align_input->timer);
  // DEBUG
  align_input->debug_flags = 0;
  align_input->debug_flags |= parameters.check_metric;
  if (parameters.check_correct) align_input->debug_flags |= ALIGN_DEBUG_CHECK_CORRECT;
  if (parameters.check_score) align_input->debug_flags |= ALIGN_DEBUG_CHECK_SCORE;
  if (parameters.check_alignments) align_input->debug_flags |= ALIGN_DEBUG_CHECK_ALIGNMENT;
  align_input->check_lineal_penalties = &parameters.lineal_penalties;
  align_input->check_affine_penalties = &parameters.affine_penalties;
  align_input->check_bandwidth = parameters.check_bandwidth;
  align_input->verbose = parameters.verbose;
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
  if (parameters.check_correct || parameters.check_score || parameters.check_alignments) {
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
  align_benchmark_configure(&align_input);
  // Read-align loop
  int seqs_processed = 0, progress = 0;
  while (true) {
    // Read input sequence-pair
    const bool input_read = align_benchmark_read_input(
        input_file,&line1,&line2,&line1_allocated,
        &line2_allocated,seqs_processed,&align_input);
    if (!input_read) break;
    // Align queries using DP
    switch (parameters.algorithm) {
      case alignment_edit_dp:
        benchmark_edit_dp(&align_input);
        break;
      case alignment_edit_dp_banded:
        benchmark_edit_dp_banded(&align_input,parameters.bandwidth);
        break;
      case alignment_gap_lineal_nw:
        benchmark_gap_lineal_nw(&align_input,&parameters.lineal_penalties);
        break;
      case alignment_gap_affine_swg:
        benchmark_gap_affine_swg(&align_input,&parameters.affine_penalties);
        break;
      case alignment_gap_affine_swg_endsfree:
        benchmark_gap_affine_swg_endsfree(
            &align_input,&parameters.affine_penalties,
            parameters.pattern_begin_free,parameters.pattern_begin_free,
            parameters.text_begin_free,parameters.text_end_free);
        break;
      case alignment_gap_affine_swg_banded:
        benchmark_gap_affine_swg_banded(&align_input,
            &parameters.affine_penalties,parameters.bandwidth);
        break;
      case alignment_gap_affine_wavefront:
        benchmark_gap_affine_wavefront(&align_input,&parameters.affine_penalties);
        break;
      case alignment_gap_affine2p_dp:
        benchmark_gap_affine2p_dp(&align_input,&parameters.affine2p_penalties);
        break;
      case alignment_gap_affine2p_wavefront:
        benchmark_gap_affine2p_wavefront(&align_input,&parameters.affine2p_penalties);
        break;
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
      "USE: ./align_benchmark -a <algorithm> -i <input>                     \n"
      "      Options::                                                      \n"
      "        [Input]                                                      \n"
      "          --algorithm|a <algorithm>                                  \n"
      "            [edit]                                                   \n"
      "              edit-dp                                                \n"
      "              edit-dp-banded                                         \n"
      "              edit-wfa                                               \n"
      "              edit-wfa-adaptive                                      \n"
      "            [gap-lineal]                                             \n"
      "              gap-lineal-nw                                          \n"
      "              gap-lineal-wfa                                         \n"
      "              gap-lineal-wfa-adaptive                                \n"
      "            [gap-affine]                                             \n"
      "              gap-affine-swg                                         \n"
      "              gap-affine-swg-banded                                  \n"
      "              gap-affine-wfa                                         \n"
      "              gap-affine-wfa-adaptive                                \n"
      "            [gap-affine-2pieces]                                     \n"
      "              gap-affine2p-dp                                        \n"
      "              gap-affine2p-wfa                                       \n"
      "              gap-affine2p-wfa-adaptive                              \n"
      "          --input|i <File>                                           \n"
      "          --output|o <File>                                          \n"
      "        [Penalties]                                                  \n"
      "          --lineal-penalties|p M,X,I,D                               \n"
      "          --affine-penalties|g M,X,O,E                               \n"
      "          --affine2p-penalties M,X,O1,E1,O2,E2                       \n"
      "        [Wavefront parameters]                                       \n"
      "          --score-only                                               \n"
      "          --minimum-wavefront-length <INT>                           \n"
      "          --maximum-difference-distance <INT>                        \n"
      "          --low-memory                                               \n"
      "          --ends-free P0,Pf,T0,Tf                                    \n"
      "        [Misc]                                                       \n"
      "          --bandwidth <INT>                                          \n"
      "          --check|c 'correct'|'score'|'alignment'                    \n"
      "          --check-distance 'edit'|'gap-lineal'|'gap-affine'          \n"
      "          --check-bandwidth <INT>                                    \n"
      "          --plot                                                     \n"
      "        [System]                                                     \n"
      "          --max-memory <bytes>                                       \n"
      "          --progress|P <integer>                                     \n"
      "          --help|h                                                   \n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* Input */
    { "algorithm", required_argument, 0, 'a' },
    { "input", required_argument, 0, 'i' },
    /* Penalties */
    { "lineal-penalties", required_argument, 0, 'p' },
    { "affine-penalties", required_argument, 0, 'g' },
    { "affine2p-penalties", required_argument, 0, 900 },
    /* Wavefront parameters */
    { "score-only", no_argument, 0, 1001 },
    { "minimum-wavefront-length", required_argument, 0, 1002 },
    { "maximum-difference-distance", required_argument, 0, 1003 },
    { "low-memory", no_argument, 0, 1004 },
    { "ends-free", required_argument, 0, 1005 },
    /* Misc */
    { "bandwidth", required_argument, 0, 2000 },
    { "check", optional_argument, 0, 'c' },
    { "check-distance", required_argument, 0, 2001 },
    { "check-bandwidth", required_argument, 0, 2002 },
    { "plot", optional_argument, 0, 2003 },
    /* System */
    { "max-memory", required_argument, 0, 3000 },
    { "progress", required_argument, 0, 'P' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  if (argc <= 1) {
    usage();
    exit(0);
  }
  while (1) {
    c=getopt_long(argc,argv,"a:i:o:p:g:P:c:vh",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /*
     * Input
     */
    case 'a': {
      /* Test bench */
      if (strcmp(optarg,"test")==0) {
        parameters.algorithm = alignment_test;
      /* Edit */
      } else if (strcmp(optarg,"edit-dp")==0) {
        parameters.algorithm = alignment_edit_dp;
      } else if (strcmp(optarg,"edit-dp-banded")==0) {
        parameters.algorithm = alignment_edit_dp_banded;
      } else if (strcmp(optarg,"edit-wfa")==0) {
        parameters.reduction_type = wavefront_reduction_none;
        parameters.algorithm = alignment_edit_wavefront;
      } else if (strcmp(optarg,"edit-wfa-adaptive")==0) {
        parameters.reduction_type = wavefront_reduction_adaptive;
        parameters.algorithm = alignment_edit_wavefront;
      /* Gap-Lineal */
      } else if (strcmp(optarg,"gap-lineal-nw")==0 ||
                 strcmp(optarg,"gap-lineal-dp")==0) {
        parameters.algorithm = alignment_gap_lineal_nw;
      } else if (strcmp(optarg,"gap-lineal-wfa")==0) {
        parameters.reduction_type = wavefront_reduction_none;
        parameters.algorithm = alignment_gap_lineal_wavefront;
      } else if (strcmp(optarg,"gap-lineal-wfa-adaptive")==0) {
        parameters.reduction_type = wavefront_reduction_adaptive;
        parameters.algorithm = alignment_gap_lineal_wavefront;
      /* Gap-Affine */
      } else if (strcmp(optarg,"gap-affine-swg")==0 ||
                 strcmp(optarg,"gap-affine-dp")==0) {
        parameters.algorithm = alignment_gap_affine_swg;
      } else if (strcmp(optarg,"gap-affine-swg-banded")==0 ||
                 strcmp(optarg,"gap-affine-dp-banded")==0) {
        parameters.algorithm = alignment_gap_affine_swg_banded;
      } else if (strcmp(optarg,"gap-affine-wfa")==0) {
        parameters.reduction_type = wavefront_reduction_none;
        parameters.algorithm = alignment_gap_affine_wavefront;
      } else if (strcmp(optarg,"gap-affine-wfa-adaptive")==0) {
        parameters.reduction_type = wavefront_reduction_adaptive;
        parameters.algorithm = alignment_gap_affine_wavefront;
      /* Gap-Affine 2-Pieces */
      } else if (strcmp(optarg,"gap-affine2p-dp")==0) {
        parameters.algorithm = alignment_gap_affine2p_dp;
      } else if (strcmp(optarg,"gap-affine2p-wfa")==0) {
        parameters.reduction_type = wavefront_reduction_none;
        parameters.algorithm = alignment_gap_affine2p_wavefront;
      } else if (strcmp(optarg,"gap-affine2p-wfa-adaptive")==0) {
        parameters.reduction_type = wavefront_reduction_adaptive;
        parameters.algorithm = alignment_gap_affine2p_wavefront;
      } else {
        fprintf(stderr,"Algorithm '%s' not recognized\n",optarg);
        exit(1);
      }
      break;
    }
    case 'i':
      parameters.input_filename = optarg;
      break;
    case 'o':
      parameters.output_filename = optarg;
      break;
    /*
     * Penalties
     */
    case 'p': { // --lineal-penalties M,X,I,D
      char* sentinel = strtok(optarg,",");
      parameters.lineal_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.lineal_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.lineal_penalties.insertion = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.lineal_penalties.deletion = atoi(sentinel);
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
    /*
     * Wavefront parameters
     */
    case 1001: // --score-only
      parameters.score_only = true;
      break;
    case 1002: // --minimum-wavefront-length
      parameters.min_wavefront_length = atoi(optarg);
      break;
    case 1003: // --maximum-difference-distance
      parameters.max_distance_threshold = atoi(optarg);
      break;
    case 1004: // --low-memory
      parameters.low_memory = true;
      break;
    case 1005: { // --ends-free P0,Pf,T0,Tf
      parameters.endsfree = true;
      char* sentinel = strtok(optarg,",");
      parameters.pattern_begin_free = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.pattern_end_free = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.text_begin_free = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.text_end_free = atoi(sentinel);
      break;
    }
    /*
     * Misc
     */
    case 2000: // --bandwidth
      parameters.bandwidth = atoi(optarg);
      break;
    case 'c':
      if (optarg ==  NULL) { // default = score
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = false;
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
    case 2001: // --check-distance
      if (strcasecmp(optarg,"edit")==0) { // default = edit
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_EDIT;
      } else if (strcasecmp(optarg,"gap-lineal")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_LINEAL;
      } else if (strcasecmp(optarg,"gap-affine")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE;
      } else {
        fprintf(stderr,"Option '--check-distance' must be in {'edit','gap-lineal','gap-affine'}\n");
        exit(1);
      }
      break;
    case 2002: // --check-bandwidth
      parameters.check_bandwidth = atoi(optarg);
      break;
    case 2003: // --plot
      parameters.plot = (optarg==NULL) ? 1000 : atoi(optarg);
      break;
    /*
     * System
     */
    case 3000:
      parameters.max_memory = atol(optarg);
      break;
    case 'P':
      parameters.progress = atoi(optarg);
      break;
    case 'v':
      parameters.verbose = true;
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
  // Checks
  if (parameters.algorithm!=alignment_test && parameters.input_filename==NULL) {
    fprintf(stderr,"Option --input is required \n");
    exit(1);
  }
  if (parameters.endsfree) {
    switch (parameters.algorithm) {
      case alignment_gap_affine_swg:
        parameters.algorithm = alignment_gap_affine_swg_endsfree;
        break;
      case alignment_gap_lineal_wavefront:
      case alignment_gap_affine_wavefront:
      case alignment_gap_affine2p_wavefront:
        break;
      default:
        fprintf(stderr,"Ends-free variant not implemented selected algorithm\n");
        exit(1);
        break;
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
    align_benchmark();
  }
}
