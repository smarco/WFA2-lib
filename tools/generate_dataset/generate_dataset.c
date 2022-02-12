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
 * DESCRIPTION: Sequence Generator for benchmarking pairwise algorithms
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

/*
 * Common
 */
#define MIN(a,b) (((a)<=(b))?(a):(b))
#define MAX(a,b) (((a)>=(b))?(a):(b))
#define ABS(a) (((a)>=0)?(a):-(a))

/*
 * DNA Alphabet
 */
//#define ALPHABET_SIZE 26
//char alphabet[] = {
//    'A','B','C','D','E',
//    'F','G','H','I','J',
//    'K','L','M','N','O',
//    'P','Q','R','S','T',
//    'U','V','W','X','Y',
//    'Z'
//};
#define ALPHABET_SIZE 4
char alphabet[] = {
    'A','C','G','T'
};

/*
 * Sequence errors
 */
typedef struct {
  // Error log
  char operation;
  int position;
  // Mismatch detail
  char base_char;
  char replacement_char;
} seq_error_t;

/*
 * Random number generator [min,max)
 */
uint64_t rand_iid(
    const uint64_t min,
    const uint64_t max) {
  const int n_rand = rand(); // [0, RAND_MAX]
  const uint64_t range = max - min;
  const uint64_t rem = RAND_MAX % range;
  const uint64_t sample = RAND_MAX / range;
  // Consider the small interval within remainder of RAND_MAX
  if (n_rand < RAND_MAX - rem) {
    return min + n_rand/sample;
  } else {
    return rand_iid(min,max);
  }
}
/*
 * Generate random sequence
 */
void sequence_generate_random(
    char* const sequence,
    const int length) {
  // Generate random characters
  int i;
  for (i=0;i<length;++i) {
    sequence[i] = alphabet[rand_iid(0,ALPHABET_SIZE)];
  }
  sequence[length] = '\0';
}
/*
 * Extract sequence
 */
int sequence_extract(
    char* const reference,
    const int reference_length,
    char* const candidate,
    const int candidate_length) {
  const int length_diff = reference_length - candidate_length;
  const int offset = rand_iid(0,length_diff+1);
  strncpy(candidate,reference+offset,candidate_length);
  // Return offset
  return offset;
}
/*
 * Debug
 */
void sequence_errors_print(
    FILE* const stream,
    seq_error_t* const errors_log,
    const int num_errors) {
  int i;
  for (i=0;i<num_errors;++i) {
    switch (errors_log[i].operation) {
      case 'M':
        fprintf(stream,"(M,%d,%c->%c)",errors_log[i].position,
            errors_log[i].base_char,errors_log[i].replacement_char);
        break;
      default:
        fprintf(stream,"(%c,%d)",errors_log[i].operation,errors_log[i].position);
        break;
    }
  }
}
/*
 * Generate sequence errors
 */
void sequence_generate_mismatch(
    char* const sequence,
    const int sequence_length,
    seq_error_t* const errors_log) {
  // Generate random mismatch
  char character;
  int position;
  do {
    position = rand_iid(0,sequence_length);
    character = alphabet[rand_iid(0,ALPHABET_SIZE)];
  } while (sequence[position] == character);
  // Log
  errors_log->operation = 'M';
  errors_log->position = position;
  errors_log->base_char = sequence[position];
  errors_log->replacement_char = character;
  // Set mismatch
  sequence[position] = character;
}
void sequence_generate_deletion(
    char* const sequence,
    int* const sequence_length,
    seq_error_t* const errors_log) {
  // Generate random deletion
  const int position = rand_iid(0,*sequence_length);
  const int new_sequence_length = *sequence_length - 1;
  int i;
  for (i=position;i<new_sequence_length;++i) {
    sequence[i] = sequence[i+1];
  }
  *sequence_length = new_sequence_length;
  // Log
  errors_log->operation = 'D';
  errors_log->position = position;
}
void sequence_generate_insertion(
    char* const sequence,
    int* const sequence_length,
    seq_error_t* const errors_log) {
  // Generate random insertion
  const int position = rand_iid(0,*sequence_length);
  const int new_sequence_length = *sequence_length + 1;
  int i;
  for (i=new_sequence_length-1;i>position;--i) {
    sequence[i] = sequence[i-1];
  }
  *sequence_length = new_sequence_length;
  // Insert random character
  sequence[position] = alphabet[rand_iid(0,ALPHABET_SIZE)];
  // Log
  errors_log->operation = 'I';
  errors_log->position = position;
}
int sequence_generate_errors(
    char* const sequence,
    const int sequence_length,
    const int num_errors,
    seq_error_t* const errors_log) {
  int length = sequence_length;
  // Generate random errors
  int i;
  for (i=0;i<num_errors;++i) {
    int error_type = rand_iid(0,3);
    switch (error_type) {
      case 0: sequence_generate_mismatch(sequence,length,errors_log+i); break;
      case 1: sequence_generate_deletion(sequence,&length,errors_log+i); break;
      default: sequence_generate_insertion(sequence,&length,errors_log+i); break;
    }
  }
  // Close sequence and return length
  sequence[length] = '\0';
  return length;
}
/*
 * Parameters
 */
typedef struct {
  int num_reads;
  char *output;
  int pattern_length;
  int text_length;
  float error_degree;
  int debug;
} countour_bench_args;
countour_bench_args parameters = {
    .num_reads = 1000,
    .output = NULL,
    .pattern_length = 100,
    .text_length = 100,
    .error_degree = 0.04,
    .debug = 0,
};
/*
 * Menu & parsing cmd options
 */
void usage() {
  fprintf(stderr, "USE: ./generate-datasets [OPTIONS]...\n"
                  "      Options::\n"
                  "        --output|o         <File>    \n"
                  "        --num-patterns|n   <Integer> \n"
                  "        --length|l         <Integer> \n"
                  "        --pattern-length|P <Integer> \n"
                  "        --text-length|T    <Integer> \n"
                  "        --error|e          <Float>   \n"
                  "        --debug|g                    \n"
                  "        --help|h                     \n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "num-patterns", required_argument, 0, 'n' },
    { "output", required_argument, 0, 'o' },
    { "length", required_argument, 0, 'l' },
    { "pattern-length", required_argument, 0, 'P' },
    { "text-length", required_argument, 0, 'T' },
    { "error", required_argument, 0, 'e' },
    { "debug", no_argument, 0, 'g' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  if (argc <= 1) {
    usage();
    exit(0);
  }
  while (1) {
    c=getopt_long(argc,argv,"n:o:l:P:T:e:gh",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
      case 'n':
        parameters.num_reads = atoi(optarg);
        break;
      case 'o':
        parameters.output = optarg;
        break;
      case 'l': {
        const int length = atoi(optarg);
        parameters.pattern_length = length;
        parameters.text_length = length;
        break;
      }
      case 'P':
        parameters.pattern_length = atoi(optarg);
        break;
      case 'T':
        parameters.text_length = atoi(optarg);
        break;
      case 'e':
        parameters.error_degree = atof(optarg);
        break;
      case 'g':
        parameters.debug = 1;
        break;
      case 'h':
        usage();
        exit(1);
      case '?': default:
        fprintf(stderr, "Option not recognized \n"); exit(1);
    }
  }
}
int main(int argc,char* argv[]) {
  // Parsing command-line options
  parse_arguments(argc,argv);
  // Open file
  FILE* const output_file = (parameters.output==NULL) ? stdout : fopen(parameters.output,"w");
  // Allocate sequences
  const int reference_length = MAX(parameters.pattern_length,parameters.text_length); // Long sequence
  const int candidate_length = MIN(parameters.pattern_length,parameters.text_length); // Short sequence
  const int num_errors = (parameters.error_degree >= 1.0) ?
      (int) parameters.error_degree :
      ceil((float)candidate_length * parameters.error_degree);
  char* const reference = malloc(reference_length+1);
  char* const candidate = malloc(candidate_length+num_errors+1);
  seq_error_t* const errors_log = malloc(num_errors*sizeof(seq_error_t));
  // Read-align loop
  srand(time(0));
  int i;
  for (i=0;i<parameters.num_reads;++i) {
    // Generate random sequence
    sequence_generate_random(reference,reference_length);
    // Extract short sequence from long
    const int offset = sequence_extract(reference,reference_length,candidate,candidate_length);
    // Add errors
    sequence_generate_errors(candidate,candidate_length,num_errors,errors_log);
    // DEBUG
    if (parameters.debug) {
      fprintf(output_file,"#DEBUG offset=%d errors=",offset);
      sequence_errors_print(output_file,errors_log,num_errors);
      fprintf(output_file,"\n");
    }
    // Print
    if (parameters.pattern_length <= parameters.text_length) {
      fprintf(output_file,">%s\n",candidate);
      fprintf(output_file,"<%s\n",reference);
    } else {
      fprintf(output_file,"<%s\n",reference);
      fprintf(output_file,">%s\n",candidate);
    }
  }
  // Close files & free
  if (parameters.output!=NULL) fclose(output_file);
  free(reference);
  free(candidate);
  free(errors_log);
}
