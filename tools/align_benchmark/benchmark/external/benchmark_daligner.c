/*
 *  Wavefront Alignment Algorithms
 *  Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of Wavefront Alignment Algorithms.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Wavefront Alignment Algorithms Benchmark
 */

#include "benchmark/benchmark_utils.h"

// Myers Wavefront
#include "external/daligner_wave/align.h"
#include "external/daligner_wave/DB.h"
#include "external/daligner_wave/QV.h"

/*
 * Benchmark Indel Myers-wavefront (Improved O(ND) implemented by himself)
 */
int enlarge_vector(
    Work_Data *work,
    int newmax);
int enlarge_points(
    Work_Data *work,
    int newmax);
int forward_wave(
    Work_Data *work,
    Align_Spec *spec,
    Alignment *align,
    Path *bpath,
    int *mind,
    int maxd,
    int mida,
    int minp,
    int maxp);
void benchmark_daligner(
    align_input_t* const align_input) {
  // Parameters
  const double identity = 0.7;
  // Init structures
  Work_Data* w = New_Work_Data();
  enlarge_vector(w, 10000 * (6 * sizeof(int) + sizeof(uint64_t)));
  enlarge_points(w, 4 * 4096 / 100 * sizeof(uint16_t) + sizeof(Path));
  float freq[4] = { 0.25, 0.25, 0.25, 0.25 };
  Align_Spec *spec = New_Align_Spec(identity, 100, freq, 1);
  Path apath, bpath;
  apath.trace = malloc(4096);
  bpath.trace = malloc(4096);
  Alignment aln;
  aln.path = &apath;
  aln.flags = 0;
  aln.aseq = (char *)align_input->pattern;
  aln.bseq = (char *)align_input->text;
  aln.alen = align_input->pattern_length;
  aln.blen = align_input->text_length;
  int low = 0;
  int high = 0;
  // Align
  timer_start(&align_input->timer);
  forward_wave(w, spec, &aln, &bpath, &low, high, 0, -INT32_MAX, INT32_MAX);
  timer_stop(&align_input->timer);
  // DEBUG
  //  fprintf(stderr, "\n(%u, %u) :: (%u, %u) -> (%u, %u)\n",
  //      aln.alen, aln.blen,
  //      aln.path->abpos, aln.path->bbpos,
  //      aln.path->aepos, aln.path->bepos);
  //  Alignment_Cartoon(stdout,&aln,1,1);
  //  Print_Alignment(stdout,&aln,w,1,1,1,1,1);
  //  Print_Reference(stdout,&aln,w,1,1,1,1,1);
  //  int INDENT    = 4;
  //  int WIDTH     = 100;
  //  int BORDER    = 10;
  //  Alignment_Cartoon(stdout,&aln,INDENT,5);
  //  //Print_Reference(stdout,&aln,w,INDENT,WIDTH,BORDER,0,5);
  //  Print_Alignment(stdout,&aln,w,INDENT,WIDTH,BORDER,0,5);
  // Free
  free(apath.trace);
  free(bpath.trace);
  Free_Align_Spec(spec);
  Free_Work_Data(w);
}
