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
 * DESCRIPTION: WaveFront-Alignment module for plot
 */

#include "wavefront_plot.h"
#include "wavefront_aligner.h"

/*
 * Setup
 */
void wavefront_plot_allocate(
    wavefront_plot_t* const wf_plot,
    const distance_metric_t distance_metric,
    const int pattern_length,
    const int text_length,
    wavefront_plot_params_t* const plot_params) {
  // Compute dimensions
  const int min_v = (plot_params->min_v == -1) ? 0 : plot_params->min_v;
  const int max_v = (plot_params->max_v == -1) ? pattern_length-1 : plot_params->max_v;
  const int min_h = (plot_params->min_h == -1) ? 0 : plot_params->min_h;
  const int max_h = (plot_params->max_h == -1) ? text_length-1 : plot_params->max_h;
  // Wavefront Components
  wf_plot->m_heatmap = heatmap_new(heatmap_min,
      min_v,max_v,min_h,max_h,plot_params->resolution_points);
  if (distance_metric == gap_affine) {
    wf_plot->i1_heatmap = heatmap_new(heatmap_min,
        min_v,max_v,min_h,max_h,plot_params->resolution_points);
    wf_plot->d1_heatmap = heatmap_new(heatmap_min,
        min_v,max_v,min_h,max_h,plot_params->resolution_points);
  } else {
    wf_plot->i1_heatmap = NULL;
    wf_plot->d1_heatmap = NULL;
  }
  if (distance_metric == gap_affine_2p) {
    wf_plot->i2_heatmap = heatmap_new(heatmap_min,
        min_v,max_v,min_h,max_h,plot_params->resolution_points);
    wf_plot->d2_heatmap = heatmap_new(heatmap_min,
        min_v,max_v,min_h,max_h,plot_params->resolution_points);
  } else {
    wf_plot->i2_heatmap = NULL;
    wf_plot->d2_heatmap = NULL;
  }
  // Behavior
  wf_plot->behavior_heatmap = heatmap_new(heatmap_value,
      min_v,max_v,min_h,max_h,plot_params->resolution_points);
}
void wavefront_plot_free(
    wavefront_plot_t* const wf_plot) {
  heatmap_delete(wf_plot->m_heatmap);
  if (wf_plot->i1_heatmap) heatmap_delete(wf_plot->i1_heatmap);
  if (wf_plot->d1_heatmap) heatmap_delete(wf_plot->d1_heatmap);
  if (wf_plot->i2_heatmap) heatmap_delete(wf_plot->i2_heatmap);
  if (wf_plot->d2_heatmap) heatmap_delete(wf_plot->d2_heatmap);
  heatmap_delete(wf_plot->behavior_heatmap);
}
/*
 * Accessors
 */
void wavefront_plot_component(
    wavefront_t* const wavefront,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score,
    heatmap_t* const wf_heatmap,
    heatmap_t* const extend_heatmap) {
  if (wavefront != NULL) {
    int k;
    for (k=wavefront->lo;k<=wavefront->hi;++k) {
      const wf_offset_t offset = wavefront->offsets[k];
      if (offset >= 0) {
        // Compute coordinates
        int v = WAVEFRONT_V(k,offset);
        int h = WAVEFRONT_H(k,offset);
        if (v>=pattern_length || h>=text_length) continue;
        heatmap_set(wf_heatmap,v,h,score);
        // Simulate extension
        if (extend_heatmap != NULL) {
          while (v<pattern_length && h<text_length && pattern[v++]==text[h++]) {
            heatmap_set(wf_heatmap,v,h,score);
            heatmap_set(extend_heatmap,v,h,10);
          }
        }
      }
    }
  }
}
void wavefront_plot(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const char* const text,
    const int score) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const int s = (wf_components->memory_modular) ? score%wf_components->max_score_scope : score;
  // Plot wavefront components
  wavefront_plot_component(
      wf_components->mwavefronts[s],
      pattern,pattern_length,text,text_length,
      score,wf_aligner->wf_plot.m_heatmap,
      wf_aligner->wf_plot.behavior_heatmap);
  if (distance_metric == gap_affine) {
    wavefront_plot_component(
        wf_components->i1wavefronts[s],
        pattern,pattern_length,text,text_length,
        score,wf_aligner->wf_plot.i1_heatmap,NULL);
    wavefront_plot_component(
        wf_components->d1wavefronts[s],
        pattern,pattern_length,text,text_length,
        score,wf_aligner->wf_plot.d1_heatmap,NULL);
  }
  if (distance_metric == gap_affine_2p) {
    wavefront_plot_component(
        wf_components->i2wavefronts[s],
        pattern,pattern_length,text,text_length,
        score,wf_aligner->wf_plot.i2_heatmap,NULL);
    wavefront_plot_component(
        wf_components->d2wavefronts[s],
        pattern,pattern_length,text,text_length,
        score,wf_aligner->wf_plot.d2_heatmap,NULL);
  }
}
//void wavefront_plot_cutoff(
//    wavefront_aligner_t* const wf_aligner,
//    const int score,
//    const int lo_base,
//    const int lo_reduced,
//    const int hi_base,
//    const int hi_reduced) {
//  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
//  const int s = (wf_components->memory_modular) ? score%wf_components->max_score_scope : score;
//  wavefront_t* const wavefront = wf_components->mwavefronts[s];
//  heatmap_t* const heatmap = wf_aligner->wf_plot.behavior_heatmap;
//  int k;
//  for (k=lo_base;k<lo_reduced;++k) {
//    const wf_offset_t offset = wavefront->offsets[k];
//    if (offset >= 0) heatmap_set(heatmap,WAVEFRONT_V(k,offset),WAVEFRONT_H(k,offset),20);
//  }
//  for (k=hi_reduced+1;k<=hi_base;++k) {
//    const wf_offset_t offset = wavefront->offsets[k];
//    if (offset >= 0) heatmap_set(heatmap,WAVEFRONT_V(k,offset),WAVEFRONT_H(k,offset),20);
//  }
//}
/*
 * Display
 */
void wavefront_plot_print_cigar(
    FILE* const stream,
    cigar_t* const cigar,
    const char target_operation) {
  int i, h=0, v=0, count=0;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    // Print point
    const char operation = cigar->operations[i];
    if (operation == target_operation) {
      if (count++ > 0) fprintf(stream,";");
      fprintf(stream,"%d,%d",h,v);
    }
    // Check operation
    switch (operation) {
      case 'M': case 'X': ++h; ++v; break;
      case 'I': ++h; break;
      case 'D': ++v; break;
      default: break;
    }
  }
}
void wavefront_plot_print(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_plot_t* const wf_plot = &wf_aligner->wf_plot;
  // Metadata
  fprintf(stream,"# PatternLength %d\n",wf_aligner->pattern_length);
  fprintf(stream,"# TextLength %d\n",wf_aligner->text_length);
  fprintf(stream,"# Penalties ");
  wavefronts_penalties_print(stream,&wf_aligner->penalties);
  fprintf(stream,"\n");
  // Alignment mode
  fprintf(stream,"# WFAMode (");
  fprintf(stream,"%s",(wf_aligner->alignment_scope==compute_score)?"S":"A");
  fprintf(stream,"%c",(wf_aligner->wf_components.bt_piggyback)?'L':'F');
  fprintf(stream,"%c",(wf_aligner->alignment_form.span==alignment_end2end)?'G':'S');
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  if (wf_heuristic->strategy != wf_heuristic_none) {
    wavefront_heuristic_print(stream,wf_heuristic);
  }
  fprintf(stream,")\n");
  // Wavefront components
  fprintf(stream,"# Heatmap M\n"); heatmap_print(stream,wf_plot->m_heatmap);
  if (distance_metric == gap_affine) {
    fprintf(stream,"# Heatmap I1\n"); heatmap_print(stream,wf_plot->i1_heatmap);
    fprintf(stream,"# Heatmap D1\n"); heatmap_print(stream,wf_plot->d1_heatmap);
  }
  if (distance_metric == gap_affine_2p) {
    fprintf(stream,"# Heatmap I2\n"); heatmap_print(stream,wf_plot->i2_heatmap);
    fprintf(stream,"# Heatmap D2\n"); heatmap_print(stream,wf_plot->d2_heatmap);
  }
  // Extend
  fprintf(stream,"# Heatmap Extend\n"); heatmap_print(stream,wf_plot->behavior_heatmap);
  // CIGAR
  if (wf_aligner->alignment_scope == compute_alignment) {
    fprintf(stream,"# List CIGAR-M ");
    wavefront_plot_print_cigar(stream,&wf_aligner->cigar,'M');
    fprintf(stream,"\n");
    fprintf(stream,"# List CIGAR-X ");
    wavefront_plot_print_cigar(stream,&wf_aligner->cigar,'X');
    fprintf(stream,"\n");
    fprintf(stream,"# List CIGAR-I ");
    wavefront_plot_print_cigar(stream,&wf_aligner->cigar,'I');
    fprintf(stream,"\n");
    fprintf(stream,"# List CIGAR-D ");
    wavefront_plot_print_cigar(stream,&wf_aligner->cigar,'D');
    fprintf(stream,"\n");
  }
}
