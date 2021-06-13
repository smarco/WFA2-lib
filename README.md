# WFA

## 1. INTRODUCTION

### 1.1 What is WFA?

The wavefront alignment (WFA) algorithm is an exact gap-affine algorithm that takes advantage of  
homologous regions between the sequences to accelerate the alignment process. As opposed to 
traditional dynamic programming algorithms that run in quadratic time, the WFA runs in time O(ns),
proportional to the read length n and the alignment score s, using O(s^2) memory. Moreover, the WFA
exhibits simple data dependencies that can be easily vectorized, even by the automatic features of 
modern compilers, for different architectures, without the need to adapt the code.

This library implements the WFA and the WFA-Adapt algorithms for gap-affine penalties. It also 
provides support functions to display and verify the results. Moreover, it implements a benchmarking
tool that serves to evaluate the performance of these two algorithms, together with other 
high-performance alignment methods (checkout branch `benchmark`). The library can be executed   
through the benchmarking tool for evaluation purposes or can be integrated into your code by calling
the WFA functions.

If you are interested in benchmarking WFA with other algorithms implemented or integrated into the
WFA library, checkout branch `benchmark`.

### 1.2 Getting started

Note: We recomend using the GCC compiler

Clone GIT and compile

```
$> git clone https://github.com/smarco/WFA.git WFA
$> cd WFA
$> make clean all
```

## 3. PROGRAMMING WITH WFA

Inside the folder `tools/examples/` the user can find two simples examples of how to program using 
the WFA library. These examples aim to illustrate how to integrate the WFA code into any tool.

### 3.1 Simple WFA example

This simple example illustrates how to align two sequences using the gap-affine WFA algorithm.
First, we need to include the WFA alignment module.

```C
#include "gap_affine/affine_wavefront_align.h"
```

Then, we prepare the text, pattern, penalties, and the memory-managed (MM) allocator. Note that the
`affine_penalties` is configured in terms of penalties. For that reason, mismatch, gap-opening, and
gap-extension are supposed to be positive values.  

```C
  // Patter & Text
  char* pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
  char* text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
  // Allocate MM
  mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  // Set penalties
  affine_penalties_t affine_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening = 6,
      .gap_extension = 2,
  };
```

Afterwards, we initialize the `affine_wavefronts` object and we align the pattern against the text
using the configured penalties.

```C
  // Init Affine-WFA
  affine_wavefronts_t* affine_wavefronts = affine_wavefronts_new_complete(
      strlen(pattern),strlen(text),&affine_penalties,NULL,mm_allocator);
  // Align
  affine_wavefronts_align(
      affine_wavefronts,pattern,strlen(pattern),text,strlen(text));
```

Finally, we can display the results of the alignment process. For example, the alignment score and
the alignment CIGAR. For this purpose, the function `edit_cigar_score_gap_affine` computes the 
CIGAR score, and the function `edit_cigar_print_pretty` prints pretty the CIGAR.

```C
  // Display alignment
  const int score = edit_cigar_score_gap_affine(
      &affine_wavefronts->edit_cigar,&affine_penalties);
  fprintf(stderr,"  PATTERN  %s\n",pattern);
  fprintf(stderr,"  TEXT     %s\n",text);
  fprintf(stderr,"  SCORE COMPUTED %d\t",score);
  edit_cigar_print_pretty(stderr,
      pattern,strlen(pattern),text,strlen(text),
      &affine_wavefronts->edit_cigar,mm_allocator);
  // Free
  affine_wavefronts_delete(affine_wavefronts);
  mm_allocator_delete(mm_allocator);
```

Compile and run:

```
$> gcc -O3 -I../.. -L../../build wfa_basic.c -o wfa_basic -lwfa
$> ./wfa_basic
```

### 3.2 WFA-Adaptive example

This example shows how to use the adaptive version of the WFA (i.e., WFA-Adaptive) to further improve
the performance of the WFA algorithm by discarding alignment paths that are unlikely to reach the 
optimal solution. This example is very similar to the previous one, we only have to include the 
parameters `minimum-wavefront-length` and `maximum-difference-distance`. 
 
 
```C
  const int min_wavefront_length = 10;
  const int max_distance_threshold = 50;
  // Init Affine-WFA
  affine_wavefronts_t* affine_wavefronts = affine_wavefronts_new_reduced(
      strlen(pattern),strlen(text),&affine_penalties,
      min_wavefront_length,max_distance_threshold,NULL,mm_allocator);
  // Align
  affine_wavefronts_align(
      affine_wavefronts,pattern,strlen(pattern),text,strlen(text));

```

In this example, we show how to access the individual elements of the CIGAR 
(i.e., 'M','X','I', and 'D') encoded using plain 8-bit ASCII.


```C
  // Count mismatches, deletions, and insertions
  int i, misms=0, ins=0, del=0;
  edit_cigar_t* const edit_cigar = &affine_wavefronts->edit_cigar;
  for (i=edit_cigar->begin_offset;i<edit_cigar->end_offset;++i) {
    switch (edit_cigar->operations[i]) {
      case 'M': break;
      case 'X': ++misms; break;
      case 'D': ++del; break;
      case 'I': ++ins; break;
    }
  }
  fprintf(stderr,
      "Alignment contains %d mismatches, %d insertions, "
      "and %d deletions\n",misms,ins,del);
```

Compile and run:

```
$> gcc -O3 -I../.. wfa_adapt.c ../../build/libwfa.a -o wfa_adapt
$> ./wfa_basic
```

### 3.3 Aligning sequences longer than 65.536 bases

By default, the WFA uses 16-bit integers to represent the alignment wavefronts. For that reason,
the maximum sequence length allowed is 2^16. In case you want to align longer sequences, you
must adjust the definitions on `gap_affine/affine_wavefront.h` and select `AFFINE_WAVEFRONT_W32`.

```
/*
 * Offset size
 */
//#define AFFINE_WAVEFRONT_W8
//#define AFFINE_WAVEFRONT_W16
#define AFFINE_WAVEFRONT_W32
```

## 4. BENCHMARKING. COMMAND-LINE AND OPTIONS

### 4.1 Introduction to benchmarking WFA. Simple tests

The WFA includes the benchmarking tool *align-benchmark* to test and compare the performance of 
several pairwise alignment implementations, including the WFA and WFA-Adapt. This tool takes as 
input a dataset containing pairs of sequences (i.e., pattern and text) to align. Patterns are 
preceded by the '>' symbol and texts by the '<' symbol. Example:

```
>ATTGGAAAATAGGATTGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTCGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTAGCTCGAAGCCCA
<GATTGGAAAATAGGATGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTGCTCGAAGCCCA
>CCGTAGAGTTAGACACTCGACCGTGGTGAATCCGCGACCACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCAGTGATTAAAC
<CCTAGAGTTAGACACTCGACCGTGGTGAATCCGCGATCTACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCGAGTGATTAAAC
[...]
```

You can either generate a custom dataset of your own, or use the *generate-dataset* tool to generate
a random dataset. For example, the following command generates a dataset named 'sample.dataset.seq' 
of 5M pairs of 100 bases with an alignment error of 5% (i.e., 5 mismatches, insertions or deletions 
per alignment).

```
$> ./bin/generate_dataset -n 5000000 -l 100 -e 0.05 -o sample.dataset.seq
```

Once you have the dataset ready, you can run the *align-benchmark* tool to benchmark the performance 
of a specific pairwise alignment method. For example, the WFA algorithm:

```
$> ./bin/align_benchmark -i sample.dataset.seq -a gap-affine-wfa
...processed 10000 reads (benchmark=125804.398 reads/s;alignment=188049.469 reads/s)
...processed 20000 reads (benchmark=117722.406 reads/s;alignment=180925.031 reads/s)
[...]
...processed 5000000 reads (benchmark=113844.039 reads/s;alignment=177325.281 reads/s)
[Benchmark]
=> Total.reads            5000000
=> Time.Benchmark        43.92 s  (    1   call,  43.92  s/call {min43.92s,Max43.92s})
  => Time.Alignment      28.20 s  ( 64.20 %) (    5 Mcalls,   5.64 us/call {min438ns,Max47.05ms})
```

The *align-benchmark* tool will finish and report overall benchmark time (including reading the 
input, setup, checking, etc.) and the time taken by the algorithm (i.e., *Time.Alignment*). If you 
want to measure the accuracy of the alignment method, you can add the option `--check` and all the
alignments will be verified. 

```
$> ./bin/align_benchmark -i sample.dataset.seq -a gap-affine-wfa --check
...processed 10000 reads (benchmark=14596.232 reads/s;alignment=201373.984 reads/s)
...processed 20000 reads (benchmark=13807.268 reads/s;alignment=194224.922 reads/s)
[...]
...processed 5000000 reads (benchmark=10625.568 reads/s;alignment=131371.703 reads/s)
[Benchmark]
=> Total.reads            5000000
=> Time.Benchmark         7.84 m  (    1   call, 470.56  s/call {min470.56s,Max470.56s})
  => Time.Alignment      28.06 s  (  5.9 %) (    5 Mcalls,   5.61 us/call {min424ns,Max73.61ms})
[Accuracy]
 => Alignments.Correct        5.00 Malg        (100.00 %) (samples=5M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
 => Score.Correct             5.00 Malg        (100.00 %) (samples=5M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => Score.Total           147.01 Mscore uds.            (samples=5M{mean29.40,min0.00,Max40.00,Var37.00,StdDev6.00)}
     => Score.Diff            0.00 score uds.  (  0.00 %) (samples=0,--n/a--)}
 => CIGAR.Correct             0.00 alg         (  0.00 %) (samples=0,--n/a--)}
   => CIGAR.Matches         484.76 Mbases      ( 96.95 %) (samples=484M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => CIGAR.Mismatches        7.77 Mbases      (  1.55 %) (samples=7M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => CIGAR.Insertions        7.47 Mbases      (  1.49 %) (samples=7M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => CIGAR.Deletions         7.47 Mbases      (  1.49 %) (samples=7M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}

```

Using the `--check` option, the tool will report *Alignments.Correct* (i.e., total alignments that 
are correct, not necessarily optimal), and *Score.Correct* (i.e., total alignments that have the 
optimal score). Note that the overall benchmark time will increase due to the overhead introduced  
by the checking routine, however the *Time.Alignment* should remain the same.

### 4.2 Generate-dataset tool (Command-line and Options)

```
        --output|o        <File>
          Filename/Path to the output dataset.
          
        --num-patterns|n  <Integer>
          Total number of pairs pattern-text to generate.
          
        --length|l        <Integer>
          Total length of the pattern.
          
        --error|e         <Float>
          Total error-rate between the pattern and the text (allowing single-base mismatches, 
          insertions and deletions). This parameter may modify the final length of the text.
          
        --help|h
          Outputs a succinct manual for the tool.
```

### 4.3 Align-benchmark tool (Command-line and Options)

Summary of algorithms/methods implemented within the benchmarking tool. If you are interested 
in benchmarking WFA with other algorithms implemented or integrated into the WFA library, 
checkout branch `benchmark`.

|      Algorithm Name        |       Code-name       | Distance Model |  Output   | Implementation | Extra Parameters                                           |
|----------------------------|:---------------------:|:--------------:|:---------:|----------------|------------------------------------------------------------|
|DP Edit                     |edit-dp                |  Edit-distace  | Alignment |WFA             |                                                            |
|DP Edit Banded              |edit-dp-banded         |  Edit-distace  | Alignment |WFA             | --bandwidth                                                |
|DP Gap-lineal               |gap-lineal-nw          |   Gap-lineal   | Alignment |WFA             |                                                            |
|DP Gap-affine               |gap-affine-swg         |   Gap-affine   | Alignment |WFA             |                                                            |
|DP Gap-affine Banded        |gap-affine-swg-banded  |   Gap-affine   | Alignment |WFA             | --bandwidth                                                |
|WFA Gap-affine              |gap-affine-wfa         |   Gap-affine   | Alignment |WFA             |                                                            |
|WFA Gap-affine Adaptive     |gap-affine-wfa-adaptive|   Gap-affine   | Alignment |WFA             | --minimum-wavefront-length / --maximum-difference-distance |

#### - Input

```
          --algorithm|a <algorithm-code-name> 
            Selects pair-wise alignment algorithm/implementation.
                                                       
          --input|i <File>
            Filename/path to the input SEQ file. That is, file containing the sequence pairs to
            align. Sequences are stored one per line, grouped by pairs where the pattern is 
            preceded by '>' and text by '<'.
```
                                     
#### - Penalties

```                                                  
          --lineal-penalties|p M,X,I,D
            Selects gap-lineal penalties for those alignment algorithms that use this penalty model.
            Example: --lineal-penalties="-1,1,2,2"
                
          --affine-penalties|g M,X,O,E
            Selects gap-affine penalties for those alignment algorithms that use this penalty model.
            Example: --affine-penalties="-1,4,2,6" 
          
```
                         
#### - Specifics

```                                                  
          --bandwidth <INT>
            Selects the bandwidth size for those algorithms that use bandwidth strategy. 
                
          --minimum-wavefront-length <INT>
            Selects the minimum wavefront length to trigger the WFA-Adapt reduction method.
            
          --maximum-difference-distance <INT>
            Selects the maximum difference distance for the WFA-Adapt reduction method.  
```
                   
#### - Misc

```                                                       
          --progress|P <integer>
            Set the progress message periodicity.
            
          --check|c 'correct'|'score'|'alignment'                    
            Activates the verification of the alignment results. 
          
          --check-distance 'edit'|'gap-lineal'|'gap-affine'
            Select the alignment-model to use for verification of the results.
          
          --check-bandwidth <INT>
            Sets a bandwidth for the simple verification functions.

          --help|h
            Outputs a succinct manual for the tool.
```


## 5. AUTHORS

  Santiago Marco-Sola \- santiagomsola@gmail.com     

## 6. REPORTING BUGS

Feedback and bug reporting it's highly appreciated. 
Please report any issue or suggestion on github, or by email to the main developer (santiagomsola@gmail.com).

## 7. LICENSE

WFA is distributed under MIT licence.

## 8. CITATION

**Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa**. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.
