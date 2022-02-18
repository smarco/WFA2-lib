# WFA

## 1. WFA TOOLS AND EXAMPLES

### 1.1 What is this?

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

## 2 GENERATE DATASET TOOL

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

## 3. ALIGNMENT BENCHMARK TOOL

### 3.1 Introduction to benchmarking WFA. Simple tests

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

### 3.2 Command-line options

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


## 4. AUTHORS

  Santiago Marco-Sola \- santiagomsola@gmail.com     

## 5. REPORTING BUGS

Feedback and bug reporting it's highly appreciated. 
Please report any issue or suggestion on github, or by email to the main developer (santiagomsola@gmail.com).

## 6. LICENSE

WFA is distributed under MIT licence.

## 7. CITATION

**Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa**. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.
