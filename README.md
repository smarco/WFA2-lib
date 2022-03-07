# WFA

## 1. INTRODUCTION

### 1.1 What is WFA?

The wavefront alignment (WFA) algorithm is an **exact** gap-affine algorithm that takes advantage of homologous regions between the sequences to accelerate the alignment process. As opposed to traditional dynamic programming algorithms that run in quadratic time, the WFA runs in time `O(ns+s^2)`, proportional to the sequence length `n` and the alignment score `s`, using `O(s^2)` memory. Moreover, the WFA algorithm exhibits simple computational patterns that can be automatically vectorized by the modern compilers for different architectures without the need to adapt the code. To intuitively illustrate why the WFA algorithm is so interesting, have a look at the following figure. On the left, it shows the cells computed by a classical dynamic programming based algorithm (like Smith-Waterman or Needleman Wunsch). On the right, the cells computed by the WFA algorithm to obtain the same result (i.e., the optimal alignment).

<p align = "center">
<img src = "img/00.wfa.vs.swg.small.png" width = "750px">
</p>

### 1.2 What is WFA2-lib?

The WFA2 library implements the WFA algorithm for different distance metrics and alignment modes. It supports a wide range of [distance functions](#wfa2.distances): indel, edit, gap-lineal, gap-affine, and dual-gap gap-affine distances. The library allows computing only the score or the complete alignment (i.e., CIGAR) (see [Alignment Scope](#wfa2.scope)). Also, the WFA2 library supports computing end-to-end alignments (a.k.a. global-alignment) and ends-free alignments (including semi-global, glocal, and extension alignment) (see [Alignment Span](#wfa2.span)). In the case of long and noise alignments, the library provides different [low-memory modes](#wfa2.mem) that significantly reduce the memory usage of the naive WFA algorithm implementation. Beyond the exact-alignment modes, the WFA2 library implements [heuristic modes](#wfa2.heuristics) that dramatically accelerate the alignment computation. Additionally, the library provides many other support functions to display and verify alignment results, control the overall memory usage, and more.

### 1.3 Getting started

Clone GIT and compile the library, tools, and examples.

```
$> git clone https://github.com/smarco/WFA2-lib
$> cd WFA2-lib
$> make clean all
```

### 1.4 Contents (where to go from here)

On the section [WFA2-lib features](#wfa2.features), you can explore the most relevant options and features of the library. Then, inside the folder [tools/](tools/README.md), you can find tools that can be used to execute and understand the WFA2 library capabilities. Additionally, inside [examples/](examples/README.md) you can find simple examples illustrating how to integrate the WFA2 code into any tool.

* [WFA2-lib Features](#wfa2.features)
    * [Distance Metrics](#wfa2.distances)
    * [Alignment Scope](#wfa2.scope)
    * [Alignment Span](#wfa2.span)
    * [Memory modes](#wfa2.mem)
    * [Heuristic modes](#wfa2.heuristics)
    * [Other options](#wfa2.other)
* [Programming with WFA2-lib](#wfa2.programming)
    * [Simple C example](#wfa2.programming.c)
    * [Simple C++ example](#wfa2.programming.cpp)
* [Reporting Bugs and Feature Request](#wfa2.complains)
* [License](#wfa2.licence)
* [Citation](#wfa2.cite) 

### 1.5 Important notes and clarifications

- The WFA algorithm is an **exact algorithm**. If no heuristic is applied (e.g., band or adaptive pruning), the core algorithm guarantees to always find the optimal solution (i.e., best alignment score). Since its first release, some authors have referenced the WFA as approximated or heuristic, and that is NOT the case.


- Given two sequences of length `n`, traditional dynamic-programming (DP) based methods (like Smith-Waterman or Needleman-Wunsch) compute the optimal alignment in `O(n^2)` time, using `O(n^2)` memory. In contrast, the WFA algorithm requires `O(ns+s^2)` and `O(s^2)` memory (being `s` the optimal alignment score). Therefore, **the memory consumption of the WFA algorithm is not intrinsically higher than that of other methods**. Most DP-based methods can use heuristics (like banded, X-drop, or Z-drop) to reduce the execution time and the memory usage at the expense of losing accuracy. Likewise, **the WFA algorithm can use heuristics too in order to dramatically reduce the execution time and the memory usage**.


- The WFA algorithm **works with plain ASCII strings**. Although we mainly focus on aligning DNA/RNA sequences, the algorithm and the WFA2-lib implementation work with any pair of strings. Moreover, these sequences doesn't have to be pre-processed (like with the Bit-Parallel Myers algorithm, used within EdLib), nor any table has to be precomputed (like the query profile, used within some Smith-Waterman implementations).


- Due to its simplicity, the WFA algorithm can be automatically vectorized for any SIMD-compliant CPU supported by the compiler. For this reason, **the WFA2-lib implementation is independent of any specific ISA or processor model**. As opposed to other hardware-dependent libraries, we aim to offer a multiplatform pairwise-alignment library that can be executed on different processors and models (e.g., SSE, AVX2, AVX512, POWER-ISA, ARM, NEON, SVE, SVE2, RISCV-RVV, ...).


- **A note for the fierce competitors.** I can understand that science and publishing have become a fierce competition these days. Many researchers want their methods to be successful and popular, seeking funding, tenure, or even fame. If you are going to benchmark the WFA using the least favourable configuration, careless programming, and a disadvantageous setup, please, go ahead. But remember, researchers like you have put a lot of effort into developing the WFA. We all joined this "competition" because we sought to find better methods that could truly help other researchers. So, try to be nice, tone down the marketing, and produce fair evaluations and honest publications.

## <a name="wfa2.features"></a> 2. WFA2-LIB FEATURES

* Support for most common and widely-used distance metrics (i.e., indel, edit, gap-lineal, gap-affine, and dual-gap gap-affine).
* Support for end-to-end (a.k.a. global) and ends-free alignment modes.
* Low memory modes to regulate the memory consumption.
* Heuristic methods to use on top of the core WFA algorithm

### <a name="wfa2.distances"></a> 2.1 Distance Metrics

The WFA2 library implements the wavefront algorithm for the most widely used distance metrics. Depending on the distance function, the practical alignment time can change, although the computational complexity always remains proportional to the alignment score or distance. The WFA2 library offers the following distance metrics or functions:

- **Indel (or LCS).** Produces alignments allowing matches, insertions, and deletions with unitary cost (i.e., {M,I,D} = {0,1,1}) but not mismatches. Also known as the longest common subsequence (LCS) problem. The LCS is defined as the longest subsequence common to both sequences, provided that the characters of the subsequence are not required to occupy consecutive positions within the original sequences. 

```
    PATTERN    A-GCTA-GTGTC--AATGGCTACT-T-T-TCAGGTCCT
               |  ||| |||||    |||||||| | | |||||||||
    TEXT       AA-CTAAGTGTCGG--TGGCTACTATATATCAGGTCCT
    ALIGNMENT  1M1I1D3M1I5M2I2D8M1I1M1I1M1I9M
```

<details><summary>Configuration Indel Alignment</summary>
<p>

```C
    // Configuration
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = indel;
```

</p>
</details>

- **Edit (a.k.a. Levenshtein).** Produces alignments allowing matches, mismatches, insertions, and deletions with unitary cost (i.e., {M,X,I,D} = {0,1,1,1}). Edit or Levenshtein distance between two sequences is the minimum number of single-character edits (i.e., insertions, deletions, or mismatches) required to transform one sequence into the other.

```
    PATTERN    AGCTA-GTGTCAATGGCTACT-T-T-TCAGGTCCT
               | ||| |||||  |||||||| | | |||||||||
    TEXT       AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT
    ALIGNMENT  1M1X3M1I5M2X8M1I1M1I1M1I9M
```

<details><summary>Configuration Edit Alignment</summary>
<p>

```C
    // Configuration
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = edit;
```

</p>
</details>

- **Gap-linear (as in Needleman-Wunsch).** Produces alignments allowing matches, mismatches, insertions, and deletions. Allows assigning a penalty (a.k.a. cost or weight) to each alignment operation. It computes the optimal alignment, minimizing the overall cost to transform one sequence into the other. Under the gap-linear model, the alignment score is computed based on {X,I}⁠, where X corresponds to the mismatch penalty and the gap penalty is expressed as the function l(N)=N·I (given the length of the gap N and the gap penalty I).

```
    PATTERN    A-GCTA-GTGTC--AATGGCTACT-T-T-TCAGGTCCT
               |  ||| |||||    |||||||| | | |||||||||
    TEXT       AA-CTAAGTGTCGG--TGGCTACTATATATCAGGTCCT
    ALIGNMENT  1M1I1D3M1I5M2I2D8M1I1M1I1M1I9M
```

<details><summary>Configuration Gap-Linear Alignment</summary>
<p>

```C
    // Configuration
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_linear;
    attributes.linear_penalties.mismatch = 6; // X > 0 
    attributes.linear_penalties.indel = 2;    // I > 0
```

</p>
</details>

- **Gap-affine (as in Smith-Waterman-Gotoh).**

Linear gap cost functions can lead to alignments populated with small gaps. In certain scenarios, like genomics or evolutionary studies, long gaps are preferred (understood as a single event). Under the gap-affine model, the alignment score is computed based on {X,O,E}⁠, where X corresponds to the mismatch penalty and the gap penalty is expressed as the function g(N)=O+N·E (given the length of the gap N, the gap opening penalty O, and the gap extension penalty E).

```
    PATTERN    AGCTA-GTGTCAATGGCTACT---TTTCAGGTCCT
               | ||| |||||  ||||||||   | |||||||||
    TEXT       AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT
    ALIGNMENT  1M1X3M1I5M2X8M3I1M1X9M
```

<details><summary>Configuration Gap-Affine Alignment</summary>
<p>

```C
    // Configuration
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
    attributes.affine_penalties.mismatch = 6;      // X > 0
    attributes.affine_penalties.gap_opening = 4;   // O >= 0
    attributes.affine_penalties.gap_extension = 2; // E > 0
```

</p>
</details>

- **Dual-cost gap-affine distances.** Also known as piece-wise gap-affine cost, this distance metric addresses some issues that the regular gap-affine distance has with long gaps. In a nutshell, the regular gap-affine distance can occasionally split long gaps by sporadic mismatches (most often, when aligning long and noisy sequences). Instead, many users would prefer to increase the gap open
cost to produce a single long gap instead. For that, the dual-cost gap-affine distance (p=2) defines two affine cost functions (i.e., for short and long gaps). Then, the alignment score is computed based on {X,O1,E1,O2,E2}⁠, where X corresponds to the mismatch penalty and the gap penalty is expressed as the function g(N)=min{O1+N·E1,O2+N·E2} (given the length of the gap N, the gap opening penalties O1 and O2, and the gap extension penalties E1 and E2).

<details><summary>Configuration Dual-Cost Gap-Affine Alignment</summary>
<p>

```C
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.mismatch = 6;       // X > 0
    attributes.affine2p_penalties.gap_opening1 = 4;   // O1 >= 0
    attributes.affine2p_penalties.gap_extension1 = 2; // E1 > 0
    attributes.affine2p_penalties.gap_opening2 = 12;  // O2 >= 0
    attributes.affine2p_penalties.gap_extension2 = 1; // E2 > 0
```

</p>
</details>

### <a name="wfa2.scope"></a> 2.2 Alignment Scope

Depending on the use case, it is often the case that an application only requires to compute the alignment score, not the complete alignment (i.e., CIGAR). As it happens with traditional dynamic programming algorithms, the WFA algorithm requires less memory (i.e., `O(s)`) to compute the alignment score. In turn, this results in slighter faster alignment executions, consuming less computational resources. For this reason, the WFA2 library implements two different modes depending on the alignment scope: score-only and full-CIGAR alignment. 

The ** score-only alignment ** mode computes only the alignment score. This mode utilizes only the front-wavefronts of the WFA algorithm to keep track of the optimal alignment score. As a result, it requires `O(s)` memory and, in practice, performs slighter faster than the standard full-CIGAR mode. 

<details><summary>Score-only configuration</summary>
<p>

```C
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.alignment_scope = compute_score;
```

</p>
</details>

The ** full-CIGAR alignment ** computes the sequence of alignment operations (i.e., {'M','X','D','I'}) that transforms one sequence into the other (i.e., alignment CIGAR). The alignment score can be obtained as a by-product of the alignment process, evaluating the score of the alignment CIGAR. This mode requires `O(s^2)` memory (using the default memory mode, wavefront_memory_high) or less (using the low-memory modes).

<details><summary>Full-CIGAR alignment configuration</summary>
<p>

```C
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.alignment_scope = compute_alignment;
```

</p>
</details>

### <a name="wfa2.span"></a> 2.3 Alignment Span

The WFA2 library allows computing alignments with different spans or shapes. Although there is certain ambiguity and confusion in the terminology, we have tried to generalized the different options available to offer flexible parameters that can capture multiple alignment scenarios. During the development of the WFA we decided to adhere to the classical approximate string matching terminology where we align a **pattern (a.k.a. query or sequence)** against a **text (a.k.a. target, database, or reference)**.

- **End-to-end alignment.** Also known as global-alignment, this alignment mode forces aligning the two sequences from the beginning to end of both.

```
    PATTERN    AATTAATTTAAGTCTAGGCTACTTTCGGTACTTTGTTCTT
               ||||    ||||||||||||||||||||||||||   |||
    TEXT       AATT----TAAGTCTAGGCTACTTTCGGTACTTT---CTT
```

<details><summary>End-to-end configuration</summary>
<p>

```C
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.alignment_form.span = alignment_end2end;
```

</p>
</details>

- **Ends-free alignment.** This alignment mode allows leading and trailing insertions or deletions for "free" (i.e., no penalty/cost on the overall alignment score). Moreover, this alignment mode allows determining the maximum gap length allowed for free at the beginning and end of the sequences. Note that this mode does not implement local alignment as it doesn't allow free insertions and deletions at the beginning/end of the sequences at the same time. However, it allows many different configurations used across different analysis, methods, and tools.

```
    PATTERN    AATTAATTTAAGTCTAGGCTACTTTCGGTACTTTGTTCTT
                   |||||||||||||||||||||||||||||| ||   
    TEXT       ----AATTTAAGTCTAGGCTACTTTCGGTACTTTCTT---
```

<details><summary>Ends-free configuration</summary>

```C
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.alignment_form.span = alignment_endsfree;
    attributes.alignment_form.pattern_begin_free = pattern_begin_free;
    attributes.alignment_form.pattern_end_free = pattern_end_free;
    attributes.alignment_form.text_begin_free = text_begin_free;
    attributes.alignment_form.text_end_free = text_end_free;
```

</details>

- **Other**

<details><summary>Glocal alignment (a.k.a. semi-global or fitting)</summary>
<p>

- **Glocal alignment (a.k.a. semi-global or fitting).** Alignment mode where the pattern is globally aligned and the text is locally aligned. Often due to the large size of one of the sequences (e.g., the text sequence being a genomic reference), this alignment mode forces one sequence (i.e., pattern) to align globally to a substring of the other (i.e., text).

```
    PATTERN    -------------AATTTAAGTCTAGGCTACTTTC---------------
                            ||||||||| ||||||||||||               
    TEXT       ACGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT
```

```C
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.alignment_form.span = alignment_endsfree;
    attributes.alignment_form.pattern_begin_free = 0;
    attributes.alignment_form.pattern_end_free = 0;
    attributes.alignment_form.text_begin_free = text_begin_free;
    attributes.alignment_form.text_end_free = text_end_free;
```

</p>
</details>

<details><summary>Extension alignment</summary>
<p>

- **Extension alignment.** Alignment mode where the start of both pattern and text sequences are forced to be aligned. However, the ends of both are free. This alignment mode is typically used within seed-and-extend algorithms.

```C
    // Right extension
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.alignment_form.span = alignment_endsfree;
    attributes.alignment_form.pattern_begin_free = 0;
    attributes.alignment_form.pattern_end_free = pattern_end_free;
    attributes.alignment_form.text_begin_free = 0;
    attributes.alignment_form.text_end_free = text_end_free;
    
    PATTERN    AATTTAAGTCTG-CTACTTTCACGCA-GCT----------
               ||||| |||||| ||||||||||| | | |          
    TEXT       AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT
```

```C
    // Left extension
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.alignment_form.span = alignment_endsfree;
    attributes.alignment_form.pattern_begin_free = pattern_begin_free;
    attributes.alignment_form.pattern_end_free = 0;
    attributes.alignment_form.text_begin_free = text_begin_free;
    attributes.alignment_form.text_end_free = 0;
    
    PATTERN    -------------AAACTTTCACGTACG-TGACAGTCTCT
                              ||||||||||||| |||||| ||||
    TEXT       AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT
```

</p>
</details>

<details><summary>Overlapped alignment</summary>
<p>

- **Overlapped alignment (a.k.a. dovetail).** 

```C
    // Overlapped (Right-Left)
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.alignment_form.span = alignment_endsfree;
    attributes.alignment_form.pattern_begin_free = pattern_begin_free;
    attributes.alignment_form.pattern_end_free = 0;
    attributes.alignment_form.text_begin_free = 0;
    attributes.alignment_form.text_end_free = text_end_free;
    
    PATTERN    ACGCGTCTGACTGACTGACTAAACTTTCATGTAC-TGACA-----------------
                                   ||||||||| |||| |||||                 
    TEXT       --------------------AAACTTTCACGTACGTGACATATAGCGATCGATGACT
```

```C
    // Overlapped (Left-Right)
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.alignment_form.span = alignment_endsfree;
    attributes.alignment_form.pattern_begin_free = 0;
    attributes.alignment_form.pattern_end_free = pattern_end_free;
    attributes.alignment_form.text_begin_free = text_begin_free;
    attributes.alignment_form.text_end_free = 0;

    PATTERN    ----------------------ACGCGTCTGACTGACTACGACTACGACTGACTAGCAT
                                     ||||||||| || ||                      
    TEXT       ACATGCATCGATCAGACTGACTACGCGTCTG-CTAAC----------------------
```

</p>
</details>

### <a name="wfa2.mem"></a> 2.4 Memory modes

The WFA2 library implements various memory modes: `wavefront_memory_high`, `wavefront_memory_med`, `wavefront_memory_low`. These modes allow regulating the overall memory consumption at the expense of execution time. The standard WFA algorithm, which stores explicitly all wavefronts in memory, correspond to the mode `wavefront_memory_high`. The other methods progressively reduce the memory usage at the expense of slightlier larger alignment times. These memory modes can be used transparently with other alignment options and generate identical results. Note that the score-only alignment is not affected by this option (using a minimal memory footprint of `O(s)`).

<details><summary>Memory-mode configuration</summary>
<p>

```C

  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.memory_mode = wavefront_memory_med;
  
```

</p>
</details>

### <a name="wfa2.heuristics"></a> 2.5 Heuristic modes

The WFA algorithm can be used together with many heuristics to reduced the alignment time and memory used. 
As it happens to other alignment methods, heuristics often result in suboptimal solutions and lost of accuracy.
Due to the 

The WFA2 library implements various of these methods.

- Adaptive-Wavefront alignment

- Banded alignment

### <a name="wfa2.other"></a> 2.7 Other options

- Pause-and-resume mode

- Parallel mode

- Maximum memory

- Lambda mode

- Plot

- Check

## <a name="wfa2.programming"></a> 3. PROGRAMMING WITH WFA2-LIB

### <a name="wfa2.programming.c"></a> 3.1 Simple C example

This simple example illustrates how to align two sequences using the WFA2 library. First, include the WFA2 alignment headers.

```C
#include "wavefront/wavefront_align.h"
```

Next, create and configure the WFA alignment object. The following example uses the defaults configuration and sets custom `gap_affine` penalties. Note that mismatch, gap-opening, and gap-extension must be positive values.  

```C
// Configure alignment attributes
wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
attributes.distance_metric = gap_affine;
attributes.affine_penalties.mismatch = 4;
attributes.affine_penalties.gap_opening = 6;
attributes.affine_penalties.gap_extension = 2;
// Initialize Wavefront Aligner
wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
```

Finally, call the `wavefront_align` function.

```C
// Patter & Text
char* pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
char* text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
// Align
wavefront_align(wf_aligner,pattern,strlen(pattern),text,strlen(text));
```

Afterwards, we can use the library to display the result of the alignment (e.g., the alignment score and CIGAR).

```C
// Display alignment
fprintf(stderr,"  Pattern  %s\n",pattern);
fprintf(stderr,"  Text     %s\n",text);
fprintf(stderr,"  WFA-Alignment Score %d\n",wf_aligner->cigar.score);
cigar_print_pretty(stderr,
    pattern,strlen(pattern),text,strlen(text),
    &wf_aligner->cigar,wf_aligner->mm_allocator);
```

At the end of the program, it is polite to release the memory used.

```C
// Free
wavefront_aligner_delete(wf_aligner);
```

To compile and run this example you need to link against the WFA library (-lwfa).

```
$> gcc -O3 wfa_example.c -o wfa_example -lwfa
$> ./wfa_example
```

**IMPORTANT.** Once an alignment object is created, **it is strongly recommended to reuse the alignment object to compute multiple alignment**. Creating and destroying the alignment object for every alignment computed can have a significant overhead. Reusing the alignment object allows repurposing internal datastructures, minimize the cost of memory allocations, avoid multiple alignment setups and precomputations, etc.

### <a name="wfa2.programming.cpp"></a> 3.2 Simple C++ example

The WFA2 library can be using from C++ code using the C++ bindings. This example is similar to the previous one but using the C++ bindings. First, include the C++ bindings and remember to use the WFA namespace.

```C
#include "bindings/cpp/WFAligner.hpp"

using namespace wfa;
```

Configure and create the WFA alignment object. In this case, gap-affine distance using custom penalties and the standard memory-usage algorithm (i.e., standard WFA algorithm).

```C++
// Create a WFAligner
WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);
```

Align two sequences (in this case, given as strings).

```C++
// Patter & Text
string pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
string text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
// Align
aligner.alignEnd2End(pattern,text);
```

Display the result of the alignment.

```C++
// Display score
cout << "WFA-Alignment returns score " << aligner.getAlignmentScore() << endl;
// Print CIGAR
string cigar = aligner.getAlignmentCigar();
cout << "PATTERN: " << pattern  << endl;
cout << "TEXT: " << text  << endl;
cout << "CIGAR: " << cigar  << endl;
```

**IMPORTANT.** Once an alignment object is created, **it is strongly recommended to reuse the alignment object to compute multiple alignment**. Creating and destroying the alignment object for every alignment computed can have a significant overhead. Reusing the alignment object allows repurposing internal datastructures, minimize the cost of memory allocations, avoid multiple alignment setups and precomputations, etc.

## <a name="wfa2.complains"></a> 4. REPORTING BUGS AND FEATURE REQUEST

Feedback and bug reporting is highly appreciated. Please report any issue or suggestion on github or via email to the main developer (santiagomsola@gmail.com).

## <a name="wfa2.licence"></a> 5. LICENSE

WFA2-lib is distributed under MIT licence.

## <a name="wfa2.cite"></a> 6. CITATION

**Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa**. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.
