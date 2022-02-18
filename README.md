# WFA

## 1. INTRODUCTION

### 1.1 What is WFA?

The wavefront alignment (WFA) algorithm is an **exact** gap-affine algorithm that takes advantage of homologous regions between the sequences to accelerate the alignment process. As opposed to traditional dynamic programming algorithms that run in quadratic time, the WFA runs in time O(ns+s^2), proportional to the read length n and the alignment score s, using O(s^2) memory.  oreover, the WFA exhibits simple data dependencies that can be easily vectorized, even by the automatic features of modern compilers, for different architectures, without the need to adapt the code.

### 1.2 What is WFA2-lib?

The WFA2-lib library implements the WFA algorithm for different distance metrics and alignment modes. 
It supports indel, edit, gap-lineal, gap-affine, and gap-affine distances, aligning end-to-end or ends-free.


Moreover it provides support functions to display and verify the results.





### 1.3 Getting started

Note: We recommend using the GCC compiler

Clone GIT and compile

```
$> git clone https://github.com/smarco/WFA.git WFA
$> cd WFA
$> make clean all
```

### 1.4 Simple C example



## 2. WFA2 LIBRARY FEATURES

### 2.1 Distance Metrics

### 2.2 Alignment Scope

### 2.3 Alignment Span

### 2.4 Memory modes

### 2.5 Heuristic modes

### 2.6 Lambda mode







## 3. OTHER OPTIONS 








## 4. AUTHORS

  Santiago Marco-Sola \- santiagomsola@gmail.com
  
## 5. REPORTING BUGS

Feedback and bug reporting it's highly appreciated. Please report any issue or suggestion on github, or by email to the main developer (santiagomsola@gmail.com).

## 6. LICENSE

WFA2-lib is distributed under MIT licence.

## 7. CITATION

**Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa**. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.
