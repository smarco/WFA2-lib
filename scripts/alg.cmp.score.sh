#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Compare alg files
# USAGE: ./alg.cmp.score.sh checkfile1.alg checkfile2.alg

# Parameters
FILE1=$1
FILE2=$1

# Compare
diff  <(awk '{print $1}' $FILE1) <(awk '{print $1}' $FILE2)
