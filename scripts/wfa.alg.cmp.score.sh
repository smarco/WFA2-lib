#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms 
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Compare alignment files (*.alg)
# USAGE: ./wfa.alg.cmp.score.sh file1.alg file2.alg

# Parameters
FILE1=$1
FILE2=$1

# Compare
diff  <(awk '{if ($1<0) print -$1; else print $1}' $FILE1) <(awk '{if ($1<0) print -$1; else print $1}' $FILE2)
