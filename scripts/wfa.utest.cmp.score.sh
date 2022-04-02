#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms (Unitary Tests)
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Compare alignment files (*.alg)
# USAGE: ./wfa.utest.cmp.score.sh file1.alg file2.alg

# Parameters
FILE1=$1
FILE2=$1

# Compare
diff  <(awk '{print $1}' $FILE1) <(awk '{print $1}' $FILE2)
