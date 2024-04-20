#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms (Unitary Tests)
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Cleanup
# USAGE: ./wfa.utest.clean.sh

# Config
OUTPUT="./tests"

# Clear
rm -f $OUTPUT/*.alg $OUTPUT/*.log* &> /dev/null

