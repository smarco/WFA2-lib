#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms (Unitary Tests)
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: WFA unitary tests (Summarise performance results)
# USAGE: ./wfa.utest.summary.sh

for FILE_ALG in $1/*.alg
do
  PREFIX=${FILE_ALG%.*}
  FILE_LOG=$PREFIX.log
  echo "[UTest::$FILE_ALG] "
  echo -en "  Time  \t"
  grep -m1 "Time.Alignment" $FILE_LOG | awk '{print $3" "$4}'
  echo -en "  Memory\t"
  grep -m1 "Maximum resident set size" $FILE_LOG | tr -d "(:)" | awk '{print $6" "$5}'
done 