#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms (Unitary Tests)
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Compares alignments (*.alg) from two folders
# USAGE: ./wfa.utest.cmp.sh wfa_results_folder_1 wfa_results_folder_2

# Parameters
FOLDER1=$1
FOLDER2=$2

echo "> Comparing $FOLDER1 vs $FOLDER2"
for FILE_ALG1 in $FOLDER1/*.alg
do
  FILENAME=$(basename -- "$FILE_ALG1")
  PREFIX=${FILENAME%.*}
  FILE_ALG2="$FOLDER2/$FILENAME"
  echo -ne "[UTest::$PREFIX]"
  if [[ ${#PREFIX} < 15 ]]; then echo -ne "   "; fi
  echo -ne "\t"
  # Check existence
  if [[ ! -f "$FILE_ALG2" ]]
  then
    echo "$FILE_ALG2 doesn't exist."
    continue
  fi
  # Check diff
  if [[ $(diff $FILE_ALG1 $FILE_ALG2) ]] 
  then
    if [[ $(diff <(awk '{if ($1<0) print -$1; else print $1}' $FILE_ALG1) <(awk '{if ($1<0) print -$1; else print $1}' $FILE_ALG2)) ]]
    then
      echo "Error"
    else
      echo "ok" # Only score
    fi
  else
    echo "OK"
  fi
done
