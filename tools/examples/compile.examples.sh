#!/bin/bash

# Check library
if [ ! -f "../../build/libwfa.a" ]; then
  echo "Library libwfa.a not found. Please compile WFA library from top folder first"
  exit 
fi

# Compile examples
WFA_INCLUDES="../.."
WFA_LIB_PATH="../../build"

gcc -O3 -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_basic.c -o wfa_basic -lwfa
gcc -O3 -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_repeated.c -o wfa_repeated -lwfa
gcc -O3 -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_adapt.c -o wfa_adapt -lwfa

WFA_CPP="../../bindings/cpp/WFAligner.cpp" 

g++ -O3 -I $WFA_INCLUDES -L $WFA_LIB_PATH $WFA_CPP wfa_bindings.cpp -o wfa_bindings -lwfa
