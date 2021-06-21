#!/bin/bash

# Check library
if [ ! -f "../../build/libwfa.a" ]; then
  echo "Library libwfa.a not found. Please compile WFA library from top folder first"
  exit 
fi

# Compile examples
WFA_INCLUDES="../.."
WFA_LIB_PATCH="../../build"

gcc -O3 -I $WFA_INCLUDES -L $WFA_LIB_PATCH wfa_basic.c -o wfa_basic -lwfa
gcc -O3 -I $WFA_INCLUDES -L $WFA_LIB_PATCH wfa_repeated.c -o wfa_repeated -lwfa
gcc -O3 -I $WFA_INCLUDES -L $WFA_LIB_PATCH wfa_adapt.c -o wfa_adapt -lwfa
