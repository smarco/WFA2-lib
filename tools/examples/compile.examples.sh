#!/bin/bash

# Check library
if [[ ! -f "../../lib/libwfa.a" || ! -f "../../lib/libwfacpp.a" ]]; then
  echo "Library libwfa.a/libwfacpp.a not found. Please compile WFA library from top folder first"
  exit 
fi

# Compile examples
WFA_INCLUDES="../.."
WFA_LIB_PATH="../../lib"
FLAGS="-O3"

gcc $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_basic.c -o wfa_basic -lwfa
gcc $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_repeated.c -o wfa_repeated -lwfa
gcc $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_adapt.c -o wfa_adapt -lwfa

g++ $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_bindings.cpp -o wfa_bindings -lwfacpp
g++ $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_lambda.cpp -o wfa_lambda -lwfacpp
