#!/bin/bash -x
# PROJECT: Wavefront Alignments Algorithms (Unitary Tests)
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: WFA unitary tests (for heuristics)
# USAGE: ./wfa.utest.heuristics.sh

SINGLE_PAIR=$1

# Clean
rm *.alg *.wfa *.png

MODE="full"
HEURISTIC=""
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot

MODE="band.10.10"
HEURISTIC="--wfa-heuristic=banded-static --wfa-heuristic-parameters=-10,10"
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot
MODE="band.10.150"
HEURISTIC="--wfa-heuristic=banded-static --wfa-heuristic-parameters=-10,150"
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot

MODE="aband.10.10"
HEURISTIC="--wfa-heuristic=banded-adaptive --wfa-heuristic-parameters=-10,10,1"
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot
MODE="aband.50.50"
HEURISTIC="--wfa-heuristic=banded-adaptive --wfa-heuristic-parameters=-50,50,1"
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot

MODE="wf-adap.10.50.1"
HEURISTIC="--wfa-heuristic=wfa-adaptive --wfa-heuristic-parameters=10,50,1"
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot
MODE="wf-adap.10.10.1"
HEURISTIC="--wfa-heuristic=wfa-adaptive --wfa-heuristic-parameters=10,10,1"
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot
MODE="wf-adap.10.50.10"
HEURISTIC="--wfa-heuristic=wfa-adaptive --wfa-heuristic-parameters=10,50,10"
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot

MODE="xdrop.20"
HEURISTIC="--wfa-heuristic=xdrop --wfa-heuristic-parameters=20,20"
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot
MODE="zdrop.20"
HEURISTIC="--wfa-heuristic=zdrop --wfa-heuristic-parameters=20,20"
./bin/align_benchmark -a gap-affine-wfa -i $SINGLE_PAIR -o $SINGLE_PAIR.$MODE.alg $HEURISTIC --plot

# Produce PNG
python3 ./scripts/wfa2png.py --compact

# Crop
# for i in *.png; do convert $i -crop 8118x6602+0+1000 ${i%.*}.cropped.png; done
