#!/bin/bash -x
# PROJECT: Wavefront Alignments Algorithms
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: WFA unitary tests (for performance & correcness)
# USAGE: ./wfa.utest.run.sh

# Clear
rm *.log *.alg

# Config
ALGORITHM="gap-affine-wfa"  
REDUCTION="--wfa-heuristic=wfa-adaptive --wfa-heuristic-parameters 10,50,1"
LOWMEMORY="--wfa-memory-mode=med"

# Utest for length=100
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100.n100K.e2.seq -o sim.l100.e2.W.alg               &> sim.l100.e2.W.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100.n100K.e2.seq -o sim.l100.e2.Wl.alg $LOWMEMORY   &> sim.l100.e2.Wl.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100.n100K.e2.seq -o sim.l100.e2.Wa.alg $REDUCTION   &> sim.l100.e2.Wa.log

# Utest for length=1K
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l1K.n10K.e20.seq -o sim.l1K.e20.W.alg                &> sim.l1K.e20.W.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l1K.n10K.e20.seq -o sim.l1K.e20.Wl.alg $LOWMEMORY    &> sim.l1K.e20.Wl.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l1K.n10K.e20.seq -o sim.l1K.e20.Wa.alg $REDUCTION    &> sim.l1K.e20.Wa.log

# Utest for length=10K
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l10K.n1K.e20.seq -o sim.l10K.e20.W.alg -P 100        &> sim.l10K.e20.W.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l10K.n1K.e20.seq -o sim.l10K.e20.Wl.alg $LOWMEMORY   &> sim.l10K.e20.Wl.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l10K.n1K.e20.seq -o sim.l10K.e20.Wa.alg $REDUCTION   &> sim.l10K.e20.Wa.log

# Utest for length=100K
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100K.n1.e10.seq -o sim.l100K.e10.Wl.alg $LOWMEMORY  &> sim.l100K.e10.Wl.log

#
# Random tests
#
# INPUT="../data/sim.l1K.n10K.e20.seq -P 500"
# ./bin/align_benchmark -a indel-wfa -i $INPUT -c score --check-distance indel -v 
# ./bin/align_benchmark -a edit-wfa -i $INPUT -c score --check-distance edit -v
# ./bin/align_benchmark -a gap-linear-wfa -i $INPUT -c score --check-distance linear -v
# ./bin/align_benchmark -a gap-linear-wfa -i $INPUT -c score --check-distance linear -v
# ./bin/align_benchmark -a gap-affine-wfa -i $INPUT -c score --check-distance affine -v
# ./bin/align_benchmark -a gap-affine-2p-wfa -i $INPUT -c score --check-distance affine2p -v