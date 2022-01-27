#!/bin/bash -x
# PROJECT: Wavefront Alignments Algorithms
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: WFA unitary tests (for performance & correcness)
# USAGE: ./wfa.utest.run.sh

rm *.log *.alg  

# Utest for length=100
\time -v ./bin/align_benchmark -a gap-affine-wfa -i ../data/sim.l100.n100K.e2.seq -o sim.l100.e2.W.alg                        &> sim.l100.e2.W.log
\time -v ./bin/align_benchmark -a gap-affine-wfa -i ../data/sim.l100.n100K.e2.seq --memory-mode=high -o sim.l100.e2.Wl.alg    &> sim.l100.e2.Wl.log
\time -v ./bin/align_benchmark -a gap-affine-wfa-adaptive -i ../data/sim.l100.n100K.e2.seq -o sim.l100.e2.Wa.alg              &> sim.l100.e2.Wa.log

# Utest for length=1K
\time -v ./bin/align_benchmark -a gap-affine-wfa -i ../data/sim.l1K.n10K.e20.seq -o sim.l1K.e20.W.alg                         &> sim.l1K.e20.W.log
\time -v ./bin/align_benchmark -a gap-affine-wfa -i ../data/sim.l1K.n10K.e20.seq --memory-mode=high -o sim.l1K.e20.Wl.alg     &> sim.l1K.e20.Wl.log
\time -v ./bin/align_benchmark -a gap-affine-wfa-adaptive -i ../data/sim.l1K.n10K.e20.seq -o sim.l1K.e20.Wa.alg               &> sim.l1K.e20.Wa.log

# Utest for length=10K
\time -v ./bin/align_benchmark -a gap-affine-wfa -i ../data/sim.l10K.n1K.e20.seq -o sim.l10K.e20.W.alg -P 100                 &> sim.l10K.e20.W.log
\time -v ./bin/align_benchmark -a gap-affine-wfa -i ../data/sim.l10K.n1K.e20.seq --memory-mode=high -o sim.l10K.e20.Wl.alg    &> sim.l10K.e20.Wl.log
\time -v ./bin/align_benchmark -a gap-affine-wfa-adaptive -i ../data/sim.l10K.n1K.e20.seq -o sim.l10K.e20.Wa.alg              &> sim.l10K.e20.Wa.log

# Other
\time -v ./bin/align_benchmark -a gap-affine-wfa   -i ../data/sim.l100K.n1.e10.seq -o sim.l100K.e10.Wl.alg --memory-mode=high &> sim.l100K.e10.Wl.log
\time -v ./bin/align_benchmark -a gap-affine2p-wfa -i ../data/sim.l1K.n10K.e20.seq -o sim.l1K.e20.2pW.alg -P 100              &> sim.l1K.e20.2pW.log

