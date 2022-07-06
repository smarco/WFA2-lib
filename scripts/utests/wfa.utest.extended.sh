#!/bin/bash -x
# PROJECT: Wavefront Alignments Algorithms (Unitary Tests)
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: WFA unitary tests (extended for different use-cases)
# USAGE: ./wfa.utest.extended.sh

# Clear
rm -rf *.log *.alg utests
mkdir utests

# Config
BIN="/usr/bin/time -v ./bin/align_benchmark"
BIWFA="--wfa-memory-mode=ultralow"
XTRAS="--check correct"

# Base data
DATA100="../data/sim.l100.n100K.e2.seq"
DATA1K="../data/sim.l1K.n10K.e20.seq"
DATA10K="../data/sim.l10K.n1K.e20.seq"

ENDSFREE1K="../data/sim.endsfree.l1K.n10K.seq"
ENDSFREE10K="../data/sim.endsfree.l10K.n1K.seq"

LONG100K="../data/sim.l100K.n10.e20.seq"

# Base (edit)
$BIN -a edit-wfa -i $DATA100 -o utests/test.l100.aedit.alg        $XTRAS &> utests/test.l100.aedit.log
$BIN -a edit-wfa -i $DATA1K  -o utests/test.l1K.aedit.alg         $XTRAS &> utests/test.l1K.aedit.log
$BIN -a edit-wfa -i $DATA10K -o utests/test.l10K.aedit.alg $BIWFA $XTRAS &> utests/test.l10K.aedit.log
$BIN -a edit-wfa -i $DATA100 -o utests/test.l100.sedit.alg --wfa-score-only $XTRAS &> utests/test.l100.sedit.log
$BIN -a edit-wfa -i $DATA1K  -o utests/test.l1K.sedit.alg  --wfa-score-only $XTRAS &> utests/test.l1K.sedit.log
$BIN -a edit-wfa -i $DATA10K -o utests/test.l10K.sedit.alg --wfa-score-only $XTRAS &> utests/test.l10K.sedit.log

# Base (gap-linear)
$BIN -a gap-linear-wfa -i $DATA100 -o utests/test.l100.alinear.alg        $XTRAS &> utests/test.l100.alinear.log
$BIN -a gap-linear-wfa -i $DATA1K  -o utests/test.l1K.alinear.alg         $XTRAS &> utests/test.l1K.alinear.log
$BIN -a gap-linear-wfa -i $DATA10K -o utests/test.l10K.alinear.alg $BIWFA $XTRAS &> utests/test.l10K.alinear.log
$BIN -a gap-linear-wfa -i $DATA100 -o utests/test.l100.slinear.alg --wfa-score-only $XTRAS &> utests/test.l100.slinear.log
$BIN -a gap-linear-wfa -i $DATA1K  -o utests/test.l1K.slinear.alg  --wfa-score-only $XTRAS &> utests/test.l1K.slinear.log
$BIN -a gap-linear-wfa -i $DATA10K -o utests/test.l10K.slinear.alg --wfa-score-only $XTRAS &> utests/test.l10K.slinear.log

# Base (gap-affine)
$BIN -a gap-affine-wfa -i $DATA100 -o utests/test.l100.aaffine.alg        $XTRAS &> utests/test.l100.aaffine.log
$BIN -a gap-affine-wfa -i $DATA1K  -o utests/test.l1K.aaffine.alg         $XTRAS &> utests/test.l1K.aaffine.log
#$BIN -a gap-affine-wfa -i $DATA10K -o utests/test.l10K.aaffine.alg $BIWFA $XTRAS &> utests/test.l10K.aaffine.log
$BIN -a gap-affine-wfa -i $DATA100 -o utests/test.l100.saffine.alg --wfa-score-only $XTRAS &> utests/test.l100.saffine.log
$BIN -a gap-affine-wfa -i $DATA1K  -o utests/test.l1K.saffine.alg  --wfa-score-only $XTRAS &> utests/test.l1K.saffine.log
#$BIN -a gap-affine-wfa -i $DATA10K -o utests/test.l10K.saffine.alg --wfa-score-only $XTRAS &> utests/test.l10K.saffine.log

# Base (gap-affine-2p)
$BIN -a gap-affine2p-wfa -i $DATA100 -o utests/test.l100.a2p.alg        $XTRAS &> utests/test.l100.a2p.log
$BIN -a gap-affine2p-wfa -i $DATA1K  -o utests/test.l1K.a2p.alg         $XTRAS &> utests/test.l1K.a2p.log
#$BIN -a gap-affine2p-wfa -i $DATA10K -o utests/test.l10K.a2p.alg $BIWFA $XTRAS &> utests/test.l10K.a2p.log
$BIN -a gap-affine2p-wfa -i $DATA100 -o utests/test.l100.s2p.alg --wfa-score-only $XTRAS &> utests/test.l100.s2p.log
$BIN -a gap-affine2p-wfa -i $DATA1K  -o utests/test.l1K.s2p.alg  --wfa-score-only $XTRAS &> utests/test.l1K.s2p.log
#$BIN -a gap-affine2p-wfa -i $DATA10K -o utests/test.l10K.s2p.alg --wfa-score-only $XTRAS &> utests/test.l10K.s2p.log

# Heuristics (gap-affine-2p alignment)
HEURISTIC="--wfa-heuristic=wfa-adaptive --wfa-heuristic-parameters=10,50,1"
$BIN -a gap-affine2p-wfa -i $DATA100 -o utests/test.adaptive.l100.alg $HEURISTIC $XTRAS &> utests/test.adaptive.l100.log
$BIN -a gap-affine2p-wfa -i $DATA1K  -o utests/test.adaptive.l1K.alg  $HEURISTIC $XTRAS &> utests/test.adaptive.l1K.log
$BIN -a gap-affine2p-wfa -i $DATA10K -o utests/test.adaptive.l10K.alg $HEURISTIC $XTRAS &> utests/test.adaptive.l10K.log

# Ends-free (length-diff=20%,e=10%,indels=2,20%)
$BIN -a gap-affine2p-wfa -i $ENDSFREE1K  -o utests/test.endsfree.l1K.alg  $BIWFA --ends-free=200,200,200,200     $XTRAS &> utests/test.endsfree.l1K.log
$BIN -a gap-affine2p-wfa -i $ENDSFREE10K -o utests/test.endsfree.l10K.alg $BIWFA --ends-free=2000,2000,2000,2000 $XTRAS &> utests/test.endsfree.l10K.log

# Very long sequences (100K)
#$BIN -a gap-affine2p-wfa -i $LONG100K -o utests/test.long.l100K.alg $BIWFA $XTRAS &> utests/test.long.l100K.log

 
