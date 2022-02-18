#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: WFA unitary tests (for performance & correcness)
# USAGE: ./wfa.utest.memcheck.sh

# Log file
echo "" > mem.asan.report

# Reduce dataset
head -n 2000 ../data/sim.l1K.n10K.e20.seq > mem.asan.data

# Valgrind memcheck
make clean all
valgrind --tool=memcheck ./bin/align_benchmark -a gap-affine-wfa -i mem.asan.data              &>> mem.asan.report
valgrind --tool=memcheck ./bin/align_benchmark -a gap-affine-wfa -i mem.asan.data --low-memory &>> mem.asan.report 
valgrind --tool=memcheck ./bin/align_benchmark -a gap-affine-wfa-adaptive -i mem.asan.data     &>> mem.asan.report

# ASAN
make clean asan
ASAN_OPTIONS=detect_leaks=1:symbolize=1 LSAN_OPTIONS=verbosity=2:log_threads=1
./bin/align_benchmark -a gap-affine-wfa -i mem.asan.data              &>> mem.asan.report
./bin/align_benchmark -a gap-affine-wfa -i mem.asan.data --low-memory &>> mem.asan.report
./bin/align_benchmark -a gap-affine-wfa-adaptive -i mem.asan.data     &>> mem.asan.report

# Clean
rm mem.asan.data
