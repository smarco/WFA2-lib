#pragma once

#include "util.hpp"

namespace genasm_gpu {
    extern bool enabled_algorithm_log;
    std::vector<Alignment_t> align_all(Genome_t &reference, std::vector<Read_t> &reads, long long* core_algorithm_ns=NULL);
    std::vector<Alignment_t> align_all(std::vector<std::string> &texts, std::vector<std::string> &queries, long long* core_algorithm_ns=NULL);
    __global__ void ascii_to_twobit_strings(int count, long long *string_lengths, char **ascii_strings, char **twobit_strings);
}
