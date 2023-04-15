#pragma once
#include "util.hpp"

namespace genasm_cpu {
    extern bool enabled_algorithm_log;
    std::vector<Alignment_t> align_all(Genome_t &reference, std::vector<Read_t> &reads, int threads=1, long long* core_algorithm_ns=NULL);
    std::vector<Alignment_t> align_all(std::vector<std::string> &texts, std::vector<std::string> &queries, int threads, long long* core_algorithm_ns=NULL);
}
