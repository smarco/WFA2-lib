#pragma once

#include <string>
#include <vector>
#include <chrono>
#include <map>

#define SEED_FILE_MAF 0
#define SEED_FILE_PAF 1

typedef struct Sequence {
    std::string description; // description line of sequence
    std::string content; // content of sequence without whitespaces
} Sequence_t;

typedef struct Genome {
    std::map<std::string, long long> chromosome_starts;
    std::string content;
} Genome_t;

typedef struct Read Read_t;
typedef struct CandidateLocation {
    std::string read_description;
    std::string chromosome;
    long long start_in_chromosome;
    long long start_in_reference;
    long long start_of_aligned_region;
    long long size_of_aligned_region;
    bool strand;
} CandidateLocation_t;

typedef struct Read {
    std::string description;
    std::string content; //content of sequence without whitespaces
    std::vector<CandidateLocation_t> locations;
} Read_t;

typedef struct Alignment {
    std::string cigar;
    long long edit_distance;
} Alignment_t;

typedef struct CigarEntry {
    uint8_t edit_count;
    char edit_type;
} CigarEntry_t;

template <typename TargetType> long long measure_ns(TargetType target){
    auto begin = std::chrono::high_resolution_clock::now();
    target();
    auto end = std::chrono::high_resolution_clock::now();

    long long ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
    return ns;
}

void remove_whitespaces(std::string &str);

std::string read_file(std::string file_path);

Sequence_t read_fasta(std::string file_path, int sequence_index);

Genome_t read_genome(std::string fasta_file_path);

void get_global_seeds(Genome_t &genome, std::vector<CandidateLocation_t> &locations);

Sequence_t read_fasta_highspeed(std::string file_path, int sequence_index);

std::vector<Read_t> read_fastq(std::string file_path);

std::vector<CandidateLocation_t> read_maf(std::string file_path);

std::vector<CandidateLocation_t> read_paf(std::string file_path);

bool ends_with(std::string const &s, std::string const &ending);

std::vector<Read_t> read_fastq_and_seed_locations(Genome_t &genome, std::string fastq_file_path, std::string seed_file_path, std::vector<Read_t> &reads);

bool cigar_char_equals(char c, char d);

#define OPT_MISSING 0
#define OPT_EXISTS 1
#define OPT_INVALID 2
int get_cmd_option(int argc, char **argv, std::string key);
int get_cmd_option(int argc, char **argv, std::string key, std::string &value);
bool check_options(int argc, char **argv, std::vector<std::string> valid_options);

std::vector<int> parse_csv_numbers(std::string csv);

std::vector<std::string> parse_csv_strings(std::string csv);
