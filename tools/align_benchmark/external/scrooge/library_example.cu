#include <string>
#include <vector>
#include <iostream>

#include "genasm_cpu.hpp"
#include "genasm_gpu.hpp"
#include "util.hpp"

using namespace std;

void cpu_string_pairs_example(){
    vector<string> texts = {"ACGTACGT"};
    vector<string> queries = {"ACGTACG"};

    int threads = 1;
    vector<Alignment_t> alignments = genasm_cpu::align_all(texts, queries, threads);

    for(Alignment_t &aln : alignments){
        cout << "edit_distance:" << aln.edit_distance << " ";
        cout << "cigar:" << aln.cigar << endl;
    }
}

void gpu_string_pairs_example(){
    vector<string> texts = {"ACGTACGT"};
    vector<string> queries = {"ACGTACG"};

    vector<Alignment_t> alignments = genasm_gpu::align_all(texts, queries);

    for(Alignment_t &aln : alignments){
        cout << "edit_distance:" << aln.edit_distance << " ";
        cout << "cigar:" << aln.cigar << endl;
    }
}

void cpu_mapping_example(){
    Genome_t reference;
    reference.content = "ACGTACGT";

    //candidate locations can be anywhere in the reference genome
    //each read can have its own set of multiple candidate locations
    //here we use only a single single set with a single candidate location at the start of the reference genome, for demonstration purposes
    CandidateLocation_t ref_begin;
    ref_begin.start_in_reference = 0; //first character index of the candidate location
    ref_begin.strand = true; //forward strand
    vector<CandidateLocation_t> ref_begin_vec(1, ref_begin);

    Read_t read;
    read.description = "example_read_id";
    read.content = "ACGTACG";
    read.locations = ref_begin_vec;
    vector<Read_t> reads(1, read);

    int threads = 1;
    vector<Alignment_t> alignments = genasm_cpu::align_all(reference, reads, threads);

    for(Alignment_t &aln : alignments){
        cout << "edit_distance:" << aln.edit_distance << " ";
        cout << "cigar:" << aln.cigar << endl;
    }
}

void gpu_mapping_example(){
    Genome_t reference;
    reference.content = "ACGTACGT";

    //candidate locations can be anywhere in the reference genome
    //each read can have its own set of multiple candidate locations
    //here we use only a single single set with a single candidate location at the start of the reference genome, for demonstration purposes
    CandidateLocation_t ref_begin;
    ref_begin.start_in_reference = 0; //first character index of the candidate location
    ref_begin.strand = true; //forward strand
    vector<CandidateLocation_t> ref_begin_vec(1, ref_begin);

    
    Read_t read;
    read.description = "example_read_id";
    read.content = "ACGTACG";
    read.locations = ref_begin_vec;
    vector<Read_t> reads(1, read);

    vector<Alignment_t> alignments = genasm_gpu::align_all(reference, reads);

    for(Alignment_t &aln : alignments){
        cout << "edit_distance:" << aln.edit_distance << " ";
        cout << "cigar:" << aln.cigar << endl;
    }
}

int main(int argc, char *argv[]){
    genasm_cpu::enabled_algorithm_log = false;
    genasm_gpu::enabled_algorithm_log = false;
    
    cpu_string_pairs_example();
    gpu_string_pairs_example();
    cpu_mapping_example();
    gpu_mapping_example();
}
