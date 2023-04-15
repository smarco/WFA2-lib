#include "util.hpp"
#include "genasm_gpu.hpp"
#include "genasm_cpu.hpp"
#include "bitvector_test.hpp"

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>

using namespace std;

#ifdef __GNUC__
    #include <experimental/filesystem>
    using namespace std::experimental;
#else
    #include <filesystem>
#endif


bool enable_log = true;
string optimization_blocker = "";

bool cigarFormatCorrect(Alignment_t &alignment){
    stringstream cigar_ss(alignment.cigar);
    cigar_ss.peek(); //set the eof bit
    while(!cigar_ss.eof()){
        unsigned int edit_count;
        char edit_type;

        cigar_ss >> edit_count;
        cigar_ss >> edit_type;
        cigar_ss.peek(); //set the eof bit

        if(cigar_ss.fail()){
            cout << "CIGAR had bad format" << endl;
            return false;
        }

        if(edit_count == 0){
            cout << "CIGAR cannot contain edits with count 0" << endl;
            return false;
        }

        if(edit_type == 'I'){
        }
        else if(edit_type == 'D'){
        }
        else if(edit_type == 'X' || edit_type == '=' || edit_type == 'M'){
        }
        else{
            cout << "CIGAR contains unknown edit type '" << edit_type << "'" << endl;
            return false;
        }
    }
    return true;
}

bool cigarInBoundsAndCoversRead(Alignment_t &alignment, CandidateLocation_t &location, Read_t &read, Genome_t &reference){
    long long i = location.start_in_reference;
    long long j = 0;

    stringstream cigar_ss(alignment.cigar);
    cigar_ss.peek(); //set the eof bit
    while(!cigar_ss.eof()){
        unsigned int edit_count;
        char edit_type;
        cigar_ss >> edit_count;
        cigar_ss >> edit_type;
        
        cigar_ss.peek(); //set the eof bit

        if(edit_type == 'I'){
            j += edit_count;
        }
        else if(edit_type == 'D'){
            i += edit_count;
        }
        else{
            i += edit_count;
            j += edit_count;
        }
    }

    if(j < read.content.size()){
        cout << "CIGAR didn't cover entire read" << endl;
        return false;
    }

    if(j > read.content.size()){
        cout << "CIGAR went out of bounds of read" << endl;
        return false;
    }

    if(i > reference.content.size()){
        cout << "CIGAR went out of bounds of reference" << endl;
        return false;
    }

    return true;
}

bool validateCigarString(Alignment_t &alignment, CandidateLocation_t &location, Read_t &read, Genome_t &reference){
    /*
     * test if the given CIGAR string is a correct transformation from reference to read
     * return true if correct, false otherwise, and print potential error messages
    */
    if(!cigarFormatCorrect(alignment)){
        cout << "CIGAR format wrong" << endl;
        return false;
    }
    if(!cigarInBoundsAndCoversRead(alignment, location, read, reference)){
        cout << "CIGAR runs out of bounds or too short to cover read" << endl;
        return false;
    }

    long long i = location.start_in_reference;
    long long j = 0;
    long long edits_in_cigar_string = 0;

    stringstream cigar_ss(alignment.cigar);
    cigar_ss.peek(); //would set the eof bit for the empty string
    while(!cigar_ss.eof()){
        unsigned int edit_count;
        char edit_type;
        cigar_ss >> edit_count;
        cigar_ss >> edit_type;
        cigar_ss.peek(); //set the eof bit

        if(edit_type == 'I'){
            j += edit_count;
            edits_in_cigar_string += edit_count;
        }
        else if(edit_type == 'D'){
            i += edit_count;
            edits_in_cigar_string += edit_count;
        }
        else{ //M, X, =
            for(int e = 0; e < edit_count; e++){
                if(edit_type == 'X' && cigar_char_equals(reference.content[i], read.content[j])){
                    cout << "CIGAR contains 'X' but reference[i] and read[j] match" << endl;
                    return false;
                }
                if(edit_type == '=' && !cigar_char_equals(reference.content[i], read.content[j])){
                    cout << "CIGAR contains '=' but reference[i] and read[j] mismatch" << endl;
                    return false;
                }
                if(edit_type == 'M' && reference.content[i] != read.content[j]){
                    edits_in_cigar_string++;
                }
                i++;
                j++;
            }
            if(edit_type == 'X'){
                edits_in_cigar_string += edit_count;
            }
        }
    }

    if(edits_in_cigar_string != alignment.edit_distance){
        cout << "CIGAR has " << edits_in_cigar_string << " edits, while the reported edit disatance is " << alignment.edit_distance << endl;
        return false;
    }

    return true;
}

void gpu_algorithm_correctness_test(){
    Genome_t reference;
    reference.content = "AAAACCCCGGGGTTTT";
    
    CandidateLocation_t ref_begin;
    ref_begin.start_in_reference = 0;
    ref_begin.start_in_chromosome = 0;
    ref_begin.strand = true;
    ref_begin.chromosome = "";
    vector<CandidateLocation_t> ref_begin_vec(1, ref_begin);

    vector<Read_t> reads;
    reads.push_back({"test_read_4d12m4i",           "CCCCGGGGTTTTAAAA",         ref_begin_vec});
    reads.push_back({"test_read_16m",               "AAAACCCCGGGGTTTT",         ref_begin_vec});
    reads.push_back({"test_read_3d7m",              "ACCCCGG",                  ref_begin_vec});
    reads.push_back({"test_read_4m4d4m4i4m",        "AAAAGGGGAAAATTTT",         ref_begin_vec});
    reads.push_back({"test_read_12s4m",             "AAAAAAAAAAAAAAAA",         ref_begin_vec});
    reads.push_back({"test_read_1m1s1i3m1s2m3i",    "ATTAACGCCTTT",             ref_begin_vec});
    reads.push_back({"test_read_oversized",         "TTTTAAAACCCCGGGGTTTTAAAA", ref_begin_vec});
    reads.push_back({"test_read_empty",             "",                         ref_begin_vec});
    reads.push_back({"test_read_len64",             "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAACCCCGGGGTTTTAAAA", ref_begin_vec});

    //vector<int> correct_edit_distances = {8, 0, 3, 8, 12, 6, 8, 0};
    vector<int> correct_edit_distances = {8, 0, 3, 8, 12, 6, 8, 0, 48};
    //vector<int> correct_edit_distances = {8};
    //vector<int> correct_edit_distances = {0};
    //vector<int> correct_edit_distances = {48};
    
    vector<Alignment_t> alignments = genasm_gpu::align_all(reference, reads);

    if(alignments.size() != correct_edit_distances.size()){
        cout << "FAILED gpu_algorithm_correctness_test: align_all() produced wrong number of alignments" << endl;
        return;
    }

    bool success = true;
    for(int i = 0; i < alignments.size(); i++){
        if(alignments[i].edit_distance != correct_edit_distances[i]){
            cout << "FAILED gpu_algorithm_correctness_test: align_all() produced distance " << alignments[i].edit_distance;
            cout << " instead of " << correct_edit_distances[i];
            cout << " for read \""<< reads[i].description << "\"" << endl;
            success = false;
        }
        if(!validateCigarString(alignments[i], ref_begin, reads[i], reference)){
            success = false;
        }

    }
    if(success){
        cout << "PASSED gpu_algorithm_correctness_test" << endl;
    }
}

void cpu_algorithm_correctness_test(){
    Genome_t reference;
    reference.content = "AAAACCCCGGGGTTTT";
    
    CandidateLocation_t ref_begin;
    ref_begin.start_in_reference = 0;
    ref_begin.start_in_chromosome = 0;
    ref_begin.strand = true;
    ref_begin.chromosome = "";
    vector<CandidateLocation_t> ref_begin_vec(1, ref_begin);

    vector<Read_t> reads;
    reads.push_back({"test_read_4d12m4i",           "CCCCGGGGTTTTAAAA",         ref_begin_vec});
    reads.push_back({"test_read_16m",               "AAAACCCCGGGGTTTT",         ref_begin_vec});
    reads.push_back({"test_read_3d7m",              "ACCCCGG",                  ref_begin_vec});
    reads.push_back({"test_read_4m4d4m4i4m",        "AAAAGGGGAAAATTTT",         ref_begin_vec});
    reads.push_back({"test_read_12s4m",             "AAAAAAAAAAAAAAAA",         ref_begin_vec});
    reads.push_back({"test_read_1m1s1i3m1s2m3i",    "ATTAACGCCTTT",             ref_begin_vec});
    reads.push_back({"test_read_oversized",         "TTTTAAAACCCCGGGGTTTTAAAA", ref_begin_vec});
    reads.push_back({"test_read_empty",             "",                         ref_begin_vec});
    reads.push_back({"test_read_len64",             "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAACCCCGGGGTTTTAAAA", ref_begin_vec});

    vector<int> correct_edit_distances = {8, 0, 3, 8, 12, 6, 8, 0, 48};

    vector<Alignment_t> alignments = genasm_cpu::align_all(reference, reads);

    if(alignments.size() != correct_edit_distances.size()){
        cout << "FAILED cpu_algorithm_correctness_test: align_all() produced wrong number of alignments" << endl;
        return;
    }

    bool success = true;
    for(int i = 0; i < alignments.size(); i++){
        if(alignments[i].edit_distance != correct_edit_distances[i]){
            cout << "FAILED cpu_algorithm_correctness_test: align_all() produced distance " << alignments[i].edit_distance;
            cout << " instead of " << correct_edit_distances[i];
            cout << " for read \""<< reads[i].description << "\"" << endl;
            success = false;
        }
        if(!validateCigarString(alignments[i], ref_begin, reads[i], reference)){
            success = false;
        }

    }
    if(success){
        cout << "PASSED cpu_algorithm_correctness_test" << endl;
    }
}

void gpu_algorithm_performance_test(string reference_file_path, string reads_file_path, string alignments_file_path, int read_length_cap=-1, int dataset_inflation=1){
    if(enable_log) cerr << "Starting performance test..." << endl;

    if(enable_log) cerr << "Reading reference sequence..." << endl;
    Genome_t reference_genome = read_genome(reference_file_path);
    
    if(enable_log) cerr << "Reading reads files (~30 seconds)..." << endl;
    vector<Read_t> reads;
    read_fastq_and_seed_locations(reference_genome, reads_file_path, alignments_file_path, reads);

    if(enable_log) cerr << "Filtering reads..." << endl;
    //filter out any reverse complement reads
    for(Read_t &read : reads){
        read.locations.erase(remove_if(
                read.locations.begin(),
                read.locations.end(),
                [](CandidateLocation_t const &l){ return l.strand==false; }
            ),
            read.locations.end()
        );
    }

    //reads.erase(reads.begin()+1, reads.end());
    //reads.erase(reads.begin()+1000, reads.end());
    //reads = vector<Read_t>(115241, reads[0]);

    if(read_length_cap >= 0){
        for(Read_t &r: reads){
            r.content = r.content.substr(0, read_length_cap);
        }
    }

    if(dataset_inflation > 1){
        int old_size = reads.size();
        reads.resize(dataset_inflation*old_size);
        for(int i = 1; i < dataset_inflation; i++){
            copy_n(reads.begin(), old_size, reads.begin()+i*old_size);
        }
    }

    if(enable_log) cerr << "Sorting reads..." << endl;
    //sort reads in descending length
    sort(reads.begin(), reads.end(), [](Read_t &a, Read_t &b){return a.content.size() > b.content.size();});

    if(enable_log) cerr << "Running alignment algorithm..." << endl;

    vector<Alignment_t> alignments;
    long long core_algorithm_ns;
    auto workload = [&](){
        alignments = genasm_gpu::align_all(reference_genome, reads, &core_algorithm_ns);
    };
    long long end_to_end_ns = measure_ns(workload); //runtime including data transfer, conversion, readout, post-processing
    long long end_to_end_alignments_per_second = alignments.size() * 1000000000 / end_to_end_ns;
    long long core_algorithm_alignments_per_second = alignments.size() * 1000000000 / core_algorithm_ns;

    if(enable_log) cerr << "Sanity checking alignments..." << endl;
    
    size_t pair_idx = 0;
    for(Read_t &read : reads){
        for(CandidateLocation_t &location : read.locations){
            if(!validateCigarString(alignments[pair_idx], location, read, reference_genome)){
                cout << "FAILED sanity check in algorithm_performance_test ";
                cout << "for alignment " << pair_idx << endl;
            }
            pair_idx++;
        }
    }

    if(enable_log) cerr << "Done" << endl;

    cout << "align_all() took " << (end_to_end_ns/1000000) << "ms (data transfers, conversion, gpu kernel and post-processing)" << endl;
    cout << "GPU kernel took " << (core_algorithm_ns/1000000) << "ms" << endl;
    cout << "GPU kernel ran at " << core_algorithm_alignments_per_second << " aligns/second" << endl;
    //cout << "ran at " << end_to_end_alignments_per_second << " aligns/second" << endl;
}

void cpu_algorithm_performance_test(string reference_file_path, string reads_file_path, string alignments_file_path, int threads, int read_length_cap=-1, int dataset_inflation=1){
    if(enable_log) cerr << "Starting performance test..." << endl;

    if(enable_log) cerr << "Reading reference sequence..." << endl;
    Genome_t reference_genome = read_genome(reference_file_path);
    if(enable_log) cerr << "Reading reads files (~30 seconds)..." << endl;
    vector<Read_t> reads;
    read_fastq_and_seed_locations(reference_genome, reads_file_path, alignments_file_path, reads);

    if(enable_log) cerr << "Filtering reads..." << endl;
    //filter out any reverse complement reads
    for(Read_t &read : reads){
        read.locations.erase(remove_if(
                read.locations.begin(),
                read.locations.end(),
                [](CandidateLocation_t const &l){ return l.strand==false; }
            ),
            read.locations.end()
        );
    }
    //reads.erase(reads.begin()+1, reads.end());
    //reads.erase(reads.begin()+1000, reads.end());
    //reads = vector<Read_t>(115241, reads[0]);
    /*for(int i = 0; i < 115241 - 3; i+=4){
        reads[i+1] = reads[i];
        reads[i+2] = reads[i];
        reads[i+3] = reads[i];
    }*/
    if(read_length_cap >= 0){
        for(Read_t &r: reads){
            r.content = r.content.substr(0, read_length_cap);
        }
    }

    if(dataset_inflation > 1){
        int old_size = reads.size();
        reads.resize(dataset_inflation*old_size);
        for(int i = 1; i < dataset_inflation; i++){
            copy_n(reads.begin(), old_size, reads.begin()+i*old_size);
        }
    }

    if(enable_log) cerr << "Sorting reads..." << endl;
    //sort reads in descending length
    sort(reads.begin(), reads.end(), [](Read_t &a, Read_t &b){return a.content.size() > b.content.size();});

    if(enable_log) cerr << "Running alignment algorithm..." << endl;

    vector<Alignment_t> alignments;
    long long core_algorithm_ns;
    auto workload = [&](){
        alignments = genasm_cpu::align_all(reference_genome, reads, threads, &core_algorithm_ns);
    };
    long long end_to_end_ns = measure_ns(workload); //runtime including data transfer, conversion, readout, post-processing
    long long end_to_end_alignments_per_second = alignments.size() * 1000000000 / end_to_end_ns;
    long long core_algorithm_alignments_per_second = alignments.size() * 1000000000 / core_algorithm_ns;

    if(enable_log) cerr << "Sanity checking alignments..." << endl;
    
    size_t pair_idx = 0;
    for(Read_t &read : reads){
        for(CandidateLocation_t &location : read.locations){
            if(!validateCigarString(alignments[pair_idx], location, read, reference_genome)){
                cout << "FAILED sanity check in cpu_algorithm_performance_test ";
                cout << "for alignment " << pair_idx << endl;
            }
            pair_idx++;
        }
    }

    if(enable_log) cerr << "Done" << endl;

    cout << "align_all() took " << (end_to_end_ns/1000000) << "ms (data transfers, conversion, cpu kernel and post-processing)" << endl;
    cout << "CPU kernel took " << (core_algorithm_ns/1000000) << "ms" << endl;
    cout << "CPU kernel ran at " << core_algorithm_alignments_per_second << " aligns/second" << endl;
    //cout << "ran at " << end_to_end_alignments_per_second << " aligns/second" << endl;
}

void read_file_performance_test(string path){
    string raw_file;
    auto workload = [&path, &raw_file](){
        raw_file = read_file(path);
    };
    long long ns = measure_ns(workload);
    long long bytes_per_second = raw_file.size() * 1000000000 / ns;
    cout << "read_file() ran at " << bytes_per_second/1000000 << "MB/s in " << ns/1000000 << "ms"<< endl;
}

void read_genome_performance_test(string path){
    Genome_t reference_genome;
    auto workload = [&path, &reference_genome](){
        reference_genome = read_genome(path);
    };
    long long ns = measure_ns(workload);
    long long bytes_per_second = reference_genome.content.size() * 1000000000 / ns;
    cout << "read_genome() ran at " << bytes_per_second/1000000 << "MB/s in " << ns/1000000 << "ms"<< endl;
}

void read_fastq_performance_test(string path){
    vector<Read_t> reads;
    auto workload = [&path, &reads](){
        reads = read_fastq(path);
    };
    long long ns = measure_ns(workload);
    
    long long reads_bytes = 0;
    for(auto it = reads.begin(); it != reads.end(); it++){
        reads_bytes += it->content.size();
    }
    long long total_bytes = filesystem::file_size(path);

    long long bytes_per_second = total_bytes * 1000000000 / ns;
    cout << "read_fastq() ran at " << bytes_per_second/1000000 << "MB/s in " << ns/1000000 << "ms"<< endl;
}

void read_maf_performance_test(string path){
    vector<CandidateLocation_t> locations;
    auto workload = [&path, &locations](){
        locations = read_maf(path);
    };
    long long ns = measure_ns(workload);

    long long total_bytes = filesystem::file_size(path);
    long long bytes_per_second = total_bytes * 1000000000 / ns;
    cout << "read_maf() ran at " << bytes_per_second/1000000 << "MB/s in " << ns/1000000 << "ms"<< endl;
}

void read_fastq_and_seed_locations_performance_test(string fastq_path, string seeds_path){
    Genome_t dummy_genome;
    vector<Read_t> reads;
    auto workload = [&fastq_path, &seeds_path, &reads, &dummy_genome](){
        read_fastq_and_seed_locations(dummy_genome, fastq_path, seeds_path, reads);
    };
    long long ns = measure_ns(workload);
    
    long long reads_bytes = 0;
    for(auto it = reads.begin(); it != reads.end(); it++){
        reads_bytes += it->content.size();
    }
    long long total_bytes = reads_bytes*4; //approximate total file size

    long long bytes_per_second = total_bytes * 1000000000 / ns;
    cout << "read_fastq_and_seed_locations() ran at " << bytes_per_second/1000000 << "MB/s in " << ns/1000000 << "ms"<< endl;
}

void io_performance_test(string reference_file_path, string reads_file_path, string alignments_file_path){
    if(enable_log) cerr << "Starting IO performance test..." << endl;

    read_file_performance_test(reads_file_path);
    read_genome_performance_test(reference_file_path);
    read_fastq_performance_test(reads_file_path);
    read_maf_performance_test(alignments_file_path);
    read_fastq_and_seed_locations_performance_test(reads_file_path, alignments_file_path);
}

#define TWOBIT_AT(I, CHARARRAY) (((CHARARRAY)[(I)>>2] >> (6 - ((I%4)<<1))) & 0x03)
void print_twobit_as_ascii(long long length, vector<char> twobit){
    char *res = (char *)malloc(length + 1);
    for(int i = 0; i < length; i++){
        char code = TWOBIT_AT(i, twobit);

        if(code == 0x00) res[i] = 'A';
        if(code == 0x01) res[i] = 'C';
        if(code == 0x02) res[i] = 'G';
        if(code == 0x03) res[i] = 'T';
    }
    res[length] = '\0';
    cout << res << endl;
    free(res);
}

void ascii_to_two_bit_correctness_test(){
    vector<string> inputs = {
        "",
        "A",
        "ACGT",
        "ACGTA",
        "AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT",
        "AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTTA"
    };

    vector<vector<char>> correct_results;
    for(int i = 0; i < inputs.size(); i++){
        vector<char> res;
        for(int quad = 0; quad*4 < inputs[i].size(); quad++){
            char twobit = 0x00;
            for(int sub_idx = 0; sub_idx < 4 && quad*4 + sub_idx < inputs[i].size(); sub_idx++){
                char c = inputs[i][quad*4 + sub_idx];
                char code;
                if(c == 'A') code = 0x00;
                if(c == 'C') code = 0x01;
                if(c == 'G') code = 0x02;
                if(c == 'T') code = 0x03;
                twobit |= code << (6 - 2*sub_idx);
            }
            res.push_back(twobit);
        }
        correct_results.push_back(res);
    }

    char **ascii_strings, **twobit_strings;
    long long *string_lengths;
    cudaMallocManaged(&ascii_strings, sizeof(char *) * inputs.size());
    cudaMallocManaged(&twobit_strings, sizeof(char *) * inputs.size());
    cudaMallocManaged(&string_lengths, sizeof(long long) * inputs.size());

    for(int i = 0; i < inputs.size(); i++){
        cudaMallocManaged(ascii_strings + i, sizeof(char)*inputs[i].size());
        cudaMallocManaged(twobit_strings + i, sizeof(char)*(inputs[i].size()+3)/4);
        string_lengths[i] = inputs[i].size();
        for(int j = 0; j < inputs[i].size(); j++){
            ascii_strings[i][j] = inputs[i][j];
        }
    }

    genasm_gpu::ascii_to_twobit_strings<<<32, 32>>>(inputs.size(), string_lengths, ascii_strings, twobit_strings);
    cudaDeviceSynchronize();

    bool success = true;
    for(int i = 0; i < inputs.size(); i++){
        bool equal = true;
        for(int j = 0; j < correct_results[i].size(); j++){
            if(twobit_strings[i][j] != correct_results[i][j]) equal = false;
        }
        if(!equal){
            cout << "FAILED ascii_to_two_bit_correctness_test: produced twobit string" << endl;
            print_twobit_as_ascii(inputs[i].size(), vector<char>(twobit_strings[i], twobit_strings[i] + correct_results[i].size()));
            cout << "instead of" << endl;
            print_twobit_as_ascii(inputs[i].size(), correct_results[i]);
            cout << " for index " << i << endl << endl;
            success = false;
        }
    }
    if(success){
        cout << "PASSED ascii_to_two_bit_correctness_test" << endl;
    }
}

void ascii_to_two_bit_performance_test(string reads_file_path){    
    vector<Read_t> reads = read_fastq(reads_file_path);
    //reads.erase(reads.begin()+1000, reads.end());


    vector<vector<char>> correct_results;
    for(int i = 0; i < reads.size(); i++){
        vector<char> res;
        for(int quad = 0; quad*4 < reads[i].content.size(); quad++){
            char twobit = 0x00;
            for(int sub_idx = 0; sub_idx < 4 && quad*4 + sub_idx < reads[i].content.size(); sub_idx++){
                char c = reads[i].content[quad*4 + sub_idx];
                char code;
                if(c == 'A') code = 0x00;
                if(c == 'C') code = 0x01;
                if(c == 'G') code = 0x02;
                if(c == 'T') code = 0x03;
                twobit |= code << (6 - 2*sub_idx);
            }
            res.push_back(twobit);
        }
        correct_results.push_back(res);
    }

    long long int total_ascii_length = 0;
    long long int total_twobit_length = 0;
    for(int i = 0; i < reads.size(); i++){
        total_ascii_length += reads[i].content.size();
        total_twobit_length += (reads[i].content.size()+3)/4;
    }

    char **ascii_strings, **twobit_strings;
    long long *string_lengths;
    char *ascii_block, *twobit_block;
    cudaMallocManaged(&ascii_strings, sizeof(char *) * reads.size());
    cudaMallocManaged(&twobit_strings, sizeof(char *) * reads.size());
    cudaMallocManaged(&string_lengths, sizeof(long long) * reads.size());
    cudaMallocManaged(&ascii_block, sizeof(char) * total_ascii_length);
    cudaMallocManaged(&twobit_block, sizeof(char) * total_twobit_length);
    char *next_ascii_start = ascii_block;
    char *next_twobit_start = twobit_block;

    for(int i = 0; i < reads.size(); i++){
        ascii_strings[i] = next_ascii_start;
        twobit_strings[i] = next_twobit_start;
        next_ascii_start += reads[i].content.size();
        next_twobit_start += (reads[i].content.size()+3)/4;
        string_lengths[i] = reads[i].content.size();
        for(int j = 0; j < reads[i].content.size(); j++){
            ascii_strings[i][j] = reads[i].content[j];
        }
    }

    auto workload = [&](){
        genasm_gpu::ascii_to_twobit_strings<<<256, 32>>>(reads.size(), string_lengths, ascii_strings, twobit_strings);
        cudaDeviceSynchronize();
    };
    long long ns = measure_ns(workload);
    long long reads_per_second = reads.size() * 1000000000 / ns;
    cout << "ascii_to_two_bit_performance_test" << endl;
    cout << "ran at " << reads_per_second << " reads/second" << endl;

    for(int i = 0; i < reads.size(); i++){
        bool equal = true;
        for(int j = 0; j < correct_results[i].size(); j++){
            if(twobit_strings[i][j] != correct_results[i][j]) equal = false;
        }
        if(!equal){
            cout << "ERROR in ascii_to_two_bit_performance_test: produced twobit string" << endl;
            print_twobit_as_ascii(reads[i].content.size(), vector<char>(twobit_strings[i], twobit_strings[i] + correct_results[i].size()));
            cout << "instead of" << endl;
            print_twobit_as_ascii(reads[i].content.size(), correct_results[i]);
            cout << " for index " << i << endl << endl;
        }
    }
}

void parse_args(int argc, char **argv, string &reference_file, string &reads_file, string &seeds_file, bool &gpu_info_only, bool &verbose, bool &unit_tests, bool &cpu_performance_test){
    //default values
    reference_file = "datasets/human_genome/pacbio-chr1-simulated-m10k-k5_0001.ref";
    reads_file =     "datasets/human_genome/pacbio-chr1-simulated-m10k-k5_0001.fastq";
    seeds_file =     "datasets/human_genome/pacbio-chr1-simulated-m10k-k5_0001.maf";


    gpu_info_only =         OPT_EXISTS == get_cmd_option(argc, argv, "--gpu_info_only");
    verbose =               OPT_EXISTS == get_cmd_option(argc, argv, "--verbose");
    unit_tests =            OPT_EXISTS == get_cmd_option(argc, argv, "--unit_tests");
    cpu_performance_test =  OPT_EXISTS == get_cmd_option(argc, argv, "--cpu_performance_test");

    bool help_and_exit = false;
    help_and_exit |= OPT_INVALID == get_cmd_option(argc, argv, "--reference", reference_file);
    help_and_exit |= OPT_INVALID == get_cmd_option(argc, argv, "--reads", reads_file);
    help_and_exit |= OPT_INVALID == get_cmd_option(argc, argv, "--seeds", seeds_file);
    help_and_exit |= OPT_INVALID == get_cmd_option(argc, argv, "--gpu_info_only");
    help_and_exit |= OPT_INVALID == get_cmd_option(argc, argv, "--verbose");
    help_and_exit |= OPT_INVALID == get_cmd_option(argc, argv, "--unit_tests");
    help_and_exit |= OPT_INVALID == get_cmd_option(argc, argv, "--cpu_performance_test");
    help_and_exit |= OPT_MISSING != get_cmd_option(argc, argv, "--help");
    help_and_exit |= !check_options(argc, argv, {"--reference", "--reads", "--seeds", "--help", "--gpu_info_only", "--verbose", "--unit_tests", "--cpu_performance_test"});

    string help_text =
        "tests[.exe] [options]\n"
        "Options:\n"
        "--reference=[path to reference FASTA] -- overide default reference data for performance test\n"
        "--reads=[path to reads FASTQ]         -- overide default reads data for performance test\n"
        "--seeds=[path to MAF or PAF]          -- overide default seeds data for performance test\n"
        "--gpu_info_only                       -- only print GPU info\n"
        "--verbose                             -- print progress to stderr. Otherwise, only test results are printed\n"
        "--unit_tests                          -- run unit tests (default: disabled)\n"
        "--cpu_performance_test                -- run cpu algorithm performance test (default: gpu)\n"
        "--help                                -- displays this information\n";

    if(help_and_exit){
        cout << help_text << flush;
        exit(0);
    }
}

void print_gpu_info(){
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    cout << nDevices << " visible GPU(s):" << endl;
    for (int i = 0; i < nDevices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        cout << "idx=" << i;
        cout << " name=\"" << prop.name  << "\"";
        cout << " SMs=" << prop.multiProcessorCount;
        cout << " smem=" << (prop.sharedMemPerMultiprocessor >> 10) << "kiB" << endl;
    }
    cout << endl;
}

int main(int argc, char **argv){
    string reference_file, reads_file, seeds_file;
    bool gpu_info_only, verbose, unit_tests, cpu_performance_test;
    parse_args(argc, argv, reference_file, reads_file, seeds_file, gpu_info_only, verbose, unit_tests, cpu_performance_test);
    if(gpu_info_only){
        print_gpu_info();
        exit(0);
    }

    genasm_cpu::enabled_algorithm_log = verbose;
    genasm_gpu::enabled_algorithm_log = verbose;
    enable_log = verbose;

    if(unit_tests){
        print_gpu_info();
        bitvector_tests();
        //io_performance_test(datasets_dir);
        ascii_to_two_bit_correctness_test();
        //ascii_to_two_bit_performance_test(reads_file);
        cpu_algorithm_correctness_test();
        gpu_algorithm_correctness_test();
    }
    else{
        if(cpu_performance_test){
            cpu_algorithm_performance_test(reference_file, reads_file, seeds_file, 32);
        }
        else{
            print_gpu_info();
            gpu_algorithm_performance_test(reference_file, reads_file, seeds_file);
        }
    }

    return 0;
}