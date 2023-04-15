#include "util.hpp"

#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits>
#include <stdexcept>
#include "filesystem.hpp"

using namespace std;

void remove_whitespaces(string &str){
    str.erase(
        remove_if(
            str.begin(),
            str.end(),
            ::isspace
        ),
        str.end()
    );
}

string read_file(string file_path){
    long long file_size = filesystem::file_size(file_path);
    string res(file_size, ' ');

    ifstream f(file_path, ios::binary);
    if(!f.is_open()){
        throw runtime_error("could not read file \"" + file_path + "\"");
    }

    f.read(&res[0], file_size);
    f.close();

    return res;
}

vector<Sequence_t> read_fasta(string file_path){
    /*
     * Reads a FASTA file to a Sequence_t
     * file_path: path to FASTA file
     * returns: vector of sequences contained in FASTA
     */
    string raw_file = read_file(file_path);

    long long i = 0;
    vector<Sequence_t> sequences;

    while(i < raw_file.size()){
        //advance to start of next sequence
        for(; i < raw_file.size(); i++){
            char c = raw_file[i];
            if(c == '>'){
                i++;
                break;
            }
        }

        if(i >= raw_file.size()) break; //no more sequences

        Sequence_t sequence;
        //read description
        for(; i < raw_file.size(); i++){
            char c = raw_file[i];
            if(c == '\n' || c == '\r'){
                break;
            }
            sequence.description.push_back(c);
        }

        //read content
        for(; i < raw_file.size(); i++){
            char c = raw_file[i];
            if(c == '\n' || c == '\r' || c == ' '){
                continue;
            }
            if(c == '>'){
                break;
            }
            sequence.content.push_back(c);
        }

        sequences.push_back(sequence);
    }

    return sequences;
}

Genome_t read_genome(string fasta_file_path){
    vector<Sequence_t> chromosomes = read_fasta(fasta_file_path);

    Genome_t genome;
    long long i = 0;
    for(Sequence_t &chromosome : chromosomes){
        genome.chromosome_starts[chromosome.description] = i;
        genome.content += chromosome.content;
        i += chromosome.content.size();
    }

    return genome;
}

vector<Read_t> read_fastq(string file_path){
    /*
     * Reads a FASTQ file to a vector of reads, start_in_reference an reverse_complement are left uninitialized
     * file_path: path to FASTQ file
     * reads: output reads vector
     * returns: true if successful
     */
    
    string raw_file = read_file(file_path);
    vector<Read_t> reads;

    long long i = 0;
    while(i < raw_file.size()){
        Read_t read;
        //search for start of next read
        for(; i < raw_file.size(); i++){
            char c = raw_file[i];
            if(c == '@'){
                i++;
                break;
            }
        }
        if(i >= raw_file.size()) break; //abort if the endof file is already reached, i.e. no read is left

        //read description
        for(; i < raw_file.size(); i++){
            char c = raw_file[i];
            if(c == '\n'){
                i++;
                break;
            }
            if(c == '\r' || c == ' ') continue;
            read.description.push_back(c);
        }

        //read content
        long long line_start = i;
        for(; i < raw_file.size(); i++){
            char c = raw_file[i];
            if(c == '\n' || c == '\r') break;
        }
        read.content = raw_file.substr(line_start, i-line_start);

        reads.push_back(read);
    }

    return reads;
}

string get_line(string &source, long long &pointer){
    long long i;
    char *s = (char *)source.c_str() + pointer;
    for(i = 0; pointer + i < source.size(); i++){
        if(s[i] == '\n'){
            i++;
            break;
        }
    }
    string res = source.substr(pointer, i);
    while(res.size() && (res.back() == '\n' || res.back() == '\r')){
        res.pop_back();
    }

    pointer += i;

    return res;
}

vector<CandidateLocation_t> read_maf(string file_path){
    string raw_file = read_file(file_path);
    vector<CandidateLocation_t> locations;

    long long i = 0;
    string line;
    while(i < raw_file.size()){
        //search for line starting with a
        line = get_line(raw_file, i);
        if(line.size() == 0) continue;
        if(line[0] == 'a'){
            CandidateLocation_t location;
            string ref_alignment, read_alignment;
            while(true){
                line = get_line(raw_file, i);
                if(line.length() == 0) break;
                if(line[0] == '\r') break;
                //cout << line[0] << endl;
                //cout << "\"" << line << "\"" << endl;
                if(line[0] != 's') continue;

                stringstream line_ss(line.substr(1));
                string src, strand, text;
                long long start, size, srcSize;
                line_ss >> src;
                line_ss >> start;
                line_ss >> size;
                line_ss >> strand;
                line_ss >> srcSize;
                line_ss >> text;
                if(src == "ref"){
                    //loacation.start_in_reference = start;
                    location.start_in_chromosome = start;
                    location.chromosome = "ref";
                    ref_alignment = text;
                }
                else{
                    location.read_description = src;
                    location.strand = (strand == "+");
                    location.start_of_aligned_region = start;
                    location.size_of_aligned_region = size;
                    read_alignment = text;
                }
            }

            long long edit_distance = 0;
            for(long long i = 0; i < read_alignment.size(); i++){
                if(read_alignment[i] != ref_alignment[i]){
                    edit_distance++;
                }
            }

            locations.push_back(location);
        }
    }
    return locations;
}

vector<CandidateLocation_t> read_paf(string file_path){
    string raw_file = read_file(file_path);
    vector<CandidateLocation_t> locations;

    long long i = 0;
    string line;
    while(i < raw_file.size()){
        //search for line starting with a
        line = get_line(raw_file, i);
        if(line.size() == 0) continue;

        stringstream line_ss(line);

        CandidateLocation_t location;
        long long qlen, qstart, qend, tlen, tstart, tend, matches, alignment_size;
        string strand;
        
        getline(line_ss, location.read_description, '\t');
        line_ss >> qlen;
        line_ss >> qstart;
        line_ss >> qend;
        line_ss >> strand;
        line_ss.get();
        getline(line_ss, location.chromosome, '\t');
        line_ss >> tlen;
        line_ss >> tstart;
        line_ss >> tend;
        line_ss >> matches;
        line_ss >> alignment_size;

        location.start_in_chromosome = tstart;
        location.strand = (strand == "+");
        location.start_of_aligned_region = qstart;
        location.size_of_aligned_region = qend - qstart;
        //location.baseline_edit_distance = alignment_size - matches;

        locations.push_back(location);
    }

    return locations;
}

bool ends_with(string const &s, string const &ending)
{
    if (ending.size() > s.size()) return false;
    return equal(ending.rbegin(), ending.rend(), s.rbegin());
}

void left_extend_locations(vector<CandidateLocation_t> &locations){
    for(CandidateLocation_t &location : locations){
        location.start_in_chromosome = max(0ll, location.start_in_chromosome - location.start_of_aligned_region);
        location.size_of_aligned_region += location.start_of_aligned_region;
        location.start_of_aligned_region = 0;
    }
}

void get_global_seeds(Genome_t &genome, vector<CandidateLocation_t> &locations){
    for(CandidateLocation_t &location : locations){
        if(genome.chromosome_starts.size() > 1){
            location.start_in_reference = genome.chromosome_starts[location.chromosome] + location.start_in_chromosome;
        }
        else{
            location.start_in_reference = location.start_in_chromosome;
        }
    }
}

vector<Read_t> read_fastq_and_seed_locations(Genome_t &genome, string fastq_file_path, string seed_file_path, vector<Read_t> &reads){
    vector<CandidateLocation_t> locations;

    if(ends_with(seed_file_path, ".paf")){
        locations = read_paf(seed_file_path);
    }
    else if(ends_with(seed_file_path, ".maf")){
        locations = read_maf(seed_file_path);
    }
    else{
        throw invalid_argument("unknown seed file ending\n");
    }
    left_extend_locations(locations);
    get_global_seeds(genome, locations);

    reads = read_fastq(fastq_file_path);

    map<string, size_t> read_desc_to_idx;
    for(size_t i = 0; i < reads.size(); i++){
        read_desc_to_idx[reads[i].description] = i;
    }

    for(CandidateLocation_t &location : locations){
        if(read_desc_to_idx.find(location.read_description) == read_desc_to_idx.end()){
            cerr << "candidate location specified unknown read \"" << location.read_description << "\"" << endl;
            exit(1);
        }

        size_t i = read_desc_to_idx[location.read_description];
        reads[i].locations.push_back(location);
    }

    return reads;
}

bool cigar_char_equals(char c, char d){
        switch(c){
        case 'A':
        case 'a':
        return d == 'A' || d == 'a';

        case 'C':
        case 'c':
        return d == 'C' || d == 'c';

        case 'G':
        case 'g':
        return d == 'G' || d == 'g';

        case 'T':
        case 't':
        return d == 'T' || d == 't';

        default: {
            cerr << "compared invalid character" << endl;
            return false;
        }
    }
}

/*
 * search for the given key in argv
 * returns: OPT_EXITS if the ket was found,
 *          OPT_MISSING if not
 */
int get_cmd_option(int argc, char **argv, string key){
    for(char **s = argv+1; s != argv+argc; s++){
        string arg(*s);
        size_t key_value_split = arg.find("=");
        string arg_key = arg.substr(0, key_value_split);
        if(arg_key == key){
            if(arg_key.size() == arg.size()){
                return OPT_EXISTS;
            }
            else{
                //had more than
                return OPT_INVALID;
            }
        }
    }
    return OPT_MISSING;
}

/*
 * search for the given key=[val] in argv, if present, retrieve val
 * returns: OPT_EXITS if the ket was found,
 *          OPT_MISSING if not
 *          OPT_INVALID if key was present, but had no value
 */
int get_cmd_option(int argc, char **argv, string key, string &value){
    for(char **s = argv+1; s != argv+argc; s++){
        string arg(*s);
        size_t key_value_split = arg.find("=");
        string arg_key = arg.substr(0, key_value_split);
        if(arg_key == key){
            if(key_value_split < arg.size()-1){
                value = arg.substr(key_value_split+1);
                return OPT_EXISTS;
            }
            else{
                return OPT_INVALID;
            }
        }
    }
    return OPT_MISSING;
}

/*
 * ensure all given option keys match one of the given valid_options
 * return true if all option keys are found in valid_options
 * return false otherwise
 */
bool check_options(int argc, char **argv, vector<string> valid_options){
    for(char **s = argv+1; s != argv+argc; s++){
        string arg(*s);
        size_t key_value_split = arg.find("=");
        string arg_key = arg.substr(0, key_value_split);

        auto valid_options_it = find(valid_options.begin(), valid_options.end(), arg_key);
        if(valid_options_it == valid_options.end()){
            return false;
        }
    }
    return true;
}

/*
 * parse the comma separates values to a vector of ints
 * e.g. "1,20,-5" -> vector<int>{1, 20, -5}
 */
vector<int> parse_csv_numbers(string csv){
    stringstream ss(csv);
    vector<int> result;

    while(ss.good())
    {
        string substr;
        getline(ss, substr, ',');
        result.push_back(stoi(substr));
    }

    return result;
}

vector<string> parse_csv_strings(string csv){
    stringstream ss(csv);
    vector<string> result;

    while(ss.good())
    {
        string substr;
        getline(ss, substr, ',');
        result.push_back(substr);
    }

    return result;
}
