#ifndef CLI_KNOBS
/**********************************************************************************************/
/*************************************** [CONFIG KNOBS] ***************************************/
/**********************************************************************************************/

    //algorithm parameters
    #define W 64
    #define K 64
    #define O 33

    //optimization toggles, comment out to disable any given optimization
    #define STORE_ENTRIES_NOT_EDGES
    //#define DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
    #define EARLY_TERMINATION

    #define DEBUG //uncomment to enable asserts in kernel
    #define DEBUG_OUTPUT //uncomment to enable error messages in kernel, requires DEBUG

/**********************************************************************************************/
/************************************* [END CONFIG KNOBS] *************************************/
/**********************************************************************************************/
#else
    #define W CLI_W
    #define K CLI_K
    #define O CLI_O
    #ifdef CLI_STORE_ENTRIES_NOT_EDGES
        #define STORE_ENTRIES_NOT_EDGES
    #endif
    #ifdef CLI_DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
        #define DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
    #endif
    #ifdef CLI_EARLY_TERMINATION
        #define EARLY_TERMINATION 
    #endif
#endif


/* internal macros */

//GPU properties
#define GPU_ID 0
#define THREAD_BLOCKS smCount(GPU_ID)*THREAD_BLOCKS_PER_SM

//warp and thread numbers for a single block
#define THREADS W
#define WARPS ((THREADS+31)/32)
#define ALL_THREADS 0xFFFFFFFF

//maximum number of text and pattern characters to trace back per window
#define TB_LIMIT (W-O)
//number of bits in bitvector and halfbitvector needed for traceback
#define TB_BITS min((size_t)(W-O+1), m)
//number of bits in bitvector not needed for traceback
#define NON_TB_BITS (m - (TB_BITS))
//index into traceback bitvector where 0==MSB, m-1==LSB, corresponding to the pattern indices
#ifdef DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
    #define TB_BIT(J) ((TB_BITS) - 1 - (J))
#else
    #define TB_BIT(J) (m - 1 - (J))
#endif

//size of R
#ifdef STORE_ENTRIES_NOT_EDGES
    #define BITVECTORS_PER_ELEMENT 1
#else
    #define BITVECTORS_PER_ELEMENT 3
#endif
#ifdef DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
    #define COLUMNS (W-O+1)
#else
    #define COLUMNS (W+1)
#endif
#define ROWS (K+1)
#define R_BITVECTORS (COLUMNS * ROWS * BITVECTORS_PER_ELEMENT)

//indexing into R
#ifdef STORE_ENTRIES_NOT_EDGES
    #define IDX(I, D) (COLUMNS*(D) + (I))
#else
    #define MAT 0
    #define INS 1
    #define DEL 2
    #define IDX(I, D, EDIT_TYPE) (3*(COLUMNS*(D) + (I)) + (EDIT_TYPE))
#endif

//codes used in twobit representation
#define A 0x00
#define C 0x01
#define G 0x02
#define T 0x03

#include <limits>
#include <omp.h>
#include <cassert>
#include <cstring>

//define bitvector types
#define BITVECTOR_NS genasm_cpu
#define BITVECTOR bitvector
#define BITVECTOR_BITS W
#include "bitvector.hpp"

#ifdef DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
    #define BITVECTOR_NS genasm_cpu
    #define BITVECTOR halfbitvector
    #if O > 0
        #define BITVECTOR_BITS (TB_LIMIT+1)
    #else
        #define BITVECTOR_BITS W
    #endif
    #include "bitvector.hpp"
#else
    namespace genasm_cpu {
        typedef bitvector halfbitvector;
    }
#endif

#include "genasm_cpu.hpp"
using namespace std;
namespace genasm_cpu {
    bool enabled_algorithm_log = false;

    typedef struct PM {
        bitvector masks[4];
    } PM_t;

    struct str {
        char *base;
        size_t len;
    };

    void print_zero_based(size_t n, char *text){
        char mapping[4];
        mapping[A] = 'A';
        mapping[C] = 'C';
        mapping[G] = 'G';
        mapping[T] = 'T';

        char *subtxt = (char*)malloc(n+1);
        memcpy(subtxt, text, n);
        for(size_t i = 0; i < n; i++){
            subtxt[i] = mapping[text[i]];
        }
        subtxt[n] = '\0';
        printf("%s\n", subtxt);
        free(subtxt);
    }

    void print_state(size_t n, char *text, size_t m, char *pattern, size_t k, halfbitvector *R){
        for(int it = min((size_t)(W-O), n); it >= 0; it--){
            #ifdef STORE_ENTRIES_NOT_EDGES
                printf("text_iteration = %d\n", it);
                for(int dt = 0; dt <= k; dt++){
                    R[IDX(it, dt)].print();
                }
            #else
                printf("text_iteration = %d insertion\n", it);
                for(int dt = 0; dt <= k; dt++){
                    R[IDX(it, dt, INS)].print();
                }
                printf("text_iteration = %d deletion\n", it);
                for(int dt = 0; dt <= k; dt++){
                    R[IDX(it, dt, DEL)].print();
                }
                printf("text_iteration = %d match\n", it);
                for(int dt = 0; dt <= k; dt++){
                    R[IDX(it, dt, MAT)].print();
                }
            #endif
        }

        printf("text   =");
        print_zero_based(n, text);
        printf("pattern=");
        print_zero_based(m, pattern);
    }

    PM_t generatePatternBitmaskACGT(int m, char *pattern){
        PM_t pm;
        pm.masks[A] = bitvector::ones();
        pm.masks[C] = bitvector::ones();
        pm.masks[G] = bitvector::ones();
        pm.masks[T] = bitvector::ones();

        for(int bit_idx = 0; bit_idx < m; bit_idx++){
            int j = m - 1 - bit_idx;
            char curChar = pattern[j];
            pm.masks[curChar] = pm.masks[curChar] & bitvector::single_zero_at(bit_idx);
        }

        //print_zero_based(m, pattern);
        //cout << "A: "; pm.masks[A].print();
        //cout << "C: "; pm.masks[C].print();
        //cout << "G: "; pm.masks[G].print();
        //cout << "T: "; pm.masks[T].print();

        return pm;
    }

    halfbitvector extract_tb_bitvector(size_t m, bitvector b){
        halfbitvector res = halfbitvector::zeros();
        for(int next_bit_offset = 0; next_bit_offset < TB_BITS; next_bit_offset+=32){
            int next_bit = NON_TB_BITS + next_bit_offset;
            uint32_t tmp = b.extract_bits(next_bit);
            res.insert_bits(next_bit_offset, tmp);
        }
        return res;
    }

    long long genasm_dc(size_t n, char *text, size_t m, char *pattern, size_t k, halfbitvector *R, bitvector *forefront){
        PM_t pm = generatePatternBitmaskACGT(m, pattern);
        long long window_edit_distance = std::numeric_limits<long long>::max();
        
        for(long long d = 0; d <= k; d++){
            bitvector center, top, right, topright;
            for(size_t i = n; i != (size_t)-1; i--){
                char curChar = text[i];
                bitvector curPm = pm.masks[curChar];
                bitvector mat, sub, ins, del;

                if(d > 0){
                    top = forefront[i];
                }

                if(d == 0 && i == n){
                    mat = bitvector::ones();
                    sub = bitvector::ones();
                    ins = bitvector::ones();
                    del = bitvector::ones();
                    center = bitvector::ones();
                }
                else if(d == 0){
                    mat = (right << 1) | curPm;
                    sub = bitvector::ones();
                    ins = bitvector::ones();
                    del = bitvector::ones();
                    center = mat;
                }
                else if(i == n){
                    mat = bitvector::ones();
                    sub = bitvector::ones();
                    ins = bitvector::ones()<<d;
                    del = bitvector::ones();
                    center = ins;
                }
                else{
                    mat = (right << 1) | curPm;
                    sub = topright << 1;
                    ins = top << 1;
                    del = topright;
                    center = mat & sub & ins & del;
                }

                right = center;
                topright = top;
                forefront[i] = center;

                #ifdef DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
                    if(i < (W-O+1)){
                        #ifdef STORE_ENTRIES_NOT_EDGES
                            R[IDX(i, d)] = extract_tb_bitvector(m, center);
                        #else
                            R[IDX(i, d, MAT)] = extract_tb_bitvector(m, mat);
                            R[IDX(i, d, INS)] = extract_tb_bitvector(m, ins);
                            R[IDX(i, d, DEL)] = extract_tb_bitvector(m, del);
                        #endif
                    }
                #else
                    #ifdef STORE_ENTRIES_NOT_EDGES
                        R[IDX(i, d)] = center;
                    #else
                        R[IDX(i, d, MAT)] = mat;
                        R[IDX(i, d, INS)] = ins;
                        R[IDX(i, d, DEL)] = del;
                    #endif
                #endif

                if(i==0 && center.has_zero_at(m-1)){
                    window_edit_distance = min(d, window_edit_distance);
                    #ifdef EARLY_TERMINATION
                        return d;
                    #endif
                }
            }
        }

        return window_edit_distance;
    }

    long long genasm_tb(size_t n, char *text, size_t m, char *pattern, long long k, halfbitvector *R, long long window_edit_distance, size_t *text_consumed, size_t *pattern_consumed, char **cigar){
        size_t i = 0;
        size_t j = 0;

        #ifdef DEBUG
            if(window_edit_distance < 0 || window_edit_distance > W){
                #ifdef DEBUG_OUTPUT
                    printf("received invalid window_edit_distance %lld\n", window_edit_distance);
                #endif
                assert(false);
            }
        #endif
        long long d = window_edit_distance;

        char current_edit_type = ' ';
        uint8_t current_edit_count = 0;

        while(j < m){ //trailing deletes are ignored, so we can terminate as soon as pattern is fully matches
            //early termination when windowing
            if(i >= TB_LIMIT) break;
            if(j >= TB_LIMIT) break;

            bool i_limit = i >= n;
            bool d_limit = d == 0;

            bool can_sub, can_ins, can_del;
            #ifdef DEBUG
                bool can_mat;
            #endif
            if(j < m-1){
                #ifdef STORE_ENTRIES_NOT_EDGES
                    can_ins = !d_limit &&             R[IDX(i  , d-1)].has_zero_at(TB_BIT(j+1));
                    can_del = !d_limit && !i_limit && R[IDX(i+1, d-1)].has_zero_at(TB_BIT(j  ));
                    can_sub = !d_limit && !i_limit && R[IDX(i+1, d-1)].has_zero_at(TB_BIT(j+1));
                    #ifdef DEBUG
                        can_mat =         !i_limit && R[IDX(i+1, d  )].has_zero_at(TB_BIT(j+1));
                    #endif
                #else
                    can_ins = R[IDX(i, d, INS)].has_zero_at(TB_BIT(j));
                    can_del = R[IDX(i, d, DEL)].has_zero_at(TB_BIT(j));
                    can_sub = R[IDX(i, d, DEL)].has_zero_at(TB_BIT(j+1));
                    #ifdef DEBUG
                        can_mat = R[IDX(i, d, MAT)].has_zero_at(TB_BIT(j));
                    #endif
                #endif
            }
            else{
                can_ins = !d_limit;
                can_del = false;
                can_sub = !d_limit && !i_limit;
                #ifdef DEBUG
                    can_mat = d==0;
                #endif
            }

            char edit_type;
            if(can_ins){
                j += 1;
                d -= 1;
                edit_type = 'I';
            }
            else if(can_del){
                i += 1;
                d -= 1;
                edit_type = 'D';
            }
            else if(can_sub){
                i += 1;
                j += 1;
                d -= 1;
                edit_type = 'X';
            }
            #ifdef DEBUG
                else if(can_mat){
            #else
                else{
            #endif
                i += 1;
                j += 1;
                edit_type = '=';
            }
            #ifdef DEBUG
                else{
                    #ifdef DEBUG_OUTPUT
                        print_state(n, text, m, pattern, k, R);
                        printf("i = %llu\n", i);
                        printf("j = %llu\n", j);
                        printf("d = %lld\n", d);
                        printf("m = %llu\n", m);
                        printf("n = %llu\n", n);
                        printf("k = %lld\n", k);
                        printf("\n");
                    #endif
                    assert(false);
                }
            #endif

            if(edit_type != current_edit_type){
                if(current_edit_count > 0){
                    int chars_used = sprintf(*cigar, "%d%c", current_edit_count, current_edit_type);
                    *cigar += chars_used;
                }
                current_edit_type = edit_type;
                current_edit_count = 1;
            }
            else{
                current_edit_count++;
            }
        }

        if(current_edit_count > 0){
            int chars_used = sprintf(*cigar, "%d%c", current_edit_count, current_edit_type);
            *cigar += chars_used;
        }

        *text_consumed = i;
        *pattern_consumed = j;
        long long edits_used = window_edit_distance - d;
        return edits_used;
    }

    long long genasm(size_t ref_len, char* ref, size_t read_len, char* read, halfbitvector* R, bitvector* forefront, char* cigar){
        size_t ref_idx = 0;
        size_t read_idx = 0;
        long long edit_distance = 0;

        while(read_idx < read_len){
            size_t n = min((size_t)W, ref_len - ref_idx);
            size_t m = min((size_t)W, read_len - read_idx);
            char *text = ref + ref_idx;
            char *pattern = read + read_idx;

            long long window_edit_distance = genasm_dc(n, text, m, pattern, K, R, forefront);

            size_t text_consumed, pattern_consumed;
            long long edits_used = genasm_tb(n, text, m, pattern, K, R, window_edit_distance, &text_consumed, &pattern_consumed, &cigar);
            /*if(threadIdx.x==0){
                edits_used = 10 + R[IDX(0, K)].has_zero_at(TB_BIT(0));
                text_consumed = min(W-O, n);
                pattern_consumed = min(W-O, m);
            }*/

            edit_distance += edits_used;
            ref_idx += text_consumed;
            read_idx += pattern_consumed;
        }

        return edit_distance;
    }

    long long genasm_kernel(vector<str> &refs, vector<str> &reads, vector<char*> &cigars, vector<long long> &edit_distances, int threads=1){
        halfbitvector *R = (halfbitvector *)malloc(sizeof(halfbitvector)*R_BITVECTORS);
        bitvector *forefront = (bitvector *)malloc(sizeof(bitvector)*(W+1));
        if (R == NULL || forefront == NULL){
            cerr << "failed to allocate working set memory" << endl;
            exit(1);
        }

        long long ns = measure_ns([&](){
          for(long long i = 0; i < refs.size(); i++){
              long long ed = genasm(refs[i].len, refs[i].base, reads[i].len, reads[i].base, R, forefront, cigars[i]);
              edit_distances[i] = ed;
          }
        });

        free(R);
        free(forefront);

        return ns;
    }

    char *ascii_to_zero_based_string(string &ascii){
        char *zero_based_string = (char*)malloc(ascii.size());
        if(zero_based_string == NULL){
            cerr << "failed to allocate zero_based_string" << endl;
            exit(1);
        }
        for(size_t i = 0; i < ascii.size(); i++){
            switch (ascii[i])
            {
            case 'A':
            case 'a':
                zero_based_string[i] = A;
                break;
            case 'C':
            case 'c':
                zero_based_string[i] = C;
                break;
            case 'G':
            case 'g':
                zero_based_string[i] = G;
                break;
            case 'T':
            case 't':
                zero_based_string[i] = T;
                break;
            default:
                assert(false);
                break;
            }
        }
        return zero_based_string;
    }

    std::vector<Alignment_t> align_all(Genome_t &reference, std::vector<Read_t> &reads, int threads, long long* core_algorithm_ns){
        size_t num_pairs = 0;
        for(Read_t &read : reads){
            num_pairs += read.locations.size();
        }
        
        vector<Alignment_t> alignments;

        vector<str> refs, rds;
        vector<char *> cigars;
        vector<long long> edit_distances;

        if(enabled_algorithm_log) cerr << "Preparing data..." << endl;
        char *zero_based_reference = ascii_to_zero_based_string(reference.content);

        for(Read_t &read : reads){
            for(CandidateLocation_t &location : read.locations){
                char* ref_base = zero_based_reference + location.start_in_reference;
                size_t ref_len = reference.content.size() - location.start_in_reference;
                refs.push_back({ref_base, ref_len});

                char* read_base = ascii_to_zero_based_string(read.content);
                size_t read_len = read.content.size();
                rds.push_back({read_base, read_len});

                char *cigar = (char *)malloc(read_len*4+1);
                if(cigar == NULL){
                    cerr << "failed to allocate cigar storage" << endl;
                    exit(1);
                }
                cigar[0] = '\0';
                cigars.push_back(cigar);
                edit_distances.push_back(0);
            }
        }

        if(enabled_algorithm_log) cerr << "Starting genasm_kernel..." << endl;

        long long ns = genasm_kernel(refs, rds, cigars, edit_distances, threads);

        if(core_algorithm_ns != NULL){
            *core_algorithm_ns = ns;
        }
        long long alignments_per_second = num_pairs * 1000000000 / ns;
        if(enabled_algorithm_log) cerr << "core algorithm ran at " << alignments_per_second << " aligns/second" << endl;

        if(enabled_algorithm_log) cerr << "post-processing results..." << endl;

        size_t pair_idx = 0;
        for(Read_t &read : reads){
            for(CandidateLocation_t &location : read.locations){
                alignments.push_back(Alignment_t{string(cigars[pair_idx]), edit_distances[pair_idx]});
                free(rds[pair_idx].base);
                pair_idx++;
            }
        }
        free(zero_based_reference);

        if(enabled_algorithm_log) cerr << "pairs=" << num_pairs << " alignments=" << alignments.size() << endl;
        return alignments;
    }

    std::vector<Alignment_t> align_all(vector<string> &texts, vector<string> &queries, int threads, long long* core_algorithm_ns){
        size_t num_pairs = texts.size();
        assert(num_pairs == queries.size());

        vector<Alignment_t> alignments;

        vector<str> zb_texts, zb_queries;
        vector<char *> cigars;
        vector<long long> edit_distances;

        if(enabled_algorithm_log) cerr << "Preparing data..." << endl;

        for(size_t i = 0; i < num_pairs; i++){
            char* query_base = ascii_to_zero_based_string(queries[i]);
            size_t query_len = queries[i].size();
            zb_queries.push_back({query_base, query_len});
            
            char* text_base = ascii_to_zero_based_string(texts[i]);
            size_t text_len = queries[i].size();
            zb_texts.push_back({text_base, text_len});

            char *cigar = (char *)malloc(query_len*4+1);
            if(cigar == NULL){
                cerr << "failed to allocate cigar storage" << endl;
                exit(1);
            }
            cigar[0] = '\0';
            cigars.push_back(cigar);
            edit_distances.push_back(0);
        }

        if(enabled_algorithm_log) cerr << "Starting genasm_kernel..." << endl;
        long long ns = measure_ns([&](){
            genasm_kernel(zb_texts, zb_queries, cigars, edit_distances, threads);
        });
        if(core_algorithm_ns != NULL){
            *core_algorithm_ns = ns;
        }
        long long alignments_per_second = num_pairs * 1000000000 / ns;
        if(enabled_algorithm_log) cerr << "core algorithm ran at " << alignments_per_second << " aligns/second" << endl;

        if(enabled_algorithm_log) cerr << "post-processing results..." << endl;

        for(size_t pair_idx = 0; pair_idx < num_pairs; pair_idx++){
                alignments.push_back(Alignment_t{string(cigars[pair_idx]), edit_distances[pair_idx]});
                free(zb_texts[pair_idx].base);
                free(zb_queries[pair_idx].base);
                pair_idx++;
        }

        if(enabled_algorithm_log) cerr << "pairs=" << num_pairs << " alignments=" << alignments.size() << endl;
        return alignments;
    }
}
