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
    #define DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
    //#define EARLY_TERMINATION

    //implementation parameters
    #define THREAD_BLOCKS_PER_SM 20
    #define CIGAR_SUBLIST_SIZE 64

    //uncomment the desired memory type for the DP table R
    #define DP_MEMORY_SHARED
    //#define DP_MEMORY_GLOBAL

    //percent of L1/smem/txcache reserved for shared memory
    //comment out to let CUDA determine it automatically
    //devies of capability 8.6 support 0, 8, 16, 32, 64 or 100kiB per SM
    //devices of capability 8.0 support 0, 8, 16, 32, 64, 100, 132 or 164kiB per SM
    #define SMEM_CARVEOUT_PERCENT 100

    //#define DEBUG //uncomment to enable asserts in kernel
    //#define DEBUG_OUTPUT //uncomment to enable error messages in kernel, requires DEBUG

/**********************************************************************************************/
/************************************* [END CONFIG KNOBS] *************************************/
/**********************************************************************************************/
#else
    #define W CLI_W
    #define K CLI_K
    #define O CLI_O
    #define THREAD_BLOCKS_PER_SM CLI_THREAD_BLOCKS_PER_SM
    #define CIGAR_SUBLIST_SIZE CLI_CIGAR_SUBLIST_SIZE

    #ifdef CLI_STORE_ENTRIES_NOT_EDGES
        #define STORE_ENTRIES_NOT_EDGES
    #endif
    #ifdef CLI_DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
        #define DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
    #endif
    #ifdef CLI_EARLY_TERMINATION
        #define EARLY_TERMINATION 
    #endif

    #ifdef CLI_DP_MEMORY_SHARED
        #define DP_MEMORY_SHARED
    #endif
    #ifdef CLI_DP_MEMORY_GLOBAL
        #define DP_MEMORY_GLOBAL
    #endif

    #ifdef CLI_SMEM_CARVEOUT_PERCENT
        #define SMEM_CARVEOUT_PERCENT CLI_SMEM_CARVEOUT_PERCENT
    #endif
#endif

/* internal macros */

#define GPU_ID 0
#define THREAD_BLOCKS smCount(GPU_ID)*THREAD_BLOCKS_PER_SM

//warp and thread numbers for a single block
#define THREADS W
#define WARPS ((THREADS+31)/32)
#define ALL_THREADS 0xFFFFFFFF

//maximum number of text and pattern characters to trace back per window
#define TB_LIMIT (W-O)
//number of bits in bitvector and halfbitvector needed for traceback
#define TB_BITS min(W-O+1, m)
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
    #define IDX(I, D) (ROWS*(I) + (D))
#else
    #define MAT 0
    #define INS 1
    #define DEL 2
    #define IDX(I, D, EDIT_TYPE) (3*(ROWS*(I) + (D)) + (EDIT_TYPE))
#endif

//codes used in twobit representation
#define A 0x00
#define C 0x01
#define G 0x02
#define T 0x03

#include "util.hpp"
#include "cuda_util.hpp"
#include "genasm_gpu.hpp"

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <assert.h>

using namespace std;

bool genasm_gpu::enabled_algorithm_log = true;

#define CUDA_LIST cigar_list
#define CUDA_LIST_SUBLIST_SIZE CIGAR_SUBLIST_SIZE
#define CUDA_LIST_CONTAINED_TYPE CigarEntry_t
#include "cuda_list.hpp"

#define BITVECTOR_NS genasm_gpu
#define BITVECTOR bitvector
#define BITVECTOR_BITS W
#include "bitvector.hpp"

#ifdef DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
    #define BITVECTOR_NS genasm_gpu
    #define BITVECTOR halfbitvector
    #if O > 0
        #define BITVECTOR_BITS (TB_LIMIT+1)
    #else
        #define BITVECTOR_BITS W
    #endif
    #include "bitvector.hpp"
#else
    namespace genasm_gpu {
        typedef bitvector halfbitvector;
    }
#endif

namespace genasm_gpu {
    typedef struct PM {
        bitvector masks[4];
    } PM_t;

    typedef struct TwoBitArray {
        const char* __restrict__ base;
        unsigned long long offset;  //offset as #ACGT characters
                                    //necessary since *base is only at 4-ACGT granularity
        unsigned long long size; //length as #ACGT characters
    } TwoBitArray_t;

    struct AlignmentRes {
        long long edit_distance;
        cigar_list cigar;
    };

    __device__ char twobit_at(unsigned int i, TwoBitArray_t twobit){
        unsigned int byte_idx = (i+twobit.offset)>>2;
        char subbyte_idx = (i+twobit.offset)%4;

        char quad_code = twobit.base[byte_idx];
        char shifted_quad_code = quad_code >> (6 - (subbyte_idx<<1));
        return shifted_quad_code & 0x3;
    }

    __host__ __device__ TwoBitArray_t twobit_add(TwoBitArray_t twobit, int x){
        TwoBitArray_t res;
        res.base = twobit.base;
        res.offset = twobit.offset + x;
        res.size = twobit.size - x;
        return res;
    }

    __device__ void print_twobit_as_ascii(long long length, TwoBitArray_t twobit){
        char *res = (char *)malloc(length + 1);
        for(int i = 0; i < length; i++){
            char code = twobit_at(i, twobit);
            if(code == A) res[i] = 'A';
            if(code == C) res[i] = 'C';
            if(code == G) res[i] = 'G';
            if(code == T) res[i] = 'T';
        }
        res[length] = '\0';
        printf("%s\n", res);
        free(res);
    }

    __device__ PM_t generatePatternBitmaskACGT(int m, TwoBitArray_t pattern){
        int bit_idx = threadIdx.x;
        int i = m - 1 - bit_idx;
        int warp_bit_idx = bit_idx%32;

        uint32_t warp_masks[4];
        warp_masks[A] = 0xFFFFFFFF;
        warp_masks[C] = 0xFFFFFFFF;
        warp_masks[G] = 0xFFFFFFFF;
        warp_masks[T] = 0xFFFFFFFF;

        if(i >= 0){
            char code = twobit_at(i, pattern);
            warp_masks[code] = ~(1u<<warp_bit_idx);
        }

        for (int offset = 16; offset > 0; offset /= 2){
            int t = threadIdx.x;
            bool is_recipient = t%32<offset;
            bool transmitter_valid = t+offset <blockDim.x;
            bool do_shuffle = transmitter_valid || !is_recipient;
            uint32_t threads_mask = __ballot_sync(ALL_THREADS, do_shuffle);
            if(!do_shuffle)
                continue;

            warp_masks[A] &= __shfl_down_sync(threads_mask, warp_masks[A], offset);
            warp_masks[C] &= __shfl_down_sync(threads_mask, warp_masks[C], offset);
            warp_masks[G] &= __shfl_down_sync(threads_mask, warp_masks[G], offset);
            warp_masks[T] &= __shfl_down_sync(threads_mask, warp_masks[T], offset);
        }

        __shared__ PM_t pm;
        if(threadIdx.x == 0){
            pm.masks[A] = bitvector::zeros();
            pm.masks[C] = bitvector::zeros();
            pm.masks[G] = bitvector::zeros();
            pm.masks[T] = bitvector::zeros();
        }
        __syncthreads();

        if(bit_idx % 32 == 0){
            pm.masks[A].insert_bits_threadsafe(bit_idx, warp_masks[A]);
            pm.masks[C].insert_bits_threadsafe(bit_idx, warp_masks[C]);
            pm.masks[G].insert_bits_threadsafe(bit_idx, warp_masks[G]);
            pm.masks[T].insert_bits_threadsafe(bit_idx, warp_masks[T]);
        }

        __syncthreads();

        return pm;
    }

    __device__ halfbitvector extract_tb_bitvector(int m, bitvector b){
        halfbitvector res = halfbitvector::zeros();
        for(int next_bit_offset = 0; next_bit_offset < TB_BITS; next_bit_offset+=32){
            int next_bit = NON_TB_BITS + next_bit_offset;
            uint32_t tmp = b.extract_bits(next_bit);
            res.insert_bits(next_bit_offset, tmp);
        }
        return res;
    }

    __device__ void print_state(int n, TwoBitArray_t text, int m, TwoBitArray_t pattern, int k, halfbitvector *R){
        for(int it = min((W-O), m); it >= 0; it--){
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
        print_twobit_as_ascii(n, text);
        printf("pattern=");
        print_twobit_as_ascii(m, pattern);
    }

    __device__ int genasm_dc(int n, TwoBitArray_t text, int m, TwoBitArray_t pattern, int k, halfbitvector* __restrict__ R
    #ifdef DEBUG
        , size_t read_number
    #endif
    ){
        PM_t pm = generatePatternBitmaskACGT(m, pattern);

        //if n is within the area TB may reach we need to initialize it as well
        #ifdef STORE_ENTRIES_NOT_EDGES
            #ifdef DISCARD_ENTRIES_NOT_USED_BY_TRACEBACK
                if(n < (W-O+1)){
                    if(threadIdx.x == 0){
                        R[IDX(n, 0)] = halfbitvector::ones();
                    }
                    
                    int d = threadIdx.x+1;
                    R[IDX(n, d)] = halfbitvector::ones() << max(0, d - NON_TB_BITS);
                }
            #else
                {
                if(threadIdx.x == 0){
                    R[IDX(n, 0)] = halfbitvector::ones();
                }
                int d = threadIdx.x+1;
                R[IDX(n, d)] = halfbitvector::ones() << d;
                }
            #endif
        #endif

        int i = n - threadIdx.x - 1;

        char curChar;
        bitvector curPm;
        if(i >= 0){
            curChar = twobit_at(i, text);
            curPm = pm.masks[curChar];
        }

        int window_edit_distance;
        window_edit_distance = m;

        #ifdef EARLY_TERMINATION
            __shared__ bool early_terminate;
            if(threadIdx.x == 0){
                early_terminate = false;
            }
        #endif

        bitvector top, left, topleft, center, diag;
        bitvector mat, sub, ins, del;
        mat = bitvector::ones();
        sub = bitvector::ones();
        ins = bitvector::ones();
        del = bitvector::ones();
        int d = 0;
        for(int cycle = 0; cycle < n+k; cycle++){
            __syncthreads();
            #ifdef EARLY_TERMINATION
                if(early_terminate) break;
            #endif
            diag = bitvector::shuffle_up(center);

            if(cycle < threadIdx.x) continue;
            if(i < 0) continue;
            if(d > k) continue;

            topleft = left;
            if(i == n-1){
                left = bitvector::ones() << d;
            }
            else{
                left = diag;
            }
            top = center;

            if(d == 0){
                mat = (left << 1) | curPm;
                center = mat;
            }
            else{
                del = topleft;
                sub = topleft << 1;
                ins = top << 1;
                mat = (left << 1) | curPm;
                center = del & sub & ins & mat;
            }

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

            if(center.has_zero_at(m-1)){
                window_edit_distance = min(d, window_edit_distance);
                #ifdef EARLY_TERMINATION
                    if(i == 0){
                        early_terminate = true;
                    }
                #endif
            }

            d++;
        }
        __shared__ int shared_window_edit_distance;
        if(i == 0){
            shared_window_edit_distance = window_edit_distance;
        }
        __syncthreads();
        return shared_window_edit_distance;
    }

    __device__ void genasm_tb(int n, TwoBitArray_t text, int m, TwoBitArray_t pattern, int k, halfbitvector* __restrict__ R, int window_edit_distance, int *text_consumed, int *pattern_consumed, int *edits_used, cigar_list *cigar_list
    #ifdef DEBUG
        , size_t read_number
    #endif
    ){
        if(threadIdx.x != 0) return;

        int i = 0;
        int j = 0;
        int d = window_edit_distance;

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
                else if(threadIdx.x == 0){
                    #ifdef DEBUG_OUTPUT
                        printf("genasm_tb got stuck in dead-end at read number %lld\n", read_number);
                        printf("i = %d\n", i);
                        printf("j = %d\n", j);
                        printf("d = %d\n", d);
                        printf("m = %d\n", m);
                        printf("n = %d\n", n);
                        printf("k = %d\n", k);
                        printf("\n");
                        print_state(n, text, m, pattern, k, R);
                    #endif
                    assert(false);
                }
            #endif

            if(edit_type != current_edit_type){
                if(current_edit_count > 0){
                    cigar_list->pushBack({current_edit_count, current_edit_type});
                }
                current_edit_type = edit_type;
                current_edit_count = 1;
            }
            else{
                current_edit_count++;
            }
        }

        if(current_edit_count > 0){
            cigar_list->pushBack({current_edit_count, current_edit_type});
        }
        *text_consumed = i;
        *pattern_consumed = j;
        *edits_used = window_edit_distance - d;
    }

    __device__ AlignmentRes genasm(TwoBitArray_t reference, TwoBitArray_t read, halfbitvector *R
    #ifdef DEBUG
        , size_t read_number
    #endif
    ){
        __shared__ AlignmentRes res;
        if(threadIdx.x == 0){
            res.cigar.init();
            res.edit_distance = 0;
        }

        size_t reference_idx = 0;
        size_t read_idx = 0;

        while(read_idx < read.size){
            int n = min((unsigned long long)W, reference.size - reference_idx);
            int m = min((unsigned long long)W, read.size - read_idx);
            TwoBitArray_t text = twobit_add(reference, reference_idx);
            TwoBitArray_t pattern = twobit_add(read, read_idx);

            #ifdef DEBUG
                int window_edit_distance = genasm_dc(n, text, m, pattern, K, R, read_number);
            #else
                int window_edit_distance = genasm_dc(n, text, m, pattern, K, R);
            #endif

            __shared__ int text_consumed, pattern_consumed, edits_used;
            #ifdef DEBUG
                genasm_tb(n, text, m, pattern, K, R, window_edit_distance, &text_consumed, &pattern_consumed, &edits_used, &res.cigar, read_number);
            #else
                genasm_tb(n, text, m, pattern, K, R, window_edit_distance, &text_consumed, &pattern_consumed, &edits_used, &res.cigar);
            #endif
            __syncthreads();

            if(threadIdx.x == 0){
                res.edit_distance += edits_used;
            }
            reference_idx += text_consumed;
            read_idx += pattern_consumed;
        }

        __syncthreads();
        return res;
    }

    __managed__ size_t next_pair_index;
    __global__ void genasm_kernel(size_t pairs_count, TwoBitArray_t *references, TwoBitArray_t *reads, AlignmentRes *results){ 
        #ifdef DP_MEMORY_SHARED
            extern __shared__ halfbitvector R[];
        #endif
        #ifdef DP_MEMORY_GLOBAL
            __shared__ halfbitvector *R;
            if(threadIdx.x == 0){
                R = (halfbitvector*)malloc(sizeof(halfbitvector)* R_BITVECTORS);
                #ifdef DEBUG
                    if(R == NULL){
                        #ifdef DEBUG_OUTPUT
                            printf("failed to allocate R for %d bitvectors of %lld bytes each (%lld total)\n", R_BITVECTORS, sizeof(halfbitvector), sizeof(halfbitvector)* R_BITVECTORS);
                        #endif
                        assert(false);
                    }
                #endif
            }
        #endif
        
        __shared__ size_t i;
        if(threadIdx.x == 0){
            i = atomicAdd((unsigned long long *)&next_pair_index, 1);
        }
        __syncthreads();

        while(i < pairs_count){
            TwoBitArray_t ref = references[i];
            TwoBitArray_t read = reads[i];

            #ifdef DEBUG
                results[i] = genasm(ref, read, R, i);
            #else
                results[i] = genasm(ref, read, R);
            #endif

            if(threadIdx.x == 0){
                i = atomicAdd((unsigned long long *)&next_pair_index, 1);
            }
            __syncthreads();
        }

        #ifdef DP_MEMORY_GLOBAL
            if(threadIdx.x==0){
                free(R);
            }
        #endif
    }

    __device__ char ascii_to_twobit_code(char ascii){
        if(ascii == 'a' || ascii == 'A') return 0x00;
        if(ascii == 'c' || ascii == 'C') return 0x01;
        if(ascii == 'g' || ascii == 'G') return 0x02;
        if(ascii == 't' || ascii == 'T') return 0x03;
        assert(false); //invalid character
        return 0xFF;
    }

    __device__ void ascii_to_twobit_string(long long length, char *ascii, char* twobit){
        long long quad_offset = threadIdx.x;
        int quad_stride = blockDim.x;

        long long quad_idx = quad_offset;
        for(quad_idx = quad_offset; quad_idx*4+3 < length; quad_idx += quad_stride){
            long long i = quad_idx * 4;

            char c0 = ascii[i+0];
            char c1 = ascii[i+1];
            char c2 = ascii[i+2];
            char c3 = ascii[i+3];
            
            char b0 = ascii_to_twobit_code(c0);        
            char b1 = ascii_to_twobit_code(c1);        
            char b2 = ascii_to_twobit_code(c2);        
            char b3 = ascii_to_twobit_code(c3);

            char b01 = (b0 << 2) | b1;
            char b23 = (b2 << 2) | b3;
            char res = (b01 << 4) | b23;

            twobit[quad_idx] = res;
        }

        if(quad_idx*4 < length){
            char res = 0x00;
            for(long long i = 0; quad_idx*4 + i < length; i++){
                char c = ascii[quad_idx*4 + i];
                res |= ascii_to_twobit_code(c) << (6 - 2*i);
            }
            twobit[quad_idx] = res;
        }
    }

    __global__ void ascii_to_twobit_strings(int count, long long *string_lengths, char **ascii_strings, char **twobit_strings){
        int offset = blockIdx.x;
        int stride = gridDim.x;
        for(int i = offset; i < count; i+=stride){
            ascii_to_twobit_string(string_lengths[i], ascii_strings[i], twobit_strings[i]);
        }
    }

    __global__ void single_ascii_to_twobit_string(long long length, char *ascii, char *twobit){
        ascii_to_twobit_string(length, ascii, twobit);
    }

    /*
     * Given a reference genome
     * convert it to twobit representation and provide list of pointers, one for each candidate location in reads
     * out_twobit_blob is the allocated storage, should later be freed with cudaFree()
     */
    void twobit_reference(Genome_t &reference, vector<Read_t> &reads, char **out_twobit_blob, TwoBitArray_t **out_twobit_arrays){
        size_t ref_len = reference.content.size();
        size_t twobit_ref_len = (ref_len+3)/4;

        size_t total_pairs = 0;
        for(size_t i = 0; i < reads.size(); i++){
            total_pairs += reads[i].locations.size();
        }

        char *reference_cudamem,
            *twobit_blob;
        TwoBitArray_t *twobit_arrays;

        CUDACHK(cudaMallocManaged(&twobit_arrays, sizeof(TwoBitArray_t) * total_pairs));
        CUDACHK(cudaMallocManaged(&reference_cudamem, ref_len));
        CUDACHK(cudaMallocManaged(&twobit_blob, twobit_ref_len));

        memcpy(reference_cudamem, reference.content.c_str(), ref_len);

        single_ascii_to_twobit_string<<<THREAD_BLOCKS, 32>>>(ref_len, reference_cudamem, twobit_blob);
        CUDACHK(cudaPeekAtLastError());
        CUDACHK(cudaDeviceSynchronize());

        // V here
        size_t pair_idx = 0;
        for(size_t i = 0; i < reads.size(); i++){
            for(CandidateLocation_t &location : reads[i].locations){
                TwoBitArray_t ref;
                ref.base = twobit_blob;
                ref.offset = location.start_in_reference;
                ref.size = reference.content.size() - location.start_in_reference;
                twobit_arrays[pair_idx] = ref;
                pair_idx++;
            }
        }
        // A 

        *out_twobit_blob = twobit_blob;
        *out_twobit_arrays = twobit_arrays;

        CUDACHK(cudaFree(reference_cudamem));

        //hint reads data as "readMostly"
        CUDACHK(cudaMemAdvise(twobit_blob, twobit_ref_len, cudaMemAdviseSetReadMostly, GPU_ID));
        CUDACHK(cudaMemAdvise(twobit_arrays, sizeof(TwoBitArray_t) * total_pairs, cudaMemAdviseSetReadMostly, GPU_ID));

        //explicitly prefetch reads data to GPU if possible on the current device
        if(canPrefetch(GPU_ID)){
            CUDACHK(cudaMemPrefetchAsync(twobit_blob, twobit_ref_len, GPU_ID));
            CUDACHK(cudaMemPrefetchAsync(twobit_arrays, sizeof(TwoBitArray_t) * total_pairs, GPU_ID));
        }
    }

    /*
     * Given a list of reads with candidate locations,
     * convert them to twobit representation and provide list of pointers, one for each candidate location
     * out_twobit_blob is the allocated storage, should later be freed with cudaFree()
     */
    void twobit_reads(vector<Read_t> &reads, char **out_twobit_blob, TwoBitArray_t **out_twobit_arrays){
        char *reads_cudamem,
            *twobit_blob;
        TwoBitArray_t *twobit_arrays;

        //count total reads characters and number of alignment pairs
        size_t total_reads_size = 0;
        size_t total_pairs = 0;
        for(size_t i = 0; i < reads.size(); i++){
            total_reads_size += reads[i].content.size();
            total_pairs += reads[i].locations.size();
        }
        size_t total_reads_size_twobit = (total_reads_size+3)/4;

        //allocate input and output for conversion kernel
        CUDACHK(cudaMallocManaged(&twobit_arrays, sizeof(TwoBitArray_t) * total_pairs));
        CUDACHK(cudaMallocManaged(&reads_cudamem, total_reads_size));
        CUDACHK(cudaMallocManaged(&twobit_blob, total_reads_size_twobit));

        //concatenate reads into kernel input
        size_t next_read_offset = 0;
        for(size_t i = 0; i < reads.size(); i++){
            memcpy(reads_cudamem + next_read_offset, reads[i].content.c_str(), reads[i].content.size());
            next_read_offset += reads[i].content.size();
        }

        //convert from byte-per-bp to twobit-per-bp
        single_ascii_to_twobit_string<<<THREAD_BLOCKS, 32>>>(total_reads_size, reads_cudamem, twobit_blob);
        CUDACHK(cudaPeekAtLastError());
        CUDACHK(cudaDeviceSynchronize());

        //initialize pointers into twobitarray, one for each alignment pair
        size_t pair_idx = 0;
        next_read_offset = 0;
        for(size_t i = 0; i < reads.size(); i++){
            TwoBitArray_t read;
            read.base = twobit_blob;
            read.offset = next_read_offset;
            read.size = reads[i].content.size();

            for(CandidateLocation_t &location : reads[i].locations){
                twobit_arrays[pair_idx] = read;
                pair_idx++;
            }

            next_read_offset += reads[i].content.size();
        }

        //write to out parameters
        *out_twobit_blob = twobit_blob;
        *out_twobit_arrays = twobit_arrays;

        //no longer need the input to conversion kernel
        CUDACHK(cudaFree(reads_cudamem));

        //hint reads data as "readMostly"
        CUDACHK(cudaMemAdvise(twobit_blob, total_reads_size_twobit, cudaMemAdviseSetReadMostly, GPU_ID));
        CUDACHK(cudaMemAdvise(twobit_arrays, sizeof(TwoBitArray_t) * total_pairs, cudaMemAdviseSetReadMostly, GPU_ID));

        //explicitly prefetch reads data to GPU if possible on the current device
        if(canPrefetch(GPU_ID)){
            CUDACHK(cudaMemPrefetchAsync(twobit_blob, total_reads_size_twobit, GPU_ID));
            CUDACHK(cudaMemPrefetchAsync(twobit_arrays, sizeof(TwoBitArray_t) * total_pairs, GPU_ID));
        }
    }

    /*
     * Given a list of strings,
     * convert them to twobit representation and provide list of pointers, one for each string
     * out_twobit_blob is the allocated storage, should later be freed with cudaFree()
     */
    void twobit_strings(vector<string> strings, char **out_twobit_blob, TwoBitArray_t **out_twobit_arrays){
        char *strings_cudamem,
            *twobit_blob;
        TwoBitArray_t *twobit_arrays;

        //count total reads characters and number of alignment pairs
        size_t total_strings_size = 0;
        for(size_t i = 0; i < strings.size(); i++){
            total_strings_size += strings[i].size();
        }
        size_t total_strings_size_twobit = (total_strings_size+3)/4;

        //allocate input and output for conversion kernel
        CUDACHK(cudaMallocManaged(&twobit_arrays, sizeof(TwoBitArray_t) * strings.size()));
        CUDACHK(cudaMallocManaged(&strings_cudamem, total_strings_size));
        CUDACHK(cudaMallocManaged(&twobit_blob, total_strings_size_twobit));

        //concatenate strings into kernel input
        size_t next_string_offset = 0;
        for(size_t i = 0; i < strings.size(); i++){
            memcpy(strings_cudamem + next_string_offset, strings[i].c_str(), strings[i].size());
            next_string_offset += strings[i].size();
        }

        //convert from byte-per-bp to twobit-per-bp
        single_ascii_to_twobit_string<<<THREAD_BLOCKS, 32>>>(total_strings_size, strings_cudamem, twobit_blob);
        CUDACHK(cudaPeekAtLastError());
        CUDACHK(cudaDeviceSynchronize());

        //initialize pointers into twobitarray, one for each alignment pair
        next_string_offset = 0;
        for(size_t i = 0; i < strings.size(); i++){
            TwoBitArray_t tba;
            tba.base = twobit_blob;
            tba.offset = next_string_offset;
            tba.size = strings[i].size();

            twobit_arrays[i] = tba;

            next_string_offset += strings[i].size();
        }

        //write to out parameters
        *out_twobit_blob = twobit_blob;
        *out_twobit_arrays = twobit_arrays;

        //no longer need the input to conversion kernel
        CUDACHK(cudaFree(strings_cudamem));

        //hint reads data as "readMostly"
        CUDACHK(cudaMemAdvise(twobit_blob, total_strings_size_twobit, cudaMemAdviseSetReadMostly, GPU_ID));
        CUDACHK(cudaMemAdvise(twobit_arrays, sizeof(TwoBitArray_t) * strings.size(), cudaMemAdviseSetReadMostly, GPU_ID));

        //explicitly prefetch reads data to GPU if possible on the current device
        if(canPrefetch(GPU_ID)){
            CUDACHK(cudaMemPrefetchAsync(twobit_blob, total_strings_size_twobit, GPU_ID));
            CUDACHK(cudaMemPrefetchAsync(twobit_arrays, sizeof(TwoBitArray_t) * strings.size(), GPU_ID));
        }
    }

    string cigarlist_to_cigar(cigar_list cl){
        stringstream cigar_ss;
        for(cigar_list_iterator it = cl.begin(); it != cl.end(); it++){
            cigar_ss << (int)it->edit_count;
            cigar_ss << it->edit_type;
        }
        return cigar_ss.str();
    }

    vector<Alignment_t> align_all(Genome_t &reference, vector<Read_t> &reads, long long* core_algorithm_ns){
        long long ref_len = reference.content.size();
        size_t num_pairs = 0;
        for(Read_t &read : reads){
            num_pairs += read.locations.size();
        }

        char *reference_twobit_blob,
             *reads_twobit_blob;
        TwoBitArray_t *reference_tbas,
                      *read_tbas;
        AlignmentRes *results;

        CUDACHK(cudaMallocManaged(&results, num_pairs * sizeof(AlignmentRes)));

        //reserve the amount of memory at most needed for output cigars
        size_t total_cigar_sublists = 0;
        for(Read_t &read : reads){
            size_t max_edits = read.content.size()*2;
            size_t max_cigar_sublists = (max_edits+CIGAR_SUBLIST_SIZE-1)/CIGAR_SUBLIST_SIZE;
            total_cigar_sublists += max((size_t)1, max_cigar_sublists) * read.locations.size();
        }
        cigar_list::backingStorageInit(total_cigar_sublists);

        twobit_reads(reads, &reads_twobit_blob, &read_tbas);
        twobit_reference(reference, reads, &reference_twobit_blob, &reference_tbas);
        
        if(enabled_algorithm_log) cerr << "Starting Kernel..." << endl;

        size_t dp_mem_per_block = R_BITVECTORS * sizeof(halfbitvector);
        if(enabled_algorithm_log) cerr << "using " << dp_mem_per_block << "B DP memory per thread block" << endl;
        #ifdef DP_MEMORY_GLOBAL
            int safety_factor = 2; //CUDA crashes if we request the exact number of bytes
            setMallocHeapLimit(R_BITVECTORS * sizeof(halfbitvector) * THREAD_BLOCKS * safety_factor);
        #endif
        #ifdef SMEM_CARVEOUT_PERCENT
            smemCarveout(SMEM_CARVEOUT_PERCENT, (void*)genasm_kernel);
        #endif
        #ifdef DP_MEMORY_SHARED
            int smem_limit = maximizeDynamicSmem((void*)genasm_kernel, GPU_ID);
            int dp_smem = R_BITVECTORS * sizeof(halfbitvector);
            if(dp_smem > smem_limit){
                cout << "R requires " << dp_smem << "B, more than the device limit " << smem_limit << "B" << endl;
                exit(1);
            }
        #else
            int dp_smem = 0;
        #endif

        next_pair_index = 0;
        long long ns = measure_ns([&](){
            genasm_kernel<<<THREAD_BLOCKS, THREADS, dp_smem>>>(num_pairs, reference_tbas, read_tbas, results);
            CUDACHK(cudaPeekAtLastError());
            CUDACHK(cudaDeviceSynchronize());
        });

        if(core_algorithm_ns != NULL){
            *core_algorithm_ns = ns;
        }

        long long alignments_per_second = num_pairs * 1000000000 / ns;
        if(enabled_algorithm_log) cerr << "core algorithm ran at " << alignments_per_second << " aligns/second" << endl;

        if(enabled_algorithm_log) cerr << "Post Processing Results..." << endl;

        if(canPrefetch(cudaCpuDeviceId)){
            CUDACHK(cudaMemPrefetchAsync(results, num_pairs*sizeof(AlignmentRes), cudaCpuDeviceId));
            cigar_list::backingStoragePrefetch(cudaCpuDeviceId);
        }

        vector<Alignment_t> alignments(num_pairs);
        size_t pair_idx = 0;
        for(Read_t &read : reads){
            for(CandidateLocation_t &location : read.locations){
                alignments[pair_idx].edit_distance = results[pair_idx].edit_distance;
                alignments[pair_idx].cigar = cigarlist_to_cigar(results[pair_idx].cigar);
                pair_idx++;
            }
        }

        CUDACHK(cudaFree(reference_twobit_blob));
        CUDACHK(cudaFree(reads_twobit_blob));
        CUDACHK(cudaFree(reference_tbas));
        CUDACHK(cudaFree(read_tbas));
        CUDACHK(cudaFree(results));

        cigar_list::backingStorageDestruct();

        if(enabled_algorithm_log) cerr << "pairs=" << num_pairs << " alignments=" << alignments.size() << endl;
        return alignments;
    }

    vector<Alignment_t> align_all(vector<string> &texts, vector<string> &queries, long long* core_algorithm_ns){
        size_t num_pairs = texts.size();
        assert(queries.size() == num_pairs);

        char *texts_twobit_blob,
             *queries_twobit_blob;
        TwoBitArray_t *texts_tbas,
                      *queries_tbas;
        AlignmentRes *results;

        CUDACHK(cudaMallocManaged(&results, num_pairs * sizeof(AlignmentRes)));

        //reserve the amount of memory at most needed for output cigars
        size_t total_cigar_sublists = 0;
        for(string &query : queries){
            size_t max_edits = query.size()*2;
            size_t max_cigar_sublists = (max_edits+CIGAR_SUBLIST_SIZE-1)/CIGAR_SUBLIST_SIZE;
            total_cigar_sublists += max((size_t)1, max_cigar_sublists);
        }
        cigar_list::backingStorageInit(total_cigar_sublists);

        twobit_strings(queries, &queries_twobit_blob, &queries_tbas);
        twobit_strings(texts, &texts_twobit_blob, &texts_tbas);
        
        if(enabled_algorithm_log) cerr << "Starting Kernel..." << endl;

        size_t dp_mem_per_block = R_BITVECTORS * sizeof(halfbitvector);
        if(enabled_algorithm_log) cerr << "using " << dp_mem_per_block << "B DP memory per thread block" << endl;
        #ifdef DP_MEMORY_GLOBAL
            int safety_factor = 2; //CUDA crashes if we request the exact number of bytes
            setMallocHeapLimit(R_BITVECTORS * sizeof(halfbitvector) * THREAD_BLOCKS * safety_factor);
        #endif
        #ifdef SMEM_CARVEOUT_PERCENT
            smemCarveout(SMEM_CARVEOUT_PERCENT, (void*)genasm_kernel);
        #endif
        #ifdef DP_MEMORY_SHARED
            int smem_limit = maximizeDynamicSmem((void*)genasm_kernel, GPU_ID);
            int dp_smem = R_BITVECTORS * sizeof(halfbitvector);
            if(dp_smem > smem_limit){
                cout << "R requires " << dp_smem << "B, more than the device limit " << smem_limit << "B" << endl;
                exit(1);
            }
        #else
            int dp_smem = 0;
        #endif

        next_pair_index = 0;
        long long ns = measure_ns([&](){
            genasm_kernel<<<THREAD_BLOCKS, THREADS, dp_smem>>>(num_pairs, texts_tbas, queries_tbas, results);
            CUDACHK(cudaPeekAtLastError());
            CUDACHK(cudaDeviceSynchronize());
        });
        
        if(core_algorithm_ns != NULL){
            *core_algorithm_ns = ns;
        }

        long long alignments_per_second = num_pairs * 1000000000 / ns;
        if(enabled_algorithm_log) cerr << "core algorithm ran at " << alignments_per_second << " aligns/second" << endl;

        if(enabled_algorithm_log) cerr << "Post Processing Results..." << endl;

        if(canPrefetch(cudaCpuDeviceId)){
            CUDACHK(cudaMemPrefetchAsync(results, num_pairs*sizeof(AlignmentRes), cudaCpuDeviceId));
            cigar_list::backingStoragePrefetch(cudaCpuDeviceId);
        }

        vector<Alignment_t> alignments(num_pairs);
        for(size_t pair_idx = 0; pair_idx < num_pairs; pair_idx++){
            alignments[pair_idx].edit_distance = results[pair_idx].edit_distance;
            alignments[pair_idx].cigar = cigarlist_to_cigar(results[pair_idx].cigar);
        }

        CUDACHK(cudaFree(texts_twobit_blob));
        CUDACHK(cudaFree(queries_twobit_blob));
        CUDACHK(cudaFree(texts_tbas));
        CUDACHK(cudaFree(queries_tbas));
        CUDACHK(cudaFree(results));

        cigar_list::backingStorageDestruct();

        if(enabled_algorithm_log) cerr << "pairs=" << num_pairs << " alignments=" << alignments.size() << endl;
        return alignments;
    }
}
