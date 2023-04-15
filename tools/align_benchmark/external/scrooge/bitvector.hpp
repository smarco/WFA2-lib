/*
 * This header implements a bitvector type of reconfigurable size and name
 * 
 * The file which includes this header should define the BITVECTOR_BITS macro
 * as the number of bits the bitvector should hold.
 *
 * Optionally, it can also define the BITVECTOR macro as the name the bitvector type will have.
 * If BITVECTOR is not defined, it will default to naming the type bitvector.
 *
 * If the shuffle_up function is used, the number of warps should be passed
 * in the WARPS macro.
 *
 * Both BITVECTOR_BITS ans BITVECTOR will be undef'd by this header.
 *
 * ***** EXAMPLE USAGE *****
 * #define BITVECTOR bitvector //type will be named bitvector
 * #define BITVECTOR_BITS 64 //bitvectors will have 64 bit size
 * #include "bitvector.hpp"
 * 
 *
 */

#include <cstdint>
#include <iostream>
#include <iomanip>
#include <sstream>

#ifndef BITVECTOR
#define BITVECTOR bitvector
#endif

#ifndef BITVECTOR_ELEMENT_TYPE
    #if BITVECTOR_BITS <= 8
        #define BITVECTOR_ELEMENT_TYPE uint8_t
        #define BITVECTOR_ELEMENT_BITS 8
    #elif BITVECTOR_BITS <= 16
        #define BITVECTOR_ELEMENT_TYPE uint16_t
        #define BITVECTOR_ELEMENT_BITS 16
    #elif BITVECTOR_BITS <= 32
        #define BITVECTOR_ELEMENT_TYPE uint32_t
        #define BITVECTOR_ELEMENT_BITS 32
    #elif BITVECTOR_BITS <= 64
        #define BITVECTOR_ELEMENT_TYPE uint64_t
        #define BITVECTOR_ELEMENT_BITS 64
    #else
        #define BITVECTOR_ELEMENT_TYPE uint32_t
        #define BITVECTOR_ELEMENT_BITS 32
    #endif
#endif

#define BITVECTOR_ELEMENTS ((BITVECTOR_BITS + BITVECTOR_ELEMENT_BITS - 1) / BITVECTOR_ELEMENT_BITS)
#define LAST_ELEMENT_BITS (BITVECTOR_BITS - BITVECTOR_ELEMENT_BITS * (BITVECTOR_ELEMENTS - 1))
#define ELEMENT_ZEROS ((BITVECTOR_ELEMENT_TYPE)0)
#define ELEMENT_ONES (~ELEMENT_ZEROS)

#ifdef __CUDACC__
    #define NVCC
#endif

#ifdef NVCC
    #define FS __device__ __host__ //when compiled with nvcc, make functions available for host and device code
#else
    #define FS //when compiled with non-nvcc, omit function specifiers
#endif

#ifdef BITVECTOR_NS
    namespace BITVECTOR_NS{
#endif

struct BITVECTOR{
    #if BITVECTOR_ELEMENTS == 1
        BITVECTOR_ELEMENT_TYPE element;
    #else
        BITVECTOR_ELEMENT_TYPE elements[BITVECTOR_ELEMENTS];
    #endif
    static const int bits;
    static const int n_elements;

    FS BITVECTOR operator|(const BITVECTOR &other) const {
        #if BITVECTOR_ELEMENTS == 1
            return {(BITVECTOR_ELEMENT_TYPE)(element | other.element)};
        #else
            BITVECTOR res;
            for(int i = 0; i < BITVECTOR_ELEMENTS; i++){
                res.elements[i] = elements[i] | other.elements[i];
            }
            return res;
        #endif
    };

    FS BITVECTOR operator&(const BITVECTOR &other) const {
        #if BITVECTOR_ELEMENTS == 1
            return {(BITVECTOR_ELEMENT_TYPE)(element & other.element)};
        #else
            BITVECTOR res;
            for(int i = 0; i < BITVECTOR_ELEMENTS; i++){
                res.elements[i] = elements[i] & other.elements[i];
            }
            return res;
        #endif
    };

    FS BITVECTOR operator~() const {
        #if BITVECTOR_ELEMENTS == 1
            return {(BITVECTOR_ELEMENT_TYPE)(~element)};
        #else
            BITVECTOR res;
            for(int i = 0; i < BITVECTOR_ELEMENTS; i++){
                res.elements[i] = ~elements[i];
            }
            return res;
        #endif
    };

    FS BITVECTOR operator<<(int n) const {
        #if BITVECTOR_ELEMENTS == 1
            if(n < BITVECTOR_ELEMENT_BITS){
                return {(BITVECTOR_ELEMENT_TYPE)(element << n)};
            }
            else{
                return {(BITVECTOR_ELEMENT_TYPE)0};;
            }
        #else
            int element_shift = n/BITVECTOR_ELEMENT_BITS;
            int bit_shift = n % BITVECTOR_ELEMENT_BITS;

            if(element_shift >= BITVECTOR_ELEMENTS){
                return BITVECTOR::zeros();
            }

            BITVECTOR res = BITVECTOR::zeros();
            res.elements[element_shift] = elements[0] << bit_shift;
            for(int i = 1; i < BITVECTOR_ELEMENTS - element_shift; i++){
                BITVECTOR_ELEMENT_TYPE upper = elements[i] << bit_shift;
                BITVECTOR_ELEMENT_TYPE lower = (bit_shift==0)? 0 : elements[i-1] >> (BITVECTOR_ELEMENT_BITS - bit_shift);
                res.elements[i+element_shift] = upper | lower;
            }
            return res;
        #endif
    };

    FS bool operator==(BITVECTOR other) const {
        #if BITVECTOR_ELEMENTS == 1
            return element==other.element;
        #else
            for(int i = 0; i < BITVECTOR_ELEMENTS-1; i++){
                if(elements[i] != other.elements[i]) return false;
            }

            #if LAST_ELEMENT_BITS == BITVECTOR_ELEMENT_BITS
                if(elements[BITVECTOR_ELEMENTS-1] != other.elements[BITVECTOR_ELEMENTS-1]) return false;
            #else
                BITVECTOR_ELEMENT_TYPE last_element_mask = ~((~(BITVECTOR_ELEMENT_TYPE)0) << LAST_ELEMENT_BITS);
                BITVECTOR_ELEMENT_TYPE last_element = elements[BITVECTOR_ELEMENTS-1]&last_element_mask;
                BITVECTOR_ELEMENT_TYPE other_last_element = other.elements[BITVECTOR_ELEMENTS-1]&last_element_mask;
                if(last_element != other_last_element) return false;
            #endif
        #endif

        return true;
    }

    FS bool has_one_at(int n) const {
        #if BITVECTOR_ELEMENTS == 1
            return element & ((BITVECTOR_ELEMENT_TYPE)1 << n);
        #else
            int element_offset = n/BITVECTOR_ELEMENT_BITS;
            int bit_offset = n%BITVECTOR_ELEMENT_BITS;
            return elements[element_offset] & ((BITVECTOR_ELEMENT_TYPE)1 << bit_offset);
        #endif
    }

    FS bool has_zero_at(int n) const {
        return !has_one_at(n);
    }

    /*
     * copy bits into the bitvector, starting at start_n (lsb), ending at start_n+31(msb)
     * start_n must be a multiple of 32
     * bitvector must have contained 0s in [start_n,start_n+31] before calling insert_bits
     * this version of insert_bits is threadsafe for threads within a single block
     */
    #ifdef NVCC
    __device__ void insert_bits_threadsafe(int start_n, uint32_t bits){
        #if BITVECTOR_ELEMENTS == 1
            #if BITVECTOR_ELEMENT_BITS <= 32
                element = (BITVECTOR_ELEMENT_TYPE)bits;
            #else
                BITVECTOR_ELEMENT_TYPE insert_bits = ((BITVECTOR_ELEMENT_TYPE)bits) << start_n;
                atomicOr_block((unsigned long long*)&element, (unsigned long long)insert_bits);
            #endif
        #else
            #if BITVECTOR_ELEMENT_BITS == 32
                int element = start_n/BITVECTOR_ELEMENT_BITS;
                elements[element] = bits;
            #elif BITVECTOR_ELEMENT_BITS > 32
                int element = start_n/BITVECTOR_ELEMENT_BITS;
                int bit_offset = start_n%BITVECTOR_ELEMENT_BITS;
                BITVECTOR_ELEMENT_TYPE insert_bits = ((BITVECTOR_ELEMENT_TYPE)bits)<<bit_offset;
                BITVECTOR_ELEMENT_TYPE old = atomicOr_block(elements+element, insert_bits);
            #else //BITVECTOR_ELEMENT_BITS < 32
                int first_element = start_n/BITVECTOR_ELEMENT_BITS;
                int remaining_elements = 32/BITVECTOR_ELEMENT_BITS;
                for(int element_offset = 0; element_offset < remaining_elements; element_offset++){
                    int element = first_element + element_offset;
                    if(element >= BITVECTOR_ELEMENTS) break;
                    BITVECTOR_ELEMENT_TYPE insert_bits = (BITVECTOR_ELEMENT_TYPE)(bits >> (element_offset*BITVECTOR_ELEMENT_BITS));
                    elements[element] = insert_bits;
                }
            #endif
        #endif
    }
    #endif
    /*
     * copy bits into the bitvector, starting at start_n (lsb), ending at start_n+31(msb)
     * start_n must be a multiple of 32
     * bitvector must have contained 0s in [start_n,start_n+31] before calling insert_bits
     */
    FS void insert_bits(int start_n, uint32_t bits){
        #if BITVECTOR_ELEMENTS == 1
            #if BITVECTOR_ELEMENT_BITS <= 32
                element = (BITVECTOR_ELEMENT_TYPE)bits;
            #else
                BITVECTOR_ELEMENT_TYPE insert_bits = ((BITVECTOR_ELEMENT_TYPE)bits) << start_n;
                element |= insert_bits;
            #endif
        #else
            #if BITVECTOR_ELEMENT_BITS == 32
                int element = start_n/BITVECTOR_ELEMENT_BITS;
                elements[element] = bits;
            #elif BITVECTOR_ELEMENT_BITS > 32
                int element = start_n/BITVECTOR_ELEMENT_BITS;
                int bit_offset = start_n%BITVECTOR_ELEMENT_BITS;
                BITVECTOR_ELEMENT_TYPE insert_bits = ((BITVECTOR_ELEMENT_TYPE)bits)<<bit_offset;
                elements[element] |= insert_bits;
            #else //BITVECTOR_ELEMENT_BITS < 32
                int first_element = start_n/BITVECTOR_ELEMENT_BITS;
                int remaining_elements = 32/BITVECTOR_ELEMENT_BITS;
                for(int element_offset = 0; element_offset < remaining_elements; element_offset++){
                    int element = first_element + element_offset;
                    if(element >= BITVECTOR_ELEMENTS) break;
                    BITVECTOR_ELEMENT_TYPE insert_bits = (BITVECTOR_ELEMENT_TYPE)(bits >> (element_offset*BITVECTOR_ELEMENT_BITS));
                    elements[element] = insert_bits;
                }
            #endif
        #endif
    }

    /*
     * extract 32 bits, starting at start_n (lsb), ending at start_n+31 (msb)
     */
    FS uint32_t extract_bits(int start_n){
        #if BITVECTOR_ELEMENTS == 1
            return (uint32_t)(element>>start_n);
        #else
            int first_element = start_n/BITVECTOR_ELEMENT_BITS;
            int first_element_shift = start_n% BITVECTOR_ELEMENT_BITS;

            BITVECTOR_ELEMENT_TYPE e = elements[first_element];
            uint32_t res = (uint32_t)(e>>first_element_shift);

            int next_bit = BITVECTOR_ELEMENT_BITS - first_element_shift;
            for(int element = first_element+1; element < BITVECTOR_ELEMENTS; element++){
                if(next_bit >=32) break;
                res |= ((uint32_t)elements[element])<<next_bit;
                next_bit += BITVECTOR_ELEMENT_BITS;
            }
            return res;
        #endif
    }

    FS void print(){
        char s[BITVECTOR_BITS+1];
        for(int i = 0; i < BITVECTOR_BITS; i++){
            s[BITVECTOR_BITS-1-i] = '0'+has_one_at(i);
        }
        s[BITVECTOR_BITS] = '\0';
        printf("%s\n", s);
    }

    FS static BITVECTOR single_one_at(int n) {
        #if BITVECTOR_ELEMENTS == 1
            return {(BITVECTOR_ELEMENT_TYPE)(((BITVECTOR_ELEMENT_TYPE)1)<<n)};
        #else
            int element_offset = n/BITVECTOR_ELEMENT_BITS;
            int bit_offset = n%BITVECTOR_ELEMENT_BITS;

            BITVECTOR res = zeros();
            res.elements[element_offset] = (BITVECTOR_ELEMENT_TYPE)1<<bit_offset;
            return res;
        #endif
    }

    FS static BITVECTOR single_zero_at(int n) {
        return ~BITVECTOR::single_one_at(n);
    }

    FS static BITVECTOR zeros(){
        return {};
    }

    #ifdef NVCC
    __device__ static BITVECTOR shuffle_up(BITVECTOR b){
        //shuffle up within a warp
        BITVECTOR res;
        #if BITVECTOR_ELEMENTS == 1
            res.element = __shfl_up_sync(0xFFFFFFFF, b.element, 1);
        #else
            for(int i = 0; i < BITVECTOR_ELEMENTS; i++){
                res.elements[i] = __shfl_up_sync(0xFFFFFFFF, b.elements[i], 1);
            }
        #endif

        #if WARPS > 1
            //shuffle across all warp borders
            __shared__ BITVECTOR carry_out[WARPS - 1];
            int warp = threadIdx.x / 32;
            int last_in_warp = threadIdx.x % 32 == 31;
            int first_in_warp = threadIdx.x % 32 == 0;
            if(warp < WARPS-1 && last_in_warp){
                carry_out[warp] = b;
            }
            __syncthreads();
            if(warp > 0 && first_in_warp){
                res = carry_out[warp-1];
            }
            //__syncthreads();
        #endif
        return res;
    }
    #endif

    FS static BITVECTOR ones(){
        return ~zeros();
    }

    static void print_defines(){
        #define STRINGIFY2(X) #X
        #define STRINGIFY(X) STRINGIFY2(X)
        #define STRINGVAR(X) #X " " STRINGIFY(X)
        #define EVALVAR(X) #X " " << X
        std::cout << STRINGVAR(BITVECTOR_ELEMENT_TYPE) << std::endl;
        std::cout << EVALVAR(BITVECTOR_ELEMENT_BITS) << std::endl;
        std::cout << EVALVAR(BITVECTOR_ELEMENTS) << std::endl;
        std::cout << EVALVAR(BITVECTOR_BITS) << std::endl;
        #undef STRINGIFY
        #undef STRINGIFY2
        #undef STRINGVAR
        #undef EVALVAR
    }
};

const int BITVECTOR::bits = BITVECTOR_BITS;
const int BITVECTOR::n_elements = BITVECTOR_ELEMENTS;

std::ostream& operator<<(std::ostream& os, const BITVECTOR x){
    os << std::hex << std::uppercase << std::setfill('0');
    #if LAST_ELEMENT_BITS < BITVECTOR_ELEMENT_BITS
        BITVECTOR_ELEMENT_TYPE last_element_mask = ~(ELEMENT_ONES << LAST_ELEMENT_BITS);
    #else
        BITVECTOR_ELEMENT_TYPE last_element_mask = ELEMENT_ONES;
    #endif

    #if BITVECTOR_ELEMENTS == 1
        os << std::setw((LAST_ELEMENT_BITS+3)/4);
        os << (x.element & last_element_mask) << " ";
    #else

        os << std::setw((LAST_ELEMENT_BITS+3)/4);
        os << (x.elements[BITVECTOR_ELEMENTS-1] & last_element_mask) << " ";

        for(int i = BITVECTOR_ELEMENTS-2; i >= 0; i--){
            os << std::setw(BITVECTOR_ELEMENT_BITS/4);
            os << x.elements[i] << " ";
        }
    #endif

    return os;
}

#ifdef BITVECTOR_NS
    }
    #undef BITVECTOR_NS
#endif

#undef BITVECTOR_BITS
#undef BITVECTOR
#undef BITVECTOR_ELEMENT_TYPE
#undef BITVECTOR_ELEMENT_BITS
#undef BITVECTOR_ELEMENTS
#undef LAST_ELEMENT_BITS
#undef ELEMENT_ZEROS
#undef ELEMENT_ONES
