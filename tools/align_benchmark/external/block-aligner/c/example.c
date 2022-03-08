#include <stdio.h>
#include <string.h>

#include "block_aligner.h"

void example1(void) {
    // global alignment
    const char* a_str = "AAAAAAAA";
    const char* b_str = "AARAAAA";
    size_t a_len = strlen(a_str);
    size_t b_len = strlen(b_str);
    SizeRange range = {.min = 32, .max = 32};
    Gaps gaps = {.open = -11, .extend = -1};

    PaddedBytes* a = block_new_padded_aa(a_len, range.max);
    PaddedBytes* b = block_new_padded_aa(b_len, range.max);
    block_set_bytes_padded_aa(a, (const uint8_t*)a_str, a_len, range.max);
    block_set_bytes_padded_aa(b, (const uint8_t*)b_str, b_len, range.max);

    BlockHandle block = block_new_aa(a_len, b_len, range.max);
    block_align_aa(block, a, b, &BLOSUM62, gaps, range, 0);
    AlignResult res = block_res_aa(block);

    printf("a: %s\nb: %s\nscore: %d\nidx: (%lu, %lu)\n",
            a_str,
            b_str,
            res.score,
            res.query_idx,
            res.reference_idx);

    block_free_aa(block);
    block_free_padded_aa(a);
    block_free_padded_aa(b);
}

void example2(void) {
    // global alignment with traceback
    const char* a_str = "AAAAAAAA";
    const char* b_str = "AARAAAA";
    size_t a_len = strlen(a_str);
    size_t b_len = strlen(b_str);
    SizeRange range = {.min = 32, .max = 32};
    Gaps gaps = {.open = -11, .extend = -1};

    PaddedBytes* a = block_new_padded_aa(a_len, range.max);
    PaddedBytes* b = block_new_padded_aa(b_len, range.max);
    block_set_bytes_padded_aa(a, (const uint8_t*)a_str, a_len, range.max);
    block_set_bytes_padded_aa(b, (const uint8_t*)b_str, b_len, range.max);

    BlockHandle block = block_new_aa_trace(a_len, b_len, range.max);
    block_align_aa_trace(block, a, b, &BLOSUM62, gaps, range, 0);
    AlignResult res = block_res_aa_trace(block);

    printf("a: %s\nb: %s\nscore: %d\nidx: (%lu, %lu)\n",
            a_str,
            b_str,
            res.score,
            res.query_idx,
            res.reference_idx);

    Cigar* cigar = block_new_cigar(res.query_idx, res.reference_idx);
    block_cigar_aa_trace(block, res.query_idx, res.reference_idx, cigar);
    size_t cigar_len = block_len_cigar(cigar);
    // Note: 'M' signals either a match or mismatch
    char ops_char[] = {' ', 'M', 'I', 'D'};
    for (int i = 0; i < cigar_len; i++) {
        OpLen o = block_get_cigar(cigar, i);
        printf("%lu%c", o.len, ops_char[o.op]);
    }
    printf("\n");

    block_free_cigar(cigar);
    block_free_aa_trace(block);
    block_free_padded_aa(a);
    block_free_padded_aa(b);
}

int main() {
    example1();
    example2();
}
