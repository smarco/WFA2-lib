
#include <iostream>
#include <iomanip>

#define BITVECTOR_BITS 65
#define BITVECTOR bv
#define BITVECTOR_ELEMENT_TYPE uint32_t
#define BITVECTOR_ELEMENT_BITS 32
#include "bitvector.hpp"

#define BITVECTOR_BITS 32
#define BITVECTOR hbv
#include "bitvector.hpp"

using namespace std;

#define TEST(CONDITION) if(!(CONDITION)) { \
    cout << "FAILED " << __func__ << " at " << __FILE__ << ":" << __LINE__ << endl; \
    passed = false; \
    }

bool shift_test(){
    bool passed = true;

    bv a = {0x89ABCDEF, 0x01234567, 0x1};
    TEST((a <<  0 == bv{0x89ABCDEF, 0x01234567, 0x1}))
    TEST((a <<  4 == bv{0x9ABCDEF0, 0x12345678, 0x0}))
    TEST((a <<  8 == bv{0xABCDEF00, 0x23456789, 0x1}))
    TEST((a << 12 == bv{0xBCDEF000, 0x3456789A, 0x0}))
    TEST((a << 16 == bv{0xCDEF0000, 0x456789AB, 0x1}))
    TEST((a << 20 == bv{0xDEF00000, 0x56789ABC, 0x0}))
    TEST((a << 24 == bv{0xEF000000, 0x6789ABCD, 0x1}))
    TEST((a << 65 == bv::zeros()))
    
    return passed;
}

bool or_test(){
    bool passed = true;

    bv a = {0xFFFF0000, 0xFFFF0000, 0x0};
    bv b = {0x0000FFFF, 0x0000FFFF, 0x1};

    bv c = {0x0F0F0F0F, 0x0F0F0F0F, 0x0};
    bv d = {0xF0F0F0F0, 0xF0F0F0F0, 0x1};

    bv e = {0x0A0A0A0A, 0x0A0A0A0A, 0x0};
    bv f = {0x05050505, 0x05050505, 0x0};

    TEST(((a | b) == bv::ones()))
    TEST(((c | d) == bv::ones()))
    TEST(((e | f) == c))
    TEST(((bv::ones() | bv::zeros()) == bv::ones()))

    return passed;
}

bool and_test(){
    bool passed = true;

    bv a = {0xFFFF0000, 0xFFFF0000, 0x0};
    bv b = {0x55550000, 0x55550000, 0x0};

    bv c = {0x5555BBBB, 0x5555BBBB, 0x1};
    bv d = {0xBBBB5555, 0xBBBB5555, 0x1};
    bv e = {0x11111111, 0x11111111, 0x1};

    TEST(((a & c) == b))
    TEST(((c & d) == e))
    TEST(((bv::ones() & bv::zeros()) == bv::zeros()))

    return passed;
}

bool not_test(){
    bool passed = true;

    bv a = {0xFFFF0000, 0xFFFF0000, 0x0};
    bv b = {0x0000FFFF, 0x0000FFFF, 0x1};

    TEST(~bv::ones() == bv::zeros())
    TEST(~a == b)

    return passed;
}

bool has_one_at_test(){
    bool passed = true;

    for(int i = 0; i < bv::bits; i++){
        TEST(bv::ones().has_one_at(i));
    }

    for(int i = 0; i < bv::bits; i++){
        TEST(!bv::zeros().has_one_at(i));
    }

    bv a = {0x11111111, 0x11111111, 0x1};
    for(int i = 0; i < bv::bits; i++){
        TEST(a.has_one_at(i) == (i%4==0))
    }

    return passed;
}

bool single_one_at_test(){
    bool passed = true;

    TEST(bv::single_one_at(0) == bv{0x1})
    TEST(bv::single_one_at(1) == bv{0x2})
    TEST(bv::single_one_at(2) == bv{0x4})
    TEST((bv::single_one_at(34) == bv{0x0, 0x4}))
    TEST((bv::single_one_at(bv::bits-1) == bv{0x0, 0x0, 0x1}))

    return passed;
}

#define BITVECTOR_BITS 64
#define BITVECTOR bv64
#include "bitvector.hpp"
bool insert_bits_test(){
    bool passed = true;

    bv64 a = bv64::zeros();
    a.insert_bits(32, 0xFF);
    a.insert_bits( 0, 0xAA);
    TEST(a == bv64{0xFF000000AA})

    return passed;
}

bool bitvector_tests(){
    bool passed = true;

    passed &= shift_test();
    passed &= or_test();
    passed &= and_test();
    passed &= not_test();
    passed &= has_one_at_test();
    passed &= single_one_at_test();
    passed &= insert_bits_test();

    if(passed){
        cout << "PASSED bitvector tests" << endl;
    }
    return passed;
}
