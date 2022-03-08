
/**
 * @file v2i64.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V2I64_H_INCLUDED
#define _V2I64_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 8bit 32cell */
typedef struct v2i64_s {
	__m128i v1;
} v2i64_t;

/* expanders (without argument) */
#define _e_x_v2i64_1(u)
#define _e_x_v2i64_2(u)

/* expanders (without immediate) */
#define _e_v_v2i64_1(a)				(a).v1
#define _e_v_v2i64_2(a)				(a).v1
#define _e_vv_v2i64_1(a, b)			(a).v1, (b).v1
#define _e_vv_v2i64_2(a, b)			(a).v1, (b).v1
#define _e_vvv_v2i64_1(a, b, c)		(a).v1, (b).v1, (c).v1
#define _e_vvv_v2i64_2(a, b, c)		(a).v1, (b).v1, (c).v1

/* expanders with immediate */
#define _e_i_v2i64_1(imm)			(imm)
#define _e_i_v2i64_2(imm)			(imm)
#define _e_vi_v2i64_1(a, imm)		(a).v1, (imm)
#define _e_vi_v2i64_2(a, imm)		(a).v1, (imm)
#define _e_vvi_v2i64_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v2i64_2(a, b, imm)	(a).v1, (b).v1, (imm)

/* address calculation macros */
#define _addr_v2i64_1(imm)			( (__m128i *)(imm) )
#define _addr_v2i64_2(imm)			( (__m128i *)(imm) )
#define _pv_v2i64(ptr)				( _addr_v2i64_1(ptr) )
/* expanders with pointers */
#define _e_p_v2i64_1(ptr)			_addr_v2i64_1(ptr)
#define _e_p_v2i64_2(ptr)			_addr_v2i64_2(ptr)
#define _e_pv_v2i64_1(ptr, a)		_addr_v2i64_1(ptr), (a).v1
#define _e_pv_v2i64_2(ptr, a)		_addr_v2i64_2(ptr), (a).v1

/* expand intrinsic name */
#define _i_v2i64(intrin) 			_mm_##intrin##_epi64
#define _i_v2i64x(intrin)			_mm_##intrin##_si128

/* apply */
#define _a_v2i64(intrin, expander, ...) ( \
	(v2i64_t) { \
		_i_v2i64(intrin)(expander##_v2i64_1(__VA_ARGS__)) \
	} \
)
#define _a_v2i64x(intrin, expander, ...) ( \
	(v2i64_t) { \
		_i_v2i64x(intrin)(expander##_v2i64_1(__VA_ARGS__)) \
	} \
)
#define _a_v2i64xv(intrin, expander, ...) { \
	_i_v2i64x(intrin)(expander##_v2i64_1(__VA_ARGS__)); \
}

/* load and store */
#define _load_v2i64(...)	_a_v2i64x(load, _e_p, __VA_ARGS__)
#define _loadu_v2i64(...)	_a_v2i64x(loadu, _e_p, __VA_ARGS__)
#define _store_v2i64(...)	_a_v2i64xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v2i64(...)	_a_v2i64xv(storeu, _e_pv, __VA_ARGS__)

/* broadcast */
// #define _set_v2i64(...)		_a_v2i64(set1, _e_i, __VA_ARGS__)
#define _set_v2i64(x)		( (v2i64_t) { _mm_set1_epi64x(x) } )
#define _zero_v2i64()		_a_v2i64x(setzero, _e_x, _unused)
#define _seta_v2i64(x, y)	( (v2i64_t) { _mm_set_epi64x(x, y) } )
#define _swap_v2i64(x) ( \
	(v2i64_t) { \
		_mm_shuffle_epi32((x).v1, 0x4e) \
	} \
)

/* logics */
#define _not_v2i64(...)		_a_v2i64x(not, _e_v, __VA_ARGS__)
#define _and_v2i64(...)		_a_v2i64x(and, _e_vv, __VA_ARGS__)
#define _or_v2i64(...)		_a_v2i64x(or, _e_vv, __VA_ARGS__)
#define _xor_v2i64(...)		_a_v2i64x(xor, _e_vv, __VA_ARGS__)
#define _andn_v2i64(...)	_a_v2i64x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v2i64(...)		_a_v2i64(add, _e_vv, __VA_ARGS__)
#define _sub_v2i64(...)		_a_v2i64(sub, _e_vv, __VA_ARGS__)
// #define _max_v2i64(...)		_a_v2i64(max, _e_vv, __VA_ARGS__)
// #define _min_v2i64(...)		_a_v2i64(min, _e_vv, __VA_ARGS__)
#define _max_v2i64(a, b)	( (v2i64_t) { _mm_max_epi32(a.v1, b.v1) } )
#define _min_v2i64(a, b)	( (v2i64_t) { _mm_min_epi32(a.v1, b.v1) } )

/* shuffle */
// #define _shuf_v2i64(...)	_a_v2i64(shuffle, _e_vv, __VA_ARGS__)

/* blend */
#define _sel_v2i64(mask, a, b) ( \
	(v2i64_t) { \
		_mm_blendv_epi8((b).v1, (a).v1, (mask).v1) \
	} \
)

/* compare */
#define _eq_v2i64(...)		_a_v2i64(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v2i64(...)		_a_v2i64(cmpgt, _e_vv, __VA_ARGS__)

/* test: take mask and test if all zero */
#define _test_v2i64(x, y)	_mm_test_all_zeros((x).v1, (y).v1)

/* insert and extract */
#define _ins_v2i64(a, val, imm) { \
	(a).v1 = _i_v2i64((a).v1, (val), (imm)); \
}
#define _ext_v2i64(a, imm) ( \
	(int64_t)_i_v2i64(extract)((a).v1, (imm)) \
)

/* shift */
#define _shlv_v2i64(a, n) ( \
	(v2i64_t) {_i_v2i64(sll)((a).v1, (n).v1)} \
)
#define _shrv_v2i64(a, n) ( \
	(v2i64_t) {_i_v2i64(srl)((a).v1, (n).v1)} \
)

/* mask */
#define _mask_v2i64(a) ( \
	(uint32_t) (_mm_movemask_epi8((a).v1)) \
)
#define V2I64_MASK_00		( 0x0000 )
#define V2I64_MASK_01		( 0x00ff )
#define V2I64_MASK_10		( 0xff00 )
#define V2I64_MASK_11		( 0xffff )

/* convert */
#define _cvt_v2i32_v2i64(a) ( \
	(v2i64_t) { \
		_mm_cvtepi32_epi64((a).v1) \
	} \
)

/* transpose */
#define _lo_v2i64(a, b) ( \
	(v2i64_t) { \
		_mm_unpacklo_epi64((a).v1, (b).v1) \
	} \
)
#define _hi_v2i64(a, b) ( \
	(v2i64_t) { \
		_mm_unpackhi_epi64((a).v1, (b).v1) \
	} \
)

/* debug print */
// #ifdef _LOG_H_INCLUDED
#define _print_v2i64(a) { \
	debug("(v2i64_t) %s(%lx, %lx)", #a, _ext_v2i64(a, 1), _ext_v2i64(a, 0)); \
}
// #else
// #define _print_v2i64(x)		;
// #endif

#endif /* _V2I64_H_INCLUDED */
/**
 * end of v2i64.h
 */
