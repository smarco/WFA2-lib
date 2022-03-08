
/**
 * @file v2i32.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V2I32_H_INCLUDED
#define _V2I32_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 8bit 32cell */
typedef struct v2i32_s {
	__m128i v1;
} v2i32_t;

/* expanders (without argument) */
#define _e_x_v2i32_1(u)
#define _e_x_v2i32_2(u)

/* expanders (without immediate) */
#define _e_v_v2i32_1(a)				(a).v1
#define _e_v_v2i32_2(a)				(a).v1
#define _e_vv_v2i32_1(a, b)			(a).v1, (b).v1
#define _e_vv_v2i32_2(a, b)			(a).v1, (b).v1
#define _e_vvv_v2i32_1(a, b, c)		(a).v1, (b).v1, (c).v1
#define _e_vvv_v2i32_2(a, b, c)		(a).v1, (b).v1, (c).v1

/* expanders with immediate */
#define _e_i_v2i32_1(imm)			(imm)
#define _e_i_v2i32_2(imm)			(imm)
#define _e_vi_v2i32_1(a, imm)		(a).v1, (imm)
#define _e_vi_v2i32_2(a, imm)		(a).v1, (imm)
#define _e_vvi_v2i32_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v2i32_2(a, b, imm)	(a).v1, (b).v1, (imm)

/* address calculation macros */
#define _addr_v2i32_1(imm)			( (__m128i *)(imm) )
#define _addr_v2i32_2(imm)			( (__m128i *)(imm) )
#define _pv_v2i32(ptr)				( _addr_v2i32_1(ptr) )
/* expanders with pointers */
#define _e_p_v2i32_1(ptr)			_addr_v2i32_1(ptr)
#define _e_p_v2i32_2(ptr)			_addr_v2i32_2(ptr)
#define _e_pv_v2i32_1(ptr, a)		_addr_v2i32_1(ptr), (a).v1
#define _e_pv_v2i32_2(ptr, a)		_addr_v2i32_2(ptr), (a).v1

/* expand intrinsic name */
#define _i_v2i32(intrin) 			_mm_##intrin##_epi32
#define _i_v2i32e(intrin)			_mm_##intrin##_epi64
#define _i_v2i32x(intrin)			_mm_##intrin##_si128

/* apply */
#define _a_v2i32(intrin, expander, ...) ( \
	(v2i32_t) { \
		_i_v2i32(intrin)(expander##_v2i32_1(__VA_ARGS__)) \
	} \
)
#define _a_v2i32e(intrin, expander, ...) ( \
	(v2i32_t) { \
		_i_v2i32e(intrin)(expander##_v2i32_1(__VA_ARGS__)) \
	} \
)
#define _a_v2i32ev(intrin, expander, ...) { \
	_i_v2i32e(intrin)(expander##_v2i32_1(__VA_ARGS__)); \
}
#define _a_v2i32x(intrin, expander, ...) ( \
	(v2i32_t) { \
		_i_v2i32x(intrin)(expander##_v2i32_1(__VA_ARGS__)) \
	} \
)
#define _a_v2i32xv(intrin, expander, ...) { \
	_i_v2i32x(intrin)(expander##_v2i32_1(__VA_ARGS__)); \
}

/* load and store */
#define _load_v2i32(...)	_a_v2i32e(loadl, _e_p, __VA_ARGS__)
#define _loadu_v2i32(...)	_a_v2i32e(loadl, _e_p, __VA_ARGS__)
#define _store_v2i32(...)	_a_v2i32ev(storel, _e_pv, __VA_ARGS__)
#define _storeu_v2i32(...)	_a_v2i32ev(storel, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v2i32(...)		_a_v2i32(set1, _e_i, __VA_ARGS__)
#define _zero_v2i32()		_a_v2i32x(setzero, _e_x, _unused)
#define _seta_v2i32(x, y) ( \
	(v2i32_t) { \
		_mm_cvtsi64_si128((((uint64_t)(x))<<32) | ((uint32_t)(y))) \
	} \
)
#define _swap_v2i32(x) ( \
	(v2i32_t) { \
		_mm_shuffle_epi32((x).v1, 0xe1) \
	} \
)

/* logics */
#define _not_v2i32(...)		_a_v2i32x(not, _e_v, __VA_ARGS__)
#define _and_v2i32(...)		_a_v2i32x(and, _e_vv, __VA_ARGS__)
#define _or_v2i32(...)		_a_v2i32x(or, _e_vv, __VA_ARGS__)
#define _xor_v2i32(...)		_a_v2i32x(xor, _e_vv, __VA_ARGS__)
#define _andn_v2i32(...)	_a_v2i32x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v2i32(...)		_a_v2i32(add, _e_vv, __VA_ARGS__)
#define _sub_v2i32(...)		_a_v2i32(sub, _e_vv, __VA_ARGS__)
#define _mul_v2i32(...)		_a_v2i32(mul, _e_vv, __VA_ARGS__)
#define _max_v2i32(...)		_a_v2i32(max, _e_vv, __VA_ARGS__)
#define _min_v2i32(...)		_a_v2i32(min, _e_vv, __VA_ARGS__)

/* blend: (mask & b) | (~mask & a) */
#define _sel_v2i32(mask, a, b) ( \
	(v2i32_t) { \
		_mm_blendv_epi8((b).v1, (a).v1, (mask).v1) \
	} \
)

/* compare */
#define _eq_v2i32(...)		_a_v2i32(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v2i32(...)		_a_v2i32(cmpgt, _e_vv, __VA_ARGS__)

/* test: take mask and test if all zero */
#define _test_v2i32(x, y)	_mm_test_all_zeros((x).v1, (y).v1)

/* insert and extract */
#define _ins_v2i32(a, val, imm) { \
	(a).v1 = _i_v2i32((a).v1, (val), (imm)); \
}
#define _ext_v2i32(a, imm) ( \
	(int32_t)_i_v2i32(extract)((a).v1, (imm)) \
)

/* shift */
#define _sal_v2i32(a, imm) ( \
	(v2i32_t) {_i_v2i32(slli)((a).v1, (imm))} \
)
#define _sar_v2i32(a, imm) ( \
	(v2i32_t) {_i_v2i32(srai)((a).v1, (imm))} \
)

/* mask */
#define _mask_v2i32(a) ( \
	(uint32_t) (0xff & _mm_movemask_epi8((a).v1)) \
)
#define V2I32_MASK_00		( 0x00 )
#define V2I32_MASK_01		( 0x0f )
#define V2I32_MASK_10		( 0xf0 )
#define V2I32_MASK_11		( 0xff )

/* convert */
typedef uint64_t v2i8_t;
#define _load_v2i8(p) ({ \
	uint8_t const *_p = (uint8_t const *)(p); \
	*((uint16_t const *)(_p)); \
})
#define _store_v2i8(p, v) { \
	uint8_t *_p = (uint8_t *)(p); \
	*((uint16_t *)_p) = (v); \
}
#define _cvt_v2i8_v2i32(a) ( \
	(v2i32_t) { _mm_cvtepi8_epi32(_mm_cvtsi64_si128(a)) } \
)
#define _cvt_v2i32_v2i8(a) ( \
	(uint16_t)_mm_cvtsi128_si64(_mm_shuffle_epi8((a).v1, _mm_cvtsi64_si128(0x0400))) \
)
#define _cvt_u64_v2i32(a) ( \
	(v2i32_t){ _mm_cvtsi64_si128(a) } \
)
#define _cvt_v2i32_u64(a) ( \
	(uint64_t)_mm_cvtsi128_si64((a).v1) \
)
#define _cvt_v2i64_v2i32(a) ( \
	(v2i32_t) { \
		_mm_shuffle_epi32((a).v1, 0xd8) \
	} \
)
#define _cvth_v2i64_v2i32(a) ( \
	(v2i32_t) { \
		_mm_shuffle_epi32((a).v1, 0x8d) \
	} \
)

/* transpose */
#define _lo_v2i32(a, b) ( \
	(v2i32_t) { \
		_mm_unpacklo_epi32((a).v1, (b).v1) \
	} \
)
#define _hi_v2i32(a, b) ( \
	(v2i32_t) { \
		_mm_shuffle_epi32(_mm_unpacklo_epi32((a).v1, (b).v1), 0x0e) \
	} \
)

/* debug print */
#ifdef _LOG_H_INCLUDED
#define _print_v2i32(a) { \
	debug("(v2i32_t) %s(%d, %d)", #a, _ext_v2i32(a, 1), _ext_v2i32(a, 0)); \
}
#else
#define _print_v2i32(x)		;
#endif

#endif /* _V2I32_H_INCLUDED */
/**
 * end of v2i32.h
 */
