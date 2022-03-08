
/**
 * @file v16i8.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V16I8_H_INCLUDED
#define _V16I8_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 8bit 32cell */
typedef struct v16i8_s {
	__m128i v1;
} v16i8_t;

/* expanders (without argument) */
#define _e_x_v16i8_1(u)
#define _e_x_v16i8_2(u)

/* expanders (without immediate) */
#define _e_v_v16i8_1(a)				(a).v1
#define _e_v_v16i8_2(a)				(a).v1
#define _e_vv_v16i8_1(a, b)			(a).v1, (b).v1
#define _e_vv_v16i8_2(a, b)			(a).v1, (b).v1
#define _e_vvv_v16i8_1(a, b, c)		(a).v1, (b).v1, (c).v1
#define _e_vvv_v16i8_2(a, b, c)		(a).v1, (b).v1, (c).v1

/* expanders with immediate */
#define _e_i_v16i8_1(imm)			(imm)
#define _e_i_v16i8_2(imm)			(imm)
#define _e_vi_v16i8_1(a, imm)		(a).v1, (imm)
#define _e_vi_v16i8_2(a, imm)		(a).v1, (imm)
#define _e_vvi_v16i8_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v16i8_2(a, b, imm)	(a).v1, (b).v1, (imm)

/* address calculation macros */
#define _addr_v16i8_1(imm)			( (__m128i *)(imm) )
#define _addr_v16i8_2(imm)			( (__m128i *)(imm) )
#define _pv_v16i8(ptr)				( _addr_v16i8_1(ptr) )
/* expanders with pointers */
#define _e_p_v16i8_1(ptr)			_addr_v16i8_1(ptr)
#define _e_p_v16i8_2(ptr)			_addr_v16i8_2(ptr)
#define _e_pv_v16i8_1(ptr, a)		_addr_v16i8_1(ptr), (a).v1
#define _e_pv_v16i8_2(ptr, a)		_addr_v16i8_2(ptr), (a).v1

/* expand intrinsic name */
#define _i_v16i8(intrin) 			_mm_##intrin##_epi8
#define _i_v16i8x(intrin)			_mm_##intrin##_si128

/* apply */
#define _a_v16i8(intrin, expander, ...) ( \
	(v16i8_t) { \
		_i_v16i8(intrin)(expander##_v16i8_1(__VA_ARGS__)) \
	} \
)
#define _a_v16i8x(intrin, expander, ...) ( \
	(v16i8_t) { \
		_i_v16i8x(intrin)(expander##_v16i8_1(__VA_ARGS__)) \
	} \
)
#define _a_v16i8xv(intrin, expander, ...) { \
	_i_v16i8x(intrin)(expander##_v16i8_1(__VA_ARGS__)); \
}

/* load and store */
#define _load_v16i8(...)	_a_v16i8x(load, _e_p, __VA_ARGS__)
#define _loadu_v16i8(...)	_a_v16i8x(loadu, _e_p, __VA_ARGS__)
#define _store_v16i8(...)	_a_v16i8xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v16i8(...)	_a_v16i8xv(storeu, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v16i8(...)		_a_v16i8(set1, _e_i, __VA_ARGS__)
#define _zero_v16i8()		_a_v16i8x(setzero, _e_x, _unused)
#define _seta_v16i8(...) ( \
	(v16i8_t) { \
		_mm_set_epi8(__VA_ARGS__) \
	} \
)

/* swap (reverse) */
#define _swap_idx_v16i8() ( \
	_mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15) \
)
#define _swap_v16i8(a) ( \
	(v16i8_t) { \
		_mm_shuffle_epi8((a).v1, _swap_idx_v16i8()) \
	} \
)
#define _swapn_idx_v16i8() ( \
	_mm_set_epi8(-16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1) \
)
#define _swapn_v16i8(a, l) ( \
	(v16i8_t) { \
		_mm_shuffle_epi8((a).v1, _mm_add_epi8(_swapn_idx_v16i8(), _mm_set1_epi8(l))) \
	} \
)

/* logics */
#define _not_v16i8(...)		_a_v16i8x(not, _e_v, __VA_ARGS__)
#define _and_v16i8(...)		_a_v16i8x(and, _e_vv, __VA_ARGS__)
#define _or_v16i8(...)		_a_v16i8x(or, _e_vv, __VA_ARGS__)
#define _xor_v16i8(...)		_a_v16i8x(xor, _e_vv, __VA_ARGS__)
#define _andn_v16i8(...)	_a_v16i8x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v16i8(...)		_a_v16i8(add, _e_vv, __VA_ARGS__)
#define _sub_v16i8(...)		_a_v16i8(sub, _e_vv, __VA_ARGS__)
#define _adds_v16i8(...)	_a_v16i8(adds, _e_vv, __VA_ARGS__)
#define _subs_v16i8(...)	_a_v16i8(subs, _e_vv, __VA_ARGS__)
#define _max_v16i8(...)		_a_v16i8(max, _e_vv, __VA_ARGS__)
#define _min_v16i8(...)		_a_v16i8(min, _e_vv, __VA_ARGS__)

/* shuffle */
#define _shuf_v16i8(...)	_a_v16i8(shuffle, _e_vv, __VA_ARGS__)

/* blend */
#define _sel_v16i8(...)		_a_v16i8(blendv, _e_vvv, __VA_ARGS__)

/* compare */
#define _eq_v16i8(...)		_a_v16i8(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v16i8(...)		_a_v16i8(cmpgt, _e_vv, __VA_ARGS__)

/* insert and extract */
#define _ins_v16i8(a, val, imm) { \
	(a).v1 = _i_v16i8(insert)((a).v1, (val), (imm)); \
}
#define _ext_v16i8(a, imm) ( \
	(int8_t)_i_v16i8(extract)((a).v1, (imm)) \
)

/* shift */
#define _bsl_v16i8(a, imm) ( \
	(v16i8_t) { \
		_i_v16i8x(slli)((a).v1, (imm)) \
	} \
)
#define _bsr_v16i8(a, imm) ( \
	(v16i8_t) { \
		_i_v16i8x(srli)((a).v1, (imm)) \
	} \
)
#define _shl_v16i8(a, imm) ( \
	(v16i8_t) { \
		_mm_slli_epi32((a).v1, (imm)) \
	} \
)
#define _shr_v16i8(a, imm) ( \
	(v16i8_t) { \
		_mm_srli_epi32((a).v1, (imm)) \
	} \
)
#define _sal_v16i8(a, imm) ( \
	(v16i8_t) { \
		_mm_slai_epi32((a).v1, (imm)) \
	} \
)
#define _sar_v16i8(a, imm) ( \
	(v16i8_t) { \
		_mm_srai_epi32((a).v1, (imm)) \
	} \
)

/* mask */
#define _mask_v16i8(a) ( \
	(v16_mask_t) { \
		.m1 = _i_v16i8(movemask)((a).v1) \
	} \
)

/* horizontal max */
#define _hmax_v16i8(a) ({ \
	__m128i _vmax = _mm_max_epi8((a).v1, \
		_mm_srli_si128((a).v1, 8)); \
	_vmax = _mm_max_epi8(_vmax, \
		_mm_srli_si128(_vmax, 4)); \
	_vmax = _mm_max_epi8(_vmax, \
		_mm_srli_si128(_vmax, 2)); \
	_vmax = _mm_max_epi8(_vmax, \
		_mm_srli_si128(_vmax, 1)); \
	(int8_t)_mm_extract_epi8(_vmax, 0); \
})

/* convert */
#define _cvt_v16i16_v16i8(a) ( \
	(v16i8_t) { \
		_mm_packs_epi16(_mm256_castsi256_si128((a).v1), _mm256_extracti128_si256((a).v1, 1)) \
	} \
)

/* debug print */
#ifdef DEBUG
#define _print_v16i8(a) { \
	int8_t _tmp[16]; v16i8_t _tmpv = (a); _storeu_v16i8(_tmp, _tmpv); \
	debug("(v16i8_t) %s(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", \
		#a, \
		_tmp[15], \
		_tmp[14], \
		_tmp[13], \
		_tmp[12], \
		_tmp[11], \
		_tmp[10], \
		_tmp[9], \
		_tmp[8], \
		_tmp[7], \
		_tmp[6], \
		_tmp[5], \
		_tmp[4], \
		_tmp[3], \
		_tmp[2], \
		_tmp[1], \
		_tmp[0]); \
}
#else
#define _print_v16i8(x)		;
#endif

#endif /* _V16I8_H_INCLUDED */
/**
 * end of v16i8.h
 */
