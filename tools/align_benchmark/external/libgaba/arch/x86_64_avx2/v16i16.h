
/**
 * @file v16i16.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V16I16_H_INCLUDED
#define _V16I16_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 8bit 32cell */
typedef struct v16i16_s {
	__m256i v1;
} v16i16_t;

/* expanders (without argument) */
#define _e_x_v16i16_1(u)
#define _e_x_v16i16_2(u)

/* expanders (without immediate) */
#define _e_v_v16i16_1(a)			(a).v1
#define _e_v_v16i16_2(a)			(a).v1
#define _e_vv_v16i16_1(a, b)		(a).v1, (b).v1
#define _e_vv_v16i16_2(a, b)		(a).v1, (b).v1
#define _e_vvv_v16i16_1(a, b, c)	(a).v1, (b).v1, (c).v1
#define _e_vvv_v16i16_2(a, b, c)	(a).v1, (b).v1, (c).v1

/* expanders with immediate */
#define _e_i_v16i16_1(imm)			(imm)
#define _e_i_v16i16_2(imm)			(imm)
#define _e_vi_v16i16_1(a, imm)		(a).v1, (imm)
#define _e_vi_v16i16_2(a, imm)		(a).v1, (imm)
#define _e_vvi_v16i16_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v16i16_2(a, b, imm)	(a).v1, (b).v1, (imm)

/* address calculation macros */
#define _addr_v16i16_1(imm)			( (__m256i *)(imm) )
#define _addr_v16i16_2(imm)			( (__m256i *)(imm) )
#define _pv_v16i16(ptr)				( _addr_v16i16_1(ptr) )
/* expanders with pointers */
#define _e_p_v16i16_1(ptr)			_addr_v16i16_1(ptr)
#define _e_p_v16i16_2(ptr)			_addr_v16i16_2(ptr)
#define _e_pv_v16i16_1(ptr, a)		_addr_v16i16_1(ptr), (a).v1
#define _e_pv_v16i16_2(ptr, a)		_addr_v16i16_2(ptr), (a).v1

/* expand intrinsic name */
#define _i_v16i16(intrin) 			_mm256_##intrin##_epi16
#define _i_v16i16x(intrin)			_mm256_##intrin##_si256

/* apply */
#define _a_v16i16(intrin, expander, ...) ( \
	(v16i16_t) { \
		_i_v16i16(intrin)(expander##_v16i16_1(__VA_ARGS__)) \
	} \
)
#define _a_v16i16x(intrin, expander, ...) ( \
	(v16i16_t) { \
		_i_v16i16x(intrin)(expander##_v16i16_1(__VA_ARGS__)) \
	} \
)
#define _a_v16i16xv(intrin, expander, ...) { \
	_i_v16i16x(intrin)(expander##_v16i16_1(__VA_ARGS__)); \
}

/* load and store */
#define _load_v16i16(...)	_a_v16i16x(load, _e_p, __VA_ARGS__)
#define _loadu_v16i16(...)	_a_v16i16x(loadu, _e_p, __VA_ARGS__)
#define _store_v16i16(...)	_a_v16i16xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v16i16(...)	_a_v16i16xv(storeu, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v16i16(...)	_a_v16i16(set1, _e_i, __VA_ARGS__)
#define _zero_v16i16()		_a_v16i16x(setzero, _e_x, _unused)
#define _seta_v16i16(...) ( \
	(v16i16_t) { \
		_mm256_set_epi16(__VA_ARGS__) \
	} \
)

/* swap (reverse) */
#define _swap_idx_v16i16() ( \
	_mm256_broadcastsi128_si256(_mm_set_epi8( \
		1, 0, 3, 2, 5, 4, 7, 6, \
		9, 8, 11, 10, 13, 12, 15, 14)) \
)
#define _swap_v16i16(a) ( \
	(v16i16_t) { \
		_mm256_permute2x128_si256( \
			_mm256_shuffle_epi8((a).v1, _swap_idx_v16i16()), \
			_mm256_shuffle_epi8((a).v1, _swap_idx_v16i16()), \
			0x01) \
	} \
)

/* logics */
#define _not_v16i16(...)	_a_v16i16x(not, _e_v, __VA_ARGS__)
#define _and_v16i16(...)	_a_v16i16x(and, _e_vv, __VA_ARGS__)
#define _or_v16i16(...)		_a_v16i16x(or, _e_vv, __VA_ARGS__)
#define _xor_v16i16(...)	_a_v16i16x(xor, _e_vv, __VA_ARGS__)
#define _andn_v16i16(...)	_a_v16i16x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v16i16(...)	_a_v16i16(add, _e_vv, __VA_ARGS__)
#define _sub_v16i16(...)	_a_v16i16(sub, _e_vv, __VA_ARGS__)
#define _max_v16i16(...)	_a_v16i16(max, _e_vv, __VA_ARGS__)
#define _min_v16i16(...)	_a_v16i16(min, _e_vv, __VA_ARGS__)

/* shuffle */
#define _shuf_v16i16(...)	_a_v16i16(shuffle, _e_vv, __VA_ARGS__)

/* blend */
// #define _sel_v16i16(...)	_a_v16i16(blendv, _e_vvv, __VA_ARGS__)

/* compare */
#define _eq_v16i16(...)		_a_v16i16(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v16i16(...)		_a_v16i16(cmpgt, _e_vv, __VA_ARGS__)

/* insert and extract */
#define _ins_v16i16(a, val, imm) { \
	(a).v1 = _i_v16i16(insert)((a).v1, (val), (imm)); \
}
#define _ext_v16i16(a, imm) ( \
	(int16_t)_i_v16i16(extract)((a).v1, (imm)) \
)

/* mask */
#define _mask_v16i16(a) ( \
	(v16_mask_t) { \
		.m1 = _mm256_movemask_epi8( \
			_mm256_packs_epi16((a).v1, \
			_mm256_castsi128_si256(_mm256_extracti128_si256((a).v1, 1)))) \
	} \
)

/* horizontal max (reduction max) */
#define _hmax_v16i16(a) ({ \
	__m128i _t = _mm_max_epi16( \
		_mm256_castsi256_si128((a).v1), \
		_mm256_extracti128_si256((a).v1, 1) \
	); \
	_t = _mm_max_epi16(_t, _mm_srli_si128(_t, 8)); \
	_t = _mm_max_epi16(_t, _mm_srli_si128(_t, 4)); \
	_t = _mm_max_epi16(_t, _mm_srli_si128(_t, 2)); \
	(int16_t)_mm_extract_epi16(_t, 0); \
})

/* convert */
#define _cvt_v16i8_v16i16(a) ( \
	(v16i16_t) { \
		_mm256_cvtepi8_epi16((a).v1) \
	} \
)

/* debug print */
#ifdef DEBUG
#define _print_v16i16(a) { \
	int16_t _tmp[16]; v16i16_t _tmpv = (a); _storeu_v16i16(_tmp, _tmpv); \
	debug("(v16i16_t) %s(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", \
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
#define _print_v16i16(x)		;
#endif

#endif /* _V16I16_H_INCLUDED */
/**
 * end of v16i16.h
 */
