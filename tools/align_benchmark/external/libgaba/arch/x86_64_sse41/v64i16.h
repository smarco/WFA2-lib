
/**
 * @file v64i16.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V64I16_H_INCLUDED
#define _V64I16_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 16bit 64cell */
typedef struct v64i16_s {
	__m128i v1;
	__m128i v2;
	__m128i v3;
	__m128i v4;
	__m128i v5;
	__m128i v6;
	__m128i v7;
	__m128i v8;
} v64i16_t;

/* expanders (without argument) */
#define _e_x_v64i16_1(u)
#define _e_x_v64i16_2(u)
#define _e_x_v64i16_3(u)
#define _e_x_v64i16_4(u)
#define _e_x_v64i16_5(u)
#define _e_x_v64i16_6(u)
#define _e_x_v64i16_7(u)
#define _e_x_v64i16_8(u)

/* expanders (without immediate) */
#define _e_v_v64i16_1(a)				(a).v1
#define _e_v_v64i16_2(a)				(a).v2
#define _e_v_v64i16_3(a)				(a).v3
#define _e_v_v64i16_4(a)				(a).v4
#define _e_v_v64i16_5(a)				(a).v5
#define _e_v_v64i16_6(a)				(a).v6
#define _e_v_v64i16_7(a)				(a).v7
#define _e_v_v64i16_8(a)				(a).v8
#define _e_vv_v64i16_1(a, b)			(a).v1, (b).v1
#define _e_vv_v64i16_2(a, b)			(a).v2, (b).v2
#define _e_vv_v64i16_3(a, b)			(a).v3, (b).v3
#define _e_vv_v64i16_4(a, b)			(a).v4, (b).v4
#define _e_vv_v64i16_5(a, b)			(a).v5, (b).v5
#define _e_vv_v64i16_6(a, b)			(a).v6, (b).v6
#define _e_vv_v64i16_7(a, b)			(a).v7, (b).v7
#define _e_vv_v64i16_8(a, b)			(a).v8, (b).v8
#define _e_vvv_v64i16_1(a, b, c)		(a).v1, (b).v1, (c).v1
#define _e_vvv_v64i16_2(a, b, c)		(a).v2, (b).v2, (c).v2
#define _e_vvv_v64i16_3(a, b, c)		(a).v3, (b).v3, (c).v3
#define _e_vvv_v64i16_4(a, b, c)		(a).v4, (b).v4, (c).v4
#define _e_vvv_v64i16_5(a, b, c)		(a).v5, (b).v5, (c).v5
#define _e_vvv_v64i16_6(a, b, c)		(a).v6, (b).v6, (c).v6
#define _e_vvv_v64i16_7(a, b, c)		(a).v7, (b).v7, (c).v7
#define _e_vvv_v64i16_8(a, b, c)		(a).v8, (b).v8, (c).v8

/* expanders with immediate */
#define _e_i_v64i16_1(imm)			(imm)
#define _e_i_v64i16_2(imm)			(imm)
#define _e_i_v64i16_3(imm)			(imm)
#define _e_i_v64i16_4(imm)			(imm)
#define _e_i_v64i16_5(imm)			(imm)
#define _e_i_v64i16_6(imm)			(imm)
#define _e_i_v64i16_7(imm)			(imm)
#define _e_i_v64i16_8(imm)			(imm)
#define _e_vi_v64i16_1(a, imm)		(a).v1, (imm)
#define _e_vi_v64i16_2(a, imm)		(a).v2, (imm)
#define _e_vi_v64i16_3(a, imm)		(a).v3, (imm)
#define _e_vi_v64i16_4(a, imm)		(a).v4, (imm)
#define _e_vi_v64i16_5(a, imm)		(a).v5, (imm)
#define _e_vi_v64i16_6(a, imm)		(a).v6, (imm)
#define _e_vi_v64i16_7(a, imm)		(a).v7, (imm)
#define _e_vi_v64i16_8(a, imm)		(a).v8, (imm)
#define _e_vvi_v64i16_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v64i16_2(a, b, imm)	(a).v2, (b).v2, (imm)
#define _e_vvi_v64i16_3(a, b, imm)	(a).v3, (b).v3, (imm)
#define _e_vvi_v64i16_4(a, b, imm)	(a).v4, (b).v4, (imm)
#define _e_vvi_v64i16_5(a, b, imm)	(a).v5, (b).v5, (imm)
#define _e_vvi_v64i16_6(a, b, imm)	(a).v6, (b).v6, (imm)
#define _e_vvi_v64i16_7(a, b, imm)	(a).v7, (b).v7, (imm)
#define _e_vvi_v64i16_8(a, b, imm)	(a).v8, (b).v8, (imm)

/* address calculation macros */
#define _addr_v64i16_1(imm)			( (__m128i *)(imm) )
#define _addr_v64i16_2(imm)			( (__m128i *)(imm) + 1 )
#define _addr_v64i16_3(imm)			( (__m128i *)(imm) + 2 )
#define _addr_v64i16_4(imm)			( (__m128i *)(imm) + 3 )
#define _addr_v64i16_5(imm)			( (__m128i *)(imm) + 4 )
#define _addr_v64i16_6(imm)			( (__m128i *)(imm) + 5 )
#define _addr_v64i16_7(imm)			( (__m128i *)(imm) + 6 )
#define _addr_v64i16_8(imm)			( (__m128i *)(imm) + 7 )
#define _pv_v64i16(ptr)				( _addr_v64i16_1(ptr) )

/* expanders with pointers */
#define _e_p_v64i16_1(ptr)			_addr_v64i16_1(ptr)
#define _e_p_v64i16_2(ptr)			_addr_v64i16_2(ptr)
#define _e_p_v64i16_3(ptr)			_addr_v64i16_3(ptr)
#define _e_p_v64i16_4(ptr)			_addr_v64i16_4(ptr)
#define _e_p_v64i16_5(ptr)			_addr_v64i16_5(ptr)
#define _e_p_v64i16_6(ptr)			_addr_v64i16_6(ptr)
#define _e_p_v64i16_7(ptr)			_addr_v64i16_7(ptr)
#define _e_p_v64i16_8(ptr)			_addr_v64i16_8(ptr)
#define _e_pv_v64i16_1(ptr, a)		_addr_v64i16_1(ptr), (a).v1
#define _e_pv_v64i16_2(ptr, a)		_addr_v64i16_2(ptr), (a).v2
#define _e_pv_v64i16_3(ptr, a)		_addr_v64i16_3(ptr), (a).v3
#define _e_pv_v64i16_4(ptr, a)		_addr_v64i16_4(ptr), (a).v4
#define _e_pv_v64i16_5(ptr, a)		_addr_v64i16_5(ptr), (a).v5
#define _e_pv_v64i16_6(ptr, a)		_addr_v64i16_6(ptr), (a).v6
#define _e_pv_v64i16_7(ptr, a)		_addr_v64i16_7(ptr), (a).v7
#define _e_pv_v64i16_8(ptr, a)		_addr_v64i16_8(ptr), (a).v8

/* expand intrinsic name */
#define _i_v64i16(intrin) 			_mm_##intrin##_epi16
#define _i_v64i16x(intrin)			_mm_##intrin##_si128

/* apply */
#define _a_v64i16(intrin, expander, ...) ( \
	(v64i16_t) { \
		_i_v64i16(intrin)(expander##_v64i16_1(__VA_ARGS__)), \
		_i_v64i16(intrin)(expander##_v64i16_2(__VA_ARGS__)), \
		_i_v64i16(intrin)(expander##_v64i16_3(__VA_ARGS__)), \
		_i_v64i16(intrin)(expander##_v64i16_4(__VA_ARGS__)), \
		_i_v64i16(intrin)(expander##_v64i16_5(__VA_ARGS__)), \
		_i_v64i16(intrin)(expander##_v64i16_6(__VA_ARGS__)), \
		_i_v64i16(intrin)(expander##_v64i16_7(__VA_ARGS__)), \
		_i_v64i16(intrin)(expander##_v64i16_8(__VA_ARGS__)) \
	} \
)
#define _a_v64i16x(intrin, expander, ...) ( \
	(v64i16_t) { \
		_i_v64i16x(intrin)(expander##_v64i16_1(__VA_ARGS__)), \
		_i_v64i16x(intrin)(expander##_v64i16_2(__VA_ARGS__)), \
		_i_v64i16x(intrin)(expander##_v64i16_3(__VA_ARGS__)), \
		_i_v64i16x(intrin)(expander##_v64i16_4(__VA_ARGS__)), \
		_i_v64i16x(intrin)(expander##_v64i16_5(__VA_ARGS__)), \
		_i_v64i16x(intrin)(expander##_v64i16_6(__VA_ARGS__)), \
		_i_v64i16x(intrin)(expander##_v64i16_7(__VA_ARGS__)), \
		_i_v64i16x(intrin)(expander##_v64i16_8(__VA_ARGS__)) \
	} \
)
#define _a_v64i16xv(intrin, expander, ...) { \
	_i_v64i16x(intrin)(expander##_v64i16_1(__VA_ARGS__)); \
	_i_v64i16x(intrin)(expander##_v64i16_2(__VA_ARGS__)); \
	_i_v64i16x(intrin)(expander##_v64i16_3(__VA_ARGS__)); \
	_i_v64i16x(intrin)(expander##_v64i16_4(__VA_ARGS__)); \
	_i_v64i16x(intrin)(expander##_v64i16_5(__VA_ARGS__)); \
	_i_v64i16x(intrin)(expander##_v64i16_6(__VA_ARGS__)); \
	_i_v64i16x(intrin)(expander##_v64i16_7(__VA_ARGS__)); \
	_i_v64i16x(intrin)(expander##_v64i16_8(__VA_ARGS__)); \
}

/* load and store */
#define _load_v64i16(...)	_a_v64i16x(load, _e_p, __VA_ARGS__)
#define _loadu_v64i16(...)	_a_v64i16x(loadu, _e_p, __VA_ARGS__)
#define _store_v64i16(...)	_a_v64i16xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v64i16(...)	_a_v64i16xv(storeu, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v64i16(...)	_a_v64i16(set1, _e_i, __VA_ARGS__)
#define _zero_v64i16()		_a_v64i16x(setzero, _e_x, _unused)

/* logics */
#define _not_v64i16(...)	_a_v64i16x(not, _e_v, __VA_ARGS__)
#define _and_v64i16(...)	_a_v64i16x(and, _e_vv, __VA_ARGS__)
#define _or_v64i16(...)		_a_v64i16x(or, _e_vv, __VA_ARGS__)
#define _xor_v64i16(...)	_a_v64i16x(xor, _e_vv, __VA_ARGS__)
#define _andn_v64i16(...)	_a_v64i16x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v64i16(...)	_a_v64i16(add, _e_vv, __VA_ARGS__)
#define _sub_v64i16(...)	_a_v64i16(sub, _e_vv, __VA_ARGS__)
#define _max_v64i16(...)	_a_v64i16(max, _e_vv, __VA_ARGS__)
#define _min_v64i16(...)	_a_v64i16(min, _e_vv, __VA_ARGS__)

/* compare */
#define _eq_v64i16(...)		_a_v64i16(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v64i16(...)		_a_v64i16(cmpgt, _e_vv, __VA_ARGS__)

/* insert and extract */
#define _V64I16_N			( sizeof(__m128i) / sizeof(int16_t) )
#define _ins_v64i16(a, val, imm) { \
	if((imm) < _V64I16_N) { \
		(a).v1 = _i_v64i8(insert)((a).v1, (val), (imm)); \
	} else if((imm) < 2 * _V64I16_N) { \
		(a).v2 = _i_v64i8(insert)((a).v2, (val), (imm) - _V64I16_N); \
	} else if((imm) < 3 * _V64I16_N) { \
		(a).v3 = _i_v64i8(insert)((a).v3, (val), (imm) - 2 * _V64I16_N); \
	} else if((imm) < 4 * _V64I16_N) { \
		(a).v4 = _i_v64i8(insert)((a).v4, (val), (imm) - 3 * _V64I16_N); \
	} else if((imm) < 5 * _V64I16_N) { \
		(a).v5 = _i_v64i8(insert)((a).v5, (val), (imm) - 4 * _V64I16_N); \
	} else if((imm) < 6 * _V64I16_N) { \
		(a).v6 = _i_v64i8(insert)((a).v6, (val), (imm) - 5 * _V64I16_N); \
	} else if((imm) < 7 * _V64I16_N) { \
		(a).v7 = _i_v64i8(insert)((a).v7, (val), (imm) - 6 * _V64I16_N); \
	} else { \
		(a).v8 = _i_v64i8(insert)((a).v8, (val), (imm) - 7 * _V64I16_N); \
	} \
}
#define _ext_v64i16(a, imm) ( \
	(int16_t)(((imm) < _V64I16_N) ? ( \
		_i_v64i16(extract)((a).v1, (imm)) \
	) : ((imm) < 2 * _V64I16_N) ? ( \
		_i_v64i16(extract)((a).v2, (imm) - _V64I16_N) \
	) : ((imm) < 3 * _V64I16_N) ? ( \
		_i_v64i16(extract)((a).v3, (imm) - 2 * _V64I16_N) \
	) : ((imm) < 4 * _V64I16_N) ? ( \
		_i_v64i16(extract)((a).v4, (imm) - 3 * _V64I16_N) \
	) : ((imm) < 5 * _V64I16_N) ? ( \
		_i_v64i16(extract)((a).v5, (imm) - 4 * _V64I16_N) \
	) : ((imm) < 6 * _V64I16_N) ? ( \
		_i_v64i16(extract)((a).v6, (imm) - 5 * _V64I16_N) \
	) : ((imm) < 7 * _V64I16_N) ? ( \
		_i_v64i16(extract)((a).v7, (imm) - 6 * _V64I16_N) \
	) : ( \
		_i_v64i16(extract)((a).v8, (imm) - 7 * _V64I16_N) \
	)) \
)

/* mask */
#define _mask_v64i16(a) ( \
	(v64_mask_t) { \
		.m1 = _mm_movemask_epi8(_mm_packs_epi16((a).v1, (a).v2)), \
		.m2 = _mm_movemask_epi8(_mm_packs_epi16((a).v3, (a).v4)), \
		.m3 = _mm_movemask_epi8(_mm_packs_epi16((a).v5, (a).v6)), \
		.m4 = _mm_movemask_epi8(_mm_packs_epi16((a).v7, (a).v8)) \
	} \
)

/* horizontal max (reduction max) */
#define _hmax_v64i16(a) ({ \
	__m128i _vmax = _mm_max_epi16( \
		_mm_max_epi16(_mm_max_epi16((a).v1, (a).v2), _mm_max_epi16((a).v3, (a).v4)), \
		_mm_max_epi16(_mm_max_epi16((a).v5, (a).v6), _mm_max_epi16((a).v7, (a).v8))); \
	_vmax = _mm_max_epi16(_vmax, \
		_mm_srli_si128(_vmax, 8)); \
	_vmax = _mm_max_epi16(_vmax, \
		_mm_srli_si128(_vmax, 4)); \
	_vmax = _mm_max_epi16(_vmax, \
		_mm_srli_si128(_vmax, 2)); \
	(int16_t)_mm_extract_epi16(_vmax, 0); \
})

#define _cvt_v64i8_v64i16(a) ( \
	(v64i16_t) { \
		_mm_cvtepi8_epi16((a).v1), \
		_mm_cvtepi8_epi16(_mm_srli_si128((a).v1, 8)), \
		_mm_cvtepi8_epi16((a).v2), \
		_mm_cvtepi8_epi16(_mm_srli_si128((a).v2, 8)), \
		_mm_cvtepi8_epi16((a).v3), \
		_mm_cvtepi8_epi16(_mm_srli_si128((a).v3, 8)), \
		_mm_cvtepi8_epi16((a).v4), \
		_mm_cvtepi8_epi16(_mm_srli_si128((a).v4, 8)) \
	} \
)

/* debug print */
#ifdef DEBUG
#define _print_v64i16(a) { \
	int16_t _tmp[64]; v64i16_t _tmpv = (a); _storeu_v64i16(_tmp, _tmpv); \
	debug("(v64i16_t) %s(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
			"%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, ", \
			"%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, ", \
			"%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", \
		#a, \
		_tmp[32 + 31], \
		_tmp[32 + 30], \
		_tmp[32 + 29], \
		_tmp[32 + 28], \
		_tmp[32 + 27], \
		_tmp[32 + 26], \
		_tmp[32 + 25], \
		_tmp[32 + 24], \
		_tmp[32 + 23], \
		_tmp[32 + 22], \
		_tmp[32 + 21], \
		_tmp[32 + 20], \
		_tmp[32 + 19], \
		_tmp[32 + 18], \
		_tmp[32 + 17], \
		_tmp[32 + 16], \
		_tmp[32 + 15], \
		_tmp[32 + 14], \
		_tmp[32 + 13], \
		_tmp[32 + 12], \
		_tmp[32 + 11], \
		_tmp[32 + 10], \
		_tmp[32 + 9], \
		_tmp[32 + 8], \
		_tmp[32 + 7], \
		_tmp[32 + 6], \
		_tmp[32 + 5], \
		_tmp[32 + 4], \
		_tmp[32 + 3], \
		_tmp[32 + 2], \
		_tmp[32 + 1], \
		_tmp[32 + 0], \
		_tmp[31], \
		_tmp[30], \
		_tmp[29], \
		_tmp[28], \
		_tmp[27], \
		_tmp[26], \
		_tmp[25], \
		_tmp[24], \
		_tmp[23], \
		_tmp[22], \
		_tmp[21], \
		_tmp[20], \
		_tmp[19], \
		_tmp[18], \
		_tmp[17], \
		_tmp[16], \
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
#define _print_v64i16(x)	;
#endif

#endif /* _V64I16_H_INCLUDED */
/**
 * end of v64i16.h
 */
