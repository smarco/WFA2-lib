
/**
 * @file arch_util.h
 *
 * @brief architecture-dependent utilities devided from util.h
 */
#ifndef _ARCH_UTIL_H_INCLUDED
#define _ARCH_UTIL_H_INCLUDED
#define MM_ARCH			"SSE4.1"

#include <x86intrin.h>
#include <stdint.h>

/**
 * misc bit operations (popcnt, tzcnt, and lzcnt)
 */

/**
 * @macro popcnt
 */
#ifdef __POPCNT__
#  define popcnt(x)		( (uint64_t)_mm_popcnt_u64(x) )
#else
	// #warning "popcnt instruction is not enabled."
	static inline
	uint64_t popcnt(uint64_t n)
	{
		uint64_t c = 0;
		c = (n & 0x5555555555555555) + ((n>>1) & 0x5555555555555555);
		c = (c & 0x3333333333333333) + ((c>>2) & 0x3333333333333333);
		c = (c & 0x0f0f0f0f0f0f0f0f) + ((c>>4) & 0x0f0f0f0f0f0f0f0f);
		c = (c & 0x00ff00ff00ff00ff) + ((c>>8) & 0x00ff00ff00ff00ff);
		c = (c & 0x0000ffff0000ffff) + ((c>>16) & 0x0000ffff0000ffff);
		c = (c & 0x00000000ffffffff) + ((c>>32) & 0x00000000ffffffff);
		return(c);
	}
#endif

/**
 * @macro ZCNT_RESULT
 * @brief workaround for a bug in gcc (<= 5), all the results of tzcnt / lzcnt macros must be modified by this label
 */
#ifndef ZCNT_RESULT
#  if defined(_ARCH_GCC_VERSION) && _ARCH_GCC_VERSION < 600
#    define ZCNT_RESULT		volatile
#  else
#    define ZCNT_RESULT
#  endif
#endif

/**
 * @macro tzcnt
 * @brief trailing zero count (count #continuous zeros from LSb)
 */
#ifdef __BMI__
/** immintrin.h is already included */
#  if defined(_ARCH_GCC_VERSION) && _ARCH_GCC_VERSION < 480
#    define tzcnt(x)		( (uint64_t)__tzcnt_u64(x) )
#  else
#    define tzcnt(x)		( (uint64_t)_tzcnt_u64(x) )
#  endif
#else
	// #warning "tzcnt instruction is not enabled."
	static inline
	uint64_t tzcnt(uint64_t n)
	{
		#ifdef __POPCNT__
			return(popcnt(~n & (n - 1)));
		#else
			if(n == 0) {
				return(64);
			} else {
				int64_t res;
				__asm__( "bsfq %1, %0" : "=r"(res) : "r"(n) );
				return(res);
			}
		#endif
	}
#endif

/**
 * @macro lzcnt
 * @brief leading zero count (count #continuous zeros from MSb)
 */
#ifdef __LZCNT__
/* __lzcnt_u64 in bmiintrin.h gcc-4.6, _lzcnt_u64 in lzcntintrin.h from gcc-4.7 */
#  if defined(_ARCH_GCC_VERSION) && _ARCH_GCC_VERSION < 470
#    define lzcnt(x)		( (uint64_t)__lzcnt_u64(x) )
#  else
#    define lzcnt(x)		( (uint64_t)_lzcnt_u64(x) )
#  endif
#else
	// #warning "lzcnt instruction is not enabled."
	static inline
	uint64_t lzcnt(uint64_t n)
	{
		if(n == 0) {
			return(64);
		} else {
			int64_t res;
			__asm__( "bsrq %1, %0" : "=r"(res) : "r"(n) );
			return(63 - res);
		}
	}
#endif

/**
 * @macro _swap_u64
 */
#ifdef __clang__
#  define _swap_u64(x)		({ uint64_t _x = (x); __asm__( "bswapq %0" : "+r"(_x) ); _x; })
#else
#  define _swap_u64(x)		( (uint64_t)_bswap64(x) )
#endif

/**
 * @macro _loadu_u64, _storeu_u64
 */
#define _loadu_u64(p)		({ uint8_t const *_p = (uint8_t const *)(p); *((uint64_t const *)_p); })
#define _storeu_u64(p, e)	{ uint8_t *_p = (uint8_t *)(p); *((uint64_t *)(_p)) = (e); }

/**
 * @macro _aligned_block_memcpy
 *
 * @brief copy size bytes from src to dst.
 *
 * @detail
 * src and dst must be aligned to 16-byte boundary.
 * copy must be multipe of 16.
 */
#define _xmm_rd_a(src, n) (xmm##n) = _mm_load_si128((__m128i *)(src) + (n))
#define _xmm_rd_u(src, n) (xmm##n) = _mm_loadu_si128((__m128i *)(src) + (n))
#define _xmm_wr_a(dst, n) _mm_store_si128((__m128i *)(dst) + (n), (xmm##n))
#define _xmm_wr_u(dst, n) _mm_storeu_si128((__m128i *)(dst) + (n), (xmm##n))
#define _memcpy_blk_intl(dst, src, size, _wr, _rd) { \
	/** duff's device */ \
	uint8_t *_src = (uint8_t *)(src), *_dst = (uint8_t *)(dst); \
	uint64_t const _nreg = 16;		/** #xmm registers == 16 */ \
	uint64_t const _tcnt = (size) / sizeof(__m128i); \
	uint64_t const _offset = ((_tcnt - 1) & (_nreg - 1)) - (_nreg - 1); \
	uint64_t _jmp = _tcnt & (_nreg - 1); \
	uint64_t _lcnt = (_tcnt + _nreg - 1) / _nreg; \
	register __m128i xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7; \
	register __m128i xmm8, xmm9, xmm10, xmm11, xmm12, xmm13, xmm14, xmm15; \
	_src += _offset * sizeof(__m128i); \
	_dst += _offset * sizeof(__m128i); \
	switch(_jmp) { \
		case 0: do { _rd(_src, 0); \
		case 15:     _rd(_src, 1); \
		case 14:     _rd(_src, 2); \
		case 13:     _rd(_src, 3); \
		case 12:     _rd(_src, 4); \
		case 11:     _rd(_src, 5); \
		case 10:     _rd(_src, 6); \
		case 9:      _rd(_src, 7); \
		case 8:      _rd(_src, 8); \
		case 7:      _rd(_src, 9); \
		case 6:      _rd(_src, 10); \
		case 5:      _rd(_src, 11); \
		case 4:      _rd(_src, 12); \
		case 3:      _rd(_src, 13); \
		case 2:      _rd(_src, 14); \
		case 1:      _rd(_src, 15); \
		switch(_jmp) { \
			case 0:  _wr(_dst, 0); \
			case 15: _wr(_dst, 1); \
			case 14: _wr(_dst, 2); \
			case 13: _wr(_dst, 3); \
			case 12: _wr(_dst, 4); \
			case 11: _wr(_dst, 5); \
			case 10: _wr(_dst, 6); \
			case 9:  _wr(_dst, 7); \
			case 8:  _wr(_dst, 8); \
			case 7:  _wr(_dst, 9); \
			case 6:  _wr(_dst, 10); \
			case 5:  _wr(_dst, 11); \
			case 4:  _wr(_dst, 12); \
			case 3:  _wr(_dst, 13); \
			case 2:  _wr(_dst, 14); \
			case 1:  _wr(_dst, 15); \
		} \
				     _src += _nreg * sizeof(__m128i); \
				     _dst += _nreg * sizeof(__m128i); \
				     _jmp = 0; \
			    } while(--_lcnt > 0); \
	} \
}
#define _memcpy_blk_aa(dst, src, len)		_memcpy_blk_intl(dst, src, len, _xmm_wr_a, _xmm_rd_a)
#define _memcpy_blk_au(dst, src, len)		_memcpy_blk_intl(dst, src, len, _xmm_wr_a, _xmm_rd_u)
#define _memcpy_blk_ua(dst, src, len)		_memcpy_blk_intl(dst, src, len, _xmm_wr_u, _xmm_rd_a)
#define _memcpy_blk_uu(dst, src, len)		_memcpy_blk_intl(dst, src, len, _xmm_wr_u, _xmm_rd_u)
#define _memset_blk_intl(dst, a, size, _wr) { \
	uint8_t *_dst = (uint8_t *)(dst); \
	__m128i const xmm0 = _mm_set1_epi8((int8_t)a); \
	uint64_t i; \
	for(i = 0; i < size / sizeof(__m128i); i++) { \
		_wr(_dst, 0); _dst += sizeof(__m128i); \
	} \
}
#define _memset_blk_a(dst, a, size)			_memset_blk_intl(dst, a, size, _xmm_wr_a)
#define _memset_blk_u(dst, a, size)			_memset_blk_intl(dst, a, size, _xmm_wr_u)

/**
 * substitution matrix abstraction
 */
/* store */
#define _store_sb(_scv, sv16)				{ _store_v16i8((_scv).v1, (sv16)); }

/* load */
#define _load_sb(scv)						( _from_v16i8_n(_load_v16i8((scv).v1)) )


/**
 * gap penalty vector abstraction macros
 */
/* store */
#define _make_gap(_e1, _e2, _e3, _e4) ( \
	(v16i8_t){ _mm_set_epi8( \
		(_e4), (_e4), (_e4), (_e4), \
		(_e3), (_e3), (_e3), (_e3), \
		(_e2), (_e2), (_e2), (_e2), \
		(_e1), (_e1), (_e1), (_e1)) \
	} \
)
#define _store_adjh(_scv, _adjh, _adjv, _ofsh, _ofsv) { \
	_store_v16i8((_scv).v2, _make_gap(_adjh, _adjv, _ofsh, _ofsv)) \
}
#define _store_adjv(_scv, _adjh, _adjv, _ofsh, _ofsv) { \
	/* nothing to do; use v2 */ \
	/* _store_v16i8((_scv).v3, _make_gap(_adjh, _adjv, _ofsh, _ofsv)) */ \
}
#define _store_ofsh(_scv, _adjh, _adjv, _ofsh, _ofsv) { \
	/* nothing to do; use v2 */ \
	/* _store_v16i8((_scv).v4, _make_gap(_adjh, _adjv, _ofsh, _ofsv)) */ \
}
#define _store_ofsv(_scv, _adjh, _adjv, _ofsh, _ofsv) { \
	/* nothing to do; use v2 */ \
	/* _store_v16i8((_scv).v5, _make_gap(_adjh, _adjv, _ofsh, _ofsv)) */ \
}

/* load */
#define _load_gap(_ptr, _idx) ( \
	(v16i8_t){ _mm_shuffle_epi32(_mm_load_si128((__m128i const *)(_ptr)), (_idx)) } \
)

#define _load_adjh(_scv)					( _from_v16i8_n(_load_gap((_scv).v2, 0x00)) )
#define _load_adjv(_scv)					( _from_v16i8_n(_load_gap((_scv).v2, 0x00)) )
#define _load_ofsh(_scv)					( _from_v16i8_n(_load_gap((_scv).v2, 0x55)) )
#define _load_ofsv(_scv)					( _from_v16i8_n(_load_gap((_scv).v2, 0x55)) )
#define _load_gfh(_scv)						( _from_v16i8_n(_load_gap((_scv).v2, 0xaa)) )
#define _load_gfv(_scv)						( _from_v16i8_n(_load_gap((_scv).v2, 0xff)) )
/*
#define _load_adjv(_scv)					( _from_v16i8_n(_load_gap((_scv).v2, 0x55)) )
#define _load_ofsv(_scv)					( _from_v16i8_n(_load_gap((_scv).v2, 0xff)) )
*/


/* cache line operation */
#define WCR_BUF_SIZE		( 128 )		/** two cache lines in x86_64 */
#define memcpy_buf(_dst, _src) { \
	register __m128i *_s = (__m128i *)(_src); \
	register __m128i *_d = (__m128i *)(_dst); \
	__m128i xmm0 = _mm_load_si128(_s); \
	__m128i xmm1 = _mm_load_si128(_s + 1); \
	__m128i xmm2 = _mm_load_si128(_s + 2); \
	__m128i xmm3 = _mm_load_si128(_s + 3); \
	__m128i xmm4 = _mm_load_si128(_s + 4); \
	__m128i xmm5 = _mm_load_si128(_s + 5); \
	__m128i xmm6 = _mm_load_si128(_s + 6); \
	__m128i xmm7 = _mm_load_si128(_s + 7); \
	_mm_stream_si128(_d, xmm0); \
	_mm_stream_si128(_d + 1, xmm1); \
	_mm_stream_si128(_d + 2, xmm2); \
	_mm_stream_si128(_d + 3, xmm3); \
	_mm_stream_si128(_d + 4, xmm4); \
	_mm_stream_si128(_d + 5, xmm5); \
	_mm_stream_si128(_d + 6, xmm6); \
	_mm_stream_si128(_d + 7, xmm7); \
}

/* 128bit register operation */
#define elem_128_t			__m128i
#define rd_128(_ptr)		( _mm_load_si128((__m128i *)(_ptr)) )
#define wr_128(_ptr, _e)	{ _mm_store_si128((__m128i *)(_ptr), (_e)); }
#define _ex_128(k, h)		( _mm_extract_epi64((elem_128_t)k, h) )
#define ex_128(k, p)		( ((((p)>>3) ? _ex_128(k, 1) : (_ex_128(k, 0))>>(((p) & 0x07)<<3)) & (WCR_OCC_SIZE-1)) )
#define p_128(v)			( _mm_cvtsi64_si128((uint64_t)(v)) )
#define e_128(v)			( (uint64_t)_mm_cvtsi128_si64((__m128i)(v)) )



/* compare and swap (cas) */
#if defined(__GNUC__)
#  if defined(_ARCH_GCC_VERSION) && _ARCH_GCC_VERSION < 470
#    define cas(ptr, cmp, val) ({ \
		uint8_t _res; \
		__asm__ volatile ("lock cmpxchg %[src], %[dst]\n\tsete %[res]" \
			: [dst]"+m"(*ptr), [res]"=a"(_res) \
			: [src]"r"(val), "a"(*cmp) \
			: "memory", "cc"); \
		_res; \
	})
#    define fence() ({ \
		__asm__ volatile ("mfence"); \
	})
#  else											/* > 4.7 */
#    define cas(ptr, cmp, val)	__atomic_compare_exchange_n(ptr, cmp, val, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)
#    define fence()				__sync_synchronize()
#  endif
#else
#  error "atomic compare-and-exchange is not supported in this version of compiler."
#endif

#endif /* #ifndef _ARCH_UTIL_H_INCLUDED */
/**
 * end of arch_util.h
 */
