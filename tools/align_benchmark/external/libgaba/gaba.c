
/**
 * @file gaba.c
 *
 * @brief libgaba (libsea3) DP routine implementation
 *
 * @author Hajime Suzuki
 * @date 2016/1/11
 * @license Apache v2
 */
// #define DEBUG
// #define DEBUG_MEM
// #define DEBUG_OVERFLOW
// #define DEBUG_ALL
/*
 * debug print configuration: -DDEBUG to enable debug print, -DDEBUG_ALL to print all the vectors, arrays, and bitmasks
 * NOTE: dumping all the vectors sometimes raises SEGV due to a stack shortage. use `ulimit -s 65536' to avoid it.
 */
#if defined(DEBUG_ALL) && !defined(DEBUG)
#  define DEBUG
#endif
#ifdef DEBUG
#  define REDEFINE_DEBUG
#endif

/* make sure POSIX APIs are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif

#if defined(__darwin__) && !defined(_BSD_SOURCE)
#  define _BSD_SOURCE
#endif

/* import general headers */
#include <stdio.h>				/* sprintf in dump_path */
#include <stdint.h>				/* uint32_t, uint64_t, ... */
#include <stddef.h>				/* offsetof */
#include <string.h>				/* memset, memcpy */
#include <inttypes.h>

#define _GABA_PARSE_EXPORT_LEVEL	static inline
#include "gaba.h"
#include "gaba_parse.h"
#include "log.h"
#include "arch/arch.h"			/* architecture dependents */


/**
 * gap penalty model configuration: choose one of the following three by a define macro
 *   linear:   g(k) =          ge * k          where              ge > 0
 *   affine:   g(k) =     gi + ge * k          where gi > 0,      ge > 0
 *   combined: g(k) = min(gi + ge * k, gf * k) where gi > 0, gf > ge > 0
 */
#define LINEAR						0
#define AFFINE						1
#define COMBINED					2

#ifdef MODEL
#  if !(MODEL == LINEAR || MODEL == AFFINE || MODEL == COMBINED)
#    error "MODEL must be LINEAR (1), AFFINE (2), or COMBINED (3)."
#  endif
#else
#  define MODEL					AFFINE
#endif

#if MODEL == LINEAR
#  define MODEL_LABEL				linear
#  define MODEL_STR					"linear"
#elif MODEL == AFFINE
#  define MODEL_LABEL				affine
#  define MODEL_STR					"affine"
#else
#  define MODEL_LABEL				combined
#  define MODEL_STR					"combined"
#endif

#ifdef BIT
#  if !(BIT == 2 || BIT == 4)
#    error "BIT must be 2 or 4."
#  endif
#else
#  define BIT						4
#endif

/* define ENABLE_FILTER to enable gapless alignment filter */
// #define ENABLE_FILTER


/* bandwidth-specific configurations aliasing vector macros */
#define BW_MAX						64
#ifndef BW
#  define BW						64
#endif

#define _W							BW
#if _W == 16
#  define _NVEC_ALIAS_PREFIX		v16i8
#  define _WVEC_ALIAS_PREFIX		v16i16
#  define DP_CTX_INDEX				2
#elif _W == 32
#  define _NVEC_ALIAS_PREFIX		v32i8
#  define _WVEC_ALIAS_PREFIX		v32i16
#  define DP_CTX_INDEX				1
#elif _W == 64
#  define _NVEC_ALIAS_PREFIX		v64i8
#  define _WVEC_ALIAS_PREFIX		v64i16
#  define DP_CTX_INDEX				0
#else
#  error "BW must be one of 16, 32, or 64."
#endif
#include "arch/vector_alias.h"

#define DP_CTX_MAX					( 3 )
#define _dp_ctx_index(_bw)			( ((_bw) == 64) ? 0 : (((_bw) == 32) ? 1 : 2) )
// _static_assert(_dp_ctx_index(BW) == DP_CTX_INDEX);


/* add suffix for gap-model- and bandwidth-wrapper (see gaba_wrap.h) */
#ifdef SUFFIX
#  define _suffix_cat3_2(a, b, c)	a##_##b##_##c
#  define _suffix_cat3(a, b, c)		_suffix_cat3_2(a, b, c)
#  define _suffix(_base)			_suffix_cat3(_base, MODEL_LABEL, BW)
#else
#  define _suffix(_base)			_base
#endif


/* add namespace for arch wrapper (see main.c) */
#ifdef NAMESPACE
#  define _export_cat(x, y)			x##_##y
#  define _export_cat2(x, y)		_export_cat(x, y)
#  define _export(_base)			_export_cat2(_suffix(_base), NAMESPACE)
#else
#  define _export(_base)			_suffix(_base)
#endif


/* import unittest */
#ifndef UNITTEST_UNIQUE_ID
#  if MODEL == LINEAR
#    if BW == 16
#      define UNITTEST_UNIQUE_ID	31
#    elif BW == 32
#      define UNITTEST_UNIQUE_ID	32
#    else
#      define UNITTEST_UNIQUE_ID	33
#    endif
#  elif MODEL == AFFINE
#    if BW == 16
#      define UNITTEST_UNIQUE_ID	34
#    elif BW == 32
#      define UNITTEST_UNIQUE_ID	35
#    else
#      define UNITTEST_UNIQUE_ID	36
#    endif
#  else
#    if BW == 16
#      define UNITTEST_UNIQUE_ID	37
#    elif BW == 32
#      define UNITTEST_UNIQUE_ID	38
#    else
#      define UNITTEST_UNIQUE_ID	39
#    endif
#  endif
#endif
#include  "unittest.h"


/* static assertion macros */
#define _sa_cat_intl(x, y)		x##y
#define _sa_cat(x, y)			_sa_cat_intl(x, y)
#define _static_assert(expr)	typedef char _sa_cat(_st, __LINE__)[(expr) ? 1 : -1]


/* internal constants */
#define BLK_BASE					( 5 )
#define BLK						( 0x01<<BLK_BASE )

#ifdef DEBUG_MEM
#  define MIN_BULK_BLOCKS			( 0 )
#  define MEM_ALIGN_SIZE			( 32 )		/* 32byte aligned for AVX2 environments */
#  define MEM_INIT_SIZE				( (uint64_t)4 * 1024 )
#  define MEM_MARGIN_SIZE			( 4096 )	/* tail margin of internal memory blocks */
#else
#  define MIN_BULK_BLOCKS			( 32 )
#  define MEM_ALIGN_SIZE			( 32 )		/* 32byte aligned for AVX2 environments */
#  define MEM_INIT_SIZE				( (uint64_t)256 * 1024 * 1024  )
#  define MEM_MARGIN_SIZE			( 4096 )	/* tail margin of internal memory blocks */
#endif

#define INIT_FETCH_APOS				( -1 )
#define INIT_FETCH_BPOS				( -1 )

/* test consistency of exported macros */
_static_assert(V2I32_MASK_01 == GABA_UPDATE_A);
_static_assert(V2I32_MASK_10 == GABA_UPDATE_B);


/**
 * @macro _likely, _unlikely
 * @brief branch prediction hint for gcc-compatible compilers
 */
#define _likely(x)					__builtin_expect(!!(x), 1)
#define _unlikely(x)				__builtin_expect(!!(x), 0)

/**
 * @macro _force_inline
 * @brief inline directive for gcc-compatible compilers
 */
#define _force_inline				inline
// #define _force_inline				/* */

/** assume 64bit little-endian system */
_static_assert(sizeof(void *) == 8);

/** check size of structs declared in gaba.h */
_static_assert(sizeof(struct gaba_params_s) == 40);
_static_assert(sizeof(struct gaba_section_s) == 16);
_static_assert(sizeof(struct gaba_fill_s) == 64);
_static_assert(sizeof(struct gaba_segment_s) == 32);
_static_assert(sizeof(struct gaba_alignment_s) == 64);
_static_assert(sizeof(nvec_masku_t) == _W / 8);

/**
 * @macro _plen
 * @brief extract plen from path_section_s
 */
#define _plen(seg)					( (seg)->alen + (seg)->blen )

/* forward declarations */
static int64_t gaba_dp_add_stack(struct gaba_dp_context_s *self, uint64_t size);
static void *gaba_dp_malloc(struct gaba_dp_context_s *self, uint64_t size);
static void gaba_dp_free(struct gaba_dp_context_s *self, void *ptr);				/* do nothing */
struct gaba_dp_context_s;


/**
 * @struct gaba_drop_s
 */
struct gaba_drop_s {
	int8_t drop[_W];					/** (32) max */
};
_static_assert(sizeof(struct gaba_drop_s) == _W);

/**
 * @struct gaba_middle_delta_s
 */
struct gaba_middle_delta_s {
	int16_t delta[_W];					/** (64) middle delta */
};
_static_assert(sizeof(struct gaba_middle_delta_s) == sizeof(int16_t) * _W);

/**
 * @struct gaba_mask_pair_s
 */
#if MODEL == LINEAR
struct gaba_mask_pair_s {
	nvec_masku_t h;						/** (4) horizontal mask vector */
	nvec_masku_t v;						/** (4) vertical mask vector */
};
_static_assert(sizeof(struct gaba_mask_pair_s) == _W / 4);
#else	/* affine and combined */
struct gaba_mask_pair_s {
	nvec_masku_t h;						/** (4) horizontal mask vector */
	nvec_masku_t v;						/** (4) vertical mask vector */
	nvec_masku_t e;						/** (4) e mask vector */
	nvec_masku_t f;						/** (4) f mask vector */
};
_static_assert(sizeof(struct gaba_mask_pair_s) == _W / 2);
#endif
#define _mask_u64(_m)				( ((nvec_masku_t){ .mask = (_m) }).all )

/**
 * @struct gaba_diff_vec_s
 */
#if MODEL == LINEAR
struct gaba_diff_vec_s {
	uint8_t dh[_W];						/** (32) dh */
	uint8_t dv[_W];						/** (32) dv */
};
_static_assert(sizeof(struct gaba_diff_vec_s) == 2 * _W);
#else	/* affine and combined gap penalty */
struct gaba_diff_vec_s {
	uint8_t dh[_W];						/** (32) dh */
	uint8_t dv[_W];						/** (32) dv */
	uint8_t de[_W];						/** (32) de */
	uint8_t df[_W];						/** (32) df */
};
_static_assert(sizeof(struct gaba_diff_vec_s) == 4 * _W);
#endif

/**
 * @struct gaba_char_vec_s
 */
struct gaba_char_vec_s {
	uint8_t w[_W];						/** (32) a in the lower 4bit, b in the higher 4bit */
};
_static_assert(sizeof(struct gaba_char_vec_s) == _W);

/**
 * @struct gaba_block_s
 * @brief a unit of banded matrix, 32 vector updates will be recorded in a single block object.
 * phantom is an alias of the block struct as a head cap of contiguous blocks.
 */
struct gaba_block_s {
	struct gaba_mask_pair_s mask[BLK];	/** (256 / 512) traceback capability flag vectors (set if transition to the ajdacent cell is possible) */
	struct gaba_diff_vec_s diff;		/** (64, 128, 256) diff variables of the last vector */
	int8_t acc, xstat;					/** (2) accumulator, and xdrop status (term detected when xstat < 0) */
	int8_t acnt, bcnt;					/** (2) forwarded lengths */
	uint32_t dir_mask;					/** (4) extension direction bit array */
	uint64_t max_mask;					/** (8) lanewise update mask (set if the lane contains the current max) */
};
struct gaba_phantom_s {
	struct gaba_diff_vec_s diff;		/** (64, 128, 256) diff variables of the last (just before the head) vector */
	int8_t acc, xstat;					/** (2) accumulator, and xdrop status (term detected when xstat < 0) */
	int8_t acnt, bcnt;					/** (4) prefetched sequence lengths (only effective at the root, otherwise zero) */
	uint32_t reserved;					/** (4) overlaps with dir_mask */
	struct gaba_block_s const *blk;		/** (8) link to the previous block (overlaps with max_mask) */
};
_static_assert(sizeof(struct gaba_block_s) % 16 == 0);
_static_assert(sizeof(struct gaba_phantom_s) % 16 == 0);
#define _last_block(x)				( (struct gaba_block_s *)(x) - 1 )
#define _last_phantom(x)			( (struct gaba_phantom_s *)(x) - 1 )
#define _phantom(x)					( _last_phantom((struct gaba_block_s *)(x) + 1) )

/**
 * @struct gaba_merge_s
 */
struct gaba_merge_s {
	uint64_t reserved1;					/** (8) keep aligned to 16byte boundary */
	struct gaba_block_s const *blk[1];	/** (8) addressed by [tail_idx[vec][q]] */
	int8_t tidx[2][_W];					/** (32, 64, 128) lanewise index array */
	struct gaba_diff_vec_s diff;		/** (64, 128, 256) diff variables of the last (just before the head) vector */
	int8_t acc, xstat;					/** (2) acc and xstat are reserved for block_s */
	uint8_t reserved2[13], qofs[1];		/** (14) displacement of vectors in the q-direction */
};
#define MERGE_TAIL_OFFSET			( BLK * sizeof(struct gaba_mask_pair_s) - 2 * _W - 2 * sizeof(void *) )
#define _merge(x)					( (struct gaba_merge_s *)((uint8_t *)(x) + MERGE_TAIL_OFFSET) )
_static_assert(sizeof(struct gaba_merge_s) % 16 == 0);
_static_assert(sizeof(struct gaba_merge_s) + MERGE_TAIL_OFFSET == sizeof(struct gaba_block_s));
_static_assert(MAX_MERGE_COUNT <= 14);

/**
 * @struct gaba_joint_tail_s
 * @brief (internal) tail cap of a contiguous matrix blocks, contains a context of the blocks
 * (band) and can be connected to the next blocks.
 */
struct gaba_joint_tail_s {
	/* char vector and delta vectors */
	struct gaba_char_vec_s ch;			/** (16, 32, 64) char vector */
	struct gaba_drop_s xd;				/** (16, 32, 64) */
	struct gaba_middle_delta_s md;		/** (32, 64, 128) */

	int16_t mdrop;						/** (2) drop from m.max (offset) */
	uint16_t istat;						/** (2) 1 if bridge */
	uint32_t pridx;						/** (4) remaining p-length */
	uint32_t aridx, bridx;				/** (8) reverse indices for the tails */
	uint32_t aadv, badv;				/** (8) advanced lengths */
	struct gaba_joint_tail_s const *tail;/** (8) the previous tail */
	uint64_t abrk, bbrk;				/** (16) breakpoint masks */
	uint8_t const *atptr, *btptr;		/** (16) tail of the current section */
	struct gaba_fill_s f;				/** (24) */
};
_static_assert((sizeof(struct gaba_joint_tail_s) % 32) == 0);
#define TAIL_BASE				( offsetof(struct gaba_joint_tail_s, f) )
#define BRIDGE_TAIL_OFFSET		( offsetof(struct gaba_joint_tail_s, mdrop) )
#define BRIDGE_TAIL_SIZE		( sizeof(struct gaba_joint_tail_s) - BRIDGE_TAIL_OFFSET )
#define _bridge(x)				( (struct gaba_joint_tail_s *)((uint8_t *)(x) - BRIDGE_TAIL_OFFSET) )
#define _tail(x)				( (struct gaba_joint_tail_s *)((uint8_t *)(x) - TAIL_BASE) )
#define _fill(x)				( (struct gaba_fill_s *)((uint8_t *)(x) + TAIL_BASE) )
#define _offset(x)				( (x)->f.max - (x)->mdrop )

#define _mem_blocks(n)				( sizeof(struct gaba_phantom_s) + (n + 1) * sizeof(struct gaba_block_s) + sizeof(struct gaba_joint_tail_s) )
#define MEM_INIT_VACANCY			( _mem_blocks(MIN_BULK_BLOCKS) )
_static_assert(2 * sizeof(struct gaba_block_s) < MEM_MARGIN_SIZE);
_static_assert(MEM_INIT_VACANCY < MEM_INIT_SIZE);

/**
 * @struct gaba_root_block_s
 */
struct gaba_root_block_s {
	uint8_t _pad1[288 - sizeof(struct gaba_phantom_s)];
	struct gaba_phantom_s blk;
	struct gaba_joint_tail_s tail;
#if _W != 64
	uint8_t _pad2[352 + 32 - sizeof(struct gaba_joint_tail_s)];
#endif
};
_static_assert(sizeof(struct gaba_root_block_s) == 640 + 32);
_static_assert(sizeof(struct gaba_root_block_s) >= sizeof(struct gaba_phantom_s) + sizeof(struct gaba_joint_tail_s));

/**
 * @struct gaba_reader_work_s
 * @brief (internal) working buffer for fill-in functions, contains sequence prefetch buffer
 * and middle/max-small deltas.
 */
struct gaba_reader_work_s {
	/** 64byte aligned */
	uint8_t bufa[BW_MAX + BLK];			/** (128) */
	uint8_t bufb[BW_MAX + BLK];			/** (128) */
	/** 256 */

	/** 64byte alidned */
	uint32_t arlim, brlim;				/** (8) asridx - aadv = aridx + arlim */
	uint32_t aid, bid;					/** (8) ids */
	uint8_t const *atptr, *btptr;		/** (16) tail of the current section */
	uint32_t pridx;						/** (4) remaining p-length (unsigned!) */
	int32_t ofsd;						/** (4) delta of large offset */
	uint32_t arem, brem;				/** (8) current ridx (redefined as rem, ridx == rem + rlim) */
	uint32_t asridx, bsridx;			/** (8) start ridx (converted to (badv, aadv) in fill_create_tail) */
	struct gaba_joint_tail_s const *tail;	/** (8) previous tail */
	/** 64 */

	/** 64byte aligned */
	struct gaba_drop_s xd;				/** (16, 32, 64) current drop from max */
#if _W != 64
	uint8_t _pad[_W == 16 ? 16 : 32];	/** padding to align to 64-byte boundary */
#endif
	struct gaba_middle_delta_s md;		/** (32, 64, 128) */
};
_static_assert((sizeof(struct gaba_reader_work_s) % 64) == 0);

/**
 * @struct gaba_merge_work_s
 */
#define MERGE_BUFFER_LENGTH		( 2 * _W )
struct gaba_merge_work_s {
	uint8_t qofs[16];					/** (16) q-offset array */
	uint32_t qw, _pad1;					/** (8) */
	uint32_t lidx, uidx;				/** (8) */
#if _W != 16
	uint8_t _pad2[32];					/** padding to align to 64-byte boundary */
#endif
	uint64_t abrk[4];					/** (32) */
	uint64_t bbrk[4];					/** (32) */
	uint8_t buf[MERGE_BUFFER_LENGTH + 7 * sizeof(int16_t) * MERGE_BUFFER_LENGTH];	/** (32, 64, 128) + (320, 640, 1280) */
};
_static_assert((sizeof(struct gaba_merge_work_s) % 64) == 0);

/**
 * @struct gaba_aln_intl_s
 * @brief internal alias of gaba_alignment_t, allocated in the local context and used as a working buffer.
 */
struct gaba_aln_intl_s {
	/* memory management */
	void *opaque;						/** (8) opaque (context pointer) */
	gaba_lfree_t lfree;					/** (8) local free */
	int64_t score;						/** (8) score (save) */
	uint32_t aicnt, bicnt;				/** (8) gap region counters */
	uint32_t aecnt, becnt;				/** (8) gap base counters */
	uint32_t dcnt;						/** (4) unused in the loop */
	uint32_t slen;						/** (4) section length (counter) */
	struct gaba_segment_s *seg;			/** (8) section ptr */
	uint32_t plen, padding;				/** (8) path length (psum; save) */
};

_static_assert(sizeof(struct gaba_alignment_s) == sizeof(struct gaba_aln_intl_s));
_static_assert(offsetof(struct gaba_alignment_s, score) == offsetof(struct gaba_aln_intl_s, score));
_static_assert(offsetof(struct gaba_alignment_s, dcnt) == offsetof(struct gaba_aln_intl_s, dcnt));
_static_assert(offsetof(struct gaba_alignment_s, slen) == offsetof(struct gaba_aln_intl_s, slen));
_static_assert(offsetof(struct gaba_alignment_s, seg) == offsetof(struct gaba_aln_intl_s, seg));
_static_assert(offsetof(struct gaba_alignment_s, plen) == offsetof(struct gaba_aln_intl_s, plen));
_static_assert(offsetof(struct gaba_alignment_s, agcnt) == offsetof(struct gaba_aln_intl_s, aecnt));
_static_assert(offsetof(struct gaba_alignment_s, bgcnt) == offsetof(struct gaba_aln_intl_s, becnt));

/**
 * @struct gaba_leaf_s
 * @brief working buffer for max score search
 */
struct gaba_leaf_s {
	struct gaba_joint_tail_s const *tail;
	struct gaba_block_s const *blk;
	uint32_t p, q;						/** (8) local p (to restore mask pointer), local q */
	uint64_t ppos;						/** (8) global p (== resulting path length) */
	uint32_t aridx, bridx;
};

/**
 * @struct gaba_writer_work_s
 * @brief working buffer for traceback (allocated in the thread-local context)
 */
struct gaba_writer_work_s {
	/** local context */
	struct gaba_aln_intl_s a;			/** (64) working buffer, copied to the result object */

	/** work */
	uint32_t afcnt, bfcnt;				/** (8) */
	uint32_t *path;						/** (8) path array pointer */
	struct gaba_block_s const *blk;		/** (8) current block */
	uint8_t p, q, ofs, state, _pad[4];	/** (8) local p, q-coordinate, [0, BW), path offset, state */

	/** save */
	uint32_t aofs, bofs;				/** (8) ofs for bridge */
	uint32_t agidx, bgidx;				/** (8) grid indices of the current trace */
	uint32_t asgidx, bsgidx;			/** (8) base indices of the current trace */
	uint32_t aid, bid;					/** (8) section ids */

	/** section info */
	struct gaba_joint_tail_s const *atail;/** (8) */
	struct gaba_joint_tail_s const *btail;/** (8) */
	struct gaba_alignment_s *aln;		/** (8) */

	struct gaba_leaf_s leaf;			/** (40) working buffer for max pos search */
	/** 64, 192 */
};
_static_assert((sizeof(struct gaba_writer_work_s) % 64) == 0);

/**
 * @struct gaba_score_vec_s
 */
struct gaba_score_vec_s {
	int8_t v1[16];
	int8_t v2[16];
	int8_t v3[16];
	int8_t v4[16];
	int8_t v5[16];
};
_static_assert(sizeof(struct gaba_score_vec_s) == 80);

/**
 * @struct gaba_mem_block_s
 */
struct gaba_mem_block_s {
	struct gaba_mem_block_s *next;
	// struct gaba_mem_block_s *prev;
	uint64_t size;
};
_static_assert(sizeof(struct gaba_mem_block_s) == 16);

/**
 * @struct gaba_stack_s
 * @brief save stack pointer
 */
struct gaba_stack_s {
	struct gaba_mem_block_s *mem;
	uint8_t *top;
	uint8_t *end;
};
_static_assert(sizeof(struct gaba_stack_s) == 24);
#define _stack_size(_s)					( (uint64_t)((_s)->end - (_s)->top) )

/**
 * @macro _init_bar, _test_bar
 * @brief memory barrier for debugging
 */
#ifdef DEBUG_MEM
#define _init_bar(_name) { \
	memset(self->_barrier_##_name, 0xa5, 256); \
}
#define _test_bar(_name) { \
	for(uint64_t __i = 0; __i < 256; __i++) { \
		if(self->_barrier_##_name[__i] != 0xa5) { \
			fprintf(stderr, "barrier broken, i(%lu), m(%u)\n", __i, self->_barrier_##_name[__i]); *((volatile uint8_t *)NULL); \
		} \
	} \
}
#define _barrier(_name)					uint8_t _barrier_##_name[256];
#else
#define _init_bar(_name)				;
#define _test_bar(_name)				;
#define _barrier(_name)
#endif

/**
 * @struct gaba_dp_context_s
 *
 * @brief (internal) container for dp implementations
 */
struct gaba_dp_context_s {
	_barrier(head);

	/** loaded on init */
	struct gaba_joint_tail_s const *root[4];	/** (32) root tail (phantom vectors) */

	/* memory management */
	struct gaba_mem_block_s mem;		/** (16) root memory block */
	struct gaba_stack_s stack;			/** (24) current stack */
	uint64_t _pad2;

	/* score constants */
	double imx, xmx;					/** (16) 1 / (M - X), X / (M - X) (precalculated constants) */
	struct gaba_score_vec_s scv;		/** (80) substitution matrix and gaps */

	/* scores */
	int8_t tx;							/** (1) xdrop threshold */
	int8_t tf;							/** (1) filter threshold */
	int8_t gi, ge, gfa, gfb;			/** (4) negative integers */
	uint8_t aflen, bflen;				/** (2) short-gap length thresholds */
	uint8_t ofs,  _pad1[7];				/** (16) */
	/** 192; 64byte aligned */

	_barrier(mid);

	/* working buffers */
	union gaba_work_s {
		struct gaba_reader_work_s r;	/** (192) */
		struct gaba_merge_work_s m;		/** (2048?) */
		struct gaba_writer_work_s l;	/** (192) */
	} w;
	/** 64byte aligned */

	_barrier(tail);
};
_static_assert((sizeof(struct gaba_dp_context_s) % 64) == 0);
#define GABA_DP_CONTEXT_LOAD_OFFSET	( 0 )
#define GABA_DP_CONTEXT_LOAD_SIZE	( offsetof(struct gaba_dp_context_s, w) )
_static_assert((GABA_DP_CONTEXT_LOAD_OFFSET % 64) == 0);
_static_assert((GABA_DP_CONTEXT_LOAD_SIZE % 64) == 0);
#define _root(_t)					( (_t)->root[_dp_ctx_index(BW)] )

/**
 * @struct gaba_opaque_s
 */
struct gaba_opaque_s {
	void *api[8];
};
#define _export_dp_context(_t) ( \
	(struct gaba_dp_context_s *)(((struct gaba_opaque_s *)(_t)) - DP_CTX_MAX + _dp_ctx_index(BW)) \
)
#define _restore_dp_context(_t) ( \
	(struct gaba_dp_context_s *)(((struct gaba_opaque_s *)(_t)) - _dp_ctx_index(BW) + DP_CTX_MAX) \
)
#define _export_dp_context_global(_t) ( \
	(struct gaba_dp_context_s *)(((struct gaba_opaque_s *)(_t)) - DP_CTX_MAX + _dp_ctx_index(BW)) \
)
#define _restore_dp_context_global(_t) ( \
	(struct gaba_dp_context_s *)(((struct gaba_opaque_s *)(_t)) - _dp_ctx_index(BW) + DP_CTX_MAX) \
)

/**
 * @struct gaba_context_s
 *
 * @brief (API) an algorithmic context.
 *
 * @sa gaba_init, gaba_close
 */
#ifdef DEBUG_MEM
#  define ROOT_BLOCK_OFFSET				( 4096 + 2048 )
#else
#  define ROOT_BLOCK_OFFSET				( 4096 )
#endif
struct gaba_context_s {
	/** opaque pointers for function dispatch */
	struct gaba_opaque_s api[4];		/** function dispatcher, used in gaba_wrap.h */
	/** 64byte aligned */

	/** templates */
	struct gaba_dp_context_s dp;		/** template of thread-local context */
	uint8_t _pad[ROOT_BLOCK_OFFSET - sizeof(struct gaba_dp_context_s) - 4 * sizeof(struct gaba_opaque_s)];
	/** 64byte aligned */

	/** phantom vectors */
	struct gaba_root_block_s ph[3];		/** (768) template of root vectors, [0] for 16-cell, [1] for 32-cell, [2] for 64-cell */
	/** 64byte aligned */
};
_static_assert(offsetof(struct gaba_context_s, ph) == ROOT_BLOCK_OFFSET);
#define _proot(_c, _bw)				( &(_c)->ph[_dp_ctx_index(_bw)] )

/**
 * @enum GABA_BLK_STATUS
 * head states and intermediate states are exclusive
 */
enum GABA_BLK_STATUS {
	/* intermediate states */
	CONT			= 0,				/* continue */
	ZERO			= 0x01,				/* internal use */
	TERM			= 0x80,				/* sign bit */
	STAT_MASK		= ZERO | TERM | CONT,
	/* head states */
	HEAD			= 0x20,
	MERGE			= 0x40,				/* merged head and the corresponding block contains no actual vector (DP cell) */
	ROOT			= HEAD | MERGE
};
_static_assert((int8_t)TERM < 0);		/* make sure TERM is recognezed as a negative value */
_static_assert((int32_t)CONT<<8 == (int32_t)GABA_CONT);
_static_assert((int32_t)TERM<<8 == (int32_t)GABA_TERM);


/**
 * coordinate conversion macros
 */
// #define _rev(pos, len)				( (len) + (uint64_t)(len) - (uint64_t)(pos) - 1 )
#define _rev(pos)					( GABA_EOU + (uint64_t)GABA_EOU - (uint64_t)(pos) - 1 )
#define _roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )

/**
 * max and min
 */
#define MAX2(x,y)		( (x) > (y) ? (x) : (y) )
#define MAX3(x,y,z)		( MAX2(x, MAX2(y, z)) )
#define MAX4(w,x,y,z)	( MAX2(MAX2(w, x), MAX2(y, z)) )

#define MIN2(x,y)		( (x) < (y) ? (x) : (y) )
#define MIN3(x,y,z)		( MIN2(x, MIN2(y, z)) )
#define MIN4(w,x,y,z)	( MIN2(MIN2(w, x), MIN2(y, z)) )


/**
 * @fn gaba_malloc, gaba_free
 * @brief a pair of malloc and free, aligned and margined.
 * any pointer created by gaba_malloc MUST be freed by gaba_free.
 */
static _force_inline
void *gaba_malloc(
	size_t size)
{
	void *ptr = NULL;

	/* roundup to align boundary, add margin at the head and tail */
	size = _roundup(size, MEM_ALIGN_SIZE);
	if(posix_memalign(&ptr, MEM_ALIGN_SIZE, size + 2 * MEM_MARGIN_SIZE) != 0) {
		debug("posix_memalign failed");
		trap(); return(NULL);
	}
	debug("posix_memalign(%p), size(%lu)", ptr, size);
	return(ptr + MEM_MARGIN_SIZE);
}
static _force_inline
void gaba_free(
	void *ptr)
{
	debug("free(%p)", ptr - MEM_MARGIN_SIZE);
	free(ptr - MEM_MARGIN_SIZE);
	return;
}


/* matrix fill functions */

/* direction macros */
/**
 * @struct gaba_dir_s
 */
struct gaba_dir_s {
	uint32_t mask;
	int32_t acc;					/* use 32bit int to avoid (sometimes inefficient) 8bit and 16bit operations on x86_64 GP registers */
};

/**
 * @macro _dir_init
 */
#define _dir_init(_blk) ((struct gaba_dir_s){ .mask = 0, .acc = (_blk)->acc })
/**
 * @macro _dir_fetch
 */
#define _dir_fetch(_d) { \
	(_d).mask <<= 1; (_d).mask |= (uint32_t)((_d).acc < 0);  \
	debug("fetched dir(%x), %s", (_d).mask, _dir_is_down(_d) ? "go down" : "go right"); \
}
/**
 * @macro _dir_update
 * @brief update direction determiner for the next band
 */
#define _dir_update(_d, _vector) { \
	(_d).acc += _ext_n(_vector, 0) - _ext_n(_vector, _W-1); \
	/*debug("acc(%d), (%d, %d)", _dir_acc, _ext_n(_vector, 0), _ext_n(_vector, _W-1));*/ \
}
/**
 * @macro _dir_adjust_remainder
 * @brief adjust direction array when termination is detected in the middle of the block
 */
#define _dir_adjust_remainder(_d, _filled_count) { \
	debug("adjust remainder, array(%x), shifted array(%x)", (_d).mask, (_d).mask<<(BLK - (_filled_count))); \
	(_d).mask <<= (BLK - (_filled_count)); \
}
/**
 * @macro _dir_is_down, _dir_is_right
 * @brief direction indicator (_dir_is_down returns true if dir == down)
 */
#define _dir_mask_is_down(_mask)	( (_mask) & 0x01 )
#define _dir_mask_is_right(_mask)	( ~(_mask) & 0x01 )
#define _dir_is_down(_d)			( _dir_mask_is_down((_d).mask) )
#define _dir_is_right(_d)			( _dir_mask_is_right((_d).mask) )
#define _dir_bcnt(_d)				( popcnt((_d).mask) )				/* count vertical transitions */
#define _dir_mask_windback(_mask)	{ (_mask) >>= 1; }
#define _dir_windback(_d)			{ (_d).mask >>= 1; }
/**
 * @macro _dir_save
 */
#define _dir_save(_blk, _d) { \
	(_blk)->dir_mask = (_d).mask;	/* store mask */ \
	(_blk)->acc = (_d).acc;			/* store accumulator */ \
}
/**
 * @macro _dir_load
 */
#define _dir_mask_load(_blk, _cnt)	( (_blk)->dir_mask>>(BLK - (_cnt)) )
#define _dir_load(_blk, _cnt) ( \
	(struct gaba_dir_s){ \
		.mask = _dir_mask_load(_blk, _cnt), \
		.acc = (_blk)->acc \
	} \
	/*debug("load dir cnt(%d), mask(%x), shifted mask(%x)", (int32_t)_filled_count, _d.mask, _d.mask>>(BLK - (_filled_count)));*/ \
)

/**
 * seqreader macros
 */
#define _rd_bufa_base(k)		( (k)->w.r.bufa + BLK + _W )
#define _rd_bufb_base(k)		( (k)->w.r.bufb )
#define _rd_bufa(k, pos, len)	( _rd_bufa_base(k) - (pos) - (len) )
#define _rd_bufb(k, pos, len)	( _rd_bufb_base(k) + (pos) )
#define _lo64(v)				_ext_v2i64(v, 0)
#define _hi64(v)				_ext_v2i64(v, 1)
#define _lo32(v)				_ext_v2i32(v, 0)
#define _hi32(v)				_ext_v2i32(v, 1)

/* scoring function and sequence encoding */

/**
 * @macro _max_match, _gap_h, _gap_v
 * @brief calculate scores (_gap: 0 for horizontal gap, 1 for vertical gap)
 */
#define _max_match(_p)				( _hmax_v16i8(_loadu_v16i8((_p)->score_matrix)) )
#define _min_match(_p)				( -_hmax_v16i8(_sub_v16i8(_zero_v16i8(), _loadu_v16i8((_p)->score_matrix))) )
#if MODEL == LINEAR
#define _gap_h(_p, _l)				( -1 * ((_p)->gi + (_p)->ge) * (_l) )
#define _gap_v(_p, _l)				( -1 * ((_p)->gi + (_p)->ge) * (_l) )
#define _gap(_p, _d, _l)			_gap_h(_p, _l)
#elif MODEL == AFFINE
#define _gap_h(_p, _l)				( -1 * ((_l) > 0) * (_p)->gi - (_p)->ge * (_l) )
#define _gap_v(_p, _l)				( -1 * ((_l) > 0) * (_p)->gi - (_p)->ge * (_l) )
#define _gap(_p, _d, _l)			_gap_h(_p, _l)
#define _gap_e(_p, _l)				( -1 * ((_l) > 0) * (_p)->gi - (_p)->ge * (_l) )
#define _gap_f(_p, _l)				( -1 * ((_l) > 0) * (_p)->gi - (_p)->ge * (_l) )
#else /* MODEL == COMBINED */
#define _gap_h(_p, _l)				( MAX2(-1 * ((_l) > 0) * (_p)->gi - (_p)->ge * (_l), -1 * (_p)->gfb * (_l)) )
#define _gap_v(_p, _l)				( MAX2(-1 * ((_l) > 0) * (_p)->gi - (_p)->ge * (_l), -1 * (_p)->gfa * (_l)) )
#define _gap(_p, _d, _l)			( MAX2(-1 * ((_l) > 0) * (_p)->gi - (_p)->ge * (_l), -1 * (&(_p)->gfb)[-(_d)] * (_l)) )
#define _gap_e(_p, _l)				( -1 * ((_l) > 0) * (_p)->gi - (_p)->ge * (_l) )
#define _gap_f(_p, _l)				( -1 * ((_l) > 0) * (_p)->gi - (_p)->ge * (_l) )
#endif
#define _ofs_h(_p)					( (_p)->gi + (_p)->ge )
#define _ofs_v(_p)					( (_p)->gi + (_p)->ge )
#define _ofs_e(_p)					( (_p)->gi )
#define _ofs_f(_p)					( (_p)->gi )

/**
 * @enum BASES
 */
#if BIT == 2
enum BASES { A = 0x00, C = 0x01, G = 0x02, T = 0x03, N = 0x04 };
#  define _max_match_base_a(_p)		( 0x0c )
#  define _max_match_base_b(_p)		( 0x03 )
#else
enum BASES { A = 0x01, C = 0x02, G = 0x04, T = 0x08, N = 0x00 };
#  define _max_match_base_a(_p)		( 0x01 )
#  define _max_match_base_b(_p)		( 0x01 )
#endif

/**
 * @macro _adjust_BLK, _comp_BLK, _match_BW
 * NOTE: _load_v16i8 and _load_v32i8 will be moved out of the loop when loading static uint8_t const[16].
 */
#if BIT == 2
/* 2bit encoding */
static uint8_t const decode_table[16] __attribute__(( aligned(16) )) = {
	'A', 'C', 'G', 'T', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};
static uint8_t const comp_mask_a[16] __attribute__(( aligned(16) )) = {
	0x03, 0x02, 0x01, 0x00, 0x04, 0x04, 0x04, 0x04,
	0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04
};
static uint8_t const shift_mask_b[16] __attribute__(( aligned(16) )) = {
	0x00, 0x04, 0x08, 0x0c, 0x02, 0x02, 0x02, 0x02,
	0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02
};
static uint8_t const compshift_mask_b[16] __attribute__(( aligned(16) )) = {
	0x0c, 0x08, 0x04, 0x00, 0x02, 0x02, 0x02, 0x02,
	0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02
};
/* q-fetch (anti-diagonal matching for vector calculation) */
#ifdef UNSAFE_FETCH
#  define _fwaq_v16i8(_v)		( _shuf_v16i8((_load_v16i8(comp_mask_a)), (_v)) )
#  define _fwaq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(comp_mask_a))), (_v)) )
#  define _rvaq_v16i8(_v)		( _swap_v16i8((_v)) )
#  define _rvaq_v32i8(_v)		( _swap_v32i8((_v)) )
#  define _fwbq_v16i8(_v)		( _shuf_v16i8((_load_v16i8(shift_mask_b)), (_v)) )
#  define _fwbq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(shift_mask_b))), (_v)) )
#  define _rvbq_v16i8(_v)		( _shuf_v16i8((_load_v16i8(compshift_mask_b)), _swap_v16i8(_v)) )
#  define _rvbq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(compshift_mask_b))), _swap_v32i8(_v)) )
#else
#  define _fwaq_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(comp_mask_a)), (_v)) )	/* _l is ignored */
#  define _fwaq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(comp_mask_a))), (_v)) )
#  define _rvaq_v16i8(_v, _l)	( _swapn_v16i8((_v), (_l)) )
#  define _rvaq_v32i8(_v)		( _swap_v32i8((_v)) )
#  define _fwbq_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(shift_mask_b)), (_v)) )	/* _l is ignored */
#  define _fwbq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(shift_mask_b))), (_v)) )
#  define _rvbq_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(compshift_mask_b)), _swapn_v16i8((_v), (_l))) )
#  define _rvbq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(compshift_mask_b))), _swap_v32i8(_v)) )
#endif

/* p-fetch (diagonal matching for alignment refinement) */
#  define _fwap_v16i8(_v, _l)	( (_v) )
#  define _rvap_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(comp_mask_a)), _swapn_v16i8((_v), (_l))) )
#  define _fwbp_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(shift_mask_b)), (_v)) )
#  define _rvbp_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(compshift_mask_b)), _swapn_v16i8((_v), (_l))) )

#  define _match_n(_a, _b)		_or_n(_a, _b)
#  define _match_v16i8(_a, _b)	_or_v16i8(_a, _b)

#else
/* 4bit encoding */
static uint8_t const decode_table[16] __attribute__(( aligned(16) )) = {
	'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
	'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
};
static uint8_t const comp_mask[16] __attribute__(( aligned(16) )) = {
	0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
	0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
};
#define comp_mask_a				comp_mask
#define compshift_mask_b		comp_mask

/* q-fetch (anti-diagonal matching for vector calculation) */
#ifdef UNSAFE_FETCH
#  define _fwaq_v16i8(_v)		( _shuf_v16i8((_load_v16i8(comp_mask)), (_v)) )
#  define _fwaq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(comp_mask))), (_v)) )
#  define _rvaq_v16i8(_v)		( _swap_v16i8((_v)) )
#  define _rvaq_v32i8(_v)		( _swap_v32i8((_v)) )
#  define _fwbq_v16i8(_v)		( (_v) )
#  define _fwbq_v32i8(_v)		( (_v) )
#  define _rvbq_v16i8(_v)		( _shuf_v16i8((_load_v16i8(comp_mask)), _swap_v16i8(_v)) )
#  define _rvbq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(comp_mask))), _swap_v32i8(_v)) )
#else
#  define _fwaq_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(comp_mask)), (_v)) )		/* _l is ignored */
#  define _fwaq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(comp_mask))), (_v)) )
#  define _rvaq_v16i8(_v, _l)	( _swapn_v16i8((_v), (_l)) )
#  define _rvaq_v32i8(_v)		( _swap_v32i8((_v)) )
#  define _fwbq_v16i8(_v, _l)	( (_v) )											/* id(x); _l is ignored */
#  define _fwbq_v32i8(_v)		( (_v) )
#  define _rvbq_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(comp_mask)), _swapn_v16i8((_v), (_l))) )
#  define _rvbq_v32i8(_v)		( _shuf_v32i8((_from_v16i8_v32i8(_load_v16i8(comp_mask))), _swap_v32i8(_v)) )
#endif

/* p-fetch (diagonal matching for alignment refinement) */
#  define _fwap_v16i8(_v, _l)	( (_v) )
#  define _rvap_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(comp_mask)), _swapn_v16i8((_v), (_l))) )
#  define _fwbp_v16i8(_v, _l)	( (_v) )
#  define _rvbp_v16i8(_v, _l)	( _shuf_v16i8((_load_v16i8(comp_mask)), _swapn_v16i8((_v), (_l))) )

#  define _match_n(_a, _b)		_and_n(_a, _b)
#  define _match_v16i8(_a, _b)	_and_v16i8(_a, _b)
#endif

/**
 * @fn fill_fetch_seq_a
 */
static _force_inline
void fill_fetch_seq_a(
	struct gaba_dp_context_s *self,
	uint8_t const *pos,
	uint64_t len)
{
	if(pos < GABA_EOU) {
		debug("reverse fetch a: pos(%p), len(%lu)", pos, len);
		/* reverse fetch: 2 * alen - (2 * alen - pos) + (len - 32) */
		#ifdef UNSAFE_FETCH
		v32i8_t ach = _loadu_v32i8(pos + (len - BLK));					/* this may touch the space before the array */
		_storeu_v32i8(_rd_bufa(self, _W, len), _rvaq_v32i8(ach));		/* reverse */
		#else
		v32i8_t ach = _loadu_v32i8(pos);								/* this will not touch the space before the array, but will touch at most 31bytes after the array */
		_storeu_v32i8(_rd_bufa(self, _W, BLK), _rvaq_v32i8(ach));	/* reverse; will not invade any buf */
		#endif
	} else {
		debug("forward fetch a: pos(%p), len(%lu), p(%p)", pos, len, _rev(pos + (len - 1)));
		/* forward fetch: 2 * alen - pos */
		#ifdef UNSAFE_FETCH
		v32i8_t ach = _loadu_v32i8(_rev(pos + (len - 1)));
		_storeu_v32i8(_rd_bufa(self, _W, len), _fwaq_v32i8(ach));		/* complement */
		#else
		v32i8_t ach = _loadu_v32i8(_rev(pos + (len - 1)));
		// _print_v32i8(ach); _print_v32i8(_fwaq_v32i8(ach));
		_storeu_v32i8(_rd_bufa(self, _W, len), _fwaq_v32i8(ach));	/* complement; will invade bufa[BLK..BLK+_W] */
		#endif
	}
	return;
}

/**
 * @fn fill_fetch_seq_a_n
 * FIXME: this function invades 16bytes before bufa.
 */
static _force_inline
void fill_fetch_seq_a_n(
	struct gaba_dp_context_s *self,
	uint64_t ofs,
	uint8_t const *pos,
	uint64_t len)
{
	#ifdef DEBUG
	if(len > _W + BLK) {
		fprintf(stderr, "invalid len(%lu)\n", len);
		*((volatile uint8_t *)NULL);
	}
	#endif

	if(pos < GABA_EOU) {
		debug("reverse fetch a: pos(%p), len(%lu)", pos, len);
		/* reverse fetch: 2 * alen - (2 * alen - pos) + (len - 32) */
		pos += len; ofs += len;		/* fetch in reverse direction */
		while(len > 0) {
			uint64_t l = MIN2(len, 16);
			#ifdef UNSAFE_FETCH
			v16i8_t ach = _loadu_v16i8(pos - 16);
			_storeu_v16i8(_rd_bufa(self, ofs - l, l), _rvaq_v16i8(ach));/* reverse */
			#else
			v16i8_t ach = _loadu_v16i8(pos - l);						/* fetch in the reverse order */
			_storeu_v16i8(_rd_bufa(self, ofs - l, l), _rvaq_v16i8(ach, l));	/* reverse; this will invade bufa[BLK..BLK+_W] */
			#endif
			len -= l; pos -= l; ofs -= l;
		}
	} else {
		debug("forward fetch a: pos(%p), len(%lu), p(%p)", pos, len, _rev(pos + len - 1));
		/* forward fetch: 2 * alen - pos */
		pos += len - 1; ofs += len;
		while(len > 0) {
			uint64_t l = MIN2(len, 16);
			#ifdef UNSAFE_FETCH
			v16i8_t ach = _loadu_v16i8(_rev(pos));
			_storeu_v16i8(_rd_bufa(self, ofs - l, l), _fwaq_v16i8(ach));/* complement */
			#else
			v16i8_t ach = _loadu_v16i8(_rev(pos));					/* fetch in the forward order */
			// _print_v16i8(ach); _print_v16i8(_fwaq_v16i8(ach, l));
			_storeu_v16i8(_rd_bufa(self, ofs - l, l), _fwaq_v16i8(ach, l));	/* complement; will invade bufa[BLK..BLK+_W] */
			#endif
			len -= l; pos -= l; ofs -= l;
		}
	}
	return;
}

/**
 * @fn fill_fetch_seq_b
 */
static _force_inline
void fill_fetch_seq_b(
	struct gaba_dp_context_s *self,
	uint8_t const *pos,
	uint64_t len)
{
	if(pos < GABA_EOU) {
		debug("forward fetch b: pos(%p), len(%lu)", pos, len);
		/* forward fetch: pos */
		v32i8_t bch = _loadu_v32i8(pos);
		_storeu_v32i8(_rd_bufb(self, _W, len), _fwbq_v32i8(bch));		/* forward; will not invade any buf */
	} else {
		debug("reverse fetch b: pos(%p), len(%lu), p(%p)", pos, len, _rev(pos + len - 1));
		/* reverse fetch: 2 * blen - pos + (len - 32) */
		#ifdef UNSAFE_FETCH
		v32i8_t bch = _loadu_v32i8(_rev(pos) - (BLK - 1));
		_storeu_v32i8(_rd_bufb(self, _W, len), _rvbq_v32i8(bch));		/* reverse complement */
		#else
		v32i8_t bch = _loadu_v32i8(_rev(pos + (len - 1)));
		// _print_v32i8(bch); _print_v32i8(_rvbq_v32i8(bch));
		_storeu_v32i8(_rd_bufb(self, _W + len - BLK, BLK), _rvbq_v32i8(bch));	/* reverse complement; not to use swapn for v32i8_t, will invade bufb[0.._W] and bufa */
		#endif
	}
	return;
}

/**
 * @fn fill_fetch_seq_b_n
 * FIXME: this function invades 16bytes after bufb.
 */
static _force_inline
void fill_fetch_seq_b_n(
	struct gaba_dp_context_s *self,
	uint64_t ofs,
	uint8_t const *pos,
	uint64_t len)
{
	#ifdef DEBUG
	if(len > _W + BLK) {
		fprintf(stderr, "invalid len(%lu)\n", len);
		*((volatile uint8_t *)NULL);
	}
	#endif

	if(pos < GABA_EOU) {
		debug("forward fetch b: pos(%p), len(%lu)", pos, len);
		/* forward fetch: pos */
		while(len > 0) {
			uint64_t l = MIN2(len, 16);									/* advance length */
			v16i8_t bch = _loadu_v16i8(pos);
			#ifdef UNSAFE_FETCH
			_storeu_v16i8(_rd_bufb(self, ofs, l), _fwbq_v16i8(bch));	/* FIXME: will invade a region after bufb, will break arlim..bid */
			#else
			_storeu_v16i8(_rd_bufb(self, ofs, l), _fwbq_v16i8(bch, l));	/* FIXME: will invade a region after bufb, will break arlim..bid */
			#endif
			len -= l; pos += l; ofs += l;
		}
	} else {
		debug("reverse fetch b: pos(%p), len(%lu), p(%p)", pos, len, _rev(pos + len - 1));
		/* reverse fetch: 2 * blen - pos + (len - 16) */
		while(len > 0) {
			uint64_t l = MIN2(len, 16);									/* advance length */
			#ifdef UNSAFE_FETCH
			v16i8_t bch = _loadu_v16i8(_rev(pos + (16 - 1)));
			_storeu_v16i8(_rd_bufb(self, ofs, l), _rvbq_v16i8(bch));
			#else
			v16i8_t bch = _loadu_v16i8(_rev(pos + (l - 1)));			/* reverse fetch */
			// _print_v16i8(bch); _print_v16i8(_rvbq_v16i8(bch, l));
			_storeu_v16i8(_rd_bufb(self, ofs, l), _rvbq_v16i8(bch, l));	/* FIXME: will invade a region after bufb */
			#endif
			len -= l; pos += l; ofs += l;
		}
	}
	return;
}

/**
 * @fn fill_fetch_core
 * @brief fetch sequences from current sections, lengths must be shorter than 32 (= BLK)
 */
static _force_inline
void fill_fetch_core(
	struct gaba_dp_context_s *self,
	uint32_t acnt,
	uint32_t alen,
	uint32_t bcnt,
	uint32_t blen)
{
	/* fetch seq a */
	nvec_t a = _loadu_n(_rd_bufa(self, acnt, _W));		/* unaligned */
	fill_fetch_seq_a(self, self->w.r.atptr - self->w.r.arem, alen);		/* will not invade bufb */
	_store_n(_rd_bufa(self, 0, _W), a);					/* always aligned */

	/* fetch seq b */
	nvec_t b = _loadu_n(_rd_bufb(self, bcnt, _W));		/* unaligned */
	fill_fetch_seq_b(self, self->w.r.btptr - self->w.r.brem, blen);		/* will invade bufa */
	_store_n(_rd_bufb(self, 0, _W), b);					/* always aligned */

	_print_n(a); _print_n(b);
	return;
}

/**
 * @fn fill_cap_fetch
 */
static _force_inline
void fill_cap_fetch(
	struct gaba_dp_context_s *self,
	struct gaba_block_s const *blk)
{
	/* fetch: len might be clipped by ridx */
	v2i32_t const rem = _load_v2i32(&self->w.r.arem);
	v2i32_t const lim = _set_v2i32(BLK);
	v2i32_t len = _min_v2i32(rem, lim);
	_print_v2i32(len);
	fill_fetch_core(self, (blk - 1)->acnt, _lo32(len), (blk - 1)->bcnt, _hi32(len));
	return;
}

/**
 * @fn fill_init_fetch
 * @brief similar to cap fetch, updating ridx and rem
 */
static _force_inline
int64_t fill_init_fetch(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk,
	v2i64_t pos)
{
	_test_bar(head); _test_bar(mid); _test_bar(tail);
	/* calc irem (length until the init_fetch breakpoint; restore remaining head margin lengths) */
	v2i32_t const adj = _seta_v2i32(0, 1);
	v2i32_t irem = _sub_v2i32(
		_seta_v2i32(INIT_FETCH_BPOS, INIT_FETCH_APOS),
		_cvt_v2i64_v2i32(pos)
	);

	/* load srem (length until the next breakpoint): no need to update srem since advance counters are always zero */
	v2i32_t srem = _load_v2i32(&self->w.r.arem);

	/* fetch, bounded by remaining sequence lengths and remaining head margins */
	v2i32_t len = _min_v2i32(
		_min_v2i32(
			irem,								/* if remaining head margin is the shortest */
			srem								/* if remaining sequence length is the shortest */
		),
		_add_v2i32(								/* if the opposite sequence length is the shortest */
			_swap_v2i32(_sub_v2i32(srem, irem)),
			_add_v2i32(adj, irem)
		)
	);
	debug("pos(%ld, %ld)", (int64_t)_hi64(pos), (int64_t)_lo64(pos));
	_print_v2i32(irem);
	_print_v2i32(srem);
	_print_v2i32(len);

	/* fetch sequence and store at (0, 0), then save newly loaded sequence lengths */
	fill_fetch_core(self, 0, _lo32(len), 0, _hi32(len));
	_print_n(_loadu_n(_rd_bufa(self, _lo32(len), _W)));
	_print_n(_loadu_n(_rd_bufb(self, _hi32(len), _W)));

	/* save fetch length for use in the next block fill / tail construction */
	_store_v2i8(&blk->acnt, _cvt_v2i32_v2i8(len));
	_store_v2i32(&self->w.r.arem, _sub_v2i32(srem, len));

	return(_hi64(pos) + _hi32(len));
}

/**
 * @fn fill_restore_fetch
 * @brief fetch sequence from an existing block
 */
static _force_inline
void fill_restore_fetch(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail,
	struct gaba_block_s const *blk,
	v2i32_t ridx)
{
	/* load segment head info */
	v2i32_t sridx = _add_v2i32(_load_v2i32(&tail->aridx), _load_v2i32(&tail->aadv));
	struct gaba_joint_tail_s const *prev_tail = tail->tail;
	_print_v2i32(sridx);

	/* calc fetch positions and lengths */
	v2i32_t dridx = _add_v2i32(ridx, _set_v2i32(_W));	/* desired fetch pos */
	v2i32_t cridx = _min_v2i32(dridx, sridx);			/* might be clipped at the head of sequence section */
	v2i32_t ofs = _sub_v2i32(dridx, cridx);				/* lengths fetched from phantom */
	v2i32_t len = _min_v2i32(
		cridx,											/* clipped at the tail of sequence section */
		_sub_v2i32(_set_v2i32(_W + BLK), ofs)			/* lengths to fill buffers */
	);
	_print_v2i32(dridx); _print_v2i32(cridx);
	_print_v2i32(ofs); _print_v2i32(len);

	/* init buffer with magic (for debugging) */
	_memset_blk_a(self->w.r.bufa, 0, 2 * (BW_MAX + BLK));

	/* fetch seq a */
	fill_fetch_seq_a_n(self, _lo32(ofs), tail->atptr - _lo32(cridx), _lo32(len));
	if(_lo32(ofs) > 0) {
		nvec_t ach = _and_n(_set_n(0x0f), _loadu_n(&prev_tail->ch));/* aligned to 16byte boundaries */
		_print_n(ach);
		_storeu_n(_rd_bufa(self, 0, _lo32(ofs)), ach);				/* invades backward */
	}
	_print_n(_loadu_n(_rd_bufa(self, 0, _W)));

	/* fetch seq b */
	if(_hi32(ofs) > 0) {
		nvec_t bch = _and_n(
			_set_n(0x0f),
			_shr_n(_loadu_n(&prev_tail->ch.w[_W - _hi32(ofs)]), 4)	/* invades backward */
		);
		_print_n(bch);
		_storeu_n(_rd_bufb(self, 0, _hi32(ofs)), bch);				/* aligned store */
	}
	fill_fetch_seq_b_n(self, _hi32(ofs), tail->btptr - _hi32(cridx), _hi32(len));
	_print_n(_loadu_n(_rd_bufb(self, 0, _W)));
	return;
}

/*
 * @fn fill_load_section
 */
static _force_inline
void fill_load_section(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail,
	v2i32_t id, v2i32_t len, v2i64_t bptr,
	uint32_t pridx)
{
	_test_bar(head); _test_bar(mid); _test_bar(tail);
	/* calc ridx */
	v2i32_t ridx = _load_v2i32(&tail->aridx);	/* (bridx, aridx), 0 for reload */
	ridx = _sel_v2i32(_eq_v2i32(ridx, _zero_v2i32()),/* if ridx is zero (occurs when section is updated) */
		len,									/* load the next section */
		ridx									/* otherwise use the same section as the previous */
	);

	/* calc #bases until the next breakpoint (FIXME: simpler calculation needed) */
	v2i32_t brem = _seta_v2i32(tzcnt(tail->bbrk) + 1, tzcnt(tail->abrk) + 1);
	v2i32_t rem = _sel_v2i32(_gt_v2i32(brem, _set_v2i32(_W)),
		ridx,									/* no merging breakpoint found for this extension */
		_min_v2i32(ridx, brem)					/* merging breakpoint found, bounded by brem */
	);
	v2i32_t rlim = _sub_v2i32(ridx, rem);		/* remaining bases after the breakpoint (rlim == 0 if no breakpoint found) */
	v2i64_t tptr = _add_v2i64(bptr,				/* calc tail pointer: tptr = bptr + rem */
		_cvt_v2i32_v2i64(_sub_v2i32(len, rlim))
	);
	_print_v2i32(brem); _print_v2i32(rem); _print_v2i32(rlim);
	_print_v2i64(tptr);

	/* save all */
	_store_v2i32(&self->w.r.arlim, rlim);
	_store_v2i32(&self->w.r.aid, id);
	_store_v2i64(&self->w.r.atptr, tptr);
	_storeu_u64(&self->w.r.pridx, pridx);		/* (pridx, ofsd) = (remaining p length, 0) */
	_store_v2i32(&self->w.r.arem, rem);
	_store_v2i32(&self->w.r.asridx, ridx);
	_print_v2i32(ridx);

	/* save tail */
	self->w.r.tail = tail;
	return;
}

/**
 * @fn fill_create_phantom
 * @brief initialize phantom block
 */
static _force_inline
struct gaba_block_s *fill_create_phantom(
	struct gaba_dp_context_s *self,
	struct gaba_block_s const *prev_blk,
	uint32_t cnt)								/* (0, 0) for the first phantom, (prev_blk->bcnt, prev_blk->acnt) for middle phantom */
{
	struct gaba_phantom_s *ph = (struct gaba_phantom_s *)self->stack.top;
	debug("start stack_top(%p), stack_end(%p)", self->stack.top, self->stack.end);

	/* head sequence buffers are already filled, continue to body fill-in (first copy phantom block) */
	_memcpy_blk_uu(&ph->diff, &prev_blk->diff, sizeof(struct gaba_diff_vec_s));
	ph->acc = prev_blk->acc;					/* just copy */
	ph->xstat = (prev_blk->xstat & ROOT) | HEAD;/* propagate root-update flag and mark head */
	_store_v2i8(&ph->acnt, cnt);				/* FIXED! */
	ph->reserved = 0;							/* overlaps with mask */
	ph->blk = prev_blk;
	debug("ph(%p), xstat(%x), prev_blk(%p), blk(%p)", ph, ph->xstat, prev_blk, (struct gaba_block_s *)(ph + 1) - 1);
	return((struct gaba_block_s *)(ph + 1) - 1);
}

/**
 * @fn fill_create_bridge
 * @brief ``bridge'' tail skips <len - ridx> bases at the head of the sections.
 */
static _force_inline
struct gaba_joint_tail_s *fill_create_bridge(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *prev_tail,
	v2i32_t id, v2i32_t len, v2i64_t bptr,
	v2i32_t adv)
{
	// struct gaba_joint_tail_s *tail = _bridge(gaba_dp_malloc(self, BRIDGE_TAIL_SIZE));
	struct gaba_joint_tail_s *tail = gaba_dp_malloc(self, sizeof(struct gaba_joint_tail_s));
	debug("create bridge, p(%p)", tail);

	/* copy ch, xd, and md */
	_storeu_n(tail->ch.w, _loadu_n(prev_tail->ch.w));
	_storeu_n(tail->xd.drop, _loadu_n(prev_tail->xd.drop));
	_storeu_w(tail->md.delta, _loadu_w(prev_tail->md.delta));

	/* vector position and scores are unchanged */
	_storeu_u64(&tail->mdrop, 0x00010000 | _loadu_u64(&prev_tail->mdrop));	/* mark istat, copy mdrop and pridx */

	/* relative section pointer is moved, set advanced length */
	_store_v2i32(&tail->aridx, _sub_v2i32(len, adv));			/* FIXME: len -> correct ridx */
	_store_v2i32(&tail->aadv, adv);

	/* copy breakpoint arrays and pointers */
	tail->tail = prev_tail;
	_store_v2i64(&tail->abrk, _load_v2i64(&prev_tail->abrk));	/* unchanged */
	_store_v2i64(&tail->atptr, _add_v2i64(bptr, _cvt_v2i32_v2i64(len)));/* correct tptr for safety (actually not required) */
	_store_v2i32(&tail->f.aid, id);								/* correct id pair is required */
	_memcpy_blk_uu(&tail->f.ascnt, &prev_tail->f.ascnt, 32);	/* just copy (unchanged) */
	tail->f.status = prev_tail->f.status;
	return(tail);
}

/**
 * @fn fill_load_vectors
 * @brief load sequences and indices from the previous tail
 */
static _force_inline
struct gaba_block_s *fill_load_vectors(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail)
{
	/* load sequence vectors */
	nvec_t const mask = _set_n(0x0f);
	nvec_t ch = _loadu_n(&tail->ch.w);
	nvec_t ach = _and_n(mask, ch);
	nvec_t bch = _and_n(mask, _shr_n(ch, 4));	/* bit 7 must be cleared not to affect shuffle mask */
	_store_n(_rd_bufa(self, 0, _W), ach);
	_store_n(_rd_bufb(self, 0, _W), bch);
	_print_n(ach); _print_n(bch);

	/* copy max and middle delta vectors */
	nvec_t xd = _loadu_n(&tail->xd);
	wvec_t md = _loadu_w(&tail->md);
	_store_n(&self->w.r.xd, xd);
	_store_w(&self->w.r.md, md);
	_print_n(xd);
	_print_w(md);

	/* extract the last block pointer, pass to fill-in loop */
	return(fill_create_phantom(self, _last_block(tail), 0));
}

/**
 * @fn fill_create_tail
 * @brief create joint_tail at the end of the blocks
 */
static _force_inline
int32_t fill_save_vectors(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s *tail,
	uint32_t cnt)
{
	/* store char vector */
	nvec_t ach = _loadu_n(_rd_bufa(self, cnt & 0xff, _W));
	nvec_t bch = _loadu_n(_rd_bufb(self, cnt>>8, _W));
	_print_n(ach); _print_n(bch);
	_storeu_n(&tail->ch, _or_n(ach, _shl_n(bch, 4)));

	/* copy delta vectors */
	nvec_t xd = _load_n(&self->w.r.xd);				/* gcc-4.8 w/ -mavx2 -mbmi2 has an optimization bug with vpmovsxbw */
	wvec_t md = _load_w(&self->w.r.md);
	_storeu_n(&tail->xd, xd);
	_storeu_w(&tail->md, md);
	_print_n(xd);
	_print_w(md);

	/* search max section */
	md = _add_w(md, _cvt_n_w(xd));					/* xd holds drop from max */
	return(_hmax_w(md));
}
static _force_inline
struct gaba_joint_tail_s *fill_save_section(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s *tail,
	uint32_t xstat,
	int32_t mdrop)
{
	v2i32_t rlim = _load_v2i32(&self->w.r.arlim);
	v2i32_t id = _load_v2i32(&self->w.r.aid);
	v2i64_t tptr = _load_v2i64(&self->w.r.atptr);
	tptr = _add_v2i64(tptr, _cvt_v2i32_v2i64(rlim));
	_print_v2i32(rlim); _print_v2i32(_load_v2i32(&self->w.r.arem));

	/* adjust reversed indices */
	tail->mdrop = mdrop;
	tail->istat = 0;
	tail->pridx = self->w.r.pridx;
	v2i32_t ridx = _add_v2i32(_load_v2i32(&self->w.r.arem), rlim);
	v2i32_t sridx = _load_v2i32(&self->w.r.asridx);
	v2i32_t adv = _sub_v2i32(sridx, ridx);
	_store_v2i32(&tail->aridx, ridx);
	_store_v2i32(&tail->aadv, adv);

	/* adjust breakpoint masks */
	struct gaba_joint_tail_s const *prev_tail = self->w.r.tail;
	v2i64_t brk = _shrv_v2i64(_loadu_v2i64(&prev_tail->abrk), _cvt_v2i32_v2i64(adv));
	tail->tail = prev_tail;
	_store_v2i64(&tail->abrk, brk);

	/* save section info */
	_store_v2i64(&tail->atptr, tptr);
	_store_v2i32(&tail->f.aid, id);
	_print_v2i32(ridx); _print_v2i32(sridx); _print_v2i32(adv);

	/* calc end-of-section flag, section counts, and base counts */
	v2i32_t update = _eq_v2i32(ridx, _zero_v2i32());
	_store_v2i32(&tail->f.ascnt, _sub_v2i32(
		_load_v2i32(&prev_tail->f.ascnt),			/* aligned section counts */
		update										/* increment by one if the pointer reached the end of the current section */
	));
	_store_v2i64(&tail->f.apos, _add_v2i64(
		_load_v2i64(&prev_tail->f.apos),
		_cvt_v2i32_v2i64(adv)
	));
	_print_v2i32(_load_v2i32(&tail->f.ascnt));
	_print_v2i64(_load_v2i64(&tail->f.apos));

	/* store max, status flag */
	tail->f.max = _offset(prev_tail) + self->w.r.ofsd + mdrop;
	tail->f.status = ((xstat & (TERM | CONT))<<8) | _mask_v2i32(update);
	debug("prev_offset(%ld), offset(%ld), max(%d, %ld)",
		_offset(prev_tail), _offset(prev_tail) + self->w.r.ofsd, mdrop, tail->f.max);
	return(tail);
}
static _force_inline
struct gaba_joint_tail_s *fill_create_tail(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	/* inspect fetched base counts */
	uint32_t cnt = _load_v2i8(&blk->acnt), xstat = blk->xstat;/* *((uint16_t const *)&blk->acnt); */

	/* create joint_tail: squash the last block if no vector was filled */
	struct gaba_joint_tail_s *tail = (struct gaba_joint_tail_s *)(blk + (cnt != 0));
	self->stack.top = (void *)(tail + 1);				/* write back stack_top */
	debug("end stack_top(%p), stack_end(%p), blk(%p)", self->stack.top, self->stack.end, blk);

	return(fill_save_section(self, tail, xstat,			/* sections */
		fill_save_vectors(self, tail, cnt)				/* vectors */
	));
}

/**
 * @macro _fill_load_context
 * @brief load vectors onto registers
 */
#if MODEL == LINEAR
#define _fill_load_context(_blk) \
	debug("blk(%p)", (_blk)); \
	/* load sequence buffer offset */ \
	uint8_t const *aptr = _rd_bufa(self, 0, _W); \
	uint8_t const *bptr = _rd_bufb(self, 0, _W); \
	/* load mask pointer */ \
	struct gaba_mask_pair_s *ptr = ((struct gaba_block_s *)(_blk))->mask; \
	/* load vector registers */ \
	register nvec_t dh = _loadu_n(((_blk) - 1)->diff.dh); \
	register nvec_t dv = _loadu_n(((_blk) - 1)->diff.dv); \
	_print_n(_add_n(dh, _load_ofsh(self->scv))); \
	_print_n(_add_n(dv, _load_ofsv(self->scv))); \
	/* load delta vectors */ \
	register nvec_t delta = _zero_n(); \
	register nvec_t drop = _load_n(self->w.r.xd.drop); \
	_print_n(drop); \
	_print_w(_add_w(_load_w(&self->w.r.md), _add_w(_cvt_n_w(delta), _set_w(_offset(self->w.r.tail) + self->w.r.ofsd - 128)))); \
	_print_w(_add_w(_add_w(_load_w(&self->w.r.md), _cvt_n_w(delta)), _add_w(_cvt_n_w(drop), _set_w(_offset(self->w.r.tail) + self->w.r.ofsd)))); \
	/* load direction determiner */ \
	struct gaba_dir_s dir = _dir_init((_blk) - 1);
#else	/* AFFINE and COMBINED */
#define _fill_load_context(_blk) \
	debug("blk(%p)", (_blk)); \
	/* load sequence buffer offset */ \
	uint8_t const *aptr = _rd_bufa(self, 0, _W); \
	uint8_t const *bptr = _rd_bufb(self, 0, _W); \
	/* load mask pointer */ \
	struct gaba_mask_pair_s *ptr = ((struct gaba_block_s *)(_blk))->mask; \
	/* load vector registers */ \
	register nvec_t dh = _loadu_n(((_blk) - 1)->diff.dh); \
	register nvec_t dv = _loadu_n(((_blk) - 1)->diff.dv); \
	register nvec_t de = _loadu_n(((_blk) - 1)->diff.de); \
	register nvec_t df = _loadu_n(((_blk) - 1)->diff.df); \
	_print_n(_sub_n(_load_ofsh(self->scv), dh)); \
	_print_n(_add_n(dv, _load_ofsv(self->scv))); \
	_print_n(_sub_n(_sub_n(de, dv), _load_adjh(self->scv))); \
	_print_n(_sub_n(_add_n(df, dh), _load_adjv(self->scv))); \
	/* load delta vectors */ \
	register nvec_t delta = _zero_n(); \
	register nvec_t drop = _load_n(self->w.r.xd.drop); \
	_print_n(drop); \
	_print_w(_add_w(_load_w(&self->w.r.md), _add_w(_cvt_n_w(delta), _set_w(_offset(self->w.r.tail) + self->w.r.ofsd - 128)))); \
	_print_w(_add_w(_add_w(_load_w(&self->w.r.md), _cvt_n_w(delta)), _add_w(_cvt_n_w(drop), _set_w(_offset(self->w.r.tail) + self->w.r.ofsd)))); \
	/* load direction determiner */ \
	struct gaba_dir_s dir = _dir_init((_blk) - 1);
#endif

/**
 * @macro _fill_body
 * @brief update vectors
 */
#if MODEL == LINEAR
#define _fill_body() { \
	register nvec_t t = _match_n(_loadu_n(aptr), _loadu_n(bptr)); \
	_print_n(_loadu_n(aptr)); _print_n(_loadu_n(bptr)); \
	t = _shuf_n(_load_sb(self->scv), t); _print_n(t); \
	t = _max_n(dh, t); \
	t = _max_n(dv, t); \
	ptr->h.mask = _mask_n(_eq_n(t, dv)); \
	ptr->v.mask = _mask_n(_eq_n(t, dh)); \
	debug("mask(%lx, %lx)", (uint64_t)ptr->h.all, (uint64_t)ptr->v.all); \
	ptr++; \
	nvec_t _dv = _sub_n(t, dh); \
	dh = _sub_n(t, dv); \
	dv = _dv; \
	_print_n(drop); \
	_print_n(_add_n(dh, _load_ofsh(self->scv))); \
	_print_n(_add_n(dv, _load_ofsv(self->scv))); \
}
#elif MODEL == AFFINE
#define _fill_body() { \
	register nvec_t t = _match_n(_loadu_n(aptr), _loadu_n(bptr)); \
	_print_n(_loadu_n(aptr)); _print_n(_loadu_n(bptr)); \
	t = _shuf_n(_load_sb(self->scv), t); _print_n(t); \
	t = _max_n(de, t); \
	t = _max_n(df, t); \
	ptr->h.mask = _mask_n(_eq_n(t, de)); \
	ptr->v.mask = _mask_n(_eq_n(t, df)); \
	/* update de and dh */ \
	de = _add_n(de, _load_adjh(self->scv)); \
	nvec_t te = _max_n(de, t); \
	ptr->e.mask = _mask_n(_eq_n(te, t)); \
	de = _add_n(te, dh); \
	dh = _add_n(dh, t); \
	/* update df and dv */ \
	df = _add_n(df, _load_adjv(self->scv)); \
	nvec_t tf = _max_n(df, t); \
	ptr->f.mask = _mask_n(_eq_n(tf, t)); \
	debug("mask(%lx, %lx, %lx, %lx)", (uint64_t)ptr->h.all, (uint64_t)ptr->v.all, (uint64_t)ptr->e.all, (uint64_t)ptr->f.all); \
	df = _sub_n(tf, dv); \
	t = _sub_n(dv, t); \
	ptr++; dv = dh; dh = t; \
	_print_n(_sub_n(_load_ofsh(self->scv), dh)); \
	_print_n(_add_n(dv, _load_ofsv(self->scv))); \
	_print_n(_sub_n(_sub_n(de, dv), _load_adjh(self->scv))); \
	_print_n(_sub_n(_add_n(df, dh), _load_adjv(self->scv))); \
}
#else /* MODEL == COMBINED */
#define _fill_body() { \
	register nvec_t t = _match_n(_loadu_n(aptr), _loadu_n(bptr)); \
	register nvec_t dfh = _add_n(dv, _load_gfh(self->scv)); \
	register nvec_t dfv = _sub_n(_load_gfv(self->scv), dh); \
	_print_n(_sub_n(_zero_n(), dh)); _print_n(dv); _print_n(de); _print_n(df); \
	_print_n(dfv); _print_n(dfh); \
	_print_n(_loadu_n(aptr)); _print_n(_loadu_n(bptr)); \
	register nvec_t s = _max_n(de, df); \
	t = _shuf_n(_load_sb(self->scv), t); _print_n(t); \
	s = _max_n(s, dfh); \
	t = _max_n(t, dfv); \
	t = _max_n(t, s); \
	_print_n(t); \
	uint64_t mask_gfh = _mask_u64(_mask_n(_eq_n(t, dfh))), mask_gh = _mask_u64(_mask_n(_eq_n(t, de))); \
	uint64_t mask_gfv = _mask_u64(_mask_n(_eq_n(t, dfv))), mask_gv = _mask_u64(_mask_n(_eq_n(t, df))); \
	debug("mask_gfh(%lx), mask_gh(%lx), mask_gfv(%lx), mask_gv(%lx)", mask_gfh, mask_gh, mask_gfv, mask_gv); \
	ptr->h.all = mask_gfh | mask_gh; mask_gh &= ~mask_gfh; \
	ptr->v.all = mask_gfv | mask_gv; mask_gv &= ~mask_gfv; \
	/* update de and dh */ \
	de = _add_n(de, _load_adjh(self->scv)); \
	nvec_t te = _max_n(de, t); \
	ptr->e.all = mask_gh | _mask_u64(_mask_n(_eq_n(te, t))); \
	de = _add_n(te, dh); \
	dh = _add_n(dh, t); \
	/* update df and dv */ \
	df = _add_n(df, _load_adjv(self->scv)); \
	nvec_t tf = _max_n(df, t); \
	ptr->f.all = mask_gv | _mask_u64(_mask_n(_eq_n(tf, t))); \
	debug("mask_ge(%lx), mask_gf(%lx), mask(%lx, %lx, %lx, %lx)", (uint64_t)_mask_u64(_mask_n(_eq_n(te, t))), (uint64_t)_mask_u64(_mask_n(_eq_n(tf, t))), (uint64_t)ptr->h.all, (uint64_t)ptr->v.all, (uint64_t)ptr->e.all, (uint64_t)ptr->f.all); \
	df = _sub_n(tf, dv); \
	t = _sub_n(dv, t); \
	ptr++; dv = dh; dh = t; \
	_print_n(_sub_n(_load_ofsh(self->scv), dh)); \
	_print_n(_add_n(dv, _load_ofsv(self->scv))); \
	_print_n(_sub_n(_sub_n(de, dv), _load_adjh(self->scv))); \
	_print_n(_sub_n(_add_n(df, dh), _load_adjv(self->scv))); \
}
#endif /* MODEL */

/**
 * @macro _fill_update_delta
 * @brief update small delta vector and max vector
 */
#define _fill_update_delta(_op_add, _vector, _ofs) { \
	nvec_t _t = _op_add(_ofs, _vector); \
	delta = _add_n(delta, _t); \
	drop = _subs_n(drop, _t); \
	_dir_update(dir, _t); \
	_print_n(delta); _print_n(drop); \
	_print_w(_add_w(_load_w(&self->w.r.md), _add_w(_cvt_n_w(delta), _set_w(_offset(self->w.r.tail) + self->w.r.ofsd - 128)))); \
	_print_w(_add_w(_add_w(_load_w(&self->w.r.md), _cvt_n_w(delta)), _add_w(_cvt_n_w(drop), _set_w(_offset(self->w.r.tail) + self->w.r.ofsd)))); \
}
/**
 * @macro _fill_right, _fill_down
 * @brief wrapper of _fill_body and _fill_update_delta
 */
#define _fill_right_update_ptr() { \
	aptr--;				/* increment sequence buffer pointer */ \
}
#define _fill_right_windback_ptr() { \
	aptr++; \
}
#if MODEL == LINEAR
#define _fill_right() { \
	dh = _bsl_n(dh, 1);	/* shift left dh */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_add_n, dh, _load_ofsh(self->scv)); \
}
#else	/* AFFINE and COMBINED */
#define _fill_right() { \
	dh = _bsl_n(dh, 1);	/* shift left dh */ \
	df = _bsl_n(df, 1);	/* shift left df */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_sub_n, dh, _load_ofsh(self->scv)); \
}
#endif /* MODEL */
#define _fill_down_update_ptr() { \
	bptr++;				/* increment sequence buffer pointer */ \
}
#define _fill_down_windback_ptr() { \
	bptr--; \
}
#if MODEL == LINEAR
#define _fill_down() { \
	dv = _bsr_n(dv, 1);	/* shift right dv */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_add_n, dv, _load_ofsv(self->scv)); \
}
#else	/* AFFINE and COMBINED */
#define _fill_down() { \
	dv = _bsr_n(dv, 1);	/* shift right dv */ \
	de = _bsr_n(de, 1);	/* shift right de */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_add_n, dv, _load_ofsv(self->scv)); \
}
#endif /* MODEL */

#if defined(DEBUG) || defined(DEBUG_OVERFLOW)
#define _check_overflow(_delta, _drop) { \
	int8_t b[_W], d[_W], flag = 0; int16_t ovf[_W], udf[_W], m1[_W], m2[_W], m3[_W]; \
	_storeu_n(b, _delta); _storeu_n(d, _drop); \
	wvec_t ov = _and_w(_set_w(0x0100), _cvt_n_w(_andn_n(_add_n(drop, delta), _and_n(drop, delta)))); \
	wvec_t uv = _and_w(_set_w(0x0100), _cvt_n_w(_or_n(_subs_n(delta, _set_n(0x40)), drop))); \
	_storeu_w(ovf, ov); _storeu_w(udf, uv); _storeu_w(m1, md); \
	wvec_t md2v = _add_w(_add_w(_add_w(md, ov), uv), _set_w(-cofs - 0x100)); _storeu_w(m2, md2v); \
	wvec_t mds1v = _add_w(md2v, _cvt_n_w(drop)); _storeu_w(m3, mds1v); \
	for(uint64_t i = 0; i < _W - 1; i++) { if(b[i + 1] > b[i] + 128) { flag |= 1; } if(b[i + 1] < b[i] - 128) { flag |= 1; } } \
	for(uint64_t i = 0; i < _W - 1; i++) { if(m2[i + 1] > m2[i] + 128) { flag |= 2; } if(m2[i + 1] < m2[i] - 128) { flag |= 2; } } \
	if(flag != 0) { \
		fprintf(stderr, "overflow detected, flag(%x)\n", flag); \
		fprintf(stderr, "delta("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%04d, ", b[i]); } fprintf(stderr, ")\n"); \
		fprintf(stderr, " drop("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%04d, ", d[i]); } fprintf(stderr, ")\n"); \
		fprintf(stderr, "   md("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%04d, ", self->w.r.md.delta[i]); } fprintf(stderr, ")\n"); \
		fprintf(stderr, "   m1("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%04d, ", m1[i]); } fprintf(stderr, ")\n"); \
		fprintf(stderr, "  ovf("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%04d, ", ovf[i]); } fprintf(stderr, ")\n"); \
		fprintf(stderr, "  udf("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%04d, ", udf[i]); } fprintf(stderr, ")\n"); \
		fprintf(stderr, "   m2("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%04d, ", m2[i]); } fprintf(stderr, ")\n"); \
		fprintf(stderr, "   m3("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%04d, ", m3[i]); } fprintf(stderr, ")\n"); \
		fprintf(stderr, "  sum("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%04d, ", (int8_t)(b[i] + d[i])); } fprintf(stderr, ")\n"); \
		fprintf(stderr, " mask("); for(uint64_t i = 0; i < _W; i++) { fprintf(stderr, "%c,    ", (int8_t)(~(b[i] + d[i]) & b[i]) < 0 ? '1' : '0'); } fprintf(stderr, ")\n"); \
	} \
}
#else
#define _check_overflow(_x, _y) {}
#endif

/**
 * @macro _fill_store_context
 * @brief store vectors at the end of the block
 */
#define _fill_store_context_intl(_blk) { \
	/* store direction array */ \
	_dir_save(_blk, dir); \
	/* update xdrop status and offsets */ \
	(_blk)->xstat = (self->tx - _ext_n(drop, _W/2)) & TERM; \
	int32_t cofs = _ext_n(delta, _W/2); \
	/* store cnt */ \
	int32_t acnt = _rd_bufa(self, 0, _W) - aptr; \
	int32_t bcnt = bptr - _rd_bufb(self, 0, _W); \
	(_blk)->acnt = acnt; (_blk)->bcnt = bcnt; \
	/* write back local working buffers */ \
	self->w.r.ofsd += cofs; self->w.r.arem -= acnt; self->w.r.brem -= bcnt; \
	/* update max and middle vectors in the working buffer */ \
	nvec_t prev_drop = _load_n(&self->w.r.xd); \
	_store_n(&self->w.r.xd, drop);		/* save max delta vector */ \
	_print_n(prev_drop); _print_n(_add_n(drop, delta)); \
	(_blk)->max_mask = _mask_u64(_mask_n(_gt_n(_add_n(drop, delta), prev_drop))); \
	debug("update_mask(%lx)", (uint64_t)(_blk)->max_mask); \
	/* update middle delta vector */ \
	wvec_t md = _load_w(&self->w.r.md); \
	md = _add_w(md, _cvt_n_w(delta)); \
	_check_overflow(delta, drop); \
	/* rescue overflow */ \
	md = _add_w(md, _and_w(_set_w(0x0100), _cvt_n_w(_andn_n(_add_n(drop, delta), _and_n(drop, delta))))); \
	/* rescue underflow */ \
	md = _add_w(md, _and_w(_set_w(0x0100), _cvt_n_w(_or_n(_subs_n(delta, _set_n(0x40)), drop)))); cofs += 0x0100; \
	md = _add_w(md, _set_w(-cofs));		/* fixup offset adjustment */ \
	_store_w(&self->w.r.md, md); \
	_print_w(md); \
}
#if MODEL == LINEAR
#define _fill_store_context(_blk) { \
	_storeu_n((_blk)->diff.dh, dh); _print_n(dh); \
	_storeu_n((_blk)->diff.dv, dv); _print_n(dv); \
	_fill_store_context_intl(_blk); \
}
#else	/* AFFINE and COMBINED */
#define _fill_store_context(_blk) { \
	_storeu_n((_blk)->diff.dh, dh); _print_n(dh); \
	_storeu_n((_blk)->diff.dv, dv); _print_n(dv); \
	_storeu_n((_blk)->diff.de, de); _print_n(de); \
	_storeu_n((_blk)->diff.df, df); _print_n(df); \
	_fill_store_context_intl(_blk); \
}
#endif

/**
 * @fn fill_bulk_test_idx
 * @brief returns negative if ij-bound (for the bulk fill) is invaded
 */
static _force_inline
int64_t fill_bulk_test_idx(
	struct gaba_dp_context_s const *self)
{
	#define _test(_label, _ofs)	( (uint64_t)self->w.r._label - (uint64_t)(_ofs) )
	debug("test(%ld, %ld, %ld), len(%u, %u, %u)",
		(int64_t)_test(arem, BLK), (int64_t)_test(brem, BLK), (int64_t)_test(pridx, BLK),
		self->w.r.arem, self->w.r.brem, self->w.r.pridx);
	return((int64_t)(_test(arem, BLK) | _test(brem, BLK) | _test(pridx, BLK)));
	#undef _test
}

/**
 * @fn fill_cap_test_idx
 * @brief returns negative if ij-bound (for the cap fill) is invaded
 */
#define _fill_cap_test_idx_init() \
	uint8_t const *alim = _rd_bufa(self, self->w.r.arem, _W); \
	uint8_t const *blim = _rd_bufb(self, self->w.r.brem, _W); \
	uint8_t const *plim = blim - (ptrdiff_t)alim + (ptrdiff_t)self->w.r.pridx;
#define _fill_cap_test_idx() ({ \
	debug("arem(%zd), brem(%zd), prem(%zd)", \
		(int64_t)aptr - (int64_t)alim, (int64_t)blim - (int64_t)bptr, (int64_t)plim - (int64_t)bptr + (int64_t)aptr); \
	((int64_t)aptr - (int64_t)alim) | ((int64_t)blim - (int64_t)bptr) | ((int64_t)plim - (int64_t)bptr + (int64_t)aptr); \
})

#ifndef DEBUG_ALL
#  undef DEBUG
#  undef _LOG_H_INCLUDED
#  include "log.h"
#endif

/**
 * @fn fill_bulk_block
 * @brief fill a block, returns negative value if xdrop detected
 */
static _force_inline
void fill_bulk_block(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	_test_bar(head); _test_bar(mid); _test_bar(tail);
	/* fetch sequence */
	fill_fetch_core(self, (blk - 1)->acnt, BLK, (blk - 1)->bcnt, BLK);

	/* load vectors onto registers */
	debug("blk(%p)", blk);
	_fill_load_context(blk);
	/**
	 * @macro _fill_block
	 * @brief unit unrolled fill-in loop
	 */
	#define _fill_block(_direction, _label, _jump_to) { \
		_dir_fetch(dir); \
		if(_unlikely(!_dir_is_##_direction(dir))) { \
			goto _fill_##_jump_to; \
		} \
		_fill_##_label: \
		_fill_##_direction##_update_ptr(); \
		_fill_##_direction(); \
		if(--i == 0) { break; } \
	}

	/* update diff vectors */
	int64_t i = BLK;
	while(1) {					/* 4x unrolled loop */
		_fill_block(down, d1, r1);
		_fill_block(right, r1, d2);
		_fill_block(down, d2, r2);
		_fill_block(right, r2, d1);
	}

	/* store vectors */
	self->w.r.pridx -= BLK;
	_fill_store_context(blk);
	return;
}

#ifdef REDEFINE_DEBUG
#  define DEBUG
#  undef _LOG_H_INCLUDED
#  include "log.h"
#endif

/**
 * @fn fill_bulk_k_blocks
 * @brief fill <cnt> contiguous blocks without ij-bound test
 */
static _force_inline
struct gaba_block_s *fill_bulk_k_blocks(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk,
	uint64_t cnt)
{
	/* bulk fill loop, first check termination */
	struct gaba_block_s *tblk = blk + cnt;
	while((blk->xstat | (ptrdiff_t)(tblk - blk)) > 0) {
		/* bulk fill */
		debug("blk(%p), aadv(%u), badv(%u), max(%ld)", blk + 1,
			self->w.r.asridx - self->w.r.arem - self->w.r.arlim,
			self->w.r.bsridx - self->w.r.brem - self->w.r.brlim,
			self->w.r.md.delta[_W/2] + self->w.r.xd.drop[_W/2] + _offset(self->w.r.tail) + self->w.r.ofsd);
		_print_w(_load_w(self->w.r.md.delta));
		fill_bulk_block(self, ++blk);
	}
	debug("return, blk(%p), xstat(%x), pridx(%u)", blk, blk->xstat, self->w.r.pridx);
	return(blk);
}

/**
 * @fn fill_bulk_seq_bounded
 * @brief fill blocks with ij-bound test
 */
static _force_inline
struct gaba_block_s *fill_bulk_seq_bounded(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	/* bulk fill loop, first check termination */
	while((blk->xstat | fill_bulk_test_idx(self)) >= 0) {
		debug("blk(%p), aadv(%u), badv(%u), max(%ld)", blk + 1,
			self->w.r.asridx - self->w.r.arem - self->w.r.arlim,
			self->w.r.bsridx - self->w.r.brem - self->w.r.brlim,
			self->w.r.md.delta[_W/2] + self->w.r.xd.drop[_W/2] + _offset(self->w.r.tail) + self->w.r.ofsd);
		fill_bulk_block(self, ++blk);
	}
	debug("return, blk(%p), xstat(%x)", blk, blk->xstat);
	return(blk);
}

#ifndef DEBUG_ALL
#  undef DEBUG
#  undef _LOG_H_INCLUDED
#  include "log.h"
#endif

/**
 * @fn fill_cap_seq_bounded
 * @brief fill blocks with cap test
 */
static _force_inline
struct gaba_block_s *fill_cap_seq_bounded(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	_test_bar(head); _test_bar(mid); _test_bar(tail);
	#define _fill_cap_seq_bounded_core(_dir) { \
		/* update sequence coordinate and then check term */ \
		_fill_##_dir##_update_ptr(); \
		if(_fill_cap_test_idx() < 0) { \
			_fill_##_dir##_windback_ptr(); \
			_dir_windback(dir); \
			break; \
		} \
		_fill_##_dir();		/* update band */ \
	}

	debug("blk(%p)", blk);
	while(blk->xstat >= 0) {
		/* fetch sequence */
		fill_cap_fetch(self, ++blk);
		_fill_cap_test_idx_init();
		_fill_load_context(blk);			/* contains ptr as struct gaba_mask_pair_s *ptr = blk->mask; */

		/* update diff vectors */
		struct gaba_mask_pair_s *tptr = &blk->mask[BLK];
		while(ptr < tptr) {					/* ptr is automatically incremented in _fill_right() or _fill_down() */
			_dir_fetch(dir);				/* determine direction */
			if(_dir_is_right(dir)) {
				_fill_cap_seq_bounded_core(right);
			} else {
				_fill_cap_seq_bounded_core(down);
			}
		}

		uint64_t i = ptr - blk->mask;		/* calc filled count */
		self->w.r.pridx -= i;				/* update remaining p-length */
		_dir_adjust_remainder(dir, i);		/* adjust dir remainder */
		_fill_store_context(blk);			/* store mask and vectors */
		if(_unlikely(i != BLK)) { break; }	/* reached the end */
	}
	debug("return, blk(%p), xstat(%x)", blk, blk->xstat);
	return(blk);
}

#ifdef REDEFINE_DEBUG
#  define DEBUG
#  undef _LOG_H_INCLUDED
#  include "log.h"
#endif

/**
 * @fn max_blocks_mem
 * @brief calculate maximum number of blocks to be filled (bounded by stack size)
 */
static _force_inline
uint64_t max_blocks_mem(
	struct gaba_dp_context_s const *self)
{
	uint64_t mem_size = _stack_size(&self->stack);
	uint64_t blk_cnt = mem_size / sizeof(struct gaba_block_s);
	debug("calc_max_block_mem, stack_top(%p), stack_end(%p), mem_size(%lu), cnt(%lu)",
		self->stack.top, self->stack.end, mem_size, (blk_cnt > 3) ? (blk_cnt - 3) : 0);
	return(((blk_cnt > 3) ? blk_cnt : 3) - 3);
}

/**
 * @fn max_blocks_idx
 * @brief calc max #expected blocks from remaining seq lengths,
 * used to determine #blocks which can be filled without mem boundary check.
 */
static _force_inline
uint64_t max_blocks_idx(
	struct gaba_dp_context_s const *self)
{
	uint64_t p = self->w.r.arem + self->w.r.brem;
	return(MIN2(p + p / 2, self->w.r.pridx) / BLK + 1);
}

/**
 * @fn min_blocks_idx
 * @brief calc min #expected blocks from remaining seq lengths,
 * #blocks filled without seq boundary check.
 */
static _force_inline
uint64_t min_blocks_idx(
	struct gaba_dp_context_s const *self)
{
	uint64_t p = MIN2(self->w.r.arem, self->w.r.brem);
	debug("aridx(%u), bridx(%u), p(%lu), pridx(%u)", self->w.r.arem, self->w.r.brem, p, self->w.r.pridx)
	return(MIN2(p, self->w.r.pridx) / BLK);
}

/**
 * @fn fill_seq_bounded
 * @brief fill blocks with seq bound tests (without mem test), adding head and tail
 */
static _force_inline
struct gaba_block_s *fill_seq_bounded(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	uint64_t cnt;						/* #blocks filled in bulk */
	debug("blk(%p), cnt(%lu)", blk, min_blocks_idx(self));
	while((cnt = min_blocks_idx(self)) > MIN_BULK_BLOCKS) {
		/* bulk fill without ij-bound test */
		debug("fill predetd blocks");
		if(((blk = fill_bulk_k_blocks(self, blk, cnt))->xstat & STAT_MASK) != CONT) {
			debug("term detected, xstat(%x)", blk->xstat);
			return(blk);				/* xdrop termination detected, skip cap */
		}
	}

	/* bulk fill with ij-bound test */
	if(((blk = fill_bulk_seq_bounded(self, blk))->xstat & STAT_MASK) != CONT) {
		debug("term detected, blk(%p), xstat(%x)", blk, blk->xstat);
		return(blk);					/* xdrop termination detected, skip cap */
	}

	/* cap fill (with p-bound test) */
	return(fill_cap_seq_bounded(self, blk));
}

/**
 * @fn fill_section_seq_bounded
 * @brief fill dp matrix inside section pairs
 */
static _force_inline
struct gaba_block_s *fill_section_seq_bounded(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	/* extra large bulk fill (with stack allocation) */
	uint64_t mem_cnt, seq_cnt;			/* #blocks filled in bulk */
	while((mem_cnt = max_blocks_mem(self)) < (seq_cnt = max_blocks_idx(self))) {
		debug("mem_cnt(%lu), seq_cnt(%lu)", mem_cnt, seq_cnt);

		/* here mem limit reaches first (seq sections too long), so fill without seq test */
		if((mem_cnt = MIN2(mem_cnt, min_blocks_idx(self))) > MIN_BULK_BLOCKS) {
			debug("mem bounded fill");
			if(((blk = fill_bulk_k_blocks(self, blk, mem_cnt))->xstat & STAT_MASK) != CONT) {
				return(blk);
			}
		}

		/* memory ran out: malloc a next stack and create a new phantom head */
		debug("add stack, blk(%p)", blk);
		if(gaba_dp_add_stack(self, _mem_blocks(seq_cnt)) != 0) { return(NULL); }
		blk = fill_create_phantom(self, blk, _load_v2i8(&blk->acnt));
	}

	/* bulk fill with seq bound check */
	return(fill_seq_bounded(self, blk));
}

/**
 * @fn gaba_dp_fill_root
 *
 * @brief build_root API
 */
struct gaba_fill_s *_export(gaba_dp_fill_root)(
	struct gaba_dp_context_s *self,
	struct gaba_section_s const *a,
	uint32_t apos,
	struct gaba_section_s const *b,
	uint32_t bpos,
	uint32_t pridx)
{
	/* restore dp context pointer by adding offset */
	self = _restore_dp_context(self);
	_init_bar(head); _init_bar(mid); _init_bar(tail);

	/* load current sections, then transpose sections to extract {id, len, base} pairs */
	v2i64_t asec = _loadu_v2i64(a), bsec = _loadu_v2i64(b);	/* tuple of (64bit ptr, 32-bit id, 32-bit len) */
	v2i64_t id_len = _lo_v2i64(asec, bsec), bptr = _hi_v2i64(asec, bsec);
	v2i32_t id = _cvt_v2i64_v2i32(id_len), len = _cvth_v2i64_v2i32(id_len);
	v2i32_t adv = _seta_v2i32(bpos, apos);		/* head offsets */
	_print_v2i32(id); _print_v2i32(len); _print_v2i64(bptr); _print_v2i32(adv);

	/* create bridge (skip (apos, bpos) at the head) */
	struct gaba_joint_tail_s *brg = fill_create_bridge(self, _root(self), id, len, bptr, adv);

	/* load sections */
	fill_load_section(self,
		brg, id, len, bptr,
		pridx == 0 ? UINT32_MAX : pridx			/* UINT32_MAX */
	);

	/* load sequences and extract the last block pointer */
	if(_stack_size(&self->stack) < MEM_INIT_VACANCY) {
		gaba_dp_add_stack(self, _mem_blocks(max_blocks_idx(self)));
	}
	struct gaba_block_s *blk = fill_load_vectors(self, _root(self));

	/* init fetch */
	if(fill_init_fetch(self, blk, _load_v2i64(&_root(self)->f.apos)) < INIT_FETCH_BPOS) {
		return(_fill(fill_create_tail(self, blk)));
	}

	/* init fetch done, issue ungapped extension here if filter is needed */
	/* fill blocks then create a tail cap */
	return(_fill(fill_create_tail(self,
		fill_section_seq_bounded(self, blk)
	)));
}

/**
 * @fn gaba_dp_fill
 *
 * @brief fill API
 */
struct gaba_fill_s *_export(gaba_dp_fill)(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *fill,
	struct gaba_section_s const *a,
	struct gaba_section_s const *b,
	uint32_t pridx)
{
	self = _restore_dp_context(self);
	_init_bar(head); _init_bar(mid); _init_bar(tail);

	/* load current sections, then transpose sections to extract {id, len, base} pairs */
	v2i64_t asec = _loadu_v2i64(a), bsec = _loadu_v2i64(b);	/* tuple of (64bit ptr, 32-bit id, 32-bit len) */
	v2i64_t id_len = _lo_v2i64(asec, bsec), bptr = _hi_v2i64(asec, bsec);
	v2i32_t id = _cvt_v2i64_v2i32(id_len), len = _cvth_v2i64_v2i32(id_len);
	_print_v2i32(id); _print_v2i32(len); _print_v2i64(bptr);

	/* load sections */
	fill_load_section(self,
		_tail(fill), id, len, bptr,
		pridx == 0 ? _tail(fill)->pridx : pridx /* UINT32_MAX */
	);

	/* load sequences and extract the last block pointer */
	_print_v2i32(_load_v2i32(&_tail(fill)->aridx));
	_print_v2i32(_load_v2i32(&_tail(fill)->aadv));
	if(_stack_size(&self->stack) < MEM_INIT_VACANCY) {
		gaba_dp_add_stack(self, _mem_blocks(max_blocks_idx(self)));
	}
	struct gaba_block_s *blk = fill_load_vectors(self, _tail(fill));

	/* check if still in the init (head) state */
	if((int64_t)_tail(fill)->f.bpos < INIT_FETCH_BPOS) {
		if(fill_init_fetch(self, blk, _load_v2i64(&_tail(fill)->f.apos)) < INIT_FETCH_BPOS) {
			return(_fill(fill_create_tail(self, blk)));
		}
		/* init fetch done, issue ungapped extension here if filter is needed */
	}

	/* fill blocks then create a tail cap */
	return(_fill(fill_create_tail(self,
		fill_section_seq_bounded(self, blk)
	)));
}


/* merge bands */
/**
 * @macro _wb, _cb, _mb, _pb
 * @brief buffer accessor macros: wb[i] -> _wb(self, i)
 */
#define _wb(_t, _i)				( &(_t)->w.m.buf[_i] )		/* char buffer is allocated at the head */
#define _cb(_t, _i)				( &(_t)->w.m.buf[1 * MERGE_BUFFER_LENGTH + sizeof(int16_t) * (_i)] )
#define _mb(_t, _i)				( &(_t)->w.m.buf[3 * MERGE_BUFFER_LENGTH + sizeof(int16_t) * (_i)] )
#define _pb(_t, _i)				( &(_t)->w.m.buf[5 * MERGE_BUFFER_LENGTH + sizeof(int16_t) * (_i)] )
#define _eb(_t, _i)				( &(_t)->w.m.buf[7 * MERGE_BUFFER_LENGTH + sizeof(int16_t) * (_i)] )
#define _fb(_t, _i)				( &(_t)->w.m.buf[9 * MERGE_BUFFER_LENGTH + sizeof(int16_t) * (_i)] )

/**
 * @fn merge_calc_qspan
 * @brief detect the toprightmost and bottomleftmost vectors.
 * read qofs array (q-coordinate displacements), extract max and min
 * (and corresponding indices), and store them to working buffers.
 */
static _force_inline
int64_t merge_calc_qspan(
	struct gaba_dp_context_s *self,
	uint8_t const *qofs,
	uint32_t cnt)
{
	/* load offset vector */
	v16i8_t qv = _loadu_v16i8(qofs);

	/* suppress invalid elements after the last element */
	static uint8_t const mask_base[16] __attribute__(( aligned(16) )) = {
		-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16
	};
	v16i8_t mask = _add_v16i8(_load_v16i8(mask_base), _set_v16i8(cnt));
	qv = _sel_v16i8(qv, _set_v16i8(-127), mask);

	/* calc max, min */
	int8_t qub = _hmax_v16i8(qv), qlb = -_hmax_v16i8(_sub_v16i8(_zero_v16i8(), qv));
	if(qub - qlb > MERGE_BUFFER_LENGTH) { return(1); }
	self->w.m.qw = qub - qlb;

	/* calc max and min indices */
	self->w.m.uidx = tzcnt(((v16i8_masku_t){
		.mask = _mask_v16i8(_eq_v16i8(qv, _set_v16i8(qub)))
	}).all);
	self->w.m.lidx = tzcnt(((v16i8_masku_t){
		.mask = _mask_v16i8(_eq_v16i8(qv, _set_v16i8(qlb)))
	}).all);

	/* save the offsetted vector */
	_store_v16i8(&self->w.m.qofs, _sub_v16i8(qv, _set_v16i8(qlb)));
	return(0);
}

/**
 * @fn merge_init_work
 * @brief clear working buffers
 */
static _force_inline
void merge_init_work(
	struct gaba_dp_context_s *self)
{
	/* clear char vector */
	_memset_blk_a(_wb(self, 0), 0, sizeof(uint8_t) * MERGE_BUFFER_LENGTH);		/* 4x vmovdqa */

	/* fill INT16_MIN */
	v16i16_t v = _set_v16i16(INT16_MIN);
	#if MODEL == LINEAR
		uint64_t init_size = 3 * MERGE_BUFFER_LENGTH;
	#elif MODEL == AFFINE
		uint64_t init_size = 5 * MERGE_BUFFER_LENGTH;
	#else
		uint64_t init_size = 5 * MERGE_BUFFER_LENGTH;
	#endif

	for(uint64_t i = 0; i < init_size; i += 64) {
		_store_v16i16(_cb(self, i     ), v);		/* 4x unrolled vmovdqa */
		_store_v16i16(_cb(self, i + 16), v);
		_store_v16i16(_cb(self, i + 32), v);
		_store_v16i16(_cb(self, i + 48), v);
	}
	return;
}

/**
 * @fn merge_paste_score_vectors
 * @brief read score vectors from the tail and paste them at position q of the working vectors
 */
static _force_inline
void merge_paste_score_vectors(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail,
	int64_t offset,
	uint64_t q)
{
	/* char (just copy) */
	_storeu_n(_wb(self, q), _loadu_n(tail->ch.w));

	/* front vector (read-modify-write) */
	wvec_t md = _loadu_w(tail->md.delta);
	wvec_t cv = _add_w(md, _set_w(_offset(tail) - offset));
	_storeu_w(_cb(self, q), _max_w(_loadu_w(_cb(self, q)), cv));

	/* max */
	nvec_t xv = _loadu_n(tail->xd.drop);
	wvec_t mv = _add_w(cv, _cvt_n_w(xv));
	_storeu_w(_mb(self, q), _max_w(_loadu_w(_mb(self, q)), mv));

	/* previous */
	struct gaba_block_s const *blk = _last_block(tail);
	#if MODEL == LINEAR
		/* dh is not negated for linear: cv - dh; offset is ignored */
		nvec_t dh = _loadu_n(blk->diff.dh);
		wvec_t pv = _sub_w(cv, _cvt_n_w(dh));
		_storeu_w(_pb(self, q), _max_w(
			_loadu_w(_pb(self, q)),
			pv
		));
	#else	/* AFFINE and COMBINED */
		/* negated for affine: cv + -dh */
		nvec_t dh = _loadu_n(blk->diff.dh);
		wvec_t pv = _add_w(cv, _cvt_n_w(dh));
		_storeu_w(_pb(self, q), _max_w(
			_loadu_w(_pb(self, q)),
			pv
		));

		/* \delta E */
		nvec_t de = _loadu_n(blk->diff.de);
		wvec_t ev = _sub_w(cv, _cvt_n_w(de));
		_storeu_w(_eb(self, q), _max_w(
			_loadu_w(_eb(self, q)),
			ev
		));

		/* \delta F */
		nvec_t df = _loadu_n(blk->diff.df);
		wvec_t fv = _sub_w(cv, _cvt_n_w(df));
		_storeu_w(_fb(self, q), _max_w(
			_loadu_w(_fb(self, q)),
			fv
		));
	#endif
	return;
}

/**
 * @fn merge_paste_mask_vectors
 * @brief the same as merge_paste_score_vectors for breakpoint masks
 */
static _force_inline
void merge_paste_mask_vectors(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail,
	uint64_t q)
{
	uint64_t abrk = tail->abrk;
	uint64_t aidx = (self->w.m.qw - q) / 64, arem = (self->w.m.qw - q) % 64;
	self->w.m.abrk[aidx    ] = self->w.m.abrk[aidx    ] | abrk<<arem;
	self->w.m.abrk[aidx + 1] = self->w.m.abrk[aidx + 1] | (abrk>>(63 - arem))>>1;

	uint64_t bbrk = tail->bbrk;
	uint64_t bidx = q / 64, brem = q % 64;
	self->w.m.bbrk[bidx    ] = self->w.m.bbrk[bidx    ] | bbrk<<brem;
	self->w.m.bbrk[bidx + 1] = self->w.m.bbrk[bidx + 1] | (bbrk>>(63 - brem))>>1;
	return;
}

/**
 * @fn merge_paste_vectors
 * @brief iterate over tail objects to paste vectors to working buffers
 */
static _force_inline
void merge_paste_vectors(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *const *fill,
	uint32_t cnt)
{
	int64_t offset = _offset(_tail(fill[0]));

	/* calc max: read-modify-write loop */
	for(uint64_t i = 0; i < cnt; i++) {
		struct gaba_joint_tail_s const *tail = _tail(fill[i]);
		uint64_t q = self->w.m.qofs[i];

		/* score vectors */
		merge_paste_score_vectors(self, tail, offset, q);

		/* breakpoint masks */
		merge_paste_mask_vectors(self, tail, q);
	}
	return;
}

/**
 * @fn merge_detect_maxpos
 * @brief find maximum-scoring cell and the corresponding index in the working score vector (the pasted forefront vector)
 */
static _force_inline
uint64_t merge_detect_maxpos(
	struct gaba_dp_context_s *self)
{
	/* extract max value from the cv vector */
	v32i16_t cmaxv = _load_v32i16(_cb(self, 0));
	for(uint64_t i = 32; i < MERGE_BUFFER_LENGTH; i += 32) {
		cmaxv = _max_v32i16(cmaxv, _load_v32i16(_cb(self, i)));
	}
	cmaxv = _set_v32i16(_hmax_v32i16(cmaxv));

	/* calc index */
	uint64_t qmax = MERGE_BUFFER_LENGTH;
	for(uint64_t i = 0; i < MERGE_BUFFER_LENGTH; i += 32) {
		uint64_t q = tzcnt(((v32i8_masku_t){
			.mask = _mask_v32i16(_eq_v32i16(cmaxv, _load_v32i16(_cb(self, i))))
		}).all);
		qmax = MIN2(qmax, i + (q >= 64 ? MERGE_BUFFER_LENGTH : q));
	}

	/* determine span to be sliced */
	uint64_t const ofs = _W / 2, ub = MERGE_BUFFER_LENGTH - _W;
	qmax = MIN2(ub, qmax - ofs);					/* shift and clip at ub */
	qmax = (qmax < ub) ? qmax : 0;					/* clip negative val */
	return(qmax);
}

/**
 * @fn merge_slice_vectors
 */
static _force_inline
int32_t merge_slice_vectors(
	struct gaba_dp_context_s *self,
	struct gaba_merge_s *mg,
	uint64_t q)
{
	struct gaba_joint_tail_s *mt = _tail(mg + 1);

	/* copy char vector */
	_storeu_n(mt->ch.w, _loadu_n(_wb(self, q)));

	/* max vector */
	wvec_t mv = _loadu_w(_mb(self, q));
	int32_t mdrop = _hmax_w(mv);					/* extract max */

	wvec_t cv = _loadu_w(_cb(self, q));
	_storeu_n(mt->xd.drop, _cvt_w_n(_sub_w(mv, cv)));
	_storeu_w(mt->md.delta, cv);

	/* calc diff vectors */
	#if MODEL == LINEAR
		/* dh is not negated for linear-gap impl */
		wvec_t dh = _sub_w(cv, _loadu_w(_pb(self, q)));
		_storeu_n(mg->diff.dh, _cvt_w_n(dh));

		wvec_t dv = _sub_w(cv, _loadu_w(_pb(self, q - 1)));	/* q - 1 to make the cells in each lane vertically aligned */
		_storeu_n(mg->diff.dv, _cvt_w_n(dv));
	#elif MODEL == AFFINE
		/* dh is negated for affine */
		wvec_t dh = _sub_w(_loadu_w(_pb(self, q)), cv);
		_storeu_n(mg->diff.dh, _cvt_w_n(dh));

		wvec_t dv = _sub_w(cv, _loadu_w(_pb(self, q - 1)));	/* q - 1 to make the cells in each lane vertically aligned */
		_storeu_n(mg->diff.dv, _cvt_w_n(dv));

		wvec_t de = _sub_w(cv, _loadu_w(_eb(self, q)));
		_storeu_n(mg->diff.de, _cvt_w_n(de));

		wvec_t df = _sub_w(cv, _loadu_w(_fb(self, q)));
		_storeu_n(mg->diff.df, _cvt_w_n(df));
	#endif
	return(mdrop);
}

/**
 * @fn merge_restore_section
 * @brief follow the section links (tail object links) to restore (section, ridx) pair for the new vector
 */
static _force_inline
struct gaba_joint_tail_s const *merge_restore_section(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s *mt,
	struct gaba_joint_tail_s const *tail,
	uint32_t q,
	uint64_t i)
{
	#define _r(_x, _idx)		( (&(_x))[(_idx)] )

	/* slice breakpoint mask vector */
	uint64_t cidx = (self->w.m.qw - q) / 64, bidx = (self->w.m.qw - q) % 64;
	uint64_t brk = (
		   _r(self->w.m.abrk, i)[cidx    ]>>bidx
		| (_r(self->w.m.abrk, i)[cidx + 1]<<(63 - bidx))<<1
	);

	/* calc ridx */
	int32_t rem = q;
	uint32_t scnt = _r(tail->f.ascnt, i);
	while(rem > (int32_t)_r(tail->aadv, i)) {
		rem -= _r(tail->aadv, i);
		scnt -= _r(tail->aridx, i) == 0;
		tail = tail->tail;
	}

	/* copy section info */
	_r(mt->atptr, i) = _r(tail->atptr, i);
	_r(mt->f.aid, i) = _r(tail->f.aid, i);
	_r(mt->abrk, i) = brk;

	/* save ridx */
	_r(mt->aridx, i) = _r(tail->aridx, i) + rem;
	_r(mt->aadv, i) = 0;							/* always zero */
	_r(mt->f.ascnt, i) = scnt;

	#undef _r
	return(tail);
}

/**
 * @fn merge_create_tail
 */
static _force_inline
struct gaba_joint_tail_s *merge_create_tail(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *const *fill,
	uint32_t cnt)
{
	/* allocate merge-tail object */
	uint64_t arr_size = _roundup(cnt * sizeof(int64_t), MEM_ALIGN_SIZE);
	struct gaba_merge_s *mg = (ptrdiff_t)arr_size + gaba_dp_malloc(self,
		arr_size + sizeof(struct gaba_merge_s) + sizeof(struct gaba_joint_tail_s)
	);
	struct gaba_joint_tail_s *mt = _tail(mg + 1);

	/* copy tail pointer array */
	uint32_t pridx = UINT32_MAX;
	int64_t max = -1;
	// int64_t ppos = -1;
	for(uint64_t i = 0; i < cnt; i++) {
		struct gaba_joint_tail_s *const tail = _tail(fill[i]);
		mg->blk[-i] = _last_block(tail);
		pridx = MIN2(pridx, tail->pridx);
		max = MAX2(max, tail->f.max);
		// ppos = MAX2(ppos, tail->f.ppos);
	}

	/* save scores */
	mt->tail = NULL;					/* always NULL */
	mt->pridx = pridx;
	mt->f.max = max;
	// mt->f.ppos = ppos;

	/* determine center cell */
	uint64_t q = merge_detect_maxpos(self);

	/* convert q vector to relative offset from qmax then store to the new tail */
	_store_v16i8(&mg->acc,
		_swap_v16i8(_sub_v16i8(
			_load_v16i8(self->w.m.qofs),
			_set_v16i8(q)
		))
	);
	mg->acc = _cb(self, q + _W - 1) - _cb(self, q);
	mg->xstat = MERGE;					/* always MERGE; the other flags are not propagated (must be handled before the merge function is called) */

	/* copy vectors from the working buffer */
	mt->mdrop = merge_slice_vectors(self, mg, q);

	/* restore section boundary info */
	merge_restore_section(self, mt, _tail(fill[self->w.m.lidx]),                q, 0);
	merge_restore_section(self, mt, _tail(fill[self->w.m.uidx]), self->w.m.qw - q, 1);
	return(mt);
}

/**
 * @fn gaba_dp_merge
 *
 * @brief merging functionality (merge multiple banded matrices aligning on the same anti-diagonal)
 */
struct gaba_fill_s *_export(gaba_dp_merge)(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *const *fill,
	uint8_t const *qofs,
	uint32_t cnt)
{
	self = _restore_dp_context(self);
	_init_bar(head); _init_bar(mid); _init_bar(tail);

	/* clear working buffer */
	if(merge_calc_qspan(self, qofs, cnt) != 0) {
		return(NULL);					/* fill contains unmergable tail object */
	}
	merge_init_work(self);

	/* paste vectors */
	merge_paste_vectors(self, fill, cnt);

	/* slice vectors and copy them to a new tail */
	return(_fill(merge_create_tail(self, fill, cnt)));
}


/* max score search functions */
/**
 * @fn leaf_load_max_mask
 */
static _force_inline
uint64_t leaf_load_max_mask(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail)
{
	debug("pos(%lu, %lu), scnt(%u, %u), offset(%ld)",
		tail->f.apos, tail->f.bpos,
		tail->f.ascnt, tail->f.bscnt,
		tail->f.max - tail->mdrop);

	/* load max vector, create mask */
	nvec_t drop = _loadu_n(tail->xd.drop);
	wvec_t md = _loadu_w(tail->md.delta);
	uint64_t max_mask = _mask_u64(_mask_w(_eq_w(
		_set_w(tail->mdrop),
		_add_w(md, _cvt_n_w(drop))
	)));
	debug("max_mask(%lx)", max_mask);
	_print_w(_set_w(tail->mdrop));
	_print_w(_add_w(md, _cvt_n_w(drop)));
	return(max_mask);
}

/**
 * @fn leaf_search_pos
 */
static _force_inline
void leaf_search_pos(
	struct gaba_dp_context_s *self,
	nvec_masku_t const *mask_arr,
	nvec_masku_t const *m,
	uint64_t max_mask)
{
	/* search max cell, probe from the tail to the head */
	while(m > mask_arr && (max_mask & ~(--m)->all) != 0) {
		debug("max_mask(%lx), m(%lx)", max_mask, (uint64_t)m->all);
		max_mask &= ~m->all;
	}
	debug("max_mask(%lx), m(%lx)", max_mask, (uint64_t)m->all);
	self->w.l.p = m - mask_arr;
	self->w.l.q = tzcnt(m->all & max_mask);
	debug("p(%u), q(%u)", self->w.l.p, self->w.l.q);
}

#ifndef DEBUG_ALL
#  undef DEBUG
#  undef _LOG_H_INCLUDED
#  include "log.h"
#endif

/**
 * @fn leaf_detect_pos
 * @brief refill block to obtain cell-wise update-masks
 */
static _force_inline
void leaf_detect_pos(
	struct gaba_dp_context_s *self,
	struct gaba_block_s const *blk,
	uint64_t max_mask)
{
	#define _fill_block_leaf(_m) { \
		_dir_fetch(dir); \
		if(_dir_is_right(dir)) { \
			_fill_right_update_ptr(); \
			_fill_right(); \
		} else { \
			_fill_down_update_ptr(); \
			_fill_down(); \
		} \
		_m++->mask = _mask_n(_gt_n(delta, max)); \
		max = _max_n(delta, max); \
		debug("mask(%lx)", (uint64_t)_mask_u64((_m - 1)->mask)); \
	}

	/* load contexts and overwrite max vector */
	nvec_masku_t mask_arr[BLK], *m = mask_arr;		/* cell-wise update-mask array */
	/* vectors on registers */ {
		_fill_load_context(blk);
		nvec_t max = delta;
		for(int64_t i = 0; i < blk->acnt + blk->bcnt; i++) {
			_fill_block_leaf(m);
		}
	}

	/* search */
	leaf_search_pos(self, mask_arr, m, max_mask);
	return;
}

#ifdef REDEFINE_DEBUG
#  undef _LOG_H_INCLUDED
#  define DEBUG
#  include "log.h"
#endif

/**
 * @fn leaf_search
 * @brief returns resulting path length
 */
static _force_inline
uint64_t leaf_search(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail)
{
	/* load mask and block pointers */
	uint64_t max_mask = leaf_load_max_mask(self, tail);
	struct gaba_block_s const *b = _last_block(tail) + 1;
	debug("max_mask(%lx)", max_mask);

	_test_bar(head); _test_bar(mid); _test_bar(tail);

	/*
	 * iteratively clear lanes with longer paths;
	 * max_mask will be zero if the block contains the maximum scoring cell with the shortest path
	 */
	v2i32_t ridx = _load_v2i32(&tail->aridx);
	_print_v2i32(ridx);
	// if((b[-1].xstat & ROOT_HEAD) == ROOT_HEAD) { debug("reached root, xstat(%x)", b[-1].xstat); return(0); }	/* actually unnecessary but placed as a sentinel */
	while(1) {
		if(((--b)->xstat & ROOT) == ROOT) { debug("reached root, xstat(%x)", b->xstat); return(0); }	/* actually unnecessary but placed as a sentinel */
		while(_unlikely(b->xstat & HEAD)) { b = _phantom(b)->blk; }	/* sometimes head chains more than one */

		/* first adjust ridx to the head of this block then test mask was updated in this block */
		v2i8_t cnt = _load_v2i8(&b->acnt);
		ridx = _add_v2i32(ridx, _cvt_v2i8_v2i32(cnt));
		_print_v2i32(_cvt_v2i8_v2i32(cnt)); _print_v2i32(ridx);
		debug("max_mask(%lx), update_mask(%lx)", max_mask, (uint64_t)b->max_mask);
		if((max_mask & ~b->max_mask) == 0) { break; }
		max_mask &= ~b->max_mask;
	}

	/* calc (p, q) coordinates from block */
	fill_restore_fetch(self, tail, b, ridx);		/* fetch from existing blocks for p-coordinate search */
	leaf_detect_pos(self, b, max_mask);				/* calc local p,q-coordinates */
	self->w.l.blk = b;								/* max detection finished and reader_work has released, save block pointer to writer_work */

	/* restore reverse indices */
	int64_t fcnt = self->w.l.p + 1;					/* #filled vectors */
	uint32_t dir_mask = _dir_load(self->w.l.blk, fcnt).mask;/* 1 for b-side extension, 0 for a-side extension */
	_print_v2i32(ridx);
	ridx = _sub_v2i32(ridx,
		_seta_v2i32(
			(0    + popcnt(dir_mask)) - (_W - self->w.l.q),	/* bcnt - bqdiff where popcnt(dir_mask) is #b-side extensions */
			(fcnt - popcnt(dir_mask)) - (1 + self->w.l.q)	/* acnt - aqdiff */
		)
	);
	_print_v2i32(ridx);

	/* increment by one to convert char index to grid index */
	v2i32_t gidx = _sub_v2i32(_set_v2i32(1), ridx);
	gidx = _add_v2i32(gidx, _load_v2i32(&tail->aridx));		/* compensate tail lengths */

	_store_v2i32(&self->w.l.agidx, gidx);			/* current gidx */
	_store_v2i32(&self->w.l.asgidx, gidx);			/* tail gidx (pos) */
	_print_v2i32(gidx);

	/* calc plen */
	v2i32_t eridx = _load_v2i32(&tail->aridx);
	v2i32_t rem = _sub_v2i32(ridx, eridx);
	uint64_t plen = tail->f.apos + tail->f.bpos - (INIT_FETCH_APOS + INIT_FETCH_BPOS) + _W - _hi32(rem) - _lo32(rem);
	_print_v2i32(eridx); _print_v2i32(rem);
	debug("path length: plen(%lu, %lu), p(%u)", plen, plen % 32, self->w.l.p);
	return(plen);
}

/**
 * @fn gaba_dp_search_max
 */
struct gaba_pos_pair_s *_export(gaba_dp_search_max)(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *fill)
{
	self = _restore_dp_context(self);
	_init_bar(head); _init_bar(mid); _init_bar(tail);
	struct gaba_joint_tail_s const *tail = _tail(fill);

	struct gaba_pos_pair_s *pos = gaba_dp_malloc(self, sizeof(struct gaba_pos_pair_s));
	pos->plen = leaf_search(self, tail);		/* may return zero */

	v2i32_t const v11 = _seta_v2i32(1, 1);
	v2i32_t gidx = _load_v2i32(&self->w.l.agidx), acc = _zero_v2i32();
	v2i32_t id = _load_v2i32(&tail->f.aid);

	debug("gidx(%d, %d), acc(%d, %d)", _hi32(gidx), _lo32(gidx), _hi32(acc), _lo32(acc));
	while(tail->tail != NULL) {
		v2i32_t update = _gt_v2i32(v11, gidx);
		if(_test_v2i32(update, v11)) { break; }

		/* accumulate advanced lengths */
		v2i32_t adv = _load_v2i32(&tail->aadv), nid = _load_v2i32(&tail->f.aid);
		acc = _add_v2i32(acc, adv);
		_print_v2i32(acc); _print_v2i32(adv); _print_v2i32(nid);

		/* fetch previous tail and calculate the last update mask */
		tail = tail->tail;
		v2i32_t ridx = _load_v2i32(&tail->aridx);
		v2i32_t mask = _and_v2i32(update, _eq_v2i32(ridx, _zero_v2i32()));

		/* add advanced lengths */
		gidx = _add_v2i32(gidx, _and_v2i32(mask, acc));
		id = _sel_v2i32(mask, nid, id);			/* reload */
		acc = _andn_v2i32(mask, acc);			/* flush if summed up to gidx */
		_print_v2i32(gidx); _print_v2i32(id);
	}

	_store_v2i32(&pos->aid, id);
	_store_v2i32(&pos->apos, gidx);
	_print_v2i32(gidx);
	return(pos);
}


/* path trace functions */
/**
 * @fn trace_reload_section
 * @brief reload section and adjust gidx, pass 0 for section a or 1 for section b
 */
static _force_inline
void trace_reload_section(
	struct gaba_dp_context_s *self,
	uint64_t i)
{
	#define _r(_x, _idx)		( (&(_x))[(_idx)] )
	debug("load section %s, idx(%d), adv(%d)",
		i == 0 ? "a" : "b", _r(self->w.l.agidx, i), _r(_r(self->w.l.atail, i)->aadv, i));

	/* load tail pointer (must be inited with leaf tail) */
	struct gaba_joint_tail_s const *tail = _r(self->w.l.atail, i), *prev_tail = tail;
	int32_t gidx = _r(self->w.l.agidx, i);
	debug("base gidx(%d)", gidx);

	while(gidx <= 0) {
		do {
			gidx += tail->istat ? 0 : _r(tail->aadv, i);
			debug("add istat(%u), ridx(%d), adv(%d), gidx(%d), stat(%x)", tail->istat, _r(tail->aridx, i), _r(tail->aadv, i), gidx, tail->f.status);
			prev_tail = tail; tail = tail->tail;
		} while(_r(tail->aridx, i) != 0);
	}

	/* reload finished, store section info */
	_r(self->w.l.atail, i) = tail;				/* FIXME: is this correct?? */
	_r(self->w.l.aid, i) = _r(prev_tail->f.aid, i);

	/* save indices, offset for bridged region adjustment */
	_r(self->w.l.aofs, i) = prev_tail->istat ? _r(prev_tail->aadv, i) : 0;
	_r(self->w.l.agidx, i) = gidx;
	_r(self->w.l.asgidx, i) = gidx;
	debug("prev_tail->istat(%u), ofs(%u)", prev_tail->istat, _r(self->w.l.aofs, i));
	return;

	#undef _r
}

/**
 * @fn trace_push_segment
 */
static _force_inline
void trace_push_segment(
	struct gaba_dp_context_s *self)
{
	/* windback pointer */
	self->w.l.a.slen++;
	self->w.l.a.seg--;

	/* calc ppos */
	uint64_t ppos = (self->w.l.path - self->w.l.aln->path) * 32 + self->w.l.ofs;
	debug("ppos(%lu), path(%p, %p), ofs(%u), seg(%p)", ppos, self->w.l.path, self->w.l.aln->path, self->w.l.ofs, self->w.l.a.seg);

	/* load section info */
	v2i32_t ofs = _load_v2i32(&self->w.l.aofs);
	v2i32_t gidx = _load_v2i32(&self->w.l.agidx);
	v2i32_t sgidx = _load_v2i32(&self->w.l.asgidx);
	v2i32_t id = _load_v2i32(&self->w.l.aid);
	_print_v2i32(gidx); _print_v2i32(sgidx); _print_v2i32(id); _print_v2i32(ofs);

	/* store section info */
	// v2i32_t mask = _eq_v2i32(gidx, _zero_v2i32());/* add bridged length if the traceback pointer has reached the head */
	// v2i32_t pos = _add_v2i32(_and_v2i32(mask, ofs), gidx), len = _sub_v2i32(sgidx, gidx);
	v2i32_t pos = _add_v2i32(ofs, gidx), len = _sub_v2i32(sgidx, gidx);		/* add bridged length, is this correct? */
	_store_v2i32(&self->w.l.a.seg->aid, id);
	_store_v2i32(&self->w.l.a.seg->apos, pos);
	_store_v2i32(&self->w.l.a.seg->alen, len);
	self->w.l.a.seg->ppos = ppos;

	/* update rsgidx */
	_store_v2i32(&self->w.l.asgidx, gidx);
	return;
}

/**
 * @enum ts_*
 */
#define TS_H							( 0x01 )
#define TS_V							( 0x02 )
#define TS_S							( 0x04 )
enum {
	ts_d  = TS_H | TS_V,
	ts_v0 = TS_V,
	ts_v1 = TS_V | TS_S,
	ts_h0 = TS_H,
	ts_h1 = TS_H | TS_S
};

/**
 * @macro _trace_inc_*
 * @brief increment gap counters
 */

#define _trace_inc_gi_h()				{ gi = _sub_v2i32(gi, v01); }
#define _trace_inc_gi_v()				{ gi = _sub_v2i32(gi, v10); }
#define _trace_inc_ge_h()				{ ge = _sub_v2i32(ge, v01); }
#define _trace_inc_ge_v()				{ ge = _sub_v2i32(ge, v10); }
#define _trace_inc_gf_h()				{ gf = _sub_v2i32(gf, v01); }
#define _trace_inc_gf_v()				{ gf = _sub_v2i32(gf, v10); }
/*
#define _trace_inc_gi()					;
#define _trace_inc_ge()					;
#define _trace_inc_gf()					;
*/

/**
 * @macro _trace_*_*_test_index
 * @brief test if path reached a section boundary
 */
#define _trace_bulk_v_test_index()		( 0 )
#define _trace_bulk_d_test_index()		( 0 )
#define _trace_bulk_h_test_index()		( 0 )
#define _trace_tail_v_test_index()		( _test_v2i32(gidx, v10) /* bgidx == 0 */ )
#define _trace_tail_d_test_index()		( !_test_v2i32(_eq_v2i32(gidx, v00), v11) /* agidx == 0 || bgidx == 0 */ )
#define _trace_tail_h_test_index()		( _test_v2i32(gidx, v01) /* agidx == 0 */ )

/**
 * @macro _trace_*_*_update_index
 * @brief update indices by one
 */
#define _trace_bulk_v_update_index()	;
#define _trace_bulk_h_update_index()	;
#define _trace_tail_v_update_index()	{ gidx = _add_v2i32(gidx, v10); _print_v2i32(gidx); }
#define _trace_tail_h_update_index()	{ gidx = _add_v2i32(gidx, v01); _print_v2i32(gidx); }

/**
 * @macro _trace_test_*
 * @brief test mask
 */
#if MODEL == LINEAR
#define _trace_test_diag_h()			( ((mask->h.all>>q) & 0x01) == 0 )
#define _trace_test_diag_v()			( ((mask->v.all>>q) & 0x01) == 0 )
#define _trace_test_gap_h()				( ((mask->h.all>>q) & 0x01) != 0 )
#define _trace_test_gap_v()				( ((mask->v.all>>q) & 0x01) != 0 )
#define _trace_test_fgap_h()			( 0 )
#define _trace_test_fgap_v()			( 0 )
#elif MODEL == AFFINE
#define _trace_test_diag_h()			( ((mask->h.all>>q) & 0x01) == 0 )
#define _trace_test_diag_v()			( ((mask->v.all>>q) & 0x01) == 0 )
#define _trace_test_gap_h()				( ((mask->e.all>>q) & 0x01) == 0 )
#define _trace_test_gap_v()				( ((mask->f.all>>q) & 0x01) == 0 )
#define _trace_test_fgap_h()			( 0 )
#define _trace_test_fgap_v()			( 0 )
#else /* MODEL == COMBINED */
#define _trace_test_diag_h()			( ((mask->h.all>>q) & 0x01) == 0 )
#define _trace_test_diag_v()			( ((mask->v.all>>q) & 0x01) == 0 )
#define _trace_test_gap_h()				( (((~mask->h.all & mask->e.all)>>q) & 0x01) == 0 )
#define _trace_test_gap_v()				( (((~mask->v.all & mask->f.all)>>q) & 0x01) == 0 )
#define _trace_test_fgap_h()			( ((mask->e.all>>q) & 0x01) == 0 )
#define _trace_test_fgap_v()			( ((mask->f.all>>q) & 0x01) == 0 )
#endif

/**
 * @macro _trace_*_update_path_q
 */
#define _trace_h_update_path_q() { \
	path_array = path_array<<1; \
	mask--; \
	q += _dir_mask_is_down(dir_mask); \
	/*fprintf(stderr, "h: q(%u)\n", q);*/ \
	_dir_mask_windback(dir_mask); \
}
#define _trace_v_update_path_q() { \
	path_array = (path_array<<1) | 0x01; \
	mask--; \
	q += _dir_mask_is_down(dir_mask) - 1; \
	/*fprintf(stderr, "v: q(%u)\n", q);*/ \
	_dir_mask_windback(dir_mask); \
}

/**
 * @macro _trace_load_block_rem
 * @brief reload mask pointer and direction mask, and adjust path offset
 */
#define _trace_load_block_rem(_cnt) ({ \
	debug("ofs(%u), cnt(%lu), path(%p), ofs(%u)", \
		ofs, (uint64_t)_cnt, path - (ofs < (_cnt)), (uint32_t)(ofs - (_cnt)) & (BLK - 1)); \
	path -= ofs < (_cnt); \
	ofs = (ofs - (_cnt)) & (BLK - 1); \
	_dir_mask_load(blk, (_cnt)); \
})

/**
 * @macro _trace_reload_tail
 * @brief reload tail, issued at each band-segment boundaries
 */
#define _trace_reload_tail(t, _vec_idx) { \
	/* store path (offset will be adjusted afterwards) */ \
	_storeu_u64(path, path_array<<ofs); \
	/* reload block pointer */ \
	blk--; do { \
		debug("reload head block, blk(%p), prev_blk(%p), head(%x), cnt(%u, %u)", blk, _phantom(blk)->blk, _phantom(blk)->blk->xstat & HEAD, _phantom(blk)->blk->acnt, _phantom(blk)->blk->bcnt); \
		blk = _phantom(blk)->blk; \
	} while((blk->xstat & HEAD) != 0); \
	while(blk->xstat & MERGE) { \
		struct gaba_merge_s const *_mg = _merge(blk); \
		blk = _mg->blk[_mg->tidx[_vec_idx][q]]; \
		q += _mg->qofs[_mg->tidx[_vec_idx][q]];	/* adjust q, restore block pointer; FIXME: tail pointer must be saved somewhere */ \
	} \
	/* reload dir and mask pointer, adjust path offset */ \
	uint64_t _cnt = blk->acnt + blk->bcnt; \
	mask = &blk->mask[_cnt - 1]; dir_mask = _trace_load_block_rem(_cnt); \
	debug("reload tail, path_array(%lx), blk(%p), idx(%lu), mask(%x)", path_array, blk, _cnt - 1, dir_mask); \
}

/**
 * @macro _trace_reload_block
 * @brief reload block for bulk trace
 */
#define _trace_reload_block() { \
	/* store path (bulk, offset does not change here) */ \
	_storeu_u64(path, path_array<<ofs); path--; \
	/* reload mask and mask pointer; always test the boundary */ \
	mask = &(--blk)->mask[BLK - 1]; dir_mask = _dir_mask_load(blk, BLK); \
	if(_unlikely((_phantom(blk)->xstat & HEAD) != 0)) { \
		do { blk = _phantom(blk)->blk; debug("reload block, cnt(%u, %u)", blk->acnt, blk->bcnt); } while((_phantom(blk)->xstat & HEAD) != 0); \
		uint64_t _cnt = blk->acnt + blk->bcnt; path++; \
		mask = &blk->mask[_cnt - 1]; dir_mask = _trace_load_block_rem(_cnt); \
	} \
	debug("reload block, path_array(%lx), blk(%p), head(%x), mask(%x), ofs(%u), cnt(%u, %u)", path_array, blk, blk->xstat & HEAD, dir_mask, ofs, blk->acnt, blk->bcnt); \
}

/**
 * @macro _trace_test_bulk
 * @brief update gidx and test if section boundary can be found in the next block.
 * adjust gidx to the next base and return 0 if bulk loop can be continued,
 * otherwise gidx will not updated. must be called just after reload_ptr or reload_tail.
 */
#if 1
#define _trace_test_bulk() ({ \
	v2i8_t _cnt = _load_v2i8(&blk->acnt); \
	v2i32_t _gidx = _sub_v2i32(gidx, _cvt_v2i8_v2i32(_cnt)); \
	debug("test bulk, xstat(%x), save(%u), _cnt(%d, %d), _gidx(%d, %d), test(%d)", \
		blk->xstat, save, \
		_hi32(_cvt_v2i8_v2i32(_cnt)), _lo32(_cvt_v2i8_v2i32(_cnt)), _hi32(_gidx), _lo32(_gidx), \
		_test_v2i32(_gt_v2i32(bw, _gidx), v11) ? 1 : 0); \
	_test_v2i32(_gt_v2i32(bw, _gidx), v11) ? (gidx = _gidx, 1) : 0;	/* agidx >= BW && bgidx >= BW */ \
})
#else
#define _trace_test_bulk()	0
#endif

/**
 * @macro _trace_*_load
 * @brief set _state 0 when in diagonal loop, otherwise pass 1
 */
#define TRACE_HEAD_CNT			( _W / BLK + (_W == 16) )
#define _trace_bulk_load_n(t, _state, _jump_to) { \
	if(_unlikely(mask < blk->mask)) { \
		_trace_reload_block(); \
		if(_unlikely(!_trace_test_bulk())) {	/* adjust gidx */ \
			if(q >= _W) { goto _trace_term; }	/* out-of-bound check */ \
			gidx = _add_v2i32(gidx, _seta_v2i32(q - save, save - q)); \
			debug("jump to %s, adjust gidx, q(%d), prev_q(%d), adj(%d, %d), gidx(%u, %u)", #_jump_to, \
				q, self->w.l.q, q - save, save - q, _hi32(gidx), _lo32(gidx)); \
			save = TRACE_HEAD_CNT; \
			goto _jump_to;						/* jump to tail loop */ \
		} \
	} \
}
#define _trace_tail_load_n(t, _state, _jump_to) { \
	if(_unlikely(mask < blk->mask)) { \
		debug("test ph(%p), xstat(%x), ppos(%lu), path(%p, %p), ofs(%u)", \
			_last_phantom(blk), _last_phantom(blk)->xstat, (path - self->w.l.aln->path) * 32 + ofs, path, self->w.l.aln->path, ofs); \
		if(_last_phantom(blk)->xstat & HEAD) {	/* head (phantom) block is marked 0x4000 */ \
			_trace_reload_tail(t, _state);		/* fetch the previous tail */ \
		} else { \
			/* load dir and update mask pointer */ \
			_trace_reload_block();				/* not reached the head yet */ \
			if(--save >= TRACE_HEAD_CNT && _trace_test_bulk()) {	/* adjust gidx, NOTE: both must be evalueted every time */ \
				debug("save q(%d)", q); \
				save = q; \
				goto _jump_to;					/* dispatch again (bulk loop) */ \
			} \
		} \
	} \
}

#ifndef DEBUG_ALL
#  undef DEBUG
#  undef _LOG_H_INCLUDED
#  include "log.h"
#endif

/**
 * @fn trace_core
 */
static _force_inline
void trace_core(
	struct gaba_dp_context_s *self)
{
	#define _pop_vector(_c, _l, _state, _jump_to) { \
		debug("go %s (%s, %s), dir(%x), mask(%lx, %lx), h(%lx, %lx, %lx), v(%lx, %lx, %lx), p(%ld), q(%d), ptr(%p), path_array(%lx)", \
			#_l, #_c, #_jump_to, dir_mask, (uint64_t)mask->h.all, (uint64_t)mask->v.all, \
			(uint64_t)_trace_test_diag_h(), (uint64_t)_trace_test_gap_h(), (uint64_t)_trace_test_fgap_h(), (uint64_t)_trace_test_diag_v(), (uint64_t)_trace_test_gap_v(), (uint64_t)_trace_test_fgap_v(), \
			(int64_t)(mask - blk->mask), (int32_t)q, mask, path_array); \
		_trace_##_c##_##_l##_update_index(); \
		_trace_##_l##_update_path_q(); \
		_trace_##_c##_load_n(t, _state, _jump_to); \
	}
	#define _trace_gap_loop(t, _c, _n, _l) { \
		_trace_##_c##_##_l##_head: \
			if(_trace_test_fgap_##_l()) { \
				if(_unlikely(_trace_##_c##_##_l##_test_index())) { \
					self->w.l.state = ts_##_l##0; goto _trace_term; \
				} \
				_trace_inc_gf_##_l(); \
				_pop_vector(_c, _l, 1, _trace_##_n##_##_l##_retd); \
				goto _trace_##_c##_##_l##_retd; \
			} \
			_trace_inc_gi_##_l(); \
			do { \
				if(_unlikely(_trace_##_c##_##_l##_test_index())) { \
					self->w.l.state = ts_##_l##1; goto _trace_term; \
				} \
				_trace_inc_ge_##_l();	/* increment #gap bases on every iter */ \
				_pop_vector(_c, _l, 1, _trace_##_n##_##_l##_tail); \
			_trace_##_c##_##_l##_tail:; \
			} while(_trace_test_gap_##_l()); \
			goto _trace_##_c##_##_l##_retd; \
	}
	#define _trace_diag_loop(t, _c, _n) { \
		while(1) { \
		_trace_##_c##_h_retd: _trace_##_c##_d_head: \
			if(!_trace_test_diag_h()) { goto _trace_##_c##_h_head; } \
			if(_unlikely(_trace_##_c##_d_test_index())) { \
				self->w.l.state = ts_d; goto _trace_term; \
			} \
			_pop_vector(_c, h, 0, _trace_##_n##_d_mid); \
		_trace_##_c##_d_mid: \
			_pop_vector(_c, v, 0, _trace_##_n##_d_tail); \
		_trace_##_c##_v_retd: _trace_##_c##_d_tail: \
			if(!_trace_test_diag_v()) { goto _trace_##_c##_v_head; } \
		} \
	}

	/* debug */
	_test_bar(head); _test_bar(mid); _test_bar(tail);

	/* constants for agidx and bgidx end detection */
	v2i32_t const v00 = _seta_v2i32(0, 0);
	v2i32_t const v01 = _seta_v2i32(0, -1);
	v2i32_t const v10 = _seta_v2i32(-1, 0);
	v2i32_t const v11 = _seta_v2i32(-1, -1);
	v2i32_t const bw = _set_v2i32(_W);

	/* gap counts; #gap regions in lower and #gap bases in higher */
	register v2i32_t gi = _load_v2i32(&self->w.l.a.aicnt);
	register v2i32_t ge = _load_v2i32(&self->w.l.a.aecnt);
	register v2i32_t gf = _load_v2i32(&self->w.l.afcnt);

	/* load path array, adjust path offset to align the head of the current block */
	uint32_t ofs = self->w.l.ofs;	/* global p-coordinate */
	uint32_t *path = self->w.l.path;
	uint64_t path_array = _loadu_u64(path)>>ofs;

	/* load pointers and coordinates */
	struct gaba_block_s const *blk = self->w.l.blk;
	struct gaba_mask_pair_s const *mask = &blk->mask[self->w.l.p];
	uint32_t q = self->w.l.q, save = TRACE_HEAD_CNT;
	uint32_t dir_mask = _trace_load_block_rem(self->w.l.p + 1);

	/* load grid indices; assigned for each of N+1 boundaries of length-N sequence */
	register v2i32_t gidx = _load_v2i32(&self->w.l.agidx);
	_print_v2i32(gidx);

	/* dispatcher: jump to the last state */
	switch(self->w.l.state) {
		case ts_d:  goto _trace_tail_d_head;
		case ts_v0: goto _trace_tail_v_head;
		case ts_v1: goto _trace_tail_v_tail;
		case ts_h0: goto _trace_tail_h_head;
		case ts_h1: goto _trace_tail_h_tail;
		default: /* trap(); */ return;			/* broken state */
	}

	/* tail loops */ {
		_trace_gap_loop(self, tail, bulk, v);	/* v is the beginning state */
		_trace_diag_loop(self, tail, bulk);
		_trace_gap_loop(self, tail, bulk, h);
	}
	/* bulk loops */ {
		_trace_gap_loop(self, bulk, tail, v);
		_trace_diag_loop(self, bulk, tail);
		_trace_gap_loop(self, bulk, tail, h);
	}

_trace_term:;
	/* reached a boundary of sections, compensate ofs */
	uint64_t rem = mask - blk->mask + 1;
	path += (ofs + rem) >= BLK;
	ofs = (ofs + rem) & (BLK - 1);
	_storeu_u64(path, path_array<<ofs);
	debug("rem(%lu), path(%p), arr(%lx), ofs(%u)", rem, path, path_array<<ofs, ofs);

	/* save gap counts */
	_store_v2i32(&self->w.l.a.aicnt, gi);
	_store_v2i32(&self->w.l.a.aecnt, ge);
	_store_v2i32(&self->w.l.afcnt, gf);

	/* store pointers and coordinates */
	self->w.l.path = path;
	self->w.l.blk = blk;
	self->w.l.p = mask - blk->mask;	/* local p-coordinate */
	self->w.l.q = q;				/* q coordinate */
	self->w.l.ofs = ofs;			/* global p-coordinate */
	_store_v2i32(&self->w.l.agidx, gidx);
	_print_v2i32(gidx);
	return;
}

#ifdef REDEFINE_DEBUG
#  undef _LOG_H_INCLUDED
#  define DEBUG
#  include "log.h"
#endif

/**
 * @fn trace_init
 */
static _force_inline
void trace_init(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail,
	struct gaba_alloc_s const *alloc,
	uint64_t plen)
{
	/* store tail pointers for sequence segment tracking */
	self->w.l.atail = tail;
	self->w.l.btail = tail;

	/* malloc container */
	uint64_t sn = tail->f.ascnt + tail->f.bscnt + 2, pn = (plen + 31) / 32 + 2;
	uint64_t size = (
		  sizeof(struct gaba_alignment_s)				/* base */
		+ sizeof(uint32_t) * _roundup(pn, 8)			/* path array and its margin */
		+ sizeof(struct gaba_segment_s) * sn			/* segment array */
	);

	/* save aln pointer and memory management stuffs to working buffer */
	self->w.l.aln = alloc->lmalloc(alloc->opaque, size);
	self->w.l.a.opaque = alloc->opaque;					/* save opaque pointer */
	self->w.l.a.lfree = alloc->lfree;

	/* use gaba_alignment_s buffer instead in the traceback loop */
	self->w.l.a.score = tail->f.max;					/* just copy */
	_store_v2i32(&self->w.l.a.aicnt, _zero_v2i32());	/* clear counters */
	_store_v2i32(&self->w.l.a.aecnt, _zero_v2i32());
	_store_v2i32(&self->w.l.afcnt, _zero_v2i32());
	self->w.l.a.dcnt = 0;

	/* section and path */
	self->w.l.a.slen = 0;
	self->w.l.a.seg = sn + (struct gaba_segment_s *)(self->w.l.aln->path + _roundup(pn, 8)),
	self->w.l.a.plen = plen;							/* save path length */
	self->w.l.a.padding = 0x40000000;

	/* store block and coordinates */
	self->w.l.ofs = plen & (32 - 1);
	self->w.l.state = ts_d;								/* clear state, the traceback always starts with a match */
	self->w.l.path = self->w.l.aln->path + plen / 32;
	debug("sn(%lu), seg(%p), pn(%lu), path(%p)", sn, self->w.l.a.seg, pn, self->w.l.path);

	/* clear array (FIXME: is head cap needed? open uint64_t head_cap; before uint32_t path[] and set 0x1<<62) if so */
	self->w.l.path[0] = 0x01<<self->w.l.ofs;
	self->w.l.path[1] = 0;
	self->w.l.path[2] = 0;								/* make sure the memory is "initialized" so that not warned by valgrind */
	return;
}

/**
 * @fn trace_body
 * @brief plen = 0 generates an alignment object with no section (not NULL object).
 * may fail when path got lost out of the band and returns NULL object
 */
static _force_inline
struct gaba_alignment_s *trace_body(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail,
	struct gaba_alloc_s const *alloc,
	uint64_t plen)
{
	/* create alignment object */
	trace_init(self, tail, alloc, plen);

	/* blockwise traceback loop, until ppos reaches the root */
	while(self->w.l.path + self->w.l.ofs > self->w.l.aln->path) {	/* !(self->w.l.path == self->w.l.aln->path && self->w.l.ofs == 0) */
		_test_bar(head); _test_bar(mid); _test_bar(tail);

		/* update section info; check the next direction */
		debug("gidx(%d, %d)", self->w.l.bgidx, self->w.l.agidx);
		if((int32_t)self->w.l.agidx < (int32_t)((self->w.l.state & TS_H) != 0)) {
			trace_reload_section(self, 0);
		}
		if((int32_t)self->w.l.bgidx < (int32_t)((self->w.l.state & TS_V) != 0)) {
			trace_reload_section(self, 1);
		}

		/* fragment trace: q must be inside [0, BW) */
		trace_core(self);
		debug("p(%d), q(%d)", self->w.l.p, self->w.l.q);
		if(_unlikely(self->w.l.q >= _W)) {
			/* out of band: abort */
			self->w.l.a.lfree(self->w.l.a.opaque, (void *)((uint8_t *)self->w.l.aln));
			return(NULL);
		}

		/* push section info to segment array */
		trace_push_segment(self);
	}

	/* estimate alignment identity */
	v2i32_t gicnt = _load_v2i32(&self->w.l.a.aicnt);
	v2i32_t gecnt = _load_v2i32(&self->w.l.a.aecnt);
	v2i32_t gfcnt = _load_v2i32(&self->w.l.afcnt);
	v2i32_t gcnt = _add_v2i32(gecnt, gfcnt);

	v2i32_t gi = _set_v2i32((int32_t)self->gi), ge = _set_v2i32((int32_t)self->ge);
	v2i32_t gf = _cvt_v2i8_v2i32(_load_v2i8(&self->gfa));
	v2i32_t g = _add_v2i32(_add_v2i32(
		_mul_v2i32(gi, gicnt),
		_mul_v2i32(ge, gecnt)),
		_mul_v2i32(gf, gfcnt)
	);

	uint64_t dlen = (self->w.l.a.plen - _hi32(gcnt) - _lo32(gcnt))>>1;
	int64_t dsc = self->w.l.a.score + _hi32(g) + _lo32(g);

	/* copy (overwrite identity and gap counts) */
	_memcpy_blk_ua(self->w.l.aln, &self->w.l.a, sizeof(struct gaba_alignment_s));
	self->w.l.aln->identity = dlen == 0 ? 0.0 : (((double)dsc / (double)dlen) * self->imx - self->xmx);
	_store_v2i32(&self->w.l.aln->agcnt, gcnt);
	self->w.l.aln->dcnt = dlen;

	/*
	fprintf(stderr, "plen(%lu), g(%d, %d, %d, %d), gic(%u, %u), gec(%u, %u), gfc(%u, %u), gc(%u, %u), dlen(%ld), sc(%ld, %ld), mc(%f), identity(%f)\n",
		plen, self->gi, self->ge, self->gfa, self->gfb,
		self->w.l.a.bicnt, self->w.l.a.aicnt,
		self->w.l.a.becnt, self->w.l.a.aecnt,
		self->w.l.bfcnt, self->w.l.afcnt,
		_hi32(gcnt), _lo32(gcnt),
		dlen, self->w.l.aln->score, dsc, dsc * self->imx - self->xmx * dlen, self->w.l.aln->identity);
	*/
	return(self->w.l.aln);
}

/**
 * @fn gaba_dp_trace
 */
struct gaba_alignment_s *_export(gaba_dp_trace)(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *fill,
	struct gaba_alloc_s const *alloc)
{
	/* restore dp context pointer by adding offset */
	self = _restore_dp_context(self);
	_init_bar(head); _init_bar(mid); _init_bar(tail);

	/* restore default alloc if NULL */
	struct gaba_alloc_s const default_alloc = {
		.opaque = (void *)self,
		.lmalloc = (gaba_lmalloc_t)gaba_dp_malloc,
		.lfree = (gaba_lfree_t)gaba_dp_free
	};
	alloc = (alloc == NULL) ? &default_alloc : alloc;

	/* search and trace */
	return(trace_body(self, _tail(fill), alloc,
		(int64_t)fill->bpos < INIT_FETCH_BPOS ? 0 : leaf_search(self, _tail(fill))
	));
}


/**
 * @fn gaba_dp_res_free
 */
void _export(gaba_dp_res_free)(
	struct gaba_dp_context_s *self,
	struct gaba_alignment_s *aln)
{
	struct gaba_aln_intl_s *a = (struct gaba_aln_intl_s *)aln;
	debug("free mem, ptr(%p), lmm(%p)", a, a->opaque);
	a->lfree(a->opaque, (void *)a);
	return;
}


/* score accumulation auxiliary macros */
/**
 * @macro _hadd_v16i8
 */
#define _hadd_v16i8(_v) ({ \
	v16i8_t v = (_v); \
	v = _add_v16i8(v, _bsr_v16i8(v, 8)); \
	v = _add_v16i8(v, _bsr_v16i8(v, 4)); \
	v = _add_v16i8(v, _bsr_v16i8(v, 2)); \
	_ext_v16i8(v, 1) + _ext_v16i8(v, 0); \
})

#if MODEL == LINEAR || MODEL == AFFINE
#define _del(_c) { \
	v2i32_t cnt = _seta_v2i32((_c) > 0, _c); gbc = _add_v2i32(gbc, cnt); \
}
#define _ins(_c) { \
	v2i32_t cnt = _seta_v2i32((_c) > 0, _c); gac = _add_v2i32(gac, cnt); \
}
#else
#define _del(_c) { \
	v2i32_t cnt = _seta_v2i32((_c) > 0, _c); \
	gbc = _add_v2i32(gbc, cnt); fbc = (_c) > self->bflen ? fbc : _add_v2i32(fbc, cnt); \
}
#define _ins(_c) { \
	v2i32_t cnt = _seta_v2i32((_c) > 0, _c); \
	gac = _add_v2i32(gac, cnt); fac = (_c) > self->aflen ? fac : _add_v2i32(fac, cnt); \
}
#endif

#define _del_f(_c) { ap += (_c); gcnt += (_c); _del(_c); }		/* horizontal gap; leaq (%r1, %r2, 2), %r1 */
#define _del_r(_c) { ap -= (_c); gcnt += (_c); _del(_c); }
#define _ins_f(_c) { bp += (_c); gcnt -= (_c); _ins(_c); }		/* vertical gap; leaq 1(%r1, %r2, 1), %r1 */
#define _ins_r(_c) { bp -= (_c); gcnt -= (_c); _ins(_c); }
#define _match_core(_rv, _qv, _c) { \
	_print_v16i8(_rv); _print_v16i8(_qv); \
	v16i8_t sp = _swapn_v16i8(_shuf_v16i8(sb, _match_v16i8((_rv), (_qv))), (_c));				/* masked score profile */ \
	debug("hadd(%d), cnt(%u), xcnt(%u)", _hadd_v16i8(sp), _c, popcnt(((v16i8_masku_t){ .mask = _mask_v16i8(sp) }).all)); \
	sacc = _add_v16i8(sacc, sp); xc += popcnt(((v16i8_masku_t){ .mask = _mask_v16i8(sp) }).all);/* count mismatches */ \
}
#define _match_ff(_c) { \
	v16i8_t sacc = _zero_v16i8(); \
	dc += (_c); ap += (_c); bp += (_c); debug("cnt(%lu)", _c); \
	for(uint64_t i = (_c); i > 0; i -= _gaba_parse_min2(i, 16)) { \
		uint64_t l = _gaba_parse_min2(i, 16); \
		v16i8_t av = _fwap_v16i8(_loadu_v16i8(ap - i), l), bv = _fwbp_v16i8(_loadu_v16i8(bp - i), l); _match_core(av, bv, l); \
	} \
	score += _hadd_v16i8(sacc); debug("score(%ld), dcnt(%lu), xcnt(%lu)", score, dc, xc); \
}
#define _match_fr(_c) { \
	v16i8_t sacc = _zero_v16i8(); \
	dc += (_c); ap += (_c); bp -= (_c); debug("cnt(%lu)", _c); \
	for(uint64_t i = (_c); i > 0; i -= _gaba_parse_min2(i, 16)) { \
		uint64_t l = _gaba_parse_min2(i, 16); \
		v16i8_t av = _fwap_v16i8(_loadu_v16i8(ap - i), l), bv = _rvbp_v16i8(_loadu_v16i8(bp + i - l), l); _match_core(av, bv, l); \
	} \
	score += _hadd_v16i8(sacc); debug("score(%ld), dcnt(%lu), xcnt(%lu)", score, dc, xc); \
}
#define _match_rf(_c) { \
	v16i8_t sacc = _zero_v16i8(); \
	dc += (_c); ap -= (_c); bp += (_c); debug("cnt(%lu)", _c); \
	for(uint64_t i = (_c); i > 0; i -= _gaba_parse_min2(i, 16)) { \
		uint64_t l = _gaba_parse_min2(i, 16); \
		v16i8_t av = _rvap_v16i8(_loadu_v16i8(ap + i - l), l), bv = _fwbp_v16i8(_loadu_v16i8(bp - i), l); _match_core(av, bv, l); \
	} \
	score += _hadd_v16i8(sacc); debug("score(%ld), dcnt(%lu), xcnt(%lu)", score, dc, xc); \
}
#define _match_rr(_c) { \
	v16i8_t sacc = _zero_v16i8(); \
	dc += (_c); ap -= (_c); bp -= (_c); debug("cnt(%lu)", _c); \
	for(uint64_t i = (_c); i > 0; i -= _gaba_parse_min2(i, 16)) { \
		uint64_t l = _gaba_parse_min2(i, 16); \
		v16i8_t av = _rvap_v16i8(_loadu_v16i8(ap + i - l), l), bv = _rvbp_v16i8(_loadu_v16i8(bp + i - l), l); _match_core(av, bv, l); \
	} \
	score += _hadd_v16i8(sacc); debug("score(%ld), dcnt(%lu), xcnt(%lu)", score, dc, xc); \
}
#define _match_end(_c) { gcnt += (_c)<<8; if((int64_t)ridx > 0) { gcnt = 0; } }

/**
 * @fn gaba_dp_calc_score
 * @brief calculate score, match count, mismatch count, and gap counts for the section
 * NOTE: this function depends on the _parser_init_* and _parser_loop_* macros defined in gaba_parse.h
 */
struct gaba_score_s *_export(gaba_dp_calc_score)(
	struct gaba_dp_context_s *self,
	uint32_t const *path,
	gaba_path_section_t const *s,
	gaba_section_t const *a,
	gaba_section_t const *b)
{
	self = _restore_dp_context(self);
	_init_bar(head); _init_bar(mid); _init_bar(tail);

	v16i8_t sb = _sub_v16i8(_to_v16i8_n(_load_sb(self->scv)), _set_v16i8(self->ofs));
	v2i32_t gac = _zero_v2i32(), fac = _zero_v2i32(), gbc = _zero_v2i32(), fbc = _zero_v2i32();
	uint64_t xc = 0, dc = 0;
	int64_t score = 0;

	uint8_t const *ap = a->base < GABA_EOU ? &a->base[s->apos] : gaba_mirror(&a->base[s->apos], 0);
	uint8_t const *bp = b->base < GABA_EOU ? &b->base[s->bpos] : gaba_mirror(&b->base[s->bpos], 0);
	_parser_init_fw(path, s->ppos, gaba_plen(s));

	uint64_t harr = gaba_parse_u64(p, lim - ridx - 2);
	int32_t gcnt = 0xb3000000<<(harr & 0x07) & 0x80000000;	/* == 0 when the last element is M */
	switch(((a->base >= GABA_EOU)<<1) | (b->base >= GABA_EOU)) {
		case 0x00: _parser_loop_fw(_del_f, _ins_f, _match_ff, _match_end); break;
		case 0x01: _parser_loop_fw(_del_f, _ins_r, _match_fr, _match_end); break;
		case 0x02: _parser_loop_fw(_del_r, _ins_f, _match_rf, _match_end); break;
		case 0x03: _parser_loop_fw(_del_r, _ins_r, _match_rr, _match_end); break;
		default: break;
	}

	struct gaba_score_s *sc = gaba_dp_malloc(self, sizeof(struct gaba_score_s));
	#if MODEL == LINEAR || MODEL == AFFINE
		v2i32_t gc = _add_v2i32(gac, gbc);
		sc->score = (score
			- (int64_t)self->gi * _hi32(gc)
			- (int64_t)self->ge * _lo32(gc)
		);
	#else /* MODEL == COMBINED */
		v2i32_t gc = _sub_v2i32(_add_v2i32(gac, gbc), _add_v2i32(fac, fbc));
		sc->score = (score
			- (int64_t)self->gi * _hi32(gc)
			- (int64_t)self->ge * _lo32(gc)
			- (int64_t)self->gfa * _lo32(fac)
			- (int64_t)self->gfb * _lo32(fbc)
		);
	#endif

	sc->identity = dc > 0 ? ((double)(dc - xc) / (double)dc) : 0.0;
	_store_v2i32(&sc->agcnt, _lo_v2i32(gbc, gac));
	sc->mcnt = dc - xc; sc->xcnt = xc;
	_store_v2i32(&sc->aicnt, _hi_v2i32(gbc, gac));
	_store_v2i32(&sc->afgcnt, _lo_v2i32(fbc, fac));			/* zero for LINEAR and AFFINE */
	_store_v2i32(&sc->aficnt, _hi_v2i32(fbc, fac));
	sc->adj = 0; sc->reserved = 0;

	#if MODEL != LINEAR
		if((uint32_t)(gcnt + 128) > 256) { return(sc); }	/* gap continues from the previous segment */

		/* compensate gap region */
		uint64_t dir = gcnt < 0 ? 1 : 0, sglen = gcnt < 0 ? -gcnt : gcnt;
		uint64_t tarr = gaba_parse_u64(p, lim), tglen = tzcnt(tarr ^ (-dir)) - dir;	/* ramaining gap length */

		/* split path by segments */
		int32_t adj = _gap(self, dir, tglen + sglen) - _gap(self, dir, sglen);
		while(tglen > 0) { uint64_t glen = (&(++s)->alen)[dir]; glen = MIN2(glen, tglen); tglen -= glen; adj -= _gap(self, dir, glen); }
		sc->adj = adj;
	#endif
	return(sc);
}

/* cleanup macros */
#undef _hadd_v16i8
#undef _del
#undef _del_f
#undef _del_r
#undef _ins
#undef _ins_f
#undef _ins_r
#undef _match_core
#undef _match_ff
#undef _match_fr
#undef _match_rf
#undef _match_rr
#undef _nop


/* context initializations */
/**
 * @fn gaba_init_restore_default
 */
static _force_inline
void gaba_init_restore_default(
	struct gaba_params_s *p)
{
	#define restore(_name, _default) { \
		p->_name = ((uint64_t)(p->_name) == 0) ? (_default) : (p->_name); \
	}
	uint32_t zm = ((v32i8_masku_t){
		.mask = _mask_v32i8(_eq_v32i8(_loadu_v32i8(p->score_matrix), _zero_v32i8()))
	}).all;
	if((zm & 0xfffff) == 0xfffff) {
		/* score matrix for erroneous long reads */
		v16i8_t sc = _seta_v16i8(
			1, -1, -1, -1,
			-1, 1, -1, -1,
			-1, -1, 1, -1,
			-1, -1, -1, 1
		);
		_storeu_v16i8(p->score_matrix, sc);
		/* affine gap penalty */
		p->gi = 1; p->ge = 1; p->gfa = 0; p->gfb = 0;
	}
	restore(xdrop,			50);
	restore(filter_thresh,	0);					/* disable filter */
	return;
}

/**
 * @fn gaba_init_check_score
 * @brief return non-zero if the score is not applicable
 */
static _force_inline
int64_t gaba_init_check_score(
	struct gaba_params_s const *p)
{
	if(_max_match(p) <= 0) { return(-1); }
	if(_max_match(p) > 6) { return(-1); }
	if(_min_match(p) >= 0) { return(-1); }
	if(_min_match(p) < -7) { return(-1); }
	if(_min_match(p) < -2 * (p->gi + p->ge)) { return(-1); }
	if(p->gfa != 0 && p->gfb != 0 && _min_match(p) <= -1 * (p->gfa + p->gfb)) { return(-1); }
	if(p->ge <= 0) { return(-1); }
	if(p->gi < 0) { return(-1); }
	if(p->gfa < 0 || (p->gfa != 0 && p->gfa <= p->ge)) { return(-1); }
	if(p->gfb < 0 || (p->gfb != 0 && p->gfb <= p->ge)) { return(-1); }
	if((p->gfa == 0) ^ (p->gfb == 0)) { return(-1); }
	for(int32_t i = 0; i < _W/2; i++) {
		int32_t t1 = _ofs_h(p) + _gap_h(p, i*2 + 1) - _gap_h(p, i*2);
		int32_t t2 = _ofs_h(p) + (_max_match(p) + _gap_v(p, i*2 + 1)) - _gap_v(p, (i + 1) * 2);
		int32_t t3 = _ofs_v(p) + (_max_match(p) + _gap_h(p, i*2 + 1)) - _gap_h(p, (i + 1) * 2);
		int32_t t4 = _ofs_v(p) + _gap_h(p, i*2 + 1) - _gap_h(p, i*2);

		if(MAX4(t1, t2, t3, t4) > 127) { return(-1); }
		if(MIN4(t2, t2, t3, t4) < 0) { return(-1); }
	}
	return(0);
}

/**
 * @fn gaba_init_score_vector
 */
static _force_inline
struct gaba_score_vec_s gaba_init_score_vector(
	struct gaba_params_s const *p)
{
	v16i8_t scv = _loadu_v16i8(p->score_matrix);
	int8_t ge = -p->ge, gi = -p->gi;			/* convert to negative values */
	struct gaba_score_vec_s sc __attribute__(( aligned(MEM_ALIGN_SIZE) ));

	/* score matrices */
	#if BIT == 4
		int8_t m = _hmax_v16i8(scv);
		int8_t x = _hmax_v16i8(_sub_v16i8(_zero_v16i8(), scv));
		scv = _add_v16i8(_bsl_v16i8(_set_v16i8(m + x), 1), _set_v16i8(-x));
	#endif
	_store_sb(sc, _add_v16i8(scv, _set_v16i8(-2 * (ge + gi))));

	/* gap penalties; adj, ofs, gfh, gfv */
	#if MODEL == LINEAR
		_store_adjh(sc, 0, ge + gi, 0, 0);
		_store_adjv(sc, 0, ge + gi, 0, 0);
		_store_ofsh(sc, 0, ge + gi, 0, 0);
		_store_ofsv(sc, 0, ge + gi, 0, 0);
	#elif MODEL == AFFINE
		_store_adjh(sc, -gi, ge + gi, 0, 0);
		_store_adjv(sc, -gi, ge + gi, 0, 0);
		_store_ofsh(sc, -gi, ge + gi, 0, 0);
		_store_ofsv(sc, -gi, ge + gi, 0, 0);
	#else	/* COMBINED */
		int8_t gfa = -p->gfa, gfb = -p->gfb;	/* convert to negative values */
		_store_adjh(sc, -gi, ge + gi, -(ge + gi) + gfb, -(ge + gi) + gfa);
		_store_adjv(sc, -gi, ge + gi, -(ge + gi) + gfb, -(ge + gi) + gfa);
		_store_ofsh(sc, -gi, ge + gi, -(ge + gi) + gfb, -(ge + gi) + gfa);
		_store_ofsv(sc, -gi, ge + gi, -(ge + gi) + gfb, -(ge + gi) + gfa);
	#endif
	return(sc);
}

/**
 * @fn gaba_init_middle_delta
 */
static _force_inline
struct gaba_middle_delta_s gaba_init_middle_delta(
	struct gaba_params_s const *p)
{
	struct gaba_middle_delta_s md __attribute__(( aligned(MEM_ALIGN_SIZE) ));
	for(int i = 0; i < _W/2; i++) {
		md.delta[_W/2 - 1 - i] = -(i + 1) * _max_match(p) + _gap_h(p, i*2 + 1);
		md.delta[_W/2     + i] = -(i + 1) * _max_match(p) + _gap_v(p, i*2 + 1);
	}
	_print_w(_load_w(&md));
	return(md);
}

/**
 * @fn gaba_init_diff_vectors
 * @detail
 * dH[i, j] = S[i, j] - S[i - 1, j]
 * dV[i, j] = S[i, j] - S[i, j - 1]
 * dE[i, j] = E[i, j] - S[i, j]
 * dF[i, j] = F[i, j] - S[i, j]
 */
static _force_inline
struct gaba_diff_vec_s gaba_init_diff_vectors(
	struct gaba_params_s const *p)
{
	struct gaba_diff_vec_s diff __attribute__(( aligned(MEM_ALIGN_SIZE) ));
	/**
	 * dh: dH[i, j] - ge
	 * dv: dV[i, j] - ge
	 * de: dE[i, j] + gi + dV[i, j] - ge
	 * df: dF[i, j] + gi + dH[i, j] - ge
	 */
	for(int i = 0; i < _W/2; i++) {
		diff.dh[_W/2 - 1 - i] = _ofs_h(p) + _gap_h(p, i*2 + 1) - _gap_h(p, i*2);
		diff.dh[_W/2     + i] = _ofs_h(p) + _max_match(p) + _gap_v(p, i*2 + 1) - _gap_v(p, (i + 1) * 2);
		diff.dv[_W/2 - 1 - i] = _ofs_v(p) + _max_match(p) + _gap_h(p, i*2 + 1) - _gap_h(p, (i + 1) * 2);
		diff.dv[_W/2     + i] = _ofs_v(p) + _gap_v(p, i*2 + 1) - _gap_v(p, i*2);
	#if MODEL == AFFINE || MODEL == COMBINED
		diff.de[_W/2 - 1 - i] = _ofs_e(p) + diff.dv[_W/2 - 1 - i] + _gap_e(p, i*2 + 1) - _gap_h(p, i*2 + 1);
		diff.de[_W/2     + i] = _ofs_e(p) + diff.dv[_W/2     + i] - p->gi;
		diff.df[_W/2 - 1 - i] = _ofs_f(p) + diff.dh[_W/2 - 1 - i] - p->gi;
		diff.df[_W/2     + i] = _ofs_f(p) + diff.dh[_W/2     + i] + _gap_f(p, i*2 + 1) - _gap_v(p, i*2 + 1);
	#endif
	}
	#if MODEL == AFFINE || MODEL == COMBINED
		/* negate dh for affine and combined */
		_store_n(&diff.dh, _sub_n(_zero_n(), _load_n(&diff.dh)));
	#endif
	return(diff);
}

/**
 * @fn gaba_init_phantom
 * @brief phantom block at root
 */
static _force_inline
void gaba_init_phantom(
	struct gaba_root_block_s *ph,
	struct gaba_params_s const *p)
{
	/* 192 bytes are reserved for phantom block */
	*_last_phantom(&ph->tail) = (struct gaba_phantom_s){
		.acc = 0,  .xstat = ROOT,
		.acnt = 0, .bcnt = 0,
		.blk = NULL,
		.diff = gaba_init_diff_vectors(p)
	};

	/*
	 * fill root tail object
	 * note on max, mdrop and xd.drop: The max-tracking variables are initialized to keep
	 * values that correspond to the initial vector at p = -1. This configuration result in
	 * the first maximum value update occurring at (i, j) = (0, 0) or the origin.
	 * The max_mask and maximum value are correctly set here and they serve as sentinels
	 * in the dp_search_max and trace_reload_section loop, especially for non-matching input
	 * sequences.
	 */
	int64_t init_max = -(_max_match(p) + _gap_h(p, 1));
	ph->tail = (struct gaba_joint_tail_s){
		/* fill object: coordinate and status */
		.f = {
			.max = init_max,
			.status = CONT | GABA_UPDATE_A | GABA_UPDATE_B,
			.ascnt = 0,    .bscnt = 0,
			.apos = -_W/2, .bpos = -_W/2,
			.aid = 0,      .bid = 0,
		},

		/* section info */
		.tail = NULL,
		.aridx = 0,    .bridx = 0,
		.aadv = 0,     .badv = 0,
		.atptr = NULL, .btptr = NULL,
		.abrk = 0,     .bbrk = 0,		/* no breakpoint set */

		/* score and vectors */
		.mdrop = init_max,
		.ch.w = { [0] = _max_match_base_a(p), [_W-1] = _max_match_base_b(p)<<4 },
		.xd.drop = { 0 /*[_W / 2] = _max_match(p) - _gap_v(p, 1)*/ },
		.md = gaba_init_middle_delta(p)
	};

	/* add offsets */
	ph->tail.mdrop -= 128;
	_store_n(&ph->tail.xd.drop,
		_add_n(_load_n(&ph->tail.xd.drop), _set_n(-128))
	);
	return;
}

/**
 * @fn gaba_init_accumulate_score
 */
static _force_inline
v2i64_t gaba_init_accumulate_score(
	struct gaba_params_s const *p)
{
	int64_t acc[2] = { 0, 0 };
	for(uint64_t i = 0; i < 16; i++) {
		acc[(i & 0x03) == (i>>2)] += (int64_t)p->score_matrix[i];
	}
	return(_loadu_v2i64(acc));
}

/**
 * @fn gaba_init_dp_context
 */
static _force_inline
void gaba_init_dp_context(
	struct gaba_context_s *ctx,
	struct gaba_params_s const *p)
{
	v2i64_t acc = gaba_init_accumulate_score(p);
	double m = (double)_hi64(acc) / 4.0, x = (double)_lo64(acc) / 12.0;
	debug("m(%f), x(%f), xmx(%f)", m, x, x / (m - x));

	/* initialize template, see also: "loaded on init" mark and GABA_DP_CONTEXT_LOAD_OFFSET */
	ctx->dp = (struct gaba_dp_context_s){
		/* score vectors */
		.scv = gaba_init_score_vector(p),
		.tx = p->xdrop - 128,
		.tf = p->filter_thresh,

		.gi = p->gi, .ge = p->ge, .gfa = p->gfa, .gfb = p->gfb,
		.imx = 1 / (m - x), .xmx = x / (m - x),
		.ofs = 2 * (p->ge + p->gi),
		.aflen = p->gi / (p->gfa - p->ge), .bflen = p->gi / (p->gfb - p->ge),

		/* pointers to root vectors */
		.root = {
			[_dp_ctx_index(16)] = &_proot(ctx, 16)->tail,
			[_dp_ctx_index(32)] = &_proot(ctx, 32)->tail,
			[_dp_ctx_index(64)] = &_proot(ctx, 64)->tail
		}
	};
	debug("g(%d, %d, %d, %d), g(%d, %d, %d, %d)",
		p->gi, p->ge, p->gfa, p->gfb,
		ctx->dp.gi, ctx->dp.ge, ctx->dp.gfa, ctx->dp.gfb);
	return;
}

/**
 * @fn gaba_init
 * @brief params can be NULL (default params for noisy long reads are loaded).
 */
gaba_t *_export(gaba_init)(
	struct gaba_params_s const *p)
{
	/* restore defaults and check sanity */
	struct gaba_params_s pi = (p != NULL		/* copy params to local stack to modify it afterwards */
		? *p
		: (struct gaba_params_s){ { 0 }, 0 }
	);

	/* malloc gaba_context_s */
	struct gaba_context_s *ctx = NULL;
	if(pi.reserved == NULL) {
		/* create new context */
		gaba_init_restore_default(&pi);			/* restore defaults */
		if(gaba_init_check_score(&pi) != 0) {	/* check the scores are applicable */
			debug("unsupported score: (%d, %d, %d, %d, %d, %d)", p->score_matrix[0], p->score_matrix[1], p->gi, p->ge, p->gfa, p->gfb);
			return(NULL);						/* cannot be instanciated with the score param */
		}

		if((ctx = gaba_malloc(sizeof(struct gaba_context_s))) == NULL) {
			return(NULL);
		}
		gaba_init_dp_context(ctx, &pi);
	} else {
		/* reuse the previous dp context: fill phantom objects of the context */
		ctx = (struct gaba_context_s *)pi.reserved;
	}

	/* load default phantom objects */
	debug("init BW(%u), proot(%p, %p), root(%p), ctx(%p)", BW, _proot(ctx, BW), &_proot(ctx, BW)->tail, ctx->dp.root[_dp_ctx_index(BW)], ctx);
	gaba_init_phantom(_proot(ctx, BW), &pi);
	return((gaba_t *)ctx);
}

/**
 * @fn gaba_clean
 */
void _export(gaba_clean)(
	struct gaba_context_s *ctx)
{
	if(ctx != NULL) { gaba_free(ctx); }
	return;
}

/**
 * @fn gaba_dp_init
 */
struct gaba_dp_context_s *_export(gaba_dp_init)(
	struct gaba_context_s const *ctx)
{
	/* malloc stack memory */
	struct gaba_dp_context_s *self = gaba_malloc(sizeof(struct gaba_dp_context_s) + MEM_INIT_SIZE);
	if(self == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* add offset */
	self = _restore_dp_context_global(self);

	/* copy template */
	_memcpy_blk_aa(
		(uint8_t *)self + GABA_DP_CONTEXT_LOAD_OFFSET,
		(uint8_t *)&ctx->dp + GABA_DP_CONTEXT_LOAD_OFFSET,
		GABA_DP_CONTEXT_LOAD_SIZE
	);

	/* init stack pointers */
	self->stack.mem = &self->mem;
	self->stack.top = (uint8_t *)(self + 1);
	self->stack.end = (uint8_t *)self + MEM_INIT_SIZE;

	/* init mem object */
	self->mem = (struct gaba_mem_block_s){
		.next = NULL,
		.size = MEM_INIT_SIZE
	};

	/* return offsetted pointer */
	return(_export_dp_context(self));
}

/**
 * @fn gaba_dp_add_stack
 * @brief returns zero when succeeded
 */
static _force_inline
int64_t gaba_dp_add_stack(
	struct gaba_dp_context_s *self,
	uint64_t size)
{
	debug("add_stack, ptr(%p)", self->stack.mem->next);
	if(self->stack.mem->next == NULL) {
		/* current stack is the tail of the memory block chain, add new block */
		size = MAX2(
			size + _roundup(sizeof(struct gaba_mem_block_s), MEM_ALIGN_SIZE),
			2 * self->stack.mem->size
		);
		struct gaba_mem_block_s *mem = gaba_malloc(size);
		debug("malloc called, mem(%p)", mem);
		if(mem == NULL) { return(-1); }

		/* link new node to the tail of the current chain */
		self->stack.mem->next = mem;

		/* initialize the new memory block */
		mem->next = NULL;
		mem->size = size;
	}
	/* follow a forward link and init stack pointers */
	self->stack.mem = self->stack.mem->next;
	self->stack.top = (uint8_t *)_roundup((uintptr_t)(self->stack.mem + 1), MEM_ALIGN_SIZE);
	self->stack.end = (uint8_t *)self->stack.mem + self->stack.mem->size;
	return(0);
}

/**
 * @fn gaba_dp_flush
 */
void _export(gaba_dp_flush)(
	struct gaba_dp_context_s *self)
{
	/* restore dp context pointer by adding offset */
	self = _restore_dp_context(self);

	self->stack.mem = &self->mem;
	self->stack.top = (uint8_t *)(self + 1);
	self->stack.end = (uint8_t *)self + MEM_INIT_SIZE;
	return;
}

/**
 * @fn gaba_dp_save_stack
 */
struct gaba_stack_s const *_export(gaba_dp_save_stack)(
	struct gaba_dp_context_s *self)
{
	self = _restore_dp_context(self);

	/* load stack content before gaba_stack_s object is allocated from it */
	struct gaba_stack_s save = self->stack;

	/* save the previously loaded content to the stack object */
	struct gaba_stack_s *stack = gaba_dp_malloc(self, sizeof(struct gaba_stack_s));
	*stack = save;
	debug("save stack(%p), (%p, %p, %p)", stack, stack->mem, stack->top, stack->end);
	return(stack);
}

/**
 * @fn gaba_dp_flush_stack
 */
void _export(gaba_dp_flush_stack)(
	struct gaba_dp_context_s *self,
	gaba_stack_t const *stack)
{
	if(stack == NULL) { return; }
	self = _restore_dp_context(self);
	self->stack = *stack;
	debug("restore stack, self(%p), mem(%p, %p, %p)", self, stack->mem, stack->top, stack->end);
	return;
}

/**
 * @fn gaba_dp_malloc
 */
static _force_inline
void *gaba_dp_malloc(
	struct gaba_dp_context_s *self,
	uint64_t size)
{
	/* roundup */
	size = _roundup(size, MEM_ALIGN_SIZE);

	/* malloc */
	debug("self(%p), stack_top(%p), size(%lu)", self, self->stack.top, size);
	if(_stack_size(&self->stack) < size) {
		if(gaba_dp_add_stack(self, size) != 0) {
			return(NULL);
		}
		debug("stack.top(%p)", self->stack.top);
	}
	self->stack.top += size;
	return((void *)(self->stack.top - size));
}

/**
 * @fn gaba_dp_free
 * @brief do nothing
 */
static _force_inline
void gaba_dp_free(
	struct gaba_dp_context_s *self,
	void *ptr)
{
	return;
}

/**
 * @fn gaba_dp_clean
 */
void _export(gaba_dp_clean)(
	struct gaba_dp_context_s *self)
{
	if(self == NULL) {
		return;
	}

	/* restore dp context pointer by adding offset */
	self = _restore_dp_context(self);

	struct gaba_mem_block_s *m = self->mem.next;
	debug("cleanup memory chain, m(%p)", m);
	while(m != NULL) {
		struct gaba_mem_block_s *mnext = m->next;
		debug("free m(%p), mnext(%p)", m, mnext);
		gaba_free(m); m = mnext;
	}
	gaba_free(_export_dp_context_global(self));
	return;
}

/* unittests */
#if UNITTEST == 1

#define UNITTEST_SEQ_MARGIN			( 8 )			/* add margin to avoid warnings in the glibc strlen */
#define UNITTEST_GABA_HEAD_MARGIN	( 0 )
#define UNITTEST_GABA_TAIL_MARGIN	( 64 )

#include <string.h>
#include <sys/types.h>
#include <unistd.h>

/**
 * @struct unittest_context_s
 */
struct unittest_context_s {
	struct gaba_params_s const *params;
	struct gaba_context_s *ctx;
	struct gaba_dp_context_s *dp;
};

/**
 * @fn unittest_build_context
 * @brief build context for unittest
 */
#if MODEL == LINEAR
static struct gaba_params_s const *unittest_default_params[8] = {
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(2, 3, 0, 6)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(1, 2, 0, 1)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(2, 4, 0, 3)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(2, 6, 0, 3)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(4, 5, 0, 8)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(4, 7, 0, 8)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(5, 3, 0, 8)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(5, 7, 0, 5))
};
#elif MODEL == AFFINE
static struct gaba_params_s const *unittest_default_params[8] = {
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(2, 3, 5, 1)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(1, 2, 2, 1)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(1, 3, 2, 1)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(2, 6, 3, 1)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(4, 5, 8, 2)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(6, 7, 5, 3)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(5, 3, 8, 2)),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(5, 7, 10, 2))
};
#else
static struct gaba_params_s const *unittest_default_params[8] = {
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(2, 3, 5, 1), .gfa = 2, .gfb = 2),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(2, 4, 5, 1), .gfa = 3, .gfb = 3),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(2, 3, 5, 1), .gfa = 4, .gfb = 2),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(2, 3, 5, 1), .gfa = 2, .gfb = 4),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(4, 5, 8, 2), .gfa = 3, .gfb = 3),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(6, 7, 5, 3), .gfa = 4, .gfb = 4),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(5, 3, 5, 2), .gfa = 3, .gfb = 3),
	GABA_PARAMS(.xdrop = 80, GABA_SCORE_SIMPLE(5, 7, 10, 2), .gfa = 4, .gfb = 4)
};
#endif
static
void *unittest_build_context(void *params)
{
	/* build context */
	struct unittest_context_s *c = malloc(sizeof(struct unittest_context_s));
	c->params = unittest_default_params[rand() % 8];
	c->ctx = _export(gaba_init)(c->params);
	c->dp = _export(gaba_dp_init)(c->ctx);

	if(c->ctx == NULL) { trap(); }
	if(c->dp == NULL) { trap(); }
	return((void *)c);
}

/**
 * @fn unittest_clean_context
 */
static
void unittest_clean_context(void *gctx)
{
	struct unittest_context_s *c = (struct unittest_context_s *)gctx;
	_export(gaba_dp_clean)(c->dp);
	_export(gaba_clean)(c->ctx);
	free(c);
	return;
}

/* global configuration of the tests */
unittest_config(
	.name = "gaba",
	.init = unittest_build_context,
	.clean = unittest_clean_context
);

/**
 * check if gaba_init and gaba_dp_init returned valid pointers to context
 */
unittest()
{
	struct unittest_context_s *c = (struct unittest_context_s *)gctx;
	assert(c != NULL, "%p", c);
	assert(c != NULL, "%p", c->params);
	assert(c != NULL, "%p", c->ctx);
	assert(c != NULL, "%p", c->dp);
}

/**
 * internal tests
 */
/* fetchers */
unittest( .name = "fetch" ) {
	struct gaba_dp_context_s dp = { { 0 } };		/* clear all */
	uint64_t lens[] = { 0, 1, 2, 12, 13, 23, 24, 25, 31, 32, 33, 40, 41, 47, 48, 49, 62, 63, 64 };
	uint64_t ofss[] = { 0, 1, 2, 7, 8, 9, 12, 15, 16, 17, 31, 32, 33, 45, 47, 48 };
	uint8_t const c[72] = {
		A, C, G, T, C, T, G, T, A,
		A, C, G, T, C, T, G, T, A,
		A, C, G, T, C, T, G, T, A,
		A, C, G, T, C, T, G, T, A,
		A, C, G, T, C, T, G, T, A,
		A, C, G, T, C, T, G, T, A,
		A, C, G, T, C, T, G, T, A,
		A, C, G, T, C, T, G, T, A
	};

	#define _clear() { \
		memset(dp.w.r.bufa, 0xff, sizeof(dp.w.r.bufa)); \
		memset(dp.w.r.bufb, 0xff, sizeof(dp.w.r.bufb)); \
	}

	#define _test_buffer(_msg, _ev, _len, _ofs) { \
		uint64_t fail = 0; \
		for(uint64_t idx = 0; idx < (_len); idx++) { if(!(_ev)) { fail++; } } \
		assert(fail == 0, "%s, len(%lu), ofs(%lu)", _msg, (_len), (_ofs)); \
		if(fail > 0) { \
			for(uint64_t i = 0; i < 96; i++) { fprintf(stderr, "%c", decode_table[dp.w.r.bufa[i] & 0x0f]); } fprintf(stderr, "\n"); \
			for(uint64_t i = 0; i < 96; i++) { fprintf(stderr, "%c", decode_table[dp.w.r.bufb[i] & 0x0f]); } fprintf(stderr, "\n"); \
		} \
	}

	/* test seq_a */
	for(uint64_t i = 0; i < sizeof(lens) / sizeof(uint64_t); i++) {
		if(lens[i] >= BLK) { continue; }
		/* fw */
		_clear(); fill_fetch_seq_a(&dp, c, lens[i]);
		_test_buffer("seq_a, fw", *_rd_bufa(&dp, BW + idx, 1) == c[idx], lens[i], 0);

		/* rv */
		_clear(); fill_fetch_seq_a(&dp, gaba_mirror(c, lens[i]), lens[i]);
		_test_buffer("seq_a, rv", *_rd_bufa(&dp, BW + idx, 1) == comp_mask_a[c[lens[i] - idx - 1]], lens[i], 0);
	}

	/* test seq_b */
	for(uint64_t i = 0; i < sizeof(lens) / sizeof(uint64_t); i++) {
		if(lens[i] >= BLK) { continue; }
		_clear(); fill_fetch_seq_b(&dp, c, lens[i]);
		_test_buffer("seq_b, fw", *_rd_bufb(&dp, BW + idx, 1) == c[idx], lens[i], 0);

		_clear(); fill_fetch_seq_b(&dp, gaba_mirror(c, lens[i]), lens[i]);
		_test_buffer("seq_b, rv", *_rd_bufb(&dp, BW + idx, 1) == compshift_mask_b[c[lens[i] - idx - 1]], lens[i], 0);
	}

	/* test seq_a_n */
	for(uint64_t k = 0; k < sizeof(ofss) / sizeof(uint64_t); k++) {
		for(uint64_t i = 0; i < sizeof(lens) / sizeof(uint64_t); i++) {
			if(lens[i] + ofss[k] >= BLK + BW) { continue; }
			_clear(); fill_fetch_seq_a_n(&dp, ofss[k], c, lens[i]);
			_test_buffer("seq_a_n, fw", *_rd_bufa(&dp, ofss[k] + idx, 1) == c[idx], lens[i], ofss[k]);

			_clear(); fill_fetch_seq_a_n(&dp, ofss[k], gaba_mirror(c, lens[i]), lens[i]);
			_test_buffer("seq_a_n, rv", *_rd_bufa(&dp, ofss[k] + idx, 1) == comp_mask_a[c[lens[i] - idx - 1]], lens[i], ofss[k]);
		}

		/* test seq_b_n */
		for(uint64_t i = 0; i < sizeof(lens) / sizeof(uint64_t); i++) {
			if(lens[i] + ofss[k] >= BLK + BW) { continue; }
			_clear(); fill_fetch_seq_b_n(&dp, ofss[k], c, lens[i]);
			_test_buffer("seq_b_n, fw", *_rd_bufb(&dp, ofss[k] + idx, 1) == c[idx], lens[i], ofss[k]);

			_clear(); fill_fetch_seq_b_n(&dp, ofss[k], gaba_mirror(c, lens[i]), lens[i]);
			_test_buffer("seq_b_n, rv", *_rd_bufb(&dp, ofss[k] + idx, 1) == compshift_mask_b[c[lens[i] - idx - 1]], lens[i], ofss[k]);
		}
	}
	#undef _clear
}

/* print_cigar test */
static
int ut_printer(
	void *pbuf,
	uint64_t len,
	char c)
{
	char *buf = *((char **)pbuf);
	len = sprintf(buf, "%" PRId64 "%c", len, c);
	*((char **)pbuf) += len;
	return(len);
}

#define _arr(...)		( (uint32_t const []){ 0, 0, __VA_ARGS__, 0, 0, 0, 0 } + 2 )
unittest()
{
	char *buf = (char *)malloc(16384);
	char *p = buf;

	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "6M4D8M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "6M4I8M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "6M5D8M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "6M5I8M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "2D5M5I8M1I5M5D8M") == 0, "%s", buf);

	free(buf);
}

unittest()
{
	char *buf = (char *)malloc(16384);
	char *p = buf;

	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "8M4D6M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "8M4I6M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "8M5D6M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "8M5I6M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	p = buf;
	gaba_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I5M2D") == 0, "%s", buf);

	free(buf);
}

/* dump_cigar test */
unittest()
{
	uint64_t const len = 16384;
	char *buf = (char *)malloc(len);

	gaba_dump_cigar_forward(buf, len, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "6M4D8M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "6M4I8M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "6M5D8M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "6M5I8M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	gaba_dump_cigar_forward(buf, len, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "2D5M5I8M1I5M5D8M") == 0, "%s", buf);

	free(buf);
}

unittest()
{
	uint64_t const len = 16384;
	char *buf = (char *)malloc(len);

	gaba_dump_cigar_reverse(buf, len, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "8M4D6M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "8M4I6M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "8M5D6M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "8M5I6M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	gaba_dump_cigar_reverse(buf, len, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I5M2D") == 0, "%s", buf);

	free(buf);
}
#undef _arr

/* cross tests */
/**
 * @struct unittest_naive_section_pair_s, unittest_naive_result_s, unittest_pos_score_s
 * @brief result container
 */
struct unittest_naive_section_pair_s {
	uint64_t aid, apos, alen, bid, bpos, blen, ppos;
};
struct unittest_naive_result_s {
	int64_t score;
	uint64_t path_length;
	uint64_t alen, blen;
	char *path;
	uint64_t scnt;
	struct unittest_naive_section_pair_s *sec;
};
struct unittest_pos_score_s {
	int64_t score;
	uint64_t apos;
	uint64_t bpos;
};
struct unittest_naive_section_work_s {
	uint64_t aid, bid;
	uint64_t const *abnd, *bbnd;
	uint64_t appos, bppos;
	struct unittest_naive_section_pair_s *sec, *ptr;
	uint64_t scnt;
};

/**
 * @fn unittest_naive_init_section, unittest_naive_push_section, unittest_naive_finalize_section
 */
static
void unittest_naive_init_section(
	struct unittest_naive_section_work_s *w,
	uint64_t apos, uint64_t bpos,
	uint64_t const *abnd,
	uint64_t const *bbnd)
{
	*w = (struct unittest_naive_section_work_s){ 0 };

	/* calc length */
	uint64_t acnt = 0, bcnt = 0;
	while(abnd[0] < apos && abnd[1] != 0) { debug("acnt(%lu), abnd(%lu, %lu, %lu), apos(%lu)", acnt, abnd[0], abnd[1], abnd[2], apos); abnd++; acnt++; } w->abnd = abnd;
	while(bbnd[0] < bpos && bbnd[1] != 0) { debug("bcnt(%lu), bbnd(%lu, %lu, %lu), bpos(%lu)", bcnt, bbnd[0], bbnd[1], bbnd[2], bpos); bbnd++; bcnt++; } w->bbnd = bbnd;

	w->aid = acnt - 1;
	w->bid = bcnt - 1;
	w->appos = apos;
	w->bppos = bpos;
	w->sec = calloc(acnt + bcnt + 1, sizeof(struct unittest_naive_section_pair_s));
	w->ptr = &w->sec[acnt + bcnt];
	return;
}
static
void unittest_naive_test_section(
	struct unittest_naive_section_work_s *w,
	struct unittest_pos_score_s pos,
	uint64_t aadv, uint64_t badv)
{
	if(pos.apos - aadv < w->abnd[-1] || pos.bpos - badv < w->bbnd[-1]) {
		debug("advance sec, a(%lu, %lu, %lu, %lu), b(%lu, %lu, %lu, %lu)", pos.apos, aadv, w->abnd[-1], w->appos, pos.bpos, badv, w->bbnd[-1], w->bppos);
		*w->ptr-- = (struct unittest_naive_section_pair_s){
			.aid = w->aid, .apos = pos.apos - w->abnd[-1], .alen = w->appos - pos.apos,
			.bid = w->bid, .bpos = pos.bpos - w->bbnd[-1], .blen = w->bppos - pos.bpos,
			.ppos = pos.apos + pos.bpos
		};
		w->scnt++;

		w->appos = pos.apos;
		w->bppos = pos.bpos;
		if(pos.apos - aadv < w->abnd[-1]) { w->abnd--; w->aid--; }
		if(pos.bpos - badv < w->bbnd[-1]) { w->bbnd--; w->bid--; }
	}
	return;
}
static
struct unittest_naive_section_pair_s *unittest_naive_finalize_section(
	struct unittest_naive_section_work_s *w,
	uint64_t *cnt)
{
	if(w->appos != 0 || w->bppos != 0) {
		*w->ptr = (struct unittest_naive_section_pair_s){
			.aid = w->aid, .apos = 0, .alen = w->appos,
			.bid = w->bid, .bpos = 0, .blen = w->bppos,
			.ppos = 0
		};
		w->scnt++;
	}
	memmove(w->sec, w->ptr, w->scnt * sizeof(struct unittest_naive_section_pair_s));
	w->sec[w->scnt] = (struct unittest_naive_section_pair_s){ 0 };
	*cnt = w->scnt;
	debug("scnt(%lu)", w->scnt);
	return(w->sec);
}

/* test if the naive implementation is sane */
#define check_naive_result(_r, _score, _path) ( \
	   (_r).score == (_score) \
	&& strcmp((_r).path, (_path)) == 0 \
	&& (_r).path_length == strlen(_path) \
)
#define print_naive_result(_r) \
	"score(%d), path(%s), len(%d)", \
	(_r).score, (_r).path, (_r).path_length

/**
 * @fn unittest_naive
 *
 * @brief naive implementation of the forward semi-global alignment algorithm
 * left-aligned gap and left-aligned deletion
 */
static
struct unittest_naive_result_s unittest_naive(
	struct gaba_params_s const *sc,
	char const *a, char const *b,
	uint64_t const *abnd, uint64_t const *bbnd)
{
	// debug("a(%s), b(%s)", a, b);
	/* utils */
	#define _a(p, q, plen)	( (q) * ((plen) + 1) + (p) )
	#define s(p, q)			_a(p, 3*(q), alen)
	#define e(p, q)			_a(p, 3*(q)+1, alen)
	#define f(p, q)			_a(p, 3*(q)+2, alen)
	#define m(p, q)			( a[(p) - 1] == b[(q) - 1] ? m : x )

	/* load gap penalties */
	v16i8_t scv = _loadu_v16i8(sc->score_matrix);
	int8_t m = _hmax_v16i8(scv);
	int8_t x = -_hmax_v16i8(_sub_v16i8(_zero_v16i8(), scv));
	int8_t gi = -sc->gi;
	int8_t ge = -sc->ge;
	int8_t gfa = sc->gfa == 0 ? gi + 2 * ge : -sc->gfa;
	int8_t gfb = sc->gfb == 0 ? gi + 2 * ge : -sc->gfb;

	/* calc lengths */
	uint64_t alen = strlen(a);
	uint64_t blen = strlen(b);

	/* calc min */
	int64_t min = INT16_MIN - x - 2*gi;

	/* malloc matrix */
	int64_t *mat = (int64_t *)malloc(
		3 * (alen + 1) * (blen + 1) * sizeof(int64_t)
	);
	if(mat == NULL) { trap(); }

	/* init */
	struct unittest_pos_score_s max = { 0, 0, 0 };

	mat[s(0, 0)] = mat[e(0, 0)] = mat[f(0, 0)] = 0;
	for(uint64_t i = 1; i < alen+1; i++) {
		mat[s(i, 0)] = mat[e(i, 0)] = MAX3(min, gi + (int64_t)i * ge, (int64_t)i * gfb);
		mat[f(i, 0)] = min;
	}
	for(uint64_t j = 1; j < blen+1; j++) {
		mat[s(0, j)] = mat[f(0, j)] = MAX3(min, gi + (int64_t)j * ge, (int64_t)j * gfa);
		mat[e(0, j)] = min;
	}

	for(uint64_t j = 1; j < blen+1; j++) {
		for(uint64_t i = 1; i < alen+1; i++) {
			int64_t score_e = mat[e(i, j)] = MAX2(
				mat[s(i - 1, j)] + gi + ge,
				mat[e(i - 1, j)] + ge
			);
			int64_t score_f = mat[f(i, j)] = MAX2(
				mat[s(i, j - 1)] + gi + ge,
				mat[f(i, j - 1)] + ge
			);
			int64_t score = mat[s(i, j)] = MAX4(min,
				mat[s(i - 1, j - 1)] + m(i, j),
				MAX2(score_e, mat[s(i - 1, j)] + gfb),
				MAX2(score_f, mat[s(i, j - 1)] + gfa)
			);
			if(score > max.score
			|| (score == max.score && (i + j) < (max.apos + max.bpos))) {
				max = (struct unittest_pos_score_s){
					score, i, j
				};
			}
		}
	}
	if(max.score == 0) {
		max = (struct unittest_pos_score_s){ 0, 0, 0 };
	}

	struct unittest_naive_result_s result = {
		.score = max.score,
		.alen = max.apos,
		.blen = max.bpos,
		.path_length = max.apos + max.bpos + 1,
		.path = (char *)malloc(max.apos + max.bpos + UNITTEST_SEQ_MARGIN)
	};
	struct unittest_naive_section_work_s w;
	unittest_naive_init_section(&w, max.apos, max.bpos, abnd, bbnd);

	struct unittest_pos_score_s curr = max;
	int64_t path_index = max.apos + max.bpos + 1;
	while(curr.apos > 0 || curr.bpos > 0) {
		/* M > I > D > X */
		if(curr.bpos > 1 && mat[s(curr.apos, curr.bpos)] == mat[s(curr.apos, curr.bpos - 1)] + gfa) {
			unittest_naive_test_section(&w, curr, 0, 1);
			curr.bpos--;
			result.path[--path_index] = 'D';
		} else if(mat[s(curr.apos, curr.bpos)] == mat[f(curr.apos, curr.bpos)]) {
			while(curr.bpos > 1 && mat[f(curr.apos, curr.bpos)] != mat[s(curr.apos, curr.bpos - 1)] + gi + ge) {
			// while(curr.bpos > 1 && mat[f(curr.apos, curr.bpos)] == mat[f(curr.apos, curr.bpos - 1)] + ge) {
				unittest_naive_test_section(&w, curr, 0, 1);
				curr.bpos--;
				result.path[--path_index] = 'D';
			}
			unittest_naive_test_section(&w, curr, 0, 1);
			curr.bpos--;
			result.path[--path_index] = 'D';
		} else if(curr.apos > 1 && mat[s(curr.apos, curr.bpos)] == mat[s(curr.apos - 1, curr.bpos)] + gfb) {
			unittest_naive_test_section(&w, curr, 1, 0);
			curr.apos--;
			result.path[--path_index] = 'R';
		} else if(mat[s(curr.apos, curr.bpos)] == mat[e(curr.apos, curr.bpos)]) {
			while(curr.apos > 1 && mat[e(curr.apos, curr.bpos)] != mat[s(curr.apos - 1, curr.bpos)] + gi + ge) {
			// while(curr.apos > 1 && mat[e(curr.apos, curr.bpos)] == mat[e(curr.apos - 1, curr.bpos)] + ge) {
				unittest_naive_test_section(&w, curr, 1, 0);
				curr.apos--;
				result.path[--path_index] = 'R';
			}
			unittest_naive_test_section(&w, curr, 1, 0);
			curr.apos--;
			result.path[--path_index] = 'R';
		} else {
			unittest_naive_test_section(&w, curr, 1, 1);
			result.path[--path_index] = 'R';
			result.path[--path_index] = 'D';
			curr.apos--;
			curr.bpos--;
		}
	}
	result.sec = unittest_naive_finalize_section(&w, &result.scnt);

	result.path_length -= path_index;
	for(uint64_t i = 0; i < result.path_length; i++) {
		result.path[i] = result.path[path_index++];
	}
	result.path[result.path_length] = '\0';
	free(mat);

	#undef _a
	#undef s
	#undef e
	#undef f
	#undef m
	return(result);
}

#if MODEL == LINEAR
unittest( .name = "naive" )
{
	struct gaba_params_s const *p = unittest_default_params[0];
	struct unittest_naive_result_s n;

	/* all matches */
	n = unittest_naive(p, "AAAA", "AAAA", (uint64_t const []){ 0, 4, 0 }, (uint64_t const []){ 0, 4, 0 });
	assert(check_naive_result(n, 8, "DRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* mismatch */
	n = unittest_naive(p, "AAAAAAAA", "TAAAAAAAA", (uint64_t const []){ 0, 8, 0 }, (uint64_t const []){ 0, 9, 0 });
	assert(check_naive_result(n, 11, "DRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	n = unittest_naive(p, "GTTTTTTTT", "TTTTTTTT", (uint64_t const []){ 0, 9, 0 }, (uint64_t const []){ 0, 8, 0 });
	assert(check_naive_result(n, 11, "DRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* with deletions */
	n = unittest_naive(p, "TTTTACGTACGT", "TTACGTACGT", (uint64_t const []){ 0, 12, 0 }, (uint64_t const []){ 0, 10, 0 });
	assert(check_naive_result(n, 8, "DRDRRRDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* with insertions */
	n = unittest_naive(p, "TTACGTACGT", "TTTTACGTACGT", (uint64_t const []){ 0, 10, 0 }, (uint64_t const []){ 0, 12, 0 });
	assert(check_naive_result(n, 8, "DRDRDDDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);
}
#elif MODEL == AFFINE
unittest( .name = "naive" )
{
	struct gaba_params_s const *p = unittest_default_params[0];
	struct unittest_naive_result_s n;

	/* all matches */
	n = unittest_naive(p, "AAAA", "AAAA", (uint64_t const []){ 0, 4, 0 }, (uint64_t const []){ 0, 4, 0 });
	assert(check_naive_result(n, 8, "DRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* mismatch */
	n = unittest_naive(p, "AAAAAAAA", "TAAAAAAAA", (uint64_t const []){ 0, 8, 0 }, (uint64_t const []){ 0, 9, 0 });
	assert(check_naive_result(n, 11, "DRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	n = unittest_naive(p, "GTTTTTTTT", "TTTTTTTT", (uint64_t const []){ 0, 9, 0 }, (uint64_t const []){ 0, 8, 0 });
	assert(check_naive_result(n, 11, "DRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* with deletions */
	n = unittest_naive(p, "TTTTACGTACGT", "TTACGTACGT", (uint64_t const []){ 0, 12, 0 }, (uint64_t const []){ 0, 10, 0 });
	assert(check_naive_result(n, 13, "DRDRRRDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* with insertions */
	n = unittest_naive(p, "TTACGTACGT", "TTTTACGTACGT", (uint64_t const []){ 0, 10, 0 }, (uint64_t const []){ 0, 12, 0 });
	assert(check_naive_result(n, 13, "DRDRDDDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* ins-match-del */
	n = unittest_naive(p, "ATGAAGCTGCGAGGC", "TGATGGCTTGCGAGGC", (uint64_t const []){ 0, 15, 0 }, (uint64_t const []){ 0, 16, 0 });
	assert(check_naive_result(n, 6, "DDDRDRDRRRDRDRDRDDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);
}
#else /* COMBINED */
unittest( .name = "naive" )
{
	struct gaba_params_s const *p = unittest_default_params[0];
	struct unittest_naive_result_s n;

	/* all matches */
	n = unittest_naive(p, "AAAA", "AAAA", (uint64_t const []){ 0, 4, 0 }, (uint64_t const []){ 0, 4, 0 });
	assert(check_naive_result(n, 8, "DRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* mismatch */
	n = unittest_naive(p, "AAAAAAAA", "TAAAAAAAA", (uint64_t const []){ 0, 8, 0 }, (uint64_t const []){ 0, 9, 0 });
	assert(check_naive_result(n, 14, "DDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	n = unittest_naive(p, "GTTTTTTTT", "TTTTTTTT", (uint64_t const []){ 0, 9, 0 }, (uint64_t const []){ 0, 8, 0 });
	assert(check_naive_result(n, 14, "RDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* with deletions */
	n = unittest_naive(p, "TTTTACGTACGT", "TTACGTACGT", (uint64_t const []){ 0, 12, 0 }, (uint64_t const []){ 0, 10, 0 });
	assert(check_naive_result(n, 16, "DRDRRRDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* with insertions */
	n = unittest_naive(p, "TTACGTACGT", "TTTTACGTACGT", (uint64_t const []){ 0, 10, 0 }, (uint64_t const []){ 0, 12, 0 });
	assert(check_naive_result(n, 16, "DRDRDDDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);

	/* ins-match-del */
	n = unittest_naive(p, "ATGAAGCTGCGAGGC", "TGATGGCTTGCGAGGC", (uint64_t const []){ 0, 15, 0 }, (uint64_t const []){ 0, 16, 0 });
	assert(check_naive_result(n, 17, "RDRDRDRDRDRDDRDRDDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path); free(n.sec);
}
#endif /* MODEL */



#define UNITTEST_MAX_SEQ_CNT		( 16 )
struct unittest_seq_pair_s {
	char const *a[UNITTEST_MAX_SEQ_CNT];
	char const *b[UNITTEST_MAX_SEQ_CNT];
};

struct unittest_sec_pair_s {
	struct gaba_section_s *a;
	struct gaba_section_s *b;
	uint32_t apos, bpos;
};

/* convert to upper case and subtract offset by 0x40 */
#define _b(x)	( (x) & 0x1f )
static
uint8_t unittest_encode_base_forward(char c)
{
	/* conversion table */
	#if BIT == 2
		enum unittest_bases { A = 0x00, C = 0x01, G = 0x02, T = 0x03 };
	#else
		enum unittest_bases { A = 0x01, C = 0x02, G = 0x04, T = 0x08 };
	#endif

	static uint8_t const table[] = {
		[_b('A')] = A,
		[_b('C')] = C,
		[_b('G')] = G,
		[_b('T')] = T,
		[_b('U')] = T,
		#if BIT == 4
			[_b('R')] = A | G,
			[_b('Y')] = C | T,
			[_b('S')] = G | C,
			[_b('W')] = A | T,
			[_b('K')] = G | T,
			[_b('M')] = A | C,
			[_b('B')] = C | G | T,
			[_b('D')] = A | G | T,
			[_b('H')] = A | C | T,
			[_b('V')] = A | C | G,
			[_b('N')] = 0,		/* treat 'N' as a gap */
		#endif
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b((uint8_t)c)]);
}
static
uint8_t unittest_encode_base_reverse(char c)
{
	/* conversion table */
	#if BIT == 2
		enum unittest_bases { A = 0x00, C = 0x01, G = 0x02, T = 0x03 };
	#else
		enum unittest_bases { A = 0x01, C = 0x02, G = 0x04, T = 0x08 };
	#endif

	static uint8_t const table[] = {
		[_b('A')] = T,
		[_b('C')] = G,
		[_b('G')] = C,
		[_b('T')] = A,
		[_b('U')] = A,
		#if BIT == 4
			[_b('R')] = T | C,
			[_b('Y')] = G | A,
			[_b('S')] = C | G,
			[_b('W')] = T | A,
			[_b('K')] = C | A,
			[_b('M')] = T | G,
			[_b('B')] = G | C | A,
			[_b('D')] = T | C | A,
			[_b('H')] = T | G | A,
			[_b('V')] = T | G | C,
			[_b('N')] = 0,		/* treat 'N' as a gap */
		#endif
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b((uint8_t)c)]);
}
#undef _b

static
char *unittest_cat_seq(char const *const *p)
{
	uint64_t len = 0;
	for(char const *const *q = p; *q != NULL; q++) {
		len += strlen(*q);
	}
	char *b = malloc(sizeof(char) * (len + 1)), *s = b;
	for(char const *const *q = p; *q != NULL; q++) {
		debug("s(%p, %s)", *q, *q);
		memcpy(s, *q, strlen(*q)); s += strlen(*q);
	}
	*s = '\0';
	return(b);
}

static
uint64_t *unittest_build_section_array(char const *const *p)
{
	uint64_t cnt = 0;
	for(char const *const *q = p; *q != NULL; q++) { cnt++; }
	uint64_t *b = malloc(sizeof(uint64_t) * (cnt + 2)), *s = b;
	uint64_t acc = 0; *s++ = acc;
	for(char const *const *q = p; *q != NULL; q++) {
		acc += strlen(*q); *s++ = acc;
	}
	*s = 0;
	return(b);
}

static
struct gaba_section_s *unittest_build_section_forward(char const *const *p, uint64_t pos)
{
	struct gaba_section_s *s = calloc(UNITTEST_MAX_SEQ_CNT + 1, sizeof(struct gaba_section_s));

	uint64_t len = pos + UNITTEST_GABA_HEAD_MARGIN + UNITTEST_GABA_TAIL_MARGIN + _W;
	for(char const *const *q = p; *q != NULL; q++) { len += strlen(*q) + 1; }
	char *a = calloc(1, len); a += pos + UNITTEST_GABA_HEAD_MARGIN;
	uint64_t i = 0;
	while(p[i] != NULL) {
		if(i == 0) {
			s[i] = gaba_build_section(i * 2, a - pos, strlen(p[i]) + pos);
		} else {
			s[i] = gaba_build_section(i * 2, a, strlen(p[i]));
		}
		for(char const *r = p[i]; *r != '\0'; r++) {
			*a++ = unittest_encode_base_forward(*r);
		}
		i++; *a++ = '\0';
	}
	s[i] = gaba_build_section(i * 2, a, _W);
	memset(a, N, _W);
	return(s);
}

static
struct gaba_section_s *unittest_build_section_reverse(char const *const *p, uint64_t pos)
{
	struct gaba_section_s *s = calloc(UNITTEST_MAX_SEQ_CNT + 1, sizeof(struct gaba_section_s));

	uint64_t len = pos + UNITTEST_GABA_HEAD_MARGIN + UNITTEST_GABA_TAIL_MARGIN + _W;
	for(char const *const *q = p; *q != NULL; q++) { len += strlen(*q) + 1; }
	char *a = calloc(1, len); a += UNITTEST_GABA_HEAD_MARGIN;
	uint64_t i = 0;
	while(p[i] != NULL) {
		if(i == 0) {
			s[i] = gaba_build_section(i * 2 + 1, gaba_mirror(a, strlen(p[i]) + pos), strlen(p[i]) + pos);
		} else {
			s[i] = gaba_build_section(i * 2 + 1, gaba_mirror(a, strlen(p[i])), strlen(p[i]));
		}
		debug("a(%p, %p)", a, s[i].base);
		for(char const *r = p[i] + strlen(p[i]); r > p[i]; r--) {
			*a++ = unittest_encode_base_reverse(r[-1]);
			debug("r[-1](%c), *a(%x)", r[-1], a[-1]);
		}
		if(i == 0) { a += pos; }
		i++; *a++ = '\0';
	}
	s[i] = gaba_build_section(i * 2 + 1, gaba_mirror(a, _W), _W);
	memset(a, N, _W);
	return(s);
}

static
void unittest_clean_section(struct unittest_sec_pair_s *s)
{
	if(s->a[0].base < GABA_EOU) {
		free((char *)s->a[0].base - UNITTEST_GABA_HEAD_MARGIN);
	} else {
		free((char *)gaba_mirror(s->a[0].base, s->a[0].len) - UNITTEST_GABA_HEAD_MARGIN);
	}
	if(s->b[0].base < GABA_EOU) {
		free((char *)s->b[0].base - UNITTEST_GABA_HEAD_MARGIN);
	} else {
		free((char *)gaba_mirror(s->b[0].base, s->b[0].len) - UNITTEST_GABA_HEAD_MARGIN);
	}
	free(s->a); free(s->b);
	free(s);
	return;
}

static
struct unittest_sec_pair_s *unittest_build_section(
	struct unittest_seq_pair_s const *p,
	struct gaba_section_s *(*build_section)(char const *const *, uint64_t))
{
	struct unittest_sec_pair_s *s = malloc(sizeof(struct unittest_sec_pair_s));
	s->apos = rand() % 64; s->a = build_section(p->a, s->apos);
	s->bpos = rand() % 64; s->b = build_section(p->b, s->bpos);
	return(s);
}

static
struct gaba_fill_s const *unittest_dp_extend(
	struct gaba_dp_context_s *dp,
	struct unittest_sec_pair_s *p)
{
	/* fill root */
	struct gaba_section_s const *a = p->a, *b = p->b;
	struct gaba_fill_s const *f = _export(gaba_dp_fill_root)(dp, a, p->apos, b, p->bpos, 0);

	gaba_fill_t const *m = f;
	while((f->status & GABA_TERM) == 0) {
		if(f->status & GABA_UPDATE_A) {
			a++;
			debug("update a(%u, %u, %p)", a->id, a->len, a->base);
		}
		if(f->status & GABA_UPDATE_B) {
			b++;
			debug("update b(%u, %u, %p)", b->id, b->len, b->base);
		}
		if(a->base == NULL || b->base == NULL) { break; }
		f = _export(gaba_dp_fill)(dp, f, a, b, 0);
		m = f->max > m->max ? f : m;
	}
	return(m);									/* never be null */
}

static
int unittest_check_maxpos(
	uint32_t id,				/* pos->id */
	uint32_t pos,				/* pos->pos */
	uint64_t ofs,				/* length before head */
	uint64_t len,				/* nr.len */
	struct gaba_section_s const *sec)
{
	int64_t acc = 0;
	debug("id(%u), pos(%u), ofs(%lu), len(%lu)", id, pos, ofs, len);
	while(acc + sec->len < ofs + len) { acc += sec++->len; }
	return((id == sec->id && ofs + len/* - 1*/ == pos + acc) ? 1 : 0);
}

static
int unittest_check_path(
	struct gaba_alignment_s const *aln,
	char const *str)
{
	int64_t plen = aln->plen, slen = strlen(str);
	uint32_t const *p = aln->path;
	debug("%s", str);

	/* first check length */
	if(plen != slen) {
		debug("plen(%ld), slen(%ld)", plen, slen);
		return(0);
	}

	/* next compare encoded string (bit string) */
	while(plen > 0) {
		uint32_t array = 0;
		for(int64_t i = 0; i < 32; i++) {
			if(plen-- == 0) {
				array = (array>>(32 - i)) | ((uint64_t)0x01<<i);
				break;
			}
			array = (array>>1) | ((*str++ == 'D') ? 0x80000000 : 0);
			// debug("%c, %x", str[-1], array);
		}
		// debug("path(%x), array(%x)", *p, array);
		if(*p++ != array) {
			return(0);
		}
	}
	return(1);
}
#define unittest_decode_path(_r) ({ \
	uint64_t plen = (_r)->plen, cnt = 0; \
	uint32_t const *path = (_r)->path; \
	uint32_t path_array = *path; \
	char *ptr = alloca(plen); \
	char *p = ptr; \
	while(plen-- > 0) { \
		*p++ = (path_array & 0x01) ? 'D' : 'R'; \
		path_array >>= 1; \
		if(++cnt == 32) { \
			path_array = *++path; \
			cnt = 0; \
		} \
	} \
	*p = '\0'; \
	ptr; \
})

static
int unittest_check_section(
	struct gaba_alignment_s const *aln,
	struct unittest_naive_section_pair_s const *sec,
	uint64_t scnt)
{
	int rc = 1;
	debug("scnt(%u, %lu)", aln->slen, scnt);
	if(aln->slen != scnt) { rc = 0; }

	for(uint64_t i = 0; i < aln->slen; i++) {

		debug("test, a(%u, %u, %u), b(%u, %u, %u), p(%lu), a(%lu, %lu, %lu), b(%lu, %lu, %lu), p(%lu)",
			aln->seg[i].aid, aln->seg[i].apos, aln->seg[i].alen, aln->seg[i].bid, aln->seg[i].bpos, aln->seg[i].blen, aln->seg[i].ppos,
			sec[i].aid, sec[i].apos, sec[i].alen, sec[i].bid, sec[i].bpos, sec[i].blen, sec[i].ppos);

		if(aln->seg[i].aid>>1 != sec[i].aid) { rc = 0; }
		if(aln->seg[i].apos != sec[i].apos) { rc = 0; }
		if(aln->seg[i].alen != sec[i].alen) { rc = 0; }
		if(aln->seg[i].bid>>1 != sec[i].bid) { rc = 0; }
		if(aln->seg[i].bpos != sec[i].bpos) { rc = 0; }
		if(aln->seg[i].blen != sec[i].blen) { rc = 0; }
		if(aln->seg[i].ppos != sec[i].ppos) { rc = 0; }
	}
	for(uint64_t i = scnt; i < aln->slen; i++) {
		debug("gaba rem, a(%u, %u, %u), b(%u, %u, %u), p(%lu)",
			aln->seg[i].aid, aln->seg[i].apos, aln->seg[i].alen, aln->seg[i].bid, aln->seg[i].bpos, aln->seg[i].blen, aln->seg[i].ppos);
	}
	for(uint64_t i = aln->slen; i < scnt; i++) {
		debug("naive rem, a(%lu, %lu, %lu), b(%lu, %lu, %lu), p(%lu)",
			sec[i].aid, sec[i].apos, sec[i].alen, sec[i].bid, sec[i].bpos, sec[i].blen, sec[i].ppos);
	}
	return(rc);
}

static
char *string_pair_diff(
	char const *a,
	char const *b)
{
	uint64_t len = 2 * (strlen(a) + strlen(b));
	char *base = malloc(len);
	char *ptr = base, *tail = base + len - 1;
	uint64_t state = 0;

	#define push(ch) { \
		*ptr++ = (ch); \
		if(ptr == tail) { \
			base = realloc(base, 2 * len); \
			ptr = base + len; \
			tail = base + 2 * len; \
			len *= 2; \
		} \
	}
	#define push_str(str) { \
		for(uint64_t i = 0; i < strlen(str); i++) { \
			push(str[i]); \
		} \
	}

	uint64_t i;
	for(i = 0; i < MIN2(strlen(a), strlen(b)); i++) {
		if(state == 0 && a[i] != b[i]) {
			push_str("\x1b[31m"); state = 1;
		} else if(state == 1 && a[i] == b[i]) {
			push_str("\x1b[39m"); state = 0;
		}
		push(a[i]);
	}
	if(state == 1) { push_str("\x1b[39m"); state = 0; }
	for(; i < strlen(a); i++) { push(a[i]); }

	push('\n');
	for(uint64_t i = 0; i < strlen(b); i++) {
		push(b[i]);
	}

	push('\0');
	return(base);
}
#define format_string_pair_diff(_a, _b) ({ \
	char *str = string_pair_diff(_a, _b); \
	char *copy = alloca(strlen(str) + 1); \
	strcpy(copy, str); \
	free(str); \
	copy; \
})
#define print_string_pair_diff(_a, _b)		"\n%s", format_string_pair_diff(_a, _b)

static
void unittest_test_pair(
	UNITTEST_ARG_DECL,
	struct gaba_params_s const *params,
	struct gaba_dp_context_s *dp,
	struct unittest_seq_pair_s const *pair,
	uint64_t dir)
{
	#define FMT			"[%s:%d:%s] (%d, %d, %d, %d, %d, %d) { .a = { \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\" }, .b = { \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\" } }"
	#define ARG			MODEL_STR, _W, dir == 0 ? "fw" : "rv", \
		params->score_matrix[0], params->score_matrix[1], params->gi, params->ge, params->gfa, params->gfb, \
		pair->a[0], pair->a[1], pair->a[2], pair->a[3], pair->a[4], pair->a[5], \
		pair->b[0], pair->b[1], pair->b[2], pair->b[3], pair->b[4], pair->b[5]

	/* prepare sequences */
	char *a = unittest_cat_seq(pair->a);
	char *b = unittest_cat_seq(pair->b);

	uint64_t *asec = unittest_build_section_array(pair->a);
	uint64_t *bsec = unittest_build_section_array(pair->b);
	struct unittest_sec_pair_s *s = unittest_build_section(
		pair,
		dir == 0 ? unittest_build_section_forward : unittest_build_section_reverse
	);

	/* test naive implementations */
	struct unittest_naive_result_s nr = unittest_naive(params, a, b, asec, bsec);
	assert(nr.score >= 0);
	assert(nr.path_length <= strlen(a) + strlen(b));
	assert(nr.path != NULL);
	assert(nr.sec != NULL);

	/* fix head segment */
	for(uint64_t i = 0; i < nr.scnt && nr.sec[i].aid == 0; i++) { nr.sec[i].apos += s->apos; }
	for(uint64_t i = 0; i < nr.scnt && nr.sec[i].bid == 0; i++) { nr.sec[i].bpos += s->bpos; }

	/* fill-in sections */
	struct gaba_fill_s const *m = unittest_dp_extend(dp, s);

	assert(m != NULL);
	assert(m->max == nr.score, FMT ", m->max(%ld), nr.score(%d)", ARG, m->max, nr.score);

	/* test max pos */
	struct gaba_pos_pair_s const *p = _export(gaba_dp_search_max)(dp, m);
	assert(p != NULL);
	assert(unittest_check_maxpos(p->aid, p->apos, s->apos, nr.alen, s->a), FMT, ARG);
	assert(unittest_check_maxpos(p->bid, p->bpos, s->bpos, nr.blen, s->b), FMT, ARG);

	/* test path */
	struct gaba_alignment_s const *r = _export(gaba_dp_trace)(dp, m, NULL);

	assert(r != NULL);
	assert(r->plen == nr.path_length, FMT ", m->plen(%lu), nr.path_length(%u)", ARG, r->plen, nr.path_length);
	assert(unittest_check_path(r, nr.path), FMT "\n%s", ARG, format_string_pair_diff(unittest_decode_path(r), nr.path));
	assert(unittest_check_section(r, nr.sec, nr.scnt), FMT, ARG);

	/* calc score */
	int64_t score = 0;
	uint64_t mcnt = 0, xcnt = 0, agcnt = 0, bgcnt = 0;
	for(uint64_t i = 0; i < r->slen; i++) {
		struct gaba_score_s const *c = _export(gaba_dp_calc_score)(dp,
			r->path, &r->seg[i], &s->a[r->seg[i].aid>>1], &s->b[r->seg[i].bid>>1]
		);
		assert(c != NULL);

		score += c->score + c->adj;
		mcnt += c->mcnt;
		xcnt += c->xcnt;
		agcnt += c->agcnt;
		bgcnt += c->bgcnt;
		assert(c->identity >= 0.0 && c->identity <= 1.0, FMT ", identity(%f)", ARG, c->identity);
	}
	assert(score == nr.score, FMT ", score(%ld, %ld, %ld)", ARG, score, m->max, nr.score);
	assert((2 * (mcnt + xcnt) + agcnt + bgcnt) == nr.path_length, FMT ", path_length(%lu), mcnt(%lu), xcnt(%lu), agcnt(%lu), bgcnt(%lu)", ARG, nr.path_length, mcnt, xcnt, agcnt, bgcnt);
	assert((mcnt + xcnt + agcnt) == nr.alen, FMT ", alen(%lu), mcnt(%lu), xcnt(%lu), agcnt(%lu), bgcnt(%lu)", ARG, nr.alen, mcnt, xcnt, agcnt, bgcnt);
	assert((mcnt + xcnt + bgcnt) == nr.blen, FMT ", blen(%lu), mcnt(%lu), xcnt(%lu), agcnt(%lu), bgcnt(%lu)", ARG, nr.blen, mcnt, xcnt, agcnt, bgcnt);

	/* FIXME: test gap counts */

	/* cleanup everything */
	unittest_clean_section(s);
	free(a);
	free(b);
	free(asec);
	free(bsec);
	free(nr.path);
	free(nr.sec);
	return;

	#undef FMT
	#undef ARG
}

unittest( .name = "base", .params = &unittest_default_params[5] )
{
	struct unittest_seq_pair_s pairs[] = {
		/**
		 * Empty segment should be supported but not for now. The logic in the sequence fetcher
		 * does not support segments with len == 0, where the reload condition (link forwarding)
		 * is detected by ridx == 0. It is possible to place compare-zero-branch barriar at the
		 * head of the fill-in functions, but not implemented yet (nor decided whether I will)
		 * because the additional branching path causes an overhead on the execution time and
		 * the code (binary) size.
		 */
/*		{
			.a = { "" },
			.b = { "" }
		},*/
		{
			.a = { "A" },
			.b = { "A" }
		},
		{
			.a = { "A", "A" },
			.b = { "A", "A" }
		},
		{
			.a = { "T", "A" },
			.b = { "T", "A" }
		},
		{
			.a = { "AA" },
			.b = { "AA" }
		},
		{
			.a = { "ACGT" },
			.b = { "ACGT" }
		},
		{
			.a = { "ACGT", "ACGT" },
			.b = { "ACGT", "ACGT" }
		},
		{
			.a = { "ACGTACGT" },
			.b = { "ACGT", "ACGT" }
		},
		{
			.a = { "ACGT", "ACGT" },
			.b = { "ACGTACGT" }
		},
		{
			.a = { "GAAAAAAAA" },
			.b = { "AAAAAAAA" }
		},
		{
			.a = { "TTTTTTTT" },
			.b = { "CTTTTTTTT" }
		},
		{
			.a = { "GACGTACGT" },
			.b = { "ACGTACGT" }
		},
		{
			.a = { "ACGTACGT" },
			.b = { "GACGTACGT" }
		},
		{
			.a = { "GACGTACGT", "ACGTACGT", "ACGTACGT", "GACGTACGT" },
			.b = { "ACGTACGT", "GACGTACGT", "GACGTACGT", "ACGTACGT" }
		},
		{
			.a = { "GACGTACGTGACGTACGT" },
			.b = { "ACGTACGT", "ACGTACGT" }
		},
		{
			.a = { "ACGTACGT", "ACGTACGT" },
			.b = { "GACGTACGTGACGTACGT" }
		},
		{
			.a = { "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT" },
			.b = { "ACGTACGT", "ACGTACGTG", "TACGTACGT", "ACGTACGT" }
		},
		{
			.a = { "AAGGGTCGCCAATTG" },
			.b = { "AAGGGTCGCCAATTG" }
		},
		/* non-matching sequences */
		{
			.a = { "TCTGCGGTAACATCCGCTAC" },
			.b = { "CACGAGGTTCATCGGGTACT" }
		},
		{
			.a = { "CCTACTGAGTTT" },
			.b = { "ACATGGCAGTTT" }
		},
		/* collection of SEGVd sequences */
		{
			.a = { "TAGGTGTACCGAG" },
			.b = { "CCTGTATACCGGC" }
		},
		{
			.a = { "TGTATGCACGGATGATCCAGCTCTAATGGATGG" },
			.b = { "TATACGGGATGGGTGTGATCCACTCTAATGGAT" }
		},
		{
			.a = { "GCAGATGGTCCCTCCTAATGAAGCTCGCTGACT" },
			.b = { "GGTAGAGGTCCCTCCTAATGAAGCTCGCTGACT" }
		},
		{
			.a = { "GAATTTGGCCGGCGGATGCTTATGGATGGTTATTTGGCATGTGACAACG" },
			.b = { "GAATTGCCTGGCGGATGTTATGATGTTTTTTGGCTTGTGACATACGATA" }
		},
		{
			.a = { "GGAGTCCCGAATTGCTAAGGCGTCATGAAATGGTAGACTTGTGCCTAGAGACAT" },
			.b = { "GGAGTCCGAGTTGTAAGGCGTCATGGAATGGAGACTTGTGCCTAGAGACATCCG" }
		},
		{
			.a = { "CGCGGTAACTACGCAACTGACCGTTAAGGTCTCGGTCGTTCTATCCTGACGCGA" },
			.b = { "CGCGGTCACCACGCAACTGACTGTTAAGGTCTCGGTCGTTCTAACTGAGCGACC" }
		},
		{
			.a = { "GGTGGGTTACTGAGCATCCTCCCTAAGTTCACTTTCCGAGATAACCAATCAATCTTG" },
			.b = { "GGGGGGTTTCCTGAGCTACTCCCTCCAGTTCCCCATTCACCGAGATAGCCAATCAAT" }
		},
		{
			.a = { "GTCTTACGAACAGAGACGT" },
			.b = { "GACTTCTCTGGACTGACACGT" }
		},
		{
			.a = { "TGGTAGTAGTCGCATACAGTACCTCTAGTCTTACGAACAGAGACGT" },
			.b = { "TGCTAGTAGTGGCATACAGTAGCCTCCTGACTTCTCTGGACTGACACGT" }
		},
		{
			.a = { "GGATCGAGCAACCTGGTAGTAGTCGCATACAGTACCTCTAGTCTTACGAACAGAGACGTGTAATAGATTTACTCAGCGACGGGAGAAAGTGGTGCAACTAATCTC" },
			.b = { "GGATCGGGCTACTGCTAGTAGTGGCATACAGTAGCCTCCTGACTTCTCTGGACTGACACGTCGGTACTAGATTTACTAGCGACGGGAGAAAAAGTCACATACTCA" }
		},
		/* debugging memory */
		{
			.a = { "GCTCTGACCGTGTTGAGAATATCGGGTTTCGGGCTCATGCCCTACGATGCGTA" },
			.b = { "GGTTGACCGTGTTGAGAATATCGGGTTTCGGGCTCATGCCCTACATGCTTTAA" }
		},
		/* debugging section boundaries; non-matching sequence (eaten as a gap) at the head */
		{
			.a = {
				"G",
				"TAGGACACTCAGAAGAGCTTGGATATTATGGGGAAGCCTAGTTGTCGCGCAAACAAATAACGGTGTCAAAAGTTATAGAGACTAGTAGGCGATCTGACGACGGCGGA",
				"AACACAATTCG",
				"GACCCCTCGCCTATAGCAAGCCGGCTGAAAGGATTGATTGACGAAGGGTTAAAGGGCAGCAATCATTAGGTGTCTATCGGAGGCACTGCCTGGT",
				"GGGGAACG",
				"GCAGTGATATAAATGCCGTTCGCTATTCCCCTTCCGTGCTCGTCCGGGAGCCAGACGTGTATGTAACTGTGTCCAGTGCTAGACCAGTTAGATATTACGGAACCCCTCAAAATCCACTTCTCCGCACGCTACTGGCTGGTGGACTCCCCAGACTTTTACCACGTTCATTGGCCGTAAACCCTCGGACGCGTCACATGCGCTCTATGCTGGATGGGAGTGTGCAGTCCCGATCCGGAATAGCAGTGAACTGACACTTGCAACTCGCTGCGGAATGTCGTGC"
			},
			.b = {
				"T",
				"GTCGGGAACTCAGAAGATCTCTGATATATTTATGGGGAAGCTAGTTGTCCGCAGACAAATAAAGGTGTCAAAAGTTATAGAGACTAGTAGGGATCTGACGACTGGCG",
				"AAACAATTCGA",
				"AACCCCTCGCCTTAAAAGCCGGCTGAAATGATTGATTGACGAAGGGTTCATGGGCAGCAATCCATTTCAGGTGTTATCGGAGCAAATGTCTGGT",
				"GTGGAAAC",
				"GCAGTGATAGAAATCATGTTCTCTATTCCCCTTCCGTGCTCGTCCGGTGCTAGACGAGTTTGTTACTGTGTCCAATGCTAGACCAGTTTGAGATTACTAAACCCCCAAACTCCACTTCTCCGACGCTAAGTGCTGGCGGACTCCCCAGACGTTTACCCCGTTCATTGGCGTAAACCCTCGTGATGCGTGCACACGCGCTCTATACTGCATGGAGTGCGCAGTCCCGAATCGGAATAGCAATGTACCTGTAACATGCAGCTCACTGCGGAATATACGTGCC"
			}
		},
		/* debugging init_fetch (BW = 64) */
		{
			.a = { "GGGGCTAGCAGCATCATTCTGTCCCGGAAA" },
			.b = { "GGGCGCGAGCAGCGAAATTCGGCCAAAAA" }
		},
		/* debugging section boundaries */
		{
			.a = {
				"T",
				"AGGCACTCCCGTACTGGG",
				"TGA",
				"CGGGGGGTTACAAGTCAAATACGGTCCTTCAAAGTGACCTCGACATTAATAAAGTTTAGTCGTCAGTCTAAAGAGTACTAATAAGTTCTGACAAGCGGGGGGGCTAATCTTCTTCCTCGGTTTTTGACAGTGCGTAGTTATGCGCTTAAC",
				"TAGGTATCACATGGGGCGCCCTGTCCGCCTAAATTTTGAATCTGTCAATGCACATGCCAT",
				"AT"
			},
			.b = {
				"T",
				"AGGCACTCTCCGTACATTGTG",
				"TGA",
				"CGCAGGGTTACAAGTCAATGCGGTTCCTTCAAAGGATTCTCACATTATAAAGGGTTGGTCGCGTCTAAAGATGTACTAATAAGCTCTGACACGCGGGGGAGCTAATCGTCTTCCTAGTGTCTTTAACAGTCGCGTAATGATGCGCTTAAC",
				"TTGGGATCACAGGGGCGCTGTCGCCTAAATTTTGAATCTGTCAATGCACATGCCACT",
				"AT"
			}
		},
		{
			.a = {
				"CTCACTGCCACACACGGTTCG",
				"ACCCACTTTAGCATATCC",
				"GCCTTTCG",
				"GCTTGAAGGTCCAAGATGTATTTCAGATTGGGTCCGGCAAAGCAGACCGAGGTCTTTTTGCCACGCTAGTATCGTTGACACTGAACATATCGAGTGTGTGCTCTCAGTCTAGGGGGAGCACTGCCCAATTGCGTAAAGTCGTTAGCTGAAGCTCTAGTCGGTTATTATCAATTTACATTATGGTTCTCAATGAGGTCTGTTGGGCGACCTGCCAAACTGCACTAACTCACATTTGACGCTCTATTATATCCTTATTGGAGTCGGAAGCGAATTAGTAGTATCTGCTCTCGGCTATTAGAGGCCCCAAGGAGCCCAGACGGTATGAGCAAGGGTTAATAGCAATCGTAGAGAAGTGGGAGGGACCTCCTCAAGATTGGTCTCCTTCCCTCTCTTGATTGGCGGTGGGCAGGTTTGTAAATCACCTTCGTTAACTTACCCCTTCTGTGTTGGCGTGGTTAGTCGAACACGCGGATCTTTGGCATAACGT",
				"GCATTATACAGGTATGCGCAAAGATGACATGGAGGTGCTTGGGCGCTATCATCGTAAACATTTACCTAGTTACCCTGCTCGAATGATTAC",
				"TAGTGTTG"
			},
			.b = {
				"GTGTCGCTGCCCTCACACCGCTCG",
				"ACCACTTTAGCATTCC",
				"GCCTTTGG",
				"GCTTGAAGGTCCAAGCTGTATTTCAGATTGGGTCGGCTAAGCAGACTGAGGTCTATTTGCACACGCTAGTATCGGTGACACTTCACATATTCGAGTGTGCGCTCTCAAGTCTAGGGGAGCACGCCGACATTGCGTAAAGTCAGTTAGATAGAAGCTCTAAGTCGGTAATTATTAAGTATCATTATGGGTCTCCAGAGGTCTGTTGGGGACCTGCCAAACTGCACTAACTCACATTTGCGCTTATTATATCACTTATTGGAATCGGAAGTGAATTAGTAGTATCTGCTCTCGGCTATTAGAGGCCCCCATGACTCAGCGCGATGAGACCAAGGTCAAGTACCAATCGTAGAGTAAGGCGAGAATCCTCCTTCAGATCGGCTCTTCCCTCTCTTGATTGGCGGTGGGCAGATTTATAATCACCTTCGTTAACCCAGCCCTTACTGTGTTGTCGATGGTAAGTCGAGACCCCGGATCTTTGGCATAAGT",
				"GCATTATACAGGGATGGCCAAGTTGACACGGAGTGCTGGGCGCATTTTCGTAAACATTTACCTAGATTACCATGCTCGAATCAGATTAAGC",
				"TAGTGGTTG"
			}
		},
		{
			.a = {
				"TCACCTCG",
				"CG",
				"CTCAACGCAC",
				"TTTGCTTAGCATGAAACGTGTGAAGCGTCTAATAGATGCCTCCCCCGTTCCTTCTCGCATACTTCGAAACCTGGAAAGTTTAGGGCTAGAATCTGGCGGCTTAGGCGGTCCTCACTCCAACATGTCACAGGTCTCTATCTATGGTACACCAGCTACCAGCGAATTGTTCATGTCAGCCTGTTAATACAATCGACGTTCTCCGTCTTTCAGTTGAACACCCCGTCTCATGTAACTCCGAGGTCGCCGGGTTTACGCAAGTGCTGTGTAAAGTTTCACTGTGCGTGAAGAGGTACCTACTACCTAAAACTGAGTGCCTCTTACGATATTCGTGTTGTTGCTACCTAGATTTGTCCTACGTGAGATCGTCAAAACCTCTTAGGTGTACTCAGTACTAGCTTCTTCAAAGTTGGCCACCCAGCTCGAAGGCCGGACTAGTGCGCTAGTCACCGGAAAAAGTGTCATGAAATGGGTCGGCTGGAAACACCCTATTCTGATGATCGCAATGCGATGCTTCACGACCTTCCTAATTCGACGTATCAGACCTTGACGGTCCAGAAGTATCACTACGCCACTGCGACTGTGTAGCATAATTCTACTCTCTAGGGTAGCTAAGTCAACTACGATGAACGTGCTCTCATTCGTGAGTCTAATGTTCCACGGACTAGGCGCCATAGATCCGCCCTACACCCATGCGGTATCCACACTAAACGGTTCTAGCCTCACGGTCGTGTTAAATAAAGGGCGTAGCCACGAGATGCCCACAACCTCCTTACAGCCAGGATCCGCAATCCTGCTTGGCAGAATCGAGGTATGCAATAAGCCCATATCGACACCAGTAAGTACGGCGCCTTGAAGTAGCGGGTTATCTTCCTGGCATTCAGATTACGATTGGGCCGATCAATTGCAGCTTCGTACTGTTCTGTCTCTAAAAATTTCTCGGCTAAACAGCTTACTCCAACAACAAATCTCTCGGGGGTGAAGAATTCATCCTGGTGGTTCAGTGACATATTCCTCAGTCCAGAGCAATTTCTCACGCCCCAAGCAAGGTTATCAGGAGGTGTGTTTCACGCCTATTACCTACGCGAAGTGGCCAAAAGACACGTCCGCGAAGCTATTAGGGATTAATAGTCGAAGGCATAAAGGTGCGGTCATATAAGAGCGCGGTATTTTATGGGATCTTGTGGTGTCTCTCGCTCCGCAGCAATGAGTTCCTAAGGTTCAGGTTCAAACCGGCGCTGATGCACCATATGCACATACTCCCGTTTCCGACCGGGTAGGAGACGCCCGTATACCACACGGCTATCGGCAACGCTAGACATGTATAAAACGTTCTGTCCACCGAGCACCGCCACCTCCAGGGCCAAATCCTGGCGCGTGACTCCAGAGTAACTCTTCTGGACAGTTTCGCCTCGACAATTCTCAAATCTCTTGACCTCTCGAATTACTAATCAACAAAGTCCGACTCAAACTCCGACTTGAGTAGCAGCGCTATTCTAGTCTGGCCGGACAAAAAGCCCGTGCTGACCCTAGGTATAGGATGCTGCGTGTGGCTCTAGACCCTAATTACCAAACAGTGGCCATGTCGTTGACTTCTGATTGAAGCAAAGGGTGCCGGTCGGGTCCACAAATGATAGTCGCTAAGGGAACCGTTCCCCGCGGAATACCTGTCCACTGTTGCCCAAGCCATGGGTGCTTGGCCTTCCGTAGGCTTGTCTCCTTTCGGGCCAATCTGTAAGCGAGCTACCAGGTCGAAAGCCGCATCGGCTCGCCCTCCAGCGAATGATGGCAAG",
				"ATTGCGTTTCAGGGGATAAAATATCTACTTGTGCGTTCTGGTCACTACACCAACTCAATATCTCGCCCTATCT",
				"GCAGG"
			},
			.b = {
				"TCACTTCA",
				"CG",
				"CTACACGCAC",
				"TTTGGTTAGGATGAAACGATGTGAAGAGTCTAATAGGTGCCTCGCCTGTTCCCTCACGATACCTCGTAGCCTGGAAGTTGAGAGCTAGAATCTGTGCGTCTTAGGCGGGCCTCACTCTAACACGCACTGCGTCTTATCTATGGTAGTCATCAGCTACAGCTGCAATTGACGCATGTGAGCCTGTTTAATGCAATCGAGGTTCTCCATCTGTCAGTTGAACACCCCGCATGCTTTAACATCCGAGGATCGCCGGGTTTACGCAAGTGCTGTGATAAGTTTCCACTTGGTGAAGAGGTTCTTACGTACCTTAAACAGAGGCCTCTTACCAAAGCTCCGTGTTGTTGCTACCTAGAGTTTCCCCTCGTGAGATCGGTCAAAACCTCTCATAAGGTGTACTCATTACTAACTTCTCAAAGGTTGGCACCCAGCTCGGGAAGGTCGGAGAGTGTACTAGCCACCGGACGAAAGTGTCCATGAAAGAGCGGGCTGGAAAACACCCTCTTCTAAGTGACTTCAATGCGATGCAGCACGACCTTCCTAATTCGACGTATCAGACCTTACGGTCCAGAAGTATCACTACGACACTGGGCTGTATCGCATAATTCTACTATTGAGGGTAGCTAAGTCAACTAGAAGATCGTGGTCTACGTTGTGAGGTCTAATGTTCACGGAGTAGGTAGGCCATAGATTCGCCCTACACCCAGGGGATCCACACTAACGGTTCTAGCCTCAAGGTCGTAAGTTAAATAAAGAGGCGTAGCCAGAGATACCCACTACTCCGTACAGTCAGGATCCTCAATCCTACTTGGCAAATGAGGTATGCAAGAAGCCCTATGACACCAGTAAATACGCTCCCTGAAGTACGCGGGTTTGATCCGTCTGGCATTCAGATTACTATTGGCCGAGCAAACCCAGCTTCGTACTGTTCTGTCTCCGAGAATTCGCGGCTAAACAGCTGATTCCACACAAAATCTCTCGGGGTGAAAAATTCATCCATAGTGGTTCATTGACAATATTCCTCAGTCCATAGCATATTTCTCACCGCCAAGCAAGGTATCCGAGAGGTGTGTTCTCACGCTCATATTCCTACGCGAAGTCGGCCAATAGAGACACGTCCGCGAGCTATCAGGGATTAATAATCGAGGGCATATAAGGTCGCGAGCAATATAAGAGCGCGGTATTTAATGGAATCTTGTCGTGTCTCTCGTCTCGCAGAAGTGAGTTCTAAGGTCCAGGCTTCAAATAGGGCTGAGCACCATATGCACATCACTGCCCGTTTTCGACTGGGTAGGAGACGCACGTATCCACACGGCTATGGTCCAACGCTGGACATGTAGAAAACGTTCTTCCCACCGAGCACCGCACCGCCAGGGCCAATCCTGGCGCTCACTCCAGAGCAAACCTTCTGGCCAGTTCGCCTCGACAATTACTCAAATCTCTTGACCTCTCGAAGAGTAAGCAACTTAAGTCCGACTCAAATCCGACTTGAGTAGCGCGCTATTCTAGTCCTCGCCGGTTCAAAAAGCCCTGATGACCCTAGGTATGGATGCCGCGTGTGGTTAGACTCTAATTTACCAAACAGTGGCCATGATCTTGACTTCTGATTAAAGCAAAGGGTGCCGGTCGGGTTCCAAAATGATGAGTCCATAGGGAACTCGGGTTCCCCGAGGAATACGTGAATCCACTGTGCCAAAGCCATGGGTGCTTGGCCTTCCTTAAGGCTTTGTCTCCTTTCGGGCCAATCTGTAAGCGAGTTACCAGGTCAAAAGCCTCATCGGCTCGCCCTCCAGCAATGATGGCAAG",
				"GTTGCGTATCAGGGATAAAACATCTACTTGATGCGTTCTGTCACTAACTTCAACTTCAATATCCGCCCTATC",
				"TCAGA"
			}
		},
		{
			.a = {
				"TGC",
				"AC",
				"GG",
				"AAAATCTATGAGTCCTCTCCAATCCTGGGCCGCCGATCCGTCCTTGGAGTGGTCTCCGTTTCCAGTCCAGAATCTGTCTCACTTGTATCAGACTATGCCTTCTTGGCCTATTCGGACTTGTGGACTTAAATACGTTAACTCAGGGAGATGGGACGCAGCAGAAAGTTGAACCATTTAGTTACGAGACGGTTCTTCGTAATCGACCTACCTCCCTCTTTGGGCTCAATTAAGCGGTTGTCTTGCACGTACCGTGCCTCAGGACTGTGGCGTCCCGCGAAGGGCCAGGTATGGTTCGGAAAGTCTATAGTAACGCTTACGAAAGATTGGTGGCCTACGATCCTTTCGGCTCCACTACGGTGACTCAAGGACTCAGCCATGAATAGGATCGGTGTGTTGCCTGAAGTCATCGTCGTTATGCCAATTTGGCTTACACATGAGTGGCACACCTGGTTGTGACTACACCCCATCGTAAAAGACTTACTTTGCATAAAACCCTCCAAAAAACGAGCAGTTCGCGGATGCTAGCTAGTAGAAGCTGTAGGTCTGGGATCTAATGCGGTCGCCTTGGGGGAAGGTAACATGTTTGGAAAAGTCTGCGATAGTAACTAGCCCTACGTTTTTTCGACACTAATTATTCGATTCCTCGCACAAATCTAGTCCAGATGTGTGGTCTAACGGCTGCTGTGGCGTTGCTCTTTTCGGGCGTGACATTCGCACTGTTCCAAGATCTATCTCTGTTTTGGCAACCTTAGACTCTTACGCATACGGCCCATTCTTTCTGCCTTAAGTATCCAACCGTGTATGAAGTTTGGAACTACGACCCTCGATAACAAAGACAAATGGAGTTTACTGTACAGCTTCCTGCCGGGGTCCCCTCCGCGCTCGAGAGCTTGGCAAAGGGCTTTCTAGCCATCGCCACTCCTATAATAGCTACGTGCTATAAGCGTTGAGACCACGAAGGGCGTACCCACACCGGACCGGAGTCGATGCGCTTTGTCTCCACGCACCCGAT",
				"GTACTTTACATCACACACGATCGTGTAATATGTTTGGGTATGCTACAATGCGTTCCGCGGCCAAATTGCGGAATTACTCAGGGCCTTTACCGGCGGACAGTTGAGGATCCTTTCCAAAACCCT",
				"TCCTTCCGAAGTGTGTAGTATAAAAACCGAACCCAAGGTTGCGAAATCGGCGGCGGGTTATTCAACCTTAGCCACGACTTTACGCTATGATCTACAGGTCGCTACAGCCCACCCTGATGGTCTTGCTATGCGTAGACACGCGTGTGAAGGGCAGAGTTGGCTAAAGACGGTTATCTCGTTAGCAACGTCACGAGTAATCTCAGCTACCGAAGGCTGTCCAGTGGCCGCCTAGTTATTCACCATAGGGCAACGTTAACGCGCGAACATACTCTAATGTCTTCACATCAATGGAATACTACCTGATGTCCAAAGACGAGTACCTCCGGGCGTGCCGCTCGCGTCATTGTCCATGTATCGCGAACGCTTCCCGTATTAGCTATTTCGTCTGGGGTTAAGTCTCTGCTCCCTCCGAATTAAGGTAAGTCGGAACCTATGCAAACCGGGGCGGTACAATTCGCGGTTCTGAATACACTGTCTCATCCTGCGTTACCTCGGTGGGGTTATCTACCAAGGTATCAGCCTAGCGGACAGATGTCGAGTAGCGCGCGGTATGAGATAAAAGAACACTATGACTTGCAATAATGCCGACGACCACGAGCACTCGGAAAACTACGTGTCGAAGGCTTTTGAATTCGCCGCCTCGCTCTTGGTGACTTCGTTGTGCAATCTAGGCATCCTAACTGCCGCGAACTTTATATATAGCCTCGAAGACGTGTCTCCTGAGCCCCACAAAGGTTATTTAACTGAAGCCCTGAATCGTTGATAGGTGCGCCGGACGATTGAGAATTGCGCTGAAAGTTTGACATCCCGACACTACGTTTCATAATAGGTGAATCTAGAGTCGACGCTCCTTCGTGCTATTGATCCGGTGAGATGCCAATCTGGCCAGACGTTGGCTATGTTCTGTAATATAACTCGTAGAGGAAACTCACATGTGGTGCTCGTGCCAACGCTACTCGGGGTGCGCTAACCTGTTAAAACCTGATTTACCGTAAAGCAGCCACATAAAACCCATCTGAATCGGCTGAAGGGAGTAGTAGTCTACGCCACCTGTCTGTTTGGCGAAGGTGAGATGCCTGCATAGTCCCCCATTCCTCTGACGGCGATTTCCTCCCATTCAAGTTTAGCTACATAGCACAAGCTTTCTGACGATCTGCCTAGCACCGGCTGCGACCTGCCGGTGCACTCGTCAACCCATGTGACGCTCTAGCGCTTCAGGGCTTCAAGAACAGCCAGACCATGGACCCCATTCGGTGTCTTGTCTACCCGGCCAACCCGGCCGTTAGCTTAGGTTCAATCGCGTGGATAGGCTAAGAGAAGAGAAGCAGACCGTAAGCTGTCTTCACGGCGAGAAGGGGCGTGAGGAATAGTCGGTGTCACTCTTCCGTTTGGCGATTCTTTGACATTTAAATCTTGACGTAACCGGACTTAGCCACATCCTCAACTTCACCACGCTAGCTCCGTGGATCCACTACAAACACGCGCGGCATCGATCTCCGGTGTGGCGGCTTTGACTGTATGCCTAGGCTAAGGTAAAGTGGGCAGGTAGTTACATCGCTGCTATGCATTTCCCGCTAATGATTCTGGGACACAACTCCGCATGCCAACATAAGTCGCAACCATGAACGCTTCGCGTCGAGCAGATAGCAGGCGCTGCAAGTTTCCCAGTCACCTCCTTGTCAATGACTAGGACCCGTGCGTCTAAACTTTGAATTGTCCTTTTACCTATAGGCGGGGACGTAGTTGCGTTCGGCGTGCTCCGTTACCGCGACGTACCGTACCGTGATAGGCTAAGACTTTGTAAGAGTGCTTGCTAGGAAACTCCGCTCAGTAGAGGGGCTCAGTTGCCTCGAATGGACACTTCCTTTCCACTCTGACTGAAACATTTTATATGCAGCTTGAACCTTCCGTCGTGTGCGGCCGATAACCACTACCCCAGGTGTCCGAATGGCTACTAGACGCCG"
			},
			.b = {
				"TGC",
				"CC",
				"GG",
				"TAAAATCTATTAGTCCTAACCCAATCCTGGGCCGCCTGATCCGCCTTCCAGTGTGTCTCCGTATCCAGTCCAGAATCATGTCTCCCCTTGTATGAACTTGCTCTGCTTGGCCTATCGGACTTGTGGGAATTAAAATTACGTTAGCTCAGGGAGCTCGGCATCCAGTAGAAATGGAACCATTTAGTTACGAACGAGTTCTTCGATTCTACTACCCCTTCTTTGGGCTCAAGTTCAAGCGGTTGTCTTGCACGTACTCCGTACTCTCAGGATTGGACGACCCGCGAAGGGACACGTATGATGCGAATGTTTATAGTAACGCTTACAGAAAGATTGGTGGCCTACGATTCTTTCGGCTGCCAATACGGCGAGCTCAAGGACTCAGCCATAATGGTTCGGTGTGAGGCCGGAAATCATCGTCCTTAGACCCATTTGGTCTTACAGGTGAGTGCCCCACCTGGTTGTGACACACCCGATCGTAACAGACTTATTCGCATAAAACCCTCCAAGATACAGAGCGGTTCGCGCGTTGCTAGCTTAGTAACGCTGTAGGTCTGAGATCTAAAGCGGTCCCTTGGGGGAAGGTCAAATGTCTGGCAAAGTCGGCGATAGTAAGCTTCGCCGCTAGTATGATGCGAACAACGGAATTATGGAGGTCCTGGCCAACTCTAGTCCAGATGTGTCGTCTCCGGCTGATGTGGCGTTGCTTTTTCGGCGTGACATCGCAGTGGTCCAAGATCTATCTCTGTTTTGCCACCTTAACTAGTTACGCATACGGCCCATTCTTTCTGCCCAAATATCCAACAGATGTATGAAGTTTGGATCTCCCGCCCCTCAGTAACAAAGACAACTTGGATTTACTGTACAGTTCTGCCGGGGTGTCCCTCCGCGCTCGAGGGCTTGGGCAAAGAGCTTTCTAGCCATCGCCATACCTTATAATAGCTGCATGGCATAAGCCAGAGACGACAAAGGGCGTGCCCACACCGGACGGAGTCGATGCGCTTTGTCTTCCACGTACCCCAT",
				"GGACTTTACATGCACACAGACTTGAGTAATCAGTTTGGCTTATGCTACAATGAGTTCCGCGGCCAAATTGCGGAATACTCGGTATCTGTACCGGCGGACAGTTGAAGAGCCTTTTCCAAACCCT",
				"ACCTTACGAAGTCGTGTATTATAAAACCGAACCAACGGTTTGCGTAACGGGCCGGCGGGTTTTTCAACCTAGCCACGACTTTCGCTATGATCTATAGGTTCGCTACATCCCCGCCTGAGGTCTTGTATGCGTAGAACGCGCGGAGGGGCATAGTTGGCAAAGAGGATTAGCTCGCATAGCGACGAACACGAGTAATCCAGGGACCGAAGCCGTCCCAGTGGAAGCCTAGTTCAGTCCCACAGGGGAACGTTAACGCGCGAACATCTCTCTAAGTCTTCACCTCGATGGAATTACACCCTGATGCCAAAGACGAGTACCTCCCGGCGTGGCGCTCGCTCATTGTCCATGTACCGCGAGCGCTTCCCTTATTAGCTATTTCGTCTGGGGATTAAGTCCCTGCTCACCTCCCGAATTAAAGTACGTCGAAACCTATGGCAACCGGGGCGATTACAATTCGCGGTTCTTAACCACTGTTCATCCAGCGTTAGCTCGGGGGGATATCTACCAAGGTATCATGCCTAGCGGACAGATGTCGAGTAGAGCGGTATGAGATAAAAGAAACTAGGACTTGCAAAACGCCGACGACCACCAGCCCTCGGAAAATACGTTCGAAGCCTTTTAACCGCCCGCCCGCTCTTGGTGACTACTATGTGCATCTAGGCATCCTAAACTGCCGCGAACTTTGTATATTACCCGAAGATGTGTTCCTGGGCCCCACCAGGATTATTTAACTGAAGCCCTGAATCGTTGATAGGTGGCGCGGCGATTGCGAATTAGCGCTGTCAATTTGACATCCCACACTACGTTTCATAATAGGTGAATCTAGAGTAGACGCGCCTTCGTGCCTAATTAGATCGGCGTCGCGCTGCCCAATCTGGACACAGACGTTGGCTCTGTCTGTAATATAAACACGTAGGAGGAAACTTCATTATGGTCCGATCAGTGCCATCGCTACTCGGGGTGCGCCAACCTGGTAAAGCCTGGTTTACCGTAAGCGGCACATAAGCCATCTGAATCGGCTGAAGGGAGAGTAGTCTACGCCACATGTCTGTTTGGCTGAGGTGAGATGCCTGCCTAGGGCCCCCCATTCCCCTGACGGCAAATTTCCTCCCATTCAAGGAGACGTACATAGCACAAGCTTTCATAGACGATCTGCCTACCACCGGCGGCGACCTACAGGTGCACTCGTAAACCCATGTGAGCTCTAGCGCTTAAGGGCTCAGGAACCAGACAGCACATAGACCCCATCCGGTGCCTTGTCTACCTCGGCCAACCGCCCGTTAGCTAGGTTCAATCGCGTGGATGGATAGGCCACGAGAAGCGTACCGTAGCTTTTCACAGCTGAAAGGGGCGTGAGGAGTAGTCGGTGTCACTCTTCCGTATGGCGATTCTCTGACATTTAAATCTGAGTAAGGGACTTAGTCACATCCACAACTTCACCACGCTAGTCCTTGATCCACTACAAAACACACTGCGAGTTTGATCTCCGGTGTGGGCGGTTTTGTCTGTATGCCTAGGCTAGGGTAGAGTGGGCCGGTAGTTGTACATCGGCGCTATGATTTCCCGCTACTGACTTCTGGACACAACCCGCATGCCAACATAAGTCGCAACCATGAAGCTTCGCGTCGAGCGATAGCAGGAAATGCAAGTTTCCCAGTCACCCCATGTCAAACGAATAGTACACGTGCGTCTTAAACTTTGAATTGTGCTATACCTATAGGCGGGGACGAGTTGCGTGGGCGTGCTCCGTTACCCGCGAACGTACCGTCCTTGATAGGCTAAGACTTTTTAAGAGTGCTGCGGTAGAAACTCCGCTCAGTAGAGGGCCTCAGTTGCTCGAATGTACCATTACTTTTCATTCTGACTGAAAAATTTTAACGCGCAGTGAACCTTCCGTCGTGTGCGGCCGATAACCACTACCCAGGTGTCACGAATGGCCACTAGACGCTTG"
			}
		},
		/* debugging reverse fetch */
		{
			.a = { "C", "TCAGCTTAC" },
			.b = { "C", "TTGGTTAC" }
		},
		{
			.a = { "GTCCTGTTACACCCCAGGCGACGGGAGT", "CGGCATAGGTTTTACACCGTTCAATGGTCCTGAGT", "CTGATCTGCATTG" },
			.b = { "GTCCTGTTACACCCCCAGGCGACCGGGAGT", "CGTCATAGCGTTTTACACCGTTCAATGTCTCTTAGGT", "CTGTTTCATG" }
		},
		/* debugging combined */
		{
			.a = {
				"CTCACTGCCACACACGGTTCG",
				"ACCCACTTTAGCATATCC",
				"GCCTTTCG",
				"GCTTGAAGGTCCAAGATGTATTTCAGATTGGGTCCGGCAAAGCAGACCGAGGTCTTTTTGCCACGCTAGTATCGTTGACACTGAACATATCGAGTGTGTGCTCTCAGTCTAGGGGGAGCACTGCCCAATTGCGTAAAGTCGTTAGCTGAAGCTCTAGTCGGTTATTATCAATTTACATTATGGTTCTCAATGAGGTCTGTTGGGCGACCTGCCAAACTGCACTAACTCACATTTGACGCTCTATTATATCCTTATTGGAGTCGGAAGCGAATTAGTAGTATCTGCTCTCGGCTATTAGAGGCCCCAAGGAGCCCAGACGGTATGAGCAAGGGTTAATAGCAATCGTAGAGAAGTGGGAGGGACCTCCTCAAGATTGGTCTCCTTCCCTCTCTTGATTGGCGGTGGGCAGGTTTGTAAATCACCTTCGTTAACTTACCCCTTCTGTGTTGGCGTGGTTAGTCGAACACGCGGATCTTTGGCATAACGT",
				"GCATTATACAGGTATGCGCAAAGATGACATGGAGGTGCTTGGGCGCTATCATCGTAAACATTTACCTAGTTACCCTGCTCGAATGATTAC",
				"TAGTGTTG"
			},
			.b = {
				"GTGTCGCTGCCCTCACACCGCTCG",
				"ACCACTTTAGCATTCC",
				"GCCTTTGG",
				"GCTTGAAGGTCCAAGCTGTATTTCAGATTGGGTCGGCTAAGCAGACTGAGGTCTATTTGCACACGCTAGTATCGGTGACACTTCACATATTCGAGTGTGCGCTCTCAAGTCTAGGGGAGCACGCCGACATTGCGTAAAGTCAGTTAGATAGAAGCTCTAAGTCGGTAATTATTAAGTATCATTATGGGTCTCCAGAGGTCTGTTGGGGACCTGCCAAACTGCACTAACTCACATTTGCGCTTATTATATCACTTATTGGAATCGGAAGTGAATTAGTAGTATCTGCTCTCGGCTATTAGAGGCCCCCATGACTCAGCGCGATGAGACCAAGGTCAAGTACCAATCGTAGAGTAAGGCGAGAATCCTCCTTCAGATCGGCTCTTCCCTCTCTTGATTGGCGGTGGGCAGATTTATAATCACCTTCGTTAACCCAGCCCTTACTGTGTTGTCGATGGTAAGTCGAGACCCCGGATCTTTGGCATAAGT",
				"GCATTATACAGGGATGGCCAAGTTGACACGGAGTGCTGGGCGCATTTTCGTAAACATTTACCTAGATTACCATGCTCGAATCAGATTAAGC",
				"TAGTGGTTG"
			}
		},
		{ .a = { "GGCATCAGGTTGCGACC", "AACGACAGATCACTGGTGTAACTTAAC", "G", "AGGCATTACTTGTAG" }, .b = { "GGCATCAGCTTGCGACC", "AAGCGCAGGATCACTGGTGTGGAAACTGTT", "G", "AGGCTGTTTACTCGTAG" } },
		{ .a = { "GTAGT", "TCAATATTAGACGATGGCTGCTTTCCTCAG", "TTCATCTGC", "GTAGGGAAGGCCCCTTAATGCGCGGATCGATATT" }, .b = { "GTAGT", "TCCATATTAGACGAAGGCTGCGGTCTCG", "TTATCCTGC", "GCAGGAGAGGCCCCTTAATGCGTGGATCGAATATG" } },
		{ .a = { "ACAAAAAGATCGGTC", "ATGGGACGTTGACAGCGGGGTGATCGGTCTTGGCAA", "GTACCATGTA", "CGGAATACCGCGCAATCCTTTGTAAATGTAGAGACTTCAATCTACTTTAGGTC" }, .b = { "ACTAAAAGTCGGAC", "ATAGACGTTCACAGCGGTGATCCGGACTTGGCAA", "GTACCATGTA", "CGGAGTACCACGCACTCCTTTATAAAGTAGAGACTTAAATCCTACTTGAATGGC" } },
		{ .a = { "CTCACTGCCACACACGGTTCG", "ACCCACTTTAGCATATCC", "GCCTTTCG", "GCTTGAAGGTCCAAGATGTATTTCAGATTGGGTCCGGCAAAGCAGACCGAGGTCTTTTTGCCACGCTAGTATCGTTGACACTGAACATATCGAGTGTGTGCTCTCAGTCTAGGGGGAGCACTGCCCAATTGCGTAAAGTCGTTAGCTGAAGCTCTAGTCGGTTATTATCAATTTACATTATGGTTCTCAATGAGGTCTGTTGGGCGACCTGCCAAACTGCACTAACTCACATTTGACGCTCTATTATATCCTTATTGGAGTCGGAAGCGAATTAGTAGTATCTGCTCTCGGCTATTAGAGGCCCCAAGGAGCCCAGACGGTATGAGCAAGGGTTAATAGCAATCGTAGAGAAGTGGGAGGGACCTCCTCAAGATTGGTCTCCTTCCCTCTCTTGATTGGCGGTGGGCAGGTTTGTAAATCACCTTCGTTAACTTACCCCTTCTGTGTTGGCGTGGTTAGTCGAACACGCGGATCTTTGGCATAACGT", "GCATTATACAGGTATGCGCAAAGATGACATGGAGGTGCTTGGGCGCTATCATCGTAAACATTTACCTAGTTACCCTGCTCGAATGATTAC", "TAGTGTTG" }, .b = { "GTGTCGCTGCCCTCACACCGCTCG", "ACCACTTTAGCATTCC", "GCCTTTGG", "GCTTGAAGGTCCAAGCTGTATTTCAGATTGGGTCGGCTAAGCAGACTGAGGTCTATTTGCACACGCTAGTATCGGTGACACTTCACATATTCGAGTGTGCGCTCTCAAGTCTAGGGGAGCACGCCGACATTGCGTAAAGTCAGTTAGATAGAAGCTCTAAGTCGGTAATTATTAAGTATCATTATGGGTCTCCAGAGGTCTGTTGGGGACCTGCCAAACTGCACTAACTCACATTTGCGCTTATTATATCACTTATTGGAATCGGAAGTGAATTAGTAGTATCTGCTCTCGGCTATTAGAGGCCCCCATGACTCAGCGCGATGAGACCAAGGTCAAGTACCAATCGTAGAGTAAGGCGAGAATCCTCCTTCAGATCGGCTCTTCCCTCTCTTGATTGGCGGTGGGCAGATTTATAATCACCTTCGTTAACCCAGCCCTTACTGTGTTGTCGATGGTAAGTCGAGACCCCGGATCTTTGGCATAAGT", "GCATTATACAGGGATGGCCAAGTTGACACGGAGTGCTGGGCGCATTTTCGTAAACATTTACCTAGATTACCATGCTCGAATCAGATTAAGC", "TAGTGGTTG" } },
		{ .a = { "GCTAGTTACAAACCCG", "AGAATCCGAAG", "TTCAGCGGAGTA", "TGTAAAAGACCCGACGGACCCCTGTCGCATTTCTTCTATTCCCTCGCTTGAATTTGCCTCAGCCCAGAAGACTCAATCTGTTTGGGATTTCAGGCTGAATGAACCCCAGATGTGACACCTAGTGACGGCGCTGGATTCTCCGTAAACACGACATAAAGTCAAAGGACCCATGTGGTGAATGACTTCTATGCTTGCCGGTAGGGCACCTCCGTCTGAAGATAGGTATCGCAGAGCCGGGCTTAGATCCATCTTTGTTCAGTAATGCAATTGGGGCTCGACCCCACTTAATGTTTGCATGAGAGAGTTTCTCGCGTGCGGCGCCCATTAGCAACTCTGAGAATGCCCTGACAGGCTGAGTTTTTAGATCCGCCCTGTCCCCGCTCGGTTCGGAGCCAGACTGGAAGGTTTAACGTGACGTGACTGTATCAAATTATCGGACGAGATCCATATTGGGCGCGGCGCAGACAGCCCCTCAAGCTTATCGCGAGAGTGAAAACATTCAAACCCAGATGTAAAGAGGGGAGGATTAGGCGTCTAGTGTGAGGATAACGGTGCCGACGGTTAGCTGTTAGCACAAAAACTTCCAGAATGGTGGCTCAGAACGCGGCGGGTCCATCTTCGCTCTCTGGTCTTGCAATGGCACTCGGGCCTTGGCGTAACTCAGTCTTACCGCTCGCTAGATCGCGAAGCCCGGAGTTCAATGCGTTTGTCGGGATCCTTGACAGCGGAAGATACCTAACTCTCTTAAATAGTCTGGCAGCGGAGTCCTCGGGCTGATGCGGTCTGTCATTCCCGTATATCGCTTACATACCGATGTGTTGTATAGTCGTATGGGCGACGAGTGTGAGCACTTCGTACAGCCCTTTACTGTTATGTCGCGGGCCAACACCGGGTCCGGTGAGCACTATGCCAATGATTAAAGTTCAA", "GTCGGTTTGTGAAACATTAGGAGGGCCTGGTACATTTGTGCCGCCTGAGGGATCTCGACAGACAAATAGGGTAAACAGGGCATACGCAGGACGCCTCATT", "CGTCTACCGAGTTGACAATGCACGCAAGTTCAGACCATGTTCGGTGTTTGGATTGATGGGGTGATTCAGTTCACTATGTGCCGAACATTTCCGTCGTGGCGCTGGGGCCGAGTGGATGTTCCAAAGGTAACGCACTGTATCACCTTCGGACTCGTGACCACGTGGACC" }, .b = { "GCTCTGTTACCAAACGCG", "AGAATCCTAAGG", "TTCAAGCGGAGTA", "TGTAAAAGCCCCGACGGATCCCCTCTCGCAGTTCTTATTCCCTTGTTAGAATTTGGCCCTCAGCCCAAACACTCATCGGTTTTGGATTTCAGACTGAATGAACCCCTGATTGTACACCTAGAGCCGGCGCTGGATTCTCCGTAAGCACACATAAAGTCAAACGACCCATGAGGTGAAAGACTTTCTAATGCTTCCGTAGGCGGACTTCCGTCTGAAGATAGGTATCGCAGGGCCGGCTTAGATCCATTTTGTTTCAGTAATCAATCGGCGCTCGCCCCCACTTATGTTTGCGTGAGAGAGTTTCTCGGGTTGCGGCGTCCATCAGCAACTCTGAGATGGCCACAAGGGTGAGTTTTAGATCCGCTCTGTGCCCGTCGTGTTCGGAGCCAACTGGAAGGTTTAAGGATCGCGATCTAGCATCTCAATTATCGGTCGAGATCGTATGGAGCGCGGCGCAGACATCCCTCAAGCTGTTATCGCCAGAATTGAAACCATACAACCCAGATGCTAAACAGGGGAGGATTAGAGCGTCTAGAGTCAGGTAACGGTGCCGCACGGGTTAGCTGTTGCACAGAAACTTCAGGAGATGTCTCTGAACGCGGCGCGTCCACTTCGCTCTTGGTTCATGTAATGGCACTTCGGGCCTTGGCGTAACCAGTCTTACACGTCGCTAGGTCGCGAGACCGGAGTTCAACTGCGTTTTGTCGGATCCTTGACAGGCGGAAGATATCTAACTATCTTAAATAGTCGTCACTGCAGCCCTCGGGCTAATGCGGTTCTTCATTTCCGTATATCGCTTACATACCGATGGTTAGTACAGTCTAGTGGGCGACGAGTGTGGGCACTCGTACAGCCCTTTAATGTAAGCCGCGGACACACCCGGCTCCGGTGAGCACTATGCAATGATTAAAGTCAA", "GTCGGGTTTATCAAACATTAGAGGGCCTGATACAATTGGTGCCGCCCGAGGGACCTGGACAGACAATCAGGTAACAGGGCATAAGCATGACGCCTCATC", "CGTATGCCGTAGTTCACAATGCACGCAAGTTCACGACCATGTTCGGTGTTTGGTTTGAGGGAGGATTCAGTTCACTATGTCGTGAACATTTCCGTGGTGGCGCTGGGGCCGGTGGCTGTCCAAAGGTAACGCGCTGCTTCACCTCGGACTGTGACGCGTGGCC" } },
		{ .a = { "AGGTCGCGGCTTAGGCGCGATGGC" }, .b = { "GGTCGCGTCTTAGGCCGATGGC" } },
		{ .a = { "CTCCCCAACCACAGAAGCTCTAAAGA" }, .b = { "ATCCTCAGACGAACATGTCGCCCTAGG" } },
		{ .a = { "GGTATTTGAATCGTGAACCCC", "TCCCGAGTTCGGCCGATGAGTGCGGGTTTAGGCTATTGTGTCCCGTCGC", "AAAGGTAGGTG", "CATCTTAGCTTGAGCGAGAATAGCTGTGCGGTGGACCAAGATCACGGCCGGATACGTCAATTTCCAGGAGAAGCCTATCCGTCCTGTGTATGTGTAGGGTAGTGAAGAATCGGGTCAGAAGTGCCGCAATTAGGTGTCGTATCGCATCGTCCGAGAGAATTGCCCTAGCTGAAACAATTTAAGAGGGGGAGCTTGATTTAGCCTGGTGCGCGGTGCCCTTCCAACTCAGAAACCTTATCCCTTAGGCATATAGACTTGTCGTCTCCTCAAGCTGCGAGGTTCTCCAAAGTCAGTATATAATTCTCTCCCTAAATCGACGGAATTTTGTGTTCGAATTCACAACCTTTTACGCCCAAAGTGCGTTTAGGCGTGGCTGAGCATTCCACCAACGCTGCGGTCTGACACAGCGAGGTTTACAACCTGCCTTTTAGCCACGAAATAGCTCATGAATGTGTACGTATCCACTGGTGCTAGGGGGAGGAACGTCTGAACCCCTATCTGGGGACAGTCACCCGTATGGTGTAATAAAAGGGAGTCAGGCTTGCGGAACTTAATCTTCATTGGTCCACAGGTCACTTCAACTACGCACGACCTTCAACCAAAACACCTGCTTAAAGCCTACACATGAAAGAAGACGCGGATCTTCACTCCTGCATCTTCTCCTTGAGTTATGCTGTTTGAACGGCCTACATGCGGTTATGAGACCGAAGGAGTTACATCCTCACTGAGGTAGGCTTTAACGAAGTATAAACTCCACTCTCAATGCTGATGGCGGTGGGGGTTGACGCATCGACCGGAAGTGACCGATCTCTGACGCAAGTCATGCGAGGGCCTGGCTCGGGGAGTGATGGTTGGCTACCTATATCGAATTGCCGTACTATCTTTCAGGACGTCCCTCCGAATCGGGTGGTTCGTGGGAAGCGGGTTCCATCAGTTCCGTTCGATCGGTTGGACATATACTATAGGTTTCGGTCGGGAAGTTAAAATATTGGTCGGTCCGTAAGCGATGCCGACGACCAAAGTATCGGTTAGTAAAGAACCCTCTCCAGGCGAAGACAGAATGTTGAGTATACGGATTCCCGTGGGGTGAAGTGGTCGCAGAATGGACGCGTTCCCTTAGAAATTGGAATAGTAGCTGGCTAAGGGTCCTTCAACTGTTTTTCGTTAGCGTATAGCGCGTCGCCGGAACTTAGCTCATAGCTGATCCCGACTGGAAAAACTAGCCCTCATTCGGGGTTTCTACCGC", "GGACGGTTTGATCAATGGCCGCGCCTTACAGTGGAAA", "TTCAGGATGGCTTTACCGCAGGAACTTCTTCTTGTCCTATCGGACTCGCGGTATTCGGTCGAACGTGTGTGTCCTGACCCTTAAGTCAGAGAAGCTGTGAGCATGTGGATTAGTACGCCAGTAGGCCACAATCATTACCAAAAGTATCGCGAAGGGAAGCATACTTATACAACATACGCTTTTAGCGGCCTCATGTTGTTTTAATTAACGTCCGCTATCCAATAACATATGTGGACCAGTGCAAACTAACCCCACCGCAGATCGCCATGGGATGATTATTCCAGCGATTACATTGCAGAGAGTAATTTTTGAATTCCTCATGATTAGCAAGGTGGTCCGATTAGGGTCTAACATTTTCCTCTCCCTAGAGTACTGTCGAAGCGGACGAGTTTGCATAACCTACAATACCAGTCAGGCAGTTCGCGTATTCAGGCTGATTGTCCCCTTGTCCGGATCTTAGCACTCACCTATGGCTTAGGATTGGTGTGCTCGTGCAGGCCTCTCGTAGCGGCGCTCTGAACCTATAACTCGGGCTAATTGGCTAGGCCACGCGCCCCCGAGAGCGCCGCAATGTTACAGCGAGACTGGAATTCCATTCCAGGTACGATGGAAGTCTGTGTTCGAGGTCTCCCTACAGGTCAACTTGGCCAACCGCAACCAGTCCTTCGCCTCCAAATTCATGCCATCCGCGCGTCCATGCGAGGAGCTAGTGTAGGCTGTAACTTGAACTCTTACTGTTGCGAAGTTTGCGTGTGCCGCGCTCATTACGTATCATTTGGGAACGATTCCCATACTTATAGGACTCCGATAGTCTCCGCGAGGTGAGCAGCTAGAATCGTCTGAACGCATTACGTTGCAGTCCATGAGTGGACACCCTGCGGCTAAGGGGACTGCTACCTATACTCTCACGGTACTATCGCGAAACTATCTTATAAATCACTGATTAAGCGTTGATATTTATGCGGCTGGTTCGTCCTCACTTAGTAGTTCTGTACCCATTTGCCGGTCTTGGGTAGTGCGTACGCGAGGTGCGTTACGCACAGTTTCGCGACCTATTGCTATCTGCAAGACGTACGGTTAAATCATCTGGGGCTGTAGTTTACGCTATGTTGATAAGCTATGTCCCAGTACGAGACCAAAAAGATGACTTGTCCATGACCAGACAAGGCTGCTACGTCTGACCGGCCAGGGGAGTAGGGTCAAACCGAGCGTTACCCGGTAAGAGCAATTAATAAAGGAAATTTTAAAGGAGGGGATTCTGCTCCAGTGGAGCTCAAGCGACCTTTTGGACACCGTCAGTGGGAGTAGCCAGGCCCCAATCCAGGGCTGATAGCAAACGGGTACACGCCGCGACCGTGGCGCCCGGCGTAACTAGTCCCTTGCTTGGGCATTCAATGGCCGGTCGCAAGCTGTACGTTAATTTAAATGGCAATTCCACCTTACGCGCAGCACCTGCTCTCGTCCGGTTAAGTCTCGCGGTAGCCCTGCACGTCAATGTGGAACCTTTCCCACTCGAAACCGAATAAGGACCACGCCCAGTCGAAATAATGCCGCAGGTAGACGCATGCTAAATCATAATTGCCAGCGAATGAAAGATAGGCGGATGTTTGTCTAGGCAGGTCGACAACCGTGCGGAACATTTCTTTCATTGTACTGCTTGGGATAGTCTTATGCGCATGAGATTCGCTCACTGGGTA" }, .b = { "GGTATTTAAAATCCGTGACTCC", "TCGAGTTCTGCCGATAGTGCGGGTTTACGCTAATGTGTCCCGGTCGC", "AAAGGTAGGTG", "TATTTAGCGTGAGCGAGAATAGCTGGCGGTGGACCAAGGTCACGGCGGCTACGTCACTATGCCAGGAGAACTTTCCTCCTGTGATGTGTACGCAGTGAAGGAATCGGGTCAGGAGTGCCGCAATTAGGTCTCGTATTCTCACGTCCGAGAGAGTGACCCTTGCTGAACAATTTAAGAGGGAGCTTGATTTGAGCCTGGTGGCGGTGCCCTTCGCACCTTCCCAGAAACCTTGATCCCTAGGCATTTAGACTTGTCGTCTCCTCAAGCATGCGAGGTTCTCCAAAGTCAGTAAATACTTCCTCCCAAGAACTCAACGGAATTTTGGTTCGATTCACACCTTTTACGCCCAAAGGTGCGATATGACGTAGGCTGAGCATTCTATCCAACGCTGCGGCTGACACAGTGAAGGTTTACCACCTGCCTTTTAGCCACGAAATAGATTCATGAATGAGCACGATCCACGTCGGTGCTAGCCGGGAGGAACGTCGGAACCCCTATCTGGGGAACAGTACAGCTGTATGTGTCATAAAAGGGGTCAGTGCTTAGCCAATGTTTAATCTTCATTGGTCCAAGGGTCACTTCACTACGCAGGACCTTCACCCAGATCCATCTGCTTAAAGCCTAATCATGAAAGAAAACGCGTGATCTGTCACTCCTGCATCTTCTCCTTGAGTTATGACTGTTGAACGCGCCTTCATGCGGTTATGAACCGAAGAAGTTACATACCTCACTGTAGGTAGGCTGTAAGGCGAAGTATAAACTCCACCCTCAATGCTGATGGCGGGTGGGGGTTGACGCTATCCCCGGAATTGACCGTTTCTGACGAGTTCATGCGAGGGCCCTGGTCGGTACGAGTGATGGGTTGGCTACTATATCGATTGCCGTACTTCTCTATCATGGCCGGCCTCCGTAATCGGGTGGTTCGTGGGAGCCGGGCCATCAGTTCCGCTCGATCCGGTTGGACATATATTTGCATTATCGTCGGGGAAGTTAAAATATTGGTCCGTCCGAAGCGTGCCGCGACGAAAGTCTCGGTTGTTAAGAACCCTCTCCAGGCGGAGGAGAATGTGCAGTACTTACGGATTCCCGTGGGTTGCAGTGTCGCAGAATCGACGCGTTCCCTATAGTATTGGAATAGTAGCTTGGCTAAGGGTCCTTCAACTGTTTTGCATTAGGGTGTAGTAGAGGCCTCGCGGACTTAGCTCATAGCTTATTCGACTGGAAGACAGCCTCATTTGGGGTTTCTACCCC", "GGGACGTTTGCGCAAGTCACGCGCGTACGTGGAA", "TTCAGGATGCTTTACTGCGGAACTTCCTTCTTGTCCTATCGGACTCGCGGTATTCGGGCGAACGTTGGATTCTGACCCTCTAAGTCAGAGAAGCTGTGAGCCAGCTGGTATTAGTACGCCAGTAGGACAAATCATTACCCAATAGACGAGAAAGGAAGCATACTTAACAATACACATTTTACGGCCTCTGTTGTTTTAATTAACGTCGCTATCTCAAAAACAATGTGGACCAGTGCAATCACTAACCCGACCGCAGACGCCAAGCGGATGATTATATGCAGCGATTACATTGGAGAGCGTAATTTTTGAATGCCTCATGATCTAGCGAGGTGGTCCGATTCGGTCTAACATTTTACTCTCCCTATGAGTCTGCCAAGGCCGCGTTTCGCATAACCTTACAATAGAATCAGGAGTTGGGCGTAATAAGCGATTGTCCCCTTGTCCGGGATCTTAAGCACTCACCCTCGCCTTAGGATTGATGTGGTCGTGCATGCCTCTCGAGCGGAGCCTGAACCTATAACTCGGGTAATTGGCTAGGCCACGCGCCCCCGAGAGCGCCGAAATGTTTCAGCGATACTGAATTCCCATTCCGGTGACGATGCGAGTCTGTGTTCGGGGTCTGCTACAGAAGTCAAATTGGCCAACCCAACCTAGACCTCGCCACCACATATCATGCCATCCAGCGCGACCATGGGAGGACTAGTGTAGGCTTTAACTTGCAACTTTACTTGTGCGAAGTTTCGTCTGCCGCGGCTCATTACGTACCATCTTGGGAACGATCTCCCATACTTATAGGTCCTCGATAGTCTCCGCTGAGGTGAGCAGCTAGAATCGTACTGCCAGCCTTACGTTGCAGTCACAATGAGTGGACACTCTGCGGTAAGGGACTGCCACCTATACATTCTTCAGGTACTATCGCGAAAATATCTAATAAATTCACATGATGAGCGTTATATTTATGCGGCTGGTTCCGTCACACTTAGTAGTTCTGTACCCATTTGCCGGTCGTGGGTAGTGCGTACGCGAGGTGCTTCGACAGTTTCGCGAACATCTATAGCTATCTGCAAGACGTACGGTCAAATCTATCTGGGGCGTGTAGTTTATGCTATCGGTGATAAGCTATTTTCAGTACGAGACCAAAAAGATTATTGTCCATGACCAGACAAGGTTGCTACAGTCCACCGGTCCAGGGGATAGGGTCAACGGACGTTACGCGTTAAGAGCAATTTAATTAAGAGAAATTTTAAAGAGGGGGATCTCTGCTCCAGTGGGGCTCAAACGGACCTTTGGACACCGTCAGTGGGAGTAGCCTAGGCCCCAATCCGGGTCTGAGTTGCCAAAAGGATACACGGTGCCGACCGTGGCGCCCGGCGTAACTAGTCCTTGCTAGGCTTTCAATGGCCGGTCGCACAAGCCTGCACGTTATAATTTTAATGGCACTCCACCTTTACCGCAGAACCTGCGTTCACGTTCCGAGTTAGTCTTGCAGGAGCCCTGCACGGCACTGTGGAACCTTTCCCACTGTAACCGAATAAGACCACGGCCCAGTCGAAATATGCGTGCAGTACACGCACTGACAAATCATAATTGCCTGCGAATCAAAGGATAGGCGTTGTTTGTCGGGCAGGTGGATAATCGTGCGGACCATTTTCTTTCATTGTACTAGCTTTGGGATAGTCTTACGCGAAAGTGAGATTCGCTCACTGAGTA" } },
		{ .a = { "GGGCATCCGCGAGTCTGACCAACTTGATTTGGCATCAA", "GG", "GA", "TATACAAGGTCATCTGTCAACCACAGTGAATAAGACGATATCAGGAACCACGGCAGTTAATTATCACAAAATCCCTCCCGGGCAACGTATTAGTTATAATCAGGCTGGCTCTCGTCTCCTAATTATTTTTTATCTCCCTGAACCGATCCACGAACATAATTTTGAGTCTTTTTACGAATGAAAACCCACCTAATGACGGCCCAAGCGGCGACGAGTCTAGATGAGCAATGCAACAGTTAACCGTCATTTTCTGATCCTCAGGAGTAAGTATACCATCTGAGTTATGCACGATAACACCGAGCCTGCGTAATGAGGTAAAGGACTCTTGACGGTTGCAACACAGCTCATCACTCTGGGCAGTGCGGTTGTATACCCAGGTTGCTTACCGTCAGACAGATAAATACCGCTAATGCGGCTTTGCAACCATGATGATGCAAGTAGGTAACCAACGCCTGGTCATCGTACAACGTGCGCCGGGGGTCAAAGTTATTCGGCTTCACAGGGTTCCGTACAAGAGCAGATTGCACACGGGAAAAACCGCCTCCCGAGGGTCGCGGTTACCGCCGCCATTCATGGAAGTGAGTGTGCGCCTGTCTCGCCCAGGTCAACGTAGCTCACGTGTGAGTATGCCTCTCCACCGTAGGCCTCTCACATGGGGCATGAAGCCGTGAGGGATATCTACAAAGGTAATTGTAACCTAATGGGTGGGTCTACTCAGAAGTAATACACAACTGACCGTACATCGGGGAGGGGGTCGACAACCTTCCCTACACAGTGACTTTCCGAGAATAAAGATTCATCCTGCCCAGGTCGAGACAACTCGTTGGGTTCAACTGGAGATCCCCTCGTCAGGGCCAAATTCTTGGCCTTCGATTGGAAGCAAAGGCGGCCCTCACTGGGCAAGTTACTATCCACTGGGAATCTACAATTGCTTTGGTTGTTTATAACGGCGCGCGTGGTCAAATTGGGCACATGTTGACACACTCTTTCGAGGACCGTTATAAGGACAATAGTGCACAATGACACGCTGATCCAATACAAAAGTGAAAAATTAATCGCCACAGGAGGATGCAGTTCTTCATATAAGCGTGTCTCTTTCAGTCGGAGTACAACATCACGGCCCGAGCTCTTAGAGACATCGTCTTAATCGGGTTATGGTTCGTTTGAAGTCATGTGGTGTGAGCTGCCAACCGACCGCCCAGCTAC", "GAAGGAAATGCCTCGATGC", "CGGGTCTGCAAACAGTAGTAAACTCTACCCGGAAATCACTAGTCGCCTTATTCATGTTTAACTAGTTTTCGATCCCGGACTTCATAGGTGTTTGACCTCCATCTCGCTCCCACGAAGGTCAGAAAACCAGGAGTATGCGTGTTTCTATCTAGCCGCCTGCGCGAAGTGTCCCCACGTCGCAGTCTCGCAAGGCAGGGCACTGTAGAGGCATCTGGGTATTAGGGGAACGAGCGCCCCTTGAGTTACCTCCGGTAACCAGTACATGCAAATGACGTCTCACGGTCTTGATCTTGTGAAGAACTCGCATACTCTGAACTTCGAAGGTAGTTTTAGATTCGTGGCGACGCCCGTCGCTCGGACCTCTAAGCTCTCGCCAAGCCAATGGGTAATCCGGAGTGTGATTGAAGAGTGGACAGTTGAGTCACGTAAGGACCTGCCGTCACCCCTCACCGAGTCAGGCAACGGAATGGTTTAAAGTCGAGCACACCGTGGAAAGAACGTACGCCGTAACAGGATGGCACTTCTTTAACGTCACGCTCGTTGGCTCACGAGGTACTTCTAAGGATACGGTGCGAGCCGGATAGGATCACAACGGCCTGACGGGAGACTTTTGTCAACGAGGATCTAGGGGTAACAGCCGAACGATTAACAGTACATCTCTCAAGCAATGGTGAGGTCGTGGCCTGGGTCAGGAGAGTCCCTATTGGGTACCTGTGACATGAGACACTTGAACAACTGTGGCACTTTAAATGACGCGTACCAGTCGAAATGGAGAACTACTTTACGCTCGATGGGTAGTACGTGATGGGGTGCTGGCTCGCGCAATTCCTCTCACCTTTGCTAGGGCTCGCTATTCGGGACTTGAAATTTCAGCGATTTCGCAGATTGCCGTCCCCCCAACGATTTATGTAAAATATGTAACTCTCGAGGCTACTAATTTGTTTTGTGCGGCAGGTTATCGACTGATGTGCCCCAATAGAGACGGTCCTCTGCACTAAGTGTATACTACCATCGGACTCCAATGAACAATTGTACTCATGCACGGTGTTAATATCCTTCTGCCCCGCGGAATACGGCAGAAGTCGCGCATCGAGACCCTTCTTGCCAGCGGGCCTTGCAGACCCACGCGGATCTCGAGAGAGACATAGACCAGGCGTATCGAATTGGAAGGACGGCCGTCTACGGTACTAGCATATCGTCCCCGTTTAGGCCTTTCGATCACGAAG" }, .b = { "GGAGCCTCGCGAGATCTGACCAAATTGATAAGGCATCAA", "GG", "GT", "CATACGAAGGGCATCTGTCAACCAAGTGATAAGACGATAACAGGACTCACGGCAATTTAATTATCACAAAATCCTTCCCGGGCACCGTATTAGTTAAGATCAGGCTGGCTCTCGTCTCCTAGTTAATTTTTTATCTCCCTGAAGCCGTCTCAACGAACATTATTTGTGAGCTATTACGAATGAAAACCCACCCTGATGCGGCCAGCCGCGACGAGTCTAGTTGAAGCATGCAACAGTTAACCTCTATTTTCTTATCTAAGGAGTAATTATACCATCTGAGTCTATGCGCGATAACAACGAGACACGCGATGAGAATAGGACCTGACGGTTGCAACACTAGCTCATCACTCTGGGCGTCGAGTTGTGTATACCCAGGTTGCTATCCGTCAGACAGATAAATACCCAATGCGGCTTTGCACTCCAGGATCATGCAGTAGGTAACCAACGCCTGGTAGCGTAGAAGTGTGTTGGGGGTCCAGTGATTCGGCTTCAAGGGTTCCGTACAAAGCTAGATTGCAGCACGGAAAAACCGCCCCCAAGGGTTCGCTGGTTAGCCGCGCCATTACTGGAGGTGGTGTGTCGCTCGTTCGCCGGGGTCCACGTAGCTCACTGTGAGTATGCCTCTCCAACGTAGCTCTCACATGGGGCATGACAGCCGTGAGGGAATATCTACAAAGGTAATGTAATCTAACGGGTGGGTGTAGTAGAGTAATCCCAATCTGACCGTACATCGGGGAGGGGGTCGCAACCTTCCTACACAGGACTTTCCGGACACAAAAGATTCTGGCTCCCAGGCGAGACAACGTCGTTGAGGTTCAACTGGATATCCCCTCGTCAGGGACCAACATTTCTTGACTTCGTTGGAAGCAAAGGGGCCTCACTGGGCAGTTCACTTGGCCTGGGAACTCTACAATTCTGTCTGCTTATTTATAACGGCGCGCGTGTAATATTGGGCACATGTGACACACTCTCTCGCAGACGTTATGAGGACAATACTGCACAATGACACGGAATACAAACCAAAAGTAAAAATTAATCGCCACAGTAGATGCTATTCTTCATAGTACGCGATGCTCGTTCAGTCGAGTACAAACCATCCCGGCTCGAGCCTTAGAGACATCGCCTGAATCGGGCTTATGGTTCGTTTGAGGTCATGTGGGGTGAGCTGCCAACCGACCGCCCTAGCTAC", "GAAGGAATGTTCCTCAAGTGC", "CGGTCCTGCAAACAGTGGTAAGCTCTACCCGGAAATCACTGAGTCGCCTTAAATTATGTTTAACTTAGCTTTCGATCCGGACTTCAAAGTGTTTGACCCCCATCCTCGCTCCCACGAAGGTCAGGAAACCAGAGATATCGGGTTACTATCTAGCCGCTGCGATAATTGCCCCCGGCTCACGTCTCTCAAGGCAGGGCATTGTAAAGGCATCTGGGTATTAAGGGGAACAGCGCCACCCTTGAGTAACCTCGTAACCAGTAAATGCAATGACGTCGTCATGGTCTTGATTCTTGGTGACGAACTCGCATAACTCTGACCTTCGAAGGCAGTTTTAGTTCGTGGCGGCGCCCGTGGTCGGACCTCTAAGCTCTCACCAAGCCAATGGGTAATCCGGAGTATGATTACGAGTGGACACGTTGGGTCACGTAAGGACCTTCTCGTCACCCCCCAACGAAGGCAGGCCACGAATGTTTAAAGTCGAGCATCTACCGTGGAAAGAACGTACGCCGTAACAGGATGGCACTTCTTTAACGTCACGCTCGTTGGCTCACGAGGTACTCACAGGATAGGGCGACCGGATAGGATCACATCGGCCTGACGGGAGACATTTGACAACGAGGATCAAACGGGTAAAACATCCGATCGGATTAACAGGATCATCTATCAGCATGTGGTGAGTCGTTGGCCGTGGGTCAGGAGAGTCCCTTTGGGTACGGTGACATGAGACACTTGAACTACTTGGCACCATAAATGACGTGAACCCAGTGAAATGGAGAACTAGTTTACGCTCGATCGGTGTAGGTGATCGGGTGCTGGCTCGCGCATGTTCCTCTCACCATTTCCTAGGGCTCGCTATCGGGACTGAAATTTCTGCGATTTCGCAGTATTCGCGGCCCCCCAACGCTTTATGTAACATATTTAACTCTCGAGGCCTACTAATTTGCTTTGTGCCGATGAACGACCGATGTGCTCTCAAATGGCACGGTCCTGCTGCAAGAAGTGTATACTACCATCGGACTCCAAGATGAACAATTGTGCACATGCACGGTGTTAATATCCTCAGCCACGCGGAAACCAGCAGCAGTCGCGCTCGACACCCTTGTTCCAGCGGGCCATGCGACCCAGGCGGATCTCGGAGAGAGTCATGGACCTAGGCGTTACGATTGGAAGGGCGGCCGCTACGCTGACTAGCATACCGTCCCCCGTTAGACTTTCGTCACGAAG" } },
		{ .a = { "ATACTCAATGAGCGCATCCGTCTGAAGCAATATACGC", "TAATGCCTGACTGAAGTAGCCGCCGCCATGGGTGCCTTCCTGC", "GCACGGGTCAATCCT", "ACTTTAATTCATTCCGTCTCTGAAGCATGTTGCGCATGTGTACTCGCATACTGGTATTTCTTAACCAATATTAAGCTTGGTTCGCCGAAGGTCGTGTATTAGGAGACCTTAGTCAAGTCTCTGGCGAAGTTGCCGACTCATCTGTTGCGTATGTCTCGATCTCTCTTGGATCGGTCCGGACTTAGCTTTTGATTTAGGGATCGCTCGCTACCTACCCAAAAGCTCATTTCAGTGTGTCTTCAACCCGGCTTCCAATATGTCGCAAAGACGTTGTATAGGCCCTTCTGCTACATGATTGTAAGTTCTGCGATTCGGGCAAACCAATATGTTTTCGGTGCAACTTTCGCTGCTTGATGATATGTTTGCGAGTAAGTTACATAGGTAGTGTGTGAGCTATGACCACTAATCCGCTCCTTAGGGCCGTGCACTATACAACCTGGGACTCAGCGCTTGCTAGTGGTCGACGCGGAAAACTCGCGAGCAAGTCAGTGCACGCTTACTACTGAAAGGAATACCTGGCCATCCCGCAGGATATCAACGGCTGGATCGCTTGGATTGATCTTATTCTTGGCTCTCCTTTAGGGAGGTAGGCCGCACGAAAGCATCCAAATAGCTGTCCTTACCGGCTCCTCTACGCTGCTCGGGCACGGCCGGGAAATCCGTCGATATGCTTTCGACTGTAACGGCTGAAAGTCCAACTGGGGGTAGCTCAACGGTTCGAGGGCTTCTACGTTACTATATAATCGGCGTTGCCATA", "ATACCTCCAACAACAAAGTGTCGGTTTTGTTTGAAATCCACTACAAGAGCTCGCTCAGATCAGAAGATACTCA", "TAACCGGAAACTCCGC" }, .b = { "ATACTCAAAGAAGCATCCTCTATGCAATATATGC", "TAATGCCTGACTGAAGTAGCCGCCGCCATGGGTGCCTTCCCGC", "GCAGGCAACCCT", "ACTTTAATTCATCCGCTCTCATGAACTGTTCGCGCATGTGTACTCGCATACTGGCTTTTCTTAACATATTAAGCTTGGTTCGCCAAAGGTCGTGTATGAGAGCTTTTAGGTCAATCTCCTGTCAAGTTGCCGACTCATCTGTTGCGTATGTCACGATCGCTCTTAGATAGTTCCCGGACTAACTTTTGATTTAGGATCGCTCGCTACCTCCCAAAACCTCATGTCAGAGTGTCCTTCAAATCCCGGCTTCCAATATGTCGCAAAGGGGTTGTATGGCCCTTCTGCGAATGATTGTAGTTCTGCGATTCGGGCAAACCAATATGTTTTCGGTGCAACTTTCCTTCTTGATGATAATTTTGCCGTTAAGTTACATAGGCTGGTCTGGTGAGCTAAGACCACTAAATCCGTCCACTTACGGCCGTGCACTGATACAACTGGCTCAGCGCTGTCTTAGTCGTCGACCGGAAATACTCGGCGAGCAGCAAAGCACTCTTACTACTCAAAGGAATACCTCGCATCCGCGGCTATCAACGGCTGGTCGTTGGATTGATACTTATTCTTGGCTCTCCTTTGGGGAGGTAGGCCGCACGAAACATCCCCGATAGCTGCCCTTTACAGGTCCTCTAACGCTGATGTGGGAGGCCGGGGAAATACGTGATATGCTTCACGCCCTTCACAGTGAAGAGTCCACGGGGGTAGCTCAACGTTCAGGGTTTCTACGTTAGTATATAATCGGAGGTTGCCATA", "ATACCCCAACAACAAAAATTCGGGCTTTGTTGAATCCACTACAAGTCGTTCAGACAGAATACTCA", "TAACCGGCAACTCGGC" } },
		{ .a = { "AAAAGCGGGA", "CCCAATATTCACTTGATTCGTCTTTA" }, .b = { "AAAAGCGG", "CGCAATAATTCACTTATTCGACTCTGA" } },
		{ .a = { "C", "TATAGCCTGTTGAAACAAACATGGATCTTGGGGCTCGAAAATTATTCTACCACTGCGTAC" }, .b = { "G", "TTCTAGCCTGTTGAAAGCAACAATGGATCTTGGGGCTCGAAATTATTACTACCATTGTGTAC" } },
		{ .a = { "C", "CAACCGAGCCCGTAATTAAA" }, .b = { "A", "TAGCCGAGCCCGAAATTAAA" } },
		{ .a = { "GACGTTT", "TCCTCCGGCTCGCATGAAATA" }, .b = { "GACGTA", "GCCTCCGGCGCGCATGAAATA" } },
		{ .a = { "AGGTCACCTGCACCCTTGC", "CCTATAGCAGCCATACTCGGTTTACATATAGAGGCCGTGGCTTGTTTGAGCTCAT" }, .b = { "AGGTCAACTGCACCCTTGA", "GAGTATAGAGCCATACTAGGTTTACATATAAGGCCGTTGGTCTTGTTTGAGCTCAT" } },
		{ .a = { "TATGCCT", "CTGTTATTGGTCACACTA" }, .b = { "TATGCA", "ATGTTATTGGGTCACACTA" } },
		{ .a = { "ATAGTCAACAACTGAGCAGAGCATATTATCCAACGTGT", "TTTAGCCAGCGTATTTAGGACTCC" }, .b = { "ATAGTCAATAATCTGACCGGAGTCATATTATATCCAACGTGG", "GAGGCGCCAGCGTATTTAGGAATTCC" } },
		{ .a = { "CTAAACTCTTGTTTT", "CCTCATCTCACCACCCACGATCGGAACAAAGCCCGGATCGT" }, .b = { "CTAAACTCTAGTAA", "CCTCATTCAACCACCCCGCGTCGAACAAGACCGGATCGT" } },

		/* fails for affine-16 due to the bandwidth shotage */
#if 0
		{
			.a = { "CCTTTACATATACGGACCAGAGATAGGTTAAATTTAAATCTTGCCGGT" },
			.b = { "CCCTAAGATACGGATTAGCTAACTAGCAACGGAGATACGTAAATGTTAAATCTTGCCGGT" }
		}
#endif
	};
	uint64_t j = 5; {
	// for(uint64_t j = 0; j < sizeof(unittest_default_params) / sizeof(struct gaba_params_s *); j++) {
		struct gaba_params_s const *p = unittest_default_params[j];
		struct gaba_context_s *g = _export(gaba_init)(p);
		struct gaba_dp_context_s *l = _export(gaba_dp_init)(g);

		for(uint64_t i = 0; i < sizeof(pairs) / sizeof(struct unittest_seq_pair_s); i++) {
			_export(gaba_dp_flush)(l);
			unittest_test_pair(UNITTEST_ARG_LIST, p, l, &pairs[i], 0);
			unittest_test_pair(UNITTEST_ARG_LIST, p, l, &pairs[i], 1);
		}

		_export(gaba_dp_clean)(l);
		_export(gaba_clean)(g);
	}
}

/**
 * @fn unittest_random_base
 */
static
char unittest_random_base(void)
{
	char const table[4] = {'A', 'C', 'G', 'T'};
	return(table[rand() % 4]);
}

/**
 * @fn unittest_generate_random_sequence
 */
static
char *unittest_generate_random_sequence(
	int64_t len)
{
	char *seq;		/** a pointer to sequence */
	seq = (char *)malloc(sizeof(char) * (len + UNITTEST_SEQ_MARGIN));

	if(seq == NULL) { return NULL; }
	for(int64_t i = 0; i < len; i++) {
		seq[i] = unittest_random_base();
	}
	seq[len] = '\0';
	return seq;
}

/**
 * @fn unittest_generate_mutated_sequence
 */
static
char *unittest_generate_mutated_sequence(
	char const *seq,
	double x,
	double d,
	int bw)
{
	if(seq == NULL) { return NULL; }

	int64_t wave = 0;			/** wave is q-coordinate of the alignment path */
	int64_t len = strlen(seq);
	char *mutated_seq = (char *)malloc(sizeof(char) * (2 * len + UNITTEST_SEQ_MARGIN));
	if(mutated_seq == NULL) { return NULL; }
	for(uint64_t i = 0, j = 0; j < len; i++, j++) {
		if(((double)rand() / (double)RAND_MAX) < x) {
			mutated_seq[i] = unittest_random_base();/** mismatch */
		} else if(((double)rand() / (double)RAND_MAX) < d) {
			if(rand() & 0x01 && wave > -bw+1) {
				mutated_seq[i] = seq[j++];			/** deletion */
				wave--;
			} else if(wave < bw-2) {
				mutated_seq[i] = unittest_random_base();
				j--; wave++;						/** insertion */
			} else {
				mutated_seq[i] = seq[j];
			}
		} else {
			mutated_seq[i] = seq[j];
		}
		mutated_seq[i + 1] = '\0';
	}
	return(mutated_seq);
}

unittest( .name = "cross" )
{
	uint64_t const cnt = 5000;

	struct unittest_context_s *c = (struct unittest_context_s *)gctx;
	struct gaba_stack_s const *stack = NULL;

	for(uint64_t i = 0; i < cnt; i++) {
		struct unittest_seq_pair_s pair = {
			.a = {
				unittest_generate_random_sequence((rand() % (_W + 10)) + 1),
				unittest_generate_random_sequence((rand() % (_W + 10)) + 1),
				unittest_generate_random_sequence((rand() % 16) + 1),
				unittest_generate_random_sequence((rand() % 2048) + 1),
				unittest_generate_random_sequence((rand() % 128) + 1),
				unittest_generate_random_sequence((rand() % 2048) + 1)
				// unittest_generate_random_sequence((rand() % 10) + TEST_LEN)
			}
		};

		for(uint64_t j = 0; pair.a[j] != NULL; j++) {
			// pair.b[j] = pair.a[j];
			pair.b[j] = unittest_generate_mutated_sequence(pair.a[j], 0.1, 0.1, _W);
		}

		/* issue flush / save-reload randomly */
		switch(rand() % 100) {
			case 0: _export(gaba_dp_flush)(c->dp);							/* fall through to save the head */
			case 1: stack = _export(gaba_dp_save_stack)(c->dp); break;		/* always overwrite */
			case 2: if(stack != NULL) { _export(gaba_dp_flush_stack)(c->dp, stack); stack = NULL; } break;
			default: break; /* do nothing */
		}

		unittest_test_pair(UNITTEST_ARG_LIST, c->params, c->dp, &pair, 0);
		unittest_test_pair(UNITTEST_ARG_LIST, c->params, c->dp, &pair, 1);

		for(uint64_t j = 0; pair.a[j] != NULL; j++) {
			free((void *)pair.a[j]);
			free((void *)pair.b[j]);
		}
	}
}

#endif /* UNITTEST */

/**
 * end of gaba.c
 */
