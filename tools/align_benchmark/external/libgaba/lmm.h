
/**
 * @fn lmm.c
 *
 * @brief malloc with local context
 */
#ifndef _LMM_H_INCLUDED
#define _LMM_H_INCLUDED

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE		200112L
#endif

#ifndef _BSD_SOURCE
#define _BSD_SOURCE
#endif

// #define LMM_DEBUG

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* roundup */
#define _lmm_cutdown(x, base)		( (x) & ~((base) - 1) )
#define _lmm_roundup(x, base)		( ((x) + (base) - 1) & ~((base) - 1) )

/* constants */
#define LMM_ALIGN_SIZE				( 16 )

#define LMM_MIN_BASE_SIZE			( 128 )
#define LMM_DEFAULT_BASE_SIZE		( 1024 )

/* max and min */
#define LMM_MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define LMM_MIN2(x,y) 		( (x) > (y) ? (y) : (x) )


/**
 * @struct lmm_s
 */
struct lmm_s {
	uint8_t need_free;
	uint8_t pad1[7];
	uint32_t head_margin, tot_margin;
	void *ptr, *lim;
};
typedef struct lmm_s lmm_t;

/**
 * @fn lmm_init_margin
 */
static inline
lmm_t *lmm_init_margin(
	void *base,
	size_t base_size,
	uint32_t head_margin,
	uint32_t tail_margin)
{
	struct lmm_s *lmm = NULL;
	uint8_t need_free = 0;
	if(base == NULL || base_size < LMM_MIN_BASE_SIZE) {
		base_size = LMM_MAX2(base_size, LMM_DEFAULT_BASE_SIZE);
		base = malloc(base_size);
		need_free = 1;
	}

	lmm = (struct lmm_s *)base;
	lmm->need_free = need_free;
	lmm->ptr = (void *)((uintptr_t)base + sizeof(struct lmm_s));
	lmm->lim = (void *)((uintptr_t)base + _lmm_cutdown(base_size, LMM_ALIGN_SIZE));
	lmm->head_margin = head_margin;
	lmm->tot_margin = head_margin + tail_margin;
	return(lmm);
}

/**
 * @fn lmm_init
 */
static inline
lmm_t *lmm_init(
	void *base,
	size_t base_size)
{
	return(lmm_init_margin(base, base_size, 0, 0));
}


/**
 * @fn lmm_clean
 */
static inline
void *lmm_clean(
	lmm_t *lmm)
{
	if(lmm != NULL && lmm->need_free == 1) {
		free((void *)lmm);
		return(NULL);
	}
	return((void *)lmm);
}

/**
 * @fn lmm_reserve_mem
 */
static inline
void *lmm_reserve_mem(
	lmm_t *lmm,
	void *ptr,
	uint64_t size)
{
	uint64_t *sp = (uint64_t *)ptr;
	size = _lmm_roundup(size, LMM_ALIGN_SIZE);
	*sp = size;
	lmm->ptr = (void *)((uintptr_t)sp + LMM_ALIGN_SIZE + size);
	return((void *)((uintptr_t)sp + LMM_ALIGN_SIZE));
}

/**
 * @fn lmm_malloc
 */
static inline
void *lmm_malloc(
	lmm_t *lmm,
	size_t size)
{
	if(lmm == NULL) { return(malloc(size)); }

	size += lmm->tot_margin;
	#ifndef LMM_DEBUG
		if(((uintptr_t)lmm->ptr + LMM_ALIGN_SIZE + size) < (uintptr_t)lmm->lim) {
			return(lmm_reserve_mem(lmm, lmm->ptr, size) + lmm->head_margin);
		}
	#endif

	void *ptr = malloc(size);
	if(ptr == NULL) { return(NULL); }
	return(ptr + lmm->head_margin);
}

/**
 * @fn lmm_realloc
 */
static inline
void *lmm_realloc(
	lmm_t *lmm,
	void *ptr,
	size_t size)
{
	if(ptr == NULL) { lmm_malloc(lmm, size); }
	if(lmm == NULL) { return(realloc(ptr, size)); }

	size += lmm->tot_margin;
	ptr -= lmm->head_margin;
	#ifndef LMM_DEBUG
		/* check if prev mem (ptr) is inside mm */
		if((void *)lmm < ptr && ptr < lmm->lim) {
			uint64_t prev_size = *((uint64_t *)((uintptr_t)ptr - LMM_ALIGN_SIZE));
			if((uintptr_t)ptr + prev_size == (uintptr_t)lmm->ptr	/* the last block */
			&& (uintptr_t)ptr + size < (uintptr_t)lmm->lim) {		/* and room for expansion */
				return(lmm_reserve_mem(lmm,
					(void *)((uintptr_t)ptr - LMM_ALIGN_SIZE), size
				) + lmm->head_margin);
			}

			/* no room realloc to outside */
			void *np = malloc(size);
			if(np == NULL) { return(NULL); }
			np += lmm->head_margin;
			memcpy(np, ptr, prev_size);
			lmm->ptr = (void *)((uintptr_t)ptr - LMM_ALIGN_SIZE);
			return(np);
		}
	#endif

	/* pass to library realloc */
	void *p = realloc(ptr, size);
	if(p == NULL) { p = ptr; }
	return(p + lmm->head_margin);
}

/**
 * @fn lmm_free
 */
static inline
void lmm_free(
	lmm_t *lmm,
	void *ptr)
{
	if(ptr == NULL) { return; }
	if(lmm == NULL) { free(ptr); }

	ptr -= lmm->head_margin;
	#ifndef LMM_DEBUG
		if((void *)lmm < ptr && ptr < lmm->lim) {
			/* no need to free */
			uint64_t prev_size = *((uint64_t *)((uintptr_t)ptr - LMM_ALIGN_SIZE));
			if((uintptr_t)ptr + prev_size == (uintptr_t)lmm->ptr) {
				lmm->ptr = (void *)((uintptr_t)ptr - LMM_ALIGN_SIZE);
			}
			return;
		}
	#endif
	free(ptr);
	return;
}

/**
 * @fn lmm_strdup
 */
static inline
char *lmm_strdup(
	lmm_t *lmm,
	char const *str)
{
	int64_t len = strlen(str);
	char *s = (char *)lmm_malloc(lmm, len + 1);
	memcpy(s, str, len + 1);
	return(s);
}


/* object pool implementation */

/**
 * @struct lmm_pool_block_s
 */
struct lmm_pool_block_s {
	struct lmm_pool_block_s *next;
	struct lmm_pool_block_s *prev;
	int64_t cnt;
	int64_t pad;
	uint8_t mem[];
};

/**
 * @struct lmm_pool_object_s
 */
struct lmm_pool_object_s {
	struct lmm_pool_object_s *next;
};

/**
 * @struct lmm_pool_s
 */
struct lmm_pool_s {
	lmm_t *lmm;
	struct lmm_pool_block_s *curr;
	struct lmm_pool_block_s *root;
	struct lmm_pool_object_s *tail;
	int64_t rem;
	int64_t object_multiplier;
};
typedef struct lmm_pool_s lmm_pool_t;

/**
 * @fn lmm_pool_init
 */
static inline
lmm_pool_t *lmm_pool_init(
	lmm_t *lmm,
	uint64_t object_size,
	uint64_t init_object_cnt)
{
	object_size = _lmm_roundup(object_size, sizeof(struct lmm_pool_object_s));
	struct lmm_pool_s *pool = (struct lmm_pool_s *)lmm_malloc(lmm,
		  sizeof(struct lmm_pool_s) + sizeof(struct lmm_pool_block_s)
		+ init_object_cnt * object_size);
	if(pool == NULL) {
		return(NULL);
	}

	struct lmm_pool_block_s *blk = (struct lmm_pool_block_s *)(pool + 1);

	/* init pool object */
	pool->lmm = lmm;
	pool->curr = blk;
	pool->root = blk;
	pool->tail = (struct lmm_pool_object_s *)blk->mem;
	pool->rem = init_object_cnt;
	pool->object_multiplier = object_size / sizeof(struct lmm_pool_object_s);

	/* init first block */
	blk->next = NULL;
	blk->prev = NULL;
	blk->cnt = init_object_cnt;

	/* init free list */
	pool->tail->next = NULL;

	return((lmm_pool_t *)pool);
}

/**
 * @fn lmm_pool_clean
 */
static inline
void lmm_pool_clean(
	lmm_pool_t *_pool)
{
	struct lmm_pool_s *pool = (struct lmm_pool_s *)_pool;
	struct lmm_pool_block_s *blk = pool->root->next;
	lmm_t *lmm = pool->lmm;

	while(blk != NULL) {
		struct lmm_pool_block_s *next_blk = blk->next;
		lmm_free(lmm, blk); blk = next_blk;
	}
	lmm_free(lmm, pool);
	return;
}

/**
 * @fn lmm_pool_flush
 */
static inline
void lmm_pool_flush(
	lmm_pool_t *_pool)
{
	struct lmm_pool_s *pool = (struct lmm_pool_s *)_pool;
	if(pool == NULL) { return; }

	pool->curr = pool->root;
	pool->tail = (struct lmm_pool_object_s *)pool->root->mem;
	pool->rem = pool->root->cnt;
	pool->tail->next = NULL;
	return;
}

/**
 * @fn lmm_pool_add_block
 */
static inline
void lmm_pool_add_block(
	struct lmm_pool_s *pool)
{
	if(pool->curr->next == NULL) {
		int64_t next_cnt = 2 * pool->curr->cnt;
		struct lmm_pool_block_s *blk = pool->curr->next =
			(struct lmm_pool_block_s *)lmm_malloc(pool->lmm,
				  sizeof(struct lmm_pool_block_s)
				+ next_cnt * sizeof(struct lmm_pool_object_s) * pool->object_multiplier);
		blk->next = NULL;
		blk->prev = pool->curr;
		blk->cnt = next_cnt;
	}

	pool->curr = pool->curr->next;
	pool->tail = (struct lmm_pool_object_s *)pool->curr->mem;
	pool->rem = pool->curr->cnt;
	return;
}

/**
 * @fn lmm_pool_create_object
 */
static inline
void *lmm_pool_create_object(
	lmm_pool_t *_pool)
{
#ifdef LMM_POOL_SEPARATE_NODE
	struct lmm_pool_s *pool = (struct lmm_pool_s *)_pool;
	return(malloc(sizeof(struct lmm_pool_object_s) * pool->object_multiplier));
#else

	struct lmm_pool_s *pool = (struct lmm_pool_s *)_pool;
	struct lmm_pool_object_s *obj = NULL;

	if(pool->tail->next != NULL) {
		obj = pool->tail->next;
		pool->tail->next = obj->next;
	} else {
		obj = pool->tail;
		pool->tail += pool->object_multiplier;

		if(--pool->rem <= 0) {
			lmm_pool_add_block(pool);
		}

		*pool->tail = *obj;
	}
	return((void *)obj);

#endif
}

/**
 * @fn lmm_pool_delete_object
 */
static inline
void lmm_pool_delete_object(
	lmm_pool_t *_pool,
	void *_obj)
{
#ifdef LMM_POOL_SEPARATE_NODE

	free(_obj);

#else

	struct lmm_pool_s *pool = (struct lmm_pool_s *)_pool;
	struct lmm_pool_object_s *obj = (struct lmm_pool_object_s *)_obj;

	obj->next = pool->tail->next;
	pool->tail->next = obj;
	return;

#endif
}

/**
 * kvec.h from https://github.com/ocxtal/kvec.h
 * the original implementation of kvec.h is found at
 * https://github.com/attractivechaos/klib
 *
 * ===== start kvec.h =====
 */

/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
  An example:

#include "kvec.h"
int main() {
	kvec_t(int) array;
	kv_init(array);
	kv_push(int, array, 10); // append
	kv_a(int, array, 20) = 5; // dynamic
	kv_A(array, 20) = 4; // static
	kv_destroy(array);
	return 0;
}
*/

/*
  
  2016-0410
    * add kv_pushm

  2016-0401

    * modify kv_init to return object
    * add kv_pusha to append arbitrary type element
    * add kv_roundup
    * change init size to 256

  2015-0307

    * add packed vector. (Hajime Suzuki)

  2008-09-22 (0.1.0):

	* The initial version.

*/
// #define kv_roundup(x) (-, 32-(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#define lmm_kv_roundup(x, base)			( (((x) + (base) - 1) / (base)) * (base) )


#define LMM_KVEC_INIT_SIZE 			( 256 )

/**
 * basic vectors (kv_*)
 */
#define lmm_kvec_t(type)		struct { uint64_t n, m; type *a; }
#define lmm_kv_init(lmm, v)		({ (v).n = 0; (v).m = LMM_KVEC_INIT_SIZE; (v).a = (type *)lmm_malloc((lmm), (v).m * sizeof(*(v).a)); (v); })
#define lmm_kv_destroy(lmm, v)	{ lmm_free((void *)(lmm), (v).a); (v).a = NULL; }
// #define lmm_kv_A(v, i)      ( (v).a[(i)] )
#define lmm_kv_pop(lmm, v)		( (v).a[--(v).n] )
#define lmm_kv_size(v)			( (v).n )
#define lmm_kv_max(v)			( (v).m )

#define lmm_kv_clear(lmm, v)	( (v).n = 0 )
#define lmm_kv_resize(lmm, v, s) ({ \
	uint64_t _size = LMM_MAX2(LMM_KVEC_INIT_SIZE, (s)); \
	(v).m = _size; \
	(v).n = LMM_MIN2((v).n, _size); \
	(v).a = (type *)lmm_realloc((lmm), (v).a, sizeof(*(v).a) * (v).m); \
})

#define lmm_kv_reserve(lmm, v, s) ( \
	(v).m > (s) ? 0 : ((v).m = (s), (v).a = (type *)lmm_realloc((lmm), (v).a, sizeof(*(v).a) * (v).m), 0) )

#define lmm_kv_copy(lmm, v1, v0) do {								\
		if ((v1).m < (v0).n) lmm_kv_resize(lmm, v1, (v0).n);			\
		(v1).n = (v0).n;									\
		memcpy((v1).a, (v0).a, sizeof(*(v).a) * (v0).n);	\
	} while (0)

#define lmm_kv_push(lmm, v, x) do {									\
		if ((v).n == (v).m) {								\
			(v).m = (v).m * 2;								\
			(v).a = (type *)lmm_realloc((lmm), (v).a, sizeof(*(v).a) * (v).m);	\
		}													\
		(v).a[(v).n++] = (x);								\
	} while (0)

#define lmm_kv_pushp(lmm, v) ( \
	((v).n == (v).m) ?									\
	((v).m = (v).m * 2,									\
	 (v).a = (type *)lmm_realloc((lmm), (v).a, sizeof(*(v).a) * (v).m), 0)	\
	: 0), ( (v).a + ((v).n++) )

/* lmm_kv_pusha will not check the alignment of elem_t */
#define lmm_kv_pusha(lmm, elem_t, v, x) do { \
		uint64_t size = lmm_kv_roundup(sizeof(elem_t), sizeof(*(v).a)); \
		if(sizeof(*(v).a) * ((v).m - (v).n) < size) { \
			(v).m = LMM_MAX2((v).m * 2, (v).n + (size));				\
			(v).a = (type *)lmm_realloc((lmm), (v).a, sizeof(*(v).a) * (v).m);	\
		} \
		*((elem_t *)&((v).a[(v).n])) = (x); \
		(v).n += size / sizeof(*(v).a); \
	} while(0)

#define lmm_kv_pushm(lmm, v, arr, size) do { \
		if(sizeof(*(v).a) * ((v).m - (v).n) < (uint64_t)(size)) { \
			(v).m = LMM_MAX2((v).m * 2, (v).n + (size));				\
			(v).a = (type *)lmm_realloc((lmm), (v).a, sizeof(*(v).a) * (v).m);	\
		} \
		for(uint64_t _i = 0; _i < (uint64_t)(size); _i++) { \
			(v).a[(v).n + _i] = (arr)[_i]; \
		} \
		(v).n += (uint64_t)(size); \
	} while(0)

#define lmm_kv_a(lmm, v, i) ( \
	((v).m <= (size_t)(i) ? \
	((v).m = (v).n = (i) + 1, lmm_kv_roundup((v).m, 32), \
	 (v).a = (type *)lmm_realloc((lmm), (v).a, sizeof(*(v).a) * (v).m), 0) \
	: (v).n <= (size_t)(i) ? (v).n = (i) + 1 \
	: 0), (v).a[(i)])

/** bound-unchecked accessor */
#define lmm_kv_at(v, i)			( (v).a[(i)] )
#define lmm_kv_ptr(v)			( (v).a )

/** heap queue : elements in v must be orderd in heap */
#define lmm_kv_hq_init(lmm, v)		{ lmm_kv_init(lmm, v); (v).n = 1; }
#define lmm_kv_hq_destroy(lmm, v)	lmm_kv_destroy(lmm, v)
#define lmm_kv_hq_size(v)			( lmm_kv_size(v) - 1 )
#define lmm_kv_hq_max(v)			( lmm_kv_max(v) - 1 )
#define lmm_kv_hq_clear(lmm, v)		( (v).n = 1 )

#define lmm_kv_hq_resize(lmm, v, s)		( lmm_kv_resize(lmm, v, (s) + 1) )
#define lmm_kv_hq_reserve(lmm, v, s)	( lmm_kv_reserve(lmm, v, (s) + 1) )

#define lmm_kv_hq_copy(lmm, v1, v0)		lmm_kv_copy(lmm, v1, v0)

#define lmm_kv_hq_n(v, i) ( *((int64_t *)&v.a[i]) )
#define lmm_kv_hq_push(lmm, v, x) { \
	lmm_kv_push(lmm, v, x); \
	uint64_t i = (v).n - 1; \
	while(i > 1 && (lmm_kv_hq_n(v, i>>1) > lmm_kv_hq_n(v, i))) { \
		(v).a[0] = (v).a[i>>1]; \
		(v).a[i>>1] = (v).a[i]; \
		(v).a[i] = (v).a[0]; \
		i >>= 1; \
	} \
}
#define lmm_kv_hq_pop(lmm, v) ({ \
	uint64_t i = 1, j = 2; \
	(v).a[0] = (v).a[i]; \
	(v).a[i] = (v).a[--(v).n]; \
	(v).a[(v).n] = (v).a[0]; \
	while(j < (v).n) { \
		uint64_t k; \
		k = (j + 1 < (v).n && lmm_kv_hq_n(v, j + 1) < lmm_kv_hq_n(v, j)) ? (j + 1) : j; \
		k = (lmm_kv_hq_n(v, k) < lmm_kv_hq_n(v, i)) ? k : 0; \
		if(k == 0) { break; } \
		(v).a[0] = (v).a[k]; \
		(v).a[k] = (v).a[i]; \
		(v).a[i] = (v).a[0]; \
		i = k; j = k<<1; \
	} \
	v.a[v.n]; \
})

/**
 * 2-bit packed vectors (lmm_kpv_*)
 * v.m must be multiple of lmm_kpv_elems(v).
 */
#define _BITS						( 2 )

/**
 * `sizeof(*((v).a)) * 8 / _BITS' is a number of packed elements in an array element.
 */
#define lmm_kpv_elems(v)			( sizeof(*((v).a)) * 8 / _BITS )
#define lmm_kpv_base(v, i)			( ((i) % lmm_kpv_elems(v)) * _BITS )
#define lmm_kpv_mask(v, e)			( (e) & ((1<<_BITS) - 1) )

#define lmm_kpvec_t(type)			struct { uint64_t n, m; type *a; }
#define lmm_kpv_init(lmm, v)		( (v).n = 0, (v).m = LMM_KVEC_INIT_SIZE * lmm_kpv_elems(v), (v).a = (type *)lmm_malloc(lmm, (v).m * sizeof(*(v).a)) )
#define lmm_kpv_destroy(lmm, v)		{ free((void *)lmm, (v).a); (v).a = NULL; }

// #define lmm_kpv_A(v, i) ( lmm_kpv_mask(v, (v).a[(i) / lmm_kpv_elems(v)]>>lmm_kpv_base(v, i)) )
#define lmm_kpv_pop(lmm, v) 		( (v).n--, lmm_kpv_mask(v, (v).a[(v).n / lmm_kpv_elems(v)]>>lmm_kpv_base(v, (v).n)) )
#define lmm_kpv_size(v)				( (v).n )
#define lmm_kpv_max(v)				( (v).m )

/**
 * the length of the array is ((v).m + lmm_kpv_elems(v) - 1) / lmm_kpv_elems(v)
 */
#define lmm_kpv_amax(v)				( ((v).m + lmm_kpv_elems(v) - 1) / lmm_kpv_elems(v) )
#define lmm_kpv_asize(v)			( ((v).n + lmm_kpv_elems(v) - 1) / lmm_kpv_elems(v) )

#define lmm_kpv_clear(lmm, v)		( (v).n = 0 )
#define lmm_kpv_resize(lmm, v, s) ({ \
	uint64_t _size = LMM_MAX2(LMM_KVEC_INIT_SIZE, (s)); \
	(v).m = _size; \
	(v).n = LMM_MIN2((v).n, _size); \
	(v).a = (type *)lmm_realloc((lmm), (v).a, sizeof(*(v).a) * lmm_kpv_amax(v)); \
})

#define lmm_kpv_reserve(lmm, v, s) ( \
	(v).m > (s) ? 0 : ((v).m = (s), (v).a = (type *)lmm_realloc(lmm, (v).a, sizeof(*(v).a) * lmm_kpv_amax(v)), 0) )

#define lmm_kpv_copy(lmm, v1, v0) do {							\
		if ((v1).m < (v0).n) lmm_kpv_resize(lmm, v1, (v0).n);	\
		(v1).n = (v0).n;								\
		memcpy((v1).a, (v0).a, lmm_kpv_amax(v));				\
	} while (0)
#define lmm_kpv_push(lmm, v, x) do {								\
		if ((v).n == (v).m) {							\
			(v).a = (type *)lmm_realloc(lmm, (v).a, 2 * sizeof(*(v).a) * lmm_kpv_amax(v)); \
			memset((v).a + lmm_kpv_amax(v), 0, lmm_kpv_amax(v));	\
			(v).m = (v).m * 2;							\
		}												\
		(v).a[(v).n / lmm_kpv_elems(v)] = 					\
			  ((v).a[(v).n / lmm_kpv_elems(v)] & ~(lmm_kpv_mask(v, 0xff)<<lmm_kpv_base(v, (v).n))) \
			| (lmm_kpv_mask(v, x)<<lmm_kpv_base(v, (v).n)); \
		(v).n++; \
	} while (0)

#define lmm_kpv_a(lmm, v, i) ( \
	((v).m <= (size_t)(i)? \
	((v).m = (v).n = (i) + 1, lmm_kv_roundup((v).m, 32), \
	 (v).a = (type *)lmm_realloc(lmm, (v).a, sizeof(*(v).a) * (v).m), 0) \
	: (v).n <= (size_t)(i)? (v).n = (i) + 1 \
	: 0), lmm_kpv_mask(v, (v).a[(i) / lmm_kpv_elems(v)]>>lmm_kpv_base(v, i)) )

/** bound-unchecked accessor */
#define lmm_kpv_at(v, i)		( lmm_kpv_mask(v, (v).a[(i) / lmm_kpv_elems(v)]>>lmm_kpv_base(v, i)) )
#define lmm_kpv_ptr(v)			( (v).a )

/**
 * ===== end of kvec.h =====
 */

#endif /* _LMM_H_INCLUDED */
/**
 * end of lmm.c
 */
