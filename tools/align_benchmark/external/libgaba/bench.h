
/**
 * @file bench.h
 *
 * @brief benchmarking utils
 *
 * @detail
 * usage:
 * bench_t b;
 * bench_init(b);	// clear accumulator
 * 
 * bench_start(b);
 * // execution time between bench_start and bench_end is accumulated
 * bench_end(b);
 *
 * printf("%lld us\n", bench_get(b));	// in us
 */
#ifndef _BENCH_H_INCLUDED
#define _BENCH_H_INCLUDED

#include <stdint.h>
#include <string.h>

/**
 * benchmark macros
 */
#ifdef BENCH
#include <sys/time.h>

/**
 * @struct _bench
 * @brief benchmark variable container
 */
struct _bench {
	struct timeval s;		/** start */
	int64_t a;				/** accumulator */
};
typedef struct _bench bench_t;

/**
 * @macro bench_init
 */
#define bench_init(b) { \
	memset(&(b).s, 0, sizeof(struct timeval)); \
	(b).a = 0; \
}

/**
 * @macro bench_start
 */
#define bench_start(b) { \
	gettimeofday(&(b).s, NULL); \
}

/**
 * @macro bench_end
 */
#define bench_end(b) { \
	struct timeval _e; \
	gettimeofday(&_e, NULL); \
	(b).a += ( (_e.tv_sec  - (b).s.tv_sec ) * 1000000000 \
	         + (_e.tv_usec - (b).s.tv_usec) * 1000); \
}

/**
 * @macro bench_get
 */
#define bench_get(b) ( \
	(b).a \
)

#else /* #ifdef BENCH */

/** disable bench */
struct _bench { uint64_t _unused; };
typedef struct _bench bench_t;
#define bench_init(b) 		;
#define bench_start(b) 		;
#define bench_end(b)		;
#define bench_get(b)		( 0LL )

#endif /* #ifdef BENCH */

#endif /* #ifndef _BENCH_H_INCLUDED */
/**
 * end of bench.h
 */
