
/**
 * @file example.c
 */
#include "gaba.h"										/* just include gaba.h */
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

int printer(void *fp, uint64_t len, char c)
{
	return(fprintf((FILE *)fp, "%" PRIu64 "%c", len, c));
}

int main(int argc, char *argv[]) {
	/* create config */
	gaba_t *ctx = gaba_init(GABA_PARAMS(
		.xdrop = 100,
		GABA_SCORE_SIMPLE(2, 3, 5, 1)					/* match award, mismatch penalty, gap open penalty (G_i), and gap extension penalty (G_e) */
	));
	char const *a = "\x01\x08\x01\x08\x01\x08";			/* 4-bit encoded "ATATAT" */
	char const *b = "\x01\x08\x01\x02\x01\x08";			/* 4-bit encoded "ATACAT" */
	char const t[64] = { 0 };							/* tail array */

	struct gaba_section_s asec = gaba_build_section(0, a, strlen(a));
	struct gaba_section_s bsec = gaba_build_section(2, b, strlen(b));
	struct gaba_section_s tail = gaba_build_section(4, t, 64);

	/* create thread-local object */
	gaba_dp_t *dp = gaba_dp_init(ctx);					/* dp[0] holds a 64-cell-wide context */
	// gaba_dp_t *dp_32 = &dp[_dp_ctx_index(32)];			/* dp[1] and dp[2] are narrower ones */
	// gaba_dp_t *dp_16 = &dp[_dp_ctx_index(16)];

	/* init section pointers */
	struct gaba_section_s const *ap = &asec, *bp = &bsec;
	struct gaba_fill_s const *f = gaba_dp_fill_root(dp,	/* dp -> &dp[_dp_ctx_index(band_width)] makes the band width selectable */
		ap, 0,											/* a-side (reference side) sequence and start position */
		bp, 0,											/* b-side (query) */
		UINT32_MAX										/* max extension length */
	);

	/* until X-drop condition is detected */
	struct gaba_fill_s const *m = f;					/* track max */
	while((f->status & GABA_TERM) == 0) {
		if(f->status & GABA_UPDATE_A) { ap = &tail; }	/* substitute the pointer by the tail section's if it reached the end */
		if(f->status & GABA_UPDATE_B) { bp = &tail; }

		f = gaba_dp_fill(dp, f, ap, bp, UINT32_MAX);	/* extend the banded matrix */
		m = f->max > m->max ? f : m;					/* swap if maximum score was updated */
	}

	struct gaba_alignment_s *r = gaba_dp_trace(dp,
		m,												/* section with the max */
		NULL											/* custom allocator: see struct gaba_alloc_s in gaba.h */
	);

	printf("score(%" PRId64 "), path length(%" PRIu32 ")\n", r->score, r->plen);
	gaba_print_cigar_forward(
		printer, (void *)stdout,						/* printer */
		r->path,										/* bit-encoded path array */
		0,												/* offset is always zero */
		r->plen											/* path length */
	);
	printf("\n");

	/* clean up */
	gaba_dp_res_free(dp, r);
	gaba_dp_clean(dp);
	gaba_clean(ctx);
	return 0;
}

/**
 * end of example.c
 */
