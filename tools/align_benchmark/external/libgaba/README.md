# libgaba

Libgaba is a fast extension alignment (semi-global alignment) library for nucleotide sequences. It uses difference DP to calculate arbitrary large score within 8-bit-wide integers. Two heuristic algorithms, adaptive banded DP and X-drop termination, are incorporated for practical use in nucleotide sequence alignment tools.

## Usage

### Build

Typing `make` will build an all-in-one archive `libgaba.a`. The library is intended to be placed and built inside the source tree of the tool that depends on it. Place the libgaba directory (as git submodule, for example) in anywhere you like and invoke make (with proper `CC`, `CFLAGS`, and `ARCHFLAGS`).


### Example Source

Include a header: `gaba.h`. `make example` will compile the following example code and link it with `libgaba.a`.

```
#include "gaba.h"										/* just include gaba.h */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* typedef int (*gaba_printer_t)(void *, uint64_t, char) */
int printer(void *fp, uint64_t len, char c)
{
	return(fprintf((FILE *)fp, "%ld%c", len, c));
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

	printf("score(%ld), path length(%lu)\n", r->score, r->plen);
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
```

### Input sequence formats

Input sequences are provided as arrays of `uint8_t` (no need to be null-terminated). Sequence encoding must be either of the following two, 2-bit encoding (A = 0x00, C = 0x01, G = 0x02, T = 0x03) or 4-bit (A = 0x01, C = 0x02, G = 0x04, T = 0x08). The sequence encoding is configured at compile time and the default encoding is 4-bit. The configuration can be changed by overwriting BIT flag; `make BIT=2` will build all-in-one binary with the 2-bit encoding setting. Ambiguous bases can be represented in the 4-bit encoding by OR-ing the base alphabets. The scoring criteria for the ambiguous bases are "match if at least one base is shared between two input letters, otherwise mismatch."

### Substitution matrix

Match awards and mismatch penalties are always represented as a 4x4-matrix. The diagonal elements (at [0][0], [1][1], ..., [3][3]) are match awards (in positive 8bit integers) for ('A', 'A'), ('C', 'C'), ..., ('T', 'T'), and others are mismatch penalties (in negative integers). In 4-bit-encoding mode (the default configuration), the match awards and mismatch penalties are respectively set uniform for all the pairs, by averaging the four matching and twelve mismatching pairs.

### Gap penalty functions

The library supports the following three gap penalty functions:

* linear:   g(k) =          ge * k          where              ge > 0
* affine:   g(k) =     gi + ge * k          where gi > 0,      ge > 0
* combined: g(k) = min(gi + ge * k, gf * k) where gi > 0, gf > ge > 0

The init function `gaba_init` determines which score model fits the most to the given parameter set and construct a dispatcher for fill-in and traceback functions: `gaba_dp_fill_root`, `gaba_dp_fill`, and `gaba_dp_trace`. The functions will internally call the actual function from the dispatcher table, which is optimized for the score parameter. Setting gi = gfa = gfb = 0 will fetch the linear-gap penalty model, gfa = gfb = 0 will fetch the affine-gap, otherwise the combined gap penalty model is selected.

### Input sequences

Libgaba interprets each of input sequences as "concatenation of subsequences," that is, alignment can be incrementally extended by multiple `gaba_dp_fill` function call. Each call of fill function (`gaba_dp_fill_root` or `gaba_dp_fill`) will forward the front vector as far as possible but always leaves a triangular region at the front that cannot be covered without projecting the vector outside of the current matrix (area marked by * in the diagram below).

(Note: The library uses anti-diagonally placed constant-wide vector to fill a banded matrix. The vector is called "forefront vector.")

```
      section A   ...
     +---------+------
     |    \    |
     |\    \   |
     | \    \  |
  B  |  \    \ |
     |   \    \| margin
     |    \   /|<-->
     |     \ /*|
     +---------+------
     |         |
```


Sometimes the forefront vector will not reach the end of the sequences. It is more likely to occur when the input sequence lengths greatly differ. Whether the vector reached the ends or not is determined by examining `tail->status & GABA_UPDATE_A` and `tail->status & GABA_UPDATE_B`. Alignment can be extended unlimitedly (until it detects the X-drop terminate condition) by iteration of call and swap pairs (see `while((f->status & GABA_TERM) == 0)` loop in the example code).

### Sections

Input subsequences are distinguished by their ids. It will be any 32-bit integer but 0xffffffff and 0xfffffffe are reserved for internal use. The `base` and `len` are a pair of pointer to a sequence and its length.

```
struct gaba_section_s {
	uint32_t id;
	uint32_t len;
	uint8_t const *base;
};
```

#### Reverse-complemented sequences

If you need to align the reverse-complement of a sequence you have on memory, you do not need to explicitly construct the reverse-complemented sequence before calling the fill-in functions. Passing a pointer mirrored at GABA_EOU (end-of-userland) will fetch the sequence at the original location in the reverse-complemented manner. `gaba_mirror(sequence, strlen(sequence))` will create a pointer to the head of the reverse-complemented sequence (or the **tail of the original sequence**).


## References and Links

* Hajime Suzuki and Masahiro Kasahara (2018), Introducing difference recurrence relations for faster semi-global alignment of long sequences. [https://doi.org/10.1186/s12859-018-2014-8](https://doi.org/10.1186/s12859-018-2014-8)

* Benchmark scripts and alternative implementations [https://github.com/ocxtal/diffbench](https://github.com/ocxtal/diffbench)
* Benchmark scripts on the adaptive banding algorithm [https://github.com/ocxtal/adaptivebandbench](https://github.com/ocxtal/adaptivebandbench)
* Older version of the library [libsea: https://bitbucket.org/suzukihajime/libsea/](https://bitbucket.org/suzukihajime/libsea/).


## License

Apache v2.

Copyright (c) 2015-2016 Hajime Suzuki