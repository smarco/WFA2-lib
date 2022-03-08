/*
The MIT License

Copyright (c) 2022-     Dana-Farber Cancer Institute

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

#include <stdint.h>
#include "lv89.h"

typedef struct {
	int32_t d, k;
} wf_diag_t;

static int32_t wf_step(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t n, wf_diag_t *a, int32_t is_global)
{
	int32_t j, m;
	wf_diag_t *b = a + n + 2; // temporary array

	// wfa_extend
#if 0 // unoptimized original version
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k, i = k + p->d;
		while (k + 1 < tl && i + 1 < ql && ts[k+1] == qs[i+1])
			++k, ++i;
		if (k + p->d == ql - 1 && (!is_global || k == tl - 1)) return -1;
		p->k = k;
	}
#else // optimized version learned from WFA
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k;
		int32_t max_k = (ql - p->d < tl? ql - p->d : tl) - 1;
		const char *ts_ = ts + 1, *qs_ = qs + p->d + 1;
		uint64_t cmp = 0;
		while (k + 7 < max_k) {
			uint64_t x = *(uint64_t*)(ts_ + k); // warning: unaligned memory access
			uint64_t y = *(uint64_t*)(qs_ + k);
			cmp = x ^ y;
			if (cmp == 0) k += 8;
			else break;
		}
		if (cmp)
			k += __builtin_ctzl(cmp) >> 3; // on x86, this is done via the BSR instruction: https://www.felixcloutier.com/x86/bsr
		else if (k + 7 >= max_k)
			while (k < max_k && *(ts_ + k) == *(qs_ + k)) // use this for generic CPUs. It is slightly faster than the unoptimized version
				++k;
		if (k + p->d == ql - 1 && (!is_global || k == tl - 1)) return -1;
		p->k = k;
	}
#endif

	// wfa_next
	b[0].d = a[0].d - 1;
	b[0].k = a[0].k + 1;
	b[1].d = a[0].d;
	b[1].k = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;
	for (j = 1; j < n - 1; ++j) {
		int32_t k = a[j-1].k;
		k = k > a[j+1].k + 1? k : a[j+1].k + 1;
		k = k > a[j].k + 1? k : a[j].k + 1;
		b[j+1].d = a[j].d, b[j+1].k = k;
	}
	if (n >= 2) {
		b[n].d = a[n-1].d;
		b[n].k = a[n-2].k > a[n-1].k + 1? a[n-2].k : a[n-1].k + 1;
	}
	b[n+1].d = a[n-1].d + 1;
	b[n+1].k = a[n-1].k;

	// drop out-of-bound cells
	for (j = 0, m = 0; j < n + 2; ++j)
		if (b[j].d + b[j].k < ql && b[j].k < tl)
			a[m++] = b[j];
	return m;
}

// mem should be at least (tl+ql)*16 long
int32_t lv_ed(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_global, uint8_t *mem)
{
	int32_t s = 0, n = 1;
	wf_diag_t *a = (wf_diag_t*)mem;
	a[0].d = 0, a[0].k = -1;
	while (1) {
		n = wf_step(tl, ts, ql, qs, n, a, is_global);
		if (n < 0) break;
		++s;
	}
	return s;
}
