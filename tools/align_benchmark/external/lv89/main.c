#include <assert.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "lv89.h"
#include "edlib.h"
#ifdef _USE_WFA
#include "gap_affine/affine_wavefront_align.h"
#endif
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c, s, is_global = 0, use_edlib = 0, use_wfa = 0;

	while ((c = ketopt(&o, argc, argv, 1, "glw", 0)) >= 0) {
		if (c == 'g') is_global = 1;
		else if (c == 'l') use_edlib = 1;
		else if (c == 'w') use_wfa = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: ed-test [options] <in1.fa> <in2.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -g    count gaps at the end of the target\n");
		fprintf(stderr, "  -l    use edlib instead\n");
#ifdef _USE_WFA
		fprintf(stderr, "  -w    use WFA\n");
#endif
		return 1;
	}
	assert(!use_edlib || !use_wfa);

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);
	kseq_read(ks1);
	kseq_read(ks2);

	if (use_edlib) {
		EdlibAlignResult rst;
		rst = edlibAlign(ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l,
				edlibNewAlignConfig(-1, is_global? EDLIB_MODE_NW : EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0));
		s = rst.editDistance;
#ifdef _USE_WFA
	} else if (use_wfa) {
		mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
		affine_penalties_t pan = { .match = 0, .mismatch = 1, .gap_opening = 1, .gap_extension = 1 }; // Init Affine-WFA
		affine_wavefronts_t *wf =
			affine_wavefronts_new_complete(ks2->seq.l, ks1->seq.l, &pan, NULL, mm_allocator);
		affine_wavefronts_align(wf, ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l);
		s = edit_cigar_score_gap_affine(&wf->edit_cigar, &pan);
#endif
	} else {
		uint8_t *mem = (uint8_t*)malloc((ks1->seq.l + ks2->seq.l) * 16);
		s = lv_ed(ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, is_global, mem);
		free(mem);
	}
	printf("%s\t%s\t%d\n", ks1->name.s, ks2->name.s, s);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
