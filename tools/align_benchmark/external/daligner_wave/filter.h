/*******************************************************************************************
 *
 *  Filter interface for the dazzler.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 ********************************************************************************************/

#ifndef _FILTER

#define _FILTER

#include "DB.h"
#include "align.h"

extern int    BIASED;
extern int    VERBOSE;
extern int    MINOVER;
extern int    HGAP_MIN;
extern int    SYMMETRIC;
extern int    IDENTITY;
extern char  *SORT_PATH;

extern uint64 MEM_LIMIT;
extern uint64 MEM_PHYSICAL;

int Set_Filter_Params(int kmer, int binshift, int suppress, int hitmin, int nthreads); 

void *Sort_Kmers(DAZZ_DB *block, int *len);

void Match_Filter(char *aname, DAZZ_DB *ablock, char *bname, DAZZ_DB *bblock,
                  void *atable, int alen, void *btable, int blen,
                  int comp, Align_Spec *asettings);

#endif
