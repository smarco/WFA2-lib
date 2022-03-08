/*******************************************************************************************
 *
 *  Fast local alignment filter for long, noisy reads based on "dumbing down" of my RECOMB 2005
 *     filter with Jens Stoye, and a "smarting up" of the k-mer matching by turning it into
 *     a threaded sort and merge paradigm using a super cache coherent radix sort.  Local
 *     alignment is accomplised with dynamically-banded O(nd) algorithm that terminates when
 *     it fails to find a e-matching patch for a significant distance, and polishes the match
 *     to the last e-prefix-positive 32-mer.
 *
 *  Author :  Gene Myers
 *  First  :  June 2013
 *  Current:  June 1, 2014
 *
 ********************************************************************************************/

//  A complete threaded code for the filter

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "DB.h"
#include "filter.h"
#include "align.h"

#undef FOR_PACBIO

#define THREAD    pthread_t

#define MAX_BIAS  2    //  In -b mode, don't consider tuples with specificity
                       //     <= 4 ^ -(kmer-MAX_BIAS)
#define MAXGRAM 10000  //  Cap on k-mer count histogram (in count_thread, merge_thread)

#define PANEL_SIZE     50000   //  Size to break up very long A-reads
#define PANEL_OVERLAP  10000   //  Overlap of A-panels

#define MATCH_CHUNK    100     //  Max expected number of hits between two reads
#define TRACE_CHUNK  20000     //  Max expected trace points in hits between two reads

#undef  TEST_LSORT
#undef  TEST_KSORT
#undef  TEST_PAIRS
#undef  TEST_CSORT
#define    HOW_MANY   3000   //  Print first HOW_MANY items for each of the TEST options above

#define DO_ALIGNMENT
#undef  TEST_GATHER
#undef  TEST_CONTAIN
#undef  SHOW_OVERLAP          //  Show the cartoon
#undef  SHOW_ALIGNMENT        //  Show the alignment
#define   ALIGN_WIDTH    80   //     Parameters for alignment
#define   ALIGN_INDENT   20
#define   ALIGN_BORDER   10

#ifdef SHOW_OVERLAP
#define NOTHREAD
#endif

#ifdef TEST_GATHER
#define NOTHREAD
#endif

#ifdef TEST_CONTAIN
#define NOTHREAD
#endif

typedef struct
  { uint64 p1;   //  The lower half
    uint64 p2;
  } Double;

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__

typedef struct
  { uint64 code;
    int    rpos;
    int    read;
  } KmerPos;

typedef struct
  { int    diag;
    int    apos;
    int    aread;
    int    bread;
  } SeedPair;

#else

typedef struct
  { uint64 code;
    int    read;
    int    rpos;
  } KmerPos;

typedef struct
  { int    apos;
    int    diag;
    int    bread;
    int    aread;
  } SeedPair;

#endif

/*******************************************************************************************
 *
 *  PARAMETER SETUP
 *
 ********************************************************************************************/

static int Kmer;
static int Hitmin;
static int Binshift;
static int Suppress;

static int    Kshift;         //  2*Kmer
static uint64 Kmask;          //  4^Kmer-1
static int    TooFrequent;    //  (Suppress != 0) ? Suppress : INT32_MAX

static int    NTHREADS;       //  Adjusted downward to nearest power of 2
static int    NSHIFT;         //  NTHREADS = 1 << NSHIFT

int Set_Filter_Params(int kmer, int binshift, int suppress, int hitmin, int nthread)
{ if (kmer <= 1)
    return (1);

  Kmer     = kmer;
  Binshift = binshift;
  Suppress = suppress;
  Hitmin   = hitmin;

  Kshift = 2*Kmer;
  if (Kmer == 32)
    Kmask = 0xffffffffffffffffllu;
  else
    Kmask = (0x1llu << Kshift) - 1;

  if (Suppress == 0)
    TooFrequent = INT32_MAX;
  else
    TooFrequent = Suppress;

  NTHREADS = 1;
  NSHIFT   = 0;
  while (2*NTHREADS <= nthread)
    { NTHREADS *= 2;
      NSHIFT   += 1;
    }

  return (0);
}


/*******************************************************************************************
 *
 *  LEXICOGRAPHIC SORT
 *
 ********************************************************************************************/

#define BMER      4
#define BSHIFT    8             //  = 2*BMER
#define BPOWR   256             //  = 2^BSHIFT
#define BMASK  0xffllu          //  = BPOWR-1

static uint64  QMASK;           //  = BMASK << NSHIFT
static int     LEX_shift;
static int64   LEX_zsize;
static int     LEX_last;
static int     LEX_next;
static Double *LEX_src;
static Double *LEX_trg;

typedef struct
  { int64  beg;
    int64  end;
    int64  tptr[BPOWR];
    int64 *sptr;
  } Lex_Arg;

static void *lex_thread(void *arg)
{ Lex_Arg    *data  = (Lex_Arg *) arg;
  int64      *sptr  = data->sptr;
  int64      *tptr  = data->tptr;
  int         shift = LEX_shift;   //  Must be a multiple of 8 in [0,120]
  int        qshift = (LEX_next - LEX_shift) - NSHIFT;
  int64       zsize = LEX_zsize;
  Double     *src   = LEX_src;
  Double     *trg   = LEX_trg;
  int64       i, n, x;
  uint64      c, b;

  n = data->end;
  if (shift >= 64)
    { shift -= 64;
      if (LEX_last)
        for (i = data->beg; i < n; i++)
          { c = src[i].p2;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p2;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
            sptr[((b >> qshift) & QMASK) + x/zsize] += 1;
          }
    }

  else if ( ! LEX_last && LEX_next >= 64)   //  && LEX_shift < 64

    { qshift = (LEX_next - 64) - NSHIFT;
      if (qshift < 0)
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
            sptr[((src[i].p2 << NSHIFT) & QMASK) + x/zsize] += 1;
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
            sptr[((src[i].p2 >> qshift) & QMASK) + x/zsize] += 1;
          }
    }

  else // LEX_last || LEX_next < 64
    if (LEX_last)
      if (shift == 0)
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            x = tptr[c&BMASK]++;
            trg[x] = src[i];
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
          }
    else
      if (shift == 0)
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            x = tptr[c&BMASK]++;
            trg[x] = src[i];
            sptr[((c >> qshift) & QMASK) + x/zsize] += 1;
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
            sptr[((b >> qshift) & QMASK) + x/zsize] += 1;
          }

  return (NULL);
}

static Double *lex_sort(int bytes[16], Double *src, Double *trg, Lex_Arg *parmx)
{ THREAD  threads[NTHREADS];

  int64   len, x, y;
  Double *xch;
  int     i, j, k, z;
  int     b, c, fb;

  len       = parmx[NTHREADS-1].end;
  LEX_zsize = (len-1)/NTHREADS + 1;
  LEX_src   = src;
  LEX_trg   = trg;
  QMASK     = (BMASK << NSHIFT);

  for (c = 0; c < 16; c++)
    if (bytes[c])
      break;
  fb = c;
  for (b = c; b < 16; b = c)
    { for (c = b+1; c < 16; c++)
        if (bytes[c])
          break;
      LEX_last  = (c >= 16);
      LEX_shift = (b << 3);
      LEX_next  = (c << 3);
 
      if (b == fb)
        { for (i = 0; i < NTHREADS; i++)
            for (z = 0; z < NTHREADS*BPOWR; z++)
              parmx[i].sptr[z] = 0;
        }
      else
        { x = 0;
          for (i = 0; i < NTHREADS; i++)
            { parmx[i].beg = x;
              x = LEX_zsize*(i+1);
              if (x > len)
                x = len;
              parmx[i].end = x;
              for (j = 0; j < BPOWR; j++)
                parmx[i].tptr[j] = 0;
            }
          parmx[NTHREADS-1].end = len;

          for (j = 0; j < BPOWR; j++)
            { k = (j << NSHIFT);
              for (z = 0; z < NTHREADS; z++)
                for (i = 0; i < NTHREADS; i++)
                  { parmx[i].tptr[j] += parmx[z].sptr[k+i];
                    parmx[z].sptr[k+i] = 0;
                  }
            }
        }

      x = 0;
      for (j = 0; j < BPOWR; j++)
        for (i = 0; i < NTHREADS; i++)
          { y = parmx[i].tptr[j];
            parmx[i].tptr[j] = x;
            x += y;
          }

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,lex_thread,parmx+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      xch     = LEX_src;
      LEX_src = LEX_trg;
      LEX_trg = xch;

#ifdef TEST_LSORT
      printf("\nLSORT %d\n",LEX_shift);
      if (LEX_shift >= 64)
        { x = (1 << ((LEX_shift-64)+BSHIFT))-1;
          for (i = 0; i < len; i++)
            { printf("%6d: %8llx %8llx %8llx %8llx : %4llx",
                     i,LEX_src[i].p2>>32,(LEX_src[i].p2)&0xffffffffll,LEX_src[i].p1>>32,
                     LEX_src[i].p1&0xffffffffll,LEX_src[i].p2&x);
              if (i > 0 && (LEX_src[i].p1 < LEX_src[i].p1 ||
                             (LEX_src[i].p1 == LEX_src[i].p1 && 
                             (LEX_src[i].p2 & x) < (LEX_src[i-1].p2 & x))))
                printf(" OO");
              printf("\n");
            }
        }
      else
        { x = (1 << (LEX_shift+BSHIFT))-1;
          for (i = 0; i < len; i++)
            { printf("%6d: %8llx %8llx %8llx %8llx : %4llx",
                     i,LEX_src[i].p2>>32,(LEX_src[i].p2)&0xffffffffll,LEX_src[i].p1>>32,
                     LEX_src[i].p1&0xffffffffll,LEX_src[i].p1&x);
              if (i > 0 && (LEX_src[i].p1 & x) < (LEX_src[i-1].p1 & x))
                printf(" OO");
              printf("\n");
            }
        }
#endif
    }

  return (LEX_src);
}


/*******************************************************************************************
 *
 *  INDEX BUILD
 *
 ********************************************************************************************/

static int *NormShift = NULL;
static int  LogNorm, LogThresh;
static int  LogBase[4];

static DAZZ_DB    *TA_block;
static KmerPos    *TA_list;
static DAZZ_TRACK *TA_track;

typedef struct
  { int    tnum;
    int64 *kptr;
    int    fill;
  } Tuple_Arg;


static void *tuple_thread(void *arg)
{ Tuple_Arg  *data  = (Tuple_Arg *) arg;
  int         tnum  = data->tnum;
  int64      *kptr  = data->kptr;
  KmerPos    *list  = TA_list;
  int         i, m, n, x, p;
  uint64      c;
  char       *s;

  c  = TA_block->nreads;
  i  = (c * tnum) >> NSHIFT;
  n  = TA_block->reads[i].boff;
  s  = ((char *) (TA_block->bases)) + n;
  n -= Kmer*i;

  if (TA_track != NULL)
    { DAZZ_READ *reads = TA_block->reads;
      int64     *anno1 = ((int64 *) (TA_track->anno)) + 1;
      int       *point = (int *) (TA_track->data);
      int64      a, b, f; 
      int        q = 0;

      f = anno1[i-1];
      for (m = (c * (tnum+1)) >> NSHIFT; i < m; i++)
        { b = f;
          f = anno1[i];
          for (a = b; a <= f; a += 2)
            { if (a == b)
                p = 0;
              else
                p = point[a-1];
              if (a == f)
                q = reads[i].rlen;
              else
                q = point[a];
              if (p+Kmer <= q)
                { c = 0;
                  for (x = 1; x < Kmer; x++)
                    c = (c << 2) | s[p++];
                  while (p < q)
                    { x = s[p];
                      c = ((c << 2) | x) & Kmask;
                      list[n].read = i;
                      list[n].rpos = p++;
                      list[n].code = c;
                      n += 1;
                      kptr[c & BMASK] += 1;
                    }
                }
            }
          s += (q+1);
        }

      m = TA_block->reads[m].boff - Kmer*m;
      kptr[BMASK] += (data->fill = m-n);
      while (n < m)
        { list[n].code = 0xffffffffffffffffllu;
          list[n].read = -1;
          list[n].rpos = -1;
          n += 1;
        }
    }

  else
    for (m = (c * (tnum+1)) >> NSHIFT; i < m; i++)
      { c = p = 0;
        for (x = 1; x < Kmer; x++)
          c = (c << 2) | s[p++];
        while ((x = s[p]) != 4)
          { c = ((c << 2) | x) & Kmask;
            list[n].read = i;
            list[n].rpos = p++;
            list[n].code = c;
            n += 1;
            kptr[c & BMASK] += 1;
          }
        s += (p+1);
      }

  return (NULL);
}

static void *biased_tuple_thread(void *arg)
{ Tuple_Arg  *data  = (Tuple_Arg *) arg;
  int         tnum  = data->tnum;
  int64      *kptr  = data->kptr;
  KmerPos    *list  = TA_list;
  int         n, i, m;
  int         x, a, k, p;
  uint64      d, c;
  char       *s, *t;

  c  = TA_block->nreads;
  i  = (c * tnum) >> NSHIFT;
  n  = TA_block->reads[i].boff;
  s  = ((char *) (TA_block->bases)) + n;
  n -= Kmer*i;

  if (TA_track != NULL)
    { DAZZ_READ *reads = TA_block->reads;
      int64     *anno1 = ((int64 *) (TA_track->anno)) + 1;
      int       *point = (int *) (TA_track->data);
      int64      j, b, f; 
      int        q = 0;

      f = anno1[i-1];
      for (m = (c * (tnum+1)) >> NSHIFT; i < m; i++)
        { b = f;
          f = anno1[i];
          t = s+1;
          for (j = b; j <= f; j += 2)
            { if (j == b)
                p = 0;
              else
                p = point[j-1];
              if (j == f)
                q = reads[i].rlen;
              else
                q = point[j];
              if (p+Kmer <= q)
                { c = 0;
                  a = 0;
                  k = 1;
                  while (p < q)
                    { x = s[p];
                      a += LogBase[x];
                      c  = ((c << 2) | x);
                      while (a < LogNorm && k < Kmer)
                        { if (++p >= q)
                            break;
                          k += 1;
                          x  = s[p];
                          a += LogBase[x];
                          c  = ((c << 2) | x);
                        }
                      while (1)
	                { int u = a-LogBase[(int) t[p-k]];
                          if (u < LogNorm) break;
                          a  = u;
                          k -= 1;
                        }
                      if (a > LogThresh)
                        { d = ((c << NormShift[k]) & Kmask);
                          list[n].read = i;
                          list[n].rpos = p;
                          list[n].code = d;
                          n += 1;
                          kptr[d & BMASK] += 1;
                        }
                      p += 1;
                      a -= LogBase[(int) s[p-k]];
                    }
                }
            }
          s += (q+1);
	}
    }

  else
    for (m = (c * (tnum+1)) >> NSHIFT; i < m; i++)
      { t = s+1;
        c = 0;
        p = a = 0;
        k = 1;
        while ((x = s[p]) != 4)
          { a += LogBase[x];
            c  = ((c << 2) | x);
            while (a < LogNorm && k < Kmer)
              { if ((x = s[++p]) == 4)
                  goto eoread2;
                k += 1;
                a += LogBase[x];
                c  = ((c << 2) | x);
              }
            while (1)
	      { int u = a-LogBase[(int) t[p-k]];
                if (u < LogNorm) break;
                a  = u;
                k -= 1;
              }
            if (a > LogThresh)
              { d = ((c << NormShift[k]) & Kmask);
                list[n].read = i;
                list[n].rpos = p;
                list[n].code = d;
                n += 1;
                kptr[d & BMASK] += 1;
              }
            p += 1;
            a -= LogBase[(int) s[p-k]];
          }
      eoread2:
        s += (p+1);
      }

  m = TA_block->reads[m].boff - Kmer*m;
  kptr[BMASK] += (data->fill = m-n);
  while (n < m)
    { list[n].code = 0xffffffffffffffffllu;
      list[n].read = -1;
      list[n].rpos = -1;
      n += 1;
    }

  return (NULL);
}

static KmerPos *FR_src;
static KmerPos *FR_trg;

typedef struct
  { int  beg;
    int  end;
    int  kept;
  } Comp_Arg;

static void *compsize_thread(void *arg)
{ Comp_Arg   *data  = (Comp_Arg *) arg;
  int         end   = data->end;
  KmerPos    *src   = FR_src;
  int         n, i, c, p;
  uint64      h, g;

  i = data->beg;
  h = src[i].code;
  n = 0;
  while (i < end)
    { p = i++;
      while ((g = src[i].code) == h)
        i += 1;
      if ((c = (i-p)) < TooFrequent)
        n += c;
      h = g;
    }

  data->kept = n;
  return (NULL);
}

static void *compress_thread(void *arg)
{ Comp_Arg   *data  = (Comp_Arg *) arg;
  int         end   = data->end;
  KmerPos    *src   = FR_src;
  KmerPos    *trg   = FR_trg;
  int         n, i, p;
  uint64      h, g;

  i = data->beg;
  h = src[i].code;
  n = data->kept;
  while (i < end)
    { p = i++;
      while ((g = src[i].code) == h)
        i += 1;
      if (i-p < TooFrequent)
        { while (p < i)
            trg[n++] = src[p++];
        }
      h = g;
    }

  return (NULL);
}

void *Sort_Kmers(DAZZ_DB *block, int *len)
{ THREAD    threads[NTHREADS];
  Tuple_Arg parmt[NTHREADS];
  Comp_Arg  parmf[NTHREADS];
  Lex_Arg   parmx[NTHREADS];
  int       mersort[16];

  KmerPos  *src, *trg, *rez;
  int       kmers, nreads;
  int       i, j, x, z;
  uint64    h;

  for (i = 0; i < NTHREADS; i++)
    parmx[i].sptr = (int64 *) alloca(NTHREADS*BPOWR*sizeof(int64));

  for (i = 0; i < 16; i++)
    mersort[i] = 0;
  for (i = 0; i < Kshift; i += 8)
    mersort[i>>3] = 1;

  if (NormShift == NULL && BIASED)
    { double scale;

      NormShift = (int *) Malloc(sizeof(int)*(Kmer+1),"Allocating Sort_Kmers bias shift");
      if (NormShift == NULL)
        exit (1);
      for (i = 0; i <= Kmer; i++)
        NormShift[i] = Kshift - 2*i;
      LogNorm = 10000 * Kmer;
      LogThresh = 10000 * (Kmer-MAX_BIAS);

      scale = -10000. / log(4.);
      for (i = 0; i < 4; i++)
        LogBase[i] = (int) ceil( scale * log(block->freq[i]) );
    }

  nreads = block->nreads;
  kmers  = block->reads[nreads].boff - Kmer * nreads;

  if (kmers <= 0)
    goto no_mers;

  if (( (Kshift-1)/BSHIFT + (TooFrequent < INT32_MAX) ) & 0x1)
    { trg = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
      src = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
    }
  else
    { src = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
      trg = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
    }
  if (src == NULL || trg == NULL)
    exit (1);

  if (VERBOSE)
    { printf("\n   Kmer count = ");
      Print_Number((int64) kmers,0,stdout);
      printf("\n   Using %.2fGb of space\n",(1. * kmers) / 33554432);
      fflush(stdout);
    }

  TA_block = block;
  TA_list  = src;
  TA_track = block->tracks;

  for (i = 0; i < NTHREADS; i++)
    { parmt[i].tnum = i;
      parmt[i].kptr = parmx[i].tptr;
      for (j = 0; j < BPOWR; j++)
        parmt[i].kptr[j] = 0;
    }

  if (BIASED)
    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,biased_tuple_thread,parmt+i);
  else
    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,tuple_thread,parmt+i);

  for (i = 0; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);

  x = 0;
  for (i = 0; i < NTHREADS; i++)
    { parmx[i].beg = x;
      j = (int) ((((int64) nreads) * (i+1)) >> NSHIFT);
      parmx[i].end = x = block->reads[j].boff - j*Kmer;
    }

  rez = (KmerPos *) lex_sort(mersort,(Double *) src,(Double *) trg,parmx);
  if (BIASED || TA_track != NULL)
    { if (Kmer%4 == 0)
        { int wedge[NTHREADS];

          for (j = 0; j < NTHREADS; j++)
            if (parmt[j].fill > 0)
              break;
          j += 1;
          if (j < NTHREADS)
            { x = kmers-1;
              for (i = NTHREADS-1; i >= j; i--)
                { x = x - parmt[i].fill;
                  z = x;
                  while (rez[x].read >= 0)
                    x -= 1;
                  wedge[i] = z-x;
                }
              x += 1;
              z = x-parmt[j-1].fill;
              for (i = j; i < NTHREADS; i++)
                { memmove(rez+z,rez+x,wedge[i]*sizeof(KmerPos));
                  x += wedge[i] + parmt[i].fill;
                  z += wedge[i];
                }
            }
        }
      for (i = 0; i < NTHREADS; i++)
        kmers -= parmt[i].fill;
    }

  if (TooFrequent < INT32_MAX && kmers > 0)
    { parmf[0].beg = 0;
      for (i = 1; i < NTHREADS; i++)
        { x = (((int64) i)*kmers) >> NSHIFT;
          h = rez[x-1].code;
          while (rez[x].code == h)
            x += 1;
          parmf[i-1].end = parmf[i].beg = x;
        }
      parmf[NTHREADS-1].end = kmers;

      if (rez[kmers-1].code == 0xffffffffffffffffllu)
        rez[kmers].code = 0;
      else
        rez[kmers].code = 0xffffffffffffffffllu;

      if (src == rez)
        { FR_src = src;
          FR_trg = rez = trg;
        }
      else
        { FR_src = trg;
          FR_trg = rez = src;
        }

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,compsize_thread,parmf+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      x = 0;
      for (i = 0; i < NTHREADS; i++)
        { z = parmf[i].kept;
          parmf[i].kept = x;
          x += z;
        }
      kmers = x;

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,compress_thread,parmf+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);
    }

  rez[kmers].code   = 0xffffffffffffffffllu;
  rez[kmers+1].code = 0;
    
  if (src != rez)
    free(src);
  else
    free(trg);

#ifdef TEST_KSORT
  { int i;

    printf("\nKMER SORT:\n");
    for (i = 0; i < HOW_MANY && i < kmers; i++)
      { KmerPos *c = rez+i;
        printf(" %9d:  %6d / %6d / %16llx\n",i,c->read,c->rpos,c->code);
      }
    fflush(stdout);
  }
#endif

  if (VERBOSE)
    { if (TooFrequent < INT32_MAX || BIASED || TA_track != NULL)
        { printf("   Revised kmer count = ");
          Print_Number((int64) kmers,0,stdout);
          printf("\n");
        }
      printf("   Index occupies %.2fGb\n",(1. * kmers) / 67108864);
      fflush(stdout);
    }

  if (kmers <= 0)
    { free(rez);
      goto no_mers;
    }

  if (kmers > (int64) (MEM_LIMIT/(4*sizeof(KmerPos))))
    { fprintf(stderr,"Warning: Block size too big, index occupies more than 1/4 of");
      if (MEM_LIMIT == MEM_PHYSICAL)
        fprintf(stderr," physical memory (%.1fGb)\n",(1.*MEM_LIMIT)/0x40000000ll);
      else
        fprintf(stderr," desired memory allocation (%.1fGb)\n",(1.*MEM_LIMIT)/0x40000000ll);
      fflush(stderr);
    }

  *len = kmers;
  return (rez);

no_mers:
  *len = 0;
  return (NULL);
}


/*******************************************************************************************
 *
 *  FILTER MATCH
 *
 ********************************************************************************************/

static int find_tuple(uint64 x, KmerPos *a, int n)
{ int l, r, m;

  // smallest k s.t. a[k].code >= x (or n if does not exist)

  l = 0;
  r = n;
  while (l < r)
    { m = ((l+r) >> 1);
      if (a[m].code < x)
        l = m+1;
      else
        r = m;
    }
  return (l);
}

  //  Determine what *will* be the size of the merged list and histogram of sizes for given cutoffs

static KmerPos  *MG_alist;
static KmerPos  *MG_blist;
static SeedPair *MG_hits;
static int       MG_comp;
static int       MG_self;

typedef struct
  { int    abeg, aend;
    int    bbeg, bend;
    int64 *kptr;
    int64  nhits;
    int    limit;
    int64  hitgram[MAXGRAM];
  } Merge_Arg;

static void *count_thread(void *arg)
{ Merge_Arg  *data  = (Merge_Arg *) arg;
  KmerPos    *asort = MG_alist;
  KmerPos    *bsort = MG_blist;
  int64      *gram  = data->hitgram;
  int64       nhits = 0;
  int         aend  = data->aend;

  int64  ct;
  int    ia, ib;
  int    jb, ja;
  uint64 ca, cb;
  uint64 da, db;
  int    ar, ap;
  int    a, b;

  ia = data->abeg;
  ca = asort[ia].code;
  ib = data->bbeg;
  cb = bsort[ib].code;
  if (MG_self)
    { while (1)
        { while (cb < ca)
            cb = bsort[++ib].code;
          while (cb > ca)
            ca = asort[++ia].code;
          if (cb == ca)
            { ja = ia++;
              while ((da = asort[ia].code) == ca)
                ia += 1;
              jb = ib++;
              while ((db = bsort[ib].code) == cb)
                ib += 1;

              if (ia > aend)
                { if (ja >= aend)
                    break;
                  da = asort[ia = aend].code;
                  db = bsort[ib = data->bend].code;
                }

              ct = 0;
              b  = jb;
              if (IDENTITY)
                for (a = ja; a < ia; a++)
                  { ar = asort[a].read;
                    if (MG_comp)
                      { while (b < ib && bsort[b].read <= ar)
                          b += 1;
                      }
                    else
                      { ap = asort[a].rpos;
                        while (b < ib && bsort[b].read < ar)
                          b += 1;
                        while (b < ib && bsort[b].read == ar && bsort[b].rpos < ap)
                          b += 1;
                      }
                    ct += (b-jb);
                  }
              else
                for (a = ja; a < ia; a++)
                  { ar = asort[a].read;
                    while (b < ib && bsort[b].read < ar)
                      b += 1;
                    ct += (b-jb);
                  }

              nhits += ct;
              ca = da;
              cb = db;

              if (ct < MAXGRAM)
                gram[ct] += 1;
            }
        }
    }
  else
    { while (1)
        { while (cb < ca)
            cb = bsort[++ib].code;
          while (cb > ca)
            ca = asort[++ia].code;
          if (cb == ca)
            { ja = ia++;
              while ((da = asort[ia].code) == ca)
                ia += 1;
              jb = ib++;
              while ((db = bsort[ib].code) == cb)
                ib += 1;

              if (ia > aend)
                { if (ja >= aend)
                    break;
                  da = asort[ia = aend].code;
                  db = bsort[ib = data->bend].code;
                }

              ct  = (ia-ja);
              ct *= (ib-jb);

              nhits += ct;
              ca = da;
              cb = db;

              if (ct < MAXGRAM)
                gram[ct] += 1;
            }
        }
    }

  data->nhits = nhits;

  return (NULL);
}

  //  Produce the merged list now that the list has been allocated and
  //    the appropriate cutoff determined.

static void *merge_thread(void *arg)
{ Merge_Arg  *data  = (Merge_Arg *) arg;
  int64      *kptr  = data->kptr;
  KmerPos    *asort = MG_alist;
  KmerPos    *bsort = MG_blist;
  SeedPair   *hits  = MG_hits;
  int64       nhits = data->nhits;
  int         aend  = data->aend;
  int         limit = data->limit;

  int64  ct;
  int    ia, ib;
  int    jb, ja;
  uint64 ca, cb;
  uint64 da, db;
  int    ar, ap;
  int    a, b, c;

  ia = data->abeg;
  ca = asort[ia].code;
  ib = data->bbeg;
  cb = bsort[ib].code;
  if (MG_self)
    { while (1)
        { while (cb < ca)
            cb = bsort[++ib].code;
          while (cb > ca)
            ca = asort[++ia].code;
          if (cb == ca)
            { ja = ia++;
              while ((da = asort[ia].code) == ca)
                ia += 1;
              jb = ib++;
              while ((db = bsort[ib].code) == cb)
                ib += 1;

              if (ia > aend)
                { if (ja >= aend)
                    break;
                  da = asort[ia = aend].code;
                  db = bsort[ib = data->bend].code;
                }

              ct = 0;
              b  = jb;
              if (IDENTITY)
                for (a = ja; a < ia; a++)
                  { ar = asort[a].read;
                    if (MG_comp)
                      { while (b < ib && bsort[b].read <= ar)
                          b += 1;
                      }
                    else
                      { ap = asort[a].rpos;
                        while (b < ib && bsort[b].read < ar)
                          b += 1;
                        while (b < ib && bsort[b].read == ar && bsort[b].rpos < ap)
                          b += 1;
                      }
                    ct += (b-jb);
                  }
              else
                for (a = ja; a < ia; a++)
                  { ar = asort[a].read;
                    while (b < ib && bsort[b].read < ar)
                      b += 1;
                    ct += (b-jb);
                  }

              if (ct < limit)
                { b = jb;
                  if (IDENTITY)
                    for (a = ja; a < ia; a++)
                      { ap = asort[a].rpos;
                        ar = asort[a].read;
                        if (MG_comp)
                          { while (b < ib && bsort[b].read <= ar)
                              b += 1;
                          }
                        else
                          { while (b < ib && bsort[b].read < ar)
                              b += 1;
                            while (b < ib && bsort[b].read == ar && bsort[b].rpos < ap)
                              b += 1;
                          }
                        if ((ct = b-jb) > 0)
                          { kptr[ap & BMASK] += ct;
                            for (c = jb; c < b; c++)
                              { hits[nhits].bread = bsort[c].read;
                                hits[nhits].aread = ar;
                                hits[nhits].apos  = ap; 
                                hits[nhits].diag  = ap - bsort[c].rpos;
                                nhits += 1;
                              }
                          }
                      }
                  else
                    for (a = ja; a < ia; a++)
                      { ap = asort[a].rpos;
                        ar = asort[a].read;
                        while (b < ib && bsort[b].read < ar)
                          b += 1;
                        if ((ct = b-jb) > 0)
                          { kptr[ap & BMASK] += ct;
                            for (c = jb; c < b; c++)
                              { hits[nhits].bread = bsort[c].read;
                                hits[nhits].aread = ar;
                                hits[nhits].apos  = ap; 
                                hits[nhits].diag  = ap - bsort[c].rpos;
                                nhits += 1;
                              }
                          }
                      }
                }
              ca = da;
              cb = db;
            }
        }
    }
  else
    { while (1)
        { while (cb < ca)
            cb = bsort[++ib].code;
          while (cb > ca)
            ca = asort[++ia].code;
          if (cb == ca)
            { if (ia >= aend) break;
              ja = ia++;
              while ((da = asort[ia].code) == ca)
                ia += 1;
              jb = ib++;
              while ((db = bsort[ib].code) == cb)
                ib += 1;

              if (ia > aend)
                { if (ja >= aend)
                    break;
                  da = asort[ia = aend].code;
                  db = bsort[ib = data->bend].code;
                }

              ct = ib-jb;
              if ((ia-ja)*ct < limit)
                { for (a = ja; a < ia; a++)
                    { ap = asort[a].rpos;
                      kptr[ap & BMASK] += ct;
                      for (b = jb; b < ib; b++)
                        { hits[nhits].bread = bsort[b].read;
                          hits[nhits].aread = asort[a].read;
                          hits[nhits].apos  = ap;
                          hits[nhits].diag  = ap - bsort[b].rpos;
                          nhits += 1;
                        }
                    }
                }
              ca = da;
              cb = db;
            }
        }
    }

  return (NULL);
}

  //  Report threads: given a segment of merged list, find all seeds and from them all alignments.

static DAZZ_DB    *MR_ablock;
static DAZZ_DB    *MR_bblock;
static SeedPair   *MR_hits;
static int         MR_two;
static Align_Spec *MR_spec;
static int         MR_tspace;

typedef struct
  { uint64   max;
    uint64   top;
    uint16  *trace;
  } Trace_Buffer;

static int Entwine(Path *jpath, Path *kpath, Trace_Buffer *tbuf, int *where)
{ int   ac, b2, y2, ae;
  int   i, j, k;
  int   num, den, min;
#ifdef SEE_ENTWINE
  int   strt = 1;
  int   iflare, oflare;
#endif

  uint16 *ktrace = tbuf->trace + (uint64) (kpath->trace);
  uint16 *jtrace = tbuf->trace + (uint64) (jpath->trace);

  min   = 10000;
  num   = 0;
  den   = 0;

#ifdef SEE_ENTWINE
  printf("\n");
#endif

  y2 = jpath->bbpos;
  j  = jpath->abpos/MR_tspace;

  b2 = kpath->bbpos;
  k  = kpath->abpos/MR_tspace;

  if (jpath->abpos == kpath->abpos)
    { min = abs(y2-b2);
      if (min == 0)
        *where = kpath->abpos;
    }

  if (j < k)
    { ac = k*MR_tspace;

      j = 1 + 2*(k-j);
      k = 1;

      for (i = 1; i < j; i += 2)
        y2 += jtrace[i];
    }
  else
    { ac = j*MR_tspace;

      k = 1 + 2*(j-k);
      j = 1;

      for (i = 1; i < k; i += 2)
        b2 += ktrace[i];
    }

  ae = jpath->aepos;
  if (ae > kpath->aepos)
    ae = kpath->aepos;

  while (1)
    { ac += MR_tspace;
      if (ac >= ae)
        break;
      y2 += jtrace[j];
      b2 += ktrace[k];
      j += 2;
      k += 2;

#ifdef SEE_ENTWINE
      printf("   @ %5d : %5d %5d = %4d\n",ac,y2,b2,abs(b2-y2));
#endif

      i = abs(y2-b2);
      if (i <= min)
        { min = i;
          if (i == 0)
            *where = ac;
        }
      num += i;
      den += 1;
#ifdef SEE_ENTWINE
      if (strt)
        { strt   = 0;
          iflare = i;
        }
      oflare = i;
#endif
    }

  if (jpath->aepos == kpath->aepos)
    { i = abs(jpath->bepos-kpath->bepos);
      if (i <= min)
        { min = i;
          if (i == 0)
            *where = kpath->aepos;
        }
    }

#ifdef SEE_ENTWINE
  if (den == 0)
    printf("Nothing\n");
  else
    printf("MINIM = %d AVERAGE = %d  IFLARE = %d  OFLARE = %d\n",min,num/den,iflare,oflare);
#endif

  if (den == 0)
    return (-1);
  else
    return (min);
}


//  Produce the concatentation of path1 and path2 where they are known to meet at
//    the trace point with coordinate ap. Place this result in a big growing buffer,
//    that gets reset when fusion is called with path1 = NULL

static void Fusion(Path *path1, int ap, Path *path2, Trace_Buffer *tbuf)
{ int     k, k1, k2;
  int     len, diff;
  uint16 *trace;

  k1 = 2 * ((ap/MR_tspace) - (path1->abpos/MR_tspace));
  k2 = 2 * ((ap/MR_tspace) - (path2->abpos/MR_tspace));

  len = k1+(path2->tlen-k2);

  if (tbuf->top + len >= tbuf->max)
    { tbuf->max = 1.2*(tbuf->top+len) + 1000;
      tbuf->trace = (uint16 *) Realloc(tbuf->trace,sizeof(uint16)*tbuf->max,"Allocating paths");
      if (tbuf->trace == NULL)
        exit (1);
    }

  trace = tbuf->trace + tbuf->top;
  tbuf->top += len;

  diff = 0;
  len  = 0;
  if (k1 > 0)
    { uint16 *t = tbuf->trace + (uint64) (path1->trace);
      for (k = 0; k < k1; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }
  if (k2 < path2->tlen)
    { uint16 *t = tbuf->trace + (uint64) (path2->trace);
      for (k = k2; k < path2->tlen; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }

  path1->aepos = path2->aepos;
  path1->bepos = path2->bepos;
  path1->diffs = diff;
  path1->trace = (void *) (trace - tbuf->trace);
  path1->tlen  = len;
}


static int Handle_Redundancies(Path *amatch, int novls, Path *bmatch, Trace_Buffer *tbuf)
{ Path *jpath, *kpath;
  int   j, k, no;
  int   dist;
  int   awhen = 0, bwhen = 0;
  int   hasB;

#ifdef TEST_CONTAIN
  for (j = 0; j < novls; j++)
    printf("  %3d: [%5d,%5d] x [%5d,%5d]\n",j,amatch[j].abpos,amatch[j].aepos,
                                              amatch[j].bbpos,amatch[j].bepos);
#endif

  hasB = (bmatch != NULL);

  for (j = 1; j < novls; j++)
    { jpath = amatch+j;
      for (k = j-1; k >= 0; k--)
        { kpath = amatch+k;

          if (kpath->abpos < 0)
            continue;

          if (jpath->abpos < kpath->abpos)

            { if (kpath->abpos <= jpath->aepos && kpath->bbpos <= jpath->bepos)
                { dist = Entwine(jpath,kpath,tbuf,&awhen);
                  if (dist == 0)
                    { if (kpath->aepos > jpath->aepos)
                        { if (hasB)
                            { if (MG_comp)
                                { dist = Entwine(bmatch+k,bmatch+j,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(jpath,awhen,kpath,tbuf);
                                  Fusion(bmatch+k,bwhen,bmatch+j,tbuf);
                                  bmatch[j] = bmatch[k];
#ifdef TEST_CONTAIN
                                  printf("  Really 1");
#endif
                                }
                              else
                                { dist = Entwine(bmatch+j,bmatch+k,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(jpath,awhen,kpath,tbuf);
                                  Fusion(bmatch+j,bwhen,bmatch+k,tbuf);
#ifdef TEST_CONTAIN
                                  printf("  Really 2");
#endif
                                }
                            }
                          else
                            { Fusion(jpath,awhen,kpath,tbuf);
#ifdef TEST_CONTAIN
                              printf("  Really 3");
#endif
                            }
                          k = j;
                        }
                      kpath->abpos = -1;
#ifdef TEST_CONTAIN
                      printf("  Fuse! A %d %d\n",j,k);
#endif
                    }
                }
            }

          else // kpath->abpos <= jpath->abpos

            { if (jpath->abpos <= kpath->aepos && jpath->bbpos <= kpath->bepos)
                { dist = Entwine(kpath,jpath,tbuf,&awhen);
                  if (dist == 0)
                    { if (kpath->abpos == jpath->abpos)
                        { if (kpath->aepos > jpath->aepos)
                            { *jpath = *kpath;
                              if (hasB)
                                bmatch[j] = bmatch[k];
                            }
                        }
                      else if (jpath->aepos > kpath->aepos)
                        { if (hasB)
                            { if (MG_comp)
                                { dist = Entwine(bmatch+j,bmatch+k,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(kpath,awhen,jpath,tbuf);
                                  *jpath = *kpath;
                                  Fusion(bmatch+j,bwhen,bmatch+k,tbuf);
#ifdef TEST_CONTAIN
                                  printf("  Really 4");
#endif
                                }
                              else
                                { dist = Entwine(bmatch+k,bmatch+j,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(kpath,awhen,jpath,tbuf);
                                  *jpath = *kpath;
                                  Fusion(bmatch+k,bwhen,bmatch+j,tbuf);
                                  bmatch[j] = bmatch[k];
#ifdef TEST_CONTAIN
                                  printf("  Really 5");
#endif
                                }
                            }
                          else
                            { Fusion(kpath,awhen,jpath,tbuf);
                              *jpath = *kpath;
#ifdef TEST_CONTAIN
                              printf("  Really 6");
#endif
                            }
                          k = j;
                        }
                      else
                        { *jpath = *kpath;
                          if (hasB)
                            bmatch[j] = bmatch[k];
                        }
                      kpath->abpos = -1;
#ifdef TEST_CONTAIN
                      printf("  Fuse! B %d %d\n",j,k);
#endif
                    }
                }
            }
        }
    }

  no = 0;
  for (j = 0; j < novls; j++)
    if (amatch[j].abpos >= 0)
      { if (hasB)
          bmatch[no] = bmatch[j];
        amatch[no++] = amatch[j];
      }
  novls = no;

#ifdef TEST_CONTAIN
  for (j = 0; j < novls; j++)
    printf("  %3d: [%5d,%5d] x [%5d,%5d]\n",j,amatch[j].abpos,amatch[j].aepos,
                                              amatch[j].bbpos,amatch[j].bepos);
#endif

  return (novls);
}

void Diagonal_Span(Path *path, int *mind, int *maxd)
{ uint16 *points;
  int     i, tlen;
  int     dd, low, hgh;

  points = path->trace;
  tlen   = path->tlen;

  dd = path->abpos - path->bbpos;
  low = hgh = dd;

  dd = path->aepos - path->bepos;
  if (dd < low)
    low = dd;
  else if (dd > hgh)
    hgh = dd;

  dd = (path->abpos/MR_tspace)*MR_tspace - path->bbpos;
  tlen -= 2;
  for (i = 1; i < tlen; i += 2)
    { dd += MR_tspace - points[i];
      if (dd < low)
        low = dd;
      else if (dd > hgh)
        hgh = dd;
    }

  *mind = (low >> Binshift)-1;
  *maxd = (hgh >> Binshift)+1;
}

typedef struct
  { int64       beg, end;
    int        *score;
    int        *lastp;
    int        *lasta;
    Work_Data  *work;
    FILE       *ofile1;
    FILE       *ofile2;
    int64       nfilt;
    int64       ncheck;
  } Report_Arg;

static void *report_thread(void *arg)
{ Report_Arg  *data   = (Report_Arg *) arg;
  SeedPair    *hits   = MR_hits;
  Double      *hitd   = (Double *) MR_hits;
  char        *aseq   = (char *) (MR_ablock->bases);
  char        *bseq   = (char *) (MR_bblock->bases);
  DAZZ_READ   *aread  = MR_ablock->reads;
  DAZZ_READ   *bread  = MR_bblock->reads;
  int         *score  = data->score;
  int         *scorp  = data->score + 1;
  int         *scorm  = data->score - 1;
  int         *lastp  = data->lastp;
  int         *lasta  = data->lasta;
  Work_Data   *work   = data->work;
  FILE        *ofile1 = data->ofile1;
  FILE        *ofile2 = data->ofile2;
  int          afirst = MR_ablock->tfirst;
  int          bfirst = MR_bblock->tfirst;
  int          maxdiag = ( MR_ablock->maxlen >> Binshift);
  int          mindiag = (-MR_bblock->maxlen >> Binshift);

  Overlap     _ovla, *ovla = &_ovla;
  Overlap     _ovlb, *ovlb = &_ovlb;
  Alignment   _align, *align = &_align;
  Path        *apath = &(ovla->path);
  Path        *bpath;
  int64        nfilt = 0;
  int64        ahits = 0;
  int64        bhits = 0;
  int          small, tbytes;

  int    AOmax, BOmax;
  int    novla, novlb;
  Path  *amatch, *bmatch;

  Trace_Buffer _tbuf, *tbuf = &_tbuf;

  Double *hitc;
  int     minhit;
  uint64  cpair;
  uint64  npair = 0;
  int64   nidx, eidx;

  //  In ovl and align roles of A and B are reversed, as the B sequence must be the
  //    complemented sequence !!

  align->flags = ovla->flags = ovlb->flags = MG_comp;
  align->path  = apath;

  if (MR_tspace <= TRACE_XOVR)
    { small  = 1;
      tbytes = sizeof(uint8);
    }
  else
    { small  = 0;
      tbytes = sizeof(uint16);
    }

  AOmax = BOmax = MATCH_CHUNK;
  amatch = Malloc(sizeof(Path)*AOmax,"Allocating match vector");
  bmatch = Malloc(sizeof(Path)*BOmax,"Allocating match vector");

  tbuf->max   = 2*TRACE_CHUNK;
  tbuf->trace = Malloc(sizeof(short)*tbuf->max,"Allocating trace vector");

  if (amatch == NULL || bmatch == NULL || tbuf->trace == NULL)
    exit (1);

  fwrite(&ahits,sizeof(int64),1,ofile1);
  fwrite(&MR_tspace,sizeof(int),1,ofile1);
  if (MR_two)
    { fwrite(&bhits,sizeof(int64),1,ofile2);
      fwrite(&MR_tspace,sizeof(int),1,ofile2);
    }

  minhit = (Hitmin-1)/Kmer + 1;
  hitc   = hitd + (minhit-1);
  eidx   = data->end - minhit;
  nidx   = data->beg;
  for (cpair = hitd[nidx].p2; nidx <= eidx; cpair = npair)
    if (hitc[nidx].p2 != cpair)
      { nidx += 1;
        while ((npair = hitd[nidx].p2) == cpair)
          nidx += 1;
      }
    else
      { int   ar, br;
        int   alen, blen;
        int   doA, doB;
        int   setaln, amark, amark2;
        int   apos, bpos, diag;
        int64 lidx, sidx;
        int64 f, h2;

        ar = hits[nidx].aread;
        br = hits[nidx].bread;
        alen = aread[ar].rlen;
        blen = bread[br].rlen;
        if (alen < HGAP_MIN && blen < HGAP_MIN)
          { nidx += 1;
            while ((npair = hitd[nidx].p2) == cpair)
              nidx += 1;
            continue;
          }

#ifdef TEST_GATHER
        printf("%5d vs %5d : %5d x %5d\n",br+bfirst,ar+afirst,blen,alen);
#endif
        setaln = 1;
        doA = doB = 0;
        amark2 = 0;
        novla  = novlb = 0;
        tbuf->top = 0;
        for (sidx = nidx; hitd[nidx].p2 == cpair; nidx = h2)
          { amark  = amark2 + PANEL_SIZE;
            amark2 = amark  - PANEL_OVERLAP;

            h2 = lidx = nidx;
            do
              { apos  = hits[nidx].apos;
                npair = hitd[++nidx].p2;
                if (apos <= amark2)
                  h2 = nidx;
              }
            while (npair == cpair && apos <= amark);

            if (nidx-lidx < minhit) continue;

            for (f = lidx; f < nidx; f++)
              { apos = hits[f].apos;
                diag = hits[f].diag >> Binshift;
                if (apos - lastp[diag] >= Kmer)
                  score[diag] += Kmer;
                else
                  score[diag] += apos - lastp[diag];
                lastp[diag] = apos;
              }

#ifdef TEST_GATHER
            printf("  %6lld upto %6d",nidx-lidx,amark);
#endif

            for (f = lidx; f < nidx; f++)
              { apos = hits[f].apos;
                diag = hits[f].diag;
                bpos = apos - diag;
                diag = diag >> Binshift;
                if (apos > lasta[diag] &&
                     (score[diag] + scorp[diag] >= Hitmin || score[diag] + scorm[diag] >= Hitmin))
                  { if (setaln)
                      { setaln = 0;
                        align->aseq = aseq + aread[ar].boff;
                        align->bseq = bseq + bread[br].boff;
                        align->alen = alen;
                        align->blen = blen;
                        ovlb->bread = ovla->aread = ar + afirst;
                        ovlb->aread = ovla->bread = br + bfirst;
#ifdef FOR_PACBIO
                        doA = 1;
                        doB = (SYMMETRIC && (ar != br || !MG_self || !MG_comp));
#else
                        doA = (alen >= HGAP_MIN);
                        doB = (SYMMETRIC && blen >= HGAP_MIN &&
                                   (ar != br || !MG_self || !MG_comp));
#endif
                      }
#ifdef TEST_GATHER
                    else
                      printf("\n                    ");

                    if (scorm[diag] > scorp[diag])
                      printf("  %5d.. x %5d.. %5d (%3d)",
                             bpos,apos,apos-bpos,score[diag]+scorm[diag]);
                    else
                      printf("  %5d.. x %5d.. %5d (%3d)",
                             bpos,apos,apos-bpos,score[diag]+scorp[diag]);
#endif
                    nfilt += 1;

#ifdef DO_ALIGNMENT
                    bpath = Local_Alignment(align,work,MR_spec,apos-bpos,apos-bpos,apos+bpos,-1,-1);

                    { int low, hgh, ae;

                      Diagonal_Span(apath,&low,&hgh);
                      if (diag < low)
                        low = diag;
                      else if (diag > hgh)
                        hgh = diag;
                      ae = apath->aepos;
                      for (diag = low; diag <= hgh; diag++)
                        if (ae > lasta[diag])
                          lasta[diag] = ae;
#ifdef TEST_GATHER
                      printf(" %d - %d @ %d",low,hgh,apath->aepos);
#endif
                    }

                    if ((apath->aepos-apath->abpos) + (apath->bepos-apath->bbpos) >= MINOVER)
                      { if (doA)
                          { if (novla >= AOmax)
                              { AOmax = 1.2*novla + MATCH_CHUNK;
                                amatch = Realloc(amatch,sizeof(Path)*AOmax,
                                                 "Reallocating match vector");
                                if (amatch == NULL)
                                  exit (1);
                              }
                            if (tbuf->top + apath->tlen > tbuf->max)
                              { tbuf->max = 1.2*(tbuf->top+apath->tlen) + TRACE_CHUNK;
                                tbuf->trace = Realloc(tbuf->trace,sizeof(short)*tbuf->max,
                                                      "Reallocating trace vector");
                                if (tbuf->trace == NULL)
                                  exit (1);
                              }
                            amatch[novla] = *apath;
                            amatch[novla].trace = (void *) (tbuf->top);
                            memmove(tbuf->trace+tbuf->top,apath->trace,sizeof(short)*apath->tlen);
                            novla += 1;
                            tbuf->top += apath->tlen;
                          }
                        if (doB)
                          { if (novlb >= BOmax)
                              { BOmax = 1.2*novlb + MATCH_CHUNK;
                                bmatch = Realloc(bmatch,sizeof(Path)*BOmax,
                                                        "Reallocating match vector");
                                if (bmatch == NULL)
                                  exit (1);
                              }
                            if (tbuf->top + bpath->tlen > tbuf->max)
                              { tbuf->max = 1.2*(tbuf->top+bpath->tlen) + TRACE_CHUNK;
                                tbuf->trace = Realloc(tbuf->trace,sizeof(short)*tbuf->max,
                                                      "Reallocating trace vector");
                                if (tbuf->trace == NULL)
                                  exit (1);
                              }
                            bmatch[novlb] = *bpath;
                            bmatch[novlb].trace = (void *) (tbuf->top);
                            memmove(tbuf->trace+tbuf->top,bpath->trace,sizeof(short)*bpath->tlen);
                            novlb += 1;
                            tbuf->top += bpath->tlen;
                          }

#ifdef TEST_GATHER
                        printf("  [%5d,%5d] x [%5d,%5d] = %4d",
                               apath->abpos,apath->aepos,apath->bbpos,apath->bepos,apath->diffs);
#endif
#ifdef SHOW_OVERLAP
                        printf("\n\n                    %d(%d) vs %d(%d)\n\n",
                               ovla->aread,ovla->alen,ovla->bread,ovla->blen);
                        Print_ACartoon(stdout,align,ALIGN_INDENT);
#ifdef SHOW_ALIGNMENT
                        Compute_Trace_ALL(align,work);
                        printf("\n                      Diff = %d\n",align->path->diffs);
                        Print_Alignment(stdout,align,work,
                                        ALIGN_INDENT,ALIGN_WIDTH,ALIGN_BORDER,0,5);
#endif
#endif // SHOW_OVERLAP

                      }
#ifdef TEST_GATHER
                    else
                      printf("  No alignment %d",
                              ((apath->aepos-apath->abpos) + (apath->bepos-apath->bbpos))/2);
#endif
#endif // DO_ALIGNMENT
                  }
              }

            for (f = lidx; f < nidx; f++)
              { diag = hits[f].diag >> Binshift;
                score[diag] = lastp[diag] = 0;
              }
#ifdef TEST_GATHER
            printf("\n");
#endif
          }

        for (f = sidx; f < nidx; f++)
          { int d;

            diag = hits[f].diag >> Binshift;
            for (d = diag; d <= maxdiag; d++)
              if (lasta[d] == 0)
                break;
              else
                lasta[d] = 0;
            for (d = diag-1; d >= mindiag; d--)
              if (lasta[d] == 0)
                break;
              else
                lasta[d] = 0;
          }

         
         { int i;

#ifdef TEST_CONTAIN
           if (novla > 1 || novlb > 1)
             printf("\n%5d vs %5d:\n",ar,br);
#endif

           if (novla > 1)
             { if (novlb > 1)
                 novla = novlb = Handle_Redundancies(amatch,novla,bmatch,tbuf);
               else
                 novla = Handle_Redundancies(amatch,novla,NULL,tbuf);
             }
           else if (novlb > 1)
             novlb = Handle_Redundancies(bmatch,novlb,NULL,tbuf);

           for (i = 0; i < novla; i++)
             { ovla->path = amatch[i];
               ovla->path.trace = tbuf->trace + (uint64) (ovla->path.trace);
               if (small)
                 Compress_TraceTo8(ovla);
               if (Write_Overlap(ofile1,ovla,tbytes))
                 { fprintf(stderr,"%s: Cannot write to %s too small?\n",SORT_PATH,Prog_Name);
                   exit (1);
                 }
             }
           for (i = 0; i < novlb; i++)
             { ovlb->path = bmatch[i];
               ovlb->path.trace = tbuf->trace + (uint64) (ovlb->path.trace);
               if (small)
                 Compress_TraceTo8(ovlb);
               if (Write_Overlap(ofile2,ovlb,tbytes))
                 { fprintf(stderr,"%s: Cannot write to %s, too small?\n",SORT_PATH,Prog_Name);
                   exit (1);
                 }
             }
           ahits += novla;
           bhits += novlb;
         }
      }

  free(tbuf->trace);
  free(bmatch);
  free(amatch);

  data->nfilt  = nfilt;
  data->ncheck = ahits + bhits;

  if (MR_two)
    { rewind(ofile2);
      fwrite(&bhits,sizeof(int64),1,ofile2);
      fclose(ofile2);
    }
  else
    ahits += bhits;

  rewind(ofile1);
  fwrite(&ahits,sizeof(int64),1,ofile1);
  fclose(ofile1);

  return (NULL);
}


/*******************************************************************************************
 *
 *  THE ALGORITHM
 *
 ********************************************************************************************/

static char *NameBuffer(char *aname, char *bname)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  len = strlen(aname) + strlen(bname) + 100;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((cat = (char *) realloc(cat,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making path name)\n",Prog_Name);
          exit (1);
        }
    }
  return (cat);
}

void Match_Filter(char *aname, DAZZ_DB *ablock, char *bname, DAZZ_DB *bblock,
                  void *vasort, int alen, void *vbsort, int blen,
                  int comp, Align_Spec *aspec)
{ THREAD     threads[NTHREADS];
  Merge_Arg  parmm[NTHREADS];
  Lex_Arg    parmx[NTHREADS];
  Report_Arg parmr[NTHREADS];
  int        pairsort[16];
  char      *fname;

  SeedPair *khit, *hhit;
  SeedPair *work1, *work2;
  int64     nhits;
  int64     nfilt, ncheck;

  KmerPos  *asort, *bsort;
  int64     atot, btot;

  asort = (KmerPos *) vasort;
  bsort = (KmerPos *) vbsort;

  atot = ablock->totlen;
  btot = bblock->totlen;

  MR_tspace = Trace_Spacing(aspec);

  { int64 powr;
    int   i, nbyte;

    for (i = 0; i < NTHREADS; i++)
      parmx[i].sptr = (int64 *) alloca(NTHREADS*BPOWR*sizeof(int64));

    for (i = 0; i < 16; i++)
      pairsort[i] = 0;

    powr = 1;
    for (nbyte = 0; powr < ablock->maxlen; nbyte += 1)
      powr <<= 8;
    for (i = 4; i < 4+nbyte; i++)
      pairsort[i] = 1;

    powr = 1;
    for (nbyte = 0; powr < ablock->nreads; nbyte += 1)
      powr <<= 8;
    for (i = 8; i < 8+nbyte; i++)
      pairsort[i] = 1;

    powr = 1;
    for (nbyte = 0; powr < bblock->nreads; nbyte += 1)
      powr <<= 8;
    for (i = 12; i < 12+nbyte; i++)
      pairsort[i] = 1;
  }

  nfilt = ncheck = nhits = 0;

  if (VERBOSE)
    { if (comp)
        printf("\nComparing %s to c(%s)\n",aname,bname);
      else
        printf("\nComparing %s to %s\n",aname,bname);
    }

  if (alen == 0 || blen == 0)
    goto zerowork;

  { int    i, j, p;
    uint64 c;
    int    limit;

    MG_alist = asort;
    MG_blist = bsort;
    MG_self  = (aname == bname);
    MG_comp  = comp;

    parmm[0].abeg = parmm[0].bbeg = 0;
    for (i = 1; i < NTHREADS; i++)
      { p = (int) ((((int64) alen) * i) >> NSHIFT);
        if (p > 0)
          { c = asort[p-1].code;
            while (asort[p].code == c)
              p += 1;
          }
        parmm[i].abeg = parmm[i-1].aend = p;
        parmm[i].bbeg = parmm[i-1].bend = find_tuple(asort[p].code,bsort,blen);
      }
    parmm[NTHREADS-1].aend = alen;
    parmm[NTHREADS-1].bend = blen;

    for (i = 0; i < NTHREADS; i++)
      for (j = 0; j < MAXGRAM; j++)
        parmm[i].hitgram[j] = 0;

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,count_thread,parmm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    if (VERBOSE)
      printf("\n");
    if (MEM_LIMIT > 0)
      { int64 histo[MAXGRAM];
        int64 tom, avail;

        for (j = 0; j < MAXGRAM; j++)
          histo[j] = parmm[0].hitgram[j];
        for (i = 1; i < NTHREADS; i++)
          for (j = 0; j < MAXGRAM; j++)
            histo[j] += parmm[i].hitgram[j];

        avail = (int64) (MEM_LIMIT - (sizeof_DB(ablock) + sizeof_DB(bblock))) / sizeof(Double);
        if (asort == bsort || avail > alen + 2*blen)
          avail = (avail - alen) / 2;
        else
          avail = avail - (alen + blen);
        avail *= .98;

        tom = 0;
        for (j = 0; j < MAXGRAM; j++)
          { tom += j*histo[j];
            if (tom > avail)
              break;
          }
        limit = j;

        if (limit <= 1)
          { fprintf(stderr,"\nError: Insufficient ");
            if (MEM_LIMIT == MEM_PHYSICAL)
              fprintf(stderr," physical memory (%.1fGb), reduce block size\n",
                             (1.*MEM_LIMIT)/0x40000000ll);
            else
              { fprintf(stderr," memory allocation (%.1fGb),",(1.*MEM_LIMIT)/0x40000000ll);
                fprintf(stderr," reduce block size or increase allocation\n");
              }
            fflush(stderr);
            exit (1);
          }
        if (limit < 10)
          { fprintf(stderr,"\nWarning: Sensitivity hampered by low ");
            if (MEM_LIMIT == MEM_PHYSICAL)
              fprintf(stderr," physical memory (%.1fGb), reduce block size\n",
                             (1.*MEM_LIMIT)/0x40000000ll);
            else
              { fprintf(stderr," memory allocation (%.1fGb),",(1.*MEM_LIMIT)/0x40000000ll);
                fprintf(stderr," reduce block size or increase allocation\n");
              }
            fflush(stderr);
          }
        if (VERBOSE)
          { printf("   Capping mutual k-mer matches over %d (effectively -t%d)\n",
                   limit,(int) sqrt(1.*limit));
            fflush(stdout);
          }

        for (i = 0; i < NTHREADS; i++)
          { parmm[i].nhits = 0;
            for (j = 1; j < limit; j++)
              parmm[i].nhits += j * parmm[i].hitgram[j];
            parmm[i].limit = limit;
          }
      }
    else
      for (i = 0; i < NTHREADS; i++)
        parmm[i].limit = INT32_MAX;

    nhits = parmm[0].nhits;
    for (i = 1; i < NTHREADS; i++)
      parmm[i].nhits = nhits += parmm[i].nhits;

    if (VERBOSE)
      { printf("   Hit count = ");
        Print_Number(nhits,0,stdout);
        if (asort == bsort || nhits >= blen)
          printf("\n   Highwater of %.2fGb space\n",
                       (1. * (alen + 2*nhits)) / 67108864);
        else
          printf("\n   Highwater of %.2fGb space\n",
                       (1. * (alen + blen + nhits)) / 67108864);
        fflush(stdout);
      }

    if (nhits == 0)
      goto zerowork;

    if (asort == bsort)
      hhit = work1 = (SeedPair *) Malloc(sizeof(SeedPair)*(nhits+1),
                                         "Allocating daligner hit vectors");
    else
      { if (nhits >= blen)
          bsort = (KmerPos *) Realloc(bsort,sizeof(SeedPair)*(nhits+1),
                                       "Reallocating daligner sort vectors");
        hhit = work1 = (SeedPair *) bsort;
      }
    khit = work2 = (SeedPair *) Malloc(sizeof(SeedPair)*(nhits+1),
                                        "Allocating daligner hit vectors");
    if (hhit == NULL || khit == NULL || bsort == NULL)
      exit (1);

    MG_blist = bsort;
    MG_hits  = khit;

    for (i = NTHREADS-1; i > 0; i--)
      parmm[i].nhits = parmm[i-1].nhits;
    parmm[0].nhits = 0;

    for (i = 0; i < NTHREADS; i++)
      { parmm[i].kptr = parmx[i].tptr;
        for (p = 0; p < BPOWR; p++)
          parmm[i].kptr[p] = 0;
      }

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,merge_thread,parmm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

#ifdef TEST_PAIRS
    printf("\nSETUP SORT:\n");
    for (i = 0; i < HOW_MANY && i < nhits; i++)
      { SeedPair *c = khit+i;
        printf(" %5d / %5d / %5d /%5d\n",c->aread,c->bread,c->apos,c->apos-c->diag);
      }
#endif
  }

  { int   i;
    int64 x;

    x = 0;
    for (i = 0; i < NTHREADS-1; i++)
      { parmx[i].beg = x;
        parmx[i].end = x = parmm[i+1].nhits;
      }
    parmx[NTHREADS-1].beg = x;
    parmx[NTHREADS-1].end = nhits;

    khit = (SeedPair *) lex_sort(pairsort,(Double *) khit,(Double *) hhit,parmx);

    khit[nhits].aread = 0x7fffffff;
    khit[nhits].bread = 0x7fffffff;
    khit[nhits].diag  = 0x7fffffff;
    khit[nhits].apos  = 0;

#ifdef TEST_CSORT
    printf("\nCROSS SORT %lld:\n",nhits);
    for (i = 0; i < HOW_MANY && i <= nhits; i++)
      { SeedPair *c = khit+i;
        printf(" %5d / %5d / %5d /%5d\n",c->aread,c->bread,c->apos,c->apos-c->diag);
      }
#endif
  }

  { int    i, w;
    int64  p;
    int    d;
    int   *counters;

    MR_ablock = ablock;
    MR_bblock = bblock;
    MR_hits   = khit;
    MR_two    = ! MG_self && SYMMETRIC;
    MR_spec   = aspec;

    parmr[0].beg = 0;
    for (i = 1; i < NTHREADS; i++)
      { p = (nhits * i) >> NSHIFT;
        if (p > 0)
          { d = khit[p-1].bread;
            while ((khit[p].bread) == d)
              p += 1;
          }
        parmr[i].beg = parmr[i-1].end = p;
      }
    parmr[NTHREADS-1].end = nhits;

    w = ((ablock->maxlen >> Binshift) - ((-bblock->maxlen) >> Binshift)) + 1;
    counters = (int *) Malloc(NTHREADS*3*w*sizeof(int),"Allocating diagonal buckets");
    if (counters == NULL)
      exit (1);

    fname = NameBuffer(aname,bname);

    for (i = 0; i < 3*w*NTHREADS; i++)
      counters[i] = 0;
    for (i = 0; i < NTHREADS; i++)
      { if (i == 0)
          parmr[i].score = counters - ((-bblock->maxlen) >> Binshift);
        else
          parmr[i].score = parmr[i-1].lasta + w;
        parmr[i].lastp = parmr[i].score + w;
        parmr[i].lasta = parmr[i].lastp + w;
        parmr[i].work  = New_Work_Data();

        sprintf(fname,"%s/%s.%s.%c%d.las",SORT_PATH,aname,bname,(comp?'C':'N'),i+1);
        parmr[i].ofile1 = Fopen(fname,"w");
        if (parmr[i].ofile1 == NULL)
          exit (1);
        if (MG_self)
          parmr[i].ofile2 = parmr[i].ofile1;
        else if (SYMMETRIC)
          { sprintf(fname,"%s/%s.%s.%c%d.las",SORT_PATH,bname,aname,(comp?'C':'N'),i+1);
            parmr[i].ofile2 = Fopen(fname,"w");
            if (parmr[i].ofile2 == NULL)
              exit (1);
          }
      }

#ifdef NOTHREAD

    for (i = 0; i < NTHREADS; i++)
      report_thread(parmr+i);

#else

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,report_thread,parmr+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

#endif

    if (VERBOSE)
      for (i = 0; i < NTHREADS; i++)
        { nfilt  += parmr[i].nfilt;
          ncheck += parmr[i].ncheck;
        }

    for (i = 0; i < NTHREADS; i++)
      Free_Work_Data(parmr[i].work);
    free(counters);
  }

  free(work2);
  free(work1);
  goto epilogue;

zerowork:
  { FILE *ofile;
    int   i;

    fname = NameBuffer(aname,bname);

    nhits  = 0;
    for (i = 0; i < NTHREADS; i++)
      { sprintf(fname,"%s/%s.%s.%c%d.las",SORT_PATH,aname,bname,(comp?'C':'N'),i+1);
        ofile = Fopen(fname,"w");
        fwrite(&nhits,sizeof(int64),1,ofile);
        fwrite(&MR_tspace,sizeof(int),1,ofile);
        fclose(ofile);
        if (! MG_self && SYMMETRIC)
          { sprintf(fname,"%s/%s.%s.%c%d.las",SORT_PATH,bname,aname,(comp?'C':'N'),i+1);
            ofile = Fopen(fname,"w");
            fwrite(&nhits,sizeof(int64),1,ofile);
            fwrite(&MR_tspace,sizeof(int),1,ofile);
            fclose(ofile);
          }
      }
  }

epilogue:

  if (VERBOSE)
    { int width;

      if (nhits <= 0)
        width = 1;
      else
        width = ((int) log10((double) nhits)) + 1;
      width += (width-1)/3;

      printf("\n     ");
      Print_Number(nhits,width,stdout);
      printf(" %d-mers (%e of matrix)\n     ",Kmer,(1.*nhits/atot)/btot);
      Print_Number(nfilt,width,stdout);
      printf(" seed hits (%e of matrix)\n     ",(1.*nfilt/atot)/btot);
      Print_Number(ncheck,width,stdout);
      printf(" confirmed hits (%e of matrix)\n",(1.*ncheck/atot)/btot);
      fflush(stdout);
    }
}
