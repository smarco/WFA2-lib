/*******************************************************************************************
 *
 *  Given a list of sorted .las files, merge them into a single sorted .las file.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"

static char *Usage = "[-va] <merge:las> <parts:las> ...";

#define MEMORY 4000   // in Mb

#undef   DEBUG

  //  Heap sort of records according to (aread,bread,COMP(flags),abpos) order

#define COMPARE(lp,rp)				\
  if (lp->aread > rp->aread)			\
    bigger = 1;					\
  else if (lp->aread < rp->aread)		\
    bigger = 0;					\
  else if (lp->bread > rp->bread)		\
    bigger = 1;					\
  else if (lp->bread < rp->bread)		\
    bigger = 0;					\
  else if (COMP(lp->flags) > COMP(rp->flags))	\
    bigger = 1;					\
  else if (COMP(lp->flags) < COMP(rp->flags))	\
    bigger = 0;					\
  else if (lp->path.abpos > rp->path.abpos)	\
    bigger = 1;					\
  else if (lp->path.abpos < rp->path.abpos)	\
    bigger = 0;					\
  else if (lp > rp)				\
    bigger = 1;					\
  else						\
    bigger = 0;

static void reheap(int s, Overlap **heap, int hsize)
{ int      c, l, r;
  int      bigger;
  Overlap *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      if (r > hsize)
        bigger = 1;
      else
        { hr = heap[r];
          COMPARE(hr,hl)
        }
      if (bigger)
        { COMPARE(hs,hl)
          if (bigger)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { COMPARE(hs,hr)
          if (bigger)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

  //  Heap sort of records according to (aread,abpos) order

#define MAPARE(lp,rp)				\
  if (lp->aread > rp->aread)			\
    bigger = 1;					\
  else if (lp->aread < rp->aread)		\
    bigger = 0;					\
  else if (lp->path.abpos > rp->path.abpos)	\
    bigger = 1;					\
  else if (lp->path.abpos < rp->path.abpos)	\
    bigger = 0;					\
  else if (lp > rp)				\
    bigger = 1;					\
  else						\
    bigger = 0;

static void maheap(int s, Overlap **heap, int hsize)
{ int      c, l, r;
  int      bigger;
  Overlap *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      if (r > hsize)
        bigger = 1;
      else
        { hr = heap[r];
          MAPARE(hr,hl)
        }
      if (bigger)
        { MAPARE(hs,hl)
          if (bigger)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { MAPARE(hs,hr)
          if (bigger)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

#ifdef DEBUG

static void showheap(Overlap **heap, int hsize)
{ int i;
  printf("\n");
  for (i = 1; i <= hsize; i++)
    printf(" %3d: %5d, %5d\n",i,heap[i]->aread,heap[i]->bread);
}

#endif

  //  Input block data structure and block fetcher

typedef struct
  { FILE   *stream;
    char   *block;
    char   *ptr;
    char   *top;
    int64   count;
  } IO_block;

static void ovl_reload(IO_block *in, int64 bsize)
{ int64 remains;

  remains = in->top - in->ptr;
  if (remains > 0)
    memmove(in->block, in->ptr, remains);
  in->ptr  = in->block;
  in->top  = in->block + remains;
  in->top += fread(in->top,1,bsize-remains,in->stream);
}

  //  The program

int main(int argc, char *argv[])
{ IO_block *in;
  int64     bsize, osize, psize;
  char     *block, *oblock;
  int       i, fway;
  Overlap **heap;
  int       hsize;
  Overlap  *ovls;
  int64     totl;
  int       tspace, tbytes;
  FILE     *output;
  char     *optr, *otop;

  int       VERBOSE;
  int       MAP_SORT;

  //  Process command line

  { int j, k;
    int flags[128];

    ARG_INIT("LAmerge")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("va") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];
    MAP_SORT = flags['a'];

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }

    fway = argc-2;
    if (fway > 252)
      { fprintf(stderr,"Exceeded maximum # of inputs and outputs (252) of merge\n");
        exit (1);
      }
  }

  //  Open all the input files and initialize their buffers

  psize  = sizeof(void *);
  osize  = sizeof(Overlap) - psize;
  bsize  = (MEMORY*1000000ll)/(fway + 1);
  block  = (char *) Malloc(bsize*(fway+1)+psize,"Allocating LAmerge blocks");
  in     = (IO_block *) Malloc(sizeof(IO_block)*fway,"Allocating LAmerge IO-reacords");
  if (block == NULL || in == NULL)
    exit (1);
  block += psize;

  totl   = 0;
  tbytes = 0;
  tspace = 0;
  for (i = 0; i < fway; i++)
    { int64  novl;
      int    mspace;
      FILE  *input;
      char  *pwd, *root;
      char  *iblock;

      pwd   = PathTo(argv[i+2]);
      root  = Root(argv[i+2],".las");
      input = Fopen(Catenate(pwd,"/",root,".las"),"r");
      if (input == NULL)
        exit (1);
      free(pwd);
      free(root);

      if (fread(&novl,sizeof(int64),1,input) != 1)
        SYSTEM_READ_ERROR
      totl += novl;
      if (fread(&mspace,sizeof(int),1,input) != 1)
        SYSTEM_READ_ERROR
      if (i == 0)
        { tspace = mspace;
          if (tspace <= TRACE_XOVR && tspace != 0)
            tbytes = sizeof(uint8);
          else
            tbytes = sizeof(uint16);
        }
      else if (tspace != mspace)
        { fprintf(stderr,"%s: PT-point spacing conflict (%d vs %d)\n",Prog_Name,tspace,mspace);
          exit (1);
        }

      in[i].stream = input;
      in[i].block  = iblock = block+i*bsize;
      in[i].ptr    = iblock;
      in[i].top    = iblock + fread(in[i].block,1,bsize,input);
      in[i].count  = 0;
    }

  //  Open the output file buffer and write (novl,tspace) header

  { char *pwd, *root;

    pwd    = PathTo(argv[1]);
    root   = Root(argv[1],".las");
    output = Fopen(Catenate(pwd,"/",root,".las"),"w");
    if (output == NULL)
      exit (1);
    free(pwd);
    free(root);

    if (fwrite(&totl,sizeof(int64),1,output) != 1)
      SYSTEM_READ_ERROR
    if (fwrite(&tspace,sizeof(int),1,output) != 1)
      SYSTEM_READ_ERROR

    oblock = block+fway*bsize;
    optr   = oblock;
    otop   = oblock + bsize;
  }

  if (VERBOSE)
    { printf("Merging %d files totalling ",fway);
      Print_Number(totl,0,stdout);
      printf(" records\n");
    }

  //  Initialize the heap

  heap = (Overlap **) Malloc(sizeof(Overlap *)*(fway+1),"Allocating heap");
  ovls = (Overlap *) Malloc(sizeof(Overlap)*fway,"Allocating heap");
  if (heap == NULL || ovls == NULL)
    exit (1);

  hsize = 0;
  for (i = 0; i < fway; i++)
    { if (in[i].ptr < in[i].top)
        { ovls[i]     = *((Overlap *) (in[i].ptr - psize));
          in[i].ptr  += osize;
          hsize      += 1;
          heap[hsize] = ovls + i;
        }
    }

  if (hsize > 3)
    { if (MAP_SORT)
        for (i = hsize/2; i > 1; i--)
          maheap(i,heap,hsize);
      else
        for (i = hsize/2; i > 1; i--)
          reheap(i,heap,hsize);
    }

  //  While the heap is not empty do

  while (hsize > 0)
    { Overlap  *ov;
      IO_block *src;
      int64     tsize, span;

      if (MAP_SORT)
        maheap(1,heap,hsize);
      else
        reheap(1,heap,hsize);

      ov  = heap[1];
      src = in + (ov - ovls);

      do
        { src->count += 1;

          tsize = ov->path.tlen*tbytes;
          span  = osize + tsize;
          if (src->ptr + span > src->top)
            ovl_reload(src,bsize);
          if (optr + span > otop)
            { if (fwrite(oblock,1,optr-oblock,output) != (size_t) (optr-oblock))
                SYSTEM_READ_ERROR
              optr = oblock;
            }

          memmove(optr,((char *) ov) + psize,osize);
          optr += osize;
          memmove(optr,src->ptr,tsize);
          optr += tsize;

          src->ptr += tsize;
          if (src->ptr >= src->top)
            { heap[1] = heap[hsize];
              hsize  -= 1;
              break;
            }
          *ov       = *((Overlap *) (src->ptr - psize));
          src->ptr += osize;
        }
      while (CHAIN_NEXT(ov->flags));
    }

  //  Flush output buffer and wind up

  if (optr > oblock)
    { if (fwrite(oblock,1,optr-oblock,output) != (size_t) (optr-oblock))
        SYSTEM_READ_ERROR
    }
  fclose(output);

  for (i = 0; i < fway; i++)
    fclose(in[i].stream);

  for (i = 0; i < fway; i++)
    totl -= in[i].count;
  if (totl != 0)
    { fprintf(stderr,"%s: Did not write all records to %s (%lld)\n",argv[0],argv[1],totl);
      exit (1);
    }

  free(ovls);
  free(heap);
  free(in);
  free(block-psize);

  exit (0);
}
