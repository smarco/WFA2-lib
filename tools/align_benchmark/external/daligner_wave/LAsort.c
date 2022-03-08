/*******************************************************************************************
 *
 *  Load a file U.las of overlaps into memory, sort them all by A,B index,
 *    and then output the result to U.S.las
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

static char *Usage = "[-va] <align:las> ...";

#define MEMORY   1000   //  How many megabytes for output buffer

static char *IBLOCK;

static int SORT_OVL(const void *x, const void *y)
{ int64 l = *((int64 *) x);
  int64 r = *((int64 *) y);

  Overlap *ol, *or;
  int      al, ar;
  int      bl, br;
  int      cl, cr;
  int      pl, pr;

  ol = (Overlap *) (IBLOCK+l);
  or = (Overlap *) (IBLOCK+r);

  al = ol->aread;
  ar = or->aread;
  if (al != ar)
    return (al-ar);

  bl = ol->bread;
  br = or->bread;
  if (bl != br)
    return (bl-br);

  cl = COMP(ol->flags);
  cr = COMP(or->flags);
  if (cl != cr)
    return (cl-cr);

  pl = ol->path.abpos;
  pr = or->path.abpos;
  if (pl != pr)
    return (pl-pr);

  if (ol < or)
    return (-1);
  else if (ol > or)
    return (1);
  else
    return (0);
}

static int SORT_MAP(const void *x, const void *y)
{ int64 l = *((int64 *) x);
  int64 r = *((int64 *) y);

  Overlap *ol, *or;
  int      al, ar;
  int      pl, pr;

  ol = (Overlap *) (IBLOCK+l);
  or = (Overlap *) (IBLOCK+r);

  al = ol->aread;
  ar = or->aread;
  if (al != ar)
    return (al-ar);

  pl = ol->path.abpos;
  pr = or->path.abpos;
  if (pl != pr)
    return (pl-pr);

  if (ol < or)
    return (-1);
  else if (ol > or)
    return (1);
  else
    return (0);
}

int main(int argc, char *argv[])
{ char     *iblock, *fblock, *iend;
  int64     isize,   osize;
  int64     ovlsize, ptrsize;
  int       tspace, tbytes;
  int       i;

  int       VERBOSE;
  int       MAP_ORDER;
 
  //  Process options

  { int j, k;
    int flags[128];

    ARG_INIT("LAsort")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("va") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];
    MAP_ORDER = flags['a'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  //  For each file do

  ptrsize = sizeof(void *);
  ovlsize = sizeof(Overlap) - ptrsize;
  isize   = 0;
  iblock  = NULL;
  osize   = MEMORY * 1000000ll;
  fblock  = Malloc(osize,"Allocating LAsort output block");

  for (i = 1; i < argc; i++)
    { int64    *perm;
      FILE     *input, *foutput;
      int64     novl, sov;

      //  Read in the entire file and output header

      { int64  size;
        struct stat info;
        char  *pwd, *root, *name;

        pwd   = PathTo(argv[i]);
        root  = Root(argv[i],".las");
        name  = Catenate(pwd,"/",root,".las");
        input = Fopen(name,"r");
        if (input == NULL)
          exit (1);

        stat(name,&info);
        size = info.st_size;

        if (fread(&novl,sizeof(int64),1,input) != 1)
          SYSTEM_READ_ERROR
        if (fread(&tspace,sizeof(int),1,input) != 1)
          SYSTEM_READ_ERROR

        if (tspace <= TRACE_XOVR && tspace != 0)
          tbytes = sizeof(uint8);
        else
          tbytes = sizeof(uint16);

        if (VERBOSE)
          { printf("  %s: ",root);
            Print_Number(novl,0,stdout);
            printf(" records ");
            Print_Number(size-novl*ovlsize,0,stdout);
            printf(" trace bytes\n");
            fflush(stdout);
          }

        foutput = Fopen(Catenate(pwd,"/",root,".S.las"),"w");
        if (foutput == NULL)
          exit (1);

        if (fwrite(&novl,sizeof(int64),1,foutput) != 1)
          SYSTEM_READ_ERROR
        if (fwrite(&tspace,sizeof(int),1,foutput) != 1)
          SYSTEM_READ_ERROR

        free(pwd);
        free(root);

        if (size > isize)
          { if (iblock == NULL)
              iblock = Malloc(size+ptrsize,"Allocating LAsort input block");
            else
              iblock = Realloc(iblock-ptrsize,size+ptrsize,"Allocating LAsort input block");
            if (iblock == NULL)
              exit (1);
            iblock += ptrsize;
            isize   = size;
          }
        size -= (sizeof(int64) + sizeof(int));
        if (size > 0)
          { if (fread(iblock,size,1,input) != 1)
              SYSTEM_READ_ERROR
          }
        fclose(input);
        iend = iblock + (size - ptrsize);
      }

      //  Set up unsorted permutation array

      perm = (int64 *) Malloc(sizeof(int64)*novl,"Allocating LAsort permutation vector");
      if (perm == NULL)
        exit (1);

      { int64 off;
        int   j;

        if (CHAIN_START(((Overlap *) (iblock-ptrsize))->flags))
          { sov = 0;
            off = -ptrsize;
            for (j = 0; j < novl; j++)
              { if (CHAIN_START(((Overlap *) (iblock+off))->flags))
                  perm[sov++] = off;
                off += ovlsize + ((Overlap *) (iblock+off))->path.tlen*tbytes;
              }
          }
        else
          { off = -ptrsize;
            for (j = 0; j < novl; j++)
              { perm[j] = off;
                off += ovlsize + ((Overlap *) (iblock+off))->path.tlen*tbytes;
              }
            sov = novl;
          }
      }

      //  Sort permutation array of ptrs to records

      IBLOCK = iblock;
      if (MAP_ORDER)
        qsort(perm,sov,sizeof(int64),SORT_MAP);
      else
        qsort(perm,sov,sizeof(int64),SORT_OVL);

      //  Output the records in sorted order

      { int      j;
        Overlap *w;
        int64    tsize, span;
        char    *fptr, *ftop, *wo;

        fptr = fblock;
        ftop = fblock + osize;
        for (j = 0; j < sov; j++)
          { w = (Overlap *) (wo = iblock+perm[j]);
            do
              { tsize = w->path.tlen*tbytes;
                span  = ovlsize + tsize;
                if (fptr + span > ftop)
                  { if (fwrite(fblock,1,fptr-fblock,foutput) != (size_t) (fptr-fblock))
                      SYSTEM_READ_ERROR
                    fptr = fblock;
                  }
                memmove(fptr,((char *) w)+ptrsize,ovlsize);
                fptr += ovlsize;
                memmove(fptr,(char *) (w+1),tsize);
                fptr += tsize;
                w = (Overlap *) (wo += span);
              }
            while (wo < iend && CHAIN_NEXT(w->flags));
          }
        if (fptr > fblock)
          { if (fwrite(fblock,1,fptr-fblock,foutput) != (size_t) (fptr-fblock))
              SYSTEM_READ_ERROR
          }
      }

      free(perm);
      fclose(foutput);
    }

  if (iblock != NULL)
    free(iblock - ptrsize);
  free(fblock);

  exit (0);
}
