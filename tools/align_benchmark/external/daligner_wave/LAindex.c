/*******************************************************************************************
 *
 *  Create an index with extension .las.idx for a .las file.
 *    Utility expects the .las file to be sorted.
 *    Header contains total # of trace points, max # of trace points for
 *    a given overlap, max # of trace points in all the overlaps for a given aread, and
 *    max # of overlaps for a given aread.  The remainder are the offsets into each pile.
 *
 *  Author:  Gene Myers
 *  Date  :  Sept 2015
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

static char *Usage = "[-v] <source:las> ...";

#define MEMORY   1000   //  How many megabytes for output buffer

int main(int argc, char *argv[])
{ char     *iblock;
  FILE     *input, *output;
  int64     novl, bsize, ovlsize, ptrsize;
  int       tspace, tbytes;
  char     *pwd, *root;
  int64     tmax, ttot;
  int64     omax, smax;
  int64     odeg, sdeg;
  int       i;

  int       VERBOSE;

  //  Process options

  { int j, k;
    int flags[128];

    ARG_INIT("LAindex")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("v") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  //  For each file do

  ptrsize = sizeof(void *);
  ovlsize = sizeof(Overlap) - ptrsize;
  bsize   = MEMORY * 1000000ll;
  iblock  = (char *) Malloc(bsize + ptrsize,"Allocating input block");
  if (iblock == NULL)
    exit (1);
  iblock += ptrsize;

  for (i = 1; i < argc; i++)
    { pwd   = PathTo(argv[i]);
      root  = Root(argv[i],".las");
      input = Fopen(Catenate(pwd,"/",root,".las"),"r");
      if (input == NULL)
        exit (1);
    
      if (fread(&novl,sizeof(int64),1,input) != 1)
        SYSTEM_READ_ERROR
      if (fread(&tspace,sizeof(int),1,input) != 1)
        SYSTEM_READ_ERROR
      if (tspace <= TRACE_XOVR && tspace != 0)
        tbytes = sizeof(uint8);
      else
        tbytes = sizeof(uint16);
    
      output = Fopen(Catenate(pwd,"/.",root,".las.idx"),"w");
      if (output == NULL)
        exit (1);
    
      free(pwd);
      free(root);

      if (VERBOSE)
        { printf("  Indexing %s: ",root);
          Print_Number(novl,0,stdout);
          printf(" records ... ");
          fflush(stdout);
        }

      fwrite(&novl,sizeof(int64),1,output);
      fwrite(&novl,sizeof(int64),1,output);
      fwrite(&novl,sizeof(int64),1,output);
      fwrite(&novl,sizeof(int64),1,output);
    
      { int         j, alst;
        Overlap    *w;
        int64       tsize;
        int64       optr;
        char       *iptr, *itop;
        int64       tlen;
    
        optr = sizeof(int64) + sizeof(int32);
        iptr = iblock;
        itop = iblock + fread(iblock,1,bsize,input);
    
        alst = -1;
        odeg = sdeg = 0;
        omax = smax = 0;
        tmax = ttot = 0;
        for (j = 0; j < novl; j++)
          { if (iptr + ovlsize > itop)
              { int64 remains = itop-iptr;
                if (remains > 0)
                  memmove(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,input);
              }
    
            w = (Overlap *) (iptr - ptrsize);
    
            tlen = w->path.tlen;
            if (alst < 0)
              { fwrite(&optr,sizeof(int64),1,output);
                alst = w->aread;
              }
            else
              while (alst < w->aread)
                { if (sdeg > smax)
                    smax = sdeg;
                  if (odeg > omax)
                    omax = odeg;
                  fwrite(&optr,sizeof(int64),1,output);
    	          odeg = sdeg = 0;
                  alst += 1;
                }
            if (tlen > tmax)
              tmax = tlen;
            ttot += tlen;
            odeg += 1;
            sdeg += tlen;
    
            iptr += ovlsize;
    
            tsize = tlen*tbytes;
            if (iptr + tsize > itop)
              { int64 remains = itop-iptr;
                if (remains > 0)
                  memmove(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,input);
              }
            
            optr += ovlsize + tsize;
            iptr += tsize;
          }
        fwrite(&optr,sizeof(int64),1,output);
      }
    
      if (sdeg > smax)
        smax = sdeg;
      if (odeg > omax)
        omax = odeg;
    
      rewind(output);
    
      fwrite(&omax,sizeof(int64),1,output);
      fwrite(&ttot,sizeof(int64),1,output);
      fwrite(&smax,sizeof(int64),1,output);
      fwrite(&tmax,sizeof(int64),1,output);

      if (VERBOSE)
        { Print_Number(ttot,0,stdout);
          printf(" trace points\n");
          fflush(stdout);
        }
    
      fclose(input);
      fclose(output);
    }

  free(iblock-ptrsize);

  exit (0);
}
