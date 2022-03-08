/*******************************************************************************************
 *
 *  Merge together in index order, overlap files <XXX>.1.las, <XXX>.2.las, ... into a
 *    single overlap file and output to the standard output
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

static char *Usage = "[-v] <source:las> > <target>.las";

#define MEMORY   1000   //  How many megabytes for output buffer

int main(int argc, char *argv[])
{ char     *iblock, *oblock;
  FILE     *input;
  int64     novl, bsize, ovlsize, ptrsize;
  int       tspace, tbytes;
  char     *pwd, *root, *root2;

  int       VERBOSE;

  //  Process options

  { int i, j, k;
    int flags[128];

    ARG_INIT("LAcat")

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

  ptrsize = sizeof(void *);
  ovlsize = sizeof(Overlap) - ptrsize;
  bsize   = MEMORY * 1000000ll;
  oblock  = (char *) Malloc(bsize,"Allocating output block");
  iblock  = (char *) Malloc(bsize + ptrsize,"Allocating input block");
  if (oblock == NULL || iblock == NULL)
    exit (1);
  iblock += ptrsize;

  pwd    = PathTo(argv[1]);
  root   = Root(argv[1],".las");

  root2 = index(root,'#');
  if (root2 == NULL)
    { fprintf(stderr,"%s: No #-sign in source name '%s'\n",Prog_Name,root);
      exit (1);
    }
  if (index(root2+1,'#') != NULL)
    { fprintf(stderr,"%s: Two or more occurences of #-sign in source name '%s'\n",Prog_Name,root);
      exit (1);
    }
  *root2++ = '\0';

  { int64    povl;
    int      i, mspace;

    novl   = 0;
    tspace = 0;
    mspace = 0;
    tbytes = 0;
    for (i = 0; 1; i++)
      { char *name = Catenate(pwd,"/",Numbered_Suffix(root,i+1,root2),".las");
        if ((input = fopen(name,"r")) == NULL) break;

        if (fread(&povl,sizeof(int64),1,input) != 1)
          SYSTEM_READ_ERROR
        novl += povl;
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

        fclose(input);
      }
    if (fwrite(&novl,sizeof(int64),1,stdout) != 1)
      SYSTEM_READ_ERROR
    if (fwrite(&tspace,sizeof(int),1,stdout) != 1)
      SYSTEM_READ_ERROR
  }

  { int      i, j;
    Overlap *w;
    int64    tsize, povl;
    int      mspace;
    char    *iptr, *itop;
    char    *optr, *otop;

    optr = oblock;
    otop = oblock + bsize;

    for (i = 0; 1; i++)
      { char *name = Catenate(pwd,"/",Numbered_Suffix(root,i+1,root2),".las");
        if ((input = fopen(name,"r")) == NULL) break;

        if (fread(&povl,sizeof(int64),1,input) != 1)
          SYSTEM_READ_ERROR
        if (fread(&mspace,sizeof(int),1,input) != 1)
          SYSTEM_READ_ERROR

        if (VERBOSE)
          fprintf(stderr,"  Concatenating %s: %lld la\'s\n",Numbered_Suffix(root,i+1,root2),povl);

        iptr = iblock;
        itop = iblock + fread(iblock,1,bsize,input);

        for (j = 0; j < povl; j++)
          { if (iptr + ovlsize > itop)
              { int64 remains = itop-iptr;
                if (remains > 0)
                  memmove(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,input);
              }

            w = (Overlap *) (iptr - ptrsize);
            tsize = w->path.tlen*tbytes;

            if (optr + ovlsize + tsize > otop)
              { if (fwrite(oblock,1,optr-oblock,stdout) != (size_t) (optr-oblock))
                  SYSTEM_READ_ERROR
                optr = oblock;
              }

            memmove(optr,iptr,ovlsize);
            optr += ovlsize;
            iptr += ovlsize;

            if (iptr + tsize > itop)
              { int64 remains = itop-iptr;
                if (remains > 0)
                  memmove(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,input);
              }
            
            memmove(optr,iptr,tsize);
            optr += tsize;
            iptr += tsize;
          }

        fclose(input);
      }

    if (optr > oblock)
      { if (fwrite(oblock,1,optr-oblock,stdout) != (size_t) (optr-oblock))
          SYSTEM_READ_ERROR
      }
  }

  if (VERBOSE)
    fprintf(stderr,"  Totalling %lld la\'s\n",novl);

  free(pwd);
  free(root);
  free(oblock);
  free(iblock-ptrsize);

  exit (0);
}
