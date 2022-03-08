/*******************************************************************************************
 *
 *  Utility for displaying the information in the overlaps of a .las file in a very
 *    simple to parse format.
 *
 *  Author:    Gene Myers
 *  Creation:  July 2013
 *  Last Mod:  Jan 2015
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"

static char *Usage =
    "[-cdtlo] <src1:db|dam> [<src2:db|dam>] <align:las> [<reads:FILE> | <reads:range> ...]";

#define LAST_READ_SYMBOL  '$'

static int ORDER(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  return (x-y);
}

int main(int argc, char *argv[])
{ DAZZ_DB   _db1, *db1 = &_db1; 
  DAZZ_DB   _db2, *db2 = &_db2; 
  Overlap   _ovl, *ovl = &_ovl;

  FILE   *input;
  int64   novl;
  int     tspace, tbytes, small;
  int     tmax;
  int     reps, *pts;
  int     input_pts;

  int     OVERLAP;
  int     DOCOORDS, DODIFFS, DOTRACE, DOLENS;
  int     ISTWO;

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("LAdump")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("ocdtl")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    OVERLAP   = flags['o'];
    DOCOORDS  = flags['c'];
    DODIFFS   = flags['d'];
    DOTRACE   = flags['t'];
    DOLENS    = flags['l'];

    if (DOTRACE)
      DOCOORDS = 1;

    if (argc <= 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      P #a #b #o #c     -");
        fprintf(stderr," (#a,#b^#o) have an LA between them where #o is 'n' or 'c' and\n");
        fprintf(stderr,"                         ");
        fprintf(stderr,"   #c is '>' (start of best chain), '+' (start of alternate chain),\n");
        fprintf(stderr,"                         ");
        fprintf(stderr,"   '-' (continuation of chain), or '.' (no chains in file).\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -c: C #ab #ae #bb #be - #a[#ab,#ae] aligns with #b^#o[#bb,#be]\n");
        fprintf(stderr,"      -d: D #               - there are # differences in the LA\n");
        fprintf(stderr,"      -t: T #n              -");
        fprintf(stderr," there are #n trace point intervals for the LA\n");
        fprintf(stderr,"           (#d #y )^#n      -");
        fprintf(stderr," there are #d difference aligning the #y bp's of B with the\n");
        fprintf(stderr,"                                 next fixed-size interval of A\n");
        fprintf(stderr,"      -l: L #la #lb         -");
        fprintf(stderr," #la is the length of the a-read and #lb that of the b-read\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -o: Output proper overlaps only\n");

        exit (1);
      }
  }

  //  Open trimmed DB or DB pair

  { int   status;
    char *pwd, *root;
    FILE *input;

    ISTWO  = 0;
    status = Open_DB(argv[1],db1);
    if (status < 0)
      exit (1);
    if (db1->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }

    if (argc > 3)
      { pwd   = PathTo(argv[3]);
        root  = Root(argv[3],".las");
        if ((input = fopen(Catenate(pwd,"/",root,".las"),"r")) != NULL)
          { ISTWO = 1;
            fclose(input);
            status = Open_DB(argv[2],db2);
            if (status < 0)
              exit (1);
            if (db2->part > 0)
              { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[2]);
                exit (1);
              }
            Trim_DB(db2);
          }
        else
          db2 = db1;
        free(root);
        free(pwd);
      }
    else
      db2 = db1;
    Trim_DB(db1);
  }

  //  Process read index arguments into a sorted list of read ranges

  input_pts = 0;
  if (argc == ISTWO+4)
    { if (argv[ISTWO+3][0] != LAST_READ_SYMBOL || argv[ISTWO+3][1] != '\0')
        { char *eptr, *fptr;
          int   b, e;

          b = strtol(argv[ISTWO+3],&eptr,10);
          if (eptr > argv[ISTWO+3] && b > 0)
            { if (*eptr == '-')
                { if (eptr[1] != LAST_READ_SYMBOL || eptr[2] != '\0')
                    { e = strtol(eptr+1,&fptr,10);
                      input_pts = (fptr <= eptr+1 || *fptr != '\0' || e <= 0);
                    }
                }
              else
                input_pts = (*eptr != '\0');
            }
          else
            input_pts = 1;
        }
    }

  if (input_pts)
    { int v, x;
      FILE *input;

      input = Fopen(argv[ISTWO+3],"r");
      if (input == NULL)
        exit (1);

      reps = 0;
      while ((v = fscanf(input," %d",&x)) != EOF)
        if (v == 0)
          { fprintf(stderr,"%s: %d'th item of input file %s is not an integer\n",
                           Prog_Name,reps+1,argv[2]);
            exit (1);
          }
        else
          reps += 1;

      reps *= 2;
      pts   = (int *) Malloc(sizeof(int)*reps,"Allocating read parameters");
      if (pts == NULL)
        exit (1);

      rewind(input);
      for (v = 0; v < reps; v += 2)
        { fscanf(input," %d",&x);
          pts[v] = pts[v+1] = x;
        }

      fclose(input);
    }

  else
    { pts  = (int *) Malloc(sizeof(int)*2*argc,"Allocating read parameters");
      if (pts == NULL)
        exit (1);

      reps = 0;
      if (argc > 3+ISTWO)
        { int   c, b, e;
          char *eptr, *fptr;

          for (c = 3+ISTWO; c < argc; c++)
            { if (argv[c][0] == LAST_READ_SYMBOL)
                { b = db1->nreads;
                  eptr = argv[c]+1;
                }
              else
                b = strtol(argv[c],&eptr,10);
              if (eptr > argv[c])
                { if (b <= 0)
                    { fprintf(stderr,"%s: %d is not a valid index\n",Prog_Name,b);
                      exit (1);
                    }
                  if (*eptr == '\0')
                    { pts[reps++] = b;
                      pts[reps++] = b;
                      continue;
                    }
                  else if (*eptr == '-')
                    { if (eptr[1] == LAST_READ_SYMBOL)
                        { e = INT32_MAX;
                          fptr = eptr+2;
                        }
                      else
                        e = strtol(eptr+1,&fptr,10);
                      if (fptr > eptr+1 && *fptr == 0 && e > 0)
                        { pts[reps++] = b;
                          pts[reps++] = e;
                          if (b > e)
                            { fprintf(stderr,"%s: Empty range '%s'\n",Prog_Name,argv[c]);
                              exit (1);
                            }
                          continue;
                        }
                    }
                }
              fprintf(stderr,"%s: argument '%s' is not an integer range\n",Prog_Name,argv[c]);
              exit (1);
            }

          qsort(pts,reps/2,sizeof(int64),ORDER);

          b = 0;
          for (c = 0; c < reps; c += 2)
            if (b > 0 && pts[b-1] >= pts[c]-1) 
              { if (pts[c+1] > pts[b-1])
                  pts[b-1] = pts[c+1];
              }
            else
              { pts[b++] = pts[c];
                pts[b++] = pts[c+1];
              }
          pts[b++] = INT32_MAX;
          reps = b;
        }
      else
        { pts[reps++] = 1;
          pts[reps++] = INT32_MAX;
        }
    }

  //  Initiate file reading and read header
  
  { char  *over, *pwd, *root;

    pwd   = PathTo(argv[2+ISTWO]);
    root  = Root(argv[2+ISTWO],".las");
    over  = Catenate(pwd,"/",root,".las");
    input = Fopen(over,"r");
    if (input == NULL)
      exit (1);

    if (fread(&novl,sizeof(int64),1,input) != 1)
      SYSTEM_READ_ERROR
    if (fread(&tspace,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR

    if (tspace <= TRACE_XOVR && tspace != 0)
      { small  = 1;
        tbytes = sizeof(uint8);
      }
    else
      { small  = 0;
        tbytes = sizeof(uint16);
      }

    free(pwd);
    free(root);
  }

  //  Scan to count sizes of things

  { int   j, al, tlen;
    int   in, npt, idx, ar;
    int64 novls, odeg, omax, sdeg, smax, ttot;

    in  = 0;
    npt = pts[0];
    idx = 1;

    //  For each record do

    novls = omax = smax = ttot = tmax = 0;
    sdeg  = odeg = 0;

    al = 0;
    for (j = 0; j < novl; j++)

       //  Read it in

      { Read_Overlap(input,ovl);
        tlen = ovl->path.tlen;
        fseeko(input,tlen*tbytes,SEEK_CUR);

        //  Determine if it should be displayed

        ar = ovl->aread+1;
        if (in)
          { while (ar > npt)
              { npt = pts[idx++];
                if (ar < npt)
                  { in = 0;
                    break;
                  }
                npt = pts[idx++];
              }
          }
        else
          { while (ar >= npt)
              { npt = pts[idx++];
                if (ar <= npt)
                  { in = 1;
                    break;
                  }
                npt = pts[idx++];
              }
          }
        if (!in)
          continue;

        //  If -o check display only overlaps

        if (OVERLAP)
          { if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
              continue;
            if (ovl->path.aepos != db1->reads[ovl->aread].rlen &&
                ovl->path.bepos != db2->reads[ovl->bread].rlen)
              continue;
          }

        if (ar != al)
          { if (sdeg > smax)
              smax = sdeg;
            if (odeg > omax)
              omax = odeg;
            sdeg = odeg = 0;
            al = ar;
          }

        novls += 1;
        odeg  += 1;
        sdeg  += tlen;
        ttot  += tlen;
        if (tlen > tmax)
          tmax = tlen;
      }

    if (sdeg > smax)
      smax = sdeg;
    if (odeg > omax)
      omax = odeg;

    printf("+ P %lld\n",novls);
    printf("%% P %lld\n",omax);
    if (DOTRACE)
      { printf("+ T %lld\n",ttot);
        printf("%% T %lld\n",smax);
        printf("@ T %d\n",tmax);
      }
  }

  //  Read the file and display selected records
  
  { int        j, k;
    uint16    *trace;
    int        in, npt, idx, ar;
    DAZZ_READ *read1, *read2;

    rewind(input);
    fread(&novl,sizeof(int64),1,input);
    fread(&tspace,sizeof(int),1,input);

    trace = (uint16 *) Malloc(sizeof(uint16)*tmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);

    read1 = db1->reads;
    read2 = db2->reads;

    in  = 0;
    npt = pts[0];
    idx = 1;

    //  For each record do

    for (j = 0; j < novl; j++)

       //  Read it in

      { Read_Overlap(input,ovl);
        ovl->path.trace = (void *) trace;
        Read_Trace(input,ovl,tbytes);

        //  Determine if it should be displayed

        ar = ovl->aread+1;
        if (in)
          { while (ar > npt)
              { npt = pts[idx++];
                if (ar < npt)
                  { in = 0;
                    break;
                  }
                npt = pts[idx++];
              }
          }
        else
          { while (ar >= npt)
              { npt = pts[idx++];
                if (ar <= npt)
                  { in = 1;
                    break;
                  }
                npt = pts[idx++];
              }
          }
        if (!in)
          continue;

        //  If -o check display only overlaps

        if (OVERLAP)
          { if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
              continue;
            if (ovl->path.aepos != db1->reads[ovl->aread].rlen &&
                ovl->path.bepos != db2->reads[ovl->bread].rlen)
              continue;
          }

        //  Display it
            
        printf("P %d %d",ovl->aread+1,ovl->bread+1);
        if (COMP(ovl->flags))
          printf(" c");
        else
          printf(" n");
        if (CHAIN_NEXT(ovl->flags))
          printf(" -");
        else if (BEST_CHAIN(ovl->flags))
          printf(" >");
        else if (CHAIN_START(ovl->flags))
          printf(" +");
        else
          printf(" .");
        printf("\n");

        if (DOLENS)
          printf("L %d %d\n",read1[ovl->aread].rlen,read2[ovl->bread].rlen);

        if (DOCOORDS)
          printf("C %d %d %d %d\n",ovl->path.abpos,ovl->path.aepos,ovl->path.bbpos,ovl->path.bepos);

        if (DODIFFS)
          printf("D %d\n",ovl->path.diffs);

        if (DOTRACE)
          { uint16 *trace = (uint16 *) ovl->path.trace;
            int     tlen  = ovl->path.tlen;

            if (small)
              Decompress_TraceTo16(ovl);
            printf("T %d\n",tlen>>1);
            for (k = 0; k < tlen; k += 2)
              printf(" %3d %3d\n",trace[k],trace[k+1]);
          }
      }

    free(trace);
  }

  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  exit (0);
}
