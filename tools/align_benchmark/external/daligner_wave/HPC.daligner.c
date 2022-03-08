/*********************************************************************************************\
 *
 *  Produce a script to compute overlaps for all block pairs of a DB, and then sort and merge
 *    them into as many .las files as their are blocks.
 *
 *  Author:  Gene Myers
 *  Date  :  June 1, 2014
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "DB.h"
#include "filter.h"

#undef  LSF  //  define if want a directly executable LSF script

static char *Usage[] =
  { "[-vbad] [-t<int>] [-w<int(6)>] [-l<int(1000)>] [-s<int(100)] [-P<dir(/tmp)>]",
    "        [-M<int>] [-B<int(4)>] [-D<int( 250)>] [-T<int(4)>] [-f<name>]",
    "      ( [-k<int(14)>] [-h<int(35)>] [-e<double(.70)>] [-H<int>]",
    "        [-k<int(20)>] [-h<int(50)>] [-e<double(.85)>] <ref:db|dam> )",
    "        [-m<track>]+ <reads:db|dam> [<first:int>[-<last:int>]]"
  };

  //  Command Options

static int    DUNIT, BUNIT;
static int    VON, BON, CON, DON;
static int    WINT, TINT, HGAP, HINT, KINT, SINT, LINT, MINT;
static int    NTHREADS;
static double EREL;
static int    MMAX, MTOP;
static char **MASK;
static char  *ONAME;
static char  *PDIR;

#define LSF_ALIGN "bsub -q medium -n 4 -o DALIGNER.out -e DALIGNER.err -R span[hosts=1] -J align#%d"
#define LSF_MERGE \
          "bsub -q short -n 12 -o MERGE%d.DAL.out -e MERGE%d.DAL.err -R span[hosts=1] -J merge#%d"
#define LSF_CHECK \
          "bsub -q short -n 12 -o CHECK%d.DAL.out -e CHECK%d.DAL.err -R span[hosts=1] -J check#%d"

void daligner_script(int argc, char *argv[])
{ int   nblocks;
  int   usepath;
  int   useblock;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  FILE *out;
  char  name[100];
  char *pwd, *root;

  //  Make sure DB exists and is partitioned, get number of blocks in partition

  pwd = PathTo(argv[1]);
  if (strcmp(argv[1]+(strlen(argv[1])-4),".dam") == 0)
    root = Root(argv[1],".dam");
  else
    root = Root(argv[1],".db");

  { int i, nfiles;
    FILE *dbvis;

    dbvis = fopen(Catenate(pwd,"/",root,".dam"),"r");
    if (dbvis == NULL)
      { dbvis = Fopen(Catenate(pwd,"/",root,".db"),"r");
        if (dbvis == NULL)
          exit (1);
      }

    if (fscanf(dbvis,"files = %d\n",&nfiles) != 1)
      SYSTEM_READ_ERROR
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];

        if (fgets(buffer,30000,dbvis) == NULL)
          SYSTEM_READ_ERROR
      }

    useblock = 1;
    if (fscanf(dbvis,"blocks = %d\n",&nblocks) != 1 || nblocks == 1)
      { useblock = 0;
        nblocks  = 1;
      }

    usepath = (strcmp(pwd,".") != 0);
  }

  //  Set range fblock-lblock checking that DB.<fblock-1>.las exists & DB.<fblock>.las does not

  { char *eptr, *fptr;
    FILE *file;

    if (argc == 3)
      { fblock = strtol(argv[2],&eptr,10);
        if (*eptr != '\0' && *eptr != '-')
          { fprintf(stderr,"%s: final argument '%s' does not start with an integer\n",
                           Prog_Name,argv[2]);
            exit (1);
          }
        useblock = 1;
        if (*eptr == '-')
          { lblock = strtol(eptr+1,&fptr,10);
            if (*fptr != '\0')
              { fprintf(stderr,"%s: second part of range '%s' is not an integer\n",
                               Prog_Name,eptr+1);
                exit (1);
              }
          }
        else
          lblock = fblock;
        if (fblock < 1 || lblock > nblocks || fblock > lblock)
          { fprintf(stderr,"%s: range %d-%d is empty or out of bounds\n",Prog_Name,fblock,lblock);
            exit (1);
          }
      }
    else
      { fblock = 1;
        lblock = nblocks;
      }

    if (fblock > 1)
      { file = fopen(Catenate(pwd,"/",root,Numbered_Suffix(".",fblock-1,".las")),"r");
        if (file == NULL)
          { if (usepath)
              fprintf(stderr,"%s: File %s/%s.%d.las should already be present!\n",
                             Prog_Name,pwd,root,fblock-1);
            else
              fprintf(stderr,"%s: File %s.%d.las should already be present!\n",
                             Prog_Name,root,fblock-1);
            exit (1);
          }
        else
          fclose(file);
      }
    if (useblock)
      file = fopen(Catenate(pwd,"/",root,Numbered_Suffix(".",fblock,".las")),"r");
    else
      file = fopen(Catenate(pwd,"/",root,".las"),"r");
    if (file != NULL)
      { if (usepath)
          if (useblock)
            fprintf(stderr,"%s: File %s/%s.%d.las should not yet exist!\n",
                           Prog_Name,pwd,root,fblock);
          else
            fprintf(stderr,"%s: File %s/%s.las should not yet exist!\n",Prog_Name,pwd,root);
        else
          if (useblock)
            fprintf(stderr,"%s: File %s.%d.las should not yet exist!\n",Prog_Name,root,fblock);
          else
            fprintf(stderr,"%s: File %s.las should not yet exist!\n",Prog_Name,root);
        exit (1);
      }

    DON = (DON && (lblock > 1));
    out = stdout;
  }

  { int level, njobs;
    int i, j, k;

    //  Create all work subdirectories if DON

    if (DON)
      { if (ONAME != NULL)
          { sprintf(name,"%s.00.MKDIR",ONAME);
            out = fopen(name,"w");
          }

        fprintf(out,"# Create work subdirectories\n");
        for (i = fblock; i <= lblock; i++)
          fprintf(out,"mkdir work%d\n",i);

        if (ONAME != NULL)
          fclose(out);
      }

    //  Produce all necessary daligner jobs

    if (ONAME != NULL)
      { sprintf(name,"%s.01.OVL",ONAME);
        out = fopen(name,"w");
      }

    njobs = 0;
    for (i = fblock; i <= lblock; i++)
      njobs += (i-1)/BUNIT+1;

    fprintf(out,"# Daligner jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      { int bits;
        int low, hgh;

        bits = (i-1)/BUNIT+1;
        low  = 1;
        for (j = 1; j <= bits; j++)
          {
#ifdef LSF
            fprintf(out,LSF_ALIGN,jobid++);
            fprintf(out," \"");
#endif
            fprintf(out,"daligner");
            if (VON)
              fprintf(out," -v");
            if (BON)
              fprintf(out," -b");
            if (KINT != 14)
              fprintf(out," -k%d",KINT);
            if (WINT != 6)
              fprintf(out," -w%d",WINT);
            if (HINT != 35)
              fprintf(out," -h%d",HINT);
            if (TINT > 0)
              fprintf(out," -t%d",TINT);
            if (HGAP > 0)
              fprintf(out," -H%d",HGAP);
            if (EREL > 0.)
              fprintf(out," -e%g",EREL);
            if (LINT != 1000)
              fprintf(out," -l%d",LINT);
            if (SINT != 100)
              fprintf(out," -s%d",SINT);
            if (MINT >= 0)
              fprintf(out," -M%d",MINT);
            if (PDIR != NULL)
              fprintf(out," -P%s",PDIR);
            if (NTHREADS != 4)
              fprintf(out," -T%d",NTHREADS);
            for (k = 0; k < MTOP; k++)
              fprintf(out," -m%s",MASK[k]);
            if (useblock)
              if (usepath)
                fprintf(out," %s/%s.%d",pwd,root,i);
              else
                fprintf(out," %s.%d",root,i);
            else
              if (usepath)
                fprintf(out," %s/%s",pwd,root);
              else
                fprintf(out," %s",root);
            hgh = (i*j)/bits + 1;
            for (k = low; k < hgh; k++)
              if (useblock)
                if (usepath)
                  fprintf(out," %s/%s.%d",pwd,root,k);
                else
                  fprintf(out," %s.%d",root,k);
              else
                if (usepath)
                  fprintf(out," %s/%s",pwd,root);
                else
                  fprintf(out," %s",root);

            if (lblock == 1)   // ==> i = 1, [low,hgh) = [1,2)
              { fprintf(out," && mv");
                if (useblock)
                  fprintf(out," %s.1.%s.1.las",root,root);
                else
                  fprintf(out," %s.%s.las",root,root);
                if (usepath)
                  fprintf(out," %s/",pwd);
                else
                  fprintf(out," ");
                if (useblock)
                  fprintf(out,"%s.1.las",root);
                else
                  fprintf(out,"%s.las",root);
              }
            else if (DON)
              { fprintf(out," && mv");
                for (k = low; k < hgh; k++)
                  fprintf(out," %s.%d.%s.%d.las",root,i,root,k);
                fprintf(out," work%d",i);
                for (k = low; k < hgh; k++)
                  if (k != i)
                    fprintf(out," && mv %s.%d.%s.%d.las work%d",root,k,root,i,k);
              }

#ifdef LSF
            fprintf(out,"\"");
#endif
            fprintf(out,"\n");
            low = hgh;
          }
      }

    //  Check .las files (optional)

    if (ONAME != NULL)
      { fclose(out);
        sprintf(name,"%s.02.CHECK.OPT",ONAME);
        out = fopen(name,"w");
      }

    fprintf(out,"# Check initial .las files jobs (%d) (optional but recommended)\n",
                 (fblock-1) * ((lblock-fblock)/(BUNIT+1) + 1) +
                 (lblock-fblock+1) * ((lblock-1)/(BUNIT+1) + 1) );

#ifdef LSF
    jobid = 1;
#endif
    for (i = 1; i <= lblock; i++)
      for (j = (i < fblock ? fblock : 1); j <= lblock; )
        { k = j+BUNIT;
          if (k > lblock)
            k = lblock;
#ifdef LSF
          fprintf(out,LSF_CHECK,0,0,jobid++);
          fprintf(out," \"");
#endif
          fprintf(out,"LAcheck -vS");
          if (usepath)
            fprintf(out," %s/%s",pwd,root);
          else
            fprintf(out," %s",root);
          while (j <= k)
            { if (lblock == 1)
                { if (usepath)
                    if (useblock)
                      fprintf(out," %s/%s.1",pwd,root);
                    else
                      fprintf(out," %s/%s",pwd,root);
                  else
                    if (useblock)
                      fprintf(out," %s.1",root);
                    else
                      fprintf(out," %s",root);
                }
              else
                { if (DON)
                    fprintf(out," work%d/%s.%d.%s.%d",i,root,i,root,j);
                  else
                    fprintf(out," %s.%d.%s.%d",root,i,root,j);
                }
              j += 1;
            }
#ifdef LSF
          fprintf(out,"\"");
#endif
          fprintf(out,"\n");
        }

    if (ONAME != NULL)
      fclose(out);

    //  Higher level merges (if lblock > 1)

    if (lblock > 1)
      { int pow, stage;

        //  Determine the number of merging levels

        stage = 3;

        pow = 1;
        for (level = 0; pow < lblock; level++)
          pow *= DUNIT;

        //  Issue the commands for each merge level

        { int  p, cnt, dnt;

          cnt = lblock;
          dnt = (lblock-fblock)+1;
          for (i = 1; i <= level; i++)
            { int bits, dits;
              int low, hgh;

              if (ONAME != NULL)
                { sprintf(name,"%s.%02d.MERGE",ONAME,stage++);
                  out = fopen(name,"w");
                }

              bits = (cnt-1)/DUNIT+1;
              dits = (dnt-1)/DUNIT+1;

              //  Incremental update merges

#ifdef LSF
              jobid = 1;
#endif
              if (dnt >= 1)
                { int last;

                  last = (dnt == 1 || i == level);
                  fprintf(out,"# Level %d merge jobs (%d)\n",
                              i,bits*((lblock-fblock)+1) + dits*(fblock-1));
                  for (j = 1; j < fblock; j++)
                    {
#ifdef LSF
                      fprintf(out,LSF_MERGE,i,i,jobid++);
                      fprintf(out," \"");
#endif
                      if (last)
                        { if (DON)
                            { if (usepath)
                                fprintf(out,"mv %s/%s.%d.las work%d/L%d.%d.0.las && ",
                                            pwd,root,j,j,i,j);
                              else
                                fprintf(out,"mv %s.%d.las work%d/L%d.%d.0.las && ",root,j,j,i,j);
                            }
                          else
                            { if (usepath)
                                fprintf(out,"mv %s/%s.%d.las L%d.%d.0.las && ",pwd,root,j,i,j);
                              else
                                fprintf(out,"mv %s.%d.las L%d.%d.0.las && ",root,j,i,j);
                            }
                        }
                      low = 1;
                      for (p = 1; p <= dits; p++)
                        { hgh = (dnt*p)/dits;
#ifdef LSF
                          if (p > 1)
                            { fprintf(out,LSF_MERGE,i,i,jobid++);
                              fprintf(out," \"");
                            }
#endif
                          fprintf(out,"LAmerge");
                          if (VON)
                            fprintf(out," -v");
                          if (CON)
                            fprintf(out," -a");
                          if (last)
                            if (DON)
                              if (usepath)
                                fprintf(out," %s/%s.%d work%d/L%d.%d.0",pwd,root,j,j,i,j);
                              else
                                fprintf(out," %s.%d work%d/L%d.%d.0",root,j,j,i,j);
                            else
                              if (usepath)
                                fprintf(out," %s/%s.%d L%d.%d.0",pwd,root,j,i,j);
                              else
                                fprintf(out," %s.%d L%d.%d.0",root,j,i,j);
                          else
                            if (DON)
                              fprintf(out," work%d/L%d.%d.%d",j,i+1,j,p);
                            else
                              fprintf(out," L%d.%d.%d",i+1,j,p);
                          for (k = low; k <= hgh; k++)
                            if (i == 1)
                              if (DON)
                                fprintf(out," work%d/%s.%d.%s.%d",j,root,j,root,k+(fblock-1));
                              else
                                fprintf(out," %s.%d.%s.%d",root,j,root,k+(fblock-1));
                            else
                              if (DON)
                                fprintf(out," work%d/L%d.%d.%d",j,i,j,k);
                              else
                                fprintf(out," L%d.%d.%d",i,j,k);
#ifdef LSF
                          fprintf(out,"\"");
#endif
                          fprintf(out,"\n");
                          low = hgh+1;
                        }
                    }
                }
              else
                fprintf(out,"# Level %d merge jobs (%d)\n",i,bits*((lblock-fblock)+1));

              //  New block merges

              for (j = fblock; j <= lblock; j++) 
                { low = 1;
                  for (p = 1; p <= bits; p++)
                    { hgh = (cnt*p)/bits;
#ifdef LSF
                      fprintf(out,LSF_MERGE,i,i,jobid++);
                      fprintf(out," \"");
#endif
                      fprintf(out,"LAmerge");
                      if (VON)
                        fprintf(out," -v");
                      if (CON)
                        fprintf(out," -a");
                      if (i == level)
                        if (usepath)
                          fprintf(out," %s/%s.%d",pwd,root,j);
                        else
                          fprintf(out," %s.%d",root,j);
                      else
                        if (DON)
                          fprintf(out," work%d/L%d.%d.%d",j,i+1,j,p);
                        else
                          fprintf(out," L%d.%d.%d",i+1,j,p);
                      for (k = low; k <= hgh; k++)
                        if (i == 1)
                          if (DON)
                            fprintf(out," work%d/%s.%d.%s.%d",j,root,j,root,k);
                          else
                            fprintf(out," %s.%d.%s.%d",root,j,root,k);
                        else
                          if (DON)
                            fprintf(out," work%d/L%d.%d.%d",j,i,j,k);
                          else
                            fprintf(out," L%d.%d.%d",i,j,k);
#ifdef LSF
                      fprintf(out,"\"");
#endif
                      fprintf(out,"\n");
                      low = hgh+1;
                    }
                }

              //  Check new .las (optional)

              if (ONAME != NULL)
                { fclose(out);
                  sprintf(name,"%s.%02d.CHECK.OPT",ONAME,stage++);
                  out = fopen(name,"w");
                }

              fprintf(out,"# Check level %d .las files jobs (%d) (optional but recommended)\n",
                          i+1,(fblock-1)*((dits-1)/(BUNIT+1)+1) +
                              (lblock-fblock+1)*((bits-1)/(BUNIT+1)+1) );

#ifdef LSF
              jobid = 1;
#endif
              if (dnt >= 1)
                { int last;

                  last = (dnt == 1 || i == level);
                  for (j = 1; j < fblock; j++)
                    for (p = 1; p <= dits;)
                      { k = p+BUNIT;
                        if (k > dits)
                          k = dits;
#ifdef LSF
                        fprintf(out,LSF_CHECK,i,i,jobid++);
                        fprintf(out," \"");
#endif
                        fprintf(out,"LAcheck -vS");
                        if (usepath)
                          fprintf(out," %s/%s",pwd,root);
                        else
                          fprintf(out," %s",root);
                        while (p <= k)
                          { if (last)
                              if (usepath)
                                fprintf(out," %s/%s.%d",pwd,root,j);
                              else
                                fprintf(out," %s.%d",root,j);
                            else
                              if (DON)
                                fprintf(out," work%d/L%d.%d.%d",j,i+1,j,p);
                              else
                                fprintf(out," L%d.%d.%d",i+1,j,p);
                            p += 1;
                          }
#ifdef LSF
                        fprintf(out,"\"");
#endif
                        fprintf(out,"\n");
                      }
                }

              for (j = fblock; j <= lblock; j++) 
                for (p = 1; p <= bits;)
                  { k = p+BUNIT;
                    if (k > bits)
                      k = bits;
#ifdef LSF
                    fprintf(out,LSF_CHECK,i,i,jobid++);
                    fprintf(out," \"");
#endif
                    fprintf(out,"LAcheck -vS");
                    if (usepath)
                      fprintf(out," %s/%s",pwd,root);
                    else
                      fprintf(out," %s",root);
                    while (p <= k)
                      { if (i == level)
                          if (usepath)
                            fprintf(out," %s/%s.%d",pwd,root,j);
                          else
                            fprintf(out," %s.%d",root,j);
                        else
                          if (DON)
                            fprintf(out," work%d/L%d.%d.%d",j,i+1,j,p);
                          else
                            fprintf(out," L%d.%d.%d",i+1,j,p);
                        p += 1;
                      }
#ifdef LSF
                    fprintf(out,"\"");
#endif
                    fprintf(out,"\n");
                  }

              //  Cleanup (optional)

              if (ONAME != NULL)
                { fclose(out);
                  if (i == 1)
                    sprintf(name,"%s.%02d.RM.OPT",ONAME,stage++);
                  else
                    sprintf(name,"%s.%02d.RM",ONAME,stage++);
                  out = fopen(name,"w");
                }
              if (i == 1)
                fprintf(out,"# Remove level %d .las files (optional)\n",i);
              else
                fprintf(out,"# Remove level %d .las files\n",i);

              if (dnt >= 1)
                { int last;

                  last = (dnt == 1 || i == level);
                  for (j = 1; j < fblock; j++)
                    { low = 1;
                      for (p = 1; p <= dits; p++)
                        { hgh = (dnt*p)/dits;
                          if (DON)
                            fprintf(out,"cd work%d; ",j);
                          fprintf(out,"rm");
                          if (last)
                            fprintf(out," L%d.%d.0.las",i,j);
                          for (k = low; k <= hgh; k++)
                            if (i == 1)
                              fprintf(out," %s.%d.%s.%d.las",root,j,root,k+(fblock-1));
                            else
                              fprintf(out," L%d.%d.%d.las",i,j,k);
                          if (DON)
                            fprintf(out,"; cd ..");
                          fprintf(out,"\n");
                          low = hgh+1;
                        }
                    }
                }

              for (j = fblock; j <= lblock; j++) 
                { low = 1;
                  for (p = 1; p <= bits; p++)
                    { hgh = (cnt*p)/bits;
                      if (DON)
                        fprintf(out,"cd work%d; ",j);
                      fprintf(out,"rm");
                      for (k = low; k <= hgh; k++)
                        if (i == 1)
                          fprintf(out," %s.%d.%s.%d.las",root,j,root,k);
                        else
                          fprintf(out," L%d.%d.%d.las",i,j,k);
                      if (DON)
                        fprintf(out,"; cd ..");
                      fprintf(out,"\n");
                      low = hgh+1;
                    }
                }

              if (ONAME != NULL)
                fclose(out);

              if (dnt >= 1)
                { if (dnt > 1)
                    dnt = dits;
                  else
                    dnt = 0;
                }
              cnt = bits;
            }
        }
    }
  }

  free(root);
  free(pwd);
}

/*********************************************************************************************\
 *
 *  Produce a script to compute overlaps for all block pairs between two DBs, and then sort
 *    and merge them into as many .las files as their are blocks of the 1st DB.
 *
 *  Author:  Gene Myers
 *  Date  :  December 31, 2014
 *
 *********************************************************************************************/
 
#define LSF_MALIGN "bsub -q medium -n 4 -o MAPALL.out -e MAPALL.err -R span[hosts=1] -J align#%d"
#define LSF_MSORT  "bsub -q short -n 12 -o SORT.ALL.out -e SORT.ALL.err -R span[hosts=1] -J sort#%d"
#define LSF_MMERGE \
            "bsub -q short -n 12 -o MERGE%d.ALL.out -e MERGE%d.ALL.err -R span[hosts=1] -J merge#%d"

void mapper_script(int argc, char *argv[])
{ int   nblocks1, nblocks2;
  int   useblock1, useblock2;
  int   usepath1, usepath2;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  FILE *out;
  char  name[100];
  char *pwd1, *root1;
  char *pwd2, *root2;

  //  Make sure DAM and DB exist and the DB is partitioned, get number of blocks in partition

  pwd1   = PathTo(argv[1]);
  if (strcmp(argv[1]+(strlen(argv[1])-4),".dam") == 0)
    root1 = Root(argv[1],".dam");
  else
    root1 = Root(argv[1],".db");

  { int i, nfiles;
    FILE *dbvis;

    dbvis = fopen(Catenate(pwd1,"/",root1,".dam"),"r");
    if (dbvis == NULL)
      { dbvis = Fopen(Catenate(pwd1,"/",root1,".db"),"r");
        if (dbvis == NULL)
          exit (1);
      }

    if (fscanf(dbvis,"files = %d\n",&nfiles) != 1)
      SYSTEM_READ_ERROR
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];

        if (fgets(buffer,30000,dbvis) == NULL)
          SYSTEM_READ_ERROR
      }

    useblock1 = 1;
    if (fscanf(dbvis,"blocks = %d\n",&nblocks1) != 1 || nblocks1 == 1)
      { useblock1 = 0;
        nblocks1  = 1;
      }

    usepath1 = (strcmp(pwd1,".") != 0);

    fclose(dbvis);
  }

  pwd2   = PathTo(argv[2]);
  if (strcmp(argv[2]+(strlen(argv[2])-4),".dam") == 0)
    root2 = Root(argv[2],".dam");
  else
    root2 = Root(argv[2],".db");

  if (strcmp(root2,root1) == 0 && strcmp(pwd1,pwd2) == 0)
    { fprintf(stderr,"%s: Comparing the same data base %s/%s against itself, use HPCdaligner\n",
                     Prog_Name,pwd1,root1);
      exit (1);
    }

  { int i, nfiles;
    FILE *dbvis;

    dbvis = fopen(Catenate(pwd2,"/",root2,".dam"),"r");
    if (dbvis == NULL)
      { dbvis = Fopen(Catenate(pwd2,"/",root2,".db"),"r");
        if (dbvis == NULL)
          exit (1);
      }

    if (fscanf(dbvis,"files = %d\n",&nfiles) != 1)
      SYSTEM_READ_ERROR
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];

        if (fgets(buffer,30000,dbvis) == NULL)
          SYSTEM_READ_ERROR
      }

    useblock2 = 1;
    if (fscanf(dbvis,"blocks = %d\n",&nblocks2) != 1 || nblocks2 == 1)
      { useblock2 = 0;
        nblocks2  = 1;
      }

    usepath2 = (strcmp(pwd2,".") != 0);

    fclose(dbvis);
  }

  //  Set range fblock-lblock checking that DB.<fblock-1>.las exists & DB.<fblock>.las does not

  { char *eptr, *fptr, *src2;
    FILE *file;

    if (argc == 4)
      { fblock = strtol(argv[3],&eptr,10);
        if ((*eptr != '\0' && *eptr != '-') || eptr <= argv[3])
          { fprintf(stderr,"%s: final argument '%s' does not start with an integer\n",
                           Prog_Name,argv[3]);
            exit (1);
          }
        useblock2 = 1;
        if (*eptr == '-')
          { lblock = strtol(eptr+1,&fptr,10);
            if (*fptr != '\0' || fptr <= eptr+1)
              { fprintf(stderr,"%s: second part of range '%s' is not an integer\n",
                               Prog_Name,eptr+1);
                exit (1);
              }
          }
        else
          lblock = fblock;
        if (fblock < 1 || lblock > nblocks2 || fblock > lblock)
          { fprintf(stderr,"%s: range %d-%d is empty or out of bounds\n",Prog_Name,fblock,lblock);
            exit (1);
          }
      }
    else
      { fblock = 1;
        lblock = nblocks2;
      }

    if (usepath2)
      src2 = Strdup(Catenate(pwd2,"/",root2,""),"Allocating small string!");
    else
      src2 = Strdup(root2,"Allocating small string!");
    if (src2 == NULL)
      exit (1);

    if (fblock > 1)
      { file = fopen(Catenate(src2,".",root1,Numbered_Suffix(".",fblock-1,".las")),"r");
        if (file == NULL)
          { fprintf(stderr,"%s: File %s.%d.%s.las should already be present!\n",
                           Prog_Name,src2,fblock-1,root1);
            exit (1);
          }
        else
          fclose(file);
      }
    if (useblock2)
      { file = fopen(Catenate(src2,".",root1,Numbered_Suffix(".",fblock,".las")),"r");
        if (file != NULL)
          { fprintf(stderr,"%s: File %s.%d.%s.las should not yet exist!\n",
                           Prog_Name,src2,fblock,root1);
            exit (1);
          }
      }
    else
      { file = fopen(Catenate(src2,".",root1,".las"),"r");
        if (file != NULL)
          { fprintf(stderr,"%s: File %s.%s.las should not yet exist!\n",
                           Prog_Name,src2,root1);
            exit (1);
          }
      }

    free(src2);

    DON = (DON && (nblocks1 > 1));
    out = stdout;
  }

  { int level, njobs;
    int i, j, k;

    //  Create all work subdirectories if DON

    if (DON)
      { if (ONAME != NULL)
          { sprintf(name,"%s.00.MKDIR",ONAME);
            out = fopen(name,"w");
          }

        fprintf(out,"# Create work subdirectories\n");
        for (i = fblock; i <= lblock; i++)
          fprintf(out,"mkdir work%d\n",i);

        if (ONAME != NULL)
          fclose(out);
      }

    //  Produce all necessary daligner jobs ...

    if (ONAME != NULL)
      { sprintf(name,"%s.01.CMP",ONAME);
        out = fopen(name,"w");
      }

    njobs = nblocks1 * ( (lblock-fblock)/BUNIT + 1);

    fprintf(out,"# Daligner jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      { int bits;
        int low, hgh;

        bits = (nblocks1-1)/BUNIT+1;
        low  = 1;
        for (j = 1; j <= bits; j++)
          {
#ifdef LSF
            fprintf(out,LSF_MALIGN,jobid++);
            fprintf(out," \"");
#endif
            fprintf(out,"daligner -A");
            if (VON)
              fprintf(out," -v");
            if (BON)
              fprintf(out," -b");
            fprintf(out," -k%d",KINT);
            if (WINT != 6)
              fprintf(out," -w%d",WINT);
            fprintf(out," -h%d",HINT);
            if (TINT > 0)
              fprintf(out," -t%d",TINT);
            if (EREL > 0.)
              fprintf(out," -e%g",EREL);
            else
              fprintf(out," -e.85");
            if (LINT != 1000)
              fprintf(out," -l%d",LINT);
            if (SINT != 100)
              fprintf(out," -s%d",SINT);
            if (NTHREADS != 4)
              fprintf(out," -T%d",NTHREADS);
            if (MINT >= 0)
              fprintf(out," -M%d",MINT);
            if (PDIR != NULL)
              fprintf(out," -P%s",PDIR);
            for (k = 0; k < MTOP; k++)
              fprintf(out," -m%s",MASK[k]);

            fprintf(out," ");
            if (usepath2)
              fprintf(out,"%s/",pwd2);
            fprintf(out,"%s",root2);
            if (useblock2)
              fprintf(out,".%d",i);

	    hgh = 1 + (nblocks1*j)/bits;
            for (k = low; k < hgh; k++)
              { fprintf(out," ");
                if (usepath1)
                  fprintf(out,"%s/",pwd1);
                fprintf(out,"%s",root1);
                if (useblock1)
                  fprintf(out,".%d",k);
              }

            if (nblocks1 == 1)
              { if (useblock1 || usepath2)
                  { fprintf(out," && mv %s",root2);
                    if (useblock2)
                      fprintf(out,".%d.las",i);
                    if (useblock1)
                      fprintf(out,".%s.1.las ",root1);
                    else
                      fprintf(out,".%s.las ",root1);
                    if (useblock1)
                      { if (usepath2)
                          fprintf(out,"%s/",pwd2);
                        fprintf(out,"%s",root2);
                        if (useblock2)
                          fprintf(out,".%d",i);
                        fprintf(out,".%s.las",root1);
                      }
                    else
                      fprintf(out,"%s",pwd2);
                  }
              }
            else if (DON)
              { fprintf(out," && mv");
                for (k = low; k < hgh; k++)
                  { fprintf(out," %s",root2);
                    if (useblock2)
                      fprintf(out,".%d",i);
                    fprintf(out,".%s.%d.las",root1,k);
                  }
                fprintf(out," work%d",i);
              }
#ifdef LSF
            fprintf(out,"\"");
#endif
            fprintf(out,"\n");
            low = hgh;
          }
      }

    //  Check .las files (optional)

    if (ONAME != NULL)
      { fclose(out);
        sprintf(name,"%s.02.CHECK.OPT",ONAME);
        out = fopen(name,"w");
      }

    fprintf(out,"# Check initial .las files jobs (%d) (optional but recommended)\n",
                 (lblock-fblock+1) * ((nblocks1-1)/(BUNIT+1) + 1) );

#ifdef LSF
    jobid = 1;
#endif
    for (j = fblock; j <= lblock; j++)
      for (i = 1; i <= nblocks1; )
        { k = i+BUNIT;
          if (k > nblocks1)
            k = nblocks1;
#ifdef LSF
          fprintf(out,LSF_CHECK,0,0,jobid++);
          fprintf(out," \"");
#endif
          fprintf(out,"LAcheck -vS");
          if (usepath2)
            fprintf(out," %s/%s",pwd2,root2);
          else
            fprintf(out," %s",root2);
          if (usepath1)
            fprintf(out," %s/%s",pwd1,root1);
          else
            fprintf(out," %s",root1);
          while (i <= k)
            { fprintf(out," ");
              if (nblocks1 == 1)
                { if (usepath2)
                    fprintf(out,"%s/",pwd2);
                  fprintf(out,"%s",root2);
                  if (useblock2)
                    fprintf(out,".%d",j);
                  fprintf(out,".%s",root1);
                }
              else
                { if (DON)
                    fprintf(out,"work%d/",j);
                  fprintf(out,"%s",root2);
                  if (useblock2)
                    fprintf(out,".%d",j);
                  fprintf(out,".%s.%d",root1,i);
                }
              i += 1;
            }
#ifdef LSF
          fprintf(out,"\"");
#endif
          fprintf(out,"\n");
        }

    if (ONAME != NULL)
      fclose(out);

    //  Higher level merges (if lblock > 1)

    if (nblocks1 > 1)
      { int pow, stage;

        //  Determine the number of merging levels

        stage = 3;

        pow = 1;
        for (level = 0; pow < nblocks1; level++)
          pow *= DUNIT;

        //  Issue the commands for each merge level

        { int  p, cnt;

          cnt = nblocks1;
          for (i = 1; i <= level; i++)
            { int bits;
              int low, hgh;

              if (ONAME != NULL)
                { sprintf(name,"%s.%02d.MERGE",ONAME,stage++);
                  out = fopen(name,"w");
                }

              bits = (cnt-1)/DUNIT+1;
              fprintf(out,"# Level %d jobs (%d)\n",i,bits*((lblock-fblock)+1));

              //  Block merges

#ifdef LSF
              jobid = 1;
#endif
              for (j = fblock; j <= lblock; j++) 
                { low = 1;
                  for (p = 1; p <= bits; p++)
                    { hgh = (cnt*p)/bits;
#ifdef LSF
                      fprintf(out,LSF_MMERGE,i,i,jobid++);
                      fprintf(out," \"");
#endif
                      fprintf(out,"LAmerge ");
                      if (VON)
                        fprintf(out,"-v ");
                      if (CON)
                        fprintf(out,"-a ");
                      if (i == level)
                        { if (usepath2)
                            fprintf(out,"%s/",pwd2);
                          fprintf(out,"%s",root2);
                          if (useblock2)
                            fprintf(out,".%d",j);
                          fprintf(out,".%s",root1);
                        }
                      else
                        { if (DON)
                            fprintf(out,"work%d/",j);
                          fprintf(out,"L%d.%d.%d",i+1,j,p);
                        }
                      for (k = low; k <= hgh; k++)
                        if (i == 1)
                          { if (DON)
                              fprintf(out," work%d/",j);
                            else
                              fprintf(out," ");
                            fprintf(out,"%s",root2);
                            if (useblock2)
                              fprintf(out,".%d",j);
                            fprintf(out,".%s.%d",root1,k);
                          }
                        else
                          if (DON)
                            fprintf(out," work%d/L%d.%d.%d",j,i,j,k);
                          else
                            fprintf(out," L%d.%d.%d",i,j,k);

#ifdef LSF
                      fprintf(out,"\"");
#endif
		      fprintf(out,"\n");
                      low = hgh+1;
                    }
                }

              //  Check new .las (optional)

              if (ONAME != NULL)
                { fclose(out);
                  sprintf(name,"%s.%02d.CHECK.OPT",ONAME,stage++);
                  out = fopen(name,"w");
                }

              fprintf(out,"# Check level %d .las files jobs (%d) (optional but recommended)\n",
                          i+1,(lblock-fblock+1)*((bits-1)/(BUNIT+1)+1) );

#ifdef LSF
              jobid = 1;
#endif
              for (j = fblock; j <= lblock; j++) 
                for (p = 1; p <= bits; )
                  { k = p+BUNIT;
                    if (k > bits)
                      k = bits;
#ifdef LSF
                    fprintf(out,LSF_CHECK,0,0,jobid++);
                    fprintf(out," \"");
#endif
                    fprintf(out,"LAcheck -vS");
                    if (usepath2)
                      fprintf(out," %s/%s",pwd2,root2);
                    else
                      fprintf(out," %s",root2);
                    if (usepath1)
                      fprintf(out," %s/%s",pwd1,root1);
                    else
                      fprintf(out," %s",root1);
                    while (p <= k)
                      { fprintf(out," ");
                        if (i == level)
                          { if (usepath2)
                              fprintf(out,"%s/",pwd2);
                            fprintf(out,"%s",root2);
                            if (useblock2)
                              fprintf(out,".%d",j);
                            fprintf(out,".%s",root1);
                          }
                        else
                          { if (DON)
                              fprintf(out,"work%d/",j);
                            fprintf(out,"L%d.%d.%d",i+1,j,p);
                          }
                        p += 1;
                      }
#ifdef LSF
                    fprintf(out,"\"");
#endif
                    fprintf(out,"\n");
                  }

              //  Cleanup (optional)

              if (ONAME != NULL)
                { fclose(out);
                  sprintf(name,"%s.%02d.RM",ONAME,stage++);
                  out = fopen(name,"w");
                }

              fprintf(out,"# Remove level %d .las files\n",i);

              for (j = fblock; j <= lblock; j++) 
                { low = 1;
                  for (p = 1; p <= bits; p++)
                    { hgh = (cnt*p)/bits;
                      if (DON)
                        fprintf(out,"cd work%d; ",j);
                      fprintf(out,"rm");
                      for (k = low; k <= hgh; k++)
                        if (i == 1)
                          { fprintf(out," %s",root2);
                            if (useblock2)
                              fprintf(out,".%d",j);
                            fprintf(out,".%s.%d.las",root1,k);
                          }
                        else
                          fprintf(out," L%d.%d.%d.las",i,j,k);
                      if (DON)
                        fprintf(out,"; cd ..");
                      fprintf(out,"\n");
                      low = hgh+1;
                    }
                }

              if (ONAME != NULL)
                fclose(out);

              cnt  = bits;
            }
        }
    }
  }

  free(root2);
  free(pwd2);
  free(root1);
  free(pwd1);

  exit (0);
}

int main(int argc, char *argv[])
{ int    i, j, k;
  int    flags[128];
  char  *eptr;
  int    mapper;

  //  Process options and decide if its a overlap or mapper script

  ARG_INIT("HPC.daligner")

  KINT  = 0;
  HINT  = 0;
  HGAP  = 0;
  EREL  = 0.;

  BUNIT = 4;
  DUNIT = 250;
  TINT  = 0;
  WINT  = 6;
  LINT  = 1000;
  SINT  = 100;
  MINT  = -1;
  PDIR  = NULL;

  MTOP = 0;
  MMAX = 10;
  MASK = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
  if (MASK == NULL)
    exit (1);
  ONAME = NULL;

  NTHREADS = 4;

  j = 1;
  for (i = 1; i < argc; i++)
    if (argv[i][0] == '-')
      switch (argv[i][1])
      { default:
          ARG_FLAGS("vbadAI");
          break;
        case 'e':
          ARG_REAL(EREL)
          if (EREL < .7 || EREL >= 1.)
            { fprintf(stderr,"%s: Average correlation must be in [.7,1.) (%g)\n",Prog_Name,EREL);
              exit (1);
            }
          break;
        case 'f':
          ONAME = argv[i]+2;
          break;
        case 'h':
          ARG_POSITIVE(HINT,"Hit threshold (in bp.s)")
          break;
        case 'k':
          ARG_POSITIVE(KINT,"K-mer length")
          if (KINT > 32)
            { fprintf(stderr,"%s: K-mer length must be 32 or less\n",Prog_Name);
              exit (1);
            }
          break;
        case 'l':
          ARG_POSITIVE(LINT,"Minimum ovlerap length")
          break;
        case 'm':
          if (MTOP >= MMAX)
            { MMAX = 1.2*MTOP + 10;
              MASK = (char **) Realloc(MASK,MMAX*sizeof(char *),"Reallocating mask track array");
              if (MASK == NULL)
                exit (1);
            }
          MASK[MTOP++] = argv[i]+2;
          break;
        case 's':
          ARG_POSITIVE(SINT,"Trace spacing")
          break;
        case 't':
          ARG_POSITIVE(TINT,"Tuple suppression frequency")
          break;
        case 'w':
          ARG_POSITIVE(WINT,"Log of bin width")
          break;
        case 'B':
          ARG_NON_NEGATIVE(BUNIT,"Blocks per command")
          break;
        case 'D':
          ARG_NON_NEGATIVE(DUNIT,"File per merge")
          if (DUNIT < 3)
            { fprintf(stderr,"%s: Files per merge must be at least 3 (%d)\n",
                             Prog_Name,DUNIT);
              exit (1);
            }
          break;
        case 'H':
          ARG_POSITIVE(HGAP,"HGAP threshold (in bp.s)")
          break;
        case 'M':
          ARG_NON_NEGATIVE(MINT,"Memory allocation (in Gb)")
          break;
        case 'P':
          PDIR = argv[i]+2;
          break;
        case 'T':
          ARG_POSITIVE(NTHREADS,"Number of threads")
          break;
      }
    else
      argv[j++] = argv[i];
  argc = j;

  VON = flags['v'];
  BON = flags['b'];
  CON = flags['a'];
  DON = flags['d'];

  if (argc < 2 || argc > 4)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[4]);
      exit (1);
    }

  if (argc == 2)
    mapper = 0;
  else if (argc == 4)
    mapper = 1;
  else
    { (void) strtol(argv[2],&eptr,10);
      if ((*eptr == '\0' || *eptr == '-') && eptr > argv[2])
        mapper = 0;
      else
        mapper = 1;
    }

  if (mapper)
    { if (HGAP > 0)
        { fprintf(stderr,"%s: Cannot use -H option in a comparison script\n",Prog_Name);
          exit (1);
        }
      if (KINT <= 0)
        KINT = 20;
      if (HINT <= 0)
        HINT = 50;
      if (EREL <= 0.)
        EREL = .85;
    }
  else
    { if (KINT <= 0)
        KINT = 14;
      if (HINT <= 0)
        HINT = 35;
    }

  for (j = 1; 2*j <= NTHREADS; j *= 2)
    ;
  NTHREADS = j;

  if (mapper)
    mapper_script(argc,argv);
  else
    daligner_script(argc,argv);

  exit (0);
}
