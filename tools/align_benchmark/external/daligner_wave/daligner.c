/*********************************************************************************************\
 *
 *  Find all local alignment between long, noisy DNA reads:
 *    Compare sequences in 'subject' database against those in the list of 'target' databases
 *    searching for local alignments of 1000bp or more (defined constant MIN_OVERLAP in
 *    filter.c).  Subject is compared in both orientations againt each target.  An output
 *    stream of 'Overlap' records (see align.h) is written in binary to the standard output,
 *    each encoding a given found local alignment between two of the sequences.  The -v
 *    option turns on a verbose reporting mode that gives statistics on each major stage.
 *
 *    The filter operates by looking for a pair of diagonal bands of width 2^'s' that contain
 *    a collection of exact matching 'k'-mers between the two sequences, such that the total
 *    number of bases covered by 'k'-mer hits is 'h'.  k cannot be larger than 32 in the
 *    current implementation.
 *
 *    Some k-mers are significantly over-represented (e.g. homopolymer runs).  These are
 *    suppressed as seed hits, with the parameter 't' -- any k-mer that occurs more than
 *    't' times in either the subject or target is not counted as a seed hit.  If the -t
 *    option is absent then no k-mer is suppressed.  Alternatively, the option -M specifies
 *    that 't' is dynamically set to the largest value such that less than -M memory is
 *    used.
 *
 *    For each subject, target pair, say XXX and YYY, the program outputs a file containing
 *    overlaps of the form XXX.YYY.[C|N]#.las where C implies that the reads in XXX were
 *    complemented and N implies they were not (both comparisons are performed), and # is
 *    the thread that detected and wrote out the collection of overlaps.  For example, if
 *    NTHREAD in the program is 4, then 8 files are output for each subject, target pair.
 *
 *  Author:  Gene Myers
 *  Date  :  June 1, 2014
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include <sys/param.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#include "DB.h"
#include "filter.h"

static char *Usage[] =
  { "[-vbAI] [-k<int(14)>] [-w<int(6)>] [-h<int(35)>] [-t<int>] [-M<int>] [-P<dir(/tmp)>]",
    "        [-e<double(.70)] [-l<int(1000)>] [-s<int(100)>] [-H<int>] [-T<int(4)>]",
    "        [-m<track>]+ <subject:db|dam> <target:db|dam> ...",
  };

int     VERBOSE;   //   Globally visible to filter.c
char   *SORT_PATH;
int     BIASED;
int     MINOVER;
int     HGAP_MIN;
int     SYMMETRIC;
int     IDENTITY;
uint64  MEM_LIMIT;
uint64  MEM_PHYSICAL;

/*  Adapted from code by David Robert Nadeau (http://NadeauSoftware.com) licensed under
 *     "Creative Commons Attribution 3.0 Unported License"
 *          (http://creativecommons.org/licenses/by/3.0/deed.en_US)
 *
 *   I removed Windows options, reformated, and return int64 instead of size_t
 */

static int64 getMemorySize( )
{
#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))

  // OSX, NetBSD, OpenBSD

  int     mib[2];
  size_t  size = 0;
  size_t  len = sizeof( size );

  mib[0] = CTL_HW;
#if defined(HW_MEMSIZE)
  mib[1] = HW_MEMSIZE;            // OSX
#elif defined(HW_PHYSMEM64)
  mib[1] = HW_PHYSMEM64;          // NetBSD, OpenBSD
#endif
  if (sysctl(mib,2,&size,&len,NULL,0) == 0)
    return ((size_t) size);
  return (0);

#elif defined(_SC_AIX_REALMEM)

  // AIX

  return ((size_t) sysconf( _SC_AIX_REALMEM ) * ((size_t) 1024L));

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  // FreeBSD, Linux, OpenBSD, & Solaris

  size_t  size = 0;

  size = (size_t) sysconf(_SC_PHYS_PAGES);
  return (size * ((size_t) sysconf(_SC_PAGESIZE)));

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)

  // ? Legacy ?

  size_t  size = 0;

  size = (size_t) sysconf(_SC_PHYS_PAGES);
  return (size * ((size_t) sysconf(_SC_PAGE_SIZE)));

#elif defined(CTL_HW) && (defined(HW_PHYSMEM) || defined(HW_REALMEM))

  // DragonFly BSD, FreeBSD, NetBSD, OpenBSD, and OSX

  int          mib[2];
  unsigned int size = 0;
  size_t       len  = sizeof( size );

  mib[0] = CTL_HW;
#if defined(HW_REALMEM)
  mib[1] = HW_REALMEM;		// FreeBSD
#elif defined(HW_PYSMEM)
  mib[1] = HW_PHYSMEM;		// Others
#endif
  if (sysctl(mib,2,&size,&len,NULL,0) == 0)
    return (size_t)size;
  return (0);

#else

  return (0);

#endif
}

typedef struct
  { int *ano;
    int *end;
    int  idx;
    int  out;
  } Event;

static void reheap(int s, Event **heap, int hsize)
{ int      c, l, r;
  Event   *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      hr = heap[r];
      if (hr->idx > hl->idx)
        { if (hs->idx > hl->idx)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { if (hs->idx > hr->idx)
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

static int64 merge_size(DAZZ_DB *block, int mtop)
{ Event       ev[mtop+1];
  Event      *heap[mtop+2];
  int         r, mhalf;
  int64       nsize;

  { DAZZ_TRACK *track;
    int         i;

    track = block->tracks;
    for (i = 0; i < mtop; i++)
      { ev[i].ano = ((int *) (track->data)) + ((int64 *) (track->anno))[0];
        ev[i].out = 1;
        heap[i+1] = ev+i;
        track = track->next;
      }
    ev[mtop].idx = INT32_MAX;
    heap[mtop+1] = ev+mtop;
  }

  mhalf = mtop/2;

  nsize = 0;
  for (r = 0; r < block->nreads; r++)
    { int         i, level, hsize;
      DAZZ_TRACK *track;

      track = block->tracks;
      for (i = 0; i < mtop; i++)
        { ev[i].end = ((int *) (track->data)) + ((int64 *) (track->anno))[r+1];
          if (ev[i].ano < ev[i].end)
            ev[i].idx = *(ev[i].ano);
          else
            ev[i].idx = INT32_MAX;
          track = track->next;
        }
      hsize = mtop;

      for (i = mhalf; i > 1; i--)
        reheap(i,heap,hsize);

      level = 0;
      while (1)
        { Event *p;

          reheap(1,heap,hsize);

          p = heap[1];
          if (p->idx == INT32_MAX) break;

          p->out = 1-p->out;
          if (p->out)
            { level -= 1;
              if (level == 0)
                nsize += 1;
            }
          else
            { if (level == 0)
                nsize += 1;
              level += 1;
            }
          p->ano += 1;
          if (p->ano >= p->end)
            p->idx = INT32_MAX;
          else
            p->idx = *(p->ano);
        }
    }

  return (nsize);
}

static DAZZ_TRACK *merge_tracks(DAZZ_DB *block, int mtop, int64 nsize)
{ DAZZ_TRACK *ntrack;
  Event       ev[mtop+1];
  Event      *heap[mtop+2];
  int         r, mhalf;
  int64      *anno;
  int        *data;

  ntrack = (DAZZ_TRACK *) Malloc(sizeof(DAZZ_TRACK),"Allocating merged track");
  if (ntrack == NULL)
    exit (1);
  ntrack->name = Strdup("merge","Allocating merged track");
  ntrack->anno = anno = (int64 *) Malloc(sizeof(int64)*(block->nreads+1),"Allocating merged track");
  ntrack->data = data = (int *) Malloc(sizeof(int)*nsize,"Allocating merged track");
  ntrack->size = sizeof(int);
  ntrack->next = NULL;
  if (anno == NULL || data == NULL || ntrack->name == NULL)
    exit (1);

  { DAZZ_TRACK *track;
    int         i;

    track = block->tracks;
    for (i = 0; i < mtop; i++)
      { ev[i].ano = ((int *) (track->data)) + ((int64 *) (track->anno))[0];
        ev[i].out = 1;
        heap[i+1] = ev+i;
        track = track->next;
      }
    ev[mtop].idx = INT32_MAX;
    heap[mtop+1] = ev+mtop;
  }

  mhalf = mtop/2;

  nsize = 0;
  for (r = 0; r < block->nreads; r++)
    { int         i, level, hsize;
      DAZZ_TRACK *track;

      anno[r] = nsize;

      track = block->tracks;
      for (i = 0; i < mtop; i++)
        { ev[i].end = ((int *) (track->data)) + ((int64 *) (track->anno))[r+1];
          if (ev[i].ano < ev[i].end)
            ev[i].idx = *(ev[i].ano);
          else
            ev[i].idx = INT32_MAX;
          track = track->next;
        }
      hsize = mtop;

      for (i = mhalf; i > 1; i--)
        reheap(i,heap,hsize);

      level = 0;
      while (1)
        { Event *p;

          reheap(1,heap,hsize);

          p = heap[1];
          if (p->idx == INT32_MAX) break;

          p->out = 1-p->out;
          if (p->out)
            { level -= 1;
              if (level == 0)
                data[nsize++] = p->idx;
            }
          else
            { if (level == 0)
                data[nsize++] = p->idx;
              level += 1;
            }
          p->ano += 1;
          if (p->ano >= p->end)
            p->idx = INT32_MAX;
          else
            p->idx = *(p->ano);
        }
    }
  anno[r] = nsize;

  return (ntrack);
}

static int read_DB(DAZZ_DB *block, char *name, char **mask, int *mstat, int mtop, int kmer)
{ int i, isdam, status, kind, stop;

  isdam = Open_DB(name,block);
  if (isdam < 0)
    exit (1);

  for (i = 0; i < mtop; i++)
    { status = Check_Track(block,mask[i],&kind);
      if (status >= 0)
        if (kind == MASK_TRACK)
          mstat[i] = 0;
        else
          { if (mstat[i] != 0)
              mstat[i] = -3;
          }
      else
        { if (mstat[i] == -2)
            mstat[i] = status;
        }
      if (status == 0 && kind == MASK_TRACK)
        Load_Track(block,mask[i]);
    }

  Trim_DB(block);

  stop = 0;
  for (i = 0; i < mtop; i++)
    { DAZZ_TRACK *track;
      int64      *anno;
      int         j;

      status = Check_Track(block,mask[i],&kind);
      if (status < 0 || kind != MASK_TRACK)
        continue;
      stop += 1;
      track = Load_Track(block,mask[i]);

      anno = (int64 *) (track->anno); 
      for (j = 0; j <= block->nreads; j++)
        anno[j] /= sizeof(int);
    }

  if (stop > 1)
    { int64       nsize;
      DAZZ_TRACK *track;

      nsize = merge_size(block,stop);
      track = merge_tracks(block,stop,nsize);

      while (block->tracks != NULL)
        Close_Track(block,block->tracks->name);

      block->tracks = track;
    }

  if (block->cutoff < kmer)
    { for (i = 0; i < block->nreads; i++)
        if (block->reads[i].rlen < kmer)
          { fprintf(stderr,"%s: Block %s contains reads < %dbp long !  Run DBsplit.\n",
                           Prog_Name,name,kmer);
            exit (1);
          }
    }

  Read_All_Sequences(block,0);

  return (isdam);
}

static void complement(char *s, int len)
{ char *t;
  int   c;

  t = s + (len-1);
  while (s < t)
    { c = *s;
      *s = (char) (3-*t);
      *t = (char) (3-c);
      s += 1;
      t -= 1;
    }
  if (s == t)
    *s = (char) (3-*s);
}

static DAZZ_DB *complement_DB(DAZZ_DB *block, int inplace)
{ static DAZZ_DB _cblock, *cblock = &_cblock;
  int            nreads;
  DAZZ_READ     *reads;
  char          *seq;
  
  nreads = block->nreads;
  reads  = block->reads;
  if (inplace)
    { seq = (char *) block->bases;
      cblock = block;
    }
  else
    { seq  = (char *) Malloc(block->reads[nreads].boff+1,"Allocating dazzler sequence block");
      if (seq == NULL)
        exit (1);
      *seq++ = 4;
      memmove(seq,block->bases,block->reads[nreads].boff);
      *cblock = *block;
      cblock->bases  = (void *) seq;
      cblock->tracks = NULL;
    }

  { int   i;
    float x;

    x = cblock->freq[0];
    cblock->freq[0] = cblock->freq[3];
    cblock->freq[3] = x;

    x = cblock->freq[1];
    cblock->freq[1] = cblock->freq[2];
    cblock->freq[2] = x;

    for (i = 0; i < nreads; i++)
      complement(seq+reads[i].boff,reads[i].rlen);
  }

  { DAZZ_TRACK *src, *trg;
    int        *data, *tata;
    int         i, x, rlen;
    int64      *tano, *anno;
    int64       j, k;

    for (src = block->tracks; src != NULL; src = src->next)
      { tano = (int64 *) src->anno;
        tata = (int   *) src->data;

        if (inplace)
          { data = tata;
            anno = tano;
            trg  = src;
          }
        else
          { data = (int *) Malloc(sizeof(int)*tano[nreads],
                                  "Allocating dazzler interval track data");
            anno = (int64 *) Malloc(sizeof(int64)*(nreads+1),
                                    "Allocating dazzler interval track index");
            trg  = (DAZZ_TRACK *) Malloc(sizeof(DAZZ_TRACK),
                                         "Allocating dazzler interval track header");
            if (data == NULL || trg == NULL || anno == NULL)
              exit (1);

            trg->name = Strdup(src->name,"Copying track name");
            if (trg->name == NULL)
              exit (1);

            trg->size = 4;
            trg->anno = (void *) anno;
            trg->data = (void *) data;
            trg->next = cblock->tracks;
            cblock->tracks = trg;
          }

        for (i = 0; i < nreads; i++)
          { rlen = reads[i].rlen;
            anno[i] = tano[i];
            j = tano[i+1]-1;
            k = tano[i];
            while (k < j)
              { x = tata[j];
                data[j--] = rlen - tata[k];
                data[k++] = rlen - x;
              }
            if (k == j)
              data[k] = rlen - tata[k];
          }
        anno[nreads] = tano[nreads];
      }
  }

  return (cblock);
}

static char *CommandBuffer(char *aname, char *bname)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  len = 2*(strlen(aname) + strlen(bname)) + 200;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((cat = (char *) realloc(cat,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making path name)\n",Prog_Name);
          exit (1);
        }
    }
  return (cat);
}

int main(int argc, char *argv[])
{ DAZZ_DB    _ablock, _bblock;
  DAZZ_DB    *ablock = &_ablock, *bblock = &_bblock;
  char       *afile,  *bfile;
  char       *aroot,  *broot;
  void       *aindex, *bindex;
  int         alen,    blen;
  Align_Spec *asettings;
  int         isdam;
  int         MMAX, MTOP, *MSTAT;
  char      **MASK;

  int    KMER_LEN;
  int    BIN_SHIFT;
  int    MAX_REPS;
  int    HIT_MIN;
  double AVE_ERROR;
  int    SPACING;
  int    NTHREADS;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;
    DIR   *dirp;

    ARG_INIT("daligner")

    KMER_LEN  = 14;
    HIT_MIN   = 35;
    BIN_SHIFT = 6;
    MAX_REPS  = 0;
    HGAP_MIN  = 0;
    AVE_ERROR = .70;
    SPACING   = 100;
    MINOVER   = 1000;    //   Globally visible to filter.c
    NTHREADS  = 4;
    SORT_PATH = "/tmp";

    MEM_PHYSICAL = getMemorySize();
    MEM_LIMIT    = MEM_PHYSICAL;
    if (MEM_PHYSICAL == 0)
      { fprintf(stderr,"\nWarning: Could not get physical memory size\n");
        fflush(stderr);
      }

    MTOP  = 0;
    MMAX  = 10;
    MASK  = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
    MSTAT = (int *) Malloc(MMAX*sizeof(int),"Allocating mask status array");
    if (MASK == NULL || MSTAT == NULL)
      exit (1);

    j    = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vbAI")
            break;
          case 'k':
            ARG_POSITIVE(KMER_LEN,"K-mer length")
            if (KMER_LEN > 32)
              { fprintf(stderr,"%s: K-mer length must be 32 or less\n",Prog_Name);
                exit (1);
              }
            break;
          case 'w':
            ARG_POSITIVE(BIN_SHIFT,"Log of bin width")
            break;
          case 'h':
            ARG_POSITIVE(HIT_MIN,"Hit threshold (in bp.s)")
            break;
          case 't':
            ARG_POSITIVE(MAX_REPS,"Tuple supression frequency")
            break;
          case 'H':
            ARG_POSITIVE(HGAP_MIN,"HGAP threshold (in bp.s)")
            break;
          case 'e':
            ARG_REAL(AVE_ERROR)
            if (AVE_ERROR < .7 || AVE_ERROR >= 1.)
              { fprintf(stderr,"%s: Average correlation must be in [.7,1.) (%g)\n",
                               Prog_Name,AVE_ERROR);
                exit (1);
              }
            break;
          case 'l':
            ARG_POSITIVE(MINOVER,"Minimum alignment length")
            break;
          case 's':
            ARG_POSITIVE(SPACING,"Trace spacing")
            break;
          case 'M':
            { int limit;

              ARG_NON_NEGATIVE(limit,"Memory allocation (in Gb)")
              MEM_LIMIT = limit * 0x40000000ll;
              break;
            }
          case 'm':
            if (MTOP >= MMAX)
              { MMAX  = 1.2*MTOP + 10;
                MASK  = (char **) Realloc(MASK,MMAX*sizeof(char *),"Reallocating mask track array");
                MSTAT = (int *) Realloc(MSTAT,MMAX*sizeof(int),"Reallocating mask status array");
                if (MASK == NULL || MSTAT == NULL)
                  exit (1);
              }
            MASK[MTOP++] = argv[i]+2;
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            if ((dirp = opendir(SORT_PATH)) == NULL)
              { fprintf(stderr,"%s: -P option: cannot open directory %s\n",Prog_Name,SORT_PATH);
                exit (1);
              }
            closedir(dirp);
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];   //  Globally declared in filter.h
    BIASED    = flags['b'];   //  Globally declared in filter.h
    SYMMETRIC = 1-flags['A'];
    IDENTITY  = flags['I'];

    if (argc <= 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        exit (1);
      }

    for (j = 0; j < MTOP; j++)
      MSTAT[j] = -2;
  }

  MINOVER *= 2;
  if (Set_Filter_Params(KMER_LEN,BIN_SHIFT,MAX_REPS,HIT_MIN,NTHREADS))
    { fprintf(stderr,"Illegal combination of filter parameters\n");
      exit (1);
    }

  /* Read in the reads in A */

  afile = argv[1];
  isdam = read_DB(ablock,afile,MASK,MSTAT,MTOP,KMER_LEN);
  if (isdam)
    aroot = Root(afile,".dam");
  else
    aroot = Root(afile,".db");

  asettings = New_Align_Spec( AVE_ERROR, SPACING, ablock->freq, 1);

  /* Compare against reads in B in both orientations */

  { int   i, j;
    char *command;

    aindex = NULL;
    broot  = NULL;
    for (i = 2; i < argc; i++)
      { bfile = argv[i];
        if (strcmp(afile,bfile) != 0)
          { isdam = read_DB(bblock,bfile,MASK,MSTAT,MTOP,KMER_LEN);
            if (isdam)
              broot = Root(bfile,".dam");
            else
              broot = Root(bfile,".db");
          }
        else
          broot = aroot;

        if (i == 2)
          { for (j = 0; j < MTOP; j++)
              { if (MSTAT[j] == -2)
                  printf("%s: Warning: -m%s option given but no track found.\n",Prog_Name,MASK[j]);
                else if (MSTAT[j] == -1)
                  printf("%s: Warning: %s track not sync'd with relevant db.\n",Prog_Name,MASK[j]);
                else if (MSTAT[j] == -3)
                  printf("%s: Warning: %s track is not a mask track.\n",Prog_Name,MASK[j]);
              }

            if (VERBOSE)
              printf("\nBuilding index for %s\n",aroot);
            aindex = Sort_Kmers(ablock,&alen);
          }

        if (aroot != broot)
          { if (VERBOSE)
              printf("\nBuilding index for %s\n",broot);
            bindex = Sort_Kmers(bblock,&blen);
            Match_Filter(aroot,ablock,broot,bblock,aindex,alen,bindex,blen,0,asettings);

            bblock = complement_DB(bblock,1);
            if (VERBOSE)
              printf("\nBuilding index for c(%s)\n",broot);
            bindex = Sort_Kmers(bblock,&blen);
            Match_Filter(aroot,ablock,broot,bblock,aindex,alen,bindex,blen,1,asettings);
          }
        else
          { Match_Filter(aroot,ablock,aroot,ablock,aindex,alen,aindex,alen,0,asettings);

            bblock = complement_DB(ablock,0);
            if (VERBOSE)
              printf("\nBuilding index for c(%s)\n",aroot);
            bindex = Sort_Kmers(bblock,&blen);
            Match_Filter(aroot,ablock,aroot,bblock,aindex,alen,bindex,blen,1,asettings);

            bblock->reads = NULL;  //  ablock & bblock share "reads" vector, don't let Close_DB
                                   //     free it !
          }

        Close_DB(bblock);

        command = CommandBuffer(aroot,broot);

        sprintf(command,"LAsort %s/%s.%s.[CN]*.las",SORT_PATH,aroot,broot);
        if (VERBOSE)
          printf("\n%s\n",command);
        system(command);
        sprintf(command,"LAmerge %s.%s.las %s/%s.%s.[CN]*.S.las",aroot,broot,SORT_PATH,aroot,broot);
        if (VERBOSE)
          printf("%s\n",command);
        system(command);
        sprintf(command,"rm %s/%s.%s.[CN]*.las",SORT_PATH,aroot,broot);
        if (VERBOSE)
          printf("%s\n",command);
        system(command);
        if (aroot != broot && SYMMETRIC)
          { sprintf(command,"LAsort %s/%s.%s.[CN]*.las",SORT_PATH,broot,aroot);
            if (VERBOSE)
              printf("%s\n",command);
            system(command);
            sprintf(command,"LAmerge %s.%s.las %s/%s.%s.[CN]*.S.las",broot,aroot,
                                                                     SORT_PATH,broot,aroot);
            if (VERBOSE)
              printf("%s\n",command);
            system(command);
            sprintf(command,"rm %s/%s.%s.[CN]*.las",SORT_PATH,broot,aroot);
            if (VERBOSE)
              printf("%s\n",command);
            system(command);
          }

        if (aroot != broot)
          free(broot);
      }
  }

  exit (0);
}
