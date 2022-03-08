/*******************************************************************************************
 *
 *  Compressed data base module.  Auxiliary routines to open and manipulate a data base for
 *    which the sequence and read information are separated into two separate files, and the
 *    sequence is compressed into 2-bits for each base.  Support for tracks of additional
 *    information, and trimming according to the current partition.  Eventually will also
 *    support compressed quality information. 
 *
 *  Author :  Gene Myers
 *  Date   :  July 2013
 *  Revised:  April 2014
 *
 ********************************************************************************************/

#ifndef _DAZZ_DB

#define _DAZZ_DB

#include <stdio.h>

#include "QV.h"

#define HIDE_FILES          //  Auxiliary DB files start with a . so they are "hidden"
                            //    Undefine if you don't want this

//  For interactive applications where it is inappropriate to simply exit with an error
//    message to standard error, define the constant INTERACTIVE.  If set, then error
//    messages are put in the global variable Ebuffer and the caller of a DB routine
//    can decide how to deal with the error.
//
//  DB, QV, or alignment routines that can encounter errors function as before in
//    non-INTERACTIVE mode by exiting after printing an error message to stderr.  In
//    INTERACTIVE mode the routines place a message at EPLACE and return an error
//    value.  For such routines that were previously void, they are now int, and
//    return 1 if an error occured, 0 otherwise.

#ifdef INTERACTIVE

#define EPRINTF sprintf
#define EPLACE  Ebuffer
#define EXIT(x) return (x)

#else // BATCH

#define EPRINTF fprintf
#define EPLACE  stderr
#define EXIT(x) exit (1)

#endif

typedef unsigned char      uint8;
typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;
typedef signed char        int8;
typedef signed short       int16;
typedef signed int         int32;
typedef signed long long   int64;
typedef float              float32;
typedef double             float64;


/*******************************************************************************************
 *
 *  COMMAND LINE INTERPRETATION MACROS
 *
 ********************************************************************************************/

extern char *Prog_Name;   //  Name of program

#ifdef INTERACTIVE

extern char Ebuffer[];

#endif

#define ARG_INIT(name)                  \
  Prog_Name = Strdup(name,"");          \
  for (i = 0; i < 128; i++)             \
    flags[i] = 0;

#define ARG_FLAGS(set)                                                                  \
  for (k = 1; argv[i][k] != '\0'; k++)                                                  \
    { if (index(set,argv[i][k]) == NULL)                                                \
        { fprintf(stderr,"%s: -%c is an illegal option\n",Prog_Name,argv[i][k]);        \
          exit (1);                                                                     \
        }                                                                               \
      flags[(int) argv[i][k]] = 1;                                                      \
    }

#define ARG_POSITIVE(var,name)                                                          \
  var = strtol(argv[i]+2,&eptr,10);                                                     \
  if (*eptr != '\0' || argv[i][2] == '\0')                                              \
    { fprintf(stderr,"%s: -%c '%s' argument is not an integer\n",			\
                     Prog_Name,argv[i][1],argv[i]+2);      				\
      exit (1);                                                                         \
    }                                                                                   \
  if (var <= 0)                                                                         \
    { fprintf(stderr,"%s: %s must be positive (%d)\n",Prog_Name,name,var);              \
      exit (1);                                                                         \
    }

#define ARG_NON_NEGATIVE(var,name)                                                      \
  var = strtol(argv[i]+2,&eptr,10);                                                     \
  if (*eptr != '\0' || argv[i][2] == '\0')                                              \
    { fprintf(stderr,"%s: -%c '%s' argument is not an integer\n",			\
                     Prog_Name,argv[i][1],argv[i]+2);      				\
      exit (1);                                                                         \
    }                                                                                   \
  if (var < 0)	                                                                        \
    { fprintf(stderr,"%s: %s must be non-negative (%d)\n",Prog_Name,name,var);          \
      exit (1);                                                                         \
    }

#define ARG_REAL(var)                                                                   \
  var = strtod(argv[i]+2,&eptr);                                                        \
  if (*eptr != '\0' || argv[i][2] == '\0')                                              \
    { fprintf(stderr,"%s: -%c '%s' argument is not a real number\n",			\
                     Prog_Name,argv[i][1],argv[i]+2);      				\
      exit (1);                                                                         \
    }


/*******************************************************************************************
 *
 *  GUARDED BATCH IO MACROS
 *
 ********************************************************************************************/

    //  Utilitieis

int Count_Args(char *arg);

#define SYSTEM_READ_ERROR						\
  { fprintf(stderr,"%s: System error, read failed!\n",Prog_Name);	\
    exit (2);								\
  }

#define SYSTEM_WRITE_ERROR						\
  { fprintf(stderr,"%s: System error, write failed!\n",Prog_Name);	\
    exit (2);								\
  }

#define SYSTEM_CLOSE_ERROR						\
  { fprintf(stderr,"%s: System error, file close failed!\n",Prog_Name);	\
    exit (2);								\
  }

    //  Output

#define FWRITE(v,s,n,file)			\
  { if (fwrite(v,s,n,file) != (size_t) n)	\
      SYSTEM_WRITE_ERROR			\
  }

#define FPRINTF(file,...)		\
  { if (fprintf(file,__VA_ARGS__) < 0)	\
      SYSTEM_WRITE_ERROR		\
  }

#define PRINTF(...)			\
  { if (printf(__VA_ARGS__) < 0)	\
      SYSTEM_WRITE_ERROR		\
  }

#define FPUTS(x,file)			\
  { if (fputs(x,file) == EOF)		\
      SYSTEM_WRITE_ERROR		\
  }

    //  Close

#define FCLOSE(file)		\
  { if (fclose(file) != 0)	\
      SYSTEM_CLOSE_ERROR	\
  }

    //  Input

#define FREAD(v,s,n,file)								\
  { if (fread(v,s,n,file) != (size_t) n)						\
      { if (ferror(file))								\
          SYSTEM_READ_ERROR								\
        else										\
          { fprintf(stderr,"%s: The file %s is corrupted\n",Prog_Name,file ## _name);	\
            exit (1);									\
          }										\
      }											\
  }

#define FSCANF(file,...)								\
  { if (fscanf(file,__VA_ARGS__) != Count_Args(#__VA_ARGS__)-1)				\
      { if (ferror(file))								\
          SYSTEM_READ_ERROR								\
        else										\
          { fprintf(stderr,"%s: The file %s is corrupted\n",Prog_Name,file ## _name);	\
            exit (1);									\
          }										\
      }											\
  }

#define FGETS(v,n,file)									\
  { if (fgets(v,n,file) == NULL)							\
      { if (ferror(file))								\
          SYSTEM_READ_ERROR								\
        else										\
          { fprintf(stderr,"%s: The file %s is corrupted\n",Prog_Name,file ## _name);	\
            exit (1);									\
          }										\
      }											\
  }

#define FSEEKO(file,p,d)	\
  { if (fseeko(file,p,d) < 0)	\
      SYSTEM_READ_ERROR		\
  }

#define FTELLO(file)		\
  ( { int x = ftello(file);	\
      if (x < 0)		\
        SYSTEM_READ_ERROR	\
      ; x; 			\
  } )

/*******************************************************************************************
 *
 *  UTILITIES
 *
 ********************************************************************************************/

//  The following general utilities return NULL if any of their input pointers are NULL, or if they
//    could not perform their function (in which case they also print an error to stderr).

void *Malloc(int64 size, char *mesg);                    //  Guarded versions of malloc, realloc
void *Realloc(void *object, int64 size, char *mesg);     //  and strdup, that output "mesg" to
char *Strdup(char *string, char *mesg);                  //  stderr if out of memory

FILE *Fopen(char *path, char *mode);     // Open file path for "mode"
char *PathTo(char *path);                // Return path portion of file name "path"
char *Root(char *path, char *suffix);    // Return the root name, excluding suffix, of "path"

// Catenate returns concatenation of path.sep.root.suffix in a *temporary* buffer
// Numbered_Suffix returns concatenation of left.<num>.right in a *temporary* buffer

char *Catenate(char *path, char *sep, char *root, char *suffix);
char *Numbered_Suffix(char *left, int num, char *right);


// DB-related utilities

void Print_Number(int64 num, int width, FILE *out);   //  Print readable big integer
int  Number_Digits(int64 num);                        //  Return # of digits in printed number

#define COMPRESSED_LEN(len)  (((len)+3) >> 2)

void   Compress_Read(int len, char *s);   //  Compress read in-place into 2-bit form
void Uncompress_Read(int len, char *s);   //  Uncompress read in-place into numeric form
void      Print_Read(char *s, int width);

void Lower_Read(char *s);     //  Convert read from numbers to lowercase letters (0-3 to acgt)
void Upper_Read(char *s);     //  Convert read from numbers to uppercase letters (0-3 to ACGT)
void Number_Read(char *s);    //  Convert read from letters to numbers

void Letter_Arrow(char *s);   //  Convert arrow pw's from numbers to uppercase letters (0-3 to 1234)
void Number_Arrow(char *s);   //  Convert arrow pw string from letters to numbers


/*******************************************************************************************
 *
 *  DB IN-CORE DATA STRUCTURES
 *
 ********************************************************************************************/

#define DB_QV   0x03ff   //  Mask for 3-digit quality value
#define DB_CSS  0x0400   //  This is the second or later of a group of reads from a given insert
#define DB_BEST 0x0800   //  This is the longest read of a given insert (may be the only 1)

#define DB_ARROW 0x2     //  DB is an arrow DB
#define DB_ALL   0x1     //  all wells are in the trimmed DB

//  Fields have different interpretations if a .db versus a .dam

typedef struct
  { int     origin; //  Well # (DB), Contig # (DAM)
    int     rlen;   //  Length of the sequence (Last pulse = fpulse + rlen)
    int     fpulse; //  First pulse (DB), left index of contig in scaffold (DAM)
    int64   boff;   //  Offset (in bytes) of compressed read in 'bases' file, or offset of
                    //    uncompressed bases in memory block
    int64   coff;   //  Offset (in bytes) of compressed quiva streams in '.qvs' file (DB),
                    //  Offset (in bytes) of scaffold header string in '.hdr' file (DAM)
                    //  4 compressed shorts containing snr info if an arrow DB.
    int     flags;  //  QV of read + flags above (DB only)
  } DAZZ_READ;

//  A track can be of 3 types:
//    data == NULL: there are nreads 'anno' records of size 'size'.
//    data != NULL && size == 4: anno is an array of nreads+1 int's and data[anno[i]..anno[i+1])
//                                    contains the variable length data
//    data != NULL && size == 8: anno is an array of nreads+1 int64's and data[anno[i]..anno[i+1])
//                                    contains the variable length data

typedef struct _track
  { struct _track *next;  //  Link to next track
    char          *name;  //  Symbolic name of track
    int            size;  //  Size in bytes of anno records
    void          *anno;  //  over [0,nreads]: read i annotation: int, int64, or 'size' records 
    void          *data;  //     data[anno[i] .. anno[i+1]-1] is data if data != NULL
  } DAZZ_TRACK;

//  The tailing part of a .anno track file can contain meta-information produced by the
//    command that produced the track.  For example, the coverage, or good/bad parameters
//    for trimming, or even say a histogram of QV values.  Each item is an array of 'nelem'
//    64-bit ints or floats ('vtype' = DB_INT or DB_REAL), has a 'name' string that
//    describes it, and an indicator as to whether the values should be equal accross all
//    block tracks, or summed accross all block tracks (by Catrack).  'value' points at the
//    array of values

#define DB_INT  0
#define DB_REAL 1

#define DB_EXACT 0
#define DB_SUM   1

typedef struct
  { int   vtype;  //  INT64 or FLOAST64
    int   nelem;  //  >= 1
    int   accum;  //  EXACT, SUM
    char *name;
    void *value;
  } DAZZ_EXTRA;

//  The information for accessing QV streams is in a DAZZ_QV record that is a "pseudo-track"
//    named ".@qvs" and is always the first track record in the list (if present).  Since normal
//    track names cannot begin with a . (this is enforced), this pseudo-track is never confused
//    with a normal track.

typedef struct
  { struct _track *next;
    char          *name;
    int            ncodes;  //  # of coding tables
    QVcoding      *coding;  //  array [0..ncodes-1] of coding schemes (see QV.h)
    uint16        *table;   //  for i in [0,db->nreads-1]: read i should be decompressed with
                            //    scheme coding[table[i]]
    FILE          *quiva;   //  the open file pointer to the .qvs file
  } DAZZ_QV;

//  The DB record holds all information about the current state of an active DB including an
//    array of DAZZ_READS, one per read, and a linked list of DAZZ_TRACKs the first of which
//    is always a DAZZ_QV pseudo-track (if the QVs have been loaded).

typedef struct
  { int         ureads;     //  Total number of reads in untrimmed DB
    int         treads;     //  Total number of reads in trimmed DB
    int         cutoff;     //  Minimum read length in block (-1 if not yet set)
    int         allarr;     //  DB_ALL | DB_ARROW
    float       freq[4];    //  frequency of A, C, G, T, respectively

    //  Set with respect to "active" part of DB (all vs block, untrimmed vs trimmed)

    int         maxlen;     //  length of maximum read (initially over all DB)
    int64       totlen;     //  total # of bases (initially over all DB)

    int         nreads;     //  # of reads in actively loaded portion of DB
    int         trimmed;    //  DB has been trimmed by cutoff/all
    int         part;       //  DB block (if > 0), total DB (if == 0)
    int         ufirst;     //  Index of first read in block (without trimming)
    int         tfirst;     //  Index of first read in block (with trimming)

       //  In order to avoid forcing users to have to rebuild all thier DBs to accommodate
       //    the addition of fields for the size of the actively loaded trimmed and untrimmed
       //    blocks, an additional read record is allocated in "reads" when a DB is loaded into
       //    memory (reads[-1]) and the two desired fields are crammed into the first two
       //    integer spaces of the record.

    char       *path;       //  Root name of DB for .bps, .qvs, and tracks
    int         loaded;     //  Are reads loaded in memory?
    void       *bases;      //  file pointer for bases file (to fetch reads from),
                            //    or memory pointer to uncompressed block of all sequences.
    DAZZ_READ  *reads;      //  Array [-1..nreads] of DAZZ_READ
    DAZZ_TRACK *tracks;     //  Linked list of loaded tracks
  } DAZZ_DB; 


/*******************************************************************************************
 *
 *  DB STUB FILE FORMAT = NFILE FDATA^nfile NBLOCK PARAMS BDATA^nblock
 *
 ********************************************************************************************/

#define MAX_NAME 10000      //  Longest file name or fasta header line

#define DB_NFILE  "files = %9d\n"   //  number of files
#define DB_FDATA  "  %9d %s %s\n"   //  last read index + 1, fasta prolog, file name
#define DB_NBLOCK "blocks = %9d\n"  //  number of blocks
#define DB_PARAMS "size = %10lld cutoff = %9d all = %1d\n"  //  block size, len cutoff, all in well
#define DB_BDATA  " %9d %9d\n"      //  First read index (untrimmed), first read index (trimmed)


/*******************************************************************************************
 *
 *  DB ROUTINES
 *
 ********************************************************************************************/

  // Suppose DB is the name of an original database.  Then there will be files .DB.idx, .DB.bps,
  //    .DB.qvs, and files .DB.<track>.anno and DB.<track>.data where <track> is a track name
  //    (not containing a . !).

  // A DAM is basically a DB except that:
  //    1. there are no QV's, instead .coff points the '\0' terminated fasta header of the read
  //          in the file .<dam>.hdr file
  //    2. .origin contains the contig # of the read within a fasta entry (assembly sequences
  //          contain N-separated contigs), and .fpulse the first base of the contig in the
  //          fasta entry

  // Open the given database or dam, "path" into the supplied DAZZ_DB record "db". If the name has
  //   a part # in it then just the part is opened.  The index array is allocated (for all or
  //   just the part) and read in.
  // Return status of routine:
  //    -1: The DB could not be opened for a reason reported by the routine to EPLACE
  //     0: Open of DB proceeded without mishap
  //     1: Open of DAM proceeded without mishap

int Open_DB(char *path, DAZZ_DB *db);

  // Trim the DB or part thereof and all loaded tracks according to the cutoff and all settings
  //   of the current DB partition.  Reallocate smaller memory blocks for the information kept
  //   for the retained reads.

void Trim_DB(DAZZ_DB *db);

  // Shut down an open 'db' by freeing all associated space, including tracks and QV structures,
  //   and any open file pointers.  The record pointed at by db however remains (the user
  //   supplied it and so should free it).

void Close_DB(DAZZ_DB *db);

  // Return the size in bytes of the given DB

int64 sizeof_DB(DAZZ_DB *db);

  // If QV pseudo track is not already in db's track list, then load it and set it up.
  //   The database must not have been trimmed yet.  -1 is returned if a .qvs file is not
  //   present, and 1 is returned if an error (reported to EPLACE) occured and INTERACTIVE
  //   is defined.  Otherwise a 0 is returned.

int Load_QVs(DAZZ_DB *db);

  // Remove the QV pseudo track, all space associated with it, and close the .qvs file.

void Close_QVs(DAZZ_DB *db);

  // Look up the file and header in the file of the indicated track.  Return:
  //     1: Track is for trimmed DB
  //     0: Track is for untrimmed DB
  //    -1: Track is not the right size of DB either trimmed or untrimmed
  //    -2: Could not find the track
  // In addition, if opened (0 or 1 returned), then kind points at an integer indicating
  //   the type of track as follows:
  //      CUSTOM  0 => a custom track
  //      MASK    1 => a mask track

#define CUSTOM_TRACK 0
#define   MASK_TRACK 1

int Check_Track(DAZZ_DB *db, char *track, int *kind);

  // If track is not already in the db's track list, then allocate all the storage for it,
  //   read it in from the appropriate file, add it to the track list, and return a pointer
  //   to the newly created DAZZ_TRACK record.  If the track does not exist or cannot be
  //   opened for some reason, then NULL is returned if INTERACTIVE is defined.  Otherwise
  //   the routine prints an error message to stderr and exits if an error occurs, and returns
  //   with NULL only if the track does not exist.

DAZZ_TRACK *Load_Track(DAZZ_DB *db, char *track);

  // Assumming file pointer for afile is correctly positioned at the start of a extra item,
  //   and aname is the name of the .anno file, decode the value present and places it in
  //   extra if extra->nelem == 0, otherwise reduce the value just read into extra according
  //   according the to the directive given by 'accum'.  Leave the read poinrt at the next
  //   extra or end-of-file.
  //   Returns:
  //      1 if at the end of file,
  //      0 if item was read and folded correctly,
  //     -1 if there was a system IO or allocation error (if interactive), and
  //     -2 if the new value could not be reduced into the currenct value of extra (interactive)

int Read_Extra(FILE *afile, char *aname, DAZZ_EXTRA *extra);

//  Write extra record to end of file afile and advance write pointer
//  If interactive, then return non-zero on error, if bash, then print
//  and halt if an error

int Write_Extra(FILE *afile, DAZZ_EXTRA *extra);

  // If track is on the db's track list, then it is removed and all storage associated with it
  //   is freed.

void Close_Track(DAZZ_DB *db, char *track);

  // Allocate and return a buffer big enough for the largest read in 'db'.
  // **NB** free(x-1) if x is the value returned as *prefix* and suffix '\0'(4)-byte
  // are needed by the alignment algorithms.  If cannot allocate memory then return NULL
  // if INTERACTIVE is defined, or print error to stderr and exit otherwise.

char *New_Read_Buffer(DAZZ_DB *db);

  // Load into 'read' the i'th read in 'db'.  As a lower case ascii string if ascii is 1, an
  //   upper case ascii string if ascii is 2, and a numeric string over 0(A), 1(C), 2(G), and 3(T)
  //   otherwise.  A '\0' (or 4) is prepended and appended to the string so it has a delimeter
  //   for traversals in either direction.  A non-zero value is returned if an error occured
  //   and INTERACTIVE is defined.

int  Load_Read(DAZZ_DB *db, int i, char *read, int ascii);

  // Exactly the same as Load_Read, save the arrow information is loaded, not the DNA sequence,
  //   and there is only a choice between numeric (0) or ascii (1);

int  Load_Arrow(DAZZ_DB *db, int i, char *read, int ascii);

  // Load into 'read' the subread [beg,end] of the i'th read in 'db' and return a pointer to the
  //   the start of the subinterval (not necessarily = to read !!! ).  As a lower case ascii
  //   string if ascii is 1, an upper case ascii string if ascii is 2, and a numeric string
  //   over 0(A), 1(C), 2(G), and 3(T) otherwise.  A '\0' (or 4) is prepended and appended to
  //   the string holding the substring so it has a delimeter for traversals in either direction.
  //   A NULL pointer is returned if an error occured and INTERACTIVE is defined.

char *Load_Subread(DAZZ_DB *db, int i, int beg, int end, char *read, int ascii);

  // Allocate a set of 5 vectors large enough to hold the longest QV stream that will occur
  //   in the database.  If cannot allocate memory then return NULL if INTERACTIVE is defined,
  //   or print error to stderr and exit otherwise.

#define DEL_QV  0   //  The deletion QVs are x[DEL_QV] if x is the buffer returned by New_QV_Buffer
#define DEL_TAG 1   //  The deleted characters
#define INS_QV  2   //  The insertion QVs
#define SUB_QV  3   //  The substitution QVs
#define MRG_QV  4   //  The merge QVs

char **New_QV_Buffer(DAZZ_DB *db);

  // Load into 'entry' the 5 QV vectors for i'th read in 'db'.  The deletion tag or characters
  //   are converted to a numeric or upper/lower case ascii string as per ascii.  Return with
  //   a zero, except when an error occurs and INTERACTIVE is defined in which case return wtih 1.

int   Load_QVentry(DAZZ_DB *db, int i, char **entry, int ascii);

  // Allocate a block big enough for all the uncompressed sequences, read them into it,
  //   reset the 'off' in each read record to be its in-memory offset, and set the
  //   bases pointer to point at the block after closing the bases file.  If ascii is
  //   1 then the reads are converted to lowercase ascii, if 2 then uppercase ascii, and
  //   otherwise the reads are left as numeric strings over 0(A), 1(C), 2(G), and 3(T).
  //   Return with a zero, except when an error occurs and INTERACTIVE is defined in which
  //   case return wtih 1.

int Read_All_Sequences(DAZZ_DB *db, int ascii);

  // For the DB or DAM "path" = "prefix/root.[db|dam]", find all the files for that DB, i.e. all
  //   those of the form "prefix/[.]root.part" and call actor with the complete path to each file
  //   pointed at by path, and the suffix of the path by extension.  The . proceeds the root
  //   name if the defined constant HIDE_FILES is set.  Always the first call is with the
  //   path "prefix/root.[db|dam]" and extension "db" or "dam".  There will always be calls for
  //   "prefix/[.]root.idx" and "prefix/[.]root.bps".  All other calls are for *tracks* and
  //   so this routine gives one a way to know all the tracks associated with a given DB.
  //   -1 is returned if the path could not be found, and 1 is returned if an error (reported
  //   to EPLACE) occured and INTERACTIVE is defined.  Otherwise a 0 is returned.

int List_DB_Files(char *path, void actor(char *path, char *extension));

#endif // _DAZZ_DB
