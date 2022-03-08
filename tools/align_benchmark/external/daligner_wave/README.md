# Daligner: The Dazzler "Overlap" Module

## _Author:  Gene Myers_
## _First:   April 10, 2016_
## _Current: April 19, 2019_

For typeset documentation, examples of use, and design philosophy please go to
my [blog](https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide).
 
### Version Numbers

v1.0 has been released, but if you need to refer to a later revision
from the stable master branch, please use ``v1.0.yyyymmdd`` where
``yyyy-mm-dd`` is the date of the commit used. This is important for
method details in scientific papers, and for software packaging
(e.g. Conda, HomeBrew, or Linux distribution packages).

The commands below permit one to find all significant local alignments between reads
encoded in Dazzler database.  The assumption is that the reads are from a PACBIO RS II
long read sequencer.  That is the reads are long and noisy, up to 15% on average.

Recall that a database has a current partition that divides it into blocks of a size
that can conveniently be handled by calling the **daligner** overlapper on all the pairs of
blocks producing a collection of .las local alignment files that can then be sorted and
merged into an ordered sequence of sorted files containing all alignments between reads
in the data set.  The alignment records are parsimonious in that they do not record an
alignment but simply a set of trace points, typically every 100bp or so, that allow the
efficient reconstruction of alignments on demand.

All programs add suffixes (e.g. .db, .las) as needed.
For the commands that take multiple .db or .las file blocks as arguments, i.e. **daligner**, **LAsort**, **LAmerge**, **LAcat**,
and **LAcheck**, one can place a @-sign in the name, which is then interpreted as the sequence of files
obtained by replacing the @-sign by 1, 2, 3, ... in sequence until a number is reached for
which no file matches.  One can also place a @-sign followed by an integer, say, i, in which
case the sequence starts at i.  Lastly, one can also use @i-j where i and j are integers, in
which case the sequence is from i to j, inclusive.

The formal UNIX command line
descriptions and options for the DALIGNER module commands are as follows:

```
1. daligner [-vaAI]
       [-k<int(16)>] [-%<int(28)>] [-h<int(50)>] [-w<int(6)>] [-t<int>] [-M<int>]
       [-e<double(.75)] [-l<int(1500)] [-s<int(100)>] [-H<int>]
       [-T<int(4)>] [-P<dir(/tmp)>] [-m<track>]+
       <subject:db|dam> <target:db|dam> ...
```

Compare sequences in the trimmed *\<subject\>* block against those in the list of *\<target\>*
blocks searching for local alignments involving at least -l base pairs (default 1000)
or more, that have an average correlation rate of -e (default 70%).  The local
alignments found will be output in a sparse encoding where a trace point on the
alignment is recorded every -s base pairs of the a-read (default 100bp). Reads are
compared in both orientations and local alignments meeting the criteria are output to
one of several created files described below.  The -v option turns on a verbose
reporting mode that gives statistics on each major step of the computation.  The
program runs with 4 threads by default, but this may be set to any positive value with
the -T option.

The options -k, -%, -h, and -w control the initial filtration search for possible matches
between reads.  Specifically, our search code looks for a pair of diagonal bands of
width 2<sup>w</sup> (default 2<sup>6</sup> = 64) that contain a collection of matching k-mers
(default 16) in the lowest %-percentifle between the two reads, such that the total number of bases covered by the k-mer hits is h (default 50). k cannot be larger than 32 in the current implementation.  *These parameters will shortly be superceded with a more intuitive interface.*

If there are one or more interval tracks specified with the -m option, then the reads
of the DB or DB's to which the mask applies are soft masked with the union of the
intervals of all the interval tracks that apply, that is any k-mers that contain any
bases in any of the masked intervals are ignored for the purposes of seeding a match.
An interval track is a track, such as the "dust" track created by DBdust, that encodes
a set of intervals over either the untrimmed or trimmed DB.

Invariably, some k-mers are significantly over-represented (e.g. homopolymer runs).
These k-mers create an excessive number of matching k-mer pairs and left unaddressed
would cause daligner to overflow the available physical memory.  One way to deal with
this is to explicitly set the -t parameter which suppresses the use of any k-mer that
occurs more than t times in either the subject or target block.  However, a better way
to handle the situation is to let the program automatically select a value of t that
meets a given memory usage limit specified (in Gb) by the -M parameter.  By default
daligner will use the amount of physical memory as the choice for -M.  If you want to
use less, say only 8Gb on a 24Gb HPC cluster node because you want to run 3 daligner
jobs on the node, then specify -M8.  Specifying -M0 basically indicates that you do not
want daligner to self adjust k-mer suppression to fit within a given amount of memory.

Each found alignment is recorded as -- a[ab,ae] x b<sup>o</sup>[bb,be] -- where a and b are the
indices (in the trimmed DB) of the reads that overlap, o indicates whether the b-read
is from the same or opposite strand, and [ab,ae] and [bb,be] are the intervals of a
and bo, respectively, that align.  For each subject, target pair of blocks, say X and Y,
the program reports alignments where the a-read is in X and the b-read is in Y, or
vice versa.  However, if the -A option is set ("A" for "asymmetric") then just overlaps
where the a-read is in X and the b-read is in Y are reported, and if X = Y then it
further reports only those overlaps where the a-read index is less than the b-read index.
In either case, if the -I option is set ("I" for "identity") then when X = Y, overlaps
between different portions of the same read will also be found and reported.  In summary,
the command `daligner -A X Y` produces a single file `X.Y.las` and `daligner X Y` produces
2 files `X.Y.las` and `Y.X.las` (unless X=Y in which case only a single file, `X.X.las`, is
produced).  The overlap records in one of these files are sorted as described for LAsort.
The -a option to daligner is passed directly through to LAsort which is actually called
as a sub-process to produce the sorted file.
In order to produce the aforementioned .las file, several temporary .las files, two for
each thread, are produce in the sub-directory /tmp by default.  You can overide this
location by specifying the directory you would like this activity to take place in with
the -P option.

By default daligner compares all overlaps between reads in the database that are
greater than the minimum cutoff set when the DB or DBs were split, typically 1 or
2 Kbp.  However, the HGAP assembly pipeline only wants to correct large reads, say
8Kbp or over, and so needs only the overlaps where the a-read is one of the large
reads.  By setting the -H parameter to say N, one alters daligner so that it only
reports overlaps where the a-read is over N base-pairs long.

While the default parameter settings are good for raw Pacbio data, daligner can be used
for efficiently finding alignments in corrected reads or other less noisy reads. For
example, for mapping applications against .dams we run `daligner -k20 -h60 -e.85` and
on corrected reads, we typically run `daligner -k25 -w5 -h60 -e.95 -s500` and at
these settings it is very fast.

```
2. LAsort [-va] <align:las> ...
```

Sort each .las alignment file specified on the command line. For each file it reads in
all the overlaps in the file and sorts them in lexicographical order of (a,b,o,ab)
assuming each alignment is recorded as a[ab,ae] x b<sup>o</sup>[bb,be]. It then writes them all
to a file named \<align\>.S.las (assuming that the input file was \<align\>.las). With the
-v option set then the program reports the number of records read and written. If the
-a option is set then it sorts LAs in lexicographical order of (a,ab) alone, which is
desired when sorting a mapping of reads to a reference.

If the .las file was produced by damapper the local alignments are organized into
chains where the LA segments of a chain are consecutive and ordered in the file.
LAsort can detects that it has been passed such a file and if so treats the chains as
a unit and sorts them on the basis of the first LA in the chain.

```
3. LAmerge [-va] [-P<dir(/tmp)>] <merge:las> <parts:las> ...
```

Merge the .las files \<parts\> into a singled sorted file \<merge\>, where it is assumed
that  the input \<parts\> files are sorted.  There are no limits to how many files can be
merged, but if there are more than 252, a typical UNIX OS limit on the number of simultaneously
open files, then the program recursively spawns sub-processes and creates temporary files
in the directory specified by the -P option, /tmp by default.
With the -v option set the program reports the number of
records read and written.  The -a option indicates the sort is as describe for LAsort
above.

If the .las file was produced by damapper the local alignments are organized into
chains where the LA segments of a chain are consecutive and ordered in the file.  When
merging such files, LAmerge treats the chains as a unit and orders them on the basis
of the first LA in the chain.

Used correctly, LAmerge and LAsort together allow one to perform an "external" sort
that produces a collection of sorted files containing in aggregate all the local
alignments found by the daligner, such that their concatenation is sorted in order of
(a,b,o,ab) (or (a,ab) if the -a option is set). In particular, this means that all the
alignments for a given a-read will be found consecutively in one of the files.  So
computations that need to look at all the alignments for a given read can operate in
simple sequential scans of these sorted files.

```
4. LAshow [-caroUF] [-i<int(4)>] [-w<int(100)>] [-b<int(10)>]
                    <src1:db|dam> [ <src2:db|dam> ]
                    <align:las> [ <reads:FILE> | <reads:range> ... ]
```

LAshow produces a printed listing of the local alignments contained in the specified
.las file, where the a- and b-reads come from src1 or from src1 and scr2, respectively.
If a file or list of read ranges is given then only the overlaps for which the a-read
is in the set specified by the file or list are displayed. See DBshow for an explanation
of how the file and list of read ranges are interpreted.  If the -F option is set then
the roles of the a- and b- reads are reversed in the display.

If the -c option is given then a cartoon rendering is displayed, and if -a or -r option
is set then an alignment of the local alignment is displayed.  The -a option puts
exactly -w columns per segment of the display, whereas the -r option puts exactly -w
a-read symbols in each segment of the display.  The -r display mode is useful when one
wants to visually compare two alignments involving the same a-read.  If a combination of
the -c, -a, and -r flags is set, then the cartoon comes first, then the -a alignment,
and lastly the -r alignment.  The -i option sets the indent for the cartoon and/or
alignment displays, if they are requested.  The -b option sets the number of symbols on
either side of the aligned segments in an alignment display, and -U specifies that
uppercase should be used for DNA sequence instead of the default lowercase.  If the
-o option is set then only alignments that are proper overlaps (a sequence end occurs
at the each end of the alignment) are displayed.  If the -F option is given then the
roles of the A- and B-reads are flipped.

When examining LAshow output it is important to keep in mind that the coordinates
describing an interval of a read are referring conceptually to positions between bases
starting at 0 for the position to the left of the first base.  That is, a coordinate c
refers to the position between the c-1'st and c'th base, and the interval [b,e] captures
the e-b bases from the b'th to the e-1'st, inclusive.  We give an example with a cartoon
and (part of an) alignment for which we will explain several additional
important points:

```
 1  1,865 c   [18,479..20,216] x [ 1,707..0>  (24,451 x 7,283 bps, 19 trace pts)

      18479              4235
  A ========+----------+======>  dif/(len1+len2) = 478/(1737+1707) = 27.76%
  B  <======+-----------
       5576

   18469 agccgcctag[tgcctcgcaaacgc-t-cggggcggcgt-gaaagcgg--
         ::::::::::[||||||||||||||*|*|||*|||*|||*||||||||**
    1717 ctcttcttta[tgcctcgcaaacgccttcggcgcg-cgttgaaagcggtt  17.9%

   18513 -ccggtgggtc--agtggcgagttctggcagtgcgctggg-ctgcgaaat
         *||||||*|||**|||||*||||*|*|*|||**|||||||*||*||||||
   1669 gccggtgcgtcgcagtgg-gagt-c-gtcag--cgctggggcttcgaaat  24.0%

        . . .
```

The display of an LA always begins with a line giving the A-read, then the B-read, then
an indication of orientation (i.e. 'n' for same strand, and 'c' for the opposite strand)
followed by the A-interval and B-interval that are aligned and in parentheses
the lengths of the two reads and the number of tracepoints in the alignment between them.
In particular,
note carefully that when the B-read is in the complement orientation (c), then the
B-interval gives the higher coordinate first, the idea being that one will align from
the highest base down to the lowest base in the descending direction on B, complement
the characters as you go.  Further note that in the alignment display the coordinates at
the start of each line follow this orientation convention and give the coordinate of the
"tick mark" just left of the first character in each line.  It is useful to know if an
interval reaches the end of read, and to signal this we use an angle-bracket \<\> instead
of a square bracket [], e.g. in the example the B-segment starts at the beginning of the
read.  Finally, observe that in the cartoon the numbers are not coordinates but rather
indicate the lengths of the unaligned bits left and right of the two aligned intervals.
Finally, observe that in the cartoon the numbers are not coordinates but rather indicate
the lengths of the unaligned bits left and right of the two aligned intervals.

With the introduction of damapper, .las files can now contain chains.  If LAshow detects
that it has been passed a file with chain information then it displays marks at the left
that reveal the chain structure, e.g.:

```
   >     117   37,630 c   [   253.. 7,980] x [   331,430..   324,027]  ~  10.5%
   +     117   37,628 n   [   253.. 7,983] x [21,493,673..21,501,079]  ~  10.6%
   +     117       57 c   [   253.. 1,086] x [ 2,008,164.. 2,007,369]  ~   9.8%
    -    117       57 c   [ 1,300.. 7,982] x [ 2,007,351.. 2,000,945]  ~  10.7%
   >     117       15 c   [ 7,992.. 8,716] x [   242,529..   241,822]  ~   7.8%
    -    117       15 c   [ 8,752..14,299] x [   241,824..   236,425]  ~  10.7%
    -    117       15 c   [14,133..14,832] x [   236,630..   235,953]  ~  12.1%
   +     117   37,628 n   [ 7,992.. 8,716] x [19,202,357..19,203,064]  ~   7.7%
    -    117   37,628 n   [ 8,752..14,832] x [19,203,062..19,208,974]  ~  10.9%
```

A chain begins with either a > or + character, where > indicates this is the highest
scoring chain and + indicates an alternate near optimal chain (controlled by the
-n parameter to damapper).  Each additional LA of a chain is marked with a - character.

```
5a. LAdump [-cdtlo] <src1:db|dam> [ <src2:db|dam> ]
                   <align:las> [ <reads:FILE> | <reads:range> ... ]

5b. dumpLA <align.las>
```

Like LAshow, LAdump allows one to display the local alignments (LAs) of a subset of the
piles in an .las file and select which information to show about them.  The difference
is that the information is written in a very simple "1-code" ASCII format that makes it
easy for one to read and parse the information for further use.  For each LA the pair of
reads is output on a line.  -c requests that one further output the coordinates of the
LA segments be output.  The -d option requests that the number of difference in the LA
be output, -t requests that the tracepoint information be output, and -l requests the
length of the two reads be output.  Finally, -o requests that only LAs that are proper
overlaps be output. 

The format is very simple.  Each requested piece of information occurs on a line.  The
first character of every line is a "1-code" character that tells you what information
to expect on the line.  The rest of the line contains information where each item is
separated by a single blank space.  The trace point line gives the number of trace
point intervals in the LA and is immediately followed by that many lines containing
a pair of integers giving the number of differences and b-displacement in each successive
trace point interval.

```
    P #a #b #o #c     - (#a,#b^#o) have an LA between them where #o is 'n' or 'c' and
                           #c is '>' (start of best chain), '+' (start of alternate chain),
                           '-' (continuation of chain), or '.' (no chains in file).
    L #la #lb         - #la is the length of the a-read and #lb that of the b-read
    C #ab #ae #bb #be - #a[#ab,#ae] aligns with #b^#o[#bb,#be]
    D #               - there are # differences in the LA
    T #n              - there are #n trace point intervals for the LA
     (#d #y )^#n      - there are #d difference aligning the #y bp's of B with the
                           next fixed-size interval of A
    + X #             - Total amount of X (X = P or T)
    % X #             - Maximum amount of X in any pile (X = P or T)
    @ T #             - Maximum number of trace points in any trace
```

1-code lines that begin with +, %, or @ are always the first lines in the output.
They give size information about what is contained in the output.  Specifically,
'+ X #' gives the total number of LAs (X=P), or the total number of trace point
intervals (X=T) in the file .  '% X #' gives the maximum number of LAs (X=P) or
the maximum number of trace point intervals (X=T) in a given *pile* (collection of
LAs all with the same a-read (applies only to sorted .las files).  A final line: '@ T #',
gives the maximum # of trace point intervals in any trace within the file.
After these lines and before the start of the lines describing alignment records is a
single line of the form 'X #' where the number is the trace point spacing for all
alignments.

The command dumpLA reads a 1-code file from the standard input and if possible produces a .las
file for it.  The 1-code file is any legitimate coding of alignments as might be produced by LAdump.
The 1-code file must contain the P-, C-, and T-lines as well as the X-line and the header lines
beginning with +, %, or @.  So for example, a 1-code file produced by LAdump with the -c and -t
options is invertible.

```
6a. LAa2b
6b. LAb2a
```

Pipes (stdin to stdout) that convert an ASCII output produced by LAdump into a compressed
binary representation (LAa2b) and vice verse (LAb2a).  The idea is to save disk space by
keeping the dumps in a more compressed format.

```
7. LAcat [-v] <source:las> ... > <target>.las
```

The sequence of \<source\> files (that can contain @-sign block ranges) are
concatenated in order
into a single .las file and pipe the result to the standard output.  The -v
option reports the files concatenated and the number of la's within them to
standard error (as the standard output receives the concatenated file).

```
8. LAsplit [-v] <target:las> (<parts:int> | <path:db|dam>) < <source>.las
```

If the second argument is an integer n, then divide the alignment file \<source\>, piped
in through the standard input, as evenly as possible into n alignment files with the
names specified by template \<target\>, subject to the restriction that all alignment
records for a given a-read are in the same file.  The name of the n files is the
string \<target\> where the single @-sign that occurs somewhere in it is replaced
by i for i in [1,n] and a .las extension is added if necessary.

If the second argument refers to a database \<path\>.db that has been partitioned, then
divide the input alignment file into block .las files where all records whose a-read is
in \<path\>.i.db are in the i'th file generated from the template \<target\>.  The -v
option reports the files produced and the number of la's within them to standard error.

```
9. LAcheck [-vaS] <src1:db|dam> [ <src2:db|dam> ] <align:las> ...
```

LAcheck checks each .las file for structural integrity, where the a- and b-sequences
come from src1 or from src1 and scr2, respectively.  That is, it makes sure each file
makes sense as a plausible .las file, e.g. values are not out of bound, the number of
records is correct, the number of trace points for a record is correct, and so on.  If
the -S option is set then it further checks that the alignments are in sorted order,
by default pile order, but if -a is also set, then map order.
If the -v option is set then a line is output for each .las file saying either the
file is OK or reporting the first error.  If the -v option is not set then the program
runs silently.  The exit status is 0 if every file is deemed good, and 1 if at least
one of the files looks corrupted.

With the introduction of damapper, LAcheck checks to see if a file has chain
information, and if it does, then it checks the validity of chains and checks the
sorting order of chains as a unit according to the -a option.

```
10. HPC.daligner [-vad] [-t<int>] [-w<int(6)>] [-l<int(1500)] [-s<int(100)] [-M<int>]
                    [-P<dir(/tmp)>] [-B<int(4)>] [-T<int(4)>] [-f<name>]
                  ( [-k<int(16)>] [-h<int(50)>] [-e<double(.75)] [-H<int>]
                    [-k<int(20)>] [-h<int(50)>] [-e<double(.85)]  <ref:db|dam>  )
                    [-m<track>]+ <reads:db|dam> [<first:int>[-<last:int>]]
```

HPC.daligner writes a UNIX shell script to the standard output or to a series of files
beginning with the prefix \<name\> if the -f option is set, that either performs an
"overlap" computation on all the blocks in a single database, or a "comparison"
computation on all pairs of blocks between two databases, depending on whether it is
given one or two DB's as arguments (\<ref\> and \<reads\>).  We describe the overlap
script first and its effect first and then later the comparison script.

An Overlap Script: consists of a sequence of commands that effectively run daligner on
all pairs of blocks of a split database and then externally sorts and merges them using
LAsort and LAmerge into a collection of alignment files with names \<path\>.#.las where #
ranges from 1 to the number of blocks the data base is split into. These sorted files
if concatenated by say LAcat would contain all the alignments in sorted order (of
a-read, then b-read, ...).  Moreover, all overlaps for a given a-read are guaranteed
to not be split across files, so one can run artifact analyzers or error correction on
each sorted file in parallel.

The data base must have been previously split by DBsplit and all the parameters, except
-a, -d, -f, -B, and -D, are passed through to the calls to daligner. The defaults for
these parameters are as for daligner. The -v and -a flags are passed to all calls to
LAsort and LAmerge. All other options are described later. For a database divided into
N sub-blocks, the calls to daligner will produce in total N<sup>2</sup> .las files,
on per block pair.
These are then merged so that there is 1 file per row of
the N x N block matrix. So at the end one has N sorted .las files, one per block of
A-reads, that when
concatenated would give a single large sorted overlap file.

The -B option (default 4) gives the desired number of block comparisons per call to
daligner. Some must contain B-1 comparisons, and the first B-2 block comparisons
even less, but the HPCdaligner "planner" does the best it can to give an average load
of -B block comparisons per command.

If the integers \<first\> and \<last\> are missing then the script produced is for every
block in the database.  If \<first\> is present then HPCdaligner produces an incremental
script that compares blocks \<first\> through \<last\> (\<last\> = \<first\> if not present)
against each other and all previous blocks 1 through \<first\>-1, and then incrementally
updates the .las files for blocks 1 through \<first\>-1, and creates the .las files for
blocks \<first\> through \<last\>.

A Comparison Script: consists of a sequence of commands that effectively maps every
read in the DB \<reads\> against a reference set of sequences in the DB \<ref\>, recording
all the found local alignments in the sequence of files \<reads\>.1.\<ref\>.las,
\<reads\>.2.\<ref\>.las, ... where \<reads\>.\<ref\>.k.las contains the alignments between all
of \<ref\> and the k'th block of \<reads\>. The parameters are exactly the same as for the
overlap script save that the -k, -h, and -e defaults are set more stringently for
mapping, and the -A, -I , and -H options make no sense as \<ref\> and \<reads\> are
expected to be distinct data sets.  If the integers \<first\> and \<last\> are missing then
the script produced is for every block in the database \<reads\>. If \<first\> is present
then HPC.daligner produces a script that compares blocks \<first\> through \<last\> (\<last\>
= \<first\> if not present) of \<reads\> against DAM \<ref\>.

The command scripts output by HPC.daligner and other HPC.\<x\> programs consists of
command blocks each of which begins with a comment line (begins with #) followed by a
potentially long list of lines each containing a shell command.  Command blocks whose
comment mentions "jobs" and gives the number of said in parenthesis, we call parallel
blocks because each command line in the block can be sent to a node in a cluster for
independent execution, i.e. none of the commands in a block depend on another in the
block.  The remaining command blocks we call house-keeping blocks because they can be
executed by the shell on the launch/server node and the commands are either checking
the integrity of .las files with LAcheck, or removing intermediate files with rm. Each
block should be performed in the order given and should complete before the next block
is performed.

If the -f option is set, then each command block is written to a file with a name of
the form \<name\>.#.\<description\> where \<name\> is specified by the user in the -f option
argument, # gives the order in which the command block in the given file is to be
performed in relation to other command block files, and \<description\> is a (very)
short symbolic reminder of what the block is doing.  For example, "HPC.daligner -fJOBS
DB" would produce the files:

```
     JOBS.01.OVL
     JOBS.02.CHECK.OPT
     JOBS.03.MERGE
     JOBS.04.RM.OPT
```

There are always 4 command blocks.  The files with the suffix .OPT are
optional and need not be executed albeit we highly recommend that one run the
CHECK block.  One should *not* run the RM block if one wants to later use
DASrealign after scrubbing.

A new -d option requests scripts that organize files into a collection of
sub-directories so as not to overwhelm the underlying OS for large genomes.  Recall
that for a DB divided into N blocks, the daligner will produce N<sup>2</sup> .las-files.
With the -d option set, N sub-directories (with respect to the directory HPC.daligner is
called in) of the form "work\<i\>" for i from 1 to N are created in an initial command
block, and then all work files are placed in those sub-directories, with a maximum
of 2N files appearing in any sub-directory at any given point in the process.

Example:

```
//  Recall G.db from the example in DAZZ_DB/README

> cat G.db
files =         1
       1862 G Sim
blocks =         2
size =        11 cutoff =         0 all = 0
         0         0
      1024      1024
      1862      1862
> HPCdaligner -mdust -t5 G | csh -v   // Run the HPCdaligner script

# Dazzler jobs (2)
dazzler -d -t5 -mdust G.1 G.1
dazzler -d -t5 -mdust G.2 G.1 G.2
# Initial sort jobs (4)
LAsort G.1.G.1.*.las && LAmerge G.L1.1.1 G.1.G.1.*.S.las && rm G.1.G.1.*.S.las
LAsort G.1.G.2.*.las && LAmerge G.L1.1.2 G.1.G.2.*.S.las && rm G.1.G.2.*.S.las
LAsort G.2.G.1.*.las && LAmerge G.L1.2.1 G.2.G.1.*.S.las && rm G.2.G.1.*.S.las
LAsort G.2.G.2.*.las && LAmerge G.L1.2.2 G.2.G.2.*.S.las && rm G.2.G.2.*.S.las
# Level 1 jobs (2)
LAmerge G.1 G.L1.1.1 G.L1.1.2 && rm G.L1.1.1.las G.L1.1.2.las
LAmerge G.2 G.L1.2.1 G.L1.2.2 && rm G.L1.2.1.las G.L1.2.2.las

> LAshow -c -a:G -w50 G.1 | more  // Take a look at the result !

G.1: 34,510 records

         1          9 c   [     0.. 1,876] x [ 9,017..10,825]  ( 18 trace pts)

                      12645
    A      ---------+====>   dif/(len1+len2) = 398/(1876+1808) = 21.61%
    B <====+---------
       9017

         1 ..........gtg-cggt--caggggtgcctgc-t-t-atcgcaatgtta
                     |||*||||**||||||||*||||*|*|*||**|*|*||||
      9008 gagaggccaagtggcggtggcaggggtg-ctgcgtcttatatccaggtta  27.5%

        35 ta-ctgggtggttaaacttagccaggaaacctgttgaaataa-acggtgg
           ||*|||||||||||||*|**|*||*|*||||||*|**|||||*|*|||||
      9057 tagctgggtggttaaa-tctg-ca-g-aacctg-t--aataacatggtgg  24.0%

        83 -ctagtggcttgccgtttacccaacagaagcataatgaaa-tttgaaagt
           *||||||||*||||||||*||**||||*|||**|||||||*||||*||||
      9100 gctagtggc-tgccgttt-ccgcacag-agc--aatgaaaatttg-aagt  20.0%

       131 ggtaggttcctgctgtct-acatacagaacgacggagcgaaaaggtaccg
           ||*|||||||||||||*|*||||*|*|*||||||||||*||||||||||*
      9144 gg-aggttcctgctgt-tcacat-c-ggacgacggagc-aaaaggtacc-  16.0%

...

> LAcat G >G.las       //  Combine G.1.las & G.2.las into a single .las file
> LAshow G G | more    //   Take another look, now at G.las

G: 62,654 records
   1    9 c   [     0.. 1,876] x [ 9,017..10,825] :   <    398 diffs  ( 18 trace pts)
   1   38 c   [     0.. 7,107] x [ 5,381..12,330] :   <  1,614 diffs  ( 71 trace pts)
   1   49 n   [ 5,493..14,521] x [     0.. 9,065] :   <  2,028 diffs  ( 91 trace pts)
   1   68 n   [12,809..14,521] x [     0.. 1,758] :   <    373 diffs  ( 17 trace pts)
   1  147 c   [     0..13,352] x [   854..14,069] :   <  2,993 diffs  (133 trace pts)
   1  231 n   [10,892..14,521] x [     0.. 3,735] :   <    816 diffs  ( 37 trace pts)
   1  292 c   [ 3,835..14,521] x [     0..10,702] :   <  2,353 diffs  (107 trace pts)
   1  335 n   [ 7,569..14,521] x [     0.. 7,033] :   <  1,544 diffs  ( 70 trace pts)
   1  377 c   [ 9,602..14,521] x [     0.. 5,009] :   <  1,104 diffs  ( 49 trace pts)
   1  414 c   [ 6,804..14,521] x [     0.. 7,812] :   <  1,745 diffs  ( 77 trace pts)
   1  415 c   [     0.. 3,613] x [ 7,685..11,224] :   <    840 diffs  ( 36 trace pts)
   1  445 c   [ 9,828..14,521] x [     0.. 4,789] :   <  1,036 diffs  ( 47 trace pts)
   1  464 n   [     0.. 1,942] x [12,416..14,281] :   <    411 diffs  ( 19 trace pts)

...
```
