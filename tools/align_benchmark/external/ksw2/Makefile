CC=			gcc
CFLAGS=		-g -Wall -Wextra -Wc++-compat -O2
CPPFLAGS=	-DHAVE_KALLOC
INCLUDES=	-I.
OBJS=		ksw2_gg.o ksw2_gg2.o ksw2_gg2_sse.o ksw2_extz.o ksw2_extz2_sse.o \
			ksw2_extd.o ksw2_extd2_sse.o ksw2_extf2_sse.o ksw2_exts2_sse.o
PROG=		ksw2-test
LIBS=		-lz

ifneq ($(gaba),) # gaba source code directory
	CPPFLAGS += -DHAVE_GABA
	INCLUDES += -I$(gaba)
	LIBS_MORE += -L$(gaba)/build -lgaba
	CFLAGS += -msse4
endif

ifneq ($(parasail),) # parasail install prefix
	CPPFLAGS += -DHAVE_PARASAIL
	INCLUDES += -I$(parasail)/include
	LIBS_MORE += $(parasail)/lib/libparasail.a # don't link against the dynamic library
endif

ifeq ($(sse2),)
	CFLAGS += -march=native
endif

ifneq ($(avx2),)
	CFLAGS += -mavx2
endif

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

ksw2-test:cli.o kalloc.o $(OBJS)
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS_MORE) $(LIBS)
		
clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

cli.o: ksw2.h kseq.h
kalloc.o: kalloc.h
ksw2_extd.o: ksw2.h
ksw2_extd2_sse.o: ksw2.h
ksw2_extf2_sse.o: ksw2.h
ksw2_extz.o: ksw2.h
ksw2_extz2_sse.o: ksw2.h
ksw2_gg.o: ksw2.h
ksw2_gg2.o: ksw2.h
ksw2_gg2_sse.o: ksw2.h
