CC = gcc
AR = ar
GIT = git
ARCHFLAGS = -march=native
BIT = 4
CFLAGS = -O3 -Wall -Wno-unused-function -Wno-unused-label -std=c99 -pipe -DBIT=$(BIT)
TARGET = libgaba.a

all: native example unittest bench

native:
	$(MAKE) -f Makefile.core CC=$(CC) CFLAGS='$(CFLAGS) $(ARCHFLAGS)' all
	$(AR) rcs $(TARGET) gaba.*.o gaba_common.o

example: example.c native
	$(CC) -o $@ $(CFLAGS) $(ARCHFLAGS) $< $(TARGET)

unittest: unittest.c native
	$(CC) -o $@ $(CFLAGS) $(ARCHFLAGS) $< gaba.*.o

bench: bench.c gaba.c
	$(CC) -o $@ $(CFLAGS) $(ARCHFLAGS) $^ -DBW=32 -DMODEL=AFFINE -DBENCH

debug: gaba.c
	$(CC) -o $@ $(CFLAGS:-O3=-g) $(ARCHFLAGS) $^ -DBW=32 -DMODEL=AFFINE -DDEBUG -DUNITTEST -DUNITTEST_ALIAS_MAIN

clean:
	rm -rf *.o $(TARGET) example unittest bench debug *~ *.a *.dSYM session*
