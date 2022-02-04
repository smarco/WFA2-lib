###############################################################################
# Flags & Folders
###############################################################################
FOLDER_BIN=bin
FOLDER_BUILD=build
FOLDER_BUILD_CPP=build/cpp
FOLDER_LIB=lib

UNAME=$(shell uname)

CC=gcc
CPP=g++

LD_FLAGS=-lm
CC_FLAGS=-Wall -g
ifeq ($(UNAME), Linux)
  LD_FLAGS+=-lrt 
endif

AR=ar
AR_FLAGS=-rsc

###############################################################################
# Compile rules
###############################################################################
LIB_WFA=$(FOLDER_LIB)/libwfa.a
LIB_WFA_CPP=$(FOLDER_LIB)/libwfacpp.a
SUBDIRS=alignment \
        benchmark \
        bindings/cpp \
        edit \
        gap_affine \
        gap_affine2p \
        gap_lineal \
        system \
        utils \
        wavefront

release: CC_FLAGS+=-O3 -march=native -flto
release: MODE=all
release: setup
release: $(SUBDIRS) lib_wfa tools

all: CC_FLAGS+=-O3 -march=native
all: MODE=all
all: setup
all: $(SUBDIRS) lib_wfa tools

debug: setup
debug: MODE=all
debug: $(SUBDIRS) lib_wfa tools

# ASAN: ASAN_OPTIONS=detect_leaks=1:symbolize=1 LSAN_OPTIONS=verbosity=2:log_threads=1
asan: CC_FLAGS+=-fsanitize=address -fno-omit-frame-pointer -fno-common
asan: MODE=all
asan: setup
asan: $(SUBDIRS) lib_wfa tools

setup:
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_BUILD_CPP) $(FOLDER_LIB)
	
lib_wfa: $(SUBDIRS)
	$(AR) $(AR_FLAGS) $(LIB_WFA) $(FOLDER_BUILD)/*.o 2> /dev/null
	$(AR) $(AR_FLAGS) $(LIB_WFA_CPP) $(FOLDER_BUILD)/*.o $(FOLDER_BUILD_CPP)/*.o 2> /dev/null

clean:
	rm -rf $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)
	
###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
	
tools: $(SUBDIRS)
	$(MAKE) --directory=$@ $(MODE)

.PHONY: $(SUBDIRS) tools

