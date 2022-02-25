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

CC_FLAGS=-Wall -g

AR=ar
AR_FLAGS=-rsc

ifndef BUILD_EXAMPLES 
BUILD_EXAMPLES=1
endif

ifndef BUILD_TOOLS 
BUILD_TOOLS=1
endif
###############################################################################
# Configuration rules
###############################################################################
LIB_WFA=$(FOLDER_LIB)/libwfa.a
LIB_WFA_CPP=$(FOLDER_LIB)/libwfacpp.a
SUBDIRS=alignment \
        bindings/cpp \
        system \
        utils \
        wavefront
ifeq ($(BUILD_TOOLS),1)        
    APPS+=tools/generate_dataset \
          tools/align_benchmark
endif
ifeq ($(BUILD_EXAMPLES),1)        
    APPS+=examples
endif

release: CC_FLAGS+=-O3 -march=native -flto
release: build

all: CC_FLAGS+=-O3 -march=native
all: build

debug: build

# ASAN: ASAN_OPTIONS=detect_leaks=1:symbolize=1 LSAN_OPTIONS=verbosity=2:log_threads=1
asan: CC_FLAGS+=-fsanitize=address -fno-omit-frame-pointer -fno-common
asan: build

###############################################################################
# Build rules
###############################################################################
build: setup
build: $(SUBDIRS) 
build: lib_wfa 
build: $(APPS)

setup:
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_BUILD_CPP) $(FOLDER_LIB)
	
lib_wfa: $(SUBDIRS)
	$(AR) $(AR_FLAGS) $(LIB_WFA) $(FOLDER_BUILD)/*.o 2> /dev/null
	$(AR) $(AR_FLAGS) $(LIB_WFA_CPP) $(FOLDER_BUILD)/*.o $(FOLDER_BUILD_CPP)/*.o 2> /dev/null

clean:
	rm -rf $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)
	$(MAKE) --directory=tools/align_benchmark clean
	$(MAKE) --directory=examples clean
	
###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
	
$(APPS):
	$(MAKE) --directory=$@ all

.PHONY: $(SUBDIRS) $(APPS)

