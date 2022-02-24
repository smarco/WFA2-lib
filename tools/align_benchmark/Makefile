###############################################################################
# Definitions
###############################################################################
FOLDER_WFA=../../..
FOLDER_ROOT=../..
FOLDER_LIB=$(FOLDER_ROOT)/lib
FOLDER_BIN=$(FOLDER_ROOT)/bin
FOLDER_BUILD=./build

ifeq ($(UNAME), Linux)
  LD_FLAGS+=-lrt 
endif

LI_FLAGS=-L$(FOLDER_LIB) -I$(FOLDER_ROOT) -I$(FOLDER_WFA) -I.
LIB_WFA=$(FOLDER_LIB)/libwfa.a
LIB_BLOCKALIGNER=external/block-aligner/target/release/
LIB_EDLIB=external/edlib/build/lib/
LIB_GABA=external/libgaba/libgaba.a
LIB_PARASAIL=external/parasail/build/libparasail.a
###############################################################################
# Rules
###############################################################################
SUBDIRS=benchmark \
        edit \
        gap_affine \
        gap_affine2p \
        gap_linear \
        indel
BIN=$(FOLDER_BIN)/align_benchmark
        
all: setup 
all: $(SUBDIRS)
all: OBJS+=$(FOLDER_BUILD)/*.o
all: OBJS+=external/build/*.o
all: LD_FLAGS+=-L$(LIB_BLOCKALIGNER) -Wl,-rpath=$(PWD)/tools/align_benchmark/$(LIB_BLOCKALIGNER) -lblock_aligner # BlockAligner
all: LD_FLAGS+=-L$(LIB_EDLIB) -Wl,-rpath=$(PWD)/tools/align_benchmark/$(LIB_EDLIB) -ledlib # Edlib
all: LD_FLAGS+=$(LIB_GABA) # Gaba
all: LD_FLAGS+=$(LIB_PARASAIL) # Parasail
all: LD_FLAGS+=-lstdc++ # SeqAn

all: align_benchmark

align_benchmark: */*.c */*.h align_benchmark.c $(LIB_WFA)
	$(CC) $(CC_FLAGS) $(LI_FLAGS) $(OBJS) align_benchmark.c -o $(BIN) $(LD_FLAGS) -lwfa -lm -fopenmp
	
setup:
	@mkdir -p $(FOLDER_BUILD)

clean: 
	rm -rf $(FOLDER_BUILD)
	
###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
	
.PHONY: $(SUBDIRS)

