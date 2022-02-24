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

