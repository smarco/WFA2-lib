###############################################################################
# Definitions
###############################################################################
FOLDER_WFA=../../..
FOLDER_ROOT=..
FOLDER_BUILD=$(FOLDER_ROOT)/build

###############################################################################
# Modules
###############################################################################
MODULES=benchmark_check \
        benchmark_edit \
        benchmark_gap_affine \
        benchmark_gap_affine2p \
        benchmark_gap_linear \
        benchmark_indel \
        benchmark_utils
        
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

###############################################################################
# Rules
###############################################################################
SUBDIRS=external

all: $(SUBDIRS) $(OBJS) 

# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -I$(FOLDER_WFA) -c $< -o $@

###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
	
.PHONY: $(SUBDIRS) 