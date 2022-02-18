###############################################################################
# Definitions
###############################################################################
FOLDER_WFA=../../../..
FOLDER_ROOT=../..
FOLDER_BUILD=$(FOLDER_ROOT)/build

###############################################################################
# Modules
###############################################################################
MODULES=benchmark_bitpal \
        benchmark_gaba \
        benchmark_ksw2 \
        benchmark_parasail \
        benchmark_seqan
        
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

LI_FLAGS=-I$(FOLDER_ROOT) -I$(FOLDER_WFA)

###############################################################################
# Rules
###############################################################################
all: $(OBJS) $(FOLDER_BUILD)/benchmark_seqan_bindings.o

$(FOLDER_BUILD)/benchmark_gaba.o : benchmark_gaba.c
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c $< -o $@

INCLUDE_SEQAN=$(FOLDER_ROOT)/external/seqan/include/
$(FOLDER_BUILD)/benchmark_seqan_bindings.o : benchmark_seqan_bindings.cpp
	$(CPP) -std=c++14 $(CC_FLAGS) $(LI_FLAGS) -I$(INCLUDE_SEQAN) -c $< -o $@

# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c $< -o $@
	
