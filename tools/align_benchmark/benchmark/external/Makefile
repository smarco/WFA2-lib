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
        benchmark_seqan_bindings \
        benchmark_seqan \
        benchmark_wfalm_bindings \
        benchmark_wfalm
        
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

LI_FLAGS=-I$(FOLDER_ROOT) -I$(FOLDER_WFA)

###############################################################################
# Rules
###############################################################################
all: $(OBJS) $(FOLDER_BUILD)/benchmark_seqan_bindings.o

$(FOLDER_BUILD)/benchmark_gaba.o : benchmark_gaba.c
	$(CC) -march=native $(CC_FLAGS) $(LI_FLAGS) -c $< -o $@

$(FOLDER_BUILD)/benchmark_seqan_bindings.o : benchmark_seqan_bindings.cpp
	$(CPP) -std=c++14 $(CC_FLAGS) $(LI_FLAGS) -I$(FOLDER_ROOT)/external/seqan/include/ -c $< -o $@
	
$(FOLDER_BUILD)/benchmark_wfalm_bindings.o : benchmark_wfalm_bindings.cpp
	$(CPP) -std=c++14 $(CC_FLAGS) -Wno-sign-compare $(LI_FLAGS) -c $< -o $@

# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c $< -o $@
	
