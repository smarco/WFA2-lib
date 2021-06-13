###############################################################################
# Definitions
###############################################################################
FOLDER_ROOT=..
FOLDER_BUILD=../build

###############################################################################
# Modules
###############################################################################
MODULES=affine2p_dp \
        affine2p_matrix \
        affine2p_penalties

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

CC_XFLAGS=-march=native -fopt-info-vec-optimized

###############################################################################
# Rules
###############################################################################
all: $(OBJS)

$(FOLDER_BUILD)/affine2p_wavefront.o : affine2p_wavefront.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
$(FOLDER_BUILD)/affine2p_wavefront_align.o : affine2p_wavefront_align.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
$(FOLDER_BUILD)/affine2p_wavefront_extend.o : affine2p_wavefront_extend.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	

