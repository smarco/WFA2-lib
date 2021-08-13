###############################################################################
# Definitions
###############################################################################
FOLDER_ROOT=..
FOLDER_BUILD=../build

###############################################################################
# Modules
###############################################################################
MODULES=wavefront_align \
        wavefront_aligner \
        wavefront_attributes \
        wavefront_backtrace_buffer \
        wavefront_backtrace \
        wavefront_components \
        wavefront_compute_affine \
        wavefront_compute_affine2p \
        wavefront_compute \
        wavefront_display \
        wavefront_extend \
        wavefront_pcigar \
        wavefront_penalties \
        wavefront_plot \
        wavefront_reduction \
        wavefront_slab \
        wavefront

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

CC_XFLAGS=-march=native -fopt-info-vec-optimized

###############################################################################
# Rules
###############################################################################
all: $(OBJS)

$(FOLDER_BUILD)/wavefront_compute_affine.o : wavefront_compute_affine.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/wavefront_compute_affine2p.o : wavefront_compute_affine2p.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
$(FOLDER_BUILD)/wavefront_extend.o : wavefront_extend.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	

