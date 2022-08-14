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
        wavefront_backtrace_offload \
        wavefront_backtrace \
        wavefront_bialign \
        wavefront_bialigner \
        wavefront_components \
        wavefront_compute_affine \
        wavefront_compute_affine2p \
        wavefront_compute_edit \
        wavefront_compute_linear \
        wavefront_compute \
        wavefront_debug \
        wavefront_display \
        wavefront_extend \
        wavefront_heuristic \
        wavefront_pcigar \
        wavefront_penalties \
        wavefront_plot \
        wavefront_slab \
        wavefront_unialign \
        wavefront

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

ifeq ($(BUILD_WFA_PARALLEL),1) 
PFLAGS=-DWFA_PARALLEL -fopenmp
endif
#CC_XFLAGS+=-fopt-info-vec-optimized

###############################################################################
# Rules
###############################################################################
all: $(OBJS)

# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(CC_FLAGS) $(PFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	

