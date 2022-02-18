###############################################################################
# Rules
###############################################################################
all: generate_dataset

generate_dataset: generate_dataset.c
	$(CC) $(CC_FLAGS) generate_dataset.c -o ../../bin/generate_dataset $(LD_FLAGS) -lm
