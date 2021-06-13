###############################################################################
# Definitions
###############################################################################
FOLDER_ROOT=..
FOLDER_BUILD=../build
FOLDER_BIN=../bin

###############################################################################
# Tools
###############################################################################
TOOLS=generate_dataset align_benchmark
TOOLS_SRC=$(addsuffix .c, $(TOOLS))

###############################################################################
# Rules
###############################################################################
OBJS=$(FOLDER_BUILD)/*.o

all: LIBS+=$(LD_FLAGS)
all: FLAGS=$(CC_FLAGS)
all: $(TOOLS)

align_benchmark: $(FOLDER_BUILD)/*.o align_benchmark.c
	$(CC) $(FLAGS) -I$(FOLDER_ROOT) align_benchmark.c $(OBJS) -o $(FOLDER_BIN)/align_benchmark $(LIBS)
	
generate_dataset: generate_dataset.c
	$(CC) $(CC_FLAGS) generate_dataset.c -o $(FOLDER_BIN)/generate_dataset $(LD_FLAGS)
