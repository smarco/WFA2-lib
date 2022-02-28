###############################################################################
# Definitions
###############################################################################
FOLDER_WFA=..
FOLDER_LIB=../lib
FOLDER_BIN=bin

###############################################################################
# Rules
###############################################################################
LIB_WFA=$(FOLDER_LIB)/libwfa.a
EXAMPLES_C=wfa_adapt \
           wfa_basic \
           wfa_repetead
EXAMPLES_CPP=wfa_bindings \
             wfa_lambda
        
all: setup $(EXAMPLES_C) $(EXAMPLES_CPP)
	
$(EXAMPLES_C): *.c $(LIB_WFA)
	$(CC) $(CC_FLAGS) -L$(FOLDER_LIB) -I$(FOLDER_WFA) wfa_basic.c -o $(FOLDER_BIN)/wfa_basic -lwfa -lm
	$(CC) $(CC_FLAGS) -L$(FOLDER_LIB) -I$(FOLDER_WFA) wfa_adapt.c -o $(FOLDER_BIN)/wfa_adapt -lwfa -lm
	$(CC) $(CC_FLAGS) -L$(FOLDER_LIB) -I$(FOLDER_WFA) wfa_repeated.c -o $(FOLDER_BIN)/wfa_repeated -lwfa -lm

$(EXAMPLES_CPP): *.cpp $(LIB_WFA)
	$(CPP) $(CC_FLAGS) -L$(FOLDER_LIB) -I$(FOLDER_WFA) wfa_bindings.cpp -o $(FOLDER_BIN)/wfa_bindings -lwfacpp
	$(CPP) $(CC_FLAGS) -L$(FOLDER_LIB) -I$(FOLDER_WFA) wfa_lambda.cpp -o $(FOLDER_BIN)/wfa_lambda -lwfacpp
	
setup:
	@mkdir -p $(FOLDER_BIN)

clean: 
	rm -rf $(FOLDER_BIN)
	
###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
	
.PHONY: $(SUBDIRS)