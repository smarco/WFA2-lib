###############################################################################
# Compile rules
###############################################################################
BUILD_RULES=bitpal_build \
            blockaligner_build \
            daligner_wave_build \
            edlib_build \
            ksw2_build \
            gaba_build \
            parasail_build \
            scrooge_build
            
CLEAN_RULES=blockaligner_clean \
            edlib_clean \
            gaba_clean \
            parasail_clean

all: prepare $(BUILD_RULES)

prepare:
	-mkdir build

clean: $(CLEAN_RULES)
	-rm -rf build

###############################################################################
# Build rules
###############################################################################
bitpal_build:
	-$(MAKE) --directory=bitpal all

blockaligner_build:
	-$(MAKE) --directory=block-aligner/c all

daligner_wave_build:
	-$(MAKE) --directory=daligner_wave all

edlib_build:
	-cd edlib; ./build.sh

ksw2_build:
	-$(MAKE) --directory=ksw2 all
	-mv ksw2/kalloc.o ksw2/ksw2_gg.o ksw2/ksw2_gg2.o ksw2/ksw2_gg2_sse.o \
	    ksw2/ksw2_extz.o ksw2/ksw2_extz2_sse.o ksw2/ksw2_extd.o ksw2/ksw2_extd2_sse.o \
	    ksw2/ksw2_extf2_sse.o ksw2/ksw2_exts2_sse.o build

gaba_build:
	-$(MAKE) --directory=libgaba all

parasail_build:
	-cd parasail; ./build.sh
	
scrooge_build:
	-$(MAKE) --directory=scrooge all

###############################################################################
# Clean rules
###############################################################################
blockaligner_clean:
	-$(MAKE) --directory=block-aligner clean

edlib_clean:
	-cd edlib; ./clean.sh;

gaba_clean:
	-$(MAKE) --directory=libgaba clean

parasail_clean:
	-cd parasail; ./clean.sh;



