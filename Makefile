# Diretories
MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
EXAMPLEDIR = $(MDIR)/example

OUTPUTDIR = $(MDIR)/output

# Directory of splinter
SPLINTERDIR = $(MDIR)/external/splinter

# Directory of Cuba
CUBADIR = $(MDIR)/external/Cuba-4.2.2

splinter: libsplinter-static-4-0.a
libsplinter-static-4-0.a:
	(cd $(SPLINTERDIR) && if ! [ -e build ]; then mkdir build ; fi; cd build CXX=g++ && cmake .. && make)

cuba: libcuba.a
libcuba.a:
	(cd $(CUBADIR) && ./configure && make)

.base-build:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; fi;
	echo "Do not remove this file or directory manually..." > build/.base-build

.base-output:
	if ! [ -e $(OUTPUTDIR) ]; then mkdir $(OUTPUTDIR) $(OUTPUTDIR)/spec $(OUTPUTDIR)/cov $(OUTPUTDIR)/fish $(OUTPUTDIR)/sample ; fi;
	echo "Do not remove this file or directory manually..." > output/.base-output

vpath %.c src
vpath %.o build

vpath .base-build build
vpath .base-output output

vpath example.c $(EXAMPLEDIR)

vpath libsplinter-static-4-0.a $(SPLINTERDIR)/build
vpath libcuba.a $(CUBADIR)

# C compiler:
CC = gcc

# Tool for creating static libraries:
AR = ar rv

# Error flags
ERRFLAG = -Wall -Wconversion -Wextra

# Optimization flag
OPTFLAG = -O3

# OpenMP flag
OMPFLAG = -fopenmp

# Other compilation flags
CCFLAG = -g -fPIC
LDFLAG = -g -fPIC

# Libraries
LIBS = -L$(FISHERLSSDIR) -L$(CUBADIR) -I$(CUBADIR) -lfisherlss -L$(SPLINTERDIR) -I$(SPLINTERDIR) -lsplinter-4-0 -lcuba -lm -lquadmath -lgsl -llapack -lblas

# Pass current working directory to the code
FISHERLSSDIR ?= $(MDIR)
CCFLAG += -D__FISHERLSSDIR__='"$(FISHERLSSDIR)"'

# Where to find include files
INCLUDES = -I../include
INCLUDESSPLINTER = -I$(SPLINTERDIR)/include/cinterface
HEADERFILES = $(wildcard ./include/*.h)

# Create object files
%.o:  %.c .base-build .base-output libcuba.a libsplinter-static-4-0.a $(HEADERFILES)
	cd $(WRKDIR);$(CC) $(ERRFLAG) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) $(INCLUDESSPLINTER) -I$(SPLINTERDIR) -I$(CUBADIR) -c ../$< -o $*.o

SOURCE = flss.o spectra.o covariance.o fisher.o print.o fiducials.o integrand.o average.o integrand-average.o kernels.o integrate.o interpolate.o sample.o shape.o dat.o mat.o misc.o

# Source files
C_SOURCE = $(addprefix src/, $(addsuffix .c,$(basename $(SOURCE))))

# All files
C_ALL = $(SOURCE)
H_ALL = $(addprefix include/, common.h $(addsuffix .h, $(basename $(notdir $(C_ALL)))))

# Create static library
all: fisher-lss example

fisher-lss: libfisherlss.a

libfisherlss.a: $(SOURCE)
	$(AR) $@ $(addprefix build/, $(SOURCE))

example: example.c libfisherlss.a
	$(CC) $(ERRFLAG) $(EXAMPLEDIR)/example.c $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $(EXAMPLEDIR)/example.o $(LIBS) $(INCLUDESSPLINTER)

# Clean
clean: clean-build clean-output clean-cuba clean-splinter

clean-cuba:
	(cd $(CUBADIR) && make distclean)

clean-splinter:
	(cd $(SPLINTERDIR)/build && make clean && rm -rf ../build)

clean-build:
	rm -rf $(WRKDIR);
	rm -f libfisherlss.a
	rm -f $(EXAMPLEDIR)/example.o

clean-output:
	rm -rf $(OUTPUTDIR);
