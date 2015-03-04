# Makefile: Makefile for PGAP
# Author: Yan Y. Liu <yanliu@illinois.edu>
# Date: 2014/08/17
# Common changes to be made to run on your platform include:
# CXX: if you use a different mpicc
# CXXFLGS: choose a different one if you want sequential or parallel
# LIBS: change if you use different random number generator
# SPRNG_HOME: environment variable to the path of SPRNG

# to build more than one exec, use NAME=VALUE at command line
# e.g., make clean;make SUARCH_MIC=-mmic HETERO=10 COMM_MODE=async EXEC_SUFFIX=10
# batch generation of execs:
# for m in async sync; do for a in -xhost -mmic; do for h in 1 5 10 15 20; do make clean;make SUARCH_MIC=$a HETERO=$h COMM_MODE=$m EXEC_SUFFIX=$h; done; done; done
# generation of general -mmic and -xhost:
# for m in async sync; do for a in -xhost -mmic; do make clean;make SUARCH_MIC=$a HETERO=10 COMM_MODE=$m ; done; done
# generation of cpu-only version
# for m in async sync; do make clean; make COMM_MODE=$m ; done
# will build ga-async-mmic10
# Architecture-related flags
# Intel MIC arch, default is empty; optoins: -mmic or -xhost
SUARCH_MIC  ?= 
#SUARCH_MIC  ?= -xhost
SUARCH_GPU  ?=
# if compile, set hetero factor
ifeq "$(SUARCH_MIC)" "-mmic"
# mig interval diff b/w host and mic
HETERO      ?= 10
CXXFLGS     += -DHETERO=$(HETERO)
# snd parallelism diff b/w host and mic. default: the same
CXXFLGS     += -DHETERO_BUFFERCAPDIFF=1
endif
EXEC_SUFFIX ?=
# stampede
#SPRNG_HOME  := $(WORK)/sprng2.0$(SUARCH_MIC)

# comm mode: async or sync
COMM_MODE   ?= async
ifeq "$(COMM_MODE)" "sync"
CXXFLGS     += -DPGAMODE_SYNC
endif

SHELL       += -x
#CXX          = gcc
CXX          = mpicc
# Default: Parallel, non-blocking 
CXXFLGS     += -g -Wall -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE -DT_PROFILING
# Sequential 
#CXXFLGS     += -g -Wall -DSPRNG
#CXXFLGS     += -g -Wall -DGSL_SPRNG -DPGAMODE -DPGA_NONBLOCK_MODE
#CXXFLGS     += -g -Wall -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE -DT_PROFILING -DNOIO
##CXXFLGS     += -g -Wall -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE -DT_PROFILING -DDEBUG_COMM
# use MPI_Ibsend(), not robust 'cause buffer policy differs on diff MPIs
#CXXFLGS     += -g -Wall -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE -DDEBUG_COMM -DPGA_USE_IBSEND
SRCC         = c
SRCH         = h
# use static libraries
#STATICLINK   = -static
STATICLINK   =
# standalone obj items
#TARGETS      = log data addr pga gsl-sprng mysprng
TARGETS      = log data addr pga mysprng randseq
#TARGETS      = log data addr
OBJS         = $(TARGETS:=.o)
DEFTARGETS = myrng
DEFS       = $(DEFTARGETS:=.h)
# executables
MAINS        = ga

# external directories, set by env vars
#EXTDIRS      = $(GSL_HOME) $(SPRNG_HOME)
EXTDIRS      = $(SPRNG_HOME)

# include files
INCPATH  = $(EXTDIRS:%=-I%/include)

# library paths
LIBPATH      = $(EXTDIRS:%=-L%/lib)
#LIBS         = -lgsl -lgslcblas -lsprng
LIBS         = -lsprng
#LIBS_DEFAULT = -lm -lpthread
LIBS_DEFAULT = -lm

all: $(MAINS) Makefile

$(MAINS): % : %.$(SRCC) $(DEFS) $(OBJS)
	@$(CXX) $(SUARCH_MIC) $(CXXFLGS) $(STATICLINK) -I. $(INCPATH) -o $@-$(COMM_MODE)$(EXEC_SUFFIX)$(SUARCH_MIC) $< $(OBJS) $(LIBPATH) $(LIBS) $(LIBS_DEFAULT)

# compile
$(DEFS):
$(OBJS): %.o: %.$(SRCC) %.$(SRCH)
	@$(CXX) $(SUARCH_MIC) $(CXXFLGS) -I. $(INCPATH) -o $@ -c $<
# link

# clean
clean:
	@rm -f $(OBJS) 
SUARCH_MICS      := -mmic -xhost
COMM_MODES       := async sync
EXECS_COMM := $(foreach c, $(COMM_MODES), $(MAINS:=-$(c)))
EXECS_COMM_MICS := $(foreach m, $(SUARCH_MICS), $(EXECS_COMM:=$(m)))
cleanall:
	@rm -f $(EXECS_COMM) $(EXECS_COMM_MICS) $(OBJS)
