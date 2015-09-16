# Location of the CUDA Toolkit
CUDA_PATH ?= /sw/cuda/7.0

OSUPPER = $(shell uname -s 2>/dev/null | tr "[:lower:]" "[:upper:]")
OSLOWER = $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")

OS_SIZE    = $(shell uname -m | sed -e "s/x86_64/64/" -e "s/armv7l/32/" -e "s/aarch64/64/")
OS_ARCH    = $(shell uname -m)
ARCH_FLAGS =

DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))
ifneq ($(DARWIN),)
	XCODE_GE_5 = $(shell expr `xcodebuild -version | grep -i xcode | awk '{print $$2}' | cut -d'.' -f1` \>= 5)
endif

# Take command line flags that override any of these settings
ifeq ($(x86_64),1)
	OS_SIZE = 64
	OS_ARCH = x86_64
endif
ifeq ($(ARMv7),1)
	OS_SIZE    = 32
	OS_ARCH    = armv7l
	ARCH_FLAGS = -target-cpu-arch ARM
endif
ifeq ($(aarch64),1)
	OS_SIZE    = 64
	OS_ARCH    = aarch64
	ARCH_FLAGS = -target-cpu-arch ARM
endif

# Common binaries
ifneq ($(DARWIN),)
ifeq ($(XCODE_GE_5),1)
  GCC ?= clang
else
  GCC ?= g++
endif
else
ifeq ($(ARMv7),1)
  GCC ?= arm-linux-gnueabihf-g++
else
  GCC ?= g++
endif
endif
# NVCC := nvcc -ccbin $(GCC)
NVCC := $(CUDA_PATH)/bin/nvcc -ccbin $(GCC)

# internal flags
NVCCFLAGS   := -m${OS_SIZE} ${ARCH_FLAGS}
CCFLAGS     :=
LDFLAGS     :=

# Extra user flags
EXTRA_NVCCFLAGS   ?=
EXTRA_LDFLAGS     ?=
EXTRA_CCFLAGS     ?=

# OS-specific build flags
ifneq ($(DARWIN),)
  LDFLAGS += -rpath $(CUDA_PATH)/lib
  CCFLAGS += -arch $(OS_ARCH)
else
  ifeq ($(OS_ARCH),armv7l)
    ifeq ($(abi),androideabi)
      NVCCFLAGS += -target-os-variant Android
    else
      ifeq ($(abi),gnueabi)
        CCFLAGS += -mfloat-abi=softfp
      else
        # default to gnueabihf
        override abi := gnueabihf
        LDFLAGS += --dynamic-linker=/lib/ld-linux-armhf.so.3
        CCFLAGS += -mfloat-abi=hard
      endif
    endif
  endif
endif

ifeq ($(ARMv7),1)
ifneq ($(TARGET_FS),)
GCCVERSIONLTEQ46 := $(shell expr `$(GCC) -dumpversion` \<= 4.6)
ifeq ($(GCCVERSIONLTEQ46),1)
CCFLAGS += --sysroot=$(TARGET_FS)
endif
LDFLAGS += --sysroot=$(TARGET_FS)
LDFLAGS += -rpath-link=$(TARGET_FS)/lib
LDFLAGS += -rpath-link=$(TARGET_FS)/usr/lib
LDFLAGS += -rpath-link=$(TARGET_FS)/usr/lib/arm-linux-$(abi)
endif
endif

ALL_CCFLAGS :=
ALL_CCFLAGS += $(NVCCFLAGS)
ALL_CCFLAGS += $(EXTRA_NVCCFLAGS)
ALL_CCFLAGS += $(addprefix -Xcompiler ,$(CCFLAGS))
ALL_CCFLAGS += $(addprefix -Xcompiler ,$(EXTRA_CCFLAGS))

ALL_LDFLAGS :=
ALL_LDFLAGS += $(ALL_CCFLAGS)
ALL_LDFLAGS += $(addprefix -Xlinker ,$(LDFLAGS))
ALL_LDFLAGS += $(addprefix -Xlinker ,$(EXTRA_LDFLAGS))

MPI_CCFLAGS :=
MPI_CCFLAGS += $(CCFLAGS)
MPI_CCFLAGS += $(EXTRA_CCFLAGS)

MPI_LDFLAGS :=
MPI_LDFLAGS += $(addprefix -Xlinker ,$(LDFLAGS))
MPI_LDFLAGS += $(addprefix -Xlinker ,$(EXTRA_LDFLAGS))

# Gdal
GDALDIR = /sw/geosoft/gdal-1.11.2-fgdb
GDALLIB := -L$(GDALDIR)/lib -lgdal

# Common includes and paths for CUDA
INCLUDES  := -I../../common/inc
LIBRARIES :=

SAMPLE_ENABLED := 1

# MPI check and binaries
MPICXX ?= $(shell which mpicxx 2>/dev/null)

ifneq ($(shell uname -m | sed -e "s/i386/i686/"), ${OS_ARCH})
      $(info -----------------------------------------------------------------------------------------------)
      $(info WARNING - attempting to detect 32-bit MPI compiler.)
      MPICXX := $(shell echo $(MPICXX) | sed -e "s/64//")
endif

ifeq ($(MPICXX),)
      $(info -----------------------------------------------------------------------------------------------)
      $(info WARNING - No MPI compiler found.)
      $(info -----------------------------------------------------------------------------------------------)
      $(info   CUDA Sample "simpleMPI" cannot be built without an MPI Compiler.)
      $(info   This will be a dry-run of the Makefile.)
      $(info   For more information on how to set up your environment to build and run this )
      $(info   sample, please refer the CUDA Samples documentation and release notes)
      $(info -----------------------------------------------------------------------------------------------)
      MPICXX=mpicxx
      SAMPLE_ENABLED := 0
else
      MPI_GCC := $(shell $(MPICXX) -v 2>&1 | grep gcc | wc -l | tr -d ' ')
ifeq ($(MPI_GCC),0)
      MPI_CCFLAGS += -stdlib=libstdc++
      MPI_LDFLAGS += -stdlib=libstdc++
endif
endif

# Gencode arguments
ifeq ($(OS_ARCH),armv7l)
SMS ?= 20 30 32 35 37 50
else
SMS ?= 20 30 35 37 50
endif

ifeq ($(SMS),)
$(info >>> WARNING - no SM architectures have been specified - waiving sample <<<)
SAMPLE_ENABLED := 0
endif

ifeq ($(GENCODE_FLAGS),)
# Generate SASS code for each SM architecture listed in $(SMS)
$(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))

# Generate PTX code from the highest SM architecture in $(SMS) to guarantee forward-compatibility
HIGHEST_SM := $(lastword $(sort $(SMS)))
ifneq ($(HIGHEST_SM),)
GENCODE_FLAGS += -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)
endif
endif

LIBSIZE :=
ifeq ($(DARWIN),)
ifeq ($(OS_SIZE),64)
LIBSIZE := 64
endif
endif

LIBRARIES += -L$(CUDA_PATH)/lib$(LIBSIZE) -lcudart

ifeq ($(SAMPLE_ENABLED),0)
EXEC ?= @echo "[@]"
endif

# Target rules
all: MultiGPUKDE

multiGPUKDE.o: multiGPUKDE.c
	$(MPICXX) $(INCLUDES) $(MPI_CCFLAGS) -o $@ -c $<

io.o: io.c
	$(GCC) -I. -I$(GDALDIR)/include -o $@ -c $<

decompose.o: decompose.c
	$(GCC) -o $@ -c $<

kde.o: kde.cu
	${NVCC} $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

MultiGPUKDE: multiGPUKDE.o io.o decompose.o kde.o
	$(MPICXX) $(MPI_LDFLAGS) -o ../$@ $+ $(LIBRARIES) $(GDALLIB)

clean:
	rm -f multiGPUKDE.o io.o decompose.o kde.o ../MultiGPUKDE
test:
	mpirun -np 1 ../MultiGPUKDE ../data/SampleData.csv ../testResult.tif -2380000 2280000 -1520000 1420000 2000 20000 1 1 


