#################################################################################################
# Makefile for CHAPSim2, by Wei Wang, July 2021                                                 #
# Updated: Mar 2026 - robust stop-on-error logging with tee, portable for bash/zsh environments #
# Usage:                                                                                        #
#       make all              to make all files with default optimization                       #
#       make cfg=gnu-g        to debug for gfortran compiler (no DEBUG_STEPS)                   #
#       make cfg=gnu-debug    to debug for gfortran compiler (with DEBUG_STEPS)                 #
#       make cfg=gnu-o3       to optimize with -O3                                              #
#       make cfg=intel        to debug for intel compiler                                       #
#       make cfg=cray         to make for cray                                                  #
#################################################################################################

#################################################################################################
# Makefile for CHAPSim2, by Wei Wang, July 2021                                                 #
# Updated: Mar 2026 - POSIX sh compatible, stop on error, no bash dependency                    #
#################################################################################################

.SUFFIXES:
SHELL := /bin/sh
.SHELLFLAGS := -ec
.DELETE_ON_ERROR:

# Program name and directories
PROGRAM = CHAPSim
DIR_SRC = ../src
DIR_BIN = ../bin
DIR_OBJ = ../obj
DIR_BLD = ../build
DIR_LIB = ../lib/fishpack4.1

# Default configuration if 'cfg' is not specified
cfg ?= default

# Detect architecture and platform
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# Set architecture-specific flags
ARCH_FLAGS :=
ifeq ($(UNAME_S),Darwin)
    ifeq ($(UNAME_M),arm64)
        ARCH_FLAGS := -arch arm64
        $(info Detected: Apple Silicon ($(UNAME_M)))
    else ifeq ($(UNAME_M),x86_64)
        ARCH_FLAGS := -arch x86_64
        $(info Detected: Intel Mac ($(UNAME_M)))
    endif
else ifeq ($(UNAME_S),Linux)
    $(info Detected: Linux ($(UNAME_M)))
else
    $(info Detected: $(UNAME_S) ($(UNAME_M)))
endif

# Common Fortran Compiler Flags
COMMON_FOPTS = $(ARCH_FLAGS) -fallow-argument-mismatch -Wno-unused -Wno-unused-dummy-argument \
               -Wno-maybe-uninitialized -fdiagnostics-color=always -Wuninitialized \
               -finit-real=snan -finit-integer=-999999 -ffpe-trap=invalid,zero,overflow \
               -Wall -Wextra -fbacktrace -fbounds-check -fimplicit-none -cpp -fopenmp \
               -finit-local-zero -ffree-line-length-512 -DDOUBLE_PREC \
               -fdefault-real-8 -fdefault-double-8

FOPTS_WARNING = -Wunused-parameter -Wunused-dummy-argument

# Compiler and flags based on configuration
ifeq ($(cfg),gnu-g)
  FC     = mpif90
  FOPTS  = -g -O0 $(COMMON_FOPTS) $(FOPTS_WARNING)
  FDEBG  =
else ifeq ($(cfg),gnu-debug)
  FC     = mpif90
  FOPTS  = -g -O0 $(COMMON_FOPTS) $(FOPTS_WARNING)
  FDEBG  = -DDEBUG_STEPS
else ifeq ($(cfg),intel)
  FC     = mpiifort
  FOPTS  = -g -assume ieee_fpe_flags -check all -check bounds -check uninit -debug all \
           -fp-stack-check fpe0 -fpe3 -fpe-all=3 -ftrapuv -ftz -warn all,nounused -fpp \
           -convert big_endian
  FDEBG  = -DDEBUG_STEPS -DDEBUG_FFT
else ifeq ($(cfg),cray)
  FC     = ftn
  FOPTS  = -cpp
  FDEBG  = -DDEBUG_STEPS -DDEBUG_FFT
else ifeq ($(cfg),gnu-o3)
  FC     = mpif90
  FOPTS  = -O3 $(COMMON_FOPTS)
  FDEBG  =
else ifeq ($(cfg),nvhpc-cpu)
  FC     = mpif90
  FOPTS  = -O3 -Kieee -Minfo=accel -target=multicore -acc -cpp -DDOUBLE_PREC
  FDEBG  =
else
  FC     = mpif90
  FOPTS  = -O3 $(COMMON_FOPTS)
  FDEBG  =
endif

# HDF5 paths
HDF5_INC =
HDF5_LIB =

# FFT library and include paths
INCLUDE = -I../lib/2decomp-fft/build/opt/include $(HDF5_INC)

# Check for library file directly
LIB_FILE_64 := $(wildcard ../lib/2decomp-fft/build/opt/lib64/libdecomp2d.a)
LIB_FILE    := $(wildcard ../lib/2decomp-fft/build/opt/lib/libdecomp2d.a)

# Determine which library to use
ifneq ($(LIB_FILE_64),)
    DECOMP_LIB = ../lib/2decomp-fft/build/opt/lib64/libdecomp2d.a
    $(info Using library: lib64/libdecomp2d.a)
else ifneq ($(LIB_FILE),)
    DECOMP_LIB = ../lib/2decomp-fft/build/opt/lib/libdecomp2d.a
    $(info Using library: lib/libdecomp2d.a)
else
    $(error ❌ Error: libdecomp2d.a not found in lib or lib64. Run library build first!)
endif

# Verify library is readable
LIB_CHECK := $(shell ar -t $(DECOMP_LIB) > /dev/null 2>&1 && echo OK || echo FAIL)
ifneq ($(LIB_CHECK),OK)
    $(error ❌ Error: Library $(DECOMP_LIB) is corrupted or unreadable. Rebuild the library!)
endif

# Use the library file directly instead of -L/-l flags
LIBS = $(DECOMP_LIB) $(HDF5_LIB)

# Object file lists
OBJS1 = $(addprefix $(DIR_OBJ)/, \
      modules.o \
      tools_general.o \
      input_thermo.o \
      bc_dirichlet.o \
      bc_ndomain_interior.o \
      bc_general.o \
      input_general.o \
      basics_algorithms.o \
      basics_operations2.o \
      tools_solver.o \
      para_conversion.o \
      geometry.o \
      io_tools.o \
      io_monitor.o \
      io_visulisation.o \
      post_statistics.o \
      domain_decomposition.o \
      poisson_fishpack.o \
      poisson_1stderivcomp_fft2d.o \
      poisson_interface.o \
			eq_continuity.o \
      bc_convective_outlet.o \
      eq_energy.o \
      eq_mhd.o \
      io_restart.o \
      eq_momentum2.o \
      initialisation.o \
      test_algrithms.o \
      chapsim.o)

OBJS2 = $(addprefix $(DIR_OBJ)/, fftpack.o)

# Compilation flags for FFT library (Fortran 77 files)
F77_OPTS = $(ARCH_FLAGS) -O2 -w -fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch

# Total number of object files for progress display
TOTAL_OBJS := $(words $(OBJS1) $(OBJS2))
COUNTER = 0

# Function to display compilation progress
define progress
    $(eval COUNTER := $(shell expr $(COUNTER) + 1))
    $(eval PERCENT := $(shell expr $(COUNTER) \* 100 / $(TOTAL_OBJS)))
    @printf "(%2d/%2d) [%3d%%] Compiling %s\n" $(COUNTER) $(TOTAL_OBJS) $(PERCENT) $(notdir $<)
endef

# Default targets
default: all

all: $(DIR_BIN)/$(PROGRAM)
	@echo ""
	@echo "========================================================================="
	@echo "  CHAPSim2 Build System"
	@echo "  Configuration: $(cfg)"
	@echo "========================================================================="
	@echo ""

# Main target for building the program
$(DIR_BIN)/$(PROGRAM): $(OBJS1) $(OBJS2)
	@echo ""
	@echo "========================================================================="
	@echo "  Linking $(PROGRAM)..."
	@echo "  Architecture: $(UNAME_M) on $(UNAME_S)"
	@echo "  Compiler: $(FC)"
	@echo "  Library: $(DECOMP_LIB)"
	@echo "========================================================================="
	@mkdir -p $(@D)
	@$(FC) $(FOPTS) $(FDEBG) $(ARCH_FLAGS) -o $@ $^ $(LIBS) >> build.log 2>&1 || \
	  { echo "❌ Link failed. Last 40 lines of build.log:"; tail -40 build.log; exit 1; }
	@echo ""
	@echo "========================================================================="
	@echo "  ✅ Successfully compiled $(PROGRAM)"
	@echo "  Binary location: $@"
	@echo "========================================================================="
	@echo ""

# Rule to compile Fortran 90 files
$(DIR_OBJ)/%.o: $(DIR_SRC)/%.f90
	@mkdir -p $(@D)
	$(call progress)
	@$(FC) $(INCLUDE) $(FOPTS) $(FDEBG) -c $< -o $@ >> build.log 2>&1 || \
	  { echo "❌ Compile failed for $<. Last 40 lines of build.log:"; tail -40 build.log; exit 1; }

# Rule to compile Fortran 77 files
$(DIR_OBJ)/%.o: $(DIR_LIB)/%.f
	@mkdir -p $(@D)
	$(call progress)
	@$(FC) $(F77_OPTS) -c $< -o $@ >> build.log 2>&1 || \
	  { echo "❌ Compile failed for $<. Last 40 lines of build.log:"; tail -40 build.log; exit 1; }

# Verify library target
verify-lib:
	@echo "Verifying 2decomp library..."
	@if [ -f "$(DECOMP_LIB)" ]; then \
		echo "✅ Library found: $(DECOMP_LIB)"; \
		echo "File type:"; \
		file $(DECOMP_LIB); \
		echo ""; \
		echo "Archive contents (first 10):"; \
		ar -t $(DECOMP_LIB) | head -10; \
		echo "..."; \
		echo "Total objects: $$(ar -t $(DECOMP_LIB) | wc -l)"; \
		if ar -t $(DECOMP_LIB) | grep -E "^/$$|^//" > /dev/null 2>&1; then \
			echo ""; \
			echo "⚠️  WARNING: Found suspicious '/' entries!"; \
			ar -t $(DECOMP_LIB) | grep -E "^/$$|^//"; \
		else \
			echo ""; \
			echo "✅ No suspicious entries found"; \
		fi; \
	else \
		echo "❌ Library not found: $(DECOMP_LIB)"; \
		exit 1; \
	fi

# Clean target
clean:
	@echo "Cleaning built files and directories..."
	@echo ""
	@rm -rf $(DIR_OBJ) $(DIR_BIN)/$(PROGRAM) build.log
	@find $(DIR_BLD) -name "*.mod" -delete -o -name "*.smod" -delete 2>/dev/null || true
	@echo "✅ Clean complete"
	@echo ""

# Help target
help:
	@echo "CHAPSim2 Makefile - Available targets:"
	@echo ""
	@echo "  make all              - Build with default settings (O3 optimization)"
	@echo "  make cfg=gnu-o3       - Build with -O3 optimization"
	@echo "  make cfg=gnu-g        - Build with debug symbols (-g -O0)"
	@echo "  make cfg=gnu-debug    - Build with debug + DEBUG_STEPS"
	@echo "  make cfg=intel        - Build with Intel compiler"
	@echo "  make cfg=cray         - Build with Cray compiler"
	@echo "  make cfg=nvhpc-cpu    - Build with NVHPC CPU config"
	@echo "  make clean            - Remove all build artifacts"
	@echo "  make verify-lib       - Verify 2decomp library integrity"
	@echo "  make help             - Show this help message"
	@echo ""
	@echo "Current system:"
	@echo "  OS: $(UNAME_S)"
	@echo "  Architecture: $(UNAME_M)"
	@echo "  Library: $(DECOMP_LIB)"
	@echo ""

# Phony targets
.PHONY: all clean default verify-lib help