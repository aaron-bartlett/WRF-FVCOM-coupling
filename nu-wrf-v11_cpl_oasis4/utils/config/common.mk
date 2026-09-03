ifndef ARCH             # Architecture, e.g., Linux
  ARCH := $(shell uname -s)
endif
ifndef MACH             # Hardware type, e.g., x86_64
  MACH := $(shell uname -m)
endif
ifndef SITE             # Site name, e.g., discover
  SITE := $(shell uname -n)
endif

ifndef NUWRF_FC
  FC = ifort
else
  FC = $(NUWRF_FC)
endif

# -----------------
#     Libraries
#
# MCB 9/19/2017 - Fix LIB_CURL script problem
# LIB_CURL had the following line:
#	$(BASELIBS)/curl/lib/lib4(lib).a) )
#                           ^
#                           |
#   '4' instead of a '$'----|
#
# -----------------
INC_CURL = $(BASELIBS)/curl/include/curl
#LIB_CURL = $(wildcard $(foreach lib,curl,\
#           $(BASELIBS)/curl/lib/lib$(lib).a) )
LIB_CURL = $(BASELIBS)/curl/lib

INC_HDF4 = $(BASELIBS)/hdf4/include
LIB_HDF4 = $(wildcard $(foreach lib,mfhdf df,\
           $(BASELIBS)/hdf4/lib/lib$(lib).a) )
LIB_HDF4 += $(wildcard $(foreach lib,jpeg,\
           $(BASELIBS)/jpeg/lib/lib$(lib).a) )

INC_HDF5 = $(BASELIBS)/hdf5/include
LIB_HDF5 = $(wildcard $(foreach lib,hdf5_hl hdf5hl_fortran hdf5_fortran hdf5,\
           $(BASELIBS)/hdf5/lib/lib$(lib).a) )

INC_NETCDF = $(BASELIBS)/netcdf4/include
LIB_NETCDF = $(wildcard $(foreach lib,netcdff netcdf,\
           $(BASELIBS)/netcdf4/lib/lib$(lib).a) )

INC_ZLIB = $(BASELIBS)/zlib/include
LIB_ZLIB = $(BASELIBS)/zlib/lib

LIBS = -L$(LIB_CURL) -lcurl $(LIB_NETCDF) $(LIB_HDF5) -L$(LIB_ZLIB) -lz \
        $(LIB_HDF4) -L$(UTILS_DIR)/shared -lshared -lz -ldl 

INCS = -I$(INC_CURL) -I$(INC_HDF5) -I$(INC_NETCDF) -I$(INC_ZLIB)

FFLAGS = $(INCS) -I$(UTILS_DIR)/shared

.SUFFIXES:
.SUFFIXES: .f90 .F90

ifeq ($(ARCH), Linux)

ifeq ($(FC), ifort)
  CC = icc
  CFLAGS = -O2
  FFLAGS_DEBUG = -g -CB -fpe0 -check uninit -ftrapuv -traceback -fpp
  FFLAGS_OPT = -O3 -ftz -ip -fp-model strict -w -fpp
endif
ifeq ($(FC), gfortran)
  CC = gcc
  CFLAGS = -O2
  FFLAGS_DEBUG = -g -fbounds-check -fcheck-array-temporaries \
          -ffpe-trap=invalid,zero,overflow -fbacktrace
  FFLAGS_OPT = -O2 -funroll-loops -ftree-vectorize -cpp -fomit-frame-pointer 
endif

ifndef UTILS_DEBUG
FFLAGS += $(FFLAGS_OPT)
else
FFLAGS += $(FFLAGS_DEBUG)
endif

endif

# WORKAROUNDS FOR SOME COMPONENTS
# -------------------------------

# Temporary workaround to deal with top level build.sh scripts
UTILS_DIR ?= $(NUWRF_DIR)/utils

# This is needed for plot_chem
export NCARG_ROOT := $(BASELIBS)/ncarg
export PATH := $(NCARG_ROOT)/bin:$(PATH)

# prep_chem_sources
ifeq ($(FC), ifort)
MAKEPSC_OPT=opt.intel.nuwrf
endif
ifeq ($(FC), gfortran)
MAKEPSC_OPT=opt.gfortran.nuwrf
endif
CHEM_TYPE=RADM_WRF_FIM
AER_TYPE=SIMPLE

