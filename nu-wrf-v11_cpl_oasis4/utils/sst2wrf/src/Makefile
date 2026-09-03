#------------------------------------------------------------------------------
# NASA/GSFC, Software Integration & Visualization Office, Code 610.3
#------------------------------------------------------------------------------
#
# Makefile for SST2WRF
#
# AUTHOR:
# Eric Kemp, NGIS
#
# 12/16 : Major mods (C. Cruz)
#------------------------------------------------------------------------------
include $(CONFIG_DIR)/common.mk
BIN_DIR=$(UTILS_DIR)/bin
SCRIPTS_DIR=$(NUWRFDIR)/scripts/python/utils

.PHONY: clean all

LIB = -L../lib -lsst
AR = ar
ARFLAGS = rc
LIBSST = libsst.a

EXE = sst2wrf.x

all:  libsst $(EXE)
	mkdir -p $(BIN_DIR)
	echo 'Installing...'$(EXE)
	cp $(EXE) $(BIN_DIR)
	cp $(SCRIPTS_DIR)/sst2wrf.py $(BIN_DIR)

libsst:
	cd $(CURDIR)/external; \
	make 'FC=$(FC)' \
	     'FFLAGS=$(FFLAGS)' \
	     'AR=$(AR)' \
	     'ARFLAGS=$(ARFLAGS)' \
	     'LIBSST=$(LIBSST)'

clean:
	make -C external clean
	rm -f *.[ox] *.mod ../lib/libsst.a


$(EXE): TimeUtil_mod.o FieldSST_RSS_mod.o \
   OutputHours_mod.o FieldWPS_mod.o FileWPS_mod.o sst2wrf.o
	$(FC) $(FFLAGS) TimeUtil_mod.o \
	FieldSST_RSS_mod.o OutputHours_mod.o FieldWPS_mod.o FileWPS_mod.o \
	sst2wrf.o \
	-o $(EXE) -L$(CURDIR)/external -lsst -L$(UTILS_DIR)/shared -lshared

testTimeUtil: TimeUtil_mod.o testTimeUtil.o
	$(FC) $(FFLAGS) testTimeUtil.o TimeUtil_mod.o -o testTimeUtil

# Source code dependencies

sst2wrf.o: sst2wrf.f90 TimeUtil_mod.o FieldSST_RSS_mod.o \
	   OutputHours_mod.o FieldWPS_mod.o FileWPS_mod.o
	$(FC) $(FFLAGS) -I$(UTILS_DIR)/shared -c sst2wrf.f90

testTimeUtil.o: testTimeUtil.f90
	$(FC) $(FFLAGS) -c testTimeUtil.f90

TimeUtil_mod.o: TimeUtil_mod.f90
	$(FC) $(FFLAGS) -c TimeUtil_mod.f90

FieldSST_RSS_mod.o: FieldSST_RSS_mod.f90 TimeUtil_mod.o
	$(FC) $(FFLAGS) -c FieldSST_RSS_mod.f90

OutputHours_mod.o:  OutputHours_mod.f90 \
                    FieldSST_RSS_mod.o 
	$(FC) $(FFLAGS) -c OutputHours_mod.f90

FieldWPS_mod.o: FieldWPS_mod.f90
	$(FC) $(FFLAGS) -c FieldWPS_mod.f90

FileWPS_mod.o: FileWPS_mod.f90 FieldWPS_mod.o
	$(FC) $(FFLAGS) -c FileWPS_mod.f90

