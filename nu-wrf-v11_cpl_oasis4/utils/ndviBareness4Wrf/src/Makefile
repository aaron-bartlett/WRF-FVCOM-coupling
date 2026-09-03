#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Sciences and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# Makefile for NDVI processing code
#
# AUTHOR:
# Eric Kemp, NASA GSFC/SSAI
#
# 12/16 : Major mods (C. Cruz)
#------------------------------------------------------------------------------
include $(CONFIG_DIR)/common.mk
BIN_DIR=$(UTILS_DIR)/bin

# Additional ndviBareness4Wrf dependencies
GDAL = $(BASELIBS)/gdal
FORTRANGIS = $(BASELIBS)/fortrangis
INCDIRS_NDVI = -I$(GDAL)/include -I$(FORTRANGIS)/include
LIBDIRS_NDVI= -L$(FORTRANGIS)/lib -L$(GDAL)/lib \
	-lfortrangis -lfortranc -lgdal -lstdc++

# Also, modify FFLAGS and LIBS
ifeq ($(FC),ifort)
   FFLAGS_NDVI = $(FFLAGS) -convert big_endian -assume byterecl
endif
ifeq ($(FC),gfortran)
   FFLAGS_NDVI = $(FFLAGS) -fconvert=big-endian -frecord-marker=4
   LIBDIRS_NDVI += -lpthread -ldl
endif

.PHONY: clean all

EXE = ndviBareness4Wrf.x

all: $(EXE)
	mkdir -p $(BIN_DIR)
	echo 'Installing...'$(EXE)
	cp $(EXE) $(BIN_DIR)

clean:
	rm -f *.[ox] *.mod

$(EXE): ndviBareness4Wrf.o CalcNdviBareness_mod.o \
                  GIMMS_MODIS_NDVI_mod.o FieldWPS_mod.o \
                  FileWPS_mod.o NDVI_mod.o \
                  SPORT_MODIS_NDVI_mod.o 
	$(FC) $(FFLAGS) $^ -o $(EXE) $(LIBDIRS_NDVI) \
        -L$(UTILS_DIR)/shared -lshared

# Source code dependencies

CalcNdviBareness_mod.o: CalcNdviBareness_mod.f90
	$(FC) $(FFLAGS_NDVI) -c $<

FieldWPS_mod.o: FieldWPS_mod.f90
	$(FC) $(FFLAGS_NDVI) -c $<

FileWPS_mod.o: FileWPS_mod.f90 FieldWPS_mod.o
	$(FC) $(FFLAGS_NDVI) -c $<

GIMMS_MODIS_NDVI_mod.o: GIMMS_MODIS_NDVI_mod.f90
	$(FC) $(FFLAGS_NDVI) -c $< $(INCDIRS_NDVI)

NDVI_mod.o: NDVI_mod.f90 SPORT_MODIS_NDVI_mod.o GIMMS_MODIS_NDVI_mod.o 
	$(FC) $(FFLAGS_NDVI) -c $<

ndviBareness4Wrf.o: ndviBareness4Wrf.f90 CalcNdviBareness_mod.o \
                    FieldWPS_mod.o FileWPS_mod.o \
                    NDVI_mod.o
	$(FC) $(FFLAGS_NDVI) -c $<

SPORT_MODIS_NDVI_mod.o: SPORT_MODIS_NDVI_mod.f90
	$(FC) $(FFLAGS_NDVI) -c $<

