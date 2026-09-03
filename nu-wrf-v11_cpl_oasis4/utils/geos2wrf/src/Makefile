#------------------------------------------------------------------------------
# NASA/GSFC, Software Systems Support Office, Code 610.3
#------------------------------------------------------------------------------
#
# Makefile fragment containing rule dependencies. Should be incorporated into
# compiler-specific Makefile via 'include Makefile.dependencies'.
#
# AUTHOR:  Eric Kemp, NASA CISTO/SSAI
#
# 12/16 : Major mods (C. Cruz)
# #------------------------------------------------------------------------------
include $(CONFIG_DIR)/common.mk 
BIN_DIR=$(UTILS_DIR)/bin
SHARED_LIB = $(UTILS_DIR)/shared/libshared.a
SCRIPTS_DIR = $(NUWRFDIR)/scripts/python/utils

.PHONY: clean all

EXE = createHGT.x createLANDSEA.x createPRESSURE.x createRH.x createSOILHGT.x \
       extrapIsobaric.x geos2wps.x merra2wrf.x splitWPS.x temporalInterpolation.x

all: $(EXE)
	mkdir -p $(BIN_DIR)
	echo 'Installing...'$(EXE)
	cp $(EXE) $(BIN_DIR)
	if [ -e $(BIN_DIR)/geos2wps.x ]; then \
	ln -sf $(BIN_DIR)/geos2wps.x $(BIN_DIR)/geos2wrf.x; fi
	cp $(SCRIPTS_DIR)/merra2wrf.py $(BIN_DIR)

clean:
	rm -f *.[ox] *.mod

createHGT.x: createHGT.o DerivedVars_mod.o FieldWPS_mod.o \
           FileWPS_mod.o WrfRH_mod.o
	$(FC) $^ -o createHGT.x $(SHARED_LIB)

createLANDSEA.x: createLANDSEA.o FieldWPS_mod.o  FileWPS_mod.o
	$(FC) $^ -o createLANDSEA.x $(SHARED_LIB)

createPRESSURE.x: createPRESSURE.o DerivedVars_mod.o FieldWPS_mod.o \
                 FileWPS_mod.o WrfRH_mod.o
	$(FC) $^ -o createPRESSURE.x $(SHARED_LIB)

createRH.x: createRH.o FieldWPS_mod.o FieldWPSList_mod.o \
          FileWPS_mod.o WrfRH_mod.o
	$(FC) $^ -o createRH.x $(SHARED_LIB)

createSOILHGT.x: createSOILHGT.o  FieldWPS_mod.o FileWPS_mod.o
	$(FC) $^ -o createSOILHGT.x $(SHARED_LIB)

extrapIsobaric.x: extrapIsobaric.o DerivedVars_mod.o FieldWPS_mod.o \
                FieldWPSList_mod.o  FileWPS_mod.o WrfRH_mod.o
	$(FC) $^ -o extrapIsobaric.x $(SHARED_LIB)

geos2wps.x: geos2wps.o NamelistGEOS2WPS_mod.o FileUtils_mod.o \
          Hdf4Utils_mod.o NetcdfUtils_mod.o FieldWPS_mod.o \
          FileWPS_mod.o GEOS2WPS_mod.o  CheckSubsetBounds_mod.o
	$(FC) $^ -o geos2wps.x $(LIBS)

merra2wrf.x: merra2wrf.o MerraData_mod.o FieldWPS_mod.o FileWPS_mod.o \
           NamelistMerra_mod.o FileUtils_mod.o Hdf4Utils_mod.o \
           NetcdfUtils_mod.o  WrfRH_mod.o CheckSubsetBounds_mod.o
	$(FC) $^ -o merra2wrf.x $(LIBS)

temporalInterpolation.x: temporalInterpolation.o DerivedVars_mod.o \
                       FieldWPS_mod.o \
                       FieldWPSList_mod.o  \
                       FileWPS_mod.o calc_epoch_sec.o WrfRH_mod.o
	$(FC) $^ -o temporalInterpolation.x $(SHARED_LIB)

splitWPS.x: splitWPS.o FieldWPS_mod.o FileWPS_mod.o 
	$(FC) $^ -o splitWPS.x  $(SHARED_LIB)

# Source code dependencies

calc_epoch_sec.o: calc_epoch_sec.c
	$(CC) $(CFLAGS) -c $<

CheckSubsetBounds_mod.o: CheckSubsetBounds_mod.f90
	$(FC) $(FFLAGS) -c $<

createHGT.o: createHGT.f90 DerivedVars_mod.o FieldWPS_mod.o \
             FileWPS_mod.o
	$(FC) $(FFLAGS) -c $<

createLANDSEA.o: createLANDSEA.f90 FieldWPS_mod.o  FileWPS_mod.o
	$(FC) $(FFLAGS) -c $<

createPRESSURE.o: createPRESSURE.f90 FieldWPS_mod.o \
                  FileWPS_mod.o DerivedVars_mod.o
	$(FC) $(FFLAGS) -c $<

createRH.o: createRH.f90 FieldWPS_mod.o FieldWPSList_mod.o \
            FileWPS_mod.o WrfRH_mod.o
	$(FC) $(FFLAGS) -c $<

createSOILHGT.o: createSOILHGT.f90  FieldWPS_mod.o FileWPS_mod.o
	$(FC) $(FFLAGS) -c $<

DerivedVars_mod.o: DerivedVars_mod.f90 FieldWPS_mod.o WrfRH_mod.o
	$(FC) $(FFLAGS) -c $<

extrapIsobaric.o: extrapIsobaric.f90 DerivedVars_mod.o FieldWPS_mod.o \
                  FieldWPSList_mod.o  FileWPS_mod.o WrfRH_mod.o
	$(FC) $(FFLAGS) -c $<

FieldWPSList_mod.o: FieldWPSList_mod.f90 FieldWPS_mod.o
	$(FC) $(FFLAGS) -c $<

FileUtils_mod.o: FileUtils_mod.f90 Hdf4Utils_mod.o NetcdfUtils_mod.o
	$(FC) $(FFLAGS) -c $< 

FileWPS_mod.o: FileWPS_mod.f90 FieldWPS_mod.o
	$(FC) $(FFLAGS) -c $<

geos2wps.o: geos2wps.f90 GEOS2WPS_mod.o NamelistGEOS2WPS_mod.o 
	$(FC) $(FFLAGS) -c $<

GEOS2WPS_mod.o: GEOS2WPS_mod.f90 FieldWPS_mod.o \
                FileUtils_mod.o FileWPS_mod.o NamelistGEOS2WPS_mod.o
	$(FC) $(FFLAGS) -c $<

Hdf4Utils_mod.o: Hdf4Utils_mod.f90 CheckSubsetBounds_mod.o
	$(FC) $(FFLAGS) -I$(INC_HDF4) -c $<

MerraData_mod.o: MerraData_mod.f90 FileWPS_mod.o FieldWPS_mod.o \
                 FileUtils_mod.o  WrfRH_mod.o
	$(FC) $(FFLAGS) -c $< 

merra2wrf.o: merra2wrf.f90 MerraData_mod.o FieldWPS_mod.o FileWPS_mod.o \
             NamelistMerra_mod.o 
	$(FC) $(FFLAGS) -c $<

NamelistGEOS2WPS_mod.o: NamelistGEOS2WPS_mod.f90 FileUtils_mod.o 
	$(FC) $(FFLAGS) -c $<

NamelistMerra_mod.o: NamelistMerra_mod.f90 FileUtils_mod.o 
	$(FC) $(FFLAGS) -c $<

NetcdfUtils_mod.o: NetcdfUtils_mod.f90 CheckSubsetBounds_mod.o
	$(FC) $(FFLAGS) -I$(INC_NETCDF) -c $<

splitWPS.o: splitWPS.f90 FieldWPS_mod.o FileWPS_mod.o 
	$(FC) $(FFLAGS) -c $<

temporalInterpolation.o: temporalInterpolation.f90 DerivedVars_mod.o \
                         FieldWPS_mod.o \
                         FieldWPSList_mod.o FileWPS_mod.o \
                         calc_epoch_sec.o WrfRH_mod.o
	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<
