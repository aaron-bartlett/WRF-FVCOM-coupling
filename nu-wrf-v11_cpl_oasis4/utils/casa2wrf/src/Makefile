#------------------------------------------------------------------------------
# NASA/GSFC, Software Integration & Visualization Office, Code 610.3
#------------------------------------------------------------------------------
#
# Makefile for CASA2WRF
#
# AUTHOR:
# Eric Kemp, NASA GSFC/Northrop Grumman
#
# 12/16 : Major mods (C. Cruz)
#------------------------------------------------------------------------------
include $(CONFIG_DIR)/common.mk
.PHONY: clean all

EXE = casa2wrf.x

all:  $(EXE)
	mkdir -p $(UTILS_DIR)/bin
	echo 'Installing...'$(EXE)
	cp $(EXE) $(UTILS_DIR)/bin

clean:
	rm -f *.[ox] *.mod

$(EXE): casa2wrf.o \
           netcdf_util_mod.o interp_util_mod.o sorted_datetimes_mod.o \
	   interpolator_mod.o \
           flux_interpolator_mod.o \
           derived_variables_mod.o \
           casa_variables_mod.o casa_emission_variables_mod.o \
           casa_flux_filenames_mod.o  
	$(FC) $(FFLAGS) $^ -o $(EXE) $(LIBS)

calc_pressure: calc_pressure.o
	$(FC) $(FFLAGS)  $^ -o calc_pressure $(LIBS)

# Source code dependencies
casa2wrf.o: casa2wrf.f90 interpolator_mod.o \
              sorted_datetimes_mod.o flux_interpolator_mod.o \
              casa_variables_mod.o \
               casa_emission_variables_mod.o casa_flux_filenames_mod.o 
	$(FC) $(FFLAGS) -c $<

interpolator_mod.o: interpolator_mod.f90 sorted_datetimes_mod.o \
	            casa_emission_variables_mod.o \
                    netcdf_util_mod.o casa_variables_mod.o \
                    derived_variables_mod.o \
                    interp_util_mod.o
	$(FC) $(FFLAGS) -c $< $(INCS)
flux_interpolator_mod.o: flux_interpolator_mod.f90 \
                    interpolator_mod.o \
                    netcdf_util_mod.o \
                    derived_variables_mod.o sorted_datetimes_mod.o \
                    interp_util_mod.o casa_flux_filenames_mod.o 
	$(FC) $(FFLAGS) -c $< $(INCS)

interp_util_mod.o: interp_util_mod.f90 sorted_datetimes_mod.o

# Generic
%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(INCS)

