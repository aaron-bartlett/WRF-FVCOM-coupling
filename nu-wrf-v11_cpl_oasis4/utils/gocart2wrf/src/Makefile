#------------------------------------------------------------------------------
# NASA/GSFC, Software Integration & Visualization Office, Code 610.3
#------------------------------------------------------------------------------
#
# Makefile for GOCART2WRF
#
# AUTHOR:
# Eric Kemp, NASA GSFC/Northrop Grumman
#
# 12/16 : Major mods (C. Cruz)
#------------------------------------------------------------------------------
include $(CONFIG_DIR)/common.mk
BIN_DIR=$(UTILS_DIR)/bin

.PHONY: clean all

EXE = gocart2wrf.x

all: $(EXE)
	mkdir -p $(BIN_DIR)
	echo 'Installing...'$(EXE)
	cp $(EXE) $(BIN_DIR)

clean:
	rm -f *.[ox] *.mod

$(EXE): gocart2wrf.o \
           netcdf_util_mod.o interp_util_mod.o interpolator_mod.o \
           derived_variables_mod.o sorted_datetimes_mod.o \
           gocart_variables_mod.o
	$(FC) $(FFLAGS) $^ -o $(EXE) $(LIBS)

calc_pressure.x: calc_pressure.o
	$(FC) $(FFLAGS)  $^ -o calc_pressure.x $(LIBS)


# Source code dependencies
gocart2wrf.o: gocart2wrf.f90 interpolator_mod.o \
              sorted_datetimes_mod.o gocart_variables_mod.o
	$(FC) $(FFLAGS) -c $<

interpolator_mod.o: interpolator_mod.f90 netcdf_util_mod.o \
                    derived_variables_mod.o sorted_datetimes_mod.o \
                    interp_util_mod.o
	$(FC) $(FFLAGS) -c $< $(INCS)

interp_util_mod.o: interp_util_mod.f90 sorted_datetimes_mod.o

# Generic
%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(INCS)

