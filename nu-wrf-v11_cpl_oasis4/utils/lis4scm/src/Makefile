#------------------------------------------------------------------------------
# NASA/GSFC, Software Integration & Visualization Office, Code 610.3
#------------------------------------------------------------------------------
#
# Makefile for lis_wrf_domain
#
# AUTHOR:
# Eric Kemp, NASA SIVO/Northrop Grumman
#
# 12/16 : Major mods (C. Cruz)
#------------------------------------------------------------------------------
include $(CONFIG_DIR)/common.mk
BIN_DIR=$(UTILS_DIR)/bin

.PHONY: clean all

EXE = lis4scm.x

all: $(EXE)
	mkdir -p $(BIN_DIR)
	echo 'Installing...'$(EXE)
	cp $(EXE) $(BIN_DIR)

clean:
	rm -f *.[ox] *.mod

$(EXE): lis4scm.o
	$(FC) $^ -o $@ $(LIBS)

lis4scm.o: lis4scm.f90
	$(FC) $(FFLAGS) -c $<

