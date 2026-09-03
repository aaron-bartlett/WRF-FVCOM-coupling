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
SCRIPTS_DIR=$(NUWRFDIR)/scripts/python/utils

.PHONY: clean all

EXE = lisWrfDomain.x

all: $(EXE)
	mkdir -p $(BIN_DIR)
	echo 'Installing...'$(EXE)
	cp $(EXE) $(BIN_DIR)
	cp $(SCRIPTS_DIR)/adjust_liswrf_domain.py $(BIN_DIR)

clean:
	rm -f *.[ox] *.mod

$(EXE): lisWrfDomain.o DomainData_mod.o
	$(FC) $^ -o $@ $(LIBS)

lisWrfDomain.o: lisWrfDomain.f90 DomainData_mod.o 
	$(FC) $(FFLAGS) -c $<

DomainData_mod.o: DomainData_mod.f90
	$(FC) $(FFLAGS) -c $<

