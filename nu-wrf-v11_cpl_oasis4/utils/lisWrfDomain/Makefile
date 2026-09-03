#------------------------------------------------------------------------------
# NASA/GSFC, Software Integration & Visualization Office, Code 610.3
#------------------------------------------------------------------------------
#
# Top level Makefile for lis_wrf_domain program.
#
# AUTHOR:  Eric Kemp, NASA SIVO/Northrop Grumman
#
# 12/16 : Major mods (C. Cruz)
#------------------------------------------------------------------------------
include $(CONFIG_DIR)/common.mk 

DIRS := src

.PHONY: all clean $(DIRS)
all: $(DIRS) 

$(DIRS):
	$(MAKE) --directory=$@

clean:
	@for d in $(DIRS) ; do \
	$(MAKE) -C $$d $@ ; done\

