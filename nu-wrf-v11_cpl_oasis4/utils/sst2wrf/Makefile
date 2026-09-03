#------------------------------------------------------------------------------
# NASA/GSFC, Software Integration & Visualization Office, Code 610.3
#------------------------------------------------------------------------------
#
# Top level Makefile for sst2wrf program
#
# AUTHOR:  Eric Kemp, NGIS
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

