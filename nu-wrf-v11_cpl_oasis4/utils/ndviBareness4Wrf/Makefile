#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Sciences and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# Top level Makefile for NDVI processing software.
#
# AUTHOR:  Eric Kemp, NASA GSFC/SSAI
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

