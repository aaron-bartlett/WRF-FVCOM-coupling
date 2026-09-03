#------------------------------------------------------------------------------
# NASA/GSFC, Software Systems Support Office, Code 610.3
#------------------------------------------------------------------------------
#
# Makefile for plot_chem
#
# AUTHOR:
# Eric Kemp, NASA SSSO/Northrop Grumman
#
# 12/16 : Minor mods (C. Cruz)
#------------------------------------------------------------------------------
include $(CONFIG_DIR)/common.mk 

DIRS := src

.PHONY: all clean $(DIRS)
all: $(DIRS) 

$(DIRS):
	$(MAKE) --directory=$@

clean:
	@for d in $(DIRS) ; do \
	$(MAKE) -C $$d $@ ; done

