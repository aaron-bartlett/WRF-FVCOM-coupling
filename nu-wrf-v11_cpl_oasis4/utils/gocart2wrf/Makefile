#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Science and Technology Office, 
# Code 606.
#------------------------------------------------------------------------------
#
# Makefile dependencies file for GOCART2WRF. Should be incorporated in 
# compiler-specific Makefile via 'include Makefile.dependencies'
#
# AUTHOR:  Eric Kemp, NASA CISTO/SSAI
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

