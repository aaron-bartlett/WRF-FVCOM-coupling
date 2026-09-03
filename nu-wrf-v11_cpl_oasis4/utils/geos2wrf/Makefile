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

