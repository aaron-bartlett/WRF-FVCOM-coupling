#
# Top Makefile for NU-WRF utils.
#

UTILS_DIR := $(CURDIR)
CONFIG_DIR := $(CURDIR)/config

export UTILS_DIR
export CONFIG_DIR

include $(CONFIG_DIR)/common.mk 

SHARED := shared
ifeq ($(COMPILER_VENDOR),intel)
DIRS := casa2wrf sst2wrf geos2wrf gocart2wrf lis4scm \
	lisWrfDomain ndviBareness4Wrf plot_chem prep_chem_sources
else
DIRS := casa2wrf sst2wrf geos2wrf gocart2wrf lis4scm \
	lisWrfDomain plot_chem prep_chem_sources
endif

.PHONY: all clean allclean $(DIRS) $(SHARED)

all: $(DIRS) 

$(DIRS) $(SHARED):
	$(MAKE) --directory=$@


ifdef DIR
clean:
	$(MAKE) -C $(DIR) clean
endif

allclean: 
	@for d in $(DIRS) $(SHARED) ; do \
	$(MAKE) -C $$d clean ; done
	rm -rf $(CURDIR)/bin

# Dependencies
$(DIRS): $(SHARED)
