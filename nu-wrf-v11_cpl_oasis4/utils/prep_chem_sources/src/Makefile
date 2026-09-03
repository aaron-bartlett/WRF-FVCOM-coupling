include $(CONFIG_DIR)/common.mk
BIN_DIR=$(UTILS_DIR)/bin
OPT ?= $(MAKEPSC_OPT)
CHEM ?= $(CHEM_TYPE)
AER ?= $(AER_TYPE)
NETCDF_INC ?= -I$(NETCDF_4_INCLUDE)
HDF5_INC ?= -I$(HDF5_INCLUDE)
NETCDF_LIBS ?=  $(LIB_NETCDF)
HDF_LIBS ?= $(LIB_HDF4) $(LIB_HDF5) -L$(LIB_ZLIB) -lz

.PHONY: clean all
# Define path and compilation 
include paths.mk
include include.mk.$(MAKEPSC_OPT)

# Compiler commands.
INCLUDES   = -I$(LIB_RAMS_PATH)/include \
	     $(NETCDF_INC) $(HDF5_INC) \
	     -I$(UTILS_DIR)/shared 
LIBS       = -L$(UTILS_DIR)/shared -lshared \
	     $(NETCDF_LIBS) $(HDF_LIBS) 
F_COMMAND  = $(F_COMP) -c $(F_OPTS) $(INCLUDES)
C_COMMAND  = $(C_COMP) -c $(C_OPTS) -D$(CMACH) -D$(CUNDERSCORE) $(INCLUDES)

# Define objects.

include objects_prepsrc.mk

# Define targets.

EXE       = prep_chem_sources_$(CHEM_TYPE)_$(AER_TYPE).x
RAMSUTILS = ramsutils.a
RAMSLIB   = ramslib.a
WPSLIB    = wpslib.a  

all:  $(RAMSUTILS) $(RAMSLIB) $(WPSLIB) $(WRFFIMUTILS) $(EXE)
	mkdir -p $(BIN_DIR)
	echo 'Installing...'$(EXE)
	cp $(EXE) $(BIN_DIR)/prep_chem_sources.x

$(RAMSUTILS): $(OBJ_RAMS_UTILS)
	$(ARCHIVE) $(RAMSUTILS) $(OBJ_RAMS_UTILS)

$(RAMSLIB): $(OBJ_RAMS)
	$(ARCHIVE) $(RAMSLIB) $(OBJ_RAMS)

$(WPSLIB): $(OBJ_WPS)
	$(ARCHIVE) $(WPSLIB) $(OBJ_WPS)

$(EXE): $(OBJ_PREPSRC) $(RAMSLIB) $(RAMSUTILS) $(WPSLIB)
	$(LOADER) $(LOADER_OPTS) -o $(EXE) $(OBJ_PREPSRC) \
        $(RAMSLIB) $(RAMSUTILS) $(WPSLIB) $(LIBS)


clean:
	rm -f *.[oa] *.mod *.oo *.ipo \
        $(LIBUTILS) $(PREPSOURCE) $(RAMSLIB) $(WPSLIB) $(WRFFIMUTILS)  $(EXE) 

# Define dependencies and rules

include depend_prepsrc.mk
