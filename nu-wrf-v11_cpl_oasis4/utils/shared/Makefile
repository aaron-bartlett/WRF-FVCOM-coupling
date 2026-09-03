include $(CONFIG_DIR)/common.mk
.PHONY : all clean

VPATH = $(NUWRF_DIR)/WRFV3/external/io_netcdf

LIB = libshared.a

OBJS = FileUnit_mod.o nml.o readh5dataset.o strtok.o
$(LIB): $(OBJS)
	-rm -f $(LIB)
	ar rcs $(LIB) $(OBJS)

all: $(LIB)

EXT_LIB = -L$(BASELIBS)/netcdf4/lib -lnetcdff -L$(BASELIBS)/hdf5/lib \
	  -L$(BASELIBS)/zlib/lib -L$(BASELIBS)/hdf4/lib -L$(BASELIBS)/jpeg/lib \
	  -L$(BASELIBS)/netcdf4/lib -lnetcdf -lnetcdf -lmfhdf -lhdf5_hl -lhdf5 -ldf -ljpeg -lz
	
%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

clean:
	rm -f *.o *.mod *.a *.x

nml.o: nml.f90 FileUnit_mod.o

