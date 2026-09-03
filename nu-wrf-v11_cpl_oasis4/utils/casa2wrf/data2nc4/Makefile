include $(CONFIG_DIR)/common.mk
.PHONY: clean all

EXE = ConvertData2Netcdf.x
all:  $(EXE)
	mkdir -p $(UTILS_DIR)/bin
	cp $(EXE) $(UTILS_DIR)/bin

%.x: %.o 
	$(FC) $^ -o $@ $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(INCS)

clean:
	rm -f *.[xo]

