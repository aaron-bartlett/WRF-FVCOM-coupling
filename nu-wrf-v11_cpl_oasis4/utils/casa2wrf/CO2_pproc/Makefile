include $(CONFIG_DIR)/common.mk
.PHONY: clean all


EXE = Read_CO2_conc.x Read_CO2_Flux.x
all:  $(EXE)
	mkdir -p $(UTILS_DIR)/bin
	cp $(EXE) $(UTILS_DIR)/bin

%.x: %.o 
	$(FC) $^ -o $@ $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(INCS)

clean:
	rm -f *.[xo]

