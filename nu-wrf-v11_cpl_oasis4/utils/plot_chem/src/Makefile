include $(CONFIG_DIR)/common.mk
BIN_DIR=$(UTILS_DIR)/bin

# Over-ride defaults
FC = ncargf90

.PHONY: all clean
EXE = plot_chem.x

all: $(EXE)
	mkdir -p $(BIN_DIR)
	echo 'Installing...'$(EXE)
	cp $(EXE) $(BIN_DIR)

clean:
	rm -f *.[ox] *.mod

$(EXE): plot_chem.o ncarg_util_mod.o grads_util_mod.o
	$(FC) $(FFLAGS) -o $(EXE) $^

# Source code dependencies

plot_chem.o: plot_chem.f90 ncarg_util_mod.o grads_util_mod.o
	$(FC) $(FFLAGS) -c $<

ncarg_util_mod.o: ncarg_util_mod.f90
	$(FC) $(FFLAGS) -c $<

grads_util_mod.o: grads_util_mod.f90
	$(FC) $(FFLAGS) -c $<


