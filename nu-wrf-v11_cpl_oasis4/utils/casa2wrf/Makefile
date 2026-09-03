include $(CONFIG_DIR)/common.mk 

DIRS := src CO2_pproc data2nc4

.PHONY: all clean $(DIRS)
all: $(DIRS) 

$(DIRS):
	$(MAKE) --directory=$@

clean:
	@for d in $(DIRS) ; do \
	$(MAKE) -C $$d $@ ; done\

