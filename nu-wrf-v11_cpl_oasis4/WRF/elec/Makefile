#

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f

MODULES =                                \
             ../elec/module_commasmpi.o \
             ../elec/module_mp_boxmgsetup.o \
             ../elec/module_mp_screen.o \
             ../elec/module_mp_discharge.o
            
LIBTARGET    =  elecphys
TARGETDIR    =  ./
$(LIBTARGET) :  MODULE DRIVERS
		$(AR) $(ARFLAGS) ../main/$(LIBWRFLIB) $(MODULES) $(OBJS)

MODULE  : $(MODULES)

DRIVERS : $(OBJS)

include ../configure.wrf

clean:
	@ echo 'use the clean script'

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

include depend.elec

            