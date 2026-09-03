SHELL=/bin/sh
EXEC =		GSDSU.x
#
# - dependency rule -
# 1) module_floating_point.o must comes first. 
# 2) module_simulator.o (and other 3D array looping module) depends on
#    module_mpi.o.  So they should comes after module_mpi.
# 3) main_SDSU.o must be at the end.  
#

#OBJSGCEM=	rdcmd.o scatdb.o module_mpi.o module_scatdb.o module_gocart_offline.o module_simulator.o module_GV.o Type_Kinds.o NESDIS_LandEM_Module.o module_telsem.o module_surface.o module_opt_micro.o simulator_micro.o module_radar.o module_lidar.o module_opt_visir.o simulator_visir.o simulator_broad.o module_opt_broad.o module_scan.o main_SDSU.o
OBJSGCEM=       module_floating_point.o module_const.o Type_Kinds.o rdcmd.o \
                scatdb.o module_scatdb.o module_tmatrix.o module_gmp_4ice.o module_psd.o \
                module_convolution.o \
                NESDIS_LandEM_Module.o module_telsem.o module_radar_power.o \
                module_mpi.o module_gocart_interact.o module_simulator.o \
                module_gas_profile.o module_sarta.o \
                module_GV.o module_surface.o module_opt_micro.o \
                simulator_micro.o module_radar.o module_lidar.o \
                module_opt_visir.o  simulator_visir.o simulator_broad.o \
                module_opt_broad.o module_scan.o main_SDSU.o

############################################################################
#                        Makefile User Modification
#                      Satellite Data Simulation Unit       
#
# CPP is C-preprocessor
# CF is fortran compiler
# CFFLAGS is fortran compiler flag
# INC_NETCDF is netcdf include file directly
#QRUN     = ./../QRUN_NCCS
#QRUN      = /path/to/qrun/directory
QRUN      = ./../QRUN

# Desk-top MAC  
CPP     = cpp -C -traditional-cpp
#CF      = gfortran # serial option
CF      = mpif90 
CFFLAGS = -O2 -funroll-loops -ftree-vectorize
CCFLAGS = -O2 -funroll-loops -ftree-vectorize
CC      = gcc
INC_NETCDF = /path/to/netcdf/include
LD_NETCDF = /path/to/netcdf/lib
LINK_MPI =   -lmpi # Make blank for non

# This is for location of define_CPP.h , and universal for any machine. 
INCLUDE_DIR = -I./

#
#
############################################################################
# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build executable.
#
all:		$(EXEC) find_thread_number
	cp GSDSU.x $(QRUN) 
	cp FIND_THREAD_NUMBER.x $(QRUN)/HOW_MANY_CPU_GSDSU

$(EXEC):	$(OBJSGCEM)
	$(CF) -o $(@) $(OBJSGCEM) -L$(LD_NETCDF) -lnetcdff -lnetcdf $(NETCDF4_DEP_LIB) $(LINK_MPI)

.SUFFIXES :	.o  .f90 .F

.F.f90 :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i > $*.f90 ; rm $*.i

.F.o :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i  >$*.f90 ; rm $*.i
	${CF} ${CFFLAGS} ${DBUG} -I$(INC_NETCDF) -c $*.f90

.f90.o :
	${CF} ${CFFLAGS} ${DBUG} -I$(INC_NETCDF) -c $*.f90

.c.o : 
	${CC} ${CCFLAGS} -c $*.c

.PHONY:  find_thread_number
find_thread_number: FIND_THREAD_NUMBER.o
	$(CF) -o FIND_THREAD_NUMBER.x $^

clean:
	-rm -f *.f90 *.o *.x *.mod PI* $(QRUN)/*.x $(QRUN)/HOW_MANY_CPU_GSDSU

