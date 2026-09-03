#------------------------------------------------------------------------------
# NASA/GSFC, Software Integration & Visualization Office, Code 610.3
#------------------------------------------------------------------------------
#
# Makefile for SST binary file read routines.  Source code provided by
# Remote Sensing Systems.
#
# AUTHOR:  Eric Kemp, NGIS
#
#------------------------------------------------------------------------------

all: read_rss_mwir_sst.o read_rss_oisst_v4.o
	$(AR) $(ARFLAGS) $(LIBSST) read_rss_mwir_sst.o read_rss_oisst_v4.o

#------------------------------------------------------------------------------

read_rss_mwir_sst.o: read_rss_mwir_sst.f90
	$(FC) $(FFLAGS) -c read_rss_mwir_sst.f90

read_rss_oisst_v4.o: read_rss_oisst_v4.f90
	$(FC) $(FFLAGS) -c read_rss_oisst_v4.f90

#------------------------------------------------------------------------------

.PHONY:  clean
clean:
	-rm -f *.o *.a 
