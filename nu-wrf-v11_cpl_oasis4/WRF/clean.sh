#!/usr/bin/env bash
#------------------------------------------------------------------------------
# NASA/GSFC, Software Integration & Visualization Office, Code 610.3           
#------------------------------------------------------------------------------
#                                                                              
# SCRIPT:  build.sh                                                            
#                                                                              
# AUTHOR:                                                                      
# Rob Burns and Eric Kemp, NGIS/NASA SIVO                                      
#                                                                              
# DESCRIPTION:                                                                 
# Wrapper script for cleaning WRFV3 component of NU-WRF.                       
#                                                                              
# REVISION HISTORY:                                                            
# 11 Apr 2011 - Updated comments.                                           
#------------------------------------------------------------------------------

echo "Cleaning WRFV3"

#------------------------------------------------------------------------------
# Grab the clean script directory and make sure all commands run from here.
#------------------------------------------------------------------------------

NU_WRF_PATH=$(cd $(dirname $0) ; pwd)
cd $NU_WRF_PATH

#------------------------------------------------------------------------------
# Source the nuwrf_util_lib.sh script to load additional functions.
#------------------------------------------------------------------------------

. $NUWRF_DIR/scripts/other/nuwrf_util_lib.sh

#------------------------------------------------------------------------------
# Make sure to blow away any previous build.
#------------------------------------------------------------------------------

./clean -a
nuwrf_check_exit "WRFV3 clean script"

#EMK...This can fail due to missing configure.lis file, so we'll create a fake
#one with 'touch' if necessary.  
cd $NU_WRF_PATH/lis/make
touch configure.lis
make clean
nuwrf_check_exit "WRFV3/lis/make make clean command"

cd $NU_WRF_PATH/lis/make/MAKDEP
make clean
nuwrf_check_exit "WRFV3/lis/make/MAKDEP make clean command"

#cd $NU_WRF_PATH/lis/lib/w3lib
#make clean
#nuwrf_check_exit "WRFV3/lis/lib/w3lib make clean command"

#cd $NU_WRF_PATH/lis/lib/grib
#make clean
#nuwrf_check_exit "WRFV3/lis/lib/grib make clean command"

#cd $NU_WRF_PATH/lis/lib/read_grib
#make clean
#nuwrf_check_exit "WRFV3/lis/lib/read_grib make clean command"

#------------------------------------------------------------------------------
# Normal exit:  Go back to the WRFV3 directory and exit.
#------------------------------------------------------------------------------

cd $NU_WRF_PATH
exit 0
