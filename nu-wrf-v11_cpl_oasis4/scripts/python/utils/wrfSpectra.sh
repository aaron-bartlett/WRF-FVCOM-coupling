#!/bin/sh
#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Science and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# PROGRAM: wrfSpectra.sh
#
# AUTHOR:
# Eric Kemp, NASA CISTO/SSAI
#
# DESCRIPTION:
# Wrapper script for ke_spectra_arw.ncl, which calculates and plots mean 
# horizontal kinetic energy spectra from WRF-ARW 3D output. Heavily influenced
# by merra2Spectra.csh and ke_spectra_dd.ncl scripts developed by NASA GMAO
# for global lat/lon data.
#
# Requires installed version of NCL. If processing NU-WRF output, NCL should
# be compiled with NetCDF4 (HDF5 compression) support.
#
# REVISION:
# 22 Jul 2015 - First version.
#
#------------------------------------------------------------------------------

# Module environment for NASA Discover supercomputer. This sets environment
# variables like PATH for running NCL.
source /usr/share/modules/init/bash
module purge
unset LD_LIBRARY_PATH
module load other/ncl-6.3.0-static
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64

# Common settings for all WRF experiments to plot
export WRFFILE="wrfout_d01_2000-02-26_00:00:00" # WRF file name

export DIRO="./" # Output directory
export PFMT="PNG" # Output format for plots
export KE_MLEVEL=21 # WRF model level to process

# Experiment specific settings
TOP_DIR_WORK="/discover/nobackup/projects/dscale/tiguchi3/runs_hasw"
TOP_DIR_PUB="/discover/nobackup/projects/dscale/tiguchi3/pub/lbc_r3135_outputs"
wrf_dirs[0]="$TOP_DIR_WORK/merra2_10years_workdirect_lbcr3135_unnudge_B24"
wrf_dirs[1]="$TOP_DIR_WORK/merra2_10years_workdirect_lbcr3135_nudge2000_B24"
wrf_dirs[2]="$TOP_DIR_WORK/merra2_10years_workdirect_lbcr3135_nudge600_B24"
wrf_dirs[3]="$TOP_DIR_WORK/merra2_10years_workdirect_lbcr3135_unnudge_B12"
wrf_dirs[4]="$TOP_DIR_PUB/B12_nudge2000km/1999-2001"
wrf_dirs[5]="$TOP_DIR_PUB/B12_nudge600km/1999-2001"
wrf_dirs[6]="$TOP_DIR_WORK/merra2_10years_workdirect_lbcr3135_unnudge_B4"
wrf_dirs[7]="$TOP_DIR_WORK/merra2_10years_workdirect_lbcr3135_nudge2000_B4"
wrf_dirs[8]="$TOP_DIR_WORK/merra2_10years_workdirect_lbcr3135_nudge600_B4"

wrf_cases[0]="B24_NONUDGE"
wrf_cases[1]="B24_NUDGE2000KM"
wrf_cases[2]="B24_NUDGE600KM"
wrf_cases[3]="B12_NONUDGE"
wrf_cases[4]="B12_NUDGE2000KM"
wrf_cases[5]="B12_NUDGE600KM"
wrf_cases[6]="B4_NONUDGE"
wrf_cases[7]="B4_NUDGE2000KM"
wrf_cases[8]="B4_NUDGE600KM"

for i in 0 1 2 3 4 5 6 7 8 ; do

    export CASE="${wrf_cases[$i]}"
    export PNAME=$CASE
    export DIRI="${wrf_dirs[$i]}"

    # Skip if WRF file doesn't exist
    file=$DIRI/$WRFFILE
    if [ ! -f $file ] ; then
	echo "WARNING, missing file $file"
	echo "Skipping to next experiment..."
	continue
    fi

    # Call NCL script to calculate and plot spectra.
    ncl 'fnamWRF="'${WRFFILE}'"' 'ptit="'${CASE}'"' ke_spectra_arw.ncl || exit 1
    sleep 1

done


# The end.
exit 0