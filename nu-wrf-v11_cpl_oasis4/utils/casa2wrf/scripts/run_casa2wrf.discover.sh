#!/bin/sh
#SBATCH -J casa2w
#SBATCH -N 1 -n 1 --ntasks-per-node=1
#SBATCH -t 1:00:00
#SBATCH -A s0942
#SBATCH -o casa2w.slurm.out
#SBATCH -p general
#------------------------------------------------------------------------------
# NASA/GSFC, Software Systems Support Office, Code 610.3           
#------------------------------------------------------------------------------
#                                                                              
# SCRIPT:  run_casa2wrf.discover.sh
#                                                                              
# AUTHOR:                                                                      
# Eric Kemp, NASA SSSO/SSAI
# Modified for casa2wrf by Jossy Jacob 
#                                                                             
# DESCRIPTION:                                                                 
# Sample batch script for running casa2wrf on NASA GSFC Discover 
# supercomputer with SLURM.
#
#------------------------------------------------------------------------------

# When a batch script is started, it starts in the user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Load config file for modules and paths.
source ./config.discover.sh || exit 1

# Location of CO2/wrfout data assumed to be in conc, flux, WRFout
# Create links to CO2/wrfout data
if [ ! -e conc ] ; then
   ln -s /path/to/CASA_CO2_conc conc
fi
if [ ! -e conc ] ; then
   ln -s /path/to/CASA_CO2_flux flux
fi
if [ ! -e WRFout ] ; then
   ln -sf /path/to/WRFout WRFout
fi

# Go to work directory and make sure namelist.casa2wrf is present.
if [ -z "$WORKDIR" ] ; then
    echo "ERROR, WORKDIR is not defined!"
    exit 1
fi
cd $WORKDIR || exit 1
if [ ! -e namelist.casa2wrf ] ; then
    echo "ERROR, namelist.casa2wrf not found!"
    exit 1
fi

# Check/create output directory
if [ ! -d "$WORKDIR/chem_flux" ] ; then
   mkdir $WORKDIR/chem_flux
fi

# Make backup copies of wrfbdy and wrfinput files, as casa2wrf will add
# new variables to them.
if [ ! -e wrfbdy_d01 ] ; then
    echo "ERROR, wrfbdy_d01 not found!"
    cp original_data/wrfbdy_d01 wrfbdy_d01
fi

if [ ! -e wrfinput_d01 ] ; then
  cd original_data
  files1=`ls wrfinput_d0[1-9]`
  cd ..
  for file in $files1 ; do
    	cp original_data/$file $file
        echo 'Copied wrfinput data'
  done
fi
# Run casa2wrf.  No MPI is used since the program is serial.
ln -fs $NUWRFDIR/utils/casa2wrf/bin/casa2wrf casa2wrf || exit 1
if [ ! -e casa2wrf ] ; then
    echo "ERROR, casa2wrf not found!"
    exit 1
fi
./casa2wrf | tee casa2wrf.out

# The end
exit 0
