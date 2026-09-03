#!/bin/sh
#SBATCH -J Preproc
#SBATCH -N 1 -n 1 --ntasks-per-node=1
#SBATCH -t 1:00:00
#SBATCH -A s0942
#SBATCH -o casa2w.pproc.out
#SBATCH -p general
#------------------------------------------------------------------------------
# NASA/GSFC, Software Systems Support Office, Code 610.3           
#------------------------------------------------------------------------------
#                                                                              
# SCRIPT:  run_casa2wrfPPROC.discover.sh
#                                                                              
# AUTHOR:                                                                      
# Jossy Jacob
# Jan 2, 2014                                                                            
# DESCRIPTION:                                                                 
# Sample batch script for running casa2wrf data preprocessor on NASA GSFC Discover 
# supercomputer with SLURM.
#
#------------------------------------------------------------------------------

# When a batch script is started, it starts in the user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi
set srcdir = /discover/nobackup/NUWRF/v6-3.4.1_CO2/utils/casa2wrf/pproc/
set workdir = $SLURM_SUBMIT_DIR 
# Load config file for modules and paths.
source ./config.discover.sh || exit 1
# Compile the code if no exectuable
if [ ! -e Read_CO2_conc.x ] then 
  cd $srcdir
  make -f Makefile_CO2_conc
  make -f Makefile_CO2_Flux
fi
cd $workdir
ln -sf $srcdir/Read_CO2_conc.x .
ln -sf $srcdir/Read_CO2_Flux.x .
#
# Run preprocessing of CO2 concentration files
# These conc files are ndays=5 days data files, with 1 hrly data.
# Output is 6 hourly
set indir = 'conc/'
set file_list1 = `ls -c1 ${indir}binary/` 
foreach file ( $file_list1 )
 Read_CO2_conc.x $file $indir 1 6 5
end
#
# Run preprocessing of CO2 Flux files
#
#Flux files are yearly data files. If there is changes, we need to change dates in the code. 
#Output is 3 hourly
#
set indir = 'flux/'
set file_list2 = `ls -c1 ${indir}binary/` 
foreach file ( $file_list2 )
 Read_CO2_Flux.x $file $indir
end
#
# The end
exit 0
