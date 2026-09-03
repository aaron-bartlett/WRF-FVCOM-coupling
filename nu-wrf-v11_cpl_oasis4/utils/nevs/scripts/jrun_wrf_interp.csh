#!/bin/csh

# This script will start several wrf_interp.sp3 jobs at once on Discover[05-08, 15-18]
# login nodes.  This script assumes:
#  - you have a customized copy of namelist.vinterp inside each of the NU-WRF output 
#    directories, such as <NEVS_root>/round0_default/namelist.vinterp
#  - you have created the output subdirectory to contain the wrf_interp output files,
#    such as  <NEVS_root>/round0_default/wrfout_INTRP/ (from namelist.vinterp)
#
# The SBATCH version of this script is not currently feasible, because compute nodes do
# not have enough job cache (only 1GB) to hold all the data required, and the job slows
# to a crawl after a couple minutes.  --Try again when SCU-14 is online
#

#For bash
#ulimit -s unlimited
#For csh:
unlimit
source /usr/share/modules/init/csh
module purge
module load comp/intel-13.1.3.192
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/lib64

set SERVERS = (discover06 discover07 discover08 discover15 discover16 discover17 discover18 discover25 discover26 discover27 discover28 discover05 discover06 discover07 discover08 )

set rootPath = $argv[1]
shift argv

foreach wrf_dir ($argv)
   echo "Server = $SERVERS[1]"
   echo "wrf_dir = $wrf_dir"
   ssh -f $SERVERS[1] "cd $rootPath/$wrf_dir && nohup /discover/nobackup/projects/nu-wrf/lib/SLES11.3/wrf_interp/wrf_interp.sp3 >& wrfinterp_log"
   shift SERVERS
end

exit 0
