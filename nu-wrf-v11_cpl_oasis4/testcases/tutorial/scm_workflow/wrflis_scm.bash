#!/bin/bash
#SBATCH -J scm
#SBATCH -o /discover/nobackup/ccruz/scratch/nu-wrf/dev_scm/results/intel-intelmpi.2019-12-13_13-31-47/wrflis_scm/scm.out
#SBATCH -e /discover/nobackup/ccruz/scratch/nu-wrf/dev_scm/results/intel-intelmpi.2019-12-13_13-31-47/wrflis_scm/scm.err
#SBATCH --account=s0942
#SBATCH --ntasks=1
#SBATCH --constraint=hasw
#SBATCH --qos=debug
#SBATCH --time=1:00:00
cd /discover/nobackup/ccruz/scratch/nu-wrf/dev_scm/results/intel-intelmpi.2019-12-13_13-31-47/wrflis_scm
# source common file for modules and paths 
source ./common.reg || exit 1 
if [ -d /discover/nobackup/projects/nu-wrf/regression_testing/data/Charney/wrflis/scm ]; then 
ln -sf /discover/nobackup/projects/nu-wrf/regression_testing/data/Charney/wrflis/scm
fi 
python /discover/nobackup/ccruz/devel/nu-wrf/code/gitlab/nu-wrf/scripts/python/regression/regression.py wrflis_scm
 
