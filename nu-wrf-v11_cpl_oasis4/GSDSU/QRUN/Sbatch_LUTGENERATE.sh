#!/bin/csh
#SBATCH -J GSDSU@LUT
#SBATCH -A s0942
#SBATCH -N 45 -n 1260 --ntasks-per-node=28:model=bro
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL
#SBATCH -o gsdsu_lut_run.txt
#SBATCH -e gsdsu_lut_err.txt

 rm gsdsu_error.txt
 rm gsdsu_run.txt
 limit stacksize unlimited

 module purge
 module load comp-intel/2016.2.181
 module load mpi-sgi/mpt.2.15r20 # up-to-date MPT module

 set EXE = 'GSDSU.x'
 set RUN_OPT_BEGIN = 'mpirun -np 1260'

 echo "${RUN_OPT_BEGIN} ./${EXE}"
 ${RUN_OPT_BEGIN} ./${EXE}

