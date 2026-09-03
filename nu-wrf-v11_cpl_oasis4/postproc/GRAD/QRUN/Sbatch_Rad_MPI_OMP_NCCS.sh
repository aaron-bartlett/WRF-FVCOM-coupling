#!/bin/csh
#SBATCH -J goddardrad
#SBATCH -A s0942
#SBATCH -N 1 -n 1 --ntasks-per-node=1 -C west
#SBATCH -t 0:03:00
#SBATCH --mail-type=ALL
#SBATCH -o sbatch_runtime.txt
#SBATCH -e sbatch_error.txt

 set EXE = 'rad.x'
 set RUN_OPT_BEGIN = 'mpirun -np 1'

 setenv OMP_STACKSIZE 512m
 setenv LD_LIBRARY_PATH  /usr/local/intel/Composer/composer_xe_2013_sp1.2.144/compiler/lib/intel64/
 setenv OMP_NUM_THREADS 12

 echo "${RUN_OPT_BEGIN} ./${EXE}"
 ${RUN_OPT_BEGIN} ./${EXE}

