#!/bin/csh
#SBATCH -J GRAD@RUN
#SBATCH -A s0942
#SBATCH -N 1 -n 1 --ntasks-per-node=1 --constraint=hasw
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL
#SBATCH -o grad_run.txt
#SBATCH -e grad_err.txt

 rm grad_run.txt
 rm grad_err.txt
 limit stacksize unlimited

 module purge
 module load other/comp/gcc-5.3-sp3
 module load comp/intel-15.0.3.187
 module load lib/mkl-15.0.3.187
 module load mpi/sgi-mpt-2.12
 module load tool/allinea-tools-5.0.1

 set EXE = 'GRAD.x'
 set RUN_OPT_BEGIN = 'mpirun -np 1'

 rm sbatch_runtime.txt
 rm sbatch_error.txt

 echo "${RUN_OPT_BEGIN} ./${EXE}"
 ${RUN_OPT_BEGIN} ./${EXE}
