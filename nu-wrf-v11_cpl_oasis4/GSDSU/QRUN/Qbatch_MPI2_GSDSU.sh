#PBS -S /bin/csh
#PBS -N gsdsu@p3
###PBS -l select=1:ncpus=37:mpiprocs=37:model=sky_ele  
#PBS -l select=2:ncpus=18:mpiprocs=18:model=has
#PBS -l walltime=00:20:00
#PBS -W group_list=s1183
#PBS -m abe
#PBS -q devel
#PBS -V
#PBS -e ./gsdsu_error.txt
#PBS -o ./gsdsu_run.txt

 rm gsdsu_error.txt
 rm gsdsu_run.txt
 limit stacksize unlimited

 module purge
 module load comp-intel/2016.2.181
 module load mpi-sgi/mpt.2.15r20 # up-to-date MPT module

 set EXE = 'GSDSU.x'
 echo $EXE
 set RUN_OPT_BEGIN = 'mpiexec -np 36'

 echo "${RUN_OPT_BEGIN} ./${EXE}"
 ${RUN_OPT_BEGIN} ./${EXE} 

