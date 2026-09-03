#!/bin/bash --login

#SBATCH -A glm200001 
#SBATCH --job-name=fvcom_test
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --time=07:55:00

cd /compass/glm200001/bart753/coupled_run/FVCOM/run

module purge
module load pnnl_proxies
module load autoconf/2.72
module load oneapi-compiler/2024.2
module load oneapi-mkl/2024.2
module load openmpi/5.0.9
module load hdf5/1.14.6
module load netcdf-c/4.9.2
module load netcdf-fortran/4.6.1


export HDF5=$(module show hdf5 2>&1 | grep "CMAKE_PREFIX_PATH" | cut -d'"' -f4 | sed 's/\/include$//')
export NETCDF=/compass/glm200001/bart753/netcdf

export WRF_DIR=/compass/glm200001/bart753/WRF
export OASIS_DIR=/compass/glm200001/bart753/oasis3-mct
export OASIS_ENV=compass_intel


./fvcom --casename=gl
#mpirun -np 4 ./fvcom --casename=gl
