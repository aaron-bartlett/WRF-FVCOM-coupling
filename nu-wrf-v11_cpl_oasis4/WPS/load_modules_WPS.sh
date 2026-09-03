#!/bin/bash --login

module purge
module load pnnl_proxies
module load autoconf/2.72
module load oneapi-compiler/2024.2
module load oneapi-mkl/2024.2
module load openmpi/5.0.9
module load hdf5/1.14.6

#module load netcdf-c/4.9.2
#module load netcdf-fortran/4.6.1

export HDF5=$(module show hdf5 2>&1 | grep "CMAKE_PREFIX_PATH" | cut -d'"' -f4 | sed 's/\/include$//')
export NETCDF=/compass/glm200001/cmu/coupled-run/libs/netcdf4
export ZLIB=/compass/glm200001/cmu/coupled-run/libs/zlib
export PATH="$NETCDF/bin:$HDF5/bin:$PATH"
export LD_LIBRARY_PATH="$NETCDF/lib:$HDF5/lib:$LD_LIBRARY_PATH"

export JASPERLIB=/compass/glm200001/cmu/coupled-run/libs/grib2/lib
export JASPERINC=/compass/glm200001/cmu/coupled-run/libs/grib2/include
export LD_LIBRARY_PATH=$JASPERLIB:$LD_LIBRARY_PATH

export NETCDF4=1
export NETCDF4_DEP=1
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export HDF5_USE_FILE_LOCKING=FALSE

export WRF_DIR=/compass/glm200001/cmu/coupled-run/nu-wrf-v11_cpl_oasis4/WRF
export OASIS_DIR=/compass/glm200001/cmu/coupled-run/libs/oasis3-mct
export OASIS_ENV=compass_intel

ulimit -s unlimited
