export LD_LIBRARY_PATH=/pexue3/chenfuh/local_libs/netcdf4.2.1.1/lib/:$LD_LIBRARY_PATH
#mpirun -n 1 ./fvcom --casename=gl
mpirun -n 4 ./fvcom --casename=gl
