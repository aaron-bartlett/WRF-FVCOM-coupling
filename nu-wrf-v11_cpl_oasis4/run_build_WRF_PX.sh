#export LIBROOT=/discover/nobackup/projects/nu-wrf/lib/intel-sgimpt-bjerknes-p5
#export LIBDIR_TAG=/discover/nobackup/projects/nu-wrf/lib/sles12/ekman/intel-intelmpi
export COMPILER_VENDOR="intel"
export MPI_VENDOR="intelmpi"
#./build.sh --config nu-wrf.cfg allclean
#./build.sh --config nu-wrf.cfg wrf

#./build.sh lis wrf
#./build.sh wps
#./build.sh rip
#./build.sh arwpost
#./build.sh utils
#./build.sh ldt

./build.sh --config nu-wrf.cfg wrf
