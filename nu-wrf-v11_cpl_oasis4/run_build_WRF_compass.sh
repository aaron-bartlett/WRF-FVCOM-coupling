source /compass/glm200001/cmu/coupled-run/load_modules.sh

#./build.sh --config nu-wrf.cfg allclean
#./build.sh --config nu-wrf.cfg wrf

#./build.sh lis wrf
#./build.sh wps
#./build.sh rip
#./build.sh arwpost
#./build.sh utils
#./build.sh ldt

./build.sh --config nu-wrf.cfg allclean
./build.sh --config nu-wrf.cfg wrfonly
nohup ./build.sh --config nu-wrf.cfg wrfonly > log.build 2>&1 &
