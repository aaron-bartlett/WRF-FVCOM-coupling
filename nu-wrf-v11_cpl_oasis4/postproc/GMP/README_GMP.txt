=========================================================
=========   Goddard Microphysics Processor V1.0 =========
=========================================================

 Oct 2018: Toshi Matsui @ NASA GSFC: Initial Version 1.0


The Goddard Microphysics Processor is developed upon the 
WRF's microphysics driver, and it essentially read the WRF
output and process one-time step of microphysics output. 
This framework is ideal for developing new code, and quickly 
testing your coding and output before plugging into the WRF/GCE. 

All directory is under GMP directory, and subdirectories are
explained as follow. 

--------- Subdirectories ---------

QRUN           : Running directory, where executable is copied in this directory. 
SRC            : Source code directory.   
INPUT          : Input directory (can be anywhere) 
OUTPUT         : Output directory (can be anywhere) 


========================================================
=========   Setting up and Running the GMP     =========
========================================================

0. Define MPI or single-CPU mode. 
  
   a. open SRC/define.h file
      either choose 0 or 2 for MPI option (default is 2)
     For quick development in earlier stage, 0 (no MPI) is fine. 

1. For NCCS Discover (NASA's super computer) users

   a. At the top of the directory (./), type
     >./build.sh -o cleanfirst gmp
      This will clean up first, and compile entire code. 

     If you modified the code, type
     >./build.sh -o rebuild gmp
      This will just re-compile the modified the code (faster for development)


   b. After compiling the code, go to QRUN
      open Config_GMP.F file, and edit run-time parameters.

   c. Then, run GMP

      For single CPU
      >./GMP.x

      For MPI run with batch job
      >sbatch Sbatch_####.sh
      Make sure you must edit Sbatch_####.sh for your own NCCS account. 

    d. Ouptut file can be dumped in the output directory with your own append name. 
       It is currently grads format. 

2. If you will work on your own computer...

   a. Go to SRC directory, and open make file.
      Modifye the following line based on your computer's compiler. 
      CPP is C-preprocessor, CF is Fortran compiler, CC is C compiler
      INC_NETCDF is NetCDF include file directory, and LD_NETCDF is 
      NetCDF library file directory. 

      CPP     =  /lib/cpp -traditional-cpp
      CF      = ifort
      CC      = icc
      INC_NETCDF = /path/to/netcdf/include
      LD_NETCDF = /path/to/netcdf/lib

   b. After editing make file, type
    >make
    for compilation. 

    >make clean
    for clean up. 

   c. After compiling the code, go to QRUN, 
      open Config_GMP.F file, and edit run-time parameters.

   d. Then, run GMP

      For single CPU, just type
      >./GMP.x

      For MPI run
      >mpirun -np 8 GMP.x

   e. Ouptut file can be dumped in the output directory with your own append name.
      It is currently grads format.

