!################################################################
!################  Define various C-preprocessor ################
!################################################################

!
! MPI OPTION (if MPI=1or2, you must have MPI library (MPICH or other), and specify MPI library lin in make file.) 
!   MPI  0 - no mpi
!   MPI  1 - file decomposition
!   MPI  2 - domain decomposition   
!
! -details of file or domain decomposition-
! MPI=1 is designed for file-loop decomposition. 
!   If you have a dozen (hudnred) of input files, it will be powerfull tool to save computational time.  
!   E.g., If you have 25 files, you can use up to 25 CPUs to gain the maximum speed. 
!
! MPI=2 is designed for domain decomposition. 
!   If you have a few of large-domain input files, it will be powerfull option to save computational time. 
!

# define MPI 0 

!
! HDF option (if HDF=1, you must have HDF library, and specify HDF library link in make file)
!    HDF 0 - does not include HDF library 
!    HDF 1 - include HDF library and its subroutine
!
! - details - 
! HDF option is linked to scan_micro or scan_radar simulator. If HDF=1 and scan_micro/scan_radar
! is turned on, HDF swath output will be generated. 
!

# define HDF 0 


