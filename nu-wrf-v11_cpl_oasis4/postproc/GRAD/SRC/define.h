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
# define MPI 2

!
! if defined some OpenMP funtion is called. If not, undefined (--> otherwise, it won't compile.)
!
# define OMP


!
! if defined, it account for GOCART aerosols. 
!
# undef WRF_CHEM

