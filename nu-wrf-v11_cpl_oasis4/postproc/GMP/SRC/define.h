!################################################################
!################  Define various C-preprocessor ################
!################################################################

!
! MPI OPTION (if MPI=1or2, you must have MPI library (MPICH or other), and specify MPI library lin in make file.)
!   MPI  0 - no mpi
!   MPI  2 - domain decomposition
!
! MPI=2 is designed for domain decomposition.
!   If you have a few of large-domain input files, it will be powerfull option to save computational time.
!
# define MPI 0

!
! if defined some OpenMP funtion is called. If not, undefined (--> otherwise, it won't compile.)
!
# define OMP


!
! if defined, it account for GOCART aerosols. 
!
# undef WRF_CHEM

# define WRF_CHEM 0
