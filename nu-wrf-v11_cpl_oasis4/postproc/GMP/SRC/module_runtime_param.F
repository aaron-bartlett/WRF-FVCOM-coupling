! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 module module_runtime_param

 implicit none

 public

 save

 character(len=200) :: CRMdir  ! Input CRM directory
 character(len=200) :: outdir  ! Input CRM directory
 character(len=100) :: out_appendix  ! output filename appendix
 integer :: ids_wrf, ide_wrf, jds_wrf, jde_wrf, kmax  ! maximum i j k bounds of WRF
 integer :: num_tiles_omp   ! the number of tiles for OMP/MIC

 integer :: mp_physics      ! microphysics scheme index
 real :: dx  !horizontal grid space [m]

!
! Choice of WRF physics scheme 
!
 integer,parameter :: & ! -> all (local)
   NUWRF4ICESCHEME = 56  ,&
   GMPUNIFIED      = 58  

 contains

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine read_runtime_param
 implicit none

 integer,parameter :: io = 101
 namelist /mp_switch/  CRMdir, outdir, out_appendix, ids_wrf, ide_wrf, jds_wrf, jde_wrf, kmax, &
                              num_tiles_omp, mp_physics, dx

!
! Read configuration file
!
 open(io, file='Config_GMP.F' )
 read(io,mp_switch)
 close(io)

 return
 end subroutine read_runtime_param

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 end module module_runtime_param

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
