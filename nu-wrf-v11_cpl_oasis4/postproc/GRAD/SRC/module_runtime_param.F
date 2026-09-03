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

 integer :: lw_physics      ! lw scheme index
 integer :: sw_physics      ! sw scheme index
 real :: dx  !horizontal grid space [m]
 integer :: goddardrad_2014_skip ! 0 (no skip) or 1 (skip)
!
! Choice of WRF physics scheme 
!
 integer,parameter :: & ! -> all (local)
   rrtmscheme=1      ,& ! rrtm lw scheme
   goddardswscheme_2011 =55  ,&
   goddardlwscheme_2011 =55  ,&
   goddardswscheme_2014 =56  ,&
   goddardlwscheme_2014 =56



 contains

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine read_runtime_param
 implicit none

 integer,parameter :: io = 101
 namelist /radiation_switch/  CRMdir, outdir, out_appendix, ids_wrf, ide_wrf, jds_wrf, jde_wrf, kmax, &
                              num_tiles_omp, lw_physics, sw_physics, dx, goddardrad_2014_skip

!
! Read configuration file
!
 open(io, file='Config_GoddardRad.F' )
 read(io,radiation_switch)
 close(io)

 return
 end subroutine read_runtime_param

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 end module module_runtime_param

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
