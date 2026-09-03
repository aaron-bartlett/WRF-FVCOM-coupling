! #############################################################################
! #####################  Configure Radiation Switch  ##########################
! #############################################################################
$mp_switch

  CRMdir = '/gpfsm/dnb04/projects/p06/data/NUWRF/TWP_Jan23/4ICE/'  ! Input WRF directory

  outdir = './../OUTPUT/' ! Output directory

  dx = 1000.0         ! horizontal grid space (m)

  ids_wrf = 300 
  ide_wrf = 333  !333 !is max

  jds_wrf = 10   ! 
  jde_wrf = 33 !333 !333 is max
 
  kmax = 60 ! maximum vertical dimension

  num_tiles_omp = 1  ! number of tile for OMP/MIC

!  out_appendix = '_NUWRF_4ICE'
  out_appendix = '_GMP_UNIFIED'  

!  mp_physics=56       ! 56 GCE 4ICE microphysics scheme index
  mp_physics=58       ! 58 Goddard Unified Microphysics (GMPUNIFIED)

$end
