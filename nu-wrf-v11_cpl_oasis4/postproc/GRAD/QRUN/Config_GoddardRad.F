! #############################################################################
! #####################  Configure Radiation Switch  ##########################
! #############################################################################

$radiation_switch

  CRMdir = './../INPUT/'  ! Input CRM directory

  outdir = './../OUTPUT/' ! Output directory

  out_appendix = '_test'  !output file name appendix

  ids_wrf = 1 
  ide_wrf = 420 !423 is max

  jds_wrf = 1   ! 
  jde_wrf = 400 !411 is max
 
  kmax = 60 ! maximum vertical dimension

  num_tiles_omp = 1  ! number of tile for OMP/MIC

! Radiation scheme
! 1-rrtm radiation scheme  (not available now)
! 55-gsfc radiation scheme (2011 version) 
! 56-gsfc radiation scheme (2014 version)
  lw_physics=56       ! lw scheme index
  sw_physics=56       ! sw scheme index

  dx = 5000.0         !  horizontal grid space (m)

  goddardrad_2014_skip = 0  ! 0 (non skip) or 1 (skip)

$end

