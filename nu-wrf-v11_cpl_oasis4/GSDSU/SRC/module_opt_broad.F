#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_opt_broad
 use module_simulator
 use module_opt_visir
 use module_surface
 use module_gmp_4ice
 use module_icescat_visir
 use module_mpi
 use module_simulator_broad
 implicit none

!---------------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!
!
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
! Comments:  
!   This module deal with broad-band single-scattering properties for particles. 
!   Single-scattering LUTs are pre-computed by integrating spectrum parameters using visir routines.
!   So, essentially assumptions and methods (particles) are consistent between visir and broadband 
!   simulators.   
! 
! History: 
!  05/2010  Toshi Matsui@NASA GSFC: Found bug in computing ebudget (corrected)
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 save     ! all module parameters will be saved

!
! Encapsulation control 
!
  private   ! encapsulate all variables and subourtines (non accessible)

!
! public (accessible) subroutines
!
  public :: makelut_broad, &          ! create broad-band single-scattering LUTs
            makelut_broad_aerosol, &  ! create broad-band single-scattering LUTs for aerosols
                opt_broad    ! compute broad-band single-scattering properties

  save 

! Goddard Radiation (Chou-Suarez) band
!
!        in the uv region :
!           index  1 for the 0.175-0.225 micron band
!           index  2 for the 0.225-0.245; 0.260-0.280 micron band (tricky here...)
!           index  3 for the 0.245-0.260 micron band
!           index  4 for the 0.280-0.295 micron band
!           index  5 for the 0.295-0.310 micron band
!           index  6 for the 0.310-0.320 micron band
!           index  7 for the 0.325-0.400 micron band
!        in the par region :
!           index  8 for the 0.400-0.700 micron band
!        in the infrared region :
!           index  9 for the 0.700-1.220 micron band
!           index 10 for the 1.220-2.270 micron band
!           index 11 for the 2.270-10.00 micron band
!
!   band     wavenumber (/cm)   absorber
!
!    1           0 - 340           h2o
!    2         340 - 540           h2o
!    3         540 - 800       h2o,cont,co2
!    4         800 - 980       h2o,cont
!                              co2,f11,f12,f22
!    5         980 - 1100      h2o,cont,o3
!                              co2,f11
!    6        1100 - 1215      h2o,cont
!                              n2o,ch4,f12,f22
!    7        1215 - 1380      h2o,cont
!                              n2o,ch4
!    8        1380 - 1900          h2o
!    9        1900 - 3000          h2o
!
! in addition, a narrow band in the 17 micrometer region is added to
!    compute flux reduction due to n2o
!
!    10        540 - 620       h2o,cont,co2,n2o
!
! band 3 (540-800/cm) is further divided into 3 sub-bands :
!
!   subband   wavenumber (/cm)
!
!    1          540 - 620
!    2          620 - 720
!    3          720 - 800
!


! bundled scattering array (mxpts_re_visir, max_spc, nband, mxisl )
 real(sdsu_fps),allocatable :: kext_q_unit_array (:,:,:,:) ! extinction cross section per unit mass [km-1 / (g/m3) ]
 real(sdsu_fps),allocatable :: omega_array       (:,:,:,:) ! signle scattering albedo [-]
 real(sdsu_fps),allocatable :: asym_array        (:,:,:,:) ! asymetry parameter (forward scattering fraction) [-]

 real(sdsu_fps),allocatable :: kext_q_unit_aero_array (:,:,:,:) ! extinction cross section per unit mass [km-1 / (g/m3) ]
 real(sdsu_fps),allocatable :: omega_aero_array       (:,:,:,:) ! signle scattering albedo [-]
 real(sdsu_fps),allocatable :: asym_aero_array        (:,:,:,:) ! asymetry parameter (forward scattering fraction) [-]


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 subroutine opt_broad
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation method.
!    
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; Revised LUT.
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 integer :: i,j,k !domain loop index
 integer :: ib,n  !band loop index
 integer :: mxband !maximum band #
 integer :: ierr  !allocation stat
 integer :: max_band, mxib

 integer :: ispc  !RAMS spicies index
 real(sdsu_fps) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation

 type ( particle_gen ) :: & !particle_gen is defined in module_simulater
     kext_gen,  & ! extinction coefficient [km-1]
     salb_gen,  & ! single scattering albedo [-]
     asym_gen     ! asymetry parameter [-]

 type ( particle_rams ) :: & !particle_rams is defined in module_simulater
     kext_rams,  & ! extinction coefficient [km-1]
     salb_rams,  & ! single scattering albedo [-]
     asym_rams     ! asymetry parameter [-]

 type ( particle_sbm ) :: & !particle_sbm is defined in module_simulater
     kext_sbm,  & ! extinction coefficient [km-1]
     salb_sbm,  & ! single scattering albedo [-]
     asym_sbm     ! asymetry parameter [-]

 type ( particle_gocart ) :: & !particle_gocart is defined in module_simulater
     kext_gocart,  & ! extinction coefficient [km-1]
     salb_gocart,  & ! single scattering albedo [-]
     asym_gocart     ! asymetry parameter [-]

 type ( particle_ucd ) :: & !particle_ucd is defined in module_simulater
     kext_ucd,  & ! extinction coefficient [km-1]
     salb_ucd,  & ! single scattering albedo [-]
     asym_ucd     ! asymetry parameter [-]


!
! single-scattering properties for condensates (cloud-precipitation)
!
 real(sdsu_fps) :: & ! 
   kextcl, & ! total extinction [km-1]
   salbcl, & ! total single scatterling albedo [-]
   asymcl    ! total asymetery parameters [-]

 real(sdsu_fps),dimension(mxlyr,nband_sw) :: & ! 
   kextcl_sw, & ! total extinction [km-1]
   salbcl_sw, & ! total single scatterling albedo [-]
   asymcl_sw    ! total asymetery parameters [-]

 real(sdsu_fps),dimension(mxlyr,nband_lw) :: & ! 
   kextcl_lw, & ! total extinction [km-1]
   salbcl_lw, & ! total single scatterling albedo [-]
   asymcl_lw    ! total asymetery parameters [-]

 real(sdsu_fps),dimension(mxlyr,nband_aw) :: & ! 
   kextcl_aw, & ! total extinction [km-1]
   salbcl_aw, & ! total single scatterling albedo [-]
   asymcl_aw    ! total asymetery parameters [-]


!
! single-scattering properties for aerosols
!
 real(sdsu_fps) :: & ! 
   kextal, & ! total extinction coef [km-1]
   salbal, & ! total single scatterling albedo [-]
   asymal    ! total asymetery parameters [-]

 real(sdsu_fps),dimension(mxlyr,nband_sw) :: & ! 
   kextal_sw, & ! total extinction [km-1]
   salbal_sw, & ! total single scatterling albedo [-]
   asymal_sw    ! total asymetery parameters [-]

 real(sdsu_fps),dimension(mxlyr,nband_lw) :: & ! 
   kextal_lw, & ! total extinction [km-1]
   salbal_lw, & ! total single scatterling albedo [-]
   asymal_lw    ! total asymetery parameters [-]

 real(sdsu_fps),dimension(mxlyr,nband_aw) :: & ! 
   kextal_aw, & ! total extinction [km-1]
   salbal_aw, & ! total single scatterling albedo [-]
   asymal_aw    ! total asymetery parameters [-]

!
! misc
!
 real(sdsu_fps),dimension(mxlyr) :: thrate_sw, thrate_lw !heating rate [K/sec]
 real(sdsu_fps),dimension(1:4) :: ERBE_sw, ERBE_lw  !radiaiton budget
!                    (LW 1-TOA LW down, 2-TOA LW up, 3-surface LW down, 4-surface LW up)
!                    (SW 1-TOA SW down, 2-TOA SW up, 3-surface SW down, 4-surface SW up)

 real(sdsu_fps) :: swdown_lsm(1:3,1:2) !1-UV, 2-PAR, 3-NIR broadband radiation flux fraction                                                   !for beam and diffuse [W/m2]


 real(sdsu_fps) :: kext_cloud  ! cloud only extinction [km-1]
 real(sdsu_fps) :: kext_cloud_par(mxlyr)


!============================================================================================

 if(masterproc) print*,'MSG opt_broad: Start BroadBand Simulator'

 mxband = max(nband_sw, nband_lw)

!
! read/initialize single-scattering look-up table
!
! call read_sslut(mxband) 

 call read_sslut_nc(mxband)

 if(account_aerosol) call read_sslut_aero_nc(mxband)
 

!
! domain loop
!
 DOM_LOOP: do j = myj_start, myj_end ; do i = myi_start, myi_end

!
! interpolate sounding climatology for a given latitude and julian date
!
 call sounding_interp(surface(i,j)%lat, INT(sdsu_julian))

!
! Calculate single-scattering properties of each grid at every SW LW band 
!

  mxib = 2  !1-sw, 2-lw

   RAD_LOOP: do ib = 1, mxib

       if(ib == 1) then !SW
         max_band = nband_sw
       else              !LW
         max_band = nband_lw
       endif

   BAND_LOOP: do n = 1, max_band


   VERTICAL_LOOP: do k = 1, mxlyr

   !
   ! cloud/rain particles
   ! 
   mic_select1: select case(trim(type_microphysics))
   case('GEN')

       !
       ! compute optical properties for different species of GCE microphysics
       !
       call opt_gen_broad( n, ib, q_gen(i,j,k), re_gen(i,j,k), kext_gen, salb_gen, asym_gen )

       !
       ! Total optical properties
       !
       call total_opt_gen_broad( kext_gen, salb_gen, asym_gen, kextcl, salbcl, asymcl , kext_cloud ) 


   case('RAMS')

       !
       ! compute optical properties for different species
       !
       call opt_rams_broad( n, ib, atmos(i,j,k)%t_air, atmos(i,j,k)%rh , q_rams(i,j,k), re_rams(i,j,k), &
                            kext_rams, salb_rams, asym_rams )

       !
       ! Total optical properties
       !
       call total_opt_rams_broad( kext_rams, salb_rams, asym_rams, kextcl, salbcl, asymcl, kext_cloud )


   case('SBM')

      call stop_sdsu(myrank,'Not yet supporting SBM')

   case default 

      kextcl  = 0.e0
      salbcl  = 0.e0
      asymcl  = 0.e0

   end select mic_select1

     if(ib == 1) then !SW
        kextcl_sw(k,n) = kextcl
        salbcl_sw(k,n) = salbcl
        asymcl_sw(k,n) = asymcl
     elseif(ib == 2) then !LW
        kextcl_lw(k,n) = kextcl
        salbcl_lw(k,n) = salbcl
        asymcl_lw(k,n) = asymcl
     endif

   !
   ! cloud (only) optical depth for visible range 
   !
   if( ib == 1 .and. n == 8) then
       kext_cloud_par(k) = kext_cloud 
   endif


   !
   ! aerosol particles
   !
      kextal  = 0.e0
      salbal  = 0.e0
      asymal  = 0.e0

   if(account_aerosol) then

    select case(trim(aerosol_microphysics) )
    case('GOCART')

       !
       ! compute optical properties for different species of GCE microphysics
       !
       call opt_gocart_broad( n, ib, q_gocart(i,j,k), atmos(i,j,k)%rh*1.e-2, &
                              kext_gocart, salb_gocart, asym_gocart )

       !
       ! Total optical properties
       !
       call total_opt_gocart_broad( kext_gocart, salb_gocart, asym_gocart, kextal, salbal, asymal )

    case('UCD')

       !
       ! compute optical properties for different species of GCE microphysics
       !
       call opt_ucd_broad( n, ib, q_ucd(i,j,k), atmos(i,j,k)%rh*1.e-2, &
                              kext_ucd, salb_ucd, asym_ucd )

       !
       ! Total optical properties
       !
       call total_opt_ucd_broad( kext_ucd, salb_ucd, asym_ucd, kextal, salbal, asymal )


    end select


   endif

     if(ib == 1) then !SW
        kextal_sw(k,n) = kextal
        salbal_sw(k,n) = salbal
        asymal_sw(k,n) = asymal
     elseif(ib == 2) then !LW
        kextal_lw(k,n) = kextal
        salbal_lw(k,n) = salbal
        asymal_lw(k,n) = asymal
     endif

   enddo VERTICAL_LOOP

   enddo BAND_LOOP

   enddo RAD_LOOP


! initialize
 ERBE_sw = 0.e0 ; ERBE_lw = 0.e0

  !
  ! 1D vertical column driver of broadband radiative transfer (SW)
  !
  call  simulator_broad( 'SW', nband_sw, mxlyr, surface(i,j)%cosz, surface(i,j)%t_skin, surface(i,j)%t_air , &
                  emiss_broad(i,j), albedo_broad(i,j,1:2,1:2), sdsu_solcon, &
                  atmos(i,j,1:mxlyr)%t_air, atmos(i,j,1:mxlyr)%press, atmos_stag(i,j,0:mxlyr)%press, &
                  atmos(i,j,1:mxlyr)%exner, atmos(i,j,1:mxlyr)%dhgt , atmos(i,j,1:mxlyr)%sh,    &
                  kextcl_sw, salbcl_sw, asymcl_sw, &
                  kextal_sw, salbal_sw, asymal_sw, &
                  ERBE_sw, thrate_sw(1:mxlyr), swdown_lsm , &
                  aod_broad(i,j), cod_broad(i,j)  ) ! output


  !
  ! 1D vertical column driver of broadband radiative transfer (LW)
  !
  call  simulator_broad( 'LW', nband_lw, mxlyr, surface(i,j)%cosz, surface(i,j)%t_skin, surface(i,j)%t_air , &
                  emiss_broad(i,j), albedo_broad(i,j,1:2,1:2), sdsu_solcon, &
                  atmos(i,j,1:mxlyr)%t_air, atmos(i,j,1:mxlyr)%press, atmos_stag(i,j,0:mxlyr)%press, &
                  atmos(i,j,1:mxlyr)%exner, atmos(i,j,1:mxlyr)%dhgt , atmos(i,j,1:mxlyr)%sh,    &
                  kextcl_lw, salbcl_lw, asymcl_lw, &
                  kextal_lw, salbal_lw, asymal_lw, &
                  ERBE_lw, thrate_lw(1:mxlyr)             ) ! output


  !
  ! SDSU output (SW)
  !
  ebudget(i,j,1) = ERBE_sw(1)  ! TOA SW downwelling flux [W/m2]
  ebudget(i,j,2) = ERBE_sw(2)  ! TOA SW upwelling flux   [W/m2]
  ebudget(i,j,3) = ERBE_sw(3)  ! Surface SW downwelling flux [W/m2]
  ebudget(i,j,4) = ERBE_sw(4)  ! Surface SW upwelling flux [W/m2]
  ebudget(i,j,5) = ERBE_sw(1)  + ERBE_sw(2)   !net toa radiation [W/m2]
  ebudget(i,j,6) = ERBE_sw(3)  + ERBE_sw(4)   !net surface radiation [W/m2]
  ebudget(i,j,7) = ebudget(i,j,5)   - ebudget(i,j,6)    !net atmosphere radiation [W/m2] 
  if(heating_rate) sw_heat(i,j,1:mxlyr) = thrate_sw(1:mxlyr) * 3600. * 24. ! [K/day] <- [K/sec]


  !
  ! SDSU output (LW)
  !
  ebudget(i,j,8)  = ERBE_lw(1)  ! TOA LW downwelling flux [W/m2]
  ebudget(i,j,9)  = ERBE_lw(2)  ! TOA LW upwelling flux   [W/m2]
  ebudget(i,j,10) = ERBE_lw(3)  ! Surface LW downwelling flux [W/m2]
  ebudget(i,j,11) = ERBE_lw(4)  ! Surface LW upwelling flux [W/m2]
  ebudget(i,j,12) = ERBE_lw(1)  + ERBE_lw(2)   !net toa radiation [W/m2]
  ebudget(i,j,13) = ERBE_lw(3)  + ERBE_lw(4)   !net surface radiation [W/m2]
  ebudget(i,j,14) = ebudget(i,j,12)  - ebudget(i,j,13)   !net atmosphere radiation [W/m2] 
  if(heating_rate) lw_heat(i,j,1:mxlyr) = thrate_lw(1:mxlyr) * 3600. * 24. ! [K/day] <- [K/sec]

  !
  ! cloud only optical depth [-]
  !
  cod_broad(i,j) =  sum( atmos(i,j,1:mxlyr)%dhgt *  kext_cloud_par(1:mxlyr) )

 enddo ; enddo DOM_LOOP


#if MPI == 2 
if(numproc_tot>1) then
! 
! Gather all outputs sub-domain from different slave processors into the master processor. 
!
  if(numproc_tot>1) call mpi_sdsu_communicate('TO_MASTER', aod_broad )
  if(numproc_tot>1) call mpi_sdsu_communicate('TO_MASTER', cod_broad )
  if(numproc_tot>1) call mpi_sdsu_communicate('TO_MASTER', ebudget )
  if(numproc_tot>1 .and. heating_rate) call mpi_sdsu_communicate('TO_MASTER', lw_heat )
  if(numproc_tot>1 .and. heating_rate) call mpi_sdsu_communicate('TO_MASTER', sw_heat )

endif
#endif


  return
 end subroutine opt_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_sslut_nc(mxband)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!   Read NetCDF LUT data.
!    
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 integer,intent(in) :: mxband !maximum band

 integer :: i,j,k !domain loop index
 integer :: ib,n  !band loop index
 integer :: ierr  !allocation stat
 integer :: max_band
 integer :: mxisl, max_spc, nband
 integer :: ncid, varid


!
! static indice
!
 mxisl = 2  !1-sw, 2-lw

 select case(trim(type_microphysics))
 case('GEN')    ! GCE (LIN_type) 1moment bulk microphsics 
   max_spc = 6
 case('RAMS')  !GCE RAMS 1- and 2-moment microphysics 
   max_spc = 16
 case('SBM')
   call stop_sdsu(myrank,'Not yet supporoting SBM')
 case default
 end select

 nband = max( nband_sw, nband_lw )

!
! allocate array parameters
!
 if( .not. allocated(kext_q_unit_array) ) then

   allocate( kext_q_unit_array ( mxpts_re_visir, max_spc, nband, mxisl ), &
             omega_array       ( mxpts_re_visir, max_spc, nband, mxisl ), &
             asym_array        ( mxpts_re_visir, max_spc, nband, mxisl ), &
             stat=ierr )
   if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_broad: allocation error -> Terminate program.')

 endif


!
! open netcdf file
!
 sdsu_io_file = trim(sdsu_dir_sslut)//'BROADBAND_CLOUD_'//trim(cloud_microphysics)//'.nc'
 if(masterproc) print*, 'MSG read_sslut_nc: reading broadband LUT:',trim(sdsu_io_file)
 
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )


!
! extinction cross section per unit mass [km-1 / (g/m3) ]
!
 call check( nf90_inq_varid(ncid, 'kext_q_unit' , varid ) )
 call check( nf90_get_var(ncid, varid, kext_q_unit_array ))

!
! signle scattering albedo [-] 
!
 call check( nf90_inq_varid(ncid, 'omega' , varid ) )
 call check( nf90_get_var(ncid, varid, omega_array ))

!
! asymetry parameter (forward scattering fraction) [-] 
!
 call check( nf90_inq_varid(ncid, 'asym' , varid ) )
 call check( nf90_get_var(ncid, varid, asym_array ))

!
! close netCDF file
!
  call check( nf90_close(ncid) )



 return
 end subroutine read_sslut_nc

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_sslut_aero_nc(mxband)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!   Read NetCDF LUT data for aerosols.
!    
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 integer,intent(in) :: mxband !maximum band

 integer :: i,j,k !domain loop index
 integer :: ib,n  !band loop index
 integer :: ierr  !allocation stat
 integer :: max_band
 integer :: mxisl, max_spc, nband
 integer :: ncid, varid


!
! static indice
!
 mxisl = 2  !1-sw, 2-lw

 select case(trim(aerosol_microphysics) )
 case('GOCART')
   max_spc = 14
 case('UCD')
   max_spc = 5
 end select

 nband = max( nband_sw, nband_lw )

!
! allocate array parameters
!
 if( .not. allocated(kext_q_unit_array) ) then

   allocate( kext_q_unit_aero_array ( mxpts_rh, max_spc, nband, mxisl ), &
             omega_aero_array       ( mxpts_rh, max_spc, nband, mxisl ), &
             asym_aero_array        ( mxpts_rh, max_spc, nband, mxisl ), &
             stat=ierr )
   if (ierr /= 0) call stop_sdsu(myrank,'MSG read_sslut_aero_nc: allocation error -> Terminate program.')

 endif


!
! open netcdf file
!
 sdsu_io_file = trim(sdsu_dir_sslut)//'BROADBAND_AEROSOL_'//trim(aerosol_microphysics)//'.nc'
 if(masterproc) print*, 'MSG read_sslut_aero_nc: reading broadband LUT:',trim(sdsu_io_file)

 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )


!
! extinction cross section per unit mass [km-1 / (g/m3) ]
!
 call check( nf90_inq_varid(ncid, 'kext_q_unit' , varid ) )
 call check( nf90_get_var(ncid, varid, kext_q_unit_aero_array ))

!
! signle scattering albedo [-] 
!
 call check( nf90_inq_varid(ncid, 'omega' , varid ) )
 call check( nf90_get_var(ncid, varid, omega_aero_array ))

!
! asymetry parameter (forward scattering fraction) [-] 
!
 call check( nf90_inq_varid(ncid, 'asym' , varid ) )
 call check( nf90_get_var(ncid, varid, asym_aero_array ))

!
! close netCDF file
!
  call check( nf90_close(ncid) )


 return
 end subroutine read_sslut_aero_nc

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mie_lut_broad( q , re, kext_lut, salb_lut, asym_lut, pts_ref, mxpts_ref, &
                           kext, salb, asym )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the extinction/scattering coefficients for
! given Re using Mie LUTs (per unit q), 
! and optal values are re-scaled to actual q (hydrometeor content).
! 
! History:
!   05/2009  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(in) :: q  ! particle mass mixing ratio [g/m3]
 real(sdsu_fps) :: re            ! particle effective radius [micron]
 real(sdsu_fps),intent(in) :: kext_lut(mxpts_ref) ! extinction coefficient per unit mass [km-1 / (g/m3) ]
 real(sdsu_fps),intent(in) :: salb_lut(mxpts_ref) ! single scattering albedo [-]
 real(sdsu_fps),intent(in) :: asym_lut(mxpts_ref) ! asymetry parameter [-]
 real(sdsu_fps),intent(in) :: pts_ref(mxpts_ref)  ! LUT of re index [micron]
 integer,intent(in) :: mxpts_ref           ! dimension of re index
 real(sdsu_fps),intent(out) :: kext   ! extinction coefficient [km-1] 
 real(sdsu_fps),intent(out) :: salb   ! single scattering albedo [-]
 real(sdsu_fps),intent(out) :: asym   ! asymetry parameter [-]

! local parameter
 integer :: j !loop index
 real(sdsu_fps) :: kext_q_unit   ! extinction coefficient per unit mass [km-1 / (g/m3) ]
 real(sdsu_fps) :: wgt1,wgt2

!
!initialise
!
 kext = 0.d0 ; salb = 0.d0 ; asym = 0.d0

 if(q <= q_min_condensate) return

 if(re < pts_ref(1)        ) re = pts_ref(1)
 if(re > pts_ref(mxpts_ref)) re = pts_ref(mxpts_ref)

!
! Linear interpolation
!
 REF_LOOP: do j = 1,mxpts_ref-1

  if( re >= pts_ref(j) .and. re <= pts_ref(j+1) ) then
      wgt2 =  ( re - pts_ref(j) )    / ( pts_ref(j+1) - pts_ref(j) )
      wgt1 = 1.0 - wgt2

      kext_q_unit = wgt1*kext_lut(j) + wgt2*kext_lut(j+1)
      salb        = wgt1*salb_lut(j) + wgt2*salb_lut(j+1)
      asym        = wgt1*asym_lut(j) + wgt2*asym_lut(j+1)

      exit REF_LOOP
  endif

 enddo REF_LOOP


!
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
! 
  kext = kext_q_unit * q 

 return
 end subroutine mie_lut_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mie_lut_broad_aerosol( q , rh_in, kext_lut, salb_lut, asym_lut, pts_rh, mxpts, &
                                  kext, salb, asym )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the extinction/scattering coefficients for
! given Re using Mie LUTs (per unit q), 
! and optal values are re-scaled to actual q (hydrometeor content).
! 
! History:
!   05/2009  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(in) :: q  ! particle mass mixing ratio [g/m3]
 real(sdsu_fps),intent(in) :: rh_in ! relative humidity [-] (0~1) 
 real(sdsu_fps),intent(in) :: kext_lut(mxpts) ! extinction coefficient per unit mass [km-1 / (g/m3) ]
 real(sdsu_fps),intent(in) :: salb_lut(mxpts) ! single scattering albedo [-]
 real(sdsu_fps),intent(in) :: asym_lut(mxpts) ! asymetry parameter [-]
 real(sdsu_fps),intent(in) :: pts_rh(mxpts)  ! LUT of re index [micron]
 integer,intent(in) :: mxpts          ! dimension of rh index
 real(sdsu_fps),intent(out) :: kext   ! extinction coefficient [km-1] 
 real(sdsu_fps),intent(out) :: salb   ! single scattering albedo [-]
 real(sdsu_fps),intent(out) :: asym   ! asymetry parameter [-]

! local parameter
 integer :: j                    ! loop index
 real(sdsu_fps) :: kext_q_unit   ! extinction coefficient per unit mass mixing ratio [km-1 / (g/m3) ]
 real(sdsu_fps) :: wgt1,wgt2
 real(sdsu_fps) :: rh            ! relative humidity [-] (0~1) 

!
!initialise
!
 kext = 0.d0 ; salb = 0.d0 ; asym = 0.d0
 rh = rh_in

 if(q <= q_min_aerosol) return

 if(rh < pts_rh(1)    ) rh = pts_rh(1)
 if(rh > pts_rh(mxpts)) rh = pts_rh(mxpts_rh)

!
! Linear interpolation
!
 RH_LOOP: do j = 1,mxpts_rh-1

  if( rh >= pts_rh(j) .and. rh <= pts_rh(j+1) ) then
      wgt2 =  ( rh - pts_rh(j) )    / ( pts_rh(j+1) - pts_rh(j) )
      wgt1 = 1.0 - wgt2

      kext_q_unit = wgt1*kext_lut(j) + wgt2*kext_lut(j+1)
      salb        = wgt1*salb_lut(j) + wgt2*salb_lut(j+1)
      asym        = wgt1*asym_lut(j) + wgt2*asym_lut(j+1)

      exit RH_LOOP
  endif

 enddo RH_LOOP


!
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
! 
  kext = kext_q_unit * q

 return
 end subroutine mie_lut_broad_aerosol

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  opt_gen_broad( n, ib, q, re, kext_gen, salb_gen, asym_gen )
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation techniqu.
!  This routine is for GCE microphysics.
!    
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; Revised completely for faster processing. 
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 implicit none
 integer,intent(in)::n  ! radiation band number
 integer,intent(in)::ib ! radiation type (1:SW or 2:LW)
 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 type ( particle_gen ),intent(out) :: & ! particle_gen is defined in module_simulater
     kext_gen,  & ! extinction coefficient [km-1]
     salb_gen,  & ! single scattering albedo [-]
     asym_gen     ! asymetry parameter [-]

 integer :: j
 integer :: idx_spc
 real(sdsu_fps) :: wgt1, wgt2
 real(sdsu_fps) :: re_in, q_in
 real(sdsu_fps) :: kext_q_unit , salb , asym

!
! initialize
! 
 call init_microphysics( kext_gen, 'zero')
 call init_microphysics( salb_gen, 'zero')
 call init_microphysics( asym_gen, 'zero')


!
! cloud species
!
 idx_spc = 1  ; q_in = q%cloud ; re_in = re%cloud

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit 
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gen%cloud = kext_q_unit * q_in ; salb_gen%cloud = salb ; asym_gen%cloud = asym  !outpout

!
! rain species
!
 idx_spc = 2  ; q_in = q%rain ; re_in = re%rain

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gen%rain = kext_q_unit * q_in ; salb_gen%rain = salb ; asym_gen%rain = asym  !outpout


!
! ice species
!
 idx_spc = 3  ; q_in = q%ice ; re_in = re%ice

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gen%ice = kext_q_unit * q_in ; salb_gen%ice = salb ; asym_gen%cloud = asym  !outpout


!
! snow species
!
 idx_spc = 4  ; q_in = q%snow ; re_in = re%snow

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gen%snow = kext_q_unit * q_in ; salb_gen%snow = salb ; asym_gen%snow = asym  !outpout

 
!
! graupel species
!
 idx_spc = 5  ; q_in = q%graupel ; re_in = re%graupel

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gen%graupel = kext_q_unit * q_in ; salb_gen%graupel = salb ; asym_gen%graupel = asym  !outpout



!
! hail species
!
 idx_spc = 1  ; q_in = q%hail ; re_in = re%hail

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gen%hail = kext_q_unit * q_in ; salb_gen%hail = salb ; asym_gen%hail = asym  !outpout


 return
 end subroutine opt_gen_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine total_opt_gen_broad( kext_gen, salb_gen, asym_gen, &
                                 kexttot, salbtot, asymtot , kext_cloud )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute total single-scattering properties from atmospheric particles.
!  This routine is for GCE microphysics.
!  
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; Add cloud-only optical depth
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     kext_gen,  & ! extinction coefficient [km-1]
     salb_gen,  & ! single scattering albedo [-]
     asym_gen     ! asymetry parameter [-]

 real(sdsu_fps),intent(out) ::  kexttot, salbtot, asymtot !total optical parameters
 real(sdsu_fps),intent(out) ::  kext_cloud !cloud component optical parameters

!
! total extinction 
!
   kexttot    = &
                kext_gen%cloud  + &
                kext_gen%rain   + &
                kext_gen%ice    + &
                kext_gen%snow   + &
                kext_gen%graupel+ &
                kext_gen%hail

   kext_cloud  = &
                kext_gen%cloud  + &
                kext_gen%ice   


!
! total single scattering albedo
!
   if ( kexttot <= 0.e0 ) then
        salbtot = 0.e0
   else
       salbtot = ( &
                  salb_gen%cloud  * kext_gen%cloud   + &
                  salb_gen%rain   * kext_gen%rain    + &
                  salb_gen%ice    * kext_gen%ice     + &
                  salb_gen%snow   * kext_gen%snow    + &
                  salb_gen%graupel* kext_gen%graupel + &
                  salb_gen%hail   * kext_gen%hail      &
                 ) / kexttot

   endif

!
! total asymetry parameter
!
   if ( salbtot <= 0.e0 ) then
       asymtot = 0.e0
   else
       asymtot = ( &
                  asym_gen%cloud  * salb_gen%cloud  * kext_gen%cloud   + &
                  asym_gen%rain   * salb_gen%rain   * kext_gen%rain    + &
                  asym_gen%ice    * salb_gen%ice    * kext_gen%ice     + &
                  asym_gen%snow   * salb_gen%snow   * kext_gen%snow    + &
                  asym_gen%graupel* salb_gen%graupel* kext_gen%graupel + &
                  asym_gen%hail   * salb_gen%hail   * kext_gen%hail      &
                 ) / ( salbtot * kexttot )

   endif

 return
 end subroutine total_opt_gen_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  opt_rams_broad( n, ib, t_air, rh, q, re, kext_rams, salb_rams, asym_rams )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation techniqu.
!  This routine is for RAMS microphysics
!    
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; Revised completely for faster processing. 
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 integer,intent(in)::n  !band #
 integer,intent(in)::ib !radiation type
 real(sdsu_fps),intent(in) :: t_air !air temperature [K]
 real(sdsu_fps),intent(in) :: rh    !relative humidity 
 type ( particle_rams ),intent(in) :: & ! particle_rams is defined in module_simulater
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 type ( particle_rams ),intent(out) :: & ! particle_rams is defined in module_simulater
     kext_rams,  & ! extinction coefficient [km-1]
     salb_rams,  & ! single scattering albedo [-]
     asym_rams     ! asymetry parameter [-]

 integer :: ns,nt !temperature and RH index
 integer :: ispc  !species index
 integer :: j
 integer :: idx_spc
 real(sdsu_fps) :: wgt1, wgt2
 real(sdsu_fps) :: re_in, q_in
 real(sdsu_fps) :: kext_q_unit , salb , asym



!
! initialize
!
 call init_microphysics( kext_rams, 'zero')
 call init_microphysics( salb_rams, 'zero')
 call init_microphysics( asym_rams, 'zero')


!
! cloud1 species
!
 idx_spc = 1  ; q_in = q%cloud1 ; re_in = re%cloud1

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_rams%cloud1 = kext_q_unit * q_in ; salb_rams%cloud1 = salb ; asym_rams%cloud1 = asym  !outpout


!
! cloud2 species
!
 idx_spc = 2  ; q_in = q%cloud2 ; re_in = re%cloud2

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_rams%cloud2 = kext_q_unit * q_in ; salb_rams%cloud2 = salb ; asym_rams%cloud2 = asym  !outpout

!
! rain species
!
 idx_spc = 3  ; q_in = q%rain ; re_in = re%rain

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_rams%rain = kext_q_unit * q_in ; salb_rams%rain = salb ; asym_rams%rain = asym  !outpout

!
! ice1 species
!
 nt = max( 1,min( 31,-nint( (t_air-273.15) )))  ! Temperature index--1 to 31 
 ns = max(79,min(100, nint( rh  )))  ! RH index--79 to 100 
 ispc = rams_jhabtab(nt,ns,1) !ice1 habit

 idx_spc = ispc  ; q_in = q%ice1 ; re_in = re%ice1

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_rams%ice1 = kext_q_unit * q_in ; salb_rams%ice1 = salb ; asym_rams%ice1 = asym  !outpout


!
! ice2 species
!
 ispc = rams_jhabtab(nt,ns,2) !ice2 habit
 idx_spc = ispc  ; q_in = q%ice2 ; re_in = re%ice2

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_rams%ice2 = kext_q_unit * q_in ; salb_rams%ice2 = salb ; asym_rams%ice2 = asym  !outpout


!
! snow species
!
 idx_spc = 14  ; q_in = q%snow ; re_in = re%snow

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_rams%snow = kext_q_unit * q_in ; salb_rams%snow = salb ; asym_rams%snow = asym  !outpout

!
! graupel species
!
 idx_spc = 15  ; q_in = q%graupel ; re_in = re%graupel

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_rams%graupel = kext_q_unit * q_in ; salb_rams%graupel = salb ; asym_rams%graupel = asym  !outpout

!
! hail species
!
 idx_spc = 16  ; q_in = q%hail ; re_in = re%hail

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    if(re_in < pts_re_visir(1)        ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 =  ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) )
         wgt1 = 1.0 - wgt2

         kext_q_unit =  wgt1*kext_q_unit_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
    enddo
 endif

! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_rams%hail = kext_q_unit * q_in ; salb_rams%hail = salb ; asym_rams%hail = asym  !outpout


 return
 end subroutine opt_rams_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine total_opt_rams_broad( kext_rams, salb_rams, asym_rams, &
                                 kexttot, salbtot, asymtot , kext_cloud)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute total single-scattering properties from atmospheric particles.
!  This routine is for RAMS microphysics.
!    
! History: 
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 type ( particle_rams ),intent(in) :: & ! particle_rams is defined in module_simulater
     kext_rams,  & ! extinction coefficient [km-1]
     salb_rams,  & ! single scattering albedo [-]
     asym_rams     ! asymetry parameter [-]

 real(sdsu_fps),intent(out) ::  kexttot, salbtot, asymtot !total optical parameters
 real(sdsu_fps),intent(out) ::  kext_cloud !cloud component optical parameters

!
! total extinction 
!
   kexttot    =  &
                kext_rams%cloud1 + &
                kext_rams%cloud2 + &
                kext_rams%rain   + &
                kext_rams%ice1   + &
                kext_rams%ice2   + &
                kext_rams%snow   + &
                kext_rams%graupel+ &
                kext_rams%hail

   kext_cloud  = &
                kext_rams%cloud1 + &
                kext_rams%cloud2 + &
                kext_rams%ice1   + &
                kext_rams%ice2  

!
! total single scattering albedo
!
   if ( kexttot <= 0.e0 ) then
        salbtot = 0.e0
   else
       salbtot = ( &
                  salb_rams%cloud1 * kext_rams%cloud1  + &
                  salb_rams%cloud2 * kext_rams%cloud2  + &
                  salb_rams%rain   * kext_rams%rain    + &
                  salb_rams%ice1   * kext_rams%ice1    + &
                  salb_rams%ice2   * kext_rams%ice2    + &
                  salb_rams%snow   * kext_rams%snow    + &
                  salb_rams%graupel* kext_rams%graupel + &
                  salb_rams%hail   * kext_rams%hail      &
                 ) / kexttot

   endif

!
! total asymetry parameter
!
   if ( salbtot <= 0.e0 ) then
       asymtot = 0.e0
   else
       asymtot = ( &
                  asym_rams%cloud1 * salb_rams%cloud1 * kext_rams%cloud1  + &
                  asym_rams%cloud2 * salb_rams%cloud2 * kext_rams%cloud2  + &
                  asym_rams%rain   * salb_rams%rain   * kext_rams%rain    + &
                  asym_rams%ice1   * salb_rams%ice1   * kext_rams%ice1    + &
                  asym_rams%ice2   * salb_rams%ice2   * kext_rams%ice2    + &
                  asym_rams%snow   * salb_rams%snow   * kext_rams%snow    + &
                  asym_rams%graupel* salb_rams%graupel* kext_rams%graupel + &
                  asym_rams%hail   * salb_rams%hail   * kext_rams%hail      &
                 ) / ( salbtot * kexttot )

   endif

 return
 end subroutine total_opt_rams_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine  opt_gocart_broad( n, ib, q, rh, kext_gocart, salb_gocart, asym_gocart )
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation techniqu.
!  This routine is for GOCART microphysics.
!        
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; Revised completely for faster processing. 
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 integer,intent(in)::n  !band #
 integer,intent(in)::ib !radiation type
 type ( particle_gocart ),intent(in) :: & ! particle_gocart is defined in module_simulater
     q   ! mixing ratio [g/m3] 
 real(sdsu_fps),intent(in) :: rh ! relative humidity [0~1]  
 type ( particle_gocart ),intent(out) :: & ! particle_gocart is defined in module_simulater
     kext_gocart,  & ! extinction coefficient [km-1]
     salb_gocart,  & ! single scattering albedo [-]
     asym_gocart     ! asymetry parameter [-]


 integer :: j
 integer :: idx_spc
 real(sdsu_fps) :: wgt1, wgt2
 real(sdsu_fps) :: rh_in, q_in
 real(sdsu_fps) :: kext_q_unit , salb , asym

!
! initialize
!
 call init_aerosol( kext_gocart, 'zero')
 call init_aerosol( salb_gocart, 'zero')
 call init_aerosol( asym_gocart, 'zero')


!
! sulfate
!
 idx_spc = 1  ; q_in = q%so4 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
   if(rh_in < pts_rh(1)       ) rh_in = pts_rh(1)
   if(rh_in > pts_rh(mxpts_rh)) rh_in = pts_rh(mxpts_rh)
   !
   ! Linear interpolation
   !
   do j = 1,mxpts_rh-1
      if( rh_in >= pts_rh(j) .and. rh_in <= pts_rh(j+1) ) then
         wgt2 =  ( rh_in - pts_rh(j) )    / ( pts_rh(j+1) - pts_rh(j) )
         wgt1 = 1.0 - wgt2
         kext_q_unit =  wgt1*kext_q_unit_aero_array(j,idx_spc,n,ib) &
                      + wgt2*kext_q_unit_aero_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
         salb        =  wgt1*omega_aero_array      (j,idx_spc,n,ib) &
                      + wgt2*omega_aero_array      (j+1,idx_spc,n,ib) !single scattering albedo []
         asym        =  wgt1*asym_aero_array       (j,idx_spc,n,ib) &
                      + wgt2*asym_aero_array       (j+1,idx_spc,n,ib) !asymetry parameter []
         exit
      endif
   enddo
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%so4 = kext_q_unit * q_in ; salb_gocart%so4 = salb ; asym_gocart%so4 = asym  !outpout


!
! black carbon
!
 idx_spc = 2  ; q_in = q%blc ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit =  wgt1*kext_q_unit_aero_array(j,idx_spc,n,ib) &
                 + wgt2*kext_q_unit_aero_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        =  wgt1*omega_aero_array      (j,idx_spc,n,ib) &
                 + wgt2*omega_aero_array      (j+1,idx_spc,n,ib) !single scattering albedo []
    asym        =  wgt1*asym_aero_array       (j,idx_spc,n,ib) &
                 + wgt2*asym_aero_array       (j+1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%blc = kext_q_unit * q_in ; salb_gocart%blc = salb ; asym_gocart%blc = asym  !outpout

!
! non hygroscopic OC
!
 idx_spc = 3  ; q_in = q%ocn; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit =  kext_q_unit_aero_array(1,idx_spc,n,ib)
    salb        =  omega_aero_array      (1,idx_spc,n,ib)
    asym        =  asym_aero_array       (1,idx_spc,n,ib)
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%ocn = kext_q_unit * q_in ; salb_gocart%ocn = salb ; asym_gocart%ocn = asym  !outpout

!
! hygroscopic OC
!
 idx_spc = 4  ; q_in = q%och ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit =  wgt1*kext_q_unit_aero_array(j,idx_spc,n,ib) &
                 + wgt2*kext_q_unit_aero_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        =  wgt1*omega_aero_array      (j,idx_spc,n,ib) &
                 + wgt2*omega_aero_array      (j+1,idx_spc,n,ib) !single scattering albedo []
    asym        =  wgt1*asym_aero_array       (j,idx_spc,n,ib) &
                 + wgt2*asym_aero_array       (j+1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%och = kext_q_unit * q_in ; salb_gocart%och = salb ; asym_gocart%och = asym  !outpout

!
! fine-mode sea salt
!
 idx_spc = 5  ; q_in = q%ssa ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit =  wgt1*kext_q_unit_aero_array(j,idx_spc,n,ib) &
                 + wgt2*kext_q_unit_aero_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        =  wgt1*omega_aero_array      (j,idx_spc,n,ib) &
                 + wgt2*omega_aero_array      (j+1,idx_spc,n,ib) !single scattering albedo []
    asym        =  wgt1*asym_aero_array       (j,idx_spc,n,ib) &
                 + wgt2*asym_aero_array       (j+1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%ssa = kext_q_unit * q_in ; salb_gocart%ssa = salb ; asym_gocart%ssa = asym  !outpout

!
! coarse-mode sea salt
!
 idx_spc = 6  ; q_in = q%ssc ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit =  wgt1*kext_q_unit_aero_array(j,idx_spc,n,ib) &
                 + wgt2*kext_q_unit_aero_array(j+1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        =  wgt1*omega_aero_array      (j,idx_spc,n,ib) &
                 + wgt2*omega_aero_array      (j+1,idx_spc,n,ib) !single scattering albedo []
    asym        =  wgt1*asym_aero_array       (j,idx_spc,n,ib) &
                 + wgt2*asym_aero_array       (j+1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%ssc = kext_q_unit * q_in ; salb_gocart%ssc = salb ; asym_gocart%ssc = asym  !outpout


!
! dust 1
!
 idx_spc = 7  ; q_in = q%du1 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%du1 = kext_q_unit * q_in ; salb_gocart%du1 = salb ; asym_gocart%du1 = asym  !outpout

!
! dust2 
!
 idx_spc = 8  ; q_in = q%du2 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%du2 = kext_q_unit * q_in ; salb_gocart%du2 = salb ; asym_gocart%du2 = asym  !outpout

!
! dust 3
!
 idx_spc = 9  ; q_in = q%du3 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%du3 = kext_q_unit * q_in ; salb_gocart%du3 = salb ; asym_gocart%du3 = asym  !outpout

!
! dust 4
!
 idx_spc = 10  ; q_in = q%du4 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%du4 = kext_q_unit * q_in ; salb_gocart%du4 = salb ; asym_gocart%du4 = asym  !outpout

!
! dust 5 
!
 idx_spc = 11  ; q_in = q%du5 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%du5 = kext_q_unit * q_in ; salb_gocart%du5 = salb ; asym_gocart%du5 = asym  !outpout

!
! dust 6
!
 idx_spc = 12  ; q_in = q%du6 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%du6 = kext_q_unit * q_in ; salb_gocart%du6 = salb ; asym_gocart%du6 = asym  !outpout

!
! dust 7
!
 idx_spc = 13  ; q_in = q%du7 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%du7 = kext_q_unit * q_in ; salb_gocart%du7 = salb ; asym_gocart%du7 = asym  !outpout

!
! dust 8
!
 idx_spc = 14  ; q_in = q%du8 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_gocart%du8 = kext_q_unit * q_in ; salb_gocart%du8 = salb ; asym_gocart%du8 = asym  !outpout


 return
 end subroutine opt_gocart_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine total_opt_gocart_broad( kext_gocart, salb_gocart, asym_gocart, &
                                 kexttot, salbtot, asymtot )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute total single-scattering properties from atmospheric particles.
!  This routine is for GCE microphysics.
!  
! History: 
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 type ( particle_gocart ),intent(in) :: & ! particle_gocart is defined in module_simulater
     kext_gocart,  & ! extinction coefficient [km-1]
     salb_gocart,  & ! single scattering albedo [-]
     asym_gocart     ! asymetry parameter [-]

 real(sdsu_fps),intent(out) ::  kexttot, salbtot, asymtot !total optical parameters

!
! total extinction 
!
   kexttot    = &
                kext_gocart%so4  + &
                kext_gocart%blc  + &
                kext_gocart%ocn  + &
                kext_gocart%och  + &
                kext_gocart%ssa  + &
                kext_gocart%ssc  + &
                kext_gocart%du1  + &
                kext_gocart%du2  + &
                kext_gocart%du3  + &
                kext_gocart%du4  + &
                kext_gocart%du5  + &
                kext_gocart%du6  + &
                kext_gocart%du7  + &
                kext_gocart%du8  

!
! total single scattering albedo
!
   if ( kexttot <= 0.e0 ) then
        salbtot = 0.e0
   else
       salbtot = ( &
                  salb_gocart%so4 * kext_gocart%so4   + &
                  salb_gocart%blc * kext_gocart%blc   + &
                  salb_gocart%ocn * kext_gocart%ocn   + &
                  salb_gocart%och * kext_gocart%och   + &
                  salb_gocart%ssa * kext_gocart%ssa   + &
                  salb_gocart%ssc * kext_gocart%ssc   + &
                  salb_gocart%du1 * kext_gocart%du1   + &
                  salb_gocart%du2 * kext_gocart%du2   + &
                  salb_gocart%du3 * kext_gocart%du3   + &
                  salb_gocart%du4 * kext_gocart%du4   + &
                  salb_gocart%du5 * kext_gocart%du5   + &
                  salb_gocart%du6 * kext_gocart%du6   + &
                  salb_gocart%du7 * kext_gocart%du7   + &
                  salb_gocart%du8 * kext_gocart%du8     &
                 ) / kexttot

   endif


!
! total asymetry parameter
!
   if ( salbtot <= 0.e0 ) then
       asymtot = 0.e0
   else
       asymtot = ( &
                  asym_gocart%so4  * salb_gocart%so4  * kext_gocart%so4   + &
                  asym_gocart%blc  * salb_gocart%blc  * kext_gocart%blc   + &
                  asym_gocart%ocn  * salb_gocart%ocn  * kext_gocart%ocn   + &
                  asym_gocart%och  * salb_gocart%och  * kext_gocart%och   + &
                  asym_gocart%ssa  * salb_gocart%ssa  * kext_gocart%ssa   + &
                  asym_gocart%ssc  * salb_gocart%ssc  * kext_gocart%ssc   + &
                  asym_gocart%du1  * salb_gocart%du1  * kext_gocart%du1   + &
                  asym_gocart%du2  * salb_gocart%du2  * kext_gocart%du2   + &
                  asym_gocart%du3  * salb_gocart%du3  * kext_gocart%du3   + &
                  asym_gocart%du4  * salb_gocart%du4  * kext_gocart%du4   + &
                  asym_gocart%du5  * salb_gocart%du5  * kext_gocart%du5   + &
                  asym_gocart%du6  * salb_gocart%du6  * kext_gocart%du6   + &
                  asym_gocart%du7  * salb_gocart%du7  * kext_gocart%du7   + &
                  asym_gocart%du8  * salb_gocart%du8  * kext_gocart%du8     &
                 ) / ( salbtot * kexttot )

   endif

 return
 end subroutine total_opt_gocart_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine  opt_ucd_broad( n, ib, q, rh, kext_ucd, salb_ucd, asym_ucd )
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation techniqu.
!  This routine is for GOCART microphysics.
!        
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 integer,intent(in)::n  !band #
 integer,intent(in)::ib !radiation type
 type ( particle_ucd ),intent(in) :: & ! particle_ucd is defined in module_simulater
     q   ! mixing ratio [g/m3] 
 real(sdsu_fps),intent(in) :: rh ! relative humidity [0~1]  
 type ( particle_ucd ),intent(out) :: & ! particle_ucd is defined in module_simulater
     kext_ucd,  & ! extinction coefficient [km-1]
     salb_ucd,  & ! single scattering albedo [-]
     asym_ucd     ! asymetry parameter [-]


 integer :: j
 integer :: idx_spc
 real(sdsu_fps) :: wgt1, wgt2
 real(sdsu_fps) :: rh_in, q_in
 real(sdsu_fps) :: kext_q_unit , salb , asym

!
! initialize
!
 call init_aerosol( kext_ucd, 'zero')
 call init_aerosol( salb_ucd, 'zero')
 call init_aerosol( asym_ucd, 'zero')


!
! dust 1
!
 idx_spc = 1  ; q_in = q%du1 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_ucd%du1 = kext_q_unit * q_in ; salb_ucd%du1 = salb ; asym_ucd%du1 = asym  !outpout

!
! dust2 
!
 idx_spc = 2  ; q_in = q%du2 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_ucd%du2 = kext_q_unit * q_in ; salb_ucd%du2 = salb ; asym_ucd%du2 = asym  !outpout

!
! dust 3
!
 idx_spc = 3  ; q_in = q%du3 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_ucd%du3 = kext_q_unit * q_in ; salb_ucd%du3 = salb ; asym_ucd%du3 = asym  !outpout

!
! dust 4
!
 idx_spc = 4  ; q_in = q%du4 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_ucd%du4 = kext_q_unit * q_in ; salb_ucd%du4 = salb ; asym_ucd%du4 = asym  !outpout

!
! dust 5 
!
 idx_spc = 5  ; q_in = q%du5 ; rh_in = rh
 if(q_in <= q_min_aerosol) then
    kext_q_unit = 0. ; salb = 0. ; asym = 0.
 else
    kext_q_unit = kext_q_unit_aero_array(1,idx_spc,n,ib) !extincion per unit q [km-1/(g/m3)]
    salb        = omega_aero_array      (1,idx_spc,n,ib) !single scattering albedo []
    asym        = asym_aero_array       (1,idx_spc,n,ib) !asymetry parameter []
 endif
! Re-scale kext from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 kext_ucd%du5 = kext_q_unit * q_in ; salb_ucd%du5 = salb ; asym_ucd%du5 = asym  !outpout


 return
 end subroutine  opt_ucd_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

subroutine total_opt_ucd_broad( kext_ucd, salb_ucd, asym_ucd, &
                                 kexttot, salbtot, asymtot )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute total single-scattering properties from atmospheric particles.
!  This routine is for GCE microphysics.
!  
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 type ( particle_ucd ),intent(in) :: & ! particle_ucd is defined in module_simulater
     kext_ucd,  & ! extinction coefficient [km-1]
     salb_ucd,  & ! single scattering albedo [-]
     asym_ucd     ! asymetry parameter [-]

 real(sdsu_fps),intent(out) ::  kexttot, salbtot, asymtot !total optical parameters

!
! total extinction 
!
   kexttot    = &
                kext_ucd%du1  + &
                kext_ucd%du2  + &
                kext_ucd%du3  + &
                kext_ucd%du4  + &
                kext_ucd%du5  

!
! total single scattering albedo
!
   if ( kexttot <= 0.e0 ) then
        salbtot = 0.e0
   else
       salbtot = ( &
                  salb_ucd%du1 * kext_ucd%du1   + &
                  salb_ucd%du2 * kext_ucd%du2   + &
                  salb_ucd%du3 * kext_ucd%du3   + &
                  salb_ucd%du4 * kext_ucd%du4   + &
                  salb_ucd%du5 * kext_ucd%du5    &
                 ) / kexttot

   endif


!
! total asymetry parameter
!
   if ( salbtot <= 0.e0 ) then
       asymtot = 0.e0
   else
       asymtot = ( &
                  asym_ucd%du1  * salb_ucd%du1  * kext_ucd%du1   + &
                  asym_ucd%du2  * salb_ucd%du2  * kext_ucd%du2   + &
                  asym_ucd%du3  * salb_ucd%du3  * kext_ucd%du3   + &
                  asym_ucd%du4  * salb_ucd%du4  * kext_ucd%du4   + &
                  asym_ucd%du5  * salb_ucd%du5  * kext_ucd%du5    &
                 ) / ( salbtot * kexttot )

   endif

 return
 end subroutine total_opt_ucd_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine makelut_broad
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Construct the look-up tables of single-scattering properties (SSP) for particles.
! For each broadband, a number of monocromatic SSP is summe with weight to compute broad-band
! SSP. So, if you want to derive precise value, numband_lut must be large (100), which initially
! takes for a while to compute SSP. But once LUTs are made, LUTs needs to be read for next time.
! Mie routine used here is ones in the VisibleIR simulator for consistency. 
! 
! History: 
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 character(len=30) :: fname_out(20)
 character(len=100):: char_wav
 integer :: max_spc ! maximum particle species
 integer :: ispc    !RAMS particle habit index
 real(sdsu_fps) :: pbck     ! backscatter phase function/(4*pi)
 real(sdsu_fps) :: frequency
 integer :: nspec,i,j,nwav,ib,isl !loop indice
 integer :: iounit
 logical :: fexist
 integer :: mxwav ! maximum # of wavelength channels 
 real(sdsu_fps) :: cfmas, pwmas ! alpha_m & beta_m 
 integer :: ns,nt !humidity and temperature index
 real(sdsu_fps) :: ws,we

 character(len=2):: rad_typ !SW or LW
 character(len=2):: char_ib !band # in character
 integer :: max_band  !max # of band

 real(sdsu_fps) :: subw( numband_lut ) !sub band for making LUT
 real(sdsu_fps) :: rhoe  !snow density [kg/m3]

 integer,save:: nang                 !! number of scaattering angles
 real(sdsu_fpd),save:: ang(knang)           !! scattering angles [degree]
 real(sdsu_fpd):: scr                  !! relative water reflactive index (real part)
 real(sdsu_fpd):: sci                  !! relative water reflactive index (imaginary part)

 real(sdsu_fpd),parameter :: q_unit = 1.0d0  !unit mass [g/m3]
 real(sdsu_fpd) :: &
     wlc      ,  &  !wavelength [cm]
     rhl      ,  &  ! RH [0~1]
     cextp_lut,  &  ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_lut,  &  ! absorption cross section [cm2/cm3] or [cm-1]
     php_lut(knang),&! volume scattering phase function [-]
     cextp_sum,  &  ! weigting sum of extinction cross section [cm2/cm3] or [cm-1]
     cabsp_sum,  &  ! weigting sum of absorption cross section [cm2/cm3] or [cm-1]
     php_sum(knang),&! weigting sum of volume scattering phase function [-]
     kext_q_unit, & ! extinction cross section per unit mass [km-1 / (g/m3) ]
     omega     , &  ! signle scattering albedo [-]
     gg(klgn1)  , & !
     phase_func(knang),& ! phase function [-]
     asym           ! asymetry parameter (forward scattering fraction) [-] 

 integer :: mxlgn1
 real(sdsu_fpd):: cs(knang),yy(knang)
 real(sdsu_fpd),parameter :: unit_depth = 1.d0 ! layer depth [cm]
 real(sdsu_fpd):: wgt,wgt_sum !weight for solar constant 
 real(sdsu_fpd):: fsol,fplnk
 real(sdsu_fpd):: sunir                !! function of solar irradiance
 real(sdsu_fpd):: bplnk
 integer :: mxisl

! for icescat data
 integer,parameter :: iuk = 120
 integer :: intvl, ipol
  real(8):: szp(kintvl)          !! size parameter
 logical :: logic_opened


 integer :: nband
 integer :: ierr

!
! NetCDF parameters
!
 integer :: ncid, id
 integer :: ref_dimid, spc_dimid, bnd_dimid, typ_dimid
 integer, parameter :: ndims4d = 4
 integer :: dimids4d(ndims4d)
 integer :: id_kext, id_omega, id_asym


!
! 
!
 if(masterproc) print*,'MSG makelut_broad: LUTs for broadband simulation - '

!
! static indice
!
   mxisl = 2  !1-sw, 2-lw

   select case(trim(type_microphysics))
     case('GEN')    ! GCE (LIN_type) 1moment bulk microphsics 
      max_spc = 6
     case('RAMS')  !GCE RAMS 1- and 2-moment microphysics 
      max_spc = 16
     case('SBM')
        call stop_sdsu(myrank,'Not yet supporoting SBM')
     case default
   end select

   nband = max( nband_sw, nband_lw )

!
! allocate array parameters
!
 if( .not. allocated(kext_q_unit_array) ) then

   allocate( kext_q_unit_array ( mxpts_re_visir, max_spc, nband, mxisl ), &
             omega_array       ( mxpts_re_visir, max_spc, nband, mxisl ), &
             asym_array        ( mxpts_re_visir, max_spc, nband, mxisl ), &
             stat=ierr )
   if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_broad: allocation error -> Terminate program.')

 endif

!
! Initialization of fundamental particle refractive index (rfi)
! outputs are stored in the saved global parameter in module paras
!
  call read_rfi

   !first derive phase angle for consistent interpolartion for icescat database.
    call open_oldfile(iuk, trim(sdsu_dir_data)//'KRNL.OUT' )
    rewind iuk
    read(iuk,*) intvl, nang, ipol  ! 73, 74, 1  
    read(iuk,*) szp(1:intvl)  ! size parameter (0.1371 ~ 938.0: Mie regime)
    read(iuk,*) ang(1:nang)   ! scattering angle [deg] (0-180deg) 0: forward 180: backscatter
    close(iuk)

!
! open NetCDF file output 
!

!
! open netcdf file
!
 sdsu_io_file = trim(sdsu_dir_sslut)//'BROADBAND_CLOUD_'//trim(cloud_microphysics)//'.nc'

 inquire(file=trim(sdsu_io_file),exist=fexist)
 if( lut_replace ) fexist = .false.
 if (fexist) then !skip
     print*,'MSG makelut_broad; ',trim(sdsu_io_file),' exist.'
     return
 endif 

 print*,'MSG makelut_broad; Generating Broad-band LUT ',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )               !32bit output

!
! Define netcdf dimensions
!
 call check( nf90_def_dim(ncid, 're_max'   , mxpts_re_visir   , ref_dimid) )
 call check( nf90_def_dim(ncid, 'spc_max'  , max_spc          , spc_dimid) )
 call check( nf90_def_dim(ncid, 'band_max' , nband            , bnd_dimid) )
 call check( nf90_def_dim(ncid, 'typ_max'  , mxisl            , typ_dimid) )

 dimids4d  = (/ ref_dimid, spc_dimid, bnd_dimid, typ_dimid /)


!
! output parameter definition
!
 call check( nf90_def_var(ncid, 'kext_q_unit', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' extinction cross section per unit mass ') )
 call check( nf90_put_att(ncid, id, 'units', '[km-1 / (g/m3)]') )
 id_kext = id

 call check( nf90_def_var(ncid, 'omega', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' single scattering albedo ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_omega = id

 call check( nf90_def_var(ncid, 'asym', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' asymetry parameter ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_asym = id

 call check( nf90_enddef(ncid) )  ! End define mode.


!
! loop 
!

  SWLW_LOOP: do isl = 1,mxisl  

    if(isl == 1) then !SW
       max_band = nband_sw
       rad_typ = 'SW'
    else              !LW
       max_band = nband_lw
       rad_typ = 'LW'
    endif

  BAND_LOOP: do ib = 1, max_band

    if(isl == 1) then !SW
       ws = wavel_min(ib)  !starting wavelength for the ib-band [micron]
       we = wavel_max(ib)  !ending wavelength for the ib-band [micron]
       call compute_subw( ws,we,subw )
    else              !LW
       ws = waven_min(ib)  !starting wavenumber for the ib-band [cm-1]
       we = waven_max(ib)  !ending wavenumber for the ib-band [cm-1]
       call compute_subw( ws,we,subw )
    endif
 
   write(char_ib,"(I2.2)") ib  !from integer to character 

   if(masterproc) print*,''
   mic_select1: select case(trim(type_microphysics))
     case('GEN')    ! GCE (LIN_type) 1moment bulk microphsics 

      fname_out(1) = 'CW-'//trim(cloud_microphysics)//'-'//rad_typ//'band'//trim(char_ib)
      fname_out(2) = 'RR-'//trim(cloud_microphysics)//'-'//rad_typ//'band'//trim(char_ib)
      fname_out(3) = 'CI-'//trim(cloud_microphysics)//'-'//rad_typ//'band'//trim(char_ib)
      fname_out(4) = 'SN-'//trim(cloud_microphysics)//'-'//rad_typ//'band'//trim(char_ib)
      fname_out(5) = 'GR-'//trim(cloud_microphysics)//'-'//rad_typ//'band'//trim(char_ib)
      fname_out(6) = 'HL-'//trim(cloud_microphysics)//'-'//rad_typ//'band'//trim(char_ib)

     case('RAMS')  !GCE RAMS 1- and 2-moment microphysics 

      fname_out(1)  = 'C1-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(2)  = 'C2-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(3)  = 'RR-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(4)  = 'I1-COL-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(5)  = 'I1-HEX-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(6)  = 'I1-DEN-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(7)  = 'I1-NDL-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(8)  = 'I1-ROS-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(9)  = 'I2-COL-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(10) = 'I2-HEX-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(11) = 'I2-DEN-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(12) = 'I2-NDL-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(13) = 'I2-ROS-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(14) = 'SN-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(15) = 'GR-RAMS-'//rad_typ//'band'//trim(char_ib)
      fname_out(16) = 'HL-RAMS-'//rad_typ//'band'//trim(char_ib)

     case('SBM')
 
        call stop_sdsu(myrank,'Not yet supporoting SBM')

     case default 
        ! do nothing
   end select mic_select1

!
! species loop 
!
 SPC_LOOP: do nspec = 1, max_spc


  print*,'MSG makelut_broad; Processing ',trim(fname_out(nspec))


 RE_LOOP: do i = 1, mxpts_re_visir
    cextp_sum=0.d0 ;cabsp_sum=0.d0 ;php_sum=0.d0 ; wgt_sum=0.d0

 WAV_LOOP : do nwav = 1 , numband_lut

!
! sub wavelength/wavenumber weight 
!
      if(isl == 1) then !SW subw is [micron]
         fsol=sunir( 1.d4/subw(nwav) )  !TOA spectrum solar constant
         wgt = fsol
         wlc = subw(nwav)*1.0d-4  ![cm]  <- [micron]
      else              !LW  subw is [cm-1]
         fplnk      =  bplnk( 1.d4/subw(nwav) , DBLE(earth_mean_temp) )  !plank function
         wgt = fplnk
         wlc = 1.d0 /subw(nwav)  ![cm] <- [cm-1]
      endif


!
! Compute single-scattering parameters 
!
       cextp_lut=0.d0 ;cabsp_lut=0.d0 ;php_lut=0.d0 
       rhl = 0.1d0 ! RH (0~1) not used for condensates (but aerosols)

       mic_select2: select case(trim(type_microphysics))
       case('GEN')      ! General 1moment bulk microphsics 

           NSPEC_GCE: if(nspec==1) then  !cloud use gamma PSD in RAMS (cloud1 class)

!             if( trim(cloud_microphysics) == 'THOM') then
!                call mie_gen_visir('qc ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%cloud), &
!                                   dble(rho_gen%cloud), wlc,rhl, & !cloud
!                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )
!             else

              ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud
              call mie_rams_visir('qc1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%cloud1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )
!             endif

            elseif(nspec==2) then
                call mie_gen_visir('qr ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%rain), &
                                  dble(rho_gen%rain), wlc,rhl, & !rain
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )
            elseif(nspec==3) then !ice use gamma PSD in RAMS (ice column)

              if (trim(cloud_microphysics)=='MORR' .or. trim(cloud_microphysics)=='MORRH' .or. &
                  trim(cloud_microphysics)=='THOM' ) then
                call mie_gen_visir('qi ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%ice), &
                                  dble(rho_gen%ice), wlc,rhl, & !ice
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )

              else

                 ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice
                 call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )
              endif

            elseif(nspec==4) then


              if( trim(cloud_microphysics) == 'GMP4ICE' ) then
               rhoe = 1.e+3 * max(0.05, min(0.9e0, 142.36 * ( pts_re_visir(i)**(-0.95) ) ) ) ! re-dependent bulk snow density [kg/m3]
              elseif( trim(cloud_microphysics) == 'THOM' ) then
               rhoe = min(900. , 0.13 / (pts_re(j)*2.*1.e-6) ) ! Kg/m^3
              else
                rhoe = rho_gen%snow
              endif


                call mie_gen_visir('qs ', q_unit, dble(pts_re_visir(i)),dble(mu_gen%snow), &
                                   dble(rhoe), wlc,rhl, & !snow
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )

            elseif(nspec==5) then
                call mie_gen_visir('qg ', q_unit, dble(pts_re_visir(i)),dble(mu_gen%graupel), &
                                    dble(rho_gen%graupel), wlc,rhl, & !graupel
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )

            elseif(nspec==6) then
                call mie_gen_visir('qh ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%hail), &
                                   dble(rho_gen%hail), wlc,rhl, & !hail
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )
            else
               call stop_sdsu(myrank, 'MSG makelut_visir: there is no such nspec')
            endif NSPEC_GCE


       case('RAMS')  !GCE RAMS 1- & 2-moment microphysics 

            NSPEC_RAMS: if(nspec==1) then
              ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud1
              call mie_rams_visir('qc1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%cloud1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==2) then
              ispc = 16; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud2
              call mie_rams_visir('qc2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%cloud2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==3) then
              ispc = 2 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! rain 
              call mie_rams_visir('qr ',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%rain), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==4) then
              ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 col
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==5) then
              ispc = 8 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 hex
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==6) then
              ispc = 9 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 den
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==7) then
              ispc = 10 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 ndl
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==8) then
              ispc = 11 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 ros
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==9) then
              ispc = 4 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 col
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==10) then
              ispc = 12 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 hex
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==11) then
              ispc = 13 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 den
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==12) then
              ispc = 14 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 ndl
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==13) then
              ispc = 15 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 ros
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==14) then
              ispc = 5 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! snow aggregate
              call mie_rams_visir('qs ',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%snow), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==15) then
              ispc = 6 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! graupel
              call mie_rams_visir('qg ',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%graupel), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==16) then
               ispc = 7 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! hail
              call mie_rams_visir('qh ',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%hail), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            else
               call stop_sdsu(myrank,'MSG makelut_visir: there is no such nspec')
            endif NSPEC_RAMS

       case('SBM') 

          call stop_sdsu(myrank,'Not yet supporting SBM')

       case default 

          ! do nothing

       end select mic_select2

!
! summing sub-band single-scattering component
!

     !
     ! avoid NaN case
     !
     if( isnan(cextp_lut) .or. cextp_lut == undefined .or. cextp_lut == 0.0 ) then 
        cextp_lut = 0.0
        cabsp_lut = 0.0
          php_lut = 0.0
     endif

     cextp_sum = cextp_sum + cextp_lut * wgt ! weigting sum of extinction cross section [cm2/cm3] or [cm-1]
     cabsp_sum = cabsp_sum + cabsp_lut * wgt ! weigting sum of absorption cross section [cm2/cm3] or [cm-1]
     php_sum   =   php_sum +   php_lut * wgt ! weigting sum of volume scattering phase function [-]
     wgt_sum   =   wgt_sum +             wgt ! weightng sum

 enddo WAV_LOOP

!
! derive band-mean values of single-scattering parameters
! Note these values are per unit mixing ratio
!
     cextp_lut = cextp_sum / wgt_sum ! broad-band extinction cross section [cm2/cm3] or [cm-1]
     cabsp_lut = cabsp_sum / wgt_sum ! broad-band absorption cross section [cm2/cm3] or [cm-1]
       php_lut =   php_sum / wgt_sum ! broad-band volume scattering phase function [-]


     kext_q_unit = cextp_lut * 1.d+5      ! [km-1] <- [cm-1]

     if( cextp_lut == 0. ) then
       omega = 0.
       asym  = 0.
     else
       omega       = (cextp_lut-cabsp_lut) / cextp_lut  ! signle scattering albedo [-]


!
! Compute aymetry parameter here
!
       phase_func  =  php_lut/(cextp_lut-cabsp_lut) !  
       cs(1:nang) = cos(ang(1:nang)*rad)
       yy(1:nang) = phase_func(1:nang)*unit_depth ! phsf(1:nang,l)
       mxlgn1 = N_stream*2 + 1  !same as nda 
       call lgndf3(mxlgn1,nang,cs,yy,gg)
       asym = gg(2)/gg(1)  ! asymetry paramer

     endif

     !
     ! avoid NaN case
     !
     if( isnan(kext_q_unit) .or. kext_q_unit == 0.0 ) then
        omega = 0.0
        asym  = 0.0
     endif

     !write(iounit,90) kext_q_unit, omega, asym

     !
     ! large array output
     !
     kext_q_unit_array(i, nspec, ib, isl) = real(kext_q_unit)
     omega_array      (i, nspec, ib, isl) = real(omega)
     asym_array       (i, nspec, ib, isl) = real(asym)


 enddo RE_LOOP

! close(iounit) !close LUT table

!   endif INQ


 enddo SPC_LOOP

 enddo BAND_LOOP


 enddo SWLW_LOOP


!
! write output in NetCDF
!
 call check( nf90_put_var(ncid, id_kext  , kext_q_unit_array ) )
 call check( nf90_put_var(ncid, id_omega , omega_array       ) )
 call check( nf90_put_var(ncid, id_asym  , asym_array        ) )

 call check( nf90_close(ncid) )  ! close nc file


 90  Format ( 1P3E14.6 )

 return
 end subroutine makelut_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine makelut_broad_aerosol
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Construct the look-up tables of single-scattering properties (SSP) for particles.
! For each broadband, a number of monocromatic SSP is summe with weight to compute broad-band
! SSP. So, if you want to derive precise value, numband_lut must be large (100), which initially
! takes for a while to compute SSP. But once LUTs are made, LUTs needs to be read for next time.
! Mie routine used here is ones in the VisibleIR simulator for consistency. 
!
! History: 
!  08/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 character(len=30) :: fname_out(20)
 character(len=100):: char_wav
 integer :: max_spc ! maximum particle species
 integer :: ispc    !RAMS particle habit index
 real(sdsu_fps) :: pbck     ! backscatter phase function/(4*pi)
 real(sdsu_fps) :: frequency
 integer :: nspec,i,j,nwav,ib,isl !loop indice
 integer :: iounit
 logical :: fexist 
 integer :: mxwav ! maximum # of wavelength channels 
 real(sdsu_fps) :: cfmas, pwmas ! alpha_m & beta_m 
 integer :: ns,nt !humidity and temperature index
 real(sdsu_fps) :: ws,we

 character(len=2):: rad_typ !SW or LW
 character(len=2):: char_ib !band # in character
 integer :: max_band  !max # of band

 real(sdsu_fps) :: subw( numband_lut ) ! sub band for making LUT

 integer,save:: nang                   ! number of scaattering angles
 real(sdsu_fpd),save:: ang(knang)      ! scattering angles [degree]
 real(sdsu_fpd):: scr                  ! relative water reflactive index (real part)
 real(sdsu_fpd):: sci                  ! relative water reflactive index (imaginary part)

 real(sdsu_fpd),parameter :: q_unit = 1.0d0   !unit mass [g/m3]
 real(sdsu_fpd) :: &
     wlc      ,  &  !wavelength [cm]
     rhl      ,  &  ! RH [0~1]
     cextp_lut,  &  ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_lut,  &  ! absorption cross section [cm2/cm3] or [cm-1]
     php_lut(knang),&! volume scattering phase function [-]
     cextp_sum,  &  ! weigting sum of extinction cross section [cm2/cm3] or [cm-1]
     cabsp_sum,  &  ! weigting sum of absorption cross section [cm2/cm3] or [cm-1]
     php_sum(knang),&! weigting sum of volume scattering phase function [-]
     kext_q_unit, & ! extinction cross section per unit mass [km-1 / (g/m3) ]
     omega     , &  ! signle scattering albedo [-]
     gg(klgn1)  , & !
     phase_func(knang),& ! phase function [-]
     asym       ,&    ! asymetry parameter (forward scattering fraction) [-] 
     kext_q_unit_previous, & ! extinction cross section per unit mass [km-1 / (g/m3) ]
     omega_previous      , &  ! signle scattering albedo [-]
     asym_previous           ! asymetry parameter (forward scattering fraction) [-] 


 integer :: mxlgn1
 real(sdsu_fpd):: cs(knang),yy(knang)
 real(sdsu_fpd),parameter :: unit_depth = 1.d0 ! layer depth [cm]
 real(sdsu_fpd):: wgt,wgt_sum !weight for solar constant 
 real(sdsu_fpd):: fsol,fplnk
 real(sdsu_fpd):: sunir                !! function of solar irradiance
 real(sdsu_fpd):: bplnk
 integer :: mxisl


 integer :: nband
 integer :: ierr

!
! NetCDF parameters
!
 integer :: ncid, id
 integer :: rh_dimid, spc_dimid, bnd_dimid, typ_dimid
 integer, parameter :: ndims4d = 4
 integer :: dimids4d(ndims4d)
 integer :: id_kext, id_omega, id_asym


!
!
!
 if(masterproc) print*,'MSG makelut_broad_aerosol: LUTs for broadband simulation - '

 mxisl = 2  !1-sw, 2-lw

 select case(trim(aerosol_microphysics) )
 case('GOCART')
   max_spc = 14
 case('UCD')
   max_spc = 5
 end select

 nband = max( nband_sw, nband_lw )

!
! allocate array parameters
!
 if( .not. allocated(kext_q_unit_aero_array) ) then

   allocate( kext_q_unit_aero_array ( mxpts_rh, max_spc, nband, mxisl ), &
             omega_aero_array       ( mxpts_rh, max_spc, nband, mxisl ), &
             asym_aero_array        ( mxpts_rh, max_spc, nband, mxisl ), &
             stat=ierr )
   if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_broad_aerosol: allocation error -> Terminate program.')

 endif


!
! initialization of fundamental particle refractive index (rfi)
! outputs are stored in the saved global parameter in module paras
!
  call read_rfi


!
! open netcdf file
!
 sdsu_io_file = trim(sdsu_dir_sslut)//'BROADBAND_AEROSOL_'//trim(aerosol_microphysics)//'.nc'

 inquire(file=trim(sdsu_io_file),exist=fexist)
 if( lut_replace ) fexist = .false.
 if (fexist) then !skip
     print*,'MSG makelut_broad_aerosol; ',trim(sdsu_io_file),' exist.'
     return
 endif

 print*,'MSG makelut_broad_aerosol; Generating Broad-band LUT ',trim(fname_out(nspec))

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )               !32bit output

!
! Define netcdf dimensions
!
 call check( nf90_def_dim(ncid, 'rh_max'   , mxpts_rh         , rh_dimid) )
 call check( nf90_def_dim(ncid, 'spc_max'  , max_spc          , spc_dimid) )
 call check( nf90_def_dim(ncid, 'band_max' , nband            , bnd_dimid) )
 call check( nf90_def_dim(ncid, 'typ_max'  , mxisl            , typ_dimid) )

 dimids4d  = (/ rh_dimid, spc_dimid, bnd_dimid, typ_dimid /)

!
! output parameter definition
!
 call check( nf90_def_var(ncid, 'kext_q_unit', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' extinction cross section per unit mass ') )
 call check( nf90_put_att(ncid, id, 'units', '[km-1 / (g/m3)]') )
 id_kext = id

 call check( nf90_def_var(ncid, 'omega', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' single scattering albedo ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_omega = id

 call check( nf90_def_var(ncid, 'asym', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' asymetry parameter ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_asym = id

 call check( nf90_enddef(ncid) )  ! End define mode.



!
! loop
!

  SWLW_LOOP: do isl = 1,mxisl  

   if(isl == 1) then !SW
      max_band = nband_sw
      rad_typ = 'SW'
   else              !LW
      max_band = nband_lw
      rad_typ = 'LW'
   endif

  BAND_LOOP: do ib = 1, max_band
 
   if(isl == 1) then !SW
      ws = wavel_min(ib)  !starting wavelength for the ib-band [micron]
      we = wavel_max(ib)  !ending wavelength for the ib-band [micron]
      call compute_subw( ws,we,subw )
   else              !LW
      ws = waven_min(ib)  !starting wavenumber for the ib-band [cm-1]
      we = waven_max(ib)  !ending wavenumber for the ib-band [cm-1]
      call compute_subw( ws,we,subw )
   endif

   write(char_ib,"(I2.2)") ib  !from integer to character 

   if(masterproc) print*,''

!
! LUT names
!
    select case(trim(aerosol_microphysics) )
    case('GOCART')

      fname_out(1)  = 'SO4-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(2)  = 'BLC-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(3)  = 'OCN-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(4)  = 'OCH-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(5)  = 'SSA-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(6)  = 'SSC-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(7)  = 'DU1-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(8)  = 'DU2-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(9)  = 'DU3-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(10) = 'DU4-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(11) = 'DU5-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(12) = 'DU6-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(13) = 'DU7-GOCART-'//rad_typ//'band'//trim(char_ib)
      fname_out(14) = 'DU8-GOCART-'//rad_typ//'band'//trim(char_ib)

    case('UCD')

      fname_out(1) = 'DU1-UCD-'//rad_typ//'band'//trim(char_ib)
      fname_out(2) = 'DU2-UCD-'//rad_typ//'band'//trim(char_ib)
      fname_out(3) = 'DU3-UCD-'//rad_typ//'band'//trim(char_ib)
      fname_out(4) = 'DU4-UCD-'//rad_typ//'band'//trim(char_ib)
      fname_out(5) = 'DU5-UCD-'//rad_typ//'band'//trim(char_ib)

    end select


!
! species loop 
!
 SPC_LOOP: do nspec = 1,  max_spc

   print*,'MSG makelut_broad_aerosol; Processing ',trim(fname_out(nspec))

!
!       Calculate scattering parameters per unit mass mixing ratio 1 [g/cm3]
!
  RH_LOOP: do i = 1, mxpts_rh
    cextp_sum=0.d0 ;cabsp_sum=0.d0 ;php_sum=0.d0 ; wgt_sum=0.d0


  WAV_LOOP : do nwav = 1 , numband_lut

!
! sub wavelength/wavenumber weight 
!
    if(isl == 1) then !SW subw is [micron]
       fsol=sunir( 1.d4/subw(nwav) )  !TOA spectrum solar constant
       wgt = fsol
       wlc = subw(nwav)*1.0d-4  ![cm]  <- [micron]
    else              !LW  subw is [cm-1]
       fplnk      =  bplnk( 1.d4/subw(nwav) , DBLE(earth_mean_temp) )  !plank function
       wgt = fplnk
       wlc = 1.d0 /subw(nwav)  ![cm] <- [cm-1]
    endif


!
! Compute single-scattering parameters per unit mass (= 1[g/m3] = 1.d-6[g/cm3] )
!
       cextp_lut=0.d0 ;cabsp_lut=0.d0 ;php_lut=0.d0
       rhl = DBLE( pts_rh(i) )   ! RH (0~1) for condensates (aerosols)

         select case(trim(aerosol_microphysics) )
         case('GOCART')

            NSPEC_GOCART: if(nspec==1) then
                call mie_gocart('so4', q_unit*1.0d-6, wlc,rhl,         &   !sulfur and its precurser
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )

            elseif(nspec==2) then
                call mie_gocart('blc', q_unit*1.0d-6, wlc,rhl,         &   !black carbon (soot)
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==3) then
                call mie_gocart('ocn', q_unit*1.0d-6, wlc,rhl,         &   !non hygroscopic OC
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==4) then
                call mie_gocart('och', q_unit*1.0d-6, wlc,rhl,         &   !hygroscopic OC
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==5) then
                call mie_gocart('ssa', q_unit*1.0d-6, wlc,rhl,         &   !sea-salt accumulation mode
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==6) then
                call mie_gocart('ssc', q_unit*1.0d-6, wlc,rhl,         &   !sea-salt coarse mode
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==7) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du1', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 1
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==8) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du2', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 2 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==9) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du3', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 3 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==10) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du4', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 4 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==11) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du5', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 5
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==12) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du6', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 6
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==13) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du7', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 7
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==14) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du8', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 8
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            else
               call stop_sdsu(myrank,'MSG makelut_broad_aerosol: there is no such nspec')
            endif NSPEC_GOCART


         case('UCD')

            if(nspec==1) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du1', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 1
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==2) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du2', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 2 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==3) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du3', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 3 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==4) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du4', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 4 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==5) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du5', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 5
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            else
               call stop_sdsu(myrank,'MSG makelut_broad_aerosol: there is no such nspec')
            endif 

         end select



!
! summing sub-band single-scattering component
!
     !
     ! avoid NaN case
     !
     if( isnan(cextp_lut) .or. cextp_lut == undefined ) then
        cextp_lut = 0.0
        cabsp_lut = 0.0
          php_lut = 0.0
     endif

     cextp_sum = cextp_sum + cextp_lut * wgt ! weigting sum of extinction cross section [cm2/cm3] or [cm-1]
     cabsp_sum = cabsp_sum + cabsp_lut * wgt ! weigting sum of absorption cross section [cm2/cm3] or [cm-1]
     php_sum   =   php_sum +   php_lut * wgt ! weigting sum of volume scattering phase function [-]
     wgt_sum   =   wgt_sum +             wgt ! weightng sum

 enddo WAV_LOOP

!
! derive band-mean values of single-scattering parameters
! Note these values are per unit mixing ratio
!
     cextp_lut = cextp_sum / wgt_sum ! broad-band extinction cross section [cm2/cm3] or [cm-1]
     cabsp_lut = cabsp_sum / wgt_sum ! broad-band absorption cross section [cm2/cm3] or [cm-1]
       php_lut =   php_sum / wgt_sum ! broad-band volume scattering phase function [-]


     kext_q_unit = cextp_lut * 1.d+5      ! [km-1] <- [cm-1]

     if( cextp_lut == 0. ) then
       omega = 0.
       asym  = 0.
     else 
       omega = (cextp_lut-cabsp_lut) / cextp_lut  ! signle scattering albedo [-]

!
! Compute aymetry parameter here
!
       phase_func  =  php_lut/(cextp_lut-cabsp_lut) !  
       cs(1:nang) = cos(ang(1:nang)*rad)
       yy(1:nang) = phase_func(1:nang)*unit_depth ! phsf(1:nang,l)
       mxlgn1 = N_stream*2 + 1  !same as nda 
       call lgndf3(mxlgn1,nang,cs,yy,gg)
       asym = gg(2)/gg(1)  ! asymetry paramer
     endif

!
! For RH index greather than 1 and dust species, SS parameters are constatnt (not RH function).
!
    if(i==1 .and. nspec>=7) then
       kext_q_unit_previous=kext_q_unit ; omega_previous=omega ; asym_previous=asym
    endif

    if(i>1 .and. nspec>=7) then
       kext_q_unit=kext_q_unit_previous ; omega=omega_previous ; asym=asym_previous
    endif

     !
     ! large array output
     !
     kext_q_unit_aero_array(i, nspec, ib, isl) = real(kext_q_unit)
     omega_aero_array      (i, nspec, ib, isl) = real(omega)
     asym_aero_array       (i, nspec, ib, isl) = real(asym)


 enddo RH_LOOP

! close(iounit) !close LUT table
! endif INQ

 enddo SPC_LOOP

 enddo BAND_LOOP

 enddo SWLW_LOOP


 90  Format ( 1P3E14.6 )

 return
 end subroutine makelut_broad_aerosol

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine compute_subw( ws,we,subwave )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute sub-wave point for making LUTs.
! History: 
!  10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: ws !starting wavelength/wavenumber
 real(sdsu_fps),intent(in) :: we !ending wavelength/wavenumber
 real(sdsu_fps),intent(out) :: subwave(numband_lut)
 integer :: i
 real(sdsu_fps) :: dsub  !interval of sub point

 if( numband_lut <= 2 ) call stop_sdsu(myrank, 'MSG compute_subw : numband_lut should be greater than 2')

 subwave(1) = ws              !starting wavelength/wavenumber
 subwave(numband_lut) = we    !ending wavelength/wavenumber

 dsub = (we - ws)/ (numband_lut -1)

 do i = 2 , numband_lut - 1
    subwave(i) = ws + dsub * REAL(i-1) 
 enddo

 return
 end subroutine compute_subw

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_opt_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

