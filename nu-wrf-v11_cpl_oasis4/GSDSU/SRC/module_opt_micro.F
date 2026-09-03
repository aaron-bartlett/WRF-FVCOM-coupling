#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_opt_micro
 use module_simulator
 use module_scatdb
 use module_tmatrix_database
 use module_gmp_4ice
 use module_mpi
 use module_openssp_ext
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
!
! Comments:  
!   This module deal with microwave particle single scattering and gaseous absorption for radar and
!   microwave simulators. Single-scattering LUTs are pre-computed by integrating spectral
!   parameters using Mie routines. Cloud/rain PSDs assumptions and density are maded consistent to visir 
!   and broadband optical modules, although subroutines are different. Aerosols are not accounted for
!   , because of their micron size and unknown refractive index with respect to this spectrum range.
! 
! History: 
!  10/2017  Toshi Matsui@NASA GSFC : Replaced LUT by new LUT array. 
!  06/2011  Toshi Matsui@NASA GSFC : Morrison two-moment class is added. 
!  06/2010  Toshi Matsui@NASA GSFC ; Significant clean up / re-organization of the routine. 
!  03/2010  Toshi Matsui@NASA GSFC ; Additional options for effective refractive index for solid particles.
!  08/2008  Toshi Matsui@NASA GSFC ; Unifield Radar & Microwve Part into this unifile module.    
!  02/2008  Toshi Matsui@NASA GSFC : Add 4 different melting options for microwave dielectric functions.
!  08/2007  Toshi Matsui@NASA GSFC ; Add GCE SBM options.
!  03/2007  Toshi Matsui@NASA GSFC ; Initial.
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
 public :: makelut_micro, & ! create microwave single-scattering properties LUTs
                opt_micro, & ! compute microwave single-scattering properties
                 watoptic, &    ! compute water dielectric constant for a given microwave frequency
               sample_psd_sbm,&  ! compute mean PSD for specif HID class
               effective_permittivity, & !compute effective ermittivity in Rayleigh approximation
               gas_absorb, iceoptic, mg_ellips, em_ellips, snow_compaction 

 public :: opt_micro_parameter

!
! 3D optical properties (public)
!
 type opt_micro_parameter 
   real(sdsu_fps) :: kext  ! total extinction coef [km-1]
   real(sdsu_fps) :: salb  ! total single scatterling albedo [-]
   real(sdsu_fps) :: asym  ! total asymetery parameters [-]
   real(sdsu_fps) :: sback ! total back scattering [km-1]
   real(sdsu_fps) :: vdop  ! Doppler velocity [cm/s]
 end type  opt_micro_parameter
 type (opt_micro_parameter),allocatable,public :: atmos_opt(:,:,:,:)

 real(sdsu_fps),allocatable,public :: t_air_micro(:,:,:) ! boundary temperature for scan simulator
 real(sdsu_fps),allocatable,public :: t_skin_micro(:,:)  ! boundary temperature for scan simulator

 integer,public :: mxfreq       !maximum # of frequency
 integer,public :: mxlyr_refine !maximum layer number for refined profile



! minimum condensate mixing ratio for microwave optical properties
  real(sdsu_fps),parameter :: q_min_micro = q_min_condensate  ! minimum detectable condensate amount from CRM [g/cm3] 

!  integer, public :: ibs,ibe,jbs,jbe !extra buffer domain bounds

! minimum extinction for stability
  real(sdsu_fps),parameter :: min_bext = 1.e-10   ! 1e-6 default


!
! universal microwave LUT
!
  real(sdsu_fps),allocatable :: &
     kext_q_unit_array (:,:,:,:), &  ! extinction coefficient per unit mass conc [1/km / (g/m3)]
     salb_array        (:,:,:,:), &  ! single-scatter albedo [-]
     asym_array        (:,:,:,:), &  ! asymmetry factor [-]
     pbck_array        (:,:,:,:), &  ! backscatter phase function/(4*pi) [-]
     vdop_array        (:,:,:,:)     ! Doppler velocity [cm/s]


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_micro(sensor,within_swath,max_inc_angle,failed_opt_micro)
 implicit none 
!---------------------------------------------------------------------------------------------------
!
! Comments:  
!  This module perpare the microwave (passive and active) optical properties
!  for a given atmosphere and condensates amounts. 
! 
! History: 
!  08/2008  Toshi Matsui@NASA GSFC ; Unifield Radar & Microwve Part into this unifile module.    
!  02/2008  Toshi Matsui@NASA GSFC : Add 4 different melting options for microwave dielectric functions.
!  08/2007  Toshi Matsui@NASA GSFC ; Add GCE SBM options.
!           
! References:
!
! Masunaga, H., and C.D. Kummerow, 2005: Combined Radar and Radiometer Analysis of 
!     Precipitation Profiles for a Parametric Retrieval Algorithm. J. Atmos. Oceanic 
!     Technol., 22, 909-929.
!
! 
!---------------------------------------------------------------------------------------------------
 character(len=5),intent(in) :: sensor  !'radar' or 'micro'
 logical,intent(in),optional :: within_swath(:,:) ! if within swath , true 
 real,intent(in),optional :: max_inc_angle   ! maximum inclination angle [deg]
 logical,intent(out),optional :: failed_opt_micro  !failure signal
 integer :: ierr !error index
 integer :: i,j,k,nf,m !looping indice
 integer :: is,ie,js,je,ks,ke,di,dj,dk
 real(sdsu_fps) :: atm_ext  !atmosphere gas extinction
 real(sdsu_fps) :: freqcy  ! frequency [GHz]
 real(sdsu_fps) :: tavg, tavg_c, pres !temperature [K], [C], and pressurej[]
 real(sdsu_fps) :: melt_frac_bulk ! parameterized melting fraction for freezing condensates of bulk microphysics
 integer :: ispc  !RAMS spicies index
 integer :: ns,nt !temperature and RH index
 real(sdsu_fps) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(sdsu_fps),allocatable :: hgt_stag(:) ![km]
 integer :: delta_ij
 real(sdsu_fps) :: incidence ! radiation incidence angle [rad]
 real :: dummy
 type (opt_micro_parameter) :: atmos_opt1d(mxlyr)  !1D column of optical properties

! type (particle_hid),allocatable :: sback_hid1d(:)         !backscatter for each HID 
! type (particle_hid),allocatable :: sback_hid1d_refine(:)  !backscatter for each HID 


 type ( particle_gen ) :: & !particle_gen is defined in module_simulater
     kext_gen,  & ! extinction coefficient [km-1]
     salb_gen,  & ! single scattering albedo [-]
     asym_gen,  & ! asymetry parameter [-]
     pbck_gen,  & ! backscattering coefficient [-]
     vdop_gen     ! Doppler velocity [m/s]

 type ( particle_rams ) :: & !particle_rams is defined in module_simulater
     kext_rams,  & ! extinction coefficient [km-1]
     salb_rams,  & ! single scattering albedo [-]
     asym_rams,  & ! asymetry parameter [-]
     pbck_rams,  & ! backscattering coefficient [-]
     vdop_rams     ! Doppler velocity [m/s]

 type ( particle_sbm ) :: & !particle_sbm is defined in module_simulater
     kext_sbm,  & ! extinction coefficient [km-1]
     salb_sbm,  & ! single scattering albedo [-]
     asym_sbm,  & ! asymetry parameter [-]
     pbck_sbm,  & ! backscattering coefficient [-]
     vdop_sbm     ! Doppler velocity 

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1

!
! Program Start
!
 if( allocated(hgt_stag) ) deallocate( hgt_stag )


  if( present(failed_opt_micro) ) failed_opt_micro = .false.

!
! determine some memory bounds
!
 select case(trim(sensor))
 case('radar')
    mxfreq = mxfreq_radar
    mxlyr_refine = mxlyr_radar
    if(masterproc) print*,'MSG opt_microwave: Compute radar backscatter '
 case('micro')
    mxfreq = mxfreq_micro
    mxlyr_refine = mxlyr_micro
    if(masterproc) print*,'MSG opt_microwave: Compute microwave single-scattering properties  '
 case default
    call stop_sdsu(myrank,'MSG opt_micowave: no such sensor')
 end select

!
! dealloacte 3D optical parameters
!
  if( allocated(atmos_opt)) deallocate( atmos_opt )

!
! allocate 3D optical parameters
!
 select case(trim(sensor))
 case('micro')

   if( scan_micro .or. slant_path_micro ) then  ! This option require domain allocation and constant profile

      !
      ! extra buffer zone to be used in slant path options
      !
      if( scan_micro ) then
        delta_ij =  int(  ( maxval( hgt_stag_micro )  / cos(max_inc_angle*const_degrad) &
                    * cos( ( 90.-max_inc_angle )*const_degrad )  * 1.25 ) / gridsize )
      else  !for slant option for non scan_micro

        call derive_incidence_angle( ona_angle_micro*const_degrad, sma_micro, 90.e0*const_degrad, &
                                     incidence,  dummy ) !90deg to overestimate incidence for this portion. 

        delta_ij =  int(  ( maxval( hgt_stag_micro )  / cos(incidence) &
                    * cos( ( 90.-incidence/const_degrad)*const_degrad )  * 1.25 ) / gridsize )
!        delta_ij =  int(  ( maxval( hgt_stag_micro )  / cos(inc_angle_micro*const_degrad) &
!                    * cos( ( 90.-inc_angle_micro )*const_degrad )  * 1.25 ) / gridsize )

      endif

      !extra domain for slant path
      call get_extra_memory_domain(delta_ij, is,ie,js,je,ibs,ibe,jbs,jbe)

      allocate( atmos_opt (ibs:ibe,jbs:jbe,mxlyr_refine,mxfreq), &  !buffer domain
                stat=ierr )
       if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_micro: allocation error atmos_opt -> Terminate program. ')
       
      if( .not. allocated(t_air_micro) ) then !domain memory air temp and t_skin
        allocate( t_air_micro (ibs:ibe,jbs:jbe,mxlyr_refine),& !buffer domain
                  t_skin_micro(ibs:ibe,jbs:jbe             ),& !surface data is domain memory
                  stat=ierr )
        t_skin_micro = undefined

        if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_micro: allocation error t_air_micro-> Terminate program. ')
      endif

      allocate( hgt_stag(0:mxlyr_refine) )
      hgt_stag(0:mxlyr_refine) =  hgt_stag_micro(0:mxlyr_refine)

   else                   ! This option allow memory allocation and CRM profile (column simulation)
      allocate( atmos_opt (is:ie,js:je,1:mxlyr,mxfreq), stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_micro: allocation error atmos_opt-> Terminate program. ')

   endif

  case ('radar')

   if( scan_radar ) then  ! This option require domain allocation and instrument profile

      !
      ! extra buffer zone to be used in slant path options
      !
      delta_ij =  max(1 , int(  ( maxval( hgt_stag_radar )  / cos(max_inc_angle*const_degrad) &
                          * cos( ( 90.-max_inc_angle )*const_degrad )  * 2. ) / gridsize ) )

      !extra domain for slant path or footprint
      call get_extra_memory_domain(delta_ij, is,ie,js,je,ibs,ibe,jbs,jbe)
      allocate( atmos_opt (ibs:ibe,jbs:jbe,mxlyr_refine,mxfreq), &  !buffer domain
                stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_micro: allocation error atmos_opt. ')

      allocate( hgt_stag(0:mxlyr_refine) )
      hgt_stag(0:mxlyr_refine) =  hgt_stag_radar(0:mxlyr_refine)

      !toshiii --> l2bundle_radar
      if( .not. allocated(L2_radar_1d) ) then !domain memory air temp and t_skin
        allocate( L2_radar_1d    (mxlyr)                        , &
                  L2_radar_refine(ibs:ibe,jbs:jbe,mxlyr_refine) , & !buffer domain
                  stat=ierr  )
       endif

   else

        allocate( atmos_opt (is:ie,js:je,1:mxlyr_refine,mxfreq), stat=ierr )
        if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_micro: allocation error atmos_opt. ')

      if( .not. allocated(L2_radar_1d) ) then !buffered sub-domain 
        allocate( L2_radar_1d    (mxlyr)                        , &
                  L2_radar_refine(ibs:ibe,jbs:jbe,mxlyr_refine) , & !buffer domain
                  L2_radar_conv  (ibs:ibe,jbs:jbe,mxlyr_refine)   , & !buffer domain
                  stat=ierr  )
      endif

!      if( .not. allocated(sback_hid1d       )) allocate(sback_hid1d       (mxlyr)) !HID
!      if( .not. allocated(sback_hid1d_refine)) allocate(sback_hid1d_refine(mxlyr_refine)) !HID

      allocate( hgt_stag(0:mxlyr_refine) )
        hgt_stag(0:mxlyr_refine) =  hgt_stag_radar(0:mxlyr_refine)

   endif
 end select


!
! allocate LUT parameters
!
  if( lut_micro ) then

   !
   ! read microwave LUT
   !
   if( trim(sensor) == 'micro') then
      call read_lut_micro_nc(micro_sensor)
   elseif( trim(sensor) == 'radar' ) then
      call read_lut_micro_nc(radar_sensor)
   endif


  endif !lut_micro


!
! Loop over the channel frequencies
!
 freq_loop: do nf = 1, mxfreq 

!
! Sensor frequency and wavelengths
!
  if(sensor=='radar') freqcy = freq_radar(nf) ! radar
  if(sensor=='micro') freqcy = freq_micro(nf) ! passive microwave
  if(verbose_SDSU .and. masterproc)  &
     print*, ' Frequency =', freqcy,' GHz -> Wavelength = ',2.997925e2/freqcy ,'mm'

!
! Calculate optical properties of each volume element of the cloud.
!

 DOM_2D_LOOP: do j = myj_start, myj_end ; do i = myi_start, myi_end

 if( present( within_swath ) .and. sensor == 'radar' ) then
     if( .not. within_swath(i,j) ) cycle   !skip for fast processing
 endif


!
! special MMF_GPROF skip 
!
 if( trim(sim_case) == 'MMF_GPROF' .and. j > max_sample_gprof ) then
     atmos_opt(i,j,:,nf)%kext  = undefined
     atmos_opt(i,j,:,nf)%salb  = undefined
     atmos_opt(i,j,:,nf)%asym  = undefined
     atmos_opt(i,j,:,nf)%sback = undefined
     atmos_opt(i,j,:,nf)%vdop  = undefined
     cycle
 endif


 VERTICAL_LOOP: do k = 1, mxlyr
    tavg = atmos(i,j,k)%t_air  ! [K]
    tavg_c = tavg - 273.16  ! [degC]
    pres = atmos(i,j,k)%press ! [hPa]

!
! Ideintify melting layer and fraction for bulk microphysics (SBM has explicit melting process.) 
!
  if( melt_opt /= 0 ) then   !debbug (toshii)
    if ( atmos_stag(i,j,k-1)%t_air == 273.16 .and. atmos_stag(i,j,k)%t_air == 273.16) then
        melt_frac_bulk = 0.  !debug found by Jiundar Chern
    elseif(atmos_stag(i,j,k-1)%t_air >= 273.16 .and. atmos_stag(i,j,k)%t_air   <= 273.16) then
        melt_frac_bulk = ( atmos_stag(i,j,k-1)%t_air - 273.16 ) / &
                         ( atmos_stag(i,j,k-1)%t_air - atmos_stag(i,j,k)%t_air )
    elseif( atmos_stag(i,j,k-1)%t_air < 273.16 ) then ! no melt
        melt_frac_bulk = 0. 
    elseif( atmos_stag(i,j,k)%t_air > 273.16 ) then
        melt_frac_bulk = 0.9
    endif
  else
     melt_frac_bulk = 0.
  endif

!
! Gas absorption code 
!
    call gas_absorb(freqcy, tavg, pres, atmos(i,j,k)%rh, 0., atm_ext)


!
! Condensate particle optical properties 
!

  mic_select1: select case(trim(type_microphysics))

  case('GEN')

       !
       ! compute optical properties for different species of Morrison microphysics
       !
       call opt_gen_micro( nf, freqcy, tavg, melt_frac_bulk, mu_gen, q_gen(i,j,k), re_gen(i,j,k), &
                           kext_gen, salb_gen, asym_gen, pbck_gen, vdop_gen )


       !
       ! Total optical properties
       !
       call total_opt_gen(atmos(i,j,k)%press, atmos(i,j,k)%w, & 
                          atm_ext, kext_gen, salb_gen, asym_gen, pbck_gen, vdop_gen, &
                          atmos_opt1d(k) )

     !
     ! radar L2 parameters
     !
!     if(  sensor == 'radar' .and. scan_radar .and. nf == 1) then
     if(  sensor == 'radar' .and. nf == 1) then
           call L2_radar_parameters_gen(mu_gen, q_gen(i,j,k), re_gen(i,j,k), &
                                        atmos(i,j,k) , L2_radar_1d(k)  )
     endif

  case('RAMS')

     nt = max( 1,min( 31,-nint(  (atmos(i,j,k)%t_air-273.15)  )))  ! Temperature index--1 to 31 T below 0 deg.C
     ns = max(79,min(100, nint( atmos(i,j,k)%rh          )))       ! RH index--79 to 100 based on integer RH
     !
     ! compute optical properties for different species of RAMS microphysics
     !
     call opt_rams_micro( nt, ns, nf, freqcy, tavg, melt_frac_bulk, q_rams(i,j,k), re_rams(i,j,k), &
                           kext_rams, salb_rams, asym_rams, pbck_rams, vdop_rams )

     !
     ! Total optical properties
     !
     call total_opt_rams(atmos(i,j,k)%rho_dair*1.e-3, atmos(i,j,k)%w, &
                         atm_ext, kext_rams, salb_rams, asym_rams, pbck_rams, vdop_rams, &
                         atmos_opt1d(k) )
     atmos_opt1d(k)%vdop     = 0.e0


  case('SBM')

     !
     ! compute optical properties for different species of HUCM SBM microphysics
     !
     call opt_sbm_micro( freqcy, tavg, q_sbm(i,j,k), n_sbm(i,j,k,1:nbin), fmelt_sbm_snow(i,j,k, 1:nbin), &
                         frime_sbm_snow(i,j,k, 1:nbin),  fmelt_sbm_graupel(i,j,k, 1:nbin), &
                         fmelt_sbm_hail(i,j,k, 1:nbin), kext_sbm, salb_sbm, asym_sbm, pbck_sbm, vdop_sbm )

     !
     ! Total optical properties
     !
     call total_opt_sbm( atmos(i,j,k)%press, atmos(i,j,k)%w, &
                         atm_ext, kext_sbm, salb_sbm, asym_sbm, pbck_sbm, vdop_sbm, &
                         atmos_opt1d(k) )

     !
     ! radar L2 parameters
     !
     if(  sensor == 'radar' .and. nf == 1) then
           call L2_radar_parameters_sbm(q_sbm(i,j,k), n_sbm(i,j,k,1:nbin)  , &
                                        atmos(i,j,k), L2_radar_1d(k)  )
     endif

     !
     ! HID function
     !
     if( sensor == 'radar' .and. radar_hid ) then
       call sback_hid_sbm( freqcy, tavg, q_sbm(i,j,k), n_sbm(i,j,k,1:nbin), fmelt_sbm_snow(i,j,k, 1:nbin), &
                           frime_sbm_snow(i,j,k, 1:nbin), fmelt_sbm_graupel(i,j,k, 1:nbin), &
                           fmelt_sbm_hail(i,j,k, 1:nbin), sback_hid1d(k) )
     endif

  case default  !for no-particle simulation (sim_case == 'SOUNDING', 'MLM' or 'LIS')

     atmos_opt1d(k)%kext  = atm_ext
     atmos_opt1d(k)%salb  = 0.e0
     atmos_opt1d(k)%asym  = 0.e0
     atmos_opt1d(k)%sback = 0.e0
     atmos_opt1d(k)%vdop  = 0.e0

  end select mic_select1


!
! check error
!
 if(trim(sim_case) == 'MMF_GPROF') then !Special output for  GPROF  
   if    ( ISNAN(atmos_opt1d(k)%kext) ) then
      print*,'MSG: opt_micro: find NAN in kex  at point (i,j,k,nf)=',i,j,k,nf
      print*,'MSG: opt_micro:Strange input file = ',sdsu_inp_name(1:efile_len)
       if( present(failed_opt_micro) ) failed_opt_micro = .true. ; return
   elseif( ISNAN(atmos_opt1d(k)%salb) ) then
      print*,'MSG: opt_micro: find NAN in salb at point (i,j,k,nf)=',i,j,k,nf
      print*,'MSG: opt_micro:Strange input file = ',sdsu_inp_name(1:efile_len)
       if( present(failed_opt_micro) ) failed_opt_micro = .true. ; return
   elseif( ISNAN(atmos_opt1d(k)%asym) ) then
      print*,'MSG: opt_micro: find NAN in asym at point (i,j,k,nf)=',i,j,k,nf
      print*,'MSG: opt_micro:Strange input file = ',sdsu_inp_name(1:efile_len)
       if( present(failed_opt_micro) )  failed_opt_micro = .true. ; return
  endif

 else
   if    ( ISNAN(atmos_opt1d(k)%kext) ) then
      print*,'MSG: opt_micro: find NAN in kex  at point (i,j,k,nf)=',i,j,k,nf ; call stop_sdsu(myrank,'')
   elseif( ISNAN(atmos_opt1d(k)%salb) ) then
      print*,'MSG: opt_micro: find NAN in salb at point (i,j,k,nf)=',i,j,k,nf ; call stop_sdsu(myrank,'')
   elseif( ISNAN(atmos_opt1d(k)%asym) ) then
      print*,'MSG: opt_micro: find NAN in asym at point (i,j,k,nf)=',i,j,k,nf ; call stop_sdsu(myrank,'')
   elseif( ISNAN(atmos_opt1d(k)%sback) ) then
      print*,'MSG: opt_micro: find NAN in sback at point (i,j,k,nf)=',i,j,k,nf ; call stop_sdsu(myrank,'')
   elseif( ISNAN(atmos_opt1d(k)%vdop    ) ) then
      print*,'MSG: opt_micro: find NAN in vdop  at point (i,j,k,nf)=',i,j,k,nf ; call stop_sdsu(myrank,'')
  endif
 endif


 enddo VERTICAL_LOOP

!
! re-sample vertical profile of optical properties with constant finer vertical layers.
!

 if( (sensor == 'micro' .and. scan_micro) .or. &
     (sensor == 'micro' .and. slant_path_micro) .or. &
     (sensor == 'radar') ) then ! re-sample vertical profile with constant layer height


    call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                           mxlyr_refine,  hgt_stag, &
                           atmos_opt1d(1:mxlyr)%kext, atmos_opt(i,j,1:mxlyr_refine,nf)%kext   )

    call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                           mxlyr_refine,  hgt_stag, &
                           atmos_opt1d(1:mxlyr)%salb, atmos_opt(i,j,1:mxlyr_refine,nf)%salb  )

    call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                           mxlyr_refine,  hgt_stag, &
                           atmos_opt1d(1:mxlyr)%asym, atmos_opt(i,j,1:mxlyr_refine,nf)%asym  )

    call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                           mxlyr_refine,  hgt_stag, &
                           atmos_opt1d(1:mxlyr)%sback, atmos_opt(i,j,1:mxlyr_refine,nf)%sback )

    call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                           mxlyr_refine,  hgt_stag, &
                           atmos_opt1d(1:mxlyr)%vdop , atmos_opt(i,j,1:mxlyr_refine,nf)%vdop  )


    !
    ! microwave anscillary parameters
    !
    if( (sensor == 'micro' .and. scan_micro       .and. nf == 1 ) .or. &
        (sensor == 'micro' .and. slant_path_micro .and. nf == 1 ) ) then 

        call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                              mxlyr_refine,  hgt_stag, &
                              atmos(i,j,1:mxlyr)%t_air, t_air_micro(i,j,1:mxlyr_refine)  )

        t_skin_micro(i,j) = surface(i,j)%t_skin   ! [K]

    endif

    !
    ! radar L2 parameters
    !
    RADAR_L2: if(  sensor == 'radar' .and. nf == 1) then

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%press, L2_radar_refine(i,j,1:mxlyr_refine)%press)

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%t_air, L2_radar_refine(i,j,1:mxlyr_refine)%t_air)

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%sh, L2_radar_refine(i,j,1:mxlyr_refine)%sh)

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%w   , L2_radar_refine(i,j,1:mxlyr_refine)%w)

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%iwc, L2_radar_refine(i,j,1:mxlyr_refine)%iwc)  

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%lwc, L2_radar_refine(i,j,1:mxlyr_refine)%lwc)

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%rwc, L2_radar_refine(i,j,1:mxlyr_refine)%rwc)

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%lfrac, L2_radar_refine(i,j,1:mxlyr_refine)%lfrac)

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%bulk_den, L2_radar_refine(i,j,1:mxlyr_refine)%bulk_den)

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 L2_radar_1d(1:mxlyr)%rain_rate, L2_radar_refine(i,j,1:mxlyr_refine)%rain_rate)

          do m = 0, 6
             call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                    mxlyr_refine,  hgt_stag, &
                                    L2_radar_1d(1:mxlyr)%mmt(m), L2_radar_refine(i,j,1:mxlyr_refine)%mmt(m) )
          enddo


          if( sensor == 'radar' .and. radar_hid ) then

          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%dz, sback_hid1d_refine(1:mxlyr_refine)%dz )
          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%rn, sback_hid1d_refine(1:mxlyr_refine)%rn )
          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%cr, sback_hid1d_refine(1:mxlyr_refine)%cr )
          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%ldg, sback_hid1d_refine(1:mxlyr_refine)%ldg )
          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%hdg, sback_hid1d_refine(1:mxlyr_refine)%hdg )
          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%vi, sback_hid1d_refine(1:mxlyr_refine)%vi )
          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%ds, sback_hid1d_refine(1:mxlyr_refine)%ds )
          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%ws, sback_hid1d_refine(1:mxlyr_refine)%ws )
          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%ha, sback_hid1d_refine(1:mxlyr_refine)%ha )
          call refine_range_bin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr)%dhgt, &
                                 mxlyr_refine,  hgt_stag, &
                                 sback_hid1d(1:mxlyr)%bd, sback_hid1d_refine(1:mxlyr_refine)%bd )

          if( sensor == 'radar' .and. radar_hid .and. trim(cloud_microphysics(1:4)) == 'HUCM') then
            call find_hid( sback_hid1d_refine(1:mxlyr_refine), L2_radar_refine(i,j,1:mxlyr_refine)%hid, &
                           L2_radar_refine(i,j,1:mxlyr_refine)%hid_conf  )  
          endif

          endif

    endif RADAR_L2


 else  !use vertical profile identical to the CRM vertical profile .   (micro column option only)

     atmos_opt(i,j,:,nf)%kext  = atmos_opt1d(:)%kext
     atmos_opt(i,j,:,nf)%salb  = atmos_opt1d(:)%salb
     atmos_opt(i,j,:,nf)%asym  = atmos_opt1d(:)%asym
     atmos_opt(i,j,:,nf)%sback = atmos_opt1d(:)%sback
     atmos_opt(i,j,:,nf)%vdop  = atmos_opt1d(:)%vdop

 endif

 enddo ; enddo DOM_2D_LOOP

 enddo freq_loop

!
! write microwave optical properties
!
 if(sensor == 'radar' .and. write_opt) call write_opt_radar

 if(sensor == 'micro' .and. write_opt) call write_opt_micro

#if MPI == 2
if(numproc_tot>1) then
!
! broadcasting to all threads for scan simulator
!
 if( (sensor == 'micro' .and. scan_micro      ) .or. &
     (sensor == 'micro' .and. slant_path_micro)  ) then

  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, atmos_opt(ibs:ibe,jbs:jbe,:,:)%kext  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, atmos_opt(ibs:ibe,jbs:jbe,:,:)%salb  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, atmos_opt(ibs:ibe,jbs:jbe,:,:)%asym  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, t_air_micro(ibs:ibe,jbs:jbe,:)       )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, t_skin_micro(ibs:ibe,jbs:jbe)        )

 endif


 if ( sensor == 'radar' .and. scan_radar ) then

  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, atmos_opt(ibs:ibe,jbs:jbe,:,:)%kext  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, atmos_opt(ibs:ibe,jbs:jbe,:,:)%sback )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, atmos_opt(ibs:ibe,jbs:jbe,:,:)%vdop     )
 endif

 if( sensor == 'radar' ) then
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%press)
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%t_air)
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%sh   )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%w    )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%iwc  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%lwc  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%rwc  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%lfrac)
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%bulk_den  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%rain_rate )
  do m = 0, 6
     call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%mmt(m) )
  enddo
  if( sensor == 'radar' .and. radar_hid ) then
     call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%hid )
     call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, L2_radar_refine(ibs:ibe,jbs:jbe,:)%hid_conf )
  endif

 endif


endif
#endif

 return
 end subroutine opt_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_lut_micro_nc(sensor_name)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Open and read the prepared Look-up table. 
! 
! History: 
!   10/2017  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 character(len=20),intent(in) :: sensor_name
 integer :: max_spc
 integer :: ierr, ncid, varid

!
! 
!
 select case(trim(type_microphysics))
 case('GEN')
   max_spc = 6
 case('RAMS')
   max_spc = 16
 case('SBM')
   call stop_sdsu(myrank,'Not yet supporoting SBM')
 case default
 end select

!
! allocate array parameters
!
 if( .not. allocated(kext_q_unit_array) ) then
   allocate( kext_q_unit_array ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             salb_array        ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             asym_array        ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             pbck_array        ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             vdop_array        ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             stat=ierr )
   if (ierr /= 0) call stop_sdsu(myrank,'MSG read_lut_micro_nc: allocation error -> Terminate program.')
 endif
   kext_q_unit_array = 0.
   salb_array = 0.
   asym_array = 0.
   pbck_array = 0.
   vdop_array = 0.

!
! open netcdf file
!
 sdsu_io_file = trim(sdsu_dir_sslut)//'MICRO_CLOUD_'//trim(cloud_microphysics)//'-'//trim(sensor_name)//'.nc'
 if(masterproc) print*, 'MSG read_lut_micro_nc: reading microwave LUT:',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! extinction coefficient per unit mass conc [km-1 / (g/m3)]
!
 call check( nf90_inq_varid(ncid, 'kext_q_unit' , varid ))
 call check( nf90_get_var(ncid, varid, kext_q_unit_array ))

!
! single scattering albedo [-]
!
 call check( nf90_inq_varid(ncid, 'salb' , varid ))
 call check( nf90_get_var(ncid, varid, salb_array ))

!
! asymetry parameter [-]
!
 call check( nf90_inq_varid(ncid, 'asym' , varid ))
 call check( nf90_get_var(ncid, varid, asym_array ))

!
! backscatter phase function/(4*pi) [-]
!
 call check( nf90_inq_varid(ncid, 'pbck' , varid ))
 call check( nf90_get_var(ncid, varid, pbck_array ))

!
! Doppler velocity [cm/s]
!
 call check( nf90_inq_varid(ncid, 'vdop' , varid ))
 call check( nf90_get_var(ncid, varid, vdop_array ))

!
! close netCDF file
!
 call check( nf90_close(ncid) )

 return
 end subroutine read_lut_micro_nc

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine interp_lut_array(nf, idx_spc, q_in, re_in, wgt1_temp , wgt2_temp, idx_temp, &
                             kext, salb , asym, pbck, vdop )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  interpolate LUT array for given index.
! 
! History: 10/2017  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 integer,intent(in) :: nf  !channel index
 integer,intent(in) :: idx_spc  ! species index
 real(sdsu_fps),intent(in) :: q_in !input mixing ratio [g/m3]
 real(sdsu_fps),intent(in) :: re_in !input effective radius [micron]
 real(sdsu_fps),intent(in) :: wgt1_temp , wgt2_temp ! temperautre weight
 integer,intent(in) :: idx_temp  ! temperature index

 real(sdsu_fps),intent(out) :: kext, salb , asym, pbck, vdop
 
 integer :: i,j
 real(sdsu_fps) :: wgt1, wgt2
 real(sdsu_fps) :: re


 if(q_in <= q_min_condensate .or. re_in == 0. ) then
   kext = 0. ; salb = 0. ; asym = 0. ; pbck = 0. ; vdop = 0.
 else
   re = re_in
   if(re_in < pts_re(1)             ) re = pts_re(1)
   if(re_in > pts_re(mxpts_re)) re = pts_re(mxpts_re)

   do j = 1, mxpts_re-1
      if( re >= pts_re(j) .and. re <= pts_re(j+1) ) then
         wgt2 =  ( re - pts_re(j) )    / ( pts_re(j+1) - pts_re(j) )
         wgt1 = 1.0 - wgt2

         kext =  q_in * (&
                        wgt1*( kext_q_unit_array ( j, idx_temp  , idx_spc, nf ) * wgt1_temp &
                              +kext_q_unit_array ( j, idx_temp+1, idx_spc, nf ) * wgt2_temp ) + &
                        wgt2*( kext_q_unit_array ( j+1, idx_temp  , idx_spc, nf ) * wgt1_temp &
                              +kext_q_unit_array ( j+1, idx_temp+1, idx_spc, nf ) * wgt2_temp ) &
                        )

         salb =  wgt1*( salb_array ( j  , idx_temp  , idx_spc, nf ) * wgt1_temp &
                       +salb_array ( j  , idx_temp+1, idx_spc, nf ) * wgt2_temp ) + &
                 wgt2*( salb_array ( j+1, idx_temp  , idx_spc, nf ) * wgt1_temp &
                       +salb_array ( j+1, idx_temp+1, idx_spc, nf ) * wgt2_temp )

         asym =  wgt1*( asym_array ( j  , idx_temp  , idx_spc, nf ) * wgt1_temp &
                       +asym_array ( j  , idx_temp+1, idx_spc, nf ) * wgt2_temp ) + &
                 wgt2*( asym_array ( j+1, idx_temp  , idx_spc, nf ) * wgt1_temp &
                       +asym_array ( j+1, idx_temp+1, idx_spc, nf ) * wgt2_temp )

         pbck =  wgt1*( pbck_array ( j  , idx_temp  , idx_spc, nf ) * wgt1_temp &
                       +pbck_array ( j  , idx_temp+1, idx_spc, nf ) * wgt2_temp ) + &
                 wgt2*( pbck_array ( j+1, idx_temp  , idx_spc, nf ) * wgt1_temp &
                       +pbck_array ( j+1, idx_temp+1, idx_spc, nf ) * wgt2_temp )

         vdop =  wgt1*( vdop_array ( j  , idx_temp  , idx_spc, nf ) * wgt1_temp &
                       +vdop_array ( j  , idx_temp+1, idx_spc, nf ) * wgt2_temp ) + &
                 wgt2*( vdop_array ( j+1, idx_temp  , idx_spc, nf ) * wgt1_temp &
                       +vdop_array ( j+1, idx_temp+1, idx_spc, nf ) * wgt2_temp )

!if( kext < 0. ) then
!    print*, j, idx_temp  , idx_spc, nf
!    print*, kext, wgt1, wgt2, wgt1_temp, wgt2_temp, j
!stop
!endif

         exit
      endif
   enddo


 endif

 return
 end subroutine interp_lut_array

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_gen_micro(nf, freq, tavg, fmelt, mu, q, re, kext_gen, salb_gen, asym_gen ,pbck_gen, vdop_gen)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation techniqu.
!  This routine is for GCE microphysics.
!    
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; New LUT array. 
!  06/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 integer,intent(in)::nf  !frequency id
 real(sdsu_fps),intent(in)::tavg  ! layer temperature
 real(sdsu_fps),intent(in)::freq  ! frequency [GHz]
 real(sdsu_fps),intent(in)::fmelt ! meltfraction 

 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     mu, & ! dispersion parameter
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]
      
 type ( particle_gen ),intent(out) :: & ! particle_gen is defined in module_simulater
     kext_gen,  & ! extinction coefficient [km-1]
     salb_gen,  & ! single scattering albedo [-]
     asym_gen,  & ! asymetry parameter [-]
     pbck_gen,  & ! backscatter [-]
     vdop_gen     ! Doppler velocity [m/s]

 integer :: ispc  !RAMS spicies index
 real(sdsu_fps) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(sdsu_fps) :: cfvt, pwvt ! alpha_vt & beta_vt

 logical,parameter :: new_method = .true.

 integer :: i,j
 integer :: idx_spc, idx_temp, i_i, i_w
 real(sdsu_fps) :: wgt1, wgt2
 real(sdsu_fps) :: re_in, q_in
 real(sdsu_fps) :: kext_q_unit , salb , asym, pbck, vdop
 real(sdsu_fps) :: tair  ! air temperature [K]
 real(sdsu_fps) :: wgt1_temp , wgt2_temp
 real(sdsu_fps) :: wgt1_wtemp, wgt2_wtemp
 real(sdsu_fps) :: wgt1_itemp, wgt2_itemp
 real(sdsu_fps) :: wgt1_re, wgt2_re

!
! Initialize
!
 call init_microphysics(kext_gen,'zero')
 call init_microphysics(salb_gen,'zero')
 call init_microphysics(asym_gen,'zero')
 call init_microphysics(pbck_gen,'zero')
 call init_microphysics(vdop_gen,'zero')


 if ( lut_micro ) then

!
! air temperature weights (liquid)
!
 tair = tavg 
 if (tavg  < pts_temp_w(1))           tair = pts_temp_w(1) 
 if (tavg  > pts_temp_w(mxpts_temp))   tair = pts_temp_w(mxpts_temp) 

 do i = 1, mxpts_temp-1
    if( tair >= pts_temp_w(i) .and. tair <= pts_temp_w(i+1) ) then
        wgt1_wtemp = ( pts_temp_w(i+1) - tair ) / ( pts_temp_w(i+1) - pts_temp_w(i) )
        wgt2_wtemp = 1.-wgt1_wtemp
        i_w = i
        if( wgt2_wtemp < 0. .or. wgt2_wtemp > 1.0 ) then
            print*, 'MSG opt_gen_micro strange wgt2_wtemp',wgt2_wtemp
            stop
        endif
        exit
    endif
 enddo

!
! air temperature weights (ice)
!
 tair = tavg
 if (tavg  < pts_temp_i(1))           tair = pts_temp_i(1)
 if (tavg  > pts_temp_i(mxpts_temp))   tair = pts_temp_i(mxpts_temp)

 do i = 1, mxpts_temp-1
    if( tair >= pts_temp_i(i) .and. tair <= pts_temp_i(i+1) ) then
        wgt1_itemp = ( pts_temp_i(i+1) - tair ) / ( pts_temp_i(i+1) - pts_temp_i(i) )
        wgt2_itemp = 1.-wgt1_itemp
        i_i = i
        if( wgt2_itemp < 0. .or. wgt2_itemp > 1.0 ) then
            print*, 'MSG opt_gen_micro strange wgt2_itemp',wgt2_itemp
            stop
        endif
        exit
    endif
 enddo

!
! cloud species
!
 call interp_lut_array(nf, 1, q%cloud, re%cloud, wgt1_wtemp, wgt2_wtemp, i_w, &
      kext_gen%cloud, salb_gen%cloud, asym_gen%cloud, pbck_gen%cloud, vdop_gen%cloud )

!
! rain species
!
 call interp_lut_array(nf, 2, q%rain, re%rain, wgt1_wtemp, wgt2_wtemp, i_w, &
      kext_gen%rain, salb_gen%rain, asym_gen%rain, pbck_gen%rain, vdop_gen%rain )

!
! ice species
!
 call interp_lut_array(nf, 3, q%ice, re%ice, wgt1_itemp, wgt2_itemp, i_i, &
      kext_gen%ice, salb_gen%ice, asym_gen%ice, pbck_gen%ice, vdop_gen%ice )

!
! snow species
!

 if(fmelt == 0.) then !snow
   call interp_lut_array(nf, 4, q%snow, re%snow, wgt1_itemp, wgt2_itemp, i_i, &
                         kext_gen%snow, salb_gen%snow, asym_gen%snow, pbck_gen%snow, vdop_gen%snow )
 else ! melting apply for snow only
   call mie_gen_micro('qs',freq,tavg,mu_gen%snow, q%snow,re%snow,rho_gen%snow , & !melting snow
                      fmelt, a_vt%snow, b_vt%snow, kext_gen%snow,salb_gen%snow, &
                      asym_gen%snow,pbck_gen%snow, vdop_gen%snow  )
 endif


!
! graupel species
!
 call interp_lut_array(nf, 5, q%graupel, re%graupel, wgt1_itemp, wgt2_itemp, i_i, &
      kext_gen%graupel, salb_gen%graupel, asym_gen%graupel, pbck_gen%graupel, vdop_gen%graupel )

!
! hail species
!
 call interp_lut_array(nf, 6, q%hail, re%hail, wgt1_itemp, wgt2_itemp, i_i, &
      kext_gen%hail, salb_gen%hail, asym_gen%hail, pbck_gen%hail, vdop_gen%hail )


 else  !no use of mie LUT (SLOW: compute extinction by integrating hydrometeor size spectra.)


     if( trim(cloud_microphysics) == 'WDM') then

        call mie_gen_micro('qc',freq,tavg,mu_gen%cloud,q%cloud,re%cloud ,rho_gen%cloud ,1.,&   ! rain
                       a_vt%cloud, b_vt%cloud, kext_gen%cloud,salb_gen%cloud,asym_gen%cloud,pbck_gen%cloud, vdop_gen%cloud )
     else

      ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud -> RAMS cloud1
                 cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
      call mie_rams_micro('cloud1',freq,tavg,q%cloud,re%cloud,gnu_rams%cloud1,cfvt,pwvt, &
                     cfmas, pwmas , 1.0,  kext_gen%cloud,salb_gen%cloud,asym_gen%cloud,pbck_gen%cloud,vdop_gen%cloud )
     endif

     call mie_gen_micro('qr',freq,tavg,mu_gen%rain,q%rain,re%rain ,rho_gen%rain ,1.,&   ! rain
                       a_vt%rain, b_vt%rain, kext_gen%rain,salb_gen%rain,asym_gen%rain,pbck_gen%rain, vdop_gen%rain )

     if( trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH' .or. &
         trim(cloud_microphysics) == 'THOM' ) then !morrison is twomoment
         call mie_gen_micro('qi',freq,tavg,mu_gen%ice, q%ice,re%ice ,rho_gen%ice ,0.,&   ! ice
                            a_vt%ice, b_vt%ice, kext_gen%ice,salb_gen%ice,asym_gen%ice,pbck_gen%ice, vdop_gen%ice )
     else
         ispc = 3  ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1
                    cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
         call mie_rams_micro('ice1', freq, tavg, q%ice, re%ice , gnu_rams%ice1, cfvt,pwvt, &  !correction...here..
                             cfmas, pwmas , 0.0,  kext_gen%ice , salb_gen%ice , asym_gen%ice , pbck_gen%ice,&
                             vdop_gen%ice )
     endif


     call mie_gen_micro('qs',freq,tavg,mu_gen%snow, q%snow,re%snow,rho_gen%snow,fmelt, & !snow
                        a_vt%snow, b_vt%snow,kext_gen%snow,salb_gen%snow,asym_gen%snow,pbck_gen%snow,vdop_gen%snow )

     call mie_gen_micro('qg',freq,tavg,mu_gen%graupel, q%graupel,re%graupel, rho_gen%graupel, 0., & !graupel
                        a_vt%graupel, b_vt%graupel,kext_gen%graupel,salb_gen%graupel,asym_gen%graupel,pbck_gen%graupel,&
                        vdop_gen%graupel  )

     call mie_gen_micro('qh', freq, tavg, mu_gen%hail, q%hail, re%hail , rho_gen%hail , 0.,& !hail
                        a_vt%hail, b_vt%hail, kext_gen%hail,salb_gen%hail,asym_gen%hail,pbck_gen%hail,vdop_gen%hail )

 endif 

 return
 end subroutine opt_gen_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine opt_rams_micro( nt, ns, nf, freq, tavg, fmelt, q, re, &
                            kext_rams, salb_rams, asym_rams ,pbck_rams,vdop_rams)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation techniqu.
!  This routine is for RAMS microphysics.
!    
! History: 
!  10/2017  Toshi Matsui@NASA GSFC ; New LUT array. 
!  06/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 integer,intent(in):: nt ! Temperature index--1 to 31 T below 0 deg.C
 integer,intent(in):: ns ! RH index--79 to 100 based on integer RH
 integer,intent(in)::nf  !frequency id
 real(sdsu_fps),intent(in)::tavg  ! layer temperature
 real(sdsu_fps),intent(in)::freq  ! frequency [GHz]
 real(sdsu_fps),intent(in)::fmelt ! meltfraction 

 type ( particle_rams ),intent(in) :: & ! particle_rams is defined in module_simulater
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 type ( particle_rams ),intent(out) :: & ! particle_rams is defined in module_simulater
     kext_rams,  & ! extinction coefficient [km-1]
     salb_rams,  & ! single scattering albedo [-]
     asym_rams,  & ! asymetry parameter [-]
     pbck_rams,  & ! backscatter
     vdop_rams     ! Doppler velocity [m/s]

 integer :: ispc  !RAMS spicies index
 real(sdsu_fps) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(sdsu_fps) :: cfvt, pwvt ! alpha_vt & beta_vt


 logical,parameter :: new_method = .true.

 integer :: i,j
 integer :: idx_spc, idx_temp, i_i, i_w
 real(sdsu_fps) :: wgt1, wgt2
 real(sdsu_fps) :: re_in, q_in
 real(sdsu_fps) :: kext_q_unit , salb , asym, pbck, vdop
 real(sdsu_fps) :: tair  ! air temperature [K]
 real(sdsu_fps) :: wgt1_temp , wgt2_temp
 real(sdsu_fps) :: wgt1_wtemp, wgt2_wtemp
 real(sdsu_fps) :: wgt1_itemp, wgt2_itemp
 real(sdsu_fps) :: wgt1_re, wgt2_re


 !
 ! initialize
 !
 call init_microphysics(kext_rams,'zero')
 call init_microphysics(salb_rams,'zero')
 call init_microphysics(asym_rams,'zero')
 call init_microphysics(pbck_rams,'zero')
 call init_microphysics(vdop_rams,'zero')


 LUT_RAMS: if ( lut_micro ) then


!
! air temperature weights (liquid)
!
 tair = tavg
 if (tavg  < pts_temp_w(1))           tair = pts_temp_w(1)
 if (tavg  > pts_temp_w(mxpts_temp))   tair = pts_temp_w(mxpts_temp)

 do i = 1, mxpts_temp-1
    if( tair >= pts_temp_w(i) .and. tair <= pts_temp_w(i+1) ) then
        wgt1_wtemp = ( pts_temp_w(i+1) - tair ) / ( pts_temp_w(i+1) - pts_temp_w(i) )
        wgt2_wtemp = 1.-wgt1_wtemp
        i_w = i
        if( wgt2_wtemp < 0. .or. wgt2_wtemp > 1.0 ) then
            print*, 'MSG opt_gen_micro strange wgt2_wtemp',wgt2_wtemp
            stop
        endif
        exit
    endif
 enddo

!
! air temperature weights (ice)
!
 tair = tavg
 if (tavg  < pts_temp_i(1))           tair = pts_temp_i(1)
 if (tavg  > pts_temp_i(mxpts_temp))   tair = pts_temp_i(mxpts_temp)

 do i = 1, mxpts_temp-1
    if( tair >= pts_temp_i(i) .and. tair <= pts_temp_i(i+1) ) then
        wgt1_itemp = ( pts_temp_i(i+1) - tair ) / ( pts_temp_i(i+1) - pts_temp_i(i) )
        wgt2_itemp = 1.-wgt1_itemp
        i_i = i
        if( wgt2_itemp < 0. .or. wgt2_itemp > 1.0 ) then
            print*, 'MSG opt_gen_micro strange wgt2_itemp',wgt2_itemp
            stop
        endif
        exit
    endif
 enddo


!
! cloud1 species
!
 call interp_lut_array(nf, 1, q%cloud1, re%cloud1, wgt1_wtemp, wgt2_wtemp, i_w, &
      kext_rams%cloud1, salb_rams%cloud1, asym_rams%cloud1, pbck_rams%cloud1, vdop_rams%cloud1 )

!
! cloud2 species
!
 call interp_lut_array(nf, 2, q%cloud2, re%cloud2, wgt1_wtemp, wgt2_wtemp, i_w, &
      kext_rams%cloud2, salb_rams%cloud2, asym_rams%cloud2, pbck_rams%cloud2, vdop_rams%cloud2 )

!
! rain species
!
 call interp_lut_array(nf, 3, q%rain, re%rain, wgt1_wtemp, wgt2_wtemp, i_w, &
      kext_rams%rain, salb_rams%rain, asym_rams%rain, pbck_rams%rain, vdop_rams%rain )

!
! ice1 species
!
 ispc = rams_jhabtab(nt,ns,1)
 call interp_lut_array(nf, ispc, q%ice1, re%ice1, wgt1_itemp, wgt2_itemp, i_i, &
      kext_rams%ice1, salb_rams%ice1, asym_rams%ice1, pbck_rams%ice1, vdop_rams%ice1 )

!
! ice2 species
!
 ispc = rams_jhabtab(nt,ns,2)
 call interp_lut_array(nf, ispc, q%ice2, re%ice2, wgt1_itemp, wgt2_itemp, i_i, &
      kext_rams%ice2, salb_rams%ice2, asym_rams%ice2, pbck_rams%ice2, vdop_rams%ice2 )

!
! snow species
!
 if(fmelt == 0.e0) then ! non-melting snow
   call interp_lut_array(nf, 14, q%snow, re%snow, wgt1_itemp, wgt2_itemp, i_i, &
         kext_rams%snow, salb_rams%snow, asym_rams%snow, pbck_rams%snow, vdop_rams%snow )
 else
    ispc = 5 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! snow aggregate
               cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
    call mie_rams_micro('snow', freq, tavg, q%snow, re%snow , gnu_rams%snow,cfvt,pwvt, &
                        cfmas, pwmas , fmelt,  kext_rams%snow, salb_rams%snow, &
                        asym_rams%snow, pbck_rams%snow, vdop_rams%snow )
 endif

!
! graupel species
!
 call interp_lut_array(nf, 15, q%graupel, re%graupel, wgt1_itemp, wgt2_itemp, i_i, &
      kext_rams%graupel, salb_rams%graupel, asym_rams%graupel, pbck_rams%graupel, vdop_rams%graupel )

!
! hail species
!
 call interp_lut_array(nf, 16, q%hail, re%hail, wgt1_itemp, wgt2_itemp, i_i, &
      kext_rams%hail, salb_rams%hail, asym_rams%hail, pbck_rams%hail, vdop_rams%hail )



 else  !no use of mie LUT (SLOW: compute extinction by integrating hydrometeor size spectra.)

   ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud1
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call mie_rams_micro('cloud1', freq, tavg, q%cloud1, re%cloud1 , gnu_rams%cloud1 ,cfvt,pwvt, &
                       cfmas, pwmas , 1.0,  kext_rams%cloud1, salb_rams%cloud1, asym_rams%cloud1, pbck_rams%cloud1,&
                       vdop_rams%cloud1 )

   ispc = 16; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud2
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call mie_rams_micro('cloud2', freq, tavg, q%cloud2, re%cloud2 , gnu_rams%cloud2 ,cfvt,pwvt, &
                    cfmas, pwmas , 1.0,  kext_rams%cloud2, salb_rams%cloud2, asym_rams%cloud2, pbck_rams%cloud2,&
                    vdop_rams%cloud2 )

   ispc = 2 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! rain 
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call mie_rams_micro('rain', freq, tavg, q%rain, re%rain , gnu_rams%rain, cfvt,pwvt,&
                    cfmas, pwmas , 1.0,  kext_rams%rain, salb_rams%rain, asym_rams%rain, pbck_rams%rain,vdop_rams%rain )

   ispc = rams_jhabtab(nt,ns,1) ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1
                                 cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call mie_rams_micro('ice1', freq, tavg, q%ice1, re%ice1 , gnu_rams%ice1 ,cfvt,pwvt, &
                     cfmas, pwmas , 0.0,  kext_rams%ice1 , salb_rams%ice1 , asym_rams%ice1 , pbck_rams%ice1,&
                      vdop_rams%ice1 )

   ispc = rams_jhabtab(nt,ns,2) ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2
                                 cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call mie_rams_micro('ice2', freq, tavg, q%ice2, re%ice2 , gnu_rams%ice2 ,cfvt,pwvt, &
                    cfmas, pwmas , 0.0,  kext_rams%ice2 , salb_rams%ice2 , asym_rams%ice2 , pbck_rams%ice2,vdop_rams%ice2 )
 
   ispc = 5 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! snow aggregate
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call mie_rams_micro('snow', freq, tavg, q%snow, re%snow , gnu_rams%snow,cfvt,pwvt, &
                    cfmas, pwmas , fmelt,  kext_rams%snow, salb_rams%snow, asym_rams%snow, pbck_rams%snow,vdop_rams%snow )
  
   ispc = 6 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! graupel
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call mie_rams_micro('graupel',freq,tavg,q%graupel,re%graupel,gnu_rams%graupel,cfvt,pwvt, &
                    cfmas, pwmas , 0.0,  kext_rams%graupel, salb_rams%graupel, asym_rams%graupel, pbck_rams%graupel,&
                    vdop_rams%graupel )

   ispc = 7 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! hail
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call mie_rams_micro( 'hail',freq, tavg, q%hail, re%hail, gnu_rams%hail ,cfvt,pwvt, &
                      cfmas, pwmas , 0.0,  kext_rams%hail ,salb_rams%hail ,asym_rams%hail, pbck_rams%hail, vdop_rams%hail )

   endif LUT_RAMS

 return
 end subroutine  opt_rams_micro


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine opt_sbm_micro( freq, tavg, q, n, fmelt, frime, fmelt_g, fmelt_h, &
                           kext_sbm, salb_sbm, asym_sbm ,pbck_sbm, vdop_sbm)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation techniqu.
!  This routine is for HUCM SBM microphysics.
!    
! History: 
!  06/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in)::tavg  ! layer temperature [K]
 real(sdsu_fps),intent(in)::freq  ! frequency [GHz]

 type ( particle_sbm ),intent(in) :: & ! particle_sbm is defined in module_simulater
     q           , & ! bulk mixing ratio [g/m3] 
     n   (1:nbin)    ! particle # concentration [1/m4]
 real(sdsu_fps),intent(in) :: fmelt(1:nbin)    ! mass melting fraction of snow [-]
 real(sdsu_fps),intent(in) :: frime(1:nbin)    ! mass rime fraction of snow [-]
 real(sdsu_fps),intent(in) :: fmelt_g(1:nbin)  ! mass melting fraction of graupel[-]
 real(sdsu_fps),intent(in) :: fmelt_h(1:nbin)  ! mass melting fraction of hail [-]

 type ( particle_sbm ),intent(out) :: & ! particle_sbm is defined in module_simulater
     kext_sbm,  & ! extinction coefficient [km-1]
     salb_sbm,  & ! single scattering albedo [-]
     asym_sbm,  & ! asymetry parameter [-]
     pbck_sbm,  & ! backscatter
     vdop_sbm     ! Doppler velocity  [m/s]

  real(sdsu_fps) :: frime_fake(1:nbin)  !fake 

!
! initialize as zero
!
 call init_microphysics(kext_sbm,'zero')
 call init_microphysics(salb_sbm,'zero')
 call init_microphysics(asym_sbm,'zero')
 call init_microphysics(pbck_sbm,'zero')
 call init_microphysics(vdop_sbm,'zero')


!
! integrate single-scattering property by size integration
!
 frime_fake(1:nbin) = 0. !fake 

 call mie_sbm_micro('liq    ',freq,tavg,den_sbm(1:nbin)%liq,rad_sbm(1:nbin)%liq,drad_sbm(1:nbin)%liq,&
                     vt_sbm(1:nbin)%liq, q%liq,n(1:nbin)%liq,fmelt(1:nbin)      ,frime_fake(1:nbin), &
                     kext_sbm%liq, salb_sbm%liq, asym_sbm%liq, pbck_sbm%liq ,vdop_sbm%liq ) 

 call mie_sbm_micro('ice_col',freq,tavg,den_sbm(1:nbin)%ice_col,rad_sbm(1:nbin)%ice_col,drad_sbm(1:nbin)%ice_col,&
                     vt_sbm(1:nbin)%ice_col, q%ice_col,n(1:nbin)%ice_col,fmelt(1:nbin)      ,frime_fake(1:nbin),&
                     kext_sbm%ice_col, salb_sbm%ice_col, asym_sbm%ice_col, pbck_sbm%ice_col,vdop_sbm%ice_col  ) 

 call mie_sbm_micro('ice_pla',freq,tavg,den_sbm(1:nbin)%ice_pla,rad_sbm(1:nbin)%ice_pla,drad_sbm(1:nbin)%ice_pla,&
                     vt_sbm(1:nbin)%ice_pla, q%ice_pla,n(1:nbin)%ice_pla,fmelt(1:nbin)     ,frime_fake(1:nbin),&
                     kext_sbm%ice_pla, salb_sbm%ice_pla, asym_sbm%ice_pla, pbck_sbm%ice_pla,vdop_sbm%ice_pla  ) 

 call mie_sbm_micro('ice_den',freq,tavg,den_sbm(1:nbin)%ice_den,rad_sbm(1:nbin)%ice_den,drad_sbm(1:nbin)%ice_den,&
                     vt_sbm(1:nbin)%ice_den, q%ice_den,n(1:nbin)%ice_den,fmelt(1:nbin)      ,frime(1:nbin),&
                     kext_sbm%ice_den, salb_sbm%ice_den, asym_sbm%ice_den, pbck_sbm%ice_den,vdop_sbm%ice_den  ) 

 call mie_sbm_micro('snow   ',freq,tavg,den_sbm(1:nbin)%snow,rad_sbm(1:nbin)%snow,drad_sbm(1:nbin)%snow,&
                     vt_sbm(1:nbin)%snow, q%snow,n(1:nbin)%snow,fmelt(1:nbin)      ,frime(1:nbin),&
                     kext_sbm%snow, salb_sbm%snow, asym_sbm%snow, pbck_sbm%snow, vdop_sbm%snow  ) 

 call mie_sbm_micro('graupel',freq,tavg,den_sbm(1:nbin)%graupel,rad_sbm(1:nbin)%graupel,drad_sbm(1:nbin)%graupel,&
                     vt_sbm(1:nbin)%graupel, q%graupel,n(1:nbin)%graupel,fmelt_g(1:nbin)     ,frime_fake(1:nbin),&
                     kext_sbm%graupel, salb_sbm%graupel, asym_sbm%graupel, pbck_sbm%graupel, vdop_sbm%graupel ) ! graupel

 call mie_sbm_micro('hail   ',freq,tavg,den_sbm(1:nbin)%hail,rad_sbm(1:nbin)%hail,drad_sbm(1:nbin)%hail,&
                     vt_sbm(1:nbin)%hail, q%hail,n(1:nbin)%hail,fmelt_h(1:nbin)    ,frime_fake(1:nbin),&
                     kext_sbm%hail, salb_sbm%hail, asym_sbm%hail, pbck_sbm%hail, vdop_sbm%hail  ) ! hail

 return
 end subroutine opt_sbm_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine sback_hid_sbm( freq, tavg, q, n, fmelt, frime, fmelt_g, fmelt_h, &
                           sback_hid )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute bask-scattering properties from spectra-bin microphysics with respect to 
!  CSU HID radar retrieval species.
!  This routine is for HUCM SBM microphysics.
!    
! History: 
!  12/2013  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in)::tavg  ! layer temperature [K]
 real(sdsu_fps),intent(in)::freq  ! frequency [GHz]

 type ( particle_sbm ),intent(in) :: & ! particle_sbm is defined in module_simulater
     q           , & ! bulk mixing ratio [g/m3] 
     n   (1:nbin)    ! particle # concentration [1/m4]
 real(sdsu_fps),intent(in) :: fmelt(1:nbin)    ! mass melting fraction of snow [-]
 real(sdsu_fps),intent(in) :: frime(1:nbin)    ! mass rime fraction of snow [-]
 real(sdsu_fps),intent(in) :: fmelt_g(1:nbin)  ! mass melting fraction of graupel[-]
 real(sdsu_fps),intent(in) :: fmelt_h(1:nbin)  ! mass melting fraction of hail [-]
 type(particle_hid),intent(out) :: sback_hid   ! backscatter
 
 real(sdsu_fps) :: sback
 real(sdsu_fps) :: & 
     kext,  & ! extinction coefficient [km-1]
     salb,  & ! single scattering albedo [-]
     asym,  & ! asymetry parameter [-]
     pbck,  & ! backscatter
     vdop     ! Doppler velocity  [m/s]

 real(sdsu_fps) :: frime_fake(1:nbin)  !fake 

 real(sdsu_fps) :: n_temp(1:nbin)    ! temporally particle # concentration [1/m4]

 real(sdsu_fps) ::  q_tot, q_melt, q_rime  !total, melt, rimed tota mass concentration. 

 real(sdsu_fps) :: hid  !hydrometero identification 
 real(sdsu_fps) :: hid_confidence ! confidence of level [%]


!
! initialize
!
 call init_microphysics(sback_hid,'zero')
 
!
! integrate single-scattering property by size integration
!
 frime_fake(1:nbin) = 0. !fake 

!
! DZ: drizzle (diameter 0.3 ~ 0.55mm of sbm liq: nbin 20 ~22) 
!
 n_temp(1:nbin) = 0.e0
 n_temp(20:22)  = n(20:22)%liq
 call mie_sbm_micro('liq    ',freq,tavg,den_sbm(1:nbin)%liq,rad_sbm(1:nbin)%liq,drad_sbm(1:nbin)%liq,&
                     vt_sbm(1:nbin)%liq, q%liq, n_temp(1:nbin), fmelt(1:nbin), frime_fake(1:nbin), &
                     kext, salb, asym, pbck ,vdop ) 

 sback_hid%dz = kext * salb * pbck  !backscatter of DZ component
 
!
! RN: rain (0.55mm~ 5mm of sbm liq: nbin 23~31)
!
 n_temp(1:nbin) = 0.e0
 n_temp(23:31)  = n(23:31)%liq
 call mie_sbm_micro('liq    ',freq,tavg,den_sbm(1:nbin)%liq,rad_sbm(1:nbin)%liq,drad_sbm(1:nbin)%liq,&
                     vt_sbm(1:nbin)%liq, q%liq,n_temp(1:nbin),fmelt(1:nbin)      ,frime_fake(1:nbin), &
                     kext, salb, asym, pbck ,vdop )

 sback_hid%rn = kext * salb * pbck  ! backscatter of RN component

!
! VI: vertically oriented ice crystal (all diameter of sbm ice column)
!

!
! CR: ice crystals (all diameter of sbm ice plate and sbm dendrite )
!
 sback = 0.
 n_temp(1:nbin) = 0.e0
 n_temp(1:nbin)  = n(1:nbin)%ice_col
 call mie_sbm_micro('ice_col',freq,tavg,den_sbm(1:nbin)%ice_col,rad_sbm(1:nbin)%ice_col,drad_sbm(1:nbin)%ice_col,&
                     vt_sbm(1:nbin)%ice_col, q%ice_col,n_temp(1:nbin),fmelt(1:nbin)      ,frime_fake(1:nbin),&
                     kext, salb, asym, pbck ,vdop  ) 

! sback_hid%vi = kext * salb * pbck   ! backscatter of VI component
 sback = kext * salb * pbck

 n_temp(1:nbin) = 0.e0
 n_temp(1:nbin)  = n(1:nbin)%ice_pla
 call mie_sbm_micro('ice_pla',freq,tavg,den_sbm(1:nbin)%ice_pla,rad_sbm(1:nbin)%ice_pla,drad_sbm(1:nbin)%ice_pla,&
                     vt_sbm(1:nbin)%ice_pla, q%ice_pla,n_temp(1:nbin),fmelt(1:nbin)     ,frime_fake(1:nbin),&
                     kext, salb, asym, pbck ,vdop ) 
 sback = kext * salb * pbck + sback

 n_temp(1:nbin) = 0.e0
 n_temp(1:nbin)  = n(1:nbin)%ice_den
 call mie_sbm_micro('ice_den',freq,tavg,den_sbm(1:nbin)%ice_den,rad_sbm(1:nbin)%ice_den,drad_sbm(1:nbin)%ice_den,&
                     vt_sbm(1:nbin)%ice_den, q%ice_den,n_temp(1:nbin),fmelt(1:nbin)      ,frime(1:nbin),&
                     kext, salb, asym, pbck ,vdop  ) 
 sback = kext * salb * pbck + sback

 n_temp(1:nbin) = 0.e0
 n_temp(1:22)  = n(1:22)%snow
 call mie_sbm_micro('snow   ',freq,tavg,den_sbm(1:nbin)%snow,rad_sbm(1:nbin)%snow,drad_sbm(1:nbin)%snow,&
                     vt_sbm(1:nbin)%snow, q%snow,n_temp(1:nbin),fmelt(1:nbin)      ,frime(1:nbin),&
                     kext, salb, asym, pbck, vdop  )
 sback = kext * salb * pbck + sback

 n_temp(1:nbin) = 0.e0
 n_temp(1:24)  = n(1:24)%graupel
 call mie_sbm_micro('graupel',freq,tavg,den_sbm(1:nbin)%graupel,rad_sbm(1:nbin)%graupel,drad_sbm(1:nbin)%graupel,&
                     vt_sbm(1:nbin)%graupel, q%graupel,n_temp(1:nbin),fmelt_g(1:nbin)     ,frime_fake(1:nbin),&
                     kext, salb, asym, pbck ,vdop ) ! graupel
 sback = kext * salb * pbck + sback

 n_temp(1:nbin) = 0.e0
 n_temp(1:24)  = n(1:24)%hail
 call mie_sbm_micro('hail   ',freq,tavg,den_sbm(1:nbin)%hail,rad_sbm(1:nbin)%hail,drad_sbm(1:nbin)%hail,&
                     vt_sbm(1:nbin)%hail, q%hail,n_temp(1:nbin),fmelt_h(1:nbin)    ,frime_fake(1:nbin),&
                     kext, salb, asym, pbck ,vdop  ) ! hail

 sback_hid%CR = kext * salb * pbck + sback   ! backscatter of IC component


!
! DS: aggregate (1mm ~ 12mm of sbm aggregate: 23 ~31bin), 
! low density graupel (with riming fraction), or wet snow (with melting fraction)
!
 n_temp(1:nbin) = 0.e0
 n_temp(23:nbin)  = n(23:nbin)%snow  ! 31 --> nbin
 call mie_sbm_micro('snow   ',freq,tavg,den_sbm(1:nbin)%snow,rad_sbm(1:nbin)%snow,drad_sbm(1:nbin)%snow,&
                     vt_sbm(1:nbin)%snow, q%snow,n_temp(1:nbin),fmelt(1:nbin)      ,frime(1:nbin),&
                     kext, salb, asym, pbck, vdop  ) 
!
! Check total riming mass and melting fraction
!
 q_tot  = sum(              x_sbm(23:31)%snow*n_temp(23:31)*drad_sbm(23:31)%snow )
 q_melt = sum( fmelt(23:31)*x_sbm(23:31)%snow*n_temp(23:31)*drad_sbm(23:31)%snow )
 q_rime = sum( frime(23:31)*x_sbm(23:31)%snow*n_temp(23:31)*drad_sbm(23:31)%snow )

 if( q_tot /= 0.e0 ) then  !use 20% threshold (--> tunable here) 
    if( q_melt / q_tot > 0.2 ) then
        sback_hid%ws = kext * salb * pbck    ! backscatter of WS component
    elseif( q_rime / q_tot > 0.3 ) then
        sback_hid%ldg = kext * salb * pbck   ! backscatter of LDG component 
    else 
        sback_hid%ds = kext * salb * pbck    ! backscatter of DS component
    endif
 endif 

!
! HDG: high density graupel (sbm graupel: 1mm ~ 10mm: 25~33bin) and solid rain (~0.5cm of sbm hail, 25~31bin)
!
 n_temp(1:nbin) = 0.e0
 n_temp(25:33)  = n(25:33)%graupel
 call mie_sbm_micro('graupel',freq,tavg,den_sbm(1:nbin)%graupel,rad_sbm(1:nbin)%graupel,drad_sbm(1:nbin)%graupel,&
                     vt_sbm(1:nbin)%graupel, q%graupel,n_temp(1:nbin),fmelt_g(1:nbin)     ,frime_fake(1:nbin),&
                     kext, salb, asym, pbck ,vdop ) ! graupel
  sback_hid%hdg = kext * salb * pbck  ! backscatter of HDG component


! n_temp(1:nbin) = 0.e0
! n_temp(25:31)  = n(25:31)%hail
! call mie_sbm_micro('hail   ',freq,tavg,den_sbm(1:nbin)%hail,rad_sbm(1:nbin)%hail,drad_sbm(1:nbin)%hail,&
!                     vt_sbm(1:nbin)%hail, q%hail,n_temp(1:nbin),fmelt_h(1:nbin)    ,frime_fake(1:nbin),&
!                     kext, salb, asym, pbck ,vdop  ) ! hail
! sback_hid%hdg = kext * salb * pbck + sback   ! backscatter of HDG component
 

!
! BD: big drops (5mm~1.1cm of large drop or melted hail: 32~35bin  )
!
 sback = 0.e0
 n_temp(1:nbin) = 0.e0
 n_temp(32:nbin)  = n(32:nbin)%liq
 call mie_sbm_micro('liq    ',freq,tavg,den_sbm(1:nbin)%liq,rad_sbm(1:nbin)%liq,drad_sbm(1:nbin)%liq,&
                     vt_sbm(1:nbin)%liq, q%liq,n_temp(1:nbin),fmelt(1:nbin)      ,frime_fake(1:nbin), &
                     kext, salb, asym, pbck ,vdop )
 sback = kext * salb * pbck  ! backscatter of RN component

 n_temp(1:nbin) = 0.e0
 n_temp(32:nbin)  = n(32:nbin)%hail
 call mie_sbm_micro('hail   ',freq,tavg,den_sbm(1:nbin)%hail,rad_sbm(1:nbin)%hail,drad_sbm(1:nbin)%hail,&
                     vt_sbm(1:nbin)%hail, q%hail,n_temp(1:nbin),fmelt_h(1:nbin)    ,frime_fake(1:nbin),&
                     kext, salb, asym, pbck ,vdop  ) ! hail
 q_tot  = sum(                  x_sbm(32:nbin)%hail*n_temp(32:nbin)*drad_sbm(32:nbin)%hail )
 q_melt = sum( fmelt_h(32:nbin)*x_sbm(32:nbin)%hail*n_temp(32:nbin)*drad_sbm(32:nbin)%hail )
 
 if( q_tot > 0.e0) then 
   if( q_melt/q_tot > 0.2 ) then !enought melt
      sback_hid%bd  = kext * salb * pbck + sback   ! backscatter of BD component
   else !no melt (--> solid big hail)
      sback = 0.
      sback = kext * salb * pbck  ! backscatter of HL component
   endif 
 else
    sback_hid%bd  = sback ! large liquid drop only
 endif


!
! HL: hail (0.5cm~5cm: 25~31bin)
!
 n_temp(1:nbin) = 0.e0
 n_temp(25:31)  = n(25:31)%hail
 call mie_sbm_micro('hail   ',freq,tavg,den_sbm(1:nbin)%hail,rad_sbm(1:nbin)%hail,drad_sbm(1:nbin)%hail,&
                     vt_sbm(1:nbin)%hail, q%hail,n_temp(1:nbin),fmelt_h(1:nbin)    ,frime_fake(1:nbin),&
                     kext, salb, asym, pbck ,vdop  ) ! hail
 sback_hid%ha = kext * salb * pbck + sback   ! backscatter of HL component


!
! find id here for diagnostic
!
  call find_hid0d( sback_hid, hid , hid_confidence)

  if( hid == 2.e0 ) then !rain

    !sample PSD

    call sample_psd_sbm( 'add', 'hail', rad_sbm(1:nbin)%hail, drad_sbm(1:nbin)%hail,&
                          q%hail,n(1:nbin)%hail,fmelt_h(1:nbin) )
         
    call sample_psd_sbm( 'add', 'liq', rad_sbm(1:nbin)%liq, drad_sbm(1:nbin)%liq,&
                          q%liq ,n(1:nbin)%liq,fmelt_h(1:nbin) )              


  endif

 return
 end subroutine sback_hid_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine makelut_micro(sensor)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! Construct the look-up tables of Mie coefficients
! for radar and microwave radiometer simulations.
!
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 character(len=5),intent(in) :: sensor  !'radar' or 'micro'
 character(len=30) :: fname_out(100)
 character(len=100):: char_freq
 integer :: max_spc ! maximum particle species
 real(sdsu_fps) :: kext  !extinction [1/km]
 real(sdsu_fps) :: salb  !single scatter albedo [-]
 real(sdsu_fps) :: asym   ! asymetry paramete [-]
 real(sdsu_fps) :: pbck   ! backscatter phase function/(4*pi) [-]
 real(sdsu_fps) :: vdop   ! Doppler velocity [m/s]
 real(sdsu_fps) :: rhoe   ! re-dependent bulk snow density [kg/m3]
 real(sdsu_fps) :: log10_nw  !log10 of Nw (or N0) [log10(1/m4)]
 real(sdsu_fps) :: n0     ! N0 [1/m4]

 real(sdsu_fps) :: frequency
 real(sdsu_fps),parameter :: q_unit = 1.e0 ! unit water content [g/m3]
 integer :: nspec,i,j,nfreq
 logical :: fexist
 integer :: mxfreq ! maximum # of frequency
 integer :: ispc   ! rams particle class
 real(sdsu_fps) :: cfmas, pwmas ! alpha_m & beta_m 
 real(sdsu_fps) :: cfvt, pwvt   ! alpha_m & beta_m 
 integer :: ierr
!
! NetCDF parameters
!
 integer :: ncid, id
 integer :: ref_dimid, temp_dimid, spc_dimid, freq_dimid
 integer, parameter :: ndims4d = 4
 integer :: dimids4d(ndims4d)
 integer :: id_kext, id_salb, id_asym, id_pbck, id_vdop


!
! program start
!

 if(sensor == 'radar') then
    mxfreq = mxfreq_radar
   write(*,*) '- MSG makelut_micro: LUT for Radar - '
 elseif(sensor == 'micro') then
    mxfreq = mxfreq_micro
    write(*,*) '- MSG makelut_micro: LUT for Microwave Radiometer- '
 else
    call stop_sdsu(myrank,'There is no such sensor')
 endif


!
! static indice
!
 select case(trim(type_microphysics))
 case('GEN')    ! GCE (LIN_type) 1moment bulk microphsics 
    max_spc = 6
 case('RAMS')  !GCE RAMS 1- and 2-moment microphysics 
    max_spc = 16
 case('SBM')
    call stop_sdsu(myrank,'Not yet supporoting SBM')
 case default
 end select


!
! allocate array parameters
!
 if( .not. allocated(kext_q_unit_array) ) then
   allocate( kext_q_unit_array ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             salb_array        ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             asym_array        ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             pbck_array        ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             vdop_array        ( mxpts_re, mxpts_temp, max_spc, max_chan ), &
             stat=ierr )
   kext_q_unit_array =  0.
   salb_array = 0.
   asym_array = 0.
   pbck_array = 0.
   vdop_array = 0.
   if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_micro: allocation error -> Terminate program.')
 endif

!
! open netcdf file
!
 if(sensor == 'radar') then
 sdsu_io_file = trim(sdsu_dir_sslut)//'MICRO_CLOUD_'//trim(cloud_microphysics)//'-'//trim(radar_sensor)//'.nc'
 elseif(sensor == 'micro') then
 sdsu_io_file = trim(sdsu_dir_sslut)//'MICRO_CLOUD_'//trim(cloud_microphysics)//'-'//trim(micro_sensor)//'.nc'
 endif

 inquire(file=trim(sdsu_io_file),exist=fexist)
 if( lut_replace ) fexist = .false.
 if (fexist) then !skip
     print*,'MSG makelut_micro; ',trim(sdsu_io_file),' exist.'
     return
 endif

 print*,'MSG makelut_micro; Generating microwave LUT ',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )               !32bit output

!
! Define netcdf dimensions
!
 call check( nf90_def_dim(ncid, 're_max'   , mxpts_re   , ref_dimid ) )
 call check( nf90_def_dim(ncid, 'temp_max' , mxpts_temp , temp_dimid) )
 call check( nf90_def_dim(ncid, 'spc_max'  , max_spc    , spc_dimid ) )
 call check( nf90_def_dim(ncid, 'freq_max' , max_chan   , freq_dimid) )
 dimids4d  = (/ ref_dimid, temp_dimid, spc_dimid, freq_dimid /)

!
! global parameter
!
 if(sensor == 'radar') then
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'Frequency[GHz]', freq_radar(1:mxfreq)  ) )
 elseif(sensor == 'micro') then
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'Frequency[GHz]', freq_micro(1:mxfreq)  ) )
 endif

!
! output parameter definition
!
 call check( nf90_def_var(ncid, 'kext_q_unit', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' extinction coefficient per unit mass conc ') )
 call check( nf90_put_att(ncid, id, 'units', '[km-1 / (g/m3)]') )
 id_kext = id

 call check( nf90_def_var(ncid, 'salb', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' single scattering albedo ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_salb = id

 call check( nf90_def_var(ncid, 'asym', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' asymetry parameter ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_asym = id

 call check( nf90_def_var(ncid, 'pbck', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' backscatter phase function/(4*pi) ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_pbck = id

 call check( nf90_def_var(ncid, 'vdop', NF90_REAL, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' Doppler velocity ') )
 call check( nf90_put_att(ncid, id, 'units', '[cm/s]') )
 id_vdop = id

 call check( nf90_enddef(ncid) )  ! End define mode.

!
! loop start
!
 freq_loop : do nfreq = 1 , mxfreq

   if(sensor == 'radar') then
      frequency = freq_radar(nfreq)
      write(*,*) 'Freq = ',freq_radar(nfreq),' GHz'
   elseif(sensor == 'micro') then
      frequency = freq_micro(nfreq)
      write(*,*) 'Freq = ',freq_micro(nfreq),' GHz'
   endif

!
! species loop 
!
 nspec_loop: do nspec = 1, max_spc


!
! Calculate scattering parameters
! toshii
 temp_loop: do i = 1, mxpts_temp
   re_loop: do j  = 1, mxpts_re
      kext = 0. ; salb = 0. ; asym = 0. ; pbck = 0. ; vdop = 0.

      mic_select2: select case(trim(type_microphysics))

      case('GEN')      ! general bulk microphsics 

           ! LUT as function of T and q (for empirical DSD)
            NSPEC_GEN: if(nspec==1) then

               if( trim(cloud_microphysics) == 'WDM') then 
                call mie_gen_micro('qc',frequency, pts_temp_w(i), mu_gen%cloud, q_unit, pts_re(j), rho_gen%cloud , 1.,& 
                                    a_vt%cloud, b_vt%cloud, kext, salb,asym, pbck, vdop)
               else
                 ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud1
                          cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)    ! 
                 call mie_rams_micro('cloud1',frequency, pts_temp_w(i), q_unit, pts_re(j) , gnu_rams%cloud1 , &
                                      cfvt,pwvt, cfmas, pwmas , 1.0,  kext, salb, asym, pbck, vdop )
               endif

            elseif(nspec==2) then

              call mie_gen_micro('qr',frequency, pts_temp_w(i), mu_gen%rain, q_unit, pts_re(j), rho_gen%rain , 1.,& !rain
                                 a_vt%rain, b_vt%rain, kext, salb,asym, pbck, vdop)

            elseif(nspec==3) then

              if( trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH' .or. &
                  trim(cloud_microphysics) == 'THOM' ) then !  twomoment
                call mie_gen_micro('qi',frequency, pts_temp_w(i), mu_gen%ice, q_unit, pts_re(j), rho_gen%ice , 0.,& !ice
                                   a_vt%ice, b_vt%ice, kext, salb,asym, pbck, vdop)

              else  
                ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice -> RAMS ice1 col
                          cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
                call mie_rams_micro('ice1', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice1 , &
                                    cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

              endif

            elseif(nspec==4) then

              if( trim(cloud_microphysics) == 'GMP4ICE' ) then
               rhoe = 1.e+3 * max(0.05, min(0.9e0, 142.36 * ( pts_re(j)**(-0.95) ) ) ) ! re-dependent bulk snow density [kg/m3]
              elseif( trim(cloud_microphysics) == 'THOM' ) then
               rhoe = min(917. , 0.13 / (pts_re(j)*2.*1.e-6) ) ! Kg/m^3  (917. to be the same as densice)
              else
                rhoe = rho_gen%snow
              endif
              call mie_gen_micro('qs',frequency, pts_temp_i(i), mu_gen%snow, q_unit, pts_re(j), rhoe , 0.,& !snow
                                a_vt%snow, b_vt%snow, kext, salb,asym, pbck, vdop)

!print*,'Snow Mie=',kext, salb,asym
              !
              ! in case of OpenSSP
              !

!              if( ss_opt_micro == 3 ) then
!                 call get_n0( mu_gen%snow, q_unit, rhoe, pts_re(j), n0 )
!                 if( n0 == 0. ) then
!                     log10_nw = 0.  
!                 else 
!                     log10_nw = log10( n0 / 8.e6 )  ! 8.e6 is reference n0
!                 endif

!                 call getsnowp2_sdsu( log10_nw-6.90309, q_unit, frequency, kext, salb, asym )
!print*,'OpenSSP Mie=', kext, salb,asym
!print*,'Re=',pts_re(j),'n0=',n0
!print*,''
!              endif

            elseif(nspec==5) then

              call mie_gen_micro('qg',frequency,pts_temp_i(i),mu_gen%graupel,q_unit,pts_re(j),rho_gen%graupel,0.,& !graupel
                                a_vt%graupel, b_vt%graupel, kext, salb,asym, pbck, vdop)

            elseif(nspec==6) then

              call mie_gen_micro('qh',frequency, pts_temp_i(i), mu_gen%hail, q_unit, pts_re(j), rho_gen%hail , 0.,& !hail
                                 a_vt%hail, b_vt%hail, kext, salb,asym, pbck, vdop)

            else
               call stop_sdsu(myrank,'MSG makelut_micro: there is no such nspec')
            endif NSPEC_GEN

       case('RAMS')  !GCE RAMS 1- & 2-moment microphysics 

            NSPEC_RAMS: if(nspec==1) then
              ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud1
                         cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('cloud1',frequency, pts_temp_w(i), q_unit, pts_re(j) , gnu_rams%cloud1 , &
                                  cfvt,pwvt,cfmas, pwmas , 1.0,  kext, salb, asym, pbck, vdop )
 
            elseif(nspec==2) then
              ispc = 16; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud2
                         cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('cloud2', frequency, pts_temp_w(i), q_unit, pts_re(j) , gnu_rams%cloud2 , &
                                  cfvt,pwvt,cfmas, pwmas , 1.0,  kext, salb, asym, pbck, vdop )

            elseif(nspec==3) then
              ispc = 2 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! rain 
                         cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('rain', frequency, pts_temp_w(i), q_unit, pts_re(j) , gnu_rams%rain, &
                                  cfvt,pwvt,cfmas, pwmas , 1.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==4) then
              ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 col
                        cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice1_col', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice1 , &
                                   cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==5) then
              ispc = 8 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 hex
                        cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice1_hex', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice1 , &
                                  cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==6) then
              ispc = 9 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 den
                        cfvt  = rams_dstprms(4,ispc) ;pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice1_den', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice1 , &
                                  cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==7) then
              ispc = 10 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 ndl
                         cfvt  = rams_dstprms(4,ispc) ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice1_ndl', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice1 , &
                                   cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==8) then
              ispc = 11 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 ros
                         cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice1_ros', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice1 , &
                                   cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==9) then
              ispc = 4 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 col
                        cfvt  = rams_dstprms(4,ispc) ;pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice2_col', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice2 , &
                                   cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==10) then
              ispc = 12 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 hex
                         cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice2_hex', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice2 , &
                                   cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==11) then
              ispc = 13 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 den
                         cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice2_den', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice2 , &
                                  cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==12) then
              ispc = 14 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 ndl
                         cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice2_ndl', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice2 , &
                                   cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==13) then
              ispc = 15 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 ros
                         cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('ice2_ros', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%ice2 , &
                                  cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==14) then 
              ispc = 5 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! snow aggregate
                         cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
              call mie_rams_micro('snow', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%snow, &   
                                  cfvt,pwvt,cfmas, pwmas ,  0.0, kext, salb, asym, pbck, vdop )
 
            elseif(nspec==15) then 
               ispc = 6 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! graupel
                          cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
               call mie_rams_micro('graupel', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%graupel, &   
                                   cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )

            elseif(nspec==16) then
               ispc = 7 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! hail
                          cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
               call mie_rams_micro('hail', frequency, pts_temp_i(i), q_unit, pts_re(j) , gnu_rams%hail , &
                                    cfvt,pwvt,cfmas, pwmas , 0.0, kext, salb, asym, pbck, vdop )
            
            else
               call stop_sdsu(myrank,'MSG makelut_micro: there is no such nspec')
            endif NSPEC_RAMS

       case('SBM')

              call stop_sdsu(myrank,'MSG makelut_micro: do not use lut option for SBM')

       case default 

       end select mic_select2

       !
       ! transfer to giant array 
       !
       kext_q_unit_array ( j,i,nspec,nfreq ) = kext
       salb_array        ( j,i,nspec,nfreq ) = salb
       asym_array        ( j,i,nspec,nfreq ) = asym
       pbck_array        ( j,i,nspec,nfreq ) = pbck
       vdop_array        ( j,i,nspec,nfreq ) = vdop

      enddo re_loop
   enddo temp_loop

 enddo nspec_loop

 enddo freq_loop



!
! write output in NetCDF
!
 call check( nf90_put_var(ncid, id_kext , kext_q_unit_array ) )
 call check( nf90_put_var(ncid, id_salb , salb_array        ) )
 call check( nf90_put_var(ncid, id_asym , asym_array        ) )
 call check( nf90_put_var(ncid, id_pbck , pbck_array        ) )
 call check( nf90_put_var(ncid, id_vdop , vdop_array        ) )

 call check( nf90_close(ncid) )  ! close nc file


 return
 end subroutine makelut_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine mie_gen_micro(spc, freq, temp, mu, wc, re, dens, frac_liq,a_vt,b_vt, &
                           ksca, asca, gsca, pbck, vdop)
  implicit none
!--------------------------------------------------------------------------------
! Purpose : Compute the extinction, absorption, asymmetry parameter and
!           backscatter for a given water content of condensates in [g/m^3], and
!           a particle size distribution n(D) with intercept n0.
!           N(D) = No * D^mu * exp( -lambda * D )
!  
!  08/2016  Toshi Matsui @ NASA GSFC : Clean up by adding effective_permittivity routine, and also add 
!                                      snow compaction for melting calculaiton consistent to RAMS and SBM.  
!  07/2015  Toshi Matsui@NASA GSFC ; Size ranges changed between precip and non-precip species.
!  01/2014  Toshi Matsui@NASA GSFC ; add Tmatrix option
!  01/2011  Toshi Matsui@NASA GSFC ; add SCATDB option
!  03/2010  Toshi Matsui@NASA GSFC ; Additional options for effective refractive index for solid particles.
!
!--------------------------------------------------------------------------------
 character(len=*),intent(in) :: spc !species character
 real(sdsu_fps),intent(in) :: freq  ! frequency of radiation [GHz]
 real(sdsu_fps),intent(in) :: temp  ! temperature of particles [K]
 real(sdsu_fps),intent(in) :: mu    ! dispersion parameter [-]
 real(sdsu_fps),intent(in) :: wc    ! water content of condensate [g/m**3]
 real(sdsu_fps),intent(in) :: re    ! drop effective radius [micron] 
 real(sdsu_fps),intent(in) :: dens  ! bulk effective density of hydrometero [kg/m**3]
 real(sdsu_fps),intent(in) :: frac_liq ! fracion of liquid (1=rain or 0=ice condensates, 0~1=melting ice)
 real(sdsu_fps),intent(in) :: a_vt     ! alpha of Vt  [-]
 real(sdsu_fps),intent(in) :: b_vt     ! beta of Vt   [-]
 real(sdsu_fps),intent(out) :: ksca  ! extinction coefficient [1/km]
 real(sdsu_fps),intent(out) :: asca  ! single-scatter albedo [-]
 real(sdsu_fps),intent(out) :: gsca  ! asymmetry factor [-]
 real(sdsu_fps),intent(out) :: pbck  ! backscatter phase function/(4*pi) [-]
 real(sdsu_fps),intent(out) :: vdop  ! Doppler velocity [cm/s]

!--------Local Parameters
 integer :: i    !looping 
 integer :: imax !max for looping
 real(sdsu_fps) :: rhoe ! bulk effective density of hydrometero [kg/m**3]

 real(sdsu_fps) :: n0   ! intercept for the exponetial DSD [1/m**4]
 real(sdsu_fps) :: wave
 real(sdsu_fps) :: rad  ! radius of particle [mm]
 real(sdsu_fps) :: lam  ! the slope of the distribution [1/m]
 real(sdsu_fps) :: num  ! a particle number density per radius increment [1/m**4]
 real(sdsu_fps) :: faa  ! fraction of air in ice 
 real(sdsu_fps) :: eicere
 real(sdsu_fps) :: eiceim
 real(sdsu_fps) :: qext
 real(sdsu_fps) :: qsca
 real(sdsu_fps) :: asym
 real(sdsu_fps) :: qbsca
 real(sdsu_fps) :: bext
 real(sdsu_fps) :: bsca
 real(sdsu_fps) :: bsym
 real(sdsu_fps) :: bq11
 real(sdsu_fps) :: vt_wsum  !weighting sum of terminal velocity 
 real(sdsu_fps) :: vt  !terminal velocity [m/s]
 real(sdsu_fps) :: fmelt


 real(sdsu_fps) :: eimag, ereal

 real(sdsu_fps) :: gfac1, gfac2, gamfac  !gamma parameter
 real(sdsu_fps) :: &
   factor     ,& ! 
   d_mm       ,& ! particle diameter [mm]
   d_m        ,& ! particle diameter [m]
   d_increment,& ! diameter increment over size loops [mm]
   d_min      ,& ! minimum diameter of particle [mm]
   d_max         ! maximum diameter of particle [mm]

 real(sdsu_fpd) :: &
      xd     ,&       ! = dble(2.*pi*rad/wave)
      freqd  ,&
      freqhzd,& ! frequency in double precision
      tempd  ,&       ! temperaure in doubpl precision
      sald   ,&
      eicered,& 
      ewatred,&
      eiceimd,&
      ewatimd,&
      qscad  ,&
      qextd  ,&
      asymd  ,&
      qbscad 

    
 complex *8 eice,ewat,eair !ice, water, and air permittivity
 complex *8 ei     !effective ice-air permittivity
 complex *8 emelt,emelt_wi,emelt_iw, emelt_em !effective water-ice-air permittivity

 complex*16 ewatd
 complex*16 eid, emeltd
 complex*16 crefd, crefd_liq, crefd_ice

 real(sdsu_fps),parameter :: dr = 0.05    ! increment bin of radius [mm]
 real(sdsu_fps),parameter :: densice = 0.917e+3  !solid ice density [kg/m**3]
 real(sdsu_fps),parameter :: densliq = 1.0e+3    !liquid density [kg/m**3]
 real(sdsu_fps),parameter :: wc_unit = 1.0  !unit water content 1 [g/m3]
 integer,save :: nshp !SCATDB shape #
 real(sdsu_fps) :: qext_ns, qsca_ns, asym_ns, qbsca_ns ! non-spherical single scattering
 real(sdsu_fps) :: d_mm_solid  !solid ice diameter [mm]
 real(sdsu_fps) :: rad_solid_micron  !solid ice radius [micron]
 real(sdsu_fps) :: area   !area of fully particle [mm2]
 logical :: success ! status 
 real(sdsu_fps) :: dens_spectrum  ! spectral effective density for each size [kg/m3]
 real(sdsu_fps),allocatable :: thom_snow_num(:)  !size distributions of thompson scheme [1/m4]
 real(sdsu_fps) :: rad_new, fmelt_new
! real(sdsu_fps) :: dcreal, dcimag  !real / imageneary component of dielectric function


!
! Begin by checking if hydrometeors of this species are present.
! If not, set scattering parameters to zero and return.
!
  if(wc .lt. q_min_micro .or. re == 0. ) then
     ksca=0. ; asca=0. ; gsca=0. ; pbck=0. ; vdop=0.
     return
  endif

!
! SCATDB option for ice crystal
!
 if( trim(spc) == 'qi' .and. ss_opt_micro == 2 )  then
    nshp = scatdb_ice_type  !determined by configuration file
    call scatdb_prep(freq, temp, nshp)
 endif

!
! preparting ice Tmatrix option  (so far does not apply to melting species or liquid species)
!
 select case( trim(spc) )
 case( 'qs','qg','qh' )
   if( ss_opt_micro == 4 ) &
     call  tmatrix_interp( i_can_say=masterproc, proc='prep', freq=freq )
 case default
 end select

!
! Assign some useful constants
!
  wave=const_ghz_mm/freq

!
! If hydrometeors are present, initialize the scattering parameters
!
  bext=0. ; bsca=0. ; bsym=0. ; bq11=0. ; vt_wsum = 0.


 if( trim(spc) == 'qc' .or. trim(spc) == 'qi' )  then  !non-precip

! define imax for a givne radius increment 20. is maximum diameter
  d_increment = 0.001  ! [mm] --> 1micron
  d_min =  0.002  ![mm]  --> 2micron
  d_max = 0.2  ! [mm]  --> 200micron

 else  !precip species

! define imax for a givne radius increment 20. is maximum diameter
  d_increment = 0.10  ! [mm]  --> 100micron
  d_min =  0.2  ![mm]  --> 200micron
  d_max = 20.  ! [mm]  --> 2cm

 endif




!
! special case (Thompson snow class starts from ice crystal size)
!toshii (remove temporally)
  if( trim(cloud_microphysics) == 'THOM' .and. &
      trim(spc) == 'qs'                  ) then  !thompson spectrum density of snow
     d_increment = 0.01  ! [mm] --> 10micron
     d_min =  0.002  ![mm]  --> 2micron
     d_max = 20.  ! [mm]  --> 2cm
!print*, 'special case (Thompson snow class starts from ice crystal size)'
  endif

 imax = nint( (d_max-d_min) /d_increment)


  if( trim(cloud_microphysics) == 'THOM' .and. &
      trim(spc) == 'qs'                  ) then  !thompson spectrum density of snow

      if( .not. allocated(thom_snow_num) ) allocate(  thom_snow_num(0:imax)   )
      call thom_snow_n (wc,temp,d_increment,d_min,d_max,imax, thom_snow_num )

  endif


!
! derive lambda and intercept of expoential DSD for unit water content
!
  call gamma_reff(mu+4.0,gfac1)
  call gamma_reff(mu+3.0,gfac2)
  gamfac = gfac1/gfac2

! derive lambda
  lam = 1. / (2.*re*1.e-6) * gamfac  !slope [1/m]

!
! WDM/WSM specific subroutine
!
  if( (trim(cloud_microphysics) == 'WDM' .or. trim(cloud_microphysics) == 'WSM') .and.  &
      trim(spc) == 'qs'                  ) then  
     lam = 1. / (2.*re*1.e-6)  !slope [1/m]  Eq. (8) in Bae et al. (2016)
  endif

!
! derive intercept
!       [g/m3]                     [m3/g]        [1/m]
  n0  = wc_unit * 6. / const_pi / (dens*1000.) * (lam**(4.+mu))  / gfac1  ! [1/m**(4+mu)]

!
! size loop
!
  do i=0,imax

     d_mm = d_min + d_increment * FLOAT(i)  !diameter [mm]   (new in V3.5)
     d_m  = d_mm * 1e-3                                 ! diameter [m]
     xd  = dble(const_pi*d_mm/wave)                  ! size parameter [-]

!
! Get N(D)
!         [1/m**(4+mu)]  [ m**mu ]      [1/m]  [m]   = [1/m4]
     num = n0           * d_m**mu * exp(-lam * d_m )

!
! Get N(D) for WDM double moment species only. 
!
  if( trim(cloud_microphysics) == 'WDM' ) then
      if( trim(spc) == 'qc' ) call wdm_derive_dsd( 'cloud', wc*1.e-3, re*1.e-6 , d_m, num ) 
      if( trim(spc) == 'qr' ) call wdm_derive_dsd( 'rain' , wc*1.e-3, re*1.e-6 , d_m, num )

  endif

!
! Get N(D) for Thompson snow number density 
!
     if( trim(cloud_microphysics) == 'THOM' .and. &
         trim(spc) == 'qs'                  ) then
       num = thom_snow_num(i) ! [1/m4]
     endif

!
! Calculate effective dielectric constant of frozen particle using different methods/assumptions. 
!
       if( trim(cloud_microphysics) == 'THOM' .and. &
           trim(spc) == 'qs'                  ) then  !thompson spectrum density of snow
          rhoe = min( densice,  (0.069 * d_m*d_m) / ( const_pi / 6. * d_m*d_m*d_m ) ) ! [kg/m3]
       else
          rhoe = dens                  ! effective density [kg/m3]
       endif


!
! Snow compatction impact (Li-Matsui scheme)
!
     fmelt = frac_liq !initialize
     if(trim(spc) == 'qs' .and. frac_liq > 0.0 .and. melt_opt /= 0) then
        call snow_compaction( 0.5*d_mm , rhoe*1e-3 , frac_liq , rad_new, fmelt_new )
        d_mm = 2.*rad_new         ! compacted diameter [mm]
        d_m  = d_mm * 1e-3        ! diameter but using different unit [m]
        xd  = dble(const_pi*d_mm/wave)  ! size parameter [-]
        fmelt = fmelt_new         ! update liquid fraction
     endif

    !
    ! derive effective permittivity (liquid, solid, mixed phases) 
    ! toshiii: found bugs..here fix

     call  effective_permittivity(freq, fmelt, rhoe, temp, crefd)


!if( trim(spc) == 'qs') print*, freq, frac_liq, dens, temp

     !
     ! Mie program (this program will be called anyway regardless of the ss_opt_micro option, 
     !                   since non-spherical routine may fail to derive output. 
     !
     if(xd <= 0. ) stop 'MSG mie_gen_micro: xd become zero or negaive'
     call mie_sphere(xd,crefd,qscad,qextd,asymd,qbscad)

     qext=REAL(qextd)
     qsca=REAL(qscad)
     asym=REAL(asymd)
     qbsca=REAL(qbscad)

!if( trim(spc) == 'qs') then
!  print*,i,qext,asym,xd,crefd
!stop
!endif

     !
     ! Tmatrix routine here (only for non-melting solid-phase particles )
     !
     select case( trim(spc) )
     case( 'qs','qg','qh' )

     if( ss_opt_micro == 4 ) then

       ! solid-sphere radius for interpolation
       rad_solid_micron = 0.5 * 1.e3 * ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![micron]

       !fluffy particle area section for estimating efficiency 
       area = (d_mm*0.5)*(d_mm*0.5)*const_pi  !area of fully particle [mm2]

       call  tmatrix_interp( i_can_say=masterproc, proc='interp', rhoe=dens*1.e-3, &
                             rad_solid=rad_solid_micron, area=area               , &
                             qext=qext_ns, qsca=qsca_ns, asym=asym_ns, qbsca=qbsca_ns, success=success  )

       if( success ) then !if interpolated successfully, use T-matrix value
          !print*,rad_solid_micron , 'Mie qext=',qext, 'Tmatrix qext=',qext_ns
          qext  = qext_ns
          qsca  = qsca_ns
          asym  = asym_ns
          qbsca = qbsca_ns
       endif

     endif
     case default
     end select


     !
     ! SCATDB option for ice crystal (this option easily go size limit, so after calculating mie.)
     !
     if( trim(spc) == 'qi' .and. (ss_opt_micro == 2 .or. ss_opt_micro == 3) )  then

       !solid-sphere radius for interpolation
       d_mm_solid = ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![mm]

       !fluffy particle area section for estimating efficiency 
       area = (d_m*0.5)*(d_m*0.5)*const_pi  !area of fully particle [m2]
       
       call scatdb_interp(nshp, d_mm_solid, area, qext_ns, qsca_ns, asym_ns, qbsca_ns, success)

       if( success ) then !if interpolated, use SCATDB value
          qext  = qext_ns
          qsca  = qsca_ns 
          asym  = asym_ns
          qbsca = qbsca_ns
       endif

     endif


    !
    ! OpenSSP_EXT option for snow aggregate
    !
    if( trim(spc) == 'qs' .and. ss_opt_micro == 3 )  then

       !solid-sphere radius for interpolation (conserving mass)
       d_mm_solid = ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![mm]

       !fluffy particle area section for estimating efficiency 
       area = (d_m*0.5)*(d_m*0.5)*const_pi  !area of fully particle [m2]

       call interp_openssp_ext(freq, area, d_mm_solid, qext_ns, qsca_ns, asym_ns, qbsca_ns, success)

       if( success ) then !if interpolated, use SCATDB value
          !print*,0.5*d_mm_solid*1.e3, 'Mie qext=',qext, 'OpenSSP qext=',qext_ns
          qext  = qext_ns
          qsca  = qsca_ns
          asym  = asym_ns
          qbsca = qbsca_ns
       endif

     endif

     
!
! terminal velocity [m/s]
!
     vt = a_vt * (d_m**b_vt)  ![m/s]

!
! integrate over particle size distribution; (Cross sectin extionction)*(cross sectional aerea)
!
!          [1/m4]*[-]*[-]         * [mm]*[mm]*[mm] * 1.e-6  --> [km-1] ( per 1 g/m3 )

     factor = num * const_pi * 0.25 * d_mm * d_mm * d_increment * 1.e-6

     bext = bext + qext      * factor   ![km-1]  ( per 1 g/m3 )
     bsca = bsca + qsca      * factor   ![km-1]
     bsym = bsym + qsca*asym * factor   ![km-1]
     bq11 = bq11 + qbsca     * factor   ![km-1]
     vt_wsum = vt_wsum + qbsca * factor * vt  !reflectivty weight of terminal velocity

  enddo !i loop

!
! check for distribution with very small extinction;
! set parameters to zero to avoid numerical problems
!

  if( bext .gt. min_bext) then
      ksca=bext      * wc  ! *wc is re-scaling to actual extinction  [km-1]
      asca=bsca/bext       ! [-]
      gsca=bsym/bsca       ! [-]
      pbck=bq11/bsca       ! [-]   
      vdop = vt_wsum / bq11 !doppler velocity [m/s]

   else
      ksca=0. ; asca=0. ; gsca=0. ; pbck=0. ; vdop=0.
   end if

! print*,'MSG mie_gen_micro', ksca, asca, gsca

!if( trim(spc) == 'qs') stop


 return
 end subroutine mie_gen_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine mie_rams_micro(spc,freq, temp, wc, re, gnu, cfvt, pwvt, cfmas, pwmas, frac_liq, &
                            ksca, asca, gsca, pbck, vdop)
  implicit none

!--------------------------------------------------------------------------------
! Purpose : Compute the extinction, absorption, asymmetry parameter and
!           backscatter for a given water content and effective radius of condensates , and
!           a particle size distribution parameter for RAMS 1moment & 2moment microphyiscs
!
!  08/2016  Toshi Matsui @ NASA GSFC : Clean up by adding effective_permittivity routine. 
!  07/2015  Toshi Matsui@NASA GSFC ; Size ranges changed between precip and non-precip species.
!  01/2014  Toshi Matsui@NASA GSFC ; add Tmatrix option
!  01/2011  Toshi Matsui@NASA GSFC ; add SCATDB option
!  03/2010  Toshi Matsui@NASA GSFC ; Additional options for effective refractive index for solid particles.
! Toshi Matsui @ NASA GSFC: Initial
!
!--------------------------------------------------------------------------------
 character*(*),intent(in) :: spc !species character
 real(sdsu_fps),intent(in) :: freq  ! frequency of radiation [GHz]
 real(sdsu_fps),intent(in) :: temp  ! temperature of particles [K]
 real(sdsu_fps),intent(in) :: wc    ! water content of condensate [g/m3]
 real(sdsu_fps),intent(in) :: re    ! drop effective radius [micron] 
 real(sdsu_fps),intent(in) :: gnu   ! PSD shape parameter for generalized gamma distribution

 real(sdsu_fps),intent(in) :: cfvt   ! alpha_vt 
 real(sdsu_fps),intent(in) :: pwvt   ! beta_vt  [-]
 real(sdsu_fps),intent(in) :: cfmas  ! alpha_m [kg /m**beta_m]
 real(sdsu_fps),intent(in) :: pwmas  ! beta_m  [-]
 real(sdsu_fps),intent(in) :: frac_liq ! fracion of liquid (1=rain or 0=ice condensates, 0~1=melting snow)
 real(sdsu_fps),intent(out) :: ksca  ! extinction coefficient [1/km]
 real(sdsu_fps),intent(out) :: asca  ! single-scatter albedo [-]
 real(sdsu_fps),intent(out) :: gsca  ! asymmetry factor [-]
 real(sdsu_fps),intent(out) :: pbck  ! backscatter phase function/(4*pi) [-]
 real(sdsu_fps),intent(out) :: vdop  ! Doppler velocity [m/s]

!--------Local Parameters
 integer :: i    !looping 
 integer :: imax !max for looping
 real(sdsu_fps) :: rhoe ! bulk effective density of hydrometero [kg/m**3]

 real(sdsu_fps) ::& 
   wave  ,& ! wavelength [mm]
   rad   ,& ! radius of particle [mm]
   d     ,& ! particle diameter 
   dens  ,& ! particle density [kg/m3]
   num   ,& ! a particle number density per radius increment [1/m4]
   faa   ,& ! fraction of air in ice 
   eicere,& ! real component of ice complex refractive index
   eiceim,& ! imaginary component of ice complex refractive index
   qext  ,& ! extinction
   qsca  ,& ! scattering
   asym  ,& ! asymetry 
   qbsca ,& ! backscattering
   bext  ,& ! extinction
   bsca  ,& ! scattering
   bsym  ,& ! asymetry
   bq11  ,& ! backscattering
   eimag ,& ! complex refractive index
   ereal ,& !
   gamfac,& ! gamma PSD factors [-]
   gfac1 ,& !
   gfac2 ,& !
   dn         ,& ! charactristic diameter [m]
   d_increment,& ! diameter increment over size loops [mm]
   d_min      ,& ! minium diameter of particle [mm]
   d_max      ,& ! maximum diameter of particle [mm]
   d_m        ,& ! particle diameter [mm]
   d_mm       ,& ! particle diameter [mm]
   rad_new    ,& ! updated particle radius [mm]
   fmelt      ,& ! melt fraction
   fmelt_new  ,& ! updated melt fraction [-]
   mean_mass  ,& ! mean mass [kg]
   fgamma     ,& ! gamma function
   ntot       ,& ! total particle number concentrations [#/m3]
   vt_wsum    ,& ! weighting sum of terminal velocity 
   vt         ,& ! terminal velocity [m/s]
   factor        ! factor = num * pi * 0.25 * d_mm * d_mm * d_increment * 1.e-6

 real(sdsu_fpd) :: &
     xd      ,&      ! size parameter [-] in double precision 
     freqd   ,&
     freqhzd ,& ! frequency in double precision [GHz]
     tempd   ,&      ! temperaure in doubpl precision
     sald    ,&      ! salinity of water
     eicered ,&
     ewatred ,& !real component of ice/water complex refractive index
     eiceimd ,&
     ewatimd ,& !imaginary component of ice/water complex refractive index
     qscad   ,&
     qextd   ,&
     asymd   ,&
     qbscad
    
 complex *8 eice,ewat,eair !ice, water, and air permittivity
 complex *8 ei     !effective ice-air permittivity
 complex *8 emelt,emelt_wi,emelt_iw, emelt_em !effective water-ice-air permittivity

 complex*16 ewatd
 complex*16 eid, emeltd
 complex*16 crefd, crefd_liq, crefd_ice

 real(sdsu_fps),parameter :: dr = 0.05          ! increment bin of radius [mm]
 real(sdsu_fps),parameter :: densice = 0.917e+3 ! solid ice density [kg/m3]
 real(sdsu_fps),parameter :: densliq = 1.0e+3   ! liquid water density [kg/m3]
 real(sdsu_fps),parameter :: wc_unit = 1.0      ! unit water content [g/m3]
 integer,save :: nshp !SCATDB shape #
 real(sdsu_fps) :: qext_ns, qsca_ns, asym_ns, qbsca_ns ! non-spherical single scattering
 real(sdsu_fps) :: d_mm_solid  !solid ice diameter [mm]
 logical :: success ! status 
 real(sdsu_fps) :: rad_solid_micron  !solid ice radius [micron]
 real(sdsu_fps) :: area   !area of fully particle [mm2]
! real(sdsu_fps) :: dcreal, dcimag  !real / imageneary component of dielectric function

!
! Begin by checking if hydrometeors of this species have enough water content.
! If not, set scattering parameters to zero and return.
!
  if(wc .lt. q_min_micro .or. re == 0.) then
     ksca=0. ; asca=0. ; gsca=0. ; pbck=0. ; vdop=0.
     return
  endif

!
! SCATDB option for ice crystal
!
 if( trim(spc(1:3)) == 'ice' .and. ss_opt_micro == 2 )  then
    if(trim(spc) == 'ice1_col') nshp = 2
    if(trim(spc) == 'ice1_hex') nshp = 4
    if(trim(spc) == 'ice1_den') nshp =10
    if(trim(spc) == 'ice1_ndl') nshp = 0
    if(trim(spc) == 'ice1_ros') nshp = 8 
    if(trim(spc) == 'ice2_col') nshp = 2
    if(trim(spc) == 'ice2_hex') nshp = 4
    if(trim(spc) == 'ice2_den') nshp =10
    if(trim(spc) == 'ice2_ndl') nshp = 0
    if(trim(spc) == 'ice2_ros') nshp = 8 
    call scatdb_prep(freq, temp, nshp)
 endif


!
! preparting ice Tmatrix option  (so far does not apply to melting species or liquid species)
!
 if( frac_liq == 0.e0 .and. ss_opt_micro == 4 ) then
   call  tmatrix_interp( i_can_say=masterproc, proc='prep', freq=freq )
 endif



!
! Assign some useful constants
!

  wave=const_ghz_mm/freq !wavelgnth [mm] <- frequency [GHz]


!
! Get characteristic diameter from effective radius (re) and PSD shape parameter (gnu)
!
  call gamma_reff(gnu+3.0,gfac1)
  call gamma_reff(gnu+2.0,gfac2)    
  gamfac = gfac1/gfac2
  dn = 2.0 * re / gamfac * 1e-6  !charactristic diameter [m]


!
! get total number concentration for a given characteristic diameter, mass mixing ratio, PSD shape parameter
! , and gamma factors.
!
  call gamma_reff(gnu+pwmas,gfac1)
  call gamma_reff(gnu      ,gfac2)
  gamfac = gfac1/gfac2
  mean_mass = cfmas * ( dn ** pwmas ) * gamfac  !mean mass [kg]
  ntot = wc * 1e-3 / mean_mass                   !total particle number concentration [#/m3]


!
! initialize the optical parameters
!
  bext=0. ; bsca=0. ; bsym=0. ; bq11=0.  ; vt_wsum = 0.0

!
!toshii: Should I extend PSD beyond 200micron for ice?
! This is related to the too small Zh values in anvil top CFADs.
!
 if( trim(spc) == 'cloud1' .or. trim(spc(1:3)) == 'ice' )  then  !non-precip

! define imax for a givne radius increment 20. is maximum diameter
  d_increment = 0.001  ! [mm] --> 1micron
  d_min =  0.002  ![mm]  --> 2micron
  d_max = 0.2  ! [mm]  --> 200micron

 else  !precip species

! define imax for a givne radius increment 20. is maximum diameter
  d_increment = 0.10  ! [mm]  --> 100micron
  d_min =  0.2  ![mm]  --> 200micron
  d_max = 20.  ! [mm]  --> 2cm

 endif

 imax = nint( (d_max-d_min) /d_increment)

!
! loop over particle size
!
  SIZE_LOOP: do i=0,imax

!
! Compute diameter, size parameter particle number density, and particle density
!
     d_mm = d_increment * 0.5 + d_increment * FLOAT(i)  !diameter [mm]
     d_m  = d_mm * 1e-3                        ! diameter but using different unit [m]
     xd  = dble(const_pi*d_mm/wave)                  ! size parameter [-]
     dens = (6.0 * cfmas / const_pi) * (d_m ** (pwmas-3.e0)) ! particle density [kg/m3]
     rhoe = dens ! particle density [kg/m3]

     fmelt = frac_liq  !melt fraction

!
! Snow compatction impact (Li-Matsui scheme)
!
     
     if(trim(spc) == 'snow' .and. frac_liq > 0.0 .and. melt_opt /= 0) then
        call snow_compaction( 0.5*d_mm , dens*1e-3 , frac_liq , rad_new, fmelt_new )
        d_mm = 2.*rad_new         ! compacted diameter [mm]
        d_m  = d_mm * 1e-3        ! diameter but using different unit [m]
        xd  = dble(const_pi*d_mm/wave)  ! size parameter [-]
        fmelt = fmelt_new         ! update liquid fraction
     endif 

!
! Gamma function
!
     call gamma_function(gnu,d_m,dn, fgamma)   ! get gamma function
     num = ntot * fgamma                       ! particle number density[1/m4]

     tempd=dble(temp)  !temperature [K] in double precision


    !
    ! derive effective permittivity (liquid, solid, mixed phases)
    !
    call  effective_permittivity(freq, fmelt, rhoe, temp, crefd )

 
     !
     ! Mie program (this program will be called anyway regardless of the ss_opt_micro option, 
     !                   since non-spherical routine may fail to derive output. 
     !
     call mie_sphere(xd,crefd,qscad,qextd,asymd,qbscad)

     qext = REAL(qextd)   !extinction efficiency 
     qsca = REAL(qscad)   !scattering 
     asym = REAl(asymd)   !asymetry
     qbsca= REAL(qbscad)  !backscattering

     !
     ! Tmatrix routine here (only for non-melting solid-phase particles )
     !
     if( frac_liq == 0.e0 .and. ss_opt_micro == 4 ) then

       ! solid-sphere radius for interpolation
       rad_solid_micron = 0.5 * 1.e3 * ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![micron]

       !fluffy particle area section for estimating efficiency 
       area = (d_mm*0.5)*(d_mm*0.5)*const_pi  !area of fully particle [mm2]

       call  tmatrix_interp( i_can_say=masterproc, proc='interp', rhoe=dens*1.e-3, &
                             rad_solid=rad_solid_micron, area=area               , &
                             qext=qext_ns, qsca=qsca_ns, asym=asym_ns, qbsca=qbsca_ns, success=success  )

       if( success ) then !if interpolated successfully, use T-matrix value

          qext  = qext_ns
          qsca  = qsca_ns
          asym  = asym_ns
          qbsca = qbsca_ns
       endif

     endif

    !
    ! SCATDB option for ice crystal
    !
    if( trim(spc(1:3)) == 'ice' .and. (ss_opt_micro == 2 .or. ss_opt_micro == 3 ) )  then

       !solid-sphere radius for interpolation
       d_mm_solid = ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![mm]

       !fluffy particle area section for estimating efficiency 
       area = (d_m*0.5)*(d_m*0.5)*const_pi  !area of fully particle [m2]

       call scatdb_interp(nshp, d_mm_solid, area, qext_ns, qsca_ns, asym_ns, qbsca_ns, success)

       if( success ) then !if interpolated, use SCATDB value

          qext  = qext_ns
          qsca  = qsca_ns
          asym  = asym_ns
          qbsca = qbsca_ns
       endif

     endif

    !
    ! OpenSSP_EXT option for snow aggregate
    !
    if( trim(spc) == 'snow' .and. ss_opt_micro == 3 )  then

       !solid-sphere radius for interpolation
       d_mm_solid = ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![mm]

       !fluffy particle area section for estimating efficiency 
       area = (d_m*0.5)*(d_m*0.5)*const_pi  !area of fully particle [m2]

       call interp_openssp_ext(freq, area, d_mm_solid, qext_ns, qsca_ns, asym_ns, qbsca_ns, success)

       if( success ) then !if interpolated, use SCATDB value
          !print*,0.5*d_mm_solid*1.e3, 'Mie asym=',asym, 'OpenSSP asym=',asym_ns
          qext  = qext_ns
          qsca  = qsca_ns
          asym  = asym_ns
          qbsca = qbsca_ns
       endif

     endif

!
! terminal velocity [m/s]
!
     vt = cfvt * (d_m**pwvt)  ![m/s]


!
! integrate over particle size distribution; (Cross sectin extionction)*(cross sectional aerea)
!
!          [1/m4]*[-]*[-]         * [mm]*[mm]*[mm] * 1.e-6  --> [km-1] ( per 1 g/m3 )

     factor = num * const_pi * 0.25 * d_mm * d_mm * d_increment * 1.e-6  

     bext = bext + qext      * factor   ![km-1]  ( per 1 g/m3 )
     bsca = bsca + qsca      * factor   ![km-1]
     bsym = bsym + qsca*asym * factor   ![km-1]
     bq11 = bq11 + qbsca     * factor   ![km-1]
     vt_wsum = vt_wsum + qbsca * factor * vt  !reflectivty weight of terminal velocity

  enddo SIZE_LOOP

!
! check for distribution with very small extinction;
! set parameters to zero to avoid numerical problems
!

  if( bext .gt. min_bext) then
      ksca=bext         ! [km-1]
      asca=bsca/bext    ! [-]
      gsca=bsym/bsca    ! [-]
      pbck=bq11/bsca    ! [-]
      vdop = vt_wsum / bq11 !doppler velocity [m/s]
   else
      ksca=0. ; asca=0. ; gsca=0. ; pbck=0. ; vdop=0.
   end if

 return
 end subroutine mie_rams_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine mie_sbm_micro(spc, freq, temp, bulkden, bulkrad, dbulkrad, vt, wc, n0, fmelt_bin,frime_bin,&
                           ksca, asca, gsca, pbck, vdop)
  implicit none
!----------------------------------------------------------------------------------------------------
!  Purpose 
!     Compute the extinction, absorption, asymmetry parameter and
!     backscatter for a given # conc. of ice species (cloud ice, snow, graupel, hail) .
!
!
! 08/2016  Toshi Matsui @ NASA GSFC : Clean up by adding effective_permittivity routine. 
! 01/2011  Toshi Matsui@NASA GSFC ; add SCATDB option
! 06/2011  Toshi Matsui @ NASA GSFC : adding rime fraction, adapting for HUCM_SBM43 scheme. 
! 03/2010  Toshi Matsui @ NASA GSFC ; Additional options for effective refractive index for solid particles.
! 08/2007  Toshi Matsui @ NASA GSFC ; initial
!----------------------------------------------------------------------------------------------------
 character(len=*),intent(in):: spc      ! index name of hydrometeor classes 
 real(sdsu_fps), intent(in) :: freq     ! frequency of radiation [GHz]
 real(sdsu_fps), intent(in) :: temp     ! temperature of particles [K]
 real(sdsu_fps), intent(in) :: bulkden(1:nbin) ! bulk density of hydrometero for each bins [g/cm3]
 real(sdsu_fps), intent(in) :: bulkrad(1:nbin) ! bulk radius of hydrometero for each bins [cm] 
 real(sdsu_fps), intent(in) :: dbulkrad(1:nbin)! binsize [cm]
 real(sdsu_fps), intent(in) :: vt(1:nbin)      ! terminal velocity [cm/s]
 real(sdsu_fps), intent(in) :: wc        ! water content of cloud ice distribution [g/m**3]
 real(sdsu_fps), intent(in) :: n0(1:nbin)! number concentration of hydrometeors [1/m4]
 real(sdsu_fps), intent(in) :: fmelt_bin(1:nbin) !mass melting fraction of hydrometoer for each bins [-] 
 real(sdsu_fps), intent(in) :: frime_bin(1:nbin) !mass rime fraction of snow only for each bin [-] 
 real(sdsu_fps), intent(out) :: ksca     ! extinction coefficient [1/km]
 real(sdsu_fps), intent(out) :: asca     ! single-scatter albedo [-]
 real(sdsu_fps), intent(out) :: gsca     ! asymmetry factor [-]
 real(sdsu_fps), intent(out) :: pbck     ! backscatter phase function/(4*pi) [-]
 real(sdsu_fps), intent(out) :: vdop     ! Doppler velocity [m/s]

!----Local variables
 real(sdsu_fps), parameter :: densice = 0.917e+3 ![kg/m3]
 real(sdsu_fps), parameter :: densliq = 1.0e+3 ![kg/m3]
 real(sdsu_fps) :: rhoe ! bulk effective density of hydrometero [kg/m**3]

 integer *4 i
 real(sdsu_fps) :: eimag, ereal

 real(sdsu_fps) :: &
      wave ,&    ! wavelength [mm]
      rad  ,&    ! radius of particle [mm]
      d_mm ,&    ! diameter [mm]
      d_m  ,&    ! diameter [m]
      rad_liq ,& ! new adjusted radius for liquid drop from snow melted portion [mm]
      dens ,&    ! effective density [kg/m3]
      fmelt ,&   ! melting fraction
      num    ,&  ! number conc. per radius bin [1/m4]
      faa    ,&  ! fraction of air in solid particles [-]
      den_diag,& ! diagnostic density [g/cm3]
      rad_diag,& ! diagnostic radisu [cm]
      eicere ,&  ! real component of ice complex refractive index
      eiceim ,&  ! imaginary component of ice complex refractive index
      ewatim ,&  
      qext   ,&  ! extinction
      qsca   ,&  ! backscattering
      asym   ,&  ! asymetry
      qbsca  ,&
      bext   ,&
      bsca   ,&
      bsym   ,&
      bq11   ,&
      vt_snow,&  !snow terminal velocity [cm/s]
      vt_wsum    !weighted sum of terminal velocity of a particle


 real(sdsu_fpd) :: &
      sald     ,&  ! salinity of water
      xd       ,&  ! size parameter [-] in double precision
      freqd    ,&
      freqhzd  ,&  ! frequency in double precision [GHz]
      tempd    ,&
      eicered  ,&
      ewatred  ,&  ! real component of ice/water complex refractive index
      eiceimd  ,&
      ewatimd  ,&  ! imaginary component of ice/water complex refractive index
      qscad    ,&
      qextd    ,&
      asymd    ,&
      qbscad
 real(sdsu_fps) :: factor
 real(sdsu_fps) :: drad

 complex *8 eice,ewat,eair !ice, water, and air permittivity
 complex *8 ei     !effective ice-air permittivity
 complex *8 emelt,emelt_wi,emelt_iw, emelt_em !effective water-ice-air permittivity

 complex *16 eid, emeltd, ewatd
 complex *16 crefd, crefd_ice, crefd_liq
 integer,save :: nshp !SCATDB shape #
 real(sdsu_fps) :: qext_ns, qsca_ns, asym_ns, qbsca_ns ! non-spherical single scattering
 real(sdsu_fps) :: d_mm_solid  !solid ice diameter [mm]
 real(sdsu_fps) :: rad_solid_micron  !solid ice radius [micron]
 real(sdsu_fps) :: area   !area of fully particle [mm2]
 logical :: success ! status 

 character(len=5),parameter :: method = 'new'
 real(sdsu_fps) :: fv_melt
 integer :: melt_temp 
! real(sdsu_fps) :: dcreal, dcimag  !real / imageneary component of dielectric function

!
!     Begin by checking if hydrometeors of this species are present.
!     If not, set scattering parameters to zero and return.
!
 ksca=0. ; asca=0. ; gsca=0. ; pbck=0. ; vdop = 0.

 if(wc .lt. q_min_micro ) then 
     return
 endif

!
! SCATDB option for ice crystal
!
 if( trim(spc(1:3)) == 'ice' .and. ss_opt_micro == 2 )  then
    if(trim(spc) == 'ice1_col') nshp = 1
    if(trim(spc) == 'ice1_pla') nshp = 3
    if(trim(spc) == 'ice1_den') nshp =10
    call scatdb_prep(freq, temp, nshp)
 endif


!
! preparting ice Tmatrix option  (so far does not apply to melting species or liquid species)
!
 if( trim(spc) /= 'liq' .and. ss_opt_micro == 4 ) then
   call  tmatrix_interp( i_can_say=masterproc, proc='prep', freq=freq )
 endif


!
!     Assign some useful constants
 wave=const_ghz_mm/freq

!
!     If hydrometeors are present, initialize the scattering parameters
!
 bext=0. ; bsca=0. ; bsym=0. ; bq11=0. ; vt_wsum = 0.

!
!     Loop over particle sizes:
 do i=1,nbin

    if(n0(i) == 0. ) cycle  !if number density is zero, skip all process


    if( trim(spc) == 'snow')   then        


      call diag_rime_melt_snow (i, frime_bin(i), fmelt_bin(i), den_sbm(1:nbin)%snow, &
                                x_sbm(1:nbin)%snow, rad_sbm(1:nbin)%snow, &
                                den_diag, rad_diag, fv_melt )

      rad  = rad_diag*10.e0    ! radius [mm]
      d_mm  = rad * 2.         ! diameter [mm]
      d_m   = d_mm * 1.e-3     ! diameter [m]
      drad = dbulkrad(i)*10.e0 ! binsize [mm]
      dens = den_diag*1.e3     ! effective density [kg/m3]
      fmelt = fv_melt          ! volumetric melting fraction [-]


    elseif( trim(spc) == 'graupel' .or. trim(spc) == 'hail')   then

      rad   = bulkrad(i)*10.e0  ! radius [mm]
      d_mm  = rad * 2.          ! diameter [mm]
      d_m   = d_mm * 1.e-3     ! diameter [m]
      drad  = dbulkrad(i)*10.e0 ! binsize [mm]
      dens  = bulkden(i) * 1.e3 ! effective density [kg/m3]
      call diag_fv_melt( fmelt_bin(i) , bulkden(i) , fv_melt ) !diagnostic volumetric liquid fraction
      fmelt = fv_melt      ! melting fraction [-]

    else

      rad   = bulkrad(i)*10.e0  ! radius [mm]
      d_mm  = rad * 2.          ! diameter [mm]
      d_m   = d_mm * 1.e-3     ! diameter [m]
      drad  = dbulkrad(i)*10.e0 ! binsize [mm]
      dens  = bulkden(i) * 1.e3 ! effective density [kg/m3]
      fmelt = 0.                ! melting fraction [-]

    endif

    rhoe = dens ! effective density [kg/m3]

    xd=dble(2.e0*const_pi*rad/wave)  ! size parameter

    num = n0(i)    ! number conc. per radius bin [1/m4]

   !
   ! derive effective permittivity (liquid, solid, mixed phases)
   !
    call  effective_permittivity(freq, fmelt, rhoe, temp, crefd)


     !
     ! Mie program (this program will be called anyway regardless of the ss_opt_micro option, 
     !                   since non-spherical routine may fail to derive output. 
     !
     call mie_sphere(xd,crefd,qscad,qextd,asymd,qbscad)

     qext=REAL(qextd)   ! extinction effciency [-]
     qsca=REAL(qscad)   ! 
     asym=REAL(asymd)
     qbsca=REAL(qbscad)

     !
     ! Tmatrix routine here (only for non-melting solid-phase particles )
     !
     if( trim(spc) /= 'liq' .and. fmelt_bin(i) == 0.e0 .and. ss_opt_micro == 4 ) then

       ! solid-sphere radius for interpolation
       rad_solid_micron = 0.5 * 1.e3 * ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![micron]

       !fluffy particle area section for estimating efficiency 
       area = (d_mm*0.5)*(d_mm*0.5)*const_pi  !area of fully particle [mm2]

       call  tmatrix_interp( i_can_say=masterproc, proc='interp', rhoe=dens*1.e-3, &
                             rad_solid=rad_solid_micron, area=area               , &
                             qext=qext_ns, qsca=qsca_ns, asym=asym_ns, qbsca=qbsca_ns, success=success  )

       if( success ) then !if interpolated successfully, use T-matrix value
          qext  = qext_ns
          qsca  = qsca_ns
          asym  = asym_ns
          qbsca = qbsca_ns
       endif

     endif



    !
    ! SCATDB option for ice crystal
    !
    if( trim(spc(1:3)) == 'ice' .and. (ss_opt_micro == 2 .and. ss_opt_micro == 3) )  then

       !solid-sphere radius for interpolation
       d_mm_solid = ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![mm]

       !fluffy particle area section for estimating efficiency 
       area = (d_m*0.5)*(d_m*0.5)*const_pi  !area of fully particle [m2]

       call scatdb_interp(nshp, d_mm_solid, area, qext_ns, qsca_ns, asym_ns, qbsca_ns, success)

       if( success ) then !if interpolated, use SCATDB value
          qext  = qext_ns
          qsca  = qsca_ns
          asym  = asym_ns
          qbsca = qbsca_ns
       endif

    endif

    !
    ! OpenSSP_EXT option for snow aggregate
    !
    if( trim(spc) == 'snow' .and. ss_opt_micro == 3 )  then

       !solid-sphere radius for interpolation
       d_mm_solid = ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![mm]

       !fluffy particle area section for estimating efficiency 
       area = (d_m*0.5)*(d_m*0.5)*const_pi  !area of fully particle [m2]

       call interp_openssp_ext(freq, area, d_mm_solid, qext_ns, qsca_ns, asym_ns, qbsca_ns, success)

       if( success ) then !if interpolated, use SCATDB value
          !print*,0.5*d_mm_solid*1.e3, 'Mie asym=',asym, 'OpenSSP asym=',asym_ns
          qext  = qext_ns
          qsca  = qsca_ns
          asym  = asym_ns
          qbsca = qbsca_ns
       endif

     endif



!
! integrate over particle size distribution;
! num*drad becomes actual number concentraions [#/m3] (new) : Mar 30, 2010
!          [1/m4]       [mm][mm] [mm]         --> [km-1] 
   factor = num*const_pi*rad*rad*drad*1.e-6

   bext=bext+ qext     * factor   ![km-1]
   bsca=bsca+ qsca     * factor   ![km-1]
   bsym=bsym+ qsca*asym* factor   ![km-1]
   bq11=bq11+ qbsca    * factor   ![km-1]

   if( spc == 'snow   ')   then
      !
      ! diagnostic snow aggregate Vt as a function of riming and meling fracitons.
      !
      call vt_rime_melting( i, nbin, frime_bin(i), fmelt_bin(i), x_sbm%snow , vt_sbm%liq , &
                            vt_sbm%snow, vt_sbm%hail, den_sbm%snow, vt_snow ) 
!multiply 0.01 to convert cm/s to m/s
      vt_wsum = vt_wsum + qbsca * factor * 0.01*vt_snow !reflectivty weight of terminal velocity [m/s]
   else !default  
      vt_wsum = vt_wsum + qbsca * factor * 0.01*vt(i)   !reflectivty weight of terminal velocity [m/s]
   endif

 enddo !bin spectra loop

!
!     check for distribution with very small extinction;
!     set parameters to zero to avoid numerical problems
!
  if( bext .gt. min_bext) then
     ksca=bext         ![km-1]
     asca=bsca/bext    ![-]
     gsca=bsym/bsca    ![-]
     pbck=bq11/bsca    ![-]
     vdop = vt_wsum / bq11 * 1.e-2 !doppler velocity [m/s]
  else
     ksca=0. ; asca=0. ; gsca=0. ; pbck=0. ; vdop = 0.
  end if

  return
 end subroutine mie_sbm_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine snow_compaction( rad , den , fmelt , rad_new , fmelt_new )
 implicit none
!--------------------------------------------------------------------------------
! Purpose : Compact snow aggregate radius as a function of melting fraction. 
! 
! Toshi Matsui @ NASA GSFC: initial
!--------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: rad        ! initial particle radius [mm]
 real(sdsu_fps),intent(in) :: den        ! particle density [g/cm3]
 real(sdsu_fps),intent(in) :: fmelt      ! 0~1 [-]
 real(sdsu_fps),intent(out) :: rad_new   ! compacted particle radius  [mm]
 real(sdsu_fps),intent(out) :: fmelt_new ! 0~1 [-]
 real(sdsu_fps) :: rad_liq ![mm]
 real(sdsu_fps),parameter :: den_liq = 1.0 !liquid density [g/cm3]


  rad_liq = ( (rad**3.)*den*fmelt / den_liq  )**(1./3.)  !radius for liquid water [mm]

  rad_new = ( ( (rad**3.)*den + (rad_liq**3.)*(den-den_liq) )/den  )**(1./3.)  ! [mm]

  fmelt_new = (rad_liq/rad_new)**3  !this is actually volumetric fraction of water within melting ice-air inclusion

  return
 end subroutine snow_compaction

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mg_ellips (fincl, ematrix, eincl, emg)
 implicit none
!--------------------------------------------------------------------------------
! Purpose :!
! Maxwell-Garnett dielctric function of 2-component media
! (elliptical inclusions) P. Bauer 1996 (core-shell mixture)
! See also Olson et al. 2001 JAM page 1153 equation (32)&(33)
! This version is closer to EM sollution, and slighly different 
! from Bohren and Battan 1980 (spherical inclusion). 
! 
! 08/2007  Toshi Matsui @ NASA GSFC ; initial
!-------------------------------------------------------------------------------

 real(sdsu_fps),intent(in) ::  fincl  !volume fraction of inclusions
 complex*8,intent(in) :: ematrix ! permittivity of matrix
 complex*8,intent(in) :: eincl   ! permittivity of inclusions
 complex*8,intent(out) ::  emg   !  effective permittivity
 complex*8 :: gamma, q


      q     = (eincl / (eincl - ematrix)) &
            * clog (eincl / ematrix) - 1.0
      gamma = 2.0 * ematrix * q / (eincl - ematrix)

      emg = ((1.0 - fincl) * ematrix + fincl * gamma * eincl) &
          / (1.0 - fincl + fincl * gamma)

    return
  end subroutine mg_ellips

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine em_ellips (fincl, ematrix, eincl, emg)
 implicit none

!-------------------------------------------------------------------------------
! Effective-Medium Debye dielectric function  (Homogeneous mixture assumption)
! See eq. (2) in Bohren and Battan 1980
!
! Toshi Matsui@NASA GSFC; Initial
!-------------------------------------------------------------------------------

 real(sdsu_fps),intent(in) ::  fincl  !volume fraction of inclusions
 complex*8,intent(in) :: ematrix ! permittivity of matrix
 complex*8,intent(in) :: eincl   ! permittivity of inclusions
 complex*8,intent(out) :: emg   !  effective permittivity
 complex*8 :: a,b,c

!
! parameters for a quadratic equation
!
 a = cmplx(-2.,0.)
 b = (3.*fincl-1.)*eincl - (3.*fincl-2.)*ematrix 
 c = ematrix*eincl

!
! quadratic solution for effective permittivity 
!
 emg = (-b - sqrt(b*b-4.*a*c) )  / (2.*a)

    return
  end subroutine em_ellips

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine debye_mixture(rho,eprsol,episol,ei)
 implicit none
!
! Debye Theory is used to determine the relative permittivity/dielectric
! of an air/ice mixture given the density of this mixture.  See Battan
! (1973) p. 40.
! LDC 3/21/98
!
 real(sdsu_fpd),intent(in) :: rho, eprsol, episol
 complex,intent(out) :: ei

 real(sdsu_fpd) :: epr, epi
 complex :: m2, m2sol, msol, kovrrho, rhosol, rhomix
 real(sdsu_fpd) :: nsol,ksol,n,k

!
! creating complex function
!
 rhomix = cmplx(rho, 0.)

!
! Densitiy of solid ice
!
  rhosol=(0.917,0.)

!
! Complex index of refraction for water and solid ice (rho=0.917 g/cm^3)
!  is obtained from subroutine ray which gives relative permittivity.
! Must convert first.
! approximate since n**2>>k**2 (n**2-k**2 ~= n**2) for ice
!
  nsol=sqrt(eprsol)
  ksol=episol/(2*nsol)
  msol=cmplx(nsol,ksol)
  m2sol=msol**2

!
! Determin ice refractive index and dielectric using Debye theory (1929)
! which postulataes that K/rho is constant for an ice/air mixture
! where rho is density and  K is (m**2-1)/(m**2+2) Battan (1973) p. 40.
!

  kovrrho=((m2sol-(1,0))/(m2sol+(2,0)))/(rhosol)
  m2=((2,0)*kovrrho*rhomix+(1,0))/((1,0)-kovrrho*rhomix)
  n=real(sqrt(m2))
  k=imag(sqrt(m2))
  epr=n**2-k**2   !real component
  epi=2*n*k       !imagenary component
  ei=cmplx(epr,epi)  !output 

  return
 end subroutine debye_mixture

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine mie_sphere (x, min, qscat, qexti, asym, qbscat)

!*
!*    Mie Routine P. Bauer
!*

      implicit none
!*
      save
      integer limitx
      parameter (limitx = 1500)
      real(sdsu_fpd) conv, pi, ghztohz
!*
      real(sdsu_fpd) :: x
      real(sdsu_fpd) :: mr, mi, n1, n2
      real(sdsu_fpd) :: qscat, qexti, qabso, asym, qbscat
!*
      real(sdsu_fpd) :: rfac1, rfac2
      real(sdsu_fpd) :: rhelp1 (2), rhelp2 (2)
!*
      complex*16 m, mx, min
      complex*16 chelp1, chelp2, cfac1, cfac2, cbscat
!*
      complex*16 dn (0:limitx), wn (-1:limitx)
      complex*16 an   (limitx), bn    (limitx)
!*
      integer nend
      integer i100, i101
!*
      equivalence (chelp1, rhelp1 (1))
      equivalence (chelp2, rhelp2 (1))
!*
      real(sdsu_fpd) ::  eps0, epshigh
!*
!*
      data pi      /  3.14159265358979323846d+00 /
      data conv    /  0.017453292d+00 /
      data ghztohz /  1.0d+09 /
!*
      data eps0     / 8.854d-12 /
      data epshigh  / 4.9d+00 /

!*
      m      = conjg (min)
      chelp1 = m
      mr     =        rhelp1 (1)
      mi     = -1.0d0 * rhelp1 (2)
!*
      mx   = m  * x
      n1   = mr * x
      n2   = mi * x
!*
      if (x .le. 20000.0) nend = int( x + 4.00 * x ** (1.0 / 3.0) + 2.0 )
      if (x .le.  4200.0) nend = int( x + 4.05 * x ** (1.0 / 3.0) + 2.0 )
      if (x .le.     8.0) nend = int( x + 4.00 * x ** (1.0 / 3.0) + 1.0 )
      if (nend .le.    5) nend = 5
      if (nend .gt. limitx) nend = limitx
!*
      rfac1      = sin  (n1) * sin  (n1) + sinh (n2) * sinh (n2)
      rhelp1 (1) = sin  (n1) * cos  (n1) / rfac1
      rhelp1 (2) = sinh (n2) * cosh (n2) / rfac1
!*
      dn (0) = chelp1
!*
      rhelp1 (1) =             cos (x)
      rhelp1 (2) = -1.0d+00 * sin (x)
      rhelp2 (1) =             sin (x)
      rhelp2 (2) =             cos (x)
!*
      wn (-1) = chelp1
      wn ( 0) = chelp2
!*
      qexti  = 0.0
      qscat  = 0.0
      qbscat = 0.0
      qabso  = 0.0
      asym   = 0.0
      cbscat = cmplx (0.0,0.0)
!*
      do 100 i100 = 1, nend
         dn (i100) = -1.0d0 * i100 / mx &
                   +  1.0d0 / (i100 / mx - dn (i100 - 1))
         wn (i100) = wn (i100 - 1) * (2.0d0 * i100 - 1.0d0) / x &
                   - wn (i100 - 2)
!*
         cfac1 = dn (i100) / m + i100 / x
         cfac2 = m * dn (i100) + i100 / x
!*
         chelp1 = wn (i100)
         chelp2 = wn (i100 - 1)
!*
         an (i100) = (cfac1 * rhelp1 (1) - rhelp2 (1)) &
                   / (cfac1 * chelp1     - chelp2    )
         bn (i100) = (cfac2 * rhelp1 (1) - rhelp2 (1)) &
                   / (cfac2 * chelp1     - chelp2    )
!*
         chelp1 = an (i100)
         chelp2 = bn (i100)
!*
         rfac1 = rhelp1 (1) + rhelp2 (1)
         rfac2 = cdabs (an (i100)) * cdabs (an (i100)) &
               + cdabs (bn (i100)) * cdabs (bn (i100))
!*
         qexti  = qexti  + (2.0d0 * i100 + 1.0d0) * rfac1
         qscat  = qscat  + (2.0d0 * i100 + 1.0d0) * rfac2
         cbscat = cbscat + (2.0d0 * i100 + 1.0d0) * (-1.0d0) ** i100 &
                * (an (i100) - bn (i100))
!*
         if (i100 .eq. 1) go to 100
!*
         chelp1 = an (i100 - 1) * conjg (an (i100)) &
                + bn (i100 - 1) * conjg (bn (i100))
         chelp2 = an (i100 - 1) * conjg (bn (i100 - 1))
!*
         i101 = i100 - 1
         rfac1  = i101 * (i101 + 2) / (i101 + 1.0d0)
         rfac2  = (2.0d0 * i101 + 1.0d0) / (i101 * (i101 + 1.0d0))
!*
         asym = asym + rfac1 * rhelp1 (1) + rfac2 * rhelp2 (1)
100   continue
!*
      qexti  = qexti * 2.0d0 / (x * x)
      qscat  = qscat * 2.0d0 / (x * x)
      asym   = asym  * 4.0d0 / (x * x * qscat)
      qbscat = cdabs (cbscat) * cdabs (cbscat) / (x * x)
      if (qscat .gt. qexti) qscat = qexti
!*
   return
 end subroutine mie_sphere

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gas_absorb(freq,temp,pres,relhum,water,kabs)
 implicit none

!----IO parameter----
      real(sdsu_fps), intent(in) :: freq   !frequency [GHz]
      real(sdsu_fps), intent(in) :: temp   !temperature [K]
      real(sdsu_fps), intent(in) :: relhum !relative humidity [%]
      real(sdsu_fps), intent(in) :: pres   !pressure [mb] 
      real(sdsu_fps), intent(in) :: water  !water content [g/m3]
      real(sdsu_fps), intent(out) :: kabs  !absorption coefficient
   
!----local parameter-------
      real(sdsu_fps),parameter:: a0=6.107799961e0
      real(sdsu_fps),parameter:: a1=4.436518521e-1
      real(sdsu_fps),parameter:: a2=1.428945805e-2
      real(sdsu_fps),parameter:: a3=2.650648471e-4
      real(sdsu_fps),parameter:: a4=3.031240396e-6
      real(sdsu_fps),parameter:: a5=2.034080948e-8
      real(sdsu_fps),parameter:: a6=6.136820929e-11

      real(sdsu_fps) :: tc !temperature [DegK]
      real(sdsu_fps) :: es !saturation papor pressure
      real(sdsu_fps) :: vapor_pressure !vapor pressure 
      real(sdsu_fps) :: rho            !water vapor density [g/m3] 
      real(sdsu_fps) :: kabs_h2o , kabs_o2 , kabs_n2 , kabs_clw  !absorption coefficient in each component. 

      tc = temp - 273.15
      es = a0+tc*(a1+tc*(a2+tc*(a3+tc*(a4+tc*(a5+a6*tc)))))
      if (es .lt. 0.) es = 0.
      vapor_pressure =  relhum*es/100.
      rho = vapor_pressure*100.*18./(8.314*temp)

      call abs_h2o(temp, pres, rho, freq, kabs_h2o)

      call abs_o2(temp, pres, rho, freq, kabs_o2)

      call abs_n2(temp, pres, freq, kabs_n2)

      call absorb_clw(freq,temp,water,kabs_clw)

      kabs = kabs_h2o + kabs_o2 + kabs_n2 + kabs_clw

 return
 end subroutine gas_absorb

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine abs_h2o(t,p,rho,f,abh2o)
      implicit none
!
!     C. Kummerow, 8/2003.  Changed function to subroutine
!     Copyright (c) 2002 Massachusetts Institute of Technology
!
!  NAME- ABH2O    LANGUAGE- FORTRAN 77
!
! PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
!
!  CALLING SEQUENCE PARAMETERS-
!    SPECIFICATIONS
      real(sdsu_fps) :: t,p,rho,f,abh2o
!      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
!      T       KELVIN    I   TEMPERATURE
!      P       MILLIBAR  I   PRESSURE              .1 TO 1000
!      RHO     G/M**3    I   WATER VAPOR DENSITY
!      F       GHZ       I   FREQUENCY             0 TO 800
!      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT
!
!   REFERENCES-
!   P.W. ROSENKRANZ, RADIO SCIENCE V.33, PP.919-928 (1998); V.34, P.1025 (1999).
!
!   LINE INTENSITIES SELECTION THRESHOLD=
!     HALF OF CONTINUUM ABSORPTION AT 1000 MB.
!   WIDTHS MEASURED AT 22,183,380 GHZ, OTHERS CALCULATED.
!     A.BAUER ET AL.ASA WORKSHOP (SEPT. 1989) (380GHz).
!     M. TRETYAKOV et al., J. MOLEC. SPEC. (2003)
!
!   REVISION HISTORY-
!    DATE- OCT.6, 1988  P.W.ROSENKRANZ - EQS AS PUBL. IN 1993.
!          OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
!                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
!          OCT. 24, 95  PWR -ADD 1 LINE.
!          JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING,
!                       REVISED CONTINUUM.
!        Aug. 28, 2002  PWR - CORRECTED LINE INTENSITIES
!        Mar. 2, 2003   PWR - LINE SHIFT
!
!   LOCAL VARIABLES:
      integer nlines,i,j
      parameter (nlines=15)
      real(sdsu_fps) :: df(2),s1(nlines),b2(nlines),w3(nlines),fl(nlines),x(nlines), &
       ws(nlines),xs(nlines),sr(nlines)
      real(sdsu_fps) :: pvap,pda,den,ti,ti2,sum,width,wsq,s,base,res,con,shift
!     LINE FREQUENCIES:
      data fl/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, &
       443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360, &
       620.7008, 752.0332, 916.1712/
!     LINE INTENSITIES AT 300K:
      data s1/ .1314e-13, .2279e-11, .8058e-13, .2701e-11, .2444e-10, &
       .2185e-11, .4637e-12, .2568e-10, .8392e-12, .3272e-11, .6676e-12, &
       .1535e-08, .1711e-10, .1014e-08, .4238e-10/
!     T COEFF. OF INTENSITIES:
      data b2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, &
       3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
!     AIR-BROADENED WIDTH PARAMETERS AT 300K:
      data w3/.00281, .00287, .0023, .00278, .00287, .0021, .00186, &
       .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/
!     T-EXPONENT OF AIR-BROADENING:
      data x/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69, &
       .71, .68, .70/
!     SELF-BROADENED WIDTH PARAMETERS AT 300K:
      data ws/.01349, .01491, .0108, .0135, .01541, .0090, .00788, &
       .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/
!     T-EXPONENT OF SELF-BROADENING:
      data xs/ .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72, &
       1.0, .68, .84, .78/
!     RATIO OF SHIFT TO WIDTH
      data sr/ 0., -.017, 13*0./
!
      if(rho.le.0.) then
        abh2o = 0.
        return
      endif
      pvap = rho*t/217.
      pda = p -pvap
      den = 3.335e16*rho ! const includes isotopic abundance
      ti = 300./t
      ti2 = ti**2.5
!
!      CONTINUUM TERMS
      con = (5.43e-10*pda*ti**3 + 1.8e-8*pvap*ti**7.5)*pvap*f*f
!
!      ADD RESONANCES
      sum = 0.
      do 30 i=1,nlines
      width = w3(i)*pda*ti**x(i) + ws(i)*pvap*ti**xs(i)
      shift = sr(i)*width  ! unknown temperature dependence
      wsq = width*width
      s = s1(i)*ti2*exp(b2(i)*(1.-ti))
      df(1) = f - fl(i) - shift
      df(2) = f + fl(i) + shift
!  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
      base = width/(562500. + wsq)
!  DO FOR POSITIVE AND NEGATIVE RESONANCES
      res = 0.
      do 20 j=1,2
      if(abs(df(j)).lt.750.) res = res + width/(df(j)**2+wsq) - base
20    continue
30    sum = sum + s*res*(f/fl(i))**2
      abh2o = .3183e-4*den*sum + con
      return
      end subroutine abs_h2o

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine abs_o2(temp,pres,vapden,freq,o2abs)
      implicit none
!
!     C. Kummerow, 8/2003.  Changed function to subroutine
!  Copyright (c) 2003 Massachusetts Institute of Technology
!
!     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
!              IN NEPERS/KM
!
!      5/1/95  P. Rosenkranz
!      11/5/97  P. Rosenkranz - 1- line modification.
!      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
!      8/21/02  pwr - revised width at 425
!      3/20/03  pwr - 1- line mixing and width revised
!
!     IMPLICIT NONE
!
!     ARGUMENTS:
    real(sdsu_fps) :: temp,pres,vapden,freq
    real(sdsu_fps),intent(out) :: o2abs
!
!     NAME    UNITS    DESCRIPTION        VALID RANGE
!
!     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
!                                          valid for atmosphere
!     PRES   MILLIBARS PRESSURE           3 TO 1000
!     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
!                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
!     FREQ    GHZ      FREQUENCY          0 TO 900
!
!     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
!     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
!      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
!     H.J. Liebe et al, JQSRT V.48, pp.629-643 (1992).
!     M.J. Schwartz, Ph.D. thesis, M.I.T. (1998).
!     A.F. Krupnov et al, J. Mol. Spect. v.215, pp.309-311 (2002).
!     M.Yu. Tretyakov et al, J. Mol. Spect. (2003 preprint).
!     SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
!
!     This version differs from Liebe's MPM92 in these significant respects:
!     1. The 1- line has the width and mixing coefficient measured by
!      Tretyakov et al.
!     2. It modifies the 1- line width temperature dependence to (1/T)**0.9
!     3. It uses the same temperature dependence (X) for submillimeter
!      line widths as in the 60 GHz band: (1/T)**0.8
!     4. The 425 GHz line width is from Krupnov et al.
!
!     Local variables:
      real(sdsu_fps) :: th,th1,b,preswv,presda,den,dens,dfnr,sum,str,y,sf1,sf2,fcen,df
      integer k
      real(sdsu_fps) :: x,wb300,w300(40),f(40),y300(40),s300(40),v(40),be(40)
      common /o2com/ x,wb300,w300,f,y300,s300,v,be
!      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      data f/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910, &
        59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002, &
        56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685, &
        55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241, &
        53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368, &
        52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7632, &
        487.2494, 715.3931, 773.8397, 834.1458/
        data s300/.2936e-14,.8079e-15, .2480e-14,.2228e-14, &
        .3351e-14,.3292e-14, .3721e-14,.3891e-14, &
        .3640e-14,.4005e-14, .3227e-14,.3715e-14, &
        .2627e-14,.3156e-14, .1982e-14,.2477e-14, &
        .1391e-14,.1808e-14, .9124e-15,.1230e-14, &
        .5603e-15,.7842e-15, .3228e-15,.4689e-15, &
        .1748e-15,.2632e-15, .8898e-16,.1389e-15, &
        .4264e-16,.6899e-16, .1924e-16,.3229e-16, &
        .8191e-17,.1423e-16, .6494e-15, .7083e-14, .3025e-14, &
        .1835e-14, .1158e-13, .3993e-14/
      data be/.009,.015, .083,.084, 2*.212, 2*.391, 2*.626, &
       2*.915, 2*1.260, 1.660,1.665, 2.119,2.115, 2.624,2.625, &
       2*3.194, 2*3.814, 2*4.484, 2*5.224, 2*6.004, 2*6.844, &
       2*7.744, .048, .044, .049, .145, .141, .145/
!      WIDTHS IN MHZ/MB
      data wb300/.56/, x/.8/
      data w300/1.67, 1.646, 1.468, 1.449, 1.382, 1.360, &
       1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171, &
       1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 2*1.05, &
       2*1.02,2*1.00,2*.97,2*.94,2*.92,2*.89, 3*1.64, 3*1.81/
      data y300/  -0.036,  0.2408, -0.3486,  0.5227, &
       -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311, &
        0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599, &
        0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246, &
        0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546, &
        0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, 6*0./
      data v/  0.0079, -0.0978,  0.0844, -0.1273, &
        0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584, &
        0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675, &
        0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590, &
        0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091, &
        0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, 6*0./
!
      th = 300./temp
      th1 = th-1.
      b = th**x
      preswv = vapden*temp/217.
      presda = pres -preswv
      den = .001*(presda*b + 1.1*preswv*th)
      dens = .001*(presda*th**.9 + 1.1*preswv*th)
      dfnr = wb300*den
      sum = 1.6e-17*freq*freq*dfnr/(th*(freq*freq + dfnr*dfnr))
      do 32 k=1,40
      if(k.eq.1) then !exception for 1- line
        df = w300(1)*dens
      else
        df = w300(k)*den
      endif
      fcen = f(k)
      y = .001*pres*b*(y300(k)+v(k)*th1)
      str = s300(k)*exp(-be(k)*th1)
      sf1 = (df + (freq-fcen)*y)/((freq-fcen)**2 + df*df)
      sf2 = (df - (freq+fcen)*y)/((freq+fcen)**2 + df*df)
32    sum = sum + str*(sf1+sf2)*(freq/f(k))**2
      o2abs = .5034e12*sum*presda*th**3/3.14159
      o2abs = amax1(o2abs,0.)
      return
      end subroutine abs_o2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

   subroutine abs_n2(t,p,f,absn2)
   implicit none

   real(sdsu_fps),intent(in) :: t ! temperature [K]
   real(sdsu_fps),intent(in) :: p ! pressure [mb]
   real(sdsu_fps),intent(in) :: f ! frequency [GHz] (valid 0 - 1000GHz)
   real(sdsu_fps),intent(out) :: absn2 
   real(sdsu_fps) :: fdepen, bf, th

!
!     C. Kummerow, 8/2003.  Changed function to subroutine
!  Copyright (c) 2002 Massachusetts Institute of Technology
!     ABSN2 = COLLISION-INDUCED ABSORPTION COEFFICIENT (NEPER/KM)
!     5/22/02 P.Rosenkranz
!
!     Equations based on:
!      Borysow, A, and L. Frommhold,
!      Astrophysical Journal, v.311, pp.1043-1057 (1986)
!     with modification of 1.29 to account for O2-O2 and O2-N2
!     collisions, as suggested by
!      J.R. Pardo, E.Serabyn, J.Cernicharo, J. Quant. Spectros.
!      Radiat. Trans. v.68, pp.419-433 (2001).
!

      th = 300./t

      fdepen = .5 + .5/(1.+(f/450.)**2)

      bf = 6.5e-14*fdepen*p*p*f*f*th**3.6

      absn2 = 1.29*bf

      return

      end subroutine abs_n2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine absorb_clw(freq,temp,water,kabs_clw)
      implicit none

      real(sdsu_fps),intent(in) :: freq !frequency [GHz]
      real(sdsu_fps),intent(in) :: temp !temperature [K]
      real(sdsu_fps),intent(in) :: water !water content [g/m3]

!     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
!     ARGUMENTS (INPUT):
!     WATER IN G/M**3
!     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
!     TEMP IN KELVIN

      real(sdsu_fps),intent(out) :: kabs_clw
!
!     REFERENCES:
!     LIEBE, HUFFORD AND MANABE, INT. J. IR & MM WAVES V.12, pp.659-675
!      (1991);  Liebe et al, AGARD Conf. Proc. 542, May 1993.
!
!     REVISION HISTORY:
!        PWR 8/3/92   original version
!        PWR 12/14/98 temp. dependence of EPS2 eliminated to agree
!                     with MPM93
!        pwr 2/27/02  use exponential dep. on T, eq. 2b instead of eq. 4a
!
    real(sdsu_fps) :: eps0, eps1, eps2
    real(sdsu_fps) :: theta1, fp, fs
    complex eps,re
      if(water.le.0.) then
       kabs_clw = 0.
       return
      endif
      theta1 = 1.-300./temp
      eps0 = 77.66 - 103.3*theta1
      eps1 = .0671*eps0
      eps2 = 3.52                 ! from MPM93
      fp = 20.1*exp(7.88*theta1)  ! from eq. 2b
      fs = 39.8*fp
      eps = (eps0-eps1)/cmplx(1.,freq/fp) + &
       (eps1-eps2)/cmplx(1.,freq/fs) +eps2
      re = (eps-1.)/(eps+2.)
      kabs_clw = -.06286*aimag(re)*freq*water
      return
      end subroutine absorb_clw

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine iceoptic (fghz, t, epsri, epsii)
      implicit none
!
!     Compute  complex refractive index of cloud ice
!
!  Hufford (1991), see Brussard and Watson (1995), p.297
!  The radio frequency may range from 0 to 1000 GHz, and the temperature from 0° to –40°C.
!
      save

      real(sdsu_fpd) :: fghz, epsri, epsii, a, b, theta, t, tk
!
      epsri = 3.15d+00

      tk = t
      if (tk .gt. 273.16) tk = 273.16  !always below freezing 

      theta = 300.0 / tk
      a     = 1.0d-04 * (50.4 + 62.0 * (theta - 1.0)) &
            * exp (-22.1 * (theta - 1.0))
      b     = 1.0d-04 * (0.633 / theta - 0.131) &
            + (7.36d-04 * theta / (theta - 0.9927)) &
            * (7.36d-04 * theta / (theta - 0.9927))
      epsii = a / fghz + b * fghz
!
      return
      end subroutine iceoptic 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine watoptic (freqhz, ktemp, salinity, epsreal, epsimag)
      implicit none
!
! complex refractive index (permittivity) of liquid water
!
      save

      real(sdsu_fpd) :: fac1
      real(sdsu_fpd) :: conv, pi, ghztohz
      real(sdsu_fpd) :: ktemp, ctemp, salinity
      real(sdsu_fpd) :: eps0, epsstat, epshigh, epsreal, epsimag
      real(sdsu_fpd) :: trelax, freqhz, omega
      data pi      /  3.141592654d+00 /
      data conv    /  0.017453292d+00 /
      data ghztohz /  1.0d+09 /
      data eps0     / 8.854d-12 /  !permittivity in vacume
      data epshigh  / 4.9d+00 /  !high-frequency permittivity


!*
      ctemp = ktemp - 273.15d+00
      omega = (pi + pi) * freqhz
!
! parameterized static permittivity (v->0)
!
      epsstat = (87.134d+00 &
              - 1.949d-01 * ctemp &
              - 1.276d-02 * ctemp * ctemp &
              + 2.491d-04 * ctemp * ctemp * ctemp) &
              * (1.0d+00 &
              + 1.613d-05 * salinity * ctemp &
              - 3.656d-03 * salinity &
              + 3.210d-05 * salinity * salinity &
              - 4.232d-07 * salinity * salinity * salinity)
!
! parameterized effective relaxation time
!
      trelax = (1.768d-11 &
             - 6.086d-13 * ctemp &
             + 1.104d-14 * ctemp * ctemp &
             - 8.111d-17 * ctemp * ctemp * ctemp) &
             * (1.0d+00 &
             + 2.282d-05  * salinity * ctemp &
             - 7.638d-04  * salinity &
             - 7.760d-06  * salinity * salinity &
             + 1.105d-08  * salinity * salinity * salinity)

      fac1    = 1.0d+00 + omega * omega * trelax * trelax

!
! real and imaginery components of permittivity
!
      epsreal = epshigh + (epsstat - epshigh) / fac1
      epsimag = ((epsstat - epshigh) * omega * trelax) / fac1

      return
      end subroutine watoptic

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine total_opt_gen(press, w, atm_ext, kext_gen, salb_gen, asym_gen, pbck_gen, vdop_gen, &
                          opt )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute total optical properties by integrating multiple hydrometeors. 
! 
! History:
! 10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: press  !atmosphere pressure [mb]
 real(sdsu_fps),intent(in) :: w  !atmosphere vertical velocity [m/s]
 real(sdsu_fps),intent(in) :: atm_ext  !atmosphere extinction 
 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     kext_gen,  & ! extinction coefficient [km-1]
     salb_gen,  & ! single scattering albedo [-]
     asym_gen,  & ! asymetry parameter [-]
     pbck_gen,  & ! backscattering coefficient [-]
     vdop_gen     ! Doppler velocity [m/s]

 type(opt_micro_parameter),intent(out) :: opt

 real(sdsu_fps) :: vdop_raw

!
! total extinction 
!
   opt%kext    = atm_ext         + &
                kext_gen%cloud  + & 
                kext_gen%rain   + &
                kext_gen%ice    + &
                kext_gen%snow   + &
                kext_gen%graupel+ &
                kext_gen%hail 

!
! total single scattering albedo
!
   if ( opt%kext <= 0.e0 ) then
        opt%salb = 0.e0
   else
       opt%salb = ( &
                  salb_gen%cloud  * kext_gen%cloud   + &
                  salb_gen%rain   * kext_gen%rain    + &
                  salb_gen%ice    * kext_gen%ice     + &
                  salb_gen%snow   * kext_gen%snow    + &
                  salb_gen%graupel* kext_gen%graupel + &
                  salb_gen%hail   * kext_gen%hail      &
                 ) / opt%kext

   endif

!
! total asymetry parameter
!
   if ( opt%salb <= 0.e0 ) then
       opt%asym = 0.e0
   else
       opt%asym = ( &
                  asym_gen%cloud  * salb_gen%cloud  * kext_gen%cloud   + &
                  asym_gen%rain   * salb_gen%rain   * kext_gen%rain    + &
                  asym_gen%ice    * salb_gen%ice    * kext_gen%ice     + &
                  asym_gen%snow   * salb_gen%snow   * kext_gen%snow    + &
                  asym_gen%graupel* salb_gen%graupel* kext_gen%graupel + &
                  asym_gen%hail   * salb_gen%hail   * kext_gen%hail      &
                 ) / ( opt%salb * opt%kext )

   endif

!
! total backscattering coef. 
!
   opt%sback = &
                  pbck_gen%cloud  * salb_gen%cloud  * kext_gen%cloud   + &
                  pbck_gen%rain   * salb_gen%rain   * kext_gen%rain    + &
                  pbck_gen%ice    * salb_gen%ice    * kext_gen%ice     + &
                  pbck_gen%snow   * salb_gen%snow   * kext_gen%snow    + &
                  pbck_gen%graupel* salb_gen%graupel* kext_gen%graupel + &
                  pbck_gen%hail   * salb_gen%hail   * kext_gen%hail      


!
! Doppler velocity of hydrometeors [m/s]
!
   if ( opt%sback <= 0.e0 ) then
     vdop_raw = 0.e0
   else
     vdop_raw = &
                ( pbck_gen%cloud  * salb_gen%cloud  * kext_gen%cloud  * vdop_gen%cloud  + &
                  pbck_gen%rain   * salb_gen%rain   * kext_gen%rain   * vdop_gen%rain   + &
                  pbck_gen%ice    * salb_gen%ice    * kext_gen%ice    * vdop_gen%ice    + &
                  pbck_gen%snow   * salb_gen%snow   * kext_gen%snow   * vdop_gen%snow   + &
                  pbck_gen%graupel* salb_gen%graupel* kext_gen%graupel* vdop_gen%graupel+ &
                  pbck_gen%hail   * salb_gen%hail   * kext_gen%hail   * vdop_gen%hail  ) / &
                  opt%sback
   endif

!
! correct it with pressure and actual atmospheric w (downward is positive)
!
  opt%vdop  = vdop_raw*sqrt(800.e0/press)  - w ! actual terminal velocity for a given pressure [m/s]


 return
 end subroutine total_opt_gen

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine total_opt_sbm( press, w, &
                           atm_ext, kext_sbm, salb_sbm, asym_sbm, pbck_sbm, vdop_sbm , &
                           opt )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute total optical properties by integrating multiple hydrometeors. 
! 
! History:
! 10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: press  !atmosphere pressure [mb]
 real(sdsu_fps),intent(in) :: w   !atmosphere vertical velocity [m/s]
 real(sdsu_fps),intent(in) :: atm_ext  !atmosphere extinction 
 type ( particle_sbm ),intent(in) :: & ! particle_sbm is defined in module_simulater
     kext_sbm,  & ! extinction coefficient
     salb_sbm,  & ! single scattering albedo [-]
     asym_sbm,  & ! asymetry parameter [-]
     pbck_sbm,  & ! backscattering coefficient
     vdop_sbm     ! Doppler velocity coefficient [m/s]
 
 type(opt_micro_parameter),intent(out) :: opt

 real(sdsu_fps) :: vdop_raw


!
! total extinction 
!
   opt%kext    = atm_ext          + &
                kext_sbm%liq    + & 
                kext_sbm%ice_col+ &
                kext_sbm%ice_pla+ &
                kext_sbm%ice_den+ &
                kext_sbm%snow   + &
                kext_sbm%graupel+ &
                kext_sbm%hail 

!
! total single scattering albedo
!
   if ( opt%kext <= 0.e0 ) then
        opt%salb = 0.e0
   else
       opt%salb = ( &
                  salb_sbm%liq    * kext_sbm%liq     + &
                  salb_sbm%ice_col* kext_sbm%ice_col + &
                  salb_sbm%ice_pla* kext_sbm%ice_pla + &
                  salb_sbm%ice_den* kext_sbm%ice_den + &
                  salb_sbm%snow   * kext_sbm%snow    + &
                  salb_sbm%graupel* kext_sbm%graupel + &
                  salb_sbm%hail   * kext_sbm%hail      &
                 ) / opt%kext

   endif

!
! total asymetry parameter
!
   if ( opt%salb <= 0.e0 ) then
       opt%asym = 0.e0
   else
       opt%asym = ( &
                  asym_sbm%liq    * salb_sbm%liq    * kext_sbm%liq     + &
                  asym_sbm%ice_col* salb_sbm%ice_col* kext_sbm%ice_col + &
                  asym_sbm%ice_pla* salb_sbm%ice_pla* kext_sbm%ice_pla + &
                  asym_sbm%ice_den* salb_sbm%ice_den* kext_sbm%ice_den + &
                  asym_sbm%snow   * salb_sbm%snow   * kext_sbm%snow    + &
                  asym_sbm%graupel* salb_sbm%graupel* kext_sbm%graupel + &
                  asym_sbm%hail   * salb_sbm%hail   * kext_sbm%hail      &
                 ) / ( opt%salb * opt%kext )

   endif

!
! total backscattering coef. 
!
   opt%sback = &
                  pbck_sbm%liq    * salb_sbm%liq    * kext_sbm%liq     + &
                  pbck_sbm%ice_col* salb_sbm%ice_col* kext_sbm%ice_col + &
                  pbck_sbm%ice_pla* salb_sbm%ice_pla* kext_sbm%ice_pla + &
                  pbck_sbm%ice_den* salb_sbm%ice_den* kext_sbm%ice_den + &
                  pbck_sbm%snow   * salb_sbm%snow   * kext_sbm%snow    + &
                  pbck_sbm%graupel* salb_sbm%graupel* kext_sbm%graupel + &
                  pbck_sbm%hail   * salb_sbm%hail   * kext_sbm%hail      


!
! Doppler velocity of hydrometeors [m/s]
!
   if ( opt%sback <= 0.e0 ) then
     vdop_raw = 0.e0
   else
     vdop_raw = &
                ( pbck_sbm%liq    * salb_sbm%liq    * kext_sbm%liq     * vdop_sbm%liq     + &
                  pbck_sbm%ice_col* salb_sbm%ice_col* kext_sbm%ice_col * vdop_sbm%ice_col + &
                  pbck_sbm%ice_pla* salb_sbm%ice_pla* kext_sbm%ice_pla * vdop_sbm%ice_pla + &
                  pbck_sbm%ice_den* salb_sbm%ice_den* kext_sbm%ice_den * vdop_sbm%ice_den + &
                  pbck_sbm%snow   * salb_sbm%snow   * kext_sbm%snow    * vdop_sbm%snow    + &
                  pbck_sbm%graupel* salb_sbm%graupel* kext_sbm%graupel * vdop_sbm%graupel + &
                  pbck_sbm%hail   * salb_sbm%hail   * kext_sbm%hail    * vdop_sbm%hail ) / &
                  opt%sback
   endif

!
! correct it with pressure and actual atmospheric w (downward is positive)
! (w: upward is positive)
  opt%vdop  = vdop_raw *sqrt(1000.e0/press) - w ! measured Doppler velocity [m/s]


 return
 end subroutine total_opt_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine total_opt_rams(rho_dair, w, atm_ext, kext_rams, salb_rams, asym_rams, pbck_rams, vdop_rams, &
                           opt )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute total optical properties by integrating multiple hydrometeors. 
! 
! History:
! 10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: rho_dair  !atmosphere dry air density  [kg/m3]
 real(sdsu_fps),intent(in) :: w  !atmosphere vertical velocity [m/s]
 real(sdsu_fps),intent(in) :: atm_ext  !atmosphere extinction 
 type ( particle_rams ),intent(in) :: & ! particle_rams is defined in module_simulater
     kext_rams,  & ! extinction coefficient
     salb_rams,  & ! single scattering albedo [-]
     asym_rams,  & ! asymetry parameter [-]
     pbck_rams,  & ! backscattering coefficient
     vdop_rams     ! Doppler velocity [m/s]

 type(opt_micro_parameter),intent(out) :: opt

 real(sdsu_fps) :: vdop_raw 

!
! total extinction 
!
   opt%kext    = atm_ext          + &
                kext_rams%cloud1 + & 
                kext_rams%cloud2 + &
                kext_rams%rain   + &
                kext_rams%ice1   + &
                kext_rams%ice2   + &
                kext_rams%snow   + &
                kext_rams%graupel+ &
                kext_rams%hail 

!
! total single scattering albedo
!
   if ( opt%kext <= 0.e0 ) then
        opt%salb = 0.e0
   else
       opt%salb = ( &
                  salb_rams%cloud1 * kext_rams%cloud1  + &
                  salb_rams%cloud2 * kext_rams%cloud2  + &
                  salb_rams%rain   * kext_rams%rain    + &
                  salb_rams%ice1   * kext_rams%ice1    + &
                  salb_rams%ice2   * kext_rams%ice2    + &
                  salb_rams%snow   * kext_rams%snow    + &
                  salb_rams%graupel* kext_rams%graupel + &
                  salb_rams%hail   * kext_rams%hail      &
                 ) / opt%kext

   endif

!
! total asymetry parameter
!
   if ( opt%salb <= 0.e0 ) then
       opt%asym = 0.e0
   else
       opt%asym = ( &
                  asym_rams%cloud1 * salb_rams%cloud1 * kext_rams%cloud1  + &
                  asym_rams%cloud2 * salb_rams%cloud2 * kext_rams%cloud2  + &
                  asym_rams%rain   * salb_rams%rain   * kext_rams%rain    + &
                  asym_rams%ice1   * salb_rams%ice1   * kext_rams%ice1    + &
                  asym_rams%ice2   * salb_rams%ice2   * kext_rams%ice2    + &
                  asym_rams%snow   * salb_rams%snow   * kext_rams%snow    + &
                  asym_rams%graupel* salb_rams%graupel* kext_rams%graupel + &
                  asym_rams%hail   * salb_rams%hail   * kext_rams%hail      &
                 ) / ( opt%salb * opt%kext )

   endif

!
! total backscattering coef. 
!
   opt%sback = &
                  pbck_rams%cloud1 * salb_rams%cloud1 * kext_rams%cloud1  + &
                  pbck_rams%cloud2 * salb_rams%cloud2 * kext_rams%cloud2  + &
                  pbck_rams%rain   * salb_rams%rain   * kext_rams%rain    + &
                  pbck_rams%ice1   * salb_rams%ice1   * kext_rams%ice1    + &
                  pbck_rams%ice2   * salb_rams%ice2   * kext_rams%ice2    + &
                  pbck_rams%snow   * salb_rams%snow   * kext_rams%snow    + &
                  pbck_rams%graupel* salb_rams%graupel* kext_rams%graupel + &
                  pbck_rams%hail   * salb_rams%hail   * kext_rams%hail      


!
! Doppler velocity of hydrometeors [cm/s]
!
   if ( opt%sback <= 0.e0 ) then
     vdop_raw = 0.e0
   else
     vdop_raw = &
                 (pbck_rams%cloud1 * salb_rams%cloud1 * kext_rams%cloud1  * vdop_rams%cloud1 + &
                  pbck_rams%cloud2 * salb_rams%cloud2 * kext_rams%cloud2  * vdop_rams%cloud2 + &
                  pbck_rams%rain   * salb_rams%rain   * kext_rams%rain    * vdop_rams%rain   + &
                  pbck_rams%ice1   * salb_rams%ice1   * kext_rams%ice1    * vdop_rams%ice1   + &
                  pbck_rams%ice2   * salb_rams%ice2   * kext_rams%ice2    * vdop_rams%ice2   + &
                  pbck_rams%snow   * salb_rams%snow   * kext_rams%snow    * vdop_rams%snow   + &
                  pbck_rams%graupel* salb_rams%graupel* kext_rams%graupel * vdop_rams%graupel+ &
                  pbck_rams%hail   * salb_rams%hail   * kext_rams%hail    * vdop_rams%hail ) / &
                  opt%sback
   endif 

!
! correct it with pressure and actual atmospheric w (downward is positive)
!
  opt%vdop  = vdop_raw*sqrt(0.7/rho_dair) - w  ! actual terminal velocity for a given air density[m/s]


 return
 end subroutine total_opt_rams

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_opt_micro
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out microwave single-scattering properties in GrADS format. 
! 
! History:
! 10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer i,j,k,l,m,irec
 integer,parameter :: io = 301
 integer :: is,ie,js,je,ks,ke
 character(len=3) :: fff

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end

! modify ke for slant path 
 if( scan_micro .or. slant_path_micro ) then
   ks = 1 ; ke = mxlyr_micro
 endif

!
! Computed Tbs in binary format
!
 sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.OPT.'//&
                trim(micro_sensor)//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG write_opt_micro ; output for microwave optical properties ',&
 ' ->',trim(sdsu_io_file)

!
! Write out binary data 
!
#if MPI == 2
 if(masterproc) &
#endif
 open( unit= io, file = sdsu_io_file, access='direct', &
       status='replace',recl=(ke-ks+1)*mxgridx*mxgridy*4 )
 irec = 1

 do l = 1,mxfreq_micro
    call dump( atmos_opt(is:ie,js:je,ks:ke,l)%kext , io, irec ) !total extinction coef [km-1]
    call dump( atmos_opt(is:ie,js:je,ks:ke,l)%salb , io, irec ) !total single scatterling albedo [-]
    call dump( atmos_opt(is:ie,js:je,ks:ke,l)%asym , io, irec ) !total asymetery parameters [-]
 enddo ! l

#if MPI == 2
 if(masterproc) &
#endif
 close(io)


!
!  write out grads control file
!
#if MPI == 2
 if(masterproc) then
#endif

 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.OPT.'//&
                  trim(micro_sensor)//trim(output_suffix)//'.ctl'
    print*,'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))

    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.OPT.'//trim(micro_sensor)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(micro_sensor) //' microwave optical properties'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize

    if( scan_micro .or. slant_path_micro ) then
       write(io,*) 'zdef   ',(ke-ks+1) ,' linear 0.',range_micro
    else
       write(io,*) 'zdef   ',(ke-ks+1) ,' linear 0. 1 '  !toshii fix this??
    endif

    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxfreq_micro*3

    do l=1,mxfreq_micro !frequency loop

       write(fff,"(I3.3)") INT(freq_micro(l))
       write(io,*) 'kext'//fff,(ke-ks+1) ,'0  Extinction '//trim(nch_micro(l))//' [km-1]'
       write(io,*) 'salb'//fff,(ke-ks+1) ,'0  SS albedo  '//trim(nch_micro(l))//' [-]'
       write(io,*) 'asym'//fff,(ke-ks+1) ,'0  Asymetry   '//trim(nch_micro(l))//' [-]'

    enddo

    write(io,*) 'endvars '
    close(io)
 endif

#if MPI == 2
 endif
#endif

 return
 end subroutine write_opt_micro


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_opt_micro_mpi
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out microwave single-scattering properties in GrADS format. 
! 
! History:
! 10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer i,j,k,l,m,irec
 integer,parameter :: io = 301
 integer :: is,ie,js,je,ks,ke
 character (len=2) :: crank
!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end

!
! Computed Tbs in binary format
!
  write(crank,"(I2.2)") myrank

 sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.OPT.'//crank//&
                trim(micro_sensor)//trim(output_suffix)//'.bin'

 print*,'MSG write_opt_micro ; output for microwave optical properties ',&
 ' ->',trim(sdsu_io_file)

!
! Write out binary data 
!
 open( unit= io, file = sdsu_io_file, access='direct', &
       status='replace',recl=mxgridx*mxgridy*4 )
 irec = 1

 do l = 1,mxfreq_micro
     write(io,rec=irec) ((atmos_opt(i,j,1,l)%kext,i=1,mxgridx),j=1,mxgridy)  !microwave Tb
 enddo ! l

 close(io)

 return
 end subroutine write_opt_micro_mpi

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_opt_radar
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out radar backscatter in GrADS format. 
! 
! History:
! 10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,k,l,m,irec
 integer,parameter :: io = 302
 integer :: is,ie,js,je,ks,ke
 character(len=3) :: fff

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end 
 ks = 1 ; ke = mxlyr_radar


!
! Computed Tbs in binary format
!
 sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.OPT.'//&
                trim(radar_sensor)//trim(output_suffix)//'.bin'
#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG write_opt_radar ; output for radar optical properties -->'&
 ,trim(sdsu_io_file)

!
! Write out binary data 
!
#if MPI == 2
 if(masterproc) &
#endif
 open( unit= io, file = sdsu_io_file, access='direct', &
       status='replace',recl=mxlyr*mxgridx*mxgridy*4)
 irec = 1

 do l = 1,mxfreq_radar
    call dump( atmos_opt(is:ie,js:je,ks:ke,l)%sback , io, irec )  !total backscatter [km-1]
    call dump( atmos_opt(is:ie,js:je,ks:ke,l)%vdop  , io, irec )  !Doppler velocity [m/s]

 enddo ! l
 
#if MPI == 2
 if(masterproc) &
#endif
 close(io)

!
! write control file
!

#if MPI == 2
 if(masterproc) then
#endif

 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.OPT.'//&
                  trim(radar_sensor)//trim(output_suffix)//'.ctl'
    print*,'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))

    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.OPT.'//trim(radar_sensor)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(radar_sensor) //' radar optical properties'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef   ',(ke-ks+1) ,' linear 0.',range_radar
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxfreq_radar*2

    do l=1,mxfreq_radar !frequency loop

       write(fff,"(I3.3)") INT(freq_radar(l))
       write(io,*) 'sback'//fff,(ke-ks+1) ,'0  Radar backscatter '//trim(nch_radar(l))//' [km-1]'
       write(io,*) 'vdop'//fff,(ke-ks+1) ,'0  Doppler velocity '//trim(nch_radar(l))//' [m/s]'

    enddo

    write(io,*) 'endvars '
    close(io)
 endif

#if MPI == 2
 endif
#endif


 return
 end subroutine write_opt_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine L2_radar_parameters_gen(mu, q, re, a, L2_radar_0d )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Write out various momemnt from SBM inputs for studying particle-size distributions.  
! 
! History:
!  03/2012  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     mu, & ! mu of gamma DSD [-]
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]
 type ( atmos_parameter ) ,intent(in) :: a
 type(particle_retrieval),intent(out) :: L2_radar_0d  !Level 2 radar parameters

 integer :: ierr,irec
 integer :: i,j,k,m
 integer :: is,ie,js,je,ks,ke,di,dj,dk

! moment parameters
 type(particle_gen) :: mmt
 real(sdsu_fps),allocatable :: wgt(:) !weight
 real(sdsu_fps) :: mth      ! moment order

 select case(trim(type_microphysics))
 case('GEN')
 case default
     if(masterproc) print*, 'MSG write_out_gen_moment: microphysics is',trim(cloud_microphysics)
     call stop_sdsu(myrank,'MSG write_out_gen_moment:termnate program')
     return
 end select 

!
! allocate parameters
! 
 if( .not. allocated(wgt) ) then
   allocate( wgt(1:nbin) )
 endif

!
! --------------- write out parameters in  netcdf file  ------------------------
!

!
! environmental parameters
!
 L2_radar_0d%press = a%press ! atmospheric pressure  [hPa]
 L2_radar_0d%t_air = a%t_air ! atmospheric temperature [K]
 L2_radar_0d%sh    = a%sh    ! specific humidity [%]
 L2_radar_0d%w     = a%w     ! vertical velocity [m/s]


!
! rainfall rate profile 
!
 L2_radar_0d%rain_rate = a%rain_rate ! rain rate [mm/hr]

!
! rain water content (rainning particles)
!
  L2_radar_0d%rwc = q%rain + q%snow + q%graupel + q%hail  ! [g/m3]

!
! liquid fraction (rainning particles)
!
 if( L2_radar_0d%rwc <= 0.) then
     L2_radar_0d%lfrac = 0.
 else
     L2_radar_0d%lfrac = q%rain / L2_radar_0d%rwc  !liquid fraction [-]
 endif


!
! cloud ice water content (non-rainning particles)
!
 L2_radar_0d%iwc =  q%ice ! [g/m3]

!
! cloud liquid water content (non-rainning particle)
!
 L2_radar_0d%lwc = q%cloud   ! cloud liquid water content [g/m3]


!
! 0th - 6th-order moments of raining particle (radius > 100um) used for gamma parameters.
!
 MOMENT_LOOP: do m = 0, 6
       mth = REAL(m)  !moment order 0~6th
       wgt(1:nbin) = 1.0e0   ! weight
       call gen_moment( mth, q%rain   , re%rain   , rho_gen%rain   ,mu%rain   , mmt%rain   )
       call gen_moment( mth, q%snow   , re%snow   , rho_gen%snow   ,mu%snow   , mmt%snow   )
       call gen_moment( mth, q%graupel, re%graupel, rho_gen%graupel,mu%graupel, mmt%graupel)
       call gen_moment( mth, q%hail   , re%hail   , rho_gen%hail   ,mu%hail   , mmt%hail   )
       L2_radar_0d%mmt(m) = mmt%rain + mmt%snow + mmt%graupel + mmt%hail   !moment [ m**mth / m3 ]
 enddo MOMENT_LOOP

!
! bulk density (raining particle)
!
 if( const_pi*L2_radar_0d%mmt(3) <= 0.) then
     L2_radar_0d%bulk_den = 0.e0
 else
     L2_radar_0d%bulk_den = L2_radar_0d%rwc / (const_pi*L2_radar_0d%mmt(3)/6.) *1.e-6 ! bulk density [g/cm3]
 endif

 return
 end subroutine L2_radar_parameters_gen

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine L2_radar_parameters_sbm(q_sbm, n_sbm, a, L2_radar_0d )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Write out various momemnt from SBM inputs for studying particle-size distributions.  
! 
! History:
!  03/2012  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 type(particle_sbm),intent(in) :: q_sbm         ! particle mixing ratio [g/m3] 
 type(particle_sbm),intent(in) :: n_sbm(1:nbin) ! particle # concentration [1/m4]
 type ( atmos_parameter ) ,intent(in) :: a

 type(particle_retrieval),intent(out) :: L2_radar_0d  !Level 2 radar parameters

 integer :: ierr,irec
 integer :: i,j,k,m
 character(len=20) :: typ
 integer,parameter :: io = 111
 integer :: is,ie,js,je,ks,ke,di,dj,dk

! moment parameters
 type(particle_sbm) :: mmt
 real(sdsu_fps),allocatable :: wgt(:) !weight
 real(sdsu_fps) :: mth !moment order
 real(sdsu_fps) :: fmelt  !melt fraction
 real(sdsu_fps) :: mass

!
! check microphysics
!
 if( trim(cloud_microphysics) == 'HUCM_SBM43' .or. &
     trim(cloud_microphysics) == 'HUCM_SBM' ) then
 else
     if(masterproc) print*, 'MSG L2_radar_parameters_sbm: microphysics must be SBM'
     return
 endif

!
! allocate parameters
! 
 if( .not. allocated(wgt) ) then
   allocate( wgt(1:nbin) )
 endif

!
! --------------- write out parameters in  netcdf file  ------------------------
!

!
! environmental parameters
!
 L2_radar_0d%press = a%press ! atmospheric pressure  [hPa]
 L2_radar_0d%t_air = a%t_air ! atmospheric temperature [K]
 L2_radar_0d%sh    = a%sh    ! specific humidity [%]
 L2_radar_0d%w     = a%w     ! vertical velocity [m/s]

!
! rainfall rate profile 
!
 L2_radar_0d%rain_rate = a%rain_rate ! rain rate [mm/hr]

!
! rain water content (rainning particles)
!
  mth = 3.e0            ! moment order 

  wgt(1:nbin) = den_sbm(1:nbin)%liq*(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
  call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%liq, drad_sbm(1:nbin)%liq,&
                  q_sbm%liq, n_sbm(1:nbin)%liq, mmt%liq )

  wgt(1:nbin) = den_sbm(1:nbin)%ice_col*(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
  call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_col, drad_sbm(1:nbin)%ice_col,&
                  q_sbm%ice_col, n_sbm(1:nbin)%ice_col, mmt%ice_col )

  wgt(1:nbin) = den_sbm(1:nbin)%ice_pla*(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
  call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_pla, drad_sbm(1:nbin)%ice_pla,&
                   q_sbm%ice_pla, n_sbm(1:nbin)%ice_pla, mmt%ice_pla )

  wgt(1:nbin) = den_sbm(1:nbin)%ice_den*(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
  call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_den, drad_sbm(1:nbin)%ice_den,&
                  q_sbm%ice_den, n_sbm(1:nbin)%ice_den, mmt%ice_den )

  wgt(1:nbin) = den_sbm(1:nbin)%snow   *(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
  call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%snow, drad_sbm(1:nbin)%snow,&
                  q_sbm%snow, n_sbm(1:nbin)%snow, mmt%snow )

  wgt(1:nbin) = den_sbm(1:nbin)%graupel*(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
  call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%graupel, drad_sbm(1:nbin)%graupel,&
                  q_sbm%graupel, n_sbm(1:nbin)%graupel, mmt%graupel )

  wgt(1:nbin) = den_sbm(1:nbin)%hail   *(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
  call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%hail, drad_sbm(1:nbin)%hail,&
                  q_sbm%hail, n_sbm(1:nbin)%hail, mmt%hail )

  L2_radar_0d%rwc =  mmt%liq + mmt%ice_col + mmt%ice_pla + mmt%ice_den &
                       + mmt%snow    + mmt%graupel  + mmt%hail ! ice water content [g/m3]

!
! liquid fraction (rainning particles)
!
 if( L2_radar_0d%rwc <= 0.) then
     L2_radar_0d%lfrac = 0.
 else
     L2_radar_0d%lfrac = mmt%liq / L2_radar_0d%rwc  !liquid fraction [-]
 endif


!
! cloud ice water content (non-rainning particles)
!
 mth = 3.e0            ! moment order 
 wgt(1:nbin) = den_sbm(1:nbin)%ice_col*(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
 call sbm_moment('non_rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_col, drad_sbm(1:nbin)%ice_col,&
                 q_sbm%ice_col, n_sbm(1:nbin)%ice_col, mmt%ice_col )

 wgt(1:nbin) = den_sbm(1:nbin)%ice_pla*(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
 call sbm_moment('non_rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_pla, drad_sbm(1:nbin)%ice_pla,&
                 q_sbm%ice_pla, n_sbm(1:nbin)%ice_pla, mmt%ice_pla )

 wgt(1:nbin) = den_sbm(1:nbin)%ice_den*(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
 call sbm_moment('non_rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_den, drad_sbm(1:nbin)%ice_den,&
                  q_sbm%ice_den, n_sbm(1:nbin)%ice_den, mmt%ice_den )

 wgt(1:nbin) = den_sbm(1:nbin)%snow   *(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
 call sbm_moment('non_rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%snow, drad_sbm(1:nbin)%snow,&
                 q_sbm%snow, n_sbm(1:nbin)%snow, mmt%snow )

 wgt(1:nbin) = den_sbm(1:nbin)%graupel*(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
 call sbm_moment('non_rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%graupel, drad_sbm(1:nbin)%graupel,&
                  q_sbm%graupel, n_sbm(1:nbin)%graupel, mmt%graupel )

 wgt(1:nbin) = den_sbm(1:nbin)%hail   *(const_pi/6.)*1.e+6   ! weight (density) [g/m3]
 call sbm_moment('non_rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%hail, drad_sbm(1:nbin)%hail,&
                  q_sbm%hail, n_sbm(1:nbin)%hail, mmt%hail )

 L2_radar_0d%iwc =  mmt%ice_col + mmt%ice_pla + mmt%ice_den &
                      + mmt%snow    + mmt%graupel  + mmt%hail ! ice water content [g/m3]


!
! cloud liquid water content (non-rainning particle)
!
 mth = 3.e0            ! moment order 

 wgt(1:nbin) = den_sbm(1:nbin)%liq * (const_pi/6.)*1.e+6   ! weight (density) [g/m3]
 call sbm_moment('non_rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%liq, drad_sbm(1:nbin)%liq,&
                 q_sbm%liq, n_sbm(1:nbin)%liq, mmt%liq )

 L2_radar_0d%lwc = mmt%liq   ! cloud liquid water content [g/m3]


!
! 0th, 3rd, 6th-order moments of raining particle (radius > 100um) used for gamma parameters.
!
 MOMENT_LOOP: do m = 0, 6

       mth = REAL(m)  !moment order 0~6th
       wgt(1:nbin) = 1.0e0   ! weight

       call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%liq, drad_sbm(1:nbin)%liq,&
                        q_sbm%liq, n_sbm(1:nbin)%liq, mmt%liq ) 
       call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_col, drad_sbm(1:nbin)%ice_col,&
                        q_sbm%ice_col, n_sbm(1:nbin)%ice_col, mmt%ice_col )
       call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_pla, drad_sbm(1:nbin)%ice_pla,&
                        q_sbm%ice_pla, n_sbm(1:nbin)%ice_pla, mmt%ice_pla )
       call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_den, drad_sbm(1:nbin)%ice_den,&
                        q_sbm%ice_den, n_sbm(1:nbin)%ice_den, mmt%ice_den )
       call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%snow, drad_sbm(1:nbin)%snow,&
                        q_sbm%snow, n_sbm(1:nbin)%snow, mmt%snow )
       call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%graupel, drad_sbm(1:nbin)%graupel,&
                        q_sbm%graupel, n_sbm(1:nbin)%graupel, mmt%graupel )
       call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%hail, drad_sbm(1:nbin)%hail,&
                        q_sbm%hail, n_sbm(1:nbin)%hail, mmt%hail )
       L2_radar_0d%mmt(m) = mmt%liq + mmt%ice_col + mmt%ice_pla + mmt%ice_den &
                              + mmt%snow    + mmt%graupel  + mmt%hail ! 0th [1/m3], 3rd [m3/m3], 6th [m6/m3]
 enddo MOMENT_LOOP


    !
    !  derive gamma parametrs from moments value. 
    !
!    if( L2_radar(i,j,k)%rwc >= q_min_condensate .and. &
!        mmt_order(1) > 0. .and. mmt_order(2) > 0. .and. mmt_order(3) > 0. ) then
!       call get_params_bisect(mmt_order(1),mmt_order(2),mmt_order(3),L2_radar(i,j,k)%psd)
!    else
!       L2_radar(i,j,k)%psd%mu     = undefined ! mu        [-]
!       L2_radar(i,j,k)%psd%lambda = undefined ! slope     [1/m]
!       L2_radar(i,j,k)%psd%n0     = undefined ! intercept [1/m3]
!       L2_radar(i,j,k)%psd%nw     = undefined ! intercept [1/m3]
!       L2_radar(i,j,k)%psd%d0     = undefined ! slope     [mm]
!    endif


!
! bulk density (raining particle)
!
 mth = 3.e0            ! moment order 

 wgt(1:nbin) = den_sbm(1:nbin)%liq       ! weight (density)[g/cm3]
 call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%liq, drad_sbm(1:nbin)%liq,&
                 q_sbm%liq, n_sbm(1:nbin)%liq, mmt%liq )
 wgt(1:nbin) = den_sbm(1:nbin)%ice_col   ! weight (density)[g/cm3]
 call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_col, drad_sbm(1:nbin)%ice_col,&
                  q_sbm%ice_col, n_sbm(1:nbin)%ice_col, mmt%ice_col )
 wgt(1:nbin) = den_sbm(1:nbin)%ice_pla   ! weight (density)[g/cm3]
 call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_pla, drad_sbm(1:nbin)%ice_pla,&
                  q_sbm%ice_pla, n_sbm(1:nbin)%ice_pla, mmt%ice_pla )
 wgt(1:nbin) = den_sbm(1:nbin)%ice_den   ! weight (density)[g/cm3]
 call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_den, drad_sbm(1:nbin)%ice_den,&
                  q_sbm%ice_den, n_sbm(1:nbin)%ice_den, mmt%ice_den )
 wgt(1:nbin) = den_sbm(1:nbin)%snow      ! weight (density)[g/cm3]
 call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%snow, drad_sbm(1:nbin)%snow,&
                 q_sbm%snow, n_sbm(1:nbin)%snow, mmt%snow )
 wgt(1:nbin) = den_sbm(1:nbin)%graupel   ! weight (density)[g/cm3]
 call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%graupel, drad_sbm(1:nbin)%graupel,&
                  q_sbm%graupel, n_sbm(1:nbin)%graupel, mmt%graupel )
 wgt(1:nbin) = den_sbm(1:nbin)%hail      ! weight (density)[g/cm3]
 call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%hail, drad_sbm(1:nbin)%hail,&
                  q_sbm%hail, n_sbm(1:nbin)%hail, mmt%hail )

 mass  = mmt%liq + mmt%ice_col + mmt%ice_pla + mmt%ice_den &
       + mmt%snow    + mmt%graupel  + mmt%hail ! [g/cm3]

 if( mass <= 0.) then
     L2_radar_0d%bulk_den = 0.e0
 else
     L2_radar_0d%bulk_den = mass / L2_radar_0d%mmt(3)  ! bulk density [g/cm3]
 endif


 return
 end subroutine L2_radar_parameters_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  find_hid( sback, hid , hid_confidence) 
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Find out hydrometeor identification by choosing HID having maximum backscatter coefficients.
! 
! History:
!  12/2013  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 type( particle_hid ),intent(in) :: sback(:) !backscatter 
 real,intent(out) :: hid(:)  !hydrometero identification 
 real,intent(out) :: hid_confidence(:) ! confidence of level [%]
 real :: sback_all(1:10)     !backscattering for all sub HID species
 integer :: k, kmax          !looping index
 integer :: iloc, n_loc(1)   !location index

 hid(:) = undefined !initialize output
 hid_confidence(:) = undefined

!
! define vertical layers
!
 kmax = size( sback(:)%dz )

!
! find out major HID for each layer
!
 do k = 1, kmax

    sback_all(1)  = sback(k)%dz
    sback_all(2)  = sback(k)%rn
    sback_all(3)  = sback(k)%cr
    sback_all(4)  = sback(k)%ds
    sback_all(5)  = sback(k)%ws
    sback_all(6)  = sback(k)%vi
    sback_all(7)  = sback(k)%ldg
    sback_all(8)  = sback(k)%hdg
    sback_all(9)  = sback(k)%ha
    sback_all(10) = sback(k)%bd

    ! find largest
    if( sum( sback_all(1:10) ) <= 0.e0 ) then
       hid(k) = undefined  ! Un-identified
       hid_confidence(k) = undefined
    else
      ! find majority
      n_loc(1:1) = maxloc( sback_all(1:10) )
      iloc = n_loc(1)
      hid(k) = real( iloc )  !hydrometeor ideintification (0~10) 

      !
      ! HID confidence of backscatter fraction
      !
      hid_confidence(k) = ( sback_all(iloc) / sum(sback_all(1:10)) )  * 100.e0  ![%]

    endif

   if( hid_confidence(k) > 100.e0 .or. hid_confidence(k) < undefined ) then
     print*,'MSG find_hid: strange hid_confidence(k) value',hid_confidence(k)
     print*,'iloc=',iloc
     print*, sback_all(iloc),  sum(sback_all(1:10))
     stop
   endif

 enddo

 return
 end subroutine find_hid

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  find_hid0d( sback, hid , hid_confidence)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Find out hydrometeor identification by choosing HID having maximum backscatter coefficients.
! 
! History:
!  12/2013  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 type( particle_hid ),intent(in) :: sback !backscatter 
 real,intent(out) :: hid  !hydrometero identification 
 real,intent(out) :: hid_confidence ! confidence of level [%]
 real :: sback_all(1:10)     !backscattering for all sub HID species
 integer :: iloc, n_loc(1)   !location index

 hid = undefined !initialize output
 hid_confidence = undefined

    sback_all(1)  = sback%dz
    sback_all(2)  = sback%rn
    sback_all(3)  = sback%cr
    sback_all(4)  = sback%ds
    sback_all(5)  = sback%ws
    sback_all(6)  = sback%vi
    sback_all(7)  = sback%ldg
    sback_all(8)  = sback%hdg
    sback_all(9)  = sback%ha
    sback_all(10) = sback%bd

    ! find largest
    if( sum( sback_all(1:10) ) <= 0.e0 ) then
       hid = undefined  ! Un-identified
       hid_confidence = undefined
    else
      ! find majority
      n_loc(1:1) = maxloc( sback_all(1:10) )
      iloc = n_loc(1)
      hid = real( iloc )  !hydrometeor ideintification (0~10) 

      !
      ! HID confidence of backscatter fraction
      !
      hid_confidence= ( sback_all(iloc) / sum(sback_all(1:10)) )  * 100.e0  ![%]

    endif

   if( hid_confidence > 100.e0 .or. hid_confidence < undefined ) then
     print*,'MSG find_hid0d: strange hid_confidence value',hid_confidence
     print*,'iloc=',iloc
     print*, sback_all(iloc),  sum(sback_all(1:10))
     stop
   endif

 return
 end subroutine find_hid0d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine sample_psd_sbm( proc, spc, bulkrad, dbulkrad, wc, n0,fmelt_bin )
 implicit none
!----------------------------------------------------------------------------------------------------
!  Purpose 
!
! 06/2014  Toshi Matsui @ NASA GSFC ; initial
!----------------------------------------------------------------------------------------------------
 character(len=*) :: proc ! process
 character(len=*),intent(in),optional :: spc      ! index name of hydrometeor classes 
 real(sdsu_fps), intent(in),optional  :: bulkrad(1:nbin) ! bulk radius of hydrometero for each bins [cm] 
 real(sdsu_fps), intent(in),optional  :: dbulkrad(1:nbin)! binsize [cm]
 real(sdsu_fps), intent(in),optional  :: wc        ! water content of hydrometeor [g/m**3]
 real(sdsu_fps), intent(in),optional  :: n0(1:nbin)! number concentration of hydrometeors [1/m4]
 real(sdsu_fps), intent(in),optional  :: fmelt_bin(1:nbin) !mass melting fraction of hydrometoer for each bins [-] 

!local parameter
 real(sdsu_fps),allocatable,dimension(:,:),save :: sum_psd, mean_psd  !sum and mean of PDF [#/m3]
 real(sdsu_fps),allocatable,dimension(:),save :: cnt_sample  !sampling count [-]
 integer :: ierr
 integer :: s
 integer :: n
 integer,parameter :: io = 699

 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

  if(masterproc) print*,'MSG sample_psd_sbm: initialize parameters'

  if( .not. allocated(sum_psd) ) then
   allocate( sum_psd (1:nbin,1:3), &
             mean_psd(1:nbin,1:3), &
             cnt_sample(1:3) , &
             stat=ierr )
   if (ierr /= 0) call stop_sdsu(myrank,'MSG sample_psd_sbm: allocation error -> Terminate program.')
   sum_psd    = 0.e0 
   mean_psd   = 0.e0
   cnt_sample = 0.e0
  endif

 case('add') !-------------------------------------------------------------------------------------

  if(wc .lt. q_min_micro) then
     return
  endif

 if( trim(spc) == 'hail' ) then

  sum_psd(1:nbin,1) = sum_psd(1:nbin,1) + n0(1:nbin) * dbulkrad(1:nbin)*1.e-2
  cnt_sample(1) = cnt_sample(1) + 1.e0
 
  sum_psd(1:nbin,2) = sum_psd(1:nbin,2) + fmelt_bin(1:nbin)
  cnt_sample(2) = cnt_sample(2) + 1.e0

 elseif( trim(spc) == 'liq' ) then  

  sum_psd(1:nbin,3) = sum_psd(1:nbin,3) + n0(1:nbin) * dbulkrad(1:nbin)*1.e-2
  cnt_sample(3) = cnt_sample(3) + 1.e0

 endif

 case('finish') !---------------------------------------------------------------------------------------

 if( .not. allocated(sum_psd) ) then
    print*,'MSG sample_psd_sbm(finish): sum_psd is not allocated '
    return
 endif

#if MPI > 0 
  call mpi_sdsu_reduce('SUM', sum_psd    )
  call mpi_sdsu_reduce('SUM', cnt_sample )

  MPI_MASTERPROC: if(masterproc) then
#endif

 do s = 1, 3
   mean_psd(1:nbin,s) = sum_psd(1:nbin,s) / cnt_sample(s)  !derive mean value
 enddo

!
! output file
!
 sdsu_io_file =trim(sdsu_dir_output)//'RN_PSD.txt'

 print*,'MSG sample_psd_sbm: dump mean PSD. See --> ',trim(sdsu_io_file)

 open( io,file=trim(sdsu_io_file) )

 do n = 1, nbin
    write(io,100) n, mean_psd(n,1:3) 
 enddo
 close(io)

 100 format(I5, 5x, E12.5, 5x, f8.3, 5x, E12.5)



#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 end select

 end subroutine sample_psd_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine effective_permittivity(freq, fmelt, dens, atemp, crefd, &
                                    dcreal, dcimag)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Calculate effective dielectric constant of frozen particle using different methods/assumptions. 
! 
! History:
! 10/2016  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(in):: freq   ! frequency [GHz]
 real(sdsu_fps),intent(in):: fmelt  ! melt fraction (0~1) [-]  0.0-solid, 1.0-liquid, 
 real(sdsu_fps),intent(in):: dens   ! particle effective density [kg/m3]
 real(sdsu_fps),intent(in):: atemp  ! air temperature [K]
 complex*16,intent(out) ::  crefd   ! effective permittivity 
 real,intent(inout),optional :: dcreal
 real,intent(inout),optional :: dcimag  !real / imageneary component of dielectric function
 real(sdsu_fpd) :: creal,rimag

 real(sdsu_fps) ::&
   faa   ,& ! fraction of air in ice 
   eicere,& ! real component of ice complex refractive index
   eiceim,& ! imaginary component of ice complex refractive index
   eimag ,& ! complex refractive index
   ereal ,& !
   fmelt_new   ! updated melt fraction [-]

 real(sdsu_fpd) :: &
     xd      ,& ! size parameter [-] in double precision 
     freqd   ,& ! frequency [GHz]
     freqhzd ,& ! frequency in double precision 
     tempd   ,& ! temperaure in doubpl precision
     sald    ,& ! salinity of water
     eicered ,& ! real component of ice complex refractive index
     ewatred ,& ! real component of water complex refractive index
     eiceimd ,& ! imaginary component of ice complex refractive index
     ewatimd    ! imaginary component of water complex refractive index

 complex *8 :: eice,ewat,eair !ice, water, and air permittivity
 complex *8 :: ei     !effective ice-air permittivity
 complex *8 :: emelt,emelt_wi,emelt_iw, emelt_em !effective water-ice-air permittivity
 complex*16 :: ewatd
 complex*16 :: eid, emeltd
 complex*16 :: crefd_liq, crefd_ice
 real(sdsu_fps),parameter :: densice = 0.917e+3  !solid ice density [kg/m**3]

!
! check input
!
 if(fmelt < 0. .or. fmelt > 1. ) then
   print*,'MSG effective_permittivity: fmelt out of range',fmelt
   call stop_sdsu(myrank,'MSG effective_permittivity: terminating program.') 
 endif
 
!
!
!
 tempd = dble(atemp)

 ICE: if(fmelt < 1.0) then ! if including solid particle, compute ice refractive index 

!
! complex refractive index of solid ice 
!
   freqd=dble(freq)  !frequency [GHz] in double precision
   call iceoptic(freqd,tempd,eicered,eiceimd)  !get complex refractive index
   eicere=REAL(eicered)      !real      component of complex refractive index
   eiceim=REAL(eiceimd)      !imaginary component of complex refractive index
   eice=cmplx(eicere,eiceim) !complex refractive index of ice
   eair=cmplx(1.0006,0.0)    !complex refractive index of air

!
! Calculate effective dielectric constant of fluffy frozen particle using different methods/assumptions. 
!
   faa =1.e0 - (dens/densice)  ! volume fraction of inclusions (air)

   ice_select: select case(ice_refraction_func)
   case(1)  ! Maxwell-Garnet air include ice 

      call mg_ellips(1.e0-faa,eair,eice,ei)
      eid=dcmplx(ei)
      crefd_ice=cdsqrt(eid)

   case(2)  ! Maxwell-Garnet ice include air

      call mg_ellips(faa,eice,eair,ei)
      eid=dcmplx(ei)
      crefd_ice=cdsqrt(eid)

   case(3)  ! Effective Medium homogeneous mixing of ice and air

      call em_ellips( 1.e0-faa, eair, eice, ei)
      eid=dcmplx(ei)
      crefd_ice=cdsqrt(eid)

   case(4)  ! Debye theory

      call debye_mixture( dens*1.d-3, eicered , eiceimd , ei )
      eid=dcmplx(ei)
      crefd_ice=cdsqrt(eid)

   case default ; call stop_sdsu(myrank,'MSG effective_permittivity: There is no such ice_refraction_func')
   end select ice_select

 endif ICE


 LIQ: if(fmelt > 0.0) then !liquid or melting particles
!
! complex refractive index of liquid water
!
   freqhzd=dble(freq*1.e+9)      !frequency
   sald=0.d0                     !salinity
   call watoptic(freqhzd,tempd,sald,ewatred,ewatimd)  !get complex refractive index 
   ewatd=dcmplx(ewatred,ewatimd) !double precision complex*16
   crefd_liq=cdsqrt(ewatd)       !square root of complex refractive index
   ewat=cmplx(ewatred,ewatimd)   !complex*8

 endif LIQ

!
! choose either liquid, ice, or melting
!
 if(fmelt == 0.) then    ! solid ice

   crefd = crefd_ice !output
   !
   ! optional output 
   !
   if(present(dcreal) ) dcreal = real(ei) ; if(present(dcimag)) dcimag = aimag(ei)

 elseif(fmelt== 1.0 ) then ! liquid condensates

   crefd = crefd_liq !output
   !
   ! optional output 
   !
   if(present(dcreal) ) dcreal=real(ewatred) ; if(present(dcimag)) dcimag=real(ewatimd)

 else                       ! mixed-phase condensates 

     melt_select: select case(melt_opt)
     case(0)   ! no melting scheme (just use ice dielectric const.)
         crefd = crefd_ice
         !
         ! optional output 
         !
         if(present(dcreal) ) dcreal = real(ei) ; if(present(dcimag)) dcimag = aimag(ei)
     case(1)   ! water include ice
        call mg_ellips( 1.-fmelt, ewat, ei, emelt_wi )
        emeltd=dcmplx(emelt_wi)
        crefd = cdsqrt(emeltd)
        !
        ! optional output 
        !
        if(present(dcreal) )  dcreal=real(emelt_wi) ; if(present(dcimag)) dcimag=aimag(emelt_wi)

     case(2)   ! ice include water     
         call mg_ellips( fmelt, ei, ewat, emelt_iw)
         emeltd=dcmplx(emelt_iw)
         crefd = cdsqrt(emeltd)
        !
        ! optional output 
        !
        if(present(dcreal) )  dcreal=dble(emelt_iw) ; if(present(dcimag)) dcimag=aimag(emelt_iw)

     case(3)   !averaging above two solutions
         call mg_ellips( 1.-fmelt, ewat, ei, emelt_wi )
         call mg_ellips( fmelt, ei, ewat, emelt_iw)
         eimag = 10.**( log10(aimag(emelt_wi))*(fmelt) + log10(aimag(emelt_iw))*(1.-fmelt)    )
         ereal = 10.**( log10(dble(emelt_wi))*fmelt + log10(dble(emelt_iw))*(1.-fmelt)    )
         emelt = cmplx(ereal,eimag)
         emeltd=dcmplx(emelt)
         crefd = cdsqrt(emeltd)
        !
        ! optional output 
        !
        if(present(dcreal) )  dcreal=real(emelt_wi) ; if(present(dcimag)) dcimag=aimag(emelt_wi)

     case(4)   ! effective medium homogeneous mixing of ice and water
         call em_ellips( 1-fmelt, ewat, ei, emelt_em)
         emeltd=dcmplx(emelt_em)  !ewatd
         crefd = cdsqrt(emeltd)
        !
        ! optional output 
        !
        if(present(dcreal) ) dcreal=real(emelt_em) ; if(present(dcimag)) dcimag=aimag(emelt_em)

     case default ; call stop_sdsu(myrank,'MSG ss_rams_micro: There is no such melt_opt')
     end select melt_select

 endif


 return
 end subroutine effective_permittivity

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_n0( mu, q, dens, re, n0 )
 implicit none
 real(sdsu_fps),intent(in) :: mu  !dispersion parameter [-]
 real(sdsu_fps),intent(in) :: q   !mass conc. [g/m3]
 real(sdsu_fps),intent(in) :: dens!bulk hydrometeor density [kg/m3]
 real(sdsu_fps),intent(in) :: re  !effective radius [micron]
 real(sdsu_fps),intent(out) :: n0 !intercept [1/m**(4+mu)]

 real(sdsu_fps) :: gfac1, gfac2, gamfac
 real(sdsu_fps) :: lam  ! slope paraameter [1/m]


 if(q <= 0. .or. re <= 0.) then
   n0 = 0. ; return
 endif

!
! derive lambda and intercept of expoential DSD for unit water content
!
  call gamma_reff(mu+4.0,gfac1)
  call gamma_reff(mu+3.0,gfac2)
  gamfac = gfac1/gfac2

! derive lambda
  lam = 1. / (2.*re*1.e-6) * gamfac  !slope [1/m]

! derive intercept
!        [g/m3]                    [m3/g]          [1/m]
  n0  =     q * 6. / const_pi / (dens*1000.) * (lam**(4.+mu))  / gfac1  ! [1/m**(4+mu)]

 end subroutine get_n0

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_opt_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
