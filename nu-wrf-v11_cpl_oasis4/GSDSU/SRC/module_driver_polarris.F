#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_driver_polarris
 use module_simulator
 use module_opt_micro
 use module_tmatrix
 use module_mueller_matrix
 use module_gmp_4ice
 use module_p3_param
 use module_mpi

 implicit none

!---------------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
!
! Comments:  
!   This module drives Muleller and T-matrix modules to calculate polarimetric radar
!   observables and retrievals.  
!
! History: 
!  02/2017  Toshi Matsui@NASA GSFC ; Created LUT subroutoines
!  01/2017  Toshi Matsui@NASA GSFC ; debbuged radial velocity.
!  09/2016  Toshi Matsui@NASA GSFC ; Initial development applied from module_opt_microwave.
!           
! References: 
!  Not yet.
!
!---------------------------------------------------------------------------------------------------
 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 private   ! encapsulate all variables and subourtines (non accessible)


!
! public (accessible) subroutines
!
 public :: read_polarris_axisorient, &   ! reading axis ratio and orienting angle table
           prep_makelut_polarris, &    ! preparing LUT call
           makelut_polarris, &         ! create microwave single-scattering properties LUTs
           opt_polarris                ! compute microwave single-scattering properties

 public :: opt_polarris_parameter  !optical properties.

!
! 3D observables parameters 
!
 type opt_polarris_parameter 
  real(sdsu_fps) :: zhh     ! reflectivity [dBZ] 
  real(sdsu_fps) :: zdr     ! differential reflectivity [dB]
  real(sdsu_fps) :: kdp     ! specific diff phase [deg/km] 
  real(sdsu_fps) :: rhohv   ! copolar (hv) correlation coefficient [-]
  real(sdsu_fps) :: vdop    ! Doppler velocity [cm/s]
  real(sdsu_fps) :: vrad    ! radial velocity [cm/s]
  real(sdsu_fps) :: elev    ! radar elevation angle [deg]
  real(sdsu_fps) :: vang    ! radar viewing angle [deg] (starting from north ward 0-359.999).
 end type  opt_polarris_parameter
 type (opt_polarris_parameter),allocatable :: atmos_optpol(:,:,:,:)


!
! 3D retrieval parameter
!
 type geo_polarris_parameter
    real(sdsu_fps) :: t_air     ! layer air temperature [degK]
    real(sdsu_fps) :: rh        ! relative humidity [%]
    real(sdsu_fps) :: u         ! u wind (eastward) [m/sec]
    real(sdsu_fps) :: v         ! v wind (northward) [m/sec]
    real(sdsu_fps) :: w         ! w wind (upward)  [m/sec]
    real(sdsu_fps) :: rain_rate ! rain rate profile [mm/hr]
    real(sdsu_fps) :: qc        ! cloud mixing ratio  [g/m3]
    real(sdsu_fps) :: qr        ! rain mixing ratio  [g/m3]
    real(sdsu_fps) :: qi        ! ice mixing ratio  [g/m3]
    real(sdsu_fps) :: qs        ! snow mixing ratio  [g/m3]
    real(sdsu_fps) :: qg        ! graupel mixing ratio  [g/m3]
    real(sdsu_fps) :: qh        ! hail mixing ratio  [g/m3]
    real(sdsu_fps) :: rer       ! rain effective radius [um]
    real(sdsu_fps) :: res       ! snow effective radius [um]
    real(sdsu_fps) :: reh       ! hail effective radius [um]
 end type  geo_polarris_parameter
 type (geo_polarris_parameter),allocatable :: atmos_geopol(:,:,:)


!
! giant POLARRIS Look-up tables
!
 real(sdsu_fpd),allocatable :: bmueller_array  (:,:,:,:,:,:) !(4,4,mxpts_re_pol,max_spc,mxpts_temp_pol,mxelev_polarris)
 real(sdsu_fps),allocatable :: fscat_real_array(:,:,:,:,:,:) !(2,2,mxpts_re_pol,max_spc,mxpts_temp_pol,mxelev_polarris)
 real(sdsu_fps),allocatable :: fscat_imag_array(:,:,:,:,:,:) !(2,2,mxpts_re_pol,max_spc,mxpts_temp_pol,mxelev_polarris)
 real(sdsu_fpd),allocatable :: vdop_array      (    :,:,:,:) !(    mxpts_re_pol,max_spc,mxpts_temp_pol,mxelev_polarris)


!
! giant POLARRIS Look-up tables for p3
!
! bmueller_array  (4,4,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
! fscat_real_array(2,2,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
! fscat_imag_array(2,2,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
! vdop_array      (    n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &

 real(sdsu_fpd),allocatable :: bmueller_array_p3  (:,:,:,:,:,:,:) !
 real(sdsu_fps),allocatable :: fscat_real_array_p3(:,:,:,:,:,:,:) !
 real(sdsu_fps),allocatable :: fscat_imag_array_p3(:,:,:,:,:,:,:) !
 real(sdsu_fpd),allocatable :: vdop_array_p3      (    :,:,:,:,:) !






! minimum condensate mixing ratio for microwave optical properties
  real(sdsu_fps),parameter :: q_min_micro = q_min_condensate  ! minimum detectable condensate amount from CRM [g/cm3] 

!  integer, public :: ibs,ibe,jbs,jbe !extra buffer domain bounds

! minimum extinction for stability
  real(sdsu_fps),parameter :: min_bext = 1.e-10   ! 1e-6 default

!  real(sdsu_fps),parameter :: const_ghz_mm = 2.997925e2

!
! for general bulk microphysics
!
 type ( particle_gen_i ) :: & 
     num_axis_gen     , & ! number of axis ratio 
     ityp_axis_gen(max_ens) , & ! type of axis ratio (-999 is undefined)
     ityp_dist_gen    , & ! particle orienting angle distribution type
     num_stdev_gen        ! number of orienting angle standard deviations

 type ( particle_gen ) :: & 
     axis_ratio_gen(max_ens)   , & ! particle axis ratio [-]
     mean_orient_gen             , & ! mean orienting angle [deg]
     stdev_orient_gen(max_ens)       ! standard deviation / theta_m of orienting angle [deg]

!
! for rams 
!
 type ( particle_rams_i ) :: &
     num_axis_rams     , & ! number of axis ratio 
     ityp_axis_rams(max_ens) , & ! type of axis ratio (-999 is undefined)
     ityp_dist_rams    , & ! particle orienting angle distribution type
     num_stdev_rams        ! number of orienting angle standard deviations

 type ( particle_rams ) :: & 
     axis_ratio_rams(max_ens)   , & ! particle axis ratio [-]
     mean_orient_rams             , & ! mean orienting angle [deg]
     stdev_orient_rams(max_ens)       ! standard deviation / theta_m of orienting angle [deg]

!
! for sbm
!
 type ( particle_sbm_i ) :: &
     num_axis_sbm     , & ! number of axis ratio 
     ityp_axis_sbm(max_ens) , & ! type of axis ratio (-999 is undefined)
     ityp_dist_sbm    , & ! particle orienting angle distribution type
     num_stdev_sbm        ! number of orienting angle standard deviations

 type ( particle_sbm ) :: & 
     axis_ratio_sbm(max_ens)   , & ! particle axis ratio [-]
     mean_orient_sbm             , & ! mean orienting angle [deg]
     stdev_orient_sbm(max_ens)       ! standard deviation / theta_m of orienting angle [deg]

!
! tag of spectrum axis ratio 
!
 integer,parameter :: &
    tag_chuang_beard_1990      = 201  , & ! (equilibrium assumption. more oblate. ref zdr=2.42)
    tag_Pruppacher_Beard_1970  = 202  , & ! (equilibrium assumption. more oblate. ref zdr=2.36)
    tag_jameson_1983           = 203  , & ! (oscilating ave . less oblate.        ref zdr=0.76) 
    tag_green_1975             = 204  , & ! (equilibrium assumption. parabolic.   ref zdr=2.33)
    tag_pruppacher_pitter_1971 = 205  , & ! (equilibrium assumption. more oblate. ref zdr=2.31)
    tag_jones_i_1959           = 206  , & ! (equilibrium assumption. more oblate. ref zdr=2.10)
    tag_jones_ii_1959          = 207  , & ! (equilibrium assumption. more oblate. ref zdr=1.84)
    tag_goddard_cherry_1984    = 208  , & ! (equilibrium assumption. more oblate. ref zdr=2.41)
    tag_brandes_2002           = 209  , & ! (oscilating ave. less oblate.         ref zdr=0.76)
    tag_ryzhkov_gh_2011        = 210  , & ! (dry graupel/hail axis ratio )
    tag_hail_knight_1986       = 211  , & ! hail model of empirical fit from knight 1986. 
    tag_snow_mc3e_hvps3        = 212      ! snow aggregate axis ratio relationship from MC3E citation HVPS-3 measurment


 real(sdsu_fpd),parameter :: dpi=3.14159265358979323846d0    ! pai 
 real(sdsu_fpd),parameter :: degrad = dpi/180.d0             ! conversion from degree to radiance
 real(sdsu_fpd),parameter :: raddeg = 180.d0/dpi             ! conversion from radiance to degree 

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_polarris_axisorient 
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Read & initialize constant parameters of axis ratio and orienting angle distributions for
!  the polarimetric radar simulator (T-mueller matrix).
! 
! History: 
!  06/2017  Toshi Matsui@NASA GSFC ; initial.
!           
! References:
!---------------------------------------------------------------------------------------------------

 integer,parameter :: io = 101

 namelist /gen_switch/ num_axis_gen, ityp_axis_gen, axis_ratio_gen, ityp_dist_gen, mean_orient_gen, &
                       num_stdev_gen, stdev_orient_gen

 namelist /rams_switch/ num_axis_rams, ityp_axis_rams, axis_ratio_rams, ityp_dist_rams, mean_orient_rams, &
                       num_stdev_rams, stdev_orient_rams

 namelist /sbm_switch/ num_axis_sbm, ityp_axis_sbm, axis_ratio_sbm, ityp_dist_sbm, mean_orient_sbm, &
                       num_stdev_sbm, stdev_orient_sbm

 if(masterproc) print*,'MSG read_polarris_axisorient: reading axis ratio and orienting angle distributions'

!
! file that specifies axis ratio and orienting angle distributions.
!

 sdsu_io_file = trim(sdsu_dir_data)//'polarris_axisorient_'//trim(polarris_axisorient_name)//'.txt'
 if(masterproc) print*,sdsu_io_file

 call open_oldfile(io, trim(sdsu_io_file) )

 mic_select: select case(trim(type_microphysics))
 case('GEN') ! general bulk microphsics 
   read(io,gen_switch )
 case('RAMS')           ! RAMS type
   read(io,rams_switch)
 case('SBM')   ! HUCM SBM type
   read(io,sbm_switch )
 case default
   call stop_sdsu(myrank,'MSG read_polarris_axisorient: There is no such cloud_microphysics')
 end select mic_select

 close(io)


 return
 end subroutine read_polarris_axisorient 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_polarris
 implicit none 
!---------------------------------------------------------------------------------------------------
! Comments:  
!  This module drives the polarimetric radar simulator (T-Mueller matrix).
!  1. CRM streching vertical layers are intepolrated into constant height levels.
!  2. Calculate radar elevation angle for each grid.
!  3. Calculate polarimetric observables for each microphysics scheme.
!
! History: 
!  09/2016  Toshi Matsui@NASA GSFC ; initial.
!           
! References:
!---------------------------------------------------------------------------------------------------
 integer :: ierr !error index
 integer :: i,j,k,ni,m !looping indice
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
 real(sdsu_fps) :: wave_m ! wavelength [m]
 character(len=100):: char_freq, fname_out
 integer :: max_spc
 integer :: mxpts_re_pol_temp !temporal index

 type (atmos_parameter),allocatable :: atmos_refine1d(:)

 type (opt_polarris_parameter),allocatable :: atmos_optpol1d(:)  !1D column of optical properties
 
 type ( particle_gen_r8 ) :: & ! particle_gen is defined in module_simulater
     bmueller_gen(4,4),  &  ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_gen                  ! Doppler velocity [m/s]
 type ( particle_gen_cmplx ) :: & ! particle_gen is defined in module_simulater
     fscat_gen(2,2)               ! 2x2 forward scattering matrix 

 type ( particle_rams_r8 ) :: & !particle_rams is defined in module_simulater
     bmueller_rams(4,4),  &  ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_rams                  ! Doppler velocity [m/s]
 type ( particle_rams_cmplx ) :: & !particle_rams is defined in module_simulater
     fscat_rams(2,2)               ! 2x2 forward scattering matrix 

 type ( particle_sbm_r8 ) :: & !particle_sbm is defined in module_simulater
     bmueller_sbm(4,4),  &  ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_sbm                  ! Doppler velocity [m/s]
 type ( particle_sbm_cmplx ) :: & !particle_sbm is defined in module_simulater
     fscat_sbm(2,2)               ! 2x2 forward scattering matrix 

 real(sdsu_fpd) :: v_inst(3), v_inst_top(3), v_grid(3)
 real(sdsu_fps) :: elev(max_ens), dist !elevation angle [deg], distance [km]
 integer :: mxlyr_refine
 real(sdsu_fps) :: elev_inst !instrumental elevation [km] 
 real(sdsu_fps),allocatable :: stag_t_air(:)


 type ( particle_gen ), allocatable, dimension(:) :: & !1D parameters 
     q_gen_refine1d,  &  ! particle mixing ratio [g/m3]
     re_gen_refine1d     ! particle effective radius [micron]

 type ( particle_rams ), allocatable, dimension(:) :: & !1D parameters 
     q_rams_refine1d,  & ! particle mixing ratio [g/m3]
     re_rams_refine1d    ! particle effective radius [micron]

 type ( particle_sbm ), allocatable, dimension(:,:) :: & 
     n_sbm_refine1d      ! particle # concentration [1/m4]

 real(sdsu_fps), allocatable, dimension(:,:) :: &
    fmelt_sbm_snow_refine1d    ,& ! mass melting fraction of snow [-]
    fmelt_sbm_graupel_refine1d ,& ! mass melting fraction of graupel [-]
    fmelt_sbm_hail_refine1d    ,& ! mass melting fraction of hail [-]
    frime_sbm_snow_refine1d       ! mass riming fraction of snow  [-]

 type ( particle_sbm ), allocatable, dimension(:) :: & 
     q_sbm_refine1d    ,& ! particle mixing ratio [g/m3]
     re_sbm_refine1d      ! particle effective radius [micron]

 real(sdsu_fps), allocatable, dimension(:) :: &
   rho_graupel_refine1d ,& ! vertically refined graupel density [kg/m3]
   rho_hail_refine1d       ! vertically refined hail density [kg/m3]

 real(sdsu_fps), allocatable, dimension(:) :: &
   p3_rhorime_refine1d ,&  ! vertically refined rimed ice density [kg/m3]
   p3_frime_refine1d   ,&  ! vertically refined rime mass fraction  [-]
   p3_mu_cloud_refine1d,&  ! cloud mu [-]
   p3_num_ice_refine1d     ! ice number conc [#/m3]



 real(sdsu_fps) :: vang_surf  !surface-level viewing angle [deg]

 real(sdsu_fps) :: liq_tot, liq_cloud, sbm_cloud_frac  !fraction of cloud component within SBM liquid category

 integer :: nens, max_ensemble
 character(len=1) :: nchar


!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1

!
! Program Start
!
 if(masterproc) print*,'MSG opt_polarris: Compute polarimetric radar observables '


!
! Determine original CRM level or instrumental level (--> mxlyr_refine )
!
  if( use_crm_level_polarris ) then
     mxlyr_refine = mxlyr             !CRM level (eta streching)
  else
     mxlyr_refine = mxlyr_polarris    !constant height level

     !instrumental stagger height level [km]
     if( allocated(hgt_stag) ) deallocate( hgt_stag )
     allocate( hgt_stag(0:mxlyr_refine) )
     hgt_stag(0:mxlyr_refine) =  hgt_stag_polarris(0:mxlyr_refine)
  endif

!
! allocate 1D/3D array 
!
  if( allocated( stag_t_air ) ) deallocate( stag_t_air ) 
  allocate( stag_t_air(0:mxlyr_refine), stat = ierr ) 

  if( allocated(atmos_refine1d)) deallocate( atmos_refine1d )
  allocate( atmos_refine1d(1:mxlyr_refine), stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_polarris: allocation error atmos_refine1d. ')

  if( allocated(atmos_optpol)) deallocate( atmos_optpol )
  allocate( atmos_optpol (is:ie,js:je,1:mxlyr_refine,mxinst_polarris), stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_polarris: allocation error atmos_optpol. ')


  if( allocated(atmos_optpol1d) ) deallocate( atmos_optpol1d )
  allocate( atmos_optpol1d(1:mxlyr_refine), stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_polarris: allocation error atmos_optpol1d. ')

  select case(trim(type_microphysics))
  case('GEN')
    if( allocated(q_gen_refine1d)) deallocate( q_gen_refine1d )
    if( allocated(re_gen_refine1d)) deallocate( re_gen_refine1d )
    allocate( q_gen_refine1d (1:mxlyr_refine) , stat=ierr ) 
    allocate( re_gen_refine1d(1:mxlyr_refine) , stat=ierr )

    if( trim(cloud_microphysics) == 'NSSL2MMT' ) then
       if( .not. allocated(rho_graupel_refine1d) ) allocate( rho_graupel_refine1d(1:mxlyr_refine) , stat=ierr )
       if( .not. allocated(rho_hail_refine1d)    ) allocate( rho_hail_refine1d   (1:mxlyr_refine) , stat=ierr )
    endif 

    if( trim(cloud_microphysics) == 'P3' ) then
       if( .not. allocated(p3_rhorime_refine1d)     ) allocate( p3_rhorime_refine1d  (1:mxlyr_refine) , stat=ierr )
       if( .not. allocated(p3_frime_refine1d)       ) allocate( p3_frime_refine1d    (1:mxlyr_refine) , stat=ierr )
       if( .not. allocated(p3_mu_cloud_refine1d)    ) allocate( p3_mu_cloud_refine1d (1:mxlyr_refine) , stat=ierr )
       if( .not. allocated(p3_num_ice_refine1d )    ) allocate( p3_num_ice_refine1d  (1:mxlyr_refine) , stat=ierr )
    endif


  case('RAMS')
    if( allocated(q_rams_refine1d)) deallocate( q_rams_refine1d )
    if( allocated(re_rams_refine1d)) deallocate( re_rams_refine1d )
    allocate( q_rams_refine1d (1:mxlyr_refine) , stat=ierr )
    allocate( re_rams_refine1d(1:mxlyr_refine) , stat=ierr )
  case('SBM')
    if( allocated(q_sbm_refine1d)) deallocate( q_sbm_refine1d )
    if( allocated(re_sbm_refine1d)) deallocate( re_sbm_refine1d )
    if( allocated(fmelt_sbm_snow_refine1d)) deallocate( fmelt_sbm_snow_refine1d )
    if( allocated(fmelt_sbm_graupel_refine1d)) deallocate( fmelt_sbm_graupel_refine1d )
    if( allocated(fmelt_sbm_hail_refine1d)) deallocate( fmelt_sbm_hail_refine1d )
    if( allocated(frime_sbm_snow_refine1d)) deallocate( frime_sbm_snow_refine1d )

    allocate( q_sbm_refine1d (1:mxlyr_refine)  , stat=ierr )
    allocate( re_sbm_refine1d(1:mxlyr_refine)  , stat=ierr )
    allocate( n_sbm_refine1d            (1:mxlyr_refine,1:nbin) , stat=ierr ) 
    allocate( fmelt_sbm_snow_refine1d   (1:mxlyr_refine,1:nbin) , stat=ierr )
    allocate( fmelt_sbm_graupel_refine1d(1:mxlyr_refine,1:nbin) , stat=ierr )
    allocate( fmelt_sbm_hail_refine1d   (1:mxlyr_refine,1:nbin) , stat=ierr )
    allocate( frime_sbm_snow_refine1d   (1:mxlyr_refine,1:nbin) , stat=ierr )
  case default
  end select


  if(polarris_crmout) then
    if( allocated(atmos_geopol)) deallocate( atmos_geopol )
    allocate( atmos_geopol (is:ie,js:je,1:mxlyr_refine), stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_polarris: allocation error atmos_geopol. ')
  endif



!
! dolan 2013 ensemble assumptions of axis ratio and orienting angle
!
 if( polarris_axisorient_typ == 4 ) then
   max_ensemble = 4
 else
   max_ensemble = 1
 endif

 if(masterproc) print*, 'MSG opt_polarris: polarris_axisorient_typ=',polarris_axisorient_typ


!
! ensemble loop
!
 ensemble_loop: do nens = 1, max_ensemble

  if( polarris_axisorient_typ == 4 ) then
      write(nchar,"(I1.1)") nens ! 1 digit decimal 
      polarris_axisorient_name = 'dolan2013_pt'//nchar
      if(masterproc) print*, nens, max_ensemble, trim(polarris_axisorient_name)

  endif

!
! initialize POLARIS parameters
!
 call read_polarris_axisorient




!
! Loop over the channel frequencies
!
 INST_LOOP: do ni = 1, mxinst_polarris 

!
! Sensor frequency and wavelengths
!
  freqcy = freq_polarris(ni) ! passive microwave
  if(verbose_SDSU .and. masterproc)  &
  print*, ' Frequency =', freqcy,' GHz -> Wavelength = ',2.997925e2/freqcy ,'mm'
     
  wave_m = 2.997925e2/freqcy * 0.001 !wavelegnth [m]


!
! allocate LUT parameters
!
  if( lut_polarris ) then
    char_freq = nch_polarris(ni)

    mic_select0: select case(trim(type_microphysics))
    case('GEN')
        max_spc = 6
        mxpts_re_pol_temp = mxpts_re_pol
        if( trim(cloud_microphysics) == 'NSSL2MMT' ) max_spc = 17
        !
        ! P3 scheme only use cloud and rain for  general LUT
        ! 
        if ( trim(cloud_microphysics) == 'P3' ) max_spc = 2
        
    case('RAMS') 
        max_spc = 16
        mxpts_re_pol_temp = mxpts_re_pol
    case('SBM')
        max_spc = 16 
        mxpts_re_pol_temp = nbin
    case default 
        stop
    end select mic_select0
    fname_out = 'POLARRIS_'//trim(cloud_microphysics)//'-'//trim(char_freq)//'_'//trim(polarris_axisorient_name)//'.nc'

    !
    ! allocate giant LUT array and initialize parameters
    !
    if( .not. allocated(bmueller_array) ) then
      allocate( bmueller_array  (4,4,mxpts_re_pol_temp,max_spc,mxpts_temp_pol,mxelev_polarris), &
                fscat_real_array(2,2,mxpts_re_pol_temp,max_spc,mxpts_temp_pol,mxelev_polarris), &
                fscat_imag_array(2,2,mxpts_re_pol_temp,max_spc,mxpts_temp_pol,mxelev_polarris), &
                vdop_array      (    mxpts_re_pol_temp,max_spc,mxpts_temp_pol,mxelev_polarris), &
                stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG opt_polarris: allocation error bundled array. ')
    endif

    !
    ! open and read polarris LUT
    !
    call open_polarris_lut(fname_out)


    !
    ! P3 ice LUT
    !
    IF_P3: if ( trim(cloud_microphysics) == 'P3' ) then

      if( .not. allocated(bmueller_array_p3) ) then
        allocate( bmueller_array_p3  (4,4,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
                  fscat_real_array_p3(2,2,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
                  fscat_imag_array_p3(2,2,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
                  vdop_array_p3      (    n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
                  stat=ierr )
      endif

      fname_out = 'POLARRIS_'//trim(cloud_microphysics)//'_ICE-'//trim(char_freq)//'_'//trim(polarris_axisorient_name)//'.nc'
      !
      ! open and read polarris LUT for P3 ICE
      !
      call open_polarris_lut_p3(fname_out)


    endif IF_P3


  endif !lut_polarris

!
! Find instrumental elevation
!
 if( int(elev_polradar(ni)) == undefined_i2 ) then

    call find_elevation(lat_polarris(ni), lon_polarris(ni), elev_inst)
    !
    ! in case cannot find....
    !
    if( elev_inst == undefined ) then
      if(masterproc) print*,'MSG opt_polarris: cannot found instrumental elevation level. '
      if(masterproc) print*,'It is located outside of domain. Force zero.'
      elev_inst = 0.e0
    endif

 else
    elev_inst = elev_polradar(ni) ![km]
 endif

!
! instrument location (ground level)
!

  call latlonhgt2cartesian( lat_polarris(ni)*degrad, lon_polarris(ni)*degrad, &
                            dble(elev_inst), v_inst )  !ground level 
  call latlonhgt2cartesian(lat_polarris(ni)*degrad, lon_polarris(ni)*degrad, 50.d0, v_inst_top) !50km AGL 


!
! Calculate optical properties of each volume element of the cloud.
!

  DOM_2D_LOOP: do j = myj_start, myj_end ; do i = myi_start, myi_end


!
! instrumental viewing angle (north is zero - clock-wise 0-360deg)
! 
 call viewing_angle( lat_polarris(ni)*degrad, lon_polarris(ni)*degrad, dble(elev_inst), &
                     surface(i,j)%lat*degrad, surface(i,j)%lon*degrad, &
                     vang_surf )

 atmos_optpol1d(:)%vang = vang_surf  !all levels


!
! refine vertical layer if needed.  
!
  if( use_crm_level_polarris ) then

     !
     ! pass atmospheric parameters
     !
     atmos_refine1d(1:mxlyr_refine)%t_air    = atmos(i,j,1:mxlyr)%t_air
     atmos_refine1d(1:mxlyr_refine)%press    = atmos(i,j,1:mxlyr)%press
     atmos_refine1d(1:mxlyr_refine)%hgt      = atmos(i,j,1:mxlyr)%hgt
     atmos_refine1d(1:mxlyr_refine)%rh       = atmos(i,j,1:mxlyr)%rh
     atmos_refine1d(1:mxlyr_refine)%u        = atmos(i,j,1:mxlyr)%u
     atmos_refine1d(1:mxlyr_refine)%v        = atmos(i,j,1:mxlyr)%v
     atmos_refine1d(1:mxlyr_refine)%w        = atmos(i,j,1:mxlyr)%w
     atmos_refine1d(1:mxlyr_refine)%rho_dair = atmos(i,j,1:mxlyr)%rho_dair
     atmos_refine1d(1:mxlyr_refine)%rain_rate = atmos(i,j,1:mxlyr)%rain_rate

     !
     ! pass microphysics parameters
     !
     select case(trim(type_microphysics))
     case('GEN')
        q_gen_refine1d (1:mxlyr_refine)  = q_gen (i,j,1:mxlyr)  
        re_gen_refine1d(1:mxlyr_refine)  = re_gen(i,j,1:mxlyr)  

        if( trim(cloud_microphysics) == 'NSSL2MMT' ) then
           rho_graupel_refine1d(1:mxlyr_refine) = rho3d_graupel(i,j,1:mxlyr)
           rho_hail_refine1d   (1:mxlyr_refine) = rho3d_hail(i,j,1:mxlyr)
        endif

        if( trim(cloud_microphysics) == 'P3' ) then
           p3_frime_refine1d  (1:mxlyr_refine) = p3_frime  (i,j,1:mxlyr)
           p3_rhorime_refine1d(1:mxlyr_refine) = p3_rhorime(i,j,1:mxlyr)
           p3_mu_cloud_refine1d(1:mxlyr_refine)= p3_mu_cloud(i,j,1:mxlyr)
           p3_num_ice_refine1d(1:mxlyr_refine) = n_gen(i,j,1:mxlyr)%ice 
        endif

     case('RAMS')
        q_rams_refine1d (1:mxlyr_refine) = q_rams (i,j,1:mxlyr)  
        re_rams_refine1d(1:mxlyr_refine) = re_rams(i,j,1:mxlyr)  
     case('SBM')
        q_sbm_refine1d (1:mxlyr_refine)  = q_sbm (i,j,1:mxlyr)  
        re_sbm_refine1d(1:mxlyr_refine)  = re_sbm(i,j,1:mxlyr)  
        n_sbm_refine1d            (1:mxlyr_refine,1:nbin) = n_sbm            (i,j,1:mxlyr,1:nbin)  
        fmelt_sbm_snow_refine1d   (1:mxlyr_refine,1:nbin) = fmelt_sbm_snow   (i,j,1:mxlyr,1:nbin)  
        fmelt_sbm_graupel_refine1d(1:mxlyr_refine,1:nbin) = fmelt_sbm_graupel(i,j,1:mxlyr,1:nbin) 
        fmelt_sbm_hail_refine1d   (1:mxlyr_refine,1:nbin) = fmelt_sbm_hail   (i,j,1:mxlyr,1:nbin)  
        frime_sbm_snow_refine1d   (1:mxlyr_refine,1:nbin) = frime_sbm_snow   (i,j,1:mxlyr,1:nbin)  
     case default  
     end select 

  else

     !
     ! interpolate atmospheric parameters in consant height level
     !
     call const_hgt_level_atmos( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr), &
                                 mxlyr_refine, hgt_stag,  atmos_refine1d(1:mxlyr_refine) ) 

     ! staggered air temperature
     do k = 2, mxlyr_refine
        stag_t_air(k) = 0.5 * ( atmos_refine1d(k-1)%t_air +  atmos_refine1d(k)%t_air )
     enddo 
     stag_t_air(0) = stag_t_air(1)  !toshii temporal

     select case(trim(type_microphysics))
     case('GEN')

        if( trim(cloud_microphysics) == 'NSSL2MMT' ) then

            call const_hgt_level_gen( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr),         &
                                      q_gen(i,j,1:mxlyr),re_gen(i,j,1:mxlyr), mxlyr_refine, hgt_stag, &
                                      q_gen_refine1d(1:mxlyr_refine), re_gen_refine1d(1:mxlyr_refine),&
                                      rho3d_graupel(i,j,1:mxlyr), rho3d_hail(i,j,1:mxlyr), &
                                      rho_graupel_refine1d(1:mxlyr_refine), rho_hail_refine1d(1:mxlyr_refine) )

        elseif( trim(cloud_microphysics) == 'P3' ) then

            call const_hgt_level_p3( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr),          &
                                      q_gen(i,j,1:mxlyr),re_gen(i,j,1:mxlyr), mxlyr_refine, hgt_stag, &
                                      q_gen_refine1d(1:mxlyr_refine), re_gen_refine1d(1:mxlyr_refine),&
                                      p3_rhorime(i,j,1:mxlyr), p3_frime(i,j,1:mxlyr), &
                                      p3_rhorime_refine1d(1:mxlyr_refine), p3_frime_refine1d(1:mxlyr_refine), &
                                      p3_mu_cloud(i,j,1:mxlyr), n_gen(i,j,1:mxlyr)%ice, &
                                      p3_mu_cloud_refine1d(1:mxlyr_refine), p3_num_ice_refine1d(1:mxlyr_refine)  )


        else  !other general bulk scheme

            call const_hgt_level_gen( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr),         &
                                      q_gen(i,j,1:mxlyr),re_gen(i,j,1:mxlyr), mxlyr_refine, hgt_stag, &
                                      q_gen_refine1d(1:mxlyr_refine), re_gen_refine1d(1:mxlyr_refine) )

        endif

     case('RAMS')

        call const_hgt_level_rams( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr),         &
                                   q_rams(i,j,1:mxlyr),re_rams(i,j,1:mxlyr), mxlyr_refine, hgt_stag, &
                                   q_rams_refine1d(1:mxlyr_refine), re_rams_refine1d(1:mxlyr_refine) )

     case('SBM')

        call const_hgt_level_sbm( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr),         &
                                   q_sbm(i,j,1:mxlyr),re_sbm(i,j,1:mxlyr), mxlyr_refine, hgt_stag, &
                                   q_sbm_refine1d(1:mxlyr_refine), re_sbm_refine1d(1:mxlyr_refine) )

        call const_hgt_level_sbmbin( mxlyr, atmos_stag(i,j,0:mxlyr)%hgt, atmos(i,j,1:mxlyr),  &
                                  n_sbm            (i,j,1:mxlyr,1:nbin) , &
                                  fmelt_sbm_snow   (i,j,1:mxlyr,1:nbin) , &
                                  fmelt_sbm_graupel(i,j,1:mxlyr,1:nbin) , &
                                  fmelt_sbm_hail   (i,j,1:mxlyr,1:nbin) , &
                                  frime_sbm_snow   (i,j,1:mxlyr,1:nbin) , &
                                  mxlyr_refine, hgt_stag, &
                                  n_sbm_refine1d            (1:mxlyr_refine,1:nbin) , &
                                  fmelt_sbm_snow_refine1d   (1:mxlyr_refine,1:nbin) , &
                                  fmelt_sbm_graupel_refine1d(1:mxlyr_refine,1:nbin) , &
                                  fmelt_sbm_hail_refine1d   (1:mxlyr_refine,1:nbin) , &
                                  frime_sbm_snow_refine1d   (1:mxlyr_refine,1:nbin) )                                   

     case default
     end select

  endif

  VERTICAL_LOOP: do k = 1, mxlyr_refine  !either CRM or instrumental levels

    !
    ! initialize output as undef
    !
    atmos_optpol1d(k)%zhh   = undefined
    atmos_optpol1d(k)%zdr   = undefined
    atmos_optpol1d(k)%kdp   = undefined
    atmos_optpol1d(k)%rhohv = undefined
    atmos_optpol1d(k)%vdop  = undefined
    atmos_optpol1d(k)%vrad  = undefined
    atmos_optpol1d(k)%elev  = undefined

    if( ni == 1 .and. polarris_crmout ) then

       !
       ! pass atmospheric parameters to output 3D array
       !
       atmos_geopol(i,j,k)%t_air =  atmos_refine1d(k)%t_air    
       atmos_geopol(i,j,k)%rh    =  atmos_refine1d(k)%rh
       atmos_geopol(i,j,k)%u     =  atmos_refine1d(k)%u
       atmos_geopol(i,j,k)%v     =  atmos_refine1d(k)%v    
       atmos_geopol(i,j,k)%w     =  atmos_refine1d(k)%w    
       atmos_geopol(i,j,k)%rain_rate =  atmos_refine1d(k)%rain_rate  

       !
       ! pass microphysics parameters to output 3D array
       !
       select case(trim(type_microphysics))
       case('GEN')  
          atmos_geopol(i,j,k)%qc = q_gen_refine1d(k)%cloud
          atmos_geopol(i,j,k)%qr = q_gen_refine1d(k)%rain
          atmos_geopol(i,j,k)%qi = q_gen_refine1d(k)%ice
          atmos_geopol(i,j,k)%qs = q_gen_refine1d(k)%snow
          atmos_geopol(i,j,k)%qg = q_gen_refine1d(k)%graupel
          atmos_geopol(i,j,k)%qh = q_gen_refine1d(k)%hail
          atmos_geopol(i,j,k)%rer= re_gen_refine1d(k)%rain
          atmos_geopol(i,j,k)%res= re_gen_refine1d(k)%snow
          atmos_geopol(i,j,k)%reh= re_gen_refine1d(k)%hail

          !toshii add NSSL and P3 output parameter here?

       case('RAMS')

          atmos_geopol(i,j,k)%qc = q_rams_refine1d(k)%cloud1 + q_rams_refine1d(k)%cloud2
          atmos_geopol(i,j,k)%qr = q_rams_refine1d(k)%rain
          atmos_geopol(i,j,k)%qi = q_rams_refine1d(k)%ice1 + q_rams_refine1d(k)%ice2
          atmos_geopol(i,j,k)%qs = q_rams_refine1d(k)%snow
          atmos_geopol(i,j,k)%qg = q_rams_refine1d(k)%graupel
          atmos_geopol(i,j,k)%qh = q_rams_refine1d(k)%hail

       case('SBM')

          liq_tot = sum( n_sbm_refine1d(k,1:nbin)%liq*drad_sbm(1:nbin)%liq*1.e-2 * x_sbm(1:nbin)%liq )
          liq_cloud = sum( n_sbm_refine1d(k,1:18)%liq*drad_sbm(1:18)%liq*1.e-2 * x_sbm(1:18)%liq )
          if( liq_tot > 0. ) then
              sbm_cloud_frac = liq_cloud / liq_tot 
          else
              sbm_cloud_frac = 0.
          endif

          atmos_geopol(i,j,k)%qc = q_sbm_refine1d(k)%liq * sbm_cloud_frac
          atmos_geopol(i,j,k)%qr = q_sbm_refine1d(k)%liq * ( 1.- sbm_cloud_frac ) 
          atmos_geopol(i,j,k)%qi = q_sbm_refine1d(k)%ice_col + q_sbm_refine1d(k)%ice_den + &
                                   q_sbm_refine1d(k)%ice_pla
          if( atmos_geopol(i,j,k)%qi < 0. ) atmos_geopol(i,j,k)%qi = undefined
          atmos_geopol(i,j,k)%qs = q_sbm_refine1d(k)%snow
          atmos_geopol(i,j,k)%qg = q_sbm_refine1d(k)%graupel
          atmos_geopol(i,j,k)%qh = q_sbm_refine1d(k)%hail
       case default
          atmos_geopol(i,j,k)%qc = undefined
          atmos_geopol(i,j,k)%qr = undefined
          atmos_geopol(i,j,k)%qi = undefined
          atmos_geopol(i,j,k)%qs = undefined
          atmos_geopol(i,j,k)%qg = undefined
          atmos_geopol(i,j,k)%qh = undefined
       end select
    endif

    !
    ! specify temperature and pressure 
    !
    tavg = atmos_refine1d(k)%t_air  ! [K]
    tavg_c = tavg - 273.16          ! [degC]
    pres = atmos_refine1d(k)%press  ! [hPa]

    !
    ! derive location of CRM grid in cartesian cordinate. 
    !
    call latlonhgt2cartesian(surface(i,j)%lat*degrad, surface(i,j)%lon*degrad, &
                             dble(atmos_refine1d(k)%hgt), v_grid ) 

    !
    ! derive elevation angle of polarimetric radar
    !
    call elevation_angle( v_inst_top , v_inst, v_grid, elev(1), dist )

    atmos_optpol1d(k)%elev  = elev(1) ! elevation angle [deg]


    !
    ! skip if out of maximum radar coverage or negative elevation angle
    !
    if( dist > coverage_polarris(ni) .or.  &
        atmos_optpol1d(k)%elev < 0. ) then 
        atmos_optpol1d(k)%elev  = undefined
        atmos_optpol1d(k)%vang  = undefined
        cycle VERTICAL_LOOP  !
    endif

    ! also if this is underground or above CRM level, skip loop.
    if( atmos_refine1d(k)%rh == undefined ) cycle VERTICAL_LOOP


    !
    ! Ideintify melting layer and fraction for bulk microphysics (SBM has explicit melting process.) 
    ! --> toshii change this for future....
    !
    if (stag_t_air(k-1) >= const_Kel2Cel .and. stag_t_air(k)<= const_Kel2Cel ) then
        melt_frac_bulk = ( stag_t_air(k-1) - const_Kel2Cel ) / &
                         ( stag_t_air(k-1) - stag_t_air(k) )
    elseif( stag_t_air(k-1) < const_Kel2Cel ) then ! no melt
        melt_frac_bulk = 0. 
    elseif( stag_t_air(k) > const_Kel2Cel   ) then
        melt_frac_bulk = 0.9
    endif

    !
    ! Gas absorption code 
    !
    call gas_absorb(freqcy, tavg, pres, atmos_refine1d(k)%rh, 0., atm_ext)


   !
   ! Condensate particle optical properties 
   !

   mic_select1: select case(trim(type_microphysics))

   case('GEN')

       !
       ! compute optical properties for different species
       !
       if( trim(cloud_microphysics) == 'NSSL2MMT' ) then

           call opt_gen_polarris( ni, freqcy, elev, tavg, melt_frac_bulk, mu_gen, q_gen_refine1d(k), re_gen_refine1d(k), &
                                 bmueller_gen, fscat_gen, vdop_gen , rho_graupel_refine1d(k), rho_hail_refine1d(k) )

       elseif( trim(cloud_microphysics) == 'P3' ) then

           call opt_gen_polarris_p3( ni, freqcy, elev, tavg, melt_frac_bulk, mu_gen, q_gen_refine1d(k), re_gen_refine1d(k), &
                                 p3_mu_cloud_refine1d(k), p3_num_ice_refine1d(k), p3_rhorime_refine1d(k), p3_frime_refine1d(k),  &
                                 bmueller_gen, fscat_gen, vdop_gen )

       else  !other general bulk scheme

           call opt_gen_polarris( ni, freqcy, elev, tavg, melt_frac_bulk, mu_gen, q_gen_refine1d(k), re_gen_refine1d(k), &
                                 bmueller_gen, fscat_gen, vdop_gen )

       endif

       !
       ! Total optical properties
       !

       call total_opt_gen(wave_m, atmos_refine1d(k)%press, atmos_refine1d(k)%w,     & 
                          atm_ext, bmueller_gen, fscat_gen, vdop_gen, &
                          atmos_optpol1d(k) )


   case('RAMS')

     nt = max( 1,min( 31,-nint(  (atmos_refine1d(k)%t_air-273.15)  )))  ! Temperature index--1 to 31 T below 0 deg.C
     ns = max(79,min(100, nint( atmos_refine1d(k)%rh          )))       ! RH index--79 to 100 based on integer RH
     !
     ! compute optical properties for different species of RAMS microphysics
     !
     call opt_rams_polarris( nt, ns, ni, freqcy, elev, tavg, melt_frac_bulk, q_rams_refine1d(k), re_rams_refine1d(k), &
                           bmueller_rams, fscat_rams, vdop_rams )

     !
     ! Total optical properties
     !
     call total_opt_rams(wave_m, atmos_refine1d(k)%rho_dair*1.e-3, atmos_refine1d(k)%w, &
                         atm_ext, bmueller_rams, fscat_rams, vdop_rams, &
                         atmos_optpol1d(k) )


   case('SBM')

     !
     ! compute optical properties for different species of HUCM SBM microphysics
     !
     call opt_sbm_polarris( freqcy, elev, tavg, q_sbm_refine1d(k), n_sbm_refine1d(k,1:nbin), &
                         fmelt_sbm_snow_refine1d(k, 1:nbin), frime_sbm_snow_refine1d(k, 1:nbin), &
                         fmelt_sbm_graupel_refine1d(k, 1:nbin),  fmelt_sbm_hail_refine1d(k, 1:nbin), &
                         bmueller_sbm, fscat_sbm, vdop_sbm )

     !
     ! Total optical properties
     !
     call total_opt_sbm( wave_m, atmos_refine1d(k)%press, atmos_refine1d(k)%w, &
                         atm_ext,  bmueller_sbm, fscat_sbm, vdop_sbm, &
                         atmos_optpol1d(k) )

   case default  !for no-particle simulation (sim_case == 'SOUNDING', 'MLM' or 'LIS')

     atmos_optpol1d(k)%zhh   = undefined
     atmos_optpol1d(k)%zdr   = undefined
     atmos_optpol1d(k)%kdp   = undefined
     atmos_optpol1d(k)%rhohv = undefined
     atmos_optpol1d(k)%vdop  = undefined
     atmos_optpol1d(k)%vrad  = undefined

   end select mic_select1

!
! derive radial velocity (u-eastward, v-northward, w (or vdop) -upward wind)
!

  call derive_radial_velocity( atmos_refine1d(k)%u, atmos_refine1d(k)%v, atmos_optpol1d(k)%vdop, &
                               atmos_optpol1d(k)%elev, &
                               i, j, k, atmos_refine1d(k)%hgt, v_grid, v_inst, atmos_optpol1d(k)%vrad )

!
! check error
!
   call Find_NaN_Inf('MSG: opt_polarris: find NAN in zhh   ',atmos_optpol1d(k)%zhh,i,j,k)
   call Find_NaN_Inf('MSG: opt_polarris: find NAN in zdr   ',atmos_optpol1d(k)%zdr,i,j,k)
   call Find_NaN_Inf('MSG: opt_polarris: find NAN in kdp   ',atmos_optpol1d(k)%kdp,i,j,k)
   call Find_NaN_Inf('MSG: opt_polarris: find NAN in rhohv ',atmos_optpol1d(k)%rhohv,i,j,k)
   call Find_NaN_Inf('MSG: opt_polarris: find NAN in vdop  ',atmos_optpol1d(k)%vdop,i,j,k)
   call Find_NaN_Inf('MSG: opt_polarris: find NAN in vrad  ',atmos_optpol1d(k)%vrad,i,j,k)


  enddo VERTICAL_LOOP

!
! passing to polarris 3D output array 
!
  atmos_optpol(i,j,:,ni)%zhh   = atmos_optpol1d(:)%zhh
  atmos_optpol(i,j,:,ni)%zdr   = atmos_optpol1d(:)%zdr
  atmos_optpol(i,j,:,ni)%kdp   = atmos_optpol1d(:)%kdp
  atmos_optpol(i,j,:,ni)%rhohv = atmos_optpol1d(:)%rhohv
  atmos_optpol(i,j,:,ni)%vdop  = atmos_optpol1d(:)%vdop
  atmos_optpol(i,j,:,ni)%vrad  = atmos_optpol1d(:)%vrad
  atmos_optpol(i,j,:,ni)%elev  = atmos_optpol1d(:)%elev
  atmos_optpol(i,j,:,ni)%vang  = atmos_optpol1d(:)%vang


 enddo ; enddo DOM_2D_LOOP

 enddo INST_LOOP


if( output_grads_polarris ) then
!
! write polarimetric radar observable
!
 call write_out_polarris_grads

endif

!
! write out bundled NetCDF output
!
 if( output_netcdf_polarris ) call write_out_polarris_netcdf


 enddo ensemble_loop

 return
 end subroutine opt_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine open_polarris_lut(fname_out)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Open and read the prepared POLARRIS LUT in NetCDF format. 
! 
! History: 02/2017  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 character(len=100) :: fname_out
 integer :: ncid, varid          !# of variables to inquire  
 integer :: i,j,n

!
! file name
!
#if MPI == 2
if(masterproc) &
#endif
 print*,'MSG open_polarris_lut: reading LUT ',trim(sdsu_dir_sslut)//trim(fname_out)

!
! open netCDF file
!
 call check( nf90_open(trim(sdsu_dir_sslut)//trim(fname_out), nf90_nowrite, ncid) )

!
! 4x4 mueller matrix [-]
!
 call check( nf90_inq_varid(ncid, 'bmueller' , varid ) )  
 call check( nf90_get_var(ncid, varid, bmueller_array ))

!
! 2x2 scattering matrix (real) [-]
!
 call check( nf90_inq_varid(ncid, 'fscat_real' , varid ) )                                          
 call check( nf90_get_var(ncid, varid, fscat_real_array ))

!
! 2x2 scattering matrix (imag) [-]
!
 call check( nf90_inq_varid(ncid, 'fscat_imag' , varid ) )
 call check( nf90_get_var(ncid, varid, fscat_imag_array ))

!
! doppler velocity [m/2]
!
 call check( nf90_inq_varid(ncid, 'vdop' , varid ) )
 call check( nf90_get_var(ncid, varid, vdop_array ))

!
! close netCDF file
!
  call check( nf90_close(ncid) )

 return
 end subroutine open_polarris_lut

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine open_polarris_lut_p3(fname_out)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Open and read the prepared POLARRIS P3 LUT in NetCDF format. 
! 
! History: 06/2019  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 character(len=100) :: fname_out
 integer :: ncid, varid          !# of variables to inquire  
 integer :: i,j,n

!
! file name
!
#if MPI == 2
if(masterproc) &
#endif
 print*,'MSG open_polarris_lut: reading LUT ',trim(sdsu_dir_sslut)//trim(fname_out)

!
! open netCDF file
!
 call check( nf90_open(trim(sdsu_dir_sslut)//trim(fname_out), nf90_nowrite, ncid) )

!
! 4x4 mueller matrix [-]
!
 call check( nf90_inq_varid(ncid, 'bmueller' , varid ) )  
 call check( nf90_get_var(ncid, varid, bmueller_array_p3 ))

!
! 2x2 scattering matrix (real) [-]
!
 call check( nf90_inq_varid(ncid, 'fscat_real' , varid ) )                                          
 call check( nf90_get_var(ncid, varid, fscat_real_array_p3 ))

!
! 2x2 scattering matrix (imag) [-]
!
 call check( nf90_inq_varid(ncid, 'fscat_imag' , varid ) )
 call check( nf90_get_var(ncid, varid, fscat_imag_array_p3 ))

!
! doppler velocity [m/2]
!
 call check( nf90_inq_varid(ncid, 'vdop' , varid ) )
 call check( nf90_get_var(ncid, varid, vdop_array_p3 ))

!
! close netCDF file
!
  call check( nf90_close(ncid) )

 return
 end subroutine open_polarris_lut_p3

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_gen_polarris( nf, freq, elev, tavg, fmelt, mu, q, re, &
                              bmueller_gen, fscat_gen, vdop_gen, rho_graupel, rho_hail )

 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Derive mueller and forward scattering matrix either using LUT or explicit computation. 
!  This routine is for general bulk microphysics microphysics.
!    
! History: 
!  09/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 integer,intent(in)::nf  !frequency id
 real(sdsu_fps),intent(in)::tavg  ! layer temperature
 real(sdsu_fps),intent(inout)::elev(max_ens)  ! elevation angle [deg]
 real(sdsu_fps),intent(in)::freq  ! frequency [GHz]
 real(sdsu_fps),intent(in)::fmelt ! meltfraction 

 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     mu, & ! dispersion parameter
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 type ( particle_gen_r8 ),intent(out) :: & ! particle_gen is defined in module_simulater
     bmueller_gen(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_gen                 ! Doppler velocity [m/s]

 type ( particle_gen_cmplx ),intent(out) :: & ! particle_gen is defined in module_simulater
     fscat_gen(2,2)             ! 2x2 forward scattering matrix 

 real(sdsu_fps),intent(in),optional :: rho_graupel, rho_hail ! effective denisty [kg/m3]


 integer,parameter :: melv = 1
 type ( particle_gen_r8 ) ::  & ! particle_gen is defined in module_simulater
     bmueller(melv,4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop(melv)             ! Doppler velocity [m/s]

 type ( particle_gen_cmplx ) :: & ! particle_gen is defined in module_simulater
     fscat(melv,2,2)             ! 2x2 forward scattering matrix 


 integer :: ispc  !RAMS spicies index
 real(sdsu_fps) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(sdsu_fps) :: cfvt, pwvt ! alpha_vt & beta_vt
 integer :: i,j,k
 real(sdsu_fps) :: rhoe !effective density [kg/m3]

!
! Initialize as zero
!
 do j = 1,4 ; do k = 1,4
  call init_microphysics(bmueller_gen(j,k),'zero')
 enddo ; enddo 

 do j = 1,2 ; do k = 1,2
  call init_microphysics(fscat_gen(j,k),'zero')
 enddo ; enddo

 call init_microphysics(vdop_gen,'zero')


   if ( lut_polarris ) then

       if( present(rho_graupel) .and. present(rho_hail) ) then !case of NSSL 2moment with explicit hail and graupel density

         !
         ! Interpolate LUT: 
         ! LUT is a function of T and Re, and re-scaled for actual q.  
         ! 
         call ss_gen_lut( tavg, elev(1), q, re, &
                          bmueller_gen, fscat_gen, vdop_gen , rho_graupel, rho_hail)


       else
         !
         ! Interpolate LUT: 
         ! LUT is a function of T and Re, and re-scaled for actual q.  
         ! 
         call ss_gen_lut( tavg, elev(1), q, re, &
                          bmueller_gen, fscat_gen, vdop_gen )

       endif

       return

   else  !no use of mie LUT (SLOW: compute extinction by integrating hydrometeor size spectra.)


     if( trim(cloud_microphysics) == 'WDM') then

        call ss_gen_polarris('cloud',freq,tavg,mu_gen%cloud,q%cloud,re%cloud ,rho_gen%cloud ,1.,&   ! rain
                           1, elev, num_axis_gen%cloud, ityp_axis_gen%cloud, axis_ratio_gen%cloud, ityp_dist_gen%cloud,  &
                           mean_orient_gen%cloud, num_stdev_gen%cloud, stdev_orient_gen%cloud, &
                           a_vt%cloud, b_vt%cloud, bmueller%cloud, fscat%cloud, vdop%cloud )
     else

      ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud -> RAMS cloud1
                 cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! but use gen orienting angles
      call ss_rams_polarris('cloud1',freq,tavg,q%cloud,re%cloud,gnu_rams%cloud1, cfvt,pwvt,  &
                          1, elev, num_axis_gen%cloud, ityp_axis_gen%cloud, axis_ratio_gen%cloud, ityp_dist_gen%cloud,  &
                          mean_orient_gen%cloud, num_stdev_gen%cloud, stdev_orient_gen%cloud, &
                          cfmas, pwmas , 1.0,  bmueller%cloud, fscat%cloud, vdop%cloud )
     endif

     call ss_gen_polarris('rain',freq,tavg,mu_gen%rain,q%rain,re%rain ,rho_gen%rain ,1.,&   ! rain
                       1, elev, num_axis_gen%rain, ityp_axis_gen%rain, axis_ratio_gen%rain, ityp_dist_gen%rain,  &
                       mean_orient_gen%rain, num_stdev_gen%rain, stdev_orient_gen%rain, &
                       a_vt%rain, b_vt%rain, bmueller%rain, fscat%rain, vdop%rain )

     if( trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH' .or. &  !
         trim(cloud_microphysics) == 'THOM' ) then !morrison is twomoment
         call ss_gen_polarris('ice',freq,tavg,mu_gen%ice, q%ice,re%ice ,rho_gen%ice ,0.,&   ! ice
                            1, elev, num_axis_gen%ice, ityp_axis_gen%ice, axis_ratio_gen%ice, ityp_dist_gen%ice,  &
                            mean_orient_gen%ice, num_stdev_gen%ice, stdev_orient_gen%ice, &
                            a_vt%ice, b_vt%ice, bmueller%ice, fscat%ice, vdop%ice )
     else
         ispc = 3  ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 --> use RAMS psd
                    cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! but using gen orienting angle.
         call ss_rams_polarris('ice1', freq, tavg, q%ice, re%ice , gnu_rams%ice1, cfvt,pwvt, &  !correction...here..
                             1, elev, num_axis_gen%ice, ityp_axis_gen%ice, axis_ratio_gen%ice, ityp_dist_gen%ice,  &
                             mean_orient_gen%ice, num_stdev_gen%ice, stdev_orient_gen%ice, &
                             cfmas, pwmas , 0.0,  bmueller%ice, fscat%ice, vdop%ice )
     endif

       if( trim(cloud_microphysics) == 'GMP4ICE' ) then
          rhoe = 1.e+3 * max(0.05, min(0.9e0, 142.36 * ( re%snow**(-0.95) ) ) ) ! re-dependent bulk snow density 
       elseif( trim(cloud_microphysics) == 'THOM' ) then
          rhoe = min(917. , 0.13 / (re%snow*2.*1.e-6) ) ! Kg/m^3  (917. to be the same as densice)
       else
          rhoe = rho_gen%snow
       endif

     call ss_gen_polarris('snow',freq,tavg,mu_gen%snow, q%snow,re%snow,rhoe,fmelt, & !snow
                        1, elev, num_axis_gen%snow, ityp_axis_gen%snow, axis_ratio_gen%snow, ityp_dist_gen%snow,  &
                        mean_orient_gen%snow, num_stdev_gen%snow, stdev_orient_gen%snow, &
                        a_vt%snow, b_vt%snow, bmueller%snow, fscat%snow, vdop%snow )

     call ss_gen_polarris('graupel',freq,tavg,mu_gen%graupel, q%graupel,re%graupel, rho_gen%graupel, 0., & !graupel
                        1, elev, num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel,  &
                        mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                        a_vt%graupel, b_vt%graupel,bmueller%graupel, fscat%graupel, vdop%graupel  )

     call ss_gen_polarris('hail', freq, tavg, mu_gen%hail, q%hail, re%hail , rho_gen%hail , 0.,& !hail
                        1, elev, num_axis_gen%hail, ityp_axis_gen%hail, axis_ratio_gen%hail, ityp_dist_gen%hail,  &
                        mean_orient_gen%hail, num_stdev_gen%hail, stdev_orient_gen%hail, &
                        a_vt%hail, b_vt%hail, bmueller%hail, fscat%hail ,vdop%hail )

   endif 


!
! output
!
 bmueller_gen(1:4,1:4) = bmueller(melv,1:4,1:4)

 fscat_gen(1:2,1:2) = fscat(melv,1:2,1:2)

 vdop_gen = vdop(melv)


  

 return
 end subroutine opt_gen_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_gen_polarris_p3( nf, freq, elev, tavg, fmelt, mu, q, re, &
                                 mu_cloud, num_ice, rhorime, frime,  &
                                 bmueller_gen, fscat_gen, vdop_gen )

 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Derive mueller and forward scattering matrix either using LUT or explicit computation. 
!  This routine is for general bulk microphysics microphysics.
!    
! History: 
!  09/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 integer,intent(in)::nf  !frequency id
 real(sdsu_fps),intent(in)::tavg  ! layer temperature
 real(sdsu_fps),intent(inout)::elev(max_ens)  ! elevation angle [deg]
 real(sdsu_fps),intent(in)::freq  ! frequency [GHz]
 real(sdsu_fps),intent(in)::fmelt ! meltfraction 

 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     mu, & ! dispersion parameter
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 real(sdsu_fps),intent(in) :: mu_cloud ! cloud mu [-]
 real(sdsu_fps),intent(in) :: num_ice  ! ice number conc. [#/m3]
 real(sdsu_fps),intent(in) :: rhorime  ! rimed ice density  [kg/m3]
 real(sdsu_fps),intent(in) :: frime    ! mass fraction [-]

 type ( particle_gen_r8 ),intent(out) :: & ! particle_gen is defined in module_simulater
     bmueller_gen(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_gen                 ! Doppler velocity [m/s]

 type ( particle_gen_cmplx ),intent(out) :: & ! particle_gen is defined in module_simulater
     fscat_gen(2,2)             ! 2x2 forward scattering matrix 

 integer,parameter :: melv = 1

 type ( particle_gen_r8 ) ::  & ! particle_gen is defined in module_simulater
     bmueller(melv,4,4),      & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop(melv)                 ! Doppler velocity [m/s]

 type ( particle_gen_cmplx ) :: & ! particle_gen is defined in module_simulater
     fscat(melv,2,2)             ! 2x2 forward scattering matrix 

 integer :: i,j,k

!
! Initialize as zero
!
 do j = 1,4 ; do k = 1,4
  call init_microphysics(bmueller_gen(j,k),'zero')
 enddo ; enddo 

 do j = 1,2 ; do k = 1,2
  call init_microphysics(fscat_gen(j,k),'zero')
 enddo ; enddo

 call init_microphysics(vdop_gen,'zero')


   if ( lut_polarris ) then

         !
         ! Interpolate LUT: 
         ! LUT is a function of T and Re, and re-scaled for actual q.  
         ! 
         call ss_gen_lut_p3( tavg, elev(1), q, re, mu_cloud, num_ice, frime, rhorime, &
                             bmueller_gen, fscat_gen, vdop_gen )


       return

   else  !no use of mie LUT (SLOW: compute extinction by integrating hydrometeor size spectra.)

!toshii make mu_cloud (2~15)
     call ss_gen_polarris('cloud',freq, tavg, mu_cloud, q%cloud, re%cloud, rho_gen%cloud , 1.,& !cloud drop
                        1, elev, num_axis_gen%cloud, ityp_axis_gen%cloud, axis_ratio_gen%cloud, ityp_dist_gen%cloud,  & 
                        mean_orient_gen%cloud, num_stdev_gen%cloud, stdev_orient_gen%cloud, &
                        a_vt%cloud, b_vt%cloud, bmueller%cloud, fscat%cloud, vdop%cloud )


     call ss_gen_polarris('rain',freq,tavg,mu_gen%rain,q%rain,re%rain ,rho_gen%rain ,1.,&   ! rain
                       1, elev, num_axis_gen%rain, ityp_axis_gen%rain, axis_ratio_gen%rain, ityp_dist_gen%rain,  &
                       mean_orient_gen%rain, num_stdev_gen%rain, stdev_orient_gen%rain, &
                       a_vt%rain, b_vt%rain, bmueller%rain, fscat%rain, vdop%rain )

     !toshii make num_ice 
     call ss_p3_polarris(freq, tavg, q%ice, num_ice, frime, rhorime , & !p3 ice
                         1, elev, bmueller%ice, fscat%ice ,vdop%ice )

   endif 


!
! output
!
 bmueller_gen(1:4,1:4) = bmueller(melv,1:4,1:4)

 fscat_gen(1:2,1:2) = fscat(melv,1:2,1:2)

 vdop_gen = vdop(melv)


  

 return
 end subroutine opt_gen_polarris_p3

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine opt_rams_polarris( nt, ns, nf, freq, elev, tavg, fmelt, q, re, &
                            bmueller_rams, fscat_rams, vdop_rams )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Derive mueller and forward scattering matrix either using LUT or explicit computation. 
!  This routine is for general RAMS microphysics microphysics.
!    
! History: 
!  09/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 integer,intent(in):: nt ! Temperature index--1 to 31 T below 0 deg.C
 integer,intent(in):: ns ! RH index--79 to 100 based on integer RH
 integer,intent(in)::nf  !frequency id
 real(sdsu_fps),intent(in)::elev(max_ens)  ! elevation angle [deg]
 real(sdsu_fps),intent(in)::tavg  ! layer temperature
 real(sdsu_fps),intent(in)::freq  ! frequency [GHz]
 real(sdsu_fps),intent(in)::fmelt ! meltfraction 

 type ( particle_rams ),intent(in) :: & ! particle_rams is defined in module_simulater
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 type ( particle_rams_r8 ),intent(out) :: & ! particle_rams is defined in module_simulater
     bmueller_rams(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_rams                 ! Doppler velocity [m/s]

 type ( particle_rams_cmplx ),intent(out) :: & ! particle_rams is defined in module_simulater
     fscat_rams(2,2)         ! 2x2 forward scattering matrix 

 integer,parameter :: melv = 1
 type ( particle_rams_r8 ) ::  & ! particle_rams is defined in module_simulater
     bmueller(melv,4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop(melv)             ! Doppler velocity [m/s]

 type ( particle_rams_cmplx ) :: & ! particle_gen is defined in module_simulater
     fscat(melv,2,2)             ! 2x2 forward scattering matrix 

 integer :: ispc  !RAMS spicies index
 real(sdsu_fps) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(sdsu_fps) :: cfvt, pwvt ! alpha_vt & beta_vt
 integer :: i,j,k

 !
 ! initialize
 !
 do j = 1,4 ; do k = 1,4
  call init_microphysics(bmueller_rams(j,k),'zero')
 enddo ; enddo

 do j = 1,2 ; do k = 1,2
  call init_microphysics(fscat_rams(j,k),'zero')
 enddo ; enddo

 call init_microphysics(vdop_rams,'zero')


 LUT_RAMS: if ( lut_polarris ) then

     stop 'rams polarris: not yet'

   else  !no use of mie LUT (SLOW: compute extinction by integrating hydrometeor size spectra.)

   ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud1
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call ss_rams_polarris('cloud1', freq, tavg, q%cloud1, re%cloud1 , gnu_rams%cloud1 ,cfvt,pwvt, &
                       1, elev, num_axis_rams%cloud1, ityp_axis_rams%cloud1, axis_ratio_rams%cloud1, ityp_dist_rams%cloud1,  &
                       mean_orient_rams%cloud1, num_stdev_rams%cloud1, stdev_orient_rams%cloud1, &
                       cfmas, pwmas , 1.0, bmueller%cloud1, fscat%cloud1, vdop%cloud1 ) 

   ispc = 16; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud2
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call ss_rams_polarris('cloud2', freq, tavg, q%cloud2, re%cloud2 , gnu_rams%cloud2 ,cfvt,pwvt, &
                       1, elev, num_axis_rams%cloud2, ityp_axis_rams%cloud2, axis_ratio_rams%cloud2, ityp_dist_rams%cloud2,  &
                       mean_orient_rams%cloud2, num_stdev_rams%cloud2, stdev_orient_rams%cloud2, &
                       cfmas, pwmas , 1.0,  bmueller%cloud2, fscat%cloud2, vdop%cloud2 )

   ispc = 2 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! rain 
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call ss_rams_polarris('rain', freq, tavg, q%rain, re%rain , gnu_rams%rain, cfvt,pwvt,&
                       1, elev, num_axis_rams%rain, ityp_axis_rams%rain, axis_ratio_rams%rain, ityp_dist_rams%rain,  &
                       mean_orient_rams%rain, num_stdev_rams%rain, stdev_orient_rams%rain, &
                       cfmas, pwmas , 1.0,  bmueller%rain, fscat%rain,vdop%rain )

   ispc = rams_jhabtab(nt,ns,1) ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1
                                 cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call ss_rams_polarris('ice1', freq, tavg, q%ice1, re%ice1 , gnu_rams%ice1 ,cfvt,pwvt, &
                       1, elev, num_axis_rams%ice1, ityp_axis_rams%ice1, axis_ratio_rams%ice1, ityp_dist_rams%ice1,  &
                       mean_orient_rams%ice1, num_stdev_rams%ice1, stdev_orient_rams%ice1, &
                       cfmas, pwmas , 0.0,  bmueller%ice1, fscat%ice1, vdop%ice1 )

   ispc = rams_jhabtab(nt,ns,2) ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2
                                 cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call ss_rams_polarris('ice2', freq, tavg, q%ice2, re%ice2 , gnu_rams%ice2 ,cfvt,pwvt, &
                       1, elev, num_axis_rams%ice2, ityp_axis_rams%ice2, axis_ratio_rams%ice2, ityp_dist_rams%ice2,  &
                       mean_orient_rams%ice2, num_stdev_rams%ice2, stdev_orient_rams%ice2, &
                       cfmas, pwmas , 0.0, bmueller%ice2, fscat%ice2, vdop%ice2 )
 
   ispc = 5 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! snow aggregate
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call ss_rams_polarris('snow', freq, tavg, q%snow, re%snow , gnu_rams%snow,cfvt,pwvt, &
                       1, elev, num_axis_rams%snow, ityp_axis_rams%snow, axis_ratio_rams%snow, ityp_dist_rams%snow,  &
                       mean_orient_rams%snow, num_stdev_rams%snow, stdev_orient_rams%snow, &
                       cfmas, pwmas, fmelt, bmueller%snow, fscat%snow, vdop%snow )
  
   ispc = 6 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! graupel
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call ss_rams_polarris('graupel',freq,tavg,q%graupel,re%graupel,gnu_rams%graupel,cfvt,pwvt, &
                       1, elev, num_axis_rams%graupel, ityp_axis_rams%graupel, axis_ratio_rams%graupel, ityp_dist_rams%graupel,  &
                       mean_orient_rams%graupel, num_stdev_rams%graupel, stdev_orient_rams%graupel, &
                       cfmas, pwmas , 0.0, bmueller%graupel, fscat%graupel, vdop%graupel )

   ispc = 7 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! hail
              cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
   call ss_rams_polarris( 'hail',freq, tavg, q%hail, re%hail, gnu_rams%hail ,cfvt,pwvt, &
                        1, elev, num_axis_rams%hail, ityp_axis_rams%hail, axis_ratio_rams%hail, ityp_dist_rams%hail,  &
                        mean_orient_rams%hail, num_stdev_rams%hail, stdev_orient_rams%hail, &
                        cfmas, pwmas , 0.0,  bmueller%hail, fscat%hail, vdop%hail )

   endif LUT_RAMS


!
! output
!
 bmueller_rams(1:4,1:4) = bmueller(melv,1:4,1:4)

 fscat_rams(1:2,1:2) = fscat(melv,1:2,1:2)

 vdop_rams = vdop(melv)




 return
 end subroutine  opt_rams_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine opt_sbm_polarris( freq, elev, tavg, q, n, fmelt, frime, fmelt_g, fmelt_h, &
                           bmueller_sbm, fscat_sbm, vdop_sbm)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Derive mueller and forward scattering matrix either using LUT or explicit computation. 
!  This routine is for HUCM spectra bin microphysics.
!    
! History: 
!  09/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(in)::tavg  ! layer temperature [K]
 real(sdsu_fps),intent(inout)::elev(max_ens)  ! elevation angle [deg]
 real(sdsu_fps),intent(in)::freq  ! frequency [GHz]

 type ( particle_sbm ),intent(in) :: & ! particle_sbm is defined in module_simulater
     q           , & ! bulk mixing ratio [g/m3] 
     n   (1:nbin)    ! particle # concentration [1/m4]
 real(sdsu_fps),intent(in) :: fmelt(1:nbin)    ! mass melting fraction of snow [-]
 real(sdsu_fps),intent(in) :: frime(1:nbin)    ! mass rime fraction of snow [-]
 real(sdsu_fps),intent(in) :: fmelt_g(1:nbin)  ! mass melting fraction of graupel[-]
 real(sdsu_fps),intent(in) :: fmelt_h(1:nbin)  ! mass melting fraction of hail [-]

 type ( particle_sbm_r8 ),intent(out) :: & ! particle_sbm is defined in module_simulater
     bmueller_sbm(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_sbm              ! Doppler velocity [m/s]

 type ( particle_sbm_cmplx ),intent(out) :: & ! particle_sbm is defined in module_simulater
     fscat_sbm(2,2)             ! 2x2 forward scattering matrix 

 integer,parameter :: melv = 1
 type ( particle_sbm_r8 ) ::  & ! particle_sbm is defined in module_simulater
     bmueller(melv,4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop(melv)             ! Doppler velocity [m/s]

 type ( particle_sbm_cmplx ) :: & ! particle_sbm is defined in module_simulater
     fscat(melv,2,2)             ! 2x2 forward scattering matrix 

 real(sdsu_fps) :: frime_fake(1:nbin)  !fake 

 integer :: j,k

!
! Initialize as zero
!
 do j = 1,4 ; do k = 1,4
  call init_microphysics(bmueller_sbm(j,k),'zero')
 enddo ; enddo

 do j = 1,2 ; do k = 1,2
  call init_microphysics(fscat_sbm(j,k),'zero')
 enddo ; enddo

 call init_microphysics(vdop_sbm,'zero')

!
! integrate single-scattering property by size integration
!
 if ( lut_polarris ) then

    !
    ! Interpolate LUT: 
    ! LUT is a function of T and Re, and re-scaled for actual q.  
    ! toshii --tavg --> 273.
    !
    call ss_sbm_lut( tavg, elev(1), q, n, fmelt, frime, fmelt_g, fmelt_h, & 
                     bmueller_sbm, fscat_sbm, vdop_sbm )

    return

 else  !diret solution of T-Mueller matrix

    frime_fake(1:nbin) = 0. !fake 

    call ss_sbm_polarris('liq    ',1,nbin,freq,tavg,den_sbm(1:nbin)%liq,rad_sbm(1:nbin)%liq,drad_sbm(1:nbin)%liq,&
                         vt_sbm(1:nbin)%liq, q%liq,n(1:nbin)%liq,fmelt(1:nbin)      ,frime_fake(1:nbin), &
                         1, elev, num_axis_sbm%liq, ityp_axis_sbm%liq, axis_ratio_sbm%liq, ityp_dist_sbm%liq,  &
                         mean_orient_sbm%liq, num_stdev_sbm%liq, stdev_orient_sbm%liq, &
                         bmueller%liq, fscat%liq, vdop%liq ) 

    call ss_sbm_polarris('ice_col',1,nbin,freq,tavg,den_sbm(1:nbin)%ice_col,rad_sbm(1:nbin)%ice_col,drad_sbm(1:nbin)%ice_col,&
                         vt_sbm(1:nbin)%ice_col, q%ice_col,n(1:nbin)%ice_col,fmelt(1:nbin)      ,frime_fake(1:nbin),&
                         1, elev, num_axis_sbm%ice_col, ityp_axis_sbm%ice_col, axis_ratio_sbm%ice_col, ityp_dist_sbm%ice_col,  &
                         mean_orient_sbm%ice_col, num_stdev_sbm%ice_col, stdev_orient_sbm%ice_col, &
                         bmueller%ice_col, fscat%ice_col, vdop%ice_col  ) 

    call ss_sbm_polarris('ice_pla',1,nbin,freq,tavg,den_sbm(1:nbin)%ice_pla,rad_sbm(1:nbin)%ice_pla,drad_sbm(1:nbin)%ice_pla,&
                         vt_sbm(1:nbin)%ice_pla, q%ice_pla,n(1:nbin)%ice_pla,fmelt(1:nbin)     ,frime_fake(1:nbin),&
                         1, elev, num_axis_sbm%ice_pla, ityp_axis_sbm%ice_pla, axis_ratio_sbm%ice_pla, ityp_dist_sbm%ice_pla,  &
                         mean_orient_sbm%ice_pla, num_stdev_sbm%ice_pla, stdev_orient_sbm%ice_pla, &
                         bmueller%ice_pla, fscat%ice_pla, vdop%ice_pla  ) 

    call ss_sbm_polarris('ice_den',1,nbin,freq,tavg,den_sbm(1:nbin)%ice_den,rad_sbm(1:nbin)%ice_den,drad_sbm(1:nbin)%ice_den,&
                         vt_sbm(1:nbin)%ice_den, q%ice_den,n(1:nbin)%ice_den,fmelt(1:nbin)      ,frime(1:nbin),&
                         1, elev, num_axis_sbm%ice_den, ityp_axis_sbm%ice_den, axis_ratio_sbm%ice_den, ityp_dist_sbm%ice_den,  &
                         mean_orient_sbm%ice_den, num_stdev_sbm%ice_den, stdev_orient_sbm%ice_den, &
                         bmueller%ice_den, fscat%ice_den, vdop%ice_den  ) 

    call ss_sbm_polarris('snow   ',1,nbin,freq,tavg,den_sbm(1:nbin)%snow,rad_sbm(1:nbin)%snow,drad_sbm(1:nbin)%snow,&
                         vt_sbm(1:nbin)%snow, q%snow,n(1:nbin)%snow,fmelt(1:nbin)      ,frime(1:nbin),&
                         1, elev,num_axis_sbm%snow, ityp_axis_sbm%snow, axis_ratio_sbm%snow, ityp_dist_sbm%snow,  &
                         mean_orient_sbm%snow, num_stdev_sbm%snow, stdev_orient_sbm%snow, &
                         bmueller%snow, fscat%snow, vdop%snow  ) 

    call ss_sbm_polarris('graupel',1,nbin,freq,tavg,den_sbm(1:nbin)%graupel,rad_sbm(1:nbin)%graupel,drad_sbm(1:nbin)%graupel,&
                         vt_sbm(1:nbin)%graupel, q%graupel,n(1:nbin)%graupel,fmelt_g(1:nbin)     ,frime_fake(1:nbin),&
                         1, elev,num_axis_sbm%graupel, ityp_axis_sbm%graupel, axis_ratio_sbm%graupel, ityp_dist_sbm%graupel,  &
                         mean_orient_sbm%graupel, num_stdev_sbm%graupel, stdev_orient_sbm%graupel, &
                         bmueller%graupel, fscat%graupel, vdop%graupel ) ! graupel

    call ss_sbm_polarris('hail   ',1,nbin,freq,tavg,den_sbm(1:nbin)%hail,rad_sbm(1:nbin)%hail,drad_sbm(1:nbin)%hail,&
                         vt_sbm(1:nbin)%hail, q%hail,n(1:nbin)%hail,fmelt_h(1:nbin)    ,frime_fake(1:nbin),&
                         1, elev,num_axis_sbm%hail, ityp_axis_sbm%hail, axis_ratio_sbm%hail, ityp_dist_sbm%hail,  &
                         mean_orient_sbm%hail, num_stdev_sbm%hail, stdev_orient_sbm%hail, &
                         bmueller%hail, fscat%hail, vdop%hail ) ! hail

    !
    ! output
    !
    bmueller_sbm(1:4,1:4) = bmueller(melv,1:4,1:4)

    fscat_sbm(1:2,1:2) = fscat(melv,1:2,1:2)

    vdop_sbm = vdop(melv)


 endif


 return
 end subroutine opt_sbm_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine ss_gen_lut( tair, elev, q, re,  &
                        bmueller_gen, fscat_gen, vdop_gen, rho_g, rho_h  )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the mueller scattering matrix for
! given elevation angle, temperature , and effective radius using POLARRIS LUTs (per unit q), 
! and values are re-scaled to actual qm (hydrometeor mass mixing ratio).
! 
! History:
!   02/2017  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in)::tair  ! layer temperature [K]
 real(sdsu_fps),intent(inout)::elev  ! elevation angle [deg]

 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]


 type ( particle_gen_r8 ),intent(out) :: & ! particle_gen is defined in module_simulater
     bmueller_gen(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_gen              ! Doppler velocity [m/s]

 type ( particle_gen_cmplx ),intent(out) :: & ! particle_gen is defined in module_simulater
     fscat_gen(2,2)             ! 2x2 forward scattering matrix 

 real(sdsu_fps),intent(in),optional :: rho_g, rho_h  !graupel and hail density [kg/m3] for NSSL 2moment

 integer :: i_elev
 real(sdsu_fps) :: wgt1_elev, wgt2_elev  !weights for elevation

 integer :: it, i_tair_liq
 real(sdsu_fps) :: wgt1_tair_liq, wgt2_tair_liq  !weights for air temperature (liquid)

 integer :: i_tair_ice
 real(sdsu_fps) :: wgt1_tair_ice, wgt2_tair_ice  !weights for air temperature  (ice)

 integer :: i_re
 real(sdsu_fps) :: wgt1_re, wgt2_re  !weights for effective radius (ice)
 real(sdsu_fps) :: ref

 integer :: ispc
 real(sdsu_fps) :: qm

!
! interpolate elevation angle
! 
 if( elev < 0. ) then
    print*,'MSG ss_gen_lut: elevatation angle is < 0. Correct --> 0.'
    elev =  0.
 elseif( elev > 90. ) then
    print*,'MSG ss_gen_lut: elevatation angle is > 90. Correct --> 90.'
    elev = 90.
 endif

 i_elev = int(elev) + 2
 wgt1_elev = ( elev_polarris( i_elev ) - elev )
 wgt2_elev = 1. - wgt1_elev

!
! interpolate temperature (liq)
!
 if( tair < pts_temp_w_pol(1) ) then
   i_tair_liq = 2
   wgt1_tair_liq = 1. ; wgt2_tair_liq = 0. 
 elseif( tair >= pts_temp_w_pol( mxpts_temp_pol ) ) then
   i_tair_liq = mxpts_temp_pol
   wgt1_tair_liq = 0. ; wgt2_tair_liq = 1.
 else  
   do it = 2, mxpts_temp_pol 
      if( pts_temp_w_pol(it-1) <= tair .and. tair < pts_temp_w_pol(it) ) then
          i_tair_liq = it
          wgt1_tair_liq = ( pts_temp_w_pol(it) - tair ) / &
                          ( pts_temp_w_pol(it) - pts_temp_w_pol(it-1) )  
          wgt2_tair_liq = 1. - wgt1_tair_liq
      endif
   enddo
 endif

!
! interpolate temperature (ice)
!
 if( tair < pts_temp_i_pol(1) ) then
   i_tair_ice = 2
   wgt1_tair_ice = 1. ; wgt2_tair_ice = 0.
 elseif( tair >= pts_temp_i_pol( mxpts_temp_pol ) ) then
   i_tair_ice = mxpts_temp_pol
   wgt1_tair_ice = 0. ; wgt2_tair_ice = 1.
 else
   do it = 2, mxpts_temp_pol
      if( pts_temp_i_pol(it-1) <= tair .and. tair < pts_temp_i_pol(it) ) then
          i_tair_ice = it
          wgt1_tair_ice = ( pts_temp_i_pol(it) - tair ) / &
                          ( pts_temp_i_pol(it) - pts_temp_i_pol(it-1) )
          wgt2_tair_ice = 1. - wgt1_tair_ice
      endif
   enddo
 endif


!
! cloud 
!
 cloud_present: if( re%cloud == 0. .or. q%cloud == 0. ) then
    bmueller_gen%cloud = 0. ;  fscat_gen%cloud = 0. ;  vdop_gen%cloud = 0. 
 else
   ispc = 1 ; ref = re%cloud ; qm = q%cloud
   if( ref < pts_re_small(1) ) then
     i_re = 2
     wgt1_re = 1. ; wgt2_re = 0.
   elseif( ref >= pts_re_small( mxpts_re_pol ) ) then
     i_re = mxpts_re_pol
     wgt1_re = 0. ; wgt2_re = 1.
   else
     do it = 2, mxpts_re_pol
        if( pts_re_small(it-1) <= ref .and. ref < pts_re_small(it) ) then
            i_re = it
            wgt1_re = ( pts_re_small(it) - ref ) / &
                      ( pts_re_small(it) - pts_re_small(it-1) )
            wgt2_re = 1. - wgt1_re
        endif
     enddo
   endif


  call interpolate_polarris_lut( qm, i_re, ispc, i_tair_liq, i_elev, &
                        wgt1_re, wgt2_re, wgt1_tair_liq, wgt2_tair_liq, wgt1_elev, wgt2_elev, &
                        bmueller_gen%cloud, fscat_gen%cloud, vdop_gen%cloud )

 endif cloud_present

!
! rain
!
 rain_present: if( re%rain == 0. .or. q%rain == 0. ) then
    bmueller_gen%rain = 0. ;  fscat_gen%rain = 0. ;  vdop_gen%rain = 0.
 else
   ispc = 2 ; ref = re%rain ; qm = q%rain
   if( ref < pts_re_large(1) ) then
     i_re = 2
     wgt1_re = 1. ; wgt2_re = 0.
   elseif( ref >= pts_re_large( mxpts_re_pol ) ) then
     i_re = mxpts_re_pol
     wgt1_re = 0. ; wgt2_re = 1.
   else
     do it = 2, mxpts_re_pol
        if( pts_re_large(it-1) <= ref .and. ref < pts_re_large(it) ) then
            i_re = it
            wgt1_re = ( pts_re_large(it) - ref ) / &
                      ( pts_re_large(it) - pts_re_large(it-1) )
            wgt2_re = 1. - wgt1_re
        endif
     enddo
   endif


  call interpolate_polarris_lut( qm, i_re, ispc, i_tair_liq, i_elev, &
                        wgt1_re, wgt2_re, wgt1_tair_liq, wgt2_tair_liq, wgt1_elev, wgt2_elev, &
                        bmueller_gen%rain, fscat_gen%rain, vdop_gen%rain )

 endif rain_present


!
! ice
!
 ice_present: if( re%ice == 0. .or. q%ice == 0. ) then
    bmueller_gen%ice = 0. ;  fscat_gen%ice = 0. ;  vdop_gen%ice = 0.
 else
   ispc = 3 ; ref = re%ice ; qm = q%ice
   if( ref < pts_re_small(1) ) then
     i_re = 2
     wgt1_re = 1. ; wgt2_re = 0.
   elseif( ref >= pts_re_small( mxpts_re_pol ) ) then
     i_re = mxpts_re_pol
     wgt1_re = 0. ; wgt2_re = 1.
   else
     do it = 2, mxpts_re_pol
        if( pts_re_small(it-1) <= ref .and. ref < pts_re_small(it) ) then
            i_re = it
            wgt1_re = ( pts_re_small(it) - ref ) / &
                      ( pts_re_small(it) - pts_re_small(it-1) )
            wgt2_re = 1. - wgt1_re
        endif
     enddo
   endif


  call interpolate_polarris_lut( qm, i_re, ispc, i_tair_ice, i_elev, &
                        wgt1_re, wgt2_re, wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                        bmueller_gen%ice, fscat_gen%ice, vdop_gen%ice )

 endif ice_present

!
! snow
!
 snow_present: if( re%snow == 0. .or. q%snow == 0. ) then
     bmueller_gen%snow = 0. ;  fscat_gen%snow = 0. ;  vdop_gen%snow = 0.
 else
   ispc = 4 ; ref = re%snow ; qm = q%snow
   if( ref < pts_re_large(1) ) then
     i_re = 2
     wgt1_re = 1. ; wgt2_re = 0.
   elseif( ref >= pts_re_large( mxpts_re_pol ) ) then
     i_re = mxpts_re_pol
     wgt1_re = 0. ; wgt2_re = 1.
   else
     do it = 2, mxpts_re_pol
        if( pts_re_large(it-1) <= ref .and. ref < pts_re_large(it) ) then
            i_re = it
            wgt1_re = ( pts_re_large(it) - ref ) / &
                      ( pts_re_large(it) - pts_re_large(it-1) )
            wgt2_re = 1. - wgt1_re
        endif
     enddo
   endif

  call interpolate_polarris_lut( qm, i_re, ispc, i_tair_ice, i_elev, &
                        wgt1_re, wgt2_re, wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                        bmueller_gen%snow, fscat_gen%snow, vdop_gen%snow )
 endif snow_present




!
! graupel
!

 graupel_present: if( re%graupel == 0. .or. q%graupel == 0. ) then
    bmueller_gen%graupel = 0. ;  fscat_gen%graupel = 0. ;  vdop_gen%graupel = 0.
 else
   ispc = 5 ; ref = re%graupel ; qm = q%graupel
   if( ref < pts_re_large(1) ) then
     i_re = 2
     wgt1_re = 1. ; wgt2_re = 0.
   elseif( ref >= pts_re_large( mxpts_re_pol ) ) then
     i_re = mxpts_re_pol
     wgt1_re = 0. ; wgt2_re = 1.
   else
     do it = 2, mxpts_re_pol
        if( pts_re_large(it-1) <= ref .and. ref < pts_re_large(it) ) then
            i_re = it
            wgt1_re = ( pts_re_large(it) - ref ) / &
                      ( pts_re_large(it) - pts_re_large(it-1) )
            wgt2_re = 1. - wgt1_re
        endif
     enddo
   endif

   if( present(rho_g) ) then !NSSL 2MMT

    call interpolate_polarris_lut_rho( 'graupel', qm, rho_g, i_re, i_tair_ice, i_elev, &
                          wgt1_re, wgt2_re, wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                          bmueller_gen%graupel, fscat_gen%graupel, vdop_gen%graupel )

   else !general scheme

    call interpolate_polarris_lut( qm, i_re, ispc, i_tair_ice, i_elev, &
                          wgt1_re, wgt2_re, wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                          bmueller_gen%graupel, fscat_gen%graupel, vdop_gen%graupel )

   endif
 endif graupel_present


!
! hail
!

 hail_present: if( re%hail  == 0. .or. q%hail == 0. ) then
    bmueller_gen%hail = 0. ;  fscat_gen%hail = 0. ;  vdop_gen%hail = 0.
 else
   ispc = 6 ; ref = re%hail ; qm = q%hail  
   if( ref < pts_re_large(1) ) then
     i_re = 2
     wgt1_re = 1. ; wgt2_re = 0.
   elseif( ref >= pts_re_large( mxpts_re_pol ) ) then
     i_re = mxpts_re_pol
     wgt1_re = 0. ; wgt2_re = 1.
   else
     do it = 2, mxpts_re_pol
        if( pts_re_large(it-1) <= ref .and. ref < pts_re_large(it) ) then
            i_re = it
            wgt1_re = ( pts_re_large(it) - ref ) / &
                      ( pts_re_large(it) - pts_re_large(it-1) )
            wgt2_re = 1. - wgt1_re
        endif
     enddo
   endif


   if( present(rho_h) ) then !NSS 2MMT

    call interpolate_polarris_lut_rho( 'hail', qm, rho_h, i_re, i_tair_ice, i_elev, &
                          wgt1_re, wgt2_re, wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                          bmueller_gen%hail, fscat_gen%hail, vdop_gen%hail )


   else !general scheme

    call interpolate_polarris_lut( qm, i_re, ispc, i_tair_ice, i_elev, &
                          wgt1_re, wgt2_re, wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                          bmueller_gen%hail, fscat_gen%hail, vdop_gen%hail )

   endif

 endif hail_present


 return
 end subroutine ss_gen_lut

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine ss_gen_lut_p3( tair, elev, q, re, mu_cloud, num_ice, frime_in, rhorime_in, &
                           bmueller_gen, fscat_gen, vdop_gen  )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the mueller scattering matrix for P3 ice microphysics with
! given elevation angle, temperature , and normzlied mass, 
! rime ice fraction, and rime ice density using POLARRIS LUTs (per unit number), 
! and scattering matrix is re-scaled to actual number conc.
! 
! History:
!   07/2019  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in)::tair  ! layer temperature [K]
 real(sdsu_fps),intent(inout)::elev  ! elevation angle [deg]

 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 real(sdsu_fps),intent(in) :: mu_cloud ! cloud mu [-]
 real(sdsu_fps),intent(in) :: num_ice  ! ice number conc. [#/m3]
 real(sdsu_fps),intent(in) :: rhorime_in  ! rimed ice density  [kg/m3]
 real(sdsu_fps),intent(in) :: frime_in    ! mass fraction [-]


 type ( particle_gen_r8 ),intent(out) :: & ! particle_gen is defined in module_simulater
     bmueller_gen(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_gen              ! Doppler velocity [m/s]

 type ( particle_gen_cmplx ),intent(out) :: & ! particle_gen is defined in module_simulater
     fscat_gen(2,2)             ! 2x2 forward scattering matrix 


 integer :: i_elev
 real(sdsu_fps) :: wgt1_elev, wgt2_elev  !weights for elevation

 integer :: it, i_tair_liq
 real(sdsu_fps) :: wgt1_tair_liq, wgt2_tair_liq  !weights for air temperature (liquid)

 integer :: i_tair_ice
 real(sdsu_fps) :: wgt1_tair_ice, wgt2_tair_ice  !weights for air temperature  (ice)

 integer :: i_re
 real(sdsu_fps) :: wgt1_re, wgt2_re  !weights for effective radius (ice)
 real(sdsu_fps) :: ref

real(sdsu_fps) :: qnorm    ! normzlied q [kg]
 real(sdsu_fps) :: rhorime  ! rimed ice density  [kg/m3]
 real(sdsu_fps) :: frime    ! mass fraction [-]

 integer :: ispc
 real(sdsu_fps) :: qm

 integer :: irho_found
 integer :: ifr_found
 integer :: iqn_found
 real(sdsu_fps) :: wgt1_rhor, wgt2_rhor
 real(sdsu_fps) :: wgt1_fr  , wgt2_fr
 real(sdsu_fps) :: wgt1_qn  , wgt2_qn
 integer            :: i_rhor          ! index for rho_rime loop                [1 .. n_rhor]
 integer            :: i_Fr            ! index for rime-mass-fraction loop      [1 .. n_Fr]
 integer            :: i_Qnorm         ! index for normalized qitot loop        [1 .. n_Qnorm]



!
! interpolate elevation angle
! 
 if( elev < 0. ) then
    print*,'MSG ss_gen_lut: elevatation angle is < 0. Correct --> 0.'
    elev =  0.
 elseif( elev > 90. ) then
    print*,'MSG ss_gen_lut: elevatation angle is > 90. Correct --> 90.'
    elev = 90.
 endif

 i_elev = int(elev) + 2
 wgt1_elev = ( elev_polarris( i_elev ) - elev )
 wgt2_elev = 1. - wgt1_elev

!
! interpolate temperature (liq)
!
 if( tair < pts_temp_w_pol(1) ) then
   i_tair_liq = 2
   wgt1_tair_liq = 1. ; wgt2_tair_liq = 0. 
 elseif( tair >= pts_temp_w_pol( mxpts_temp_pol ) ) then
   i_tair_liq = mxpts_temp_pol
   wgt1_tair_liq = 0. ; wgt2_tair_liq = 1.
 else  
   do it = 2, mxpts_temp_pol 
      if( pts_temp_w_pol(it-1) <= tair .and. tair < pts_temp_w_pol(it) ) then
          i_tair_liq = it
          wgt1_tair_liq = ( pts_temp_w_pol(it) - tair ) / &
                          ( pts_temp_w_pol(it) - pts_temp_w_pol(it-1) )  
          wgt2_tair_liq = 1. - wgt1_tair_liq
      endif
   enddo
 endif

!
! interpolate temperature (ice)
!
 if( tair < pts_temp_i_pol(1) ) then
   i_tair_ice = 2
   wgt1_tair_ice = 1. ; wgt2_tair_ice = 0.
 elseif( tair >= pts_temp_i_pol( mxpts_temp_pol ) ) then
   i_tair_ice = mxpts_temp_pol
   wgt1_tair_ice = 0. ; wgt2_tair_ice = 1.
 else
   do it = 2, mxpts_temp_pol
      if( pts_temp_i_pol(it-1) <= tair .and. tair < pts_temp_i_pol(it) ) then
          i_tair_ice = it
          wgt1_tair_ice = ( pts_temp_i_pol(it) - tair ) / &
                          ( pts_temp_i_pol(it) - pts_temp_i_pol(it-1) )
          wgt2_tair_ice = 1. - wgt1_tair_ice
      endif
   enddo
 endif


!
! cloud 
!
 cloud_present: if( re%cloud == 0. .or. q%cloud == 0. ) then
    bmueller_gen%cloud = 0. ;  fscat_gen%cloud = 0. ;  vdop_gen%cloud = 0. 
 else
   ispc = 1 ; ref = re%cloud ; qm = q%cloud
   if( ref < pts_re_small(1) ) then
     i_re = 2
     wgt1_re = 1. ; wgt2_re = 0.
   elseif( ref >= pts_re_small( mxpts_re_pol ) ) then
     i_re = mxpts_re_pol
     wgt1_re = 0. ; wgt2_re = 1.
   else
     do it = 2, mxpts_re_pol
        if( pts_re_small(it-1) <= ref .and. ref < pts_re_small(it) ) then
            i_re = it
            wgt1_re = ( pts_re_small(it) - ref ) / &
                      ( pts_re_small(it) - pts_re_small(it-1) )
            wgt2_re = 1. - wgt1_re
        endif
     enddo
   endif


  call interpolate_polarris_lut( qm, i_re, ispc, i_tair_liq, i_elev, &
                        wgt1_re, wgt2_re, wgt1_tair_liq, wgt2_tair_liq, wgt1_elev, wgt2_elev, &
                        bmueller_gen%cloud, fscat_gen%cloud, vdop_gen%cloud )

 endif cloud_present

!
! rain
!
 rain_present: if( re%rain == 0. .or. q%rain == 0. ) then
    bmueller_gen%rain = 0. ;  fscat_gen%rain = 0. ;  vdop_gen%rain = 0.
 else
   ispc = 2 ; ref = re%rain ; qm = q%rain
   if( ref < pts_re_large(1) ) then
     i_re = 2
     wgt1_re = 1. ; wgt2_re = 0.
   elseif( ref >= pts_re_large( mxpts_re_pol ) ) then
     i_re = mxpts_re_pol
     wgt1_re = 0. ; wgt2_re = 1.
   else
     do it = 2, mxpts_re_pol
        if( pts_re_large(it-1) <= ref .and. ref < pts_re_large(it) ) then
            i_re = it
            wgt1_re = ( pts_re_large(it) - ref ) / &
                      ( pts_re_large(it) - pts_re_large(it-1) )
            wgt2_re = 1. - wgt1_re
        endif
     enddo
   endif


  call interpolate_polarris_lut( qm, i_re, ispc, i_tair_liq, i_elev, &
                        wgt1_re, wgt2_re, wgt1_tair_liq, wgt2_tair_liq, wgt1_elev, wgt2_elev, &
                        bmueller_gen%rain, fscat_gen%rain, vdop_gen%rain )

 endif rain_present


!
! ice
!
 ice_present: if( q%ice == 0. .or. num_ice == 0. ) then
      bmueller_gen%ice=0. ; fscat_gen%ice=0. ; vdop_gen%ice=0.
 else
  
  qnorm  = q%ice * 1.e-3 / num_ice  ! [kg]

   !
   ! qnorm weight
   !
   if( qnorm <=  Qn(1) ) then ! Qtot/Ntot [kg]
      wgt1_qn = 1. ; wgt2_qn = 0. ; iqn_found = 1
   elseif( qnorm > Qn(n_Qnorm) )   then
      wgt1_qn = 0. ; wgt2_qn = 1. ; iqn_found = n_Qnorm - 1
   else
      qn_loop: do i_Qnorm = 1, n_Qnorm-1
         wgt1_qn = 1. ; wgt2_qn = 0. ; iqn_found = 1
         if( qnorm > Qn( i_Qnorm ) .and. qnorm <= Qn( i_Qnorm+1 ) ) then
             wgt1_qn = ( Qn( i_Qnorm+1 ) - qnorm ) / ( Qn( i_Qnorm+1 ) - Qn( i_Qnorm ) )
             wgt2_qn = 1. - wgt1_qn
             iqn_found = i_Qnorm
             if( wgt1_qn > 1.0 .or. wgt1_qn < 0. ) stop 'wrong wgt1_qn'
             exit qn_loop
         endif
      enddo qn_loop
   endif

   !
   ! density weight
   !
   rhorime = min( max(rhorime_in, 50.), 900.)
   if( rhorime == 50. ) then
      wgt1_rhor = 1. ; wgt2_rhor = 0. ; irho_found = 1
   else
      rhor_loop: do i_rhor = 1, n_rhor-1
         wgt1_rhor = 1. ; wgt2_rhor = 0. ; irho_found = 1
         if( rhorime > rhor( i_rhor ) .and. rhorime <= rhor( i_rhor+1 ) ) then
             wgt1_rhor = ( rhor( i_rhor+1 ) - rhorime ) / ( rhor( i_rhor+1 ) - rhor( i_rhor ) )
             wgt2_rhor = 1. - wgt1_rhor
             irho_found = i_rhor
            if( wgt1_rhor > 1.0 .or. wgt1_rhor < 0. ) stop 'wrong wgt1_rhor'
             exit rhor_loop
         endif
      enddo rhor_loop
   endif

   !
   ! fraction weight
   !
   frime = min( max(frime_in, 0.), 1.)
   if( frime == 0. ) then
      wgt1_fr = 1. ; wgt2_fr = 0. ; ifr_found = 1
   else
      fr_loop: do i_Fr = 1, n_Fr-1
         wgt1_fr = 1. ; wgt2_fr = 0. ; ifr_found = 1
         if( frime > Fr( i_Fr ) .and. frime <= Fr( i_Fr+1 ) ) then
             wgt1_fr = ( Fr( i_Fr+1 ) - frime ) / ( Fr( i_Fr+1 ) - Fr( i_Fr ) )
             wgt2_fr = 1. - wgt1_fr
             ifr_found = i_Fr
             exit fr_loop
            if( wgt1_fr > 1.0 .or. wgt1_fr < 0. ) stop 'wrong wgt1_fr'
         endif
      enddo fr_loop
   endif

  call interpolate_polarris_lut_p3( num_ice, i_tair_ice, i_elev, iqn_found, irho_found, ifr_found, &
                                    wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                                    wgt1_qn, wgt2_qn, wgt1_rhor, wgt2_rhor, wgt1_fr, wgt2_fr, &
                                    bmueller_gen%ice, fscat_gen%ice, vdop_gen%ice )

 endif ice_present

 return
 end subroutine ss_gen_lut_p3

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine interpolate_polarris_lut( qm, i_re, ispc, i_tair, i_elev, &
                        wgt1_re, wgt2_re, wgt1_tair, wgt2_tair, wgt1_elev, wgt2_elev, &
                        bmueller0d, fscat0d, vdop0d )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the mueller scattering matrix for
! given elevation angle, temperature , and effective radius using POLARRIS LUTs (per unit q), 
! and values are re-scaled to actual qm (hydrometeor mass mixing ratio).
! 
! History:
!   02/2017  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: qm  !mass mixing ratio [g/m3]
 integer,intent(in) :: i_re
 integer,intent(in) :: ispc
 integer,intent(in) :: i_tair
 integer,intent(in) :: i_elev
 real(sdsu_fps),intent(in) :: wgt1_re, wgt2_re  !weights for effective radius (ice)
 real(sdsu_fps),intent(in) :: wgt1_tair, wgt2_tair  !weights for air temperature (liquid)
 real(sdsu_fps),intent(in) :: wgt1_elev, wgt2_elev  !weights for elevation
 
 real(sdsu_fpd),intent(out) :: & 
     bmueller0d(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop0d              ! Doppler velocity [m/s]

 complex,intent(out) :: & ! particle_gen is defined in module_simulater
     fscat0d(2,2)             ! 2x2 forward scattering matrix 

 real(sdsu_fps) :: fscat_real(2,2), fscat_imag(2,2)
 real(sdsu_fps) :: wgt_tot

!
! total weight
!
 wgt_tot =   wgt1_re * wgt1_elev * wgt1_tair &
           + wgt2_re * wgt1_elev * wgt1_tair &
           + wgt1_re * wgt2_elev * wgt1_tair &
           + wgt2_re * wgt2_elev * wgt1_tair &
           + wgt1_re * wgt1_elev * wgt2_tair &
           + wgt2_re * wgt1_elev * wgt2_tair &
           + wgt1_re * wgt2_elev * wgt2_tair &
           + wgt2_re * wgt2_elev * wgt2_tair

! if( wgt_tot /= 1.0 ) print*,'MSG interpolate_polarris_lut: wgt_tot is not ONE =',wgt_tot
!
! 4x4 mueller matrix, need qm weight
!
  bmueller0d(:,:) = qm * (  &
  ((   bmueller_array  (:,:,i_re-1,ispc,i_tair-1,i_elev-1) * wgt1_re &
     + bmueller_array  (:,:,i_re  ,ispc,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   bmueller_array  (:,:,i_re-1,ispc,i_tair-1,i_elev  ) * wgt1_re  &
     + bmueller_array  (:,:,i_re  ,ispc,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   bmueller_array  (:,:,i_re-1,ispc,i_tair,i_elev-1  ) * wgt1_re &
     + bmueller_array  (:,:,i_re  ,ispc,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   bmueller_array  (:,:,i_re-1,ispc,i_tair,i_elev    ) * wgt1_re  &
     + bmueller_array  (:,:,i_re  ,ispc,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    )

!
! 2x2 scattering matrix (real), need qm weight
!
  fscat_real(:,:) = qm * (  &
  ((   fscat_real_array(:,:,i_re-1,ispc,i_tair-1,i_elev-1) * wgt1_re &
     + fscat_real_array(:,:,i_re  ,ispc,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_real_array(:,:,i_re-1,ispc,i_tair-1,i_elev  ) * wgt1_re  &
     + fscat_real_array(:,:,i_re  ,ispc,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   fscat_real_array(:,:,i_re-1,ispc,i_tair,i_elev-1  ) * wgt1_re &
     + fscat_real_array(:,:,i_re  ,ispc,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_real_array(:,:,i_re-1,ispc,i_tair,i_elev    ) * wgt1_re  &
     + fscat_real_array(:,:,i_re  ,ispc,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    )

!
! 2x2 scattering matrix (imag), need qm weight
!
  fscat_imag(:,:) = qm * (  &
  ((   fscat_imag_array(:,:,i_re-1,ispc,i_tair-1,i_elev-1) * wgt1_re &
     + fscat_imag_array(:,:,i_re  ,ispc,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_imag_array(:,:,i_re-1,ispc,i_tair-1,i_elev  ) * wgt1_re  &
     + fscat_imag_array(:,:,i_re  ,ispc,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   fscat_imag_array(:,:,i_re-1,ispc,i_tair,i_elev-1  ) * wgt1_re &
     + fscat_imag_array(:,:,i_re  ,ispc,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_imag_array(:,:,i_re-1,ispc,i_tair,i_elev    ) * wgt1_re  &
     + fscat_imag_array(:,:,i_re  ,ispc,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    )

  fscat0d = cmplx( fscat_real, fscat_imag)

!
! doppler velocity (not qm weight)
!
  vdop0d = ( &
  ((   vdop_array(i_re-1,ispc,i_tair-1,i_elev-1) * wgt1_re &
     + vdop_array(i_re  ,ispc,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   vdop_array(i_re-1,ispc,i_tair-1,i_elev  ) * wgt1_re  &
     + vdop_array(i_re  ,ispc,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   vdop_array(i_re-1,ispc,i_tair,i_elev-1  ) * wgt1_re &
     + vdop_array(i_re  ,ispc,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   vdop_array(i_re-1,ispc,i_tair,i_elev    ) * wgt1_re  &
     + vdop_array(i_re  ,ispc,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    )


 return
 end subroutine interpolate_polarris_lut

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine interpolate_polarris_lut_rho( spc, qm, rhoe, i_re, i_tair, i_elev, &
                        wgt1_re, wgt2_re, wgt1_tair, wgt2_tair, wgt1_elev, wgt2_elev, &
                        bmueller0d, fscat0d, vdop0d )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the mueller scattering matrix for
! given elevation angle, temperature , and effective radius using POLARRIS LUTs (per unit q), 
! and values are re-scaled to actual qm (hydrometeor mass mixing ratio).
! 
! History:
!   02/2017  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: spc ! species
 real(sdsu_fps),intent(in) :: qm   !mass mixing ratio [g/m3]
 real(sdsu_fps),intent(in) :: rhoe !effective density [kg/m3]
 integer,intent(in) :: i_re
 integer,intent(in) :: i_tair
 integer,intent(in) :: i_elev
 real(sdsu_fps),intent(in) :: wgt1_re, wgt2_re  !weights for effective radius (ice)
 real(sdsu_fps),intent(in) :: wgt1_tair, wgt2_tair  !weights for air temperature (liquid)
 real(sdsu_fps),intent(in) :: wgt1_elev, wgt2_elev  !weights for elevation
 
 real(sdsu_fpd),intent(out) :: & 
     bmueller0d(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop0d              ! Doppler velocity [m/s]

 complex,intent(out) :: & ! particle_gen is defined in module_simulater
     fscat0d(2,2)             ! 2x2 forward scattering matrix 

 real(sdsu_fps) :: fscat_real(2,2), fscat_imag(2,2)
 real(sdsu_fps) :: wgt_tot, wgt1_rho, wgt2_rho

 integer :: ispc1, ispc2

!
! total weight (debbuging purpose)
!
 wgt_tot =   wgt1_re * wgt1_elev * wgt1_tair &
           + wgt2_re * wgt1_elev * wgt1_tair &
           + wgt1_re * wgt2_elev * wgt1_tair &
           + wgt2_re * wgt2_elev * wgt1_tair &
           + wgt1_re * wgt1_elev * wgt2_tair &
           + wgt2_re * wgt1_elev * wgt2_tair &
           + wgt1_re * wgt2_elev * wgt2_tair &
           + wgt2_re * wgt2_elev * wgt2_tair

! if( wgt_tot /= 1.0 ) print*,'MSG interpolate_polarris_lut: wgt_tot is not ONE =',wgt_tot

!
! density weight and index
!

if(trim(spc) == 'graupel') then

  if    ( rhoe >= 0.17e+3 .and. rhoe < 0.3e+3 ) then
    wgt2_rho = (rhoe - 0.17e+3 ) / ( 0.3e+3 - 0.17e+3 ) 
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 5 ; ispc2 = ispc1+1
  elseif( rhoe >= 0.3e+3 .and. rhoe < 0.4e+3 ) then
    wgt2_rho = (rhoe - 0.3e+3 ) / ( 0.4e+3 - 0.3e+3 )       
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 6 ; ispc2 = ispc1+1
  elseif( rhoe >= 0.4e+3 .and. rhoe < 0.5e+3 ) then
    wgt2_rho = (rhoe - 0.4e+3 ) / ( 0.5e+3 - 0.4e+3 )       
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 7 ; ispc2 = ispc1+1
  elseif( rhoe >= 0.5e+3 .and. rhoe < 0.6e+3 ) then
    wgt2_rho = (rhoe - 0.5e+3 ) / ( 0.6e+3 - 0.5e+3 )       
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 8 ; ispc2 = ispc1+1
  elseif( rhoe >= 0.6e+3 .and. rhoe < 0.7e+3 ) then
    wgt2_rho = (rhoe - 0.6e+3 ) / ( 0.7e+3 - 0.6e+3 )       
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 9 ; ispc2 = ispc1+1
  elseif( rhoe >= 0.7e+3 .and. rhoe < 0.8e+3 ) then
    wgt2_rho = (rhoe - 0.7e+3 ) / ( 0.8e+3 - 0.7e+3 )       
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 10 ; ispc2 = ispc1+1
  elseif( rhoe >= 0.8e+3 .and. rhoe <=0.9e+3 ) then
    wgt2_rho = (rhoe - 0.8e+3 ) / ( 0.9e+3 - 0.8e+3 )       
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 11 ; ispc2 = ispc1+1
  endif

elseif( trim(spc) == 'hail') then

  if    ( rhoe >= 0.5e+3 .and. rhoe < 0.6e+3 ) then
    wgt2_rho = (rhoe - 0.5e+3 ) / ( 0.6e+3 - 0.5e+3 )
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 13 ; ispc2 = ispc1+1
  elseif( rhoe >= 0.6e+3 .and. rhoe < 0.7e+3 ) then
    wgt2_rho = (rhoe - 0.6e+3 ) / ( 0.7e+3 - 0.6e+3 )
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 14 ; ispc2 = ispc1+1
  elseif( rhoe >= 0.7e+3 .and. rhoe < 0.8e+3 ) then
    wgt2_rho = (rhoe - 0.7e+3 ) / ( 0.8e+3 - 0.7e+3 )
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 15 ; ispc2 = ispc1+1
  elseif( rhoe >= 0.8e+3 .and. rhoe <=0.9e+3 ) then
    wgt2_rho = (rhoe - 0.8e+3 ) / ( 0.9e+3 - 0.8e+3 )
    wgt1_rho = 1. - wgt2_rho
    ispc1 = 16 ; ispc2 = ispc1+1
  endif

else 
 stop 'MSG: interpolate_polarris_lut_rho: no such spc'
endif

if( wgt1_rho > 1.0 .or. wgt1_rho < 0. ) stop 'MSG strange wgt1_rho'


!
! 4x4 mueller matrix, need qm weight
!
  bmueller0d(:,:) = &
  qm * (  &
  ((   bmueller_array  (:,:,i_re-1,ispc1,i_tair-1,i_elev-1) * wgt1_re &
     + bmueller_array  (:,:,i_re  ,ispc1,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   bmueller_array  (:,:,i_re-1,ispc1,i_tair-1,i_elev  ) * wgt1_re  &
     + bmueller_array  (:,:,i_re  ,ispc1,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   bmueller_array  (:,:,i_re-1,ispc1,i_tair,i_elev-1  ) * wgt1_re &
     + bmueller_array  (:,:,i_re  ,ispc1,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   bmueller_array  (:,:,i_re-1,ispc1,i_tair,i_elev    ) * wgt1_re  &
     + bmueller_array  (:,:,i_re  ,ispc1,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    ) * wgt1_rho + &
  qm * (  &
  ((   bmueller_array  (:,:,i_re-1,ispc2,i_tair-1,i_elev-1) * wgt1_re &
     + bmueller_array  (:,:,i_re  ,ispc2,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   bmueller_array  (:,:,i_re-1,ispc2,i_tair-1,i_elev  ) * wgt1_re  &
     + bmueller_array  (:,:,i_re  ,ispc2,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   bmueller_array  (:,:,i_re-1,ispc2,i_tair,i_elev-1  ) * wgt1_re &
     + bmueller_array  (:,:,i_re  ,ispc2,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   bmueller_array  (:,:,i_re-1,ispc2,i_tair,i_elev    ) * wgt1_re  &
     + bmueller_array  (:,:,i_re  ,ispc2,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    ) * wgt2_rho


!
! 2x2 scattering matrix (real), need qm weight
!
  fscat_real(:,:) = &
  qm * (  &
  ((   fscat_real_array(:,:,i_re-1,ispc1,i_tair-1,i_elev-1) * wgt1_re &
     + fscat_real_array(:,:,i_re  ,ispc1,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_real_array(:,:,i_re-1,ispc1,i_tair-1,i_elev  ) * wgt1_re  &
     + fscat_real_array(:,:,i_re  ,ispc1,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   fscat_real_array(:,:,i_re-1,ispc1,i_tair,i_elev-1  ) * wgt1_re &
     + fscat_real_array(:,:,i_re  ,ispc1,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_real_array(:,:,i_re-1,ispc1,i_tair,i_elev    ) * wgt1_re  &
     + fscat_real_array(:,:,i_re  ,ispc1,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    ) * wgt1_rho + &
  qm * (  &
  ((   fscat_real_array(:,:,i_re-1,ispc2,i_tair-1,i_elev-1) * wgt1_re &
     + fscat_real_array(:,:,i_re  ,ispc2,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_real_array(:,:,i_re-1,ispc2,i_tair-1,i_elev  ) * wgt1_re  &
     + fscat_real_array(:,:,i_re  ,ispc2,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   fscat_real_array(:,:,i_re-1,ispc2,i_tair,i_elev-1  ) * wgt1_re &
     + fscat_real_array(:,:,i_re  ,ispc2,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_real_array(:,:,i_re-1,ispc2,i_tair,i_elev    ) * wgt1_re  &
     + fscat_real_array(:,:,i_re  ,ispc2,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    ) * wgt2_rho 


!
! 2x2 scattering matrix (imag), need qm weight
!
  fscat_imag(:,:) = &
  qm * (  &
  ((   fscat_imag_array(:,:,i_re-1,ispc1,i_tair-1,i_elev-1) * wgt1_re &
     + fscat_imag_array(:,:,i_re  ,ispc1,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_imag_array(:,:,i_re-1,ispc1,i_tair-1,i_elev  ) * wgt1_re  &
     + fscat_imag_array(:,:,i_re  ,ispc1,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   fscat_imag_array(:,:,i_re-1,ispc1,i_tair,i_elev-1  ) * wgt1_re &
     + fscat_imag_array(:,:,i_re  ,ispc1,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_imag_array(:,:,i_re-1,ispc1,i_tair,i_elev    ) * wgt1_re  &
     + fscat_imag_array(:,:,i_re  ,ispc1,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    ) * wgt1_rho + &
  qm * (  &
  ((   fscat_imag_array(:,:,i_re-1,ispc2,i_tair-1,i_elev-1) * wgt1_re &
     + fscat_imag_array(:,:,i_re  ,ispc2,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_imag_array(:,:,i_re-1,ispc2,i_tair-1,i_elev  ) * wgt1_re  &
     + fscat_imag_array(:,:,i_re  ,ispc2,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   fscat_imag_array(:,:,i_re-1,ispc2,i_tair,i_elev-1  ) * wgt1_re &
     + fscat_imag_array(:,:,i_re  ,ispc2,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   fscat_imag_array(:,:,i_re-1,ispc2,i_tair,i_elev    ) * wgt1_re  &
     + fscat_imag_array(:,:,i_re  ,ispc2,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    ) * wgt2_rho 


  fscat0d = cmplx( fscat_real, fscat_imag)

!
! doppler velocity (not qm weight)
!
  vdop0d = &
  ( &
  ((   vdop_array(i_re-1,ispc1,i_tair-1,i_elev-1) * wgt1_re &
     + vdop_array(i_re  ,ispc1,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   vdop_array(i_re-1,ispc1,i_tair-1,i_elev  ) * wgt1_re  &
     + vdop_array(i_re  ,ispc1,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   vdop_array(i_re-1,ispc1,i_tair,i_elev-1  ) * wgt1_re &
     + vdop_array(i_re  ,ispc1,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   vdop_array(i_re-1,ispc1,i_tair,i_elev    ) * wgt1_re  &
     + vdop_array(i_re  ,ispc1,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    ) * wgt1_rho + &
  ( &
  ((   vdop_array(i_re-1,ispc2,i_tair-1,i_elev-1) * wgt1_re &
     + vdop_array(i_re  ,ispc2,i_tair-1,i_elev-1) * wgt2_re ) * wgt1_elev  &
       + &
   (   vdop_array(i_re-1,ispc2,i_tair-1,i_elev  ) * wgt1_re  &
     + vdop_array(i_re  ,ispc2,i_tair-1,i_elev  ) * wgt2_re ) * wgt2_elev ) * wgt1_tair &
            + &
  ((   vdop_array(i_re-1,ispc2,i_tair,i_elev-1  ) * wgt1_re &
     + vdop_array(i_re  ,ispc2,i_tair,i_elev-1  ) * wgt2_re ) * wgt1_elev  &
       + &
   (   vdop_array(i_re-1,ispc2,i_tair,i_elev    ) * wgt1_re  &
     + vdop_array(i_re  ,ispc2,i_tair,i_elev    ) * wgt2_re ) * wgt2_elev ) * wgt2_tair &
    ) * wgt2_rho 

 return
 end subroutine interpolate_polarris_lut_rho

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine interpolate_polarris_lut_p3( num, i_tair, i_elev, i_qn, i_rho, i_fr, &
                                         wgt1_tair, wgt2_tair, wgt1_elev, wgt2_elev, &
                                         wgt1_qn, wgt2_qn, wgt1_rhor, wgt2_rhor, wgt1_fr, wgt2_fr, &
                                         bmueller0d, fscat0d, vdop0d )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the mueller scattering matrix for
! given elevation angle, temperature , and effective radius using POLARRIS LUTs (per unit q), 
! and values are re-scaled to actual qm (hydrometeor mass mixing ratio).
! 
! History:
!   07/2019  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: num  !number concentration [#/m3] 
 integer,intent(in) :: i_tair
 integer,intent(in) :: i_elev
 integer,intent(in) :: i_qn
 integer,intent(in) :: i_rho
 integer,intent(in) :: i_fr

 real(sdsu_fps),intent(in) :: wgt1_tair, wgt2_tair  !weights for air temperature (liquid)
 real(sdsu_fps),intent(in) :: wgt1_elev, wgt2_elev  !weights for elevation
 real(sdsu_fps),intent(in) :: wgt1_qn, wgt2_qn      !weights for normalized q 
 real(sdsu_fps),intent(in) :: wgt1_rhor, wgt2_rhor  !weights for rime density
 real(sdsu_fps),intent(in) :: wgt1_fr, wgt2_fr      !weights for rime fraction
 
 real(sdsu_fpd),intent(out) :: & 
     bmueller0d(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop0d              ! Doppler velocity [m/s]

 complex,intent(out) :: & ! particle_gen is defined in module_simulater
     fscat0d(2,2)             ! 2x2 forward scattering matrix 

 real(sdsu_fps) :: fscat_real(2,2), fscat_imag(2,2)
 real(sdsu_fps) :: wgt_tot


! if( wgt_tot /= 1.0 ) print*,'MSG interpolate_polarris_lut_p3: wgt_tot is not ONE =',wgt_tot
!
! 4x4 mueller matrix, need qm weight
!
  bmueller0d(:,:) = num * (  &
((((   bmueller_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair-1, i_elev-1 ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + & 
   (   bmueller_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair-1, i_elev-1 ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   bmueller_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair-1, i_elev-1 ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   bmueller_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair-1, i_elev-1 ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt1_tair &
   + & 
 (((   bmueller_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair, i_elev-1 ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   bmueller_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair, i_elev-1 ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   bmueller_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair, i_elev-1 ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   bmueller_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair, i_elev-1 ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt2_tair ) * wgt1_elev &
  + &
((((   bmueller_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair-1, i_elev ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair-1, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   bmueller_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair-1, i_elev ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair-1, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   bmueller_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair-1, i_elev ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair-1, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   bmueller_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair-1, i_elev ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair-1, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt1_tair &
   + &
 (((   bmueller_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair, i_elev ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   bmueller_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair, i_elev ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   bmueller_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair, i_elev ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   bmueller_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair, i_elev ) * wgt1_rhor &
     + bmueller_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt2_tair ) * wgt2_elev &
    )



!
! 2x2 scattering matrix (real), need qm weight
!
  fscat_real(:,:) = num * (  &
((((   fscat_real_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair-1, i_elev-1 ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_real_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair-1, i_elev-1 ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   fscat_real_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair-1, i_elev-1 ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_real_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair-1, i_elev-1 ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt1_tair &
   + &
 (((   fscat_real_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair, i_elev-1 ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_real_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair, i_elev-1 ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   fscat_real_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair, i_elev-1 ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_real_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair, i_elev-1 ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt2_tair ) * wgt1_elev &
  + &
((((   fscat_real_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair-1, i_elev ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair-1, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_real_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair-1, i_elev ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair-1, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   fscat_real_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair-1, i_elev ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair-1, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_real_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair-1, i_elev ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair-1, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt1_tair &
   + &
 (((   fscat_real_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair, i_elev ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_real_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair, i_elev ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   fscat_real_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair, i_elev ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_real_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair, i_elev ) * wgt1_rhor &
     + fscat_real_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt2_tair ) * wgt2_elev &
    )


!
! 2x2 scattering matrix (imag), need qm weight
!
  fscat_imag(:,:) = num * (  &
((((   fscat_imag_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair-1, i_elev-1 ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_imag_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair-1, i_elev-1 ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   fscat_imag_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair-1, i_elev-1 ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_imag_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair-1, i_elev-1 ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt1_tair &
   + &
 (((   fscat_imag_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair, i_elev-1 ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_imag_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair, i_elev-1 ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   fscat_imag_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair, i_elev-1 ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_imag_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair, i_elev-1 ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt2_tair ) * wgt1_elev &
  + &
((((   fscat_imag_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair-1, i_elev ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair-1, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_imag_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair-1, i_elev ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair-1, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   fscat_imag_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair-1, i_elev ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair-1, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_imag_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair-1, i_elev ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair-1, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt1_tair &
   + &
 (((   fscat_imag_array_p3  ( :,:,i_rho  , i_fr  , i_qn  , i_tair, i_elev ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr  , i_qn  , i_tair, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_imag_array_p3  ( :,:,i_rho  , i_fr+1, i_qn  , i_tair, i_elev ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn  , i_tair, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   fscat_imag_array_p3  ( :,:,i_rho  , i_fr  , i_qn+1, i_tair, i_elev ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr  , i_qn+1, i_tair, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   fscat_imag_array_p3  ( :,:,i_rho  , i_fr+1, i_qn+1, i_tair, i_elev ) * wgt1_rhor &
     + fscat_imag_array_p3  ( :,:,i_rho+1, i_fr+1, i_qn+1, i_tair, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt2_tair ) * wgt2_elev &
    )

  ! merge real and imagenary component
  fscat0d = cmplx( fscat_real, fscat_imag)

!
! doppler velocity (not qm weight)
!
  vdop0d = ( &
((((   vdop_array_p3  ( i_rho  , i_fr  , i_qn  , i_tair-1, i_elev-1 ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr  , i_qn  , i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   vdop_array_p3  ( i_rho  , i_fr+1, i_qn  , i_tair-1, i_elev-1 ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr+1, i_qn  , i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   vdop_array_p3  ( i_rho  , i_fr  , i_qn+1, i_tair-1, i_elev-1 ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr  , i_qn+1, i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   vdop_array_p3  ( i_rho  , i_fr+1, i_qn+1, i_tair-1, i_elev-1 ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr+1, i_qn+1, i_tair-1, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt1_tair &
   + &
 (((   vdop_array_p3  ( i_rho  , i_fr  , i_qn  , i_tair, i_elev-1 ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr  , i_qn  , i_tair, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   vdop_array_p3  ( i_rho  , i_fr+1, i_qn  , i_tair, i_elev-1 ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr+1, i_qn  , i_tair, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   vdop_array_p3  ( i_rho  , i_fr  , i_qn+1, i_tair, i_elev-1 ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr  , i_qn+1, i_tair, i_elev-1 ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   vdop_array_p3  ( i_rho  , i_fr+1, i_qn+1, i_tair, i_elev-1 ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr+1, i_qn+1, i_tair, i_elev-1 ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt2_tair ) * wgt1_elev &
  + &
((((   vdop_array_p3  ( i_rho  , i_fr  , i_qn  , i_tair-1, i_elev ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr  , i_qn  , i_tair-1, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   vdop_array_p3  ( i_rho  , i_fr+1, i_qn  , i_tair-1, i_elev ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr+1, i_qn  , i_tair-1, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   vdop_array_p3  ( i_rho  , i_fr  , i_qn+1, i_tair-1, i_elev ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr  , i_qn+1, i_tair-1, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   vdop_array_p3  ( i_rho  , i_fr+1, i_qn+1, i_tair-1, i_elev ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr+1, i_qn+1, i_tair-1, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt1_tair &
   + &
 (((   vdop_array_p3  ( i_rho  , i_fr  , i_qn  , i_tair, i_elev ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr  , i_qn  , i_tair, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   vdop_array_p3  ( i_rho  , i_fr+1, i_qn  , i_tair, i_elev ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr+1, i_qn  , i_tair, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn  &
    + &
  ((   vdop_array_p3  ( i_rho  , i_fr  , i_qn+1, i_tair, i_elev ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr  , i_qn+1, i_tair, i_elev ) * wgt2_rhor ) * wgt1_fr  &
     + &
   (   vdop_array_p3  ( i_rho  , i_fr+1, i_qn+1, i_tair, i_elev ) * wgt1_rhor &
     + vdop_array_p3  ( i_rho+1, i_fr+1, i_qn+1, i_tair, i_elev ) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn ) * wgt2_tair ) * wgt2_elev &
    )

 return
 end subroutine interpolate_polarris_lut_p3

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine ss_sbm_lut( tair, elev, q, n, fmelt, frime, fmelt_g, fmelt_h, &
                        bmueller_sbm, fscat_sbm, vdop_sbm )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the mueller scattering matrix for
! given elevation angle and temperature for SBM bin using POLARRIS LUTs (per unit number density), 
! and values are re-scaled to actual qm (hydrometeor mass mixing ratio).
! 
! History:
!   04/2017  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in)::tair  ! layer temperature [K]
 real(sdsu_fps),intent(inout)::elev  ! elevation angle [deg]

 type ( particle_sbm ),intent(in) :: & ! particle_sbm is defined in module_simulater
     q           , & ! bulk mixing ratio [g/m3] 
     n   (1:nbin)    ! particle # concentration [1/m4]
 real(sdsu_fps),intent(in) :: fmelt(1:nbin)    ! mass melting fraction of snow [-]
 real(sdsu_fps),intent(in) :: frime(1:nbin)    ! mass rime fraction of snow [-]
 real(sdsu_fps),intent(in) :: fmelt_g(1:nbin)  ! mass melting fraction of graupel[-]
 real(sdsu_fps),intent(in) :: fmelt_h(1:nbin)  ! mass melting fraction of hail [-]

 type ( particle_sbm_r8 ),intent(out) :: & ! particle_sbm is defined in module_simulater
     bmueller_sbm(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_sbm              ! Doppler velocity [m/s]

 type ( particle_sbm_cmplx ),intent(out) :: & ! particle_sbm is defined in module_simulater
     fscat_sbm(2,2)             ! 2x2 forward scattering matrix 

 real(sdsu_fps) :: n_temp(1:nbin)  ! particle # concentration [1/m4]

 integer :: i_elev
 real(sdsu_fps) :: wgt1_elev, wgt2_elev  !weights for elevation

 integer :: it, i_tair_liq
 real(sdsu_fps) :: wgt1_tair_liq, wgt2_tair_liq  !weights for air temperature (liquid)

 integer :: i_tair_ice
 real(sdsu_fps) :: wgt1_tair_ice, wgt2_tair_ice  !weights for air temperature  (ice)

 integer :: i_re
 real(sdsu_fps) :: wgt1_re, wgt2_re  !weights for effective radius (ice)
 real(sdsu_fps) :: ref

 integer :: ispc
 real(sdsu_fps) :: qm

!
! interpolate elevation angle
! 
 if( elev < 0. ) then
    print*,'MSG ss_sbm_lut: elevatation angle is < 0. Correct --> 0.'
    elev =  0.
 elseif( elev > 90. ) then
    print*,'MSG ss_sbm_lut: elevatation angle is > 90. Correct --> 90.'
    elev = 90.
 endif

 i_elev = int(elev) + 2
 wgt1_elev = ( elev_polarris( i_elev ) - elev )
 wgt2_elev = 1. - wgt1_elev

!
! interpolate temperature (liq)
!
 if( tair < pts_temp_w_pol(1) ) then
   i_tair_liq = 2
   wgt1_tair_liq = 1. ; wgt2_tair_liq = 0.
 elseif( tair >= pts_temp_w_pol( mxpts_temp_pol ) ) then
   i_tair_liq = mxpts_temp_pol
   wgt1_tair_liq = 0. ; wgt2_tair_liq = 1.
 else
   do it = 2, mxpts_temp_pol
      if( pts_temp_w_pol(it-1) <= tair .and. tair < pts_temp_w_pol(it) ) then
          i_tair_liq = it
          wgt1_tair_liq = ( pts_temp_w_pol(it) - tair ) / &
                          ( pts_temp_w_pol(it) - pts_temp_w_pol(it-1) )
          wgt2_tair_liq = 1. - wgt1_tair_liq
      endif
   enddo
 endif

!
! interpolate temperature (ice)
!
 if( tair < pts_temp_i_pol(1) ) then
   i_tair_ice = 2
   wgt1_tair_ice = 1. ; wgt2_tair_ice = 0.
 elseif( tair >= pts_temp_i_pol( mxpts_temp_pol ) ) then
   i_tair_ice = mxpts_temp_pol
   wgt1_tair_ice = 0. ; wgt2_tair_ice = 1.
 else
   do it = 2, mxpts_temp_pol
      if( pts_temp_i_pol(it-1) <= tair .and. tair < pts_temp_i_pol(it) ) then
          i_tair_ice = it
          wgt1_tair_ice = ( pts_temp_i_pol(it) - tair ) / &
                          ( pts_temp_i_pol(it) - pts_temp_i_pol(it-1) )
          wgt2_tair_ice = 1. - wgt1_tair_ice
      endif
   enddo
 endif


!
! liq
!
 liq_present: if( q%liq == 0. ) then
    bmueller_sbm%liq = 0. ;  fscat_sbm%liq = 0. ;  vdop_sbm%liq = 0.
 else
   ispc = 1 ; n_temp(1:nbin) = n(1:nbin)%liq

   call interpolate_polarris_lut_sbm( n_temp, ispc, i_tair_liq, i_elev, &
                                      wgt1_tair_liq, wgt2_tair_liq, wgt1_elev, wgt2_elev, &
                                      bmueller_sbm%liq, fscat_sbm%liq, vdop_sbm%liq )

 endif liq_present


!
! ice column
!
 ice_col_present: if( q%ice_col == 0. ) then
    bmueller_sbm%ice_col = 0. ;  fscat_sbm%ice_col = 0. ;  vdop_sbm%ice_col = 0.
 else
   ispc = 2 ; n_temp(1:nbin) = n(1:nbin)%ice_col

   call interpolate_polarris_lut_sbm( n_temp, ispc, i_tair_ice, i_elev, &
                                      wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                                      bmueller_sbm%ice_col, fscat_sbm%ice_col, vdop_sbm%ice_col )

 endif ice_col_present


!
! ice plate
!
 ice_pla_present: if( q%ice_pla == 0. ) then
    bmueller_sbm%ice_pla = 0. ;  fscat_sbm%ice_pla = 0. ;  vdop_sbm%ice_pla = 0.
 else
   ispc = 3 ; n_temp(1:nbin) = n(1:nbin)%ice_pla

   call interpolate_polarris_lut_sbm( n_temp, ispc, i_tair_ice, i_elev, &
                                      wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                                      bmueller_sbm%ice_pla, fscat_sbm%ice_pla, vdop_sbm%ice_pla )

 endif ice_pla_present

!
! ice dendrite
!
 ice_den_present: if( q%ice_den == 0. ) then
    bmueller_sbm%ice_den = 0. ;  fscat_sbm%ice_den = 0. ;  vdop_sbm%ice_den = 0.
 else
   ispc = 4 ; n_temp(1:nbin) = n(1:nbin)%ice_den

   call interpolate_polarris_lut_sbm( n_temp, ispc, i_tair_ice, i_elev, &
                                      wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                                      bmueller_sbm%ice_den, fscat_sbm%ice_den, vdop_sbm%ice_den )

 endif ice_den_present


!
! graupel
!
 graupel_present: if( q%graupel == 0. ) then
    bmueller_sbm%graupel = 0. ;  fscat_sbm%graupel = 0. ;  vdop_sbm%graupel = 0.
 else
   ispc = 5 ; n_temp(1:nbin) = n(1:nbin)%graupel

   call interpolate_polarris_lut_sbm( n_temp, ispc, i_tair_ice, i_elev, &
                                      wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                                      bmueller_sbm%graupel, fscat_sbm%graupel, vdop_sbm%graupel )

 endif graupel_present


!
! hail
!
 hail_present: if( q%hail == 0. ) then
    bmueller_sbm%hail = 0. ;  fscat_sbm%hail = 0. ;  vdop_sbm%hail = 0.
 else
   ispc = 6 ; n_temp(1:nbin) = n(1:nbin)%hail

   call interpolate_polarris_lut_sbm( n_temp, ispc, i_tair_ice, i_elev, &
                                      wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                                      bmueller_sbm%hail, fscat_sbm%hail, vdop_sbm%hail )

 endif hail_present


!
! snow + riming 
!
 snow_present: if( q%snow == 0. ) then
    bmueller_sbm%snow = 0. ;  fscat_sbm%snow = 0. ;  vdop_sbm%snow = 0.
 else
   n_temp(1:nbin) = n(1:nbin)%snow

   call interpolate_polarris_lut_sbm_snow( n_temp, frime, i_tair_ice, i_elev, &
                                           wgt1_tair_ice, wgt2_tair_ice, wgt1_elev, wgt2_elev, &
                                           bmueller_sbm%snow, fscat_sbm%snow, vdop_sbm%snow )

 endif snow_present



 return
 end subroutine ss_sbm_lut

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine interpolate_polarris_lut_sbm( n_temp, ispc, i_tair, i_elev, &
                                          wgt1_tair, wgt2_tair, wgt1_elev, wgt2_elev, &
                                          bmueller0d, fscat0d, vdop0d )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the mueller scattering matrix for
! given elevation angle, temperature for each SBM bin using POLARRIS LUTs (per unit number density), 
! and values are re-scaled to actual qm (hydrometeor mass mixing ratio).
! 
! History:
!   04/2017  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: n_temp(1:nbin)  !number density [1/m4]
 integer,intent(in) :: ispc
 integer,intent(in) :: i_tair
 integer,intent(in) :: i_elev
 real(sdsu_fps),intent(in) :: wgt1_tair, wgt2_tair  !weights for air temperature (liquid)
 real(sdsu_fps),intent(in) :: wgt1_elev, wgt2_elev  !weights for elevation


 real(sdsu_fpd),intent(out) :: &
     bmueller0d(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop0d              ! Doppler velocity [m/s]

 complex,intent(out) :: & ! particle_gen is defined in module_simulater
     fscat0d(2,2)             ! 2x2 forward scattering matrix 

 real(sdsu_fps) :: fscat_real(2,2), fscat_imag(2,2)
 real(sdsu_fps) :: wgt_tot

 real(sdsu_fpd) :: bmueller_temp(4,4)
 real(sdsu_fpd) :: bscat      ! backscatter coeff
 real(sdsu_fpd) ::  vdop_sum 
 real(sdsu_fpd) :: wgt_sum  
 integer :: i


!
! total weight
!
 wgt_tot =   wgt1_elev * wgt1_tair &
           + wgt2_elev * wgt1_tair &
           + wgt1_elev * wgt2_tair &
           + wgt2_elev * wgt2_tair 


!
! initialize
!
 bmueller0d(:,:) = 0.  ; fscat_real(:,:) = 0. ; fscat_imag(:,:) = 0. 
 vdop_sum = 0. ; wgt_sum = 0.  

 do i = 1, nbin

!
! 4x4 mueller matrix, need qm weight
!
 
  bmueller_temp(:,:) = n_temp(i) * (  &
                    bmueller_array  (:,:,i,ispc,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
                  + bmueller_array  (:,:,i,ispc,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
                  + bmueller_array  (:,:,i,ispc,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
                  + bmueller_array  (:,:,i,ispc,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
                    )

  bmueller0d(:,:) = bmueller0d(:,:) + bmueller_temp(:,:)

!
! 2x2 scattering matrix (real), need qm weight
!
  fscat_real(:,:) = fscat_real(:,:) + n_temp(i) * (  &
                    fscat_real_array (:,:,i,ispc,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
                  + fscat_real_array (:,:,i,ispc,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
                  + fscat_real_array (:,:,i,ispc,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
                  + fscat_real_array (:,:,i,ispc,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
                    )

!
! 2x2 scattering matrix (imag), need qm weight
!
  fscat_imag(:,:) = fscat_imag(:,:) + n_temp(i) * (  &
                    fscat_imag_array (:,:,i,ispc,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
                  + fscat_imag_array (:,:,i,ispc,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
                  + fscat_imag_array (:,:,i,ispc,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
                  + fscat_imag_array (:,:,i,ispc,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
                    )

!
! doppler velocity (not qm weight)
!
  call bscat_mueller( bmueller_temp ( 1:4,1:4 ), bscat )

  vdop_sum = vdop_sum + bscat * ( &
             vdop_array (i,ispc,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
           + vdop_array (i,ispc,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
           + vdop_array (i,ispc,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
           + vdop_array (i,ispc,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
             )

  wgt_sum  = wgt_sum + bscat
 
 enddo

!
! output
!
 fscat0d = cmplx( fscat_real, fscat_imag)

 if( wgt_sum > 0. ) then
    vdop0d =  vdop_sum / wgt_sum
 else
    vdop0d = 0. 
 endif

 
 return
 end subroutine interpolate_polarris_lut_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine interpolate_polarris_lut_sbm_snow( n_temp, frime, i_tair, i_elev, &
                                               wgt1_tair, wgt2_tair, wgt1_elev, wgt2_elev, &
                                               bmueller0d, fscat0d, vdop0d )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the mueller scattering matrix for
! given elevation angle, temperature for each SBM bin using POLARRIS LUTs (per unit number density), 
! and values are re-scaled to actual qm (hydrometeor mass mixing ratio).
! 
! History:
!   04/2017  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: n_temp(1:nbin)  !number density [1/m4]
 real(sdsu_fps),intent(in) :: frime (1:nbin)  !riming fraction [-]

 integer,intent(in) :: i_tair
 integer,intent(in) :: i_elev
 real(sdsu_fps),intent(in) :: wgt1_tair, wgt2_tair  !weights for air temperature (liquid)
 real(sdsu_fps),intent(in) :: wgt1_elev, wgt2_elev  !weights for elevation


 real(sdsu_fpd),intent(out) :: &
     bmueller0d(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop0d              ! Doppler velocity [m/s]
          
 complex,intent(out) :: & ! particle_gen is defined in module_simulater
     fscat0d(2,2)             ! 2x2 forward scattering matrix 
         
 real(sdsu_fps) :: fscat_real(2,2), fscat_imag(2,2)
 real(sdsu_fps) :: wgt_tot
             
 real(sdsu_fpd) :: bmueller_temp(4,4)
 real(sdsu_fpd) :: bscat      ! backscatter coeff
 real(sdsu_fpd) ::  vdop_sum 
 real(sdsu_fpd) :: wgt_sum  
 integer :: i

 integer :: ispc
 real(sdsu_fps) :: wgt1_frime, wgt2_frime
 real(sdsu_fps) :: frime_temp

!   
! total weight 
!   
 wgt_tot =   wgt1_elev * wgt1_tair &
           + wgt2_elev * wgt1_tair &
           + wgt1_elev * wgt2_tair &
           + wgt2_elev * wgt2_tair


!
! initialize
!
 bmueller0d(:,:) = 0.  ; fscat_real(:,:) = 0. ; fscat_imag(:,:) = 0.
 vdop_sum = 0. ; wgt_sum = 0.

 do i = 1, nbin

!
! ispc interpolation weight
!
  frime_temp = max(0., min(0.9, frime(i)) )
  if( frime_temp == 0. ) then
    ispc = 8
    wgt2_frime = 0.
    wgt1_frime = 1. - wgt2_frime
  elseif( frime_temp == 0.9 ) then
    ispc = 16
    wgt2_frime = 1.
    wgt1_frime = 1. - wgt2_frime
  else
    ispc = int( frime_temp * 10. ) + 8
    wgt2_frime =  10. * ( frime_temp - real( int( frime_temp * 10. )    ) / 10. )
    wgt1_frime = 1. - wgt2_frime
  endif

! if( n_temp(i) > 0. ) &
!  print*, i, frime_temp, ispc-1, ispc, wgt1_frime, wgt2_frime, int( frime_temp * 10. ), frime_temp * 10.

!
! 4x4 mueller matrix, need qm weight
!

  bmueller_temp(:,:) = n_temp(i) * (  &
                   ( &
                    bmueller_array  (:,:,i,ispc-1,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
                  + bmueller_array  (:,:,i,ispc-1,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
                  + bmueller_array  (:,:,i,ispc-1,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
                  + bmueller_array  (:,:,i,ispc-1,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
                   ) * wgt1_frime &
                  + &
                   ( &
                    bmueller_array  (:,:,i,ispc,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
                  + bmueller_array  (:,:,i,ispc,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
                  + bmueller_array  (:,:,i,ispc,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
                  + bmueller_array  (:,:,i,ispc,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
                   ) * wgt2_frime &
                    )

  bmueller0d(:,:) = bmueller0d(:,:) + bmueller_temp(:,:)

!
! 2x2 scattering matrix (real), need qm weight
!
  fscat_real(:,:) = fscat_real(:,:) + n_temp(i) * (  &
                   ( &
                    fscat_real_array (:,:,i,ispc-1,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
                  + fscat_real_array (:,:,i,ispc-1,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
                  + fscat_real_array (:,:,i,ispc-1,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
                  + fscat_real_array (:,:,i,ispc-1,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
                   ) * wgt1_frime &
                  + &
                   ( &
                    fscat_real_array (:,:,i,ispc,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
                  + fscat_real_array (:,:,i,ispc,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
                  + fscat_real_array (:,:,i,ispc,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
                  + fscat_real_array (:,:,i,ispc,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
                   ) * wgt2_frime &
                    )

!
! 2x2 scattering matrix (imag), need qm weight
!
  fscat_imag(:,:) = fscat_imag(:,:) + n_temp(i) * (  &
                   ( &
                    fscat_imag_array (:,:,i,ispc-1,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
                  + fscat_imag_array (:,:,i,ispc-1,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
                  + fscat_imag_array (:,:,i,ispc-1,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
                  + fscat_imag_array (:,:,i,ispc-1,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
                   ) * wgt1_frime &
                  + &
                   ( &
                    fscat_imag_array (:,:,i,ispc,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
                  + fscat_imag_array (:,:,i,ispc,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
                  + fscat_imag_array (:,:,i,ispc,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
                  + fscat_imag_array (:,:,i,ispc,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
                    ) &
                    )

!
! doppler velocity (not qm weight)
!
  call bscat_mueller( bmueller_temp ( 1:4,1:4 ), bscat )

  vdop_sum = vdop_sum + bscat * ( &
            ( &
             vdop_array (i,ispc-1,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
           + vdop_array (i,ispc-1,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
           + vdop_array (i,ispc-1,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
           + vdop_array (i,ispc-1,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
            ) * wgt1_frime &
           + &
            ( &
             vdop_array (i,ispc,i_tair-1,i_elev-1) * wgt1_tair * wgt1_elev &
           + vdop_array (i,ispc,i_tair  ,i_elev-1) * wgt2_tair * wgt1_elev &
           + vdop_array (i,ispc,i_tair-1,i_elev  ) * wgt1_tair * wgt2_elev &
           + vdop_array (i,ispc,i_tair  ,i_elev  ) * wgt2_tair * wgt2_elev &
            ) * wgt2_frime &
             )

  wgt_sum  = wgt_sum + bscat

 enddo

!
! output
!
 fscat0d = cmplx( fscat_real, fscat_imag )

 if( wgt_sum > 0. ) then
    vdop0d =  vdop_sum / wgt_sum
 else
    vdop0d = 0.
 endif

 return
 end subroutine interpolate_polarris_lut_sbm_snow

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine prep_makelut_polarris
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! Prepare MPI parameters for constructing the look-up tables of Mueller and forward scattering matrix
! for polarimetric radar simulations.
!
! History:
! 01/2017  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 integer :: k, nspec, ntemp, nre  !loop indice
 integer :: max_spc
 integer :: cnt_loop
 integer :: nloop_max
 logical :: exit_loop
 integer :: mxpts_re_pol_temp

 integer,allocatable,dimension(:) :: &
      loop_nspec, loop_ntemp, loop_nre , &
      loop_rhor , loop_fr   , loop_qnorm

 integer            :: i_rhor          ! index for rho_rime loop                [1 .. n_rhor]
 integer            :: i_Fr            ! index for rime-mass-fraction loop      [1 .. n_Fr]
 integer            :: i_Qnorm         ! index for normalized qitot loop        [1 .. n_Qnorm]


!
! P3 scheme only use cloud and rain for makelut routine
! 

 P3_IF: if ( trim(cloud_microphysics) == 'P3' ) then  !toshiii debbug

    cnt_loop = 0
    do ntemp = 1, mxpts_temp_pol
        do i_rhor = 1, n_rhor
           do i_Fr = 1, n_Fr
               do i_Qnorm = 1, n_Qnorm
                  cnt_loop = 1 + cnt_loop
               enddo
           enddo
       enddo
    enddo

    !
    ! allocate idx array
    !
    nloop_max = cnt_loop
    if( allocated(loop_ntemp) ) deallocate(loop_ntemp)
    if( allocated(loop_rhor ) ) deallocate(loop_rhor )
    if( allocated(loop_fr   ) ) deallocate(loop_fr   )
    if( allocated(loop_qnorm) ) deallocate(loop_qnorm)
    

    allocate( loop_ntemp( nloop_max ), loop_rhor( nloop_max ), loop_fr( nloop_max ), loop_qnorm( nloop_max ) )

    !
    ! asign each parameter in giant loop
    !
    cnt_loop = 0
    do ntemp = 1, mxpts_temp_pol
        do i_rhor = 1, n_rhor
           do i_Fr = 1, n_Fr
               do i_Qnorm = 1, n_Qnorm
                  cnt_loop = 1 + cnt_loop
                  loop_ntemp(cnt_loop) = ntemp
                  loop_rhor (cnt_loop) = i_rhor
                  loop_fr   (cnt_loop) = i_Fr
                  loop_qnorm(cnt_loop) = i_Qnorm
               enddo
           enddo
        enddo
    enddo


   !
   ! special MPI subroutine call here
   !
   exit_loop = .false. !initialize as false

#if MPI > 0

    call mpi_sdsu_init  ! initialize MPI
    call mpi_sdsu_loop(nloop_max, exit_loop)  ! soring input files for each processor

#else

    ! Default values of MPI parameters for single CPU run 
    masterproc = .true.
    myrank = 0
    myl_start = 1  ; myl_end   = nloop_max
    print*, 'MSG: prep_makelut_polarris, warning you are about to generate POLARRIS LUT with single CPUs. '
    print*, 'It takes for a month.... You should use MPI=2 option with at least 100CPUs. '
    print*,'Best Scalable number of CPUs is ',nloop_max
!  stop    'Stop the program'

#endif

    !
    ! generating POLARRIS LUT
    !
    call makelut_polarris_p3( nloop_max, loop_ntemp, loop_rhor, loop_fr, loop_qnorm, exit_loop)

#if MPI > 0
  call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
#endif

 endif P3_IF


!
! choice of microphysics to specify microphysics constant
!
 select case(trim(type_microphysics))
   case('GEN')      ! general bulk microphsics 
     call bulk_DSD
     max_spc = 6 
     mxpts_re_pol_temp = mxpts_re_pol 
     
     ! 
     ! NSSL scheme has explicit density for graupel (+7 species) and hail (+4 species)
     !
     if ( trim(cloud_microphysics) == 'NSSL2MMT' ) max_spc = 6 + 7 + 4
  
     !
     ! P3 scheme only use cloud and rain for makelut routine
     ! 
     if ( trim(cloud_microphysics) == 'P3' ) max_spc = 2  
   
   case('RAMS')  !GCE RAMS 1- and 2-moment microphysics 
     max_spc = 16
     mxpts_re_pol_temp = mxpts_re_pol
   case('SBM')
     call read_SBM_bin
     max_spc = 16
     mxpts_re_pol_temp = nbin
   case default
     stop 'MSG prep_makelut_polarris: no such microphysics options'
 end select
 call re_rams_gamma('init')  !call rams for all cases.


 cnt_loop = 0

 do nspec = 1, max_spc
    do ntemp = 1, mxpts_temp_pol
       do nre  = 1, mxpts_re_pol_temp
          cnt_loop = 1 + cnt_loop
       enddo
    enddo
 enddo


!
! allocate idx array
!
 nloop_max = cnt_loop 
 if( allocated(loop_nspec) ) deallocate( loop_nspec )
 if( allocated(loop_ntemp) ) deallocate(loop_ntemp  )
 if( allocated(loop_nre  ) ) deallocate(loop_nre    )

 allocate( loop_nspec( nloop_max ), loop_ntemp( nloop_max ), loop_nre( nloop_max ) )

!
! asign each parameter in giant loop
!
 cnt_loop = 0

 do nspec = 1, max_spc
    do ntemp = 1, mxpts_temp_pol
       do nre  = 1, mxpts_re_pol_temp
          cnt_loop = 1 + cnt_loop
          loop_nspec(cnt_loop) = nspec
          loop_ntemp(cnt_loop) = ntemp
          loop_nre  (cnt_loop) = nre
       enddo
    enddo
 enddo

!
! special MPI subroutine call here
!
exit_loop = .false. !initialize as false

#if MPI > 0

 if ( trim(cloud_microphysics) /= 'P3' ) call mpi_sdsu_init  ! P3 scheme already initailized this...
 call mpi_sdsu_loop(nloop_max, exit_loop)  ! soring input files for each processor

#else

  ! Default values of MPI parameters for single CPU run 
  masterproc = .true.
  myrank = 0
  myl_start = 1  ; myl_end   = nloop_max
  print*, 'MSG: prep_makelut_polarris, warning you are about to generate POLARRIS LUT with single CPUs. '
  print*, 'It takes for a month.... You should use MPI=2 option with at least 100CPUs. '
  print*,'Best Scalable number of CPUs is ',nloop_max
!  stop    'Stop the program'

#endif

!
! generating POLARRIS LUT
!
 call makelut_polarris(nloop_max, loop_nspec, loop_ntemp, loop_nre, exit_loop, mxpts_re_pol_temp, max_spc)



#if MPI > 0
!
! terminate mpi process
!
 call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
 call mpi_sdsu_finale  !terminate mpi process

#endif

 stop

 return
 end subroutine prep_makelut_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine makelut_polarris(nloop_max, loop_nspec, loop_ntemp, loop_nre, exit_loop, mxpts_re_pol_temp, max_spc)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! Construct the look-up tables of Mueller and forward scattering matrix
! for polarimetric radar simulations. This is x8000 faster than original implementation. 
!
! History:
! 06/2017  Toshi Matsui@NASA GSFC ; Add polarimetric output in LUT.nc.
! 03/2017  Toshi Matsui@NASA GSFC ; Add SBM options.
! 01/2017  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: nloop_max !loop maxima
 integer,intent(in) :: loop_nspec(nloop_max)
 integer,intent(in) :: loop_ntemp(nloop_max)
 integer,intent(in) :: loop_nre  (nloop_max)
 logical,intent(in) :: exit_loop  !loop exit 
 integer,intent(in) :: mxpts_re_pol_temp !re index 
 integer,intent(in) :: max_spc  ! maximum particle species

!
! local parameters
!
 character(len=100) :: fname_out
 character(len=100):: char_freq
 real(sdsu_fps) :: rhoe       ! re-dependent bulk snow density [kg/m3]
 real(sdsu_fps) :: frequency
 real(sdsu_fps),parameter :: q_unit = 1.e0 ! unit water content [g/m3]
 integer :: nspec,i,j,l,m,nfreq, im, jm, nf
 logical :: fexist
 integer :: ispc   ! rams particle class
 real(sdsu_fps) :: cfmas, pwmas ! alpha_m & beta_m 
 real(sdsu_fps) :: cfvt, pwvt   ! alpha_m & beta_m 
 integer,parameter :: io = 2
 integer :: ierr, n

 real(sdsu_fps) :: frime_fake(1:nbin) , frime_temp(1:nbin)
 real(sdsu_fps) :: fmelt_fake(1:nbin)   ! fake 
 real(sdsu_fps) :: n_unit(1:nbin)       ! number density [1/m4]

! local
 real(sdsu_fpd) :: bmueller_tot(4,4)  !total mueller matrix
 complex :: fscat_tot(2,2)            !total forward scattering 
 type( polradar_parameter ) :: polout  !this is defined in module_mueller_matrix.F
 real(sdsu_fps) :: wave_m  !wavelength [m]
 real(sdsu_fps),allocatable,dimension(:,:,:,:) :: array_zhh, array_zdr, array_kdp

!
! LUT output
!
 real(sdsu_fpd) :: bmueller_mat(mxelev_polarris,4,4) ! 4x4 mueller matrix (backscatering component only)
 complex        :: fscat_mat   (mxelev_polarris,2,2) ! 2x2 forward scattering matrix 
 real(sdsu_fpd) :: vdop        (mxelev_polarris    )  ! Doppler velocity [m/s]

!
! NetCDF parameters
!
 integer :: ncid
 integer :: i_dimid, j_dimid, n_dimid, m_dimid, im_dimid, jm_dimid, is_dimid, js_dimid

 integer, parameter :: ndims6d = 6, ndims4d = 4
 integer :: dimids6d_mueller(ndims6d), dimids6d_scat(ndims6d), dimids4d_vdop(ndims4d)

 integer :: id, id_bmueller, id_fscat_real, id_fscat_imag, id_vdop, id_zhh, id_zdr, id_kdp, id_rhohv

 integer :: nens, max_ensemble
 character(len=1) :: nchar


!
! program start
!
 if(masterproc) print*, '- MSG makelut_polarris: LUT for POLARRIS '

!
! allocate giant LUT array and initialize parameters
!
 if( allocated(bmueller_array) ) then
       deallocate(bmueller_array  , &
                  fscat_real_array, &
                  fscat_imag_array, &
                  vdop_array      )
 endif
 allocate( bmueller_array  (4,4,mxpts_re_pol_temp,max_spc,mxpts_temp_pol,mxelev_polarris), &
           fscat_real_array(2,2,mxpts_re_pol_temp,max_spc,mxpts_temp_pol,mxelev_polarris), &
           fscat_imag_array(2,2,mxpts_re_pol_temp,max_spc,mxpts_temp_pol,mxelev_polarris), &
           vdop_array      (    mxpts_re_pol_temp,max_spc,mxpts_temp_pol,mxelev_polarris), &
           stat=ierr )
 if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_polarris: allocation error bunded array. ')
 bmueller_array   = undefined_dble
 fscat_real_array = undefined
 fscat_imag_array = undefined
 vdop_array       = undefined_dble

 if( masterproc ) then
   allocate( array_zhh(mxpts_re_pol_temp, max_spc,mxpts_temp_pol,mxelev_polarris), &
             array_zdr(mxpts_re_pol_temp, max_spc,mxpts_temp_pol,mxelev_polarris), &
             array_kdp(mxpts_re_pol_temp, max_spc,mxpts_temp_pol,mxelev_polarris), &
            stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_polarris: allocation error polout array. ')
  array_zhh = undefined
  array_zdr = undefined
  array_kdp = undefined
 endif

!
! dolan 2013 ensemble assumptions of axis ratio and orienting angle
!
if( polarris_axisorient_name == 'dolan2013' ) then
  max_ensemble = 4
else
  max_ensemble = 1
endif


!
! ensemble loop
!
 ensemble_loop: do nens = 1, max_ensemble

  if( polarris_axisorient_typ == 4 ) then
      write(nchar,"(I1.1)") nens ! 1 digit decimal 
      polarris_axisorient_name = 'dolan2013_pt'//nchar
  endif

!
! initialize POLARIS parameters
!
 call read_polarris_axisorient

!
! initialize 
!
 frime_fake(1:nbin) = 0.
 fmelt_fake(1:nbin) = 0.  !fake 
 n_unit(1:nbin) = 1.0  ! number density [1/m4]

!
! frequency loop
!
 freq_loop : do nfreq = 1 , mxinst_polarris

   !
   ! if same frequencies in freq loop, cycle loop
   !
   if(nfreq >= 2) then
     do nf = 1, nfreq-1
        if( freq_polarris(nfreq) == freq_polarris(nf) ) cycle freq_loop
     enddo
   endif

   frequency = freq_polarris(nfreq)
   if(masterproc) write(*,*) 'Freq = ',freq_polarris(nfreq),' GHz'
   if(masterproc) write(*,*) 'Wavelength = ',0.1*const_ghz_mm/freq_polarris(nfreq),' cm'

   wave_m = 0.1*const_ghz_mm/freq_polarris(nfreq)/100.  !wavelength [m]

   char_freq = nch_polarris(nfreq)

 if(masterproc) then 

!
! open netcdf file
!
 fname_out = 'POLARRIS_'//trim(cloud_microphysics)//'-'//trim(char_freq)//'_'//trim(polarris_axisorient_name)//'.nc'
 sdsu_io_file = trim(sdsu_dir_sslut)//trim(fname_out)

 print*,'MGS makelut_polarris: start creating POLARRIS LUT',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )               !32bit output

!
! Define netcdf dimensions
!
 call check( nf90_def_dim(ncid, 'nspec' , max_spc           , n_dimid) )
 call check( nf90_def_dim(ncid, 'itemp' , mxpts_temp_pol    , i_dimid) )
 call check( nf90_def_dim(ncid, 'jre'   , mxpts_re_pol_temp , j_dimid) )
 call check( nf90_def_dim(ncid, 'melev' , mxelev_polarris, m_dimid) )
 call check( nf90_def_dim(ncid, 'imueller' , 4, im_dimid) )
 call check( nf90_def_dim(ncid, 'jmueller' , 4, jm_dimid) )
 call check( nf90_def_dim(ncid, 'iscat' , 2, is_dimid) )
 call check( nf90_def_dim(ncid, 'jscat' , 2, js_dimid) )

 dimids6d_mueller = (/ im_dimid, jm_dimid, j_dimid, n_dimid, i_dimid, m_dimid /)
 dimids6d_scat    = (/ is_dimid, js_dimid, j_dimid, n_dimid, i_dimid, m_dimid /)
 dimids4d_vdop    = (/                     j_dimid, n_dimid, i_dimid, m_dimid /)



!
! global attribute (need to add case statement)
!

 mic_select: select case(trim(type_microphysics))
 case('GEN')      ! general bulk microphsics 

!
! cloud axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_cloud', &
                        ityp_axis_gen(1:num_axis_gen%cloud)%cloud  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_cloud', &
                        axis_ratio_gen(1:num_axis_gen%cloud)%cloud  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_cloud', ityp_dist_gen%cloud  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_cloud', mean_orient_gen%cloud  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_cloud', &
                       stdev_orient_gen(1:num_stdev_gen%cloud)%cloud  ) )


!
! rain axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_rain', &
                        ityp_axis_gen(1:num_axis_gen%rain)%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_rain', &
                        axis_ratio_gen(1:num_axis_gen%rain)%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_rain', ityp_dist_gen%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_rain', mean_orient_gen%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_rain', &
                       stdev_orient_gen(1:num_stdev_gen%rain)%rain  ) )

!
! ice axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_ice', &
                        ityp_axis_gen(1:num_axis_gen%ice)%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_ice', &
                        axis_ratio_gen(1:num_axis_gen%ice)%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_ice', ityp_dist_gen%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_ice', mean_orient_gen%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice', &
                       stdev_orient_gen(1:num_stdev_gen%ice)%ice  ) )

!
! snow axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_snow', &
                        ityp_axis_gen(1:num_axis_gen%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_snow', &
                        axis_ratio_gen(1:num_axis_gen%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_snow', ityp_dist_gen%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_snow', mean_orient_gen%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_gen_snow', &
                       stdev_orient_gen(1:num_stdev_gen%snow)%snow  ) )

!
! graupel axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_graupel', &
                        ityp_axis_gen(1:num_axis_gen%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_graupel', &
                        axis_ratio_gen(1:num_axis_gen%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_graupel', ityp_dist_gen%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_graupel', mean_orient_gen%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_graupel', &
                       stdev_orient_gen(1:num_stdev_gen%graupel)%graupel  ) )

!
! hail axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_hail', &
                        ityp_axis_gen(1:num_axis_gen%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_hail', &
                        axis_ratio_gen(1:num_axis_gen%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_hail', ityp_dist_gen%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_hail', mean_orient_gen%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_hail', &
                       stdev_orient_gen(1:num_stdev_gen%hail)%hail  ) )


 case('RAMS')  !RAMS type

!
! cloud axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_cloud1', &
                        ityp_axis_rams(1:num_axis_rams%cloud1)%cloud1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_cloud1', &
                        axis_ratio_rams(1:num_axis_rams%cloud1)%cloud1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_cloud1', ityp_dist_rams%cloud1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_cloud1', mean_orient_rams%cloud1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_cloud1', &
                       stdev_orient_rams(1:num_stdev_rams%cloud1)%cloud1  ) )

 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_cloud2', &
                        ityp_axis_rams(1:num_axis_rams%cloud2)%cloud2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_cloud2', &
                        axis_ratio_rams(1:num_axis_rams%cloud2)%cloud2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_cloud2', ityp_dist_rams%cloud2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_cloud2', mean_orient_rams%cloud2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_cloud2', &
                       stdev_orient_rams(1:num_stdev_rams%cloud2)%cloud2  ) )

!
! rain axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_rain', &
                        ityp_axis_rams(1:num_axis_rams%rain)%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_rain', &
                        axis_ratio_rams(1:num_axis_rams%rain)%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_rain', ityp_dist_rams%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_rain', mean_orient_rams%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_rain', &
                       stdev_orient_rams(1:num_stdev_rams%rain)%rain  ) )

!
! ice axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_ice1', &
                        ityp_axis_rams(1:num_axis_rams%ice1)%ice1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_ice1', &
                        axis_ratio_rams(1:num_axis_rams%ice1)%ice1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_ice1', ityp_dist_rams%ice1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_ice1', mean_orient_rams%ice1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice1', &
                       stdev_orient_rams(1:num_stdev_rams%ice1)%ice1  ) )

 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_ice2', &
                        ityp_axis_rams(1:num_axis_rams%ice2)%ice2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_ice2', &
                        axis_ratio_rams(1:num_axis_rams%ice2)%ice2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_ice2', ityp_dist_rams%ice2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_ice2', mean_orient_rams%ice2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice2', &
                       stdev_orient_rams(1:num_stdev_rams%ice2)%ice2  ) )

!
! snow axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_snow', &
                        ityp_axis_rams(1:num_axis_rams%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_snow', &
                        axis_ratio_rams(1:num_axis_rams%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_snow', ityp_dist_rams%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_snow', mean_orient_rams%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_rams_snow', &
                       stdev_orient_rams(1:num_stdev_rams%snow)%snow  ) )

!
! graupel axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_graupel', &
                        ityp_axis_rams(1:num_axis_rams%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_graupel', &
                        axis_ratio_rams(1:num_axis_rams%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_graupel', ityp_dist_rams%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_graupel', mean_orient_rams%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_graupel', &
                       stdev_orient_rams(1:num_stdev_rams%graupel)%graupel  ) )

!
! hail axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_hail', &
                        ityp_axis_rams(1:num_axis_rams%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_hail', &
                        axis_ratio_rams(1:num_axis_rams%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_hail', ityp_dist_rams%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_hail', mean_orient_rams%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_hail', &
                       stdev_orient_rams(1:num_stdev_rams%hail)%hail  ) )



 case('SBM') !HUCM SBM type

!
! liquid axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_liq', &
                        ityp_axis_sbm(1:num_axis_sbm%liq)%liq  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_liq', &
                        axis_ratio_sbm(1:num_axis_sbm%liq)%liq  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_liq', ityp_dist_sbm%liq  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_liq', mean_orient_sbm%liq  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_liq', &
                       stdev_orient_sbm(1:num_stdev_sbm%liq)%liq  ) )


!
! ice column axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_ice_col', &
                        ityp_axis_sbm(1:num_axis_sbm%ice_col)%ice_col  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_ice_col', &
                        axis_ratio_sbm(1:num_axis_sbm%ice_col)%ice_col  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_ice_col', ityp_dist_sbm%ice_col  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_ice_col', mean_orient_sbm%ice_col  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice_col', &
                       stdev_orient_sbm(1:num_stdev_sbm%ice_col)%ice_col  ) )


!
! ice plate axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_ice_pla', &
                        ityp_axis_sbm(1:num_axis_sbm%ice_pla)%ice_pla  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_ice_pla', &
                        axis_ratio_sbm(1:num_axis_sbm%ice_pla)%ice_pla  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_ice_pla', ityp_dist_sbm%ice_pla  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_ice_pla', mean_orient_sbm%ice_pla  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice_pla', &
                       stdev_orient_sbm(1:num_stdev_sbm%ice_pla)%ice_pla  ) )

!
! ice dendrite axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_ice_den', &
                        ityp_axis_sbm(1:num_axis_sbm%ice_den)%ice_den  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_ice_den', &
                        axis_ratio_sbm(1:num_axis_sbm%ice_den)%ice_den  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_ice_den', ityp_dist_sbm%ice_den  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_ice_den', mean_orient_sbm%ice_den  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice_den', &
                       stdev_orient_sbm(1:num_stdev_sbm%ice_den)%ice_den  ) )


!
! snow axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_snow', &
                        ityp_axis_sbm(1:num_axis_sbm%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_snow', &
                        axis_ratio_sbm(1:num_axis_sbm%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_snow', ityp_dist_sbm%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_snow', mean_orient_sbm%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_snow', &
                       stdev_orient_sbm(1:num_stdev_sbm%snow)%snow  ) )

!
! graupel axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_graupel', &
                        ityp_axis_sbm(1:num_axis_sbm%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_graupel', &
                        axis_ratio_sbm(1:num_axis_sbm%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_graupel', ityp_dist_sbm%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_graupel', mean_orient_sbm%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_graupel', &
                       stdev_orient_sbm(1:num_stdev_sbm%graupel)%graupel  ) )

!
! hail axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_hail', &
                        ityp_axis_sbm(1:num_axis_sbm%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_hail', &
                        axis_ratio_sbm(1:num_axis_sbm%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_hail', ityp_dist_sbm%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_hail', mean_orient_sbm%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_hail', &
                       stdev_orient_sbm(1:num_stdev_sbm%hail)%hail  ) )


 case default
   call stop_sdsu(myrank,'MSG makelut_polarris: There is no such cloud_microphysics')
 end select mic_select



!
! define netcdf output parameters
!
 call check( nf90_def_var(ncid, 'bmueller', NF90_DOUBLE, dimids6d_mueller, id) )
 call check( nf90_put_att(ncid, id, 'description',' 4x4 mueller matrix ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_bmueller = id

 call check( nf90_def_var(ncid, 'fscat_real', NF90_FLOAT, dimids6d_scat, id) )
 call check( nf90_put_att(ncid, id, 'description',' 2x2 scattering matrix (real)') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_fscat_real = id

 call check( nf90_def_var(ncid, 'fscat_imag', NF90_FLOAT, dimids6d_scat, id) )
 call check( nf90_put_att(ncid, id, 'description',' 2x2 scattering matrix (imagenary)') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_fscat_imag = id

 call check( nf90_def_var(ncid, 'vdop', NF90_DOUBLE, dimids4d_vdop, id) )
 call check( nf90_put_att(ncid, id, 'description',' doppler velocity ') )
 call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
 id_vdop = id


 call check( nf90_def_var(ncid, 'zhh', NF90_FLOAT, dimids4d_vdop, id) )
 call check( nf90_put_att(ncid, id, 'description',' reflectivity at H pol ') )
 call check( nf90_put_att(ncid, id, 'units', '[dBZ]') )
 id_zhh = id

 call check( nf90_def_var(ncid, 'zdr', NF90_FLOAT, dimids4d_vdop, id) )
 call check( nf90_put_att(ncid, id, 'description',' differential reflectivity ') )
 call check( nf90_put_att(ncid, id, 'units', '[dB]') )
 id_zdr = id

 call check( nf90_def_var(ncid, 'kdp', NF90_FLOAT, dimids4d_vdop, id) )
 call check( nf90_put_att(ncid, id, 'description',' specific diff phase ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg/km]') )
 id_kdp = id



 call check( nf90_enddef(ncid) )  ! End define mode.

 endif



#if MPI > 0
 call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
#endif


!
! giant loop
!
   giant_loop: do l = myl_start, myl_end   

     if( exit_loop ) cycle  !this ensure to cycle proces, when processor is not assigned for task in MPI=1 options

     !
     ! asign original nspec, i, j parameters from giant loop index
     !
     nspec = loop_nspec(l)  ! species index
     i     = loop_ntemp(l)  ! temperature index
     j     = loop_nre  (l)  ! effective radius index (or bin for SBM)

     ! zero out output parameters
     bmueller_mat = 0. ; fscat_mat = 0. ; vdop = 0. 

      mic_select2: select case(trim(type_microphysics))

      case('GEN')      ! general bulk microphsics 

        NSSL_IF: if ( trim(cloud_microphysics) == 'NSSL2MMT' ) then  !special class for explicit riming rate
          if(nspec==1) then

           call ss_gen_polarris('cloud',frequency, pts_temp_w_pol(i), mu_gen%cloud, q_unit, pts_re_small(j), rho_gen%cloud , 1.,& !cloud drop
                                mxelev_polarris, elev_polarris, &
                                num_axis_gen%cloud, ityp_axis_gen%cloud, axis_ratio_gen%cloud, ityp_dist_gen%cloud,  &
                                mean_orient_gen%cloud, num_stdev_gen%cloud, stdev_orient_gen%cloud, &
                                a_vt%cloud, b_vt%cloud, bmueller_mat, fscat_mat , vdop)

          elseif(nspec==2) then

          call ss_gen_polarris('rain',frequency, pts_temp_w_pol(i), mu_gen%rain, q_unit, pts_re_large(j), rho_gen%rain , 1.,& !rain
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%rain, ityp_axis_gen%rain, axis_ratio_gen%rain, ityp_dist_gen%rain,  &
                             mean_orient_gen%rain, num_stdev_gen%rain, stdev_orient_gen%rain, &
                             a_vt%rain, b_vt%rain, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==3) then

          call ss_gen_polarris('ice',frequency, pts_temp_i_pol(i), mu_gen%ice, q_unit, pts_re_small(j), rho_gen%ice , 0.,& !ice
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%ice, ityp_axis_gen%ice, axis_ratio_gen%ice, ityp_dist_gen%ice, &
                             mean_orient_gen%ice, num_stdev_gen%ice, stdev_orient_gen%ice, &
                             a_vt%ice, b_vt%ice, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==4) then

          call ss_gen_polarris('snow',frequency, pts_temp_i_pol(i), mu_gen%snow, q_unit, pts_re_large(j), rho_gen%snow, 0.,& !snow
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%snow, ityp_axis_gen%snow, axis_ratio_gen%snow, ityp_dist_gen%snow, &
                             mean_orient_gen%snow, num_stdev_gen%snow, stdev_orient_gen%snow, &
                             a_vt%snow, b_vt%snow, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==5) then

          rhoe = 0.17e+3  !graupel density [kg/m3]
          call ss_gen_polarris('graupel_017',frequency,pts_temp_i_pol(i),mu_gen%graupel,q_unit,pts_re_large(j),rhoe,0.,& !graupel
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel, &
                             mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                             a_vt%graupel, b_vt%graupel, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==6) then

          rhoe = 0.30e+3  !graupel density [kg/m3]
          call ss_gen_polarris('graupel_030',frequency,pts_temp_i_pol(i),mu_gen%graupel,q_unit,pts_re_large(j),rhoe,0.,& !graupel
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel, &
                             mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                             a_vt%graupel, b_vt%graupel, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==7) then
          
          rhoe = 0.40e+3  !graupel density [kg/m3]
          call ss_gen_polarris('graupel_040',frequency,pts_temp_i_pol(i),mu_gen%graupel,q_unit,pts_re_large(j),rhoe,0.,& !graupel
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel, &
                             mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                             a_vt%graupel, b_vt%graupel, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==8) then

          rhoe = 0.50e+3  !graupel density [kg/m3]
          call ss_gen_polarris('graupel_050',frequency,pts_temp_i_pol(i),mu_gen%graupel,q_unit,pts_re_large(j),rhoe,0.,& !graupel
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel, &
                             mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                             a_vt%graupel, b_vt%graupel, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==9) then

          rhoe = 0.60e+3  !graupel density [kg/m3]
          call ss_gen_polarris('graupel_060',frequency,pts_temp_i_pol(i),mu_gen%graupel,q_unit,pts_re_large(j),rhoe,0.,& !graupel
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel, &
                             mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                             a_vt%graupel, b_vt%graupel, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==10) then

          rhoe = 0.70e+3  !graupel density [kg/m3]
          call ss_gen_polarris('graupel_070',frequency,pts_temp_i_pol(i),mu_gen%graupel,q_unit,pts_re_large(j),rhoe,0.,& !graupel
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel, &
                             mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                             a_vt%graupel, b_vt%graupel, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==11) then

          rhoe = 0.80e+3  !graupel density [kg/m3]
          call ss_gen_polarris('graupel_080',frequency,pts_temp_i_pol(i),mu_gen%graupel,q_unit,pts_re_large(j),rhoe,0.,& !graupel
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel, &
                             mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                             a_vt%graupel, b_vt%graupel, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==12) then

          rhoe = 0.90e+3  !graupel density [kg/m3]
          call ss_gen_polarris('graupel_090',frequency,pts_temp_i_pol(i),mu_gen%graupel,q_unit,pts_re_large(j),rhoe,0.,& !graupel
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel, &
                             mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                             a_vt%graupel, b_vt%graupel, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==13) then

          rhoe = 0.50e+3  !hail density [kg/m3]
          call ss_gen_polarris('hail_050',frequency, pts_temp_i_pol(i), mu_gen%hail, q_unit, pts_re_large(j), rhoe , 0.,& !hail
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%hail, ityp_axis_gen%hail, axis_ratio_gen%hail, ityp_dist_gen%hail, &
                             mean_orient_gen%hail, num_stdev_gen%hail, stdev_orient_gen%hail, &
                             a_vt%hail, b_vt%hail, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==14) then

          rhoe = 0.60e+3  !hail density [kg/m3]
          call ss_gen_polarris('hail_060',frequency, pts_temp_i_pol(i), mu_gen%hail, q_unit, pts_re_large(j), rhoe , 0.,& !hail
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%hail, ityp_axis_gen%hail, axis_ratio_gen%hail, ityp_dist_gen%hail, &
                             mean_orient_gen%hail, num_stdev_gen%hail, stdev_orient_gen%hail, &
                             a_vt%hail, b_vt%hail, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==15) then

          rhoe = 0.70e+3  !hail density [kg/m3]
          call ss_gen_polarris('hail_070',frequency, pts_temp_i_pol(i), mu_gen%hail, q_unit, pts_re_large(j), rhoe , 0.,& !hail
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%hail, ityp_axis_gen%hail, axis_ratio_gen%hail, ityp_dist_gen%hail, &
                             mean_orient_gen%hail, num_stdev_gen%hail, stdev_orient_gen%hail, &
                             a_vt%hail, b_vt%hail, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==16) then

          rhoe = 0.80e+3  !hail density [kg/m3]
          call ss_gen_polarris('hail_080',frequency, pts_temp_i_pol(i), mu_gen%hail, q_unit, pts_re_large(j), rhoe , 0.,& !hail
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%hail, ityp_axis_gen%hail, axis_ratio_gen%hail, ityp_dist_gen%hail, &
                             mean_orient_gen%hail, num_stdev_gen%hail, stdev_orient_gen%hail, &
                             a_vt%hail, b_vt%hail, bmueller_mat, fscat_mat, vdop)

          elseif(nspec==17) then

          rhoe = 0.90e+3  !hail density [kg/m3]
          call ss_gen_polarris('hail_090',frequency, pts_temp_i_pol(i), mu_gen%hail, q_unit, pts_re_large(j), rhoe , 0.,& !hail
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%hail, ityp_axis_gen%hail, axis_ratio_gen%hail, ityp_dist_gen%hail, &
                             mean_orient_gen%hail, num_stdev_gen%hail, stdev_orient_gen%hail, &
                             a_vt%hail, b_vt%hail, bmueller_mat, fscat_mat, vdop)

           endif

        else  !NSSL_IF

        ! LUT as function of T and q (for empirical DSD)
        NSPEC_GEN: if(nspec==1) then

          if( trim(cloud_microphysics) == 'WDM') then 

             call ss_gen_polarris('cloud',frequency, pts_temp_w_pol(i), mu_gen%cloud, q_unit, pts_re_small(j), rho_gen%cloud , 1.,& !cloud drop
                                mxelev_polarris, elev_polarris, &
                                num_axis_gen%cloud, ityp_axis_gen%cloud, axis_ratio_gen%cloud, ityp_dist_gen%cloud,  &
                                mean_orient_gen%cloud, num_stdev_gen%cloud, stdev_orient_gen%cloud, & 
                                a_vt%cloud, b_vt%cloud, bmueller_mat, fscat_mat , vdop)

          elseif( trim(cloud_microphysics) == 'P3') then

             !
             ! toshii currently LUT specify cloud mu as 5 (2~15). 
             ! for weather radar, cloud mu does not affect radar variables anyway. 
             ! may need to change in future. 
             !
             call ss_gen_polarris('cloud',frequency, pts_temp_w_pol(i), 5., q_unit, pts_re_small(j), rho_gen%cloud , 1.,& !cloud drop
                                mxelev_polarris, elev_polarris, &
                                num_axis_gen%cloud, ityp_axis_gen%cloud, axis_ratio_gen%cloud, ityp_dist_gen%cloud,  &
                                mean_orient_gen%cloud, num_stdev_gen%cloud, stdev_orient_gen%cloud, &
                                a_vt%cloud, b_vt%cloud, bmueller_mat, fscat_mat , vdop)


          else

             ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud1 --> uisng RAMS psd
                        cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)    ! but using gen orienting angle
             call ss_rams_polarris('cloud1',frequency, pts_temp_w_pol(i), q_unit, pts_re_small(j) , gnu_rams%cloud1 , cfvt,pwvt,&
                                 mxelev_polarris, elev_polarris, &
                                 num_axis_gen%cloud, ityp_axis_gen%cloud, axis_ratio_gen%cloud, ityp_dist_gen%cloud,  &
                                 mean_orient_gen%cloud, num_stdev_gen%cloud, stdev_orient_gen%cloud, &
                                 cfmas, pwmas , 1.0,  bmueller_mat, fscat_mat, vdop )


          endif

        elseif(nspec==2) then

          call ss_gen_polarris('rain',frequency, pts_temp_w_pol(i), mu_gen%rain, q_unit, pts_re_large(j), rho_gen%rain , 1.,& !rain
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%rain, ityp_axis_gen%rain, axis_ratio_gen%rain, ityp_dist_gen%rain,  &
                             mean_orient_gen%rain, num_stdev_gen%rain, stdev_orient_gen%rain, &
                             a_vt%rain, b_vt%rain, bmueller_mat, fscat_mat, vdop)

        elseif(nspec==3) then

         if( trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH' .or.&
             trim(cloud_microphysics) == 'THOM' ) then !  twomoment
          call ss_gen_polarris('ice',frequency, pts_temp_i_pol(i), mu_gen%ice, q_unit, pts_re_small(j), rho_gen%ice , 0.,& !ice
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%ice, ityp_axis_gen%ice, axis_ratio_gen%ice, ityp_dist_gen%ice, &
                             mean_orient_gen%ice, num_stdev_gen%ice, stdev_orient_gen%ice, &
                             a_vt%ice, b_vt%ice, bmueller_mat, fscat_mat, vdop)

         else  
                ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice -> RAMS ice1 col
                          cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! but using gen orienting angle
          call ss_rams_polarris('ice1', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice1 , cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_gen%ice, ityp_axis_gen%ice, axis_ratio_gen%ice, ityp_dist_gen%ice, &
                              mean_orient_gen%ice, num_stdev_gen%ice, stdev_orient_gen%ice, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

         endif

        elseif(nspec==4) then

          if( trim(cloud_microphysics) == 'GMP4ICE' ) then
              rhoe = 1.e+3 * max(0.05, min(0.9e0, 142.36 * ( pts_re_large(j)**(-0.95) ) ) ) ! re-dependent bulk snow density [kg/m3]
          elseif( trim(cloud_microphysics) == 'THOM' ) then
              rhoe = min(917. , 0.13 / (pts_re_large(j)*2.*1.e-6) ) ! Kg/m^3  (917. to be the same as densice)
          else
              rhoe = rho_gen%snow
          endif

          call ss_gen_polarris('snow',frequency, pts_temp_i_pol(i), mu_gen%snow, q_unit, pts_re_large(j), rhoe , 0.,& !snow
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%snow, ityp_axis_gen%snow, axis_ratio_gen%snow, ityp_dist_gen%snow, &
                             mean_orient_gen%snow, num_stdev_gen%snow, stdev_orient_gen%snow, &
                             a_vt%snow, b_vt%snow, bmueller_mat, fscat_mat, vdop)


        elseif(nspec==5) then

          call ss_gen_polarris('graupel',frequency,pts_temp_i_pol(i),mu_gen%graupel,q_unit,pts_re_large(j),rho_gen%graupel,0.,& !graupel
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%graupel, ityp_axis_gen%graupel, axis_ratio_gen%graupel, ityp_dist_gen%graupel, &
                             mean_orient_gen%graupel, num_stdev_gen%graupel, stdev_orient_gen%graupel, &
                             a_vt%graupel, b_vt%graupel, bmueller_mat, fscat_mat, vdop)

        elseif(nspec==6) then

          call ss_gen_polarris('hail',frequency, pts_temp_i_pol(i), mu_gen%hail, q_unit, pts_re_large(j), rho_gen%hail , 0.,& !hail
                             mxelev_polarris, elev_polarris, &
                             num_axis_gen%hail, ityp_axis_gen%hail, axis_ratio_gen%hail, ityp_dist_gen%hail, &
                             mean_orient_gen%hail, num_stdev_gen%hail, stdev_orient_gen%hail, &
                             a_vt%hail, b_vt%hail, bmueller_mat, fscat_mat, vdop)

        else
          call stop_sdsu(myrank,'MSG makelut_polarris: there is no such nspec')
        endif NSPEC_GEN

        !special case for NSSL scheme for 

        endif NSSL_IF



      case('RAMS')  !GCE RAMS 1- & 2-moment microphysics 

        NSPEC_RAMS: if(nspec==1) then
          ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud1
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('cloud1',frequency, pts_temp_w_pol(i), q_unit, pts_re_small(j) , gnu_rams%cloud1 , cfvt,pwvt,&
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%cloud1, ityp_axis_rams%cloud1, axis_ratio_rams%cloud1, ityp_dist_rams%cloud1, &
                              mean_orient_rams%cloud1, num_stdev_rams%cloud1, stdev_orient_rams%cloud1, &
                              cfmas, pwmas , 1.0,  bmueller_mat, fscat_mat, vdop )
 
        elseif(nspec==2) then
          ispc = 16; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud2
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('cloud2', frequency, pts_temp_w_pol(i), q_unit, pts_re_small(j) , gnu_rams%cloud2 , cfvt,pwvt,&
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%cloud2, ityp_axis_rams%cloud2, axis_ratio_rams%cloud2, ityp_dist_rams%cloud2, &
                              mean_orient_rams%cloud2, num_stdev_rams%cloud2, stdev_orient_rams%cloud2, &
                              cfmas, pwmas , 1.0,  bmueller_mat, fscat_mat, vdop )

        elseif(nspec==3) then
          ispc = 2 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! rain 
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('rain', frequency, pts_temp_w_pol(i), q_unit, pts_re_large(j) , gnu_rams%rain, cfvt,pwvt,&
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%rain, ityp_axis_rams%rain, axis_ratio_rams%rain, ityp_dist_rams%rain, &
                              mean_orient_rams%rain, num_stdev_rams%rain, stdev_orient_rams%rain, &
                              cfmas, pwmas , 1.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==4) then
          ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 col
                    cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice1_col', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice1 , cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice1, ityp_axis_rams%ice1, axis_ratio_rams%ice1, ityp_dist_rams%ice1, &
                              mean_orient_rams%ice1, num_stdev_rams%ice1, stdev_orient_rams%ice1, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==5) then
          ispc = 8 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 hex
                    cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice1_hex', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice1 ,cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice1, ityp_axis_rams%ice1, axis_ratio_rams%ice1, ityp_dist_rams%ice1, &
                              mean_orient_rams%ice1, num_stdev_rams%ice1, stdev_orient_rams%ice1, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==6) then
          ispc = 9 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 den
                    cfvt  = rams_dstprms(4,ispc) ;pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice1_den', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice1 , cfvt,pwvt,&
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice1, ityp_axis_rams%ice1, axis_ratio_rams%ice1, ityp_dist_rams%ice1, &
                              mean_orient_rams%ice1, num_stdev_rams%ice1, stdev_orient_rams%ice1, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==7) then
          ispc = 10 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 ndl
                     cfvt  = rams_dstprms(4,ispc) ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice1_ndl', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice1 ,cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice1, ityp_axis_rams%ice1, axis_ratio_rams%ice1, ityp_dist_rams%ice1, &
                              mean_orient_rams%ice1, num_stdev_rams%ice1, stdev_orient_rams%ice1, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==8) then
          ispc = 11 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 ros
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice1_ros', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice1 ,cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice1, ityp_axis_rams%ice1, axis_ratio_rams%ice1, ityp_dist_rams%ice1, &
                              mean_orient_rams%ice1, num_stdev_rams%ice1, stdev_orient_rams%ice1, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==9) then
          ispc = 4 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 col
                    cfvt  = rams_dstprms(4,ispc) ;pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice2_col', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice2 , cfvt,pwvt,&
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice2, ityp_axis_rams%ice2, axis_ratio_rams%ice2, ityp_dist_rams%ice2, &
                              mean_orient_rams%ice2, num_stdev_rams%ice2, stdev_orient_rams%ice2, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==10) then
          ispc = 12 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 hex
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice2_hex', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice2 , cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice2, ityp_axis_rams%ice2, axis_ratio_rams%ice2, ityp_dist_rams%ice2, &
                              mean_orient_rams%ice2, num_stdev_rams%ice2, stdev_orient_rams%ice2, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==11) then
          ispc = 13 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 den
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice2_den', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice2 ,cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice2, ityp_axis_rams%ice2, axis_ratio_rams%ice2, ityp_dist_rams%ice2, &
                              mean_orient_rams%ice2, num_stdev_rams%ice2, stdev_orient_rams%ice2, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==12) then
          ispc = 14 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 ndl
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice2_ndl', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice2 ,cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice2, ityp_axis_rams%ice2, axis_ratio_rams%ice2, ityp_dist_rams%ice2, &
                              mean_orient_rams%ice2, num_stdev_rams%ice2, stdev_orient_rams%ice2, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==13) then
          ispc = 15 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 ros
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('ice2_ros', frequency, pts_temp_i_pol(i), q_unit, pts_re_small(j) , gnu_rams%ice2 ,cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%ice2, ityp_axis_rams%ice2, axis_ratio_rams%ice2, ityp_dist_rams%ice2, &
                              mean_orient_rams%ice2, num_stdev_rams%ice2, stdev_orient_rams%ice2, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==14) then 
          ispc = 5 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! snow aggregate
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('snow', frequency, pts_temp_i_pol(i), q_unit, pts_re_large(j) , gnu_rams%snow,cfvt,pwvt, &   
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%snow, ityp_axis_rams%snow, axis_ratio_rams%snow, ityp_dist_rams%snow, &
                              mean_orient_rams%snow, num_stdev_rams%snow, stdev_orient_rams%snow, &
                              cfmas, pwmas ,  0.0, bmueller_mat, fscat_mat, vdop )
 
        elseif(nspec==15) then 
          ispc = 6 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! graupel
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('graupel', frequency, pts_temp_i_pol(i), q_unit, pts_re_large(j) , gnu_rams%graupel,cfvt,pwvt, &   
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%graupel, ityp_axis_rams%graupel, axis_ratio_rams%graupel, ityp_dist_rams%graupel, &
                              mean_orient_rams%graupel, num_stdev_rams%graupel, stdev_orient_rams%graupel, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )

        elseif(nspec==16) then
          ispc = 7 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! hail
                     cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! 
          call ss_rams_polarris('hail', frequency, pts_temp_i_pol(i), q_unit, pts_re_large(j) , gnu_rams%hail , cfvt,pwvt, &
                              mxelev_polarris, elev_polarris, &
                              num_axis_rams%hail, ityp_axis_rams%hail, axis_ratio_rams%hail, ityp_dist_rams%hail, &
                              mean_orient_rams%hail, num_stdev_rams%hail, stdev_orient_rams%hail, &
                              cfmas, pwmas , 0.0, bmueller_mat, fscat_mat, vdop )
            
        else
          call stop_sdsu(myrank,'MSG makelut_polarris: there is no such nspec for RAMS microphysics')
        endif NSPEC_RAMS

      case('SBM')


        NSPEC_SBM: if(nspec==1) then

          call ss_sbm_polarris('liq    ',j,j,frequency,pts_temp_w_pol(i),&
                               den_sbm(1:nbin)%liq,rad_sbm(1:nbin)%liq,drad_sbm(1:nbin)%liq,&
                               vt_sbm(1:nbin)%liq, q_unit, n_unit(j:j), fmelt_fake(j:j), frime_fake(j:j), &
                               mxelev_polarris, elev_polarris, &
                               num_axis_sbm%liq, ityp_axis_sbm%liq, axis_ratio_sbm%liq, ityp_dist_sbm%liq,  &
                                mean_orient_sbm%liq, num_stdev_sbm%liq, stdev_orient_sbm%liq, &
                                bmueller_mat, fscat_mat, vdop )

        elseif(nspec==2) then

           call ss_sbm_polarris('ice_col',j,j,frequency,pts_temp_i_pol(i),&
                                den_sbm(1:nbin)%ice_col,rad_sbm(1:nbin)%ice_col,drad_sbm(1:nbin)%ice_col,&
                                vt_sbm(1:nbin)%ice_col, q_unit, n_unit(j:j), fmelt_fake(j:j), frime_fake(j:j), &
                                mxelev_polarris, elev_polarris, &
                                num_axis_sbm%ice_col, ityp_axis_sbm%ice_col, axis_ratio_sbm%ice_col, ityp_dist_sbm%ice_col,  &
                                mean_orient_sbm%ice_col, num_stdev_sbm%ice_col, stdev_orient_sbm%ice_col, &
                                bmueller_mat, fscat_mat, vdop  )

        elseif(nspec==3) then

           call ss_sbm_polarris('ice_pla',j,j,frequency,pts_temp_i_pol(i),&
                                den_sbm(1:nbin)%ice_pla,rad_sbm(1:nbin)%ice_pla,drad_sbm(1:nbin)%ice_pla,&
                                vt_sbm(1:nbin)%ice_pla, q_unit, n_unit(j:j), fmelt_fake(j:j), frime_fake(j:j), &
                                mxelev_polarris, elev_polarris, &
                                num_axis_sbm%ice_pla, ityp_axis_sbm%ice_pla, axis_ratio_sbm%ice_pla, ityp_dist_sbm%ice_pla,  &
                                mean_orient_sbm%ice_pla, num_stdev_sbm%ice_pla, stdev_orient_sbm%ice_pla, &
                                bmueller_mat, fscat_mat, vdop )

        elseif(nspec==4) then

           call ss_sbm_polarris('ice_den',j,j,frequency,pts_temp_i_pol(i),&
                                den_sbm(1:nbin)%ice_den,rad_sbm(1:nbin)%ice_den,drad_sbm(1:nbin)%ice_den,&
                                vt_sbm(1:nbin)%ice_den, q_unit, n_unit(j:j), fmelt_fake(j:j), frime_fake(j:j), &
                                mxelev_polarris, elev_polarris, &
                                num_axis_sbm%ice_den, ityp_axis_sbm%ice_den, axis_ratio_sbm%ice_den, ityp_dist_sbm%ice_den,  &
                                mean_orient_sbm%ice_den, num_stdev_sbm%ice_den, stdev_orient_sbm%ice_den, &
                                bmueller_mat, fscat_mat, vdop )


        elseif(nspec==5) then

           call ss_sbm_polarris('graupel',j,j,frequency,pts_temp_i_pol(i),&
                                den_sbm(1:nbin)%graupel,rad_sbm(1:nbin)%graupel,drad_sbm(1:nbin)%graupel,&
                                vt_sbm(1:nbin)%graupel, q_unit, n_unit(j:j), fmelt_fake(j:j), frime_fake(j:j), &
                                mxelev_polarris, elev_polarris, &
                                num_axis_sbm%graupel, ityp_axis_sbm%graupel, axis_ratio_sbm%graupel, ityp_dist_sbm%graupel,  &
                                mean_orient_sbm%graupel, num_stdev_sbm%graupel, stdev_orient_sbm%graupel, &
                                bmueller_mat, fscat_mat, vdop ) ! graupel

        elseif(nspec==6) then

           call ss_sbm_polarris('hail   ',j,j,frequency,pts_temp_i_pol(i),&
                                den_sbm(1:nbin)%hail,rad_sbm(1:nbin)%hail,drad_sbm(1:nbin)%hail,&
                                vt_sbm(1:nbin)%hail, q_unit, n_unit(j:j), fmelt_fake(j:j), frime_fake(j:j), &
                                mxelev_polarris, elev_polarris, &
                                num_axis_sbm%hail, ityp_axis_sbm%hail, axis_ratio_sbm%hail, ityp_dist_sbm%hail,  &
                                mean_orient_sbm%hail, num_stdev_sbm%hail, stdev_orient_sbm%hail, &
                                bmueller_mat, fscat_mat, vdop ) ! hail

        elseif(nspec >= 7 .and. nspec <= 16) then

           !varying riming fraction from 0., 0.1 ~ 0.9
           frime_temp(j:j) = real(nspec - 7) * 0.1  !7: frime=0.0 , 8: frime=0.1, 9: frime=0.2..., 16: frime=0.9
           call ss_sbm_polarris('snow   ',j,j,frequency,pts_temp_i_pol(i),&
                                den_sbm(1:nbin)%snow,rad_sbm(1:nbin)%snow,drad_sbm(1:nbin)%snow,&
                                vt_sbm(1:nbin)%snow, q_unit, n_unit(j:j), fmelt_fake(j:j), frime_temp(j:j), &
                                mxelev_polarris, elev_polarris, &
                                num_axis_sbm%snow, ityp_axis_sbm%snow, axis_ratio_sbm%snow, ityp_dist_sbm%snow,  &
                                mean_orient_sbm%snow, num_stdev_sbm%snow, stdev_orient_sbm%snow, &
                                bmueller_mat, fscat_mat, vdop  )


        else
           call stop_sdsu(myrank,'MSG makelut_polarris: there is no such nspec for SBM')
        endif NSPEC_SBM


      case default 
        call stop_sdsu(myrank,'MSG makelut_polarris: no such microphysics option.')
      end select mic_select2

      !
      !debbug toshii
      !
      !call derive_output_polradar( bmueller_mat(1,1:4,1:4), fscat_mat(1,1:2,1:2), dble(wave_m), polout )
      !print*,j,'zhh=',polout%zhh,' zdr=',polout%zdr , 'kdp =', polout%kdp



     !
     ! store output in the giant array
     !
     do m = 1, mxelev_polarris
        bmueller_array  (1:4,1:4,j,nspec,i,m) = bmueller_mat(m,1:4,1:4)

!if(nspec==6 .and. m==1 .and. i==mxpts_temp_pol) print*, 'jre=',j,nspec,i,m, 'hail bmueller_array  (1,1, j, 6, mxpts_temp_pol, 1)', &
!      bmueller_array  (1,1, j, nspec, mxpts_temp_pol, m), 'bmueller_mat(m,1,1)=',bmueller_mat(m,1,1)

        fscat_real_array(1:2,1:2,j,nspec,i,m) = real( fscat_mat(m,1:2,1:2) )
        fscat_imag_array(1:2,1:2,j,nspec,i,m) = aimag( fscat_mat(m,1:2,1:2) )
        vdop_array      (        j,nspec,i,m) = vdop     (m)
     enddo


 enddo giant_loop  ! nspec-species, i-temperature, j-effective radius loop

!
! MPI gather parameters to mater proc.
!
#if MPI > 0
 call mpi_wait_for_master !this process hold slave processors to wait for master proc. 

 if(numproc_tot>1) then

  do im = 1, 4 ; do jm = 1, 4
    call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined_dble, bmueller_array(im,jm,:,:,:,:) )
  enddo ; enddo

  do im = 1, 2 ; do jm = 1, 2
    call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, fscat_real_array(im,jm,:,:,:,:) )
    call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, fscat_imag_array(im,jm,:,:,:,:) )
  enddo ; enddo

  call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined_dble, vdop_array(:,:,:,:) )

 endif

 call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
#endif

!
! dump output
!
 if(masterproc) then

print*,'bmueller_array after MPI gather'
do n = 1, mxpts_re_pol_temp
  print*, bmueller_array  (1,1, n, 2, mxpts_temp_pol, 1), bmueller_array  (1,1, n, 6, mxpts_temp_pol, 1)
enddo

   call check( nf90_put_var(ncid, id_bmueller   , bmueller_array   ) )
   call check( nf90_put_var(ncid, id_fscat_real , fscat_real_array ) )
   call check( nf90_put_var(ncid, id_fscat_imag , fscat_imag_array ) ) !debug found here
   call check( nf90_put_var(ncid, id_vdop       , vdop_array       ) )



   !
   ! adding polarimetic variables
   !

   do m = 1, mxelev_polarris ; do i = 1, mxpts_temp_pol ; do nspec = 1,max_spc ; do j = 1, mxpts_re_pol_temp

      bmueller_tot(1:4,1:4) = bmueller_array      (1:4,1:4, j,nspec,i,m)
      fscat_tot(1:2,1:2) = cmplx( fscat_real_array(1:2,1:2, j,nspec,i,m), &
                                  fscat_imag_array(1:2,1:2, j,nspec,i,m)  )
      !
      ! derive polarimetric variable from mueller scattering matrix. 
      !
      call derive_output_polradar( bmueller_tot, fscat_tot, dble(wave_m), polout )

      array_zhh(j,nspec,i,m) = real( polout%zhh )   ! reflectivity [dBZ] 
      array_zdr(j,nspec,i,m) = real( polout%zdr )   ! differential reflectivity [dB]
      array_kdp(j,nspec,i,m) = real( polout%kdp )   ! specific diff phase [deg/km] 

   enddo ; enddo ; enddo ; enddo

   call check( nf90_put_var(ncid, id_zhh       , array_zhh  ) )
   call check( nf90_put_var(ncid, id_zdr       , array_zdr  ) )
   call check( nf90_put_var(ncid, id_kdp       , array_kdp  ) )

   call check( nf90_close(ncid) )  ! close nc file

   print*,'MGS makelut_polarris: Congraturation! Finished generating POLARRIS LUT'

 endif

 enddo freq_loop


 enddo ensemble_loop


 return
 end subroutine makelut_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine makelut_polarris_p3(nloop_max, loop_ntemp, loop_rhor, loop_fr, loop_qnorm, &
                                exit_loop )
 implicit none

 integer,intent(in) :: nloop_max !loop maxima
 integer,intent(in) :: loop_ntemp(nloop_max)
 integer,intent(in) :: loop_rhor (nloop_max)
 integer,intent(in) :: loop_fr   (nloop_max)
 integer,intent(in) :: loop_qnorm(nloop_max)

 logical,intent(in) :: exit_loop  !loop exit 


!
! local parameters
!
 character(len=100) :: fname_out
 character(len=100):: char_freq
 real(sdsu_fps) :: rhoe       ! re-dependent bulk snow density [kg/m3]
 real(sdsu_fps) :: frequency
 integer :: nspec,i,j,l,m,nfreq, im, jm, nf
 logical :: fexist
 integer,parameter :: io = 2
 integer :: ierr, n
 integer :: i_temp, i_rhor,i_Fr, i_Qnorm

! local
 real(sdsu_fpd) :: bmueller_tot(4,4)  !total mueller matrix
 complex :: fscat_tot(2,2)            !total forward scattering 
 type( polradar_parameter ) :: polout  !this is defined in module_mueller_matrix.F
 real(sdsu_fps) :: wave_m  !wavelength [m]
 real(sdsu_fps),allocatable,dimension(:,:,:,:,:) :: array_zhh, array_zdr, array_kdp

!
! LUT output
!
 real(sdsu_fpd) :: bmueller_mat(mxelev_polarris,4,4) ! 4x4 mueller matrix (backscatering component only)
 complex        :: fscat_mat   (mxelev_polarris,2,2) ! 2x2 forward scattering matrix 
 real(sdsu_fpd) :: vdop        (mxelev_polarris    )  ! Doppler velocity [m/s]

!
! NetCDF parameters
!
 integer :: ncid
 integer :: i_dimid, j_dimid, n_dimid, m_dimid, im_dimid, jm_dimid, is_dimid, js_dimid, &
            nr_dimid, nf_dimid, nq_dimid

 integer, parameter :: ndims7d = 7, ndims5d = 5
 integer :: dimids7d_mueller(ndims7d), dimids7d_scat(ndims7d), dimids5d_vdop(ndims5d)

 integer :: id, id_bmueller, id_fscat_real, id_fscat_imag, id_vdop, id_zhh, id_zdr, id_kdp, id_rhohv

 integer :: nens, max_ensemble
 character(len=1) :: nchar

 real :: n_unit = 1.  !unit number [#/m3]
 real :: qnorm ![kg]

!
! program start
!
 if(masterproc) print*, '- MSG makelut_polarris_p3: LUT for POLARRIS (P3 ice microphysics)'

!
! allocate giant LUT array and initialize parameters
!
 if( allocated(bmueller_array_p3) ) then
       deallocate(bmueller_array_p3  , &
                  fscat_real_array_p3, &
                  fscat_imag_array_p3, &
                  vdop_array_p3      )
 endif
 allocate( bmueller_array_p3  (4,4,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
           fscat_real_array_p3(2,2,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
           fscat_imag_array_p3(2,2,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
           vdop_array_p3      (    n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
           stat=ierr )
 if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_polarris_p3: allocation error bunded array. ')
 bmueller_array_p3   = undefined_dble
 fscat_real_array_p3 = undefined
 fscat_imag_array_p3 = undefined
 vdop_array_p3       = undefined_dble

 if( masterproc ) then
   allocate( array_zhh(n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
             array_zdr(n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
             array_kdp(n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &
            stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_polarris_p3: allocation error polout array. ')
  array_zhh = undefined
  array_zdr = undefined
  array_kdp = undefined
 endif


!
! dolan 2013 ensemble assumptions of axis ratio and orienting angle
!
if( polarris_axisorient_name == 'dolan2013' ) then
  max_ensemble = 4
else
  max_ensemble = 1
endif


!
! ensemble loop
!
 ensemble_loop: do nens = 1, max_ensemble

  if( polarris_axisorient_typ == 4 ) then
      write(nchar,"(I1.1)") nens ! 1 digit decimal 
      polarris_axisorient_name = 'dolan2013_pt'//nchar
  endif

!
! initialize POLARIS parameters
!
 call read_polarris_axisorient


!
! frequency loop
!
 freq_loop : do nfreq = 1 , mxinst_polarris

   !
   ! if same frequencies in freq loop, cycle loop
   !
   if(nfreq >= 2) then
     do nf = 1, nfreq-1
        if( freq_polarris(nfreq) == freq_polarris(nf) ) cycle freq_loop
     enddo
   endif

   frequency = freq_polarris(nfreq)
   if(masterproc) write(*,*) 'Freq = ',freq_polarris(nfreq),' GHz'
   if(masterproc) write(*,*) 'Wavelength = ',0.1*const_ghz_mm/freq_polarris(nfreq),' cm'

   wave_m = 0.1*const_ghz_mm/freq_polarris(nfreq)/100.  !wavelength [m]

   char_freq = nch_polarris(nfreq)


 if_masterproc: if(masterproc) then

!
! open netcdf file
!
 fname_out = 'POLARRIS_'//trim(cloud_microphysics)//'_ICE-'//trim(char_freq)//'_'//trim(polarris_axisorient_name)//'.nc'
 sdsu_io_file = trim(sdsu_dir_sslut)//trim(fname_out)

 print*,'MGS makelut_polarris_p3: start creating POLARRIS LUT for P3 ice',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )               !32bit output

!
! Define netcdf dimensions
!  allocate( bmueller_array_p3  (4,4,n_rhor,n_Fr,n_Qnorm,mxpts_temp_pol,mxelev_polarris), &

 call check( nf90_def_dim(ncid, 'melev' , mxelev_polarris, m_dimid) )
 call check( nf90_def_dim(ncid, 'itemp' , mxpts_temp_pol , i_dimid) )
 call check( nf90_def_dim(ncid, 'nqnorm', n_Qnorm        , nq_dimid) )
 call check( nf90_def_dim(ncid, 'nfr'   , n_Fr           , nf_dimid) )
 call check( nf90_def_dim(ncid, 'nrhor' , n_rhor         , nr_dimid) )
 call check( nf90_def_dim(ncid, 'imueller' , 4, im_dimid) )
 call check( nf90_def_dim(ncid, 'jmueller' , 4, jm_dimid) )
 call check( nf90_def_dim(ncid, 'iscat' , 2, is_dimid) )
 call check( nf90_def_dim(ncid, 'jscat' , 2, js_dimid) )

 dimids7d_mueller = (/ im_dimid, jm_dimid, nr_dimid, nf_dimid, nq_dimid, i_dimid, m_dimid /)
 dimids7d_scat    = (/ is_dimid, js_dimid, nr_dimid, nf_dimid, nq_dimid, i_dimid, m_dimid /)
 dimids5d_vdop    = (/                     nr_dimid, nf_dimid, nq_dimid, i_dimid, m_dimid /)


!
! global attribute 
!

!
! ice axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_ice', &
                        ityp_axis_gen(1:num_axis_gen%ice)%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_ice', &
                        axis_ratio_gen(1:num_axis_gen%ice)%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_ice', ityp_dist_gen%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_ice', mean_orient_gen%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice', &
                       stdev_orient_gen(1:num_stdev_gen%ice)%ice  ) )

!
! snow axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_snow', &
                        ityp_axis_gen(1:num_axis_gen%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_snow', &
                        axis_ratio_gen(1:num_axis_gen%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_snow', ityp_dist_gen%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_snow', mean_orient_gen%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_gen_snow', &
                       stdev_orient_gen(1:num_stdev_gen%snow)%snow  ) )

!
! graupel axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_graupel', &
                        ityp_axis_gen(1:num_axis_gen%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_graupel', &
                        axis_ratio_gen(1:num_axis_gen%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_graupel', ityp_dist_gen%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_graupel', mean_orient_gen%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_graupel', &
                       stdev_orient_gen(1:num_stdev_gen%graupel)%graupel  ) )

!
! hail axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_hail', &
                        ityp_axis_gen(1:num_axis_gen%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_hail', &
                        axis_ratio_gen(1:num_axis_gen%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_hail', ityp_dist_gen%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_hail', mean_orient_gen%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_hail', &
                       stdev_orient_gen(1:num_stdev_gen%hail)%hail  ) )


!
! define netcdf output parameters
!
 call check( nf90_def_var(ncid, 'bmueller', NF90_DOUBLE, dimids7d_mueller, id) )
 call check( nf90_put_att(ncid, id, 'description',' 4x4 mueller matrix ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_bmueller = id

 call check( nf90_def_var(ncid, 'fscat_real', NF90_FLOAT, dimids7d_scat, id) )
 call check( nf90_put_att(ncid, id, 'description',' 2x2 scattering matrix (real)') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_fscat_real = id

 call check( nf90_def_var(ncid, 'fscat_imag', NF90_FLOAT, dimids7d_scat, id) )
 call check( nf90_put_att(ncid, id, 'description',' 2x2 scattering matrix (imagenary)') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_fscat_imag = id

 call check( nf90_def_var(ncid, 'vdop', NF90_DOUBLE, dimids5d_vdop, id) )
 call check( nf90_put_att(ncid, id, 'description',' doppler velocity ') )
 call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
 id_vdop = id


 call check( nf90_def_var(ncid, 'zhh', NF90_FLOAT, dimids5d_vdop, id) )
 call check( nf90_put_att(ncid, id, 'description',' reflectivity at H pol ') )
 call check( nf90_put_att(ncid, id, 'units', '[dBZ]') )
 id_zhh = id

 call check( nf90_def_var(ncid, 'zdr', NF90_FLOAT, dimids5d_vdop, id) )
 call check( nf90_put_att(ncid, id, 'description',' differential reflectivity ') )
 call check( nf90_put_att(ncid, id, 'units', '[dB]') )
 id_zdr = id

 call check( nf90_def_var(ncid, 'kdp', NF90_FLOAT, dimids5d_vdop, id) )
 call check( nf90_put_att(ncid, id, 'description',' specific diff phase ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg/km]') )
 id_kdp = id


 call check( nf90_enddef(ncid) )  ! End define mode.

 endif if_masterproc

!
! giant loop
!
   giant_loop: do l = myl_start, myl_end

     if( exit_loop ) cycle  !this ensure to cycle proces, when processor is not assigned for task in MPI=1 options

     !
     ! asign original nspec, i, j parameters from giant loop index
     !
     i       = loop_ntemp(l)  ! temperature index
     i_rhor  = loop_rhor (l)
     i_Fr    = loop_fr   (l) 
     i_Qnorm = loop_qnorm(l) 


     ! zero out output parameters
     bmueller_mat = 0. ; fscat_mat = 0. ; vdop = 0.


     qnorm = 261.7**(real(i_Qnorm+10)*0.1)*1.e-18 ! Qtot/Ntot [kg]
 
     call ss_p3_polarris(frequency, pts_temp_i_pol(i), qnorm, n_unit, Fr(i_Fr), rhor(i_rhor) , & !p3 ice
                         mxelev_polarris, elev_polarris, &
                         bmueller_mat, fscat_mat, vdop)

     !
     !debbug toshii
     !
     !call derive_output_polradar( bmueller_mat(1,1:4,1:4), fscat_mat(1,1:2,1:2), dble(wave_m), polout )
     !print*,j,'zhh, zdr, kdp =', polout%zhh, polout%zdr, polout%kdp


     !
     ! store output in the giant array (this is normzlied by ntot)
     !
     do m = 1, mxelev_polarris
        bmueller_array_p3  (1:4,1:4,i_rhor,i_Fr,i_Qnorm,i,m) =     bmueller_mat(m,1:4,1:4)
        fscat_real_array_p3(1:2,1:2,i_rhor,i_Fr,i_Qnorm,i,m) = real ( fscat_mat(m,1:2,1:2) )
        fscat_imag_array_p3(1:2,1:2,i_rhor,i_Fr,i_Qnorm,i,m) = aimag( fscat_mat(m,1:2,1:2) )
        vdop_array_p3      (        i_rhor,i_Fr,i_Qnorm,i,m) =             vdop(m)

!if(nspec==6 .and. m==1 .and. i==mxpts_temp_pol) print*, 'jre=',j,nspec,i,m, 'hail bmueller_array  (1,1, j, 6, mxpts_temp_pol, 1)', &
!      bmueller_array  (1,1, j, nspec, mxpts_temp_pol, m), 'bmueller_mat(m,1,1)=',bmueller_mat(m,1,1)
     enddo

 enddo giant_loop  ! nspec-species, i-temperature, j-effective radius loop


!
! MPI gather parameters to mater proc.
!
#if MPI > 0
 call mpi_wait_for_master !this process hold slave processors to wait for master proc. 

 if(numproc_tot>1) then

  do im = 1, 4 ; do jm = 1, 4; do i_rhor = 1, n_rhor
    call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined_dble, bmueller_array_p3(im,jm,i_rhor,:,:,:,:) )
  enddo ; enddo ; enddo

  do im = 1, 2 ; do jm = 1, 2 ; do i_rhor = 1, n_rhor
    call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, fscat_real_array_p3(im,jm,i_rhor,:,:,:,:) )
    call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, fscat_imag_array_p3(im,jm,i_rhor,:,:,:,:) )
  enddo ; enddo ; enddo

  do i_rhor = 1, n_rhor
    call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined_dble, vdop_array_p3(i_rhor,:,:,:,:) )
  enddo 
 endif

 call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
#endif



!
! dump output
!
 if(masterproc) then


 print*,'bmueller_array_p3 after MPI gather'
 do i_rhor = 1, n_rhor ; do i_Fr = 1, n_Fr ; do i_Qnorm = 1, n_Qnorm
    print*, bmueller_array_p3  (1,1,i_rhor,i_Fr,i_Qnorm, mxpts_temp_pol, 1)
 enddo ; enddo ; enddo


   call check( nf90_put_var(ncid, id_bmueller   , bmueller_array_p3   ) )
   call check( nf90_put_var(ncid, id_fscat_real , fscat_real_array_p3 ) )
   call check( nf90_put_var(ncid, id_fscat_imag , fscat_imag_array_p3 ) ) 
   call check( nf90_put_var(ncid, id_vdop       , vdop_array_p3       ) )

   !
   ! deriving polarimetic variables for visualizing LUT purpose
   !

   do m = 1, mxelev_polarris ; do i = 1, mxpts_temp_pol ; do i_Qnorm = 1, n_Qnorm ; do i_Fr = 1, n_Fr ; do i_rhor = 1, n_rhor

      bmueller_tot(1:4,1:4) = bmueller_array_p3      (1:4,1:4, i_rhor,i_Fr,i_Qnorm,i,m)
      fscat_tot(1:2,1:2) = cmplx( fscat_real_array_p3(1:2,1:2, i_rhor,i_Fr,i_Qnorm,i,m), &
                                  fscat_imag_array_p3(1:2,1:2, i_rhor,i_Fr,i_Qnorm,i,m)  )
      !
      ! derive polarimetric variable from mueller scattering matrix. 
      !
      call derive_output_polradar( bmueller_tot, fscat_tot, dble(wave_m), polout )

      array_zhh(i_rhor,i_Fr,i_Qnorm,i,m) = real( polout%zhh )   ! reflectivity [dBZ] 
      array_zdr(i_rhor,i_Fr,i_Qnorm,i,m) = real( polout%zdr )   ! differential reflectivity [dB]
      array_kdp(i_rhor,i_Fr,i_Qnorm,i,m) = real( polout%kdp )   ! specific diff phase [deg/km] 

   enddo ; enddo ; enddo ; enddo ; enddo

   call check( nf90_put_var(ncid, id_zhh       , array_zhh  ) )
   call check( nf90_put_var(ncid, id_zdr       , array_zdr  ) )
   call check( nf90_put_var(ncid, id_kdp       , array_kdp  ) )

   call check( nf90_close(ncid) )  ! close nc file

   print*,'MGS makelut_polarris_p3: Congraturation! Finished generating POLARRIS LUT for P3 ice'

 endif

 enddo freq_loop

 enddo ensemble_loop


 return
 end subroutine makelut_polarris_p3

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


 subroutine ss_gen_polarris(spc, freq, temp, mu, wc, re, dens, frac_liq, &
                            num_elev, elev, num_axis, ityp_axis, axis_ratio, ityp_dist, & 
                            mean_orient, num_stdev, stdev_orient, &
                            a_vt,b_vt, bmueller_mat_out, fscat_mat_out, vdop_out )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Integrate over size distributions to derive total Mueller and forward scattering matrix for a given
!  size distributions, effective density, and axis ratio of particles of general bulk microphysics.
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: spc !species character
 real(sdsu_fps),intent(in) :: freq  ! frequency of radiation [GHz]
 real(sdsu_fps),intent(in) :: temp  ! temperature of particles [K]
 real(sdsu_fps),intent(in) :: mu    ! dispersion parameter [-]
 real(sdsu_fps),intent(in) :: wc    ! water content of condensate [g/m**3]
 real(sdsu_fps),intent(in) :: re    ! drop effective radius [micron] 
 real(sdsu_fps),intent(in) :: dens  ! bulk effective density of hydrometero [kg/m**3]
 real(sdsu_fps),intent(in) :: frac_liq ! fracion of liquid (1=rain or 0=ice condensates, 0~1=melting ice)

 integer,intent(in) :: num_elev   ! number of elevation angle
 real(sdsu_fps),intent(in) :: elev(num_elev)  ! elevaation angle [deg]
 integer,intent(in) :: num_axis   ! number of axis ratio 
 integer,intent(in) :: ityp_axis(num_axis)  ! particle axis ratio model type [-] 
 real(sdsu_fps),intent(in) :: axis_ratio(num_axis)  ! particle axis ratio [-]
 integer,intent(in) :: ityp_dist  ! particle orienting angle distribution type
 real(sdsu_fps),intent(in) :: mean_orient  ! mean orienting angle [deg]
 integer,intent(in) :: num_stdev  ! number of orienting angle standard deviations
 real(sdsu_fps),intent(in) :: stdev_orient(num_stdev)  ! standard deviation / theta_m of orienting angle [deg]

 real(sdsu_fps),intent(in) :: a_vt     ! alpha of Vt  [-]
 real(sdsu_fps),intent(in) :: b_vt     ! beta of Vt   [-]
 real(sdsu_fpd),intent(out) :: vdop_out(num_elev)! Doppler velocity [cm/s]

!--------Local Parameters
 integer :: i    ! looping 
 integer :: imax ! max for looping

 real(sdsu_fps) :: rhoe ! bulk effective density of hydrometero [kg/m**3]
 real(sdsu_fps) :: n0   ! intercept for the exponetial DSD [1/m**4]
 real(sdsu_fps) :: wave
 real(sdsu_fps) :: rad  ! radius of particle [mm]
 real(sdsu_fps) :: lam  ! the slope of the distribution [1/m]
 real(sdsu_fps) :: num  ! a particle number density per radius increment [1/m4]
!NSSL start
 real(sdsu_fps) :: num_v ! particle number density per volume increment [1/m6]
 real(sdsu_fps) :: ntot !total number conc [#/m3]
 real(sdsu_fps) :: d0   ! characteristic diameter [m] 
 real(sdsu_fps) :: dbar ! characteristics diameter [m]
 real(sdsu_fps) :: v0   ! characteristic volume [m3]
 real(sdsu_fps) :: v_increment ! volume increment [m3]
 real(sdsu_fps) :: v_m3 ! particle volume [m3]
 real(sdsu_fps) :: cd   ! drag coefficients
 real(sdsu_fps),parameter :: rhoair_800mb = 1.007  !dry air desnti of US standar air at 800mb [kg/m3]
!NSSL end
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
 real(sdsu_fps) :: fmelt
 real(sdsu_fps) :: vt_wsum  !weighting sum of terminal velocity 
 real(sdsu_fps) :: vt       !terminal velocity [m/s]


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
 integer,save :: nshp !SCATDB shape #
 real(sdsu_fps) :: qext_ns, qsca_ns, asym_ns, qbsca_ns ! non-spherical single scattering
 real(sdsu_fps) :: d_mm_solid  !solid ice diameter [mm]
 real(sdsu_fps) :: rad_solid_micron  !solid ice radius [micron]
 real(sdsu_fps) :: area   !area of fully particle [mm2]
 logical :: success ! status 
 real(sdsu_fps),allocatable :: thom_snow_num(:)  !size distributions of thompson scheme [1/m4]
 real(sdsu_fps) :: rad_new, fmelt_new

! tmatrix
 integer :: nrank, nm
 real(sdsu_fps) :: creal, cimag  !real / imageneary component of dielectric function
 real(sdsu_fpd) :: diam_sph        !spherical diaemter [cm] 
 real(sdsu_fpd) :: dwave           !wavelength of radar band [cm]
 real(sdsu_fpd) :: daxis_ratio, conk, aovrb, sigma_temp
 complex*8      :: tmat_out(20,40,40)   ! maybe dynamic allocation (save is debbug)
 real(sdsu_fps) :: axis
 integer :: naxis
 real(sdsu_fps) :: axis_temp

!
! mueller matrix
!
 integer :: nstdev , melev 
 integer,parameter :: nforward = 1 , nback = 5  !index of forward and backscatter angle in mueller matrix
 real(sdsu_fpd)    :: mueller_mat          (max_ens,max_elv,4,4) ! 4x4 mueller matrix  (backscatter component)
 real(sdsu_fpd) :: mueller_mat_sum         (max_ens,max_elv,4,4) ! summation of 4x4 mueller matrix 
 real(sdsu_fpd) :: mueller_mat_full(max_ens,max_ens,max_elv,4,4) ! full bundle of 4x4 mueller matrix 
 complex :: fscat_mat              (max_ens,max_elv,2,2) ! 2x2 forward scattering matrix 
 complex :: fscat_mat_sum          (max_ens,max_elv,2,2) ! summation of 2x2 forward scattering matrix 
 complex :: fscat_mat_full (max_ens,max_ens,max_elv,2,2) ! full bundle of 2x2 forward scattering matrix 
 real(sdsu_fpd),intent(out) :: bmueller_mat_out(num_elev,4,4) ! 4x4 mueller matrix (backscatering component only)
 complex,intent(out) :: fscat_mat_out (num_elev,2,2)    ! 2x2 forward scattering matrix 
 type ( polradar_parameter ) :: output ! polarimeterica radar measureable

 real(sdsu_fpd) :: bscat            ( max_ens,max_elv ) ! backscatter coeff
 real(sdsu_fpd) :: vt_sum           ( max_ens,max_elv ) ! sum of vt
 real(sdsu_fpd) :: vt_wgt           ( max_ens,max_elv ) ! weight of vt  
 real(sdsu_fpd) :: vdop_full( max_ens,max_ens,max_elv ) ! doppler velocity [m/s]

 type( polradar_parameter ) :: polout  !this is defined in module_mueller_matrix.F


 real(sdsu_fpd) :: dn0, dlam, dmu, dnum


!
! Begin by checking if hydrometeors of this species are present.
! If not, set scattering parameters to zero and return.
!
  if(wc .lt. q_min_micro) then
     bmueller_mat_out=0. ; fscat_mat_out=cmplx(0.,0.) ; vdop_out=0.
     return
  endif


!
! Assign some useful constants
!
  wave=const_ghz_mm/freq  !wavelength [mm]

!
! If hydrometeors are present, initialize the scattering parameters
!
 mueller_mat = 0. ; mueller_mat_sum = 0. ; mueller_mat_full = 0.
 fscat_mat  = cmplx(0.,0.) ; fscat_mat_sum  = cmplx(0.,0.) ; fscat_mat_full  = cmplx(0.,0.)
 vdop_full = 0. ;  vt_sum = 0.  ;  vt_wgt = 0.   

 if( trim(spc) == 'cloud' .or. trim(spc) == 'ice' )  then  !non-precip

! define imax for a givne radius increment 20. is maximum diameter
  d_increment = 0.01  ! [mm] --> 1micron
  d_min =  0.002  ![mm]  --> 2micron
  d_max = 0.2  ! [mm]  --> 200micron

 else  !precip species

! define imax for a givne radius increment 20. is maximum diameter
  d_increment = 1.0  ! [mm]  --> 100micron
  d_min =  0.2  ![mm]  --> 200micron
  d_max = 20.  ! [mm]  --> 2cm

 endif

!
! special case (Thompson snow class starts from ice crystal size)
!
  if( trim(cloud_microphysics) == 'THOM' .and. &
      trim(spc) == 'snow'  ) then  !thompson spectrum density of snow

     d_increment = 0.01  ! [mm] --> 10micron
     d_min =  0.002  ![mm]  --> 2micron
     d_max = 20.  ! [mm]  --> 2cm

  endif

 imax = nint( (d_max-d_min) /d_increment)

  if( trim(cloud_microphysics) == 'THOM' .and. &
      trim(spc) == 'snow'    ) then  !thompson spectrum density of snow

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
      trim(spc) == 'snow'  ) then  
     lam = 1. / (2.*re*1.e-6)  !slope [1/m]  Eq. (8) in Bae et al. (2016)
  endif

!
! derive intercept (became double precision)
!       [g/m3]                     [m3/g]        [1/m]
!  n0  = wc       * 6. / const_pi / (dens*1000.) * (lam**(4.+mu))  / gfac1  ! [1/m**(4+mu)]
  dn0 = dble(wc)       * 6. / dble(const_pi) / dble(dens*1000.) * (dble(lam)**(4.+dble(mu)))  / dble(gfac1)  ! 


!
! NSSL 2moment scheme
!
 if( trim(cloud_microphysics) == 'NSSL2MMT' ) then
    if( trim(spc) == 'cloud' .or. trim(spc) == 'ice' .or. trim(spc) == 'snow' ) then

      ! volume-based DSD
      call gamma_reff(mu+1.0,gfac1)
      call gamma_reff(mu+5.0/3.0,gfac2)
      lam = 0.5 * (1.+mu)*gfac1/gfac2/re * 1.e+6  !slope [1/m]
      dbar = 1./lam                   ! characteristics diameter [m]
      v0 = const_pi / 6. * dbar**3.   ! characteristic volume [m3]

      call gamma_reff(mu+2.0,gfac2)
      ntot = 6./const_pi*( lam**3 )/dens*gfac1/gfac2*(wc*1.e-3) ! total number conc. [#/m3]

    else !rain, graupel, and hail

      ! unique gamma DSD Ziegler (1985)
      call gamma_reff(mu+4.0,gfac1)
      call gamma_reff(mu+1.0,gfac2)
      gamfac = gfac2/gfac1   

      ntot =  (lam**3)*6.*wc/const_pi/(dens*1000.)*gamfac  ! total number conc [#/m3]
      d0 = 1./lam                                          ! characteristic diameter [m] 

    endif 
 endif


!
! loop for axix ratio
!
 loop_axis: do naxis = 1, num_axis
   axis_temp = axis_ratio(naxis) 


!
! size loop
!

   loop_size: do i=0,imax

     d_mm = d_min + d_increment * FLOAT(i)  !diameter [mm]   (new in V3.5)
     d_m  = d_mm * 1e-3                     ! diameter [m]
     xd  = dble(const_pi*d_mm/wave)         ! size parameter [-]

!
! Get N(D)
!         [1/m**(4+mu)]  [ m**mu ]      [1/m]  [m]   = [1/m4]
!     num = n0           * d_m**mu * exp(-lam * d_m )

    dnum = dn0          *dble(d_m)**dble(mu) * exp(-dble(lam) * dble(d_m) )

    num = real(dnum) !back to single precision


!
! Get N(D) for WDM double moment species only. 
!
     if( trim(cloud_microphysics) == 'WDM' ) then
        if( trim(spc) == 'cloud' ) call wdm_derive_dsd( 'cloud', wc*1.e-3, re*1.e-6 , d_m, num ) 
        if( trim(spc) == 'rain'  ) call wdm_derive_dsd( 'rain' , wc*1.e-3, re*1.e-6 , d_m, num )
     endif

!
! Get N(D) for Thompson snow number density 
!
     if( trim(cloud_microphysics) == 'THOM' .and. &
         trim(spc) == 'snow'     ) then
       num = thom_snow_num(i) ! [1/m4]
     endif

!
! NSSL 2moment scheme
!
 if( trim(cloud_microphysics) == 'NSSL2MMT' ) then
    if( trim(spc) == 'cloud' .or. trim(spc) == 'ice' .or. trim(spc) == 'snow' ) then

     v_m3 = const_pi / 6. * d_m**3  !particle volume [m3]
! volume-based number density [1/m6]
!    [1/m6]   [1/m3]                                                                                   [1/m3]
      num_v = ntot * ( (mu+1.)**(mu+1.) ) * ( ( v_m3 / v0 )**mu ) * exp( -1.*(mu+1.)*v_m3/v0 ) / gfac1 / v0

    else !rain, graupel, and hail


!    [1/m4] [1/m3]   [1/m^1+mu]            [m^mu]        
      num =  ntot  / d0**(1.+mu) / gfac2 * d_m**mu * exp( -1.*d_m/d0 )  !number density [1/m4]

    endif
 endif


!
! upadte meling fraction
!
     fmelt = frac_liq

!
! Calculate effective dielectric constant of frozen particle using different methods/assumptions. 
!
     if( trim(cloud_microphysics) == 'THOM' .and. &
          trim(spc) == 'snow'    ) then  !thompson spectrum density of snow
          rhoe = min( densice,  (0.069 * d_m*d_m) / ( const_pi / 6. * d_m*d_m*d_m ) ) ! [kg/m3]
     else
          rhoe = dens                  ! effective density [kg/m3]
     endif


!
! Snow compatction impact (Li-Matsui scheme)
!
     if(trim(spc) == 'snow' .and. frac_liq > 0.0 .and. melt_opt /= 0) then
        call snow_compaction( 0.5*d_mm , rhoe*1e-3 , frac_liq , rad_new, fmelt_new )
        d_mm = 2.*rad_new         ! compacted diameter [mm]
        d_m  = d_mm * 1e-3        ! diameter but using different unit [m]
        xd  = dble(const_pi*d_mm/wave)  ! size parameter [-]
        fmelt = fmelt_new         ! update liquid fraction
     endif

     !
     ! derive effective permittivity (liquid, solid, mixed phases)
     !
     call  effective_permittivity(freq, fmelt, rhoe, temp, crefd, &
                                  creal,cimag )

     !
     ! derive size-dependent axis ratio, if specified function id.
     !
     call spectrum_axis_ratio( trim(spc) , ityp_axis(naxis), d_mm*0.1, axis_temp )
     daxis_ratio = dble(axis_temp)

     !
     ! drive tmatrix :primal output is tmat_out(20,40,40)
     !
     diam_sph = dble(d_mm) * 0.1d0 ! spherical diaemter [cm]
     dwave  = dble(wave) * 0.1d0   ! wavelength [cm]


     call driver_tmatrix (daxis_ratio, dble(creal), dble(cimag), diam_sph, dwave, &   ! input
                          tmat_out, nrank, nm, conk, aovrb, sigma_temp )  ! output


     !
     ! drive Mueller matrix to derive 4x4 mueller matrix and 2x2 forward scattering matrix
     ! We also need extinciton and backscatter from this.....
     !
     call driver_mueller( daxis_ratio, dble(creal), dble(cimag), diam_sph, dwave, &   ! input
                          ityp_dist, dble(mean_orient), num_stdev, dble(stdev_orient) , &
                          tmat_out, conk, sigma_temp, num_elev, dble(elev(1:num_elev)) , &
                          mueller_mat, fscat_mat  )

      !
      !debbug
      !
      !call derive_output_polradar( mueller_mat(1,1,1:4,1:4), fscat_mat(1,1,1:2,1:2), dble(wave*1.e3), polout )
      !print*,daxis_ratio,creal,cimag,diam_sph,dwave
      !print*,trim(spc),i,'d_mm=',d_mm,'zhh, zdr, kdp =', polout%zhh, polout%zdr, polout%kdp


     !
     ! derive backscattering component 
     !
     do nstdev = 1, num_stdev ; do melev = 1, num_elev
        call bscat_mueller( mueller_mat ( nstdev, melev, 1:4,1:4 ), bscat(nstdev, melev) )
     enddo ; enddo

     !
     ! terminal velocity [m/s]
     !
     vt = a_vt * (d_m**b_vt)  ![m/s]  toshii --> need upper limit of vt_max to inhibit large values. 

     !
     ! P3 scheme
     !
     if( trim(cloud_microphysics) == 'P3' ) then
        call p3_vt_liq(spc, d_m, vt) 
     endif

     !
     ! NSSL 2moment scheme (revise above vt according to the NSSL definition)
     !
     if( trim(cloud_microphysics) == 'NSSL2MMT' ) then
         select case( trim(spc) ) 
         case('rain')

           vt = 10.*(1.-exp(-516.575*d_m))   !rain Vt (Mansel --> interpretation of Gunn and Kinzer data) [m/s]

         case('graupel_017','graupel_030','graupel_040','graupel_050',&
              'graupel_060','graupel_070','graupel_080','graupel_090')

           cd = max(0.45, min(1.2, 0.45 + 0.55*(800. - max( 170., min( 800., rhoe ) ) )/(800. - 170.) ) )
           vt = sqrt(4.0*rhoe*const_g/(3.0*cd*rhoair_800mb)) * d_m**0.5  ! Wisner form [m/s]

         case('hail_050','hail_060','hail_070','hail_080','hail_090') 

           cd = max(0.45, min(1.2, 0.45 + 0.55*(800. - max( 500., min( 800., rhoe ) ) )/(800. - 170.) ) )
           vt = sqrt(4.0*rhoe*const_g/(3.0*cd*rhoair_800mb)) * d_m**0.5  ! Wisner form [m/s]

         end select 

     endif

!
! integrate over particle size distribution; (Cross sectin extionction)*(cross sectional aerea)
!
!    [1/m3] =[1/m4]*[mm]        *[m/mm] --> [1/m3] 
     factor =  num * d_increment * 1.e-3


    if( trim(cloud_microphysics) == 'NSSL2MMT' ) then
      if( trim(spc) == 'cloud' .or. trim(spc) == 'ice' .or. trim(spc) == 'snow' ) then

        v_increment = const_pi / 6. * ( d_increment * 1.e-3 )**3.  !volume increment [m3]

        ![1/m3] =  [1/m6]* [m3]         
         factor =  num_v * v_increment 
 
      endif
    endif


     mueller_mat_sum = mueller_mat_sum + mueller_mat * factor !accumulating 4x4 mueller matrix

     fscat_mat_sum  = fscat_mat_sum  + fscat_mat  * factor !accumulating 2x2 forward scattering matrix

     do nstdev = 1, num_stdev ; do melev = 1, num_elev
        vt_sum( nstdev, melev ) = vt_sum( nstdev, melev ) +  vt * factor * bscat(nstdev, melev)
        vt_wgt( nstdev, melev ) = vt_wgt( nstdev, melev ) +       factor * bscat(nstdev, melev)
!        print*, 'MSG ss_gen_micro:', i, melev, factor, bscat(nstdev, melev)
     enddo ; enddo
!print*, 'mueller_mat=',mueller_mat(1,1,1:4,1:4)
!print*,'factor',factor, num, d_increment
!print*, 'mueller_mat_sum=',mueller_mat_sum (1,1,1:4,1:4)  ! mueller_mat_sum         (max_ens,max_elv,4,4)


   enddo loop_size

!
! for each axis ratio loop
!
   mueller_mat_full(naxis,1:num_stdev,1:num_elev,1:4,1:4 ) = &
             mueller_mat_sum       (1:num_stdev,1:num_elev,1:4,1:4 ) !output 4x4 mueller matrix


   fscat_mat_full(naxis,1:num_stdev,1:num_elev,1:2,1:2 )  = &
             fscat_mat_sum       (1:num_stdev,1:num_elev,1:2,1:2 ) !output 2x2 forward scattering matrix

   do nstdev = 1, num_stdev ; do melev = 1, num_elev 
      if( vt_wgt( nstdev, melev ) <= 0.0 ) then
          vdop_full(naxis,nstdev, melev) = 0.
      else    
          vdop_full(naxis,nstdev, melev) = vt_sum( nstdev, melev ) / vt_wgt( nstdev, melev ) !doppler velocity [m/s] 
      endif
   enddo ; enddo

 enddo loop_axis

!==================================================================================
! toshii --> eventually we must derive ensemble mean, max, minima component or stdev
! right now, just 1st memmber of ensemble
!

 do melev = 1, num_elev

   !
   !forward and backscatter component of mueller matrix (so far no ensembles)
   !
   bmueller_mat_out(melev,1:4,1:4) =  mueller_mat_full(1,1,melev,  1:4,1:4 ) 

   !
   ! forward scattering matrix output 
   !
   fscat_mat_out (melev,1:2,1:2) =  fscat_mat_full(1,1,melev,     1:2,1:2 )  

   !
   ! doppler velocity 
   !
   vdop_out(melev) = vdop_full(1,1,melev) 

   !
   ! derive polarimetric parameters (testing only)
   !
   ! call derive_output_polradar( bmueller_mat_out(1:4,1:4), fscat_mat_out (1:2,1:2), wave*0.001d0, &
   !                              output  )

 enddo

 return
 end subroutine ss_gen_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine ss_p3_polarris( freq, temp, qtot, ntot, frim, rhorim, &
                            num_elev, elev, &
                            bmueller_mat_out, fscat_mat_out, vdop_out )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Integrate over size distributions to derive total Mueller and forward scattering matrix for a given
!  size distributions, effective density, and axis ratio of particles of general bulk microphysics.
!
! History:
! 06/2019  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: freq  ! frequency of radiation [GHz]
 real(sdsu_fps),intent(in) :: temp  ! temperature of particles [K]
 real(sdsu_fps),intent(in) :: qtot  ! total mass conc. [g/m3]
 real(sdsu_fps),intent(in) :: ntot  ! total number conc. [#/m3] 
 real(sdsu_fps),intent(in) :: frim   ! rime mass fraction (Qrim/Qtot=Fr) [-]
 real(sdsu_fps),intent(in) :: rhorim ! rime density (Qrim/Brim=RHOr) [kg/m3]

 integer,intent(in) :: num_elev   ! number of elevation angle
 real(sdsu_fps),intent(in) :: elev(num_elev)  ! elevaation angle [deg]

 real(sdsu_fpd),intent(out) :: bmueller_mat_out(num_elev,4,4) ! 4x4 mueller matrix (backscatering component only)
 complex,intent(out) :: fscat_mat_out (num_elev,2,2)    ! 2x2 forward scattering matrix 
 real(sdsu_fpd),intent(out) :: vdop_out(num_elev)! Doppler velocity [cm/s]

!--------Local Parameters
 integer :: i    ! looping 
 integer :: imax ! max for looping

 real(sdsu_fps) :: rhoe ! bulk effective density of hydrometero [kg/m**3]
 real(sdsu_fps) :: n0   ! intercept for the exponetial DSD [1/m**4]
 real(sdsu_fps) :: wave
 real(sdsu_fps) :: rad  ! radius of particle [mm]
 real(sdsu_fps) :: lam  ! the slope of the distribution [1/m]
 real(sdsu_fps) :: num  ! a particle number density per radius increment [1/m4]
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
 real(sdsu_fps) :: fmelt
 real(sdsu_fps) :: vt_wsum  !weighting sum of terminal velocity 
 real(sdsu_fps) :: vt       !terminal velocity [m/s]


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
 integer,save :: nshp !SCATDB shape #
 real(sdsu_fps) :: qext_ns, qsca_ns, asym_ns, qbsca_ns ! non-spherical single scattering
 real(sdsu_fps) :: d_mm_solid  !solid ice diameter [mm]
 real(sdsu_fps) :: rad_solid_micron  !solid ice radius [micron]
 real(sdsu_fps) :: area   !area of fully particle [mm2]
 logical :: success ! status 
 real(sdsu_fps),allocatable :: thom_snow_num(:)  !size distributions of thompson scheme [1/m4]
 real(sdsu_fps) :: rad_new, fmelt_new
 real(sdsu_fps) :: coef ! coefficient for power law


 real(sdsu_fps) :: axis_ratio  ! particle axis ratio [-]
 integer :: ityp_dist  ! particle orienting angle distribution type
 real(sdsu_fps) :: mean_orient  ! mean orienting angle [deg]
 integer :: num_stdev  ! number of orienting angle standard deviations
 real(sdsu_fps) :: stdev_orient(max_ens)  ! standard deviation / theta_m of orienting angle [deg]
 integer :: p3_ice_regime  !P3 ice regime (1-4)


! tmatrix
 integer :: nrank, nm
 real(sdsu_fps) :: creal, cimag  !real / imageneary component of dielectric function
 real(sdsu_fpd) :: diam_sph        !spherical diaemter [cm] 
 real(sdsu_fpd) :: dwave           !wavelength of radar band [cm]
 real(sdsu_fpd) :: daxis_ratio, conk, aovrb, sigma_temp
 complex*8      :: tmat_out(20,40,40)   ! maybe dynamic allocation (save is debbug)
 real(sdsu_fps) :: axis
 integer :: naxis
 real(sdsu_fps) :: axis_temp

!
! mueller matrix
!
 integer :: nstdev , melev 
 integer,parameter :: nforward = 1 , nback = 5  !index of forward and backscatter angle in mueller matrix
 real(sdsu_fpd)    :: mueller_mat          (max_ens,max_elv,4,4) ! 4x4 mueller matrix  (backscatter component)
 real(sdsu_fpd) :: mueller_mat_sum         (max_ens,max_elv,4,4) ! summation of 4x4 mueller matrix 
 real(sdsu_fpd) :: mueller_mat_full(max_ens,max_ens,max_elv,4,4) ! full bundle of 4x4 mueller matrix 
 complex :: fscat_mat              (max_ens,max_elv,2,2) ! 2x2 forward scattering matrix 
 complex :: fscat_mat_sum          (max_ens,max_elv,2,2) ! summation of 2x2 forward scattering matrix 
 complex :: fscat_mat_full (max_ens,max_ens,max_elv,2,2) ! full bundle of 2x2 forward scattering matrix 
 type ( polradar_parameter ) :: output ! polarimeterica radar measureable

 real(sdsu_fpd) :: bscat            ( max_ens,max_elv ) ! backscatter coeff
 real(sdsu_fpd) :: vt_sum           ( max_ens,max_elv ) ! sum of vt
 real(sdsu_fpd) :: vt_wgt           ( max_ens,max_elv ) ! weight of vt  
 real(sdsu_fpd) :: vdop_full( max_ens,max_ens,max_elv ) ! doppler velocity [m/s]

 type( polradar_parameter ) :: polout  !this is defined in module_mueller_matrix.F

!
! p3 ice output
!
 real(sdsu_fpd) :: &
  qnorm, &       ! normalized mass (Qtot/Ntot) [kg]
  Dth, &         ! critical size separating spherical/nonspherical ice [m]
  Dgr , &        ! size of equal mass for graupel and unrimed ice [m]
  Dcr , &        ! size of equal mass for graupel and partially rimed ice [m]
  alpha_md(4), & ! alpha for m-D relationship ( mass = alpha * D ** beta )
  beta_md(4),  & ! beta for m-D relationship ( mass = alpha * D ** beta )
  alpha_ad(4), & ! alpha for a-D relationshi (aera = alpha * D ** beta)
  beta_ad(4),  & ! beta for a-D relationshi (aera = alpha * D ** beta)
  n0norm   ,  &  ! normalized intercept of gamma PSD [1/m**(4+mu)]
  lambda ,  &    ! slope of gamma PSD [1/m]
  mu     ,  &    ! shape param of gamma PSD [-]
  re             ! effective radius of p3 ice [micron]


real :: q_min

!
! Begin by checking if hydrometeors of this species are present.
! If not, set scattering parameters to zero and return.
!
  q_min = 261.7**(real(1+10)*0.1)*1.e-18 ! Qtot/Ntot [kg]

  if(qtot < q_min .or. ntot == 0. ) then
     bmueller_mat_out=0. ; fscat_mat_out=cmplx(0.,0.) ; vdop_out=0.
     return
  endif


!
! Assign some useful constants
!
  wave=const_ghz_mm/freq  !wavelength [mm]

!
! If hydrometeors are present, initialize the scattering parameters
!
 mueller_mat = 0. ; mueller_mat_sum = 0. ; mueller_mat_full = 0.
 fscat_mat  = cmplx(0.,0.) ; fscat_mat_sum  = cmplx(0.,0.) ; fscat_mat_full  = cmplx(0.,0.)
 vdop_full = 0. ;  vt_sum = 0.  ;  vt_wgt = 0.   


!
! size increment: power of 10 base
!     D = 10.**(real(i)*coef)  ![mm]
!
 imax = 22 !1000  ! maximum size loop index
 d_max = 20. ! max diameter [mm]  --> 2cm
 coef = 0.5759 !log10(d_max)/real(imax)  !coefficient for power law

 !
 ! derive P3 PSD parameters and m-D and a-D relationship
 !
 qnorm  = dble( qtot*1.e-3/ntot ) ! [kg]
 call p3_psd_init( qnorm, dble(frim), dble(rhorim), Dth, Dgr, Dcr, &
                   alpha_md, beta_md, alpha_ad, beta_ad, n0norm, lambda, mu , re)
!print*, 'Dth, Dgr, Dcr',Dth, Dgr, Dcr


!
! loop for axix ratio
!
 loop_axis: do naxis = 1, 1  !toshii num_axis

!
! size loop
!

   loop_size: do i=1,imax

     d_mm = 0.5572*exp(coef*real(i)) * 1.e-3
     d_increment = ( 0.5572*exp((real(i)+0.5)*coef) - .5572*exp((real(i)-0.5)*coef) ) * 1.e-3

!     d_mm = 10.**(real(i)*coef)  !diameter [mm]
!     d_increment = 10.**((real(i)+0.5)*coef) - 10.**((real(i)-0.5)*coef)  ! dD in [mm]

     d_m  = d_mm * 1e-3                     ! diameter [m]
     xd  = dble(const_pi*d_mm/wave)         ! size parameter [-]



!
! Get N(D), density, area, vt
!         [1/m**(4+mu)]  [ m**mu ]      [1/m]  [m]   = [1/m4]
!     num = n0           * d_m**mu * exp(-lam * d_m )
     call p3_psd_param( dble(d_m), Dth, Dgr, Dcr, alpha_md, beta_md, alpha_ad, beta_ad, n0norm*ntot, lambda, mu ,&
                        num, rhoe, area, vt , p3_ice_regime)

!print*, d_m*1.e6, num, rhoe, vt , p3_ice_regime

   !
   ! assign axis ratio and orientation angle distributions for different P3 ice regime
   !
   select case(p3_ice_regime)
   case(1)  !small spherical ice
     axis_ratio = 1. ! sphere
     ityp_dist  = 0  !ramdom orientation 
     mean_orient= 0. ! [deg]
     num_stdev  = 1  !
     stdev_orient(1) = 0. ! [deg] does not account for random

   case(2)  ! dendrite
     axis_ratio = axis_ratio_gen(naxis)%ice
     ityp_dist  = ityp_dist_gen%ice
     mean_orient= mean_orient_gen%ice
     num_stdev  = num_stdev_gen%ice
     stdev_orient(1:num_stdev) = stdev_orient_gen(1:num_stdev)%ice

   case(3)  ! graupel
     axis_ratio = axis_ratio_gen(naxis)%graupel
     ityp_dist  = ityp_dist_gen%graupel
     mean_orient= mean_orient_gen%graupel
     num_stdev  = num_stdev_gen%graupel
     stdev_orient(1:num_stdev) = stdev_orient_gen(1:num_stdev)%graupel

   case(4)  ! partially rimed crystals (snow-grapuel)
     axis_ratio = axis_ratio_gen(naxis)%snow
     ityp_dist  = ityp_dist_gen%snow
     mean_orient= mean_orient_gen%snow
     num_stdev  = num_stdev_gen%snow
     stdev_orient(1:num_stdev) = stdev_orient_gen(1:num_stdev)%snow

   end select 

     !
     ! derive effective permittivity (liquid, solid, mixed phases)
     !
     call  effective_permittivity(freq, 0., rhoe, temp, crefd, &
                                  creal,cimag )

     !
     ! derive size-dependent axis ratio, if specified function id.
     ! 
     !call spectrum_axis_ratio( trim(spc) , ityp_axis(naxis), d_mm*0.1, axis_temp )
     daxis_ratio = dble(axis_ratio)

     !
     ! drive tmatrix :primal output is tmat_out(20,40,40)
     !
     diam_sph = dble(d_mm) * 0.1d0 ! spherical diaemter [cm]
     dwave  = dble(wave) * 0.1d0   ! wavelength [cm]

     call driver_tmatrix (daxis_ratio, dble(creal), dble(cimag), diam_sph, dwave, &   ! input
                          tmat_out, nrank, nm, conk, aovrb, sigma_temp )  ! output


     !
     ! drive Mueller matrix to derive 4x4 mueller matrix and 2x2 forward scattering matrix
     ! We also need extinciton and backscatter from this.....
     !
     call driver_mueller( daxis_ratio, dble(creal), dble(cimag), diam_sph, dwave, &   ! input
                          ityp_dist, dble(mean_orient), num_stdev, dble(stdev_orient) , &
                          tmat_out, conk, sigma_temp, num_elev, dble(elev(1:num_elev)) , &
                          mueller_mat, fscat_mat  )

      !
      !debbug toshii (testing only)
      !
      !call derive_output_polradar( mueller_mat(1,1,1:4,1:4), fscat_mat(1,1,1:2,1:2), dble(wave*1.e3), polout )
      !print*,daxis_ratio,creal,cimag,diam_sph,dwave
      !print*,'zhh, zdr, kdp =', polout%zhh, polout%zdr, polout%kdp


     !
     ! derive backscattering component 
     !
     do nstdev = 1, num_stdev ; do melev = 1, num_elev
        call bscat_mueller( mueller_mat ( nstdev, melev, 1:4,1:4 ), bscat(nstdev, melev) )
     enddo ; enddo


!
! integrate over particle size distribution; (Cross sectin extionction)*(cross sectional aerea)
!
!    [1/m3] =[1/m4]*[mm]        *[m/mm] --> [1/m3] 
     factor =  num * d_increment * 1.e-3

     mueller_mat_sum = mueller_mat_sum + mueller_mat * factor !accumulating 4x4 mueller matrix

     fscat_mat_sum  = fscat_mat_sum  + fscat_mat  * factor !accumulating 2x2 forward scattering matrix

     do nstdev = 1, num_stdev ; do melev = 1, num_elev
        vt_sum( nstdev, melev ) = vt_sum( nstdev, melev ) +  vt * factor * bscat(nstdev, melev)
        vt_wgt( nstdev, melev ) = vt_wgt( nstdev, melev ) +       factor * bscat(nstdev, melev)
!        print*, 'MSG ss_gen_micro:', i, melev, factor, bscat(nstdev, melev)
     enddo ; enddo

   enddo loop_size

!
! for each axis ratio loop
!
   mueller_mat_full(naxis,1:num_stdev,1:num_elev,1:4,1:4 ) = &
             mueller_mat_sum       (1:num_stdev,1:num_elev,1:4,1:4 ) !output 4x4 mueller matrix


   fscat_mat_full(naxis,1:num_stdev,1:num_elev,1:2,1:2 )  = &
             fscat_mat_sum       (1:num_stdev,1:num_elev,1:2,1:2 ) !output 2x2 forward scattering matrix

   do nstdev = 1, num_stdev ; do melev = 1, num_elev 
      if( vt_wgt( nstdev, melev ) <= 0.0 ) then
          vdop_full(naxis,nstdev, melev) = 0.
      else    
          vdop_full(naxis,nstdev, melev) = vt_sum( nstdev, melev ) / vt_wgt( nstdev, melev ) !doppler velocity [m/s] 
      endif
   enddo ; enddo

 enddo loop_axis

!==================================================================================
! toshii --> eventually we must derive ensemble mean, max, minima component or stdev
! right now, just 1st memmber of ensemble
!

 do melev = 1, num_elev

   !
   !forward and backscatter component of mueller matrix (so far no ensembles)
   !
   bmueller_mat_out(melev,1:4,1:4) =  mueller_mat_full(1,1,melev,  1:4,1:4 ) 

   !
   ! forward scattering matrix output 
   !
   fscat_mat_out (melev,1:2,1:2) =  fscat_mat_full(1,1,melev,     1:2,1:2 )  

   !
   ! doppler velocity 
   !
   vdop_out(melev) = vdop_full(1,1,melev) 

   !
   ! derive polarimetric parameters (testing only)
   !
   ! call derive_output_polradar( bmueller_mat_out(1:4,1:4), fscat_mat_out (1:2,1:2), wave*0.001d0, &
   !                              output  )

 enddo

 return
 end subroutine ss_p3_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine ss_rams_polarris(spc,freq, temp, wc, re, gnu, cfvt, pwvt, &
                          num_elev, elev, num_axis, ityp_axis, axis_ratio, ityp_dist, &
                          mean_orient, num_stdev, stdev_orient, &
                          cfmas, pwmas, frac_liq, &
                          bmueller_mat_out, fscat_mat_out, vdop_out)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Integrate over size distributions to derive total Mueller and forward scattering matrix for a given
!  size distributions, effective density, and axis ratio of particles of general RAMS microphysics.
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 character*(*),intent(in) :: spc !species character
 real(sdsu_fps),intent(in) :: freq  ! frequency of radiation [GHz]
 real(sdsu_fps),intent(in) :: temp  ! temperature of particles [K]
 real(sdsu_fps),intent(in) :: wc    ! water content of condensate [g/m3]
 real(sdsu_fps),intent(in) :: re    ! drop effective radius [micron] 
 real(sdsu_fps),intent(in) :: gnu   ! PSD shape parameter for generalized gamma distribution
 real(sdsu_fps),intent(in) :: cfvt   ! alpha_vt 
 real(sdsu_fps),intent(in) :: pwvt   ! beta_vt  [-]

 integer,intent(in) :: num_elev   ! number of elevation angle
 real(sdsu_fps),intent(in) :: elev(num_elev)  ! particle axis ratio [-]
 integer,intent(in) :: num_axis   ! number of axis ratio 
 integer,intent(in) :: ityp_axis(num_axis)  ! particle axis ratio model type [-] 
 real(sdsu_fps),intent(in) :: axis_ratio(num_axis)  ! particle axis ratio [-]
 integer,intent(in) :: ityp_dist  ! particle orienting angle distribution type
 real(sdsu_fps),intent(in) :: mean_orient   ! mean orienting angle [deg]
 integer,intent(in) :: num_stdev  ! number of orienting angle standard deviations
 real(sdsu_fps),intent(in) :: stdev_orient(num_stdev)  ! standard deviation / theta_m of orienting angle [deg]

 real(sdsu_fps),intent(in) :: cfmas  ! alpha_m [kg /m**beta_m]
 real(sdsu_fps),intent(in) :: pwmas  ! beta_m  [-]
 real(sdsu_fps),intent(in) :: frac_liq ! fracion of liquid (1=rain or 0=ice condensates, 0~1=melting snow)
 real(sdsu_fpd),intent(out) :: vdop_out(num_elev)  ! Doppler velocity [m/s]

!--------Local Parameters
 real(sdsu_fps) :: rhoe ! bulk effective density of hydrometero [kg/m**3]

 integer :: i    !looping 
 integer :: imax !max for looping

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
 integer,save :: nshp !SCATDB shape #
 real(sdsu_fps) :: qext_ns, qsca_ns, asym_ns, qbsca_ns ! non-spherical single scattering
 real(sdsu_fps) :: d_mm_solid  !solid ice diameter [mm]
 logical :: success ! status 
 real(sdsu_fps) :: rad_solid_micron  !solid ice radius [micron]
 real(sdsu_fps) :: area   !area of fully particle [mm2]

! tmatrix
 integer :: nrank, nm
 real(sdsu_fps) :: creal, cimag  !real / imageneary component of dielectric function
 real(sdsu_fpd) :: diam_sph        !spherical diaemter [cm] 
 real(sdsu_fpd) :: dwave           !wavelength of radar band [cm]
 real(sdsu_fpd) :: daxis_ratio, conk, aovrb, sigma_temp
 complex*8 :: tmat_out(20,40,40)   ! maybe dynamic allocation
 real(sdsu_fps) :: axis
 integer :: naxis
 real(sdsu_fps) :: axis_temp

!
! mueller matrix
!
 integer :: nstdev , melev
 integer,parameter :: nforward = 1 , nback = 5  !index of forward and backscatter angle in mueller matrix
 real(sdsu_fpd)    :: mueller_mat          (max_ens,max_elv,4,4) ! 4x4 mueller matrix  (backscatter component)
 real(sdsu_fpd) :: mueller_mat_sum         (max_ens,max_elv,4,4) ! summation of 4x4 mueller matrix 
 real(sdsu_fpd) :: mueller_mat_full(max_ens,max_ens,max_elv,4,4) ! full bundle of 4x4 mueller matrix 
 complex :: fscat_mat              (max_ens,max_elv,2,2) ! 2x2 forward scattering matrix 
 complex :: fscat_mat_sum          (max_ens,max_elv,2,2) ! summation of 2x2 forward scattering matrix 
 complex :: fscat_mat_full (max_ens,max_ens,max_elv,2,2) ! full bundle of 2x2 forward scattering matrix 
 real(sdsu_fpd),intent(out) :: bmueller_mat_out(num_elev,4,4) ! 4x4 mueller matrix (forward and backscatering component only)
 complex,intent(out) :: fscat_mat_out (num_elev,2,2)    ! 2x2 forward scattering matrix 

 real(sdsu_fpd) :: bscat            ( max_ens,max_elv ) ! backscatter coeff
 real(sdsu_fpd) :: vt_sum           ( max_ens,max_elv ) ! sum of vt
 real(sdsu_fpd) :: vt_wgt           ( max_ens,max_elv ) ! weight of vt  
 real(sdsu_fpd) :: vdop_full( max_ens,max_ens,max_elv ) ! doppler velocity [m/s]


!
! Begin by checking if hydrometeors of this species have enough water content.
! If not, set scattering parameters to zero and return.
!

 if(wc .lt. q_min_micro) then
    bmueller_mat_out=0. ; fscat_mat_out=cmplx(0.,0.) ; vdop_out=0.
    return
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
 ntot = wc * 1e-3 / mean_mass                  !total particle number concentration [#/m3]


!
! initialize the optical parameters
!
 mueller_mat = 0. ; mueller_mat_sum = 0. ; mueller_mat_full = 0.
 fscat_mat = cmplx(0.,0.) ; fscat_mat_sum = cmplx(0.,0.) ; fscat_mat_full = cmplx(0.,0.)
 vdop_full = 0. ;  vt_sum = 0.  ;  vt_wgt = 0.


 if( trim(spc) == 'cloud1' .or. trim(spc(1:3)) == 'ice' )  then  !non-precip

! define imax for a givne radius increment 20. is maximum diameter
  d_increment = 0.01  ! [mm] --> 1micron
  d_min =  0.002  ![mm]  --> 2micron
  d_max = 0.2  ! [mm]  --> 200micron

 else  !precip species

! define imax for a givne radius increment 20. is maximum diameter
  d_increment = 1.0  ! [mm]  --> 100micron
  d_min =  0.2  ![mm]  --> 200micron
  d_max = 20.  ! [mm]  --> 2cm

 endif

 imax = nint( (d_max-d_min) /d_increment)

!
! loop for axix ratio
!
 loop_axis: do naxis = 1, num_axis
   axis_temp = axis_ratio(naxis)

!
! loop over particle size
!
   size_loop: do i=0,imax

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
      call  effective_permittivity(freq, fmelt, rhoe, temp, crefd,  creal, cimag)

     !     
     ! derive size-dependent axis ratio, if specified function id.
     !     
     call spectrum_axis_ratio( trim(spc) , ityp_axis(naxis), d_mm*0.1, axis_temp )
     daxis_ratio = dble(axis_temp)


     !        
     ! drive tmatrix :primal output is tmat_out
     !
     diam_sph = dble(d_mm) * 0.1d0 ! spherical diaemter [cm]
     dwave  = dble(wave) * 0.1d0   ! wavelength [cm]
     call driver_tmatrix(daxis_ratio, dble(creal), dble(cimag), diam_sph, dwave, &              ! input
                         tmat_out, nrank, nm, conk, aovrb, sigma_temp )  ! output


     !
     ! drive Mueller matrix to derive 4x4 mueller matrix and 2x2 forward scattering matrix
     ! We also need extinciton and backscatter from this.....
     !

     call driver_mueller( daxis_ratio, dble(creal), dble(cimag), diam_sph, dwave, &   ! input
                          ityp_dist, dble(mean_orient), num_stdev, dble(stdev_orient) , &
                          tmat_out, conk, sigma_temp, num_elev, dble(elev(1:num_elev)) , &
                          mueller_mat, fscat_mat  )



     !
     ! derive backscattering component 
     !
     do nstdev = 1, num_stdev ; do melev = 1, num_elev
        call bscat_mueller( mueller_mat ( nstdev, melev, 1:4,1:4 ), bscat(nstdev, melev) )
     enddo ; enddo


     !
     ! terminal velocity [m/s]
     !
      vt = cfvt * (d_m**pwvt)  ![m/s]

!
! integrate over particle size distribution; (Cross sectin extionction)*(cross sectional aerea)
!
!    [1/m3] =[1/m4]*[mm]        *[m/mm] --> [1/m3] 
     factor =  num * d_increment * 1.e-3

     mueller_mat_sum = mueller_mat_sum + mueller_mat * factor !accumulating 4x4 mueller matrix

     fscat_mat_sum  = fscat_mat_sum  + fscat_mat  * factor !accumulating 2x2 forward scattering matrix

     do nstdev = 1, num_stdev ; do melev = 1, num_elev
        vt_sum( nstdev, melev ) = vt_sum( nstdev, melev ) +  vt * factor * bscat(nstdev, melev)
        vt_wgt( nstdev, melev ) = vt_wgt( nstdev, melev ) +       factor * bscat(nstdev, melev)
     enddo ; enddo

  enddo size_loop 

  !
  ! output for each axis ratio
  !
   mueller_mat_full(naxis,1:num_stdev,1:num_elev,1:4,1:4 ) = &
             mueller_mat_sum       (1:num_stdev,1:num_elev,1:4,1:4 ) !output 4x4 mueller matrix


   fscat_mat_full(naxis,1:num_stdev,1:num_elev,     1:2,1:2 )  = &
             fscat_mat_sum       (1:num_stdev,1:num_elev,     1:2,1:2 ) !output 2x2 forward scattering matrix

   do nstdev = 1, num_stdev ; do melev = 1, num_elev
      if( vt_wgt( nstdev, melev ) <= 0.0 ) then
          vdop_full(naxis,nstdev, melev) = 0.
      else
          vdop_full(naxis,nstdev, melev) = vt_sum( nstdev, melev ) / vt_wgt( nstdev, melev ) !doppler velocity [m/s] 
      endif
   enddo ; enddo


 enddo loop_axis

 do melev = 1, num_elev

   !
   !forward and backscatter component of mueller matrix (so far no ensembles)
   !
   bmueller_mat_out(melev,1:4,1:4) =  mueller_mat_full(1,1,melev,  1:4,1:4 )

   !
   ! forward scattering matrix output 
   !
   fscat_mat_out (melev,1:2,1:2) =  fscat_mat_full(1,1,melev,     1:2,1:2 )

   !
   ! doppler velocity 
   !
   vdop_out(melev) = vdop_full(1,1,melev)

   !
   ! derive polarimetric parameters (testing only)
   !
   ! call derive_output_polradar( bmueller_mat_out(1:4,1:4), fscat_mat_out (1:2,1:2), wave*0.001d0, &
   !                              output  )

 enddo

 return
 end subroutine ss_rams_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine ss_sbm_polarris(spc, sbin, ebin, freq, temp, bulkden, bulkrad, dbulkrad, vt, wc, n0, fmelt_bin,frime_bin,&
                            num_elev, elev, num_axis, ityp_axis, axis_ratio, ityp_dist,  &
                            mean_orient, num_stdev, stdev_orient, &
                            bmueller_mat_out, fscat_mat_out, vdop_out )
  implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Integrate over size distributions to derive total Mueller and forward scattering matrix for a given
!  size distributions, effective density, and axis ratio of particles of HUCM SBM. 
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in):: spc      ! index name of hydrometeor classes 
 integer,intent(in)         :: sbin,ebin! starting and ending bin for size integration
 real(sdsu_fps), intent(in) :: freq     ! frequency of radiation [GHz]
 real(sdsu_fps), intent(in) :: temp     ! temperature of particles [K]
 real(sdsu_fps), intent(in) :: bulkden(1:nbin) ! bulk density of hydrometero for each bins [g/cm3]
 real(sdsu_fps), intent(in) :: bulkrad(1:nbin) ! bulk radius of hydrometero for each bins [cm] 
 real(sdsu_fps), intent(in) :: dbulkrad(1:nbin)! binsize [cm]
 real(sdsu_fps), intent(in) :: vt(1:nbin)      ! terminal velocity [cm/s]
 real(sdsu_fps), intent(in) :: wc        ! water content of cloud ice distribution [g/m**3]
 real(sdsu_fps), intent(in) :: n0(sbin:ebin)        ! number concentration of hydrometeors [1/m4]
 real(sdsu_fps), intent(in) :: fmelt_bin(sbin:ebin) !mass melting fraction of hydrometoer for each bins [-] 
 real(sdsu_fps), intent(in) :: frime_bin(sbin:ebin) !mass rime fraction of snow only for each bin [-] 

 integer,intent(in) :: num_elev   ! number of elevation angle
 real(sdsu_fps),intent(in) :: elev(num_elev)  ! particle axis ratio [-]
 integer,intent(in) :: num_axis   ! number of axis ratio 
 integer,intent(in) :: ityp_axis(num_axis)  ! particle axis ratio model type [-] 
 real(sdsu_fps),intent(in) :: axis_ratio(num_axis)  ! particle axis ratio [-]
 integer,intent(in) :: ityp_dist  ! particle orienting angle distribution type
 real(sdsu_fps),intent(in) :: mean_orient   ! mean orienting angle [deg]
 integer,intent(in) :: num_stdev  ! number of orienting angle standard deviations
 real(sdsu_fps),intent(in) :: stdev_orient(num_stdev)  ! standard deviation / theta_m of orienting angle [deg]

 real(sdsu_fpd), intent(out) :: vdop_out(num_elev)   ! Doppler velocity [m/s]

!----Local variables
 real(sdsu_fps) :: morient   ! mean orienting angle [deg]

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

 real(sdsu_fps) :: fv_melt
 integer :: melt_temp 

! tmatrix
 integer :: nrank, nm
 real(sdsu_fps) :: creal, cimag  !real / imageneary component of dielectric function
 real(sdsu_fpd) :: diam_sph        !spherical diaemter [cm] 
 real(sdsu_fpd) :: dwave           !wavelength of radar band [cm]
 real(sdsu_fpd) :: daxis_ratio, conk, aovrb, sigma_temp
 complex*8 :: tmat_out(20,40,40)   ! maybe dynamic allocation
 real(sdsu_fps) :: axis
 integer :: naxis
 real(sdsu_fps) :: axis_temp

!
! mueller matrix
!
 integer :: nstdev , melev
 integer,parameter :: nforward = 1 , nback = 5  !index of forward and backscatter angle in mueller matrix
 real(sdsu_fpd)    :: mueller_mat          (max_ens,max_elv,4,4) ! 4x4 mueller matrix  (backscatter component)
 real(sdsu_fpd) :: mueller_mat_sum         (max_ens,max_elv,4,4) ! summation of 4x4 mueller matrix 
 real(sdsu_fpd) :: mueller_mat_full(max_ens,max_ens,max_elv,4,4) ! full bundle of 4x4 mueller matrix 
 complex :: fscat_mat              (max_ens,max_elv,2,2) ! 2x2 forward scattering matrix 
 complex :: fscat_mat_sum          (max_ens,max_elv,2,2) ! summation of 2x2 forward scattering matrix 
 complex :: fscat_mat_full (max_ens,max_ens,max_elv,2,2) ! full bundle of 2x2 forward scattering matrix 
 real(sdsu_fpd),intent(out) :: bmueller_mat_out(num_elev,4,4) ! 4x4 mueller matrix (backscatering component only)
 complex,intent(out) :: fscat_mat_out (num_elev,2,2)    ! 2x2 forward scattering matrix 
 type ( polradar_parameter ) :: output ! polarimeterica radar measureable

 real(sdsu_fpd) :: bscat            ( max_ens,max_elv ) ! backscatter coeff
 real(sdsu_fpd) :: vt_sum           ( max_ens,max_elv ) ! sum of vt
 real(sdsu_fpd) :: vt_wgt           ( max_ens,max_elv ) ! weight of vt  
 real(sdsu_fpd) :: vdop_full( max_ens,max_ens,max_elv ) ! doppler velocity [m/s]


 type( polradar_parameter ) :: polout  !this is defined in module_mueller_matrix.F

!
! Begin by checking if hydrometeors of this species are present.
! If not, set scattering parameters to zero and return.
!
  if(wc .lt. q_min_micro) then
     bmueller_mat_out=0. ; fscat_mat_out=cmplx(0.,0.) ; vdop_out=0.
     return
  endif


!
!     Assign some useful constants
 wave=const_ghz_mm/freq

!
!     If hydrometeors are present, initialize the scattering parameters
!
 mueller_mat = 0. ; mueller_mat_sum = 0. ; mueller_mat_full = 0.
 fscat_mat  = cmplx(0.,0.) ; fscat_mat_sum  = cmplx(0.,0.) ; fscat_mat_full  = cmplx(0.,0.)
 vdop_full = 0. ;  vt_sum = 0.  ;  vt_wgt = 0.


!
! loop for axix ratio
!
 loop_axis: do naxis = 1, num_axis
   axis_temp = axis_ratio(naxis)


!
!     Loop over particle sizes:
!
  loop_size: do i=sbin, ebin  !1,nbin

    if(n0(i) == 0. ) cycle  !if number density is zero, skip all process

    !
    ! estimate diagnostic radius and density based on melting and riming fraction
    !
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
      fmelt = fv_melt           ! melting fraction [-]

    else

      rad   = bulkrad(i)*10.e0  ! radius [mm]
      d_mm  = rad * 2.          ! diameter [mm]
      d_m   = d_mm * 1.e-3      ! diameter [m]
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
    call  effective_permittivity(freq, fmelt, rhoe, temp, crefd, creal, cimag)

    !     
    ! derive axis ratio
    !     
    call spectrum_axis_ratio( trim(spc) , ityp_axis(naxis), d_mm*0.1, axis_temp )
    daxis_ratio = dble(axis_temp)


    !        
    ! drive tmatrix :primal output is tmat_out
    !
!     print*, axis_ratio 
    diam_sph = dble(d_mm) * 0.1d0 ! spherical diaemter [cm]
    dwave  = dble(wave) * 0.1d0   ! wavelength [cm]
    call driver_tmatrix(daxis_ratio, dble(creal), dble(cimag), diam_sph, dwave, &              ! input
                         tmat_out, nrank, nm, conk, aovrb, sigma_temp )  ! output

!toshii debbugg
!    if( trim(spc) == 'hail')   then
!     if( d_mm < 28.) then
!        morient = 0.
!     else  !for hail stone
!        morient =  90.
!     endif
!    else
     morient = mean_orient
!    endif

    !
    ! drive Mueller matrix to derive 4x4 mueller matrix and 2x2 forward scattering matrix
    ! We also need extinciton and backscatter from this.....
    !
    call driver_mueller( daxis_ratio, dble(creal), dble(cimag), diam_sph, dwave, &   ! input
                          ityp_dist, dble(morient), num_stdev, dble(stdev_orient) , &
                          tmat_out, conk, sigma_temp, num_elev, dble(elev(1:num_elev)) , &
                          mueller_mat, fscat_mat  )



    !
    ! derive backscattering component 
    !
    do nstdev = 1, num_stdev ; do melev = 1, num_elev
       call bscat_mueller( mueller_mat ( nstdev, melev, 1:4,1:4 ), bscat(nstdev, melev) )
    enddo ; enddo



!
! integrate over particle size distribution;
!    [1/m3] =[1/m4]*[mm]  *[m/mm] --> [1/m3] 
     factor =  num * drad * 1.e-3

     mueller_mat_sum = mueller_mat_sum + mueller_mat * factor !accumulating 4x4 mueller matrix

     fscat_mat_sum  = fscat_mat_sum  + fscat_mat  * factor !accumulating 2x2 forward scattering matrix

      !
      !debbug toshii
      !
      !call derive_output_polradar( mueller_mat(1,1,1:4,1:4), fscat_mat(1,1,1:2,1:2), dble(wave*1.e3), polout )
      !print*,i,dwave,'d_mm=',d_mm,'zhh, zdr, kdp =', polout%zhh, polout%zdr, polout%kdp

    if( spc == 'snow   ')   then

      !
      ! diagnostic snow aggregate Vt as a function of riming and meling fracitons.
      !
      call vt_rime_melting( i, nbin, frime_bin(i), fmelt_bin(i), x_sbm%snow , vt_sbm%liq , &
                            vt_sbm%snow, vt_sbm%hail, den_sbm%snow, vt_snow ) 

     do nstdev = 1, num_stdev ; do melev = 1, num_elev
        vt_sum( nstdev, melev ) = vt_sum( nstdev, melev ) + 0.01*vt_snow * factor * bscat(nstdev, melev)
        vt_wgt( nstdev, melev ) = vt_wgt( nstdev, melev ) + factor * bscat(nstdev, melev)
     enddo ; enddo

    else !default (all other species)

     do nstdev = 1, num_stdev ; do melev = 1, num_elev
        vt_sum( nstdev, melev ) = vt_sum( nstdev, melev ) + 0.01*vt(i) * factor * bscat(nstdev, melev)
        vt_wgt( nstdev, melev ) = vt_wgt( nstdev, melev ) + factor * bscat(nstdev, melev)
     enddo ; enddo

    endif

  enddo loop_size !bin spectra loop

!
! for each axis ratio loop
!
   mueller_mat_full(naxis,1:num_stdev,1:num_elev,1:4,1:4 ) = &
             mueller_mat_sum       (1:num_stdev,1:num_elev,1:4,1:4 ) !output 4x4 mueller matrix


   fscat_mat_full(naxis,1:num_stdev,1:num_elev,     1:2,1:2 )  = &
             fscat_mat_sum       (1:num_stdev,1:num_elev,     1:2,1:2 ) !output 2x2 forward scattering matrix

   do nstdev = 1, num_stdev ; do melev = 1, num_elev
      if( vt_wgt( nstdev, melev ) <= 0.0 ) then
          vdop_full(naxis,nstdev, melev) = 0.
      else
          vdop_full(naxis,nstdev, melev) = vt_sum( nstdev, melev ) / vt_wgt( nstdev, melev ) !doppler velocity [m/s] 
      endif
   enddo ; enddo

 enddo loop_axis

!==================================================================================
! toshii --> eventually we must derive ensemble mean, max, minima component or stdev
! right now, just 1st memmber of ensemble
!

 do melev = 1, num_elev

   !
   !forward and backscatter component of mueller matrix (so far no ensembles)
   !
   bmueller_mat_out(melev,1:4,1:4) =  mueller_mat_full(1,1,melev,  1:4,1:4 )

   !
   ! forward scattering matrix output 
   !
   fscat_mat_out (melev,1:2,1:2) =  fscat_mat_full(1,1,melev,     1:2,1:2 )

   !
   ! doppler velocity 
   !
   vdop_out(melev) = vdop_full(1,1,melev)

   !
   ! derive polarimetric parameters (testing only)
   !
   ! call derive_output_polradar( bmueller_mat_out(1:4,1:4), fscat_mat_out (1:2,1:2), wave*0.001d0, &
   !                              output  )

 enddo

 return
 end subroutine ss_sbm_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine total_opt_gen(wave_m, press, w, atm_ext, bmueller_gen, fscat_gen, vdop_gen, &
                           opt )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Integrate over various hydrometeor to derive total Mueller and forward scattering matrix for 
!  general bulk microphyscs. 
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(in) :: wave_m  !wavelength [m]
 real(sdsu_fps),intent(in) :: press  !atmosphere pressure [mb]
 real(sdsu_fps),intent(in) :: w      !atmosphere vertical velocity [m/s]
 real(sdsu_fps),intent(in) :: atm_ext  !atmosphere extinction 
 type ( particle_gen_r8 ),intent(in) :: & ! particle_gen is defined in module_simulater
     bmueller_gen(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_gen                 ! Doppler velocity [m/s]
 type ( particle_gen_cmplx ),intent(in) :: & ! particle_gen is defined in module_simulater
     fscat_gen(2,2)             ! 2x2 forward scattering matrix 

 type(opt_polarris_parameter),intent(out) :: opt

! local
 real(sdsu_fpd) :: bmueller_tot(4,4)  !total mueller matrix
 real(sdsu_fpd) :: vdop_tot           !total doppler velocity [m/s]
 real(sdsu_fpd) :: bscat_tot          !total backscatter
 complex :: fscat_tot(2,2)            !total forward scattering 
 type ( particle_gen_r8 ) :: bscat ! backscatterig coeff 
 real(sdsu_fpd) :: vdop_raw
 type( polradar_parameter ) :: polout  !this is defined in module_mueller_matrix.F

!
! total mueller matrix 
!
 bmueller_tot = bmueller_gen%cloud  + &
                bmueller_gen%rain   + &
                bmueller_gen%ice    + &
                bmueller_gen%snow   + &
                bmueller_gen%graupel+ &
                bmueller_gen%hail 

!
! if all zero, return.
!
 if( maxval( bmueller_tot(:,:) ) == 0.d0 ) then
   opt%zhh   = undefined   ! reflectivity [dBZ] 
   opt%zdr   = undefined   ! differential reflectivity [dB]
   opt%kdp   = undefined   ! specific diff phase [deg/km] 
   opt%rhohv = undefined   ! copolar (hv) correlation coefficient [-]
   opt%vdop  = undefined   ! doppler velocity [m/s]
   return
 endif 


!
! total forward scattering matrix 
!
 fscat_tot    = fscat_gen%cloud  + &
                fscat_gen%rain   + &
                fscat_gen%ice    + &
                fscat_gen%snow   + &
                fscat_gen%graupel+ &
                fscat_gen%hail

!
! derive polarimetric observables
!
 call derive_output_polradar( bmueller_tot, fscat_tot, dble(wave_m), polout  )

 opt%zhh     = real( polout%zhh )   ! reflectivity [dBZ] 
 opt%zdr     = real( polout%zdr )   ! differential reflectivity [dB]
 opt%kdp     = real( polout%kdp )   ! specific diff phase [deg/km] 
 opt%rhohv   = real( polout%rhohv ) ! copolar (hv) correlation coefficient [-]

!
! derive bakscatter
!
 call bscat_mueller( bmueller_gen%cloud, bscat%cloud )
 call bscat_mueller( bmueller_gen%rain , bscat%rain  )
 call bscat_mueller( bmueller_gen%ice  , bscat%ice   )
 call bscat_mueller( bmueller_gen%snow , bscat%snow  )
 call bscat_mueller( bmueller_gen%graupel, bscat%graupel )
 call bscat_mueller( bmueller_gen%hail , bscat%hail  )

 bscat_tot = bscat%cloud + bscat%rain + bscat%ice + bscat%snow + bscat%graupel + bscat%hail

!
! Doppler velocity of hydrometeors [m/s]
!
 vdop_raw = &
                ( bscat%cloud * vdop_gen%cloud  + &
                  bscat%rain  * vdop_gen%rain   + &
                  bscat%ice   * vdop_gen%ice    + &
                  bscat%snow  * vdop_gen%snow   + &
                  bscat%graupel * vdop_gen%graupel  + &
                  bscat%hail  * vdop_gen%hail ) /  bscat_tot

!
! correct it with pressure and actual atmospheric w (downward is positive)
!
  opt%vdop  = real(vdop_raw)*sqrt(800.e0/press)  - w ! actual terminal velocity for a given pressure [m/s]


 return
 end subroutine total_opt_gen

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine total_opt_sbm( wave_m, press, w, atm_ext, bmueller_sbm, fscat_sbm,  vdop_sbm , &
                           opt )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Integrate over various hydrometeor to derive total Mueller and forward scattering matrix for 
!  HUCM SBM.
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: wave_m  !wavelength [m]
 real(sdsu_fps),intent(in) :: press  !atmosphere pressure [mb]
 real(sdsu_fps),intent(in) :: w      !atmosphere vertical velocity [m/s]
 real(sdsu_fps),intent(in) :: atm_ext  !atmosphere extinction 
 type ( particle_sbm_r8 ),intent(in) :: & ! particle_sbm is defined in module_simulater
     bmueller_sbm(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_sbm                 ! Doppler velocity [m/s]
 type ( particle_sbm_cmplx ),intent(in) :: & ! particle_sbm is defined in module_simulater
     fscat_sbm(2,2)             ! 2x2 forward scattering matrix 
 type(opt_polarris_parameter),intent(out) :: opt

! local
 real(sdsu_fpd) :: bmueller_tot(4,4)  !total mueller matrix
 real(sdsu_fpd) :: vdop_tot           !total doppler velocity [m/s]
 real(sdsu_fpd) :: bscat_tot          !total backscatter
 complex :: fscat_tot(2,2)            !total forward scattering 
 type ( particle_sbm_r8 ) :: bscat ! backscatterig coeff 
 real(sdsu_fpd) :: vdop_raw
 type( polradar_parameter ) :: polout  !this is defined in module_mueller_matrix.F

!
! total mueller matrix 
!
 bmueller_tot = bmueller_sbm%liq     + &
                bmueller_sbm%ice_col + &
                bmueller_sbm%ice_pla + &
                bmueller_sbm%ice_den + &
                bmueller_sbm%snow    + &
                bmueller_sbm%graupel + &
                bmueller_sbm%hail

!
! if all zero, return.
!
 if( maxval( bmueller_tot(:,:) ) == 0.d0 ) then
   opt%zhh   = undefined   ! reflectivity [dBZ] 
   opt%zdr   = undefined   ! differential reflectivity [dB]
   opt%kdp   = undefined   ! specific diff phase [deg/km] 
   opt%rhohv = undefined   ! copolar (hv) correlation coefficient [-]
   opt%vdop  = undefined   ! doppler velocity [m/s]
   return
 endif


!
! total forward scattering matrix 
!
 fscat_tot    = fscat_sbm%liq     + &
                fscat_sbm%ice_col + &
                fscat_sbm%ice_pla + &
                fscat_sbm%ice_den + &
                fscat_sbm%snow    + &
                fscat_sbm%graupel + &
                fscat_sbm%hail

!
! derive polarimetric observables
!
 call derive_output_polradar( bmueller_tot, fscat_tot, dble(wave_m), polout  )

 opt%zhh     = real( polout%zhh )   ! reflectivity [dBZ] 
 opt%zdr     = real( polout%zdr )   ! differential reflectivity [dB]
 opt%kdp     = real( polout%kdp )   ! specific diff phase [deg/km] 
 opt%rhohv   = real( polout%rhohv ) ! copolar (hv) correlation coefficient [-]

!
! derive bakscatter
!
 call bscat_mueller( bmueller_sbm%liq     , bscat%liq )
 call bscat_mueller( bmueller_sbm%ice_col , bscat%ice_col  )
 call bscat_mueller( bmueller_sbm%ice_pla , bscat%ice_pla  )
 call bscat_mueller( bmueller_sbm%ice_den , bscat%ice_den  )
 call bscat_mueller( bmueller_sbm%snow    , bscat%snow  )
 call bscat_mueller( bmueller_sbm%graupel , bscat%graupel )
 call bscat_mueller( bmueller_sbm%hail    , bscat%hail  )

 bscat_tot = bscat%liq  + bscat%ice_col + bscat%ice_pla + bscat%ice_den &
           + bscat%snow + bscat%graupel + bscat%hail

!
! Doppler velocity of hydrometeors [m/s]
!
 vdop_raw = &
                ( bscat%liq     * vdop_sbm%liq     + &
                  bscat%ice_col * vdop_sbm%ice_col + &
                  bscat%ice_pla * vdop_sbm%ice_pla + &
                  bscat%ice_den * vdop_sbm%ice_den + &
                  bscat%snow    * vdop_sbm%snow   + &
                  bscat%graupel * vdop_sbm%graupel  + &
                  bscat%hail    * vdop_sbm%hail ) /  bscat_tot


!
! correct it with pressure and actual atmospheric w (downward is positive)
! (w: upward is positive)
  opt%vdop  = real(vdop_raw) *sqrt(1000.e0/press) - w ! measured Doppler velocity [m/s]


 return
 end subroutine total_opt_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine total_opt_rams(wave_m, rho_dair, w, atm_ext, bmueller_rams, fscat_rams, vdop_rams, &
                           opt )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Integrate over various hydrometeor to derive total Mueller and forward scattering matrix for 
!  general RAMS microphysics. 
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(in) :: wave_m  !wavelength [m]
 real(sdsu_fps),intent(in) :: rho_dair  !atmosphere dry air density  [kg/m3]
 real(sdsu_fps),intent(in) :: w  !atmosphere vertical velocity [m/s]
 real(sdsu_fps),intent(in) :: atm_ext  !atmosphere extinction 
 type ( particle_rams_r8 ),intent(in) :: & ! particle_rams is defined in module_simulater
     bmueller_rams(4,4),  & ! 4x4 forward and backward component of 4x4 mueller matrix 
     vdop_rams                 ! Doppler velocity [m/s]
 type ( particle_rams_cmplx ),intent(in) :: & ! particle_rams is defined in module_simulater
     fscat_rams(2,2)             ! 2x2 forward scattering matrix 

 type(opt_polarris_parameter),intent(out) :: opt

! local
 real(sdsu_fpd) :: bmueller_tot(4,4)  !total mueller matrix
 real(sdsu_fpd) :: vdop_tot           !total doppler velocity [m/s]
 real(sdsu_fpd) :: bscat_tot          !total backscatter
 complex :: fscat_tot(2,2)            !total forward scattering 
 type ( particle_rams_r8 ) :: bscat ! backscatterig coeff 
 real(sdsu_fpd) :: vdop_raw
 type( polradar_parameter ) :: polout  !this is defined in module_mueller_matrix.F


!
! total mueller matrix 
!
 bmueller_tot = bmueller_rams%cloud1 + &
                bmueller_rams%cloud2 + &
                bmueller_rams%rain   + &
                bmueller_rams%ice1   + &
                bmueller_rams%ice2   + &
                bmueller_rams%snow   + &
                bmueller_rams%graupel+ &
                bmueller_rams%hail

!
! if all zero, return.
!
 if( maxval( bmueller_tot(:,:) ) == 0.d0 ) then
   opt%zhh   = undefined   ! reflectivity [dBZ] 
   opt%zdr   = undefined   ! differential reflectivity [dB]
   opt%kdp   = undefined   ! specific diff phase [deg/km] 
   opt%rhohv = undefined   ! copolar (hv) correlation coefficient [-]
   opt%vdop  = undefined   ! doppler velocity [m/s]
   return
 endif


!
! total forward scattering matrix 
!
 fscat_tot = fscat_rams%cloud1 + &
             fscat_rams%cloud2 + &
             fscat_rams%rain   + &
             fscat_rams%ice1   + &
             fscat_rams%ice2   + &
             fscat_rams%snow   + &
             fscat_rams%graupel+ &
             fscat_rams%hail

!
! derive polarimetric observables
!
 call derive_output_polradar( bmueller_tot, fscat_tot, dble(wave_m), polout  )

 opt%zhh     = real( polout%zhh )   ! reflectivity [dBZ] 
 opt%zdr     = real( polout%zdr )   ! differential reflectivity [dB]
 opt%kdp     = real( polout%kdp )   ! specific diff phase [deg/km] 
 opt%rhohv   = real( polout%rhohv ) ! copolar (hv) correlation coefficient [-]

!
! derive bakscatter
!
 call bscat_mueller( bmueller_rams%cloud1, bscat%cloud1 )
 call bscat_mueller( bmueller_rams%cloud2, bscat%cloud2 )
 call bscat_mueller( bmueller_rams%rain  , bscat%rain   )
 call bscat_mueller( bmueller_rams%ice1  , bscat%ice1   )
 call bscat_mueller( bmueller_rams%ice2  , bscat%ice2   )
 call bscat_mueller( bmueller_rams%snow  , bscat%snow   )
 call bscat_mueller( bmueller_rams%graupel, bscat%graupel )
 call bscat_mueller( bmueller_rams%hail  , bscat%hail   )

 bscat_tot = bscat%cloud1 + bscat%cloud2 + bscat%rain + bscat%ice1 + &
             bscat%ice2   + bscat%snow   + bscat%graupel + bscat%hail

!
! Doppler velocity of hydrometeors [m/s]
!
 vdop_raw = &
                ( bscat%cloud1 * vdop_rams%cloud1  + &
                  bscat%cloud2 * vdop_rams%cloud2  + &
                  bscat%rain   * vdop_rams%rain   + &
                  bscat%ice1   * vdop_rams%ice1    + &
                  bscat%ice2   * vdop_rams%ice2    + &
                  bscat%snow   * vdop_rams%snow   + &
                  bscat%graupel * vdop_rams%graupel  + &
                  bscat%hail  * vdop_rams%hail ) /  bscat_tot


!
! correct it with pressure and actual atmospheric w (downward is positive)
!
  opt%vdop  = real(vdop_raw)*sqrt(0.7/rho_dair) - w  ! actual terminal velocity for a given air density[m/s]


 return
 end subroutine total_opt_rams

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine spectrum_axis_ratio( spc , ityp, d_cm , axis_ratio )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Derive axis ratio for a given particle size only if specified the tag ID. 
!
! axis ratio > 1.0 : prolate, like Almondo (hail, crystals)
! axis ratio < 1.0 : oblate , like mochi   (snow, rain)
!
! History:
! 05/2017  Toshi Matsui@NASA GSFC ; Add Ryzhkov2011 option for graupel/hail
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in):: spc     ! name of hydrometeor classes 
 integer,intent(in) :: ityp            ! axis ratio model type ID
 real, intent(in) :: d_cm              ! particle diameter [cm]
 real, intent(inout) :: axis_ratio     ! axix ratio [-]

 real,parameter :: cloud_thresh_cm = 0.02  !threshold to seperate cloud and precip [cm]

 real,parameter :: min_axis_rain = 0.4, max_axis_rain = 0.999

 real :: d_mm

 d_mm = 10.* d_cm

 if( ityp == tag_chuang_beard_1990 ) then

   axis_ratio=1.01668-0.09806*d_cm-2.52686*d_cm**2+ &
              +3.75061*d_cm**3-1.68692*d_cm**4
   if( axis_ratio < min_axis_rain ) axis_ratio = min_axis_rain
   if( axis_ratio > max_axis_rain ) axis_ratio = max_axis_rain
   if(  d_cm < cloud_thresh_cm    ) axis_ratio = 1.0

 elseif( ityp == tag_Pruppacher_Beard_1970 ) then

   axis_ratio =1.03 - 0.62*d_mm
   if( axis_ratio < min_axis_rain ) axis_ratio = min_axis_rain
   if( axis_ratio > max_axis_rain ) axis_ratio = max_axis_rain
   if(  d_cm < cloud_thresh_cm ) axis_ratio = 1.00

 elseif( ityp == tag_jameson_1983 )  then

   axis_ratio =1.01-0.23*d_mm
   if( axis_ratio < min_axis_rain ) axis_ratio = min_axis_rain
   if( axis_ratio > max_axis_rain ) axis_ratio = max_axis_rain
   if(  d_cm < cloud_thresh_cm ) axis_ratio = 1.00

 elseif( ityp == tag_green_1975 ) then

    axis_ratio =1.6993*d_cm**5-5.3728*d_cm**4+6.4728*d_cm**3- &
                +3.3406*d_cm**2+0.033*d_cm+1.0001
   if( axis_ratio < min_axis_rain ) axis_ratio = min_axis_rain
   if( axis_ratio > max_axis_rain ) axis_ratio = max_axis_rain
    if(  d_cm < cloud_thresh_cm ) axis_ratio = 1.00

 elseif( ityp == tag_pruppacher_pitter_1971 ) then

     axis_ratio = -4.1901*d_cm**4+8.1542*d_cm**3-4.9659*d_cm**2+ &
                  +0.3924*d_cm+0.9888
   if( axis_ratio < min_axis_rain ) axis_ratio = min_axis_rain
   if( axis_ratio > max_axis_rain ) axis_ratio = max_axis_rain
     if(  d_cm < cloud_thresh_cm ) axis_ratio = 1.00

 elseif( ityp == tag_jones_i_1959 ) then

! Jones (1959) I:  oblate raindrops a/b vs. D best fit line similar to his eyeball fit
     axis_ratio =1.03-0.57*d_mm
   if( axis_ratio < min_axis_rain ) axis_ratio = min_axis_rain
   if( axis_ratio > max_axis_rain ) axis_ratio = max_axis_rain
     if(  d_cm < cloud_thresh_cm ) axis_ratio = 1.00

 elseif( ityp == tag_jones_ii_1959 ) then

!  Jones (1959) II:  ALL drops a/b vs. D best fit line similar to his eyeball fit
     axis_ratio =1.13-0.69*d_mm
     if( axis_ratio < min_axis_rain ) axis_ratio = min_axis_rain
    if( axis_ratio > max_axis_rain ) axis_ratio = max_axis_rain
     if(  d_cm < cloud_thresh_cm ) axis_ratio = 1.00

 elseif( ityp == tag_goddard_cherry_1984 ) then

!     Goddard and Cherry (1984) 3rd order polynomial fit that is
!     only valid to about 3.2 mm.  Composite Goddard and Cherry
!     with Chuang and Beard (1990) in order to utilize.

     if(d_cm.lt.0.32)then
          axis_ratio =2.3434*d_cm**3-3.3039*d_cm**2+0.2771*d_cm+1.0036
     else
          axis_ratio =1.01668-0.09806*d_cm-2.52686*d_cm**2+&
                3.75061*d_cm**3-1.68692*d_cm**4
     endif
     if( axis_ratio < min_axis_rain ) axis_ratio = min_axis_rain
     if( axis_ratio > max_axis_rain ) axis_ratio = max_axis_rain
     if(  d_cm < cloud_thresh_cm ) axis_ratio = 1.00

 elseif( ityp == tag_brandes_2002 ) then

     !  Brandes (2002)
     axis_ratio =0.9951+0.02510*d_mm-0.03644*d_mm**2+&
           0.005303*d_mm**3-0.0002492*d_mm**4

     if( axis_ratio < min_axis_rain ) axis_ratio = min_axis_rain
     if( axis_ratio > max_axis_rain ) axis_ratio = max_axis_rain
     if(  d_cm < cloud_thresh_cm    ) axis_ratio = 1.00

 elseif( ityp == tag_ryzhkov_gh_2011 ) then

     axis_ratio = max(0.8, 1.0 - 0.02 * d_mm)  !miniam value is 0.8

 elseif( ityp == tag_hail_knight_1986 ) then

   !oblate hail from weighted mean values of knight 1986 (all three fields values are averaged)
    axis_ratio = max(0.725, 0.897 - 0.0008 * d_mm - 0.0002* d_mm**2 )   

 elseif( ityp == tag_snow_mc3e_hvps3 ) then

    !snow aggregate axis ratio relationship from MC3E citation HVPS-3 measurment

    if( d_mm > 10. ) d_mm = 10.  !stay diameter less than 10mm 

    axis_ratio = 0.003 * d_mm**2 - 0.05 * d_mm + 0.70 

 else

  ! do nothing

 endif

 return
 end subroutine spectrum_axis_ratio

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine bscat_mueller(mat,bscat) 
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Compute backscattering coefficients from Mueller matrix.  
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: mat(4,4)
 real(sdsu_fpd),intent(out) ::  bscat
 bscat = mat(1,1) - mat(1,2) - mat(2,1) + mat(2,2) ! backscattering coeff
 end subroutine bscat_mueller

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine latlonhgt2cartesian(lat,lon,hgt_km,v)
 implicit none
!-----------------------------------------------------------------------------------
! Comments: Using latitude and longitude, find out vector point on reference ellipsoid
!           in ECEF cartesian coordinate. 
!
!           Note that (x,y,z) = (rs, 0, 0) => lat= 0deg, lon= 0deg
!                     (x,y,z) = ( 0,rs, 0) => lat= 0deg, lon=90deg
!                     (x,y,z) = ( 0, 0,rs) => lat=90deg, lon= 0deg
!
! History:
! 12/2010  Toshi Matsui@NASA GSFC : Initial. 
!
! References:
!-----------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: lat, lon !earth's latitude and longitude [rad] 
 real(sdsu_fpd),intent(in) :: hgt_km   !hgt above elipsoid [km]
 real(sdsu_fpd),intent(out) :: v(3) !cartesian cordinate position [km]
 real(sdsu_fpd) :: re !distance from earth's core to ellipsoidal surface  [km]

!
! Get distance from the earth's core to the ellipoidal surface. 
!
 call earth_radius(lat, re)

!
! Derive vector x,y,z
!
 v(1) = (re+hgt_km)*cos(lat)*cos(lon)  ! x [km]

 v(2) = (re+hgt_km)*cos(lat)*sin(lon)  ! y [km]

 v(3) = (re+hgt_km)*sin(lat)           ! z [km]


 return
 end subroutine latlonhgt2cartesian

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine earth_radius( lat, re )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Compute earth's radius as spherical oblate (for a given latitude). 
!           Namely, from the earth's core of earth to reference ellipsoid based ono WGS 1984. 
! History:
! 10/2010  Toshi Matsui@NASA GSFC : Initial
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: lat ! latitude [rad]
 real(sdsu_fpd),intent(out) :: re ! distance between earth's core and ellipsoid [km]

 real(sdsu_fpd) :: a,b,c,d
 real(sdsu_fpd),parameter :: req=6378.1370d0     ! earth's equatorial radius [km]
 real(sdsu_fpd),parameter :: rpl=6356.7523142d0  ! earth's polar radius [km]

 a = req*req*dcos(lat)*req*req*dcos(lat)
 b = rpl*rpl*dsin(lat)*rpl*rpl*dsin(lat)
 c = req*dcos(lat)*req*dcos(lat)
 d = rpl*dsin(lat)*rpl*dsin(lat)

 re = dsqrt( (a+b)/(c+d) )   ! distance between earth's core and ellipsoid [km]

 return
 end subroutine earth_radius

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_distance( v1, v2, dist)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Derive distance between v1 and v2 points in ECEF cordinate.
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: v1(3)  !3D vector 1
 real(sdsu_fpd),intent(in) :: v2(3)  !3D vector 2
 real(sdsu_fpd),intent(out) :: dist

!
! distance [unit is same as v1 and v2 
!
 dist = sqrt( (v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 + (v1(3)-v2(3))**2 )

 return
 end subroutine get_distance

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  elevation_angle( v_inst_top , v_inst_grnd, v_grid, elev, dist )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For given ground instrument, instrument normal direction, and CRM grid locaitons, 
!            this subroutine compute elevation angle and sitance. 
!
!  length a = CRM grid ~ earth core
!         b = satellite ~ earth core 
!         c = satellite  ~ CRM grid
!
!  angle  alpha for a  (sensor viewing angle)
!         beta  for b  (90 - incidence angle)
!         gam   for c  (earth core angle)
!  
! History:
! 09/2016  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fpd), intent(in) :: v_inst_top (1:3)  ! vector of instrumental top [km]
 real(sdsu_fpd), intent(in) :: v_inst_grnd(1:3)  ! vector of instrumet grid [km]
 real(sdsu_fpd), intent(in) :: v_grid     (1:3)  ! vector of grid point [km]

 real(sdsu_fps), intent(out) :: elev !elevation angle [deg]
 real(sdsu_fps), intent(out) :: dist !distance between instrument and grid [km]

 real(sdsu_fpd) :: alpha,beta,gam ![rad]
 real(sdsu_fpd) :: a,b,c  ![km]


 call get_distance( v_inst_grnd , v_grid, a)
 call get_distance( v_inst_top  , v_grid, b)
 call get_distance( v_inst_top  , v_inst_grnd, c)

!
! check if this is on nadir
!
 if( (a+b) - c < 1.d-10 )  then

   alpha =   0.d0 * degrad
   beta  =   0.d0 * degrad
   gam   = 180.d0 * degrad

 elseif( (a+b) - c < 1.d-10 )  then

   alpha = 180.d0 * degrad
   beta  =   0.d0 * degrad
   gam   =   0.d0 * degrad

 else

   alpha = dacos( ( b*b + c*c - a*a ) / (2.*b*c) )  ![rad]
   beta  = dacos( ( a*a + c*c - b*b ) / (2.*a*c) )  ![rad]  --> elevation angle = (90-beta*raddeg)
   gam   = dacos( ( a*a + b*b - c*c ) / (2.*a*b) )  ![rad]

 endif

 elev = real( 90.d0 - beta*raddeg )  !elevation angle [deg]
 dist = real( a )                    !distance [km]

 return
 end subroutine elevation_angle 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  angle_to_inst( v_grid_next , v_inst_grnd, v_grid, angle )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For given ground instrument, CRM grid, and neiboriing CRM grid, 
!            this subroutine compute angles to compute radial velocity. 
!
!  length a = CRM grid ~ earth core
!         b = satellite ~ earth core 
!         c = satellite  ~ CRM grid
!
!  angle  alpha for a  (sensor viewing angle)
!         beta  for b  (90 - incidence angle)
!         gam   for c  (earth core angle)
!  
! History:
! 09/2016  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fpd), intent(in) :: v_grid_next(1:3)  ! vector of neibor grid [km]
 real(sdsu_fpd), intent(in) :: v_inst_grnd(1:3)  ! vector of instrumet grid [km]
 real(sdsu_fpd), intent(in) :: v_grid     (1:3)  ! vector of grid point [km]

 real(sdsu_fps), intent(out) :: angle !angle to instrumental direction  [rad]

 real(sdsu_fpd) :: alpha,beta,gam ![rad]
 real(sdsu_fpd) :: a,b,c  ![km]


 call get_distance( v_inst_grnd , v_grid, a)
 call get_distance( v_grid_next , v_grid, b)
 call get_distance( v_grid_next , v_inst_grnd, c)

!
! check if extreme angle
!
 if( a == (b+c)  )  then

   alpha = 180.d0 * degrad
   beta  =   0.d0 * degrad
   gam   =   0.d0 * degrad

 elseif ( (c+b) == a  )  then

   alpha =   0.d0 * degrad
   beta  =   0.d0 * degrad
   gam   = 180.d0 * degrad

 else

   alpha = dacos( ( b*b + c*c - a*a ) / (2.*b*c) )  ![rad]
   beta  = dacos( ( a*a + c*c - b*b ) / (2.*a*c) )  ![rad]  --> elevation angle = (90-beta*raddeg)
   gam   = dacos( ( a*a + b*b - c*c ) / (2.*a*b) )  ![rad]

 endif

!
! output angle [rad]
!
 angle = real( gam )


 return
 end subroutine angle_to_inst

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine const_hgt_level_atmos( mxlyr_crm, hgt_stag_crm, atmos_in, &
                             mxlyr_inst, hgt_stag_inst,  atmos_out )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Interpolate CRM streching vertical level into constant height levels
!  for atmospheric parameters.  
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 integer,intent(in)        :: mxlyr_crm 
 real(sdsu_fps),intent(in) :: hgt_stag_crm(0:mxlyr_crm)    ! height leve at CRM interface [km] 
 type ( atmos_parameter ),intent(in)  :: atmos_in(1:mxlyr_crm) !atmos parameter
 integer,intent(in)        :: mxlyr_inst                   ! maximum layer for radar instrument.
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)  ! instrument measurement height at interface [km]
 type ( atmos_parameter ),intent(out)  :: atmos_out(1:mxlyr_inst)  !atmos parameter at instrumental level 

 integer,parameter :: np_max = 30
 real(sdsu_fps) :: bundle_crm (np_max,mxlyr_crm)
 real(sdsu_fps) :: bundle_inst(np_max,mxlyr_inst)

 integer :: k, np


 np = 1
 bundle_crm(np,1:mxlyr_crm) = atmos_in(1:mxlyr_crm)%t_air     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = atmos_in(1:mxlyr_crm)%press     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = atmos_in(1:mxlyr_crm)%rh      ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = atmos_in(1:mxlyr_crm)%u      ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = atmos_in(1:mxlyr_crm)%v      ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = atmos_in(1:mxlyr_crm)%w         ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = atmos_in(1:mxlyr_crm)%rho_dair  ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = atmos_in(1:mxlyr_crm)%rain_rate ; np = np + 1


!
! refine vertical layer
!
 call refine_range_bin( mxlyr_crm, hgt_stag_crm(0:mxlyr_crm), atmos_in(1:mxlyr_crm)%dhgt, &
                        mxlyr_inst,  hgt_stag_inst, &
                        np-1, bundle_crm(1:np-1,1:mxlyr_crm), bundle_inst(1:np-1,1:mxlyr_inst) )


 np = 1
 atmos_out(1:mxlyr_inst)%t_air = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 atmos_out(1:mxlyr_inst)%press = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 atmos_out(1:mxlyr_inst)%rh    = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 atmos_out(1:mxlyr_inst)%u     = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 atmos_out(1:mxlyr_inst)%v     = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 atmos_out(1:mxlyr_inst)%w     = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 atmos_out(1:mxlyr_inst)%rho_dair = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 atmos_out(1:mxlyr_inst)%rain_rate = bundle_inst(np,1:mxlyr_inst)  ; np = np + 1

 do k = 1, mxlyr_inst
    atmos_out(k)%hgt = 0.5 * ( hgt_stag_inst(k-1) + hgt_stag_inst(k) ) 
 enddo

 return
 end subroutine const_hgt_level_atmos

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine const_hgt_level_gen( mxlyr_crm, hgt_stag_crm, atmos_in, q_gen_in, re_gen_in,  &
                                 mxlyr_inst, hgt_stag_inst,  q_gen_out, re_gen_out , &
                                 rho_graupel_in, rho_hail_in, rho_graupel_out, rho_hail_out )

 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Interpolate CRM streching vertical level into constant height levels
!  for generalized bulk microphysics.   
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 integer,intent(in)        :: mxlyr_crm
 real(sdsu_fps),intent(in) :: hgt_stag_crm(0:mxlyr_crm)    ! height leve at CRM interface [km] 
 type ( atmos_parameter ),intent(in)  :: atmos_in(1:mxlyr_crm) !atmos parameter
 type ( particle_gen ),intent(in)  :: q_gen_in(1:mxlyr_crm) ! mixing ratio [g/m3]
 type ( particle_gen ),intent(in)  :: re_gen_in(1:mxlyr_crm) ! effectie radius [um]

 integer,intent(in)        :: mxlyr_inst                   ! maximum layer for radar instrument.
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)  ! instrument measurement height at interface [km]
 type ( particle_gen ),intent(inout) :: q_gen_out(1:mxlyr_inst)  !mixing ratio at instrumental level [g/m3]
 type ( particle_gen ),intent(inout) :: re_gen_out(1:mxlyr_inst) !effectie radius at instrumental level [um]

 real(sdsu_fps),intent(in),optional  :: rho_graupel_in(1:mxlyr_crm)  ! graupel density [kg/m3]
 real(sdsu_fps),intent(in),optional  :: rho_hail_in(1:mxlyr_crm)     ! hail density [kg/m3]
 real(sdsu_fps),intent(out),optional :: rho_graupel_out(1:mxlyr_inst) 
 real(sdsu_fps),intent(out),optional :: rho_hail_out(1:mxlyr_inst)

 integer,parameter :: np_max = 30 
 real(sdsu_fps) :: bundle_crm (np_max,mxlyr_crm)
 real(sdsu_fps) :: bundle_inst(np_max,mxlyr_inst)

 integer :: k, np


 np = 1
 bundle_crm(np,1:mxlyr_crm) = q_gen_in(1:mxlyr_crm)%cloud    ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_gen_in(1:mxlyr_crm)%rain     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_gen_in(1:mxlyr_crm)%ice      ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_gen_in(1:mxlyr_crm)%snow     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_gen_in(1:mxlyr_crm)%graupel  ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_gen_in(1:mxlyr_crm)%hail     ; np = np + 1

 bundle_crm(np,1:mxlyr_crm) = re_gen_in(1:mxlyr_crm)%cloud    ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_gen_in(1:mxlyr_crm)%rain     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_gen_in(1:mxlyr_crm)%ice      ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_gen_in(1:mxlyr_crm)%snow     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_gen_in(1:mxlyr_crm)%graupel  ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_gen_in(1:mxlyr_crm)%hail     ; np = np + 1

 if( present(rho_graupel_in) ) bundle_crm(np,1:mxlyr_crm) = rho_graupel_in(1:mxlyr_crm)   ; np = np + 1
 if( present(rho_hail_in)    ) bundle_crm(np,1:mxlyr_crm) = rho_hail_in   (1:mxlyr_crm)   ; np = np + 1

!
! refine vertical layer
!
 call refine_range_bin( mxlyr_crm, hgt_stag_crm(0:mxlyr_crm), atmos_in(1:mxlyr_crm)%dhgt, &
                        mxlyr_inst,  hgt_stag_inst, &
                        np-1, bundle_crm(1:np-1,1:mxlyr_crm), bundle_inst(1:np-1,1:mxlyr_inst) )

 np = 1
 q_gen_out(1:mxlyr_inst)%cloud = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_gen_out(1:mxlyr_inst)%rain  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_gen_out(1:mxlyr_inst)%ice   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_gen_out(1:mxlyr_inst)%snow  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_gen_out(1:mxlyr_inst)%graupel=bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_gen_out(1:mxlyr_inst)%hail  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1

 re_gen_out(1:mxlyr_inst)%cloud = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_gen_out(1:mxlyr_inst)%rain  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_gen_out(1:mxlyr_inst)%ice   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_gen_out(1:mxlyr_inst)%snow  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_gen_out(1:mxlyr_inst)%graupel=bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_gen_out(1:mxlyr_inst)%hail  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1

!after interpoation, predicted rhoe become out of boundary, --> correct here. 
 if( present(rho_graupel_out) ) then
     rho_graupel_out(1:mxlyr_inst)  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
     do k = 1, mxlyr_inst
        rho_graupel_out(k) = max(170., min(900., rho_graupel_out(k) ) )  !redfine boundary
     enddo
 endif

 if( present(rho_hail_out)    ) then
     rho_hail_out   (1:mxlyr_inst)  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
     do k = 1, mxlyr_inst
        rho_hail_out(k) = max(500. ,min(900., rho_hail_out(k) ) )  !redfine boundary
     enddo
 endif

 return
 end subroutine const_hgt_level_gen

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine const_hgt_level_p3( mxlyr_crm, hgt_stag_crm, atmos_in, q_gen_in, re_gen_in,  &
                                 mxlyr_inst, hgt_stag_inst,  q_gen_out, re_gen_out , &
                                 rhorime_in, frime_in, rhorime_out, frime_out , &
                                 mu_cloud_in, n_ice_in, mu_cloud_out, n_ice_out )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Interpolate CRM streching vertical level into constant height levels
!  for P3 bulk microphysics.   
!
! History:
! 07/2019  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 integer,intent(in)        :: mxlyr_crm
 real(sdsu_fps),intent(in) :: hgt_stag_crm(0:mxlyr_crm)    ! height leve at CRM interface [km] 
 type ( atmos_parameter ),intent(in)  :: atmos_in(1:mxlyr_crm) !atmos parameter
 type ( particle_gen ),intent(in)  :: q_gen_in(1:mxlyr_crm) ! mixing ratio [g/m3]
 type ( particle_gen ),intent(in)  :: re_gen_in(1:mxlyr_crm) ! effectie radius [um]

 integer,intent(in)        :: mxlyr_inst                   ! maximum layer for radar instrument.
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)  ! instrument measurement height at interface [km]
 type ( particle_gen ),intent(inout) :: q_gen_out(1:mxlyr_inst)  !mixing ratio at instrumental level [g/m3]
 type ( particle_gen ),intent(inout) :: re_gen_out(1:mxlyr_inst) !effectie radius at instrumental level [um]

 real(sdsu_fps),intent(in)  :: rhorime_in(1:mxlyr_crm)  ! graupel density [kg/m3]
 real(sdsu_fps),intent(in)  :: frime_in(1:mxlyr_crm)     ! hail density [kg/m3]
 real(sdsu_fps),intent(out) :: rhorime_out(1:mxlyr_inst) 
 real(sdsu_fps),intent(out) :: frime_out  (1:mxlyr_inst)

 real(sdsu_fps),intent(in)  :: mu_cloud_in (1:mxlyr_crm)   ! cloud mu [-]
 real(sdsu_fps),intent(in)  :: n_ice_in    (1:mxlyr_crm)   ! ice number conc.  [#/m3]
 real(sdsu_fps),intent(out) :: mu_cloud_out(1:mxlyr_inst)
 real(sdsu_fps),intent(out) :: n_ice_out   (1:mxlyr_inst)

 integer,parameter :: np_max = 30 
 real(sdsu_fps) :: bundle_crm (np_max,mxlyr_crm)
 real(sdsu_fps) :: bundle_inst(np_max,mxlyr_inst)

 integer :: k, np

 np = 1
 bundle_crm(np,1:mxlyr_crm) = q_gen_in(1:mxlyr_crm)%cloud    ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_gen_in(1:mxlyr_crm)%rain     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_gen_in(1:mxlyr_crm)%ice      ; np = np + 1

 bundle_crm(np,1:mxlyr_crm) = re_gen_in(1:mxlyr_crm)%cloud    ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_gen_in(1:mxlyr_crm)%rain     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_gen_in(1:mxlyr_crm)%ice      ; np = np + 1

 bundle_crm(np,1:mxlyr_crm) = rhorime_in(1:mxlyr_crm)         ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = frime_in  (1:mxlyr_crm)         ; np = np + 1

 bundle_crm(np,1:mxlyr_crm) = mu_cloud_in(1:mxlyr_crm)        ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = n_ice_in   (1:mxlyr_crm)        ; np = np + 1

!
! refine vertical layer
!
 call refine_range_bin( mxlyr_crm, hgt_stag_crm(0:mxlyr_crm), atmos_in(1:mxlyr_crm)%dhgt, &
                        mxlyr_inst,  hgt_stag_inst, &
                        np-1, bundle_crm(1:np-1,1:mxlyr_crm), bundle_inst(1:np-1,1:mxlyr_inst) )

 np = 1
 q_gen_out(1:mxlyr_inst)%cloud = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_gen_out(1:mxlyr_inst)%rain  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_gen_out(1:mxlyr_inst)%ice   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1

 re_gen_out(1:mxlyr_inst)%cloud = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_gen_out(1:mxlyr_inst)%rain  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_gen_out(1:mxlyr_inst)%ice   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1

!after interpoation, predicted rhoe become out of boundary, --> correct here. 
 rhorime_out(1:mxlyr_inst)  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 do k = 1, mxlyr_inst
    rhorime_out(k) = max(50., min(900., rhorime_out(k) ) )  !redfine boundary
 enddo

 frime_out   (1:mxlyr_inst)  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 do k = 1, mxlyr_inst
   frime_out(k) = max(0. ,min(1., frime_out(k) ) )  !redfine boundary
 enddo

 mu_cloud_out   (1:mxlyr_inst)  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 do k = 1, mxlyr_inst
   mu_cloud_out(k) = max(2. ,min(15., mu_cloud_out(k) ) )  !redfine boundary
 enddo

 n_ice_out      (1:mxlyr_inst)  = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 do k = 1, mxlyr_inst
   n_ice_out(k) = max(0. , n_ice_out(k) )  !redfine boundary
 enddo
 
 return
 end subroutine const_hgt_level_p3

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine const_hgt_level_rams( mxlyr_crm, hgt_stag_crm, atmos_in, q_rams_in, re_rams_in,  &
                                  mxlyr_inst, hgt_stag_inst,  q_rams_out, re_rams_out )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Interpolate CRM streching vertical level into constant height levels
!  for generalized RAMS microphysics.   
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 integer,intent(in)        :: mxlyr_crm
 real(sdsu_fps),intent(in) :: hgt_stag_crm(0:mxlyr_crm)    ! height leve at CRM interface [km] 
 type ( atmos_parameter ),intent(in)  :: atmos_in(1:mxlyr_crm) !atmos parameter
 type ( particle_rams ),intent(in)  :: q_rams_in(1:mxlyr_crm) ! mixing ratio [g/m3]
 type ( particle_rams ),intent(in)  :: re_rams_in(1:mxlyr_crm) ! effectie radius [um]

 integer,intent(in)        :: mxlyr_inst                   ! maximum layer for radar instrument.
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)  ! instrument measurement height at interface [km]
 type ( particle_rams ),intent(inout) :: q_rams_out(1:mxlyr_inst)  !mixing ratio at instrumental level [g/m3]
 type ( particle_rams ),intent(inout) :: re_rams_out(1:mxlyr_inst) !effectie radius at instrumental level [um]

 integer,parameter :: np_max = 30 
 real(sdsu_fps) :: bundle_crm (np_max,mxlyr_crm)
 real(sdsu_fps) :: bundle_inst(np_max,mxlyr_inst)

 integer :: k, np

!
! bundle input
!
 np = 1
 bundle_crm(np,1:mxlyr_crm) = q_rams_in(1:mxlyr_crm)%cloud1   ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_rams_in(1:mxlyr_crm)%cloud2   ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_rams_in(1:mxlyr_crm)%rain     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_rams_in(1:mxlyr_crm)%ice1     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_rams_in(1:mxlyr_crm)%ice2     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_rams_in(1:mxlyr_crm)%snow     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_rams_in(1:mxlyr_crm)%graupel  ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_rams_in(1:mxlyr_crm)%hail     ; np = np + 1

 bundle_crm(np,1:mxlyr_crm) = re_rams_in(1:mxlyr_crm)%cloud1   ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_rams_in(1:mxlyr_crm)%cloud2   ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_rams_in(1:mxlyr_crm)%rain     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_rams_in(1:mxlyr_crm)%ice1     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_rams_in(1:mxlyr_crm)%ice2     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_rams_in(1:mxlyr_crm)%snow     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_rams_in(1:mxlyr_crm)%graupel  ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_rams_in(1:mxlyr_crm)%hail     ; np = np + 1


!
! refine vertical layer
!
 call refine_range_bin( mxlyr_crm, hgt_stag_crm(0:mxlyr_crm), atmos_in(1:mxlyr_crm)%dhgt, &
                        mxlyr_inst,  hgt_stag_inst, &
                        np-1, bundle_crm(1:np-1,1:mxlyr_crm), bundle_inst(1:np-1,1:mxlyr_inst) )

!
! release output
!
 np = 1
 q_rams_out(1:mxlyr_inst)%cloud1 = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_rams_out(1:mxlyr_inst)%cloud2 = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_rams_out(1:mxlyr_inst)%rain   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_rams_out(1:mxlyr_inst)%ice1   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_rams_out(1:mxlyr_inst)%ice2   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_rams_out(1:mxlyr_inst)%snow   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_rams_out(1:mxlyr_inst)%graupel= bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_rams_out(1:mxlyr_inst)%hail   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1

 re_rams_out(1:mxlyr_inst)%cloud1 = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_rams_out(1:mxlyr_inst)%cloud2 = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_rams_out(1:mxlyr_inst)%rain   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_rams_out(1:mxlyr_inst)%ice1   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_rams_out(1:mxlyr_inst)%ice2   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_rams_out(1:mxlyr_inst)%snow   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_rams_out(1:mxlyr_inst)%graupel= bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_rams_out(1:mxlyr_inst)%hail   = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1


 return
 end subroutine const_hgt_level_rams

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine const_hgt_level_sbm( mxlyr_crm, hgt_stag_crm, atmos_in, q_sbm_in, re_sbm_in,  &
                             mxlyr_inst, hgt_stag_inst,  q_sbm_out, re_sbm_out )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Interpolate CRM streching vertical level into constant height levels
!  for HUCM SBM bulk components.   
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 integer,intent(in)        :: mxlyr_crm
 real(sdsu_fps),intent(in) :: hgt_stag_crm(0:mxlyr_crm)    ! height leve at CRM interface [km] 
 type ( atmos_parameter ),intent(in)  :: atmos_in(1:mxlyr_crm) !atmos parameter
 type ( particle_sbm ),intent(in)  :: q_sbm_in(1:mxlyr_crm) ! mixing ratio [g/m3]
 type ( particle_sbm ),intent(in)  :: re_sbm_in(1:mxlyr_crm) ! effectie radius [um]

 integer,intent(in)        :: mxlyr_inst                   ! maximum layer for radar instrument.
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)  ! instrument measurement height at interface [km]
 type ( particle_sbm ),intent(inout) :: q_sbm_out(1:mxlyr_inst)  !mixing ratio at instrumental level [g/m3]
 type ( particle_sbm ),intent(inout) :: re_sbm_out(1:mxlyr_inst) !effectie radius at instrumental level [um]

 integer,parameter :: np_max = 30 
 real(sdsu_fps) :: bundle_crm (np_max,mxlyr_crm)
 real(sdsu_fps) :: bundle_inst(np_max,mxlyr_inst)

 integer :: k, np


 np = 1
 bundle_crm(np,1:mxlyr_crm) = q_sbm_in(1:mxlyr_crm)%liq      ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_sbm_in(1:mxlyr_crm)%ice_col  ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_sbm_in(1:mxlyr_crm)%ice_pla  ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_sbm_in(1:mxlyr_crm)%ice_den  ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_sbm_in(1:mxlyr_crm)%snow     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_sbm_in(1:mxlyr_crm)%graupel  ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = q_sbm_in(1:mxlyr_crm)%hail     ; np = np + 1

 bundle_crm(np,1:mxlyr_crm) = re_sbm_in(1:mxlyr_crm)%liq     ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_sbm_in(1:mxlyr_crm)%ice_col ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_sbm_in(1:mxlyr_crm)%ice_pla ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_sbm_in(1:mxlyr_crm)%ice_den ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_sbm_in(1:mxlyr_crm)%snow    ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_sbm_in(1:mxlyr_crm)%graupel ; np = np + 1
 bundle_crm(np,1:mxlyr_crm) = re_sbm_in(1:mxlyr_crm)%hail    ; np = np + 1

!
! refine vertical layer
!
 call refine_range_bin( mxlyr_crm, hgt_stag_crm(0:mxlyr_crm), atmos_in(1:mxlyr_crm)%dhgt, &
                        mxlyr_inst,  hgt_stag_inst, &
                        np-1, bundle_crm(1:np-1,1:mxlyr_crm), bundle_inst(1:np-1,1:mxlyr_inst) )


 np = 1
 q_sbm_out(1:mxlyr_inst)%liq     = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_sbm_out(1:mxlyr_inst)%ice_col = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_sbm_out(1:mxlyr_inst)%ice_pla = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_sbm_out(1:mxlyr_inst)%ice_den = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_sbm_out(1:mxlyr_inst)%snow    = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_sbm_out(1:mxlyr_inst)%graupel = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 q_sbm_out(1:mxlyr_inst)%hail    = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1

 re_sbm_out(1:mxlyr_inst)%liq     = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_sbm_out(1:mxlyr_inst)%ice_col = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_sbm_out(1:mxlyr_inst)%ice_pla = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_sbm_out(1:mxlyr_inst)%ice_den = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_sbm_out(1:mxlyr_inst)%snow    = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_sbm_out(1:mxlyr_inst)%graupel = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1
 re_sbm_out(1:mxlyr_inst)%hail    = bundle_inst(np,1:mxlyr_inst)   ; np = np + 1


 return
 end subroutine const_hgt_level_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine const_hgt_level_sbmbin( mxlyr_crm, hgt_stag_crm, atmos_in, n_sbm_in ,  &
                                 fmelt_sbm_snow_in, fmelt_sbm_graupel_in, &
                                 fmelt_sbm_hail_in, frime_sbm_snow_in, &
                                 mxlyr_inst, hgt_stag_inst, &
                                 n_sbm_out ,  &
                                 fmelt_sbm_snow_out, fmelt_sbm_graupel_out, &
                                 fmelt_sbm_hail_out, frime_sbm_snow_out )

 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Interpolate CRM streching vertical level into constant height levels
!  for HUCM SBM bin components.    
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 integer,intent(in)        :: mxlyr_crm
 real(sdsu_fps),intent(in) :: hgt_stag_crm(0:mxlyr_crm)    ! height leve at CRM interface [km] 
 type ( atmos_parameter ),intent(in)  :: atmos_in(1:mxlyr_crm) !atmos parameter
 
 type ( particle_sbm ),intent(in)  :: n_sbm_in(1:mxlyr_crm,1:nbin) ! mixing ratio [g/m3]

 real(sdsu_fps),intent(in),dimension(1:mxlyr_crm,1:nbin) :: &
         fmelt_sbm_snow_in     , &  ! melting fraction of snow [-]
         fmelt_sbm_graupel_in  , &  ! melting fraction of graupel [-]
         fmelt_sbm_hail_in     , &  ! melting fraction of hail [-]
         frime_sbm_snow_in          ! riming fraction of snow [-]

 integer,intent(in)        :: mxlyr_inst                   ! maximum layer for radar instrument.
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)  ! instrument measurement height at interface [km]

 type ( particle_sbm ),intent(out)  :: n_sbm_out(1:mxlyr_inst,1:nbin) ! mixing ratio [g/m3]
 real(sdsu_fps),intent(out),dimension(1:mxlyr_inst,1:nbin) :: &  
         fmelt_sbm_snow_out     , &  ! melting fraction of snow [-]
         fmelt_sbm_graupel_out  , &  ! melting fraction of graupel [-]
         fmelt_sbm_hail_out     , &  ! melting fraction of hail [-]
         frime_sbm_snow_out          ! riming fraction of snow [-]

 integer,parameter :: np_max = 30 
 real(sdsu_fps) :: bundle_crm (mxlyr_crm, nbin, np_max)
 real(sdsu_fps) :: bundle_inst(mxlyr_inst, nbin, np_max)

 integer :: k, np, n

!
! 
!
 np = 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = n_sbm_in(1:mxlyr_crm,1:nbin)%liq      ; np = np + 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = n_sbm_in(1:mxlyr_crm,1:nbin)%ice_col  ; np = np + 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = n_sbm_in(1:mxlyr_crm,1:nbin)%ice_pla  ; np = np + 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = n_sbm_in(1:mxlyr_crm,1:nbin)%ice_den  ; np = np + 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = n_sbm_in(1:mxlyr_crm,1:nbin)%snow     ; np = np + 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = n_sbm_in(1:mxlyr_crm,1:nbin)%graupel  ; np = np + 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = n_sbm_in(1:mxlyr_crm,1:nbin)%hail     ; np = np + 1

! convert fraction into number conc.
 bundle_crm(1:mxlyr_crm,1:nbin,np) = fmelt_sbm_snow_in(1:mxlyr_crm,1:nbin) &
                                     * n_sbm_in(1:mxlyr_crm,1:nbin)%snow     ; np = np + 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = fmelt_sbm_graupel_in(1:mxlyr_crm,1:nbin) &
                                     * n_sbm_in(1:mxlyr_crm,1:nbin)%graupel  ; np = np + 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = fmelt_sbm_hail_in(1:mxlyr_crm,1:nbin) &
                                     * n_sbm_in(1:mxlyr_crm,1:nbin)%hail     ; np = np + 1
 bundle_crm(1:mxlyr_crm,1:nbin,np) = frime_sbm_snow_in(1:mxlyr_crm,1:nbin) &
                                     * n_sbm_in(1:mxlyr_crm,1:nbin)%snow     ; np = np + 1

!
! refine vertical layer
!
 call refine_range_bin( mxlyr_crm, hgt_stag_crm(0:mxlyr_crm), atmos_in(1:mxlyr_crm)%dhgt, &
                        mxlyr_inst,  hgt_stag_inst, np-1, nbin, &
                        bundle_crm(1:mxlyr_crm,1:nbin,1:np-1), bundle_inst(1:mxlyr_inst,1:nbin,1:np-1) )


 np = 1
 n_sbm_out(1:mxlyr_inst,1:nbin)%liq     = bundle_inst(1:mxlyr_inst,1:nbin,np)   ; np = np + 1
 n_sbm_out(1:mxlyr_inst,1:nbin)%ice_col = bundle_inst(1:mxlyr_inst,1:nbin,np)   ; np = np + 1
 n_sbm_out(1:mxlyr_inst,1:nbin)%ice_pla = bundle_inst(1:mxlyr_inst,1:nbin,np)   ; np = np + 1
 n_sbm_out(1:mxlyr_inst,1:nbin)%ice_den = bundle_inst(1:mxlyr_inst,1:nbin,np)   ; np = np + 1
 n_sbm_out(1:mxlyr_inst,1:nbin)%snow    = bundle_inst(1:mxlyr_inst,1:nbin,np)   ; np = np + 1
 n_sbm_out(1:mxlyr_inst,1:nbin)%graupel = bundle_inst(1:mxlyr_inst,1:nbin,np)   ; np = np + 1
 n_sbm_out(1:mxlyr_inst,1:nbin)%hail    = bundle_inst(1:mxlyr_inst,1:nbin,np)   ; np = np + 1

!
! convert number conc. into fraction.
!
 fmelt_sbm_snow_out  = 0. ; fmelt_sbm_graupel_out = 0. 
 fmelt_sbm_hail_out  = 0. ; frime_sbm_snow_out    = 0.

 do n = 1, nbin ; do k = 1, mxlyr_inst
    if( bundle_inst(k,n,5) > 0. ) &
       fmelt_sbm_snow_out(k,n)    = bundle_inst(k,n,8) / bundle_inst(k,n,5)
    if( bundle_inst(k,n,6) > 0. ) &
       fmelt_sbm_graupel_out(k,n) = bundle_inst(k,n,9)   / bundle_inst(k,n,6)
    if( bundle_inst(k,n,7) > 0. ) &
        fmelt_sbm_hail_out(k,n)    = bundle_inst(k,n,10)  / bundle_inst(k,n,7)
    if( bundle_inst(k,n,5) > 0. ) &
        frime_sbm_snow_out(k,n)    = bundle_inst(k,n,11)  / bundle_inst(k,n,5)

 enddo ; enddo

 return
 end subroutine const_hgt_level_sbmbin

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine derive_radial_velocity( u_in,v_in, w_in, elev, &
                                    i_in, j_in, k_in, hgt, v_grid, v_inst , vrad)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Derive radial velocity for Polaritric radar instrument. 
!
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: u_in, v_in, w_in     ! u,v,w wind [m/s]
 real(sdsu_fps),intent(in) :: elev                 ! [deg]
 integer,intent(in) :: i_in, j_in, k_in            ! i, j location
 real(sdsu_fps),intent(in) :: hgt                  ! height [km]
 real(sdsu_fpd),intent(in) :: v_inst(3), v_grid(3) ! cartesian cordinate of instrument and CRM grid.
 real(sdsu_fps),intent(out) :: vrad     ! radial velocity [m/s]

 real(sdsu_fps) :: vrad_temp     ! radial velocity [m/s]

 real(sdsu_fpd) :: v_egrid(3) , v_ngrid(3) , v_ugrid(3) ! cartesian cordinate of instrument and CRM grid.

 real(sdsu_fps) :: lat_egrid, lon_egrid !lat/lon [deg] of east side of grid
 real(sdsu_fps) :: lat_ngrid, lon_ngrid !lat/lon [deg] of north side of grid
 real(sdsu_fps) :: lat_ugrid, lon_ugrid !lat/lon [deg] of north side of grid

 real(sdsu_fps) :: eangle, nangle, uangle  !angle between instrument toward east, north, and up cell [deg]

 if( mxgridx == 1 ) call stop_sdsu(myrank, 'MSG derive_radial_velocity: mxgridx must be greater than one.')
 if( mxgridy == 1 ) call stop_sdsu(myrank, 'MSG derive_radial_velocity: mxgridy must be greater than one.')

! print*,'MSG derive_radial_velocity u,v,w=',u_in,v_in,w_in

 if( u_in == undefined .or.  &
     v_in == undefined .or.  &
     w_in == undefined ) then
    vrad = undefined
    return
 endif

 
!
! cartesian location of east grid (note this does not work for single grid)
!
 if( i_in == mxgridx ) then
   lat_egrid = lat_domain(i_in, j_in)
   lon_egrid = 2. * lon_domain(i_in, j_in) - lon_domain(i_in-1, j_in) 
 else
   lat_egrid = lat_domain(i_in+1, j_in)
   lon_egrid = lon_domain(i_in+1, j_in) 
 endif

!print*, 'grid',lat_domain(i_in, j_in), lon_domain(i_in, j_in)
!print*, 'egrid',lat_egrid, lon_egrid

  call latlonhgt2cartesian( lat_egrid*degrad, lon_egrid*degrad, &
                            dble(hgt), v_egrid )

!
! angle between grid-uwind and grid-instrument 
!
  call  angle_to_inst( v_egrid , v_inst, v_grid, eangle )  


!
! cartesian loaction of north grid
!
 if( j_in == mxgridy ) then
   lat_ngrid = 2. * lat_domain(i_in, j_in) - lat_domain(i_in, j_in-1) 
   lon_ngrid = lon_domain(i_in, j_in)
 else
   lat_ngrid = lat_domain(i_in, j_in+1) 
   lon_ngrid = lon_domain(i_in, j_in+1)
 endif

  call latlonhgt2cartesian( lat_ngrid*degrad, lon_ngrid*degrad, &
                            dble(hgt), v_ngrid )

!
! angle between grid-vwind and grid-instrument 
!
  call  angle_to_inst( v_ngrid , v_inst, v_grid, nangle )  

!print*,'nangle=',nangle*raddeg

!
! cartesian loaction of upward grid
!
  lat_ugrid = lat_domain(i_in,j_in) 
  lon_ugrid = lon_domain(i_in,j_in)

  call latlonhgt2cartesian( lat_ugrid*degrad, lon_ugrid*degrad, &
                            dble(hgt+1.0), v_ugrid )
! if(j_in == 1 .and. k_in == 10) print*, i_in, lat_ugrid, lon_ugrid, dble(hgt+1.0), v_ugrid

!
! angle between grid-wwind and grid-instrument 
!
  call  angle_to_inst( v_ugrid , v_inst, v_grid, uangle )  
!print*,'uangle=',uangle*raddeg

!
! derive radial velocity  
!
 vrad_temp = u_in*cos(eangle) + v_in*cos(nangle) + w_in*cos(uangle) 

 vrad = -1. * vrad_temp        ! negative toward instrument


!print*, 'vrad=',vrad, u_in*cos(eangle), v_in*cos(nangle), w_in*cos(uangle)

 return
 end subroutine derive_radial_velocity

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine viewing_angle( lat_inst, lon_inst, delev_inst, &
                           lat_grid, lon_grid , vang ) 
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Estimate radar viewing angle at surafce level. North direction is zero degree, and 
!  rotate clock wise from top-down view ranging from 0 to 359.999deg. 
! 
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: lat_inst, lon_inst   ! lat/long of instrument [rad]
 real(sdsu_fpd),intent(in) :: delev_inst            ! elevation instrument [km]
 real(sdsu_fpd),intent(in) :: lat_grid, lon_grid   ! lat/long of CRM grid [rad]
 real(sdsu_fps),intent(out) :: vang

 real(sdsu_fpd) :: v_inst(3), v_grid(3), v_north(3)  ! cartesian cordinate position [km]
 real(sdsu_fpd) :: lat_north  ! north of instrument [rad]

 real(sdsu_fps) :: angle_rad
!
! instrumental location
!
 call latlonhgt2cartesian( lat_inst, lon_inst, delev_inst, v_inst )

!
! grid location
!
 call latlonhgt2cartesian( lat_grid, lon_grid, delev_inst, v_grid )

!
! 0.1deg north of instrumental location (So, this algorithm does not work
! for North Pole.)
!
 lat_north = lat_inst + 0.1*degrad ![rad]
 if( lat_north*raddeg > 89.90 ) then
   call stop_sdsu(myrank,'MSG viewing_angle: instrument location cannot be 89.9deg North.')
 endif
 call latlonhgt2cartesian( lat_north, lon_inst, delev_inst, v_north )


!
! angle
!
 call  angle_to_inst( v_north , v_grid, v_inst, angle_rad )


!
! output angle
! 
 if( lon_inst < lon_grid ) then
    vang = angle_rad * real(raddeg )   !viewing angle
 else
    vang = 360.0 - angle_rad * real(raddeg )   !viewing angle
 endif

 return
 end subroutine viewing_angle

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_polarris_grads
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out POLARRIS 3D outpu in GrADS format. 
! 
! History:
! 09/2016  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: ierr
 integer :: i,j,k,irec
 character(len=20) :: typ
 integer,parameter :: io = 110
 integer :: is,ie,js,je,ks,ke, nf
 character(len=2) :: ff

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end

 if( use_crm_level_polarris ) then
   ks=myk_start ; ke=mxlyr
 else
   ks=myk_start ; ke=mxlyr_polarris
 endif

!
! write out CRM surface input for grads binarry format
!
  typ = '.POLARRIS'

!
! write out POLARRIS output in grads binarry format
!
  sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//trim(typ)//&
                 '.'//trim(polarris_axisorient_name)//trim(output_suffix)//'.bin'


#if MPI == 2
 if(masterproc) &
#endif
  print("(a)"),'MSG write_out_polarris_grads: write POLARRIS 3D output file in grads format ',&
  ' ->',trim(sdsu_io_file)

!
! open grads file
!
#if MPI == 2
 if(masterproc) &
#endif
  open(unit= io, file = sdsu_io_file, access='direct', status='replace', &
        recl=mxgridx*mxgridy*(ke-ks+1)*4 )
  irec = 1

!
! write binary 3D files 
!
 do nf = 1, mxinst_polarris 
  call dump( atmos_optpol(is:ie,js:je,ks:ke,nf)%zhh  , io, irec )
  call dump( atmos_optpol(is:ie,js:je,ks:ke,nf)%zdr  , io, irec )
  call dump( atmos_optpol(is:ie,js:je,ks:ke,nf)%kdp  , io, irec )
  call dump( atmos_optpol(is:ie,js:je,ks:ke,nf)%rhohv, io, irec )
  call dump( atmos_optpol(is:ie,js:je,ks:ke,nf)%vdop , io, irec )
  call dump( atmos_optpol(is:ie,js:je,ks:ke,nf)%vrad , io, irec )
  call dump( atmos_optpol(is:ie,js:je,ks:ke,nf)%elev , io, irec )
  call dump( atmos_optpol(is:ie,js:je,ks:ke,nf)%vang , io, irec )

 enddo 

if ( polarris_crmout ) then

!
! write binary 3D files 
!
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%t_air , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%rh    , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%u     , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%v     , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%w     , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%rain_rate , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%qc    , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%qr    , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%qi    , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%qs    , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%qg    , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%qh    , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%rer   , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%res   , io, irec )
  call dump( atmos_geopol(is:ie,js:je,ks:ke)%reh   , io, irec )

endif



#if MPI == 2
 if(masterproc) &
#endif
 close(io)

 if(masterproc) print*,''

!
!  write out grads control file
!
#if MPI == 2
 if(masterproc) then
#endif

 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//trim(typ)//&
                  '.'//trim(polarris_axisorient_name)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//trim(typ)//&
                  '.'//trim(polarris_axisorient_name)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  POLARRIS output file '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef   ',ke      ,' linear 0. ',range_polarris
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',irec-1

    do nf = 1, mxinst_polarris

    write(ff,"(I2.2)") nf ! five digit decimal 
    write(io,*) 'zhh'//ff   ,ke,'  0  reflectivity  [dBZ] ('//trim(name_polarris(nf))//')'
    write(io,*) 'zdr'//ff   ,ke,'  0  differential reflectivity [dB] ('//trim(name_polarris(nf))//')'
    write(io,*) 'kdp'//ff   ,ke,'  0  specific diff phase [deg/km] ('//trim(name_polarris(nf))//')'
    write(io,*) 'rhohv'//ff ,ke,'  0  copolar (hv) correlation coefficient [-] ('//trim(name_polarris(nf))//')'
    write(io,*) 'vdop'//ff  ,ke,'  0  doppler velocity [m/s] ('//trim(name_polarris(nf))//')'
    write(io,*) 'vrad'//ff  ,ke,'  0  radial velocity [m/s] ('//trim(name_polarris(nf))//')'
    write(io,*) 'elev'//ff  ,ke,'  0  elevation angle [deg] ('//trim(name_polarris(nf))//')'
    write(io,*) 'vang'//ff  ,ke,'  0  instrumental view angle [deg] ('//trim(name_polarris(nf))//')'

    enddo


    if ( polarris_crmout ) then
      write(io,*) 't_air'  ,ke,'  0  air temperature [K]'
      write(io,*) 'rh   '  ,ke,'  0  RH [-]'
      write(io,*) 'u   '   ,ke,'  0  u wind [m/s]'
      write(io,*) 'v   '   ,ke,'  0  v wind [m/s]'
      write(io,*) 'w   '   ,ke,'  0  w wind [m/s]'
      write(io,*) 'precip' ,ke,'  0  precipitation rate [mm/hr]'
      write(io,*) 'qc    ' ,ke,'  0  cloud mixing ratio [g/m3]'
      write(io,*) 'qr    ' ,ke,'  0  rain mixing ratio [g/m3]'
      write(io,*) 'qi    ' ,ke,'  0  ice mixing ratio [g/m3]'
      write(io,*) 'qs    ' ,ke,'  0  snow mixing ratio [g/m3]'
      write(io,*) 'qg    ' ,ke,'  0  graupel mixing ratio [g/m3]'
      write(io,*) 'qh    ' ,ke,'  0  hail mixing ratio [g/m3]'
      write(io,*) 'rer   ' ,ke,'  0  rain effective radius [um]'
      write(io,*) 'res   ' ,ke,'  0  snow effective radius [um]'
      write(io,*) 'reh   ' ,ke,'  0  hail effective radius [um]'
    endif

    write(io,*) 'endvars '
    close(io)
 endif

#if MPI == 2
 endif
#endif

 return
 end subroutine write_out_polarris_grads

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_polarris_netcdf
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out POLARRIS 3D outpu in NetCDF format. 
! 
! History:
! 02/2016  Toshi Matsui@NASA GSFC ; Add instrument and orienting angle information in global attribute. 
! 09/2016  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer :: ierr,irec
 integer :: i,j,k,m, nf
 character(len=2) :: ff
 character(len=20) :: typ
 integer,parameter :: io = 111
 integer :: is,ie,js,je,ks,ke,di,dj,dk, ns,ne,dn

!netcdf parameters
 integer :: ncid
 integer :: i_dimid, j_dimid, k_dimid, inst_dimid, char_dimid

 integer, parameter :: ndims3d = 3, ndims2d = 2, ndims1d=1, ndims_inst=2
 integer :: dimids1d(ndims1d), dimids2d(ndims2d), dimids3d(ndims3d), dimids_inst(ndims_inst)

 integer :: id
 integer :: id_rname, id_rfreq, id_rlat, id_rlon, id_relev, id_rcov
 integer :: id_hgt
 integer :: id_lat, id_lon
 integer,dimension(abs_mxinst_polarris) :: id_zhh, id_zdr, id_kdp, id_rhohv, id_vdop, id_vrad, id_elev, id_vang
 integer :: id_t_air, id_rh, id_u, id_v, id_w, id_precip, &
            id_qc, id_qr, id_qi, id_qs, id_qg, id_qh

 integer :: mxlyr_refine
 character(len=20) :: name_temp

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end 
 if( use_crm_level_polarris ) then
     mxlyr_refine = mxlyr             !CRM level (eta streching)
 else
     mxlyr_refine = mxlyr_polarris    !constant height level
 endif
 ks=myk_start ; ke=mxlyr_refine
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1
 ns = 1; ne=mxinst_polarris ; dn = ne-ns+1


!
! write out CRM surface input for grads binarry format
!
  typ = '.POLARRIS'

!
! write out POLARRIS output in grads binarry format
!
  sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//trim(typ)//&
                 '.'//trim(polarris_axisorient_name)//trim(output_suffix)//'.nc'

#if MPI == 2
 if(masterproc)  then
#endif

  print("(a)"),'MSG write_out_polarris_netcdf: write POLARRIS output in NetCDF format. ',&
  ' ->',trim(sdsu_io_file)

! if ( real(di*dj*dk) > 5.e7 ) then  !for big data
!   call check( nf90_create(trim(sdsu_io_file), IOR(4,nf90_64bit_offset) , ncid) )  !64bit output
   call check( nf90_create(trim(sdsu_io_file), NF90_64BIT_OFFSET , ncid) )  !64bit output

! else
!   call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )               !32bit output
! endif


!
! Define 3D dimensions
!
 call check( nf90_def_dim(ncid, 'x' , mxgridx, i_dimid) )
 call check( nf90_def_dim(ncid, 'y' , mxgridy, j_dimid) )
 call check( nf90_def_dim(ncid, 'z' , mxlyr_refine, k_dimid) )
 dimids3d = (/ i_dimid, j_dimid, k_dimid /)
 dimids2d = (/ i_dimid, j_dimid /)
 dimids1d = (/ k_dimid /)

!
! instrumental dimension
!
 call check( nf90_def_dim(ncid, 'mxinst_polarris'  ,dn , inst_dimid) )
 call check( nf90_def_dim(ncid, 'char20'  ,20 , char_dimid) )

 dimids_inst = (/ char_dimid, inst_dimid /)

!
! global attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'Number_of_Instruments', mxinst_polarris) )
! call check( nf90_put_att(ncid, NF90_GLOBAL, 'Instrument name ', name_polarris(1:mxinst_polarris) ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'Frequenices_GHz', freq_polarris(1:mxinst_polarris) ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'Latitude_deg' , lat_polarris(1:mxinst_polarris) ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'Longitude_deg', lon_polarris(1:mxinst_polarris) ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'Elevation_km', elev_polradar(1:mxinst_polarris) ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'Coverage_km', coverage_polarris(1:mxinst_polarris) ) )


!
! global attribute (need to add case statement)
!

 mic_select: select case(trim(type_microphysics))
 case('GEN')      ! general bulk microphsics 

!
! cloud axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_cloud', &
                        ityp_axis_gen(1:num_axis_gen%cloud)%cloud  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_cloud', &
                        axis_ratio_gen(1:num_axis_gen%cloud)%cloud  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_cloud', ityp_dist_gen%cloud  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_cloud', mean_orient_gen%cloud  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_cloud', &
                       stdev_orient_gen(1:num_stdev_gen%cloud)%cloud  ) )


!
! rain axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_rain', &
                        ityp_axis_gen(1:num_axis_gen%rain)%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_rain', &
                        axis_ratio_gen(1:num_axis_gen%rain)%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_rain', ityp_dist_gen%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_rain', mean_orient_gen%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_rain', &
                       stdev_orient_gen(1:num_stdev_gen%rain)%rain  ) )

!
! ice axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_ice', &
                        ityp_axis_gen(1:num_axis_gen%ice)%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_ice', &
                        axis_ratio_gen(1:num_axis_gen%ice)%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_ice', ityp_dist_gen%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_ice', mean_orient_gen%ice  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice', &
                       stdev_orient_gen(1:num_stdev_gen%ice)%ice  ) )

!
! snow axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_snow', &
                        ityp_axis_gen(1:num_axis_gen%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_snow', &
                        axis_ratio_gen(1:num_axis_gen%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_snow', ityp_dist_gen%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_snow', mean_orient_gen%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_gen_snow', &
                       stdev_orient_gen(1:num_stdev_gen%snow)%snow  ) )

!
! graupel axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_graupel', &
                        ityp_axis_gen(1:num_axis_gen%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_graupel', &
                        axis_ratio_gen(1:num_axis_gen%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_graupel', ityp_dist_gen%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_graupel', mean_orient_gen%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_graupel', &
                       stdev_orient_gen(1:num_stdev_gen%graupel)%graupel  ) )

!
! hail axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_gen_hail', &
                        ityp_axis_gen(1:num_axis_gen%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_gen_hail', &
                        axis_ratio_gen(1:num_axis_gen%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_gen_hail', ityp_dist_gen%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_gen_hail', mean_orient_gen%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_hail', &
                       stdev_orient_gen(1:num_stdev_gen%hail)%hail  ) )


 case('RAMS')  !RAMS type

!
! cloud axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_cloud1', &
                        ityp_axis_rams(1:num_axis_rams%cloud1)%cloud1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_cloud1', &
                        axis_ratio_rams(1:num_axis_rams%cloud1)%cloud1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_cloud1', ityp_dist_rams%cloud1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_cloud1', mean_orient_rams%cloud1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_cloud1', &
                       stdev_orient_rams(1:num_stdev_rams%cloud1)%cloud1  ) )

 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_cloud2', &
                        ityp_axis_rams(1:num_axis_rams%cloud2)%cloud2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_cloud2', &
                        axis_ratio_rams(1:num_axis_rams%cloud2)%cloud2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_cloud2', ityp_dist_rams%cloud2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_cloud2', mean_orient_rams%cloud2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_cloud2', &
                       stdev_orient_rams(1:num_stdev_rams%cloud2)%cloud2  ) )

!
! rain axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_rain', &
                        ityp_axis_rams(1:num_axis_rams%rain)%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_rain', &
                        axis_ratio_rams(1:num_axis_rams%rain)%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_rain', ityp_dist_rams%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_rain', mean_orient_rams%rain  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_rain', &
                       stdev_orient_rams(1:num_stdev_rams%rain)%rain  ) )

!
! ice axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_ice1', &
                        ityp_axis_rams(1:num_axis_rams%ice1)%ice1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_ice1', &
                        axis_ratio_rams(1:num_axis_rams%ice1)%ice1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_ice1', ityp_dist_rams%ice1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_ice1', mean_orient_rams%ice1  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice1', &
                       stdev_orient_rams(1:num_stdev_rams%ice1)%ice1  ) )

 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_ice2', &
                        ityp_axis_rams(1:num_axis_rams%ice2)%ice2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_ice2', &
                        axis_ratio_rams(1:num_axis_rams%ice2)%ice2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_ice2', ityp_dist_rams%ice2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_ice2', mean_orient_rams%ice2  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice2', &
                       stdev_orient_rams(1:num_stdev_rams%ice2)%ice2  ) )

!
! snow axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_snow', &
                        ityp_axis_rams(1:num_axis_rams%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_snow', &
                        axis_ratio_rams(1:num_axis_rams%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_snow', ityp_dist_rams%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_snow', mean_orient_rams%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_rams_snow', &
                       stdev_orient_rams(1:num_stdev_rams%snow)%snow  ) )

!
! graupel axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_graupel', &
                        ityp_axis_rams(1:num_axis_rams%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_graupel', &
                        axis_ratio_rams(1:num_axis_rams%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_graupel', ityp_dist_rams%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_graupel', mean_orient_rams%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_graupel', &
                       stdev_orient_rams(1:num_stdev_rams%graupel)%graupel  ) )

!
! hail axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_rams_hail', &
                        ityp_axis_rams(1:num_axis_rams%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_rams_hail', &
                        axis_ratio_rams(1:num_axis_rams%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_rams_hail', ityp_dist_rams%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_rams_hail', mean_orient_rams%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_hail', &
                       stdev_orient_rams(1:num_stdev_rams%hail)%hail  ) )



 case('SBM') !HUCM SBM type

!
! liquid axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_liq', &
                        ityp_axis_sbm(1:num_axis_sbm%liq)%liq  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_liq', &
                        axis_ratio_sbm(1:num_axis_sbm%liq)%liq  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_liq', ityp_dist_sbm%liq  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_liq', mean_orient_sbm%liq  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_liq', &
                       stdev_orient_sbm(1:num_stdev_sbm%liq)%liq  ) )


!
! ice column axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_ice_col', &
                        ityp_axis_sbm(1:num_axis_sbm%ice_col)%ice_col  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_ice_col', &
                        axis_ratio_sbm(1:num_axis_sbm%ice_col)%ice_col  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_ice_col', ityp_dist_sbm%ice_col  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_ice_col', mean_orient_sbm%ice_col  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice_col', &
                       stdev_orient_sbm(1:num_stdev_sbm%ice_col)%ice_col  ) )


!
! ice plate axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_ice_pla', &
                        ityp_axis_sbm(1:num_axis_sbm%ice_pla)%ice_pla  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_ice_pla', &
                        axis_ratio_sbm(1:num_axis_sbm%ice_pla)%ice_pla  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_ice_pla', ityp_dist_sbm%ice_pla  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_ice_pla', mean_orient_sbm%ice_pla  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice_pla', &
                       stdev_orient_sbm(1:num_stdev_sbm%ice_pla)%ice_pla  ) )

!
! ice dendrite axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_ice_den', &
                        ityp_axis_sbm(1:num_axis_sbm%ice_den)%ice_den  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_ice_den', &
                        axis_ratio_sbm(1:num_axis_sbm%ice_den)%ice_den  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_ice_den', ityp_dist_sbm%ice_den  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_ice_den', mean_orient_sbm%ice_den  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_ice_den', &
                       stdev_orient_sbm(1:num_stdev_sbm%ice_den)%ice_den  ) )


!
! snow axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_snow', &
                        ityp_axis_sbm(1:num_axis_sbm%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_snow', &
                        axis_ratio_sbm(1:num_axis_sbm%snow)%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_snow', ityp_dist_sbm%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_snow', mean_orient_sbm%snow  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_snow', &
                       stdev_orient_sbm(1:num_stdev_sbm%snow)%snow  ) )

!
! graupel axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_graupel', &
                        ityp_axis_sbm(1:num_axis_sbm%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_graupel', &
                        axis_ratio_sbm(1:num_axis_sbm%graupel)%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_graupel', ityp_dist_sbm%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_graupel', mean_orient_sbm%graupel  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_graupel', &
                       stdev_orient_sbm(1:num_stdev_sbm%graupel)%graupel  ) )

!
! hail axis ratio and orienting angle attribute
!
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_axis_sbm_hail', &
                        ityp_axis_sbm(1:num_axis_sbm%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'axis_ratio_sbm_hail', &
                        axis_ratio_sbm(1:num_axis_sbm%hail)%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'ityp_dist_sbm_hail', ityp_dist_sbm%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'mean_orient_sbm_hail', mean_orient_sbm%hail  ) )
 call check( nf90_put_att(ncid, NF90_GLOBAL, 'stdev_orient_ge_hail', &
                       stdev_orient_sbm(1:num_stdev_sbm%hail)%hail  ) )


 case default
   call stop_sdsu(myrank,'MSG write_out_polarris_netcdf: There is no such cloud_microphysics')
 end select mic_select



!
! instrument parameters
!
 call check( nf90_def_var(ncid, 'name_polarris', NF90_CHAR, dimids_inst, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(mxinst_polarris)') )
 call check( nf90_put_att(ncid, id, 'description',' radar name or ID ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_rname = id

 call check( nf90_def_var(ncid, 'freq_polarris', NF90_REAL, dimids_inst, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(mxinst_polarris)') )
 call check( nf90_put_att(ncid, id, 'description',' radar frequency ') )
 call check( nf90_put_att(ncid, id, 'units', '[GHz]') )
 id_rfreq = id

 call check( nf90_def_var(ncid, 'lat_radar', NF90_REAL, dimids_inst, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(mxinst_polarris)') )
 call check( nf90_put_att(ncid, id, 'description',' radar latitude ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_rlat = id

 call check( nf90_def_var(ncid, 'lon_radar', NF90_REAL, dimids_inst, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(mxinst_polarris)') )
 call check( nf90_put_att(ncid, id, 'description',' radar longitude ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_rlon = id

 call check( nf90_def_var(ncid, 'elev_radar', NF90_REAL, dimids_inst, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(mxinst_polarris)') )
 call check( nf90_put_att(ncid, id, 'description',' radar elevation (above sea level)') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_relev = id

 call check( nf90_def_var(ncid, 'coverage_radar', NF90_REAL, dimids_inst, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(mxinst_polarris)') )
 call check( nf90_put_att(ncid, id, 'description',' radar coverage ') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_rcov = id


!
! 1D parameter ID
!
   if( use_crm_level_polarris ) then
      call check( nf90_def_var(ncid, 'hgt', NF90_REAL, dimids3d, id) )
      call check( nf90_put_att(ncid, id, 'MemoryOrder','(z)') )
      call check( nf90_put_att(ncid, id, 'description',' height above sea level ') )
      call check( nf90_put_att(ncid, id, 'units', '[km]') )
      id_hgt = id
   else
      call check( nf90_def_var(ncid, 'hgt', NF90_REAL, dimids1d, id) )
      call check( nf90_put_att(ncid, id, 'MemoryOrder','(z)') )
      call check( nf90_put_att(ncid, id, 'description',' height above sea level ') )
      call check( nf90_put_att(ncid, id, 'units', '[km]') )
      id_hgt = id
   endif


!
! 2D parameter ID
!
 call check( nf90_def_var(ncid, 'latitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y)') )
 call check( nf90_put_att(ncid, id, 'description',' CRM latitude )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_lat = id

 call check( nf90_def_var(ncid, 'longitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y)') )
 call check( nf90_put_att(ncid, id, 'description',' CRM longitude )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_lon = id


!
! 3D radar parameter ID
!
 do nf = 1, mxinst_polarris

 write(ff,"(I2.2)") nf ! two digit decimal 

 call check( nf90_def_var(ncid, 'zhh'//ff, NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
 call check( nf90_put_att(ncid, id, 'description',' reflectivity at H-polarization ('//trim(name_polarris(nf))//')') )
 call check( nf90_put_att(ncid, id, 'units', '[dBZ]') )
 id_zhh(nf) = id

 call check( nf90_def_var(ncid, 'zdr'//ff, NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
 call check( nf90_put_att(ncid, id, 'description',' differential reflectivity ('//trim(name_polarris(nf))//')') )
 call check( nf90_put_att(ncid, id, 'units', '[dB]') )
 id_zdr(nf) = id

 call check( nf90_def_var(ncid, 'kdp'//ff, NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
 call check( nf90_put_att(ncid, id, 'description',' specific differential phase ('//trim(name_polarris(nf))//')') )
 call check( nf90_put_att(ncid, id, 'units', '[deg/km]') )
 id_kdp(nf) = id

 call check( nf90_def_var(ncid, 'rhohv'//ff, NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
 call check( nf90_put_att(ncid, id, 'description',' copolar (hv) correlation coefficient ('//trim(name_polarris(nf))//')') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_rhohv(nf) = id

 call check( nf90_def_var(ncid, 'vdop'//ff, NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
 call check( nf90_put_att(ncid, id, 'description',' doppler velocity ('//trim(name_polarris(nf))//')') )
 call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
 id_vdop(nf) = id

 call check( nf90_def_var(ncid, 'vrad'//ff, NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
 call check( nf90_put_att(ncid, id, 'description',' radial velocity ('//trim(name_polarris(nf))//')') )
 call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
 id_vrad(nf) = id

 call check( nf90_def_var(ncid, 'elev'//ff, NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
 call check( nf90_put_att(ncid, id, 'description',' elevation angle ('//trim(name_polarris(nf))//')') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_elev(nf) = id

 call check( nf90_def_var(ncid, 'vang'//ff, NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
 call check( nf90_put_att(ncid, id, 'description',' instrumental view angle ('//trim(name_polarris(nf))//')') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_vang(nf) = id

 enddo

!
! CRM 3D parameter ID
!
 CRMOUT_IF: if( polarris_crmout )  then

   call check( nf90_def_var(ncid, 't_air', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' air temperature )') )
   call check( nf90_put_att(ncid, id, 'units', '[degK]') )
   id_t_air = id

   call check( nf90_def_var(ncid, 'rh', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' relative humidity )') )
   call check( nf90_put_att(ncid, id, 'units', '[-]') )
   id_rh = id

   call check( nf90_def_var(ncid, 'u', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' u wind )') )
   call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
   id_u = id

   call check( nf90_def_var(ncid, 'v', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' v wind )') )
   call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
   id_v = id

   call check( nf90_def_var(ncid, 'w', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' w wind )') )
   call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
   id_w = id

   call check( nf90_def_var(ncid, 'precip', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' precipitation )') )
   call check( nf90_put_att(ncid, id, 'units', '[mm/hr]') )
   id_precip = id

   call check( nf90_def_var(ncid, 'qc', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' cloud mass concentration )') )
   call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
   id_qc = id

   call check( nf90_def_var(ncid, 'qr', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' rain mass concentration )') )
   call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
   id_qr = id

   call check( nf90_def_var(ncid, 'qi', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' ice mass concentration )') )
   call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
   id_qi = id

   call check( nf90_def_var(ncid, 'qs', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' snow mass concentration )') )
   call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
   id_qs = id

   call check( nf90_def_var(ncid, 'qg', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' graupel mass concentration )') )
   call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
   id_qg = id

   call check( nf90_def_var(ncid, 'qh', NF90_REAL, dimids3d, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(x-y-z)') )
   call check( nf90_put_att(ncid, id, 'description',' hail mass concentration )') )
   call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
   id_qh = id

 endif CRMOUT_IF

 call check( nf90_enddef(ncid) )  ! End define mode.

#if MPI == 2
 endif
#endif

!
! dump radar information
!

#if MPI == 2
 if(masterproc) then
#endif
     do nf = 1, mxinst_polarris
       call check( nf90_put_var(ncid, id_rname, name_polarris(nf) , start=(/1, nf/), count=(/20,1/) ) ) 
     enddo
     call check( nf90_put_var(ncid, id_rfreq, freq_polarris    (ns:ne)) )
     call check( nf90_put_var(ncid, id_rlat , lat_polarris     (ns:ne)) )
     call check( nf90_put_var(ncid, id_rlon , lon_polarris     (ns:ne)) )
     call check( nf90_put_var(ncid, id_relev, elev_polradar    (ns:ne)) )
     call check( nf90_put_var(ncid, id_rcov , coverage_polarris(ns:ne)) )
#if MPI == 2
 endif
#endif

!
! dump out 1D binary data
!

#if MPI == 2
 if(masterproc) then
#endif

   if( use_crm_level_polarris ) then
      call check( nf90_put_var(ncid, id_hgt, atmos(is:ie,js:je,ks:ke)%hgt ))
   else
      call check( nf90_put_var(ncid, id_hgt, hgt_polarris(ks:ke)) )
   endif

#if MPI == 2
 endif
#endif

!
! dump out 2D binary data
!
 call dumpnc( surface(is:ie,js:je)%lat, id_lat, ncid )
 call dumpnc( surface(is:ie,js:je)%lon, id_lon, ncid )

!
! dump out 3D binary data
!
 do nf = 1, mxinst_polarris

    call dumpnc( atmos_optpol(is:ie,js:je,ks:ke,nf)%zhh, id_zhh(nf), ncid )
    call dumpnc( atmos_optpol(is:ie,js:je,ks:ke,nf)%zdr, id_zdr(nf), ncid )
    call dumpnc( atmos_optpol(is:ie,js:je,ks:ke,nf)%kdp, id_kdp(nf), ncid )
    call dumpnc( atmos_optpol(is:ie,js:je,ks:ke,nf)%rhohv, id_rhohv(nf), ncid )
    call dumpnc( atmos_optpol(is:ie,js:je,ks:ke,nf)%vdop, id_vdop(nf), ncid )
    call dumpnc( atmos_optpol(is:ie,js:je,ks:ke,nf)%vrad, id_vrad(nf), ncid )
    call dumpnc( atmos_optpol(is:ie,js:je,ks:ke,nf)%elev, id_elev(nf), ncid )
    call dumpnc( atmos_optpol(is:ie,js:je,ks:ke,nf)%vang, id_vang(nf), ncid )

 enddo

!
! dump out 3D CRM binary data
!
 if( polarris_crmout )  then

  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%t_air , id_t_air, ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%rh    , id_rh   , ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%u     , id_u    , ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%v     , id_v    , ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%w     , id_w    , ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%rain_rate , id_precip, ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%qc    , id_qc   , ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%qr    , id_qr   , ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%qi    , id_qi   , ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%qs    , id_qs   , ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%qg    , id_qg   , ncid )
  call dumpnc( atmos_geopol(is:ie,js:je,ks:ke)%qh    , id_qh   , ncid )

 endif


#if MPI == 2
   if(masterproc) &
#endif
  call check( nf90_close(ncid) )  ! close nc file


 return
 end subroutine write_out_polarris_netcdf

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_driver_polarris

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
