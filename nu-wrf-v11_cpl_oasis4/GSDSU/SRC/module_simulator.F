#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_simulator
 use module_psd
 use module_floating_point
 use module_const
 use netcdf
 use module_mpi 
 use module_thermo
 use module_gocart_interact
 use module_convolution
 use module_gmp_4ice
 use module_thompson_mp
 use module_wdm_mp
 use module_nssl_2mmt 
 use module_p3_param

 implicit none

!--------------------------------------------------------------------------------------------------
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
!  This is an IO interface module for Satellite Data Simulation Unit (SDSU)
!  All the IO processes are controled in this module.
!
! subroutines
!
! History:
! 09/2020  Toshi Matsui@NASA GSFC : Microwave mxpts_re has been increased from 37 to 41
! 08/2016  Toshi Matsui@NASA GSFC : Add POLARRIS option.
! 06/2015  Toshi Matsui@NASA GSFC : Add zonal / diurnal statistics package. 
! 08/2011  Toshi Matsui@NASA GSFC : GOCART input options is added in WRF input. 
! 08/2011  Toshi Matsui@NASA GSFC : GCE-SBM 3D option is included.
! 07/2011  Toshi Matsui@NASA GSFC : SBM moment output routines was added. 
! 07/2011  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 06/2011  Toshi Matsui@NASA GSFC : Morrison two-moment class is added. 
! 06/2011  Toshi Matsui@NASA GSFC : SBM 43bin option is added. 
! 06/2011  Toshi Matsui@NASA GSFC : Range resolution can specify in radar and lidar simulator. 
! 01/2010  Toshi Matsui@NASA GSFC : Memory decomposition for MPI=2 option. 
! 12/2009  Toshi Matsui@NASA GSFC : Module run-parameter constant are replaced by namelist input.
! 12/2009  Toshi Matsui@NASA GSFC : Add global kind control
! 05/2009  Toshi Matsui@NASA GSFC : Add GCE_RAMS options
! 08/2008  Toshi Matsui@NASA GSFC : Add Re for Mie LUT to account 2-moment 
! 06/2008  Toshi Matsui@NASA GSFC : Add MPI options via c-preprocessor
! 04/2008  Toshi Matsui@NASA GSFC : Add 5 different WRF microphysics (DSD) options.
! 02/2008  Toshi Matsui@NASA GSFC : Add 4 different melting options for microwave dielectric functions.
! 01/2008  Toshi Matsui@NASA GSFC : Add Mie diagnostic routine. 
! 01/2008  Toshi Matsui@NASA GSFC : Add Goddard SW/LW scheme (Broadband & ISCCP-like simulators)
! 12/2007  Toshi Matsui@NASA GSFC ; Add MMF input options
! 11/2007  Toshi Matsui@NASA GSFC ; Input files are unified in a NetCDF format. 
! 10/2007  Toshi Matsui@NASA GSFC ; Microwave/radar/visir simulation account FOV
! 09/2007  Toshi Matsui@NASA GSFC ; Add Lidar options. 
! 08/2007  Toshi Matsui@NASA GSFC ; Add GCE SBM options.
! 07/2007  Toshi Matsui@NASA GSFC ; Add dynamic allocation and better-memory coding.  
! 05/2007  Toshi Matsui@NASA GSFC ; adding downwelling tb for output. 
! 03/2007  Toshi Matsui@NASA GSFC ; Initial
!
!
! References:
!
! Matsui, T., J. Santanello, J. J. Shi, W.-K. Tao, D. Wu, C. Peters-Lidard, E. Kemp, M. Chin, 
!   D. Starr, M. Sekiguchi, and F. Aires, (2014): Introducing multisensor satellite radiance-based 
!   evaluation for regional Earth System modeling, Journal of Geophysical Research, 119, 8450-8475, 
!   doi:10.1002/2013JD021424.
!
!-----------------------------------------------------------------------------------------------------

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 public   ! ALL variables and subourtines are accessible in module_simulator.


! #############################################################################
! ##########################  Run-time Parameters #############################
! #############################################################################

!
! Simulator options (simulator_switch)
!
  logical :: micro       ! microwave simulator switch; on when .true.
  logical :: radar       ! radar simulator switch; on  when .true. 
  logical :: visir       ! visble/IR simulator switch; on when .true. 
  logical :: lidar       ! Lidar simulator switch; on when .true.
  logical :: broad       ! Lidar simulator switch; on when .true.
  logical :: GV          ! GV simulator switch; on when .true. 
  logical :: sarta       ! SARTA (AIRS) simulator switch; on when .true. 
  logical :: polarris    ! !POLLLARIS switch; on when .true.

!
! input/output options (io_options)
!
  character(len=500) :: sdsu_dir_sslut   ! Directory path for the Single-Scattering LUTs
  character(len=500) :: sdsu_dir_data    ! Directory path for the various datafiles 
  character(len=500) :: sdsu_io_name     ! Model input-file-list file under sdsu_dir_input 
  logical :: verbose_SDSU    !if true, print out more comments during run. 
  logical :: write_surface   !if true, write out surface opt properties
  logical :: write_opt       !if true, write out optical (single scattering) properties
  logical :: write_CRM3D     !if true, write out CRM 3D input file in GrADS format.
  logical :: write_CRM2D     !if true, write out CRM 2D input file in GrADS format.
  logical :: write_grads_ctl !if true, write grads control files (logical)
  logical :: write_psdmmt    !if true, write out PSD moment value from HUCM SBM
  character(len=50) :: output_suffix

  logical :: write_cftext  =  .false.   ! if true, write climate-format text output for Super Cloud Library (logical) 
  character(len=50) :: cf_init_time  = '2011-11-23_00:00:00' ! if write_cftext is true, you must specify initial time in this format.

  real(sdsu_fps) :: minlat     ! min longitude of lat-lon output [deg]
  real(sdsu_fps) :: maxlat     ! max latitude  of lat-lon output [deg]
  real(sdsu_fps) :: minlon     ! min longitude of lat-lon output [deg] 
  real(sdsu_fps) :: maxlon     ! max longitude of lat-lon output [deg]
  real(sdsu_fps) :: res_latlon ! grid spacing  of lat-lon output [deg]


  character(len=10) :: mask_method   ! 'NONE' does not use any mask
                                     ! 'GRID' use grid number
                                     ! 'LATLON' use latitude and longitude 

  integer :: is_diag               !statistics starting domain i index
  integer :: ie_diag               !statistics ending domain i index
  integer :: js_diag               !statistics starting domain j index
  integer :: je_diag               !statistics ending domain j index

  integer :: minlat_diag           ! min longitude of statistical diagnosis ( -90 ~  90) [deg]
  integer :: maxlat_diag           ! max latitude  of statistical diagnosis ( -90 ~  90) [deg]
  integer :: minlon_diag           ! min longitude of statistical diagnosis (-180 ~ 180) [deg]
  integer :: maxlon_diag           ! max longitude of statistical diagnosis (-180 ~ 180) [deg]

  logical :: diag_crm_column_mean   ! time-averaged column profiles (MPI=0 or 2)
  logical :: diag_crm_spatial_mean  ! time-averaged gridded values (MPI=0 or 2)
  logical :: diag_crm_spatial_mean_daily ! time-averaged gridded values (MPI=0 or 2)
  logical :: diag_crm_spatial_mean_hourly ! time-averaged gridded values (MPI=0 or 2)
  logical :: diag_crm_diurnal_mean  ! time-averaged diurnal cycle (MPI=0 or 2)
  logical :: diag_crm_diurnal_daily ! daily diurnal cycle (MPI=0 or 2)

!
! MMF database optimzation (io_opt_prps) 
!
  real :: dqr_frac, dqs_frac   ! fraction of qr and qs [%]


!
! Input CRM type, dimension, grid size and directries (crm_options)
!                          
  character(len=20)  :: sim_case   ! CRM type
  character(len=300) :: sdsu_dir_input  ! input CRM directory
  character(len=300) :: sdsu_dir_output  ! output directory
  integer  :: mxgridx   ! max grid # in horizontal x direction 
  integer  :: mxgridy   ! max grid # in horizontal y direction  
  integer  :: mxlyr     ! max grid # in vertical direction 
  real(sdsu_fps) :: gridsize  ! horizontal grid spacings [km]

  character(len=20),parameter :: mmf_case ='atrain_1pm'
                                 ! grid selection for MMF
                                 ! 'global' global domain 
                                 ! 'trmm' latitude is within 37deg domain
                                 ! 'atrain_1pm' longitude is 1pm (1 strip of longtitude)
                                 ! 'atrain_1am' latitude is 1am (1strip of longitude)

  character(len=20) :: cloud_microphysics   ! Cloud Microphyiscs Type (character*20)
                             ! GOD: Goddard bulk 1-mmt scheme [Tao et al. 2003] 
                             ! GOD10: Goddard bulk 1-mmt scheme 2010 [Lang et al. 2010]
                             ! GMP4ICE: Goddard bulk 1-mmt 4ICE scheme [Lang et al. 2015]
                             ! LIN: LIN bulk 1-mmt scheme  [Lin et al. ]
                             ! WSM: WRF-Single-Moment 6-Class Scheme [Hong et al. 2004]
                             ! WDM: WRF-Double-Moment 6-Class Scheme [Hong et al. 2010]
                             ! MORR/MORRH: Morrison two-moment scheme graupel/hail option [Morrison et al. 2009]
                             ! THOMPSON: Thompson hybrid moment scheme [Thompson et al. 2008] 
                             ! NSSL2MMT: National Severe Storm Lab 2-moment scheme [Mansel et al. 2010]
                             ! P3: Predicted Particle Properties (P3) scheme [Morrison and Milbrandt 2015]
                             ! RAMS1: RAMS 1-mmt scheme [Walko et al., 1995]
                             ! RAMS2: RAMS 2-mmt scheme [Meyers et al., 1997] 
                             ! HUCM_SBM: HUCM spectra-bin microphysics scheme [Khain et al. 2007]
                             ! HUCM_SBM43: HUCM spectra-bin microphysics 43 bin scheme [Khain et al. 2010] 

 character(len=10) :: type_microphysics  !type of microphysics (GEN, RAMS, or SBM)

 logical :: gsfc_hail        ! gsf hail options 

 logical :: clear_sky_scene  ! if .true., zero out all condensates. 

 logical :: account_aerosol   !if true,  account aerosol particles
 character(len=20) :: aerosol_microphysics ! Aerosol Microphyiscs Type (character*20)
                                ! GOCART: GOCART 1moment bulk aerosol scheme [Chin et al. ] 
                                ! UCD: UC Davis dust scheme

 logical :: nudge_gocart_on   !if true, read global GOCART data (you must prepare seperately)
 
 logical :: uniform_surface  ! if true  (sim_case='GCE' must be always .true.), it assign uniform
                             ! surface character over entire domain. 

!
! Options for single-scattering computation (single_scatter_options)
!

 logical :: lut_micro         ! Particle single-scattering LUT options for micro/radar simulator.
 logical :: lut_visir         ! Particle single-scattering LUT option for visir simulator.
 logical :: lut_polarris      ! Polarimetric Radar Observable LUT option for POLARRIS 
 logical :: lut_replace       ! Replace existing LUT, if you modify single-scattering routines. 

 integer :: ss_opt_micro  ! single scattering calculation options for Microwave/Radar simulator (integer)
                          ! 1 - Mie (fluffy sphere)
                          ! 2 - SCATDB (non-spherical ice crystals)
                          ! 3 - SnowFake (non-spherical ice crystals) , notavailable

 integer :: scatdb_ice_type ! SCATDB ice crystal shape index (if ss_opt_micro = 2, you must choose)
                            ! 0-hexl, 1-hexs, 2-hexb, 3-hexf, 4-hexp, 5-ros3,
                            ! 6-ros4, 7-ros5, 8-ros6, 9-sstr, 10-sden


 integer :: ice_refraction_func !Effective refraction functions for frozen particles for Microwave/Radar simulator
                          ! 1: Oblique Maxwell-Garnett function that assumes ice inclusion within air matrix.
                          ! 2: Oblique Maxwell-Garnett function that assumes air inclusion within ice matrix.
                          ! 3: Effective-Medium function that assumes homogeneous mixing.

 integer :: melt_opt   ! Methods for computing dielectic constant functions of melting ice particles (Microwave/Radar) 
                       ! 0: Does not account melting particles. 
                       ! 1: Maxwell-Garnett function (ice inclusion within water matrix)
                       ! 2: Maxwell-Garnett function (water inclusion within ice matrix)
                       ! 3: Maxwell-Garnett function (average) ; RECOMMENDED
                       ! 4: Effective-Medium function (homogeneous mixing)

 integer :: ss_opt_visir ! Single scattering calculation options for visir/lidar simulator (integer)
                         ! If you change this option, you must set lut_replace = .true.
                         ! 1 - Mie      (default in R-start)
                         ! 2 - Bryan Baum Non-spherical Ice Scattering Database (column, mix, aggregate)
                         ! 3 - MODIS team Non-spherical Ice Scattering Database (not available yet )

 integer :: land_emiss_micro  ! 0: simple scheme
                              ! 1: NESDIS land emissivity model (V1)
                              ! 2: TELSEM emissivity database

 character(len=200) :: input_bf_emiss  ! name input for IR spectrum emissivity for visir simulator


 integer,parameter :: max_chan = 20  ! maximum channel array size
 integer,parameter :: max_chan_radar = 20  ! maximum channel array size

!
! Microwave sensor options (micro_options)
!
 logical :: slant_path_micro  ! if true it account for slant path option.

 character(len=20) :: micro_sensor   !sensor name 
 logical :: ground_micro             !=.true. for ground based; =.false. for satellite based
 integer :: mxfreq_micro             ! The number of microwave-radiometer channels
 real(sdsu_fps) :: ona_angle_micro        ! off-nadir sensor viewing angle [deg]
 real(sdsu_fps) :: freq_micro(max_chan)   ! Channel frequencies [GHz]
 character(len=20) :: nch_micro(max_chan) ! lut character that is consistent to freq
 real(sdsu_fps) :: fov_ct_micro(max_chan) ! Spatial resolution for cross-track FOV     
 real(sdsu_fps) :: fov_dt_micro(max_chan) ! Spatial resolution for down-track FOV
 real(sdsu_fps) :: beamwidth(max_chan)    ! beam width 
 real(sdsu_fps) :: sma_micro              ! semi-major axis of satellite [km]
 integer        :: nf_L2    ! channel of FOV to be used in L2 parameters and deconvolution


! Orbit / Scan Simulator set up. 
 logical :: scan_micro       ! if true, scan / orbit simulator will run. It slow but more accurate. (logical)
 character(len=20) :: satellite_micro     ! available satellite type (character*20)
 character(len=20) :: scan_type_micro     ! available instrumental type 
 character(len=100) :: inpfile_overpass_micro ! input overpassing informations

 real(sdsu_fps),parameter :: mxhgt_micro = 20.0   ! maximum height of measurement (above sea level) [km] (toshii)
 real(sdsu_fps),parameter :: range_micro = 0.25    ! measurement range resolution [km] (toshii)
 integer :: mxlyr_micro
 real(sdsu_fps),allocatable :: hgt_stag_micro(:)
 real(sdsu_fps),allocatable :: hgt_micro(:)

!
! Radar sensor options (radar_options)
!
 logical :: radar_output_L2  !if .true. it also output radar-range-sampled geophysical parameter. 
                             !but output size could be large
 character(len=20) :: radar_sensor  !sensor name 
 logical :: ground_radar       !=.true. for ground-based sensor; =.false. for satellite-based sensor
 integer :: mxfreq_radar       !The number of channels
 real(sdsu_fps) :: min_echo    !minimal_detactable echo [dBZ]  17 & 20 for pre- and post-boost
 real(sdsu_fps) :: inc_angle_radar     !incidence angle [deg] 12.13 is derived from mean of 1/mu (0 ~ 17) 
 real(sdsu_fps) :: k2(max_chan_radar)  !dielectric constant |k^2| defaults (if not known -> -999.)
 real(sdsu_fps) :: freq_radar(max_chan_radar) !Channel frequencies [GHz] 
 character(len=20) :: nch_radar(max_chan_radar)   !lut character that is consistent to freq_radar
 real(sdsu_fps) :: fov_ct_radar(max_chan_radar)   !Spatial resolution for cross-track FOV
 real(sdsu_fps) :: fov_dt_radar(max_chan_radar)   !Spatial resolution for down-track FOV
 real(sdsu_fps) :: beamwidth_radar(max_chan_radar)!beamwidth [deg]
 real(sdsu_fps) :: mxhgt_radar  ! maximum height of measurement (above sea level) [km]
 real(sdsu_fps) :: minhgt_radar ! minimum height of measurement [km]
 real(sdsu_fps) :: range_radar  ! radar measurement range distance [km]
 real(sdsu_fps) :: range_res_radar ! radar measurement range resolution [km]
 logical :: convolution_radar_on 

! Orbit / Scan Simulator set up.
 logical :: scan_radar  ! if true, scan / orbit simulator will run. It slow but more accurate. (logical)
 character(len=20) :: satellite_radar  ! available satellite type (character*20)
 character(len=20) :: scan_type_radar  ! available instrumental type 
 character(len=100) :: inpfile_overpass_radar ! input overpassing informations

! (following parameters will be computed in allocate_all)
 integer :: mxlyr_radar  ! maximum layer for radar simulator 
 real(sdsu_fps),allocatable :: hgt_stag_radar(:) ! radar measurement height at interface [km]
 real(sdsu_fps),allocatable :: dhgt_radar(:)     ! layer depth of radar height  [km]

!
! VisIR sensor options (visir_options)
!
  integer,parameter :: N_stream = 3  ! Stream number of vis-IR simulator (N_stream=1 -> 1x2 = two stream)
                                     ! More stream -> more accurate but more CPU time 
                                     ! I did not put this options in Configure_SDSU.F

  logical,parameter :: lut_grnd_replace = .false. ! If true, it replace lut_grnd.(takes a few min) 
                                                  ! If you changed N_stream, it must be .true.
                                                  ! Default is false.

  character(len=20) :: visir_sensor  ! sensor name 
  real(sdsu_fps) :: znth_slr     ! solar zenith angle [deg] (if -999. coszen depends on model time.)
  real(sdsu_fps) :: znth_obs     ! veiwing zenith angle [deg]
  real(sdsu_fps) :: azmth        ! azimuth angle between the sun and sensor [deg]
  integer :: mxwavel             ! The number of channels
  real(sdsu_fps) :: wavel(max_chan)          ! Channel wavelengths  [micron]
  character(len=20) :: nch_wavel(max_chan)   ! lut charactere that is consistent to wavel
  real(sdsu_fps) :: fov_ct_visir(max_chan)   ! Spatial resolution for cross-track FOV
  real(sdsu_fps) :: fov_dt_visir(max_chan)   ! Spatial resolution for down-track FOV

! Orbit / Scan Simulator set up. 
 logical :: scan_visir       ! if true, scan / orbit simulator will run. It slow but more accurate. (logical)
 character(len=20) :: satellite_visir     ! available satellite type (character*20)
 character(len=20) :: scan_type_visir     ! available instrumental type 
 character(len=100) :: inpfile_overpass_visir ! input overpassing informations

!
! Lidar Simulator options (lidar_options)
!
  character(len=20) :: lidar_sensor  ! sensor name in three character
  logical :: ground_lidar            !=.true. for ground-based sensor; =.false. for satellite-based sensor
  real(sdsu_fps)  :: MS_Correct      ! multiple scttering correction factor
  integer :: mxwavel_lidar           ! The number of channels
  real(sdsu_fps) :: wavel_lidar(max_chan)         ! Channel wavelengths  [micron]
  character(len=20) :: nch_wavel_lidar(max_chan)  ! lut charactere that is consistent to wavel
  logical :: inst_profile_lidar       ! = .true. for instrument-defined profile (must define mxhgt_radar,dght_radar)
  real(sdsu_fps) :: mxhgt_lidar       ! maximum height of measuremenit (above sea level) [km]
  real(sdsu_fps) :: range_lidar       ! radar measurement range resolution [km]

  ! (following parameters will be computed in allocate_all)
  integer :: mxlyr_lidar  ! maximum layer for lidar simulator 
  real(sdsu_fps),allocatable :: hgt_stag_lidar(:) ! lidar measurement height at interface [km]
  real(sdsu_fps),allocatable :: dhgt_lidar(:)     ! layer depth of lidar height  [km]

! Orbit / Scan Simulator set up. 
 logical :: scan_lidar       ! if true, scan / orbit simulator will run. It slow but more accurate. (logical)
 character(len=20) :: satellite_lidar     ! available satellite type (character*20)
 character(len=20) :: scan_type_lidar     ! available instrumental type 
 character(len=100) :: inpfile_overpass_lidar ! input overpassing informations



!
! Broadband simulator (broad_options)
!
  character(len=20) :: broad_scheme  ! broad band radiation scheme 

  logical :: heating_rate  ! write out 3D braodband SW/LW heating rate [K/day]
                           ! in addition to the default energy budget output.
  real(sdsu_fps) :: fov_ct_broad   ! Spatial resolution for cross-track FOV (CERES)
  real(sdsu_fps) :: fov_dt_broad   ! Spatial resolution for down-track FOV  (CERES)

! Orbit / Scan Simulator set up. 
 logical :: scan_broad       ! if true, scan / orbit simulator will run. It slow but more accurate. (logical)
 character(len=20) :: satellite_broad     ! available satellite type (character*20)
 character(len=20) :: scan_type_broad     ! available instrumental type 
 character(len=100) :: inpfile_overpass_broad ! input overpassing informations


!
! GV simulator (GV_options)
!

  logical :: aircraft_on      ! if true, simulate aircraft 2D-probe measurables parameters.
  logical ::  aircraft_ice    ! if true, it only sample ice particle
                              ! if .false. it only sample liquid particle
  logical :: parsivel_on      ! if ture, simulate ground-based Parsivel measurable parameters.
  logical :: parsivel_liq_cutoff ! if True, parisvel parameters are only accounted for diameter less than 8mm. 
  integer :: particle_shape   ! 0 - sphere
                              ! 1 - non sphere (from SnowFake and 2DVD measurements)
  logical :: dump_psd         ! if true, dump output of full 33-bin PSDs.
  logical :: zonal_sampling_on! true - sample particule zone, and create statistical composite for given area. 
                              ! if zonal_sampling_on is true, define sampling zone below 

  real(sdsu_fps) :: minlat_gv   ! minimum latitude [deg] (Aircraf and Parsivel)
  real(sdsu_fps) :: maxlat_gv   ! maximum latitudfe [deg] (Aircraf and Parsivel)
  real(sdsu_fps) :: minlon_gv   ! minimum longitude [deg] (Aircraf and Parsivel)
  real(sdsu_fps) :: maxlon_gv   ! maximum longitude [deg] (Aircraf and Parsivel)                   
  real(sdsu_fps) :: minalt_gv   ! minimum altitude [km]
  real(sdsu_fps) :: maxalt_gv   ! maximum altitude [km]

!
! SARTA (AIRS) simulator (sarta_options)
!
  integer,parameter :: mxchan_airs = 2834   !total airs channels
  logical :: all_channel_airs   ! if true, it simulates default all AIRS (2784) channels.
                                ! if false, you must specify specific channels. 
  integer :: nchan_airs         ! number of airs channels
  integer :: ichan_airs(mxchan_airs) ! specific airs channel ID numbers
  real(sdsu_fps) ::  sol_zen_airs    ! solar zenith angle [deg] (if -999. coszen depends on model time.)
  real(sdsu_fps) :: view_ang_airs    ! veiwing zenith angle [deg]
  character(len=5) :: airs_unit      !'tb' gives output in brightness temperature.
                                     !'rad' gives output in radiance.

!
! POLARRIS ,polarimetric radar simulator (polarris_options)
!
 logical :: generate_lut_polarris   ! if true, it generate new scattering LUTs
                                    ! , which usually requires 1140 CPUs in MPI task.
 integer,save :: polarris_axisorient_typ  !types of axis ratio and orienting angle assumptions.
                                   ! 0:  User manual 
                                   ! 1:  Ryzhkov et al. 2011 JAMC
                                   ! 2:  Putnam et al. 2017 MWR
                                   ! 3:  MC3E 
 character(len=20) :: polarris_axisorient_name  !certain type name 
 
  integer,parameter :: abs_mxinst_polarris = 50  !maximum array for radar insruments.
  logical :: polarris_crmout        ! if true, dump output of CRM interpolated geophysical parameter (logical)
  integer :: mxinst_polarris        ! maximum instrument number 
  character(len=20),dimension(abs_mxinst_polarris) :: name_polarris  ! sensor name in three character
  real(sdsu_fps) :: freq_polarris(abs_mxinst_polarris)     ! polarimetric radar frequencies [GHz]
  real(sdsu_fps) :: lat_polarris(abs_mxinst_polarris)      ! polarimetric radar location (lat) [deg]
  real(sdsu_fps) :: lon_polarris(abs_mxinst_polarris)      ! polarimetric radar location (lat) [deg]
  real(sdsu_fps) :: elev_polradar(abs_mxinst_polarris)     ! polarimetric radar elevation [km]
  real(sdsu_fps) :: coverage_polarris(abs_mxinst_polarris) ! polarimetric radar coverage [km]


  logical :: use_crm_level_polarris ! if true, use CRM levels, otherwise define vertical levels
  real(sdsu_fps) :: mxhgt_polarris  ! maximum height of measuremenit (above sea level) [km]    
  real(sdsu_fps) :: range_polarris  ! radar measurement range resolution [km]
  integer,parameter :: abs_mxelev_polarris = 20        ! Absolute maximum dimension for elevation angle.

  logical :: output_grads_polarris  ! output control
  logical :: output_netcdf_polarris  ! output control

  integer,parameter :: max_chan_polarris = 100  ! maximum channel array size
  character(len=20) :: nch_polarris(max_chan_radar)   !lut character that is consistent to freq_polarris

  integer :: mxlyr_polarris !maximum height
  real(sdsu_fps),allocatable :: hgt_stag_polarris(:) ! radar measurement height at interface [km]
  real(sdsu_fps),allocatable :: dhgt_polarris(:)     ! layer depth of radar height  [km]
  real(sdsu_fps),allocatable :: hgt_polarris(:)      ! layer height   [km]



! ################################################################################
! ###################    Single-Scattering LUT parameter    ######################
! ################################################################################

!
! LUT bins for microwave and radar (as a function of temperature and effective radius)
!
  integer, parameter :: mxpts_temp  = 12  ! LUT dimension for temperature
  real(sdsu_fps),parameter :: pts_temp_w(mxpts_temp) = & ! LUT temperature vector [K] for liquid hydrometeors
                       (/ 255., 260., 265., 270., 275., 280., 285., &
                          290., 295., 300., 305., 310. /)
  real(sdsu_fps), parameter :: pts_temp_i(mxpts_temp) = & ! LUT temperature vector [K] for frozen hydrometeors
                        (/ 180., 220., 230., 240., 245., 250., 255., &
                           260., 265., 270., 275., 280. /)

  integer, parameter :: mxpts_re  = 41  ! LUT dimension for effective radius [micron] 
  real(sdsu_fps), parameter :: pts_re(mxpts_re) = & ! LUT re vector [micron] for frozen hydrometeors
              (/ 1e-0, 2e-0, 3e-0, 4e-0, 5e-0, 6e-0, 7e-0, 8e-0, 9e-0,&
                 1e+1, 2e+1, 3e+1, 4e+1, 5e+1, 6e+1, 7e+1, 8e+1, 9e+1,&
                 1e+2, 2e+2, 3e+2, 4e+2, 5e+2, 6e+2, 7e+2, 8e+2, 9e+2,&
                 1e+3, 1.2e+3, 1.4e+3, 1.6e+3,1.8e+3, &
                       2e+3, 3e+3, 4e+3, 5e+3, 6e+3, 7e+3, 8e+3, 9e+3,&
                 1e+4 /)

!
! LUT bins for polairmetric radar (POLARRIS)
! toshii-debbug
!  integer, parameter :: mxelev_polarris =  1                   ! The number of elevatation angle 
!  real(sdsu_fps),parameter :: elev_polarris(mxelev_polarris) = & ! Elevatation angle [deg]
!                       (/0./)

  integer, parameter :: mxelev_polarris =  91                  ! The number of elevatation angle 
  real(sdsu_fps),parameter :: elev_polarris(mxelev_polarris) = & ! Elevatation angle [deg]
                       (/0., 1., 2., 3., 4., 5., 6., 7., 8., 9., & 
                        10.,11.,12.,13.,14.,15.,16.,17.,18.,19., &
                        20.,21.,22.,23.,24.,25.,26.,27.,28.,29., &
                        30.,31.,32.,33.,34.,35.,36.,37.,38.,39., &
                        40.,41.,42.,43.,44.,45.,46.,47.,48.,49., &
                        50.,51.,52.,53.,54.,55.,56.,57.,58.,59., &
                        60.,61.,62.,63.,64.,65.,66.,67.,68.,69., &
                        70.,71.,72.,73.,74.,75.,76.,77.,78.,79., &
                        80.,81.,82.,83.,84.,85.,86.,87.,88.,89., &
                        90./)

!toshii debbug
!  integer, parameter :: mxpts_temp_pol  = 1  ! LUT dimension for temperature
!  real(sdsu_fps),parameter :: pts_temp_w_pol(mxpts_temp_pol) = & ! LUT temperature vector [K] for liquid hydrometeors
!                       (/ 255. /)
!  real(sdsu_fps),parameter :: pts_temp_i_pol(mxpts_temp_pol) = & ! LUT temperature vector [K] for ice hydrometeors
!                        (/ 180. /)

  integer, parameter :: mxpts_temp_pol  = 12  ! LUT dimension for temperature
  real(sdsu_fps),parameter :: pts_temp_w_pol(mxpts_temp_pol) = & ! LUT temperature vector [K] for liquid hydrometeors
                       (/ 255., 260., 265., 270., 275., 280., 285., &
                          290., 295., 300., 305., 310. /)
  real(sdsu_fps),parameter :: pts_temp_i_pol(mxpts_temp_pol) = & ! LUT temperature vector [K] for ice hydrometeors
                        (/ 180., 220., 230., 240., 245., 250., 255., &
                           260., 265., 270., 275., 280. /)


!toshii debbug
!  integer, parameter :: mxpts_re_pol  = 1  ! LUT dimension for effective radius [micron] 
!  real(sdsu_fps), parameter :: pts_re_small(mxpts_re_pol) = & ! LUT re vector [micron] for frozen hydrometeors
!                (/ 1.0e+1/)
!  real(sdsu_fps), parameter :: pts_re_large(mxpts_re_pol) = & ! LUT re vector [micron] for frozen hydrometeors
!                (/ 1.0e+3 /)

  integer, parameter :: mxpts_re_pol  = 70 !65  ! LUT dimension for effective radius [micron] 
  real(sdsu_fps), parameter :: pts_re_small(mxpts_re_pol) = & ! LUT re vector [micron] for frozen hydrometeors
                (/ 1.0e-3, 1.0e-2, 1.0e-1, & 
                   1.0e+0, 1.2e+0, 1.4e+0, 1.6e+0, 1.8e+0, &
                   2.0e+0, 2.2e+0, 2.4e+0, 2.6e+0, 2.8e+0, &
                   3.0e+0, 3.2e+0, 3.4e+0, 3.6e+0, 3.8e+0, &
                   4.0e+0, 4.2e+0, 4.4e+0, 4.6e+0, 4.8e+0, &
                   5.0e+0, 5.2e+0, 5.4e+0, 5.6e+0, 5.8e+0, &
                   6.0e+0, 6.2e+0, 6.4e+0, 6.6e+0, 6.8e+0, &
                   7.0e+0, 7.2e+0, 7.4e+0, 7.6e+0, 7.8e+0, &
                   8.0e+0, 8.2e+0, 8.4e+0, 8.6e+0, 8.8e+0, &
                   9.0e+0, 9.2e+0, 9.4e+0, 9.6e+0, 9.8e+0, &
                   1.0e+1, 1.1e+1, 1.2e+1, 1.3e+1, 1.4e+1, & 
                   1.5e+1, 1.6e+1, 1.7e+1, 1.8e+1, 1.9e+1, &
                   2.0e+1, 2.1e+1, 2.2e+1, 2.3e+1, 2.4e+1, &
                   2.5e+1, 2.6e+1, 2.7e+1, 2.8e+1, 2.9e+1, &
                   3.0e+1, 1.0e+2 /)

  real(sdsu_fps), parameter :: pts_re_large(mxpts_re_pol) = & ! LUT re vector [micron] for frozen hydrometeors
                 (/ 1.0e-1, 1.0e+0, 1.0e+1, &   
                   1.0e+2, 1.2e+2, 1.4e+2, 1.6e+2, 1.8e+2, &
                   2.0e+2, 2.2e+2, 2.4e+2, 2.6e+2, 2.8e+2, &
                   3.0e+2, 3.2e+2, 3.4e+2, 3.6e+2, 3.8e+2, &
                   4.0e+2, 4.2e+2, 4.4e+2, 4.6e+2, 4.8e+2, &
                   5.0e+2, 5.2e+2, 5.4e+2, 5.6e+2, 5.8e+2, &
                   6.0e+2, 6.2e+2, 6.4e+2, 6.6e+2, 6.8e+2, &
                   7.0e+2, 7.2e+2, 7.4e+2, 7.6e+2, 7.8e+2, &
                   8.0e+2, 8.2e+2, 8.4e+2, 8.6e+2, 8.8e+2, &
                   9.0e+2, 9.2e+2, 9.4e+2, 9.6e+2, 9.8e+2, &
                   1.0e+3, 1.1e+3, 1.2e+3, 1.3e+3, 1.4e+3, &
                   1.5e+3, 1.6e+3, 1.7e+3, 1.8e+3, 1.9e+3, &
                   2.0e+3, 2.1e+3, 2.2e+3, 2.3e+3, 2.4e+3, &
                   2.5e+3, 2.6e+3, 2.7e+3, 2.8e+3, 2.9e+3, &
                   3.0e+3, 1.0e+4 /)


!
! LUT bins for visible-IR and lidar simulators
!
  integer, parameter :: mxpts_re_visir  = 109  ! LUT dimension for effective radius [micron] 
  real(sdsu_fps), parameter :: pts_re_visir(mxpts_re_visir) = & ! LUT re vector [micron] for frozen hydrometeors
              (/ 1.0e-2, 1.5e-2, 2.0e-2, 2.5e-2, 3.0e-2, 3.5e-2, 4.0e-2, 4.5e-2, 5.0e-2, &
                 5.5e-2, 6.0e-2, 6.5e-2, 7.0e-2, 7.5e-2, 8.0e-2, 8.5e-2, 9.0e-2, 9.5e-2, &
                 1.0e-1, 1.5e-1, 2.0e-1, 2.5e-1, 3.0e-1, 3.5e-1, 4.0e-1, 4.5e-1, 5.0e-1, &
                 5.5e-1, 6.0e-1, 6.5e-1, 7.0e-1, 7.5e-1, 8.0e-1, 8.5e-1, 9.0e-1, 9.5e-1, &
                 1.0e-0, 1.5e-0, 2.0e-0, 2.5e-0, 3.0e-0, 3.5e-0, 4.0e-0, 4.5e-0, 5.0e-0, &
                 5.5e-0, 6.0e-0, 6.5e-0, 7.0e-0, 7.5e-0, 8.0e-0, 8.5e-0, 9.0e-0, 9.5e-0, &
                 1.0e+1, 1.5e+1, 2.0e+1, 2.5e+1, 3.0e+1, 3.5e+1, 4.0e+1, 4.5e+1, 5.0e+1, &
                 5.5e+1, 6.0e+1, 6.5e+1, 7.0e+1, 7.5e+1, 8.0e+1, 8.5e+1, 9.0e+1, 9.5e+1, &
                 1.0e+2, 1.5e+2, 2.0e+2, 2.5e+2, 3.0e+2, 3.5e+2, 4.0e+2, 4.5e+2, 5.0e+2, &
                 5.5e+2, 6.0e+2, 6.5e+2, 7.0e+2, 7.5e+2, 8.0e+2, 8.5e+2, 9.0e+2, 9.5e+2, &
                 1.0e+3, 1.5e+3, 2.0e+3, 2.5e+3, 3.0e+3, 3.5e+3, 4.0e+3, 4.5e+3, 5.0e+3, &
                 5.5e+3, 6.0e+3, 6.5e+3, 7.0e+3, 7.5e+3, 8.0e+3, 8.5e+3, 9.0e+3, 9.5e+3, &
                 1e+4 /)


  integer, parameter :: mxpts_rh  = 36       ! LUT dimension for relative humidity [fraction] 
  real(sdsu_fps), parameter :: pts_rh(mxpts_rh) = & ! LUT rh vector for aerosol particles
              (/ 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, &
                 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.81, 0.82, 0.83, &
                 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, &
                 0.94, 0.95, 0.96, 0.97, 0.98, 0.99 /)


!
! Goddard SW/LW radiation broadband information (note actuall minimam LW wavenumber starts from 0cm-1
!                                                 but R-START (visir simulator) account upto 1000micron = 10cm-1) 
!
  integer,parameter :: nband_lw = 10      ! # of LW radiation band
  real(sdsu_fps),parameter :: waven_min(nband_lw) = (/  10.,  340.,  540.,  800.,  980., &
                                            1100., 1215., 1380., 1900.,  540./) !starting wavenumber [cm-1]
  real(sdsu_fps),parameter :: waven_max(nband_lw) = (/ 340.,  540.,  800.,  980., 1100., &
                                            1215., 1380., 1900., 3000.,  620./) !ending wavenumber [cm-1]

  integer,parameter :: nband_sw = 11       !# of SW radiative band
  real(sdsu_fps),parameter :: wavel_min(nband_sw) = (/0.175, 0.225, 0.245, 0.280, 0.295, 0.310, &
                                            0.325, 0.400, 0.700, 1.220, 2.270/) !starting wavelength [micron]
  real(sdsu_fps),parameter :: wavel_max(nband_sw) = (/0.225, 0.280, 0.260, 0.295, 0.310, 0.320, &
                                            0.400, 0.700, 1.220, 2.270, 10.00/) !ending wavelength [micron]

  integer,parameter :: numband_lut = 5   ! sub band number for LUT making (must be greater than 2) 
                                         ! larger values -> more accurate -> need more time to construct LUTs.
                                         ! Once LUTs are maded, small/large numband_lut does not affect CPU time.  
                                         ! (100 is used for LUT of goddard radiatioin for NU-WRF)

  real(sdsu_fps),parameter :: earth_mean_temp = 300. !earth mean temperature for plank-function weight. 

!
! MSTRNX broadband radiation information 
!
  integer,parameter :: nband_aw = 29      ! # of SW and LW radiation band
  real(sdsu_fps),parameter :: waven_mstrn(nband_aw+1) = (/ &
          10,  250,  400, 530,  610,  670,  750,   820, &
          980, 1175, 1225, 1325, 1400, 2000, 2500, 3300, 3800, 4700, 5200,&
          6000, 10000, 12750, 13250,14750, 23000, 30000, 33500, 36000, &
          43500, 50000/)

! ################################################################################
! ##############    Input Model (Atmosphere Land Time) parameter    ##############
! ################################################################################

!
! IO
!
 logical :: bad_sdsu_io   !if input data does not present --> .true.

!
! Time related parameter
!
  real(sdsu_fps) :: sdsu_yyyy   ! year
  real(sdsu_fps) :: sdsu_mm     ! month
  real(sdsu_fps) :: sdsu_dd     ! day
  real(sdsu_fps) :: sdsu_julian ! julian day
  real(sdsu_fps) :: sdsu_hh     ! hour
  real(sdsu_fps) :: sdsu_nn     ! minutes
  real(sdsu_fps) :: sdsu_ss     ! secounds
  real(sdsu_fps) :: sdsu_gmt    ! GMT hour 
  real(sdsu_fps) :: umu         ! Direction cosine (assigned automatically)
  real(sdsu_fps) :: sdsu_solcon ! Broad-band solar constant at TOA [W/m2]

  character(len=15) :: grads_time  !grads time character (e.g., 00:00Z11jul1985
  integer :: efile_len ! enf-file-name length (excluding suffix of sdsu_inp_name)

  integer :: sdsu_nmax_file !maximum file number 
  character(len=100),allocatable,dimension(:) :: sdsu_inp_list ! input file names
  character(len=100) :: sdsu_inp_name   !Input file names
  character(len=200) :: sdsu_io_file  !output pathfile character

  
!
! Atmospere parameters at layer level, dimension(is:ie,js:je,ks:ke)
! ks=1 is bottom of atmosphere, ks=mxlyr is top of atmosphere
 type atmos_parameter 
    real(sdsu_fps) :: dhgt      ! thickness of layer [km] (derived in get_others)
    real(sdsu_fps) :: press     ! layer pressure [hPa] 
    real(sdsu_fps) :: t_air     ! layer air temperature [degK]
    real(sdsu_fps) :: rho_dair  ! dry air density [g/m3]
    real(sdsu_fps) :: rh        ! relative humidity [%]
    real(sdsu_fps) :: qv        ! water vapor mixing ratio [qv]
    real(sdsu_fps) :: sh        ! specific humidity [g/g] (derived in get_others)
    real(sdsu_fps) :: exner     ! exner function [-] (derived in get_others)
    real(sdsu_fps) :: hgt       ! height [km] (derived in get_others)
    real(sdsu_fps) :: latent_heat ! latent heat [K/day] 
    real(sdsu_fps) :: cloud_nuc ! cloud nucleation rate [kg/ks/s] 
    real(sdsu_fps) :: lw_heat   ! LW radiative heat [K/day] 
    real(sdsu_fps) :: sw_heat   ! SW radiative heat [K/day] 
    real(sdsu_fps) :: rain_rate ! rain rate profile [mm/hr]
    real(sdsu_fps) :: tot_cloud ! total cloud condensate (cloud+ice) [g/m3]
    real(sdsu_fps) :: ccn       ! cloud condensation nuclei [#/cm3]
    real(sdsu_fps) :: icn       ! ice nuclei conc [#/cm3]
    real(sdsu_fps) :: u         ! u wind  [m/sec]
    real(sdsu_fps) :: v         ! v wind  [m/sec]
    real(sdsu_fps) :: w         ! w wind (vertical velocity)  [m/sec]
    real(sdsu_fps) :: sctot     ! total space charge density [C/m3]
 end type atmos_parameter
 type ( atmos_parameter ), allocatable, dimension(:,:,:) :: atmos !atmos parameter

!
! Atmosphere parameters at staggered level, dimension(is:ie,js:je,ks-1:ke) 
!
 type atmos_stagger_parameter 
    real(sdsu_fps) :: hgt         ! 3D height at interface [km] 
    real(sdsu_fps) :: press       ! pressure level at interface [hPa]
    real(sdsu_fps) :: t_air       ! air temperature at interface [degK]
    real(sdsu_fps) :: w           ! vertical velocity [m/s]
 end type atmos_stagger_parameter
 type ( atmos_stagger_parameter ), allocatable, dimension(:,:,:) :: atmos_stag !atmos parameter

 real(sdsu_fps),allocatable,dimension(:) ::   &! (0:mxlyr)             !1D 
     hgt_lev       !1D height at interface [km] (for slantpath)

!
! Surface parameters, dimension(is:ie,js:je)
!
 type surface_parameter
    integer :: diag_mask       !1-diagnosisis 0-skip
    integer :: iland           ! 1-land,  2-water
    integer :: igbp_typ        ! IGBP land-cover type (dominant vegetation type )
    real(sdsu_fps) :: cosz     ! cosine of solar zenith angle [-]
    real(sdsu_fps) :: solar_azimuth ! solar azimuth angle [deg]
    real(sdsu_fps) :: t_skin   ! surface skin temperaure [K] 
    real(sdsu_fps) :: t_air    ! surface air temperature [K]
    real(sdsu_fps) :: rain_rate! rainrate_sfc! surface rainfall rate [mm/hr] (diagnostic purpose)
    real(sdsu_fps) :: u10m     ! horizontal wind speed at 10m AGL [m/s]
    real(sdsu_fps) :: u1km     ! u-wind speed at 1km AGL [m/s] (for GPROF)
    real(sdsu_fps) :: v1km     ! v-wind speed at 1km AGL [m/s] (for GPROF)
    real(sdsu_fps) :: lat      ! latitude [deg] 
    real(sdsu_fps) :: lon      ! lon      [deg] 
    real(sdsu_fps) :: frac_veg ! vegetation fraction [%] (optional for WRF)
    real(sdsu_fps) :: albedo   ! surface SW albedo [-]
    real(sdsu_fps) :: h2o_snow ! snow water equivalent [kg m-2] 
    real(sdsu_fps) :: h2o_col  ! column water vapor [kg m-2] 
    real(sdsu_fps) :: t_soil   ! top-soil temperature [K]
    real(sdsu_fps) :: dhgt_snow! snow depth [m]
    real(sdsu_fps) :: h2o_soil ! soil moisture fraction [0-1]
    real(sdsu_fps) :: elev     ! surface elevation [m]
    real(sdsu_fps) :: hgt_0degC ! freezing height [m]
    real(sdsu_fps) :: path_cloud_liq ! cloud liquid water path [kg/m2]
    real(sdsu_fps) :: path_rain_liq  ! rain liquid water path [kg/m2]
    real(sdsu_fps) :: path_cloud_ice ! cloud ice water path [kg/m2]
    real(sdsu_fps) :: path_rain_ice  ! rain ice water path [kg/m2]
    real(sdsu_fps) :: rain_rate_liq ! rainrate_sfc! surface rainfall rate [mm/hr] (diagnostic purpose)
    real(sdsu_fps) :: rain_rate_ice ! rainrate_sfc! surface rainfall rate [mm/hr] (diagnostic purpose)
    real(sdsu_fps) :: toa_lw_up !TOA LW upwelling   radiation [W/m2] 
    real(sdsu_fps) :: toa_lw_dn !TOA LW downwelling radiation [W/m2] 
    real(sdsu_fps) :: toa_sw_up !TOA SW upwelling   radiation [W/m2] 
    real(sdsu_fps) :: toa_sw_dn !TOA SW downwelling radiation [W/m2] 
    real(sdsu_fps) :: boa_lw_up !BOA LW upwelling   radiation [W/m2] 
    real(sdsu_fps) :: boa_lw_dn !BOA LW downwelling radiation [W/m2] 
    real(sdsu_fps) :: boa_sw_up !BOA SW upwelling   radiation [W/m2] 
    real(sdsu_fps) :: boa_sw_dn !BOA SW downwelling radiation [W/m2] 
    real(sdsu_fps) :: cape      !Convective Available Potential Energy [J/kg] 
    real(sdsu_fps) :: ncape     !Normalized Convective Available Potential Energy [J/kg/m] or [m/s2]
    real(sdsu_fps) :: lcl       !lifting condensation level [m]
    real(sdsu_fps) :: wcd       !warm cloud depth [m]
    integer :: iclass    !0-undefined, 1-low, 2-mid, 3-high, 4-convective, 5-cirrus overlapp with shallow
 end type surface_parameter
 type ( surface_parameter ), allocatable, dimension(:,:) :: surface !surface parameter
 type ( surface_parameter ) :: idealized_surface  !uniform idealied surface


!
! domain statistics
!

!
! domain statistics
!
 type domain_statistics
   real(sdsu_fps) :: maxlat  !domain maximum latitudue [deg]
   real(sdsu_fps) :: maxlon  !domain maximum longitude [deg]
   real(sdsu_fps) :: minlat  !domain minimum latitudue [deg]
   real(sdsu_fps) :: minlon  !domain minimum longitude [deg]
 end type domain_statistics
 type (domain_statistics) :: domain_stat


! ################################################################################
! ########################     Microphysics parameters     #######################
! Goddard SDSU define abstract class for each different microphysics scheme, 
! including 
! * General-type bulk one-moment 6-class microphysics (Goddard, Lin, WSM,WDM,Thompson,Morrison, and etc.)
! * RAMS bulk one/two-moment 8-class microphysics
! * HUCM spectra-bin one-moment 7-class microphysics
! * GOCART bulk 14-class aerosol microphysics
! 2D or 3D arrays associated with above scheme will be allocated only if these 
! scheme is selected in Configure_SDSU.F file. This will minimize memory and 
! increas computation speed, while developer needs to create their own abstract types
! and associated single-scattering routines, when new scheme is introduced in the 
! Goddard SDSU. In this case, contact to Toshi Matsui@NASA GSFC. 
! ################################################################################

  real(sdsu_fps), parameter :: r_cld  = 9.      !radius of cloud water [micron] (consistent to mie_clw_bs)
  real(sdsu_fps), parameter :: r_ice  = 20.     !radius of cloud ice [micron] (fixed)

  real(sdsu_fps),parameter :: Nt_c_thom = 100.e6 ! Thompson cloud number concentration [#/m3]

!
! generalized 1- or 2-moment 6-class microphysics parameters with generalized gamma distribution 
!
 type particle_gen !gen paticle type
    real(sdsu_fps):: cloud,rain,ice,snow,graupel,hail
 end type particle_gen

 type particle_gen_r8 !gen paticle type in double precision
    real(sdsu_fpd) :: cloud,rain,ice,snow,graupel,hail
 end type particle_gen_r8

 type particle_gen_i !gen paticle type in integer
    integer :: cloud,rain,ice,snow,graupel,hail
 end type particle_gen_i

 type particle_gen_cmplx !gen paticle type in integer
    complex :: cloud,rain,ice,snow,graupel,hail
 end type particle_gen_cmplx

 type ( particle_gen ), allocatable, dimension(:,:,:) :: & !3D parameters (is:ie,js:je,ks:ke)
     q_gen,  & ! particle mixing ratio [g/m3]
     n_gen,  & ! particle number conc. [#/m3]
     re_gen    ! particle effective radius [micron]

 type ( particle_gen ), allocatable, dimension(:,:) :: & !2D parameters (is:ie,js:je)
     qcol_gen  ! column integrated (equivalent water path) particle amount [kg/m2]

 type ( particle_gen ) :: & !scalar  parameters
     n0_gen     , & ! intercept for exponential PSD [1/m4]
     lambda_gen , & ! lambda
     mu_gen     , & ! dispersion parameter
     rho_gen    , & ! bulk density [kg/m3]
     a_vt       , & ! alpha of Vt  [-]
     b_vt       , & ! beta of Vt   [-]
     vt_max         ! maximum limit of terminal velocity [m/s]

 real(sdsu_fps),allocatable, dimension(:,:,:) :: & ! special parameters for the NSSL2MMT scheme (is:ie,js:je,ks:ke)
    rho3d_graupel    ,& ! 3D graupel effective density [kg/m3]
    rho3d_hail          ! 3D hail effective density [kg/m3]

 real(sdsu_fps),allocatable, dimension(:,:,:) :: & ! special parameters for the P3 scheme (is:ie,js:je,ks:ke)
    p3_mu_cloud ,& ! dispersion parameters of P3 cloud [-]
    p3_frime    ,& ! Rime ice mass fraction [-]
    p3_rhorime     ! Rime ice density [kg/m3]

!
! RAMS 1- and 2-momemtn microphysics parameters
!
 type particle_rams !rams paticle type
    real(sdsu_fps):: cloud1,cloud2,rain,ice1,ice2,snow,graupel,hail
 end type particle_rams

 type particle_rams_r8 !rams paticle type in double precision
    real(sdsu_fpd):: cloud1,cloud2,rain,ice1,ice2,snow,graupel,hail
 end type particle_rams_r8

 type particle_rams_i !rams paticle type in integer
    integer :: cloud1,cloud2,rain,ice1,ice2,snow,graupel,hail
 end type particle_rams_i

 type particle_rams_cmplx !rams paticle type in complex
    complex :: cloud1,cloud2,rain,ice1,ice2,snow,graupel,hail
 end type particle_rams_cmplx


 type ( particle_rams ), allocatable, dimension(:,:,:) :: & !3D parameters (is:ie,js:je,ks:ke)
     q_rams,  & ! particle mixing ratio [g/m3]
     n_rams,  & ! particle number concentration [#/m3]
     re_rams    ! particle effective radius [micron]

 type ( particle_rams ), allocatable, dimension(:,:) :: & !2D parameters (is:ie,js:je)
     qcol_rams  ! column integrated (equivalent water path) particle amount [kg/m2]

 type ( particle_rams ) :: & !scalar  parameters
     gnu_rams  ! shape parameter for RAMS generalized gamma distirbution

 real(sdsu_fps),dimension(7,16) :: rams_dstprms !RAMS habit table

 integer, dimension(1:31,79:100,1:2) :: rams_jhabtab  ! ice mode1 & mode2 habitat index 
                                           ! as a function of T(1~31) and RH(79~100) 
                                           ! 1st dimension: Temperature (1~31)
                                           ! 2nd dimension: RH (79~100)
                                           ! 3rd dimension: 1- ice mode1 ; 2- ice mode2


!
! HUCM Spectra-Bin Microphysics (SBM) parameters
!
 integer :: nbin ! number of particle size bin of SBM (33 or 43)

 type particle_sbm !SBM paticle type
    real(sdsu_fps) :: liq,ice_col,ice_pla,ice_den,snow,graupel,hail
 end type particle_sbm

 type particle_sbm_r8 !SBM paticle type in double precision
    real(sdsu_fpd) :: liq,ice_col,ice_pla,ice_den,snow,graupel,hail
 end type particle_sbm_r8

 type particle_sbm_i !SBM paticle type in integer
    integer :: liq,ice_col,ice_pla,ice_den,snow,graupel,hail
 end type particle_sbm_i

 type particle_sbm_cmplx !SBM paticle type in complex
    complex :: liq,ice_col,ice_pla,ice_den,snow,graupel,hail
 end type particle_sbm_cmplx


 type ( particle_sbm ), allocatable, dimension(:,:,:,:) :: & !4D parameters (is:ie,js:je,ks:ke,1:nbin)
     n_sbm        ! particle # concentration [1/m4]

 real(sdsu_fps), allocatable, dimension(:,:,:,:) :: &
    fmelt_sbm_snow    ,& !mass melting fraction of snow [-]
    fmelt_sbm_graupel ,& !mass melting fraction of graupel [-]
    fmelt_sbm_hail    ,& !mass melting fraction of hail [-]
    frime_sbm_snow       !mass riming fraction of snow  [-]

 type ( particle_sbm ), allocatable, dimension(:,:,:) :: & !3D parameters (is:ie,js:je,ks:ke)
     q_sbm    ,& ! particle mixing ratio [g/m3]
     re_sbm      ! particle effective radius [micron]

 type ( particle_sbm ), allocatable, dimension(:,:) :: & !2D parameters (is:ie,js:je)
     qcol_sbm  ! column integrated (equivalent water path) particle amount [kg/m2]

 type ( particle_sbm ), allocatable, dimension(:) :: & !size bin parameters (1:nbin)
     x_sbm   ,& ! mass per particle [g]
     den_sbm ,& ! density per particle [g/cm3]
     rad_sbm ,& ! radius of particle [cm]
     drad_sbm,& ! d (radius) [cm]
     vt_sbm     ! terminal velocity [cm/s]

 type ( particle_sbm ), allocatable, dimension(:) :: & !size bin parameters (0:nbin)
     brad_sbm  ! boundary of radius bin [cm]

!
! retrieval hydrometeor class 
!
 type particle_retrieval
  real(sdsu_fps) :: press     ! layer pressure [hPa] 
  real(sdsu_fps) :: t_air     ! layer air temperature [degK]
  real(sdsu_fps) :: sh        ! specific humidity [g/g] (derived in get_others)
  real(sdsu_fps) :: w         ! vertical velocity [m/s]
  real(sdsu_fps) :: iwc       ! cloud ice water content (non-raining) [g/m3]
  real(sdsu_fps) :: lwc       ! cloud liquid water content (non-raining) [g/m3]
  real(sdsu_fps) :: rwc       ! rain water content (raining) [g/m3]
  real(sdsu_fps) :: lfrac     ! liquid fraction (raining) [-]
  real(sdsu_fps) :: bulk_den  ! bulk density (raining) [g/m3]
  real(sdsu_fps) :: rain_rate ! rain rate (raining) [mm/hr]
  real(sdsu_fps) :: mmt(0:6)  ! moments [m^i / m3]  (i=0~6)
  real(sdsu_fps) :: hid       ! hydrometeor identification [-]
  real(sdsu_fps) :: hid_conf  ! hydrometeor identification confidence [%]
 end type particle_retrieval

 type (particle_retrieval) , allocatable , dimension(:)     :: L2_radar_1d     ! radar Level2 Parameters
 type (particle_retrieval) , allocatable , dimension(:,:,:) :: L2_radar_refine ! radar Level2 Parameters
 type (particle_retrieval) , allocatable , dimension(:,:,:) :: L2_radar_conv   ! radar Level2 Parameters
 type (particle_retrieval) , allocatable , dimension(:,:,:) :: L2_radar   !DPR Level2 Parameters


!
! CSU Hydrometeor Identification Definition 
!
 logical,parameter :: radar_hid = .true.  !if true, it will estimate HID from radar sub-fraction.
 type particle_hid !SBM paticle type
    real(sdsu_fps) :: dz,rn,cr,ds,ws,vi,ldg,hdg,ha,bd
 end type particle_hid

 type (particle_hid),allocatable :: sback_hid1d(:)         !backscatter for each HID 
 type (particle_hid),allocatable :: sback_hid1d_refine(:)  !backscatter for each HID 


 

!
! GOCART 1-moment aerosol microphysics parameters
!
 type particle_gocart !GOCART paticle type
    real(sdsu_fps) :: so4,blc,ocn,och,ssa,ssc,du1,du2,du3,du4,du5,du6,du7,du8
 end type particle_gocart

 type particle_gocart_r8 !GOCART paticle type
    real(sdsu_fpd) :: so4,blc,ocn,och,ssa,ssc,du1,du2,du3,du4,du5,du6,du7,du8
 end type particle_gocart_r8

 type ( particle_gocart ), allocatable, dimension(:,:,:) :: & !3D parameters (is:ie,js:je,ks:ke)
     q_gocart  ! particle mixing ratio [g/m3]

 integer,parameter :: mxspc_gocart = 14  !total GOCART aerosol species 

  real(8),dimension(1:6,0:99) :: re_gocart, rm_gocart !for SSLUT index
                             ! effective and mode radius [um] for different RH index (0:99)
                             ! 1~6 index is for gocart spc
                             ! 1- "suso" :sulfate and its precursors
                             ! 2- "waso" : water soluble organic carbon
                             ! 3- "soot" : black carbon
                             ! 4- "ssam" : sea salt (accumulation,fine mode)
                             ! 5- "sscm" : sea salt (coarse mode)
                             ! 6- "dust" : dust

 real,parameter :: dust_frac(4)=(/ 0.01053,0.08421,0.25263,0.65263 /) !fraction for fine dust


!
! UCD dust model
!
 type particle_ucd !GOCART paticle type
    real(sdsu_fps) :: du1,du2,du3,du4,du5
 end type particle_ucd

 type particle_ucd_r8 !GOCART paticle type
    real(sdsu_fpd) :: du1,du2,du3,du4,du5
 end type particle_ucd_r8

 type ( particle_ucd ), allocatable, dimension(:,:,:) :: & !3D parameters (is:ie,js:je,ks:ke)
     q_ucd  ! particle mixing ratio [g/m3]

 integer,parameter :: mxspc_ucd = 5  !total species 

!
! microphysics io funcition
!
 logical :: wdm_re_cloud_exist, wdm_re_ice_exist, wdm_re_snow_exist  !IO parameters for WDM scheme


! ################################################################################
! #######################    Output Satellite Sinals   ###########################
! ################################################################################

! microwave tb -> dimension(mxgridx,mxgridy,mxfreq_micro,1:2)
!                 The 4th dimension (1:2) defines horizontal and vertical polarizations in this order.
  real(sdsu_fps),allocatable :: tb_out(:,:,:,:)     ! Microwave brightness temperature [K]
  real(sdsu_fps),allocatable :: tb_out_fov(:,:,:,:) ! tb_out for FOV [K]





! radar Z and PIA -> dimension(mxgridx,mxgridy,mxlyr,mxfreq_radar) & (mxgridx,mxgridy,mxlyr,mxfreq_radar)
 type radar_output_parameter
  real(sdsu_fps) :: Zt     ! non-attenuated equiv. reflectivity factor [Z] = [mm6/m3] = [um3]
  real(sdsu_fps) :: Zm     ! attenuated     equiv. reflectivity factor [Z] = [mm6/m3] = [um3]
  real(sdsu_fps) :: Vdop   ! Doppler velocity [m/s]
  real(sdsu_fps) :: tau    ! path integrated optical depth [-]
 end type radar_output_parameter
 type( radar_output_parameter ),allocatable :: radar_out(:,:,:,:),radar_out_conv(:,:,:,:)

 type radar_output_2d_parameter
  real(sdsu_fps) :: Zt_max   ! composite (column maximum) non.atte. reflectivity [Z] = [mm6/m3] = [um3] 
  real(sdsu_fps) :: Zm_max   ! composite (column maximum) atte      reflectivity [Z] = [mm6/m3] = [um3] 
  real(sdsu_fps) :: Het      ! echo-top height [km] 
  real(sdsu_fps) :: pia      ! 2way path integrated attenuation (dBZt-dBZm)|bottom [dB]
 end type radar_output_2d_parameter
 type( radar_output_2d_parameter ),allocatable :: radar_out_2d_conv(:,:,:)


! vis/IR radiance -> dimension(mxgridx,mxgridy,mxwavel)
  real(sdsu_fps),allocatable :: radiance_out(:,:,:)     !Radiance [W/m2/str/micron] for visible and NIR band
                                            !Brightness temperature [K] for IR band (>10 micron)
  real(sdsu_fps),allocatable :: radiance_out_fov(:,:,:) ! radiance_out for FOV

! lidar -> dimension (mxgridx,mxgridy,mxlyer,mxfreq_radar) 
  real(sdsu_fps),allocatable :: att_B(:,:,:,:)       !total attenuated backscattering coefficient [m-1 str-1] 
  real(sdsu_fps),allocatable :: att_p22_B(:,:,:,:)   !perpendicular component [m-1 str-1] 
  real(sdsu_fps),allocatable :: att_clear_B(:,:,:,:)       !total attenuated backscattering coefficient [m-1 str-1] 
  real(sdsu_fps),allocatable :: lidar_sback(:,:,:,:) !total backscattering coeff
  real(sdsu_fps),allocatable :: lidar_ratio(:,:,:,:) !lidar ratio 
  real(sdsu_fps),allocatable :: optical_depth(:,:,:,:) !optical depth  [-]

! ISCCP simulator  -> dimension(mxgridx,mxgridy)
 real(sdsu_fps),allocatable,dimension(:,:) :: &  
    cld_opt  ,& !column-cloud optical depth [-]
    cldtop_p ,& !cloud-top pressure [hPa]
    aerosol_opt  ,& !column-aerosol optical depth [-]
    aerotop_p    !aerosol-top pressure [hPa]

! Broadband energy budget -> dimension (mxgridx,mxgridy,14)
! Last dimension 14 is 
! 1  TOA SW downwelling flux [W/m2]
! 2  TOA SW upwelling flux   [W/m2]
! 3  Surface SW downwelling flux [W/m2]
! 4  Surface SW upwelling flux [W/m2]
! 5  Net toa SW radiation [W/m2]
! 6  Net surface SW radiation [W/m2]
! 7  Net atmosphere SW radiation [W/m2] 
! 8  TOA LW downwelling flux [W/m2]
! 9  TOA LW upwelling flux   [W/m2]
! 10 Surface LW downwelling flux [W/m2]
! 11 Surface LW upwelling flux [W/m2]
! 12 Net toa LW radiation [W/m2]
! 13 Net surface LW radiation [W/m2]
! 14 Net atmosphere LW radiation [W/m2] 
  real(sdsu_fps),allocatable,dimension(:,:,:) :: ebudget, ebudget_fov 

  real(sdsu_fps),allocatable,dimension(:,:) :: aod_broad,cod_broad !aerosol and cloud column optical depths at visible broadband

! Braodband heating rate -> dimension(mxgridx,mxgridy,mxlyr)
  real(sdsu_fps),allocatable :: sw_heat(:,:,:)  !shortwave heating rate [K/day]
  real(sdsu_fps),allocatable :: lw_heat(:,:,:)  !longwave heating rate [K/day]

!
! SARTA (AIRS) simulator
!
 integer :: nchan_sarta              ! number of airs channels
 integer :: ichan_sarta(mxchan_airs) ! specific airs channel ID numbers

 real(sdsu_fps),allocatable :: sarta_out(:,:,:)  ! AIRS output radiance [Watts/m^2/cm^-1] or brightness temperature [K]
 integer,allocatable        :: sarta_chan_id(:)      ! AIRS channl ID [-] 
 real(sdsu_fps),allocatable :: sarta_chan_freq(:)    ! AIRS channel frequency [cm-1]

!
! temporal array used for dumpting output in MPI domain decomposition 
!
 real(sdsu_fps),allocatable :: out_domain2d(:,:)   !temporal 2d domain (mxgridx,mxgridy)
 real(sdsu_fps),allocatable :: out_domain3d(:,:,:) !temporal 3d domain (mxgridx,mxgridy,mxlyr)


!
! latlon output
!
 integer :: imax_latlon 
 integer :: jmax_latlon
 type latlon_grid_parameter
  real(sdsu_fps) :: var  !real variables
  real(sdsu_fps) :: vsum !sum of real variables
  integer        :: cnt  !count 
  real(sdsu_fps) :: area !area [km2]
 end type latlon_grid_parameter
 type( latlon_grid_parameter ) ,allocatable,save :: latlon_grid(:,:)
 real(sdsu_fps),allocatable,save :: lat_grid(:), lon_grid(:) ! latitude and longitude of latlon grid[deg]

 real(sdsu_fps),allocatable :: lat_domain(:,:)   !domain latitude (mxgridx,mxgridy)  [deg]
 real(sdsu_fps),allocatable :: lon_domain(:,:)   !domain longitude (mxgridx,mxgridy) [deg]
 real(sdsu_fps),allocatable :: elev_domain(:,:)   !domain elevation (mxgridx,mxgridy) [km]


! ################################################################################
! ########################    Miscelaneous parameter    ##########################
! ################################################################################

!
! Lidar parameter
!
  real(sdsu_fps), allocatable,dimension(:,:) :: overpass  !lidar overpass file (1:pass 0:no pass)

!
! MMF parameters  (2x2.5 grid + 64gce grid)
!
  logical,parameter :: mmf_remove_cirrus = .false. ! if .true. it remove cirrus deck from MMF run (cheating). 

  integer ,parameter :: mxgridx_geos = 144  ! GEOS x grid
  integer ,parameter :: mxgridy_geos =  91  ! GEOS y grid
  integer ,parameter :: mxgridx_gce  =  64  ! GCE  x grid
  integer ,parameter :: mxgridy_gce  =   1  ! GCE  y grid
  integer ,parameter :: mxlyr_gce    =  44 !28  ! GCE  z grid

  integer :: ys_geos    ! start GEOS y grid for TRMM domain (partial)
  integer :: ye_geos    ! ending GEOS y grid for TRMM domain (partial)
  integer :: xs_geos    ! start GEOS x grid for target domain (partial)
  integer :: xe_geos    ! ending GEOS x grid for target domain (partial)

!
! PRPS parameters (The Precipitation Retrieval and Profiling Scheme)
!
 integer :: max_sample_gprof  !maximum sampling number for GPROF MMF database
 integer,parameter :: min_scan_angle = -55.  !minimum scan angle [deg]
 integer,parameter :: max_scan_angle =  55.  !maximum scan angle [deg]
 integer :: min_index_angle
 integer :: max_index_angle
 real(sdsu_fps),allocatable :: scan_angle_prps(:) ! off-nadir sensor viewing angle (scan angle) [deg]

 real(sdsu_fps),allocatable :: tb_prps_out(:,:,:,:,:)  ! Microwave brightness temperature for PRPS [K]
 real(sdsu_fps),allocatable :: tb_prps_out_fov(:,:,:,:,:) ! tb_prps_out for FOV [K]

 real(sdsu_fps),allocatable :: gf_prps(:,:,:,:)      ! normalized gain function [-]
 real(sdsu_fps),allocatable :: inc_angle_prps(:,:,:) ! incidence angle for PRPS [deg]

!
! MMF parameters  (1x1.25 grid + 16gce grid)
!
!  integer ,parameter :: mxgridx_geos = 288  ! GEOS x grid
!  integer ,parameter :: mxgridy_geos = 181  ! GEOS y grid
!  integer ,parameter :: mxgridx_gce  =  16  ! GCE  x grid
!  integer ,parameter :: mxgridy_gce  =   1  ! GCE  y grid
!  integer ,parameter :: mxlyr_gce    =  28  ! GCE  z grid
!  integer, parameter :: mm =  51  ! start GEOS y grid for TRMM domain (partial)
!  integer, parameter :: ye_geos = 131  ! start GEOS y grid for TRMM domain (partial)

! MMF parameters  (2x2.5 grid + 16gce grid)
!  integer ,parameter :: mxgridx_geos = 144  ! GEOS x grid
!  integer ,parameter :: mxgridy_geos =  91  ! GEOS y grid
!  integer ,parameter :: mxgridx_gce  =  16  ! GCE  x grid
!  integer ,parameter :: mxgridy_gce  =   1  ! GCE  y grid
!  integer ,parameter :: mxlyr_gce    =  28  ! GCE  z grid
!  integer, parameter :: ys_geos =  26  ! start GEOS y grid for TRMM domain (partial)
!  integer, parameter :: ye_geos =  66  ! start GEOS y grid for TRMM domain (partial)

!
! GEOS5 options
!
 type geos5_grid
  integer :: imax   ! maximum i dimension of GEOS5
  integer :: jmax   ! maximum j dimension of GEOS5
  integer :: kmax   ! maximum k dimension of GEOS5
  integer :: ids, ide, di
  integer :: jds, jde, dj
  integer :: kds, kde, dk
 end type geos5_grid
 type( geos5_grid ) :: geos  !geos5 grid group

 real(sdsu_fpd),allocatable :: lon_geos5(:)    ! full longitude GEOS5 [deg]
 real(sdsu_fpd),allocatable :: lat_geos5(:)    ! full latitude GEOS5 [deg]
 real(sdsu_fpd),allocatable :: lev_geos5(:)    ! full eta vertical level [-]
 real(sdsu_fps) :: dlon_geos5                  ! longtidue grid space [deg]
 real(sdsu_fps) :: dlat_geos5                  ! latitdue grid space [deg]



!
! ISCCP-like simulator parameters
!

 real(sdsu_fps),parameter :: isccp_tau(7) = (/0.1, 1.3, 3.6, 9.4, 23., 60., 200./)  !cloud optical depath boundary
 real(sdsu_fps),parameter :: isccp_pre(8) = (/1100., 800., 680., 560., 440., 310., 180., 0./) !cloud-top pressure boundary

!
! AIRS parameters
!
 integer,parameter :: nw_airs = 1238

!
! satellite overpass input parameters
!
 type overpass_parameter
    real(sdsu_fps) :: lat       ! latitude [deg] 
    real(sdsu_fps) :: lon       ! longitude [deg]
    logical        :: ascend    ! if true, ascending satellite overpass
    character(len=100) :: tag   ! orbi tag 
 end type overpass_parameter

 type ( overpass_parameter ), allocatable, dimension(:) :: &
     overpass_micro ,&! 
     overpass_radar ,&!
     overpass_lidar ,&!
     overpass_visir ,&!
     overpass_broad   !


!
! useful swtich
!
 logical :: activate_zonal_freq_cloudsat_mask  ! activation of zonal_freq_cloudsat_mask



! ################################################################################
! ############################    Module Interface     ###########################
! ################################################################################

 interface dump
    module procedure dump2d
    module procedure dump3d
 end interface dump

 interface dump_cf
    module procedure dump2d_cf
    module procedure dump3d_cf
 end interface dump_cf

 interface dump_gprof
   module procedure dump_gprof_2d
   module procedure dump_gprof_3d
 end interface dump_gprof

 interface dump_prps
   module procedure dump_prps_2d
   module procedure dump_prps_3d
 end interface dump_prps

 interface dumpnc
    module procedure dump2d_netcdf
    module procedure dump3d_netcdf
 end interface dumpnc

 interface bulk_rain
!   module procedure gce_rain
   module procedure gen_rain
 end interface bulk_rain

 interface refine_range_bin
    module procedure refine_range_bin_micro
    module procedure refine_range_bin_lidar
    module procedure refine_range_bin_bundle
    module procedure refine_range_bin_bundle_sbm
 end interface refine_range_bin

 interface Find_NaN_Inf
    module procedure Find_NaN_Inf_Single
    module procedure Find_NaN_Inf_Double
 end interface Find_NaN_Inf

 interface init_microphysics
    module procedure init_gen
    module procedure init_gen_r8
    module procedure init_gen_cmplx
    module procedure init_rams
    module procedure init_rams_r8
    module procedure init_rams_cmplx
    module procedure init_sbm
    module procedure init_sbm_r8
    module procedure init_sbm_cmplx
    module procedure init_hid
 end interface init_microphysics

 interface init_aerosol
    module procedure init_gocart
    module procedure init_gocart_r8
    module procedure init_ucd
    module procedure init_ucd_r8
 end interface init_aerosol


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_config
 implicit none
!--------------------------------------------------------------------------------------------
! Comments: 
!  This subroutine read SDSU configuration file (run-time parameters).
!
! History:
! 09/2017  Toshi Matsui@NASA GSFC : Add new bundled microphysics type.
! 12/2009  Toshi Matsui@NASA GSFC ; Initial   
!        
! References: 
!-----------------------------------------------------------------------------------------------------

 integer,parameter :: io = 101

 namelist /simulator_switch/ micro, radar, visir, lidar, broad, GV, sarta, polarris

 namelist /io_options/ sdsu_dir_sslut, sdsu_dir_data, sdsu_io_name , verbose_SDSU, write_surface, &
                       write_opt, write_CRM3D, write_CRM2D, write_grads_ctl, &
                       write_psdmmt, output_suffix,&
                       minlat, maxlat, minlon, maxlon, res_latlon

 namelist /diag_options/ diag_crm_column_mean, diag_crm_spatial_mean, diag_crm_diurnal_mean, &
                         diag_crm_diurnal_daily, diag_crm_spatial_mean_daily, diag_crm_spatial_mean_hourly, &
                         mask_method, is_diag, ie_diag, js_diag, je_diag, &
                         minlat_diag, maxlat_diag, minlon_diag, maxlon_diag

 namelist /crm_options/ sim_case, sdsu_dir_input, sdsu_dir_output, mxgridx, mxgridy, mxlyr, gridsize, &
                       cloud_microphysics, gsfc_hail, clear_sky_scene, account_aerosol, aerosol_microphysics, &
                       nudge_gocart_on, uniform_surface, idealized_surface

 namelist /single_scatter_options/ lut_micro, lut_visir, lut_polarris, lut_replace, ss_opt_micro, scatdb_ice_type, &
                                   ice_refraction_func, melt_opt, ss_opt_visir, land_emiss_micro, input_bf_emiss

 namelist /micro_options/ slant_path_micro, micro_sensor, ground_micro, ona_angle_micro, &
                          mxfreq_micro, freq_micro, fov_ct_micro, fov_dt_micro,&
                          beamwidth, sma_micro, nf_L2 , &                   
                          scan_micro, satellite_micro, scan_type_micro, &
                          inpfile_overpass_micro

 namelist /radar_options/ radar_output_L2, radar_sensor, ground_radar, mxfreq_radar, min_echo, inc_angle_radar, &
                          k2, freq_radar, fov_ct_radar,  fov_dt_radar,  &
                          mxhgt_radar, range_radar, &
                          scan_radar, satellite_radar, scan_type_radar, &
                          inpfile_overpass_radar

 namelist /visir_options/ visir_sensor, znth_slr, znth_obs, azmth, &
                          mxwavel, wavel, fov_ct_visir, fov_dt_visir, &
                          scan_visir, satellite_visir, scan_type_visir, &
                          inpfile_overpass_visir

 namelist /lidar_options/ lidar_sensor, ground_lidar, MS_Correct, mxwavel_lidar, &
                          wavel_lidar, inst_profile_lidar, mxhgt_lidar, range_lidar, &
                          scan_lidar, satellite_lidar, scan_type_lidar, &
                          inpfile_overpass_lidar

 namelist /broad_options/ broad_scheme, heating_rate, fov_ct_broad, fov_dt_broad, &
                          scan_broad, satellite_broad, scan_type_broad, &
                          inpfile_overpass_broad

 namelist /GV_options/ aircraft_on, aircraft_ice, parsivel_on, parsivel_liq_cutoff, &
                       particle_shape, dump_psd, zonal_sampling_on, &
                       minlat_gv, maxlat_gv, minlon_gv, maxlon_gv, minalt_gv, maxalt_gv

 namelist /sarta_options/ all_channel_airs, nchan_airs, ichan_airs, sol_zen_airs, view_ang_airs, airs_unit

 namelist /polarris_options/ generate_lut_polarris, polarris_axisorient_typ, polarris_crmout, &
                             mxinst_polarris, name_polarris, freq_polarris,  &
                             lat_polarris, lon_polarris, elev_polradar, coverage_polarris, &
                             use_crm_level_polarris, mxhgt_polarris, range_polarris, &
                             output_grads_polarris, output_netcdf_polarris
 namelist /io_opt_prps/ dqr_frac, dqs_frac


!
! Open and read Configure File. 
!
 call open_oldfile(io, 'Configure_SDSU.F' )

 read(io,simulator_switch)
 read(io,diag_options)
 read(io,io_options)
 read(io,crm_options)
 if( sim_case == 'MMF_GPROF' ) read(io,io_opt_prps)
 read(io,single_scatter_options)
 read(io,micro_options)
 read(io,radar_options)
 read(io,visir_options)
 read(io,lidar_options)
 read(io,broad_options)
 read(io,GV_options)
 read(io,sarta_options)
 read(io,polarris_options)


 close(io)

!
! set name for POLARRIS axis ratio and orienting angle assumption talble
!
if( polarris_axisorient_typ == 0 ) polarris_axisorient_name = 'manual'
if( polarris_axisorient_typ == 1 ) polarris_axisorient_name = 'ryzhkov2011'
if( polarris_axisorient_typ == 2 ) polarris_axisorient_name = 'putnam2017'
if( polarris_axisorient_typ == 3 ) polarris_axisorient_name = 'matsui2018'
if( polarris_axisorient_typ == 4 ) polarris_axisorient_name = 'dolan2013'


!
! check configuration parameter
!
 if( sim_case == 'GCE' ) then
     print*,sim_case
     if( .not. uniform_surface ) then
        call stop_sdsu(myrank,'MSG read_config: GCE simulation must specify uniform_surface')
     endif
 endif

!
! turn off LUT options for SBM
!
 if(cloud_microphysics == 'HUCM_SBM'  .or. &
    cloud_microphysics == 'HUCM_SBM43' ) then
    lut_micro    = .false.
    lut_visir    = .false.
 endif

!
! let's make various switch here after
!
 if( radar .and. scan_radar .and. scan_type_radar == 'CPR' ) then
     activate_zonal_freq_cloudsat_mask = .true.
 else
     activate_zonal_freq_cloudsat_mask = .false.
 endif


!
! determine bundled microphysics type
!
 mic_select: select case(trim(cloud_microphysics))
 case('GOD','GOD10','GMP4ICE','LIN','WSM','WDM','MORR','MORRH','THOM','NTUCLR','NSSL2MMT','P3') ! general bulk microphsics 
  type_microphysics = 'GEN'
 case('RAMS1','RAMS2')  !RAMS type
  type_microphysics = 'RAMS'
 case('HUCM_SBM','HUCM_SBM43') !HUCM SBM type
  type_microphysics = 'SBM'
 case default
   call stop_sdsu(myrank,'MSG read_config: There is no such cloud_microphysics')
 end select mic_select




 return
 end subroutine read_config

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine reconfigure_grid
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!    Read GEOS5 (Fortuna) Forecast output in NetCDF4 (HDF5) format.
! 
! History:
!   11/2012  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer :: i,j,k,l  !looping
 integer :: ncid                !# of variables to inquire  
 integer :: varid               !# of variable id


!
! Open netCDF file
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_list(1)) 
! print*,'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! longitude
!
 call check( nf90_inq_dimid(ncid, 'lon',  varid) )
 call check( nf90_inquire_dimension( ncid, varid, len=geos%imax ) )
 allocate( lon_geos5(geos%imax) )
 call check( nf90_get_var(ncid, varid, lon_geos5 ) )  !longitude [deg]
 dlon_geos5 = real( lon_geos5(2) - lon_geos5(1) )     ! delta longitude [deg]

!
! latitude
! 
 call check( nf90_inq_dimid(ncid, 'lat',  varid) )
 call check( nf90_inquire_dimension( ncid, varid, len=geos%jmax ) )
 allocate( lat_geos5(geos%jmax) )               
 call check( nf90_get_var(ncid, varid, lat_geos5 ) )  ! latitude [deg]
 dlat_geos5 = real( lat_geos5(2) - lat_geos5(1) )     ! delta latitude [deg]

!
! eta vertical level
!
 call check( nf90_inq_dimid(ncid, 'lev',  varid) )
 call check( nf90_inquire_dimension( ncid, varid, len=geos%kmax ) )
 allocate( lev_geos5(geos%kmax) )
 call check( nf90_get_var(ncid, varid, lev_geos5 ) )  ! eta vertical level [-]


!
! Focused domain
!
 do i = 1, geos%imax
    if( lon_geos5(i) >= REAL(minlon)) then
        geos%ids = i ; exit
    endif
 enddo
 do i = 1, geos%imax
    if( lon_geos5(i) >= REAL(maxlon) ) then
        geos%ide = i ; exit
    endif
 enddo
 geos%di = geos%ide-geos%ids + 1

 do j = 1, geos%jmax
    if( lat_geos5(j) >= REAL(minlat) ) then
        geos%jds = j ; exit
    endif
 enddo
 do j = 1, geos%jmax
    if( lat_geos5(j) >= REAL(maxlat) ) then
        geos%jde = j ; exit
    endif
 enddo
 geos%dj = geos%jde-geos%jds + 1

 geos%kds = 1 ; geos%kde = geos%kmax
 geos%dk  = geos%kde - geos%kds + 1

!
! Reconfigure domain size
!

  mxgridx   = geos%di    ! max grid # in horizontal x direction  (integer)
  mxgridy   = geos%dj     ! max grid # in horizontal y direction (integer) 
  mxlyr     = geos%dk     ! max grid # in vertical direction  (integer) Toshi- check this lator

  gridsize  = 33.e0   ! horizontal grid spacings [km] (real) toshii change later

!
! close netCDF file
!
 call check( nf90_close(ncid) )

 return
 end subroutine reconfigure_grid

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_scan_sensor
 implicit none
!-----------------------------------------------------------------------------------------------------
! Comments: 
!  This subroutine re-setup simulator run-time parameters for orbit/scan simulator. 
!
! History:
! 12/2009  Toshi Matsui@NASA GSFC ; Initial   
!        
! References: 
!-----------------------------------------------------------------------------------------------------

!
! re-set up microwave simulator run-time parameters
!
 if( micro .and. scan_micro ) then

    select case(trim(scan_type_micro))

    case('GMI_LF')
  micro_sensor = trim(scan_type_micro)   !sensor name (GPM Microwave Imager)
  ground_micro = .false.    !=.true. for ground based; =.false. for satellite based
  mxfreq_micro = 5                 ! The number of microwave-radiometer channels
  freq_micro(1:mxfreq_micro) = (/10.65, 18.7, 23.8, 36.5, 89.0 /)       ! Channel frequencies [GHz]
  fov_ct_micro(1:mxfreq_micro) = (/ 19.4, 11.2, 9.2, 8.6, 4.4   /)     ! Spatial resolution for cross-track FOV
  fov_dt_micro(1:mxfreq_micro) = (/ 32.2, 18.3, 15.0, 14.4, 7.3  /)    ! Spatial resolution for down-track FOV
  beamwidth(1:mxfreq_micro) = (/ 1.75, 1.00, 0.90, 0.90, 0.40 /)     !!beamwidth of half-power [deg]
  nf_L2  = 4  !channel of FOV to be used in L2 parameters 

    case('GMI_LF37')
  micro_sensor = trim(scan_type_micro)   !sensor name (GPM Microwave Imager)
  ground_micro = .false.    !=.true. for ground based; =.false. for satellite based
  mxfreq_micro = 1                 ! The number of microwave-radiometer channels
  freq_micro(1:mxfreq_micro) = (/36.5/)       ! Channel frequencies [GHz]
  fov_ct_micro(1:mxfreq_micro) = (/ 8.6  /)     ! Spatial resolution for cross-track FOV
  fov_dt_micro(1:mxfreq_micro) = (/ 14.4  /)    ! Spatial resolution for down-track FOV
  beamwidth(1:mxfreq_micro) = (/ 0.90 /)     !!beamwidth of half-power [deg]
  nf_L2  = 1  !channel of FOV to be used in L2 parameters 

    case('GMI_HF')
  micro_sensor = trim(scan_type_micro)   !sensor name (GPM Microwave Imager)
  ground_micro = .false.    !=.true. for ground based; =.false. for satellite based
  mxfreq_micro = 5               ! The number of microwave-radiometer channels
  freq_micro(1:mxfreq_micro) = (/166., 176.31, 180.31, 186.31, 190.31  /)          ! Channel frequencies [GHz]
  fov_ct_micro(1:mxfreq_micro) = (/ 4.4, 4.4, 4.4, 4.4, 4.4  /)      ! Spatial resolution for cross-track FOV
  fov_dt_micro(1:mxfreq_micro) = (/ 7.3, 7.3, 7.3, 7.3, 7.3   /)     ! Spatial resolution for down-track FOV
  beamwidth(1:mxfreq_micro) = (/ 0.4, 0.4, 0.4, 0.4, 0.4 /)     !!beamwidth of half-power [deg]
  nf_L2  = 1  !channel of FOV to be used in L2 parameters 

    case('TMI_LF')
  micro_sensor = trim(scan_type_micro)   !sensor name (TRMM Microwave Imager)
  ground_micro = .false.    ! =.true. for ground based; =.false. for satellite based
  mxfreq_micro = 4          ! The number of microwave-radiometer channels
  freq_micro(1:mxfreq_micro) = (/ 10.65, 19.35, 21.3, 37.0  /)                    ! Channel frequencies [GHz]
  fov_ct_micro(1:mxfreq_micro) =  (/36.0, 18.3, 16.5, 9.7  /)               ! Spatial resolution for cross-track FOV
  fov_dt_micro(1:mxfreq_micro) = (/ 59.5, 30.3, 27.2, 16.0   /)              ! Spatial resolution for down-track FOV
  beamwidth(1:mxfreq_micro) = (/ 1.75, 1.00, 0.90, 0.90 /)     !!beamwidth of half-power [deg]
  nf_L2  = 4  !channel of FOV to be used in L2 parameters 

    case('TMI_HF')
  micro_sensor = trim(scan_type_micro)  !sensor name (TRMM Microwave Imager)
  ground_micro = .false.    !=.true. for ground based; =.false. for satellite based
  mxfreq_micro = 1          ! The number of microwave-radiometer channels
  freq_micro(1:mxfreq_micro) =  85.5          ! Channel frequencies [GHz]
  fov_ct_micro(1:mxfreq_micro) =  4.2                 ! Spatial resolution for cross-track FOV
  fov_dt_micro(1:mxfreq_micro) =  6.8                ! Spatial resolution for down-track FOV
  beamwidth(1:mxfreq_micro) = (/ 0.40 /)     !!beamwidth of half-power [deg]
  nf_L2  = 1  !channel of FOV to be used in L2 parameters 

    case('AMSR_E_LF')
  micro_sensor = trim(scan_type_micro) !sensor name (AMSR-E)
  ground_micro = .false.    !=.true. for ground based; =.false. for satellite based
  mxfreq_micro = 5                ! The number of microwave-radiometer channels
  freq_micro(1:mxfreq_micro) = (/  6.925, 10.65, 18.7, 23.8, 36.5 /) ! Channel frequencies [GHz]
  fov_ct_micro(1:mxfreq_micro) = (/ 43.0, 30.0, 16.0, 18.0,  8.0 /)      ! Spatial resolution for cross-track FOV
  fov_dt_micro(1:mxfreq_micro) = (/ 74.0, 51.0, 27.0, 31.0, 14.0 /)     ! Spatial resolution for down-track FOV
  beamwidth(1:mxfreq_micro) = (/ 2.2, 1.4, 0.8, 0.90, 0.40/)     !!beamwidth of half-power [deg]
  nf_L2  = 5  !channel of FOV to be used in L2 parameters 

    case('AMSR_E_HF')
  micro_sensor = trim(scan_type_micro) !sensor name (AMSR-E)
  ground_micro = .false.    !=.true. for ground based; =.false. for satellite based
  mxfreq_micro = 1                ! The number of microwave-radiometer channels
  freq_micro(1:mxfreq_micro)   = (/  89.0  /) ! Channel frequencies [GHz]
  fov_ct_micro(1:mxfreq_micro) = (/ 4.0  /)      ! Spatial resolution for cross-track FOV
  fov_dt_micro(1:mxfreq_micro) = (/ 6.0  /)     ! Spatial resolution for down-track FOV
  beamwidth(1:mxfreq_micro)    = (/ 0.18 /)     !!beamwidth of half-power [deg]
  nf_L2  = 1  !channel of FOV to be used in L2 parameters 

    case('SSMIS')
  micro_sensor = trim(scan_type_micro) !sensor name (SSMIS slected channels)
  ground_micro = .false.    !=.true. for ground based; =.false. for satellite based
  mxfreq_micro = 9                ! The number of microwave-radiometer channels
  freq_micro(1:mxfreq_micro) = (/  19.35 , 22.235, 37.0, 91.65 , &
                                   150.0 , 176.71, 180.31, 186.31, 189.91 /) ! Channel frequencies [GHz]
  fov_ct_micro(1:mxfreq_micro) = (/ 46.5, 46.5, 31.2, 13.2, &
                                    13.2, 13.2, 13.2, 13.2, 13.2  /)      ! Spatial resolution for cross-track FOV
  fov_dt_micro(1:mxfreq_micro) = (/ 73.6, 73.6, 45.0, 15.5, &
                                    15.5, 15.5, 15.5, 15.5, 15.5  /)     ! Spatial resolution for down-track FOV
  beamwidth(1:mxfreq_micro) =    (/ 2.46, 2.46, 1.51, 0.52, &
                                    0.52, 0.52, 0.52, 0.52, 0.52/)     !!beamwidth of half-power [deg]
  nf_L2  = 3  !channel of FOV to be used in L2 parameters 

    case default

      print("(a)"),'MSG init_scan_sensor: You put sensor name = ',trim(scan_type_micro)
      print*,' It is not available yet '
      print*,'Available sensor names are GMI_LF,GMI_HF,TMI_LF,TMI_HF, AMSR_E, SSMIS_LF, SSMIS_HF'
      call stop_sdsu(myrank,'MSG init_scan_sensor: terminate program')

    end select

 endif

!
! re-set up radar simulator run-time parameters
!
 if( radar .and. scan_radar ) then

    select case(trim(scan_type_radar))

    case('SDPR_Ku')
  radar_sensor = trim(scan_type_radar)   !sensor name (Dual-frequency Precipitation Radar)
  ground_radar = .false.    !=.true. for ground-based sensor; =.false. for satellite-based sensor
  mxfreq_radar = 1          ! The number of channels
  min_echo = 17.            ! minimal_detactable echo [dBZ]
  inc_angle_radar = 17.     ! incidence angle [deg]
  k2 = 0.925                ! dielectric constant |k^2| defaults (if not known -> -999.)
  freq_radar = 13.6         ! Channel frequencies [GHz]  
  fov_ct_radar =  0.1       ! Spatial resolution for cross-track FOV
  fov_dt_radar =  0.1       ! Spatial resolution for down-track FOV 
  beamwidth_radar = 0.01    ! beamwidth of half-power [deg]
  mxhgt_radar  = 22.0       ! maximum height of measuremenit (above sea level) [km]    
  minhgt_radar = -5.0       ! minimum height of measuremenit (above sea level) [km]    
  range_radar= 0.125        ! radar measurement range distance [km]
  range_res_radar= 0.25     ! radar measurement range resolution [km]

    case('DPR_Ku')
  radar_sensor = trim(scan_type_radar)   !sensor name (Dual-frequency Precipitation Radar)
  ground_radar = .false.    !=.true. for ground-based sensor; =.false. for satellite-based sensor
  mxfreq_radar = 1          ! The number of channels
  min_echo = 17.            ! minimal_detactable echo [dBZ]
  inc_angle_radar = 17.     ! incidence angle [deg]
  k2 = 0.925                ! dielectric constant |k^2| defaults (if not known -> -999.)
  freq_radar = 13.6         ! Channel frequencies [GHz]  
  fov_ct_radar =  5.0       ! Spatial resolution for cross-track FOV
  fov_dt_radar =  5.0       ! Spatial resolution for down-track FOV 
  beamwidth_radar = 0.71    ! beamwidth of half-power [deg]
  mxhgt_radar  = 22.0       ! maximum height of measuremenit (above sea level) [km]    
  minhgt_radar = -5.0       ! minimum height of measuremenit (above sea level) [km]    
  range_radar= 0.125        ! radar measurement range distance [km]
  range_res_radar= 0.25     ! radar measurement range resolution [km]

    case('DPR_Ka')
  radar_sensor = trim(scan_type_radar)  !sensor name (Dual-frequency Precipitation Radar)
  ground_radar = .false.    !=.true. for ground-based sensor; =.false. for satellite-based sensor
  mxfreq_radar = 1              !The number of channels
  min_echo = 17.                     !minimal_detactable echo [dBZ]
  inc_angle_radar = 17.            !incidence angle [deg]
  k2 = 0.925                 !dielectric constant |k^2| defaults (if not known -> -999.)
  freq_radar = 35.5    !Channel frequencies [GHz]  
  fov_ct_radar =  5.0       ! Spatial resolution for cross-track FOV
  fov_dt_radar =  5.0       ! Spatial resolution for down-track FOV 
  beamwidth_radar = 0.71    ! beamwidth of half-power [deg]
  mxhgt_radar  = 22.0       ! maximum height of measuremenit (above sea level) [km]    
  minhgt_radar = -5.0       ! minimum height of measuremenit (above sea level) [km]    
  range_radar= 0.125        ! radar measurement range distance [km]
  range_res_radar= 0.25     ! radar measurement range resolution [km]

    case('PR')
  radar_sensor = trim(scan_type_radar)       !sensor name (TRMM Precipitation Radar) (character*20)
  ground_radar = .false.    !=.true. for ground-based sensor; =.false. for satellite-based sensor (logical)
  mxfreq_radar = 1          !The number of channels (integer)
  min_echo = 17.            !minimal_detactable echo [dBZ]  17 & 20 for pre- and post-boost (real)
  inc_angle_radar = 17.     !incidence angle [deg] 12.13 is derived from mean of 1/mu (0 ~ 17) (real)
  k2 = 0.925                !Radar constant |k^2| defaults (if not known -> -999.) (real,dimension(mxfreq_radar))
  freq_radar = 13.8         !Channel frequencies [GHz] (real,dimension(mxfreq_radar))
  fov_ct_radar = 5.0        !Spatial resolution for cross-track FOV (real,dimension(mxfreq_radar))
  fov_dt_radar = 5.0        !Spatial resolution for down-track FOV (real,dimension(mxfreq_radar))
  beamwidth_radar = 0.71    ! beamwidth of half-power [deg] ??
  mxhgt_radar  = 15.0       ! maximum height of measuremenit (above sea level) [km]    
  minhgt_radar =  0.0       ! minimum height of measuremenit (above sea level) [km]    
  range_radar = 0.25        ! radar measurement range resolution [km] ??
  range_res_radar= 0.25     ! radar measurement range resolution [km] ??

   case('CPR')
  radar_sensor = trim(scan_type_radar)   ! sensor name (Cloud Profile Radar or CloudSat)
  ground_radar = .false.    ! =.true. for ground-based sensor; =.false. for satellite-based sensor
  mxfreq_radar = 1          ! The number of channels
  min_echo =     -28.           ! minimal_detactable echo [dBZ]
  inc_angle_radar = 0.16   ! incidence angle [deg]
  k2 = 0.75                 ! dielectric constant |k^2| defaults (if not known -> -999.)
  freq_radar = 94.15        ! Channel frequencies [GHz]
  fov_ct_radar =  1.4       ! Spatial resolution for cross-track FOV
  fov_dt_radar =  2.5       ! Spatial resolution for down-track FOV 
  beamwidth_radar = 0.1    ! beamwidth of half-power [deg] ??
  mxhgt_radar  =  20.0      ! maximum height of measuremenit (above sea level) [km]  
  minhgt_radar =   0.0       ! minimum height of measuremenit (above sea level) [km]    
  range_radar  =  0.2398    ! radar measurement range resolution [km]
  range_res_radar= 0.2398   ! radar measurement range resolution [km] ??

 case default

  print*,'MSG init_scan_sensor: You put sensor name = ', trim(scan_type_radar)
  print*,' It is not available yet.'
  print*,'Available sensor names are PR, DPR_Ku, DPR_Ka '
  call stop_sdsu(myrank,'MSG init_scan_sensor: terminate program') 

 end select

 endif

!
! re-set up radar simulator run-time parameters
!
 if( visir .and. scan_visir ) then

    select case(trim(scan_type_visir))

    case('MODIS_IR')

  visir_sensor = trim(scan_type_visir)   ! sensor name  (character*20)
  znth_slr = 0.      ! solar zenith angle [deg] (if -999. coszen depends on model time.) (real)
  znth_obs = 12.13   ! viewing zenith angle [deg] (real)
  azmth    = 0.      ! azimuth angle between the sun and sensor [deg] (real)
  mxwavel = 3        ! The number of channels (real)
                        !WV    Cloud Cloud
  wavel(1:mxwavel) = (/  7.32, 11.0, 12.0  /)    ! Channel wavelengths  [micron] (real,dimension(mxwavel))
  fov_ct_visir(1:mxwavel) =   (/  1.0, 1.0, 1.0  /)  
! Spatial resolution for cross-track FOV (adjusted for PR FOV) (real,dimension(mxwavel))
  fov_dt_visir(1:mxwavel) =   (/  1.0, 1.0, 1.0  /)     ! Spatial resolution for down-track FOV  (adjusted for PR FOV) (real,dimension(mxwavel))
    end select
 endif


!
! re-set up radar simulator run-time parameters for lidar simulator
!
 if( lidar .and. scan_lidar ) then

    select case(trim(scan_type_lidar))

    case('CALIOP')

  lidar_sensor = 'CALIOP'         ! sensor name in three character (character*20)
  ground_lidar = .false.          !=.true. for ground-based sensor; =.false. for satellite-based sensor (logical)
  MS_Correct = 1.0 !0.7                ! multiple scttering correction factor (real)
  mxwavel_lidar = 2               ! The number of channels (integer)
  wavel_lidar(1:mxwavel_lidar)  = (/ 0.532,1.032 /) ! Channel wavelengths  [micron] (real,dimension(mxwavel_lidar))
  inst_profile_lidar = .true.     ! = .true. for instrument-defined profile (must define mxhgt_lidar,range_lidar)
  mxhgt_lidar  =  20.0            ! maximum height of measuremenit (above sea level) [km]  
  range_lidar  =  0.2398          ! lidar measurement range resolution [km]

    end select
 endif




 return
 end subroutine init_scan_sensor

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_input_filename   
 implicit none
!--------------------------------------------------------------------------------------------
! Comments: 
!  This subroutine read CRM file list.
!
! History:
! 01/2008  Toshi Matsui@NASA GSFC ; Initial   
!        
! References: 
!-----------------------------------------------------------------------------------------------------

 integer,parameter :: inp = 111 !IO unit
 integer :: n !looping indice
 integer :: ioer,ierr   !IO stat
 character(len=100) :: dummy  


!
! Open input-file namelist
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_io_name)
!print*, trim(sdsu_io_file)
 open(unit=inp,file=trim(sdsu_io_file),status='old',iostat=ioer)

 if(ioer /= 0) then
    sdsu_io_file = trim(sdsu_io_name)  ! test local inpfile exist or not

    ! try to open file
    open(unit=inp,file=trim(sdsu_io_file),status='old',iostat=ioer)

    if(ioer /= 0) then
        print("(a)"), 'MSG read_input_filename; cannot find ',trim(sdsu_io_file)
        call stop_sdsu(myrank,'MSG read_input_filename: terminate program') 
    endif

 endif

!
! Estimate sdsu_nmax_file  
!
 n = 0
 forever: do
   n = n + 1
   read(unit=inp,FMT='(A100)',END=999) dummy  !Read CRM input file name (global parameters)
 enddo forever
 999 sdsu_nmax_file = n - 1 !sdsu_nmax_file is global parameter
 close(inp)

!
! Memory allocation
!
 if( .not. allocated(sdsu_inp_list) ) then
    allocate( sdsu_inp_list(1:sdsu_nmax_file), stat=ierr ) 
    if (ierr /= 0) call stop_sdsu(myrank,'MSG read_input_filename: allocation error -> Terminate program.')
 endif


!
! Open input-file namelist
!
 open(unit=inp,file=trim(sdsu_io_file),status='old',iostat=ioer)

!
! Read CRM file name and count file number. 
!
 do n = 1, sdsu_nmax_file
   read(unit=inp,FMT='(A100)',END=999) sdsu_inp_list(n)  !Read CRM input file name (global parameters)
 enddo 
 close(inp)


 return
 end subroutine read_input_filename

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine read_overpass(sim_typ)
 implicit none
!--------------------------------------------------------------------------------------------
! Comments: 
!  This subroutine read satellite orbit information for corresponding CRM file input file.
!
! History:
! 01/2011  Toshi Matsui@NASA GSFC ; Initial   
!        
! References: 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: sim_typ
 integer,parameter :: inp = 111 ! IO unit
 integer :: n                   ! looping indice
 integer :: ioer,ierr           ! IO stat
 character(len=200) :: dummy
 character(len=1) :: a_or_d

!
! Memory allocation
!
 select case (trim(sim_typ))
 case ('micro')
   if( .not. allocated(overpass_micro) ) then
      allocate( overpass_micro(1:sdsu_nmax_file), stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG read_overpass: allocation error -> Terminate program.')
   endif
 case('radar')
   if( .not. allocated(overpass_radar) ) then
      allocate( overpass_radar(1:sdsu_nmax_file), stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG read_overpass: allocation error -> Terminate program.')
   endif
 case('visir')
   if( .not. allocated(overpass_visir) ) then
      allocate( overpass_visir(1:sdsu_nmax_file), stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG read_overpass: allocation error -> Terminate program.')
   endif
 case('lidar')
   if( .not. allocated(overpass_lidar) ) then
      allocate( overpass_lidar(1:sdsu_nmax_file), stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG read_overpass: allocation error -> Terminate program.')
   endif
 case('broad')
   if( .not. allocated(overpass_broad) ) then
      allocate( overpass_broad(1:sdsu_nmax_file), stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG read_overpass: allocation error -> Terminate program.')
   endif

 case default ; call stop_sdsu(myrank,'MSG read_overpass: There is no such sim_typ')
 end select

!
! input overpass file name
!
 select case (trim(sim_typ))
 case ('micro')
   sdsu_io_file = trim(inpfile_overpass_micro)  
 case('radar')
   sdsu_io_file = trim(inpfile_overpass_radar)  
 case('visir')
   sdsu_io_file = trim(inpfile_overpass_visir)
 case('lidar')
   sdsu_io_file = trim(inpfile_overpass_lidar)
 case('broad')
   sdsu_io_file = trim(inpfile_overpass_broad)

 case default ; call stop_sdsu(myrank,'MSG read_overpass: There is no such sim_typ')
 end select

 open(unit=inp,file=trim(sdsu_io_file),status='old',iostat=ioer)

 if(ioer /= 0) then
    print("(a)"), 'MSG read_overpass; cannot find ',trim(sdsu_io_file)
    call stop_sdsu(myrank,'MSG read_overpass: terminate program')
 endif


!
! Read satellite overpass information for each corresponding CRM input file 
!
 do n = 1, sdsu_nmax_file

    select case (trim(sim_typ))
    case ('micro')
      read(inp,*,iostat=ioer) dummy, overpass_micro(n)%lat,       overpass_micro(n)%lon, &
                                  a_or_d,  overpass_micro(n)%tag
      if( a_or_d == 'A') then
          overpass_micro(n)%ascend = .true. 
      elseif( a_or_d == 'D') then
          overpass_micro(n)%ascend = .false.
      else
         call stop_sdsu(myrank,'MSG read_overpass: There is no such a_or_d')
      endif

    case('radar')
      read(inp,*,iostat=ioer) dummy, overpass_radar(n)%lat   , overpass_radar(n)%lon, &
                              a_or_d, overpass_radar(n)%tag
      if( a_or_d == 'A') then
          overpass_radar(n)%ascend = .true.
      elseif( a_or_d == 'D') then
          overpass_radar(n)%ascend = .false.
      else
         call stop_sdsu(myrank,'MSG read_overpass: There is no such a_or_d')
      endif

    case('visir')
      read(inp,*,iostat=ioer) dummy, overpass_visir(n)%lat   , overpass_visir(n)%lon, &
                              a_or_d, overpass_visir(n)%tag
      if( a_or_d == 'A') then
          overpass_visir(n)%ascend = .true.
      elseif( a_or_d == 'D') then
          overpass_visir(n)%ascend = .false.
      else
         call stop_sdsu(myrank,'MSG read_overpass: There is no such a_or_d')
      endif

    case('lidar')
      read(inp,*,iostat=ioer) dummy, overpass_lidar(n)%lat   , overpass_lidar(n)%lon, &
                              a_or_d, overpass_lidar(n)%tag
      if( a_or_d == 'A') then
          overpass_lidar(n)%ascend = .true.
      elseif( a_or_d == 'D') then
          overpass_lidar(n)%ascend = .false.
      else
         call stop_sdsu(myrank,'MSG read_overpass: There is no such a_or_d')
      endif

    case('broad')
      read(inp,*,iostat=ioer) dummy, overpass_broad(n)%lat   , overpass_broad(n)%lon, &
                              a_or_d, overpass_broad(n)%tag
      if( a_or_d == 'A') then
          overpass_broad(n)%ascend = .true.
      elseif( a_or_d == 'D') then
          overpass_broad(n)%ascend = .false.
      else
         call stop_sdsu(myrank,'MSG read_overpass: There is no such a_or_d')
      endif

    case default ; call stop_sdsu(myrank,'MSG read_overpass: There is no such sim_typ')
    end select 

    if(ioer /= 0) then
       call stop_sdsu(myrank,'MSG read_overpass; overpass info is strange')
    endif

    if( trim(dummy) /=  trim(sdsu_inp_list(n)) ) then
       call stop_sdsu(myrank,'MSG read_overpass: CRM inpfile and overpass inpfile are not compatible.')
    endif

 enddo
 close(inp)

 return
 end subroutine read_overpass

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_channel_character
 implicit none
!--------------------------------------------------------------------------------------------
!
! Comments: 
!  This subroutine is convert numeric parameters of sensor frequencies/wavelength into 
!  character format. 
!
! History:
! 11/2013  Toshi Matsui@NASA GSFC ; Initial   
!        
! References: 
!-----------------------------------------------------------------------------------------------------

 if(micro ) call convert_channel_character( 'GHz', mxfreq_micro, freq_micro(1:mxfreq_micro), &
                                             nch_micro(1:mxfreq_micro) )   

 if(radar ) call convert_channel_character( 'GHz', mxfreq_radar, freq_radar(1:mxfreq_radar), &
                                             nch_radar(1:mxfreq_radar) )

 if(visir ) call convert_channel_character( 'micron', mxwavel, wavel     (1:mxwavel), &
                                             nch_wavel(1:mxwavel) )

 if(lidar ) call convert_channel_character( 'micron', mxwavel_lidar, wavel_lidar(1:mxwavel_lidar), &
                                            nch_wavel_lidar(1:mxwavel_lidar) )

 if(polarris ) call convert_channel_character( 'GHz', mxinst_polarris, freq_polarris(1:mxinst_polarris), &
                                             nch_polarris(1:mxinst_polarris) )

 end subroutine init_channel_character

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine convert_channel_character( app, nch, nch_float, nch_char )
 implicit none
!--------------------------------------------------------------------------------------------
!
! Comments: 
!  This subroutine is convert numeric parameters of sensor frequencies/wavelength into 
!  character format. 
!
! History:
! 12/2013  Toshi Matsui@NASA GSFC ; Fixed bug (nch --> n)
! 11/2013  Toshi Matsui@NASA GSFC ; Initial   
!        
! References: 
!-----------------------------------------------------------------------------------------------------

 character(len=*),intent(in) :: app  !output append
 integer :: nch !channel number
 real(sdsu_fps),dimension(nch), intent(in) :: nch_float   ! channel freq/wavelength in float 
 character(len=20),dimension(nch),intent(out) :: nch_char ! channel freq/wavelength in character
 integer :: n
 character(len=5) :: iiiii  !integer portion of character
 character(len=5) :: ddddd  !decimal portion of character
 character(len=1) :: char_int1
 integer :: ith_non0
 integer :: i
 character(len=20) :: d_app


 do n = 1, nch

    !
    ! decimal portion
    !
    write(ddddd,"(I5.5)") NINT(1.e5*( nch_float(n) - REAL(INT(nch_float(n) )))) ! five digit decimal 
     ith_non0=1 !initialize
    find_non0: do i = 5, 1, -1
       if( ddddd(i:i) /= '0' ) then
           ith_non0 = i
           exit find_non0
       endif
    enddo find_non0
    d_app = '.'//ddddd(1:ith_non0)//trim(app)

    !
    ! integer portion
    !
    write(iiiii,"(I5.5)") INT( nch_float(n) )  !five digit integer
    ith_non0=1 !initialize
    find_non0_again: do i = 1,5
       if( iiiii(i:i) /= '0' ) then
           ith_non0 = i
           exit find_non0_again
       endif
    enddo find_non0_again


    !
    ! output
    !
    nch_char(n) =  iiiii(ith_non0:5)//trim(d_app)  ! found bug here

 enddo

 end subroutine convert_channel_character

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine allocate_all
  implicit none
!--------------------------------------------------------------------------------------------
!              = Satellite Data Simulation Unit =
!           Dynamic allocation for global parameters
!
! Comments: 
!  This subroutine is allocating memory for module-global parameters.
!  Also initialize the value as undefined. 
!
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Add GCE SBM options.
! 07/2007  Toshi Matsui@NASA GSFC ; Initial   
!        
! References: 
!-----------------------------------------------------------------------------------------------------
 integer :: ierr !error index
 integer :: nw   !# of wavelengh
 integer :: i,j,k, n
 integer :: is,ie,js,je,ks,ke,di,dj,dk !memory domain parameter
 integer :: k_inst, k_slant
 real(sdsu_fps) :: area  !area [km2]
 integer :: delta_ij
!
! simplify loop index name
!
#if MPI == 2
 if(trim(sim_case) == 'WRF' .or. &
    trim(sim_case) == 'GCESBM' .or. &
    trim(sim_case) ==  'GEOS5_forecast' ) then  !memory decomposition is available for WRF case. 
!domain decomp
    is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr
 else
!entire domain
    is=1 ; ie=mxgridx ; js=1 ; je=mxgridy ; ks=1 ; ke=mxlyr
 endif
#else
!entire domain
    is=1 ; ie=mxgridx ; js=1 ; je=mxgridy ; ks=1 ; ke=mxlyr
#endif

 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1

  if(masterproc) print*,''
  if(masterproc) print*,'MSG allocate_all : allocate memory for global parameters.'
  if(masterproc) print*,''

!
! environmental parameters 
!
  allocate( &
            hgt_lev(0:mxlyr)                  ,&
            stat=ierr )


  if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: CRM input allocation error -> Terminate program. ')

  allocate( atmos     (is:ie,js:je,  mxlyr)    ,& !atmosphere layer parameters
            atmos_stag(is:ie,js:je,0:mxlyr)    ,& !atmosphere staggered-level parameters
            surface   (is:ie,js:je        )    ,& !surface parameters
            stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: CRM input allocation error -> Terminate program. ')

  !initialize
  atmos%dhgt = undefined ; atmos%press = undefined ; atmos%t_air = undefined 
  atmos%rh   = undefined ; atmos%sh    = undefined ; atmos%exner = undefined
  atmos%hgt  = undefined ; atmos%w      = undefined; atmos%u = undefined ; atmos%v = undefined
  atmos%tot_cloud = undefined  ; atmos%latent_heat = undefined ; atmos%lw_heat  = undefined
  atmos%sw_heat = undefined ; atmos%rain_rate = undefined  ; atmos%ccn = undefined
  atmos%icn = undefined  ; atmos%sctot = undefined

  atmos_stag%hgt   = undefined ; atmos_stag%press = undefined 
  atmos_stag%t_air = undefined ; atmos_stag%w     = undefined

  surface%diag_mask = undefined_i2
  surface%iland     = undefined_i2 ; surface%igbp_typ = undefined_i2 ; surface%h2o_soil  = undefined
  surface%cosz      = undefined    ; surface%t_skin   = undefined    ; surface%t_air     = undefined
  surface%rain_rate = undefined    ; surface%u10m     = undefined    ; surface%lat       = undefined
  surface%lon       = undefined    ; surface%frac_veg = undefined    ; surface%albedo    = undefined
  surface%h2o_snow  = undefined    ; surface%t_soil   = undefined    ; surface%dhgt_snow = undefined
  surface%elev      = undefined    ; surface%hgt_0degC= undefined    ; surface%path_cloud_liq = undefined
  surface%rain_rate = undefined    ; surface%h2o_col  = undefined    ; surface%dhgt_snow = undefined
  surface%path_rain_liq = undefined  ; surface%path_cloud_ice = undefined 
  surface%path_rain_ice = undefined  ; surface%rain_rate_liq = undefined 
  surface%rain_rate_ice = undefined


!
! Condensate parameters 
!
 mic_select: select case(trim(type_microphysics))

 case('GEN')

  allocate( q_gen     (is:ie,js:je,mxlyr),& ! particle mixing ratio [g/m3]
            re_gen    (is:ie,js:je,mxlyr),& ! particle effective radius [micron]
            qcol_gen  (is:ie,js:je)      ,& ! column integrated (equivalent water path) particle amount [kg/m2]
            stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: General 6-class microphysics allocation error -> Terminate program. ')

  !initialize
  q_gen%cloud = undefined  ; q_gen%rain    = undefined  ; q_gen%ice   = undefined
  q_gen%snow  = undefined  ; q_gen%graupel = undefined  ; q_gen%hail  = undefined

  re_gen%cloud = undefined  ; re_gen%rain    = undefined  ; re_gen%ice   = undefined
  re_gen%snow  = undefined  ; re_gen%graupel = undefined  ; re_gen%hail  = undefined

  qcol_gen%cloud = undefined  ; qcol_gen%rain    = undefined  ; qcol_gen%ice   = undefined
  qcol_gen%snow  = undefined  ; qcol_gen%graupel = undefined  ; qcol_gen%hail  = undefined

  if( trim(cloud_microphysics) == 'MORR'   .or. &
      trim(cloud_microphysics) == 'MORRH'  .or. &
      trim(cloud_microphysics) == 'NSSL2MMT' .or. &
      trim(cloud_microphysics) == 'P3' .or. &
      trim(cloud_microphysics) == 'THOM'   .or. &
      trim(cloud_microphysics) == 'NTUCLR' .or. &
      trim(cloud_microphysics) == 'WDM' ) then  ! only two-moment
    allocate( n_gen     (is:ie,js:je,mxlyr),&   ! particle number conc [#/m3]
              stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: General 6-class microphysics allocation error -> Terminate program. ')

    n_gen%cloud = undefined  ; n_gen%rain    = undefined  ; n_gen%ice   = undefined
    n_gen%snow  = undefined  ; n_gen%graupel = undefined  ; n_gen%hail  = undefined

  endif             

  if( trim(cloud_microphysics) == 'NSSL2MMT' ) then ! 
    allocate( rho3d_graupel     (is:ie,js:je,mxlyr),&   ! 3D graupel effectie density [kg/m3]
              rho3d_hail        (is:ie,js:je,mxlyr),&   ! 3D hail effective density   [kg/m3] 
              stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: NSSL density array allocation error -> Terminate program. ')

    rho3d_graupel = undefined  ; rho3d_hail  = undefined

  endif

  if( trim(cloud_microphysics) == 'P3' ) then ! 
    allocate( p3_mu_cloud (is:ie,js:je,mxlyr),&   ! dispersion parameter for P3 cloud [-]
              p3_frime    (is:ie,js:je,mxlyr),&   ! Rime ice mass fraction [-]
              p3_rhorime  (is:ie,js:je,mxlyr),&   ! Rime ice density [kg/m3]
              stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: P3 density array allocation error -> Terminate program. ')

    p3_mu_cloud = undefined ; p3_frime = undefined  ; p3_rhorime = undefined

  endif


 case('SBM')
!
! Spectra-bin microphysics 
!
  allocate( n_sbm     (is:ie,js:je,mxlyr,nbin),& ! particle # concentration [1/m4]
            fmelt_sbm_snow(is:ie,js:je,mxlyr,nbin),& ! mass melting fraction of snow aggregate [-]
            fmelt_sbm_graupel(is:ie,js:je,mxlyr,nbin),& ! mass melting fraction of graupel [-]
            fmelt_sbm_hail(is:ie,js:je,mxlyr,nbin),& ! mass melting fraction of hail [-]
            frime_sbm_snow(is:ie,js:je,mxlyr,nbin),& ! rime fraction of snow aggregate [-]
            q_sbm     (is:ie,js:je,mxlyr)     ,& ! particle mixing ratio [g/m3]
            re_sbm    (is:ie,js:je,mxlyr)     ,& ! particle effective radius [micron]
            qcol_sbm  (is:ie,js:je)           ,& ! column integrated (equivalent water path) particle amount [kg/m2]
            stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: HUCM_SBM HUCM_SBM43 allocation error -> Terminate program.')
  !initialize
  n_sbm%liq =undefined ; n_sbm%ice_col=undefined ; n_sbm%ice_pla=undefined ; n_sbm%ice_den=undefined
  n_sbm%snow=undefined ; n_sbm%graupel=undefined ; n_sbm%hail   =undefined

!  fmelt_sbm%liq =undefined ; fmelt_sbm%ice_col=undefined ; fmelt_sbm%ice_pla=undefined ; fmelt_sbm%ice_den=undefined
!  fmelt_sbm%snow=undefined ; fmelt_sbm%graupel=undefined ; fmelt_sbm%hail   =undefined

  q_sbm%liq =undefined ; q_sbm%ice_col=undefined ; q_sbm%ice_pla=undefined ; q_sbm%ice_den=undefined
  q_sbm%snow=undefined ; q_sbm%graupel=undefined ; q_sbm%hail   =undefined

  re_sbm%liq =undefined ; re_sbm%ice_col=undefined ; re_sbm%ice_pla=undefined ; re_sbm%ice_den=undefined
  re_sbm%snow=undefined ; re_sbm%graupel=undefined ; re_sbm%hail   =undefined

  qcol_sbm%liq =undefined ; qcol_sbm%ice_col=undefined ; qcol_sbm%ice_pla=undefined ; qcol_sbm%ice_den=undefined
  qcol_sbm%snow=undefined ; qcol_sbm%graupel=undefined ; qcol_sbm%hail   =undefined

  frime_sbm_snow    =undefined ; fmelt_sbm_snow = undefined ; fmelt_sbm_graupel = undefined ; fmelt_sbm_hail = undefined 

 case('RAMS')
!
! RAMS microphyiscs
!
  allocate( q_rams     (is:ie,js:je,mxlyr),& ! particle mixing ratio [g/m3]
            n_rams     (is:ie,js:je,mxlyr),& ! particle number concentration [#/m3]
            re_rams    (is:ie,js:je,mxlyr),& ! particle effective radius [micron]
            qcol_rams  (is:ie,js:je)      ,& ! column integrated (equivalent water path) particle amount [kg/m2]
            stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: RAMS allocation error -> Terminate program.')

  !initialize
  q_rams%cloud1=undefined ; q_rams%cloud2=undefined ; q_rams%rain=undefined    ; q_rams%ice1=undefined 
  q_rams%ice2  =undefined ; q_rams%snow  =undefined ; q_rams%graupel=undefined ; q_rams%hail=undefined

  n_rams%cloud1=undefined ; n_rams%cloud2=undefined ; n_rams%rain=undefined    ; n_rams%ice1=undefined
  n_rams%ice2  =undefined ; n_rams%snow  =undefined ; n_rams%graupel=undefined ; n_rams%hail=undefined

  re_rams%cloud1=undefined ; re_rams%cloud2=undefined ; re_rams%rain=undefined    ; re_rams%ice1=undefined
  re_rams%ice2  =undefined ; re_rams%snow  =undefined ; re_rams%graupel=undefined ; re_rams%hail=undefined

  qcol_rams%cloud1=undefined ; qcol_rams%cloud2=undefined ; qcol_rams%rain=undefined    ; qcol_rams%ice1=undefined
  qcol_rams%ice2  =undefined ; qcol_rams%snow  =undefined ; qcol_rams%graupel=undefined ; qcol_rams%hail=undefined

 case default 
  ! MLM or LIS has no atmos particles. 
 end select mic_select

!
! Aerosol particles (GOCART only so far)
!
 if(account_aerosol) then

    select case(trim(aerosol_microphysics) )
    case('GOCART') 
      allocate( q_gocart (is:ie,js:je,mxlyr),& ! particle mixing ratio [g/m3]
                stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: GOCART allocation error -> Terminate program.')

      !initialize
      q_gocart%so4=undefined ; q_gocart%blc=undefined ; q_gocart%ocn=undefined ; q_gocart%och=undefined
      q_gocart%ssa=undefined ; q_gocart%ssc=undefined ; q_gocart%du1=undefined ; q_gocart%du2=undefined
      q_gocart%du3=undefined ; q_gocart%du4=undefined ; q_gocart%du5=undefined ; q_gocart%du6=undefined
      q_gocart%du7=undefined ; q_gocart%du8=undefined 
    case('UCD')
      allocate( q_ucd (is:ie,js:je,mxlyr),& ! particle mixing ratio [g/m3]
                stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: UCD allocation error -> Terminate program.')

      !initialize
      q_ucd%du1=undefined ; q_ucd%du2=undefined
      q_ucd%du3=undefined ; q_ucd%du4=undefined ; q_ucd%du5=undefined

    case default
      call stop_sdsu(myrank,'MSG allocate_all: There is no such aerosol_microphysics')
    end select
 endif

!
! microwave tb output-> The dimension (1:2) defines horizontal and vertical polarizations in this order.
!
  if(micro) then

! new optionis for slant path
    if( scan_micro .or. slant_path_micro ) then
       mxlyr_micro =   nint( mxhgt_micro/range_micro )  ! new mxlyr with constant range (maximum height ~15km)
       allocate( hgt_stag_micro (0:mxlyr_micro), hgt_micro(1:mxlyr_micro) )
       hgt_stag_micro(0) = 0.  !always zero [km]
       do k_slant = 1, mxlyr_micro
          hgt_stag_micro(k_slant) = real(k_slant) * range_micro   ![km]
          hgt_micro(k_slant) = real(k_slant) * range_micro - 0.5 * range_micro
       enddo
    endif

    if(trim(sim_case) == 'MMF_GPROF') then

      ! for PRPS
      allocate(tb_prps_out(mxgridx,mxgridy,mxfreq_micro,1:2,min_index_angle:max_index_angle), &  !Tb [K]
               tb_prps_out_fov(mxgridx,mxgridy,mxfreq_micro,1:2,min_index_angle:max_index_angle), &
               gf_prps( myi_start:myi_end, myj_start:myj_end, mxfreq_micro, min_index_angle:max_index_angle ), &
               inc_angle_prps( myi_start:myi_end, myj_start:myj_end, min_index_angle:max_index_angle ), &
               stat=ierr ) !Microwave brightness temperature [K]
      if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: PRPS Tb output alocation error -> Terminate program.')
      tb_prps_out = undefined     ; tb_prps_out_fov = undefined
      gf_prps     = undefined     ; inc_angle_prps  = undefined

    else

      allocate(tb_out(mxgridx,mxgridy,mxfreq_micro,1:2), &  !Microwave brightness temperature [K]
               tb_out_fov(mxgridx,mxgridy,mxfreq_micro,1:2), stat=ierr ) !Microwave brightness temperature [K]
      if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: microw output allocation error -> Terminate program.')
      tb_out = undefined     ; tb_out_fov = undefined

    endif


  endif

!
! radar echo output
!
  if(radar) then


     mxlyr_radar = NINT( mxhgt_radar / range_radar ) !radar instrumental maximum height level 
     if(masterproc) print*,'MSG allocateall: radar instrumental max level is =',mxlyr_radar
     allocate( hgt_stag_radar(0:mxlyr_radar) ,&
               dhgt_radar    (1:mxlyr_radar)  &
                )

     ! compute radar range-based height at interface
     hgt_stag_radar(0) = 0.  !sea level 0km
     do k_inst = 1, mxlyr_radar
        hgt_stag_radar(k_inst)  = REAL(k_inst) * range_radar  ![km]
     enddo
     dhgt_radar = range_radar  ![km]
   
     if( radar_hid ) then
        allocate(sback_hid1d       (mxlyr)) !HID
        allocate(sback_hid1d_refine(mxlyr_radar)) !HID
     endif

     !
     ! column simulation (traditional)
     !
     if( .not. scan_radar ) then 

        if( (maxval( fov_ct_radar ) <= gridsize) .or. &
            (mxgridx == 1 .and. mxgridy == 1)  ) then
            convolution_radar_on = .false.
            delta_ij = 0
        else 
           convolution_radar_on = .true.
           !
           ! convolution pixels
           !
           delta_ij =  int(  1.5* ( maxval( fov_ct_radar )  /  gridsize )  )
        endif

        !extra domain for slant path
        call get_extra_memory_domain(delta_ij, is,ie,js,je,ibs,ibe,jbs,jbe)

        allocate(radar_out        (ibs:ibe,jbs:jbe,mxlyr_radar,mxfreq_radar), &
                 radar_out_conv   (ibs:ibe,jbs:jbe,mxlyr_radar,mxfreq_radar), &
                 radar_out_2d_conv(ibs:ibe,jbs:jbe,mxfreq_radar)            , &
                 stat=ierr )   ! 
        if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: radar output allocation error -> Terminate program.')

        !
        ! initialize as undefined
        !	
        radar_out%Zt = undefined ; radar_out%Zm = undefined ; radar_out%Vdop = undefined 
        radar_out_conv%Zt = undefined ; radar_out_conv%Zm = undefined ; radar_out_conv%Vdop = undefined
        radar_out_2d_conv%Zt_max = undefined ; radar_out_2d_conv%Zm_max = undefined
        radar_out_2d_conv%Het    = undefined ; radar_out_2d_conv%pia    = undefined
     endif


  endif
!
! Radiance output [W/m2/str/micron] VIS and NIR band & tb [K] for IR band (>10 micron)
!
  if(visir) then
     if (visir_sensor == 'AIRS' ) then  !SPECIAL AIRS case

       nw = nw_airs ! total channel #
       if(masterproc) print*,'MSG allocate_all: Special AIRS-MLM case'
       if(masterproc) print*,'total channel number is ',nw
       allocate(radiance_out(mxgridx,mxgridy,nw),&
                radiance_out_fov(mxgridx,mxgridy,nw), stat=ierr )
       if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: visir output allocation error -> Terminate program.')
     
     else !normal

       allocate(radiance_out(mxgridx,mxgridy,mxwavel),&
                radiance_out_fov(mxgridx,mxgridy,mxwavel), stat=ierr ) 
       if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: visir output allocation error -> Terminate program.')

     endif
     radiance_out = undefined ; radiance_out_fov = undefined
  endif   

!
! Lidar total attenuating backscattering coef [m-1 str-1]
!
  if(lidar) then

    if( inst_profile_lidar ) then  !instrumental vertical cordinate

        mxlyr_lidar = NINT( mxhgt_lidar / range_lidar ) !lidar instrumental maximum height level 
        if(masterproc) print*,'MSG allocateall: lidar instrumental max level is =',mxlyr_lidar
        allocate( hgt_stag_lidar(0:mxlyr_lidar) ,&
                  dhgt_lidar    (1:mxlyr_lidar)  &
                )

        ! compute lidar range-based height at interface
        hgt_stag_lidar(0) = 0.  !sea level 0km
        do k_inst = 1, mxlyr_lidar
            hgt_stag_lidar(k_inst)  = REAL(k_inst) * range_lidar  ![km]
        enddo
        dhgt_lidar = range_lidar  ![km]

        allocate(att_B      (mxgridx,mxgridy,mxlyr_lidar,mxwavel_lidar), &
                 att_p22_B  (mxgridx,mxgridy,mxlyr_lidar,mxwavel_lidar), &
                 att_clear_B(mxgridx,mxgridy,mxlyr_lidar,mxwavel_lidar), &
                 lidar_sback(mxgridx,mxgridy,mxlyr_lidar,mxwavel_lidar), &
                 lidar_ratio(mxgridx,mxgridy,mxlyr_lidar,mxwavel_lidar), &
                 optical_depth(mxgridx,mxgridy,mxlyr_lidar,mxwavel_lidar), &
                 stat=ierr )
        if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: lidar output allocation error -> Terminate program.')

   else  !CRM level profile

     allocate(att_B      (mxgridx,mxgridy,mxlyr,mxwavel_lidar), &
              att_p22_B  (mxgridx,mxgridy,mxlyr,mxwavel_lidar), &
              att_clear_B(mxgridx,mxgridy,mxlyr,mxwavel_lidar), &
              lidar_sback(mxgridx,mxgridy,mxlyr,mxwavel_lidar), &
              lidar_ratio(mxgridx,mxgridy,mxlyr,mxwavel_lidar), &
              optical_depth(mxgridx,mxgridy,mxlyr,mxwavel_lidar), &
              stat=ierr )
     if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: lidar output allocation error -> Terminate program. ')

   endif

   !
   ! initialize as undefined
   !
   att_B = undefined ; lidar_sback = undefined ; lidar_ratio = undefined  ; optical_depth = undefined
   att_clear_B = undefined ; att_p22_B = undefined
 

  endif

!
! ISCCP output
!
! if(isccp) then
!    allocate(cld_opt(mxgridx,mxgridy) , & ! column-cloud optical depth [-]
!             cldtop_p(mxgridx,mxgridy), & ! cloud-top pressure [hPa]
!             aerosol_opt(mxgridx,mxgridy) , & ! column-aerosol optical depth [-]
!             aerotop_p(mxgridx,mxgridy), & ! aerosol-top pressure [hPa]
!             stat=ierr)
!    if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: ISCCP output allocation error -> Terminate program.')
!    cld_opt = undefined ; cldtop_p = undefined ; aerosol_opt = undefined ; aerotop_p = undefined
! endif !isccp

!
! Broadband output 
!
 if(broad) then

    allocate(ebudget(mxgridx,mxgridy,14)    , & ! broadband SW/IR energy budget
             ebudget_fov(mxgridx,mxgridy,14), & ! broadband SW/IR energy budget
             aod_broad(mxgridx,mxgridy) , cod_broad(mxgridx,mxgridy), & ! aerosol and cloud column tau [-]
             stat=ierr)
    if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: broadband output allocation error -> Terminate program.')
    ebudget   = undefined ; ebudget_fov = undefined 
    aod_broad = undefined ; cod_broad   = undefined

    if(heating_rate) then
       allocate(sw_heat(mxgridx,mxgridy,mxlyr),lw_heat(mxgridx,mxgridy,mxlyr),stat=ierr)  !heating rate [K/day]
       if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: broadband output allocation error -> Terminate program.') 
       sw_heat = undefined ; lw_heat = undefined 
    endif !heating_rate


 endif !broad



!
! SARTA (AIRS) simulator
!
 if(sarta) then            

    if(all_channel_airs) then

      nchan_sarta = 2378  !default all-chanell AIRS simulation.
      do n = 1, nchan_sarta
         ichan_sarta(n) = n
      enddo

    else 

      nchan_sarta  = nchan_airs      ! number of airs channels
      ichan_sarta(1:mxchan_airs) = ichan_airs(1:mxchan_airs) ! specific airs channel ID numbers
 
    endif

    ! sarta output is domain tile due to size limitation. 

    allocate( sarta_out    (is:ie,js:je,nchan_sarta)    ,&
              sarta_chan_id(            nchan_sarta)    ,&
              sarta_chan_freq(          nchan_sarta)    ,&
              stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_all: sarta output allocation error -> Terminate program.')
    sarta_out       = undefined
    sarta_chan_id   = undefined_i2
    sarta_chan_freq = undefined

 endif


!
! radar echo output
!
 if(polarris) then

   if( .not. use_crm_level_polarris ) then

     mxlyr_polarris = NINT( mxhgt_polarris / range_polarris ) !radar instrumental maximum height level 
     if(masterproc) print*,'MSG allocateall: polarimetric radar instrumental max level is =',mxlyr_polarris
     allocate( hgt_stag_polarris(0:mxlyr_polarris) ,&
               dhgt_polarris    (1:mxlyr_polarris) ,&
                hgt_polarris    (1:mxlyr_polarris)  &
                )

     ! compute radar range-based height at interface
     hgt_stag_polarris(0) = 0.  !sea level 0km
     do k_inst = 1, mxlyr_polarris
        hgt_stag_polarris(k_inst)  = REAL(k_inst) * range_polarris  ![km]
        hgt_polarris(k_inst)  = REAL(k_inst) * range_polarris - 0.5*range_polarris
     enddo
     dhgt_polarris = range_polarris  ![km]

   endif



 !
 endif



!
! 2d and 3D domain array used for dumpting output in MPI domain decomposition 
! array size must be (1:mxgridx,1:mxgridy        ) & (1:mxgridx,1:mxgridy,1:mxlyr)
!
 allocate( out_domain2d(1:mxgridx,1:mxgridy        )    ,& !temporal 2d array for entire domain
           out_domain3d(1:mxgridx,1:mxgridy,1:mxlyr)    ,& !temporal 3d array for entire domain
           stat=ierr )


!
! output latlon array
!
#if MPI == 2
 if( masterproc) then
#endif
   imax_latlon = int( ( maxlon - minlon ) / res_latlon )
   jmax_latlon = int( ( maxlat - minlat ) / res_latlon )

   if( imax_latlon < 1 ) call stop_sdsu(myrank,'MSG allocate_all: maxlon is smaller than minlon')
   if( jmax_latlon < 1 ) call stop_sdsu(myrank,'MSG allocate_all: maxlat is smaller than minlat')

   allocate( latlon_grid(imax_latlon,jmax_latlon), &
             lat_grid(jmax_latlon)          , &
             lon_grid(imax_latlon)          , &
             stat=ierr )

   ! initialized
   do j = 1, jmax_latlon
      lat_grid(j) = minlat + res_latlon * real(j-1) ![deg]
   enddo
   do i = 1, imax_latlon
      lon_grid(i) = minlon + res_latlon * real(i-1) ![deg]
   enddo
   latlon_grid%var = undefined ; latlon_grid%vsum = undefined ;  latlon_grid%cnt = 0

   !derive lat-lon grid area per pixel. 
   do j = 1, jmax_latlon
      call derive_area_latlon(res_latlon, lat_grid(j), area)
      latlon_grid(1:imax_latlon,j)%area = area !area of latlon grid [km2]
   enddo


 

#if MPI == 2
 endif
#endif

!
! now allocated for each processor
!

   allocate(  lat_domain(1:mxgridx,1:mxgridy),&   !domain latitude (mxgridx,mxgridy)  [deg]
              lon_domain(1:mxgridx,1:mxgridy),&   !domain longitude (mxgridx,mxgridy) [deg]
              elev_domain(1:mxgridx,1:mxgridy),&   !domain elevation (mxgridx,mxgridy) [km]
             stat=ierr )



  return
 end subroutine allocate_all

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  derive_area_latlon(res, lat_in, area)
 implicit none
 real(sdsu_fps),intent(in) :: res     !resolution [deg]
 real(sdsu_fps),intent(in) :: lat_in  !resolution [deg]
 real(sdsu_fps),intent(out) :: area   !area of lat lon cordinate [km2] 
 real(sdsu_fps) :: r_earth !earth's radius [km]

 real(sdsu_fps) :: theta_s, theta_e ! angle from south pole [rad]

 theta_s = (90.e0 + lat_in - 0.5*res) * const_degrad !start theta [rad] (south pole is zero radian)
 theta_e = (90.e0 + lat_in + 0.5*res) * const_degrad !start theta [rad] (south pole is zero radian)

 call earth_radius_single( lat_in*const_degrad, r_earth )

!
! derie area per lat-lon pixel [km2]
!
 area = res*const_degrad *r_earth*r_earth * (cos( theta_s ) - cos( theta_e ))  

 
 return
 end subroutine  derive_area_latlon

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine earth_radius_single( lat, re )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Compute earth's radius as spherical oblate (for a given latitude). 
!           Namely, from the earth's core of earth to reference ellipsoid based ono WGS 1984. 
! History:
! 10/2010  Toshi Matsui@NASA GSFC : Initial
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: lat ! latitude [rad]
 real(sdsu_fps),intent(out) :: re ! distance between earth's core and ellipsoid [km]

 real(sdsu_fps) :: a,b,c,d
 real(sdsu_fps),parameter :: req=6378.1370e0     ! earth's equatorial radius [km]
 real(sdsu_fps),parameter :: rpl=6356.7523142e0  ! earth's polar radius [km]

 a = req*req*cos(lat)*req*req*cos(lat)
 b = rpl*rpl*sin(lat)*rpl*rpl*sin(lat)
 c = req*cos(lat)*req*cos(lat)
 d = rpl*sin(lat)*rpl*sin(lat)

 re = sqrt( (a+b)/(c+d) )   ! distance between earth's core and ellipsoid [km]

 return
 end subroutine earth_radius_single

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine LIS_filename_convert( file_name , LIS_file_name )
 implicit none
 character(len=100),intent(in)  :: file_name
 character(len=100),intent(out) :: LIS_file_name

 LIS_file_name = 'lisout_'//file_name(14:16)//'_'//&
 file_name(1:4)//'-'//file_name(5:6)//'-'//file_name(7:8)//'_'//file_name(9:10)//':'//file_name(11:12)//':00'

 return
 end subroutine LIS_filename_convert

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine GEOS_filename_convert( file_name , GEOS_file_name )
 implicit none
 character(len=100),intent(in)  :: file_name
 character(len=100),intent(out) :: GEOS_file_name

 character(len=100) :: file_name_head
 character(len=4) :: yyyy
 character(len=2) :: mm, dd , hh
 integer :: file_len  !total file length
 integer :: i !looping
 integer :: ystart_len  !character length that start four digit number
 integer :: start_len !character start length

!
! file file length
!
 file_len = len(trim(file_name))  !length of character file name

!
! find the position of four digit number from the sdsu_inp_name 
!
 do i = file_len, 1, -1
    if( file_name(i:i) == '+' ) then
        ystart_len = i + 1  ; exit
    endif 
 enddo

 file_name_head = trim(file_name(1:ystart_len-1))
 
!
! get year ~ hour
!
 yyyy=file_name(ystart_len:(ystart_len+3) )
   mm=file_name(ystart_len+4:ystart_len+5)
   dd=file_name(ystart_len+6:ystart_len+7)
   hh=file_name(ystart_len+9:ystart_len+10)

!
! Converted GEOS input file name
!
  GEOS_file_name = trim(file_name_head)//yyyy//'-'//mm//'-'//dd//'_'//hh//':00:00.nc'

 return
 end subroutine GEOS_filename_convert

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_gmt( file_name )
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute GMT from the sdsu_inp_name(characters). GMT will be used in computing local solar
! zenith angle in simulator_broad. Also this routine filter and check the file name compatible
! to SDSU.  
! 
! History:
! 12/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!---------------------------------------------------------------------------------------------------
 character(len=100),intent(in) :: file_name   !Input file names

 integer :: file_len  !total file length
 integer :: i !looping
 integer :: cnt !count
 integer :: ystart_len  !character length that start four digit number
 integer :: start_len !character start length
 integer :: digit  !digit 
 logical :: num   ! true if input is numeric character
 integer :: n    !integer number


!
! initialization
!
  sdsu_yyyy=0. ; sdsu_mm=0. ; sdsu_dd=0. ; sdsu_julian=0. ; sdsu_hh=0. ; sdsu_nn=0. ; sdsu_ss=0. ; sdsu_gmt=0. 
  efile_len = 0

  ystart_len = undefined_i2
!
! file file length
!
 file_len = len(trim(file_name))  !length of character file name

!
! find the position of four digit number from the sdsu_inp_name 
!
 cnt=0
 do i = file_len, 1, -1
    call is_this_num( file_name(i:i), num )
    if(num) then ; cnt = cnt+1 
    else         ; cnt = 0      ; endif 

    if(cnt == 4) then
       ystart_len = i   ; exit
    endif 
       ystart_len = undefined_i2
 enddo !i

!
! action based upon the previous filter.
!
 if( ystart_len == undefined_i2 ) then
    print("(a)"),'MSG get_gmt: File name appers wrong',trim(file_name)
    call stop_sdsu(myrank,'MSG get_gmt: Cannot find year from the file list -> Modify file name as instructed')
 else

!
! Derive year
!
    start_len = ystart_len + 0 ; digit = 3
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_yyyy = sdsu_yyyy + n*(10.e0**(digit-i+start_len)) 
    enddo

    if(sdsu_yyyy < 1977.) then
       print("(a)"),'MSG get_gmt: Operational satellites did not exist in this year ', sdsu_yyyy 
       call stop_sdsu(myrank,'MSG get_gmt: Check input file name -> terminate simulation')
    endif



!
! Derive month
!
    start_len = ystart_len + 5 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_mm = sdsu_mm + n*(10.e0**(digit-i+start_len))
    enddo

    if(sdsu_mm > 12.e0) then
       print*,'MSG get_gmt: Strange month ', sdsu_mm 
       call stop_sdsu(myrank,'MSG get_gmt: Check input file name -> terminate simulation')
    endif

!
! Derive day
!
    start_len = ystart_len + 8 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_dd = sdsu_dd + n*(10.e0**(digit-i+start_len))
    enddo

    if(sdsu_dd > 31.e0) then
       print*,'MSG get_gmt: Strange day ', sdsu_dd 
       call stop_sdsu(myrank,'MSG get_gmt: Check input file name -> terminate simulation')
    endif

!
! Derive hour
!
    start_len = ystart_len + 11 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_hh = sdsu_hh + n*(10.e0**(digit-i+start_len))
    enddo

    if(sdsu_hh > 24.e0) then
       print*,'MSG get_gmt: Strange hour', sdsu_hh 
       call stop_sdsu(myrank,'MSG get_gmt: Check input file name -> terminate simulation')
    endif


!
! special case of MMF GPROF case
! 
if(trim(sim_case) == 'MMF_GPROF') then

 sdsu_nn=0. ; sdsu_ss=0.

 else


!
! Derive minutes
!
    start_len = ystart_len + 14 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_nn = sdsu_nn + n*(10.e0**(digit-i+start_len))
    enddo

    if(sdsu_nn > 60.e0) then
       print*,'MSG get_gmt: Strange minutes ', sdsu_nn
       call stop_sdsu(myrank,'MSG get_gmt: Check input file name -> terminate simulation')
    endif


!
! Derive secounds
!
    start_len = ystart_len + 17 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_ss = sdsu_ss + n*(10.e0**(digit-i+start_len))
    enddo


    if(sdsu_mm > 60.e0) then
       print*,'MSG get_gmt: Strange seconds ', sdsu_ss 
       call stop_sdsu(myrank,'MSG get_gmt: Check input file name -> terminate simulation')
    endif

 
  endif

!
! Derive gmt [hr]
!
  sdsu_gmt = sdsu_hh + sdsu_nn/60.e0 + sdsu_ss/3600.e0

!
! Derive Julian day
!
  call get_julian(sdsu_yyyy, sdsu_mm ,sdsu_dd, sdsu_julian)

!
! Get the end-file length (excluding suffix of sdsu_inp_name)
!
  efile_len = ystart_len + 18


 endif


  return 
 end subroutine get_gmt

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_gmt_geos(file_name)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute GMT from the sdsu_inp_name(characters) for GEOS5 forecast outpout. 
! GMT will be used in computing local solar
! zenith angle in simulator_broad. Also this routine filter and check the file name compatible
! to SDSU.  
! 
! History:
! 12/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!---------------------------------------------------------------------------------------------------
 character(len=100),intent(in) :: file_name   !Input file names

 integer :: file_len  !total file length
 integer :: i !looping
 integer :: cnt !count
 integer :: ystart_len  !character length that start four digit number
 integer :: start_len !character start length
 integer :: digit  !digit 
 logical :: num   ! true if input is numeric character
 integer :: n    !integer number


!
! initialization
!
  sdsu_yyyy=0. ; sdsu_mm=0. ; sdsu_dd=0. ; sdsu_julian=0. ; sdsu_hh=0. ; sdsu_nn=0. ; sdsu_ss=0. ; sdsu_gmt=0. 
  efile_len = 0

  ystart_len = undefined_i2
!
! file file length
!
 file_len = len(trim(file_name))  !length of character file name

!
! find the position of four digit number from the sdsu_inp_name 
! in case of GEOS5 forecast data it is right after "+" sign. 
!

 cnt=0
 do i = file_len, 1, -1
    if( file_name(i:i) == '+' ) then
        ystart_len = i+1 ; exit
    endif
 enddo

!
! action based upon the previous filter.
!
 if( ystart_len == undefined_i2 ) then
    call stop_sdsu(myrank,'Cannot find year from the file list -> Modify file name as instructed')
 else

!
! Derive year
!
    start_len = ystart_len + 0 ; digit = 3
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_yyyy = sdsu_yyyy + n*(10.e0**(digit-i+start_len)) 
    enddo

    if(sdsu_yyyy < 1977.) then
       print*,'Operational satellites did not exist in this year ', sdsu_yyyy 
       call stop_sdsu(myrank,'Check input file name -> terminate simulation')
    endif

!
! Derive month
!
    start_len = ystart_len + 4 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_mm = sdsu_mm + n*(10.e0**(digit-i+start_len))
    enddo

    if(sdsu_mm > 12.e0) then
       print*,'Strange month ', sdsu_mm 
       call stop_sdsu(myrank,'Check input file name -> terminate simulation')
    endif

!
! Derive day
!
    start_len = ystart_len + 6 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_dd = sdsu_dd + n*(10.e0**(digit-i+start_len))
    enddo

    if(sdsu_dd > 31.e0) then
       print*,'Strange day ', sdsu_dd 
       call stop_sdsu(myrank,'Check input file name -> terminate simulation')
    endif

!
! Derive hour
!
    start_len = ystart_len + 9 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_hh = sdsu_hh + n*(10.e0**(digit-i+start_len))
    enddo

    if(sdsu_hh > 24.e0) then
       print*,'Strange hour', sdsu_hh 
       call stop_sdsu(myrank,'Check input file name -> terminate simulation')
    endif

!
! Derive minutes
!
    start_len = ystart_len + 11 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       sdsu_nn = sdsu_nn + n*(10.e0**(digit-i+start_len))
    enddo

    if(sdsu_nn > 60.e0) then
       print*,'Strange minutes ', sdsu_nn
       call stop_sdsu(myrank,'Check input file name -> terminate simulation')
    endif
!
! Derive secounds
!
    sdsu_ss = 0.0

!
! Derive gmt [hr]
!
  sdsu_gmt = sdsu_hh + sdsu_nn/60.e0 + sdsu_ss/3600.e0

!
! Derive Julian day
!
  call get_julian(sdsu_yyyy, sdsu_mm ,sdsu_dd, sdsu_julian)

!
! Get the end-file length (excluding suffix of sdsu_inp_name)
!
 do i = file_len, 1, -1
    if( file_name(i:i) == '.' ) then
        efile_len = i-1 ; exit
    endif
 enddo

 endif

 return
 end subroutine get_gmt_geos

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine is_this_num( char_num, num )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  This subroutine check input character is numerical or not.  
! 
! History:
! 12/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 character(len=1),intent(in) :: char_num !numeric number input as character 
 logical,intent(out) :: num !logical output if input is numeric number 

 num = .false.
 if( char_num == '0') num=.true.
 if( char_num == '1') num=.true.
 if( char_num == '2') num=.true.
 if( char_num == '3') num=.true.
 if( char_num == '4') num=.true.
 if( char_num == '5') num=.true.
 if( char_num == '6') num=.true.
 if( char_num == '7') num=.true.
 if( char_num == '8') num=.true.
 if( char_num == '9') num=.true.
 return
 end subroutine is_this_num

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
 subroutine char_to_int( char_num, num )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  This subroutine convert character number to integer number. . 
! 
! History:
! 12/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 character(len=1),intent(in) :: char_num !character number
 integer,intent(out) :: num              !integer number

 num = undefined_i2
 if( char_num == '0') num=0
 if( char_num == '1') num=1
 if( char_num == '2') num=2
 if( char_num == '3') num=3
 if( char_num == '4') num=4
 if( char_num == '5') num=5
 if( char_num == '6') num=6
 if( char_num == '7') num=7
 if( char_num == '8') num=8
 if( char_num == '9') num=9

 if(num == undefined_i2) then
  print*,char_num,'MSG char_to_int; wrong input of char_num'
  call stop_sdsu(myrank,'MSG char_to_int; wrong input of char_num')
 endif

 end subroutine char_to_int

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine get_julian(yyyy,mm,dd,julian)
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
! Compute julian date from Month, Day, and Year.  
! 
! History:
! 12/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: yyyy,mm,dd !year,month,day
 real(sdsu_fps),intent(out) :: julian    !julian day
 real(sdsu_fps),parameter :: day(12) = (/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./) !for normal year
 real(sdsu_fps),parameter :: dayl(12)= (/31.,29.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./) !for leap year  

  if( mod(yyyy,4.) == 0. ) then !leap year
     julian=sum( dayl(1:(int(mm)-1) ) ) + dd
  else                           ! normal year
     julian=sum( day(1:(int(mm)-1) ) ) + dd
  endif

 end subroutine get_julian


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_grads_time
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!   derive time character for grads control file. 
! 
! History:
! 07/2011  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 character :: &
    hh*2 ,&
    mm*2 ,&
    dd*2 ,&
    mon*3,&
    yyyy*4


 if(masterproc) print*,'MSG get_grads_time: derive grads time character'

 write(hh,"(I2.2)") INT(sdsu_hh)
 write(mm,"(I2.2)") INT(sdsu_nn)
 write(dd,"(I2.2)") INT(sdsu_dd)

 if( sdsu_mm == 1. )  mon = 'jan'
 if( sdsu_mm == 2. )  mon = 'feb' 
 if( sdsu_mm == 3. )  mon = 'mar' 
 if( sdsu_mm == 4. )  mon = 'apr' 
 if( sdsu_mm == 5. )  mon = 'may' 
 if( sdsu_mm == 6. )  mon = 'jun' 
 if( sdsu_mm == 7. )  mon = 'jul' 
 if( sdsu_mm == 8. )  mon = 'aug' 
 if( sdsu_mm == 9. )  mon = 'sep' 
 if( sdsu_mm == 10.)  mon = 'oct' 
 if( sdsu_mm == 11.)  mon = 'nov' 
 if( sdsu_mm == 12.)  mon = 'dec' 

 write(yyyy,"(I4.4)") INT(sdsu_yyyy)

!
! get grads time character
!
 grads_time = hh//':'//mm//'Z'//dd//mon//yyyy

 return
 end subroutine get_grads_time

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine bulk_DSD
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
! Set constants for gamma size distributions
! for different microphysics scheme.
! 
! History:
! 06/2011  Toshi Matsui@NASA GSFC : Added Morrison 2-moment
! 04/2008  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps) :: re !sude re (not used in this routine)


    if(trim(cloud_microphysics) == 'GOD') then
       if(masterproc) print*,'MSG bulk_DSD ;  Exponential DSD prescribed parameters for Goddard Microphysics'

       n0_gen%rain = 0.08e+8   ! intercept for rain [1/m4] (=0.08[1/cm4])
       n0_gen%snow = 0.16e+8   ! intercept for snow [1/m4] (=0.16[1/cm4])
       n0_gen%graupel= 0.04e+8 ! intercept for graupel [1/m4] (=0.04[1/cm4])
       n0_gen%hail = 0.002e+8  ! intercept for hail [1/m4] (=0.002[1/cm4]) 

       rho_gen%rain = 1.0e+3  ! density of rain [kg/m3] (=1.0[g/cm3])
       rho_gen%snow = 0.1e+3  ! density of snow [kg/m3] (=0.1[g/cm3]) 
       rho_gen%graupel= 0.4e+3  ! density of graupel [kg/m3] (=0.4[g/cm3])
       rho_gen%hail = 0.9e+3  ! density of hail [kg/m3] (=0.9[g/cm3])

       call re_LUT_Heymsfield_Platt_1984('init',270.,1., re)

       mu_gen%cloud = 0.   !
       mu_gen%rain  = 0.   !
       mu_gen%ice   = 0.   !
       mu_gen%snow  = 0.   !
       mu_gen%graupel = 0. !
       mu_gen%hail  = 0.   !

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%ice     = 0.0
       b_vt%ice     = 0.0
       a_vt%snow    = 1.30493
       b_vt%snow    = 0.11
       a_vt%graupel = 19.3
       b_vt%graupel = 0.37
       a_vt%rain    = 841.99667
       b_vt%rain    = 0.8
       a_vt%hail    = 841.99667
       b_vt%hail    = 0.8
       vt_max%cloud = 0.  ![m/s]
       vt_max%rain = 9.1  ![m/s]
       vt_max%ice = 0.    ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 9.1  ![m/s]


    elseif(trim(cloud_microphysics) == 'GOD10') then
       if(masterproc) print("(a)"),'MSG bulk_DSD ;  Exponential DSD prescribed parameters for Goddard Microphysics 2010'

       n0_gen%rain = 0.08e+8   ! intercept for rain [1/m4] (=0.08[1/cm4])
       n0_gen%snow = 0.16e+8   ! intercept for snow [1/m4] (=0.16[1/cm4]) !this is default (actually T-q function)
       n0_gen%graupel= 0.04e+8   ! intercept for graupel [1/m4] (=0.04[1/cm4]) !his is default (actually T-q function)
       n0_gen%hail = 0.002e+8  ! intercept for hail [1/m4] (=0.002[1/cm4]) 

       rho_gen%rain = 1.0e+3  ! density of rain [kg/m3] (=1.0[g/cm3])
       rho_gen%snow = 0.05e+3  ! density of snow [kg/m3] (=0.05[g/cm3]) (Note it's reduced from GOD.)
       rho_gen%graupel= 0.3e+3  ! density of graupel [kg/m3] (=0.2[g/cm3]) (Note it's reduced from GOD.)
       rho_gen%hail = 0.9e+3  ! density of hail [kg/m3] (=0.9[g/cm3])

       call re_LUT_Heymsfield_Platt_1984('init',270.,1., re)

       mu_gen%cloud = 0.   !
       mu_gen%rain  = 0.   !
       mu_gen%ice   = 0.   !
       mu_gen%snow  = 0.   !
       mu_gen%graupel = 0. !
       mu_gen%hail  = 0.   !

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%ice     = 0.0
       b_vt%ice     = 0.0
       a_vt%snow    = 1.30493
       b_vt%snow    = 0.11
       a_vt%graupel = 19.3
       b_vt%graupel = 0.37
       a_vt%rain    = 841.99667
       b_vt%rain    = 0.8
       a_vt%hail    = 841.99667
       b_vt%hail    = 0.8

       vt_max%cloud = 0.  ![m/s]
       vt_max%rain = 9.1  ![m/s]
       vt_max%ice = 0.    ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 9.1  ![m/s]

    elseif(trim(cloud_microphysics) == 'GMP4ICE') then
       if(masterproc) print("(a)"),'MSG bulk_DSD ;  Exponential DSD prescribed parameters for Goddard Microphysics 4ICE'

       n0_gen%rain = 0.08e+8   ! intercept for rain [1/m4] (=0.08[1/cm4])
       n0_gen%snow = 0.10e+8   ! intercept for snow [1/m4] (=0.10[1/cm4]) !this is default (actually T-q function)
       n0_gen%graupel= 0.04e+8   ! intercept for graupel [1/m4] (=0.04[1/cm4]) !his is default (actually T-q function)
       n0_gen%hail = 0.002e+8  ! intercept for hail [1/m4] (=0.002[1/cm4]) 

       rho_gen%rain = 1.0e+3  ! density of rain [kg/m3] (=1.0[g/cm3])
       rho_gen%snow = 0.05e+3  ! density of snow [kg/m3] (=0.05[g/cm3]) (its T relationship )
       rho_gen%graupel= 0.3e+3  ! density of graupel [kg/m3] (=0.2[g/cm3]) 
       rho_gen%hail = 0.9e+3  ! density of hail [kg/m3] (=0.9[g/cm3])

       call re_LUT_Heymsfield_Platt_1984('init',270.,1., re)

       mu_gen%cloud = 0.   !
       mu_gen%rain  = 0.   !
       mu_gen%ice   = 0.   !
       mu_gen%snow  = 0.   !
       mu_gen%graupel = 0. !
       mu_gen%hail  = 0.   !

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%ice     = 0.0       
       b_vt%ice     = 0.0
       a_vt%snow    = 1.30493
       b_vt%snow    = 0.11
       a_vt%graupel = 19.3
       b_vt%graupel = 0.37
       a_vt%rain    = 841.99667
       b_vt%rain    = 0.8
       a_vt%hail    = 841.99667
       b_vt%hail    = 0.8

       vt_max%cloud = 0.  ![m/s]
       vt_max%rain = 9.1  ![m/s]
       vt_max%ice = 0.    ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 9.1  ![m/s]


    elseif(trim(cloud_microphysics) == 'LIN') then
       if(masterproc) print*,'MSG bulk_DSD ;  Exponential DSD prescribed parameters for LIN scheme'
       n0_gen%rain = 0.08e+8   ! intercept for rain [1/m4] (=0.08[1/cm4])
       n0_gen%snow = 0.03e+8   ! intercept for snow [1/m4] (=0.03[1/cm4])
       n0_gen%graupel= 0.04e+8   ! intercept for graupel [1/m4] (=0.04[1/cm4])
       n0_gen%hail = 0.002e+8  ! does not exist

       rho_gen%rain = 1.0e+3  ! density of rain [kg/m3] (=1.0[g/cm3])
       rho_gen%snow = 0.1e+3  ! density of snow [kg/m3] (=0.1[g/cm3]) 
       rho_gen%graupel= 0.4e+3  ! density of graupel [kg/m3] (=0.4[g/cm3])
       rho_gen%hail = 0.9e+3  ! does not exist

       call re_LUT_Heymsfield_Platt_1984('init',270.,1., re)
       mu_gen%cloud = 0.   !
       mu_gen%rain  = 0.   !
       mu_gen%ice   = 0.   !
       mu_gen%snow  = 0.   !
       mu_gen%graupel = 0. !
       mu_gen%hail  = 0.   !

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%ice     = 0.0
       b_vt%ice     = 0.0
       a_vt%snow    = 11.72
       b_vt%snow    = 0.41
       a_vt%graupel = 19.3
       b_vt%graupel = 0.37
       a_vt%rain    = 841.99667
       b_vt%rain    = 0.8
       a_vt%hail    = 841.99667
       b_vt%hail    = 0.8

       vt_max%cloud = 0.  ![m/s]
       vt_max%rain = 9.1  ![m/s]
       vt_max%ice = 0.    ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 9.1  ![m/s]

    elseif(trim(cloud_microphysics) == 'WSM' .or. trim(cloud_microphysics) == 'WDM' ) then
       if(masterproc) print*,'MSG bulk_DSD ; Exponential DSD prescribed parameters of WSM/WDM'

!toshiwdm: modify here
       n0_gen%rain = 0.08e+8   ! intercept for rain [1/m4] (=0.08[1/cm4])
       n0_gen%snow = 0.02e+8   ! -> temperature dependent (fake)
       n0_gen%graupel= 0.04e+8   ! intercept for graupel [1/m4] (=0.04[1/cm4])
       n0_gen%hail = 0.002e+8  ! does not exist  

       rho_gen%cloud  = 0.997e+3  ! density of cloud   [kg/m3] (=1.0[g/cm3])
       rho_gen%rain = 1.0e+3  ! density of rain [kg/m3] (=1.0[g/cm3])
       rho_gen%snow = 0.1e+3  ! density of snow [kg/m3] (=0.1[g/cm3]) 
       rho_gen%graupel= 0.5e+3  ! density of graupel [kg/m3] (=0.4[g/cm3])
       rho_gen%hail = 0.9e+3  ! does not exist

       call re_LUT_Heymsfield_Platt_1984('init',270.,1., re)

       mu_gen%cloud = 0.   !
       mu_gen%rain  = 0.   !
       mu_gen%ice   = 0.   !
       mu_gen%snow  = 0.   !
       mu_gen%graupel = 0. !
       mu_gen%hail  = 0.   !

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%ice     = 0.0
       b_vt%ice     = 0.0
       a_vt%snow    = 11.72
       b_vt%snow    = 0.41
       a_vt%graupel = 19.3
       b_vt%graupel = 0.37
       a_vt%rain    = 841.99667
       b_vt%rain    = 0.8
       a_vt%hail    = 841.99667
       b_vt%hail    = 0.8

       vt_max%cloud = 0.  ![m/s]
       vt_max%rain = 9.1  ![m/s]
       vt_max%ice = 0.    ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 9.1  ![m/s]


    elseif(trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH') then

! modify --> Mei Han
       mu_gen%cloud = 2.   !
       mu_gen%rain  = 0.   !
       mu_gen%ice   = 0.   !
       mu_gen%snow  = 0.   !
       mu_gen%graupel = 0. !
       mu_gen%hail  = 0.   !

       rho_gen%cloud  = 0.997e+3  ! density of cloud   [kg/m3] (=1.0[g/cm3])
       rho_gen%rain   = 0.997e+3  ! density of rain    [kg/m3] (=1.0[g/cm3])
       rho_gen%ice    = 0.1e+3  ! density of ice     [kg/m3] (=0.1[g/cm3]) 
       rho_gen%snow   = 0.1e+3  ! density of snow    [kg/m3] (=0.1[g/cm3]) 
       rho_gen%graupel= 0.4e+3  ! density of graupel [kg/m3] (=0.4[g/cm3])
       rho_gen%hail   = 0.9e+3  ! density of hail (does not exit, but spefify for)

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%rain    = 841.99667
       b_vt%rain    = 0.8
       a_vt%ice     = 0.0
       b_vt%ice     = 0.0
       a_vt%snow    = 11.72
       b_vt%snow    = 0.41
       a_vt%graupel = 19.3
       b_vt%graupel = 0.37
       a_vt%hail    = 114.5
       b_vt%hail    = 0.5

       vt_max%cloud = 0.  ![m/s]
       vt_max%rain = 9.1  ![m/s]
       vt_max%ice = 0.    ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 9.1  ![m/s]

    elseif(trim(cloud_microphysics) == 'P3' ) then

       !these are fake for P3 scheme
       mu_gen%cloud = 2.   !
       mu_gen%rain  = 0.   !
       mu_gen%ice   = 0.   !
       mu_gen%snow  = 0.   !
       mu_gen%graupel = 0. !
       mu_gen%hail  = 0.   !

       rho_gen%cloud  = 0.997e+3  ! density of cloud   [kg/m3] (=1.0[g/cm3])
       rho_gen%rain   = 0.997e+3  ! density of rain    [kg/m3] (=1.0[g/cm3])
       rho_gen%ice    = 0.1e+3  ! density of ice     [kg/m3] (=0.1[g/cm3]) 
       rho_gen%snow   = 0.1e+3  ! density of snow    [kg/m3] (=0.1[g/cm3]) 
       rho_gen%graupel= 0.4e+3  ! density of graupel [kg/m3] (=0.4[g/cm3])
       rho_gen%hail   = 0.9e+3  ! density of hail (does not exit, but spefify for)

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%rain    = 841.99667
       b_vt%rain    = 0.8
       a_vt%ice     = 0.0
       b_vt%ice     = 0.0
       a_vt%snow    = 11.72
       b_vt%snow    = 0.41
       a_vt%graupel = 19.3
       b_vt%graupel = 0.37
       a_vt%hail    = 114.5
       b_vt%hail    = 0.5

       vt_max%cloud = 0.  ![m/s]
       vt_max%rain = 9.1  ![m/s]
       vt_max%ice = 0.    ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 9.1  ![m/s]


    elseif(trim(cloud_microphysics) == 'NSSL2MMT') then

       !note that mu values are diffrent than an ordinal gamma DSD in NSSL2MMT scheme
       mu_gen%cloud = 0.   ! cnu  N(v)
       mu_gen%ice   = 0.   ! cinu N(v)
       mu_gen%snow  = -0.8 ! snu  N(v)
       mu_gen%rain  = 0.   !  alphar N(d)
       mu_gen%graupel = 0. !  alphah N(d)
       mu_gen%hail  = 2.   !  alphahl N(d)

       rho_gen%cloud  = 1.e+3     ! density of cloud   [kg/m3] (=1.0[g/cm3])
       rho_gen%rain   = 1.e+3     ! density of rain    [kg/m3] (=1.0[g/cm3])
       rho_gen%ice    = 0.917e+3  ! density of ice     [kg/m3] (=0.1[g/cm3]) 
       rho_gen%snow   = 0.1e+3    ! density of snow    [kg/m3] (=0.1[g/cm3]) 
       !density of graupel and hail are predicted  in NSSL2MMT scheme

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%rain    = 841.9    !this is fake
       b_vt%rain    = 0.8      !this is fake
       a_vt%ice     = 42.2999  !modified Ferrier (1994)
       b_vt%ice     = 0.55     !modified Ferrier (1994)
       a_vt%snow    = 8.96613 !Ferrier (1994)
       b_vt%snow    = 0.42    !Ferrier (1994)
       a_vt%graupel = 19.3    !this is fake
       b_vt%graupel = 0.37    !this is fake
       a_vt%hail    = 114.5   !this is fake
       b_vt%hail    = 0.5     !this is fake

       vt_max%rain = 9.1  ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 20.  ![m/s]


    elseif(trim(cloud_microphysics) == 'THOM') then

! modify --> Mei Han

!Cloud water: 
!               m-D relationship and density are the same as rain, but with a non-zero mu.
!               It assumes the number concentration (not the intercept) for maritime case as Nt_c as below:
!               Nt_c = 100.E6  ![/m^3]
!               mu_c = MIN(15., (1000.E6/Nt_c + 2.))
       mu_gen%cloud = min(15., 1000.e6/Nt_c_thom + 2. )   !
       mu_gen%rain  = 0.   !
       mu_gen%ice   = 0.   !
       mu_gen%snow  = 0.   !
       mu_gen%graupel = 0. !
       mu_gen%hail  = 0.   !

       rho_gen%cloud  = 0.997e+3  ! density of cloud   [kg/m3] (=1.0[g/cm3])
       rho_gen%rain   = 0.997e+3  ! density of rain    [kg/m3] (=1.0[g/cm3])
       rho_gen%ice    = 0.89e+3  ! density of ice     [kg/m3] (=0.1[g/cm3]) 
       rho_gen%snow   = 0.1e+3  ! density of snow    [kg/m3] (=0.1[g/cm3]) 
       rho_gen%graupel= 0.4e+3  ! density of graupel [kg/m3] (=0.4[g/cm3])
       rho_gen%hail   = 0.9e+3  ! density of hail (does not exit, but spefify for)

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%ice     = 0.0
       b_vt%ice     = 0.0
       a_vt%snow    = 11.72
       b_vt%snow    = 0.41
       a_vt%graupel = 19.3
       b_vt%graupel = 0.37
       a_vt%rain    = 841.99667
       b_vt%rain    = 0.8
       a_vt%hail    = 841.99667
       b_vt%hail    = 0.8

       vt_max%cloud = 0.  ![m/s]
       vt_max%rain = 9.1  ![m/s]
       vt_max%ice = 0.    ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 9.1  ![m/s]

    elseif(trim(cloud_microphysics) == 'NTUCLR') then

! modify --> 

       n0_gen%snow    = 2.e+7   ! intercept for snow [1/m4] (=0.16[1/cm4]) 
       n0_gen%graupel = 4.e+6   ! intercept for graupel [1/m4] (=0.04[1/cm4])


       mu_gen%cloud = 0.   !
       mu_gen%rain  = 0.   !
       mu_gen%ice   = 0.   !
       mu_gen%snow  = 0.   !
       mu_gen%graupel = 0. !
       mu_gen%hail  = 0.   !

       rho_gen%cloud  = 0.9998396e+3  ! density of cloud   [kg/m3] (=1.0[g/cm3])
       rho_gen%rain   = 0.9998396e+3  ! density of rain    [kg/m3] (=1.0[g/cm3])
       rho_gen%ice    = 0.5e+3  ! density of ice     [kg/m3] (=0.1[g/cm3]) 
       rho_gen%snow   = 0.1e+3  ! density of snow    [kg/m3] (=0.1[g/cm3]) 
       rho_gen%graupel= 0.4e+3  ! density of graupel [kg/m3] (=0.4[g/cm3])
       rho_gen%hail   = 0.9e+3  ! density of hail (does not exit, but spefify for)

       !parameters for Vt
       a_vt%cloud   = 0.0
       b_vt%cloud   = 0.0
       a_vt%ice     = 700.0
       b_vt%ice     = 1.0
       a_vt%snow    = 11.72
       b_vt%snow    = 0.41
       a_vt%graupel = 19.3
       b_vt%graupel = 0.37
       a_vt%rain    = 841.99667
       b_vt%rain    = 0.8
       a_vt%hail    = 841.99667
       b_vt%hail    = 0.8

       vt_max%cloud = 0.  ![m/s]
       vt_max%rain = 9.1  ![m/s]
       vt_max%ice = 0.    ![m/s]
       vt_max%snow = 1.2  ![m/s]
       vt_max%graupel = 20. ![m/s]
       vt_max%hail = 9.1  ![m/s]


   endif 

 if(masterproc) print*,''

 end subroutine bulk_DSD

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_SBM_bin
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
! Read SBM particle size bin tables in ascii format.
! 
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
! Li, X., W.-K. Tao, A. Khain, J. Simpson and D. Johnson, 2008: Sensitivity of a cloud-resolving 
! model to bulk and explicit-bin microphysics schemes: Part I: Comparisons.  J. Atmos. Sci., (accepted).
!
! Li, X., W.-K. Tao, A. Khain, J. Simpson and D. Johnson, 2008: Sensitivity of a cloud-resolving 
! model to bulk and explicit-bin microphysics schemes:: Part II: Cloud microphysics and storm 
! dynamics interactions.  J. Atmos. Sci., (accepted).
!-----------------------------------------------------------------------------------------------------
 integer :: i,j ! loop indices
 integer :: ierr !error index

!test
 real :: fmelt_bin, frime_bin
 real :: den_diag, rad_diag, rad, rad_liq, fmelt

 if(masterproc) print*,'MSG read_SBM_bin: Reading bulk masses and density of SBM scheme'

!
! set up mass-bin size
!
 if( trim(cloud_microphysics)  == 'HUCM_SBM'   ) nbin = 33  !33 mass bin 
 if( trim(cloud_microphysics)  == 'HUCM_SBM43' ) nbin = 43  !43 mass bin

!
! allocate mass bin parameter
!
 allocate( &
          x_sbm   (1:nbin) ,& ! mass per particle [g]
          den_sbm (1:nbin) ,& ! density per particle [g/cm3]
          rad_sbm (1:nbin) ,& ! radius of particle [cm]
          drad_sbm(1:nbin) ,& ! d (radius) [cm]
          vt_sbm  (1:nbin) ,& ! terminal velocity [cm/s]
          brad_sbm(0:nbin) ,& ! boundary of radius bin [cm]
          stat=ierr )


!
! open/read bulk masses of SBM scheme
!
 if( trim(cloud_microphysics)  == 'HUCM_SBM'   ) sdsu_io_file = trim(sdsu_dir_data)//'masses_sbm.asc'
 if( trim(cloud_microphysics)  == 'HUCM_SBM43' ) sdsu_io_file = trim(sdsu_dir_data)//'masses_sbm43.asc'

! open(1,file=trim(sdsu_io_file),status='old')
 call open_oldfile(1, trim(sdsu_io_file) )
 read(1,900) x_sbm%liq, x_sbm%ice_col, x_sbm%ice_pla, x_sbm%ice_den, x_sbm%snow, x_sbm%graupel, x_sbm%hail
 close(1)

!
! open/read bulk density of SBM scheme (HUCM SBM assume equivalent-volume sphere for all species)
! , which is different from maximum diameter in field measurements. 
!
 if( trim(cloud_microphysics)  == 'HUCM_SBM'   ) sdsu_io_file = trim(sdsu_dir_data)//'bulkdens_sbm.asc'
 if( trim(cloud_microphysics)  == 'HUCM_SBM43' ) sdsu_io_file = trim(sdsu_dir_data)//'bulkdens_sbm43.asc'

! open(1,file=trim(sdsu_io_file),status='old')
 call open_oldfile(1, trim(sdsu_io_file) )
 read(1,900) den_sbm%liq, den_sbm%ice_col, den_sbm%ice_pla, den_sbm%ice_den, den_sbm%snow, den_sbm%graupel, den_sbm%hail
 close(1)


!
! open/read terminal velocity of SBM scheme
!
 if( trim(cloud_microphysics)  == 'HUCM_SBM'   ) sdsu_io_file = trim(sdsu_dir_data)//'termvels_sbm.asc'
 if( trim(cloud_microphysics)  == 'HUCM_SBM43' ) sdsu_io_file = trim(sdsu_dir_data)//'termvels_sbm43.asc'

! open(1,file=trim(sdsu_io_file),status='old')
 call open_oldfile(1, trim(sdsu_io_file) )
 read(1,900) vt_sbm%liq, vt_sbm%ice_col, vt_sbm%ice_pla, vt_sbm%ice_den, vt_sbm%snow, vt_sbm%graupel, vt_sbm%hail
 close(1)


 900 format(6e13.5)

!
! compute bulk radius [cm]
!
 do i = 1, nbin
    rad_sbm(i)%liq     = ((3.e0*x_sbm(i)%liq    /4.e0/const_pi/den_sbm(i)%liq    )**(1.e0/3.e0))
    rad_sbm(i)%ice_col = ((3.e0*x_sbm(i)%ice_col/4.e0/const_pi/den_sbm(i)%ice_col)**(1.e0/3.e0))
    rad_sbm(i)%ice_pla = ((3.e0*x_sbm(i)%ice_pla/4.e0/const_pi/den_sbm(i)%ice_pla)**(1.e0/3.e0))
    rad_sbm(i)%ice_den = ((3.e0*x_sbm(i)%ice_den/4.e0/const_pi/den_sbm(i)%ice_den)**(1.e0/3.e0))
    rad_sbm(i)%snow    = ((3.e0*x_sbm(i)%snow   /4.e0/const_pi/den_sbm(i)%snow   )**(1.e0/3.e0))
    rad_sbm(i)%graupel = ((3.e0*x_sbm(i)%graupel/4.e0/const_pi/den_sbm(i)%graupel)**(1.e0/3.e0))
    rad_sbm(i)%hail    = ((3.e0*x_sbm(i)%hail   /4.e0/const_pi/den_sbm(i)%hail   )**(1.e0/3.e0))
 enddo


!
! compute boundary of size bin [cm]
!

 do i = 0, nbin
      if(i==0) then
         brad_sbm(i)%liq     = rad_sbm(1)%liq     - ( rad_sbm(2)%liq     - rad_sbm(1)%liq     )/2.e0 
         brad_sbm(i)%ice_col = rad_sbm(1)%ice_col - ( rad_sbm(2)%ice_col - rad_sbm(1)%ice_col )/2.e0
         brad_sbm(i)%ice_pla = rad_sbm(1)%ice_pla - ( rad_sbm(2)%ice_pla - rad_sbm(1)%ice_pla )/2.e0
         brad_sbm(i)%ice_den = rad_sbm(1)%ice_den - ( rad_sbm(2)%ice_den - rad_sbm(1)%ice_den )/2.e0
         brad_sbm(i)%snow    = rad_sbm(1)%snow    - ( rad_sbm(2)%snow    - rad_sbm(1)%snow    )/2.e0
         brad_sbm(i)%graupel = rad_sbm(1)%graupel - ( rad_sbm(2)%graupel - rad_sbm(1)%graupel )/2.e0
         brad_sbm(i)%hail    = rad_sbm(1)%hail    - ( rad_sbm(2)%hail    - rad_sbm(1)%hail    )/2.e0
     elseif(i<nbin) then
         brad_sbm(i)%liq     = ( rad_sbm(i)%liq     + rad_sbm(i+1)%liq     )/2.e0 
         brad_sbm(i)%ice_col = ( rad_sbm(i)%ice_col + rad_sbm(i+1)%ice_col )/2.e0
         brad_sbm(i)%ice_pla = ( rad_sbm(i)%ice_pla + rad_sbm(i+1)%ice_pla )/2.e0
         brad_sbm(i)%ice_den = ( rad_sbm(i)%ice_den + rad_sbm(i+1)%ice_den )/2.e0
         brad_sbm(i)%snow    = ( rad_sbm(i)%snow    + rad_sbm(i+1)%snow    )/2.e0
         brad_sbm(i)%graupel = ( rad_sbm(i)%graupel + rad_sbm(i+1)%graupel )/2.e0
         brad_sbm(i)%hail    = ( rad_sbm(i)%hail    + rad_sbm(i+1)%hail    )/2.e0
      elseif(i==nbin) then
         brad_sbm(i)%liq     = rad_sbm(nbin)%liq     + ( rad_sbm(nbin)%liq     - rad_sbm(nbin-1)%liq     )/2.e0 
         brad_sbm(i)%ice_col = rad_sbm(nbin)%ice_col + ( rad_sbm(nbin)%ice_col - rad_sbm(nbin-1)%ice_col )/2.e0
         brad_sbm(i)%ice_pla = rad_sbm(nbin)%ice_pla + ( rad_sbm(nbin)%ice_pla - rad_sbm(nbin-1)%ice_pla )/2.e0
         brad_sbm(i)%ice_den = rad_sbm(nbin)%ice_den + ( rad_sbm(nbin)%ice_den - rad_sbm(nbin-1)%ice_den )/2.e0
         brad_sbm(i)%snow    = rad_sbm(nbin)%snow    + ( rad_sbm(nbin)%snow    - rad_sbm(nbin-1)%snow    )/2.e0
         brad_sbm(i)%graupel = rad_sbm(nbin)%graupel + ( rad_sbm(nbin)%graupel - rad_sbm(nbin-1)%graupel )/2.e0
         brad_sbm(i)%hail    = rad_sbm(nbin)%hail    + ( rad_sbm(nbin)%hail    - rad_sbm(nbin-1)%hail    )/2.e0
      endif
 enddo !nbin

!
! compute width of seize bins [cm]
!
 do i = 1, nbin
         drad_sbm(i)%liq     =  brad_sbm(i)%liq     - brad_sbm(i-1)%liq     
         drad_sbm(i)%ice_col =  brad_sbm(i)%ice_col - brad_sbm(i-1)%ice_col 
         drad_sbm(i)%ice_pla =  brad_sbm(i)%ice_pla - brad_sbm(i-1)%ice_pla    
         drad_sbm(i)%ice_den =  brad_sbm(i)%ice_den - brad_sbm(i-1)%ice_den    
         drad_sbm(i)%snow    =  brad_sbm(i)%snow    - brad_sbm(i-1)%snow    
         drad_sbm(i)%graupel =  brad_sbm(i)%graupel - brad_sbm(i-1)%graupel 
         drad_sbm(i)%hail    =  brad_sbm(i)%hail    - brad_sbm(i-1)%hail    
 enddo !nbin

 return
 end subroutine read_SBM_bin

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine diag_rime_snow (ibin, fr, dens, xs, rs, &
                           dens_rime, rs_rime )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Diagnose bulk density and radius for rimed aggregate particle. 
! 
! History:
!   06/2011  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)        :: ibin  ! bin #
 real(sdsu_fps),intent(in) :: fr    ! riming fraction [-]
 real(sdsu_fps),intent(in) :: dens(1:nbin) ! snow bulk density from LUT (without riming) [g/cm3] 
 real(sdsu_fps),intent(in) :: xs(1:nbin)   ! snow bulk mass per particlefrom LUT (without riming) [g]
 real(sdsu_fps),intent(in) :: rs(1:nbin)   ! snow bulk radius per particle from LUT (without riming) [cm]

 real(sdsu_fps),intent(out) :: dens_rime ! diagnostic bulk (effective) density of snow witm riming.  [g/cm3]
 real(sdsu_fps),intent(out) :: rs_rime   ! bulk radius of snow taking with riming [cm]

! local params
 real(sdsu_fps) :: xs_dry   !snow bulk mass without riming [g]
 real(sdsu_fps) :: dens_dry !snow bulk density without riming [g/cm3]
 real(sdsu_fps),parameter :: rho_rime = 0.9e0  !rimed density  [g/cm3]

 integer :: i !looping index

!
! Take account for riming effect upon bulk snow aggregate density and radius. 
!
 if     (fr == 0.e0) then

   dens_rime = dens(ibin)     ! rime fraction 0--> dry snow agg

 elseif (fr == 1.e0) then

   dens_rime = rho_rime          ! rime fraction 1--> 0.9g/cm3 

 elseif (fr > 0.e0 .and. fr < 1.e0) then  ! rime fraction 0~1

   xs_dry = xs(ibin) * ( 1.e0 - fr )   ! non-rimed dry snow mass per particle [g]

   if ( xs_dry < xs(1)) then  ! less than the smallest mass bin

       dens_dry = dens(1)  ![g/cm3]

   else ! interpolate

      do i = 1, nbin-1  !bin loop
         if ( xs_dry >= xs(i) .and. xs_dry <=  xs(i+1)) then
            dens_dry = dens(i) + (dens(i+1) - dens(i)) / &
                        (xs(i+1) - xs(i)) * (xs_dry - xs(i))  !density of dry snow component
            exit
         endif
      enddo

   endif

   !
   ! get total bulk diagnostic density considering riming fraction from dens_dry
   !
   dens_rime = dens_dry * rho_rime / (rho_rime*(1.e0 - fr ) + dens_dry * fr ) ! [g/cm3] 

 else

     print*,'MSG: diag_snow_SBM, strange riming fraction frime=',fr
     call stop_sdsu(myrank,'MSG: diag_snow_SBM: terminate program') 

 endif 

!
! update bulk radius considering riming fraction
!
  rs_rime = rs(ibin) * (dens(ibin)/dens_rime)**(1.e0/3.e0)  ! [cm]


 return
 end subroutine diag_rime_snow

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine diag_rime_melt_snow (ibin, fr, fl, dens, xs, rs, &
                                 rho_dry_rime, r_bulk, fv_liq  )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Diagnose bulk density and radius for rimed aggregate particle. 
! 
! History:
!   09/2017  Toshi Matsui@NASA GSFC : bug fixed for completely melted snow aggregate.
!   06/2011  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)        :: ibin  ! bin #
 real(sdsu_fps),intent(in) :: fr    ! riming mass fraction [-]
 real(sdsu_fps),intent(in) :: fl    ! melting mass fraction [-]
 real(sdsu_fps),intent(in) :: dens(1:nbin) ! snow bulk density from LUT (without riming) [g/cm3] 
 real(sdsu_fps),intent(in) :: xs(1:nbin)   ! snow bulk mass per particle from LUT (without riming) [g]
 real(sdsu_fps),intent(in) :: rs(1:nbin)   ! snow bulk radius per particle from LUT (without riming) [cm]

 real(sdsu_fps),intent(out) :: rho_dry_rime ! diagnostic bulk (effective) density of snow witm dry and riming. [g/cm3]
 real(sdsu_fps),intent(out) :: r_bulk  ! bulk radius [cm]
 real(sdsu_fps),intent(out) :: fv_liq  ! volume fraction of melted component [-]


! local params
 real(sdsu_fps) :: xs_dry   !snow bulk mass without riming and melting [g]
 real(sdsu_fps) :: rho_dry !snow bulk density without riming and melting [g/cm3]
 real(sdsu_fps),parameter :: rho_rime = 0.9e0    !rimed density  [g/cm3]
 real(sdsu_fps),parameter :: rho_liq  = 1.0e0    !liquid density  [g/cm3]
 real(sdsu_fps),parameter :: rho_ice  = 0.917e0  !pure ice density  [g/cm3]

 real(sdsu_fps) :: v_dry    !volume of dry component [cm3]
 real(sdsu_fps) :: v_rime   !volume of rimed component [cm3]
 real(sdsu_fps) :: v_liq    !volume of melted component [cm3]
 real(sdsu_fps) :: fd       !mass dry fraction [-]

 integer :: i !looping index

!
! initial check of values
!
 if( fr > 1.0 .or. fr < 0. ) then
   print*,'MSG: diag_rime_melt_snow, strange riming fraction frime=',fr
   stop
 elseif( fl > 1.0 .or. fl < 0. ) then
   print*,'MSG: diag_rime_melt_snow, strange melting fraction fmelt=',fl
   stop
 elseif( (fl + fr) > 1.0 ) then
   print*,'MSG: diag_rime_melt_snow, total melting and riming fraction > 1 ',fr+fl
   stop
 endif

 fd = max(0., 1.-fr-fl)  !dry component. 

!
! if completely melted....
!
 if( fl == 1. ) then  !completely melted
   v_liq  = xs(ibin) / rho_liq  !volume of melted component [cm3]
   rho_dry_rime = 1.0   !dry and rimed bulk density [g/cm3]
   r_bulk =  ( 3.*v_liq / (4.*const_pi) )**(1./3.) !bulk radius [cm]
   fv_liq = 1.0  
   return
 endif
 

!
! Determine dry snow aggregate density  (rho_dry)
! by taking account for riming and melting effect upon bulk snow aggregate density and radius. 
!
   xs_dry = xs(ibin) * fd   ! non-rimed and non-melted dry snow mass per particle [g]

   if ( xs_dry < xs(1)) then  ! less than the smallest mass bin

       rho_dry = dens(1)  ![g/cm3]

   else ! interpolate

      do i = 1, nbin-1  !bin loop
         if ( xs_dry >= xs(i) .and. xs_dry <=  xs(i+1)) then
            rho_dry = dens(i) + (dens(i+1) - dens(i)) / &
                        (xs(i+1) - xs(i)) * (xs_dry - xs(i))  !density of dry snow component [g/cm3]
            exit
         endif
      enddo

   endif


!
! get total bulk diagnostic density considering riming fraction from rho_dry
!

   rho_dry_rime = rho_dry * rho_rime / ( rho_rime * fd + rho_dry * fr ) ! [g/cm3] 

   v_dry  = fd*xs(ibin) / rho_dry  !volume of dry component [cm3]
   v_rime = fr*xs(ibin) / rho_rime !volume of rimed component [cm3]
   v_liq  = fl*xs(ibin) / rho_liq  !volume of melted component [cm3]

   rho_dry_rime = xs(ibin) * ( fd + fr ) / ( v_dry + v_rime )  !dry and rimed bulk density [g/cm3]


!
! output
!
   fv_liq  = v_liq  / (v_dry+v_rime+v_liq)  ! volumetric fraction of melted component [-]

   r_bulk = 0.5*( ( 6.e0/3.14159265 * xs(ibin) * (fd/rho_dry+fl/rho_liq+fr/rho_rime)  )**(1.e0/3.e0) )  ! bulk radius[cm]

   if(fv_liq < 0. .or. fv_liq > 1.0) call stop_sdsu(myrank,'MSG diag_rime_melt_snow: strange fv_liq')
   if( r_bulk < 0. )        call stop_sdsu(myrank,'MSG diag_rime_melt_snow: strange r_bulk' )
   if( rho_dry_rime < 0. )  call stop_sdsu(myrank,'MSG diag_rime_melt_snow: strange rho_dry_rime')


 return
 end subroutine diag_rime_melt_snow

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine diag_fv_melt( fl , rho_dry , fv_liq )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Diagnose volumetric fraction of liquid in grauple or hail. 
! 
! History:
!   08/2013  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: fl    ! melting mass fraction [-]
 real(sdsu_fps),intent(in) :: rho_dry ! density of dry graupel/hail [g/cm3] 
 real(sdsu_fps),intent(out) :: fv_liq  ! volume fraction of melted component [-]

! local params
 real(sdsu_fps),parameter :: rho_liq  = 1.0e0    !liquid density  [g/cm3]
 real(sdsu_fps) :: fd       !mass dry fraction [-]

!
! check input
!
 if( fl > 1.0 .or. fl < 0. ) then
   print*,'MSG: diag_fv_melt, strange melting fraction fl=',fl
   call stop_sdsu(myrank,'')
 endif
 fd = 1.e0 - fl  !dry fraction

 fv_liq =( fl/rho_liq ) / ( fl/rho_liq + fd/rho_dry )

 if(fv_liq < 0. .or. fv_liq > 1.0) call  stop_sdsu(myrank, 'MSG diag_fv_melt: strange fv_liq' )


 return
 end subroutine diag_fv_melt

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine rd_CRM_WRF
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Read WRF output (wrfout) in in netCDF format. All microphysics schemes are one-moment bulk.
!
! History:
! 11/2009  Toshi Matsui@NASA GSFC : Modified into paralell input. 
! 05/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
! WRF team, User's guide for Advanced Research WRF (ARW) modeling system version 2.2.
!            
!---------------------------------------------------------------------------------------------------
  
 integer :: i,j,k,n  !loop
 integer :: ncid,ncidsbm          !# of variables to inquire  
 integer :: varid               !# of variable id
 integer :: ierr                !allocation stat
 real(sdsu_fps)    :: press_top                              !model top pressure at interface [Pa] 
 real(sdsu_fps)    :: e, es  !unsaturated and saturated water vapor pressure [kPa]
 integer :: itemp  
 integer :: inet
 real(sdsu_fps) :: net
 integer :: is,ie,js,je,ks,ke,di,dj,dk
 real(sdsu_fps),allocatable :: wind_u(:,:)   !surface u wind speed [m/s]
 real(sdsu_fps),allocatable :: wind_v(:,:)   !surface v wind speed [m/s]
 real(sdsu_fps),allocatable :: rho(:,:,:)    !dry air density [kg/m3]
 real(sdsu_fps),allocatable :: net2d(:,:)
 real(sdsu_fps),allocatable :: net3d(:,:,:)
 real(sdsu_fps),allocatable :: temp3d(:,:,:)
 real(sdsu_fps),allocatable :: net3d_stag(:,:,:)
 integer,allocatable :: inet2d(:,:)
 character(len=2) :: char_bin, tag_char
 character(len=50) :: para_char   
 real(sdsu_fps) :: x,r
 real(sdsu_fps) :: q_bin  !mixing ratio per one bin [g/m3]
 real(sdsu_fps) :: n_bin !drop size ditributions [1/m3]

 real(sdsu_fps) :: dD, D, Nt, dr, V
 real(sdsu_fps) :: volume  ! particle volume [cm3/m3]
 real(sdsu_fps) :: reduce
 real(sdsu_fps) :: noise
 logical :: missing_para
 real(sdsu_fps) :: rime_mass


!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1 

!
! memory allocation
!
 if( .not. allocated(net2d) ) then
 allocate( &
           net2d     (is:ie,js:je          ), &
           inet2d    (is:ie,js:je          ), &
           net3d     (is:ie,js:je,1:mxlyr  ), &
           temp3d    (is:ie,js:je,1:mxlyr  ), &
           net3d_stag(is:ie,js:je,1:mxlyr+1), &
           wind_u    (is:ie,js:je          ), &
           wind_v    (is:ie,js:je          ), &
           rho       (is:ie,js:je,1:mxlyr  ), &
           stat=ierr )

 if(ierr /= 0) call stop_sdsu(myrank,'MSG rd_CRM_WRF: cannot allocate the variable')
 endif

!
! file name
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name)
 if(masterproc) print*,'MSG rd_CRM_WRF: ++++++++++++++++++++ READING NEW FILE ++++++++++++++++++++ '
#if MPI == 2
if(masterproc) &
#endif
 print*,'Input ->',trim(sdsu_io_file)

!
! open netCDF file
!
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! land mask
!
 call check( nf90_inq_varid(ncid, 'XLAND' , varid ) )  ! land-water mask 
 call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = myj_start, myj_end ; do i = myi_start, myi_end
    surface(i,j)%iland = INT( net2d(i,j) ) ! 1-land surface 2-ocean
    if( surface(i,j)%iland == 2 .or. surface(i,j)%iland == 1 ) then
    else
       surface(i,j)%iland = 2  !
    endif
 enddo ; enddo

!
! height and all sorts
!
 call check( nf90_inq_varid(ncid, 'PHB' , varid ) ) ! geopotential [m/s2]  (staggered height) 
 call check( nf90_get_var(ncid, varid, net3d_stag(is:ie,js:je,ks:ke+1), start=(/is,js,ks/), count=(/di,dj,dk+1/)))
 do k = 0, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos_stag(i,j,k)%hgt = net3d_stag(i,j,k+1)/const_g/1.e3  ![km] <- [m] (hydrostatic assumption)
 enddo ; enddo ; enddo

 do k = 0, mxlyr
    hgt_lev(k) = atmos_stag(myi_start,myj_start,k)%hgt      !toshi temporal [km] (z=0: bottom, z=mxlyr: top)
 enddo

!
! wind and all sorts
!
 call check( nf90_inq_varid(ncid, 'U10' , varid ) ) !surface u wind [m/s]
 call check( nf90_get_var(ncid, varid, wind_u(is:ie,js:je), start=(/is,js/), count=(/di,dj/))) 
 call check( nf90_inq_varid(ncid, 'V10' , varid ) ) !surface v wind [m/s]
 call check( nf90_get_var(ncid, varid, wind_v(is:ie,js:je), start=(/is,js/), count=(/di,dj/))) 
 do j = myj_start, myj_end ; do i = myi_start, myi_end
    surface(i,j)%u10m = sqrt(wind_u(i,j)**2. + wind_v(i,j)**2.)
 enddo ; enddo


 call check( nf90_inq_varid(ncid, 'W' , varid ) )   !vertical wind [m/s] (staggered height)
 call check( nf90_get_var(ncid, varid, net3d_stag(is:ie,js:je,ks:ke+1), start=(/is,js,ks/), count=(/di,dj,dk+1/))) 
 do k = 0, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos_stag(i,j,k)%w  = net3d_stag(i,j,k+1)
 enddo ; enddo ; enddo

!actually this is staggered, but reading as un-staggered....

!        float U(Time, bottom_top, south_north, west_east_stag) ;
 call check( nf90_inq_varid(ncid, 'U' , varid ) )   !u wind [m/s] (staggered height)
 call check( nf90_get_var(ncid, varid, atmos(is:ie,js:je,ks:ke)%u, start=(/is,js,ks/), count=(/di,dj,dk/)))

!        float V(Time, bottom_top, south_north_stag, west_east) ;
 call check( nf90_inq_varid(ncid, 'V' , varid ) )   !v wind [m/s] (staggered height)
 call check( nf90_get_var(ncid, varid, atmos(is:ie,js:je,ks:ke)%v, start=(/is,js,ks/), count=(/di,dj,dk/)))

!
! pressure and all sorts
!
 call check( nf90_inq_varid(ncid, 'PSFC' , varid ) ) !surface pressure [Pa]
 call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos_stag(i,j,0)%press = net2d(i,j) / 100.e0     ! [hPa] <- [Pa]
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'P_TOP' , varid ) ) !model top pressure [Pa]
 call check( nf90_get_var(ncid, varid, press_top ) )
 do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos_stag(i,j,mxlyr)%press = press_top/100.e0  ! [hPa] <- [Pa]
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'PB' , varid ) ) ! Use base pressure [Pa]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos(i,j,k)%press = net3d(i,j,k) / 100.e0   ! [hPa] <- [Pa]
 enddo ; enddo ; enddo

 do k = 1, mxlyr-1 ; do j = myj_start, myj_end ; do i = myi_start, myi_end  !toshii may change later...
   atmos_stag(i,j,k)%press = ( atmos(i,j,k)%press + atmos(i,j,k+1)%press ) * 0.5e0 !interafce <- average [hPa]
 enddo ; enddo ; enddo

!
! temperature and all sorts
!
 call check( nf90_inq_varid(ncid, 'TSK' , varid ) )  ! surface skin temperature [K]
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%t_skin, start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'T2' , varid ) )   ! 2-m air temp [K]
 call check( nf90_get_var(ncid, varid, atmos_stag(is:ie,js:je,0)%t_air, start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'T' , varid ) )    !purturvation potential temperature [K]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos(i,j,k)%t_air = (net3d(i,j,k)+300.e0) / (1.e3/(atmos(i,j,k)%press))**0.286e0 !air temperature [K]
 enddo ; enddo ; enddo

 do k = 1, mxlyr-1 ; do j = myj_start, myj_end ; do i = myi_start, myi_end  !toshii may change later
    atmos_stag(i,j,k)%t_air = 0.5e0 * ( atmos(i,j,k)%t_air + atmos(i,j,k+1)%t_air ) !interface <- average
 enddo ; enddo ; enddo

 do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos_stag(i,j,mxlyr)%t_air = atmos_stag(i,j,mxlyr-1)%t_air &
                       - (atmos_stag(i,j,mxlyr-1)%t_air-atmos(i,j,mxlyr)%t_air)*2.e0  !air temperature for top level
 enddo ; enddo

 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     rho(i,j,k) = (atmos(i,j,k)%press*100.e0) / (const_Rd*atmos(i,j,k)%t_air)  !dry air density [kg/m3]
     atmos(i,j,k)%rho_dair = rho(i,j,k) * 1.e3 ! [g/m3]
 enddo ; enddo ; enddo 

!
! humidity and all sorts
!
 call check( nf90_inq_varid(ncid, 'QVAPOR' , varid ) ) !water vapor mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))

 atmos%qv = net3d  ! water vapor mixing ratio [kg/kg]

!print*, net3d(1,1,:)
!stop
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    es = 0.611e0 * exp( const_Lv_Rv*( 1.e0/const_Kel2Cel - 1.e0/atmos(i,j,k)%t_air ) )  ! sat vapor pressure [kPa]
    e = (net3d(i,j,k) / (net3d(i,j,k)+const_Rd_Rv)) * atmos(i,j,k)%press * 0.1e0  ! vapor pressure [mb] -> [kPa]
    atmos(i,j,k)%rh = max( 1e-5, e/es) *100.e0     ! relative humidity [%]
 enddo ; enddo ; enddo

!
! hydrometeors
!
 mic_select0: select case(trim(cloud_microphysics))

 case('GOD','GOD10','GMP4ICE','LIN','WSM','WDM','MORR','MORRH','THOM','NTUCLR','NSSL2MMT')

 !mixing ratio
 call check( nf90_inq_varid(ncid, 'QCLOUD' , varid ) ,'QCLOUD',missing_para) !cloud water mixing ratio[kg/kg]
 if(missing_para) then
   q_gen%cloud = 0.e0
 else
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gen(i,j,k)%cloud = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
 enddo ; enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'QRAIN' , varid ) ,'QRAIN',missing_para) !rain mixing ratio[kg/kg]
 if(missing_para) then
   q_gen%rain = 0.e0
 else
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gen(i,j,k)%rain = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
 enddo ; enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'QICE' , varid ) ,'QICE',missing_para) !cloud ice mixing ratio[kg/kg]
 if(missing_para) then
   q_gen%ice = 0.e0
 else
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gen(i,j,k)%ice = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
 enddo ; enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'QSNOW' , varid ) ,'QSNOW',missing_para) !snow mixing ratio[kg/kg]
 if(missing_para) then
   q_gen%snow = 0.e0
 else
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gen(i,j,k)%snow = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
 enddo ; enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'QGRAUP' , varid ) ,'QGRAUP',missing_para) !graupel mixing ratio[kg/kg]
 if(missing_para) then
   q_gen%graupel = 0.e0
 else

 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    if(gsfc_hail .or. trim(cloud_microphysics) == 'MORRH') then  ! hail option --> read graupel as hail
       q_gen(i,j,k)%hail    = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
       q_gen(i,j,k)%graupel = 0.e0
       if(net3d(i,j,k)> 0.) print*,i,j,k,'reading hail mass', q_gen(i,j,k)%hail
    else  ! reading graupel as graupel
       q_gen(i,j,k)%graupel = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
       q_gen(i,j,k)%hail = 0.e0
    endif
 enddo ; enddo ; enddo
 endif

 if(gsfc_hail .or. trim(cloud_microphysics) == 'MORRH') then  ! hail option --> read graupel as hail
 else

   call check( nf90_inq_varid(ncid, 'QHAIL' , varid ), 'QHAIL',missing_para) !hail mixing ratio[kg/kg]
   if(missing_para) then
     q_gen%hail = 0.e0
   else
     call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
     do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
        q_gen(i,j,k)%hail    = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
     enddo ; enddo ; enddo
    endif
 endif

!
! NSSL electric microphysics only parameter 
!
 call check( nf90_inq_varid(ncid, 'SCTOT' , varid ) ,'SCTOT',missing_para) !total space charge density [C/m3]
 if(missing_para) then
   atmos%sctot = 0.e0
 else
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos(i,j,k)%sctot = net3d(i,j,k) ![C/m3]
 enddo ; enddo ; enddo
 endif


!
! No hydrometero option (toshii)
!
! q_gen%cloud = 0. 
! q_gen%rain  = 0.
! q_gen%ice = 0. 
! q_gen%snow = 0.  
! q_gen%graupel = 0. 
! q_gen%hail = 0. 
!


 WDM_IF: if ( trim(cloud_microphysics) == 'WDM' ) then

  call check( nf90_inq_varid(ncid, 'QNRAIN' , varid ) ) !rain number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%rain = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNCLOUD' , varid ) ) !cloud number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%cloud = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  n_gen(:,:,:)%ice     = undefined ![g/m3] < - [kg/kg]
  n_gen(:,:,:)%snow    = undefined ![g/m3] < - [kg/kg]
  n_gen(:,:,:)%graupel = undefined ![g/m3] < - [kg/kg]
  n_gen(:,:,:)%hail    = undefined ![g/m3] < - [kg/kg]


 call check( nf90_inq_varid(ncid, 'RE_CLOUD' , varid ), 'RE_CLOUD',missing_para) !cloud effective radius [m]
 if( missing_para ) then
    wdm_re_cloud_exist = .false.
   re_gen%cloud = undefined
 else
    wdm_re_cloud_exist = .true.

  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     re_gen%cloud = net3d(i,j,k) * 1.e+6 ![um] < - [m]
  enddo ; enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'RE_ICE' , varid ), 'RE_ICE',missing_para) !ice effective radius [m]
 if( missing_para ) then
    wdm_re_ice_exist = .false.
    re_gen%ice = undefined
 else
    wdm_re_ice_exist = .true.
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     re_gen(i,j,k)%ice = net3d(i,j,k) * 1.e+6 ![um] < - [m]
  enddo ; enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'RE_SNOW' , varid ), 'RE_SNOW',missing_para) !snow effective radius [m]
 if( missing_para ) then
    re_gen%snow = undefined
    wdm_re_snow_exist = .false.
 else
    wdm_re_snow_exist = .true.

  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     re_gen(i,j,k)%snow = net3d(i,j,k) * 1.e+6 ![um] < - [m]
  enddo ; enddo ; enddo
 endif


 endif WDM_IF

 THOM_IF: if ( trim(cloud_microphysics) == 'THOM' ) then

  call check( nf90_inq_varid(ncid, 'QNRAIN' , varid ) ) !rain number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%rain = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNICE' , varid ) ) !cloud ice number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%ice = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  n_gen(:,:,:)%cloud   = undefined ![g/m3] < - [kg/kg]
  n_gen(:,:,:)%snow    = undefined ![g/m3] < - [kg/kg]
  n_gen(:,:,:)%graupel = undefined ![g/m3] < - [kg/kg]
  n_gen(:,:,:)%hail    = undefined ![g/m3] < - [kg/kg]

 endif THOM_IF


 MORR_IF: if ( trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH' ) then

  !number conc
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end  !cloud number conc [#/kg]
     n_gen(i,j,k)%cloud = undefined ![g/m3] < - [kg/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNRAIN' , varid ) ) !rain number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%rain = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNICE' , varid ) ) !cloud ice number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%ice = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNSNOW' , varid ) ) !snow number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%snow = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNGRAUPEL' , varid ) ) !graupel number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     if( trim(cloud_microphysics) == 'MORRH' ) then
        n_gen(i,j,k)%hail = net3d(i,j,k)*rho(i,j,k) ! hail number conc [#/m3]
        n_gen(i,j,k)%graupel = 0.
        if(net3d(i,j,k)> 0.) print*,i,j,k,'reading hail conc',n_gen(i,j,k)%hail
     else
        n_gen(i,j,k)%graupel = net3d(i,j,k)*rho(i,j,k) ! graupel number conc [#/m3]
        n_gen(i,j,k)%hail = 0.
     endif
  enddo ; enddo ; enddo

 endif  MORR_IF

 NTUCLR_IF: if ( trim(cloud_microphysics) == 'NTUCLR' ) then


  call check( nf90_inq_varid(ncid, 'QNRAIN' , varid ) ) !rain number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%rain = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNCLOUD' , varid ) ) !cloud number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%cloud = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNICE' , varid ) ) !ice number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%ice = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  n_gen(:,:,:)%snow    = undefined ![g/m3] < - [kg/kg]
  n_gen(:,:,:)%graupel = undefined ![g/m3] < - [kg/kg]
  n_gen(:,:,:)%hail    = undefined ![g/m3] < - [kg/kg]

 endif NTUCLR_IF

 NSSL_IF: if ( trim(cloud_microphysics) == 'NSSL2MMT' ) then

  call check( nf90_inq_varid(ncid, 'QNDROP' , varid ) ) !cloud number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%cloud = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNRAIN' , varid ) ) !rain number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%rain = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNICE' , varid ) ) !ice number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%ice = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNSNOW' , varid ) ) !snow number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%snow = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNGRAUPEL' , varid ) ) !graupel number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%graupel = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNHAIL' , varid ) ) !hail number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%hail = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QVGRAUPEL' , varid ) ) !graupel volume conc [m3/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     if( net3d(i,j,k) > 0. ) then
         rho3d_graupel(i,j,k) = q_gen(i,j,k)%graupel / (net3d(i,j,k) * rho(i,j,k)) / 1000. !3D graupel density [kg/m3] 
         rho3d_graupel(i,j,k) = min(900., max(170., rho3d_graupel(i,j,k) ) )  !confine min and max range
     else
         !make sure q, n, and rho are all zero. 
         rho3d_graupel(i,j,k) = 0. 
         q_gen(i,j,k)%graupel = 0. 
         n_gen(i,j,k)%graupel = 0.
     endif
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QVHAIL' , varid ) ) !hail volume conc [m3/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     if( net3d(i,j,k) > 0. ) then
         rho3d_hail(i,j,k) = q_gen(i,j,k)%hail / (net3d(i,j,k) * rho(i,j,k)) / 1000. !3D hail density [kg/m3] 
         rho3d_hail(i,j,k) = min(900., max(500., rho3d_hail(i,j,k) ))  !confine min and max range
     else
         !make sure q, n, and rho are all zero.
         rho3d_hail(i,j,k) = 0.
         q_gen(i,j,k)%hail = 0. 
         n_gen(i,j,k)%hail = 0.
     endif

  enddo ; enddo ; enddo


 endif NSSL_IF

 case('P3') 

 !mixing ratio
 call check( nf90_inq_varid(ncid, 'QCLOUD' , varid ) ,'QCLOUD',missing_para) !cloud water mixing ratio[kg/kg]
 if(missing_para) then
   q_gen%cloud = 0.e0
 else
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gen(i,j,k)%cloud = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
 enddo ; enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'QRAIN' , varid ) ,'QRAIN',missing_para) !rain mixing ratio[kg/kg]
 if(missing_para) then
   q_gen%rain = 0.e0
 else
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gen(i,j,k)%rain = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
 enddo ; enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'QICE' , varid ) ,'QICE',missing_para) !cloud ice mixing ratio[kg/kg]
 if(missing_para) then
   q_gen%ice = 0.e0
 else
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gen(i,j,k)%ice = net3d(i,j,k) * rho(i,j,k) * 1.e3 ![g/m3] < - [kg/kg]
 enddo ; enddo ; enddo
 endif

!toshii debbug for P3
!    q_gen%cloud  = 0.

 q_gen%snow = 0. ; q_gen%graupel = 0. ; q_gen%hail = 0. 

  call check( nf90_inq_varid(ncid, 'QNCLOUD' , varid ) ) !cloud number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%cloud = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNRAIN' , varid ) ) !rain number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%rain = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QNICE' , varid ) ) !ice number conc [#/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     n_gen(i,j,k)%ice = net3d(i,j,k) * rho(i,j,k) ![#/m3] < - [#/kg]
  enddo ; enddo ; enddo

  n_gen%snow = 0. ; n_gen%graupel = 0. ; n_gen%hail = 0.


  call check( nf90_inq_varid(ncid, 'QIR' , varid ) ) !Rime ice mass-1 mixing ratio [kg/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     if( q_gen(i,j,k)%ice == 0. ) then
        p3_frime(i,j,k) = 0. 
     else
    !     [-]                      [kg/m3]                        [kg/m3]
        p3_frime(i,j,k) = (net3d(i,j,k) * rho(i,j,k))/( q_gen(i,j,k)%ice*1.e-3) ![kg/m3] < - [kg/kg]
     endif
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'QIB' , varid ) ) !Rime ice volume-1 mixing ratio [m3/kg]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     if ( net3d(i,j,k) == 0. ) then
        p3_rhorime(i,j,k) = 0. 
     else
        rime_mass = p3_frime(i,j,k) * ( q_gen(i,j,k)%ice*1.e-3)  ! [kg/m3]
   !        [kg/m3]          [kg/m3]               [m3/m3] 
!        p3_rhorime(i,j,k) = rime_mass / (net3d(i,j,k) * rho(i,j,k)) ! got division by zero...
        p3_rhorime(i,j,k) = rime_mass / max(const_smallest_real,  (net3d(i,j,k) * rho(i,j,k)) ) !got division by zero...
     endif
  enddo ; enddo ; enddo


 case('RAMS1','RAMS2')
   call stop_sdsu(myrank, 'MSG rd_CRM_WRF: Not yet supporting RAMS microphysics in WRF input')
 case('HUCM_SBM')

 frime_sbm_snow(:,:,:,:)     = 0.e0
 fmelt_sbm_snow(:,:,:,:)     = 0.e0
 fmelt_sbm_graupel(:,:,:,:)     = 0.e0
 fmelt_sbm_hail(:,:,:,:)     = 0.e0

 !
 ! liquid  (cloud and rain)
 !
 n_sbm%liq = 0.e0 ; q_sbm%liq = 0.e0
 tag_char = 'dr'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist2'//trim(sdsu_inp_name(efile_len-23:efile_len))
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char  
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD function [/ (cm3 g)]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function(x_sbm(n)%liq, drad_sbm(n)%liq*1.e-2, net3d(i,j,k), &
                         n_sbm(i,j,k,n)%liq, q_sbm(i,j,k)%liq )
    enddo ; enddo ; enddo 
 enddo   
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! ice column
 !
 n_sbm%ice_col = 0.e0 ; q_sbm%ice_col = 0.e0
 tag_char = 'ic'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist3'//trim(sdsu_inp_name(efile_len-23:efile_len))
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char  
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD function [/ (cm3 g)]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function(x_sbm(n)%ice_col, drad_sbm(n)%ice_col*1.e-2, net3d(i,j,k), &
                         n_sbm(i,j,k,n)%ice_col, q_sbm(i,j,k)%ice_col )
    enddo ; enddo ; enddo
 enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! ice plate
 !
 n_sbm%ice_pla = 0.e0 ; q_sbm%ice_pla = 0.e0
 tag_char = 'ip'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist4'//trim(sdsu_inp_name(efile_len-23:efile_len))
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char  
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD function [/ (cm3 g)]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function(x_sbm(n)%ice_pla, drad_sbm(n)%ice_pla*1.e-2, net3d(i,j,k), &
                         n_sbm(i,j,k,n)%ice_pla, q_sbm(i,j,k)%ice_pla )
    enddo ; enddo ; enddo
 enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! ice dendrite
 !
 n_sbm%ice_den = 0.e0 ; q_sbm%ice_den = 0.e0
 tag_char = 'id'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist5'//trim(sdsu_inp_name(efile_len-23:efile_len))
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char 
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD function [/ (cm3 g)]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function(x_sbm(n)%ice_den, drad_sbm(n)%ice_den*1.e-2, net3d(i,j,k), &
                         n_sbm(i,j,k,n)%ice_den, q_sbm(i,j,k)%ice_den )
    enddo ; enddo ; enddo
 enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! snow aggregate 
 !
 n_sbm%snow = 0.e0 ; q_sbm%snow = 0.e0
 tag_char = 'sn'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist7'//trim(sdsu_inp_name(efile_len-23:efile_len))
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char 
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD function [/ (cm3 g)]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function(x_sbm(n)%snow, drad_sbm(n)%snow*1.e-2, net3d(i,j,k), &
                         n_sbm(i,j,k,n)%snow, q_sbm(i,j,k)%snow )
   enddo ; enddo ; enddo
 enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! graupel
 !
 n_sbm%graupel = 0.e0 ; q_sbm%graupel = 0.e0
 tag_char = 'gr'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist8'//trim(sdsu_inp_name(efile_len-23:efile_len))
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char 
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD function [/ (cm3 g)] 
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function(x_sbm(n)%graupel, drad_sbm(n)%graupel*1.e-2, net3d(i,j,k), &
                         n_sbm(i,j,k,n)%graupel, q_sbm(i,j,k)%graupel )
    enddo ; enddo ; enddo
 enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! hail
 !
 n_sbm%hail = 0.e0 ; q_sbm%hail = 0.e0
 tag_char = 'ha'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist9'//trim(sdsu_inp_name(efile_len-23:efile_len))
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char 
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD function [/ (cm3 g)]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function(x_sbm(n)%hail, drad_sbm(n)%hail*1.e-2, net3d(i,j,k), &
                         n_sbm(i,j,k,n)%hail, q_sbm(i,j,k)%hail )
    enddo ; enddo ; enddo
 enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 case('HUCM_SBM43')

 !
 ! liquid  (cloud and rain)
 !
 n_sbm%liq = 0.e0 ; q_sbm%liq = 0.e0
 tag_char = 'dr'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist2'//trim(sdsu_inp_name(efile_len-23:efile_len))
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input (liquid) ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char  
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )    ! SBM PSD mass mixing ratio [kg/kg]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function43(rho(i,j,k), x_sbm(n)%liq, drad_sbm(n)%liq*1.e-2, net3d(i,j,k), &
                           n_sbm(i,j,k,n)%liq, q_sbm(i,j,k)%liq )
    enddo ; enddo ; enddo 
 enddo   
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! ice column
 !
 n_sbm%ice_col = 0.e0 ; q_sbm%ice_col = 0.e0
 tag_char = 'ic'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist3'//trim(sdsu_inp_name(efile_len-23:efile_len))
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input (ice column)->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char  
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg] 
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function43(rho(i,j,k), x_sbm(n)%ice_col, drad_sbm(n)%ice_col*1.e-2, net3d(i,j,k), &
                           n_sbm(i,j,k,n)%ice_col, q_sbm(i,j,k)%ice_col )
    enddo ; enddo ; enddo
 enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! ice plate
 !
 n_sbm%ice_pla = 0.e0 ; q_sbm%ice_pla = 0.e0
 tag_char = 'ip'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist4'//trim(sdsu_inp_name(efile_len-23:efile_len))
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input (ice plate)->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char  
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg] 
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function43(rho(i,j,k), x_sbm(n)%ice_pla, drad_sbm(n)%ice_pla*1.e-2, net3d(i,j,k), &
                           n_sbm(i,j,k,n)%ice_pla, q_sbm(i,j,k)%ice_pla )
    enddo ; enddo ; enddo
 enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! ice dendrite
 !
 n_sbm%ice_den = 0.e0 ; q_sbm%ice_den = 0.e0
 tag_char = 'id'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist5'//trim(sdsu_inp_name(efile_len-23:efile_len))
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input (ice dendrite) ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char 
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function43(rho(i,j,k), x_sbm(n)%ice_den, drad_sbm(n)%ice_den*1.e-2, net3d(i,j,k), &
                           n_sbm(i,j,k,n)%ice_den, q_sbm(i,j,k)%ice_den )
    enddo ; enddo ; enddo
 enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! snow aggregate 
 !
 n_sbm%snow = 0.e0 ; q_sbm%snow = 0.e0
 tag_char = 'sn'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist7'//trim(sdsu_inp_name(efile_len-23:efile_len))
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input (snow aggregate) ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char   ! total component
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function43(rho(i,j,k), x_sbm(n)%snow, drad_sbm(n)%snow*1.e-2, net3d(i,j,k), &
                           n_sbm(i,j,k,n)%snow, q_sbm(i,j,k)%snow )
    enddo ; enddo ; enddo

    write(char_bin,"(I2.2)") n  ; para_char = 'flc'//char_bin//tag_char   ! liquid (melted) component
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg]
    call check( nf90_get_var(ncidsbm, varid, temp3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       if( isnan(temp3d(i,j,k)) ) print*,'find nan in melt',i,j,k,n
       if( temp3d(i,j,k) > 0.e0 .and. net3d(i,j,k) > 0.e0 ) then
           fmelt_sbm_snow(i,j,k,n) = MAX( 0.e0, MIN( 1.0e0, temp3d(i,j,k)/net3d(i,j,k) ) )  !melt fraction [-]  (0~1)
       else
           fmelt_sbm_snow(i,j,k,n) = 0.e0
       endif
    enddo ; enddo ; enddo

    write(char_bin,"(I2.2)") n  ; para_char = 'rf3'//char_bin//tag_char   ! rime component
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg]
    call check( nf90_get_var(ncidsbm, varid, temp3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       if( isnan(temp3d(i,j,k)) ) print*,'find nan in rime frac',i,j,k,n
       if( temp3d(i,j,k) > 0.e0 .and. net3d(i,j,k) > 0.e0) then
           frime_sbm_snow(i,j,k,n) = MAX( 0.e0, MIN( 0.9e0, temp3d(i,j,k)/net3d(i,j,k)) ) !rime fraction [-]  (0~1)
       else
           frime_sbm_snow(i,j,k,n) = 0.e0
       endif
    enddo ; enddo ; enddo

    !
    ! adjust riming and melting fraction to remove numerical noise
    !
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       if( ( fmelt_sbm_snow(i,j,k,n) + frime_sbm_snow(i,j,k,n) ) > 1.0 ) then 
           noise = ( fmelt_sbm_snow(i,j,k,n) + frime_sbm_snow(i,j,k,n) ) - 1.e0
           frime_sbm_snow(i,j,k,n) = MAX( 0.e0, MIN( 0.9e0, frime_sbm_snow(i,j,k,n) - noise  ))  ! take out from melting....
       endif
    enddo ; enddo ; enddo


 enddo

 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! graupel
 !
 n_sbm%graupel = 0.e0 ; q_sbm%graupel = 0.e0
 tag_char = 'gr'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist8'//trim(sdsu_inp_name(efile_len-23:efile_len))
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input (graupel) ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

 ! total
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char  !total
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function43(rho(i,j,k), x_sbm(n)%graupel, drad_sbm(n)%graupel*1.e-2, net3d(i,j,k), &
                           n_sbm(i,j,k,n)%graupel, q_sbm(i,j,k)%graupel )
    enddo ; enddo ; enddo

    write(char_bin,"(I2.2)") n  ; para_char = 'flc'//char_bin//tag_char  !liquid (melted) component
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg]
    call check( nf90_get_var(ncidsbm, varid, temp3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       if( temp3d(i,j,k) > 0.e0 .and. net3d(i,j,k) > 0.e0) then 
           fmelt_sbm_graupel(i,j,k,n) = MAX( 0.e0, MIN( 1.0e0, temp3d(i,j,k)/net3d(i,j,k) ))  !melt fraction [-]  (0~1)
       else 
           fmelt_sbm_graupel(i,j,k,n) = 0.e0 
       endif
    enddo ; enddo ; enddo

 enddo

 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! hail
 !
 n_sbm%hail = 0.e0 ; q_sbm%hail = 0.e0
 tag_char = 'ha'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist9'//trim(sdsu_inp_name(efile_len-23:efile_len))
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input (hail) ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 do n = 1, nbin
    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char  !total portion
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg]
    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       call sbm_function43(rho(i,j,k), x_sbm(n)%hail, drad_sbm(n)%hail*1.e-2, net3d(i,j,k), &
                           n_sbm(i,j,k,n)%hail, q_sbm(i,j,k)%hail )
    enddo ; enddo ; enddo

    write(char_bin,"(I2.2)") n  ; para_char = 'flc'//char_bin//tag_char  ! liquid (melted) component
    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD mass mixing ratio [kg/kg]
    call check( nf90_get_var(ncidsbm, varid, temp3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
       if( temp3d(i,j,k) > 0.e0 .and. net3d(i,j,k) > 0.e0) then
           fmelt_sbm_hail(i,j,k,n) = MAX( 0.e0, MIN( 1.0e0, temp3d(i,j,k)/net3d(i,j,k) ) ) !melt fraction [-]  (0~1)
       else
           fmelt_sbm_hail(i,j,k,n) = 0.e0
       endif
    enddo ; enddo ; enddo

 enddo

! No hydrometero option (toshii)
 if(masterproc) print*,' ------------ remove hydrometero -----------------------'
! n_sbm%liq = 0.e0     ; q_sbm%liq = 0.e0
! n_sbm%ice_col = 0.e0 ; q_sbm%ice_col = 0.e0
! n_sbm%ice_pla = 0.e0 ; q_sbm%ice_pla = 0.e0
! n_sbm%ice_den = 0.e0 ; q_sbm%ice_den = 0.e0
! n_sbm%snow = 0.e0    ; q_sbm%snow = 0.e0
! n_sbm%graupel = 0.e0 ; q_sbm%graupel = 0.e0
! n_sbm%hail = 0.e0    ; q_sbm%hail = 0.e0

 call check( nf90_close(ncidsbm) )  ! close nc file

 case default
    ! do nothing
 end select mic_select0

!
! surface precipitation and all sorts
!
! call check( nf90_inq_varid(ncid, 'RAINNCV' , varid ) ) !time-step precip [mm]
! call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
! do j = myj_start, myj_end ; do i = myi_start, myi_end
!    surface(i,j)%rain_rate = net2d(i,j) / 3.33333333e0 * 3600.e0   ! [mm/hr] <- [mm]
! enddo ; enddo

! call check( nf90_inq_varid(ncid, 'RE_SNOW' , varid ), 'RE_SNOW',missing_para) !snow effective radius [m]

! call check( nf90_inq_varid(ncid, 'RE_SNOW' , varid ), 'RE_SNOW',missing_para) !snow effective radius [m]
! if( missing_para ) then
!    re_gen%snow = undefined
!    wdm_re_snow_exist = .false.
! else
!    wdm_re_snow_exist = .true.

!  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
!  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
!     re_gen(i,j,k)%snow = net3d(i,j,k) * 1.e+6 ![um] < - [m]
!  enddo ; enddo ; enddo
! endif

!
! surface parameters and all sorts
!
 call check( nf90_inq_varid(ncid, 'IVGTYP' , varid ), 'IVGTYP', missing_para )  ! dominant vegetation type
 if( missing_para ) then
     surface%igbp_typ = 0
 else
  call check( nf90_get_var(ncid, varid, inet2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
  do j = myj_start, myj_end ; do i = myi_start, myi_end
     call convert_usgs_igbp( inet2d(i,j) , surface(i,j)%igbp_typ)
  enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'VEGFRA' , varid ), 'VEGFRA' , missing_para )  ! vegetation fraction [%]
 if( missing_para ) then
     surface%frac_veg = 0.
 else
     call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%frac_veg, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'SNOW' , varid ), 'SNOW', missing_para )  ! snow water equivalent [kg m-2]
 if( missing_para ) then
     surface%h2o_snow = 0.
 else
    call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%h2o_snow, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'SNOWH' , varid ), 'SNOWH', missing_para )  ! snow depth [m]
 if( missing_para ) then
    surface%dhgt_snow = 0.
 else
   call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%dhgt_snow, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'XLAT' , varid ) )  ! latitude [deg]
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%lat, start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'XLONG' , varid ) )  ! longitude [deg]
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%lon, start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'ALBEDO' , varid ), 'ALBEDO' , missing_para )  ! BROADBAND ALBEDO  
 if( missing_para ) then
   surface%albedo = 0.05
 else
   call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%albedo, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'SMOIS' , varid ),  'SMOIS', missing_para )  ! soil moisture [m3 m-3]
 if( missing_para ) then
   surface%h2o_soil= 0.0
 else
   call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%h2o_soil, start=(/is,js,1/), count=(/di,dj,1/)))
 endif

 call check( nf90_inq_varid(ncid, 'TSLB' , varid ), 'TSLB', missing_para )  ! soil temperature [K]
 if( missing_para ) then
   surface%t_soil= atmos(:,:,1)%t_air
 else
   call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%t_soil, start=(/is,js,1/), count=(/di,dj,1/)))
 endif

 call check( nf90_inq_varid(ncid, 'HGT' , varid ),'HGT', missing_para )  ! surface elevation [m]  
 if( missing_para ) then
   surface%elev = 0.
 else
   call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%elev, start=(/is,js/), count=(/di,dj/)))
 endif


!
! GOCART aerosols
!
 if(account_aerosol) then

 select case(trim(aerosol_microphysics) )
 case('GOCART')
 
 call check( nf90_inq_varid(ncid, 'sulf' , varid ) ) ! SO4 [ppmv]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%so4 = net3d(i,j,k) *1.0e-6*(atmos(i,j,k)%press*100.e0)*96.0 &
                          /(8.314*atmos(i,j,k)%t_air ) ! SO4 [g/m3] 
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'BC1' , varid ) ) ! Black Carbon 1 [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%blc = net3d(i,j,k) * rho(i,j,k) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'BC2' , varid ) ) ! Black Carbon 2 [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%blc = q_gocart(i,j,k)%blc + net3d(i,j,k) * rho(i,j,k) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'OC1' , varid ) ) ! Organic Carbon 1 (non-hygroscopic) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%ocn = net3d(i,j,k) * rho(i,j,k) * 1.0e-6 * 1.4  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'OC2' , varid ) ) ! Organic Carbon 2 (hygroscopic) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%och = net3d(i,j,k) * rho(i,j,k) * 1.0e-6 * 1.4  ![g/m3]
 enddo ; enddo ; enddo


 call check( nf90_inq_varid(ncid, 'SEAS_1' , varid ) ) ! Sea Salt (fine mode: 0.3re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%ssa = net3d(i,j,k) * rho(i,j,k) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SEAS_2' , varid ) ) ! Sea Salt (coarse mode: 1.0re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%ssc = net3d(i,j,k) * rho(i,j,k) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SEAS_3' , varid ) ) ! Sea Salt (coarse mode: 3.2re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%ssc = q_gocart(i,j,k)%ssc + net3d(i,j,k) * rho(i,j,k) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SEAS_4' , varid ) ) ! Sea Salt (coarse mode: 7.5re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%ssc = q_gocart(i,j,k)%ssc + net3d(i,j,k) * rho(i,j,k) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 reduce = 1.  !sensitivity factor (normal is one.)

 call check( nf90_inq_varid(ncid, 'DUST_1' , varid ) ) ! Dust mode 1~4 (0.5re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%du1 = net3d(i,j,k) * rho(i,j,k) * 1.0e-6 * dust_frac(1)*reduce ![g/m3] 
    q_gocart(i,j,k)%du2 = net3d(i,j,k) * rho(i,j,k) * 1.0e-6 * dust_frac(2)*reduce ![g/m3] 
    q_gocart(i,j,k)%du3 = net3d(i,j,k) * rho(i,j,k) * 1.0e-6 * dust_frac(3)*reduce ![g/m3] 
    q_gocart(i,j,k)%du4 = net3d(i,j,k) * rho(i,j,k) * 1.0e-6 * dust_frac(4)*reduce ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DUST_2' , varid ) ) ! Dust mode 5 (1.4re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%du5 = net3d(i,j,k) * rho(i,j,k) * 1.0e-6 *reduce ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DUST_3' , varid ) ) ! Dust mode 6 (2.4re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%du6 = net3d(i,j,k) * rho(i,j,k) * 1.0e-6 *reduce ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DUST_4' , varid ) ) ! Dust mode 7 (4.5re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%du7 = net3d(i,j,k) * rho(i,j,k) * 1.0e-6 *reduce ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DUST_5' , varid ) ) ! Dust mode 8 (8.0re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    q_gocart(i,j,k)%du8 = net3d(i,j,k) * rho(i,j,k) * 1.0e-6  *reduce ![g/m3]
 enddo ; enddo ; enddo

! mask (temporal)
!    q_gocart%du1 = 0.
!    q_gocart%du2 = 0.
!    q_gocart%du3 = 0.
!    q_gocart%du4 = 0.
!    q_gocart%du5 = 0.
!    q_gocart%du6 = 0.
!    q_gocart%du7 = 0.
!    q_gocart%du8 = 0.

 case('UCD')

  call check( nf90_inq_varid(ncid, 'DUST_1' , varid ) ) ! Dust mode 1  [kg kg-1]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     q_ucd(i,j,k)%du1 = net3d(i,j,k) * rho(i,j,k) * 1.0e+3   ![g/m3]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'DUST_2' , varid ) ) ! Dust mode 2  [kg kg-1]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     q_ucd(i,j,k)%du2 = net3d(i,j,k) * rho(i,j,k) * 1.0e+3   ![g/m3]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'DUST_3' , varid ) ) ! Dust mode 3  [kg kg-1]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     q_ucd(i,j,k)%du3 = net3d(i,j,k) * rho(i,j,k) * 1.0e+3   ![g/m3]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'DUST_4' , varid ) ) ! Dust mode 4  [kg kg-1]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     q_ucd(i,j,k)%du4 = net3d(i,j,k) * rho(i,j,k) * 1.0e+3   ![g/m3]
  enddo ; enddo ; enddo

  call check( nf90_inq_varid(ncid, 'DUST_5' , varid ) ) ! Dust mode 5  [kg kg-1]
  call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
     q_ucd(i,j,k)%du5 = net3d(i,j,k) * rho(i,j,k) * 1.0e+3   ![g/m3]
  enddo ; enddo ; enddo

 case default
  call stop_sdsu(myrank,'MSG rd_CRM_WRF: There is no such aerosol_microphysics')
 end select

 endif  !account aerosols

!
!----------------  additional diagnostic parameters, if available -------------------------------------
!

!
! total atmospheric latent heat
!
 if ( trim(cloud_microphysics) == 'P3' ) then
   call check( nf90_inq_varid(ncid, 'LATHEAT' , varid ), 'LATHEAT',missing_para) !latent heat 
 else
   call check( nf90_inq_varid(ncid, 'PHYS_TOT' , varid ), 'PHYS_TOT',missing_para) !latent heat 
 endif
 if(missing_para) then
   atmos%latent_heat = undefined
 else
   call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
   do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
      atmos(i,j,k)%latent_heat    = net3d(i,j,k) * 24. * 3600. ![K/day] < - [degC/s]
   enddo ; enddo ; enddo
 endif

!
! cloud nucleation rate
!
 if ( trim(cloud_microphysics) == 'P3' ) then
   call check( nf90_inq_varid(ncid, 'CLDNUC' , varid ), 'CLDNUC',missing_para) !latent heat 

 if(missing_para) then
   atmos%cloud_nuc = undefined
 else
   call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
   do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
      atmos(i,j,k)%cloud_nuc = net3d(i,j,k)  ![kg/kg/s]
   enddo ; enddo ; enddo
 endif

 endif



!
! LW radiative heating
!
 call check( nf90_inq_varid(ncid, 'RTHRATLW' , varid ), 'RTHRATLW',missing_para) !LW radiative heating
 if(missing_para) then
   atmos%lw_heat = undefined
 else
   call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
   do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
      atmos(i,j,k)%lw_heat  = net3d(i,j,k) * 24. * 3600. ![K/day] < - [K/s]
   enddo ; enddo ; enddo
 endif

!
! SW radiative heating
!
 call check( nf90_inq_varid(ncid, 'RTHRATSW' , varid ), 'RTHRATSW',missing_para) !SW radiative heating
 if(missing_para) then
   atmos%sw_heat = undefined
 else
   call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
   do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
      atmos(i,j,k)%sw_heat  = net3d(i,j,k) * 24. * 3600. ![K/day] < - [K/s]
   enddo ; enddo ; enddo
 endif

!
! Energy budget
!
 call check( nf90_inq_varid(ncid, 'TLWUP' , varid ), 'TLWUP',missing_para) !TOA LW Up [W/m2]
 if(missing_para) then 
    surface%toa_lw_up = undefined 
 else
    call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%toa_lw_up, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'TLWDN' , varid ), 'TLWDN',missing_para) !TOA LW down [W/m2]
 if(missing_para) then   
    surface%toa_lw_dn = undefined
 else
    call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%toa_lw_dn, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'TSWUP' , varid ), 'TSWUP',missing_para) !TOA SW Up [W/m2]
 if(missing_para) then   
    surface%toa_sw_up = undefined
 else
    call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%toa_sw_up, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'TSWDN' , varid ), 'TSWDN',missing_para) !TOA SW down [W/m2]
 if(missing_para) then
    surface%toa_sw_dn = undefined
 else
    call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%toa_sw_dn, start=(/is,js/), count=(/di,dj/)))
 endif


 call check( nf90_inq_varid(ncid, 'SLWUP' , varid ), 'SLWUP',missing_para) !BOA LW Up [W/m2]
 if(missing_para) then   
    surface%boa_lw_up = undefined
 else
    call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%boa_lw_up, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'SLWDN' , varid ), 'SLWDN',missing_para) !BOA LW down [W/m2]
 if(missing_para) then
    surface%boa_lw_dn = undefined
 else
    call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%boa_lw_dn, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'SSWUP' , varid ), 'SSWUP',missing_para) !BOA SW up [W/m2]
 if(missing_para) then
    surface%boa_sw_up = undefined
 else
    call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%boa_sw_up, start=(/is,js/), count=(/di,dj/)))
 endif

 call check( nf90_inq_varid(ncid, 'SSWDN' , varid ), 'SSWDN',missing_para) !BOA SW down [W/m2]
 if(missing_para) then
    surface%boa_sw_dn = undefined
 else
    call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%boa_sw_dn, start=(/is,js/), count=(/di,dj/)))
 endif


!
! close netCDF file
!
  call check( nf90_close(ncid) )

 if(masterproc) print*,''

 return
 end subroutine rd_CRM_WRF
 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine check(status,para,missing_para)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Check status of reading process of netCDF. 
!
! History:
! 05/2014  Toshi Matsui@NASA GSFC ; included optional statement for skiping
! 05/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!---------------------------------------------------------------------------------------------------
  integer, intent (in) :: status
  character(len=*),intent(in),optional :: para  !optional parameter for missing, but skippable parameters.
  logical,intent(out),optional :: missing_para  !optional parameter for missing, but skippable parameters.

  if( present(missing_para) )  missing_para = .false.

   
  if(status /= nf90_noerr) then
    if(masterproc)  print *, trim(nf90_strerror(status))

    !
    !for optional parameters, keep going on program.
    !
    if( present(para) ) then
          if(masterproc) print*,'Missing ',trim(para)
          if( present(missing_para) ) missing_para = .true.
          return
    endif

    !
    !for required parameter, stop program
    ! 
    call stop_sdsu(myrank, "Stopped missing parameters!"  )
  end if
 end subroutine check


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine check_io(status,bad_io)
  implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Check status of reading process of netCDF. 
!
! History:
! 05/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!---------------------------------------------------------------------------------------------------
  integer, intent (in) :: status
  logical, intent (out) :: bad_io
  bad_io = .false.
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    bad_io = .true.
  end if
 end subroutine check_io

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine rd_CRM_GCESBM
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Read unified GCE output in in netCDF format, adapeted for domain decomposition.  
!
! History:
! 07/2014  Toshi Matsui@NASA GSFC ; Add Morrison MP
! 07/2011  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!
! Tao, W.K., 2003:  Goddard Cumulus Ensemble (GCE) Model: Application for Understanding 
! Precipitation Processes. Meteor. Monogr., 29, 107.
!
!---------------------------------------------------------------------------------------------------

 integer :: i,j,k,l,n    ! looping
 integer :: is,ie,js,je,ks,ke,di,dj,dk
 integer :: ncid, ncidsbm     !netcdf file index 
 integer :: varid    !parameter id
 integer :: ierr     !allocation stat

 real(sdsu_fps) :: dsd_conc_net !drop size ditributions in l,z,y,x order [1/m3]
 real(sdsu_fps),allocatable :: hgt(:), hgt_stag(:) !staggered GCE height level [km]
 real(sdsu_fps),allocatable :: p0(:) !initial pressure [hPa]
 real(sdsu_fps),allocatable :: net2d(:,:)   !2d parameter 
 real(sdsu_fps),allocatable :: net3d(:,:,:) !3d parameter 
 real(sdsu_fps),allocatable :: net4d(:,:,:,:)!4d parameter 
 real(sdsu_fps),allocatable :: rho_air(:,:,:)  !dry air density [kg/m3]
 real(sdsu_fps),allocatable :: u(:,:,:)   
 real(sdsu_fps),allocatable :: v(:,:,:)

 real :: rho(mxlyr)  !dry air density

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1


!
! Memory allocation 
!
 if( .not. allocated(hgt) ) then
 allocate( &
           hgt      (1:               mxlyr  ) ,&
           hgt_stag (0:               mxlyr  ) ,&
           p0       (                   mxlyr) ,&                       
           net2d    (is:ie,js:je             ) ,&          
           net3d    (is:ie,js:je,ks:ke       ) ,&      
           net4d    (is:ie,js:je,ks:ke,1:nbin) ,&  
           rho_air  (is:ie,js:je,ks:ke       ) ,&  
           u        (is:ie,js:je,ks:ke       ) ,&      
           v        (is:ie,js:je,ks:ke       ) ,&        
           stat=ierr )
 if(ierr /= 0) call stop_sdsu(myrank,'MSG rd_CRM_GCESBM: cannot allocate the variable')
 endif

!
! open netCDF file
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name(1:efile_len))//'.cdf'
 if(masterproc) print("(a)"),'MSG rd_CRM_GCESBM: ++++++++++++++++++++ READING NEW FILE ++++++++++++++++++++'
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)

!
! open netCDF file
!
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )


!
! latitude and longitude (new) 
!
 call check( nf90_inq_varid(ncid, 'latitude' , varid ) )  !latitude [deg]
 call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = js, je ; do i = is, ie
    surface(i,j)%lat = net2d(i,j)    ! latitude [deg]
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'longitude' , varid ) )  !longitude [deg]
 call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = js, je ; do i = is, ie
    surface(i,j)%lon = net2d(i,j)    ! longitude [deg]
 enddo ; enddo

!
! staggered height level and all sorts
!
 call check( nf90_inq_varid(ncid, 'z' , varid ) )  ! height in staggered level [m]
 call check( nf90_get_var(ncid, varid, hgt ) )

 call check( nf90_inq_varid(ncid, 'z2' , varid ) )  ! height in staggered level [m]
 call check( nf90_get_var(ncid, varid, hgt_stag ) )


 do k = 0,mxlyr
    hgt_lev(k)            = hgt_stag(k)  * 1.e-3  !height level 1D [km]
    atmos_stag(:,:,k)%hgt = hgt_stag(k)  * 1.e-3  !height level 3D [km]
 enddo


!
! temperature and all sorts
!
 call check( nf90_inq_varid(ncid, 'T' , varid ) )  ! air temperrature [K]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    atmos(i,j,k)%t_air = net3d(i,j,k)  ! air temperature [K]
 enddo ; enddo ; enddo

 do k = ks, ke-1 ; do j = js, je ; do i = is, ie
    atmos_stag(i,j,k)%t_air = 0.5e0 * (atmos(i,j,k)%t_air + atmos(i,j,k+1)%t_air) ![K]
 enddo ; enddo ; enddo

 atmos_stag(:,:,mxlyr)%t_air = atmos(:,:,mxlyr)%t_air  !top level [K]

 call check( nf90_inq_varid(ncid, 'Ts' , varid ) )  !surface skin temperature [K]
 call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = js, je ; do i = is, ie
    surface(i,j)%t_skin = net2d(i,j)    ! skin temp [K]
 enddo ; enddo
 surface%t_soil = surface%t_skin !top-soil temperature [K]

 atmos_stag(:,:,0)%t_air     = surface(:,:)%t_skin !surface air temp = skin temp [K]



!
! pressure and all sorts
!
 call check( nf90_inq_varid(ncid, 'p0' , varid ) )  !initial pressure[hPa]
 call check( nf90_get_var(ncid, varid, p0 ) )

 call check( nf90_inq_varid(ncid, 'psfc' , varid ) )  ! surface pressure[hPa]
 call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = js, je ; do i = is, ie
    atmos_stag(i,j,0)%press = net2d(i,j)  !surface pressure [hPa]
 enddo ; enddo

 do k = ks, ke
    atmos(:,:,k)%press = p0(k) ![hPa] 
 enddo

 do k = ks, ke-1
    atmos_stag(:,:,k)%press = 0.5e0* (p0(k) + p0(k+1)) ![hPa]
 enddo

 atmos_stag(:,:,mxlyr)%press = 0.5e0* (p0(mxlyr) + 0.)     ! top layer

!
! relative humidity
!
 call check( nf90_inq_varid(ncid, 'rh' , varid ) )  !relative humidity [%]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    atmos(i,j,k)%rh = net3d(i,j,k)  !  relative humidity [%]
 enddo ; enddo ; enddo

 do k = ks, ke ; do j = js, je ; do i = is, ie
    rho_air(i,j,k) = atmos(i,j,k)%press / ( atmos(i,j,k)%t_air*2.87e0)  ! dry air density  [kg/m3]
    atmos(i,j,k)%rho_dair = rho_air(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

!
! mixing ratio of condensates
!

 mic_select0: select case(trim(cloud_microphysics))
 case('GOD','GOD10','GMP4ICE','LIN','MORR','MORRH')

 call check( nf90_inq_varid(ncid, 'Qc' , varid ) )  ! cloud water mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_gen(i,j,k)%cloud = net3d(i,j,k)*rho_air(i,j,k) ! cloud water mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qr' , varid ) )     ! rain water mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_gen(i,j,k)%rain = net3d(i,j,k)*rho_air(i,j,k) ! rain water mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qs' , varid ) )  !snow mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_gen(i,j,k)%snow  = net3d(i,j,k)*rho_air(i,j,k)  ! snow mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qg' , varid ) )  !graupel mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    if( trim(cloud_microphysics) == 'MORRH' ) then
      q_gen(i,j,k)%hail = net3d(i,j,k)*rho_air(i,j,k) ! hail mixing ratio [g/m3]
      q_gen(i,j,k)%graupel = 0.
    else
       q_gen(i,j,k)%graupel = net3d(i,j,k)*rho_air(i,j,k) ! graupel mixing ratio [g/m3]
       q_gen(i,j,k)%hail = 0.
    endif
 enddo ; enddo ; enddo

 if( trim(cloud_microphysics) == 'GMP4ICE' ) then
    call check( nf90_inq_varid(ncid, 'Qh' , varid ) )  !graupel mixing ratio [g/kg]
    call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
    do k = ks, ke ; do j = js, je ; do i = is, ie
       q_gen(i,j,k)%hail = net3d(i,j,k)*rho_air(i,j,k) ! hail mixing ratio [g/m3]
    enddo ; enddo ; enddo
 endif

 call check( nf90_inq_varid(ncid, 'Qi' , varid ) )  !ice  mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_gen(i,j,k)%ice =  net3d(i,j,k)*rho_air(i,j,k) ! cloud ice mixing ratio [g/m3]
 enddo ; enddo ; enddo

 if( trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH' ) then

 call check( nf90_inq_varid(ncid, 'Nc' , varid ) )  ! cloud water number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_gen(i,j,k)%cloud = net3d(i,j,k)*rho_air(i,j,k) ! cloud water number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Nr' , varid ) )     ! rain water number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_gen(i,j,k)%rain = net3d(i,j,k)*rho_air(i,j,k) ! rain water number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Ns' , varid ) )  !snow number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_gen(i,j,k)%snow  = net3d(i,j,k)*rho_air(i,j,k)  ! snow number conc [#/m3]
 enddo ; enddo ; enddo


 call check( nf90_inq_varid(ncid, 'Ng' , varid ) )  !graupel number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
 if( trim(cloud_microphysics) == 'MORRH' ) then
    n_gen(i,j,k)%hail = net3d(i,j,k)*rho_air(i,j,k) ! hail number conc [#/m3]
    n_gen(i,j,k)%graupel = 0.
 else
    n_gen(i,j,k)%graupel = net3d(i,j,k)*rho_air(i,j,k) ! graupel number conc [#/m3]
    n_gen(i,j,k)%hail = 0.
 endif
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Ni' , varid ) )  !ice number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_gen(i,j,k)%ice =  net3d(i,j,k)*rho_air(i,j,k) ! cloud ice number conc [#/m3]
 enddo ; enddo ; enddo


 endif

 case('WSM','WDM')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GCE: WSM or WDM schemes do not exist in GCE input')

 case('THOM')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GCE: Thompson scheme does not exist in GCE input')

 case('NTUCLR')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GCE: NTUCLR scheme does not exist in GCE input')


 case('RAMS1','RAMS2') !-------------------------------------------

 call check( nf90_inq_varid(ncid, 'Qc' , varid ) )  ! cloud mode1 (small) mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_rams(i,j,k)%cloud1 = net3d(i,j,k)*rho_air(i,j,k) ! cloud mode1 (small)  mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qc2' , varid ) )  ! cloud mode2 (large) mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_rams(i,j,k)%cloud2 = net3d(i,j,k)*rho_air(i,j,k) ! cloud mode2 (large)  mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qr' , varid ) )     ! rain mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_rams(i,j,k)%rain = net3d(i,j,k)*rho_air(i,j,k) ! rain mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qi' , varid ) )  ! ice mode1 (small) mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_rams(i,j,k)%ice1 = net3d(i,j,k)*rho_air(i,j,k)! ice mode1 (small)  mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qs' , varid ) )   ! ice mode2 (large) mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_rams(i,j,k)%ice2 = net3d(i,j,k)*rho_air(i,j,k) ! ice mode2 (large)  mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qa' , varid ) )  !snow aggregate mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_rams(i,j,k)%snow = net3d(i,j,k)*rho_air(i,j,k)  ! snow aggregate mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qg' , varid ) )  !graupel mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_rams(i,j,k)%graupel = net3d(i,j,k)*rho_air(i,j,k) ! graupel mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qh' , varid ) )  !hail mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    q_rams(i,j,k)%hail = net3d(i,j,k)*rho_air(i,j,k)   ! hail mixing ratio [g/m3]
 enddo ; enddo ; enddo

  n_rams(:,:,:)%cloud1  = 0.e0  !initialize
  n_rams(:,:,:)%cloud2  = 0.e0
  n_rams(:,:,:)%rain    = 0.e0
  n_rams(:,:,:)%ice1    = 0.e0
  n_rams(:,:,:)%ice2    = 0.e0
  n_rams(:,:,:)%snow    = 0.e0
  n_rams(:,:,:)%graupel = 0.e0
  n_rams(:,:,:)%hail    = 0.e0

 if( trim(cloud_microphysics) == 'RAMS2' ) then !RAMS 2-moment only ---------------------------------

 call check( nf90_inq_varid(ncid, 'Nc' , varid ) )  ! cloud mode1 (small) total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_rams(i,j,k)%cloud1 = net3d(i,j,k)*rho_air(i,j,k) ! cloud mode1 (small) total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Nc2' , varid ) )  ! cloud mode2 (large) total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_rams(i,j,k)%cloud2 = net3d(i,j,k)*rho_air(i,j,k) ! cloud mode2 (large) total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Nr' , varid ) )     ! rain total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_rams(i,j,k)%rain = net3d(i,j,k)*rho_air(i,j,k) ! rain total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Ni' , varid ) )  ! ice mode1 (small) total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_rams(i,j,k)%ice1 = net3d(i,j,k)*rho_air(i,j,k)! ice mode1 (small) total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Ns' , varid ) )   ! ice mode2 (large) total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_rams(i,j,k)%ice2 = net3d(i,j,k)*rho_air(i,j,k) ! ice mode2 (large) total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Na' , varid ) )  !snow aggregate total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_rams(i,j,k)%snow = net3d(i,j,k)*rho_air(i,j,k)  ! snow aggregate total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Ng' , varid ) )  !graupel total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_rams(i,j,k)%graupel = net3d(i,j,k)*rho_air(i,j,k) ! graupel total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Nh' , varid ) )  !hail total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    n_rams(i,j,k)%hail = net3d(i,j,k)*rho_air(i,j,k)   ! hail total number conc [#/m3]
 enddo ; enddo ; enddo


 endif !RAMS2


 case('HUCM_SBM') !--------------------------------------------


!
! particle size distribution of condensates
!
 frime_sbm_snow(:,:,:,:)     = 0.e0

 !
 ! liquid  (cloud and rain)
 !
 n_sbm%liq = 0.e0 ; q_sbm%liq = 0.e0
 sdsu_io_file =trim(sdsu_dir_input)//trim(sdsu_inp_name(1:efile_len))//'water.cdf'
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

 call check( nf90_inq_varid(ncidsbm, 'fc' , varid ) )  !cloud-rain drop size distribution [g/g]
 call check( nf90_get_var(ncidsbm, varid, net4d(is:ie,js:je,ks:ke,1:nbin), &
             start=(/is,js,ks,1/), count=(/di,dj,dk,nbin/)))
 do n = 1, nbin ; do k = ks, ke ; do j = js, je ; do i = is, ie
    call sbm_function43(rho_air(i,j,k), x_sbm(n)%liq, drad_sbm(n)%liq*1.e-2, net4d(i,j,k,n), &
                        n_sbm(i,j,k,n)%liq, q_sbm(i,j,k)%liq )
 enddo ; enddo ; enddo ; enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! ice plate
 !
 n_sbm%ice_pla = 0.e0 ; q_sbm%ice_pla = 0.e0
 sdsu_io_file =trim(sdsu_dir_input)//trim(sdsu_inp_name(1:efile_len))//'ice_p.cdf'
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

 call check( nf90_inq_varid(ncidsbm, 'fip' , varid ) )  !ice plate drop size distribution [g/g]
 call check( nf90_get_var(ncidsbm, varid, net4d(is:ie,js:je,ks:ke,1:nbin), &
             start=(/is,js,ks,1/), count=(/di,dj,dk,nbin/)))
 do n = 1, nbin ; do k = ks, ke ; do j = js, je ; do i = is, ie
    call sbm_function43(rho_air(i,j,k), x_sbm(n)%ice_pla, drad_sbm(n)%ice_pla*1.e-2, net4d(i,j,k,n), &
                        n_sbm(i,j,k,n)%ice_pla, q_sbm(i,j,k)%ice_pla )
 enddo ; enddo ; enddo ; enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

 !
 ! ice column
 !
 n_sbm%ice_col = 0.e0 ; q_sbm%ice_col = 0.e0
 sdsu_io_file =trim(sdsu_dir_input)//trim(sdsu_inp_name(1:efile_len))//'ice_c.cdf'
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

 call check( nf90_inq_varid(ncidsbm, 'fic' , varid ) )  !ice column drop size distribution [g/g]
 call check( nf90_get_var(ncidsbm, varid, net4d(is:ie,js:je,ks:ke,1:nbin), &
             start=(/is,js,ks,1/), count=(/di,dj,dk,nbin/)))
 do n = 1, nbin ; do k = ks, ke ; do j = js, je ; do i = is, ie
    call sbm_function43(rho_air(i,j,k), x_sbm(n)%ice_col, drad_sbm(n)%ice_col*1.e-2, net4d(i,j,k,n), &
                        n_sbm(i,j,k,n)%ice_col, q_sbm(i,j,k)%ice_col )
 enddo ; enddo ; enddo ; enddo
 call check( nf90_close(ncidsbm) )  ! close nc file


 !
 ! ice dendrite
 !
 n_sbm%ice_den = 0.e0 ; q_sbm%ice_den = 0.e0
 sdsu_io_file =trim(sdsu_dir_input)//trim(sdsu_inp_name(1:efile_len))//'ice_d.cdf'
#if MPI == 2 
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file
 
 call check( nf90_inq_varid(ncidsbm, 'fid' , varid ) )  !ice dendrite drop size distribution [g/g]
 call check( nf90_get_var(ncidsbm, varid, net4d(is:ie,js:je,ks:ke,1:nbin), &
             start=(/is,js,ks,1/), count=(/di,dj,dk,nbin/)))
 do n = 1, nbin ; do k = ks, ke ; do j = js, je ; do i = is, ie 
    call sbm_function43(rho_air(i,j,k), x_sbm(n)%ice_den, drad_sbm(n)%ice_den*1.e-2, net4d(i,j,k,n), &
                        n_sbm(i,j,k,n)%ice_den, q_sbm(i,j,k)%ice_den )
 enddo ; enddo ; enddo ; enddo 
 call check( nf90_close(ncidsbm) )  ! close nc file


!
! snow aggregate
!

 n_sbm%snow = 0.e0 ; q_sbm%snow = 0.e0
 sdsu_io_file =trim(sdsu_dir_input)//trim(sdsu_inp_name(1:efile_len))//'snow.cdf'
#if MPI == 2 
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

 call check( nf90_inq_varid(ncidsbm, 'fs' , varid ) )  !snow drop size distribution [g/g]
 call check( nf90_get_var(ncidsbm, varid, net4d(is:ie,js:je,ks:ke,1:nbin), &
             start=(/is,js,ks,1/), count=(/di,dj,dk,nbin/)))
 do n = 1, nbin ; do k = ks, ke ; do j = js, je ; do i = is, ie
    call sbm_function43(rho_air(i,j,k), x_sbm(n)%snow, drad_sbm(n)%snow*1.e-2, net4d(i,j,k,n), &
!    call sbm_function43(rho(k)*1.e+3, x_sbm(n)%snow, drad_sbm(n)%snow*1.e-2, net4d(i,j,k,n), &
                        n_sbm(i,j,k,n)%snow, q_sbm(i,j,k)%snow )
 enddo ; enddo ; enddo ; enddo

 call check( nf90_close(ncidsbm) )  ! close nc file


 fmelt_sbm_snow = 0.e0 
 sdsu_io_file =trim(sdsu_dir_input)//trim(sdsu_inp_name(1:efile_len))//'snow_melt.cdf'
#if MPI == 2 
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

 call check( nf90_inq_varid(ncidsbm, 'fls' , varid ) )  !melting fraction  [-]
 call check( nf90_get_var(ncidsbm, varid, net4d(is:ie,js:je,ks:ke,1:nbin), &
             start=(/is,js,ks,1/), count=(/di,dj,dk,nbin/)))
 do n = 1, nbin ; do k = ks, ke ; do j = js, je ; do i = is, ie
    fmelt_sbm_snow(i,j,k,n) = net4d(i,j,k,n)
 enddo ; enddo ; enddo ; enddo
 call check( nf90_close(ncidsbm) )  ! close nc file


!
! graupel
!
 n_sbm%graupel = 0.e0 ; q_sbm%graupel = 0.e0
 sdsu_io_file =trim(sdsu_dir_input)//trim(sdsu_inp_name(1:efile_len))//'graupel.cdf'
#if MPI == 2 
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

 call check( nf90_inq_varid(ncidsbm, 'fg' , varid ) )  !graupel drop size distribution [g/g]
 call check( nf90_get_var(ncidsbm, varid, net4d(is:ie,js:je,ks:ke,1:nbin), &
             start=(/is,js,ks,1/), count=(/di,dj,dk,nbin/)))
 do n = 1, nbin ; do k = ks, ke ; do j = js, je ; do i = is, ie
    call sbm_function43(rho_air(i,j,k), x_sbm(n)%graupel, drad_sbm(n)%graupel*1.e-2, net4d(i,j,k,n), &
                        n_sbm(i,j,k,n)%graupel, q_sbm(i,j,k)%graupel )
 enddo ; enddo ; enddo ; enddo
 call check( nf90_close(ncidsbm) )  ! close nc file

!
! hail
!
 n_sbm%hail = 0.e0 ; q_sbm%hail = 0.e0
 sdsu_io_file =trim(sdsu_dir_input)//trim(sdsu_inp_name(1:efile_len))//'hail.cdf'
#if MPI == 2 
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

 call check( nf90_inq_varid(ncidsbm, 'fh' , varid ) )  !hail drop size distribution [g/g]
 call check( nf90_get_var(ncidsbm, varid, net4d(is:ie,js:je,ks:ke,1:nbin), &
             start=(/is,js,ks,1/), count=(/di,dj,dk,nbin/)))
 do n = 1, nbin ; do k = ks, ke ; do j = js, je ; do i = is, ie
    call sbm_function43(rho_air(i,j,k), x_sbm(n)%hail, drad_sbm(n)%hail*1.e-2, net4d(i,j,k,n), &
                        n_sbm(i,j,k,n)%hail, q_sbm(i,j,k)%hail )
 enddo ; enddo ; enddo ; enddo
 call check( nf90_close(ncidsbm) )  ! close nc file


 case('HUCM_SBM43')
 frime_sbm_snow(:,:,:,:)     = 0.e0
 fmelt_sbm_snow(:,:,:,:)     = 0.e0

   call stop_sdsu(myrank,'MSG: rd_CRM_GCESBM: new option HUCM_SBM43 under construction')

 case default
    ! do nothing
 end select mic_select0

!
! vertical velocity
!
 call check( nf90_inq_varid(ncid, 'w' , varid ) )  ! w-component wind [m/s]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do k = ks, ke ; do j = js, je ; do i = is, ie
    atmos(i,j,k)%w = net3d(i,j,k)  ! vertical velocity [m/s]
 enddo ; enddo ; enddo

 !stagger level
 do k = ks, ke-1
    atmos_stag(:,:,k)%w = 0.5e0*(atmos(:,:,k)%w + atmos(:,:,k+1)%w)
 enddo
 atmos_stag(:,:, 0)%w = 0. ! bottom-level vertical velocity [m/s]
 atmos_stag(:,:,ke)%w = 0. ! top-level vertical velocity    [m/s]





!
! surface wind speed
!
 call check( nf90_inq_varid(ncid, 'u' , varid ) )  ! u-component wind [m/s]
 call check( nf90_get_var(ncid, varid, u(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 call check( nf90_inq_varid(ncid, 'v' , varid ) )  ! v-component wind [m/s]
 call check( nf90_get_var(ncid, varid, v(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do j = js, je ; do i = is, ie
   surface(i,j)%u10m = sqrt( u(i,j,1)**2 + v(i,j,1)**2  )  ! surface wind [m/s] 
 enddo ; enddo

!
! close netCDF file
!
 call check( nf90_close(ncid) )

 return
 end subroutine rd_CRM_GCESBM

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine rd_CRM_GCE_unified
  implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Read unified GCE output in in netCDF format. 
!
! History:
! 05/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!
! Tao, W.K., 2003:  Goddard Cumulus Ensemble (GCE) Model: Application for Understanding 
! Precipitation Processes. Meteor. Monogr., 29, 107.
!
! Zeng, X., W.K. Tao, S. Lang, A. Y. Hou, M. Zhang, and J. Simpson, 2008:  On the sensitivity 
! of atmospheric ensemble states to cloud microphysics in long-term cloud-resolving model 
! simulations, Journal of the Meteorological Society of Japan, (In press).
!
!---------------------------------------------------------------------------------------------------

 integer :: i,j,k,l    !looping
 integer :: is,ie,js,je,ks,ke,di,dj,dk
 integer :: ncid     !# of variables to inquire  
 integer :: varid    !# of variable id
 integer :: ierr     !allocation stat

 real(sdsu_fps) :: dsd_conc_net !drop size ditributions in l,z,y,x order [1/m3]
 real(sdsu_fps),allocatable :: hgt2_sbm(:) !staggered GCE height level (SBM only remove lator) toshii [km]
 real(sdsu_fps),allocatable :: hgt2(:) !staggered GCE height level [km]
 real(sdsu_fps),allocatable :: temp_net(:,:,:) !layer averaged temperaure in z,y,x order [K]
 real(sdsu_fps),allocatable :: p0(:) !initial pressure [hPa]
 real(sdsu_fps),allocatable :: net2d(:,:)   !2d parameter in y,x order
 real(sdsu_fps),allocatable :: net3d(:,:,:) !3d parameter in z,y,x order
 real(sdsu_fps),allocatable :: net4d(:,:,:,:)!4d parameter in l,z,y,x order
 real(sdsu_fps),allocatable :: rho_air(:,:,:)  !moist air density [kg/m3]
 real(sdsu_fps),allocatable :: u(:,:,:)   
 real(sdsu_fps),allocatable :: v(:,:,:)


!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1

!
! Memory allocation 
!
 if( .not. allocated(hgt2) ) then
 allocate( &
 hgt2(mxlyr) ,&                       ! staggered GCE height level [m]
 hgt2_sbm(mxlyr+2) ,&!staggered GCE height level (SBM only remove lator) toshii [km]
 temp_net(mxlyr,mxgridy,mxgridx) ,&   ! layer averaged temperaure in z,y,x order [C]
 p0(mxlyr) ,&                         ! initial pressure [hPa] 
 net2d(mxgridy,mxgridx)   ,&          ! 2d parameer in y,x order
 net3d(mxlyr,mxgridy,mxgridx) ,&      ! 3d parameters in z,y,x order
 net4d(nbin,mxlyr,mxgridy,mxgridx) ,& ! 4d parameters in z,y,x order
 rho_air(mxlyr,mxgridy,mxgridx) ,&    ! density of dry air [kg/m3]
 u(mxlyr,mxgridy,mxgridx) ,&          ! u wind
 v(mxlyr,mxgridy,mxgridx) ,&          ! v wind
 stat=ierr )
 if(ierr /= 0) call stop_sdsu(myrank,'MSG rd_CRM_GCE_unified: cannot allocate the variable')
 endif

!
! open netCDF file
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name)
 if(masterproc) print("(a)"),'MSG rd_CRM_GCE_unified: ++++++++++++++++++++ READING NEW FILE ++++++++++++++++++++'
 if(masterproc) print("(a)"),'Input ->',trim(sdsu_io_file)

 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! staggered height level and all sorts
!
 !for SBM option: make sure that z2 must be consistent to unified GCE output. toshii

 call check( nf90_inq_varid(ncid, 'z2' , varid ) )  ! height in staggered level [m]
 call check( nf90_get_var(ncid, varid, hgt2 ) )

 do j = js, je ; do i = is, ie
 do k = 0,mxlyr-1 
    hgt_lev(k)       = hgt2(k+1)/1.e3  !height 1D [km]
    atmos_stag(i,j,k)%hgt = hgt2(k+1)/1.e3  !height 3D [km]
 enddo
 enddo ; enddo

 do j = js, je ; do i = is, ie

 hgt_lev(mxlyr)       = hgt_lev(mxlyr-1)       + hgt_lev(mxlyr-1)       - hgt_lev(mxlyr-2)
 atmos_stag(i,j,mxlyr)%hgt = atmos_stag(i,j,mxlyr-1)%hgt + atmos_stag(i,j,mxlyr-1)%hgt &
                             - atmos_stag(i,j,mxlyr-2)%hgt
 enddo ; enddo

!
! temperature and all sorts
!
 call check( nf90_inq_varid(ncid, 'T' , varid ) )  ! air temperrature [Celvius]
 call check( nf90_get_var(ncid, varid, temp_net ) )

 do i = is, ie ; do j = js, je ; do k = ks, ke
    atmos(i,j,k)%t_air = temp_net(k,j,i) + const_Kel2Cel  ![K]
 enddo ; enddo ; enddo

 do i = is, ie ; do j = js, je ; do k = ks, ke
    atmos_stag(i,j,k)%t_air = 0.5e0 * (atmos(i,j,k)%t_air + atmos(i,j,k+1)%t_air) ![K]
 enddo ; enddo ; enddo

 do j = js, je ; do i = is, ie 
 atmos_stag(i,j,0)%t_air     = atmos(i,j,1)%t_air      !surface level [K]
 atmos_stag(i,j,mxlyr)%t_air = atmos(i,j,mxlyr)%t_air  !top level [K]
 enddo  ; enddo

 call check( nf90_inq_varid(ncid, 'Ts' , varid ) )  !surface skin temperature [Celcius]
 call check( nf90_get_var(ncid, varid, net2d ) )

 do j = js, je ; do i = is, ie
    surface(i,j)%t_skin = net2d(j,i) + const_Kel2Cel    ! skin temp [K]
    surface(i,j)%t_soil = surface(i,j)%t_skin !top-soil temperature [K]
 enddo ; enddo

!
! pressure and all sorts
!
 call check( nf90_inq_varid(ncid, 'p0' , varid ) )  !initial pressure[hPa]
 call check( nf90_get_var(ncid, varid, p0 ) )

 do j = js, je ; do i = is, ie
 do k = 1, mxlyr
    atmos(i,j,k)%press = p0(k) ![hPa] 
 enddo
 enddo ; enddo

 do j = js, je ; do i = is, ie
 do k = 1, mxlyr-1
    atmos_stag(i,j,k)%press = 0.5e0* (p0(k) + p0(k+1)) ![hPa]
 enddo
 enddo ; enddo

 do j = js, je ; do i = is, ie
    atmos_stag(i,j,0)%press = p0(1) +  0.5e0*(p0(1) - p0(2))  ! surface  
    atmos_stag(i,j,mxlyr)%press = 0.5e0* (p0(mxlyr) + 0.)     ! top layer
 enddo ; enddo

!
! relative humidity
!
 call check( nf90_inq_varid(ncid, 'f' , varid ) )  !relative humidity [%]
 call check( nf90_get_var(ncid, varid, net3d ) )

 do i = is, ie ; do j = js, je ; do k = ks, ke
    atmos(i,j,k)%rh = net3d(k,j,i)     ! relative humidity [%]
 enddo ; enddo ; enddo

!
! mixing ratio of condensates
!

 mic_select0: select case(trim(cloud_microphysics))
 case('GOD','GOD10','LIN','GMP4ICE')

 call check( nf90_inq_varid(ncid, 'Qv' , varid ) )  ! water vapor mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )

 do i = is, ie ; do j = js, je ; do k = ks, ke
    rho_air(k,j,i) = p0(k) / ((temp_net(k,j,i)+const_Kel2Cel)*2.87e0)  !dry air density [kg/m3] 
    atmos(i,j,k)%rho_dair = rho_air(k,j,i) * 1.e3  ! [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qc' , varid ) )  ! cloud water mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_gen(i,j,k)%cloud = net3d(k,j,i)*rho_air(k,j,i) ! cloud water mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qr' , varid ) )     ! rain water mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_gen(i,j,k)%rain = net3d(k,j,i)*rho_air(k,j,i) ! rain water mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qs' , varid ) )  !snow mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_gen(i,j,k)%snow  = net3d(k,j,i)*rho_air(k,j,i)  ! snow mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qg' , varid ) )  !graupel mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_gen(i,j,k)%graupel = net3d(k,j,i)*rho_air(k,j,i) ! graupel mixing ratio [g/m3]
 enddo ; enddo ; enddo

 if( trim(cloud_microphysics) == 'GMP4ICE' ) then
    call check( nf90_inq_varid(ncid, 'Qh' , varid ) )  !hail mixing ratio [g/kg]
    call check( nf90_get_var(ncid, varid, net3d ) )
    do i = is, ie ; do j = js, je ; do k = ks, ke
       q_gen(i,j,k)%hail = net3d(k,j,i)*rho_air(k,j,i) ! hail mixing ratio [g/m3]
    enddo ; enddo ; enddo
 else
    q_gen%hail = 0.e0 
 endif
 
 call check( nf90_inq_varid(ncid, 'Qi' , varid ) )  !ice  mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_gen(i,j,k)%ice = net3d(k,j,i)*rho_air(k,j,i)! cloud ice mixing ratio [g/m3]
 enddo ; enddo ; enddo

 case('WSM','WDM','MORR','MORRH','THOM','NTUCLR','NSSL2MMT','P3')

   call stop_sdsu(myrank,'MSG rd_CRM_GCE_unified: Scheme does not exist in GCE')

 case('RAMS1','RAMS2') !-------------------------------------------

 call check( nf90_inq_varid(ncid, 'Qv' , varid ) )  ! water vapor mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )

 do i = is, ie ; do j = js, je ; do k = ks, ke
    rho_air(k,j,i) = p0(k) / ((temp_net(k,j,i)+const_Kel2Cel)*2.87e0)  !dry air density [kg/m3] 
    atmos(i,j,k)%rho_dair = rho_air(k,j,i) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qc' , varid ) )  ! cloud mode1 (small) mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_rams(i,j,k)%cloud1 = net3d(k,j,i)*rho_air(k,j,i) ! cloud mode1 (small)  mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qc2' , varid ) )  ! cloud mode2 (large) mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_rams(i,j,k)%cloud2 = net3d(k,j,i)*rho_air(k,j,i) ! cloud mode2 (large)  mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qr' , varid ) )     ! rain mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_rams(i,j,k)%rain = net3d(k,j,i)*rho_air(k,j,i) ! rain mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qi' , varid ) )  ! ice mode1 (small) mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_rams(i,j,k)%ice1 = net3d(k,j,i)*rho_air(k,j,i)! ice mode1 (small)  mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qs' , varid ) )   ! ice mode2 (large) mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )         ! Note that Qs is deinfed as snow in RAMS.
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_rams(i,j,k)%ice2 = net3d(k,j,i)*rho_air(k,j,i) ! ice mode2 (large)  mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qa' , varid ) )  !snow aggregate mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_rams(i,j,k)%snow = net3d(k,j,i)*rho_air(k,j,i)  ! snow aggregate mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qg' , varid ) )  !graupel mixing ratio [g/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_rams(i,j,k)%graupel = net3d(k,j,i)*rho_air(k,j,i) ! graupel mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Qh' , varid ) )  !hail mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_rams(i,j,k)%hail = net3d(k,j,i)*rho_air(k,j,i)   ! hail mixing ratio [g/m3]
 enddo ; enddo ; enddo

  n_rams(:,:,:)%cloud1  = 0.e0  !initialize
  n_rams(:,:,:)%cloud2  = 0.e0
  n_rams(:,:,:)%rain    = 0.e0
  n_rams(:,:,:)%ice1    = 0.e0
  n_rams(:,:,:)%ice2    = 0.e0
  n_rams(:,:,:)%snow    = 0.e0
  n_rams(:,:,:)%graupel = 0.e0
  n_rams(:,:,:)%hail    = 0.e0

 if( trim(cloud_microphysics) == 'RAMS2' ) then !RAMS 2-moment only ---------------------------------

 call check( nf90_inq_varid(ncid, 'Nc' , varid ) )  ! cloud mode1 (small) total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    n_rams(i,j,k)%cloud1 = net3d(k,j,i)*rho_air(k,j,i) ! cloud mode1 (small) total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Nc2' , varid ) )  ! cloud mode2 (large) total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    n_rams(i,j,k)%cloud2 = net3d(k,j,i)*rho_air(k,j,i) ! cloud mode2 (large) total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Nr' , varid ) )     ! rain total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    n_rams(i,j,k)%rain = net3d(k,j,i)*rho_air(k,j,i) ! rain total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Ni' , varid ) )  ! ice mode1 (small) total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    n_rams(i,j,k)%ice1 = net3d(k,j,i)*rho_air(k,j,i)! ice mode1 (small) total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Ns' , varid ) )   ! ice mode2 (large) total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )         ! Note that Qs is deinfed as snow in RAMS.
 do i = is, ie ; do j = js, je ; do k = ks, ke
    n_rams(i,j,k)%ice2 = net3d(k,j,i)*rho_air(k,j,i) ! ice mode2 (large) total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Na' , varid ) )  !snow aggregate total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    n_rams(i,j,k)%snow = net3d(k,j,i)*rho_air(k,j,i)  ! snow aggregate total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Ng' , varid ) )  !graupel total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    n_rams(i,j,k)%graupel = net3d(k,j,i)*rho_air(k,j,i) ! graupel total number conc [#/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'Nh' , varid ) )  !hail total number conc [#/kg]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    n_rams(i,j,k)%hail = net3d(k,j,i)*rho_air(k,j,i)   ! hail total number conc [#/m3]
 enddo ; enddo ; enddo


 endif !RAMS2


 case('HUCM_SBM') !--------------------------------------------


!
! height level and all sorts (temporal using old 2D GCE SBM result from Xiaowen), remove lator toshii
!
 call check( nf90_inq_varid(ncid, 'z2' , varid ) )  ! height in staggered level
 call check( nf90_get_var(ncid, varid, hgt2 ) )

 do i = is, ie ; do j = js, je 
 do k = 0,mxlyr
    hgt_lev(k)       = hgt2(k+2)  !height 1D [km]
    atmos_stag(i,j,k)%hgt = hgt2(k+2)  !height 3D [km]
 enddo
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'qcl' , varid ) )  !cloud water mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_sbm(i,j,k)%liq = net3d(k,j,i)     ! cloud water mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'qrn' , varid ) )  !rain water mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_sbm(i,j,k)%liq = q_sbm(i,j,k)%liq + net3d(k,j,i)     ! rain water mixing ratio [g/m3] (cloud+rain=liq)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'qcic' , varid ) )  !column ice  mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_sbm(i,j,k)%ice_col = net3d(k,j,i)              ! cloud ice mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'qcip' , varid ) )  !plate ice  mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_sbm(i,j,k)%ice_pla = net3d(k,j,i)     ! cloud ice mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'qcid' , varid ) )  !dendrite ice  mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_sbm(i,j,k)%ice_den = net3d(k,j,i)     ! cloud ice mixing ratio [g/m3]
 enddo ; enddo ; enddo


 call check( nf90_inq_varid(ncid, 'qcs' , varid ) )  !snow mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_sbm(i,j,k)%snow = net3d(k,j,i)     ! snow mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'qcg' , varid ) )  !graupel mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_sbm(i,j,k)%graupel = net3d(k,j,i)     ! graupel mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'qch' , varid ) )  !hail mixing ratio [g/m3]
 call check( nf90_get_var(ncid, varid, net3d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke
    q_sbm(i,j,k)%hail = net3d(k,j,i)     ! hail mixing ratio [g/m3]
 enddo ; enddo ; enddo

!
! size distribution of condensates
!
 call check( nf90_inq_varid(ncid, 'fc' , varid ) )  !cloud-rain drop size distribution [g/cm3]
 call check( nf90_get_var(ncid, varid, net4d ) )
  do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    dsd_conc_net = 3.e0*net4d(l,k,j,i)*1.e8/x_sbm(l)%liq/rad_sbm(l)%liq  ! [1/m4] <- [g/m3]  
    if( dsd_conc_net <= 1.e0 ) dsd_conc_net = 0.e0
    n_sbm(i,j,k,l)%liq       = dsd_conc_net    ! liquid (cloud-rain) particle size density [1/m4]
 enddo ; enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'fic' , varid ) )  !column shape ice drop size distribution [g/cm3]
 call check( nf90_get_var(ncid, varid, net4d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    dsd_conc_net = 3.e0*net4d(l,k,j,i)*1.e8/x_sbm(l)%ice_col/rad_sbm(l)%ice_col  ! [1/m4] <- [g/m3]
    if( dsd_conc_net <= 1.e0) dsd_conc_net = 0.e0
    n_sbm(i,j,k,l)%ice_col   = dsd_conc_net     ! cloud ice column particle size density [1/m4]
 enddo ; enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'fip' , varid ) )  !plate shape ice drop size distribution [g/cm3]
 call check( nf90_get_var(ncid, varid, net4d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    dsd_conc_net = 3.e0*net4d(l,k,j,i)*1.e8/x_sbm(l)%ice_pla/rad_sbm(l)%ice_pla  ! [1/m4] <- [g/m3]
    if( dsd_conc_net <= 1.e0) dsd_conc_net = 0.e0
    n_sbm(i,j,k,l)%ice_pla    =   dsd_conc_net ! cloud ice plate particle size density [1/m4]
 enddo ; enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'fid' , varid ) )  !dendrite shape ice drop size distribution [g/cm3]
 call check( nf90_get_var(ncid, varid, net4d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    dsd_conc_net = 3.e0*net4d(l,k,j,i)*1.e8/x_sbm(l)%ice_den/rad_sbm(l)%ice_den  ! [1/m4] <- [g/m3]
    if( dsd_conc_net <= 1.e0) dsd_conc_net = 0.e0
    n_sbm(i,j,k,l)%ice_den = dsd_conc_net ! cloud ice dendride particle size density [1/m4]
 enddo ; enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'fs' , varid ) )  !snow drop size distribution [g/cm3]
 call check( nf90_get_var(ncid, varid, net4d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    dsd_conc_net = 3.e0*net4d(l,k,j,i)*1.e8/x_sbm(l)%snow/rad_sbm(l)%snow ! [1/m4] <- [g/m3]
    if( dsd_conc_net <= 1.e0) dsd_conc_net = 0.e0
    n_sbm(i,j,k,l)%snow      = dsd_conc_net     ! snow particle size density [1/m4]
 enddo ; enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'fg' , varid ) )  !grauple drop size distribution [g/cm3]
 call check( nf90_get_var(ncid, varid, net4d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    dsd_conc_net = 3.e0*net4d(l,k,j,i)*1.e8/x_sbm(l)%graupel/rad_sbm(l)%graupel  ! [1/m4] <- [g/m3]
    if( dsd_conc_net <= 1.e0) dsd_conc_net = 0.e0
    n_sbm(i,j,k,l)%graupel   = dsd_conc_net     ! graupel particle size density [1/m4]
 enddo ; enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'fh' , varid ) )  !hail drop size distribution [g/cm3]
 call check( nf90_get_var(ncid, varid, net4d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    dsd_conc_net = 3.e0*net4d(l,k,j,i)*1.e8/x_sbm(l)%hail/rad_sbm(l)%hail  ! [1/m4] <- [g/m3]
    if( dsd_conc_net <= 1.e0) dsd_conc_net = 0.e0
    n_sbm(i,j,k,l)%hail       = dsd_conc_net     ! hail particle size density [1/m4]
 enddo ; enddo ; enddo ; enddo

!
! melting portion of hydrometeor
!

 call check( nf90_inq_varid(ncid, 'fls' , varid ) )  !snow melting portion  [-]
 call check( nf90_get_var(ncid, varid, net4d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    fmelt_sbm_snow(i,j,k,l)       = net4d(l,k,j,i)  ! snow melting fraction  [-]
 enddo ; enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'flg' , varid ) )   ! graupel melting portion  [-]
 call check( nf90_get_var(ncid, varid, net4d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    fmelt_sbm_graupel(i,j,k,l)   = net4d(l,k,j,i)     ! graupel melting portion  [-]
 enddo ; enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'flh' , varid ) )  !hail melting portion  [-]
 call check( nf90_get_var(ncid, varid, net4d ) )
 do i = is, ie ; do j = js, je ; do k = ks, ke ; do l = 1, nbin
    fmelt_sbm_hail(i,j,k,l)       = net4d(l,k,j,i)     ! hail melting portion  [-]
 enddo ; enddo ; enddo ; enddo

 case('HUCM_SBM43')
   frime_sbm_snow(:,:,:,:)     = 0.e0
   fmelt_sbm_snow(:,:,:,:)     = 0.e0

   call stop_sdsu(myrank,'MSG: rd_CRM_GCE_unified: new option HUCM_SBM43 under construction')

 case default
    ! do nothing
 end select mic_select0
 
!
! winds
!
 call check( nf90_inq_varid(ncid, 'w' , varid ) )  ! w-component wind [m/s]
 call check( nf90_get_var(ncid, varid, net3d ) )

 do i = is, ie ; do j = js, je ; do k = ks, ke 
    atmos_stag(i,j,k)%w = net3d(k,j,i)  ! vertical velocity [m/s]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'u' , varid ) )  ! u-component wind [m/s]
 call check( nf90_get_var(ncid, varid, u ) )
 call check( nf90_inq_varid(ncid, 'v' , varid ) )  ! v-component wind [m/s]
 call check( nf90_get_var(ncid, varid, v ) )

 do i = is, ie ; do j = js, je 
   surface(i,j)%u10m = sqrt( u(1,j,i)**2 + v(1,j,i)**2  )  ! surface wind [m/2] 
 enddo ; enddo

!
! close netCDF file
!
  call check( nf90_close(ncid) )

 return
 end subroutine rd_CRM_GCE_unified

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine sbm_function(x,dr,f,n,q)
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Convert HUCM SBM (33bin version) PSD function to common mass mixing ration (or number 
!  concentrations) per volume. 
!
! History:
!  12/2010  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: x      ! mass per particle [g]
 real(sdsu_fps),intent(in) :: dr     ! SBM size-bin width [m]
 real(sdsu_fps),intent(in) :: f      ! SBM PSD function [/ (cm3 g)]
 real(sdsu_fps),intent(out) :: n     ! particle size density [1/m4]
 real(sdsu_fps),intent(inout) :: q   ! total mixing ratio [g/m3] 
 real(sdsu_fps) :: q_bin  ! mixing ratio per one bin [g/m3]
 real(sdsu_fps) :: n_bin  ! drop size ditributions [1/m3]

!
! initialize
!
 n_bin = 0.e0 ; q_bin = 0.e0

 q_bin = f*3.0e0*0.23105e0*x*x*1.0e6 ! mass conc per bin [g/m3]

 n_bin = q_bin / x / dr              ! # conc per bin [#/m4]

!
! clean numerical noise
!
 if( n_bin <= 0.e0 ) then
     n_bin = 0.e0  ; q_bin = 0.e0
 endif

!
! update PSD and total mixing ratio
!
 n = n_bin     ! particle size density [1/m4]
 q = q + q_bin ! total mixing ratio [g/m3]


 return
 end subroutine sbm_function

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine sbm_function43(rho,x,dr,m,n,q)
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Convert HUCM SBM (43bin version) PSD mass mixing ratio to common mass mixing (or number 
!  concentrations) per volume. 
!
! History:
!  12/2010  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: rho    ! dry air density [kg/m3]
 real(sdsu_fps),intent(in) :: x      ! mass per particle [g]
 real(sdsu_fps),intent(in) :: dr     ! SBM size-bin width [m]
 real(sdsu_fps),intent(in) :: m      ! SBM PSD mass mixing ratio [kg/kg]
 real(sdsu_fps),intent(out) :: n     ! particle size density [1/m4]
 real(sdsu_fps),intent(inout) :: q   ! total mixing ratio [g/m3] 
 real(sdsu_fps) :: q_bin  ! mixing ratio per one bin [g/m3]
 real(sdsu_fps) :: n_bin  ! drop size ditributions [1/m3]

!
! initialize
!
 n_bin = 0.e0 ; q_bin = 0.e0

 q_bin = m * rho * 1.e3      ! mass conc per bin [g/m3]

 n_bin = q_bin / x / dr      ! # conc per bin [#/m4]

 if( isnan(m) ) then
     print*,'MSG sbm_function43, find NaN in m' 
 endif

!
! clean numerical noise
!
 if( n_bin <= 0.e0 ) then
     n_bin = 0.e0  ; q_bin = 0.e0
 endif

!
! update PSD and total mixing ratio
!
 n = n_bin     ! particle size density [1/m4]
 q = q + q_bin ! total mixing ratio [g/m3]


 return
 end subroutine sbm_function43

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine rd_CRM_MMF_GPROF
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
! Read MMF output in netCDF format, specially designed for GPROF database. 
! It first reads as MMF structure (4D), and transfer them into usual CRM strucure (3D) 
! 
! History:
! 04/2012  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------

 integer :: xg,yg,j_crm,i,j,k,l  !looping
 integer :: ncid                !# of variables to inquire  
 integer :: varid               !# of variable id
 integer :: ierr                !allocation stat
 integer :: dimid,tdim
 real(sdsu_fps),allocatable :: hgt2(:) !staggered GCE height level [km]
 real(sdsu_fps),allocatable :: p0(:) !initial pressure [hPa]
 real(sdsu_fps),allocatable :: rho(:,:,:)  !moist air density [kg/m3]

 real(sdsu_fps),allocatable :: var_5d(:,:,:,:,:), var_4d(:,:,:,:),var_3d(:,:,:),var_2d(:,:),var_1d(:)
 real(sdsu_fps) :: lon_geos(mxgridx_geos),lat_geos(mxgridy_geos)
 logical,parameter :: temporal = .false.

!new
 integer :: is,ie,js,je,ks,ke,di,dj,dk
 integer :: se  !end sample index
 integer, dimension(nf90_max_var_dims) :: dimIDs

 real(sdsu_fps),parameter :: biass_correction = 0.5e0  !reduce scattering

 real(sdsu_fps) :: q_orig, q_tune !original and tuned mass [g/m3]
 real(sdsu_fps) :: wgt1, wgt2

!
!  open netCDF file
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name)
 if(masterproc) print("(a)"),'MSG rd_CRM_MMF: ++++++++++++++++++++ READING NEW FILE ++++++++++++++++++++'
 print("(a)"),'Input ->',trim(sdsu_io_file)

 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! read file-dependent sample number
!
 call check( nf90_inq_varid (ncid, 'lon' , varid ) )                       ! find varid of 'lon'
 call check( nf90_inquire_variable ( ncid, varid, dimids = dimIDs ) )      ! get dimIDs of 'lon'
 call check( nf90_inquire_dimension( ncid, dimIDs(1), len = max_sample_gprof ) ) ! get maximum sample number of 'lon'


!
! For Testing PURPOSE !!!!!
!
! max_sample_gprof = 5 !for testing/development (toshii)...

 if( max_sample_gprof > mxgridy ) max_sample_gprof = mxgridy  
 ! correct max_sample number in case of maximum domain bounds
 print*,'MSG rd_CRM_MMF_GPROF: max_sample_gprof = ',max_sample_gprof  !global parameters

 se = max_sample_gprof  !simplify index

 is=myi_start ; ie=myi_end ; js=1 ; je=max_sample_gprof ; ks=1 ; ke=mxlyr
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1


!----------------allocation-------------------------------
 allocate( var_1d(js:je            )  ,& !
           var_2d(js:je,      ks:ke)  ,& !
           var_3d(js:je,is:ie,ks:ke)  ,& !
           stat=ierr )
 if(ierr /= 0) call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: cannot allocate the variable')

!
! longitude
!
 call check( nf90_inq_varid(ncid, 'lon' , varid ) )
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
   surface(is:ie,j)%lon = var_1d(j)  ! longitude [deg]  toshiii change this for more heterogeneous emissivty 
 enddo

!
! latitude
!
 call check( nf90_inq_varid(ncid, 'lat' , varid ) )
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
   surface(is:ie,j)%lat = var_1d(j)  ! latidue [deg]
 enddo


!
! surface height (terrain height)
!
 call check( nf90_inq_varid(ncid, 'ht_sfc' , varid ) )  ! DEM - surface height [m]
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
   atmos_stag(is:ie,j,0)%hgt = var_1d(j)*1e-3  !surface height [km]
   !print*,'hgt',j,atmos_stag(is,j,0)%hgt
 enddo

!
! CRM height (atmospheric height), which do not accounted DEM.
!
 call check( nf90_inq_varid(ncid, 'crm_ht' , varid ) ) ! CRM height (middle) above sea level [m]
 call check( nf90_get_var(ncid, varid, var_2d(js:je,ks:ke), start=(/js,ks/), count=(/dj,dk/) ) )
 do j = js, je
    do k = ks,ke-1
          atmos_stag(is:ie,j,k)%hgt = 0.5 * ( var_2d(j,k) + var_2d(j,k+1) ) * 1e-3  &
                                     + atmos_stag(is:ie,j,0)%hgt   ! [km]
    enddo
!    print*,'stag hgt', atmos_stag(is,j,ks)%hgt
    atmos_stag(is:ie,j,ke)%hgt = atmos_stag(is:ie,j,ke-1)%hgt + &
                                 (var_2d(j,ke) - var_2d(j,ke-1))*1e-3 + &
                                 atmos_stag(is:ie,j,0)%hgt !top level [km]
!    print*,'stag hgt',j, atmos_stag(is,j,0)%hgt, atmos_stag(is,j,38)%hgt
 enddo

!stop

!
! middle pressure
!
 call check( nf90_inq_varid(ncid, 'crm_pres' , varid ) ) ! pressure of CRM [Pa]
 call check( nf90_get_var(ncid, varid, var_2d(js:je,ks:ke), start=(/js,ks/), count=(/dj,dk/) ) )
 do j = js, je
    do k = ks, ke
         atmos(is:ie,j,k)%press = var_2d(j,k)*1e-2  ! [hPa] <- [Pa]  
    enddo
 enddo

!
! interface pressure
! 
 do k = ks, ke-1
   atmos_stag(:,:,k)%press = 0.5e0*(atmos(:,:,k)%press + atmos(:,:,k+1)%press) !interafce [hPa]
 enddo
 atmos_stag(:,:,0)%press = atmos(:,:,1)%press + atmos(:,:,1)%press - atmos_stag(:,:,1)%press   !surface level [hPa]
 atmos_stag(:,:,mxlyr)%press = 0.5e0*(atmos(:,:,mxlyr)%press + 0.) !top pressure [hPa]

!
! skin temperature
! 
 call check( nf90_inq_varid(ncid, 'tskin' , varid ) )  ! surface skin temperature [K]
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
   surface(is:ie,j)%t_skin = var_1d(j)  !surface  skin temperature [K]
 enddo

!
! air temperature (middle level)
!
 call check( nf90_inq_varid(ncid, 'crm_ta' , varid ) ) ! temperature of [K]
 call check( nf90_get_var(ncid, varid, var_3d(js:je,is:ie,ks:ke), start=(/js,is,ks/), count=(/dj,di,dk/) ) )
 do k = ks,ke ; do j = js,je ; do i = is,ie
    atmos(i,j,k)%t_air = var_3d(j,i,k)   ! temperature  [K]
 enddo ; enddo ; enddo

!
! air temperature (interface)
!
 do k = 1, mxlyr-1
    atmos_stag(:,:,k)%t_air = 0.5e0 * (atmos(:,:,k)%t_air + atmos(:,:,k+1)%t_air) ![K]
 enddo
 atmos_stag(:,:,0)%t_air     = 0.5e0* (surface(:,:)%t_skin+atmos_stag(:,:,1)%t_air)      !surface level [K]
 atmos_stag(:,:,mxlyr)%t_air = atmos(:,:,mxlyr)%t_air  !top level [K]


!
! dry air density
!
 atmos%rho_dair = 1.e+3*atmos%press*100.e0 / (const_Rd*atmos%t_air)  !dry air density [g/m3]

!
! RH
!
 call check( nf90_inq_varid(ncid, 'crm_rh' , varid ) ) !relative humidity [%]
 call check( nf90_get_var(ncid, varid, var_3d(js:je,is:ie,ks:ke), start=(/js,is,ks/), count=(/dj,di,dk/) ) )
 do k = ks,ke ; do j = js,je ; do i = is,ie
    atmos(i,j,k)%rh = var_3d(j,i,k)   ! relative humidiy [%]
 enddo ; enddo ; enddo


!
! hydrometeors
!
 mic_select0: select case(trim(cloud_microphysics))
 case('GOD','GOD10','GMP4ICE','LIN')

 call check( nf90_inq_varid(ncid, 'crm_qc' , varid ) ) !cloud water mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_3d(js:je,is:ie,ks:ke), start=(/js,is,ks/), count=(/dj,di,dk/) ) )
 do k = ks,ke ; do j = js,je ; do i = is,ie
    q_gen(i,j,k)%cloud = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3               ! cloud water mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'crm_qr' , varid ) ) !rain mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_3d(js:je,is:ie,ks:ke), start=(/js,is,ks/), count=(/dj,di,dk/) ) )
 do k = ks,ke ; do j = js,je ; do i = is,ie
    q_gen(i,j,k)%rain = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3                 ! rain water mixing ratio [g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'crm_qi' , varid ) ) !cloud ice mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_3d(js:je,is:ie,ks:ke), start=(/js,is,ks/), count=(/dj,di,dk/) ) )
 do k = ks,ke ; do j = js,je ; do i = is,ie
    q_orig = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3
    q_tune = (var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3*dqr_frac/100.) ** (dqs_frac/100.)
!    wgt1 = min( 1.0, max( 0.0, (3.0 - q_gen(i,j,1)%rain)/2.7 ) ) ; wgt2 = 1.-wgt1 !6th
    wgt1 = 1.0 ; if(q_gen(i,j,1)%rain >= 0.5) wgt1 = 0.0  ; wgt2 = 1.-wgt1 
    q_gen(i,j,k)%ice = q_orig*wgt1 + q_tune*wgt2

!    q_gen(i,j,k)%ice = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3    ! cloud ice mixing ratio [g/m3]
!     q_gen(i,j,k)%ice = (var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3*dqr_frac/100.) ** (dqs_frac/100.)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'crm_qs' , varid ) ) !snow mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_3d(js:je,is:ie,ks:ke), start=(/js,is,ks/), count=(/dj,di,dk/) ) )
 do k = ks,ke ; do j = js,je ; do i = is,ie
    q_orig = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3
    q_tune = (var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3*dqr_frac/100.) ** (dqs_frac/100.)
!    wgt1 = min( 1.0, max( 0.0, (3.0 - q_gen(i,j,1)%rain)/2.7 ) ) ; wgt2 = 1.-wgt1
    wgt1 = 1.0 ; if(q_gen(i,j,1)%rain >= 0.5) wgt1 = 0.0  ; wgt2 = 1.-wgt1
    q_gen(i,j,k)%snow = q_orig*wgt1 + q_tune*wgt2

!    q_gen(i,j,k)%snow = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3   ! snow mixing ratio [g/m3]
!    q_gen(i,j,k)%snow = (var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3*dqr_frac/100.) ** (dqs_frac/100.)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'crm_qg' , varid ) ) !graupel mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_3d(js:je,is:ie,ks:ke), start=(/js,is,ks/), count=(/dj,di,dk/) ) )
 do k = ks,ke ; do j = js,je ; do i = is,ie
    q_orig = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3
    q_tune = (var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3*dqr_frac/100.) ** (dqs_frac/100.)
!    wgt1 = min( 1.0, max( 0.0, (3.0 - q_gen(i,j,1)%rain)/2.7 ) ) ; wgt2 = 1.-wgt1
    wgt1 = 1.0 ; if(q_gen(i,j,1)%rain >= 0.5) wgt1 = 0.0  ; wgt2 = 1.-wgt1
    q_gen(i,j,k)%graupel = q_orig*wgt1 + q_tune*wgt2

!    q_gen(i,j,k)%graupel = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3   ! graupel mixing ratio [g/m3]
!    q_gen(i,j,k)%graupel = (var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3*dqr_frac/100.) ** (dqs_frac/100.)
 enddo ; enddo ; enddo

 if( trim(cloud_microphysics) == 'GMP4ICE' ) then
   call check( nf90_inq_varid(ncid, 'crm_qh' , varid ) ) !hail mixing ratio[g/kg]
   call check( nf90_get_var(ncid, varid, var_3d(js:je,is:ie,ks:ke), start=(/js,is,ks/), count=(/dj,di,dk/) ) )
   do k = ks,ke ; do j = js,je ; do i = is,ie
    q_orig = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3
    q_tune = (var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3*dqr_frac/100.) ** (dqs_frac/100.)
!    wgt1 = min( 1.0, max( 0.0, (3.0 - q_gen(i,j,1)%rain)/2.7 ) ) ; wgt2 = 1.-wgt1
    wgt1 = 1.0 ; if(q_gen(i,j,1)%rain >= 0.5) wgt1 = 0.0  ; wgt2 = 1.-wgt1
    q_gen(i,j,k)%hail = q_orig*wgt1 + q_tune*wgt2

!      q_gen(i,j,k)%hail = var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3  ! hail mixing ratio [g/m3]
!      q_gen(i,j,k)%hail = (var_3d(j,i,k)*atmos(i,j,k)%rho_dair*1.e-3*dqr_frac/100.) ** (dqs_frac/100.)
   enddo ; enddo ; enddo
 else
   q_gen%hail = 0.  
 endif


 case('WSM','WDM')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: WSM or WDM scheme does not exist in MMF')
 case('MORR','MORRH')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: Morrison scheme does not exist in MMF input')
 case('THOM')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: Thompson scheme does not exist in MMF input')
 case('NTUCLR')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: NTUCLR scheme does not exist in MMF input')
 case('NSSL2MMT') 
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: NSSL2MMT scheme does not exist in MMF input')
 case('P3')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: P3 scheme does not exist in MMF input')
 case('RAMS1','RAMS2')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: Not yet supporting RAMS microphysics in MMF input')
 case('HUCM_SBM','HUCM_SBM43')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: Not yet supporting SBM microphysics in MMF input')
 case default
    ! do nothing
 end select mic_select0

!
! vertical velocity 
!
 call check( nf90_inq_varid(ncid, 'crm_wz' , varid ) ) ! vertical velocity [m/s]
 call check( nf90_get_var(ncid, varid, var_3d(js:je,is:ie,ks:ke), start=(/js,is,ks/), count=(/dj,di,dk/) ) )
 do k = ks,ke ; do j = js,je ; do i = is,ie
    atmos_stag(i,j,k)%w = var_3d(j,i,k) ! vertical velocity [m/s]
 enddo ; enddo ; enddo
 atmos_stag(is:ie,js:je,0)%w = 0.

!
! surface mask
!
 call check( nf90_inq_varid(ncid, 'mask' , varid ) )  !
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )

 do j = js, je
    if(     INT(var_1d(j)) == 1 ) then
        surface(:,j)%iland = 1  ! GEOS land -> land
    elseif( INT(var_1d(j)) == 2 ) then
        surface(:,j)%iland = 2  ! GEOS ocean -> water
    elseif( INT(var_1d(j)) == 3 ) then
        surface(:,j)%iland = 2  ! GEOS ice sheet ->  (change this lator)
    else
        call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: strange mask value')
    endif
 enddo

!
! soil temperature
!
 surface%t_soil = surface%t_skin !

!
! snow depth
!
 call check( nf90_inq_varid(ncid, 'snow_depth' , varid ) )  ! snow depth [m]
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
    surface(:,j)%dhgt_snow = var_1d(j) !snow depth [m]
 enddo

!
! surface albedo 
!
 call check( nf90_inq_varid(ncid, 'albedo' , varid ) )  ! surface albedo [-]
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
    surface(:,j)%albedo = var_1d(j) ! surface albedo [-]
 enddo

!
! surface wind velocity 
!
 call check( nf90_inq_varid(ncid, 'wind_sfc' , varid ) )  ! surface wind [m/s]
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
    surface(:,j)%u10m = var_1d(j) ! surface wind [m/s]
 enddo

!
! u-wind 1km AGL
!
 call check( nf90_inq_varid(ncid, 'crm_u1km' , varid ) )  ! surface wind [m/s]
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
    surface(:,j)%u1km = var_1d(j) ! u-wind 1km AGL [m/s]
 enddo

!
! v-wind 1km AGL
!
 call check( nf90_inq_varid(ncid, 'crm_v1km' , varid ) )  ! surface wind [m/s]
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
    surface(:,j)%v1km = var_1d(j) ! v-wind 1km AGL [m/s]
 enddo

!
! soil moisture fraction 
!
 call check( nf90_inq_varid(ncid, 'soilH2O' , varid ) )  !soil moisture fraction [-]
 call check( nf90_get_var( ncid, varid, var_1d(js:je), start=(/js/), count=(/dj/) ) )
 do j = js, je
    surface(:,j)%h2o_soil = var_1d(j) ! soil moisture fraction [-]
 enddo

!
! vegetation fraction
!
 surface%frac_veg = 0.


!
!   close netCDF file
!
 call check( nf90_close(ncid) )


 return
 end subroutine rd_CRM_MMF_GPROF

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine rd_CRM_MMF
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
! Read MMF output in netCDF format.
! It first reads as MMF structure (4D), and transfer them into usual CRM strucure (3D) 
! 
! History:
! 12/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------

 integer :: xg,yg,j_crm,i,j,k,l  !looping
 integer :: ncid                !# of variables to inquire  
 integer :: varid               !# of variable id
 integer :: ierr                !allocation stat
 integer :: dimid,tdim
 real(sdsu_fps),allocatable :: hgt2(:) !staggered GCE height level [km]
 real(sdsu_fps),allocatable :: p0(:) !initial pressure [hPa]
 real(sdsu_fps),allocatable :: rho(:,:,:)  !moist air density [kg/m3]

 real(sdsu_fps),allocatable :: var_5d(:,:,:,:,:), var_4d(:,:,:,:),var_3d(:,:,:),var_2d(:,:),var_1d(:)
 real(sdsu_fps) :: lon_geos(mxgridx_geos),lat_geos(mxgridy_geos)

 integer :: xg_1pm, xg_1am

!----------------allocation-------------------------------
 allocate( var_5d(mxgridx_geos,mxgridy_geos,mxgridx_gce,mxgridy_gce,mxlyr_gce),& !
           var_3d(mxgridx_geos,mxgridy_geos,mxlyr_gce),& !
           var_2d(mxgridx_geos,mxgridy_geos),  & !
           rho(mxgridx,mxgridy,mxlyr),             & !
           stat=ierr )
 if(ierr /= 0) call stop_sdsu(myrank,'MSG rd_CRM_MMF: cannot allocate the variable')



!----------------------- open netCDF file
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name)
 if(masterproc) print("(a)"),'MSG rd_CRM_MMF: ++++++++++++++++++++ READING NEW FILE ++++++++++++++++++++'
 print("(a)"),'Input ->',trim(sdsu_io_file)

 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!------------- lat and lon ----------------------------------------------
 call check( nf90_inq_varid(ncid, 'lon' , varid ) ) 
 call check( nf90_get_var(ncid, varid, lon_geos ) )

!
! find the longitude of local 1:30pm
!
 xg_1pm = -999
 do xg = 1, mxgridx_geos
  if( mod(sdsu_gmt + lon_geos(xg)/15., 24.) == 13. ) then
      xg_1pm = xg ; exit
  endif
 enddo
 if(xg_1pm == -999) stop 'xg_1pm is not defined'

 xg_1am = -999
 do xg = 1, mxgridx_geos
  if( mod(sdsu_gmt + lon_geos(xg)/15., 24.) == 1. ) then  !1 am
      xg_1am = xg ; exit
  endif
 enddo
 if(xg_1am == -999) stop 'xg_1am is not defined'
 
 select case(trim(mmf_case))
 case('global')
   xs_geos = 1    ! start GEOS x grid for target domain (partial)
   xe_geos = 144  ! ending GEOS x grid for target domain (partial)
   ys_geos = 1   ! start GEOS y grid for TRMM domain (partial)
   ye_geos = 91   ! ending GEOS y grid for TRMM domain (partial)
 case('trmm')
   xs_geos = 1    ! start GEOS x grid for target domain (partial)
   xe_geos = 144  ! ending GEOS x grid for target domain (partial)
   ys_geos = 26   ! start GEOS y grid for TRMM domain (partial)
   ye_geos = 66   ! ending GEOS y grid for TRMM domain (partial)
 case('atrain_1pm')
   xs_geos = xg_1pm
   xe_geos = xg_1pm
   ys_geos = 1   ! start GEOS y grid for TRMM domain (partial)
   ye_geos = 91   ! ending GEOS y grid for TRMM domain (partial)
 case('atrain_1am')
   xs_geos = xg_1am
   xe_geos = xg_1am 
   ys_geos = 1   ! start GEOS y grid for TRMM domain (partial)
   ye_geos = 91   ! ending GEOS y grid for TRMM domain (partial)
 case default
   call stop_sdsu(myrank,'MSG rd_CRM_MMF: strange mmf_case')
 end select

 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    if( lon_geos(xg) > 180. ) then
      surface(:,j_crm)%lon =  lon_geos(xg) - 360.  !geos longitude  [deg]
    else
      surface(:,j_crm)%lon =  lon_geos(xg)  !geos longitude  [deg]
    endif
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

 call check( nf90_inq_varid(ncid, 'lat' , varid ) ) 
 call check( nf90_get_var(ncid, varid, lat_geos ) )

 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
          surface(:,j_crm)%lat =  lat_geos(yg)  !geos latitude [deg]
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg



!--------------- height & thickness --------------------------------
 call check( nf90_inq_varid(ncid, 'ht_sfc' , varid ) )  ! DEM - surface height (GEOS surface) [m]
 call check( nf90_get_var(ncid, varid, var_2d ) )

 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
          atmos_stag(:,j_crm,0)%hgt = var_2d(xg,yg)*1e-3 !surface height [km]
          surface(:,j_crm)%elev = var_2d(xg,yg) !surface height [m]
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

 call check( nf90_inq_varid(ncid, 'crm_ht' , varid ) ) !CRM height (middle) above ground [m]
 call check( nf90_get_var(ncid, varid, var_3d ) )

! toshii: change this algorithm lator with surface pressure
! find bug here. 
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    do k = 1, mxlyr-1
          atmos_stag(:,j_crm,k)%hgt = atmos_stag(:,j_crm,0)%hgt + &
                                      0.5 * ( var_3d(xg,yg,k) + var_3d(xg,yg,k+1) ) * 1e-3 ![km] 
    enddo ! k 
          atmos_stag(:,j_crm,mxlyr)%hgt = atmos_stag(:,j_crm,mxlyr-1)%hgt + &
                                     (var_3d(xg,yg,mxlyr) - var_3d(xg,yg,mxlyr-1))*1e-3 !top level [km]
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

!--------------- pressure -------------------------------

!
! middle pressure
!
 call check( nf90_inq_varid(ncid, 'crm_pres' , varid ) ) !pressure of CRM [Pa]
 call check( nf90_get_var(ncid, varid, var_3d ) )

 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    do i = 1, mxgridx
          atmos(i,j_crm,1:mxlyr_gce)%press = var_3d(xg,yg,1:mxlyr_gce)*1e-2 ! [hPa] <- [Pa]  ]
    enddo ! i 
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

!
! interface pressure
! 
 do k = 1, mxlyr-1
   atmos_stag(:,:,k)%press = 0.5e0*(atmos(:,:,k)%press + atmos(:,:,k+1)%press) !interafce [hPa]
 enddo
 atmos_stag(:,:,0)%press = atmos(:,:,1)%press + atmos(:,:,1)%press - atmos_stag(:,:,1)%press   !surface level [hPa]
 atmos_stag(:,:,mxlyr)%press = 0.5e0*(atmos(:,:,mxlyr)%press + 0.) !top pressure [hPa]

!-------------- temperaure --------------------------------

!
! skin temperature
! 
 call check( nf90_inq_varid(ncid, 'tskin' , varid ) )  ! surface skin temperature [K]
 call check( nf90_get_var(ncid, varid, var_2d ) )

 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
          surface(:,j_crm)%t_skin = var_2d(xg,yg) !surface skin temperature [K]
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

!
! air temperature (middle level)
!
 call check( nf90_inq_varid(ncid, 'crm_ta' , varid ) ) !temperature of [K]
 call check( nf90_get_var(ncid, varid, var_5d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg =xs_geos, xe_geos
    do i = 1, mxgridx
          atmos(i,j_crm,1:mxlyr)%t_air =  var_5d(xg,yg,i,1,1:mxlyr)  !temperature  [K]
    enddo ! i 
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

!
! air temperature (interface)
!
 do k = 1, mxlyr-1
    atmos_stag(:,:,k)%t_air = 0.5e0 * (atmos(:,:,k)%t_air + atmos(:,:,k+1)%t_air) ![K]
 enddo
 atmos_stag(:,:,0)%t_air     = 0.5e0* (surface(:,:)%t_skin+atmos_stag(:,:,1)%t_air)      !surface level [K]
 atmos_stag(:,:,mxlyr)%t_air = atmos(:,:,mxlyr)%t_air  !top level [K]


!-------------- air density---------------------------------------
 rho = atmos%press*100.e0 / (const_Rd*atmos%t_air)  !dry air density [kg/m3]
 atmos%rho_dair = 1.e+3*atmos%press*100.e0 / (const_Rd*atmos%t_air)  !dry air density [g/m3]

! do k = 1,mxlyr ; do j = 1, mxgridy ; do i = 1, mxgridx
!    if( atmos(i,j,k)%t_air < 0. ) then
!      print*,'press is negative',  atmos(i,j,k)%t_air  ,i,j,k
!    endif
! enddo ; enddo ; enddo

!-------------- RH ---------------------------------------
 call check( nf90_inq_varid(ncid, 'crm_rh' , varid ) ) !relative humidity [%]
 call check( nf90_get_var(ncid, varid, var_5d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    do i = 1, mxgridx
          atmos(i,j_crm,1:mxlyr)%rh =  var_5d(xg,yg,i,1,1:mxlyr)  !relative humidiy [%]
    enddo ! i 
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg



!----------------hydrometeors---------------------------------------
 mic_select0: select case(trim(cloud_microphysics))
 case('GOD','GOD10','GMP4ICE','LIN')

 call check( nf90_inq_varid(ncid, 'crm_qc' , varid ) ) !cloud water mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_5d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    do i = 1, mxgridx
          q_gen(i,j_crm,:)%cloud =  var_5d(xg,yg,i,1,:)*rho(i,j_crm,:)  !cloud water mixing ratio [g/m3]
    enddo ! i
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg


 call check( nf90_inq_varid(ncid, 'crm_qr' , varid ) ) !rain mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_5d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    do i = 1, mxgridx
          q_gen(i,j_crm,:)%rain =  var_5d(xg,yg,i,1,:)*rho(i,j_crm,:)  !rain water mixing ratio [g/m3]
    enddo ! i 
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg


 call check( nf90_inq_varid(ncid, 'crm_qi' , varid ) ) !cloud ice mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_5d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    do i = 1, mxgridx
          q_gen(i,j_crm,:)%ice =  var_5d(xg,yg,i,1,:)*rho(i,j_crm,:)  !cloud ice mixing ratio [g/m3]
    enddo ! i 
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg




 call check( nf90_inq_varid(ncid, 'crm_qs' , varid ) ) !snow mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_5d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    do i = 1, mxgridx
          q_gen(i,j_crm,:)%snow =  var_5d(xg,yg,i,1,:)*rho(i,j_crm,:)   !snow mixing ratio [g/m3]
    enddo ! i 
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg


 call check( nf90_inq_varid(ncid, 'crm_qg' , varid ) ) !graupel mixing ratio[g/kg]
 call check( nf90_get_var(ncid, varid, var_5d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    do i = 1, mxgridx
          q_gen(i,j_crm,:)%graupel  =  var_5d(xg,yg,i,1,:)*rho(i,j_crm,:)  !graupel mixing ratio [g/m3]
    enddo ! i 
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

 if( trim(cloud_microphysics) == 'GMP4ICE' ) then
   call check( nf90_inq_varid(ncid, 'crm_qh' , varid ) ) !hail mixing ratio[g/kg]
   call check( nf90_get_var(ncid, varid, var_5d ) )
   j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
      do i = 1, mxgridx
          q_gen(i,j_crm,:)%hail  =  var_5d(xg,yg,i,1,:)*rho(i,j_crm,:)  !hail mixing ratio [g/m3]
      enddo ! i 
   j_crm = j_crm + 1 ; enddo ; enddo !xg yg
 else
   q_gen%hail = 0.  !currently zero 
 endif


 case('WSM','WDM')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF: WSM or WDM scheme does not exist in MMF')
 case('MORR','MORRH')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF: Morrison scheme does not exist in MMF input')
 case('NSSL2MMT') 
   call stop_sdsu(myrank,'MSG rd_CRM_MMF: NSSL2MMT scheme does not exist in MMF input')
 case('P3')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF: P3 scheme does not exist in MMF input')
 case('THOM')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: Thompson scheme does not exist in MMF input')
 case('NTUCLR')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF_GPROF: NTUCLR scheme does not exist in MMF input')
 case('RAMS1','RAMS2')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF: Not yet supporting RAMS microphysics in MMF input')
 case('HUCM_SBM','HUCM_SBM43')
   call stop_sdsu(myrank,'MSG rd_CRM_MMF: Not yet supporting SBM microphysics in MMF input')
 case default
    ! do nothing
 end select mic_select0



!---------------- Misc parameters--------------------------------------------


!
! land-ocean mask
!
 call check( nf90_inq_varid(ncid, 'mask' , varid ) )  !
 call check( nf90_get_var(ncid, varid, var_2d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
          if( INT(var_2d(xg,yg)) == 1 ) then
                 surface(:,j_crm)%iland = 1  ! GEOS land -> land
          elseif( INT(var_2d(xg,yg)) == 2 ) then 
                 surface(:,j_crm)%iland = 2  ! GEOS ocean -> water
          elseif( INT(var_2d(xg,yg)) == 3 ) then
              surface(:,j_crm)%iland = 2  ! GEOS ice sheet ->  (change this lator)
          else
              print*,xg,yg,INT(var_2d(xg,yg))
              call stop_sdsu(myrank,'MSG rd_CRM_MMF: strange mask value')
          endif 
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

!
! soil temperature
!
 surface%t_soil = surface%t_skin !
! call check( nf90_inq_varid(ncid, 'tsoil' , varid ) )  ! surface skin temperature [K]
! call check( nf90_get_var(ncid, varid, var_2d ) )
! j_crm = 1 ; do xg = 1,mxgridx_geos ; do yg = 1,mxgridy_geos
!       if(iland(1,j_crm) == 1) then
!          surface(:,j_crm)%t_soil = var_2d(xg,yg) !top-soil temperature [K]
!       else
!          surface(:,j_crm)%t_soil = surface(:,j_crm)%t_skin !skin temp 
!       endif
! j_crm = j_crm + 1 ; enddo ; enddo !xg yg


!
! snow depth
!
 call check( nf90_inq_varid(ncid, 'snow_depth' , varid ) )  ! snow depth [m]
 call check( nf90_get_var(ncid, varid, var_2d ) )

 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
          surface(:,j_crm)%dhgt_snow = var_2d(xg,yg) !snow depth [m]
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

!
! surface albedo 
!
 call check( nf90_inq_varid(ncid, 'albedo' , varid ) )  ! surface albedo [-]
 call check( nf90_get_var(ncid, varid, var_2d ) )

 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
          surface(:,j_crm)%albedo = var_2d(xg,yg) !surface albedo [-]
          if(surface(1,j_crm)%albedo > 1.)  surface(:,j_crm)%albedo = undefined !night time..
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

!
! surface wind velocity 
!
 call check( nf90_inq_varid(ncid, 'wind_sfc' , varid ) )  ! surface wind [m/s]
 call check( nf90_get_var(ncid, varid, var_2d ) )

 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
          surface(:,j_crm)%u10m = var_2d(xg,yg) !surface wind [m/s]
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg


!
! soil moisture fraction 
!
 call check( nf90_inq_varid(ncid, 'soilH2O' , varid ) )  !soil moisture fraction [-]
 call check( nf90_get_var(ncid, varid, var_2d ) )

 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
          surface(:,j_crm)%h2o_soil = var_2d(xg,yg) ! soil moisture fraction [-]
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

!
! Vegetation Fraction  
!
! call check( nf90_inq_varid(ncid, 'lai' , varid ) )  !leaf area index [-]
! call check( nf90_get_var(ncid, varid, var_2d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
         if(surface(1,j_crm)%iland == 1) then !land 

!          surface(:,j_crm)%frac_veg = 1.- exp(-var_2d(xg,yg)) ! vegetation fraction  [-]  (Temporal)
          !
          ! consistent to NESDIS algorithm 
          !
!          surface(:,j_crm)%frac_veg = MAX( 2.e0, ( var_2d(xg,yg) - 0.5e0 ) / 3.e0 * 100.e0) 
              surface(:,j_crm)%frac_veg = 0.5  !toshii (temporal)
         else !ocean or ice-sheet
              surface(:,j_crm)%frac_veg = 0.
         endif
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

!
! vertical velocity 
!
 atmos_stag%w = 5. !toshi

 call check( nf90_inq_varid(ncid, 'crm_wz' , varid ) ) !vertical velocity [m/s]
 call check( nf90_get_var(ncid, varid, var_5d ) )
 j_crm = 1 ; do yg = ys_geos, ye_geos ; do xg = xs_geos, xe_geos
    do i = 1, mxgridx
          atmos(i,j_crm,:)%w =  var_5d(xg,yg,i,1,:) !vertical velocity [m/s]
    enddo ! i
 j_crm = j_crm + 1 ; enddo ; enddo !xg yg

 do j = 1, mxgridy ; do i = 1, mxgridx
    atmos_stag(i,j,0)%w = 0.
    do k = 1,mxlyr
       atmos_stag(i,j,0)%w = ( atmos_stag(i,j,k)%w + atmos_stag(i,j,k-1)%w ) * 0.5e0 ![m/s]
    enddo
 enddo ; enddo

!-------------------------- close netCDF file

  call check( nf90_close(ncid) )


 return
 end subroutine rd_CRM_MMF

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine rd_CRM_GEOS5
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
! Read GEOS5 output in netCDF3 format.
! This is tempral from Oreste's run....
! 
! History:
! 05/2011  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer :: i,j,k,l  !looping
 integer :: ncid                !# of variables to inquire  
 integer :: varid               !# of variable id
 integer :: ierr                !allocation stat
 integer :: dimid,tdim

 real(sdsu_fps),allocatable :: hgt2(:)    ! staggered GCE height level [km]
 real(sdsu_fps),allocatable :: rho(:,:,:) ! dry air density [kg/m3]

 integer,parameter :: mxgridx_geos5 = 1080
 integer,parameter :: mxgridy_geos5 =  721
 real(sdsu_fpd) :: lon_geos5(mxgridx_geos5),lat_geos5(mxgridy_geos5)
 integer,parameter :: mxlyr_p_geos5 =   36
 real(sdsu_fpd) :: p0(mxlyr_p_geos5)  ! initial pressure [hPa]

 integer :: ids, ide, di
 integer :: jds, jde, dj
 integer :: kds, kde, dk
 integer :: ig, jg, kg

 real(sdsu_fps) :: lat_min =   0.
 real(sdsu_fps) :: lat_max =  30.
 real(sdsu_fps) :: lon_min = -40.
 real(sdsu_fps) :: lon_max =  20.

 real(sdsu_fps),allocatable :: net2d(:,:)
 real(sdsu_fps),allocatable :: temp2d(:,:)
 real(sdsu_fps),allocatable :: net3d(:,:,:)
 real(sdsu_fps),allocatable,dimension(:,:,:) :: airdens, du, ple, dumass

 integer :: ilen
 character(len=2) :: hh



!
! Open netCDF file
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name)
 if(masterproc) print("(a)"),'MSG rd_CRM_GEOS5: ++++++++++++++++++++ READING 3D FILE ++++++++++++++++++++'
 print("(a)"),'Input ->',trim(sdsu_io_file)

 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )


 call check( nf90_inq_varid(ncid, 'lon' , varid ) )
 call check( nf90_get_var(ncid, varid, lon_geos5 ) )

 call check( nf90_inq_varid(ncid, 'lat' , varid ) )
 call check( nf90_get_var(ncid, varid, lat_geos5 ) )

 call check( nf90_inq_varid(ncid, 'lev' , varid ) )
 call check( nf90_get_var(ncid, varid, p0 ) )

!
! NAMMA domain
!
 do i = 1, mxgridx_geos5
    if( lon_geos5(i) >= REAL(lon_min)) then 
        ids = i ; exit
    endif
 enddo
 do i = 1, mxgridx_geos5
    if( lon_geos5(i) >= REAL(lon_max) ) then
        ide = i ; exit
    endif
 enddo
 print*,ids,ide, 'mxgridx =', ide-ids + 1
 di = ide-ids + 1

 do j = 1, mxgridy_geos5
    if( lat_geos5(j) >= REAL(lat_min) ) then
        jds = j ; exit
    endif
 enddo
 do j = 1, mxgridy_geos5
    if( lat_geos5(j) >= REAL(lat_max) ) then
        jde = j ; exit
    endif
 enddo
 print*,jds,jde, 'mxgridy =', jde-jds + 1
 dj = jde-jds + 1


 kds=1 ; kde=mxlyr
 dk = kde-kds + 1

!
! get lat lon for SDSU
!
 do i = 1, mxgridx
    surface(i,:)%lon = REAL(lon_geos5(i+ids-1))
 enddo 

 do j = 1, mxgridy
    surface(:,j)%lat = REAL(lat_geos5(j+jds-1))
 enddo 

!
! memory allocation
!
 allocate( &
           net2d     (ids:ide,jds:jde          ), &
           temp2d    (ids:ide,jds:jde          ), &
           net3d     (ids:ide,jds:jde,kds:kde  ), &
           rho       (ids:ide,jds:jde,kds:kde  ), &
           stat=ierr )

 do j = 1, mxgridy ; do i = 1, mxgridx
      atmos(i,j,1:mxlyr)%press =  p0(1:mxlyr)  ![hPa]
 enddo ; enddo

!
! Surface Pressure
!
 call check( nf90_inq_varid(ncid, 'PS' , varid ) )  ! Pa
 call check( nf90_get_var(ncid, varid, net2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))

 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    atmos_stag(i,j,0)%press = net2d(ig,jg) * 0.01  ![hPa] <-- [Pa]
    do k = 1,mxlyr-1
       atmos_stag(i,j,k)%press = 0.5*( atmos(i,j,k)%press + atmos(i,j,k+1)%press )  ![hPA]
    enddo
    atmos_stag(i,j,mxlyr)%press = 0.1  !TOA [hPa]
 enddo ; enddo
 
!
! Geopotential height
!
 call check( nf90_inq_varid(ncid, 'H' , varid ) ) !  [m]
 call check( nf90_get_var(ncid, varid, net3d(ids:ide,jds:jde,kds:kde), start=(/ids,jds,kds/), count=(/di,dj,dk/)))

 do k = 1, mxlyr ; do j = 1, mxgridy ; do i = 1, mxgridx
    kg = k ; jg = j+jds-1 ; ig = i+ids-1 
    atmos(i,j,k)%hgt = net3d(ig,jg,kg)*1.e-3  ![km] 
 enddo ; enddo ; enddo

!
! Air temperature
!
 call check( nf90_inq_varid(ncid, 'T' , varid ) ) !  [K]
 call check( nf90_get_var(ncid, varid, net3d(ids:ide,jds:jde,kds:kde), start=(/ids,jds,kds/), count=(/di,dj,dk/)))

 do k = 1, mxlyr ; do j = 1, mxgridy ; do i = 1, mxgridx
    kg = k ; jg = j+jds-1 ; ig = i+ids-1
    atmos(i,j,k)%t_air = net3d(ig,jg,kg)  ![K] 
 enddo ; enddo ; enddo

!
! Dry air density
!
 do k = 1, mxlyr ; do j = 1, mxgridy ; do i = 1, mxgridx
     kg = k ; jg = j+jds-1 ; ig = i+ids-1
     rho(ig,jg,kg) = (atmos(i,j,k)%press*100.e0) / (const_Rd*atmos(i,j,k)%t_air)  !dry air density [kg/m3]
 enddo ; enddo ; enddo
 atmos%rho_dair = 1.e+3*atmos%press*100.e0 / (const_Rd*atmos%t_air)  !dry air density [g/m3]

!
! layer thickness
!



!
! Relative humidity 
!
 call check( nf90_inq_varid(ncid, 'RH' , varid ) ) !  [%]
 call check( nf90_get_var(ncid, varid, net3d(ids:ide,jds:jde,kds:kde), start=(/ids,jds,kds/), count=(/di,dj,dk/)))

 do k = 1, mxlyr ; do j = 1, mxgridy ; do i = 1, mxgridx
    kg = k ; jg = j+jds-1 ; ig = i+ids-1
    atmos(i,j,k)%rh = net3d(ig,jg,kg)  ![%] 
 enddo ; enddo ; enddo


!
! Specific humidity 
!
 call check( nf90_inq_varid(ncid, 'QV' , varid ) ) !  [g/g/]
 call check( nf90_get_var(ncid, varid, net3d(ids:ide,jds:jde,kds:kde), start=(/ids,jds,kds/), count=(/di,dj,dk/)))
 atmos%qv = net3d  ! water vapor mixing ratio [kg/kg]

 do k = 1, mxlyr ; do j = 1, mxgridy ; do i = 1, mxgridx
    kg = k ; jg = j+jds-1 ; ig = i+ids-1
    atmos(i,j,k)%sh = net3d(ig,jg,kg)  ![g/g] 
 enddo ; enddo ; enddo

!
! total condensate --> cloud and ice
!
 call check( nf90_inq_varid(ncid, 'QC' , varid ) ) !  [g/g]
 call check( nf90_get_var(ncid, varid, net3d(ids:ide,jds:jde,kds:kde), start=(/ids,jds,kds/), count=(/di,dj,dk/)))

 do k = 1, mxlyr ; do j = 1, mxgridy ; do i = 1, mxgridx
    kg = k ; jg = j+jds-1 ; ig = i+ids-1

    if( atmos(i,j,k)%t_air > 273.14 ) then !warm phase
       q_gen(i,j,k)%cloud = net3d(ig,jg,kg) * rho(ig,jg,kg) * 1.e3
       q_gen(i,j,k)%ice   = 0.e0
    else !cold phase
       q_gen(i,j,k)%cloud = 0.e0 
       q_gen(i,j,k)%ice   = net3d(ig,jg,kg) * rho(ig,jg,kg) * 1.e3
    endif

 enddo ; enddo ; enddo

! Note that raining species are always zero for GEOS5. 
 q_gen%rain    = 0.e0 
 q_gen%snow    = 0.e0 
 q_gen%graupel = 0.e0
 q_gen%hail    = 0.e0

!
! close netCDF file
!
  call check( nf90_close(ncid) )


 ilen= len(trim(sdsu_inp_name))

 sdsu_io_file=trim(sdsu_dir_input)//sdsu_inp_name(1:18)//'sfc'//sdsu_inp_name(31:ilen)
 if(masterproc) print("(a)"),'MSG rd_CRM_GEOS5: ++++++++++++++++++++ READING SFC FILE ++++++++++++++++++++'
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )


!
! Surface Geopotential height
!
 call check( nf90_inq_varid(ncid, 'PHIS' , varid ) )  ! m2/s2
 call check( nf90_get_var(ncid, varid, net2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))

 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    atmos_stag(i,j,0)%hgt = net2d(ig,jg)/const_g/1.e3  ! [km] 
    surface(i,j)%elev     = net2d(ig,jg)/const_g       ! [m]

    do k = 1, mxlyr-1
       atmos_stag(i,j,k)%hgt = 0.5 * ( atmos(i,j,k)%hgt + atmos(i,j,k+1)%hgt ) ![km]
    enddo
    atmos_stag(i,j,mxlyr)%hgt =  2.*atmos(i,j,mxlyr)%hgt - atmos_stag(i,j,mxlyr-1)%hgt  ![km]

    do k = 1,mxlyr
       atmos(i,j,k)%dhgt = atmos_stag(i,j,k)%hgt - atmos_stag(i,j,k-1)%hgt ![km]
       if( atmos(i,j,k)%dhgt < 0.e0 ) print*,'dhgt negative',i,j,k,atmos(i,j,k)%dhgt
    enddo

 enddo ; enddo


!
! Surface parameters
!
 call check( nf90_inq_varid(ncid, 'TSKIN' , varid ) )  ! surface skin temperature [K]
 call check( nf90_get_var(ncid, varid, net2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))
 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    surface(i,j)%t_skin  = net2d(ig,jg)   ! [K]
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'TA' , varid ) )  ! surface air temperature [K]
 call check( nf90_get_var(ncid, varid, net2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))
 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    surface(i,j)%t_air  = net2d(ig,jg)   ! [K]
    atmos_stag(i,j,0)%t_air = net2d(ig,jg)   ! [K]

    do k = 1, mxlyr-1
       atmos_stag(i,j,k)%t_air = 0.5*( atmos(i,j,k)%t_air + atmos(i,j,k+1)%t_air )
    enddo
    atmos_stag(i,j,mxlyr)%t_air = atmos(i,j,mxlyr)%t_air   

 enddo ; enddo


 call check( nf90_inq_varid(ncid, 'VA' , varid ) )  ! surface v wind [m/s]
 call check( nf90_get_var(ncid, varid, net2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'UA' , varid ) )  ! surface u wind [m/s]
 call check( nf90_get_var(ncid, varid, temp2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))
 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    surface(i,j)%u10m  = sqrt( net2d(ig,jg)**2 + temp2d(ig,jg)**2 )   ! [m/s]
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SNOWDP' , varid ) )  ! snow water equivalent depth [mm]
 call check( nf90_get_var(ncid, varid, net2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))
 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    if( net2d(ig,jg) /= 1.e+15 ) then
      surface(i,j)%h2o_snow  = net2d(ig,jg)   !snow water equivalent [kg/m2] == [ mm ]
      surface(i,j)%dhgt_snow = surface(i,j)%h2o_snow / 0.3 / 1000.0  !snow depth [m] 
    else
      surface(i,j)%h2o_snow  = 0.   !snow water equivalent [kg/m2] == [ mm ]
      surface(i,j)%dhgt_snow = 0.  !snow depth [m] 
    endif
 enddo ; enddo


 call check( nf90_inq_varid(ncid, 'GWETTOP' , varid ) )  ! top-soil wetness fraction [-]
 call check( nf90_get_var(ncid, varid, net2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))
 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    if( net2d(ig,jg) /= 1.e+15 ) then
      surface(i,j)%h2o_soil   = net2d(ig,jg)   !soil moisture frac [-]
    else
      surface(i,j)%h2o_soil   = 0.   !soil moisture frac [-]
    endif
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'TSOIL1' , varid ) )  ! top-soil temperature [K]
 call check( nf90_get_var(ncid, varid, net2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))
 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    if( net2d(ig,jg) /= 1.e+15 ) then
      surface(i,j)%t_soil   = net2d(ig,jg)   !top-soil temperature [K]
    else
      surface(i,j)%t_soil   = undefined   !top-soil temperature [K]
    endif
 enddo ; enddo


 call check( nf90_inq_varid(ncid, 'FRLAND' , varid ) )  ! land fraction [-0]
 call check( nf90_get_var(ncid, varid, net2d(ids:ide,jds:jde), start=(/ids,jds/), count=(/di,dj/)))
 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    if( net2d(ig,jg) > 0.5 ) then
       surface(i,j)%iland    = 1   !land index [-]
       surface(i,j)%igbp_typ = 9   !land --> savanna (modify lator)
    else
       surface(i,j)%iland    = 2   !ocean index [-]
       surface(i,j)%igbp_typ = 0   !water body (modify lator)
    endif
 enddo ; enddo

! other parameters
  surface%albedo = 0.1
  surface%frac_veg = 0.1
  surface%rain_rate = 0.


 deallocate( &
           net2d     , &
           temp2d    , &
           net3d     , &
           rho       , &
           stat=ierr )


 if(account_aerosol) then

    q_gocart%so4=0. ; q_gocart%blc=0. ; q_gocart%ocn=0. ; q_gocart%och=0.
    q_gocart%ssa=0. ; q_gocart%ssc=0. ; q_gocart%du1=0. ; q_gocart%du2=0.
    q_gocart%du3=0. ; q_gocart%du4=0. ; q_gocart%du5=0. ; q_gocart%du6=0.
    q_gocart%du7=0. ; q_gocart%du8=0.


 ilen= len(trim(sdsu_inp_name))
 write(hh,"(I2.2)") INT(sdsu_hh + 1.)
 sdsu_io_file=trim(sdsu_dir_input)//sdsu_inp_name(1:13)//'tavg3d_aer_p'//sdsu_inp_name(31:ilen-10) &
              //hh//'30z.nc4.nc'
 if(masterproc) print("(a)"),'MSG rd_CRM_GEOS5: ++++++++++++++++++++ READING AER3D FILE ++++++++++++++++++++'
 print("(a)"),'Input ->',trim(sdsu_io_file)
 call check_io( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) ,bad_sdsu_io)

!
! If aerosol data does not exist return
!
 if(bad_sdsu_io) return

!
! If present, memory allocation for aerosol parameters
!
 kds = 1 ; kde = 26
 dk = kde-kds + 1

 allocate( &
           airdens   (ids:ide,jds:jde,kds:kde  ), &
           du        (ids:ide,jds:jde,kds:kde  ), &
           ple       (ids:ide,jds:jde,kds:kde  ), &
           dumass    (ids:ide,jds:jde,kds:kde  ), &
           stat=ierr )

 call check( nf90_inq_varid(ncid, 'AIRDENS' , varid ) ) !  [kg/m3]
 call check( nf90_get_var(ncid, varid, airdens(ids:ide,jds:jde,kds:kde), start=(/ids,jds,kds/), count=(/di,dj,dk/)))
 call check( nf90_inq_varid(ncid, 'DU' , varid ) ) ! dust [kg/kg]
 call check( nf90_get_var(ncid, varid, du(ids:ide,jds:jde,kds:kde), start=(/ids,jds,kds/), count=(/di,dj,dk/)))
 call check( nf90_inq_varid(ncid, 'PLE' , varid ) ) !  [Pa]
 call check( nf90_get_var(ncid, varid, ple(ids:ide,jds:jde,kds:kde), start=(/ids,jds,kds/), count=(/di,dj,dk/)))

!
! create interpolation here....
!
 dumass = airdens * du * 1.e3  ![g/m3]
 do j = 1, mxgridy ; do i = 1, mxgridx
    jg = j+jds-1 ; ig = i+ids-1
    do k = 1, mxlyr
       do kg = kds, kde
          if( INT(atmos(i,j,k)%press) == INT(ple(ig,jg,kg)/100.) ) then
             q_gocart(i,j,k)%du6 = dumass(ig,jg,kg)
             exit
          endif
       enddo
    enddo
    q_gocart(i,j, 6)%du6 = 0.5*(dumass(ig,jg,5) + dumass(ig,jg,6) )
    q_gocart(i,j, 8)%du6 = 0.5*(dumass(ig,jg,6) + dumass(ig,jg,7) )
    q_gocart(i,j,26)%du6 = 0.5*(dumass(ig,jg,23) + dumass(ig,jg,24) )
 enddo ; enddo


 deallocate( &
           airdens   , &
           du        , &
           ple       , &
           dumass    , &
           stat=ierr )


 endif !account aerosol


 return
 end subroutine rd_CRM_GEOS5

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine rd_CRM_GEOS5_forecast_latlon
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!    Read GEOS5 (Fortuna) Forecast output in NetCDF4 (HDF5) format.
! 
! History:
!   11/2012  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer :: i,j,k,l  !looping
 integer :: ncid                !# of variables to inquire  
 integer :: varid               !# of variable id
 integer :: ierr                !allocation stat
 integer :: dimid,tdim
 integer :: ig, jg, kg
 real(sdsu_fps),allocatable :: rho(:,:,:) ! dry air density [kg/m3]
 real(sdsu_fps),allocatable :: temp2d(:,:)
 real(sdsu_fps),allocatable :: temp3d(:,:,:)
 real(sdsu_fps),allocatable,dimension(:,:,:) ::  du, ple, dumass
 integer :: ilen
 character(len=2) :: hh
 integer :: is,ie,js,je,ks,ke,di,dj,dk
 real(sdsu_fps),parameter :: geos_undef = 1.e+15
 integer :: igs,ige,jgs,jge,kgs,kge,dig,djg,dkg
 character(len=200) :: aerosol_inp_name

!
! GEOS5 options
!
! type geos5_grid
!  integer :: imax   ! maximum i dimension of GEOS5
!  integer :: jmax   ! maximum j dimension of GEOS5
!  integer :: kmax   ! maximum k dimension of GEOS5
!  integer :: ids, ide, di
!  integer :: jds, jde, dj
!  integer :: kds, kde, dk
! end type geos5_grid
! type( geos5_grid ) :: geos  !geos5 grid group

! real(sdsu_fpd),allocatable :: lon_geos5(:)    ! full longitude GEOS5 [deg]
! real(sdsu_fpd),allocatable :: lat_geos5(:)    ! full latitude GEOS5 [deg]
! real(sdsu_fpd),allocatable :: lev_geos5(:)    ! full eta vertical level [-]
! real(sdsu_fps) :: dlon_geos5                  ! longtidue grid space [deg]
! real(sdsu_fps) :: dlat_geos5                  ! latitdue grid space [deg]


!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1

!
! define GEOS tile index
!
 igs = is + geos%ids - 1  ;  ige = ie + geos%ids - 1
 jgs = js + geos%jds - 1  ;  jge = je + geos%jds - 1
 kgs = geos%kds           ;  kge = geos%kde

 if( .not. allocated(temp2d) ) then
    allocate( temp2d(is:ie,js:je      ) , &
              temp3d(is:ie,js:je,ks:ke) , &
               rho  (is:ie,js:je,ks:ke) , &
              stat=ierr )
 endif




!
! Open netCDF file
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name)
 if(masterproc) print("(a)"),'MSG rd_CRM_GEOS5_forecast_latlon: ++++++++++++++++++++ READING 3D FILE ++++++++++++++++++++'
 if(masterproc) print("(a)"),'Input ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! set up lat lon
!
 do i = is , ie
    ig = i + geos%ids - 1
    surface(i,js:je)%lon = real( lon_geos5(ig) )
 enddo

 do j = js , je
    jg = j + geos%jds - 1
    surface(is:ie,j)%lat = real( lat_geos5(jg) )
 enddo

!
! pressure
!
 call check( nf90_inq_varid(ncid, 'PL', varid ) ) ! Mid-level pressure [Pa]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    atmos(i,j,k)%press = temp3d(i,j,kg)*0.01         ![hPa]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'PS', varid ) ) ! Surface pressure [Pa]
 call check( nf90_get_var(ncid, varid, temp2d(is:ie,js:je), start=(/igs,jgs/), count=(/di,dj/)))
 do j = js,je ; do i = is,ie
    atmos_stag(i,j,0)%press = temp2d(i,j)*0.01         ![hPa]
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DELP', varid ) ) ! Pressure thickness [Pa]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    atmos_stag(i,j,k)%press = atmos_stag(i,j,0)%press - SUM( temp3d(i,j,kg:kge) )*0.01
 enddo ; enddo ; enddo

 !must clean up for stratosphere
 do j = js,je ; do i = is,ie
    do k = ks+1,ke
       if( atmos_stag(i,j,k-1)%press <= atmos(i,j,k)%press ) then
           atmos_stag(i,j,k-1)%press = 0.5*( atmos(i,j,k)%press + atmos(i,j,k-1)%press )
       endif
    enddo
 enddo ; enddo


!
! Geopotential height
!
 call check( nf90_inq_varid(ncid, 'H' , varid ) ) !  [m]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    atmos(i,j,k)%hgt = temp3d(i,j,kg)*1.e-3   ![km]
 enddo ; enddo ; enddo


!
! Air temperature
!
 call check( nf90_inq_varid(ncid, 'T' , varid ) ) ![K]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    atmos(i,j,k)%t_air = temp3d(i,j,kg)         ![K]
 enddo ; enddo ; enddo

 atmos_stag(:,:,0)%t_air = atmos(:,:,1)%t_air   !lowest level
 do k = ks+1,ke
    atmos_stag(:,:,k)%t_air = 0.5*(atmos(:,:,k-1)%t_air + atmos(:,:,k)%t_air)
 enddo

!
! Relative Humidity 
!
 call check( nf90_inq_varid(ncid, 'RH' , varid ) ) ![-]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    atmos(i,j,k)%rh = temp3d(i,j,kg)*100.         ![%]
 enddo ; enddo ; enddo

!
! Specific Humidity
!
 call check( nf90_inq_varid(ncid, 'QV' , varid ) ) ![g/g]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 atmos%qv = temp3d  ! water vapor mixing ratio [kg/kg]

 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    atmos(i,j,k)%sh = temp3d(i,j,kg)               ![g/g]
 enddo ; enddo ; enddo

!
! Dry air density
!
 do k = ks,ke ; do j = js,je ; do i = is,ie
     rho(i,j,k) = (atmos(i,j,k)%press*100.e0) / (const_Rd*atmos(i,j,k)%t_air)  !dry air density [kg/m3]
 enddo ; enddo ; enddo
  atmos%rho_dair = 1.e+3*atmos%press*100.e0 / (const_Rd*atmos%t_air)  !dry air density [g/m3]



!
! cloud liquid water content
!
 call check( nf90_inq_varid(ncid, 'QL' , varid ) )           ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gen(i,j,k)%cloud = temp3d(i,j,kg) * rho(i,j,k) * 1.e3  ![g/m3]
 enddo ; enddo ; enddo

!
! cloud ice water content
!
 call check( nf90_inq_varid(ncid, 'QI' , varid ) )           ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gen(i,j,k)%ice = temp3d(i,j,kg) * rho(i,j,k) * 1.e3  ![g/m3]
 enddo ; enddo ; enddo

!
! other spceices is zero. 
!
 q_gen%snow = 0.e0 ; q_gen%graupel = 0.e0  ; q_gen%hail = 0.e0 

!
! Layer thickness
!
 do k = ks,ke ; do j = js,je ; do i = is,ie
     atmos(i,j,k)%dhgt = 100.e0 * (atmos_stag(i,j,k-1)%press  -  atmos_stag(i,j,k)%press) &
                        / ( rho(i,j,k) * const_g ) * 1.e-3  ![km]
 enddo ; enddo ; enddo


!
! Surface Geopotential height
!
 call check( nf90_inq_varid(ncid, 'PHIS' , varid ) )  ! m2/s2
 call check( nf90_get_var(ncid, varid, temp2d(is:ie,js:je), start=(/igs,jgs/), count=(/di,dj/)))
 do j = js,je ; do i = is,ie
    atmos_stag(i,j,0)%hgt = temp2d(i,j)/const_g/1.e3  ! [km] 
    surface(i,j)%elev     = temp2d(i,j)/const_g       ! [m]
 enddo ; enddo

!
! Stagger height
!
 do j = js,je ; do i = is,ie
    do k = ks,ke
       atmos_stag(i,j,k)%hgt = atmos_stag(i,j,0)%hgt + SUM( atmos(i,j,ks:k)%dhgt )
    enddo
 enddo ; enddo

 do j = js,je ; do i = is,ie  !filter
    do k = ks,ke-1
      if ( atmos(i,j,k)%hgt >= atmos_stag(i,j,k)%hgt ) then
         atmos_stag(i,j,k)%hgt = 0.5*(atmos(i,j,k)%hgt+atmos(i,j,k+1)%hgt)
      endif
    enddo   
 enddo ; enddo


!
!==================== Surface Field ======================================
!
! modify when you get surface data. 
!
 surface(:,:)%t_skin = atmos(:,:,1)%t_air
 surface(:,:)%h2o_snow  = 0.   !snow water equivalent [kg/m2] == [ mm ]
 surface(:,:)%dhgt_snow = 0.   !snow depth [m] 
 surface(:,:)%h2o_soil  = 0.
 surface(:,:)%t_soil    = atmos(:,:,1)%t_air
 surface(:,:)%albedo    = 0.1
 surface(:,:)%frac_veg  = 0.1
 surface(:,:)%rain_rate = 0.
 surface(:,:)%u10m      = 5.

 do j = js,je ; do i = is,ie  
    if( surface(i,j)%elev == 0. ) then
        surface(i,j)%iland    = 2 !ocean
        surface(i,j)%igbp_typ = 0 !water body
    else
        surface(i,j)%iland    = 1 !land
        surface(i,j)%igbp_typ = 9 !land --> savanna (modify lator)
    endif
 enddo ; enddo



!
! close netCDF file
!
  call check( nf90_close(ncid) )

!
! =================== Aerosol Fields ====================================
!

 if(account_aerosol) then

    q_gocart%so4=0. ; q_gocart%blc=0. ; q_gocart%ocn=0. ; q_gocart%och=0.
    q_gocart%ssa=0. ; q_gocart%ssc=0. ; q_gocart%du1=0. ; q_gocart%du2=0.
    q_gocart%du3=0. ; q_gocart%du4=0. ; q_gocart%du5=0. ; q_gocart%du6=0.
    q_gocart%du7=0. ; q_gocart%du8=0.

!
! Open netCDF file
!
 call get_aerosolfile_name( trim(sdsu_inp_name), aerosol_inp_name )
 sdsu_io_file = trim(sdsu_dir_input)//trim(aerosol_inp_name)
 if(masterproc) print("(a)"),'Open aerosol file ->',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )


!
! aerosols 
!
 call check( nf90_inq_varid(ncid, 'DU001' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%du1 = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 * dust_frac(1)  ![g/m3]
    q_gocart(i,j,k)%du2 = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 * dust_frac(2)  ![g/m3]
    q_gocart(i,j,k)%du3 = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 * dust_frac(3)  ![g/m3]
    q_gocart(i,j,k)%du4 = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 * dust_frac(4)  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DU002' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%du5 = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DU003' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%du6 = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DU004' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%du7 = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DU005' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%du8 = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo


 call check( nf90_inq_varid(ncid, 'SS001' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%ssa = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SS002' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%ssc = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SS003' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%ssc = q_gocart(i,j,k)%ssc + temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SS004' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%ssc = q_gocart(i,j,k)%ssc + temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo
!toshii --> only GEOS5 version of GOCART has SS mode 5. 
 call check( nf90_inq_varid(ncid, 'SS005' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%ssc = q_gocart(i,j,k)%ssc + temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo


 call check( nf90_inq_varid(ncid, 'SO4' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%so4 = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo


 call check( nf90_inq_varid(ncid, 'BCPHOBIC' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%blc = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'BCPHILIC' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%blc = q_gocart(i,j,k)%blc + temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'OCPHOBIC' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%ocn = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'OCPHILIC' , varid ) )      ![kg/kg]
 call check( nf90_get_var(ncid, varid, temp3d(is:ie,js:je,ks:ke), start=(/igs,jgs,kgs/), count=(/di,dj,dk/)))
 do k = ks,ke ; do j = js,je ; do i = is,ie
    kg = kge - k + 1  !reverse level
    q_gocart(i,j,k)%och = temp3d(i,j,kg) * rho(i,j,k) * 1.e3 ![g/m3]
 enddo ; enddo ; enddo

!
! close netCDF file
!
  call check( nf90_close(ncid) )


 endif



 return
 end subroutine rd_CRM_GEOS5_forecast_latlon
 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine get_aerosolfile_name( inp_name, out_name )
 implicit none
 character(len=*),intent(in)  :: inp_name
 character(len=*),intent(out) :: out_name

 integer :: i         ! looping
 integer :: file_len  ! total file length

!
! get file length
!
 file_len = len(trim(inp_name))  !length of character file name

!
! Find asm, and replace by aer
!  
 do i = 1, file_len-3
    if( inp_name(i:i+2) == 'asm' ) then
        out_name = inp_name(1:i-1)//'aer'//inp_name(i+3:file_len)
        exit 
    endif
 enddo

 return
 end subroutine get_aerosolfile_name

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine rd_CRM_SOUNDING
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Read Sounding data in in netCDF format for SARTA simualtor. 
! You must specify surface properties. 
!         
! History:
! 05/2014  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!---------------------------------------------------------------------------------------------------

 integer :: i,j,k,n  !loop
 integer :: ncid,ncidsbm          !# of variables to inquire  
 integer :: varid               !# of variable id
 integer :: ierr                !allocation stat
 integer :: is,ie,js,je,ks,ke,di,dj,dk
 integer :: nsample
 integer :: nskip ! skip number for sampling vertical level of sounding data. 
 integer :: k_sounding
 real(sdsu_fps),allocatable :: net1d(:)

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1

!
! memory allocation
!

!
! file name
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name)
 if(masterproc) print("(a)"),'MSG rd_CRM_SOUNDING: ++++++++++++++++++++ READING NEW FILE ++++++++++++++++++++ '
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)

!
! open netCDF file
!
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! derive sampling number 
!
 call check( nf90_inq_dimid(ncid, 'time',  varid) )
 call check( nf90_inquire_dimension( ncid, varid, len=nsample ) )

 nskip = nsample/mxlyr

 print*,'total sample',nsample,nsample/mxlyr

!
!
!
 if( .not. allocated(net1d) ) then
 allocate( &
           net1d     (nsample), &
           stat=ierr )
 if(ierr /= 0) call stop_sdsu(myrank,'MSG rd_CRM_SOUNDER: cannot allocate the variable')
 endif

!
! IMPORTANT NOTE 
! Sounding starts from 300m, so for simualtor purpose, we need surface anchor point at 0m altitude. 
!

!
! single point (always)
!
 i = 1 ; j = 1
!
! height and all sorts
!
 call check( nf90_inq_varid(ncid, 'alt' , varid ) ) ! geopotential [m/s2]  (staggered height) 
 call check( nf90_get_var(ncid, varid, net1d))

 atmos_stag(i,j,0)%hgt = 0.e0 !surface anchor point

 surface(i,j)%elev     = atmos_stag(i,j,k)%hgt
 do k = 1, mxlyr !convert
    k_sounding = k*nskip
    atmos_stag(i,j,k)%hgt = net1d(k_sounding)*1.e-3  !level altitude [km] <- [m]
 !   print*,k, atmos_stag(i,j,k)%hgt
 enddo
 do k = 1, mxlyr  !get layer value
    atmos(i,j,k)%dhgt = 0.5e0*(atmos_stag(i,j,k)%hgt-atmos_stag(i,j,k-1)%hgt) !thickness [km]
    atmos(i,j,k)%hgt  = 0.5e0*(atmos_stag(i,j,k)%hgt+atmos_stag(i,j,k-1)%hgt) !layer height [km]
    print*,k, atmos(i,j,k)%dhgt, atmos(i,j,k)%hgt
 enddo

!
! pressure 
!
 call check( nf90_inq_varid(ncid, 'pres' , varid ) ) ! pressure [m/s2]  (staggered height) 
 call check( nf90_get_var(ncid, varid, net1d))
 atmos_stag(i,j,0)%press = 1010.e0 ! surface pressure  (toshii: change later)
 do k = 1, mxlyr
    k_sounding = k*nskip
    atmos_stag(i,j,k)%press = net1d(k_sounding)  !pressure  [hPa] 
 enddo
 do k = 1, mxlyr  !get layer value
    atmos(i,j,k)%press  = 0.5e0*(atmos_stag(i,j,k)%press+atmos_stag(i,j,k-1)%press) !layer height [km]
 enddo

!
! temperature
!
 call check( nf90_inq_varid(ncid, 'tdry' , varid ) ) ! dry air temperature [degC]  (staggered height) 
 call check( nf90_get_var(ncid, varid, net1d))
 do k = 1, mxlyr
    k_sounding = k*nskip
    atmos_stag(i,j,k)%t_air = net1d(k_sounding) + 273.15 !air temp [K] <-- [degC]
    !print*,k, atmos_stag(i,j,k)%press
 enddo
 atmos_stag(i,j,0)%t_air = 325.387  ! surface temperature from WRF (toshii modify later,)
 surface(i,j)%t_skin = atmos_stag(i,j,0)%t_air      ! may get MODIS LST

 do k = 1, mxlyr  !get layer value
    atmos(i,j,k)%t_air  = 0.5e0*(atmos_stag(i,j,k)%t_air+atmos_stag(i,j,k-1)%t_air) !layer temp [K]
    !print*,k, atmos(i,j,k)%press
 enddo

!
! dry air density
!
 do k = 1, mxlyr
     atmos(i,j,k)%rho_dair = (atmos(i,j,k)%press*100.e0) / (const_Rd*atmos(i,j,k)%t_air)  * 1.e3 ! [g/m3]
 enddo

!
! RH
!
 call check( nf90_inq_varid(ncid, 'rh' , varid ) ) ! RH [%]  (staggered height) 
 call check( nf90_get_var(ncid, varid, net1d))
 do k = 1, mxlyr
    k_sounding = k*nskip
    atmos(i,j,k)%rh = net1d(k_sounding)  !RH [%]
 enddo


 call check( nf90_inq_varid(ncid, 'lat' , varid ) )  ! latitude [deg]
 call check( nf90_get_var(ncid, varid, surface%lat))

 call check( nf90_inq_varid(ncid, 'lon' , varid ) )  ! longitude [deg]
 call check( nf90_get_var(ncid, varid, surface%lon))

!
! close netCDF file
!
  call check( nf90_close(ncid) )

 return
 end subroutine rd_CRM_SOUNDING

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine rd_CRM_MLM
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
!    Read MLM output from Joe Santanello.
!  
! History:
! 04/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!
!----------------------------------------------------------------------------------------------------

 integer :: i,j,k,l  !looping
 integer :: ncid                !# of variables to inquire  
 integer :: varid               !# of variable id
 integer :: ierr                !allocation stat

 real(sdsu_fps),dimension(mxlyr) :: geohgt
 real(sdsu_fps),dimension(mxlyr) :: tpot
 real(sdsu_fps),dimension(mxlyr) :: p0
 real(sdsu_fps),dimension(mxlyr) :: q
 real(sdsu_fps) :: ts
 real(sdsu_fps) :: psurf

 real(sdsu_fps) :: esat,ewat,mix
!
! Open netCDF file
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name)
 if(masterproc) print("(a)"),'MSG rd_CRM_MLM: ++++++++++++++++++++ READING NEW FILE ++++++++++++++++++++'
 print("(a)"),'Input ->',trim(sdsu_io_file)

 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! height level and all sorts
!
 call check( nf90_inq_varid(ncid, 'Geoht' , varid ) )  ! geopotential height [m]
 call check( nf90_get_var(ncid, varid, geohgt ) )

 do k = 1,mxlyr - 1
    hgt_lev(k)       =    0.5e0*(geohgt(mxlyr-k+1) + geohgt(mxlyr-k))*1e-3   !1D height [km]
    atmos_stag(:,:,k)%hgt = 0.5e0*(geohgt(mxlyr-k+1) + geohgt(mxlyr-k))*1e-3  !3D height [km] 
 enddo

 hgt_lev(0)       =  0.  ! 1D height [km]
 atmos_stag(:,:,0)%hgt =  0.  ! 3D height [km] 

 hgt_lev(mxlyr)       = hgt_lev(mxlyr-1)       + hgt_lev(mxlyr-1)       - hgt_lev(mxlyr-2) !top height 1D [km]
 atmos_stag(:,:,mxlyr)%hgt = atmos_stag(:,:,mxlyr-1)%hgt + atmos_stag(:,:,mxlyr-1)%hgt &
                             - atmos_stag(:,:,mxlyr-2)%hgt !top height 3D [km]

!
! pressure and all sorts
!
 call check( nf90_inq_varid(ncid, 'p0' , varid ) )  !background pressure [mb] = [hPa]
 call check( nf90_get_var(ncid, varid, p0 ) )

 call check( nf90_inq_varid(ncid, 'Psurf' , varid ) )  !background pressure [mb] = [hPa]
 call check( nf90_get_var(ncid, varid, psurf ) )

 do k = 1,mxlyr
    atmos(:,:,k)%press = p0(mxlyr-k+1)
 enddo

 do k = 1, mxlyr-1
    atmos_stag(:,:,k)%press = 0.5e0* (p0(mxlyr-k+1) + p0(mxlyr-k)) ![hPa]
 enddo

 atmos_stag(:,:,0)%press = psurf  ! surface
 atmos_stag(:,:,mxlyr)%press = 0.5* (p0(1) + 0.)     ! top layer

!
! temperature and all sorts
!
 call check( nf90_inq_varid(ncid, 'T' , varid ) )  ! potential temperrature [K]
 call check( nf90_get_var(ncid, varid, tpot ) )

 call check( nf90_inq_varid(ncid, 'Ts' , varid ) )  ! potential temperrature [K]
 call check( nf90_get_var(ncid, varid, ts ) )

 do k = 1,mxlyr
    atmos(:,:,k)%t_air =   tpot(mxlyr-k+1)/( (1.e3/p0(mxlyr-k+1))**0.281 )
 enddo 

 do k = 1, mxlyr-1
    atmos_stag(:,:,k)%t_air = 0.5 * (atmos(:,:,k)%t_air + atmos(:,:,k+1)%t_air) ![K]
 enddo
 atmos_stag(:,:,0)%t_air     = atmos(:,:,1)%t_air      !surface level [K]
 atmos_stag(:,:,mxlyr)%t_air = atmos(:,:,mxlyr)%t_air  !top level [K]

 surface(:,:)%t_skin = ts


!
! relative humidity
!
 call check( nf90_inq_varid(ncid, 'q' , varid ) )  ! specific humidity [g/kg]
 call check( nf90_get_var(ncid, varid, q ) )


 do i = 1, mxgridx ; do j = 1, mxgridy ; do k = 1, mxlyr

    esat = 0.611 * exp (const_Lv_Rv* ( 1./const_Kel2Cel - 1./atmos(i,j,k)%t_air ) ) *10.e0  ![hPa]

    ewat =  p0(mxlyr-k+1) * q(mxlyr-k+1) / const_Rd_Rv / 1.e3 ! H2O partial pressure [hPa]

    atmos(i,j,k)%rh = ewat/esat * 100.  !RH [%]

 enddo ; enddo ; enddo


 atmos%rho_dair = 1.e+3*atmos%press*100.e0 / (const_Rd*atmos%t_air)  !dry air density [g/m3]


!
! Misc parameter not spcefified in MLM
!
  atmos_stag%w    = 1.
  surface%u10m = 1.

!
! close netCDF
!

 call check( nf90_close(ncid) )



 return
 end subroutine rd_CRM_MLM

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine rd_CRM_LIS
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Read WRF output (wrfout) in in netCDF format. All microphysics schemes are one-moment bulk.
!         
! History:
! 05/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
! WRF team, User's guide for Advanced Research WRF (ARW) modeling system version 2.2.
!            
!---------------------------------------------------------------------------------------------------

 integer :: i,j,k,n  !loop
 integer :: ncid,ncidsbm          !# of variables to inquire  
 integer :: varid               !# of variable id
 integer :: ierr                !allocation stat
 real(sdsu_fps)    :: press_top                              !model top pressure at interface [Pa] 
 real(sdsu_fps)    :: e, es     ! unsaturated and saturated water vapor pressure [kPa]
 integer :: itemp
 integer :: inet
 real(sdsu_fps) :: net
 integer :: is,ie,js,je,ks,ke,di,dj,dk
 real(sdsu_fps),allocatable :: rho(:,:,:)    !dry air density [kg/m3]
 real(sdsu_fps),allocatable :: net2d(:,:)
 real(sdsu_fps),allocatable :: net3d(:,:,:)
 real(sdsu_fps),allocatable :: net3d_stag(:,:,:)
 integer,allocatable :: inet2d(:,:)
 character(len=2) :: char_bin, tag_char
 character(len=50) :: para_char
 real(sdsu_fps) :: esat,ewat


!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1


!
! memory allocation
!
 if( .not. allocated(net2d) ) then
 allocate( &
           net2d     (is:ie,js:je          ), &
           inet2d    (is:ie,js:je          ), &
           net3d     (is:ie,js:je,1:mxlyr  ), &
           net3d_stag(is:ie,js:je,1:mxlyr+1), &
           rho       (is:ie,js:je,1:mxlyr  ), &
           stat=ierr )

 if(ierr /= 0) call stop_sdsu(myrank,'MSG rd_CRM_LIS: cannot allocate the variable')
 endif

!
! file name
!
 sdsu_io_file = trim(sdsu_dir_input)//trim(sdsu_inp_name)
 if(masterproc) print("(a)"),'MSG rd_CRM_LIS: ++++++++++++++++++++ READING NEW FILE ++++++++++++++++++++ '
#if MPI == 2
if(masterproc) &
#endif
 print("(a)"),'Input ->',trim(sdsu_io_file)

!
! open netCDF file
!
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )


!
! land mask
!
 call check( nf90_inq_varid(ncid, 'Landmask' , varid ) )  ! land-water mask 
 call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = myj_start, myj_end ; do i = myi_start, myi_end
    surface(i,j)%iland = INT( net2d(i,j) ) ! 1-land surface 2-ocean
 enddo ; enddo

!
! temperature and all sorts
!
 call check( nf90_inq_varid(ncid, 'AvgSurfT' , varid ) )  ! surface skin temperature [K]
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%t_skin, start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'Tair_f' , varid ) )    !air temperature [K]
 call check( nf90_get_var(ncid, varid, atmos(is:ie,js:je,ks:ke)%t_air, start=(/is,js,ks/), count=(/di,dj,dk/)))

 do j = myj_start, myj_end ; do i = myi_start, myi_end  !eoshii replace this by LIS input
    atmos_stag(i,j,0)%t_air = 0.5e0*(surface(i,j)%t_skin+atmos(i,j,1)%t_air) !interface <- average
 enddo ; enddo

 do k = 1, mxlyr-1 ; do j = myj_start, myj_end ; do i = myi_start, myi_end  
    atmos_stag(i,j,k)%t_air = 0.5e0 * ( atmos(i,j,k)%t_air + atmos(i,j,k+1)%t_air ) !interface <- average
 enddo ; enddo ; enddo

 do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos_stag(i,j,mxlyr)%t_air = atmos_stag(i,j,mxlyr-1)%t_air &
                       - (atmos_stag(i,j,mxlyr-1)%t_air-atmos(i,j,mxlyr)%t_air)*2.e0  !air temperature for top level
 enddo ; enddo 

!
! pressure and all sorts
!
 call check( nf90_inq_varid(ncid, 'Psurf_f' , varid ) ) ! Use base pressure [Pa]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))

 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
    atmos(i,j,k)%press = net3d(i,j,k) / 100.e0   ! [hPa] <- [Pa]
 enddo ; enddo ; enddo

! put staggred pressure here.....



!
! relative humidity
!
 call check( nf90_inq_varid(ncid, 'Qair_f' , varid ) ) ! specific humidity [kg/kg]
 call check( nf90_get_var(ncid, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))

 do i = 1, mxgridx ; do j = 1, mxgridy ; do k = 1, mxlyr
    esat = 0.611e0 * exp (const_Lv_Rv* ( 1.e0/const_Kel2Cel - 1.e0/atmos(i,j,k)%t_air ) ) *10.e0  ![hPa]
    ewat =  atmos(i,j,k)%press * (net3d(i,j,k)*1.e3) / const_Rd_Rv / 1.e3 ! H2O partial pressure [hPa]
    atmos(i,j,k)%rh = ewat/esat * 100.e0  !RH [%]
 enddo ; enddo ; enddo


  atmos%rho_dair = 1.e+3*atmos%press*100.e0 / (const_Rd*atmos%t_air)  !dry air density [g/m3]


!
! surface wind 
!
 call check( nf90_inq_varid(ncid, 'Wind_f' , varid ) ) !surface u wind [m/s]
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%u10m, start=(/is,js/), count=(/di,dj/)))


!
! surface precipitation 
!
 call check( nf90_inq_varid(ncid, 'Rainf_f' , varid ) ) ! precip [mm/s]
 call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = myj_start, myj_end ; do i = myi_start, myi_end
    surface(i,j)%rain_rate = net2d(i,j)  * 3600.e0   ! [mm/hr] <- [mm/s]
 enddo ; enddo

!
! surface parameters and all sorts
!
 call check( nf90_inq_varid(ncid, 'Landcover' , varid ) )  ! dominant vegetation type
 call check( nf90_get_var(ncid, varid, net2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = myj_start, myj_end ; do i = myi_start, myi_end
    surface(i,j)%igbp_typ =  INT(net2d(i,j))  
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SWE' , varid ) )  ! snow water equivalent [m]
 call check( nf90_get_var(ncid, varid, net2D(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 do j = myj_start, myj_end ; do i = myi_start, myi_end
    surface(i,j)%h2o_snow =  net2d(i,j) * 1.e3   ! [kg/m2] <- [m]
 enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SnowDepth' , varid ) )  ! snow depth [m]
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%dhgt_snow, start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'Greenness' , varid ) )  ! vegetation fraction [%]
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%frac_veg, start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'Albedo' , varid ) )  ! BROADBAND ALBEDO  
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%albedo, start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'SoilMoist' , varid ) )  ! soil moisture [m3 m-3]
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%h2o_soil, start=(/is,js,1/), count=(/di,dj,1/)))
        
 call check( nf90_inq_varid(ncid, 'SoilTemp' , varid ) )  ! soil temperature [K]
 call check( nf90_get_var(ncid, varid, surface(is:ie,js:je)%t_soil, start=(/is,js,1/), count=(/di,dj,1/)))




!
! close netCDF file
!
  call check( nf90_close(ncid) )


 return
 end subroutine rd_CRM_LIS

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine write_out_micro
  implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out microwave brightness temperature (V & H polarization) in GrADS format. 
! 
! History:
! 07/2011  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 10/2007  Toshi Matsui@NASA GSFC ; Microwave/radar/visir simulation account FOV
! 05/2007  Toshi Matsui@NASA GSFC ; adding downwelling tb for output. 
! 03/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer i,j,l,m,irec
 integer,parameter :: io = 103
 character(len=3) :: fff 
 character(len=6) :: tb_char

!
! compute Tb weighting averaged over FOV
!
 if( trim(sim_case) == 'MMF') then
     tb_out_fov = tb_out
 else
 do m = 1,2     ! vertical & horizontal loop
    do l = 1, mxfreq_micro ! frequency loop
        call fov2d(fov_ct_micro(l), fov_dt_micro(l), tb_out(1:mxgridx,1:mxgridy,l,m), &
                 tb_out_fov(1:mxgridx,1:mxgridy,l,m) )
    enddo
 enddo
 endif
!
! Computed Tbs in binary format
!
  sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                 trim(micro_sensor)//trim(output_suffix)//'.bin'
  print("(a)"),' ->',trim(sdsu_io_file)

  open(unit= io, file = sdsu_io_file, access = 'direct', &
       status = 'replace', recl = mxgridx*mxgridy*4)
  irec = 1

  do l=1,mxfreq_micro !frequency loop
     do m=1,2  !polarization loop
        write(io,rec=irec) ((tb_out_fov(i,j,l,m),i=1,mxgridx),j=1,mxgridy)  !microwave Tb
        irec = irec + 1
     enddo
  enddo
 
 close (io)

 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(micro_sensor)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.'//trim(micro_sensor)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(micro_sensor) //' microwave Tb'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxfreq_micro*2

    do l=1,mxfreq_micro !frequency loop
       do m=1,2  !polarization loop

          write(fff,"(I3.3)") INT(freq_micro(l))
          if(m==1)  then
              tb_char = 'tb'//fff//'h'
              write(io,*) tb_char,' 0  0  Tb '//trim(nch_micro(l))//'(H) [K]'
          elseif(m==2) then
              tb_char = 'tb'//fff//'v'
              write(io,*) tb_char,' 0  0  Tb '//trim(nch_micro(l))//'(V) [K]'
          endif

       enddo
    enddo
    write(io,*) 'endvars '
    close(io)
 endif


 return
 end subroutine write_out_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_radar
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out radar reflectivity [dBZ] in GrADS format. 
! 
! History:
! 03/2013  Toshi Matsui@NASA GSFC : MPI IO option.
! 07/2011  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 10/2007  Toshi Matsui@NASA GSFC ; FOV options
! 03/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,k,l  !for looping
 integer :: irec  !irec
 integer,parameter :: io = 104
 integer :: kmax
 character(len=3) :: fff
 integer :: is,ie,js,je,ks,ke
 integer :: nf
 real,allocatable :: out2d(:,:)

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr_radar
 kmax = mxlyr_radar
 if(allocated(out2d)) deallocate(out2d) ; allocate(out2d(is:ie,js:je))


!
! Computed radar echoes in binary format
!
  sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                trim(radar_sensor)//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
   print("(a)"),'MSG write_out_radar; dump out Zt, Zm, Vdop ',&
   ' ->',trim(sdsu_io_file)

!
! Write out binary data 
!
#if MPI == 2
 if(masterproc) &
#endif
  open(unit= io, file = sdsu_io_file, access='direct', &
       status='replace',recl=kmax*mxgridx*mxgridy*4)
  irec = 1 

  do nf = 1,mxfreq_radar
     call dump( radar_out_conv(is:ie,js:je,ks:ke,nf)%Zt    , io, irec )  !non-attenuated (true) Z [mm6/m3]
     call dump( radar_out_conv(is:ie,js:je,ks:ke,nf)%Zm    , io, irec )  !attenuated (measured) Z [mm6/m3]
     call dump( radar_out_conv(is:ie,js:je,ks:ke,nf)%Vdop  , io, irec )  !Doppler velocity [m/s]
  enddo ! l

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
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(radar_sensor)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//trim(radar_sensor)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(radar_sensor)//' reflectivity '
    write(io,*) 'undef  ',undefined 
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize  
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef   ',kmax  ,' linear 0. ',range_radar 
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxfreq_radar*3

    do l = 1, mxfreq_radar
       write(fff,"(I3.3)") INT(freq_radar(l))
       write(io,*) 'zt'//fff    ,kmax, ' 0  non-attenuated reflectivity '//trim(nch_radar(l))//' [mm6/m3]'
       write(io,*) 'zm'//fff    ,kmax, ' 0  attenuated reflectivty '//trim(nch_radar(l))//' [mm6/m3]'
       write(io,*) 'vdop'//fff  ,kmax, ' 0  Doppler velocity '//trim(nch_radar(l))//' [m/s]'
    enddo


    write(io,*) 'endvars '

    close(io)
 endif

#if MPI == 2
 endif
#endif



!
!  -------------------------     2D output     --------------------------------
!

  sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                trim(radar_sensor)//'_2D'//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
   print("(a)"),'MSG write_out_radar; dump out radar 2D ouptut  ',&
   ' ->',trim(sdsu_io_file)

!
! Write out binary data 
!
#if MPI == 2
 if(masterproc) &
#endif
  open(unit= io, file = sdsu_io_file, access='direct', &
       status='replace',recl=mxgridx*mxgridy*4)
  irec = 1

  do nf = 1,mxfreq_radar
     call dump( radar_out_2d_conv(is:ie,js:je,nf)%Het    , io, irec ) !echo-top height (AGL) [km]
     call dump( radar_out_2d_conv(is:ie,js:je,nf)%Zt_max , io, irec ) !Zt maxima [mm6/m3]
     call dump( radar_out_2d_conv(is:ie,js:je,nf)%Zm_max , io, irec ) !Zm maxima [mm6/m3]
     call dump( radar_out_2d_conv(is:ie,js:je,nf)%pia    , io, irec ) !PIA [dB]
  enddo ! l

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
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(radar_sensor)//'_2D'//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//trim(radar_sensor)//'_2D'//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(radar_sensor)//' reflectivity 2D composite'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef   1    linear 0. 1'
    write(io,*) 'tdef   1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxfreq_radar*4

    do l = 1, mxfreq_radar
       write(fff,"(I3.3)") INT(freq_radar(l))
       write(io,*) 'het'//fff  , '    0  0  echo-top height '//trim(nch_radar(l))//' [mm6/m3]'
       write(io,*) 'mxzt'//fff  ,'    0  0  maximum Zt '//trim(nch_radar(l))//' [mm6/m3]'
       write(io,*) 'mxzm'//fff  ,'    0  0  maximum Zm '//trim(nch_radar(l))//' [mm6/m3]'
       write(io,*) 'pia'//fff  , '    0  0  PIA '//trim(nch_radar(l))//' [dB]'
    enddo

    write(io,*) 'endvars '

    close(io)
 endif

#if MPI == 2
 endif
#endif





 return
 end subroutine write_out_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_visir
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out visible IR imager radiance (or brightness temperautre for wave > 10micon) in GrADS format. 
! 
! History:
! 07/2011  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 10/2007  Toshi Matsui@NASA GSFC ; FOV options
! 03/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,iw,nw,irec
 real(sdsu_fps) :: wave_number    !ending wavenumber [1/cm]
 integer :: iw_r 
 real(sdsu_fps) :: wl 
 integer,parameter :: io = 105
 character(len=2) :: ww

 if (visir_sensor == 'AIRS' ) then  !SPECIAL AIRS case

  sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(visir_sensor)//'.txt'
  print("(a)"),' ->',trim(sdsu_io_file)
  open(unit=io,file=sdsu_io_file,status='replace')
  write(io,*) 'id       [micron]      [cm-1]          [K] '

   nw = nw_airs

   do iw=1,nw !compute wavelength for every 1cm of wavenumberr
      iw_r = 1626 + iw
      wave_number = 10.0 ** (2.0 + REAL(iw_r)/2000.0)     ! [1/cm]
      wl = 1 / wave_number * 1.e+4                    ! wavelength [micron]
      write(io,78) iw, wl, 1/(wl*1e-4), radiance_out(1,1,iw)
 78 format(i4, 5x, f8.3, 5x, f8.2, 5x, f8.3)
   enddo

   close(io)
   
   return
 endif


!
! compute radiance/Tb wgt averaged over FOV
!

 do nw = 1, mxwavel
    call fov2d(fov_ct_visir(nw), fov_dt_visir(nw), radiance_out(1:mxgridx,1:mxgridy,nw), &
             radiance_out_fov(1:mxgridx,1:mxgridy,nw) )
 enddo


!
! write radiances in binary format
!
  sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                trim(visir_sensor)//trim(output_suffix)//'.bin'
  print("(a)"),' ->',trim(sdsu_io_file)
  open( unit= io, file = sdsu_io_file, access='direct', status='replace', recl=mxgridx*mxgridy*4)

  irec = 1 
  do nw = 1, mxwavel
     write(io,rec=irec) ((radiance_out_fov(i,j,nw),i=1,mxgridx), j=1,mxgridy)  ! Visible Radiance or IR Tb
     irec = irec + 1
  enddo

  close (io)


!
! write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(visir_sensor)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.'//trim(visir_sensor)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(visir_sensor) //' radiance/Tb'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxwavel

    do nw=1,mxwavel !frequency loop
       write(ww,"(I2.2)") nw
       if(wavel(nw) < tb_or_rad) then ! radiance output
          write(io,*) 'cha'//ww//' 0  0  radiance at '//trim(nch_wavel(nw))//' [W/m2/str/micron]'
       else ! Tb output
          write(io,*) 'cha'//ww//' 0  0  Tb at '//trim(nch_wavel(nw))//' [K]'
       endif
    enddo
    write(io,*) 'endvars '
    close(io)
  endif

  return
 end subroutine write_out_visir 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_visir_latlon
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out visible IR imager radiance (or brightness temperautre for wave > 10micon) in GrADS format. 
! 
! History:
! 03/2012  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,iw,nw,irec
 real(sdsu_fps) :: wave_number    !ending wavenumber [1/cm]
 integer :: iw_r 
 real(sdsu_fps) :: wl 
 integer,parameter :: io = 105
 character(len=2) :: ww

!
! compute radiance/Tb wgt averaged over FOV
!

 do nw = 1, mxwavel
    call fov2d(fov_ct_visir(nw), fov_dt_visir(nw), radiance_out(1:mxgridx,1:mxgridy,nw), &
             radiance_out_fov(1:mxgridx,1:mxgridy,nw) )
 enddo


!
! write radiances in binary format
!
  sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                trim(visir_sensor)//trim(output_suffix)//'latlon.bin'
  print("(a)"),' ->',trim(sdsu_io_file)
  open( unit= io, file = sdsu_io_file, access='direct', status='replace', recl=mxgridx*mxgridy*4)

  irec = 1 
  do nw = 1, mxwavel

    call dump2d_latlon( radiance_out(1:mxgridx,1:mxgridy,nw) ,&
                        lat_domain(1:mxgridx,1:mxgridy)      ,&
                        lon_domain(1:mxgridx,1:mxgridy)      ,&
                        io, irec )
  enddo

  close (io)


!
! write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(visir_sensor)//trim(output_suffix)//'latlon.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.'//trim(visir_sensor)//trim(output_suffix)//'latlon.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(visir_sensor) //' radiance/Tb'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',imax_latlon,' linear ',minlon, res_latlon
    write(io,*) 'ydef   ',jmax_latlon,' linear ',minlat, res_latlon
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxwavel

    do nw=1,mxwavel !frequency loop
       write(ww,"(I2.2)") nw
       if(wavel(nw) < tb_or_rad) then ! radiance output
          write(io,*) 'cha'//ww//' 0  0  radiance at '//trim(nch_wavel(nw))//' [W/m2/str/micron]'
       else ! Tb output
          write(io,*) 'cha'//ww//' 0  0  Tb at '//trim(nch_wavel(nw))//' [K]'
       endif
    enddo
    write(io,*) 'endvars '
    close(io)
  endif

  return
 end subroutine write_out_visir_latlon

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_lidar
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out lidar attenuating extinction coeficient in GrADS format. 
! 
! History:
! 07/2011  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 06/2011  Toshi Matsui@NASA GSFC : Options for either CRM or instrument level
! 10/2007  Toshi Matsui@NASA GSFC ; Initial 
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,k,nw
 integer :: irec
 integer :: kmax
 integer,parameter :: io = 106
 character(len=2) :: ww

 if( inst_profile_lidar ) then  !instrumental vertical cordinate
    kmax = mxlyr_lidar
 else
    kmax = mxlyr 
 endif
 
!
! write lidar attenuating backscatter in binary format
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(lidar_sensor)//trim(output_suffix)//'.bin'
 print("(a)"),' ->',trim(sdsu_io_file)

 open(unit= io, file = sdsu_io_file, access='direct', &
      status='replace',recl=kmax*mxgridx*mxgridy*4)

 irec = 1 
 do nw=1,mxwavel_lidar

    ! total attenuating backscattering [km-1 str-1]
    write(io,rec=irec) (((att_B      (i,j,k,nw),i=1,mxgridx),j=1,mxgridy),k=1,kmax) 
    irec = irec + 1

    ! perpendicular polarizationi component  [km-1 str-1]
    write(io,rec=irec) (((att_p22_B  (i,j,k,nw),i=1,mxgridx),j=1,mxgridy),k=1,kmax)
    irec = irec + 1

    ! total backscattering  [km-1 str-1]
    write(io,rec=irec) (((lidar_sback(i,j,k,nw),i=1,mxgridx),j=1,mxgridy),k=1,kmax) 
    irec = irec + 1

    ! lidar ratio
    write(io,rec=irec) (((lidar_ratio(i,j,k,nw),i=1,mxgridx),j=1,mxgridy),k=1,kmax) 
    irec = irec + 1 

 enddo


 close(io)


!
!  write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(lidar_sensor)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//trim(lidar_sensor)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(lidar_sensor)//' signal '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize

    if ( inst_profile_lidar ) then
       write(io,*) 'zdef   ',kmax  ,' linear 0. ',range_lidar
    else
       write(io,*) 'zdef   ',kmax  ,' linear 1. 1'
    endif

    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',irec-1

    do nw=1,mxwavel_lidar

       write(ww,"(I2.2)") nw
       write(io,*) 'att'//ww   ,kmax, ' 0  atten backscatter '//trim(nch_wavel_lidar(nw))//' [km-1 str-1]'
       write(io,*) 'ptt'//ww   ,kmax, ' 0  p22 atten backscatter '//trim(nch_wavel_lidar(nw))//' [km-1 str-1]'
       write(io,*) 'ntt'//ww   ,kmax, ' 0  non-atten backscatter '//trim(nch_wavel_lidar(nw))//'  [km-1 str-1]'
       write(io,*) 'ratio'//ww ,kmax, ' 0  lidar ratio '//trim(nch_wavel_lidar(nw))//'  [-]'
 
    enddo


    write(io,*) 'endvars '
    close(io)
 endif


 return
 end subroutine write_out_lidar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_broad
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out broadband energy budget in GrADS format. 
! 
! History:
! 07/2011  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 12/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,k,n
 integer,parameter :: io = 107

!
! FOV convolution
!
   call fov3d(fov_ct_broad, fov_dt_broad, ebudget(:,:,:), &
             ebudget_fov(:,:,:) )

!
! Energy budget in binary format
!
  sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.BROAD'//trim(output_suffix)//'.bin'
  print("(a)"),' ->',trim(sdsu_io_file)
 
  open( unit= io, file = sdsu_io_file, access='direct', status = 'replace', recl=mxgridx*mxgridy*4)

  do n = 1, 14
     write(io,rec=n) ((ebudget_fov(i,j,n) ,i=1,mxgridx), j=1,mxgridy)  !energy budget
  enddo 
  write(io,rec=15) ((cod_broad(i,j) ,i=1,mxgridx), j=1,mxgridy)  !cloud optical depths
  write(io,rec=16) ((aod_broad(i,j) ,i=1,mxgridx), j=1,mxgridy)  !aerosol optical depths
  close (io)


!
! write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.BROAD'//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.BROAD'//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  Energy budgets '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',14+2
    write(io,*) 'swdt   0 0  TOA SW downwelling flux [W/m2]'
    write(io,*) 'swut   0 0  TOA SW upwelling flux   [W/m2]'
    write(io,*) 'swds   0 0  Surface SW downwelling flux [W/m2]'
    write(io,*) 'swus   0 0  Surface SW upwelling   flux [W/m2]'
    write(io,*) 'swnt   0 0  Net toa        SW radiation [W/m2]'
    write(io,*) 'swns   0 0  Net surface    SW radiation [W/m2]'
    write(io,*) 'swna   0 0  Net atmosphere SW radiation [W/m2]'
    write(io,*) 'lwdt   0 0  TOA LW downwelling flux [W/m2]'
    write(io,*) 'lwut   0 0  TOA LW upwelling   flux [W/m2]'
    write(io,*) 'lwds   0 0  Surface LW downwelling flux [W/m2]'
    write(io,*) 'lwus   0 0  Surface LW upwelling   flux [W/m2]'
    write(io,*) 'lwnt   0 0  Net        toa LW radiation [W/m2]'
    write(io,*) 'lwns   0 0  Net    surface LW radiation [W/m2]'
    write(io,*) 'lwna   0 0  Net atmosphere LW radiation [W/m2]'
    write(io,*) 'cod    0 0  cloud optical depth (vis) [-]'
    write(io,*) 'aod    0 0  aerosol optical depth [-]'
    write(io,*) 'endvars '
    close(io)
 endif

!
!  Heating rate
!
  if(heating_rate) then

      sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.HEAT'//trim(output_suffix)//'.bin'
      open( unit= io, file = sdsu_io_file, access='direct', status = 'replace', recl=mxgridx*mxgridy*mxlyr*4)
      write(io,rec=1) (((sw_heat(i,j,k), i=1,mxgridx), j=1,mxgridy), k=1,mxlyr)  ! SW Heating Rate [K/day]
      write(io,rec=2) (((lw_heat(i,j,k), i=1,mxgridx), j=1,mxgridy), k=1,mxlyr)  ! LW Heating Rate [K/day]
      close (io)

    if( write_grads_ctl ) then
       sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.HEAT'//trim(output_suffix)//'.ctl'
       print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
       open(io,file=trim(sdsu_io_file))
       sdsu_io_file = sdsu_inp_name(1:efile_len)//'.HEAT'//trim(output_suffix)//'.bin'

       write(io,"(a)") 'dset  '//trim(sdsu_io_file)
       write(io,*) 'title  SW/LW heating rate '
       write(io,*) 'undef  ',undefined
       write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
       write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
       write(io,*) 'zdef   ',mxlyr,  ' linear 0. 1 '
       write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
       write(io,*) 'vars ',2
       write(io,*) 'swhr ',mxlyr,  ' 0  shortwave heating rate [K/day]'
       write(io,*) 'lwhr ',mxlyr,  ' 0  longwave heating rate  [K/day]'
       write(io,*) 'endvars '
       close(io)
    endif

  endif

  return
 end subroutine write_out_broad

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_isccp
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out ISCCP-product consistent (tau and top pressure) output in GrADS format. 
! 
! History:
! 07/2011  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 12/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j
 integer,parameter :: io = 108

!
! Write out radiances in binary format
!
  sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.ISCCP'//trim(output_suffix)//'.bin'
  print("(a)"),' ->',trim(sdsu_io_file)

  open( unit= io, file = sdsu_io_file, access='direct', status = 'replace', recl=mxgridx*mxgridy*4)
  write(io,rec=1) ((cld_opt(i,j) ,i=1,mxgridx), j=1,mxgridy)     ! column cloud optical depth [-]
  write(io,rec=2) ((cldtop_p(i,j)   ,i=1,mxgridx), j=1,mxgridy)  ! cloud-top pressure [mb]

!     write(22,rec=3) ((aerosol_opt(i,j) ,i=1,mxgridx), j=1,mxgridy)  ! column aerosol optical depth [-]
!     write(22,rec=4) ((aerotop_p(i,j)   ,i=1,mxgridx), j=1,mxgridy)  ! aerosol-top pressure [mb]

  close (io)


 return
 end subroutine write_out_isccp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_isccp_color
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out ISCCP color diagrams in GrADS format. 
! 
! History:
! 12/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: i,j,ind
 integer,parameter :: imax=599, jmax=699
 real(sdsu_fps),dimension(0:imax,0:jmax) :: tau, pres
 real(sdsu_fps) :: w1,w2

!
! ISCCP Color diagrams for Grads
!

  do i = 0,imax
    ind = i/100+1
    w1  = REAL(REAL(ind)*100.-REAL(i))/100.
    w2  = 1.-w1

    if( ind == 1) tau(i,:) = isccp_tau(ind)*w1 + isccp_tau(ind+1)*w2
    if( ind == 2) tau(i,:) = isccp_tau(ind)*w1 + isccp_tau(ind+1)*w2 
    if( ind == 3) tau(i,:) = isccp_tau(ind)*w1 + isccp_tau(ind+1)*w2
    if( ind == 4) tau(i,:) = isccp_tau(ind)*w1 + isccp_tau(ind+1)*w2
    if( ind == 5) tau(i,:) = isccp_tau(ind)*w1 + isccp_tau(ind+1)*w2
    if( ind == 6) tau(i,:) = isccp_tau(ind)*w1 + isccp_tau(ind+1)*w2
  enddo

  do j = 0,jmax
    ind = j/100+1
    w1  = REAL(REAL(ind)*100.-REAL(j))/100.
    w2  = 1.-w1

    if( ind == 1) pres(:,j) = isccp_pre(ind)*w1 + isccp_pre(ind+1)*w2
    if( ind == 2) pres(:,j) = isccp_pre(ind)*w1 + isccp_pre(ind+1)*w2
    if( ind == 3) pres(:,j) = isccp_pre(ind)*w1 + isccp_pre(ind+1)*w2
    if( ind == 4) pres(:,j) = isccp_pre(ind)*w1 + isccp_pre(ind+1)*w2
    if( ind == 5) pres(:,j) = isccp_pre(ind)*w1 + isccp_pre(ind+1)*w2
    if( ind == 6) pres(:,j) = isccp_pre(ind)*w1 + isccp_pre(ind+1)*w2
    if( ind == 7) pres(:,j) = isccp_pre(ind)*w1 + isccp_pre(ind+1)*w2
  enddo

!
! Write out ISCCP color diagrams for Grads
!
  sdsu_io_file =trim(sdsu_dir_output)//'ISCCP_color.bin'
  print("(a)"),'ISCCP_color diagram for GrADS ->',trim(sdsu_io_file)

  open( unit= 22, file = sdsu_io_file, access='direct', status = 'replace', recl=(1+imax)*(1+jmax)*4)
  write(22,rec=1) ((tau (i,j) ,i=0,imax), j=0,jmax)   
  write(22,rec=2) ((pres(i,j) ,i=0,imax), j=0,jmax)   
  close (22)

  return
 end subroutine write_out_isccp_color

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine fov3d(fov_ct,fov_dt,var,var_fov)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute FOV convolution via Gaussian beam filling.  
! 
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
! Masunaga, H., and C.D. Kummerow, 2005: Combined Radar and Radiometer Analysis of 
!      Precipitation Profiles for a Parametric Retrieval Algorithm. J. Atmos. Oceanic 
!      Technol., 22, 909-929.
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: fov_ct, fov_dt       ! FOV size in cross track and downtrack [km]
 real(sdsu_fps),intent(in) :: var(:,:,:)  ! input variables (-999 is no data)
 real(sdsu_fps),intent(out) :: var_fov(:,:,:)  ! FOV variables

 integer :: n_ct, n_dt           ! grid # for cross/down-tracking path
 integer :: ii,jj,i,j,i_dt,j_ct  ! looping
 real(sdsu_fps) :: sum_wgt
 real(sdsu_fps),allocatable :: sum_var(:) ! sum of weight and variables
 real(sdsu_fps) :: p,wgt           ! weight and its function
 real(sdsu_fps) :: xa,ya           ! physical distance in x and y direction
 integer :: bnd(3)  ! 3D array size


 if(mxgridx <= 10 .or. mxgridy <= 10 ) then !test
    var_fov = var
    return
 endif

 bnd = UBOUND(var)
 allocate( sum_var(bnd(3)) )

  n_ct = nint(fov_ct/(2.*gridsize)) + 2
  n_dt = nint(fov_dt/(2.*gridsize)) + 2

  if( fov_ct > REAL(mxgridx)*gridsize .or. fov_dt > REAL(mxgridx)*gridsize   ) then
      n_dt = 0 ; n_ct = 0
  endif

  if( fov_ct > REAL(mxgridy)*gridsize .or. fov_dt > REAL(mxgridy)*gridsize   ) then
      n_dt = 0 ; n_ct = 0
  endif

!
! temporal for TRMM
!
  if(mxgridx == 1 .or. mxgridy == 1) then
        n_dt = 0 ; n_ct = 0
  endif


  do j = 1, mxgridy
     do i = 1, mxgridx
        sum_wgt = 0
        sum_var = 0
        do ii = i - n_dt, i + n_dt
!
!          code uses cyclic boundaries in order to fill the FOV
!          near the edge of the scan.  
!
           i_dt = ii
           if (i_dt < 1) i_dt = i_dt + mxgridx
           if (i_dt > mxgridx) i_dt = i_dt - mxgridx
           do jj = j - n_ct, j + n_ct
              j_ct = jj
              if (j_ct < 1) j_ct = j_ct + mxgridy
              if (j_ct > mxgridy) j_ct = j_ct - mxgridy
              ! Compute pixel wgt based on distance from center point
              xa = gridsize*(ii-i)  ![km]
              ya = gridsize*(jj-j)  ![km]
              p = xa*xa/((fov_ct**2)/4.0) + ya*ya/((fov_dt**2)/4.0)
              wgt = exp(-0.5*(1.1774*p)**2)
              if(var(i_dt,j_ct,1) /= undefined ) then
                 sum_wgt = sum_wgt + wgt
                 sum_var(:) = sum_var(:) + wgt*var(i_dt,j_ct,:)
              endif

           enddo !jj
        enddo !ii
        if(sum_wgt > 0.e0) then
           var_fov(i,j,:) = sum_var(:)/sum_wgt
        else
           var_fov(i,j,:) = undefined
        endif
     enddo !i
  enddo !j

 deallocate(sum_var)

 return
 end subroutine fov3d


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine fov3d_radar(kmax, fov_ct,fov_dt,var3d,var3d_fov)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute FOV convolution via Gaussian beam filling.  
! 
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
! Masunaga, H., and C.D. Kummerow, 2005: Combined Radar and Radiometer Analysis of 
!      Precipitation Profiles for a Parametric Retrieval Algorithm. J. Atmos. Oceanic 
!      Technol., 22, 909-929.
!----------------------------------------------------------------------------------------------------
 integer,intent(in) :: kmax
 real(sdsu_fps),intent(in) :: fov_ct, fov_dt       ! FOV size in cross track and downtrack [km]
 real(sdsu_fps),intent(in) :: var3d(1:,1:,1:)  ! input variables (-999 is no data)
 real(sdsu_fps),intent(out) :: var3d_fov(1:,1:,1:)  ! FOV variables

 integer :: bnd(3)  ! 3D array size
 integer :: ics,ice,jcs,jce
 integer :: n_ct, n_dt           ! grid # for cross/down-tracking path
 integer :: k, ii,jj,i,j,i_dt,j_ct  ! looping
 real(sdsu_fps) :: p,wgt,mxwgt     ! weight and its function
 real(sdsu_fps) :: xa,ya           ! physical distance in x and y direction
 real(sdsu_fps),allocatable :: sum_var(:),sum_wgt(:) ! sum of weight and variables

!
! get dimension
!
 bnd = UBOUND(var3d)

 if(allocated(sum_var)) deallocate(sum_var,sum_wgt)
 allocate( sum_var(1:bnd(3)), sum_wgt(1:bnd(3)) ) 

  n_ct = nint(fov_ct/(2.*gridsize)) 
  n_dt = nint(fov_dt/(2.*gridsize)) 

  if( fov_ct > REAL(mxgridx)*gridsize .or. fov_dt > REAL(mxgridx)*gridsize   ) then
      n_dt = 0 ; n_ct = 0
  endif

  if( fov_ct > REAL(mxgridy)*gridsize .or. fov_dt > REAL(mxgridy)*gridsize   ) then
      n_dt = 0 ; n_ct = 0
  endif



 do j = 1, bnd(2) ; do i = 1, bnd(1)

    sum_wgt = 0.
    sum_var = 0.

    !
    ! convolution loop bounds
    !
    ics = max(1, i - n_dt) ; ice = min( bnd(1), i + n_dt )
    jcs = max(1, j - n_ct) ; jce = min( bnd(2), j + n_ct )

    sub_loop: do jj = jcs, jce ; do ii = ics, ice

       !
       ! Compute pixel wgt based on distance from center point
       !
       xa = gridsize* (ii-i)  ![km]
       ya = gridsize* (jj-j)  ![km]
       p = xa*xa/((fov_ct**2)/4.0) + ya*ya/((fov_dt**2)/4.0)
       wgt = exp(-0.5*(1.1774*p)**2)

       do k = 1, bnd(3)
              
          if(var3d(ii,jj,k) /= undefined ) then

             sum_wgt(k) = sum_wgt(k) + wgt
             sum_var(k) = sum_var(k) + wgt*var3d(ii,jj,k)

          endif

       enddo !k

       !
       ! if one of parameter is undefine convoluted FOV value must be undefined and exit. 
       !

    enddo ; enddo sub_loop

    mxwgt = maxval( sum_wgt(:) )
    do k = 1, bnd(3)
       if(sum_wgt(k) >= mxwgt*0.95 ) then
          var3d_fov(i,j,k) = sum_var(k)/sum_wgt(k)
       else
          var3d_fov(i,j,k) = undefined
       endif
    enddo

  enddo ; enddo 


 return
 end subroutine fov3d_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine fov2d(fov_ct,fov_dt,var,var_fov)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute FOV convolution via Gaussian beam filling.  
! 
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
! Masunaga, H., and C.D. Kummerow, 2005: Combined Radar and Radiometer Analysis of 
!      Precipitation Profiles for a Parametric Retrieval Algorithm. J. Atmos. Oceanic 
!      Technol., 22, 909-929.
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: fov_ct, fov_dt       ! FOV size in cross track and downtrack [km]
 real(sdsu_fps),intent(in) :: var(mxgridx,mxgridy)  ! input variables (-999 is no data)
 real(sdsu_fps),intent(out) :: var_fov(mxgridx,mxgridy)  ! FOV variables

 integer :: n_ct, n_dt           ! grid # for cross/down-tracking path
 integer :: ii,jj,i,j,i_dt,j_ct  ! looping
 real(sdsu_fps) :: sum_wgt,sum_var ! sum of weight and variables
 real(sdsu_fps) :: p,wgt           ! weight and its function
 real(sdsu_fps) :: xa,ya           ! physical distance in x and y direction

 if(mxgridx <= 10 .or. mxgridy <= 10 ) then !test
    var_fov = var
    return
 endif

  n_ct = nint(fov_ct/(2.*gridsize)) + 2
  n_dt = nint(fov_dt/(2.*gridsize)) + 2


  if( fov_ct > REAL(mxgridx)*gridsize .or. fov_dt > REAL(mxgridx)*gridsize   ) then
      n_dt = 0 ; n_ct = 0
  endif

  if( fov_ct > REAL(mxgridy)*gridsize .or. fov_dt > REAL(mxgridy)*gridsize   ) then
      n_dt = 0 ; n_ct = 0
  endif

!
! temporal for TRMM
!
  if(mxgridx == 1 .or. mxgridy == 1) then 
        n_dt = 0 ; n_ct = 0
  endif


  do j = 1, mxgridy
     do i = 1, mxgridx
        sum_wgt = 0
        sum_var = 0
        do ii = i - n_dt, i + n_dt
!
!          code uses cyclic boundaries in order to fill the FOV
!          near the edge of the scan.  
!
           i_dt = ii
           if (i_dt < 1) i_dt = i_dt + mxgridx
           if (i_dt > mxgridx) i_dt = i_dt - mxgridx
           do jj = j - n_ct, j + n_ct
              j_ct = jj
              if (j_ct < 1) j_ct = j_ct + mxgridy
              if (j_ct > mxgridy) j_ct = j_ct - mxgridy
              ! Compute pixel wgt based on distance from center point
              xa = gridsize*(ii-i)  ![km]
              ya = gridsize*(jj-j)  ![km]
              p = xa*xa/((fov_ct**2)/4.0) + ya*ya/((fov_dt**2)/4.0)
              wgt = exp(-0.5*(1.1774*p)**2)
              if(var(i_dt,j_ct) /= undefined ) then
                 sum_wgt = sum_wgt + wgt
                 sum_var = sum_var + wgt*var(i_dt,j_ct)
              endif

           enddo !jj
        enddo !ii
        if(sum_wgt > 0.) then
           var_fov(i,j) = sum_var/sum_wgt
        else
           var_fov(i,j) = undefined
        endif
     enddo !i
  enddo !j

 return
 end subroutine fov2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_CRM2D
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out CRM 2D input (surface data) in GrADS format. 
! 
! History:
! 07/2011  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: i,j,k,irec
 character(len=20) :: typ
 integer,parameter :: io = 109
 integer :: is,ie,js,je,ks,ke

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end

!
! write out CRM surface input for grads binarry format
!
  typ = '.CRM2D'

!
! output path + file name
!
 sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//trim(typ)//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif 
 print("(a)"),'MSG write_out_CRM2D: write CRM 2D input file in grads format',&
 ' ->',trim(sdsu_io_file)

!
! open output file
!
#if MPI == 2
 if(masterproc) &
#endif
 open(unit= io, file = sdsu_io_file, access='direct', status='replace',recl = mxgridx*mxgridy*4 )
 irec = 1
!
! write column-integrated condensate amount [kg/m2]
!
 mic_select1: select case(trim(type_microphysics))
 case('GEN')
  call dump( qcol_gen(is:ie,js:je)%cloud  , io, irec ) ! cloud water 
  call dump( qcol_gen(is:ie,js:je)%rain   , io, irec ) ! rain 
  call dump( qcol_gen(is:ie,js:je)%ice    , io, irec ) ! cloud ice 
  call dump( qcol_gen(is:ie,js:je)%snow   , io, irec ) ! snow 
  call dump( qcol_gen(is:ie,js:je)%graupel, io, irec ) ! graupel 
  call dump( qcol_gen(is:ie,js:je)%hail   , io, irec ) ! hail 

 case('RAMS')
  call dump( qcol_rams(is:ie,js:je)%cloud1 , io, irec ) ! cloud mode1 (small)
  call dump( qcol_rams(is:ie,js:je)%cloud2 , io, irec ) ! cloud mode2 (large)
  call dump( qcol_rams(is:ie,js:je)%rain   , io, irec ) ! rain 
  call dump( qcol_rams(is:ie,js:je)%ice1   , io, irec ) ! ice mode1 (small) 
  call dump( qcol_rams(is:ie,js:je)%ice2   , io, irec ) ! ice mode2 (large) 
  call dump( qcol_rams(is:ie,js:je)%snow   , io, irec ) ! snow aggregate 
  call dump( qcol_rams(is:ie,js:je)%graupel, io, irec ) ! graupel 
  call dump( qcol_rams(is:ie,js:je)%hail   , io, irec ) ! hail 
 case('SBM')
  call dump( qcol_sbm(is:ie,js:je)%liq    , io, irec ) ! liquid (cloud+rain) 
  call dump( qcol_sbm(is:ie,js:je)%ice_col, io, irec ) ! ice column 
  call dump( qcol_sbm(is:ie,js:je)%ice_pla, io, irec ) ! ice plate 
  call dump( qcol_sbm(is:ie,js:je)%ice_den, io, irec ) ! ice dendride 
  call dump( qcol_sbm(is:ie,js:je)%snow   , io, irec ) ! snow 
  call dump( qcol_sbm(is:ie,js:je)%graupel, io, irec ) ! graupel 
  call dump( qcol_sbm(is:ie,js:je)%hail   , io, irec ) ! graupel 
 end select mic_select1

!
! write land-surface properties
!
  call dump( surface(is:ie,js:je)%h2o_col       , io, irec )  ! column water vapor [kg/m2]
  call dump( surface(is:ie,js:je)%lat           , io, irec )  ! latitude [deg]
  call dump( surface(is:ie,js:je)%lon           , io, irec )  ! longitude [deg]
  call dump( REAL(surface(is:ie,js:je)%igbp_typ), io, irec )  ! IGBP vegetation type 
  call dump( surface(is:ie,js:je)%h2o_soil      , io, irec )  ! soil moisture frac [m3 m-3]
  call dump( surface(is:ie,js:je)%dhgt_snow     , io, irec )  ! snow depth [m]
  call dump( surface(is:ie,js:je)%frac_veg      , io, irec )  ! vegetation fraction  [%]
  call dump( surface(is:ie,js:je)%t_skin        , io, irec )  ! skin temperature [K]
  call dump( surface(is:ie,js:je)%rain_rate     , io, irec )  ! surface rainfall [mm/hr]
  call dump( surface(is:ie,js:je)%elev          , io, irec )  ! surface elevation [m]
  call dump( surface(is:ie,js:je)%cape          , io, irec )  ! CAPE [J/kg]
  call dump( surface(is:ie,js:je)%ncape         , io, irec )  ! NCAPE [J/kg/m]
  call dump( surface(is:ie,js:je)%lcl           , io, irec )  ! LCL [m]
  call dump( surface(is:ie,js:je)%wcd           , io, irec )  ! WCD [m]
  call dump( REAL(surface(is:ie,js:je)%iclass)  , io, irec )  ! cloud class [-]


#if MPI == 2
 if(masterproc) &
#endif
 close (io)

 if(masterproc) print*,''


!
!  write out grads control file
!
#if MPI == 2
 if(masterproc) then
#endif

 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//trim(typ)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//trim(typ)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  CRM 2D file for diagnostic'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',irec-1

 mic_select2: select case(trim(type_microphysics))
 case('GEN')
    write(io,*) 'qc     0  0  cloud water  [kg/m2]'
    write(io,*) 'qr     0  0  rain water  [kg/m2]'
    write(io,*) 'qi     0  0  ice [kg/m2]'
    write(io,*) 'qs     0  0  snow  [kg/m2]'
    write(io,*) 'qg     0  0  graupel  [kg/m2]'
    write(io,*) 'qh     0  0  hail  [kg/m2]'
 case('RAMS')
    write(io,*) 'qc1     0  0  cloud mode1 [kg/m2]'
    write(io,*) 'qc2     0  0  cloud mode2 [kg/m2]'
    write(io,*) 'qr      0  0  rain water  [kg/m2]'
    write(io,*) 'qi1     0  0  ice mode1 [kg/m2]'
    write(io,*) 'qi2     0  0  ice mode2 [kg/m2]'
    write(io,*) 'qs      0  0  snow  [kg/m2]'
    write(io,*) 'qg      0  0  graupel  [kg/m2]'
    write(io,*) 'qh      0  0  hail  [kg/m2]'
 case('SBM')
    write(io,*) 'ql     0  0  liquid water path  [kg/m2]'
    write(io,*) 'qic    0  0  ice column path [kg/m2]'
    write(io,*) 'qip    0  0  ice plate path [kg/m2]'
    write(io,*) 'qid    0  0  ice dendrate path [kg/m2]'
    write(io,*) 'qs     0  0  snow path [kg/m2]'
    write(io,*) 'qg     0  0  graupel path [kg/m2]'
    write(io,*) 'qh     0  0  hail path [kg/m2]'
 end select mic_select2

 
    write(io,*) 'h2ocol    0  0  column water vapro [kg/m2]'
    write(io,*) 'lati      0  0  latitude  [deg]'
    write(io,*) 'long      0  0  longitude [deg]'
    write(io,*) 'igbp      0  0  IGBP vegetation type'
    write(io,*) 'soilm     0  0  soil moisture frac [m3 m-3]'
    write(io,*) 'snowd     0  0  snow depth [m]'
    write(io,*) 'fveg      0  0  vegetation fraction  [%]'
    write(io,*) 'tsk       0  0  skin temperature [K]'
    write(io,*) 'srain     0  0  surface rainfall [mm/hr]'
    write(io,*) 'elev      0  0  surface elevation [m]'
    write(io,*) 'cape      0  0  CAPE [J/kg]'
    write(io,*) 'ncape     0  0  NCAPE [J/kg/m]'
    write(io,*) 'lcl       0  0  LCL [m]'
    write(io,*) 'wcd       0  0  WCD [m]'
    write(io,*) 'cldtyp    0  0  1-shallow 2-mid, 3-high, 4 convective'

    write(io,*) 'endvars '
    close(io)
 endif

#if MPI == 2
 endif
#endif

 return
 end subroutine write_out_CRM2D

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_CRM3D
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out CRM 3D input (profile) in GrADS format. 
! 
! History:
! 07/2011  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: ierr
 integer :: i,j,k,irec
 character(len=20) :: typ
 integer,parameter :: io = 110
 integer :: is,ie,js,je,ks,ke

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end

!
! write out CRM surface input for grads binarry format
!
  typ = '.CRM3D'

!
! write out CRM input for grads binarry format
!
  sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//trim(typ)//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
  print("(a)"),'MSG write_out_CRM3D: write CRM 3D input file in grads format ',&
  ' ->',trim(sdsu_io_file)

!
! open CRM3D output file
!
#if MPI == 2
 if(masterproc) &
#endif
  open(unit= io, file = sdsu_io_file, access='direct', status='replace', recl=mxgridx*mxgridy*mxlyr*4 )
  irec = 1

!
! write atmospheric parameter
!
  call dump( atmos(is:ie,js:je,ks:ke)%t_air    , io, irec )     
  call dump( atmos(is:ie,js:je,ks:ke)%w        , io, irec )
  call dump( atmos(is:ie,js:je,ks:ke)%rh       , io, irec )
  call dump( atmos(is:ie,js:je,ks:ke)%dhgt*1.e3, io, irec )
  call dump( atmos(is:ie,js:je,ks:ke)%hgt *1.e3, io, irec )
  call dump( atmos(is:ie,js:je,ks:ke)%rain_rate , io, irec )  !3D rain rate

!
! write condensate parameter
! 
 mic_select1: select case(trim(type_microphysics))
 case('GEN') 
  call dump( q_gen(is:ie,js:je,ks:ke)%cloud  , io, irec )
  call dump( q_gen(is:ie,js:je,ks:ke)%rain   , io, irec ) 
  call dump( q_gen(is:ie,js:je,ks:ke)%ice    , io, irec )
  call dump( q_gen(is:ie,js:je,ks:ke)%snow   , io, irec )
  call dump( q_gen(is:ie,js:je,ks:ke)%graupel, io, irec )
  call dump( q_gen(is:ie,js:je,ks:ke)%hail   , io, irec )
  call dump( re_gen(is:ie,js:je,ks:ke)%cloud  , io, irec )
  call dump( re_gen(is:ie,js:je,ks:ke)%rain   , io, irec )
  call dump( re_gen(is:ie,js:je,ks:ke)%ice    , io, irec )
  call dump( re_gen(is:ie,js:je,ks:ke)%snow   , io, irec )
  call dump( re_gen(is:ie,js:je,ks:ke)%graupel, io, irec )
  call dump( re_gen(is:ie,js:je,ks:ke)%hail   , io, irec )
 case('RAMS')
  call dump( q_rams(is:ie,js:je,ks:ke)%cloud1  , io, irec )
  call dump( q_rams(is:ie,js:je,ks:ke)%cloud2  , io, irec )
  call dump( q_rams(is:ie,js:je,ks:ke)%rain    , io, irec )
  call dump( q_rams(is:ie,js:je,ks:ke)%ice1    , io, irec )
  call dump( q_rams(is:ie,js:je,ks:ke)%ice2    , io, irec )
  call dump( q_rams(is:ie,js:je,ks:ke)%snow    , io, irec )
  call dump( q_rams(is:ie,js:je,ks:ke)%graupel , io, irec )
  call dump( q_rams(is:ie,js:je,ks:ke)%hail    , io, irec )
  call dump( re_rams(is:ie,js:je,ks:ke)%cloud1  , io, irec )
  call dump( re_rams(is:ie,js:je,ks:ke)%cloud2  , io, irec )
  call dump( re_rams(is:ie,js:je,ks:ke)%rain    , io, irec )
  call dump( re_rams(is:ie,js:je,ks:ke)%ice1    , io, irec )
  call dump( re_rams(is:ie,js:je,ks:ke)%ice2    , io, irec )
  call dump( re_rams(is:ie,js:je,ks:ke)%snow    , io, irec )
  call dump( re_rams(is:ie,js:je,ks:ke)%graupel , io, irec )
  call dump( re_rams(is:ie,js:je,ks:ke)%hail    , io, irec )
 case('SBM')
  call dump( q_sbm(is:ie,js:je,ks:ke)%liq    , io, irec )
  call dump( q_sbm(is:ie,js:je,ks:ke)%ice_col, io, irec )
  call dump( q_sbm(is:ie,js:je,ks:ke)%ice_pla, io, irec )
  call dump( q_sbm(is:ie,js:je,ks:ke)%ice_den, io, irec )
  call dump( q_sbm(is:ie,js:je,ks:ke)%snow   , io, irec )
  call dump( q_sbm(is:ie,js:je,ks:ke)%graupel, io, irec )
  call dump( q_sbm(is:ie,js:je,ks:ke)%hail   , io, irec )
  call dump( re_sbm(is:ie,js:je,ks:ke)%liq    , io, irec )
  call dump( re_sbm(is:ie,js:je,ks:ke)%ice_col, io, irec )
  call dump( re_sbm(is:ie,js:je,ks:ke)%ice_pla, io, irec )
  call dump( re_sbm(is:ie,js:je,ks:ke)%ice_den, io, irec )
  call dump( re_sbm(is:ie,js:je,ks:ke)%snow   , io, irec )
  call dump( re_sbm(is:ie,js:je,ks:ke)%graupel, io, irec )
  call dump( re_sbm(is:ie,js:je,ks:ke)%hail   , io, irec )
 case default 
   ! do nothing
 end select mic_select1



 if( account_aerosol ) then
!   call dump( atmos(is:ie,js:je,ks:ke)%ccn    , io, irec )
!   call dump( atmos(is:ie,js:je,ks:ke)%icn    , io, irec )

 select case(trim(aerosol_microphysics) )
 case('GOCART')
  call dump( q_gocart(is:ie,js:je,ks:ke)%so4 * 1.e+6 , io, irec )
  call dump( q_gocart(is:ie,js:je,ks:ke)%blc * 1.e+6 , io, irec )
  call dump( (q_gocart(is:ie,js:je,ks:ke)%ocn &
             +  q_gocart(is:ie,js:je,ks:ke)%och ) * 1.e+6 , io, irec )
  call dump( (q_gocart(is:ie,js:je,ks:ke)%ssa &
             +  q_gocart(is:ie,js:je,ks:ke)%ssc ) * 1.e+6 , io, irec )
  call dump( (q_gocart(is:ie,js:je,ks:ke)%du1 &
           + q_gocart(is:ie,js:je,ks:ke)%du2 &
           + q_gocart(is:ie,js:je,ks:ke)%du3 &
           + q_gocart(is:ie,js:je,ks:ke)%du4 &
           + q_gocart(is:ie,js:je,ks:ke)%du5 &
           + q_gocart(is:ie,js:je,ks:ke)%du6 &
           + q_gocart(is:ie,js:je,ks:ke)%du7 &
           + q_gocart(is:ie,js:je,ks:ke)%du8 ) * 1.e+6  , io, irec )

 case('UCD')
  call dump( q_ucd(is:ie,js:je,ks:ke)%du1 * 1.e+6 , io, irec )
  call dump( q_ucd(is:ie,js:je,ks:ke)%du2 * 1.e+6 , io, irec )
  call dump( q_ucd(is:ie,js:je,ks:ke)%du3 * 1.e+6 , io, irec )
  call dump( q_ucd(is:ie,js:je,ks:ke)%du4 * 1.e+6 , io, irec )
  call dump( q_ucd(is:ie,js:je,ks:ke)%du5 * 1.e+6 , io, irec )

 end select
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
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//trim(typ)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//trim(typ)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  CRM 2D file for diagnostic'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef   ',mxlyr  ,' linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',irec-1

    write(io,*) 'tair   ',mxlyr,'  0  air temp  [K]'
    write(io,*) 'w      ',mxlyr,'  0  vertical velocity  [m/s]'
    write(io,*) 'rh     ',mxlyr,'  0  relative humidity  [%]'
    write(io,*) 'dhgt   ',mxlyr,'  0  layer thickness  [m]'
    write(io,*) 'hgt    ',mxlyr,'  0  layer height  [m]'
    write(io,*) 'precip ',mxlyr,'  0  precipitation [mm/hr]'

 mic_select3: select case(trim(type_microphysics))
 case('GEN')
    write(io,*) 'qc    ',mxlyr,'  0  cloud water  [g/m3]'
    write(io,*) 'qr    ',mxlyr,'  0  rain water  [g/m3]'
    write(io,*) 'qi    ',mxlyr,'  0  ice [g/m3]'
    write(io,*) 'qs    ',mxlyr,'  0  snow  [g/m3]'
    write(io,*) 'qg    ',mxlyr,'  0  graupel  [g/m3]'
    write(io,*) 'qh    ',mxlyr,'  0  hail  [g/m3]'
    write(io,*) 'rec   ',mxlyr,'  0  cloud water re [um]'
    write(io,*) 'rer   ',mxlyr,'  0  rain water re [um]'
    write(io,*) 'rei   ',mxlyr,'  0  ice re [um]'
    write(io,*) 'res   ',mxlyr,'  0  snow re [um]'
    write(io,*) 'reg   ',mxlyr,'  0  graupel re [um]'
    write(io,*) 'reh   ',mxlyr,'  0  hail re [um]'
 case('RAMS')
    write(io,*) 'qc1   ',mxlyr,'  0  cloud mode1 [g/m3]'
    write(io,*) 'qc2   ',mxlyr,'  0  cloud mode2 [g/m3]'
    write(io,*) 'qr    ',mxlyr,'  0  rain water  [g/m3]'
    write(io,*) 'qi1   ',mxlyr,'  0  ice mode1 [g/m3]'
    write(io,*) 'qi2   ',mxlyr,'  0  ice mode2 [g/m3]'
    write(io,*) 'qs    ',mxlyr,'  0  snow  [g/m3]'
    write(io,*) 'qg    ',mxlyr,'  0  graupel  [g/m3]'
    write(io,*) 'qh    ',mxlyr,'  0  hail  [g/m3]'
    write(io,*) 'rec1   ',mxlyr,'  0  cloud mode1 re [um]'
    write(io,*) 'rec2   ',mxlyr,'  0  cloud mode2 re [um]'
    write(io,*) 'rer    ',mxlyr,'  0  rain water  re [um]'
    write(io,*) 'rei1   ',mxlyr,'  0  ice mode1 re [um]'
    write(io,*) 'rei2   ',mxlyr,'  0  ice mode2 re [um]'
    write(io,*) 'res    ',mxlyr,'  0  snow re [um]'
    write(io,*) 'reg    ',mxlyr,'  0  graupel re [um]'
    write(io,*) 'reh    ',mxlyr,'  0  hail re [um]'
 case('SBM')
    write(io,*) 'ql     ',mxlyr,'  0  liquid water [g/m3]'
    write(io,*) 'qic    ',mxlyr,'  0  ice column [g/m3]'
    write(io,*) 'qip    ',mxlyr,'  0  ice plate [g/m3]'
    write(io,*) 'qid    ',mxlyr,'  0  ice dendrate [g/m3]'
    write(io,*) 'qs     ',mxlyr,'  0  snow [g/m3]'
    write(io,*) 'qg     ',mxlyr,'  0  graupel [g/m3]'
    write(io,*) 'qh     ',mxlyr,'  0  hail [g/m3]'
    write(io,*) 'rel    ',mxlyr,'  0  liquid water re [um]'
    write(io,*) 'reic   ',mxlyr,'  0  ice column re [um]'
    write(io,*) 'reip   ',mxlyr,'  0  ice plate re [um]'
    write(io,*) 'reid   ',mxlyr,'  0  ice dendrate re [um]'
    write(io,*) 'res    ',mxlyr,'  0  snow re [um]'
    write(io,*) 'reg    ',mxlyr,'  0  graupel re [um]'
    write(io,*) 'reh    ',mxlyr,'  0  hail re [um]'
 end select mic_select3


 if( account_aerosol ) then
 select case(trim(aerosol_microphysics) )
 case('GOCART')
     write(io,*) 'so4  ',mxlyr,' 0  so4 mixing ratio [ug/m3]'
     write(io,*) 'blc  ',mxlyr,' 0  black carbon mixing ratio [ug/m3]'
     write(io,*) 'oc   ',mxlyr,' 0  organic carbon mixing ratio [ug/m3]'
     write(io,*) 'ss   ',mxlyr,' 0  sea salt mixing ratio [ug/m3]'
     write(io,*) 'du   ',mxlyr,' 0  dust mixing ratio [ug/m3]'
 case('UCD')
    write(io,*) 'du1   ',mxlyr,' 0  mineral dust1 [ug/m3]'
    write(io,*) 'du2   ',mxlyr,' 0  mineral dust2 [ug/m3]'
    write(io,*) 'du3   ',mxlyr,' 0  mineral dust3 [ug/m3]'
    write(io,*) 'du4   ',mxlyr,' 0  mineral dust4 [ug/m3]'
    write(io,*) 'du5   ',mxlyr,' 0  mineral dust5 [ug/m3]'
 end select 
 endif

    write(io,*) 'endvars '
    close(io)
 endif

#if MPI == 2
 endif
#endif

 return
 end subroutine write_out_CRM3D

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_L2_radar
 implicit none

 integer :: ierr
 integer :: i,j,k,m,irec
 character(len=30) :: typ
 integer,parameter :: io = 110
 integer :: is,ie,js,je,ks,ke
 logical,parameter :: need_moments = .false.
!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr_radar

!
! write out CRM surface input for grads binarry format
!
 typ = '.'//trim(radar_sensor)//'.L2_RADAR'

!
! write out CRM input for grads binarry format
!
  sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//trim(typ)//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
  print("(a)"),'MSG write_out_L2_radar: write L2_radar parameters in grads format ',&
  ' ->',trim(sdsu_io_file)

!
! open CRM3D output file
!
#if MPI == 2
 if(masterproc) &
#endif
  open(unit= io, file = sdsu_io_file, access='direct', status='replace', recl=mxgridx*mxgridy*mxlyr_radar*4 )
  irec = 1

!
! write atmospheric parameter
!

  call dump( L2_radar_conv(is:ie,js:je,:)%press    , io, irec )
  call dump( L2_radar_conv(is:ie,js:je,:)%t_air    , io, irec )
  call dump( L2_radar_conv(is:ie,js:je,:)%sh       , io, irec )
  call dump( L2_radar_conv(is:ie,js:je,:)%w        , io, irec )
  call dump( L2_radar_conv(is:ie,js:je,:)%iwc      , io, irec )     
  call dump( L2_radar_conv(is:ie,js:je,:)%lwc      , io, irec )
  call dump( L2_radar_conv(is:ie,js:je,:)%rwc      , io, irec )
  call dump( L2_radar_conv(is:ie,js:je,:)%lfrac    , io, irec )
  call dump( L2_radar_conv(is:ie,js:je,:)%bulk_den , io, irec )
  call dump( L2_radar_conv(is:ie,js:je,:)%rain_rate, io, irec ) 
  if( radar_hid) then
     call dump( L2_radar_conv(is:ie,js:je,:)%hid      , io, irec )
     call dump( L2_radar_conv(is:ie,js:je,:)%hid_conf , io, irec )
  endif
  if(need_moments) then
     do m = 0,6 
     call dump( L2_radar_conv(is:ie,js:je,:)%mmt(m), io, irec )
     enddo
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
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//trim(typ)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//trim(typ)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  Radar Level 2 parameter in CRM grid'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef   ',mxlyr_radar  ,' linear 0. ',range_radar
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
  if(need_moments) then    
    if( radar_hid ) then
        write(io,*) 'vars 19'
    else 
        write(io,*) 'vars 17'
    endif
  else
    if( radar_hid ) then 
       write(io,*) 'vars 12'
    else
       write(io,*) 'vars 10'
    endif
  endif

    write(io,*) 'press   ',mxlyr_radar,'  0  layer pressure [hPa] '
    write(io,*) 't_air   ',mxlyr_radar,'  0  layer air temperature [degK]'
    write(io,*) 'sh      ',mxlyr_radar,'  0  specific humidity [g/g] (derived in get_others)'
    write(io,*) 'w       ',mxlyr_radar,'  0  vertical velocity (m/s)'
    write(io,*) 'iwc     ',mxlyr_radar,'  0  cloud ice water content [g/m3]'
    write(io,*) 'lwc     ',mxlyr_radar,'  0  cloud liquid water content [g/m3]'
    write(io,*) 'rwc     ',mxlyr_radar,'  0  rain water content (raining) [g/m3]'
    write(io,*) 'lfrac   ',mxlyr_radar,'  0  liquid fraction of rwc [-]'
    write(io,*) 'brho    ',mxlyr_radar,'  0  bulk density (raining) [g/cm3] '
    write(io,*) 'rain    ',mxlyr_radar,'  0  rain rate  [mm/hr] '
  if( radar_hid )  then
     write(io,*) 'hid    ',mxlyr_radar,'  0  HID (0~10) [-] '
     write(io,*) '* 0-UN, 1-DZ, 2-RN, 3-CR, 4-DS, 5-WS '
     write(io,*) '*    6-VI, 7-LDG, 8-HDG, 9-HA, 10-BD '
     write(io,*) 'hidcon ',mxlyr_radar,'  0  HID confidence [%] '
  endif

  if(need_moments) then
    write(io,*) 'mmt0    ',mxlyr_radar,'  0  0th moment [m0/m3] '
    write(io,*) 'mmt1    ',mxlyr_radar,'  0  1st moment [m1/m3] '
    write(io,*) 'mmt2    ',mxlyr_radar,'  0  2nd moment [m2/m3] '
    write(io,*) 'mmt3    ',mxlyr_radar,'  0  3rd moment [m3/m3] '
    write(io,*) 'mmt4    ',mxlyr_radar,'  0  4th moment [m4/m3] '
    write(io,*) 'mmt5    ',mxlyr_radar,'  0  5th moment [m5/m3] '
    write(io,*) 'mmt6    ',mxlyr_radar,'  0  6th moment [m6/m3] '
  endif
    write(io,*) 'endvars '
    close(io)
 endif

#if MPI == 2
 endif
#endif

 return
 end subroutine write_out_L2_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump2d( var2d, io, irec)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out 2D real array in GrADS (binary) format. 
! If size of var2d is smaller than domain size, var2d will be gathered into the domain2d.
!
! History:
! 11/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 real(sdsu_fps),dimension(myi_start:myi_end, myj_start:myj_end), intent(in)  :: var2d  ! 
 integer,intent(in)    :: io   ! io index
 integer,intent(inout) :: irec ! binary record #
 integer :: i,j



 if( myi_start==1 .and. myi_end==mxgridx .and. &
     myj_start==1 .and. myj_end==mxgridy )  then

     write(io,rec=irec) ((var2d(i,j),i=1,mxgridx),j=1,mxgridy)
     irec=irec+1

 else !domain decomposition case

#if MPI == 2
     call mpi_sdsu_collect_tile( var2d(myi_start:myi_end, myj_start:myj_end), &
                                 out_domain2d(1:mxgridx,1:mxgridy) )
     if(masterproc) write(io,rec=irec) ((out_domain2d(i,j),i=1,mxgridx),j=1,mxgridy)
     irec=irec+1
#endif

 endif

 return
 end subroutine dump2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
    
 subroutine dump3d( var3d, io, irec)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out 3D real array in GrADS (binary) format. 
! If size of var3d is smaller than domain size, var3d will be gathered into the domain3d via mpi library. 
! Vertical level can be either CRM (1~mxlyr) or satellite (1~variable) levels. 
!
! History:
! 11/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),dimension(myi_start:,myj_start:,1:), intent(in)  :: var3d  ! assumed 3D array input 
                                                                           !  (specify starting index)
 integer,intent(in)    :: io   ! io index
 integer,intent(inout) :: irec ! binary record #
 integer :: i,j,k
 integer :: bnd(3)      ! upper memory bound for 3D array
 integer :: ks, ke      ! starting and ending k index
 real,allocatable :: out_satellite_domain3d(:,:,:) !temporally 3D array for output purpose (satellite output)

!
! 1st, and 2nd array should be always myi_start:myi_end, myj_start:myj_end.
! otherwise stop program
! 3rd array could be CRM vertycal level (mxlyr) or satellite level (> mxlyr).
!
 bnd = UBOUND(var3d)  !check bounds 
 if( bnd(1) /= myi_end .or. bnd(2) /= myj_end ) then
    print*,'MSG dump3d: 1st, and 2nd array of var3d should be always myi_start:myi_end, myj_start:myj_end.'
    print*, 'bnd(1)=',bnd(1) , 'myi_end', myi_end , 'bnd(2)=', bnd(2) , 'myj_end=',myj_end
    call stop_sdsu(myrank,'MSG dump3d: terminate program') 
 endif
 ks = 1 ; ke = bnd(3) ! starting and ending k index


#if MPI < 2  

!    
! In this case var3d can be CRM or Satellite vertical levels 
!
    write(io,rec=irec) (((var3d(i,j,k),i=1,mxgridx),j=1,mxgridy),k=ks,ke)


#else

 if( ke /= mxlyr ) then  ! satellite vertical level (k can 1 ~ some variables ) 

   !
   ! allocate satellite 3D array for output purpose 
   !
   if( allocated(out_satellite_domain3d) ) then !if allocated check dimension bounds
       bnd = UBOUND(out_satellite_domain3d)  !check bounds
       if(bnd(3) /= ke) then  !if vertical level is different from previous step --> re-allocate
          deallocate( out_satellite_domain3d )
          allocate( out_satellite_domain3d(1:mxgridx,1:mxgridy,ks:ke) )
       endif
   else !if not allocated, allocate it.  
       allocate( out_satellite_domain3d(1:mxgridx,1:mxgridy,ks:ke) )
   endif

   !
   ! gather slave memory into master memory
   !
   call mpi_sdsu_collect_tile( ke, var3d(myi_start:myi_end, myj_start:myj_end, ks:ke),&
                               out_satellite_domain3d(1:mxgridx,1:mxgridy,ks:ke)  )
   !
   ! write satellite-dimension output
   !
   if(masterproc) write(io,rec=irec) (((out_satellite_domain3d(i,j,k),i=1,mxgridx),j=1,mxgridy),k=ks,ke)

 else  ! CRM vertical level  (k is always 1~mxlyr)

   !
   ! gather slave memory into master memory
   !
   call mpi_sdsu_collect_tile( ke, var3d(myi_start:myi_end, myj_start:myj_end, ks:ke),&
                                out_domain3d (1:mxgridx,1:mxgridy,1:mxlyr)  )
   !
   ! write CRM-dimension output
   !
   if(masterproc) write(io,rec=irec) (((out_domain3d(i,j,k),i=1,mxgridx),j=1,mxgridy),k=1,mxlyr)

 endif

#endif

 !
 ! add record # for next output
 !
 irec=irec+1

 return
 end subroutine dump3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump2d_latlon( var2d, lat_in, lon_in, io, irec)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out 2D real array in GrADS lat lon (binary) format. 
! 1st and 2nd dimension of var2d is treated as i and j horizontal loop.
! Input var2d, lat_in, lon_in must be in the same grid system. 
!
! History:
! 01/2011  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),dimension(:,:), intent(in)  :: var2d   ! assumed shape 2D array
 real(sdsu_fps),dimension(:,:), intent(in)  :: lat_in  ! assumed shape 2D array [deg]
 real(sdsu_fps),dimension(:,:), intent(in)  :: lon_in  ! assumed shape 2D array [deg]
 integer,intent(in)    :: io   ! io index
 integer,intent(inout) :: irec ! binary record #
 integer :: bnd(2)             ! upper memory bound for 3D array
 integer :: i,j                ! loop indice
 integer :: i_latlon,j_latlon

!
! initialize latlon_grid
!
 latlon_grid%var = undefined ; latlon_grid%vsum = 0. ;  latlon_grid%cnt = 0

!
! find upper bounds from assumped 2D array
!
 bnd = UBOUND(var2d)

!
! Fill value on lat lon grid
!
 do j = 1, bnd(2) ; do i = 1,bnd(1)

    if(var2d(i,j) == undefined) cycle

    wthin_domain: if( lat_in(i,j) >= minlat .and.  lat_in(i,j) <= maxlat ) then
      if( lon_in(i,j) >= minlon .and.  lon_in(i,j) <= maxlon ) then

       !
       ! loop for latlon_grid domain. 
       !
       do j_latlon = 1, jmax_latlon

         ! note it suppose to have *0.5 for res_latlon, but for visual reazon, I removed. it. (Toshi Matsui)
         within_grid: if( lat_in(i,j) > (lat_grid(j_latlon)-res_latlon) .and. &
                          lat_in(i,j) <= (lat_grid(j_latlon)+res_latlon) ) then
         do i_latlon = 1, imax_latlon

           if( lon_in(i,j) > (lon_grid(i_latlon)-res_latlon) .and. &
               lon_in(i,j) <= (lon_grid(i_latlon)+res_latlon) ) then

               latlon_grid(i_latlon,j_latlon)%vsum = latlon_grid(i_latlon,j_latlon)%vsum + var2d(i,j)
               latlon_grid(i_latlon,j_latlon)%cnt  = latlon_grid(i_latlon,j_latlon)%cnt + 1

           endif

         enddo
         endif within_grid
       enddo

      endif
    endif wthin_domain
 enddo ; enddo
   
!
! get mean value on lat-lon cordinate grid. 
!
 do j_latlon = 1, jmax_latlon
    do i_latlon = 1, imax_latlon
       if( latlon_grid(i_latlon,j_latlon)%cnt > 0 ) then
           latlon_grid(i_latlon,j_latlon)%var = latlon_grid(i_latlon,j_latlon)%vsum / &
                                                real(latlon_grid(i_latlon,j_latlon)%cnt)
       endif

       if( ISNAN( latlon_grid(i_latlon,j_latlon)%var ) ) then !clean up
           latlon_grid(i_latlon,j_latlon)%var = undefined
       endif

    enddo
 enddo

!
! write binary data
!
 write(io,rec=irec) ((latlon_grid(i_latlon,j_latlon)%var,i_latlon=1,imax_latlon),j_latlon=1,jmax_latlon)
 irec=irec+1  !update record length


 return
 end subroutine dump2d_latlon

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  fill_latlon_grid( lat, lon, var )
 implicit none
 real,intent(in) :: lat,lon ![deg]
 real,intent(in) :: var  !variable
 integer :: i_latlon,j_latlon

!
! loop for latlon_grid domain. 
!
 do j_latlon = 1, jmax_latlon

! note it suppose to have *0.5 for res_latlon, but for visual reazon, I removed. it. (Toshi Matsui)
    if( lat > (lat_grid(j_latlon)-res_latlon) .and. lat <= (lat_grid(j_latlon)+res_latlon) ) then
    do i_latlon = 1, imax_latlon

       if( lon > (lon_grid(i_latlon)-res_latlon) .and. lon <= (lon_grid(i_latlon)+res_latlon) ) then

           latlon_grid(i_latlon,j_latlon)%vsum = latlon_grid(i_latlon,j_latlon)%vsum + var
           latlon_grid(i_latlon,j_latlon)%cnt  = latlon_grid(i_latlon,j_latlon)%cnt + 1

       endif

    enddo
 endif
 enddo

 return
 end subroutine  fill_latlon_grid

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine dump2d_netcdf( var2d, varid, ncid)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out 2D real array in NetCDF format. 
! If size of var2d is smaller than domain size, var2d will be gathered into the domain2d via mpi library. 
!
! History:
! 11/2013  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),dimension(myi_start:myi_end, myj_start:myj_end), intent(in)  :: var2d  ! 
 integer,intent(in)    :: varid   ! var index
 integer,intent(inout) :: ncid    ! netcdf file id 
 integer :: i,j

! if( myi_start==1 .and. myi_end==mxgridx .and. &
!     myj_start==1 .and. myj_end==mxgridy )  then
#if MPI < 2  

     call check( nf90_put_var(ncid, varid, var2d) )

! else !domain decomposition case

#else 
   call mpi_sdsu_collect_tile(  var2d(myi_start:myi_end, myj_start:myj_end),&
                                out_domain2d (1:mxgridx,1:mxgridy)  )

     if(masterproc) call check( nf90_put_var(ncid, varid, out_domain2d) )
#endif

! endif


 return
 end subroutine dump2d_netcdf

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine dump3d_netcdf( var3d, varid, ncid)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out 3D real array in NetCDF format. 
! If size of var3d is smaller than domain size, var3d will be gathered into the domain3d via mpi library. 
! Vertical level can be either CRM (1~mxlyr) or satellite (1~variable) levels. 
!
! History:
! 11/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),dimension(myi_start:, myj_start:, 1:), intent(in)  :: var3d  ! assumed 3D array input
 integer,intent(in)    :: varid   ! var index
 integer,intent(inout) :: ncid    ! netcdf file id 

 integer :: i,j,k
 integer :: bnd(3)      ! upper memory bound for 3D array
 integer :: ks, ke      ! starting and ending k index
 real,allocatable :: out_satellite_domain3d(:,:,:) !temporally 3D array for output purpose (satellite output)


!
! 1st, and 2nd array should be always myi_start:myi_end, myj_start:myj_end.
! otherwise stop program
! 3rd array could be CRM vertycal level (mxlyr) or satellite level (> mxlyr).
!
 bnd = UBOUND(var3d)  !check bounds 
 if( bnd(1) /= myi_end .or. bnd(2) /= myj_end ) then
    print*,'MSG dump3d: 1st, and 2nd array of var3d should be always myi_start:myi_end, myj_start:myj_end.'
    print*, 'bnd(1)=',bnd(1) , 'myi_end', myi_end , 'bnd(2)=', bnd(2) , 'myj_end=',myj_end
    call stop_sdsu(myrank,'MSG dump3d: terminate program')
 endif
 ks = 1 ; ke = bnd(3) ! starting and ending k index


#if MPI < 2  

!    
! In this case var3d can be CRM or Satellite vertical levels 
!
     call check( nf90_put_var(ncid, varid, var3d) )

#else

 if( ke /= mxlyr ) then  ! satellite vertical level (k can 1 ~ some variables ) 

   !
   ! allocate satellite 3D array for output purpose 
   !
   if( allocated(out_satellite_domain3d) ) then !if allocated check dimension bounds
       bnd = UBOUND(out_satellite_domain3d)  !check bounds
       if(bnd(3) /= ke) then  !if vertical level is different from previous step --> re-allocate
          deallocate( out_satellite_domain3d )
          allocate( out_satellite_domain3d(1:mxgridx,1:mxgridy,ks:ke) )
       endif
   else !if not allocated, allocate it.  
       allocate( out_satellite_domain3d(1:mxgridx,1:mxgridy,ks:ke) )
   endif

   !
   ! gather slave memory into master memory
   !
   call mpi_sdsu_collect_tile( ke, var3d(myi_start:myi_end, myj_start:myj_end, ks:ke),&
                               out_satellite_domain3d(1:mxgridx,1:mxgridy,ks:ke)  )

   !
   ! write intrument-dimension output
   !
   if(masterproc) call check( nf90_put_var(ncid, varid, out_satellite_domain3d) )


 else  ! CRM vertical level  (k is always 1~mxlyr)

   !
   ! gather slave memory into master memory
   !
   call mpi_sdsu_collect_tile( ke, var3d(myi_start:myi_end, myj_start:myj_end, ks:ke),&
                                out_domain3d (1:mxgridx,1:mxgridy,1:mxlyr)  )
   !
   ! write CRM-dimension output
   !
   if(masterproc) call check( nf90_put_var(ncid, varid, out_domain3d ) )



 endif

#endif



 return
 end subroutine dump3d_netcdf

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine dump3d_cf( min_since, var3d, io )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:
! Write out 3D real array in GrADS (binary) format.
! If size of var3d is smaller than domain size, var3d will be gathered into the domain3d via mpi library.
! Vertical level can be either CRM (1~mxlyr) or satellite (1~variable) levels.
!
! History:
! 11/2009 Toshi Matsui@NASA GSFC ; Initial
!
! References:
!----------------------------------------------------------------------------------------------------
 integer,intent(in) :: min_since  ![min]

 real(sdsu_fps),dimension(myi_start:,myj_start:,1:), intent(in) :: var3d ! assumed 3D array input
                                                                           ! (specify starting index)
 integer,intent(in) :: io ! io index
 integer :: i,j,k
 integer :: bnd(3) ! upper memory bound for 3D array
 integer :: ks, ke ! starting and ending k index
 real,allocatable :: out_satellite_domain3d(:,:,:) !temporally 3D array for output purpose (satellite output)
 integer(8) :: id

!
! 1st, and 2nd array should be always myi_start:myi_end, myj_start:myj_end.
! otherwise stop program
! 3rd array could be CRM vertycal level (mxlyr) or satellite level (> mxlyr).
!
 bnd = UBOUND(var3d) !check bounds

 if( bnd(1) /= myi_end .or. bnd(2) /= myj_end ) then
    print*,'MSG dump3d_cf: 1st, and 2nd array of var3d should be always myi_start:myi_end, myj_start:myj_end.'
    print*, 'bnd(1)=',bnd(1) , 'myi_end', myi_end , 'bnd(2)=', bnd(2) , 'myj_end=',myj_end
    call stop_sdsu(myrank,'MSG dump3d_cf: terminate program')
 endif

 ks = 1 ; ke = bnd(3) ! starting and ending k index


#if MPI < 2  

!    
! In this case var3d can be CRM or Satellite vertical levels 
!
   id = 0
      do k=ks,ke ; do j = 1,mxgridy ; do i = 1,mxgridx
         id = id + 1
         write(io,400) id, min_since, k, lat_domain(i,j),lon_domain(i,j),var3d(i,j,k)
      enddo ; enddo ; enddo


#else



 if( ke /= mxlyr ) then ! satellite vertical level (k can 1 ~ some variables )

   !
   ! allocate satellite 3D array for output purpose
   !
   if( allocated(out_satellite_domain3d) ) then !if allocated check dimension bounds
       bnd = UBOUND(out_satellite_domain3d) !check bounds
       if(bnd(3) /= ke) then !if vertical level is different from previous step --> re-allocate
          deallocate( out_satellite_domain3d )
          allocate( out_satellite_domain3d(1:mxgridx,1:mxgridy,ks:ke) )
       endif
   else !if not allocated, allocate it.
       allocate( out_satellite_domain3d(1:mxgridx,1:mxgridy,ks:ke) )
   endif

   !
   ! gather slave memory into master memory
   !
   call mpi_sdsu_collect_tile( ke, var3d(myi_start:myi_end, myj_start:myj_end, ks:ke),&
                               out_satellite_domain3d(1:mxgridx,1:mxgridy,ks:ke) )
   !
   ! write satellite-dimension output in text file format
   !
   if(masterproc) then
      do k = ks, ke ; do j = 1,mxgridy ; do i = 1,mxgridx
         write(io,*) id, min_since, k, lat_domain(i,j),lon_domain(i,j), out_satellite_domain3d(i,j,k)
      enddo ; enddo ; enddo
   endif

 else ! CRM vertical level (k is always 1~mxlyr)


   !
   ! gather slave memory into master memory
   !
   call mpi_sdsu_collect_tile( ke, var3d(myi_start:myi_end, myj_start:myj_end, ks:ke),&
                                out_domain3d (1:mxgridx,1:mxgridy,1:mxlyr) )
   !
   ! write CRM-dimension output
   !

   id = 0
   if(masterproc) then
      do k = 1,mxlyr ; do j = 1,mxgridy ; do i = 1,mxgridx
         id = id + 1
         write(io,400) id, min_since, k, lat_domain(i,j),lon_domain(i,j),out_domain3d(i,j,k)
      enddo ; enddo ; enddo
   endif


 endif

#endif

 400 format(I22,1x,I12,1x,I4,1x,F10.5,1x,F10.5,1x,E14.8)

 return
 end subroutine dump3d_cf


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump2d_cf
 implicit none


 end subroutine dump2d_cf

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_sbm_moment
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Write out various momemnt from SBM inputs for studying particle-size distributions.  
! 
! History:
!  07/2011  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: ierr,irec
 integer :: i,j,k,m
 character(len=20) :: typ
 integer,parameter :: io = 111
 integer :: is,ie,js,je,ks,ke,di,dj,dk

!netcdf parameters
 integer :: ncid
 integer :: i_dimid, j_dimid, k_dimid 
 integer, parameter :: ndims3d = 3
 integer :: dimids(ndims3d)
 integer :: id
 type (particle_gen_i) :: &
     varid_mmt(1:3)  ,&   ! netcdf var ID for moments
     varid_rho       ,&   ! netcdf var ID for density papameters
     varid_fmelt          ! netcdf var ID for melt fraction

 character(len=3) :: char_mth
 character(len=7) :: char_unit

! moment parameters
 real(sdsu_fps),allocatable :: wgt(:) !weight
 real(sdsu_fps) :: mth !moment order
 real(sdsu_fps),allocatable  :: moment3d(:,:,:)
 real(sdsu_fps) :: &
    mmt     , &! moment [1/m3], [m3/m3], [m6/m3]
    mmt_col , &! moment [1/m3], [m3/m3], [m6/m3]
    mmt_pla , &! moment [1/m3], [m3/m3], [m6/m3]
    mmt_den    ! moment [1/m3], [m3/m3], [m6/m3]
 real(sdsu_fps) :: &
    rho     , &! bulk density [g/cm3]
    rho_col , &! bulk density [g/cm3]
    rho_pla , &! bulk density [g/cm3]
    rho_den    ! bulk density [g/cm3]
 real(sdsu_fps) :: fmelt

!
! check microphysics
!
 if( trim(cloud_microphysics) == 'HUCM_SBM43' .or. &
     trim(cloud_microphysics) == 'HUCM_SBM' ) then 
 else
     if(masterproc) print*, 'MSG write_out_sbm_moment: microphysics must be SBM'
     return
 endif
!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1

!
! allocate parameters
! 
 if( .not. allocated(wgt) ) then
   allocate( wgt(1:nbin) ) 
 endif

 if( .not. allocated(moment3d) ) then
    allocate( moment3d(is:ie,js:je,ks:ke) &
              , stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG write_out_sbm_moment: allocation error -> Terminate program.') 
 endif

!
! --------------- open netcdf output file and define metafile  ------------------------
!
#if MPI == 2
 if(masterproc) then 
#endif
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.PSD_MOMENT.nc'
 print("(a)"),'Moment output ->',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )

!
! Define 3D dimensions
!
 call check( nf90_def_dim(ncid, 'west_east'  , mxgridx, i_dimid) )
 call check( nf90_def_dim(ncid, 'south_north', mxgridy, j_dimid) )
 call check( nf90_def_dim(ncid, 'bottom_top' , mxlyr  , k_dimid) )
 dimids = (/ i_dimid, j_dimid, k_dimid /)

!
! moments parameteers
!
 MOMENT_LOOP1: do m = 1, 3
    mth = REAL(m-1)*3.e0  !moment order 
    if(m==1) then
      char_mth = '0th' ; char_unit = '[m0/m3]'
    elseif(m==2) then
      char_mth = '3rd' ; char_unit = '[m3/m3]'
    elseif(m==3) then
      char_mth = '6th' ; char_unit = '[m6/m3]'
    endif

    call check( nf90_def_var(ncid, 'm'//char_mth(1:1)//'_cloud', NF90_REAL, dimids, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','xyz') )
    call check( nf90_put_att(ncid, id, 'description',char_mth//'-moment of PSD (cloud)') )
    call check( nf90_put_att(ncid, id, 'units', char_unit) )
    varid_mmt(m)%cloud = id

    call check( nf90_def_var(ncid, 'm'//char_mth(1:1)//'_rain', NF90_REAL, dimids, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','xyz') )
    call check( nf90_put_att(ncid, id, 'description',char_mth//'-moment of PSD (rain)') )
    call check( nf90_put_att(ncid, id, 'units', char_unit) )
    varid_mmt(m)%rain = id

    call check( nf90_def_var(ncid, 'm'//char_mth(1:1)//'_ice', NF90_REAL, dimids, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','xyz') )
    call check( nf90_put_att(ncid, id, 'description',char_mth//'-moment of PSD (ice)') )
    call check( nf90_put_att(ncid, id, 'units', char_unit) )
    varid_mmt(m)%ice = id

    call check( nf90_def_var(ncid, 'm'//char_mth(1:1)//'_snow', NF90_REAL, dimids, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','xyz') )
    call check( nf90_put_att(ncid, id, 'description',char_mth//'-moment of PSD (snow)') )
    call check( nf90_put_att(ncid, id, 'units', char_unit) )
    varid_mmt(m)%snow = id

    call check( nf90_def_var(ncid, 'm'//char_mth(1:1)//'_graupel', NF90_REAL, dimids, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','xyz') )
    call check( nf90_put_att(ncid, id, 'description',char_mth//'-moment of PSD (graupel)') )
    call check( nf90_put_att(ncid, id, 'units', char_unit) )
    varid_mmt(m)%graupel = id

    call check( nf90_def_var(ncid, 'm'//char_mth(1:1)//'_hail', NF90_REAL, dimids, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','xyz') )
    call check( nf90_put_att(ncid, id, 'description',char_mth//'-moment of PSD (hail)') )
    call check( nf90_put_att(ncid, id, 'units', char_unit) )
    varid_mmt(m)%hail = id

 enddo MOMENT_LOOP1

!
! density parameters 
!
    call check( nf90_def_var(ncid, 'rho_ice', NF90_REAL, dimids, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','xyz') )
    call check( nf90_put_att(ncid, id, 'description','bulk density (ice)') )
    call check( nf90_put_att(ncid, id, 'units','[g/cm3]') )
    varid_rho%ice = id

    call check( nf90_def_var(ncid, 'rho_snow', NF90_REAL, dimids, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','xyz') )
    call check( nf90_put_att(ncid, id, 'description','bulk density (snow)') )
    call check( nf90_put_att(ncid, id, 'units', '[g/cm3]') )
    varid_rho%snow = id

!
! bulk melt fraction
!
    call check( nf90_def_var(ncid, 'fmelt_snow', NF90_REAL, dimids, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','xyz') )
    call check( nf90_put_att(ncid, id, 'description','bulk melt fraction (snow)') )
    call check( nf90_put_att(ncid, id, 'units', '[-]') )
    varid_fmelt%snow = id

  call check( nf90_enddef(ncid) )  ! End define mode.

#if MPI == 2
 endif
#endif

!
! --------------- write out parameters in  netcdf file  ------------------------
!

!
! 0th, 3rd, 6th-order moments
!
 MOMENT_LOOP: do m = 1, 3

    mth = REAL(m-1)*3.e0  ! moment order 
    wgt(1:nbin) = 1.0e0   ! weight

    do k = ks, ke ; do j = js, je ; do i = is, ie
       call sbm_moment('cloud', mth, wgt(1:nbin), rad_sbm(1:nbin)%liq, drad_sbm(1:nbin)%liq,&
                        q_sbm(i,j,k)%liq, n_sbm(i,j,k,1:nbin)%liq, mmt ) 
       moment3d(i,j,k)= mmt 
    enddo ; enddo ; enddo
    call dumpnc( moment3d(is:ie,js:je,ks:ke) , varid_mmt(m)%cloud, ncid )

    do k = ks, ke ; do j = js, je ; do i = is, ie
       call sbm_moment('rain', mth, wgt(1:nbin), rad_sbm(1:nbin)%liq, drad_sbm(1:nbin)%liq,&
                       q_sbm(i,j,k)%liq, n_sbm(i,j,k,1:nbin)%liq, mmt )
       moment3d(i,j,k)= mmt 
    enddo ; enddo ; enddo
    call dumpnc( moment3d(is:ie,js:je,ks:ke) , varid_mmt(m)%rain, ncid )

    do k = ks, ke ; do j = js, je ; do i = is, ie
       call sbm_moment('ice', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_col, drad_sbm(1:nbin)%ice_col,&
                        q_sbm(i,j,k)%ice_col, n_sbm(i,j,k,1:nbin)%ice_col, mmt_col )
       call sbm_moment('ice', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_pla, drad_sbm(1:nbin)%ice_pla,&
                        q_sbm(i,j,k)%ice_pla, n_sbm(i,j,k,1:nbin)%ice_pla, mmt_pla )
       call sbm_moment('ice', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_den, drad_sbm(1:nbin)%ice_den,&
                        q_sbm(i,j,k)%ice_den, n_sbm(i,j,k,1:nbin)%ice_den, mmt_den )
       moment3d(i,j,k)= mmt_col + mmt_pla + mmt_den
    enddo ; enddo ; enddo
    call dumpnc( moment3d(is:ie,js:je,ks:ke) , varid_mmt(m)%ice, ncid )

    do k = ks, ke ; do j = js, je ; do i = is, ie
       call sbm_moment('snow', mth, wgt(1:nbin), rad_sbm(1:nbin)%snow, drad_sbm(1:nbin)%snow,&
                       q_sbm(i,j,k)%snow, n_sbm(i,j,k,1:nbin)%snow, mmt )
       moment3d(i,j,k)= mmt 
    enddo ; enddo ; enddo
    call dumpnc( moment3d(is:ie,js:je,ks:ke) , varid_mmt(m)%snow, ncid )

    do k = ks, ke ; do j = js, je ; do i = is, ie
       call sbm_moment('graupel', mth, wgt(1:nbin), rad_sbm(1:nbin)%graupel, drad_sbm(1:nbin)%graupel,&
                        q_sbm(i,j,k)%graupel, n_sbm(i,j,k,1:nbin)%graupel, mmt )
       moment3d(i,j,k)= mmt 
    enddo ; enddo ; enddo
    call dumpnc( moment3d(is:ie,js:je,ks:ke) , varid_mmt(m)%graupel, ncid )

    do k = ks, ke ; do j = js, je ; do i = is, ie
       call sbm_moment('hail', mth, wgt(1:nbin), rad_sbm(1:nbin)%hail, drad_sbm(1:nbin)%hail,&
                        q_sbm(i,j,k)%hail, n_sbm(i,j,k,1:nbin)%hail, mmt )
       moment3d(i,j,k)= mmt 
    enddo ; enddo ; enddo
    call dumpnc( moment3d(is:ie,js:je,ks:ke) , varid_mmt(m)%hail, ncid )

 enddo MOMENT_LOOP


!
! bulk density for ice and snow
!
    do k = ks, ke ; do j = js, je ; do i = is, ie
       mth = 3.e0            ! moment order 
       wgt(1:nbin) = 1.0e0   ! weight
       call sbm_moment('ice', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_col, drad_sbm(1:nbin)%ice_col,&
                        q_sbm(i,j,k)%ice_col, n_sbm(i,j,k,1:nbin)%ice_col, mmt_col )
       call sbm_moment('ice', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_pla, drad_sbm(1:nbin)%ice_pla,&
                        q_sbm(i,j,k)%ice_pla, n_sbm(i,j,k,1:nbin)%ice_pla, mmt_pla )
       call sbm_moment('ice', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_den, drad_sbm(1:nbin)%ice_den,&
                        q_sbm(i,j,k)%ice_den, n_sbm(i,j,k,1:nbin)%ice_den, mmt_den )
       mmt = mmt_col + mmt_pla + mmt_den

       wgt(1:nbin) = den_sbm(1:nbin)%ice_col   ! weight (density)
       call sbm_moment('ice', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_col, drad_sbm(1:nbin)%ice_col,&
                        q_sbm(i,j,k)%ice_col, n_sbm(i,j,k,1:nbin)%ice_col, rho_col )
       wgt(1:nbin) = den_sbm(1:nbin)%ice_pla   ! weight (density)
       call sbm_moment('ice', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_pla, drad_sbm(1:nbin)%ice_pla,&
                        q_sbm(i,j,k)%ice_pla, n_sbm(i,j,k,1:nbin)%ice_pla, rho_pla )
       wgt(1:nbin) = den_sbm(1:nbin)%ice_den   ! weight (density)
       call sbm_moment('ice', mth, wgt(1:nbin), rad_sbm(1:nbin)%ice_den, drad_sbm(1:nbin)%ice_den,&
                        q_sbm(i,j,k)%ice_den, n_sbm(i,j,k,1:nbin)%ice_den, rho_den )
       rho = rho_col + rho_pla + rho_den

       if( mmt <= 0.) then
         moment3d(i,j,k) = 0.e0 
       else    
         moment3d(i,j,k) = rho / mmt !bulk ice density [g/cm3]
       endif

    enddo ; enddo ; enddo
    call dumpnc( moment3d(is:ie,js:je,ks:ke) , varid_rho%ice, ncid )

    do k = ks, ke ; do j = js, je ; do i = is, ie
       mth = 3.e0    ! moment order 
       wgt(1:nbin) = 1.0e0   ! weight
       call sbm_moment('snow', mth, wgt(1:nbin), rad_sbm(1:nbin)%snow, drad_sbm(1:nbin)%snow,&
                       q_sbm(i,j,k)%snow, n_sbm(i,j,k,1:nbin)%snow, mmt )

       wgt(1:nbin) = den_sbm(1:nbin)%snow ! weight (density)
       call sbm_moment('snow', mth, wgt(1:nbin), rad_sbm(1:nbin)%snow, drad_sbm(1:nbin)%snow,&
                       q_sbm(i,j,k)%snow, n_sbm(i,j,k,1:nbin)%snow, rho )

       if( mmt <= 0.) then
         moment3d(i,j,k) = 0.e0
       else
         moment3d(i,j,k) = rho / mmt !bulk snowdensity [g/cm3]
       endif

    enddo ; enddo ; enddo
    call dumpnc( moment3d(is:ie,js:je,ks:ke) , varid_rho%snow, ncid )

!
! bulk melt fraction
!
    do k = ks, ke ; do j = js, je ; do i = is, ie
       mth = 3.e0    ! moment order 
       wgt(1:nbin) = den_sbm(1:nbin)%snow * fmelt_sbm_snow(i,j,k,1:nbin)   ! weight
       call sbm_moment('snow', mth, wgt(1:nbin), rad_sbm(1:nbin)%snow, drad_sbm(1:nbin)%snow,&
                       q_sbm(i,j,k)%snow, n_sbm(i,j,k,1:nbin)%snow, fmelt )

       wgt(1:nbin) = den_sbm(1:nbin)%snow ! weight (density)
       call sbm_moment('snow', mth, wgt(1:nbin), rad_sbm(1:nbin)%snow, drad_sbm(1:nbin)%snow,&
                       q_sbm(i,j,k)%snow, n_sbm(i,j,k,1:nbin)%snow, rho )

       if( rho <= 0.) then
         moment3d(i,j,k) = 0.e0
       else
         moment3d(i,j,k) = fmelt / rho  !bulk melt fraction for snowy [-]
       endif

    enddo ; enddo ; enddo
    call dumpnc( moment3d(is:ie,js:je,ks:ke) , varid_fmelt%snow, ncid )


!
! close netcdf file
!
#if MPI == 2
 if(masterproc) &
#endif
  call check( nf90_close(ncid) )

 return
 end subroutine write_out_sbm_moment

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

! subroutine gce_moment( mth, q, re, rho, moment)
! implicit none
!----------------------------------------------------------------------------------------------------
! Purpose 
!   compute mth moment from SBM bin input. 
!
! 03/2012  Toshi Matsui @ NASA GSFC ; initial
!----------------------------------------------------------------------------------------------------
! real(sdsu_fps), intent(in) :: mth              ! moment order
! real(sdsu_fps), intent(in) :: q                ! mass mixing ratio [g/m3]
! real(sdsu_fps), intent(in) :: re               ! effective radius [micron]
! real(sdsu_fps), intent(in) :: rho              ! effectie density of hydrometero [kg/m3]
! real(sdsu_fps), intent(out) :: moment          ! ouput moments

! real(sdsu_fps) :: lam  ! slope [1/m]
! real(sdsu_fps) :: n0   ! intercept [1/m4]
! real(sdsu_fps) :: mu   ! mu of gamma PSD [-]
! real(sdsu_fps) :: gfac1  ! gamma factor [-] 
! real(sdsu_fps) :: order 

! if(q < q_min_condensate .or. re == 0. ) then
!    moment = 0.
!    return
! endif

! lam = 3./2./(re*1e-6)  ![1/m]

 !n0 = (lam**4)  * q * 1e-3 / (const_pi*rho) ![1/m4]

! mu = 0.
! order = mu + mth + 1

! call gamma_reff(order,gfac1)

! moment = n0 * gfac1 / ( lam**order )  !mth-order moment


! return
! end subroutine gce_moment 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gen_moment( mth, q, re, rho, mu, moment)
 implicit none
!----------------------------------------------------------------------------------------------------
! Purpose 
!   compute mth moment from generalized gamma distribution. 
!
! 03/2013  Toshi Matsui @ NASA GSFC ; initial
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps), intent(in) :: mth              ! moment order
 real(sdsu_fps), intent(in) :: q                ! mass mixing ratio [g/m3]
 real(sdsu_fps), intent(in) :: re               ! effective radius [micron]
 real(sdsu_fps), intent(in) :: rho              ! effectie density of hydrometero [kg/m3]
 real(sdsu_fps), intent(in) :: mu               ! mu of gamma PSD [-]
 real(sdsu_fps), intent(out) :: moment          ! ouput moments

 real(sdsu_fps) :: lam  ! slope [1/m]
 real(sdsu_fps) :: n0   ! intercept [1/m4]
 real(sdsu_fps) :: gfac1,gfac2,gamfac  ! gamma factor [-] 
 real(sdsu_fps) :: order

 if(q < q_min_condensate .or. re == 0. ) then
    moment = 0.
    return
 endif

!
! derive lambda and intercept of expoential DSD for unit water content
!
  call gamma_reff(mu+4.0,gfac1)
  call gamma_reff(mu+3.0,gfac2)
  gamfac = gfac1/gfac2

!
! derive lambda
!
  lam = 1. / (2.*re*1.e-6) * gamfac  !slope [1/m]

!
! derive intercept
!       [g/m3]              [m3/g]        [1/m]
  n0  = q * 6. / const_pi / (rho*1000.) * (lam**(4.+mu))  / gfac1  ! [1/m**(4+mu)]

 order = mu + mth + 1

 call gamma_reff(order,gfac1)

!
! compute m-th moment 
!
 moment = n0 * gfac1 / ( lam**order )  !mth-order moment [m**(mth-3)]


 return
 end subroutine gen_moment

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine sbm_moment(spec, mth, wgt , bulkrad, dbulkrad, q, n0, &
                       moment)
 implicit none
!----------------------------------------------------------------------------------------------------
! Purpose 
!   compute mth moment from SBM bin input. 
!
! 08/2007  Toshi Matsui @ NASA GSFC ; initial
!----------------------------------------------------------------------------------------------------
 character(len=*),intent(in):: spec             ! index name of hydrometeor classes 
 real(sdsu_fps), intent(in) :: mth              ! moment order
 real(sdsu_fps), intent(in) :: wgt     (1:nbin) ! weight (density or )
 real(sdsu_fps), intent(in) :: bulkrad (1:nbin) ! bulk radius of hydrometero for each bins [cm] 
 real(sdsu_fps), intent(in) :: dbulkrad(1:nbin) ! binsize [cm]
 real(sdsu_fps), intent(in) :: q                ! total water content [g/m3]
 real(sdsu_fps), intent(in) :: n0(1:nbin)       ! particle number density [1/m4]
 real(sdsu_fps), intent(out) :: moment          ! moment output   [m**(mth-3)]

 integer :: ib, ibs,ibe  !looping index
 real(sdsu_fps) :: Di    !bulk diameter [m]
 real(sdsu_fps) :: Ni    !number concentration [#/m3]
 real(sdsu_fps) :: Wi    !weight (could be one, density, melt frac, etc..)
 real(sdsu_fps) :: sum_moment

!
! if no mixing ratio, zero out and return. 
!
 if( q <= 0.e0 ) then
     moment = 0.  
     return
 endif

!
! set up bin bounds
!
 do ib = 1, nbin-1
    if( bulkrad(ib)*1.e4 <= 100. .and. bulkrad(ib+1)*1.e4 > 100. ) then  !100um is threshold of rain radius
        if( trim(spec) == 'rain') then
            ibs = ib+1  ; ibe = nbin
        elseif( trim(spec) == 'non_rain') then
            ibs = 1     ; ibe = ib
        endif
        exit
    endif
 enddo


!
! integrating size bin
!
 sum_moment = 0.e0

 do ib = ibs, ibe

    Ni  = 0.01e0 * dbulkrad(ib) * n0(ib)  ! number conc [#/m3]
    Di  = 2.e0   *  bulkrad(ib) * 0.01e0  ! bulk diameter [m]
    Wi  = wgt(ib)                         ! unit could be different

    sum_moment = sum_moment + Wi * ( Di ** mth ) * Ni   !summing N-th moment

 enddo

!
! output Mth moment
!
 moment = sum_moment  !  [m**(mth-3)]

 return
 end subroutine sbm_moment

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine re_all
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Choose subrouutines for different microphysics schme in order to compute effective radius of 
! condensates.  
! 
! History:
! 05/2017  Toshi Matsui@NASA GSFC ; fixed bugs in WDM call in units. (suggested from Jue-Hye Kim)
! 08/2016  Toshi Matsui@NASA GSFC ; added new WSM/WDM routine.
! 08/2008  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: i,j,k !loop indice
 real(sdsu_fps) :: n0 !intercept [1/m**4]
 real(sdsu_fps) :: rhoe !effective density 
 real(sdsu_fps) :: Nt_temp

 real(sdsu_fps) :: mu_c, lam_c ! temp

!
!
!


 if(masterproc) print*,'MSG re_all: compute particle effective radius'
 if(masterproc) print*,' '


 RAMS: if( trim(cloud_microphysics) == 'RAMS1' .or. trim(cloud_microphysics) == 'RAMS2' ) then

    call re_rams_gamma('re')

 else

    do k = 1,mxlyr ; do j=myj_start, myj_end ; do i=myi_start, myi_end

       mic_select: select case(trim(cloud_microphysics))
       case('GOD') ! Exponential DSD prescribed parameters for Goddard Microphysics

          re_gen(i,j,k)%cloud = r_cld  !mono distribution    
          if(q_gen(i,j,k)%cloud == 0.) re_gen(i,j,k)%cloud = 0.

          call re_LUT_Heymsfield_Platt_1984('proc', atmos(i,j,k)%t_air, q_gen(i,j,k)%ice, re_gen(i,j,k)%ice)

          call re_gamma_1mmt(rho_gen%rain  , n0_gen%rain  , q_gen(i,j,k)%rain , re_gen(i,j,k)%rain )

          call re_gamma_1mmt(rho_gen%snow  , n0_gen%snow  , q_gen(i,j,k)%snow , re_gen(i,j,k)%snow )


          call re_gamma_1mmt(rho_gen%graupel , n0_gen%graupel , q_gen(i,j,k)%graupel , re_gen(i,j,k)%graupel )
          call re_gamma_1mmt(rho_gen%hail  , n0_gen%hail  , q_gen(i,j,k)%hail , re_gen(i,j,k)%hail )

       case('GOD10') ! Exponential DSD prescribed parameters for Goddard Microphysics 2010

          re_gen(i,j,k)%cloud = r_cld  !mono distribution    
          if(q_gen(i,j,k)%cloud == 0.) re_gen(i,j,k)%cloud = 0.

          call re_LUT_Heymsfield_Platt_1984('proc', atmos(i,j,k)%t_air, q_gen(i,j,k)%ice, re_gen(i,j,k)%ice)

          call re_gamma_1mmt(rho_gen%rain  , n0_gen%rain  , q_gen(i,j,k)%rain , re_gen(i,j,k)%rain )

!          call re_gamma_1mmt(rho_gen%snow  , n0_gen%snow  , q_gen(i,j,k)%snow , re_gen(i,j,k)%snow )

          call re_bulk_god10_snow   ( atmos(i,j,k)%t_air, rho_gen%snow, q_gen(i,j,k)%snow, re_gen(i,j,k)%snow )

          call re_bulk_god10_graupel( atmos(i,j,k)%t_air, rho_gen%graupel, q_gen(i,j,k)%graupel, re_gen(i,j,k)%graupel )

          call re_gamma_1mmt(rho_gen%hail  , n0_gen%hail  , q_gen(i,j,k)%hail , re_gen(i,j,k)%hail )

       case('GMP4ICE')  ! Goddard Microphysics 4ICE 2014

          re_gen(i,j,k)%cloud = r_cld  !mono distribution    
          if(q_gen(i,j,k)%cloud == 0.) re_gen(i,j,k)%cloud = 0.

          call re_LUT_Heymsfield_Platt_1984('proc', atmos(i,j,k)%t_air, q_gen(i,j,k)%ice, re_gen(i,j,k)%ice)

          call re_gamma_1mmt(rho_gen%rain  , n0_gen%rain  , q_gen(i,j,k)%rain , re_gen(i,j,k)%rain )
          !bin correction for rain drop
          call gmp13_bin_corretion( atmos(i,j,k)%t_air, atmos(i,j,k)%rh, q_gen(i,j,k)%cloud, &
                                    q_gen(i,j,k)%rain , atmos(i,j,k)%rho_dair, re_gen(i,j,k)%rain )

          !4ICE diagnostic size-density mapping for snow
          call gmp13_size_rho_diag ('snow', atmos(i,j,k)%t_air, q_gen(i,j,k)%snow, re_gen(i,j,k)%snow, rhoe)
          !4ICE diagnostic size mapping for graupel
          call gmp13_size_rho_diag ('graupel', atmos(i,j,k)%t_air, q_gen(i,j,k)%graupel, re_gen(i,j,k)%graupel, rhoe)

          call re_gamma_1mmt(rho_gen%hail  , n0_gen%hail  , q_gen(i,j,k)%hail , re_gen(i,j,k)%hail )

       case('TED') ! Exponential DSD prescribed parameters for Goddard Microphysics plus TEDD scheme

          re_gen(i,j,k)%cloud = r_cld  !mono distribution    
          if(q_gen(i,j,k)%cloud == 0.) re_gen(i,j,k)%cloud = 0.
          call re_LUT_Heymsfield_Platt_1984('proc', atmos(i,j,k)%t_air, q_gen(i,j,k)%ice,  re_gen(i,j,k)%ice )
          call re_gamma_1mmt(rho_gen%rain  , n0_gen%rain  , q_gen(i,j,k)%rain  , re_gen(i,j,k)%rain )

!                                                                  Tctl,  a                  , b
          call re_bulk_tedd(atmos(i,j,k)%t_air, q_gen(i,j,k)%snow, 273.15, 3000.*qcol_gen(i,j)%snow, 200.       , &
                            3000. , re_gen(i,j,k)%snow    )
          call re_bulk_tedd(atmos(i,j,k)%t_air, q_gen(i,j,k)%graupel, 273.15, 5000. , 15.*qcol_gen(i,j)%graupel+200., &
                            3000. , re_gen(i,j,k)%graupel ) 
          call re_bulk_tedd(atmos(i,j,k)%t_air, q_gen(i,j,k)%hail   , 273.15, 5000.   , 30.*qcol_gen(i,j)%hail,      &
                            5000. , re_gen(i,j,k)%hail    ) 


       case('LIN') ! Exponential DSD prescribed parameters for LIN scheme 
      
 
          re_gen(i,j,k)%cloud = r_cld  !mono distribution    
          if(q_gen(i,j,k)%cloud == 0.) re_gen(i,j,k)%cloud = 0.
          call re_LUT_Heymsfield_Platt_1984('proc', atmos(i,j,k)%t_air, q_gen(i,j,k)%ice, re_gen(i,j,k)%ice )
          call re_gamma_1mmt(rho_gen%rain  , n0_gen%rain  , q_gen(i,j,k)%rain , re_gen(i,j,k)%rain )
          call re_gamma_1mmt(rho_gen%snow  , n0_gen%snow  , q_gen(i,j,k)%snow , re_gen(i,j,k)%snow )
          call re_gamma_1mmt(rho_gen%graupel , n0_gen%graupel , q_gen(i,j,k)%graupel , re_gen(i,j,k)%graupel )
          call re_gamma_1mmt(rho_gen%hail  , n0_gen%hail  , q_gen(i,j,k)%hail , re_gen(i,j,k)%hail )

       case('WSM') ! Exponential DSD prescribed parameters for WSM scheme 

          ! subroutine for wdm cloud, ice, and snow aggregate.
          call effectRad_wdm6_0d ( .false., .false., .false., &
                                   atmos(i,j,k)%t_air, atmos(i,j,k)%rho_dair*1.e-3, &
                                   q_gen(i,j,k)%cloud / ( atmos(i,j,k)%rho_dair ) , &
                                   3.E8 , &  !constant for single WSM
                                   q_gen(i,j,k)%ice   / ( atmos(i,j,k)%rho_dair ) , &
                                   q_gen(i,j,k)%snow / ( atmos(i,j,k)%rho_dair ) , &
                                   re_gen(i,j,k)%cloud, re_gen(i,j,k)%ice, re_gen(i,j,k)%snow )

          call re_gamma_1mmt(rho_gen%rain  , n0_gen%rain  , q_gen(i,j,k)%rain  , re_gen(i,j,k)%rain )
          call re_gamma_1mmt(rho_gen%graupel , n0_gen%graupel , q_gen(i,j,k)%graupel , re_gen(i,j,k)%graupel )
          call re_gamma_1mmt(rho_gen%hail  , n0_gen%hail  , q_gen(i,j,k)%hail , re_gen(i,j,k)%hail )

       case('WDM') ! Exponential DSD prescribed parameters for WSM scheme 

          ! subroutine for wdm cloud, ice, and snow aggregate.
          call effectRad_wdm6_0d ( wdm_re_cloud_exist, wdm_re_ice_exist, wdm_re_snow_exist, &
                                   atmos(i,j,k)%t_air, atmos(i,j,k)%rho_dair*1.e-3, &
                                   q_gen(i,j,k)%cloud / ( atmos(i,j,k)%rho_dair ) , &
                                   n_gen(i,j,k)%cloud / ( atmos(i,j,k)%rho_dair ) , &
                                   q_gen(i,j,k)%ice   / ( atmos(i,j,k)%rho_dair ) , &
                                   q_gen(i,j,k)%snow / ( atmos(i,j,k)%rho_dair ) , &
                                   re_gen(i,j,k)%cloud, re_gen(i,j,k)%ice, re_gen(i,j,k)%snow )

          call re_gamma_2mmt( 'qr',mu_gen%rain, rho_gen%rain, q_gen(i,j,k)%rain, n_gen(i,j,k)%rain, re_gen(i,j,k)%rain )
          call re_gamma_1mmt(rho_gen%graupel , n0_gen%graupel , q_gen(i,j,k)%graupel , re_gen(i,j,k)%graupel )
          call re_gamma_1mmt(rho_gen%hail  , n0_gen%hail  , q_gen(i,j,k)%hail , re_gen(i,j,k)%hail )

        case('THOM')

          !cloud
          Nt_temp = Nt_c_thom
          call re_gamma_2mmt( 'qc',mu_gen%cloud , rho_gen%cloud , q_gen(i,j,k)%cloud , Nt_temp , re_gen(i,j,k)%cloud )

          !precip species
          call thom_re( q_gen(i,j,k)%rain, n_gen(i,j,k)%rain, q_gen(i,j,k)%snow, q_gen(i,j,k)%graupel, &
                        atmos(i,j,k)%t_air, re_gen(i,j,k)%rain, re_gen(i,j,k)%snow, re_gen(i,j,k)%graupel)

          !ice crystal
          call re_gamma_2mmt( 'qi',mu_gen%ice , rho_gen%ice , q_gen(i,j,k)%ice , n_gen(i,j,k)%ice , re_gen(i,j,k)%ice  )

          !snow, graupel

        case('MORR','MORRH')
      
          Nt_temp = Nt_c_thom !tempoeral fix
          call re_gamma_2mmt( 'qc',mu_gen%cloud , rho_gen%cloud , q_gen(i,j,k)%cloud , Nt_temp , re_gen(i,j,k)%cloud )
          call re_gamma_2mmt( 'qr',mu_gen%rain, rho_gen%rain, q_gen(i,j,k)%rain, n_gen(i,j,k)%rain, re_gen(i,j,k)%rain )
          call re_gamma_2mmt( 'qi',mu_gen%ice , rho_gen%ice , q_gen(i,j,k)%ice , n_gen(i,j,k)%ice , re_gen(i,j,k)%ice  )
          call re_gamma_2mmt( 'qs',mu_gen%snow, rho_gen%snow, q_gen(i,j,k)%snow, n_gen(i,j,k)%snow, re_gen(i,j,k)%snow )
          call re_gamma_2mmt( 'qg',mu_gen%graupel, rho_gen%graupel, q_gen(i,j,k)%graupel, n_gen(i,j,k)%graupel, &
                              re_gen(i,j,k)%graupel )
          call re_gamma_2mmt( 'qh',mu_gen%hail, rho_gen%hail, q_gen(i,j,k)%hail, n_gen(i,j,k)%hail, re_gen(i,j,k)%hail )

        case('P3')

          call re_cloud_p3( q_gen(i,j,k)%cloud, n_gen(i,j,k)%cloud ,atmos(i,j,k)%rho_dair, p3_mu_cloud(i,j,k), re_gen(i,j,k)%cloud ) 
          call re_rain_p3 ( q_gen(i,j,k)%rain , n_gen(i,j,k)%rain, atmos(i,j,k)%rho_dair, re_gen(i,j,k)%rain )
          call re_ice_p3  ( q_gen(i,j,k)%ice , n_gen(i,j,k)%ice, p3_frime(i,j,k), p3_rhorime(i,j,k), re_gen(i,j,k)%ice ) 

        case('NSSL2MMT')

          !cloud, ice, and snow species use N(volume) PSD.  
          call re_nssl_volume( q_gen(i,j,k)%cloud, q_gen(i,j,k)%ice, q_gen(i,j,k)%snow, &
                               n_gen(i,j,k)%cloud, n_gen(i,j,k)%ice, n_gen(i,j,k)%snow, &
                               mu_gen%cloud,       mu_gen%ice,       mu_gen%snow,       &
                               rho_gen%cloud,      rho_gen%ice,      rho_gen%snow,      &
                               re_gen(i,j,k)%cloud,re_gen(i,j,k)%ice,re_gen(i,j,k)%snow )

          !use general gamma assumption for re (need to make sure to Ted)
          call re_gamma_2mmt( 'qr',mu_gen%rain, rho_gen%rain, q_gen(i,j,k)%rain, n_gen(i,j,k)%rain, re_gen(i,j,k)%rain )
          call re_gamma_2mmt( 'qg',mu_gen%graupel, rho3d_graupel(i,j,k), q_gen(i,j,k)%graupel, n_gen(i,j,k)%graupel, &
                              re_gen(i,j,k)%graupel )
          call re_gamma_2mmt( 'qh',mu_gen%hail, rho3d_hail(i,j,k), q_gen(i,j,k)%hail, n_gen(i,j,k)%hail, re_gen(i,j,k)%hail )


        case('NTUCLR')

          call re_jpchen_2mmt( q_gen(i,j,k)%cloud, n_gen(i,j,k)%cloud, rho_gen%cloud, atmos(i,j,k)%rho_dair, &
                               -1.2560e+0, -1.7904e-2, 8.5536e-1 , re_gen(i,j,k)%cloud )

          call re_jpchen_2mmt( q_gen(i,j,k)%rain , n_gen(i,j,k)%rain , rho_gen%rain, atmos(i,j,k)%rho_dair, &
                               -9.9216e-2,  2.9490e-2, 9.9238e-1 , re_gen(i,j,k)%rain )

          call re_gamma_2mmt( 'qi',mu_gen%ice , rho_gen%ice , q_gen(i,j,k)%ice , n_gen(i,j,k)%ice , re_gen(i,j,k)%ice  )

          call re_gamma_1mmt(rho_gen%snow    , n0_gen%snow    , q_gen(i,j,k)%snow    , re_gen(i,j,k)%snow )

          call re_gamma_1mmt(rho_gen%graupel , n0_gen%graupel , q_gen(i,j,k)%graupel , re_gen(i,j,k)%graupel )


        case('HUCM_SBM','HUCM_SBM43') 

          call re_sbm_bin( q_sbm(i,j,k)%liq     ,    rad_sbm(1:nbin)%liq      , n_sbm(i,j,k,1:nbin)%liq &
                          ,drad_sbm(1:nbin)%liq ,re_sbm(i,j,k)%liq       )   !liquid (cloud+rain) 

          call re_sbm_bin( q_sbm (i,j,k)%ice_col,   rad_sbm(1:nbin)%ice_col  , n_sbm(i,j,k,1:nbin)%ice_col &
                          ,drad_sbm(1:nbin)%ice_col ,re_sbm(i,j,k)%ice_col   )   !ice column 

          call re_sbm_bin( q_sbm (i,j,k)%ice_pla,   rad_sbm(1:nbin)%ice_pla  , n_sbm(i,j,k,1:nbin)%ice_pla &
                          ,drad_sbm(1:nbin)%ice_pla ,re_sbm(i,j,k)%ice_pla   )   !ice plate

          call re_sbm_bin( q_sbm (i,j,k)%ice_den,   rad_sbm(1:nbin)%ice_den  , n_sbm(i,j,k,1:nbin)%ice_den &
                          ,drad_sbm(1:nbin)%ice_den ,re_sbm(i,j,k)%ice_den   )   !ice dendride

          call re_sbm_bin( q_sbm (i,j,k)%snow   ,   rad_sbm(1:nbin)%snow     , n_sbm(i,j,k,1:nbin)%snow &
                          ,drad_sbm(1:nbin)%snow ,re_sbm(i,j,k)%snow   )      !snow

          call re_sbm_bin( q_sbm (i,j,k)%graupel,   rad_sbm(1:nbin)%graupel  , n_sbm(i,j,k,1:nbin)%graupel &
                          ,drad_sbm(1:nbin)%graupel ,re_sbm(i,j,k)%graupel   )   !graupel

          call re_sbm_bin( q_sbm (i,j,k)%hail   ,   rad_sbm(1:nbin)%hail     , n_sbm(i,j,k,1:nbin)%hail    &
                          ,drad_sbm(1:nbin)%hail ,re_sbm(i,j,k)%hail   )       !hail


         case('RAMS1','RAMS2')

          !do nothing

       case default
          call stop_sdsu(myrank,'MSG re_all: There is no such cloud_microphysics')
       end select mic_select

     enddo ; enddo ; enddo

 endif RAMS

 end subroutine re_all 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine re_gamma_2mmt( hclass, mu, rho, q, nt,  re )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute effective radius for generaized gamma distribution: N(D) = No * (D**mu) * exp(-lambda*D)
!  from 2-moment microphysics. 
! 
! History:
! 01/2012  Toshi Matsui@NASA GSFC : incorporate Mei Han's add. 
! 06/2011  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 character(len=2),intent(in) :: hclass  !hydrometeor class 
 real(sdsu_fps),intent(in) :: mu  ! spectral shape parameter [-]
 real(sdsu_fps),intent(in) :: rho ! bulk particle density [kg/m3]
 real(sdsu_fps),intent(inout) :: q   ! mixing ratio [g/m3]
 real(sdsu_fps),intent(inout) :: nt  ! total number concentration [#/m3]
 real(sdsu_fps),intent(out) :: re  !effective radius [micron]

 real(sdsu_fps) :: gfac1, gfac2, gamfac  !gamma factor
 real(sdsu_fps) :: lambda !slope parameter [1/m]

 real(sdsu_fps),parameter :: dcs = 125.e-6    ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
 real(sdsu_fps),parameter :: lammaxi = 1./1.e-6
 real(sdsu_fps),parameter :: lammini = 1./(2.*dcs+100.e-6)
 real(sdsu_fps),parameter :: lammaxr = 1./20.e-6
 real(sdsu_fps),parameter :: lamminr = 1./500.e-6    !1./2800.E-6 is from WRFV3.2; lamminr = 1./500.E-6 is in WRFV3.1
 real(sdsu_fps),parameter :: lammaxs = 1./10.e-6
 real(sdsu_fps),parameter :: lammins = 1./2000.e-6
 real(sdsu_fps),parameter :: lammaxg = 1./20.e-6
 real(sdsu_fps),parameter :: lamming = 1./2000.e-6

 logical,parameter :: limitted = .false. 
 
!
! for no particles. 
!
 if( q < q_min_condensate .or. nt <= 0.) then

    re = 0.e0 ; q = 0.e0  ; nt = 0.e0

 else

!
! gamma functions
!
   call gamma_reff(mu+4.e0 ,gfac1)
   call gamma_reff(mu+1.e0 ,gfac2)
   gamfac = (gfac1/gfac2)

!
! derive slope paramter
!                          [g/m3]     [1/m3]    [m3/g]  
   lambda = (const_pi/6. * (rho*1000.) *  nt   /  q     * gamfac    )**(1./3.)  ![1/m]

!
! limits for lambda for morrison scheme
!
 IF_MORR: if(trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH') then
  select case(hclass)
  case('qi') 
    lambda = min( max( lambda, lammini ), lammaxi )
  case('qr') 
    lambda = min( max( lambda, lamminr ), lammaxr ) 
  case('qs') 
    lambda = min( max( lambda, lammins ), lammaxs )
  case('qg') 
    lambda = min( max( lambda, lamming ), lammaxg )
  end select

 endif IF_MORR

!
! derive effective radius
!
   call gamma_reff(mu+4.e0 ,gfac1)
   call gamma_reff(mu+3.e0 ,gfac2)
   gamfac = (gfac1/gfac2)

   re = 1./(2.*lambda) * gamfac * 1.e+6  ![micron] <- [m]

 endif

 return
 end subroutine re_gamma_2mmt

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine re_jpchen_2mmt( q, num, rhoe, rho_dair, a, b, c, re_out )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute effective radius from JP-Chen's microphysics, which is SBM-empirical fitting. 
! 
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(inout) :: q   ! water content [g/m3]
 real(sdsu_fps),intent(inout) :: num ! number concentration [#/m3]
 real(sdsu_fps),intent(in) :: rhoe ! particle effective density [kg/m3]

 real(sdsu_fps),intent(in) :: rho_dair ! dry air mass [g/m3]

 real(sdsu_fps),intent(in) :: a,b,c  !fitting coefficients

 real(sdsu_fps),intent(out) :: re_out ! effective radius [um]

 real(sdsu_fps) :: Nc  ! number mixing ratio  [#/kg]
 real(sdsu_fps) :: Qc  ! mass mixing ratio [kg/kg]
 real(sdsu_fps) :: rc  ! mean volume radius [m]

!
! filter bad values
!
 if( q <= q_min_condensate .or. num <= 0. ) then
    q = 0. ; num = 0. ; re_out = 0.  ; return
 endif

 Nc = num / ( rho_dair * 1.e-3 )  ! number mixing ratio  [#/kg]

 Qc = q / rho_dair  ! mass mixing ratio [kg/kg]

 rc = ( 3.*Qc / (4.*const_pi*rhoe*Nc) )**(1./3.)  !mean volume radius [m]

!
! filter bad values
! 
 if( Nc <= 0. .or. rc <= 0. ) then
    q = 0. ; num = 0. ; re_out = 0.  ; return
 endif

!
! Chen and Liu (2004) Table 3 equation 6a and 7
!
 re_out = exp( a + b * log(Nc) + c * log(rc) ) * 1.e6  ! effecive radius [um]


 return
 end subroutine re_jpchen_2mmt

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine re_sbm_bin(wc, rad, n0 , drad , re , frime )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute effective radius from SBM output. 
! 
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(inout) :: wc   !water content [g/m3]
 real(sdsu_fps),intent(in),dimension(nbin)   :: rad  !bulk radius [cm]
 real(sdsu_fps),intent(inout),dimension(nbin)   :: n0   !number concentration per bin DSD [1/m4]
 real(sdsu_fps),intent(in),dimension(nbin)   :: drad  !bin width [cm]
 real(sdsu_fps),intent(out) :: re  !effective radius [micron] 
 real(sdsu_fps),intent(in),optional :: frime(nbin) !riming fraction of snow aggregate [-]

 integer :: n  !local looping indice
 real(sdsu_fps) :: third_moment , sec_moment  !3rd and 2nd moment of droplet size distribution
 real(sdsu_fps) :: r_m   ! particle radius [m]
 real(sdsu_fps) :: dr_m  ! bin width [m]
 real(sdsu_fps) :: den_diag, rad_diag 

!
! for no particles. 
!
 if( wc < q_min_condensate ) then 
    re = 0.e0 ; wc = 0.e0 ; n0 = 0.e0
    return
 endif

!
! initialization
!
 third_moment = 0.e0 ; sec_moment = 0.e0

!
! For no particles. 
!
 if( sum(n0(:)) == 0.e0) then
    re = 0.e0 ; wc = 0.e0 
    return
 endif

!
! n0(n) * (rad(n)/100./3.) is tricky but based on the relationship 
! G(m)/X = n0(r)*r/3  (SBM scheme)
! where G(m) is [g/m3], X is [g], n0(r) is [#/m4]
!
 do n = 1, nbin

    ! third_moment = third_moment + ( (rad(n)/100.e0)**3) * n0(n) * (rad(n)/100.e0/3.) ![m3]  
    ! sec_moment = sec_moment + ( (rad(n)/100.e0)**2) * n0(n) * (rad(n)/100.e0/3.)     ![m2]

    if( present(frime) ) then
      call diag_rime_snow (n, frime(n), den_sbm(1:nbin)%snow, &
                           x_sbm(1:nbin)%snow, rad_sbm(1:nbin)%snow, &
                           den_diag, rad_diag )
      r_m = rad_diag * 1.e-2
    else
      r_m = rad(n) * 1.e-2  ! particle radius [m]
    endif

    dr_m = drad(n) * 1.e-2  ! bin width [m]
    third_moment = third_moment + ( r_m**3) * n0(n) * dr_m ![m3/m3]  
    sec_moment   = sec_moment   + ( r_m**2) * n0(n) * dr_m ![m2/m3]

 enddo

  re = third_moment / sec_moment *1.0e+6   ! [micron] <- [m]

 return
 end subroutine re_sbm_bin

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine re_gamma_1mmt(den, n0 , wc, re )
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute drop effective radius for one-moment bulk scheme using analytic solution (fast). 
!  
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(in) :: den  !density [kg/m3]
 real(sdsu_fps),intent(in) :: n0   !intercept of exponential DSD [1/m4]
 real(sdsu_fps),intent(inout) :: wc   !liquid/ice water content [g/m3] 
 real(sdsu_fps),intent(out) :: re  !effective radius [micron] 
 real(sdsu_fps) :: lam   !intercept [1/m]

!
! for no particles. 
!
 if( wc < q_min_condensate ) then
    re = 0.e0 ; wc = 0.e0
 else

!
! compute drop effective radius for exponential distribution N(D) = N0*exp(-lam*D)
!

 lam = (n0*const_pi*den/(wc*(1.e-3)))**(0.25e0)  ![1/m]
 
 re = 3.e0 / 2.e0 / lam *1.0e+6  ![micron] <- [m]


!
! gamma functions
!
!   call gamma_reff(mu+4.e0 ,gfac1)
!   call gamma_reff(mu+1.e0 ,gfac2)
!   gamfac = (gfac1/gfac2)
!
! derive intercept
!       [g/m3]                     [m3/g]        [1/m]
!   lam = (n0 * gfac1 * (dens*1000.) * cost_pi /( wc * 6.) ) ** (-4.-mu)  ! [1/m**(4+mu)]
!
! derive effective radius
!
!   call gamma_reff(mu+4.e0 ,gfac1)
!   call gamma_reff(mu+3.e0 ,gfac2)
!   gamfac = (gfac1/gfac2)
!   re = 1./(2.*lambda) * gamfac * 1.e+6  ![micron] <- [m]

 endif

 return
 end subroutine re_gamma_1mmt

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine re_bulk_tedd(temp, wc, Tctl, a, b, re_max, re) 

 implicit none
!----------------------------------------------------------------------
! Purpose: Compute new intercept from the empirical relationship between re and temperature
!
! History:
! 12/2007  Toshi Matsui@NASA GSFC ; Initial
!----------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: temp   ! temperature [K]
 real(sdsu_fps),intent(in) :: wc     ! water content [g/m3]
 real(sdsu_fps),intent(in) :: Tctl   ! control temperature [K]
 real(sdsu_fps),intent(in) :: a,b    ! tuning parameter [-]
 real(sdsu_fps),intent(in) :: re_max ! maximum drop effectve radius [micon]
 real(sdsu_fps),intent(out):: re     ! drplet effective radius [micron]

!
! for no particles. 
!
 if( wc == 0.) then
    re = 0.
    return
 endif

!
!
!
 if(temp > Tctl) then !below melting level
   re = re_max ![micron]
 else  !above melting level
   re = min(  re_max ,  a/(Tctl-temp)+b  )  ![micron]
 endif


 return
 end subroutine re_bulk_tedd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine re_bulk_tedd_lut(spc, qcol,tair, wc, re )
 implicit none

!--------------------------------------------------------------------------------------------------
! Comments: 
!  This program compute droplet effective radius as a function of air temperature and 
!  condensates water path of each condensate species. 
!  LUT is derived from the result of the GCE with spectra-bin microphysics (by X. Li). 
!   
! History:
! 07/2008  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: spc      ! spicies index: 1-snow, 2-graupel, 3-hail
 real(sdsu_fps),intent(in)    :: qcol     ! water path [kg/m2]
 real(sdsu_fps),intent(in)    :: tair     ! air temperature [K]
 real(sdsu_fps) :: re_sbm_temp_mass !effective radius [mm]
 real(sdsu_fps) :: w1,w2 !weights for interpolation
 integer :: k,m  !looping indices
!
! Look-Up Table parameter
!
 integer,parameter :: spcmax= 3 !maximum # of spicies
 integer,parameter :: kmax=18 !maximum # of temperature bin
 integer,parameter :: mmax=10 !maximum # of mass bin
  real(sdsu_fps) :: temp_sbm(spcmax,kmax)   ! temp bin for re LUT [K]
  real(sdsu_fps) :: mass_sbm(spcmax,mmax)   ! mass bin for re LUT [kg/m2]
  real(sdsu_fps) :: re_sbm(spcmax,kmax,mmax)! re [mm] LUT as a functino of temperature (18) and mass (10)
  real(sdsu_fps) :: re_sbm_temp(mmax) !re_sbm interpolated for temperature [mm]

 real(sdsu_fps),intent(in) :: wc   !liquid/ice water content [g/m3] 
 real(sdsu_fps),intent(out) :: re  !effective radius [micron] 

 ! 
 ! Snow parameter
 ! 
 data (temp_sbm(1,k),k=1,kmax)/ &
  283.62,&
  282.31,&
  280.36,&
  277.48,&
  275.54,&
  271.08,&
  267.65,&
  264.69,&
  260.80,&
  256.91,&
  252.55,&
  247.15,&
  241.14,&
  234.43,&
  226.27,&
  218.06,&
  209.73,&
  201.95/
 data (mass_sbm(1,m),m=1,mmax)/ &
 .5950E+00,.1585E+01,.2575E+01,.3565E+01,.4555E+01,.5545E+01,.6535E+01,.7525E+01,.8515E+01,.9505E+01/
 data ((re_sbm(1,k,m),m=1,mmax),k=1,kmax)/ &
 .8897E+01,.7589E+01,.7346E+01,.7129E+01,.7528E+01,.8031E+01,.8326E+01,.8572E+01,.8542E+01,.8983E+01,&
 .2734E+01,.3041E+01,.3023E+01,.3525E+01,.4743E+01,.6134E+01,.6422E+01,.6422E+01,.6633E+01,.6773E+01,&
 .1242E+01,.8845E+00,.1053E+01,.1782E+01,.2799E+01,.3660E+01,.3829E+01,.3861E+01,.3812E+01,.3860E+01,&
 .5873E+00,.4417E+00,.5955E+00,.7478E+00,.1124E+01,.1626E+01,.1752E+01,.1799E+01,.1721E+01,.1684E+01,&
 .3777E+00,.3577E+00,.4227E+00,.4738E+00,.6305E+00,.9038E+00,.1065E+01,.1164E+01,.1176E+01,.1173E+01,&
 .2845E+00,.3246E+00,.3786E+00,.4278E+00,.5940E+00,.8226E+00,.9036E+00,.9633E+00,.9504E+00,.9563E+00,&
 .2654E+00,.2994E+00,.3081E+00,.3196E+00,.3729E+00,.4591E+00,.5243E+00,.6116E+00,.6436E+00,.6764E+00,&
 .2582E+00,.2708E+00,.2702E+00,.2785E+00,.3076E+00,.3478E+00,.3920E+00,.4534E+00,.4799E+00,.5197E+00,&
 .2446E+00,.2478E+00,.2486E+00,.2573E+00,.2710E+00,.2981E+00,.3279E+00,.3628E+00,.3848E+00,.4151E+00,&
 .2268E+00,.2343E+00,.2367E+00,.2379E+00,.2405E+00,.2608E+00,.2796E+00,.3010E+00,.3146E+00,.3235E+00,&
 .2036E+00,.2276E+00,.2291E+00,.2227E+00,.2197E+00,.2310E+00,.2436E+00,.2576E+00,.2658E+00,.2679E+00,&
 .1787E+00,.2286E+00,.2257E+00,.2134E+00,.2050E+00,.2099E+00,.2224E+00,.2336E+00,.2377E+00,.2381E+00,&
 .1590E+00,.2257E+00,.2193E+00,.2046E+00,.1917E+00,.1920E+00,.2026E+00,.2144E+00,.2184E+00,.2162E+00,&
 .1483E+00,.2145E+00,.2087E+00,.1944E+00,.1789E+00,.1746E+00,.1817E+00,.1908E+00,.1999E+00,.2073E+00,&
 .1481E+00,.2040E+00,.1961E+00,.1824E+00,.1672E+00,.1611E+00,.1658E+00,.1758E+00,.1801E+00,.1985E+00,&
 .1379E+00,.1880E+00,.1765E+00,.1642E+00,.1525E+00,.1461E+00,.1495E+00,.1564E+00,.1582E+00,.1669E+00,&
 .1200E+00,.1631E+00,.1506E+00,.1390E+00,.1305E+00,.1245E+00,.1274E+00,.1326E+00,.1330E+00,.1338E+00,&
 .9719E-01,.1240E+00,.1140E+00,.1037E+00,.9496E-01,.8972E-01,.9407E-01,.1017E+00,.1063E+00,.1051E+00/


 ! 
 ! Graupel parameter
 ! 
 data (temp_sbm(2,k),k=1,kmax)/ &
  284.16,&
  281.92,&
  279.69,&
  277.60,&
  275.46,&
  271.06,&
  267.52,&
  264.59,&
  260.81,&
  257.02,&
  252.72,&
  247.34,&
  241.27,&
  234.49,&
  226.33,&
  218.11,&
  209.76,&
  201.98/
 data (mass_sbm(2,m),m=1,mmax)/ &
 .1845E+01,.5335E+01,.8825E+01,.1232E+02,.1581E+02,.1930E+02,.2279E+02,.2628E+02,.2977E+02,.3326E+02/
 data ((re_sbm(2,k,m),m=1,mmax),k=1,kmax)/ &
 .3886E+01,.3268E+01,.3324E+01,.2924E+01,.2774E+01,.2779E+01,.2675E+01,.2613E+01,.2572E+01,.2615E+01,&
 .2556E+01,.2200E+01,.2382E+01,.2273E+01,.2250E+01,.2130E+01,.1956E+01,.1920E+01,.1802E+01,.1643E+01,&
 .1076E+01,.1388E+01,.1612E+01,.1576E+01,.1519E+01,.1336E+01,.1115E+01,.1100E+01,.9941E+00,.9445E+00,&
 .3444E+00,.4616E+00,.6489E+00,.7234E+00,.6757E+00,.4937E+00,.4137E+00,.4407E+00,.4614E+00,.4799E+00,&
 .2214E+00,.2303E+00,.3031E+00,.3352E+00,.4045E+00,.3663E+00,.3636E+00,.3865E+00,.4225E+00,.4409E+00,&
 .2071E+00,.1977E+00,.2412E+00,.2708E+00,.3392E+00,.3346E+00,.3608E+00,.3815E+00,.4207E+00,.4437E+00,&
 .1913E+00,.1856E+00,.2258E+00,.2601E+00,.3197E+00,.3230E+00,.3601E+00,.3773E+00,.4164E+00,.4427E+00,&
 .1834E+00,.1792E+00,.2160E+00,.2506E+00,.3065E+00,.3106E+00,.3566E+00,.3713E+00,.4040E+00,.4299E+00,&
 .1772E+00,.1754E+00,.2101E+00,.2446E+00,.2960E+00,.3017E+00,.3519E+00,.3621E+00,.3872E+00,.4141E+00,&
 .1717E+00,.1719E+00,.2051E+00,.2399E+00,.2889E+00,.2988E+00,.3470E+00,.3517E+00,.3668E+00,.3861E+00,&
 .1657E+00,.1677E+00,.2007E+00,.2372E+00,.2836E+00,.2945E+00,.3435E+00,.3453E+00,.3477E+00,.3457E+00,&
 .1611E+00,.1633E+00,.1953E+00,.2340E+00,.2784E+00,.2917E+00,.3422E+00,.3425E+00,.3329E+00,.3284E+00,&
 .1579E+00,.1583E+00,.1894E+00,.2294E+00,.2728E+00,.2891E+00,.3383E+00,.3346E+00,.3199E+00,.3251E+00,&
 .1546E+00,.1534E+00,.1825E+00,.2220E+00,.2658E+00,.2841E+00,.3324E+00,.3295E+00,.3157E+00,.3255E+00,&
 .1499E+00,.1483E+00,.1744E+00,.2103E+00,.2574E+00,.2755E+00,.3237E+00,.3251E+00,.3120E+00,.3202E+00,&
 .1433E+00,.1424E+00,.1656E+00,.1970E+00,.2432E+00,.2565E+00,.3112E+00,.3202E+00,.3103E+00,.3012E+00,&
 .1343E+00,.1353E+00,.1555E+00,.1822E+00,.2257E+00,.2352E+00,.2861E+00,.3052E+00,.3080E+00,.2954E+00,&
 .1216E+00,.1241E+00,.1424E+00,.1663E+00,.2049E+00,.2122E+00,.2682E+00,.2738E+00,.2918E+00,.3000E+00/


 ! 
 ! Hail parameter
 ! 
 data (temp_sbm(3,k),k=1,kmax)/ &
  284.03,&
  281.83,&
  279.43,&
  277.14,&
  275.14,&
  270.66,&
  267.35,&
  264.66,&
  261.13,&
  257.43,&
  253.08,&
  247.59,&
  241.41,&
  234.52,&
  226.40,&
  218.16,&
  209.85,&
  202.09/
 data (mass_sbm(3,m),m=1,mmax)/ &
 .7950E+00,.2185E+01,.3575E+01,.4965E+01,.6355E+01,.7745E+01,.9135E+01,.1052E+02,.1192E+02,.1331E+02/
 data ((re_sbm(3,k,m),m=1,mmax),k=1,kmax)/ &
 .2300E+01,.1954E+01,.2162E+01,.2250E+01,.2304E+01,.2210E+01,.2112E+01,.2068E+01,.2063E+01,.1987E+01,&
 .1221E+01,.1440E+01,.1590E+01,.1692E+01,.1726E+01,.1683E+01,.1627E+01,.1576E+01,.1665E+01,.1556E+01,&
 .8192E+00,.9854E+00,.1090E+01,.1163E+01,.1184E+01,.1170E+01,.1133E+01,.1145E+01,.1202E+01,.1062E+01,&
 .2741E+00,.4490E+00,.5047E+00,.5507E+00,.5486E+00,.5927E+00,.6179E+00,.6235E+00,.6672E+00,.5525E+00,&
 .1660E+00,.2362E+00,.2780E+00,.3315E+00,.3481E+00,.3931E+00,.4640E+00,.4885E+00,.5098E+00,.4615E+00,&
 .1506E+00,.2043E+00,.2291E+00,.2692E+00,.2958E+00,.3393E+00,.4171E+00,.4532E+00,.4341E+00,.4156E+00,&
 .1400E+00,.1853E+00,.2120E+00,.2469E+00,.2799E+00,.3273E+00,.3958E+00,.4422E+00,.4187E+00,.4272E+00,&
 .1348E+00,.1761E+00,.2007E+00,.2314E+00,.2693E+00,.3222E+00,.3930E+00,.4418E+00,.4189E+00,.4350E+00,&
 .1294E+00,.1688E+00,.1897E+00,.2155E+00,.2555E+00,.3079E+00,.3892E+00,.4408E+00,.4182E+00,.4444E+00,&
 .1243E+00,.1600E+00,.1779E+00,.2015E+00,.2418E+00,.2935E+00,.3819E+00,.4387E+00,.4155E+00,.4491E+00,&
 .1206E+00,.1512E+00,.1670E+00,.1868E+00,.2274E+00,.2821E+00,.3712E+00,.4321E+00,.4094E+00,.4509E+00,&
 .1174E+00,.1429E+00,.1563E+00,.1728E+00,.2102E+00,.2699E+00,.3574E+00,.4201E+00,.3971E+00,.4492E+00,&
 .1150E+00,.1363E+00,.1466E+00,.1602E+00,.1961E+00,.2535E+00,.3423E+00,.4034E+00,.3777E+00,.4336E+00,&
 .1125E+00,.1305E+00,.1377E+00,.1492E+00,.1842E+00,.2379E+00,.3256E+00,.3813E+00,.3552E+00,.4067E+00,&
 .1093E+00,.1246E+00,.1288E+00,.1402E+00,.1702E+00,.2222E+00,.3062E+00,.3605E+00,.3307E+00,.3756E+00,&
 .1054E+00,.1180E+00,.1194E+00,.1302E+00,.1551E+00,.2030E+00,.2790E+00,.3380E+00,.3116E+00,.3440E+00,&
 .1009E+00,.1112E+00,.1110E+00,.1199E+00,.1405E+00,.1827E+00,.2472E+00,.3016E+00,.2837E+00,.3131E+00,&
 .9614E-01,.1037E+00,.1013E+00,.1111E+00,.1280E+00,.1593E+00,.2162E+00,.2651E+00,.2500E+00,.2811E+00/
!------------------------------------------------------------------------------------------------------

!
! for no particles. 
!
 if( wc == 0.) then
    re = 0.
    return
 endif

!
! interpolate for temperature
!
  if( tair <= temp_sbm(spc,kmax) ) then !too cold 
    re_sbm_temp(1:mmax) = re_sbm(spc,kmax,1:mmax)
  elseif( tair >= temp_sbm(spc,1) ) then !too hot
    re_sbm_temp(1:mmax) = re_sbm(spc,1,1:mmax)
  else  !linear interpolation

    do k = 1, kmax-1
       if( tair <= temp_sbm(spc,k) .and. tair > temp_sbm(spc,k+1) ) then
           w1 = ( tair - temp_sbm(spc,k+1) ) / ( temp_sbm(spc,k) - temp_sbm(spc,k+1) )
           w2 = 1.-w1
           re_sbm_temp(1:mmax) = w1*re_sbm(3,k,1:mmax) + w2*re_sbm(3,k+1,1:mmax)
           exit
       endif
    enddo

  endif

!
! intepolate for mass
!
 if(qcol <= mass_sbm(spc,1) ) then         !too small
    re_sbm_temp_mass = re_sbm_temp(1)
 elseif( qcol >= mass_sbm(spc,mmax) ) then !too large 
    re_sbm_temp_mass = re_sbm_temp(mmax)
 else
    do m = 1, mmax-1
       if( qcol >= mass_sbm(spc,m) .and. qcol < mass_sbm(spc,m+1) )then
           w1 = ( mass_sbm(spc,m+1) - qcol) / (mass_sbm(spc,m+1)-mass_sbm(spc,m))
           w2 = 1.- w1
           re_sbm_temp_mass = w1*re_sbm_temp(m) + w2*re_sbm_temp(m+1)
           exit
       endif

    enddo
 endif

!
! unit conversion
!
 re = re_sbm_temp_mass *1e+3 ![micron] <- [mm]

 return
 end subroutine re_bulk_tedd_lut

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine re_LUT_Heymsfield_Platt_1984(proc,temp,wc, re)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute drop effective radius from observed DSD in Heymsfield Platt 1984]  
! And prepared it as LUT for faster computation. 
! 
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!
!----------------------------------------------------------------------------------------------------

  character(len=4),intent(in) :: proc !'init' or 'proc'
  real(sdsu_fps),intent(in)  :: temp  !temperature [K]
  real(sdsu_fps),intent(in)  :: wc    ! water content [g/m3]
  real(sdsu_fps),intent(out) :: re    !effective radius [micron]

  integer :: t ! loop indice
  integer,parameter :: tmin = 210 , tmax = 260
  real(sdsu_fps),save :: re_lut(tmin:tmax)  !LUT of re as a function of temperature [micron]
  real(sdsu_fps) :: wgt1,wgt2 !weighting function


  proc_select: select case(proc)
  case('init')

!
! creating LUT of re 
!
    do t = tmin, tmax
       call re_Heymsfield_Platt_1984(REAL(t),1., re)
       re_lut(t) = re
    enddo !t

  case('proc')

   if(wc == 0.) then
      re = 0.
      return
   endif
!
! interpolate re from LUT
!
   if(temp <= REAL(tmin) ) then  !temp is colder than tmin
      re = re_lut(tmin) !micron 
   elseif(temp >= REAL(tmax) ) then !temp is warmer than tmax
      re = re_lut(tmax)
   else !interpolate
      wgt2 = temp - int(temp)
      wgt1 = 1.-wgt2
      re  = wgt1*re_lut(int(temp)) + wgt2*re_lut(int(temp)+1)
   endif

  case default
       call stop_sdsu(myrank,'MSG re_LUT_Heymsfield_Platt_1984: the option does not exist: proc ')
  end select proc_select


 return
 end subroutine re_LUT_Heymsfield_Platt_1984

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine re_Heymsfield_Platt_1984(temp,lwc, re)
!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute drop effective radius from observed DSD in Heymsfield Platt 1984]  
!  
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; initial
!           
! References: 
!   Heymsfield, A.J., and C. Platt, 1984: A Parameterization of the Particle Size Spectrum of Ice 
!     Clouds in Terms of the Ambient Temperature and the Ice Water Content. J. Atmos. Sci., 41, 846-855.
!
!----------------------------------------------------------------------------------------------------

  implicit none
  real(sdsu_fps),intent(in)  :: temp  !temperature [K]
  real(sdsu_fps),intent(in)  :: lwc   !liquid water content [g/m3]
  real(sdsu_fps),intent(out) :: re    !effective radius [micron]

  integer :: imax
  integer :: i
  integer :: nopt
  real(sdsu_fps) :: rad
  real(sdsu_fps) :: densice
  real(sdsu_fps) :: densliq
  real(sdsu_fps) :: density
  real(sdsu_fps) :: densi
  real(sdsu_fps) :: iwctest
  real(sdsu_fps) :: norm
  real(sdsu_fps) :: num
  real(sdsu_fps) :: faa
  real(sdsu_fps) :: dr

  data densliq/1.0e+3/
  data densice/0.917e+3/
  real(sdsu_fps) :: third_moment  ![m3]  
  real(sdsu_fps) :: sec_moment     ![m2]

  third_moment = 0.  ![m3]  
  sec_moment    =0.  ![m2]

!
!     Begin by checking if hydrometeors of this species are present.
!     If not, set scattering parameters to zero and return.
!
  if(lwc .lt. q_min_condensate) then
     return
  endif

!
!     Loop over particle sizes:

!     increments of particle radius are 0.005 mm; the particle
!     size distribution is expressed as a particle number density,
!     num, per radius increment.  This distribution is given
!     by the fit to observed cloud ice distributions by
!     Heymsfield and Platt (1984)
!     two options are available:
!     nopt=0:   ice particle mass distributed in spherical
!               volume with diameter equal to maximum crystal
!               dimension (l).
!     nopt=1:   ice particle described as pure ice sphere
!               with same mass as elongated crystal.
  nopt=1
  dr = 0.005e0
  imax = nint(2.5e0/dr)
!
!     first compute normalization factor for particle size distribution
  norm=1.
  iwctest=0.

 do i=0,imax
    rad=dr*0.5e0+dr*float(i)
    call heymplatt(nopt,temp,lwc,rad,densi,num,norm)  !num [1/m**4]
    density=densi
    iwctest=iwctest+ &
    num*(1.e-3)*density*4.e0*const_pi*((rad*1.e0-1.e0)**3.)/3.e0*.005e0
 enddo

 norm=iwctest/lwc
 do i=0,imax
    rad=dr*0.5e0+dr*float(i)

    call heymplatt(nopt,temp,lwc,rad,densi,num,norm) !num [1/m**4]

!
! compute 2nd 3rd moment
!
    third_moment = third_moment + ( rad**3) * num *dr ! [mm3]  
    sec_moment   = sec_moment   + ( rad**2) * num *dr ! [mm2]

 end do

!
! compute effective radius
!
 re = third_moment / sec_moment *1e+3


 return
 end subroutine re_Heymsfield_Platt_1984

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine heymplatt(nopt,t,iwc,r,den,numden,norm)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Return particle density and number density
! of ice crystals, given temperature, ice water content of
! ice crystal distribution, and maximum particle
! dimension.  Follows the empirical relation of
! Heymsfield and Platt (1984).
!
! Input
! Two options are available:
! nopt=0:   ice particle mass distributed in spherical
!             volume with diameter equal to maximum crystal
! nopt=1:      ice particle described as pure ice sphere
!            with same mass as elongated crystal.
! t            temperature [K]
! iwc      equivalent water content of cloud ice [g/m**3]
! r            particle radius [mm]
! norm     normalization factor required to contrain
!            integrated distribution to equal the equivalent
!            ice water content.  heymplatt must be called
!            first with norm set equal to 1 to give
!           proper normalization
!
! Output
! den      particle density [g/cm**3]
! numden   particle number density [1/m**4]
!
!
! History:
! 03/2008  Toshi Matsui@NASA GSFC : Addapted to SDSU (contribution from Bill Olson@NASA GSFC)
!
! References: 
!   Heymsfield, A.J., and C. Platt, 1984: A Parameterization of the Particle Size Spectrum of Ice 
!     Clouds in Terms of the Ambient Temperature and the Ice Water Content. J. Atmos. Sci., 41, 846-855.
!---------------------------------------------------------------------------------------------

 integer :: n
 integer :: nopt
 real(sdsu_fps):: t
 real(sdsu_fps):: iwc
 real(sdsu_fps):: r
 real(sdsu_fps):: den
 real(sdsu_fps):: numden
 real(sdsu_fps):: a(0:4)
 real(sdsu_fps):: b(0:4)
 real(sdsu_fps):: c(0:4)
 real(sdsu_fps):: d(0:4)
 real(sdsu_fps):: e(0:4)
 real(sdsu_fps):: suma
 real(sdsu_fps):: sumb
 real(sdsu_fps):: sumc
 real(sdsu_fps):: sumd
 real(sdsu_fps):: sume
 real(sdsu_fps):: tc
 real(sdsu_fps):: mass
 real(sdsu_fps):: l
 real(sdsu_fps):: l0
 real(sdsu_fps):: b1
 real(sdsu_fps):: b2
 real(sdsu_fps):: a1
 real(sdsu_fps):: a2
 real(sdsu_fps):: n100diwc
 real(sdsu_fps):: n1000diwc
 real(sdsu_fps):: lc
 real(sdsu_fps):: dc
 real(sdsu_fps):: densice
 real(sdsu_fps):: norm
 real(sdsu_fps):: dldr

 data a/-1.1430e+1,-7.3892e-1,-1.8647e-2,-1.4045e-4,0.0e+0/
 data b/1.8940e+1,2.7658e0,1.2833e-1,2.7750e-3,2.2994e-5/
 data c/-1.0159e+1,-1.4538e0,-1.3511e-2,1.1318e-3,2.2360e-5/
 data d/1.6764e+1,-1.5072e-1,-1.9713e-2,-3.5051e-4,-1.6727e-6/
 data e/1.5508e+2,1.8377e+1,8.5312e-1,1.6879e-2,1.1873e-4/
 data densice/0.917/


! Temperature [C]
  tc=t-const_Kel2Cel

! Note: empirical formulae only apply to range
! -60 < tc < -20
  if(tc .gt. -20.) tc=-20.
  if(tc .lt. -60.) tc=-60.

! Calculate maximum particle dimension [microns]
  if(nopt .eq. 0) then
     l=2.*r*1.e+3
  else
     mass=densice*4.*const_pi*((r*1.e-1)**3.)/3.
    if(2.*r .le. 0.3) then
       l=((mass/(densice*(3.*sqrt(3.)/2.)* &
         ((0.5/2.)**2.)*1.e-3))**.33333)* &
        1.e+3
       dldr=(4./3.)*(1.e-3)*densice*const_pi*(mass**(-.66666))*r*r/ &
            ((densice*(3.*sqrt(3.)/2.)*((0.5/2.)**2.)*1.e-3)**.33333)
    else
       l=((mass/(densice*(3.*sqrt(3.)/2.)*((.2/2.)**2.)*1.e-3)) &
         **.55249)*1.e+3
       dldr=(4./1.82)*(1.e-3)*densice*const_pi*(mass**(1./1.82-1.))*r*r/ &
            ((densice*(3.*sqrt(3.)/2.)*((0.2/2.)**2.)*1.e-3)**(1./1.82))
    endif
 endif

 suma=0.
 sumb=0.
 sumc=0.
 sumd=0.
 sume=0.
 do n=0,4
    suma=suma+a(n)*(tc**n)
    sumb=sumb+b(n)*(tc**n)
    sumc=sumc+c(n)*(tc**n)
    sumd=sumd+d(n)*(tc**n)
    sume=sume+e(n)*(tc**n)
 enddo
 b1=suma
 b2=sumb

!  Liou's fit of b2 fails at low temperature;
!  since Heymsfield and Platt (1984) data indicate
!  a nearly constant value of -4., we use it here.
 b2=-4.

 if(tc .ge. -37.5) then
    n100diwc=exp(sumc)
 else
    n100diwc=exp(sumd)
 endif
 n1000diwc=sume

 a1=n100diwc/(100.**b1)
 a2=n1000diwc/(1000.**b2)
 l0=(a2/a1)**(1./(b1-b2))

 if(nopt .eq. 0) then
    if(l .le. l0) then
       numden=(1.e+6)*2.*a1*(l**b1)*iwc/norm
    else
       numden=(1.e+6)*2.*a2*(l**b2)*iwc/norm
    endif
 else
    if(l .le. l0) then
       numden=(1.e+6)*dldr*a1*(l**b1)*iwc/norm
    else
       numden=(1.e+6)*dldr*a2*(l**b2)*iwc/norm
    endif
 endif

! Compute particle density
! Assume randomly oriented hexagonal column with
! width/length relationship determined
! from Heymsfield's empirical relation

 if(nopt .eq. 0) then
    if(l*1.e-3 .le. 0.3) then
       dc=0.5*(l*1.e-3)*1.e-1
       lc=l*1.e-4
    else
       dc=0.2*((l*1.e-3)**0.41)*1.e-1
       lc=l*1.e-4
    endif
       mass=densice*(3.*sqrt(3.)/2.)*((dc/2.)**2.)*lc
       den=mass/(4.*const_pi*((r*1.e-1)**3.)/3.)
 else
    den=0.917
 end if

 return
 end subroutine heymplatt

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine re_bulk_god10_snow(tair, rho, q_snow, re)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!   Compute slope of exponential PSD as a function of temperature and mixing ratio for snow.
!  
! History:
! 03/2009  Toshi Matsui@NASA GSFC ; Adapted from GCE. 
!           
! References: 
!
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: tair    ! air temperature [k]
 real(sdsu_fps),intent(in) :: rho     ! density [kg/m3]
 real(sdsu_fps),intent(in) :: q_snow  ! snow mass mixing ratio [g/m3]
 real(sdsu_fps),intent(out) :: re     ! effective radius [micron]

!local parameters
 real(sdsu_fps) :: slope  ! slope of exponential PSD [1/cm]
 real(sdsu_fps) :: roqs 
 real(sdsu_fps) :: tslopes, cpi, dsnomin
 real(sdsu_fps) :: xs, sno11, sno00, dsno11, dsno00, sexp11, sexp00
 real(sdsu_fps) :: slim11, slim00, slit1, slit0, slexp, stt, stexp, sbase 
 real(sdsu_fps) :: a2, tairc, tairc1 ,ftnsT, sno1, dsno1, sexp1, slim1, tairc11
 real(sdsu_fps) :: sdl,  xx
 real(sdsu_fps) :: fexp, fftnsT, ftnsQ, fftnsTQ, tnsmax
 real(sdsu_fps) :: lambda   ! lambda [1/cm]

 real,parameter :: tns = 0.10e0  ! default snow intercept [cm-4]

!
! too small mixing ratio -> return
!
 if( q_snow < q_min_condensate ) then
    re = 0.e0 ; return 
 endif 

 roqs= rho * 1.e-3  !Snow density [g/cm3] <- [kg/m3]

! roqs=.05    !SNOW DENSITY (G/CM^3)

 tslopes=0.03218876   ! increase tns by  5 from  0 to -50C
 cpi=4.*atan(1.)
 dsnomin=0.0185 !minimum snow diameter (cm)
 dsnomin=dsnomin**4.*roqs*cpi

!
! tunable parameters
!
 xs=0.97
 sno11=0.40
 sno00=0.10
 dsno11=0.95
 dsno00=0.50
 sexp11=0.9
 sexp00=0.9
 slim11=1.0
 slim00=1.0
 slit1=-60.
 slit0=-20.
 slexp=3.3
 stt=-20.
 stexp=0.5
 sbase=0.00167


 a2= q_snow  ! [g/m3]

 fftnstq=1.0

 tairc=min(0.,max(-50.,tair-273.16)+0.0)
 tairc1=tair-273.16

 ftnsT=exp(-1.*tslopes*tairc)
 sno1=sno11
 dsno1=dsno11
 sexp1=sexp11
 slim1=slim11
 if(tairc1.gt.slit1.and.tairc1.lt.slit0)then
    tairc11=tairc1-slit0
    sdl=slit1-slit0
    slim1=slim00-(slim00-slim11)*(tairc11/sdl)**slexp
 endif
 if(tairc1.ge.slit0) slim1=slim00
 if(tairc.gt.stt)then
     sno1=sno00-(sno00-sno11)*(tairc/(stt))**stexp
     dsno1=dsno00-(dsno00-dsno11)*(tairc/(stt))**stexp
     sexp1=sexp00-(sexp00-sexp11)*(tairc/(stt))**stexp
 endif
 xx=xs-xs*min(slim1,max(0.0,(a2-sno1)/dsno1)**sexp1)


 ftnsT=ftnsT**xx
 fftnsT=ftnsT*tns
 fexp=xx
 ftnsQ=1.0
 if(a2.gt.0.) ftnsQ=(a2/sbase)**fexp
     fftnsTQ=fftnsT*ftnsQ
     tnsmax=(q_snow*1.e-6)/dsnomin   ! intercept [1/cm4]
!             [g/cm3]     / [g cm]

 if(fftnsTQ.gt.tnsmax.and.a2.gt.1.e-20) fftnsTQ=tnsmax  ! intercept [1/cm4]

 lambda = ( fftnsTQ * const_pi * roqs / (q_snow*1.e-6) ) ** 0.25e0  ! lambda [1/cm]

 re = 3.e0 / 2.e0 / lambda * 1.e+4   !effective radius [micron]

 return
 end subroutine re_bulk_god10_snow

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine re_bulk_god10_graupel(tair, rho, q_graupel,re)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!   Compute slope of exponential PSD as a function of temperature and mixing ratio for graupel.
!  
! History:
! 03/2009  Toshi Matsui@NASA GSFC ; Adapted from GCE. 
!           
! References: 
!
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: tair      ! air temperature [k]
 real(sdsu_fps),intent(in) :: rho       ! density [kg/m3]
 real(sdsu_fps),intent(in) :: q_graupel ! graupel mass mixing ratio [g/m3]
 real(sdsu_fps),intent(out) :: re       ! effective radius [micron]

!local parameters
 real(sdsu_fps) :: slope    ! slope of exponential PSD [1/cm]
 real(sdsu_fps) :: roqg
 real(sdsu_fps) :: tslopeg, cpi, dgrpmin, xg, grp11, grp00, dgrp11, dgrp00, gexp11, gexp00
 real(sdsu_fps) :: gtt, gtexp, gbase, a3, fftngtq
 real(sdsu_fps) :: tairc, tairc1, ftngT, grp1, dgrp1, gexp1, xx
 real(sdsu_fps) :: fexp, fftngT, ftngQ, tngmax
 real(sdsu_fps),parameter :: tng = 0.04e0  ! default snow intercept [cm-4]
 real(sdsu_fps) :: lambda   ! lambda [1/cm]

!
! too small mixing ratio -> return
!
 if( q_graupel < q_min_condensate ) then
    re = 0.e0 ; return
 endif

 roqg= rho * 1.e-3  !Graupel density [g/cm3] <- [kg/m3]
! roqg=.2     !graupel DENSITY (G/CM^3)

 tslopeg=0.02197225   ! increase tng by  3 from  0 to -50C

 cpi=4.*atan(1.)
 dgrpmin=0.0225 !minimum graupel diameter (cm)
 dgrpmin=dgrpmin**4.*roqg*cpi

!
! tunable parameters
!
 xg=0.93
 grp11=0.55
 grp00=0.45
 dgrp11=3.90
 dgrp00=0.25
 gexp11=0.5
 gexp00=0.7
 gtt=-20.
 gtexp=0.50
 gbase=0.0095

 a3=q_graupel  ! [g/m3]
 fftngtq=1.0

 tairc=min(0.,max(-50.,tair-273.16)+0.0)
 tairc1=tair-273.16

 tairc=min(0.,max(-50.,tair-273.16)+0.0)
 ftngT=exp(-1.*tslopeg*tairc)
 grp1=grp11
 dgrp1=dgrp11
 gexp1=gexp11
 if(tairc.gt.gtt)then
    grp1=grp00-(grp00-grp11)*(tairc/(gtt))**gtexp
    dgrp1=dgrp00-(dgrp00-dgrp11)*(tairc/(gtt))**gtexp
    gexp1=gexp00-(gexp00-gexp11)*(tairc/(gtt))**gtexp
 endif
 xx=xg-xg*min(1.,max(0.0,(a3-grp1)/dgrp1)**gexp1)
 ftngT=ftngT**xx
 fftngT=ftngT*tng
 fexp=xx
 ftngQ=1.0
 if(a3.gt.0.) ftngQ=(a3/gbase)**fexp
 fftngTQ=fftngT*ftngQ

 tngmax=(q_graupel*1.e-6)/dgrpmin  ! intercept [1/cm4]

 if(fftngTQ.gt.tngmax.and.a3.gt.1.e-20) fftngTQ=tngmax ! intercept [1/cm4]

 lambda = ( fftngTQ * const_pi * roqg / (q_graupel*1.e-6) ) ** 0.25e0  ! lambda [1/cm]

 re = 3.e0 / 2.e0 / lambda * 1.e+4   !effective radius [micron]

 return
 end subroutine re_bulk_god10_graupel

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine get_solar_parameter
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!   Compute parameters related to sun.  
!  
! History:
! 01/2011  Toshi Matsui@NASA GSFC : add azimuth angle (and clean ups)
! 03/2009  Toshi Matsui@NASA GSFC ; Initial.
!           
! References: 
!    Paltridge, G. W., Platt, C. M. R., 1976:Radiative Process in Meteorology and Climatology. New York: Elsevier.
!----------------------------------------------------------------------------------------------------
 integer :: i,j
 real(sdsu_fps) :: julian  !julian date
 real(sdsu_fps) :: declin   !solar delination angle []
 real(sdsu_fps) :: obecl,sinob,sxlong,arg,  &
                   djul,rjul,eccfac, lat_rad, alpha
 real(sdsu_fps) :: tloctm, hrang, xxlat
 real(sdsu_fps),parameter :: dpd = 360.e0/365.e0       ! degrees per day for earth's orbital position (deg/day)

 real(sdsu_fps) :: sol_zen ! solar zenith angle [rad]
 real(sdsu_fps) :: hour_loc     ! local [hr]


 declin=0.e0 ; sdsu_solcon =0.e0 ; eccfac=0.e0

 if(masterproc) print("(a)"),'MSG get_solar_parameter: compute cosine of solar zenith angle and solar constant'
 if(masterproc) print*,' '

!
! from global parameter
!
 julian = (sdsu_julian-1.e0) + sdsu_gmt/24.e0 !julian date

!
! Earth's obliquity = 23.5 degree.
!
 obecl=23.5e0*const_degrad  !earth's oblique angle [rad]
 sinob=sin(obecl)           !sine of oblique angle [-]

!
! Calculate solar declination angle
!
 if(julian >= 80.e0) sxlong = dpd * (julian - 80.e0)  ![rad]
 if(julian <  80.e0) sxlong = dpd * (julian + 285.e0) ![rad]
 arg    = sinob  * sin(sxlong * const_degrad )
 declin = asin(arg)    !solar delination angle [rad]

!
! Solar constant eccentricity factor (paltridge and platt 1976)
!
 djul = julian * dpd
 rjul = djul   * const_degrad
 eccfac = 1.000110e0 + 0.034221e0 * cos(rjul) + 0.001280e0 * sin(rjul) + 0.000719e0 * &
          cos( 2.e0 * rjul ) + 0.000077e0 * sin( 2.e0 * rjul )


!
! solar constant (to be used in broadband simulator)
!
 sdsu_solcon = 1370.e0 * eccfac  ! solar constant [W/m2]



 do j = myj_start, myj_end ; do i = myi_start, myi_end

    !
    ! Cosine of solar zenith angle (used for broadband and visir simulator)
    !
    tloctm = sdsu_gmt + surface(i,j)%lon / 15.e0
    hrang  = 15.e0 * (tloctm - 12.e0) * const_degrad
    xxlat  = surface(i,j)%lat * const_degrad
    surface(i,j)%cosz = sin(xxlat) * sin(declin) + &
                        cos(xxlat) * cos(declin) * cos(hrang)

    !
    ! solar azimuth angle (used for visir simulator)
    !
    sol_zen = acos( surface(i,j)%cosz ) ! solar zenith angle [rad]

    if( visir ) then

     lat_rad = surface(i,j)%lat*const_degrad

     alpha = acos(  (sin(declin) - sin(lat_rad)*cos(sol_zen)) / (cos(lat_rad)*sin(sol_zen)) )

     hour_loc = sdsu_gmt + surface(i,j)%lon/15.
     if( hour_loc < 0. ) hour_loc = hour_loc + 24.
     if( hour_loc > 24.) hour_loc = hour_loc -24. 


     if( hour_loc < 12. ) then ! local AM
       surface(i,j)%solar_azimuth = alpha / const_degrad       ![deg]
     else                      ! local PM
       surface(i,j)%solar_azimuth = ( 2.*const_pi - alpha ) / const_degrad  ![deg]
     endif

    else
      surface(i,j)%solar_azimuth = undefined
    endif

 enddo ; enddo


 return
 end subroutine get_solar_parameter

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU  

 subroutine Find_NaN_Inf_Single(Warning_MSG, real_input, i_in,j_in,k_in)
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  Find out infinity of NaN values for Debugging  
!  
! History:
! 10/2008  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(inout) :: real_input  !anykind of Non-dimensional single-precision Real parameters
 integer,intent(in) :: i_in, j_in, k_in
 character*(*),intent(in) :: Warning_MSG

!
! Find Infinity
!
 if( (real_input > const_largest_real) .or. &
     (real_input < -1*const_largest_real)     ) then
    print*,'MSG Find_NaN_Inf: '//Warning_MSG//'Infinity at',i_in,j_in,k_in
    real_input = undefined
!    call stop_sdsu(myrank,'MSG Find_NaN_Inf: terminate program')
    return
 endif

!
! Find NaN
! 
 if( isnan(real_input) ) then 
    print*,'MSG Find_NaN_Inf: '//Warning_MSG//'NaN at',i_in,j_in,k_in
    real_input = undefined
!    call stop_sdsu(myrank,'MSG Find_NaN_Inf: terminate program') 
    return
 endif

 return
 end subroutine Find_NaN_Inf_Single

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine Find_NaN_Inf_Double(Warning_MSG, dble_input, i_in,j_in,k_in)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  Find out infinity of NaN values for Debugging  
!  
! History:
! 10/2008  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 real(sdsu_fpd),intent(inout) :: dble_input  !anykind of Non-dimensional input Real parameters
 integer,intent(in) :: i_in, j_in, k_in
 character*(*),intent(in) :: Warning_MSG

!
! Find Infinity
!

 if( (dble_input >    const_largest_dble) .or. &
     (dble_input < -1*const_largest_dble)     ) then
    print*,'MSG Find_NaN_Inf: '//Warning_MSG//'Infinity at',i_in,j_in,k_in
    dble_input = undefined_dble
   return
 endif

!
! Find NaN
! 
 if( isnan(dble_input) ) then
    print*,'MSG Find_NaN_Inf: '//Warning_MSG//'NaN at',i_in,j_in,k_in
    dble_input = undefined_dble
    return
 endif

 return
 end subroutine Find_NaN_Inf_Double

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine re_rams_gamma(proc)
 implicit none
!--------------------------------------------------------------------------------------
! Comments:
! This f90 subroutine computes the effective radius for each species of condensate
! within the RAMS microphysics. It operates on a column basis for ease in future 
! parallelism. The basis for the computation is the definition of effective
! radius for a modified gamma distribution (from Miller et al. 2001, JGR):
!
!     Reff = (1/2) * [gamma(nu+3)/gamma(nu+2)] * Dn
!
! Here, "gamma" is the gamma function, "nu" is the width parameter in the
! gamma distribution, and "Dn" is the characteristic diameter.
!
! The gamma function and log gamma function are provided as a subroutine 
! and a function (respectively) in this code and are called as "gamma_reff"
! and "gammln_reff" to avoid conflict with other gamma functions that may be 
! compiled in the driver.
!
! Dn is computed as it is in the RAMS microphysics according to:
! 
!     dn = dnfac(ihcat) * ( emb(k,lcat) ** pwmasi(ihcat) ) * 1.e6
!
! "dnfac" is a precomputed constant that is a function of the mass-diameter
! relationship for each microphsical species and habit. This is compute as
! 
!     dnfac(lhcat) = (cfmasi * exp(glg - glgm)) ** pwmasi(lhcat)
!
! where "lhcat" is an index from 1 to 16 through the possible species and
! habits. Here, g1g and g1gm are
! 
!     glg = gammln(gnu(lcat))
!     glgm = gammln(gnu(lcat) + pwmas(lhcat))
!
! where gammln computes the log of the gamma function
!
! ----------------------------------------------------------------------
!                            COMPUTING EMB
! ----------------------------------------------------------------------
!
! emb needs some explanation--the computation of this depends on whether
! 1 or 2 moment microphysics is being used. Specifically:
!
! In mic_init.F emb is first computed as:
! 
! do lcat = 1,8
!    lhcat = lcat0(lcat)
!    if (jnmb(lcat) .eq. 2) then
!       do k = 2,m1-1
!          emb(k,lcat) = cfmas(lhcat) * parm(lcat) ** pwmas(lhcat)
!       enddo
!    endif
!    do k = 2,m1-1 
!       jhcat(k,lcat) = lhcat
!    enddo
! enddo
! 
! BUT note that this only applies to the situation in which jnmb = 2. This is true for 1-moment 
! microphysics for species 2, 4, 5, 6, 7 (rain, snow, aggregates, graupel, and hail)
! species:          cloud     rain   pristine   snow  agg    grau   hail   cloud2
! jnmb:                 4       2       5        2      2      2      2       4
! irams_flags:          1       1       5        1      1      1      1       1   
! rams_parms:    300000000.0  0.001  100000.0  0.001  0.001  0.001  0.003  100000.0
! 
! Now, if jnmb = 4 (cloud1 and cloud2) then
! 
!    parmi = 1. / parm(lcat)
!       emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat) * parmi))
! 
! where emb0 and emb1 are computed as:
! 
! do lcat=1,ncat
!    lhcat = lcat0(lcat)
!    emb0 (lcat) = cfmas(lhcat) * dstprms(6,lhcat) ** pwmas(lhcat)
!    emb1 (lcat) = cfmas(lhcat) * dstprms(7,lhcat) ** pwmas(lhcat)
! enddo
! 
! Here, dstprms(6,lhcat) is "dmb0" and dstprms(7,lhcat) is "dmb1"--set constant.
! 
! In the case of 2-moment microphysics, emb is computed as:
! 
!       emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat)  &
!          / max(1.e-12,cx(k,lcat))))
! 
! Note that in this case, we need the number concentration (cx is Nc of each species)
!
! ----------------------------------------------------------------------
!
! Procedure:
!
! 1. Read in rams.config file to set parameters and options
! 2. Precompute dnfac (prior to loop over all grid points)
! 3. Loop over all i, j, k
! 4. Look up the habit for pristine ice and snow based on T and RH
! 5. Compute emb
! 6. Compute Dn
! 7. Compute effective radius
! 
! History:
! 04/2009  Toshi Matsui@NASA GSFC ; RAMS routine was Adapted for SDSU (contribution from D. Posselt)
!
! References:
! Walko et al., 1995. R. Walko, W.R. Cotton, M.P. Meyers and J.Y. Harrington, New RAMS cloud 
!    microphysics parameterization. Part I: The single-moment scheme. Atmos. Res. 38 (1995), 
!    pp. 29-62.
!
! Meyers M. P., R. L. Walko, J. Y. Harrington, and W. R. Cotton, 1997: New RAMS cloud microphysics 
!    parameterization. Part II: The two-moment scheme. Atmos. Res., 45, 3-39
!
!-----------------------------------------------------------------------
  character*(*),intent(in) :: proc

  ! Local variables
  integer, parameter :: unit_config=555
  integer, parameter :: nhcat = 16  ! # of fundamental habitat classes
  integer, parameter :: ncat = 8   ! # of simulated condensates classes
  integer :: lhcat !loop index for habitat classes
  integer :: lcat  !loop index for simulated condensate class
  integer :: ihcat
  integer, dimension(ncat) :: jnmb
  real(sdsu_fps),    dimension(ncat) :: emb
  integer, dimension(ncat) :: jhcat
  real(sdsu_fps), dimension(ncat)  :: emb0,emb1,gnu,parm, gamfac
  real(sdsu_fps), dimension(nhcat) :: cfmas,pwmas,dnfac,pwmasi


  integer :: icloud, irain, ipris, isnow, iaggr, igraup, ihail, icloud2
  real(sdsu_fps),dimension(ncat) :: rx, cx
  real(sdsu_fps) :: parmi, cfmasi
  real(sdsu_fps) :: gfac1, gfac2, glg, glgm
  real(sdsu_fps) :: dn_native, dn !characteristic diameter [micron]

  integer, parameter :: lcat0(ncat) = (/1,2,3,4,5,6,7,16/) ! lcat corressponding to lhcat

  integer, parameter :: lhcat0(nhcat) = &  !lhcat corressponding to lcat
                        (/1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/) 

  integer, dimension(1:31,79:100,1:2) :: jhabtab  ! ice mode1 & mode2 habitat index 
                                           ! as a function of T(1~31) and RH(1~100) 
                                           ! 1st dimension: Temperature (1~31)
                                           ! 2nd dimension: RH (1~100)
                                           ! 3rd dimension: 1- ice mode1 ; 2- ice mode2

! Data statement--contains necessary mass parameters
  real(sdsu_fps), dimension(7,nhcat) :: dstprms
  data dstprms/ &
!  --------------------------------------------------------------------------------------
!   shape     cfmas   pwmas     cfvt    pwvt     dmb0      dmb1     lhcat - habit name
!  --------------------------------------------------------------------------------------
!      .5,      524.,     3.,    3173.,     2.,   2.e-6,   40.e-6,  & !1 -cloud
!      .5,      524.,     3.,     144.,   .497,   .1e-3,    5.e-3,  & !2 -rain
!    .179,     110.8,   2.91,    1538.,   1.00,  15.e-6,  125.e-6,  & !3 -pris col
!    .179,  2.739e-3,   1.74,     27.7,   .484,   .1e-3,   10.e-3,  & !4 -snow col
!      .5,      .496,    2.4,     16.1,   .416,   .1e-3,   10.e-3,  & !5 -aggreg
!      .5,      157.,     3.,     332.,   .786,   .1e-3,    5.e-3,  & !6 -graup
!      .5,      471.,     3.,    152.1,   .497,   .8e-3,   10.e-3,  & !7 -hail 
!    .429,     .8854,    2.5,   20801.,  1.377,      00,       00,  & !8 -pris hex
!   .3183,   .377e-2,     2.,     56.4,   .695,      00,       00,  & !9 -pris den
!   .1803,   1.23e-3,    1.8,   1617.9,   .983,      00,       00,  & !10-pris ndl
!      .5,     .1001,  2.256,    6239.,   1.24,      00,       00,  & !11-pris ros
!    .429,     .8854,    2.5,    30.08,   .563,      00,       00,  & !12-snow hex
!   .3183,   .377e-2,     2.,     3.39,   .302,      00,       00,  & !13-snow den
!   .1803,   1.23e-3,    1.8,     44.6,   .522,      00,       00,  & !14-snow ndl
!      .5,     .1001,  2.256,    125.7,   .716,      00,       00,  & !15-snow ros
!      .5,      524.,     3.,   1.26e7,   1.91,  65.e-6,  100.e-6/    !16-cloud2
! --------------------------------------------------------------------------------------
! New table (correct Vt coefficient)
!  --------------------------------------------------------------------------------------
      .5,      524.,     3.,      0.,   .497,   2.e-6,   40.e-6,  & !1 -cloud
      .5,      524.,     3.,     144.,   .497,   .1e-3,    5.e-3,  & !2 -rain
    .179,     110.8,   2.91,     0.0,   .484,  15.e-6,  125.e-6,  & !3 -pris col
    .179,  2.739e-3,   1.74,     27.7,   .484,   .1e-3,   10.e-3,  & !4 -snow col
      .5,      .496,    2.4,     16.1,   .416,   .1e-3,   10.e-3,  & !5 -aggreg
      .5,      157.,     3.,     332.,   .786,   .1e-3,    5.e-3,  & !6 -graup
      .5,      471.,     3.,    152.1,   .497,   .8e-3,   10.e-3,  & !7 -hail 
    .429,     .8854,    2.5,     0.,   .484,      00,       00,  & !8 -pris hex
   .3183,   .377e-2,     2.,     0.,   .484,      00,       00,  & !9 -pris den
   .1803,   1.23e-3,    1.8,     0.,   .484,      00,       00,  & !10-pris ndl
      .5,     .1001,  2.256,     0.,   .484,     00,       00,  & !11-pris ros
    .429,     .8854,    2.5,     27.7,   .484,      00,       00,  & !12-snow hex
   .3183,   .377e-2,     2.,     27.7,   .484,      00,       00,  & !13-snow den
   .1803,   1.23e-3,    1.8,     27.7,   .484,      00,       00,  & !14-snow ndl
      .5,     .1001,  2.256,     27.7,   .484,      00,       00,  & !15-snow ros
      .5,      524.,     3.,      0.,   .497,  65.e-6,  100.e-6/    !16-cloud2



! Set min and max characteristic diameter and effective radius.
! Thest are approximate, and it is up to the user to decide the values and
! to decide whether the range is applied to effective radius or to characteristic diameter...

! Min/Max for characteristic diameter...
  real(sdsu_fps) dnmin(8),dnmax(8)
!              cld,  rain, pris,   snow,    agg,   grau,   hail, cld2
  data dnmin /  1.,    80.,   1.,   125.,    10.,    10.,    10., 40./
  data dnmax / 40., 10000., 125., 10000., 10000., 10000., 10000., 80./

  ! Min/Max for effective radii...
  real(sdsu_fps) remin(8),remax(8)
!              cld,  rain, pris,   snow,    agg,   grau,   hail, cld2
  data remin /  5.,    80.,   5.,   125.,    10.,    10.,    10., 40./
  data remax / 40., 10000., 125., 10000., 10000., 10000., 10000., 80./

  real(sdsu_fps) :: rho_air !air density [g/m3]

  integer :: i, j, k, it, is, nt, ns ! Loop index
  real(sdsu_fps) :: gammln_reff      ! function 
  real(sdsu_fps) :: mean_mass        ! [kg]
  real(sdsu_fps) :: ntot             ! total particle number coonc 

  save

! -------------------------- PROGRAMS START ------------------------------------

 proc_select: select case(proc)
 case('init')

!
! Set up habit tables as a function of temperature and humidity
! (will be used later to diagnose variable crystal habit
! nt is temp, ns = satur (liq)
!
! we can reduce this is  79 ~ 100
  do is = 79,100 ; do it = 1,31
          if (it .ge. 0 .and. it .le. 2) then
              if (is .le. 95) then
                jhabtab(it,is,1) = 3 !pris col
                jhabtab(it,is,2) = 4 !snow col
              else
                jhabtab(it,is,1) = 8  !pris hex
                jhabtab(it,is,2) = 12 !snow hex
              endif
          else if(it .gt. 2 .and. it .le. 4) then
              if (is .lt. 90) then
                jhabtab(it,is,1) = 3  !pris col
                jhabtab(it,is,2) = 4  !snow col
              else
                jhabtab(it,is,1) = 8  !pris hex
                jhabtab(it,is,2) = 12 !snow hex
              endif
          else if(it .gt. 4 .and. it .le. 6) then
              if (is .lt. 85) then
                jhabtab(it,is,1) = 3  !pris col
                jhabtab(it,is,2) = 4  !snow col
              else
                jhabtab(it,is,1) = 10 !pris ndl
                jhabtab(it,is,2) = 14 !snow ndl
              endif
          else if(it .gt. 6 .and. it .le. 9) then
              if (is .lt. 90) then
                jhabtab(it,is,1) = 3  !pris col
                jhabtab(it,is,2) = 4  !snow col
              else
                jhabtab(it,is,1) = 10 !pris ndl
                jhabtab(it,is,2) = 14 !snow ndl
              endif
          else if(it .gt. 9 .and. it .le. 22) then
              if (is .lt. 90) then
                jhabtab(it,is,1) = 8  !pris hex
                jhabtab(it,is,2) = 12 !snow hex
              else
                jhabtab(it,is,1) = 9  !pris den
                jhabtab(it,is,2) = 13 !snow den
              endif
          elseif(it .gt. 22 .and. it .le. 30) then
              if (is .lt. 80) then
                jhabtab(it,is,1) = 3  !pris col
                jhabtab(it,is,2) = 4  !snow col
              else
                jhabtab(it,is,1) = 10 !pris ndl
                jhabtab(it,is,2) = 14 !snow ndl
              endif
          elseif(it .gt. 30) then
              if (is .lt. 90) then
                jhabtab(it,is,1) = 3  !pris col
                jhabtab(it,is,2) = 4  !snow hex
              else
                jhabtab(it,is,1) = 11 !pris ros
                jhabtab(it,is,2) = 15 !snow ros
              endif
          endif

  enddo; enddo

!
! Choose 1moment or 2moment
!
 if( trim(cloud_microphysics) == 'RAMS2') then
    sdsu_io_file = trim(sdsu_dir_data)//'rams.config.2moment'
 else
    sdsu_io_file = trim(sdsu_dir_data)//'rams.config.1moment'
 endif


!
! Read RAMS microphysics configuration file
!
 if(masterproc) print("(a)"), 'MSG re_rams_gamma: Reading RAMS configuration file =>',trim(sdsu_io_file)
 if(masterproc) print*,''
  open ( unit_config, file=trim(sdsu_io_file), form='formatted', action='read', status='old' )

  read ( unit_config, * ) icloud
  read ( unit_config, * ) irain
  read ( unit_config, * ) ipris
  read ( unit_config, * ) isnow
  read ( unit_config, * ) iaggr
  read ( unit_config, * ) igraup
  read ( unit_config, * ) ihail
  read ( unit_config, * ) icloud2

  jnmb(1) = icloud
  jnmb(2) = irain
  jnmb(3) = ipris
  jnmb(4) = isnow
  jnmb(5) = iaggr
  jnmb(6) = igraup
  jnmb(7) = ihail
  jnmb(8) = icloud2   

  ! Fill jnmb according to whether we are using 1 or 2 moment microphysics
! 1moment case
  if (icloud == 1) jnmb(1) = 4
  if (irain  == 1) jnmb(2) = 2
  if (ipris  >= 1) jnmb(3) = 5
  if (isnow  == 1) jnmb(4) = 2
  if (iaggr  == 1) jnmb(5) = 2
  if (igraup == 1) jnmb(6) = 2
  if (ihail  == 1) jnmb(7) = 2
  if (icloud2== 1) jnmb(8) = 4

! 2moment case
  if (irain == 5 .or. isnow == 5 .or. iaggr == 5 .or.  &
      igraup == 5 .or. ihail == 5) then
    if (irain  >= 1) jnmb(2) = 5
    if (isnow  >= 1) jnmb(4) = 5
    if (iaggr  >= 1) jnmb(5) = 5
    if (igraup >= 1) jnmb(6) = 5
    if (ihail  >= 1) jnmb(7) = 5
  endif

  ! Skip past flags that determine whether we are using new power laws and bin riming
  read ( unit_config, * )
  read ( unit_config, * )
  read ( unit_config, * )

  ! Read parameters (characteristic diameter [m])
  read ( unit_config, * ) parm(1)
  read ( unit_config, * ) parm(2)
  read ( unit_config, * ) parm(3)
  read ( unit_config, * ) parm(4)
  read ( unit_config, * ) parm(5)
  read ( unit_config, * ) parm(6)
  read ( unit_config, * ) parm(7)
  read ( unit_config, * ) parm(8)

!1moment options
  if (icloud <= 1) parm(1) = .3e9  !# per kg ~ m^3
  if (irain  == 1) parm(2) = .1e-2
  if (isnow  == 1) parm(4) = .1e-2
  if (iaggr  == 1) parm(5) = .1e-2
  if (igraup == 1) parm(6) = .1e-2
  if (ihail  == 1) parm(7) = .3e-2
  if (icloud2== 1) parm(8) = .1e6  !# per kg ~ m^3 

  ! Read past CCN, GCCN parameters
  read ( unit_config, * )
  read ( unit_config, * )

  ! Read gamma width parameter
  read ( unit_config, * ) gnu(1)
  read ( unit_config, * ) gnu(2)
  read ( unit_config, * ) gnu(3)
  read ( unit_config, * ) gnu(4)
  read ( unit_config, * ) gnu(5)
  read ( unit_config, * ) gnu(6)
  read ( unit_config, * ) gnu(7)
  read ( unit_config, * ) gnu(8)

  close ( unit_config )

!
! Precompute gamma factors as a function of gnu for condensate class
!
  do lcat=1,ncat
    call gamma_reff(gnu(lcat)+3,gfac1)
    call gamma_reff(gnu(lcat)+2,gfac2)
    gamfac(lcat) = gfac1/gfac2
  enddo

!
! Get alpha_m (cfmas) and beta_m (pwmas) for mass-Diameter powerlow relationship for each habit class
!
  do lhcat=1,nhcat
    cfmas(lhcat) = dstprms(2,lhcat)   !alpha m 
    pwmas(lhcat) = dstprms(3,lhcat)   !beta  m
  enddo

!
! Compute maximum and minimum emb for each condensate class
!
  do lcat=1,ncat
    lhcat = lcat0(lcat) !from habit class to simulated condensate class
    emb0 (lcat) = cfmas(lhcat) * dstprms(6,lhcat) ** pwmas(lhcat)  ! minimum emb
    emb1 (lcat) = cfmas(lhcat) * dstprms(7,lhcat) ** pwmas(lhcat)  ! maximum emb
  enddo

!
! Compute constant pwmasi & dnfac for each habit
!
  do lhcat = 1,nhcat   
    lcat=lhcat0(lhcat)  !from simulated condensate class to habit class
    cfmasi = 1. / cfmas(lhcat)
    pwmasi(lhcat) = 1. / pwmas(lhcat)
    glg  = gammln_reff(gnu(lcat))
    glgm = gammln_reff(gnu(lcat) + pwmas(lhcat))
    dnfac(lhcat) = (cfmasi * exp(glg - glgm)) ** pwmasi(lhcat)
  enddo

!
! Passing some parameters to global module_simulater parameter
!
   gnu_rams%cloud1 = gnu(1)   
   gnu_rams%rain   = gnu(2) 
   gnu_rams%ice1   = gnu(3) 
   gnu_rams%ice2   = gnu(4) 
   gnu_rams%snow   = gnu(5) 
   gnu_rams%graupel= gnu(6) 
   gnu_rams%hail   = gnu(7) 
   gnu_rams%cloud2 = gnu(8) 

   rams_dstprms(1:7,1:nhcat) = dstprms(1:7,1:nhcat)
   rams_jhabtab(1:31,79:100,1:2) = jhabtab(1:31,79:100,1:2)


 case('re')

!
!----------------- Start loop over all grid points --------------------
!
  do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end

!
! Get air density for unit converesion
!
     rho_air = 1.e3 * atmos(i,j,k)%press / (atmos(i,j,k)%t_air*2.87)  !dry air density [g/m3]

!
! Compute emb--depends on variables read in from rams.config
!
    LCAT_LOOP: do lcat = 1,8

!
! Particle mixing ratio (rx) [g/g] 
!
      if ( lcat .eq. 1 ) rx(lcat) = q_rams(i,j,k)%cloud1 / rho_air    
      if ( lcat .eq. 2 ) rx(lcat) = q_rams(i,j,k)%rain   / rho_air 
      if ( lcat .eq. 3 ) rx(lcat) = q_rams(i,j,k)%ice1   / rho_air 
      if ( lcat .eq. 4 ) rx(lcat) = q_rams(i,j,k)%ice2   / rho_air 
      if ( lcat .eq. 5 ) rx(lcat) = q_rams(i,j,k)%snow   / rho_air 
      if ( lcat .eq. 6 ) rx(lcat) = q_rams(i,j,k)%graupel/ rho_air 
      if ( lcat .eq. 7 ) rx(lcat) = q_rams(i,j,k)%hail   / rho_air 
      if ( lcat .eq. 8 ) rx(lcat) = q_rams(i,j,k)%cloud2 / rho_air 

!
! Particle number concentration (cx) [#/kg] 
!
      if ( lcat .eq. 1 ) cx(lcat) = n_rams(i,j,k)%cloud1 / rho_air * 1e+3 
      if ( lcat .eq. 2 ) cx(lcat) = n_rams(i,j,k)%rain   / rho_air * 1e+3  
      if ( lcat .eq. 3 ) cx(lcat) = n_rams(i,j,k)%ice1   / rho_air * 1e+3 
      if ( lcat .eq. 4 ) cx(lcat) = n_rams(i,j,k)%ice2   / rho_air * 1e+3 
      if ( lcat .eq. 5 ) cx(lcat) = n_rams(i,j,k)%snow   / rho_air * 1e+3 
      if ( lcat .eq. 6 ) cx(lcat) = n_rams(i,j,k)%graupel/ rho_air * 1e+3 
      if ( lcat .eq. 7 ) cx(lcat) = n_rams(i,j,k)%hail   / rho_air * 1e+3 
      if ( lcat .eq. 8 ) cx(lcat) = n_rams(i,j,k)%cloud2 / rho_air * 1e+3  

      lhcat = lcat0(lcat)  !from condensate class to habit class

!
! Estimate effective mass (emb) for 1moment and 2moment scheme
!
      if (jnmb(lcat) .eq. 2) then     ! 1moment case (rain,snow,aggregate,graupel,hail)

        emb(lcat) = cfmas(lhcat) * parm(lcat) ** pwmas(lhcat)  ![kg]

      elseif (jnmb(lcat) .eq. 4) then ! 1moment case (cloud1 & cloud2)

        parmi = 1. / parm(lcat)
        emb(lcat) = max(emb0(lcat),min(emb1(lcat),rx(lcat) * parmi)) ![g/g]*/[#/kg] -> [kg]

      elseif (jnmb(lcat) >= 5) then   ! 2moment case

        emb(lcat) = max(emb0(lcat),min(emb1(lcat),rx(lcat) / max(1.e-12,cx(lcat))))

      endif

      jhcat(lcat) = lhcat ! habit class for lcat loop

    enddo LCAT_LOOP

!
! adjusting jhcat() for ice1 (pristine ice) and ice2 (snow) classes as a function of T and RH
!
    nt = max(1,min(31,-nint(  (atmos(i,j,k)%t_air-const_Kel2Cel)  )))  ! Temperature index
    ns = max(79,min(100,nint( atmos(i,j,k)%rh )))          ! RH index
    jhcat(3) = jhabtab(nt,ns,1)          ! Habit of pristine ice
    jhcat(4) = jhabtab(nt,ns,2)          ! Habit of snow


    ! Now, compute effective radius for each species
    do lcat = 1,8

      ! Set the index based on condensate species and habit
      ihcat = jhcat(lcat)  !

      ! Compute the characteristic diameter (microns)
      dn_native = dnfac(ihcat) * ( emb(lcat) ** pwmasi(ihcat) ) * 1.e6 ![micron] 

      ! Constrain characteristic diameter for given max and min values (as is done in RAMS)
       dn = max(  dnmin(lcat) , min( dnmax(lcat), dn_native) )   ![micron]

      ! Compute effective radius (output) [micron]
      if ( lcat == 1 ) re_rams(i,j,k)%cloud1  = 0.5 * gamfac(1) * dn
      if ( lcat == 2 ) re_rams(i,j,k)%rain    = 0.5 * gamfac(2) * dn
      if ( lcat == 3 ) re_rams(i,j,k)%ice1    = 0.5 * gamfac(3) * dn
      if ( lcat == 4 ) re_rams(i,j,k)%ice2    = 0.5 * gamfac(4) * dn
      if ( lcat == 5 ) re_rams(i,j,k)%snow    = 0.5 * gamfac(5) * dn
      if ( lcat == 6 ) re_rams(i,j,k)%graupel = 0.5 * gamfac(6) * dn
      if ( lcat == 7 ) re_rams(i,j,k)%hail    = 0.5 * gamfac(7) * dn
      if ( lcat == 8 ) re_rams(i,j,k)%cloud2  = 0.5 * gamfac(8) * dn

!
! get total number concentration for a given characteristic diameter, mass mixing ratio, PSD shape parameter
! , and gamma factors.
!
      call gamma_reff(gnu(lcat)+pwmas(ihcat),gfac1)
      call gamma_reff(gnu(lcat)             ,gfac2)
      mean_mass = cfmas(ihcat) * ( (dn*1e-6) ** pwmas(ihcat) ) * (gfac1/gfac2)  !mean mass [kg]
      ntot = rx(lcat) * 1e-3 * rho_air / mean_mass     ! re-compute total particle number concentrations [#/m3]

      if( ntot==0. .or. ntot>0. .or. ntot<0. .or. ntot>=0. .or. ntot<=0. ) then  
      else ; ntot = 0. ; endif  !filter NaN -> zero

      MOMENT2: if( trim(cloud_microphysics) == 'RAMS1' ) then 
      !adjust number conoc (output)
      if ( lcat == 1 ) n_rams(i,j,k)%cloud1  = ntot
      if ( lcat == 2 ) n_rams(i,j,k)%rain    = ntot
      if ( lcat == 3 ) n_rams(i,j,k)%ice1    = ntot
      if ( lcat == 4 ) n_rams(i,j,k)%ice2    = ntot
      if ( lcat == 5 ) n_rams(i,j,k)%snow    = ntot
      if ( lcat == 6 ) n_rams(i,j,k)%graupel = ntot
      if ( lcat == 7 ) n_rams(i,j,k)%hail    = ntot
      if ( lcat == 8 ) n_rams(i,j,k)%cloud2  = ntot
      endif MOMENT2

    enddo !lcat

!
! Filter numerically insignificant data 
!
   if (q_rams(i,j,k)%cloud1 <q_min_condensate) then 
       q_rams(i,j,k)%cloud1=0.e0;re_rams(i,j,k)%cloud1=0.e0;n_rams(i,j,k)%cloud1=0.e0
   endif
   if (q_rams(i,j,k)%cloud2 <q_min_condensate) then
     q_rams(i,j,k)%cloud2=0.e0;re_rams(i,j,k)%cloud2=0.e0;n_rams(i,j,k)%cloud2=0.e0
   endif
   if (q_rams(i,j,k)%rain   <q_min_condensate) then
     q_rams(i,j,k)%rain  =0.e0;re_rams(i,j,k)%rain  =0.e0;n_rams(i,j,k)%rain  =0.e0
   endif
   if (q_rams(i,j,k)%ice1   <q_min_condensate) then
     q_rams(i,j,k)%ice1  =0.e0;re_rams(i,j,k)%ice1  =0.e0;n_rams(i,j,k)%ice1  =0.e0
   endif
   if (q_rams(i,j,k)%ice2   <q_min_condensate) then
     q_rams(i,j,k)%ice2  =0.e0;re_rams(i,j,k)%ice2  =0.e0;n_rams(i,j,k)%ice2  =0.e0
   endif
   if (q_rams(i,j,k)%snow   <q_min_condensate) then
     q_rams(i,j,k)%snow  =0.e0;re_rams(i,j,k)%snow  =0.e0;n_rams(i,j,k)%snow  =0.e0
   endif
   if (q_rams(i,j,k)%graupel<q_min_condensate) then
     q_rams(i,j,k)%graupel=0.e0;re_rams(i,j,k)%graupel=0.e0;n_rams(i,j,k)%graupel=0.e0
   endif
   if (q_rams(i,j,k)%hail   <q_min_condensate) then
     q_rams(i,j,k)%hail  =0.e0;re_rams(i,j,k)%hail  =0.e0;n_rams(i,j,k)%hail  =0.e0
   endif

  enddo ; enddo; enddo !mxgridx & mxgridy & mxlyr

!
!----------------- End of 3D loops ----------------------------
!


 case default ; call stop_sdsu(myrank,'MSG re_rams_gamma: There is no such proc')
 end select proc_select

  return
 end subroutine re_rams_gamma

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine idealized_scene
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   Set up idealized scene by modulating environmental parameters.   
!                                       
! History:
! 06/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 if(uniform_surface) then
   if(masterproc) print*,'MSG idealized_scene: set up uniform surface conditions.'
   if(masterproc) print*,''

   !
   ! For lat lon is idealized, if they are not specified. 
   !
   if(surface(myi_start,myj_start)%lat == undefined ) then
      surface(:,:)%lat = idealized_surface%lat       
      surface(:,:)%lon = idealized_surface%lon       
   endif

   surface(:,:)%iland     = idealized_surface%iland  
   surface(:,:)%igbp_typ  = idealized_surface%igbp_typ   
   surface(:,:)%elev      = idealized_surface%elev   
   surface(:,:)%frac_veg  = idealized_surface%frac_veg
   surface(:,:)%albedo    = idealized_surface%albedo    
   surface(:,:)%h2o_soil  = idealized_surface%h2o_soil  
   surface(:,:)%h2o_snow  = idealized_surface%h2o_snow   
   surface(:,:)%dhgt_snow = idealized_surface%dhgt_snow

 endif


 if(clear_sky_scene) then
   if(masterproc) print*,'MSG idealized_scene: zero out condensates.'
   if(masterproc) print*,''
   !
   ! zero out condensates 
   !
   mic_select0: select case(type_microphysics)
   case('GEN')
    q_gen(:,:,:)%cloud    = 0.e0
    q_gen(:,:,:)%rain     = 0.e0
    q_gen(:,:,:)%ice      = 0.e0
    q_gen(:,:,:)%snow     = 0.e0
    q_gen(:,:,:)%graupel  = 0.e0
    q_gen(:,:,:)%hail     = 0.e0
   case('RAMS')
    q_rams(:,:,:)%cloud1  = 0.e0 
    q_rams(:,:,:)%cloud2  = 0.e0 
    q_rams(:,:,:)%rain    = 0.e0 
    q_rams(:,:,:)%ice1    = 0.e0 
    q_rams(:,:,:)%ice2    = 0.e0 
    q_rams(:,:,:)%snow    = 0.e0 
    q_rams(:,:,:)%graupel = 0.e0 
    q_rams(:,:,:)%hail    = 0.e0 
   case('SBM')
    q_sbm(:,:,:)%liq     = 0.e0 
    q_sbm(:,:,:)%ice_col = 0.e0 
    q_sbm(:,:,:)%ice_pla = 0.e0 
    q_sbm(:,:,:)%ice_den = 0.e0 
    q_sbm(:,:,:)%snow    = 0.e0 
    q_sbm(:,:,:)%graupel = 0.e0 
    q_sbm(:,:,:)%hail    = 0.e0 
   case default

   end select mic_select0
 endif


 return
 end subroutine idealized_scene

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_others
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:                             
!   compute other parameters from input CRM data  
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: i,j,k
 real(sdsu_fps) :: esat, ewat
 real(sdsu_fps) :: aero(mxlyr,1:mxspc_gocart)   ! bundled aeros mixing ratio [g/m3]
 real(sdsu_fps) :: s   !super saturation [%]
 real(sdsu_fps) :: ccn_out !CCN conc [#/cm3]
 real(sdsu_fps) :: in_out  !IN conc [#/Litter]
 real(sdsu_fps) :: rain_liq, rain_ice ! [mm/hr]
 real(sdsu_fps) :: water_path(4) !1-cloud liq, 2-rain liq, 3-cloud ice, 4-rain ice path [kg/m2]
 integer :: k_0degC

 real(sdsu_fps) :: qgh_col, wmax_col

 if(masterproc) print*,'MSG get_others: calculate other important parameters'
 if(masterproc) print*,''

!
! get other 3D parameters
!

 do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end

   if( trim(sim_case) == 'MMF_GPROF' .and. j > max_sample_gprof ) cycle

   !
   ! get specific humidity from RH and pressure
   !
   esat = 0.611e0 * exp (const_Lv_Rv* ( 1.e0/const_Kel2Cel - 1.e0/atmos(i,j,k)%t_air ) ) *10.e0  ! [hPa]
   ewat = esat * atmos(i,j,k)%rh * 1.e-2                                    ! water vapor pressure [hPa]
   atmos(i,j,k)%sh = ewat * const_Rd_Rv / (atmos(i,j,k)%press-0.378e0*ewat) !specific humitiy [g/g]

   !
   ! exner function
   !
   atmos(i,j,k)%exner = (atmos(i,j,k)%press / 1.e3)**(const_Rd/1004.e0)

   !
   ! layer thickness
   !
   atmos(i,j,k)%dhgt=atmos_stag(i,j,k)%hgt-atmos_stag(i,j,k-1)%hgt  ! thickness [km]

   !
   ! layer height 
   !
   atmos(i,j,k)%hgt=0.5e0*( atmos_stag(i,j,k)%hgt+atmos_stag(i,j,k-1)%hgt)  !layer height [km]

   !
   ! vertical velocity [m/s]
   !
   atmos(i,j,k)%w=0.5e0*( atmos_stag(i,j,k)%w+atmos_stag(i,j,k-1)%w)  !


!tot_cloud
   select case(type_microphysics)
   case('GEN')
     atmos(i,j,k)%tot_cloud = q_gen(i,j,k)%cloud + q_gen(i,j,k)%ice
   case('RAMS')
     atmos(i,j,k)%tot_cloud = q_rams(i,j,k)%cloud1 + q_rams(i,j,k)%cloud2 + q_rams(i,j,k)%ice1 + q_rams(i,j,k)%ice2
   case('SBM')
     atmos(i,j,k)%tot_cloud = q_sbm(i,j,k)%liq + q_sbm(i,j,k)%ice_col + q_sbm(i,j,k)%ice_pla + q_sbm(i,j,k)%ice_den
   case default
   end select

   !
   ! aerosols
   !
   if(account_aerosol .and. trim(aerosol_microphysics) == 'GOCART') then

    aero(k,1) =q_gocart(i,j,k)%so4  !all unit [g/m3]
    aero(k,2) =q_gocart(i,j,k)%blc
    aero(k,3) =q_gocart(i,j,k)%ocn
    aero(k,4) =q_gocart(i,j,k)%och
    aero(k,5) =q_gocart(i,j,k)%ssa
    aero(k,6) =q_gocart(i,j,k)%ssc
    aero(k,7) =q_gocart(i,j,k)%du1
    aero(k,8) =q_gocart(i,j,k)%du2
    aero(k,9) =q_gocart(i,j,k)%du3
    aero(k,10)=q_gocart(i,j,k)%du4
    aero(k,11)=q_gocart(i,j,k)%du5
    aero(k,12)=q_gocart(i,j,k)%du6
    aero(k,13)=q_gocart(i,j,k)%du7
    aero(k,14)=q_gocart(i,j,k)%du8

     s = max(0.e0, atmos(i,j,k)%rh - 100.e0) !super saturation [%]

     call mass2ccn(atmos(i,j,k)%t_air,s,aero(k,:),ccn_out)
     atmos(i,j,k)%ccn = ccn_out  !CCN conc [#/cm3]

     call mass2icn (atmos(i,j,k)%press,atmos(i,j,k)%t_air,aero(k,:),in_out)
     atmos(i,j,k)%icn = in_out   !IN conc [#/Litter]

     endif


 enddo ; enddo ; enddo

!
! get other 2D parameters
!
 do j = myj_start, myj_end ; do i = myi_start, myi_end

   if( trim(sim_case) == 'MMF_GPROF' .and. j > max_sample_gprof ) cycle
   !
   ! near surface air temperature
   !
   surface(i,j)%t_air = atmos_stag(i,j,0)%t_air  !air temperaute at 10m AGL 

   !
   ! column water vapor [kg/m2]
   !
    surface(i,j)%h2o_col = sum( atmos(i,j,1:mxlyr)%sh * atmos(i,j,1:mxlyr)%rho_dair &
                               *  atmos(i,j,1:mxlyr)%dhgt )


!print*,i,j, surface(i,j)%h2o_col
!if(j == 106 .or. j == 107) then
!    print*, 'sh',sum(atmos(i,j,1:mxlyr)%sh)
!    print*, 'rho_dair',sum(atmos(i,j,1:mxlyr)%rho_dair)
!    print*, 'dhgt',sum(atmos(i,j,1:mxlyr)%dhgt)
!endif


   if(surface(i,j)%h2o_col <= 0.) then
      print*, 'sh=', atmos(i,j,1:mxlyr)%sh
      print*, 'rho=', atmos(i,j,1:mxlyr)%rho_dair
      print*, 'dhgt=', atmos(i,j,1:mxlyr)%dhgt
      print*, 'hgt=', atmos_stag(i,j,0:mxlyr)%hgt
      print*, '(i,j)=',i,j
      stop 'find negative h2o_col. Check input.'
   endif

   !
   ! column-integrated condensate  [kg/m2] and instanteneous surface rate [mm/hr]
   !
   mic_select0: select case(type_microphysics)

   case('GEN')
    qcol_gen(i,j)%cloud   = sum(q_gen(i,j,1:mxlyr)%cloud    * atmos(i,j,1:mxlyr)%dhgt)
    qcol_gen(i,j)%rain    = sum(q_gen(i,j,1:mxlyr)%rain     * atmos(i,j,1:mxlyr)%dhgt)
    qcol_gen(i,j)%ice     = sum(q_gen(i,j,1:mxlyr)%ice      * atmos(i,j,1:mxlyr)%dhgt)
    qcol_gen(i,j)%snow    = sum(q_gen(i,j,1:mxlyr)%snow     * atmos(i,j,1:mxlyr)%dhgt)
    qcol_gen(i,j)%graupel = sum(q_gen(i,j,1:mxlyr)%graupel  * atmos(i,j,1:mxlyr)%dhgt)
    qcol_gen(i,j)%hail    = sum(q_gen(i,j,1:mxlyr)%hail     * atmos(i,j,1:mxlyr)%dhgt)

!toshii (temporal for paper revision)
! temporal [km] (z=0: bottom, z=mxlyr: top)
!k_0degC = 1
!do k = 1, mxlyr
!   if( atmos(i,j,k)%t_air < 273.) then
!!       k_0degC = k ; exit 
!endif
!enddo
!qcol_gen(i,j)%cloud   = sum(q_gen(i,j,1:k_0degC)%cloud    * atmos(i,j,1:k_0degC)%dhgt)  !cloud liquid in warm-temp only.

    surface(i,j)%path_cloud_liq = qcol_gen(i,j)%cloud
    surface(i,j)%path_rain_liq  = qcol_gen(i,j)%rain
    surface(i,j)%path_cloud_ice = qcol_gen(i,j)%ice    !+ qcol_gen(i,j)%snow !toshii (temporal to add snow)
    surface(i,j)%path_rain_ice  = qcol_gen(i,j)%snow + qcol_gen(i,j)%graupel + qcol_gen(i,j)%hail

!
! special routine to remove cirrus deck from MMF
!
    if( trim(sim_case) == 'MMF' .or. trim(sim_case) == 'MMF_GPROF' ) then
       if(mmf_remove_cirrus) then
          call remove_cirrus_deck( qcol_gen(i,j), q_gen(i,j,1:mxlyr) )
          !if modified, then re-calculate path amount
          qcol_gen(i,j)%ice     = sum(q_gen(i,j,1:mxlyr)%ice      * atmos(i,j,1:mxlyr)%dhgt)
          qcol_gen(i,j)%snow    = sum(q_gen(i,j,1:mxlyr)%snow     * atmos(i,j,1:mxlyr)%dhgt)
          qcol_gen(i,j)%graupel = sum(q_gen(i,j,1:mxlyr)%graupel  * atmos(i,j,1:mxlyr)%dhgt)
          qcol_gen(i,j)%hail    = sum(q_gen(i,j,1:mxlyr)%hail     * atmos(i,j,1:mxlyr)%dhgt)
          surface(i,j)%path_cloud_ice = qcol_gen(i,j)%ice
          surface(i,j)%path_rain_ice  = qcol_gen(i,j)%snow + qcol_gen(i,j)%graupel + qcol_gen(i,j)%hail

       endif
    endif 


    call bulk_rain( atmos(i,j,1)%press, q_gen(i,j,1), re_gen(i,j,1), &
                    rain_liq, rain_ice   )
    surface(i,j)%rain_rate_liq = rain_liq !surface liquid rain rate [mm/hr]
    surface(i,j)%rain_rate_ice = rain_ice !surface liquid rain rate [mm/hr]
    surface(i,j)%rain_rate = surface(i,j)%rain_rate_liq + surface(i,j)%rain_rate_ice !total surface rain rate

    do k = 1, mxlyr
       call bulk_rain( atmos(i,j,k)%press, q_gen(i,j,k), re_gen(i,j,k), &
                       rain_liq, rain_ice   )
       atmos(i,j,k)%rain_rate = rain_liq + rain_ice  ! rain rate profile [mm/hr]
    enddo




   case('RAMS')
    qcol_rams(i,j)%cloud1  = sum( q_rams(i,j,1:mxlyr)%cloud1  * atmos(i,j,1:mxlyr)%dhgt )
    qcol_rams(i,j)%cloud2  = sum( q_rams(i,j,1:mxlyr)%cloud2  * atmos(i,j,1:mxlyr)%dhgt )
    qcol_rams(i,j)%rain    = sum( q_rams(i,j,1:mxlyr)%rain    * atmos(i,j,1:mxlyr)%dhgt )
    qcol_rams(i,j)%ice1    = sum( q_rams(i,j,1:mxlyr)%ice1    * atmos(i,j,1:mxlyr)%dhgt )
    qcol_rams(i,j)%ice2    = sum( q_rams(i,j,1:mxlyr)%ice2    * atmos(i,j,1:mxlyr)%dhgt )
    qcol_rams(i,j)%snow    = sum( q_rams(i,j,1:mxlyr)%snow    * atmos(i,j,1:mxlyr)%dhgt )
    qcol_rams(i,j)%graupel = sum( q_rams(i,j,1:mxlyr)%graupel * atmos(i,j,1:mxlyr)%dhgt )
    qcol_rams(i,j)%hail    = sum( q_rams(i,j,1:mxlyr)%hail    * atmos(i,j,1:mxlyr)%dhgt )

    surface(i,j)%path_cloud_liq = qcol_rams(i,j)%cloud1 + qcol_rams(i,j)%cloud2
    surface(i,j)%path_rain_liq  = qcol_rams(i,j)%rain
    surface(i,j)%path_cloud_ice = qcol_rams(i,j)%ice1 
    surface(i,j)%path_rain_ice  = qcol_rams(i,j)%ice2 + qcol_rams(i,j)%snow &
                                + qcol_rams(i,j)%graupel + qcol_rams(i,j)%hail

   ! has not yet coded for RAMS microphysics rain rate

   case('SBM')
    qcol_sbm(i,j)%liq     = sum(q_sbm(i,j,1:mxlyr)%liq     * atmos(i,j,1:mxlyr)%dhgt)
    qcol_sbm(i,j)%ice_col = sum(q_sbm(i,j,1:mxlyr)%ice_col * atmos(i,j,1:mxlyr)%dhgt)
    qcol_sbm(i,j)%ice_pla = sum(q_sbm(i,j,1:mxlyr)%ice_pla * atmos(i,j,1:mxlyr)%dhgt)
    qcol_sbm(i,j)%ice_den = sum(q_sbm(i,j,1:mxlyr)%ice_den * atmos(i,j,1:mxlyr)%dhgt)
    qcol_sbm(i,j)%snow    = sum(q_sbm(i,j,1:mxlyr)%snow    * atmos(i,j,1:mxlyr)%dhgt)
    qcol_sbm(i,j)%graupel = sum(q_sbm(i,j,1:mxlyr)%graupel * atmos(i,j,1:mxlyr)%dhgt)
    qcol_sbm(i,j)%hail    = sum(q_sbm(i,j,1:mxlyr)%hail    * atmos(i,j,1:mxlyr)%dhgt)

    call sbm_path(  atmos(i,j,:)%dhgt, n_sbm(i,j,:,:), water_path(1:4) )
    surface(i,j)%path_cloud_liq = water_path(1)  !cloud liquid [kg/m2]
    surface(i,j)%path_rain_liq  = water_path(2)  !rain liquid  [kg/m2]
    surface(i,j)%path_cloud_ice = water_path(3)  !cloud ice    [kg/m2]
    surface(i,j)%path_rain_ice  = water_path(4)  !rain liquid  [kg/m2]

    call sbm_rain(  atmos(i,j,1)%press, n_sbm(i,j,1,:), frime_sbm_snow(i,j,1,:), &
                   rain_liq, rain_ice   )    
    surface(i,j)%rain_rate_liq = rain_liq !surface liquid rain rate [mm/hr]
    surface(i,j)%rain_rate_ice = rain_ice !surface liquid rain rate [mm/hr]
    surface(i,j)%rain_rate = surface(i,j)%rain_rate_liq + surface(i,j)%rain_rate_ice !total surface rain rate

    do k = 1, mxlyr
       call sbm_rain(  atmos(i,j,k)%press, n_sbm(i,j,k,:), frime_sbm_snow(i,j,k,:), &
                       rain_liq, rain_ice   )
       atmos(i,j,k)%rain_rate = rain_liq + rain_ice  ! rain rate profile [mm/hr]
    enddo

   case default 

   end select mic_select0


   if( trim(sim_case) == 'MMF_GPROF' ) then
       surface(i,j)%hgt_0degC = undefined
       surface(i,j)%wcd       = undefined
       surface(i,j)%cape      = undefined
       surface(i,j)%ncape     = undefined
       surface(i,j)%lcl       = undefined
   else
   !
   ! freezing level AGL [m]
   !
   call freezing_level( mxlyr, atmos_stag(i,j,0:mxlyr)%t_air , atmos_stag(i,j,0:mxlyr)%hgt, &
                         surface(i,j)%hgt_0degC )

   !
   ! derive various thermodynamic parameters [km]
   !
   call compute_cape( undefined, mxlyr, atmos(i,j,:)%qv, atmos(i,j,:)%hgt*1000., atmos(i,j,:)%dhgt*1000., &
                      atmos(i,j,:)%t_air, atmos(i,j,:)%rh/100., atmos(i,j,:)%press, &
                      surface(i,j)%cape, surface(i,j)%ncape, surface(i,j)%lcl )

   surface(i,j)%wcd = max( 0., surface(i,j)%hgt_0degC - surface(i,j)%lcl )  !warm cloud depth

   
   endif


   !
   ! find the column class (iclass)
   !

   select case(type_microphysics)
   case('GEN')
     qgh_col = qcol_gen(i,j)%graupel + qcol_gen(i,j)%hail
   case('RAMS')
     qgh_col = qcol_rams(i,j)%graupel + qcol_rams(i,j)%hail
   case('SBM')
     qgh_col = qcol_sbm(i,j)%graupel + qcol_sbm(i,j)%hail
   case default
   end select

   wmax_col = maxval( atmos(i,j,:)%w )  !maximum vertical velocity 

   call column_category( surface(i,j)%path_rain_liq,qgh_col, wmax_col,  atmos(i,j,:)%tot_cloud, atmos(i,j,:)%hgt, atmos(i,j,:)%press, &
                         surface(i,j)%lcl, surface(i,j)%iclass )

 enddo ; enddo


!
! domain statistics
!
 domain_stat%maxlat = maxval(surface(:,:)%lat)   !domain maximum latitudue [deg]
 domain_stat%maxlon = maxval(surface(:,:)%lon)   !domain maximum longitude [deg]
 domain_stat%minlat = minval(surface(:,:)%lat)   !domain minimum latitudue [deg]
 domain_stat%minlon = minval(surface(:,:)%lon)   !domain minimum longitude [deg]

#if MPI == 2

  call mpi_sdsu_allreduce('MAX', domain_stat%maxlat )
  call mpi_sdsu_allreduce('MAX', domain_stat%maxlon )
  call mpi_sdsu_allreduce('MIN', domain_stat%minlat )
  call mpi_sdsu_allreduce('MIN', domain_stat%minlon )

#endif

! print*,myrank,'domain_stat',domain_stat 
! call stop_sdsu(myrank,'test')

!
! broadcast domain lat lon
!
#if MPI == 2

  do j = myj_start, myj_end ; do i = myi_start, myi_end
     lat_domain(i,j)  = surface(i,j)%lat
     lon_domain(i,j)  = surface(i,j)%lon
     elev_domain(i,j) = surface(i,j)%elev
  enddo ; enddo  

!
! let's send each tile into the all processor. 
!
   call mpi_sdsu_communicate( 'TO_ALL', lat_domain )
   call mpi_sdsu_communicate( 'TO_ALL', lon_domain )
   call mpi_sdsu_communicate( 'TO_ALL', elev_domain )

#else

    lat_domain = surface%lat  !latitude  [deg]
    lon_domain = surface%lon  !longitude [deg] 
    elev_domain = surface%elev

#endif
   
 return
 end subroutine get_others

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine freezing_level( kmax, tair, hgt, hgt_0degC )
 implicit none
 integer,intent(in) :: kmax
 real(sdsu_fps), intent(in) :: tair(0:kmax)  ! air temperature [K]
 real(sdsu_fps), intent(in) :: hgt(0:kmax)   ! height [km]
 real(sdsu_fps), intent(out) :: hgt_0degC ! freezing height  [m]

 integer :: k
 real(sdsu_fps) :: wgt1, wgt2

 if( minval(tair(:) ) >  const_Kel2Cel ) then
      hgt_0degC = undefined
      return
 endif
 if (maxval( tair(:) ) < const_Kel2Cel ) then
       hgt_0degC = hgt(0)*1000.e0
       return
 endif
  
 do k = 0, kmax-1

    if( tair(k) > const_Kel2Cel .and. tair(k+1) <= const_Kel2Cel ) then 
       wgt2 = (tair(k) - const_Kel2Cel)  / (tair(k) - tair(k+1))  ! 
       wgt1 = 1.-wgt2  !0.33
       hgt_0degC = (wgt1 * hgt(k) + wgt2 * hgt(k+1) ) * 1000.e0
       exit
    endif

 enddo

 return
 end subroutine freezing_level

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine convert_usgs_igbp(usgs_typ, igbp_typ)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   Convert USGU land-cover type into IGBP land-cover type 
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer,intent(in) :: usgs_typ
 integer,intent(out):: igbp_typ

!--------------IGBP LULC type--------------
!water body =                    0
!evergreen needleleaf forest =   1
!evergreen broadleaf forest =    2
!deciduous needleleaf forest =   3
!deciduous broadleaf forest =    4
!mixed forests =                 5
!closed shrubland =              6
!open shrublands =               7
!woody savannas =                8 
!savannas =                      9
!grasslands =                   10 
!permanent wetlands =           11  
!croplands =                    12
!urban and built-up =           13
!cropland/natural vegetation mosaic =  14
!snow and ice =                 15
!barren or sparsely vegetated = 16
!-------------------------------------------- 

 usgs_select: select case(usgs_typ)
     case(1) ; igbp_typ = 13    ! (USGS) Urban and Built-up Land
     case(2) ; igbp_typ = 12    ! (USGS) Dryland Cropland and Pasture
     case(3) ; igbp_typ = 12    ! (USGS) Irrigated Cropland and Pasture
     case(4) ; igbp_typ = 12    ! (USGS) Mixed Dryland/Irrigated Cropland and Pasture
     case(5) ; igbp_typ = 14    ! (USGS) Cropland/Grassland Mosaic
     case(6) ; igbp_typ = 14    ! (USGS) Cropland/Woodland Mosaic
     case(7) ; igbp_typ = 10    ! (USGS) Grassland
     case(8) ; igbp_typ = 6    ! (USGS) Shrubland
     case(9) ; igbp_typ = 7    ! (USGS) Mixed Shrubland/Grassland
     case(10); igbp_typ = 9    ! (USGS) Savanna
     case(11); igbp_typ = 4    ! (USGS) Deciduous Broadleaf Forest
     case(12); igbp_typ = 3    ! (USGS) Deciduous Needleleaf Forest
     case(13); igbp_typ = 2    ! (USGS) Evergreen Broadleaf
     case(14); igbp_typ = 1    ! (USGS) Evergreen Needleleaf
     case(15); igbp_typ = 5    ! (USGS) Mixed Forest
     case(16); igbp_typ = 0    ! (USGS) Water Bodies
     case(17); igbp_typ = 11    ! (USGS) Herbaceous Wetland
     case(18); igbp_typ = 11    ! (USGS) Wooden Wetland
     case(19); igbp_typ = 16    ! (USGS) Barren or Sparsely Vegetated
     case(20); igbp_typ = 7    ! (USGS) Herbaceous Tundra
     case(21); igbp_typ = 7    ! (USGS) Wooded Tundra
     case(22); igbp_typ = 7    ! (USGS) Mixed Tundra
     case(23); igbp_typ = 10    ! (USGS) Bare Ground Tundra
     case(24); igbp_typ = 15    ! (USGS) Snow or Ice
 
     case default     
       print*,'MSG convert_usgs_igbp: There is no such usgs_typ',usgs_typ
       call stop_sdsu(myrank,'MSG convert_usgs_igbp: terminate program') 
 end select usgs_select

 end subroutine convert_usgs_igbp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_gocart_psd
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   Read GOCART's particle-size distribution functions from LUT.  
!                                       
! History:
! 01/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: t,k,nh
 integer :: nrs
 character :: aa*39   ! dummy character
 character :: at*4    ! aerosol types
 character*4,dimension(6),parameter :: mie_tbl = &
       (/'suso', 'waso', 'soot',  'ssam', 'sscm', 'dust'/)
!--------------local variables-----------------------------
 integer,parameter :: nmt = 6    !# of mie tables
 integer,parameter :: nrh  = 36  !# of column (RH) for mie table
 integer,parameter :: nrmx = 99
 real(8) :: mRHr(nmt,nrh)      ! relative humidity [-]
 real(8) :: mRMr(nmt,nrh)      ! relative humidity [-]
 real(8) :: nRMr(mxspc_gocart,0:nrmx)     ! reasigned relative humidity [-]
 real(8) :: mREr(nmt,nrh)      ! relative humidity [-]
 real(8) :: nREr(mxspc_gocart,0:nrmx)     ! reasigned relative humidity [-]
 real(8) :: w1, w2             !weight for pressure interpolation

!
! Mie table loop
! 1- "suso" : sulfate and its precursors
! 2- "waso" : water soluble organic carbon
! 3- "soot" : black carbon
! 4- "ssam" : sea salt (accumulation,fine mode)
! 5- "sscm" : sea salt (coarse mode)
! 6- "dust" : dust

  if(masterproc .and. verbose_SDSU) print*,'MSG read_gocart_psd; read mie tables'

  do t = 1, nmt  !# of mie tables
     at=mie_tbl(t)
     nrs = nrh
     if(at == "dust") nrs = 8
     open(1,file=trim(sdsu_dir_data)//'mie.'//at//'.clirad',status='old')
     read(1,*)aa  !read comments
     read(1,'(6x,36f10.2)')(mRHr(t,k),k=1,nrs)
     read(1,'(6x,36f10.2)')(mRMr(t,k),k=1,nrs)
     read(1,'(6x,36f10.2)')(mREr(t,k),k=1,nrs)
     close(1)
  enddo

!
!interpolate for RH 0~79%
!
  do t = 1, nmt-1  !spc loop (excepting dust)
     do nh = 0, 79  !RH 0~79%
       do k = 1, 16
         if(nh == nint(mRHr(t,k)*100.) ) then
          nRMr(t,nh) = mRMr(t,k)
          nREr(t,nh) = mREr(t,k)
          exit
        elseif( nh > nint(mRHr(t,k)*100.) .and. &
                nh < nint(mRHr(t,k+1)*100.)    ) then
          w2 = ( real(nh)/100. - mRHr(t,k) ) / &
               (mRHr(t,k+1)-mRHr(t,k) )
          w1 = max(min( 1.-w2,1. ),0.)
          nRMr(t,nh) = w1*mRMr(t,k)+w2*mRMr(t,k+1)
          nREr(t,nh) = w1*mREr(t,k)+w2*mREr(t,k+1)
          exit
         endif
       enddo
     enddo
  enddo

!
!assign for RH 80~99%
!
  do t = 1, nmt-1  !spc loop (excepting dust)
     do nh = 80, nrmx  !RH 90~99%
        nRMr(t,nh) = mRMr(t,nh-63)
        nREr(t,nh) = mREr(t,nh-63)
     enddo
  enddo

!
! passing to the module global parameters
!

!
! non dust case (for each RH index)
!
  re_gocart(1:nmt-1,0:nrmx) = nREr(1:nmt-1,0:nrmx)
  rm_gocart(1:nmt-1,0:nrmx) = nRMr(1:nmt-1,0:nrmx)

!
! dust case (for 8 size mode)
!
  re_gocart(nmt,1:8) = mREr(nmt,1:8)  !effective radius [um]
  rm_gocart(nmt,1:8) = mRMr(nmt,1:8)  !mode radius [um]

  return
 end subroutine read_gocart_psd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine nudge_gocart_offline
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!  Read GOCART simulation result from netcdf file.   
!                                       
! History:
! 01/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!  GOCART aerosol species indice
!     index                                  visir index
!  1 = SO4        !sulfur and its precure      "so4"
!  2 = BC1+BC2    !black carbon (soot)         "blc"
!  3 = OC1        !non hygroscopic OC          "ocn"
!  4 = OC2        !hygroscopic OC              "och"
!  5 = SS1        !sea-salt accumulation mode  "ssa"
!  6 = SS2+SS3+SS4!sea-salt coarse mode        "ssc"
!  7 = DU1        ! dust mode 1                "du1"
!  8 = DU1        ! dust mode 2                "du2"
!  9 = DU1        ! dust mode 3                "du3"
! 10 = DU1        ! dust mode 4                "du4"
! 11 = DU2        ! dust mode 5                "du5"
! 12 = DU3        ! dust mode 6                "du6"
! 13 = DU4        ! dust mode 7                "du7"
! 14 = DU5        ! dust mode 8                "du8"
!----------------------------------------------------------------------------------------------------
 integer :: i,j,k
 integer :: is,ie,js,je,ks,ke !memory domain parameter
 integer :: julday    ! julian day 
 integer :: julyr     ! year
 integer,save :: julday_prev  !previous julidan day 
 real(sdsu_fps) :: julian       ! julian day in real  
 real(sdsu_fps) :: aero(mxlyr,1:mxspc_gocart)   ! bundled aeros mixing ratio [g/m3]
! real :: lat(1,1), lon(1,1)  !test for GCE


!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=1 ; ke=mxlyr

!
! Day/time parameter needed for gocart_read subroutine from global module parameter
!
 julday = int(sdsu_julian)
 julian = (sdsu_julian-1.e0) + sdsu_gmt/24.e0
 julyr  = int(sdsu_yyyy)

!
! GOCART IO PROCESS (this must be called whenever julidan date is changed. )
! Everyday, it open/read/store new GOCART global aerosol files
 if( julday /= julday_prev ) then
    call gocart_read( is,ie, js,je, surface(is:ie,js:je)%lat, surface(is:ie,js:je)%lon,&
                      julian, julyr )  !new routine to read subset of gocart
    julday_prev = julday
 endif


!    
! Spatially Temporally Interpolate Aerosol Mass Conc to get 1D aero(1:mxlyr,1:mxspc_gocart) [g/m3] 
! aerosol type index (1 ~ 14)
!
 do j = js, je ; do i = is, ie 

    call aero_interp(mxlyr, atmos(i,j,1:mxlyr)%press, atmos(i,j,1:mxlyr)%t_air, aero(1:mxlyr,1:mxspc_gocart), &
                     surface(i,j)%lat, surface(i,j)%lon, sdsu_gmt)

    q_gocart(i,j,1:mxlyr)%so4=aero(1:mxlyr,1)  !all unit [g/m3]
    q_gocart(i,j,1:mxlyr)%blc=aero(1:mxlyr,2)
    q_gocart(i,j,1:mxlyr)%ocn=aero(1:mxlyr,3)
    q_gocart(i,j,1:mxlyr)%och=aero(1:mxlyr,4)
    q_gocart(i,j,1:mxlyr)%ssa=aero(1:mxlyr,5)
    q_gocart(i,j,1:mxlyr)%ssc=aero(1:mxlyr,6)
    q_gocart(i,j,1:mxlyr)%du1=aero(1:mxlyr,7)
    q_gocart(i,j,1:mxlyr)%du2=aero(1:mxlyr,8)
    q_gocart(i,j,1:mxlyr)%du3=aero(1:mxlyr,9)
    q_gocart(i,j,1:mxlyr)%du4=aero(1:mxlyr,10)
    q_gocart(i,j,1:mxlyr)%du5=aero(1:mxlyr,11)
    q_gocart(i,j,1:mxlyr)%du6=aero(1:mxlyr,12)
    q_gocart(i,j,1:mxlyr)%du7=aero(1:mxlyr,13)
    q_gocart(i,j,1:mxlyr)%du8=aero(1:mxlyr,14)

 enddo ; enddo 

 return
 end subroutine nudge_gocart_offline

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine refine_range_bin_bundle( mxlyr_crm, hgt_stag_crm, dhgt_crm, mxlyr_inst, hgt_stag_inst, &
                                     np, var1d_crm, var1d_inst ) 
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!   This routine modify variables from 
!   CRM level to instrumental level.    (for single precision parameter) 
!   Interface of instrument range bin start from earth's ellipsoid.
!
! History:
!   05/2011  Toshi Matsui@NASA GSFC :: Initial
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)        :: mxlyr_crm                  ! maximum vertical layer of CRM
 real(sdsu_fps),intent(in) :: hgt_stag_crm(0:mxlyr_crm)  ! height leve at CRM interface [km] 
 real(sdsu_fps),intent(in) :: dhgt_crm(1:mxlyr_crm)      ! layer depth of CRM vertical height [km] 
 integer,intent(in)        :: mxlyr_inst                 ! maximum layer for radar instrument. 
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)! instrument measurement height at interface [km]
 integer,intent(in)        :: np
 real(sdsu_fps),intent(in) :: var1d_crm(1:np,1:mxlyr_crm)     ! parameter in CRM level []
 real(sdsu_fps),intent(out) :: var1d_inst(1:np,1:mxlyr_inst)  ! parameter in instrument level  []

 integer :: k_crm, k_inst !vertical loop indice
 integer :: k_crm_low, k_crm_hig ! 
 real(sdsu_fps) :: hgt_btm_inst, hgt_top_inst !layer interface height [km]
 real(sdsu_fps) :: range_bin ! instrument range bin height [km]
 real(sdsu_fps) :: dhgt_low, dhgt_hig
 real(sdsu_fps) :: dhgt_surf  !surface residual layer depth [km]

 integer :: n

!
! Looping of instrumental range
!
 INSTR_RANGE_LOOP: do k_inst = 1, mxlyr_inst

   hgt_btm_inst = hgt_stag_inst(k_inst-1)  ![km]
   hgt_top_inst = hgt_stag_inst(k_inst)    ![km]

   !
   ! range bin completely under the ground surface  --> all undefine
   !
! old code
!   if(hgt_btm_inst < hgt_stag_crm(0) .or. hgt_top_inst < hgt_stag_crm(0) ) then
!      var1d_inst (k_inst) = undefined
!      cycle INSTR_RANGE_LOOP
!   endif

! new code 
   !
   ! range bin completely under the ground surface  --> all undefine
   !
   if(hgt_btm_inst < hgt_stag_crm(0) .and. hgt_top_inst < hgt_stag_crm(0) ) then
      var1d_inst (1:np,k_inst) = undefined
      cycle INSTR_RANGE_LOOP
   endif

   !
   ! range bin completely above the CRM-defined top layer  --> all undefine
   !
   if(hgt_btm_inst >= hgt_stag_crm(mxlyr_crm) .or. hgt_top_inst >= hgt_stag_crm(mxlyr_crm) ) then
      var1d_inst (1:np,k_inst) = undefined
      cycle INSTR_RANGE_LOOP
   endif

   !
   ! range bin partially under the ground  --> surface index
   !
   if(hgt_btm_inst <= hgt_stag_crm(0) .and. hgt_top_inst >= hgt_stag_crm(0) ) then
      dhgt_surf = hgt_top_inst - hgt_stag_crm(0)  !surface residual layer depth [km]

      if( dhgt_surf/(hgt_top_inst-hgt_btm_inst) < 0.5e0 ) then !if residual layer < 50% of instr layer
            var1d_inst (1:np,k_inst) = undefined 
            cycle INSTR_RANGE_LOOP
      endif

      k_crm_low = 1
      dhgt_low  = hgt_stag_crm(1) - hgt_stag_crm(0)

   else
   !
   ! otherwise interpolate/extrapolate CRM-level value into instrumental range bin
   !

   ! get k_crm_low index
     do k_crm = 1, mxlyr_crm
        if( hgt_btm_inst >= hgt_stag_crm(k_crm-1) .and. hgt_btm_inst < hgt_stag_crm(k_crm) ) then
            k_crm_low = k_crm
            dhgt_low = hgt_stag_crm(k_crm) - hgt_btm_inst  ! [km] 
            exit
        endif
     enddo
   endif

   ! get k_crm_hig index
   do k_crm = 1, mxlyr_crm
      if( hgt_top_inst >= hgt_stag_crm(k_crm-1) .and. hgt_top_inst < hgt_stag_crm(k_crm) ) then
          k_crm_hig = k_crm
          dhgt_hig = hgt_top_inst - hgt_stag_crm(k_crm-1)    ! [km] 
          exit
      endif
   enddo

   !interplate
   if( k_crm_low == k_crm_hig ) then          !within same CRM layer

       var1d_inst(1:np,k_inst) = var1d_crm(1:np,k_crm_low)

   elseif( k_crm_low + 1 == k_crm_hig ) then  !just 1 vertical index lag

       range_bin = hgt_stag_inst(k_inst) - hgt_stag_inst(k_inst-1) ![km]

       var1d_inst(1:np,k_inst) = (dhgt_low * var1d_crm(1:np,k_crm_low) + dhgt_hig * var1d_crm(1:np,k_crm_hig) ) &
                               / (dhgt_low+dhgt_hig )   !dhgt_low+dhgt_hig= range_bin

   else  ! More than 1 vertical index lag. 

       range_bin = hgt_stag_inst(k_inst) - hgt_stag_inst(k_inst-1) ![km]

       do n = 1, np
          var1d_inst(n,k_inst) = ( dhgt_low * var1d_crm(n,k_crm_low) + dhgt_hig * var1d_crm(n,k_crm_hig)  &
                                +SUM(dhgt_crm(k_crm_low+1:k_crm_hig-1)*var1d_crm(n,k_crm_low+1:k_crm_hig-1)) ) &
                               / (dhgt_low+dhgt_hig+SUM(dhgt_crm(k_crm_low+1:k_crm_hig-1)))   !   range_bin
       enddo

   endif

 enddo INSTR_RANGE_LOOP


 return
 end subroutine refine_range_bin_bundle

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine refine_range_bin_bundle_sbm( mxlyr_crm, hgt_stag_crm, dhgt_crm, mxlyr_inst, hgt_stag_inst, &
                                         np,nbin, var1d_crm, var1d_inst ) 
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!   This routine modify radar exitinction and backscatter profile from 
!   CRM level to instrumental level.    (for single precision parameter) 
!   Interface of instrument range bin start from earth's ellipsoid.
!
! History:
!   05/2011  Toshi Matsui@NASA GSFC :: Initial
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)        :: mxlyr_crm                  ! maximum vertical layer of CRM
 real(sdsu_fps),intent(in) :: hgt_stag_crm(0:mxlyr_crm)  ! height leve at CRM interface [km] 
 real(sdsu_fps),intent(in) :: dhgt_crm(1:mxlyr_crm)      ! layer depth of CRM vertical height [km] 
 integer,intent(in)        :: mxlyr_inst                 ! maximum layer for radar instrument. 
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)! instrument measurement height at interface [km]
 integer,intent(in)        :: np,nbin
 real(sdsu_fps),intent(in) :: var1d_crm(1:mxlyr_crm,1:nbin,1:np)     ! parameter in CRM level []
 real(sdsu_fps),intent(out) :: var1d_inst(1:mxlyr_inst,1:nbin,1:np)  ! parameter in instrument level  []

 integer :: k_crm, k_inst !vertical loop indice
 integer :: k_crm_low, k_crm_hig ! 
 real(sdsu_fps) :: hgt_btm_inst, hgt_top_inst !layer interface height [km]
 real(sdsu_fps) :: range_bin ! instrument range bin height [km]
 real(sdsu_fps) :: dhgt_low, dhgt_hig
 real(sdsu_fps) :: dhgt_surf  !surface residual layer depth [km]

 integer :: n,l

!
! Looping of instrumental range
!
 INSTR_RANGE_LOOP: do k_inst = 1, mxlyr_inst

   hgt_btm_inst = hgt_stag_inst(k_inst-1)  ![km]
   hgt_top_inst = hgt_stag_inst(k_inst)    ![km]

   !
   ! range bin completely under the ground surface  --> all undefine
   !
! old code
!   if(hgt_btm_inst < hgt_stag_crm(0) .or. hgt_top_inst < hgt_stag_crm(0) ) then
!      var1d_inst (k_inst) = undefined
!      cycle INSTR_RANGE_LOOP
!   endif

! new code 
   !
   ! range bin completely under the ground surface  --> all undefine
   !
   if(hgt_btm_inst < hgt_stag_crm(0) .and. hgt_top_inst < hgt_stag_crm(0) ) then
      var1d_inst (k_inst,1:nbin,1:np) = undefined
      cycle INSTR_RANGE_LOOP
   endif

   !
   ! range bin completely above the CRM-defined top layer  --> all undefine
   !
   if(hgt_btm_inst >= hgt_stag_crm(mxlyr_crm) .or. hgt_top_inst >= hgt_stag_crm(mxlyr_crm) ) then
      var1d_inst (k_inst,1:nbin,1:np) = undefined
      cycle INSTR_RANGE_LOOP
   endif

   !
   ! range bin partially under the ground  --> surface index
   !
   if(hgt_btm_inst <= hgt_stag_crm(0) .and. hgt_top_inst >= hgt_stag_crm(0) ) then
      dhgt_surf = hgt_top_inst - hgt_stag_crm(0)  !surface residual layer depth [km]

      if( dhgt_surf/(hgt_top_inst-hgt_btm_inst) < 0.5e0 ) then !if residual layer < 50% of instr layer
            var1d_inst (k_inst,1:nbin,1:np) = undefined 
            cycle INSTR_RANGE_LOOP
      endif

      k_crm_low = 1
      dhgt_low  = hgt_stag_crm(1) - hgt_stag_crm(0)

   else
   !
   ! otherwise interpolate/extrapolate CRM-level value into instrumental range bin
   !

   ! get k_crm_low index
     do k_crm = 1, mxlyr_crm
        if( hgt_btm_inst >= hgt_stag_crm(k_crm-1) .and. hgt_btm_inst < hgt_stag_crm(k_crm) ) then
            k_crm_low = k_crm
            dhgt_low = hgt_stag_crm(k_crm) - hgt_btm_inst  ! [km] 
            exit
        endif
     enddo
   endif

   ! get k_crm_hig index
   do k_crm = 1, mxlyr_crm
      if( hgt_top_inst >= hgt_stag_crm(k_crm-1) .and. hgt_top_inst < hgt_stag_crm(k_crm) ) then
          k_crm_hig = k_crm
          dhgt_hig = hgt_top_inst - hgt_stag_crm(k_crm-1)    ! [km] 
          exit
      endif
   enddo

   !interplate
   if( k_crm_low == k_crm_hig ) then          !within same CRM layer

       var1d_inst(k_inst,1:nbin,1:np) = var1d_crm(k_crm_low,1:nbin,1:np)

   elseif( k_crm_low + 1 == k_crm_hig ) then  !just 1 vertical index lag

       range_bin = hgt_stag_inst(k_inst) - hgt_stag_inst(k_inst-1) ![km]

       var1d_inst(k_inst,1:nbin,1:np) = (dhgt_low * var1d_crm(k_crm_low,1:nbin,1:np) + dhgt_hig * var1d_crm(k_crm_hig,1:nbin,1:np) ) &
                               / (dhgt_low+dhgt_hig )   !dhgt_low+dhgt_hig= range_bin

   else  ! More than 1 vertical index lag. 

       range_bin = hgt_stag_inst(k_inst) - hgt_stag_inst(k_inst-1) ![km]

       do n = 1, np ; do l = 1, nbin
          var1d_inst(k_inst,l,n) = ( dhgt_low * var1d_crm(k_crm_low,l,n) + dhgt_hig * var1d_crm(k_crm_hig,l,n)  &
                                +SUM(dhgt_crm(k_crm_low+1:k_crm_hig-1)*var1d_crm( k_crm_low+1:k_crm_hig-1,l,n )) ) &
                               / (dhgt_low+dhgt_hig+SUM(dhgt_crm(k_crm_low+1:k_crm_hig-1)))   !   range_bin
       enddo  ; enddo

   endif

 enddo INSTR_RANGE_LOOP


 return
 end subroutine refine_range_bin_bundle_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine refine_range_bin_micro( mxlyr_crm, hgt_stag_crm, dhgt_crm, mxlyr_inst, hgt_stag_inst, &
                                     var1d_crm, var1d_inst ) 
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!   This routine modify radar exitinction and backscatter profile from 
!   CRM level to instrumental level.    (for single precision parameter) 
!   Interface of instrument range bin start from earth's ellipsoid.
!
! History:
!   05/2011  Toshi Matsui@NASA GSFC :: Initial
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)        :: mxlyr_crm                  ! maximum vertical layer of CRM
 real(sdsu_fps),intent(in) :: hgt_stag_crm(0:mxlyr_crm)  ! height leve at CRM interface [km] 
 real(sdsu_fps),intent(in) :: dhgt_crm(1:mxlyr_crm)      ! layer depth of CRM vertical height [km] 
 integer,intent(in)        :: mxlyr_inst                 ! maximum layer for radar instrument. 
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)! instrument measurement height at interface [km]
 real(sdsu_fps),intent(in) :: var1d_crm(1:mxlyr_crm)     ! parameter in CRM level []
 real(sdsu_fps),intent(out) :: var1d_inst(1:mxlyr_inst)  ! parameter in instrument level  []

 integer :: k_crm, k_inst !vertical loop indice
 integer :: k_crm_low, k_crm_hig ! 
 real(sdsu_fps) :: hgt_btm_inst, hgt_top_inst !layer interface height [km]
 real(sdsu_fps) :: range_bin ! instrument range bin height [km]
 real(sdsu_fps) :: dhgt_low, dhgt_hig
 real(sdsu_fps) :: dhgt_surf  !surface residual layer depth [km]

!
! Looping of instrumental range
!
 INSTR_RANGE_LOOP: do k_inst = 1, mxlyr_inst

   hgt_btm_inst = hgt_stag_inst(k_inst-1)  ![km]
   hgt_top_inst = hgt_stag_inst(k_inst)    ![km]

   !
   ! range bin completely under the ground surface  --> all undefine
   !
! old code
!   if(hgt_btm_inst < hgt_stag_crm(0) .or. hgt_top_inst < hgt_stag_crm(0) ) then
!      var1d_inst (k_inst) = undefined
!      cycle INSTR_RANGE_LOOP
!   endif

! new code 
   !
   ! range bin completely under the ground surface  --> all undefine
   !
   if(hgt_btm_inst < hgt_stag_crm(0) .and. hgt_top_inst < hgt_stag_crm(0) ) then
      var1d_inst (k_inst) = undefined
      cycle INSTR_RANGE_LOOP
   endif

   !
   ! range bin completely above the CRM-defined top layer  --> all undefine
   !
   if(hgt_btm_inst >= hgt_stag_crm(mxlyr_crm) .or. hgt_top_inst >= hgt_stag_crm(mxlyr_crm) ) then
      var1d_inst (k_inst) = undefined
      cycle INSTR_RANGE_LOOP
   endif

   !
   ! range bin partially under the ground  --> surface index
   !
   if(hgt_btm_inst <= hgt_stag_crm(0) .and. hgt_top_inst >= hgt_stag_crm(0) ) then
      dhgt_surf = hgt_top_inst - hgt_stag_crm(0)  !surface residual layer depth [km]

      if( dhgt_surf/(hgt_top_inst-hgt_btm_inst) < 0.5e0 ) then !if residual layer < 50% of instr layer
            var1d_inst (k_inst) = undefined 
            cycle INSTR_RANGE_LOOP
      endif

      k_crm_low = 1
      dhgt_low  = hgt_stag_crm(1) - hgt_stag_crm(0)

   else
   !
   ! otherwise interpolate/extrapolate CRM-level value into instrumental range bin
   !

   ! get k_crm_low index
     do k_crm = 1, mxlyr_crm
        if( hgt_btm_inst >= hgt_stag_crm(k_crm-1) .and. hgt_btm_inst < hgt_stag_crm(k_crm) ) then
            k_crm_low = k_crm
            dhgt_low = hgt_stag_crm(k_crm) - hgt_btm_inst  ! [km] 
            exit
        endif
     enddo
   endif

   ! get k_crm_hig index
   do k_crm = 1, mxlyr_crm
      if( hgt_top_inst >= hgt_stag_crm(k_crm-1) .and. hgt_top_inst < hgt_stag_crm(k_crm) ) then
          k_crm_hig = k_crm
          dhgt_hig = hgt_top_inst - hgt_stag_crm(k_crm-1)    ! [km] 
          exit
      endif
   enddo

   !interplate
   if( k_crm_low == k_crm_hig ) then          !within same CRM layer

       var1d_inst(k_inst) = var1d_crm(k_crm_low)

   elseif( k_crm_low + 1 == k_crm_hig ) then  !just 1 vertical index lag

       range_bin = hgt_stag_inst(k_inst) - hgt_stag_inst(k_inst-1) ![km]

       var1d_inst(k_inst) = (dhgt_low * var1d_crm(k_crm_low) + dhgt_hig * var1d_crm(k_crm_hig) ) &
                               / (dhgt_low+dhgt_hig )   !dhgt_low+dhgt_hig= range_bin

   else  ! More than 1 vertical index lag. 

       range_bin = hgt_stag_inst(k_inst) - hgt_stag_inst(k_inst-1) ![km]

       var1d_inst(k_inst) = ( dhgt_low * var1d_crm(k_crm_low) + dhgt_hig * var1d_crm(k_crm_hig)  &
                                +SUM(dhgt_crm(k_crm_low+1:k_crm_hig-1)*var1d_crm(k_crm_low+1:k_crm_hig-1)) ) &
                               / (dhgt_low+dhgt_hig+SUM(dhgt_crm(k_crm_low+1:k_crm_hig-1)))   !   range_bin

   endif

 enddo INSTR_RANGE_LOOP


 return
 end subroutine refine_range_bin_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine refine_range_bin_lidar( nln, hgt_stag_visir, dhgt_visir, mxlyr_inst, hgt_stag_inst, &
                                     var1d_visir, var1d_inst_visir ) 
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!   This routine modify lidar exitinction and backscatter profile from 
!   CRM level to instrumental level.    (for doubple precision input) 
!   level mxlyr_crm
!   Interface of instrument range bin start from earth's ellipsoid.
!
! History:
!   05/2011  Toshi Matsui@NASA GSFC :: Initial
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)        :: nln                        ! maximum vertical layer of CRM
 real(sdsu_fpd),intent(in) :: hgt_stag_visir(1:nln+1)      ! height leve at CRM interface [km] 
 real(sdsu_fpd),intent(in) :: dhgt_visir(1:nln)        ! layer depth of CRM vertical height [km] 
 integer,intent(in)        :: mxlyr_inst                 ! maximum layer for radar instrument. 
 real(sdsu_fps),intent(in) :: hgt_stag_inst(0:mxlyr_inst)! instrument measurement height at interface [km]
 real(sdsu_fpd),intent(in) :: var1d_visir(1:nln)           ! parameter in CRM level []
 real(sdsu_fpd),intent(out) :: var1d_inst_visir(1:mxlyr_inst)  ! parameter in instrument level  []

 integer        :: mxlyr_crm                  ! maximum vertical layer of CRM
 real(sdsu_fps) :: hgt_stag_crm(0:nln)  ! height leve at CRM interface [km] 
 real(sdsu_fps) :: dhgt_crm(1:nln)      ! layer depth of CRM vertical height [km] 
 real(sdsu_fpd) :: var1d_crm(1:nln)           ! parameter in CRM level []
 real(sdsu_fpd) :: var1d_inst(1:mxlyr_inst)  ! parameter in instrument level  []

 integer :: n
 integer :: k_crm, k_inst, k_rev !vertical loop indice
 integer :: k_crm_low, k_crm_hig ! 
 real(sdsu_fps) :: hgt_btm_inst, hgt_top_inst !layer interface height [km]
 real(sdsu_fps) :: range_bin ! instrument range bin height [km]
 real(sdsu_fps) :: dhgt_low, dhgt_hig


!
! from opt_visir, vertical level is opposite. 1 is TOA, nln is BOA
!
 mxlyr_crm = nln

!
! reverse vertical level
!
 do k_crm = 1 , mxlyr_crm
    n = nln - k_crm + 1
    dhgt_crm(k_crm) = dhgt_visir(n)
    var1d_crm(k_crm) = var1d_visir(n) 
enddo

 do k_crm = 0 , mxlyr_crm
    n = nln - k_crm + 1
    hgt_stag_crm(k_crm) = hgt_stag_visir(n)
 enddo


!
! Looping of instrumental range
!
 INSTR_RANGE_LOOP: do k_inst = 1, mxlyr_inst

   hgt_btm_inst = hgt_stag_inst(k_inst-1)
   hgt_top_inst = hgt_stag_inst(k_inst)

   !
   ! rang bin under the ground  --> all undefine
   !
   if(hgt_btm_inst < hgt_stag_crm(0) .or. hgt_top_inst < hgt_stag_crm(0) ) then
      var1d_inst (k_inst) = undefined
      cycle INSTR_RANGE_LOOP
   endif

   !
   ! rang bin above the CRM-defiined TOA  --> all undefine
   !
   if(hgt_btm_inst >= hgt_stag_crm(mxlyr_crm) .or. hgt_top_inst >= hgt_stag_crm(mxlyr_crm) ) then
      var1d_inst (k_inst) = undefined
      cycle INSTR_RANGE_LOOP
   endif

   !
   ! otherwise interpolate/extrapolate CRM-level value into instrumental range bin
   !

   ! get k_crm_low index
   do k_crm = 1, mxlyr_crm
      if( hgt_btm_inst >= hgt_stag_crm(k_crm-1) .and. hgt_btm_inst < hgt_stag_crm(k_crm) ) then
          k_crm_low = k_crm
          dhgt_low = hgt_stag_crm(k_crm) - hgt_btm_inst  ! [km] 
          exit
      endif
   enddo

   ! get k_crm_hig index
   do k_crm = 1, mxlyr_crm
      if( hgt_top_inst >= hgt_stag_crm(k_crm-1) .and. hgt_top_inst < hgt_stag_crm(k_crm) ) then
          k_crm_hig = k_crm
          dhgt_hig = hgt_top_inst - hgt_stag_crm(k_crm-1)    ! [km] 
          exit
      endif
   enddo

   !interplate
   !interplate
   if( k_crm_low == k_crm_hig ) then          !within same CRM layer

       var1d_inst(k_inst) = var1d_crm(k_crm_low)

   elseif( k_crm_low + 1 == k_crm_hig ) then  !just 1 vertical index lag

       range_bin = hgt_stag_inst(k_inst) - hgt_stag_inst(k_inst-1) ![km]

       var1d_inst(k_inst) = (DBLE(dhgt_low) * var1d_crm(k_crm_low) + DBLE(dhgt_hig) * var1d_crm(k_crm_hig) ) &
                               / DBLE(range_bin)

   else  ! More than 1 vertical index lag.   (toshii verify lator)

       range_bin = hgt_stag_inst(k_inst) - hgt_stag_inst(k_inst-1) ![km]

       var1d_inst(k_inst) = ( DBLE(dhgt_low) * var1d_crm(k_crm_low) + DBLE(dhgt_hig) * var1d_crm(k_crm_hig)  &
                                +SUM(DBLE(dhgt_crm(k_crm_low+1:k_crm_hig-1))*var1d_crm(k_crm_low+1:k_crm_hig-1)) ) &
                               / DBLE(range_bin)

   endif

 enddo INSTR_RANGE_LOOP


!
! reverse vertical level again consistent to opt_visir
!
 do k_inst = 1, mxlyr_inst
    k_rev = mxlyr_inst - k_inst + 1

    var1d_inst_visir(k_rev) = var1d_inst(k_inst) 
 enddo

 return
 end subroutine refine_range_bin_lidar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine sbm_path( dhgt, n_bin, water_path_out )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Derive cloud/rain ice/liquid water path linking to L2 data
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real,intent(in) :: dhgt(mxlyr)               ! thickness of layer [km]
 type( particle_sbm ),intent(in) :: n_bin(mxlyr, nbin) ! number density [#/m4]
 real,intent(out) :: water_path_out(4)  ! 1-cloud liq, 2-rain liq, 3-cloud ice, 4-rain ice path [kg/m2]

 real :: water_path(4) !1-cloud liq, 2-rain liq, 3-cloud ice, 4-rain ice path [kg/m2]

 type( particle_sbm) :: path(1:2) !rain rate [mm/sec]
 integer :: k

 water_path = 0.

!
! integrate vertically
!
 do k = 1, mxlyr

   path(1:2)%liq     = 0. ; path(1:2)%ice_col = 0. ; path(1:2)%ice_pla = 0.
   path(1:2)%ice_den = 0. ; path(1:2)%snow    = 0. ; path(1:2)%graupel = 0. 
   path(1:2)%hail    = 0. 

 call sbm_path_layer(dhgt(k), n_bin(k,1:nbin)%liq,     rad_sbm%liq,     drad_sbm%liq,     x_sbm%liq,     path(1:2)%liq     )
 call sbm_path_layer(dhgt(k), n_bin(k,1:nbin)%ice_col, rad_sbm%ice_col, drad_sbm%ice_col, x_sbm%ice_col, path(1:2)%ice_col )
 call sbm_path_layer(dhgt(k), n_bin(k,1:nbin)%ice_pla, rad_sbm%ice_pla, drad_sbm%ice_pla, x_sbm%ice_pla, path(1:2)%ice_pla )
 call sbm_path_layer(dhgt(k), n_bin(k,1:nbin)%ice_den, rad_sbm%ice_den, drad_sbm%ice_den, x_sbm%ice_den, path(1:2)%ice_den )
 call sbm_path_layer(dhgt(k), n_bin(k,1:nbin)%snow,    rad_sbm%snow,    drad_sbm%snow,    x_sbm%snow,    path(1:2)%snow    )
 call sbm_path_layer(dhgt(k), n_bin(k,1:nbin)%graupel, rad_sbm%graupel, drad_sbm%graupel, x_sbm%graupel, path(1:2)%graupel )
 call sbm_path_layer(dhgt(k), n_bin(k,1:nbin)%hail,    rad_sbm%hail,    drad_sbm%hail,    x_sbm%hail,    path(1:2)%hail    )

 water_path(1) = water_path(1) + path(1)%liq  !cloud water path [kg/m2]
 water_path(2) = water_path(2) + path(2)%liq  !rain water path [kg/m2]
 water_path(3) = water_path(3) + path(1)%ice_col + path(1)%ice_pla + path(1)%ice_den &  !cloud ice path [kg/m2]
                               + path(1)%snow    + path(1)%graupel + path(1)%hail
 water_path(4) = water_path(4) + path(2)%ice_col + path(2)%ice_pla + path(2)%ice_den &  !rain ice path [kg/m2]
                               + path(2)%snow    + path(2)%graupel + path(2)%hail

 enddo

 water_path_out(1:4) = water_path(1:4)

 return
 end subroutine sbm_path

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine sbm_path_layer( dhgt,  ns, r, dr, x, path_out )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Derive various L2 2dsurface data to be convolved in later.   
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real,intent(in) :: dhgt       ! layer depth [km]
 real,intent(in) :: ns (nbin)  ! number density [#/m4]
 real,intent(in) :: r (nbin)   ! particle radius [cm]
 real,intent(in) :: dr (nbin)  ! bin width [cm]
 real,intent(in) :: x  (nbin)  ! mass per particle [g]
 real,intent(out) :: path_out(2)  ! water path within layer[kg/m2]

 integer :: n       ! looping
 real :: path(1:2)  ! accumulating rain rate [mm/hr]
 real :: q          ! mixing ratio [g/m3]

 path = 0.     ! initialize

 do n = 1, nbin

    q = ns(n) * dr(n)*1.e-2 * x(n)    ! mixing ratio [g/m3] 

    if( r(n)*1e+4 < 100. ) then ! cloud particle (radius < 100 um)
       path(1) = path(1) + q * 1.e-3 * dhgt * 1.e+3 ! [kg/m2]
    else              ! rain particle
       path(2) = path(2) + q * 1.e-3 * dhgt * 1.e+3  ! [kg/m2]
    endif

 enddo

 path_out(1:2) = path(1:2)  ! [kg/m2]

 return
 end subroutine sbm_path_layer

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

! subroutine gce_rain( press, q, re, rain_liq, rain_ice )
! implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation techniqu.
!  This routine is for GCE microphysics.
!    
! History: 
!  06/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
! real(sdsu_fps),intent(in)::press  ! layer pressure [hPa]

 !type ( particle_gce ),intent(in) :: & ! particle_gen is defined in module_simulater
 !    q,  & ! mixing ratio [g/m3] 
 !    re    ! effective radius  [micron]

! real(sdsu_fps),intent(out) :: rain_liq, rain_ice ! [mm/hr]
! type( particle_gce) :: rain !rain rate [mm/sec]


! call bulk_rain_rate(press, q%rain, re%rain, mu_gen%rain, a_vt%rain, b_vt%rain, vt_max%rain, rho_gce%rain, rain%rain )

! call bulk_rain_rate(press, q%snow, re%snow, mu_gen%snow, a_vt%snow, b_vt%snow, vt_max%snow, rho_gce%snow, rain%snow )

! call bulk_rain_rate(press, q%graupel, re%graupel, mu_gen%graupel, a_vt%graupel, b_vt%graupel, vt_max%graupel, &
!                      rho_gce%graupel, rain%graupel )

! call bulk_rain_rate(press, q%hail, re%hail, mu_gen%hail, a_vt%hail, b_vt%hail, vt_max%hail, rho_gce%hail, rain%hail )

! rain_liq = rain%rain
! rain_ice = rain%snow + rain%graupel + rain%hail


! return
! end subroutine gce_rain

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gen_rain( press, q, re, rain_liq, rain_ice )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:
!  compute instanteneous rain rate profile for generalized bulk microphysics scheme.   
!    
! History: 
!  06/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in)::press  ! layer pressure [hPa]

 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 real(sdsu_fps),intent(out) :: rain_liq, rain_ice ! [mm/hr]

 type( particle_gen) :: rain !rain rate [mm/sec]


 call bulk_rain_rate(press, q%rain, re%rain, mu_gen%rain, a_vt%rain, b_vt%rain, vt_max%rain, rho_gen%rain, rain%rain )

 call bulk_rain_rate(press, q%snow, re%snow, mu_gen%snow, a_vt%snow, b_vt%snow, vt_max%snow, rho_gen%snow, rain%snow )

 call bulk_rain_rate(press, q%graupel, re%graupel, mu_gen%graupel, a_vt%graupel, b_vt%graupel, vt_max%graupel, &
                      rho_gen%graupel, rain%graupel )

 call bulk_rain_rate(press, q%hail, re%hail, mu_gen%hail, a_vt%hail, b_vt%hail, vt_max%hail, rho_gen%hail, rain%hail )

 rain_liq = rain%rain
 rain_ice = rain%snow + rain%graupel + rain%hail

 return
 end subroutine gen_rain

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine bulk_rain_rate(p, q, re, mu, a_vt, b_vt, vt_max, dens, rain_out)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute single-scattering properties from atmospheric particles using LUTs-interpolation techniqu.
!  This routine is for GCE microphysics.
!    
! History: 
!  06/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!---------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(in) :: p      ! pressure [mb]
 real(sdsu_fps),intent(in) :: q      ! mass mixing ratio [g/m3]
 real(sdsu_fps),intent(in) :: re     ! effective radius  [micron]
 real(sdsu_fps),intent(in) :: mu     ! dispersion parameters [-]
 real(sdsu_fps),intent(in) :: a_vt   ! alpha of Vt [-]
 real(sdsu_fps),intent(in) :: b_vt   ! beta of Vt [-]
 real(sdsu_fps),intent(in) :: vt_max ! maximum limit of Vt [m/s]
 real(sdsu_fps),intent(in) :: dens   ! density of particle [kg/m3]

 real(sdsu_fps),intent(out) :: rain_out ! rain rate [mm/hr]

 integer :: i    !looping 
 integer :: imax !max for looping

 real(sdsu_fps) :: n0   ! intercept for the exponetial DSD [1/m**4]
 real(sdsu_fps) :: rad  ! radius of particle [mm]
 real(sdsu_fps) :: lam  ! the slope of the distribution [1/m]
 real(sdsu_fps) :: num  ! a particle number density per radius increment [1/m**4]
 real(sdsu_fps) :: gfac1, gfac2, gamfac  !gamma parameter
 real(sdsu_fps) :: &
   q_inc      ,& ! mass mixing ratio [kg/m3]
   rain_inc   ,& ! rainfall rate [mm/hr]
   factor     ,& ! pressure/density factor
   vt         ,& ! terminal velocity [m/s]
   d_mm       ,& ! particle diameter [mm]
   d_m        ,& ! particle diameter [m]
   d_increment,& ! diameter increment over size loops [mm]
   d_max         ! maximum diameter of particle


 real(sdsu_fps),parameter :: dr = 0.05    ! increment bin of radius [mm]
 real(sdsu_fps),parameter :: densice = 0.917e+3  !solid ice density
 real(sdsu_fps),parameter :: densliq = 1.0e+3    !liquid density
 real(sdsu_fps),parameter :: wc_unit = 1.0  !unit water content 1 [g/m3]


 if(q .lt. q_min_condensate) then
     rain_out = 0. 
     return
 endif

 if( re == 0. ) then
    rain_out = 0. ; return
 endif

 factor = sqrt(800.e0/p)  !cetered on 800mb pressure 


 d_increment = 0.10  ! [mm]
 d_max = 20.
! define imax for a givne radius increment 20. is maximum diameter
 imax = nint(d_max/d_increment)

!
! derive lambda and intercept of expoential DSD for unit water content
!
  call gamma_reff(mu+4.0,gfac1)
  call gamma_reff(mu+3.0,gfac2)
  gamfac = gfac1/gfac2

!
! derive lambda
!
  lam = 1. / (2.*re*1.e-6) * gamfac  !slope [1/m]

  if( lam <= const_smallest_real**(1./(4.+mu)) .or. dens ==0. ) then
      rain_out = 0. ; return
  endif
 
!
! derive intercept
!       [g/m3]                     [m3/g]        [1/m]
  n0  = wc_unit * 6. / const_pi / (dens*1000.) * (lam**(4.+mu))  / gfac1  ! [1/m**(4+mu)]

!
! size loop
!
  rain_inc = 0.
  do i=0,imax

     d_mm = d_increment * 0.5 + d_increment * FLOAT(i)  !diameter [mm]
     d_m  = d_mm * 1e-3                        ! diameter but using different unit [m]

!
! Get N(D)  [1/m4]
!         [1/m**(4+mu)]  [ m**mu ]      [1/m]  [m]   = [1/m4]
     num = n0           * (d_m**mu) * exp(-lam * d_m )


!
! derive mixing ratio [kg/m3]
!
     q_inc = num * (d_increment * 1.e-3)  * const_pi / 6.e0 * (d_m*d_m*d_m) * (dens)
!           [1/m4]  * [m]                                      [m3]       *  [kg/m3]    = [kg/m3]

!
! terminal velocity [m/s]
!
     vt    = min( vt_max, a_vt * (d_m**b_vt) * factor )

!
! rainfall rate
!
     rain_inc = rain_inc + q_inc * vt * 3600.e0  !rain rate [mm/hr] per bin per uniq water content
  enddo

!
! ouptut total rainfall rate [mm/hr]
!
 rain_out = rain_inc  * q

 return
 end subroutine bulk_rain_rate

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine sbm_rain( press, n_bin, frime ,rain_liq, rain_ice)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Derive L2 surface ice/liquid rain.   
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real,intent(in) :: press               ! pressure [hPa]
 type( particle_sbm ),intent(in) :: n_bin(1:nbin) ! number density [#/m4]
 real,intent(in) :: frime(1:nbin)       ! riming fraction for snow [-]
 real,intent(out) :: rain_liq, rain_ice ! [mm/hr]
 type( particle_sbm) :: rain !rain rate [mm/sec]

!
! get rainfall for each species
!
 call sbm_rain_rate(press, vt_sbm%liq,     n_bin%liq,     drad_sbm%liq,     x_sbm%liq,     rain%liq    )

 call sbm_rain_rate(press, vt_sbm%ice_col, n_bin%ice_col, drad_sbm%ice_col, x_sbm%ice_col, rain%ice_col)

 call sbm_rain_rate(press, vt_sbm%ice_pla, n_bin%ice_pla, drad_sbm%ice_pla, x_sbm%ice_pla, rain%ice_pla)

 call sbm_rain_rate(press, vt_sbm%ice_den, n_bin%ice_den, drad_sbm%ice_den, x_sbm%ice_den, rain%ice_den)

 call sbm_rain_rate(press, vt_sbm%snow,    n_bin%snow,    drad_sbm%snow,    x_sbm%snow,    rain%snow   )

 call sbm_rain_rate(press, vt_sbm%graupel, n_bin%graupel, drad_sbm%graupel, x_sbm%graupel, rain%graupel)

 call sbm_rain_rate(press, vt_sbm%hail,    n_bin%hail,    drad_sbm%hail,    x_sbm%hail,    rain%hail   )

!
! output
!
 rain_liq = rain%liq        ! liquid rain [mm/hr]
 rain_ice = rain%ice_col+rain%ice_pla+rain%ice_den+rain%snow+rain%graupel+rain%hail !ice rain [mm/hr]

 return
 end subroutine sbm_rain

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine sbm_rain_rate( p, vts, ns, dr, x, rain_out )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Derive various L2 2dsurface data to be convolved in later.   
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real,intent(in) :: p          ! pressure [mb]
 real,intent(in) :: vts(nbin)  ! terminal velcoty [m/s]
 real,intent(in) :: ns (nbin)  ! number density [#/m4]
 real,intent(in) :: dr (nbin)  ! bin width [cm]
 real,intent(in) :: x  (nbin)  ! mass per particle [g]
 real,intent(out) :: rain_out  ! rain rate per hour[mm/hr]

 integer :: n  ! looping
 real :: vt    ! actual terminal velocity for a given pressure [m/s]
 real :: rain  ! accumulating rain rate [mm/hr]
 real :: q     ! mixing ratio [g/m3]

 rain = 0.     ! initialize

 do n = 1, nbin
    vt      = 1.e-2*vts(n)*sqrt(1000.e0/p)  ! actual terminal velocity for a given pressure [m/s]
    q       = ns(n) * dr(n)*1.e-2 * x(n)    ! mixing ratio [g/m3] 
    rain    = rain + q * 1.e-3 * vt * 3600.e0  !rain rate [mm/hr]
 enddo

 rain_out = rain  ![mm/hr]

 return
 end subroutine sbm_rain_rate

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_extra_memory_domain(dij, is,ie,js,je,ibs,ibe,jbs,jbe)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  This routine estimate estimate extra-buffer memory domain larger than actual memory domain. 
!  This buffer zone is required for slant-path radiative transfer options for scan_simulator. 
!
! History:
! 01/2012  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)  :: dij             ! extra buffere
 integer,intent(in)  :: is,ie,js,je     ! memory bounds
 integer,intent(out) :: ibs,ibe,jbs,jbe !extended buffer bounds
 integer :: ierr
 integer :: irank, imyrank
 integer :: totcnt

 if( is==1 .and. ie==mxgridx .and. &
     js==1 .and. je==mxgridy )  then

     !
     ! single processor case
     !
     ibs=1 ; ibe=mxgridx
     jbs=1 ; jbe=mxgridy
     totcnt = (ibe-ibs+1)*(jbe-jbs+1)
 
 else !domain decomposition case

#if MPI == 2

    ibs = max(       1, is - dij )
    ibe = min( mxgridx, ie + dij )
    jbs = max(       1, js - dij )
    jbe = min( mxgridy, je + dij )
    totcnt = (ibe-ibs+1)*(jbe-jbs+1)

!
! distribute for all processor
!
    if(.not. allocated(buffer)) allocate(buffer(0:numproc_tot-1))

!
! distribute defined-only vector size to other memory. 
!
    call MPI_ALLGATHER(      ibs, 1, mpi_integer,                  &
                       buffer%is, 1, mpi_integer, mpi_comm_world, ierr)

    call MPI_ALLGATHER(      ibe, 1, mpi_integer,                  &
                       buffer%ie, 1, mpi_integer, mpi_comm_world, ierr)

    call MPI_ALLGATHER(      jbs, 1, mpi_integer,                  &
                       buffer%js, 1, mpi_integer, mpi_comm_world, ierr)

    call MPI_ALLGATHER(      jbe, 1, mpi_integer,                  &
                       buffer%je, 1, mpi_integer, mpi_comm_world, ierr)

    call MPI_ALLGATHER(   totcnt, 1, mpi_integer,                  &
                       buffer%totnum, 1, mpi_integer, mpi_comm_world, ierr)


!
! Calculate overlap (memory and buffer) domain bounds and total grid numbers
! (these info will be used in shapring info across processors)
   if(.not. allocated(overlap))  allocate( overlap(0:numproc_tot-1, 0:numproc_tot-1 ) )

   do imyrank = 0, numproc_tot-1
     do irank = 0, numproc_tot-1

        call buffer_overlapp( memory(imyrank), buffer(irank) , overlap(imyrank,irank) )
!if(masterproc) print*, 'imyrank=', imyrank, 'irank=',irank
!if(masterproc) print*, 'm=',memory(imyrank)
!if(masterproc) print*, 'b=',buffer(irank)
!if(masterproc) print*, 'o=',overlap(imyrank,irank)

     enddo
   enddo

!call stop_sdsu(myrank,'')



#endif

 endif

 end subroutine get_extra_memory_domain

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_sarta
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  write out SARTA-simulated AIRS Tb in netcdf format
!
! History:
! 11/2012  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer :: i,j,k,irec,nch
 character(len=20) :: typ
 integer,parameter :: io = 121
 integer :: is,ie,js,je,ks,ke

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end

!
! write out CRM surface input for grads binarry format
!

!
! =============open output file name, and prepare header ==========================
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.SARTA_AIRS'&
               //trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
  print("(a)"),'MSG write_out_sarta: write SARTA (AIRS) output file in grads format ',&
  ' ->',trim(sdsu_io_file)

!
! Open output file
!
#if MPI == 2
 if(masterproc) &
#endif
 open(unit= io, file = sdsu_io_file, access='direct', status='replace',recl = mxgridx*mxgridy*4 )
 irec = 1

 do nch = 1, nchan_sarta
    call dump( sarta_out(is:ie,js:je,nch) , io, irec ) ! AIRS Tb [K]
 enddo

#if MPI == 2
 if(masterproc) &
#endif
 close (io)

 if(masterproc) print*,''
    
!
!  write out grads control file
!
#if MPI == 2
 if(masterproc) then
#endif

 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.SARTA_AIRS'//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.SARTA_AIRS'//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  AIRS Tb (Zdimension is used for channel dimension)'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef    ',nchan_sarta,' linear 0. 1'
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',1

    if( airs_unit == 'tb' ) then
       write(io,*) 'tb   ',nchan_sarta,'  0  AIRS brightness temperature [K]'
    elseif(  airs_unit == 'rad' ) then
       write(io,*) 'rad   ',nchan_sarta,'  0  AIRS radiance [W/m2/cm-1]'
    endif

    write(io,*) 'endvars '
    close(io)
 endif


!
! write out sarta AIRS channel ID and frequncy 
!
 if(masterproc) print*,''
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.SARTA_AIRS_CHANNEL'//trim(output_suffix)//'.txt'
 print("(a)"),'see AIRS CHANNEL ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 write(io,*) ' Chan ID    Chan Frquency [cm-1]'
 do nch = 1, nchan_sarta
   write(io,FMT='(I9,7x,F9.3)') sarta_chan_id(nch), sarta_chan_freq(nch)
 enddo
 close(io)



#if MPI == 2
 endif
#endif


 return
 end subroutine write_out_sarta

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_sarta_nc
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  write out SARTA-simulated AIRS Tb in netcdf format
!
! History:
! 11/2012  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer :: i,j,k,nf,nch
 integer :: is,ie,js,je,ks,ke, im

!netcdf
 integer :: ncid, id, ich,np
 integer :: i_dimid, j_dimid, k_dimid, n_dimid
 integer, parameter :: ndims3d = 3, ndims2d = 2, ndims1d=1
 integer :: dimids0dchan(ndims1d), dimids2d(ndims2d), dimids3d(ndims3d)
 integer :: id_chan, id_freq, id_lon, id_lat, id_tskin, id_tair, id_wv  ! var IDs 
 character(len=100) :: description

 integer,allocatable :: id_out(:)
 real(sdsu_fps),allocatable :: var3d(:,:,:)
 character(len=4) :: cccc !channel id 
 character(len=4) :: ffff !channel frequency 
 character(len=3) :: ddd  !channel digit

!
! initialize array
!
 is=myi_start ; ie=myi_end ;  js=myj_start ; je=myj_end  ; ks=1 ; ke=mxlyr


 if( .not. allocated(id_out) ) then
     allocate( id_out(1:nchan_sarta),&
               var3d(is:ie,js:je,ks:ke) )
 endif

!
! =============open output file name, and prepare header ==========================
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.SARTA_AIRS'&
               //trim(output_suffix)//'.nc'

#if MPI == 2
 if(masterproc)  then
#endif

  print("(a)"),'MSG write_out_sarta_nc: write SARTA (AIRS) output file in netcdf format ',&
  ' ->',trim(sdsu_io_file)



! call check( nf90_create(trim(sdsu_io_file), IOR(4,nf90_64bit_offset) , ncid) )  !64bit output (big output data

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )


!
! Define dimensions
!
 call check( nf90_def_dim(ncid, 'west_east'  , mxgridx, i_dimid) )
 call check( nf90_def_dim(ncid, 'south_north', mxgridy, j_dimid) )
 call check( nf90_def_dim(ncid, 'bottom_top' , mxlyr  , k_dimid) )
 call check( nf90_def_dim(ncid, 'channel' , nchan_sarta , n_dimid) )

 dimids0dchan = (/ n_dimid /)
 dimids2d     = (/ i_dimid, j_dimid /)
 dimids3d     = (/ i_dimid, j_dimid, k_dimid /)

 call check( nf90_def_var(ncid, 'Channel_ID', NF90_INT, dimids0dchan, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(channel)') )
 call check( nf90_put_att(ncid, id, 'description',' Channel ID ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_chan = id

 call check( nf90_def_var(ncid, 'Frequency', NF90_REAL, dimids0dchan, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(channel)') )
 call check( nf90_put_att(ncid, id, 'description',' Channel Center Frequency ') )
 call check( nf90_put_att(ncid, id, 'units', '[cm-1]') )
 id_freq = id

 call check( nf90_def_var(ncid, 'lat', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(i,j)') )
 call check( nf90_put_att(ncid, id, 'description',' Latitude  ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_lat = id

 call check( nf90_def_var(ncid, 'lon', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(i,j)') )
 call check( nf90_put_att(ncid, id, 'description',' Longitude  ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_lon = id

 call check( nf90_def_var(ncid, 'skintemp', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(i,j)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Skin Temperature ') )
 call check( nf90_put_att(ncid, id, 'units', '[K]') )
 id_tskin = id

 call check( nf90_def_var(ncid, 'airtemp', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(i,j,k)') )
 call check( nf90_put_att(ncid, id, 'description',' Air Temperature ') )
 call check( nf90_put_att(ncid, id, 'units', '[K]') )
 id_tair = id

 call check( nf90_def_var(ncid, 'watervapor', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(i,j,k)') )
 call check( nf90_put_att(ncid, id, 'description',' Water vapor mass mixing ratio ') )
 call check( nf90_put_att(ncid, id, 'units', '[gm/kg dry air]') )
 id_wv = id


 do nch = 1, nchan_sarta
    write(cccc,"(I4.4)") sarta_chan_id(nch)  !channel ID
    write(ffff,"(I4.4)") INT(sarta_chan_freq(nch)) !channel center frequency four integer digit.
    write(ddd,"(I3.3)") INT(1.e3*(sarta_chan_freq(nch) - REAL(INT(sarta_chan_freq(nch))))) ! three decimal digit
 
    if( airs_unit == 'tb' ) then
      call check( nf90_def_var(ncid, 'tb_'//cccc, NF90_REAL, dimids2d, id) )
      call check( nf90_put_att(ncid, id, 'MemoryOrder','(i,j)') )
      call check( nf90_put_att(ncid, id, 'description',' AIRS Tb at '//ffff//'.'//ddd//' [cm-1]' ) )
      call check( nf90_put_att(ncid, id, 'units', '[K]') )
    elseif(  airs_unit == 'rad' ) then
      call check( nf90_def_var(ncid, 'rad_'//cccc, NF90_REAL, dimids2d, id) )
      call check( nf90_put_att(ncid, id, 'MemoryOrder','(i,j)') )
      call check( nf90_put_att(ncid, id, 'description',' AIRS radiance at '//ffff//'.'//ddd//' [cm-1]' )  )
      call check( nf90_put_att(ncid, id, 'units', '[W/m2/cm-1]') )
    else
      stop 'MSG write_out_sarta_nc: no such airs_unit'
    endif

    id_out(nch) = id

 enddo


 call check( nf90_enddef(ncid) )  ! End define mode.

!
! ----------------------------- DUMP DATA ------------------------------------------
!

!
! dump channel information
!
  call check( nf90_put_var(ncid, id_chan , sarta_chan_id  (1:nchan_sarta)) )
  call check( nf90_put_var(ncid, id_freq , sarta_chan_freq(1:nchan_sarta)) )

#if MPI == 2
 endif
#endif

!
! following parameter are tile domain in MPI=2 case. 
!
   call dumpnc( surface(is:ie,js:je)%lat    , id_lat  , ncid )
   call dumpnc( surface(is:ie,js:je)%lon    , id_lon  , ncid )
   call dumpnc( surface(is:ie,js:je)%t_skin , id_tskin, ncid )
   call dumpnc( atmos(is:ie,js:je,ks:ke)%t_air , id_tair, ncid )

   do k = ks,ke ; do j = js,je ; do i = is,ie 
      var3d(i,j,k) = atmos(i,j,k)%sh / (1.e0 - atmos(i,j,k)%sh ) * 1.e-3 !water vapor mass mixing ratio [gm/kg]
   enddo ; enddo ; enddo
   call dumpnc( var3d(is:ie,js:je,ks:ke) , id_wv, ncid )

   do nch = 1, nchan_sarta
      call dumpnc( sarta_out(is:ie,js:je,nch), id_out(nch), ncid )
   enddo


#if MPI == 2
   if(masterproc) &
#endif
  call check( nf90_close(ncid) )  ! close nc file


 return
 end subroutine write_out_sarta_nc

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_gprof
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!
! History:
! 04/2012  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer :: i,j,k,nf
 integer :: is,ie,js,je,ks,ke, im
 real(sdsu_fps) :: bw_deg
 real(sdsu_fps),allocatable,save :: bw_array(:) ! beam width [deg]
 real(sdsu_fps),allocatable,save :: wgt(:,:)    ! weight 
 real(sdsu_fps),allocatable,save :: gf(:,:,:)     ! gain function [-]
 logical,allocatable,save :: in_ifov(:,:,:)     !within IFOV or not
 real(sdsu_fps) :: dx, dx_intercept
 real(sdsu_fpd) :: wgt_dble

!netcdf
 integer :: ncid, id, ich,np
 integer :: i_dimid, j_dimid, k_dimid, k_radar_dimid
 integer, parameter :: ndims2d = 2, ndims1d=1
 integer :: dimids1d(ndims1d), dimids2d(ndims2d), dimids2d_radar(ndims2d)
 integer :: id_data, id_lat, id_lon !var IDs 
 integer :: id_year, id_month, id_day, id_hour, id_min, id_sec !var IDs (time)
 integer :: id_rain, id_rf, id_snow, id_snowd, id_tb(max_chan), id_tb37(max_chan)
 integer :: id_wind, id_skin, id_t2m, id_tcwv, id_u1km, id_v1km, id_elev
 integer :: id_rwc, id_lwc, id_mwc, id_iwc, id_lh, id_dhgt
 integer :: id_zeku, id_zeka, id_zew1,id_zew2,id_zew3, id_piaw1,id_piaw2,id_piaw3
 character(len=100) :: description
 character :: pol*3, cc*2

 real(sdsu_fps) :: incidence_ang ! incidence angle [rad]
 real(sdsu_fps) :: dist_sat_fov ! distance between satellite and FOV [km]

 real(sdsu_fps),allocatable :: var1d(:) , var2d(:,:), tb(:,:), temp2d(:,:), tb1d(:), precip1d(:)  !temporal array
 integer,allocatable :: ivar1d(:)  !temporal array
 logical :: cycle_this
 integer :: cnt_total , cnt_rain
 logical,parameter :: wgt_2d = .false.  ! true - weight is calculated in 9x7 x-y 2d grid. (CSU order)
                                       ! false- weight is calcualted in 64  x   1d grid. (native GCE grid)

 integer,parameter :: mxlyr_20km = 38  ! GCE level of 20km above sea level 
 real :: phi  !polarization angle [rad]
 real :: a, b ! coefficients 
!
!===========  First, get gain function for all MMF samples. ===========================
!

 if( .not. allocated(bw_array) ) then
    allocate( bw_array(myi_start:myi_end) , &
              tb1d(myj_start:myj_end)     , &
              precip1d(myj_start:myj_end) , &
              wgt(myi_start:myi_end,mxfreq_micro), &
              gf( myi_start:myi_end, myj_start:myj_end, mxfreq_micro ) ,&
              in_ifov( myi_start:myi_end, myj_start:myj_end, mxfreq_micro )   )
 endif

!
! center i grid point is considered to be FOV center. 
! so dx must be adjusted by dx_intercept
!
 dx_intercept = real(myi_end - myi_start + 1) * gridsize / 2.e0  ![km]


 J_LOOP: do j = myj_start, myj_end

    if( trim(sim_case) == 'MMF_GPROF' .and. j > max_sample_gprof ) then
        gf(:,j,:) = undefined
        cycle J_LOOP
    endif

    if( wgt_2d ) then

       !
       ! derive incidence angle
       !
       call derive_incidence_angle( ona_angle_micro*const_degrad, sma_micro, surface(1,j)%lat*const_degrad, &
                                    incidence_ang,  dist_sat_fov )

       !
       ! derive beam width (transform 64 x-grid to 9x7x-y grid) <-- Requirement from C. Kummerow. 
       !
       call beam_width_9x7plane( gridsize, incidence_ang, dist_sat_fov, bw_array )

    else

       !
       ! estimate angle from center pixel to neighbor  surface(i,j)%lat
       !
       do i = myi_start, myi_end
          dx =  real(i) * gridsize - dx_intercept  !distance from center grid [km]
          call beam_width_dx(surface(i,j)%lat , sma_micro , ona_angle_micro, dx, bw_deg ) !derive angle [deg]
          bw_array(i) = bw_deg  !beam width from center grid  [deg]
       enddo

    endif

    !
    ! get antenna gain function  (Gaussian weight based on half-power beamwidth)
    !
    do nf = 1, mxfreq_micro ; do i = myi_start, myi_end
       call gain_function_gaussian( dble(beamwidth(nf)*const_degrad), dble(bw_array(i)*const_degrad), wgt_dble  )
       wgt(i,nf) = real( wgt_dble )  ! single <- double
    enddo ; enddo

    !
    ! normalized gain function (sum == 1)
    !
    do nf = 1, mxfreq_micro ; do i = myi_start, myi_end
       gf(i,j,nf) = wgt(i,nf) / sum( wgt(:,nf) )  !normalized gain functi
    enddo ; enddo
    !
    ! logical (in IFOV) 
    !
    do nf = 1, mxfreq_micro ; do i = myi_start, myi_end
       if( bw_array(i) <= beamwidth(nf)*0.5 ) then
         in_ifov(i,j,nf) = .true.
       else
         in_ifov(i,j,nf) = .false.
       endif
    enddo ; enddo

    ! Note that i=64th grid should not be used in 9x7 sampling. 
    in_ifov(64,:,:) = .false.

    if(j==myj_start .and. masterproc .and. wgt_2d) then
    do nf = 1, mxfreq_micro
       print*,''
       print*,'FOV footprint of Microwave Tb at '//trim(nch_micro(nf))//' in 9x7 2D plane'
       print*,in_ifov(55:63,j,nf)
       print*,in_ifov(46:54,j,nf)
       print*,in_ifov(37:45,j,nf)
       print*,in_ifov(28:36,j,nf)   ! 31,32,33 for CloudSat W
       print*,in_ifov(19:27,j,nf)
       print*,in_ifov(10:18,j,nf)
       print*,in_ifov(1 : 9,j,nf)
    enddo
    endif

 enddo J_LOOP


!
! =============open output file name, and prepare header ==========================
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(micro_sensor)//&
  '.GPROF'//trim(output_suffix)//'.nc'
 print("(a)"),'MSG write_out_gprof: output ->',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )

 is=myi_start ; ie=myi_end ; js=1 ; je=max_sample_gprof ; ks=1

!
! Now MMF GCE is only stored at level of 20km ASL. 
! (actual mxlyr is 44, but upper is just sponge layers).  
!
 ke = mxlyr_20km


 if(ie == is ) then
   im = is
 else
   im = (ie-is+1) / 2
 endif


!
! Define dimensions
!
 call check( nf90_def_dim(ncid, 'sample'  , je-js+1 , j_dimid) )
 call check( nf90_def_dim(ncid, 'mxlyr'   , ke-ks+1 , k_dimid) )
 call check( nf90_def_dim(ncid, 'mxlyr_radar', mxlyr_radar-ks+1 , k_radar_dimid) )

 dimids2d_radar = (/ j_dimid, k_radar_dimid /)
 dimids2d       = (/ j_dimid, k_dimid /)
 dimids1d       = (/ j_dimid /)


 call check( nf90_def_var(ncid, 'Datasetnum', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Dataset number ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_data = id

 call check( nf90_def_var(ncid, 'Latitude', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' latitude (-90:90) ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_lat = id

 call check( nf90_def_var(ncid, 'Longitude', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' longitude (-180:180) ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_lon = id

!
! time header & ID
!
 call check( nf90_def_var(ncid, 'year', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Year) )') )
 call check( nf90_put_att(ncid, id, 'units', '[year]') )
 id_year = id

 call check( nf90_def_var(ncid, 'month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[month]') )
 id_month = id

 call check( nf90_def_var(ncid, 'day_of_month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Days of Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[day]') )
 id_day = id

 call check( nf90_def_var(ncid, 'hour', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (hour) )') )
 call check( nf90_put_att(ncid, id, 'units', '[hour]') )
 id_hour = id

 call check( nf90_def_var(ncid, 'minute', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (minute) )') )
 call check( nf90_put_att(ncid, id, 'units', '[min]') )
 id_min = id
 call check( nf90_def_var(ncid, 'second', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (second) )') )
 call check( nf90_put_att(ncid, id, 'units', '[sec]') )
 id_sec = id

 call check( nf90_def_var(ncid, 'SurfaceRain', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Rain ') )
 call check( nf90_put_att(ncid, id, 'units', '[mm/hr]') )
 id_rain = id

 call check( nf90_def_var(ncid, 'FrozenPrecip', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description','  Frozen Precip ') )
 call check( nf90_put_att(ncid, id, 'units', '[mm/hr]') )
 id_snow = id

 call check( nf90_def_var(ncid, 'Rain_fraction', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description','  Rain fraction (>1mm/hr) ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_rf = id


!
! L1B brightness temperature
!
 ich = 0
 do nf=1,mxfreq_micro !frequency loop
    do np=2,1,-1  !polarization loop
       if(np==1) pol = '(H)'
       if(np==2) pol = '(V)'

       !default
       description = 'Microwave Tb at '//trim(nch_micro(nf))//'Hz'//pol

       call skip_channel(trim(micro_sensor), nf, np, cycle_this, 'Microwave Tb ', description)
       if(cycle_this) cycle

       ich = ich + 1

       write(cc,"(I2.2)") ich
       call check( nf90_def_var(ncid, 'Tb_'//cc , NF90_REAL, dimids1d, id) )
       call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
       call check( nf90_put_att(ncid, id, 'description',trim(description) ) )
       call check( nf90_put_att(ncid, id, 'units', '[K]') )
       id_tb(ich) = id

    enddo
 enddo

!
! De-convoluted L1 brightness temperature in 37Ghz- FOV
!

 ich = 0
 do nf=1,mxfreq_micro !frequency loop
    do np=2,1,-1  !polarization loop
       if(np==1) pol = '(H)'
       if(np==2) pol = '(V)'

       !default
       description = 'Deconvoluted Microwave Tb at '//trim(nch_micro(nf))//'Hz'//pol

       call skip_channel(trim(micro_sensor), nf, np, cycle_this, 'Deconvoluted Microwave Tb ', description)
       if(cycle_this) cycle

       ich = ich + 1


       write(cc,"(I2.2)") ich
!       call check( nf90_def_var(ncid, 'Tb37_'//cc , NF90_REAL, dimids1d, id) )
       call check( nf90_def_var(ncid, 'TbDeconv_'//cc , NF90_REAL, dimids1d, id) )
       call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
       call check( nf90_put_att(ncid, id, 'description',trim(description) ) )
       call check( nf90_put_att(ncid, id, 'units', '[K]') )
       id_tb37(ich) = id

       !toshii I found the very strange bugs. dimids2d_radar is (1,3) in default, but
       !around the end of this loop, it becomes (57, 58), and screw up the dimesion of radar output. 
       !this may be related to the NetCDF library. 
       ! here is brute force correction in next line. 
       dimids2d_radar = (/ j_dimid, k_radar_dimid /)

    enddo
 enddo

 call check( nf90_def_var(ncid, 'SurfWndSpd', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Wind speed ') )
 call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
 id_wind = id

 call check( nf90_def_var(ncid, 'SkinTemp', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Skin Temperature ') )
 call check( nf90_put_att(ncid, id, 'units', '[K]') )
 id_skin = id

 call check( nf90_def_var(ncid, 'T2m', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' 2m Air Temperature ') )
 call check( nf90_put_att(ncid, id, 'units', '[K]') )
 id_t2m = id

 call check( nf90_def_var(ncid, 'SnowDepth', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Snow Depth ') )
 call check( nf90_put_att(ncid, id, 'units', '[m]') )
 id_snowd = id

 call check( nf90_def_var(ncid, 'ELEV', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Elevation ') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_elev = id

 call check( nf90_def_var(ncid, 'TCWV', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Total Column Water Vapor ') )
 call check( nf90_put_att(ncid, id, 'units', '[mm]') )
 id_tcwv = id

 call check( nf90_def_var(ncid, 'U1km', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' 1km above surface U component winds ') )
 call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
 id_u1km = id

 call check( nf90_def_var(ncid, 'V1km', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' 1km above surface V component winds ') )
 call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
 id_v1km = id

 call check( nf90_def_var(ncid, 'RWC', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,mxlyr)') )
 call check( nf90_put_att(ncid, id, 'description',' Rain water content ') )
 call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
 id_rwc = id

 call check( nf90_def_var(ncid, 'LWC', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,mxlyr)') )
 call check( nf90_put_att(ncid, id, 'description',' Cloud liquid water content ') )
 call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
 id_lwc = id

 call check( nf90_def_var(ncid, 'MWC', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,mxlyr)') )
 call check( nf90_put_att(ncid, id, 'description',' Mixed-phase water content ') )
 call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
 id_mwc = id

 call check( nf90_def_var(ncid, 'IWC', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,mxlyr)') )
 call check( nf90_put_att(ncid, id, 'description',' Cloud ice water content ') )
 call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
 id_iwc = id

 call check( nf90_def_var(ncid, 'LH', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,mxlyr)') )
 call check( nf90_put_att(ncid, id, 'description',' Latent Heating  ') )
 call check( nf90_put_att(ncid, id, 'units', '[K/day]') )
 id_lh = id

 call check( nf90_def_var(ncid, 'THICKNESS', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,mxlyr)') )
 call check( nf90_put_att(ncid, id, 'description',' layer thickness ') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_dhgt = id

 call check( nf90_enddef(ncid) )  ! End define mode.

!
! =============open output file name, and prepare header ==========================
!

 if( allocated(var1d) ) deallocate( var1d, var2d, ivar1d, temp2d )
 allocate( var1d(js:je), var2d(js:je,ks:ke), &
           tb(is:ie,js:je), ivar1d(js:je), temp2d(is:ie,js:je) )


  select case(trim(micro_sensor))
  case('TMI_PRE')
    ivar1d(js:je) =  9
  case('AMSR-E')
    ivar1d(js:je) = 10
  case('TMI_POST')
    ivar1d(js:je) = 11
  case('SSMIS')
    ivar1d(js:je) = 12
  case('GMI')
    ivar1d(js:je) = 13
  case('MADRAS')
    ivar1d(js:je) = 14
  case('MHS')
    ivar1d(js:je) = 15
  case('SAPHIR')
    ivar1d(js:je) = 16
  case('ATMS')
    ivar1d(js:je) = 17
  case('MICROMAS')
    ivar1d(js:je) = 18
  case('MWHS_II_LIGHT')
    ivar1d(js:je) = 19
  case default
    call stop_sdsu(myrank,'MSG write_out_gprof: no such micro_sensor type')
      stop
  end select
  call check( nf90_put_var(ncid, id_data, ivar1d(js:je) )  ) !data number


  var1d(js:je) = surface(im,js:je)%lat
  call check( nf90_put_var(ncid, id_lat, var1d(js:je) )  )  !latitude

  do j = js, je
     if( surface(im,j)%lon > 180. ) then
        var1d(j) = surface(im,j)%lon - 360.
     else
        var1d(j) = surface(im,j)%lon
     endif
  enddo
  call check( nf90_put_var(ncid, id_lon, var1d(js:je) )  )  !longitude

  ivar1d(js:je) = int(sdsu_yyyy)
  call check( nf90_put_var(ncid, id_year, ivar1d(js:je) )  )  !year
  ivar1d(js:je) = int(sdsu_mm)
  call check( nf90_put_var(ncid, id_month, ivar1d(js:je) )  ) !month
  ivar1d(js:je) = int(sdsu_dd)
  call check( nf90_put_var(ncid, id_day, ivar1d(js:je) )  )   !day of month
  ivar1d(js:je) = int(sdsu_hh)
  call check( nf90_put_var(ncid, id_hour, ivar1d(js:je) )  )  !hour
  ivar1d(js:je) = int(sdsu_nn)
  call check( nf90_put_var(ncid, id_min, ivar1d(js:je) )  )   !minute
  ivar1d(js:je) = int(sdsu_ss)
  call check( nf90_put_var(ncid, id_sec, ivar1d(js:je) )  )   !second

  call dump_gprof( surface(is:ie,js:je)%rain_rate_liq, gf(is:ie,js:je,nf_L2), id_rain, ncid ) !liquid rain
  call dump_gprof( surface(is:ie,js:je)%rain_rate_ice, gf(is:ie,js:je,nf_L2), id_snow, ncid ) !solid precip

!
! rain fraction 
!
 do j = js, je

    cnt_total = 0 ; cnt_rain = 0
    do i = is, ie
       if(in_ifov(i,j,nf_L2)) then
          cnt_total = cnt_total + 1
          if(surface(i,j)%rain_rate_liq+surface(i,j)%rain_rate_ice >= 1.e0) then
             cnt_rain = cnt_rain + 1
          endif
       endif
    enddo

    if( cnt_total == 0 ) then
       var1d(j) = 0.
    else
       var1d(j) = real(cnt_rain)/real(cnt_total)  !rain fraction. 
    endif
    call check( nf90_put_var(ncid, id_rf, var1d(js:je) )  )
 enddo


!
! L1 brightness temperature 
!
 ich = 0
 do nf=1,mxfreq_micro !frequency loop
    do np=2,1,-1  !polarization loop

       call skip_channel(trim(micro_sensor), nf, np, cycle_this)
       if(cycle_this) cycle

       ! default value
       tb(is:ie,js:je) = tb_out(is:ie,js:je,nf,np)

       if(trim(micro_sensor)=='SSMIS') then
         if(nf==11) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,11,1) + tb_out(is:ie,js:je,16,1) )
         if(nf==12) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,12,1) + tb_out(is:ie,js:je,15,1) )
         if(nf==13) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,13,1) + tb_out(is:ie,js:je,14,1) )
       endif

       if(trim(micro_sensor)=='GMI') then
         if(nf==7) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,7,2) + tb_out(is:ie,js:je,10,2) )
         if(nf==8) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,8,2) + tb_out(is:ie,js:je, 9,2) )
       endif

       if(trim(micro_sensor)=='MHS') then
         if(nf==3) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,3,np) + tb_out(is:ie,js:je, 6,np) )
         if(nf==4) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,4,np) + tb_out(is:ie,js:je, 5,np) )
       endif

       if(trim(micro_sensor)=='SAPHIR') then
         if(nf==1) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,1,np) + tb_out(is:ie,js:je,12,np) )
         if(nf==2) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,2,np) + tb_out(is:ie,js:je,11,np) )
         if(nf==3) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,3,np) + tb_out(is:ie,js:je,10,np) )
         if(nf==4) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,4,np) + tb_out(is:ie,js:je, 9,np) )
         if(nf==5) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,5,np) + tb_out(is:ie,js:je, 8,np) )
         if(nf==6) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,6,np) + tb_out(is:ie,js:je, 7,np) )
       endif

       if(trim(micro_sensor)=='ATMS') then
         if(nf==5) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,5,np) + tb_out(is:ie,js:je,14,np) )
         if(nf==6) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,6,np) + tb_out(is:ie,js:je,13,np) )
         if(nf==7) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,7,np) + tb_out(is:ie,js:je,12,np) )
         if(nf==8) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,8,np) + tb_out(is:ie,js:je,11,np) )
         if(nf==9) tb(is:ie,js:je) = 0.5*( tb_out(is:ie,js:je,9,np) + tb_out(is:ie,js:je,10,np) )
       endif

       if(trim(micro_sensor)=='MICROMAS') then

         if( nf <= 8) phi = 70. * const_degrad   ! based on MICROMAS ATBD
         if( nf >= 9) phi = -20. * const_degrad  ! based on MICROMAS ATBD  
         a = sin(phi)*sin(phi) ; b = cos(phi)*cos(phi)
         !fixed polarization degree 
         tb(is:ie,js:je) = a*tb_out(is:ie,js:je,nf,2) + b*tb_out(is:ie,js:je,nf,1)

       endif



       ich = ich + 1
       call dump_gprof( tb(is:ie,js:je), gf(is:ie,js:je,np)   , id_tb(ich) , ncid ) ! Tb

       ! can be any channel depending on specification of nf_L2
       call dump_gprof( tb(is:ie,js:je), gf(is:ie,js:je,nf_L2), id_tb37(ich), ncid ) ! convoluted Tb 


       !
       ! accumulate statistical composite between Tb and precipitation for validation.
       !
       if( .not. allocated(tb1d)     ) allocate( tb1d    (js:je) )
       if( .not. allocated(precip1d) ) allocate( precip1d(js:je) )
       do j = js , je
          tb1d    (j) = sum(      tb(is:ie,j)          *gf(is:ie,j,np) )
          precip1d(j) = sum( surface(is:ie,j)%rain_rate*gf(is:ie,j,np) )
       enddo
       call tb_precip_sample( 'add', ich, tb1d(js:je), precip1d(js:je) )


    enddo
 enddo

! 1D output
  call dump_gprof( surface(is:ie,js:je)%u10m   , gf(is:ie,js:je,nf_L2), id_wind, ncid )  ! surface wind
  call dump_gprof( surface(is:ie,js:je)%t_skin , gf(is:ie,js:je,nf_L2), id_skin, ncid )  ! skin temp
  call dump_gprof( surface(is:ie,js:je)%dhgt_snow , gf(is:ie,js:je,nf_L2), id_snowd, ncid ) ! snow depth
  call dump_gprof( surface(is:ie,js:je)%h2o_col, gf(is:ie,js:je,nf_L2), id_tcwv, ncid )  ! TCWV
  call dump_gprof( surface(is:ie,js:je)%u1km   , gf(is:ie,js:je,nf_L2), id_u1km, ncid )  ! U-wind 1km AGL
  call dump_gprof( surface(is:ie,js:je)%v1km   , gf(is:ie,js:je,nf_L2), id_v1km, ncid )  ! V-wind 1km AGL

  var1d(js:je) = atmos_stag(im,js:je,0)%hgt   ! surface elevation [km]
  call check( nf90_put_var(ncid, id_elev, var1d(js:je) )  )
! 2D output
  call dump_gprof( q_gen(is:ie,js:je,ks:ke)%rain , gf(is:ie,js:je,nf_L2), id_rwc, ncid ) ! RWC (make sure...)
  call dump_gprof( q_gen(is:ie,js:je,ks:ke)%cloud, gf(is:ie,js:je,nf_L2), id_lwc, ncid ) ! LWC 
  var2d(js:je,ks:ke) = undefined
  call check( nf90_put_var(ncid, id_mwc, var2d(js:je,ks:ke) )  )                       ! MWC (undefined)
  call dump_gprof( q_gen(is:ie,js:je,ks:ke)%ice, gf(is:ie,js:je,nf_L2), id_iwc, ncid ) ! IWC (make sure...)
  var2d(js:je,ks:ke) = undefined
  call check( nf90_put_var(ncid, id_lh, var2d(js:je,ks:ke) )  )                        ! LH (undefined)

  var2d(js:je,ks:ke) = atmos(im,js:je,ks:ke)%dhgt    ! thickness of atmospheric layer [km]
  call check( nf90_put_var(ncid, id_dhgt, var2d(js:je,ks:ke) )  )                      ! THICKNESS


!
! 2m air temp (temporally defined in this way)
!
  do j = js, je ; do i = is,ie
     temp2d(i,j) = 0.5*( surface(i,j)%t_skin + atmos(i,j,1)%t_air )  ! 2m air temperature [k]
  enddo ; enddo
  call dump_gprof( temp2d(is:ie,js:je) , gf(is:ie,js:je,nf_L2), id_t2m, ncid )  ! 2m air temperature


!
! accumulate statisitcal sample for GPROF verification.
!
 call tpw_tskin_precip_sample('add',  surface(is:ie,js:je)%t_skin   , &
                                      surface(is:ie,js:je)%h2o_col  , &
                                      surface(is:ie,js:je)%rain_rate, &
                                      gf     (is:ie,js:je,nf_L2)    , &
                                      surface(is:ie,js:je)%lat      , &
                                      surface(is:ie,js:je)%iland      &
                                        )



!
! close file
!
 call check( nf90_close(ncid) )  ! close nc file


 return
 end subroutine write_out_gprof

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_prps
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!
! History:
! 04/2012  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer :: i,j,k,nf,n
 integer :: is,ie,js,je,im, ns,ne

!netcdf
 integer :: ncid, id, ich
 integer :: i_dimid, j_dimid, k_dimid, n_dimid
 integer, parameter :: ndims2d = 2, ndims1d=1
 integer :: dimids1d(ndims1d), dimids2d(ndims2d), dimids1d_scan(ndims1d)
 integer :: id_data, id_lat, id_lon, id_scan, id_inc !var IDs 
 integer :: id_year, id_month, id_day, id_hour, id_min, id_sec !var IDs (time)
 integer :: id_rain, id_rf, id_snow, id_snowd, id_tb(max_chan), id_tb37(max_chan)
 integer :: id_wind, id_skin, id_t2m, id_tcwv, id_u1km, id_v1km, id_elev
 integer :: id_rwc, id_lwc, id_mwc, id_iwc, id_lh, id_dhgt
 integer :: id_zeku, id_zeka, id_zew1,id_zew2,id_zew3, id_piaw1,id_piaw2,id_piaw3
 character(len=100) :: description
 character :: pol*3, cc*2

 real(sdsu_fps) :: incidence_ang ! incidence angle [rad]
 real(sdsu_fps) :: dist_sat_fov ! distance between satellite and FOV [km]

 real(sdsu_fps),allocatable :: var1d(:) , var2d(:,:), tb(:,:,:,:), tb_vh(:,:,:), &
                               temp2d(:,:), tb2d(:,:), precip2d(:,:)  !temporal array
 integer,allocatable :: ivar1d(:)  !temporal array
 logical :: cycle_this
 integer :: cnt_total , cnt_rain
 logical,parameter :: wgt_2d = .false.  ! true - weight is calculated in 9x7 x-y 2d grid. (CSU order)
                                       ! false- weight is calcualted in 64  x   1d grid. (native GCE grid)

 integer,parameter :: mxlyr_20km = 38  ! GCE level of 20km above sea level 
 real :: phi  ! sensor viewing angle [rad]
 real :: wgt_h, wgt_v ! coefficients for H and V pol. 

 logical,parameter :: output_prps = .false.   !usually true 


!
! =============open output file name, and prepare header ==========================
!

if(output_prps) then

 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(micro_sensor)//&
  '.GPROF'//trim(output_suffix)//'.nc'
 print("(a)"),'MSG write_out_prps: output ->',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )

endif

 is=myi_start ; ie=myi_end ; js=1 ; je=max_sample_gprof 
 ns=min_index_angle ; ne=max_index_angle

 if(ie == is ) then
   im = is
 else
   im = (ie-is+1) / 2
 endif


if(output_prps) then
!
! Define dimensions
!
 call check( nf90_def_dim(ncid, 'sample'  , je-js+1 , j_dimid) )
 call check( nf90_def_dim(ncid, 'scan'    , ne-ns+1 , n_dimid ) )
 
 dimids1d_scan  = (/ n_dimid /)
 dimids2d       = (/ j_dimid, n_dimid /)
 dimids1d       = (/ j_dimid /)


 call check( nf90_def_var(ncid, 'scan_angle', NF90_REAL, dimids1d_scan, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(scan)') )
 call check( nf90_put_att(ncid, id, 'description',' sensor scan angles ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_scan = id

 call check( nf90_def_var(ncid, 'IncidenceAngle', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Incidence Angle ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_inc = id

 call check( nf90_def_var(ncid, 'Datasetnum', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Dataset number ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_data = id

 call check( nf90_def_var(ncid, 'Latitude', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' latitude (-90:90) ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_lat = id

 call check( nf90_def_var(ncid, 'Longitude', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' longitude (-180:180) ') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_lon = id

!
! time header & ID
!
 call check( nf90_def_var(ncid, 'year', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Year) )') )
 call check( nf90_put_att(ncid, id, 'units', '[year]') )
 id_year = id

 call check( nf90_def_var(ncid, 'month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[month]') )
 id_month = id

 call check( nf90_def_var(ncid, 'day_of_month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Days of Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[day]') )
 id_day = id

 call check( nf90_def_var(ncid, 'hour', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (hour) )') )
 call check( nf90_put_att(ncid, id, 'units', '[hour]') )
 id_hour = id

 call check( nf90_def_var(ncid, 'minute', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (minute) )') )
 call check( nf90_put_att(ncid, id, 'units', '[min]') )
 id_min = id
 call check( nf90_def_var(ncid, 'second', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (second) )') )
 call check( nf90_put_att(ncid, id, 'units', '[sec]') )
 id_sec = id


!
! L1B brightness temperature
!
 ich = 0
 do nf=1,mxfreq_micro !frequency loop
    pol = '(VH)'

    !default
    description = 'Microwave Tb at '//trim(nch_micro(nf))//'GHz'//pol

    call skip_channel_prps(trim(micro_sensor), nf, cycle_this, 'Microwave Tb ', description)
    if(cycle_this) cycle

    ich = ich + 1

    write(cc,"(I2.2)") ich
    call check( nf90_def_var(ncid, 'Tb_'//cc , NF90_REAL, dimids2d, id) )
    call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
    call check( nf90_put_att(ncid, id, 'description',trim(description) ) )
    call check( nf90_put_att(ncid, id, 'units', '[K]') )
    id_tb(ich) = id

 enddo



 call check( nf90_def_var(ncid, 'SurfaceRain', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Rain ') )
 call check( nf90_put_att(ncid, id, 'units', '[mm/hr]') )
 id_rain = id

 call check( nf90_def_var(ncid, 'FrozenPrecip', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
 call check( nf90_put_att(ncid, id, 'description','  Frozen Precip ') )
 call check( nf90_put_att(ncid, id, 'units', '[mm/hr]') )
 id_snow = id

 call check( nf90_def_var(ncid, 'SurfWndSpd', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Wind speed ') )
 call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
 id_wind = id

 call check( nf90_def_var(ncid, 'SkinTemp', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Skin Temperature ') )
 call check( nf90_put_att(ncid, id, 'units', '[K]') )
 id_skin = id

 call check( nf90_def_var(ncid, 'T2m', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
 call check( nf90_put_att(ncid, id, 'description',' 2m Air Temperature ') )
 call check( nf90_put_att(ncid, id, 'units', '[K]') )
 id_t2m = id

 call check( nf90_def_var(ncid, 'SnowDepth', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Snow Depth ') )
 call check( nf90_put_att(ncid, id, 'units', '[m]') )
 id_snowd = id

 call check( nf90_def_var(ncid, 'ELEV', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface Elevation ') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_elev = id

 call check( nf90_def_var(ncid, 'TCWV', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(sample,scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Total Column Water Vapor ') )
 call check( nf90_put_att(ncid, id, 'units', '[mm]') )
 id_tcwv = id

 call check( nf90_enddef(ncid) )  ! End define mode.

endif

!
! =============write data  ==========================
!

 if( allocated(var1d) ) deallocate( var1d, var2d, ivar1d, temp2d )
 allocate( var1d(js:je), var2d(js:je,ns:ne), &
           tb(is:ie,js:je,1:2,ns:ne), tb_vh(is:ie,js:je,ns:ne), ivar1d(js:je), temp2d(is:ie,js:je) )


if(output_prps) then
  select case(trim(micro_sensor))
  case('TMI_PRE')
    ivar1d(js:je) =  9
  case('AMSR-E')
    ivar1d(js:je) = 10
  case('TMI_POST')
    ivar1d(js:je) = 11
  case('SSMIS')
    ivar1d(js:je) = 12
  case('GMI')
    ivar1d(js:je) = 13
  case('MADRAS')
    ivar1d(js:je) = 14
  case('MHS')
    ivar1d(js:je) = 15
  case('SAPHIR')
    ivar1d(js:je) = 16
  case('ATMS')
    ivar1d(js:je) = 17
  case('MICROMAS')
    ivar1d(js:je) = 18
  case('MWHS_II_LIGHT')
    ivar1d(js:je) = 19
  case default
    call stop_sdsu(myrank,'MSG write_out_prps: no such micro_sensor type')
      stop
  end select


!
! Writing index date/location/time info
!

  call check( nf90_put_var(ncid, id_scan, scan_angle_prps(ns:ne) )  ) !scan angle

  call check( nf90_put_var(ncid, id_inc, inc_angle_prps(im,js:je,ns:ne) ) ) !incidence angle

  call check( nf90_put_var(ncid, id_data, ivar1d(js:je) )  ) !data number


  var1d(js:je) = surface(im,js:je)%lat
  call check( nf90_put_var(ncid, id_lat, var1d(js:je) )  )  !latitude

  do j = js, je
     if( surface(im,j)%lon > 180. ) then
        var1d(j) = surface(im,j)%lon - 360.
     else
        var1d(j) = surface(im,j)%lon
     endif
  enddo
  call check( nf90_put_var(ncid, id_lon, var1d(js:je) )  )  !longitude

  ivar1d(js:je) = int(sdsu_yyyy)
  call check( nf90_put_var(ncid, id_year, ivar1d(js:je) )  )  !year
  ivar1d(js:je) = int(sdsu_mm)
  call check( nf90_put_var(ncid, id_month, ivar1d(js:je) )  ) !month
  ivar1d(js:je) = int(sdsu_dd)
  call check( nf90_put_var(ncid, id_day, ivar1d(js:je) )  )   !day of month
  ivar1d(js:je) = int(sdsu_hh)
  call check( nf90_put_var(ncid, id_hour, ivar1d(js:je) )  )  !hour
  ivar1d(js:je) = int(sdsu_nn)
  call check( nf90_put_var(ncid, id_min, ivar1d(js:je) )  )   !minute
  ivar1d(js:je) = int(sdsu_ss)
  call check( nf90_put_var(ncid, id_sec, ivar1d(js:je) )  )   !second
endif

!
! L1 brightness temperature
!
 ich = 0
 do nf=1,mxfreq_micro !frequency loop

   call skip_channel_prps(trim(micro_sensor), nf, cycle_this)
   if(cycle_this) cycle

    tb(is:ie,js:je,1:2,ns:ne) = tb_prps_out(is:ie,js:je,nf,1:2,ns:ne)

   ! for sounder channel (plus-minus frequency)
   select case(trim(micro_sensor))
   case('MHS')
     if(nf==3) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,3,1:2,ns:ne) + tb_prps_out(is:ie,js:je, 6,1:2,ns:ne) )
     if(nf==4) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,4,1:2,ns:ne) + tb_prps_out(is:ie,js:je, 5,1:2,ns:ne) )
     !VH polarization mixture for MHS sensor
     do n = ns, ne
        phi = scan_angle_prps(n) * const_degrad  ! scan angle [rad]
        if(nf==1 .or. nf==2 .or. nf==7) then ! standard channel    V                                       H 
           tb_vh(is:ie,js:je,n) = cos(phi)*cos(phi)*tb(is:ie,js:je,2,n) + sin(phi)*sin(phi)*tb(is:ie,js:je,1,n)
        else                       ! for window   channel          H                                       V 
           tb_vh(is:ie,js:je,n) = cos(phi)*cos(phi)*tb(is:ie,js:je,1,n) + sin(phi)*sin(phi)*tb(is:ie,js:je,2,n)    
        endif
     enddo
   case('SAPHIR') 
     if(nf==1) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,1,1:2,ns:ne) + tb_prps_out(is:ie,js:je,12,1:2,ns:ne) )
     if(nf==2) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,2,1:2,ns:ne) + tb_prps_out(is:ie,js:je,11,1:2,ns:ne) )
     if(nf==3) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,3,1:2,ns:ne) + tb_prps_out(is:ie,js:je,10,1:2,ns:ne) )
     if(nf==4) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,4,1:2,ns:ne) + tb_prps_out(is:ie,js:je, 9,1:2,ns:ne) )
     if(nf==5) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,5,1:2,ns:ne) + tb_prps_out(is:ie,js:je, 8,1:2,ns:ne) )
     if(nf==6) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,6,1:2,ns:ne) + tb_prps_out(is:ie,js:je, 7,1:2,ns:ne) )
     do n = ns, ne
        phi = scan_angle_prps(n) * const_degrad  ! scan angle [rad]
        wgt_v = 0. ; wgt_h = 1.   !toshii fix this later..... 
        tb_vh(is:ie,js:je,n) = wgt_v*tb(is:ie,js:je,2,n) + wgt_h*tb(is:ie,js:je,1,n)
     enddo
   case('ATMS')
     if(nf==5) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,5,1:2,ns:ne) + tb_prps_out(is:ie,js:je,14,1:2,ns:ne) )
     if(nf==6) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,6,1:2,ns:ne) + tb_prps_out(is:ie,js:je,13,1:2,ns:ne) )
     if(nf==7) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,7,1:2,ns:ne) + tb_prps_out(is:ie,js:je,12,1:2,ns:ne) )
     if(nf==8) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,8,1:2,ns:ne) + tb_prps_out(is:ie,js:je,11,1:2,ns:ne) )
     if(nf==9) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,9,1:2,ns:ne) + tb_prps_out(is:ie,js:je,10,1:2,ns:ne) )

     !VH polarization mixture for ATMS sensor
     do n = ns, ne
        phi = scan_angle_prps(n) * const_degrad  ! scan angle [rad]
        if(nf >  3) then           ! for window channel            H                                       V 
           tb_vh(is:ie,js:je,n) = cos(phi)*cos(phi)*tb(is:ie,js:je,1,n) + sin(phi)*sin(phi)*tb(is:ie,js:je,2,n)
        else                       ! for standard channel(1,2,3)   V                                       H 
           tb_vh(is:ie,js:je,n) = cos(phi)*cos(phi)*tb(is:ie,js:je,2,n) + sin(phi)*sin(phi)*tb(is:ie,js:je,1,n)
        endif
     enddo


   case('MICROMAS')
    if( nf <= 8) phi =  70. * const_degrad   ! based on MICROMAS ATBD
    if( nf >= 9) phi = -20. * const_degrad   ! based on MICROMAS ATBD
    !fixed polarization degree
    tb_vh(is:ie,js:je,ns:ne) = sin(phi)*sin(phi)*tb(is:ie,js:je,2,ns:ne) + cos(phi)*cos(phi)*tb(is:ie,js:je,1,ns:ne)

   case('MWHS_II_LIGHT')

    if(nf==2 ) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,2,1:2,ns:ne) + tb_prps_out(is:ie,js:je,7,1:2,ns:ne) )
    if(nf==3 ) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,3,1:2,ns:ne) + tb_prps_out(is:ie,js:je,6,1:2,ns:ne) ) 
    if(nf==4 ) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,4,1:2,ns:ne) + tb_prps_out(is:ie,js:je,5,1:2,ns:ne) )
    if(nf==9 ) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,9,1:2,ns:ne) + tb_prps_out(is:ie,js:je,14,1:2,ns:ne) ) 
    if(nf==10) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,10,1:2,ns:ne) + tb_prps_out(is:ie,js:je,13,1:2,ns:ne) )
    if(nf==11) tb(is:ie,js:je,1:2,ns:ne) = 0.5*( tb_prps_out(is:ie,js:je,11,1:2,ns:ne) + tb_prps_out(is:ie,js:je,12,1:2,ns:ne) )

     !VH polarization mixture for MWHS sensor
     do n = ns, ne
        phi = scan_angle_prps(n) * const_degrad  ! scan angle [rad]
        if(nf == 1 .or. nf == 8) then  ! for standard channel      H                                       V 
           tb_vh(is:ie,js:je,n) = cos(phi)*cos(phi)*tb(is:ie,js:je,1,n) + sin(phi)*sin(phi)*tb(is:ie,js:je,2,n)
        else                       ! for window channel            V                                       H 
           tb_vh(is:ie,js:je,n) = cos(phi)*cos(phi)*tb(is:ie,js:je,2,n) + sin(phi)*sin(phi)*tb(is:ie,js:je,1,n)
        endif
     enddo



   case default

   end select


   ich = ich + 1

if(output_prps) then
   call dump_prps( tb_vh(is:ie,js:je,ns:ne), gf_prps(is:ie,js:je,nf,ns:ne), id_tb(ich), ncid ) ! Tb
endif

   !
   ! accumulate statistical composite between Tb and precipitation for validation.
   !
   if( .not. allocated(tb2d)     ) allocate( tb2d    (js:je,ns:ne) )
   if( .not. allocated(precip2d) ) allocate( precip2d(js:je,ns:ne) )
   do n = ns, ne ; do j = js , je
      tb2d    (j,n) = sum(  tb_vh(is:ie,j,n)          *gf_prps(is:ie,j,nf,n) )
      precip2d(j,n) = sum( surface(is:ie,j)%rain_rate*gf_prps(is:ie,j,nf,n) )
   enddo ; enddo

   call tb_precip_prps( 'add', ich, surface(im,js:je)%lat, surface(im,js:je)%iland, &
                        tb2d(js:je,ns:ne), precip2d(js:je,ns:ne) )

 enddo

if(output_prps) then
!
! geophysical parameters
!
 call dump_prps( surface(is:ie,js:je)%rain_rate_liq, gf_prps(is:ie,js:je,nf_L2,ns:ne), id_rain, ncid ) ! liquid rain
 call dump_prps( surface(is:ie,js:je)%rain_rate_ice, gf_prps(is:ie,js:je,nf_L2,ns:ne), id_snow, ncid ) ! solid precip
 call dump_prps( surface(is:ie,js:je)%u10m         , gf_prps(is:ie,js:je,nf_L2,ns:ne), id_wind, ncid ) ! surface wind
 call dump_prps( surface(is:ie,js:je)%t_skin       , gf_prps(is:ie,js:je,nf_L2,ns:ne), id_skin, ncid ) ! skin temp
 call dump_prps( surface(is:ie,js:je)%dhgt_snow    , gf_prps(is:ie,js:je,nf_L2,ns:ne), id_snowd,ncid ) ! snow depth
 call dump_prps( surface(is:ie,js:je)%h2o_col      , gf_prps(is:ie,js:je,nf_L2,ns:ne), id_tcwv, ncid ) ! TCWV
 call dump_prps( atmos_stag(is:ie,js:je,0)%hgt     , gf_prps(is:ie,js:je,nf_L2,ns:ne), id_elev, ncid ) ! elevation

!
! 2m air temp (temporally defined in this way)
!
  do j = js, je ; do i = is,ie
     temp2d(i,j) = 0.5*( surface(i,j)%t_skin + atmos(i,j,1)%t_air )  ! 2m air temperature [k]
  enddo ; enddo
  call dump_prps( temp2d(is:ie,js:je) , gf_prps(is:ie,js:je,nf_L2,ns:ne), id_t2m, ncid )  ! 2m air temperature

!
! accumulate statisitcal sample for GPROF verification.
!
 call tpw_tskin_precip_sample('add',  surface(is:ie,js:je)%t_skin   , &
                                      surface(is:ie,js:je)%h2o_col  , &
                                      surface(is:ie,js:je)%rain_rate, &
                                      gf_prps(is:ie,js:je,nf_L2,0)  , &
                                      surface(is:ie,js:je)%lat      , &
                                      surface(is:ie,js:je)%iland      &
                                        )

!
! close file
!
 call check( nf90_close(ncid) )  ! close nc file
endif


 return
 end subroutine write_out_prps

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine skip_channel(sensor_name, nf, np, skip_this, header, description)
 implicit none
 character(*),intent(in) :: sensor_name
 integer,intent(in) :: nf
 integer,intent(in) :: np
 logical,intent(out) :: skip_this
 character(*),intent(in),optional :: header
 character(*),intent(inout),optional :: description

!
! initialize as false
!
 skip_this = .false.

 select case(trim(micro_sensor))
 case('TMI','TMI_POST','TMI_PRE')
    if(nf==3 .and. np==1) skip_this = .true. !no 24G horizontal channel

 case('GMI')
    if(nf==3 .and. np==1) skip_this = .true. !no 24G horizontal channel
    if(nf>=7 .and. np==1) skip_this = .true. !no H channel in 183GHz
    if(nf==7) description = trim(header)//'at 183.3plus-minus7GHz (V)'
    if(nf==8) description = trim(header)//'at 183.3plus-minus3GHz (V)'
    if(nf>=9) skip_this = .true.

 case('SSMIS')
    if(nf==2 .and. np==1) skip_this = .true.  !no 22.235G horizontal channel (V only)
    if( (nf>=4 .and. nf<=8) .and. np==2 )  skip_this = .true. !no O3-band vertical channel (H only)
    if(nf>=10 .and. np==2) skip_this = .true.  !no V in sounding channles
    if(nf==11) description = trim(header)//'at 183.3plus-minus6.6GHz (H)'
    if(nf==12) description = trim(header)//'at 183.3plus-minus3.0GHz (H)'
    if(nf==13) description = trim(header)//'at 183.3plus-minus1.0GHz (H)'
    if(nf>=14) skip_this = .true.
 case('MADRAS')
    if(nf==2 .and. np==1) skip_this = .true. !no 24G horizontal channel
 case('MHS')
    if(nf==3 .and. np==1) description = trim(header)//'at 183.3plus-minus3.0GHz (H)'
    if(nf==3 .and. np==2) description = trim(header)//'at 183.3plus-minus3.0GHz (V)'
    if(nf==4 .and. np==1) description = trim(header)//'at 183.3plus-minus1.0GHz (H)'
    if(nf==4 .and. np==2) description = trim(header)//'at 183.3plus-minus1.0GHz (V)'
    if(nf==5) skip_this = .true.
    if(nf==6) skip_this = .true.
 case('SAPHIR')
    if(np == 2) skip_this = .true.  !no V chan in sounding
    if(nf >= 7) skip_this = .true.
    if(nf==1 .and. np==1) description = trim(header)//'at 183.3plus-minus11.0GHz (H)'
    if(nf==2 .and. np==1) description = trim(header)//'at 183.3plus-minus6.8GHz (H)'
    if(nf==3 .and. np==1) description = trim(header)//'at 183.3plus-minus4.2GHz (H)'
    if(nf==4 .and. np==1) description = trim(header)//'at 183.3plus-minus2.8GHz (H)'
    if(nf==5 .and. np==1) description = trim(header)//'at 183.3plus-minus1.1GHz (H)'
    if(nf==6 .and. np==1) description = trim(header)//'at 183.3plus-minus0.2GHz (H)'
 case('ATMS')
    if(nf==5 .and. np==1) description = trim(header)//'at 183.3plus-minus7.0GHz (H)'
    if(nf==5 .and. np==2) description = trim(header)//'at 183.3plus-minus7.0GHz (V)'
    if(nf==6 .and. np==1) description = trim(header)//'at 183.3plus-minus4.5GHz (H)'
    if(nf==6 .and. np==2) description = trim(header)//'at 183.3plus-minus4.5GHz (V)'
    if(nf==7 .and. np==1) description = trim(header)//'at 183.3plus-minus3.0GHz (H)'
    if(nf==7 .and. np==2) description = trim(header)//'at 183.3plus-minus3.0GHz (V)'
    if(nf==8 .and. np==1) description = trim(header)//'at 183.3plus-minus1.8GHz (H)'
    if(nf==8 .and. np==2) description = trim(header)//'at 183.3plus-minus1.8GHz (V)'
    if(nf==9 .and. np==1) description = trim(header)//'at 183.3plus-minus1.0GHz (H)'
    if(nf==9 .and. np==2) description = trim(header)//'at 183.3plus-minus1.0GHz (V)'
    if(nf>=10) skip_this = .true.

 case('MICROMAS')
    if( np == 1 ) skip_this = .true.  !skip H channel in sounding 
    if(nf==1 .and. np==2) description = trim(header)//'at 91.655GHz (VH)'
    if(nf==2 .and. np==2) description = trim(header)//'at 114.5GHz (VH)'
    if(nf==3 .and. np==2) description = trim(header)//'at 115.95GHz (VH)'
    if(nf==4 .and. np==2) description = trim(header)//'at 116.65GHz (VH)'
    if(nf==5 .and. np==2) description = trim(header)//'at 117.25GHz (VH)'
    if(nf==6 .and. np==2) description = trim(header)//'at 117.80GHz (VH)'
    if(nf==7 .and. np==2) description = trim(header)//'at 118.24GHz (VH)'
    if(nf==8 .and. np==2) description = trim(header)//'at 118.58GHz (VH)'
    if(nf==9 .and. np==2) description = trim(header)//'at 184.41GHz (VH)'
    if(nf==10 .and. np==2) description = trim(header)//'at 186.51GHz (VH)'
    if(nf==11 .and. np==2) description = trim(header)//'at 190.31GHz (VH)'
    if(nf==12 .and. np==2) description = trim(header)//'at 204.8GHz (VH)'
 case default
    !do nothing
 end select

 return
 end subroutine skip_channel

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine skip_channel_prps(sensor_name, nf, skip_this, header, description)
 implicit none
 character(*),intent(in) :: sensor_name
 integer,intent(in) :: nf
 logical,intent(out) :: skip_this
 character(*),intent(in),optional :: header
 character(*),intent(inout),optional :: description

!
! initialize as false
!
 skip_this = .false.

 select case(trim(micro_sensor))
 case('MHS')
    if(nf==3) description = trim(header)//' at 183.3plus-minus3.0GHz (VH)'
    if(nf==4) description = trim(header)//' at 183.3plus-minus1.0GHz (VH)'
    if(nf==5) skip_this = .true.
    if(nf==6) skip_this = .true.
 case('SAPHIR')
    if(nf >= 7) skip_this = .true.
    if(nf==1 ) description = trim(header)//' at 183.3plus-minus11.0GHz (VH)'
    if(nf==2 ) description = trim(header)//' at 183.3plus-minus6.8GHz (VH)'
    if(nf==3 ) description = trim(header)//' at 183.3plus-minus4.2GHz (VH)'
    if(nf==4 ) description = trim(header)//' at 183.3plus-minus2.8GHz (VH)'
    if(nf==5 ) description = trim(header)//' at 183.3plus-minus1.1GHz (VH)'
    if(nf==6 ) description = trim(header)//' at 183.3plus-minus0.2GHz (VH)'
 case('ATMS' )
    if(nf==5 ) description = trim(header)//' at 183.3plus-minus7.0GHz (VH)'
    if(nf==6 ) description = trim(header)//' at 183.3plus-minus4.5GHz (VH)'
    if(nf==7 ) description = trim(header)//' at 183.3plus-minus3.0GHz (VH)'
    if(nf==8 ) description = trim(header)//' at 183.3plus-minus1.8GHz (VH)'
    if(nf==9 ) description = trim(header)//' at 183.3plus-minus1.0GHz (VH)'
    if(nf>=10) skip_this = .true.
 case('MICROMAS')
    if(nf==1 ) description = trim(header)//' at 91.655GHz (VH)'
    if(nf==2 ) description = trim(header)//' at 114.5GHz (VH)'
    if(nf==3 ) description = trim(header)//' at 115.95GHz (VH)'
    if(nf==4 ) description = trim(header)//' at 116.65GHz (VH)'
    if(nf==5 ) description = trim(header)//' at 117.25GHz (VH)'
    if(nf==6 ) description = trim(header)//' at 117.80GHz (VH)'
    if(nf==7 ) description = trim(header)//' at 118.24GHz (VH)'
    if(nf==8 ) description = trim(header)//' at 118.58GHz (VH)'
    if(nf==9 ) description = trim(header)//' at 184.41GHz (VH)'
    if(nf==10) description = trim(header)//' at 186.51GHz (VH)'
    if(nf==11) description = trim(header)//' at 190.31GHz (VH)'
    if(nf==12) description = trim(header)//' at 204.8GHz (VH)'

 case('MWHS_II_LIGHT')
    skip_this = .true.
    if(nf==1 ) then ; description = trim(header)//' at 89.0GHz (VH)' ; skip_this = .false. ; endif
    if(nf==2 ) then ; description = trim(header)//' at 118.75plus-minus5.0GHz (VH)' ; skip_this = .false. ; endif
    if(nf==3 ) then ; description = trim(header)//' at 118.75plus-minus3.0GHz (VH)' ; skip_this = .false. ; endif
    if(nf==4 ) then ; description = trim(header)//' at 118.75plus-minus1.1GHz (VH)' ; skip_this = .false. ; endif
    if(nf==8 ) then ; description = trim(header)//' at 150.0GHz (VH)'               ; skip_this = .false. ; endif
    if(nf==9 ) then ; description = trim(header)//' at 183.3plus-minus7.0GHz (VH)'  ; skip_this = .false. ; endif
    if(nf==10) then ; description = trim(header)//' at 183.3plus-minus3.0GHz (VH)'  ; skip_this = .false. ; endif
    if(nf==11) then ; description = trim(header)//' at 183.3plus-minus1.0GHz (VH)'  ; skip_this = .false. ; endif


 case default
    !do nothing
 end select

 return
 end subroutine skip_channel_prps

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump_gprof_2d( var2d, wgt2d, varid, ncid)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!
! History:
! 04/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),dimension(:,:), intent(in)  :: var2d  ! 
 real(sdsu_fps),dimension(:,:), intent(in)  :: wgt2d  ! 
 integer,intent(in)    :: varid   ! var index
 integer,intent(inout) :: ncid    ! netcdf file id 
 integer :: i,j
 integer :: is,ie,js,je
 integer :: bnd(2)  !upper memory bound for 3D array
 real(sdsu_fps),allocatable :: out1d(:) 

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(var2d)
 is = 1 ; ie = bnd(1) ; js = 1 ; je = bnd(2)
 allocate( out1d(js:je) )

 do j = js , je
    out1d(j) = sum( var2d(is:ie,j)*wgt2d(is:ie,j) )        
 enddo

 call check( nf90_put_var(ncid, varid, out1d) )


 deallocate( out1d )
 return
 end subroutine dump_gprof_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump_gprof_3d( var3d, wgt2d, varid, ncid)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!
! History:
! 04/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),dimension(:,:,:), intent(in)  :: var3d  ! 
 real(sdsu_fps),dimension(:,:), intent(in)  :: wgt2d  ! 
 integer,intent(in)    :: varid   ! var index
 integer,intent(inout) :: ncid    ! netcdf file id 
 integer :: i,j,k
 integer :: is,ie,js,je, ks, ke
 integer :: bnd(3)  !upper memory bound for 3D array
 real(sdsu_fps),allocatable :: out2d(:,:)

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(var3d)
 is = 1 ; ie = bnd(1) ; js = 1 ; je = bnd(2) ; ks = 1 ; ke = bnd(3)
 allocate( out2d(js:je,ks:ke) )

 do j = js , je
    do k = ks, ke
       out2d(j,k) = sum( var3d(is:ie,j,k)*wgt2d(is:ie,j) )
    enddo
 enddo

 call check( nf90_put_var(ncid, varid, out2d) )


 deallocate( out2d )
 return
 end subroutine dump_gprof_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump_prps_2d( var2d, wgt3d, varid, ncid)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!
! History:
! 04/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),dimension(:,:), intent(in)  :: var2d  ! 
 real(sdsu_fps),dimension(:,:,:), intent(in)  :: wgt3d  ! 
 integer,intent(in)    :: varid   ! var index
 integer,intent(inout) :: ncid    ! netcdf file id 
 integer :: i,j,k,n
 integer :: is,ie,js,je, ns, ne
 integer :: bnd(3)  !upper memory bound for 3D array
 real(sdsu_fps),allocatable :: out2d(:,:)

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(wgt3d)
 is = 1 ; ie = bnd(1) ; js = 1 ; je = bnd(2) ; ns = 1 ; ne = bnd(3)
 allocate( out2d(js:je,ns:ne) )

 do n = ns, ne
    do j = js , je

       out2d(j,n) = sum( var2d(is:ie,j)*wgt3d(is:ie,j,n) )

    enddo
 enddo

 call check( nf90_put_var(ncid, varid, out2d) )

 deallocate( out2d )
 return
 end subroutine dump_prps_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump_prps_3d( var3d, wgt3d, varid, ncid)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!
! History:
! 04/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),dimension(:,:,:), intent(in)  :: var3d  ! 
 real(sdsu_fps),dimension(:,:,:), intent(in)  :: wgt3d  ! 
 integer,intent(in)    :: varid   ! var index
 integer,intent(inout) :: ncid    ! netcdf file id 
 integer :: i,j,k,n
 integer :: is,ie,js,je, ns, ne
 integer :: bnd(3)  !upper memory bound for 3D array
 real(sdsu_fps),allocatable :: out2d(:,:)

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(wgt3d)
 is = 1 ; ie = bnd(1) ; js = 1 ; je = bnd(2) ; ns = 1 ; ne = bnd(3)
 allocate( out2d(js:je,ns:ne) )

 do n = ns, ne
    do j = js , je

       out2d(j,n) = sum( var3d(is:ie,j,n)*wgt3d(is:ie,j,n) )

    enddo
 enddo

 call check( nf90_put_var(ncid, varid, out2d) )

 deallocate( out2d )
 return
 end subroutine dump_prps_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine open_oldfile(io_unit,io_file)
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  This subroutine open existing file more carefully. It actually try 10 times before it getting failure. 
! 
! History:
! 09/2012  Toshi Matsui@NASA GSFC : Initial. 
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: io_unit
 character(*),intent(in) :: io_file  !io file name and directory


 integer :: open_stat 
 integer :: ntries  !number of tries to open the file
 character(len=200) :: error_msg

 ntries = 0   ! number of tries to open the file

 forever: do

    open (unit=io_unit, file=trim(io_file), status='old', iostat=open_stat)

    if (open_stat /= 0) then
       ntries = ntries + 1
       if (ntries > 10) then  !if more than 10 sec
           error_msg = 'MSG open_oldfile: Cannot open file '//trim(io_file)
           call stop_sdsu(myrank,trim(error_msg))
       endif
       call sleep(1)   ! Fail openiing a file,  so wait 1 second before retrying
    else 
       exit forever   !success opening a file
    endif

 enddo forever

 return
 end subroutine open_oldfile

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine init_gen(param,proc)
 implicit none
 type( particle_gen ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%cloud    =0.0
   param%rain     =0.0
   param%ice      =0.0
   param%snow     =0.0
   param%graupel  =0.0
   param%hail     =0.0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_gen

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_gen_r8(param,proc)
 implicit none
 type( particle_gen_r8 ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%cloud    =0.d0
   param%rain     =0.d0
   param%ice      =0.d0
   param%snow     =0.d0
   param%graupel  =0.d0
   param%hail     =0.d0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_gen_r8

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine init_gen_cmplx(param,proc)
 implicit none
 type( particle_gen_cmplx ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%cloud    = cmplx(0.,0.)
   param%rain     = cmplx(0.,0.)
   param%ice      = cmplx(0.,0.)
   param%snow     = cmplx(0.,0.)
   param%graupel  = cmplx(0.,0.)
   param%hail     = cmplx(0.,0.)
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_gen_cmplx

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_rams(param,proc)
 implicit none
 type( particle_rams ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%cloud1   =0.0
   param%cloud2   =0.0
   param%rain     =0.0
   param%ice1     =0.0
   param%ice2     =0.0
   param%snow     =0.0
   param%graupel  =0.0
   param%hail     =0.0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_rams

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_rams_r8(param,proc)
 implicit none
 type( particle_rams_r8 ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc 
 
 select case(trim(proc))
 case('zero')
   param%cloud1   =0.0
   param%cloud2   =0.0
   param%rain     =0.0
   param%ice1     =0.0
   param%ice2     =0.0
   param%snow     =0.0
   param%graupel  =0.0
   param%hail     =0.0
 case('nan')
   
 case('inf')

 end select

 return
 end subroutine init_rams_r8

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_rams_cmplx(param,proc)
 implicit none
 type( particle_rams_cmplx ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%cloud1   =cmplx(0.,0.)
   param%cloud2   =cmplx(0.,0.)
   param%rain     =cmplx(0.,0.)
   param%ice1     =cmplx(0.,0.)
   param%ice2     =cmplx(0.,0.)
   param%snow     =cmplx(0.,0.)
   param%graupel  =cmplx(0.,0.)
   param%hail     =cmplx(0.,0.)
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_rams_cmplx

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_sbm(param,proc)
 implicit none
 type( particle_sbm ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%liq    =0.0
   param%ice_col=0.0
   param%ice_pla=0.0
   param%ice_den=0.0
   param%snow   =0.0
   param%graupel=0.0
   param%hail   =0.0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_sbm_r8(param,proc)
 implicit none
 type( particle_sbm_r8 ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%liq    =0.d0
   param%ice_col=0.d0
   param%ice_pla=0.d0
   param%ice_den=0.d0
   param%snow   =0.d0
   param%graupel=0.d0
   param%hail  =0.d0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_sbm_r8

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_sbm_cmplx(param,proc)
 implicit none
 type( particle_sbm_cmplx ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%liq    = cmplx(0.,0)
   param%ice_col= cmplx(0.,0)
   param%ice_pla= cmplx(0.,0)
   param%ice_den= cmplx(0.,0)
   param%snow   = cmplx(0.,0)
   param%graupel= cmplx(0.,0)
   param%hail   = cmplx(0.,0)
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_sbm_cmplx

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_gocart(param,proc)
 implicit none
 type( particle_gocart ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%so4   =0.0
   param%blc   =0.0
   param%ocn   =0.0
   param%och   =0.0
   param%ssa   =0.0
   param%ssc   =0.0
   param%du1   =0.0
   param%du2   =0.0
   param%du3   =0.0
   param%du4   =0.0
   param%du5   =0.0
   param%du6   =0.0
   param%du7   =0.0
   param%du8   =0.0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_gocart

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_gocart_r8(param,proc)
 implicit none
 type( particle_gocart_r8 ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%so4   =0.0
   param%blc   =0.0
   param%ocn   =0.0
   param%och   =0.0
   param%ssa   =0.0
   param%ssc   =0.0
   param%du1   =0.0
   param%du2   =0.0
   param%du3   =0.0
   param%du4   =0.0
   param%du5   =0.0
   param%du6   =0.0
   param%du7   =0.0
   param%du8   =0.0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_gocart_r8

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_ucd(param,proc)
 implicit none
 type( particle_ucd ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%du1   =0.0
   param%du2   =0.0
   param%du3   =0.0
   param%du4   =0.0
   param%du5   =0.0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_ucd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_ucd_r8(param,proc)
 implicit none
 type( particle_ucd_r8 ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%du1   =0.0
   param%du2   =0.0
   param%du3   =0.0
   param%du4   =0.0
   param%du5   =0.0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_ucd_r8


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_hid(param,proc)
 implicit none
 type( particle_hid ) ,intent(inout) :: param
 character(len=*),intent(in) :: proc

 select case(trim(proc))
 case('zero')
   param%dz  =0.0
   param%rn  =0.0
   param%cr  =0.0
   param%ldg =0.0
   param%hdg =0.0
   param%vi  =0.0
   param%ds  =0.0
   param%ws  =0.0
   param%ha  =0.0
   param%bd  =0.0
 case('nan')

 case('inf')

 end select

 return
 end subroutine init_hid

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine vt_rime_melting(ibin, nbin, rf, fl, ms, vt_dr, vt_sn, vt_ha, rho_sn, vt_out)

  implicit none

  integer, intent(in ) :: ibin
  integer, intent(in ) :: nbin
  real(sdsu_fps), intent(in ) :: rf ! mass riming fraction
  real(sdsu_fps),intent(in )  :: fl ! mass melting fraction
  real(sdsu_fps),intent(in )  :: ms(nbin)    !mass [g]
  real(sdsu_fps),intent(in )  :: vt_dr(nbin) !vt of liquid [cm/s]
  real(sdsu_fps),intent(in )  :: vt_sn(nbin) !vt of snow [cm/s]
  real(sdsu_fps),intent(in )  :: vt_ha(nbin) !vt of hail [cm/s]
  real(sdsu_fps),intent(in )  :: rho_sn(nbin)! [g/cm3]
  real(sdsu_fps),intent(out)  :: vt_out      !vt [cms]

  real :: ms_sn_dry, rho_sn_int, rho_sn_dry, rho_sn_blk
  real :: a, b, c
  integer :: kr


  if (fl > 0.d0) then
    ms_sn_dry = ms(ibin) * (1.d0 - fl)
    if (ms_sn_dry <= ms(1)) then
      rho_sn_int = rho_sn(1)
    else
      do kr = 1, nbin-1
        if (ms_sn_dry >= ms(kr) .and. ms_sn_dry <= ms(kr+1)) then
          exit
        endif
      enddo
      rho_sn_int = rho_sn(kr) + (rho_sn(kr+1) - rho_sn(kr)) / &
                  (ms(kr+1) - ms(kr)) * (ms_sn_dry - ms(kr))
    endif
    rho_sn_dry = rho_sn_int * 0.9d0 / (0.9d0*(1.d0 - rf) + rho_sn_int * rf)
    a = rho_sn_dry * 1.0d0
    b = 1.d0 * (1.d0 - fl)
    c = rho_sn_dry * fl
    rho_sn_blk = a / (b+c)
    b = 0.9d0 - rho_sn(ibin)
    c = rho_sn_blk 
    if (b == 0.d0) then
      a = 0.d0
    else
      a = (c - rho_sn(ibin)) / b
    endif
    vt_out = vt_sn(ibin) * (1.d0 - a) + a * (vt_ha(ibin) * (1.d0 - fl) + vt_dr(ibin) * fl)
  else
    if (rf <= 0.d0) then
      rho_sn_blk = rho_sn(ibin)
    elseif (rf >= 1.d0) then
      rho_sn_blk = 0.9d0
    else
      ms_sn_dry = ms(ibin) * (1.d0 - rf)
      if (ms_sn_dry <= ms(1)) then
        rho_sn_int = rho_sn(1)
      else
        do kr = 1, nbin-1
          if (ms_sn_dry >= ms(kr) .and. ms_sn_dry <= ms(kr+1)) then
            exit
          endif
        enddo

        rho_sn_int = rho_sn(kr) + (rho_sn(kr+1) - rho_sn(kr)) / &
                    (ms(kr+1) - ms(kr)) * (ms_sn_dry - ms(kr))
      endif
      rho_sn_blk = rho_sn_int * 0.9d0 / (0.9d0*(1.d0 - rf) + rho_sn_int * rf)
    endif
    b = 0.9d0 - rho_sn(ibin)
    c = rho_sn_blk 
    if (b == 0.d0) then
      a = 0.d0
    else
      a = (c - rho_sn(ibin)) / b
    endif
    vt_out = vt_sn(ibin) * (1.d0 - a) + a * vt_ha(ibin)
  endif


  return
  end subroutine vt_rime_melting

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine  remove_cirrus_deck( qcol, q1d )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  This subroutine remove cirrus deck from MMF simuations.
! 
! History:
! 11/2013  Toshi Matsui@NASA GSFC : Initial. 
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 type ( particle_gen ), intent(in) :: qcol ! collumn hydrometeor mass [kg/m2]
 type ( particle_gen ), dimension(mxlyr),intent(inout) :: q1d
 integer :: k

!
! if rain present and cirrus deck present, but without snow+grapuel+hail amount
!
 if( qcol%ice+qcol%snow < 0.1e0 ) then !could be thin cirrus 
    q1d(:)%ice     = 0. !zero out
    q1d(:)%snow    = 0. !zero out
    q1d(:)%graupel = 0. !zero out
    q1d(:)%hail    = 0. !zero out 
 endif

 return
 end subroutine remove_cirrus_deck 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine tb_precip_sample(proc,ich,tb,precip)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!
! History:
! 07/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 character*(*) :: proc ! process
 integer,intent(in),optional :: ich
 real(sdsu_fps),intent(in),optional ::  tb(:)
 real(sdsu_fps),intent(in),optional ::  precip(:)

 integer :: j,ip,ic,io      ! loop indice

 integer,parameter :: maxbin_precip = 200
 integer,parameter :: maxbin_chan   = 12
 integer :: bnd(1)  !upper memory bound for 3D array

 type matrix_stat
   real(sdsu_fps) :: tb_sum
   real(sdsu_fps) :: cnt_sum
   real(sdsu_fps) :: tb_ave
 end type matrix_stat
 type( matrix_stat ), dimension(1:maxbin_precip,1:maxbin_chan),save :: matrix  !sum of PDF 


 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------
!print*,'MSG tb_precip_sample: init'

   matrix(:,:)%tb_sum  = 0.
   matrix(:,:)%cnt_sum = 0.
   matrix(:,:)%tb_ave  = 0.

 case('add') !-------------------------------------------------------------------------------------
!print*,'MSG tb_precip_sample: add'

 if(present(ich) .and. present(tb) .and. present(precip) ) then

    bnd = UBOUND(tb)

    do j = 1, bnd(1)
       ip = nint( precip(j)*10. )
       !print*,ich,j,precip(j),ip,tb(j)
       if( ip >= 1 .and. ip <= maxbin_precip ) then

           if( tb(j) > 50.e0 .and. tb(j) < 400.e0 ) then
               matrix(ip,ich)%tb_sum  = matrix(ip,ich)%tb_sum  + tb(j)
               matrix(ip,ich)%cnt_sum = matrix(ip,ich)%cnt_sum +   1.e0
           endif

       endif
    enddo

 endif

 case('finish') !---------------------------------------------------------------------------------------

!print*,'MSG tb_precip_sample: finish'
#if MPI > 0 
  call mpi_sdsu_reduce( 'SUM', matrix%tb_sum  )
  call mpi_sdsu_reduce( 'SUM', matrix%cnt_sum )

  MPI_MASTERPROC: if(masterproc) then
#endif

!print*,'cnt',matrix(:,:)%cnt_sum

 do ip = 1, maxbin_precip ; do ic = 1, maxbin_chan
    if( matrix(ip,ic)%cnt_sum > 0.) then
      matrix(ip,ic)%tb_ave = matrix(ip,ic)%tb_sum / matrix(ip,ic)%cnt_sum
    else
      matrix(ip,ic)%tb_ave = undefined
    endif
 enddo ; enddo

!print*,'tb_ave',matrix(:,:)%tb_ave


!
! Output in binary format
!
  sdsu_io_file =trim(sdsu_dir_output)//'TB_PRECIP_STAT'//trim(output_suffix)//'.bin'
  print("(a)"),' ->',trim(sdsu_io_file)

  open( unit= io, file = sdsu_io_file, access='direct', status = 'replace', recl=maxbin_precip*4)
  do ic = 1, maxbin_chan
     write(io,rec=ic) (matrix(ip,ic)%tb_ave ,ip=1, maxbin_precip)
  enddo
  close (io)


!
! write out grads control file
!
    sdsu_io_file =trim(sdsu_dir_output)//'TB_PRECIP_STAT'//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = 'TB_PRECIP_STAT'//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  Tb Precip diagram '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',maxbin_precip ,' linear 1. 1'
    write(io,*) 'ydef                  1   linear 1. 1'
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'

if( trim(micro_sensor) == 'MICROMAS' ) then
    write(io,*) 'vars  12'
    write(io,*) 'tb1    0  0  Tb at 91.0GHzHz(VH) [K]'
    write(io,*) 'tb2    0  0  Tb at 114.5GHzHz(VH) [K]'
    write(io,*) 'tb3    0  0  Tb at 115.95GHzHz(VH) [K]'
    write(io,*) 'tb4    0  0  Tb at 116.65GHzHz(VH) [K]'
    write(io,*) 'tb5    0  0  Tb at 117.25GHz(VH)[K]'
    write(io,*) 'tb6    0  0  Tb at 117.80GHz(VH)[K]'
    write(io,*) 'tb7    0  0  Tb at 118.24GHz(VH) [K]'
    write(io,*) 'tb8    0  0  Tb at 118.58GHz(VH) [K]'
    write(io,*) 'tb9    0  0  Tb at 184.41GHz(VH) [K]'
    write(io,*) 'tb10    0  0  Tb at 186.51GHz(VH) [K]'
    write(io,*) 'tb11    0  0  Tb at 190.31GHz(VH) [K]'
    write(io,*) 'tb12    0  0  Tb at 204.8GHz(VH) [K]'
else
    write(io,*) 'vars  12'
    write(io,*) 'tb1    0  0  Tb at 89.0GHzHz(V) [K]'
    write(io,*) 'tb2    0  0  Tb at 89.0GHzHz(H) [K]'
    write(io,*) 'tb3    0  0  Tb at 157.0GHzHz(V) [K]'
    write(io,*) 'tb4    0  0  Tb at 157.0GHzHz(H) [K]'
    write(io,*) 'tb5    0  0  Tb at 183.3plus-minus3.0GHz (V)[K]'
    write(io,*) 'tb6    0  0  Tb at 183.3plus-minus3.0GHz (H)[K]'
    write(io,*) 'tb7    0  0  Tbat 183.3plus-minus1.0GHz (V) [K]'
    write(io,*) 'tb8    0  0  Tbat 183.3plus-minus1.0GHz (H) [K]'
    write(io,*) 'tb9    0  0  Tb at 190.31GHzHz(V) [K]'
    write(io,*) 'tb10    0  0  Tb at 190.31GHzHz(H) [K]'
    write(io,*) 'tb11    0  0  Tb for channel 11 [K]'
    write(io,*) 'tb12    0  0  Tb for channel 12 [K]'
endif
    write(io,*) 'endvars '
    close(io)

 print*,'MSG tb_precip_sample: finished writing output data'

#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 end select

 end subroutine tb_precip_sample

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine tb_precip_prps(proc, ich, lat1d, iland, tb2d, precip2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!
! History:
! 07/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 character*(*) :: proc ! process
 integer,intent(in),optional :: ich  !channel ID
 real(sdsu_fps),intent(in),optional ::  lat1d    (:)                 ! latitude [deg]
 integer       ,intent(in),optional ::  iland    (:)                 ! 1-land,  2-water
 real(sdsu_fps),intent(in),optional ::  tb2d    (:,min_index_angle:) ! brightness temperature (sample,scan) [K]
 real(sdsu_fps),intent(in),optional ::  precip2d(:,min_index_angle:) ! surface precipitation (sample,scan) [mm/hr]

 integer :: n, j,ip,ic,io,nf      ! loop indice
 logical :: cycle_this

 integer,parameter :: sbin = -41 , ebin = 60 ! precipitation bin
 integer :: il, js, je, ns ,ne !scan bin

 real(sdsu_fps) :: precip_bin(sbin:ebin)
 integer :: izone  !0,1,2 
 character(len=5) ::  ltyp  !surface type
 character(len=5) ::  ztyp  !zone type


 type matrix_stat
   real(sdsu_fps) :: tb_sum
   real(sdsu_fps) :: cnt_sum
   real(sdsu_fps) :: tb_ave
 end type matrix_stat
 type( matrix_stat ),allocatable, save :: &
   matrix(:,:,:,:,:)  !sum of PDF (0:2, 0:2, sbin:ebin, 1:max_chan, ns:ne )
                    !           (surface, zone, precip, channel, view angle)
!netcdf
 integer :: ncid, id
 integer :: ip_dimid, n_dimid
 integer, parameter :: ndims2d = 2, ndims1d=1
 integer :: dimids1d_scan(ndims1d), dimids2d(ndims2d), dimids1d_precip(ndims1d)
 integer :: id_scan, id_precip, id_tb(max_chan) !var IDs 
 character(len=100) :: description
 character :: pol*3, cc*2


 select case(trim(proc))

 case('add') !-------------------------------------------------------------------------------------
!print*,'MSG tb_precip_prps: add'

 ns=min_index_angle ; ne= max_index_angle  !scan bin
 js=1               ; je=max_sample_gprof

 if( .not. allocated(matrix) ) then
    allocate( matrix( 0:2, 0:2, sbin:ebin, 1:max_chan, ns:ne ) )
    matrix%tb_sum  = 0.
    matrix%cnt_sum = 0.
    matrix%tb_ave  = 0.
 endif


 if(present(ich) .and. present(tb2d) .and. present(precip2d) ) then

    do n = ns,ne ; do j = js,je
     izone = -1
     if( lat1d(j) > -60.  .and. lat1d(j) < -30. ) izone = 0  !Southern Midlat
     if( lat1d(j) > -30.  .and. lat1d(j) <  30. ) izone = 1  !Tropics
     if( lat1d(j) >  30.  .and. lat1d(j) <  60. ) izone = 2  !Northern Midlat
     if( izone == -1 ) cycle
       !
       ! find precip bin
       !
       if( precip2d(j,n) > 0. ) then
          ip = max( sbin, min( ebin, nint((log10(precip2d(j,n)*100.)-2.0)*20.) ) )
       else
          ip = sbin
       endif

       !
       ! add Tb if it is in the reasonable range. 
       !
       if( tb2d(j,n) > 50.e0 .and. tb2d(j,n) < 400.e0 ) then 

           ! total (land+ocean)
           matrix(0,izone, ip,ich,n)%tb_sum  = matrix(0,izone, ip,ich,n)%tb_sum  + tb2d(j,n)
           matrix(0,izone, ip,ich,n)%cnt_sum = matrix(0,izone, ip,ich,n)%cnt_sum +   1.e0

           ! land-ocean sepeartely 
           il = iland(j)  !1-land 2-ocean
           matrix(il,izone, ip,ich,n)%tb_sum  = matrix(il,izone, ip,ich,n)%tb_sum  + tb2d(j,n)
           matrix(il,izone, ip,ich,n)%cnt_sum = matrix(il,izone, ip,ich,n)%cnt_sum +   1.e0
           
       endif

    enddo ; enddo

 endif

 case('finish') !---------------------------------------------------------------------------------------


! if at this point, matrix is not allocated, no files are processed. So return 
! if( .not. allocated(matrix ) ) return  

! if( .not. allocated(matrix) ) then
!    allocate( matrix( 0:2, 0:2, sbin:ebin, 1:max_chan, ns:ne ) )
!    matrix%tb_sum  = 0.
!    matrix%cnt_sum = 0.
!    matrix%tb_ave  = 0.
! endif



 ns=min_index_angle ; ne= max_index_angle  !scan bin
 
!print*,'MSG tb_precip_prps: finish'
#if MPI > 0 
  call mpi_sdsu_reduce( 'SUM', matrix(:,:,:,:,:)%tb_sum  )
  call mpi_sdsu_reduce( 'SUM', matrix(:,:,:,:,:)%cnt_sum )
#endif


MPI_MASTERPROC: if(masterproc) then

 do il = 0,2 ; do izone = 0, 2; do n = ns, ne ; do ip = sbin, ebin ; do ic = 1, max_chan
    if( matrix(il,izone,ip,ic,n)%cnt_sum > 0.) then
      matrix(il,izone, ip,ic,n)%tb_ave = matrix(il,izone, ip,ic,n)%tb_sum / matrix(il,izone, ip,ic,n)%cnt_sum
    else
      matrix(il,izone, ip,ic,n)%tb_ave = undefined
    endif
 enddo ; enddo ; enddo ; enddo ; enddo

 do ip = sbin, ebin
    precip_bin(ip) = ( 10.**(real(ip)/20.+2) ) / 100. !precip rate of bin [mm/hr]
 enddo


!
! Output in netcdf format
!
 surface_loop: do il = 0,2  !0-total, 1-land, 2-water


    select case(il)
    case(0) ; ltyp = 'TOTAL'
    case(1) ; ltyp = 'LAND '
    case(2) ; ltyp = 'WATER'
    end select

 zone_loop: do izone = 0, 2 !0-south midlatidue, 1-tropics, 2-north midlatitude

    select case(izone)
    case(0) ; ztyp = 'SMID'
    case(1) ; ztyp = 'TROP'
    case(2) ; ztyp = 'NMID'
    end select

 
    sdsu_io_file = trim(sdsu_dir_output)//'TB_PRECIP_STAT_'//trim(micro_sensor)//'_'&
                   //trim(ltyp)//'_'//trim(ztyp)//trim(output_suffix)//'.nc'
    print("(a)"),'Write ->',trim(sdsu_io_file)

    call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )

   !
   ! Define dimensions
   !
   call check( nf90_def_dim(ncid, 'precip_bin', ebin-sbin+1   , ip_dimid) )
   call check( nf90_def_dim(ncid, 'scan_bin'  , ne-ns+1       , n_dimid ) )

   dimids2d          = (/ ip_dimid, n_dimid /)
   dimids1d_precip   = (/ ip_dimid /)
   dimids1d_scan     = (/ n_dimid /)


   call check( nf90_def_var(ncid, 'scan_angle', NF90_REAL, dimids1d_scan, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(scan_bin)') )
   call check( nf90_put_att(ncid, id, 'description',' sensor scan angles ') )
   call check( nf90_put_att(ncid, id, 'units', '[deg]') )
   id_scan = id

   call check( nf90_def_var(ncid, 'precip', NF90_REAL, dimids1d_precip, id) )
   call check( nf90_put_att(ncid, id, 'MemoryOrder','(precip_bin)') )
   call check( nf90_put_att(ncid, id, 'description',' surface precipitation rate ') )
   call check( nf90_put_att(ncid, id, 'units', '[mm/hr]') )
   id_precip = id

   !
   ! L1B mean brightness temperature
   !
   ic = 0
   do nf=1,mxfreq_micro !frequency loop
      pol = '(VH)'

      !default
      description = 'Microwave Tb at '//trim(nch_micro(nf))//'GHz'//pol

      call skip_channel_prps(trim(micro_sensor), nf, cycle_this, 'Microwave Tb ', description)
      if(cycle_this) cycle

      ic = ic + 1

      write(cc,"(I2.2)") ic
      call check( nf90_def_var(ncid, 'Tb_ave_'//cc , NF90_REAL, dimids2d, id) )
      call check( nf90_put_att(ncid, id, 'MemoryOrder','(precip_bin,scan_bin)') )
      call check( nf90_put_att(ncid, id, 'description',trim(description) ) )
      call check( nf90_put_att(ncid, id, 'units', '[K]') )
      id_tb(ic) = id

   enddo

   call check( nf90_enddef(ncid) )  ! End define mode


  !
  ! Writing index date/location/time info
  !
   call check( nf90_put_var(ncid, id_precip, precip_bin(sbin:ebin)  ) ) ! precip bin (toshii)

   call check( nf90_put_var(ncid, id_scan  , scan_angle_prps(ns:ne) ) ) ! scan angle

  !
  ! L1B mean brightness temperature
  !
  ic = 0
  do nf=1,mxfreq_micro !frequency loop

     call skip_channel_prps(trim(micro_sensor), nf, cycle_this, 'Microwave Tb ', description)
     if(cycle_this) cycle

     ic = ic + 1
     call check( nf90_put_var(ncid, id_tb(ic) , matrix(il,izone, :,ic,:)%tb_ave  )  )
  
  enddo

  !
  ! close file
  !
  call check( nf90_close(ncid) )  ! close nc file

 enddo zone_loop

 enddo surface_loop

 endif MPI_MASTERPROC


 end select  !-------------------------------------------------------------------------------------


 return
 end subroutine tb_precip_prps

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine tpw_tskin_precip_sample(proc, tskin,tpw, precip, wgt2d, lat, iland)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: This suborutine construct statistics output for verifying GPROF database.  
!
! History:
! 09/2014  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 character*(*) :: proc ! process
 real(sdsu_fps),dimension(:,:),intent(in),optional :: tskin  !surface skin temperature[K]
 real(sdsu_fps),dimension(:,:),intent(in),optional :: tpw    !total preciptable water [mm]
 real(sdsu_fps),dimension(:,:),intent(in),optional :: precip !surface rain rate [mm/hr]
 real(sdsu_fps),dimension(:,:),intent(in),optional :: wgt2d  !anntenna weight [-] 
 real(sdsu_fps),dimension(:,:),intent(in),optional :: lat    !latitude [deg] 
 integer,       dimension(:,:),intent(in),optional :: iland  !land mask 1:land, 2:ocean [-] 

 integer :: i,j,n      ! loop indice
 integer :: is,ie,js,je
 integer,parameter :: io=988
 integer :: bnd(2)  !upper memory bound for 3D array
 real(sdsu_fps) :: precip0d, tskin0d, tpw0d

 integer,parameter :: minbin_tskin = 180, maxbin_tskin = 350
 integer,parameter :: minbin_tpw   =   0, maxbin_tpw   = 100

 integer :: i_tskin, j_tpw !bin id
 type matrix_stat
   real(sdsu_fps) :: precip_sum
   real(sdsu_fps) :: cnt_sum
   real(sdsu_fps) :: precip_ave
   real(sdsu_fps) :: precip_gmi !GMI-derive precipirate
 end type matrix_stat
 type( matrix_stat ), dimension(minbin_tskin:maxbin_tskin,minbin_tpw:maxbin_tpw), save :: &
     matrix_o   ,&!matrix for ocean
     matrix_l   ,&!matrix for land
     matrix_a     !matrix for land and ocean


 real :: table(172,102) !table value from GMI_MATRIX.txt data


 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------
   ! initialize matrix
   matrix_o(:,:)%precip_sum  = 0.
   matrix_o(:,:)%cnt_sum     = 0.
   matrix_o(:,:)%precip_ave  = 0.

   matrix_l(:,:)%precip_sum  = 0.
   matrix_l(:,:)%cnt_sum     = 0.
   matrix_l(:,:)%precip_ave  = 0.

   matrix_a(:,:)%precip_sum  = 0.
   matrix_a(:,:)%cnt_sum     = 0.
   matrix_a(:,:)%precip_ave  = 0.


   !
   ! derived observed one
   !
   sdsu_io_file = trim(sdsu_dir_data)//'GMI_MATRIX.txt'
   open(io,file=trim(sdsu_io_file))
   do n = 1, 102
      read(io,FMT=*) table(1:172,n)
   enddo
   close(io)

   !
   ! derive GMI precipitation for comparison purpose
   ! 
   do n = 102, 2, -1
      j = 102 - n
      matrix_a(minbin_tskin:maxbin_tskin, j  )%precip_gmi = table(2:172, n)
   enddo

 case('add') !-------------------------------------------------------------------------------------


 if(present(tskin) .and. present(tpw) .and. present(precip) .and. present(wgt2d) ) then

    !
    ! find upper bound from assumped 3D array
    !
    bnd = UBOUND(tskin)
    is = 1 ; ie = bnd(1) ; js = 1 ; je = bnd(2)

    J_LOOP: do j = js , je

       !
       ! gaussian antenna weighted parameters
       !
       precip0d = sum( precip(is:ie,j)*wgt2d(is:ie,j) )  ! [mm/hr]
       tskin0d  = sum( tskin (is:ie,j)*wgt2d(is:ie,j) )  ! [K]
       tpw0d    = sum( tpw   (is:ie,j)*wgt2d(is:ie,j) )  ! [mm]

       !
       ! derive indices for i,j of matrix
       !
       i_tskin = NINT( tskin0d )
       j_tpw   = NINT( tpw0d   )
       if( i_tskin < minbin_tskin .or. i_tskin > maxbin_tskin ) cycle J_LOOP
       if( j_tpw   < minbin_tpw   .or. j_tpw   > maxbin_tpw   ) cycle J_LOOP


       !
       ! tropical accumulated precip and sampling count. 
       !
       if( lat(1,j) >= -37.e0 .and. lat(1,j) <= 37.e0 ) then

           !
           ! global accumulated precip and sampling count. 
           !
           ! for both ocean and land
              matrix_a(i_tskin,j_tpw)%precip_sum  = matrix_o(i_tskin,j_tpw)%precip_sum  + precip0d  !sum of precip
              matrix_a(i_tskin,j_tpw)%cnt_sum     = matrix_o(i_tskin,j_tpw)%cnt_sum     +   1.e0    !sum of count
           if   ( iland(1,j) == 2 ) then !ocean
              matrix_o(i_tskin,j_tpw)%precip_sum  = matrix_o(i_tskin,j_tpw)%precip_sum  + precip0d  !sum of precip
              matrix_o(i_tskin,j_tpw)%cnt_sum     = matrix_o(i_tskin,j_tpw)%cnt_sum     +   1.e0    !sum of count
           elseif( iland(1,j) == 1) then !land
              matrix_l(i_tskin,j_tpw)%precip_sum  = matrix_l(i_tskin,j_tpw)%precip_sum  + precip0d  !sum of precip
              matrix_l(i_tskin,j_tpw)%cnt_sum     = matrix_l(i_tskin,j_tpw)%cnt_sum     +   1.e0    !sum of count
           endif

       endif

    enddo J_LOOP

 endif


 case('finish') !---------------------------------------------------------------------------------------

!print*,'MSG precip_precip_sample: finish'
#if MPI > 0 
  call mpi_sdsu_reduce( 'SUM', matrix_o%precip_sum  )
  call mpi_sdsu_reduce( 'SUM', matrix_o%cnt_sum )
  call mpi_sdsu_reduce( 'SUM', matrix_l%precip_sum  )
  call mpi_sdsu_reduce( 'SUM', matrix_l%cnt_sum )
  call mpi_sdsu_reduce( 'SUM', matrix_a%precip_sum  )
  call mpi_sdsu_reduce( 'SUM', matrix_a%cnt_sum )

  MPI_MASTERPROC: if(masterproc) then
#endif

 !
 ! derive average precipitation per bin
 !
 do i = minbin_tskin,maxbin_tskin ; do j = minbin_tpw, maxbin_tpw

    ! land case
    if( matrix_l(i,j)%cnt_sum > 0.) then
      matrix_l(i,j)%precip_ave = matrix_l(i,j)%precip_sum / matrix_l(i,j)%cnt_sum
    else
      matrix_l(i,j)%precip_ave = undefined
    endif

    ! ocean case
    if( matrix_o(i,j)%cnt_sum > 0.) then
      matrix_o(i,j)%precip_ave = matrix_o(i,j)%precip_sum / matrix_o(i,j)%cnt_sum
    else
      matrix_o(i,j)%precip_ave = undefined
    endif

    ! ocean case
    if( matrix_a(i,j)%cnt_sum > 0.) then
      matrix_a(i,j)%precip_ave = matrix_a(i,j)%precip_sum / matrix_a(i,j)%cnt_sum
    else
      matrix_a(i,j)%precip_ave = undefined
    endif



 enddo ; enddo

!
! Output in binary format
!

  sdsu_io_file =trim(sdsu_dir_output)//'TPW_TSKIN_PRECIP'//trim(output_suffix)//'.bin'
  print("(a)"),' ->',trim(sdsu_io_file)

  open( unit= io, file = sdsu_io_file, access='direct', status = 'replace', &
        recl=(maxbin_tskin-minbin_tskin+1)*(maxbin_tpw-minbin_tpw+1)*4 )
  write(io,rec=1) (( matrix_o(i,j)%precip_ave ,i=minbin_tskin,maxbin_tskin), j=minbin_tpw, maxbin_tpw)
  write(io,rec=2) (( matrix_l(i,j)%precip_ave ,i=minbin_tskin,maxbin_tskin), j=minbin_tpw, maxbin_tpw)
  write(io,rec=3) (( matrix_a(i,j)%precip_ave ,i=minbin_tskin,maxbin_tskin), j=minbin_tpw, maxbin_tpw)
  write(io,rec=4) (( matrix_a(i,j)%precip_gmi ,i=minbin_tskin,maxbin_tskin), j=minbin_tpw, maxbin_tpw)
  close (io)

!
! write out grads control file
!
    sdsu_io_file =trim(sdsu_dir_output)//'TPW_TSKIN_PRECIP'//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = 'TPW_TSKIN_PRECIP'//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  Tropical Precip diagram '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',(maxbin_tskin-minbin_tskin+1),' linear ',minbin_tskin,' 1'
    write(io,*) 'ydef   ',(maxbin_tpw-minbin_tpw+1    ),' linear ',minbin_tpw  ,' 1'
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars  4'
    write(io,*) 'oprecip   0  0  mean ocean precipitation  [mm/hr]'
    write(io,*) 'lprecip   0  0  mean land precipitation   [mm/hr]'
    write(io,*) 'aprecip   0  0  mean land-ocean precipitation   [mm/hr]'
    write(io,*) 'gmiprecip 0  0  GMI mean precipitation    [mm/hr]'
    write(io,*) 'endvars '
    close(io)

    print*,'MSG tpw_tskin_precip_sample: finished writing output data'

#if MPI > 0 
  endif MPI_MASTERPROC
#endif
 end select

 end subroutine tpw_tskin_precip_sample

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_out_cc_mask
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  This subroutine create cloudsat-CALIPSO-defined cloud fractoin. 
! 
! History:
! 07/2015  Toshi Matsui@NASA GSFC : Initial. 
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer i,j,k,l,nw
 integer :: irec
 integer :: kmax
 integer,parameter :: io = 106
 character(len=2) :: ww
 real,allocatable :: cloud_mask(:,:), cloud_fraction(:,:,:)

!
! fileter case situation
!
 check_radar: if( radar .and. trim(radar_sensor) == 'CPR' ) then

 check_lidar: if( lidar .and. trim(lidar_sensor) == 'CALIOP' ) then

 check_sim: if( trim(sim_case) == 'MMF' ) then

 !
 ! vertical profile
 ! 
 if( inst_profile_lidar ) then  !instrumental vertical levels
    kmax = mxlyr_lidar
 else  ! model native vertical levels
    kmax = mxlyr
 endif

!
! allocate array
!
 if( .not. allocated(cloud_mask) ) then
   allocate( cloud_mask(mxgridy,3) , & 
             cloud_fraction(mxgridy,kmax,3) )
                     
 endif  

 do k=1,kmax ; do j=1,mxgridy 

    do i=1,mxgridx

       cloud_mask(i,1:3) = 0.  !no cloud
       if( radar_out_conv(i,j,k,1)%Zm /= undefined ) then
         if( 10.*log10( radar_out_conv(i,j,k,1)%Zm ) >= -30. ) then
               cloud_mask(i,1) = 1.  !consistento COSP
               cloud_mask(i,2) = 1.  !consistento COSP
         endif
       endif

       ! total attenuating backscattering [km-1 str-1]
       if( att_B(i,j,k,1) /= undefined .or. att_clear_B(i,j,k,1) /= undefined ) then
          if( att_B(i,j,k,1) >= att_clear_B(i,j,k,1)*5. ) then
               cloud_mask(i,1) = 1.  !
               cloud_mask(i,3) = 1.  !consistento COSP CALIOP
          endif
       endif
    enddo 

    ! calculate cloud fraction 
    do l = 1, 3
       cloud_fraction(j,k,l) = sum(cloud_mask(1:mxgridx,l)) / real(mxgridx)
    enddo

 enddo ; enddo 

!
! Write out binary data 
!
  sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.CC_MASK'//&
                trim(output_suffix)//'.bin'
  open(unit= io, file = sdsu_io_file, access='direct', &
       status='replace',recl=kmax*1*mxgridy*4)
  irec = 1
  do l = 1, 3
    write(io,rec=irec) ((cloud_fraction(j,k,l),j=1,mxgridy),k=1,kmax)
   irec = irec + 1
  enddo
  close(io)

!
!  write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.CC_MASK'//&
                  trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.CC_MASK'//trim(output_suffix)//'.bin'
    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(radar_sensor)//' reflectivity '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   1  linear ',surface(1,1)%lon,' 1'
    write(io,*) ' ydef             91  linear -90.    2.0'
    write(io,*) 'zdef   ',kmax  ,' linear 0. ',range_radar
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars 3'
    write(io,*) 'ccmask' ,kmax,    ' 0  Cloud Mask [0~1]'
    write(io,*) 'cloudmask' ,kmax, ' 0  Cloud Mask CloudSat only [0~1]'
    write(io,*) 'calmask' ,kmax,   ' 0  Cloud Mask CALIPSO only [0~1]'
    write(io,*) 'endvars '
    close(io)
 endif
 
 endif check_sim
 endif check_lidar
 endif check_radar

 return
 end subroutine write_out_cc_mask

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_cf_all( min_since )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  This subroutine call routine to write out climate-format text output.
! 
! History:
! 05/2015  Toshi Matsui@NASA GSFC : Initial. 
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: min_since


 call write_cf_file('TEMP', min_since)
 call write_cf_file('PRESS', min_since)
 call write_cf_file('RH', min_since)
 call write_cf_file('HGT', min_since)
 call write_cf_file('QC', min_since)
 call write_cf_file('QR', min_since)
 call write_cf_file('QI', min_since)
 call write_cf_file('QS', min_since)
 call write_cf_file('QG', min_since)
 call write_cf_file('QH', min_since)

 return
 end subroutine write_cf_all

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_cf_file(typ, min_since)
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  This subroutine write out climate-format text output.
! 
! History:
! 05/2015  Toshi Matsui@NASA GSFC : Initial. 
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: typ
 integer,intent(in) :: min_since  ![min]

 integer :: i,j,k,irec
 integer,parameter :: io = 201 , io_met = 202
 integer :: is,ie,js,je,ks,ke
 character(len=200) :: hadoop_name

 character :: &
    hh*2 ,&
    mm*2 ,&
    dd*2 ,&
    ss*2 ,&
    mon*3,&
    yyyy*4


 write(hh,"(I2.2)") INT(sdsu_hh)
 write(mm,"(I2.2)") INT(sdsu_nn)
 write(dd,"(I2.2)") INT(sdsu_dd)
 write(ss,"(I2.2)") INT(sdsu_ss)

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end


!
! hadoop CF name
!
 hadoop_name = sdsu_inp_name(1:efile_len-6)//'.'//mm//'.'//ss//'.CF.'


#if MPI == 2
 if(masterproc) then
#endif

   sdsu_io_file = trim(sdsu_dir_output)//trim(hadoop_name)//trim(typ)//trim(output_suffix)//'.txt'
   print("(a)"),'MSG write_cf_file: write climate-format file in text format',&
              ' ->',trim(sdsu_io_file)
   open(unit= io, file = sdsu_io_file )

!
! open metafile
!
    sdsu_io_file = trim(sdsu_dir_output)//trim(hadoop_name)//trim(typ)//trim(output_suffix)//'.met'
    open(unit= io_met, file = sdsu_io_file )
    write(io_met,"(a)") '=== Ascii Climate Format File of NU-WRF output ==='
    write(io_met,"(a)") 'Order: id, time[min], level,latitude[deg],longitude[deg],geophysical parameters'
    write(io_met,"(a)") ''
    write(io_met,*)  'Current Time: ',sdsu_inp_name(efile_len-18:efile_len)
    write(io_met,*)  min_since,'[min] since ',trim(cf_init_time)
    write(io_met,"(a)") ''

#if MPI == 2
 endif
#endif


!
! dump 2D or 3D fields here
!
 select case(trim(typ))

 case('TEMP')
    call dump_cf( min_since, atmos(is:ie,js:je,ks:ke)%t_air , io )
    if(masterproc) write(io_met,*) 'Parameter: Air temperature [Kelvin]'

 case('PRESS')
    call dump_cf( min_since, atmos(is:ie,js:je,ks:ke)%press , io )
    if(masterproc) write(io_met,*) 'Parameter: Air pressure [hPa]'

 case('RH')
    call dump_cf( min_since, atmos(is:ie,js:je,ks:ke)%rh , io )
    if(masterproc) write(io_met,*) 'Parameter: Relative humidity [%]'

 case('HGT')
    call dump_cf( min_since, atmos(is:ie,js:je,ks:ke)%hgt*1000.0 , io )
    if(masterproc) write(io_met,*) 'Parameter: Layer height [m]'

 case('QC')
    call dump_cf( min_since, q_gen(is:ie,js:je,ks:ke)%cloud*1.e-3/atmos(is:ie,js:je,ks:ke)%rho_dair , io )
    if(masterproc) write(io_met,*) 'Parameter: Cloud mixing ratio  [g/g]'

 case('QR')
    call dump_cf( min_since, q_gen(is:ie,js:je,ks:ke)%rain*1.e-3/atmos(is:ie,js:je,ks:ke)%rho_dair , io )
    if(masterproc) write(io_met,*) 'Parameter: Rain mixing ratio  [g/g]'

 case('QI')
    call dump_cf( min_since, q_gen(is:ie,js:je,ks:ke)%ice*1.e-3/atmos(is:ie,js:je,ks:ke)%rho_dair , io )
    if(masterproc) write(io_met,*) 'Parameter: Ice mixing ratio  [g/g]'

 case('QS')
    call dump_cf( min_since, q_gen(is:ie,js:je,ks:ke)%snow*1.e-3/atmos(is:ie,js:je,ks:ke)%rho_dair , io )
    if(masterproc) write(io_met,*) 'Parameter: Snow mixing ratio  [g/g]'

 case('QG')
    call dump_cf( min_since, q_gen(is:ie,js:je,ks:ke)%graupel*1.e-3/atmos(is:ie,js:je,ks:ke)%rho_dair , io )
    if(masterproc) write(io_met,*) 'Parameter: Graupel mixing ratio  [g/g]'

 case('QH')
    call dump_cf( min_since, q_gen(is:ie,js:je,ks:ke)%hail*1.e-3/atmos(is:ie,js:je,ks:ke)%rho_dair , io )
    if(masterproc) write(io_met,*) 'Parameter: Hail mixing ratio  [g/g]'

 case default
   stop 'MSG write_cf_file: no such typ'
 end select

#if MPI == 2
 if(masterproc) then
#endif
   close(io)
   close(io_met)
#if MPI == 2
 endif
#endif 

 return
 end subroutine write_cf_file

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine find_elevation(lat_in, lon_in, elev_km)
 implicit none

 real(sdsu_fps),intent(in) :: lat_in   !latitude [deg] 
 real(sdsu_fps),intent(in) :: lon_in   !longitude [deg] 
 real(sdsu_fps),intent(out) :: elev_km !elevation [km]

 integer :: i,j
 integer :: i_find, j_find
 real(sdsu_fps) :: lat_stag1, lat_stag2, lon_stag1, lon_stag2

!
! default output, if locaiton is not found...
!
 elev_km = undefined
 
 j_loop: do j = 2, mxgridy-1 
    i_loop: do i = 2, mxgridx-1

    lat_stag1 = 0.5*( lat_domain(i,j-1) + lat_domain(i,j  ) ) 
    lat_stag2 = 0.5*( lat_domain(i,j  ) + lat_domain(i,j+1) )  

    lon_stag1 = 0.5*( lon_domain(i-1,j) + lon_domain(i,j  ) )
    lon_stag2 = 0.5*( lon_domain(i,j  ) + lon_domain(i+1,j) )

!    print*,i,j,lat_stag1,lat_stag2,lon_stag1,lon_stag2
    if( lat_in > lat_stag1 .and. lat_in <= lat_stag2 ) then
       if( lon_in > lon_stag1 .and. lon_in <= lon_stag2 ) then

           ! find 
           elev_km = elev_domain(i, j) * 0.001 ! output elevation
           if(masterproc) print*,i,j,'MSG: find_elevation radar instrument is at ',elev_km,'km height', elev_domain(i, j)
           exit j_loop

       endif
    endif
    enddo i_loop
 enddo j_loop


 return
 end subroutine find_elevation

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine scan_angle_bin
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Set up viewing angle of scan bin for PRPS database. 
! 
! History:
! 09/2019  Toshi Matsui@NASA GSFC : Initial. 
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer :: i,j,nf,n
 real(sdsu_fps) :: incidence_ang,  dist_sat_fov


 if(masterproc) print*,'MSG scan_angle_bin: set up PRPS scan bin'

!
! make max and minum scan index
! 
 min_index_angle = INT( min_scan_angle )
 max_index_angle = INT( max_scan_angle )

!
! allocate 
!
 if( .not. allocated( scan_angle_prps ) ) allocate( scan_angle_prps( min_index_angle:max_index_angle ) )

!
! set up scan angle for PRPS database
!
 do n = min_index_angle, max_index_angle

    scan_angle_prps(n) =  real( n )  ! off-nadir sensor viewing angle (scan angle) [deg]

 enddo 

 return
 end subroutine scan_angle_bin

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine prep_scan_array_prps
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Set up viewing angle of scan bin for PRPS database. 
! 
! History:
! 09/2019  Toshi Matsui@NASA GSFC : Initial. 
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer :: i,j,nf,n
 real(sdsu_fps) :: incidence !incidence angle [rad]

 real(sdsu_fps) :: bw_deg !beam width [deg]
 real(sdsu_fps),allocatable :: bw_array(:,:) ! beam width [deg]
 real(sdsu_fps),allocatable :: wgt(:,:,:)    ! weight  [-]
 real(sdsu_fps) :: dx, dx_intercept
 real(sdsu_fpd) :: wgt_dble
 real(sdsu_fps) :: dummy

!
! calculate sensor incidence angle and viewing angles. 
!
 if( .not. allocated(bw_array) ) allocate( bw_array(myi_start:myi_end, min_index_angle:max_index_angle ) )
 if( .not. allocated(wgt) ) allocate( wgt(myi_start:myi_end,mxfreq_micro, min_index_angle:max_index_angle ) )

!         
! center i grid point is considered to be FOV center. 
! so dx must be adjusted by dx_intercept
!        
 dx_intercept = real(myi_end - myi_start + 1) * gridsize / 2.e0  ![km]

             
 j_loop: do j = myj_start, myj_end

    if( j > max_sample_gprof ) then
        gf_prps(:,j,:,:) = undefined  
        cycle j_loop
    endif
    
    !
    ! estimate angle from center pixel to neighbor  surface(i,j)%lat
    !
     do i = myi_start, myi_end  ; do n =  min_index_angle, max_index_angle
        dx =  real(i) * gridsize - dx_intercept  !distance from center grid [km]
        call beam_width_dx(surface(i,j)%lat , sma_micro , scan_angle_prps(n), dx, bw_deg ) !derive angle [deg]
        bw_array(i,n) = bw_deg  !beam width from center grid  [deg]
     enddo ; enddo


    !
    ! get antenna gain function  (Gaussian weight based on half-power beamwidth)
    !
    do nf = 1, mxfreq_micro ; do i = myi_start, myi_end ; do n =  min_index_angle, max_index_angle

       if( bw_array(i,n) > beamwidth(nf)*0.5 ) then !outside of half-power beam width 
          wgt_dble = 0.
       else
          call gain_function_gaussian( dble(beamwidth(nf)*const_degrad), dble(bw_array(i,n)*const_degrad), wgt_dble  )
       endif
       wgt(i,nf,n) = real( wgt_dble )  ! single <- double

    enddo ; enddo ; enddo

    !
    ! normalized gain function (sum == 1)
    !
    do nf = 1, mxfreq_micro ; do i = myi_start, myi_end ; do n =  min_index_angle, max_index_angle
       gf_prps(i,j,nf,n) = wgt(i,nf,n) / sum( wgt(:,nf,n) )  !normalized gain function [-]
    enddo ; enddo ; enddo

    !
    ! incidence angle
    !
    do i = myi_start, myi_end ; do n =  min_index_angle, max_index_angle

       call derive_incidence_angle( abs(scan_angle_prps(n))*const_degrad, sma_micro, &
                                    surface(i,j)%lat*const_degrad, incidence,  dummy )

       inc_angle_prps(i,j,n) = incidence / const_degrad   ! incidence angle for PRPS [deg]
    enddo ; enddo

 enddo j_loop

 return
 end subroutine prep_scan_array_prps

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine column_category(qr_col, qgh_col, wmax_col, qtot, hgt, press, hgt_lcl, iclass)
 implicit none
 real,intent(in) :: qr_col        !column integrated qr [kg/m2]
 real,intent(in) :: qgh_col       !column integrated qg+qh [kg/m2]
 real,intent(in) :: wmax_col      !column integrated maximum vertical velocity [m/s]

 real,intent(in) :: qtot(mxlyr)   !qc+qi [g/m3]
 real,intent(in) :: hgt(mxlyr)    !atmospheric height [km]
 real,intent(in) :: press(mxlyr)  !atmospheric pressure [mb]

 real,intent(in) :: hgt_lcl
 integer,intent(out) :: iclass

 integer :: bnd(1)    ! upper memory bound for 2D array

 integer :: cnt
 integer :: dk, k, k_ctop, k_lcl
 real :: hgt_ctop     ! cloud top height [km]
 real :: press_ctop   !cloud top pressure [km]
 real :: press_lcl    !cloud-base pressure [mb]
 real,parameter :: press_680 = 680.  !seperation for low and mid cloud
 real,parameter :: press_440 = 440.  !seperation for mid and high cloud
 integer :: icld(mxlyr)
 integer :: icld_sum

!
! Seperate column into several class
! 0-undefined, 1-low, 2-mid, 3-high, 4-convective, 5-deep stratiform
! using ice+cloud condensate + cloud-base height + height information
!
 iclass = 0   !undefined

!
! check if cloudy column or not
!
 cnt = 0 
 do k = 1, mxlyr
    if( qtot(k) > q_min_condensate ) cnt = cnt + 1
 enddo
 

 if( cnt == 0 ) then
    !clear sky (no cloud)
     iclass = 0 ; return
 endif

!
! cloudy mask
!
 do k = 1, mxlyr
    if( qtot(k) > q_min_condensate ) then
       icld(k) = 1
    else
       icld(k) = 0
    endif
 enddo
 icld_sum = sum( icld(:) )

!
! find cloud top level
!
 k_ctop = 1  ; press_ctop = press(1)
 ctop_loop: do k = mxlyr, 2, -1   !from top to bottom
    if( qtot(k) > q_min_condensate .and. qtot(k-1) > q_min_condensate ) then  ! at least two layer of cloud
        k_ctop = k  
        hgt_ctop = hgt(k)  !cloud top height (km)
        press_ctop = press(k) !clout top pressure [mb]
        exit ctop_loop
    endif
 enddo ctop_loop

 k_lcl = 1 ; press_lcl = press(1) 
 lcl_loop: do k = mxlyr, 2, -1
    if( hgt(k-1) < hgt_lcl .and. hgt_lcl <= hgt(k) ) then
       k_lcl = k
       press_lcl = press(k)
       exit lcl_loop
    endif
 enddo lcl_loop

 dk = k_ctop - k_lcl  !number of pixels between cloud top and base

!
! asign cloud category
!
 if( press_ctop > press_680 ) then
     iclass = 1  !shallow cloud

 elseif( press_680 > press_ctop .and. press_ctop > press_440  ) then
     iclass = 2  !mid cloud

!     if( real(icld_sum) >=  real(dk)*0.8 ) then  !90% pixel is filled down to cloud base
       if( qr_col > 10. ) iclass = 4  !convective, with intense rain
       if( wmax_col > 3.    ) iclass = 4  !convective, with maximum updraft 
       if( qgh_col > 0.5 ) iclass = 4  !convective, with maximum updraft  
!     endif

 elseif( press_ctop < press_440 ) then
     iclass = 3  !heigh cloud

     if( qr_col > 0.1 ) iclass = 5  !deep stratiform with minimum rain

!     if( real(icld_sum) >=  real(dk)*0.8 ) then  !90% pixel is filled down to cloud base
       if( qr_col > 10. ) iclass = 4  !convective, with intense rain
       if( wmax_col > 3.    ) iclass = 4  !convective, with maximum updraft 
       if( qgh_col > 0.5 ) iclass = 4  !convective, with presence of hail and graupel
!     endif


 endif


 return
 end subroutine column_category

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_simulator

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

