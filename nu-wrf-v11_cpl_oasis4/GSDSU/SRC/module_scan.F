#include <define_CPP.h>
#undef FULL_L2_RADAR
!
! if you deinfe FULL_L2_radar, full of Level 2 geophysical parameters are defined. 
!

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_scan
 use module_simulator
 use module_opt_micro
 use module_radar
 use module_gmp_4ice
 use module_radar_power
 use module_surface
 use module_mpi

#if HDF >= 1
 use HDF5
 USE ISO_C_BINDING
#endif

 implicit none
!--------------------------------------------------------------------------------------------------
!                   = Goddard Satellite Data Simulator Unit =
!                      - Satellite Orbit & Scanning Module -
!
! Comments:  This module computes satellite orbit and scanning points.
!            Orbit assume Kepplarian orbit, thus several Kepler parameters
!            are required to initialize satellite orbit pattern. Scanning
!            pattern are all unique to each instruments. Scan point technically
!            means that center of FOV location on Earth's reference ellipsoid. Earth's terrain
!            is not taken account for.  3D vector is based on ECEF cordinate. 
!
!  ECEF cordinate (3D vector):
!     ECEF stands for Earth-Centered, Earth-Fixed, and is a Cartesian coordinate system. 
!     ECEF represents positions as an X, Y, and Z coordinate. The point (0,0,0) is defined 
!     as the center of mass of the earth (earth's core). Its axes are aligned with the International 
!     Reference Pole (IRP) and International Reference Meridian (IRM) that are fixed with respect 
!     to the surface of the Earth.
!
!     The z-axis is pointing towards the north but it does not coincide exactly with the 
!     instantaneous Earth rotational axis. 
!     The x-axis intersects the sphere of the earth at 0° latitud(Equator) and 0° longitude Greenwich). 
!     This means that ECEF rotates with the earth and therefore, coordinates of a point fixed on the 
!     surface of the earth do not change.
!
! scan_simulator( principle driver, public access)
!  |
!  |- init_scan_simulator  (initialize various constant parameters)
!  |    |- init_orbit   (initialize orbit/keppler parameters)
!  |    |- init_sensor  (initialize sensor scanning parameters)
!  |    |- init_time_one_cycle  (initialize orbit time parameters for one cycle)
!  |    |- init_overpass (initialize orbit location for a given overpassing point)
!  |    |- allocate_scan (allocate scanning parameters )
!  |
!  |- orbit_scan_simulator (propagate satellite and scan points in time steps)
!  |    |- nadir_angle_ecross (predict nadir angle of cross-trackng instrument)
!  |    |- right_ascension (estimate earth's rotation angle)
!  |    |- solve_true_anomalies (estimate true anomaly from mean anomaly)
!  |    |- cartesian2latlon (get latitude and longitude of satellite)
!  |    |- satellite_vector (estimate satellite vector position)
!  |    |    |- rotation_matrix (position matrix rotation)
!  |    |- earth_radius (estimate distance from earth core to reference ellipsoid)
!  |    |- scan_vector (estimate fov vector position)
!  |    |- cartesian2latlon (get latitude and longitude of fov)
!  |    |- write_orbit_map (write out scan orbit output)
!  |
!  |- scan_statistics (derive various statiscs from orbit/scan outputs)
!  |    |- get_arc_distance (estimate distance between two point on reference ellipoid)
!  |
!  |- scan_interpolate (Interpolate sat-fov vector upon CRM surface grid.)
!       |- latlon2cartesian (get vector of surface CRM grid point)
!       |- get_arc_distance (estimate distance between two point on reference ellipoid)
!       |- triangle_parameter (estimate angle and length of triangle )
!
!
!
!
! drive_micro_simulator_scan(driver of microwave simulator of scanning geometry)
!  |
!  |- driver_slant_path_micro (driver of 1D slant-path microwave simulator)
!  |    |- slant_path1d_micro (Find slant-path up-down profile)
!  |    |    |- vector_increment (get inremental slant vector)
!  |    |    |- find_nearest_grid (find nearest grid for slant vector)
!  |    |    |- find_mirror_position (find mirror position of satellite)
!  |    |- eddington_2way (solve microwave radiative transfer)
!  |
!  |- convolution_micro(convolve microwave Tb within beamwidth)
!  !
!  |- call write_orbital_micro_netcdf (write out netcdf orbital output )
!  |
!  |- call write_orbital_micro_latlon (write out binary orbital output on lat-lon cordinate grid)
!
!
!
!
! History:
! 02/2011  Toshi Matsui@NASA GSFC : Adapted to G-SDSU lidar-visir modules. 
! 09/2011  Toshi Matsui@NASA GSFC : Adapted to G-SDSU microwave optical module. 
! 12/2010  Toshi Matsui@NASA GSFC : Initial. 
!
!
! References:
!
!  Brooks, D. R., 1977: An introduction to orbit dynamics and its application to satellite 
!       based earth monitoring missions. NASA Ref. Publication 1009 Washington, DC, 80 pp. 
!       [NTIS No. N78-12113]
!
!  Kidder, S. Q., and T. H. Vonder Haar, 1995: Satellite Meteorology: An Introduction, 
!         (Chapter 2: Orbits and Navigation), Academic Press, San Diego, 466 pp. 
!
!----------------------------------------------------------------------

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 private   ! encapsulate all parameters and subourtines (non accessible)

!
! public (accessible) subroutines
!
 public :: scan_simulator, & ! driver of orbit-scan simulator
           zonal_freq_cloudsat_mask, & ! create zonal cloud frequencies from CloudSat CPR
           drive_micro_simulator_scan, &
           drive_radar_simulator_scan, &
           drive_visir_simulator_scan, &
           drive_lidar_simulator_scan



! #############################################################################
! ##########################  Run-time Parameters #############################
! #############################################################################
 logical,parameter :: verbous_scan = .true.  !If true, it gives alot of comments.  
 

 character(len=50) :: satellite_name  ! GPM       : Global Precipitation Measurement 
                                      ! TRMM_PRE  : TRMM satellite pre-boost 
                                      ! TRMM_POST : TRMM satellite post boost 
                                      ! AQUA      : Aqua satellite 
                                      ! TERRA     : Terra satellite
                                      ! CALIPSO   : CALIPSO satellite
                                      ! CLOUDSAT  : CLOUDSAT satellite


 character(len=50) :: sensor_name   ! GMI_LF : GPM Microwave Imager low-frequency channel (10~89)
                                    ! GMI_HF : GPM Microwave Imager high-Frequency channel (166-183)
                                    ! TMI_LF : TRMM Microwave Imager low-frequency channel (10~37GHz)
                                    ! TMI_HF : TRMM Microwave Imager high-Frequency channel (85GHz)
                                    ! AMSR_E_LF : AMSR-E low-freq channels  (6.9~37GHz)
                                    ! AMSR_E_HF : AMSR-E high-freq channels (89GHz)
                                    ! DPR_Ku : GPM Dual-Frequency Radar Ku band (13.6GHz)
                                    ! DPR_Ka : GPM Dual Frequency Radar Ka band (35.5GHz)
                                    ! PR     : TRMM Precipitation Radar (13.6GHz)

 real(sdsu_fpd) :: overpass_lat  ! overpassing latitude [rad]
 real(sdsu_fpd) :: overpass_lon  ! overpassing longitude [rad]
 logical :: ascend_overpass         ! if true, overpass in ascending orbit 
 character(len=100) :: overpass_tag ! 

!
! output options
!
 logical,parameter :: write_grid_data    = .true.  ! if true, write out grads-gridded data
 logical,parameter :: write_station_data = .false. ! if true, write out station-format scan data in text format.


! #############################################################################
! ###################### End of Run-time Parameters ###########################
! #############################################################################

 character(len=5) :: simulator_typ  !'radar' or 'micro'

!
! Time parameters
!
 type gregorian_time
    integer :: yyyy  !year 
    integer :: mm    !month
    integer :: dd    !day
    integer :: hh    !hour
    integer :: mn    !minute
    real(sdsu_fpd) :: ss    !second
    integer :: dayofyear
    real(sdsu_fpd) :: secofday
 end type gregorian_time
 type (gregorian_time) :: &
    epoch    ,&  !epoch time
    start    ,&  !sampling start time
    finish       !sampling finish time 

 type julian_day
   real(sdsu_fpd) :: epoch   ! satellite epoch julian day [day]
   real(sdsu_fpd) :: start   ! sampling start julian day [day]
   real(sdsu_fpd) :: finish  ! sampling finish julian day [day]
 end type julian_day
 type( julian_day ) :: julian  

 real(sdsu_fpd) :: dt_one_cycle ! time for one orbit cycle [sec]
 real(sdsu_fpd) :: dt_half_cycle ! time for half orbit cycle [sec]

 
 character(len=200) :: io_file_name  !io file name
 character(len=200) :: ctl_file_name


!
! Keplarian orbit parameters
!
 type orbit_parameter
   real(sdsu_fpd) :: a          ! semi-major axis [km]
   real(sdsu_fpd) :: e          ! eccentricity [-]
   real(sdsu_fpd) :: i          ! orbit inclination angle [rad]
   real(sdsu_fpd) :: f          ! true anomalies [rad]
   real(sdsu_fpd) :: omega      ! argument of perigee [rad]
   real(sdsu_fpd) :: ohm        ! right ascension of the ascending node [rad]
   real(sdsu_fpd) :: M          ! mean anomalies [rad]
   real(sdsu_fpd) :: M_dot      ! perturbed mean motion [rad/sec]
   real(sdsu_fpd) :: M0         ! initial value of mean anomalies [rad]
   real(sdsu_fpd) :: omega0     ! initial value of omega [rad]
   real(sdsu_fpd) :: omega_dot  ! perturbation term in omega [rad/sec]
   real(sdsu_fpd) :: ohm0       ! initial value of ohm [rad]
   real(sdsu_fpd) :: ohm_dot    ! perturbation term in ohm [rad/sec] 
   real(sdsu_fpd) :: ohm_earth  ! earth's rotation angle [rad]
 end type orbit_parameter
 type ( orbit_parameter ) :: orbit

!
! satllite sensor parameters
!
 type sensor_parameter
   character(len=50) :: scan_typ !'conical' , 'ecross' , 'cross'
   character(len=50) :: direc   ! 'F' (forward looking) or 'B' (backward looking)
   integer :: mx_sample   ! maximum sample number per one scan
   real(sdsu_fpd) :: ona         ! off-nadir angle [rad]
   real(sdsu_fpd) :: evs         ! earth viewing sector (0~180deg)
   real(sdsu_fpd) :: scan_sec    ! scanning per sec  [scan/sec]
   real(sdsu_fpd) :: samp_dt     ! sensor sampling interval [sec]
 end type sensor_parameter
 type ( sensor_parameter ) :: sensor  

 real(sdsu_fpd),save :: dona_dpr ! bin angle for DPR and PR scanning [deg]

!
! output scan grid information 
!
 integer :: imax_scan    ! sample number per one scan 
 integer :: jmax_scan    ! scan number per total orbit time 
 integer :: js_scan      ! start j_scan for do-loop search
 integer :: je_scan      ! end j_scan for do-loop search
 integer :: js_sample    
 integer :: je_sample 
 type scan_parameter
   real(sdsu_fpd) :: vector_fov(3)  !fov center ECEF vector x,y,z [km]
   real(sdsu_fpd) :: vector_sat(3)  !satellite position ECEF vector x,y,z [km]
   real(sdsu_fpd) :: vel_sat(3)     !satellite velocity ECEF vector x,y,z [km/sec]
   real(sdsu_fpd) :: latlon_fov(2)  ![rad]
   real(sdsu_fpd) :: latlon_sat(2)  ![rad]
   real(sdsu_fpd) :: inci_ang       ! fov incidence angle [rad]
   real(sdsu_fpd) :: alt_sat        ! satellite altidue above reference ellipoid [km]
   integer :: i_crm, j_crm          ! collocated CRM grid i-j index
   type (gregorian_time) :: time     !gregorian time parameter
   real(sdsu_fpd) :: GPSatomic      !GPS atmoic time [sec]
 end type scan_parameter
 type(scan_parameter),allocatable :: scans(:,:)  ! output scan parameter

!
! output scan grid for microwave simulator 
!
 type scan_parameter_micro
   real(sdsu_fps) :: tb_fov(max_chan,1:2)  !microwave Tb within fov
   real(sdsu_fps) :: emis_fov(max_chan,1:2)  !microwave surface emissivity within fov
   real(sdsu_fps) :: wgt_fov(max_chan,1:2) !total weight wihtin fov (sample# x gain function)
 end type scan_parameter_micro
 type(scan_parameter_micro),allocatable :: scans_micro(:,:)  ! output scan parameter for micro 

!
! output scan grid for radar simulator 
!



!
! output scan grid for L2 parameters
!
 type scan_L2_parameter
  real(sdsu_fps) :: params(1:20)  !bundled L2 parameters
 end type scan_L2_parameter
 type(scan_L2_parameter),allocatable :: scans_L2(:,:)

!
! statistics 
!
 public :: stat_parameter
 type stat_parameter
  real(sdsu_fpd) :: mean_scan_dist    !swath-mean distance between scan lines [km]
  real(sdsu_fpd) :: max_scan_dist     !swath-max  distance between scan lines [km]
  real(sdsu_fpd) :: mean_fov_dist     !swath-mean distance between one fov sample to next. [km]
  real(sdsu_fpd) :: max_fov_dist      !swath-max  distance between one fov sample to next. [km]
  real(sdsu_fpd) :: max_sat_fov_dist  !maximum satellite-fov distance [km]
  real(sdsu_fpd) :: min_sat_fov_dist  !maximum satellite-fov distance [km]
  real(sdsu_fpd) :: mean_swath_width  !swath-mean swath length[km]
  real(sdsu_fpd) :: mean_scan_angle   !swath-mean scan angle [rad]
  real(sdsu_fpd) :: mean_inc_angle    !swath-mean incident angle [rad]
  real(sdsu_fpd) :: max_inc_angle     !swath-max incident angle [rad]
  logical        :: grid_gt_fov       !.true. if CRM grid is larger than FOV
 end type stat_parameter
 type( stat_parameter ),public :: stat  

!
! surface parameter for CRM domains 
!
!  length a = CRM grid  ~ earth core
!         b = satellite ~ earth core 
!         c = satellite ~ CRM grid
!
!  angle  alpha for a
!         beta  for b
!         gam   for c
!
 public :: surface_CRM_parameter
 type surface_CRM_parameter
    real(sdsu_fpd) :: vector_grid(3)    ! ECEF vector position of CRM surface grid  [km]
    real(sdsu_fpd) :: vector_sat(3)     ! ECEF vector position of satellite [km]
    real(sdsu_fpd) :: latlon_grid(2)    ! latitude and longitude of CRM surface grid [rad]
    real(sdsu_fpd) :: latlon_sat(2)     ! latitude and longitude of interpolated satellite position [rad]
    real(sdsu_fpd) :: a,b,c             ! length of triangle  [km]
    real(sdsu_fpd) :: alpha,beta,gam    ! angle of triangle [rad]
    logical        :: within_swath      ! if within swath, .true.
    real(sdsu_fps) :: umu_local         ! cosine of local incident angle [-]
 end type surface_CRM_parameter
 type ( surface_CRM_parameter ), allocatable, dimension(:,:),save, public :: surface_CRM !surface parameter

!
!
!
 type location_slant
  real :: ir
  real :: jr
 end type location_slant
 type ( location_slant ), allocatable, dimension(:,:) :: bot_slant !bottom locatio

!
! grads orbit map output 
!
 real(sdsu_fpd),parameter :: grid_space=0.25d0 ! grid spacing of output [deg]
 integer,parameter :: io_station = 200  ! output index (text output)
 integer :: j_grid, i_grid              ! latitude position
 integer :: jmax,imax                   ! output array bounds
 integer,allocatable :: cnt_scan(:,:)   ! count of satellite scanning (output)
 integer,allocatable :: cnt_sat(:,:)    ! count of satellite overpassing (output)
 real(sdsu_fpd),allocatable :: time_scan(:,:)  ! timing of satellite scanning [sec]
 real(sdsu_fpd),allocatable :: time_sat(:,:)   ! timing of satellite overpassing [sec]

!
! some universal constants
!
 real(sdsu_fpd),parameter :: dpi=3.14159265358979323846d0  ! pai 
 real(sdsu_fpd),parameter :: degrad=dpi/180.d0             ! conversion from degree to radiance
 real(sdsu_fpd),parameter :: day_sec = 86400.d0           ! total daylength in sec [sec]
 real(sdsu_fpd),parameter :: earth_rotate_1day = 360.d0 + 360.d0/365.2422d0  ! [deg]
 real(sdsu_fpd),parameter :: vector_ec(3) = (/0.d0,0.d0,0.d0/) ! ECEF position of earth's core. [km]
 real(sdsu_fpd),parameter :: undefined_scan = -9999.d-10

!
! radar orbital output array
!
 type  radar_orbital_output
   real(sdsu_fps) :: echo       !  radar echo [dBZ]
 end type radar_orbital_output
 type( radar_orbital_output ),allocatable :: radar_orbital(:,:,:)

!
! radar range parameter
!
 integer :: irange_max            ! absolute range bin maximum 

 type radar_range_parameter
    real(sdsu_fpd)   :: hgt       ! above-ellipsoid height normal to ellipsoid  [m] 
    real(sdsu_fpd)   :: distance  ! range distance [m]
    real(sdsu_fpd)   :: dhgt      ! range length [m] 
 end type radar_range_parameter
 type ( radar_range_parameter ), allocatable, dimension(:,:,:) :: radar_range3d !atmos parameter

 type radar_bundle_parameter
  real(sdsu_fps) :: gain      ! normalized gain function (0~1)       
  real(sdsu_fps) :: Zt        ! non-attenuating equiv. reflectivity factor [Z] = [mm6/m3] = [um3]
  real(sdsu_fps) :: Zm        ! attenuating     equiv. reflectivity factor [Z] = [mm6/m3] = [um3]
  real(sdsu_fps) :: tau       ! path-integrated optical depth [-]
  real(sdsu_fps) :: iwc       ! cloud ice water content (non-raining) [g/m3]
  real(sdsu_fps) :: lwc       ! cloud liquid water content (non-raining) [g/m3]
  real(sdsu_fps) :: rwc       ! rain water content (raining) [g/m3]
  real(sdsu_fps) :: rain_rate ! rain rate (raining) [mm/hr]

#ifdef FULL_L2_RADAR
  real(sdsu_fps) :: press     ! layer pressure [hPa] 
  real(sdsu_fps) :: t_air     ! layer air temperature [degK]
  real(sdsu_fps) :: sh        ! specific humidity [g/g] (derived in get_others)
  real(sdsu_fps) :: omega     ! vertical velocity (m/s)
  real(sdsu_fps) :: lfrac     ! liquid fraction (raining) [-]
  real(sdsu_fps) :: bulk_den  ! bulk density (raining) [g/m3]
  real(sdsu_fps) :: mmt(0:6)  ! moments [m^i / m3]  (i=0~6)
#endif

  real(sdsu_fpd) :: hgt       ! above-ellipsoid height normal to ellipsoid  [m] 
  real(sdsu_fpd) :: distance  ! range distance [m]
!  real(sdsu_fpd) :: dhgt      ! range length [m] 
 end type radar_bundle_parameter
 type( radar_bundle_parameter ), allocatable :: radar_bundle3d(:,:,:)
 type( radar_bundle_parameter ), allocatable :: scans_radar(:,:,:)  ! output scan parameter for radar



 type radar_bounds
  integer :: ks    !starting k index (TOA)
  integer :: ke    !ending   k index (at reference ellipsoid)
  integer :: k_srt !SRT level of k index
  integer :: k_piae! Effective PIA level of k index 
 end type radar_bounds
 type (radar_bounds), allocatable :: bounds(:,:) ! radar vertical bounds


!
! lidar orbital output array
!
 type  lidar_orbital_output 
   real(sdsu_fps) :: att_B       ! total attenuating backscattering [km-1 str-1]
   real(sdsu_fps) :: att_p22_B   ! perpendicular polarization component [km-1 str-1]
   real(sdsu_fps) :: lidar_sback ! total backscattering  [km-1 str-1]
   real(sdsu_fps) :: lidar_ratio ! lidar ratio
   real(sdsu_fps) :: optical_depth ! lidar optical depth [-]
 end type lidar_orbital_output
 type( lidar_orbital_output ),allocatable :: lidar_orbital(:,:,:)
 real(sdsu_fps),allocatable :: lidar_l2bundle(:,:,:)


!
! Bundled Level-2 output 
!
 real(sdsu_fps),allocatable :: L2bundle2d(:,:,:)
 integer,parameter :: nparam_L2 = 11  !maximum parameter number in L2 output
 real(sdsu_fps),allocatable :: L2_out_fov(:,:,:)  !bundled L2 parameters

 integer,parameter :: nparam_L2_lidar = 9
 real(sdsu_fps),allocatable :: L2bundle3d(:,:,:,:)

 real(sdsu_fps),allocatable :: domain_stag_hgt(:,:,:) 
 real(sdsu_fps),allocatable :: domain_dhgt    (:,:,:)


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  scan_simulator(typ, sat, inst, overpass_info)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  driver satellite orbit and sensor scanning subroutines.
!
! History:
! 01/2012  Toshi Matsui@NASA GSFC   : separate as preparation called before module_opt
! 11/2011  Toshi Matsui@NASA GSFC   : Adapted to G-SDSU
! 10/2010  Toshi Matsui@NASA GSFC   : initial
!
! References:
!----------------------------------------------------------------------
 character(len=*),intent(in) :: typ
 character(len=20),intent(in) :: sat  !satellite name
 character(len=20),intent(in) :: inst  !instrument name
 type(overpass_parameter),intent(in) :: overpass_info 

 if(masterproc) print*,''
 if(masterproc) print*,'MSG scan_simulator: Start orbit and scan simulator'

 !
 ! transfer G-SDSU configured input parameters to module parameters.
 !
 simulator_typ  = trim(typ)
 satellite_name = sat
 sensor_name    = inst

 overpass_lat = DBLE(overpass_info%lat)*degrad  ![rad]
 overpass_lon = DBLE(overpass_info%lon)*degrad  ![rad]
 ascend_overpass = overpass_info%ascend  !logical  
 overpass_tag    = overpass_info%tag     !overpass tag (character)


 !
 ! initialize atellite-sensor scan simulator
 !
 call init_scan_simulator

 !
 ! running satellite-sensor orbit-scan simulator 
 ! 
 call orbit_scan_simulator

 !
 ! satellite-sensor scan statistics
 !
 call scan_statistics

 !
 ! intepolate satellite scan sample onto CRM domain 
 !
 call scan_interpolate

 return
 end subroutine scan_simulator

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine drive_micro_simulator_scan
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  driver satellite orbit and sensor scanning subroutines.
!
! History:
! 01/2012  Toshi Matsui@NASA GSFC   : separated from scan simulator. 
!
! References:
!----------------------------------------------------------------------

!
! compute microwave Tb along 2way slant path
!
 call driver_slant_path_micro

!
! convolve Tb withi fov beamwidth
! 
 call convolution_micro

#if MPI == 2
 if(masterproc) then
#endif

 !
 ! write out orbital data in netcdf format
 !
 call write_orbital_micro_netcdf

 !
 ! writge out orbital data in latlon grads format
 !
 call write_orbital_micro_latlon

 !
 ! writge out full-orbit geolocation in netcdf
 !
 call write_geolocation_netcdf ( trim(scan_type_micro), 1, imax_scan )

 !
 ! write out FOV images in binary file.  (for presentation purpose)
 ! 
 call write_out_micro_scan_image

!     call write_2A12_HDF5         !HDF5 format (orbital grid)
!     call write_1BGMI_HDF5        !HDF5 format  (orbital grid)

#if MPI == 2
 endif
#endif

 !
 ! dellocate dynamic array
 !
 call deallocate_all_scan

 return
 end subroutine drive_micro_simulator_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine drive_radar_simulator_scan
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  driver satellite orbit and sensor scanning subroutines.
!
! History:
! 01/2016  Toshi Matsui@NASA GSFC   : add cloud zonal freqency
! 01/2012  Toshi Matsui@NASA GSFC   : separated from scan simulator. 
!
! References:
!----------------------------------------------------------------------

!
! compute radar power along slant path
!
 call driver_slant_path_radar

!
! convolve Z within fov beamwidth
! 
 call convolution_radar


#if MPI == 2
 if(masterproc) then
#endif

 !
 ! find bounds of the range bin
 !
 call find_range_bounds

 !
 ! write out orbital radar reflectivity and height in GrADS binary format. 
 ! & in NetCDF format. 
 !
 if( trim(scan_type_radar) == 'DPR_Ka' ) then !special DPR Ka band

     ! write GrADS binary data
     call write_orbital_radar( 'DPR_Ka_MA',  1, 25 )  !Ka Matched data
     call write_orbital_radar( 'DPR_Ka_HS', 26, 49 )  !Ka (interlaced) High Sensitivity data

     ! write NetCDF data
     call write_orbital_radar_netcdf( 'DPR_Ka_MA',  1, 25 ) ! Ka Matched data
     call write_orbital_radar_netcdf( 'DPR_Ka_HS', 26, 49 ) ! Ka (interlaced) High Sensitivity data
     call write_geolocation_netcdf ( 'DPR_Ka_MA',  1, 25 )
     call write_geolocation_netcdf ( 'DPR_Ka_HS', 26, 49 )
 
 else !general case

     ! write GrADS binary data 
     call write_orbital_radar( trim(scan_type_radar), 1, imax_scan )

     if( activate_zonal_freq_cloudsat_mask ) then
         call write_composite_radar
         call zonal_freq_cloudsat_mask('proc') !derive zonal cloud mask 
     endif

     ! write NetCDF data
     call write_orbital_radar_netcdf( trim(scan_type_radar), 1, imax_scan )
     call write_geolocation_netcdf ( trim(scan_type_radar), 1, imax_scan )

 endif



#if MPI == 2
 endif
#endif

 !
 ! dellocate dynamic array
 !
 call deallocate_all_scan

 return
 end subroutine drive_radar_simulator_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine drive_visir_simulator_scan
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  driver satellite orbit and sensor scanning subroutines.
!
! History:
! 01/2012  Toshi Matsui@NASA GSFC   : separated from scan simulator. 
!
! References:
!----------------------------------------------------------------------

!
! get L2 parameters
!
 call get_bundle_2d_visir


#if MPI == 2
 if(masterproc) then
#endif

 !
 ! writge out orbital data in latlon grads format
 !
 call write_orbital_visir_latlon

 !
 ! write ou CRM radiance in CRM grid format. 
 !
 call write_orbital_visir 


#if MPI == 2
 endif
#endif

 !
 ! dellocate dynamic array
 !
 call deallocate_all_scan

 return
 end subroutine drive_visir_simulator_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine drive_lidar_simulator_scan
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  driver satellite orbit and sensor scanning subroutines for lidar.
!            Unlike radar and microwave simulator, this does not do slant-path approach. 
! History:
! 01/2012  Toshi Matsui@NASA GSFC   : separated from scan simulator. 
!
! References:
!----------------------------------------------------------------------

!
! get bundled 3D parameters
!
 call get_bundle_3d_lidar

#if MPI == 2
 if(masterproc) then
#endif

 !
 ! convolute lidar output upon orbital grid 
 !
 call orbital_sample_lidar

 !
 ! write out lidar backscatter on orbital output.  
 !
 call write_orbital_lidar

 !
 ! write out lidar backscatter composite on orbital output.  
 !
 call write_composite_lidar




#if MPI == 2
 endif
#endif

 !
 ! dellocate dynamic array
 !
 call deallocate_all_scan


 return
 end subroutine drive_lidar_simulator_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


 subroutine  init_scan_simulator
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Initialize various subroutine to set up initial condition for sensor and satellite orbit. 
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC   : Adapted to G-SDSU
! 10/2010  Toshi Matsui@NASA GSFC   : initial
!
! References:
!----------------------------------------------------------------------

 !
 ! initialize orbit parameters
 !
 call init_orbit

 !
 ! initialize sensor parameters
 !
 call init_sensor

 !
 ! initialize time parameters
 !
 ! call init_time
 call init_time_one_cycle


 !
 ! set up initial position upon the overpassiing point.
 !
 call init_overpass


 !
 ! allocate scan sample
 !
 call allocate_scan


 return
 end subroutine init_scan_simulator

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_orbit
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Initialize some parameters to set up satellite orbit pattern. 
!
! History:
! 10/2010  Toshi Matsui@NASA GSFC         : initial
!
! References:
!
!----------------------------------------------------------------------

 select case(trim(satellite_name))
 case('GPM')
   orbit%a         = 6776.14d0     ! semi-major axis [km]
   orbit%e         = 0.0d0         ! eccentricity [-]  
   orbit%i         = 65.0d0*degrad ! orbit inclination angle [rad]
 case('TRMM_PRE')   
   orbit%a         = 6729.d0       ! semi-major axis [km]
   orbit%e         = 1.923d-4      ! eccentricity [-]
   orbit%i         = 34.98d0*degrad! orbit inclination angle [rad]
 case('TRMM_POST')  
   orbit%a         = 6780.5d0      ! semi-major axis [km]
   orbit%e         = 1.923d-4      ! eccentricity [-]
   orbit%i         = 34.98d0*degrad! orbit inclination angle [rad]
 case('AQUA')
   orbit%a         = 7077.75d0     ! semi-major axis [km]
   orbit%e         = 0.001203d0    ! eccentricity [-]
   orbit%i         = 98.14d0*degrad! orbit inclination angle [rad]
 case('TERRA')
   orbit%a         = 7085.d0     ! semi-major axis [km]
   orbit%e         = 0.0015d0    ! eccentricity [-]
   orbit%i         = 98.14d0*degrad! orbit inclination angle [rad]
 case('DMSP_F16')
   orbit%a         = 7235.14d0     ! semi-major axis [km]
   orbit%e         = 6.9194572d-4  ! eccentricity [-]
   orbit%i         = 98.9d0*degrad! orbit inclination angle [rad]
 case('DMSP_F17')
   orbit%a         = 7233.14d0     ! semi-major axis [km]
   orbit%e         = 0.00096       ! eccentricity [-]
   orbit%i         = 98.8d0*degrad ! orbit inclination angle [rad]
 case('CALIPSO')
   orbit%a         = 7075.2d0      ! semi-major axis [km] actual value is 7080.8, but made it identical to CLOUDSAT
   orbit%e         = 0.001203d0    ! eccentricity [-]
   orbit%i         = 98.2d0*degrad ! orbit inclination angle [rad]
 case('CLOUDSAT')
   orbit%a         = 7075.2d0      ! semi-major axis [km]
   orbit%e         = 0.001203d0    ! eccentricity [-]
   orbit%i         = 98.2d0*degrad! orbit inclination angle [rad]

 case default 
    print("(a)"),'MSG init_orbit: You put satellite name = ',trim(satellite_name)
    print("(a)"),' It is not available yet in module_scan.F'
    if(simulator_typ == 'micro')  print*,'Available satellite names are GPM,TRMM_PRE,TRMM_POST,AQUA'
    if(simulator_typ == 'radar')  print*,'Available satellite names are GPM,TRMM_PRE,TRMM_POST'
    if(simulator_typ == 'visir')  print*,'Available satellite names are TERRA, AQUA, TRMM'
    call stop_sdsu(myrank,'') 

 end select

 if(masterproc .and. verbous_scan ) then  !for domain decomposition case, print by masterproc only.
    print*, ' ------------- Satellite Orbit Spec ------------------'
    print*,  'satellite name     ', trim(satellite_name)
    print*,  'semi-major axis    ', orbit%a            ,'[km]'
    print*,  'inclination angle  ', orbit%i/degrad     ,'[deg]'
    print*, '        '
 endif

 if ( dabs(overpass_lat) > dabs(orbit%i) ) then
   print*,'MSG init_orbit: Overpassing latitude is larger than satellite inclination angle'
   print*,'Your overpassing latitude is ',overpass_lat/degrad,'[deg]'
   call stop_sdsu(myrank,'') 
 endif

!
! initialize orbit motion constant  (computing J2 purturbation terms)
!
 call init_motion_parameters( orbit )


 return
 end subroutine init_orbit

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_motion_parameters( orb )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Set up initial constant parameters for satellite orbit motion including
!           1st-order orbit purturbation term (J2 term). 
!
! History:
! 11/2010  Toshi Matsui@NASA GSFC : Initial. 
!
! References:
!  Brooks, D. R., 1977: An introduction to orbit dynamics and its application to satellite 
!       based earth monitoring missions. NASA Ref. Publication 1009 Washington, DC, 80 pp. 
!       [NTIS No. N78-12113]
!
!--------------------------------------------------------------------------------------------------
 type( orbit_parameter ),intent(inout) :: orb  !keplarian orbit parameters
 real(sdsu_fpd) :: Mo_dot ! unperturbed (Keplerian) mean motion [rad/sec]
 real(sdsu_fpd) :: p, J2_r2_p2
 real(sdsu_fpd),parameter :: mu = 398601.2d0    ! gravitational constant [km3/sec2]
 real(sdsu_fpd),parameter :: J2 = 1.6238235d-3  ! 3/2*J2 term
 real(sdsu_fpd),parameter :: re_eq = 6378.145d0 ! earth's mean equatorial radius [km]

!
! Mean unpurturbed motion
!
 Mo_dot = 1.d0 / (orb%a*sqrt(orb%a/mu))       ! eq6 [rad/sec]

 p = orb%a*(1.d0-orb%e*orb%e)  ![km]

 J2_r2_p2 = J2 * (re_eq*re_eq)/(p*p)  ! J2 bundled term

!
! addint 1st-order purturbation term into Mean Motion
!
 orb%M_dot = Mo_dot * (1.d0 + J2_r2_p2*dsqrt(1.d0-orb%e*orb%e)*(1.d0-1.5d0*dsin(orb%i)))  ! eq8  [rad/sec]

!
! Purturbation terms in omega 
!
 orb%omega_dot = J2_r2_p2 * orb%M_dot * ( 2.d0 - 2.5d0 * dsin(orb%i)*dsin(orb%i) )  ! eq5a [rad/sec]

!
! Purturbation terms in ohm
!
 orb%ohm_dot = -1.d0 * J2_r2_p2 * orb%M_dot * dcos(orb%i)    ! eq5b [rad/sec]

!
! time for one/half orbit cycle
!
 dt_one_cycle = 2.e0 * dpi / orb%M_dot   ! [sec]

 dt_half_cycle = 1.e0 * dpi / orb%M_dot   ! [sec]
 

 return
 end subroutine init_motion_parameters

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_sensor
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Initialize some parameters to set up satellite instrumental scan. 
!
! History:
! 10/2010  Toshi Matsui@NASA GSFC         : initial
!
! References:
!
!----------------------------------------------------------------------

!
! allocate L2_out_fov parameters
!
 if( allocated(L2_out_fov) ) deallocate(L2_out_fov)
 allocate( L2_out_fov(mxgridx,mxgridy,1:nparam_L2) )
 L2_out_fov = undefined

!
!
! 
 select case(trim(sensor_name))
 case('GMI_LF','GMI_LF37')
   sensor%scan_typ    = 'conical'  ! scan type 'conical' or 'ecross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         = 48.5*degrad   ! off-nadir angle [rad]
!   sensor%evs         = 140.0*degrad  ! earth viewing sector (0~180deg)  (pre-launch)
   sensor%evs         = 153.0*degrad  ! earth viewing sector (0~180deg)
   sensor%scan_sec    = 32.d0/60.d0! rotation per sec [ss]  (32 rpm)
   sensor%samp_dt     = 3.6d-3     ! sensor sampling interval [sec] (3.6ms)

 case('GMI_HF')
   sensor%scan_typ    = 'conical'  ! scan type 'conical' or 'ecross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         = 45.36*degrad  ! off-nadir angle [rad]
!   sensor%evs         = 140.0*degrad  ! earth viewing sector (0~180deg)  (pre-launch)
   sensor%evs         = 153.0*degrad  ! earth viewing sector (0~180deg) (post -launch)
   sensor%scan_sec    = 32.d0/60.d0! rotation per sec [ss]  (32 rpm)
   sensor%samp_dt     = 3.6d-3     ! sensor sampling interval [sec] (3.6 ms)
 case('TMI_LF')
   sensor%scan_typ    = 'conical'  ! scan type 'conical' or 'ecross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         = 49.00*degrad  ! off-nadir angle [rad]
   sensor%evs         = 130.0*degrad  ! earth viewing sector (0~180deg)
   sensor%scan_sec    = 31.6d0/60.d0! rotation per sec [ss]  (31.6 rpm)
   sensor%samp_dt     = 6.6d-3     ! sensor sampling interval [sec] (6.6 ms)
 case('TMI_HF')
   sensor%scan_typ    = 'conical'  ! scan type 'conical' or 'ecross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         = 49.00*degrad  ! off-nadir angle [rad]
   sensor%evs         = 130.0*degrad  ! earth viewing sector (0~180deg)
   sensor%scan_sec    = 31.6d0/60.d0! rotation per sec [ss]  (31.6 rpm)
   sensor%samp_dt     = 3.3d-3     ! sensor sampling interval [sec] (3.3 ms)
 case('AMSR_E_LF')
   sensor%scan_typ    = 'conical'  ! scan type 'conical' or 'ecross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         = 47.4d0*degrad  ! off-nadir angle [rad]
   sensor%evs         = 122.0*degrad  ! earth viewing sector (0~180deg)
   sensor%scan_sec    = 40.d0/60.d0! rotation per sec [ss]  (40.0 rpm)
   sensor%samp_dt     = 2.6d-3     ! sensor sampling interval [sec] (2.6 ms)
 case('AMSR_E_HF')
   sensor%scan_typ    = 'conical'  ! scan type 'conical' or 'ecross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         = 47.4d0*degrad  ! off-nadir angle [rad]
   sensor%evs         = 122.0*degrad    ! earth viewing sector (0~180deg)
   sensor%scan_sec    = 40.d0/60.d0! rotation per sec [ss]  (40.0 rpm)
   sensor%samp_dt     = 1.3d-3     ! sensor sampling interval [sec] (1.3 ms)
 case('SSMIS')
   sensor%scan_typ    = 'conical'  ! scan type 'conical' or 'ecross'
   sensor%direc       = 'B'        ! F- forward looking, B- backward looking  
   sensor%ona         = 45.d0*degrad  ! off-nadir angle [rad]
   sensor%evs         = (143.2+0.8)*degrad  ! earth viewing sector (0~180deg)
   sensor%scan_sec    = 31.6d0/60.d0  ! rotation per sec [ss]  (31.6 rpm)
   sensor%samp_dt     = 0.8d0/(sensor%scan_sec*360.d0) ! sensor sampling interval [sec] (4.2194 ms)
 case('SDPR_Ku')
! SDPR_Ku is 5x better resolution than DPR_Ku, only operated at Supercomputer. 
! right now assuming bin angle is dona_dpr/5 degree, scan per sec is 1/(0.6/5) 
! and 249 sampling at this point
   dona_dpr           = 0.71d0/5.d0    !bin angle [deg]
   sensor%mx_sample   = 25*2*5-1       ! maximum sample number per scan
   sensor%scan_typ    = 'ecross'        ! scan type 'conical' or 'cross'
   sensor%ona         = ((dona_dpr)*dble(sensor%mx_sample))/2.d0*degrad  ! off-nadir angle [rad] 
   sensor%evs         =  0.0*degrad    ! earth viewing sector (0~180deg) (does not use)
   sensor%scan_sec    = 1.d0/(0.6/5.d0)    ! scan per sec [ss] 
   sensor%samp_dt     = (0.6d0/5.d0)/dble(sensor%mx_sample) ! sensor sampling interval [sec] 

 case('DPR_Ku')
   dona_dpr = 0.71d0
   sensor%mx_sample   = 49         ! maximum sample number per scan
   sensor%scan_typ    = 'ecross'    ! scan type 'conical' or 'ecross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         = (dona_dpr*dble(sensor%mx_sample-1))/2.d0*degrad  ! off-nadir angle [rad] 
   sensor%evs         =  0.0*degrad    ! earth viewing sector (0~180deg) (does not use)
   sensor%scan_sec    = 1.d0/0.6    ! rotation per sec [ss] 
   sensor%samp_dt     = 0.6d0/dble(sensor%mx_sample) ! sensor sampling interval [sec] 
 case('DPR_Ka')
   dona_dpr = 0.71d0
   sensor%mx_sample   = 49         ! maximum sample number per scan
   sensor%scan_typ    = 'ecross'    ! scan type 'conical' or 'ecross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         = (dona_dpr*dble(sensor%mx_sample-1))/2.d0*degrad  ! off-nadir angle [rad]  (This is fake for Ka band)
   sensor%evs         =  0.0*degrad    ! earth viewing sector (0~180deg) (does not use)
   sensor%scan_sec    = 1.d0/0.6    ! rotation per sec [ss] 
   sensor%samp_dt     = 0.6d0/dble(sensor%mx_sample) ! sensor sampling interval [sec] 
 case('PR')  ! TRMM precipitation radar
   dona_dpr = 0.71d0
   sensor%mx_sample   = 49         ! maximum sample number per scan   
   sensor%scan_typ    = 'ecross'    ! scan type 'conical' or 'ecross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         = (dona_dpr*dble(sensor%mx_sample-1))/2.d0*degrad  ! off-nadir angle [rad] 
   sensor%evs         =  0.0*degrad    ! earth viewing sector (0~180deg) (does not use)
   sensor%scan_sec    = 1.d0/0.6    !  rotation per sec [ss] 
   sensor%samp_dt     = 0.6d0/dble(sensor%mx_sample) ! sensor sampling interval [sec] 

 case('MODIS_IR')  ! MODIS IR channel 
! Note original sampling per mirror image is 22 + 1
! but mirror image is further disaggregated for 1km, 500m, 250m fov. 
! so mx_sample rate is artificially modified to x10 for 1km resolution scan.
   sensor%mx_sample   = 22+1         ! maximum sample number per scan   
   sensor%scan_typ    = 'cross'    ! scan type 'conical' , 'ecross' , 'cross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         =  55.d0*degrad  ! off-nadir angle [rad] 
   sensor%evs         =  0.0*degrad    ! earth viewing sector (0~180deg) (does not use)
! Note the original rotation speed is 20.3 rpm. 
! MODIS use double-sided mirror so x 2. 
! and mirror image is further disaggregated for 1km, 500m, 250m fov.
! So scan_sec is artificially modified to x10 for 1km resolution scan.
   sensor%scan_sec    = 2.*20.3d0/60.d0  ! rotation per sec [ss]  (20.3 rpm)
   sensor%samp_dt     = (2.d0*sensor%ona/(2.d0 * dpi * sensor%scan_sec)) /dble(sensor%mx_sample-1) ! sensor sampling interval [sec] 

 case('CALIOP')  ! CALIPSO CALIOP 
   sensor%mx_sample   = 1         ! maximum sample number per scan   
   sensor%scan_typ    = 'ecross'    ! scan type 'conical' , 'ecross' , 'cross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         =  0.d0*degrad  ! off-nadir angle [rad] 
   sensor%evs         =  0.0*degrad    ! earth viewing sector (0~180deg) (does not use)
   sensor%samp_dt     = 0.16 ! sensor sampling interv (made same to CloudSat CPR)
   sensor%scan_sec    = 1./sensor%samp_dt  ! rotation per sec [ss]

 case('CPR')  ! CloudSat CPR 
   sensor%mx_sample   = 1         ! maximum sample number per scan   
   sensor%scan_typ    = 'ecross'    ! scan type 'conical' , 'ecross' , 'cross'
   sensor%direc       = 'F'        ! F- forward looking, B- backward looking  
   sensor%ona         =  0.d0*degrad  ! off-nadir angle [rad] 
   sensor%evs         =  0.0*degrad    ! earth viewing sector (0~180deg) (does not use)
   sensor%samp_dt     = 0.16 ! sensor sampling interv 
   sensor%scan_sec    = 1./sensor%samp_dt  ! rotation per sec [ss]

 case default 
    if(masterproc) then
    print*,'MSG init_sensor: You put sensor name = ',trim(satellite_name)
    print*,' It is not available yet in module_scan.F'
    if(simulator_typ == 'micro')  print*,'Available sensor names are GMI_LF,GMI_HF,TMI_LF,TMI_HF, AMSR_E, SSMIS'
    if(simulator_typ == 'radar')  print*,'Available sensor names are PR, DPR_Ku, DPR_Ka , CPR'
    endif
    call stop_sdsu(myrank,'') 

 end select

 if(masterproc .and. verbous_scan ) then
 print*,  ' ------------- Sensor Spec ------------------'
 print*,  'sensor_name         ', trim(sensor_name)
 print*,  'scan type           ', trim(sensor%scan_typ)
 print*,  'off-nadir angle     ', sensor%ona/degrad  ,'[deg]' 
 print*,  'earth viewing sector', sensor%evs/degrad  ,'[deg]'
 print*,  'rotation per sec        ', sensor%scan_sec ,'[ss]' 
 print*,  'sampling interval   ', sensor%samp_dt  ,'[sec]'
 print*,  '  '
 endif

 return
 end subroutine init_sensor

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_time
 implicit none
!--------------------------------------------------------------------------------------------------
!
! Comments: Initialize various time parameters.  
! 
! History:
! 12/2010  Toshi Matsui@NASA GSFC : Initial.  
!
! References:
!-----------------------------------------------------------------------------------------------------

!
! satellite epoch time
!
  epoch%yyyy = 1998 !year 
  epoch%mm   =    1 !month
  epoch%dd   =    1 !day
  epoch%hh   =    0 !hour
  epoch%mn   =    0 !minute
  epoch%ss   = 0.d0 !second

!
! sampling start time
!
  start%yyyy = 1998 !year 
  start%mm   =    1 !month
  start%dd   =    1 !day
  start%hh   =    0 !hour
  start%mn   =    0 !minute
  start%ss   = 0.d0 !second

!
! sampling finish time
!
  finish%yyyy = 1998 !year 
  finish%mm   =    1 !month
  finish%dd   =    1 !day
  finish%hh   =    1 !hour
  finish%mn   =    0 !minute
  finish%ss   = 0.d0 !second


!
! Convert Gregorian day to Juliay day
!
  call cal2jd( epoch,  julian%epoch  )    ! for epoch time

  call cal2jd( start,  julian%start  )    ! for sampling start time

  call cal2jd( finish, julian%finish )    ! for sampling finish time

 return
 end subroutine init_time

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine derive_time_cf ( file_name , min_since  )
 implicit none

!--------------------------------------------------------------------------------------------------
! Comments: Compute minutes since initial CF time specified in the Configure file.
! 
! History:
! 05/2015  Toshi Matsui@NASA GSFC : Initial.  
!
! References:
!-----------------------------------------------------------------------------------------------------
 character(len=50),intent(in) :: file_name   !Input file names
 integer,intent(out) :: min_since  ! minutes since starting time [min]
 integer :: file_len  !total file length
 integer :: i !looping
 integer :: cnt !count
 integer :: ystart_len  !character length that start four digit number
 integer :: start_len !character start length
 integer :: digit  !digit 
 logical :: num   ! true if input is numeric character
 integer :: n    !integer number

 real(sdsu_fps) :: init_yyyy   ! year
 real(sdsu_fps) :: init_mm     ! month
 real(sdsu_fps) :: init_dd     ! day
 real(sdsu_fps) :: init_julian ! julian day
 real(sdsu_fps) :: init_hh     ! hour
 real(sdsu_fps) :: init_nn     ! minutes
 real(sdsu_fps) :: init_ss     ! secounds


!
! initialization
!
  init_yyyy=0. ; init_mm=0. ; init_dd=0. ; init_julian=0. ; init_hh=0. ; init_nn=0. ; init_ss=0.


 ystart_len = 1

!
! Derive year
!
 start_len = ystart_len + 0 ; digit = 3
  do i = start_len, start_len+digit
      call char_to_int( file_name(i:i), n )
      init_yyyy = init_yyyy + n*(10.e0**(digit-i+start_len))
  enddo


!
! Derive month
!
    start_len = ystart_len + 5 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       init_mm = init_mm + n*(10.e0**(digit-i+start_len))
    enddo

!
! Derive day
!
    start_len = ystart_len + 8 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       init_dd = init_dd + n*(10.e0**(digit-i+start_len))
    enddo

!
! Derive hour
!
    start_len = ystart_len + 11 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       init_hh = init_hh + n*(10.e0**(digit-i+start_len))
    enddo

!
! Derive minutes
!
    start_len = ystart_len + 14 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       init_nn = init_nn + n*(10.e0**(digit-i+start_len))
    enddo



!
! Derive secounds
!
    start_len = ystart_len + 17 ; digit = 1
    do i = start_len, start_len+digit
       call char_to_int( file_name(i:i), n )
       init_ss = init_ss + n*(10.e0**(digit-i+start_len))
    enddo

!
! CF format initial time 
!
  epoch%yyyy = INT(init_yyyy) ! year 
  epoch%mm   = INT(init_mm)   ! month
  epoch%dd   = INT(init_dd)   ! day
  epoch%hh   = INT(init_hh)   ! hour
  epoch%mn   = INT(init_nn)   !minute
  epoch%ss   = DBLE(init_ss)  !second


!
! current time
!
  start%yyyy = INT(sdsu_yyyy) ! year 
  start%mm   = INT(sdsu_mm)   ! month
  start%dd   = INT(sdsu_dd)   ! day
  start%hh   = INT(sdsu_hh)   ! hour
  start%mn   = INT(sdsu_nn)   !minute
  start%ss   = DBLE(sdsu_ss)  !second


!
! Convert Gregorian day to Juliay day
!
  call cal2jd( epoch,  julian%epoch  )    ! for epoch time

  call cal2jd( start,  julian%start  )    ! for sampling start time

!
! output time in minutes since initial CF time
!
  min_since =  nint( (julian%start - julian%epoch)*24.d0*60.d0 )  !difference [min]

 end subroutine derive_time_cf

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_time_one_cycle
 implicit none
!--------------------------------------------------------------------------------------------------
!
! Comments: Initialize various time parameters based on SDSU input time, and one-orbit time. 
! 
! History:
! 11/2011  Toshi Matsui@NASA GSFC : Adapted to G-SDSU. 
! 12/2010  Toshi Matsui@NASA GSFC : Initial.  
!
! References:
!-----------------------------------------------------------------------------------------------------

!
! satellite epoch time
!
  epoch%yyyy = INT(sdsu_yyyy) ! year 
  epoch%mm   = INT(sdsu_mm)   ! month
  epoch%dd   = INT(sdsu_dd)   ! day
  epoch%hh   = INT(sdsu_hh)   ! hour
  epoch%mn   = INT(sdsu_nn)   !minute
  epoch%ss   = DBLE(sdsu_ss)  !second

!
! sampling start time
!
  start%yyyy = epoch%yyyy !year 
  start%mm   = epoch%mm   !month
  start%dd   = epoch%dd   !day
  start%hh   = epoch%hh   !hour
  start%mn   = epoch%mn   !minute
  start%ss   = epoch%ss   !second

!
! Convert Gregorian day to Juliay day
!
  call cal2jd( epoch,  julian%epoch  )    ! for epoch time

  call cal2jd( start,  julian%start  )    ! for sampling start time

  julian%finish = julian%start + dt_one_cycle/60.d0/60.d0/24.d0  ![day]

 return
 end subroutine init_time_one_cycle

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_overpass
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Estiamte initial value of right ascension of the ascending node (ohm0) toward
!            overpassing points assuming circular orbit. 
!            Note that for this routine, 
!            * satellite epoch time and sampling start time must be same. 
!
! History:
! 09/2011  Toshi Matsui@NASA GSFC : Initial
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fpd) :: re              ! distance from earth's core to elliptical surface [km]
 real(sdsu_fpd) :: ohm_overpass  ! ohm modification toward overpassing point [rad]
 real(sdsu_fpd) :: ohm_adjust    !adjusted ohm [rad]
 real(sdsu_fpd) :: a, b,theta,theta_dash, dhr,alpha,adpi
 character(len=10) :: start_point ! 'north' 'equater' 'south'
 real(sdsu_fpd) :: M0_start !mean anomalies value of starting point. [rad]
 
!
! specify some parameters depending on orbit starting poit. 
!
 start_point = 'south'

 if    (trim(start_point) == 'north') then
    M0_start = 90.d0*degrad  ![deg]
    adpi     = 0.d0*dpi
 elseif(trim(start_point) == 'equater') then
    M0_start = 0.d0*degrad  ![deg]
    adpi     = 1.5d0*dpi
 elseif(trim(start_point) == 'south') then
    M0_start = -90.d0*degrad  ![deg]
    adpi     = 1.0d0*dpi
 else
   call stop_sdsu(myrank,'MSG init_overpass: no such start_point') 
 endif


!
! overpassing  (hook these three parameter with SDSU configuration file. 
!
 if( dabs(overpass_lat) > dabs(orbit%i) ) then
     call stop_sdsu(myrank,'MSG init_overpass: overpass_lat is greater than satellite orbit inclination angle.')
 endif

! derive earth's radius
 call earth_radius( overpass_lat , re )

 a = re * dsin(overpass_lat) / dsin(orbit%i)

 theta = dacos( (a*dcos(orbit%i)) / ( re*dcos(overpass_lat))  )

!
! initial guess of ohm modification without earth's rotation 
!
  if(ascend_overpass ) then
    alpha = (adpi-theta)   ![rad]
  else  !case of descending pass
    alpha = (adpi+theta)   ![rad]
  endif

  ohm_overpass = dmod(  overpass_lon - alpha , 2.d0*dpi)  ![rad]

!
! however during the orbit of (dpi-theta) angle, earth rotate 15deg per hour, 
! so it needs further adjustment. 
!
 b = re * dcos(overpass_lat) * sin(theta)
 theta_dash = dasin( b / re )
  if(ascend_overpass ) then
    alpha = (adpi-theta_dash)   ![rad]
  else  !case of descending pass
    alpha = (adpi+theta_dash)   ![rad]
  endif
 dhr = dt_one_cycle * alpha/(2.d0*dpi) / 60.d0 / 60.d0  ![hr]
 ohm_adjust = dhr/24.d0 * (earth_rotate_1day*degrad) !adjusting ohm [rad]

!
! earth rotation angle at julian start time. 
!
  call right_ascension ( julian%start, orbit%ohm_earth )  ! earth's rotation angle [rad]

!
! This is the ohm0 that overpass specified overpassing lat,lon point by ascending or descending pass. 
! ohm0: initial value of right ascension of the ascending node
!
 orbit%ohm0   = -1.d0*M0_start + orbit%ohm_earth + ohm_overpass +  ohm_adjust    !  [rad]

!
! initial value of motion constant
!
 orbit%M0     = M0_start      ! initial value of mean anomalies [rad]
 orbit%omega0 = 0.d0          ! initial value of argument of perigee [rad]

 return
 end subroutine init_overpass

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine allocate_scan
 implicit none
 integer :: ierr
 real(sdsu_fpd) :: dt_tot   !total simulation time. 


!
! total simulation time
!
 dt_tot =    (julian%finish - julian%start)*3600.d0*24.d0  ![sec]

 if( sensor%scan_typ  == 'conical' ) then

    !
    ! sample number per one scan (scan only within earth-viewing sector)
    !
    imax_scan = DINT(  (1.d0/ sensor%scan_sec) * sensor%evs / (2.d0*dpi)  /  sensor%samp_dt  )

    !
    ! scan number per total orbit time 
    !
    jmax_scan  =  DINT(  dt_tot  * sensor%scan_sec )

 elseif( sensor%scan_typ  == 'ecross' .or. &
         sensor%scan_typ  == 'cross' ) then

    !
    ! sample number per one scan ( electri scan )
    !
    imax_scan = sensor%mx_sample

    !
    ! scan number per total orbit time 
    !
    jmax_scan  =  DINT(  dt_tot  * sensor%scan_sec )

 else
  
  call stop_sdsu(myrank,'MSG: allocate_scan. No such scan_typ')

 endif

  
 !
 ! allocate output parameter
 !
 if(masterproc .and. verbous_scan ) then
 print*,'MSG alloacte_scan: allocate saple grid'
 print*,'    imax_scan=',imax_scan
 print*,'    jmax_scan=',jmax_scan
 endif

 if( .not. allocated(scans) ) then
   allocate( scans(imax_scan,jmax_scan) , &  
             stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_scan: scans allocation error -> Terminate program. ')
 endif


!
! estimate output size and allocate output array 
!
 imax = INT(360.d0/grid_space)
 jmax = INT(180.d0/grid_space) 

 if(  write_grid_data ) then
     if( .not. allocated(cnt_scan) ) then
     allocate( cnt_scan(0:imax-1,0:jmax-1), &
               cnt_sat (0:imax-1,0:jmax-1), &
              time_scan(0:imax-1,0:jmax-1),&
              time_sat (0:imax-1,0:jmax-1)  ) 
      cnt_scan(:,:) = -999 
      cnt_sat (:,:) = -999 
      time_sat      = -999.d0
      time_scan     = -999.d0
      endif
 endif

 return
 end subroutine allocate_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine cal2jd( t, jd )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Convert Gregorian day to Julian day number of day connencing at greenwich noon. 
!           Initial julian day (0.d0 day) is on Jan 1 1900 at noon at Greenwitch.  
! History:
! 11/2010  Toshi Matsui@NASA GSFC : Initial
!
! References:
!  http://www.tondering.dk/claus/cal/node3.html#SECTION003161000000000000000
!
!  Brooks, D. R., 1977: An introduction to orbit dynamics and its application to satellite 
!       based earth monitoring missions. NASA Ref. Publication 1009 Washington, DC, 80 pp. 
!       [NTIS No. N78-12113]
!-----------------------------------------------------------------------------------------------------
 type (gregorian_time),intent(in) :: t
 real(sdsu_fpd),intent(out) :: jd    ! julian day  in full (including UT) component ( from Jan 0.5 1900 ) [day]

 integer :: a,y,m    ! 
 real(sdsu_fpd) :: ut       ! universal time [day]
 integer :: jd_int   !
 real(sdsu_fpd) :: jd_day   ! julian days without UT components. 

!
! Derive Julian days from Jan 0.5 1900
!
 a = (14-t%mm)/12
 y = t%yyyy + 4800 - a
 m = t%mm + 12*a -3

 jd_int = t%dd + (153*m+2)/5 + 365*y + y/4 - y/100 + y/400 - 32045   !julian day

 jd_day = dble(jd_int) - 0.5d0  !julian day from Jan 0.5 1900

 ut  =  dble(t%hh)/24.d0 + dble(t%mn)/(24.d0*60.d0) + t%ss / (24.d0*60.d0*60.d0)  ! julian days in UT component

 jd = jd_day + ut  ! julian day  in full (including UT) component ( from Jan 0.5 1900 ) [day]

 return
 end subroutine cal2jd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine jd2cal(jd,t)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Convert Julian day to Gregorian day connencing at greenwich noon. 
! 
! History:
! 11/2010  Toshi Matsui@NASA GSFC : Initial
!
! References:
!  http://www.tondering.dk/claus/cal/node3.html#SECTION003161000000000000000
!
!  Brooks, D. R., 1977: An introduction to orbit dynamics and its application to satellite 
!       based earth monitoring missions. NASA Ref. Publication 1009 Washington, DC, 80 pp. 
!       [NTIS No. N78-12113]
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: jd        ! julian day  in full (including UT) component ( from Jan 0.5 1900 ) [day]
 type (gregorian_time),intent(out) :: t

 integer  :: a,b,c,d,e,m    ! 
 real(sdsu_fpd) :: ut       ! universal time [day]
 integer :: jd_int   !
 real(sdsu_fpd) :: jd_day   ! julian days without UT components. 
 real(sdsu_fpd) :: jd_0z
 type (gregorian_time) :: ny
 real(sdsu_fpd) :: jd_ny

 jd_0z  = jd + 0.5d0

 jd_int = INT( jd_0z )

 ut = jd_0z - INT(jd_0z)
 t%hh = INT(ut*24.d0)
 t%mn = INT( (ut-DBLE(t%hh)/24.d0)*24.d0*60.d0 )
 t%ss =  (ut-DBLE(t%hh)/24.d0-DBLE(t%mn)/(24.d0*60.d0))*24.d0*3600.d0

 t%secofday = ut*24.d0*3600.d0

 a = jd_int + 32044
! b = (4*a + 3)/146097
! c = a + (146097*b)/4
 b = 0
 c = jd_int + 32082 + 13  !don't know why, but it works....

 d = (4*c + 3) / 1461
 e = c - (1461*d)/4
 m = (5*e+2)/153


 t%yyyy = 100*b + d - 4800 + (m/10)

 t%dd = e-(153*m+2)/5 + 1
 t%mm = m + 3 - 12 * (m/10)

!
! new year
! 
 ny%yyyy = t%yyyy
 ny%mm   =    1 !month
 ny%dd   =    1 !day
 ny%hh   =    0 !hour
 ny%mn   =    0 !minute
 ny%ss   = 0.d0 !second
 call cal2jd( ny, jd_ny )
 t%dayofyear = INT( jd - jd_ny + 1.d0 )

 return
 end subroutine jd2cal

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine jd2GPS_atomic_time(jd, GPSatomic)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Convert Julian day to GPS atomic time [sec] since 00UTC, 6, Jan 1980. 
! 
! History:
! 11/2010  Toshi Matsui@NASA GSFC : Initial
!
! References:
!  http://www.tondering.dk/claus/cal/node3.html#SECTION003161000000000000000
!
!  Brooks, D. R., 1977: An introduction to orbit dynamics and its application to satellite 
!       based earth monitoring missions. NASA Ref. Publication 1009 Washington, DC, 80 pp. 
!       [NTIS No. N78-12113]
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: jd         ! julian day  in full (including UT) component ( from Jan 0.5 1900 ) [day]
 real(sdsu_fpd),intent(out) :: GPSatomic ! GPS atomic time since 00UTC, 6, Jan 1980 [sec]

 type (gregorian_time) :: GPS ! gregoria time of initial GPM atomic time
 real(sdsu_fpd) :: jd_GPSatomic ! julidan day of initial GPM atomic time

!
! get julian day of 00UTC, 6, Jan 1980
!
 GPS%yyyy = 1980 !year
 GPS%mm   =    1 !month
 GPS%dd   =    6 !day
 GPS%hh   =    0 !hour
 GPS%mn   =    0 !minute
 GPS%ss   = 0.d0 !second
 call cal2jd( GPS, jd_GPSatomic )

!
! output 
!

 GPSatomic = (jd - jd_GPSatomic)*24.d0*3600.d0  ! [sec]


 return
 end subroutine jd2GPS_atomic_time

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine orbit_scan_simulator
 implicit none
!--------------------------------------------------------------------------------------------------
!
! Comments: This routine simulate satellite orbit and scanning point in terms of lat and lon cordinate. 
! 
! History:
! 12/2010  Toshi Matsui@NASA GSFC : Initial.  
!
! References:
!  Brooks, D. R., 1977: An introduction to orbit dynamics and its application to satellite 
!       based earth monitoring missions. NASA Ref. Publication 1009 Washington, DC, 80 pp. 
!       [NTIS No. N78-12113]
!
!  Kidder, S. Q., and T. H. Vonder Haar, 1995: Satellite Meteorology: An Introduction, 
!         (Chapter 2: Orbits and Navigation), Academic Press, San Diego, 466 pp. 
!
!-----------------------------------------------------------------------------------------------------
 integer :: i,j               ! looping
 integer :: j_scan_old         ! previous count of sanning
 integer :: i_scan, j_scan

 real(sdsu_fpd) :: lon_sat,lat_sat   ! satellite position in lat lon [rad]
 real(sdsu_fpd) :: lat_fov,lon_fov ! scanning fov positin in lat lon [rad]
 real(sdsu_fpd) :: conical_angle     ! angle of conical scanning angle [rad]
 real(sdsu_fpd) :: scan_speed        ! scan angular velocity [rad/sec] 
 real(sdsu_fpd) :: nadir_angle       ! nadir angle for cross tracking scanning [rad]
 real(sdsu_fpd) :: jd                ! scanning time in Julian Day [day]
 real(sdsu_fpd) :: dt                ! progress time since jd_epoch [sec]
 real(sdsu_fpd) :: re                ! distance between reference ellipsoid and earth's core [km]
 real(sdsu_fpd) :: rs                ! distance between satellite and earth's core [km]
 real(sdsu_fpd) :: beta              ! scan incidence angle [rad]
 real(sdsu_fpd) :: v_sat(3)          ! satellite position vector (x,y,z) in ECEF cartesizn cordinate [km]
 real(sdsu_fpd) :: v_fov(3)          ! scanng position vector (on ellipsoid) (x,y,z) in ECEF cartesizn cordinate [km]
 real(sdsu_fpd) :: vel_sat(3)        ! satellite volocity vector (x,y,z) in ECEF cartesizn cordinate [km/sec]
 real(sdsu_fpd) :: d_angle           ! residual angle [rad]
 real(sdsu_fpd) :: t_offset          ! offset time [sec]
 real(sdsu_fpd) :: ang_thresh        ! threshold angle [rad]


 if(masterproc) print*,'MSG orbit_scan_simulator: Compute orbit track and ifov point.'

!
! open station data output
!
 if(write_station_data) then
    io_file_name = trim(sdsu_dir_output)//trim(satellite_name)//'_'//trim(sensor_name)//'.txt'
    open(io_station, file= trim(io_file_name)) 
    print*,'writing text data -->',trim(io_file_name)
 endif

!
! estimate scan angular velocity
!
  if ( sensor%scan_typ == 'conical' ) then

     scan_speed = 2.d0 * dpi         * sensor%scan_sec   ! scan angular velocity [rad/sec]
     if( sensor%direc == 'F' ) ang_thresh = 270.d0*degrad - sensor%evs/2.d0   !forward scanning  
     if( sensor%direc == 'B' ) ang_thresh =  90.d0*degrad - sensor%evs/2.d0   !forward scanning  

  elseif( sensor%scan_typ == 'ecross' ) then

     scan_speed = 2.d0 * sensor%ona  * sensor%scan_sec   ! scan angular velocity [rad/sec]

  elseif( sensor%scan_typ == 'cross' ) then

     scan_speed = 2.d0 * dpi         * sensor%scan_sec  !scan angular velocity [rad/sec]
     ang_thresh =  sensor%ona 

  endif

!
! initialize sample counting
!
  i_scan = 0  ; j_scan_old = 0  ; j_scan = 1

!
! propagating satellite 
!
 jd = julian%start ! initial time for observation [mjd]

 DO_PROPAGATE: do  !------------------------------------
    
   !
   ! update satellite scanning time by djd 
   !
   jd  = jd + sensor%samp_dt / day_sec   ! scanning time in Julian Day [day]              
   dt  = (jd - julian%epoch)*day_sec     ! progress time since julian%epoch  [sec]

    select case(trim(sensor%scan_typ))

    case('conical') ! (counter clock wise rotation for GMI&TMI, assuming starting from dpi (when t-julian%epoch=0)

      j_scan = 1 + dint( scan_speed * dt / (2.d0*dpi) )   !
      if(j_scan /= j_scan_old)  i_scan = 0                ! initialize i_scan for new nscan
      conical_angle = dmod( -1.d0*scan_speed * dt , 2.d0*dpi       )  ! conical angle [rad] counter clockwise
!      conical_angle = dmod( scan_speed * dt , 2.d0*dpi       )  ! conical angle [rad] clockwise

    case('ecross')

      if( i_scan ==  sensor%mx_sample ) then       ! maximum sample number per scan
          j_scan = 1 + j_scan
          i_scan = 0
      endif
      conical_angle = 0.d0

    case('cross')
  
      j_scan = 1 + dint( scan_speed * dt / (2.d0*dpi) )
      if(j_scan /= j_scan_old)  i_scan = 0                ! initialize i_scan for new nscan
      nadir_angle = dmod( 1.d0*scan_speed * dt - 0.5d0 *dpi, 2.d0*dpi       )  ! nadir angle [rad]
      if(nadir_angle > dpi ) nadir_angle =  nadir_angle  - 2.d0 * dpi    ! adjust nadir angle [rad] 
      conical_angle = 0.d0

    case default ; call stop_sdsu(myrank,'MSG orbit_scan_simulator : There is no such sensor%scan_typ')
    end select

   !
   ! make sure exist, when scan line exceed maximum value
   !
   if( j_scan > jmax_scan ) exit DO_PROPAGATE   

   WITHIN_VIEW_SECTOR: if( (sensor%scan_typ == 'ecross')   &
                                                 .or.   &
                           ( (sensor%scan_typ == 'cross') .and.  &
                              DABS(nadir_angle) <= ang_thresh) &
                                                 .or.   &
                           ( (sensor%scan_typ == 'conical')  .and.  &
                              DABS(conical_angle) >= ang_thresh ) ) then  ! b

      i_scan = i_scan + 1  !count of sampling

      !
      ! adjust conical angle and time for 1st sample
      !
      if((sensor%scan_typ == 'conical')  .and. i_scan == 1) then
          d_angle = DABS(conical_angle) - ang_thresh  ! residual angle [rad]
          t_offset = d_angle/scan_speed               ! offset time [sec]
          jd = jd - t_offset/day_sec                  ! adjust scanning time in Julian Day using offset time [day]
          dt  = (jd - julian%epoch)*day_sec           ! re-calculate progress time since julian%epoch  [sec]
          conical_angle = dmod( -1.d0*scan_speed * dt , 2.d0*dpi       )  ! adjust conical angle [rad]
      endif

      if((sensor%scan_typ == 'cross')  .and. i_scan == 1) then

          d_angle  = nadir_angle + DABS(ang_thresh)   ! residual angle [rad]
          t_offset = d_angle/scan_speed               ! offset time [sec]
          jd       = jd - t_offset/day_sec            ! adjust scanning time in Julian Day using offset time [day]
          dt       = (jd - julian%epoch)*day_sec      ! re-calculate progress time since julian%epoch  [sec]
          nadir_angle = dmod( 1.d0*scan_speed * dt - 0.5d0 *dpi, 2.d0*dpi       )  ! nadir angle [rad]
          if(nadir_angle > dpi ) nadir_angle =  nadir_angle  - 2.d0 * dpi    ! adjust nadir angle [rad] 

      endif

      !
      ! make sure if i_scan exceed sampling number per scan, skip rest of process. 
      !
      if(i_scan > imax_scan) cycle DO_PROPAGATE


      !
      ! ecross tracking instrument nadir angle
      !
      if ( sensor%scan_typ == 'ecross' ) then
         call nadir_angle_ecross(i_scan, scan_speed, nadir_angle)
      endif

      !
      ! Derive orbit parameter for a given time proceeding after epoch time.
      !
      call right_ascension ( jd, orbit%ohm_earth )                    ! earth's rotation angle [rad]
      orbit%omega = orbit%omega0 + orbit%omega_dot  * dt              ! argument of perigee [rad]
      orbit%ohm   = orbit%ohm0 - orbit%ohm_earth + orbit%ohm_dot * dt ! right ascension of the ascending node [rad]
      orbit%M     = ( orbit%M0 + orbit%M_dot * dt )                   ! mean anomalies [rad]

      call solve_true_anomalie(orbit%e, orbit%M, orbit%f)       ! true anomalies, f [rad]

      !
      ! estimate satellite position vector through matrix rotation method
      !
      call satellite_vector(orbit, v_sat)

      !
      ! estimate satellite velocity vector 
      !
      call satellite_velocity(orbit, vel_sat)

      !
      ! derive satellite altitude from earth's core
      !
      call get_distance( v_sat , vector_ec, rs)

      !
      ! satellite position in lat-lon cordinate
      !
      call cartesian2latlon(v_sat,lat_sat,lon_sat)

      !
      ! store into gridded array for output purpose. 
      !
      if(  write_grid_data ) then
           call grid_point(lat_sat,lon_sat,j_grid,i_grid) ! satellite position
           cnt_sat(i_grid, j_grid) = 1               ! satellite overpassing count for grid [-]
           time_sat(i_grid, j_grid) = dt             ! satellite overpassing time in progress [sec]
      endif

      !
      ! derive earth's radius as oblate spheroid
      !
      call earth_radius( lat_sat, re )

      !
      ! scanning perturbation
      !
      call scan_vector(orbit, re, conical_angle, nadir_angle, v_fov, beta)

      !
      ! scan position in lat-lon cordinate
      !
      call cartesian2latlon(v_fov,lat_fov,lon_fov)

      !
      ! store in global map for output 
      !
      if(  write_grid_data ) then
         call grid_point(lat_fov,lon_fov,j_grid,i_grid) ! scan position
         cnt_scan(i_grid, j_grid) = 1                   ! scan count for grid [-]
         time_scan(i_grid, j_grid) = dt                 ! scan time in progress [sec]
      endif

      !
      ! output
      !
      scans(i_scan,j_scan)%vector_fov(1:3) = v_fov  (1:3)![km]
      scans(i_scan,j_scan)%vector_sat(1:3) = v_sat  (1:3)![km]
      scans(i_scan,j_scan)%vel_sat   (1:3) = vel_sat(1:3)![km/sec]
      scans(i_scan,j_scan)%latlon_fov(1) = lat_fov       ![rad]
      scans(i_scan,j_scan)%latlon_fov(2) = lon_fov       ![rad]
      scans(i_scan,j_scan)%latlon_sat(1) = lat_sat       ![rad]
      scans(i_scan,j_scan)%latlon_sat(2) = lon_sat       ![rad]
      scans(i_scan,j_scan)%inci_ang      = DABS(beta)    ![rad]
      scans(i_scan,j_scan)%alt_sat       = rs - re       ![km]
      call jd2cal( jd, scans(i_scan,j_scan)%time )       !gregorian time [yyyymmddhhmnss]
      call jd2GPS_atomic_time(jd, scans(i_scan,j_scan)%GPSatomic) !GPS atmoic time [sec]

      !
      ! initialize other parameters
      !
      scans(i_scan,j_scan)%i_crm = undefined_i2       
      scans(i_scan,j_scan)%j_crm = undefined_i2       

      !
      ! output text data
      !
      if( write_station_data ) &
      write(io_station,FMT='(2I10, 2F11.5)') i_scan, j_scan, real(lon_fov/degrad), real(lat_fov/degrad)

   endif WITHIN_VIEW_SECTOR

   j_scan_old = j_scan  !preivious scan number

 enddo DO_PROPAGATE !-----------------------------------------------

!
! write gridded data in Grads binary format. 
!

 if(  write_grid_data ) then
#if MPI == 2
     if(masterproc) &  !for domain decomposition case, write by masterproc only.
#endif
     call write_orbit_map
 endif

 return
 end subroutine orbit_scan_simulator

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_orbit_map
 implicit none
 integer :: i,j
 integer,parameter :: io_ctl     = 706  ! output idnex (grads control file)
 integer,parameter :: io_grid    = 705  ! output index (grads gridded data)
 logical,parameter :: write_satellite_position = .false.

!
! scan position
!
    print*,'MSG: write_orbit_map. Write GrADS gridded data'
    io_file_name = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.FOV_'// &
                   trim(satellite_name)//'_'//trim(sensor_name)//trim(output_suffix)//'.bin'
    print("(a)"),'  dump--> ',trim(io_file_name)
    open (io_grid,file=trim(io_file_name),access='direct', status='replace', recl=jmax*imax*4)
    write(io_grid,rec=1) (( real(time_scan(i,j)) , i=0,imax-1), j=0,jmax-1)
    close(io_grid)

    io_file_name  = sdsu_inp_name(1:efile_len)//'.FOV_'// &
                   trim(satellite_name)//'_'//trim(sensor_name)//trim(output_suffix)//'.bin'
    ctl_file_name = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.FOV_'// &
                   trim(satellite_name)//'_'//trim(sensor_name)//trim(output_suffix)//'.ctl'
    open (io_ctl,file=trim(ctl_file_name)) 
    write(io_ctl,"(a)") 'dset   '//trim(io_file_name)
    write(io_ctl,*) 'title  '//trim(sensor_name)//' fov points'
    write(io_ctl,*) 'undef  -999.'
    write(io_ctl,*) 'xdef   ',imax,' linear 0.00',grid_space
    write(io_ctl,*) 'ydef   ',jmax,' linear ',-90.+grid_space/2.,grid_space
    write(io_ctl,*) 'zdef    1 linear 0. 1 '
    write(io_ctl,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io_ctl,*) 'vars 1'
    write(io_ctl,*) 'time  0 0  scan time [sec]'
    write(io_ctl,*) 'endvars'
    close(io_ctl)

 if(write_satellite_position) then
    io_file_name = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.ORBIT_'//&
                   trim(satellite_name)//'_'//trim(sensor_name)//trim(output_suffix)//'.bin'
    print("(a)"),'  dump--> ',trim(io_file_name)
    open (io_grid,file=trim(io_file_name),access='direct', status='replace', recl=jmax*imax*4)
    write(io_grid,rec=1) (( real(time_sat(i,j)) , i=0,imax-1), j=0,jmax-1)
    close(io_grid)

    io_file_name = sdsu_inp_name(1:efile_len)//'.ORBIT_'//&
                   trim(satellite_name)//'_'//trim(sensor_name)//trim(output_suffix)//'.bin'
    ctl_file_name = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.sat_position_'//&
                   trim(satellite_name)//'_'//trim(sensor_name)//trim(output_suffix)//'.ctl'
    open (io_ctl,file=trim(ctl_file_name))
    write(io_ctl,"(a)") 'dset   '//trim(io_file_name)
    write(io_ctl,*) 'title  '//trim(sensor_name)//' satellite points'
    write(io_ctl,*) 'undef  -999.'
    write(io_ctl,*) 'xdef   ',imax,' linear 0.00',grid_space
    write(io_ctl,*) 'ydef   ',jmax,' linear ',-90.+grid_space/2.,grid_space
    write(io_ctl,*) 'zdef    1 linear 0. 1 '
    write(io_ctl,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io_ctl,*) 'vars 1'
    write(io_ctl,*) 'time  0 0  scan time [sec]'
    write(io_ctl,*) 'endvars'
    close(io_ctl)
 endif

 return
 end subroutine write_orbit_map

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine solve_true_anomalie(ec, M, f_true)
 implicit none
!----------------------------------------------------------------------------------
! Comments:  Estimate true anomalies for a given eccentricity and mean anomalies through
!            either iterative method or analytic solution.  
!
! History:
! 12/2010  Toshi Matsui@NASA GSFC   : Initial.
!
! References:
! Meeus, Jean, 1991, Astronomical Algorithms, Willmann-Bell, 
!  1st English edition, 1991, ISBN 0-943396-35-2
!
!----------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: ec      ! eccentricity    [-]
 real(sdsu_fpd),intent(in) :: M       ! mean anomalies [rad]
 real(sdsu_fpd),intent(out) :: f_true ! true anomalies [rad]

 integer :: i   ! looping
 real(sdsu_fpd) :: E   !eccentric anomalies [rad]
 real(sdsu_fpd) :: dE  ! 
 real(sdsu_fpd),parameter :: degrad = 2.d0*3.1415926535d0/360.d0

!
! choose solution 
! true --> iterative solution (good for all ec value) DEFAULT
! false --> solution in series expansion (only for small ec value)
!
 logical,parameter :: iterative = .true.

 if(iterative) then

   !
   ! initial value of eccentric anomalies [rad]
   ! 
    if(ec<0.8d0) then
       E=M
    else
       E=180.d0*degrad
    endif

   !
   ! iteratively estimate eccentric anomalie
   !
   dE = E - ec*dsin(M) - M
    i = 0

    do while( dabs(dE)>1.d-10 )
      E  = E - dE/(1.d0-ec*dcos(E))
      dE = E - ec*dsin(E) - M
      i  = i + 1
      if(i>500) call stop_sdsu(myrank,'MSG newton_solution_true_anomalie: too many iteration')
    enddo


   !
   ! find true anomalies (ranges -180deg ~ 180deg)
   !
    f_true = datan2( dsqrt(1.d0-ec*ec)*dsin(E) , dcos(E)-ec )  !true anomalies [rad]


 else  

!
! From Meeus 1991 (ch3) solution through series expansion (good enough for small eccentricity)
! Maybe useful to check the accuracy of your numerical iterative process for small eccentricity. 
!
    f_true = M  + (2.d0*ec-1.d0/4.d0*ec*ec*ec+5.d0/96.d0*ec*ec*ec*ec*ec)*dsin(M) &
                + (5.d0/4.d0*ec*ec-11.d0/24.d0*ec*ec*ec*ec)*dsin(2.d0*M) &
                + (13.d0/12.d0*ec*ec*ec-43.d0/64.d0*ec*ec*ec*ec*ec)*dsin(3.d0*M) &
                + 103.d0/96.d0*ec*ec*ec*ec*dsin(4.d0*M) &
                + 1097.d0/960.d0*ec*ec*ec*ec*ec*dsin(5.d0*M)

 endif


 return
 end subroutine solve_true_anomalie

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine satellite_velocity(o, vel_sat)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Compute satellite volocity vector in ECEF cartezian cordinate (x,y,z).
!
! History:
! 10/2011  Toshi Matsui@NASA GSFC : Initial. 
!
! References:
!  Kidder, S. Q., and T. H. Vonder Haar, 1995: Satellite Meteorology: An Introduction, 
!         (Chapter 2: Orbits and Navigation), Academic Press, San Diego, 466 pp. 
!
!-------------------------------------------------------------------------------------------------
 type( orbit_parameter ),intent(in) :: o  ! keplarian orbit parameters
 real(sdsu_fpd),intent(out) :: vel_sat(3) ! satellite velocity vector [km/sec]
 real(sdsu_fpd),parameter :: mu = 398601.2d0    ! gravitational constant [km3/sec2]
 real(sdsu_fpd) :: p,coef 

!
! Mean unpurturbed motion
!
 p = o%a*(1.d0-o%e*o%e)  ![km]
 coef = dsqrt(mu/p)      ![km/sec]


!
! x direction
!
 vel_sat(1) = coef * ( (dcos(o%f)+o%e)*(-dsin(o%omega)*dcos(o%ohm) - dcos(o%i)*dsin(o%ohm)*dcos(o%omega)) &
                       - dsin(o%f) * ( dcos(o%omega)*dcos(o%ohm) - dcos(o%i)*dsin(o%ohm)*dsin(o%omega) ) & 
                     )

 vel_sat(2) = coef * ( (dcos(o%f)+o%e)*(-dsin(o%omega)*dcos(o%ohm) + dcos(o%i)*dcos(o%ohm)*dcos(o%omega)) &
                       - dsin(o%f) * ( dcos(o%omega)*dsin(o%ohm) - dcos(o%i)*dcos(o%ohm)*dsin(o%omega) ) &
                     )

 vel_sat(3) = coef * ( (dcos(o%f)+o%e)*(dsin(o%i)*dcos(o%omega) - dsin(o%f)*dsin(o%i)*dsin(o%omega) ) )


 return
 end subroutine satellite_velocity

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


 subroutine satellite_vector(orb, v_sat)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Compute satellite position vector in ECEF cartezian cordinate (x,y,z)
!           through matrix rotation method. 
!
! History:
! 12/2010  Toshi Matsui@NASA GSFC : Initial. 
!
! References:
!  Kidder, S. Q., and T. H. Vonder Haar, 1995: Satellite Meteorology: An Introduction, 
!         (Chapter 2: Orbits and Navigation), Academic Press, San Diego, 466 pp. 
!
!-------------------------------------------------------------------------------------------------
 type( orbit_parameter ),intent(in) :: orb  ! keplarian orbit parameters
 real(sdsu_fpd),intent(out) :: v_sat(3) ! satellite position vector [km]

 real(sdsu_fpd) :: rs       ! distance from earth's core to satellite [km]
 real(sdsu_fpd),dimension(3) :: v, v_omega, v_omega_i ! position vector 

!
! distance from earth's core to satellite
!
 rs = orb%a*(1.d0-orb%e*orb%e)/(1.d0+orb%e*dcos(orb%f))  ! [km]

!
! initial position vector (the satellite orbital ellipse is assumed to lie in x-y plane)
!
 v(1) =  rs*dcos(orb%f)
 v(2) =  rs*dsin(orb%f)
 v(3) =  0.d0


!
! Rotate v about the z axis throught the argument of perigee (omega)
!
 call rotation_matrix ( 'z', orb%omega, v       , v_omega        )

!
! Rotate v_omega about the x axis throught the inclination angle (i)
!
 call rotation_matrix ( 'x',     orb%i, v_omega , v_omega_i      )

!
! Rotate v_omega_i about the z axis through the right ascension of ascending node (ohm)
!
 call rotation_matrix ( 'z',   orb%ohm, v_omega_i, v_sat         )

 return
 end subroutine satellite_vector

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine rotation_matrix ( axis, a, v , v_rot   )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Rotate position vector (x,y,z) in certain angle about x/y/z axis, and
!           convert into new (x',y',z') position vector. 
!
! History:
! 12/2010  Toshi Matsui@NASA GSFC : Initial. 
!
! References:
!-------------------------------------------------------------------------------------------------
 character(len=1),intent(in) :: axis  ! rotation axis x, y, or z
 real(sdsu_fpd),intent(in) :: a              ! rotation angle [rad]
 real(sdsu_fpd),intent(in) :: v(3)           ! original position vector (x,y,z)
 real(sdsu_fpd),intent(out) :: v_rot(3)      ! rotated position vector (x',y',z')

 real(sdsu_fpd) :: Rot(3,3)  !rotation matrix


 select case(axis)
 case('x')
! Rotation Matrix about x-axis
!           1      2        3
!
!  1        1      0        0
!  2        0  cos(a)  -sin(a) 
!  3        0  sin(a)   cos(a) 

 Rot = reshape( (/1.d0,          0.d0,      0.d0,   &
                  0.d0,       dcos(a),   dsin(a),   &
                  0.d0, -1.d0*dsin(a),   dcos(a)/), &
                  (/3,3/) )

 case('y')
! Rotation Matrix about y-axis
! Ry        1      2        3
!
!  1    cos(a)     0    sin(a) 
!  2        0      1        0 
!  3   -sin(a)     0    cos(a) 

 Rot = reshape( (/dcos(a),       0.d0,  -1.d0*dsin(a) , &
                     0.d0,       1.d0,           0.d0 , &
                  dsin(a),       0.d0,       dcos(a)/), &
                  (/3,3/) )

 case('z')
! Rotation Matrix about z-axis
! Rz        1       2       3
!
!  1    cos(a) -sin(a)      0 
!  2    sin(a)  cos(a)      0 
!  3        0       0       1 

 Rot = reshape( (/      dcos(a),    dsin(a),           0.d0 , &
                  -1.d0*dsin(a),    dcos(a),           0.d0 , &
                           0.d0,       0.d0,           1.d0/), &
                  (/3,3/) )

 case default ; call stop_sdsu(myrank,'MSG rotation_matrix: There is no such axis. Stop.')
 end select 


!
! Get rotated new position vector
!
 v_rot = matmul( Rot , v )  ! matrix cross operator (Rot x v)

 return
 end subroutine rotation_matrix

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine cartesian2latlon(v,lat,lon)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Convert cartezian cordinate (x,y,z) into lat-lon cordinate.
!           Note that (x,y,z) = (rs, 0, 0) => lat= 0deg, lon= 0deg
!                     (x,y,z) = ( 0,rs, 0) => lat= 0deg, lon=90deg
!                     (x,y,z) = ( 0, 0,rs) => lat=90deg, lon= 0deg
!
! History:
! 12/2010  Toshi Matsui@NASA GSFC : Initial. 
!
! References:
!-------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: v(3)     !3D position vector in ECEF cartesian cardinate [km]
 real(sdsu_fpd),intent(out) :: lat      ! latitude (-90deg ~ 90deg) [rad]
 real(sdsu_fpd),intent(out) :: lon      ! longitude (-180deg ~ 180deg) [rad]
 real(sdsu_fpd) :: r  !distance from earth's core [km]
 real(sdsu_fpd) :: rx,ry,rz

 rx = v(1)  ; ry = v(2) ; rz = v(3) 

 r = dsqrt(rx*rx+ry*ry+rz*rz)

 if(r == 0.d0) call stop_sdsu(myrank,'MSG cartesian2latlon: stop program distance between earth core and satellite is zero')
!
! latitude [rad]
!
 lat = dasin(rz/r)  

!
! longitude [rad]
!

!special case
 if( ry==0.d0 .and. rx==0.d0) then  
     lon = 0.d0
     return
 endif
     
!usual case
 if    ( ry >= 0.d0 .and. rx >= 0.d0 ) then
     lon =             datan(ry/rx) 
 elseif( ry >= 0.d0 .and. rx <  0.d0 ) then
     lon = dpi      + datan(ry/rx) 
 elseif( ry <  0.d0 .and. rx <  0.d0 ) then
     lon = datan(ry/rx)  - dpi
 elseif( ry <  0.d0 .and. rx >= 0.d0 ) then
     lon =  datan(ry/rx) 
 else 
  call stop_sdsu(myrank,'MSG cartesian2latlon: no such ry,rx combinations')
 endif

 return
 end subroutine cartesian2latlon

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine latlon2cartesian(lat,lon,v)
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
 real(sdsu_fpd),intent(out) :: v(3) !cartesian cordinate position [km]
 real(sdsu_fpd) :: re !distance from earth's core to ellipsoidal surface  [km]

!
! Get distance from the earth's core to the ellipoidal surface. 
!
 call earth_radius(lat, re) 

!
! Derive vector x,y,z
!
 v(1) = re*cos(lat)*cos(lon)  ! x [km]

 v(2) = re*cos(lat)*sin(lon)  ! y [km]

 v(3) = re*sin(lat)           ! z [km]


 return
 end subroutine latlon2cartesian

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_distance( v1, v2, dist)
 implicit none
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

 subroutine get_arc_distance( v1, v2, arc_dist)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Estimate arc distance between two vector points on earth's surface. 
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fpd),intent(in)  :: v1(3),v2(3)
 real(sdsu_fpd),intent(out) :: arc_dist
 real(sdsu_fpd) :: a,b,c  ![km]
 real(sdsu_fpd) :: alpha ![rad]


 call get_distance( v1 ,        v2, a)
 call get_distance( v1 , vector_ec, b)
 call get_distance( v2 , vector_ec, c)

 alpha = dacos( ( b*b + c*c - a*a ) / (2.*b*c) )  ![rad]

 arc_dist = alpha * 0.5*(b+c)

 return
 end subroutine get_arc_distance

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine nadir_angle_ecross(nsample, scan_speed, nadir_angle)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Estimate viewing nadir angle for electric cross tracking scanning system.  
!
! History:
! 10/2010  Toshi Matsui@NASA GSFC         : initial
!
! References:
!
!----------------------------------------------------------------------
 integer,intent(in) :: nsample             ! count of sampling per one scan
 real(sdsu_fpd),intent(in) :: scan_speed   ! scan angular velocity [rad/sec] 
 real(sdsu_fpd),intent(out) :: nadir_angle ! nadir angle for cross tracking scanning [rad]

 real(sdsu_fpd) :: dna  ! delta nadir angle per one sampling time [rad]
 real(sdsu_fpd) :: sna  ! starting nadir angle [rad]

!
!
!
 dna = dona_dpr * degrad   ! delta nadir angle per one sampling time [rad]
 sna = sensor%ona         ! starting nadir angle [rad]
 nadir_angle = dna * dble(nsample-1) - sna ! nadir angle [rad]

!
! special correction for matched ka band
!
 if( trim(sensor_name) == 'DPR_Ka' ) then 
   if(nsample >= 1 .and. nsample <= 25) then   
      sna = dona_dpr*degrad * 12.d0 
      nadir_angle = dna * dble(nsample-1) - sna ! nadir angle for Ka_MA [rad]
   elseif( nsample >= 26 .and. nsample <= 49 ) then
      sna = dona_dpr*degrad * 11.5d0 
      nadir_angle = dna * dble(nsample-26) - sna ! nadir angle for Ka_MA [rad]
   else
      print*,'MSG nadir_angle_ecross: nsample is too high',nsample
      call stop_sdsu(myrank,'') 
   endif

! print*,'nadir angle',nsample,nadir_angle/degrad

 endif

 return
 end subroutine nadir_angle_ecross

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine scan_vector(orb, re, delta, alpha , v_fov, beta)  
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Estimate the maximum angle (ganma_scan) of observation perpendicular to the nadir track 
! by an instrument on the satellite for a given beta, and evs
!
! History:
! 10/2010  Toshi Matsui@NASA GSFC         : Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 type( orbit_parameter ), intent(in) :: orb  ! keplarian orbit parameters
 real(sdsu_fpd),intent(in) :: re        ! distance between ellipsoid and earth's core [km]
 real(sdsu_fpd),intent(in) :: delta     ! conical sensor scanning angle [rad]
 real(sdsu_fpd),intent(in) :: alpha     ! cross-track sensor scanning angle [rad]
 real(sdsu_fpd),intent(out) :: v_fov(3) ! satellite position vector [km]
 real(sdsu_fpd),intent(out) :: beta      ! earth's incidence angle  [rad]

 real(sdsu_fpd) :: theta     ! earth's off-nadir angle [rad]
 real(sdsu_fpd) :: M_mod     ! modifined mean anomalies [rad]
 real(sdsu_fpd) :: f_mod     ! modified true anomalies [rad]
 real(sdsu_fpd) :: dxy,dx,dy,dz   ! 
 real(sdsu_fpd),dimension(3) :: v, v_omega, v_omega_i ! position vector  [km]


 select case ( trim(sensor%scan_typ) ) 
 case('ecross','cross')

   beta = asin(orb%a/re*sin(alpha))  ! earth's incidence angle [rad] 
   theta = beta - alpha     ! angle of scan point - earth's core - satellite position 
                            !<-- dpi = theta + (dpi-beta) + alpha 

   !
   ! initial scan position vector (the satellite orbital ellipse is assumed to lie in x-y plane)
   !
   v(1) =  re*dcos(theta)*dcos(orb%f)
   v(2) =  re*dcos(theta)*dsin(orb%f)
   v(3) =  -1.d0*re*dsin(theta)    ! -1 mean left-to-right scan

   !
   ! special correction for DPR matched ka band 
   ! (adding 12 sampling-time lag to match up with DPR Ku)
   !
   if( trim(sensor_name) == 'DPR_Ka' ) then  

      M_mod = orb%M + orb%M_dot  * ( 12.d0 * sensor%samp_dt)  !modified mean anomalies [rad]

      call solve_true_anomalie( orb%e, M_mod, f_mod )

      v(1) =  re*dcos(theta)*dcos(f_mod)
      v(2) =  re*dcos(theta)*dsin(f_mod)
      v(3) =  -1.d0*re*dsin(theta)   !-1 mean left-to-right scan

   endif

 
 case('conical')
  
   beta = asin(orb%a/re*sin(sensor%ona))  ! earth's incidence angle based on off-nadir angle [rad] 
   theta = beta - sensor%ona              ! dpi = earth_core_angle + (dpi-beta) + na

    dxy = re*dsin(theta)*dsin(delta)
    dx = dxy*dcos(orb%f+dpi/2.d0)
    dy = dxy*dsin(orb%f+dpi/2.d0)
    dz = re*dsin(theta)*dcos(delta)
   !
   ! initial scan position vector (the satellite orbital ellipse is assumed to lie in x-y plane)
   !
   v(1) =  re*dcos(theta)*dcos(orb%f) + dx
   v(2) =  re*dcos(theta)*dsin(orb%f) + dy
   v(3) =  1.d0*dz      ! +1 mean right-to-left scan

 case default ; call stop_sdsu(myrank, 'MSG scan_vector: There is no such sensor%scan_typ')
 end select


!
! Rotate v about the z axis throught the argument of perigee (omega)
!
 call rotation_matrix ( 'z', orb%omega, v       , v_omega        )

!
! Rotate v_omega about the x axis throught the inclination angle (i)
!
 call rotation_matrix ( 'x',     orb%i, v_omega , v_omega_i      )

!
! Rotate v_omega_i about the z axis through the right ascension of ascending node (ohm)
!
 call rotation_matrix ( 'z',   orb%ohm, v_omega_i, v_fov        )


 return
 end subroutine scan_vector

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

 subroutine right_ascension ( jd, ohm_earth )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Compute earth's rotation angle (right ascensioin) using Julian days and UT.  
! 
! History:
! 11/2010  Toshi Matsui@NASA GSFC : Initial
!
! References:
!
!  Brooks, D. R., 1977: An introduction to orbit dynamics and its application to satellite 
!       based earth monitoring missions. NASA Ref. Publication 1009 Washington, DC, 80 pp. 
!       [NTIS No. N78-12113]
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: jd         ! julian day (including UT component) ( from Jan 0.5 1900 )
 real(sdsu_fpd),intent(out) :: ohm_earth ! right ascension termr of earth's rotation [rad]

 real(sdsu_fpd) :: ut        ! julian days in UT component
 real(sdsu_fpd) :: jd_wout   ! julian day from Jan 0.5 1900
 real(sdsu_fpd) :: jc        ! Julian century from Jan 0.5 1900
 real(sdsu_fpd) :: alpha_g0  ! Right ascension of the Greenwich meridian at 12 midnight (0hr UT) [deg]
 real(sdsu_fpd) :: alpha     ! Right ascension of the Greenwich meridian for a given UT [deg]

!
! estimate UT 
!
 ut = (jd - 0.5d0) - dble(int(jd - 0.5d0)) ! universal time [day] (0.5 account for Julday 0.5)

!
! julian day without UT from Jan 0.5 1900
!
 jd_wout = jd - ut  ! [day]

!
! Julian century from Jan 0.5 1900
!
 jc = ( jd_wout - 2415020.d0 ) / 36525.d0   !eq2 in Books [1977]

!
! Right ascension (local sidereal time) of the Greenwich meridian at 12 midnight (0hr UT)
!
 alpha_g0 = 99.6909833d0 + 36000.7689d0*jc + 0.00038708d0*jc*jc  ![deg] eq1 in Books [1977]

 alpha = dmod( alpha_g0 + earth_rotate_1day * ut , 360.d0 ) !eq3 in Books [1977]

 ohm_earth = alpha*degrad  ! right ascension termr of earth's rotation [rad]

 return
 end subroutine right_ascension

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine grid_point(lat,lon,j_grid,i_grid)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Estimate grid point for given lat and lon cordinate.   
!
! History:
! 10/2010  Toshi Matsui@NASA GSFC         : initial
!
! References:
!
!----------------------------------------------------------------------

 real(sdsu_fpd),intent(in)  :: lat,lon
 integer,intent(out) :: j_grid,i_grid

 integer :: nlat,nlon
 real(sdsu_fpd),parameter :: rdels=1.d0/grid_space
 real(sdsu_fpd),parameter :: cyclei=360.d0/grid_space
 real(sdsu_fpd) :: lon_k, lat_k 
 real(sdsu_fpd),parameter:: xlon= 0.d0+grid_space*0.5d0   ! starting longitude [deg]
 real(sdsu_fpd),parameter:: xlat=90.d0-grid_space*0.5d0   ! starting latitude [deg]

 lon_k = 1.d0+ rdels*(lon /degrad - xlon)
 lat_k = 1.d0- rdels*(lat /degrad - xlat)

 if(lon_k.lt.xlon) lon_k=lon_k+cyclei

 if(lon_k.ge.cyclei+xlon) lon_k=lon_k-cyclei

 nlon = nint(lon_k) ! make it interger in order to apply 1:1 for the number of array. 
 nlat = nint(lat_k) ! 

 if( nlon == dint(360.d0/grid_space) ) nlon = 1
 if( nlat == dint(180.d0/grid_space) ) nlon = dint(180.d0/grid_space) - 1

 i_grid = nlon
 j_grid = dint(180.d0/grid_space) - nlat

 return
 end subroutine grid_point

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine scan_statistics
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Derive various scan orbit statistics.   
!
! History:
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer :: i,j,i_scan, j_scan
 integer :: im,ie,jm,jm_cnt, ierr
 integer :: n, num
 real(sdsu_fpd) :: distance, distance_max, distance_min ![km] 
 real(sdsu_fpd) :: sum_value
 real(sdsu_fpd) :: a,b,c  ![km]
 real(sdsu_fpd) :: alpha,beta,gam ![rad]
 real(sdsu_fpd) :: swath ![km]
 real(sdsu_fpd) :: arc_distance ![km]
 real(sdsu_fpd) :: crit_distance ![km]
 real(sdsu_fpd) :: vector_op(3)
 real(sdsu_fpd) :: scan_speed
 real(sdsu_fpd) :: mean_angle_fov  ![rad]
 real(sdsu_fpd) :: mean_angle_scan ![rad]
 real(sdsu_fps) :: lat_fov, lon_fov
 integer :: cnt_point(0:jmax_scan)
 integer :: jmaxloc(1:1)

 if(masterproc .and. verbous_scan ) then
     print*,''
     print*,'MSG scan_statistics: derive orbit-scan statistics'
 endif
 im = imax_scan / 2
 if(imax_scan == 1) im = 1

!
! distance between one scan to next scan
!
 sum_value = 0.  ; num = 0 ; distance_max = 0.
 do j_scan = 1, jmax_scan-1
    call get_arc_distance( scans(im,j_scan)%vector_fov, scans(im,j_scan+1)%vector_fov,&
                           distance ) 
    num = num + 1
    sum_value = sum_value + distance  ![km]  

    if(distance > distance_max) distance_max = distance

 enddo

 stat%mean_scan_dist = sum_value / dble(num)  !mean distance between one scan and next scan.  [km]
 stat%max_scan_dist  = distance_max

if(verbous_scan) then
#if MPI == 2
        if(masterproc) &  !for domain decomposition case, write by masterproc only.
#endif
 print*,'  Mean distance between scan line is ', sum_value / dble(num) ,'[km]'
endif

!
! distance between fov in a line
!
 if(im==1) then
  stat%mean_fov_dist = 0.
  stat%max_fov_dist  = 0. 
 else
  sum_value = 0.  ; num = 0  ; distance_max = 0.

  do j_scan = 1, jmax_scan
     do i_scan = 1, imax_scan - 1
     call get_arc_distance( scans(i_scan,j_scan)%vector_fov, scans(i_scan+1,j_scan)%vector_fov,&
                           distance )
     num = num + 1
     sum_value = sum_value + distance  ![km]  
     if(distance > distance_max) distance_max = distance
     enddo
  enddo
 
  stat%mean_fov_dist = sum_value / dble(num)  !mean distance between fov in line
  stat%max_fov_dist  = distance_max
 endif

if(verbous_scan) then
#if MPI == 2
        if(masterproc) &  !for domain decomposition case, write by masterproc only.
#endif
 print*,'  Mean distance between fov within one scan is ', sum_value / dble(num) ,'[km]'
endif

!
! swath width (on elliptical surface )
!
 if(im==1) then
  stat%mean_swath_width = 0.
 else
  sum_value = 0.  ; num = 0
  do j_scan = 1, jmax_scan
     ie = imax_scan
     if( trim(sensor_name) == 'DPR_Ka' ) ie = 25
     call get_arc_distance( scans(1,j_scan)%vector_fov , scans(ie,j_scan)%vector_fov, swath )
     num = num + 1
     sum_value = sum_value + swath ![km]  
  enddo
  stat%mean_swath_width = sum_value / dble(num)  ! mean swath width [km]
 endif

 if(verbous_scan) then
#if MPI == 2
    if(masterproc) &  !for domain decomposition case, write by masterproc only.
#endif
    print*,'  Mean swath width is ', sum_value / dble(num) ,'[km]'
 endif


!
! mean incidence angle 
!
 if(im==1) then

  sum_value = 0.  ; num = 0
  do j_scan = 1, jmax_scan
     num = num + 1
     sum_value = sum_value + scans(im,j_scan)%inci_ang  ![rad]  
  enddo
  stat%mean_inc_angle = sum_value / dble(num)  !mean incident angle [rad]

  stat%max_inc_angle = maxval( scans(im,:)%inci_ang )  !max incident angle [rad]

 else
  sum_value = 0.  ; num = 0
  do j_scan = 1, jmax_scan
     do i_scan = 1, imax_scan - 1
     num = num + 1
     sum_value = sum_value + scans(i_scan,j_scan)%inci_ang  ![rad]  
     enddo
  enddo

  stat%mean_inc_angle = sum_value / dble(num)  !mean incident angle [rad]

  stat%max_inc_angle = maxval( scans(:,:)%inci_ang )  !max incident angle [rad]

 endif

if(verbous_scan) then
#if MPI == 2
        if(masterproc) &  !for domain decomposition case, write by masterproc only.
#endif
 print*,'  Mean incident angle ', stat%mean_inc_angle/degrad ,'[deg]'

#if MPI == 2
        if(masterproc) &  !for domain decomposition case, write by masterproc only.
#endif
 print*,'  Max  incident angle ', stat%max_inc_angle/degrad ,'[deg]'

endif


!
! maximum/minimum distance between FOV and satellite (used to estimate range bin in radar)
!
 distance_max = 0. ; distance_min = 1.d+10
 do j_scan = 1, jmax_scan ; do i_scan = 1, imax_scan

    call get_distance( scans(i_scan,j_scan)%vector_fov , scans(i_scan,j_scan)%vector_sat, &
                       distance )

    if(distance > distance_max) distance_max = distance
    if(distance < distance_min) distance_min = distance


 enddo ; enddo

 stat%max_sat_fov_dist  = distance_max  ! [km]
 stat%min_sat_fov_dist  = distance_min  ! [km]


!
! Find out start and end index for j_scan index (js_scan and je_scan).
! Reason --> to reduce orbital output size (memory).
! This is needed to determin js_scan and je_scan to minimize do loop in subroutine scan_interpolate.
!
  js_scan = undefined_i2
  forward_loop: do j_scan = 1, jmax_scan 
     do i_scan = 1, imax_scan

        lat_fov = real(scans(i_scan,j_scan)%latlon_fov(1)/degrad) !fov latitude [deg]
        lon_fov = real(scans(i_scan,j_scan)%latlon_fov(2)/degrad) !fov longitude [deg]

        !
        ! within domain ?
        !
!        if( lat_fov < domain_stat%maxlat .and. lat_fov > domain_stat%minlat ) then !for GPM orbital testbed V1.1
!        if( lon_fov < domain_stat%maxlon .and. lon_fov > domain_stat%minlon ) then
        if( lat_fov < maxlat .and. lat_fov > minlat ) then
        if( lon_fov < maxlon .and. lon_fov > minlon ) then

           js_scan = j_scan
           exit forward_loop

        endif
        endif

     enddo
  enddo forward_loop
  !
  ! stop program here, if overpass is not close enough to the CRM domain.
  !
   if( js_scan == undefined_i2  ) then
      call stop_sdsu(myrank, 'MSG scan_statistics: This overpassing point will not overpass CRM domain.')
   endif

  je_scan = undefined_i2
  backward_loop: do j_scan = jmax_scan, 1, -1
     do i_scan = 1, imax_scan

        lat_fov = real(scans(i_scan,j_scan)%latlon_fov(1)/degrad) !fov latitude [deg]
        lon_fov = real(scans(i_scan,j_scan)%latlon_fov(2)/degrad) !fov longitude [deg]

        !
        ! within domain ?
        !
!        if( lat_fov < domain_stat%maxlat .and. lat_fov > domain_stat%minlat ) then  !for GPM orbital testbed V1.1
!        if( lon_fov < domain_stat%maxlon .and. lon_fov > domain_stat%minlon ) then
        if( lat_fov < maxlat .and. lat_fov > minlat ) then
        if( lon_fov < maxlon .and. lon_fov > minlon ) then


           je_scan = j_scan
           exit backward_loop

        endif
        endif

     enddo
  enddo backward_loop

 if(masterproc) print*,'within-domain j_scan are', js_scan, je_scan

!
! mean angle between one scan and next scan.
!
 sum_value = 0.  ; num = 0
 do j_scan = 1, jmax_scan-1

    call triangle_parameter_scan( scans(im,j_scan)%vector_sat, scans(im,j_scan)%vector_fov, &
                                  scans(im,j_scan+1)%vector_fov,alpha, beta, gam, a, b, c )
    num = num + 1
    sum_value = sum_value + alpha  ![rad]  
 enddo

 mean_angle_scan = sum_value / dble(num)    ![rad]

!
! mean angle between one fov sample and next fov sample
!
if(im==1) then
 stat%mean_scan_angle = 0.
else
 sum_value = 0.  ; num = 0
 mean_angle_fov = mean_angle_scan
 if( imax_scan > 3 ) then
   do j_scan = 1, jmax_scan-1

      call triangle_parameter_scan( scans(im,j_scan)%vector_sat, scans(im,j_scan)%vector_fov, &
                                    scans(im+1,j_scan)%vector_fov,alpha, beta, gam, a, b, c )
      num = num + 1
      sum_value = sum_value + alpha  ![rad]  
   enddo
   mean_angle_fov = sum_value / dble(num) ![rad]

 endif

  ! use smaller angle for output fov beamwidth 

  stat%mean_scan_angle =  min( mean_angle_scan , mean_angle_fov )   ![rad]
endif


!
! allocate other scan array
!
 if( .not. allocated(scans_L2) ) then
   allocate( scans_L2(imax_scan,js_scan:je_scan) , &
             stat=ierr )
   if (ierr /= 0) call stop_sdsu(myrank,'MSG allocate_scan: scans_L2 allocation error -> Terminate program. ')
 endif

 do n = 1, nparam_L2
     scans_L2%params(n) = undefined  !this initialization MUST be required for MPI_defined communication 
 enddo

 return
 end subroutine scan_statistics

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine scan_interpolate
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Interpolate satellite-fov vector upon CRM surface grid.  
!            This routine apply to any type of scanning sensors. 
!            output is surface_CRM(:,:)%###
!
! History:
! 10/2012  Toshi Matsui@NASA GSFC         : Add option for coarse-resolution mesoscale model. 
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer :: i,j   !loop for CRM domain
 integer :: i_scan, j_scan !loop for scan sample domain
 integer :: ierr, io
 integer :: il, jl
 integer :: ij_loc(1:2) 
 integer :: i_rank, irank_max
 real(sdsu_fpd),allocatable :: distance_grid2fov(:,:)
 real(sdsu_fpd) :: max_dt,max_ct,crit_distance !critical distance
 real(sdsu_fpd) :: distance,lat_test,lon_test, crit_collocate
 real(sdsu_fpd) :: sum_dist , sum_vec(1:3)
 real(sdsu_fpd) :: dgrid  !diagonal distance of CRM grid
 logical :: large_grid    ! 
!
! critical ccollocation
!
 crit_collocate = 0.5d0*sqrt( gridsize**2+gridsize**2  )
 large_grid = .false.
!
! simplify loop index name
!
 if(masterproc .and. verbous_scan ) then
   print*,''
   print*,'MSG scan_interpolate: interpolate scan sample to CRM grid'
 endif

!
! allocate surface_CRM parameters 
!
 if( .not. allocated(surface_CRM) ) then
   allocate( surface_CRM(mxgridx,mxgridy) , stat=ierr )
  if (ierr /= 0) call stop_sdsu(myrank,'MSG scan_interpolate: allocation error -> Terminate program. ')
 endif
 surface_CRM(:,:)%within_swath = .false.  !logical parameter CRM grid within sensor swath or not. 


!
! allocate grid-FOV distance parameter
!
 if( allocated(distance_grid2fov) ) deallocate(distance_grid2fov)
 allocate( distance_grid2fov( 1:imax_scan, js_scan:je_scan ) , stat=ierr )
 if (ierr /= 0) call stop_sdsu(myrank, 'MSG scan_interpolate: allocation error -> Terminate program. ')
 distance_grid2fov = 1.d5 !initialize as large value

!
! get critical searching distance (toshii : modify lator by better alogirhtm)
! modify this for different sensors
!
 irank_max = 4  !default searching 

 select case( trim(simulator_typ) )
 case('micro') 
   max_dt = dble( maxval(fov_dt_micro(1:mxfreq_micro)) )
   max_ct = dble( maxval(fov_ct_micro(1:mxfreq_micro)) )
   dgrid = sqrt( gridsize**2+gridsize**2  ) * 1.5  ! 1.2 x (diagonal distance of CRM grid)
   stat%grid_gt_fov = .false.
   if( dgrid >= max_dt .and. dgrid >= max_ct ) then
        max_dt = dgrid  !for large grid size application. 
!        irank_max = 1 ! have many missing data
        stat%grid_gt_fov = .true.
   endif
 case('radar')
   max_dt = dble( 2.*maxval(fov_dt_radar(1:mxfreq_radar)) )  !make sure 2x fov size
   max_ct = dble( 2.*maxval(fov_ct_radar(1:mxfreq_radar)) )
   dgrid = sqrt( gridsize**2+gridsize**2  ) * 1.5    ! 1.2 x (diagonal distance of CRM grid)
   stat%grid_gt_fov = .false.

  if( dgrid > max_dt .and. dgrid > max_ct ) then
       max_dt = dgrid  !for large grid size application. 
       irank_max = sensor%mx_sample * (INT(crit_collocate/stat%max_scan_dist) + 4) !
       stat%grid_gt_fov = .true.
       if(masterproc) print*,'MSG scan_interpolate: Coarse Resolution Interpolation'
  endif

 case('visir')

  ! visir sensor tends to have large sampling distance than FOV size. 
  ! so it changed as follow. 
   max_dt =   stat%max_fov_dist   
   max_ct =   stat%max_scan_dist
   stat%grid_gt_fov = .false.

 case('lidar')  
  ! LIDAR FOV is generally much smaller than CSRM. So use grid size criteria. 
   max_dt = sqrt( gridsize**2+gridsize**2  ) * 1.5
   max_ct = sqrt( gridsize**2+gridsize**2  ) * 1.5
   irank_max = sensor%mx_sample * (INT(crit_collocate/stat%max_scan_dist) + 4) !
   stat%grid_gt_fov = .true.
 case('broad')
   max_dt = dble( fov_dt_broad ) 
   max_ct = dble( fov_ct_broad ) 
   stat%grid_gt_fov = .false.
 end select

!
! critical distance 
!
 crit_distance = max(max_dt, max_ct)


!
! CRM doamin spatial 2D loop start here (actually memory loop here)
!
  do j = myj_start, myj_end ; do i = myi_start, myi_end

     call latlon2cartesian(surface(i,j)%lat*degrad , surface(i,j)%lon*degrad ,&
                           surface_CRM(i,j)%vector_grid )

     !
     ! scan sample loop to get arc distance betweeen CRM grid and FOV
     !
     do j_scan = js_scan , je_scan  ; do i_scan = 1, imax_scan
        call get_arc_distance( surface_CRM(i,j)%vector_grid, scans(i_scan,j_scan)%vector_fov,&
                               distance_grid2fov(i_scan,j_scan) ) 
     enddo ; enddo 

     !
     ! get the closest irank_max scan points from one CRM grid to each scan_fov locations
     !
     sum_dist = 0.  ; sum_vec(:) = 0.
     do i_rank = 1,irank_max
        ij_loc(1:2) = minloc( distance_grid2fov(1:imax_scan,js_scan:je_scan) )
        il = ij_loc(1) ; jl = js_scan+ij_loc(2)-1  !adjust it since minloc assum ubound is 1. 

        !
        ! sums up to get get mean satellite vector position
        !
        if(stat%grid_gt_fov) then !coarse mode 
           if ( distance_grid2fov(il,jl) < crit_collocate ) then
               sum_dist = sum_dist + distance_grid2fov(il,jl)   !sum of weight
               sum_vec(1:3)  = sum_vec(1:3) +  scans(il,jl)%vector_sat(1:3)*distance_grid2fov(il,jl) !weighting sum
           endif
        else
           if ( distance_grid2fov(il,jl) < crit_distance ) then
               sum_dist = sum_dist + distance_grid2fov(il,jl)   !sum of weight
               sum_vec(1:3)  = sum_vec(1:3) +  scans(il,jl)%vector_sat(1:3)*distance_grid2fov(il,jl) !weighting sum
           endif
        endif

        !
        ! center FOV and CRM grid are collocated or not 
        !
        if ( distance_grid2fov(il,jl) <= crit_collocate ) then
           scans(il,jl)%i_crm = i
           scans(il,jl)%j_crm = j
        endif

        ! then remove the closest points
        distance_grid2fov( il , jl ) = 1.d9  !remove out 

     enddo    

     !
     ! output to CRM grid --> surface_CRM(i,j)
     !
     if ( sum_dist == 0. ) then !bad CRM grid
         surface_CRM(i,j)%vector_sat(1:3) = undefined_scan
         surface_CRM(i,j)%latlon_sat(1:2) = undefined_scan
         surface_CRM(i,j)%latlon_grid(1:2)= undefined_scan
         surface_CRM(i,j)%a               = undefined_scan
         surface_CRM(i,j)%b               = undefined_scan
         surface_CRM(i,j)%c               = undefined_scan
         surface_CRM(i,j)%alpha           = undefined_scan
         surface_CRM(i,j)%beta            = undefined_scan
         surface_CRM(i,j)%gam             = undefined_scan
         surface_CRM(i,j)%within_swath    = .false.
         surface_CRM(i,j)%umu_local       = undefined
     else
         surface_CRM(i,j)%vector_sat(1:3) = sum_vec(1:3)/sum_dist  ! satellite vector [km]
         call cartesian2latlon( surface_CRM(i,j)%vector_sat, surface_CRM(i,j)%latlon_sat(1),&
                                surface_CRM(i,j)%latlon_sat(2) )   ! lat lon [rad]
         call cartesian2latlon( surface_CRM(i,j)%vector_grid, surface_CRM(i,j)%latlon_grid(1),&
                                surface_CRM(i,j)%latlon_grid(2) )   ! lat lon [rad]

         call triangle_parameter( surface_CRM(i,j)%vector_sat, surface_CRM(i,j)%vector_grid, &
                                  surface_CRM(i,j)%alpha , surface_CRM(i,j)%beta, surface_CRM(i,j)%gam,    &  
                                  surface_CRM(i,j)%a, surface_CRM(i,j)%b, surface_CRM(i,j)%c  ) !
         surface_CRM(i,j)%umu_local = real( dcos( dpi - surface_CRM(i,j)%beta ) )  ! cosine of local incident angle 
         surface_CRM(i,j)%within_swath    = .true.

         ! visir sensor scan is a bit tricky. 
         if( trim(simulator_typ) == 'visir' .and. surface_CRM(i,j)%alpha > sensor%ona ) then !off the swath
            surface_CRM(i,j)%vector_sat(1:3) = undefined_scan
            surface_CRM(i,j)%latlon_sat(1:2) = undefined_scan
            surface_CRM(i,j)%latlon_grid(1:2)= undefined_scan
            surface_CRM(i,j)%a               = undefined_scan
            surface_CRM(i,j)%b               = undefined_scan
            surface_CRM(i,j)%c               = undefined_scan
            surface_CRM(i,j)%alpha           = undefined_scan
            surface_CRM(i,j)%beta            = undefined_scan
            surface_CRM(i,j)%gam             = undefined_scan
            surface_CRM(i,j)%within_swath    = .false.
            surface_CRM(i,j)%umu_local       = undefined
         endif

     endif     

  enddo ; enddo  !end of CRM domain 2D loop

!toshii
!  io = 10
!  open(unit= io, file = 'surface.bin', access='direct', status='replace', recl=mxgridx*mxgridy*4 )
!  write(io,rec=1) ((surface_CRM(i,j)%umu_local,i=1,mxgridx),j=1,mxgridy)
!  close(io)



#if MPI == 2
if(numproc_tot>1) then
!
! broadcasting vector_grid to all threads
!
  call mpi_sdsu_communicate( 'TO_ALL', surface_CRM(:,:)%vector_grid(1)  )
  call mpi_sdsu_communicate( 'TO_ALL', surface_CRM(:,:)%vector_grid(2)  )
  call mpi_sdsu_communicate( 'TO_ALL', surface_CRM(:,:)%vector_grid(3)  )

  call mpi_sdsu_communicate( 'TO_MASTER', surface_CRM(:,:)%latlon_grid(1) )
  call mpi_sdsu_communicate( 'TO_MASTER', surface_CRM(:,:)%latlon_grid(2) )
  call mpi_sdsu_communicate( 'TO_MASTER', surface_CRM(:,:)%umu_local      )

endif
#endif


 return
 end subroutine scan_interpolate

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  triangle_parameter( v_sat, v_grid, alpha, beta, gam, a, b, c  )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For given satellite and surface grid position, this subroutine compute angles and length
!            of triangular originating from the earth core v=(0,0,0).   
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
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fpd), intent(in) :: v_sat(1:3)   ! vector of satellite [km]
 real(sdsu_fpd), intent(in) :: v_grid(1:3)  ! vector of bottom grid [km]
 real(sdsu_fpd), intent(out) :: alpha,beta,gam ![rad]
 real(sdsu_fpd), intent(out) :: a,b,c  ![km]


 call get_distance( v_grid , vector_ec, a)
 call get_distance( v_sat  , vector_ec, b)
 call get_distance( v_sat  ,    v_grid, c)

!
! check if this is on nadir
!
 if( (a+c) - b < 1.d-10 )  then

   alpha = 0.     * degrad
   beta  = 180.d0 * degrad
   gam   = 0.     * degrad

 else
 
   alpha = dacos( ( b*b + c*c - a*a ) / (2.*b*c) )  ![rad]
   beta  = dacos( ( a*a + c*c - b*b ) / (2.*a*c) )  ![rad]
   gam   = dacos( ( a*a + b*b - c*c ) / (2.*a*b) )  ![rad]
 endif

 return
 end subroutine  triangle_parameter

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  triangle_parameter_scan( v_sat, v_grid1, v_grid2, alpha, beta, gam, a, b, c  )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For given satellite and surface grid position, this subroutine compute angles and length
!            of triangular of three arbitrally vectors.    
!
!  length a = CRM grid 1 ~ CRM grid 2
!         b = satellite ~ CRM grid 2 
!         c = satellite  ~ CRM grid 1
!
!  angle  alpha for a
!         beta  for b
!         gam   for c
!  
! History:
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fpd), intent(in) :: v_sat(1:3)   ! vector of satellite [km]
 real(sdsu_fpd), intent(in) :: v_grid1(1:3), v_grid2(1:3)  ! vector of bottom grid [km]
 real(sdsu_fpd), intent(out) :: alpha,beta,gam ![rad]
 real(sdsu_fpd), intent(out) :: a,b,c  ![km]

 call get_distance( v_grid1 ,   v_grid2, a )
 call get_distance( v_sat  ,   v_grid2, b )
 call get_distance( v_sat  ,    v_grid1, c )

 alpha = dacos( ( b*b + c*c - a*a ) / (2.*b*c) )  ![rad]
 beta  = dacos( ( a*a + c*c - b*b ) / (2.*a*c) )  ![rad]
 gam   = dacos( ( a*a + b*b - c*c ) / (2.*a*b) )  ![rad]

 return
 end subroutine  triangle_parameter_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine driver_slant_path_micro
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  At each CRM surface grid, this subroutine first interpolate optical paorperties and 
!            temperature along computed grid-satellite vector (slant path) for both downwelling and
!            upwelling directions. Then, it drives microwave radiative transfer to compute brightness
!            temperature.  
!
! History:
! 03/2014  Toshi Matsui@NASA GSFC         : bug fixed for k_top and up/downwelling index.
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------

 integer :: i,j,k,nf,np                 ! loop indice
 integer :: iup,idown                   
 integer :: k_bot,k_top(2),kmax            ! vertical memory bounds
 logical :: good_slant,good_slant_down  ! logic if slant path exist or not.
 type (opt_micro_parameter) :: slant_opt1d ( mxlyr_refine,mxfreq,2 ) !optical properties in slant paths
 real(sdsu_fps) :: slant_t_air(mxlyr_refine,2) !layer slant-path air temperature [K]
 real(sdsu_fps) :: t_air_stag(0:mxlyr_refine,2) !staggered slant-path air tempereature [K]
 real(sdsu_fps) :: slant_t_skin          !skin temperature in slant path [K]
 real(sdsu_fps) :: slant_emiss(mxfreq,2) ! microwave emissivity in slant path [-]
 real(sdsu_fps) :: slant_ebar (mxfreq)   ! ebar in slant path [-]
 real(sdsu_fps),parameter :: tb_dummy = 0. !dummy Tb
 real(sdsu_fps) :: tb_down_surface       ! boundary conditions for upwelling path
 real(sdsu_fps) :: umu_local             ! cosine of local incident angle

 if( masterproc ) then
   print*,''
   print*,'MSG driver_slant_path_micro: extract optical properties and Tb along sensor slant path'
 endif

 if( .not. allocated(bot_slant) ) allocate(bot_slant(mxgridx,mxgridy))

!
! DOMAIN_2D_LOOP
!
 do j = myj_start, myj_end ; do i = myi_start, myi_end
    !
    ! if CRM grid is out of satellite sensor swath, skip process. . 
    !
    if( .not. surface_CRM(i,j)%within_swath ) then
        tb_out(i,j,:,:) = undefined
        cycle
    endif

    umu_local = surface_CRM(i,j)%umu_local  ! cosine of local incident angle 

    !
    ! interpolate optical properties along slant path.
    !
    call slant_path1d_micro( i,j, surface_CRM(i,j) , good_slant , good_slant_down , k_bot, k_top, slant_opt1d , &
                             slant_t_air, slant_t_skin, slant_emiss, slant_ebar , bot_slant(i,j) )

    !
    ! if CRM grid does not have good slant-path
    ! , skip the processes. 
    !
    if( .not. good_slant ) then
        tb_out(i,j,:,:) = undefined
        cycle
    endif !bad slant path

!    print*,'k_top',k_top(1:2),'k_bot',k_bot, 'good_slant',good_slant,'good_slant_down',good_slant_down

    !
    ! Drive micriwave radiative transfer model along one/two-way slant paths. 
    !
    POL_LOOP: do np = 1, 2  !vertical and horizontal loop
    FREQ_LOOP: do nf = 1, mxfreq  !frequency loop

       if    ( good_slant .and. .not.good_slant_down) then  !one-way slant path


          !
          ! upwelling path is used here.
          !
          kmax = k_top(1) - k_bot + 1  !upwelling path only
          t_air_stag(k_bot-1,1) = slant_t_air(k_bot,1)
          t_air_stag(k_top(1)  ,1) = slant_t_air(k_top(1),1)
          do k = k_bot , k_top(1)-1
             t_air_stag(k,1) = 0.5 * ( slant_t_air(k,1) + slant_t_air(k+1,1) )
          enddo

          call eddington_2way('downup', kmax , umu_local, ground_micro, slant_t_skin         ,  &
                              t_air_stag(k_bot-1:k_top(1),1), hgt_stag_micro(k_bot-1:k_top(1)),  &
                              slant_opt1d (k_bot:k_top(1),nf,1)%kext,  &
                              slant_opt1d (k_bot:k_top(1),nf,1)%salb,  &
                              slant_opt1d (k_bot:k_top(1),nf,1)%asym,  &
                              slant_emiss(nf,np), slant_ebar(nf),      &
                              tb_dummy, tb_out(i,j,nf,np)    )

        elseif( good_slant .and. good_slant_down) then  !two-way slant path

          ! 
          ! downwelling path first
          !
          idown = 2
          kmax = k_top(idown) - k_bot + 1  !downwelling path 

          ! compute stagger-level air temperature
          t_air_stag(k_bot-1,idown) = slant_t_air(k_bot,idown)
          t_air_stag(k_top(idown),idown) = slant_t_air(k_top(idown),idown)
          do k = k_bot , k_top(idown)-1
             t_air_stag(k,idown) = 0.5 * ( slant_t_air(k,idown) + slant_t_air(k+1,idown) )
          enddo

          call eddington_2way('down', kmax , umu_local, ground_micro, slant_t_skin ,&
                              t_air_stag(k_bot-1:k_top(idown),idown),  &
                              hgt_stag_micro(k_bot-1:k_top(idown)),  &
                              slant_opt1d (k_bot:k_top(idown),nf,idown)%kext,  &
                              slant_opt1d (k_bot:k_top(idown),nf,idown)%salb,  &
                              slant_opt1d (k_bot:k_top(idown),nf,idown)%asym,  &
                              slant_emiss(nf,np), slant_ebar(nf),      &
                              tb_dummy, tb_down_surface   )

          ! 
          ! upwelling path
          !
          iup = 1
          kmax = k_top(iup) - k_bot + 1  !upwelling path
          t_air_stag(k_bot-1,iup) = slant_t_air(k_bot,iup)
          t_air_stag(k_top(iup),iup) = slant_t_air(k_top(iup),iup)
          do k = k_bot , k_top(iup)-1
             t_air_stag(k,iup) = 0.5 * ( slant_t_air(k,iup) + slant_t_air(k+1,iup) )
          enddo

          call eddington_2way('up', kmax , umu_local, ground_micro, slant_t_skin  ,  &
                              t_air_stag(k_bot-1:k_top(iup),iup), &
                              hgt_stag_micro(k_bot-1:k_top(iup)),  &
                              slant_opt1d (k_bot:k_top(iup),nf,iup)%kext,  &
                              slant_opt1d (k_bot:k_top(iup),nf,iup)%salb,  &
                              slant_opt1d (k_bot:k_top(iup),nf,iup)%asym,  &
                              slant_emiss(nf,np), slant_ebar(nf),      &
                              tb_down_surface, tb_out(i,j,nf,np)    )

        endif

!        print*,'TB OUTPUT',i,j,nf,np,tb_out(i,j,nf,np)

    enddo FREQ_LOOP
    enddo POL_LOOP

 enddo ; enddo


#if MPI == 2
if(numproc_tot>1) then
!
! let's send each tile into the all processor. 
!
   call mpi_sdsu_communicate( 'TO_ALL', tb_out )
   call mpi_sdsu_communicate( 'TO_ALL', bot_slant(:,:)%ir )
   call mpi_sdsu_communicate( 'TO_ALL', bot_slant(:,:)%jr )

endif
#endif



 return
 end subroutine driver_slant_path_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine driver_slant_path_radar
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  At each CRM surface grid, this subroutine first interpolate optical paorperties and 
!            temperature along computed grid-satellite vector (slant path) for both downwelling and
!            upwelling directions. Then, it compuate radar reflectivities.
!
! History:
! 02/2012  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer :: i,j,k,nf,np                 ! loop indice
 integer :: k_bot,k_top,kmax            ! vertical memory bounds
 logical :: good_slant,good_slant_down  ! logic if slant path exist or not.
 type (opt_micro_parameter) :: slant_opt1d ( mxlyr_refine,mxfreq) !optical properties in slant paths
 type (particle_retrieval) :: L2_radar_slant1d ( mxlyr_refine ) !optical properties in slant paths
 real(sdsu_fpd) :: dist_from_sat( mxlyr_refine )
 real(sdsu_fps) :: umu_local             ! cosine of local incident angle

! output from interp_radar_range 
 type (opt_micro_parameter)  ,allocatable :: range_opt1d  ( :, : ) ! optical properties along radar range
 type (particle_retrieval)   ,allocatable :: range_L2     ( : )    ! L2 parameter along radar range. 
 type (radar_range_parameter),allocatable :: range_info1d ( : )   ! atmos parameter

! output from radar return
 type (radar_output_parameter) ,allocatable :: radar_out1d( : )

 integer :: m
 integer :: ierr,io

 real(sdsu_fpd) :: range_max, range_min  !maximum and minimum range [km]
 real(sdsu_fpd) :: freqhzd, epsreal, epsimag  ! dielectric constant parameters
 real(sdsu_fps) :: kradar     ! radar constant |K^2| 
 real(sdsu_fps) :: lambda_cm  ! wavelength [cm]
 real(sdsu_fps) :: loss ! 
 real(sdsu_fps) :: pt ! peak transmitt power [W]
 real(sdsu_fps) :: binres  ! bin resolution [m]
 real(sdsu_fps) :: sigma0  ! surface return factor [dB]
 logical,parameter :: estimate_power = .false.
 real(sdsu_fps),allocatable :: temp3d(:,:,:)


 if( masterproc ) then
   print*,''
   print*,'MSG driver_slant_path_radar: extract optical properties and dBZ along sensor slant path'
 endif

 if( .not. allocated(bot_slant) ) allocate(bot_slant(mxgridx,mxgridy))

!
! Prepration ( define basic radar range distance )
!    Based on the maximum satellite-FOV distance (off nadir 17deg),
!    define maximum range distance from satellite to targe point adjusted by
!    additional range distance below ground.
!
 range_max = stat%max_sat_fov_dist - minhgt_radar  !maximum radar range [km]
 range_min = stat%min_sat_fov_dist - mxhgt_radar   !minimum radar range  [km]

!
! define absolute maximum number of radar range (irange_max)
! toshii --> this can be reduced from ellipsoidal to specific level.(can reduce memory size 70%). 
!
 irange_max =  int( (range_max - range_min) / range_radar )  !define the number of range

 if( masterproc ) then
   print*,'Radar range number (absolute) is', irange_max
 endif

!
! allocate 1D output parameters
!
 if( allocated(range_opt1d) ) deallocate(range_opt1d)
 allocate( range_opt1d  ( irange_max, mxfreq ) ) ! optical properties along radar range

 if( allocated(range_L2)    ) deallocate(range_L2)
 allocate( range_L2     ( irange_max )         ) ! L2 parameter along radar range. 

 if( allocated(range_info1d)) deallocate(range_info1d)
 allocate( range_info1d ( irange_max )         ) ! radar range parameter

 if( allocated(radar_out1d) ) deallocate(radar_out1d)
 allocate( radar_out1d ( irange_max )        )  ! radar power outoupt


!
! allocate/initialize 3D output parameter
!
 if( allocated(radar_bundle3d) ) deallocate(radar_bundle3d)
 allocate( radar_bundle3d(ibs:ibe,jbs:jbe,irange_max) , &   !single frequency only to save memory
           stat=ierr  )

 do k = 1,irange_max ; do j = jbs,jbe ; do i = ibs,ibe
    call undefine_radar_bundle( radar_bundle3d(i,j,k) )
 enddo ; enddo ; enddo


!
! DOMAIN_2D_LOOP
!
 do j = myj_start, myj_end ; do i = myi_start, myi_end


    !
    ! if CRM grid is out of satellite sensor swath, skip process. . 
    !
    if( .not. surface_CRM(i,j)%within_swath  ) then
        cycle
    endif

    umu_local = surface_CRM(i,j)%umu_local  ! cosine of local incident angle 

    !
    ! interpolate optical properties along slant path.
    !
    call slant_path1d_radar( i,j, surface_CRM(i,j), good_slant, k_bot, slant_opt1d , &
                             L2_radar_slant1d , dist_from_sat, bot_slant(i,j) )


    !
    ! if CRM grid does not have good slant-path
    ! , skip the processes. 
    !
    if( .not. good_slant ) then
        cycle
    endif !bad slant path

    !
    ! interpolate slant path into radar range bin
    !
    call interp_radar_range( dble(umu_local), irange_max, dble(range_radar), range_min, &
                             slant_opt1d, L2_radar_slant1d , dist_from_sat, &
                             range_opt1d, range_L2, range_info1d  )
!    print*,'L2_radar_slant1d(:)%t_air=', L2_radar_slant1d(:)%t_air


    FREQ_LOOP: do nf = 1, mxfreq  !frequency loop

       freqhzd=dble(freq_radar(nf)*1.e+9)
       call watoptic(freqhzd, dble(273.16e0), dble(0.e0), epsreal, epsimag)
       kradar = ((epsreal-1.e0)**2 + epsimag**2)/ ((epsreal+2.e0)**2 + epsimag**2)
       if( k2(nf) /= undefined ) kradar = k2(nf)  !radar constant


       if( estimate_power ) then

         lambda_cm = 2.997925e2/freq_radar(nf) * 1.e-1  !wavelength [cm]
         loss = 10.**(-0./10.)   ! input 0dB
         pt = 1012.000
         binres = 125.
         sigma0 = 10.**(10./10.)  !replace this from sima0 map and use  bot_slant(i,j)
        !
        ! Estimate radar power and equivalent reflectivity factor in slant paths. 
        !
!      print*,'call return power',umu_local
!        call return_radar_power( irange_max, real(range_info1d%distance), range_opt1d(:,nf)%kext, &
!                                range_opt1d(:,nf)%sback, &
!                                lambda_cm, kradar, loss, pt, &
!                                binres, sigma0, radar_power1d)

      else

          call radar_return_scan( irange_max, real(range_info1d%dhgt)*1.e-3 , range_opt1d(:,nf)%kext , &
                                  range_opt1d(:,nf)%sback ,  &
                                  2.997925e2/freq_radar(nf), kradar, radar_out1d%Zt, radar_out1d%Zm, radar_out1d%tau )

          !
          ! transfer parameter to radar_bundle3d
          !
          radar_bundle3d(i,j,:)%Zt = radar_out1d(:)%Zt  !non-attenuating equiv reflectivity factor [Z]
          radar_bundle3d(i,j,:)%Zm = radar_out1d(:)%Zm  !non-attenuating equiv reflectivity factor [Z]
          radar_bundle3d(i,j,:)%tau = radar_out1d(:)%tau !path-integrated optical depth [-]
          radar_bundle3d(i,j,:)%iwc       = range_L2(:)%iwc       
          radar_bundle3d(i,j,:)%lwc       = range_L2(:)%lwc       
          radar_bundle3d(i,j,:)%rwc       = range_L2(:)%rwc      
          radar_bundle3d(i,j,:)%rain_rate = range_L2(:)%rain_rate
#ifdef FULL_L2_RADAR
          radar_bundle3d(i,j,:)%press     = range_L2(:)%press
          radar_bundle3d(i,j,:)%t_air     = range_L2(:)%t_air
          radar_bundle3d(i,j,:)%sh        = range_L2(:)%sh
          radar_bundle3d(i,j,:)%w         = range_L2(:)%w
          radar_bundle3d(i,j,:)%lfrac     = range_L2(:)%lfrac     
          radar_bundle3d(i,j,:)%bulk_den  = range_L2(:)%bulk_den  
          do m = 0,6
             radar_bundle3d(i,j,:)%mmt(m)  = range_L2(:)%mmt(m)
          enddo
#endif

      endif

     !
     ! range smoothing all bundled parameters
     !
     call smooth_radar_bundle( radar_bundle3d(i,j,:) )

    enddo FREQ_LOOP

 enddo ; enddo

#if MPI == 2
if(numproc_tot>1) then
!
! let's send each tile into the all processor. 
!
  allocate(temp3d(ibs:ibe,jbs:jbe,1:irange_max))

  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%Zt  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%Zm  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%tau  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%iwc )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%lwc )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%rwc )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%rain_rate)

#ifdef FULL_L2_RADAR
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%press )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%t_air )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%sh  )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%w   )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%lfrac    )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%bulk_den )
  do m = 0,6
     call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_bundle3d(ibs:ibe,jbs:jbe,:)%mmt(m) )
  enddo
#endif


  call mpi_sdsu_communicate( 'TO_ALL', bot_slant(:,:)%ir )
  call mpi_sdsu_communicate( 'TO_ALL', bot_slant(:,:)%jr )

endif
#endif

 return
 end subroutine driver_slant_path_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine interp_radar_range( umu, irange_max, range_resolution, range_min,  &
                                slant_opt1d, L2_radar_slant1d, dist_from_sat,  &
                                range_opt1d, range_L2, range_info1d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: interpolate CRM slant-path parameters into radar-range resolution.
!
! History:
! 11/2012  Toshi Matsui@NASA GSFC   : Add nearest neiboubor and linear interpolation options 
! 03/2012  Toshi Matsui@NASA GSFC   : Initial. 
!
! References:
!-------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: umu    !cosine of incident angle
 integer,intent(in) :: irange_max    !maximum range number
 real(sdsu_fpd),intent(in) :: range_resolution   ! [km]
 real(sdsu_fpd),intent(in) :: range_min          ! [km]
 type (opt_micro_parameter),intent(in) :: slant_opt1d      (mxlyr_refine, mxfreq)  !optical properties in slant paths
 type (particle_retrieval) ,intent(in) :: L2_radar_slant1d (mxlyr_refine) ! L2 parameter in slant paths
 real(sdsu_fpd)            ,intent(in) :: dist_from_sat    (mxlyr_refine)
 type (opt_micro_parameter),intent(out) :: range_opt1d  (irange_max, mxfreq)
 type (particle_retrieval),intent(out)  :: range_L2     (irange_max) ! optical properties in slant paths
 type ( radar_range_parameter ), dimension(irange_max), intent(inout) :: range_info1d ! atmos parameter

 real(sdsu_fpd) :: range_km
 real(sdsu_fpd) :: dist1, dist2  ![km]
 integer :: m
 integer :: irange
 integer :: islant, ifound
 character(len=20),parameter :: method = 'INTERP'  !'NN' is Nearest Neighbor 
                                                   !'INTERP' is Linear Interpolation

!
! Find nearest slant path parameter, and fill radar range 1D array
!
 do irange = 1, irange_max  ! range integration loop  1=top  ,irange_max=bottom

    range_info1d(irange)%dhgt = range_resolution * 1.d+3  !range length [m]

    range_info1d(irange)%distance = ( range_min + dble(irange-1) * range_resolution )*1.d+3  ![m]

    range_km = range_info1d(irange)%distance * 1.d-3 !range distance in [km]

    range_info1d(irange)%hgt = (dist_from_sat(1) - range_km ) * umu *1.d+3 !above ellipsoid height nomral to ellipsoid [m]

    ifound = undefined_i2

    if( range_km > dist_from_sat(1) ) then
!       print*,' below ellipsoid ',range_km, dist_from_sat(1)
       
    elseif( range_km < dist_from_sat(mxlyr_refine) ) then
!       print*,' above CRM top', range_km , dist_from_sat(mxlyr_refine)
 
    else ! find nearest

       slant_loop: do islant = 1, mxlyr_refine-1

          if( range_km >= dist_from_sat(islant+1) .and. &
              range_km <= dist_from_sat(islant)          ) then

              dist1 = dabs( range_km - dist_from_sat(islant+1) )  ![km]
              dist2 = dabs( range_km - dist_from_sat(islant)   )  ![km] 
              select case(trim(method))
              case('NN')
                if( dist1 < dist2 ) ifound = islant+1  
                if( dist2 < dist1 ) ifound = islant
              case('INTERP')
                ifound = islant
              end select       

              exit slant_loop
           endif

       enddo slant_loop

    endif


    if( ifound == undefined_i2 ) then

       range_opt1d(irange,:)%kext  = undefined
       range_opt1d(irange,:)%salb  = undefined
       range_opt1d(irange,:)%asym  = undefined
       range_opt1d(irange,:)%sback = undefined
       range_L2(irange)%press     = undefined
       range_L2(irange)%t_air     = undefined
       range_L2(irange)%sh        = undefined
       range_L2(irange)%w         = undefined
       range_L2(irange)%iwc       = undefined
       range_L2(irange)%lwc       = undefined
       range_L2(irange)%rwc       = undefined
       range_L2(irange)%lfrac     = undefined
       range_L2(irange)%bulk_den  = undefined
       range_L2(irange)%rain_rate = undefined
       range_L2(irange)%mmt(0:6)  = undefined

    else
       select case(trim(method))
       case('NN')
       range_opt1d(irange,:)%kext  = slant_opt1d(ifound,:)%kext
       range_opt1d(irange,:)%salb  = slant_opt1d(ifound,:)%salb
       range_opt1d(irange,:)%asym  = slant_opt1d(ifound,:)%asym
       range_opt1d(irange,:)%sback = slant_opt1d(ifound,:)%sback

       range_L2(irange)%press     = L2_radar_slant1d(ifound)%press
       range_L2(irange)%t_air     = L2_radar_slant1d(ifound)%t_air
       range_L2(irange)%sh        = L2_radar_slant1d(ifound)%sh
       range_L2(irange)%w         = L2_radar_slant1d(ifound)%w
       range_L2(irange)%iwc       = L2_radar_slant1d(ifound)%iwc
       range_L2(irange)%lwc       = L2_radar_slant1d(ifound)%lwc
       range_L2(irange)%rwc       = L2_radar_slant1d(ifound)%rwc
       range_L2(irange)%lfrac     = L2_radar_slant1d(ifound)%lfrac
       range_L2(irange)%bulk_den  = L2_radar_slant1d(ifound)%bulk_den
       range_L2(irange)%rain_rate = L2_radar_slant1d(ifound)%rain_rate
       range_L2(irange)%mmt(0:6)  = L2_radar_slant1d(ifound)%mmt(0:6) 

       case('INTERP')
       if( L2_radar_slant1d(ifound  )%press == undefined .or. &
           L2_radar_slant1d(ifound+1)%press == undefined ) then
       range_opt1d(irange,:)%kext  = undefined
       range_opt1d(irange,:)%salb  = undefined
       range_opt1d(irange,:)%asym  = undefined
       range_opt1d(irange,:)%sback = undefined
       range_L2(irange)%press     = undefined
       range_L2(irange)%t_air     = undefined
       range_L2(irange)%sh        = undefined
       range_L2(irange)%w         = undefined
       range_L2(irange)%iwc       = undefined
       range_L2(irange)%lwc       = undefined
       range_L2(irange)%rwc       = undefined
       range_L2(irange)%lfrac     = undefined
       range_L2(irange)%bulk_den  = undefined
       range_L2(irange)%rain_rate = undefined
       range_L2(irange)%mmt(0:6)  = undefined

       else
       range_opt1d(irange,:)%kext  = (dist1*slant_opt1d(ifound,:)%kext+&
                                         dist2*slant_opt1d(ifound+1,:)%kext)&
                                        / (dist1+dist2)
       range_opt1d(irange,:)%salb  = (dist1*slant_opt1d(ifound,:)%salb+&                                         
                                         dist2*slant_opt1d(ifound+1,:)%salb)&
                                        / (dist1+dist2)
       range_opt1d(irange,:)%asym  = (dist1*slant_opt1d(ifound,:)%asym+&
                                         dist2*slant_opt1d(ifound+1,:)%asym)&
                                        / (dist1+dist2)
       range_opt1d(irange,:)%sback = (dist1*slant_opt1d(ifound,:)%sback+&
                                         dist2*slant_opt1d(ifound+1,:)%sback)&
                                        / (dist1+dist2)

       range_L2(irange)%press     = (dist1*L2_radar_slant1d(ifound)%press+&
                                     dist2*L2_radar_slant1d(ifound+1)%press)&
                                        / (dist1+dist2)
       range_L2(irange)%t_air     = (dist1*L2_radar_slant1d(ifound)%t_air+&
                                     dist2*L2_radar_slant1d(ifound+1)%t_air)&
                                        / (dist1+dist2)
       range_L2(irange)%sh        = (dist1*L2_radar_slant1d(ifound)%sh+&
                                     dist2*L2_radar_slant1d(ifound+1)%sh)&
                                        / (dist1+dist2)
       range_L2(irange)%w         = (dist1*L2_radar_slant1d(ifound)%w+&
                                     dist2*L2_radar_slant1d(ifound+1)%w)&
                                        / (dist1+dist2)
       range_L2(irange)%iwc       = (dist1*L2_radar_slant1d(ifound)%iwc+&
                                     dist2*L2_radar_slant1d(ifound+1)%iwc)&
                                        / (dist1+dist2)
       range_L2(irange)%lwc       = (dist1*L2_radar_slant1d(ifound)%lwc+&
                                     dist2*L2_radar_slant1d(ifound+1)%lwc)&
                                        / (dist1+dist2)
       range_L2(irange)%rwc       = (dist1*L2_radar_slant1d(ifound)%rwc+&       
                                     dist2*L2_radar_slant1d(ifound+1)%rwc)&
                                        / (dist1+dist2)
       range_L2(irange)%lfrac     = (dist1*L2_radar_slant1d(ifound)%lfrac+&
                                     dist2*L2_radar_slant1d(ifound+1)%lfrac)&
                                        / (dist1+dist2)       
       range_L2(irange)%bulk_den  = (dist1*L2_radar_slant1d(ifound)%bulk_den+&
                                     dist2*L2_radar_slant1d(ifound+1)%bulk_den)&
                                     / (dist1+dist2)
       range_L2(irange)%rain_rate  = (dist1*L2_radar_slant1d(ifound)%rain_rate+& 
                                      dist2*L2_radar_slant1d(ifound+1)%rain_rate)&
                                     / (dist1+dist2)
       do m = 0, 6
          range_L2(irange)%mmt(m) = (dist1*L2_radar_slant1d(ifound)%mmt(m)+&
                                     dist2*L2_radar_slant1d(ifound+1)%mmt(m) )&
                                     / (dist1+dist2)
       enddo
       endif

       end select

    endif

 enddo

 return
 end subroutine interp_radar_range

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine slant_path1d_micro(i_ref, j_ref, crm , good_slant, good_slant_down, k_bot, k_top, slant_opt1d , &
                         slant_t_air, slant_t_skin, slant_emiss, slant_ebar , bot)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Find optical paorperties, temperature, and surface properties 
!            along the computed grid-satellite vector (slant path) for both downwelling and
!            upwelling directions. 
!
! History:
! 04/2014  Toshi Matsui@NASA GSFC : modification to include upper stratosphere layers.
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------

 integer,intent(in) :: i_ref, j_ref ! i,j point at reference ellisoid
 type(surface_CRM_parameter),intent(in) :: crm 
 logical,intent(out) :: good_slant
 logical,intent(out) :: good_slant_down

 integer,intent(out) :: k_bot, k_top(1:2) !k_top has two value for upwelling and downwelling.
 type (opt_micro_parameter),intent(inout) :: slant_opt1d(mxlyr_refine,mxfreq,2)

 integer :: slant_terrain    ! slant terrain index (0-air  1-terrain)
 real(sdsu_fps) :: slant_t_air(mxlyr_refine,2)
 real(sdsu_fps),intent(out) :: slant_t_skin
 real(sdsu_fps),intent(out) :: slant_emiss(mxfreq,2)
 real(sdsu_fps),intent(out) :: slant_ebar (mxfreq)
 type( location_slant ),intent(inout) :: bot
 
 integer :: k  !vertical loop
 integer :: i_found, j_found
 integer :: i_bot,j_bot
 real(sdsu_fpd) ::  vector_inc(3)  !incremental vector [km]
 real(sdsu_fpd) ::  vector_bot(3)  !bottom vector [km]
 real(sdsu_fpd) ::  vector_mirror(3) !vector of mirror satellite position [km]
 type(surface_CRM_parameter) :: mirror !mirror position parameters

 logical :: logic_found_grid
 logical :: it_was_terrain
 integer ::  kmax_up,  kmax_down
 logical,parameter :: comment = .false.

!
! initialize bottom and found grid index 
!
 i_found = i_ref 
 j_found = j_ref
 i_bot   = i_found
 j_bot   = j_found
 k_bot   = 1    
 k_top(1:2) = undefined_i2

 bot%ir = real(i_bot)  !bottom i index in real format
 bot%jr = real(j_bot)  !bottom j index in real format

 good_slant      = .true.  ! initialize logical value
 good_slant_down = .true.  ! initialize logical value
 it_was_terrain = .true.   ! start from terrain

 slant_t_air = undefined


!
! initialize surface parameter as reference point.
!
 slant_t_skin      = t_skin_micro(i_bot,j_bot    )
 slant_emiss(:,:)  = emis_micro  (i_bot,j_bot,:,:)
 slant_ebar(:)     = ebar_micro  (i_bot,j_bot,:  )

!
! ---------derive upwelling-slant (from bottom to top) properties -----------
!
 UP_LOOP: do k = 1, mxlyr_refine

    slant_terrain = 0 !initially assume atmosphere ( 0-air  1-terrain )
    !
    ! get vector increment of slant profile at every range_micro height
    !
    call vector_increment( range_micro,  k, crm, vector_inc )

    !
    ! for k==1 or previous k is terrain, update bottom vector.
    ! bottom vector is defined in the middle point of bottom atmospheric layer. 
    !
    if(k==1 .or. slant_terrain == 1) then 
       vector_bot = vector_inc  ! [km]
    endif

    !
    ! find out if increment vector is within the CRM grid points or not. 
    !
    call find_nearest_grid( k, i_found, j_found, vector_inc , logic_found_grid, &
                            slant_opt1d(k,:,1) , slant_t_air(k,1), slant_terrain )
    if( .not. logic_found_grid  ) then
        good_slant = .false.
        return 
    endif

    !
    ! if incremented point is terrain, update surafce position
    ! (make sure previous grid is not atmosphere)
    !
    if( it_was_terrain .and. slant_terrain == 1 ) then 
        k_bot             = k+1
        i_bot             = i_found
        j_bot             = j_found
        slant_t_skin      = t_skin_micro(i_bot,j_bot    )
        slant_emiss(:,:)  = emis_micro  (i_bot,j_bot,:,:)
        slant_ebar (:)    = ebar_micro  (i_bot,j_bot,:  )
        bot%ir = real(i_bot)
        bot%jr = real(j_bot)
        it_was_terrain = .true.
        cycle
    endif

    it_was_terrain = .false. 

 enddo UP_LOOP

!  print*,'slant_t_air(:,1)', slant_t_air(:,1)
 
 if( k_bot < 1 .or. k_bot > mxlyr_refine ) then
        if(comment) print*, 'MSG slant_path1d_micro: strange k_bot value =',k_bot
        good_slant = .false.
        return
 endif

!
! find out top layer height for upwelling path
!
 find_ktop: do k = mxlyr_refine, k_bot, -1
    if( slant_t_air(k,1) /= undefined ) then
       k_top(1) = k
       exit find_ktop
    endif
 enddo find_ktop

! strange k_top value 
 if( k_top(1) < 1 .or. k_top(1) > mxlyr_refine )  then
     if(comment)print*,'MSG slant_path1d_micro: strange k_top(1) value =',k_top(1)
     good_slant = .false.
     return
 endif 


!
! if there is terrain within layer, also .false.
!
 do k = k_bot, k_top(1)
    if(  slant_t_air(k,1) == undefined ) then
     good_slant = .false.
     return
    endif
 enddo

 
!
! Find out mirror satellite position based on v_bot and v_sat
!
 call find_mirror_position( crm%vector_sat, vector_bot, mirror)


!
! initial location of search location is bottom i,j
! 
 i_found = i_bot
 j_found = j_bot
 it_was_terrain = .true.   ! start from terrain

!
! ---------derive downwelling-slant (from bottom to top) properties -----------
!
 DOWN_LOOP: do k = k_bot, mxlyr_refine
    slant_terrain = 0 !assume atmosphere
    !
    ! get vector increment of slant profile
    !
    call vector_increment( range_micro, k, mirror, vector_inc)

    !
    ! find out if increment vector is within the CRM grid points or not. 
    !
    call find_nearest_grid( k, i_found, j_found, vector_inc , logic_found_grid, &
                            slant_opt1d(k,:,2) , slant_t_air(k,2), slant_terrain )

    if( .not. logic_found_grid  ) then
        good_slant_down = .false.
        if(comment)print*,'MSG slant_path1d_micro: downwelling path hits by bounds'
        return
    endif

!    if( slant_terrain == 1 ) then
!        good_slant_down = .false.
!        it_was_terrain = .true.   ! start from terrain
!        return
!    endif


 enddo DOWN_LOOP


!
! find out top layer height for downwelling path
!
 find_ktop2: do k = mxlyr_refine, k_bot, -1
    if( slant_t_air(k,2) /= undefined ) then
       k_top(2) = k
       if(comment)print*,'find k_top(2)',k,slant_t_air(k,2)
       exit find_ktop2
    endif
 enddo find_ktop2

 if( k_top(2) < 1 .or. k_top(2) > mxlyr_refine )  then
     if(comment)print*,'MSG slant_path1d_micro: strange k_top(2) value =',k_top(2)
     good_slant_down = .false.
     return
 endif 

!
! if downwelling layers are too shorter than upwelling layer, 
!  it may be stackeck by real steep terrain. 
!
 kmax_up   = k_top(1) - k_bot
 kmax_down = k_top(2) - k_bot
 if( kmax_down <  kmax_up/2 ) then
     if(comment) print*,'MSG slant_path1d_micro: downwelling layers are too shorter than upwelling layer'
     good_slant_down = .false.
     return
 endif

!
! if there is terrain within layer, also .false.
!
 do k = k_bot, k_top(2) 
    if(  slant_t_air(k,2) == undefined ) then
     if(comment) print*, 'MSG slant_path1d_micro: there is terrain within layer'
     good_slant_down = .false.
     return
    endif
 enddo

 return
 end subroutine slant_path1d_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine slant_path1d_radar(i_ref, j_ref, crm , good_slant, k_bot, slant_opt1d , &
                               L2_radar_slant1d , dist_from_sat, bot)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Find optical paorperties, temperature, and surface properties 
!            along the computed grid-satellite vector (slant path) for both downwelling and
!            upwelling directions. 
!
! History:
! 03/2012  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer,intent(in) :: i_ref, j_ref ! i,j point at reference ellisoid
 type(surface_CRM_parameter),intent(in) :: crm 
 logical,intent(out) :: good_slant
 integer,intent(out) :: k_bot
 type (opt_micro_parameter),intent(inout) :: slant_opt1d      ( mxlyr_refine, mxfreq )
 type (particle_retrieval),intent(inout)  :: L2_radar_slant1d ( mxlyr_refine ) ! optical properties in slant paths
 real(sdsu_fpd),intent(inout) :: dist_from_sat(mxlyr_refine) !distance from satellite [km]

 type( location_slant ),intent(inout) :: bot

 integer :: slant_terrain    ! slant terrain index (0-air  1-terrain)
 integer :: m,k                    ! loop
 integer :: i_found, j_found
 integer :: i_bot,j_bot
 real(sdsu_fpd) ::  vector_inc(3)  ! incremental vector [km]
 real(sdsu_fpd) ::  vector_bot(3)  ! bottom vector [km]
 logical :: logic_found_grid
 
!
! initialize bottom and found grid index 
!
 i_found = i_ref 
 j_found = j_ref
 i_bot   = i_found
 j_bot   = j_found
 k_bot   = 1    

 bot%ir = real(i_bot)  !bottom i index in real format
 bot%jr = real(j_bot)  !bottom j index in real format

 good_slant      = .true.  ! initialize logical value

!
! ---------derive upwelling-slant (from bottom to top) properties -----------
!
 UP_LOOP: do k = 1, mxlyr_refine

    slant_terrain = 0 !initially assume atmosphere ( 0-air  1-terrain )

    !
    ! get vector increment of slant profile at every range_micro height
    !
    call vector_increment( range_radar, k, crm, vector_inc )

    !
    ! distance from satellite to 
    !
    call get_distance( vector_inc, crm%vector_sat, dist_from_sat(k) )

    !
    ! For k==1 or previous k is terrain, update bottom vector.
    ! bottom vector is defined in the middle point of bottom atmospheric layer. 
    !
    if(k==1 .or. slant_terrain == 1) then 
       vector_bot = vector_inc  ! [km]
    endif

    !
    ! find out if increment vector is within the CRM grid points or not. 
    !
    call find_nearest_point( k, i_found, j_found, vector_inc , logic_found_grid, &
                             slant_terrain )
 
    IF_GOOD_SLANT: if( .not. logic_found_grid  ) then
        good_slant = .false.
        return 
    else

     if( atmos_opt(i_found,j_found,k,1)%kext == undefined ) then
         slant_terrain           = 1 ! terrain
         slant_opt1d(k,:)%kext  = undefined
         slant_opt1d(k,:)%salb  = undefined
         slant_opt1d(k,:)%asym  = undefined
         slant_opt1d(k,:)%sback = undefined

         L2_radar_slant1d(k)%press     = undefined
         L2_radar_slant1d(k)%t_air     = undefined
         L2_radar_slant1d(k)%sh        = undefined
         L2_radar_slant1d(k)%w         = undefined
         L2_radar_slant1d(k)%iwc       = undefined
         L2_radar_slant1d(k)%lwc       = undefined
         L2_radar_slant1d(k)%rwc       = undefined
         L2_radar_slant1d(k)%lfrac     = undefined
         L2_radar_slant1d(k)%bulk_den  = undefined
         L2_radar_slant1d(k)%rain_rate = undefined

         do m = 0, 6
            L2_radar_slant1d(k)%mmt(m) = undefined  
         enddo

     else

         slant_terrain           = 0 ! atmosphere
         slant_opt1d(k,:)%kext  = atmos_opt(i_found,j_found,k,:)%kext
         slant_opt1d(k,:)%salb  = atmos_opt(i_found,j_found,k,:)%salb
         slant_opt1d(k,:)%asym  = atmos_opt(i_found,j_found,k,:)%asym
         slant_opt1d(k,:)%sback = atmos_opt(i_found,j_found,k,:)%sback
         L2_radar_slant1d(k)%press     = L2_radar_refine(i_found,j_found,k)%press
         L2_radar_slant1d(k)%t_air     = L2_radar_refine(i_found,j_found,k)%t_air
         L2_radar_slant1d(k)%sh        = L2_radar_refine(i_found,j_found,k)%sh
         L2_radar_slant1d(k)%w         = L2_radar_refine(i_found,j_found,k)%w
         L2_radar_slant1d(k)%iwc       = L2_radar_refine(i_found,j_found,k)%iwc
         L2_radar_slant1d(k)%lwc       = L2_radar_refine(i_found,j_found,k)%lwc
         L2_radar_slant1d(k)%rwc       = L2_radar_refine(i_found,j_found,k)%rwc
         L2_radar_slant1d(k)%lfrac     = L2_radar_refine(i_found,j_found,k)%lfrac
         L2_radar_slant1d(k)%bulk_den  = L2_radar_refine(i_found,j_found,k)%bulk_den
         L2_radar_slant1d(k)%rain_rate = L2_radar_refine(i_found,j_found,k)%rain_rate

         do m = 0, 6
            L2_radar_slant1d(k)%mmt(m) = L2_radar_refine(i_found,j_found,k)%mmt(m)
         enddo

     endif
 
    endif IF_GOOD_SLANT

    !
    ! if incremented point is terrain, update surafce parameters
    !
    if( slant_terrain == 1 ) then 
        k_bot             = k+1
        i_bot             = i_found
        j_bot             = j_found

!
!  put sigma 0 here
!

!        slant_t_skin      = t_skin_micro(i_bot,j_bot    )  
!        slant_emiss(:,:)  = emis_micro  (i_bot,j_bot,:,:)
!        slant_ebar (:)    = ebar_micro  (i_bot,j_bot,:  )
        bot%ir = real(i_bot)
        bot%jr = real(j_bot)
    endif

 enddo UP_LOOP


 return
 end subroutine slant_path1d_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine find_mirror_position( v_sat, v_bot, mirror)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For a given satellite and bottom vectors, it find mirror position in contrast to the 
!            satellite vector.  
!
! History:
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: v_sat(3)  ! satellite vector  [km]
 real(sdsu_fpd),intent(in) :: v_bot(3)  ! bottom vector  [km]
 type(surface_CRM_parameter),intent(out) :: mirror

 real(sdsu_fpd) :: alpha, beta, gam
 real(sdsu_fpd) :: a,b,c 
 real(sdsu_fpd) :: a_dash, c_dash
 real(sdsu_fpd) :: v_center(3)  ! center of v_sat and v_mirror  [km]
 real(sdsu_fpd) :: v_mirror(3)  ! bottom vector  [km]
 
!
! get triangle parameters
!
 call triangle_parameter( v_sat, v_bot, alpha, beta, gam, a, b, c  )

 a_dash = b*dsin(0.5*dpi-gam) / dsin(0.5*dpi)  ![km]
 c_dash = b*dsin(gam)        / dsin(0.5*dpi)  ![km]

 v_center(1:3) = v_bot(1:3) * a_dash / a  !center of v_sat and v_mirror [km]

 v_mirror(1:3) =  v_center(1:3) - (v_sat(1:3) - v_center(1:3)) !vector of satellite mirror porision
                                                               ! relative to bottom vector [km]


!
! output to mirror data to the original surface_CRM data
!
   mirror%vector_grid(1:3) = v_bot(1:3)    ! vector position of CRM surface grid  [km]
   mirror%vector_sat(1:3)  = v_mirror(1:3) ! vector position of satellite [km]
   call cartesian2latlon( v_mirror(1:3), mirror%latlon_sat(1), mirror%latlon_sat(2) )
   mirror%a  = a
   mirror%b  = b
   mirror%c  = c
   mirror%alpha = alpha
   mirror%beta  = beta
   mirror%gam   = gam

 return
 end subroutine find_mirror_position

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 subroutine find_nearest_point( k, i_found, j_found, v_inc, logic_found_grid, &
                                slant_terrain )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Find nearest optical paorperties, temperature, and surface properties 
!            along the computed grid-satellite vector (slant path).
!
! History:
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer,intent(in) :: k
 integer,intent(inout) :: i_found, j_found  !i j grid point that match in preivous step
 real(sdsu_fpd),intent(in) :: v_inc(3)  !vector inrement [km]
 logical,intent(out) :: logic_found_grid
 integer,intent(inout) :: slant_terrain  !terrain index 0 or 1

 integer :: inc
 integer :: n,i,j, is,ie,js,je
 integer :: n_loc(1:1),n_min,j_min,i_min
 integer,parameter :: dbox = 1  !toshii this could be linked to off-nadir angle and range_micro
 integer,parameter :: nmax = (1+dbox*2)**2

 real(sdsu_fpd) :: v_inc_ellip(3)  !modified vector inrement at reference surface [km]
 real(sdsu_fpd) :: lat_inc, lon_inc , dist, dist_min, dist_crit

 type neighbor_grid
   integer :: i
   integer :: j
   real(sdsu_fpd) :: dist ! physical distance [km]
 end type neighbor_grid
 type( neighbor_grid ) :: neighbor(nmax)

 !
 ! This is tricky, but...here is the implementation. 
 ! 1. get lat lon of incremental vector
 call cartesian2latlon( v_inc,lat_inc,lon_inc )
 ! 2. get vector (at reference ellipoid) from lat lon above
 call latlon2cartesian(lat_inc,lon_inc,v_inc_ellip)
 ! 3. this is required since distance between incremental vector and neiboring grid 
 !    will be calculated based on surface reference ellipsoid. 

!
! critical distance for searching criteria
!
 dist_crit = dsqrt( dble(gridsize**2. + gridsize**2) )

 ! search i,j boundary
 is = min(mxgridx, max(1, i_found - dbox))
 ie = min(mxgridx, max(1, i_found + dbox))

 js = min(mxgridy, max(1, j_found - dbox))
 je = min(mxgridy, max(1, j_found + dbox))

 !
 ! search loop
 ! 
 n = 0  ; neighbor%dist = 1.d6
 SEARCH_LOOP: do j = js,je ; do i = is,ie
 n = n+1
   call get_distance(surface_CRM(i,j)%vector_grid, v_inc_ellip, dist )
   neighbor(n)%i    = i 
   neighbor(n)%j    = j
   neighbor(n)%dist = dist
 enddo ; enddo SEARCH_LOOP

!
! find the shortest distance grid and its index
!
 n_loc(1:1) = minloc( neighbor(1:n)%dist )
 n_min = n_loc(1)
 i_min = neighbor(n_min)%i
 j_min = neighbor(n_min)%j
 dist_min = neighbor(n_min)%dist

!
! nearest neighbor grid is close enough or not
! if not, incrementa vector goes to out of domain --> bad slant
!
 if( dist_min < dist_crit ) then !find good neighbor

     logic_found_grid        = .true.
    ! update i_found and j_found 
    i_found = i_min
    j_found = j_min

 else ! this means all slant path parameters will be undefined
    logic_found_grid        = .false.
 endif

 return
 end  subroutine find_nearest_point

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine find_nearest_grid( k, i_found, j_found, v_inc, logic_found_grid, &
                               slant_opt0d , slant_t_air, slant_terrain )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Find nearest optical paorperties, temperature, and surface properties 
!            along the computed grid-satellite vector (slant path).
!
! History:
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer,intent(in) :: k
 integer,intent(inout) :: i_found, j_found  !i j grid point that match in preivous step
 real(sdsu_fpd),intent(in) :: v_inc(3)  !vector inrement [km]
 logical,intent(out) :: logic_found_grid
 type (opt_micro_parameter),intent(inout) :: slant_opt0d(:)  !interpolated optical propeties along slant path. 
 real(sdsu_fps),intent(inout) :: slant_t_air !air temperature in slant path [K]
 integer,intent(inout) :: slant_terrain  !terrain index 0 or 1

 integer :: inc
 integer :: n,i,j, is,ie,js,je
 integer :: n_loc(1:1),n_min,j_min,i_min
 integer,parameter :: dbox = 1  !toshii this could be linked to off-nadir angle and range_micro
 integer,parameter :: nmax = (1+dbox*2)**2

 real(sdsu_fpd) :: v_inc_ellip(3)  !modified vector inrement at reference surface [km]
 real(sdsu_fpd) :: lat_inc, lon_inc , dist, dist_min, dist_crit

 type neighbor_grid
   integer :: i
   integer :: j
   real(sdsu_fpd) :: dist ! physical distance [km]
 end type neighbor_grid
 type( neighbor_grid ) :: neighbor(nmax)

 !
 ! This is tricky, but...here is the implementation. 
 ! 1. get lat lon of incremental vector
 call cartesian2latlon( v_inc,lat_inc,lon_inc )
 ! 2. get vector (at reference ellipoid) from lat lon above
 call latlon2cartesian(lat_inc,lon_inc,v_inc_ellip)
 ! 3. this is required since distance between incremental vector and neiboring grid 
 !    will be calculated based on surface reference ellipsoid. 

!
! critical distance for searching criteria
!
 dist_crit = dsqrt( dble(gridsize**2. + gridsize**2) )

 ! search i,j boundary
 is = min(mxgridx, max(1, i_found - dbox))
 ie = min(mxgridx, max(1, i_found + dbox))

 js = min(mxgridy, max(1, j_found - dbox))
 je = min(mxgridy, max(1, j_found + dbox))

 !
 ! search loop
 ! 
 n = 0  ; neighbor%dist = 1.d6
 SEARCH_LOOP: do j = js,je ; do i = is,ie
 n = n+1
   call get_distance(surface_CRM(i,j)%vector_grid, v_inc_ellip, dist )
   neighbor(n)%i    = i 
   neighbor(n)%j    = j
   neighbor(n)%dist = dist
 enddo ; enddo SEARCH_LOOP

!
! find the shortest distance grid and its index
!
 n_loc(1:1) = minloc( neighbor(1:n)%dist )
 n_min = n_loc(1)
 i_min = neighbor(n_min)%i
 j_min = neighbor(n_min)%j
 dist_min = neighbor(n_min)%dist

!
! nearest neighbor grid is close enough or not
! if not, incrementa vector goes to out of domain --> bad slant
!
 if( dist_min < dist_crit ) then !find good neighbor

     logic_found_grid        = .true.
     if( atmos_opt(i_min,j_min,k,1)%kext == undefined ) then
         slant_terrain           = 1 ! terrain or upper stratosphere
         slant_t_air             = undefined
         slant_opt0d(:)%kext  = undefined
         slant_opt0d(:)%salb  = undefined
         slant_opt0d(:)%asym  = undefined
         slant_opt0d(:)%sback = undefined
     else
         slant_terrain           = 0 ! atmosphere
         slant_t_air             = t_air_micro(i_min,j_min,k)
         slant_opt0d(:)%kext  = atmos_opt(i_min,j_min,k,:)%kext
         slant_opt0d(:)%salb  = atmos_opt(i_min,j_min,k,:)%salb
         slant_opt0d(:)%asym  = atmos_opt(i_min,j_min,k,:)%asym
         slant_opt0d(:)%sback = atmos_opt(i_min,j_min,k,:)%sback
     endif
   ! update i_found and j_found 
   i_found = i_min
   j_found = j_min

 else ! this means all slant path parameters will be undefined
    logic_found_grid        = .false.
 endif

 return
 end  subroutine find_nearest_grid

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine vector_increment(range_hgt, k, crm ,v_inc)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Compute incremental vector for satellite-grid slant path.  
!
! History:
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: range_hgt
 integer,intent(in) :: k         !vertical level index       
 type(surface_CRM_parameter),intent(in) :: crm
 real(sdsu_fpd) :: a_dash,c_dash, dc
 real(sdsu_fpd) :: beta_dash, gam_dash
 real(sdsu_fpd) :: dv(3), v_inc(3)  !vector [km]
 real(sdsu_fpd) :: lat_inc, lon_inc


!
! c_dash is the inrement length from reference ground point
! toward satellite position for every range_micro height. 
!
 a_dash = crm%a + dble(k-1)*range_hgt + 0.5d0 * range_hgt

 beta_dash = dpi -  dasin( crm%b * dsin( crm%alpha ) / a_dash )  !corresponding new beta [rad]

 gam_dash = dpi - beta_dash - crm%alpha !corresnponding new gamma [rad]

 c_dash = crm%b * dsin(gam_dash) / dsin(beta_dash) !corresponding new c [km]


 !
 ! this is increment of c
 !
 dc = crm%c-c_dash

 !
 ! get inremental vector point of dc
 !
 dv(1:3) = crm%vector_sat(1:3) - crm%vector_grid(1:3)  ![km]

 v_inc(1:3) = crm%vector_grid(1:3) + dv(1:3) * dc / crm%c ![km]

 return
 end subroutine vector_increment

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine convolution_micro
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Convolve CRM-grid microwave Tb on sensor-level Tb, taking account of antenna gain 
!            patterns of the passive-microwave instrument.   
!
! History:
! 03/2014  Toshi Matsui@NASA GSFC         : Bug fixed in ifov gain function input.
! 10/2012  Toshi Matsui@NASA GSFC         : Coarse CRM grid option was addded. 
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 
 integer :: i_scan , j_scan, i,j,l,m
 integer :: i_crm  , j_crm
 integer :: is_crm, ie_crm, js_crm, je_crm
 real(sdsu_fpd) :: bw  !beamwidth between FOV center and  [deg]
 real(sdsu_fps) :: tb_sum(2), emis_sum(2)
 integer :: cnt
 integer :: nf,np
 integer :: ierr
 integer :: iter
 integer :: ib,jb !i,j, CRM index of bottom of slant-path
 integer :: dscan 
 real(sdsu_fpd) :: wgt, wgt_prev  ! Gaussian gain weight [-]
 real(sdsu_fps) :: wgt_sum, wgt_max, wgt_per
 integer :: iter_efov  ! maximum iteration time
 real(sdsu_fpd),allocatable :: v_sat_iter(:,:), v_fov_iter(:,:)
 real(sdsu_fps) :: L2bundle(20)
 character(len=100),parameter :: convolution_pattern = 'efov_gaussian'  
                                                       !ifov_gaussian
                                                       !efov_gaussian

 if( masterproc ) then
   print*,''
   print*,'MSG convolution_micro: convolve Tb within fov beamwidth '
 endif


!
! allocate scan L1B output
!
 if( allocated(scans_micro) ) deallocate( scans_micro )
 allocate( scans_micro( imax_scan , js_scan:je_scan) , &
           stat=ierr )
 do np = 1, 2 ; do nf = 1, mxfreq_micro
    scans_micro%tb_fov  (nf,np) = undefined
    scans_micro%emis_fov(nf,np) = undefined
    scans_micro%wgt_fov (nf,np) = undefined
 enddo ; enddo

!
! If L2 simulator, get bundle 2D parameters. 
!
 call get_bundle_2d_micro


!
! initialize output
!
 tb_out_fov = undefined 
 L2_out_fov = undefined

!
! derive iteration time for EFOV gaussian option
!
 if( trim(convolution_pattern) == 'efov_gaussian' ) then 
    iter_efov = max(1, INT(stat%mean_fov_dist/dble(gridsize) ) ) ![km]
    allocate( v_sat_iter(3,1:iter_efov), v_fov_iter(3,1:iter_efov) )
 endif


!
! gain function independent on channel frequencies
!
 FREQ_LOOP: do nf = 1, mxfreq_micro

 !
 ! determin how many surrounding pixcels must be searched. 
 ! This algorithm search area twie large as default IFOV. 
 !
 dscan = NINT(  max( fov_dt_micro(nf) , fov_ct_micro(nf) ) /gridsize  ) + 1  !+1 added toshii

!
! scan sample loop to get arc distance betweeen CRM grid and each sample FOV
!
 SCAN_LOOP: do j_scan = js_scan , je_scan  ; do i_scan = 1, imax_scan  

    if( scans(i_scan,j_scan)%i_crm /= undefined_i2 ) then

       !
       ! determin search bounds within domain memory
       !
       is_crm = min(mxgridx, max(1,scans(i_scan,j_scan)%i_crm - dscan) )
       ie_crm = min(mxgridx, max(1,scans(i_scan,j_scan)%i_crm + dscan) )
       js_crm = min(mxgridy, max(1,scans(i_scan,j_scan)%j_crm - dscan) )
       je_crm = min(mxgridy, max(1,scans(i_scan,j_scan)%j_crm + dscan) )

       wgt_prev = 0.
       tb_sum = 0. ; emis_sum = 0. ; wgt_sum = 0. ; L2bundle = 0.
       do j_crm = js_crm,je_crm ; do i_crm = is_crm,ie_crm
          if( tb_out(i_crm,j_crm,nf,1) /= undefined) then

             select case(trim(convolution_pattern))
             case('ifov_flat','ifov_gaussian') ! flat gain function within IFOV        

                !
                ! derive angle between FOV center and surrounding CRM grid
                !
                call get_beamwidth( scans(i_scan,j_scan)%vector_sat, surface_CRM(i_crm,j_crm)%vector_grid, &
                                    scans(i_scan,j_scan)%vector_fov, bw )

                if( bw <= dble(beamwidth(nf)) ) then ! within IFOV

                    wgt = 1.  ! flat weight
                    if(trim(convolution_pattern)=='ifov_gaussian') &
                    call ifov_gaussian( dble(beamwidth(nf)), bw*degrad, wgt) ! gaussian gain function
                    tb_sum(1:2) = tb_sum(1:2) + tb_out(i_crm,j_crm,nf,1:2)*real(wgt)
                    emis_sum(1:2) = emis_sum(1:2) + emis_micro(i_crm,j_crm,nf,1:2)*real(wgt)
                    wgt_sum = wgt_sum + real(wgt)
                    ib = int(bot_slant(i_crm,j_crm)%ir)
                    jb = int(bot_slant(i_crm,j_crm)%jr)
                    L2bundle(1:nparam_L2) = L2bundle(1:nparam_L2) + L2bundle2d(ib,jb,1:nparam_L2) *real(wgt)

                endif

             case('efov_gaussian') ! gaussian gain function within EFOV 

               if( stat%grid_gt_fov ) then !course resolution CRM case (does not work for EFOV) 

                  !
                  ! derive angle between FOV center and surrounding CRM grid
                  !
                  call get_beamwidth( scans(i_scan,j_scan)%vector_sat, surface_CRM(i_crm,j_crm)%vector_grid, &
                                      scans(i_scan,j_scan)%vector_fov, bw )
                  call ifov_gaussian( dble(beamwidth(nf))*degrad, bw*degrad, wgt) ! gaussian gain function

                     if( wgt > wgt_prev ) then ! replace
                       tb_sum = 0. ; emis_sum = 0. ; wgt_sum = 0. ; L2bundle = 0.
                       tb_sum(1:2) = tb_out(i_crm,j_crm,nf,1:2)*real(wgt)
                       emis_sum(1:2) = emis_micro(i_crm,j_crm,nf,1:2)*real(wgt)
                       wgt_sum = real(wgt)
                       ib = int(bot_slant(i_crm,j_crm)%ir)
                       jb = int(bot_slant(i_crm,j_crm)%jr)
                       L2bundle(1:nparam_L2) = L2bundle2d(ib,jb,1:nparam_L2) *real(wgt)
                     endif
                     wgt_prev = wgt  !previous weight


               else ! fine-resolution CRM case  (get radar values from multiple CRM grids)
               
                 !
                 ! get aditional vector for EFOV iteration
                 !
                 call efov_vector_iter(i_scan,j_scan, iter_efov, v_sat_iter, v_fov_iter ) 

                 do iter = 1, iter_efov ! EFOV vector iteration
                    !
                    ! derive angle between FOV center and surrounding CRM grid
                    !
                    call get_beamwidth( v_sat_iter(1:3,iter), surface_CRM(i_crm,j_crm)%vector_grid, &
                                        v_fov_iter(1:3,iter), bw )
                    call ifov_gaussian( dble(beamwidth(nf))*degrad, bw*degrad, wgt) !gaussian gain function

                    tb_sum(1:2) = tb_sum(1:2) + tb_out(i_crm,j_crm,nf,1:2)*real(wgt)
                    emis_sum(1:2) = emis_sum(1:2) + emis_micro(i_crm,j_crm,nf,1:2)*real(wgt)
                    wgt_sum = wgt_sum + real(wgt)

                    ib = int(bot_slant(i_crm,j_crm)%ir)
                    jb = int(bot_slant(i_crm,j_crm)%jr)
                    L2bundle(1:nparam_L2) = L2bundle(1:nparam_L2) + L2bundle2d(ib,jb,1:nparam_L2) *real(wgt)

                 enddo 

               endif

             case default
                call stop_sdsu(myrank, 'MSG convolution_micro: no such convolution type')
             end select

          endif
       enddo ; enddo

       !
       ! get convolved Tb for tb_out_fov grid for presentation purpose
       ! (return to tb on CRM grid)
       !
       if( wgt_sum > 0. ) then
          scans_micro(i_scan,j_scan)%tb_fov(nf,1:2)   = tb_sum(1:2) / wgt_sum
          scans_micro(i_scan,j_scan)%emis_fov(nf,1:2) = emis_sum(1:2) / wgt_sum
          scans_micro(i_scan,j_scan)%wgt_fov(nf,1:2)  = wgt_sum  ! sum of weigting function

           if( nf == nf_L2 ) then
              scans_L2(i_scan,j_scan)%params(1:nparam_L2) = L2bundle(1:nparam_L2) /  wgt_sum
           endif

          !
          ! output purpose paint FOV Tb on CRM grid
          !
          do j_crm = js_crm,je_crm ; do i_crm = is_crm,ie_crm
              call get_beamwidth( scans(i_scan,j_scan)%vector_sat, surface_CRM(i_crm,j_crm)%vector_grid, &
                                  scans(i_scan,j_scan)%vector_fov, bw )

              if( bw <= 0.45*stat%mean_scan_angle/degrad  ) then ! within small FOV < mean scan angle
                 tb_out_fov(i_crm,j_crm,nf,1:2) = scans_micro(i_scan,j_scan)%tb_fov(nf,1:2)
                 if( nf == nf_L2 ) then
                    L2_out_fov(i_crm,j_crm,1:nparam_L2) = scans_L2(i_scan,j_scan)%params(1:nparam_L2)
                 endif
              endif
          enddo ; enddo

       else

          scans_micro(i_scan,j_scan)%tb_fov (nf,1:2) = undefined
          scans_micro(i_scan,j_scan)%wgt_fov(nf,1:2) = undefined
           if( nf == nf_L2 ) then
              scans_L2(i_scan,j_scan)%params(1:nparam_L2) = undefined
           endif
       endif


    endif ! i_undefined

 enddo ; enddo SCAN_LOOP

 enddo FREQ_LOOP


!
! deallocate allocated params
!
 if( trim(convolution_pattern) == 'efov_gaussian' ) then
    deallocate( v_sat_iter, v_fov_iter )
 endif


#if MPI == 2
if(numproc_tot>1) then
!
! distribute parameters across different CPU memory. 
!
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, L2_out_fov(:,:,:) )
  
   do np = 1, nparam_L2
      call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                      scans_L2(1:imax_scan,js_scan:je_scan)%params(np) )
   enddo

   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, tb_out_fov(:,:,:,:) )

   do np = 1, 2 ; do nf = 1, mxfreq_micro
      call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                      scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov   (nf,np) )
      call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                      scans_micro(1:imax_scan,js_scan:je_scan)%emis_fov (nf,np) )
      call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                      scans_micro(1:imax_scan,js_scan:je_scan)%wgt_fov  (nf,np) )
   enddo ; enddo

endif
#endif

!
! re-scale wgt_fov [%] (0~100) --> sampling quality flag (sample#*gain_function)
! 100 % : good sampling
!   0 % : poor sampling due to edge of CRM domains
!
 do np = 1, 2 ; do nf = 1, mxfreq_micro

   wgt_max = maxval( scans_micro(1:imax_scan,js_scan:je_scan)%wgt_fov(nf,np) ) !max weight

   do j_scan = js_scan , je_scan  
      
      do i_scan = 1, imax_scan

         wgt_per = scans_micro(i_scan,j_scan)%wgt_fov (nf,np) / wgt_max * 100. ! weight [%]  

         scans_micro(i_scan,j_scan)%wgt_fov (nf,np) = wgt_per  !weight [%]

!         if( wgt_per < 50. ) then  !bad quality
!             scans_micro(i_scan,j_scan)%tb_fov(nf,np) = undefined
!         endif 

       enddo

      !
      ! 2x1averaging for SSMIS sensor
      !
      if( trim(scan_type_micro) == 'SSMIS' .and. nf <= 3 ) then 
         do i_scan = 1, imax_scan-1
            if    ( scans_micro(i_scan,j_scan)%tb_fov(nf,np) /= undefined .and. &
                    scans_micro(i_scan+1,j_scan)%tb_fov(nf,np) /= undefined  ) then

                scans_micro(i_scan,j_scan)%tb_fov(nf,np) = 0.5*(&
                      scans_micro(i_scan  ,j_scan)%tb_fov(nf,np) +&
                      scans_micro(i_scan+1,j_scan)%tb_fov(nf,np))

            elseif( scans_micro(i_scan,j_scan)%tb_fov(nf,np) /= undefined .and. &
                    scans_micro(i_scan+1,j_scan)%tb_fov(nf,np) == undefined  ) then

                scans_micro(i_scan,j_scan)%tb_fov(nf,np) = scans_micro(i_scan  ,j_scan)%tb_fov(nf,np)

            elseif( scans_micro(i_scan,j_scan)%tb_fov(nf,np) == undefined .and. &
                    scans_micro(i_scan+1,j_scan)%tb_fov(nf,np) /= undefined  ) then

                scans_micro(i_scan,j_scan)%tb_fov(nf,np) = scans_micro(i_scan+1,j_scan)%tb_fov(nf,np)

            else

                scans_micro(i_scan,j_scan)%tb_fov(nf,np) = undefined

            endif

         enddo
      endif
 
   enddo 

 enddo ; enddo

 return
 end subroutine convolution_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine convolution_radar
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Convolve CRM-grid radar parameters on sensor-level footprint, taking account of antenna gain 
!            patterns of the radar instrument.   
!
! History:
! 01/2016  Toshi Matsui@NASA GSFC         : bug fix for Coarse CRM grid option
! 10/2012  Toshi Matsui@NASA GSFC         : Coarse CRM grid option was addded. 
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 
 integer :: i_scan , j_scan, i,j,l,m
 integer :: ic, jc
 integer :: i_crm  , j_crm, k_range
 integer :: is_crm, ie_crm, js_crm, je_crm
 real(sdsu_fpd) :: bw  !beamwidth between FOV center and  [deg]
 real(sdsu_fps) :: tb_sum(2), emis_sum(2)
 integer :: cnt
 integer :: nf,np
 integer :: ierr
 integer :: iter
 integer :: ib,jb !i,j, CRM index of bottom of slant-path
 integer :: dscan 
 real(sdsu_fpd) :: wgt, wgt_prev  ! Gaussian gain weight
 real(sdsu_fps) :: wgt_sum, wgt_max, wgt_per
 integer :: iter_efov  ! maximum iteration time
 real(sdsu_fpd),allocatable :: v_sat_iter(:,:), v_fov_iter(:,:)
 real(sdsu_fps) :: L2bundle(20)
 type( radar_bundle_parameter ), allocatable :: radar_wgt_sum1d(:)
 real(sdsu_fps), allocatable :: wgt_sum1d(:)
 character(len=100) :: convolution_pattern  ! ifov_gaussian' 
                                            ! efov_gaussian

 
 if( masterproc ) then
   print*,''
   print*,'MSG convolution_radar: convolve Z and L2 parameters within FOV'
 endif

!
! EFOV is used for CloudSat CPR only
!
! if( trim(scan_type_radar)  == 'CPR' ) then
!   convolution_pattern = 'efov_gaussian'
! else
   convolution_pattern = 'ifov_gaussian'
! endif

!
! derive iteration time for EFOV gaussian option (for CloudSat CPR)
!
 if( trim(convolution_pattern) == 'efov_gaussian' ) then 
    iter_efov = max(3, INT(stat%mean_scan_dist/dble(gridsize) ) ) ![km]
    allocate( v_sat_iter(3,1:iter_efov), v_fov_iter(3,1:iter_efov) )
 endif


!
! allocate/initialize scans_radar output
!
 if( allocated(radar_wgt_sum1d) ) deallocate(radar_wgt_sum1d)
 allocate( radar_wgt_sum1d(irange_max) )

 if( allocated(wgt_sum1d) ) deallocate(wgt_sum1d)
 allocate( wgt_sum1d(irange_max) )

!
! ISSUE. This parameters eat alot of memory ..... (toshii)
!       It often causes memory issue for bigger WRF domains... how to decompose?
!       It has not resolved yet. But it must be for GEOS-5 GMMS application. 
!
 if( allocated(scans_radar) ) deallocate( scans_radar )
 allocate( scans_radar( imax_scan , js_scan:je_scan, irange_max ) , &
           stat=ierr )
 ! EMK NUWRF
 if (ierr .ne. 0) then
    call stop_sdsu(myrank, 'MSG convolution_radar: ERROR allocating scans_radar')
 endif
 ! EMK NUWRF END

 do k_range = 1, irange_max ; do j_scan = js_scan , je_scan  ; do i_scan = 1, imax_scan
    call undefine_radar_bundle( scans_radar(i_scan,j_scan,k_range) )  !all parameters are undefined
 enddo ; enddo ; enddo


 FREQ_LOOP: do nf = 1, mxfreq_radar   ! mxfreq_radar must be one

 !
 ! determin how many surrounding pixcels must be searched. 
 ! This algorithm search area twice large as default IFOV to consider side-robe impact. 
 !
 dscan = NINT(  max( fov_dt_radar(nf) , fov_ct_radar(nf) ) / gridsize  ) + 1  !+1 added for stability (6 for DPR)

 !
 ! scan sample loop to get arc distance betweeen CRM grid and each sample FOV
 !
 SCAN_LOOP: do j_scan = js_scan , je_scan  ; do i_scan = 1, imax_scan  

    if( scans(i_scan,j_scan)%i_crm /= undefined_i2 ) then

       !
       ! determin search bounds within sub-domain bounds
       !
       is_crm = min(ibe, max(ibs,scans(i_scan,j_scan)%i_crm - dscan) )
       ie_crm = min(ibe, max(ibs,scans(i_scan,j_scan)%i_crm + dscan) )
       js_crm = min(jbe, max(jbs,scans(i_scan,j_scan)%j_crm - dscan) )
       je_crm = min(jbe, max(jbs,scans(i_scan,j_scan)%j_crm + dscan) )

       !initialize weight
       wgt_sum1d = 0. 
       do k_range = 1, irange_max
          call zero_radar_bundle( radar_wgt_sum1d(k_range) )
       enddo
  
       wgt_prev = 0. 
       CRM_PIX_LOOP: do j_crm = js_crm,je_crm ; do i_crm = is_crm,ie_crm

          if( maxval(radar_bundle3d(i_crm,j_crm,:)%tau) /= undefined) then  !toshii (check again )

             select case(trim(convolution_pattern))
             case('ifov_gaussian') ! flat gain function within IFOV  (TRMM PR, GPM DPR )

                !
                ! derive angle between FOV center and surrounding CRM grid
                !
                call get_beamwidth( scans(i_scan,j_scan)%vector_sat, surface_CRM(i_crm,j_crm)%vector_grid, &
                                    scans(i_scan,j_scan)%vector_fov, bw )

       
                if( stat%grid_gt_fov ) then !course resolution CRM case (get radar value from single CRM grid)       

!
! for coarse grid, radar beamwidth must be artificially enlarged otherwise wgt become zero
! by exceeding 2xbeamwidth in gaussina weighting function. So I put 10degree here, but 
! it does not really matter, since coarse mode only pick up nearest neighbor.
! 
                call ifov_gaussian( dble(10.)*degrad, bw*degrad, wgt) ! gaussian gain function


                   if( wgt > wgt_prev ) then ! replace
                      wgt_sum1d = 0.
                      do k_range = 1, irange_max
                         call zero_radar_bundle( radar_wgt_sum1d(k_range) )
                         call  sum_radar_bundle( real(wgt), radar_bundle3d(i_crm,j_crm,k_range), &
                                                 radar_wgt_sum1d(k_range), wgt_sum1d(k_range) )
                      enddo
                   endif
                   wgt_prev = wgt  !previous weight

                else ! fine-resolution CRM case  (get radar values from multiple CRM grids)

                call ifov_gaussian( dble(beamwidth_radar(nf))*degrad, bw*degrad, wgt) ! gaussian gain function
                  !
                  ! weighting sum
                  !

                  do k_range = 1, irange_max
                     call  sum_radar_bundle( real(wgt), radar_bundle3d(i_crm,j_crm,k_range), &
                                             radar_wgt_sum1d(k_range), wgt_sum1d(k_range) )
                  enddo

                endif

             case('efov_gaussian') ! gaussian gain function within EFOV  --> ( CloudSat CPR )
 
               !
               ! get aditional vector for EFOV iteration
               !
               call efov_vector_iter_CPR(i_scan,j_scan, iter_efov, v_sat_iter, v_fov_iter ) 

               do iter = 1, iter_efov ! EFOV vector iteration

                  !
                  ! derive angle between FOV center and surrounding CRM grid
                  !
                  call get_beamwidth( v_sat_iter(1:3,iter), surface_CRM(i_crm,j_crm)%vector_grid, &
                                      v_fov_iter(1:3,iter), bw )
                  call ifov_gaussian( dble(beamwidth_radar(nf))*degrad, bw*degrad, wgt)    ! gaussian gain function


                  if( stat%grid_gt_fov ) then !course resolution CRM case (get radar value from single CRM grid)       

                     if( wgt > wgt_prev ) then ! replace
                        do k_range = 1, irange_max
                           call zero_radar_bundle( radar_wgt_sum1d(k_range) )
                           call  sum_radar_bundle( real(wgt), radar_bundle3d(i_crm,j_crm,k_range), &
                                                   radar_wgt_sum1d(k_range), wgt_sum1d(k_range) )
                        enddo
                     endif
                     wgt_prev = wgt  !previous weight

                  else ! fine-resolution CRM case  (get radar values from multiple CRM grids)

                    !
                    ! weighting sum
                    !
                    do k_range = 1, irange_max
                       call  sum_radar_bundle( real(wgt), radar_bundle3d(i_crm,j_crm,k_range), &
                                               radar_wgt_sum1d(k_range), wgt_sum1d(k_range) )
                    enddo
                  endif
               enddo 

             case default
                call stop_sdsu(myrank, 'MSG convolution_radar: no such convolution type')
             end select

          endif
       enddo ; enddo  CRM_PIX_LOOP

!       endif


       !
       ! get convolved radar parameters 
       !
       wgt_max = maxval( wgt_sum1d(:) ) !maximum weight within radar range
       do k_range = 1, irange_max
          call ave_radar_bundle( radar_wgt_sum1d(k_range), wgt_sum1d(k_range), &
                                 wgt_max, scans_radar(i_scan,j_scan,k_range) )
       enddo 

    endif ! i_undefined

 enddo ; enddo SCAN_LOOP

 enddo FREQ_LOOP


!
! deallocate allocated params
!
 if( trim(convolution_pattern) == 'efov_gaussian' ) then
    deallocate( v_sat_iter, v_fov_iter )
 endif


#if MPI == 2
if(numproc_tot>1) then
!
! distribute defined parameters across different thread memory 
!
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%gain )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%Zt )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%Zm )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%tau )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%iwc )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%lwc )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%rwc )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%rain_rate )

#ifdef FULL_L2_RADAR

   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%press )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%t_air )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%sh )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%w  )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%lfrac )
   call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%bulk_den )
   do m = 0,6
      call mpi_sdsu_communicate_defined( 'TO_MASTER', undefined, &
                                       scans_radar(1:imax_scan,js_scan:je_scan,:)%mmt(m) )
   enddo

#endif

endif
#endif


 return
 end subroutine convolution_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine smooth_radar_bundle( prm_col )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! Both PR, DPR, and CPR does oversampling: radar sample resolution is twice as large as range resolution. 
! In order to acount it, here smoothing is goind. 
!
! History:
! 09/2012  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 type( radar_bundle_parameter ), intent(inout) :: prm_col(:) !sum of weighted parameter for each parameters
 integer :: m

 call range_smoothing(  prm_col(:)%gain      )
 call range_smoothing(  prm_col(:)%Zt        )
 call range_smoothing(  prm_col(:)%Zm        )
 call range_smoothing(  prm_col(:)%tau       )
 call range_smoothing(  prm_col(:)%iwc       )
 call range_smoothing(  prm_col(:)%lwc       )
 call range_smoothing(  prm_col(:)%rwc       )
 call range_smoothing(  prm_col(:)%rain_rate )

#ifdef FULL_L2_RADAR
 call range_smoothing(  prm_col(:)%press     )
 call range_smoothing(  prm_col(:)%t_air     )
 call range_smoothing(  prm_col(:)%sh        )
 call range_smoothing(  prm_col(:)%w         )
 call range_smoothing(  prm_col(:)%lfrac     )
 call range_smoothing(  prm_col(:)%bulk_den  )
  do m = 0,6
     call range_smoothing( prm_col(:)%mmt(m) )
  enddo
#endif

 return
 end subroutine smooth_radar_bundle

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine range_smoothing(var)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  vertical weighting average. 
!
! History:
! 09/2012  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------

 real(sdsu_fps),intent(inout) :: var(:)

 integer :: bnd(1)    ! upper memory bound for 3D array
 integer :: k
 real(sdsu_fps),allocatable :: tmp(:)
 real(sdsu_fps),parameter :: w1=0.25,w2=0.5,w3=0.25

! define upper bounds
 bnd = UBOUND(var)

 if( allocated(tmp) ) deallocate(tmp)
 allocate( tmp( 1:bnd(1) ) )
 tmp = undefined

!
! smoothing
!
 do k = 2, bnd(1)-1

    if( var(k-1) /= undefined .and. var(k) /= undefined .and. var(k+1) /= undefined )  then
         tmp(k) = ( w1*var(k-1)+w2*var(k)+w3*var(k+1) ) / (w1+w2+w3) 
    elseif(  var(k-1) == undefined .and. var(k) /= undefined .and. var(k+1) /= undefined ) then
         tmp(k) = ( w2*var(k)+w3*var(k+1) ) / (w2+w3)
    elseif(  var(k-1) /= undefined .and. var(k) == undefined .and. var(k+1) /= undefined ) then
         tmp(k) = ( w1*var(k-1)+w3*var(k+1) ) / (w1+w3)
    elseif(  var(k-1) /= undefined .and. var(k) /= undefined .and. var(k+1) == undefined ) then
         tmp(k) = ( w1*var(k-1)+w2*var(k) ) / (w1+w2)
    elseif(  var(k-1) == undefined .and. var(k) /= undefined .and. var(k+1) == undefined ) then
         tmp(k) = ( w2*var(k) ) / (w2)
    elseif(  var(k-1) == undefined .and. var(k) == undefined .and. var(k+1) /= undefined ) then
         tmp(k) = ( w3*var(k+1) ) / (w3)
    elseif(  var(k-1) /= undefined .and. var(k) == undefined .and. var(k+1) == undefined ) then
         tmp(k) = ( w1*var(k-1) ) / (w1)
    elseif(  var(k-1) == undefined .and. var(k) == undefined .and. var(k+1) == undefined ) then
         tmp(k) = undefined
    endif

 enddo

!
! output
!
 do k = 2, bnd(1)-1
     var(k) = tmp(k)
 enddo


 return
 end subroutine range_smoothing 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine undefine_radar_bundle( prm_sum ) 
 implicit none
 type( radar_bundle_parameter ), intent(inout) :: prm_sum   !sum of weighted parameter for each parameters
 integer :: m

  prm_sum%gain      = undefined
  prm_sum%Zt        = undefined
  prm_sum%Zm        = undefined
  prm_sum%tau       = undefined
  prm_sum%iwc       = undefined
  prm_sum%lwc       = undefined
  prm_sum%rwc       = undefined
  prm_sum%rain_rate = undefined

#ifdef FULL_L2_RADAR
  prm_sum%press     = undefined
  prm_sum%t_air     = undefined
  prm_sum%sh        = undefined
  prm_sum%w         = undefined
  prm_sum%lfrac     = undefined
  prm_sum%bulk_den  = undefined
  do m = 0,6
     prm_sum%mmt(m) = undefined
  enddo
#endif

 return
 end subroutine undefine_radar_bundle

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine zero_radar_bundle( prm_sum )
 implicit none
 type( radar_bundle_parameter ), intent(inout) :: prm_sum   !sum of weighted parameter for each parameters
 integer :: m

  prm_sum%gain      = 0.e0
  prm_sum%Zt        = 0.e0
  prm_sum%Zm        = 0.e0
  prm_sum%tau       = 0.e0
  prm_sum%iwc       = 0.e0
  prm_sum%lwc       = 0.e0
  prm_sum%rwc       = 0.e0
  prm_sum%rain_rate = 0.e0

#ifdef FULL_L2_RADAR
  prm_sum%press     = 0.e0
  prm_sum%t_air     = 0.e0
  prm_sum%sh        = 0.e0
  prm_sum%w         = 0.e0
  prm_sum%lfrac     = 0.e0
  prm_sum%bulk_den  = 0.e0
  do m = 0,6
     prm_sum%mmt(m) = 0.e0
  enddo
#endif

 return
 end subroutine zero_radar_bundle

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine sum_radar_bundle( wgt, prm_in, prm_sum, wgt_sum  )
 implicit none
 real(sdsu_fps),intent(in) :: wgt  !weight
 type( radar_bundle_parameter ), intent(in)    :: prm_in  ! input parameters
 type( radar_bundle_parameter ), intent(inout) :: prm_sum ! sum of weighted parameter for each parameters
 real(sdsu_fps), intent(inout) :: wgt_sum                 ! sum of weight for each param
 integer :: m

 if( prm_in%Zt == undefined ) then
! if( prm_in%t_air == undefined ) then

   return !do nothing

 else

   wgt_sum              = wgt_sum              + wgt

   prm_sum%Zt        = prm_sum%Zt        + prm_in%Zt        * wgt
   prm_sum%Zm        = prm_sum%Zm        + prm_in%Zm        * wgt
   prm_sum%tau       = prm_sum%tau       + prm_in%tau       * wgt
   prm_sum%iwc       = prm_sum%iwc       + prm_in%iwc       * wgt
   prm_sum%lwc       = prm_sum%lwc       + prm_in%lwc       * wgt
   prm_sum%rwc       = prm_sum%rwc       + prm_in%rwc       * wgt
   prm_sum%rain_rate = prm_sum%rain_rate + prm_in%rain_rate * wgt

#ifdef FULL_L2_RADAR
   prm_sum%press     = prm_sum%press     + prm_in%press     * wgt
   prm_sum%t_air     = prm_sum%t_air     + prm_in%t_air     * wgt
   prm_sum%sh        = prm_sum%sh        + prm_in%sh        * wgt
   prm_sum%w         = prm_sum%w         + prm_in%w         * wgt
   prm_sum%lfrac     = prm_sum%lfrac     + prm_in%lfrac     * wgt
   prm_sum%bulk_den  = prm_sum%bulk_den  + prm_in%bulk_den  * wgt
   do m = 0,6
      prm_sum%mmt(m) = prm_sum%mmt(m) + prm_in%mmt(m) * wgt
   enddo
#endif

 endif

 return
 end subroutine sum_radar_bundle

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine ave_radar_bundle( prm_sum, wgt_sum, wgt_max, scans )
 implicit none
 type( radar_bundle_parameter ), intent(in) :: prm_sum   !sum of weighted parameter for each parameters
 real(sdsu_fps), intent(in) :: wgt_sum   !sum of weight for each param
 real(sdsu_fps), intent(in) :: wgt_max   !sum of weight for each param
 type( radar_bundle_parameter ), intent(out)    :: scans   !input parameters
 integer :: m

 if( wgt_sum <= 0. ) then   !absolutely no sampling due to surface (all --> undefined)

   scans%gain      = undefined
   scans%Zt        = undefined
   scans%Zm        = undefined
   scans%tau       = undefined
   scans%iwc       = undefined
   scans%lwc       = undefined
   scans%rwc       = undefined
   scans%rain_rate = undefined

#ifdef FULL_L2_RADAR
   scans%press     = undefined
   scans%t_air     = undefined
   scans%sh        = undefined
   scans%w         = undefined
   scans%lfrac     = undefined
   scans%bulk_den  = undefined
   do m = 0,6
      scans%mmt(m) = undefined
   enddo
#endif

   return  

 else !full sampling (all defined)

   scans%gain      = wgt_sum           / wgt_max  !normalized gain function 
   scans%Zt        = prm_sum%Zt        / wgt_sum
   scans%Zm        = prm_sum%Zm        / wgt_sum
   scans%tau       = prm_sum%tau       / wgt_sum
   scans%iwc       = prm_sum%iwc       / wgt_sum
   scans%lwc       = prm_sum%lwc       / wgt_sum
   scans%rwc       = prm_sum%rwc       / wgt_sum
   scans%rain_rate = prm_sum%rain_rate / wgt_sum

#ifdef FULL_L2_RADAR
   scans%press     = prm_sum%press     / wgt_sum
   scans%t_air     = prm_sum%t_air     / wgt_sum
   scans%sh        = prm_sum%sh        / wgt_sum
   scans%w         = prm_sum%w         / wgt_sum
   scans%lfrac     = prm_sum%lfrac     / wgt_sum
   scans%bulk_den  = prm_sum%bulk_den  / wgt_sum
   do m = 0,6
      scans%mmt(m) = prm_sum%mmt(m)    / wgt_sum
   enddo
#endif

 endif

 end subroutine ave_radar_bundle

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine find_range_bounds
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Estimate low and high bounds of radar range bin for orbital output. 
!            This also estimate PIAsrt and PIAe level. 
! History:
! 10/2012  Toshi Matsui@NASA GSFC         : PIA_SRT and PIA_effective levels are added. 
! 03/2012  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer :: i_scan, j_scan, k_range
 real(sdsu_fpd) :: range_resolution, range_km
 real(sdsu_fpd) :: range_max, range_min  !maximum and minimum range [km]
 real(sdsu_fpd) :: dist_fov_sat          !distance between FOV and satellite [km]
 real(sdsu_fpd) :: umu_local  !cosine of local incident angle [-]
 real(sdsu_fpd) :: alpha , beta, gam, a, b, c !triangle parameter

 if(masterproc)  print*,'MSG find_range_bounds: find bounds of radar range for output'

!
! Prepration ( define basic radar range distance )
!    Based on the maximum satellite-FOV distance (off nadir 17deg),
!    define maximum range distance from satellite to targe point adjusted by
!    additional range distance below ground.
!
 range_max = stat%max_sat_fov_dist - minhgt_radar  ! maximum radar range [km]
 range_min = stat%min_sat_fov_dist - mxhgt_radar   ! minimum radar range  [km]

!
! define absolute maximum number of radar range (irange_max)
!
 irange_max =  int( (range_max - range_min)/range_radar )  !define the number of range

!
! allocate bounds parameters
!
 if(allocated(bounds)) deallocate(bounds)
 allocate( bounds(1:imax_scan,js_scan:je_scan) ) 


 SCAN_LOOP: do j_scan = js_scan , je_scan  ; do i_scan = 1, imax_scan

   !
   ! radar range parameters for scan array
   !
   do k_range = 1, irange_max

      ! range resolution
!      scans_radar(i_scan,j_scan,k_range)%dhgt = range_radar * 1.d+3   ! [m]

      ! range distance 
      scans_radar(i_scan,j_scan,k_range)%distance = ( range_min + dble(k_range-1) * range_radar )*1.d+3  ![m]
      range_km = scans_radar(i_scan,j_scan,k_range)%distance * 1.d-3  ![km]

      call triangle_parameter(  scans(i_scan,j_scan)%vector_sat,  scans(i_scan,j_scan)%vector_fov, &
                               alpha , beta, gam, a, b, c )
      umu_local = real( dcos( dpi - beta ) )  ! cosine of local incident angle 

      ! get FOV-satellite distance [km]
      call get_distance( scans(i_scan,j_scan)%vector_fov, scans(i_scan,j_scan)%vector_sat, dist_fov_sat )

      ! above ellipsoid height nomral to ellipsoid [m]
      scans_radar(i_scan,j_scan,k_range)%hgt = (dist_fov_sat - range_km ) * umu_local * 1.d+3 ! [m]

   enddo

   !
   ! find bounds of start and ending range 
   !
   BOUNDS_LOOP: do k_range = irange_max, 1, -1
      if( scans_radar(i_scan,j_scan,k_range)%hgt > 0. ) then
         bounds(i_scan,j_scan)%ke = k_range
         bounds(i_scan,j_scan)%ks = k_range -int(mxhgt_radar/range_radar) + 1  
         exit BOUNDS_LOOP
      endif
   enddo BOUNDS_LOOP

   !
   ! find bounds for PIAsrt level (1st range bin of 0.5 of gain )
   !
   bounds(i_scan,j_scan)%k_srt = bounds(i_scan,j_scan)%ke !initialize
   SRT_LOOP: do k_range = bounds(i_scan,j_scan)%ke, bounds(i_scan,j_scan)%ks, -1
      if( scans_radar(i_scan,j_scan,k_range)%gain >= 0.5 ) then
         bounds(i_scan,j_scan)%k_srt = k_range
         exit SRT_LOOP
      endif
   enddo SRT_LOOP

   !
   ! find bounds for PIAe level (1st range bin of 1.0 of gain )
   !
   bounds(i_scan,j_scan)%k_piae = bounds(i_scan,j_scan)%ke !initialize
   PIA_LOOP: do k_range =  bounds(i_scan,j_scan)%k_srt, bounds(i_scan,j_scan)%ks, -1
      if( scans_radar(i_scan,j_scan,k_range)%gain >= 1.0 ) then
         bounds(i_scan,j_scan)%k_piae = k_range
         exit PIA_LOOP
      endif
   enddo PIA_LOOP

 enddo ; enddo SCAN_LOOP

 return
 end subroutine find_range_bounds

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine efov_vector_iter(i_scan, j_scan, iter_efov, v_sat_iter, v_fov_iter )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For given satellite and surface grid position, this subroutine compute angles and length
!            of triangular originating from the earth core.   
!
! History:
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer,intent(in) :: i_scan,j_scan 
 integer,intent(in) :: iter_efov
 real(sdsu_fpd), intent(out) :: v_sat_iter(1:3,1:iter_efov)  ! vector of satellite [km]
 real(sdsu_fpd), intent(out) :: v_fov_iter(1:3,1:iter_efov)  ! vector of satellite [km]

 real(sdsu_fpd),dimension(3) :: vs_sat , ve_sat, vs_fov, ve_fov ! vector of satellite [km]
 real(sdsu_fpd),dimension(3) :: v_sat_inc , v_fov_inc !incremental vector  [km]
 integer :: iter

 !
 !  intepolate vector between FOV (half vector)
 !
 if(i_scan == 1 ) then

   ve_sat = 0.5d0*(scans(1,j_scan)%vector_sat + scans(2,j_scan)%vector_sat) !end vector [km]
   vs_sat = 2.d0*scans(1,j_scan)%vector_sat - ve_sat !start vector [km]

   ve_fov = 0.5d0*(scans(1,j_scan)%vector_fov + scans(2+1,j_scan)%vector_fov) !end vector [km]
   vs_fov = 2.d0*scans(1,j_scan)%vector_fov - ve_fov !start vector [km]

 elseif(i_scan == imax_scan ) then

   vs_sat = 0.5d0*(scans(imax_scan,j_scan)%vector_sat + scans(imax_scan-1,j_scan)%vector_sat) !start vector [km]
   ve_sat = 2.d0*scans(imax_scan,j_scan)%vector_sat - vs_sat   !end vector [km]

   vs_fov = 0.5d0*(scans(imax_scan,j_scan)%vector_fov + scans(imax_scan-1,j_scan)%vector_fov) !start vector [km]
   ve_fov = 2.d0*scans(imax_scan,j_scan)%vector_fov - vs_fov  !end vector [km]


 else  !interpolate

   vs_sat = 0.5d0*(scans(i_scan,j_scan)%vector_sat + scans(i_scan-1,j_scan)%vector_sat) !start vector [km]
   ve_sat = 0.5d0*(scans(i_scan,j_scan)%vector_sat + scans(i_scan+1,j_scan)%vector_sat) !end vector [km]

   vs_fov = 0.5d0*(scans(i_scan,j_scan)%vector_fov + scans(i_scan-1,j_scan)%vector_fov) !start vector [km]
   ve_fov = 0.5d0*(scans(i_scan,j_scan)%vector_fov + scans(i_scan+1,j_scan)%vector_fov) !end vector [km]

 endif 

!
! start and endpoint
!
  v_sat_iter(1:3,1)         = vs_sat(1:3)
  v_sat_iter(1:3,iter_efov) = ve_sat(1:3)
  v_fov_iter(1:3,1)         = vs_fov(1:3)
  v_fov_iter(1:3,iter_efov) = ve_fov(1:3)
 
!
! interpolate middle point
! 
 if( iter_efov > 2 ) then
   v_sat_inc = (ve_sat - vs_sat) / dble( iter_efov -1 )
   v_fov_inc = (ve_fov - vs_fov) / dble( iter_efov -1 )
   do iter = 2, iter_efov-1
     v_sat_iter(1:3,iter) = vs_sat(1:3) + v_sat_inc(1:3)*dble(iter-1) ![km]
     v_fov_iter(1:3,iter) = vs_fov(1:3) + v_fov_inc(1:3)*dble(iter-1) ![km]
   enddo
 endif

 return
 end subroutine efov_vector_iter

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine efov_vector_iter_CPR(i_scan, j_scan, iter_efov, v_sat_iter, v_fov_iter )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For given satellite and surface grid position, this subroutine compute angles and length
!            of triangular originating from the earth core.   
!
! History:
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer,intent(in) :: i_scan,j_scan 
 integer,intent(in) :: iter_efov
 real(sdsu_fpd), intent(out) :: v_sat_iter(1:3,1:iter_efov)  ! vector of satellite [km]
 real(sdsu_fpd), intent(out) :: v_fov_iter(1:3,1:iter_efov)  ! vector of satellite [km]

 real(sdsu_fpd),dimension(3) :: vs_sat , ve_sat, vs_fov, ve_fov ! vector of satellite [km]
 real(sdsu_fpd),dimension(3) :: v_sat_inc , v_fov_inc !incremental vector  [km]
 integer :: iter

 !
 !  intepolate vector between FOV (half vector)
 !

 vs_sat = 0.5d0*(scans(i_scan,j_scan)%vector_sat + scans(i_scan,j_scan-1)%vector_sat) !start vector [km]
 ve_sat = 0.5d0*(scans(i_scan,j_scan)%vector_sat + scans(i_scan,j_scan+1)%vector_sat) !end vector [km]

 vs_fov = 0.5d0*(scans(i_scan,j_scan)%vector_fov + scans(i_scan,j_scan-1)%vector_fov) !start vector [km]
 ve_fov = 0.5d0*(scans(i_scan,j_scan)%vector_fov + scans(i_scan,j_scan+1)%vector_fov) !end vector [km]

!
! start and endpoint
!
  v_sat_iter(1:3,1)         = vs_sat(1:3)
  v_sat_iter(1:3,iter_efov) = ve_sat(1:3)
  v_fov_iter(1:3,1)         = vs_fov(1:3)
  v_fov_iter(1:3,iter_efov) = ve_fov(1:3)
 
!
! interpolate middle point 
! 
 if( iter_efov > 2 ) then
   v_sat_inc = (ve_sat - vs_sat) / dble( iter_efov -1 )
   v_fov_inc = (ve_fov - vs_fov) / dble( iter_efov -1 )
   do iter = 2, iter_efov-1
     v_sat_iter(1:3,iter) = vs_sat(1:3) + v_sat_inc(1:3)*dble(iter-1) ![km]
     v_fov_iter(1:3,iter) = vs_fov(1:3) + v_fov_inc(1:3)*dble(iter-1) ![km]
   enddo
 endif

 return
 end subroutine efov_vector_iter_CPR

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  get_beamwidth( v_sat, v_grid, v_fov, beamwidth )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For given satellite and surface grid position, this subroutine compute angles and length
!            of triangular originating from the earth core.   
!
! History:
! 09/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fpd), intent(in) :: v_sat(1:3)  ! vector of satellite [km]
 real(sdsu_fpd), intent(in) :: v_grid(1:3) ! vector of bottom grid [km]
 real(sdsu_fpd), intent(in) :: v_fov(1:3)  ! vector of fov [km]
 real(sdsu_fpd), intent(out) :: beamwidth  !  beam width angle [deg]

 real(sdsu_fpd) :: gam    ! [rad]
 real(sdsu_fpd) :: a,b,c  !  [km]

 call get_distance( v_sat  , v_grid, a)
 call get_distance( v_sat  , v_fov , b)
 call get_distance( v_grid , v_fov , c)

 gam   = dacos( ( a*a + b*b - c*c ) / (2.*a*b) )  ![rad]

 beamwidth = gam/degrad  ! [deg]

 return
 end subroutine  get_beamwidth

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine ifov_gaussian(beamwidth3db, alpha, f_gaussian)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For a given 3dB (half power) beam width of specific sensor, this routine give IFOV
!            gain function for various beam angle. Gain function is assumed in Gaussian distribution. 
!
! History:
! 04/2014  Toshi Matsui@NASA GSFC         : Add 2xbeamwidth limit (not significant change though)
! 10/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: beamwidth3db  ! sensor 3dB (half-power) beamwidth [rad]
 real(sdsu_fpd),intent(in) :: alpha         ! angle between FOV center and surround grid.  [rad]
 real(sdsu_fpd),intent(out) :: f_gaussian   ! gain function [dB]
 real(sdsu_fpd),parameter :: sqrt2 = 1.4142135381698608   ! sqrt(2.)
 real(sdsu_fpd),parameter :: sqrt2pi = 2.5066282749176025 ! sqrt(2.*acos(-1.0) )
 real(sdsu_fpd),parameter :: bunbo = 2.3548200130462646   ! 2.0*sqrt(2.0*log(2.0) ) 
 real(sdsu_fpd) :: sigma

!
! gain function must be considered within side robe 
! So, if measured angle is greater than 2xbeamwidth, gain function is zero.
!
 if( alpha > 2.d0 * beamwidth3db ) then
    f_gaussian = 0.d0
    return
 endif       

!
! variation (sigma) for a given beam width
!
 sigma = beamwidth3db / bunbo  !variation

!
! Gaussian gain function, mimicing [dB]
!
 f_gaussian = ( exp( -(alpha/(sqrt2*sigma))**2 ) )/(sigma*sqrt2pi)



 return
 end subroutine ifov_gaussian

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_orbital_micro_netcdf
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  This routine output satellite orbital data in NetCDF format. 
!
! History:
! 10/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 integer :: ncid
 integer :: i,j, is,ie,js,je !simplified index for i j bounds
 integer :: i_dimid, j_dimid, k_dimid
 integer, parameter :: ndims2d = 2, ndims1d=1
 integer :: dimids1d(ndims1d), dimids2d(ndims2d)
 integer :: id
 integer :: n ,ich, nf, np
 integer :: im  !middle of sample
 integer :: id_month, id_day, id_hour, id_min, id_sec, id_gps !var IDs (time)
 integer :: id_sat_alt, id_ecef_x, id_ecef_y, id_ecef_z, &
                        id_vel_x , id_vel_y , id_vel_z      ! var IDs (nagivation)
 integer :: id_fov_lat, id_fov_lon, id_sat_lat, id_sat_lon , id_inc ! var IDs (location)
 integer :: id_tb(max_chan), id_emis(max_chan)
 integer :: id_l2(nparam_L2) ! var IDs (L2 geophysical parameters)

 character(len=100) :: description
 character :: pol*3, cc*2
 real,allocatable :: temp2d(:,:)
 logical,parameter :: big_file  = .false.

 
!
! output file name
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(micro_sensor)//&
  '.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.nc'
 print("(a)"),'MSG write_orbital_micro_netcdf: output ->',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), NF90_64BIT_OFFSET, ncid) )

!
! set up array bounds
!
 im    = imax_scan/2 !
 if( imax_scan == 1 ) im = 1

! if(big_file) then
! for all scans
!    is = 1 ; ie = imax_scan ; js = 1 ; je = jmax_scan
! else
  ! for CRM domains-overlapped scans only. 
   is = 1 ; ie = imax_scan ; js = js_scan ; je = je_scan
! endif

 if( .not. allocated(temp2d) ) allocate( temp2d(is:ie,js:je) )

!
!  -------------------- define header file --------------------------
!

!
! Define dimensions
!
 call check( nf90_def_dim(ncid, 'n_sample' , ie-is+1 , i_dimid) )
 call check( nf90_def_dim(ncid, 'n_scan'   , je-js+1 , j_dimid) )
 dimids2d = (/ i_dimid, j_dimid /)
 dimids1d = (/ j_dimid /)

!
! time header & ID
!
 call check( nf90_def_var(ncid, 'month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[month]') )
 id_month = id

 call check( nf90_def_var(ncid, 'day_of_month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Days of Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[day]') )
 id_day = id

 call check( nf90_def_var(ncid, 'hour', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (hour) )') )
 call check( nf90_put_att(ncid, id, 'units', '[hour]') )
 id_hour = id

 call check( nf90_def_var(ncid, 'minute', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (minute) )') )
 call check( nf90_put_att(ncid, id, 'units', '[min]') )
 id_min = id

 call check( nf90_def_var(ncid, 'second', NF90_DOUBLE, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (second) )') )
 call check( nf90_put_att(ncid, id, 'units', '[sec]') )
 id_sec = id

 call check( nf90_def_var(ncid, 'GPS_atomic_time', NF90_DOUBLE, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' GPS atomic time (second) )') )
 call check( nf90_put_att(ncid, id, 'units', '[sec]') )
 id_gps = id


!
! satllite navigation
!
 call check( nf90_def_var(ncid, 'sat_altitude', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite altitude above Earth reference )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_sat_alt = id

 call check( nf90_def_var(ncid, 'sat_ecef_x', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite ECEF position (x-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_ecef_x = id

 call check( nf90_def_var(ncid, 'sat_ecef_y', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite ECEF position (y-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_ecef_y = id

 call check( nf90_def_var(ncid, 'sat_ecef_z', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite ECEF position (z-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_ecef_z = id

 call check( nf90_def_var(ncid, 'sat_vel_x', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite velocity (x-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km/sec]') )
 id_vel_x = id

 call check( nf90_def_var(ncid, 'sat_vel_y', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite velocity (y-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km/sec]') )
 id_vel_y = id

 call check( nf90_def_var(ncid, 'sat_vel_z', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite velocity (z-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km/sec]') )
 id_vel_z = id


!
! geolocation headher & ID
!
 call check( nf90_def_var(ncid, 'sat_latitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite geolocation (lat) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_sat_lat = id

 call check( nf90_def_var(ncid, 'sat_longitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite geolocation (lon) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_sat_lon = id

 call check( nf90_def_var(ncid, 'fov_latitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' center efov geolocation (lat) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_fov_lat = id

 call check( nf90_def_var(ncid, 'fov_longitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' center efov geolocation (lon) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_fov_lon = id

 call check( nf90_def_var(ncid, 'incident_angle', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' sensor incident angle') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_inc = id


!
! L1 brightness temperature 
!
 ich = 0
 do nf=1,mxfreq_micro !frequency loop
    do np=2,1,-1  !polarization loop
       if(np==1) pol = '(H)'
       if(np==2) pol = '(V)'

       !default
       description = 'Microwave Tb at '//trim(nch_micro(nf))//'Hz'//pol       
      
       select case(trim(scan_type_micro))
       case('GMI_LF','TMI_LF')
         if(nf==3 .and. np==1) cycle  !no 24G horizontal channel

       case('GMI_HF' )
         if(nf>=2 .and. np==1) cycle
         if(nf==2) description = 'Microwave Tb at 183.3plus-minus7GHz'//pol
         if(nf==3) description = 'Microwave Tb at 183.3plus-minus3GHz'//pol 
         if(nf>=4) cycle
 
       case('SSMIS')
          if(nf==2 .and. np==1) cycle  !no 22.235G horizontal channel (vertical only)
          if(nf==5 .and. np==2) cycle  !no 150G vertical channel
          if(nf>=6 .and. np==2) cycle  !no vertical sounder channel
          if(nf==6) description = 'Microwave Tb at 183.3plus-minus6.6GHz'//pol
          if(nf==7) description = 'Microwave Tb at 183.3plus-minus3.0GHz'//pol
          if(nf>=8) cycle

       case default
       end select

       ich = ich + 1
       write(cc,"(I2.2)") ich
       call check( nf90_def_var(ncid, 'channel_'//cc , NF90_REAL, dimids2d, id) )
       call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
       call check( nf90_put_att(ncid, id, 'description',trim(description) ) )
       call check( nf90_put_att(ncid, id, 'units', '[K]') )
       id_tb(ich) = id

    enddo
 enddo

!
! L1 emissivity
!
 ich = 0
 do nf=1,mxfreq_micro !frequency loop
    do np=2,1,-1  !polarization loop
       if(np==1) pol = '(H)'
       if(np==2) pol = '(V)'

       !default
       description = 'Surface emissivity at '//trim(nch_micro(nf))//'Hz'//pol

       select case(trim(scan_type_micro))
       case('GMI_LF','TMI_LF')
         if(nf==3 .and. np==1) cycle  !no 24G horizontal channel

       case('GMI_HF' )
         if(nf>=2 .and. np==1) cycle
         if(nf==2) description = 'Surface emissivity at 183.3plus-minus7GHz'//pol
         if(nf==3) description = 'Surface emissivity at 183.3plus-minus3GHz'//pol
         if(nf>=4) cycle

       case('SSMIS')
          if(nf==2 .and. np==1) cycle  !no 22.235G horizontal channel (vertical only)
          if(nf==5 .and. np==2) cycle  !no 150G vertical channel
          if(nf>=6 .and. np==2) cycle  !no vertical sounder channel
          if(nf==6) description = 'Surface emissivity at 183.3plus-minus6.6GHz'//pol
          if(nf==7) description = 'Surface emissivity at 183.3plus-minus3.0GHz'//pol
          if(nf>=8) cycle

       case default
       end select

       ich = ich + 1
       write(cc,"(I2.2)") ich
       call check( nf90_def_var(ncid, 'emis_'//cc , NF90_REAL, dimids2d, id) )
       call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
       call check( nf90_put_att(ncid, id, 'description',trim(description) ) )
       call check( nf90_put_att(ncid, id, 'units', '[-]') )
       id_emis(ich) = id

    enddo
 enddo


!
! L2 output 
!
 call check( nf90_def_var(ncid, 'rain_rate', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Surface rainfall rate ') )
 call check( nf90_put_att(ncid, id, 'units', '[mm/hr]') )
 id_l2(1) = id

 call check( nf90_def_var(ncid, 'clwp', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' cloud (non-raining particle only) water path') )
 call check( nf90_put_att(ncid, id, 'units', '[kg/m2]') )
 id_l2(2) = id

 call check( nf90_def_var(ncid, 'rwp', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' rain (raining particle only) water path') )
 call check( nf90_put_att(ncid, id, 'units', '[kg/m2]') )
 id_l2(3) = id

 call check( nf90_def_var(ncid, 'iwp', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' ice (ice particle only) water path') )
 call check( nf90_put_att(ncid, id, 'units', '[kg/m2]') )
 id_l2(4) = id

 call check( nf90_def_var(ncid, 'freez_hgt', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description','0deg-C isotherm height') )
 call check( nf90_put_att(ncid, id, 'units', '[m]') )
 id_l2(5) = id

 call check( nf90_def_var(ncid, 'frac_land', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' land fraction within FOV (0~1) ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_l2(6) = id

 call check( nf90_def_var(ncid, 'skin', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' surface skin temperature') )
 call check( nf90_put_att(ncid, id, 'units', '[K]') )
 id_l2(7) = id

 call check( nf90_def_var(ncid, 'frac_veg', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' vegetation fraction (0~1)') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_l2(8) = id

 call check( nf90_def_var(ncid, 'soilm', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' soil moisture fraction') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_l2(9) = id

 call check( nf90_def_var(ncid, 'swind', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' surface 10m wind speed') )
 call check( nf90_put_att(ncid, id, 'units', '[m/sec]') )
 id_l2(10) = id

 call check( nf90_def_var(ncid, 'h2ocol', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' column water vapor') )
 call check( nf90_put_att(ncid, id, 'units', '[kg/m2]') )
 id_l2(11) = id


  call check( nf90_enddef(ncid) )  ! End define mode.

!
! ----------------------------- DUMP DATA ------------------------------------------
!


!
! dump time group
!
  call check( nf90_put_var(ncid, id_month, scans(im,js:je)%time%mm)  )
  call check( nf90_put_var(ncid, id_day  , scans(im,js:je)%time%dd)  )
  call check( nf90_put_var(ncid, id_hour , scans(im,js:je)%time%hh)  )
  call check( nf90_put_var(ncid, id_min  , scans(im,js:je)%time%mn)  )
  call check( nf90_put_var(ncid, id_sec  , scans(im,js:je)%time%ss)  )
  call check( nf90_put_var(ncid, id_gps  , scans(im,js:je)%GPSatomic))

!
! dump satellite navigation group
!
  call check( nf90_put_var(ncid, id_sat_alt, real(scans(im,js:je)%alt_sat)      ) )
  call check( nf90_put_var(ncid, id_ecef_x , real(scans(im,js:je)%vector_sat(1))) )
  call check( nf90_put_var(ncid, id_ecef_y , real(scans(im,js:je)%vector_sat(2))) )
  call check( nf90_put_var(ncid, id_ecef_z , real(scans(im,js:je)%vector_sat(3))) )
  call check( nf90_put_var(ncid, id_vel_x  , real(scans(im,js:je)%vel_sat(1))   ) )
  call check( nf90_put_var(ncid, id_vel_y  , real(scans(im,js:je)%vel_sat(2))   ) )
  call check( nf90_put_var(ncid, id_vel_z  , real(scans(im,js:je)%vel_sat(3))   ) )

!
! dump satellite FOV position and incident angle [deg]
!
  call check( nf90_put_var(ncid, id_sat_lat, real(scans(is:ie,js:je)%latlon_sat(1)/degrad)) )
  call check( nf90_put_var(ncid, id_sat_lon, real(scans(is:ie,js:je)%latlon_sat(2)/degrad)) )
  call check( nf90_put_var(ncid, id_fov_lat, real(scans(is:ie,js:je)%latlon_fov(1)/degrad)) )
  call check( nf90_put_var(ncid, id_fov_lon, real(scans(is:ie,js:je)%latlon_fov(2)/degrad)) )
  call check( nf90_put_var(ncid, id_inc    , real(scans(is:ie,js:je)%inci_ang     /degrad)) )


!
! L1 brightness temperature 
!
 ich = 0
 do nf=1,mxfreq_micro !frequency loop
    do np=2,1,-1  !polarization loop

       ! default value
       temp2d(is:ie,js:je) = scans_micro(is:ie,js:je)%tb_fov(nf,np)

       select case(trim(scan_type_micro))
       case('GMI_LF','TMI_LF')
         if(nf==3 .and. np==1) cycle  !no 24G horizontal channel
       case('GMI_HF' )
         if(nf>=2 .and. np==1) cycle  !no 183.31 horizontal channel
         if(nf==2) then
            temp2d(is:ie,js:je) = 0.5*( scans_micro(is:ie,js:je)%tb_fov(nf  ,np) &
                                      + scans_micro(is:ie,js:je)%tb_fov(nf+3,np) ) 
         endif
         if(nf==3) then
            temp2d(is:ie,js:je) = 0.5*( scans_micro(is:ie,js:je)%tb_fov(nf  ,np) &
                                      + scans_micro(is:ie,js:je)%tb_fov(nf+1,np) )
         endif
         if(nf>=4) cycle

       case('SSMIS')
          if(nf==2 .and. np==1) cycle  !no 22.235G horizontal channel (vertical only)
          if(nf==5 .and. np==2) cycle  !no 150G vertical channel
          if(nf>=6 .and. np==2) cycle  !no vertical sounder channel
          if(nf==6) then
             temp2d(is:ie,js:je) = 0.5*( scans_micro(is:ie,js:je)%tb_fov(nf  ,np) &
                                       + scans_micro(is:ie,js:je)%tb_fov(nf+3,np) )
          endif
          if(nf==7) then
             temp2d = 0.5*( scans_micro(is:ie,js:je)%tb_fov(nf  ,np) &
                          + scans_micro(is:ie,js:je)%tb_fov(nf+1,np) )
          endif
          if(nf>=8) cycle

       end select

       do j = js,je ; do i = is,ie
          if ( ISNAN( temp2d(i,j) ) )  temp2d(i,j) = undefined  !final cleanup
       enddo ; enddo

       ich = ich + 1
       call check( nf90_put_var(ncid, id_tb(ich)    , temp2d(is:ie,js:je) ) )

    enddo
 enddo


!
! L1 surface emissivity 
!
 ich = 0
 do nf=1,mxfreq_micro !frequency loop
    do np=2,1,-1  !polarization loop

       ! default value
       temp2d(is:ie,js:je) = scans_micro(is:ie,js:je)%emis_fov(nf,np)

       select case(trim(scan_type_micro))
       case('GMI_LF','TMI_LF')
         if(nf==3 .and. np==1) cycle  !no 24G horizontal channel
       case('GMI_HF' )
         if(nf>=2 .and. np==1) cycle
         if(nf==2) then
            temp2d(is:ie,js:je) = 0.5*( scans_micro(is:ie,js:je)%emis_fov(nf,np) &
                                   + scans_micro(is:ie,js:je)%emis_fov(nf+3,np) ) 
         endif
         if(nf==3) then
            temp2d(is:ie,js:je) = 0.5*( scans_micro(is:ie,js:je)%emis_fov(nf,np) &
                                       + scans_micro(is:ie,js:je)%emis_fov(nf+1,np) )
         endif
         if(nf>=4) cycle

       case('SSMIS')
          if(nf==2 .and. np==1) cycle  !no 22.235G horizontal channel (vertical only)
          if(nf==5 .and. np==2) cycle  !no 150G vertical channel
          if(nf>=6 .and. np==2) cycle  !no vertical sounder channel
          if(nf==6) then
             temp2d(is:ie,js:je) = 0.5*( scans_micro(is:ie,js:je)%emis_fov(nf,np) &
                                        + scans_micro(is:ie,js:je)%emis_fov(nf+3,np) )
          endif
          if(nf==7) then
             temp2d = 0.5*( scans_micro(is:ie,js:je)%emis_fov(nf,np) &
                           + scans_micro(is:ie,js:je)%emis_fov(nf+1,np) )
          endif
          if(nf>=8) cycle

       end select

       do j = js,je ; do i = is,ie
          if ( ISNAN( temp2d(i,j) ) )  temp2d(i,j) = undefined  !final cleanup
       enddo ; enddo

       ich = ich + 1
       call check( nf90_put_var(ncid, id_emis(ich)    , temp2d(is:ie,js:je) ) )

    enddo
 enddo

!
! dump geophysical parameter convolved within FOV  (L2 parameters)
!
 do n = 1, nparam_L2
    call check( nf90_put_var(ncid, id_l2(n), scans_L2(is:ie,js:je)%params(n) ) )
 enddo

 call check( nf90_close(ncid) )  ! close nc file


 return
 end subroutine write_orbital_micro_netcdf

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_1BGMI_HDF5
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  This routine output GPIm 1BGMI-like "GMI brightness temperature" in HDF5 format.  
!
! History:
! 10/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------

!
! Dimension definitions
!
 integer :: nscan                 ! Number of scans in the granule
 integer,parameter :: nchan1 = 9  ! Number of channels in Swath 1
 integer,parameter :: nchan2 = 4  ! Number of channels in Swath 2
 integer :: nchan
 integer,parameter :: npix1 = 202  ! Number of channels in Swath 1
 integer,parameter :: npix2 = 202  ! Number of channels in Swath 2
 integer :: npix ! Number of pixels in Swath 1 or 2
 integer,parameter :: ncolds1 = 9 ! Maximum number of cold samples in Swath 1
 integer,parameter :: ncolds2 = 9 ! Maximum number of cold samples in Swath 2
 integer,parameter :: nhots1  = 30 ! Maximum number of hot samples in Swath 1
 integer,parameter :: nhots2  = 25 ! Maximum number of hot samples in Swath 2
 integer,parameter :: ntherm  = 11 ! Number of hot load thermisters

#if HDF >= 1

!
! HDF5 parameters
!
 character(len=100) :: groupdata ! group/data character name
 integer(hid_t) :: file_id       ! File identifier
 integer(hid_t) :: dataset_id    ! Dataset identifier
 integer(hid_t) :: dataspace_id  ! Data space identifier
 integer(hid_t) :: group_1st_id, group_2nd_id, group_3rd_id  ! Group identifier
 integer(HSIZE_T),dimension(1) :: dims1d
 integer(HSIZE_T),dimension(2) :: dims2d
 integer(HSIZE_T),dimension(3) :: dims3d   ! Datasets dimensions
  INTEGER(hid_t) :: h5_kind_type_i ! HDF type corresponding to the specified KIND
  INTEGER(HID_T) :: attr_id       ! Attribute identifier
  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  INTEGER     ::   arank = 1                      ! Attribure rank
  INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
  CHARACTER(LEN=80), DIMENSION(2) ::  attr_data  ! Attribute data
  CHARACTER(LEN=9) :: aname   ! Attribute name


 integer ::i, j
 integer :: error ! Error flag


!special byte parameter
! integer,parameter :: one_byte  = selected_int_kind(4) ! one Byte integer
! integer(kind=one_byte),allocatable :: temp1d_int1(:)
 integer,allocatable :: temp1d_int1(:)
 real,allocatable :: temp2d_real(:,:)
 integer ::   rank  ! Datasets rank
 character(len=100) :: group_1st, group_2nd, group_3rd !gorup name
 character(len=100) :: dsetname  !dataset name
 integer :: im       ! middle scan index

 real(sdsu_fps), allocatable :: tb_orbit(:,:,:)  ! microwave brightness temperature
 integer :: irec
 integer :: ns, ne  !nscan start and end index
 integer :: np, nf
 character(len=2) :: tag
 logical :: fexist


!
! dimension
!
 select case(trim(scan_type_micro))
 case('GMI_LF')
  ns = 1 ; ne = jmax_scan - 4
  nscan = ne-ns+1 
  nchan = nchan1
  npix  = npix1  
  tag   = 'S1'
 case('GMI_HF')
  ns = 4 ; ne = jmax_scan  !4-scan backward
  nscan = ne-ns+1 
  nchan = nchan2
  npix  = npix2  
  tag   = 'S2'
 case default
   call stop_sdsu(myrank, 'MSG write_1BGMI_HDF5:  This routine is for GMI channels only')
 end select

 im    = imax_scan/2 !
 if( imax_scan == 1 ) im = 1

!
! Initialize FORTRAN interface.
!
 call h5open_f(error)

!
! Create a new output file (if present replace)
! 
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//&
 '.1BGMI.'//tag //'.'//trim(overpass_tag)//trim(output_suffix)//'.h5'
 print("(a)"),'MSG write_1BGMI_HDF5: GPM 1BGMI output ->',trim(sdsu_io_file)

! inquire(file=sdsu_io_file,exist=fexist) !inquire file exist or not
! if (fexist) then !open file
!   call h5fopen_f(trim(sdsu_io_file), H5F_ACC_RDWR_F, file_id, error)
! else !create file
   call h5fcreate_f(trim(sdsu_io_file), H5F_ACC_TRUNC_F, file_id, error)
! endif

 !
 ! Create swath (1st layer) group  in the root group using absolute name.
 !
 group_1st = '/'//tag
 call h5gcreate_f(file_id, trim(group_1st), group_1st_id, error)

 !
 ! Store ScanTime group
 !
 groupdata = tag//'_ScanTime'
 call h5gcreate_f(group_1st_id, trim(groupdata), group_2nd_id, error) ! Create group (S#/S#_ScanTime)

  rank = 1  ; dims1d(1) = nscan 
  call h5screate_simple_f(rank, dims1d, dataspace_id, error)  ! Create the data space (nscan)

  dsetname = tag//'_Year'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%yyyy, dims1d, error)  

!  CALL h5screate_simple_f(arank, adims, aspace_id, error)
!  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error) 
!  aname = 'unit year'
!  call h5acreate_f(dataset_id, aname, atype_id, aspace_id, attr_id, error)
!  call h5awrite_f(attr_id, atype_id, attr_data, dims1d, error)
!  call h5aclose_f(attr_id, error)

  call h5dclose_f(dataset_id, error)


!following does not work --> one byte integer does not much h5write_f
  if( allocated(temp1d_int1) ) deallocate(temp1d_int1)
  allocate( temp1d_int1(ns:ne) )
  temp1d_int1(ns:ne) =  scans(im,ns:ne)%time%mm

  dsetname = tag//'_Month'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%mm , dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_DayOfMonth'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%dd, dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_Hour'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%hh, dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_Minute'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%mn, dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_Second'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, scans(im,ns:ne)%time%ss, dims1d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_DayOfYear'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%dayofyear, dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_SecondOfDay'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, scans(im,ns:ne)%time%secofday, dims1d, error)  
  call h5dclose_f(dataset_id, error)

 call h5sclose_f(dataspace_id, error)  ! close the data space (nscan)
 call h5gclose_f(group_2nd_id, error)  ! close the group (S#/S#_ScanTime)

!
! Store FOV lat lon
!
 rank = 2  ; dims2d(1) = npix ;  dims2d(2) = nscan
 call h5screate_simple_f(rank, dims2d, dataspace_id, error)  ! Create the data space (npix x nscan)

  dsetname = tag//'_Latitude'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(:,ns:ne)%latlon_fov(1)/degrad), dims2d, error)  
  call h5dclose_f(dataset_id, error)
 
  dsetname = tag//'_Longitude'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(:,ns:ne)%latlon_fov(2)/degrad), dims2d, error) 
  call h5dclose_f(dataset_id, error)

 call h5sclose_f(dataspace_id, error)  ! close the data space (npix x nscan)

 !
 ! navigation group
 !
 groupdata = tag//'_navigation'
 call h5gcreate_f(group_1st_id, trim(groupdata), group_2nd_id, error) ! Create group (/S#/S#_navigation)

  rank = 2  ; dims2d(1) = 3 ;  dims2d(2) = nscan
  call h5screate_simple_f(rank, dims2d, dataspace_id, error)  ! Create the data space (3xnscan)

  if( allocated(temp2d_real) ) deallocate(temp2d_real)
  allocate( temp2d_real(3,ns:ne) ) 
  do i = 1,3 ; do j = ns,ne
     temp2d_real(i,j) = real(scans(im,j)%vector_sat(i))*1000. ![m]
  enddo ; enddo

  dsetname = tag//'_scPos'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, temp2d_real, dims2d, error)  
  call h5dclose_f(dataset_id, error)

  do i = 1,3 ; do j = ns,ne
     temp2d_real(i,j) = real(scans(im,j)%vel_sat(i))*1000. ![m/sec]
  enddo ; enddo
  dsetname = tag//'_scVel'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, temp2d_real, dims2d, error)  
  call h5dclose_f(dataset_id, error)

  call h5sclose_f(dataspace_id, error)  ! close the data space (3xnscan)


  rank = 1  ; dims1d(1) = nscan
  call h5screate_simple_f(rank, dims1d, dataspace_id, error)  ! Create the data space (nscan)

  dsetname = tag//'_scLat'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(im,ns:ne)%latlon_sat(1)/degrad), dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_scLon'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(im,ns:ne)%latlon_sat(2)/degrad), dims1d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_scAlt'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(im,ns:ne)%alt_sat*1.d3), dims1d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_timeMidScan'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, scans(im,ns:ne)%GPSatomic, dims1d, error)
  call h5dclose_f(dataset_id, error)

 call h5sclose_f(dataspace_id, error)  ! close the data space (nscan)
 call h5gclose_f(group_2nd_id, error)  ! close the group (/S#/S#_navigation)

!
! incidence angle
!
 rank = 2  ; dims2d(1) = npix ;  dims2d(2) = nscan
 call h5screate_simple_f(rank, dims2d, dataspace_id, error)  ! Create the data space (npix x nscan)

  dsetname = tag//'_incidenceAngle'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(:,ns:ne)%inci_ang/degrad), dims2d, error)
  call h5dclose_f(dataset_id, error)

 call h5sclose_f(dataspace_id, error)  ! close the data space (npix x nscan)


!
! brightness temperature
!
 rank = 3  ; dims3d(1) = nchan ; dims3d(2) = npix ;  dims3d(3) = nscan
 call h5screate_simple_f(rank, dims3d, dataspace_id, error)  ! Create the data space (nchan x npix x nscan)

 if( allocated(tb_orbit) ) deallocate(tb_orbit)

 select case(trim(scan_type_micro))
 case('GMI_LF')
  allocate( tb_orbit(nchan,npix,ns:ne) )

  irec = 1
  do nf=1,mxfreq_micro ! frequency loop
     do np=2,1,-1      ! polarization loop

        if(nf==3 .and. np==1) cycle  ! no 24G(H) channel
        do j = ns,ne ; do i = 1, imax_scan
           tb_orbit(irec,i,j) = scans_micro(i,j)%tb_fov (nf,np)
        enddo ; enddo
        irec = irec + 1

     enddo
  enddo

 case('GMI_HF')

  allocate( tb_orbit(nchan,npix,ns:ne) )

  do j = ns,ne ; do i = 1, imax_scan
     tb_orbit(1,i,j) = scans_micro(i,j)%tb_fov (1,2)  !166(V)
     tb_orbit(2,i,j) = scans_micro(i,j)%tb_fov (1,1)  !166(H)
     tb_orbit(3,i,j) = 0.5 * (scans_micro(i,j)%tb_fov(3,2) + scans_micro(i,j)%tb_fov(4,2) )  !183pm3 (V)
     tb_orbit(4,i,j) = 0.5 * (scans_micro(i,j)%tb_fov(2,2) + scans_micro(i,j)%tb_fov(5,2) )  !183pm7 (V)
  enddo ; enddo

 case default  
  call stop_sdsu(myrank,'MSG write_1BGMI_HDF5:  This routine is for GMI channels only')
 end select

  dsetname = tag//'_Tb'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, tb_orbit, dims3d, error)
  call h5dclose_f(dataset_id, error)

!
! close 1st group
!
 call h5gclose_f(group_1st_id, error)  ! close the group (S#)

!
! Terminate access to the file.
!
 call h5fclose_f(file_id, error)
!
! Close FORTRAN interface.
!
 call h5close_f(error)

#else

 if(masterproc) print*,'MSG write_1BGMI_HDF5: no HDF option' 

#endif

 return
 end subroutine write_1BGMI_HDF5

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_2A12_HDF5
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  This routine output GPM 2A12-like precipitaion and related parameter outputs in HDF5 format.  
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------

!
! Dimension definitions
!
 integer :: nscan                 ! Number of scans in the granule
 integer,parameter :: nchan1 = 9  ! Number of channels in Swath 1
 integer,parameter :: nchan2 = 4  ! Number of channels in Swath 2
 integer :: nchan
 integer,parameter :: npix1 = 202  ! Number of channels in Swath 1
 integer,parameter :: npix2 = 202  ! Number of channels in Swath 2
 integer :: npix ! Number of pixels in Swath 1 or 2
 integer,parameter :: ncolds1 = 9 ! Maximum number of cold samples in Swath 1
 integer,parameter :: ncolds2 = 9 ! Maximum number of cold samples in Swath 2
 integer,parameter :: nhots1  = 30 ! Maximum number of hot samples in Swath 1
 integer,parameter :: nhots2  = 25 ! Maximum number of hot samples in Swath 2
 integer,parameter :: ntherm  = 11 ! Number of hot load thermisters

#if HDF >= 1

!
! HDF5 parameters
!
 character(len=100) :: groupdata ! group/data character name
 integer(hid_t) :: file_id       ! File identifier
 integer(hid_t) :: dataset_id    ! Dataset identifier
 integer(hid_t) :: dataspace_id  ! Data space identifier
 integer(hid_t) :: group_1st_id, group_2nd_id, group_3rd_id  ! Group identifier
 integer(HSIZE_T),dimension(1) :: dims1d
 integer(HSIZE_T),dimension(2) :: dims2d
 integer(HSIZE_T),dimension(3) :: dims3d   ! Datasets dimensions
  INTEGER(hid_t) :: h5_kind_type_i ! HDF type corresponding to the specified KIND
  INTEGER(HID_T) :: attr_id       ! Attribute identifier
  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  INTEGER     ::   arank = 1                      ! Attribure rank
  INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string
  CHARACTER(LEN=80), DIMENSION(2) ::  attr_data  ! Attribute data
  CHARACTER(LEN=9) :: aname   ! Attribute name


 integer ::i, j
 integer :: error ! Error flag


!special byte parameter
! integer,parameter :: one_byte  = selected_int_kind(4) ! one Byte integer
! integer(kind=one_byte),allocatable :: temp1d_int1(:)
 integer,allocatable :: temp1d_int1(:)
 real,allocatable :: temp2d_real(:,:)
 integer ::   rank  ! Datasets rank
 character(len=100) :: group_1st, group_2nd, group_3rd !gorup name
 character(len=100) :: dsetname  !dataset name
 integer :: im       ! middle scan index

 real(sdsu_fps), allocatable :: tb_orbit(:,:,:)  ! microwave brightness temperature
 integer :: irec
 integer :: ns, ne  !nscan start and end index
 integer :: np, nf
 character(len=2) :: tag
 logical :: fexist


!this does notwork (F2003 feature)
!  h5_kind_type_i = h5kind_to_type(one_byte,H5_INTEGER_KIND)

!
! dimension/bounds
!
 ns = 1 ; ne = jmax_scan - 4
 nscan = ne-ns+1 
 nchan = nchan1
 npix  = npix1  
 tag   = 'S1'
 im    = imax_scan/2 !
 if( imax_scan == 1 ) im = 1

!
! Initialize FORTRAN interface.
!
 call h5open_f(error)

!
! Create a new output file (if present replace)
! 
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//&
 '.2A25_GPM.'//trim(overpass_tag)//trim(output_suffix)//'.h5'
 print("(a)"),'MSG write_2A25_HDF5: GPM 2A25-like output ->',trim(sdsu_io_file)
 call h5fcreate_f(trim(sdsu_io_file), H5F_ACC_TRUNC_F, file_id, error)

 !
 ! Create swath (1st layer) group in the root group using absolute name.
 !
 group_1st = '/Swath'
 call h5gcreate_f(file_id, trim(group_1st), group_1st_id, error)

 !
 ! Store ScanTime group
 !
 groupdata = tag//'_ScanTime'
 call h5gcreate_f(group_1st_id, trim(groupdata), group_2nd_id, error) ! Create group (S#/S#_ScanTime)

  rank = 1  ; dims1d(1) = nscan 
  call h5screate_simple_f(rank, dims1d, dataspace_id, error)  ! Create the data space (nscan)

  dsetname = tag//'_Year'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%yyyy, dims1d, error)  

!  CALL h5screate_simple_f(arank, adims, aspace_id, error)
!  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error) 
!  aname = 'unit year'
!  call h5acreate_f(dataset_id, aname, atype_id, aspace_id, attr_id, error)
!  call h5awrite_f(attr_id, atype_id, attr_data, dims1d, error)
!  call h5aclose_f(attr_id, error)

  call h5dclose_f(dataset_id, error)


!following does not work --> one byte integer does not much h5write_f
  if( allocated(temp1d_int1) ) deallocate(temp1d_int1)
  allocate( temp1d_int1(ns:ne) )
  temp1d_int1(ns:ne) =  scans(im,ns:ne)%time%mm

  dsetname = tag//'_Month'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%mm , dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_DayOfMonth'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%dd, dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_Hour'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%hh, dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_Minute'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%mn, dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_Second'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, scans(im,ns:ne)%time%ss, dims1d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_DayOfYear'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, scans(im,ns:ne)%time%dayofyear, dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_SecondOfDay'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, scans(im,ns:ne)%time%secofday, dims1d, error)  
  call h5dclose_f(dataset_id, error)

 call h5sclose_f(dataspace_id, error)  ! close the data space (nscan)
 call h5gclose_f(group_2nd_id, error)  ! close the group (S#/S#_ScanTime)

!
! Store FOV lat lon
!
 rank = 2  ; dims2d(1) = npix ;  dims2d(2) = nscan
 call h5screate_simple_f(rank, dims2d, dataspace_id, error)  ! Create the data space (npix x nscan)

  dsetname = tag//'_Latitude'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(:,ns:ne)%latlon_fov(1)/degrad), dims2d, error)  
  call h5dclose_f(dataset_id, error)
 
  dsetname = tag//'_Longitude'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(:,ns:ne)%latlon_fov(2)/degrad), dims2d, error) 
  call h5dclose_f(dataset_id, error)

 call h5sclose_f(dataspace_id, error)  ! close the data space (npix x nscan)

 !
 ! navigation group
 !
 groupdata = tag//'_navigation'
 call h5gcreate_f(group_1st_id, trim(groupdata), group_2nd_id, error) ! Create group (/S#/S#_navigation)

  rank = 2  ; dims2d(1) = 3 ;  dims2d(2) = nscan
  call h5screate_simple_f(rank, dims2d, dataspace_id, error)  ! Create the data space (3xnscan)

  if( allocated(temp2d_real) ) deallocate(temp2d_real)
  allocate( temp2d_real(3,ns:ne) ) 
  do i = 1,3 ; do j = ns,ne
     temp2d_real(i,j) = real(scans(im,j)%vector_sat(i))*1000. ![m]
  enddo ; enddo

  dsetname = tag//'_scPos'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, temp2d_real, dims2d, error)  
  call h5dclose_f(dataset_id, error)

  do i = 1,3 ; do j = ns,ne
     temp2d_real(i,j) = real(scans(im,j)%vel_sat(i))*1000. ![m/sec]
  enddo ; enddo
  dsetname = tag//'_scVel'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, temp2d_real, dims2d, error)  
  call h5dclose_f(dataset_id, error)

  call h5sclose_f(dataspace_id, error)  ! close the data space (3xnscan)


  rank = 1  ; dims1d(1) = nscan
  call h5screate_simple_f(rank, dims1d, dataspace_id, error)  ! Create the data space (nscan)

  dsetname = tag//'_scLat'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(im,ns:ne)%latlon_sat(1)/degrad), dims1d, error)  
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_scLon'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(im,ns:ne)%latlon_sat(2)/degrad), dims1d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_scAlt'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, real(scans(im,ns:ne)%alt_sat*1.d3), dims1d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = tag//'_timeMidScan'
  call h5dcreate_f(group_2nd_id, trim(dsetname), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, scans(im,ns:ne)%GPSatomic, dims1d, error)
  call h5dclose_f(dataset_id, error)

 call h5sclose_f(dataspace_id, error)  ! close the data space (nscan)
 call h5gclose_f(group_2nd_id, error)  ! close the group (/S#/S#_navigation)

!
! various parameters
!
! 1 - surface rainfall rate [mm/hr]
! 2 - cloud liquid water path [kg/m2]
! 3 - toal rain water path [kg/m2]
! 4 - cloud ice water path [kg/m2]
! 5 - freezing height [m]
! 6 - land fraciotn [-]
! 7 - skin temp  [K]
! 8 - vegetation cover [-]
! 9 - soil moisture [-]
! 10 - wind speed [m/s]

 rank = 2  ; dims2d(1) = npix ;  dims2d(2) = nscan
 call h5screate_simple_f(rank, dims2d, dataspace_id, error)  ! Create the data space (npix x nscan)

  dsetname = 'SurfacePrecipitation'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(1), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'CloudWaterPath'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(2), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'RainWaterPath'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(3), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'IceWaterPath'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(4), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'freezingHeight'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(5), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'landFraction'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(6), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'surfaceTemperature'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(7), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'vegetationCover'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(8), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'soilMoisture'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(9), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'windSpeed'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(10), dims2d, error)
  call h5dclose_f(dataset_id, error)

  dsetname = 'h2ocolumn'
  call h5dcreate_f(group_1st_id, trim(dsetname), H5T_NATIVE_REAL, dataspace_id, dataset_id, error)
  call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, scans_L2(:,ns:ne)%params(11), dims2d, error)
  call h5dclose_f(dataset_id, error)



 call h5sclose_f(dataspace_id, error)  ! close the data space (npix x nscan)

!
! close 1st group
!
 call h5gclose_f(group_1st_id, error)  ! close the group (S#)

!
! Terminate access to the file.
!
 call h5fclose_f(file_id, error)
!
! Close FORTRAN interface.
!
 call h5close_f(error)

#else

 if(masterproc) print*,'MSG write_2A12_HDF5: no HDF library'

#endif

 return
 end subroutine write_2A12_HDF5

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_bundle_2d_micro
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Bundle L2 2D parameter into one 2D file.   
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------

 integer :: i,j 

!
!  - L2bundle2d index -
!
! 1 - surface rainfall rate [mm/hr]
! 2 - cloud liquid water path [kg/m2]
! 3 - toal rain water path [kg/m2]
! 4 - cloud ice water path [kg/m2]
! 5 - freezing height [m]
! 6 - land fraciotn [-]
! 7 - skin temp  [K]
! 8 - vegetation cover [-]
! 9 - soil moisture [-]
! 10 - wind speed [m/s]
!

 if( .not. allocated(L2bundle2d) ) allocate( L2bundle2d(mxgridx,mxgridy,nparam_L2) )


 do j = myj_start, myj_end ; do i = myi_start, myi_end

   L2bundle2d(i,j,1) = surface(i,j)%rain_rate_liq + surface(i,j)%rain_rate_ice    ! total surface rain [mm/hr]

   L2bundle2d(i,j,2) = surface(i,j)%path_cloud_liq  !cloud liquid [kg/m2]
   L2bundle2d(i,j,3) = surface(i,j)%path_rain_liq + surface(i,j)%path_rain_ice  !rain liquid  [kg/m2]
   L2bundle2d(i,j,4) = surface(i,j)%path_cloud_ice  !cloud ice    [kg/m2]
   L2bundle2d(i,j,5) = surface(i,j)%hgt_0degC       !0degC isotherm layer height [m]

   if(surface(i,j)%iland == 1) L2bundle2d(i,j,6) = 1.e0 !land
   if(surface(i,j)%iland == 2) L2bundle2d(i,j,6) = 0.e0 !ocean

   L2bundle2d(i,j,7)  = surface(i,j)%t_skin   ! skin temp  [K]
   L2bundle2d(i,j,8)  = surface(i,j)%frac_veg ! vegetation cover [-]
   L2bundle2d(i,j,9)  = surface(i,j)%h2o_soil ! soil moisture [-]
   L2bundle2d(i,j,10) = surface(i,j)%u10m     ! wind speed [m/s]
   L2bundle2d(i,j,11) = surface(i,j)%h2o_col  ! column water vapor [kg/m2]

 enddo ; enddo


#if MPI == 2
if(numproc_tot>1) then
  call mpi_sdsu_communicate( 'TO_ALL', L2bundle2d  )
endif
#endif


 return
 end subroutine get_bundle_2d_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine get_bundle_2d_visir
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Bundle L2 2D parameter into one 2D file for visir simulator. 
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------

 integer :: i,j,k
 real :: tot_path

!
!  - L2bundle2d index (visible-ir simulator) -
!

 if( .not. allocated(L2bundle2d) ) allocate( L2bundle2d(mxgridx,mxgridy,nparam_L2) )


 do j = myj_start, myj_end ; do i = myi_start, myi_end

   if( surface_CRM(i,j)%within_swath ) then 

   L2bundle2d(i,j,1) = surface(i,j)%rain_rate_liq + surface(i,j)%rain_rate_ice    ! total surface rain [mm/hr]
   L2bundle2d(i,j,2) = surface(i,j)%path_cloud_liq  !cloud liquid path [kg/m2]
   L2bundle2d(i,j,3) = surface(i,j)%path_rain_liq + surface(i,j)%path_rain_ice  !rain liquid  [kg/m2]
   L2bundle2d(i,j,4) = surface(i,j)%path_cloud_ice  !cloud ice    path [kg/m2]

   tot_path = surface(i,j)%path_rain_liq + surface(i,j)%path_rain_ice + &
              surface(i,j)%path_cloud_liq + surface(i,j)%path_cloud_ice
   if(tot_path > 1.e-5) then
           L2bundle2d(i,j,5) = 1.e0 ! cloud mask
   else
           L2bundle2d(i,j,5) = 0.
   endif

   if(surface(i,j)%iland == 1) L2bundle2d(i,j,6) = 1.e0 !land
   if(surface(i,j)%iland == 2) L2bundle2d(i,j,6) = 0.e0 !ocean

!toshi- temporal for publication....
!   do k = mxlyr, 1, -1
!      if(q_gen(i,j,k)%cloud + q_gen(i,j,k)%ice  > 1.e-5 ) then
!         L2bundle2d(i,j,7)  = atmos(i,j,k)%t_air
!         exit 
!      endif
!      L2bundle2d(i,j,7)  = surface(i,j)%t_skin 
!   enddo

!   do k = mxlyr, 1, -1
!      if(q_gen(i,j,k)%cloud + q_gen(i,j,k)%ice  > 1.e-4 ) then
!         L2bundle2d(i,j,8)  = atmos(i,j,k)%t_air
!         exit 
!      endif
!      L2bundle2d(i,j,8)  = surface(i,j)%t_skin
!   enddo

!   do k = mxlyr, 1, -1
!      if(q_gen(i,j,k)%cloud + q_gen(i,j,k)%ice  > 1.e-3 ) then
!         L2bundle2d(i,j,9)  = atmos(i,j,k)%t_air 
!          exit 
!      endif
!      L2bundle2d(i,j,9)  = surface(i,j)%t_skin
!   enddo

!   do k = mxlyr, 1, -1
!      if(q_gen(i,j,k)%cloud + q_gen(i,j,k)%ice  > 1.e-2 ) then
!         L2bundle2d(i,j,10)  = atmos(i,j,k)%t_air
!         exit 
!      endif
!      L2bundle2d(i,j,10)  = surface(i,j)%t_skin
!   enddo

   L2bundle2d(i,j,7)  = surface(i,j)%t_skin   ! skin temp  [K]
   L2bundle2d(i,j,8)  = surface(i,j)%frac_veg ! vegetation cover [-] 
   L2bundle2d(i,j,9)  = surface(i,j)%albedo   ! broadband albedo [-]
   L2bundle2d(i,j,10) = surface_CRM(i,j)%umu_local   ! cosine of local incident angle [-]
   L2bundle2d(i,j,11) = surface(i,j)%h2o_col  ! column water vapor [kg/m2]

   else
    L2bundle2d(i,j,:)  = undefined
   endif

 enddo ; enddo


#if MPI == 2
if(numproc_tot>1) then
  call mpi_sdsu_communicate( 'TO_ALL', L2bundle2d  )
endif
#endif


 return
 end subroutine get_bundle_2d_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine get_bundle_3d_lidar
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  Bundle L2 parameter into one 3D file.   
!
! History:
! 02/2012  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------

 integer :: i,j 

!
!  - L2bundle3d index (LIDAR simulator) -
!
! 1 - air temperature [K]
! 2 - relative humidity [%]
! 3 - cloud liquid mixing ratio [g/m3]
! 4 - cloud ice mixing ratio [g/m3]
! 5 - sulfate mixing ratio [ug/m3]
! 6 - black carbon mixing ratio [ug/m3]
! 7 - organic carbon mixing ratio [ug/m3]
! 8 - sea-salt mixing ratio [ug/m3]
! 9 - dust mixing ratio [ug/m3]

 if( .not. allocated(L2bundle3d) ) allocate( L2bundle3d(mxgridx,mxgridy,mxlyr, nparam_L2_lidar) )
 if( .not. allocated(domain_stag_hgt) ) allocate( domain_stag_hgt(mxgridx,mxgridy,0:mxlyr) )
 if( .not. allocated(domain_dhgt    ) ) allocate( domain_dhgt    (mxgridx,mxgridy,1:mxlyr) )


 do j = myj_start, myj_end ; do i = myi_start, myi_end

   domain_stag_hgt(i,j,0:mxlyr) = atmos_stag(i,j,0:mxlyr)%hgt

   domain_dhgt    (i,j,1:mxlyr) = atmos(i,j,1:mxlyr)%dhgt

   L2bundle3d(i,j,:,1) = atmos(i,j,:)%t_air   ! air temperature [K]

   L2bundle3d(i,j,:,2) = atmos(i,j,:)%rh  !relative humidity [%]

   !
   ! column-integrated condensate  [kg/m2] and instanteneous surface rate [mm/hr]
   !
   mic_select0: select case(type_microphysics)
   case('GEN')
    L2bundle3d(i,j,:,3)   = q_gen(i,j,:)%cloud   
    L2bundle3d(i,j,:,4)   = q_gen(i,j,:)%ice     
   case('RAMS')
    L2bundle3d(i,j,:,3)   = q_rams(i,j,:)%cloud1 + q_rams(i,j,:)%cloud2
    L2bundle3d(i,j,:,4)   = q_rams(i,j,:)%ice1  + q_rams(i,j,:)%ice2
   case('SBM')
    L2bundle3d(i,j,:,3)     = q_sbm(i,j,:)%liq    
    L2bundle3d(i,j,:,4)     = q_sbm(i,j,:)%ice_col &
                            + q_sbm(i,j,:)%ice_pla &
                            + q_sbm(i,j,:)%ice_den
   case default
   end select mic_select0
   if(account_aerosol) then
    select case(trim(aerosol_microphysics))
    case('GOCART')
    L2bundle3d(i,j,:,5)  = q_gocart(i,j,:)%so4 * 1.e+6  !sulfate [ug/m3]
    L2bundle3d(i,j,:,6)  = q_gocart(i,j,:)%blc * 1.e+6  !black carbon [ug/m3]
    L2bundle3d(i,j,:,7)  = (q_gocart(i,j,:)%ocn + q_gocart(i,j,:)%och) * 1.e+6  !organic carbon [ug/m3]
    L2bundle3d(i,j,:,8)  = (q_gocart(i,j,:)%ssa + q_gocart(i,j,:)%ssc) * 1.e+6  !sea salt [ug/m3]
    L2bundle3d(i,j,:,9) = (   q_gocart(i,j,:)%du1 &
                            + q_gocart(i,j,:)%du2 &
                            + q_gocart(i,j,:)%du3 &
                            + q_gocart(i,j,:)%du4 &
                            + q_gocart(i,j,:)%du5 &
                            + q_gocart(i,j,:)%du6 &
                            + q_gocart(i,j,:)%du7 &
                            + q_gocart(i,j,:)%du8 ) * 1.e+6 !dust [ug/m3]

    case('UCD')

    L2bundle3d(i,j,:,5) = q_ucd(i,j,:)%du1 * 1.e+6 !unit conversion [ug/m3]
    L2bundle3d(i,j,:,6) = q_ucd(i,j,:)%du2 * 1.e+6 
    L2bundle3d(i,j,:,7) = q_ucd(i,j,:)%du3 * 1.e+6
    L2bundle3d(i,j,:,8) = q_ucd(i,j,:)%du4 * 1.e+6
    L2bundle3d(i,j,:,9) = q_ucd(i,j,:)%du5 * 1.e+6

    case default
    end select

  else
    L2bundle3d(i,j,:,5:9) = 0.
  endif
 enddo ; enddo


#if MPI == 2
if(numproc_tot>1) then

  call mpi_sdsu_communicate( 'TO_MASTER', L2bundle3d  )

  call mpi_sdsu_communicate( 'TO_MASTER', domain_stag_hgt(:,:,0:mxlyr)  )

  call mpi_sdsu_communicate( 'TO_MASTER', domain_dhgt )
endif
#endif

!
! deallocate L2 bundle in slave  
!
 if( .not. masterproc ) then
     deallocate( L2bundle3d, domain_stag_hgt, domain_dhgt )
 endif


 return
 end subroutine get_bundle_3d_lidar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine write_out_micro_scan_image
  implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out microwave EFOV-shape brightness temperature image (V & H polarization) in GrADS format. 
! for visualization purpose. 
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
 integer,parameter :: io = 701
 character(len=3) :: fff
 character(len=6) :: tb_char

!
! Computed Tbs in binary format
!
  sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                 trim(micro_sensor)//'.'//trim(overpass_tag)//trim(output_suffix)//'_fovimage.bin'
  print("(a)"),' ->',trim(sdsu_io_file)

  open(unit= io, file = sdsu_io_file, access = 'direct', &
       status = 'replace', recl = mxgridx*mxgridy*4)

 irec = 1
 select case(trim(scan_type_micro))
 case('GMI_LF','TMI_LF')

  do l=1,mxfreq_micro !frequency loop
     do m=2,1,-1  !polarization loop
       if(l==3 .and. m==1) cycle  !no 24G horizontal channel
        write(io,rec=irec) ((tb_out_fov(i,j,l,m),i=1,mxgridx),j=1,mxgridy)  !microwave Tb

        irec = irec + 1
     enddo
  enddo

 case('GMI_HF')

    write(io,rec=1) ((tb_out_fov(i,j,1,2),i=1,mxgridx),j=1,mxgridy)  !166 V
    write(io,rec=2) ((tb_out_fov(i,j,1,1),i=1,mxgridx),j=1,mxgridy)  !166 H

    out_domain2d(:,:) =0.5e0 * ( tb_out_fov(:,:,3,2) + tb_out_fov(:,:,4,2) )
    write(io,rec=3) ((out_domain2d(i,j),i=1,mxgridx),j=1,mxgridy)  !183pm3 (V)

    out_domain2d(:,:) =0.5e0 * ( tb_out_fov(:,:,2,2) + tb_out_fov(:,:,5,2) )
    write(io,rec=4) ((out_domain2d(i,j),i=1,mxgridx),j=1,mxgridy)  !183pm7 (V)

 case('SSMIS')

  do l=1,5 !frequency loop
     do m=2,1,-1  !polarization loop
       if(l==3 .and. m==1) cycle  !no 24G horizontal channel
       if(l==5 .and. m==2) cycle  !no 150G vertical channel

        write(io,rec=irec) ((tb_out_fov(i,j,l,m),i=1,mxgridx),j=1,mxgridy)  !microwave Tb

        irec = irec + 1
     enddo
  enddo

  out_domain2d(:,:) =0.5e0 * ( tb_out_fov(:,:,7,1) + tb_out_fov(:,:,8,1) )
  write(io,rec=irec) ((out_domain2d(i,j),i=1,mxgridx),j=1,mxgridy)  !183.3pm3 (H)

  irec = irec + 1

  out_domain2d(:,:) =0.5e0 * ( tb_out_fov(:,:,6,1) + tb_out_fov(:,:,9,1) )
  write(io,rec=irec) ((out_domain2d(i,j),i=1,mxgridx),j=1,mxgridy)  !183.3pm6.6 (H)

 case default  ! TMI_HF, AMSR_E, GMI_LF_37, SSMIS_LF

  do l=1,mxfreq_micro !frequency loop
     do m=2,1,-1  !polarization loop
        write(io,rec=irec) ((tb_out_fov(i,j,l,m),i=1,mxgridx),j=1,mxgridy)  !microwave Tb
        irec = irec + 1
     enddo
  enddo

 end select

 close (io)


!
! Write grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                   trim(micro_sensor)//'.'//trim(overpass_tag)//trim(output_suffix)//'_fovimage.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.'//trim(micro_sensor)// &
                   '.'//trim(overpass_tag)//trim(output_suffix)//'_fovimage.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(micro_sensor) //' microwave Tb'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'

 select case(trim(scan_type_micro))
 case('GMI_LF','TMI_LF')
    write(io,*) 'vars ',mxfreq_micro*2-1

    do l=1,mxfreq_micro !frequency loop
       do m=2,1,-1  !polarization loop
       if(l==3 .and. m==1) cycle  !no 24G horizontal channel

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

 case('SSMIS')

    write(io,*) 'vars ', 10

    do l=1,5 !frequency loop
       do m=2,1,-1  !polarization loop
       if(l==2 .and. m==1) cycle  !no 22.235G horizontal channel (vertical only)
       if(l==5 .and. m==2) cycle  !no 150G vertical channel

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
    write(io,*) 'tb183pm3  0  0  Tb 183.3pm3(H) [K]'
    write(io,*) 'tb183pm6  0  0  Tb 183.3pm6.6(H) [K]'


 case('GMI_HF')
    write(io,*) 'vars ',4

     write(io,*) 'tb166v    0  0  Tb 165(V) [K]'
     write(io,*) 'tb166h    0  0  Tb 165(H) [K]'
     write(io,*) 'tb183pm3  0  0  Tb 183pm3(V) [K]'
     write(io,*) 'tb183pm7  0  0  Tb 183pm7(V) [K]'

 case default  ! TMI_HF, AMSR_E

! default
    write(io,*) 'vars ',mxfreq_micro*2

    do l=1,mxfreq_micro !frequency loop
       do m=2,1,-1  !polarization loop

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

 end select

    write(io,*) 'endvars '
    close(io)

 endif

 return
 end subroutine write_out_micro_scan_image

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  subroutine write_orbital_micro_latlon
  implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out microwave brightness temperature (V & H polarization) in latlon format in GrADS format. 
! 
! History:
! 01/2012  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,l,m,irec
 integer,parameter :: io = 702
 character(len=3) :: fff
 character(len=6) :: tb_char
 integer :: nf, np
 real,allocatable :: temp_scans(:,:)

!
! allocate temporal array
!
 if( .not. allocated(temp_scans) )  allocate( temp_scans(1:imax_scan,js_scan:je_scan) )


!
! set up output file name
!
 sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                trim(micro_sensor)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.latlon.bin'
 print("(a)"),'output ->',trim(sdsu_io_file)

!
! open output file
!
 open(unit= io, file = sdsu_io_file, access = 'direct', &
      status = 'replace', recl = imax_latlon*jmax_latlon*4)

 irec = 1

 select case(trim(scan_type_micro))
 case('GMI_LF','TMI_LF')

  do nf=1,mxfreq_micro !frequency loop
     do np=2,1,-1  !polarization loop
        if(nf==3 .and. np==1) cycle  ! no 24G(H) channel

           call dump2d_latlon( scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (nf,np) ,&
                               real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                               real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                               io, irec ) 

     enddo
  enddo

 case('GMI_HF')

    call dump2d_latlon( scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (1,2) ,&  !166 V
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                        io, irec )

    call dump2d_latlon( scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (1,1) ,&   !166 H 
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                        io, irec )

    temp_scans(1:imax_scan,js_scan:je_scan) = &  ! 183pm3 (V)
                        0.5e0 * (  scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (3,2) &
                                 + scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (4,2) )
    call dump2d_latlon( temp_scans(1:imax_scan,js_scan:je_scan) ,&  
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                        io, irec )

    temp_scans(1:imax_scan,js_scan:je_scan) = &  ! 183pm7 (V)
                        0.5e0 * (  scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (2,2) &
                                 + scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (5,2) )
    call dump2d_latlon( temp_scans(1:imax_scan,js_scan:je_scan) ,&   
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                        io, irec )

 case('SSMIS')

  do nf=1,5 !frequency loop
     do np=2,1,-1  !polarization loop
       if(nf==2 .and. np==1) cycle  !no 22.235G horizontal channel
       if(nf==5 .and. np==2) cycle  !no 150G vertical channel

           call dump2d_latlon( scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (nf,np) ,&
                               real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                               real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                               io, irec )

     enddo
  enddo

    temp_scans(1:imax_scan,js_scan:je_scan) = &  ! 183pm3 (H)
                        0.5e0 * (  scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (7,1) &
                                 + scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (8,1) )
    call dump2d_latlon( temp_scans(1:imax_scan,js_scan:je_scan) ,&
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                        io, irec )

    temp_scans(1:imax_scan,js_scan:je_scan) = &  ! 183pm6.6 (H)
                        0.5e0 * (  scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (6,1) &
                                 + scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (9,1) )
    call dump2d_latlon( temp_scans(1:imax_scan,js_scan:je_scan) ,&
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                        real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                        io, irec )


 case default  ! TMI_HF, AMSR_E, GMI_LF_37, SSMIS_LF

  do nf=1,mxfreq_micro !frequency loop
     do np=2,1,-1  !polarization loop

           call dump2d_latlon( scans_micro(1:imax_scan,js_scan:je_scan)%tb_fov (nf,np) ,&
                               real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                               real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                               io, irec )

     enddo
  enddo

 end select

!
! L2 parameters
!

  do np=1,nparam_L2 !param Loop

     call dump2d_latlon( scans_L2(1:imax_scan,js_scan:je_scan)%params(np) ,&
                         real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(1)/degrad) ,&
                         real(scans(1:imax_scan,js_scan:je_scan)%latlon_fov(2)/degrad) ,&
                         io, irec )

  enddo

 write(io,rec=irec) latlon_grid%area

 close (io)  !close binary file


!
! Write grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                   trim(micro_sensor)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.latlon.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.'//trim(micro_sensor)// &
                   '.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.latlon.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(micro_sensor) //' microwave orbital L1 and L2 data'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',imax_latlon,' linear ',minlon, res_latlon
    write(io,*) 'ydef   ',jmax_latlon,' linear ',minlat, res_latlon
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'

 select case(trim(scan_type_micro))
 case('GMI_LF','TMI_LF')
    write(io,*) 'vars ',mxfreq_micro*2-1 + nparam_L2 + 1

    do l=1,mxfreq_micro !frequency loop
       do m=2,1,-1  !polarization loop
       if(l==3 .and. m==1) cycle  !no 24G horizontal channel

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

 case('SSMIS')

    write(io,*) 'vars ', 10+nparam_L2 + 1

    do l=1,5 !frequency loop
       do m=2,1,-1  !polarization loop
       if(l==2 .and. m==1) cycle  !no 22.235G horizontal channel (vertical only)
       if(l==5 .and. m==2) cycle  !no 150G vertical channel

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
    write(io,*) 'tb183pm3  0  0  Tb 183.3pm3(H) [K]'
    write(io,*) 'tb183pm6  0  0  Tb 183.3pm6.6(H) [K]'


 case('GMI_HF')
    write(io,*) 'vars ',4+nparam_L2 +1

     write(io,*) 'tb166v    0  0  Tb 165(V) [K]'
     write(io,*) 'tb166h    0  0  Tb 165(H) [K]'
     write(io,*) 'tb183pm3  0  0  Tb 183pm3(V) [K]'
     write(io,*) 'tb183pm7  0  0  Tb 183pm7(V) [K]'

 case default  ! TMI_HF, AMSR_E

! default
    write(io,*) 'vars ',mxfreq_micro*2+nparam_L2 + 1

    do l=1,mxfreq_micro !frequency loop
       do m=2,1,-1  !polarization loop

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

 end select


    write(io,*) 'rain    0  0  surface rainfall rate [mm/hr]'
    write(io,*) 'clwp    0  0  cloud liquid water path [kg/m2]'
    write(io,*) 'trwp    0  0  total rain water path [kg/m2]'
    write(io,*) 'ciwp    0  0  cloud ice water path [kg/m2]'
    write(io,*) 'fhgt    0  0  freezing height [m]'
    write(io,*) 'lfrac   0  0  land fraction [-]'
    write(io,*) 'skin    0  0  skin temperature  [K]'
    write(io,*) 'vfrac   0  0  vegetation cover [%]'
    write(io,*) 'smo     0  0  soil moisture [-]'
    write(io,*) 'ws      0  0  wind speed [m/s]'
    write(io,*) 'h2ocol  0  0  column water vapor [kg/m2]'
    write(io,*) 'area    0  0  area of pixel [km2]'
    write(io,*) 'endvars '
    close(io)

 endif

 return
 end subroutine write_orbital_micro_latlon

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_orbital_visir_latlon
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out orbital data of visible IR imager radiance (or brightness temperautre for wave > 10micon)
! in lat-lon cordinate in GrADS format. 
! 
! History:
! 02/2012  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,iw,nw,irec
 real(sdsu_fps) :: wave_number    !ending wavenumber [1/cm]
 integer :: iw_r, np
 real(sdsu_fps) :: wl
 integer,parameter :: io = 105
 real(sdsu_fps),parameter :: tb_or_rad = 10.e0
 character(len=2) :: ww

 print*,'MSG write_orbital_visir_latlon: start dumping output'


!
! write radiances in binary format
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
               trim(visir_sensor)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.latlon.bin'
 print("(a)"),' ->',trim(sdsu_io_file)
 open( unit= io, file = sdsu_io_file, access='direct', status='replace', recl=imax_latlon*jmax_latlon*4 )

 irec = 1
 do nw = 1, mxwavel

    call dump2d_latlon( radiance_out(1:mxgridx,1:mxgridy,nw)  ,&
                        real(surface_CRM(1:mxgridx,1:mxgridy)%latlon_grid(1)/degrad) , &
                        real(surface_CRM(1:mxgridx,1:mxgridy)%latlon_grid(2)/degrad) , &
                        io, irec )

 enddo



  do np=1,nparam_L2 !param Loop

     call dump2d_latlon( L2bundle2d(1:mxgridx,1:mxgridy,np) ,&
                        real(surface_CRM(1:mxgridx,1:mxgridy)%latlon_grid(1)/degrad) , &
                        real(surface_CRM(1:mxgridx,1:mxgridy)%latlon_grid(2)/degrad) , &
                        io, irec )

  enddo


 write(io,rec=irec) latlon_grid%area
 close (io)

!
! write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(visir_sensor)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.latlon.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//trim(visir_sensor) &
                  //'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.latlon.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(visir_sensor) //' radiance/Tb'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',imax_latlon,' linear ',minlon, res_latlon
    write(io,*) 'ydef   ',jmax_latlon,' linear ',minlat, res_latlon
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxwavel + nparam_L2 + 1

    do nw=1,mxwavel !frequency loop
       write(ww,"(I2.2)") nw
       if(wavel(nw) < tb_or_rad) then ! radiance output
          write(io,*) 'cha'//ww//' 0  0  radiance at '//trim(nch_wavel(nw))//' [W/m2/str/micron]'
       else ! Tb output
          write(io,*) 'cha'//ww//' 0  0  Tb at '//trim(nch_wavel(nw))//' [K]'
       endif
    enddo

    write(io,*) 'rain    0  0  surface rainfall rate [mm/hr]'
    write(io,*) 'clwp    0  0  cloud liquid water path [kg/m2]'
    write(io,*) 'trwp    0  0  total rain water path [kg/m2]'
    write(io,*) 'ciwp    0  0  cloud ice water path [kg/m2]'
    write(io,*) 'cmask   0  0  cloud mask [-]'
    write(io,*) 'lfrac   0  0  land fraction [-]'
    write(io,*) 'skin    0  0  skin temperature  [K]'
    write(io,*) 'vfrac   0  0  vegetation cover [%]'
    write(io,*) 'alb     0  0  broadband albedo [-]'
    write(io,*) 'umu     0  0  cosine of incident angle [-]'
    write(io,*) 'h2ocol  0  0  column water vapor [kg/m2]'
    write(io,*) 'area    0  0  area of pixel [km2]'
    write(io,*) 'endvars '
    close(io)
 endif



 return
 end subroutine write_orbital_visir_latlon

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_orbital_visir
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
 integer i,j,iw,nw,irec,np
 real(sdsu_fps) :: wave_number    !ending wavenumber [1/cm]
 integer :: iw_r
 real(sdsu_fps) :: wl
 integer,parameter :: io = 105
 real(sdsu_fps),parameter :: tb_or_rad = 10.e0
 character(len=2) :: ww

!
! write radiances in binary format
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
               trim(visir_sensor)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.bin'
  print("(a)"),' ->',trim(sdsu_io_file)
  open( unit= io, file = sdsu_io_file, access='direct', status='replace', recl=mxgridx*mxgridy*4)

  irec = 1
  do nw = 1, mxwavel
     write(io,rec=irec) ((radiance_out(i,j,nw),i=1,mxgridx), j=1,mxgridy)  ! Visible Radiance or IR Tb
     irec = irec + 1
  enddo

  do np=1,nparam_L2 !param Loop
     write(io,rec=irec) ((L2bundle2d(i,j,np),i=1,mxgridx), j=1,mxgridy)  ! L2 parameter output
     irec = irec + 1
  enddo

  close (io)


!
! write out grads control file
!
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//trim(visir_sensor) &
                  //'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.latlon.bin'


 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(visir_sensor)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//trim(visir_sensor) &
                  //'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.bin'
    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(visir_sensor) //' radiance/Tb'
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxwavel + nparam_L2 + 1

    do nw=1,mxwavel !frequency loop
       write(ww,"(I2.2)") nw
       if(wavel(nw) < tb_or_rad) then ! radiance output
          write(io,*) 'cha'//ww//' 0  0  radiance at '//trim(nch_wavel(nw))//' [W/m2/str/micron]'
       else ! Tb output
          write(io,*) 'cha'//ww//' 0  0  Tb at '//trim(nch_wavel(nw))//' [K]'
       endif
    enddo

    write(io,*) 'rain    0  0  surface rainfall rate [mm/hr]'
    write(io,*) 'clwp    0  0  cloud liquid water path [kg/m2]'
    write(io,*) 'trwp    0  0  total rain water path [kg/m2]'
    write(io,*) 'ciwp    0  0  cloud ice water path [kg/m2]'
    write(io,*) 'cmask   0  0  cloud mask [-]'
    write(io,*) 'lfrac   0  0  land fraction [-]'
    write(io,*) 'skin    0  0  skin temperature  [K]'
    write(io,*) 'vfrac   0  0  vegetation cover [%]'
    write(io,*) 'alb     0  0  broadband albedo [-]'
    write(io,*) 'umu     0  0  cosine of incident angle [-]'
    write(io,*) 'h2ocol  0  0  column water vapor [kg/m2]'
    write(io,*) 'area    0  0  area of pixel [km2]'
    write(io,*) 'endvars '
    close(io)
  endif

  return
 end subroutine write_orbital_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine orbital_sample_lidar
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Since Lidar instrument has much smaller foot print size (~30m) than usual CRM, so this subroutine
! essentially sample CRM-grid lidar valules for orbital sampling. 
! 
! History:
! 04/2012  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: i,j ,im ,kl,nw,ip !loop for CRM domain
 integer :: i_scan, j_scan !loop for scan sample domain
 integer :: ierr
 integer :: cnt_closest
 integer :: sample
 integer :: ic, jc
 integer :: cnt_loc(1:1)
 integer , allocatable :: i_crm(:), j_crm(:), j_scan_lidar(:)
 real(sdsu_fpd),allocatable :: distance_grid2fov(:)
 real(sdsu_fpd) :: max_dt,max_ct,crit_distance !critical distance
 real(sdsu_fpd) :: distance,lat_test,lon_test, crit_collocate
 real(sdsu_fpd) :: sum_dist , sum_vec(1:3)
 integer :: cnt, tot_cnt
 real(sdsu_fps),allocatable :: var_lidar(:)

 if(masterproc) print("(a)"),'MSG orbital_sample_lidar: start re-sample LIDAR signals along orbital array.'
!
! allocate orbital output array over a CRM domain
!
 if( .not. allocated(var_lidar) ) allocate( var_lidar(1:mxlyr_lidar) ) 

 if( allocated(lidar_orbital) ) deallocate(lidar_orbital, lidar_l2bundle )
 allocate( lidar_orbital   (js_scan:je_scan, mxlyr_lidar ,mxwavel_lidar),&
           lidar_l2bundle  (js_scan:je_scan, mxlyr_lidar ,nparam_L2_lidar ) )
 lidar_orbital%att_B       = undefined
 lidar_orbital%att_p22_B   = undefined
 lidar_orbital%lidar_sback = undefined
 lidar_orbital%lidar_ratio = undefined
 lidar_orbital%optical_depth = undefined
 lidar_l2bundle          = undefined


!
! crtical length for collocation
!
 crit_collocate = 0.5d0*sqrt( gridsize**2+gridsize**2  )  ![km]

 !
 ! CRM domain loop to calculate total count of on-nadir pixels. 
 !
 cnt = 0
  do j = 1, mxgridy ; do i = 1, mxgridx
    if( surface_CRM(i,j)%latlon_grid(1) == undefined_scan ) cycle   !this is domain memory by masterproc 
    cnt = cnt + 1
  enddo ; enddo
 tot_cnt = cnt

!
! if there is not lidar-overlapped CRM pixel, return
!
 if(tot_cnt == 0) return 

!
! allcoate distance grid
!
 if( allocated(distance_grid2fov) ) deallocate(distance_grid2fov)
 if( allocated(i_crm)             ) deallocate(i_crm)
 if( allocated(j_crm)             ) deallocate(j_crm)

 allocate( distance_grid2fov( 1:tot_cnt ) , stat=ierr )
 if (ierr /= 0) call stop_sdsu(myrank, 'MSG orbital_sample_lidar: allocation error -> Terminate program. ')
 distance_grid2fov = 1.d5 !initialize as large value

 allocate( i_crm(1:tot_cnt) ) ; i_crm = undefined_i2
 allocate( j_crm(1:tot_cnt) ) ; j_crm = undefined_i2

 !
 ! scan sample loop to get arc distance betweeen CRM grid and each sample FOV
 !
 sample =  0
 scan_loop: do j_scan = js_scan , je_scan
 im = 1

   !
   ! CRM memory loop
   !
   cnt = 0
   crm_loop: do j = 1, mxgridy ; do i = 1, mxgridx

      if( surface_CRM(i,j)%latlon_grid(1) == undefined_scan ) cycle
      cnt = cnt + 1

!      print*,'CRM memory loop cnt=',i,j,cnt, surface_CRM(i,j)%latlon_grid(1)
      call get_arc_distance( surface_CRM(i,j)%vector_grid, scans(im,j_scan)%vector_fov,&
                            distance_grid2fov(cnt) )
      i_crm(cnt) = i
      j_crm(cnt) = j 
   enddo ; enddo crm_loop

   !
   ! get the closest single CRM grid from one FOV 
   !
   cnt_loc(1:1) = minloc( distance_grid2fov(:) )
   cnt_closest = cnt_loc(1)  

   !
   ! if the closest grid cell is withing the critical distance (grid spacging), assigne value here. 
   !
   if( distance_grid2fov(cnt_closest) <= crit_collocate ) then

       ic= i_crm(cnt_closest) ; jc= j_crm(cnt_closest)  !closest CRM grid point

       lidar_orbital(j_scan,:,:)%att_B       = att_B      (ic,jc,:,:)  !attenuating backscatter [km-1 str-1]
       lidar_orbital(j_scan,:,:)%att_p22_B   = att_p22_B  (ic,jc,:,:)  !perpendicular polarization [km-1 str-1]
       lidar_orbital(j_scan,:,:)%lidar_sback = lidar_sback(ic,jc,:,:)  !non-attenuating backscattering  [km-1 str-1]
       lidar_orbital(j_scan,:,:)%lidar_ratio = lidar_ratio(ic,jc,:,:)  !lidar ratio
       lidar_orbital(j_scan,:,:)%optical_depth = optical_depth(ic,jc,:,:) !optical depth [-]

       scans(im,j_scan)%i_crm = ic
       scans(im,j_scan)%j_crm = jc

    do ip = 1, nparam_L2_lidar

       ! resample CRM-vertical L2 parameters to CALIPSO vertical level 
       call refine_range_bin( mxlyr, domain_stag_hgt(ic,jc,0:mxlyr), domain_dhgt(ic,jc,1:mxlyr), &
                              mxlyr_lidar, hgt_stag_lidar, &
                              L2bundle3d(ic,jc,1:mxlyr,ip), var_lidar(1:mxlyr_lidar) )

       do kl = 1, mxlyr_lidar
          lidar_l2bundle(j_scan,kl,ip) = var_lidar( kl ) !L2 parameters
       enddo

    enddo

    !print*,j_scan, lidar_l2bundle(j_scan,:,1)

    !
    ! vertical smoothing (~500m range regoslution) to match up with Merged_CloudSat_CALIPSO data. 
    !
    do nw=1,mxwavel_lidar
       call range_smoothing(  lidar_orbital(j_scan,:,nw)%att_B       )
       call range_smoothing(  lidar_orbital(j_scan,:,nw)%att_p22_B   )
       call range_smoothing(  lidar_orbital(j_scan,:,nw)%lidar_sback )
       call range_smoothing(  lidar_orbital(j_scan,:,nw)%lidar_ratio )
       call range_smoothing(  lidar_orbital(j_scan,:,nw)%optical_depth )
    enddo

    do ip = 1, nparam_L2_lidar
       call range_smoothing( lidar_l2bundle(j_scan,:,ip)   )
    enddo

   endif

 enddo scan_loop

 return
 end subroutine orbital_sample_lidar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_orbital_lidar
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out orbital lidar simulator output in GrADS format. 
! 
! History:
! 05/2012  Toshi Matsui@NASA GSFC ; Initial 
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,k,nw,ip
 integer :: irec
 integer :: kmax
 integer,parameter :: io = 106
 character(len=2) :: ww
 real :: od

!
! define kmax
!
 kmax = mxlyr_lidar

!
! write lidar attenuating backscatter in binary format
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(lidar_sensor)//&
               '.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.bin'
 print("(a)"),' ->',trim(sdsu_io_file)

 open(unit= io, file = sdsu_io_file, access='direct', &
      status='replace',recl=kmax*(je_scan-js_scan+1)*4)

 irec = 1
 do nw=1,mxwavel_lidar

    ! total attenuating backscattering [km-1 str-1]
    write(io,rec=irec)  ((( lidar_orbital(j,k,nw)%att_B),j=js_scan,je_scan),k=1,kmax) 
    irec = irec + 1

    ! perpendicular polariztion component  [km-1 str-1]
    write(io,rec=irec)  ((( lidar_orbital(j,k,nw)%att_p22_B),j=js_scan,je_scan),k=1,kmax)
    irec = irec + 1

    ! total backscattering  [km-1 str-1]
    write(io,rec=irec)  ((( lidar_orbital(j,k,nw)%lidar_sback),j=js_scan,je_scan),k=1,kmax)
    irec = irec + 1

    ! lidar ratio [-]
    write(io,rec=irec)  ((( lidar_orbital(j,k,nw)%lidar_ratio),j=js_scan,je_scan),k=1,kmax) 
    irec = irec + 1

    ! optical depth [-]
    write(io,rec=irec)  ((( lidar_orbital(j,k,nw)%optical_depth),j=js_scan,je_scan),k=1,kmax)
    irec = irec + 1

    ! column optical depth [-]
    do j = js_scan,je_scan
       od = 0.
       do k = 1, kmax 
          if( lidar_orbital(j,k,nw)%optical_depth /= undefined ) od = od + lidar_orbital(j,k,nw)%optical_depth
       enddo      
       lidar_orbital(j,1:kmax,nw)%optical_depth = od  !replace    
    enddo
    write(io,rec=irec)  ((( lidar_orbital(j,k,nw)%optical_depth),j=js_scan,je_scan),k=1,kmax)
    irec = irec + 1

 enddo

 do ip = 1, nparam_L2_lidar  ! Level 2 geophysical parameters
    write(io,rec=irec)  (( lidar_l2bundle(j,k,ip),j=js_scan,je_scan),k=1,kmax)    
    irec = irec + 1
 enddo


 close(io)

!
!  write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(lidar_sensor)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//&
                 trim(lidar_sensor)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(lidar_sensor)//' signal '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef    1      linear 1. 1. '
    write(io,*) 'ydef   ',je_scan-js_scan+1,' linear 1. ', stat%mean_scan_dist
    write(io,*) 'zdef   ',kmax  ,' linear 0. ',range_lidar
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',irec-1

    do nw=1,mxwavel_lidar

       write(ww,"(I2.2)") nw
       write(io,*) 'att'//ww   ,kmax, ' 0  attenuating backscatter '//trim(nch_wavel_lidar(nw))//' [km-1 str-1]'
       write(io,*) 'ptt'//ww   ,kmax, ' 0  p22 attenuating backscatter '//trim(nch_wavel_lidar(nw))//' [km-1 str-1]'
       write(io,*) 'ntt'//ww   ,kmax, ' 0  non-attenuating backscatter '//trim(nch_wavel_lidar(nw))//'  [km-1 str-1]'
       write(io,*) 'ratio'//ww ,kmax, ' 0  lidar ratio '//trim(nch_wavel_lidar(nw))//'  [-]'
       write(io,*) 'od'//ww ,kmax, ' 0  optical depth '//trim(nch_wavel_lidar(nw))//'  [-]'
       write(io,*) 'odcol'//ww ,kmax, ' 0  column optical depth '//trim(nch_wavel_lidar(nw))//'  [-]'

    enddo

    write(io,*) 'temp ',kmax,' 0  air temperature [K]'
    write(io,*) 'rh   ',kmax,' 0  relative humidity [%]'
    write(io,*) 'qc   ',kmax,' 0  cloud liq mixing ratio [g/m3]'
    write(io,*) 'qi   ',kmax,' 0  cloud ice mixing ratio [g/m3]'

    select case(trim(aerosol_microphysics))
    case('GOCART')
     write(io,*) 'so4  ',kmax,' 0  so4 mixing ratio [ug/m3]'
     write(io,*) 'blc  ',kmax,' 0  black carbon mixing ratio [ug/m3]'
     write(io,*) 'oc   ',kmax,' 0  organic carbon mixing ratio [ug/m3]'
     write(io,*) 'ss   ',kmax,' 0  sea salt mixing ratio [ug/m3]'
     write(io,*) 'du   ',kmax,' 0  dust mixing ratio [ug/m3]'
    case('UCD')
     write(io,*) 'du1  ',kmax,' 0  dust1 mixing ratio [ug/m3]'
     write(io,*) 'du2  ',kmax,' 0  dust2 mixing ratio [ug/m3]'
     write(io,*) 'du3  ',kmax,' 0  dust3 mixing ratio [ug/m3]'
     write(io,*) 'du4  ',kmax,' 0  dust4 mixing ratio [ug/m3]'
     write(io,*) 'du5  ',kmax,' 0  dust5 mixing ratio [ug/m3]'

    case default
    end select

    write(io,*) 'endvars '
    close(io)
 endif

 return
 end subroutine write_orbital_lidar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_composite_lidar 
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!   Write out lidar backscatter-color_ratio diagram in GrADS format. 
! 
! History:
! 09/2012  Toshi Matsui@NASA GSFC ; Initial 
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: j,k
 integer :: kmax
 integer,parameter :: io = 107
 integer :: irec       !record number for binary output

 real :: color_ratio   ! color ratio (1012/532) [-] 
 real :: back_scatt    ! lidar backscatter [km-1 str-1]

 integer,parameter :: maxbin_color = 15         !maximum bin for color ratio bin
 real,parameter :: bin_color(maxbin_color) = &  !color ratio bin boundary
  (/0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6/)
 integer,parameter :: maxbin_back = 14        !maximum bin for backscatter bin
 real,parameter :: bin_back(maxbin_back) = &  !backscatter bin boundary
  (/5e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2, &
    2e-2, 3e-2, 5e-2, 1e-1, 2e-1, 3e-1, 5e-1, 1e-0/)

 real,dimension(maxbin_back) :: bin_back1, bin_back2
 integer :: ibin_col, ibin_back  !histogram bin indice
 real,dimension(maxbin_back,maxbin_color) :: cnt, pd  !histrogram count and probability density

!
! define kmax
!
 kmax = mxlyr_lidar

!
! set up bin for backscatter coeff
!
 bin_back1(1) = bin_back(1)
 do ibin_back = 2, maxbin_back
     bin_back1(ibin_back) = 0.5 * (bin_back(ibin_back-1) + bin_back(ibin_back) )
 enddo

 bin_back2(maxbin_back) = bin_back(maxbin_back)
 do ibin_back = 1, maxbin_back-1
     bin_back2(ibin_back) = 0.5 * (bin_back(ibin_back+1) + bin_back(ibin_back) )
 enddo

!
! now estimate joint historgram
!
 cnt = 0.
 do j=js_scan,je_scan
! do n = ns_domain, ne_domain
    k_loop: do k = 1, kmax
       if( lidar_orbital(j,k,1)%att_B <= tiny(1.0) ) cycle 
       color_ratio =  lidar_orbital(j,k,2)%att_B / lidar_orbital(j,k,1)%att_B   !color ratio [-]
       back_scatt  =  lidar_orbital(j,k,1)%att_B   ! backscatter
       if( back_scatt < 1.e-3 ) cycle k_loop

       do ibin_back = 1, maxbin_back
          if( back_scatt > bin_back1(ibin_back) .and. back_scatt <= bin_back2(ibin_back) ) then
             do ibin_col = 1, maxbin_color
                if( color_ratio > bin_color(ibin_col) - 0.1 .and.  color_ratio <= bin_color(ibin_col) + 0.1 ) then
                    cnt(ibin_back,ibin_col)= cnt(ibin_back,ibin_col) + 1.
                endif
             enddo
          endif
       enddo

    enddo k_loop
 enddo

 pd = cnt / SUM(cnt(:,:) )  !get normalized frequency (probability density)


!
! write output 
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(lidar_sensor)//&
               '.Joint_PDF.'//trim(overpass_tag)//trim(output_suffix)//'.bin'
 print("(a)"),'write_composite_lidar: output binary ->',trim(sdsu_io_file)

 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=maxbin_back*maxbin_color*4)
 irec = 1
 write(io,rec=irec) ((pd(ibin_back,ibin_col) ,ibin_back=1,maxbin_back),ibin_col=1,maxbin_color)
 close(io)

!
! write grads control file
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(lidar_sensor)//&
               '.Joint_PDF.'//trim(overpass_tag)//trim(output_suffix)//'.ctl'
 print("(a)"),'MSG write_composite_lidar: output control file->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//trim(lidar_sensor)//&
               '.Joint_PDF.'//trim(overpass_tag)//trim(output_suffix)//'.bin'

 write(io,"(a)") 'dset  '//trim(sdsu_io_file)
 write(io,*) 'title  CALIPSO color ratio - backscatter joint PDF '
 write(io,*) 'undef   -999. '
 write(io,*) 'xdef   ',maxbin_back  ,'linear 1 1'
 write(io,*) 'ydef   ',maxbin_color ,'linear 1 1'
 write(io,*) 'zdef    1 linear 0. 1 '
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars 1 '
 write(io,*) 'pd    0 0  Frequency [-]'
 write(io,*) 'endvars  '

 return
 end subroutine write_composite_lidar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_orbital_radar( sensor_name, is, ie )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out radar reflectivity [dBZ] in GrADS format. 
! 
! History:
! 03/2012  Toshi Matsui@NASA GSFC : grads control file is automatically writen with output file. 
! 10/2007  Toshi Matsui@NASA GSFC ; FOV options
! 03/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 character(*),intent(in) :: sensor_name  !sensor name
 integer,intent(in) :: is, ie
 integer :: js,je,ks,ke,k_index
 integer :: i,j,k,l,m  !for looping
 integer :: irec  !irec
 integer,parameter :: io = 104
 integer :: kmax
 real(sdsu_fps),allocatable :: temp3d(:,:,:)
 character(len=3) :: fff
 character(len=50) :: sn  !sensor name


 js = js_scan ; je = je_scan
 ks = 1       ; ke = irange_max
  
!
! find out start vertical index and end vertical index (20bin below ellipsoid, 130bin above) 
!
 if( trim(sensor_name) == 'DPR_Ka_HS' ) then
   kmax = int( mxhgt_radar / range_radar )/2  ! DPR Ka High Sensivity (250m range bin)
 else 
   kmax = int( mxhgt_radar / range_radar )    ! bounds(imax_scan,je_scan)%ke - bounds(imax_scan,je_scan)%ks + 1 !toshii
 endif

! kmax = irange_max

 if(allocated(temp3d)) deallocate(temp3d)
 allocate( temp3d(is:ie,js:je,ks:ke) )

!
! Computed radar echoes in binary format
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                trim(sensor_name)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.bin'
 print("(a)"),' ->',trim(sdsu_io_file)

!
! Write out binary data 
!
  open( unit= io, file = sdsu_io_file, access='direct', &
        status='replace',recl=kmax*(ie-is+1)*(je-js+1)*4)
  irec = 1
  sn = sensor_name
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%gain     , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%Zt       , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%Zm       , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, real(scans_radar(is:ie,js:je,ks:ke)%hgt), bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%iwc      , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%lwc      , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%rwc      , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%rain_rate, bounds(is:ie,js:je), io, irec )

#ifdef FULL_L2_RADAR
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%press    , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%t_air    , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%sh       , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%w        , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%lfrac    , bounds(is:ie,js:je), io, irec )
  call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%bulk_den , bounds(is:ie,js:je), io, irec )
  do m = 0, 6
     call dump_radar_range( sn, scans_radar(is:ie,js:je,ks:ke)%mmt(m), bounds(is:ie,js:je), io, irec )
  enddo
#endif

  do j = js,je ; do i = is,ie
     temp3d(i,j,:) =  real(scans(i,j)%latlon_fov(1)/degrad)  !latitude [deg]
  enddo ; enddo
  call dump_radar_range( sn, temp3d(is:ie,js:je,ks:ke) , bounds(is:ie,js:je), io, irec )
  do j = js,je ; do i = is,ie
     temp3d(i,j,:) =  real(scans(i,j)%latlon_fov(2)/degrad)  !longitude [deg]
  enddo ; enddo
  call dump_radar_range( sn, temp3d(is:ie,js:je,ks:ke) , bounds(is:ie,js:je), io, irec )

  !PIA at SRT level
  do j = js,je ; do i = is,ie
     k_index = bounds(i,j)%k_srt
     if( scans_radar(i,j, k_index )%Zm == undefined ) then
        temp3d(i,j,ks:ke) = undefined
     elseif( scans_radar(i,j, k_index )%Zm == 0. ) then !if no return use tau
        temp3d(i,j,ks:ke) = 20.* scans_radar(i,j,k_index)%tau * log10(2.71828183) 
     else !if echo is significant, use dBZt-dBZm
        temp3d(i,j,ks:ke) =  10.*log10( scans_radar(i,j, k_index )%Zt/scans_radar(i,j, k_index )%Zm)
     endif
  enddo ; enddo

  call dump_radar_range( sn, temp3d(is:ie,js:je,ks:ke) , bounds(is:ie,js:je), io, irec )

  !PIA at effective level
  do j = js,je ; do i = is,ie
     k_index = bounds(i,j)%k_piae
     if( scans_radar(i,j, k_index )%Zm == undefined ) then
        temp3d(i,j,ks:ke) = undefined
     elseif( scans_radar(i,j, k_index )%Zm == 0. ) then !if no return use tau
        temp3d(i,j,ks:ke) = 20.* scans_radar(i,j,k_index)%tau * log10(2.71828183) 
     else  !if echo is significant, use dBZt-dBZm
        temp3d(i,j,ks:ke) =  10.*log10( scans_radar(i,j, k_index )%Zt/scans_radar(i,j, k_index )%Zm)
     endif

  enddo ; enddo
  call dump_radar_range( sn, temp3d(is:ie,js:je,ks:ke) , bounds(is:ie,js:je), io, irec )


  close(io)

!
!  write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
                  trim(sensor_name)//'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//trim(sensor_name)&
                  //'.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(sensor_name)//' echo '
    write(io,*) 'undef  ',undefined
    if( sensor%mx_sample == 1 ) then ! on-nadir sampling only (CPR)
        write(io,*) 'xdef   ',(ie-is+1)    ,' linear 1. 1,'
    else
       write(io,*) 'xdef   ',(ie-is+1)    ,' linear 1. ',stat%mean_fov_dist
    endif

    write(io,*) 'ydef   ',(je-js+1)    ,' linear 1. ',stat%mean_scan_dist 

    if( trim(sensor_name) == 'DPR_Ka_HS' ) then
       write(io,*) 'zdef   ',kmax  ,' linear 0. ',range_radar*2.
    else
       write(io,*) 'zdef   ',kmax  ,' linear 0. ',range_radar
    endif

    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ', irec-1
    l = 1
    write(fff,"(I3.3)") INT(freq_radar(l))
    write(io,*) 'gain'     ,kmax, ' 0  normalized gain function (<1 is bad sampling)'
    write(io,*) 'Zt'       ,kmax, ' 0  non-attenuating echo at'//trim(nch_radar(l))//' [Z]'
    write(io,*) 'Zm'       ,kmax, ' 0  attenuating echo at '//trim(nch_radar(l))//' [Z]'
    write(io,*) 'hgt'      ,kmax, ' 0  echo height [m]'
    write(io,*) 'iwc'      ,kmax, ' 0  ice water content [g/m3]'
    write(io,*) 'lwc'      ,kmax, ' 0  liquid water content [g/m3]'
    write(io,*) 'rwc'      ,kmax, ' 0  rain water content [g/m3]'
    write(io,*) 'rain'     ,kmax, ' 0  rain rate  [mm/hr]'

#ifdef FULL_L2_RADAR
    write(io,*) 'press'    ,kmax, ' 0  air pressure [hPa]'
    write(io,*) 'tair'     ,kmax, ' 0  air temperature [K]'
    write(io,*) 'sh'       ,kmax, ' 0  specific humidity [g/g]'
    write(io,*) 'w    '    ,kmax, ' 0  vertical velocity [m/s]'
    write(io,*) 'lfrac'    ,kmax, ' 0  liquid fraction [-]'
    write(io,*) 'rhoe'     ,kmax, ' 0  bulk density [g/m3]'
    write(io,*) 'mmt0'     ,kmax, ' 0  0th moment [m0/m3]'
    write(io,*) 'mmt1'     ,kmax, ' 0  1st moment [m1/m3]'
    write(io,*) 'mmt2'     ,kmax, ' 0  2nd moment [m2/m3]'
    write(io,*) 'mmt3'     ,kmax, ' 0  3rd moment [m3/m3]'
    write(io,*) 'mmt4'     ,kmax, ' 0  4th moment [m4/m3]'
    write(io,*) 'mmt5'     ,kmax, ' 0  5th moment [m5/m3]'
    write(io,*) 'mmt6'     ,kmax, ' 0  6th moment [m6/m3]'
#endif

    write(io,*) 'rlat'     ,kmax, ' 0  latitude [deg]'
    write(io,*) 'rlon'     ,kmax, ' 0  longitude [deg]'
    write(io,*) 'PIAsrt'   ,kmax, ' 0  PIAsrt  [dB]'
    write(io,*) 'PIAe'     ,kmax, ' 0  PIAe    [dB]'

    write(io,*) 'endvars '
    close(io)
 endif

 return
 end subroutine write_orbital_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 subroutine write_composite_radar
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out ISCCP-like composite from CPR reflectivity [dBZ] in GrADS format. 
! 
! History:
! 09/2012  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer :: iscan, jscan, k, n
 integer :: is,ie,js,je,ks,ke,k_index
 character(len=200) :: outputfile, outputfile_ctl
 integer,parameter :: io = 102
 integer :: k_echo_top , k_echo_botm
 real(sdsu_fps) :: hgt_echotop,hgt_echobotm ![km]
 real(sdsu_fps) :: echo_thick
 integer ::i,j, ii,jj                   !diagram bin index
 integer,parameter :: binmax = 18
 real,dimension(0:binmax,0:binmax) :: isccp_cnt,isccp_pdf !(thick, hgt)
 real(sdsu_fps),allocatable :: dBZ(:)  !dB of Z

!
! initialize
!
 isccp_cnt = 0. ; isccp_pdf = 0.0
 is = 1 ; ie = imax_scan
 js = js_scan ; je = je_scan
 ks = 1       ; ke = irange_max

 if(allocated(dBZ)) deallocate(dBZ)
 allocate(dBZ(ks:ke))

!
! now estimate joint historgram
!
 do jscan = js_scan , je_scan
 do iscan = 1, 1

   !
   ! convert Z into dBZ
   !
   dBZ(:) = undefined
   do k = bounds(iscan,jscan)%ke, bounds(iscan,jscan)%ks,-1
      if( scans_radar(iscan,jscan,k)%Zm /= undefined) then
          dBZ(k) =  10*log10(scans_radar(iscan,jscan,k)%Zm)
      endif
   enddo

   !
   ! determin echo top , bottom and thickness (1st-layer cloud )
   !
    hgt_echotop = 0.e0  ; k_echo_top = 0 !init   --------(1st cloud layer)
    LAYER_TOP_1st: do k = bounds(iscan,jscan)%ks, bounds(iscan,jscan)%ke-3
 
       if( dBZ(k)   >= -28.e0 .and. &  ! significant echo for sucessive 3 layers.
           dBZ(k+1) >= -28.e0 .and. &
           dBZ(k+2) >= -28.e0 ) then

           hgt_echotop = real(scans_radar(iscan,jscan,k)%hgt)*1.e-3   !echo top height [km]
           k_echo_top = k
           exit LAYER_TOP_1st
       endif
    enddo LAYER_TOP_1st

    if( k_echo_top == 0 ) cycle !no cloud

    hgt_echobotm = 0.e0 ;  k_echo_botm = bounds(iscan,jscan)%ke
    LAYER_BOT_1st: do k = k_echo_top, bounds(iscan,jscan)%ke
       if( dBZ(k) < -28.e0 ) then ! loosing significant echo
           hgt_echobotm  = real(scans_radar(iscan,jscan,k)%hgt)*1.e-3   !echo bottom height [km]
           k_echo_botm = k
           exit LAYER_BOT_1st
       endif
    enddo LAYER_BOT_1st

    echo_thick = hgt_echotop - hgt_echobotm  !physical echo thickness [km]

    !
    ! count statistics for 1st layer clouds
    !
    if( echo_thick > 0.5) then
        jj = int(hgt_echotop) ; ii = int(echo_thick)
        ! if bin index is within the maximum range, count statistics
        if( ii <= binmax .and. jj <= binmax ) then  
            isccp_cnt(ii,jj) = 1.e0 + isccp_cnt(ii,jj)
        endif
    endif
   if( k_echo_botm >= bounds(iscan,jscan)%ke-3 ) cycle  !if near surface cycle

   !
   ! determin echo top , bottom and thickness (2nd-layer cloud )
   !
    hgt_echotop = k_echo_botm  ; k_echo_top = k_echo_botm !init   --------(2nd cloud layer)
    LAYER_TOP_2nd: do k = k_echo_botm, bounds(iscan,jscan)%ke-3

       if( dBZ(k  ) >= -28.e0 .and. &  ! significant echo for sucessive 3 layers.
           dBZ(k+1) >= -28.e0 .and. &
           dBZ(k+2) >= -28.e0 ) then

           hgt_echotop = real(scans_radar(iscan,jscan,k)%hgt)*1.e-3   !echo top height [km]
           k_echo_top = k
           exit LAYER_TOP_2nd
       endif
    enddo LAYER_TOP_2nd

    if( k_echo_top == k_echo_botm  ) cycle !no cloud

    hgt_echobotm = 0.e0 ;  k_echo_botm = bounds(iscan,jscan)%ke
    LAYER_BOT_2nd: do k = k_echo_top, bounds(iscan,jscan)%ke
       if( dBZ(k) < -28.e0 ) then ! loosing significant echo
           hgt_echobotm  = real(scans_radar(iscan,jscan,k)%hgt)*1.e-3   !echo bottom height [km]
           k_echo_botm = k
           exit LAYER_BOT_2nd
       endif
    enddo LAYER_BOT_2nd

    echo_thick = hgt_echotop - hgt_echobotm  !physical echo thickness [km]

    !
    ! count statistics for 2nd layer clouds
    !
    if( echo_thick > 0.5) then
        jj = int(hgt_echotop) ; ii = int(echo_thick)
        ! if bin index is within the maximum range, count statistics
        if( ii <= binmax .and. jj <= binmax ) then
            isccp_cnt(ii,jj) = 1.e0 + isccp_cnt(ii,jj)
        endif
    endif
    if( k_echo_botm >= bounds(iscan,jscan)%ke-3 ) cycle  !if near surface cycle

 enddo
 enddo

!
! create joint isccp-like diagram
!
 isccp_pdf(:,:) = isccp_cnt(:,:) / sum( isccp_cnt(:,:) )  ! relative frequency

!
! write output 
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
               trim(sensor_name)//'.ISCCP.'//trim(overpass_tag)//trim(output_suffix)//'.bin'
 print("(a)"),'see output.ctl ->',trim(sdsu_io_file)

 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=(binmax+1)*(binmax+1)*4)
 write(io,rec=1) ((isccp_pdf(i,j),i=0,binmax),j=0,binmax)
 close(io)

!
! write grads control file
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//&
               trim(sensor_name)//'.ISCCP.'//trim(overpass_tag)//trim(output_suffix)//'.ctl'
 print("(a)"),'see output.ctl ->',trim(sdsu_io_file)

 open(io,file=trim(sdsu_io_file))
 sdsu_io_file =sdsu_inp_name(1:efile_len)//'.'//&
               trim(sensor_name)//'.ISCCP.'//trim(overpass_tag)//trim(output_suffix)//'.bin'
 write(io,"(a)") 'dset  '//trim(sdsu_io_file)
 write(io,*) 'title   CloudSat CPR ISCCP diagram '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    ',binmax+1,' linear  0     1'
 write(io,*) 'ydef    ',binmax+1,' linear  0     1'
 write(io,*) 'zdef      1 linear 0. 0.01'
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars 1                    '
 write(io,*) 'pd     0  0  probability       [-]'
 write(io,*) 'endvars                   '
 close(io)


 return
 end subroutine write_composite_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine zonal_freq_cloudsat_mask(proc)
 implicit none
 character(len=*),intent(in) :: proc

 logical,allocatable,save :: cloud_mask(:,:)
 real(sdsu_fps),allocatable,save :: cloud_count(:,:), &
                                    lat_cloud_count(:)         !deg
 integer,save :: jmax_mask, kmax
 integer :: j, jscan, iscan, k ,lev, k_normal  !looping
 real(sdsu_fps),allocatable,save :: dBZ(:)  !dB of Z
 integer,parameter :: io = 102


 select case(proc) !--------------------------------------------------------------------------

 case('init')
  jmax_mask  =INT( (maxlat - minlat) / res_latlon )
  kmax  = int( mxhgt_radar / range_radar )

  allocate( cloud_mask( jmax_mask, kmax), cloud_count( jmax_mask, kmax ), lat_cloud_count( jmax_mask ) )
  cloud_count = 0. !initialization

! assign latitude bin
  do j = 1, jmax_mask
     lat_cloud_count(j) = minlat + res_latlon * real(j) - 0.5 * res_latlon
  enddo
 if(allocated(dBZ)) deallocate(dBZ)
 allocate(dBZ( kmax ))

! print*, 'MSG zonal_freq_cloudsat_mask: kmax=',kmax, mxhgt_radar, range_radar

 case('proc') !--------------------------------------------------------------------------------

!
!
!
! print*,'initialize cloud mask'
 do j = 1, jmax_mask ;  do lev = 1, kmax
    cloud_mask(j,lev) = .false. ! initialization 
 enddo ; enddo

!
! now estimate joint historgram
!
 do jscan = js_scan , je_scan
 do iscan = 1, 1

   do j = 1, jmax_mask

      if( real(scans(iscan,jscan)%latlon_fov(1)/degrad) >  lat_cloud_count(j) - 0.5 * res_latlon .and. &
          real(scans(iscan,jscan)%latlon_fov(1)/degrad) <= lat_cloud_count(j) + 0.5 * res_latlon ) then
             do k = bounds(iscan,jscan)%ks, bounds(iscan,jscan)%ke
                lev = k - bounds(iscan,jscan)%ks + 1
                if( scans_radar(iscan,jscan,k)%Zm > 0. ) then
                   if( real(10*log10(scans_radar(iscan,jscan,k)%Zm))  >= -28.) then
                        ! print*,'fill it',jscan, k,lev
                        cloud_mask(j,lev) = .true.
                   endif
                endif
             enddo
      endif

   enddo


 enddo !iscan
 enddo !jscan

!
! accumulate cloud counting
!
 do j = 1, jmax_mask ;  do lev = 1, kmax

    if( cloud_mask(j,lev) .eqv. .true. ) then
         cloud_count(j,lev) = cloud_count(j,lev) + 1.
    endif

 enddo ; enddo


 case('fine') !----------------------------------------------------------------------------

 print*, 'MSG zonal_freq_cloudsat_mask: finish proces'
!
! write out clouf frequency zonal profile
!
 sdsu_io_file = trim(sdsu_dir_output)//'zonal_cloud_freq.grads'
 print*,'see grads output file ->',trim(sdsu_io_file)

 open(10,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=jmax_mask*kmax*4)
 print*,'see output ->',trim(sdsu_io_file)
 write(10,rec=1) (( cloud_count(j,lev),j=1,jmax_mask),lev=kmax,1,-1)

 close(10)

!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'zonal_cloud_freq.ctl'
 print*,'see output.ctl ->',trim(sdsu_io_file)
 open(10,file=trim(sdsu_io_file))

 sdsu_io_file = './zonal_cloud_freq.grads'
 write(10,*) 'dset  ',trim(sdsu_io_file)
 write(10,*) 'title   cloud frequency'
 write(10,*) 'undef   -999.             '
 write(10,*) 'xdef    1 linear 1. 1'
 write(10,*) 'ydef    ',jmax_mask,' linear ',minlat,' ',res_latlon
 write(10,*) 'zdef    ',kmax,' linear 0. 0.239'
 write(10,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(10,*) 'vars 1                    '
 write(10,*) 'cnt   ',kmax,' 0  cloud count [-]'
 write(10,*) 'endvars                   '
 close(10)

 end select

 end subroutine zonal_freq_cloudsat_mask

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump_radar_range( sensor_name, var3d , bounds, io, irec )
 implicit none
 character(len=*),intent(in) :: sensor_name
 real(sdsu_fps),intent(in)     :: var3d(:, :, :)
 type(radar_bounds),intent(in) :: bounds( :, : )
 integer,intent(in) :: io
 integer,intent(inout) :: irec
 integer :: i, j, k, kmax, k_inc, ks_inc
 real(sdsu_fps), allocatable :: temp3d(:,:,:)

 integer :: is,ie,js,je,ks,ke
 integer :: bnd(3)  !upper memory bound for 3D array

 bnd = UBOUND(var3d)

!
! define bounds index
!

 is = 1 ; ie=bnd(1) ; js=1 ; je=bnd(2) ; ks = 1 ; ke = bnd(3)


 if( trim(sensor_name) == 'DPR_Ka_HS' ) then

    kmax = ( bounds(is,js)%ke - bounds(is,js)%ks + 1 ) / 2
    if(allocated(temp3d) ) deallocate(temp3d)
    allocate( temp3d(is:ie,js:je,1:kmax ) )

    do j = js, je ; do i = is, ie
       kmax = ( bounds(i,j)%ke - bounds(i,j)%ks + 1 ) / 2
       ks_inc = bounds(i,j)%ks

       do k = 1, kmax

          k_inc = ks_inc + (k-1) * 2
          if    ( var3d( i,j,k_inc ) /= undefined .and. var3d( i,j,k_inc+1 ) /= undefined) then
              temp3d( i,j,k ) = 0.5*( var3d( i,j,k_inc ) + var3d( i,j,k_inc+1 ) )
          elseif( var3d( i,j,k_inc ) == undefined .and. var3d( i,j,k_inc+1 ) /= undefined) then
              temp3d( i,j,k ) =                            var3d( i,j,k_inc+1 )
          elseif( var3d( i,j,k_inc ) /= undefined .and. var3d( i,j,k_inc+1 ) == undefined) then
              temp3d( i,j,k ) =       var3d( i,j,k_inc )
          elseif( var3d( i,j,k_inc ) == undefined .and. var3d( i,j,k_inc+1 ) == undefined) then
              temp3d( i,j,k ) = undefined
          endif

       enddo
    enddo ; enddo
    write(io,rec=irec) (((temp3d(i,j,k), i=is,ie), j=js,je), k= kmax,1,-1)

 else !default case

    kmax = bounds(is,js)%ke - bounds(is,js)%ks + 1
    if(allocated(temp3d) ) deallocate(temp3d)
    allocate( temp3d(is:ie,js:je,1:kmax ) )
    temp3d = undefined

    do j = js, je ; do i = is, ie
       kmax = bounds(i,j)%ke - bounds(i,j)%ks + 1
       k_inc = bounds(i,j)%ks

       do k = 1, kmax
          temp3d( i,j,k ) = var3d(i,j,k_inc)
          k_inc = k_inc + 1
       enddo 

    enddo ; enddo
    write(io,rec=irec) (((temp3d(i,j,k), i=is,ie), j=js,je), k= kmax,1,-1)

 endif

 irec=irec+1

 return
 end subroutine dump_radar_range

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_geolocation_netcdf ( sensor_name, is, ie )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  This routine output satellite orbital data in NetCDF format. 
!
! History:
! 4/2012  Toshi Matsui@NASA GSFC    : initial
!
! References:
!----------------------------------------------------------------------
 character(*),intent(in) :: sensor_name  !sensor name
 integer,intent(in) :: is, ie

 character(len=50) :: sn  !sensor name
 integer :: i,j,k,m, js,je,ks,ke !simplified index for i j bounds
 integer :: kmax

 integer :: ncid
 integer :: i_dimid, j_dimid, k_dimid
 integer, parameter :: ndims3d = 3, ndims2d = 2, ndims1d=1
 integer :: dimids1d(ndims1d), dimids2d(ndims2d), dimids3d(ndims3d)
 integer :: id
 integer :: n ,ich, nf, np
 integer :: im  !middle of sample
 integer :: id_month, id_day, id_hour, id_min, id_sec, id_gps !var IDs (time)
 integer :: id_sat_alt, id_ecef_x, id_ecef_y, id_ecef_z, &
                        id_vel_x , id_vel_y , id_vel_z      ! var IDs (nagivation)
 integer :: id_fov_lat, id_fov_lon, id_sat_lat, id_sat_lon , id_inc ! var IDs (location)
 
!
! output file name
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(sensor_name)//&
  '.GEOLOCATION.'//trim(overpass_tag)//trim(output_suffix)//'.nc'
 print("(a)"),'MSG write_geolocation_netcdf: output ->',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), NF90_64BIT_OFFSET, ncid) )

!
! set up array bounds
!
 im    = (is+ie)/2 !
 if( is == ie ) im = is

 js = 1 ; je = jmax_scan  !full orbit

!
!  -------------------- define header file --------------------------
!

!
! find out start vertical index and end vertical index (20bin below ellipsoid, 130bin above) 
!
! if( trim(sensor_name) == 'DPR_Ka_HS' ) then
!   kmax = int( mxhgt_radar / range_radar )/2  ! DPR Ka High Sensivity (250m range bin)
! else
!   kmax = int( mxhgt_radar / range_radar )    ! bounds(imax_scan,je_scan)%ke - bounds(imax_scan,je_scan)%ks + 1 !toshii
! endif

!
! Define dimensions
!
 call check( nf90_def_dim(ncid, 'n_sample' , ie-is+1 , i_dimid) )
 call check( nf90_def_dim(ncid, 'n_scan'   , je-js+1 , j_dimid) )
! call check( nf90_def_dim(ncid, 'n_range'  , kmax    , k_dimid) )

! dimids3d = (/ i_dimid, j_dimid, k_dimid /)
 dimids2d = (/ i_dimid, j_dimid /)
 dimids1d = (/ j_dimid /)

!
! time header & ID
!
 call check( nf90_def_var(ncid, 'month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[month]') )
 id_month = id

 call check( nf90_def_var(ncid, 'day_of_month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Days of Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[day]') )
 id_day = id

 call check( nf90_def_var(ncid, 'hour', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (hour) )') )
 call check( nf90_put_att(ncid, id, 'units', '[hour]') )
 id_hour = id

 call check( nf90_def_var(ncid, 'minute', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (minute) )') )
 call check( nf90_put_att(ncid, id, 'units', '[min]') )
 id_min = id

 call check( nf90_def_var(ncid, 'second', NF90_DOUBLE, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (second) )') )
 call check( nf90_put_att(ncid, id, 'units', '[sec]') )
 id_sec = id

 call check( nf90_def_var(ncid, 'GPS_atomic_time', NF90_DOUBLE, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' GPS atomic time (second) )') )
 call check( nf90_put_att(ncid, id, 'units', '[sec]') )
 id_gps = id


!
! satllite navigation
!
 call check( nf90_def_var(ncid, 'sat_altitude', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite altitude above Earth reference )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_sat_alt = id

 call check( nf90_def_var(ncid, 'sat_ecef_x', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite ECEF position (x-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_ecef_x = id

 call check( nf90_def_var(ncid, 'sat_ecef_y', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite ECEF position (y-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_ecef_y = id

 call check( nf90_def_var(ncid, 'sat_ecef_z', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite ECEF position (z-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_ecef_z = id

 call check( nf90_def_var(ncid, 'sat_vel_x', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite velocity (x-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km/sec]') )
 id_vel_x = id

 call check( nf90_def_var(ncid, 'sat_vel_y', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite velocity (y-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km/sec]') )
 id_vel_y = id

 call check( nf90_def_var(ncid, 'sat_vel_z', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite velocity (z-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km/sec]') )
 id_vel_z = id


!
! geolocation headher & ID
!
 call check( nf90_def_var(ncid, 'sat_latitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite geolocation (lat) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_sat_lat = id

 call check( nf90_def_var(ncid, 'sat_longitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite geolocation (lon) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_sat_lon = id

 call check( nf90_def_var(ncid, 'fov_latitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' center efov geolocation (lat) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_fov_lat = id

 call check( nf90_def_var(ncid, 'fov_longitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' center efov geolocation (lon) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_fov_lon = id

 call check( nf90_def_var(ncid, 'incident_angle', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' sensor incident angle') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_inc = id

 call check( nf90_enddef(ncid) )  ! End define mode.

!
! ----------------------------- DUMP DATA ------------------------------------------
!


!
! dump time group
!
  call check( nf90_put_var(ncid, id_month, scans(im,js:je)%time%mm)  )
  call check( nf90_put_var(ncid, id_day  , scans(im,js:je)%time%dd)  )
  call check( nf90_put_var(ncid, id_hour , scans(im,js:je)%time%hh)  )
  call check( nf90_put_var(ncid, id_min  , scans(im,js:je)%time%mn)  )
  call check( nf90_put_var(ncid, id_sec  , scans(im,js:je)%time%ss)  )
  call check( nf90_put_var(ncid, id_gps  , scans(im,js:je)%GPSatomic))

!
! dump satellite navigation group
!
  call check( nf90_put_var(ncid, id_sat_alt, real(scans(im,js:je)%alt_sat)      ) )
  call check( nf90_put_var(ncid, id_ecef_x , real(scans(im,js:je)%vector_sat(1))) )
  call check( nf90_put_var(ncid, id_ecef_y , real(scans(im,js:je)%vector_sat(2))) )
  call check( nf90_put_var(ncid, id_ecef_z , real(scans(im,js:je)%vector_sat(3))) )
  call check( nf90_put_var(ncid, id_vel_x  , real(scans(im,js:je)%vel_sat(1))   ) )
  call check( nf90_put_var(ncid, id_vel_y  , real(scans(im,js:je)%vel_sat(2))   ) )
  call check( nf90_put_var(ncid, id_vel_z  , real(scans(im,js:je)%vel_sat(3))   ) )

!
! dump satellite FOV position and incident angle [deg]
!
  call check( nf90_put_var(ncid, id_sat_lat, real(scans(is:ie,js:je)%latlon_sat(1)/degrad)) )
  call check( nf90_put_var(ncid, id_sat_lon, real(scans(is:ie,js:je)%latlon_sat(2)/degrad)) )
  call check( nf90_put_var(ncid, id_fov_lat, real(scans(is:ie,js:je)%latlon_fov(1)/degrad)) )
  call check( nf90_put_var(ncid, id_fov_lon, real(scans(is:ie,js:je)%latlon_fov(2)/degrad)) )
  call check( nf90_put_var(ncid, id_inc    , real(scans(is:ie,js:je)%inci_ang     /degrad)) )


 call check( nf90_close(ncid) )  ! close nc file


 return
 end subroutine write_geolocation_netcdf

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_orbital_radar_netcdf ( sensor_name, is, ie )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  This routine output satellite orbital data in NetCDF format. 
!
! History:
! 4/2012  Toshi Matsui@NASA GSFC    : initial
!
! References:
!----------------------------------------------------------------------
 character(*),intent(in) :: sensor_name  !sensor name
 integer,intent(in) :: is, ie

 character(len=50) :: sn  !sensor name
 integer :: i,j,k,m, js,je,ks,ke !simplified index for i j bounds
 integer :: kmax, k_index

 integer :: ncid
 integer :: i_dimid, j_dimid, k_dimid
 integer, parameter :: ndims3d = 3, ndims2d = 2, ndims1d=1
 integer :: dimids1d(ndims1d), dimids2d(ndims2d), dimids3d(ndims3d)
 integer :: id
 integer :: n ,ich, nf, np
 integer :: im  !middle of sample
 integer :: id_month, id_day, id_hour, id_min, id_sec, id_gps !var IDs (time)
 integer :: id_sat_alt, id_ecef_x, id_ecef_y, id_ecef_z, &
                        id_vel_x , id_vel_y , id_vel_z      ! var IDs (nagivation)
 integer :: id_fov_lat, id_fov_lon, id_sat_lat, id_sat_lon , id_inc ! var IDs (location)
 integer :: id_gain, id_zt, id_zm, id_piasrt,id_piae ! var IDs ( radar params )
 integer :: id_hgt, id_press, id_tair, id_sh, id_w, & ! var IDs ( L2 params )
            id_iwc, id_lwc, id_rwc, id_lfrac, &
            id_rhoe, id_rain, id_mmt(0:6)

 real(sdsu_fps) , allocatable :: temp2d(:,:)
 logical,parameter  :: big_file = .false.  ! .true. write 64bit to write large NetCDF (>2Gb) 
 
!
! output file name
!
 sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(sensor_name)//&
  '.ORBITAL.'//trim(overpass_tag)//trim(output_suffix)//'.nc'
 print("(a)"),'MSG write_orbital_radar_netcdf: output ->',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), NF90_64BIT_OFFSET , ncid) )  !64bit output

!
! set up array bounds
!
 im    = (is+ie)/2 !
 if( is == ie ) im = is

 js = js_scan ; je = je_scan



!
!  -------------------- define header file --------------------------
!

!
! find out start vertical index and end vertical index (20bin below ellipsoid, 130bin above) 
!
 if( trim(sensor_name) == 'DPR_Ka_HS' ) then
   kmax = int( mxhgt_radar / range_radar )/2  ! DPR Ka High Sensivity (250m range bin)
 else
   kmax = int( mxhgt_radar / range_radar )    ! bounds(imax_scan,je_scan)%ke - bounds(imax_scan,je_scan)%ks + 1 !toshii
 endif

!
! Define dimensions
!
 call check( nf90_def_dim(ncid, 'n_sample' , ie-is+1 , i_dimid) )
 call check( nf90_def_dim(ncid, 'n_scan'   , je-js+1 , j_dimid) )
 call check( nf90_def_dim(ncid, 'n_range'  , kmax    , k_dimid) )

 dimids3d = (/ i_dimid, j_dimid, k_dimid /)
 dimids2d = (/ i_dimid, j_dimid /)
 dimids1d = (/ j_dimid /)

!
! time header & ID
!
 call check( nf90_def_var(ncid, 'month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[month]') )
 id_month = id

 call check( nf90_def_var(ncid, 'day_of_month', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Gregorian Calendar (Days of Month) )') )
 call check( nf90_put_att(ncid, id, 'units', '[day]') )
 id_day = id

 call check( nf90_def_var(ncid, 'hour', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (hour) )') )
 call check( nf90_put_att(ncid, id, 'units', '[hour]') )
 id_hour = id

 call check( nf90_def_var(ncid, 'minute', NF90_INT, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (minute) )') )
 call check( nf90_put_att(ncid, id, 'units', '[min]') )
 id_min = id

 call check( nf90_def_var(ncid, 'second', NF90_DOUBLE, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' UTC Time (second) )') )
 call check( nf90_put_att(ncid, id, 'units', '[sec]') )
 id_sec = id

 call check( nf90_def_var(ncid, 'GPS_atomic_time', NF90_DOUBLE, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' GPS atomic time (second) )') )
 call check( nf90_put_att(ncid, id, 'units', '[sec]') )
 id_gps = id


!
! satllite navigation
!
 call check( nf90_def_var(ncid, 'sat_altitude', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite altitude above Earth reference )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_sat_alt = id

 call check( nf90_def_var(ncid, 'sat_ecef_x', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite ECEF position (x-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_ecef_x = id

 call check( nf90_def_var(ncid, 'sat_ecef_y', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite ECEF position (y-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_ecef_y = id

 call check( nf90_def_var(ncid, 'sat_ecef_z', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite ECEF position (z-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km]') )
 id_ecef_z = id

 call check( nf90_def_var(ncid, 'sat_vel_x', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite velocity (x-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km/sec]') )
 id_vel_x = id

 call check( nf90_def_var(ncid, 'sat_vel_y', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite velocity (y-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km/sec]') )
 id_vel_y = id

 call check( nf90_def_var(ncid, 'sat_vel_z', NF90_REAL, dimids1d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite velocity (z-compoment) )') )
 call check( nf90_put_att(ncid, id, 'units', '[km/sec]') )
 id_vel_z = id


!
! geolocation headher & ID
!
 call check( nf90_def_var(ncid, 'sat_latitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite geolocation (lat) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_sat_lat = id

 call check( nf90_def_var(ncid, 'sat_longitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' satellite geolocation (lon) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_sat_lon = id

 call check( nf90_def_var(ncid, 'fov_latitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' center efov geolocation (lat) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_fov_lat = id

 call check( nf90_def_var(ncid, 'fov_longitude', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' center efov geolocation (lon) )') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_fov_lon = id

 call check( nf90_def_var(ncid, 'incident_angle', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' sensor incident angle') )
 call check( nf90_put_att(ncid, id, 'units', '[deg]') )
 id_inc = id



!
! L1B/L2 radar reflectivity 
!
 call check( nf90_def_var(ncid, 'gain', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' Normalized gain function (<1 is bad sampling)') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_gain = id


 call check( nf90_def_var(ncid, 'Zt', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' non-attenuated equivalent reflectivity factor') )
 call check( nf90_put_att(ncid, id, 'units', '[mm6/m3]') )
 id_zt = id

 call check( nf90_def_var(ncid, 'Zm', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' attenuated equivalent reflectivity factor') )
 call check( nf90_put_att(ncid, id, 'units', '[mm6/m3]') )
 id_zm = id

 call check( nf90_def_var(ncid, 'PIAsrt', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' path-integrated attenuation at SRT level') )
 call check( nf90_put_att(ncid, id, 'units', '[dB]') )
 id_piasrt = id

 call check( nf90_def_var(ncid, 'PIAe', NF90_REAL, dimids2d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' effective path-integrated attenuation above clutter level') )
 call check( nf90_put_att(ncid, id, 'units', '[dB]') )
 id_piae = id

!
! ancillary data or L2 output 
!
 call check( nf90_def_var(ncid, 'hgt', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' range height above ellipsoid ') )
 call check( nf90_put_att(ncid, id, 'units', '[m]') )
 id_hgt = id

 call check( nf90_def_var(ncid, 'press', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' atmospheric pressure ') )
 call check( nf90_put_att(ncid, id, 'units', '[hPa]') )
 id_press = id

 call check( nf90_def_var(ncid, 'tair', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' air temperature ') )
 call check( nf90_put_att(ncid, id, 'units', '[degK]') )
 id_tair = id

 call check( nf90_def_var(ncid, 'sh', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' specific humidity ') )
 call check( nf90_put_att(ncid, id, 'units', '[g/g]') )
 id_sh = id

 call check( nf90_def_var(ncid, 'w', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' vertical velocity') )
 call check( nf90_put_att(ncid, id, 'units', '[m/s]') )
 id_w = id

 call check( nf90_def_var(ncid, 'iwc', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' cloud ice water content ') )
 call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
 id_iwc = id

 call check( nf90_def_var(ncid, 'lwc', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' cloud liquid water content ') )
 call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
 id_lwc = id

 call check( nf90_def_var(ncid, 'rwc', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' rain water content ') )
 call check( nf90_put_att(ncid, id, 'units', '[g/m3]') )
 id_rwc = id

 call check( nf90_def_var(ncid, 'lfrac', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' liquid fraction of rain ') )
 call check( nf90_put_att(ncid, id, 'units', '[-]') )
 id_lfrac = id

 call check( nf90_def_var(ncid, 'rhoe', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' bulk effective density of rain ') )
 call check( nf90_put_att(ncid, id, 'units', '[g/cm3]') )
 id_rhoe = id

 call check( nf90_def_var(ncid, 'rain', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' rainfall rate ') )
 call check( nf90_put_att(ncid, id, 'units', '[mm/hr]') )
 id_rain = id

 call check( nf90_def_var(ncid, 'mmt0', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' 0th moment ') )
 call check( nf90_put_att(ncid, id, 'units', '[m0/m3]') )
 id_mmt(0) = id

 call check( nf90_def_var(ncid, 'mmt1', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' 1st moment ') )
 call check( nf90_put_att(ncid, id, 'units', '[m1/m3]') )
 id_mmt(1) = id

 call check( nf90_def_var(ncid, 'mmt2', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' 2nd moment ') )
 call check( nf90_put_att(ncid, id, 'units', '[m2/m3]') )
 id_mmt(2) = id

 call check( nf90_def_var(ncid, 'mmt3', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' 3rd moment ') )
 call check( nf90_put_att(ncid, id, 'units', '[m3/m3]') )
 id_mmt(3) = id

 call check( nf90_def_var(ncid, 'mmt4', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' 4th moment ') )
 call check( nf90_put_att(ncid, id, 'units', '[m4/m3]') )
 id_mmt(4) = id

 call check( nf90_def_var(ncid, 'mmt5', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' 5th moment ') )
 call check( nf90_put_att(ncid, id, 'units', '[m5/m3]') )
 id_mmt(5) = id

 call check( nf90_def_var(ncid, 'mmt6', NF90_REAL, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'MemoryOrder','(n_range,n_sample,n_scan)') )
 call check( nf90_put_att(ncid, id, 'description',' 6th moment ') )
 call check( nf90_put_att(ncid, id, 'units', '[m6/m3]') )
 id_mmt(6) = id


 call check( nf90_enddef(ncid) )  ! End define mode.


!
! ----------------------------- DUMP DATA ------------------------------------------
!


!
! dump time group
!
  call check( nf90_put_var(ncid, id_month, scans(im,js:je)%time%mm)  )
  call check( nf90_put_var(ncid, id_day  , scans(im,js:je)%time%dd)  )
  call check( nf90_put_var(ncid, id_hour , scans(im,js:je)%time%hh)  )
  call check( nf90_put_var(ncid, id_min  , scans(im,js:je)%time%mn)  )
  call check( nf90_put_var(ncid, id_sec  , scans(im,js:je)%time%ss)  )
  call check( nf90_put_var(ncid, id_gps  , scans(im,js:je)%GPSatomic))

!
! dump satellite navigation group
!
  call check( nf90_put_var(ncid, id_sat_alt, real(scans(im,js:je)%alt_sat)      ) )
  call check( nf90_put_var(ncid, id_ecef_x , real(scans(im,js:je)%vector_sat(1))) )
  call check( nf90_put_var(ncid, id_ecef_y , real(scans(im,js:je)%vector_sat(2))) )
  call check( nf90_put_var(ncid, id_ecef_z , real(scans(im,js:je)%vector_sat(3))) )
  call check( nf90_put_var(ncid, id_vel_x  , real(scans(im,js:je)%vel_sat(1))   ) )
  call check( nf90_put_var(ncid, id_vel_y  , real(scans(im,js:je)%vel_sat(2))   ) )
  call check( nf90_put_var(ncid, id_vel_z  , real(scans(im,js:je)%vel_sat(3))   ) )

!
! dump satellite FOV position and incident angle [deg]
!
  call check( nf90_put_var(ncid, id_sat_lat, real(scans(is:ie,js:je)%latlon_sat(1)/degrad)) )
  call check( nf90_put_var(ncid, id_sat_lon, real(scans(is:ie,js:je)%latlon_sat(2)/degrad)) )
  call check( nf90_put_var(ncid, id_fov_lat, real(scans(is:ie,js:je)%latlon_fov(1)/degrad)) )
  call check( nf90_put_var(ncid, id_fov_lon, real(scans(is:ie,js:je)%latlon_fov(2)/degrad)) )
  call check( nf90_put_var(ncid, id_inc    , real(scans(is:ie,js:je)%inci_ang     /degrad)) )


  ks = 1 ; ke = irange_max
  sn = sensor_name
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%gain     , bounds(is:ie,js:je), ncid, id_gain )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%Zt       , bounds(is:ie,js:je), ncid, id_zt )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%Zm       , bounds(is:ie,js:je), ncid, id_zm )
  call dump_radar_range_nc( sn, real(scans_radar(is:ie,js:je,ks:ke)%hgt), bounds(is:ie,js:je), ncid, id_hgt )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%iwc      , bounds(is:ie,js:je), ncid, id_iwc )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%lwc      , bounds(is:ie,js:je), ncid, id_lwc )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%rwc      , bounds(is:ie,js:je), ncid, id_rwc )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%rain_rate, bounds(is:ie,js:je), ncid, id_rain )

#ifdef FULL_L2_RADAR
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%press    , bounds(is:ie,js:je), ncid, id_press )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%t_air    , bounds(is:ie,js:je), ncid, id_tair )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%sh       , bounds(is:ie,js:je), ncid, id_sh )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%w        , bounds(is:ie,js:je), ncid, id_w  )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%lfrac    , bounds(is:ie,js:je), ncid, id_lfrac )
  call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%bulk_den , bounds(is:ie,js:je), ncid, id_rhoe )
  do m = 0, 6
     call dump_radar_range_nc( sn, scans_radar(is:ie,js:je,ks:ke)%mmt(m), bounds(is:ie,js:je), ncid, id_mmt(m) )
  enddo
#endif

!
! PIA here defined as dB(Zt) - dB(Zm)
!
 if(allocated(temp2d)) deallocate(temp2d)
 allocate( temp2d(is:ie,js:je) )

  do j = js,je ; do i = is,ie
     k_index = bounds(i,j)%k_srt
     if( scans_radar(i,j, k_index )%Zm == undefined ) then
        temp2d(i,j) = undefined
     elseif( scans_radar(i,j, k_index )%Zm == 0. ) then !if no return, use tau
        temp2d(i,j) = 20.* scans_radar(i,j,k_index)%tau * log10(2.71828183)
     else  !if echo is significant, use dBZt-dBZm
        temp2d(i,j) =  10.*log10( scans_radar(i,j, k_index )%Zt/scans_radar(i,j, k_index )%Zm)
     endif
  enddo ; enddo
  call check( nf90_put_var(ncid, id_piasrt , temp2d(is:ie,js:je) ) ) 

  do j = js,je ; do i = is,ie
     k_index = bounds(i,j)%k_piae
     if( scans_radar(i,j, k_index )%Zm == undefined ) then
        temp2d(i,j) = undefined
     elseif( scans_radar(i,j, k_index )%Zm == 0. ) then !if no return, use tau
        temp2d(i,j) = 20.* scans_radar(i,j,k_index)%tau * log10(2.71828183)
     else  !if echo is significant, use dBZt-dBZm
        temp2d(i,j) = 10.*log10( scans_radar(i,j, k_index )%Zt/scans_radar(i,j, k_index )%Zm)
     endif
  enddo ; enddo
  call check( nf90_put_var(ncid, id_piae , temp2d(is:ie,js:je) ) )

 call check( nf90_close(ncid) )  ! close nc file

 return
 end subroutine write_orbital_radar_netcdf

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump_radar_range_nc( sensor_name, var3d , bounds, ncid, id )
 implicit none
 character(len=*),intent(in) :: sensor_name
 real(sdsu_fps),intent(in)     :: var3d(:, :, :)
 type(radar_bounds),intent(in) :: bounds( :, : )
 integer,intent(in) :: ncid
 integer,intent(in) :: id

 integer :: i, j, k, kmax, k_inc, ks_inc
 real(sdsu_fps), allocatable :: temp3d(:,:,:)

 integer :: is,ie,js,je,ks,ke
 integer :: bnd(3)  !upper memory bound for 3D array

 bnd = UBOUND(var3d)

!
! define bounds index
!

 is = 1 ; ie=bnd(1) ; js=1 ; je=bnd(2) ; ks = 1 ; ke = bnd(3)


 if( trim(sensor_name) == 'DPR_Ka_HS' ) then

    kmax = ( bounds(is,js)%ke - bounds(is,js)%ks + 1 ) / 2
    if(allocated(temp3d) ) deallocate(temp3d)
    allocate( temp3d(is:ie,js:je,1:kmax ) )

    do j = js, je ; do i = is, ie
       kmax = ( bounds(i,j)%ke - bounds(i,j)%ks + 1 ) / 2
       ks_inc = bounds(i,j)%ks

       do k = 1, kmax

          k_inc = ks_inc + (k-1) * 2
          if    ( var3d( i,j,k_inc ) /= undefined .and. var3d( i,j,k_inc+1 ) /= undefined) then
              temp3d( i,j,k ) = 0.5*( var3d( i,j,k_inc ) + var3d( i,j,k_inc+1 ) )
          elseif( var3d( i,j,k_inc ) == undefined .and. var3d( i,j,k_inc+1 ) /= undefined) then
              temp3d( i,j,k ) =                            var3d( i,j,k_inc+1 )
          elseif( var3d( i,j,k_inc ) /= undefined .and. var3d( i,j,k_inc+1 ) == undefined) then
              temp3d( i,j,k ) =       var3d( i,j,k_inc )
          elseif( var3d( i,j,k_inc ) == undefined .and. var3d( i,j,k_inc+1 ) == undefined) then
              temp3d( i,j,k ) = undefined
          endif

       enddo
    enddo ; enddo
!    write(io,rec=irec) (((temp3d(i,j,k), i=is,ie), j=js,je), k= kmax,1,-1)

    call check( nf90_put_var( ncid, id, temp3d(is:ie,js:je,1:kmax) ) )


 else !default case

    kmax = bounds(is,js)%ke - bounds(is,js)%ks + 1
    if(allocated(temp3d) ) deallocate(temp3d)
    allocate( temp3d(is:ie,js:je,1:kmax ) )
    temp3d = undefined

    do j = js, je ; do i = is, ie
       kmax = bounds(i,j)%ke - bounds(i,j)%ks + 1
       k_inc = bounds(i,j)%ks

       do k = 1, kmax
          temp3d( i,j,k ) = var3d(i,j,k_inc)
          k_inc = k_inc + 1
       enddo 

    enddo ; enddo
    call check( nf90_put_var( ncid, id, temp3d(is:ie,js:je,1:kmax) ) )

 endif

 return
 end subroutine dump_radar_range_nc

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine deallocate_all_scan
 implicit none
 integer :: ierr

 deallocate( scans    , &
             scans_L2 , &
             surface_CRM,&
             stat=ierr  )
 if (ierr /= 0) call stop_sdsu(myrank, 'MSG deallocate_all_scan: deallocation error -> Terminate program. ')

 if(  write_grid_data ) then
     deallocate( cnt_scan, &
                 cnt_sat , &
                 time_scan,&
                 time_sat ,&
                 stat=ierr  )
 endif

 return
 end subroutine deallocate_all_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
