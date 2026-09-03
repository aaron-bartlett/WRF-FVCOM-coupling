#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 program main_SDSU
 use module_psd
 use module_tmatrix_database
 use module_simulator
 use module_gas_profile
 use module_sarta
 use module_surface
 use module_opt_micro
 use module_opt_visir
 use module_opt_broad
 use module_GV
 use module_mpi
 use module_radar
 use module_scan
 use module_icescat_yang_visir
 use module_gmp_4ice
 use module_p3_param
 use module_diag
 use module_driver_polarris
 use module_openssp_ext
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
!
! Comments: 
!  This program is the main driver of the G-SDSU, which controls all procedure according to the 
!  user-defined run-time options in Configure_SDSU.F. 
!
! History:
! 08/2016  Toshi Matsui@NASA GSFC : Add POLARRIS (polarimetric radar simulator)
! 01/2016  Toshi Matsui@NASA GSFC : Add zonal cloud frequency 
! 02/2012  Toshi Matsui@NASA GSFC : Add scan simulator for A-Train satellite. 
! 09/2011  Toshi Matsui@NASA GSFC : Add scan simulator for the GPM mission. 
! 05/2011  Toshi Matsui@NASA GSFC : Add preliminarl version of GEOS5. 
! 04/2011  Toshi Matsui@NASA GSFC : Add new HUCM_SBM43 option. 
! 06/2008  Toshi Matsui@NASA GSFC : Add MPI options via c-preprocessor. 
! 04/2008  Toshi Matsui@NASA GSFC : Add 5 different WRF microphysics (DSD) options.
! 01/2008  Toshi Matsui@NASA GSFC : Add Mie diagnostic routine. 
! 01/2008  Toshi Matsui@NASA GSFC : Add Goddard SW/LW scheme (Broadband & ISCCP-like simulators)
! 12/2007  Toshi Matsui@NASA GSFC ; Add MMF input options.
! 09/2007  Toshi Matsui@NASA GSFC ; Add Lidar options. 
! 08/2007  Toshi Matsui@NASA GSFC ; Add GCE SBM options.
! 03/2007  Toshi Matsui@NASA GSFC ; Initial
!
! References:
! 
! Matsui, T., J. Santanello, J. J. Shi, W.-K. Tao, D. Wu, C. Peters-Lidard, E. Kemp, M. Chin, 
!   D. Starr, M. Sekiguchi, and F. Aires, (2014): Introducing multisensor satellite radiance-based 
!   evaluation for regional Earth System modeling, Journal of Geophysical Research, 119, 8450-8475, 
!   doi:10.1002/2013JD021424. 
!
!-----------------------------------------------------------------------------------------------------
 integer :: n ,nproc   ! looping indice
 real(sdsu_fpd) :: ts, te  ! starting & ending wall-clock time
 character(len=100) :: LIS_file_name, GEOS_file_name
 logical :: exit_proc  ! if true exit entire rest of the processes. 
 integer :: min_since  ! [min]
 logical :: fexist, failed_opt_micro
!
! read sdsu.config file & input file name list & configu file
!

 call read_config
 call read_input_filename
 call init_scan_sensor  !optional for orbit simulator
 call init_channel_character 

!
! special routines to prepare scattering LUT for POLARRIS-F. 
! (Needed to build here to differentiate MPI tasks other than file/domain decompositions)
!
 if(polarris .and. generate_lut_polarris) call prep_makelut_polarris


!
! initialize major/minor gaseous profiles
!
 if( visir .or. lidar ) call init_gas_profile

!
! re-define maximum grid for GCM Lat Lon cordinate
!
 if( trim(sim_case) == 'GEOS5_forecast' )  then
   call reconfigure_grid
 endif

!
! if scan simulator is active, read satelite overpass info 
!
 if(micro .and. scan_micro) call read_overpass('micro')
 if(radar .and. scan_radar) call read_overpass('radar')
 if(visir .and. scan_visir) call read_overpass('visir')
 if(lidar .and. scan_lidar) call read_overpass('lidar')
 if(broad .and. scan_broad) call read_overpass('broad')

!
! Define MPI decomposition parametrs 
!
  exit_proc = .false. !initialize as false
#if MPI == 1

 call mpi_sdsu_init  ! initialize MPI
 call mpi_sdsu_file(sdsu_nmax_file,exit_proc)  ! soring input files for each processor
  myi_start = 1  ; myi_end   = mxgridx
  myj_start = 1  ; myj_end   = mxgridy
  myk_start = 1  ; myk_end   = mxlyr

#elif MPI == 2

 call mpi_sdsu_init  ! initialize MPI
 call mpi_sdsu_tile(1,mxgridx,1,mxgridy,1,mxlyr) ! sorting domain into tiles for each processor
  myn_start = 1  ; myn_end   = sdsu_nmax_file

#else

  ! Default values of MPI parameters for single CPU run 
  masterproc = .true.
  myrank = 0
  myn_start = 1  ; myn_end   = sdsu_nmax_file
  myi_start = 1  ; myi_end   = mxgridx
  myj_start = 1  ; myj_end   = mxgridy
  myk_start = 1  ; myk_end   = mxlyr

#endif


#if MPI > 0
  call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
#endif

 if( trim(sim_case) == 'GEOS5_forecast' )  then
   if(masterproc) then
     print*,''
     print*,'Model grid numbers are re-configured for GEOS5_forecast case.'
     print*,'mxgridx   =',mxgridx
     print*,'mxgridy   =',mxgridy
     print*,'mxlyr     =',mxlyr 
     print*,''
   endif
 endif


!
! prepare Tmatrix database
!
 if( (micro .or. radar) .and. ss_opt_micro == 4 ) &
   call tmatrix_prep( masterproc, trim(sdsu_dir_data) )

!
! prepare icescat database for visbir-lidar simulators
!
 if( (visir .or. lidar .or. broad) .and. ss_opt_visir == 1 ) then
    sdsu_io_file = trim(sdsu_dir_data)//'ICESCAT_YANG/'  !data directory
    call prep_icescat_yang_visir( trim(sdsu_io_file), ss_opt_visir )
 endif

!
! prepare OpenSSP database for microwave-radar simulator
! 
 if( (micro .or. radar) .and. ss_opt_micro == 3 ) then
!    sdsu_io_file = trim(sdsu_dir_data)//'OpenSSP/'  !data directory
!    call readTablesLiang2(masterproc, trim(sdsu_io_file), 5,8)  !5 mus (2,2,2,2,2) 8 frequencies (10,19,22,36,89,166,183.3+3,183.3+7)
    sdsu_io_file = trim(sdsu_dir_data)//'OpenSSP_EXT/'  !data directory
    call read_openssp_ext(masterproc, trim(sdsu_io_file) )  !5 mus (2,2,2,2,2) 8 frequencies (10,19,22,36,89,166,183.3+3,183.3+7)


 endif


!
! initialize polarris constant
!
! if( polarris ) call read_polarris_axiscant

!
! initialize GV simulator
!
 if ( GV ) call GV_init 
 if ( GV .and. aircraft_on .and. zonal_sampling_on) call aircraft_zonal_sampling_psd('init')
 if ( GV .and. parsivel_on .and. zonal_sampling_on) call parsivel_zonal_sampling_psd('init')

!
! initialize global statistical routine for MMF_GPROF evaluation
!
 if(trim(sim_case) == 'MMF_GPROF') then
   call scan_angle_bin
   call tpw_tskin_precip_sample('init')
 endif


  if(masterproc) print*,'####################################################### '
  if(masterproc) print*,'##################### START G-SDSU #################### '
  if(masterproc) print*,'####################################################### '
  if(masterproc) print*,''
  if(masterproc) print*,'MSG SDSU: G-SDSU will processes ', sdsu_nmax_file ,'files'
  if(masterproc) print*,''

#if MPI > 0
  call mpi_wait_for_master  !this process hold slave processors to wait for master proc. 
#endif


!
! Set up PSD-density parameters for different microphysics scheme
!
    mic_select0: select case(trim(type_microphysics))
     case('GEN') ; call bulk_DSD      ! DSD for exponential DSDs
     case('SBM')   ; call read_SBM_bin  ! SBM input from X.Li
     case('RAMS')           ; 
     case default                  
    end select mic_select0
    call re_rams_gamma('init')  ! RAMS PSD is called for any simulation case.

    if( trim(cloud_microphysics) == 'P3' ) call makelut_p3_re_ice(masterproc)  !LUT generator for P3 ice


    if(account_aerosol) call read_gocart_psd  ! read GOCART PSD functions and RH grouth parameters

#if MPI > 0
  call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
#endif

!
! Create Single-Scattering LUT for faster simulation
! 
  lut_if: if(masterproc) then

    if(trim(sim_case) == 'MLM' .or. trim(sim_case) == 'SOUNDING' ) then
    else
       if( broad                       ) call makelut_broad        ! LUT for broadband simulator
       if( broad .and. account_aerosol ) call makelut_broad_aerosol! LUT for broadband simulator
       if( lut_micro .and. micro ) call makelut_micro    ('micro') ! LUT for microwave radiometer simultion
       if( lut_micro .and. radar ) call makelut_micro    ('radar') ! LUT for radar simultion
       if( lut_visir .and. visir ) call makelut_visir    ('visir') ! LUT for visir imager simulation
       if( lut_visir .and. lidar ) call makelut_visir    ('lidar') ! LUT for lidar simulation
       if( lut_visir .and. visir .and. account_aerosol)  call makelut_visir_aerosol('visir') ! LUT for aerosol
       if( lut_visir .and. lidar .and. account_aerosol)  call makelut_visir_aerosol('lidar') ! LUT for aerosol
       if(                 visir ) call makelut_visir_grnd         ! surface LUT for visir imager
    endif !MLM


  endif lut_if


  if( account_aerosol ) call makelut_ccn_icn

#if MPI > 0
  call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
#endif


!
! dynamically allocate memory for the CRM parameters
!
  call allocate_all

!
! initialize time series statistics
! (following diagnostic subouroutine DOES NOT work for MPI=1 option (file decomposition))
#if MPI == 1
#else
 call diag_domain_mask  !create diag mask surface%diag_mask

 if( diag_crm_column_mean ) call dmain_mean_profile('init',mxlyr)  !basic domain mean statistics
 if( activate_zonal_freq_cloudsat_mask ) call zonal_freq_cloudsat_mask('init')
 if( diag_crm_spatial_mean ) call spatial_mean_2d('init')
 if( diag_crm_spatial_mean ) call spatial_mean_3d('init')
 if( diag_crm_diurnal_mean ) call diurnal_mean_2d('init')
 if( diag_crm_diurnal_mean ) call diurnal_mean_3d('init')
#endif 

!
! initialize SARTA (AIRS) simulator
!
 if ( sarta ) call SARTA_DRIVER('init')

!
! initialize PSD sampler from SBM
!
 if( radar .and. radar_hid .and. trim(cloud_microphysics(1:4)) == 'HUCM') then
    call sample_psd_sbm('init')
 endif

!--------------------------------- BEGIN FILE LOOP -------------------------------------------------

  file_loop: do n = myn_start, myn_end

    if( exit_proc ) cycle  !this ensure to cycle proces, when processor is not assigned for task in MPI=1 options

    sdsu_inp_name = trim(sdsu_inp_list(n))

!
! compute GMT time from sdsu_inp_name 
!
   if ( trim(sim_case) == 'LIS' ) then
     call LIS_filename_convert( sdsu_inp_list(n) , LIS_file_name )  
     call get_gmt( LIS_file_name )
   elseif ( trim(sim_case) == 'GEOS5' ) then
     call GEOS_filename_convert( sdsu_inp_list(n) , GEOS_file_name )
     call get_gmt( GEOS_file_name )
   elseif ( trim(sim_case) == 'GEOS5_forecast' ) then
     call get_gmt_geos( sdsu_inp_list(n) )
   else
     call get_gmt( sdsu_inp_name )
   endif

!
! get grads time character for output purpose
!
   if(write_grads_ctl)   call get_grads_time


!
! MMF_GPROF skipping option if output is available (we may get rid of this in future.....)
!
   if(trim(sim_case) == 'MMF_GPROF') then !Special output for  GPROF  
       sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.'//trim(micro_sensor)//&
                     '.GPROF'//trim(output_suffix)//'.nc'
       inquire(file=sdsu_io_file,exist=fexist)
       if (fexist) then !skip
           print("(a)"),'MSG output exisit ->',trim(sdsu_io_file)
           print("(a)"),'skip'
           cycle
       endif
   endif


!
! Read input CRM parameters 
!
   call cpu_time(ts)
   crm_select1: select case(trim(sim_case))
!     case('GCE2D')    ; call rd_CRM_GCE2D      ! 2D GCE input from X.Li
     case('GCE')      ; call rd_CRM_GCE_unified! Goddard Cumulus Ensemble (GCE) model
     case('GCESBM')   ; call rd_CRM_GCESBM     ! Goddard Cumulus Ensemble (GCE) SBM model
     case('WRF')      ; call rd_CRM_WRF        ! Weather Research and Forecast (WRF) model
     case('MMF')      ; call rd_CRM_MMF        ! Goddard Multi-scale Modeling Framework (MMF) 
     case('MMF_GPROF'); call rd_CRM_MMF_GPROF  ! Goddard Multi-scale Modeling Framework (MMF) 
     case('GEOS5')    ; call rd_CRM_GEOS5      ! GEOS5 (preliminarl) 
     case('GEOS5_forecast') ; call rd_CRM_GEOS5_forecast_latlon  ! GEOS5 forecast in classic lat lon format
     case('SOUNDING') ; call rd_CRM_SOUNDING   ! Sounding data (SOUNDING)
     case('MLM')      ; call rd_CRM_MLM        ! Mixed-Layer Model (MLM) 
     case('LIS')      ; call rd_CRM_LIS        ! Land Information Systems + Reanalysis (LIS) 
     case default     ; call stop_sdsu(myrank,'MSG main_SDSU: There is no such sim_case')
   end select crm_select1
   call cpu_time(te)
   if(masterproc) print*,'total CPU time of reading input file is',te-ts  ; if(masterproc) print*,''


 if(bad_sdsu_io) cycle !If input data does not exixt skip the rest of routines. ...

!
! Set constrains environmental parameter for an idealized scene. 
!
  call idealized_scene 

!
! read aerosol from GLOBAL GOCART simulation
!
 if(account_aerosol .and. nudge_gocart_on) call nudge_gocart_offline

!
! compute various parameters for simulator or diagnostics purpose. 
!

  call re_all ! compute drop effective radius

  call get_others  !compute miscellaneous environmental/particle parameters

  if(.not. write_cftext) then
    call get_solar_parameter ! compute variou solar parameters 
  endif

!
! special PSD output from HUCM_SBM
!
 if(write_psdmmt) then
    call cpu_time(ts)
    call write_out_sbm_moment
    call cpu_time(te)
    if(masterproc) print*,'total CPU time of writing SBM moment file is',te-ts  ; if(masterproc) print*,''
 endif


!
! #######################   POLARRIS   #############################
!
  if ( polarris ) then
     call cpu_time(ts)
     call opt_polarris
     call cpu_time(te)
     if(masterproc) print*,'total CPU time of POLARRIS is',te-ts  ; if(masterproc) print*,''
  endif

!
! #######################  GV simulator  #############################
!
  if ( GV ) then
     call cpu_time(ts)
     if(aircraft_on) call GV_aircraft
     if(aircraft_on .and. zonal_sampling_on) call aircraft_zonal_sampling_psd('add')
     if(parsivel_on) call GV_parsivel
     if(parsivel_on .and. zonal_sampling_on) call parsivel_zonal_sampling_psd('add')
     call cpu_time(te) 
     if(masterproc) print*,'total CPU time of GV Simulator is',te-ts  ; if(masterproc) print*,''
  endif 

!
! ######################  Radar Simulator ###############################
!
  if ( radar ) then
     call cpu_time(ts)
     call surface_radar

     if( scan_radar ) then
         call scan_simulator('radar',satellite_radar,scan_type_radar,&
                             overpass_radar(n) )  !simulate radar scan vector
         call opt_micro('radar',surface_CRM(:,:)%within_swath,&
                        real(stat%max_inc_angle)/const_degrad ) !compute radar optical properties
         call drive_radar_simulator_scan
     else
        call opt_micro('radar')  !compute microwave optical properties
        call simulator_radar()   !simulate radar reflectivity
!        if(trim(sim_case) /= 'MMF_GPROF') then
!           call write_out_radar    ! Write out simulated radar reflectivity. 
!           if(radar_output_L2) call write_out_L2_radar ! Write L2 geophysical parameter
!        endif
     endif
     call cpu_time(te) 
     if(masterproc) print*,'total CPU time of Radar Simulator is',te-ts  ; if(masterproc) print*,''
  endif





!
! #################   Passive Microwave Simulator   #######################
!
  if ( micro ) then
     call cpu_time(ts)

     if( scan_micro ) then  

         call scan_simulator('micro',satellite_micro,scan_type_micro,&
                              overpass_micro(n) )                !simulate satellite scan vector
         call surface_micro( surface_CRM(:,:)%umu_local )        !compute surface emissivity
         call opt_micro('micro', surface_CRM(:,:)%within_swath,&
                         real(stat%max_inc_angle)/const_degrad)  !compute microwave optical properties
         call drive_micro_simulator_scan                         !drive slant-path microwave simulator

     else 

        if(trim(sim_case) == 'MMF_GPROF') then !Special output for  GPROF  

           call prep_scan_array_prps  !calculate incidence angle and gain function

           call surface_micro_prps( )  ! compute surface emissivity

           call opt_micro( sensor='micro',failed_opt_micro=failed_opt_micro) ! compute optical properties of microwave range
           if( failed_opt_micro ) cycle

           if(slant_path_micro) then
             call simulator_slant_micro   ! compute microwave Tb on slant path.
           else
             call simulator_micro_prps   ! compute microwave Tb on CRM grid.
           endif
           call write_out_prps   ! write out GPROF database
                   
        else  

           call surface_micro( )     !compute surface emissivity
           call opt_micro('micro')  !compute optical properties of microwave range
           if(slant_path_micro) then
             call simulator_slant_micro   ! compute microwave Tb on slant path.
           else
             call simulator_micro   ! compute microwave Tb on CRM grid.
           endif
#if MPI == 2
           if(masterproc) &  !for domain decomposition case, write by masterproc only.
#endif
           call write_out_micro ! - Write out Simulated brightness temperature. 

        endif

     endif

     call cpu_time(te) 
     if(masterproc) print*,'total CPU time of Microwave Simulator is',te-ts ; if(masterproc) print*,''
  endif


!
! #####################   Broadband SW/LW simulator   ############################
!
  if (broad) then
      call cpu_time(ts)
      call surface_broad  !compute broadband surface albedo
      call opt_broad      !compute optical properties and drive broadband flux simulator
#if MPI == 2
      if(masterproc) &  !for domain decomposition case, write by masterproc only.
#endif
      call write_out_broad  !write energy budget and heating ratio
      call cpu_time(te) 
      if(masterproc) print*,'total CPU time of Broadband Simulator is',te-ts ; if(masterproc) print*,''
  endif

!
! #########################    Lidar Simulator   ##################################
!
  if ( lidar ) then
      call cpu_time(ts)

      if( scan_lidar ) then !orbit + scan simulator

          call scan_simulator('lidar',satellite_lidar,scan_type_lidar,&
                               overpass_lidar(n) ) !simulate lidar scan vector
          call opt_visir('lidar',surface_CRM(:,:)%umu_local) !compute lidar optics and backscatter

          call drive_lidar_simulator_scan  !write out orbital lidar output
      else

        call opt_visir('lidar') 
#if MPI == 2
        if(masterproc) &  !for domain decomposition case, write by masterproc only.
#endif
        call write_out_lidar  ! - write out simulated lidar backscattering coefficient

      endif

      call cpu_time(te) 
      if(masterproc) print*,'total CPU time of Lidar Simulator is',te-ts ; if(masterproc) print*,''
  endif

!
! ########################   Passive vis-IR simulator  ###############################
!
  if ( visir ) then
      call cpu_time(ts)
      call surface_visir !compute surface spectrum albedo/emissivity

      if( scan_visir ) then  !orbit+scan simulator 

          call scan_simulator('visir',satellite_visir,scan_type_visir,&
                               overpass_visir(n) )  !compute satelltie sensor scan vector
          call opt_visir('visir',surface_CRM(:,:)%umu_local) !compute visir optical properties and radiance.

          call drive_visir_simulator_scan ! write out orbital output

      else ! normal column simulator

        call opt_visir('visir')  
#if MPI == 2
        if(masterproc) then  !for domain decomposition case, write by masterproc only.
#endif
          call write_out_visir ! - Write out Simulated visible/IR radiance in native CRM domain
!          if( visir_sensor /= 'AIRS' )&
!          call write_out_visir_latlon ! write out in defined lat-lon domain.
#if MPI == 2
         endif 
#endif

      endif

      call cpu_time(te) 
      if(masterproc) print*,'total CPU time of VisibleIR simulator is',te-ts ; if(masterproc) print*,''
  endif

!
! ########################   SARTA (ARIS)  simulator  ###############################
!
 if ( sarta ) then 
     call cpu_time(ts)
     call surface_sarta        ! compute surface spectrum emissivity

     call SARTA_DRIVER('calc') ! compute AIRS radiance/Tb

!     call write_out_sarta_nc  ! in netcdf format (does not work for larger output size)
     call write_out_sarta

     call cpu_time(te)
     if(masterproc) print*,'total CPU time of SARTA simulator is',te-ts ; if(masterproc) print*,''
 endif

!
! ########### CALIPSO-CloudSat Mask ###########
!
! this only works for special case( MMF with CloudSat CPR and CALIPSO CALIOP simuators)
  call write_out_cc_mask 

!
! ##################### Write diagnostic ouput ###########################
!
  if(write_CRM3D) then 
    call cpu_time(ts)
    call write_out_CRM3D !write oute default 3D output
    call cpu_time(te)
    if(masterproc) print*,'total CPU time of writing CRM3D file is',te-ts ; if(masterproc) print*,''
  endif

  if(write_CRM2D)  then
    call cpu_time(ts)
    call write_out_CRM2D !write out default 2D output
    call cpu_time(te)
    if(masterproc) print*,'total CPU time of writing CRM2D file is',te-ts ; if(masterproc) print*,''
  endif


!
! following diagnostic subouroutine DOES NOT work for MPI=1 option (file decomposition)
!
#if MPI == 1
#else
  if( diag_crm_column_mean ) call dmain_mean_profile('add',mxlyr, q_gen, n_gen, atmos%w, atmos_stag%hgt, atmos%dhgt, &
                                                     atmos%sctot, atmos%latent_heat, atmos%cloud_nuc, &
                                                     atmos%rho_dair, surface%diag_mask, surface%iclass )  
  if( diag_crm_spatial_mean ) call spatial_mean_2d('add',surface)
  if( diag_crm_spatial_mean ) call spatial_mean_3d('add',atmos, surface)
  if( diag_crm_diurnal_mean ) call diurnal_mean_2d('add',surface)
  if( diag_crm_diurnal_mean ) call diurnal_mean_3d('add',atmos, atmos_stag%press, surface)

  if( diag_crm_spatial_mean_daily ) then

      if( sdsu_hh == 0. .or. n == 1 ) call spatial_mean_daily_2d('init')
      call spatial_mean_daily_2d('add',surface)
      if( sdsu_hh == 23.) call spatial_mean_daily_2d('finish')

      if( sdsu_hh == 0. .or. n == 1 ) call spatial_mean_daily_3d('init')
      call spatial_mean_daily_3d('add',atmos,surface)
      if( sdsu_hh == 23.) call spatial_mean_daily_3d('finish')

  endif

  if( diag_crm_diurnal_daily ) then 
      if( sdsu_hh == 0. .or. n == 1 ) call diurnal_daily_2d('init')
      call diurnal_daily_2d('add',surface)
      if( sdsu_hh == 23.) call diurnal_daily_2d('finish')
  endif

! modify for sub-hourly sample later (toshii)
  if( diag_crm_spatial_mean_hourly ) then

      call spatial_mean_hourly_2d('init')
      call spatial_mean_hourly_2d('add',surface)
      call spatial_mean_hourly_2d('finish')

      call spatial_mean_hourly_3d('init')
      call spatial_mean_hourly_3d('add',atmos,surface)
      call spatial_mean_hourly_3d('finish')

  endif

#endif

 enddo file_loop

!------------------------------------- END OF FILE LOOP -------------------------------------

#if MPI > 0
  call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
#endif


!
! finalize GV zonal sampling
!
  if(GV .and. aircraft_on .and. zonal_sampling_on) call aircraft_zonal_sampling_psd('finish')  !Finalize domain sampling
  if(GV .and. parsivel_on .and. zonal_sampling_on) call parsivel_zonal_sampling_psd('finish')  !Finalize domain sampling


 if( radar .and. radar_hid .and. trim(cloud_microphysics(1:4)) == 'HUCM') then
    call sample_psd_sbm('finish')  !finialize PDF sampling and dump ouput
 endif



!
! finalize time-series statistics (does not work for MPI=1 option)
!
#if MPI == 1
#else
  if( activate_zonal_freq_cloudsat_mask )  call zonal_freq_cloudsat_mask('fine') !derive zonal cloud mask 
  if( diag_crm_column_mean  ) call dmain_mean_profile('finish',mxlyr)  !basic domain mean statistics
  if( diag_crm_spatial_mean ) call spatial_mean_2d('finish')
  if( diag_crm_spatial_mean ) call spatial_mean_3d('finish')
  if( diag_crm_diurnal_mean ) call diurnal_mean_2d('finish')
  if( diag_crm_diurnal_mean ) call diurnal_mean_3d('finish')
#endif


!
! finalize GPROF statistics
!
 if(trim(sim_case) == 'MMF_GPROF') then
   call tb_precip_prps('finish')
   call tpw_tskin_precip_sample('finish')
 endif

  if(masterproc) print*,'####################################################### '
  if(masterproc) print*,'#################### FINISH G-SDSU #################### '
  if(masterproc) print*,'####################################################### '

#if MPI > 0
!
! terminate mpi process
!
 call mpi_sdsu_finale  !terminate mpi process
#endif

 end program main_SDSU 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
