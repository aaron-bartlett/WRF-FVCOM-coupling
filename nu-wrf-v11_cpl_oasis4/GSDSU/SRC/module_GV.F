#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_GV
 use netcdf
 use module_simulator
 use module_gmp_4ice
#if MPI > 0 
 use module_mpi
#endif
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
!  This module is designed to create GV-observable parameters from WRF SBM.  
!  Global Precipitation Measurement (GPM) Ground Validation (GV) is a series of 
!  field campaigns that support rainfall-algorithm developement. 
!
! History:
! 09/2017  Toshi Matsui@NASA GSFC ; Major update of SBM 43bin + including bulk microphysics. 
! 03/2010  Toshi Matsui@NASA GSFC ; Initial.
! 11/2011  Toshi Matsui@NASA GSFC ; update Parsivel simulator (simplified only account for surface). 
! References:
!-----------------------------------------------------------------------------------------------------
 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 private   ! ALL variables and subourtines are non accessible in module_GV.

 public :: GV_init, GV_aircraft, GV_parsivel, &
           aircraft_zonal_sampling_psd, &
           parsivel_zonal_sampling_psd    ! public (accessible) subroutines

!
! ---------------------- Aircraft Parameters --------------------------------------------------------------
! Size-bin (based on maximum diameter) of particle size distribution (PSD) analysis 
! used by Andy Heymsfield@NCAR from aircraft 2D-P & 2D-C measurements. 
!

! integer,parameter :: maxbin_2D = 33  !maximum bin number of aircraft 2D analysis
! real(sdsu_fps), parameter :: psd_endbin_2D(1:maxbin_2D+1) = & ! PSD end bin in aircraft 2D analysis [micron] 
!    (/   15.e0,   45.e0,   75.e0,  105.e0,  135.e0,  195.e0,  255.e0,  315.e0,  375.e0,  435.e0,&
!        495.e0,  550.e0,  650.e0,  750.e0,  850.e0,  950.e0, 1400.e0, 1800.e0, 2200.e0, 2600.e0,&
!       3000.e0, 3400.e0, 3800.e0, 4600.e0, 5400.e0, 6200.e0, 7000.e0, 9000.e0,11000.e0,13000.e0,&
!      15000.e0,19000.e0,23000.e0,27000.e0 /)


 integer,parameter :: maxbin_2D = 39  !maximum bin number of aircraft 2D analysis
 real(sdsu_fps), parameter :: psd_endbin_2D(1:maxbin_2D+1) = & ! PSD end bin in aircraft 2D analysis [micron] 
    (/   15.e0,   45.e0,   75.e0,  105.e0,  135.e0,  195.e0,  255.e0,  315.e0,  375.e0,  435.e0,&
        495.e0,  550.e0,  650.e0,  750.e0,  850.e0,  950.e0, 1400.e0, 1800.e0, 2200.e0, 2600.e0,&
       3000.e0, 3400.e0, 3800.e0, 4600.e0, 5400.e0, 6200.e0, 7000.e0, 9000.e0,11000.e0,13000.e0,&
      15000.e0,19000.e0,23000.e0,27000.e0,31000.e0,35000.e0,39000.e0,43000.e0,47000.e0,51000.e0 /)



!
! GV aircraft parameters dimension(mxgridx,mxgridy,mxlyr)
!
 type aircraft_parameter
    real(sdsu_fps) :: bulk_re      ! bulk effective radius compabile to aircraft 2D measurements [micron]
    real(sdsu_fps) :: bulk_wc      ! bulk water content compabile to aircraft CVI measurements [g/m3]
    real(sdsu_fps) :: bulk_vol     ! bulk particle volume  [cm3/m3]
    real(sdsu_fps) :: bulk_rho     ! bulk effective density compabile to 2D-CVI analysis [g/cm3]
    real(sdsu_fps) :: Ntot         ! total number conc [#/m3]
    real(sdsu_fps) :: bulk_psd(maxbin_2D) ! bulk PSD density [#/m4]
 end type aircraft_parameter
 type ( aircraft_parameter ), allocatable, dimension(:,:,:) :: aircraft !gv aircraft parameter

!
! ---------------------- Parsivel Parameters --------------------------------------------------------------
! Size-bin (based on maximum diameter) of particle size distribution (PSD) analysis 
! used by Parsivel measurements/analysis. 
!
! Note from Ali Tokay: 
! 1. Parsivel tends to underestimate particle less than 0.8mm (800micron). 
! 2. Parsivel tends to overestimate particle with very large drops.
! 3. For liquid particle, Parsivel account up to 23rd bin (~8mm diameter, which is maximum liquid drop size.)
!

 real(sdsu_fps),parameter :: dt_parsivel = 10.e0 ! Parsivel sampling interval [sec] 

 integer,parameter :: maxbin_parsivel = 32  !maximum bin number of parsivel analysis
 real(sdsu_fps), parameter :: psd_endbin_parsivel(maxbin_parsivel+1) = & ! PSD end bin in parsivel analysis [micron] 
    (/    0.0e0,  128.5e0,  256.5e0,  385.5e0,  514.5e0,  643.5e0,  771.5e0,  900.5e0, 1029.5e0, 1158.5e0,&
       1287.5e0, 1545.5e0, 1802.5e0, 2060.0e0, 2317.5e0, 2574.5e0, 3089.5e0, 3604.5e0, 4120.5e0, 4634.5e0,&
       5150.0e0, 6180.0e0, 7210.0e0, 8240.0e0, 9270.0e0,10300.0e0,12360.0e0,14420.0e0,16480.0e0,18580.0e0,&
      20600.0e0,23690.0e0,26780.0e0 /)

!
! GV parsivel parameters dimension(mxgridx,mxgridy)
!
 type parsivel_parameter
    real(sdsu_fps) :: parsivel_rain    ! rain rate measured by Parsivel (volumed-based ) [mm/hr]
    real(sdsu_fps) :: true_rain        ! true rain rate fall [mm/hr]
    real(sdsu_fps) :: bulk_rho         ! bulk effective density compabile to Parsivel-Geonor analysis [g/cm3]
    real(sdsu_fps) :: bulk_re          ! bulk effective radius compabile to Parsivel measurements [micron]
    real(sdsu_fps) :: Ntot             ! total number concentration [#/m3]
    real(sdsu_fps) :: bulk_psd(maxbin_parsivel) ! bulk PSD compabile to Parsivel measurements [#/m4]
 end type parsivel_parameter
 type ( parsivel_parameter ), allocatable, dimension(:,:) :: parsivel !gv parsivel parameter

!
! Non-spherical particle shape parameters for Aircraft 2D probe (Dmax = Wideth enclosing box)
!
 integer,parameter :: nmax_ratio = 30  ! maximum number of Dmax-Dsphere ratio bin
 integer,parameter :: nmax_dsph  = 10  ! maximum number of Dsphere-size bin

 type ( particle_sbm ),dimension(1:nmax_ratio,1:nmax_dsph) :: & 
     pd_sbm_2D   ,& ! Probability density of Dmax-Dsphere ratio for Aircraft 2D probe [0~1]
     pd_sbm_parsivel! Probability density of WHD-Dsphere ratio for PARSIVEL [0~1]

 type ( particle_sbm ),dimension(1:nmax_ratio) :: & 
     ratio_sbm    ! Ratio between Dmax and Dsphere [-] (X-axis)  

 type ( particle_sbm ),dimension(1:nmax_dsph) :: &  
     dsph_sbm     ! Dsphere [mm] (Y-axis)


!
! psd profile parameters
!
 real(sdsu_fps),allocatable,save :: psd_profile_cnt(:)
 real(sdsu_fps),allocatable,save :: psd_profile_sum(:,:)
 real(sdsu_fps),allocatable,save :: psd_profile_mean(:,:)

!
! psd parameters
!
 real(sdsu_fps),save :: psd_cnt
 real(sdsu_fps),allocatable,save :: psd_sum(:)
 real(sdsu_fps),allocatable,save :: psd_mean(:)


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine GV_init
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Read Dmax_Dsph ratio and Dmax_WHD tables.
!  Essentially, HUCM_SBM (and most of other scheme) assume sphere shape of particle, 
!  while a) an Aircraft 2D probe measures edge-to-edge length (width of enclosing box)of randomlly oriented
!  solid particle), and b) PARSIVEL measures widest horizontal dimension (WHD).  
!  So, in order to take account for irregular snow aggregate shape, 
!  we use probability density tables of Dmax-Dsph ratio a) measured from in-situ 2DVD (snow aggregation
!  and b) simulated from SnowFake model (ice crystal habit) for Aircraft 2D.  
!  Also, we use WHD_Dsph ratio measured from in-situ 2DVD. 
!  
!
! History:
! 06/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer,parameter :: io = 600
 integer :: i,n
 integer :: i_mode
 real(sdsu_fps) :: dummy
 real(sdsu_fps) :: residual

!
! open/read Dmax-Dsphere ratio PD tables.  -----------------------------------------------------------
!

 open(io,file=trim(sdsu_dir_data)//'Dmax_Dsph_Snow.asc',status='old')

 do n = nmax_ratio, 1 , -1
    read(io, *) ratio_sbm(n)%snow, pd_sbm_2D(n,1:nmax_dsph)%snow  !Y-axis, and probability
 enddo
 read(io, *) dummy, dsph_sbm(1:nmax_dsph)%snow  !X-axis
 close(io) 

!
! Correction. Sum of PD for a given Dsph index must be unity.
!
 do n = 1,nmax_dsph
    
    do i = 1, nmax_ratio 
       if( pd_sbm_2D(i,n)%snow == MAXVAL( pd_sbm_2D(:,n)%snow ) ) then
          i_mode = i ; exit
       endif
    enddo

    residual = 1.e0 - SUM(pd_sbm_2D(:,n)%snow)  

    pd_sbm_2D(i_mode,n)%snow = pd_sbm_2D(i_mode,n)%snow + residual  !put residual value on mode bin

 enddo

!
! Inaddition, PD of first two Dsph bins measured from 2DVD are not reliable.
! So, let's replace them by using 3rd bin.  
!
 pd_sbm_2D(:,1)%snow = pd_sbm_2D(:,3)%snow
 pd_sbm_2D(:,2)%snow = pd_sbm_2D(:,3)%snow


!
! open/read WHD-Dsphere ratio PD tables.  ---------------------------------------------------------
!
 open(io,file=trim(sdsu_dir_data)//'WHD_Dsph_Snow.asc',status='old')

 do n = nmax_ratio, 1 , -1
    read(io, *) ratio_sbm(n)%snow, pd_sbm_parsivel(n,1:nmax_dsph)%snow  !Y-axis, and probability
 enddo
 read(io, *) dummy, dsph_sbm(1:nmax_dsph)%snow  !X-axis
 close(io)

!
! Correction. Sum of PD for a given Dsph index must be unity.
!
 do n = 1,nmax_dsph

    do i = 1, nmax_ratio
       if( pd_sbm_parsivel(i,n)%snow == MAXVAL( pd_sbm_parsivel(:,n)%snow ) ) then
          i_mode = i ; exit
       endif
    enddo

    residual = 1.e0 - SUM(pd_sbm_parsivel(:,n)%snow)

    pd_sbm_parsivel(i_mode,n)%snow = pd_sbm_parsivel(i_mode,n)%snow + residual  !put residual value on mode bin
 enddo

!
! Inaddition, PD of first two Dsph bins measured from 2DVD are not reliable.
! So, let's replace them by using 3rd bin.  
!
 pd_sbm_parsivel(:,1)%snow = pd_sbm_parsivel(:,3)%snow
 pd_sbm_parsivel(:,2)%snow = pd_sbm_parsivel(:,3)%snow


 return
 end subroutine GV_init
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine GV_aircraft
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine deal with domain loop to compute GV parameters.
!
! subroutines
!
! History:
! 06/2018  Toshi Matsui@NASA GSFC ; Add bulk microphysics 
! 03/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps) :: bulk_psd(maxbin_2D)     ! bulk particle size distributions [#/m3]
 integer :: i,j,k, n, nb
 integer :: ierr
 integer :: is,ie,js,je,ks,ke
 real(sdsu_fps) :: psd_dbin_2D(maxbin_2D)  ! [m]

!
! Estimate width of size bin
!
 do nb = 1, maxbin_2D
    psd_dbin_2D(nb) = (psd_endbin_2D(nb+1) - psd_endbin_2D(nb) ) * 1.e-6  ![m]
 enddo


!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end


 if(masterproc) print*,'MSG GV_aircraft: start GV simulator. '


!
! parameter allocation
!
 if( .not. allocated(aircraft) ) then
    allocate( aircraft(is:ie,js:je,ks:ke), &
    stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG GV_aircraft: aircraft allocation error -> Terminate program.')
 endif
 aircraft%bulk_re = 0.e0 ; aircraft%bulk_wc = 0.e0 ; aircraft%bulk_vol = 0.e0 
 aircraft%bulk_rho = 0.e0 

!
! Estimate bulk parameters for each grid compatible to GV aircraft measurements. 
! Size-bin (based on maximum diameter) of particle size distribution (PSD) analysis 
! used by Andy Heymsfield@NCAR from aircraft 2D-P & 2D-C measurements.
!

 DOM_LOOP: do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end

   bulk_psd = 0.e0

   !
   ! Derive bulk water content from all species compatible to aircraft CVI measurements. 
   !
   select case(trim(type_microphysics))
   case('GEN')
       call bulk_water_content_gen ( q_gen(i,j,k) , aircraft(i,j,k)%bulk_wc )
   case('RAMS')
       call bulk_water_content_rams( q_rams(i,j,k), aircraft(i,j,k)%bulk_wc )
   case('SBM')
       call bulk_water_content_sbm ( q_sbm(i,j,k) , aircraft(i,j,k)%bulk_wc )
   end select   

   if( aircraft(i,j,k)%bulk_wc > 0.e0) then

      !
      ! Derive bulk particle size distribution compatible to aircraft 2D analysis.
      !
      select case(trim(type_microphysics))
      case('GEN')
         call bulk_psd_2D_gen( atmos(i,j,k)%t_air, q_gen(i,j,k), re_gen(i,j,k), mu_gen, rho_gen, &
                               maxbin_2D, psd_endbin_2D, bulk_psd )
      case('RAMS')
       stop  !add later toshii
      case('SBM')
         call bulk_psd_2D_sbm( n_sbm(i,j,k,1:nbin),frime_sbm_snow(i,j,k, 1:nbin), &
                               fmelt_sbm_snow(i,j,k, 1:nbin), &
                               maxbin_2D, psd_endbin_2D, bulk_psd )
      end select
     
      !
      ! total number conc [#/m3]
      !
      aircraft(i,j,k)%Ntot = SUM( bulk_psd(:) )      

      !
      ! Derive bulk effective radius and bulk particle volume compatible to aircraft 2D analysis.
      !
      call bulk_re_vol( maxbin_2D, psd_endbin_2D, bulk_psd, aircraft(i,j,k)%bulk_re, aircraft(i,j,k)%bulk_vol ) 

      !
      ! Derive bulk density compatible to aircraft 2D-CVI analysis.
      !
      if( aircraft(i,j,k)%bulk_vol == 0. ) then
        aircraft(i,j,k)%bulk_wc  = undefined
        aircraft(i,j,k)%bulk_vol = undefined
        aircraft(i,j,k)%bulk_re  = undefined
        aircraft(i,j,k)%bulk_rho = undefined
        aircraft(i,j,k)%Ntot     = undefined
        aircraft(i,j,k)%bulk_psd(:) = undefined

      else

        aircraft(i,j,k)%bulk_rho = aircraft(i,j,k)%bulk_wc / aircraft(i,j,k)%bulk_vol  ! [g/cm3]

        !
        ! PSD density [#/m4]
        !
        aircraft(i,j,k)%bulk_psd(:) = bulk_psd(:)/psd_dbin_2D(:) 
      endif


   else

      aircraft(i,j,k)%bulk_wc  = undefined
      aircraft(i,j,k)%bulk_vol = undefined
      aircraft(i,j,k)%bulk_re  = undefined
      aircraft(i,j,k)%bulk_rho = undefined
      aircraft(i,j,k)%Ntot     = undefined
      aircraft(i,j,k)%bulk_psd(:) = undefined

   endif

 enddo ; enddo ; enddo DOM_LOOP

!
! Write GV aircraft parameters 
!
 call write_GV_aircraft



 return
 end subroutine GV_aircraft

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine aircraft_zonal_sampling_psd(proc)
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

 integer :: i,j,k, kk, p, n, r             ! loop indice
 real(sdsu_fps) :: Ntot         ! total number concentration [#/m3]
 real(sdsu_fps) :: pd(maxbin_2D)  ! PSD [#/m3]


 integer,parameter :: io = 601
 integer,parameter :: method = 1  ! 1 (normalized PSD)  or 2 (PSD)

 real(sdsu_fps) :: psd_dbin_2D(maxbin_2D)    ! 2D bin width [m]
 real(sdsu_fps) :: psd_midbin_2D(maxbin_2D)  ! 2D mid bin size [micron]


 integer,parameter :: play = 20, plev = 21
 real(sdsu_fps),parameter :: pressure_lay(play) =    (/1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,500.,&
                                                     450.,400.,350.,300.,250.,200.,150.,100.,50./)
 real(sdsu_fps),parameter :: pressure_lev(plev) = (/1025., 975.,925.,875.,825.,775.,725.,675.,625.,575.,525.,475.,&
                                                        425.,375.,325.,275.,225.,175.,125.,75.,25./)

 real(sdsu_fps),dimension(1:play,1:maxbin_2D),save :: matrix, mean_pd  !sum of PDF 



 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

!
! allocate
!
 if( .not. allocated(psd_profile_cnt) ) then

    allocate( psd_profile_cnt(mxlyr), &
              psd_profile_sum(mxlyr,maxbin_2D), &
              psd_profile_mean(mxlyr,maxbin_2D) )

    allocate( psd_sum (maxbin_2D), &
              psd_mean(maxbin_2D) )

 endif

!
!initialize zero
!
 psd_profile_cnt = 0.
 psd_profile_sum = 0.
 psd_profile_mean = 0.

 psd_cnt = 0.
 psd_sum = 0.
 psd_mean = 0.

!
! get middle bin values
!
 do n = 1, maxbin_2D
    psd_midbin_2D(n) = 0.5*( psd_endbin_2D(n) + psd_endbin_2D(n+1) )
 enddo

 case('add') !-------------------------------------------------------------------------------------


!
! estimate normalized PSD 
!
 DOM_LOOP: do j = myj_start, myj_end ; do i = myi_start, myi_end

    !
    ! sampling good values within defined zone
    !
    if( surface(i,j)%lat >= minlat_gv .and. surface(i,j)%lat <= maxlat_gv .and.&
        surface(i,j)%lon >= minlon_gv .and. surface(i,j)%lon <= maxlon_gv ) then

        k_loop: do k = 1, mxlyr 
           if( aircraft(i,j,k)%bulk_wc > 0.e0 .and. surface(i,j)%rain_rate > 15.) then  !toshi convective only

             kk_loop: do kk = 1, play
                if( atmos(i,j,k)%press <= pressure_lev(kk) .and. &
                    atmos(i,j,k)%press >  pressure_lev(kk+1)  ) then
                    psd_profile_cnt(kk) = psd_profile_cnt(kk) + 1.0
                    psd_profile_sum(kk,1:maxbin_2D) = psd_profile_sum(kk,1:maxbin_2D) + aircraft(i,j,k)%bulk_psd(1:maxbin_2D)
                endif
             enddo kk_loop

             if( atmos(i,j,k)%hgt >= minalt_gv .and. atmos(i,j,k)%hgt <= maxalt_gv ) then
                psd_cnt = psd_cnt + 1.0
                psd_sum(1:maxbin_2D) = psd_sum(1:maxbin_2D) + aircraft(i,j,k)%bulk_psd(1:maxbin_2D)
             endif
           endif
        enddo k_loop

    endif

 enddo ; enddo DOM_LOOP


 case('finish') !---------------------------------------------------------------------------------------


#if MPI > 0 
  call mpi_sdsu_reduce('SUM', psd_profile_cnt )
  call mpi_sdsu_reduce('SUM', psd_profile_sum ) 
  call mpi_sdsu_reduce('SUM', psd_cnt )
  call mpi_sdsu_reduce('SUM', psd_sum )

  MPI_MASTERPROC: if(masterproc) then
#endif

!
! compute mean value
!
 do kk = 1, play
    if( psd_profile_cnt(kk) > 0. ) then
        psd_profile_mean(kk,1:maxbin_2D) = psd_profile_sum(kk,1:maxbin_2D) / psd_profile_cnt(kk)
    else
        psd_profile_mean(kk,1:maxbin_2D) = 0.
    endif
 enddo

 if( psd_cnt > 0. ) then
     psd_mean(1:maxbin_2D) = psd_sum(1:maxbin_2D) / psd_cnt
 else
     psd_mean(1:maxbin_2D) = 0.
 endif



!
! write profile output
!

 if(aircraft_ice) then
    sdsu_io_file =trim(sdsu_dir_output)//'AIRCRAFT_PSD_PROFILE_ICE'//trim(output_suffix)//'.txt'
 else
    sdsu_io_file =trim(sdsu_dir_output)//'AIRCRAFT_PSD_PROFILE_LIQ'//trim(output_suffix)//'.txt'
 endif


 print*,'MSG aircraft_zonal_sampling_psd: dump PSD profile composite. See --> ',trim(sdsu_io_file)
 open( io,file=trim(sdsu_io_file) )

 write(io,101) 0.0, psd_midbin_2D(1:maxbin_2D)
 do kk = 1, play
    write(io,100) pressure_lay(kk), psd_profile_mean(kk,1:maxbin_2D)
 enddo
 close(io)


!
! write mean output
!
 if(aircraft_ice) then
    sdsu_io_file =trim(sdsu_dir_output)//'AIRCRAFT_PSD_ICE'//trim(output_suffix)//'.txt'
 else
    sdsu_io_file =trim(sdsu_dir_output)//'AIRCRAFT_PSD_LIQ'//trim(output_suffix)//'.txt'
 endif

 print*,'MSG aircraft_zonal_sampling_psd: dump PSD mean composite. See --> ',trim(sdsu_io_file)
 open( io,file=trim(sdsu_io_file) )

 write(io,103) psd_midbin_2D(1:maxbin_2D)
 write(io,102) psd_mean(1:maxbin_2D)
 close(io)

 100 format(F10.3,39E10.3)
 101 format(40F10.3)

 102 format(39E10.3)
 103 format(39F10.3)

#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 end select 

 return
 end subroutine aircraft_zonal_sampling_psd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine aircraft_zonal_sampling_psd_re(proc)
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

 integer :: i,j,k,p, n, r             ! loop indice
 integer,parameter :: maxbin_re=10    ! 
 real(sdsu_fps) :: Ntot         ! total number concentration [#/m3]
 real(sdsu_fps) :: pd(maxbin_2D)  ! PSD [#/m3]
 real(sdsu_fps),dimension(1:maxbin_re,1:maxbin_2D),save :: matrix, mean_pd  !sum of PDF 
 real(sdsu_fps),save :: cnt(1:maxbin_re)  !sampling count for each re bin
 integer,parameter :: io = 601
 integer,parameter :: method = 1  ! 1 (normalized PSD)  or 2 (PSD)
 real(sdsu_fps) :: psd_dbin_2D(maxbin_2D)    ! 2D bin width [m]
 real(sdsu_fps) :: psd_midbin_2D(maxbin_2D)  ! 2D mid bin size [micron]

 real(sdsu_fps) :: re_bin(1:maxbin_re) ! re bin (center) [micron]
 real(sdsu_fps),parameter :: re_bin_end(0:maxbin_re) = & ! re bin (boundary) [micron]
                (/ 100.e0, 400., 500., 600., 650., 700., 900., 1200., &
                   2000., 4000., 20000. /)  


 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

!
! set up re bin size
!
 do r = 1, maxbin_re 
    re_bin(r) = 0.5e0 * ( re_bin_end(r-1) + re_bin_end(r) )
 enddo

 matrix = 0.e0
 cnt = 0.e0

!
! Estimate width of size bin
!
 do p = 1, maxbin_2D
    psd_dbin_2D(p) = (psd_endbin_2D(p+1) - psd_endbin_2D(p) ) * 1.e-6  !bin width [m]
 enddo

 do p = 1, maxbin_2D
    psd_midbin_2D(p) = 0.5e0*(psd_endbin_2D(p+1) + psd_endbin_2D(p) ) ! mid bin size [micron]
 enddo

 case('add') !-------------------------------------------------------------------------------------

!
! estimate normalized PSD 
!
 DOM_LOOP: do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end

    !
    ! sampling good values within defined zone
    !
    if( surface(i,j)%lat >= minlat_gv .and. surface(i,j)%lat <= maxlat_gv .and.&
        surface(i,j)%lon >= minlon_gv .and. surface(i,j)%lon <= maxlon_gv .and.&
        atmos(i,j,k)%hgt >= minalt_gv .and. atmos(i,j,k)%hgt <= maxalt_gv .and.&
        aircraft(i,j,k)%bulk_wc > 0.e0         ) then !selected zone

       Ntot = 0. 
       do p = 1, maxbin_2D
          Ntot = Ntot +   aircraft(i,j,k)%bulk_psd(p) * psd_dbin_2D(p)   ![#/m3]  
       enddo

       pd = 0.e0
       do p = 1, maxbin_2D
          pd(p) = aircraft(i,j,k)%bulk_psd(p) * psd_dbin_2D(p) / Ntot  ! fraction [-]  
       enddo

       re_loop: do r = 1, maxbin_re
          if( aircraft(i,j,k)%bulk_re >= re_bin_end(r-1) .and. aircraft(i,j,k)%bulk_re < re_bin_end(r) ) then
          
              if(method == 1)  matrix(r,1:maxbin_2D) = matrix(r,1:maxbin_2D) + pd(1:maxbin_2D)  !normalized PSD [-] 
              if(method == 2)  matrix(r,1:maxbin_2D) = matrix(r,1:maxbin_2D) + aircraft(i,j,k)%bulk_psd(1:maxbin_2D)  !PSD [1/m4]
              cnt(r)    = cnt(r) + 1.e0  !sampling count
              exit re_loop

          endif
       enddo re_loop

    endif

 enddo ; enddo ; enddo DOM_LOOP


 case('finish') !---------------------------------------------------------------------------------------


#if MPI > 0 
  call mpi_sdsu_reduce('SUM', matrix ) 
  call mpi_sdsu_reduce('SUM', cnt    )

  MPI_MASTERPROC: if(masterproc) then
#endif

!
! compute mean value
!
 do r = 1, maxbin_re 
    mean_pd(r,1:maxbin_2D) = matrix(r,1:maxbin_2D) / cnt(r)
 enddo

!
! write output
!

 sdsu_io_file =trim(sdsu_dir_output)//'AIRCRAFT_PSD_COMPOSITE'//trim(output_suffix)//'.txt'

 print*,'MSG aircraft_zonal_sampling_psd_re: dump PSD composite. See --> ',trim(sdsu_io_file)
 open( io,file=trim(sdsu_io_file) )

 write(io,101) 0.0, psd_midbin_2D(1:maxbin_2D)
 do r = 1, maxbin_re 
    if(method == 1) write(io,100) re_bin(r), mean_pd(r,1:maxbin_2D)/psd_dbin_2D(1:maxbin_2D)  ![/m]
    if(method == 2) write(io,100) re_bin(r), mean_pd(r,1:maxbin_2D)
 enddo
 close(io)

 100 format(F10.3,33E10.3)
 101 format(34F10.3)


#if MPI > 0 
  endif MPI_MASTERPROC
#endif


 end select 

 return
 end subroutine aircraft_zonal_sampling_psd_re

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine bulk_psd_2D_gen( temp, q, re, mu, rho, maxbin, psd_endbin, bulk_psd )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine re-assign SBM PSD from 7 species into the bulk PSD measurable from 
! aircraft 2D probe. This routine assumes SBM particles are either sphere or non-sphere (choise). 
!
! History:
! 03/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: temp  ! temperature [K]
 type ( particle_gen ),intent(in) :: q  ! particle mass concentration [g/m3]
 type ( particle_gen ),intent(in) :: re ! particle effective radius [um]
 type ( particle_gen ),intent(in) :: mu ! dispersion parameter of gamm PSD [-]
 type ( particle_gen ),intent(in) :: rho! bulk effective density [kg/m3] 

 integer,intent(in) :: maxbin                      ! maximum bin number
 real(sdsu_fps),intent(in) :: psd_endbin(maxbin+1) ! psd endbin size [micron]
 real(sdsu_fps),intent(out) :: bulk_psd(maxbin)    ! bulk particle size distributions [1/m3]

 integer :: n_2d, i, n_sbm, n_pdf      ! loop for size bin
 real(sdsu_fps),dimension(maxbin) :: num    !number concentration per volume for a given size  [#/m3]
 real(sdsu_fps),dimension(maxbin) :: n_tot  ! number concentration per volume for a given size  [#/m3]
 real(sdsu_fps) :: rs_cm, re_cm ! starting,ending,and center radius of bulk PSD bin [cm]
 real(sdsu_fps),dimension(maxbin) :: rad_cm, drad_cm

 integer :: ispc  !RAMS spicies index
 real(sdsu_fps) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(sdsu_fps) :: cfvt, pwvt ! alpha_vt & beta_vt

!
! zero out 
!
 n_tot = 0.e0  
 bulk_psd = 0.e0

!
! create size bin
!
 DO_2D: do n_2d = 1, maxbin  !2D probe (max diameter) size-bin loop
    rs_cm = 0.5e0* psd_endbin( n_2d   ) * 1.e-4 
    re_cm = 0.5e0* psd_endbin( n_2d+1 ) * 1.e-4  !starting, ending radius
    rad_cm (n_2d) = 0.5*(rs_cm+re_cm)  !center radius
    drad_cm(n_2d) = re_cm - rs_cm  !bin width [cm]
 enddo DO_2D


! select case(particle_shape)
! case(0) !sphere

 if( aircraft_ice ) then

   if( trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH' .or. &  !toshii found bug
       trim(cloud_microphysics) == 'THOM' ) then !morrison is twomoment
     call num_gen('ice', temp, maxbin, rad_cm, drad_cm, q%cloud, re%cloud, mu%cloud, rho%cloud, num)
   else
     ispc = 3  ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 --> use RAMS psd
                cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! but using gen canting angle.
     call num_rams( 'ice1', maxbin, rad_cm, drad_cm, q%ice, re%ice , gnu_rams%ice1, &
                     cfvt, pwvt, cfmas, pwmas, num )  
   endif
   n_tot = n_tot + num  ! [#/m3]

   call num_gen('snow',    temp, maxbin, rad_cm, drad_cm, q%snow, re%snow, mu%snow, rho%snow, num)
   n_tot = n_tot + num  ! [#/m3]
   call num_gen('graupel', temp, maxbin, rad_cm, drad_cm, q%graupel, re%graupel, mu%graupel, rho%graupel, num)
   n_tot = n_tot + num  ! [#/m3]
   call num_gen('hail',    temp, maxbin, rad_cm, drad_cm, q%hail, re%hail, mu%hail, rho%hail, num)
   n_tot = n_tot + num  ! [#/m3]
   
 else

   if( trim(cloud_microphysics) == 'WDM') then
     call num_gen('cloud', temp, maxbin, rad_cm, drad_cm, q%cloud, re%cloud, mu%cloud, rho%cloud, num) 
     n_tot = n_tot + num  ! [#/m3]
   else
     ispc = 1  ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 --> use RAMS psd
                cfvt  = rams_dstprms(4,ispc)  ; pwvt  = rams_dstprms(5,ispc)  ! but using gen canting angle.
     call num_rams( 'cloud1', maxbin, rad_cm, drad_cm, q%cloud, re%cloud , gnu_rams%cloud1, &
                     cfvt, pwvt, cfmas, pwmas, num )
   endif

   call num_gen('rain', temp, maxbin, rad_cm, drad_cm, q%rain, re%rain, mu%rain, rho%rain, num)
   n_tot = n_tot + num  ! [#/m3]

 endif

! case(1) !irregular assumption from SnowFake and 2DVD
! case default
!   call stop_sdsu(myrank,'MSG bulk_psd_2D_gen: There is no such particle_shape')
! end select

!
! output
!
  bulk_psd =  n_tot ! bulk particle size distributions [#/m3] 

 return
 end subroutine bulk_psd_2D_gen

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine bulk_psd_2D_sbm( n, frime, fmelt, maxbin, psd_endbin, bulk_psd )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine re-assign SBM PSD from 7 species into the bulk PSD measurable from 
! aircraft 2D probe. This routine assumes SBM particles are either sphere or non-sphere (choise). 
!
! History:
! 03/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 type ( particle_sbm ),intent(in) :: n(nbin)       ! particle # concentration [1/m4]
 integer,intent(in) :: maxbin                      ! maximum bin number
 real(sdsu_fps),intent(in) :: psd_endbin(maxbin+1) ! psd endbin size [micron]
 real(sdsu_fps),intent(in) :: frime(1:nbin)  !rime fraction for snow [-]
 real(sdsu_fps),intent(in) :: fmelt(1:nbin)  !melt fraction for snow [-]

 real(sdsu_fps),intent(out) :: bulk_psd(maxbin)    ! bulk particle size distributions [1/m3]

 integer :: n_2d, i, n_sbm, n_pdf      ! loop for size bin
 real(sdsu_fps) :: n_tot        ! number concentration per volume for a given size  [#/m3]
 real(sdsu_fps) :: rs_cm, re_cm ! starting and ending radius of bulk PSD bin [cm]

 real(sdsu_fps) :: den_diag ! diagnostic density [g/cm3]
 real(sdsu_fps) :: rad_diag ! diagnostic radisu [cm]
 real(sdsu_fps) :: fv_melt  ! volumetric fraction of melting medium [-]
 type ( particle_sbm ) :: num(nbin)  ! particle # concentration [1/m4]

!
! filter ice or liquid sampling
!
 if( aircraft_ice ) then
  num(:)%liq     = 0.
  num(:)%ice_col = n(:)%ice_col
  num(:)%ice_pla = n(:)%ice_pla
  num(:)%ice_den = n(:)%ice_den
  num(:)%snow    = n(:)%snow
  num(:)%graupel = n(:)%graupel
  num(:)%hail    = n(:)%hail
 else
  num(:)%liq     = n(:)%liq
  num(:)%ice_col = 0.
  num(:)%ice_pla = 0.
  num(:)%ice_den = 0.
  num(:)%snow    = 0.
  num(:)%graupel = 0.
  num(:)%hail    = 0.
 endif
 

 bulk_psd = 0.e0
 DO_2D: do n_2d = 1, maxbin  !2D probe (max diameter) size-bin loop

    rs_cm = 0.5e0* psd_endbin( n_2d ) * 1.e-4 ; re_cm = 0.5e0* psd_endbin( n_2d+1 ) * 1.e-4  !starting, ending radius

    n_tot = 0.e0  !initialize

    select case(particle_shape) 
    case(0) !sphere

       DO_SBM: do n_sbm = 1, nbin  !SBM size-bin (bulk radius) loop
          if( rad_sbm(n_sbm)%liq >= rs_cm .and. rad_sbm(n_sbm)%liq < re_cm )         &        
              n_tot = n_tot + num(n_sbm)%liq * drad_sbm(n_sbm)%liq*1.e-2               ! liquid [#/m3]
          if( rad_sbm(n_sbm)%ice_col >= rs_cm .and. rad_sbm(n_sbm)%ice_col < re_cm ) &
              n_tot = n_tot + num(n_sbm)%ice_col * drad_sbm(n_sbm)%ice_col*1.e-2       ! ice column [#/m3]
          if( rad_sbm(n_sbm)%ice_pla >= rs_cm .and. rad_sbm(n_sbm)%ice_pla < re_cm ) &
              n_tot = n_tot + num(n_sbm)%ice_pla * drad_sbm(n_sbm)%ice_pla*1.e-2       ! ice plate [#/m3]
          if( rad_sbm(n_sbm)%ice_den >= rs_cm .and. rad_sbm(n_sbm)%ice_den < re_cm ) &
              n_tot = n_tot + num(n_sbm)%ice_den * drad_sbm(n_sbm)%ice_den*1.e-2       ! ice dendrite [#/m3]

          if(frime(n_sbm) > 0.01 .or. fmelt(n_sbm) > 0.01 ) then  !at least 1%...
             !
             ! Modify bulk radius and density for a given rimed fraction
             !
             call diag_rime_melt_snow (n_sbm, frime(n_sbm), fmelt(n_sbm), den_sbm(1:nbin)%snow, &
                                       x_sbm(1:nbin)%snow, rad_sbm(1:nbin)%snow, &
                                       den_diag, rad_diag, fv_melt )

             !call diag_rime_snow (n_sbm, frime(n_sbm), den_sbm(1:nbin)%snow, &
             !                     x_sbm(1:nbin)%snow, rad_sbm(1:nbin)%snow, &
             !                     den_diag, rad_diag )

             if( rad_diag >= rs_cm .and. rad_diag < re_cm )       &
                 n_tot = n_tot + num(n_sbm)%snow * drad_sbm(n_sbm)%snow*1.e-2             ! rimed snow [#/m3]
          else
             if( rad_sbm(n_sbm)%snow >= rs_cm .and. rad_sbm(n_sbm)%snow < re_cm )       &
                 n_tot = n_tot + num(n_sbm)%snow * drad_sbm(n_sbm)%snow*1.e-2             ! non-rimed snow [#/m3]
          endif


          if( rad_sbm(n_sbm)%graupel >= rs_cm .and. rad_sbm(n_sbm)%graupel < re_cm ) &
              n_tot = n_tot + num(n_sbm)%graupel * drad_sbm(n_sbm)%graupel*1.e-2       ! graupel [#/m3]
          if( rad_sbm(n_sbm)%hail >= rs_cm .and. rad_sbm(n_sbm)%hail < re_cm )       &
              n_tot = n_tot + num(n_sbm)%hail * drad_sbm(n_sbm)%hail*1.e-2             ! hail [#/m3]
       enddo DO_SBM

    case(1) !irregular assumption from SnowFake and 2DVD

       DO_SBM2: do n_sbm = 1, nbin  !SBM size-bin (bulk radius) loop 

          if(frime(n_sbm) > 0.0 ) then
             !
             ! Modify bulk radius and density for a given rimed fraction
             !
             call diag_rime_melt_snow (n_sbm, frime(n_sbm), fmelt(n_sbm), den_sbm(1:nbin)%snow, &
                                       x_sbm(1:nbin)%snow, rad_sbm(1:nbin)%snow, &
                                       den_diag, rad_diag, fv_melt )
!             call diag_rime_snow (n_sbm, frime(n_sbm), den_sbm(1:nbin)%snow, &
!                                  x_sbm(1:nbin)%snow, rad_sbm(1:nbin)%snow, &
!                                  den_diag, rad_diag )
          else
             rad_diag = rad_sbm(n_sbm)%snow ! [cm]
          endif

          !
          ! Non-sphere speicies rely on probability density distribution of Dmax / Rsphere ratio
          ! to estimate Dmax from SBM's sphere particles. 
          !

          call Dmax_Dsphere_PD( rs_cm,re_cm, num(n_sbm)%snow, rad_diag, drad_sbm(n_sbm)%snow, &
                                ratio_sbm%snow, pd_sbm_2D%snow, dsph_sbm%snow, n_tot )   ! Snow aggregate


          !
          ! Sphere species
          !
          if( rad_sbm(n_sbm)%liq >= rs_cm .and. rad_sbm(n_sbm)%liq < re_cm )         &
              n_tot = n_tot + num(n_sbm)%liq * drad_sbm(n_sbm)%liq*1.e-2               ! liquid [#/m3]
          if( rad_sbm(n_sbm)%graupel >= rs_cm .and. rad_sbm(n_sbm)%graupel < re_cm ) &
              n_tot = n_tot + num(n_sbm)%graupel * drad_sbm(n_sbm)%graupel*1.e-2       ! graupel [#/m3]
          if( rad_sbm(n_sbm)%hail >= rs_cm .and. rad_sbm(n_sbm)%hail < re_cm )       &
              n_tot = n_tot + num(n_sbm)%hail * drad_sbm(n_sbm)%hail*1.e-2             ! hail [#/m3]

          ! These species will be added lator into non-spherical category 
          if( rad_sbm(n_sbm)%ice_col >= rs_cm .and. rad_sbm(n_sbm)%ice_col < re_cm ) &
              n_tot = n_tot + num(n_sbm)%ice_col * drad_sbm(n_sbm)%ice_col*1.e-2       ! ice column [#/m3]
          if( rad_sbm(n_sbm)%ice_pla >= rs_cm .and. rad_sbm(n_sbm)%ice_pla < re_cm ) &
              n_tot = n_tot + num(n_sbm)%ice_pla * drad_sbm(n_sbm)%ice_pla*1.e-2       ! ice plate [#/m3]
          if( rad_sbm(n_sbm)%ice_den >= rs_cm .and. rad_sbm(n_sbm)%ice_den < re_cm ) &
              n_tot = n_tot + num(n_sbm)%ice_den * drad_sbm(n_sbm)%ice_den*1.e-2       ! ice dendrite [#/m3]

       enddo DO_SBM2


    case default
       call stop_sdsu(myrank,'MSG bulk_psd_2D_sbm: There is no such particle_shape')
    end select 
     
   bulk_psd( n_2d ) =  n_tot ! bulk particle size distributions [#/m3] 

 enddo DO_2D


 return
 end subroutine bulk_psd_2D_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine Dmax_Dsphere_PD( rs_cm,re_cm, n, rad, drad, &  ! Snow aggregate
                             ratio, pd, dsph, n_tot )   
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Essentially, HUCM_SBM (and most of other scheme) assumes sphere shape of particle, 
!  while an in-situ 2D probe measures edge-to-edge length (maximum diameter of randomlly oriented
!  solid particle). So, in order to take account for irregular shape of ice-phas particles, 
!  we use probability density tables of Dmax-Dsph ratio a) measured from in-situ 2DVD (snow aggregation
!  and b) simulated from SnowFake model (ice crystal habit).
!  
!
! History:
! 06/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: rs_cm,re_cm              ! starting and ending radius of bulk PSD bin [cm]
 real(sdsu_fps),intent(in) :: n                        ! particle # concentration [1/m4]
 real(sdsu_fps),intent(in) :: rad                      ! sphere radius of particle [cm]
 real(sdsu_fps),intent(in) :: drad                     ! width of radius bin [cm]
 real(sdsu_fps),intent(in) :: ratio(nmax_ratio)        ! ratio between Dmax and Dsphere [-] 
 real(sdsu_fps),intent(in) :: pd(nmax_ratio,nmax_dsph) ! probability density of Dmax-Dsphere ratio [0~1]
 real(sdsu_fps),intent(in) :: dsph(nmax_dsph)          ! Dsphere [mm]
 real(sdsu_fps),intent(inout) :: n_tot                 ! number concentration per volume for a given size  [#/m3] 

 integer :: index_dsph     ! index of Dsph
 integer :: l_ratio,l_rad  ! for looping
 real(sdsu_fps) :: n_sub   ! subcomponent of particle # concentration [1/m3]
 real(sdsu_fps) :: pdf     ! probability  (0~1)
 real(sdsu_fps) :: Dsphere ! sphere diameter [mm]
 real(sdsu_fps) :: Dmax    ! edge-to-edge maximum diameter of randomlly oriented solid particle.  [mm]
 real(sdsu_fps) :: ds_mm,de_mm ! starting and ending diameter size bin [mm]

!
! if particle density is zero return
!
! if( n == 0.e0 ) return

!
! unit and radious-diameter conversion
!
 ds_mm = 2.e0*rs_cm*10.e0   ; de_mm = 2.e0*re_cm*10.e0  ![mm]


!
! determine index for Dsphere bin
!

 if    ( 2.e0*rad*10.e0 <= (dsph(1)+dsph(2))*0.5e0 ) then
     index_dsph =  1
 elseif( 2.e0*rad*10.e0 >  (dsph(nmax_dsph-1)+dsph(nmax_dsph))*0.5e0 ) then
     index_dsph = 10
 else

   do l_rad = 2, nmax_dsph-1
      if( 2.e0*rad*10.e0 >  (dsph(l_rad-1)+dsph(l_rad  ))*0.5e0 .and. &
          2.e0*rad*10.e0 <= (dsph(l_rad  )+dsph(l_rad+1))*0.5e0 ) then
          index_dsph = l_rad
          exit
      endif     
   enddo 

 endif

 
!
! Loop over probability density of Dmax-Dsphere ratio. 
!
 do l_ratio = 1, nmax_ratio

    pdf =  pd(l_ratio, index_dsph) !probability density for a given Dsphere index
    if(pdf == 0.e0) cycle  

    Dsphere = 2.e0*rad*10.e0  ![mm]
    Dmax    = Dsphere * ratio(l_ratio) ![mm]
    if( Dmax >= ds_mm .and. Dmax < de_mm )  then
        n_sub = n * drad*1.e-2  * pdf   ! Nconc for a given pd  [#/m3]
        n_tot = n_tot + n_sub           ! adding particle conc [#/m3] 
    endif

 enddo  

 return
 end subroutine Dmax_Dsphere_PD

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


 subroutine WHD_Dsphere_PD( rs_cm,re_cm, n, rad, &         ! Snow aggregate
                            ratio, pd, dsph, n_tot )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Essentially, HUCM_SBM (and most of other scheme) assumes sphere shape of particle, 
!  while an in-situ 2D probe measures edge-to-edge length (maximum diameter of randomlly oriented
!  solid particle). So, in order to take account for irregular shape of ice-phas particles, 
!  we use probability density tables of Dmax-Dsph ratio a) measured from in-situ 2DVD (snow aggregation
!  and b) simulated from SnowFake model (ice crystal habit).
!  
!
! History:
! 06/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: rs_cm,re_cm              ! starting and ending radius of bulk PSD bin [cm]
 real(sdsu_fps),intent(in) :: n                        ! particle # concentration [#/m3] (different from aircraft)
 real(sdsu_fps),intent(in) :: rad                      ! sphere radius of particle [cm]
 real(sdsu_fps),intent(in) :: ratio(nmax_ratio)        ! ratio between Dmax and Dsphere [-] 
 real(sdsu_fps),intent(in) :: pd(nmax_ratio,nmax_dsph) ! probability density of Dmax-Dsphere ratio [0~1]
 real(sdsu_fps),intent(in) :: dsph(nmax_dsph)          ! Dsphere [mm]
 real(sdsu_fps),intent(inout) :: n_tot  ![#/m3]

 integer :: index_dsph     ! index of Dsph
 integer :: l_ratio,l_rad  ! for looping
 real(sdsu_fps) :: pdf     ! probability  (0~1)
 real(sdsu_fps) :: Dsphere ! sphere diameter [mm]
 real(sdsu_fps) :: Dmax    ! edge-to-edge maximum diameter of randomlly oriented solid particle.  [mm]
 real(sdsu_fps) :: ds_mm,de_mm ! starting and ending diameter size bin [mm]

!
! if particle density is zero return
!
! if( n == 0.e0 ) return

!
! unit and radious-diameter conversion
!
 ds_mm = 2.e0*rs_cm*10.e0   ; de_mm = 2.e0*re_cm*10.e0  ![mm]


!
! determine index for Dsphere bin
!

 if    ( 2.e0*rad*10.e0 <= (dsph(1)+dsph(2))*0.5e0 ) then
     index_dsph =  1
 elseif( 2.e0*rad*10.e0 >  (dsph(nmax_dsph-1)+dsph(nmax_dsph))*0.5e0 ) then
     index_dsph = 10
 else

   do l_rad = 2, nmax_dsph-1
      if( 2.e0*rad*10.e0 >  (dsph(l_rad-1)+dsph(l_rad  ))*0.5e0 .and. &
          2.e0*rad*10.e0 <= (dsph(l_rad  )+dsph(l_rad+1))*0.5e0 ) then
          index_dsph = l_rad
          exit
      endif
   enddo

 endif


!
! Loop over probability density of Dmax-Dsphere ratio. 
!
 do l_ratio = 1, nmax_ratio

    pdf =  pd(l_ratio, index_dsph) !probability density for a given Dsphere index
    if(pdf == 0.e0) cycle

    Dsphere = 2.e0*rad*10.e0  ![mm]
    Dmax    = Dsphere * ratio(l_ratio) ![mm]
    if( Dmax >= ds_mm .and. Dmax < de_mm )  then
        n_tot = n_tot + n * pdf  ! adding particle conc for a given pd [#/m3] 
    endif


 enddo

 return

 end subroutine WHD_Dsphere_PD

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine irregular_volume( n, rad, ratio, pd, dsph, vol )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Essentially, HUCM_SBM (and most of other scheme) assumes sphere shape of particle, 
!  while an in-situ 2D probe measures edge-to-edge length (maximum diameter of randomlly oriented
!  solid particle). So, in order to take account for irregular shape of ice-phase particles, 
!  we use probability density tables of Dmax-Dsph ratio a) measured from in-situ 2DVD (snow aggregation
!  and b) simulated from SnowFake model (ice crystal habit).
!  
!
! History:
! 06/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: n                        ! particle # concentration [#/m3] 
 real(sdsu_fps),intent(in) :: rad                      ! sphere radius of particle [cm]
 real(sdsu_fps),intent(in) :: ratio(nmax_ratio)        ! ratio between Dmax and Dsphere [-] 
 real(sdsu_fps),intent(in) :: pd(nmax_ratio,nmax_dsph) ! probability density of Dmax-Dsphere ratio [0~1]
 real(sdsu_fps),intent(in) :: dsph(nmax_dsph)          ! Dsphere [mm]
 real(sdsu_fps),intent(out) :: vol  ! irregular-shape volume [cm3 / m3]

 integer :: index_dsph     ! index of Dsph
 integer :: l_ratio,l_rad  ! for looping
 real(sdsu_fps) :: pdf     ! probability  (0~1)
 real(sdsu_fps) :: Dsphere ! sphere diameter [mm]
 real(sdsu_fps) :: Dmax    ! edge-to-edge maximum diameter of randomlly oriented solid particle.  [mm]

!
! determine index for Dsphere bin
!

 if    ( 2.e0*rad*10.e0 <= (dsph(1)+dsph(2))*0.5e0 ) then
     index_dsph =  1
 elseif( 2.e0*rad*10.e0 >  (dsph(nmax_dsph-1)+dsph(nmax_dsph))*0.5e0 ) then
     index_dsph = 10
 else

   do l_rad = 2, nmax_dsph-1
      if( 2.e0*rad*10.e0 >  (dsph(l_rad-1)+dsph(l_rad  ))*0.5e0 .and. &
          2.e0*rad*10.e0 <= (dsph(l_rad  )+dsph(l_rad+1))*0.5e0 ) then
          index_dsph = l_rad
          exit
      endif
   enddo

 endif


!
! Loop over probability density of Dmax-Dsphere ratio. 
!
 vol = 0.e0 

 do l_ratio = 1, nmax_ratio

    pdf =  pd(l_ratio, index_dsph) !probability density for a given Dsphere index
    if(pdf == 0.e0) cycle

    Dsphere = 2.e0*rad*10.e0  ![mm]
    Dmax    = Dsphere * ratio(l_ratio) ![mm]


    !In PARSIVEL lowest two size bins are not used at all due to low signla-to-noise ratio. 
    if(Dmax >= psd_endbin_parsivel(3)*1.e-3 ) then
       vol = vol +  const_pi/6.e0 * ( (Dmax*1.e-1)**3) *  n * pdf  ! accumulating particle volume [cm3 / m3]
    endif
 
 enddo



 return

 end subroutine irregular_volume

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine bulk_re_vol( maxbin, psd_endbin, bulk_psd, bulk_re, bulk_vol )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine computes bulk total particle volume and effective radius, using 
! size bins of aircraft 2D analysis.  
!
! History:
! 03/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 integer,intent(in) :: maxbin                      ! maximum bin number
 real(sdsu_fps),intent(in) :: psd_endbin(maxbin+1) ! psd endbin size [micron]
 real(sdsu_fps),intent(in) :: bulk_psd(maxbin)     ! bulk particle size distributions [1/m3]
 real(sdsu_fps),intent(out) :: bulk_re   ! bulk effective radius [micron]
 real(sdsu_fps),intent(out) :: bulk_vol  ! bulk particle volume  [cm3/m3]
 integer :: p  !loop
 real(sdsu_fps) :: r, dD, Nt, V
 real(sdsu_fps) :: volume  ! particle volume [cm3/m3]
 real(sdsu_fps) :: third_moment, sec_moment  ! [micron2/m3]

 volume = 0.e0 ; third_moment = 0.e0 ; sec_moment = 0.e0

 do p = 1, maxbin
    dD = (psd_endbin(p+1) - psd_endbin(p))*1.e-6           ! width of diameter size bin [micron]
     r = 0.5e0* (psd_endbin(p+1) + psd_endbin(p))/2.e0     ! radius of max-diameter size bin [micron]
    Nt = bulk_psd(p)                                       ! bulk particle size distributions [1/m3]
     V = 4.e0*const_pi/3.e0 * ( (r*1.e-4)**3)              ! particle volume using a half of maximum diameter [cm3]

            volume =    volume    + V       * Nt   ! [cm3/m3]

      third_moment = third_moment + ( r**3) * Nt   ! [micron3/m3]
        sec_moment =   sec_moment + ( r**2) * Nt   ! [micron2/m3]

 enddo

 bulk_vol = volume  ! bulk particle volume [cm3/m3]

 if( sec_moment <= 0. ) then
     bulk_re = 0.
 else
     bulk_re  = third_moment / sec_moment ! full-size-bin effective radius [micron]
 endif

 return
 end subroutine bulk_re_vol

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine bulk_water_content_gen( q, bulk_wc )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine compute bulk (total) water content compatible to CVI measurements from SBM simulation. 
! Aircraft CVI (Counterflow Virtial Impactor) measures total water content from all types of particles. 
! Thus, bulk microphysics' all species are added for each grid box to get CVI-observable water content. 
!
! History:
! 08/2017  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 type ( particle_gen ),intent(in) :: q ! mixing ratio [g/m3] 
 real(sdsu_fps),intent(out) :: bulk_wc   !bulk water content [g/m3] 

 bulk_wc = 0.e0

 if( aircraft_ice ) then
    bulk_wc = q%ice + q%snow + q%graupel + q%hail
 else
    bulk_wc = q%cloud + q%rain
 endif

 return
 end subroutine bulk_water_content_gen

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine bulk_water_content_rams( q, bulk_wc )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine compute bulk (total) water content compatible to CVI measurements from SBM simulation. 
! Aircraft CVI (Counterflow Virtial Impactor) measures total water content from all types of particles. 
! Thus, RAMS's all species are added for each grid box to get CVI-observable water content. 
!
! History:
! 09/2017  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 type ( particle_rams ),intent(in) :: q ! mixing ratio [g/m3] 
 real(sdsu_fps),intent(out) :: bulk_wc   !bulk water content [g/m3] 

 bulk_wc = 0.e0

 if( aircraft_ice ) then
    bulk_wc = q%ice1 + q%ice2 + q%snow + q%graupel + q%hail
 else
    bulk_wc = q%cloud1 + q%cloud2 + q%rain
 endif

 return
 end subroutine bulk_water_content_rams

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine bulk_water_content_sbm( q, bulk_wc )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine compute bulk (total) water content compatible to CVI measurements from SBM simulation. 
! Aircraft CVI (Counterflow Virtial Impactor) measures total water content from all types of particles. 
! Thus, HUCM_SBM's all species are added for each grid box to get CVI-observable water content. 
!
! History:
! 03/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 type ( particle_sbm ),intent(in) :: q ! mixing ratio [g/m3] 
 real(sdsu_fps),intent(out) :: bulk_wc   !bulk water content [g/m3] 

 bulk_wc = 0.e0

 if( aircraft_ice ) then
    bulk_wc = q%ice_col + q%ice_pla + q%ice_den + q%snow + q%graupel + q%hail
 else
    bulk_wc = q%liq
 endif


 return
 end subroutine bulk_water_content_sbm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_GV_aircraft
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out GV aircraft 3D parameters in GrADS format. 
! 
! History:
! 03/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer i,j,k,l,m,irec,n
 integer,parameter :: io = 602
 integer :: is,ie,js,je,ks,ke

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end

!
!         --------------- GV aircraft parameter  ------------------
!
 if(aircraft_ice) then
   sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.PSD_ICE'//trim(output_suffix)//'.bin'
 else
   sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.PSD_LIQ'//trim(output_suffix)//'.bin'
 endif


#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG write_GV_aircraft ; output for GV aircraft parameters ',&
 ' ->',trim(sdsu_io_file)

!
! Write out binary data 
!
#if MPI == 2
 if(masterproc) &
#endif
 open( unit= io, file = sdsu_io_file, access='direct', &
       status='replace',recl=mxlyr*mxgridx*mxgridy*4 )
 irec = 1

 ! bulk water content comparable to aircraft CVI [g/m3]
 call dump( aircraft(is:ie,js:je,ks:ke)%bulk_wc  , io, irec ) 
 ! bulk effective radius comparable to aircraft 2D probe [micron]
 call dump( aircraft(is:ie,js:je,ks:ke)%bulk_re  , io, irec ) 
 ! bulk particle volume comparable to aircraft 2D probe [cm3/m3]
 call dump( aircraft(is:ie,js:je,ks:ke)%bulk_vol , io, irec ) 
 ! bulk effective density compatile to 2D-CVI analysis [g/cm3]
 call dump( aircraft(is:ie,js:je,ks:ke)%bulk_rho , io, irec ) 
 ! total number conc [#/m3]
 call dump( aircraft(is:ie,js:je,ks:ke)%Ntot     , io, irec )

#if MPI == 2
 if(masterproc) &
#endif
 close(io)

!
! write out grads control file
!
#if MPI == 2
 if(masterproc) then
#endif
 if( write_grads_ctl  ) then

    if(aircraft_ice) then
      sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.PSD_ICE'//trim(output_suffix)//'.ctl'
    else
      sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.PSD_LIQ'//trim(output_suffix)//'.ctl'
    endif

    print*,'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))

    if(aircraft_ice) then
       sdsu_io_file = sdsu_inp_name(1:efile_len)//'.PSD_ICE'//trim(output_suffix)//'.bin'
    else
       sdsu_io_file = sdsu_inp_name(1:efile_len)//'.PSD_LIQ'//trim(output_suffix)//'.bin'
    endif

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  aircraft psd data '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef   ',mxlyr,  ' linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',irec-1
    write(io,*) 'bwc  ',mxlyr,  ' 0  bulk water content [g/m3]'
    write(io,*) 'bre  ',mxlyr,  ' 0  bulk effective radius [micron]'
    write(io,*) 'bvl  ',mxlyr,  ' 0  bulk particle volume [cm3/m3]' 
    write(io,*) 'brho ',mxlyr,  ' 0  bulk effective density [g/cm3]'
    write(io,*) 'ntot ',mxlyr,  ' 0  total number conc [#/m3]'
    write(io,*) 'endvars '
    close(io)
 endif

#if MPI == 2
 endif
#endif


 return
 end subroutine write_GV_aircraft

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine GV_parsivel
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine deal with domain loop to compute GV Parsivel 2D parameters.
!
! subroutines
!
! History:
! 05/2010  Toshi Matsui@NASA GSFC ; Initial.
! 11/2011  Toshi Matsui@NASA GSFC ; update Parsivel simulator (simplified only account for surface). 
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(sdsu_fps) :: bulk_psd(maxbin_parsivel)     ! bulk particle size distributions [#/m3]
 integer :: i,j,k,n,nb
 integer :: ierr
 integer :: is,ie,js,je
 real :: bulk_vol ! not used 
 real(sdsu_fps) :: psd_dbin_parsivel(maxbin_parsivel)  ! [m]

!
! estimate width of size bin
!
 do nb = 1, maxbin_parsivel
    psd_dbin_parsivel(nb) = (psd_endbin_parsivel(nb+1) - psd_endbin_parsivel(nb) ) * 1.e-6  ![m]
 enddo

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end 


 if(masterproc) print*,'MSG GV_parsivel: start GV simulator. '

!
! Initial check 
! 
 if ( cloud_microphysics == 'HUCM_SBM'   .or. &
      cloud_microphysics == 'HUCM_SBM43' ) then
 else
      call stop_sdsu(myrank,'MSG GV_parsivel: cloud_microphysics must be HUCM_SBM or HUCM_SBM43')
 endif

!
! parameter allocation
!
 if( .not. allocated(parsivel) ) then
    allocate( parsivel(is:ie,js:je), &
    stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG GV_parsivel: parsivel allocation error -> Terminate program.')
 endif
 parsivel%parsivel_rain = 0.e0   ;  parsivel%true_rain = 0.e0  ; parsivel%bulk_rho = 0.e0
 parsivel%bulk_re = 0.e0         


 HOR_LOOP: do j = myj_start, myj_end ; do i = myi_start, myi_end

   bulk_psd = 0.e0

   !
   ! Estimate bulk parameters for each grid column, compatible to GV parsivel measurements. 
   ! Size-bin (based on maximum diameter) of particle size distribution (PSD) analysis 
   ! used by Parsivel measurements/analysis.
   !
   call parsivel_surface( 'liq', atmos(i,j,1)%press, q_sbm(i,j,1)%liq, n_sbm(i,j,1,:)%liq, &
                         frime_sbm_snow(i,j,1,:), x_sbm%liq, den_sbm%liq, rad_sbm%liq, &
                         drad_sbm%liq, vt_sbm%liq, parsivel(i,j) , bulk_psd)         ! liquid 

   call parsivel_surface( 'ice_col', atmos(i,j,1)%press, q_sbm(i,j,1)%ice_col, n_sbm(i,j,1,:)%ice_col, &
                         frime_sbm_snow(i,j,1,:),x_sbm%ice_col,den_sbm%ice_col,rad_sbm%ice_col, &
                         drad_sbm%ice_col, vt_sbm%ice_col, parsivel(i,j) , bulk_psd) ! ice column

   call parsivel_surface( 'ice_pla', atmos(i,j,1)%press, q_sbm(i,j,1)%ice_pla, n_sbm(i,j,1,:)%ice_pla,&
                         frime_sbm_snow(i,j,1,:),x_sbm%ice_pla,den_sbm%ice_pla,rad_sbm%ice_pla, &
                         drad_sbm%ice_pla, vt_sbm%ice_pla, parsivel(i,j) , bulk_psd) ! ice plate

   call parsivel_surface( 'ice_den', atmos(i,j,1)%press, q_sbm(i,j,1)%ice_den, n_sbm(i,j,1,:)%ice_den,&
                         frime_sbm_snow(i,j,1,:),x_sbm%ice_den,den_sbm%ice_den,rad_sbm%ice_den, &
                         drad_sbm%ice_den, vt_sbm%ice_den, parsivel(i,j) , bulk_psd) ! ice dendride

   call parsivel_surface( 'snow', atmos(i,j,1)%press, q_sbm(i,j,1)%snow, n_sbm(i,j,1,:)%snow,&
                         frime_sbm_snow(i,j,1,:),x_sbm%snow, den_sbm%snow, rad_sbm%snow, &
                         drad_sbm%snow, vt_sbm%snow , parsivel(i,j) , bulk_psd )     ! snow

   call parsivel_surface( 'graupel', atmos(i,j,1)%press, q_sbm(i,j,1)%graupel, n_sbm(i,j,1,:)%graupel,&
                         frime_sbm_snow(i,j,1,:),x_sbm%graupel,den_sbm%graupel,rad_sbm%graupel, &
                         drad_sbm%graupel, vt_sbm%graupel, parsivel(i,j) , bulk_psd) ! graupel

   call parsivel_surface( 'hail', atmos(i,j,1)%press, q_sbm(i,j,1)%hail, n_sbm(i,j,1,:)%hail,&
                         frime_sbm_snow(i,j,1,:),x_sbm%hail,den_sbm%hail, rad_sbm%hail, &
                         drad_sbm%hail, vt_sbm%hail, parsivel(i,j) , bulk_psd)       ! hail


   if( SUM( bulk_psd(:) ) > 0.e0 ) then


       !
       ! true rain rate [mm/hr]
       !
       parsivel(i,j)%true_rain   = surface(i,j)%rain_rate_liq + surface(i,j)%rain_rate_ice 

       !
       ! total number concentration [#/m3]
       !
       parsivel(i,j)%Ntot = SUM( bulk_psd(:) )

       !
       ! Derive bulk effective radius [micron]
       !
       call bulk_re_vol( maxbin_parsivel, psd_endbin_parsivel, bulk_psd, parsivel(i,j)%bulk_re, bulk_vol )

       !
       ! Estimate bulk effective density
       !
       parsivel(i,j)%bulk_rho = parsivel(i,j)%true_rain / parsivel(i,j)%parsivel_rain

       !
       ! PSD density [#/m4]
       !
       parsivel(i,j)%bulk_psd(:) = bulk_psd(:) /  psd_dbin_parsivel(:) 
       

    else

       parsivel(i,j)%parsivel_rain = undefined   
       parsivel(i,j)%true_rain     = undefined
       parsivel(i,j)%bulk_rho      = undefined 
       parsivel(i,j)%bulk_re       = undefined
       parsivel(i,j)%bulk_psd      = undefined
       parsivel(i,j)%Ntot          = undefined 

    endif

 enddo ; enddo HOR_LOOP

!
! Write GV parsivel parameters 
!
 call write_GV_parsivel

 return
 end subroutine GV_parsivel

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine parsivel_surface( spc, press, q, &
                             num, frime, x, den, rad, &
                             drad, vt  , par , bulk_psd )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Estimate volume flux of condensates at surface through Lagrangian scheme. 
!
! History:
! 05/2010  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 character*(*),intent(in) :: spc            ! species id
 real(sdsu_fps),intent(in) :: press         ! pressure [hPa]
 real(sdsu_fps),intent(in) ::  q            ! total mixing ratio [g/m3]
 real(sdsu_fps),intent(in) :: num  (1:nbin) ! particle # concentration [1/m4]
 real(sdsu_fps),intent(in) :: frime(1:nbin) ! rime fraction for snow [-]
 real(sdsu_fps),intent(in) :: x   (1:nbin)  ! mass per particle [g]
 real(sdsu_fps),intent(in) :: den (1:nbin)  ! density per particle [g/cm3]
 real(sdsu_fps),intent(in) :: rad (1:nbin)  ! radius of particle [cm]
 real(sdsu_fps),intent(in) :: drad(1:nbin)  ! d (radius) [cm]
 real(sdsu_fps),intent(in) :: vt  (1:nbin)  ! terminal velocity LUT [cm/s]
 type ( parsivel_parameter ),intent(inout) :: par
 real(sdsu_fps),intent(inout) :: bulk_psd(maxbin_parsivel)  ! bulk particle size distributions [#/m3]

 integer :: i,j,k,n, np ! for looping
 real(sdsu_fps),parameter :: den_water = 1.e0  ! density of water [g/cm3]
 real(sdsu_fps) :: termvel    ! actual terminal velocity [m/s]
 real(sdsu_fps) :: number_conc
 real(sdsu_fps) :: volume
 real(sdsu_fps) :: rain_fall , number_fall
 real(sdsu_fps) :: number_tot
 real(sdsu_fps) :: D_micron  ![micron]
 real(sdsu_fps) :: rs_cm, re_cm ! starting and ending radius of bulk PSD bin [cm]
 real(sdsu_fps) :: den_diag ! diagnostic density [g/cm3]
 real(sdsu_fps) :: rad_diag ! diagnostic radisu [cm]
 real(sdsu_fps) :: num_snow ! number conc [#/m3]
 real(sdsu_fps) :: wgt      !weight...
 real :: sph_vol

!
! Fiilter small value
!
 if( q <= q_min_condensate ) return


!
! PSD bin loop
!
 SBMBIN_LOOP: do n = 1, nbin  


    !
    ! Modify bulk radius and density for a given rimed fraction
    !
    if(spc=='snow' .and. frime(n) > 0.0) then
      call diag_rime_snow (n, frime(n), den_sbm(1:nbin)%snow, &
                           x_sbm(1:nbin)%snow, rad_sbm(1:nbin)%snow, &
                           den_diag, rad_diag )
      D_micron =   2.e0 * rad_diag * 1.e4  !particle diameter (riming sphere) [micron] 
    else
      D_micron =   2.e0 * rad (n) * 1.e4  !particle diameter (sphere) [micron]  
    endif

    !
    ! Actual particle terminal velocity for a given pressure
    !
    termvel = 1.e-2*vt(n)*sqrt(1000.e0/press)  ! actual terminal velocity for a given pressure [m/s]

    !
    ! Minimum detactable bin is 3rd bin due to low signal-to-noise ratio. 
    !
    if( D_micron >= psd_endbin_parsivel(3) ) then
        number_conc = num (n) * drad(n)*1.e-2                           ! number concentrations [#/m3]
        volume      = 4.e0*const_pi/3.e0 * ( rad(n)**3)  * number_conc  ! sphere volume [cm3 / m3]
    else 
        number_conc = 0.e0
        volume      = 0.e0
    endif

   !
   ! For those that account for liquid particle only, all number bigger than 8mm diameter is 
   ! elliminated (--> Ali Tokay)
   ! 
    if (parsivel_liq_cutoff .and. ( D_micron > psd_endbin_parsivel(23)  )    ) then
       number_conc = 0.e0
       volume      = 0.e0
    endif

    !
    ! Modify volume(k) that take account for irregular shape of snow
    ! 
    IRREGULAR: if(particle_shape == 1 .and. spc=='snow') then !irregular snow case from SnowFake and 2DVD
         sph_vol = volume
         call irregular_volume( number_conc, D_micron/2.e4, ratio_sbm%snow, pd_sbm_parsivel%snow, dsph_sbm%snow, &
                                volume ) 
    endif IRREGULAR

    !
    ! swe_fall is actual rain rate (namely rain rate in liquid phase) that can be measured by Geonor bucket.
    ! rain_fall is rain rate measured by Parsivel. 
    !
    rain_fall = termvel * (volume * den_water * 1.e-3) * 3600. ! PARSIVEL snow volume  [kg/m2] = [mm/hr]
             !   [m/s]   *( [cm3/m3]  *  [g/cm3]  * [kg/g]) * [s]

    number_fall = number_conc  ! number conc [#/m3]
       

    !
    ! Derive bulk PSD and Re compatible to Parsivel
    !
    IF_IRREGULAR: if(particle_shape == 1 .and. spc=='snow') then !irregular case for snow from SnowFake and 2DVD

      ! account for irregular shape using 2DVD statistics from C3VP. 
      do np = 3, maxbin_parsivel  !starting from 3rd bin (since 1st and 2nd bins are always zero in Parsivel)
         rs_cm = 0.5e0* psd_endbin_parsivel( np   ) * 1.e-4  ! starting radius
         re_cm = 0.5e0* psd_endbin_parsivel( np+1 ) * 1.e-4  ! ending radius
         call WHD_Dsphere_PD( rs_cm,re_cm, number_fall, D_micron/2.e4, &
                              ratio_sbm%snow, pd_sbm_parsivel%snow, dsph_sbm%snow, bulk_psd(np) ) ! Snow
      enddo

    else  !sphere case for other species. 

      !
      ! update parsivel_psd (sphere case D_micron is same for all level so just D_micron(1)
      !
      do np = 3, maxbin_parsivel  !starting from 3rd bin (since 1st and 2nd bins are always zero in Parsivel)
         if( D_micron >= psd_endbin_parsivel(np) .and. D_micron < psd_endbin_parsivel(np+1) ) then 
            bulk_psd(np) = bulk_psd(np) + number_fall  ! PSD measured on ground [#/m3]     
         endif
      enddo 

    endif IF_IRREGULAR

    !
    ! output
    !
    par%parsivel_rain = par%parsivel_rain + rain_fall   ! PARSIVEL-observable snow rate (volumetric) [mm]


 enddo SBMBIN_LOOP !SBM

 end subroutine parsivel_surface

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_GV_parsivel
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out GV parsivel 2D parameters in GrADS format. 
! 
! History:
! 05/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer i,j,k,l,m,n,irec
 integer,parameter :: io = 603
 integer :: is,ie,js,je

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end

!
! --------------- Parisvel 2D parameter ------------------ 
!
 sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.GV2D'//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG write_GV_parsivel ; output for GV parsivel parameters ',&
 ' ->',trim(sdsu_io_file)

!
! Write out binary data 
!
#if MPI == 2
 if(masterproc) &
#endif
 open( unit= io, file = sdsu_io_file, access='direct', &
       status='replace',recl=mxgridx*mxgridy*4 )
 irec = 1

 ! rainfall measured by Parsivel (density is assumed to 1.0g/m3) [mm/hr] 
 call dump( parsivel(is:ie,js:je)%parsivel_rain , io, irec )

 ! true surface rainfall (true rainfall rate in liquid phase) [mm/hr]
 call dump( parsivel(is:ie,js:je)%true_rain   , io, irec )

 ! bulk effective density compatile to Parsivel-gauge analysis [g/cm3]
 call dump( parsivel(is:ie,js:je)%bulk_rho      , io, irec )

 ! bulk effective radius compatile to Parsivel measurements [micron]
 call dump( parsivel(is:ie,js:je)%bulk_re       , io, irec )

 ! total number concentration compatile to Parsivel measurements [micron]
 call dump( parsivel(is:ie,js:je)%Ntot       , io, irec )



#if MPI == 2
 if(masterproc) &
#endif
 close(io)

!
! write out grads control file
!
#if MPI == 2
 if(masterproc) then
#endif
 if( write_grads_ctl  ) then

    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.GV2D'//trim(output_suffix)//'.ctl'
    print*,'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.GV2D'//trim(output_suffix)//'.bin'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  GV PARSIVEL data '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef   1 linear 0. 1 '
    write(io,*) 'tdef   1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',5
    write(io,*) 'prain  0 0  Parsivel rainfall [mm/hr]'
    write(io,*) 'train  0 0  true rainfall [mm/hr]'
    write(io,*) 'brho   0 0  bulk effective density [g/cm3]'
    write(io,*) 'bre    0 0  bulk effective radius [micron]'
    write(io,*) 'ntot   0 0  total number conc [#/m3]'
    write(io,*) 'endvars '
    close(io)
 endif

#if MPI == 2
 endif
#endif

!
!          -------------------- 2D PSD -----------------------
!
 if(dump_psd) then

 sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.GV2D_PSD'//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG write_GV_parsivel ; output for GV parsivel PSD parameters ',&
 ' ->',trim(sdsu_io_file)

!
! Write out binary data 
!
#if MPI == 2
 if(masterproc) &
#endif
 open( unit= io, file = sdsu_io_file, access='direct', &
       status='replace',recl=mxgridx*mxgridy*4 )
 irec = 1

 ! bulk PSD compatible to Parsivel measurements [#/m4]
 do n = 1, maxbin_parsivel 
    call dump( parsivel(is:ie,js:je)%bulk_psd(n), io, irec )
 enddo

#if MPI == 2
 if(masterproc) &
#endif
 close(io)

 endif

 return
 end subroutine write_GV_parsivel

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine parsivel_zonal_sampling_psd(proc)
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

 integer :: i,j,k,p, n, r             ! loop indice
 integer,parameter :: maxbin_re=10
 real(sdsu_fps) :: Ntot         ! total number concentration [#/m3]
 real(sdsu_fps) :: pd(maxbin_parsivel)  ! PSD [#/m3]
 real(sdsu_fps),dimension(1:maxbin_re,1:maxbin_parsivel),save :: matrix, mean_pd  !sum of PDF 
 real(sdsu_fps),dimension(1:maxbin_parsivel),save :: matrix_tot, mean_pd_tot  !sum of PDF for total re bin

 real(sdsu_fps),allocatable,save :: matrix_sbm(:), mean_pd_sbm(:)

 real(sdsu_fps),save :: cnt(1:maxbin_re) , cnt_tot !sampling count for each re bin
 
 integer,parameter :: io = 604
 integer,parameter :: method = 2  ! 1 (normalized PSD)  or 2 (PSD)
 real(sdsu_fps) :: psd_dbin_parsivel(maxbin_parsivel)    ! parsivel bin width [m]
 real(sdsu_fps) :: psd_midbin_parsivel(maxbin_parsivel)  ! parsivel mid bin size [micron]

 real(sdsu_fps) :: re_bin(1:maxbin_re) ! re bin (center) [micron]
 real(sdsu_fps),parameter :: re_bin_end(0:maxbin_re) = & ! re bin (boundary) [micron]
                (/ 100.e0, 400., 500., 600., 650., 700., 900., 1200., &
                   2000., 4000., 20000. /)  


 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

!
! set up re bin size
!
 do r = 1, maxbin_re 
    re_bin(r) = 0.5e0 * ( re_bin_end(r-1) + re_bin_end(r) )
 enddo

 matrix = 0.e0
 cnt = 0.e0
 matrix_tot = 0.e0
 cnt_tot = 0. 

!
! Estimate width of size bin
!
 do p = 1, maxbin_parsivel
    psd_dbin_parsivel(p) = (psd_endbin_parsivel(p+1) - psd_endbin_parsivel(p) ) * 1.e-6  !bin width [m]
 enddo

 do p = 1, maxbin_parsivel
    psd_midbin_parsivel(p) = 0.5e0*(psd_endbin_parsivel(p+1) + psd_endbin_parsivel(p) ) ! mid bin size [micron]
 enddo

 nbin = 43  !toshii needed to define here, since nbin is defined after this initialization. 
 if( .not. allocated(matrix_sbm)  ) allocate( matrix_sbm(nbin) )
 if( .not. allocated(mean_pd_sbm) ) allocate( mean_pd_sbm(nbin) )   !

 matrix_sbm = 0.  ; mean_pd_sbm  = 0. 



 case('add') !-------------------------------------------------------------------------------------

!
! estimate normalized PSD 
!
 DOM_LOOP: do j = myj_start, myj_end ; do i = myi_start, myi_end

    !
    ! sampling good values within defined zone
    !
    if( surface(i,j)%lat >= minlat_gv .and. surface(i,j)%lat <= maxlat_gv .and.&
        surface(i,j)%lon >= minlon_gv .and. surface(i,j)%lon <= maxlon_gv .and.&
        parsivel(i,j)%parsivel_rain > 15.e0         ) then !selected zone (toshii convective >15mm/hr

       Ntot = 0. 
       do p = 1, maxbin_parsivel
          Ntot = Ntot +   parsivel(i,j)%bulk_psd(p) * psd_dbin_parsivel(p)   ![#/m3]  
       enddo

       pd = 0.e0
       do p = 1, maxbin_parsivel
          pd(p) = parsivel(i,j)%bulk_psd(p) * psd_dbin_parsivel(p) / Ntot  ! fraction [-]  
       enddo

       re_loop: do r = 1, maxbin_re 
          if( parsivel(i,j)%bulk_re >= re_bin_end(r-1) .and. parsivel(i,j)%bulk_re < re_bin_end(r) ) then
          
              if(method == 1)  matrix(r,1:maxbin_parsivel) = matrix(r,1:maxbin_parsivel) + pd(1:maxbin_parsivel)  !normalized PSD [-] 
              if(method == 2)  matrix(r,1:maxbin_parsivel) = matrix(r,1:maxbin_parsivel) + parsivel(i,j)%bulk_psd(1:maxbin_parsivel)  !PSD [1/m4]
              cnt(r)    = cnt(r) + 1.e0  !sampling count
              exit re_loop

          endif
       enddo re_loop

       !
       ! total
       !
       if(method == 1)  matrix_tot(1:maxbin_parsivel) = matrix_tot(1:maxbin_parsivel) + pd(1:maxbin_parsivel)  !normalized PSD [-] 
       if(method == 2)  matrix_tot(1:maxbin_parsivel) = matrix_tot(1:maxbin_parsivel) + parsivel(i,j)%bulk_psd(1:maxbin_parsivel)  !PSD [1/m4]
       cnt_tot    = cnt_tot + 1.e0  !sampling count

       !
       ! direct sbm output (rain and hail only) [1/m4] 
       !
       if( trim(cloud_microphysics)  == 'HUCM_SBM43' ) &
       matrix_sbm(1:nbin) = matrix_sbm(1:nbin) + n_sbm(i,j,1,1:nbin)%liq &
                                               + n_sbm(i,j,1,1:nbin)%hail             
 

    endif

 enddo ; enddo DOM_LOOP


 case('finish') !---------------------------------------------------------------------------------------


#if MPI > 0 
  call mpi_sdsu_reduce('SUM', matrix  ) 
  call mpi_sdsu_reduce('SUM', cnt    )
  call mpi_sdsu_reduce('SUM', matrix_tot )
  call mpi_sdsu_reduce('SUM', cnt_tot   )

  if( trim(cloud_microphysics)  == 'HUCM_SBM43' ) &
  call mpi_sdsu_reduce('SUM', matrix_sbm )


  MPI_MASTERPROC: if(masterproc) then
#endif

!
! compute mean value for each re bin
!
 do r = 1, maxbin_re
    mean_pd(r,1:maxbin_parsivel) = matrix(r,1:maxbin_parsivel) / cnt(r)
 enddo

 mean_pd_tot(1:maxbin_parsivel) = matrix_tot(1:maxbin_parsivel) / cnt_tot

 mean_pd_sbm(1:nbin) = matrix_sbm(1:nbin) / cnt_tot   ![1/m4]

 print*,'total sampling number: cnt_tot=',cnt_tot

!
! write output for each re bin
!

 sdsu_io_file =trim(sdsu_dir_output)//'PARSIVEL_PSD_COMPOSITE'//trim(output_suffix)//'.txt'

 print*,'MSG parsivel_zonal_sampling_psd: dump PSD composite. See --> ',trim(sdsu_io_file)
 open( io,file=trim(sdsu_io_file) )

 write(io,101) 0.0, psd_midbin_parsivel(1:maxbin_parsivel)
 do r = 1, maxbin_re
    if(method == 1) write(io,100) re_bin(r), mean_pd(r,1:maxbin_parsivel)/psd_dbin_parsivel(1:maxbin_parsivel)  ![/m]
    if(method == 2) write(io,100) re_bin(r), mean_pd(r,1:maxbin_parsivel)  ! [1/m4]
 enddo
 close(io)

 100 format(F10.3,33E10.3)
 101 format(34F10.3)


!
! write output for all re bin
!

 sdsu_io_file =trim(sdsu_dir_output)//'PARSIVEL_PSD'//trim(output_suffix)//'.txt'

 print*,'MSG parsivel_zonal_sampling_psd: dump PSD . See --> ',trim(sdsu_io_file)
 open( io,file=trim(sdsu_io_file) )

 write(io,102) psd_midbin_parsivel(1:maxbin_parsivel)

 if(method == 1) write(io,103) mean_pd_tot(1:maxbin_parsivel)/psd_dbin_parsivel(1:maxbin_parsivel)  ![/m]
 if(method == 2) write(io,103) mean_pd_tot(1:maxbin_parsivel)  ![1/m4]
 close(io)

 102 format(33F10.3)
 103 format(33E10.3)



!
! write output of direct SBM bin
!
 if( trim(cloud_microphysics)  == 'HUCM_SBM43' ) then

 sdsu_io_file =trim(sdsu_dir_output)//'PARSIVEL_SBM_PSD'//trim(output_suffix)//'.txt'

 print*,'MSG parsivel_zonal_sampling_psd: dump SBM PSD . See --> ',trim(sdsu_io_file)
 open( io,file=trim(sdsu_io_file) )

 write(io,104) rad_sbm(1:nbin)%liq * 2.e+1  !diameter [mm]

 write(io,105) mean_pd_sbm(1:nbin)  ![1/m4]
 close(io)

 104 format(43F10.3)
 105 format(43E10.3)

 endif

#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 end select 

 return
 end subroutine parsivel_zonal_sampling_psd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine num_gen(spc, temp, maxbin, rad_cm, drad_cm, wc, re, mu, dens, num)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Derive number concentration for generalized bulk microphysics. 
!
! History:
! 09/2017  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 character(len=*),intent(in) ::spc 
 real(sdsu_fps),intent(in) :: temp !temperature [K]
 integer,intent(in) :: maxbin
 real(sdsu_fps),intent(in) :: rad_cm(maxbin) !radius [cm]
 real(sdsu_fps),intent(in) :: drad_cm(maxbin) !radius [cm]
 real(sdsu_fps),intent(in) :: wc !mass concentration [g/m3]
 real(sdsu_fps),intent(in) :: re !effective radius [um]
 real(sdsu_fps),intent(in) :: mu !dispersion parameter [-]
 real(sdsu_fps),intent(in) :: dens !effective density [kg/m3]
 real(sdsu_fps),intent(out) :: num(maxbin) !  [#/m3]

 real(sdsu_fps) :: gfac1, gfac2, gamfac  !gamma parameter
 real(sdsu_fps) :: lam  ! the slope of the distribution [1/m]
 real(sdsu_fps) :: n0   ! intercept for the exponetial DSD [1/m**4]

 real(sdsu_fps) :: rad_um  !radius [um]
 real(sdsu_fps) :: d_m     !diameter  [um]
 real(sdsu_fps) :: dinc_m  !diameter bin increment [m]
 real(sdsu_fps) :: numd    !number density [#/m4]
 real(sdsu_fps) :: thom_snow_num(maxbin)  !size distributions of thompson scheme [1/m4]

 real(sdsu_fps) :: rhoe !effective density [kg/m3]

 integer :: i

!
! initialize outpuot
!
 num = 0. 

!
! screen bad value
!
 if( wc <= 0. .or. re <= 0. ) return


 if( trim(cloud_microphysics) == 'THOM' .and. &
    trim(spc) == 'snow'    ) then  !thompson spectrum density of snow

    call thom_snow_n_part2 ( wc, temp, maxbin, rad_cm, thom_snow_num )

    do i = 1, maxbin
       dinc_m = drad_cm(i) * 1.e-2 * 2.  ! bin width in diameter [m]
       num(i) = thom_snow_num(i) * dinc_m  ![1/m3]
    enddo
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
! WDM/WSM scheme
!
 if( (trim(cloud_microphysics) == 'WDM' .or. trim(cloud_microphysics) == 'WSM') .and.  &
      trim(spc) == 'snow'  ) then  
      lam = 1. / (2.*re*1.e-6)  !slope [1/m]  Eq. (8) in Bae et al. (2016)
 endif

 rhoe = dens  !effective density 
 if( trim(cloud_microphysics) == 'GMP4ICE' ) &
  rhoe = 1.e+3 * max(0.05, min(0.9e0, 142.36 * ( re**(-0.95) ) ) ) ! re-dependent bulk snow density [kg/m3]
 
!
! derive intercept
!       [g/m3]                     [m3/g]        [1/m]
 n0  = wc       * 6. / const_pi / (rhoe*1000.) * (lam**(4.+mu))  / gfac1  ! [1/m**(4+mu)]

!
! size loop
!
  loop_size: do i=1,maxbin

    rad_um = rad_cm(i) * 1.e4         ! radius [micron] 
    d_m  = rad_cm(i) * 2. * 1e-2      ! diameter [m]
    dinc_m = drad_cm(i) * 1.e-2 * 2.  ! bin width in diameter [m]

!
! chenck size and species
!
    if( trim(spc) == 'cloud' .or. trim(spc) == 'ice' ) then  !non-precipitation
      if( rad_um > 100. ) cycle loop_size
    else  !precipitation
      if( rad_um <= 100. ) cycle loop_size
    endif

!
! Get N(D) : number density 
!         [1/m**(4+mu)]  [ m**mu ]      [1/m]  [m]   = [1/m4]
    numd = n0           * d_m**mu * exp(-lam * d_m )

!
! Get N(D) for WDM double moment species only. 
!
    if( trim(cloud_microphysics) == 'WDM' ) then
       if( trim(spc) == 'cloud' ) call wdm_derive_dsd( 'cloud', wc*1.e-3, re*1.e-6 , d_m, numd )
       if( trim(spc) == 'rain'  ) call wdm_derive_dsd( 'rain' , wc*1.e-3, re*1.e-6 , d_m, numd )
    endif

!
! output     [1/m4] [m]
!
    num(i) = numd * dinc_m  ![1/m3]

 enddo loop_size

 return
 end subroutine num_gen 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine num_rams(spc, maxbin, rad_cm, drad_cm, wc, re, gnu, &
                     cfvt, pwvt, cfmas, pwmas, num)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Derive number concentration for RAMS microphysics. 
!
! History:
! 09/2017  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

 character(len=*),intent(in) ::spc 
 integer,intent(in) :: maxbin
 real(sdsu_fps),intent(in) :: rad_cm(maxbin) !radius [cm]
 real(sdsu_fps),intent(in) :: drad_cm(maxbin) !radius [cm]
 real(sdsu_fps),intent(in) :: wc !mass concentration [g/m3]
 real(sdsu_fps),intent(in) :: re !effective radius [um]
 real(sdsu_fps),intent(in) :: gnu !dispersion parameter [-]
 real(sdsu_fps),intent(in) :: cfvt   ! alpha_vt 
 real(sdsu_fps),intent(in) :: pwvt   ! beta_vt  [-]
 real(sdsu_fps),intent(in) :: cfmas  ! alpha_m [kg /m**beta_m]
 real(sdsu_fps),intent(in) :: pwmas  ! beta_m  [-]
 real(sdsu_fps),intent(out) :: num(maxbin) !  [#/m3]

 real(sdsu_fps) :: gfac1, gfac2, gamfac  !gamma parameter
 real(sdsu_fps) :: lam  ! the slope of the distribution [1/m]
 real(sdsu_fps) :: n0   ! intercept for the exponetial DSD [1/m**4]

 real(sdsu_fps) :: rad_um  !radius [um]
 real(sdsu_fps) :: d_m     !diameter  [um]
 real(sdsu_fps) :: dinc_m  !diameter bin increment [m]
 real(sdsu_fps) :: numd    !number density [#/m4]
 real(sdsu_fps) :: thom_snow_num(maxbin)  !size distributions of thompson scheme [1/m4]

 real(sdsu_fps) :: rhoe !effective density [kg/m3]

 integer :: i

 real(sdsu_fps) ::  &
   dn         ,& ! charactristic diameter [m]
   mean_mass  ,& ! mean mass [kg]
   ntot       ,& ! total particle number concentrations [#/m3]
   fgamma        ! gamma function


!
! initialize outpuot
!
 num = 0. 

!
! screen bad value
!
 if( wc <= 0. .or. re <= 0. ) return


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
! size loop
!
 loop_size: do i=1,maxbin

     d_m  = rad_cm(i) * 2. * 1e-2      ! diameter [m]
     rad_um = rad_cm(i) * 1.e4         ! radius [micron] 
     dinc_m = drad_cm(i) * 1.e-2 * 2.  ! bin width in diameter [m]
     rhoe = (6.0 * cfmas / const_pi) * (d_m ** (pwmas-3.e0)) ! particle density [kg/m3]

!
! chenck size and species
!
    if( trim(spc) == 'cloud1' .or. trim(spc) == 'ice1' ) then  !non-precipitation
      if( rad_um > 100. ) cycle loop_size
    else  !precipitation
      if( rad_um <= 100. ) cycle loop_size
    endif

!
! Gamma function
!
    call gamma_function(gnu,d_m,dn, fgamma)   ! get gamma function
    numd = ntot * fgamma                       ! particle number density[1/m4]


!
! output     [1/m4] [m]
!
    num(i) = numd * dinc_m  ![1/m3]

 enddo loop_size

 return
 end subroutine num_rams

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_GV

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

