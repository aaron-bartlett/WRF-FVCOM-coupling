#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_surface
 use module_simulator
 use module_telsem
 use module_bf_emiss
 use module_mpi
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
!  This module perpares surface optical (albedo/emissivity) properties for variety of active/passive sensor
!  The wavelength including from ultraviolet to microwave. 
! 
! History: 
!  08/2009  Toshi Matsui@NASA GSFC ; Initial    
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
  public :: surface_micro  ,& ! compute microwave surface emissivity
            surface_micro_prps,&! compute microwave surface emissivity for PRPS database
            surface_radar  ,& !         surface sigma zero
            surface_visir  ,& ! compute albedo or emissivity at visible-IR spectrum
            surface_broad  ,&  ! compute albedo or emissivity at each SW/LW broad band 
            surface_sarta    ! compute spectrum emissivty for sarta

!
! module parameters for surface_micro
!
 logical,allocatable,public ::  lambert(:,:)         ! lambertian surface or not
 real(sdsu_fps),allocatable,public ::  emis_micro(:,:,:,:)  ! surface emissivity []
 real(sdsu_fps),allocatable,public ::  ebar_micro(:,:,:)    ! e-bar [] (emissivity averaged for diffuse radiation)

!
! module parameters for surface_micro_prps
!
 real(sdsu_fps),allocatable,public ::  emis_micro_prps(:,:,:,:,:)  ! surface emissivity []
 real(sdsu_fps),allocatable,public ::  ebar_micro_prps(:,:,:,:)    ! e-bar [] (emissivity averaged for diffuse radiation)




!
! visible-IR surface properties
!
 real(sdsu_fps),allocatable,public :: albedo_visir(:,:,:) !spectrum surface albedo
 real(sdsu_fps),allocatable,public ::  emiss_visir(:,:,:) !spectrum surface emissivity
 real(sdsu_fps),allocatable,public ::  emiss_sarta(:,:,:) !spectrum surface emissivity for sarta

!
! broad-band surface properties
!
 real(sdsu_fps),allocatable,public :: albedo_broad(:,:,:,:) !broad-band surface albedo
 real(sdsu_fps),allocatable,public ::  emiss_broad(:,:) !broad-band surface emissivity

!
! MODIS albedo LUT parameters
!
 integer,parameter :: imax_lulc = 16
 integer,parameter :: imax_wave = 7
 integer,parameter :: imax_ndvi = 20
 real(sdsu_fps) :: modis_albedo_lut( 1:imax_lulc , 1:imax_wave , 1:imax_ndvi )  !MODIS-derived albedo lut
 real(sdsu_fps),parameter :: modis_albedo_chan(imax_wave) = &
  (/0.47,0.555,0.659,0.858,1.24,1.64,2.13/)  !wavelengths of MODIS channels [micron]

 real(sdsu_fps) :: modis_albedo_broad( 1:imax_lulc , 1:imax_ndvi ,1:2) !MODIS albedo at visible and NIR


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine surface_micro(umu_grid)
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute microwave surface emissivity.
!  
! History:
! 09/2012  Toshi Matsui@NASA GSFC ; Incidence angle is latitude dependent. 
! 12/2011  Toshi Matsui@NASA GSFC ; add TELSEM emissivity option  
! 08/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in),optional :: umu_grid(:,:) ! cosine of local incident angle [-]
 integer :: i, j, nf, np !looping indice
 integer :: ierr  !allocation stat
 type(atlas_emis_mw) :: atlas
 real(sdsu_fps) :: umu_local
 real(sdsu_fps) :: emis_water(1:2) !water emissivity 1:horizontal 2:vertical [-]
 real(sdsu_fps) :: ebar_water      !water ebar [-]
 real(sdsu_fps) :: srain !surface rain rate [mm/hr]
 real(sdsu_fps) :: vwgt  !vegetation weight [-]
 real(sdsu_fps) :: incidence !incidence angle [rad]
 real :: dummy

 if(masterproc) print*,'MSG surface_micro: Compute surface emissivity'


!
! allocate parameters
!
 if( .not. allocated(lambert)) then
 allocate( lambert(mxgridx,mxgridy), &
           emis_micro(mxgridx,mxgridy,mxfreq_micro,1:2), &
           ebar_micro(mxgridx,mxgridy,mxfreq_micro), &
           stat=ierr )
 if (ierr /= 0) call stop_sdsu(myrank,'MSG surface_micro: allocation error -> Terminate program. ')
 endif
 emis_micro = undefined
 ebar_micro = undefined

!
! open/read TELSEM database
!
 if( land_emiss_micro == 2) then

     sdsu_io_file = trim(sdsu_dir_data)//'telsem_atlas/'

     if(masterproc) print*,'MSG surface_micro: Read TELSEM Emissivity Atlas map of month =',int( sdsu_mm )

! month 7 has problem
     call rttov_readmw_atlas( sdsu_io_file, int( sdsu_mm ) , atlas, ierr ) 
!     call rttov_readmw_atlas( sdsu_io_file, 7 , atlas, ierr ) 

 endif


!
! Compute all the emissivities needed as a function of
! polarization and surface parameter (wind or soil wetness)
!

! frequency loop
 do nf = 1, mxfreq_micro
   
    if(verbose_SDSU .and. masterproc)  &
     print*, ' Frequency = ', freq_micro(nf),' GHz -> Wavelength = ',2.997925e2/freq_micro(nf) ,'mm'

    do j = myj_start, myj_end 

       if( trim(sim_case) == 'MMF_GPROF' .and. j > max_sample_gprof ) cycle

       do i = myi_start, myi_end

          ! 
          ! determine cosine of incidence angle
          !
          if( present( umu_grid ) ) then
             if( umu_grid(i,j) == undefined ) cycle
             umu_local = umu_grid(i,j)     ! cosine of local incident angle 
          else

             !
             ! derive cosine of incidence angle
             !
             if(trim(micro_sensor)== 'GMI' .and. freq_micro(nf) >= 100. ) then
               call derive_incidence_angle( 45.36*const_degrad, sma_micro, surface(i,j)%lat*const_degrad, &
                                          incidence,  dummy )
               umu_local = cos(incidence)  !cosine of sensor incidence angle
             else
               call derive_incidence_angle( ona_angle_micro*const_degrad, sma_micro, surface(i,j)%lat*const_degrad, &
                                          incidence,  dummy )   
               umu_local = cos(incidence)  !cosine of sensor incidence angle
             endif 
          endif



          if(surface(i,j)%iland == 2) then   ! over water

             lambert(i,j)=.false.
             call emit_water_micro(umu_local, freq_micro(nf), surface(i,j)%t_skin, surface(i,j)%u10m, &
                                   emis_micro(i,j,nf,1), emis_micro(i,j,nf,2), ebar_micro(i,j,nf) )

          elseif(surface(i,j)%iland == 1) then ! over land

             select case( land_emiss_micro)
             case(1)  ! NESDIS microwave land-emissivity model (use specular assmption) 

                lambert(i,j)=.false.
                call emit_nesdis_micro( umu_local, freq_micro(nf), &
                           surface(i,j)%h2o_soil, surface(i,j)%frac_veg/100.e0, surface(i,j)%t_soil, &
                           surface(i,j)%t_skin, surface(i,j)%dhgt_snow, & 
                           emis_micro(i,j,nf,1), emis_micro(i,j,nf,2), ebar_micro(i,j,nf))

                !
                ! Special treatment for surface standing water especially for desert surface ... wet surface
                ! Infiltration rate 1~5mm/hr for clay --> use 1mm/hr
                ! so whenever surface rain rate exceed this threhold, 
                ! one may assume standing water exist.
                !
                srain = surface(i,j)%rain_rate_liq + surface(i,j)%rain_rate_liq 
                if( srain > 1.0 ) then
                   call emit_water_micro(umu_local, freq_micro(nf), surface(i,j)%t_skin, surface(i,j)%u10m, &
                                        emis_water(1), emis_water(2), ebar_water )
                    vwgt = surface(i,j)%frac_veg/100.e0  
                    emis_micro(i,j,nf,1:2)= (1.-vwgt)*emis_water(1:2) + vwgt*emis_micro(i,j,nf,1:2)
                    ebar_micro(i,j,nf)    = (1.-vwgt)*ebar_water      + vwgt*ebar_micro(i,j,nf)

!                    call emit_water_micro(umu_local, freq_micro(nf), surface(i,j)%t_skin, surface(i,j)%u10m, &
!                                    emis_micro(i,j,nf,1), emis_micro(i,j,nf,2), ebar_micro(i,j,nf) )
                endif

             case(2)  ! TELSEM emissivity database (use specular assmption) 

                 call emit_telsem_micro( umu_local, freq_micro(nf), surface(i,j)%lat, surface(i,j)%lon,0.01, atlas, & 
                                          emis_micro(i,j,nf,1), emis_micro(i,j,nf,2), ebar_micro(i,j,nf) )

                 if( emis_micro(i,j,nf,1) == 0. ) then !if TELSEM give zero, call NESDIS
                    call emit_nesdis_micro( umu_local, freq_micro(nf), &
                             surface(i,j)%h2o_soil, surface(i,j)%frac_veg/100.e0, surface(i,j)%t_soil, &
                             surface(i,j)%t_skin, surface(i,j)%dhgt_snow, &
                             emis_micro(i,j,nf,1), emis_micro(i,j,nf,2), ebar_micro(i,j,nf) )
                 endif

                !
                ! Special treatment for surface standing water especially for desert surface ... wet surface
                ! Infiltration rate 1~5mm/hr for clay --> use 1mm/hr
                ! so whenever surface rain rate exceed this threhold, 
                ! one may assume standing water exist.
                !
                srain = surface(i,j)%rain_rate_liq + surface(i,j)%rain_rate_liq
                if( srain > 1.0 ) then
                   call emit_water_micro(umu_local, freq_micro(nf), surface(i,j)%t_skin, surface(i,j)%u10m, &
                                        emis_water(1), emis_water(2), ebar_water )
                    vwgt = surface(i,j)%frac_veg/100.e0  
                    emis_micro(i,j,nf,1:2)= (1.-vwgt)*emis_water(1:2) + vwgt*emis_micro(i,j,nf,1:2)
                    ebar_micro(i,j,nf)    = (1.-vwgt)*ebar_water      + vwgt*ebar_micro(i,j,nf)

!                    call emit_water_micro(umu_local, freq_micro(nf), surface(i,j)%t_skin, surface(i,j)%u10m, &
!                                    emis_micro(i,j,nf,1), emis_micro(i,j,nf,2), ebar_micro(i,j,nf) )
                endif



             case(0) ! simple land emissivity scheme (lambertian assumption)

                lambert(i,j)=.true.
                call emit_water_micro(umu_local, freq_micro(nf), surface(i,j)%t_skin, &
                            0., emis_micro(i,j,nf,1), emis_micro(i,j,nf,2), ebar_micro(i,j,nf))

                do np = 1,2 !polarizatin loop
                   emis_micro(i,j,nf,np) = (1.-surface(i,j)%h2o_soil)*0.9 + &
                                               surface(i,j)%h2o_soil*emis_micro(i,j,nf,np)
                enddo !np
                ebar_micro(i,j,nf) = (1.-surface(i,j)%h2o_soil)*0.9+ &
                                       surface(i,j)%h2o_soil*ebar_micro(i,j,nf)

             case default

                print*,'MSG surface_micro: No such land_emiss_micro=', land_emiss_micro
                print*,' It must be 0 ~ 2'
                call stop_sdsu(myrank,'MSG surface_micro: terminate program')

             end select

          else
            print*,surface(i,j)%iland
            call stop_sdsu(myrank,'MSG surface_micro; no such background surface')

          end if

       end do !i
    end do !j
  enddo !nf

!
! close TELSEM database
!
 if( land_emiss_micro == 2) then
    call rttov_closemw_atlas(atlas)
 endif


!
! write out surface emissivity
!
  if(write_surface) call write_surface_micro

#if MPI == 2
if(numproc_tot>1) then
!
! broadcasting to all threads
!
 if( scan_micro .or. scan_radar ) then
  call mpi_sdsu_communicate( 'TO_ALL', emis_micro )
  call mpi_sdsu_communicate( 'TO_ALL', ebar_micro )
 endif
endif
#endif

 return
 end subroutine surface_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine surface_micro_prps
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute microwave surface emissivity.
!  
! History:
! 09/2019  Toshi Matsui@NASA GSFC ; surface_micro is converted for PRPS database. 
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: i, j, nf, np, n !looping indice
 integer :: ierr  !allocation stat
 type(atlas_emis_mw) :: atlas
 real(sdsu_fps),allocatable :: umu_local(:)
 real(sdsu_fps) :: emis_water(1:2) !water emissivity 1:horizontal 2:vertical [-]
 real(sdsu_fps) :: ebar_water      !water ebar [-]
 real(sdsu_fps) :: srain !surface rain rate [mm/hr]
 real(sdsu_fps) :: vwgt  !vegetation weight [-]
 real(sdsu_fps) :: incidence !incidence angle [rad]
 real :: dummy

 if(masterproc) print*,'MSG surface_micro: Compute surface emissivity'

!
! allocate parameters
!
 if( .not. allocated(lambert)) allocate( lambert(mxgridx,mxgridy),stat=ierr )
 if( .not. allocated(emis_micro_prps) ) &
       allocate( emis_micro_prps(mxgridx,mxgridy,mxfreq_micro,1:2,min_index_angle:max_index_angle),stat=ierr ) 
 if( .not. allocated(ebar_micro_prps) ) &
       allocate( ebar_micro_prps(mxgridx,mxgridy,mxfreq_micro,min_index_angle:max_index_angle),stat=ierr )
 if( .not. allocated(umu_local) ) allocate( umu_local      ( min_index_angle:max_index_angle ),stat=ierr )

 if (ierr /= 0) call stop_sdsu(myrank,'MSG surface_micro_prps: allocation error -> Terminate program. ')

 emis_micro_prps = undefined
 ebar_micro_prps = undefined
 umu_local = undefined

!
! open/read TELSEM database
!
 if( land_emiss_micro == 2) then

     sdsu_io_file = trim(sdsu_dir_data)//'telsem_atlas/'

     if(masterproc) print*,'MSG surface_micro: Read TELSEM Emissivity Atlas map of month =',int( sdsu_mm )

     call rttov_readmw_atlas( sdsu_io_file, int( sdsu_mm ) , atlas, ierr ) 

 endif

!
! Compute all the emissivities needed as a function of
! polarization and surface parameter (wind or soil wetness) for all viewing angles
!

! frequency loop
 do nf = 1, mxfreq_micro
   
    if(verbose_SDSU .and. masterproc)  &
     print*, ' Frequency = ', freq_micro(nf),' GHz -> Wavelength = ',2.997925e2/freq_micro(nf) ,'mm'

    do j = myj_start, myj_end  !lat index

       if( j > max_sample_gprof ) cycle

       do i = myi_start, myi_end !lon index

          ! 
          ! determine cosine of incidence angle
          !
          do n = min_index_angle, max_index_angle

             umu_local(n) = cos( inc_angle_prps(i,j,n)*const_degrad )  !cosine of sensor incidence angle

          enddo 

          !
          ! derive emisivity
          !
          if(surface(i,j)%iland == 2) then   ! over water

             lambert(i,j)=.false.
            
             do n = min_index_angle, max_index_angle 

                if( gf_prps(i,j,nf,n) > 0. ) then !compute just within FOV

                   call emit_water_micro(umu_local(n), freq_micro(nf), surface(i,j)%t_skin, surface(i,j)%u10m, &
                              emis_micro_prps(i,j,nf,1,n), emis_micro_prps(i,j,nf,2,n), ebar_micro_prps(i,j,nf,n) )

                endif

             enddo

          elseif(surface(i,j)%iland == 1) then ! over land

             select case( land_emiss_micro)
             case(1)  ! NESDIS microwave land-emissivity model (use specular assmption) 

                lambert(i,j)=.false.

                do n = min_index_angle, max_index_angle

                   if( gf_prps(i,j,nf,n) > 0. ) then !compute just within FOV
                       call emit_nesdis_micro( umu_local(n), freq_micro(nf), &
                                surface(i,j)%h2o_soil, surface(i,j)%frac_veg/100.e0, surface(i,j)%t_soil, &
                                surface(i,j)%t_skin, surface(i,j)%dhgt_snow, & 
                                emis_micro_prps(i,j,nf,1,n), emis_micro_prps(i,j,nf,2,n), ebar_micro_prps(i,j,nf,n))
                   endif
                enddo

                !
                ! Special treatment for surface standing water especially for desert surface ... wet surface
                ! Infiltration rate 1~5mm/hr for clay --> use 1mm/hr
                ! so whenever surface rain rate exceed this threhold, 
                ! one may assume standing water exist.
                !
                srain = surface(i,j)%rain_rate_liq + surface(i,j)%rain_rate_liq 
                if( srain > 1.0 ) then

                   do n = min_index_angle, max_index_angle

                      if( gf_prps(i,j,nf,n) > 0. ) then !compute just within FOV

                        call emit_water_micro(umu_local(n), freq_micro(nf), surface(i,j)%t_skin, surface(i,j)%u10m, &
                                            emis_water(1), emis_water(2), ebar_water )
                        vwgt = surface(i,j)%frac_veg/100.e0  
                        emis_micro_prps(i,j,nf,1:2,n) = (1.-vwgt)*emis_water(1:2) + vwgt*emis_micro_prps(i,j,nf,1:2,n)
                        ebar_micro_prps(i,j,nf,n)     = (1.-vwgt)*ebar_water      + vwgt*ebar_micro_prps(i,j,nf,n)
                      endif
                   enddo

                endif

             case(2)  ! TELSEM emissivity database (use specular assmption) 

                 do n = min_index_angle, max_index_angle
                    if( gf_prps(i,j,nf,n) > 0. ) then !compute just within FOV

                      call emit_telsem_micro( umu_local(n), freq_micro(nf), surface(i,j)%lat, surface(i,j)%lon,0.01, atlas, & 
                                          emis_micro_prps(i,j,nf,1,n), emis_micro_prps(i,j,nf,2,n), ebar_micro_prps(i,j,nf,n) )

                      if( emis_micro(i,j,nf,1) == 0. ) then !if TELSEM give zero, call NESDIS
                          call emit_nesdis_micro( umu_local(n), freq_micro(nf), &
                                surface(i,j)%h2o_soil, surface(i,j)%frac_veg/100.e0, surface(i,j)%t_soil, &
                                surface(i,j)%t_skin, surface(i,j)%dhgt_snow, &
                                emis_micro_prps(i,j,nf,1,n), emis_micro_prps(i,j,nf,2,n), ebar_micro_prps(i,j,nf,n) )
                      endif
                    endif
                 enddo

                !
                ! Special treatment for surface standing water especially for desert surface ... wet surface
                ! Infiltration rate 1~5mm/hr for clay --> use 1mm/hr
                ! so whenever surface rain rate exceed this threhold, 
                ! one may assume standing water exist.
                !
                srain = surface(i,j)%rain_rate_liq + surface(i,j)%rain_rate_liq
                if( srain > 1.0 ) then

                   do n = min_index_angle, max_index_angle
                      if( gf_prps(i,j,nf,n) > 0. ) then !compute just within FOV
                         call emit_water_micro(umu_local(n), freq_micro(nf), surface(i,j)%t_skin, surface(i,j)%u10m, &
                                               emis_water(1), emis_water(2), ebar_water )
                         vwgt = surface(i,j)%frac_veg/100.e0
                         emis_micro_prps(i,j,nf,1:2,n) = (1.-vwgt)*emis_water(1:2) + vwgt*emis_micro_prps(i,j,nf,1:2,n)
                         ebar_micro_prps(i,j,nf,n)     = (1.-vwgt)*ebar_water      + vwgt*ebar_micro_prps(i,j,nf,n)
                      endif
                   enddo

                endif


             case(0) ! simple land emissivity scheme (lambertian assumption)

                lambert(i,j)=.true.

                do n = min_index_angle, max_index_angle
                   if( gf_prps(i,j,nf,n) > 0. ) then !compute just within FOV

                       call emit_water_micro(umu_local(n), freq_micro(nf), surface(i,j)%t_skin, &
                                    0., emis_micro_prps(i,j,nf,1,n), emis_micro_prps(i,j,nf,2,n), ebar_micro_prps(i,j,nf,n))

                      do np = 1,2 !polarizatin loop
                          emis_micro_prps(i,j,nf,np,n) = (1.-surface(i,j)%h2o_soil)*0.9 + &
                                                  surface(i,j)%h2o_soil*emis_micro_prps(i,j,nf,np,n)
                      enddo !np
                      ebar_micro_prps(i,j,nf,n) = (1.-surface(i,j)%h2o_soil)*0.9+ &
                                            surface(i,j)%h2o_soil*ebar_micro_prps(i,j,nf,n)
                   endif
                enddo

             case default

                print*,'MSG surface_micro_prps: No such land_emiss_micro=', land_emiss_micro
                print*,' It must be 0 ~ 2'
                call stop_sdsu(myrank,'MSG surface_micro_prps: terminate program')

             end select

          else
            print*,surface(i,j)%iland
            call stop_sdsu(myrank,'MSG surface_micro_prps; no such background surface')

          end if

       end do !i (lon index)
    end do !j (lat index)
  enddo !nf (freq)

!
! close TELSEM database
!
 if( land_emiss_micro == 2) then
    call rttov_closemw_atlas(atlas)
 endif


!
! write out surface emissivity
!

!#if MPI == 2
!if(numproc_tot>1) then
!
! broadcasting to all threads
!
! if( scan_micro .or. scan_radar ) then
!  call mpi_sdsu_communicate( 'TO_ALL', emis_micro_prps )
!  call mpi_sdsu_communicate( 'TO_ALL', ebar_micro_prps )
! endif
!endif
!#endif

 return
 end subroutine surface_micro_prps

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  emit_water_micro( umu_local, f, ts, w, emish_out, emisv_out , ebar_out )
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute emissivity of water surface at various angles. 
!  
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Adapted to SDSU.
!    ~~    Chris Kummerow@CSU     ; initial
!           
! References: 
!
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: umu_local !consine of local incident angle [-] 
 real(sdsu_fps),intent(in) :: f  ! frequency [GHz]
 real(sdsu_fps),intent(in) :: ts ! surface skin temperature [K]
 real(sdsu_fps),intent(in) :: w  ! surface wind speed [m/s]
 real(sdsu_fps),intent(out) :: emish_out ! emissivity at a given angle (horizontal)
 real(sdsu_fps),intent(out) :: emisv_out ! emissivity at a given angle (vertical)
 real(sdsu_fps),intent(out) :: ebar_out  ! mean emissivity for all angle (horizonta+vertical)

!
! Local parameter
!
 integer :: i !looping
 integer,parameter:: nang = 21 !number of ray angle
 real(sdsu_fps) :: ang(nang), mu(nang), esum(nang)
 real(sdsu_fps) :: dmu, avmu
 real(sdsu_fps) :: sume
 real(sdsu_fps) :: eh, ev ! emissivity for horizontal- and vertical-polarzation
 real(sdsu_fps) :: eavg
 real(sdsu_fps) :: pi
 real(sdsu_fps) :: s !Salinity [PPM]
 real(sdsu_fps) :: angle ! incident angle [degree]
 real(sdsu_fps) :: emish,emisv


 pi = 2.*asin(1.0)
 s = 35.                ! SALINITY IN PPM
! salinity ranges from 30 ~40 g/kg over globe.

 angle = acos(umu_local)*180./pi  !incident angle
!
!     CALCULATE EMIS AT GIVEN ANGLE
!
 call emiss_micro (f,angle,s,ts,w,ev,eh,emish,emisv)

 emish_out = eh
 emisv_out = ev

!
!     CALCULATE EMIS AT VARIOUS ANGLES
!
  do i = 1,nang
     ang(i) = 4.*( i - 1 )
     mu(i) = cos(ang(i)*const_pi/180.)
     call emiss_micro (f,ang(i),s,ts,w,ev,eh,emish,emisv)
     esum(i) =  ev + eh
 enddo 

!
!     CALCULATE EBAR
!
 sume = 0.0
 do i = 1,nang-1
    eavg = 0.5*( esum(i) + esum(i+1) )
    dmu = mu(i) - mu(i+1)
    avmu = 0.5*( mu(i) + mu(i+1) )
    sume = sume + eavg*avmu*dmu
 enddo
 ebar_out = sume
 return

 end subroutine  emit_water_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  emiss_micro (f,angle,s,ts,w,ev,eh,emish,emisv)
  implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute emissivity of water surface. 
!  
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Adapted to SDSU.
!    ~~    Chris Kummerow@CSU     ; initial
!           
! References: 
!
!----------------------------------------------------------------------------------------------------

  real(sdsu_fps),intent(in) :: f     ! frequency [GHz]
  real(sdsu_fps),intent(in) :: angle ! incident angle [deg]
  real(sdsu_fps),intent(in) :: s     ! salinity [PPM]
  real(sdsu_fps),intent(in) :: ts    ! skin temperateur [K]
  real(sdsu_fps),intent(in) :: w     !
  real(sdsu_fps),intent(out) :: eh,ev,emish,emisv

  real(sdsu_fps) :: es,t,e1,e2,theta,cmhu,cmhu2,fact1,b,b2,arg,psiz,cospz,f11,psiy,g,cospy,f22
  real(sdsu_fps) :: rsh,rsv,foam,gh,gv
  real(sdsu_fps) :: a1,rfv,rfh,y,sqrtf,corrv,corrh,rrv,rrh,rv,rh
  real(sdsu_fps),parameter :: dtr =  0.01745329252 


  t = ts-273.16
  theta = angle*dtr
!
  call diecon(s,t,f,e1,e2)
!
  cmhu = cos (theta)
  cmhu2 = cmhu*cmhu
  fact1 = e1+cmhu2-1.
  b = (fact1*fact1+e2*e2)**0.25
  b2 = b*b
  arg = e2/fact1
  psiz = 0.5*atan(arg)
  cospz = cos(psiz)*cmhu*b
  f11 = cmhu2+b2+2.*cospz
  psiy = atan(e2/e1)
  g = cmhu*sqrt(e1*e1+e2*e2)
  cospy = cos(psiy-psiz)*g*b
  f22 = g*g+b2+2.*cospy
!
!     FOR SPECULAR SURFACES THE HORIZONTAL AND VERTICAL EMISSIVITY
!     CALCULATION
!
  emish = 4.*cospz/f11
  emisv = 4.*cospy/f22
!
!     FROM HERE THE EFFECT OF SURFACE ROUGHNESS AND FOAM ARE INCLUDED
!     BASED ON HOLLINGER MODEL AND FOAM MODEL OF STOGRYN
!
  rsh = 1.-emish
  rsv = 1.-emisv
!
!     P0 = 1.707476E-2+8.560329E-4*F+1.120024E-5*F*F
!     P1 = -1.500792E-2+1.820672E-3*F-4.633806E-5*F*F
!     P2 = 2.442217E-4-2.282022E-6*F+4.194212E-7*F*F
!
!     FOAM = (P0+P1*W+P2*W*W)
!

  foam = 7.751e-06 * w ** 3.231
!
  gh = 1.-1.748e-3*angle-7.336e-5*angle**2+1.044e-7*angle**3
  gv = 1.-9.946e-4*angle+3.218e-5*angle**2-1.187e-6*angle**3 &
      +7.0e-20*angle**10
!


  a1 = (208.0+1.29*f)/ts
!
!     RFV = 1.-A1*GV-0.005*F
!     RFH = 1.-A1*GH-0.005*F
!
  rfv = 1. - a1 * gv
  rfh = 1. - a1 * gh
!
  y = 7.32e-02*angle
!
!     TS SURFACE TEMP IS IN DEGREE KELVIN
!
  sqrtf = sqrt(f)
!
!     CORRV = (W*(1.17E-01-2.09E-03*EXP(Y))*SQRTF/TS)-0.00065*F
!     CORRH = (W*(1.15E-01+3.80E-05*ANGLE**2)*SQRTF/TS)-0.00065*F
!
  corrv = (w*(1.17e-01-2.09e-03*exp(y))*sqrtf/ts)
  corrh = (w*(1.15e-01+3.80e-05*angle**2)*sqrtf/ts)
!
  rrv = rsv-corrv
  rrh = rsh-corrh
!
  rv = rrv*(1.-foam)+rfv*foam
  rh = rrh*(1.-foam)+rfh*foam
!
  eh = 1.-rh
  ev = 1.-rv
!
  return
 end subroutine  emiss_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  diecon(s,t,freq,e1,e2)
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute dielectic constant of saline water 
!  
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Adapted to SDSU.
!    ~~    Chris Kummerow@CSU     ; initial
!           
! References: 
!
!----------------------------------------------------------------------------------------------------


 real(sdsu_fps),intent(in) :: freq  ! frequency [GHz]
 real(sdsu_fps),intent(in) :: s     ! salinity [PPM]
 real(sdsu_fps),intent(in) :: t     ! skin temperateur [Cel]

 real(sdsu_fps),intent(out) :: e1,e2

 real(sdsu_fps), save :: st, s2, t2, sst, stt, sstt, es, tau, sigma
 real(sdsu_fps) :: znu,omega,den
 real(sdsu_fps),parameter :: twopi = 6.283185307
 real(sdsu_fps),parameter :: einf  = 4.9
 real(sdsu_fps),save :: sold = 0.0
 real(sdsu_fps),save :: told = -99.

  if (s .eq. sold .and. t .eq. told) go to 10

  st = s*t
  s2 = s*s
  t2 = t*t
  sst = s2*t
  stt = t2*s
  sstt = s2*t2
  es = 88.-4.339e-01*s+1.71e-03*s2-4.035e-01*t+8.065e-04*t2+6.170e-03 &
        * st-8.910e-05*sst-6.934e-05*stt+1.439e-06*sstt
!
  tau = (18.70-7.924e-02*s+6.35e-04*s2-5.489e-01*t+5.758e-03*t2+ &
        1.889e-03*st-7.209e-06*sst-5.299e-07*stt-2.101e-07*sstt)*1.0e-12
!
  sigma = (7.788e-03*s-1.672e-06*s2-8.570e-15*t+2.996e-16*t2+4.059e-04 &
           * st-3.215e-06*sst-1.423e-06*stt+3.229e-08*sstt)*1.0e11
!
10 znu = freq*1.e09
   omega = twopi*znu
   den = 1. + (omega * tau) ** 2
   e1 = (es-einf)/den+einf
   e2 = (es-einf)*omega*tau/den+2.*sigma/znu
   sold = s
   told = t
!
  return
 end subroutine  diecon

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  emit_nesdis_micro( umu_local, f, soilH2O, vfrac, tsoil, tskin, &
                                snow_depth, emish_out, emisv_out , ebar_out )
 use NESDIS_LandEM_Module
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute emissivity of land-surface via NESDIS model.
!  
! History:
! 09/2008  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: umu_local !consine of local incident angle [-] 
!
! input & output parameter
!
  real(sdsu_fps),intent(in) :: &
     f       , & !frequency
     soilH2O , & !volumetric soil moisture [-]
     vfrac   , & !vegetation fraction [-]
     tsoil   , & !top-soil temperature [K]
     tskin   , & !surface skin temperature (including snow)
     snow_depth  !physical snow deph [m] 

  real(sdsu_fps),intent(out) :: &
    emish_out , &! emissivity at given angle (horizontal) [-]
    emisv_out , &! emissivity at given angle (vertical) [-]
    ebar_out    ! mean emissivity averaged over different angles [-]

!
! Local parameter
!
  integer :: i !looping
  integer,parameter:: nang = 21
  real(sdsu_fps) :: ang(nang), mu(nang), esum(nang)
  real(sdsu_fps) :: dmu, avmu
  real(sdsu_fps) :: eh, ev ! emissivity for horizontal- and vertical-polarzation
  real(sdsu_fps) :: eavg
  real(sdsu_fps) :: tot
  real(sdsu_fps) :: s !Salinity [PPM]
  real(sdsu_fps) :: angle !incident angle [degree]

  s = 35.                             ! SALINITY IN PPM
  angle = acos(umu_local)*180./const_pi     ! incident angle [deg]

!
!  CALCULATE EMIS AT GIVEN ANGLE
!


 call NESDIS_LandEM( angle , f , & !input
                     soilH2O, vfrac, tsoil, tskin, snow_depth*1000., & !input   
                     eh, ev ) !output


 emish_out = eh !horizontal polarization
 emisv_out = ev !vertical polarization 

!
!     CALCULATE EMIS AT VARIOUS ANGLES
!
 do i = 1,nang
    ang(i) = 4.*REAL( i - 1 )
    mu(i) = cos(ang(i)*const_pi/180.)
    angle = ang(i)
    call NESDIS_LandEM(angle , f , & !input
                       soilH2O, vfrac, tsoil, tskin, snow_depth*1000., & !input   
                       eh, ev ) !output
    esum(i) =  ev + eh
 enddo

!
!     CALCULATE EBAR
!
 tot = 0.0
 do i = 1,nang-1
    eavg = 0.5*( esum(i) + esum(i+1) )
    dmu = mu(i) - mu(i+1)
    avmu = 0.5*( mu(i) + mu(i+1) )
    tot = tot + eavg*avmu*dmu
 enddo
 ebar_out = tot
 return
 end subroutine  emit_nesdis_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  emit_telsem_micro( umu_local, f, lat, lon, res, atlas, emish_out, emisv_out , ebar_out )

!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute emissivity of land-surface via TELSEM emissiity atlas database.
!  
! History:
! 12/2011  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: umu_local !consine of local incident angle [-] 
  real(sdsu_fps),intent(in) :: &
     f       , & ! frequency [GHz]
     lat     , & ! latutide [deg]
     lon     , & ! longitude [deg]
     res         ! resolution (horizontal grid spcing [deg] 

  type(atlas_emis_mw),intent(in) :: atlas

  real(sdsu_fps),intent(out) :: &
    emish_out , &! emissivity at given angle (horizontal) [-]
    emisv_out , &! emissivity at given angle (vertical) [-]
    ebar_out    ! mean emissivity averaged over different angles [-]

!
! Local parameter
!
  integer :: i !looping
  integer,parameter:: nang = 21
  real(sdsu_fps) :: ang(nang), mu(nang), esum(nang)
  real(sdsu_fps) :: dmu, avmu
  real(sdsu_fps) :: eh, ev ! emissivity for horizontal- and vertical-polarzation
  real(sdsu_fps) :: eavg
  real(sdsu_fps) :: tot
  real(sdsu_fps) :: s !Salinity [PPM]
  real(sdsu_fps) :: angle !incident angle [degree]
  real(sdsu_fps) :: stdv, stdh
  real(sdsu_fps) :: lat_in
  integer :: verb  !=1 for TRUE and 0 for FALSE
  logical :: bilinear = .false.  !if true -> spatially disaggregate TELSEM 0.25 emissivity map

  verb=0 !0 for none, 1 for full

  angle = acos(umu_local)*180./const_pi     ! incident angle [deg]

  lat_in = max(-89.5, min(89.5, lat) )  !adjust latitude 

!
!     CALCULATE EMIS AT GIVEN ANGLE
!

 if(bilinear) then
    call emis_interp_ind_sing_bilinear(lat_in,lon,angle,f, atlas, ev, eh,stdv,stdh,verb)
 else
    call emis_interp_ind_sing(lat_in,lon,angle,f, atlas, ev, eh,stdv,stdh,verb)
 endif

 emish_out = eh !horizontal polarization
 emisv_out = ev !vertical polarization 

!
!     CALCULATE EMIS AT VARIOUS ANGLES
!
 do i = 1,nang
    ang(i) = 4.*REAL( i - 1 )
    mu(i) = cos(ang(i)*const_pi/180.)
    angle = ang(i)

    if(bilinear) then
      call emis_interp_ind_sing_bilinear(lat_in,lon,angle,f, atlas, ev, eh,stdv,stdh,verb)
    else
      call emis_interp_ind_sing(lat_in,lon,angle,f, atlas, ev, eh,stdv,stdh,verb)
    endif

    esum(i) =  ev + eh
 enddo

!
!     CALCULATE EBAR
!
 tot = 0.0
 do i = 1,nang-1
    eavg = 0.5*( esum(i) + esum(i+1) )
    dmu = mu(i) - mu(i+1)
    avmu = 0.5*( mu(i) + mu(i+1) )
    tot = tot + eavg*avmu*dmu
 enddo
 ebar_out = tot
 return
 end subroutine  emit_telsem_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_surface_micro
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
! Write out microwave surface emissivity (V & H polarization) in GrADS format. 
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: i,j,l,m
 integer :: irec
 integer,parameter :: io = 201
 integer :: is,ie,js,je,ks,ke
 character(len=3) :: fff

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end

!
! Computed Tbs in binary format
!
  sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.EMIS.'//&
                 trim(micro_sensor)//trim(output_suffix)//'.bin'
#if MPI == 2
 if(masterproc) &
#endif
  print*,'MSG write_surface_micro; output for surface microwave emissivity',&
  ' ->',trim(sdsu_io_file)


#if MPI == 2
 if(masterproc) &
#endif
  open( unit = io, file = sdsu_io_file, access = 'direct', &
        status = 'replace', recl = mxgridx*mxgridy*4  )
 irec = 1
 do l = 1,mxfreq_micro
    do m = 1,2
       call dump(  emis_micro(is:ie,js:je,l,m) , io, irec )
    enddo
 enddo

#if MPI == 2
 if(masterproc) &
#endif
 close (io)

!
! write grads ctl file
!
#if MPI == 2
 if(masterproc) then
#endif

 if( write_grads_ctl  ) then
    sdsu_io_file =trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.EMIS.'//&
                  trim(micro_sensor)//trim(output_suffix)//'.ctl'
    print*,'see output.ctl ->',trim(sdsu_io_file)

    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.EMIS.'//trim(micro_sensor)//trim(output_suffix)//'.bin'

    write(io,*) 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(micro_sensor) //' microwave surface emissivity'
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
              write(io,*) 'emis'//fff//'h',' 0  0  emissivity '//trim(nch_micro(l))//'(H) [-]'
          elseif(m==2) then
              write(io,*) 'emis'//fff//'v',' 0  0  emissivity '//trim(nch_micro(l))//'(V) [-]'
          endif

       enddo
    enddo
    write(io,*) 'endvars '

    close(io)
 endif

#if MPI == 2
 endif
#endif

 return
 end subroutine write_surface_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine surface_radar
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!  
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

! if(masterproc) print*,'MSG surface_radar: Compute surface sigma zero'
!
! sensor angle (umu is global parameter)
!
 umu = cos(inc_angle_radar*const_pi/180.) !sensor incident angle


 return
 end subroutine surface_radar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine surface_visir
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   Compute surface albedo or emissivity at visible (monocromatic) wavelength.
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: ierr
 integer :: i,j,nw !looping indice


!
! initialize albedo LUT
!
 if( .not. allocated(albedo_visir) ) then

    allocate( &
             albedo_visir(1:mxgridx,1:mxgridy,1:mxwavel), &
             stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG surface_visir: allocation error -> Terminate program. ')

    !
    ! read MODIS albedo LUT
    !
    call read_modis_albedo_lut

 endif


!
! initialize emissivity 
!

 if( .not. allocated(emiss_visir) ) then

    allocate( &
             emiss_visir(1:mxgridx,1:mxgridy,1:mxwavel), &
             stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG surface_visir: allocation error -> Terminate program.')

    !
    ! input IR emissivity data file and directory
    !
    sdsu_io_file = trim(sdsu_dir_data)//'FILLED_IR_EMISS/DATABASE/'//trim(input_bf_emiss)

    call read_bf_emiss( trim(sdsu_io_file) , surface%lat, surface%lon  )

 endif


 if(masterproc) print*,'MSG surface_visir: Compute surface albedo or emissivity at visir spectrum'

 WAVE_LOOP: do nw = 1, mxwavel

    DOM_LOOP: do j = myj_start, myj_end ; do i = myi_start, myi_end

!
! compute spectrum land-surface albedo
!
       if( wavel(nw) < 3.3e0 ) then  !solar radiation spectrum

           if(surface(i,j)%iland == 2) then  !water

             albedo_visir(i,j,nw) = -999.  !will be computed in other subroutines (simulator_visir.F)

           elseif(surface(i,j)%iland == 1) then !land 

            call interp_modis_albedo_lut( wavel(nw), surface(i,j)%frac_veg, surface(i,j)%igbp_typ, &
                                        albedo_visir(i,j,nw)  )	
           endif
        else
          albedo_visir(i,j,nw) = undefined  !IR range
        endif

!
! put emissivity computation toshii (under construction)
!

        if( wavel(nw) >= 3.3e0 ) then  !IR radiation spectrum

            if(surface(i,j)%iland == 2) then  !water

               emiss_visir(i,j,nw) = -999.  !will be computed in other subroutines (simulator_visir.F)

            elseif(surface(i,j)%iland == 1) then !land 


               call interp_bf_emiss( wavel(nw), surface(i,j)%lat, surface(i,j)%lon, emiss_visir(i,j,nw)   )

!               emiss_visir(i,j,nw) = 0.99  !constant so far

            endif
        else
            emiss_visir(i,j,nw) = undefined  !negligible IR spectrum put unity just in case
        endif

   enddo ; enddo DOM_LOOP

 enddo WAVE_LOOP

!
! write out surface albedo in GrADS format
!
  if(write_surface) call write_surface_visir

 return
 end subroutine surface_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_modis_albedo_lut
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   Read MODIS-derived albedo LUT 
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer,parameter :: io_albedo=4
 character(len=200) :: file_albedo  ! file and directory
 character(len=100) :: dummy_character 
 integer :: id1, id2, id3, id4

 integer :: index_ndvi, index_lulc, index_wave !index

!--------------IGBP LULC type--------------
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

!
! open file
!
  file_albedo = trim(sdsu_dir_data)//'modis_albedo_lut.asc'
  open(io_albedo,file=file_albedo,status='old')

  read(io_albedo,*) dummy_character 

!
! read LUT as a function of LULC, wavelength, NDVI index
!
  do index_lulc = 1, imax_lulc
    do index_wave = 1, imax_wave
       do index_ndvi = 1, imax_ndvi
          read(io_albedo,*) id1, id2, id3, id4, modis_albedo_lut(index_lulc, index_wave, index_ndvi)
       enddo
    enddo
  enddo
  close(unit=io_albedo)


 return
 end subroutine read_modis_albedo_lut

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine vis_nir_albedo_lut
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!  Using MODIS-derived albedo LUT, this subroutine compute visible and NIR albedo.  
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer :: index_ndvi, index_lulc, index_wave !index
 real(sdsu_fps):: wave !wavelength [micron]

 real(sdsu_fpd) :: fplnk,fsol
 real(sdsu_fpd) :: bplnk,sunir            !! function of plank function
 real(sdsu_fps) :: wgt,sumwgt  !plank function weight and its sum
 real(sdsu_fps) :: sumalb !weighted sum of surface albedo
 real(sdsu_fps) :: alb_vis, alb_nir  !visible and NIR albedo
 

 do index_lulc = 1, imax_lulc
    do index_ndvi = 1, imax_ndvi

!
! visible albedo
!
        sumwgt = 0. ; sumalb = 0.
        do index_wave = 1, 3
           wave = modis_albedo_chan( index_wave ) !wavelength [micron]
           fsol=sunir( 1.d4/wave )
           wgt = REAL(fsol)
           sumwgt = sumwgt + wgt
           sumalb = sumalb + wgt * modis_albedo_lut(index_lulc, index_wave, index_ndvi)
        enddo
        alb_vis = sumalb / sumwgt

!
! NIR albedo
!
        sumwgt = 0. ; sumalb = 0.
        do index_wave = 4, imax_wave
           wave = modis_albedo_chan( index_wave ) !wavelength [micron]
           fsol=sunir( 1.d4/wave )
           wgt = REAL(fsol)
           sumwgt = sumwgt + wgt
           sumalb = sumalb + wgt * modis_albedo_lut(index_lulc, index_wave, index_ndvi)
        enddo
        alb_nir = sumalb / sumwgt

!
! compute vis-NIR albodo ratio
!
        modis_albedo_broad( index_lulc , index_ndvi , 1) = alb_vis
        modis_albedo_broad( index_lulc , index_ndvi , 2) = alb_nir
    enddo
 enddo

 return
 end subroutine vis_nir_albedo_lut

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine interp_modis_albedo_lut( wave, Fg, igbp, albedo )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
! Interpolate spectrum land-surface albedo from MODIS LUT. 
!
! Vegetation Fraction (Fg) is define here (Gutman and Ignatov 1998).
!
! Fg = [ NDVI - NDVI(Fg=0) ] / [ NDVI(Fg=1) - NDVI(Fg=0) ]
!
! Now we assume NDVI(Fg=1) = 1.0 & NDVI(Fg=0) = 0.0, then
!
! Fg = NDVI  (ranging 0~1)
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),intent(in) :: wave !wavelength [micron]
 real(sdsu_fps),intent(in) :: Fg  !vegetation fraction [%]
 integer,intent(in) :: igbp !IGBP land-cover type
 real(sdsu_fps) :: albedo !albedo type
 real(sdsu_fps) :: ndvi   !NDVI (0~1)
 real(sdsu_fps) :: wgt1_wav, wgt2_wav , wgt1_ndvi, wgt2_ndvi
 real(sdsu_fps) :: modis_albedo_wave(1:imax_wave)
 integer :: index_ndvi,iw

!
! Convert Fg to NDVI
! 
 ndvi = Fg * 1e-2
 index_ndvi = max( 1,  min( imax_ndvi , INT(ndvi*1e+3)/50 + 1 )) !NDVI index


!
! interpolate MODIS albedo LUT for NDVI index
!
 if(index_ndvi == 1 ) then
    
   modis_albedo_wave(1:imax_wave) = modis_albedo_lut(igbp, 1:imax_wave , 1) 
  
 elseif (index_ndvi == imax_ndvi) then

   modis_albedo_wave(1:imax_wave) = modis_albedo_lut(igbp, 1:imax_wave , imax_ndvi)
 
 else

   wgt1_ndvi = ( REAL(index_ndvi)*0.05 - ndvi ) / 0.05
   wgt2_ndvi = 1.e0 - wgt1_ndvi

   modis_albedo_wave(1:imax_wave) = wgt1_ndvi * modis_albedo_lut(igbp, 1:imax_wave , index_ndvi-1) &
                                   +wgt2_ndvi * modis_albedo_lut(igbp, 1:imax_wave , index_ndvi  )

 endif


!
! interpolate MODIS albeod LUT for wavelength index
!
 if     ( wave < modis_albedo_chan(1) ) then !less than MODIS blue channel

    albedo = modis_albedo_wave(1) 

 elseif ( wave > modis_albedo_chan(imax_wave) ) then !greater than MODIS MIR channel

    albedo = modis_albedo_wave(imax_wave)

 else !Linear interpolate  

    iw_loop: do iw = 1, imax_wave - 1
       if( wave >= modis_albedo_chan(iw) .and. wave <= modis_albedo_chan(iw+1) ) then
           wgt1_wav = ( modis_albedo_chan(iw+1) - wave ) / &
                      ( modis_albedo_chan(iw+1) - modis_albedo_chan(iw) )
           wgt2_wav = 1.- wgt1_wav

           albedo =   wgt1_wav * modis_albedo_wave( iw   )      &
                    + wgt2_wav * modis_albedo_wave( iw+1 )    

           exit iw_loop
       endif
    enddo iw_loop

 endif


 return
 end subroutine interp_modis_albedo_lut

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_surface_visir
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
! Write out visir surface albedo (white-sky) in GrADS format. 
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: i,j,nw
 integer :: irec
 integer,parameter :: io = 202
 integer :: is,ie,js,je,ks,ke
 character(len=2) :: ww

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end

!
! Computed Tbs in binary format
!
 sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.LAND_ALBEDO_EMISS.'//&
                trim(visir_sensor)//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG write_surface_visir; output for surface visible-IR albedo',&
 ' ->',trim(sdsu_io_file)

#if MPI == 2
 if(masterproc) &
#endif
 open( unit= io, file = sdsu_io_file, access = 'direct', &
       status = 'replace', recl = mxgridx*mxgridy*4  )

 irec=1
 do nw=1,mxwavel
    if( wavel(nw) < tb_or_rad ) then  !visible albedo
      call dump( albedo_visir(is:ie,js:je,nw) , io, irec )
    else !IR emissivity 
      call dump( emiss_visir(is:ie,js:je,nw) , io, irec )
    endif
 enddo

#if MPI == 2
 if(masterproc) &
#endif
 close (io)


#if MPI == 2
 if(masterproc) then
#endif

!
! write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.LAND_ALBEDO_EMISS.'//&
                   trim(visir_sensor)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.LAND_ALBEDO_EMISS.'//&
                   trim(visir_sensor)//trim(output_suffix)//'.bin'
    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(visir_sensor) //' albedo or emissivity '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',mxwavel

    do nw=1,mxwavel !frequency loop
       write(ww,"(I2.2)") nw
       if(wavel(nw) < tb_or_rad) then 
          write(io,*) 'alb'//ww//' 0  0  albedo at '//trim(nch_wavel(nw))//' [-]'
       else 
          write(io,*) 'emis'//ww//' 0  0  emissivity at '//trim(nch_wavel(nw))//' [-]'
       endif
    enddo
    write(io,*) 'endvars '
    close(io)
  endif

#if MPI == 2
 endif
#endif


 return
 end subroutine write_surface_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine surface_broad
 implicit none 
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   Compute surface albedo or emissivity at each broad band.
!   For SW albedo, it essentially computes VIS and NIR components from total broad-band albedo. 
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 integer :: ierr
 integer :: i,j,nw !looping indice
 real(sdsu_fps) :: alb_vis, alb_nir !albedo for visible and NIR band


 if(masterproc) print*,'MSG surface_broad: Compute surface albedo or emissivity at broadbands '

!
! initialize albedo
!
 if( .not. allocated(albedo_broad) ) then

    allocate( albedo_broad(1:mxgridx,1:mxgridy,1:2,1:2),stat=ierr ) !1-direct , 2-diffuse 
    if (ierr /= 0) call stop_sdsu(myrank,'MSG surface_broad: allocation error -> Terminate program. ')

   !
   ! read MODIS albedo LUT
   !
    call read_modis_albedo_lut

   !
   ! derive albedo ration betweeen visible and NIR spectrum
   !
    call vis_nir_albedo_lut  

 endif

!
! initialize emissivity
!
 if( .not. allocated(emiss_broad)) then

   allocate(emiss_broad(1:mxgridx,1:mxgridy),stat=ierr)  
   if (ierr /= 0) call stop_sdsu(myrank,'MSG surface_broad: allocation error -> Terminate program. ')

 endif

 DOM_LOOP: do j = myj_start, myj_end ; do i = myi_start, myi_end

     !
     ! compute spectrum albedo
     !
     if(surface(i,j)%iland == 2) then  !water

       !
       ! PUT CERES empirical ocean albedo later!!!!
       !
        albedo_broad(i,j,1:2,1:2) = surface(i,j)%albedo
       
     elseif(surface(i,j)%iland == 1) then !land  (simple vis-NIR separation)

        call modis_vis_nir_albedo( surface(i,j)%frac_veg, surface(i,j)%igbp_typ,&
                                   alb_vis, alb_nir )        

        !
        ! at present, separation happens for two broadband. (direct & diffuse values are same so far)
        !
        albedo_broad(i,j,1, 1:2) = alb_vis  !VISIBLE
        albedo_broad(i,j,2, 1:2) = alb_nir  !NIR
        
  
     endif

    !
    ! compute spectrum emissivity
    !

     emiss_broad(i,j) = 0.99

 enddo ; enddo DOM_LOOP

 return
 end subroutine surface_broad 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine modis_vis_nir_albedo( Fg, igbp, alb_vis, alb_nir )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!  Get broad-band albedo for  visible and NIR components, based MODIS LUT as
!  function of land-cover type and vegetation fraction 
!  THis is for land-surface only.
!
! Vegetation Fraction (Fg) is define here (Gutman and Ignatov 1998).
!
! Fg = [ NDVI - NDVI(Fg=0) ] / [ NDVI(Fg=1) - NDVI(Fg=0) ]
!
! Now we assume NDVI(Fg=1) = 1.0 & NDVI(Fg=0) = 0.0, then
!
! Fg = NDVI  (ranging 0~1)
!                                       
! History:
! 10/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 real(sdsu_fps),intent(in) :: Fg  !vegetation fraction [%]
 integer,intent(in) :: igbp !IGBP land-cover type
 real(sdsu_fps),intent(out) :: alb_vis, alb_nir  !albedo
 real(sdsu_fps) :: ndvi   !NDVI (0~1)
 real(sdsu_fps) :: wgt1_ndvi, wgt2_ndvi
 real(sdsu_fps) :: modis_albedo_wave(1:imax_wave)
 integer :: index_ndvi,iw
 real(sdsu_fps) :: ratio

!
! Convert Fg to NDVI
! 
 ndvi = Fg * 1e-2
 index_ndvi = max( 1,  min( imax_ndvi , INT(ndvi*1e+3)/50 + 1 )) !NDVI index

!
! interpolate MODIS albedo LUT for NDVI index
!
 if(index_ndvi == 1 ) then
    
   alb_vis = modis_albedo_broad( igbp , index_ndvi ,1)
   alb_nir = modis_albedo_broad( igbp , index_ndvi ,2) 
 
 elseif (index_ndvi == imax_ndvi) then

   alb_vis = modis_albedo_broad( igbp , index_ndvi ,1)
   alb_nir = modis_albedo_broad( igbp , index_ndvi ,2)

 else

   wgt1_ndvi = ( REAL(index_ndvi)*0.05 - ndvi ) / 0.05
   wgt2_ndvi = 1.e0 - wgt1_ndvi

   alb_vis =  wgt1_ndvi * modis_albedo_broad(igbp, index_ndvi-1, 1) &
             +wgt2_ndvi * modis_albedo_broad(igbp, index_ndvi  , 1)
   alb_nir =  wgt1_ndvi * modis_albedo_broad(igbp, index_ndvi-1, 2) &
             +wgt2_ndvi * modis_albedo_broad(igbp, index_ndvi  , 2)

 endif


 return
 end subroutine modis_vis_nir_albedo

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine surface_sarta
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   Compute surface emissivity at IR wavelength.
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: ierr
 integer :: i,j,nw !looping indice

!
! initialize emissivity 
!

 if( .not. allocated(emiss_sarta) ) then

    allocate( &
             emiss_sarta(1:mxgridx,1:mxgridy,1:anchor_bs_emiss), &
             stat=ierr )
    if (ierr /= 0) call stop_sdsu(myrank,'MSG surface_sarta: allocation error -> Terminate program.')

    !
    ! input IR emissivity data file and directory
    !
    sdsu_io_file = trim(sdsu_dir_data)//'FILLED_IR_EMISS/DATABASE/'//trim(input_bf_emiss)

    call read_bf_emiss( trim(sdsu_io_file) , surface%lat, surface%lon  )

 endif


 if(masterproc) print*,'MSG surface_sarta: Compute land-surface emissivity at IR spectrum'

 DOM_LOOP: do j = myj_start, myj_end ; do i = myi_start, myi_end

!
! put emissivity computation toshii (under construction)
!

    if(surface(i,j)%iland == 2) then  !water

       emiss_sarta(i,j,:) = undefined  

    elseif(surface(i,j)%iland == 1) then !land 

       call interp_bf_emiss_all( surface(i,j)%lat, surface(i,j)%lon, emiss_sarta(i,j,:)   )

    endif

 enddo ; enddo DOM_LOOP

!
! write out surface albedo in GrADS format
!
  if(write_surface) call write_surface_sarta

 return
 end subroutine surface_sarta

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_surface_sarta
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:                             
! Write out visir surface albedo (white-sky) in GrADS format. 
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 integer :: i,j,nw
 integer :: irec
 integer,parameter :: io = 204
 integer :: is,ie,js,je,ks,ke
 character(len=2) :: ww

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end ; ks=myk_start ; ke=myk_end

!
! Computed Tbs in binary format
!
 sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.LAND_SARTA_EMISS.'//&
                trim(visir_sensor)//trim(output_suffix)//'.bin'

#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG write_surface_visir; output for land-surface IR emissivity ',&
 ' ->',trim(sdsu_io_file)

#if MPI == 2
 if(masterproc) &
#endif
 open( unit= io, file = sdsu_io_file, access = 'direct', &
       status = 'replace', recl = mxgridx*mxgridy*4  )

 irec=1
 do nw=1,anchor_bs_emiss
    call dump( emiss_sarta(is:ie,js:je,nw) , io, irec )
 enddo

#if MPI == 2
 if(masterproc) &
#endif
 close (io)


#if MPI == 2
 if(masterproc) then
#endif

!
! write out grads control file
!
 if( write_grads_ctl  ) then
    sdsu_io_file = trim(sdsu_dir_output)//sdsu_inp_name(1:efile_len)//'.LAND_SARTA_EMISS.'//&
                   trim(visir_sensor)//trim(output_suffix)//'.ctl'
    print("(a)"),'see output.ctl ->',trim(sdsu_io_file)
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = sdsu_inp_name(1:efile_len)//'.LAND_SARTA_EMISS.'//&
                   trim(visir_sensor)//trim(output_suffix)//'.bin'
    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  '//trim(visir_sensor) //' baseline-fitted emissivity '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   ',mxgridx,' linear 1. ',gridsize
    write(io,*) 'ydef   ',mxgridy,' linear 1. ',gridsize
    write(io,*) 'zdef    1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ', anchor_bs_emiss

    do nw=1, anchor_bs_emiss !frequency loop
       write(ww,"(I2.2)") nw
       write(io,*) 'emis'//ww//' 0  0  emissivity at '//trim(wave_bs_emiss_cha(nw))//' [-]'
    enddo
    write(io,*) 'endvars '
    close(io)
  endif

#if MPI == 2
 endif
#endif


 return
 end subroutine write_surface_sarta

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 end module module_surface

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
