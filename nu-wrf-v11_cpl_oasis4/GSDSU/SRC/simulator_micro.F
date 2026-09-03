#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine simulator_micro
 use module_simulator
 use module_opt_micro
 use module_surface
 use module_mpi
 implicit none
!--------------------------------------------------------------------------------------------
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
!  Microwave brightness temperature is computed using an Eddington radiative
!  transfer code developed by Kummerow (1993). The geophysical parameters that
!  affect microwave brightness temperature are summarized in the following.
!
!  Surface emissivity is computed from 
!  surface type (ocean or land), surface tempterature in [K] (only for ocean),
!  near-surface wind in [m/s] (only for ocean), and soil moisture (0-1)
!  (only for land).
!
!  Output variables are horizontally and vertically polarized brightness
!  temperature in [K] at each frequency.
!
! History:
!   5/2017  Toshi Matsui@NASA GSFC ; output filter threshold has been changed for ground_micro (suggested by Mario Montopoli). 
!  09/2012  Toshi Matsui@NASA GSFC ; Incident angle are as a function of latitude. 
!  08/2008  Toshi Matsui@NASA GSFC ; Remove optical properties part into module_opt_micowave
!  07/2007  Toshi Matsui@NASA GSFC ; Add dynamic allocation and better-memory coding.  
!  05/2007  Toshi Matsui@NASA GSFC ; adding downwelling tb for output.
!  03/2007    Toshi Matsui@NASA GSFC ; Significant revision for f90  
!  2002~2006  Hirohiko Masunaga@CSU ; Adapted for SDSU
!  ~~~2000 Chris Kummerow & Willam Olson @ NASA GSFC  ; Initial version.  
!
! References:
!
!  Kummerow, C., 1993:  On the accuracy of the Eddington approximation for radiative transfer in 
!   the microwave frequencies. J. Geophys. Res., 98, 2757-2765.
!--------------------------------------------------------------------------------------------
 integer :: i,j,nz,nf,np,ng !looping indice
 integer :: ierr  !alloc stat
 integer :: k_bot , k_top , kmax !bottom  and top of k level
 real :: tb_down_surface
 real :: umu_local
 logical :: cycle_this
 real :: incidence !incidence angle [rad]
 real :: dummy

!
! Program Start
!
 if(masterproc) print*,'MSG simulator_micro: Compute microwave radiative transfer '

! frequency loop
  do nf = 1, mxfreq_micro

! polarization loop
    POL_LOOP: do np = 1,2

          if(verbose_SDSU .and. np == 0 .and. masterproc)  &
          print*, ' Frequency = ', freq_micro(nf),' GHz -> Wavelength = ',2.997925e2/freq_micro(nf) ,'mm' 

       
       ! skip radiative transfer for special channles.
       call skip_transfer(trim(micro_sensor), nf, np, cycle_this)
       if(cycle_this) then
           tb_out(:,:,nf,np) = undefined
            cycle POL_LOOP
       endif 

       do j = myj_start, myj_end ; do i = myi_start, myi_end

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

          if( trim(sim_case) == 'MMF_GPROF' .and. j > max_sample_gprof ) then
              tb_out(i,j,nf,np) = undefined 
              cycle
          endif

          !
          ! Apply Eddington calculation to vertical columns in
          ! model domain. -> now 1D driver
          !
          k_bot = 1 ; k_top = mxlyr  ; kmax = k_top - k_bot + 1
          tb_down_surface = 0.  

          call eddington_2way('downup', kmax , umu_local, ground_micro, surface(i,j)%t_skin, &
                              atmos_stag(i,j,k_bot-1:k_top)%t_air,   &
                              atmos_stag(i,j,k_bot-1:k_top)%hgt,     &
                              atmos_opt(i,j,k_bot:k_top,nf)%kext, &
                              atmos_opt(i,j,k_bot:k_top,nf)%salb, &
                              atmos_opt(i,j,k_bot:k_top,nf)%asym, &
                              emis_micro(i,j,nf,np), ebar_micro(i,j,nf),  &
                              tb_down_surface, tb_out(i,j,nf,np) )

          if( ground_micro ) then

           if ( (tb_out(i,j,nf,np) <   0.) .or. &
                (tb_out(i,j,nf,np) > 500.)        ) then
                 print*,'MSG simulator_micro: tb out of bounds i=',i,'j=',j,'nf=',nf,'np=',np &
                        ,'tb_out=',tb_out(i,j,nf,np),'[K]'
                 tb_out(i,j,nf,np) = undefined
           endif

          else

           if ( (tb_out(i,j,nf,np) <  50.) .or. &
                (tb_out(i,j,nf,np) > 350.)        ) then
                 print*,'MSG simulator_micro: tb out of bounds i=',i,'j=',j,'nf=',nf,'np=',np &
                        ,'tb_out=',tb_out(i,j,nf,np),'[K]'
                 tb_out(i,j,nf,np) = undefined
           endif

          endif


       enddo ; enddo

    enddo POL_LOOP  ! np

   enddo !nf


#if MPI == 2
if(numproc_tot>1) then
!
! let's send each tile into the master processor. 
!
   if(numproc_tot>1) call mpi_sdsu_communicate( 'TO_MASTER', tb_out )
endif
#endif



 return
 end subroutine simulator_micro 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine simulator_micro_prps
 use module_simulator
 use module_opt_micro
 use module_surface
 use module_mpi
 implicit none
!--------------------------------------------------------------------------------------------
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
!  Microwave brightness temperature is computed using an Eddington radiative
!  transfer code developed by Kummerow (1993). The geophysical parameters that
!  affect microwave brightness temperature are summarized in the following.
!
!  Surface emissivity is computed from 
!  surface type (ocean or land), surface tempterature in [K] (only for ocean),
!  near-surface wind in [m/s] (only for ocean), and soil moisture (0-1)
!  (only for land).
!
!  Output variables are horizontally and vertically polarized brightness
!  temperature in [K] at each frequency.
!
! History:
!   5/2019  Toshi Matsui@NASA GSFC ; modified simulator_micro for PRPS. 
!
! References:
!
!  Kummerow, C., 1993:  On the accuracy of the Eddington approximation for radiative transfer in 
!   the microwave frequencies. J. Geophys. Res., 98, 2757-2765.
!--------------------------------------------------------------------------------------------
 integer :: i,j,nz,nf,np,ng,n !looping indice
 integer :: ierr  !alloc stat
 integer :: k_bot , k_top , kmax !bottom  and top of k level
 real :: tb_down_surface
 real :: umu_local
 logical :: cycle_this
 real :: incidence !incidence angle [rad]
 real :: dummy

!
! Program Start
!
 if(masterproc) print*,'MSG simulator_micro_prps: Compute microwave radiative transfer '


! viewing angle loop
 do n = min_index_angle, max_index_angle

! frequency loop
  do nf = 1, mxfreq_micro

! polarization loop
    POL_LOOP: do np = 1,2

       if(verbose_SDSU .and. np == 0 .and. n == 0 .and. masterproc)  &
          print*, ' Frequency = ', freq_micro(nf),' GHz -> Wavelength = ',2.997925e2/freq_micro(nf) ,'mm' 

       !
       ! skip radiative transfer for special channles.
       !
       call skip_transfer(trim(micro_sensor), nf, np, cycle_this)
       if(cycle_this) then
           tb_prps_out(:,:,nf,np,:) = undefined
           cycle POL_LOOP
       endif 

       j_loop: do j = myj_start, myj_end 

          if( j > max_sample_gprof ) then
              tb_prps_out(:,j,nf,np,:) = undefined 
              cycle j_loop
          endif

          i_loop: do i = myi_start, myi_end

            !
            ! for faster calculation, calcualte RT within beamwidth
            !
            if( gf_prps(i,j,nf,n) == 0. ) then  !out of FOV
               tb_prps_out(i,j,nf,np,n) = undefined
               cycle i_loop 
            endif

            !
            ! derive cosine of incidence angle
            !
            umu_local = cos( inc_angle_prps(i,j,n)*const_degrad )  !cosine of sensor incidence angle

            !
            ! Apply Eddington calculation to vertical columns in
            ! model domain. -> now 1D driver
            !
            k_bot = 1 ; k_top = mxlyr  ; kmax = k_top - k_bot + 1
            tb_down_surface = 0.  

            call eddington_2way( 'downup', kmax , umu_local, ground_micro, surface(i,j)%t_skin, &
                              atmos_stag(i,j,k_bot-1:k_top)%t_air,   &
                              atmos_stag(i,j,k_bot-1:k_top)%hgt,     &
                              atmos_opt(i,j,k_bot:k_top,nf)%kext, &
                              atmos_opt(i,j,k_bot:k_top,nf)%salb, &
                              atmos_opt(i,j,k_bot:k_top,nf)%asym, &
                              emis_micro_prps(i,j,nf,np,n), ebar_micro_prps(i,j,nf,n),  &
                              tb_down_surface, tb_prps_out(i,j,nf,np,n) )

            if ( (tb_prps_out(i,j,nf,np,n) <  50.) .or. &
                (tb_prps_out(i,j,nf,np,n) > 350.)        ) then
                print*,'MSG simulator_micro_prps: tb out of bounds i=',i,'j=',j,'nf=',nf,'np=',np,'n=',n &
                      ,'tb_prps_out=',tb_prps_out(i,j,nf,np,n),'[K]'
                tb_prps_out(i,j,nf,np,n) = undefined
print*, kmax, umu_local, surface(i,j)%t_skin, emis_micro_prps(i,j,nf,np,n), ebar_micro_prps(i,j,nf,n)
stop

            endif

         enddo i_loop
       enddo j_loop

    enddo POL_LOOP  ! np

   enddo !nf

 enddo !n

!#if MPI == 2
!if(numproc_tot>1) then
!
! let's send each tile into the master processor. 
!
!   if(numproc_tot>1) call mpi_sdsu_communicate( 'TO_MASTER', tb_prps_out )
!endif
!#endif

 return
 end subroutine simulator_micro_prps 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine simulator_slant_micro
 use module_simulator
 use module_opt_micro
 use module_surface
 use module_mpi
 implicit none
!--------------------------------------------------------------------------------------------
!
! Comments: 
!  Microwave brightness temperature is computed using an Eddington radiative
!  transfer code through slant path, developed by Kummerow (1993) and Olson (2001). 
!  The geophysical parameters that
!  affect microwave brightness temperature are summarized in the following.
!
!  Surface emissivity is computed from 
!  surface type (ocean or land), surface tempterature in [K] (only for ocean),
!  near-surface wind in [m/s] (only for ocean), and soil moisture (0-1)
!  (only for land).
!
!  Output variables are horizontally and vertically polarized brightness
!  temperature in [K] at each frequency.
!
! History:
!   5/2017  Toshi Matsui@NASA GSFC ; output filter threshold has been changed for ground_micro (suggested by Mario Montopoli). 
!   3/2013  Toshi Matsui@NASA GSFC ; Bug fixed for k_top value. 
!  11/2013  Toshi Matsui@NASA GSFC ; Init.
!
! References:
!
!  Olson W.S., P. Bauer, C. D. Kummerow, Y. Hong, W.-K. Tao, 2001: A melting-layer model for 
!   passive/active microwave remote sensing applications. Part II: simulation of TRMM observations, 
!   J. Appl. Meteor, 40, 1164-1179.
!--------------------------------------------------------------------------------------------
 integer :: i,j,k, nz,nf,np,ng !looping indice
 integer :: iup, idown  
 integer :: ierr  !alloc stat
 integer :: k_bot , k_top , kmax !bottom  and top of k level
 real :: tb_down_surface
 real :: umu_local
 logical :: cycle_this
 real :: incidence !incidence angle [rad]
 real :: dummy
 type (opt_micro_parameter), allocatable :: slant_opt1d(:,:,:)  !slant-path optical properties. 
 real :: tair_slant1d( 1:mxlyr_refine , 1:2 )  !slant-path air temperature
 real :: t_air_stag  ( 0:mxlyr_refine , 1:2 )  !slant-path air temperature (staggered) 
 real(sdsu_fps),parameter :: tb_zero = 0. !dummy Tb

!
! Program Start
!
 if(masterproc) print*,'MSG simulator_slant_micro: Compute microwave radiative transfer '

 if( .not. allocated(slant_opt1d) ) then
    allocate( slant_opt1d(1:mxlyr_refine,1:mxfreq,1:2) , &
              stat=ierr)
    if (ierr /= 0) call stop_sdsu(myrank,'MSG simulator_slant_micro: allocation error -> Terminate program. ')
 endif

!
! Domain loop
!
 do j = myj_start, myj_end ; do i = myi_start, myi_end

    ! initialize output as undefined 
     tb_out(i,j,:,:) = undefined

    !special skip for MMF_GPROF cases
    if( trim(sim_case) == 'MMF_GPROF' .and. j > max_sample_gprof ) then
        cycle
    endif

    !
    ! derive cosine of incidence angle
    !
    call derive_incidence_angle( ona_angle_micro*const_degrad, sma_micro, surface(i,j)%lat*const_degrad, &
                                 incidence,  dummy )
    umu_local = cos(incidence)  !cosine of sensor incidence angle

    !
    ! serch slant-path column 
    !
    call search_slant_path_2d ( i, ibs, ibe, gridsize, incidence, hgt_stag_micro(0:mxlyr_micro) ,&
                                atmos_opt  ( ibs:ibe, j, 1:mxlyr_refine, 1:mxfreq )   , &
                                t_air_micro( ibs:ibe, j, 1:mxlyr_refine  )            , &
                                slant_opt1d( 1:mxlyr_refine,1:mxfreq,1:2 )   , &
                                tair_slant1d( 1:mxlyr_refine, 1:2), k_bot , k_top   )

    !
    ! skip strange column
    !
    if( k_bot == undefined_i2 ) then
      print*,'MSG simulator_slant_micro: k_bot is undefined at (i,j)',i,j 
      cycle
    endif

    !
    ! initialize memory bounds 
    !
    kmax = k_top - k_bot + 1

    !
    ! extract stagger air temperature
    !
    t_air_stag(k_bot-1,1:2) = tair_slant1d(k_bot,1:2)
    t_air_stag(k_top  ,1:2) = tair_slant1d(k_top,1:2)
    do k = k_bot , k_top-1
       t_air_stag(k,1:2) = 0.5 * ( tair_slant1d(k,1:2) + tair_slant1d(k+1,1:2) )
    enddo


    ! frequency loop
    do nf = 1, mxfreq_micro

       !
       ! Modify cosine of incidence angle for the case of GMI high frequency channels
       ! This small modification is only radiation path length.
       !
       if(trim(micro_sensor)== 'GMI' .and. freq_micro(nf) >= 100. ) then
          call derive_incidence_angle( 45.36*const_degrad, sma_micro, surface(i,j)%lat*const_degrad, &
                                       incidence,  dummy )
          umu_local = cos(incidence)  !cosine of sensor incidence angle
       endif

       ! polarization loop
       POL_LOOP: do np = 1,2

          ! skip radiative transfer for special channles.
          call skip_transfer(trim(micro_sensor), nf, np, cycle_this)
          if(cycle_this) then
              tb_out(i,j,nf,np) = undefined
              cycle POL_LOOP
          endif

          !
          ! Apply Eddington calculation to vertical columns in
          ! model domain. -> now 1D driver
          !
          tb_down_surface = 0.   !initialize

          idown = 2
          call eddington_2way('down', kmax , umu_local, ground_micro, surface(i,j)%t_skin  ,  &
                            t_air_stag(k_bot-1:k_top,idown), hgt_stag_micro(k_bot-1:k_top),  &
                            slant_opt1d (k_bot:k_top,nf,idown)%kext,  &
                            slant_opt1d (k_bot:k_top,nf,idown)%salb,  &
                            slant_opt1d (k_bot:k_top,nf,idown)%asym,  &
                            emis_micro(i,j,nf,np), ebar_micro(i,j,nf),      &
                            tb_zero, tb_down_surface   )
          iup = 1
          call eddington_2way('up', kmax , umu_local, ground_micro, surface(i,j)%t_skin  ,  &
                            t_air_stag(k_bot-1:k_top,iup), hgt_stag_micro(k_bot-1:k_top),  &
                            slant_opt1d (k_bot:k_top,nf,iup)%kext,  &
                            slant_opt1d (k_bot:k_top,nf,iup)%salb,  &
                            slant_opt1d (k_bot:k_top,nf,iup)%asym,  &
                            emis_micro(i,j,nf,np), ebar_micro(i,j,nf), &
                            tb_down_surface, tb_out(i,j,nf,np)    )

          !
          ! Filter bad Tb
          !
          if( ground_micro ) then

            if ( (tb_out(i,j,nf,np) <  0.) .or. &
                 (tb_out(i,j,nf,np) > 500.)        ) then
                 print*,'MSG simulator_slant_micro: tb out of bounds i=',i,'j=',j,'nf=',nf,'np=',np &
                       ,'tb_out=',tb_out(i,j,nf,np),'[K]'
                 tb_out(i,j,nf,np) = undefined
            endif


          else

            if ( (tb_out(i,j,nf,np) <  50.) .or. &
                 (tb_out(i,j,nf,np) > 350.)        ) then
                 print*,'MSG simulator_slant_micro: tb out of bounds i=',i,'j=',j,'nf=',nf,'np=',np &
                       ,'tb_out=',tb_out(i,j,nf,np),'[K]'
                 tb_out(i,j,nf,np) = undefined
            endif

          endif

       enddo POL_LOOP  ! np

    enddo !nf

 enddo ; enddo ! domain loop



#if MPI == 2
if(numproc_tot>1) then
!
! let's send each tile into the master processor. 
!
   if(numproc_tot>1) call mpi_sdsu_communicate( 'TO_MASTER', tb_out )
endif
#endif



 return
 end subroutine simulator_slant_micro 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine search_slant_path_2d ( i_bot, ibs, ibe, gridsize, incidence, hgt_stag_micro, &
                                   atmos_opt2d, tair_2d, slant_opt1d , tair_slant1d , k_bot, k_top )
 use module_opt_micro
 use module_simulator,only : sdsu_fps, stop_sdsu, mxlyr_micro
 use module_const
 implicit none
!---------------------------------------------------------------------------------------------------
! Note treating 2D i-k domain as cartesian cordinate. 
! * Downwelling radiation path is always from west. 
! * Upwelling radiaiton path is toward east. 
! * There is no interpolation, just sample nearest grid. 
!
!  down     upward
!    \      /
!     \    /
!      \  /
!       \/
!===================
!       i_bot 
!
! History:
!  05/2014  Eric Kemp  : Bug fixed (removed mxlyr_refine and mxfreq from module_simulator)
!   3/2013  Toshi Matsui@NASA GSFC ; Bug fixed for k_top value. 
!  11/2013  Toshi Matsui@NASA GSFC ; Init.
!
! References:
!---------------------------------------------------------------------------------------------------

 integer,intent(in) :: i_bot               ! i potisiotin of bottom layer
 integer,intent(in) :: ibs, ibe            ! i direction starting and ending index
 real,intent(in) :: gridsize  !horizontal grid spacing [km]
 real,intent(in) :: incidence !incidence angle [rad]
 real,intent(in) :: hgt_stag_micro(mxlyr_micro)
 type(opt_micro_parameter),intent(in)  :: atmos_opt2d ( ibs:ibe, 1:mxlyr_refine, 1:mxfreq )
 real,intent(in) :: tair_2d( ibs:ibe, 1:mxlyr_refine )
 type(opt_micro_parameter),intent(out) :: slant_opt1d( mxlyr_refine,mxfreq ,2 )
 real,intent(out)    :: tair_slant1d( 1:mxlyr_refine , 1:2 )
 integer,intent(out) :: k_bot !bottomo index
 integer,intent(out) :: k_top !top index
 integer :: i, k
 integer :: i_slante, i_slantw
 real :: agl_km  ! above ground level [km]

!
! initialize as undefined
!
 slant_opt1d%kext = undefined
 slant_opt1d%salb = undefined
 slant_opt1d%asym = undefined
 tair_slant1d     = undefined


!
! determine bottom layer considering terrain
!
 k_bot = undefined_i2
 do k = 1, mxlyr_refine
    if( atmos_opt2d ( i_bot, k, 1 )%kext /= undefined ) then !find bottom layer (just above terrain)
       k_bot = k
       exit 
    endif 
 enddo

 if(k_bot == undefined_i2) then
    return
 endif

!
! re-sample optical properties and air temperature in slant path
!
 do k = k_bot, mxlyr_refine
   
    !
    ! find i position based on incidence angle
    !
    agl_km  = hgt_stag_micro(k) - hgt_stag_micro(k_bot)  ! above ground level [km]

    i_slantw = max(ibs, i_bot - nint( agl_km * incidence/gridsize ) )  ! slant i index (from west )
    i_slante = min(ibe, i_bot + nint( agl_km * incidence/gridsize ) )  ! slant i index (to east )

    slant_opt1d( k ,1:mxfreq ,1) = atmos_opt2d ( i_slantw, k , 1:mxfreq )  ! upwelling portion (from west)
    slant_opt1d( k ,1:mxfreq ,2) = atmos_opt2d ( i_slante, k , 1:mxfreq )  ! downwelling portion (to east)

    tair_slant1d( k , 1 )  =  tair_2d( i_slantw, k ) ! upwelling portion (from west)
    tair_slant1d( k , 2 )  =  tair_2d( i_slante, k ) ! downwelling portion (to east)

 enddo


!
! find out top layer height
!
 k_top = undefined_i2
 find_ktop: do k = mxlyr_refine, k_bot, -1 
    if( tair_slant1d( k , 1 ) /= undefined .and. tair_slant1d( k , 2 ) /= undefined ) then
       k_top = k
       exit
    endif
 enddo find_ktop


 return
 end subroutine search_slant_path_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine skip_transfer(sensor_name, nf, np, skip_this)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!          Skip certain channels for faster processing. 
!
! History:
!  09/2012 Toshi Matsui    @NASA GSFC : Init
!           
! References: 
!---------------------------------------------------------------------------------------------------
 character(*),intent(in) :: sensor_name
 integer,intent(in) :: nf
 integer,intent(in) :: np
 logical,intent(out) :: skip_this

!
! initialize as false
!
 skip_this = .false.

 select case(trim(sensor_name))
 case('TMI')
    if(nf==3 .and. np==1) skip_this = .true. !no 24G horizontal channel

 case('GMI')
    if(nf==3 .and. np==1) skip_this = .true. !no 24G horizontal channel
    if(nf>=7 .and. np==1) skip_this = .true. !no H channel in 183GHz

 case('SSMIS')
    if(nf==2 .and. np==1) skip_this = .true.  !no 22.235G horizontal channel (V only)
    if( (nf>=4 .and. nf<=8) .and. np==2 )  skip_this = .true. !no O3-band vertical channel (H only)
    if(nf>=10 .and. np==2) skip_this = .true.  !no V in sounding channles
 case('MADRAS')
    if(nf==2 .and. np==1) skip_this = .true. !no 24G horizontal channel
 case default
    !do nothing
 end select

 return
 end subroutine skip_transfer

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

!     Plane parallel Eddington computes Tb for one vertical hydrometeor
!     column at a time.  

 subroutine eddington_2way(path, kmax, umu, ground_based,  tskin, lyrtemp, lyrhgt, kexttot, salbtot, asymtot, &
                          emis, ebar, tb_down_surface, tb_out )
 use module_mpi
 use module_simulator, only : sdsu_fps, stop_sdsu
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
!   Compute downwelling and upwelling microwave brightness temperature (Tb) on the assumption of 
!   specular surface throgh analytical Eddington solution.  This routine can be used eith one-way
!   or two-way plane-paralell slant-path assumptions. 
!   
!   A. one-way assumption. 
!     1.  Call eddingto_2way by setting path = "downup'. (tb_down_surface won't be used.)
!
!   B. two-way assumption. 
!     1. Call eddington_2way, by setting path = "down" and use atmosphrofile from downwelling slant path, 
!      and get tb_out, which is downwelling surface Tb. 
!     2.  Second, call eddington_2way again, but set path = "up" using atmospheric profiles from upwelling 
!        slant path. also use tb_out as input for tb_down_surface.  
!
! History:
!  09/2011 Toshi Matsui    @NASA GSFC : Modifcaiton for f90 + specular 2way assumption 
!  02/1995 Bill Olson      @NASA GSFC : Include asymptotic expression (term(3:5) for small tau)
!     1993 Chris Kummerow  @NASA GSFC : Initial
!           
! References: 
! Kummerow, C., 1993:  On the accuracy of the Eddington approximation for radiative transfer in 
!                      the microwave frequencies. J. Geophys. Res., 98, 2757-2765.
! Weinman, J., and R. Davies, 1978: Thermal Microwave Radiances From Horizontally Finite Clouds 
!                      of Hydrometeors, J. Geophys. Res., 83(C6), 3099-3107.
!---------------------------------------------------------------------------------------------------
 character*(*),intent(in) :: path            ! up, down, updown
 integer,intent(in) :: kmax                  ! array bounds
 logical,intent(in) :: ground_based          ! ground- or satellite-based 
 real(sdsu_fps),intent(in) :: umu            ! cosine of incidence angle [-]
 real(sdsu_fps),intent(in) :: lyrtemp(0:kmax) ! air temperature [K]
 real(sdsu_fps),intent(in) :: lyrhgt (0:kmax) ! layer height [km]
 real(sdsu_fps),intent(in) :: kexttot(kmax) ! total extinction coef [km-1]
 real(sdsu_fps),intent(in) :: salbtot(kmax) ! total single scatterling albedo [-]
 real(sdsu_fps),intent(in) :: asymtot(kmax) ! total asymetery parameters [-]
 real(sdsu_fps),intent(in) :: tskin                  ! surface skin temperature [K]
 real(sdsu_fps),intent(in) :: emis                   ! surface emissivity [-]
 real(sdsu_fps),intent(in) :: ebar                   ! 
  real(sdsu_fps),intent(in) ::tb_down_surface        ! downwelling surface Tb from other path [K}
 real(sdsu_fps),intent(out) :: tb_out                ! output microwave Tb [K]

!
! local parameters
!
 integer :: i,j,k ! loop
 real(sdsu_fps) :: btemp !bottom temp [K]
 real(sdsu_fps) :: dz  ! layer thickness [km]
 real(sdsu_fps) :: term(5)  !various terms
 real(sdsu_fps) :: mu, nu  !mu parameters
 real(sdsu_fps) :: xnu,xa,xb,xc,xd,ya,yb,yc,rcond !mix parameter
 real(sdsu_fps) :: w (2*kmax,2*kmax)
 real(sdsu_fps) :: bb(2*kmax)
 real(sdsu_fps) :: z (0:kmax)
 real(sdsu_fps) :: tb_up(0:kmax)       !upwelling microwave Tb [K]
 real(sdsu_fps) :: tb_down(1:kmax+1)   !downwelling microwave Tb [K]
 real(sdsu_fps),parameter :: tb_cosmic= 2.725e0  ! Cosmic microwave background radiation [K]

 real(sdsu_fps),dimension(kmax) :: &
   kext       ,&  ! extinction coef [km-1]
   salb       ,&  ! single scatterling albedo [-]
   asym       ,&  ! asymetery parameters [-]
   l          ,&  ! 
   h          ,&  !
   b0         ,&  ! air temperature at stagger height [K]
   b1         ,&  ! dT/dZ [K/km]
   dp         ,&  !
   dm             !


!
! Pass input parameters to local parameters 
!
 z(0) = lyrhgt(0) ! height at surface [km]
 btemp    = tskin           ! bottom temperature [K]

 do k = 1, kmax
    z(k)  = lyrhgt(k)     !layer height (stagger) [km]
    b0(k) = lyrtemp(k-1)  !air temperature at stagger height [K]
    b1(k) = (lyrtemp(k) - lyrtemp(k-1))/ &
            (z(k) - z(k-1))  ! dT/dZ  [K/km]
    kext(k) = kexttot(k)  !extinction [km-1]
    salb(k) = salbtot(k)  !single scattering albedo [-]
    asym(k) = asymtot(k)  !asymetry parameter
    l(k)    = sqrt(3.*kext(k)*kext(k)*(1. - salb(k))* &
              (1. - salb(k)*asym(k)))         ! lambda in eq (6))  through alaysitcal Eddington approximation.
    h(k)    = 1.5*kext(k)*(1. - salb(k)*asym(k))  ! extinction except forward scatter
 enddo 

!
!  FILL IN THE MATRIX ELEMENTS WHICH FORM THE BOUNDARY CONDITIONS
!  AT THE TOP, BOTTOM AND LAYER INTERFACES OF THE CLOUD.  THERE ARE
!  TWO QUANTITIES, D+ "DP", AND D- "DM" AT EACH BOUNDARY, SO THE
!  MATRIX HAS DIMENSIONS  2*NLYR BY 2*NLYR.
!  ORDER OF D'S:  D+(1),D-(1),D+(2),D-(2), .... , D+(NLYR),D-(NLYR)
!

 w(:,:) = 0.0e0 ! Set all matrix elements to zero. 

!
! Fill in the non-zero matrix elements
!
 w(1,1)   = ((ebar - 2.)*l(1)/h(1)) + ebar
 w(1,2)   = ((2. - ebar)*l(1)/h(1)) + ebar

 do i = 2,2*(kmax-1),2
    w(i,i-1)   =  (1. - l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2)- z(i/2-1)))
    w(i,i  )   =  (1. + l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2)- z(i/2-1)))
    w(i,i+1)   = -(1. - l(i/2+1)/h(i/2+1))
    w(i,i+2)   = -(1. + l(i/2+1)/h(i/2+1))
    w(i+1,i-1) =  (1. + l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2)- z(i/2-1)))
    w(i+1,i)   =  (1. - l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2)- z(i/2-1)))
    w(i+1,i+1) = -(1. + l(i/2+1)/h(i/2+1))
    w(i+1,i+2) = -(1. - l(i/2+1)/h(i/2+1))
 enddo !i

 w(2*kmax,2*kmax-1) =  (1. + l(kmax)/h(kmax))*exp(+l(kmax)* &
                          (z(kmax)-z(kmax-1)))
 w(2*kmax,2*kmax)   =  (1. - l(kmax)/h(kmax))*exp(-l(kmax)* &
                          (z(kmax)-z(kmax-1)))

!
! Fill in the row of constants in the linear equations
!
 bb(1)    = ebar *btemp - ebar*b0(1) -  (ebar - 2.)*b1(1)/h(1)
 do i = 2,2*(kmax-1),2
    bb(i)   =  + b1(i/2)/h(i/2) - b1(i/2+1)/h(i/2+1)
    bb(i+1) =  - b1(i/2)/h(i/2) + b1(i/2+1)/h(i/2+1)
 enddo 

 bb(2*kmax)  =  tb_cosmic - b0(kmax) - b1(kmax)*(z(kmax) - &
                 z(kmax-1) + 1/h(kmax))

!
! Matrix inversion 
!
 call linpak2( kmax, w, bb, rcond)

 do i = 1,kmax
    dp(i) = bb(2*i-1)
    dm(i) = bb(2*i)
 enddo

!
! constant 
!
 mu = umu
 nu = -umu
 tb_down(kmax+1) = tb_cosmic  !cosmic background micorwave Tb [K]


if( path == 'down' .or. path == 'downup' ) then
!
! Compute downwelling radiance from top of layer k  at angle mu only (specular-surface assumption)
!
 DOWN_LOOP: do k = kmax, 1 ,-1  !loop from top to bottom

    xa = b0(k) - 1.5*salb(k)*asym(k)*nu*b1(k)/h(k)
    xb = b1(k)
    xc = salb(k)*dp(k)*(1. - 1.5*asym(k)*nu*l(k)/h(k))
    xd = salb(k)*dm(k)*(1. + 1.5*asym(k)*nu*l(k)/h(k))
    ya = kext(k)/nu
    yb = ya + l(k)
    yc = ya - l(k)
    dz = z(k) - z(k-1)

    term(1) = tb_down(k+1)*exp(ya*dz)

    term(2) = xa*(1. - exp(ya*dz))

    if(abs(ya*dz) .lt. 1.e-5) then
       term(3) = -xb*ya*dz*dz
    else
       term(3) = xb/ya*(exp(ya*dz)*(1. - ya*dz) - 1.)
    end if

    if(abs(yb*dz) .lt. 1.e-5) then
       term(4) = -xc*ya*dz
    else
       term(4) = xc*ya/yb*(1. - exp(yb*dz))
    end if

    if(abs(yc*dz) .lt. 1.e-5) then
       term(5) = -xd*ya*dz
    else
       term(5) = xd*ya/yc*(1. - exp(yc*dz))
    end if

    tb_down(k) = sum( term(1:5) ) ! downwelling Tb [K]

 enddo DOWN_LOOP

endif  !down or downup

if( path == 'downup' ) then
   !
   ! bottom boundary using the same path
   !
   tb_up(0) = emis*btemp + (1. - emis)*tb_down(1)  !bottom of upwellling Tb [K]

elseif( path == 'up') then
   !
   ! bottom bound using different path
   !
   tb_up(0) = emis*btemp + (1. - emis)*tb_down_surface
elseif( path == 'down') then
   !
   ! output tb and return
   !
   tb_out = tb_down(1)
   return
else
 call stop_sdsu(myrank,'MSG: eddington_2way: no such path')
endif



!
! Compute upwelling radiance from top of layer k  at angle mu only (specular-surface assumption)
!
 UP_LOOP: do k = 1,kmax

    xa = b0(k) - 1.5*salb(k)*asym(k)*mu*b1(k)/h(k)
    xb = b1(k)
    xc = salb(k)*dp(k)*(1. - 1.5*asym(k)*mu*l(k)/h(k))
    xd = salb(k)*dm(k)*(1. + 1.5*asym(k)*mu*l(k)/h(k))
    ya = kext(k)/mu
    yb = ya + l(k)
    yc = ya - l(k)
    dz = z(k) - z(k-1)

    term(1) = tb_up(k-1)*exp(-ya*dz)

    term(2) = xa*(1. - exp(-ya*dz))

    if(abs(ya*dz) .lt. 1.e-5) then
       term(3) = 0.
    else
       term(3) = xb/ya*(ya*dz - 1. + exp(-ya*dz))
    end if

    if(abs(yb*dz) .lt. 1.e-5) then
        term(4) = xc*ya*dz*exp(-ya*dz)
    else
        term(4) = xc*ya/yb*(exp( (yb-ya)*dz ) - exp(-ya*dz) )
    end if

    if(abs(yc*dz) .lt. 1.e-5) then
        term(5) = xd*ya*dz*exp(-ya*dz)
    else
        term(5) = xd*ya/yc*exp(-ya*dz)*(exp(yc*dz) - 1.)
    end if

    tb_up(k) = sum(term(1:5) ) ! upwelling Tb [K]

 enddo UP_LOOP


!
! Output
!
 if(ground_based) then
    !
    ! If ground-based microwave radiometer, output is downwellling 
    ! microwave Tb at surface level
    !
    tb_out = tb_down(1)
 else
    !
    ! If satellite-based microwave radiometer, output is upwelling
    ! microwave Tb at TOA
    !
    tb_out = tb_up(kmax)
 endif
    

 return
 end subroutine eddington_2way

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine linpak2(kmax, a, b, rcond)
 use module_simulator, only : sdsu_fps
 implicit none
!
! TEST OUT THE MATRIX INVERSION PROCEDURE
!
 integer,intent(in) :: kmax
 real(sdsu_fps) :: a(2*kmax, 2*kmax), b(2*kmax), z(2*kmax), rcond
 integer lda, n, ipvt(2*kmax)

 lda = 2*kmax
 n   = 2*kmax

 call sgeco(a,lda,n,ipvt,rcond,z)
 call sgesl(a,lda,n,ipvt,b    ,0)

 return
 end subroutine linpak2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!     Plane parallel Eddington computes Tb for one vertical hydrometeor
!     column at a time.  

      subroutine eddington_pp(tskin, lyrtemp, lyrhgt, kexttot, salbtot, asymtot, &
                              emis, ebar, lambert, tb_out )
      use module_simulator, only : umu, mxlyr,ground_micro,sdsu_fps
!
!     CHRIS KUMMEROW
!     INCLUDES ASYMPTOTIC EXPRESSIONS FOR TERM3, TERM4, AND TERM5 IN
!     THE LIMIT OF SMALL EFFECTIVE OPTICAL DEPTHS; BILL OLSON FEB, 1995.
!
      implicit none
!
! IO parameters
!
      logical,intent(in) :: lambert
      real(sdsu_fps),intent(in) :: lyrtemp(0:mxlyr)
      real(sdsu_fps),intent(in) :: lyrhgt(0:mxlyr)
      real(sdsu_fps),intent(in) :: tskin
      real(sdsu_fps),intent(in) :: emis
      real(sdsu_fps),intent(in) :: ebar
      real(sdsu_fps),intent(in) :: kexttot(mxlyr)
      real(sdsu_fps),intent(in) :: salbtot(mxlyr)
      real(sdsu_fps),intent(in) :: asymtot(mxlyr)

      real(sdsu_fps) :: tb     !upwelling TOA Tb
      real(sdsu_fps) :: tb_bot !downwelling surface Tb 
      real(sdsu_fps) :: tb_out !output
!
! local parameters
!
      integer :: i,j,nn ! loop
      integer :: nang   ! stream number

      real(sdsu_fps) :: btemp !bottom temp [K]
      real(sdsu_fps) :: xnu,xa,xb,xc,xd,ya,yb,yc,rcond
      real(sdsu_fps) :: dz  !thickness
      real(sdsu_fps) :: term1, term2, term3, term4, term5
      real(sdsu_fps) :: xiup  !total donwelling radiation

      real(sdsu_fps) :: kext(mxlyr), salb(mxlyr), asym(mxlyr)
      real(sdsu_fps) :: l(mxlyr), h(mxlyr), b0(mxlyr), b1(mxlyr)
      real(sdsu_fps) :: w(2*mxlyr,2*mxlyr), bb(2*mxlyr), dp(mxlyr), dm(mxlyr)
      real(sdsu_fps) :: z(0:mxlyr), iout(0:mxlyr), i_in(mxlyr+1,100), mu, nu
!
      real(sdsu_fps),parameter :: fisot= 2.725e0 !cosmic background Tb [K]

!
! Pass input parameters to local parameters 
!
          z(0) = lyrhgt(0)
          btemp = tskin
          do j = 1,mxlyr
            z(j)  = lyrhgt(j)
            b0(j) = lyrtemp(j-1)
            b1(j) = (lyrtemp(j) - lyrtemp(j-1))/ &
                    (z(j) - z(j-1))
            kext(j) = kexttot(j)
            salb(j) = salbtot(j)
            asym(j) = asymtot(j)
            l(j) = sqrt(3.*kext(j)*kext(j)*(1. - salb(j))* &
                  (1. - salb(j)*asym(j)))
            h(j) = 1.5*kext(j)*(1. - salb(j)*asym(j))
          enddo !j

!
!         FILL IN THE MATRIX ELEMENTS WHICH FORM THE BOUNDARY CONDITIONS
!         AT THE TOP, BOTTOM AND LAYER INTERFACES OF THE CLOUD.  THERE ARE
!         TWO QUANTITIES, D+ "DP", AND D- "DM" AT EACH BOUNDARY, SO THE
!         MATRIX HAS DIMENSIONS  2*NLYR BY 2*NLYR.
!         ORDER OF D'S:  D+(1),D-(1),D+(2),D-(2), .... , D+(NLYR),D-(NLYR)

!         SET ALL MATRIX ELEMENTS TO ZERO
          do i = 1,2*mxlyr
            do j = 1,2*mxlyr
              w(i,j) = 0.0
            enddo !j
          enddo !i

!         FILL IN THE NON-ZERO MATRIX ELEMENTS
          w(1,1)   = ((ebar - 2.)*l(1)/h(1)) + ebar
          w(1,2)   = ((2. - ebar)*l(1)/h(1)) + ebar
          do i = 2,2*(mxlyr-1),2
            w(i,i-1)   =  (1. - l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i,i  )   =  (1. + l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i,i+1)   = -(1. - l(i/2+1)/h(i/2+1))
            w(i,i+2)   = -(1. + l(i/2+1)/h(i/2+1))

            w(i+1,i-1) =  (1. + l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i+1,i)   =  (1. - l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i+1,i+1) = -(1. + l(i/2+1)/h(i/2+1))
            w(i+1,i+2) = -(1. - l(i/2+1)/h(i/2+1))
          enddo !i

          w(2*mxlyr,2*mxlyr-1) =  (1. + l(mxlyr)/h(mxlyr))*exp(+l(mxlyr)* &
                                (z(mxlyr)-z(mxlyr-1)))
          w(2*mxlyr,2*mxlyr)   =  (1. - l(mxlyr)/h(mxlyr))*exp(-l(mxlyr)* &
                                (z(mxlyr)-z(mxlyr-1)))

!         FILL IN THE ROW OF CONSTANTS IN THE LINEAR EQUATIONS
          bb(1)    = ebar *btemp - ebar*b0(1) - &
                     (ebar - 2.)*b1(1)/h(1)
          do i = 2,2*(mxlyr-1),2
            bb(i)   =  + b1(i/2)/h(i/2) - b1(i/2+1)/h(i/2+1)
            bb(i+1) =  - b1(i/2)/h(i/2) + b1(i/2+1)/h(i/2+1)
          enddo 

          bb(2*mxlyr)  =  fisot - b0(mxlyr) - b1(mxlyr)*(z(mxlyr) - &
                         z(mxlyr-1) + 1/h(mxlyr))

!         MATRIX INVERSION IN DONE IN SUBROUTINE LINPAK
          call linpak( w, bb, rcond)

!
          do i = 1,mxlyr
            dp(i) = bb(2*i-1)
            dm(i) = bb(2*i)
          enddo

!         AFTER D'S ARE KNOWN, CALCULATE SURFACE RADIANCE

          mu = umu
          nu = -umu

!         FOR THE FOLLOWING CALCULATIONS, REFER TO APPENDIX B OF THESIS
!         *************************************************************
          if ( lambert ) then   !Lambertian (land background)
!           CALCULATE THE DOWNWELLING FLUX AT 81 ANGLES
            nang = 81
            do nn = 1,nang
              xnu = -(2.*nn - 1.)/(nang*2.)
              i_in(mxlyr+1,nn) = fisot
!             LOOP THROUGH THE REMAINING LAYERS
              do j = mxlyr,1,-1
!               CALCULATE RADIANCE FROM TOP OF LAYER "J"
                xa = b0(j) - 1.5*salb(j)*asym(j)*xnu*b1(j)/h(j)
                xb = b1(j)
                xc = salb(j)*dp(j)*(1. - 1.5*asym(j)*xnu*l(j)/h(j))

                xd = salb(j)*dm(j)*(1. + 1.5*asym(j)*xnu*l(j)/h(j))
                ya = kext(j)/xnu
                yb = ya + l(j)
                yc = ya - l(j)
                dz = z(j) - z(j-1)

                term1 = i_in(j+1,nn)*exp(ya*dz)
                term2 = xa*(1. - exp(ya*dz))
                if(abs(ya*dz) .lt. 1.e-5) then
                 term3=-xb*ya*dz*dz
                else
                 term3 = xb/ya*(exp(ya*dz)*(1. - ya*dz) - 1.)
                endif
                if(abs(yb*dz) .lt. 1.e-5) then
                 term4=-xc*ya*dz
                else
                 term4 = xc*ya/yb*(1. - exp(yb*dz))
                end if
                if(abs(yc*dz) .lt. 1.e-5) then
                 term5=-xd*ya*dz
                else
                 term5 = xd*ya/yc*(1. - exp(yc*dz))
                end if
                i_in(j,nn) = term1 + term2 + term3 + term4 + term5
              enddo ! j loop
            enddo ! nn loop
!
!           CALCULATE THE TOTAL DOWNWELLING FLUX REACHING THE SURFACE
            xiup = 0.
            do nn = 1,nang
              xiup = xiup + i_in(1,nn)*(1./nang)*(2.*nn-1.)/(2.*nang)
            enddo !nn
            xiup = 2.*xiup

          else  !non lambert (water backgroud)

!           CALCULATE THE DOWNWELLING FLUX AT ANGLE MU ONLY
            nn = 22            ! THIS IS A DUMMY INDEX FOR I_IN
            i_in(mxlyr+1,nn) = fisot
!           LOOP THROUGH THE REMAINING LAYERS
            do j = mxlyr,1,-1
!             CALCULATE RADIANCE FROM TOP OF LAYER "J"
              xa = b0(j) - 1.5*salb(j)*asym(j)*nu*b1(j)/h(j)
              xb = b1(j)
              xc = salb(j)*dp(j)*(1. - 1.5*asym(j)*nu*l(j)/h(j))
              xd = salb(j)*dm(j)*(1. + 1.5*asym(j)*nu*l(j)/h(j))
              ya = kext(j)/nu
              yb = ya + l(j)
              yc = ya - l(j)
              dz = z(j) - z(j-1)

              term1 = i_in(j+1,nn)*exp(ya*dz)
              term2 = xa*(1. - exp(ya*dz))
              if(abs(ya*dz) .lt. 1.e-5) then
               term3=-xb*ya*dz*dz
              else
               term3 = xb/ya*(exp(ya*dz)*(1. - ya*dz) - 1.)
              end if
              if(abs(yb*dz) .lt. 1.e-5) then
               term4=-xc*ya*dz
              else
               term4 = xc*ya/yb*(1. - exp(yb*dz))
              end if
              if(abs(yc*dz) .lt. 1.e-5) then
               term5=-xd*ya*dz
              else
               term5 = xd*ya/yc*(1. - exp(yc*dz))
              end if
              i_in(j,nn) = term1 + term2 + term3 + term4 + term5
            enddo !j loop

            xiup = i_in(1,22)
!
          endif  !lambert (land or ocean background)

! bottom boundary
          iout(0) = emis*btemp + (1. - emis)*xiup
          tb_bot = xiup


!           CALCULATE THE UPWELLING RADIANCES AT THE TOP OF EACH LAYER J
          do j = 1,mxlyr
            xa = b0(j) - 1.5*salb(j)*asym(j)*mu*b1(j)/h(j)
            xb = b1(j)
            xc = salb(j)*dp(j)*(1. - 1.5*asym(j)*mu*l(j)/h(j))
            xd = salb(j)*dm(j)*(1. + 1.5*asym(j)*mu*l(j)/h(j))
            ya = kext(j)/mu
            yb = ya + l(j)
            yc = ya - l(j)
            dz = z(j) - z(j-1)
            term1 = iout(j-1)*exp(-ya*dz)
            term2 = xa*(1. - exp(-ya*dz))
            if(abs(ya*dz) .lt. 1.e-5) then
             term3=0.
            else
             term3 = xb/ya*(ya*dz - 1. + exp(-ya*dz))
            end if
            if(abs(yb*dz) .lt. 1.e-5) then
             term4=xc*ya*dz*exp(-ya*dz)
            else
             term4 = xc*ya/yb*(exp( (yb-ya)*dz ) - exp(-ya*dz) )
            end if
            if(abs(yc*dz) .lt. 1.e-5) then
             term5=xd*ya*dz*exp(-ya*dz)
            else
             term5 = xd*ya/yc*exp(-ya*dz)*(exp(yc*dz) - 1.)
            end if
            iout(j) = term1 + term2 + term3 + term4 + term5
          enddo !j loop
!
          tb = iout(mxlyr)

      
      if(ground_micro) then
         tb_out = tb_bot
      else
         tb_out = tb
      endif
    

      return
      end subroutine eddington_pp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine eddington_sp( tskin, lyrtemp, lyrhgt, kexttot, salbtot, asymtot, &
                               emis, ebar, lambert, tb_out )
      use module_simulator, only :umu, gridsize, mxgridx, mxgridy, mxlyr,ground_micro,sdsu_fps

!
!     CHRIS KUMMEROW
!     INCLUDES ASYMPTOTIC EXPRESSIONS FOR TERM3, TERM4, AND TERM5 IN
!     THE LIMIT OF SMALL EFFECTIVE OPTICAL DEPTHS; BILL OLSON FEB, 1995.
!
      logical lambert(mxgridx,mxgridy)
      real(sdsu_fps) ::  lyrtemp(mxgridx,mxgridy,0:mxlyr)
      real(sdsu_fps) ::  lyrhgt(0:mxlyr) !sp uses 1D hgt for simplicity (Toshi: temporal )
      real(sdsu_fps) ::  tskin(mxgridx,mxgridy)
      real(sdsu_fps) ::  emis(mxgridx,mxgridy)
      real(sdsu_fps) ::  ebar(mxgridx,mxgridy)
      real(sdsu_fps) ::  kexttot(mxgridx,mxgridy,mxlyr)
      real(sdsu_fps) ::  salbtot(mxgridx,mxgridy,mxlyr)
      real(sdsu_fps) ::  asymtot(mxgridx,mxgridy,mxlyr)
      real(sdsu_fps) ::  kext(mxlyr), &
            salb(mxlyr), asym(mxlyr), l(mxlyr), h(mxlyr), &
            b0(mxlyr), b1(mxlyr), w(2*mxlyr,2*mxlyr), &
            bb(2*mxlyr), dp(mxlyr), dm(mxlyr), z(0:mxlyr)
      real(sdsu_fps) :: iout(0:mxlyr), i_in(mxlyr+1,100), mu, nu
      real(sdsu_fps) ::  tb(mxgridx,mxgridy)
      real(sdsu_fps) ::  tb_bot(mxgridx,mxgridy) !downwelling TB at the bottome
      real(sdsu_fps),intent(out) :: tb_out(mxgridx,mxgridy)


      integer nymid

      real(sdsu_fps) :: zmid
      integer *4 jin(0:mxlyr)
      integer *4 jout(0:mxlyr)

!
      data fisot / 2.725e0 /
!
!     Slant-Path Eddington computes Tb for one oblique radiance
!     path at a time.  The coordinate i,j is the gridpoint where
!     the radiance path intersects the earth's surface.

!     Calculate the mid-layer coordinates of the slant path for
!     both incoming and outgoing paths.


      do nz=0,mxlyr

        ! Find altitude of midlayer (or surface)
        if(nz .eq. 0) then
          zmid = lyrhgt(nz)
        else
          zmid = (lyrhgt(nz-1)+lyrhgt(nz))/2.
        end if

        ! Find nearest gridpoint to slant path intersection of midlayer.

        jin(nz) = nint(zmid*sqrt(1.-umu*umu)/(umu*gridsize))
        jout(nz) = nint(-zmid*sqrt(1.-umu*umu)/(umu*gridsize))

      enddo



!     Loop over gridpoints.

      do nx = 1, mxgridx
        do ny = 1, mxgridy


    ! Construct 1-D "column" along incoming ray path

          z(0) = lyrhgt(0)
          btemp = tskin(nx,ny)

          do j = 1,mxlyr

    ! Index of path intersection of midlayer; grid is assumed
    ! periodic in y.
            nymid = mod(2*mxgridy + ny + jin(j) - 1, mxgridy) + 1

            z(j)  = lyrhgt(j)
            b0(j) = lyrtemp(nx,nymid,j-1)
            b1(j) = (lyrtemp(nx,nymid,j) - &
                    lyrtemp(nx,nymid,j-1))/ &
                    (z(j) - z(j-1))
            kext(j) = kexttot(nx,nymid,j)
            salb(j) = salbtot(nx,nymid,j)
            asym(j) = asymtot(nx,nymid,j)
            l(j) = sqrt(3.*kext(j)*kext(j)*(1. - salb(j))* &
                  (1. - salb(j)*asym(j)))
            h(j) = 1.5*kext(j)*(1. - salb(j)*asym(j))

          enddo ! j loop
!
!         FILL IN THE MATRIX ELEMENTS WHICH FORM THE BOUNDARY CONDITIONS
!         AT THE TOP, BOTTOM AND LAYER INTERFACES OF THE CLOUD.  THERE ARE
!         TWO QUANTITIES, D+ "DP", AND D- "DM" AT EACH BOUNDARY, SO THE
!         MATRIX HAS DIMENSIONS  2*NLYR BY 2*NLYR.
!         ORDER OF D'S:  D+(1),D-(1),D+(2),D-(2), .... , D+(NLYR),D-(NLYR)

!         SET ALL MATRIZ ELEMENTS TO ZERO
          do i = 1,2*mxlyr
            do j = 1,2*mxlyr
              w(i,j) = 0.0
            enddo !j
          enddo !i

!         FILL IN THE NON-ZERO MATRIX ELEMENTS
          w(1,1)   = ((ebar(nx,ny) - 2.)*l(1)/h(1)) + ebar(nx,ny)
          w(1,2)   = ((2. - ebar(nx,ny))*l(1)/h(1)) + ebar(nx,ny)
          do i = 2,2*(mxlyr-1),2
            w(i,i-1)   =  (1. - l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i,i  )   =  (1. + l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i,i+1)   = -(1. - l(i/2+1)/h(i/2+1))
            w(i,i+2)   = -(1. + l(i/2+1)/h(i/2+1))

            w(i+1,i-1) =  (1. + l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i+1,i)   =  (1. - l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i+1,i+1) = -(1. + l(i/2+1)/h(i/2+1))
            w(i+1,i+2) = -(1. - l(i/2+1)/h(i/2+1))
          enddo !i
 
          w(2*mxlyr,2*mxlyr-1) =  (1. + l(mxlyr)/h(mxlyr))*exp(+l(mxlyr)* &
                                (z(mxlyr)-z(mxlyr-1)))
          w(2*mxlyr,2*mxlyr)   =  (1. - l(mxlyr)/h(mxlyr))*exp(-l(mxlyr)* &
                                (z(mxlyr)-z(mxlyr-1)))

!         FILL IN THE ROW OF CONSTANTS IN THE LINEAR EQUATIONS
          bb(1)    = ebar(nx,ny)*btemp - ebar(nx,ny)*b0(1) - &
                     (ebar(nx,ny) - 2.)*b1(1)/h(1)
          do i = 2,2*(mxlyr-1),2
            bb(i)   =  + b1(i/2)/h(i/2) - b1(i/2+1)/h(i/2+1)
            bb(i+1) =  - b1(i/2)/h(i/2) + b1(i/2+1)/h(i/2+1)
          enddo !i loop
 
          bb(2*mxlyr)  =  fisot - b0(mxlyr) - b1(mxlyr)*(z(mxlyr) - &
                         z(mxlyr-1) + 1/h(mxlyr))
!
!         MATRIX INVERSION IN DONE IN SUBROUTINE LINPAK
          call linpak( w, bb, rcond)
!
          do i = 1,mxlyr
            dp(i) = bb(2*i-1)
            dm(i) = bb(2*i)
          enddo !i 

!         AFTER D'S ARE KNOWN, CALCULATE SURFACE RADIANCE

          mu = umu
          nu = -umu

!         FOR THE FOLLOWING CALCULATIONS, REFER TO APPENDIX B OF THESIS
!         *************************************************************

          if ( lambert(nx,ny) ) then
!           CALCULATE THE DOWNWELLING FLUX AT 81 ANGLES
            nang = 81
            do nn = 1,nang
              xnu = -(2.*nn - 1.)/(nang*2.)
              i_in(mxlyr+1,nn) = fisot
!             LOOP THROUGH THE REMAINING LAYERS
              do j = mxlyr,1,-1
!               CALCULATE RADIANCE FROM TOP OF LAYER "J"
                xa = b0(j) - 1.5*salb(j)*asym(j)*xnu*b1(j)/h(j)
                xb = b1(j)
                xc = salb(j)*dp(j)*(1. - 1.5*asym(j)*xnu*l(j)/h(j))
                xd = salb(j)*dm(j)*(1. + 1.5*asym(j)*xnu*l(j)/h(j))
                ya = kext(j)/xnu
                yb = ya + l(j)
                yc = ya - l(j)
                dz = z(j) - z(j-1)

                term1 = i_in(j+1,nn)*exp(ya*dz)
                term2 = xa*(1. - exp(ya*dz))
                if(abs(ya*dz) .lt. 1.e-5) then
                 term3=-xb*ya*dz*dz
                else
                 term3 = xb/ya*(exp(ya*dz)*(1. - ya*dz) - 1.)
                endif
                if(abs(yb*dz) .lt. 1.e-5) then
                 term4=-xc*ya*dz
                else
                 term4 = xc*ya/yb*(1. - exp(yb*dz))
                end if
                if(abs(yc*dz) .lt. 1.e-5) then
                 term5=-xd*ya*dz
                else
                 term5 = xd*ya/yc*(1. - exp(yc*dz))
                end if
                i_in(j,nn) = term1 + term2 + term3 + term4 + term5

              enddo !j loop
           enddo ! nn 
!
!           CALCULATE THE TOTAL DOWNWELLING FLUX REACHING THE SURFACE
            xiup = 0.
            do nn = 1,nang
              xiup = xiup + i_in(1,nn)*(1./nang)*(2.*nn-1.)/(2.*nang)
            enddo !nn loop
 
            xiup = 2.*xiup

          else

!           CALCULATE THE DOWNWELLING FLUX AT ANGLE MU ONLY
            nn = 22            ! THIS IS A DUMMY INDEX FOR I_IN
            i_in(mxlyr+1,nn) = fisot
!           LOOP THROUGH THE REMAINING LAYERS
            do j = mxlyr,1,-1
!             CALCULATE RADIANCE FROM TOP OF LAYER "J"
              xa = b0(j) - 1.5*salb(j)*asym(j)*nu*b1(j)/h(j)
              xb = b1(j)
              xc = salb(j)*dp(j)*(1. - 1.5*asym(j)*nu*l(j)/h(j))
              xd = salb(j)*dm(j)*(1. + 1.5*asym(j)*nu*l(j)/h(j))
              ya = kext(j)/nu
              yb = ya + l(j)
              yc = ya - l(j)
              dz = z(j) - z(j-1)

              term1 = i_in(j+1,nn)*exp(ya*dz)
              term2 = xa*(1. - exp(ya*dz))
              if(abs(ya*dz) .lt. 1.e-5) then
               term3=-xb*ya*dz*dz
              else
               term3 = xb/ya*(exp(ya*dz)*(1. - ya*dz) - 1.)
              end if
              if(abs(yb*dz) .lt. 1.e-5) then
               term4=-xc*ya*dz
              else
               term4 = xc*ya/yb*(1. - exp(yb*dz))
              end if
              if(abs(yc*dz) .lt. 1.e-5) then
               term5=-xd*ya*dz
              else
               term5 = xd*ya/yc*(1. - exp(yc*dz))
              end if
              i_in(j,nn) = term1 + term2 + term3 + term4 + term5
            enddo !j loop 
            xiup = i_in(1,22)
!
          endif

!  bottom boundary
          iout(0) = emis(nx,ny)*btemp + (1. - emis(nx,ny))*xiup
          tb_bot(nx,ny) = xiup

!    Recalculate flux conditions along outgoing path.

    ! Construct 1-D "column" along outgoing ray path

          z(0) = lyrhgt(0)
          btemp = tskin(nx,ny)

          do j = 1,mxlyr

    ! Index of path intersection of midlayer; grid is assumed
    ! periodic in y.
            nymid = mod(2*mxgridy + ny + jout(j) - 1, mxgridy) + 1

            z(j)  = lyrhgt(j)
            b0(j) = lyrtemp(nx,nymid,j-1)
            b1(j) = (lyrtemp(nx,nymid,j) - lyrtemp(nx,nymid,j-1))/ &
                    (z(j) - z(j-1))
            kext(j) = kexttot(nx,nymid,j)
            salb(j) = salbtot(nx,nymid,j)
            asym(j) = asymtot(nx,nymid,j)
            l(j) = sqrt(3.*kext(j)*kext(j)*(1. - salb(j))* &
                  (1. - salb(j)*asym(j)))
            h(j) = 1.5*kext(j)*(1. - salb(j)*asym(j))

          enddo !j loop 
!
!         FILL IN THE MATRIX ELEMENTS WHICH FORM THE BOUNDARY CONDITIONS
!         AT THE TOP, BOTTOM AND LAYER INTERFACES OF THE CLOUD.  THERE ARE
!         TWO QUANTITIES, D+ "DP", AND D- "DM" AT EACH BOUNDARY, SO THE
!         MATRIX HAS DIMENSIONS  2*NLYR BY 2*NLYR.
!         ORDER OF D'S:  D+(1),D-(1),D+(2),D-(2), .... , D+(NLYR),D-(NLYR)

!         SET ALL MATRIZ ELEMENTS TO ZERO
          do i = 1,2*mxlyr
            do j = 1,2*mxlyr
              w(i,j) = 0.0
            enddo !j
          enddo !i

!         FILL IN THE NON-ZERO MATRIX ELEMENTS
          w(1,1)   = ((ebar(nx,ny) - 2.)*l(1)/h(1)) + ebar(nx,ny)
          w(1,2)   = ((2. - ebar(nx,ny))*l(1)/h(1)) + ebar(nx,ny)
          do i = 2,2*(mxlyr-1),2
            w(i,i-1)   =  (1. - l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i,i  )   =  (1. + l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i,i+1)   = -(1. - l(i/2+1)/h(i/2+1))
            w(i,i+2)   = -(1. + l(i/2+1)/h(i/2+1))

            w(i+1,i-1) =  (1. + l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i+1,i)   =  (1. - l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2)- &
                          z(i/2-1)))
            w(i+1,i+1) = -(1. + l(i/2+1)/h(i/2+1))
            w(i+1,i+2) = -(1. - l(i/2+1)/h(i/2+1))
          enddo !i loop 
          w(2*mxlyr,2*mxlyr-1) =  (1. + l(mxlyr)/h(mxlyr))*exp(+l(mxlyr)* &
                                (z(mxlyr)-z(mxlyr-1)))
          w(2*mxlyr,2*mxlyr)   =  (1. - l(mxlyr)/h(mxlyr))*exp(-l(mxlyr)* &
                                (z(mxlyr)-z(mxlyr-1)))

!         FILL IN THE ROW OF CONSTANTS IN THE LINEAR EQUATIONS
          bb(1)    = ebar(nx,ny)*btemp - ebar(nx,ny)*b0(1) - &
                     (ebar(nx,ny) - 2.)*b1(1)/h(1)
          do i = 2,2*(mxlyr-1),2
            bb(i)   =  + b1(i/2)/h(i/2) - b1(i/2+1)/h(i/2+1)
            bb(i+1) =  - b1(i/2)/h(i/2) + b1(i/2+1)/h(i/2+1)
          enddo !i loop
 
          bb(2*mxlyr)  =  fisot - b0(mxlyr) - b1(mxlyr)*(z(mxlyr) - &
                         z(mxlyr-1) + 1/h(mxlyr))
!
!         MATRIX INVERSION IN DONE IN SUBROUTINE LINPAK
          call linpak( w, bb, rcond)
!
          do i = 1,mxlyr
            dp(i) = bb(2*i-1)
            dm(i) = bb(2*i)
          enddo !i loop 

          do j = 1,mxlyr
!           CALCULATE THE UPWELLING RADIANCES AT THE TOP OF EACH LAYER J
            xa = b0(j) - 1.5*salb(j)*asym(j)*mu*b1(j)/h(j)
            xb = b1(j)
            xc = salb(j)*dp(j)*(1. - 1.5*asym(j)*mu*l(j)/h(j))
            xd = salb(j)*dm(j)*(1. + 1.5*asym(j)*mu*l(j)/h(j))
            ya = kext(j)/mu
            yb = ya + l(j)
            yc = ya - l(j)
            dz = z(j) - z(j-1)
            term1 = iout(j-1)*exp(-ya*dz)
            term2 = xa*(1. - exp(-ya*dz))
            if(abs(ya*dz) .lt. 1.e-5) then
             term3=0.
            else
             term3 = xb/ya*(ya*dz - 1. + exp(-ya*dz))
            end if
            if(abs(yb*dz) .lt. 1.e-5) then
             term4=xc*ya*dz*exp(-ya*dz)
            else
             term4 = xc*ya/yb*(exp( (yb-ya)*dz ) - exp(-ya*dz) )
            end if
            if(abs(yc*dz) .lt. 1.e-5) then
             term5=xd*ya*dz*exp(-ya*dz)
            else
             term5 = xd*ya/yc*exp(-ya*dz)*(exp(yc*dz) - 1.)
            end if
            iout(j) = term1 + term2 + term3 + term4 + term5
          enddo !j loop 
!
          tb(nx,ny) = iout(mxlyr)

        end do !nx
      end do !ny
      !$OMP END PARALLEL DO


      if(ground_micro) then
         tb_out = tb_bot
      else
         tb_out = tb 
      endif
 

      return
      end subroutine eddington_sp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine linpak(a, b, rcond)
!
!     TEST OUT THE MATRIX INVERSION PROCEDURE
      use module_simulator, only : mxlyr,sdsu_fps


!      PARAMETER ( MXLYR = 41 ) !28

      real(sdsu_fps) :: a(2*mxlyr, 2*mxlyr), b(2*mxlyr), z(2*mxlyr)
      integer lda, n, ipvt(2*mxlyr)

      lda = 2*mxlyr
      n   = 2*mxlyr

      call sgeco(a,lda,n,ipvt,rcond,z)
      call sgesl(a,lda,n,ipvt,b,0)
      return
      end subroutine linpak

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine sgeco(a,lda,n,ipvt,rcond,z)
      use module_simulator, only : sdsu_fps
!***BEGIN PROLOGUE  SGECO
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  CONDITION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION AND ESTIMATES
!            THE CONDITION NUMBER OF THE MATRIX.
!***DESCRIPTION
!
!     SGECO FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
!     AND ESTIMATES THE CONDITION OF THE MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW SGECO BY SGESL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW SGECO BY SGEDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW SGECO BY SGEDI.
!
!     ON ENTRY
!
!        A       REAL(LDA, N)
!                THE MATRIX TO BE FACTORED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        RCOND   REAL
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!
!        Z       REAL(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK SGEFA
!     BLAS SAXPY,SDOT,SSCAL,SASUM
!     FORTRAN ABS,AMAX1,SIGN
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  SASUM,SAXPY,SDOT,SGEFA,SSCAL
!***END PROLOGUE  SGECO
      integer lda,n,ipvt(1)
      real(sdsu_fps) :: a(lda,1),z(1)
      real(sdsu_fps) :: rcond
!
      real(sdsu_fps) :: sdot,ek,t,wk,wkm
      real(sdsu_fps) ::  anorm,s,sasum,sm,ynorm
      integer info,j,k,kb,kp1,l
!
!     COMPUTE 1-NORM OF A
!
!***FIRST EXECUTABLE STATEMENT  SGECO
      anorm = 0.0e0
      do 10 j = 1, n
         anorm = amax1(anorm,sasum(n,a(1,j),1))
   10 continue
!
!     FACTOR
!
      call sgefa(a,lda,n,ipvt,info)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.
!
!     SOLVE TRANS(U)*W = E
!
      ek = 1.0e0
      do 20 j = 1, n
         z(j) = 0.0e0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))
         if (abs(ek-z(k)) .le. abs(a(k,k))) go to 30
            s = abs(a(k,k))/abs(ek-z(k))
            call sscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = abs(wk)
         sm = abs(wkm)
         if (a(k,k) .eq. 0.0e0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0e0
            wkm = 1.0e0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + abs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + abs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
!
!     SOLVE TRANS(L)*Y = W
!
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + sdot(n-k,a(k+1,k),1,z(k+1),1)
         if (abs(z(k)) .le. 1.0e0) go to 110
            s = 1.0e0/abs(z(k))
            call sscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
!
      ynorm = 1.0e0
!
!     SOLVE L*V = Y
!
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call saxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (abs(z(k)) .le. 1.0e0) go to 130
            s = 1.0e0/abs(z(k))
            call sscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
!
!     SOLVE  U*Z = V
!
      do 160 kb = 1, n
         k = n + 1 - kb
         if (abs(z(k)) .le. abs(a(k,k))) go to 150
            s = abs(a(k,k))/abs(z(k))
            call sscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0e0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0e0) z(k) = 1.0e0
         t = -z(k)
         call saxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
!     MAKE ZNORM = 1.0
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
!
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end subroutine sgeco

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine sgesl(a,lda,n,ipvt,b,job)
      use module_simulator, only : sdsu_fps
!***BEGIN PROLOGUE  SGESL
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  SOLVES THE REAL SYSTEM A*X=B OR TRANS(A)*X=B
!            USING THE FACTORS OF SGECO OR SGEFA
!***DESCRIPTION
!
!     SGESL SOLVES THE REAL SYSTEM
!     A * X = B  OR  TRANS(A) * X = B
!     USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
!
!     ON ENTRY
!
!        A       REAL(LDA, N)
!                THE OUTPUT FROM SGECO OR SGEFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SGECO OR SGEFA.
!
!        B       REAL(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
!                            TRANS(A)  IS THE TRANSPOSE.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
!        OR SGEFA HAS SET INFO .EQ. 0 .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS SAXPY,SDOT
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  SAXPY,SDOT
!***END PROLOGUE  SGESL
      integer lda,n,ipvt(1),job
      real(sdsu_fps) :: a(lda,1),b(1)
!
      real(sdsu_fps) :: sdot,t
      integer k,kb,l,nm1
!***FIRST EXECUTABLE STATEMENT  SGESL
      nm1 = n - 1
      if (job .ne. 0) go to 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B
!
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
!
!        NOW SOLVE  U*X = Y
!
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call saxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
         do 60 k = 1, n
            t = sdot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
!
!        NOW SOLVE TRANS(L)*X = Y
!
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end subroutine sgesl

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      real function sasum(n,sx,incx)
      use module_simulator, only : sdsu_fps
!
!***PURPOSE  SUM OF MAGNITUDES OF S.P VECTOR COMPONENTS
!
!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
!
!     --OUTPUT--
!    SASUM  SINGLE PRECISION RESULT (ZERO IF N .LE. 0)
!
!     RETURNS SUM OF MAGNITUDES OF SINGLE PRECISION SX.
!     SASUM = SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
!
      real(sdsu_fps) :: sx(1)
!
      sasum = 0.0e0
      if(n.le.0)return
      if(incx.eq.1)goto 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      ns = n*incx
          do 10 i=1,ns,incx
          sasum = sasum + abs(sx(i))
   10     continue
      return
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
!
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sasum = sasum + abs(sx(i))
   30 continue
      if( n .lt. 6 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        sasum = sasum + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2)) &
        + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
   50 continue
      return
      end function sasum 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine saxpy(n,sa,sx,incx,sy,incy)
      use module_simulator, only : sdsu_fps
!
!***PURPOSE  S.P. COMPUTATION Y = A*X + Y
!
!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SA  SINGLE PRECISION SCALAR MULTIPLIER
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
!       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY
!
!     --OUTPUT--
!       SY  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)
!
!     OVERWRITE SINGLE PRECISION SY WITH SINGLE PRECISION SA*SX +SY.
!     FOR I = 0 TO N-1, REPLACE  SY(LY+I*INCY) WITH SA*SX(LX+I*INCX) +
!       SY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N
!       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
!
      real(sdsu_fps) :: sx(1),sy(1),sa
!
      if(n.le.0.or.sa.eq.0.e0) return

      if(incx.eq.incy) if(incx-1) 5,20,60
      !if(incx.eq.incy .and. incx-1 < 0) then
    5 continue
!
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return

!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          sy(i) = sa*sx(i) + sy(i)
   70     continue
      return



      end subroutine saxpy

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      real function sdot(n,sx,incx,sy,incy)
      use module_simulator, only : sdsu_fps
!
!***PURPOSE  S.P. INNER PRODUCT OF S.P. VECTORS
!
!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
!       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY
!
!     --OUTPUT--
!     SDOT  SINGLE PRECISION DOT PRODUCT (ZERO IF N .LE. 0)
!
!     RETURNS THE DOT PRODUCT OF SINGLE PRECISION SX AND SY.
!     SDOT = SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
      real(sdsu_fps) :: sx(1),sy(1)
!
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1)5,20,60
    5 continue
!
!        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sdot = sdot + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sdot = sdot + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sdot = sdot + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) + &
         sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
      return
!
!        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
!
   60 continue
      ns=n*incx
      do 70 i=1,ns,incx
        sdot = sdot + sx(i)*sy(i)
   70   continue
      return
      end function sdot

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine sscal(n,sa,sx,incx)
      use module_simulator, only : sdsu_fps
!
!***PURPOSE  S.P. VECTOR SCALE X = A*X
!
!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SA  SINGLE PRECISION SCALE FACTOR
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
!
!     --OUTPUT--
!       SX  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)
!
!     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX.
!     FOR I = 0 TO N-1, REPLACE SX(1+I*INCX) WITH  SA * SX(1+I*INCX)
!
      real(sdsu_fps) :: sa,sx(1)
!
      if(n.le.0)return
      if(incx.eq.1)goto 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      ns = n*incx
          do 10 i = 1,ns,incx
          sx(i) = sa*sx(i)
   10     continue
      return
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end subroutine sscal

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      integer function isamax(n,sx,incx)
      use module_simulator, only : sdsu_fps
!
!***PURPOSE  FIND LARGEST COMPONENT OF S.P. VECTOR
!
!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
!
!     --OUTPUT--
!   ISAMAX  SMALLEST INDEX (ZERO IF N .LE. 0)
!
!     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF SINGLE PRECISION SX.
!     ISAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(SX(1-INCX+I*INCX)
!
      real(sdsu_fps) :: sx(1),smax,xmag
!
      isamax = 0
      if(n.le.0) return
      isamax = 1
      if(n.le.1)return
      if(incx.eq.1)goto 20
!
!        CODE FOR INCREMENTS NOT EQUAL TO 1.
!
      smax = abs(sx(1))
      ns = n*incx
      ii = 1
          do 10 i=1,ns,incx
          xmag = abs(sx(i))
          if(xmag.le.smax) go to 5
          isamax = ii
          smax = xmag
    5     ii = ii + 1
   10     continue
      return
!
!        CODE FOR INCREMENTS EQUAL TO 1.
!
   20 smax = abs(sx(1))
      do 30 i = 2,n
         xmag = abs(sx(i))
         if(xmag.le.smax) go to 30
         isamax = i
         smax = xmag
   30 continue
      return
      end function isamax

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

      subroutine sgefa(a,lda,n,ipvt,info)
      use module_simulator, only : sdsu_fps
!***BEGIN PROLOGUE  SGEFA
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
!***DESCRIPTION
!
!     SGEFA FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
!
!     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
!
!     ON ENTRY
!
!        A       REAL(LDA, N)
!                THE MATRIX TO BE FACTORED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
!                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS SAXPY,SSCAL,ISAMAX
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  ISAMAX,SAXPY,SSCAL
!***END PROLOGUE  SGEFA
      integer lda,n,ipvt(1),info
      real(sdsu_fps) :: a(lda,n)
!
      real(sdsu_fps) :: t
      integer isamax,j,k,kp1,l,nm1
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
!***FIRST EXECUTABLE STATEMENT  SGEFA
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!
!        FIND L = PIVOT INDEX
!
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
         if (a(l,k) .eq. 0.0e0) go to 40
!
!           INTERCHANGE IF NECESSARY
!
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
!
!           COMPUTE MULTIPLIERS
!
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end subroutine sgefa

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
