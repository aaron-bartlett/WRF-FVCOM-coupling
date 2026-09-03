#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_radar
 use module_simulator
 use module_opt_micro
 use module_surface
 use module_mpi

 implicit none
!--------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!                             Radar Module
!
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
!
! Comments: 
!    Radar module that compute radar echo at various profiles. 
!
! History:
!  05/2011  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
! Matsui, T., J. Santanello, J. J. Shi, W.-K. Tao, D. Wu, C. Peters-Lidard, E. Kemp, M. Chin, 
!   D. Starr, M. Sekiguchi, and F. Aires, (2014): Introducing multisensor satellite radiance-based 
!   evaluation for regional Earth System modeling, Journal of Geophysical Research, 119, 8450-8475, 
!   doi:10.1002/2013JD021424.
!-----------------------------------------------------------------------------------------------------

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 private   ! ALL variables and subourtines are not accessible in module_radar

!            except...
 public :: simulator_radar, radar_return_scan, radar_return


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine simulator_radar(within_swath)
 implicit none
!--------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!
! Comments: 
!  The methodology to simulate radar echoes and path-integrated attenuation
!  is described in the Appendix of Masunaga and Kummerow (2005). The geophysical
!  parameters that affect radar simulation are summarized in the following.
!
!
! History:
!  06/2017  Toshi Matsui@NASA GSFC : PIA revised for both space and ground based.
!  03/2013  Toshi Matsui@NASA GSFC : Reinvent entire module and subroutines. 
!  08/2008  Toshi Matsui@NASA GSFC ; Remove optical properties part into module_opt_micowave
!  07/2007  Toshi Matsui@NASA GSFC ; Add dynamic allocation and better-memory coding.  
!  03/2007  Toshi Matsui@NASA GSFC ; Add options for attenuation and ground-based radar options. 
!  03/2007  Toshi Matsui@NASA GSFC ; Adapted for Goddard SDSU.
!
! References:
!  Masunaga, H., and C.D. Kummerow, 2005: Combined Radar and Radiometer Analysis of 
!    Precipitation Profiles for a Parametric Retrieval Algorithm. J. Atmos. Oceanic 
!    Technol., 22, 909-929.
!
!-----------------------------------------------------------------------------------------------------
 logical,intent(in),optional :: within_swath(:,:)
 integer :: i,j,k,n,nf,m ! looping indice
 integer :: ierr             ! alloc statistics

 real(sdsu_fpd) :: freqhzd, epsreal, epsimag  !dielectric constant parameters
 real(sdsu_fps) :: kradar ! radar constant |K^2|
 real(sdsu_fps) :: umu_radar !consine of incident angle [-]
 integer :: is,ie,js,je,ks,ke !memory domain parameter
 integer :: k_bot, k_top
 real(sdsu_fps) :: dBZ_profile(mxlyr_radar)  ! [dBZ]

!
! simplify loop index name
!
 is=myi_start ; ie=myi_end ; js=myj_start ; je=myj_end

!
! Program start
!
 if(masterproc) print*, 'MSG simulator_radar: Compute radar reflectivity'


 umu_radar = cos(inc_angle_radar*const_degrad) !sensor incident angle

 FREQ_LOOP: do nf = 1, mxfreq_radar

    if(verbose_SDSU .and. masterproc) &
    print*, ' Frequency = ', freq_radar(nf),' GHz -> Wavelength = ',2.997925e2/freq_radar(nf) ,'mm' 


    !
    ! Derive radar constant for different frequency, if not known.  
    !
    freqhzd=dble(freq_radar(nf)*1.e+9)
    call watoptic(freqhzd, dble(273.16e0), dble(0.e0), epsreal, epsimag)
    kradar = ((epsreal-1.e0)**2 + epsimag**2)/ ((epsreal+2.e0)**2 + epsimag**2)
    if( k2(nf) /= undefined ) kradar = k2(nf)
    if(kradar == 0.) call stop_sdsu(myrank,'MSG simulator_radar; kradar should not be zero. Check config file')


    !
    ! horizontal loop
    !
    do j = js,je  ; do i = is,ie

       if( trim(sim_case) == 'MMF_GPROF' .and. j > max_sample_gprof ) cycle

       if( present( within_swath ) ) then
          if( .not. within_swath(i,j) ) cycle !if outside of swath skip rest of all processes. 
       endif

       !
       ! Compute radar echo profile
       !          
       call radar_return( mxlyr_radar, umu_radar, dhgt_radar, atmos_opt(i,j,1:mxlyr_radar,nf) , &
                          2.997925e2/freq_radar(nf), kradar, radar_out(i,j,1:mxlyr_radar,nf)%Zt, &
                          radar_out(i,j,1:mxlyr_radar,nf)%Zm, radar_out(i,j,1:mxlyr_radar,nf)%Vdop, &
                          radar_out(i,j,1:mxlyr_radar,nf)%tau )

    enddo ; enddo

 enddo FREQ_LOOP


#if MPI == 2
if(numproc_tot>1) then

  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_out(ibs:ibe,jbs:jbe,:,:)%Zt )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_out(ibs:ibe,jbs:jbe,:,:)%Zm )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_out(ibs:ibe,jbs:jbe,:,:)%Vdop )
  call mpi_sdsu_communicate_buf( ibs,ibe,jbs,jbe, radar_out(ibs:ibe,jbs:jbe,:,:)%tau  )

endif
#endif

 !
 ! convolution
 !
 if( convolution_radar_on ) then 

   FREQ_LOOP2: do nf = 1, mxfreq_radar

      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), radar_out(:,:,:,nf)%Zt, &
                       radar_out_conv(:,:,:,nf)%Zt   )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), radar_out(:,:,:,nf)%Zm, &
                       radar_out_conv(:,:,:,nf)%Zm   )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), radar_out(:,:,:,nf)%Vdop, &
                       radar_out_conv(:,:,:,nf)%Vdop )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), radar_out(:,:,:,nf)%tau , &
                       radar_out_conv(:,:,:,nf)%tau )

      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%press, &
                       L2_radar_conv(:,:,:)%press  )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%t_air, &
                       L2_radar_conv(:,:,:)%t_air  )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%sh, &
                       L2_radar_conv(:,:,:)%sh  )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%w  , &
                       L2_radar_conv(:,:,:)%w  )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%iwc, &
                       L2_radar_conv(:,:,:)%iwc  )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%lwc, &
                       L2_radar_conv(:,:,:)%lwc  )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%rwc, &
                       L2_radar_conv(:,:,:)%rwc  )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%lfrac, &
                       L2_radar_conv(:,:,:)%lfrac  )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%bulk_den, &
                       L2_radar_conv(:,:,:)%bulk_den  )
      call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%rain_rate, &
                       L2_radar_conv(:,:,:)%rain_rate  )
      do m = 0, 6
         call fov3d_radar(mxlyr_radar, fov_ct_radar(nf), fov_dt_radar(nf), L2_radar_refine(:,:,:)%mmt(m), &
                          L2_radar_conv(:,:,:)%mmt(m) )
      enddo
      

   enddo FREQ_LOOP2

 else !no convolution

   radar_out_conv(is:ie,js:je,:,:)%Zt   = radar_out(is:ie,js:je,:,:)%Zt
   radar_out_conv(is:ie,js:je,:,:)%Zm   = radar_out(is:ie,js:je,:,:)%Zm
   radar_out_conv(is:ie,js:je,:,:)%Vdop = radar_out(is:ie,js:je,:,:)%Vdop
   radar_out_conv(is:ie,js:je,:,:)%tau  = radar_out(is:ie,js:je,:,:)%tau

   L2_radar_conv(is:ie,js:je,:)%press    =L2_radar_refine(is:ie,js:je,:)%press
   L2_radar_conv(is:ie,js:je,:)%t_air    =L2_radar_refine(is:ie,js:je,:)%t_air
   L2_radar_conv(is:ie,js:je,:)%sh       =L2_radar_refine(is:ie,js:je,:)%sh   
   L2_radar_conv(is:ie,js:je,:)%w        =L2_radar_refine(is:ie,js:je,:)%w
   L2_radar_conv(is:ie,js:je,:)%iwc      =L2_radar_refine(is:ie,js:je,:)%iwc  
   L2_radar_conv(is:ie,js:je,:)%lwc      =L2_radar_refine(is:ie,js:je,:)%lwc  
   L2_radar_conv(is:ie,js:je,:)%rwc      =L2_radar_refine(is:ie,js:je,:)%rwc  
   L2_radar_conv(is:ie,js:je,:)%lfrac    =L2_radar_refine(is:ie,js:je,:)%lfrac
   L2_radar_conv(is:ie,js:je,:)%bulk_den =L2_radar_refine(is:ie,js:je,:)%bulk_den  
   L2_radar_conv(is:ie,js:je,:)%rain_rate=L2_radar_refine(is:ie,js:je,:)%rain_rate 

   do m = 0, 6
      L2_radar_conv(is:ie,js:je,:)%mmt(m) = L2_radar_refine(is:ie,js:je,:)%mmt(m) 
   enddo

 endif

! no convolution for HID
   if( radar_hid ) then
       L2_radar_conv(is:ie,js:je,:)%hid      = L2_radar_refine(is:ie,js:je,:)%hid
       L2_radar_conv(is:ie,js:je,:)%hid_conf = L2_radar_refine(is:ie,js:je,:)%hid_conf
   endif


!
!   ---------------------     derive 2D radar output parameter     -------------------------
!
  !
  !PIA at effective level
  !
  do j = js,je ; do i = is,ie

   if( ground_radar ) then

     !
     ! find top layer index (ground based )
     !
     k_top = int(undefined)
     find_k1: do k =  1, mxlyr_radar  !find from bottom
          if( radar_out_conv(i,j,k,1)%tau /= undefined ) then
              k_top = k
              exit find_k1
          endif
     enddo find_k1

     !if cannot find k_top, skip 
     if( k_top == int(undefined) ) then
          radar_out_2d_conv(i,j,:)%pia = undefined
     else

        !
        ! Derive PIA for clear-sky and cloudy sky
        !
        do nf = 1, mxfreq_radar
           if( radar_out_conv(i,j, k_top,nf )%Zm == 0. ) then !if no return use tau
              radar_out_2d_conv(i,j,nf)%pia = 20.* radar_out_conv(i,j,k_top,nf)%tau * log10(2.71828183)
           else  !if echo is significant, use dBZt-dBZm
              radar_out_2d_conv(i,j,nf)%pia = &
                 10.* log10( radar_out_conv(i,j,k_top,nf)%Zt / radar_out_conv(i,j,k_top,nf)%Zm )
           endif
        enddo

     endif


   else

     !
     ! find bottom layer index (space born)
     !
     k_bot = int(undefined)
     find_k2: do k =  1, mxlyr_radar  !find from bottom
          if( radar_out_conv(i,j,k,1)%tau /= undefined ) then
              k_bot = k
              exit find_k2
          endif
     enddo find_k2

     !if cannot find k_bot, skip 
     if( k_bot == int(undefined) ) then
          radar_out_2d_conv(i,j,:)%pia = undefined
     else 

        !
        ! Derive PIA for clear-sky and cloudy sky
        !
        do nf = 1, mxfreq_radar
           if( radar_out_conv(i,j, k_bot,nf )%Zm == 0. ) then !if no return use tau
              radar_out_2d_conv(i,j,nf)%pia = 20.* radar_out_conv(i,j,k_bot,nf)%tau * log10(2.71828183)
           else  !if echo is significant, use dBZt-dBZm
              radar_out_2d_conv(i,j,nf)%pia = &
                 10.* log10( radar_out_conv(i,j,k_bot,nf)%Zt / radar_out_conv(i,j,k_bot,nf)%Zm )
           endif
        enddo

     endif

   endif

  enddo ; enddo


!
! derive maximum composite
!
 do nf = 1, mxfreq_radar ; do j = js,je ; do i = is,ie 

     radar_out_2d_conv(i,j,nf)%Zt_max = maxval( radar_out_conv(i,j,:,nf)%Zt ) ![Z]
     radar_out_2d_conv(i,j,nf)%Zm_max = maxval( radar_out_conv(i,j,:,nf)%Zm ) ![Z]

 enddo ; enddo ; enddo


!
! derive echo-top height , based on Zt (above ground level)
!
 do nf = 1, mxfreq_radar ; do j = js,je ; do i = is,ie

     do k =  1, mxlyr_radar
        if( radar_out_conv(i,j,k,nf)%Zt > 0.e0 ) then
            dBZ_profile(k) = 10.e0 * log10( radar_out_conv(i,j,k,nf)%Zt )   ![dBZ]
        else
            dBZ_profile(k) = undefined 
        endif
     enddo

     !
     ! find echo-top
     !
     k_top = int(undefined)
     find_top: do k =  mxlyr_radar,3,-1
        if( dBZ_profile(k  ) >= min_echo .and. &
            dBZ_profile(k-1) >= min_echo .and. &
            dBZ_profile(k-2) >= min_echo ) then
            k_top = k
            exit find_top
        endif
     enddo find_top

     if( k_top == int(undefined) ) then
         
         radar_out_2d_conv(i,j,nf)%Het = undefined 

     else

         !
         ! correct echo-top height from above-sea-level to above-ground-level
         !
         radar_out_2d_conv(i,j,nf)%Het = hgt_stag_radar( k_top ) - atmos_stag(i,j,0)%hgt  ! [km]

     endif

 enddo ; enddo ; enddo


 do nf = 1, mxfreq_radar ; do k =  1, mxlyr_radar ; do j = js,je ; do i = is,ie
    if( radar_out_conv(i,j,k,nf)%Zt > 0.e0 ) then
        dBZ_profile(k) = 10.e0 * log10( radar_out_conv(i,j,k,nf)%Zt )   ![dBZ]
        if( dBZ_profile(k) < min_echo ) then
            radar_out_conv(i,j,k,nf)%Zt = 0.e0
            radar_out_conv(i,j,k,nf)%Zm = 0.e0
            radar_out_conv(i,j,k,nf)%Vdop = undefined
        endif
    endif
 enddo ; enddo ; enddo ; enddo


 return
 end subroutine simulator_radar 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine radar_return_scan( kmax, dhgt, kext, sback, lambda, &
                          kradar, Zt, Zm, tau )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Simulating attenuating and non-attenuating radar reflectivity factor for a space-based instrument. 
!  * Unlike radar_echo routine, vertical index k = 1 is TOA, and k = kmax is BOA. 
!  * This does not account multiple scattering. 
!  * Z == zero, if sback == zero. 
!  * Z == undefined, if sback == undefined (above CRM TOA, or below CRM BOA). 
!
! 
! History:
! 11/2012  Toshi Matsui@NASA GSFC ; output modified from Zm to tau. 
! 04/2012  Toshi Matsui@NASA GSFC ; adapted from radar_echo routine          
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)        :: kmax         !maximum vertical index
 real(sdsu_fps),intent(in) :: dhgt(kmax)   !layer thickness [km] 
 real(sdsu_fps),intent(in) :: kext(kmax)   !total extinction [km-1] 
 real(sdsu_fps),intent(in) :: sback(kmax)  !total backscatter cross section [km-1]  (<- [m2/m3] with unit change )
 real(sdsu_fps),intent(in) :: lambda       !Radar echo wavelength [mm]
 real(sdsu_fps),intent(in) :: kradar       !radar constant |k^2|
 real(sdsu_fps),intent(out) :: Zt(kmax)    !non-attenuating radar reflectivity [mm6/m3] (below TOA and BOA = undefined)
 real(sdsu_fps),intent(out) :: Zm(kmax)    !attenuating radar reflectivity [mm6/m3] (below TOA and BOA = undefined)
 real(sdsu_fps),intent(out) :: tau(kmax)  ! path-integrated ptical depth [-]

 integer :: k      ! vertical loop index
 integer :: ks, ke  !starting/ending k index
 real(sdsu_fps) :: dtau

!
! initialize output
!
  Zt(1:kmax) = undefined ; Zm(1:kmax) = undefined  ; tau(1:kmax) = undefined

!
! bad input
!
 if( maxval( kext(:) ) <= 0.e0 ) then
   return
 endif 



 if(ground_radar) then  !radar reflectivity measured from ground
   !
   ! findout starting and ending k index
   !
   ks = undefined_i2  ; ke = undefined_i2

   do k = 1,kmax
      if( kext(k) >= 0. ) then 
         ke = k ; exit
      endif
   enddo
   do k = kmax, 1, -1
      if( kext(k) >= 0. ) then  
         ks = k ; exit
      endif
   enddo
  !
  ! all undefined --> return
  !
   if(ks == undefined_i2 .or. ke == undefined_i2 ) then
      return
   endif


   do k = ke, ks, -1  !ground-based remote sensing

       !
       ! path-integrated ptical depth [-]
       !
       tau(k) = SUM(kext(k:ke)* dhgt(k:ke)) - 0.5*kext(k)* dhgt(k)  

       !
       ! yet defined top and bottom layer, there are still strange tau value (very very few)
       ! in the middle of layer. So remove this noise. 
       !
       if( tau(k) < 0.e0 )  tau(k) = 0.e0

       if (sback(k) <= 0.e0) then
           Zt(k) = 0.e0  ; Zm(k) = 0.e0
       else

           !
           ! non-attenuating Z [mm6/m3]
           !
           Zt(k) = ( (lambda**4.)*sback(k)*1000. )/ ( kradar * (const_pi**5) )   

           !
           ! Attenuating Z [mm6/m3]
           !
           Zm(k) = Zt(k) * exp(-2.*tau(k))

       endif !sback >0
    enddo !k

 else !radar reflectivity from space
    !
    ! findout starting and ending k index
    !
    ks = undefined_i2  ; ke = undefined_i2

    do k = 1,kmax
       if(kext(k) >= 0. ) then 
          ks = k ; exit
       endif
    enddo
    do k = kmax, 1, -1
       if(kext(k) >= 0. ) then 
         ke = k ; exit
       endif
    enddo
  !
  ! all undefined --> return
  !
   if(ks == undefined_i2 .or. ke == undefined_i2 ) then
      return
   endif

    do k = ks, ke, 1 !space-born satellite

       !
       ! path-integrated ptical depth [-]
       !
       tau(k) = SUM(kext(ks:k)* dhgt(ks:k)) - 0.5*kext(k)* dhgt(k) 

       !
       ! yet defined top and bottom layer, there are still strange tau value (very very few)
       ! in the middle of layer. So remove this noise. 
       !
       if( tau(k) <= 0.e0 ) tau(k) = 0.

       if (sback(k) <= 0.e0) then
           Zt(k) = 0.e0  ; Zm(k) = 0.e0
       else

           !
           ! non-attenuating Z [mm6/m3]
           !
           Zt(k) = ( (lambda**4.)*sback(k)*1000. )/ ( kradar * (const_pi**5) )   

           !
           ! Attenuating Z [mm6/m3]
           !
           Zm(k) = Zt(k) * exp(-2.*tau(k))

       endif !sback > 0
    enddo !k

   endif !ground_radar

 return
 end subroutine radar_return_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine radar_return( kmax, umu_radar, dhgt, opt1d, lambda, &
                          kradar, Zt, Zm, Vdop, tau )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Simulating attenuating and non-attenuating radar reflectivity factor for a space-based instrument. 
!  * Unlike radar_echo routine, vertical index k = 1 is TOA, and k = kmax is BOA. 
!  * This does not account multiple scattering. 
!  * Z == zero, if sback == zero. 
!  * Z == undefined, if sback == undefined (above CRM TOA, or below CRM BOA). 
!
! 
! History:
! 11/2012  Toshi Matsui@NASA GSFC ; output modified from Zm to tau. 
! 04/2012  Toshi Matsui@NASA GSFC ; adapted from radar_echo routine          
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)        :: kmax         ! maximum vertical index
 real(sdsu_fps),intent(in) :: umu_radar    ! cosine of incident angle
 real(sdsu_fps),intent(in) :: dhgt(kmax)   ! layer thickness [km] 
 type(opt_micro_parameter),intent(in) :: opt1d(kmax) !optical properties profiles
 real(sdsu_fps),intent(in)  :: lambda      ! Radar echo wavelength [mm]
 real(sdsu_fps),intent(in)  :: kradar      ! radar constant |k^2|
 real(sdsu_fps),intent(out) :: Zt(kmax)    ! non-attenuating radar reflectivity [mm6/m3] (below TOA and BOA = undefined)
 real(sdsu_fps),intent(out) :: Zm(kmax)    ! attenuating radar reflectivity [mm6/m3] (below TOA and BOA = undefined)
 real(sdsu_fps),intent(out) :: Vdop(kmax)  ! doppler velocity [m/s]
 real(sdsu_fps),intent(out) :: tau(kmax)   ! path-integrated optical depth [-]


 integer :: k      ! vertical loop index
 integer :: ks, ke  !starting/ending k index
 real(sdsu_fps) :: dtau
 real(sdsu_fps) :: pia         ! path-integrated ptical depth [-]

!
! initialize
!
  Zt = undefined ; Zm = undefined  ; tau = undefined ; Vdop = undefined

 if(ground_radar) then  !radar reflectivity measured from ground
   !
   ! findout starting and ending k index
   !
   ks = undefined_i2  ; ke = undefined_i2

   do k = 1,kmax
      if( opt1d(k)%kext >= 0. ) then 
         ks = k ; exit
      endif
   enddo
   do k = kmax, 1, -1
      if( opt1d(k)%kext >= 0. ) then  
         ke = k ; exit
      endif
   enddo
  !
  ! all undefined --> return
  !
   if(ks == undefined_i2 .or. ke == undefined_i2 ) then
      return
   endif


   do k = ks, ke  !ground-based remote sensing

       !
       ! path-integrated ptical depth [-]
       !
       tau(k) = (SUM(opt1d(ks:k)%kext* dhgt(ks:k)) - 0.5*opt1d(k)%kext* dhgt(k)) / umu_radar

       if (opt1d(k)%sback <= 0.) then
           Zt(k) = 0.e0  ; Zm(k) = 0.e0 ; Vdop(k) = 0.e0
       else

           !
           ! non-attenuating Z [mm6/m3]
           !
           Zt(k) = ( (lambda**4.)*opt1d(k)%sback*1000. )/ ( kradar * (const_pi**5.) )   

           !
           ! Attenuating Z [mm6/m3]
           !
           Zm(k) = Zt(k) * exp(-2.*tau(k))

           !
           ! Doppler velocity
           !
           Vdop(k) = opt1d(k)%vdop  ![m/s]

       endif !sback >0

    enddo !k
    pia = 20.* tau(ke) * log10(2.71828183)

 else !radar reflectivity from space

    !
    ! findout starting and ending k index
    !
    ks = undefined_i2  ; ke = undefined_i2

    do k = 1,kmax
       if(opt1d(k)%kext >= 0. ) then 
          ks = k ; exit
       endif
    enddo
    do k = kmax, 1, -1
       if(opt1d(k)%kext >= 0. ) then 
         ke = k ; exit
       endif
    enddo
  !
  ! all undefined --> return
  !
   if(ks == undefined_i2 .or. ke == undefined_i2 ) then
      return
   endif

    do k = ke, ks, -1 !space-born satellite

       !
       ! path-integrated ptical depth [-]
       !
       tau(k) = (SUM( opt1d(k:ke)%kext * dhgt(k:ke)) - 0.5*opt1d(k)%kext * dhgt(k) ) / umu_radar 

       if( tau(k) <= 0. ) then
           print*,'MSG radar_return: tau become negative.',k, tau(k)
!           print*,'kext(ks:k)=',kext(ks:k)
!           print*,'dhgt(ks:k)=',dhgt(ks:k)
           tau(k) = 0.
       endif

       if (opt1d(k)%sback <= 0.e0) then
           Zt(k) = 0.e0  ; Zm(k) = 0.e0 ; Vdop(k) = 0.
       else

           !
           ! non-attenuating Z [mm6/m3]
           !
           Zt(k) = ( (lambda**4.)*opt1d(k)%sback*1000. )/ ( kradar * (const_pi**5.) )   

           !
           ! Attenuating Z [mm6/m3]
           !
           Zm(k) = Zt(k) * exp(-2.*tau(k))

           !
           ! Doppler velocity
           !
           Vdop(k) = opt1d(k)%vdop  ![m/s]

       endif !sback > 0
    enddo !k
    pia = 20.* tau(ke) * log10(2.71828183)


   endif !ground_radar

 return
 end subroutine radar_return

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine attenuate_z (kmax, Zt, tau, PIA, Zm )
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:
!  Simulate attenuation of radar reflectivity for a given optical depth profile.   
!  This routine must be called after convolution of radar reflecitivty (Z). 
! 
! History:
! 11/2012  Toshi Matsui@NASA GSFC ; Init. 
! References: 
!-----------------------------------------------------------------------------------------------------

 integer,intent(in) :: kmax
 real(sdsu_fps),intent(in) :: Zt(kmax)        ! non-attenuating equiv. reflectivity factor [Z] = [mm6/m3] = [um3]
 real(sdsu_fps),intent(in) :: tau(kmax)       ! optical depth [-]
 real(sdsu_fps),intent(out) :: Zm(kmax)        ! attenuating     equiv. reflectivity factor [Z] = [mm6/m3] = [um3]
 real(sdsu_fps),intent(out) :: PIA(kmax)       ! 2-way path integrated attenuation [dB]

 integer :: k, ks, ke
 real(sdsu_fps) :: tau_tot

!
! initialize output. 
!
 Zm(:) = undefined ; PIA(:) = undefined

 !
 ! findout starting and ending k index
 !
 ks = undefined_i2  ; ke = undefined_i2

 do k = 1,kmax
    if(Zt(k) /= undefined ) then
       ks = k ; exit
    endif
 enddo
 do k = kmax, 1, -1
    if(Zt(k) /= undefined ) then
       ke = k ; exit
    endif
 enddo

 !
 ! all undefined --> return
 !
 if(ks == undefined_i2 .or. ke == undefined_i2 ) then
    return
 endif


 !
 ! derive PIA and Zm profile
 !
 do k = ks, ke

    !tau_tot = SUM(tau(ks:k)) - 0.5*tau(k)  
    tau_tot = tau(k)  !toshii

    !
    ! 2-way path integrated attentuation profile [dB]
    !
    PIA(k) = 20.0*tau_tot*log10(2.71828183)

    !
    ! attenuating reflectivity profile [mm6/m3]
    !
    Zm(k) = Zt(k) * exp(-2.*tau_tot )

 enddo

 return
 end subroutine attenuate_z

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 function dB(inval) result(value)
     real        :: inval, value
     if( inval <= 0. ) then
       value = undefined
     else
       value = 10. * alog10( max(inval,1.0e-35 ) )
     endif
 end function dB

 end module module_radar

