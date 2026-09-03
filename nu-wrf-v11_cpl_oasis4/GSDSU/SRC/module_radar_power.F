#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_radar_power
 use module_floating_point
 use module_const
 implicit none
!
!
!!A generalized surface ehoc radar equation for down-looking pencil beam radar, T. Kozu, IEICE Trans, Coomun. 
!1995
!TRMM Precipitatioon Radar algorithm Instruction Manual For Version 6. 
!
!
!
! The pulse repetition interval is adjusted according to the satellite algtitude and the angle of 
! observations. As a result, the number of indpeendent samples changes slightly as a function of the
! scan angle. However, for simulator, we define the distance of starting ragne and ending range identical
! to all scan angle, which makes easy for IFOV convolution. 
!
!

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 private   ! ALL variables and subourtines are non-accessible, except....

 public :: return_radar_power  !public subroutine

 public :: radar_power_out

!
! output of radar power routine
!
 type radar_power_out
    real :: totk      ! path-integrated attenuation (2way) [dB]
    real :: avgk      ! attenuation [dB/km] 
    real :: Pru       ! non-attenuating radar power [W]
    real :: Pr        ! attenuating radar power [W]
    real :: Zt        ! non-attenuating reflectivity  [Z] = [mm6/m3] = [um3]
    real :: Zm        ! attenuating reflectivity  [Z] = [mm6/m3] = [um3]
    real :: Psu       ! non-attenuating radar power from surface [W]
    real :: Ps        ! attenuating radar power from surface [W]
!    real :: p_noise   ! noise power [W]
!    real :: p_total   ! total power [W]
!    real :: p_total_u ! non-attenuating total power [W]
!    real :: p_hat     ! total power wihtout noise [W]
 end type radar_power_out

! real, parameter       :: c  = 2.997925E+08*1.0E-09  ! ! The speed of light is given in meter*GHz

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine return_radar_power( irange_max, range_distance, kext, sback, zone, lambda_cm, k2_radar, loss, ptp, &
                         range_resolution, sigma0, range1d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:   
!            * too week surface power return --> must be -50 ~ -90dBm
!            * too week atmos power return ??
!            * too strong molecular attenuation??
! History:
! 04/2012  Toshi Matsui@NASA GSFC    : initial (adapted from Bob/Hyukyong's code)
!
! References:
! R. Meneghini, 1995: Analysis of Radar and Radar-Radiometer Methods For Spaceborne
!      Measurments of Precipitation, Ph.D. Dissertation, Kyoto University. 
!-------------------------------------------------------------------------------------------------
 integer,intent(in) :: irange_max      ! maximum range number
 real,intent(in) :: range_distance(irange_max)    ! radar range [m]
 real,intent(in) :: kext(irange_max)   ! extinction [1/km]
 real,intent(in) :: sback(irange_max)  ! backscatter coefficient [1/km]
 character(len=3),intent(in) :: zone(irange_max) !'top' 'mid' 'bot'
 real,intent(in) :: lambda_cm        ! wavelength [cm]
 real,intent(in) :: k2_radar         ! radar constant |k^2|
 real,intent(in) :: loss             ! loss of power [W]
 real,intent(in) :: ptp              ! peak transmit power [W] (1012.0Watts for KuPR 146.5Watts for KaPR)
 real,intent(in) :: range_resolution ! range resolution  [m]  (DPR Ku)  must be km??
 real,intent(in) :: sigma0           ! normalized scattering cross section of the surface [dB] (area per unit area)
 type ( radar_power_out ),intent(out) :: range1d(irange_max) !output parameters

 integer :: irange       ! range loop idex
 real :: lambda_m        ! wavelength [m]
 real :: const_atmos     ! conversion from power to Z <-- radar equation
 real :: const_surface   ! conversion constant factor
 real :: attenuation     ! attenuation  [dB/m] 
 real :: k_prev,k_cur    ! attenuation  [dB/m] 
 real :: kbar            ! path-integrated attenuation  [dB] 
 real :: decay           ! decay due to 2way attenutation 
 real :: range_m          ! range distance [m]
 real :: z_atmos         ! equivalent reflectivity factor of atmosphere [Z] = [mm6/m3] = [um3]
 real :: z_surf          ! equivalent reflelctivity factor of surface   [Z] = [mm6/m3] = [um3]
 real :: range_volume   ! range volume [m3]
 real :: range_area   ! range area [m2]

! real,intent(in) :: fnoise
! real,intent(in) :: band
! integer,intent(in) :: nbin ! ??
! real,intent(in) :: ag      ! ??
! real :: fr, fn   ! ramdome function ???
! real :: pn      ! noise power [W]

!
! initialize all output parameters 
!
  range1d%totk      = 0. ! path-integrated attenuation (2way) [dB]
  range1d%avgk      = 0. ! attenuation [dB/km] 
  range1d%pru       = 0. ! non-attenuating radar power [W]
  range1d%pr        = 0. ! attenuating radar power [W]
  range1d%zt        = 0. ! non-attenuating reflectivity  [Z] = [mm6/m3] = [um3]
  range1d%zm        = 0. ! attenuating reflectivity  [Z] = [mm6/m3] = [um3]
  range1d%psu       = 0. ! non-attenuating radar power from surface [W]
  range1d%ps        = 0. ! attenuating radar power from surface [W]


 lambda_m = lambda_cm / 100.      !wavelength  [m]   <-- [cm]

!
! 1.e-18 is conversion from [um3] to [m3].
!
 range_volume = 6.1683460E+18  !const_pi * 2500. * 2500. * range_resolution ![m3]
! print*,range_volume, dB( sqrt( range_volume ) ), (10**(93.95084/10.))**2
! stop

 const_atmos =  range_volume * loss * k2_radar * (const_pi/lambda_m)**2 * ptp / 64. * 1.0e-18  ! from eq 8 & 9, A3.4

 range_area   = const_pi * 5000. * 5000. ![m2]    
 const_surface  = range_area * loss * (lambda_m**2) * ptp / (4.*const_pi)**3  ! [m4] from eq 12 & 13, A3.5


 kbar = 0.
 do irange = 1 , irange_max       ! radial integral

    range_m = range_distance(irange)    ! radar range distance [m]

    select case(zone(irange)) !'top' 'mid' 'bot'

    case('top','mid') !above CRM height

      ! 
      ! derive reflectivity (z_atmos = [Z] = [mm6/m3] = [um3] )
      ! and specific attenuation (k= [dB/m] )
      !
      if( kext(irange) == undefined) then
        z_atmos = 0. ; attenuation = 0.
      else
        call get_Z(lambda_cm, k2_radar, sback(irange)*1.e-5, kext(irange)*1.e-5, z_atmos, attenuation)
      endif

      k_cur = attenuation                                 ! attenuation in current step [dB/m]
      kbar = kbar + 0.5 * (k_prev+k_cur) * range_resolution  ! path-integrated attenuation  [dB] 
      k_prev = k_cur                                      ! attenuation in previous step [dB/m]

      range1d(irange)%totk = kbar                         ! path-integrated attenuation  [dB]
      range1d(irange)%avgk = attenuation * 1.e+3          ! specific attenuation [dB/km] 

      !
      ! atmopsheric power output:  dBm(power) ranges from -150mW to -50mW.
      !
      decay = exp( -0.2*log(10.0) * kbar  )  ! decay due to 2way path-integrated attenuation
      range1d(irange)%pru  = const_atmos         * z_atmos / (range_m**4)  ! non-attenuating radar power [W]
      range1d(irange)%pr   = const_atmos * decay * z_atmos / (range_m**4)  ! attenuating radar power [W]

      !
      ! equivalent reflectivity factor [Z] = [mm6/m3] = [um3]
      !
      range1d(irange)%zt = (range_m**4) * range1d(irange)%pru / const_atmos  ! non-attenuating reflectivity
      range1d(irange)%zm = (range_m**4) * range1d(irange)%pr  / const_atmos  ! attenuating reflectivity 

      print*, 'mid',irange, dBm (range1d(irange)%pru ), dB( range1d(irange)%zt ) , decay, kbar

    case('bot') ! in ground 

!      attenuation = sigma0/range_resolution              ! attenuation due to surface return [dB/m]
!      k_cur = attenuation                                ! attenuation in current step [dB/m]
!      kbar = kbar + 0.5 * (k_prev+k_cur) * range_resolution  ! PIA [dB]
!      k_prev = kbar                               ! attenuation in previous step [dB/m]
!      range1d(irange)%totk = kbar                 ! path-integrated attenuation  [dB]
!      range1d(irange)%avgk = attenuation * 1.e+3  ! specific attenuation [dB/km] 
   
      !
      ! surface signals from from eq (12), A3.4
      !
      decay = exp( -0.2*log(10.0) * kbar  )  ! decay due to 2way path-integrated attenuation
      range1d(irange)%psu = const_surface         * sigma0 / (range_m**4.) ! non-attenuating surafce return power [W]
      range1d(irange)%ps  = const_surface * decay * sigma0 / (range_m**4.) ! attenuating     surface return power [W]

      z_surf = (range_m**4) * range1d(irange)%ps  / const_atmos

      print*,'bot', irange, dBm ( range1d(irange)%ps ), dB( z_surf ), decay, kbar
      exit

    case default 
       stop
    end select


    !
    ! add noise and get total return power (this should be implemented after convolution)
    !
     ! call noise_power(fnoise, band, pn)  !get noise power

    ! ramdome power
    !    call randomization_power(fr,fn)

    !
    ! final output with random power
    !
!    range1d(ibin)%p_noise   = pn * fn
!    range1d(ibin)%p_total   = ( range1d(ibin)%pr + range1d(ibin)%ps + pn ) * fr
!    range1d(ibin)%p_total_u = ( range1d(ibin)%pru+ range1d(ibin)%psu + pn ) * fr
!    range1d(ibin)%p_hat     = range1d(ibin)%p_total - range1d(ibin)%p_noise

!    write(unit=pout,fmt='(4f10.5,11f12.5)') ran(j)/1000., xst(j)/1000., yst(j)/1000., zst(j)/1000.,&
!dBm(pru(j)),dBm(pr(j)),dBm(psu(j)),dBm(ps(j)),dBm(pn),fr,dBm(p_total_u(j)),dBm(p_total(j)),fn,dBm(p_noise(j)),dBm(p_hat(j))

!    write(unit=out,fmt='(4f10.5,18f12.5)') ran(j)/1000., xst(j)/1000., yst(j)/1000., zst(j)/1000.,&
!    dB(zm(j)),dB(zt(j)),avgk(j),totk(j),avgrr(j),avgtemp(j),avgpres(j),(avgmass(j,i),i=1,8),(avgD0(j,i),i=1,3)

!    if(zst(j)/1000. .gt. 3. .and. zst(j)/1000. .lt. 3.1) then
!       write(*,fmt='(a20,f7.2,i5,7f12.5)') 'binres=',binres,j,ran(j)/1000., xst(j)/1000., yst(j)/1000., zst(j)/1000.,&
!       dBm(pr(j)),dBm(pn),dBm(pn+pr(j))
!    endif

 enddo


 return
 end subroutine return_radar_power

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 function deg2rad(inval) result(value)
    real        :: inval, value
    value = inval *const_pi/180.
 end function deg2rad

 function rad2deg(inval) result(value)
    real        :: inval, value
    value = inval *180./const_pi
 end function rad2deg

 function dBm(inval) result(value)
     real        :: inval, value
!     if( inval <= 0. ) then
!       value = undefined
!     else
       value = 10. * alog10( max(inval,1.0e-35 ) )+30.
!     endif

 end function dBm

 function dB(inval) result(value)
     real        :: inval, value
!     if( inval <= 0. ) then
!       value = undefined
!     else
       value = 10. * alog10( max(inval,1.0e-35 ) )
!     endif
 end function dB

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 
 subroutine get_Z ( lambda_cm, kradar, sback_cm, kext_cm, reflect, atten )
 real,intent(in) :: lambda_cm! wavelength [cm]
 real,intent(in) :: kradar   ! radar constant |K|^2 
 real,intent(in) :: sback_cm ! backscatter coefficient [1/cm]
 real,intent(in) :: kext_cm  ! extinction coefficient  [1/cm]

 real,intent(out) :: reflect  ! equivalent reflectivity factor [mm6/m3]
 real,intent(out) :: atten    ! specific attenuation  [dB/m]

 real,parameter :: cm3_to_micron3 =  1.0e+12       !constatnt for unit conversion from [cm3] to [micron3]
 real,parameter :: conK = 4.343*1.e+06  ! eq. (6) A3.3
 real,parameter :: lne = 0.43429448e0 ! = log10(2.71828183)


 reflect = cm3_to_micron3 * (lambda_cm**4)/(kradar* (const_pi**5) ) * sback_cm  ! [Z] = [mm6/m3] = [micron3]

 atten   = conK * kext_cm * 1.e-3  ! specific attenuation  [dB/m]
                                   ! (1e+6*0.433 from 1e+6/km*4.33 dB/neper) <-- ??
                                   ! conK * kext_cm = [dB/km]

!toshii attenuation is too large?? 

 return

 end subroutine get_Z

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine noise_power(fnoise,band,pn)
 implicit none
 real,intent(in) :: fnoise
 real,intent(in) :: band
 real,intent(out) :: pn

 real, parameter :: const_boltzmann = 1.380658e-23
 real, parameter :: receiver_temp   = 290.
 real            :: fpower

 fpower  = 10.**(fnoise/10.)    ! fnoise, band comes from namelist.input
 pn = const_boltzmann * receiver_temp * fpower * band  !receiver_temp=290K

 end subroutine noise_power

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


 subroutine randomization_power(fr,fn)
 implicit none

 real,intent(out) :: fr, fn

 real :: z1,z2,mu1,mu2,std1,std2,a,b
 integer :: num1,num2,ix
 integer,parameter :: nrain=100,nnoise=400
 real, dimension(nrain) :: u
 real, dimension(nnoise) :: v

 mu1=1.
 std1=1./sqrt(1.*nrain)
 mu2=1.
 std2=1./sqrt(1.*nnoise)
 do ix=1,nrain
    call random_normal(z1)
    u(ix)=mu1+std1*z1
 enddo

 do ix=1,nnoise
    call random_normal(z2)
    v(ix)=mu2+std2*z2
 enddo

 call random_number(a)
 num1=1+nint(a*99.)

 call random_number(b)
 num2=1+nint(b*399.)
 fr=u(num1)
 fn=v(num2)
 return
 end subroutine randomization_power

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine random_normal(fn_val)
 implicit none
!
! adapted from the following fortran 77 code
! algorithm 712, collected algorithms from acm.
! this work published in transactions on mathematical software,
! vol. 18, no. 4, december, 1992, pp. 434-435.
! the function random_normal() returns a normally distributed pseudo-random
! number with zero mean and unit variance.
! the algorithm uses the ratio of uniforms method of a.j. kinderman
! and j.f. monahan augmented with quadratic bounding curves.
!
 real :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0
 real :: vsmall = tiny(1.0), vlarge = huge(1.0)
 integer, parameter :: dp = selected_real_kind(12, 60)
 real :: fn_val

! local variables
 real ::  s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472
 real :: r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

 do
  call random_number(u)
  call random_number(v)
  v = 1.7156 * (v - half)

! evaluate the quadratic form
  x = u - s
  y = abs(v) - t
  q = x**2 + y*(a*y - b*x)
  if (q < r1) exit                  !     accept p if inside inner ellipse
  if (q > r2) cycle                 ! reject p if outside outer ellipse
  if (v**2 < -4.0*log(u)*u**2) exit ! reject p if outside acceptance region
 end do

 fn_val = v/u                        ! return ratio of p's coordinates as the normal deviate
 return
 end subroutine random_normal

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_radar_power

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
