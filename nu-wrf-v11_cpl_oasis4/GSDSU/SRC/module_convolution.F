#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_convolution
 use module_floating_point
 use module_const
 implicit none

!--------------------------------------------------------------------------------------------------
! Comments:                             
!   Various subroutines to estimate antenna gain function of satellite for non-scan simulation. 
!   
!  History:
!  04/2012  Toshi Matsui@NASA GSFC : Initial developed for MMF_GPROF database.  
!
! References:
!-----------------------------------------------------------------------------------------------------

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 private   ! ALL variables and subourtines are non-accessible.

!
! Accessible public subroutines
!
 public :: beam_width_dx         ,& ! 
           gain_function_gaussian,& ! calculate antenna gain function under Ganssian assumption. 
           derive_incidence_angle,& ! derive incidence angle and sat-fov distance
           beam_width_9x7plane      ! estimate beam width in 9x78 2D plance


 real(sdsu_fpd),parameter :: dpi=3.14159265358979323846d0  ! pai 

 contains 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine derive_incidence_angle( sva, sma, lat, incidence_ang,  dist_sat_fov )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For given satellite lat, semi-major axix, sensor viewing angles, this routines
!            estiamte incidence angle and satellite-FOV distance. 
!
! History:
! 09/2012  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fps), intent(in) :: sva ! sensor viewing angle [rad]
 real(sdsu_fps), intent(in) :: sma   ! semi-major axis of satellite [km]
 real(sdsu_fps), intent(in) :: lat   ! latitude of satellite [rad]
 real(sdsu_fps), intent(out) :: incidence_ang ! incidence angle [rad]
 real(sdsu_fps), intent(out) :: dist_sat_fov ! distance between satellite and FOV [km]

 real(sdsu_fpd) :: a,b,c                   !triangle length [km]
 real(sdsu_fpd) :: alpha, beta, gam, theta  !triangle angle [rad]

 alpha = dble(sva) !sensor view angle [rad]
 b = dble(sma)     !semi-major axis of satellite [km]

! Earth radius 
!
 call earth_radius( dble(lat), a )

 if( alpha == 0.d0 ) then

   c  = b - a ![km]
   theta = 0.d0 ![rad]

 else

   beta = dpi - asin(b/a*sin(alpha))  ! For Fotran it must be dpi -  ! [rad]

   theta = dpi - beta ! incidence angle [rad] 

   gam = dpi - alpha - beta  !earth-core angle  [rad]

   c   = dsin(gam)/dsin(beta) * b  ! distance between satellite and fov [km]

 endif

!
! output
!
 incidence_ang = real(theta)  ! incidence angle [rad]

 dist_sat_fov  = real(c)      ! distance between satellite and FOV [km]

 return
 end  subroutine derive_incidence_angle

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine beam_width_9x7plane( dgrid, incidence_ang, dist_sat_fov, bw_deg )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Estimate unique beam width on 9x7 2D plance.  
!
! History:
! 09/2012  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fps), intent(in) :: dgrid         ! grid spacing [km]
 real(sdsu_fps), intent(in) :: incidence_ang ! incidence angle [rad]
 real(sdsu_fps), intent(in) :: dist_sat_fov  ! distance between satellite and FOV [km]
 real(sdsu_fps), intent(out) :: bw_deg(64)       ! beam width [deg]

 integer :: i,j, i_1d
 real(sdsu_fps) :: dist_2d(9,7) !distance from the center of 2D plane

 real(sdsu_fpd),dimension(3) :: v_sat, v_fov, v_grid ! 3D position vector [km]

 real(sdsu_fpd) :: a,b,c                   !triangle length [km]
 real(sdsu_fpd) :: alpha, beta, gam, theta  !triangle angle [rad]

!
! initilaize output
!
 bw_deg(:) = undefined

!
! FOV center is treated as center of 3D Cartesian grid.
!
 v_fov(1:3) = 0.d0  ! [km]

!
! Satellite position
!
 v_sat(1) = dble( dist_sat_fov * sin(incidence_ang) )  ![km]
 v_sat(2) = 0.d0 !y-component is zero [km]
 v_sat(3) = dble( dist_sat_fov * cos(incidence_ang) )  ![km]

!
! Note: Center of FOV is located in the center of 8x8 2D plane.
!
 i_1d = 0
 do j = 1, 7 ; do i = 1, 9 
    v_grid(1) = dble( abs( 5. - real(i) ) * dgrid ) !x-distance [km]
    v_grid(2) = dble( abs( 4. - real(j) ) * dgrid ) !y-distance [km]
    v_grid(3) = 0.d0

    !
    ! At this point, we know all components of v_sat, v_fov, v_grid. 
    ! so get triangle parameters. 
    !
    call triangle_parameter_scan(v_sat, v_fov, v_grid, alpha, beta, gam, a, b, c )

    i_1d = i_1d + 1
    
    ! output
    bw_deg(i_1d) = real(alpha)/const_degrad  ! beam width [deg]

 enddo ; enddo


 return
 end subroutine beam_width_9x7plane

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
 call get_distance( v_sat   ,   v_grid2, b )
 call get_distance( v_sat   ,   v_grid1, c )

 alpha = dacos( ( b*b + c*c - a*a ) / (2.*b*c) )  ![rad]
 beta  = dacos( ( a*a + c*c - b*b ) / (2.*a*c) )  ![rad]
 gam   = dacos( ( a*a + b*b - c*c ) / (2.*a*b) )  ![rad]

 return
 end subroutine  triangle_parameter_scan

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

 subroutine beam_width_dx(lat_deg, sma , ona_deg, delta_x, bw_deg )
 implicit none

 real(sdsu_fps), intent(in) :: lat_deg  !latitude [deg]
 real(sdsu_fps), intent(in) :: sma      !semi-major axis [km]
 real(sdsu_fps), intent(in) :: ona_deg  !off-nadir angle [deg]
 real(sdsu_fps), intent(in) :: delta_x  !grid space increment [km]
 real(sdsu_fps), intent(out) :: bw_deg  !beam width [deg]

 real(sdsu_fpd) :: a,b,c,l,dx                    !triangle length [km]
 real(sdsu_fpd) :: alpha, beta, gam, eta, theta  !triangle angle [rad]

!
! in case
!
 if( delta_x == 0. ) then
   bw_deg = 0.
   return
 endif

!
! initialize a,b,alpha into double precision
! 
 call earth_radius( dble(lat_deg*const_degrad), a ) !Earth Core - FOV [km]
 b = dble(sma)                                      !Earth core - Satellite  [km]
 alpha = dble( ona_deg*const_degrad )               !Sensor off-nadir angle  [rad]


 if( alpha == 0. ) then !special case

   c    = b-a
   dx   = dble( delta_x )
   l    = dsqrt( c*c + dx*dx )
   theta = dacos( c / l )
   bw_deg = real(theta) / const_degrad  ! beam width from center pixel to dx-shifted pixel [deg]
   return

 else

   beta = dpi - dasin( b / a * dsin(alpha) )   !caution here ... [rad]
   gam  = dpi - (alpha + beta)           !earth core angle [rad]
   c    = b * dsin( gam ) / dsin( beta ) ! FOV - Satellite [km] 

   eta  = beta - dpi/2.                  !incident angle [rad] 

! print*,'abc',a,b,c
! print*,'angle', alpha/const_degrad, beta/const_degrad, gam/const_degrad , (dpi-beta)/const_degrad, eta/const_degrad

 if( delta_x > 0. ) then
     dx    = dble( delta_x )                            !increment [dx]
      l    = dsqrt( c*c + dx*dx - 2.d0 * c * dx * dcos( eta ) ) ! [km]
     theta = dacos( ( c*c + l*l - dx*dx ) / ( 2.d0 * c * l ) )   ! [rad]
 else
     dx    = dble( -1.*delta_x )                            !increment [dx]
      l    = dsqrt( c*c + dx*dx - 2.d0 * c * dx * dcos( dpi - eta ) ) ! [km]
     theta = dacos( ( c*c + l*l - dx*dx ) / ( 2.d0 * c * l ) )   ! [rad]
 endif

 bw_deg = real(theta) / const_degrad  ! beam width from center pixel to dx-shifted pixel [deg]

 return

 endif
 end subroutine beam_width_dx
 
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


 subroutine gain_function_gaussian(beamwidth3db, alpha, f_gaussian)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments:  For a given 3dB (half power) beam width of specific sensor, this routine give IFOV
!            gain function for various beam angle. Gain function is assumed in Gaussian distribution. 
!
! History:
! 10/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: beamwidth3db  ! sensor half-power beamwidth [rad]
 real(sdsu_fpd),intent(in) :: alpha         ! angle between FOV center and surround grid.  [rad]
 real(sdsu_fpd),intent(out) :: f_gaussian   ! gain function [dB]
 real(sdsu_fpd),parameter :: sqrt2 = 1.4142135381698608   ! sqrt(2.)
 real(sdsu_fpd),parameter :: sqrt2pi = 2.5066282749176025 ! sqrt(2.*acos(-1.0) )
 real(sdsu_fpd),parameter :: bunbo = 2.3548200130462646   ! 2.0*sqrt(2.0*log(2.0) ) 
 real(sdsu_fpd) :: sigma

 if(alpha == -999.) then
   f_gaussian = 0.
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
 end subroutine gain_function_gaussian

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_convolution

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

