#include <define_MPI.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_simulator_micro
 use module_simulator
 use module_opt_micro
 use module_surface
 use module_scan
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
!   Microwave radiometer simulator that compute microwave brightness temperature. 
!
! History:
!  09/2011  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 private   ! ALL variables and subourtines are not accessible in module_radar

!            except...
 public simulator_micro_scan


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine simulator_micro_scan
 implicit none

!--------------------------------------------------------------------------------------------
!
! Comments:  This routine interpolate optical properties along slant path.  
!
! History:
!  09/2011  Toshi Matsui@NASA GSFC ;  Initial. 
!
! References:
!-----------------------------------------------------------------------------------------------------

 integer :: i,j,k


!
! DOMAIN_3D_LOOP
!
 do j = myj_start, myj_end ; do i = myi_start, myi_end

!
! conver double precision (in module_scan) to single precision. 
!
    call get_slant_path( i,j, real( surface_CRM(i,j)%vector_grid  ), &
                              real( surface_CRM(i,j)%vector_sat   ), &
                              real( surface_CRM(i,j)%alpha        ), &
                              real( surface_CRM(i,j)%re           ), &
                              real( surface_CRM(i,j)%dist_sat_grid) &
                                  )

 enddo ; enddo

stop
 return
 end subroutine simulator_micro_scan

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_slant_path(i_ref, j_ref, v_grid, v_sat , alpha, re, c_orig )
 implicit none
 integer,intent(in) :: i_ref, j_ref ! i,j point at reference ellisoid
 real,intent(in) :: v_grid(1:3)  ! grid reference vector position [km]
 real,intent(in) :: v_sat(1:3)   ! satellite vector position [km]
 real,intent(in) :: alpha        ! off-nadir angle
 real,intent(in) :: re          ! distance from earth core to reference ellipsoid (under satellite position)
 real,intent(in) :: c_orig          ! distance from satellite to grid

 integer :: k
 

 do k = 1, mxlyr_micro
   print*,k
    call get_latlon_nextpoint( k, v_grid, v_sat , alpha, re , c_orig)

 enddo  
stop

 return
 end subroutine get_slant_path

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_latlon_nextpoint(k, v_grid, v_sat , alpha, re , c_orig )
 implicit none 
 integer,intent(in) :: k         !vertical level index       
 real,intent(in) :: v_grid(1:3)  ! grid reference vector position [km]
 real,intent(in) :: v_sat(1:3)   ! satellite vector position [km]
 real,intent(in) :: alpha        ! off-nadir angle [rad]
 real,intent(in) :: re          ! distance from earth core to reference ellipsoid (under satellite position)
 real,intent(in) :: c_orig       ! distance from satellite to grid [km]
 real :: a,b,c
 real :: beta, beta_dash, alpha_dash, gamma_dash, gamma_orig

! make triangle

  call get_distance(v_grid, v_sat, c) 
  b = sqrt( sum( v_sat**2 ) )                         ![km]
  a = sqrt( sum( v_grid**2) )

!  a = re + range_micro * real(k) + 0.5 * range_micro  ![km]
! c is unknown

 gamma_orig = acos( (b*b+a*a-c_orig*c_orig)/(2.*b*a) )
 beta = const_pi - gamma_orig  - alpha

print*, alpha/const_degrad, beta/const_degrad, gamma_orig/const_degrad

print*,a/sin(alpha), b/sin(beta), c / sin(gamma_orig)
stop

print*,2./sin(90.*const_degrad),1./sin(30.*const_degrad), sqrt(3.)/sin(60.*const_degrad)
stop

  alpha_dash = const_pi - beta

  beta_dash = asin( b * sin(alpha) / re ) 

  gamma_dash = const_pi - alpha_dash - beta_dash


print*,alpha/const_degrad, beta/const_degrad, alpha_dash/const_degrad, beta_dash/const_degrad, gamma_dash/const_degrad
stop

!
! get c
!


 return
 end subroutine get_latlon_nextpoint

 subroutine get_distance( v1, v2, dist)
 implicit none
 real(sdsu_fps),intent(in) :: v1(3)  !3D vector 1
 real(sdsu_fps),intent(in) :: v2(3)  !3D vector 2
 real(sdsu_fps),intent(out) :: dist

!
! distance [unit is same as v1 and v2 
!
 dist = sqrt( (v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 + (v1(3)-v2(3))**2 )

 return
 end subroutine get_distance


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine within_latlon ( lat_fov, lon_fov, lat_cen, lon_cen, logic_within_fov )
 implicit none
 real,intent(in) :: lat_fov(1:4) , lon_fov(1:4)  !lat lon of FOV 4 corners. 
 real,intent(in) :: lat_cen, lon_cen             !lat lon of center fov

 logical,intent(out) :: logic_within_fov

 real :: mlon_fov(1:4), mlon_cen  ![deg]

 real :: high_lat, low_lat, left_lon, right_lon  ![deg]

 logic_within_fov = .false.

 if( maxval(lon_fov(1:4)) > 100.e0 .and. minval(lon_fov(1:4)) < -100.e0 ) then !pixcel crossing lon=-180 

    mlon_fov(1:4) = lon_fov(1:4)
    if( lon_fov(1) < 0.e0 ) mlon_fov(1) = lon_fov(1) + 360.e0
    if( lon_fov(2) < 0.e0 ) mlon_fov(2) = lon_fov(2) + 360.e0
    if( lon_fov(3) < 0.e0 ) mlon_fov(3) = lon_fov(3) + 360.e0
    if( lon_fov(4) < 0.e0 ) mlon_fov(4) = lon_fov(4) + 360.e0

    mlon_cen = lon_cen 
    if( lon_cen < 0.e0 ) mlon_cen = lon_cen + 360.e0
    !print*,'SPECIAL CASE'

 else  !normal case 

    mlon_fov(1:4) = lon_fov(1:4)
    mlon_cen = lon_cen

 endif

 high_lat = maxval( lat_fov(1:4) )   ! [deg]
 low_lat  = minval( lat_fov(1:4) )   ! [deg]

 left_lon  = minval( mlon_fov(1:4) ) ! [deg]
 right_lon = maxval( mlon_fov(1:4) ) ! [deg] 

!
! 1st test....
!
 if(  lat_cen >  low_lat .and.  lat_cen <  high_lat .and. &
     mlon_cen > left_lon .and. mlon_cen < right_lon ) then

    ! 2nd test ...
    call within_fov( lat_fov, mlon_fov, lat_cen, mlon_cen, logic_within_fov )    

 else

     return  

 endif 


 return
 end subroutine within_latlon

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine within_fov ( lat_fov, lon_fov, lat_cen, lon_cen, logic_within_fov )
 implicit none
 real,intent(in) :: lat_fov(1:4) , lon_fov(1:4)  ! lat lon of FOV 4 corners. 
 real,intent(in) :: lat_cen, lon_cen             ! lat lon of center fov
 logical,intent(out) :: logic_within_fov         ! logic

 integer :: i
 integer :: i_1, i_2, i_3, i_4

 real :: y, slope, intercept

 logic_within_fov = .true.

!
! determin grid point here (as long as rectangular this algorithm works to 
! find the pixel within OMI corner or not. 
!  
!  |-------------------2----------|
!  |                   |          |
!  |      case 1       |  case 2  |
!  |       |-----------|----------3 
!  1-------|-----------|          |
!  | case4 |        case 3        |
!  |       |                      |
!  |-------4----------------------|

 i_1 = 0; i_2 = 0; i_3=0; i_4=0

 do i = 1, 4
   if( lon_fov(i) == minval(lon_fov(:)) )then
      i_1 = i ; exit
   endif
 enddo

 do i = 1, 4
   if( lat_fov(i) == maxval(lat_fov(:)) ) then
      i_2 = i ; exit
   endif 
 enddo

 do i = 1, 4
   if( lon_fov(i) == maxval(lon_fov(:)) )then
      i_3 = i ; exit
   endif
 enddo

 do i = 1, 4
   if( lat_fov(i) == minval(lat_fov(:)) )then
      i_4 = i ; exit
   endif
 enddo

 if( lat_cen > lat_fov(i_1) .and. lat_cen < lat_fov(i_2) .and. &
    lon_cen > lon_fov(i_1) .and. lon_cen < lon_fov(i_2) ) then     !case 1

    slope = (lat_fov(i_2)-lat_fov(i_1))/(lon_fov(i_2)-lon_fov(i_1))
    intercept = lat_fov(i_1) 
    y = slope*(lon_cen-lon_fov(i_1) ) + intercept
    if( lat_cen >= y ) logic_within_fov = .false.
!print*,'case 1 cal lat',lat_cen, 'y=', y,slope*(lon_cen-lon_fov(i_1) ), intercept
 endif

 if( lat_cen > lat_fov(i_3) .and. lat_cen < lat_fov(i_2) .and. &
     lon_cen > lon_fov(i_2) .and. lon_cen < lon_fov(i_3) ) then   !case 2

    slope = (lat_fov(i_3)-lat_fov(i_2))/(lon_fov(i_3)-lon_fov(i_2))
    intercept = lat_fov(i_2) 
    y = slope*(lon_cen-lon_fov(i_2) ) + intercept
    if( lat_cen >= y ) logic_within_fov = .false.
!print*,'case 2 cal lat',lat_cen, 'y=', y,slope*(lon_cen-lon_fov(i_1) ), intercept
 endif

 if( lat_cen > lat_fov(i_4) .and. lat_cen < lat_fov(i_3) .and. &
        lon_cen > lon_fov(i_4) .and. lon_cen < lon_fov(i_3) ) then   !case 3

    slope = (lat_fov(i_3)-lat_fov(i_4))/(lon_fov(i_3)-lon_fov(i_4))
    intercept = lat_fov(i_4)
    y = slope*(lon_cen-lon_fov(i_4) ) + intercept
    if( lat_cen <= y ) logic_within_fov = .false.
!print*,'case 3 cal lat',lat_cen, 'y=', y,slope*(lon_cen-lon_fov(i_1) ), intercept
 endif

 if( lat_cen > lat_fov(i_4) .and. lat_cen < lat_fov(i_1) .and. &
        lon_cen > lon_fov(i_1) .and. lon_cen < lon_fov(i_4) ) then   !case 4

    slope = (lat_fov(i_4)-lat_fov(i_1))/(lon_fov(i_4)-lon_fov(i_1))
    intercept = lat_fov(i_1)
    y = slope*(lon_cen-lon_fov(i_1) ) + intercept
    if( lat_cen <= y ) logic_within_fov = .false.
!print*,'case 4 cal lat',lat_cen, 'y=', y,slope*(lon_cen-lon_fov(i_1) ), intercept
 endif


 return
 end subroutine within_fov

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


 end module module_simulator_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
