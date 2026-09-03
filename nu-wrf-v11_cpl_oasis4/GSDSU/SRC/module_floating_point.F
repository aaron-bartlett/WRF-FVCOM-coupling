#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_floating_point
 implicit none

!--------------------------------------------------------------------------------------------------
! Comments:                             
!      This module specify floating point in the G-SDSU. 
!
!  History:
!  03/2012  Toshi Matsui@NASA GSFC : Init 
!
! References:
!-----------------------------------------------------------------------------------------------------

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 public   ! ALL variables and subourtines are accessible in module_simulator.


!
! SDSU floating point control
! 
 integer,parameter :: sdsu_fps = SELECTED_REAL_KIND( 6, 37)  !single precision precision kind
 integer,parameter :: sdsu_fpd = SELECTED_REAL_KIND(15,307)  !double precision precision kind


 end module module_floating_point

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
