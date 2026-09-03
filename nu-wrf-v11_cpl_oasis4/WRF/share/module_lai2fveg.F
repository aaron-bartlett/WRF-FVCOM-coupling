!NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF 
!NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF 

 module module_lai2fveg
 implicit none
!--------------------------------------------------------------------------------------------
!             -    LAI to Vegetation Fraction     -
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
! Comments: 
!    This module and subroutine convert LAI to vegetation fraction through Beer's law. 
!
! History:
!  03/2014  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 save
 private
 public :: lai2fveg 

 contains

!NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF 
!NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF 

 subroutine lai2fveg( lai, fveg )
 implicit none
 real,intent(in) :: lai    ! leaf area index [0~7]
 real,intent(out) :: fveg  ! vegetation fraction [0~1]

 real,parameter :: k = 0.5  ! k reflect geometric influence associated wih the 
                            ! angle between leaves and the sun.
                            ! 0.5 is used for planar leaves, which makes theoretical limit of 
                            ! plane-paralell LAI about 6. 

 if ( lai < 0.0 ) then
    print*,'MSG lai2fveg: input LAI is negative',lai
    print*,'MSG lai2fveg: output fveg is forced zero.'
    fveg = 0.e0
    return
 elseif ( isnan(lai) ) then
    print*,'MSG lai2fveg: input LAI is NaN',lai
    print*,'MSG lai2fveg: output fveg is forced zero.'
    fveg = 0.e0
    return 
 elseif ( 1.0e+10/lai == 0.e0 ) then
    print*,'MSG lai2fveg: input LAI is infinity',lai
    print*,'MSG lai2fveg: output fveg is forced one.'
    fveg = 1.e0
    return
 else

    !
    ! Compute vegetation fraction through Lambert-Beer's Law
    ! , which assume plane-paralell canopy. 
    ! So transmission is used as bare fraction.
    !
    fveg = 1.0 - exp ( -1.0 * k * lai ) 

 endif

 if( isnan(fveg) ) then
   print*,'MSG lai2fveg: computed fveg become NaN --> Input LAI is wrong',LAI 
   print*,'forced to zero'
   fveg = 0.0
 elseif( fveg < 0.0 .or. fveg > 1.e0 ) then
   print*,'MSG lai2fveg: computed fveg is out of range [0~1] ',fveg
   print*,'forced to zero'
   fveg = 0.0
 endif


 return
 end subroutine lai2fveg

!NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF 
!NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF 

 end module  module_lai2fveg

!NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF 
!NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF 


