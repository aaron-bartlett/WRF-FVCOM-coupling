!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE: CalcNdviBareness_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Contains code for deriving surface bareness from Normalized Difference
! Vegetation Index (NDVI) field using two thresholds.
!
! REVISION HISTORY:
! 02 Jun 2014 - Initial version
! 16 Aug 2016 - Replaced ERROR with FATAL in error messages where the program
!               halts.
! 18 Aug 2016 - Introduced missing value for bareness for obviously invalid
!               NDVI.
!
!------------------------------------------------------------------------------

module CalcNdviBareness_mod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: calcNdviBareness

   ! Constants
   integer, parameter :: MISSING_BARENESS = -9999.

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calcNdviBareness
   !
   ! DESCRIPTION: Derives surface bareness from NDVI field.
   !
   !---------------------------------------------------------------------------

   subroutine calcNdviBareness(nx,ny,ndvi_2d,ndviBarenessThreshold1, &
        ndviBarenessThreshold2,bareness_2d)

      ! Defaults
      implicit none

      ! Arguments
      integer,intent(in) :: nx
      integer,intent(in) :: ny
      real,intent(in) :: ndvi_2d(nx,ny)
      real,intent(in) :: ndviBarenessThreshold1
      real,intent(in) :: ndviBarenessThreshold2
      real,allocatable,intent(out) :: bareness_2d(:,:)

      ! Local variables
      real :: minThreshold,maxThreshold
      integer :: i,j

      ! Sanity check NDVI thresholds before finding max and min.
      if (ndviBarenessThreshold1 < -0.2.or. &
          ndviBarenessThreshold1 > 1) then
         print*,'FATAL, bad value for ndviBarenessThreshold1'
         print*,'Must be between -0.2 and 1'
         print*,'Found ',ndviBarenessThreshold1
         stop 1
      end if
      if (ndviBarenessThreshold2 < -0.2 .or. &
          ndviBarenessThreshold2 > 1) then
         print*,'FATAL, bad value for ndviBarenessThreshold2'
         print*,'Must be between -0.2 and 1'
         print*,'Found ',ndviBarenessThreshold2
         stop 1
      end if
      if (ndviBarenessThreshold1 .eq. ndviBarenessThreshold2) then
         print*,'FATAL, two different thresholds needed for NDVI bareness!'
         print*,'ndviBarenessThreshold1 = ',ndviBarenessThreshold1
         print*,'ndviBarenessThreshold2 = ',ndviBarenessThreshold2
      end if
      minThreshold = min(ndviBarenessThreshold1,ndviBarenessThreshold2)
      maxThreshold = max(ndviBarenessThreshold1,ndviBarenessThreshold2)

      allocate(bareness_2d(nx,ny))

      ! Soil marked as "bare" if it falls between two NDVI thresholds.
      ! NDVI below lower threshold is probably water, while NDVI above
      ! higher threshold has vegetation.  EXCEPTION:  Obviously low 
      ! NDVI values indicate missing or bad data, and the bareness should
      ! be flagged accordingly.
      do j = 1,ny
         do i = 1,nx
            if (ndvi_2d(i,j) > minThreshold .and. &
                ndvi_2d(i,j) < maxThreshold) then
               bareness_2d(i,j) = 1 
            else if (ndvi_2d(i,j) < -9990.) then
               bareness_2d(i,j) = MISSING_BARENESS
            else
               bareness_2d(i,j) = 0
            end if
         end do
      end do

   end subroutine calcNdviBareness
end module CalcNdviBareness_mod
