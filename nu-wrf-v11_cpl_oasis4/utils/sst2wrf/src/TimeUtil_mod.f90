!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  TimeUtil_mod
!
! AUTHOR:
! Eric Kemp, NGIS 
!
! DESCRIPTION:
! Contains date and time routines.
!
! REVISION HISTORY:
! 05 Apr 2010 - Initial version
!
!------------------------------------------------------------------------------

module TimeUtil_mod

  ! Force explicit variable declarations
  implicit none

  ! Force explicit public declarations
  private

  ! Public routines
  public :: calcMonthDayOfMonth

contains

  !----------------------------------------------------------------------------
  ! ROUTINE:
  ! calcMonthDayOfMonth
  !
  ! DESCRIPTION:
  ! Calculates month and day-of-month from year and day-of-year.
  !----------------------------------------------------------------------------

  subroutine calcMonthDayOfMonth(year,dayOfYear,month,dayOfMonth)

    ! Force explicit variable declarations
    implicit none

    ! Arguments
    integer,intent(in) :: year
    integer,intent(in) :: dayOfYear
    integer,intent(out) :: month
    integer,intent(out) :: dayOfMonth

    ! Local variables
    integer :: daysInMonth(12)
    integer :: extraDayForLeapYear
    integer :: sum
    integer :: i
    
    ! Leap year determination
    extraDayForLeapYear = 0
    if (mod(year,4).eq.0) then
       extraDayForLeapYear = 1
       if (mod(year,100).eq.0) then
          extraDayForLeapYear = 0
          if (mod(year,400).eq.0) then
             extraDayForLeapYear = 1
          endif
       endif
    endif
    
    ! Find number of days in each month
    daysInMonth(1)  = 31
    daysInMonth(2)  = 28 + extraDayForLeapYear
    daysInMonth(3)  = 31
    daysInMonth(4)  = 30
    daysInMonth(5)  = 31
    daysInMonth(6)  = 30
    daysInMonth(7)  = 31
    daysInMonth(8)  = 31
    daysInMonth(9)  = 30
    daysInMonth(10) = 31
    daysInMonth(11) = 30
    daysInMonth(12) = 31
       
    ! Find current month
    sum = 0
    month = 0
    do i = 1,12
       sum = sum + daysInMonth(i)
       if (sum >= dayOfYear) then
          month = i
          exit ! Break out of do loop
       end if
    end do
    if (month == 0) then
       print*,'ERROR, could not find month from dayOfYear!'
       print*,'dayOfYear = ',dayOfYear
       stop
    end if

    ! Find current day of month
    dayOfMonth = dayOfYear - (sum - daysInMonth(month))

    return
  end subroutine calcMonthDayOfMonth
end module TimeUtil_mod
