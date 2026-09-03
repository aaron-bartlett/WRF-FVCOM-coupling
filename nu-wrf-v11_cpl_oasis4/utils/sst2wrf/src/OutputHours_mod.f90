!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  OutputHours_mod
!
! AUTHOR:
! Eric Kemp, NGIS
!
! DESCRIPTION:
! Provides data type and methods for building list of output hours.
!
!------------------------------------------------------------------------------

module OutputHours_mod

  ! Use external modules
  use FieldSST_RSS_mod, only:  FieldSST_RSS

  ! Force explicit declarations
  implicit none

  ! Force explicit public declarations
  private

  ! Public data type
  public :: OutputHours
  type OutputHours
     integer :: numHours
     integer,allocatable :: hours(:)
  end type OutputHours

  ! Public methods
  public :: createOutputHours
  public :: destroyOutputHours

contains

  !----------------------------------------------------------------------------
  ! ROUTINE:
  ! createOutputHours
  !
  ! DESCRIPTION:
  ! Public "constructor method" for OutputHours data type.
  !----------------------------------------------------------------------------
  
  function createOutputHours(numFakeHours, fakeHours, sstData) result (this)

    ! Force explicit variable declarations
    implicit none

    ! Arguments
    integer, intent(in) :: numFakeHours
    integer :: fakeHours(:)
    type(FieldSST_RSS),intent(in) :: sstData

    ! Return variable
    type(OutputHours) :: this

    ! Local variables
    integer :: count
    integer :: curHour
    integer :: i,j,k
    logical :: unique

    ! Build count of unique valid hours 
    count = 0
    do i = 1,numFakeHours
       curHour = fakeHours(i)
       unique = .true. 
       do j = 1,i
          if (i == j) cycle
          if (curHour == fakeHours(j)) then
             unique = .false.
             exit ! Get out of j do loop
          end if
       end do
       if (unique) count = count + 1
    end do 
    curHour = sstData%hour
    unique = .true.
    do i = 1,numFakeHours
       if (curHour == fakeHours(i)) then
          unique = .false.
          exit ! Get out of i loop
       end if
    end do
    if (unique) count = count + 1
    this%numHours = count

    ! Allocate memory
    allocate(this%hours(count))

    ! Now create list of unique valid hours
    i = 1
    this%hours(i) = sstData%hour
    do j = 1,numFakeHours
       curHour = fakeHours(j)
       unique = .true.
       if (curHour == sstData%hour) then
          unique = .false.
          cycle
       end if
       do k = 1,j
          if (j == k) cycle
          if (curHour == fakeHours(k)) then
             unique = .false.
             exit ! Get out of k loop
          end if
       end do
       if (unique) then
          i = i + 1
          this%hours(i) = curHour
       end if
    end do

    return
  end function createOutputHours

  !----------------------------------------------------------------------------
  ! ROUTINE:
  ! destroyOutputHours
  !
  ! DESCRIPTION:
  ! Public "destructor method" for OutputHours data type.
  !----------------------------------------------------------------------------

  subroutine destroyOutputHours(this)

    ! Force explicit declarations
    implicit none

    ! Arguments
    type(OutputHours),intent(inout) :: this

    ! Clean up
    this%numHours = 0
    deallocate(this%hours)

    return
  end subroutine destroyOutputHours
end module OutputHours_mod
