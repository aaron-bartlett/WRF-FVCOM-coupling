!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  FileUnit_mod
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Contains utility routines and parameters for using files.
!
! HISTORY:
! 16 Aug 2016 - Forked from FileUnit_mod.  Added PATH_MAX public parameter
!               to specify maximum file path.  Changed ERROR to FATAL in
!               error messages where the program halts. Minor tweak to
!               code for selecting available file unit (now checks if unit
!               exists).
!------------------------------------------------------------------------------

module FileUnit_mod

   ! Change default behavior
   implicit none
   private
   save

   ! Public routines
   public :: select_file_unit
   public :: close_file_unit

   ! Maximum guaranteed path length (in characters) by POSIX standard
   integer, parameter, public :: PATH_MAX = 255
   
contains

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! select_file_unit
   !
   ! DESCRIPTION:
   ! Selects available file unit.
   !---------------------------------------------------------------------------

   function select_file_unit() result (file_unit)

      ! Return variable
      integer :: file_unit

      ! Local variables
      logical :: connected
      integer :: i

      ! Determine available file unit
      file_unit = 0
      connected=.true. ! Assume already in use
      do i = 1,100
         if (i == 5) cycle ! Skip standard input
         if (i == 6) cycle ! Skip standard output
         inquire(unit=i,opened=connected)
         if (.not. connected) then
            file_unit=i
            exit ! Break out of do loop
         end if
      end do
      if (file_unit == 0) then
         write(0,*)'ERROR, cannot find available file unit number!'
         stop
      end if

   end function select_file_unit

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! close_file_unit
   !
   ! DESCRIPTION:
   ! Closes file and handles error.
   !---------------------------------------------------------------------------

   subroutine close_file_unit(file_unit)

      ! Arguments
      integer,intent(in) :: file_unit

      ! Local variables
      integer :: io_status

      ! Close the file
      close(file_unit,iostat=io_status)
      if (io_status /= 0) then
         write(0,*)'ERROR, cannot close file unit ',file_unit
         stop
      end if

   end subroutine close_file_unit

end module FileUnit_mod
