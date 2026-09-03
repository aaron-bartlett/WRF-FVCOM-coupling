!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  SPORT_MODIS_NDVI_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Contains code to read NASA MSFC SPoRT daily MODIS-based NDVI product in 
! GrADS format.
!
! REVISION HISTORY:
! 02 Jun 2014 - Initial version.
! 22 Jul 2016 - Renamed SPORT_NDVI_mod
! 28 Jul 2016 - Added status return flag
! 29 Jul 2016 - Renamed to emphasize MODIS as source of NDVI product.  Bug
!               fix for passing root directory with SPORT data file.  Also 
!               changed logic to expect YYYY and MM subdirectories for SPORT 
!               data files.
! 16 Aug 2016 - Now uses FileUtil_mod to specify maximum length of file path.
!               Replaced ERROR with FATAL in error messages where the program
!               halts.
!
!------------------------------------------------------------------------------

module SPORT_MODIS_NDVI_mod

   ! Defaults
   implicit none
   private

   ! Public method
   public :: readNDVI

   ! SPoRT grid specifications for North American grid
   integer, parameter :: NX = 6300
   integer, parameter :: NY = 2900
   real, parameter :: SW_LAT =   23.00 !   23.00 N
   real, parameter :: SW_LON = -128.00 ! -128.00 E
   real, parameter :: NE_LAT =   52.00 !   52.00 N
   real, parameter :: NE_LON =  -65.00 !  -65.00 E
   real, parameter :: D_LATLON = 0.01  ! Resolution of data
   real, parameter :: NDVI_MIN = 0.05
   real, parameter :: NDVI_MISSING = -9999.

contains
   
   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNDVI
   !
   ! DESCRIPTION:  Routine for reading 2D NDVI field from GrADS file.
   !
   !---------------------------------------------------------------------------

   subroutine readNDVI(topInputDir,year,month,mday, &
        dim1,dim2,ndvi,swlat,swlon,dlatlon,istatus)

      ! Modules
      use FileUnit_mod

      ! Defaults
      implicit none

      ! Arguments
      character(len=*),intent(in) :: topInputDir
      integer,intent(in) :: year
      integer,intent(in) :: month
      integer,intent(in) :: mday
      integer,intent(out) :: dim1
      integer,intent(out) :: dim2
      real,allocatable,intent(out) :: ndvi(:,:)
      real,intent(out) :: swlat
      real,intent(out) :: swlon
      real,intent(out) :: dlatlon
      integer,intent(out) :: istatus

      ! Local variables
      character(len=PATH_MAX) :: filename
      integer :: iunit
      integer :: i,j

      ! Get filename
      call create_sport_modis_filename(topInputDir,year,month,mday,filename)

      ! Open the file
      iunit = select_file_unit()
      print*,'Opening ',trim(filename)
      open(unit=iunit, file=trim(filename), status='old', form='unformatted', &
           access='direct', recl=NX*4)

      ! Allocate memory for the NDVI data
      dim1 = NX
      dim2 = NY
      allocate(ndvi(dim1,dim2))

      ! Read the NDVI field
      do j = 1,NY
         read(unit=iunit,rec=j) ndvi(:,j)
      end do

      ! Close the file.
      call close_file_unit(iunit)

      ! Adjust non-missing NDVI to ensure consistent floor value.
      ! Adjust negative NDVI to missing flag.
      do j = 1,NY
         do i = 1,NX
            if (ndvi(i,j) > 0 .and. ndvi(i,j) < NDVI_MIN) then
               ndvi(i,j) = NDVI_MIN
            else if (ndvi(i,j) < 0) then
               ndvi(i,j) = NDVI_MISSING
            end if
         end do
      end do

      swlat = SW_LAT
      swlon = SW_LON
      dlatlon = D_LATLON 

      istatus = 0
   end subroutine readNDVI

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  create_sport_modis_filename
   !
   ! DESCRIPTION:  Assembles file name for SPoRT MODIS NDVI file based on valid
   ! date.
   !
   !---------------------------------------------------------------------------

   subroutine create_sport_modis_filename(idir,year,month,mday,filename)

      ! Modules
      use FileUnit_mod, only: PATH_MAX

      ! Defaults
      implicit none

      ! Arguments
      character(len=*),intent(in)     :: idir      
      integer,intent(in)              :: year
      integer,intent(in)              :: month
      integer,intent(in)              :: mday
      character(len=PATH_MAX),intent(out)    :: filename

      ! Local variables
      character(len=4)  :: syear
      character(len=2)  :: smonth
      character(len=2)  :: sday
      character(len=37) :: basename
      integer :: sum_lens

      ! Construct subdirectory and base file name
      write(unit=syear, fmt='(i4.4)') year
      write(unit=smonth, fmt='(i2.2)') month
      write(unit=sday, fmt='(i2.2)') mday
      basename = syear//'/'//smonth// &
           '/ndvi1km_lismask_'//syear//smonth//sday//'.bin'
      
      ! Sanity check length of full file path                      
      sum_lens = len_trim(idir) + len_trim(basename)
      if (sum_lens > PATH_MAX) then
         print*,'FATAL, full path of SPORT file is too long!'
         print*,'Shorten the directory name and try again!'
         stop 1
      end if

      ! Create full file path
      filename = trim(idir)//trim(basename)
      
   end subroutine create_sport_modis_filename

end module SPORT_MODIS_NDVI_mod
