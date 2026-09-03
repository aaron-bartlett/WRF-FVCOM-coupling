!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  MODIS_NDVI_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Contains code to read NASA MODIS daily NDVI product in GrADS format.
!
! REVISION HISTORY:
! 2 June 2014 - Initial version.
!------------------------------------------------------------------------------

module MODIS_NDVI_mod

   ! New defaults
   implicit none
   private

   ! Public method
   public :: readNDVI
   public :: readNDVI_test

   ! Public constants
   public :: SW_LAT, SW_LON
   public :: D_LATLON

   ! MODIS grid specifications for global grid
   !integer, parameter :: NX = 43200
   !integer, parameter :: NY = 21600
   !real, parameter :: SW_LAT =  -89.99583 !
   !real, parameter :: SW_LON = -179.99583 !
   !real, parameter :: NE_LAT =   89.9960  !
   !real, parameter :: NE_LON =  179.952   !
   !real, parameter :: D_LATLON = 0.00833333  ! Resolution of data
   !real, parameter :: NDVI_MIN = -0.2        ! 0.05
   !real, parameter :: NDVI_MISSING = -9999.

   ! MODIS grid specifications for ACMAP-asia (eas) grid
   integer, parameter :: NX = 15604          !(42006-26403)
   integer, parameter :: NY = 9600           !(18000-8401)
   real, parameter :: SW_LAT = -19.9955
   real, parameter :: SW_LON =  40.0002
   real, parameter :: NE_LAT =  60.0015
   real, parameter :: NE_LON = 170.0040 
   real, parameter :: D_LATLON = 0.00833333  ! Resolution of data
   real, parameter :: NDVI_MIN = -0.2        ! 0.05
   real, parameter :: NDVI_MISSING = -9999.

  ! MODIS grid specifications for North American grid
   !integer, parameter :: NX = 7562
   !integer, parameter :: NY = 3481
   !real, parameter :: SW_LAT =  23.0040
   !real, parameter :: SW_LON = -127.992
   !real, parameter :: NE_LAT =  52.0081
   !real, parameter :: NE_LON = -64.9989
   !real, parameter :: D_LATLON = 0.00833333  ! Resolution of data
   !real, parameter :: NDVI_MIN = 0.00        ! 0.05
   !real, parameter :: NDVI_MISSING = -9999.

contains
   
   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNDVI
   !
   ! DESCRIPTION:  Routine for reading 2D NDVI field from GrADS file.
   !
   !---------------------------------------------------------------------------

   subroutine readNDVI(filename,dim1,dim2,ndvi)

      ! Import modules
      use FileUnit_mod

      implicit none

      ! Arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: dim1
      integer,intent(out) :: dim2
      !real,allocatable,intent(out) :: ndvi(:,:)
      real*4,allocatable,intent(out) :: ndvi(:,:)

      ! Local variables
      integer :: iunit
      integer :: i,j
      
      ! Open the file
      iunit = select_file_unit()
      print*,'Opening ',trim(filename)
      open(unit=iunit, file=trim(filename), status='old', form='unformatted')

      ! Allocate memory for the NDVI data
      dim1 = NX
      dim2 = NY
      allocate(ndvi(dim1,dim2))

      ! Read the NDVI field
      read(unit=iunit) ndvi(:,:)
   
      !print*,'test(2000,1000:1100)'
      !print*,ndvi(2000,1000:1100)

      ! Close the file.
      call close_file_unit(iunit)

      ! Adjust non-missing NDVI to ensure consistent floor value.
      ! Adjust negative NDVI to missing flag.
      do j = 1,NY
         do i = 1,NX
            if (ndvi(i,j) > 0 .and. ndvi(i,j) < NDVI_MIN) then
               !dkim6, comment out
               !ndvi(i,j) = NDVI_MIN
            else if (ndvi(i,j) < 0) then
               ndvi(i,j) = NDVI_MISSING
            else if (ndvi(i,j) > 1) then  !also exclude large NDVI > 1
               ndvi(i,j) = NDVI_MISSING
            end if
         end do
      end do

      return
   end subroutine readNDVI

   subroutine readNDVI_test
    real*4 :: ndvi1(7562,3481)  ! -- NAM 
    character(len=100) :: fname
    ndvi1(:,:)=1.0
    fname='/discover/nobackup/dkim6/WRFDATA/satellite/modis_vi/00ncl/nam_MOD13A1_20110712_30s.bin'
    open(unit=1,file=trim(fname),status='old',form='unformatted')
    read(1)ndvi1
    close(1)
    print*,'ndvi1'
    print*,maxval(ndvi1),minval(ndvi1)
    print*,ndvi1(2000,1000:1100)
   end subroutine readNDVI_test

end module MODIS_NDVI_mod
