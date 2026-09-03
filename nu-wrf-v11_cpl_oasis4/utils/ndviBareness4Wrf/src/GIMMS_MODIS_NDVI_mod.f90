!------------------------------------------------------------------------------
! NASA/GSFC, Computational and Information Sciences and Technology Office,
! Code 606
!------------------------------------------------------------------------------
!
! MODULE: GIMMS_MODIS_NDVI_mod
!
! AUTHOR:
! Eric Kemp, NASA GSFC/SSAI
!
! DESCRIPTION:
! Contains code to read NASA GIMMS MODIS NDVI data in GeoTIFF format.
!
! REVISION HISTORY:
! 22 July 2016 - Initial version.  Borrows heavily from code from LIS 
!                developers.
! 28 July 2016 - Changed code to only process and return a single tile.
!                This is due to memory constraints that prevents METGRID
!                from processing a mosaic of NDVI.
! 29 Jul 2016 - Renamed to emphasize MODIS as source of data.
! 03 Aug 2016 - Now discriminates between Aqua and Terra.
! 16 Aug 2016 - Now uses FileUtil_mod to specify maximum length of file paths.
!               Replaced ERROR with FATAL in error messages where the program
!               halts.  Removed redundant addition of "/" at end of directory.
!
!------------------------------------------------------------------------------

module GIMMS_MODIS_NDVI_mod
   
   ! Defaults
   implicit none
   private

   ! Public methods
   public :: readNDVI
   public :: calcTileLimits
   public :: updatePrefixForTile

   ! Public constants
   integer, parameter, public :: GIMMS_AQUA_MODIS = 1
   integer, parameter, public :: GIMMS_TERRA_MODIS = 2

   ! GIMMS regional tile size
   integer, parameter :: TILE_NC = 4000
   integer, parameter :: TILE_NR = 4000

   real, parameter :: NDVI_MISSING = -9999.

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNDVI
   !
   ! DESCRIPTION:  Routine for reading 2D NDVI field from GIMMS MODIS GeoTIFF 
   ! files.
   !
   ! This borrows heavily from LIS program read_gimms_modisndvi; however,
   ! only a single tile is processed.
   !
   !---------------------------------------------------------------------------

   subroutine readNDVI(isatellite,topInputDir,year,yday,tile_x,tile_y, &
        dim1,dim2,ndvi,sw_lat,sw_lon,dlatlon,istatus)

      ! Modules
      use, intrinsic:: iso_c_binding ! Standard Fortran 2003
      use FileUnit_mod, only: PATH_MAX
      use fortranc
      use gdal

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: isatellite
      character(len=*),intent(in) :: topInputDir
      integer,intent(in) :: year
      integer,intent(in) :: yday
      integer,intent(in) :: tile_x
      integer,intent(in) :: tile_y
      integer,intent(out) :: dim1
      integer,intent(out) :: dim2
      real,allocatable,intent(out) :: ndvi(:,:)
      real,intent(out) :: sw_lat, sw_lon
      real,intent(out) :: dlatlon
      integer,intent(out) :: istatus

      ! Local variables
      real, allocatable :: ndvi_tile_uint(:,:)
      real, allocatable :: ndvi_tile(:,:)
      character(len=PATH_MAX) :: filename
      logical :: file_exists
      real :: dres
      integer :: i,r,c,j

      ! GDAL related local variables
      type(gdaldriverh)      :: driver
      type(gdaldataseth)     :: ds
      integer(kind=c_int)    :: ierr
      real(kind=c_double)    :: gt(6)
      integer(kind=c_int)    :: xsize,ysize
      real(kind=c_double)    :: x1, y1, x2, y2
      type(gdalrasterbandh)  :: band

      ! Initialize error flag
      istatus = 1 ! Assume error

      ! Initialize GDAL library
      call GDALAllRegister()

      ! Allocate and initialize array
      allocate(ndvi_tile(TILE_NC, TILE_NR))
      ndvi_tile = NDVI_MISSING

      ! What is the filename for the current tile?
      call create_gimms_modis_filename(isatellite, topInputDir,&
           year, yday, tile_x, tile_y, filename )

      ! Make sure the file exists
      inquire(file=filename,exist=file_exists)
      if (.not. file_exists) then
         print*,'WARNING, ',trim(filename),' does not exist!'
         return
      end if

      ! Use GDAL routines to open the TIFF file
      driver = gdalgetdriverbyname('Tif'//CHAR(0))
      ds = gdalopen( trim(filename)//CHAR(0), GA_ReadOnly)
      if (.not. gdalassociated(ds)) then
         print*,'FATAL, cannot open dataset in ',trim(filename)
         stop 1
      end if

      ! Retrieve coefficients for transforming pixel/line raster space
      ! and projection coordinate (Xp,Yp) space
      ierr = gdalgetgeotransform(ds, gt)
      dres = abs(gt(6))

      ! Get x- and y-resolutions of TIFF file being read.
      xsize = gdalgetrasterxsize(ds)
      ysize = gdalgetrasterysize(ds)
      
      ! Convert pixel,line tuple to georeferenced location
      call gdalapplygeotransform(gt, 0.5_c_double, 0.5_c_double, x1, y1)
      y2 = y1 - (ysize-1)*dres ! Lower left latitude

      ! Get NDVI band/layer (layer=1)
      band = gdalgetrasterband(ds, 1)
      if (.not. gdalassociated(band)) then
         print*,'FATAL, cannot get NDVI band from file ',trim(filename)
         stop 1
      end if

      ! Get NDVI layer input.  Data are unsigned 8-bit integers
      allocate(ndvi_tile_uint(TILE_NC, TILE_NR))
      ierr = gdalrasterio_f( band, GF_Read, 0, 0, ndvi_tile_uint)
      if (ierr .ne. 0) then
         print*,'FATAL, problem reading data from file ',trim(filename)
         stop 1
      end if

      ! Convert data to regular floating point.  Note the following:
      ! Valid range of NDVI is 0 to 250.
      ! 251 and 252 are reserved mask values.
      ! 253 indicates invalid land (out of range NDVI)
      ! 254 indicates water.
      ! 255 indicates no data (unfilled, cloudy, or snow contaminated)
      do j = 1, TILE_NR
         do i = 1, TILE_NC
            if (ndvi_tile_uint(i,j) .le. 250) then
               ndvi_tile(i,j) = ndvi_tile_uint(i,j) * 0.004
            else
               ndvi_tile(i,j) = NDVI_MISSING
            end if
         end do
      end do
      deallocate(ndvi_tile_uint)

      ! Close file
      call gdalclose(ds)

      ! Now reverse the Y-rows.  GIMMS NDVI data start from N and go S, but
      ! WRF wants to start S and go N.
      dim1 = TILE_NC
      dim2 = TILE_NR
      allocate(ndvi(dim1,dim2))
      ndvi = NDVI_MISSING
      i = 0
      do r = dim2,1,-1
         i = i + 1
         do c = 1, dim1
            ndvi(c,i) = ndvi_tile(c,r)
         end do
      end do

      ! Copy map information
      sw_lat = y2
      sw_lon = x1
      dlatlon = dres
      
      ! Clean up
      deallocate(ndvi_tile)

      istatus = 0 ! No error

   end subroutine readNDVI

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calcTileLimits
   !
   ! DESCRIPTION:  Calculates range of tiles that contain NDVI data for 
   ! specified latitude/longitude bounds.
   !
   !---------------------------------------------------------------------------

   subroutine calcTileLimits(lowerLeftLat,lowerLeftLon, &
        upperRightLat,upperRightLon, &
        south_y,west_x,north_y,east_x)

      ! Defaults
      implicit none

      ! Arguments
      real,intent(in) :: lowerLeftLat
      real,intent(in) :: lowerLeftLon
      real,intent(in) :: upperRightLat
      real,intent(in) :: upperRightLon
      integer,intent(out) :: south_y
      integer,intent(out) :: west_x
      integer,intent(out) :: north_y
      integer,intent(out) :: east_x

      ! Calculate tile limits
      west_x  = floor((180 + lowerLeftLon ) / 9)
      south_y = floor(( 90 - lowerLeftLat ) / 9)
      east_x  = floor((180 + upperRightLon) / 9)
      north_y = floor(( 90 - upperRightLat) / 9)

   end subroutine calcTileLimits

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  create_gimms_modis_filename
   !
   ! DESCRIPTION:  Routine for assembling GIMMS GeoTIFF tile filename.  
   !
   ! This borrows heavily from LIS program read_gimms_modisndvi, developed by
   ! Code 617.
   !
   !---------------------------------------------------------------------------

   subroutine create_gimms_modis_filename(isatellite,idir,yr,doy,x,y,filename)

      ! Modules
      use FileUnit_mod, only: PATH_MAX

      ! Defaults
      implicit none

      ! Arguments
      integer,intent(in)              :: isatellite
      character(len=*),intent(in)     :: idir      
      integer,intent(in)              :: yr
      integer,intent(in)              :: doy
      integer,intent(in)              :: x
      integer,intent(in)              :: y
      character(len=PATH_MAX),intent(out)    :: filename

      ! Local variables
      character(len=4)          :: fyr
      character(len=3)          :: fdoy
      character(len=2)          :: fx, fy
      character(len=58)         :: basename
      integer :: sum_lens

      ! Construct subdirectory and base file name
      write(unit=fyr, fmt='(i4.4)')  yr
      write(unit=fdoy, fmt='(i3.3)') doy 
      write(unit=fx, fmt='(i2.2)') x
      write(unit=fy, fmt='(i2.2)') y
      if (isatellite .eq. GIMMS_AQUA_MODIS) then
         basename = fyr//'/'//fdoy//&
              '/GMYD09Q1.A'//fyr//fdoy//'.08d.latlon.x'//fx//'y'//fy//&
              '.6v1.NDVI.tif'
      else if (isatellite .eq. GIMMS_TERRA_MODIS) then
         basename = fyr//'/'//fdoy//&
              '/GMOD09Q1.A'//fyr//fdoy//'.08d.latlon.x'//fx//'y'//fy//&
              '.6v1.NDVI.tif'
      else
         print*,'FATAL, invalid satellite selected for GIMMS MODIS NDVI!'
         stop 1
      end if

      ! Sanity check length of full file path
      sum_lens = len_trim(idir) + len_trim(basename)
      if (sum_lens > PATH_MAX) then
         print*,'FATAL, full path of GIMMS file is too long!'
         print*,'Shorten the directory name and try again!'
         stop 1
      end if

      ! Create full file path
      filename = trim(idir)//trim(basename)
         
   end subroutine create_gimms_modis_filename

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  updatePrefixForTile
   !
   ! DESCRIPTION:  Replaces the WPS file prefix with one that indicate the 
   ! GIMMS tile numbers.
   !
   !---------------------------------------------------------------------------

   subroutine updatePrefixForTile(isatellite,tile_x,tile_y,outputFilePrefix)

      ! Modules
      use FileUnit_mod, only: PATH_MAX
      
      ! Defaults
      implicit none

      ! Arguments
      integer,intent(in) :: isatellite
      integer,intent(in) :: tile_x
      integer,intent(in) :: tile_y
      character(len=PATH_MAX),intent(out) :: outputFilePrefix

      ! Replace the prefix
      if (isatellite .eq. GIMMS_AQUA_MODIS) then
         write(outputFilePrefix,'(A,I2.2,A,I2.2)') &
              'AQUA_X',tile_x,'Y',tile_y
      else if (isatellite .eq. GIMMS_TERRA_MODIS) then
         write(outputFilePrefix,'(A,I2.2,A,I2.2)') &
              'TERRA_X',tile_x,'Y',tile_y
      else
         print*,'FATAL, invalid satellite for MODIS NDVI!'
         stop 1
      end if
   end subroutine updatePrefixForTile

end module GIMMS_MODIS_NDVI_mod
