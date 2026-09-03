!------------------------------------------------------------------------------
! NASA/GSFC, Computational and Information Sciences and Technology Office,
! Code 606
!------------------------------------------------------------------------------
!
! MODULE:  NDVI_mod
!
! AUTHOR:
! Eric Kemp, NASA GSFC/SSAI
!
! DESCRIPTION:
! Contains code to read NDVI data.  Several formats are supported, each by
! separate modules called under the hood.
!
! REVISION HISTORY:
! 22 Jul 2016 - Initial version.
! 29 Jul 2016 - Revised dependencies.
! 03 Aug 2016 - Now discriminates between GIMMS Aqua and Terra data. 
! 16 Aug 2016 - Now uses FileUtil_mod to specify maximum length of file paths.
!               Replaced ERROR with FATAL in error messages where the program
!               halts.
!
!------------------------------------------------------------------------------

module NDVI_mod

   ! Defaults
   implicit none
   private

   ! Public methods
   public :: readNDVI
   public :: calcTileLimits
   public :: updatePrefixForTile

   ! Public constants
   integer, parameter, public :: SPORT_MODIS_GRADS = 1
   integer, parameter, public :: GIMMS_MODIS_GEOTIFF_AQUA  = 10
   integer, parameter, public :: GIMMS_MODIS_GEOTIFF_TERRA = 11

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNDVI
   !
   ! DESCRIPTION:  Routine for reading 2D NDVI field from supported file 
   ! format.
   !
   !---------------------------------------------------------------------------

   subroutine readNDVI(iformat,inputDirectory,tile_x,tile_y, &
        year,month,mday,yday,dim1,dim2,ndvi,map_source,desc, &
        swlat,swlon,dlatlon,istatus)
      
      ! Modules
      use GIMMS_MODIS_NDVI_mod, only: GIMMS_AQUA_MODIS
      use GIMMS_MODIS_NDVI_mod, only: GIMMS_TERRA_MODIS
      use GIMMS_MODIS_NDVI_mod, only: GIMMS_MODIS_readNDVI => readNDVI
      use SPORT_MODIS_NDVI_mod, only: SPORT_MODIS_readNDVI => readNDVI

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: iformat
      character(len=*),intent(in) :: inputDirectory
      integer,intent(in) :: tile_x,tile_y
      integer,intent(in) :: year, month, mday, yday
      integer,intent(out) :: dim1
      integer,intent(out) :: dim2
      real,allocatable,intent(out) :: ndvi(:,:)
      character(len=32),intent(out) :: map_source
      character(len=46),intent(out) :: desc
      real,intent(out) :: swlat
      real,intent(out) :: swlon
      real,intent(out) :: dlatlon      
      integer,intent(out) :: istatus

      ! Local variables
      integer :: isatellite

      ! Call appropriate read routine
      if (iformat .eq. SPORT_MODIS_GRADS) then
         map_source = 'NASA SPoRT'
         desc = 'NASA SPoRT Daily MODIS NDVI Composite'
         call SPORT_MODIS_readNDVI(trim(inputDirectory),year,month,mday,&
              dim1,dim2,ndvi,swlat,swlon,dlatlon,istatus)
      else if (iformat .eq. GIMMS_MODIS_GEOTIFF_AQUA) then
         map_source = 'NASA GIMMS'
         desc = 'NASA GIMMS Aqua MODIS NDVI'
         isatellite = GIMMS_AQUA_MODIS
         call GIMMS_MODIS_readNDVI(isatellite,trim(inputDirectory),year,yday, &
              tile_x,tile_y, &
              dim1,dim2,ndvi,swlat,swlon,dlatlon,istatus)
      else if (iformat .eq. GIMMS_MODIS_GEOTIFF_TERRA) then
         map_source = 'NASA GIMMS'
         desc = 'NASA GIMMS Terra MODIS NDVI'
         isatellite = GIMMS_TERRA_MODIS
         call GIMMS_MODIS_readNDVI(isatellite,trim(inputDirectory),year,yday, &
              tile_x,tile_y, &
              dim1,dim2,ndvi,swlat,swlon,dlatlon,istatus)
      else
         print*,'FATAL, invalid format for NDVI data!'
         stop 1
      end if
         
   end subroutine readNDVI

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calcTileLimits
   !
   ! DESCRIPTION:  Calculate range of tiles for specified latitude/longitude
   ! bounds.  Only relevant for GIMMS data.
   !
   !---------------------------------------------------------------------------

   subroutine calcTileLimits(ndviFormat, &
        lowerLeftLat,lowerLeftLon, &
        upperRightLat,upperRightLon, &
        south_y,west_x,north_y,east_x)

      ! Modules
      use GIMMS_MODIS_NDVI_mod, only: &
           GIMMS_MODIS_calcTileLimits => calcTileLimits

      ! Defaults
      implicit none

      ! Arguments
      integer,intent(in) :: ndviFormat
      real,intent(in) :: lowerLeftLat
      real,intent(in) :: lowerLeftLon
      real,intent(in) :: upperRightLat
      real,intent(in) :: upperRightLon
      integer,intent(out) :: south_y
      integer,intent(out) :: west_x
      integer,intent(out) :: north_y
      integer,intent(out) :: east_x

      ! Initialize tile range
      north_y = 1
      south_y = 1
      east_x = 1
      west_x = 1

      ! Currently only GIMMS MODIS data are divided into tiles.
      if (ndviFormat .eq. GIMMS_MODIS_GEOTIFF_AQUA .or. &
          ndviFormat .eq. GIMMS_MODIS_GEOTIFF_TERRA) then
         call GIMMS_MODIS_calcTileLimits(lowerLeftLat,lowerLeftLon, &
              upperRightLat,upperRightLon, &
              south_y,west_x,north_y,east_x)
      end if

   end subroutine calcTileLimits

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  updatePrefixForTile
   !
   ! DESCRIPTION:  Replaces the WPS prefix with one that indicates the tile
   ! number.  Only relevant for GIMMS data.
   !
   !---------------------------------------------------------------------------

   subroutine updatePrefixForTile(ndviFormat,tile_x,tile_y,outputFilePrefix)
      
      ! Modules
      use GIMMS_MODIS_NDVI_mod, only: GIMMS_AQUA_MODIS
      use GIMMS_MODIS_NDVI_mod, only: GIMMS_TERRA_MODIS
      use GIMMS_MODIS_NDVI_mod, only: &
           GIMMS_MODIS_updatePrefixForTile => updatePrefixForTile
      use FileUnit_mod, only: PATH_MAX

      ! Defaults
      implicit none

      ! Arguments
      integer,intent(in) :: ndviFormat
      integer,intent(in) :: tile_x
      integer,intent(in) :: tile_y
      character(len=PATH_MAX),intent(inout) :: outputFilePrefix

      ! Local variables
      integer :: isatellite

      ! Only GIMMS MODIS data needs to change the prefix since these data are
      ! split into tiles.
      if (ndviFormat .eq. GIMMS_MODIS_GEOTIFF_AQUA) then
         isatellite = GIMMS_AQUA_MODIS
         call GIMMS_MODIS_updatePrefixForTile(isatellite,tile_x,tile_y, &
              outputFilePrefix)
      else if (ndviFormat .eq. GIMMS_MODIS_GEOTIFF_TERRA) then
         isatellite = GIMMS_TERRA_MODIS
         call GIMMS_MODIS_updatePrefixForTile(isatellite,tile_x,tile_y, &
              outputFilePrefix)
      end if

   end subroutine updatePrefixForTile

end module NDVI_mod
