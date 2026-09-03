!------------------------------------------------------------------------------
! NASA/GSFC, Computational and Information Sciences and Technology Office,
! Code 606
!------------------------------------------------------------------------------
!
! PROGRAM:  ndviBareness4Wrf
!
! AUTHOR:
! Eric Kemp, NASA GSFC/SSAI
!
! DESCRIPTION:
! Top level program for reading NDVI data, deriving surface
! bareness using NDVI thresholds, and writing bareness in WPS intermediate
! binary format for METGRID.
!
! REVISION HISTORY:
! 02 Jun 2014 - Initial version
! 22 Jul 2016 - Revised to add support for GIMMS MODIS NDVI data in GeoTIFF
!               format
! 28 Jul 2016 - Dropped support for copying NDVI data to multiple times
!               (better to just use symbolic links).  GIMMS data now processed
!               by tile due to memory limitations.
! 03 Aug 2016 - Now discriminates between GIMMS Aqua and Terra files.
! 16 Aug 2016 - Now uses FileUtil_mod to specify maximum length of file paths.
!               Replaced ERROR with FATAL in error messages where the program
!               halts.
!
!------------------------------------------------------------------------------

program ndviBareness4Wrf

   ! Modules
   use CalcNdviBareness_mod
   use FieldWPS_mod
   use FileUnit_mod
   use FileWPS_mod
   use NDVI_mod 
   use nml_mod
   
   ! Defaults
   implicit none

   ! Local variables
   integer :: ndviFormat
   real :: lowerleftLat
   real :: lowerleftLon
   real :: upperrightLat
   real :: upperrightLon
   real :: ndviBarenessThreshold1
   real :: ndviBarenessThreshold2
   character(len=PATH_MAX) :: inputDirectory
   integer :: year
   integer :: month
   integer :: dayOfMonth
   integer :: dayOfYear
   character(len=PATH_MAX) :: outputDirectory
   character(len=PATH_MAX) :: outputFilePrefix

   integer :: nx
   integer :: ny
   real,allocatable :: ndvi_2d(:,:)
   character(len=32) :: ndvi_map_source
   character(len=46) :: ndvi_desc
   real :: swlat,swlon,dlatlon
   real,allocatable :: bareness_2d(:,:)
   integer :: hour

   integer :: fileUnit
   type(FileWPS) :: wpsFileData
   type(FieldWPS) :: wpsFieldData
   real :: earth_radius 
   integer :: istatus
   integer :: west_x, south_y, east_x, north_y
   integer :: tile_x,tile_y
   logical :: foundData
   character(len=256) :: filename
   type(nml) :: nml_ndvi
   
   foundData = .false.

   ! Process namelist file.
   filename = 'namelist.ndviBareness4Wrf'
   nml_ndvi = nml_new(filename)
   call nml_ndvi%read("settings", "ndviFormat", &
        ndviFormat)      
   call nml_ndvi%read("settings", "lowerLeftLat", &
        lowerleftLat)      
   call nml_ndvi%read("settings", "lowerLeftLon", &
        lowerleftLon) 
   call nml_ndvi%read("settings", "upperRightLat", &
        upperrightLat)      
   call nml_ndvi%read("settings", "upperRightLon", &
        upperrightLon)      
   call nml_ndvi%read("settings", "ndviBarenessThreshold1", &
        ndviBarenessThreshold1)      
   call nml_ndvi%read("settings", "ndviBarenessThreshold2", &
        ndviBarenessThreshold2)      
   call nml_ndvi%read("settings", "inputDirectory", &
        inputDirectory)      
   call nml_ndvi%read("settings", "year", year)      
   call nml_ndvi%read("settings", "month", month)      
   call nml_ndvi%read("settings", "dayOfMonth", &
        dayOfMonth)      
   call nml_ndvi%read("settings", "dayOfYear", &
        dayOfYear)      
   call nml_ndvi%read("settings", "outputDirectory", &
        outputDirectory)      
   call nml_ndvi%read("settings", "outputFilePrefix", &
        outputFilePrefix)      

   ! Calculate earth radius in km
   earth_radius = 111.111111*180./(4.*atan(1.))

   ! Read the NDVI data.  For very large datasets, we need to process into
   ! tiles.
   call calcTileLimits(ndviFormat, &
        lowerLeftLat,lowerLeftLon,upperRightLat,upperRightLon, &
        south_y,west_x,north_y,east_x)

   ! Loop through tiles (or whole dataset)
   do tile_y = north_y, south_y
      do tile_x = west_x, east_x

         ! Fetch the NDVI data
         call readNDVI(ndviFormat,trim(inputDirectory), &
              tile_x,tile_y, &
              year,month,dayOfMonth,dayOfYear, &
              nx,ny,ndvi_2d,ndvi_map_source,ndvi_desc, &
              swlat,swlon,dlatlon,istatus)
         
         if (istatus .ne. 0) cycle
         foundData = .true.

         ! Derive bareness from the NDVI
         call calcNdviBareness(nx,ny,ndvi_2d,ndviBarenessThreshold1, &
              ndviBarenessThreshold2,bareness_2d)
      
         ! For convenience, assume NDVI and bareness are valid at 00 UTC.
         hour = 0
         
         ! If individual tiles are processed, we must change the WPS file
         ! prefix to reflect the tile numbers.
         call updatePrefixForTile(ndviFormat,tile_x,tile_y,outputFilePrefix)

         ! Create the WPS file
         fileUnit = select_file_unit()
         wpsFileData = createFileWPS(fileUnit=fileUnit, &
              outputDirectory=trim(outputDirectory), &
              prefix=trim(outputFilePrefix), &
              year=year,month=month,day=dayOfMonth,hour=hour)
         
         ! Write the bareness.
         wpsFieldData = createFieldWPS(year=year,month=month, day=dayOfMonth, &
              hour=hour, xfcst=0., &
              map_source='NASA GSFC', &
              field ='BARE_DYN', units='dimensionless', &
              desc = 'Bareness based on NDVI', &
              xlvl = 200100., & ! "Ground level" WPS code
              nx=nx,ny=ny, iproj=IPROJ_LATLON, startloc=SWCORNER, &
              startlat=swlat, startlon=swlon, earth_radius=earth_radius, &
              is_wind_grid_rel=.false. , &
              slab=bareness_2d, deltaLat=dlatlon, deltaLon=dlatlon)
         call writeFileWPS(wpsFileData,wpsFieldData)
         call destroyFieldWPS(wpsFieldData)

         ! Write the NDVI field.
         wpsFieldData = createFieldWPS(year=year,month=month, day=dayOfMonth, &
              hour=hour, xfcst=0., &
              map_source=trim(ndvi_map_source), &
              field ='NDVI', units='dimensionless', &
              desc = trim(ndvi_desc), &
              xlvl = 200100., & ! "Ground level" WPS code
              nx=nx,ny=ny, iproj=IPROJ_LATLON, startloc=SWCORNER, &
              startlat=swlat, startlon=swlon, earth_radius=earth_radius, &
              is_wind_grid_rel=.false. , &
              slab=ndvi_2d, deltaLat=dlatlon, deltaLon=dlatlon)
         call writeFileWPS(wpsFileData,wpsFieldData)
         call destroyFieldWPS(wpsFieldData)
         
         ! Clean up
         call destroyFileWPS(wpsFileData)
         call close_file_unit(fileUnit)  
         deallocate(bareness_2d)
         deallocate(ndvi_2d)      

      end do
   end do

   ! Report error if necessary
   if (.not. foundData) then
      print*,'FATAL, no NDVI files were processed!'
      stop 1
   end if

end program ndviBareness4Wrf
