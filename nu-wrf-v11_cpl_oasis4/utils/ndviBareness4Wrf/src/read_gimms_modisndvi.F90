program read_gimms_modisndvi

  use, INTRINSIC:: iso_c_binding
  use fortranc
  use gdal

  implicit none
!
! DATASET
!
!        Description     GIMMS MODIS Terra & Aqua NDVI 8-bit 9x9 degree Tiles
!        Coverage        Global 180W to 180E, 80N to 60S
!        Projection      Lat/Lon grid
!        Spatial Res.    0.00225 x 0.00225 degrees (~250-meter)
!        Temporal Res.   8-day composite period
!        File Format     TIFF
!        Download        http://gimms.gsfc.nasa.gov/download/MODIS/
!                        ftp://gimms.gsfc.nasa.gov/MODIS/
!
!  NEAR REAL-TIME AND STANDARD PROCESSING
!
!        Forward near real-time (nrt) processing of this dataset is typically
!        available within 12 hours from the last day of the 8-day composite
!        period.  Science quality, standard (std) processing is typically
!        available one or more days after the last day of the 8-day composite
!        period. The two processing modes differ in their upstream data source.
!
!        Processing          Data Source
!        ----------          -----------
!        nrt                 LANCE near real-time
!        std                 MODAPS science quality
!
!        * GIMMS nrt datasets are only available for 30 days after production.
!
!   Qualifier           Value       Description
!   ---------           -----       -----------
!   Dataset             GMOD09Q1    GIMMS MODIS Terra (MOD=Terra MYD=Aqua)
!   Start date          A2010001    Starting year 2010, day-of-year 001
!   Composite period    08d         8-day composite
!   Projection          latlon      Lat/Lon grid
!   9x9 Tile index      x00y02      Column 00, Row 02
!   Versions            5v3         MODAPS collection 5, GIMMS version 3
!   Layer name          NDVI        Normalized Vegetation Index
!   File format         tif         Tagged Image File
!
! PROJECTION
!
!   This dataset is mapped to a global 180W to 180E, 90N to 90S Lat/Lon grid
!    and divided into 9x9 degree tiles.
!
!   Parameter           Value       Units
!   ---------           -----       -----
!   Datum               WGS84       -
!   Upper left Lat      90.0        deg
!   Upper left Lon      -180.0      deg
!   Pixel size          0.00225     deg
!   Grid x size         160000      pixel
!   Grid y size         80000       pixel
!
!   Coordinates are mapped and measured from the upper left corner of the
!     grid cell.
!
!   Since the size of a global grid file may be too large for users,
!     the grid is divided into 40 columns and 20 rows (starting at the
!     upper left corner at 180W 90N) to create 9x9 degree tiles of
!     4000 x 4000 pixels.
!
!     Parameter                            Value
!     ---------                            -----
!     Tile size             9 / 0.00225 =  4000
!     No. columns (x)     160000 / 4000 =  40
!     No. rows (y)         80000 / 4000 =  20
!
! GIMMS-MODIS NDVI tile files:
!- To calculate the upper left corner {Lon,Lat} for tile {x,y}:
!   UL-Lon = -180 + (x * 9)
!   UL-Lat =   90 - (y * 9)
!
! ____________________________________________________________________

  integer       :: year
  integer       :: doy
  character*100 :: odir
  character*140 :: filename
  logical       :: file_exists
  integer       :: tile_x, west_x, east_x, x_cnt
  integer       :: tile_y, south_y, north_y, y_cnt
  integer       :: Xdir_numtiles, Ydir_numtiles
  real          :: lat, lon
  real          :: ll_lat, ll_lon
  real          :: ur_lat, ur_lon
  real          :: dres

  integer       :: nc, nr
  real          :: cornerlat1, cornerlat2
  real          :: cornerlon1, cornerlon2
  real          :: gridDesc(20) 

  integer, parameter :: tile_nc=4000, tile_nr=4000
  integer       :: col_cnt1, col_cnt2, row_cnt1, row_cnt2
  integer       :: ire

!#if (defined USE_GDAL)
  TYPE(gdaldriverh)      :: driver
  TYPE(gdaldataseth)     :: ds
  TYPE(gdalrasterbandh)  :: band
  INTEGER(kind=c_int)    :: numbands
  INTEGER(kind=c_int)    :: xsize, ysize
  REAL(kind=c_double)    :: x1, y1, x2, y2, gt(6)
  INTEGER(kind=c_int)    :: i1, j1, k1, i, j, k, ierr, c,r

  real, allocatable      :: read_ndvi_uint(:,:)   ! for unsigned integer*8
  real, allocatable      :: read_ndvi(:,:,:)   
  real, allocatable      :: mosaic_ndvi(:,:)
  real, allocatable      :: yrev_ndvi(:,:)

!#endif

! __________________________________________________


  write(*,*) " Reading in GIMMS MODIS NDVI *.tif files ... "

!#if (defined USE_GDAL)
  call GDALAllRegister()

! Set output directory and initialize filename:
  odir = "../NDVI/"
  filename = "none"

! Initialize year and day-of-year "doy":
  year = 2015
  doy  = 001

! Regional lat/lon box selected by user:
! Example run domain:
!  Map projection of the LIS domain:    latlon
!  Run domain lower left lat:         -11.750 #south
!  Run domain lower left lon:          22.050 #west
!  Run domain upper right lat:         22.950 #north
!  Run domain upper right lon:         51.350
!  Run domain resolution (dx):          0.10
!  Run domain resolution (dy):          0.10

   ll_lon = 22.050
   ll_lat = -11.750 
   ur_lon = 51.350
   ur_lat = 22.950 
!
! - USE THIS APPROACH FOR SELECTING THE TILED FILES TO OPEN -
! - To calculate the {x,y} tile index for a given {Lon,Lat}:
!   x = floor((180 + Lon) / 9)
!   y = floor(( 90 - Lat) / 9)
!
   west_x  = floor((180 + ll_lon) / 9)
   south_y = floor(( 90 - ll_lat) / 9)
   east_x  = floor((180 + ur_lon) / 9)
   north_y = floor(( 90 - ur_lat) / 9)
   Xdir_numtiles = (east_x-west_x)+1
   Ydir_numtiles = (south_y-north_y)+1

!   print *, ll_lon, ur_lon, ur_lat, ll_lat
!   print *, west_x, east_x, north_y, south_y
!   print *, " west_x, east_x, north_y, south_y"
!   print *, " Number of Xdir_numtiles :: ", Xdir_numtiles
!   print *, " Number of Ydir_numtiles :: ", Ydir_numtiles

   allocate( read_ndvi(tile_nc, tile_nr, Xdir_numtiles) )
   allocate( mosaic_ndvi(Xdir_numtiles*tile_nc, Ydir_numtiles*tile_nr) )
   read_ndvi = -9999.
   mosaic_ndvi = -9999.

   ! Derive full mosaicked GIMMS tiles domain corner points == 
   !  used for setting the parameter array for reprojection/interp
   !  to LVT eval or LIS runnin domain:
   ! Initialize lat/lon points:
   cornerlat1 =  9999. ! [LL lat of mosaicked tile domain]
   cornerlat2 = -9999. ! [UR lat of mosaicked tile domain]
   cornerlon1 =  9999. ! [LL lon of mosaicked tile domain]
   cornerlon2 = -9999. ! [UR lon of mosaicked tile domain]

!- Loop over relevant tiles and mosaic them together:
   y_cnt = 0
   do tile_y = north_y, south_y
      y_cnt = y_cnt + 1
      row_cnt1 = (tile_y-north_y)*tile_nr+1
      row_cnt2 = row_cnt1+tile_nr-1

      x_cnt = 0
      do tile_x = west_x, east_x
         x_cnt = x_cnt + 1 
         col_cnt1 = (tile_x-west_x)*tile_nc+1
         col_cnt2 = col_cnt1+tile_nc-1
         print *, tile_y, row_cnt1, row_cnt2, tile_x, col_cnt1, col_cnt2

       ! Put together GIMMS MODIS filename to be read in:
         call create_gimms_modis_filename( odir,&
              year, doy, tile_x, tile_y, filename )

         inquire( file=filename, exist=file_exists )
         if( file_exists ) then
            write(*,*) "[INFO] Reading GIMMS MODIS file ... ",trim(filename)
         else
            write(*,*) " ## GIMMS MODIS file DOES NOT exist ... ",trim(filename)
         endif

       ! Use GDAL routines to open the tiff files:
         driver = gdalgetdriverbyname('Tif'//CHAR(0))
         ds = gdalopen( trim(filename)//CHAR(0), GA_ReadOnly)

         if( .not.gdalassociated(ds) ) then
            write(*,*) "[ERR] opening dataset on file ",trim(filename)," failed ..."
!            call LVT_endrun()
            stop
         end if

       ! (GDAL) Retrieves the coefficients for transforming pixel/line raster
       !   space and projection coordinate space (Xp,Yp) space:

         ierr = gdalgetgeotransform(ds, gt)

         print *, "UL_long = gt(1): ",gt(1)
         print *, "Xres(+) = gt(2): ",gt(2)
         print *, "UL_Lat: = gt(4): ",gt(4)
         print *, "Yres(-) = gt(6): ",gt(6)
         dres = abs(gt(6))
         print *, 'Resolution = ', dres

        ! Note: Upper-left most geographic point (Xp1,Yp1):
        ! UL_long = gt(1) = -171.000000000000
        ! xres    = gt(2) =    2.250000000000000E-003  (+ direction from W to E)
        ! UL_lat  = gt(4) =   72.0000000000000
        ! yres    = gt(6) =   -2.250000000000000E-003  (- direction from N to S)

       ! (GDAL) Check for x- and y-resolutions of *Tif file being read:
         xsize = gdalgetrasterxsize(ds)
         ysize = gdalgetrasterysize(ds)
         print *, "Number of xpts, ypts: ",xsize, ysize, ierr

       ! (GDAL) Note: GDALapplygeotransform is used to convert
       !   a (pixel,line) coordinate into a georeferenced (geo_x,geo_y) location:

         CALL gdalapplygeotransform(gt, 0.5_c_double, 0.5_c_double, x1, y1)

         x2 = (xsize-1)*dres + x1
         y2 = y1 - (ysize-1)*dres

         print *, 'UL_lat(y1) = ', y1
         print *, 'LR_lat(y2) = ', y2
         print *, 'UL_lon(x1) = ', x1
         print *, 'LR_lon(x2) = ', x2

         cornerlat1 = min(cornerlat1, y2)  ! [LL lat of mosaicked tile domain]
         cornerlat2 = max(cornerlat2, y1)  ! [UR lat of mosaicked tile domain]
         cornerlon1 = min(cornerlon1, x1)  ! [LL lon of mosaicked tile domain]
         cornerlon2 = max(cornerlon2, x2)  ! [UR lon of mosaicked tile domain]

      ! Check number of NDVI file layers (or bands):
!        numbands = GDALGetRasterCount(ds)

      ! Retrieve NDVI band/layer (layer=1) --
        band = gdalgetrasterband(ds, 1)

      ! Check if NDVI band is present:
        if (.NOT.gdalassociated(band)) THEN
           write(*,*) "[ERR] failed getting NDVI band from TIFF dataset in file, ",&
                      TRIM(filename)
           stop
!           call LVT_endrun()
        endif

      ! Read in NDVI layer input:
        allocate( read_ndvi_uint(tile_nc, tile_nr) )
        ierr = gdalrasterio_f( band, GF_Read, 0, 0, read_ndvi_uint )
        if (ierr /= 0) then
           write(*,*) '[ERR] reading data from TIFF dataset on file ',TRIM(filename)
           stop
!           call LVT_endrun()
        endif

     !- Convert data format (UINT8) to regular floating point::      
        do j = 1, tile_nr
           do i = 1, tile_nc
           !- Convert the UINT8 values to reals::
              if ( read_ndvi_uint(i,j) .le. 250 ) then
                 read_ndvi(i,j,x_cnt) = read_ndvi_uint(i,j) * 0.004
              else
                 read_ndvi(i,j,x_cnt) = read_ndvi_uint(i,j) 
              endif
           end do
         end do
         deallocate( read_ndvi_uint )

      !- Mosaic NDVI tiles together:
         mosaic_ndvi(col_cnt1:col_cnt2, row_cnt1:row_cnt2) = read_ndvi(:,:,x_cnt)

         call gdalclose(ds)
      enddo
   enddo
   deallocate( read_ndvi )

!#ENDIF

 ! Reverse-Y rows of read-in input file:
   allocate( yrev_ndvi(Xdir_numtiles*tile_nc, Ydir_numtiles*tile_nr) )
   yrev_ndvi = -9999.
   i = 0
   do r = Ydir_numtiles*tile_nr, 1, -1
      i = i + 1
      do c = 1, Xdir_numtiles*tile_nc
         yrev_ndvi(c,i) = mosaic_ndvi(c,r)
      end do
   end do
   deallocate( mosaic_ndvi )

   print *, "Mosaicked domain Lat range: ",cornerlat1, cornerlat2
   print *, "Mosaicked domain Lon range: ",cornerlon1, cornerlon2

   ! Filling the items needed by the interpolation library
     nr = nint((cornerlat2 - cornerlat1)/dres)+1
     nc = nint((cornerlon2 - cornerlon1)/dres)+1

     gridDesc(:) = 0
     gridDesc(1) = 0 
     gridDesc(2) = nc
     gridDesc(3) = nr
     gridDesc(4) = cornerlat1
     gridDesc(5) = cornerlon1
     gridDesc(7) = cornerlat2
     gridDesc(8) = cornerlon2
     gridDesc(9) = dres
     gridDesc(10) = dres
     gridDesc(6) = 128
     gridDesc(20) = 64

   deallocate( yrev_ndvi) 

   write(*,*) " -- Completed reading in the GIMMS MODIS NDVI files -- "

end program read_gimms_modisndvi


!BOP
! 
! !ROUTINE: create_gimms_modis_filename
! \label{create_gimms_modis_filename}
!
! !INTERFACE: 

subroutine create_gimms_modis_filename(odir,yr,doy,x,y,filename)
! 
! !USES:   
  implicit none
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for GIMMS_MODIS data files 
! based on the given date (year, doy) and x,y tile.
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GIMMS base directory
!   \item[yr]        Year of data
!   \item[doy]       Day of year 
!   \item[x]         X-direction tile index
!   \item[y]         Y-direction tile index
!   \item[filename]  Name of the GIMMS_MODIS file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
! !ARGUMENTS: 
  character(len=*)     :: odir
!  integer              :: use_anomaly
  integer              :: yr
  integer              :: doy
  integer              :: x
  integer              :: y
  character(len=*)     :: filename

  character*4          :: fyr
  character*3          :: fdoy
  character*2          :: fx, fy

  write(unit=fyr, fmt='(i4.4)')  yr
  write(unit=fdoy, fmt='(i3.3)') doy 
  write(unit=fx, fmt='(i2.2)') x
  write(unit=fy, fmt='(i2.2)') y


! ./GIMMS_MODIS/NDVI/2015/001/GMOD09Q1.A2015001.08d.latlon.x21y09.5v3.NDVI.tif

  filename = trim(odir)//'/'//fyr//'/'//fdoy//&
             '/GMOD09Q1.A'//fyr//fdoy//'.08d.latlon.x'//fx//'y'//fy//'.5v3.NDVI.tif'


end subroutine create_gimms_modis_filename
