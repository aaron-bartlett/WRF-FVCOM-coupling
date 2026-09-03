!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  GEOS2WPS_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! Defines subroutines for reading GEOS variables and writing out in WPS
! format.
!
! REVISION:
! 16 Apr 2012 - First version
!  3 Nov 2014 - Added subsetting option.
!
!------------------------------------------------------------------------------

module GEOS2WPS_mod

   ! Import modules
   use FieldWPS_mod
   use FileUnit_mod
   use FileUtils_mod
   use FileWPS_mod
   use NamelistGEOS2WPS_mod

   ! Reset defaults
   implicit none
   private

   ! Public routines
   public :: driver

   ! Private data type to simplify internal subroutine calls
   type internalData
      integer :: fileID
      integer :: longitudeDimension
      integer :: latitudeDimension
      real :: deltaLongitude
      real :: deltaLatitude
      real :: southwestLongitude
      real :: southwestLatitude
      real :: variableIndex
      real :: earthRadius
   end type internalData

contains
   
   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  driver
   !
   ! DESCRIPTION:  High-level driver for reading GEOS variables and writing
   ! to WPS format.
   !
   !--------------------------------------------------------------------------

   subroutine driver(namelist)

      ! Arguments
      type(namelistGEOS2WPS), intent(in) :: namelist

      ! Local variables
      double precision, allocatable :: longitudes(:)
      double precision, allocatable :: latitudes(:)
      type(internalData) :: internal
      integer :: iv

      ! Calculate earth radius in km
      internal%earthRadius = 111.111111*180./(4.*atan(1.)) 

      ! Open the GEOS file
      internal%fileID = openReadFile(namelist%geosFileFormat, &
           trim(namelist%geosFileName))

      ! Read the longitudes
      if (namelist%subset) then
         call readDoubleArray1d(namelist%geosFileFormat, &
              internal%fileID, trim(namelist%longitudeName), &
              internal%longitudeDimension, longitudes, &
              namelist%subset, &
              namelist%iLonMin, &
              namelist%iLonMax)
      else
         call readDoubleArray1d(namelist%geosFileFormat, &
              internal%fileID, trim(namelist%longitudeName), &
              internal%longitudeDimension, longitudes)
      end if
      internal%deltaLongitude = longitudes(2) - longitudes(1)
      internal%southwestLongitude = longitudes(1)
      deallocate(longitudes)

      ! Read the latitudes
      if (namelist%subset) then
         call readDoubleArray1d(namelist%geosFileFormat, &
              internal%fileID, trim(namelist%latitudeName), &
              internal%latitudeDimension, latitudes, &
              namelist%subset, &
              namelist%jLatMin, &
              namelist%jLatMax)
      else
         call readDoubleArray1d(namelist%geosFileFormat, &
              internal%fileID, trim(namelist%latitudeName), &
              internal%latitudeDimension, latitudes)
      end if
      internal%deltaLatitude = latitudes(2) - latitudes(1)
      internal%southwestLatitude = latitudes(1)
      deallocate(latitudes)
      
      ! Loop through the requested variables
      do iv = 1, namelist%numberOfVariables
         if (namelist%variableRanks(iv) == 3) then
            internal%variableIndex=iv
            call process3dGeosVar(namelist,internal)            
         else if (namelist%variableRanks(iv) == 4) then            
            internal%variableIndex=iv
            call process4dGeosVar(namelist,internal)
         end if
      end do

   end subroutine driver

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  process3dGeosVar
   !
   ! DESCRIPTION:  Read 3d array from GEOS HDF4 or netCDF file and outputs
   ! slices to WPS intermediate format file.
   !
   !--------------------------------------------------------------------------

   subroutine process3dGeosVar(namelist,internal)
      
      ! Arguments
      type(NamelistGEOS2WPS), intent(in) :: namelist
      type(InternalData), intent(in) :: internal

      ! Local variables
      type(FileWPS) :: wpsFileData
      integer :: timeDimension
      integer :: dim1, dim2, dim3
      integer :: year,month,day,hour,minute,second
      real, allocatable :: tmp3d(:,:,:)
      character(len=132) :: prefix
      real :: xlvl
      integer :: iv,in

      iv = internal%variableIndex

      ! Select xlvl value for WPS slab.
      select case (namelist%variableLevelTypes(iv))
      case (GROUND_LEVEL, MEAN_SEA_LEVEL, TWO_METERS_ABOVE_GROUND_LEVEL, &
           TEN_METERS_ABOVE_GROUND_LEVEL)
         xlvl = selectXlvl(namelist%variableLevelTypes(iv))
      case default
         print*,'ERROR, invalid level type for 3D GEOS array!'
         print*,'Value is ',namelist%variableLevelTypes(iv)
         print*,'Valid values:'
         print*,'    Ground level:  ',GROUND_LEVEL
         print*,'    2-meters AGL: ',TWO_METERS_ABOVE_GROUND_LEVEL
         print*,'    10-meters AGL: ',TEN_METERS_ABOVE_GROUND_LEVEL
         print*,'    Mean sea level: ',MEAN_SEA_LEVEL
         stop 1
      end select

      ! Write file prefix
      select case (namelist%variableLevelTypes(iv))
      case (GROUND_LEVEL)
         prefix=trim(namelist%variableNamesOut(iv))//'_GROUND_LEVEL'
      case (TWO_METERS_ABOVE_GROUND_LEVEL)
         prefix=trim(namelist%variableNamesOut(iv))//&
              '_2M_ABOVE_GROUND_LEVEL'
      case (TEN_METERS_ABOVE_GROUND_LEVEL)
         prefix=trim(namelist%variableNamesOut(iv))//&
              '_10M_ABOVE_GROUND_LEVEL'
      case (MEAN_SEA_LEVEL)
         prefix=trim(namelist%variableNamesOut(iv))//&
              '_MEAN_SEA_LEVEL'
      end select

      ! If we are subsetting the data, read and process each time
      ! slice separately.
      if (namelist%subset) then
         do in = 1, namelist%numberOfTimes
            call readRealArray3d(namelist%geosFileFormat,internal%fileID, &
                 trim(namelist%variableNamesIn(iv)),dim1,dim2,dim3,tmp3d, &
                 namelist%subset, &
                 namelist%iLonMin,namelist%jLatMin,&
                 namelist%timeIndices(in), &
                 namelist%iLonMax,namelist%jLatMax, &
                 namelist%timeIndices(in))
            if (dim1 .ne. internal%longitudeDimension .or. &
                dim2 .ne. internal%latitudeDimension) then
               print*,'ERROR, dimension mismatch for ', &
                    trim(namelist%variableNamesIn(iv))
               print*,'longitude dimension was ',dim1,', expected ', &
                    internal%longitudeDimension
               print*,'latitude dimension was ',dim2,', expected ', &
                    internal%latitudeDimension
               stop 1
            end if
            call openSlabFile(namelist,prefix,in,wpsFileData, &
                 year,month,day,hour,minute,second)
            call writeSlab(namelist,iv,in,xlvl,internal, &
                 tmp3d(1,1,1), &
                 internal%longitudeDimension, &
                 internal%latitudeDimension, & 
                 wpsFileData, &
                 year,month,day,hour,minute,second)
            call closeSlabFile(wpsFileData)

            deallocate(tmp3d)
         end do
         
      else
         ! If processing entire grid, read all at once and then process
         ! each time slice.
         call readRealArray3d(namelist%geosFileFormat,internal%fileID, &
              trim(namelist%variableNamesIn(iv)),dim1,dim2,dim3,tmp3d)
         if (dim1 .ne. internal%longitudeDimension .or. &
             dim2 .ne. internal%latitudeDimension) then
            print*,'ERROR, dimension mismatch for ', &
                 trim(namelist%variableNamesIn(iv))
            print*,'longitude dimension was ',dim1,', expected ', &
                 internal%longitudeDimension
            print*,'latitude dimension was ',dim2,', expected ', &
                 internal%latitudeDimension
            stop 1
         end if
         timeDimension = dim3
         do in = 1, namelist%numberOfTimes
            if (namelist%timeIndices(in) < 1 .or. &
                namelist%timeIndices(in) > timeDimension) then
               print*,'ERROR, invalid time index requested!'
               print*,'Requested time = ',in
               print*,'Valid range is 1 to ',timeDimension
               stop 1
            end if
            call openSlabFile(namelist,prefix,in,wpsFileData, &
                 year,month,day,hour,minute,second)
            call writeSlab(namelist,iv,in,xlvl,internal, &
                 tmp3d(1,1,namelist%timeIndices(in)), &
                 internal%longitudeDimension, &
                 internal%latitudeDimension, & 
                 wpsFileData, &
                 year,month,day,hour,minute,second)
            call closeSlabFile(wpsFileData)

         end do
         deallocate(tmp3d)
      end if

      return

   end subroutine process3dGeosVar

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  process4dGeosVar
   !
   ! DESCRIPTION:  Read 4d array from GEOS HDF4 or netCDF file and outputs
   ! slices to WPS intermediate format file.
   !
   !--------------------------------------------------------------------------

   subroutine process4dGeosVar(namelist,internal)
      
      ! Arguments
      type(NamelistGEOS2WPS), intent(in) :: namelist
      type(InternalData), intent(in) :: internal

      ! Local variables
      double precision, allocatable :: vertical(:)
      type(FileWPS) :: wpsFileData
      integer :: timeDimension
      integer :: verticalDimension
      integer :: dim1, dim2, dim3, dim4
      integer :: year,month,day,hour,minute,second
      real, allocatable :: tmp4d(:,:,:,:)
      real, allocatable :: tmp1d(:)
      logical :: needToFlip
      character(len=132) :: prefix
      real :: xlvl
      integer :: iv,in
      integer :: k,kk

      iv = internal%variableIndex

      ! Read the vertical coordinate
      if ( .not. namelist%hasVerticalDimension) then
         print*,'ERROR, trying to read 4D array'
         print*,'But namelist set hasVerticalDimension to false!'
         stop 1
      end if
      if (namelist%subset) then
         call readDoubleArray1d(namelist%geosFileFormat,internal%fileID, &
              trim(namelist%verticalName),verticalDimension,vertical, &
              namelist%subset, &
              namelist%kVertMin, &
              namelist%kVertMax)
      else
         call readDoubleArray1d(namelist%geosFileFormat,internal%fileID, &
              trim(namelist%verticalName),verticalDimension,vertical)
      end if
      needToFlip=.false. ! First guess
      if (vertical(2) - vertical(1) < 0 .and. &
           namelist%variableLevelTypes(iv) == MODEL_LEVEL) then
         needToFlip=.true.
         allocate(tmp1d(verticalDimension))
         do k = 1, verticalDimension
            tmp1d(k) = vertical(verticalDimension - k + 1)
         end do
         do k = 1, verticalDimension
            vertical(k) = tmp1d(k)
         end do
         deallocate(tmp1d)
      end if

      if (.not. namelist%variableLevelTypes(iv) .eq. MODEL_LEVEL .and. &
          .not. namelist%variableLevelTypes(iv) .eq. ISOBARIC_LEVEL) then
         print*,"ERROR, invalid level type for rank 4 data!"
         print*,'Selected level type is ', &
              namelist%variableLevelTypes(iv)
         print*,'Valid options: '
         print*,'    Model level: ',MODEL_LEVEL
         print*,'    Pressure level: ',ISOBARIC_LEVEL
         stop 1
      end if
      
      ! Write file prefix
      select case (namelist%variableLevelTypes(iv))
      case (MODEL_LEVEL)
         prefix=trim(namelist%variableNamesOut(iv))//'_MODEL_LEVEL'
      case (ISOBARIC_LEVEL)
         prefix=trim(namelist%variableNamesOut(iv))//'_ISOBARIC_LEVEL'
      end select

      if (namelist%subset) then

         ! For subset, read and process each time slice separately.
         do in = 1, namelist%numberOfTimes

            call readRealArray4d(namelist%geosFileFormat,internal%fileID, &
                 trim(namelist%variableNamesIn(iv)),dim1,dim2,dim3,dim4, &
                 tmp4d, &
                 namelist%subset, &
                 namelist%iLonMin,namelist%jLatMin,namelist%kVertMin, &
                 namelist%timeIndices(in), &
                 namelist%iLonMax,namelist%jLatMax,namelist%kVertMax, &
                 namelist%timeIndices(in))
            if (dim1 .ne. internal%longitudeDimension .or. &
                dim2 .ne. internal%latitudeDimension .or. &
                dim3 .ne. verticalDimension) then
               print*,'ERROR, dimension mismatch for ', &
                    trim(namelist%variableNamesIn(iv))
               print*,'longitude dimension was ',dim1,', expected ', &
                    internal%longitudeDimension
               print*,'latitude dimension was ',dim2,', expected ', &
                    internal%latitudeDimension
               print*,'vertical dimension was ',dim3,', expected ', &
                    verticalDimension
               stop 1
            end if

            call openSlabFile(namelist,prefix,in,wpsFileData, &
                 year,month,day,hour,minute,second)
            
            do k = 1, verticalDimension
               if (needToFlip) then
                  kk = verticalDimension - k + 1
               else
                  kk = k
               end if

               ! Select xlvl value for WPS slab.
               ! FIXME: Move to new subroutine.
               select case (namelist%variableLevelTypes(iv))
               case (MODEL_LEVEL)
                  ! Need to figure out which GEOS model level was used
                  xlvl = selectXlvl(namelist%variableLevelTypes(iv), &
                       level=real(namelist%kVertMin + k - 1))
               case (ISOBARIC_LEVEL)
                  xlvl = selectXlvl(namelist%variableLevelTypes(iv), &
                       level=real(vertical(k)*100.))
               case default
                  print*,'ERROR, invalid level type for 4D GEOS array!'
                  print*,'Value is ',namelist%variableLevelTypes(iv)
                  print*,'Valid values:'
                  print*,'    Model level: ',MODEL_LEVEL
                  print*,'    Isobaric level: ',ISOBARIC_LEVEL
                  stop 1
               end select

               call writeSlab(namelist,iv,in,xlvl,internal, &
                    tmp4d(1,1,kk,1), &
                    internal%longitudeDimension, &
                    internal%latitudeDimension, &
                    wpsFileData, &
                    year,month,day,hour,minute,second)
            end do ! End vertical loop

            call closeSlabFile(wpsFileData)

         end do ! End time loop
      else
         ! If processing whole grid, read entire array first and then
         ! process each time slice.
         call readRealArray4d(namelist%geosFileFormat,internal%fileID, &
              trim(namelist%variableNamesIn(iv)),dim1,dim2,dim3,dim4,tmp4d)
         
         if (dim1 .ne. internal%longitudeDimension .or. &
             dim2 .ne. internal%latitudeDimension .or. &
             dim3 .ne. verticalDimension) then
            print*,'ERROR, dimension mismatch for ', &
                 trim(namelist%variableNamesIn(iv))
            print*,'longitude dimension was ',dim1,', expected ', &
                 internal%longitudeDimension
            print*,'latitude dimension was ',dim2,', expected ', &
                 internal%latitudeDimension
            print*,'vertical dimension was ',dim3,', expected ', &
                 verticalDimension
            stop 1
         end if
         
         timeDimension = dim4

         do in = 1, namelist%numberOfTimes
            if (namelist%timeIndices(in) < 1 .or. &
                namelist%timeIndices(in) > timeDimension) then
               print*,'ERROR, invalid time index requested!'
               print*,'Requested time = ',in
               print*,'Valid range is 1 to ',timeDimension
               stop 1
            end if

            call openSlabFile(namelist,prefix,in,wpsFileData, &
                 year,month,day,hour,minute,second)
                           
            do k = 1, verticalDimension
               if (needToFlip) then
                  kk = verticalDimension - k + 1
               else
                  kk = k
               end if

               ! Select xlvl value for WPS slab.
               ! FIXME: Put in subroutine
               select case (namelist%variableLevelTypes(iv))
               case (MODEL_LEVEL)
                  ! Need to figure out which GEOS model level was used
                  xlvl = selectXlvl(namelist%variableLevelTypes(iv), &
                       level=real(kk))
               case (ISOBARIC_LEVEL)
                  xlvl = selectXlvl(namelist%variableLevelTypes(iv), &
                       level=real(vertical(kk)*100.))
               case default
                  print*,'ERROR, invalid level type for 4D GEOS array!'
                  print*,'Value is ',namelist%variableLevelTypes(iv)
                  print*,'Valid values:'
                  print*,'    Model level: ',MODEL_LEVEL
                  print*,'    Isobaric level: ',ISOBARIC_LEVEL
                  stop 1
               end select
               
               call writeSlab(namelist,iv,in,xlvl,internal, &
                    tmp4d(1,1,kk,namelist%timeIndices(in)), &
                    internal%longitudeDimension, &
                    internal%latitudeDimension, &
                    wpsFileData, &
                    year,month,day,hour,minute,second)

            end do ! Vertical loop

            call closeSlabFile(wpsFileData)

         end do ! Time loop

      end if ! subset?
         
      deallocate(vertical)
      return
   end subroutine process4dGeosVar

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  selectXlvl
   !
   ! DESCRIPTION:  Creates xlvl variable specifying level of 2D slab written
   ! to WPS intermediate format.
   !
   !--------------------------------------------------------------------------

   function selectXlvl(variableLevelType,level) result (xlvl)

      ! Arguments
      integer, intent(in) :: variableLevelType
      real, intent(in), optional :: level ! Pressure or model level

      ! Return variable
      real :: xlvl

      select case (variableLevelType)
      case (GROUND_LEVEL)
         xlvl = 200100.
      case (MEAN_SEA_LEVEL)
         xlvl = 201300.
      case (TWO_METERS_ABOVE_GROUND_LEVEL)
         xlvl = 200100.
      case (TEN_METERS_ABOVE_GROUND_LEVEL)
         xlvl = 200100.
      case (MODEL_LEVEL, ISOBARIC_LEVEL)
         if (.not. present(level)) then
            print*,'ERROR, missing level value!'
            stop 1
         end if
         xlvl = level
         if (xlvl == 200100. .or. &
             xlvl == 201300.) then
            print*,'ERROR, input level inconsistent with level type!'
            print*,'Input level is ',level
            if (xlvl == 200100.) then
               print*,'WPS reserves this for ground or near ground level.'
            else if (xlvl == 201300.) then
               print*,'WPS reserves this for mean sea level.'
            end if
            stop 1
         end if
      case default
         print*,'Invalid variable level type!'
         print*,'Value is ',variableLevelType
         print*,'Valid values:'
         print*,'    Ground level:  ',GROUND_LEVEL
         print*,'    2-meters AGL: ',TWO_METERS_ABOVE_GROUND_LEVEL
         print*,'    10-meters AGL: ',TEN_METERS_ABOVE_GROUND_LEVEL
         print*,'    Mean sea level: ',MEAN_SEA_LEVEL
         print*,'    Model level: ',MODEL_LEVEL
         print*,'    Pressure level: ',ISOBARIC_LEVEL
         stop 1
      end select

      return
   end function selectXlvl

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  openSlabFile
   !
   ! DESCRIPTION:  Opens WPS file for 2d slabs.
   !
   !--------------------------------------------------------------------------

   subroutine openSlabFile(namelist,prefix,in,wpsFileData, &
        year,month,day,hour,minute,second)

      implicit none

      ! Arguments
      type(namelistGEOS2WPS), intent(in) :: namelist
      character(len=*),intent(in) :: prefix
      integer,intent(in) :: in
      type(FileWPS),intent(out) :: wpsFileData
      integer,intent(out) :: year,month,day,hour,minute,second

      ! Local variables
      integer :: fileUnit
      logical :: includeMinute,includeSecond

      ! Write to the WPS file
      fileUnit=select_file_unit()
      includeMinute=.false.
      includeSecond=.false.
      minute = 0
      second = 0
      if (len_trim(namelist%validTimes(in)) .eq. 13) then
         read(namelist%validTimes(in),'(i4.4,1x,i2.2,1x,i2.2,1x,i2.2)') &
              year,month,day,hour
      else if (len_trim(namelist%validTimes(in)) .eq. 16) then
         read(namelist%validTimes(in), &
              '(i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2)') &
              year,month,day,hour,minute
         includeMinute=.true.
      else if (len_trim(namelist%validTimes(in)) .eq. 19) then
         read(namelist%validTimes(in), &
              '(i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2)') &
              year,month,day,hour,minute,second
         includeMinute=.true.
         includeSecond=.true.
      else
         print*,'ERROR parsing validTimes!'
         stop 1
      end if

      wpsFileData = createFileWPS(fileUnit=fileUnit, &
           outputDirectory=trim(namelist%outputDirectory), &
           prefix=trim(prefix), &
           year=year,month=month,day=day,hour=hour, &
           minute=minute,second=second, &
           includeMinute=includeMinute, &
           includeSecond=includeSecond)

      return
   end subroutine openSlabFile

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  closeSlabFile
   !
   ! DESCRIPTION:  Closes WPS file for 2d slabs.
   !
   !--------------------------------------------------------------------------
   
   subroutine closeSlabFile(wpsFileData)
      implicit none
      type(FileWPS),intent(inout) :: wpsFileData
      call destroyFileWPS(wpsFileData)
      return
   end subroutine closeSlabFile

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  writeSlab
   !
   ! DESCRIPTION:  Writes 2d slab to WPS file.
   !
   !--------------------------------------------------------------------------

   subroutine writeSlab(namelist,iv,in,xlvl,internal, &
        tmp2d,nx,ny,wpsFileData, &
        year,month,day,hour,minute,second)

      implicit none

      ! Arguments
      type(namelistGEOS2WPS), intent(in) :: namelist
      integer,intent(in) :: iv
      integer,intent(in) :: in
      real,intent(in) :: xlvl
      type(internalData),intent(in) :: internal
      integer,intent(in) :: nx,ny
      real,intent(in) :: tmp2d(nx,ny)
      type(FileWPS),intent(in) :: wpsFileData
      integer,intent(in) :: year,month,day,hour,minute,second

      ! Local variables
      type(FieldWPS) :: wpsFieldData
      logical :: includeMinute,includeSecond

      includeMinute=.true.
      includeSecond=.true.

      wpsFieldData = createFieldWPS(year=year, month=month, day=day, &
           hour=hour,minute=minute, second=second, &
           includeMinute=includeMinute, &
           includeSecond=includeSecond, &
           xfcst=namelist%forecastHours(in), &
           map_source='GEOS/NASA GMAO', &
           field=trim(namelist%variableNamesOut(iv)), &
           units=trim(namelist%variableUnits(iv)), &
           desc=trim(namelist%variableDescriptions(iv)), &
           xlvl=xlvl, &
           nx=internal%longitudeDimension, &
           ny=internal%latitudeDimension, &
           iproj=IPROJ_LATLON, startloc=SWCORNER, &
           startlat=internal%southwestLatitude, &
           startlon=internal%southwestLongitude, &
           earth_radius=internal%earthRadius, &
           is_wind_grid_rel=.false., &
           slab=tmp2d, deltaLat=internal%deltaLatitude, &
           deltaLon=internal%deltaLongitude)
      call writeFileWPS(wpsFileData,wpsFieldData)
      call destroyFieldWPS(wpsFieldData)

   end subroutine writeSlab

end module GEOS2WPS_mod

