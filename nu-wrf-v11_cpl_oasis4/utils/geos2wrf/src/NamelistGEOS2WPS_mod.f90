!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  NamelistGEOS2WPS_mod
!
! AUTHOR:
! Eric Kemp, NASA/Northrop Grumman
!
! DESCRIPTION:
! Defines data structure and routines for reading namelist input for 
! GEOS to WPS converter program.
!
! REVISION HISTORY:
! 9 Apr 2012 - Initial version
! 17 Apr 2012 - Renamed module.
! 3 Nov 2014 - Added support for subsetting.
!
!------------------------------------------------------------------------------

module NamelistGEOS2WPS_mod
   
   ! Import modules
   use FileUnit_mod
   use FileUtils_mod, only : HDF4_FORMAT
   use FileUtils_mod, only : NETCDF_FORMAT
   use FileUtils_mod, only : HDFEOS2_FORMAT

   ! Change defaults
   implicit none
   private

   ! Data structure for namelist entries
   public :: NamelistGEOS2WPS
   type NamelistGEOS2WPS
      integer :: geosFileFormat
      character(len=256) :: geosFileName
      character(len=256) :: outputDirectory
      character(len=64) :: longitudeName
      character(len=64) :: latitudeName
      logical :: hasVerticalDimension
      character(len=64) :: verticalName
      integer :: numberOfTimes
      integer,allocatable :: timeIndices(:)
      character(len=19),allocatable :: validTimes(:)
      real,allocatable :: forecastHours(:)
      integer :: numberOfVariables
      integer,allocatable :: variableRanks(:)
      integer,allocatable :: variableLevelTypes(:)
      character(len=64), allocatable :: variableNamesIn(:)
      character(len=9), allocatable :: variableNamesOut(:)
      character(len=25),allocatable :: variableUnits(:)
      character(len=46),allocatable :: variableDescriptions(:)
      logical :: Subset
      integer :: iLonMin
      integer :: iLonMax
      integer :: jLatMin
      integer :: jLatMax
      integer :: kVertMin
      integer :: kVertMax
      integer :: mTimeMin
      integer :: mTimeMax
   end type NamelistGEOS2WPS

   ! Public methods
   public :: createNamelistGEOS2WPS
   public :: destroyNamelistGEOS2WPS

   ! Public constants
   integer, parameter, public :: GROUND_LEVEL=1
   integer, parameter, public :: TWO_METERS_ABOVE_GROUND_LEVEL=2
   integer, parameter, public :: TEN_METERS_ABOVE_GROUND_LEVEL=3
   integer, parameter, public :: MEAN_SEA_LEVEL=4

   integer, parameter, public :: MODEL_LEVEL=11
   integer, parameter, public :: ISOBARIC_LEVEL=12

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  createNamelistGEOS2WPS
   !
   ! DESCRIPTION:  Public constructor for NamelistGEOS2WPS data type.
   !
   !---------------------------------------------------------------------------

   function createNamelistGEOS2WPS(filename) result (this)

      ! Arguments
      character(len=*),intent(in) :: filename

      ! Return variable
      type(NamelistGEOS2WPS) :: this

      ! Local variables
      integer,parameter :: MAX_NUMBER_OF_TIMES = 24
      integer,parameter :: MAX_NUMBER_OF_VARIABLES = 30

      integer :: geosFileFormat
      character(len=256) :: geosFileName
      character(len=256) :: outputDirectory
      namelist /files/ geosFileFormat, geosFileName, outputDirectory

      character(len=64) :: longitudeName
      character(len=64) :: latitudeName
      logical :: hasVerticalDimension
      character(len=64) :: verticalName
      namelist /coordinates/ longitudeName,latitudeName,hasVerticalDimension,&
           verticalName

      integer :: numberOfTimes
      character(len=19) :: validTimes(MAX_NUMBER_OF_TIMES)
      integer :: timeIndices(MAX_NUMBER_OF_TIMES)
      real :: forecastHours(MAX_NUMBER_OF_TIMES)
      namelist /forecast/ numberOfTimes,validTimes,timeIndices,forecastHours

      integer :: numberOfVariables
      integer :: variableRanks(MAX_NUMBER_OF_VARIABLES)
      integer :: variableLevelTypes(MAX_NUMBER_OF_VARIABLES)
      character(len=64) :: variableNamesIn(MAX_NUMBER_OF_VARIABLES)
      character(len=9)  :: variableNamesOut(MAX_NUMBER_OF_VARIABLES)
      character(len=25) :: variableUnits(MAX_NUMBER_OF_VARIABLES)
      character(len=46) :: variableDescriptions(MAX_NUMBER_OF_VARIABLES)
      namelist /variables/ numberOfVariables, &
           variableRanks,variableLevelTypes,variableNamesIn,variableNamesOut, &
           variableUnits, variableDescriptions

      logical :: subset
      integer :: iLonMin
      integer :: iLonMax
      integer :: jLatMin
      integer :: jLatMax
      integer :: kVertMin
      integer :: kVertMax
      integer :: mTimeMin
      integer :: mTimeMax
      namelist /subsetData/ subset,iLonMin,iLonMax,jLatMin,jLatMax, &
           kVertMin,kVertMax,mTimeMin,mTimeMax

      integer :: fileUnit
      integer :: year,month,day,hour,minute,second
      integer :: i

      ! Initialize data structure
      call initNamelistGEOS2WPS(this)

      ! Initialize namelist values
      geosFileFormat=0
      geosFileName='NULL'
      outputDirectory='NULL'
      longitudeName='NULL'
      latitudeName='NULL'
      hasVerticalDimension=.false.
      verticalName='NULL'
      numberOfTimes=0
      validTimes(:) = 'NULL'
      timeIndices(:)=0
      forecastHours(:) = -9999
      numberOfVariables=0
      variableRanks(:) = 0
      variableLevelTypes(:) = 0
      variableNamesIn(:) = 'NULL'
      variableNamesOut(:) = 'NULL'
      variableUnits(:) = 'NULL'
      variableDescriptions(:) = 'NULL'

      subset=.false.
      iLonMin=0
      iLonMax=0
      jLatMin=0
      jLatMax=0
      kVertMin=0
      kVertMax=0
      mTimeMin=0
      mTimeMax=0

      ! Read from namelist file
      fileUnit=select_file_unit()
      open(unit=fileUnit,file=trim(filename),delim='APOSTROPHE')

      read(unit=fileUnit,nml=files)
      if (geosFileFormat .ne. HDF4_FORMAT .and. &
          geosFileFormat .ne. NETCDF_FORMAT .and. &
          geosFileFormat .ne. HDFEOS2_FORMAT) then
         print*,'ERROR, invalid format selected for GEOS file!'
         print*,'Read in ',geosFileFormat
         print*,'Use 1 for HDF4, 2 for netCDF'
         stop 1
      end if
      if (trim(geosFileName) == 'NULL') then
         print*,'ERROR assigning name to geosFileName!'
         print*,'Read in ',trim(geosFileName)
         stop 1
      end if
      if (trim(outputDirectory) == 'NULL') then
         print*,'ERROR assigning name to outputDirectory!'
         print*,'Read in ',trim(outputDirectory)
         stop 1
      end if
      this%geosFileFormat = geosFileFormat
      this%geosFileName = trim(geosFileName)
      this%outputDirectory = trim(outputDirectory)

      read(unit=fileUnit,nml=coordinates)
      if (trim(longitudeName) == 'NULL') then
         print*,'ERROR assigning name to longitude variable!'
         stop 1
      end if
      if (trim(latitudeName) == 'NULL') then
         print*,'ERROR assigning name to latitude variable!'
         stop 1
      end if
      this%longitudeName=trim(longitudeName)
      this%latitudeName=trim(latitudeName)
      if (hasVerticalDimension) then
         if (trim(verticalName) == 'NULL') then
            print*,'ERROR assigning name to vertical variable!'
            stop 1
         end if
      else
         if (trim(verticalName) .ne. 'NULL') then
            print*,'ERROR, vertical variable assigned, but '
            print*,'  hasVerticalDimension is .false.'
            print*,'  verticalName:  ',trim(verticalName)
            stop 1
         end if
      end if
      this%hasVerticalDimension = hasVerticalDimension
      this%verticalName = trim(verticalName)

      read(unit=fileUnit,nml=forecast)
      if (numberOfTimes <= 0 .or. numberOfTimes > MAX_NUMBER_OF_TIMES) then
         print*,'ERROR, invalid value for numberOfTimes!'
         print*,'Read in ',numberOfTimes
         print*,'Must be 0 < numberOfTimes <= ',MAX_NUMBER_OF_TIMES
         stop 1
      end if
      do i = 1, numberOfTimes
         if (trim(validTimes(i)) == 'NULL') then
            print*,'ERROR assigning validTimes for i = ',i
            print*,'Value is ',trim(validTimes(i))
            stop 1
         end if
         print*,'Will attempt to extract year/month/day/hour'
         print*,'from validTimes(i) = ',validTimes(i)
         ! See if minute and second are available
         if (len_trim(validTimes(i)) .eq. 13) then 
            read(validTimes(i),'(i4.4,1x,i2.2,1x,i2.2,1x,i2.2)') &
                 year,month,day,hour
            minute = 0
            second = 0
         else if (len_trim(validTimes(i)) .eq. 16) then 
            read(validTimes(i),'(i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2)') &
                 year,month,day,hour,minute
            second = 0
         else if (len_trim(validTimes(i)) .eq. 19) then
              read(validTimes(i), &
                   '(i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2)') &
              year,month,day,hour,minute,second
         else
            print*,'ERROR parsing validTimes!'
            stop 1
         end if

      end do
      do i = 1, numberOfTimes
         if (timeIndices(i) == 0) then
            print*,'ERROR, invalid time index given for i= ',i
            print*,'Read in ',timeIndices(i)
            print*,'Must be > 0!'
            stop 1
         end if
      end do
      do i = 1, numberOfTimes
         if (forecastHours(i) < 0) then
            print*,'ERROR assigning forecastHours for i=',i
            print*,'Value is ',forecastHours(i)
            stop 1
         end if
      end do
      this%numberOfTimes=numberOfTimes
      allocate(this%validTimes(numberOfTimes))
      do i = 1, numberOfTimes
         this%validTimes(i) = validTimes(i)
      end do
      allocate(this%timeIndices(numberOfTimes))
      do i = 1, numberOfTimes
         this%timeIndices(i) = timeIndices(i)
      end do
      allocate(this%forecastHours(i))
      do i = 1, numberOfTimes
         this%forecastHours(i) = forecastHours(i)
      end do

      read(unit=fileUnit,nml=variables)
      if (numberOfVariables <= 0 .or. &
          numberOfVariables > MAX_NUMBER_OF_VARIABLES) then
         print*,'ERROR, invalid value for numberOfVariables!'
         print*,'Read in ',numberOfVariables
         print*,'Valid range is 0 <= numberOfVariables < ', &
              MAX_NUMBER_OF_VARIABLES
         stop 1
      end if
      do i = 1, numberOfVariables
         if (variableRanks(i) .ne. 3 .and. variableRanks(i) .ne. 4) then
            print*,'ERROR, bad value of variableRanks for i = ',i
            print*,'Read in ',variableRanks(i)
            print*,'Must be either 3 or 4'
            stop 1
         end if
      end do
      do i = 1, numberOfVariables
         if (variableLevelTypes(i) .ne. GROUND_LEVEL .and. &
             variableLevelTypes(i) .ne. ISOBARIC_LEVEL .and. &
             variableLevelTypes(i) .ne. MEAN_SEA_LEVEL .and. &
             variableLevelTypes(i) .ne. TWO_METERS_ABOVE_GROUND_LEVEL .and. &
             variableLevelTypes(i) .ne. TEN_METERS_ABOVE_GROUND_LEVEL .and. &
             variableLevelTypes(i) .ne. MODEL_LEVEL) then
            print*,'ERROR, bad value of variableLevelTypes for i = ',i
            print*,'Read in ',variableLevelTypes(i)
            print*,'Valid values are: '
            print*,'    For ground level: ',GROUND_LEVEL
            print*,'    For 2 meters AGL: ',TWO_METERS_ABOVE_GROUND_LEVEL
            print*,'    For 10 meters AGL: ',TEN_METERS_ABOVE_GROUND_LEVEL
            print*,'    For mean sea level: ',MEAN_SEA_LEVEL
            print*,'    For model level: ',MODEL_LEVEL
            print*,'    For isobaric level: ',ISOBARIC_LEVEL  
            stop 1
         end if
      end do
      do i = 1, numberOfVariables
         if (trim(variableNamesIn(i)) .eq. 'NULL') then
            print*,'ERROR assigning value to variableNamesIn for i=',i
            print*,'Read in ',variableNamesIn(i)
            stop 1
         end if
      end do
      do i = 1, numberOfVariables
         if (trim(variableNamesOut(i)) .eq. 'NULL') then
            print*,'ERROR assigning value to variableNamesOut for i=',i
            print*,'Read in ',variableNamesOut(i)
            stop 1
         end if
      end do
      do i = 1, numberOfVariables
         if (trim(variableUnits(i)) .eq. 'NULL') then
            print*,'ERROR assigning value to variableUnits for i=',i
            print*,'Read in ',variableUnits(i)
            stop 1
         end if
      end do
      do i = 1, numberOfVariables
         if (trim(variableDescriptions(i)) .eq. 'NULL') then
            print*,'ERROR assigning value to variableDescriptions for i=',i
            print*,'Read in ',variableDescriptions(i)
            stop 1
         end if
      end do
      this%numberOfVariables = numberOfVariables
      allocate(this%variableLevelTypes(numberOfVariables))
      do i = 1, numberOfVariables
         this%variableLevelTypes(i) = variableLevelTypes(i)
      end do
      allocate(this%variableRanks(numberOfVariables))
      do i = 1, numberOfVariables
         this%variableRanks(i) = variableRanks(i)
      end do
      allocate(this%variableNamesIn(numberOfVariables))
      do i = 1, numberOfVariables
         this%variableNamesIn(i) = variableNamesIn(i)
      end do
      allocate(this%variableNamesOut(numberOfVariables))
      do i = 1, numberOfVariables
         this%variableNamesOut(i) = variableNamesOut(i)
      end do
      allocate(this%variableUnits(numberOfVariables))
      do i = 1, numberOfVariables
         this%variableUnits(i) = variableUnits(i)
      end do
      allocate(this%variableDescriptions(numberOfVariables))
      do i = 1, numberOfVariables
         this%variableDescriptions(i) = variableDescriptions(i)
      end do

      read(unit=fileUnit,nml=subsetData)
      this%subset = subset
      this%iLonMin  = iLonMin
      this%iLonMax  = iLonMax
      this%jLatMin  = jLatMin
      this%jLatMax  = jLatMax
      this%kVertMin = kVertMin
      this%kVertMax = kVertMax
      this%mTimeMin = mTimeMin
      this%mTimeMax = mTimeMax

      call close_file_unit(fileUnit)

      return
   end function createNamelistGEOS2WPS

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  destroyNamelistGEOS2WPS
   !
   ! DESCRIPTION:  Public destructor for NamelistGEOS2WPS data type.
   !
   !---------------------------------------------------------------------------

   subroutine destroyNamelistGEOS2WPS(this)

      ! Arguments
      type(NamelistGEOS2WPS), intent(inout) :: this

      call initNamelistGEOS2WPS(this)

      return

   end subroutine destroyNamelistGEOS2WPS

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  initNamelistGEOS2WPS
   !
   ! DESCRIPTION:  Private method for initializing data members and 
   ! deallocating member arrays for NamelistGEOS2WPS data type.
   !
   !---------------------------------------------------------------------------

   subroutine initNamelistGEOS2WPS(this)

      ! Arguments
      type(NamelistGEOS2WPS),intent(inout) :: this

      this%geosFileFormat=0
      this%geosFileName='NULL'
      this%outputDirectory='NULL'
      this%longitudeName='NULL'
      this%latitudeName='NULL'
      this%hasVerticalDimension=.false.
      this%verticalName='NULL'
      this%numberOfTimes=0
      if (allocated(this%validTimes)) deallocate(this%validTimes)
      if (allocated(this%timeIndices)) deallocate(this%timeIndices)
      if (allocated(this%forecastHours)) deallocate(this%forecastHours)
      this%numberOfVariables = 0
      if (allocated(this%variableRanks)) deallocate(this%variableRanks)
      if (allocated(this%variableLevelTypes)) &
           deallocate(this%variableLevelTypes)
      if (allocated(this%variableNamesIn)) deallocate(this%variableNamesIn)
      if (allocated(this%variableNamesOut)) deallocate(this%variableNamesOut)
      if (allocated(this%variableUnits)) deallocate(this%variableUnits)
      if (allocated(this%variableDescriptions)) &
           deallocate(this%variableDescriptions)

      this%subset=.false.
      this%iLonMin  = 0
      this%iLonMax  = 0
      this%jLatMin  = 0
      this%jLatMax  = 0
      this%kVertMin = 0
      this%kVertMax = 0
      this%mTimeMin = 0
      this%mTimeMax = 0

      return
   end subroutine initNamelistGEOS2WPS
end module NamelistGEOS2WPS_mod
