!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  createPRESSURE
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Back-end utility program for GEOS2WRF package.  Reads in WPS intermediate
! format file with layer pressure thicknesses, calculates 3D layer pressures
! (at mid-points of the layers), and write out in a new WPS file.
!
! REVISION:
! 3 May 2012 - First version
!
!------------------------------------------------------------------------------

program createPRESSURE

   ! Import modules
   use DerivedVars_mod
   use FieldWPS_mod
   use FileUnit_mod
   use FileWPS_mod

   ! Change defaults
   implicit none

   ! Local variables
   type(FieldWPS) :: field 
   type(FieldWPS), allocatable :: layerPressureThicknessSlabs(:)
   type(FileWPS) :: inputFile, outputFile
   integer :: layerCount, numberOfSlabs
   logical :: endOfFile
   integer :: fileUnit1,fileUnit2
   real, allocatable :: layerPressureThicknesses(:,:,:)
   real, allocatable :: edgePressures(:,:,:)
   real, allocatable :: layerPressures(:,:,:)
   integer :: longitudeDimension, latitudeDimension, verticalDimension
   logical :: isCorrectTime
   integer :: i,j,k
   
   ! Namelist variables
   character(len=132) :: directory
   character(len=64) :: prefix
   integer :: year,month,day,hour,minute,second
   logical :: includeMinute,includeSecond
   character(len=64) :: layerPressureThicknessName
   real :: modelTopPressure
   namelist /input/ directory,prefix,year,month,day,hour,minute,second,  &
        includeMinute,includeSecond, &
        layerPressureThicknessName,modelTopPressure

   ! Read namelist file
   directory='NULL'
   prefix='NULL'
   year=0000
   month=00
   day=00
   hour=00
   minute=00
   second=00
   includeMinute=.false.
   includeSecond=.false.
   layerPressureThicknessName='NULL'
   modelTopPressure=1. ! Default for GEOS-5 is 0.01 hPa = 1 Pa

   fileUnit1=select_file_unit()
   open(unit=fileUnit1,file='namelist.createPRESSURE',delim='APOSTROPHE')
   read(unit=fileUnit1,nml=input)
   call close_file_unit(fileUnit1)

   if (trim(directory) .eq. 'NULL') then
      print*,'ERROR assigning value to directory in namelist!'
      stop 1
   end if
   if (trim(prefix) .eq. 'NULL') then
      print*,'ERROR assigning value to prefix in namelist!'
      stop 1
   end if
   if (year .eq. 0) then
      print*,'ERROR assigning value to year in namelist!'
      stop 1
   end if
   if (month .eq. 0) then
      print*,'ERROR assigning value to month in namelist!'
      stop 1
   end if
   if (day .eq. 0) then
      print*,'ERROR assigning value to day in namelist!'
      stop 1
   end if
   if (hour < 0 .or. hour > 23) then
      print*,'ERROR assigning valid value to hour in namelist!'
      stop 1
   end if

   if (minute < 0 .or. minute > 59) then
      print*,'ERROR assigning valid value to minute in namelist!'
      stop 1
   end if
   if (second < 0 .or. second > 59) then
      print*,'ERROR assigning valid value to second in namelist!'
      stop 1
   end if
   if (trim(layerPressureThicknessName) .eq. 'NULL') then
      print*,'ERROR assigning value to layerPressureThicknessName in namelist!'
      stop 1
   end if
   if ( .not. (modelTopPressure > 0)) then
      print*,'ERROR, non-positive modelTopPressure in namelist!'
      stop 1
   end if

   ! Open the GEOS WPS file
   fileUnit1=select_file_unit()
   inputFile=createFileWPS(fileUnit1,trim(directory),trim(prefix),year,month, &
        day,hour,minute,second,includeMinute,includeSecond,preserve=.true.)

   ! Loop through the GEOS WPS file looking for pressureThicknesses
   layerCount=0
   do
      ! Read the next slab
      call readFileWPS(inputFile,field,endOfFile)
      if (endOfFile) exit

      ! Skip if wrong date/time
      isCorrectTime = checkDateTimeFieldWPS(field,year,month,day,hour, &
           minute,second)
      if (.not. isCorrectTIme) cycle

      if (trim(field%field) .eq. trim(layerPressureThicknessName)) then
         layerCount=layerCount + 1
      end if
      call destroyFieldWPS(field)
   end do
   call destroyFileWPS(inputFile)
   call close_file_unit(fileUnit1)

   ! Check count
   if (layerCount .eq. 0) then
      print*,'ERROR, found no data for layerPressureThicknesses!'
      print*,'Searched for ',trim(layerPressureThicknessName)
      stop 1
   end if

   ! Allocate array for slabs
   numberOfSlabs = layerCount
   allocate(layerPressureThicknessSlabs(numberOfSlabs))
   
   ! Reopen the GEOS WPS file
   fileUnit1=select_file_unit()
   inputFile=createFileWPS(fileUnit1,trim(directory),trim(prefix),year,month, &
        day,hour,minute,second,includeMinute,includeSecond,preserve=.true.)

   ! Loop through each slab and copy layer pressure thicknesses to array.
   layerCount=0
   do
      ! Read the next slab
      call readFileWPS(inputFile,field,endOfFile)
      if (endOfFile) exit

      ! Skip if wrong date/time
      isCorrectTime = checkDateTimeFieldWPS(field,year,month,day,hour, &
           minute,second)
      if (.not. isCorrectTIme) cycle

      if (trim(field%field) .eq. trim(layerPressureThicknessName)) then
         layerCount=layerCount + 1
         layerPressureThicknessSlabs(layerCount) = copyFieldWPS(field)
      end if
      call destroyFieldWPS(field)
   end do
   call destroyFileWPS(inputFile)
   call close_file_unit(fileUnit1)

   ! Make sure thickness slabs are all on same grid.
   do i = 2, numberOfSlabs
      call compareGridsFieldWPS(layerPressureThicknessSlabs(1), &
                                layerPressureThicknessSlabs(i))
   end do

   ! Check for duplicate slab levels
   call checkDuplicatesXlvls(numberOfSlabs,layerPressureThicknessSlabs)

   ! Make sure pressure thicknesses are sorted by xlvl
   call bubbleSortXlvl(numberOfSlabs,layerPressureThicknessSlabs)

   ! Allocate 3D arrays for calculating pressures.  We'll deal with FieldWPS
   ! data types after the values are calculated.
   longitudeDimension=layerPressureThicknessSlabs(1)%nx
   latitudeDimension=layerPressureThicknessSlabs(1)%ny
   verticalDimension=numberOfSlabs

   ! Calculate edge pressures
   call copySlabsToArray3d(verticalDimension,layerPressureThicknessSlabs, &
        layerPressureThicknesses)
   call calcEdgePressures(longitudeDimension,latitudeDimension, &
        verticalDimension,modelTopPressure,layerPressureThicknesses, &
        edgePressures)

   ! Calculate layer pressure
   call calcLayerPressures(longitudeDimension,latitudeDimension, &
        verticalDimension,edgePressures,layerPressures)

   ! Open output WPS file
   fileUnit2=select_file_unit()
   outputFile=createFileWPS(fileUnit2,trim(directory),'PRESSURE_MODEL_LEVEL',&
        year,month,day,hour,minute,second,includeMinute,includeSecond)

   ! Loop through layerPressures, create slabs, write to file
   do k = 1, verticalDimension
      field=copyFieldWPS(layerPressureThicknessSlabs(k))
      field%field="PRESSURE"
      field%units="Pa"
      field%desc="Pressure"
      do j = 1, latitudeDimension
         do i = 1, longitudeDimension
            field%slab(i,j)=layerPressures(i,j,k)
         end do
      end do
      call writeFileWPS(outputFile,field)
      call destroyFieldWPS(field)
   end do

   ! Clean up
   do k = 1, verticalDimension
      call destroyFieldWPS(layerPressureThicknessSlabs(k))
   end do
   deallocate(layerPressureThicknessSlabs)
   deallocate(layerPressureThicknesses)
   deallocate(layerPressures)
   deallocate(edgePressures)
   
end program createPRESSURE
