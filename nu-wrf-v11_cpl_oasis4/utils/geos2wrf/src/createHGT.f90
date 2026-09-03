!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  createHGT
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Back-end utility program for GEOS2WRF package.  Reads in WPS intermediate
! format file with layer pressure thicknesses, temperature, and specific
! humidity, calculates 3D geopotential heights (at mid-points of the layers),
! and write out in a new WPS file.
!
! REVISION:
! 9 May 2012 - First version
!
!------------------------------------------------------------------------------

program createHGT

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
   type(FieldWPS), allocatable :: layerTemperatureSlabs(:)
   type(FieldWPS), allocatable :: layerSpecificHumiditySlabs(:)
   type(FieldWPS) :: soilHeightSlab
   type(FileWPS) :: inputFile, outputFile
   integer :: numberOfSlabs
   integer :: pressureThicknessCount, temperatureCount, specificHumidityCount
   integer :: soilHeightCount
   logical :: endOfFile
   integer :: fileUnit1,fileUnit2
   real, allocatable :: layerPressureThicknesses(:,:,:)
   real, allocatable :: layerTemperatures(:,:,:)
   real, allocatable :: layerSpecificHumidities(:,:,:)
   real, allocatable :: edgePressures(:,:,:)
   real, allocatable :: layerPressures(:,:,:)
   real, allocatable :: soilHeights(:,:)
   real, allocatable :: edgeGeopotentialHeights(:,:,:)
   real, allocatable :: layerGeopotentialHeights(:,:,:)
   integer :: longitudeDimension, latitudeDimension, verticalDimension
   logical :: isCorrectTime
   integer :: i,j,k
   
   ! Namelist variables
   character(len=132) :: directory
   character(len=64) :: prefix
   integer :: year,month,day,hour,minute,second
   logical :: includeMinute,includeSecond
   character(len=64) :: layerPressureThicknessName
   character(len=64) :: layerTemperatureName   
   character(len=64) :: layerSpecificHumidityName
   character(len=64) :: soilHeightName
   real :: modelTopPressure
   namelist /input/ directory,prefix,year,month,day,hour,minute,second, &
        includeMinute,includeSecond, &
        layerPressureThicknessName,layerTemperatureName, &
        layerSpecificHumidityName,soilHeightName,modelTopPressure

   real, parameter :: SURFACE_XLVL = 200100.

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
   layerTemperatureName='NULL'
   layerSpecificHumidityName='NULL'
   soilHeightName='NULL'
   modelTopPressure=1. ! Default for GEOS-5 is 0.01 hPa = 1 Pa

   fileUnit1=select_file_unit()
   open(unit=fileUnit1,file='namelist.createHGT',delim='APOSTROPHE')
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
   if (trim(layerTemperatureName) .eq. 'NULL') then
      print*,'ERROR assigning value to layerTemperatureName in namelist!'
      stop 1
   end if
   if (trim(layerSpecificHumidityName) .eq. 'NULL') then
      print*,'ERROR assigning value to layerSpecificHumidityName in namelist!'
      stop 1
   end if
   if (trim(soilHeightName) .eq. 'NULL') then
      print*,'ERROR assigning value to soilHeightName in namelist!'
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

   ! Loop through the GEOS WPS file looking for pressureThicknesses, 
   ! temperatures, and specific humidities
   pressureThicknessCount=0
   temperatureCount=0
   specificHumidityCount=0
   soilHeightCount=0
   do
      ! Read the next slab
      call readFileWPS(inputFile,field,endOfFile)
      if (endOfFile) exit

      ! Skip if wrong date/time
      isCorrectTime = checkDateTimeFieldWPS(field,year,month,day,hour, &
           minute,second)
      if (.not. isCorrectTIme) cycle

      ! Increment counts
      if (trim(field%field) .eq. trim(layerPressureThicknessName)) then
         pressureThicknessCount = pressureThicknessCount + 1
      else if (trim(field%field) .eq. trim(layerTemperatureName)) then
         if (field%xlvl .ne. SURFACE_XLVL) then
            temperatureCount = temperatureCount + 1
         end if
      else if (trim(field%field) .eq. trim(layerSpecificHumidityName)) then
         if (field%xlvl .ne. SURFACE_XLVL) then
            specificHumidityCount = specificHumidityCount + 1
         end if
      else if (trim(field%field) .eq. trim(soilHeightName)) then
         soilHeightCount = soilHeightCount + 1
      end if
      call destroyFieldWPS(field)
   end do
   call destroyFileWPS(inputFile)
   call close_file_unit(fileUnit1)
   
   ! Check counts
   if (pressureThicknessCount .ne. temperatureCount .or. &
       pressureThicknessCount .ne. specificHumidityCount) then
      print*,'ERROR, mismatched counts of pressure thicknesses, temperature,'
      print*,'and specific humidity slabs!'
      print*,'    Pressure thickness slabs:  ',pressureThicknessCount
      print*,'    Temperature slabs:  ',temperatureCount
      print*,'    Specific humidity slabs:  ',specificHumidityCount
      stop 1
   end if
   if (pressureThicknessCount .eq. 0) then
      print*,'ERROR, no pressure thickness data found in WPS file!'
      stop 1
   end if
   if (soilHeightCount .ne. 1) then
      print*,'ERROR, expected single slab of soil height data in WPS file!'
      print*,'Instead found ',soilHeightCount
      stop 1
   end if

   ! Allocate arrays from slabs
   numberOfSlabs=pressureThicknessCount
   allocate(layerPressureThicknessSlabs(numberOfSlabs))
   allocate(layerTemperatureSlabs(numberOfSlabs))
   allocate(layerSpecificHumiditySlabs(numberOfSlabs))

   ! Reopen the GEOS WPS file
   fileUnit1=select_file_unit()
   inputFile=createFileWPS(fileUnit1,trim(directory),trim(prefix),year,month, &
        day,hour,minute,second,includeMinute,includeSecond,preserve=.true.)

   ! Loop through each slab and copy pressureThicknesses, temperatures, and
   ! specific humidities
   pressureThicknessCount=0
   temperatureCount=0
   specificHumidityCount=0
   do
      ! Read the next slab
      call readFileWPS(inputFile,field,endOfFile)
      if (endOfFile) exit

      ! Skip if wrong date/time
      isCorrectTime = checkDateTimeFieldWPS(field,year,month,day,hour, &
           minute,second)
      if (.not. isCorrectTIme) cycle

      if (trim(field%field) .eq. trim(layerPressureThicknessName)) then
         pressureThicknessCount=pressureThicknessCount + 1
         layerPressureThicknessSlabs(pressureThicknessCount) = &
              copyFieldWPS(field)
      else if (trim(field%field) .eq. trim(layerTemperatureName)) then
         if (field%xlvl .ne. SURFACE_XLVL) then
            temperatureCount = temperatureCount + 1
            layerTemperatureSlabs(temperatureCount) = copyFieldWPS(field)
         end if
      else if (trim(field%field) .eq. trim(layerSpecificHumidityName)) then
         if (field%xlvl .ne. SURFACE_XLVL) then
            specificHumidityCount = specificHumidityCount + 1
            layerSpecificHumiditySlabs(specificHumidityCount) = &
                 copyFieldWPS(field)
         end if
      else if (trim(field%field) .eq. trim(soilHeightName)) then
         soilHeightSlab = copyFieldWPS(field)
      end if
      call destroyFieldWPS(field)
   end do
   call destroyFileWPS(inputFile)
   call close_file_unit(fileUnit1)

   ! Make sure slabs are all on the same grid
   do i = 2, numberOfSlabs
      call compareGridsFieldWPS(layerPressureThicknessSlabs(1), &
                                layerPressureThicknessSlabs(i))
   end do
   do i = 1, numberOfSlabs
      call compareGridsFieldWPS(layerPressureThicknessSlabs(1), &
                                layerTemperatureSlabs(i))
   end do
   do i = 1, numberOfSlabs
      call compareGridsFieldWPS(layerPressureThicknessSlabs(1), &
                                layerSpecificHumiditySlabs(i))
   end do
   call compareGridsFieldWPS(layerPressureThicknessSlabs(1), &
                        soilHeightSlab)
   
   ! Check for duplicate slab levels
   call checkDuplicatesXlvls(numberOfSlabs,layerPressureThicknessSlabs)
   call checkDuplicatesXlvls(numberOfSlabs,layerTemperatureSlabs)
   call checkDuplicatesXlvls(numberOfSlabs,layerSpecificHumiditySlabs)

   ! Make sure slabs are sorted by xlvl
   call bubbleSortXlvl(numberOfSlabs,layerPressureThicknessSlabs)
   call bubbleSortXlvl(numberOfSlabs,layerTemperatureSlabs)
   call bubbleSortXlvl(numberOfSlabs,layerSpecificHumiditySlabs)

   ! Allocate 3D arrays for calculating pressures.  
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

   ! Calculate edge geopotential heights
   call copySlabToArray2d(soilHeightSlab,soilHeights)
   call copySlabsToArray3d(numberOfSlabs,layerTemperatureSlabs, &
        layerTemperatures)
   call copySlabsToArray3d(numberOfSlabs,layerSpecificHumiditySlabs, &
        layerSpecificHumidities)
   call calcEdgeGeopotentialHeights(longitudeDimension,latitudeDimension, &
        verticalDimension,edgePressures,layerTemperatures, &
        layerSpecificHumidities,soilHeights,edgeGeopotentialHeights)

   ! Calculate layer geopotential heights
   call calcLayerGeopotentialHeights(longitudeDimension,latitudeDimension, &
        verticalDimension,edgePressures,layerPressures, &
        layerTemperatures,layerSpecificHumidities,edgeGeopotentialHeights, &
        layerGeopotentialHeights)
   
   ! Open output WPS file
   fileUnit2=select_file_unit()
   outputFile=createFileWPS(fileUnit2,trim(directory),'HGT_MODEL_LEVEL', &
        year,month,day,hour,minute,second,includeMinute,includeSecond)

   ! Loop through layerGeopotentialHeights, create slabs, and write to file
   do k = 1, verticalDimension
      field=copyFieldWPS(layerPressureThicknessSlabs(k))
      field%field="HGT"
      field%units="m"
      field%desc="Height"
      do j = 1, latitudeDimension
         do i = 1, longitudeDimension
            field%slab(i,j)=layerGeopotentialHeights(i,j,k)
         end do
      end do
      call writeFileWPS(outputFile,field)
      call destroyFieldWPS(field)
   end do

   ! Close file
   call destroyFileWPS(outputFile)
   call close_file_unit(fileUnit2)

   ! Clean up
   do k = 1, verticalDimension
      call destroyFieldWPS(layerPressureThicknessSlabs(k))
      call destroyFieldWPS(layerTemperatureSlabs(k))
      call destroyFielDWPS(layerSpecificHumiditySlabs(k))
   end do
   call destroyFieldWPS(soilHeightSlab)
   deallocate(layerPressureThicknessSlabs)
   deallocate(layerTemperatureSlabs)
   deallocate(layerSpecificHumiditySlabs)
   deallocate(layerPressureThicknesses)
   deallocate(layerTemperatures)
   deallocate(layerSpecificHumidities)
   deallocate(soilHeights)
   deallocate(edgePressures)
   deallocate(layerPressures)
   deallocate(edgeGeopotentialHeights)
   deallocate(layerGeopotentialHeights)

end program createHGT
