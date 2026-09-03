!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  createRH
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Back-end utility program for GEOS2WRF package.  Reads in WPS intermediate
! format file with 3D pressure, 3D temperature, and 3D specific humidity,
! and calculates 3D relative humidity w.r.t. liquid.  This new field is then 
! written to a new WPS file.
!
! REVISION:
! 02 May 2012 - First version
! 17 May 2012 - Revisions to speed up code.
!------------------------------------------------------------------------------

program createRH
   
   ! Import modules
   use FieldWPS_mod
   use FieldWPSList_mod
   use FileUnit_mod
   use FileWPS_mod
   use WrfRH_mod

   ! Defaults
   implicit none

   ! Local variables
   type(FieldWPS) :: field, pressure, temperature, specificHumidity, &
        relativeHumidity
   type(FieldWPSList) :: pressureList, temperatureList, specificHumidityList
   type(FileWPS) :: inputFile, outputFile
   real, allocatable :: pressureXLvls(:)
   type(FieldWPS), allocatable :: relativeHumidities(:)
   integer :: pressureCount, temperatureCount, specificHumidityCount
   integer :: numberOfSlabs
   logical :: endOfFile
   integer :: fileUnit1,fileUnit2
   real :: temperatureK, pressurePa, specificHumidityKgKg, &
        satSpecificHumidityKgKg
   real :: relativeHumidityPct
   real :: mixingRatioKgKg, satMixingRatioKgKg
   logical :: isCorrectTime
   integer :: i,j

   ! Namelist variables
   character(len=132) :: directory
   character(len=64) :: prefix
   integer :: year,month,day,hour,minute,second
   logical :: includeMinute,includeSecond
   logical :: processSurfacePressure
   logical :: onIsobaricLevels
   character(len=64) :: surfacePressureName
   character(len=64) :: pressureName
   character(len=64) :: temperatureName
   character(len=64) :: specificHumidityName
   namelist /input/ directory,prefix,year,month,day,hour,minute,second, &
        includeMinute,includeSecond, &
        processSurfacePressure,surfacePressureName, &
        pressureName,temperatureName,specificHumidityName, &
        onIsobaricLevels
   
   real,parameter :: SURFACE_XLVL = 200100.

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
   processSurfacePressure = .true.
   pressureName='NULL'
   surfacePressureName='NULL'
   temperatureName='NULL'
   specificHumidityName='NULL'
   onIsobaricLevels=.false.
   fileUnit1=select_file_unit()
   open(unit=fileUnit1,file='namelist.createRH',delim='APOSTROPHE')
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
   if (year .eq. 0000) then
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
   if (trim(pressureName) .eq. 'NULL') then
      print*,'ERROR assigning value to pressureName in namelist!'
      stop 1
   end if
   if (processSurfacePressure) then
      if (trim(surfacePressureName) .eq. 'NULL') then
       print*,'ERROR assigning value to surfacePressureName in namelist!'
         print*,'processSurfacePressure is set to .true.'
         stop 1
      end if
   end if
   if (trim(temperatureName) .eq. 'NULL') then
      print*,'ERROR assigning value to temperatureName in namelist!'
      stop 1
   end if
   if (trim(specificHumidityName) .eq. 'NULL') then
      print*,'ERROR assigning value to specificHumidityName in namelist!'
      stop 1
   end if

   ! Open the GEOS WPS file
   fileUnit1=select_file_unit()
   inputFile=createFileWPS(fileUnit1,trim(directory),trim(prefix),year,month, &
        day,hour,minute,second,includeMinute,includeSecond,preserve=.true.)

   ! Loop through the GEOS WPS file looking for pressure, temperature, and 
   ! specific humidity slabs.  We will add these slabs to linked lists
   ! as they are encountered, and count them.
   pressureCount = 0
   temperatureCount = 0
   specificHumidityCount = 0
   pressureList=createFieldWPSList()
   temperatureList=createFieldWPSList()
   specificHumidityList=createFieldWPSList()
   do
      ! Read the next slab
      call readFileWPS(inputFile,field,endOfFile)      
      if (endOfFile) exit

      ! Skip if wrong date/time
      isCorrectTime = checkDateTimeFieldWPS(field,year,month,day,hour, &
           minute,second)
      if (.not. isCorrectTime) then
         call destroyFieldWPS(field)
         cycle
      end if

      if (trim(field%field) .eq. trim(surfacePressureName)) then
         if (processSurfacePressure) then
            pressureCount = pressureCount + 1
            call appendFieldWPSList(pressureList,field)
         end if
      else if (trim(field%field) .eq. trim(pressureName) .and. &
               .not. onIsobaricLevels) then
         pressureCount = pressureCount + 1
         call appendFieldWPSList(pressureList,field)
      else if (trim(field%field) .eq. trim(temperatureName)) then
         temperatureCount = temperatureCount + 1
         call appendFieldWPSList(temperatureList,field)
         if (onIsobaricLevels .and. field%xlvl .ne. SURFACE_XLVL) then
            pressureCount = pressureCount + 1
            field%field=trim(pressureName)
            field%units='Pa'
            field%desc='Pressure'
            do j = 1, field%ny
               do i = 1, field%nx
                  field%slab(i,j) = field%xlvl
               end do
            end do
            call appendFieldWPSList(pressureList,field)
         end if
      else if (trim(field%field) .eq. trim(specificHumidityName)) then
         specificHumidityCount = specificHumidityCount + 1
         call appendFieldWPSList(specificHumidityList,field)
      end if

      call destroyFieldWPS(field)

   end do

   call destroyFileWPS(inputFile)
   call close_file_unit(fileUnit1)
   
   ! Check counts
   if (pressureCount .ne. temperatureCount .or. &
       pressureCount .ne. specificHumidityCount) then
      print*,'ERROR, mismatched counts of pressure, temperature, and '
      print*,'       specific humidity slabs!'
      print*,'    Pressure slabs:  ',pressureCount
      print*,'    Temperature slabs:  ',temperatureCount
      print*,'    Specific humidity slabs:  ',specificHumidityCount
      if (processSurfacePressure) then
         print*,'NOTE:  processSurfacePressure is set to .true.'
      else
         print*,'NOTE:  processSurfacePressure is set to .false.'
      end if
      stop 1
   end if
   if (pressureCount .eq. 0) then
      print*,'ERROR, no pressure data found in WPS file!'
      stop 1
   end if
   numberOfSlabs=pressureCount

   ! Collect xlvls for pressure, temperature, and specific humidity
   allocate(pressureXlvls(numberOfSlabs))
   call getXlvlsFieldWPSList(pressureList,numberOfSlabs,pressureXlvls)
   
   ! Loop through pressures, find matching temperatures and specific 
   ! humidities, and calculate relative humidity.
   pressureCount=1
   allocate(relativeHumidities(numberOfSlabs))
   do
      if (pressureCount > numberOfSlabs) exit

      if (pressureXlvls(pressureCount) .eq. SURFACE_XLVL .and. &
          (processSurfacePressure)) then
         pressure=findFieldWPSList(pressureList,trim(surfacePressureName), &
              pressureXlvls(pressureCount))
      else 
         pressure=findFieldWPSList(pressureList,trim(pressureName), &
              pressureXlvls(pressureCount))
      end if

      ! First search for corresponding temperature slab.
      temperature=findFieldWPSList(temperatureList,trim(temperatureName),&
           pressureXlvls(pressureCount))

      ! Next, search for corresponding specific humidity slab
      specificHumidity=findFieldWPSList(specificHumidityList, &
           trim(specificHumidityName),pressureXlvls(pressureCount))

      ! Make sure data are on same grid.
      call compareGridsFieldWPS(pressure,temperature)
      call compareGridsFieldWPS(pressure,specificHumidity)

      ! Next, calculate relative humidity
      relativeHumidity = copyFieldWPS(specificHumidity)
      relativeHumidity%field="RH"
      relativeHumidity%units="%"
      if (relativeHumidity%xlvl .eq. SURFACE_XLVL) then
         relativeHumidity%desc="Relative Humidity at 2 m"
      else
         relativeHumidity%desc="Relative Humidity"
      end if

      do j = 1, pressure%ny
         do i = 1, pressure%nx
            temperatureK = temperature%slab(i,j)

            ! If this is isobaric data, underground values will be garbage.
            ! Check and set to a missing value, which will be corrected
            ! by the extrapIsobaric utility.
            if (temperatureK > 9.9999999e+10) then
               relativeHumidity%slab(i,j) = 1.e+15
               cycle
            end if

            pressurePa = pressure%slab(i,j)
            
            specificHumidityKgKg = specificHumidity%slab(i,j)
            if (specificHumidityKgKg > 9.9999999e+10) then
               relativeHumidity%slab(i,j) = 1.e+15
               cycle
            end if

            mixingRatioKgKg = calcMixingRatio(specificHumidityKgKg)
            satSpecificHumidityKgKg = &
                 calcSatSpecificHumidity(temperatureK, pressurePa)
            satMixingRatioKgKg = calcMixingRatio(satSpecificHumidityKgKg)
            
            relativeHumidityPct = mixingRatioKgKg/satMixingRatioKgKg*100.
            relativeHumidityPct = max(0.,relativeHumidityPct)
            relativeHumidityPct = min(100.,relativeHumidityPct)

            relativeHumidity%slab(i,j) = relativeHumidityPct

         end do
      end do

      relativeHumidities(pressureCount) = copyFieldWPS(relativeHumidity)

      ! Clean up slabs
      call destroyFieldWPS(relativeHumidity)
      call destroyFieldWPS(specificHumidity)
      call destroyFieldWPS(temperature)
      call destroyFieldWPS(pressure)

      ! Go to next press
      pressureCount=pressureCount + 1

   end do

   ! Handle 2-m relative humidity
   if (processSurfacePressure) then

      ! Open new WPS file
      fileUnit2=select_file_unit()
      outputFile=createFileWPS(fileUnit2,trim(directory), &
           'RH_2M_ABOVE_GROUND_LEVEL',year,month,day,hour,minute,second, &
           includeMinute,includeSecond)
      
      ! Write 2M relative humidity to file
      pressureCount=1
      do
         if (pressureCount > numberOfSlabs) exit
         if (relativeHumidities(pressureCount)%xlvl == SURFACE_XLVL) then
            call writeFileWPS(outputFile,relativeHumidities(pressureCount))
         end if
         pressureCount = pressureCount + 1
      end do
      
      ! Close 2M file
      call destroyFileWPS(outputFile)
      call close_file_unit(fileUnit2)

   end if

   ! Open new WPS file
   fileUnit2=select_file_unit()
   if (onIsobaricLevels) then
      outputFile=createFileWPS(fileUnit2,trim(directory), &
           'RH_ISOBARIC_LEVEL',year,month,day,hour,minute,second, &
           includeMinute,includeSecond)
   else
      outputFile=createFileWPS(fileUnit2,trim(directory), &
           'RH_MODEL_LEVEL',year,month,day,hour,minute,second, &
           includeMinute,includeSecond)
   end if

   ! Write relative humidities to new file
   pressureCount=1
   do
      if (pressureCount > numberOfSlabs) exit
      if (relativeHumidities(pressureCount)%xlvl .ne. SURFACE_XLVL) then
         call writeFileWPS(outputFile,relativeHumidities(pressureCount))
      end if
      pressureCount = pressureCount + 1
   end do

   ! Clean up
   call destroyFileWPS(outputFile)
   call close_file_unit(fileUnit2)

   do i = 1, numberOfSlabs
      call destroyFieldWPS(relativeHumidities(i))
   end do
   deallocate(relativeHumidities)
   call destroyFieldWPSList(pressureList)
   call destroyFieldWPSList(temperatureList)
   call destroyFieldWPSList(specificHumidityList)

end program createRH
