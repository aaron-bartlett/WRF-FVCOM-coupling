!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  createLANDSEA
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Back-end utility program for GEOS2WRF package.  Reads in WPS intermediate
! format file with lake and ocean fractions, and calculates land-sea mask.
! This new field is then written to a new WPS file.
!
! REVISION:
! 1 May 2012 - First version.
!
!------------------------------------------------------------------------------

program createLANDSEA

   ! Import modules
   use FieldWPS_mod
   use FileUnit_mod
   use FileWPS_mod

   ! Change defaults
   implicit none

   ! Local variables
   type(FieldWPS) :: field, lakeFraction, oceanFraction, landSea
   type(FileWPS) :: inputFile, outputFile
   logical :: isCorrectTime
   logical :: foundLakeFraction, foundOceanFraction
   logical :: endOfFile
   integer :: fileUnit1,fileUnit2
   real :: waterFraction
   integer :: i,j

    ! Namelist variables
   character(len=132) :: directory
   character(len=64) :: prefix
   integer :: year,month,day,hour,minute,second
   logical :: includeMinute,includeSecond
   character(len=64) :: lakeFractionName
   character(len=64) :: oceanFractionName
   namelist /input/ directory,prefix,year,month,day,hour,minute,second, &
        includeMinute,includeSecond, &
        lakeFractionName,oceanFractionName

   ! Internal constants
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
   lakeFractionName='NULL'
   oceanFractionName='NULL'
   fileUnit1=select_file_unit()
   open(unit=fileUnit1,file='namelist.createLANDSEA',delim='APOSTROPHE')
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
   if (trim(lakeFractionName) .eq. 'NULL') then
      print*,'ERROR assigning value to lakeFractionName in namelist!'
      stop 1
   end if
   if (trim(oceanFractionName) .eq. 'NULL') then
      print*,'ERROR assigning value to oceanFractionName in namelist!'
      stop 1
   end if
   
   ! Open the WPS file
   fileUnit1=select_file_unit()
   inputFile=createFileWPS(fileUnit1,trim(directory),trim(prefix),year,month, &
        day,hour,minute,second,includeMinute,includeSecond,preserve=.true.)
     
   ! Loop through each slab looking for lake and ocean fractions.
   foundLakeFraction=.false.
   foundOceanFraction=.false.
   do
      ! Read the next slab
      call readFileWPS(inputFile,field,endOfFile)
      if (endOfFile) exit

      ! Skip if wrong time
      isCorrectTime = checkDateTimeFieldWPS(field,year,month,day,hour, &
           minute,second)
      if (.not. isCorrectTime) then
         cycle
      end if

      ! Skip if not surface layer (or 2-m AGL)
      if (field%xlvl .ne. SURFACE_XLVL) then
         call destroyFieldWPS(field)
         cycle
      end if

      if (trim(field%field) .eq. trim(lakeFractionName)) then
         lakeFraction = copyFieldWPS(field)
         foundLakeFraction = .true.
      else if (trim(field%field) .eq. trim(oceanFractionName)) then
         oceanFraction = copyFieldWPS(field)
         foundOceanFraction = .true.
      end if
      call destroyFieldWPS(field)

      ! All done.
      if (foundLakeFraction .and. &
          foundOceanFraction) then
         exit
      end if

   end do
   call destroyFileWPS(inputFile)
   call close_file_unit(fileUnit1)

   ! Sanity checks
   if (.not. foundLakeFraction) then
      print*,'ERROR, did not find lake fraction!'
      print*,'Looked for variable ',trim(lakeFractionName)
      print*,'and level ',SURFACE_XLVL
      stop 1
   end if
   if (.not. foundOceanFraction) then
      print*,'ERROR, did not find ocean fraction!'
      print*,'Looked for variable ',trim(oceanFractionName)
      print*,'and level ',SURFACE_XLVL
      stop 1
   end if
   
   ! Make sure slabs are for the correct time
!   call checkDateTimeFieldWPS(lakeFraction,year,month,day,hour,minute,second)
!   call checkDateTimeFieldWPS(oceanFraction,year,month,day,hour,minute,second)

   ! Make sure grids match for the two input fields
   call compareGridsFieldWPS(lakeFraction,oceanFraction)

   ! Create landsea mask field
   landSea = copyFieldWPS(lakeFraction)
   landSea%field="LANDSEA"
   landSea%units="proprtn"
   landSea%desc="Land/Sea flag (1=land, 0 or 2=sea)"
   do j = 1, landSea%ny
      do i = 1, landSea%nx
         waterFraction = lakeFraction%slab(i,j) + oceanFraction%slab(i,j)
         if (waterFraction < 0.5) then
            landSea%slab(i,j) = 1.
         else
            landSea%slab(i,j) = 0.
         end if
      end do
   end do

   ! Clean up input data
   call destroyFieldWPS(lakeFraction)
   call destroyFieldWPS(oceanFraction)
   
   ! Write to new WPS file
   fileUnit2=select_file_unit()
   outputFile=createFileWPS(fileUnit2,trim(directory),'LANDSEA_GROUND_LEVEL',&
        year,month,day,hour,minute,second,includeMinute,includeSecond)
   call writeFileWPS(outputFile,landSea)
   call destroyFieldWPS(landSea)
   call destroyFileWPS(outputFile)
   call close_file_unit(fileUnit2)

end program createLANDSEA
