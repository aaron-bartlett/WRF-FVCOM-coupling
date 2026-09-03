!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  createSOILHGT
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! Back-end utility program for GEOS2WRF package.  Reads in WPS intermediate
! format file with surface geopotential field, derives terrain field, and
! writes terrain in new WPS file.
!
! REVISION:
! 18 Apr 2012 - First version.
!
!------------------------------------------------------------------------------

program createSOILHGT

   ! Import modules
   use FieldWPS_mod
   use FileUnit_mod
   use FileWPS_mod

   ! Change defaults
   implicit none

   ! Local variables
   type(FieldWPS) :: field,soilhgt
   type(FileWPS) :: inputFile, outputFile
   logical :: done
   logical :: endOfFile
   integer :: fileUnit1,fileUnit2
   character(len=132) :: directory
   character(len=64) :: prefix
   integer :: year,month,day,hour,minute,second
   character(len=64) :: surfaceGeopotentialName
   logical :: includeMinute,includeSecond
   namelist /input/ directory,prefix,year,month,day,hour,minute,second, &
        includeMinute,includeSecond,surfaceGeopotentialName
   
   integer :: i,j

   real, parameter :: INVERSE_G_0 = 1./9.80665 ! Gravity at mean sea level
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
   surfaceGeopotentialName='NULL'
   fileUnit1=select_file_unit()
   open(unit=fileUnit1,file='namelist.createSOILHGT',delim='APOSTROPHE')
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

   if (trim(surfaceGeopotentialName) .eq. 'NULL') then
      print*,'ERROR assigning value to surfaceGeopotentialName in namelist!'
      stop 1
   end if

   ! Open the file
   fileUnit1=select_file_unit()
   inputFile=createFileWPS(fileUnit1,trim(directory),trim(prefix),year,month, &
        day,hour,minute,second,includeMinute,includeSecond,preserve=.true.)

   ! Loop through each slab looking for surface geopotential
   done = .false.
   do
      ! Read the next slab
      call readFileWPS(inputFile,field,endOfFile)
      if (endOfFile) exit

      ! Check the variable and level
      if (trim(field%field) .ne. trim(surfaceGeopotentialName) .or. &
          field%xlvl .ne. SURFACE_XLVL) then
         call destroyFieldWPS(field)
         cycle
      end if

      ! Copy the field, then update the values, units, etc.
      soilhgt = copyFieldWPS(field)
      do j = 1, field%ny
         do i = 1, field%nx
            soilhgt%slab(i,j) = INVERSE_G_0 * field%slab(i,j)
         end do
      end do
      soilhgt%field = "SOILHGT"
      soilhgt%units = "m"
      soilhgt%desc = "Terrain field of source analysis"

      call destroyFieldWPS(field)
      done = .true.
      exit

   end do
   call destroyFileWPS(inputFile)
   call close_file_unit(fileUnit1)

   ! Sanity check
   if (.not. done) then
      print*,'ERROR, did not find surface geopotential!'
      print*,'Looked for variable ',trim(surfaceGeopotentialName)
      print*,'and level ',SURFACE_XLVL
      stop 1
   end if

   ! Write to new WPS file
   fileUnit2=select_file_unit()
   outputFile=createFileWPS(fileUnit2,trim(directory),'SOILHGT_GROUND_LEVEL',&
        year,month,day,hour,minute,second,includeMinute,includeSecond)
   call writeFileWPS(outputFile,soilhgt)
   call destroyFieldWPS(soilhgt)
   call destroyFileWPS(outputFile)
   call close_file_unit(fileUnit2)

end program createSOILHGT
