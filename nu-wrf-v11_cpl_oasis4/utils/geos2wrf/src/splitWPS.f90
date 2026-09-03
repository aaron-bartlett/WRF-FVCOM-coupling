!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  splitWPS
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! Reads in WPS intermediate file, loops through each 2D slab, and writes each
! slab in new files.
!
! REVISION:
! 17 Apr 2012 - First version.
!------------------------------------------------------------------------------

program splitWPS

   ! Import modules
   use FieldWPS_mod
   use FileUnit_mod
   use FileWPS_mod
   
   ! Change defaults
   implicit none

   ! Local variables
   integer :: counter
   logical :: endOfFile
   integer :: fileUnit1,fileUnit2
   type(FieldWPS) :: field
   type(FileWPS) :: file1,file2
   integer :: int_xlvl

   ! Namelist variables
   character(len=132) :: directory
   character(len=64) :: prefix
   integer :: year,month,day,hour
   namelist /input/ directory,prefix,year,month,day,hour

   ! Read namelist file
   directory='NULL'
   prefix='NULL'
   year=0000
   month=00
   day=00
   hour=00
   fileUnit1=select_file_unit()
   open(unit=fileUnit1,file='namelist.splitWPS',delim='APOSTROPHE')
   read(unit=fileUnit1,nml=input)
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
      print*,'ERROR assigning value to year in namelist!'
      stop 1
   end if
   if (day .eq. 0) then
      print*,'ERROR assigning value to year in namelist!'
      stop 1
   end if

   ! Open the file
   fileUnit1=select_file_unit()
   file1=createFileWPS(fileUnit1,trim(directory),trim(prefix),year,month, &
        day,hour,preserve=.true.)

   ! Loop through each slab and write to new files
   counter = 0
   do 
      ! Read the next slab
      call readFileWPS(file1,field,endOfFile)
      if (endOfFile) exit

      ! Find the day and time
      read(field%hdate,'(i4.4,1x,i2.2,1x,i2.2,1x,i2.2)') year,month,day,hour
      
      ! Write new single-slab WPS file
      counter = counter + 1
!      write(prefix,'(A,i4.4)') 'MERRA_',counter
      int_xlvl = field%xlvl
      write(prefix,'(A,A,I6.6)') trim(field%field),'_',int_xlvl
      fileUnit2=select_file_unit()
      file2=createFileWPS(fileUnit2,trim(directory),trim(prefix), &
           year,month,day,hour)
      call writeFileWPS(file2,field)
      call destroyFieldWPS(field)
      call destroyFileWPS(file2)
      call close_file_unit(fileUnit2)
   end do

   ! Clean up
   call destroyFileWPS(file1)
   call close_file_unit(fileUnit1)

end program splitWPS
