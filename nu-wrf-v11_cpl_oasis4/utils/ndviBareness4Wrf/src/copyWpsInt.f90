!------------------------------------------------------------------------------
! NASA/GSFC, Computational and Information Sciences and Technology Office,
! Code 606
!------------------------------------------------------------------------------
!
! PROGRAM:  copyWpsInt
!
! AUTHOR:
! Eric Kemp, NASA CISTO/SSAI
!
! DESCRIPTION:
! Copies data in a WPS intermediate file to another WPS intermediate file with
! different date/time.
!
! REVISION:
! 27 July 2016 - First version.
!
!------------------------------------------------------------------------------

program copyWpsInt
   
   ! Modules
   use FieldWPS_mod
   use FileUnit_mod
   use FileWPS_mod

   ! Defaults
   implicit none

   ! Local variables
   character(len=105) :: dataDirectory
   character(len=14) :: wpsPrefix
   character(len=13) :: inputDateTime
   character(len=13) :: outputDateTime
   namelist /settings/ dataDirectory,wpsPrefix,inputDateTime,outputDateTime

   type(FileWPS) :: inputWpsFile,outputWpsFile
   type(FieldWPS) :: field
   integer :: inputFileUnit,outputFileUnit
   integer :: inputYear,inputMonth,inputDay,inputHour
   integer :: outputYear,outputMonth,outputDay,outputHour
   logical :: endOfFile
   character(len=24) :: hdate
   
   ! Initialize namelist variables
   dataDirectory = 'NULL'
   wpsPrefix = 'NULL'
   inputDateTime  = '0000-00-00_00'
   outputDateTime = '0000-00-00_00'

   ! Read namelist file
   inputFileUnit=selectFileUnit()
   open(unit=inputFileUnit,file='namelist.copyWpsInt', &
        delim='APOSTROPHE')
   read(unit=inputFileUnit,nml=settings)
   call closeFileUnit(inputFileUnit)

   ! Sanity check namelist settings
   if (trim(dataDirectory) .eq. 'NULL') then
      print*,'FATAL, no value assigned to dataDirectory in namelist!'
      stop 1
   end if
   if (trim(wpsPrefix) .eq. 'NULL') then
      print*,'FATAL, no value assigned to wpsPrefix in namelist!'
      stop 1
   end if
   if (trim(inputDateTime) .eq. '0000-00-00_00') then
      print*,'FATAL, no value assigned to inputDateTime in namelist!'
      stop 1
   end if
   if (trim(outputDateTime) .eq. '0000-00-00_00') then
      print*,'FATAL, no value assigned to outputDateTime in namelist!'
      stop 1
   end if
   
   ! Open input WPS file without overwriting it.
   read(inputDateTime,'(I4.4,X,I2.2,X,I2.2,X,I2.2)') &
        inputYear,inputMonth,inputDay,inputHour
   inputFileUnit=selectFileUnit()
   inputWpsFile=createFileWPS(inputFileUnit,trim(dataDirectory), &
        trim(wpsPrefix),inputYear,inputMonth,inputDay,inputHour, &
        preserve=.true.)

   ! Open new output WPS file.
   read(outputDateTime,'(I4.4,X,I2.2,X,I2.2,X,I2.2)') &
        outputYear,outputMonth,outputDay,outputHour
   outputFileUnit=selectFileUnit()
   outputWpsFile=createFileWPS(outputFileUnit,trim(dataDirectory), &
        trim(wpsPrefix),outputYear,outputMonth,outputDay,outputHour, &
        preserve=.false.)

   ! Loop through each field in the input WPS file, and copy to output WPS
   ! file with changed date/time.
   write(hdate,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A)') &
        outputYear,':',outputMonth,':',outputDay,'_',outputHour,':00:00'
   do
      call readFileWPS(inputWpsFile,field,endOfFile)
      if (endOfFile) exit      
      field%hdate = hdate
      call writeFileWPS(outputWpsFile,field)      
      call destroyFieldWPS(field)
   end do

   ! Clean up
   call destroyFileWPS(outputWpsFile)
   call closeFileUnit(outputFileUnit)
   call destroyFileWPS(inputWpsFile)
   call closeFileUnit(inputFileUnit)
   call destroyFieldWPS(field)
   
end program copyWpsInt
