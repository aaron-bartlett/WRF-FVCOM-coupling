!------------------------------------------------------------------------------
! NASA/GSFC, Computational and Information Science and Technology Office,
! Code 606
!------------------------------------------------------------------------------
!
! PROGRAM:  temporalInterpolation
!
! AUTHOR:
! Eric Kemp, NASA CISTO/SSAI
!
! DESCRIPTION:
! Back-end utility program for GEOS2WRF package.  Reads in two WPS intermediate
! format files valid at different times, performs linear temporal
! interpolation, and writes new WPS file.
!
! REVISION:
! 20 May 2015 - First version
!
!------------------------------------------------------------------------------

program temporalInterpolation

   ! Import modules
   use, intrinsic :: iso_c_binding, only: c_int
   use DerivedVars_mod, only: bubbleSort
   use FieldWPS_mod
   use FieldWPSList_mod
   use FileUnit_mod
   use FileWPS_mod

   ! Defaults
   implicit none

   ! Wrapper for C function
   interface
      function calc_epoch_sec(year,month,day,hour,minute,second) bind (c)
         import :: c_int
         integer(kind=c_int) :: year,month,day
         integer(kind=c_int) :: hour,minute,second
         integer(kind=c_int) :: calc_epoch_sec ! Return variable
      end function calc_epoch_sec
   end interface
      
   ! Local variables
   integer :: minute,second
   type(FieldWPSList) :: var1List, var2List
   real, allocatable :: var1(:,:,:),var2(:,:,:),varOutput(:,:,:)
   character(len=24) :: hdateOutput
   type(FieldWPS) :: firstField1,firstField2
   type(FileWPS) :: inputFile1, inputFile2, outputFile
   integer :: var1Count,var2Count
   integer :: fileUnit1, fileUnit2
   integer :: epochSec1,epochSec2,epochSecOutput
   real,allocatable :: xlvls1(:),xlvls2(:)
   real :: x,x1,x2,weight
   integer :: i,j,k

   ! Common namelist
   character(len=9) :: fieldName
   namelist /all/ fieldName

   ! Namelist for first input file
   character(len=132) :: directory1
   character(len=64) :: prefix1
   integer :: year1,month1,day1,hour1
   namelist /input1/ directory1,prefix1,year1,month1,day1, &
        hour1

   ! Namelist for second input file
   character(len=132) :: directory2
   character(len=64) :: prefix2
   integer :: year2,month2,day2,hour2
   namelist /input2/ directory2,prefix2,year2,month2,day2, &
        hour2

   ! Namelist for output file
   character(len=132) :: directoryOutput
   integer :: yearOutput,monthOutput,dayOutput, &
        hourOutput
   namelist /output/ directoryOutput, &
        yearOutput,monthOutput,dayOutput, &
        hourOutput
   
   ! Initialize namelist variables
   fieldName='NULL'
   directory1='NULL' ; directory2='NULL' ; directoryOutput='NULL'
   prefix1='NULL'    ; prefix2='NULL'    ; 
   year1=0000        ; year2=0000        ; yearOutput=0000
   month1=00         ; month2=00         ; monthOutput=00
   day1=00           ; day2=00           ; dayOutput=00
   hour1=-1          ; hour2=-1          ; hourOutput=-1


   ! Read namelist file
   fileUnit1=select_file_unit()
   open(unit=fileUnit1,file='namelist.temporalInterpolation', &
        delim='APOSTROPHE')
   read(unit=fileUnit1,nml=all)
   read(unit=fileUnit1,nml=input1)
   read(unit=fileUnit1,nml=input2)
   read(unit=fileUnit1,nml=output)
   call close_file_unit(fileUnit1)

   ! Sanity check namelist settings
   if (trim(fieldName) .eq. 'NULL') then
      print*,'ERROR assigning value to fieldName in namelist!'
      stop 1
   end if

   if (trim(directory1) .eq. 'NULL') then
      print*,'ERROR assigning value to directory1 in namelist!'
      stop 1
   end if
   if (trim(directory2) .eq. 'NULL') then
      print*,'ERROR assigning value to directory2 in namelist!'
      stop 1
   end if
   if (trim(directoryOutput) .eq. 'NULL') then
      print*,'ERROR assigning value to directoryOutput in namelist!'
      stop 1
   end if

   if (trim(prefix1) .eq. 'NULL') then
      print*,'ERROR assigning value to prefix1 in namelist!'
      stop 1
   end if
   if (trim(prefix2) .eq. 'NULL') then
      print*,'ERROR assigning value to prefix2 in namelist!'
      stop 1
   end if

   if (year1 .eq. 0000 .or. year1 < 0) then
      print*,'ERROR assigning valid value to year1 in namelist!'
      stop 1
   end if
   if (year2 .eq. 0000 .or. year2 < 0) then
      print*,'ERROR assigning valid value to year2 in namelist!'
      stop 1
   end if
   if (yearOutput .eq. 0000 .or. yearOutput < 0) then
      print*,'ERROR assigning valid value to yearOutput in namelist!'
      stop 1
   end if

   if (month1 .eq. 00 .or. month1 < 0 .or. month1 > 12) then
      print*,'ERROR assigning valid value to month1 in namelist!'
      stop 1
   end if
   if (month2 .eq. 00 .or. month2 < 0 .or. month2 > 12) then
      print*,'ERROR assigning valid value to month2 in namelist!'
      stop 1
   end if
   if (monthOutput .eq. 00 .or. monthOutput < 0 .or. monthOutput > 12) then
      print*,'ERROR assigning valid value to monthOutput in namelist!'
      stop 1
   end if

   if (day1 .eq. 00 .or. day1 < 0 .or. day1 > 31) then
      print*,'ERROR assigning valid value to day1 in namelist!'
      stop 1
   end if
   if (day2 .eq. 00 .or. day2 < 0 .or. day2 > 31) then
      print*,'ERROR assigning valid value to day2 in namelist!'
      stop 1
   end if
   if (dayOutput .eq. 00 .or. dayOutput < 0 .or. dayOutput > 31) then
      print*,'ERROR assigning valid value to dayOutput in namelist!'
      stop 1
   end if

   if (hour1 .eq. -1 .or. hour1 < 0 .or. hour1 > 23) then
      print*,'ERROR assigning valid value to hour1 in namelist!'
      print*,'hour1 = ',hour1
      stop 1
   end if
   if (hour2 .eq. -1 .or. hour2 < 0 .or. hour2 > 23) then
      print*,'ERROR assigning valid value to hour2 in namelist!'
      stop 1
   end if
   if (hourOutput .eq. -1 .or. hourOutput < 0 .or. hourOutput > 23) then
      print*,'ERROR assigning valid value to hourOutput in namelist!'
      stop 1
   end if

   ! Get epoch seconds of two input date/times and that of output
   ! WPS intermediate files assume data valid at top of hour.
   minute=0
   second=0
   epochSec1 = calc_epoch_sec(year1,month1,day1,hour1,minute,second)
   epochSec2 = calc_epoch_sec(year2,month2,day2,hour2,minute,second)
   epochSecOutput = calc_epoch_sec(yearOutput,monthOutput, &
        dayOutput,hourOutput,minute,second)

   ! Check limits of epochSec values
   if (.not. (epochSec1 .le. epochSec2)) then
      print*,'ERROR, date/time for first input file is not before second file!'
      stop 1
   end if
   if (.not. (epochSec1 .le. epochSecOutput) .or. &
       .not. (epochSecOutput .le. epochSec2)) then
      print*,'ERROR, date/time of output file not between input files!'
      stop 1
   end if

   ! Create time string for interpolated data.
   write(hdateOutput,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
        yearOutput,":",monthOutput,":",dayOutput,"_", &
        hourOutput,":",minute,":",second

   ! Open the WPS input files.
   fileUnit1=select_file_unit()
   inputFile1=createFileWPS(fileUnit1,trim(directory1),trim(prefix1), &
        year1,month1,day1,hour1,minute,second,preserve=.true.)
   fileUnit2=select_file_unit()
   inputFile2=createFileWPS(fileUnit2,trim(directory2),trim(prefix2), &
        year2,month2,day2,hour2,minute,second,preserve=.true.)
   
   ! Get the slabs of the specified variable, plus a count of the valid slabs 
   ! and sample first slabs. Note that the getListOfSlabs subroutine will also
   ! check against the requested valid time, and check for differences in map 
   ! projection or grid dimensions.
   call getListOfSlabs(inputFile1,fieldName, &
        year1,month1,day1,hour1,minute,second, &
        var1List,var1Count,firstField1)
   call destroyFileWPS(inputFile1)
   call close_file_unit(fileUnit1)

   call getListOfSlabs(inputFile2,fieldName, &
        year2,month2,day2,hour2,minute,second, &
        var2List,var2Count,firstField2)
   call destroyFileWPS(inputFile2)
   call close_file_unit(fileUnit2)
        
   if (var1Count .ne. var2Count) then
      print*,'ERROR, total number of slabs differ between input files!'
      print*,'First has ',var1Count
      print*,'Second has ',var2Count
      stop 1
   end if

   ! Collect and sort slab levels
   allocate(xlvls1(var1Count))
   call getXlvlsFieldWPSList(var1List,var1Count,xlvls1)
   call bubbleSort(var1Count,xlvls1)
   allocate(xlvls2(var2Count))
   call getXlvlsFieldWPSList(var2List,var2Count,xlvls2)
   call bubbleSort(var2Count,xlvls2)
   
   ! Make sure slab levels are identical
   do i = 1,var1Count
      if (xlvls1(i) .ne. xlvls2(i)) then
         print*,'ERROR, slab levels do not match between input files!'
         print*,'i, xlvl1, xlvl2: ',i,xlvls1(i),xlvls2(i)
         stop 1
      end if
   end do

   ! Make sure grids in first input file match those of second.
   call compareGridsFieldWPS(firstField1,firstField2)

   ! Collect the slabs into arrays.
   allocate(var1(firstField1%nx,firstField1%ny,var1Count))
   call copyListToArray(firstField1%nx,firstField1%ny,var1Count, &
        var1List,fieldName,xlvls1,var1)
   allocate(var2(firstField2%nx,firstField2%ny,var2Count))
   call copyListToArray(firstField2%nx,firstField2%ny,var2Count, &
        var2List,fieldName,xlvls1,var2)

   ! Now linearly interpolate in time.
   x1 = real(epochSec1)
   x2 = real(epochSec2)
   x  = real(epochSecOutput)
   weight = (x - x1)/(x2 - x1)
   allocate(varOutput(firstField1%nx,firstField1%ny,var1Count))
   do k = 1,var1Count
      do j = 1,firstField1%ny
         do i = 1,firstField1%nx
            varOutput(i,j,k) = var1(i,j,k) + &
                 (var2(i,j,k) - var1(i,j,k))*weight
         end do
      end do
   end do
   deallocate(var2)
   deallocate(var1)

   ! Now write out new slabs
   fileUnit1 = select_file_unit()
   outputFile=createFileWPS(fileUnit1,trim(directoryOutput), &
        trim(fieldName),yearOutput,monthOutput,dayOutput, &
        hourOutput,minute,second)
   call writeNewSlabs(firstField1%nx,firstField1%ny,var1Count, &
        var1List,trim(fieldName),xlvls1,varOutput,hdateOutput, &
        outputFile)

   ! Close output file
   call destroyFileWPS(outputFile)
   call close_file_unit(fileUnit1)

   ! Clean up
   call destroyFieldWPS(firstField1)
   call destroyFieldWPS(firstField2)
   call destroyFieldWPSList(var1List)
   call destroyFieldWPSList(var2List)
   deallocate(varOutput)
   deallocate(xlvls2)
   deallocate(xlvls1)

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  getListofSlabs
   !
   ! DESCRIPTION:  Internal subroutine to pull slab fields from file and put
   ! slabs for a particular variable into a linked list. Also returns count of
   ! fields in linked list, and copy of the first field put in the list. Also
   ! checks fields as they are read to ensure valid time is correct, and that
   ! the grid dimensions and map projections are consistent.
   !
   !---------------------------------------------------------------------------

   subroutine getListOfSlabs(inputFile,fieldName, &
        year,month,day,hour,minute,second, &
        varList,varCount,firstField) 

      ! Defaults
      implicit none

      ! Arguments
      type(FileWPS),intent(in) :: inputFile
      character(len=9),intent(in) :: fieldName
      integer,intent(in) :: year,month,day,hour,minute,second
      type(FieldWPSList),intent(out) :: varList
      integer,intent(out) :: varCount
      type(fieldWPS),intent(out) :: firstField

      ! Local variables
      logical :: firstSlab
      logical :: endOfFile
      logical :: isCorrectTime
      type(FieldWPS) :: field

      firstSlab=.true.
      varCount=0
      do         
         ! Read the next slab
         call readFileWPS(inputFile,field,endOfFile)
         if (endOfFile) exit

         ! Skip if wrong variable
         if (trim(field%field) .ne. trim(fieldName)) then
            call destroyFieldWPS(field)
            cycle
         end if

         ! Skip if wrong date/time
         isCorrectTime=checkDateTimeFieldWPS(field,year,month,day,hour, &
              minute,second)
         if (.not. isCorrectTime) then
            call destroyFieldWPS(field)
            cycle
         end if

         ! If this is the first slab, make a copy.  Otherwise, compare the 
         ! grids between the current slab and the first one saved.
         if (firstSlab) then
            firstField=copyFieldWPS(field)
         else
            call compareGridsFieldWPS(firstField,field)
         end if

         varCount = varCount + 1
         call appendFieldWPSList(varList,field)
         call destroyFieldWPS(field)
      end do

      ! Sanity check
      if (varCount .eq. 0) then
         print*,'ERROR, no slabs found for variable ',trim(fieldName)
         stop 1
      end if

      return
   end subroutine getListOfSlabs

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  copyListToArray
   !
   ! DESCRIPTION:  Internal subroutine to pull slabs from linked list and 
   ! copy gridded values to 3D array.
   !
   !---------------------------------------------------------------------------

   subroutine copyListToArray(nx,ny,nz,fieldList,fieldName,xlvls,array3d)

      ! Arguments
      integer,intent(in) :: nx,ny,nz
      type(FieldWPSList),intent(in) :: fieldList
      character(len=*),intent(in) :: fieldName
      real,intent(in) :: xlvls(nz)
      real,intent(out) :: array3d(nx,ny,nz)

      ! Local variables
      type(FieldWPS) :: field
      integer :: i,j,k

      do k = 1, nz
         field=findFieldWPSList(fieldList, &
              trim(fieldName),xlvls(k))
         do j = 1, ny
            do i = 1, nx
               array3d(i,j,k) = field%slab(i,j)
            end do
         end do
         call destroyFieldWPS(field)
      end do
      
      return
   end subroutine copyListToArray
   
   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  writeNewSlabs
   !
   ! DESCRIPTION:  Internal subroutine to write new set of slabs in an 
   ! already opened WPS intermediate file.
   !
   !---------------------------------------------------------------------------

   subroutine writeNewSlabs(nx,ny,numberOfSlabs,fieldList,fieldName, &
        xlvls,array,hdateOutput,outputFile)

      implicit none

      ! Arguments
      integer,intent(in) :: nx,ny,numberOfSlabs
      type(FieldWPSList),intent(in) :: fieldList
      character(len=*),intent(in) :: fieldName
      real,intent(in) :: xlvls(numberOfSlabs)
      real,intent(in) :: array(nx,ny,numberOfSlabs)
      character(len=*),intent(in) :: hdateOutput
      type(FileWPS),intent(in) :: outputFile

      ! Local variables
      type(FieldWPS) :: field
      integer :: i,j,k

      do k = 1,numberOfSlabs
         field = findFieldWPSList(fieldList,trim(fieldName),xlvls(k))
         field%hdate = trim(hdateOutput)
         field%field = trim(fieldName)
         do j = 1,ny
            do i = 1,nx
               field%slab(i,j) = array(i,j,k)
            end do
         end do
         call writeFileWPS(outputFile,field)
         call destroyFieldWPS(field)
      end do

      return
   end subroutine writeNewSlabs
end program temporalInterpolation
