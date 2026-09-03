!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  extrapIsobaric
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Back-end utility program for GEOS2WRF package.  Reads in WPS intermediate
! format file with 3D geopotential height, temperature, relative humidity,
! U and V winds, and extrapolate data to levels that are underground.
! A new WPS file with the updated data will then be written.
!
! REVISION:
! 18 May 2012 - First version
!
!------------------------------------------------------------------------------

program extrapIsobaric

   ! Import modules
   use DerivedVars_mod
   use FieldWPS_mod
   use FieldWPSList_mod
   use FileUnit_mod
   use FileWPS_mod
   use WrfRH_mod

   ! Defaults
   implicit none

   ! Local variables
   type(FieldWPSList) :: geopotentialHeightList, temperatureList, &
        relativeHumidityList, uList, vList
   real, allocatable :: heights(:,:,:), temperatures(:,:,:), &
        relativeHumidities(:,:,:), uWinds(:,:,:), vWinds(:,:,:)
   logical, allocatable :: missing3d(:,:,:), missing1d(:)
   real, allocatable :: pressureXlvls(:)
   type(FieldWPS) :: fieldFirst, field
   type(FileWPS) :: inputFile, outputFile
   integer :: geopotentialHeightCount, temperatureCount, &
        relativeHumidityCount,uCount,vCount
   integer :: numberOfSlabs
   logical :: endOfFile
   integer :: fileUnit1,fileUnit2
   logical :: isCorrectTime
   logical :: firstSlab
   integer :: i,j,k
   real :: zAbove,pAbove,tAbove
   real :: pBelow
   real :: lnPRatio
   real :: RHAbove,qsatAbove,qvAbove
   real :: tvAbove,tvBelow,meanTv
   real :: part, gamma_s

   ! Namelist variables
   character(len=132) :: directory
   character(len=64) :: prefix
   integer :: year,month,day,hour,minute,second
   character(len=64) :: geopotentialHeightName
   character(len=64) :: temperatureName
   character(len=64) :: relativeHumidityName
   character(len=64) :: uName
   character(len=64) :: vName
   namelist /input/ directory,prefix,year,month,day,hour,minute,second, &
        geopotentialHeightName,temperatureName,relativeHumidityName, &
        uName, vName

   ! Constants
   real,parameter :: MISSING=9.9999999e+10
   real,parameter :: SURFACE_XLVL = 200100.
   real,parameter :: MEAN_SEA_LEVEL_XLVL = 201300.
   real,parameter :: HALF_GAMMA = 6.5*0.001/2. ! Half assumed lapse rate
   real,parameter :: GAMMA = 6.5*0.001 ! 6.5 K per km
   real,parameter :: R_OVR_G = 287. / 9.81 ! Dry gas constant over gravity
   real,parameter :: G_OVR_R = 9.81/287.
   real,parameter :: ZSHUL=75.
   real,parameter :: TVSHUL=290.66

   ! Read namelist file
   directory='NULL'
   prefix='NULL'
   year=0000
   month=00
   day=00
   hour=00
   minute=00
   second=00
   geopotentialHeightName='NULL'
   temperatureName='NULL'
   relativeHumidityName='NULL'
   uName='NULL'
   vName='NULL'
   fileUnit1=select_file_unit()
   open(unit=fileUnit1,file='namelist.extrapIsobaric',delim='APOSTROPHE')
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
   if (minute < 0 .or. minute > 23) then
      print*,'ERROR assigning valid value to minute in namelist!'
      stop 1
   end if
   if (second < 0 .or. second > 23) then
      print*,'ERROR assigning valid value to second in namelist!'
      stop 1
   end if
   if (trim(geopotentialHeightName) .eq. 'NULL') then
      print*,'ERROR assigning value to geopotentialHeightName in namelist!'
      stop 1
   end if
   if (trim(temperatureName) .eq. 'NULL') then
      print*,'ERROR assigning value to temperatureName in namelist!'
      stop 1
   end if
   if (trim(relativeHumidityName) .eq. 'NULL') then
      print*,'ERROR assigning value to relativeHumidityName in namelist!'
      stop 1
   end if
   if (trim(uName) .eq. 'NULL') then
      print*,'ERROR assigning value to uName in namelist!'
      stop 1
   end if
   if (trim(vName) .eq. 'NULL') then
      print*,'ERROR assigning value to vName in namelist!'
      stop 1
   end if

   ! Open the GEOS WPS file
   fileUnit1=select_file_unit()
   inputFile=createFileWPS(fileUnit1,trim(directory),trim(prefix),year,month, &
        day,hour,minute,second,preserve=.true.)

   ! Loop through the GEOS WPS file looking for geopotential height, 
   ! temperature, relative humidity, u and v slabs.  We will add these slabs 
   ! to linked lists as they are encountered, and count them.
   geopotentialHeightCount=0
   temperatureCount=0
   relativeHumidityCount=0
   uCount=0
   vCount=0
   geopotentialHeightList=createFieldWPSList()
   temperatureList=createFieldWPSList()
   relativeHumidityList=createFieldWPSList()
   uList=createFieldWPSList()
   vList=createFieldWPSList()
   firstSlab=.true.
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

      ! Skip mean sea level and (near) ground level data.
      if (field%xlvl .eq. SURFACE_XLVL .or. &
          field%xlvl .eq. MEAN_SEA_LEVEL_XLVL) then
         call destroyFieldWPS(field)
         cycle
      end if

      ! If this is the first slab, make a copy.  Otherwise, compare the grids
      ! between the current slab and the first one saved.
      if (firstSlab) then
         fieldFirst=copyFieldWPS(field)
      else
         call compareGridsFieldWPS(fieldFirst,field)
      end if

      if (trim(field%field) .eq. trim(geopotentialHeightName)) then
         geopotentialHeightCount = geopotentialHeightCount + 1
         call appendFieldWPSList(geopotentialHeightList,field)
      else if (trim(field%field) .eq. trim(temperatureName)) then
         temperatureCount = temperatureCount + 1
         call appendFieldWPSList(temperatureList,field)
      else if (trim(field%field) .eq. trim(relativeHumidityName)) then
         relativeHumidityCount = relativeHumidityCount + 1
         call appendFieldWPSList(relativeHumidityList,field)
      else if (trim(field%field) .eq. trim(uName)) then
         uCount = uCount + 1
         call appendFieldWPSList(uList,field)
      else if (trim(field%field) .eq. trim(vName)) then
         vCount = vCount + 1
         call appendFieldWPSList(vList,field)
      end if

      call destroyFieldWPS(field)
   end do

   call destroyFileWPS(inputFile)
   call close_file_unit(fileUnit1)

   ! Check counts
   if (geopotentialHeightCount .ne. temperatureCount .or. &
       geopotentialHeightCount .ne. relativeHumidityCount .or. &
       geopotentialHeightCount .ne. uCount .or. &
       geopotentialHeightCount .ne. vCount) then
      print*,'ERROR, mismatched counts of geopotential height, temperature,'
      print*,'       relative humidity, u, and v slabs!'
      print*,'    Geopotential height slabs:  ',geopotentialHeightCount
      print*,'    Temperature slabs:  ',temperatureCount
      print*,'    Relative humidity slabs:  ',relativeHumidityCount
      print*,'    U slabs:  ',uCount
      print*,'    V slabs:  ',vCount
      stop 1
   end if
   if (geopotentialHeightCount .eq. 0) then
      print*,'ERROR, no geopotential height data found in WPS file!'
      stop 1
   end if
   numberOfSlabs=geopotentialHeightCount

   ! Collect pressures.
   allocate(pressureXlvls(numberOfSlabs))
   call getXlvlsFieldWPSList(geopotentialHeightList,numberOfSlabs, &
        pressureXlvls)
   call bubbleSort(numberOfSlabs,pressureXlvls)

   ! Now collect the heights, temperatures, relative humidity, u and v into
   ! arrays sorted by xlvl.
   call copyListToArray(fieldFirst%nx,fieldFirst%ny,numberOfSlabs,&
        geopotentialHeightList,trim(geopotentialHeightName), &
        pressureXlvls,heights)
   call copyListToArray(fieldFirst%nx,fieldFirst%ny,numberOfSlabs,&
        temperatureList,trim(temperatureName), &
        pressureXlvls,temperatures)
   call copyListToArray(fieldFirst%nx,fieldFirst%ny,numberOfSlabs,&
        relativeHumidityList,trim(relativeHumidityName), &
        pressureXlvls,relativeHumidities)
   call copyListToArray(fieldFirst%nx,fieldFirst%ny,numberOfSlabs,&
        uList,trim(uName), &
        pressureXlvls,uWinds)
   call copyListToArray(fieldFirst%nx,fieldFirst%ny,numberOfSlabs,&
        vList,trim(vName), &
        pressureXlvls,vWinds)

   ! Now create 3D missing mask
   allocate(missing3d(fieldFirst%nx,fieldFirst%ny,numberOfSlabs))
   allocate(missing1d(numberOfSlabs))
   do k = 1, numberOfSlabs
      missing1d(k) = .false.
      do j = 1, fieldFirst%ny
         do i = 1, fieldFirst%nx
            if (temperatures(i,j,k)       > MISSING .or. &
                relativeHumidities(i,j,k) > MISSING .or. &
                heights(i,j,k)            > MISSING .or. &
                uWinds(i,j,k)             > MISSING .or. &
                vWinds(i,j,k)             > MISSING) then
               missing3d(i,j,k) = .true.
               missing1d(k) = .true.
            else
               missing3d(i,j,k) = .false.
            end if
         end do
      end do
   end do

   ! Now loop through each column, going from model top to model bottom, find
   ! missing values, and extrapolate.
   do j = 1, fieldFirst%ny
      do i = 1, fieldFirst%nx
         zAbove=heights(i,j,1)
         tAbove=temperatures(i,j,1)
         pAbove=pressureXlvls(1)
         RHAbove=relativeHumidities(i,j,1)
         qsatAbove=calcSatSpecificHumidity(tAbove,pAbove)
         qsatAbove=calcMixingRatio(qsatAbove)
         qvAbove=0.01*RHAbove*qsatAbove
         tvAbove=tAbove*(1. + 0.61*qvAbove)

         do k = 2,numberOfSlabs

            if (missing3d(i,j,k)) then

               pBelow=pressureXlvls(k)

               if (temperatures(i,j,k) > MISSING .and. &
                   heights(i,j,k)      > MISSING) then
                  
                  ! General case if both geopotential height and temperature
                  ! are missing.  Extrapolate temperature assuming a 
                  ! constant virtual temperature lapse rate (GAMMA) using 
                  ! an expression derived from the hydrostatic equation.  
                  ! We're cheating a little bit here by using the mixing
                  ! ratio of the upper pressure level to convert to
                  ! temperature, but that is consistent with NCEP's
                  ! procedure for NAM.
                  lnPRatio=log(pBelow/pAbove)
                  tvBelow=tvAbove*exp(lnPRatio*R_OVR_G*GAMMA)
                  temperatures(i,j,k)=tvBelow/(1. + 0.61*qvAbove)
                  
                  ! Use the hypsometric equation to find the geopotential 
                  ! height.
                  meanTv = (tvBelow + tvAbove) * 0.5
                  heights(i,j,k) = zAbove - (R_OVR_G*meanTv*lnPRatio)

                  ! Shuell method from NCEP, used with GFS
                  lnPRatio=log(pBelow/pAbove)
                  if (zAbove > ZSHUL) then
                     tvBelow = tvAbove + (GAMMA*zAbove)
                     if (tvBelow > TVSHUL) then
                        if (tvAbove > TVSHUL) then
                           tvBelow = TVSHUL - 5.e-3*(tvAbove-TVSHUL)**2
                        else
                           tvBelow = TVSHUL
                        end if
                     end if
                     gamma_s=(tvAbove-tvBelow)/zAbove
                  else
                     gamma_s=0.
                  end if
                  part=R_OVR_G*lnPRatio
                  heights(i,j,k) = zAbove - tvAbove*part/(1.+0.5*gamma_s*part)
                  temperatures(i,j,k) = tAbove - GAMMA*(heights(i,j,k)-zAbove)

               else if (temperatures(i,j,k) > MISSING) then
                  ! Situation where heights are available but temperature
                  ! is missing.

                  ! Using hydrostatic approximation and constant lapse rate.
                  lnPRatio=log(pBelow/pAbove)
                  tvBelow=tvAbove*exp(lnPRatio*R_OVR_G*GAMMA)
                  temperatures(i,j,k)=tvBelow/(1. + 0.61*qvAbove)
                  
                  ! From NCEP Shuell method
!                  temperatures(i,j,k) = tAbove - GAMMA*(heights(i,j,k)-zAbove)
                  
               end if
               
               if (relativeHumidities(i,j,k) > MISSING) then
                  relativeHumidities(i,j,k) = relativeHumidities(i,j,k-1)
               end if
               if (uWinds(i,j,k) > MISSING) then
                  uWinds(i,j,k) = uWinds(i,j,k-1)
               end if
               if (vWinds(i,j,k) > MISSING) then
                  vWinds(i,j,k) = vWinds(i,j,k-1)
               end if
               
            else ! No missing data

               zAbove=heights(i,j,k)
               tAbove=temperatures(i,j,k)
               pAbove=pressureXlvls(k)
               RHAbove=relativeHumidities(i,j,k)
               qsatAbove=calcSatSpecificHumidity(tAbove,pAbove)
               qsatAbove=calcMixingRatio(qsatAbove)
               qvAbove=0.01*RHAbove*qsatAbove
               tvAbove=tAbove*(1. + 0.61*qvAbove)

            end if
               
         end do
      end do
   end do

   ! Smooth the fields where data was extrapolated
   do k = 1, numberOfSlabs
      if (missing1d(k)) then
         call smooth9p(fieldFirst%nx,fieldFirst%ny, &
              heights(1,1,k))
         call smooth9p(fieldFirst%nx,fieldFirst%ny, &
              temperatures(1,1,k))
         call smooth9p(fieldFirst%nx,fieldFirst%ny, &
              relativeHumidities(1,1,k))
         call smooth9p(fieldFirst%nx,fieldFirst%ny, &
              uWinds(1,1,k))
         call smooth9p(fieldFirst%nx,fieldFirst%ny, &
              vWinds(1,1,k))         
      end if
   end do

   ! Open new WPS file
   fileUnit2=select_file_unit()
   outputFile=createFileWPS(fileUnit2,trim(directory), &
        'ISOBARIC',year,month,day,hour,minute,second)
   
   ! Write the updated slabs to new file
   call writeIsobaricSlabs(fieldFirst%nx,fieldFirst%ny,numberOfSlabs, &
        geopotentialHeightList,trim(geopotentialHeightName), &
        pressureXlvls,missing1d,heights,outputFile)
   call writeIsobaricSlabs(fieldFirst%nx,fieldFirst%ny,numberOfSlabs, &
        temperatureList,trim(temperatureName), &
        pressureXlvls,missing1d,temperatures,outputFile)
   call writeIsobaricSlabs(fieldFirst%nx,fieldFirst%ny,numberOfSlabs, &
        relativeHumidityList,trim(relativeHumidityName), &
        pressureXlvls,missing1d,relativeHumidities,outputFile)
   call writeIsobaricSlabs(fieldFirst%nx,fieldFirst%ny,numberOfSlabs, &
        uList,trim(uName), &
        pressureXlvls,missing1d,uWinds,outputFile)
   call writeIsobaricSlabs(fieldFirst%nx,fieldFirst%ny,numberOfSlabs, &
        vList,trim(vName), &
        pressureXlvls,missing1d,vWinds,outputFile)

   ! Close output file
   call destroyFileWPS(outputFile)
   call close_file_unit(fileUnit2)

   ! Clean up
   call destroyFieldWPS(fieldFirst)
   call destroyFieldWPSList(geopotentialHeightList)
   call destroyFieldWPSList(temperatureList)
   call destroyFieldWPSList(relativeHumidityList)
   call destroyFieldWPSList(uList)
   call destroyFieldWPSList(vList)
   deallocate(missing3d)
   deallocate(missing1d)
   deallocate(heights)
   deallocate(temperatures)
   deallocate(relativeHumidities)
   deallocate(uWinds)
   deallocate(vWinds)
   deallocate(pressureXlvls)

contains
   
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
      real,allocatable,intent(inout) :: array3d(:,:,:)

      ! Local variables
      type(FieldWPS) :: field
      integer :: i,j,k

      allocate(array3d(nx,ny,nz))

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
   ! ROUTINE:  copyListToArray
   !
   ! DESCRIPTION:  Internal subroutine to employ 9-point smoother on 2-D
   ! array.
   !
   !---------------------------------------------------------------------------

   subroutine smooth9p(nx,ny,array2d)

      ! Arguments
      integer,intent(in) :: nx
      integer,intent(in) :: ny
      real,intent(inout) :: array2d(nx,ny)

      ! Local variables
      real :: tem(nx,ny)
      integer :: i,j
      
      ! Inverse sum of weights for 9-point filter
      real,parameter :: NORM16 = 1.0/16.0
      real,parameter :: NORM12 = 1.0/12.0
      real,parameter :: NORM9 = 1.0/9.0

      ! General-case 9-point smoother
      do j = 2,ny-1
         do i = 2,nx-1
            tem(i,j) =    &
                 array2d(i-1,j+1) + 2.*array2d(i ,j+1) +    array2d(i+1,j+1) &
            + 2.*array2d(i-1,j  ) + 4.*array2d(i ,j  ) + 2.*array2d(i+1,j  ) &
            +    array2d(i-1,j-1) + 2.*array2d(i ,j-1) +    array2d(i+1,j-1)
          
            tem(i,j) = NORM16*tem(i,j)
         end do
      end do

      ! Along west edge
      do j = 2,ny-1
         do i = 1,1
            tem(i,j) =    &
                 2.*array2d(i ,j+1) +    array2d(i+1,j+1) &
               + 4.*array2d(i ,j  ) + 2.*array2d(i+1,j  ) &
               + 2.*array2d(i ,j-1) +    array2d(i+1,j-1)

            tem(i,j) = NORM12*tem(i,j)
         end do
      end do

      ! Along east edge
      do j = 2,ny-1
         do i = nx,nx
            tem(i,j) =   &
                 array2d(i-1,j+1) + 2.*array2d(i ,j+1)                   &
            + 2.*array2d(i-1,j  ) + 4.*array2d(i ,j  )                   &
            +    array2d(i-1,j-1) + 2.*array2d(i ,j-1)
            tem(i,j) = NORM12*tem(i,j)
         end do
      end do

     ! Along south edge
      do j = 1,1
         do i = 2,nx-1      
            tem(i,j) =  &
                 array2d(i-1,j+1) + 2.*array2d(i ,j+1) +    array2d(i+1,j+1) &
          +   2.*array2d(i-1,j  ) + 4.*array2d(i ,j  ) + 2.*array2d(i+1,j  ) 
            tem(i,j) = NORM12*tem(i,j)
         end do
      end do

      ! Along north edge
      do j = ny,ny
         do i = 2,nx-1
            tem(i,j) =  &
               2.*array2d(i-1,j  ) + 4.*array2d(i ,j  ) + 2.*array2d(i+1,j  ) &
             +    array2d(i-1,j-1) + 2.*array2d(i ,j-1) +    array2d(i+1,j-1)
            tem(i,j) = NORM12*tem(i,j)
         end do
      end do

      ! Southwest corner
      do j = 1,1
         do i = 1,1
            tem(i,j) =  &
                 2.*array2d(i ,j+1) +    array2d(i+1,j+1) &
               + 4.*array2d(i ,j  ) + 2.*array2d(i+1,j  ) 
            tem(i,j) = NORM9*tem(i,j)
         end do
      end do

      ! Southeast corner
      do j = 1,1
         do i = nx,nx
            tem(i,j) =   &
                 array2d(i-1,j+1) + 2.*array2d(i ,j+1)                   &
          +   2.*array2d(i-1,j  ) + 4.*array2d(i ,j  )                   
            tem(i,j) = NORM9*tem(i,j)
         end do
      end do

      ! Northwest corner
      do j = ny,ny
         do i = 1,1
            tem(i,j) =   &
                 4.*array2d(i ,j  ) + 2.*array2d(i+1,j  ) &
               + 2.*array2d(i ,j-1) +    array2d(i+1,j-1)
            tem(i,j) = NORM9*tem(i,j)
         end do
      end do

      ! Northeast corner
      do j = ny,ny
         do i = nx,nx
            tem(i,j) =   &
                 2.*array2d(i-1,j  ) + 4.*array2d(i ,j  )                   &
               +    array2d(i-1,j-1) + 2.*array2d(i ,j-1)
            tem(i,j) = NORM9*tem(i,j)
          end do
      end do

      ! Overwrite original slab.
      do j = 1,ny
         do i = 1,nx
            array2d(i,j) = tem(i,j)
         end do
      end do

      return
   end subroutine smooth9p

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  writeIsobaricSlabs
   !
   ! DESCRIPTION:  Internal subroutine to write slabs with updated extrapolated
   ! data to new WPS file.
   !
   !---------------------------------------------------------------------------

   subroutine writeIsobaricSlabs(nx,ny,numberOfSlabs,fieldList,fieldName, &
        pressureXlvls,missing1d,array,outputFile)

      ! Arguments
      integer,intent(in) :: nx,ny,numberOfSlabs
      type(FieldWPSList),intent(in) :: fieldList
      character(len=*),intent(in) :: fieldName
      real,intent(in) :: pressureXlvls(numberOfSlabs)
      logical,intent(in) :: missing1d(numberOfSlabs)
      real,intent(in) :: array(nx,ny,numberOfSlabs)
      type(FileWPS),intent(in) :: outputFile

      ! Local variables
      type(FieldWPS) :: field
      integer :: i,j,k

      do k = 1, numberOfSlabs
         field=findFieldWPSList(fieldList,trim(fieldName),pressureXlvls(k))
         if (missing1d(k)) then
            do j = 1, ny
               do i = 1, nx
                  field%slab(i,j) = array(i,j,k)
               end do
            end do
         end if
         call writeFileWPS(outputFile,field)
         call destroyFieldWPS(field)
      end do      

   end subroutine writeIsobaricSlabs

end program extrapIsobaric
