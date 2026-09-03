!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  DomainData_mod
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Provides routines to capture WRF grid information from namelist.wps file and
! from netCDF output files from geogrid.exe.
!
! REVISION HISTORY:
! 9 May 2011 - Initial version (processed metgrid.exe files)
! 16 Jul 2012 - New version for geogrid.exe files.
! 19 Aug 2014 - Removed checks for reference lat/lon and x/y point. The
!               reference lat/lon are not sent back to LIS and LDT, and
!               insisting that the reference lat/lon are the same as the
!               center lat/lon unnecessarily reduces the utility of this
!               program.
! 20 Oct 2014 - Added multiple nests to output. Added code to convert
!               string to lower case characters (more robust handling of
!               map projection setting.
! 26 Aug 2016 - Added map projection line to domain_data output file.
!------------------------------------------------------------------------------

module DomainData_mod

   ! Change default behavior
   implicit none
   private

   ! Public parameters
   integer,parameter :: DATE_LENGTH=19
   integer,parameter :: MAX_DOM_LIMIT=10
   integer,parameter :: MAP_PROJ_LENGTH=8
   integer,parameter :: GEOGRID_FILENAME_LENGTH=13

   real,parameter :: MISSING = -9999

   character(len=MAP_PROJ_LENGTH),parameter :: LAMBERT="lambert"
   character(len=MAP_PROJ_LENGTH),parameter :: POLAR="polar"
   character(len=MAP_PROJ_LENGTH),parameter :: MERCATOR="mercator"
   character(len=DATE_LENGTH),parameter :: &
        DEFAULT_DATE = 'YYYY-MM-DD_HH:mm:ss'

   ! Public data type
   public :: DomainData
   type DomainData
      
      ! Encapsulate attributes
      private
      
      ! Maximum number of WRF/WPS domains
      integer :: maxDomains
      
      ! Number of west-east grid points for each WRF domain.
      integer,allocatable :: endWestEastGridPoints(:) 
                                     
      ! Number of south-north grid points for each WRF domain.
      integer,allocatable :: endSouthNorthGridPoints(:) 
         
      ! Grid resolutions in west-east direction (m)
      real, allocatable :: dx(:)

      ! Grid resolutions in south-north direction (m)
      real, allocatable :: dy(:)

      ! Map projection 
      character(len=MAP_PROJ_LENGTH) :: mapProjection

!      ! Reference (center) latitude (deg N) of map projection
!      real :: refLatitude

!      ! Reference (center) longitude (deg E) of map projection
!      real :: refLongitude

      ! First true latitude (Lambert projection), or only true latitude
      ! (Mercator and Polar Stereographic projections)
      real :: trueLatitude1

      ! Second true latitude (Lambert projection only)
      real :: trueLatitude2

      ! Standard longitude (Lambert and Polar Stereographic projections)
      real :: standLongitude

      ! Corner latitude (deg N) of each grid
      real, allocatable :: cornerLatitudes(:)

      ! Corner longitude (deg E) of each grid
      real, allocatable :: cornerLongitudes(:)

   end type DomainData

   ! Public methods
   public :: createDomainData
   public :: destroyDomainData
   public :: writeDomainData

contains

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! createDomainData
   !
   ! DESCRIPTION:
   ! Public "constructor method" for DomainData data type.
   !---------------------------------------------------------------------------

   function createDomainData() result(this)
      
      ! Return variable
      type(DomainData) :: this

      ! Local variables
      character(len=GEOGRID_FILENAME_LENGTH),allocatable :: geogridFiles(:)
      integer :: fileUnit
      integer :: i

      ! Initialize data structure
      call destroyDomainData(this)

      ! Open namelist.wps and return file unit number
      fileUnit = openNamelistWps()

      ! Process share namelist
      call processShareNamelistDomainData(this,fileUnit)

      ! Process geogrid namelist
      call processGeogridNamelistDomainData(this,fileUnit)

      ! Close namelist.wps
      call closeNamelistWps(fileUnit)

      ! Get list of geogrid files to process
      allocate(geogridFiles(this%maxDomains))
      call selectGeogridFilesDomainData(this,geogridFiles)

      ! Process each geogrid file
      allocate(this%cornerLatitudes(this%maxDomains))
      allocate(this%cornerLongitudes(this%maxDomains))
      do i = 1,this%maxDomains
         call readGeogridFileDomainData(this,i,geogridFiles(i))
      end do

      ! Clean up
      deallocate(geogridFiles)

   end function createDomainData

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! destroyDomainData
   !
   ! DESCRIPTION:
   ! Public "destructor method" for DomainData data type.
   !---------------------------------------------------------------------------

   subroutine destroyDomainData(this)

      ! Arguments
      type(DomainData),intent(inout) :: this

      ! Clean up data structure
      this%maxDomains=0
      if (allocated(this%endWestEastGridPoints)) then
         deallocate(this%endWestEastGridPoints)
      end if
      if (allocated(this%endSouthNorthGridPoints)) then
         deallocate(this%endSouthNorthGridPoints)
      end if
      if (allocated(this%dx)) deallocate(this%dx)
      if (allocated(this%dy)) deallocate(this%dy)
      if (allocated(this%cornerLatitudes)) then
         deallocate(this%cornerLatitudes)
      end if
      if (allocated(this%cornerLongitudes)) then
         deallocate(this%cornerLongitudes)
      end if
      this%mapProjection="MISSING"
!      this%refLatitude = MISSING
!      this%refLongitude = MISSING
      this%trueLatitude1 = MISSING
      this%trueLatitude2 = MISSING
      this%standLongitude = MISSING

   end subroutine destroyDomainData

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! openNamelistWps
   !
   ! DESCRIPTION:
   ! Internal subroutine for opening namelist.wps file.
   !---------------------------------------------------------------------------

   function openNamelistWps() result (fileUnit)

      ! Import module
      use FileUnit_mod

      ! Return variable
      integer :: fileUnit

      ! Local variables
      integer :: ioStatus

      ! Determine valid file unit for namelist file
      fileUnit = select_file_unit()

      ! Open the namelist.wps file
      open(unit=fileUnit,file="namelist.wps",status="old", &
           iostat=ioStatus, delim='apostrophe')
      if (ioStatus /= 0) then
         write(0,*)'ERROR, cannot open file namelist.wps!'
         stop
      end if

   end function openNamelistWps

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! processShareNamelistDomainData
   !
   ! DESCRIPTION:
   ! Internal subroutine for reading share namelist from namelist.wps file and
   ! saving grid information.
   !---------------------------------------------------------------------------

   subroutine processShareNamelistDomainData(this,fileUnit)
      
      ! Arguments
      type(DomainData),intent(inout) :: this
      integer, intent(in) :: fileUnit

      ! WPS share namelist
      character(len=3) :: wrf_core
      integer :: max_dom
      integer :: start_year(MAX_DOM_LIMIT)
      integer :: start_month(MAX_DOM_LIMIT)
      integer :: start_day(MAX_DOM_LIMIT)
      integer :: start_hour(MAX_DOM_LIMIT)
      integer :: end_year(MAX_DOM_LIMIT)
      integer :: end_month(MAX_DOM_LIMIT)
      integer :: end_day(MAX_DOM_LIMIT)
      integer :: end_hour(MAX_DOM_LIMIT)
      character(len=DATE_LENGTH) :: start_date(MAX_DOM_LIMIT) 
      character(len=DATE_LENGTH) :: end_date(MAX_DOM_LIMIT) 
      integer :: interval_seconds
      integer :: io_form_geogrid
      character(len=132) :: opt_output_from_geogrid_path
      integer :: debug_level
      namelist /share/ wrf_core,max_dom,start_year,start_month,start_day, &
           start_hour,end_year,end_month,end_day,end_hour,start_date, &
           end_date,interval_seconds, io_form_geogrid, &
           opt_output_from_geogrid_path,debug_level

      ! Read share namelist
      read(unit=fileUnit,nml=share)
      
      ! Copy number of domains to structure.
      if (max_dom < 1) then
         write(0,*) &
              'ERROR, max_dom must be positive!  Value read in is ',max_dom
         stop
      endif
      if (max_dom > MAX_DOM_LIMIT) then
         write(0,*) 'ERROR, max_dom exceeds MAX_DOM_LIMIT!'
         write(0,*) 'max_dom = ',max_dom
         write(0,*) 'MAX_DOM_LIMIT = ',MAX_DOM_LIMIT
         write(0,*) &
              'Change MAX_DOM_LIMIT in DomainData_mod.f90 and recompile!'
         stop
      endif
      this%maxDomains = max_dom

   end subroutine processShareNamelistDomainData

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! processGeogridNamelistDomainData
   !
   ! DESCRIPTION:
   ! Internal subroutine for reading geogrid namelist from namelist.wps file 
   ! and saving grid information.
   !---------------------------------------------------------------------------

   subroutine processGeogridNamelistDomainData(this,fileUnit)
      
      ! Arguments
      type(DomainData),intent(inout) :: this
      integer,intent(in) :: fileUnit

      ! WPS geogrid namelist
      integer :: parent_id(MAX_DOM_LIMIT)
      integer :: parent_grid_ratio(MAX_DOM_LIMIT)
      integer :: i_parent_start(MAX_DOM_LIMIT)
      integer :: j_parent_start(MAX_DOM_LIMIT)
      integer :: s_we(MAX_DOM_LIMIT)
      integer :: e_we(MAX_DOM_LIMIT)
      integer :: s_sn(MAX_DOM_LIMIT)
      integer :: e_sn(MAX_DOM_LIMIT)
      character(len=4) :: geog_data_res(MAX_DOM_LIMIT)
      real :: dx
      real :: dy
      character(len=MAP_PROJ_LENGTH) :: map_proj
      real :: ref_lat
      real :: ref_lon
      real :: ref_x
      real :: ref_y
      real :: truelat1
      real :: truelat2
      real :: stand_lon
      real :: pole_lat
      real :: pole_lon
      character(len=132) :: geog_data_path
      character(len=132) :: opt_geogrid_tbl_path
      namelist /geogrid/ parent_id,parent_grid_ratio,i_parent_start, &
           j_parent_start, s_we, e_we, s_sn, e_sn, geog_data_res, dx, dy, &
           map_proj, ref_lat, ref_lon, ref_x, ref_y, truelat1, truelat2, &
           stand_lon, pole_lat, pole_lon, geog_data_path, opt_geogrid_tbl_path

      ! Other local variables
      integer :: max_dom
      integer :: i

      ! Initialize ref_x and ref_y for later checks
      ref_x = MISSING
      ref_y = MISSING

      ! Allocate portions of data structure
      max_dom = this%maxDomains
      allocate(this%endWestEastGridPoints(max_dom), &
               this%endSouthNorthGridPoints(max_dom), &
               this%dx(max_dom), &
               this%dy(max_dom))     

      ! Process geogrid namelist
      read(unit=fileUnit,nml=geogrid)

      do i = 1,max_dom
         if (e_we(i) < 1) then
            write(0,*) 'ERROR, negative value of e_we read in!'
            write(0,*) 'e_we = ',e_we
            stop
         end if
         this%endWestEastGridPoints(i) = e_we(i)

         if (e_sn(i) < 1) then
            write(0,*) 'ERROR, negative value of e_sn read in!'
            write(0,*) 'e_sn = ',e_sn
            stop
         end if
         this%endSouthNorthGridPoints(i) = e_sn(i)
      end do

      if (.not. dx > 0) then
         write(0,*) 'ERROR, dx is supposed to be positive!'
         write(0,*) 'dx = ',dx
         stop
      end if
      this%dx(1) = dx
      
      if (.not. dy > 0) then
         write(0,*) 'ERROR, dy is supposed to be positive!'
         write(0,*) 'dy = ',dy
         stop
      end if
      this%dy(1) = dy

      do i = 2,max_dom
         if (parent_grid_ratio(i) <= 2) then
            write(0,*)'ERROR, parent_grid_ratio must be > 2!'
            write(0,*)'parent_grid_ratio = ',parent_grid_ratio
            stop
         end if
         this%dx(i) = this%dx(i-1) / parent_grid_ratio(i)
         this%dy(i) = this%dy(i-1) / parent_grid_ratio(i)
      end do
      
      ! Make sure ref_x and ref_y are set to the center of the first domain;
      ! otherwise, ref_lat and ref_lon will be misinterpreted.
!       if (ref_x /= MISSING) then
!          if ( abs(ref_x - (e_we(1)/2.)) > tiny(ref_x)) then
!             write(0,*) 'ERROR, ref_x is not set at center of the domain!'
!             write(0,*) 'ref_x = ',ref_x
!             write(0,*) '(e_we(1)/2.) = ',(e_we(1)/2.)
!             stop
!          end if
!       end if
!       if (ref_y /= MISSING) then
!          if ( abs(ref_y - (e_sn(1)/2.)) > tiny(ref_y)) then
!             write(0,*) 'ERROR, ref_y is not set at center of the domain!'
!             write(0,*) 'ref_y = ',ref_y
!             write(0,*) '(e_sn(1)/2.) = ',(e_sn(1)/2.)
!             stop
!          end if
!       end if

      this%mapProjection = map_proj
      ! Make sure all characters are lower case.
      call lowercase(this%mapProjection) 
!      this%refLatitude = ref_lat
!      this%refLongitude = ref_lon
      this%trueLatitude1 = truelat1
      this%trueLatitude2 = truelat2
      this%standLongitude = stand_lon
      
   end subroutine processGeogridNamelistDomainData

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! closeNamelistWps
   !
   ! DESCRIPTION:
   ! Internal subroutine for closing namelist.wps file.
   !---------------------------------------------------------------------------

   subroutine closeNamelistWps(fileUnit)

      ! Import modules
      use FileUnit_mod

      ! Arguments
      integer,intent(in) :: fileUnit

      ! Close file unit.
      call close_file_unit(fileUnit)

   end subroutine closeNamelistWps

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! selectGeogridFilesDomainData
   !
   ! DESCRIPTION:
   ! Internal subroutine for assembling geogrid filenames to read.
   !---------------------------------------------------------------------------

   subroutine selectGeogridFilesDomainData(this,geogridFiles)

      ! Arguments
      type(DomainData),intent(in) :: this
      character(len=GEOGRID_FILENAME_LENGTH) :: geogridFiles(this%maxDomains)

      ! Local variables
      character(len=GEOGRID_FILENAME_LENGTH) :: tmpString
      integer :: i

      ! Create list of file names
      do i = 1,this%maxDomains
         write(tmpString,'(A,I2.2,A)') &
              'geo_em.d',i,'.nc'
         geogridFiles(i) = tmpString
      end do

   end subroutine selectGeogridFilesDomainData

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! readGeogridFileDomainData
   !
   ! DESCRIPTION:
   ! Internal subroutine for reading geogrid netCDF file for corner lat/lon.
   !---------------------------------------------------------------------------

   subroutine readGeogridFileDomainData(this,domainID,geogridFile)

      ! Arguments
      type(DomainData),intent(inout) :: this
      integer,intent(in) :: domainID
      character(len=GEOGRID_FILENAME_LENGTH),intent(in) :: geogridFile

      ! Local variables
      integer :: ncid
      real :: cornerLatitudes(16)
      real :: cornerLongitudes(16)

      ! Include netCDF header file
      include 'netcdf.inc'

      ! Open the netCDF file
      call checkNetcdfStatus(nf_open(geogridFile,NF_NOWRITE,ncid))
      
      ! Fetch corner_lats
      call checkNetcdfStatus(nf_get_att_real(ncid,NF_GLOBAL,"corner_lats", &
           cornerLatitudes))

      ! Fetch corner lons
      call checkNetcdfStatus(nf_get_att_real(ncid,NF_GLOBAL,"corner_lons", &
           cornerLongitudes))

      ! Close the netCDF file
      call checkNetcdfStatus(nf_close(ncid))

      ! Copy the first corner lat/lon
      this%cornerLatitudes(domainID) = cornerLatitudes(1)
      this%cornerLongitudes(domainID) = cornerLongitudes(1)

   end subroutine readGeogridFileDomainData

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! checkNetcdfStatus
   !
   ! DESCRIPTION:
   ! Internal subroutine for reading netcdf status results.
   !---------------------------------------------------------------------------

   subroutine checkNetcdfStatus(status)
      
      ! Arguments
      integer,intent(in) :: status

      ! Include netCDF header file
      include 'netcdf.inc'

      ! Handle bad status
      if (status /= NF_NOERR) then
         write(0,*) 'ERROR returned from netCDF!'
         write(0,*) trim(nf_strerror(status))
         stop
      end if

   end subroutine checkNetcdfStatus

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! writeDomainData
   !
   ! DESCRIPTION:
   ! Public method for writing WRF domain.
   !---------------------------------------------------------------------------

   subroutine writeDomainData(this)

      ! Import modules
      use FileUnit_mod

      ! Arguments
      type(DomainData),intent(in) :: this

      ! Local variables
      character(len=11),parameter :: outputFile = 'domain_data'
      integer :: fileUnit
      integer :: ioStatus
      character(len=132) :: tmpString, tmpStringLong
      integer :: i

      ! Determine valid file unit for output file
      fileUnit = select_file_unit()

      ! Open the output file
      open(unit=fileUnit,file=outputFile,status="unknown", &
           iostat=ioStatus)
      if (ioStatus /= 0) then
         write(0,*)'ERROR, cannot open file ',outputFile
         stop
      end if

      ! Write out the number of nests
      ! First, for ldt.config
      tmpStringLong = 'LIS number of nests: '
      write(tmpString,'(I3)') this%maxDomains
      tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
      write(fileUnit,*) trim(tmpStringLong)
      ! Second, for list.config
      tmpStringLong = 'Number of nests: '
      write(tmpString,'(I3)') this%maxDomains
      tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
      write(fileUnit,*) trim(tmpStringLong)

      ! Write out the domain data
      if ( trim(this%mapProjection) == LAMBERT .or. &
           trim(this%mapProjection) == POLAR .or. &
           trim(this%mapProjection) == MERCATOR ) then

         tmpStringLong = 'Map projection of the LIS domain:'
         if (trim(this%mapProjection) == LAMBERT) then
            write(tmpString,'(A)') 'lambert'
         else if (trim(this%mapProjection) == POLAR) then
            write(tmpString,'(A)') 'polar'
         else if (trim(this%mapProjection) == MERCATOR) then
            write(tmpString,'(A)') 'mercator'
         end if
         tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
         write(fileUnit,*) trim(tmpStringLong)

         tmpStringLong = 'Run domain lower left lat:  '
         do i = 1,this%maxDomains
            write(tmpString,'(F10.5)') this%cornerLatitudes(i)
            tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
         end do
         write(fileUnit,*) trim(tmpStringLong)

         tmpStringLong = 'Run domain lower left lon:  '
         do i = 1,this%maxDomains
            write(tmpString,'(F10.5)') this%cornerLongitudes(i)
            tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
         end do
         write(fileUnit,*) trim(tmpStringLong)

         if ( trim(this%mapProjection) == LAMBERT ) then

            tmpStringLong = 'Run domain true lat1:  '
            do i = 1,this%maxDomains
               write(tmpString,'(F10.5)') this%trueLatitude1
               tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
            end do
            write(fileUnit,*) trim(tmpStringLong)

            tmpStringLong = 'Run domain true lat2:  '
            do i = 1,this%maxDomains
               write(tmpString,'(F10.5)') this%trueLatitude2
               tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
            end do
            write(fileUnit,*) trim(tmpStringLong)

         else

            tmpStringLong = 'Run domain true lat:  '
            do i = 1,this%maxDomains
               write(tmpString,'(F10.5)') this%trueLatitude1
               tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
            end do
            write(fileUnit,*) trim(tmpStringLong)

         end if

         tmpStringLong = 'Run domain standard lon:  '
         do i = 1,this%maxDomains
            write(tmpString,'(F10.5)') this%standLongitude
            tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
         end do
         write(fileUnit,*) trim(tmpStringLong)

         if ( trim(this%mapProjection) == POLAR ) then
            tmpStringLong = 'Run domain orientation:  '
            do i = 1,this%maxDomains
               write(tmpString,'(F10.5)') 0.0
               tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
            end do
            write(fileUnit,*) trim(tmpStringLong)
         end if

         tmpStringLong = 'Run domain resolution:  '
         do i = 1,this%maxDomains
            write(tmpString,'(F10.5)') this%dx(i)*0.001
            tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
         end do
         write(fileUnit,*) trim(tmpStringLong)

         tmpStringLong = 'Run domain x-dimension size:  '
         do i = 1,this%maxDomains
            write(tmpString,'(I5)') this%endWestEastGridPoints(i) - 1
            tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
         end do
         write(fileUnit,*) trim(tmpStringLong)

         tmpStringLong = 'Run domain y-dimension size:  '
         do i = 1,this%maxDomains
            write(tmpString,'(I5)') this%endSouthNorthGridPoints(i) - 1
            tmpStringLong = trim(tmpStringLong)//' '//trim(tmpString)
         end do
         write(fileUnit,*) trim(tmpStringLong)

      else
         write(6,*)'ERROR, invalid map projection read in!'
         write(6,*)'Read in ',trim(this%mapProjection)
         write(6,*)'Valid WRF-LIS options are ',trim(LAMBERT), &
              ', ',trim(MERCATOR),', ',trim(POLAR)
         stop 1
      end if
      
      ! Close file
      call close_file_unit(fileUnit)
   end subroutine writeDomainData

   !---------------------------------------------------------------------------
   ! ROUTINE:
   ! lowercase
   !
   ! DESCRIPTION:
   ! Changes characters in string to all lower case.
   !---------------------------------------------------------------------------

   subroutine lowercase(string)
      implicit none
      ! Arguments
      character(len=*),intent(inout) :: string
      ! Internal variables
      integer :: i,imax
      ! Internal constants
      integer,parameter :: UPPER_TO_LOWER = iachar("a") - iachar("A")
      imax = len_trim(string)
      do i = 1,imax
         if(string(i:i) .ge. "A" .and. string(i:i) .le. "Z") then
            string(i:i) = achar(iachar(string(i:i)) + UPPER_TO_LOWER)
         end if
      end do
   end subroutine lowercase

end module DomainData_mod
