#ifdef WRF_LIS
!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  ReconcileWithLis_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Contains code to update WRF land use variables in real.exe with LDT/LIS 7 
! values.
!
!------------------------------------------------------------------------------

module ReconcileWithLis_mod

   ! Import modules   
   use module_date_time                          ! From WRF
   use module_domain                             ! From WRF
   use module_wrf_error                          ! From WRF
   use module_configure                          ! From WRF
   use NetcdfUtils_mod                           ! New code

   ! Change defaults
   implicit none
   private

   ! Include netCDF constants and routines. Keep private.
   ! TODO: Rewrite this module to use netCDF Fortran 90 calls directly,
   ! and remove NetcdfUtils_mod.
   include 'netcdf.inc'

   ! Public methods
   public :: destroyLisData
   public :: createLisData
   public :: setRealData
   public :: checkTsk
   public :: checkTslb
   public :: compareHtWithLIS

   ! Collection of LIS and LDT data
   type LisData
      private
      integer :: soilWater  ! ID for "water" soil category
      integer :: numSoilCategories ! Number of soil categories
      integer :: luWater    ! ID for "water" land use category
      integer :: luSnowIce  ! ID for "snow/ice" land use category
      integer :: luUrban    ! ID for "urban" land use category
      integer :: numLUCategories  ! Number of land use categories
      character(len=256) :: luClassSystem ! Name of land use category system 
                                          ! (e.g., "USGS")
      integer :: is, ie, js, je ! WRF process boundaries, i.e., the WRF memory
                                ! boundaries constrained by the WRF domain
                                ! boundaries.
      integer :: numSoilLayers    ! Number of soil layers
 
      real,allocatable :: Landmask(:,:)
      real,allocatable :: Landcover(:,:)
      real,allocatable :: Soiltype(:,:)
      real,allocatable :: Elevation(:,:)
      real,allocatable :: Shdmax(:,:)
      real,allocatable :: Shdmin(:,:)
      real,allocatable :: Snowcover(:,:)
      real,allocatable :: Snowdepth(:,:)
      real,allocatable :: SWE(:,:)
      real,allocatable :: CanopInt(:,:)
      real,allocatable :: AvgSurfT(:,:)
      real,allocatable :: Albedo_inst(:,:) ! Instantaneous
      real,allocatable :: Albedo(:,:,:)    ! Monthly
      real,allocatable :: Mxsnowalbedo(:,:)
      real,allocatable :: Tempbot(:,:)
      real,allocatable :: Greenness_inst(:,:) ! Instantaneous
      real,allocatable :: Landusef(:,:,:)
      real,allocatable :: Greenness(:,:,:) ! Monthly climo
      real,allocatable :: SoilMoist(:,:,:) ! nx,ny,soillevel
      real,allocatable :: SoilTemp(:,:,:)  ! nx,ny,soillevel
      real,allocatable :: SoilWet(:,:,:)   ! nx,ny,soillevel
      real,allocatable :: RelSMC(:,:,:)    ! nx,ny,soillevel

      real,allocatable :: LAI_inst(:,:)
      real,allocatable :: zs(:) ! Soil layer depth
      real,allocatable :: dzs(:) ! Soil layer thickness
   end type LisData
   public :: LisData

   ! Landuse table options
   integer,parameter :: LIS_LANDUSE_USGS  = 1
   integer,parameter :: LIS_LANDUSE_MODIS = 2
   integer,parameter :: LIS_LANDUSE_UMD   = 3

   ! Landuse water flags
   integer,parameter :: WATER_LANDUSE_USGS  = 16
   integer,parameter :: WATER_LANDUSE_MODIS = 17
   integer,parameter :: WATER_LANDUSE_UMD   = 14
   integer,parameter :: LIS_WATER_LANDUSE   = -9999

   ! Soiltype water flag
   integer,parameter :: WATER_SOILTYPE = 14

   ! Landuse snow/ice flag
   integer,parameter :: SNOWICE_LANDUSE_USGS  = 24
   integer,parameter :: SNOWICE_LANDUSE_MODIS = 15
   integer,parameter :: SNOWICE_LANDUSE_UMD   = -1  ! No such category

   ! Landuse urban
   integer,parameter :: URBAN_LANDUSE_USGS  = 1
   integer,parameter :: URBAN_LANDUSE_MODIS = 13
   integer,parameter :: URBAN_LANDUSE_UMD   = 13

   ! Landuse systems. Names vary between WRF and LDT.
   character(len=256),parameter :: LANDUSE_USGS_IN_WRF = 'USGS'   
   character(len=256),parameter :: LANDUSE_MODIS_IN_WRF = &
        'MODIFIED_IGBP_MODIS_NOAH'
   character(len=256),parameter :: LANDUSE_UMD_IN_WRF = 'UMD'


   character(len=256),parameter :: LANDUSE_USGS_IN_LDT = 'USGS'   
   character(len=256),parameter :: LANDUSE_MODIS_IN_LDT = 'IGBPNCEP'
   character(len=256),parameter :: LANDUSE_UMD_IN_LDT = 'UMD'

   ! Number of land use categories
   integer,parameter :: NUM_CAT_LANDUSE_USGS = 24
   integer,parameter :: NUM_CAT_LANDUSE_MODIS = 20
   integer,parameter :: NUM_CAT_LANDUSE_UMD = 14

   ! Number of soil categories
   integer,parameter :: NUM_CAT_SOIL = 16

   ! Max message length for calling wrf_error_fatal
   integer,parameter :: MAX_MESSAGE_LEN = 256
   integer,parameter :: MAX_LANDCOVER_SCHEME_LEN = 80

contains

   !---------------------------------------------------------------------------
   ! Public destructor
   !---------------------------------------------------------------------------

   subroutine destroyLisData(this)
      implicit none

      ! Arguments
      type(LisData), intent(inout) :: this

      ! Local variables
      integer :: status
      character(len=MAX_MESSAGE_LEN) :: message

      message = 'Deallocation failure'

      if (allocated(this%Landmask)) then
         deallocate(this%Landmask,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Landcover)) then
         deallocate(this%Landcover,stat=status)          
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Soiltype)) then
         deallocate(this%Soiltype,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Elevation)) then
         deallocate(this%Elevation,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Shdmax)) then
         deallocate(this%Shdmax,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Shdmin)) then
         deallocate(this%Shdmin,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%SnowCover)) then
         deallocate(this%SnowCover,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%SnowDepth)) then
         deallocate(this%SnowDepth,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%SWE)) then
         deallocate(this%SWE,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%CanopInt)) then
         deallocate(this%CanopInt,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%AvgSurfT)) then
         deallocate(this%AvgSurfT,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Albedo_inst)) then
         deallocate(this%Albedo_inst,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Albedo)) then
         deallocate(this%Albedo,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Mxsnowalbedo)) then
         deallocate(this%Mxsnowalbedo,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Tempbot)) then
         deallocate(this%Tempbot,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Greenness_inst)) then
         deallocate(this%Greenness_inst,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Greenness)) then
         deallocate(this%Greenness,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%Landusef)) then
         deallocate(this%Landusef,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%SoilMoist)) then
         deallocate(this%SoilMoist,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%SoilTemp)) then
         deallocate(this%SoilTemp,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%SoilWet)) then
         deallocate(this%SoilWet,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%RelSMC)) then
         deallocate(this%RelSMC,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%LAI_inst)) then
         deallocate(this%LAI_inst,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%zs)) then
         deallocate(this%zs,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if
      if (allocated(this%dzs)) then   
         deallocate(this%dzs,stat=status) 
         if (status .ne. 0) call wrf_error_fatal(trim(message))
      end if

      this%is = 0 ; this%ie = 0 ; this%js = 0 ; this%je = 0
      this%numSoilLayers = 0
      
      this%soilWater = 0
      this%numSoilCategories = 0
      this%luWater = 0
      this%luSnowIce = 0
      this%luUrban = 0
      this%numLUCategories = 0
      this%luClassSystem = 'NULL'
   end subroutine destroyLisData

   !---------------------------------------------------------------------------
   ! Public constructor
   !---------------------------------------------------------------------------

   function createLisData(grid,ims,ime,jms,jme,ids,ide,jds,jde) result (this)
      
      ! Change defaults
      implicit none

      ! Arguments
      type(domain),intent(in) :: grid
      integer, intent(in) :: ims, ime, jms, jme ! WRF horizontal memory limits
      integer, intent(in) :: ids,ide,jds,jde ! Max WRF domain limits

      ! Return variable
      type(LisData) :: this

      ! Local variables
      character(len=256) :: lisFilename
      integer :: ncID
      character(len=MAX_MESSAGE_LEN) :: message
      character(len=MAX_LANDCOVER_SCHEME_LEN) :: lis_landcover_scheme
      integer :: status
      integer :: numSoilLayers
      integer :: k
      
      ! Initialize the data structure
      call destroyLisData(this)

      ! Get unstaggered horizontal dimensions
      this%is = max(ims,ids)
      this%js = max(jms,jds)
      this%ie = min(ime,ide-1)
      this%je = min(jme,jde-1)

      ! Fetch namelist settings
      call nl_get_lis_filename(grid%id,lisFilename)

      ! Fetch the LIS and LDT data
      ncID = openNetcdfReadfile(trim(lisFilename))

      ! Get soil layer thickness and depths
      call readNetcdfRealGlobalAttribute(ncID,'soil_layer_thicknesses', &
           this%numSoilLayers,this%dzs)

      allocate(this%zs(this%numSoilLayers),stat=status)
      if (status .ne. 0) then
         write(message,*) 'Error allocating this%zs'
         call wrf_error_fatal(message)
      end if
      this%zs(1) = 0
      this%zs(2) = this%dzs(2)
      do k = 3,this%numSoilLayers
         this%zs(k) = this%zs(k) + this%dzs(k)
      end do

      ! Read 2D data for each WRF memory patch.
      call readLisRealArray2d(ncID,'Landmask',this%is,this%ie, &
           this%js,this%je,this%Landmask)
      call readLisRealArray2d(ncID,'Landcover',this%is,this%ie, &
           this%js,this%je,this%Landcover)
      call readLisRealArray2d(ncID,'Soiltype',this%is,this%ie, &
           this%js,this%je,this%Soiltype)
      call readLisRealArray2d(ncID,'Elevation',this%is,this%ie, &
           this%js,this%je,this%Elevation)
      call readLisRealArray2d(ncID,'Shdmax',this%is,this%ie, &
           this%js,this%je,this%Shdmax)
      call readLisRealArray2d(ncID,'Shdmin',this%is,this%ie, &
           this%js,this%je,this%Shdmin)
      call readLisRealArray2d(ncID,'Snowcover',this%is,this%ie, &
           this%js,this%je,this%Snowcover)
      call readLisRealArray2d(ncID,'Snowdepth',this%is,this%ie, &
           this%js,this%je,this%Snowdepth)
      call readLisRealArray2d(ncID,'SWE',this%is,this%ie, &
           this%js,this%je,this%SWE)
      call readLisRealArray2d(ncID,'CanopInt',this%is,this%ie, &
           this%js,this%je,this%CanopInt)
      call readLisRealArray2d(ncID,'AvgSurfT',this%is,this%ie, &
           this%js,this%je,this%AvgSurfT)
      call readLisRealArray2d(ncID,'Albedo_inst',this%is,this%ie, &
           this%js,this%je,this%Albedo_inst)
      call readLisRealArray2d(ncID,'Mxsnowalbedo',this%is,this%ie, &
           this%js,this%je,this%Mxsnowalbedo)
      call readLisRealArray2d(ncID,'Tempbot',this%is,this%ie, &
           this%js,this%je,this%Tempbot)
      call readLisRealArray2d(ncID,'Greenness_inst',this%is,this%ie, &
           this%js,this%je,this%Greenness_inst)
!      call readLisRealArray2d(ncID,'LAI_inst',this%is,this%ie, &
!      this%js,this%je,this%LAI_inst)

      ! Read 3d data for each WRF memory patch.
      ! These routines only check the horizontal dimensions. We'll check the 
      ! vertical dimension further below.
      call readLisRealArray3d(ncID,'SoilTemp',this%is,this%ie, &
           this%js,this%je,this%SoilTemp)
      call readLisRealArray3d(ncID,'SoilMoist',this%is,this%ie, &
           this%js,this%je,this%SoilMoist)
      call readLisRealArray3d(ncID,'SoilWet',this%is,this%ie, &
           this%js,this%je,this%SoilWet)
      call readLisRealArray3d(ncID,'RelSMC',this%is,this%ie, &
           this%js,this%je,this%RelSMC)

      ! The number of LIS soil layers may differ from METGRID. That's okay, 
      ! we'll use the LIS value.  But we need to ensure the soil layer number
      ! is identical for all LIS values we read in.
      ! TODO: Read in metadata listing the LSM that LIS used and the number 
      ! of soil layers.
      if (this%numSoilLayers .ne. size(this%SoilMoist,3)) then
         write(message,*) 'Soil layer dimension mismatch for ', &
              ' SoilMoist, expected ' &
              ,this%numSoilLayers,' found ',size(this%SoilMoist,3)
         call wrf_error_fatal(message)
      end if
      if (this%numSoilLayers .ne. size(this%SoilWet,3)) then
         write(message,*) 'Soil layer dimension mismatch for ', &
              ' SoilWet, expected ' &
              ,this%numSoilLayers,' found ',size(this%SoilWet,3)
         call wrf_error_fatal(message)
      end if
      if (this%numSoilLayers .ne. size(this%RelSMC,3)) then
         write(message,*) 'Soil layer dimension mismatch for ', &
              ' RelSMC, expected ' &
              ,this%numSoilLayers,' found ',size(this%RelSMC,3)
         call wrf_error_fatal(message)
      end if

      ! Get fractions of each land use (or "vegetation") type
      call readLisRealArray3d(ncID,'Landusef',this%is,this%ie, &
           this%js,this%je,this%Landusef)
      this%numLUCategories = size(this%Landusef,3)

      ! Get time-varying greennness
      call readLisRealArray3d(ncID,'Greenness',this%is,this%ie, &
           this%js,this%je,this%Greenness)
      ! FIXME Account for case with real-time greenness
      if (size(this%Greenness,3) .ne. 12) then
         write(message,*) "Time-varying greenness isn't monthly!"
         call wrf_error_fatal(message)
      end if

      ! Get time-varying albedo
      call readLisRealArray3d(ncID,'Albedo',this%is,this%ie, &
           this%js,this%je,this%Albedo)
      if (size(this%Albedo,3) .ne. 12) then
         write(message,*) "Time-varying albedo isn't monthly!"
         call wrf_error_fatal(message)
      end if

      ! Get remaining land use info
      call readNetcdfTextGlobalAttribute(ncID,'landcover_scheme', &
           MAX_LANDCOVER_SCHEME_LEN,lis_landcover_scheme)
      call setLisLandSoilCodes(this,trim(lis_landcover_scheme))

      call closeNetcdfFile(ncID)

   end function createLisData

   !---------------------------------------------------------------------------
   ! Public subroutine for updating METGRID data with LIS/LDT data.
   !---------------------------------------------------------------------------

   subroutine setRealData(this,grid,its,ite,jts,jte,ide,jde, &
        we_have_tsk,we_have_tavgsfc,internal_time_loop, &
        flag_sst, flag_pinterp, flag_ptheta, &
        flag_icefrac,flag_icepct,config_flags)

      implicit none

      ! Arguments
      type(LisData),intent(in) :: this
      type(domain),intent(inout) :: grid
      integer,intent(in) :: its,ite,jts,jte
      integer,intent(in) :: ide,jde
      logical,intent(in) :: we_have_tsk
      logical,intent(in) :: we_have_tavgsfc
      integer,intent(in) :: internal_time_loop
      integer,intent(in) :: flag_sst
      integer,intent(in) :: flag_pinterp
      integer,intent(in) :: flag_ptheta
      integer,intent(in) :: flag_icefrac
      integer,intent(in) :: flag_icepct
      type(grid_config_rec_type),intent(in) :: config_flags

      ! Local variables
      integer :: i,j

      ! Reallocate the REAL/WRF soil layer arrays
      call setRealSoilDepths(this,grid)

      ! If we have "optional seaice", use that for xice (from REAL code).
      if (flag_icefrac .eq. 1) then
         do j = jts,min(jde-1,jte)
            do i = its, min(ide-1,ite)
               grid%xice(i,j) = grid%icefrac_gc(i,j)
            end do
         end do
      end if

      ! If we have "optional seaice percent" from AFWA, use that for xice
      ! (from REAL code).
      if (flag_icepct .eq. 1) then
         do j = jts,min(jde-1,jte)
            do i = its, min(ide-1,ite)
               grid%xice(i,j) = grid%icepct(i,j)*0.01
            end do
         end do
      end if

      ! Modify points where LIS has land
      call setRealDataLandPoints(this,grid,its,ite,jts,jte,ide,jde, &
           internal_time_loop)

      ! Modify points where LIS has water
      call setRealDataWaterPoints(this,grid,its,ite,jts,jte,ide,jde, &
           we_have_tsk,we_have_tavgsfc,internal_time_loop, &
           flag_sst,flag_pinterp,flag_ptheta, &
           config_flags)

      ! Change the fractional land use in each soil layer. This requires
      ! reallocating the array to change the number of soil layers.
      call setRealLandUseFractions(this,grid,its,ite,jts,jte,ide,jde)

      ! Update the land use metadata
      call setRealLandUseFlags(this,grid)

   end subroutine setRealData

   !---------------------------------------------------------------------------

   subroutine setRealDataLandPoints(this,grid,its,ite,jts,jte,ide,jde, &
        internal_time_loop)

      ! Change defaults
      implicit none

      ! Arguments
      type(LisData),intent(in) :: this
      type(domain),intent(inout) :: grid
      integer,intent(in) :: its,ite,jts,jte
      integer,intent(in) :: ide,jde
      integer,intent(in) :: internal_time_loop

      ! Local variables
      integer :: i,j,k
      character(len=MAX_MESSAGE_LEN) :: message
      character(len=2) :: mon,day15
      integer :: middle(0:13)
      integer :: julyr,julday,target_julyr,target_julday,target_date
      integer :: int_month,month1,month2
      integer :: k1,k2
      real :: gmt
      
      ! Temporally interpolate background albedo and (after first time level)
      ! the veg fraction (e.g., for wrflowinp)
      write(day15,fmt='(i2.2)') 15

      do k = 1,12
         write(mon,fmt='(i2.2)') k
         call get_julgmt( &
              current_date(1:4)//'-'//mon//'-'//day15//'_'//'00:00:00.0000', &
              julyr,julday,gmt)
         middle(k) = julyr*1000 + julday
      end do
      middle( 0) = middle( 1) - 31
      middle(13) = middle(12) + 31

      call get_julgmt(current_date, target_julyr, target_julday, gmt)
      target_date = target_julyr * 1000 + target_julday
      do k = 0,12
         if ( (middle(k  ) .lt. target_date) .and. &
              (middle(k+1) .ge. target_date) ) then
            int_month = k
            k1 = k
            k2 = k+1
            if ( (int_month .eq. 0) .or. (int_month .eq. 12) ) then
               month1 = 12
               month2 = 1
            else
               month1 = int_month
               month2 = month1 + 1
            end if
            exit ! Get out of k loop
         end if
      end do

      ! Update 2D land data at Land points
      do j = jts,min(jde-1,jte)
         do i = its,min(ide-1,ite)

            if (this%Landmask(i,j) .lt. 0.5) cycle ! Skip water points

            grid%sst(i,j) = 0

            grid%landmask(i,j) = 1
            grid%xland(i,j) = 1

            if (internal_time_loop .eq. 1) then
               grid%vegfra(i,j) = this%Greenness_inst(i,j)*100.
            else
               grid%vegfra(i,j) = &
             ( this%Greenness(i,j,month2) * (target_date - middle(k1)) + &
               this%Greenness(i,j,month1) * (middle(k2) - target_date) ) / &
                 (middle(k2) - middle(k1))
               grid%vegfra(i,j) = grid%vegfra(i,j) * 100.
            end if

            grid%shdmax(i,j) = this%Shdmax(i,j)*100.
            grid%shdmin(i,j) = this%Shdmin(i,j)*100.
            grid%canwat(i,j) = this%CanopInt(i,j)
            grid%soilcat(i,j) = this%Soiltype(i,j)
            grid%isltyp(i,j) = this%Soiltype(i,j)
            grid%vegcat(i,j) = this%Landcover(i,j)
            grid%ivgtyp(i,j) = this%Landcover(i,j)
            grid%lu_index(i,j) = this%Landcover(i,j)
            ! FIXME Get soil type fractions from LIS/LDT
            grid%soilctop(i,:,j) = 0
            grid%soilcbot(i,:,j) = 0
            grid%soilctop(i,int(this%Soiltype(i,j)),j) = 1
            grid%soilcbot(i,int(this%Soiltype(i,j)),j) = 1

            grid%snow(i,j) = this%SWE(i,j)
            grid%snowh(i,j) = this%Snowdepth(i,j)
            grid%snowc(i,j) = this%Snowcover(i,j)

            grid%tavgsfc(i,j) = this%AvgSurfT(i,j)
            grid%tsk(i,j) = this%AvgSurfT(i,j)
            grid%ht(i,j) = this%Elevation(i,j)

            ! Handle current and time-varying albedo
            grid%albedo(i,j) = this%Albedo_inst(i,j)
            grid%albbck(i,j) = &
             ( this%Albedo(i,j,month2) * (target_date - middle(k1)) + &
               this%Albedo(i,j,month1) * (middle(k2) - target_date) ) / &
                 (middle(k2) - middle(k1))

            grid%snoalb(i,j) = this%Mxsnowalbedo(i,j)

            ! Just use the LIS deep soil temperature.
            grid%tmn(i,j) = this%Tempbot(i,j)

            ! Remove positive "sea ice fraction" over land
            grid%xice(i,j) = 0.
         end do
      end do

      ! Update 3D soil layer data at Land points
      if (internal_time_loop .eq. 1) then
         do k = 1, this%numSoilLayers
            do j = jts,min(jde-1,jte)
               do i = its,min(ide-1,ite)
                  
                  if (this%Landmask(i,j) .lt. 0.5) cycle ! Skip water points
                  
                  grid%tslb(i,k,j) = this%SoilTemp(i,j,k)
                  if (grid%tslb(i,k,j) .lt. 170 .or. &
                       grid%tslb(i,k,j) .gt. 400) then
                     write(message,'(A,I4.4,A,I4.4,F5.2)') &
                          'update grid%tslb unreasonable, i = ',i,' j = ',j, &
                          grid%tslb(i,k,j)
                     call wrf_error_fatal(trim(message))
                  end if
                  
                  grid%smois(i,k,j) = this%SoilMoist(i,j,k)
                  grid%sh2o(i,k,j) = this%SoilWet(i,j,k)
!                  grid%smcrel(i,k,j) = this%RelSMC(i,j,k)
                  
               end do
            end do
         end do
      end if

   end subroutine setRealDataLandPoints

   !---------------------------------------------------------------------------

   subroutine setRealDataWaterPoints(this,grid,its,ite,jts,jte,ide,jde, &
        we_have_tsk,we_have_tavgsfc,internal_time_loop, &
        flag_sst,flag_pinterp,flag_ptheta,config_flags)
      implicit none

      ! Arguments
      type(LisData),intent(in) :: this
      type(domain),intent(inout) :: grid
      integer,intent(in) :: its,ite,jts,jte
      integer,intent(in) :: ide,jde
      logical,intent(in) :: we_have_tsk
      logical,intent(in) :: we_have_tavgsfc
      integer,intent(in) :: internal_time_loop
      integer,intent(in) :: flag_sst
      integer,intent(in) :: flag_pinterp
      integer,intent(in) :: flag_ptheta
      type(grid_config_rec_type),intent(in) :: config_flags

      ! Local variables
      integer :: i,j,k
      character(len=MAX_MESSAGE_LEN) :: message
      real :: xice_threshold
      real :: tsk_save
      integer :: ns
      real :: total_depth, mid_point_depth
      integer :: loop
      real, allocatable :: t2(:,:),t_gc(:,:,:)
      integer :: its_local,ite_local,jts_local,jte_local
      integer :: status

      ! Adjust 2-meter temperature based on REAL logic. We'll use internal
      ! arrays so REAL can repeat the procedure for all points.
      its_local = its
      jts_local = jts
      ite_local = min(ide-1,ite)
      jte_local = min(jde-1,jte)
      allocate(t2(its_local:ite_local,jts_local:jte_local),stat=status)
      if (status .ne. 0) then
         write(message,*) 'Error allocating t2'
         call wrf_error_fatal(message)
      end if
      allocate(t_gc(its_local:ite_local,1:grid%num_metgrid_levels, &
           jts_local:jte_local),stat=status)
      if (status .ne. 0) then
         write(message,*) 'Error allocating t_gc'
         call wrf_error_fatal(message)
      end if

      ! Local array copies
      do j = jts,min(jde-1,jte)
         do k = 1,grid%num_metgrid_levels
            do i = its,min(ide-1,ite)
               t_gc(i,k,j) = grid%t_gc(i,k,j)
            end do
         end do
      end do
      do j = jts,min(jde-1,jte)
         do i = its,min(ide-1,ite)
            t2(i,j) = grid%t2(i,j)
         end do
      end do

      ! Correction for data from P_INTERP
      if (flag_pinterp .eq. 1) then
         do j = jts,min(jde-1,jte)
            do i = its,min(ide-1,ite)
               t_gc(i,1,j) = t2(i,j)
            end do
         end do
      end if

      ! Correction for theta level data from UM
      if (flag_ptheta .eq. 1) then
         do j = jts,min(jde-1,jte)
            do i = its,min(ide-1,ite)
               t_gc(i,grid%num_metgrid_levels,j) = &
                    (t_gc(i,1,j) + t_gc(i,grid%num_metgrid_levels-1,j)) * 0.5
            end do
         end do
      end if

      ! Semi-fix missing surface field
      if (grid%p_gc(its,grid%num_metgrid_levels,jts) .lt. &
           grid%p_gc(its,2,jts)) then
         k = 2
      else
         k = grid%num_metgrid_levels
      end if
      if (t_gc(its,1,jts) .eq. -1.E30) then
         do j = jts,min(jde-1,jte)
            do i = its,min(ide-1,ite)
               t_gc(i,1,j) = t_gc(i,k,j)
            end do
         end do
      end if

      ! Final copy back to t2 array.
      do j = jts,min(jde-1,jte)
         do i = its,min(ide-1,ite)
            t2(i,j) = t_gc(i,1,j)
         end do
      end do
      
      ! Set xice_threshold
      if (config_flags%fractional_seaice .eq. 0) then
         xice_threshold = 0.5
      elseif (config_flags%fractional_seaice .eq. 1) then
         xice_threshold = 0.02
      else
         call wrf_error_fatal('Illegal value for fractional seaice!')
      end if

      ! Update 2D land data at Water points, except water temperatures
      do j = jts,min(jde-1,jte)
         do i = its,min(ide-1,ite)

            if (this%Landmask(i,j) .ge. 0.5) cycle ! Skip land points

            grid%landmask(i,j) = 0
            grid%xland(i,j) = 2
            grid%vegfra(i,j) = 0
            grid%shdmax(i,j) = 0
            grid%shdmin(i,j) = 0
            grid%canwat(i,j) = 0
            grid%soilcat(i,j) = this%soilWater
            grid%isltyp(i,j) =  this%soilWater
            grid%vegcat(i,j) = this%luWater
            grid%ivgtyp(i,j) = this%luWater
            grid%lu_index(i,j) = this%luWater
            ! FIXME See if LIS can provide fractional soil type fractions
            grid%soilctop(i,:,j) = 0
            grid%soilcbot(i,:,j) = 0
            grid%soilctop(i,this%soilWater,j) = 1
            grid%soilcbot(i,this%soilWater,j) = 1

            grid%snow(i,j) = 0
            grid%snowh(i,j) = 0
            grid%snowc(i,j) = 0

            ! Use LDT terrain if available.  This only happens if LDT was
            ! run with "Elevation fill option" set to "none"; otherwise,
            ! LDT sticks -9999 for the value.
            if (this%Elevation(i,j) > -9998.) then
               grid%ht(i,j) = this%Elevation(i,j)
            end if

            ! For albedo, use REAL default for water points
            grid%albbck(i,j) = 0.08
            grid%albedo(i,j) = 0.08

            ! For snow albedo, use REAL defaults for icy and non-icy
            ! points. Unfortunately LIS can't help us here, so we'll use 
            ! the XICE field from METGRID and hope for the best.
            if (grid%xice(i,j) > 0.5) then
               grid%snoalb(i,j) = 0.75
            else
               grid%snoalb(i,j) = 0.08
            end if
                        
         end do
      end do

      ! Set "soil layer" moisture values at water points
      do j = jts,min(jde-1,jte)
         do k = 1, this%numSoilLayers
            do i = its,min(ide-1,ite)
               
               if (this%Landmask(i,j) .ge. 0.5) cycle ! Skip land points
               
               grid%smois(i,k,j) = 1
               grid%sh2o(i,k,j) = 1
!               grid%smcrel(i,k,j) = 1
               
            end do
         end do
      end do

      ! Now handle water temperatures and sea ice. We'll mimic REAL, which is 
      ! very convoluted.
      CALL nl_get_iswater ( grid%id , grid%iswater )
      CALL nl_get_islake  ( grid%id , grid%islake )

      do j = jts,min(jde-1,jte)
         do i = its,min(ide-1,ite)

            if (this%Landmask(i,j) .ge. 0.5) cycle ! Skip land points

            ! Option to use average surface temperature for skin temperature.
            ! Typically this is only needed if skin temperature is missing.
            if (config_flags%use_tavg_for_tsk) then
               if (.not. we_have_tsk .and. .not. we_have_tavgsfc) then
                  call wrf_error_fatal( &
                       'We need either TSK or TAVGSFC over water from WPS')
               end if              
                if (we_have_tavgsfc) then
                  grid%tsk(i,j) = grid%tavgsfc(i,j)
               end if
            end if

            ! Replace water skin temperature, deep temperature, and
            ! elevation with METGRID values. LIS can't help us here.
            ! This should be run every time level, as METGRID fields
            ! are stored in the *_gc arrays when read in.
            grid%tsk(i,j) = grid%tsk_gc(i,j)
            grid%tmn(i,j) = grid%tmn_gc(i,j)
            grid%ht(i,j) = grid%ht_gc(i,j)  
            
            ! We now handle SST for inland lakes.  Typically the land use data
            ! doesn't discriminate between "lake" and "ocean", so we'll just
            ! update all water points. If SST is bad or missing, try average 
            ! surface temperature, skin temperature, and 2-meter temperature 
            ! in that order.
            if (grid%islake < 0) then
               grid%lakeflag=0

               if (grid%sst(i,j) .lt. 150 .or. &
                    grid%sst(i,j) .gt. 400) then
                  if (we_have_tavgsfc) then
                     grid%sst(i,j) = grid%tavgsfc(i,j)
                  end if
                  if (grid%sst(i,j) .lt. 150 .or. &
                       grid%sst(i,j) .gt. 400) then
                     grid%sst(i,j) = grid%tsk(i,j)
                  end if
                  if (grid%sst(i,j) .lt. 150 .or. &
                       grid%sst(i,j) .gt. 400) then
                     grid%sst(i,j) = t2(i,j) ! 2-meter temperature
                  end if
               end if
            else
               grid%lakeflag=1
               ! Here we do discriminate between "lake" and "ocean".  For 
               ! lake points, use average surface temperature, and fall back
               ! on 2-meter air temperature if necessary.
               if ( we_have_tavgsfc ) then
                  if ( grid%landusef(i,grid%islake,j) >= 0.5  .or. &
                       grid%lu_index(i,j) == grid%islake) then
                     ! EMK NUWRF...Logic change. Only do this if we lack
                     ! good SST.  That way we preserve, e.g., RSS SST for the
                     ! Great Lakes.
                     if (flag_sst .eq. 0 .or. &
                          (grid%sst(i,j) .lt. 150 .or. &
                           grid%sst(i,j) .gt. 400)) then
                        grid%sst(i,j) = grid%tavgsfc(i,j)
                        grid%tsk(i,j) = grid%tavgsfc(i,j)
                     end if
                     if (grid%sst(i,j) .lt. 150 .or. &
                          grid%sst(i,j) .gt. 400) then
                        grid%sst(i,j) = t2(i,j) ! 2-meter temperature
                     end if
                  end if
               else 
                  ! No average surface temperature available.  We redesignate
                  ! all lake points as ocean.
                  if ( grid%landusef(i,grid%islake,j) >= 0.5  .or. &
                       grid%lu_index(i,j) == grid%islake) then
                     grid%landusef(i,grid%iswater,j) = &
                          grid%landusef(i,grid%iswater,j) + &
                          grid%landusef(i,grid%islake,j)
                     grid%landusef(i,grid%islake,j) = 0.
                  end if
               end if
               ! We should still sanity check SST for ocean points.  
               ! Use skin temperature and 2-meter temperature in that order.
               if (grid%sst(i,j) .lt. 150 .or. &
                    grid%sst(i,j) .gt. 400) then
                  grid%sst(i,j) = grid%tsk(i,j)
               end if
               if (grid%sst(i,j) .lt. 150 .or. &
                    grid%sst(i,j) .gt. 400) then
                  grid%sst(i,j) = t2(i,j) ! 2-meter temperature
               end if
            end if
            
            ! Next, protect against bad skin temperature over water.
            tsk_save = grid%tsk(i,j) ! Save for sea ice adjustments
            if (flag_sst .eq. 1 .and. &
                 grid%sst(i,j) .gt. 170. .and. &
                 grid%sst(i,j) .lt. 400.) then
               grid%tsk(i,j) = grid%sst(i,j)
            end if

            ! Set tslb values following logic in process_soil_real and
            ! subordinate routines.
            do k = 1,this%numSoilLayers
               if (flag_sst .eq. 1) then
                  grid%tslb(i,k,j) = grid%sst(i,j)
               else
                  grid%tslb(i,k,j) = grid%tsk(i,j)
               end if
            end do

            ! Skip "pre" seaice adjustment step, as REAL only changes land
            ! points here (for Noah LSM).            

            ! Next, perform updates to "deep soil" and skin temperature.
            ! (In REAL, this occurs after adjust_for_seaice_pre)
            if ( (flag_sst .eq. 1) .and. &
                 (grid%sst(i,j) .gt. 170.) .and. &
                 (grid%sst(i,j) .lt. 400.) ) then
               grid%tmn(i,j) = grid%sst(i,j)
               grid%tsk(i,j) = grid%sst(i,j) 
            else
               grid%tmn(i,j) = grid%tsk(i,j)
            end if

            ! Adjust skin temperature again.
            ! Use 2-m temperature
            if (internal_time_loop .ne. 1) then
               if (grid%tsk(i,j) .lt. 170 .or. &
                   grid%tsk(i,j) .gt. 400.) then
                  ! FIXME...REAL uses the 2-meter temperature interpolated to
                  ! the model terrain. However, that interpolation occurs
                  ! after this LIS processing occurs. As a substitute, we'll
                  ! use the 2-meter temperature read in from the GRIB
                  ! file.
!                  grid%tsk(i,j) = grid%t_2(i,1,j)
                  grid%tsk(i,j) = t2(i,j)
               end if
            else
               if (grid%tsk(i,j) .lt. 170. .or. &
                   grid%tsk(i,j) .gt. 400.) then
                  if (grid%tmn(i,j) .gt. 170. .and. &
                       grid%tmn(i,j) .lt. 400.) then
                     grid%tsk(i,j)=grid%tmn(i,j)
                  else if (grid%sst(i,j).gt.170. .and. &
                           grid%sst(i,j).lt.400.) then
                     grid%tsk(i,j)=grid%sst(i,j)
                  else
                     CALL wrf_error_fatal ( 'grid%tsk unreasonable' )
                  end if
               end if
            end if

            ! Skip "tmn reasonable" check as REAL only applies it to land.
            
            ! Next, check tslb. Substitute with interpolated 2-m temperature
            if (internal_time_loop .ne. 1) then
               do ns = 1,this%numSoilLayers
                  ! FIXME...REAL uses the 2-meter temperature interpolated to
                  ! the model terrain. However, that interpolation occurs
                  ! after this LIS processing occurs. As a substitute, we'll
                  ! use the 2-meter temperature read in from the GRIB
                  ! file.
                  if (grid%tslb(i,ns,j) .LT. 170 .or. &
                      grid%tslb(i,ns,j) .GT. 400. ) THEN
!                     grid%tslb(i,ns,j) = grid%t_2(i,1,j)
                     grid%tslb(i,ns,j) = t2(i,j)
                  end if
               end do
            end if

            ! Another seaice adjustment
            ! Follow approach for NOAH LSM, minus land use, landmask, and
            ! soil type changes.
            if ( grid%xice(i,j) .gt. 200.0) grid%xice(i,j) = 0.0
            if ( (tsk_save .lt. 170 .or. tsk_save .gt. 400) .and. &
                 (grid%tsk(i,j) .gt. 170 .and. grid%tsk(i,j) .lt. 400) ) then
               tsk_save = grid%tsk(i,j)
            end if
            if ( (tsk_save .lt. 170 .or. tsk_save .gt. 400) .and. &
                 (grid%tsk(i,j) .lt. 170 .or. &
                 grid%tsk(i,j) .gt. 400) ) then
               call wrf_error_fatal('Bad skin temperature for sea ice!')
            else if (grid%xice(i,j) .ge. xice_threshold .or. &
                     tsk_save .lt. config_flags%seaice_threshold) then
               if (config_flags%fractional_seaice .eq. 0) grid%xice(i,j) = 1.0
               grid%tmn(i,j) = 271.4
               grid%tsk(i,j) = tsk_save
               total_depth = 3. ! Assume ice is 3 m deep and equispaced layers
               do loop = 1,this%numSoilLayers
                  mid_point_depth=(total_depth/this%numSoilLayers)*0.5 + &
                       (loop-1)*(total_depth/this%numSoilLayers)
                  grid%tslb(i,loop,j) = &
                       ( (total_depth-mid_point_depth)*tsk_save + &
                          mid_point_depth*grid%tmn(i,j) ) / total_depth
                  grid%smois(i,loop,j) = 1.0
                  grid%sh2o(i,loop,j) = 0.0
               end do
               grid%snoalb(i,j) = 0.75
            else if (grid%xice(i,j) .lt. xice_threshold) then
               grid%xice(i,j) = 0.
            end if

            ! Finally, Another adjustment for SST
            if (flag_sst .ne. 1) then
               grid%sst(i,j) = grid%tsk(i,j)
            end if

         end do
      end do

      message = "Deallocation error"
      deallocate(t2,stat=status)
      if (status .ne. 0) call wrf_error_fatal(trim(message))
      deallocate(t_gc,stat=status)
      if (status .ne. 0) call wrf_error_fatal(trim(message))

   end subroutine setRealDataWaterPoints

   !---------------------------------------------------------------------------
   ! Public
   !---------------------------------------------------------------------------

   subroutine checkTsk(grid,its,ite,jts,jte,ide,jde)

      implicit none

      ! Arguments
      type(domain),intent(inout) :: grid
      integer,intent(in) :: its,ite,jts,jte
      integer,intent(in) :: ide,jde

      ! Local variables
      integer :: i,j
      character(len=MAX_MESSAGE_LEN) :: message

      do j = jts,MIN(jde-1,jte)
         do i = its,MIN(ide-1,ite)
            if (grid%tsk(i,j) .lt. 170 .or. grid%tsk(i,j) .gt. 400) then
               write(message,'(A,I4.4,A,I4.4,A,F5.2)') &
                    'grid%tsk unreasonable, i = ',i,' j = ',j, &
                    'tsk = ',grid%tsk(i,j)
               call wrf_error_fatal(trim(message))
            end if
         end do
      end do

      return
   end subroutine checkTsk

   !---------------------------------------------------------------------------
   ! Public
   !---------------------------------------------------------------------------

   subroutine checkTslb(grid,its,ite,jts,jte,ide,jde,numSoilLayers)

      implicit none

      ! Arguments
      type(domain),intent(inout) :: grid
      integer,intent(in) :: its,ite,jts,jte
      integer,intent(in) :: ide,jde
      integer,intent(in) :: numSoilLayers

      ! Local variables
      integer :: i,j,k
      character(len=MAX_MESSAGE_LEN) :: message

      do j = jts,MIN(jde-1,jte)
         do k = 1,numSoilLayers
            do i = its,MIN(ide-1,ite)
               if (grid%tslb(i,k,j) .lt. 170 .or. &
                    grid%tslb(i,k,j) .gt. 400) then
                  write(message,'(A,I4.4,A,I4.4,A,I4.4,F5.2)') &
                       'grid%tslb unreasonable, i = ',i,' j = ',j, &
                       ' k = ',k,grid%tslb(i,k,j)
                  call wrf_error_fatal(trim(message))
               end if
            end do
         end do
      end do

   end subroutine checkTslb

   !---------------------------------------------------------------------------
   ! Public
   !---------------------------------------------------------------------------

   subroutine compareHtWithLIS(this,grid,its,ite,jts,jte,ide,jde)

      implicit none

      ! Arguments
      type(LisData),intent(in) :: this
      type(domain),intent(in) :: grid
      integer,intent(in) :: its,ite,jts,jte
      integer,intent(in) :: ide,jde

      ! Local variables
      real :: diff
      integer :: i,j

      ! Compare REAL terrain with LIS and land points
      do j = jts,MIN(jde-1,jte)
         do i = its,MIN(ide-1,ite)
            if (this%Landmask(i,j) .ge. 0.5) then
               diff = grid%ht(i,j) - this%Elevation(i,j)
               if (abs(diff) > 0) then
                  print*,'WARNING: i,j,ht,elevation = ', &
                       i,j,grid%ht(i,j),this%Elevation(i,j)
               end if
            endif
         end do
      end do
   end subroutine compareHtWithLIS

   !---------------------------------------------------------------------------
   ! Private subroutine for reading 2D array from LIS/LDT
   !---------------------------------------------------------------------------

   subroutine readLisRealArray2d(ncID,variableName,ims,ime,jms,jme, &
        array)

      implicit none

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(in) :: ims,ime,jms,jme ! WRF horizontal memory limits
      real, allocatable, intent(out) :: array(:,:)

      ! Local variables
      character(len=MAX_MESSAGE_LEN) :: message
      integer :: varID
      integer :: start(2)
      integer :: count(2)
      integer :: status

      ! Get variable ID
      varID = findVarId(ncID,trim(variableName),2,NF_FLOAT)

      ! Allocate array for WRF memory patch
      allocate(array(ims:ime,jms:jme),stat=status)
      if (status .ne. 0) then
         write(message,*)'ERROR allocating variable ',trim(variableName)
         call wrf_error_fatal(trim(message))
      end if

      ! Read array memory patch from file.
      start = (/ims,jms/)
      count = (/ime-ims+1, jme-jms+1/)
      status = nf_get_vara_real(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*)'ERROR reading variable ',trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

   end subroutine readLisRealArray2d

   !---------------------------------------------------------------------------
   ! Private subroutine for reading a 3D variable from LIS/LDT
   !---------------------------------------------------------------------------

   subroutine readLisRealArray3d(ncID,variableName,ims,ime,jms,jme, &
        array)

      implicit none

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(in) :: ims,ime,jms,jme ! WRF horizontal memory limits
      real, allocatable, intent(out) :: array(:,:,:)

      ! Local variables
      character(len=MAX_MESSAGE_LEN) :: message
      integer :: status
      integer :: varID
      integer :: dimIDS(3)
      integer :: ndims
      integer :: start(3)
      integer :: count(3)
      integer :: dim3

      ! Get variable ID
      varID = findVarID(ncID,variableName,3,NF_FLOAT)

      ! Get dimension IDs for variable. 
      status = nf_inq_vardimid(ncID,varID,dimIDs)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving dimension IDs for variable ', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
      
      ! Get length of third dimension.
      status = nf_inq_dimlen(ncID,dimIDs(3),dim3)
      if (status .ne. NF_NOERR) then
         write(message,*)'ERROR inquiring dimension, ', &
              trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Allocate array for WRF memory patch.
      allocate(array(ims:ime,jms:jme,dim3),stat=status)
      if (status .ne. 0) then
         write(message,*)'ERROR allocating variable ',trim(variableName)
         call wrf_error_fatal(trim(message))
      end if

      ! Read array from file.
      start=(/ims,jms,1/)
      count=(/ime-ims+1,jme-jms+1,dim3/)
      status = nf_get_vara_real(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*)'ERROR reading variable ',trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

   end subroutine readLisRealArray3d

   !---------------------------------------------------------------------------

   subroutine setRealLandUseFractions(this,grid,its,ite,jts,jte,ide,jde)

      ! Arguments
      type(LisData),intent(in) :: this
      type(domain),intent(inout) :: grid
      integer,intent(in) :: its,ite,jts,jte
      integer,intent(in) :: ide,jde

      ! Local variables
      integer :: status
      character(len=MAX_MESSAGE_LEN) :: message
      integer :: sm31,em31,sm33,em33
      integer :: i,j,k
      type(fieldlist),pointer :: p

      ! Update WRF/REAL landusef field.  We need to reallocate the array
      ! to match the number of land use categories provided by LIS.  Also,
      ! LIS only provides dominant categories, so the fractions here will be
      ! either zero or one.

      p => grid%head_statevars%next
      do while ( associated(p))
         if (trim(p%dataname) .eq. 'LANDUSEF') then
            call wrf_debug(300,'Changing vertical dimensions for landusef')
            p%ed2 = this%numLUCategories
            p%em2 = this%numLUCategories
            p%ep2 = this%numLUCategories

            deallocate(grid%landusef,stat=status)
            if (status .ne. 0) then
               message = 'Deallocation failure'
               call wrf_error_fatal(trim(message))
            end if

            sm31 = p%sm1
            em31 = p%em1
            sm33 = p%sm3
            em33 = p%em3

            allocate(grid%landusef(sm31:em31,1:this%numLUCategories, &
                 sm33:em33),  &
                 stat=status )
            if (status .ne. 0) then
               message = 'Allocation failure'
               call wrf_error_fatal(trim(message))
            end if

            p%rfield_3d => grid%landusef
            exit ! Get out of loop
         end if
         p => p%next
      end do

      do j = jts,MIN(jde-1,jte)
         do i = its,MIN(ide-1,ite)
            if (this%Landmask(i,j) .gt. 0.5) then ! Land points
               do k = 1, this%numLUCategories
                  grid%landusef(i,k,j) = this%Landusef(i,j,k)
               end do
            else ! Water points
               ! FIXME See if LIS can help us at water points here
               grid%landusef(i,:,j) = 0
               grid%landusef(i,this%luWater,j) = 1
            endif
         end do
      end do
   end subroutine setRealLandUseFractions

   !---------------------------------------------------------------------------

   subroutine reallocRealSoilLayers(this,grid,its,ite,jts,jte,ide,jde)

      ! Arguments
      type(LisData),intent(in) :: this
      type(domain),intent(inout) :: grid
      integer,intent(in) :: its,ite,jts,jte
      integer,intent(in) :: ide,jde

      ! Local variables
      integer :: status
      character(len=MAX_MESSAGE_LEN) :: message
      integer :: sm31,em31,sm33,em33
      integer :: i,j,k
      type(fieldlist),pointer :: p

      ! Reallocate soil layer arrays to have same number of layers as LIS.

      p => grid%head_statevars%next
      do while ( associated(p))
         if (trim(p%dataname) .eq. 'TSLB' .or. &
             trim(p%dataname) .eq. 'SMOIS' .or. &
             trim(p%dataname) .eq. 'SH2O' .or. &
             trim(p%dataname) .eq. 'SMCREL' .or. &
             trim(p%dataname) .eq. 'SOILCTOP' .or. & 
             trim(p%dataname) .eq. 'SOILCBOT' ) then
            call wrf_debug(300, &
                 'Changing vertical dimensions for '//trim(p%dataname))
            p%ed2 = this%numSoilLayers
            p%em2 = this%numSoilLayers
            p%ep2 = this%numSoilLayers

            if (trim(p%dataname) .eq. 'TSLB') then
               deallocate(grid%tslb,stat=status)
            else if (trim(p%dataname) .eq. 'SMOIS') then
               deallocate(grid%smois,stat=status)
            else if (trim(p%dataname) .eq. 'SH2O') then
               deallocate(grid%sh2o,stat=status)
            else if (trim(p%dataname) .eq. 'SMCREL') then
               deallocate(grid%smcrel,stat=status)
            else if (trim(p%dataname) .eq. 'SOILCTOP') then
               deallocate(grid%soilctop,stat=status)
            else if (trim(p%dataname) .eq. 'SOILCBOT') then
               deallocate(grid%soilcbot,stat=status)
            end if
            if (status .ne. 0) then
               message = 'Deallocation failure'
               call wrf_error_fatal(trim(message))
            end if

            sm31 = p%sm1
            em31 = p%em1
            sm33 = p%sm3
            em33 = p%em3
            
            if (trim(p%dataname) .eq. 'TSLB') then
               allocate(grid%tslb(sm31:em31,1:this%numSoilLayers,sm33:em33),  &
                 stat=status )
            else if (trim(p%dataname) .eq. 'SMOIS') then
               allocate(grid%smois(sm31:em31,1:this%numSoilLayers,sm33:em33),  &
                 stat=status )
            else if (trim(p%dataname) .eq. 'SH2O') then
               allocate(grid%sh2o(sm31:em31,1:this%numSoilLayers,sm33:em33),  &
                 stat=status )
            else if (trim(p%dataname) .eq. 'SMCREL') then
               allocate(grid%smcrel(sm31:em31,1:this%numSoilLayers,sm33:em33),  &
                 stat=status )
            else if (trim(p%dataname) .eq. 'SOILCTOP') then
!               allocate(grid%soilctop(sm31:em31,1:this%ns,sm33:em33),  &
!                 stat=status )
               allocate(grid%soilctop(sm31:em31,1:this%numSoilCategories, &
                    sm33:em33),stat=status )

            else if (trim(p%dataname) .eq. 'SOILCBOT') then
!               allocate(grid%soilcbot(sm31:em31,1:this%ns,sm33:em33),  &
!                 stat=status )
               allocate(grid%soilcbot(sm31:em31,1:this%numSoilCategories, &
                    sm33:em33),stat=status )

            end if

            if (status .ne. 0) then
               message = 'Allocation failure'
               call wrf_error_fatal(trim(message))
            end if

            if (trim(p%dataname) .eq. 'TSLB') then
               p%rfield_3d => grid%tslb
            else if (trim(p%dataname) .eq. 'SMOIS') then
               p%rfield_3d => grid%smois
            else if (trim(p%dataname) .eq. 'SH2O') then
               p%rfield_3d => grid%sh2o
            else if (trim(p%dataname) .eq. 'SMCREL') then
               p%rfield_3d => grid%smcrel
            else if (trim(p%dataname) .eq. 'SOILCTOP') then
               p%rfield_3d => grid%soilctop
            else if (trim(p%dataname) .eq. 'SOILCBOT') then
               p%rfield_3d => grid%soilcbot
            end if            
         end if
         p => p%next
      end do

      ! Just initialize the arrays. We'll substitute values elsewhere.
      do j = jts,min(jde-1,jte)
         do k = 1, this%numSoilLayers
            do i = 1,its,min(ide-1,ite)
               grid%tslb(i,k,j) = 0
               grid%smois(i,k,j) = 0
               grid%sh2o(i,k,j) = 0
!               grid%smcrel(i,k,j) = 0
               grid%soilctop(i,k,j) = 0
               grid%soilcbot(i,k,j) = 0
            end do
         end do
      end do
   end subroutine reallocRealSoilLayers

   !---------------------------------------------------------------------------

   subroutine setRealSoilDepths(this,grid)

      ! Arguments
      type(LisData),intent(in) :: this
      type(domain),intent(inout) :: grid

      ! Local variables
      integer :: status
      character(len=MAX_MESSAGE_LEN) :: message
      type(fieldlist),pointer :: p
      integer :: k

      p => grid%head_statevars%next
      do while ( associated(p))
         if (trim(p%dataname) .eq. 'ZS' .or. &
             trim(p%dataname) .eq. 'DZS') then
            call wrf_debug(300,'Changing vertical dimensions for soil layers')
            p%ed1 = this%numSoilLayers
            p%em1 = this%numSoilLayers
            p%ep1 = this%numSoilLayers

            if (trim(p%dataname) .eq. 'ZS') then
               deallocate(grid%zs,stat=status)
            else if (trim(p%dataname) .eq. 'DZS') then
               deallocate(grid%dzs,stat=status)
            end if
            if (status .ne. 0) then
               message = 'Deallocation failure'
               call wrf_error_fatal(trim(message))
            end if

            if (trim(p%dataname) .eq. 'ZS') then
               allocate(grid%zs(1:this%numSoilLayers),stat=status )
            else if (trim(p%dataname) .eq. 'DZS') then
               allocate(grid%dzs(1:this%numSoilLayers),stat=status )
            end if
            if (status .ne. 0) then
               message = 'Allocation failure'
               call wrf_error_fatal(trim(message))
            end if

            if (trim(p%dataname) .eq. 'ZS') then
               p%rfield_1d => grid%zs
            else if (trim(p%dataname) .eq. 'DZS') then
               p%rfield_1d => grid%dzs
            end if
         end if
         p => p%next
      end do

      do k = 1, this%numSoilLayers
         grid%zs(k) = this%zs(k)
         grid%dzs(k) = this%dzs(k)
      end do
   end subroutine setRealSoilDepths

   !---------------------------------------------------------------------------

   ! FIXME Have these values included in the LDT netCDF file
   subroutine setLisLandSoilCodes(this,lisLandcoverScheme)
      implicit none

      ! Arguments 
      type(LisData),intent(inout) :: this
      character(len=*),intent(in) :: lisLandcoverScheme

      ! Local variables
      character(len=MAX_MESSAGE_LEN) :: message

      this%soilWater = WATER_SOILTYPE
      this%numSoilCategories = NUM_CAT_SOIL

      if (trim(lisLandcoverScheme) .eq. LANDUSE_USGS_IN_LDT) then
         this%luWater = WATER_LANDUSE_USGS
         this%luSnowIce = SNOWICE_LANDUSE_USGS
         this%luUrban = URBAN_LANDUSE_USGS
         this%luClassSystem = LANDUSE_USGS_IN_WRF
         if (this%numLUCategories .ne. NUM_CAT_LANDUSE_USGS) then
            write(message,*) &
                 "ERROR, bad number of LIS veg types. Expected ", &
                 NUM_CAT_LANDUSE_USGS,', found ',this%numLUCategories
         end if
      else if (trim(lisLandcoverScheme) .eq. LANDUSE_MODIS_IN_LDT) then
         this%luWater = WATER_LANDUSE_MODIS
         this%luSnowIce = SNOWICE_LANDUSE_MODIS
         this%luUrban = URBAN_LANDUSE_MODIS
         this%luClassSystem = LANDUSE_MODIS_IN_WRF
         if (this%numLUCategories .ne. NUM_CAT_LANDUSE_MODIS) then
            write(message,*) &
                 "ERROR, bad number of LIS veg types. Expected ", &
                 NUM_CAT_LANDUSE_USGS,', found ',this%numLUCategories
         end if
      else if (trim(lisLandcoverScheme) .eq. LANDUSE_UMD_IN_LDT) then
         this%luWater = WATER_LANDUSE_UMD
         this%luSnowIce = SNOWICE_LANDUSE_UMD
         this%luUrban = URBAN_LANDUSE_UMD
         this%luClassSystem = LANDUSE_UMD_IN_WRF
         if (this%numLUCategories .ne. NUM_CAT_LANDUSE_UMD) then
            write(message,*) &
                 "ERROR, bad number of LIS veg types. Expected ", &
                 NUM_CAT_LANDUSE_USGS,', found ',this%numLUCategories
         end if
      else
         write(message,'(A,A)') &
              "Invalid value of lis_landcover_scheme:  ", &
              trim(lisLandcoverScheme)
         call wrf_error_fatal(message)
      end if

   end subroutine setLisLandSoilCodes

   !---------------------------------------------------------------------------

   subroutine setRealLandUseFlags(this,grid)

      implicit none

      ! Arguments
      type(LisData),intent(in) :: this
      type(domain),intent(inout) :: grid

      ! Reset flag for water
      call nl_set_iswater(grid%id,this%luWater)
      grid%iswater = this%luWater

      ! Reset islake flag to undefined. LIS doesn't have this info.
      call nl_set_islake(grid%id,-1)
      grid%islake = -1

      ! Reset flag for snow/ice
      call nl_set_isice(grid%id,this%luSnowIce)
      grid%isice = this%luSnowIce

      ! Reset flag for urban
      call nl_set_isurban(grid%id,this%luUrban)
      grid%isurban = this%luUrban

      ! Reset string for land classification system
      call nl_set_mminlu(grid%id,this%luClassSystem)
      grid%mminlu = this%luClassSystem

      ! Reset the number of land categories
      call nl_set_num_land_cat(grid%id,this%numLUCategories)
      grid%num_land_cat = this%numLUCategories

   end subroutine setRealLandUseFlags

end module ReconcileWithLis_mod
#else
module ReconcileWithLis_mod
end module ReconcileWithLis_mod
#endif
