!------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group / CISTO 606
!------------------------------------------------------------------------------
!
! MODULE:  interpolator_mod
!
! AUTHOR:
! Jossy P. Jacob (SSAI)
! May 2015
! Modified April 2016 Jossy Jacob
! 12/16 : Major mods (C. Cruz)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This module is for interpolating the CO2 fluxes to NU-WRF grid
! in relation to NU-WRF model state and to provide QC (Quality Control) to the 
! CO2 Flux output based on terrestrial CO2 regional flux computation. Components
! of CO2 fluxes from NPP (net production, monthly mean), RESP (respiratory flux, 
! monthly mean), Fire flux (daily), Fossil Fuel (annual mean) and Ocean CO2 
! (monthly mean) are combined to get the total flux. NPP and RESP components 
! are computed in relation to SW and T2 from NU-WRF state. 
!------------------------------------------------------------------------------

module flux_interpolator_mod

   use netcdf_util_mod
   use derived_variables_mod
   use sorted_datetimes_mod
   use interp_util_mod
   use casa_flux_filenames_mod
   use interpolator_mod

   implicit none
   private
   
   public :: new_flux_interpolator
   public :: destroy_flux_interpolator
   public :: run_flux_interpolator

   type wrfinput_grid1
      private
      integer :: nc_id
      integer :: wrf_x_dim
      integer :: wrf_y_dim
      integer :: wrf_z_dim
      real :: dx
      real :: dy
      character(len=MAX_DATE_LEN) :: datetime
      real,allocatable :: wrf_latitudes(:,:)
      real,allocatable :: wrf_longitudes(:,:)
   end type wrfinput_grid1
   type wrfbdy_grid1
      private
      integer :: nc_id
      integer :: wrf_x_dim
      integer :: wrf_y_dim
      integer :: wrf_z_dim
      integer :: wrf_t_dim
      integer :: wrf_w_dim
      character(len=MAX_DATE_LEN),allocatable :: datetimes(:)
      character(len=MAX_DATE_LEN),allocatable :: thisbdytime(:)
      character(len=MAX_DATE_LEN),allocatable :: nextbdytime(:)
      integer :: delta_time
   end type wrfbdy_grid1

! Define private data type for storing the CASA CO2 flux data
   public :: flux_interpolator
   type flux_interpolator
      private
      integer :: num_wrf_domains
      type(wrfinput_grid1), allocatable :: wrfinput(:)
      type(wrfbdy_grid1) :: wrfbdy
      character(len=MAX_DIR_LEN) :: wrf_dir
      character(len=MAX_DIR_LEN) :: flux_dir
      character(len=MAX_DIR_LEN) :: flux_prefix
      integer :: num_datetimes
      character(len=MAX_DATE_LEN),allocatable :: datetimes(:)
      logical,allocatable :: for_wrfbdy(:)
      logical,allocatable :: for_wrfinput(:,:)
      integer :: casa_lon_dim
      integer :: casa_lat_dim
      integer :: casa_lev_dim
      real,allocatable :: casa_longitudes(:)
      real,allocatable :: casa_latitudes(:)
      integer :: num_aerosols
      character(len=80),allocatable :: aerosol_names(:)
      character(len=80),allocatable :: aerosol_units(:)
      character(len=80),allocatable :: aerosol_descriptions(:)

      real, allocatable :: npp(:,:,:)  
      real, allocatable :: resp(:,:,:)  
      real, allocatable :: bfuel(:,:,:)  
      real, allocatable :: OceanCO2(:,:,:)  
      real, allocatable :: Fossilfuel(:,:,:)  
      real, allocatable :: fire(:,:,:,:)  
      real, allocatable :: CO2flux(:,:,:)  
   
   end type flux_interpolator
   type fileDetails
      !private
      character(len=132) :: filename
      character (len = 6) :: varname
      character (len = 4) :: year
      integer :: ndata
      integer :: nx
      integer :: ny 
      real, allocatable  :: lon(:)
      real, allocatable  :: lat(:)
      real, allocatable :: data1(:,:,:)
      character(len=MAX_DATE_LEN),allocatable :: data_times(:)
   end type fileDetails

   type nuwrf_grid
      character(len=MAX_DATE_LEN),allocatable :: wrf_datetimes(:)
      integer :: nx, ny, ntimes
      real :: dx, dy
      real, allocatable :: lon(:,:), lat(:,:), area(:,:),mask(:,:)
   end type nuwrf_grid
   type qc 
      real, allocatable :: flux_orig(:)
      real, allocatable :: flux_interp(:)
      real, allocatable :: flux_interp_time(:)
      real :: limit(4)
   end type qc
         
   integer, dimension(12) :: mday = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
   character(len=40) :: outfilename 
   real, parameter :: undfn = 1.e14
   
   ! Default moving average window
   integer :: window = 1
   ! Location of NU-WRF output
   character(len=MAX_DIR_LEN) :: wrf_out_dir
   
contains

!---------------------------------------------------------------------------
!
! ROUTINE:  new_flux_interpolator
!
! DESCRIPTION:  Public constructor method for interpolator data type.
!
!---------------------------------------------------------------------------
   function new_flux_interpolator(nml_c2w)
      use nml_mod
      use casa_variables_mod
      use casa_emission_variables_mod
      type (nml) :: nml_c2w
      type(flux_interpolator) :: new_flux_interpolator
      ! Local variables
      integer :: i
      character(len=256) :: wrf_dir
      integer :: num_domains
      character(len=256) :: flux_dir
      character(len=32) :: flux_prefix

      call nml_c2w%read("wrf","wrf_dir",wrf_dir)
      call nml_c2w%read("wrf","max_dom",num_domains)
      call nml_c2w%read("wrf","wrf_out_dir",wrf_out_dir,init=.true.)
      if (trim(wrf_out_dir) .eq. "MISSING") then
         wrf_out_dir = '.'
      end if
      call nml_c2w%read("wrf","moving_average_window",window)
      call nml_c2w%read("casa","flux_dir",flux_dir)
      call nml_c2w%read("casa","flux_prefix",flux_prefix)

      ! Initialize the data type
      call init_flux_interpolator(new_flux_interpolator)

      new_flux_interpolator%num_aerosols = NUM_CASA_eVARIABLES
      allocate(new_flux_interpolator%aerosol_names(NUM_CASA_eVARIABLES))
      allocate(new_flux_interpolator%aerosol_units(NUM_CASA_eVARIABLES))
      allocate(new_flux_interpolator%aerosol_descriptions(NUM_CASA_eVARIABLES))
      do i = 1, NUM_CASA_eVARIABLES
         new_flux_interpolator%aerosol_names(i) = trim(casa_wrf_variable_names(i))
         new_flux_interpolator%aerosol_units(i) = trim(casa_wrf_variable_units(i))
         new_flux_interpolator%aerosol_descriptions(i) = &
              trim(casa_wrf_variable_descriptions(i))
      end do
      
      ! Assemble space and time coordinate data from wrfinput and wrfbdy
      ! files.
      call prep_from_wrf_input(new_flux_interpolator, wrf_dir, num_domains) 
      call prep_sorted_datetimes(new_flux_interpolator, flux_dir, flux_prefix)

   end function new_flux_interpolator

!------------------------------------------------------------------------------
   subroutine init_flux_interpolator(this)
!------------------------------------------------------------------------------
      ! Arguments
      type(flux_interpolator), intent(inout) :: this

      ! Initialize immediate members
      this%num_wrf_domains = 0
      this%wrf_dir = trim('NULL')
      this%flux_dir = trim('NULL')
      this%flux_prefix = trim('NULL')
      this%num_datetimes = 0
      this%casa_lon_dim = 0
      this%casa_lat_dim = 0
      this%casa_lev_dim = 0
      this%num_aerosols = 0

      ! Initialize wrfbdy members
      this%wrfbdy%nc_id = 0
      this%wrfbdy%wrf_x_dim = 0
      this%wrfbdy%wrf_y_dim = 0
      this%wrfbdy%wrf_z_dim = 0
      this%wrfbdy%wrf_t_dim = 0
      this%wrfbdy%wrf_w_dim = 0
      this%wrfbdy%delta_time = 0

   end subroutine init_flux_interpolator

!---------------------------------------------------------------------------
!
! ROUTINE:  destroy_flux_interpolator
!
! DESCRIPTION:  Public destructor method for flux_interpolator data type.
!
!---------------------------------------------------------------------------
   subroutine destroy_flux_interpolator(this)

      ! Arguments
      type(flux_interpolator), intent(inout) :: this

      ! Local variables
      integer :: i

      ! Free memory associated with wrfinput data
      do i = 1, this%num_wrf_domains
         if (allocated(this%wrfinput(i)%wrf_latitudes)) &
              deallocate(this%wrfinput(i)%wrf_latitudes)
         if (allocated(this%wrfinput(i)%wrf_longitudes)) &
              deallocate(this%wrfinput(i)%wrf_longitudes)

         ! If netCDF file is open, close it.
         call close_netcdf_file_if_open(this%wrfinput(i)%nc_id)
      end do
      if (allocated(this%wrfinput)) deallocate(this%wrfinput)

      ! If netCDF file is open, close it.
      call close_netcdf_file_if_open(this%wrfbdy%nc_id)

      ! Free memory associated with wrfbdy data
      if (allocated(this%wrfbdy%datetimes)) &
           deallocate(this%wrfbdy%datetimes)
      if (allocated(this%wrfbdy%thisbdytime)) &
           deallocate(this%wrfbdy%thisbdytime)
      if (allocated(this%wrfbdy%nextbdytime)) &
           deallocate(this%wrfbdy%nextbdytime)
      
      ! Miscellaneous
      if (allocated(this%datetimes)) deallocate(this%datetimes)
      if (allocated(this%for_wrfbdy)) deallocate(this%for_wrfbdy)
      if (allocated(this%for_wrfinput)) deallocate(this%for_wrfinput)
      if (allocated(this%casa_longitudes)) &
           deallocate(this%casa_longitudes)
      if (allocated(this%casa_latitudes)) &
           deallocate(this%casa_latitudes)

      if (allocated(this%aerosol_names)) deallocate(this%aerosol_names)
      if (allocated(this%aerosol_units)) deallocate(this%aerosol_units)
      if (allocated(this%aerosol_descriptions)) &
           deallocate(this%aerosol_descriptions)

      ! Take care of non-allocatable data members
      call init_flux_interpolator(this)

   end subroutine destroy_flux_interpolator

!---------------------------------------------------------------------------
!
! ROUTINE:  prep_from_wrf_input
!
! DESCRIPTION:  Private method used to assemble time and space coordinates
! from wrfinput and wrfbdy files.
!
!---------------------------------------------------------------------------
   
   subroutine prep_from_wrf_input(this, wrf_dir, num_wrf_domains)

      ! Arguments
      type(flux_interpolator),intent(inout) :: this
      ! Directory with wrfinput and wrfbdy files
      character(len=*),intent(in) :: wrf_dir 
      ! Number of wrf domains to process
      integer, intent(in) :: num_wrf_domains


      character(len=MAX_DATE_LEN) :: datetimes1(1)

      character(len=132) :: filename
      character(len=2) :: char_int
      integer :: ncid
      integer :: wrf_x_dim,wrf_y_dim,wrf_z_dim,wrf_t_dim,wrf_w_dim
      integer :: wrf_domain_id
      integer :: i,j
      real :: ddx,ddy
      real,allocatable :: xlat(:,:,:)
      real,allocatable :: xlong(:,:,:)

      write(6,'(1x,A,I2.2,A)') '[FLUX_INTERP] SUB prep_from_wrf_input'

      !  Get data from wrfinput files
      if (num_wrf_domains < 1) then
         print*,'ERROR, illegal number of WRF domains!'
         print*,'num_wrf_domains should be positive, but equals ', &
              num_wrf_domains
         stop
      end if
      this%num_wrf_domains = num_wrf_domains
      allocate(this%wrfinput(num_wrf_domains))

      this%wrf_dir = trim(wrf_dir)

      do wrf_domain_id = 1, num_wrf_domains

         write(6,'(1x,A,I2.2,A)') &
              '[FLUX_INTERP]    preparing wrfinput_d',wrf_domain_id

         ! Open the wrfinput netCDF file in .  
         write(char_int,"(I2.2)") wrf_domain_id
         filename = trim(wrf_dir)//"/wrfinput_d"//trim(char_int)
         ncid = open_netcdf_readfile(filename)
         this%wrfinput(wrf_domain_id)%nc_id = ncid

         ! Get required dimensions from file
         wrf_x_dim = read_netcdf_dimension(ncid,"west_east")
         wrf_y_dim = read_netcdf_dimension(ncid,"south_north")
         wrf_z_dim = read_netcdf_dimension(ncid,"bottom_top")
         this%wrfinput(wrf_domain_id)%wrf_x_dim = wrf_x_dim
         this%wrfinput(wrf_domain_id)%wrf_y_dim = wrf_y_dim
         this%wrfinput(wrf_domain_id)%wrf_z_dim = wrf_z_dim
        ! ! Get wrfinput dx and dy in km
         call read_netcdf_real_global_attribute(ncid,"DX",ddx)
         call read_netcdf_real_global_attribute(ncid, "DY",ddy)
         this%wrfinput(wrf_domain_id)%dx = ddx/1000.0  !m to km
         this%wrfinput(wrf_domain_id)%dy = ddy/1000.0  !meter to Km
         ! Copy latitude and longitude data from wrfinput
         !call copy_latlon_from_wrfinput(this,ncid,wrf_domain_id)
         allocate(xlat(wrf_x_dim,wrf_y_dim,1))      ! At base of column
         allocate(xlong(wrf_x_dim,wrf_y_dim,1))     ! At base of column
         call read_netcdf_real_array_3d(ncid,"XLAT",wrf_x_dim,wrf_y_dim,1,xlat)
         call read_netcdf_real_array_3d(ncid,"XLONG",wrf_x_dim,wrf_y_dim,1,xlong)
         allocate(this%wrfinput(wrf_domain_id)%wrf_latitudes(wrf_x_dim,wrf_y_dim))
         allocate(this%wrfinput(wrf_domain_id)%wrf_longitudes(wrf_x_dim, &
           wrf_y_dim))
         do j = 1, wrf_y_dim
         do i = 1, wrf_x_dim
            this%wrfinput(wrf_domain_id)%wrf_latitudes(i,j) = xlat(i,j,1)
            this%wrfinput(wrf_domain_id)%wrf_longitudes(i,j) = xlong(i,j,1)
         end do
         end do
         deallocate(xlat)
         deallocate(xlong)

         ! Read datetime from wrfinput
         call read_netcdf_character_array_1d(ncid,"Times",MAX_DATE_LEN,1, &
              datetimes1)
         this%wrfinput(wrf_domain_id)%datetime = trim(datetimes1(1))
          ! We only have to process wrfbdy for domain 1, so if we're working
         ! with a different domain, we're done.  Deallocate the remaining 
         ! temporary arrays and move on to the next domain.
         if (wrf_domain_id .ne. 1) then
            cycle
         end if

         write(6,'(1x,A)')'[FLUX_INTERP]    preparing wrfbdy_d01'

         ! Open the wrfbdy file
         filename = trim(wrf_dir)//"/wrfbdy_d"//trim(char_int)
         ncid = open_netcdf_readfile(filename)
         this%wrfbdy%nc_id = ncid
         
         ! Get the required dimensions from the file
         wrf_x_dim = read_netcdf_dimension(ncid,"west_east")
         wrf_y_dim = read_netcdf_dimension(ncid,"south_north")
         wrf_z_dim = read_netcdf_dimension(ncid,"bottom_top")
         wrf_t_dim = read_netcdf_dimension(ncid,"Time")
         wrf_w_dim = read_netcdf_dimension(ncid,"bdy_width")
         write(6,'(1x,A,1x,i5)')'[FLUX_INTERP]    wrf_t_dim for bdy files:',wrf_t_dim
         this%wrfbdy%wrf_x_dim = wrf_x_dim
         this%wrfbdy%wrf_y_dim = wrf_y_dim
         this%wrfbdy%wrf_z_dim = wrf_z_dim
         this%wrfbdy%wrf_t_dim = wrf_t_dim
         this%wrfbdy%wrf_w_dim = wrf_w_dim
         ! Sanity check the x, y, and z-dimensions
         !call compare_wrfinput_wrfbdy_dims(this)
 ! Sanity check the x, y, and z-dimensions
         if (this%wrfbdy%wrf_x_dim .ne. this%wrfinput(1)%wrf_x_dim) then
           print*,'ERROR, x-dimension mismatch between wrfbdy and wrfinput!'
           print*,'wrfinput x-dimension:  ',this%wrfinput(1)%wrf_x_dim
           print*,'wrfbdy x-dimension:  ',this%wrfbdy%wrf_x_dim
         end if
         if (this%wrfbdy%wrf_y_dim .ne. this%wrfinput(1)%wrf_y_dim) then
           print*,'ERROR, y-dimension mismatch between wrfbdy and wrfinput!'
           print*,'wrfinput y-dimension:  ',this%wrfinput(1)%wrf_y_dim
           print*,'wrfbdy y-dimension:  ',this%wrfbdy%wrf_y_dim
         end if
         if (this%wrfbdy%wrf_z_dim .ne. this%wrfinput(1)%wrf_z_dim) then
           print*,'ERROR, z-dimension mismatch between wrfbdy and wrfinput!'
           print*,'wrfinput z-dimension:  ',this%wrfinput(1)%wrf_z_dim
           print*,'wrfbdy z-dimension:  ',this%wrfbdy%wrf_z_dim
         end if

         ! Copy datetimes from wrfbdy
         allocate(this%wrfbdy%datetimes(wrf_t_dim))

         call read_netcdf_character_array_1d(ncid,"Times",MAX_DATE_LEN, &
              wrf_t_dim, this%wrfbdy%datetimes)
         allocate(this%wrfbdy%thisbdytime(wrf_t_dim))
         call read_netcdf_character_array_1d(ncid, &
              "md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", &
              MAX_DATE_LEN, wrf_t_dim, this%wrfbdy%thisbdytime)
         allocate(this%wrfbdy%nextbdytime(wrf_t_dim))
         call read_netcdf_character_array_1d(ncid, &
              "md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", &
              MAX_DATE_LEN, wrf_t_dim, this%wrfbdy%nextbdytime)
         ! Sanity check datetimes
        ! call compare_wrfbdy_wrfinput_time(this)
         if (trim(this%wrfbdy%datetimes(1)) .ne. &
           trim(this%wrfinput(1)%datetime)) then
         print*,'ERROR, time stamps differ between wrfinput and wrfbdy!'
         print*,'wrfinput                    ', &
              trim(this%wrfinput(1)%datetime)
         print*,'wrfbdy (first time level):  ', &
              trim(this%wrfbdy%datetimes(1))
         stop
         end if

         ! Calculate seconds between wrfbdy time levels
         this%wrfbdy%delta_time = &
              calc_delta_time_wrfbdy(this%wrfbdy%thisbdytime(1), &
              this%wrfbdy%nextbdytime(1))
      end do
      write(6,'(1x,A,I2.2,A)') '[FLUX_INTERP] END SUB prep_from_wrf_input'
      
   end subroutine prep_from_wrf_input

!---------------------------------------------------------------------------
!
! ROUTINE:  prep_sorted_datetimes (for flux emission data)
!
! DESCRIPTION: Private method for assembling horizontally interpolated
! CASA pressure values on WRF grid.  Targets GEOS-5 netCDF4 files.
!
!---------------------------------------------------------------------------

   subroutine prep_sorted_datetimes(this,flux_dir,flux_prefix)

      ! Arguments
      type(flux_interpolator),intent(inout) :: this
      ! Directory with CASA netCDF files
      character(len=*),intent(in) :: flux_dir
      ! Prefix for CASA netCDF filenames
      character(len=*),intent(in) :: flux_prefix

      ! Local variables
      type(sorted_datetimes) :: datetime_object
      integer :: ii

      this%flux_dir = trim(flux_dir)
      this%flux_prefix = trim(flux_prefix)

      ! Assemble sorted list of datetimes.  Save information on which
      ! wrf files need which datetimes.
      datetime_object = create_sorted_datetimes(this%num_wrf_domains)
      do ii = 1, this%wrfbdy%wrf_t_dim
         call add_to_sorted_datetimes(datetime_object, &
              this%wrfbdy%datetimes(ii),for_wrfbdy=.true., &
              for_wrfinput=.false., wrf_domain_id=1)
      end do
      ! Need one more datetime to calculate wrfbdy tendencies
      call add_to_sorted_datetimes(datetime_object, &
           this%wrfbdy%nextbdytime(this%wrfbdy%wrf_t_dim),for_wrfbdy=.true., &
           for_wrfinput=.false., wrf_domain_id=1)
      do ii = 1, this%num_wrf_domains
         call add_to_sorted_datetimes(datetime_object, &
              this%wrfinput(ii)%datetime, for_wrfbdy=.false., &
              for_wrfinput=.true., wrf_domain_id=ii)
      end do
      call get_sorted_datetimes(datetime_object, &
           this%num_datetimes,this%num_wrf_domains, &
           this%datetimes, this%for_wrfbdy, this%for_wrfinput)
      call destroy_sorted_datetimes(datetime_object)

!!$      write(6,'(1x,a,1x,i5)')'[FLUX_INTERP] number of datetimes', &
!!$           this%num_datetimes
!!$      do jj = 1, this%num_datetimes
!!$         print*,'Datetime:  ',this%datetimes(jj)
!!$      !   write(6,'(A,L1)') 'For wrfbdy_d01:  ',this%for_wrfbdy(jj)
!!$      end do

   end subroutine prep_sorted_datetimes

!-----------------------------------------------------------------------
!
! ROUTINE:  run_flux_interpolator  for emissions data
!
! DESCRIPTION:  Using the interpolator data type and information on
! variable names and conversion factors, reads/interpolates/writes
! an aerosol/chemistry species from CASA emission to cheminput file.
!
!---------------------------------------------------------------------------

   subroutine run_flux_interpolator(this, casa_varname, wrf_varname, &
        flux_dt)

      ! Arguments
      type(flux_interpolator),intent(inout) :: this
      character(len=*),intent(in) :: casa_varname
      character(len=*),intent(in) :: wrf_varname
      real :: conversion_factor, &
              conversion_factor_day, &
              conversion_factor_mon, &
              conversion_factor_year
      integer             :: dt_day, num_yeardays

      ! Local variables
      character(len=256) :: flux_outfile
      character(len=MAX_DATE_LEN):: wrf_datetime
      character(len=MAX_DATE_LEN):: wrf_datetime1
      character(len=MAX_DATE_LEN):: wrf_datetime2
      real, allocatable :: tmp_flux(:,:,:)
      real, allocatable :: resp_hr(:,:)
      real, allocatable :: npp_hr(:,:)
      real, allocatable :: bfuel_hr(:,:)
      real, allocatable :: fire_hr(:,:)
      real, allocatable :: Fossil_hr(:,:)
      real, allocatable :: OceanCO2_hr(:,:)
      real, allocatable :: flux_hr(:,:),veg_hr(:,:),gpp_hr(:,:)
      real, allocatable :: flux_hr_last(:,:)
      real, allocatable :: veg_flux_last(:,:), &
                           fire_flux_last(:,:), &
                           oceanCO2_hr_last(:,:), &
                           Fossil_hr_last(:,:)
!!$      real,allocatable :: Finalflux(:,:,:)
      real,allocatable :: Finalflux(:,:,:),veg_flux(:,:),fire_flux(:,:)
      real,allocatable :: E_flux(:,:), ET_flux(:,:)
 
      real,allocatable :: npp_interp(:,:,:) 
      real,allocatable :: resp_interp(:,:,:) 
      real,allocatable :: fire_interp(:,:,:,:) 
      real,allocatable :: BFuel_interp(:,:,:) 
      real,allocatable :: FossilFuel_interp(:,:,:) 
      real,allocatable :: OceanCO2_interp(:,:,:) 

      integer :: lat_dim,lon_dim,time_dim 
      character(len=MAX_DATE_LEN),allocatable :: wrf_times(:)
      character(len=MAX_DATE_LEN) :: wrf_times_last
      real, allocatable :: wrf_temp(:,:,:)
      real, allocatable :: wrf_lon(:,:)
      real, allocatable :: wrf_lat(:,:)
      real, allocatable :: wrf_sw(:,:,:),sw_mon(:,:,:),sw(:,:)
      real, allocatable :: wrf_t2(:,:,:),t2(:,:)
      real, allocatable :: wrf_q10(:,:,:),q10_mon(:,:,:),q10(:,:)
      real, allocatable :: tmp(:,:)
      real, allocatable :: t2_aver(:,:,:), sw_aver(:,:,:)
      real, allocatable :: sum_sw(:,:), sum_t2(:,:)

      real :: flux_dt
      integer :: wrf_x_dim, wrf_y_dim
      integer :: ncid
      integer :: ii,nn,i,j,nmonth,it,nj1,itt
      character(len=132) :: tmp_file
      character(len=2) :: dom_id
      integer :: pyear, pmonth, pday, phour

      type (fileDetails) :: NPPfile, FIREfile, RESPfile,&
           BFUELfile, OceanCO2file, Fossilfile 
      type (nuwrf_grid) :: wrf_grid1
! these are for computing the total flux
      type (qc) :: qc_npp, qc_resp, qc_fossil, qc_fire, qc_bfuel, qc_oco2

      character(len=4) :: cyear
      character(len=2) :: cmonth,cday,chour,cminute,csecond
      real :: fac_day
      character(len=MAX_DATE_LEN):: mydate
      integer :: num_fluxtimes, num_dailyfluxtimes, ntime_year
      integer :: year, month, day, hour,qcval, qtime
      real ::  sum_area_wrf, wrf_param(5)
!!$      real, allocatable :: qc_total_flux(:), qc_veg_flux(:)

      character(len=:), allocatable :: date_id
      character(len=:), allocatable :: CO2_datetimes(:)
      character(len=:), allocatable :: CO2_dailydatetimes(:)
      character(len=:), allocatable :: year_datetimes(:)

      write(6,'(1x,a)') '[FLUX_INTERP] SUB run_flux_interpolator'
      
      ! Get the year
      wrf_datetime = this%datetimes(1) 
      read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear,cmonth,cday,chour,cminute,csecond
      
! Conversion Factors
      conversion_factor = 3.0e11 !!Kg/m2/sec to mol/Km^2 /Sec 
      fac_day=86400. !seconds      
      conversion_factor_day = 1.e-3/fac_day
      
      do ii = 1, this%num_wrf_domains !ii = 1 ! for outer domain only
         
         if(.not. allocated(wrf_grid1%wrf_datetimes)) &
            allocate(wrf_grid1%wrf_datetimes(this%num_datetimes)) 
         wrf_grid1%wrf_datetimes = this%datetimes
         wrf_grid1%nx = this%wrfinput(ii)%wrf_x_dim
         wrf_grid1%ny = this%wrfinput(ii)%wrf_y_dim 
         wrf_grid1%dx = this%wrfinput(ii)%dx
         wrf_grid1%dy = this%wrfinput(ii)%dy 

         wrf_x_dim = this%wrfinput(ii)%wrf_x_dim
         wrf_y_dim =  this%wrfinput(ii)%wrf_y_dim
         
         if(.not. allocated(wrf_grid1%lon)) allocate(wrf_grid1%lon(wrf_grid1%nx,wrf_grid1%ny))
         if(.not. allocated(wrf_grid1%lat)) allocate(wrf_grid1%lat(wrf_grid1%nx,wrf_grid1%ny))
         if(.not. allocated(wrf_grid1%area)) allocate(wrf_grid1%area(wrf_grid1%nx, wrf_grid1%ny)) 
         if(.not. allocated(wrf_grid1%mask)) allocate(wrf_grid1%mask(wrf_grid1%nx, wrf_grid1%ny)) 
         wrf_grid1%lon = this%wrfinput(ii)%wrf_longitudes 
         wrf_grid1%lat = this%wrfinput(ii)%wrf_latitudes
         ! WRF grid 27km equal area grids 
         wrf_grid1%area = wrf_grid1%DX*wrf_grid1%DY 

         if(.not. allocated(tmp)) allocate(tmp(wrf_grid1%nx, wrf_grid1%ny)) 
         tmp = 0
         qcval = 0
         do i= 1, wrf_grid1%nx
            do j = 1, wrf_grid1%ny
               if ((wrf_grid1%lon(i,j) >= qclimit(3)) .and. &
                    (wrf_grid1%lon(i,j) <= qclimit(4)) .and. &
                    (wrf_grid1%lat(i,j) >= qclimit(1)) .and. &
                    (wrf_grid1%lat(i,j) <= qclimit(2))) then 
                  tmp(i,j)=1
                  qcval = qcval +1 
               end if
            end do
         end do      
         wrf_grid1%mask = tmp
         sum_area_wrf = sum(wrf_grid1%area, mask = (wrf_grid1%mask .eq. 1))

! Read the NPP data file NPP_2010_monthly.nc    
         NPPfile%filename = NPPname
         qc_npp%limit(:) = qclimit(:)
         if(allocated(this%npp)) deallocate(this%npp)
         call read_interp_myfile3d(NPPfile, this%flux_dir, this%npp, &
              wrf_grid1, qc_npp)

! Read the RESP data file RESP_2010_monthly.nc 
         RESPfile%filename = RESPname
         qc_resp%limit(:) = qclimit(:)
         if(allocated(this%resp)) deallocate(this%resp)
         call read_interp_myfile3d(RESPfile, this%flux_dir, this%resp, &
              wrf_grid1, qc_resp)
         time_dim = RESPfile%ndata

! Read the BioFUEL data file BioFuel_2010_monthly.nc
         BFUELfile%filename = BFUELname 
         qc_bfuel%limit(:) = qclimit(:)
         if(allocated(this%bfuel)) deallocate(this%bfuel)
         call read_interp_myfile3d(BFUELfile, this%flux_dir, this%bfuel, &
              wrf_grid1, qc_bfuel)
         time_dim = BFUELfile%ndata

! Read the OCEANCO2 data file Ocean_CO2_2010_monthly.nc  
         OceanCO2file%filename = OceanCO2name
         qc_oco2%limit(:) = qclimit(:)
         if(allocated(this%OceanCO2)) deallocate(this%OceanCO2)
         call read_interp_myfile3d(OceanCO2file,this%flux_dir, this%OceanCO2, &
              wrf_grid1, qc_oco2)
         time_dim = OceanCO2file%ndata

! Read the Fossil Fuel data file FossilFuel_2010_clim.nc
! (Fossil Fuel is on different grid)    
         Fossilfile%filename = FossilFuelname 
         qc_fossil%limit(:) = qclimit(:)
         if(allocated(this%FossilFuel)) deallocate(this%FossilFuel)
         call read_interp_myfile3d(Fossilfile,this%flux_dir, this%Fossilfuel, &
              wrf_grid1, qc_fossil)
         time_dim = Fossilfile%ndata

! Read the Wild FIRE data file FIRE_2010_daily_01.nc
! Daily data/ 1 month data in one file 
         nmonth = 12
         Firefile%filename = FIREname 
         qc_fire%limit(:) = qclimit(:)
         if(allocated(this%fire)) deallocate(this%fire)
         call read_interp_myfile4d(FIREfile,this%flux_dir, this%fire, &
              wrf_grid1, qc_fire)
      
         write(6,'(1x,a)')'[FLUX_INTERP]    write qc_flux.dat'
         open(unit=11, file='qc_flux.dat', form='unformatted',convert='big_endian')
         write(11) qc_npp%flux_orig
         write(11) qc_npp%flux_interp
         write(11) qc_resp%flux_orig
         write(11) qc_resp%flux_interp
         write(11) qc_bfuel%flux_orig
         write(11) qc_bfuel%flux_interp
         write(11) qc_oco2%flux_orig
         write(11) qc_oco2%flux_interp
         write(11) qc_fossil%flux_orig
         write(11) qc_fossil%flux_interp 
         write(11) qc_fire%flux_orig
         write(11) qc_fire%flux_interp
         close(11)
         
! Write the interpolated fields
         write(6,'(1x,a)')'[FLUX_INTERP]    write the Interp_file_(s)'
         outfilename = 'Interp_file_npp.nc'
         call write_test_netcdf_file(wrf_x_dim,wrf_y_dim, &
              12,1, this%npp, "npp",trim(outfilename))
         outfilename = 'Interp_file_resp.nc'
         call write_test_netcdf_file(wrf_x_dim,wrf_y_dim, &
              12,1, this%resp, "resp",trim(outfilename))
         outfilename = 'Interp_file_Fossil_fuel.nc'
         call write_test_netcdf_file(wrf_x_dim,wrf_y_dim, &
              1,1, this%FossilFuel, "FossilFuel",trim(outfilename))
         outfilename = 'Interp_file_Bio_fuel.nc'
         call write_test_netcdf_file(wrf_x_dim,wrf_y_dim, &
              12,1, this%BFuel, "BioFuel",trim(outfilename))
         outfilename = 'Interp_file_OceanCO2.nc'
         call write_test_netcdf_file(wrf_x_dim,wrf_y_dim, &
              12,1, this%OceanCO2, "OceanCO2",trim(outfilename))
         outfilename = 'Interp_file_Fire.nc'
         call write_test_netcdf_file(wrf_x_dim,wrf_y_dim, &
              31,12, this%fire, "WildFire",trim(outfilename))

! Now read the WRF hourly output from daily files
!
! First create a list of datetimes from first to last wrfout filename
!      
         write(6,'(1x,a)')'[FLUX_INTERP]    read WRF hourly output from daily files'
         wrf_datetime1 = this%datetimes(1) 
         wrf_datetime2 = this%datetimes(this%num_datetimes)
         write(6,'(1x,4(a))')'[FLUX_INTERP]    start-end ',&
              trim(wrf_datetime1),' - ',trim(wrf_datetime2)

         dt_day = 1 !daily fluxdata output
         date_id = 'hourly' 
         ! Looks like CO2_datetimes are not used
!!$         call get_datetimes(wrf_datetime, wrf_datetime2, date_id, int(flux_dt), &
!!$              num_fluxtimes, CO2_datetimes)
!!$         write(6,'(1x,a,1x,i6)') &
!!$              '[FLUX_INTERP]    num_fluxtimes: ',num_fluxtimes
         date_id = 'daily' 
         call get_datetimes(wrf_datetime1, wrf_datetime2, date_id, dt_day, &
              num_dailyfluxtimes, CO2_dailydatetimes)
         write(6,'(1x,a,1x,i6)') &
              '[FLUX_INTERP]    num_dailyfluxtimes: ',num_dailyfluxtimes
         date_id = 'yearly' 
         call get_datetimes(wrf_datetime1, wrf_datetime2, date_id, dt_day, &
              ntime_year, year_datetimes)
         write(6,'(1x,a,1x,i6)') &
              '[FLUX_INTERP]    ntime_year: ',ntime_year
       
         ! construct a filename
         wrf_datetime = CO2_dailydatetimes(1)
         write(6,'(1x,a,1x,a)') &
              '[FLUX_INTERP]   wrf_datetime : ',wrf_datetime
         read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") & 
              cyear,cmonth,cday,chour,cminute,csecond
         write(dom_id,"(I2.2)") ii
         tmp_file = trim(wrf_out_dir)//'/wrfout_d'//dom_id//'_'//wrf_datetime 
         write(6,'(1x,3(a))')'[FLUX_INTERP]    get dims from: <',trim(tmp_file),'>'
         
         ncid = open_netcdf_readfile(trim(tmp_file))
         time_dim = read_netcdf_dimension(ncid,"Time")
         lon_dim = read_netcdf_dimension(ncid,"west_east")
         lat_dim = read_netcdf_dimension(ncid,"south_north")
         write(6,'(1x,a,3(1x,i6))')'[FLUX_INTERP]    t, x, y: ',&
              time_dim, lon_dim, lat_dim
         
         if(.not. allocated(wrf_temp)) allocate(wrf_temp(lon_dim,lat_dim,1))
         if(.not. allocated(wrf_lon)) allocate(wrf_lon(lon_dim,lat_dim))
         if(.not. allocated(wrf_lat)) allocate(wrf_lat(lon_dim,lat_dim))
         if(.not. allocated(wrf_times)) allocate(wrf_times(time_dim))
         wrf_lon(:,:) = wrf_temp(:,:,1); 
         wrf_lat(:,:) = wrf_temp(:,:,1); 
         call close_netcdf_file(ncid)
!
! Compute the monthly sums for SW and Q10 from T2m from NUWRF hourly output
!
         allocate(q10_mon(lon_dim,lat_dim,12))
         allocate(sw_mon(lon_dim,lat_dim,12))
         write(6,'(1x,a)')&
              '[FLUX_INTERP]    compute monthly sums for SW and Q10 from WRFOUT files'

         call get_sw_q10_mon(year_datetimes, dom_id, sw_mon, q10_mon, ntime_year)
!
! Now read the wrfout files for 1 year and compute the flux
!
         allocate (wrf_sw(lon_dim,lat_dim,time_dim))
         allocate (wrf_t2(lon_dim,lat_dim,time_dim))
         allocate (wrf_q10(lon_dim,lat_dim,time_dim))
         allocate (npp_hr(lon_dim,lat_dim))
         allocate (gpp_hr(lon_dim,lat_dim))
         allocate (resp_hr(lon_dim,lat_dim))
         allocate (bfuel_hr(lon_dim,lat_dim))
         allocate (Fossil_hr(lon_dim,lat_dim))
         allocate (OceanCO2_hr(lon_dim,lat_dim))
         
         allocate (q10(lon_dim,lat_dim))
         allocate (sw(lon_dim,lat_dim))
         allocate (t2(lon_dim,lat_dim))
         allocate (t2_aver(lon_dim,lat_dim,time_dim))
         allocate (sw_aver(lon_dim,lat_dim,time_dim))
         allocate (sum_sw(lon_dim,lat_dim))
         allocate (sum_t2(lon_dim,lat_dim))
         allocate (fire_hr(lon_dim,lat_dim))
         allocate (flux_hr(lon_dim,lat_dim))
         allocate (flux_hr_last(lon_dim,lat_dim))
         allocate (veg_flux_last(lon_dim,lat_dim))
         allocate (fire_flux_last(lon_dim,lat_dim))
         allocate (OceanCO2_hr_last(lon_dim,lat_dim))
         allocate (Fossil_hr_last(lon_dim,lat_dim))
         allocate (E_flux(lon_dim,lat_dim))
         allocate (ET_flux(lon_dim,lat_dim))
         allocate (veg_hr(lon_dim,lat_dim))
         allocate (veg_flux(lon_dim,lat_dim))
         allocate (fire_flux(lon_dim,lat_dim))     
         allocate (tmp_flux(lon_dim, lat_dim, time_dim))
         
         wrf_param(1)=wrf_grid1%dx * 1000.
         wrf_param(2)=wrf_grid1%dy * 1000.
         flux_hr_last=0.
         veg_flux_last=0.
         fire_flux_last=0.
         OceanCO2_hr_last=0.
         Fossil_hr_last=0.
         nj1=0

         write(6,'(1x,a,1x,i3,a)', advance='no') &
              '[FLUX_INTERP]    Loop over ', num_dailyfluxtimes,' days...'
         call flush(6)
         do nn = 1, num_dailyfluxtimes  ! desired daily data
           
            wrf_datetime = CO2_dailydatetimes(nn)
            read(wrf_datetime,"(i4,X,i2,X,i2,x,i2)") year,month,day,hour

            tmp_file = trim(wrf_out_dir)//'/wrfout_d'//dom_id//'_'//wrf_datetime
            
            if (nn < num_dailyfluxtimes) then
               write(6,'(1x,i3)', advance='no')nn
            else
               write(6,'(1x,i3)')nn
            end if
            call flush(6)
            
            ncid = open_netcdf_readfile(trim(tmp_file))
            time_dim = read_netcdf_dimension(ncid,"Time")
            call read_netcdf_character_array_1d(ncid,"Times", &
                 MAX_DATE_LEN, time_dim, wrf_times)
            call read_netcdf_real_array_3d(ncid, "T2", &
                 lon_dim, lat_dim, time_dim,wrf_t2)
            call read_netcdf_real_array_3d(ncid, "SWDOWN", &
                 lon_dim, lat_dim, time_dim,wrf_sw) 
            call close_netcdf_file(ncid)
!
!     Compute the Total flux
!
            do it = 1, time_dim  - window + 1 ! time_dim = 24 hours
               
               nj1 = nj1 + 1         
               mydate = wrf_times(it)
               read(mydate,"(i4,X,i2,X,i2,x,i2)") pyear,pmonth,pday,phour

               num_yeardays = 365
                 if (mod(pyear,400) == 0 .or. (mod(pyear,4) == 0 &
                    .and. (mod(pyear,100) .ne. 0))) then          
                  mday(2) = 29
                  num_yeardays = 366
               end if
               
               conversion_factor_mon = 1.e-3/(real(mday(pmonth))*fac_day)      
               conversion_factor_year = 1.e-3/(real(num_yeardays)*fac_day)  
               npp_hr = this%npp(:,:,pmonth) * conversion_factor_mon
               resp_hr = (this%npp(:,:,pmonth)+this%resp(:,:,pmonth)) &
                    * conversion_factor_mon
               bfuel_hr = this%bfuel(:,:,pmonth) * conversion_factor_mon
               
               ! Fossil Fuel data can be monthly or yearly. 
               if (size(this%FossilFuel,3) == 1) then  
                  Fossil_hr = this%Fossilfuel(:,:,1) *conversion_factor
               else
                  Fossil_hr = this%Fossilfuel(:,:,pmonth)*conversion_factor
               end if
               OceanCO2_hr = this%OceanCO2(:,:,pmonth) *conversion_factor
               fire_hr = this%fire(:,:,day,pmonth)* conversion_factor_day
               
               sw = 0.0
               q10 = 0.0
               
               if (window > 1) then ! window is set in namelist
                  call compute_sw_q10_using_aver()
               else
                  call compute_sw_q10_using_inst()
               end if
                              
               sw = sw * 24.* mday(pmonth)
               q10 = q10 * 24. * mday(pmonth)
               gpp_hr = 2.* sw * npp_hr 
               resp_hr = q10 * resp_hr
               veg_flux = (resp_hr - gpp_hr)*conversion_factor
               fire_flux = (bfuel_hr + fire_hr)*conversion_factor
               flux_hr = veg_flux + fire_flux + Fossil_hr + OceanCO2_hr 
               tmp_flux(:,:,it) = flux_hr 
               E_flux = flux_hr_last
               ET_flux = 0.
               
               ! ET_flux computes for every time except first
               if ( it > 1 ) then
                  !Flux tendency per hr
                  ET_flux=(flux_hr(:,:)-flux_hr_last(:,:))/flux_dt  
               end if

               ! write the Total flux data hourly
               if ( it > 1 ) then   
                  flux_outfile = &
                       trim(casa_varname)//'_d'//dom_id//'_'//wrf_times_last
                  call write_allflux_netcdf_file3d(wrf_x_dim,wrf_y_dim, &
                       1,wrf_times_last,flux_hr_last,ET_flux, &
                       veg_flux_last, fire_flux_last, oceanCO2_hr_last, &
                       Fossil_hr_last, &
                       trim(wrf_varname),trim(this%aerosol_units(1)), &
                       trim(this%aerosol_descriptions(1)), &
                       'ET_CASACO2','Flux Tendency', &
                       'Veg_flux', 'Vegetation Flux=npp+resp+fire+biofuel', &
                       'Fire_flux', 'Wild fire + bio fuel', &
                       'OceanCO2', 'Ocean CO2', &
                       'FossilFuel', 'Fossil Fuel', &
                       flux_outfile,wrf_param)  
               end if
               
               flux_hr_last=flux_hr
               veg_flux_last=veg_flux
               fire_flux_last=fire_flux
               OceanCO2_hr_last=OceanCO2_hr
               Fossil_hr_last=Fossil_hr
               wrf_times_last = wrf_times(it)
               
            end do ! it
            
         end do ! nn
         
         deallocate(wrf_sw)
         deallocate(wrf_t2)
         deallocate(wrf_q10)
         deallocate (npp_hr)
         deallocate (gpp_hr)
         deallocate (resp_hr)
         deallocate (bfuel_hr)
         deallocate (Fossil_hr)
         deallocate (OceanCO2_hr)
         deallocate (q10)
         deallocate (sw)
         deallocate (t2)
         deallocate (t2_aver)
         deallocate (sw_aver)
         deallocate (sum_sw)
         deallocate (sum_t2)
         deallocate (fire_hr)
         deallocate (veg_hr)
         deallocate (flux_hr)
         deallocate (flux_hr_last)
         deallocate (veg_flux_last)
         deallocate (fire_flux_last)
         deallocate (oceanCO2_hr_last)
         deallocate (Fossil_hr_last)
         deallocate (veg_flux)
         deallocate (fire_flux)
         deallocate (E_flux)
         deallocate (ET_flux)
         deallocate (tmp_flux)
         deallocate (q10_mon)
         deallocate (sw_mon)
         
      end do ! ii for domain
      
      write(6,'(1x,a)') '[FLUX_INTERP] END SUB run_flux_interpolator'
      
   contains

      subroutine compute_sw_q10_using_inst()
                          
         t2 = wrf_t2(:,:,it)-273.
         where(t2 > -30.)
            q10 = 1.5**((t2-30.)/10.)
         end where
               
         where ((abs(sw_mon(:,:,pmonth)) > 0.) &
              .and. (abs(wrf_sw(:,:,it)) < undfn))
            sw = wrf_sw(:,:,it)/sw_mon(:,:,pmonth)
         end where
               
         where (q10_mon(:,:,pmonth) > 0.)
            q10 = q10/q10_mon(:,:,pmonth)
         end where
                  
 
      end subroutine compute_sw_q10_using_inst
     
      subroutine compute_sw_q10_using_aver()
         integer :: iw
         
         iw=(window-1)/2
         q10 = 0.
         t2_aver = 0.
         sw_aver = 0.
         sum_sw = 0.
         sum_t2 = 0.

         ! --- compute centered moving average ---
         ! We align each moving average with the midpoint of the data
         ! that it averages. The edges of the the series are handled
         ! separately by averaging neighboring pairs.
         if (it < iw+1) then
            sw_aver(:,:,it) = (wrf_sw(:,:,it)+wrf_sw(:,:,it+1))/2.
            t2_aver(:,:,it) = (wrf_t2(:,:,it)+wrf_t2(:,:,it+1))/2.
         else if (it > time_dim-window) then
            sw_aver(:,:,it) = (wrf_sw(:,:,it)+wrf_sw(:,:,it-1))/2.
            t2_aver(:,:,it) = (wrf_t2(:,:,it)+wrf_t2(:,:,it-1))/2.
         else
            do itt = it-iw, it-iw + window - 1
               sum_sw = sum_sw + wrf_sw(:,:,itt)
               sum_t2 = sum_t2 + wrf_t2(:,:,itt)
            end do
            sw_aver(:,:,it) = sum_sw/window
            t2_aver(:,:,it) = sum_t2/window
         end if

!!$         ! --- compute moving average ---   
!!$         do itt = it, it + window - 1
!!$            sum_sw = sum_sw + wrf_sw(:,:,itt)
!!$            sum_t2 = sum_t2 + wrf_t2(:,:,itt)
!!$         end do
!!$         sw_aver(:,:,it) = sum_sw/window
!!$         t2_aver(:,:,it) = sum_t2/window
           
         ! --- update sw, t2 using moving average --
         t2 = t2_aver(:,:,it)-273.
         where(t2 > -30.)
            q10 = 1.5**((t2-30.)/10.)
         end where
               
         where ((abs(sw_mon(:,:,pmonth)) > 0.) &
              .and. (abs(sw_aver(:,:,it)) < undfn))
            sw = sw_aver(:,:,it)/sw_mon(:,:,pmonth)
         end where
               
         where (q10_mon(:,:,pmonth) > 0.)
            q10 = q10/q10_mon(:,:,pmonth)
         end where
                          
      end subroutine compute_sw_q10_using_aver
      
   end subroutine run_flux_interpolator
   
!------------------------------------------------------------------------------
   subroutine get_sw_q10_mon(year_datetimes, dom_id, sw_mon, q10_mon, ntimes)
!------------------------------------------------------------------------------
      character(len=*), intent(in) :: year_datetimes(:)
      character(len=2), intent(in) :: dom_id
      real, intent(inout) :: sw_mon(:,:,:), q10_mon(:,:,:)
      integer, intent(in) :: ntimes
!      
      integer :: imon, time_dim,lat_dim, lon_dim, it, itt
      integer :: year, month, day, hour,ncid, iy, nt
      character(len=MAX_DATE_LEN):: wrf_startdate, wrf_datetime
      character(len=256) :: tmp_file
      character(len=MAX_DATE_LEN), allocatable :: wrf_times(:) 
      real, allocatable :: wrf_sw(:,:,:),wrf_t2(:,:,:)
      real, allocatable :: t2(:,:), q10(:,:)
      real, allocatable :: t2_aver(:,:,:), sw_aver(:,:,:)
      real, allocatable :: sum_sw(:,:), sum_t2(:,:)

      write(6,'(1x,a,a)')'[FLUX_INTERP] SUB get_sw_q10_mon'
      if (window > 1) then
         write(6,'(1x,a,1x,i3)')&
              '[FLUX_INTERP]    computing moving average, window = ',window
      end if
      
      wrf_startdate = year_datetimes(1)
      ! Use first-date file to get dimensions
      tmp_file = trim(wrf_out_dir)//'/wrfout_d'//dom_id//'_'//wrf_startdate 
      write(6,'(1x,a,a,a)')'[FLUX_INTERP]    read start date <',trim(tmp_file),'>'
      call flush(6)
      
      ncid = open_netcdf_readfile(trim(tmp_file))
      time_dim = read_netcdf_dimension(ncid,"Time")
      lon_dim = read_netcdf_dimension(ncid,"west_east")
      lat_dim = read_netcdf_dimension(ncid,"south_north")     
      call close_netcdf_file(ncid)
      
! Read the wrfout files for 1 year and compute the sw_mon and q10
      allocate(wrf_times(time_dim))
      allocate(wrf_sw(lon_dim,lat_dim,time_dim))
      allocate(wrf_t2(lon_dim,lat_dim,time_dim))
      allocate(q10(lon_dim,lat_dim))
      allocate(t2(lon_dim,lat_dim))
      allocate(t2_aver(lon_dim,lat_dim,time_dim))
      allocate(sw_aver(lon_dim,lat_dim,time_dim))
      allocate(sum_sw(lon_dim,lat_dim))
      allocate(sum_t2(lon_dim,lat_dim))
      
      write(6,'(1x,a)',advance='no')'[FLUX_INTERP]    Loop over months...'
      call flush(6)

      ! CC: This is inefficient. We should loop only over the "needed" data
      ! (while keeping in mind the moving average window).
      do imon = 1, 12
         
         if (imon<12) then
            write(6,'(1x,i3)',advance='no')imon
         else
            write(6,'(1x,i3)')imon
         end if
         call flush(6)
         
         sw_mon(:,:,imon) = 0.
         q10_mon(:,:,imon) = 0.
         
         do iy = 1, ntimes
            
            wrf_datetime = year_datetimes(iy)

            read(wrf_datetime,"(i4,x,i2,x,i2,x,i2,6x)") &                    
                 year,month,day,hour
            if (month == imon ) then 
               tmp_file = &
                    trim(wrf_out_dir)//'/wrfout_d'//dom_id//'_'//wrf_datetime 
               ncid = open_netcdf_readfile(trim(tmp_file))
               call read_netcdf_character_array_1d(ncid,"Times",MAX_DATE_LEN, &
                    time_dim, wrf_times)
               call read_netcdf_real_array_3d(ncid, "T2", &
                    lon_dim, lat_dim, time_dim, wrf_t2)
               call read_netcdf_real_array_3d(ncid, "SWDOWN", &
                    lon_dim, lat_dim, time_dim, wrf_sw) 
               call close_netcdf_file(ncid)
               
               if (window > 1) then ! window is set in namelist
                  call compute_sums_using_aver()
               else
                  call compute_sums_using_inst()
               end if
                              
            end if
            
         end do ! iy
      end do ! imon
      
      deallocate(wrf_times)
      deallocate(wrf_sw)
      deallocate(wrf_t2)
      deallocate(q10)
      deallocate(t2)
      deallocate(t2_aver)
      deallocate(sw_aver)
      deallocate(sum_sw)
      deallocate(sum_t2)
      
      write(6,'(1x,a,a)')'[FLUX_INTERP] END SUB get_sw_q10_mon'
      
   contains

      subroutine compute_sums_using_inst()
         do it = 1, time_dim  ! nu-wrf time dim (hourly)
                  
            q10 = 0.
            where ( abs(wrf_sw(:,:,it)) < undfn)
               sw_mon(:,:,imon) = sw_mon(:,:,imon) + wrf_sw(:,:,it)
            end where
            t2 = wrf_t2(:,:,it)-273.
            
            where(t2 > -30.)
               q10 = 1.5**((t2-30.)/10.)
            end where
            q10_mon(:,:,imon) = q10_mon(:,:,imon) + q10(:,:)
                  
         end do
      end subroutine compute_sums_using_inst
      
      subroutine compute_sums_using_aver()
         integer :: iw
         
         do it = 1, time_dim
                  
            iw=(window-1)/2
            q10 = 0.
            t2_aver = 0.
            sw_aver = 0.
            sum_sw = 0.
            sum_t2 = 0.
            ! --- compute centered moving average ---
            ! We align each moving average with the midpoint of the data
            ! that it averages. The edges of the the series are handled
            ! separately by averaging neighboring pairs.
            if (it < iw+1) then
               sw_aver(:,:,it) = (wrf_sw(:,:,it)+wrf_sw(:,:,it+1))/2.
               t2_aver(:,:,it) = (wrf_t2(:,:,it)+wrf_t2(:,:,it+1))/2.
            else if (it > time_dim-window) then
               sw_aver(:,:,it) = (wrf_sw(:,:,it)+wrf_sw(:,:,it-1))/2.
               t2_aver(:,:,it) = (wrf_t2(:,:,it)+wrf_t2(:,:,it-1))/2.
            else
               do itt = it-iw, it-iw + window - 1
                  sum_sw = sum_sw + wrf_sw(:,:,itt)
                  sum_t2 = sum_t2 + wrf_t2(:,:,itt)
               end do
               sw_aver(:,:,it) = sum_sw/window
               t2_aver(:,:,it) = sum_t2/window
            end if

            ! --- update sw_mon, t2 using moving average --
            where ( abs(sw_aver(:,:,it)) < undfn)
               sw_mon(:,:,imon) = sw_mon(:,:,imon) + sw_aver(:,:,it)
            end where
            t2 = t2_aver(:,:,it)-273.
            
            where(t2 > -30.)
               q10 = 1.5**((t2-30.)/10.)
            end where
            q10_mon(:,:,imon) = q10_mon(:,:,imon) + q10(:,:)
                 
         end do
         
      end subroutine compute_sums_using_aver
         
   end subroutine get_sw_q10_mon
   
!------------------------------------------------------------------------------
   subroutine get_datetimes(wrf_startdate, wrf_enddate, id, dt, &
                            ntime, date_times)
!------------------------------------------------------------------------------
      character(len=*), intent(in) :: wrf_startdate, wrf_enddate 
      character(len=*), intent(in) :: id
      integer, intent(in) :: dt
      integer, intent(out) :: ntime
      character*(:), allocatable,intent(out):: date_times(:)
!      
      character(len=256) :: stime(1024)
      integer :: iyear, imon, iday, ihr, it,hour,i
      integer :: year1, month1, day1, hour1,minute1, second1
      integer :: year2, month2, day2, hour2,minute2, second2
      integer :: year, month, day
      integer :: nt, istart, iend
      integer, dimension(12) :: mday = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
      character(len=2) :: cmonth, cday, chour
      character(len=2) :: cmonth1, cday1, chour1,cminute1, csecond1
      character(len=2) :: cmonth2, cday2, chour2,cminute2, csecond2
      character(len=4) :: cyear1, cyear2, cyear

      read(wrf_startdate,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear1,cmonth1,cday1,chour1,cminute1,csecond1
      read(wrf_startdate,"(i4,X,i2,X,i2,x,i2,x,i2,X,i2)") &                    
           year1,month1,day1,hour1,minute1,second1
      read(wrf_enddate,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear2,cmonth2,cday2,chour2,cminute2,csecond2
      read(wrf_enddate,"(i4,X,i2,X,i2,x,i2,x,i2,X,i2)") &                    
           year2,month2,day2,hour2,minute2,second2
      if (year1 /= year2 ) then 
        print *, '*** ERROR, Year1 and Year2 are different', year1, year2
        return
     endif
     
! compute for one year and then select the ones between start time and end time. 
!   
!  Assign the number of days in Feb for Leap Year      
      iyear = year1
      if (mod(iyear,400) == 0 .or. &
         (mod(iyear,4) == 0 .and. &
         (mod(iyear,100) .ne. 0))) then          
         mday(2) = 29
      endif
      iyear = year1
      ntime = 0
      select case (id)
      case('hourly')
         call get_hourly()
      case('daily')
         call get_daily()
      case('yearly')
         call get_yearly()
      end select
      
   contains
      
!------------------------------------------------------------------------------
      subroutine get_hourly()
         do imon = month1,  month2
            do iday = 1, mday(imon) 
               do ihr = 1,24,dt
                  ntime = ntime + 1
               end do
            end do
         end do
         !allocate (character(MAX_DATE_LEN)::stime(ntime))
         nt=1
         do imon = month1,  month2
            do iday = 1, mday(imon) 
               do ihr = 1,24,dt
                  hour = ihr -1 
                  it = iyear*1000000 + imon *10000 + iday * 100 + hour 
                  write(cyear,'(I4)') iyear
                  write(cmonth,'(I2.2)') imon
                  write(cday,'(I2.2)') iday
                  write(chour,'(I2.2)') hour
                  stime(nt) = cyear//'-'//cmonth//'-'//cday//'_'//chour//':00:00'
                  nt = nt+1
               end do
            end do
         end do
         ! Now select the time you need based on WRFbdy file
         do i = 1, ntime
            if (wrf_startdate == stime(i)) istart = i
            if (wrf_enddate == stime(i)) iend = i
         end do
         ntime = iend-istart +1 
         if(.not. allocated(date_times)) &
                 allocate (character(ntime)::date_times(ntime)) 
         date_times = stime(istart:iend)
      end subroutine get_hourly
   
!------------------------------------------------------------------------------
      subroutine get_daily()
         do imon = month1,  month2
            do iday = 1, mday(imon), dt
               ntime = ntime + 1
            end do
         end do
         !allocate (character(MAX_DATE_LEN)::stime(ntime))
         nt=1
         do imon = month1,  month2
            do iday = 1, mday(imon),dt
               hour = 1 
               it = iyear*1000000 + imon *10000 + iday * 100 + hour 
               write(cyear,'(I4)') iyear
               write(cmonth,'(I2.2)') imon
               write(cday,'(I2.2)') iday
               write(chour,'(I2.2)') hour
               stime(nt) = cyear//'-'//cmonth//'-'//cday//'_'//chour//':00:00'
               nt = nt+1
            end do
         end do
         ! Now select the time you need based on WRFbdy file
         do i = 1, ntime
            read(stime(i),"(i4,X,i2,X,i2)")year,month,day
            if ((month == month1) .and. (day == day1)) istart = i
            if ((month == month2) .and. (day == day2)) iend = i
         end do
         if (hour == 1) istart = istart -1 
         ntime = iend-istart +1 
         if(.not. allocated(date_times)) &
                 allocate (character(ntime)::date_times(ntime)) 
         date_times = stime(istart:iend)
      end subroutine get_daily
      
!------------------------------------------------------------------------------
      subroutine get_yearly()
         do imon = 1,  12
            do iday = 1, mday(imon) 
               ntime = ntime + 1
            end do
         end do
         !allocate (character(MAX_DATE_LEN)::stime(ntime))
         nt=1
         hour=1
         do imon = 1,  12
            do iday = 1, mday(imon), dt 
               it = iyear*1000000 + imon *10000 + iday * 100 + hour 
               write(cyear,'(I4)') iyear
               write(cmonth,'(I2.2)') imon
               write(cday,'(I2.2)') iday
               write(chour,'(I2.2)') hour
               stime(nt) = cyear//'-'//cmonth//'-'//cday//'_'//chour//':00:00'
               nt = nt+1
            end do
         end do
         if(.not. allocated(date_times)) &
                 allocate (character(ntime)::date_times(ntime)) 
         date_times = stime
      end subroutine get_yearly
      
   end subroutine get_datetimes


!!---------------------------------------------------------------------------
   subroutine  read_interp_myfile3d(this, vardir, var_interp,wrf_grid1,qc1)
!!---------------------------------------------------------------------------
      type (fileDetails), intent(inout) :: this 
      type (nuwrf_grid), intent(in) :: wrf_grid1
      type (qc), intent(inout) :: qc1
      character(len=*), intent(in) ::  vardir
      real, allocatable, intent(out) :: var_interp(:,:,:)
      integer :: lon_dim, lat_dim, time_dim
      integer :: nn  
      character(len=4) :: cyear
      character(len=2) :: cmonth,cday,chour,cminute,csecond
      real, allocatable :: horiz_interp(:,:,:,:), &
                           var_emission(:,:,:)
      real, allocatable :: temp_lon(:), temp_lon2(:)
      real, allocatable :: temp_var(:,:,:), temp_var2(:,:)
      real :: inverse_delta_lon, inverse_delta_lat,lon0,dx,lat0,dy
      character(len=MAX_DATE_LEN):: wrf_datetime !, wrf_datetime2 
      character(len = 20) :: varname
      character(len = 30) :: myfilename
      character(len = 132) :: tmp_file
      integer :: i,j, wrf_x_dim, wrf_y_dim,nt,n1,ncid
      integer :: iq1, iq2, jq1, jq2, nval
      integer, allocatable :: tmp(:,:),temp(:,:),mask1(:,:)
      real, allocatable :: area_orig(:,:), area_wrf(:,:),tmp1(:,:), tmp2(:,:)
      real :: sum_area_wrf, sum_area_orig
   
      
      ! Get the year 
      wrf_datetime = wrf_grid1%wrf_datetimes(1)
      read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear,cmonth,cday,chour,cminute,csecond 
      myfilename = this%filename
      tmp_file = trim(vardir)//"/"//myfilename(1:5)//cyear//myfilename(10:30)
      write(6,'(1x,3(a))') &
           '[FLUX_INTERP] SUB read_interp_myfile3d <',trim(tmp_file),'>'

      ncid = open_netcdf_readfile(trim(tmp_file))
      time_dim = read_netcdf_dimension(ncid,"time")
      lon_dim = read_netcdf_dimension(ncid,'lon')
      lat_dim = read_netcdf_dimension(ncid,'lat')
      allocate(this%data_times(time_dim))
      allocate(this%lon(lon_dim))
      allocate(this%lat(lat_dim))
      allocate(this%data1(lon_dim,lat_dim,time_dim))
      this%ndata = time_dim
      this%varname = myfilename(1:4)
      this%year = cyear
      this%nx = lon_dim
      this%ny = lat_dim
      
      !dy = 0.5 ! Hard coded for emission data 0.5x0.5deg resolution
      !dx = 0.5
      !lat0 = -89.75
      if (myfilename(1:3) == 'npp' .or. myfilename(1:3) == 'NPP' ) then 
          varname = "NPP" 
      elseif (myfilename(1:3) == 'res' .or. myfilename(1:3) == 'RES' ) then 
          varname = "RespCO2" 
      elseif (myfilename(1:3) == 'bfu' .or. myfilename(1:3) == 'BFU' ) then 
          varname = "BioFuel" 
      elseif (myfilename(1:3) == 'ffu' .or. myfilename(1:3) == 'FFU' ) then 
          varname = "FossilFuel" 
         ! dx = 1.25 ! Hard coded for emission data 1.25 x1.0 deg resolution
         ! dy = 1.0
         ! lat0 = -89.5
      elseif (myfilename(1:3) == 'oco' .or. myfilename(1:3) == 'OCO' ) then 
          varname = "OceanCO2"
          !dx = 1.25 ! Hard coded for emission data 1.25 x1.0 deg resolution
         ! dy = 1.0 
          !lat0 = -89.5
      elseif (myfilename(1:3) == 'fir' .or. myfilename(1:3) == 'FIR' ) then 
          varname = "WildFire" 
      endif

      call read_netcdf_character_array_1d(ncid,"Times",MAX_DATE_LEN, &
              time_dim, this%data_times)
      call read_netcdf_real_array_1d(ncid,'lon',lon_dim,this%lon)
      call read_netcdf_real_array_1d(ncid,'lat',lat_dim,this%lat)
      call read_netcdf_real_array_3d(ncid,trim(varname),lon_dim,lat_dim, &
           time_dim,this%data1)
      dx=this%lon(2)-this%lon(1)
      dy=this%lat(2)-this%lat(1)
      lon0=this%lon(1)
      lat0=this%lat(1)
!!$      write(6,'(1x,a,1x,4(f12.4))')'[FLUX_INTERP]    dx, dy, lon0, lat0: ',&
!!$           dx, dy, lon0, lat0
      call close_netcdf_file(ncid)
      
!! Rotate the data if the longitude 0 to 360 (Expected longitude -180 to +180)
      if (this%lon(1) >= 0 ) then 
         allocate(temp_lon(lon_dim))
         allocate(temp_lon2(lon_dim))
         allocate(temp_var(lon_dim,lat_dim,time_dim)) 
         allocate(temp_var2(lon_dim,lat_dim)) 
         temp_lon = this%lon 
         n1 = lon_dim/2 +1 
         dx = temp_lon(2) - temp_lon(1)
         lon0 = temp_lon(n1) 
         do i = 1, lon_dim
            temp_lon2(i) = lon0 + dx * (i-1) - 360. 
         enddo
         this%lon(:) = temp_lon2(:) 
         temp_var = this%data1
         do nt = 1, time_dim
            temp_var2(1:n1-1,:) = temp_var(n1:lon_dim,:,nt) 
            temp_var2(n1:lon_dim,:) = temp_var(1:n1-1,:,nt) 
            this%data1(:,:,nt)= temp_var2(:,:)
         enddo
         deallocate(temp_lon)
         deallocate(temp_lon2)
         deallocate(temp_var) 
         deallocate(temp_var2) 
      endif

      if (.not. allocated(qc1%flux_orig)) allocate(qc1%flux_orig(time_dim))
      if (.not. allocated(qc1%flux_interp)) allocate(qc1%flux_interp(time_dim))
      do i = 1, lon_dim-1
         if ((this%lon(i) <= qc1%limit(3)) .and. (this%lon(i+1) >= qc1%limit(3))) iq1 = i
         if ((this%lon(i) <= qc1%limit(4)) .and. (this%lon(i+1) >= qc1%limit(4))) iq2 = i
      enddo
      do i = 1, lat_dim-1
         if ((this%lat(i) <= qc1%limit(1)) .and. (this%lat(i+1) >= qc1%limit(1))) jq1 = i
         if ((this%lat(i) <= qc1%limit(2)) .and. (this%lat(i+1) >= qc1%limit(2))) jq2 = i
      enddo 
      !make lon and lat 2d fields
      allocate(mask1(lon_dim,lat_dim))
      allocate(tmp1(lon_dim,lat_dim))
      allocate(tmp2(lon_dim,lat_dim))
      ! Compute the area of grids in sq. km 
      allocate(area_orig(lon_dim, lat_dim))
      call grid_area(this%lon, this%lat, area_orig) 
      nval = (iq2-iq1+1) * (jq2-jq1+1)
      mask1(iq1:iq2, jq1:jq2)=1
      sum_area_orig = sum(area_orig,mask=(mask1.eq.1))
      mask1 = 0
      nval=0
      do j = 1, lat_dim
         do i = 1, lon_dim
            if ((this%lon(i) >= qc1%limit(3)) .and. (this%lon(i)<= qc1%limit(4)) &
                 .and. (this%lat(j) >= qc1%limit(1)) .and. (this%lat(j) <= qc1%limit(2))) then 
               mask1(i,j) = 1
               nval = nval+1
            endif
         enddo
      enddo
      
      !Get the total emission = sum(area * emission)
      sum_area_orig = sum(area_orig,mask=(mask1.eq.1))
      do nt = 1, time_dim
         tmp1 = this%data1(:,:,nt) 
         tmp2= tmp1*area_orig
         qc1%flux_orig(nt) = sum(tmp2,mask=(mask1.eq.1))/sum_area_orig
      enddo
      
      allocate(tmp(wrf_grid1%nx, wrf_grid1%ny)) 
      tmp = 0
      nval = 0
      do j = 1, wrf_grid1%ny
         do i= 1, wrf_grid1%nx
            if ((wrf_grid1%lon(i,j) >= qc1%limit(3)-dx/2.) .and. &
                 (wrf_grid1%lon(i,j) <= qc1%limit(4)+dx/2.) .and. &
                 (wrf_grid1%lat(i,j) >= qc1%limit(1)-dy/2.) .and. &
                 (wrf_grid1%lat(i,j) <= qc1%limit(2)+dy/2.)) then 
               tmp(i,j)=1
               nval = nval +1 
            endif
         enddo
      enddo

      ! Compute the area of wrf grid in sq. km
      allocate(area_wrf(wrf_grid1%nx,wrf_grid1%ny))
      area_wrf = wrf_grid1%area 
      sum_area_wrf = sum(area_wrf, mask = (tmp .eq. 1))

      inverse_delta_lon = &
           1./(this%lon(2) - this%lon(1))
      inverse_delta_lat = &
           1./(this%lat(2) - this%lat(1))
      allocate(var_emission(lon_dim,lat_dim,1))
 
      write(6,'(1x,a,2(i4),a,2(i4))')'[FLUX_INTERP]    interp data domain dims: ',&
           lon_dim, lat_dim, ' to ', wrf_grid1%nx, wrf_grid1%ny
!! Horizontally interpolate these datasets to WRF domain
!! Now do it for outer domain only
      wrf_x_dim = wrf_grid1%nx
      wrf_y_dim = wrf_grid1%ny
      allocate(horiz_interp(wrf_x_dim,wrf_y_dim,1,1))
      if (.not. allocated(var_interp)) allocate(var_interp(wrf_x_dim,wrf_y_dim,time_dim))
      allocate(temp(wrf_x_dim,wrf_y_dim))
      do nn = 1, time_dim 
         var_emission(:,:,1) = this%data1(:,:,nn)
         call horiz_interpolation(wrf_x_dim, wrf_y_dim, &
              wrf_grid1%lat, wrf_grid1%lon,     &
              this%nx, this%ny, this%lon, this%lat,1, &
              var_emission, &
              inverse_delta_lon, inverse_delta_lat, &
              horiz_interp)
         var_interp(:,:,nn)=horiz_interp(:,:,1,1)
      enddo  ! end of months
!!$      print *, 'myfilename', myfilename
!!$      write(6,'(1x,a,4x,a,4x,a)')' nt ','flux_orig','flux_interp'
!!$      print *, ' ===================================================='
      do nt = 1, time_dim        
         qc1%flux_interp(nt) = sum((var_interp(:,:, nt)*area_wrf(:,:)),mask= (tmp .eq. 1))/sum_area_wrf
!!$         write(6,'(1x,i4,4x,f12.4,4x,f12.4)')&
!!$              nt, qc1%flux_orig(nt), qc1%flux_interp(nt)
      enddo

      deallocate(area_wrf)    
      deallocate(area_orig)    
      deallocate(var_emission)    
      deallocate (horiz_interp)
      deallocate (mask1)
      deallocate (tmp1)
      deallocate (tmp2)
      deallocate (tmp)
      deallocate (temp)

      deallocate(this%data_times)
      deallocate(this%lon)
      deallocate(this%lat)
      deallocate(this%data1)
      write(6,'(1x,a,a,a)')'[FLUX_INTERP] END SUB read_interp_myfile3d'
      
   end subroutine read_interp_myfile3d
   
!------------------------------------------------------------------------------
   subroutine  read_interp_myfile4d(this, vardir, var_interp,wrf_grid1,qc2)
!------------------------------------------------------------------------------
      type (fileDetails), intent(inout) :: this 
      type (nuwrf_grid), intent(in) :: wrf_grid1
      character(len=*), intent(in) :: vardir
      real, allocatable, intent(out) :: var_interp(:,:,:,:)
      type (qc), intent(inout) :: qc2
      integer :: lon_dim, lat_dim, time_dim
      integer :: icount,ncid, nn
      character(len=4) :: cyear
      character(len=2) :: cmonth,cday,chour,cminute,csecond
      real, allocatable :: horiz_interp(:,:,:,:), &
                           var_emission(:,:,:)
      real, allocatable :: temp_var(:,:,:)
      real :: inverse_delta_lon, inverse_delta_lat,lon0,dx,lat0,dy
      character(len=MAX_DATE_LEN):: wrf_datetime !, wrf_datetime2 
      character(len=MAX_DATE_LEN), allocatable ::file_times(:) 
      character(len = 30) :: varname, myfilename
      character(len = 132) :: tmp_file
      integer :: i,j, wrf_x_dim, wrf_y_dim,nt,ndays, nmonth,ydays,jj,myyear      
      integer :: iq1, iq2, jq1, jq2, iq, nval
      integer, allocatable :: tmp(:,:),mask1(:,:)
      real, allocatable :: area_orig(:,:), area_wrf(:,:),tmp1(:,:), tmp2(:,:)
      real :: sum_area_wrf, sum_area_orig
!    
!    
      write(6,'(1x,a,a,a)')'[FLUX_INTERP] SUB read_interp_myfile4d'
      ! Get the year 
      nmonth = 12
      ndays = 31

      if (.not. allocated(var_interp)) &
          allocate(var_interp(wrf_grid1%nx,wrf_grid1%ny,ndays,nmonth))
      wrf_datetime = wrf_grid1%wrf_datetimes(1)
      read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear,cmonth,cday,chour,cminute,csecond 
      myfilename = this%filename
      tmp_file = trim(vardir)//"/"//myfilename(1:5)//cyear//myfilename(10:30)
      write(6,'(1x,a,a,a)')'[FLUX_INTERP]    read <', trim(tmp_file),'>'
      ydays = 365
      read(wrf_datetime,"(i4)")myyear
 
      if(mod(myyear,400) == 0 .or. (mod(myyear,4) == 0 .and. (mod(myyear,100) .ne. 0))) then
           ydays = 366
      endif

      !dy = 0.5 ! Hard coded for emission data 0.5x0.5deg resolution
      !dx = 0.5
      !lat0 = -89.75

      if (myfilename(1:3) == 'npp' .or. myfilename(1:3) == 'NPP' ) then 
          varname = "NPP" 
      elseif (myfilename(1:3) == 'res' .or. myfilename(1:3) == 'RES' ) then 
          varname = "RespCO2" 
      elseif (myfilename(1:3) == 'bfu' .or. myfilename(1:3) == 'BFU' ) then 
          varname = "BioFuel" 
      elseif (myfilename(1:3) == 'ffu' .or. myfilename(1:3) == 'FFU' ) then 
          varname = "FossilFuel" 
      elseif (myfilename(1:3) == 'oco' .or. myfilename(1:3) == 'OCO' ) then 
          varname = "OceanCO2" 
      elseif (myfilename(1:3) == 'fir' .or. myfilename(1:3) == 'FIR' ) then 
          varname = "WildFire" 
      endif 


      ncid = open_netcdf_readfile(trim(tmp_file))
      time_dim = read_netcdf_dimension(ncid,"time")
      lon_dim = read_netcdf_dimension(ncid,'lon')
      lat_dim = read_netcdf_dimension(ncid,'lat') 
      allocate(this%data_times(ydays))
      allocate(this%lon(lon_dim))
      allocate(this%lat(lat_dim))
      allocate(this%data1(lon_dim,lat_dim,ydays))
      allocate(tmp1(lon_dim,lat_dim))
      allocate(tmp2(lon_dim,lat_dim))
      allocate(mask1(lon_dim,lat_dim))
      allocate(area_orig(lon_dim, lat_dim))
      call read_netcdf_real_array_1d(ncid,'lon',lon_dim,this%lon)
      call read_netcdf_real_array_1d(ncid,'lat',lat_dim,this%lat)
      dx=this%lon(2)-this%lon(1)
      dy=this%lat(2)-this%lat(1)
      lon0=this%lon(1)
      lat0=this%lat(1)
!!$      write(6,'(1x,a,1x,4(f12.4))')'[FLUX_INTERP]    dx, dy, lon0, lat0: ',&
!!$           dx, dy, lon0, lat0

!      this%filename = trim(myfilename)
      this%ndata = ydays
      this%varname = myfilename(1:4)
      this%year = cyear
      this%nx = lon_dim
      this%ny = lat_dim
      call close_netcdf_file(ncid)
!!
!! Remove this after fixing the data for lat start from lat0 
!!Jossy
!        do j = 1, lat_dim
!          this%lat(j) = lat0 + dy * (j-1) 
!        enddo 
!!
      do i = 1, lon_dim-1
         if ((this%lon(i) <= qc2%limit(3)) .and. (this%lon(i+1) >= qc2%limit(3))) iq1 = i
         if ((this%lon(i) <= qc2%limit(4)) .and. (this%lon(i+1) >= qc2%limit(4))) iq2 = i
      enddo
      do i = 1, lat_dim-1
         if ((this%lat(i) <= qc2%limit(1)) .and. (this%lat(i+1) >= qc2%limit(1))) jq1 = i
         if ((this%lat(i) <= qc2%limit(2)) .and. (this%lat(i+1) >= qc2%limit(2))) jq2 = i
      enddo

      !make lon and lat 2d fields
      do j = 1,lat_dim
        tmp1(:,j) = this%lon(:) 
      end do
      do i = 1,lon_dim
        tmp2(i,:) = this%lat(:) 
      end do 
      ! Compute the area of grids in sq. km
      call grid_area(this%lon, this%lat, area_orig) 
      !Get the area weighted total emission = sum(area * emission)
      nval = (iq2-iq1+1) * (jq2-jq1+1)
      mask1 = 0
      mask1(iq1:iq2, jq1:jq2)=1
      sum_area_orig = sum(area_orig,mask=(mask1.eq.1)) 

      if (.not. allocated(qc2%flux_orig)) allocate(qc2%flux_orig(ydays))
      if (.not. allocated(qc2%flux_interp)) allocate(qc2%flux_interp(ydays))
      allocate(tmp(wrf_grid1%nx, wrf_grid1%ny)) 
      tmp = 0
      nval = 0
      do j = 1, wrf_grid1%ny
         do i= 1, wrf_grid1%nx
            if  ((wrf_grid1%lon(i,j) >= qc2%limit(3)-dx/2.) .and. &
                 (wrf_grid1%lon(i,j) <= qc2%limit(4)+dx/2.) .and. &
                 (wrf_grid1%lat(i,j) >= qc2%limit(1)-dy/2.) .and. &
                 (wrf_grid1%lat(i,j) <= qc2%limit(2)+dy/2.)) then 
               tmp(i,j)=1 
               nval = nval +1 
            endif
         enddo
      enddo

      ! Compute the area of wrf grid in sq. km
      allocate(area_wrf(wrf_grid1%nx,wrf_grid1%ny))
      area_wrf = wrf_grid1%area  
      sum_area_wrf = sum(area_wrf, mask = (tmp .eq. 1))
      iq = 0
      icount = 0
      do jj = 1, nmonth ! for 12 months
         write(cmonth,"(I2.2)") jj
         tmp_file = trim(vardir)//'/'//FIREname(1:5)//cyear//FIREname(10:16)//cmonth//'.nc'
         ncid = open_netcdf_readfile(trim(tmp_file))
         time_dim = read_netcdf_dimension(ncid,"time")
         lon_dim = read_netcdf_dimension(ncid,'lon')
         lat_dim = read_netcdf_dimension(ncid,'lat')
         allocate(temp_var(lon_dim, lat_dim, time_dim))
         allocate(file_times(time_dim))

         call read_netcdf_character_array_1d(ncid,"Times",MAX_DATE_LEN, &
              time_dim, file_times)
         call read_netcdf_real_array_1d(ncid,'lon',lon_dim,this%lon)
         call read_netcdf_real_array_1d(ncid,'lat',lat_dim,this%lat)
         call read_netcdf_real_array_3d(ncid,trim(varname),lon_dim,lat_dim, &
              time_dim,temp_var)
         call close_netcdf_file(ncid) 
         do nt = 1, time_dim
            iq = iq +1
            tmp1 = temp_var(:,:,nt) 
            tmp2= tmp1*area_orig
            qc2%flux_orig(iq) = sum(tmp2,mask=(mask1.eq.1))/sum_area_orig
         enddo

         inverse_delta_lon = &
              1./(this%lon(2) - this%lon(1))
         inverse_delta_lat = &
              1./(this%lat(2) - this%lat(1))
         allocate(var_emission(lon_dim,lat_dim,1))
         
         if (jj==1) then
            write(6,'(1x,a,2(i4),a,2(i4))')'[FLUX_INTERP]    interp data domain dims: ',&
           lon_dim, lat_dim, ' to ', wrf_grid1%nx, wrf_grid1%ny
         end if
      
!! Horizontally interpolate these datasets to WRF domain
!! Now do it for outer domain only
         wrf_x_dim = wrf_grid1%nx
         wrf_y_dim = wrf_grid1%ny
         do nn = 1, time_dim
            icount = icount +1 
            this%data_times(icount)=file_times(nn)
            allocate(horiz_interp(wrf_x_dim,wrf_y_dim,1,1))
            var_emission(:,:,1) = temp_var(:,:,nn)
            call horiz_interpolation(wrf_x_dim, wrf_y_dim, &
                 wrf_grid1%lat, wrf_grid1%lon,     &
                 this%nx, this%ny, this%lon, this%lat,1, &
                 var_emission, &
                 inverse_delta_lon, inverse_delta_lat, &
                 horiz_interp)
            var_interp(:,:,nn,jj)=horiz_interp(:,:,1,1) 
            qc2%flux_interp(icount) = sum((horiz_interp(:,:,1,1)*area_wrf(:,:)), &
                 mask=(tmp.eq.1))/sum_area_wrf
            deallocate(horiz_interp)
         enddo  ! end of 1 datafile (for 1 day)
         deallocate(file_times)
         deallocate(temp_var)
         deallocate(var_emission)  
      enddo ! end number of months
      deallocate(this%data_times)
      deallocate(this%lon)
      deallocate(this%lat)
      deallocate(this%data1)
      write(6,'(1x,a,a,a)')'[FLUX_INTERP] END SUB read_interp_myfile4d'

   end subroutine read_interp_myfile4d

!---------------------------------------------------------------------------
!
! ROUTINE:  write_flux_netcdf_file
!
! DESCRIPTION:  Writes a 3D real array to a new netCDF file. 
!
!---------------------------------------------------------------------------
   subroutine write_allflux_netcdf_file3d(dim1,dim2,dim3,STIME,& 
        var3d, var3d_0,var3d_1,var3d_2, var3d_3,var3d_4, &
        varname, varunits, vardescription, &
        varname0, vardescription0,  &
        varname1, vardescription1,  &
        varname2, vardescription2, &
        varname3, vardescription3, &
        varname4, vardescription4, filename,wrf_param)

      ! Arguments
      integer,intent(in) :: dim1
      integer,intent(in) :: dim2
      integer,intent(in) :: dim3
      integer, parameter :: dim4=1
      integer, parameter :: field_type(1) = 104
      real, intent(in) :: wrf_param(5)
      real,intent(in) :: var3d(dim1,dim2,dim3)
      real,intent(in) :: var3d_0(dim1,dim2,dim3)
      real,intent(in) :: var3d_1(dim1,dim2,dim3)
      real,intent(in) :: var3d_2(dim1,dim2,dim3)
      real,intent(in) :: var3d_3(dim1,dim2,dim3)
      real,intent(in) :: var3d_4(dim1,dim2,dim3)


      real :: var4d(dim1,dim2,dim3,1)
      character(len=*),intent(in) :: varname
      character(len=*),intent(in) :: filename
      character(len=*),intent(in) :: varunits
      character(len=*),intent(in) :: vardescription
      character(len=*),intent(in) :: vardescription0
      character(len=*),intent(in) :: vardescription1
      character(len=*),intent(in) :: vardescription2
      character(len=*),intent(in) :: vardescription3
      character(len=*),intent(in) :: vardescription4
      character(len=*),intent(in) :: varname0
      character(len=*),intent(in) :: varname1
      character(len=*),intent(in) :: varname2
      character(len=*),intent(in) :: varname3
      character(len=*),intent(in) :: varname4

      character(len=10),parameter :: STIME_NAME = "Times"
      character(len=MAX_DATE_LEN), intent(in)::STIME(dim3)
      ! Local variables
      integer :: nc_id,var_id,var_id0,var_id1, var_id2, var_id3, var_id4
      integer :: dim1_id, dim2_id, dim3_id,Timevar_id,Time2varid,dim4_id
      character (len = *), parameter :: TIME_UNITS2 = "YYYY-MM-DD_HH:MI:SS"
      character (len=256) ::attrname
      integer :: attr_int(1)
      real :: attr_real(1) 

      ! Open the file
      nc_id = open_netcdf_newfile(trim(filename))
      !call set_netcdf_define_mode(nc_id)
      
      ! Define the dimensions
      dim1_id = write_netcdf_dimension(nc_id,"west_east",dim1)
      dim2_id = write_netcdf_dimension(nc_id,"south_north",dim2)
      Timevar_id = write_netcdf_dimension(nc_id,'DateStrLen',MAX_DATE_LEN)  
      dim4_id = write_netcdf_dimension(nc_id,"emissions_zdim_stag",dim4)
      dim3_id = write_netcdf_dimension(nc_id,"Time",dim3)
      ! Define the arrays
      Time2varid = check_define_netcdf_char_array_2d(nc_id,trim(STIME_NAME), &
           "DateStrLen","Time")

      var_id = check_define_netcdf_real_array_4d(nc_id,trim(varname), &
           "west_east","south_north","emissions_zdim_stag","Time")
      call write_netcdf_integer_attribute(nc_id,var_id,"FieldType", &
             1,field_type)
      call write_netcdf_text_attribute(nc_id,var_id,"MemoryOrder", &
           3,'XYZ')
      call write_netcdf_text_attribute(nc_id,var_id,"descrition", &
           len_trim(vardescription),trim(vardescription))   
      call write_netcdf_text_attribute(nc_id,var_id,"units", &
           len_trim(varunits),trim(varunits))   
      call write_netcdf_text_attribute(nc_id,var_id,"coordinates", &
           10,"XLONG XLAT")   
      call write_netcdf_text_attribute(nc_id,var_id,"stagger",1," ")


      var_id0 = check_define_netcdf_real_array_4d(nc_id,trim(varname0), &
           "west_east","south_north","emissions_zdim_stag","Time")
      call write_netcdf_integer_attribute(nc_id,var_id0,"FieldType", &
           1,field_type)
      call write_netcdf_text_attribute(nc_id,var_id0,"MemoryOrder", &
           3,'XYZ')
      call write_netcdf_text_attribute(nc_id,var_id0,"descrition", &
           len_trim(vardescription0),trim(vardescription0))   
      call write_netcdf_text_attribute(nc_id,var_id0,"units", &
           len_trim(varunits),trim(varunits))   
      call write_netcdf_text_attribute(nc_id,var_id0,"coordinates", &
           10,"XLONG XLAT")   
      call write_netcdf_text_attribute(nc_id,var_id0,"stagger",1," ")

      var_id1 = check_define_netcdf_real_array_4d(nc_id,trim(varname1), &
           "west_east","south_north","emissions_zdim_stag","Time")
      call write_netcdf_integer_attribute(nc_id,var_id1,"FieldType", 1,field_type)
      call write_netcdf_text_attribute(nc_id,var_id1,"MemoryOrder", 3,'XYZ')
      call write_netcdf_text_attribute(nc_id,var_id1,"description", &
           len_trim(vardescription1),trim(vardescription1))   
      call write_netcdf_text_attribute(nc_id,var_id1,"units", &
           len_trim(varunits),trim(varunits))   
      call write_netcdf_text_attribute(nc_id,var_id1,"coordinates", 10,"XLONG XLAT")   
      call write_netcdf_text_attribute(nc_id,var_id1,"stagger", 1," ")

      var_id2 = check_define_netcdf_real_array_4d(nc_id,trim(varname2), &
           "west_east","south_north","emissions_zdim_stag","Time")
      call write_netcdf_integer_attribute(nc_id,var_id2,"FieldType", 1,field_type)
      call write_netcdf_text_attribute(nc_id,var_id2,"MemoryOrder", 3,'XYZ')
      call write_netcdf_text_attribute(nc_id,var_id2,"descrition", &
           len_trim(vardescription2),trim(vardescription2))   
      call write_netcdf_text_attribute(nc_id,var_id2,"units", &
           len_trim(varunits),trim(varunits))   
      call write_netcdf_text_attribute(nc_id,var_id2,"coordinates", 10,"XLONG XLAT")   
      call write_netcdf_text_attribute(nc_id,var_id2,"stagger", 1," ")

      var_id3 = check_define_netcdf_real_array_4d(nc_id,trim(varname3), &
           "west_east","south_north","emissions_zdim_stag","Time")
      call write_netcdf_integer_attribute(nc_id,var_id3,"FieldType", 1,field_type)
      call write_netcdf_text_attribute(nc_id,var_id3,"MemoryOrder", 3,'XYZ')
      call write_netcdf_text_attribute(nc_id,var_id3,"descrition", &
           len_trim(vardescription3),trim(vardescription3))   
      call write_netcdf_text_attribute(nc_id,var_id3,"units", &
           len_trim(varunits),trim(varunits))   
      call write_netcdf_text_attribute(nc_id,var_id3,"coordinates", 10,"XLONG XLAT")   
      call write_netcdf_text_attribute(nc_id,var_id3,"stagger", 1," ")

      var_id4 = check_define_netcdf_real_array_4d(nc_id,trim(varname4), &
           "west_east","south_north","emissions_zdim_stag","Time")
      call write_netcdf_integer_attribute(nc_id,var_id4,"FieldType", 1,field_type)
      call write_netcdf_text_attribute(nc_id,var_id4,"MemoryOrder", 3,'XYZ')
      call write_netcdf_text_attribute(nc_id,var_id4,"descrition", &
           len_trim(vardescription4),trim(vardescription4))   
      call write_netcdf_text_attribute(nc_id,var_id4,"units", &
           len_trim(varunits),trim(varunits))   
      call write_netcdf_text_attribute(nc_id,var_id4,"coordinates", 10,"XLONG XLAT")   
      call write_netcdf_text_attribute(nc_id,var_id4,"stagger", 1," ")



!
! Write the global attributes same as in wrfinput
!
      call write_netcdf_text_global_attribute(nc_id,"MMINLU",4,"USGS") 
      attrname = 'WEST-EAST_GRID_DIMENSION'; attr_int(1) = dim2+1
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'NORTH-SOUTH_GRID_DIMENSION'; attr_int(1) = dim1+1
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'BOTTOM-TOP_GRID_DIMENSION'; attr_int(1) = 41
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'DX'; attr_real(1) = wrf_param(1)
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      attrname = 'DY'; attr_real(1) = wrf_param(2)
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      call write_netcdf_text_global_attribute(nc_id,"GRIDTYPE",1,"C") 
      attrname = 'DIFF_OPT'; attr_int(1) = 1
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'MAP_PROJ'; attr_int(1) = 1
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'NUM_LAND_CAT'; attr_int(1) = 24
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'GRID_ID'; attr_int(1) = 1
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'PARENT_ID'; attr_int(1) = 0
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'I_PARENT_START'; attr_int(1) = 1
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'J_PARENT_START'; attr_int(1) = 1
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'PARENT_GRID_RATIO'; attr_int(1) = 1
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'ISWATER'; attr_int(1) = 16
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'ISICE'; attr_int(1) = 24
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'ISOILWATER'; attr_int(1) = 14
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'ISLAKE'; attr_int(1) = -1
      call write_netcdf_integer_global_attribute(nc_id,attrname,1,attr_int) 
      attrname = 'CEN_LAT'; attr_real(1) = 43.
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      attrname = 'CEN_LON'; attr_real(1) = -98.
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      attrname = 'TRUELAT1'; attr_real(1) = 37.
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      attrname = 'TRUELAT2'; attr_real(1) = 49.
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      attrname = 'MOAD_CEN_LAT'; attr_real(1) = 43.
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      attrname = 'STAND_LON'; attr_real(1) = -98.
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      attrname = 'POLE_LAT'; attr_real(1) = 90.
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      attrname = 'POLE_LON'; attr_real(1) = 0.
      call write_netcdf_real_global_attribute(nc_id,attrname,1,attr_real) 
      ! Exit define mode
      call unset_netcdf_define_mode(nc_id)
      !Write times 
      call write_netcdf_character_array_1d(nc_id,Time2varid, & 
           MAX_DATE_LEN,dim3,STIME)
      ! Write the array
      var4d(:,:,:,1)=var3d(:,:,:)
!!$      print *, dim1, dim2, dim3, size(var3d,1), size(var3d,2), size(var3d,3)
      !call write_netcdf_real_array_3d(nc_id,var_id,dim1,dim2,dim3,var3d)
      call write_netcdf_real_array_4d(nc_id,var_id,dim1,dim2,dim3,1,var4d)
      var4d(:,:,:,1)=var3d_0(:,:,:)
!!$      print *, dim1, dim2, dim3, size(var3d_0,1), size(var3d_0,2), size(var3d_0,3)
      call write_netcdf_real_array_4d(nc_id,var_id0,dim1,dim2,dim3,1,var4d)
      var4d(:,:,:,1)=var3d_1(:,:,:)
      !call write_netcdf_real_array_3d(nc_id,var_id1,dim1,dim2,dim3,var3d_1)
      call write_netcdf_real_array_4d(nc_id,var_id1,dim1,dim2,dim3,1,var4d)
      var4d(:,:,:,1)=var3d_2(:,:,:)
      call write_netcdf_real_array_4d(nc_id,var_id2,dim1,dim2,dim3,1,var4d)
      var4d(:,:,:,1)=var3d_3(:,:,:)
      call write_netcdf_real_array_4d(nc_id,var_id3,dim1,dim2,dim3,1,var4d)
      var4d(:,:,:,1)=var3d_4(:,:,:)
      call write_netcdf_real_array_4d(nc_id,var_id4,dim1,dim2,dim3,1,var4d)
      ! enddo
      ! Close the file
      call close_netcdf_file(nc_id)

   end subroutine write_allflux_netcdf_file3d

!!---------------------------------------------------------------------------
   subroutine wrf_grid_area(lon,lat,area)
!!---------------------------------------------------------------------------
       integer :: nx1, ny1
       real, intent(in) :: lon(:,:), lat(:,:)
       real, intent(inout) :: area(:,:)
       real, allocatable :: ddx(:,:), ddy(:,:),ddx1(:,:)
       real :: deg2km, deg2rad
       integer :: i, j
       
       nx1 = size(lon,1)
       ny1 = size(lat, 2)
       allocate(ddx(nx1, ny1))
       allocate(ddy(nx1, ny1))
       allocate(ddx1(nx1, ny1))
       deg2km = 111.
       do j = 1, ny1
          do i = 1, nx1-1
             ddx(i,j)= deg2km*abs(lon(i+1,j)-lon(i,j))
          end do
       end do
       do i = 1, nx1
          do j = 1, ny1-1
             ddy(i,j)= deg2km*abs(lat(i,j+1)-lat(i,j))
          end do
       end do
       area = 0. 
       deg2rad = atan(1.) *4./180.  ! pi/180 
       do i = 1, nx1-1
          do j = 1, ny1-1 
             ddx1(i,j) =  ddx(i,j)* cos(deg2rad*(lat(i,j)+lat(i,j+1))/2.)
             area(i,j) = ddx1(i,j) * ddy(i,j) 
          end do
       end do
       deallocate (ddx)
       deallocate (ddx1)
       deallocate (ddy)
       
   end subroutine wrf_grid_area

!!---------------------------------------------------------------------------
   subroutine grid_area(lon,lat,area)
!!---------------------------------------------------------------------------
      integer :: nx1, ny1
      real, intent(in) :: lon(:), lat(:)
      real, intent(inout) :: area(:,:)
      real, allocatable :: ddx(:), ddy(:),ddx1(:)
      real :: deg2km, deg2rad
      integer :: i, j

      nx1 = size(lon,1)
      ny1 = size(lat, 1)
      allocate(ddx(nx1))
      allocate(ddy( ny1))
      allocate(ddx1(nx1))
      deg2km = 111.
      area = 0.0
      ddx = 0.0
      ddx1=0.0
      ddy = 0.0
      do i = 1, nx1-1
         ddx(i)= deg2km*abs(lon(i+1)-lon(i))
      enddo
      do j = 1, ny1-1
         ddy(j)= deg2km*abs(lat(j+1)-lat(j))
      enddo
      area = 0. 
      deg2rad = atan(1.) *4./180.  ! pi/180 
      do i = 1, nx1
         do j = 1, ny1 
            ddx1(i) =  ddx(i)* cos(deg2rad*(lat(j)))
            area(i,j) = ddx1(i) * ddy(j) 
         enddo
      enddo
      deallocate (ddx)
      deallocate (ddx1)
      deallocate (ddy)

   end subroutine grid_area
    
end module flux_interpolator_mod

