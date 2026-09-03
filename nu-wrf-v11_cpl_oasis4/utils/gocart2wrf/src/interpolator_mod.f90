!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration and Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  interpolator_mod
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Defines a data type and attendant methods for interpolating gocart data
! to wrfinput and wrfbdy files.
!------------------------------------------------------------------------------

module interpolator_mod

   ! Import modules
   use netcdf_util_mod
   use derived_variables_mod
   use sorted_datetimes_mod
   use interp_util_mod

   ! Reset defaults
   implicit none
   private

   ! Private data types to collect data specific to wrfinput and wrfbdy grids
   type wrfinput_grid
      private
      integer :: nc_id
      integer :: wrf_x_dim
      integer :: wrf_y_dim
      integer :: wrf_z_dim
      character(len=MAX_DATE_LEN) :: datetime
      real,allocatable :: wrf_latitudes(:,:)
      real,allocatable :: wrf_longitudes(:,:)
      real,allocatable :: wrf_pressures(:,:,:)
      real,allocatable :: gocart_horiz_interp_pressures(:,:,:) 
   end type wrfinput_grid
   type wrfbdy_grid
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
      real,allocatable :: wrf_latitudes_bxs(:,:)
      real,allocatable :: wrf_latitudes_bxe(:,:)
      real,allocatable :: wrf_latitudes_bys(:,:)
      real,allocatable :: wrf_latitudes_bye(:,:)
      real,allocatable :: wrf_longitudes_bxs(:,:)
      real,allocatable :: wrf_longitudes_bxe(:,:)
      real,allocatable :: wrf_longitudes_bys(:,:)
      real,allocatable :: wrf_longitudes_bye(:,:)
      real,allocatable :: wrf_pressures_bxs(:,:,:,:)
      real,allocatable :: wrf_pressures_bxe(:,:,:,:)
      real,allocatable :: wrf_pressures_bys(:,:,:,:)
      real,allocatable :: wrf_pressures_bye(:,:,:,:)
      real,allocatable :: gocart_horiz_interp_pressures_bxs(:,:,:,:)
      real,allocatable :: gocart_horiz_interp_pressures_bxe(:,:,:,:)
      real,allocatable :: gocart_horiz_interp_pressures_bys(:,:,:,:)
      real,allocatable :: gocart_horiz_interp_pressures_bye(:,:,:,:)
   end type wrfbdy_grid

   ! Define public data type for storing wrf grid data
   public :: interpolator
   type interpolator
      private
      integer :: num_wrf_domains
      type(wrfinput_grid), allocatable :: wrfinput(:)
      type(wrfbdy_grid) :: wrfbdy
      character(len=MAX_DIR_LEN) :: wrf_dir
      character(len=MAX_DIR_LEN) :: gocart_dir
      character(len=MAX_DIR_LEN) :: gocart_prefix
      character(len=MAX_SOURCE_LEN) :: gocart_source
      integer :: num_datetimes
      character(len=MAX_DATE_LEN),allocatable :: datetimes(:)
      logical,allocatable :: for_wrfbdy(:)
      logical,allocatable :: for_wrfinput(:,:)
      integer :: gocart_lon_dim
      integer :: gocart_lat_dim
      integer :: gocart_lev_dim
      real,allocatable :: gocart_longitudes(:)
      real,allocatable :: gocart_latitudes(:)
      integer :: num_aerosols
      character(len=80),allocatable :: aerosol_names(:)
      character(len=80),allocatable :: aerosol_units(:)
      character(len=80),allocatable :: aerosol_descriptions(:)
   end type interpolator

   ! Public methods for interpolator data type
   public :: create_interpolator
   public :: destroy_interpolator
   public :: run_interpolator

   ! Public constants
   character(len=*),parameter,public :: MERRA2 = 'MERRA2'
   character(len=*),parameter,public :: MERRAERO = 'MERRAERO'
   character(len=*),parameter,public :: GEOS5 = 'GEOS5'
   character(len=*),parameter,public :: OFFLINE = 'OFFLINE'

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  init_interpolator
   !
   ! DESCRIPTION:  Private method to initialize non-allocatable data members
   ! in interpolator type.
   !
   !---------------------------------------------------------------------------

   subroutine init_interpolator(this)

      ! Arguments
      type(interpolator), intent(inout) :: this

      ! Initialize immediate members
      this%num_wrf_domains = 0
      this%wrf_dir = trim('NULL')
      this%gocart_dir = trim('NULL')
      this%gocart_prefix = trim('NULL')
      this%gocart_source = trim('NULL')
      this%num_datetimes = 0
      this%gocart_lon_dim = 0
      this%gocart_lat_dim = 0
      this%gocart_lev_dim = 0
      this%num_aerosols = 0

      ! Initialize wrfbdy members
      this%wrfbdy%nc_id = 0
      this%wrfbdy%wrf_x_dim = 0
      this%wrfbdy%wrf_y_dim = 0
      this%wrfbdy%wrf_z_dim = 0
      this%wrfbdy%wrf_t_dim = 0
      this%wrfbdy%wrf_w_dim = 0
      this%wrfbdy%delta_time = 0

   end subroutine init_interpolator

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  create_interpolator
   !
   ! DESCRIPTION:  Public constructor method for interpolator data type.
   !
   !---------------------------------------------------------------------------

   function create_interpolator(wrf_dir, num_wrf_domains, gocart_dir, &
        gocart_prefix,num_aerosols,aerosol_names,aerosol_units, &
        aerosol_descriptions, gocart_source) result (this)

      ! Arguments
      ! Directory with wrfinput and wrfbdy files
      character(len=*),intent(in) :: wrf_dir 
      ! Number of wrf domains to process
      integer, intent(in) :: num_wrf_domains
      ! Directory with GOCART netCDF files
      character(len=*),intent(in) :: gocart_dir
      ! Prefix for GOCART netCDF filenames
      character(len=*),intent(in) :: gocart_prefix
      ! Number of aerosols to be interpolated
      integer,intent(in) :: num_aerosols
      ! Aerosol names
      character(len=*),intent(in) :: aerosol_names(num_aerosols)
      ! Aerosol units
      character(len=*),intent(in) :: aerosol_units(num_aerosols)
      ! Aerosol descriptions
      character(len=*),intent(in) :: aerosol_descriptions(num_aerosols)
      ! Source of GOCART data
      character(len=*),intent(in) :: gocart_source

      ! Return variable
      type(interpolator) :: this

      ! Local variables
      integer :: i

      ! Initialize the data type
      call init_interpolator(this)

      ! Assign the aerosol names to the data structure.
      ! TODO:  Put in private method.
      if (num_aerosols < 1) then
         print*,'ERROR, non-positive number of aerosols selected!'
         print*,'num_aerosols = ',num_aerosols
         stop
      end if
      this%num_aerosols = num_aerosols
      allocate(this%aerosol_names(num_aerosols))
      allocate(this%aerosol_units(num_aerosols))
      allocate(this%aerosol_descriptions(num_aerosols))
      do i = 1, num_aerosols
         this%aerosol_names(i) = trim(aerosol_names(i))
         this%aerosol_units(i) = trim(aerosol_units(i))
         this%aerosol_descriptions(i) = trim(aerosol_descriptions(i))
      end do

      ! Assemble space and time coordinate data from wrfinput and wrfbdy
      ! files.
      call create_interpolator_wrf(this,wrf_dir,num_wrf_domains)

      ! Assemble horizontally interpolated pressure values for gocart data
      call create_interpolator_gocart(this,gocart_dir,gocart_prefix, &
           gocart_source)

      return

   end function create_interpolator

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  create_interpolator_wrf
   !
   ! DESCRIPTION:  Private method used to assemble time and space coordinates
   ! from wrfinput and wrfbdy files.
   !
   !---------------------------------------------------------------------------

   subroutine create_interpolator_wrf(this,wrf_dir,num_wrf_domains)

      ! Arguments
      type(interpolator),intent(inout) :: this
      ! Directory with wrfinput and wrfbdy files
      character(len=*),intent(in) :: wrf_dir 
      ! Number of wrf domains to process
      integer, intent(in) :: num_wrf_domains

      ! Local variables
      real,allocatable :: phb(:,:,:,:)
      real,allocatable :: mub(:,:,:)
      real,allocatable :: rdnw(:,:)
      real,allocatable :: ph(:,:,:,:)
      real,allocatable :: t(:,:,:,:)
      real,allocatable :: qvapor(:,:,:,:)
      real,allocatable :: mu(:,:,:)

      character(len=MAX_DATE_LEN) :: datetimes(1)

      character(len=132) :: filename
      character(len=2) :: char_int
      integer :: ncid
      integer :: wrf_x_dim,wrf_y_dim,wrf_z_dim,wrf_t_dim,wrf_w_dim
      integer :: wrf_domain_id

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

         write(6,'(A,I2.2,A)') &
              '***Preparing wrfinput_d',wrf_domain_id,'***'

         ! Open the wrfinput netCDF file.  Since we will ultimately insert 
         ! aerosol data into it, open in write mode.
         write(char_int,"(I2.2)") wrf_domain_id
         filename = trim(wrf_dir)//"/wrfinput_d"//trim(char_int)
         ncid = open_netcdf_writefile(filename)
         this%wrfinput(wrf_domain_id)%nc_id = ncid

         ! Get required dimensions from file
         wrf_x_dim = read_netcdf_dimension(ncid,"west_east")
         wrf_y_dim = read_netcdf_dimension(ncid,"south_north")
         wrf_z_dim = read_netcdf_dimension(ncid,"bottom_top")
         this%wrfinput(wrf_domain_id)%wrf_x_dim = wrf_x_dim
         this%wrfinput(wrf_domain_id)%wrf_y_dim = wrf_y_dim
         this%wrfinput(wrf_domain_id)%wrf_z_dim = wrf_z_dim

         ! Copy latitude and longitude data from wrfinput
         call copy_latlon_from_wrfinput(this,ncid,wrf_domain_id)

         ! Read datetime from wrfinput
         call read_netcdf_character_array_1d(ncid,"Times",MAX_DATE_LEN,1, &
              datetimes)
         this%wrfinput(wrf_domain_id)%datetime = trim(datetimes(1))

         ! Read intermediate variables for calculating pressure
         allocate(ph(wrf_x_dim,wrf_y_dim,wrf_z_dim+1,1))   ! At cell interfaces
         allocate(t(wrf_x_dim,wrf_y_dim,wrf_z_dim,1))  ! At cell midpoints
         allocate(qvapor(wrf_x_dim,wrf_y_dim,wrf_z_dim,1)) ! At cell midpoints
         allocate(mu(wrf_x_dim,wrf_y_dim,1))        ! At base of column
         allocate(phb(wrf_x_dim,wrf_y_dim,wrf_z_dim+1,1))  ! At cell interfaces
         allocate(mub(wrf_x_dim,wrf_y_dim,1))       ! At base of column
         allocate(rdnw(wrf_z_dim,1))         ! Cell vertical thicknesses
         call read_netcdf_real_array_4d(ncid,"PH",wrf_x_dim,wrf_y_dim, &
              wrf_z_dim+1,1,ph)
         call read_netcdf_real_array_4d(ncid,"T",wrf_x_dim,wrf_y_dim, &
              wrf_z_dim,1,t)
         call read_netcdf_real_array_4d(ncid,"QVAPOR",wrf_x_dim,wrf_y_dim, &
              wrf_z_dim,1,qvapor)
         call read_netcdf_real_array_3d(ncid,"MU",wrf_x_dim,wrf_y_dim,1,mu)
         call read_netcdf_real_array_4d(ncid,"PHB",wrf_x_dim,wrf_y_dim, &
              wrf_z_dim+1,1,phb)
         call read_netcdf_real_array_3d(ncid,"MUB",wrf_x_dim,wrf_y_dim,1,mub)
         call read_netcdf_real_array_2d(ncid,"RDNW",wrf_z_dim,1,rdnw)

         ! Modify the wrfinput netCDF header to define the new aerosol 
         ! variables.  Doing so all at once reduces internal final copying
         ! and saves time.  Note that the variables will *not* be written at 
         ! this point.
         call update_wrfinput_header(this,ncid)

         ! Derive wrfinput pressure
         ! Note:  Allocates member array.
         call calc_wrfinput_pressure(this,wrf_domain_id,wrf_x_dim,wrf_y_dim, &
              wrf_z_dim,ph,phb,mu,mub,t,qvapor, rdnw)

         ! These variables aren't needed for wrfbdy, so we'll free up some 
         ! memory.
         deallocate(ph)   
         deallocate(t)
         deallocate(qvapor)
         deallocate(mu)

         ! We only have to process wrfbdy for domain 1, so if we're working
         ! with a different domain, we're done.  Deallocate the remaining 
         ! temporary arrays and move on to the next domain.
         if (wrf_domain_id .ne. 1) then
            deallocate(phb)  
            deallocate(mub)
            deallocate(rdnw)
            cycle
         end if

         print*,'***Preparing wrfbdy_d01***'

         ! Open the wrfbdy file
         filename = trim(wrf_dir)//"/wrfbdy_d"//trim(char_int)
         ncid = open_netcdf_writefile(filename)
         this%wrfbdy%nc_id = ncid

         ! Get the required dimensions from the file
         wrf_x_dim = read_netcdf_dimension(ncid,"west_east")
         wrf_y_dim = read_netcdf_dimension(ncid,"south_north")
         wrf_z_dim = read_netcdf_dimension(ncid,"bottom_top")
         wrf_t_dim = read_netcdf_dimension(ncid,"Time")
         wrf_w_dim = read_netcdf_dimension(ncid,"bdy_width")

         this%wrfbdy%wrf_x_dim = wrf_x_dim
         this%wrfbdy%wrf_y_dim = wrf_y_dim
         this%wrfbdy%wrf_z_dim = wrf_z_dim
         this%wrfbdy%wrf_t_dim = wrf_t_dim
         this%wrfbdy%wrf_w_dim = wrf_w_dim

         ! Sanity check the x, y, and z-dimensions
         call compare_wrfinput_wrfbdy_dims(this)

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
         call compare_wrfbdy_wrfinput_time(this)

         ! Calculate seconds between wrfbdy time levels
         this%wrfbdy%delta_time = &
              calc_delta_time_wrfbdy(this%wrfbdy%thisbdytime(1), &
              this%wrfbdy%nextbdytime(1))

         ! Handle pressure calculation for wrfbdy
         ! Note:  Allocates member arrays  
         call calc_wrfbdy_pressures_single_zone(this,phb,mub,rdnw,"BXS", &
              ncid)
         call calc_wrfbdy_pressures_single_zone(this,phb,mub,rdnw,"BXE", &
              ncid)
         call calc_wrfbdy_pressures_single_zone(this,phb,mub,rdnw,"BYS", &
              ncid)
         call calc_wrfbdy_pressures_single_zone(this,phb,mub,rdnw,"BYE", &
              ncid)

         ! Modify the wrfbdy netCDF header to define the new aerosol variables.
         ! Doing so all at once reduces internal final copying
         ! and saves time.  Note that the variables will *not* be written at 
         ! this point.
         call update_wrfbdy_header(this,ncid)

         ! Copy latitudes and longitudes for boundary zones
         ! Note:  Allocates member arrays.
         call copy_latlon_to_wrfbdy(this)

         ! Free memory
         deallocate(phb)  
         deallocate(mub)
         deallocate(rdnw)

      end do
   end subroutine create_interpolator_wrf

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  copy_latlon_from_wrfinput
   !
   ! DESCRIPTION:  Private method for copying latitude and longitude data
   ! from wrfinput file and into interpolator data type.
   !
   !---------------------------------------------------------------------------

   subroutine copy_latlon_from_wrfinput(this,ncid,wrf_domain_id)

      ! Arguments
      type(interpolator),intent(inout) :: this
      integer, intent(in) :: ncid
      integer, intent(in) :: wrf_domain_id

      ! Local variables
      real,allocatable :: xlat(:,:,:)
      real,allocatable :: xlong(:,:,:)
      integer :: wrf_x_dim, wrf_y_dim
      integer :: i,j

      ! Copy latitude and longitude data from wrfinput
      wrf_x_dim = this%wrfinput(wrf_domain_id)%wrf_x_dim
      wrf_y_dim = this%wrfinput(wrf_domain_id)%wrf_y_dim
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

   end subroutine copy_latlon_from_wrfinput

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calc_wrfinput_pressure
   !
   ! DESCRIPTION:  Private method to calculate pressure for wrfinput grid
   !               points.
   !
   !---------------------------------------------------------------------------

   subroutine calc_wrfinput_pressure(this,wrf_domain_id,wrf_x_dim,wrf_y_dim, &
        wrf_z_dim,ph,phb,mu,mub,t,qvapor,rdnw)

      ! Arguments
      type(interpolator),intent(inout) :: this
      integer,intent(in) :: wrf_domain_id
      integer,intent(in) :: wrf_x_dim
      integer,intent(in) :: wrf_y_dim
      integer,intent(in) :: wrf_z_dim
      real,intent(in) :: ph(wrf_x_dim,wrf_y_dim,wrf_z_dim+1,1)
      real,intent(in) :: phb(wrf_x_dim,wrf_y_dim,wrf_z_dim+1,1)
      real,intent(in) :: mu(wrf_x_dim,wrf_y_dim,1)
      real,intent(in) :: mub(wrf_x_dim,wrf_y_dim,1)
      real,intent(in) :: t(wrf_x_dim,wrf_y_dim,wrf_z_dim,1)
      real,intent(in) :: qvapor(wrf_x_dim,wrf_y_dim,wrf_z_dim,1)
      real,intent(in) :: rdnw(wrf_z_dim,1)

      real :: geopotential_above, geopotential_below
      real :: theta_moist, alpha_dry
      integer :: i,j,k

      allocate(this%wrfinput(wrf_domain_id)%wrf_pressures( &
           wrf_x_dim,wrf_y_dim,wrf_z_dim))
      do k = 1,wrf_z_dim
         do j = 1,wrf_y_dim
            do i = 1, wrf_x_dim
               geopotential_above = &
                    calc_geopotential(ph(i,j,k+1,1),phb(i,j,k+1,1))
               geopotential_below = &
                    calc_geopotential(ph(i,j,k,1),phb(i,j,k,1))
               alpha_dry = calc_inverse_dryair_density(geopotential_above, &
                    geopotential_below, rdnw(k,1), &
                    mu(i,j,1), mub(i,j,1))
               theta_moist = &
                    calc_moist_potential_temperature(t(i,j,k,1), &
                    qvapor(i,j,k,1))
               this%wrfinput(wrf_domain_id)%wrf_pressures(i,j,k) = &
                    calc_full_pressure(theta_moist,alpha_dry)
            end do
         end do
      end do

   end subroutine calc_wrfinput_pressure

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  compare_wrfbdy_wrfinput_time
   !
   ! DESCRIPTION:  Private method to ensure first time level in wrfbdy 
   !               matches the wrfinput time.
   !
   !---------------------------------------------------------------------------

   subroutine compare_wrfbdy_wrfinput_time(this)

      ! Arguments
      type(interpolator),intent(in) :: this

      ! Compare with first wrfbdy time level
      if (trim(this%wrfbdy%datetimes(1)) .ne. &
           trim(this%wrfinput(1)%datetime)) then
         print*,'ERROR, time stamps differ between wrfinput and wrfbdy!'
         print*,'wrfinput                    ', &
              trim(this%wrfinput(1)%datetime)
         print*,'wrfbdy (first time level):  ', &
              trim(this%wrfbdy%datetimes(1))
         stop
      end if

   end subroutine compare_wrfbdy_wrfinput_time

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  copy_latlon_to_wrfbdy
   !
   ! DESCRIPTION:  Private method for copying boundary zone latitudes and
   !               longitudes from wrfinput.
   !
   !---------------------------------------------------------------------------

   subroutine copy_latlon_to_wrfbdy(this)

      ! Arguments
      type(interpolator),intent(inout) :: this

      ! Local variables
      integer :: wrf_x_dim,wrf_y_dim,wrf_w_dim
      integer :: i,j,w

      ! Copy domain 1 wrfinput latitudes and longitudes to boundary zones
      wrf_x_dim = this%wrfbdy%wrf_x_dim
      wrf_y_dim = this%wrfbdy%wrf_y_dim
      wrf_w_dim = this%wrfbdy%wrf_w_dim
      allocate(this%wrfbdy%wrf_latitudes_bxs(wrf_y_dim,wrf_w_dim))
      allocate(this%wrfbdy%wrf_latitudes_bxe(wrf_y_dim,wrf_w_dim))
      allocate(this%wrfbdy%wrf_latitudes_bys(wrf_x_dim,wrf_w_dim))
      allocate(this%wrfbdy%wrf_latitudes_bye(wrf_x_dim,wrf_w_dim))
      allocate(this%wrfbdy%wrf_longitudes_bxs(wrf_y_dim,wrf_w_dim))
      allocate(this%wrfbdy%wrf_longitudes_bxe(wrf_y_dim,wrf_w_dim))
      allocate(this%wrfbdy%wrf_longitudes_bys(wrf_x_dim,wrf_w_dim))
      allocate(this%wrfbdy%wrf_longitudes_bye(wrf_x_dim,wrf_w_dim))
      do w = 1, wrf_w_dim
         do j = 1, wrf_y_dim
            i = w
            this%wrfbdy%wrf_latitudes_bxs(j,w) = &
                 this%wrfinput(1)%wrf_latitudes(i,j)
            this%wrfbdy%wrf_longitudes_bxs(j,w) = &
                 this%wrfinput(1)%wrf_longitudes(i,j)
         end do
      end do
      do w = 1, wrf_w_dim
         do j = 1, wrf_y_dim
            i = wrf_x_dim - w + 1
            this%wrfbdy%wrf_latitudes_bxe(j,w) = &
                 this%wrfinput(1)%wrf_latitudes(i,j)
            this%wrfbdy%wrf_longitudes_bxe(j,w) = &
                 this%wrfinput(1)%wrf_longitudes(i,j)
         end do
      end do
      do w = 1, wrf_w_dim
         do i = 1, wrf_x_dim
            j = w
            this%wrfbdy%wrf_latitudes_bys(i,w) = &
                 this%wrfinput(1)%wrf_latitudes(i,j)
            this%wrfbdy%wrf_longitudes_bys(i,w) = &
                 this%wrfinput(1)%wrf_longitudes(i,j)
         end do
      end do
      do w = 1, wrf_w_dim
         do i = 1, wrf_x_dim
            j = wrf_y_dim - w + 1
            this%wrfbdy%wrf_latitudes_bye(i,w) = &
                 this%wrfinput(1)%wrf_latitudes(i,j)
            this%wrfbdy%wrf_longitudes_bye(i,w) = &
                 this%wrfinput(1)%wrf_longitudes(i,j)
         end do
      end do

   end subroutine copy_latlon_to_wrfbdy

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  compare_wrfinput_wrfbdy_dims
   !
   ! DESCRIPTION:  Private method for comparing dimensions from wrfbdy and
   !               wrfinput files.
   !
   !---------------------------------------------------------------------------

   subroutine compare_wrfinput_wrfbdy_dims(this)

      ! Arguments
      type(interpolator), intent(in) :: this

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

   end subroutine compare_wrfinput_wrfbdy_dims

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  create_interpolator_gocart
   !
   ! DESCRIPTION: Private method for assembling horizontally interpolated
   ! GOCART pressure values on WRF grid.  Targets GEOS-5 netCDF4 files.
   !
   !---------------------------------------------------------------------------

   subroutine create_interpolator_gocart(this,gocart_dir,gocart_prefix, &
        gocart_source)

      ! Arguments
      type(interpolator),intent(inout) :: this
      ! Directory with GOCART netCDF files
      character(len=*),intent(in) :: gocart_dir
      ! Prefix for GOCART netCDF filenames
      character(len=*),intent(in) :: gocart_prefix
      ! Source of GOCART data
      character(len=*),intent(in) :: gocart_source

      ! Local variables
      type(sorted_datetimes) :: datetime_object
      integer :: ii,jj

      ! Local constants
      character(len=*), parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
      character(len=*), parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      this%gocart_dir = trim(gocart_dir)
      this%gocart_prefix = trim(gocart_prefix)

      ! Make sure gocart_source is uppercase
      this%gocart_source = trim(gocart_source)
      do jj = 1, len(this%gocart_source)
         ii = index(LOWER_CASE,this%gocart_source(jj:jj))
         if (ii .ne. 0) then
            this%gocart_source(jj:jj) = UPPER_CASE(ii:ii)
         end if
      end do

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

      print*,'***Summary of required datetimes***'
      do jj = 1, this%num_datetimes
         print*,'Datetime:  ',this%datetimes(jj)
         write(6,'(A,L1)') 'For wrfbdy_d01:  ',this%for_wrfbdy(jj)
         do ii = 1, this%num_wrf_domains
            write(6,'(A,I2.2,A,L1)') &
                 'For wrfinput_d',ii,':  ',this%for_wrfinput(ii,jj)
         end do
      end do

      ! Loop through each datetime and interpolate gocart data
      call interpolate_geos5_gocart_pressures(this)

   end subroutine create_interpolator_gocart

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  interpolate_geos5_gocart_pressures
   !
   ! DESCRIPTION:  Subroutine for reading and interpolating GEOS5 GOCART data
   ! pressures to WRF grid (wrfinput and/or wrfbdy as appropriate).
   !
   !---------------------------------------------------------------------------

   subroutine interpolate_geos5_gocart_pressures(this)

      ! Arguments
      type(interpolator),intent(inout) :: this

      ! Local variables
      character(len=256) :: gocart_filename
      integer :: gocart_lat_dim,gocart_lon_dim,gocart_lev_dim,gocart_time_dim
      double precision, allocatable :: gocart_longitudes_dp(:)
      double precision, allocatable :: gocart_latitudes_dp(:)
      real, allocatable :: gocart_longitudes(:)
      real, allocatable :: gocart_latitudes(:)
      real, allocatable :: gocart_ps(:,:,:) ! Surface pressure
      real, allocatable :: gocart_delp(:,:,:,:) ! Pressure thickness
      real, allocatable :: gocart_pressures(:,:,:,:) ! Total pressure 
      real, allocatable :: gocart_a(:) ! Hybrid A coefficient
      real, allocatable :: gocart_b(:) ! Hybrid B coefficient
      real :: gocart_p0 ! Reference pressure
      integer :: wrf_x_dim, wrf_y_dim, wrf_w_dim, wrf_t_dim
      integer :: ncid
      integer :: nn,dd
      integer :: time_level, ihour

      ! Preliminaries
      wrf_x_dim = this%wrfbdy%wrf_x_dim
      wrf_y_dim = this%wrfbdy%wrf_y_dim
      wrf_w_dim = this%wrfbdy%wrf_w_dim
      wrf_t_dim = this%wrfbdy%wrf_t_dim

      ! Loop through each datetime and interpolate gocart data
      do nn = 1, this%num_datetimes
         if (trim(this%gocart_source) .eq. MERRA2) then
            gocart_filename = build_merra2_filename(this%datetimes(nn), &
                 this%gocart_dir) 
            read(this%datetimes(nn),"(11X,I2,6X)") ihour
            time_level = (ihour/3) + 1
         else if (trim(this%gocart_source) .eq. OFFLINE) then
            ! Offline GOCART are split into collections with different
            ! aerosol species. But pressure data should be identical between
            ! files, so just pick one.
            gocart_filename = build_offline_filename(this%datetimes(nn), &
                 this%gocart_dir, this%gocart_prefix, 'bc1')
            read(this%datetimes(nn),"(11X,I2,6X)") ihour
            time_level = (ihour/3) + 1
         else if (trim(this%gocart_source) .eq. MERRAERO) then
            gocart_filename = build_merraero_filename(this%datetimes(nn), &
                 this%gocart_dir)
         else
            gocart_filename = build_gocart_filename(this%datetimes(nn), &
                 this%gocart_dir, &
                 this%gocart_prefix)
         end if

         print*,'Getting lat/lon/pressure from ',trim(gocart_filename)

         ncid = open_netcdf_readfile(gocart_filename)

         gocart_lon_dim = read_netcdf_dimension(ncid,'lon')
         gocart_lat_dim = read_netcdf_dimension(ncid,'lat')
         gocart_lev_dim = read_netcdf_dimension(ncid,'lev')
         gocart_time_dim = read_netcdf_dimension(ncid,'time')

         this%gocart_lon_dim = gocart_lon_dim
         this%gocart_lat_dim = gocart_lat_dim
         this%gocart_lev_dim = gocart_lev_dim

         ! Get longitudes
         if (trim(this%gocart_source) .eq. OFFLINE) then
            allocate(gocart_longitudes(gocart_lon_dim))
            call read_netcdf_real_array_1d(ncid,'lon',gocart_lon_dim, &
                 gocart_longitudes)
         else
            allocate(gocart_longitudes_dp(gocart_lon_dim))
            call read_netcdf_double_array_1d(ncid,'lon',gocart_lon_dim, &
                 gocart_longitudes_dp)
            allocate(gocart_longitudes(gocart_lon_dim))
            gocart_longitudes(:) = gocart_longitudes_dp(:)
            deallocate(gocart_longitudes_dp)
         end if
         if (.not. allocated(this%gocart_longitudes)) then
            allocate(this%gocart_longitudes(gocart_lon_dim))
            this%gocart_longitudes(:) = gocart_longitudes(:)
         end if

         ! Get latitudes
         if (trim(this%gocart_source) .eq. OFFLINE) then
            allocate(gocart_latitudes(gocart_lat_dim))
            call read_netcdf_real_array_1d(ncid,'lat',gocart_lat_dim, &
                 gocart_latitudes)
         else
            allocate(gocart_latitudes_dp(gocart_lat_dim))
            call read_netcdf_double_array_1d(ncid,'lat',gocart_lat_dim, &
                 gocart_latitudes_dp)
            allocate(gocart_latitudes(gocart_lat_dim))
            gocart_latitudes(:) = gocart_latitudes_dp(:)
            deallocate(gocart_latitudes_dp)
         end if
         if (.not. allocated(this%gocart_latitudes)) then
            allocate(this%gocart_latitudes(gocart_lat_dim))
            this%gocart_latitudes(:) = gocart_latitudes(:)
         end if

         ! Get pressure data
         if (trim(this%gocart_source) .eq. MERRA2) then
            ! Need single time slice
            allocate(gocart_ps(gocart_lon_dim,gocart_lat_dim,1))
            call read_netcdf_real_array_3d(ncid,'PS',gocart_lon_dim, &
                 gocart_lat_dim,1,gocart_ps, &
                 startarg=(/1, 1, time_level/), &
                 countarg=(/gocart_lon_dim,gocart_lat_dim,1/))
            allocate(gocart_delp(gocart_lon_dim,gocart_lat_dim, &
                 gocart_lev_dim,1))
            call read_netcdf_real_array_4d(ncid,'DELP',gocart_lon_dim, &
                 gocart_lat_dim,gocart_lev_dim,1,gocart_delp, &
                 startarg=(/1, 1, 1, time_level/), &
                 countarg=(/gocart_lon_dim,gocart_lat_dim,gocart_lev_dim,1/))
         else if (trim(this%gocart_source) .eq. OFFLINE) then
            ! Need single time slice
            allocate(gocart_ps(gocart_lon_dim,gocart_lat_dim,1))
            call read_netcdf_real_array_3d(ncid,'ps',gocart_lon_dim, &
                 gocart_lat_dim,1,gocart_ps, &
                 startarg=(/1, 1, time_level/), &
                 countarg=(/gocart_lon_dim,gocart_lat_dim,1/))
            gocart_ps(:,:,:) = gocart_ps*100. ! Convert to Pa
            call read_netcdf_real_variable(ncid,'p0',gocart_p0)
            allocate(gocart_a(gocart_lev_dim))
            gocart_p0 = gocart_p0*100. ! Convert to Pa
            call read_netcdf_real_array_1d(ncid,'hyam',gocart_lev_dim,gocart_a)
            allocate(gocart_b(gocart_lev_dim))
            call read_netcdf_real_array_1d(ncid,'hybm',gocart_lev_dim,gocart_b)
         else
            allocate(gocart_ps(gocart_lon_dim,gocart_lat_dim,gocart_time_dim))
            call read_netcdf_real_array_3d(ncid,'ps',gocart_lon_dim, &
                 gocart_lat_dim,gocart_time_dim,gocart_ps)
            allocate(gocart_delp(gocart_lon_dim,gocart_lat_dim, &
                 gocart_lev_dim,gocart_time_dim))
            call read_netcdf_real_array_4d(ncid,'delp',gocart_lon_dim, &
                 gocart_lat_dim,gocart_lev_dim,gocart_time_dim,gocart_delp)
         end if

         ! Close the GOCART netCDF file.
         call close_netcdf_file(ncid)

         ! Calculate total pressure at each GEOS mid-layer.  Note that 
         ! GEOS-5 grid convention is vertically flipped (k=1 is model top)
         ! compared to WRF (where k=1 is model terrain), so we'll flip it 
         ! back here. But offline GOCART doesn't have this issue.
         if (trim(this%gocart_source) .eq. MERRA2) then
            call flip_vertical_levels_4d(gocart_lon_dim,gocart_lat_dim,&
                 gocart_lev_dim,1,gocart_delp)
            allocate(gocart_pressures(gocart_lon_dim,gocart_lat_dim, &
                 gocart_lev_dim,1))
            call calc_pressure_geos5_gocart(gocart_lon_dim,gocart_lat_dim, &
                 gocart_lev_dim,1,gocart_ps,gocart_delp, &
                 gocart_pressures)
         else if (trim(this%gocart_source) .eq. OFFLINE) then
            allocate(gocart_pressures(gocart_lon_dim,gocart_lat_dim, &
                 gocart_lev_dim,1))
            call calc_pressure_offline_gocart(gocart_lon_dim,gocart_lat_dim, &
                 gocart_lev_dim,1,gocart_ps,gocart_p0,gocart_a,gocart_b, &
                 gocart_pressures)
         else
            call flip_vertical_levels_4d(gocart_lon_dim,gocart_lat_dim,&
                 gocart_lev_dim,gocart_time_dim,gocart_delp)
            allocate(gocart_pressures(gocart_lon_dim,gocart_lat_dim, &
                 gocart_lev_dim,gocart_time_dim))
            call calc_pressure_geos5_gocart(gocart_lon_dim,gocart_lat_dim, &
                 gocart_lev_dim,gocart_time_dim,gocart_ps,gocart_delp, &
                 gocart_pressures)
         end if
         if (allocated(gocart_ps)) deallocate(gocart_ps)
         if (allocated(gocart_delp)) deallocate(gocart_delp)
         if (allocated(gocart_a)) deallocate(gocart_a)
         if (allocated(gocart_b)) deallocate(gocart_b)

         ! Interpolate to wrfbdy boundary zones.
         if (this%for_wrfbdy(nn)) then

            ! Make sure these dims are for wrfbdy, and haven't been changed
            ! to wrfinput by code further down.
            wrf_x_dim = this%wrfbdy%wrf_x_dim
            wrf_y_dim = this%wrfbdy%wrf_y_dim
            wrf_w_dim = this%wrfbdy%wrf_w_dim
            wrf_t_dim = this%wrfbdy%wrf_t_dim

            ! Make sure interpolated gocart pressure arrays are allocated.
            ! (We can't do this until we read the number of GOCART levels 
            ! from the GEOS5 file, which we just did above).  Make sure
            ! we have an extra time level for calculating tendencies.
            if (.not. &
                 allocated(this%wrfbdy%gocart_horiz_interp_pressures_bxs)) then
               allocate(this%wrfbdy%gocart_horiz_interp_pressures_bxs( &
                    wrf_y_dim, gocart_lev_dim, wrf_w_dim, wrf_t_dim+1))
            end if
            if (.not. &
                 allocated(this%wrfbdy%gocart_horiz_interp_pressures_bxe)) then
               allocate(this%wrfbdy%gocart_horiz_interp_pressures_bxe( &
                    wrf_y_dim, gocart_lev_dim, wrf_w_dim, wrf_t_dim+1))
            end if
            if (.not. &
                 allocated(this%wrfbdy%gocart_horiz_interp_pressures_bys)) then
               allocate(this%wrfbdy%gocart_horiz_interp_pressures_bys( &
                    wrf_x_dim, gocart_lev_dim, wrf_w_dim, wrf_t_dim+1))
            end if
            if (.not. &
                 allocated(this%wrfbdy%gocart_horiz_interp_pressures_bye)) then
               allocate(this%wrfbdy%gocart_horiz_interp_pressures_bye( &
                    wrf_x_dim, gocart_lev_dim, wrf_w_dim, wrf_t_dim+1))
            end if

            ! Handle bxs region
            call horiz_interpolation_wrfbdy(wrf_y_dim, wrf_w_dim, &
                 this%wrfbdy%wrf_latitudes_bxs, &
                 this%wrfbdy%wrf_longitudes_bxs, &
                 gocart_lon_dim,gocart_lat_dim, &
                 gocart_longitudes, gocart_latitudes, &
                 gocart_lev_dim, gocart_pressures(1,1,1,1), &
                 this%wrfbdy%gocart_horiz_interp_pressures_bxs(1,1,1,nn))
            ! Handle bxe region
            call horiz_interpolation_wrfbdy(wrf_y_dim, wrf_w_dim, &
                 this%wrfbdy%wrf_latitudes_bxe, &
                 this%wrfbdy%wrf_longitudes_bxe, &
                 gocart_lon_dim,gocart_lat_dim, &
                 gocart_longitudes, gocart_latitudes, &
                 gocart_lev_dim, gocart_pressures(1,1,1,1), &
                 this%wrfbdy%gocart_horiz_interp_pressures_bxe(1,1,1,nn))
            ! Handle bys region
            call horiz_interpolation_wrfbdy(wrf_x_dim, wrf_w_dim, &
                 this%wrfbdy%wrf_latitudes_bys, &
                 this%wrfbdy%wrf_longitudes_bys, &
                 gocart_lon_dim,gocart_lat_dim, &
                 gocart_longitudes, gocart_latitudes, &
                 gocart_lev_dim, gocart_pressures(1,1,1,1), &
                 this%wrfbdy%gocart_horiz_interp_pressures_bys(1,1,1,nn))
            ! Handle bye region
            call horiz_interpolation_wrfbdy(wrf_x_dim, wrf_w_dim, &
                 this%wrfbdy%wrf_latitudes_bye, &
                 this%wrfbdy%wrf_longitudes_bye, &
                 gocart_lon_dim,gocart_lat_dim, &
                 gocart_longitudes, gocart_latitudes, &
                 gocart_lev_dim, gocart_pressures(1,1,1,1), &
                 this%wrfbdy%gocart_horiz_interp_pressures_bye(1,1,1,nn))

         end if ! If for_wrfbdy

         ! Interpolate to wrfinput grids
         do dd = 1, this%num_wrf_domains
            if (this%for_wrfinput(dd,nn)) then

               ! Update dimensions to reflect current wrfinput domain
               wrf_x_dim = this%wrfinput(dd)%wrf_x_dim
               wrf_y_dim = this%wrfinput(dd)%wrf_y_dim               
               allocate(this%wrfinput(dd)%gocart_horiz_interp_pressures( &
                    wrf_x_dim, wrf_y_dim, gocart_lev_dim))

               call horiz_interpolation(wrf_x_dim, wrf_y_dim, &
                    this%wrfinput(dd)%wrf_latitudes, &
                    this%wrfinput(dd)%wrf_longitudes, &
                    gocart_lon_dim,gocart_lat_dim, &
                    gocart_longitudes, gocart_latitudes, &
                    gocart_lev_dim, gocart_pressures(1,1,1,1), &
                    this%wrfinput(dd)%gocart_horiz_interp_pressures)

            end if ! If for_wrfinput
         end do ! Loop through wrfinput domains

         ! Clean up
         deallocate(gocart_latitudes)
         deallocate(gocart_longitudes)
         deallocate(gocart_pressures)

      end do ! Loop through times
   end subroutine interpolate_geos5_gocart_pressures

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  run_interpolator
   !
   ! DESCRIPTION:  Using the interpolator data type and information on
   ! variable names and conversion factors, reads/interpolates/writes
   ! an aerosol/chemistry species from GEOS5 GOCART to WRF.
   !
   !---------------------------------------------------------------------------

   subroutine run_interpolator(this,gocart_varname,wrf_varname, &
        conversion_factor)

      ! Arguments
      type(interpolator),intent(in) :: this
      character(len=*),intent(in) :: gocart_varname
      character(len=*),intent(in) :: wrf_varname
      real,intent(in) :: conversion_factor

      ! Local variables
      character(len=256) :: gocart_filename
      real, allocatable :: gocart_variables(:,:,:,:)
      real, allocatable :: gocart_variables_horiz_interp(:,:,:,:)
      real, allocatable :: gocart_variables_final_interp(:,:,:,:)
      integer :: wrf_x_dim, wrf_y_dim
      integer :: ncid
      integer :: ii,jj
      character(len=2) :: char_int
      integer :: ihour, time_level
      character(len=256) :: gocart_filename_prev
      
      ! Allocate temporary memory for GOCART data
      allocate(gocart_variables(this%gocart_lon_dim, &
           this%gocart_lat_dim, this%gocart_lev_dim, &
           this%num_datetimes))

      print*,'Reading ',trim(gocart_varname),' from GOCART files...'

      gocart_filename_prev = 'NULL'
      do ii = 1, this%num_datetimes

         if (trim(this%gocart_source) .eq. MERRA2 ) then
            gocart_filename = build_merra2_filename(this%datetimes(ii), &
                 this%gocart_dir) 
            read(this%datetimes(ii),"(11X,I2,6X)") ihour
            time_level = (ihour/3) + 1
         else if (trim(this%gocart_source) .eq. OFFLINE) then
            gocart_filename = build_offline_filename(this%datetimes(ii), &
                 this%gocart_dir,this%gocart_prefix,trim(gocart_varname))
            read(this%datetimes(ii),"(11X,I2,6X)") ihour
            time_level = (ihour/3) + 1
         else if (trim(this%gocart_source) .eq. MERRAERO) then
            gocart_filename = build_merraero_filename(this%datetimes(ii), &
                 this%gocart_dir) 
         else
            gocart_filename = build_gocart_filename(this%datetimes(ii), &
                 this%gocart_dir,this%gocart_prefix) 
         end if

         write(6,*)'Reading from ',trim(gocart_filename)

         if (ii .eq. 1) then
            ncid = open_netcdf_readfile(trim(gocart_filename))
            gocart_filename_prev = trim(gocart_filename)
         else if (ii .gt. 1 .and. &
              (trim(gocart_filename) .ne. trim(gocart_filename_prev))) then
            call close_netcdf_file(ncid)
            ncid = open_netcdf_readfile(trim(gocart_filename))
            gocart_filename_prev = trim(gocart_filename)
         end if

         if (trim(this%gocart_source) .eq. MERRA2 .or. &
              trim(this%gocart_source) .eq. OFFLINE) then
            ! Need single time slice
            call read_netcdf_real_array_4d(ncid,trim(gocart_varname), &
                 this%gocart_lon_dim,this%gocart_lat_dim,this%gocart_lev_dim, &
                 1,gocart_variables(1,1,1,ii), &
                 startarg=(/1,1,1,time_level/), &
                 countarg=(/this%gocart_lon_dim,this%gocart_lat_dim, &
                 this%gocart_lev_dim,1/))            
         else
            call read_netcdf_real_array_4d(ncid,trim(gocart_varname), &
                 this%gocart_lon_dim,this%gocart_lat_dim,this%gocart_lev_dim, &
                 1,gocart_variables(1,1,1,ii))
         end if

      end do

      ! Close the last GOCART netCDF file.
      call close_netcdf_file(ncid)

      ! Convert the units from GEOS-5 to WRF.
      gocart_variables(:,:,:,:) = conversion_factor*gocart_variables(:,:,:,:)

      ! Flip the data if necessary (GEOS-5 k=1 is model top, but WRF k=1 is 
      ! terrain)
      if (.not. trim(this%gocart_source) .eq. OFFLINE) then
         do ii = 1, this%num_datetimes
            call flip_vertical_levels_3d(this%gocart_lon_dim, &
                 this%gocart_lat_dim, this%gocart_lev_dim, &
                 gocart_variables(:,:,:,ii))
         enddo
      end if

      ! For debugging
#ifdef DEBUG
      call write_test_netcdf_file(this%gocart_lon_dim, this%gocart_lat_dim, &
           this%gocart_lev_dim,this%num_datetimes,gocart_variables, &
           trim(gocart_varname)//"_flipped", &
           trim(gocart_varname)//"_flipped.nc")
#endif

      ! Loop through wrfinput domains and interpolate as needed.
      do jj = 1, this%num_datetimes
         do ii = 1, this%num_wrf_domains

            if ( .not. this%for_wrfinput(ii,jj)) cycle
            wrf_x_dim = this%wrfinput(ii)%wrf_x_dim
            wrf_y_dim = this%wrfinput(ii)%wrf_y_dim
            allocate(gocart_variables_horiz_interp(wrf_x_dim,wrf_y_dim, &
                 this%gocart_lev_dim,1))

            call horiz_interpolation(wrf_x_dim, wrf_y_dim, &
                 this%wrfinput(ii)%wrf_latitudes, &
                 this%wrfinput(ii)%wrf_longitudes, &
                 this%gocart_lon_dim, this%gocart_lat_dim, &
                 this%gocart_longitudes, this%gocart_latitudes, &
                 this%gocart_lev_dim, gocart_variables(1,1,1,jj), &
                 gocart_variables_horiz_interp)

            write(char_int,"(I2.2)") ii

            ! For debugging
#ifdef DEBUG
            call write_test_netcdf_file(wrf_x_dim,wrf_y_dim, &
                 this%gocart_lev_dim,1,gocart_variables_horiz_interp, &
                 trim(gocart_varname)//"_horiz_interp_wrfinput_d"//char_int, &
                 trim(gocart_varname)//"_horiz_interp_wrfinput_d"//char_int//".nc")
#endif
            allocate(gocart_variables_final_interp(wrf_x_dim,wrf_y_dim, &
                 this%wrfinput(ii)%wrf_z_dim,1))

            call vert_interpolation(wrf_x_dim,wrf_y_dim, &
                 this%wrfinput(ii)%wrf_z_dim,1, &
                 this%gocart_lev_dim,this%wrfinput(ii)%wrf_pressures, &
                 this%wrfinput(ii)%gocart_horiz_interp_pressures, &
                 gocart_variables_horiz_interp,gocart_variables_final_interp)
            deallocate(gocart_variables_horiz_interp)

            ! For debugging
#ifdef DEBUG
            call write_test_netcdf_file(wrf_x_dim,wrf_y_dim, &
                 this%wrfinput(ii)%wrf_z_dim,1,gocart_variables_final_interp, &
                 trim(gocart_varname)//"_final_interp_wrfinput_d"//char_int, &
                 trim(gocart_varname)//"_final_interp_wrfinput_d"//char_int//".nc")
#endif
            ! Write to the wrfinput netCDF file
            write(6,'(A,A,A,I2.2)') &
                 'Writing ',trim(wrf_varname),' to wrfinput_d',ii
            call write_to_wrfinput(this,ii, &
                 gocart_variables_final_interp, trim(wrf_varname))

            ! Clean up
            deallocate(gocart_variables_final_interp)
         end do
      end do

      ! Interpolate to wrfbdy
      call run_interpolator_wrfbdy(this,gocart_variables,trim(wrf_varname))

      ! The end
      deallocate(gocart_variables)

   end subroutine run_interpolator

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  run_interpolator_wrfbdy
   !
   ! DESCRIPTION:  Private method for interpolating GOCART data to wrfbdy 
   ! zones.
   !
   !---------------------------------------------------------------------------

   subroutine run_interpolator_wrfbdy(this,gocart_variables_for_wrfbdy, &
        wrf_varname)

      ! Arguments
      type(interpolator),intent(in) :: this
      real,intent(in) :: gocart_variables_for_wrfbdy(this%gocart_lon_dim, &
           this%gocart_lat_dim,this%gocart_lev_dim,this%num_datetimes)
      character(len=*),intent(in) :: wrf_varname

      ! Local variables
      !real,allocatable :: gocart_variables_for_wrfbdy(:,:,:,:)
      integer :: num_for_wrfbdy
      integer :: wrfbdy_dim_1,wrfbdy_w_dim,wrfbdy_z_dim,wrfbdy_t_dim_plus_one
      integer :: nn,nn_for_wrfbdy

      ! We need to subset the GOCART data to those time levels actually 
      ! required by wrfbdy.  Use the same dimension order as for wrfinput.
      num_for_wrfbdy = 0
      do nn = 1, this%num_datetimes
         if (this%for_wrfbdy(nn)) then
            num_for_wrfbdy = num_for_wrfbdy + 1
         end if
      end do

      !allocate(gocart_variables_for_wrfbdy(this%gocart_lon_dim, &
      !     this%gocart_lat_dim, this%gocart_lev_dim, num_for_wrfbdy))
      nn_for_wrfbdy = 0
      do nn = 1, this%num_datetimes
         if (this%for_wrfbdy(nn)) then
            nn_for_wrfbdy = nn_for_wrfbdy + 1
            !do kk = 1, this%gocart_lev_dim
            !   do jj = 1, this%gocart_lat_dim
            !      do ii = 1, this%gocart_lon_dim
            ! gocart_variables_for_wrfbdy(ii,jj,kk,nn_for_wrfbdy) = &
            !     gocart_variables(ii,jj,kk,nn)
            !      end do
            !   end do
            !end do
         end if
      end do

      ! Now interpolate to bxs zone
      wrfbdy_dim_1 = this%wrfbdy%wrf_y_dim
      wrfbdy_z_dim = this%wrfbdy%wrf_z_dim
      wrfbdy_w_dim = this%wrfbdy%wrf_w_dim
      wrfbdy_t_dim_plus_one = num_for_wrfbdy
      call run_interpolator_wrfbdy_single_zone(this, &
           wrfbdy_dim_1,wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim_plus_one, &
           this%wrfbdy%wrf_latitudes_bxs,this%wrfbdy%wrf_longitudes_bxs, &
           this%wrfbdy%wrf_pressures_bxs, &
           this%wrfbdy%gocart_horiz_interp_pressures_bxs, &
           gocart_variables_for_wrfbdy,trim(wrf_varname),"BXS","BTXS")

      wrfbdy_dim_1 = this%wrfbdy%wrf_y_dim
      wrfbdy_z_dim = this%wrfbdy%wrf_z_dim
      wrfbdy_w_dim = this%wrfbdy%wrf_w_dim
      wrfbdy_t_dim_plus_one = num_for_wrfbdy
      call run_interpolator_wrfbdy_single_zone(this, &
           wrfbdy_dim_1,wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim_plus_one, &
           this%wrfbdy%wrf_latitudes_bxe,this%wrfbdy%wrf_longitudes_bxe, &
           this%wrfbdy%wrf_pressures_bxe, &
           this%wrfbdy%gocart_horiz_interp_pressures_bxe, &
           gocart_variables_for_wrfbdy,trim(wrf_varname),"BXE","BTXE")

      wrfbdy_dim_1 = this%wrfbdy%wrf_x_dim
      wrfbdy_z_dim = this%wrfbdy%wrf_z_dim
      wrfbdy_w_dim = this%wrfbdy%wrf_w_dim
      wrfbdy_t_dim_plus_one = num_for_wrfbdy
      call run_interpolator_wrfbdy_single_zone(this, &
           wrfbdy_dim_1,wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim_plus_one, &
           this%wrfbdy%wrf_latitudes_bys,this%wrfbdy%wrf_longitudes_bys, &
           this%wrfbdy%wrf_pressures_bys, &
           this%wrfbdy%gocart_horiz_interp_pressures_bys, &
           gocart_variables_for_wrfbdy,trim(wrf_varname),"BYS","BTYS")

      wrfbdy_dim_1 = this%wrfbdy%wrf_x_dim
      wrfbdy_z_dim = this%wrfbdy%wrf_z_dim
      wrfbdy_w_dim = this%wrfbdy%wrf_w_dim
      wrfbdy_t_dim_plus_one = num_for_wrfbdy
      call run_interpolator_wrfbdy_single_zone(this, &
           wrfbdy_dim_1,wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim_plus_one, &
           this%wrfbdy%wrf_latitudes_bye,this%wrfbdy%wrf_longitudes_bye, &
           this%wrfbdy%wrf_pressures_bye, &
           this%wrfbdy%gocart_horiz_interp_pressures_bye, &
           gocart_variables_for_wrfbdy,trim(wrf_varname),"BYE","BTYE")

      ! Clean up
      !deallocate(gocart_variables_for_wrfbdy)

   end subroutine run_interpolator_wrfbdy

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  run_interpolator_wrfbdy_single_zone
   !
   ! DESCRIPTION:  Private method for interpolating GOCART data to a 
   ! single wrfbdy zone.  This is a generic subroutine intended to be used
   ! by a higher-level driver.
   !
   !---------------------------------------------------------------------------

   subroutine run_interpolator_wrfbdy_single_zone(this, &
        wrfbdy_dim_1,wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim_plus_one, &
        wrf_latitudes_bdy,wrf_longitudes_bdy, wrf_pressures_bdy, &
        gocart_horiz_interp_pressures_bdy,gocart_variables_for_wrfbdy, &
        wrf_varname,bdy_suffix,bdy_t_suffix)

      ! Arguments
      type(interpolator),intent(in) :: this
      integer,intent(in) :: wrfbdy_dim_1
      integer,intent(in) :: wrfbdy_w_dim
      integer,intent(in) :: wrfbdy_z_dim
      integer,intent(in) :: wrfbdy_t_dim_plus_one
      real,intent(in) :: wrf_latitudes_bdy(wrfbdy_dim_1,wrfbdy_w_dim)
      real,intent(in) :: wrf_longitudes_bdy(wrfbdy_dim_1,wrfbdy_w_dim)
      real,intent(in) :: wrf_pressures_bdy(wrfbdy_dim_1, &
           wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim_plus_one)
      real,intent(in) :: gocart_horiz_interp_pressures_bdy(wrfbdy_dim_1, &
           wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim_plus_one)
      real,intent(in) :: gocart_variables_for_wrfbdy(this%gocart_lon_dim, &
           this%gocart_lat_dim,this%gocart_lev_dim,wrfbdy_t_dim_plus_one)
      character(len=*),intent(in) :: wrf_varname
      character(len=*),intent(in) :: bdy_suffix
      character(len=*),intent(in) :: bdy_t_suffix

      ! Local variables
      real,allocatable :: gocart_variables_horiz_interp(:,:,:,:)
      real,allocatable :: gocart_variables_final_interp(:,:,:,:)
      real,allocatable :: gocart_variable_tendencies(:,:,:,:)
      real :: inverse_wrfbdy_delta_time
      character(len=132) :: wrf_varname_bdy
      integer :: nn,kk,jj,ii

      inverse_wrfbdy_delta_time = 1. / this%wrfbdy%delta_time

      allocate(gocart_variables_horiz_interp(wrfbdy_dim_1, &
           this%gocart_lev_dim, &
           wrfbdy_w_dim,wrfbdy_t_dim_plus_one))

      do nn = 1, wrfbdy_t_dim_plus_one
         call horiz_interpolation_wrfbdy(wrfbdy_dim_1, &
              wrfbdy_w_dim, &
              wrf_latitudes_bdy, &
              wrf_longitudes_bdy, &
              this%gocart_lon_dim, this%gocart_lat_dim, &
              this%gocart_longitudes, this%gocart_latitudes, &
              this%gocart_lev_dim, gocart_variables_for_wrfbdy(1,1,1,nn), &
              gocart_variables_horiz_interp(1,1,1,nn))
      end do

      ! EMK For debugging
#ifdef DEBUG
      call write_test_netcdf_file(wrfbdy_dim_1,this%gocart_lev_dim, &
           wrfbdy_w_dim,wrfbdy_t_dim_plus_one-1, &
           gocart_variables_horiz_interp, &
           trim(wrf_varname)//"_"//trim(bdy_suffix)//"_horiz_interp_wrfbdy_d01", &
           trim(wrf_varname)//"_"//trim(bdy_suffix)//"_horiz_interp_wrfbdy_d01.nc")
#endif

      allocate(gocart_variables_final_interp(wrfbdy_dim_1, &
           wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim_plus_one))

      call vert_interpolation_wrfbdy(wrfbdy_dim_1, &
           wrfbdy_z_dim, &
           wrfbdy_w_dim, wrfbdy_t_dim_plus_one, &
           this%gocart_lev_dim,wrf_pressures_bdy, &
           gocart_horiz_interp_pressures_bdy, &
           gocart_variables_horiz_interp,gocart_variables_final_interp)
      deallocate(gocart_variables_horiz_interp)

      ! EMK For debugging
#ifdef DEBUG
      call write_test_netcdf_file(wrfbdy_dim_1,wrfbdy_z_dim, &
           wrfbdy_w_dim,wrfbdy_t_dim_plus_one-1, &
           gocart_variables_final_interp, &
           trim(wrf_varname)//"_"//trim(bdy_suffix)//"_final_interp_wrfbdy_d01", &
           trim(wrf_varname)//"_"//trim(bdy_suffix)//"_final_interp_wrfbdy_d01.nc")
#endif

      ! Write to wrfbdy
      wrf_varname_bdy=trim(wrf_varname)//"_"//trim(bdy_suffix)
      write(6,'(A,A,A)') 'Writing ',trim(wrf_varname_bdy), ' to wrfbdy_d01'
      call write_to_wrfbdy(this,wrfbdy_dim_1,wrfbdy_z_dim, &
           wrfbdy_w_dim, wrfbdy_t_dim_plus_one-1, &
           gocart_variables_final_interp, &
           trim(wrf_varname_bdy))

      ! Calculate tendencies.  Note that the number of tendencies is one less
      ! than the number of GOCART time levels, because each tendency is the
      ! difference between two time levels.
      allocate(gocart_variable_tendencies(wrfbdy_dim_1, wrfbdy_z_dim, &
           wrfbdy_w_dim,wrfbdy_t_dim_plus_one-1))
      do nn = 1, wrfbdy_t_dim_plus_one-1
         do kk = 1, wrfbdy_w_dim
            do jj = 1, wrfbdy_z_dim
               do ii = 1, wrfbdy_dim_1
                  gocart_variable_tendencies(ii,jj,kk,nn) = &
                       inverse_wrfbdy_delta_time * &
                       ( gocart_variables_final_interp(ii,jj,kk,nn+1) - &
                       gocart_variables_final_interp(ii,jj,kk,nn  ) ) 
               end do
            end do
         end do
      end do

      ! EMK For debugging
#ifdef DEBUG
      call write_test_netcdf_file(wrfbdy_dim_1,wrfbdy_z_dim, &
           wrfbdy_w_dim,wrfbdy_t_dim_plus_one-1, &
           gocart_variable_tendencies, &
           trim(wrf_varname)//"_"//trim(bdy_t_suffix)//"_final_interp_wrfbdy_d01", &
           trim(wrf_varname)//"_"//trim(bdy_t_suffix)//"_final_interp_wrfbdy_d01.nc")
#endif

      ! Write to wrfbdy
      wrf_varname_bdy=trim(wrf_varname)//"_"//trim(bdy_t_suffix)
      write(6,'(A,A,A)') 'Writing ',trim(wrf_varname_bdy),' to wrfbdy_d01'
      call write_to_wrfbdy(this,wrfbdy_dim_1,wrfbdy_z_dim,wrfbdy_w_dim, &
           wrfbdy_t_dim_plus_one-1,gocart_variable_tendencies, &
           trim(wrf_varname_bdy))

      ! Clean up
      deallocate(gocart_variable_tendencies)
      deallocate(gocart_variables_final_interp)

   end subroutine run_interpolator_wrfbdy_single_zone

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calc_wrfbdy_pressures_single_zone
   !
   ! DESCRIPTION:  A generic subroutine for calculating WRF-based pressures
   ! in a user-requested wrfbdy boundary zone.
   !
   !---------------------------------------------------------------------------

   subroutine calc_wrfbdy_pressures_single_zone(this,phb,mub,rdnw,bdy_zone, &
        ncid)

      ! Arguments
      type(interpolator),intent(inout) :: this
      real,intent(in) :: phb(this%wrfinput(1)%wrf_x_dim, &
           this%wrfinput(1)%wrf_y_dim, &
           this%wrfinput(1)%wrf_z_dim+1, &
           1)
      real,intent(in) :: mub(this%wrfinput(1)%wrf_x_dim, &
           this%wrfinput(1)%wrf_y_dim, &
           1)
      real,intent(in) :: rdnw(this%wrfinput(1)%wrf_z_dim,1)
      character(len=3),intent(in) :: bdy_zone
      integer,intent(in) :: ncid

      ! Local variables
      real,allocatable :: phb_bdy(:,:,:)
      real,allocatable :: mub_bdy(:,:)
      real,allocatable :: mu_bdy(:,:,:), mu_bdy_t(:,:,:)
      real,allocatable :: ph_bdy(:,:,:,:), ph_bdy_t(:,:,:,:)
      real,allocatable :: t_bdy(:,:,:,:), t_bdy_t(:,:,:,:)
      real,allocatable :: qvapor_bdy(:,:,:,:), qvapor_bdy_t(:,:,:,:)
      character(len=4) :: bdy_zone_tendency
      integer :: wrf_x_dim,wrf_y_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim
      integer :: bdy_first_dim
      integer :: ii,jj,kk,ww
      real :: inverse_delta_time

      inverse_delta_time = 1. / this%wrfbdy%delta_time 

      wrf_x_dim = this%wrfbdy%wrf_x_dim
      wrf_y_dim = this%wrfbdy%wrf_y_dim
      wrf_z_dim = this%wrfbdy%wrf_z_dim
      wrf_w_dim = this%wrfbdy%wrf_w_dim
      wrf_t_dim = this%wrfbdy%wrf_t_dim

      ! General logic applied to a boundary zone:
      ! 1.  Copy phb from wrfinput array to wrfbdy array
      ! 2.  Copy mub from wrfinput array to wrfbdy array
      ! 3.  Copy mu from wrfbdy file
      ! 4.  Copy coupled ph from wrfbdy and uncouple (divide by mu + mub)
      ! 5.  Copy coupled t from wrfdy and uncouple (divide by mu + mub)
      ! 6.  Copy coupled qvapor from wrfbdy and uncouple (divide by mu + mub)
      ! 7.  Calculate total pressure

      ! Steps 1 and 2 require special code for each boundary zone.  Also,
      ! it is convenient to copy the "length" dimension for a given boundary
      ! zone to a new variable that will be the first array dimension.

      if (bdy_zone == "BXS") then

         bdy_zone_tendency = "BTXS"

         ! Special handling for bxs boundary zone
         allocate(phb_bdy(wrf_y_dim,wrf_z_dim+1,wrf_w_dim))
         do ww = 1, wrf_w_dim
            do kk = 1, wrf_z_dim+1
               do jj = 1, wrf_y_dim
                  ii = ww
                  phb_bdy(jj,kk,ww) = phb(ii,jj,kk,1)
               end do
            end do
         end do
         allocate(mub_bdy(wrf_y_dim,wrf_w_dim))
         do ww = 1, wrf_w_dim
            do jj = 1, wrf_y_dim
               ii = ww
               mub_bdy(jj,ww) = mub(ii,jj,1)
            end do
         end do
         bdy_first_dim = wrf_y_dim

      else if (bdy_zone == "BXE") then

         bdy_zone_tendency = "BTXE"

         ! Special handling for bxe boundary zone
         allocate(phb_bdy(wrf_y_dim,wrf_z_dim+1,wrf_w_dim))
         do ww = 1, wrf_w_dim
            do kk = 1, wrf_z_dim+1
               do jj = 1, wrf_y_dim
                  ii = wrf_x_dim - ww + 1
                  phb_bdy(jj,kk,ww) = phb(ii,jj,kk,1)
               end do
            end do
         end do
         allocate(mub_bdy(wrf_y_dim,wrf_w_dim))
         do ww = 1, wrf_w_dim
            do jj = 1, wrf_y_dim
               ii = wrf_x_dim - ww + 1
               mub_bdy(jj,ww) = mub(ii,jj,1)
            end do
         end do
         bdy_first_dim = wrf_y_dim

      else if (bdy_zone == "BYS") then

         bdy_zone_tendency = "BTYS" 

         ! Special handling for bys boundary zone
         allocate(phb_bdy(wrf_x_dim,wrf_z_dim+1,wrf_w_dim))
         do ww = 1, wrf_w_dim
            do kk = 1, wrf_z_dim+1
               do ii = 1, wrf_x_dim
                  jj = ww
                  phb_bdy(ii,kk,ww) = phb(ii,jj,kk,1)
               end do
            end do
         end do
         allocate(mub_bdy(wrf_x_dim,wrf_w_dim))
         do ww = 1, wrf_w_dim
            do ii = 1, wrf_x_dim
               jj = ww
               mub_bdy(ii,ww) = mub(ii,jj,1)
            end do
         end do
         bdy_first_dim = wrf_x_dim

      else if (bdy_zone == "BYE") then

         bdy_zone_tendency = "BTYE"

         ! Special handling for bye boundary zone
         allocate(phb_bdy(wrf_x_dim,wrf_z_dim+1,wrf_w_dim))
         do ww = 1, wrf_w_dim
            do kk = 1, wrf_z_dim+1
               do ii = 1, wrf_x_dim
                  jj = wrf_y_dim - ww + 1
                  phb_bdy(ii,kk,ww) = phb(ii,jj,kk,1)
               end do
            end do
         end do
         allocate(mub_bdy(wrf_x_dim,wrf_w_dim))
         do ww = 1, wrf_w_dim
            do ii = 1, wrf_x_dim
               jj = wrf_y_dim - ww + 1
               mub_bdy(ii,ww) = mub(ii,jj,1)
            end do
         end do
         bdy_first_dim = wrf_x_dim

      else

         ! Internal error
         print*,'ERROR, illegal bdy_zone value entered!'
         print*,'Expected bxs, bxe, bys, or bye'
         print*,'Received ',bdy_zone
         stop

      end if

      ! Step 3...Generic code.
      allocate(mu_bdy(bdy_first_dim,wrf_w_dim,wrf_t_dim+1))
      call read_netcdf_real_array_3d(ncid,"MU_"//bdy_zone,bdy_first_dim, &
           wrf_w_dim,wrf_t_dim,mu_bdy)
      allocate(mu_bdy_t(bdy_first_dim,wrf_w_dim,wrf_t_dim))
      call read_netcdf_real_array_3d(ncid,"MU_"//bdy_zone_tendency, &
           bdy_first_dim, wrf_w_dim,wrf_t_dim,mu_bdy_t)
      do ww = 1, wrf_w_dim
         do ii = 1, bdy_first_dim
            mu_bdy(ii,ww,wrf_t_dim+1) = &
                 inverse_delta_time*mu_bdy_t(ii,ww,wrf_t_dim) + &
                 mu_bdy(ii,ww,wrf_t_dim)
         end do
      end do
      deallocate(mu_bdy_t)

      ! Step 4...Generic code.
      allocate(ph_bdy(bdy_first_dim,wrf_z_dim+1,wrf_w_dim,wrf_t_dim+1))
      call read_netcdf_real_array_4d(ncid,"PH_"//bdy_zone,bdy_first_dim, &
           wrf_z_dim+1, wrf_w_dim, wrf_t_dim, ph_bdy)
      allocate(ph_bdy_t(bdy_first_dim,wrf_z_dim+1,wrf_w_dim,wrf_t_dim))
      call read_netcdf_real_array_4d(ncid,"PH_"//bdy_zone_tendency, &
           bdy_first_dim, wrf_z_dim+1, wrf_w_dim, wrf_t_dim, ph_bdy_t)
      do ww = 1, wrf_w_dim
         do kk = 1, wrf_z_dim+1
            do ii = 1, bdy_first_dim
               ph_bdy(ii,kk,ww,wrf_t_dim+1) = &
                    inverse_delta_time*ph_bdy_t(ii,kk,ww,wrf_t_dim) + &
                    ph_bdy(ii,kk,ww,wrf_t_dim)
            end do
         end do
      end do
      deallocate(ph_bdy_t)
      call decouple_wrfbdy(bdy_first_dim, wrf_z_dim+1, wrf_w_dim, &
           wrf_t_dim+1, mu_bdy, mub_bdy, ph_bdy)

      ! Step 5...Generic code.
      allocate(t_bdy(bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim+1))
      call read_netcdf_real_array_4d(ncid,"T_"//bdy_zone,bdy_first_dim, &
           wrf_z_dim, wrf_w_dim,wrf_t_dim,t_bdy)
      allocate(t_bdy_t(bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim))
      call read_netcdf_real_array_4d(ncid,"T_"//bdy_zone_tendency, &
           bdy_first_dim,wrf_z_dim, wrf_w_dim,wrf_t_dim,t_bdy_t)
      do ww = 1, wrf_w_dim
         do kk = 1, wrf_z_dim
            do ii = 1, bdy_first_dim
               t_bdy(ii,kk,ww,wrf_t_dim+1) = &
                    inverse_delta_time*t_bdy_t(ii,kk,ww,wrf_t_dim) + &
                    t_bdy(ii,kk,ww,wrf_t_dim)
            end do
         end do
      end do
      deallocate(t_bdy_t)
      call decouple_wrfbdy(bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim+1, &
           mu_bdy, mub_bdy,t_bdy)

      ! Step 6...Generic code.
      allocate(qvapor_bdy(bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim+1))
      call read_netcdf_real_array_4d(ncid,"QVAPOR_"//bdy_zone,bdy_first_dim, &
           wrf_z_dim, wrf_w_dim,wrf_t_dim,qvapor_bdy)
      allocate(qvapor_bdy_t(bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim))
      call read_netcdf_real_array_4d(ncid,"QVAPOR_"//bdy_zone_tendency, &
           bdy_first_dim,wrf_z_dim, wrf_w_dim,wrf_t_dim,qvapor_bdy_t)
      do ww = 1, wrf_w_dim
         do kk = 1, wrf_z_dim
            do ii = 1, bdy_first_dim
               qvapor_bdy(ii,kk,ww,wrf_t_dim+1) = &
                    inverse_delta_time*qvapor_bdy_t(ii,kk,ww,wrf_t_dim) + &
                    qvapor_bdy(ii,kk,ww,wrf_t_dim)
            end do
         end do
      end do
      deallocate(qvapor_bdy_t)
      call decouple_wrfbdy(bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim+1, &
           mu_bdy, mub_bdy,qvapor_bdy)

      ! Step 7 requires special handling again, since the pressures will
      ! be stored to separate arrays (one for each boundary zone).
      if (bdy_zone == "BXS") then
         allocate(this%wrfbdy%wrf_pressures_bxs( &
              bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim+1))
         call calc_wrfbdy_pressure_internal(bdy_first_dim, wrf_z_dim, &
              wrf_w_dim, wrf_t_dim, ph_bdy, phb_bdy, rdnw, &
              mu_bdy, mub_bdy, t_bdy, qvapor_bdy, &
              this%wrfbdy%wrf_pressures_bxs)
      else if (bdy_zone == "BXE") then
         allocate(this%wrfbdy%wrf_pressures_bxe( &
              bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim+1))
         call calc_wrfbdy_pressure_internal(bdy_first_dim, wrf_z_dim, &
              wrf_w_dim, wrf_t_dim, ph_bdy, phb_bdy, rdnw, &
              mu_bdy, mub_bdy, t_bdy, qvapor_bdy, &
              this%wrfbdy%wrf_pressures_bxe)
      else if (bdy_zone == "BYS") then
         allocate(this%wrfbdy%wrf_pressures_bys( &
              bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim+1))
         call calc_wrfbdy_pressure_internal(bdy_first_dim, wrf_z_dim, &
              wrf_w_dim, wrf_t_dim, ph_bdy, phb_bdy, rdnw, &
              mu_bdy, mub_bdy, t_bdy, qvapor_bdy, &
              this%wrfbdy%wrf_pressures_bys)
      else
         allocate(this%wrfbdy%wrf_pressures_bye( &
              bdy_first_dim,wrf_z_dim,wrf_w_dim,wrf_t_dim+1))
         call calc_wrfbdy_pressure_internal(bdy_first_dim, wrf_z_dim, &
              wrf_w_dim, wrf_t_dim, ph_bdy, phb_bdy, rdnw, &
              mu_bdy, mub_bdy, t_bdy, qvapor_bdy, &
              this%wrfbdy%wrf_pressures_bye)
      end if

      ! Clean up
      deallocate(phb_bdy)
      deallocate(ph_bdy)
      deallocate(t_bdy)
      deallocate(mu_bdy)
      deallocate(mub_bdy)
      deallocate(qvapor_bdy)

   contains

      ! An internal generic subroutine to do the final pressure calculation.
      subroutine calc_wrfbdy_pressure_internal(bdy_first_dim, wrf_z_dim, &
           wrf_w_dim, wrf_t_dim, ph_bdy, phb_bdy, rdnw, &
           mu_bdy, mub_bdy, t_bdy, qvapor_bdy, wrf_pressures_bdy)

         ! Arguments
         integer,intent(in) :: bdy_first_dim
         integer,intent(in) :: wrf_z_dim
         integer,intent(in) :: wrf_w_dim
         integer,intent(in) :: wrf_t_dim
         real,intent(in) :: ph_bdy(bdy_first_dim,wrf_z_dim+1,wrf_w_dim, &
              wrf_t_dim+1)
         real,intent(in) :: phb_bdy(bdy_first_dim,wrf_z_dim+1,wrf_w_dim)
         real,intent(in) :: rdnw(wrf_z_dim,1)
         real,intent(in) :: mu_bdy(bdy_first_dim,wrf_w_dim,wrf_t_dim+1)
         real,intent(in) :: mub_bdy(bdy_first_dim,wrf_w_dim)
         real,intent(in) :: t_bdy(bdy_first_dim,wrf_z_dim,wrf_w_dim, &
              wrf_t_dim+1)
         real,intent(in) :: qvapor_bdy(bdy_first_dim,wrf_z_dim,wrf_w_dim, &
              wrf_t_dim+1)
         real,intent(inout) :: wrf_pressures_bdy(bdy_first_dim,wrf_z_dim, &
              wrf_w_dim, wrf_t_dim+1)

         ! Local variables
         real :: geopotential_above, geopotential_below, alpha_dry,theta_moist
         integer :: i,j,k,n

         do n = 1, wrf_t_dim+1
            do j = 1, wrf_w_dim
               do k = 1, wrf_z_dim
                  do i = 1, bdy_first_dim
                     geopotential_above = &
                          calc_geopotential(ph_bdy(i,k+1,j,n),phb_bdy(i,k+1,j))
                     geopotential_below = calc_geopotential(ph_bdy(i,k,j,n), &
                          phb_bdy(i,k,j))
                     alpha_dry = &
                          calc_inverse_dryair_density(geopotential_above, &
                          geopotential_below, rdnw(k,1), mu_bdy(i,j,n), &
                          mub_bdy(i,j))
                     theta_moist = &
                          calc_moist_potential_temperature(t_bdy(i,k,j,n), &
                          qvapor_bdy(i,k,j,n))
                     wrf_pressures_bdy(i,k,j,n) = &
                          calc_full_pressure(theta_moist,alpha_dry)
                  end do
               end do
            end do
         end do
      end subroutine calc_wrfbdy_pressure_internal
   end subroutine calc_wrfbdy_pressures_single_zone

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  write_to_wrfinput
   !
   ! DESCRIPTION:  Writes interpolated GOCART variable to wrfinput file.
   !
   !---------------------------------------------------------------------------

   subroutine write_to_wrfinput(this,wrf_domain_id, &
        gocart_variables_final_interp, var_name)

      ! Arguments
      type(interpolator),intent(in) :: this
      integer,intent(in) :: wrf_domain_id
      real,intent(in) :: gocart_variables_final_interp( &
           this%wrfinput(wrf_domain_id)%wrf_x_dim, &
           this%wrfinput(wrf_domain_id)%wrf_y_dim, & 
           this%wrfinput(wrf_domain_id)%wrf_z_dim, &
           1)
      character(len=*),intent(in) :: var_name

      ! Local variables
      integer :: ncid, var_id
      integer :: wrf_x_dim, wrf_y_dim,wrf_z_dim

      ! For convenience, copy certain data members to local variables
      wrf_x_dim = this%wrfinput(wrf_domain_id)%wrf_x_dim
      wrf_y_dim = this%wrfinput(wrf_domain_id)%wrf_y_dim
      wrf_z_dim = this%wrfinput(wrf_domain_id)%wrf_z_dim
      ncid = this%wrfinput(wrf_domain_id)%nc_id

      ! Fetch the ID of the defined new netCDF variable, and write the
      ! data.
      var_id = read_netcdf_var_id(ncid,trim(var_name))
      call write_netcdf_real_array_4d(ncid,var_id,wrf_x_dim, wrf_y_dim, &
           wrf_z_dim, 1, gocart_variables_final_interp)

   end subroutine write_to_wrfinput

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  write_to_wrfbdy
   !
   ! DESCRIPTION:  Writes interpolated GOCART variable to wrfbdy file.
   !
   !---------------------------------------------------------------------------

   subroutine write_to_wrfbdy(this,wrf_dim_1,wrf_dim_2,wrf_dim_3,wrf_dim_4, &
        gocart_variables_final_interp,wrf_varname)

      ! Arguments
      type(interpolator),intent(in) :: this
      integer,intent(in) :: wrf_dim_1
      integer,intent(in) :: wrf_dim_2
      integer,intent(in) :: wrf_dim_3
      integer,intent(in) :: wrf_dim_4
      real,intent(in) :: gocart_variables_final_interp(wrf_dim_1,wrf_dim_2, &
           wrf_dim_3,wrf_dim_4)
      character(len=*),intent(in) :: wrf_varname

      ! Local variables
      integer :: ncid
      integer :: var_id

      ! For convenience, copy data member to local variable
      ncid = this%wrfbdy%nc_id

      ! Fetch the ID of the defined new netCDF variable, and write the data.
      var_id = read_netcdf_var_id(ncid,trim(wrf_varname))
      call write_netcdf_real_array_4d(ncid,var_id,wrf_dim_1,wrf_dim_2, &
           wrf_dim_3,wrf_dim_4,gocart_variables_final_interp)

   end subroutine write_to_wrfbdy

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  destroy_interpolator
   !
   ! DESCRIPTION:  Public destructor method for interpolator data type.
   !
   !---------------------------------------------------------------------------

   subroutine destroy_interpolator(this)

      ! Arguments
      type(interpolator), intent(inout) :: this

      ! Local variables
      integer :: i

      ! Free memory associated with wrfinput data
      do i = 1, this%num_wrf_domains
         if (allocated(this%wrfinput(i)%wrf_latitudes)) &
              deallocate(this%wrfinput(i)%wrf_latitudes)
         if (allocated(this%wrfinput(i)%wrf_longitudes)) &
              deallocate(this%wrfinput(i)%wrf_longitudes)
         if (allocated(this%wrfinput(i)%wrf_pressures)) &
              deallocate(this%wrfinput(i)%wrf_pressures)
         if (allocated(this%wrfinput(i)%gocart_horiz_interp_pressures)) &
              deallocate(this%wrfinput(i)%gocart_horiz_interp_pressures)

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

      if (allocated(this%wrfbdy%wrf_latitudes_bxs)) &
           deallocate(this%wrfbdy%wrf_latitudes_bxs)
      if (allocated(this%wrfbdy%wrf_latitudes_bxe)) &
           deallocate(this%wrfbdy%wrf_latitudes_bxe)
      if (allocated(this%wrfbdy%wrf_latitudes_bys)) &
           deallocate(this%wrfbdy%wrf_latitudes_bys)
      if (allocated(this%wrfbdy%wrf_latitudes_bye)) &
           deallocate(this%wrfbdy%wrf_latitudes_bye)

      if (allocated(this%wrfbdy%wrf_longitudes_bxs)) &
           deallocate(this%wrfbdy%wrf_longitudes_bxs)
      if (allocated(this%wrfbdy%wrf_longitudes_bxe)) &
           deallocate(this%wrfbdy%wrf_longitudes_bxe)
      if (allocated(this%wrfbdy%wrf_longitudes_bys)) &
           deallocate(this%wrfbdy%wrf_longitudes_bys)
      if (allocated(this%wrfbdy%wrf_longitudes_bye)) &
           deallocate(this%wrfbdy%wrf_longitudes_bye)

      if (allocated(this%wrfbdy%wrf_pressures_bxs)) &
           deallocate(this%wrfbdy%wrf_pressures_bxs)
      if (allocated(this%wrfbdy%wrf_pressures_bxe)) &
           deallocate(this%wrfbdy%wrf_pressures_bxe)
      if (allocated(this%wrfbdy%wrf_pressures_bys)) &
           deallocate(this%wrfbdy%wrf_pressures_bys)
      if (allocated(this%wrfbdy%wrf_pressures_bye)) &
           deallocate(this%wrfbdy%wrf_pressures_bye)

      if (allocated(this%wrfbdy%gocart_horiz_interp_pressures_bxs)) &
           deallocate(this%wrfbdy%gocart_horiz_interp_pressures_bxs)
      if (allocated(this%wrfbdy%gocart_horiz_interp_pressures_bxe)) &
           deallocate(this%wrfbdy%gocart_horiz_interp_pressures_bxe)
      if (allocated(this%wrfbdy%gocart_horiz_interp_pressures_bys)) &
           deallocate(this%wrfbdy%gocart_horiz_interp_pressures_bys)
      if (allocated(this%wrfbdy%gocart_horiz_interp_pressures_bye)) &
           deallocate(this%wrfbdy%gocart_horiz_interp_pressures_bye)

      ! Miscellaneous
      if (allocated(this%datetimes)) deallocate(this%datetimes)
      if (allocated(this%for_wrfbdy)) deallocate(this%for_wrfbdy)
      if (allocated(this%for_wrfinput)) deallocate(this%for_wrfinput)
      if (allocated(this%gocart_longitudes)) &
           deallocate(this%gocart_longitudes)
      if (allocated(this%gocart_latitudes)) &
           deallocate(this%gocart_latitudes)

      if (allocated(this%aerosol_names)) deallocate(this%aerosol_names)
      if (allocated(this%aerosol_units)) deallocate(this%aerosol_units)
      if (allocated(this%aerosol_descriptions)) &
           deallocate(this%aerosol_descriptions)

      ! Take care of non-allocatable data members
      call init_interpolator(this)

   end subroutine destroy_interpolator

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  update_wrfinput_header
   !
   ! DESCRIPTION:  Modifies wrfinput header to define new GOCART variables.
   !
   !---------------------------------------------------------------------------

   subroutine update_wrfinput_header(this,ncid)

      ! Arguments
      type(interpolator),intent(in) :: this
      integer,intent(in) :: ncid

      ! Local variables
      integer,dimension(1) :: field_type
      character(len=80) :: units,description,name
      logical :: is_netcdf4
      integer :: storage
      integer :: chunksizes(4)
      integer :: shuffle,deflate,deflate_level
      integer :: var_id
      integer :: ii

      ! Get netCDF4 chunking and compression information from file.
      call get_real_array_4d_netcdf4_settings(ncid,'T',is_netcdf4, &
           storage,chunksizes,shuffle,deflate,deflate_level)

      ! Go into define mode
      call set_netcdf_define_mode(ncid)
      field_type(1) = 104

      ! Loop through aerosol species
      do ii = 1, this%num_aerosols
         name=trim(this%aerosol_names(ii))
         units=trim(this%aerosol_units(ii))
         description=trim(this%aerosol_descriptions(ii))
         print*,'Defining ',trim(name)
         !         var_id = define_netcdf_real_array_4d(ncid,trim(name), &
         !              "west_east","south_north","bottom_top","Time")
         ! Allow for cases where array already exists with same dimensions
         ! and data type.
         var_id = check_define_netcdf_real_array_4d(ncid,trim(name), &
              "west_east","south_north","bottom_top","Time")
         if (is_netcdf4) then
            call define_real_array_4d_netcdf4_settings(ncid,var_id, &
                 storage,chunksizes,shuffle,deflate,deflate_level)
         end if
         call write_netcdf_integer_attribute(ncid,var_id,"FieldType", &
              1,field_type)
         call write_netcdf_text_attribute(ncid,var_id,"MemoryOrder", &
              3,"XYZ")
         call write_netcdf_text_attribute(ncid,var_id,"description", &
              len_trim(description),trim(description))
         call write_netcdf_text_attribute(ncid,var_id,"units", &
              len_trim(units),trim(units))
         call write_netcdf_text_attribute(ncid,var_id,"stagger", &
              1," ")
         call write_netcdf_text_attribute(ncid,var_id,"coordinates", &
              10,"XLONG XLAT")
      end do

      ! Exit define mode and update header
      call unset_netcdf_define_mode(ncid)

   end subroutine update_wrfinput_header

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  update_wrfbdy_header
   !
   ! DESCRIPTION:  Modifies wrfbdy header to define new GOCART variables.
   !
   !---------------------------------------------------------------------------

   subroutine update_wrfbdy_header(this,ncid)

      ! Arguments
      type(interpolator),intent(in) :: this
      integer,intent(in) :: ncid

      ! Local variables
      character(len=80) :: units,description
      character(len=132) :: varname
      integer :: field_type(1) = 104
      character(len=5) :: suffix(8)
      character(len=3) :: memory_order(8)
      integer :: var_id
      logical :: is_netcdf4
      integer :: storage
      integer :: chunksizes(4)
      integer :: shuffle,deflate,deflate_level
      integer :: i,j

      ! Go into netCDF define mode
      call set_netcdf_define_mode(ncid)

      ! Loop through each aerosol species and boundary zone, and define
      ! the new variable.
      suffix(1) = "BXS"  ; memory_order(1) = "XSZ"
      suffix(2) = "BXE"  ; memory_order(2) = "XEZ"
      suffix(3) = "BTXS" ; memory_order(3) = "XSZ"
      suffix(4) = "BTXE" ; memory_order(4) = "XEZ"
      suffix(5) = "BYS"  ; memory_order(5) = "YSZ"
      suffix(6) = "BYE"  ; memory_order(6) = "YEZ"
      suffix(7) = "BTYS" ; memory_order(7) = "YSZ"
      suffix(8) = "BTYE" ; memory_order(8) = "YEZ"
      do i = 1, this%num_aerosols
         units=trim(this%aerosol_units(i))
         description=trim(this%aerosol_descriptions(i))
         do j = 1, 8

            ! Get netCDF4 chunking and compression information from file.
            call get_real_array_4d_netcdf4_settings(ncid, &
                 'U_'//trim(suffix(j)),is_netcdf4, &
                 storage,chunksizes,shuffle,deflate,deflate_level)

            varname = trim(this%aerosol_names(i))//"_"//trim(suffix(j))
            print*,'Defining ',trim(varname)
            if (j < 5) then
               !               var_id = define_netcdf_real_array_4d(ncid,&
               !                    trim(varname),&
               !                    "south_north","bottom_top","bdy_width","Time")
               ! Allow for cases where array already exists with correct
               ! dimensions and data type.
               var_id = check_define_netcdf_real_array_4d(ncid,&
                    trim(varname),&
                    "south_north","bottom_top","bdy_width","Time")
            else
               !               var_id = define_netcdf_real_array_4d(ncid,&
               !                    trim(varname), &
               !                    "west_east","bottom_top","bdy_width","Time")
               ! Allow for cases where array already exists with correct
               ! dimensions and data type.
               var_id = check_define_netcdf_real_array_4d(ncid,&
                    trim(varname), &
                    "west_east","bottom_top","bdy_width","Time")
            end if
            if (is_netcdf4) then
               call define_real_array_4d_netcdf4_settings(ncid,var_id, &
                    storage,chunksizes,shuffle,deflate,deflate_level)
            end if
            call write_netcdf_integer_attribute(ncid,var_id,"FieldType", &
                 1,field_type)
            call write_netcdf_text_attribute(ncid,var_id,"MemoryOrder", &
                 3,trim(memory_order(j)))
            call write_netcdf_text_attribute(ncid,var_id,"description", &
                 len_trim(description),trim(description))
            call write_netcdf_text_attribute(ncid,var_id,"units", &
                 len_trim(units),trim(units))
            call write_netcdf_text_attribute(ncid,var_id,"stagger", &
                 1," ")
         end do
      end do

      ! Exit netCDF define mode and update the header
      call unset_netcdf_define_mode(ncid)

   end subroutine update_wrfbdy_header

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  write_test_netcdf_file
   !
   ! DESCRIPTION:  Writes a 4D real array to a new netCDF file.  Intended
   ! for debugging.
   !
   !---------------------------------------------------------------------------

   subroutine write_test_netcdf_file(dim1,dim2,dim3,dim4,var4d,varname, &
        filename)

      ! Arguments
      integer,intent(in) :: dim1
      integer,intent(in) :: dim2
      integer,intent(in) :: dim3
      integer,intent(in) :: dim4
      real,intent(in) :: var4d(dim1,dim2,dim3,dim4)
      character(len=*),intent(in) :: varname
      character(len=*),intent(in) :: filename

      ! Local variables
      integer :: nc_id, var_id
      integer :: dim1_id, dim2_id, dim3_id, dim4_id

      ! Open the file
      nc_id = open_netcdf_newfile(trim(filename))

      ! Define the dimensions
      dim1_id = write_netcdf_dimension(nc_id,"dim1",dim1)
      dim2_id = write_netcdf_dimension(nc_id,"dim2",dim2)
      dim3_id = write_netcdf_dimension(nc_id,"dim3",dim3)
      dim4_id = write_netcdf_dimension(nc_id,"dim4",dim4)

      ! Define the array
      var_id = define_netcdf_real_array_4d(nc_id,trim(varname), &
           "dim1","dim2","dim3","dim4")

      ! Exit define mode
      call unset_netcdf_define_mode(nc_id)

      ! Write the array
      call write_netcdf_real_array_4d(nc_id,var_id,dim1,dim2,dim3,dim4,var4d)

      ! Close the file
      call close_netcdf_file(nc_id)

   end subroutine write_test_netcdf_file

end module interpolator_mod

