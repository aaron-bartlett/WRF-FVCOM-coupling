!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration and Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  interp_util_mod
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Contains utility functions for interpolating data to wrfinput and wrfbdy
! grids.
!------------------------------------------------------------------------------

module interp_util_mod

   ! Import modules
   use sorted_datetimes_mod

   ! Reset defaults
   implicit none
   private

   ! Public parameters
   integer,parameter,public :: MAX_DIR_LEN = 132
   integer,parameter,public :: MAX_NAME_LEN = 10
   integer,parameter,public :: MAX_SOURCE_LEN = 8

   ! Public routines
   public :: build_gocart_filename
   public :: build_offline_filename
   public :: build_merraero_filename
   public :: build_merra2_filename
   public :: flip_vertical_levels_4d
   public :: flip_vertical_levels_3d
   public :: calc_pressure_geos5_gocart
   public :: calc_pressure_offline_gocart
   public :: find_gocart_latlon_indices
   public :: bilinear_interp
   public :: horiz_interpolation
   public :: horiz_interpolation_wrfbdy
   public :: vert_interpolation
   public :: vert_interpolation_wrfbdy
   public :: decouple_wrfbdy
   public :: calc_delta_time_wrfbdy

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  build_gocart_filename
   !
   ! DESCRIPTION:  Function for building a gocart (GEOS-5 netCDF4) filename
   ! based on requested WRF date/time, file directory, and file prefix.
   !
   !---------------------------------------------------------------------------

   function build_gocart_filename(wrf_datetime,gocart_dir,gocart_prefix) &
        result(gocart_filename)

      ! Arguments
      character(len=MAX_DATE_LEN), intent(in) :: wrf_datetime
      character(len=*), intent(in) :: gocart_dir
      character(len=*), intent(in) :: gocart_prefix

      ! Return variable
      character(len=256) :: gocart_filename

      ! Local variables
      character(len=4) :: cyear
      character(len=2) :: cmonth,cday,chour,cminute,csecond

      ! Build gocart filename 
      read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear,cmonth,cday,chour,cminute,csecond

      gocart_filename = trim(gocart_dir)//"/"//trim(gocart_prefix)//"." &
           //cyear//cmonth//cday//"_"//chour//cminute//"z.nc4"

      return
   end function build_gocart_filename

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  build_offline_filename
   !
   ! DESCRIPTION:  Function for building an offline GOCART netCDF filename
   ! based on requested WRF date/time, file directory, file prefix, and
   ! collection name.
   !
   !---------------------------------------------------------------------------

   function build_offline_filename(wrf_datetime,gocart_dir,gocart_prefix, &
        gocart_varname) result(gocart_filename)

      ! Arguments
      character(len=MAX_DATE_LEN), intent(in) :: wrf_datetime
      character(len=*), intent(in) :: gocart_dir
      character(len=*), intent(in) :: gocart_prefix
      character(len=3), intent(in) :: gocart_varname

      ! Return variable
      character(len=256) :: gocart_filename

      ! Local variables
      character(len=2) :: collection_name
      character(len=4) :: cyear
      character(len=2) :: cmonth,cday,chour,cminute,csecond
      integer :: iyear,imonth,iday,ihour,iminute,isecond

      ! Get date and time info
      read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear,cmonth,cday,chour,cminute,csecond
      read(wrf_datetime,"(I4,X,I2,X,I2,x,I2,x,I2,X,I2)") &                    
           iyear,imonth,iday,ihour,iminute,isecond

      ! Make sure requested time is valid.
      if (iminute .ne. 0 .or. isecond .ne. 0 .or. &
          (ihour .ne. 0 .and. mod(ihour,3) .ne. 0)) then
         write(6,*)'ERROR, invalid time requested!'
         write(6,*)'Offline GOCART only available 3-hourly from 00Z!'
         write(6,*)'Requested ',chour//":"//cminute//":"//csecond
         stop
      end if

      ! Select appropriate collection
      select case (gocart_varname)
         case ("bc1","oc1","bc2","oc2")
            collection_name = "cc"
         case ("du1","du2","du3","du4","du5")
            collection_name = "du"
         case ("ss1","ss2","ss3","ss4")
            collection_name = "ss"
         case ("dms","so2","so4","msa")
            collection_name = "su"
         case default
            write(6,*)'ERROR, unknown GOCART variable ',trim(gocart_varname)
            stop
      end select
      
      ! Construct filename
      gocart_filename = trim(gocart_dir)//"/" &
           //trim(gocart_prefix)//"_"//collection_name//"_tracer_3HINST_" &
              //cyear//cmonth//cday//".nc"

   end function build_offline_filename

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  build_merraero_filename
   !
   ! DESCRIPTION:  Function for building a GOCART MERRAero netCDF4 filename
   ! based on requested WRF date/time and file directory. The file prefix
   ! is standardized in this routine.
   !
   !---------------------------------------------------------------------------

   function build_merraero_filename(wrf_datetime,gocart_dir) &
        result(gocart_filename)

      ! Arguments
      character(len=MAX_DATE_LEN), intent(in) :: wrf_datetime
      character(len=*), intent(in) :: gocart_dir

      ! Return variable
      character(len=256) :: gocart_filename

      ! Local variables
      character(len=4) :: cyear
      character(len=2) :: cmonth,cday,chour,cminute,csecond
      integer :: iyear,imonth,iday,ihour,iminute,isecond

      ! Extract date and time info.
      read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear,cmonth,cday,chour,cminute,csecond
      read(wrf_datetime,"(I4,X,I2,X,I2,x,I2,x,I2,X,I2)") &                    
           iyear,imonth,iday,ihour,iminute,isecond

      ! Make sure requested time is valid.
      if (iminute .ne. 0 .or. isecond .ne. 0 .or. &
          (ihour .ne. 0 .and. mod(ihour,3) .ne. 0)) then
         write(6,*)'ERROR, invalid time requested!'
         write(6,*)'MERRAero only available 3-hourly from 00Z!'
         write(6,*)'Requested ',chour//":"//cminute//":"//csecond
         stop
      end if

      ! Create the file name.
      gocart_filename = trim(gocart_dir)//"/Y"//cyear//'/M'//cmonth//'/' &
           //"dR_MERRA-AA-r2.inst3d_aer_v." &
              //cyear//cmonth//cday//"_"//chour//cminute//"z.nc4"

      return
   end function build_merraero_filename

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  build_merra2_filename
   !
   ! DESCRIPTION:  Function for building a GOCART MERRA2 netCDF4 filename
   ! based on requested WRF date/time and file directory. The file prefix
   ! is standardized in this routine.
   !
   !---------------------------------------------------------------------------

   function build_merra2_filename(wrf_datetime,gocart_dir) &
        result(gocart_filename)

      ! Arguments
      character(len=MAX_DATE_LEN), intent(in) :: wrf_datetime
      character(len=*), intent(in) :: gocart_dir

      ! Return variable
      character(len=256) :: gocart_filename

      ! Local variables
      character(len=10) :: merra2_stream
      character(len=256) :: myfile
      character(len=4) :: cyear
      character(len=2) :: cmonth,cday,chour,cminute,csecond
      integer :: iyear,imonth,iday,ihour,iminute,isecond

      ! Extract date and time info.
      read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear,cmonth,cday,chour,cminute,csecond
      read(wrf_datetime,"(I4,X,I2,X,I2,x,I2,x,I2,X,I2)") &                    
           iyear,imonth,iday,ihour,iminute,isecond

      ! Make sure requested time is valid.
      if (iminute .ne. 0 .or. isecond .ne. 0 .or. &
          (ihour .ne. 0 .and. mod(ihour,3) .ne. 0)) then
         write(6,*)'ERROR, invalid time requested!'
         write(6,*)'MERRA2 only available 3-hourly from 00Z!'
         write(6,*)'Requested ',chour//":"//cminute//":"//csecond
         stop
      end if

      ! Find the correct MERRA2 Stream
      if ( iyear < 1992 ) then 
         merra2_stream = 'MERRA2_100'
      elseif ( iyear < 2001) then 
         merra2_stream = 'MERRA2_200'
      elseif ( iyear < 2011 ) then 
         merra2_stream = 'MERRA2_300'
      else
         merra2_stream = 'MERRA2_400'
      endif

      ! Create the file name
      myfile = merra2_stream//'.inst3_3d_aer_Nv.' &
           //cyear//cmonth//cday//'.nc4'
      gocart_filename = trim(gocart_dir)//'/'//merra2_stream//'/' &
           //'/stage/Y'//cyear//'/M'//cmonth//'/'//trim(myfile)

      return
   end function build_merra2_filename

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  flip_vertical_levels_3d
   !
   ! DESCRIPTION:  Private method to flip GEOS5 data.  GEOS5 output is
   ! written to file with model top at k=1.  We'll switch this to match
   ! the WRF convention (k=1 is ground level).
   !
   !---------------------------------------------------------------------------

   subroutine flip_vertical_levels_3d(dim1,dim2,dim3,array3d)

      ! Arguments
      integer, intent(in) :: dim1, dim2, dim3
      real,dimension(dim1,dim2,dim3),intent(inout) :: array3d

      ! Local variables
      real, allocatable :: tmp3d(:,:,:)
      integer :: i,j,k
      allocate(tmp3d(dim1,dim2,dim3))
      tmp3d = 0.0
      do k = 1, dim3
         do j = 1, dim2
            do i = 1, dim1
               tmp3d(i,j,dim3-k+1) = array3d(i,j,k)
            end do
         end do
      end do
      do k = 1, dim3
         do j = 1, dim2
            do i = 1, dim1
               array3d(i,j,k) = tmp3d(i,j,k)
            end do
         end do
      end do
      deallocate(tmp3d)
   end subroutine flip_vertical_levels_3d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  flip_vertical_levels_4d
   !
   ! DESCRIPTION:  Private method to flip GEOS5 data.  GEOS5 output is
   ! written to file with model top at k=1.  We'll switch this to match
   ! the WRF convention (k=1 is ground level).
   !
   !---------------------------------------------------------------------------

   subroutine flip_vertical_levels_4d(dim1,dim2,dim3,dim4,array4d)

      ! Arguments
      integer, intent(in) :: dim1, dim2, dim3, dim4
      real,dimension(dim1,dim2,dim3,dim4),intent(inout) :: array4d

      ! Local variables
      real, allocatable :: tmp4d(:,:,:,:)
      integer :: i,j,k,n
      allocate(tmp4d(dim1,dim2,dim3,dim4))
      tmp4d = 0.0
      do n = 1, dim4
         do k = 1, dim3
            do j = 1, dim2
               do i = 1, dim1
                  tmp4d(i,j,dim3-k+1,n) = array4d(i,j,k,n)
               end do
            end do
         end do
      end do
      do n = 1, dim4
         do k = 1, dim3
            do j = 1, dim2
               do i = 1, dim1
                  array4d(i,j,k,n) = tmp4d(i,j,k,n)
               end do
            end do
         end do
      end do
      deallocate(tmp4d)
   end subroutine flip_vertical_levels_4d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calc_pressure_geos5_gocart
   !
   ! DESCRIPTION:  Private method for calculating 3D pressure field from 
   ! GEOS5 data.
   !
   !---------------------------------------------------------------------------

   subroutine calc_pressure_geos5_gocart(gocart_lon_dim,gocart_lat_dim, &
        gocart_lev_dim,gocart_time_dim,gocart_ps,gocart_delp,gocart_pressure)

      ! Arguments
      integer, intent(in) :: gocart_lon_dim
      integer, intent(in) :: gocart_lat_dim 
      integer, intent(in) :: gocart_lev_dim
      integer, intent(in) :: gocart_time_dim
      real,intent(in) ::   gocart_ps(gocart_lon_dim,gocart_lat_dim, &
           gocart_time_dim)
      real,intent(in) :: gocart_delp(gocart_lon_dim,gocart_lat_dim, &
           gocart_lev_dim,gocart_time_dim)
      real,intent(inout) :: gocart_pressure(gocart_lon_dim,gocart_lat_dim, &
           gocart_lev_dim,gocart_time_dim)

      ! Local variables
      real :: p_bottom, p_top
      integer :: ii,jj,kk,nn

      do nn = 1, gocart_time_dim
         do jj = 1, gocart_lat_dim
            do ii = 1, gocart_lon_dim
               p_bottom = gocart_ps(ii,jj,nn)
               do kk = 1, gocart_lev_dim
                  p_top = p_bottom - gocart_delp(ii,jj,kk,nn)
                  gocart_pressure(ii,jj,kk,nn) = (p_bottom + p_top) * 0.5
                  p_bottom = p_top
               end do
            end do
         end do
      end do

   end subroutine calc_pressure_geos5_gocart

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calc_pressure_offline_gocart
   !
   ! DESCRIPTION:  Private method for calculating 3D pressure field from 
   ! offline GOCART data
   !
   !---------------------------------------------------------------------------

   subroutine calc_pressure_offline_gocart(gocart_lon_dim,gocart_lat_dim, &
        gocart_lev_dim,gocart_time_dim,gocart_ps,gocart_p0,gocart_a,gocart_b, &
        gocart_pressures)

      ! Arguments
      integer, intent(in) :: gocart_lon_dim
      integer, intent(in) :: gocart_lat_dim 
      integer, intent(in) :: gocart_lev_dim
      integer, intent(in) :: gocart_time_dim
      real,intent(in) ::   gocart_ps(gocart_lon_dim,gocart_lat_dim, &
           gocart_time_dim)
      real,intent(in) :: gocart_p0
      real,intent(in) :: gocart_a(gocart_lev_dim)
      real,intent(in) :: gocart_b(gocart_lev_dim)
      real,intent(inout) :: gocart_pressures(gocart_lon_dim,gocart_lat_dim, &
           gocart_lev_dim,gocart_time_dim)

      ! Local variables
      integer :: ii,jj,kk,nn

      do nn = 1, gocart_time_dim
         do kk = 1, gocart_lev_dim
            do jj = 1, gocart_lat_dim
               do ii = 1, gocart_lon_dim               
                  gocart_pressures(ii,jj,kk,nn) = &
                       (gocart_p0*gocart_a(kk)) + &
                       (gocart_ps(ii,jj,nn)*gocart_b(kk))
               end do
            end do
         end do
      end do

   end subroutine calc_pressure_offline_gocart

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  find_gocart_latlon_indices
   !
   ! DESCRIPTION:  Given arrays of GOCART latitudes and longitudes, find
   ! array indices matching the southwest corner of four points surrounding
   ! the WRF grid column latitude/longitude.
   !
   !---------------------------------------------------------------------------

   subroutine find_gocart_latlon_indices(gocart_lon_dim, gocart_lat_dim, &
        gocart_longitudes, gocart_latitudes, wrf_longitude, wrf_latitude, &
        i_gocart_longitude_w, j_gocart_latitude_s)

      ! Arguments
      integer, intent(in) :: gocart_lon_dim
      integer, intent(in) :: gocart_lat_dim
      real,intent(in) :: gocart_longitudes(gocart_lon_dim)
      real,intent(in) :: gocart_latitudes(gocart_lat_dim)
      real,intent(in) :: wrf_longitude
      real,intent(in) :: wrf_latitude
      integer,intent(inout) :: i_gocart_longitude_w
      integer,intent(inout) :: j_gocart_latitude_s

      ! Local variables
      integer :: ii,jj

      i_gocart_longitude_w = gocart_lon_dim
      do ii = 1, gocart_lon_dim-1
         if (gocart_longitudes(ii  ) .le. wrf_longitude .and. &
              gocart_longitudes(ii+1) .gt. wrf_longitude) then
            i_gocart_longitude_w = ii
         end if
      end do

      j_gocart_latitude_s = gocart_lat_dim
      do jj = 1, gocart_lat_dim-1
         if (gocart_latitudes(jj  ) .le. wrf_latitude .and. &
              gocart_latitudes(jj+1) .gt. wrf_latitude) then
            j_gocart_latitude_s = jj
         end if
      end do

   end subroutine find_gocart_latlon_indices

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  bilinear_interp
   !
   ! DESCRIPTION:  Performs bilinear interpolation from a regular grid 
   ! (e.g., GOCART) to an arbitrary grid point (e.g., WRF grid column).
   ! See http://en.wikipedia.org/wiki/Bilinear_interpolation for an
   ! in-depth description.
   !
   ! NOTE:  To eliminate divisions, this function uses an alternative 
   ! but equivalent formula that assumes a unit square (see nonlinear
   ! section on Wikipedia article):
   ! 
   ! Interpolant = b1 + b2*x + b3*y + b4*x*y
   ! where:  b1 = f(0,0),
   !         b2 = f(1,0) - f(0,0), 
   !         b3 = f(0,1) - f(0,0),
   !         b4 = f(0,0) - f(1,0) - f(0,1) + f(1,1)
   !
   ! The function will convert the input x,y coordinates to the unit square
   ! convention by using x1, y1, and the inverse of the delta-x and delta-y
   ! between corners.  The *inverses* are used to eliminate division 
   ! entirely in the function.
   !
   !---------------------------------------------------------------------------

   function bilinear_interp(x1,inverse_delta_x,x, &
        y1,inverse_delta_y,y,q11,q21,q12,q22) result(q)

      ! Arguments
      ! First corner x-coordinate
      real,intent(in) :: x1 
      ! Inverse distance to second corner x-coordinate
      real,intent(in) :: inverse_delta_x 
      ! Interpolation point x-coordinate
      real,intent(in) :: x  
      ! First corner y-coordinate
      real,intent(in) :: y1 
      ! Inverse distance to second cornder y-coordinate
      real,intent(in) :: inverse_delta_y
      ! Interpolation point y-coordinate
      real,intent(in) :: y  
      ! Quantity at x1,y1
      real,intent(in) :: q11 
      ! Quantity at x2,y1
      real,intent(in) :: q21 
      ! Quantity at x1,y2
      real,intent(in) :: q12 
      ! Quantity at x2,y2
      real,intent(in) :: q22 

      ! Return variable
      real :: q ! Interpolated value at x,y

      ! Local variables
      real :: x_relative,y_relative
      real :: b1,b2,b3,b4

      x_relative = (x - x1)*inverse_delta_x ! Recasts x as in range [0,1]
      y_relative = (y - y1)*inverse_delta_y ! Recasts y as in range [0,1]

      b1 = q11
      b2 = q21 - q11
      b3 = q12 - q11
      b4 = q11 - q21 - q12 + q22

      q = b1 + &
           b2*x_relative + &
           b3*y_relative + &
           b4*x_relative*y_relative

   end function bilinear_interp

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  horiz_interpolation
   !
   ! DESCRIPTION:  Oversees horizontal interpolation of GOCART data to
   ! WRF grid.
   !
   ! NOTE:  Interpolation is actually done via a lower-level bilinear
   ! interpolation scheme. 
   !
   !---------------------------------------------------------------------------

   subroutine horiz_interpolation(wrf_dim_1, wrf_dim_2, &
        wrf_latitudes, wrf_longitudes, &
        gocart_lon_dim,gocart_lat_dim, &
        gocart_longitudes, gocart_latitudes, &
        gocart_lev_dim, gocart_variable, &
         gocart_interp_variable)

      ! Arguments
      integer,intent(in) :: wrf_dim_1
      integer,intent(in) :: wrf_dim_2
      real,intent(in) :: wrf_latitudes(wrf_dim_1,wrf_dim_2)
      real,intent(in) :: wrf_longitudes(wrf_dim_1,wrf_dim_2)
      integer,intent(in) :: gocart_lon_dim
      integer,intent(in) :: gocart_lat_dim
      real,intent(in) :: gocart_longitudes(gocart_lon_dim)
      real,intent(in) :: gocart_latitudes(gocart_lat_dim)
      integer,intent(in) :: gocart_lev_dim
      real,intent(in) :: gocart_variable(gocart_lon_dim,gocart_lat_dim, &
           gocart_lev_dim)
      real,intent(inout) :: gocart_interp_variable(wrf_dim_1,wrf_dim_2, &
           gocart_lev_dim)

      ! Local variables
      real :: inverse_gocart_delta_lon
      real :: inverse_gocart_delta_lat
      real :: wrf_latitude, wrf_longitude
      real :: gocart_longitude_w, gocart_latitude_s
      real :: q_sw,q_se,q_nw,q_ne,q_interp
      integer :: i_gocart_longitude_w, j_gocart_latitude_s
      integer :: ii_w, ii_e, jj_s, jj_n
      integer :: ii,jj,kk

      do jj = 1, wrf_dim_2
         do ii = 1, wrf_dim_1
            wrf_latitude = wrf_latitudes(ii,jj)
            wrf_longitude = wrf_longitudes(ii,jj)

            call find_gocart_latlon_indices(gocart_lon_dim, &
                 gocart_lat_dim, gocart_longitudes, gocart_latitudes, &
                 wrf_longitude, wrf_latitude, &
                 i_gocart_longitude_w, j_gocart_latitude_s)

            ii_w = i_gocart_longitude_w
            gocart_longitude_w = gocart_longitudes(ii_w)
            if (i_gocart_longitude_w == gocart_lon_dim) then
               ii_e = 1
            else
               ii_e = ii_w + 1
            end if

            jj_s = j_gocart_latitude_s
            jj_n = jj_s + 1
            gocart_latitude_s = gocart_latitudes(jj_s)

            inverse_gocart_delta_lon = &
                 1./(gocart_longitudes(ii_e) - gocart_longitudes(ii_w))
            inverse_gocart_delta_lat = &
                 1./(gocart_latitudes(jj_n) - gocart_latitudes(jj_s))

            ! Loop through each gocart level and interpolate pressure
            ! to WRF grid
            do kk = 1, gocart_lev_dim
               q_sw = gocart_variable(ii_w,jj_s,kk)
               q_se = gocart_variable(ii_e,jj_s,kk)
               q_nw = gocart_variable(ii_w,jj_n,kk)
               q_ne = gocart_variable(ii_e,jj_n,kk)

               q_interp = bilinear_interp(gocart_longitude_w, &
                    inverse_gocart_delta_lon, &
                    wrf_longitude, &
                    gocart_latitude_s, &
                    inverse_gocart_delta_lat,wrf_latitude, &
                    q_sw,q_se,q_nw,q_ne)

               gocart_interp_variable(ii,jj,kk) = q_interp

            end do ! kk
         end do ! jj
      end do ! ii

   end subroutine horiz_interpolation

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  horiz_interpolation_wrfbdy
   !
   ! DESCRIPTION:  Oversees horizontal interpolation of GOCART data to
   ! WRF boundary zone.  
   !
   ! NOTE:  Interpolation is actually done via a lower-level bilinear
   ! interpolation scheme. 
   !
   !---------------------------------------------------------------------------

   subroutine horiz_interpolation_wrfbdy(wrfbdy_dim_1, wrfbdy_w_dim, &
        wrfbdy_latitudes, wrfbdy_longitudes, &
        gocart_lon_dim, gocart_lat_dim, &  
        gocart_longitudes, gocart_latitudes, &
        gocart_lev_dim, gocart_variable, &
        gocart_interp_variable)

      ! Arguments
      integer,intent(in) :: wrfbdy_dim_1
      integer,intent(in) :: wrfbdy_w_dim
      real,intent(in) :: wrfbdy_latitudes(wrfbdy_dim_1,wrfbdy_w_dim)
      real,intent(in) :: wrfbdy_longitudes(wrfbdy_dim_1,wrfbdy_w_dim)
      integer,intent(in) :: gocart_lon_dim
      integer,intent(in) :: gocart_lat_dim
      real,intent(in) :: gocart_longitudes(gocart_lon_dim)
      real,intent(in) :: gocart_latitudes(gocart_lat_dim)
      integer,intent(in) :: gocart_lev_dim
      real,intent(in) :: gocart_variable(gocart_lon_dim,gocart_lat_dim, &
           gocart_lev_dim)
      real,intent(inout) :: gocart_interp_variable(wrfbdy_dim_1, &
           gocart_lev_dim, wrfbdy_w_dim)

      ! Local variables
      real :: inverse_gocart_delta_lon
      real :: inverse_gocart_delta_lat
      real :: wrfbdy_latitude, wrfbdy_longitude
      real :: gocart_longitude_w, gocart_latitude_s
      real :: q_sw,q_se,q_nw,q_ne,q_interp
      integer :: i_gocart_longitude_w, j_gocart_latitude_s
      integer :: ii_w, ii_e, jj_s, jj_n
      integer :: ii,jj,kk

      do jj = 1, wrfbdy_w_dim
         do ii = 1, wrfbdy_dim_1
            wrfbdy_latitude = wrfbdy_latitudes(ii,jj)
            wrfbdy_longitude = wrfbdy_longitudes(ii,jj)

            call find_gocart_latlon_indices(gocart_lon_dim, &
                 gocart_lat_dim, gocart_longitudes, gocart_latitudes, &
                 wrfbdy_longitude, wrfbdy_latitude, &
                 i_gocart_longitude_w, j_gocart_latitude_s)

            ii_w = i_gocart_longitude_w
            gocart_longitude_w = gocart_longitudes(ii_w)
            if (i_gocart_longitude_w == gocart_lon_dim) then
               ii_e = 1
            else
               ii_e = ii_w + 1
            end if

            jj_s = j_gocart_latitude_s
            jj_n = jj_s + 1
            gocart_latitude_s = gocart_latitudes(jj_s)

            inverse_gocart_delta_lon = &
                 1./(gocart_longitudes(ii_e) - gocart_longitudes(ii_w))
            inverse_gocart_delta_lat = &
                 1./(gocart_latitudes(jj_n) - gocart_latitudes(jj_s))

            ! Loop through each gocart level and interpolate pressure
            ! to WRF grid
            do kk = 1, gocart_lev_dim
               q_sw = gocart_variable(ii_w,jj_s,kk)
               q_se = gocart_variable(ii_e,jj_s,kk)
               q_nw = gocart_variable(ii_w,jj_n,kk)
               q_ne = gocart_variable(ii_e,jj_n,kk)

               q_interp = bilinear_interp(gocart_longitude_w, &
                    inverse_gocart_delta_lon, &
                    wrfbdy_longitude, &
                    gocart_latitude_s, &
                    inverse_gocart_delta_lat,wrfbdy_latitude, &
                    q_sw,q_se,q_nw,q_ne)

               gocart_interp_variable(ii,kk,jj) = q_interp
            end do ! kk
         end do ! jj
      end do ! ii

   end subroutine horiz_interpolation_wrfbdy

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  vert_interpolation
   !
   ! DESCRIPTION:  Private method for vertically interpolating GOCART data
   ! to WRF grid levels.  Interpolation is linear in log(p)
   !
   !---------------------------------------------------------------------------

   subroutine vert_interpolation(wrf_x_dim,wrf_y_dim,wrf_z_dim,wrf_t_dim, &
        gocart_lev_dim,wrf_pressures,gocart_horiz_interp_pressures, &
        gocart_horiz_interp_var,gocart_final_interp_var)

      ! Arguments
      integer, intent(in) :: wrf_x_dim
      integer, intent(in) :: wrf_y_dim
      integer, intent(in) :: wrf_z_dim
      integer, intent(in) :: wrf_t_dim
      integer, intent(in) :: gocart_lev_dim
      real,intent(in) :: wrf_pressures(wrf_x_dim,wrf_y_dim,wrf_z_dim,wrf_t_dim)
      real,intent(in) :: gocart_horiz_interp_pressures(wrf_x_dim,wrf_y_dim, &
           gocart_lev_dim,wrf_t_dim)
      real,intent(in) :: gocart_horiz_interp_var(wrf_x_dim,wrf_y_dim, &
           gocart_lev_dim,wrf_t_dim)
      real,intent(out) :: gocart_final_interp_var(wrf_x_dim,wrf_y_dim, &
           wrf_z_dim,wrf_t_dim)

      ! Local variables
      real :: gocart_p_base, gocart_log_p_base
      real :: gocart_p_top, gocart_log_p_top
      real :: gocart_var_base, gocart_var_top, gocart_var_interp
      real :: wrf_p, wrf_log_p
      logical :: done
      integer :: ii,jj,kk_wrf,kk_gocart,nn

      do nn = 1, wrf_t_dim
         do jj = 1, wrf_y_dim
            do ii = 1, wrf_x_dim

               do kk_wrf = 1, wrf_z_dim

                  wrf_p = wrf_pressures(ii,jj,kk_wrf,nn)

                  ! Special case:  WRF level is physically below GOCART
                  ! data.  Just use lowest GOCART level.
                  gocart_p_base = gocart_horiz_interp_pressures(ii,jj,1,nn)

                  if (wrf_p .ge. gocart_p_base) then
                     gocart_final_interp_var(ii,jj,kk_wrf,nn) = &
                          gocart_horiz_interp_var(ii,jj,1,nn)
                     cycle
                  end if

                  ! Special case:  WRF level is physically above GOCART
                  ! data.  Just use highest GOCART level.
                  gocart_p_top = &
                       gocart_horiz_interp_pressures(ii,jj,gocart_lev_dim,nn)

                  if (wrf_p .le. gocart_p_top) then
                     gocart_final_interp_var(ii,jj,kk_wrf,nn) = &
                          gocart_horiz_interp_var(ii,jj,gocart_lev_dim,&
                          nn)
                     cycle
                  end if

                  ! General case:  Linearly interpolate in log(p).
                  done = .false.
                  do kk_gocart = 1, gocart_lev_dim-1
                     gocart_p_base = &
                          gocart_horiz_interp_pressures(ii,jj,kk_gocart,nn)
                     gocart_p_top = &
                          gocart_horiz_interp_pressures(ii,jj,kk_gocart+1,nn)

                     if (gocart_p_base .ge. wrf_p .and. &
                          gocart_p_top .lt. wrf_p) then
                        done = .true.
                        exit
                     end if
                  end do
                  if (.not. done) then
                     print*,'ERROR, cannot find wrf level in gocart!'
                     print*,'gocart_p_base = ',gocart_p_base
                     print*,'gocart_p_top = ',gocart_p_top
                     print*,'wrf_p = ',wrf_p
                     print*,'ii,jj,kk_wrf,nn = ',ii,jj,kk_wrf,nn
                     stop
                  end if

                  gocart_log_p_base = log(gocart_p_base)
                  gocart_log_p_top = log(gocart_p_top)
                  wrf_log_p = log(wrf_p)

                  gocart_var_base = &
                       gocart_horiz_interp_var(ii,jj,kk_gocart,nn)
                  gocart_var_top = &
                       gocart_horiz_interp_var(ii,jj,kk_gocart+1,nn)
                  gocart_var_interp = &
                       gocart_var_base + &
                       ( (gocart_var_top - gocart_var_base) * &
                       (wrf_log_p - gocart_log_p_base) / &
                       (gocart_log_p_top - gocart_log_p_base) )

                  gocart_final_interp_var(ii,jj,kk_wrf,nn) = gocart_var_interp

               end do
            end do
         end do
      end do
   end subroutine vert_interpolation

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  vert_interpolation_wrfbdy
   !
   ! DESCRIPTION:  Private method for vertically interpolating GOCART data
   ! to WRF grid levels.  Interpolation is linear in log(p).  For use with
   ! wrfbdy arrays (vertical is second dimension).
   !
   !---------------------------------------------------------------------------

   subroutine vert_interpolation_wrfbdy(wrfbdy_dim_1,wrfbdy_z_dim, &
        wrfbdy_w_dim,wrfbdy_t_dim,gocart_lev_dim,wrf_pressures, &
        gocart_horiz_interp_pressures, &
        gocart_horiz_interp_var,gocart_final_interp_var)

      ! Arguments
      integer, intent(in) :: wrfbdy_dim_1
      integer, intent(in) :: wrfbdy_z_dim
      integer, intent(in) :: wrfbdy_w_dim
      integer, intent(in) :: wrfbdy_t_dim
      integer, intent(in) :: gocart_lev_dim
      real,intent(in) :: wrf_pressures(wrfbdy_dim_1, &
           wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim)
      real,intent(in) :: gocart_horiz_interp_pressures(wrfbdy_dim_1, &
           gocart_lev_dim,wrfbdy_w_dim,wrfbdy_t_dim)
      real,intent(in) :: gocart_horiz_interp_var(wrfbdy_dim_1, &
           gocart_lev_dim,wrfbdy_w_dim,wrfbdy_t_dim)
      real,intent(out) :: gocart_final_interp_var(wrfbdy_dim_1, &
           wrfbdy_z_dim, wrfbdy_w_dim,wrfbdy_t_dim)

      ! Local variables
      real :: gocart_p_base, gocart_log_p_base
      real :: gocart_p_top, gocart_log_p_top
      real :: gocart_var_base, gocart_var_top, gocart_var_interp
      real :: wrf_p, wrf_log_p
      logical :: done
      integer :: ii,jj_wrf,jj_gocart,kk,nn

      do nn = 1, wrfbdy_t_dim
         do kk = 1, wrfbdy_w_dim
            do ii = 1, wrfbdy_dim_1

               do jj_wrf = 1, wrfbdy_z_dim

                  wrf_p = wrf_pressures(ii,jj_wrf,kk,nn)

                  ! Special case:  WRF level is physically below GOCART
                  ! data.  Just use lowest GOCART level.
                  gocart_p_base = gocart_horiz_interp_pressures(ii,1,kk,nn)
                  if (wrf_p .ge. gocart_p_base) then
                     gocart_final_interp_var(ii,jj_wrf,kk,nn) = &
                          gocart_horiz_interp_var(ii,1,kk,nn)
                     cycle
                  end if

                  ! Special case:  WRF level is physically above GOCART
                  ! data.  Just use highest GOCART level.
                  gocart_p_top = &
                       gocart_horiz_interp_pressures(ii,gocart_lev_dim,kk,nn)
                  if (wrf_p .le. gocart_p_top) then
                     gocart_final_interp_var(ii,jj_wrf,kk,nn) = &
                          gocart_horiz_interp_var(ii,gocart_lev_dim,kk,&
                          nn)
                     cycle
                  end if

                  ! General case:  Linearly interpolate in log(p).
                  done = .false.
                  do jj_gocart = 1, gocart_lev_dim-1
                     gocart_p_base = &
                          gocart_horiz_interp_pressures(ii,jj_gocart  ,kk,nn)
                     gocart_p_top = &
                          gocart_horiz_interp_pressures(ii,jj_gocart+1,kk,nn)

                     if (gocart_p_base .ge. wrf_p .and. &
                          gocart_p_top .lt. wrf_p) then
                        done = .true.
                        exit
                     end if
                  end do
                  if (.not. done) then
                     print*,'ERROR, cannot find wrf level in gocart!'
                     print*,'gocart_p_base = ',gocart_p_base
                     print*,'gocart_p_top = ',gocart_p_top
                     print*,'wrf_p = ',wrf_p
                     print*,'ii,jj_wrf,kk,nn = ',ii,jj_wrf,kk,nn
                     stop
                  end if

                  gocart_log_p_base = log(gocart_p_base)
                  gocart_log_p_top = log(gocart_p_top)
                  wrf_log_p = log(wrf_p)

                  gocart_var_base = &
                       gocart_horiz_interp_var(ii,jj_gocart  ,kk,nn)
                  gocart_var_top = &
                       gocart_horiz_interp_var(ii,jj_gocart+1,kk,nn)
                  gocart_var_interp = &
                       gocart_var_base + &
                       ( (gocart_var_top - gocart_var_base) * &
                       (wrf_log_p - gocart_log_p_base) / &
                       (gocart_log_p_top - gocart_log_p_base) )

                  gocart_final_interp_var(ii,jj_wrf,kk,nn) = gocart_var_interp

               end do
            end do
         end do
      end do
   end subroutine vert_interpolation_wrfbdy

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  decouple_wrfbdy
   !
   ! DESCRIPTION:  A generic subroutine to "decouple" a wrfbdy variable--i.e.,
   ! divide out the total dry air mass.
   !
   !---------------------------------------------------------------------------

   subroutine decouple_wrfbdy(length_dim, height_dim, width_dim, time_dim, &
        mu, mub, data)

      ! Arguments
      integer,intent(in) :: length_dim
      integer,intent(in) :: height_dim
      integer,intent(in) :: width_dim
      integer,intent(in) :: time_dim
      real,intent(in) :: mu(length_dim,width_dim,time_dim)
      real,intent(in) :: mub(length_dim,width_dim)
      real,intent(inout) :: data(length_dim,height_dim,width_dim,time_dim)

      ! Local variables
      integer :: i,j,k,n

      do n = 1, time_dim
         do j = 1, width_dim
            do k = 1, height_dim
               do i = 1, length_dim
                  data(i,k,j,n) = data(i,k,j,n) / &
                       (mu(i,j,n) + mub(i,j))
               end do
            end do
         end do
      end do
   end subroutine decouple_wrfbdy

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calc_delta_time_wrfbdy
   !
   ! DESCRIPTION:  A function for calculating the time (in seconds) between
   ! wrfbdy time levels, for use in calculating tendencies.
   !
   !---------------------------------------------------------------------------

   function calc_delta_time_wrfbdy(datetime1, datetime2) result (delta_time)

      ! Arguments
      character(len=MAX_DATE_LEN),intent(in) :: datetime1
      character(len=MAX_DATE_LEN),intent(in) :: datetime2

      ! Return variable
      integer :: delta_time

      ! Local variables
      integer :: order
      integer :: hours(2)

      ! Pull out the two-digit hours and convert to integer
      order = compare_datetimes(datetime1,datetime2)
      if (order > 0) then
         read(datetime1(12:13),'(i2)') hours(1)
         read(datetime2(12:13),'(i2)') hours(2)
      else if (order < 0) then
         read(datetime2(12:13),'(i2)') hours(1)
         read(datetime1(12:13),'(i2)') hours(2)
      else if (order == 0) then
         print*,'ERROR, datetimes are the same!'
         print*,'datetime1 = ',datetime1
         stop
      end if

      delta_time = hours(2) - hours(1)
      if (delta_time < 0) then
         delta_time = (hours(2) - 24) - hours(1)
      end if
      delta_time = delta_time*3600

      return
   end function calc_delta_time_wrfbdy

end module interp_util_mod
