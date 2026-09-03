!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration and Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  interp_util_mod
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
! Modified for CASA CO2 Jossy P. Jacob, SSSO
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

   ! Public routines
   public :: build_casa_filename
   public :: build_casa_emission_filename
   public :: flip_vertical_levels_4d
   public :: calc_pressure_casa
   public :: find_casa_latlon_indices
   public :: bilinear_interp
   public :: horiz_interpolation
   public :: horiz_interpolation2
   public :: horiz_interpolation_wrfbdy
   public :: vert_interpolation
   public :: vert_interpolation_wrfbdy
   public :: decouple_wrfbdy
   public :: calc_delta_time_wrfbdy

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  build_casa_filename
   !
   ! DESCRIPTION:  Function for building a casa (GEOS-5 netCDF4) filename
   ! based on requested WRF date/time, file directory, and file prefix.
   !
   !---------------------------------------------------------------------------

   function build_casa_filename(wrf_datetime,casa_dir,casa_prefix) &
        result(casa_filename)

      ! Arguments
      character(len=MAX_DATE_LEN), intent(in) :: wrf_datetime
      character(len=*), intent(in) :: casa_dir
      character(len=*), intent(in) :: casa_prefix

      ! Return variable
      character(len=256) :: casa_filename

      ! Local variables
      character(len=4) :: cyear
      character(len=2) :: cmonth,cday,chour,cminute,csecond
      
      ! Build casa filename 
      read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear,cmonth,cday,chour,cminute,csecond

      !casa_filename = trim(casa_dir)//"/"//trim(casa_prefix)//"." &
      !     //cyear//cmonth//cday//'_'//chour//".nc"
      casa_filename = trim(casa_dir)//"/"//trim(casa_prefix)// &
           wrf_datetime//".nc"

      return
   end function build_casa_filename

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  build_casa_emission_filename
   !
   ! DESCRIPTION:  Function for building a casa (GEOS-5 netCDF4) filename
   ! based on requested WRF date/time, file directory, and file prefix.
   !
   !---------------------------------------------------------------------------

   function build_casa_emission_filename(wrf_datetime,casa_dir,casa_prefix) &
        result(casa_filename)

      ! Arguments
      character(len=MAX_DATE_LEN), intent(in) :: wrf_datetime
      character(len=*), intent(in) :: casa_dir
      character(len=*), intent(in) :: casa_prefix

      ! Return variable
      character(len=256) :: casa_filename

      ! Local variables
      character(len=4) :: cyear
      character(len=2) :: cmonth,cday,chour,cminute,csecond
      
      ! Build casa filename 
      read(wrf_datetime,"(A4,X,A2,X,A2,x,A2,x,A2,X,A2)") &                    
           cyear,cmonth,cday,chour,cminute,csecond
      !casa_filename = trim(casa_dir)//"/"//trim(casa_prefix)//"." &
      !     //cyear//cmonth//cday//'_'//chour//".nc"
      casa_filename = trim(casa_dir)//"/"//trim(casa_prefix)// &
           cyear//".nc"

      return
   end function build_casa_emission_filename


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
   ! ROUTINE:  calc_pressure_casa
   !
   ! DESCRIPTION:  Private method for calculating 3D pressure field from 
   ! GEOS5 data.
   !
   !---------------------------------------------------------------------------

   subroutine calc_pressure_casa(casa_lon_dim,casa_lat_dim, &
        casa_lev_dim,casa_time_dim,casa_ps,casa_pressure)
      
      ! Arguments
      integer, intent(in) :: casa_lon_dim
      integer, intent(in) :: casa_lat_dim 
      integer, intent(in) :: casa_lev_dim
      integer, intent(in) :: casa_time_dim
      real,intent(in) ::   casa_ps(casa_lon_dim,casa_lat_dim, &
           casa_time_dim)
      !real,intent(in) :: casa_delp(casa_lon_dim,casa_lat_dim, &
      !     casa_lev_dim,casa_time_dim)
      real,intent(inout) :: casa_pressure(casa_lon_dim,casa_lat_dim, &
           casa_lev_dim,casa_time_dim)

      ! Local variables
      real :: p_bottom, p_top
      integer :: ii,jj,kk,nn
      real ap1(1:casa_lev_dim+1),bp1(1:casa_lev_dim+1)
      real ap(1:casa_lev_dim+1),bp(1:casa_lev_dim+1)

      ap1=(/ &
        0.100000E-01, 0.208497E+01, 0.262021E+01, 0.327643E+01, &
        0.407657E+01, 0.504680E+01, 0.621680E+01, 0.761984E+01, &
        0.929294E+01, 0.112769E+02, 0.136434E+02, 0.164571E+02, &
        0.197916E+02, 0.237304E+02, 0.283678E+02, 0.338100E+02, &
        0.401754E+02, 0.476439E+02, 0.563879E+02, 0.666034E+02, &
        0.785123E+02, 0.923657E+02, 0.108663E+03, 0.127837E+03, &
        0.150393E+03, 0.176930E+03, 0.201192E+03, 0.216865E+03, &
        0.224363E+03, 0.223898E+03, 0.218776E+03, 0.212150E+03, &
        0.203259E+03, 0.193097E+03, 0.181619E+03, 0.169609E+03, &
        0.156260E+03, 0.142910E+03, 0.128696E+03, 0.118959E+03, &
        0.109182E+03, 0.993652E+02, 0.890999E+02, 0.788342E+02, &
        0.706220E+02, 0.643626E+02, 0.580532E+02, 0.516961E+02, &
        0.453390E+02, 0.389820E+02, 0.325708E+02, 0.260920E+02, &
        0.196131E+02, 0.131348E+02, 0.659375E+01, 0.480483E-01, &
        0.000000E+00/)

      bp1=(/  &
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, & 
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
        0.000000E+00, 0.817541E-08, 0.696003E-02, 0.280100E-01, &
        0.637201E-01, 0.113602E+00, 0.156224E+00, 0.200350E+00, &
        0.246741E+00, 0.294403E+00, 0.343381E+00, 0.392891E+00, &
        0.443740E+00, 0.494590E+00, 0.546304E+00, 0.581042E+00, &
        0.615818E+00, 0.650635E+00, 0.685900E+00, 0.721166E+00, &
        0.749378E+00, 0.770638E+00, 0.791947E+00, 0.813304E+00, &
        0.834661E+00, 0.856018E+00, 0.877429E+00, 0.898908E+00, &
        0.920387E+00, 0.941865E+00, 0.963406E+00, 0.984952E+00, &
        0.100000E+01 /)

! Flip ap and bp
       nn=casa_lev_dim+1
      do kk =1, casa_lev_dim+1
         ap(kk)=ap1(nn)
         bp(kk)=bp1(nn)
         nn=nn-1
      enddo

      do nn = 1, casa_time_dim
         do jj = 1, casa_lat_dim
            do ii = 1, casa_lon_dim
               p_bottom = casa_ps(ii,jj,nn)
               do kk = 1, casa_lev_dim
                  p_top = ap(kk+1) + bp(kk+1) * casa_ps(ii,jj,nn)
                  !p_top = ap(kk+1) + bp(kk+1) * p_bottom
                  casa_pressure(ii,jj,kk,nn) = (p_bottom + p_top) * 0.5
                  p_bottom = p_top
               end do
            end do
         end do
      end do
      !print *, 'CASA pressure', casa_pressure(120,120,:,1)
   end subroutine calc_pressure_casa

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  find_casa_latlon_indices
   !
   ! DESCRIPTION:  Given arrays of GOCART latitudes and longitudes, find
   ! array indices matching the southwest corner of four points surrounding
   ! the WRF grid column latitude/longitude.
   !
   !---------------------------------------------------------------------------

   subroutine find_casa_latlon_indices(casa_lon_dim, casa_lat_dim, &
        casa_longitudes, casa_latitudes, wrf_longitude, wrf_latitude, &
        i_casa_longitude_w, j_casa_latitude_s)
      
      ! Arguments
      integer, intent(in) :: casa_lon_dim
      integer, intent(in) :: casa_lat_dim
      real,intent(in) :: casa_longitudes(casa_lon_dim)
      real,intent(in) :: casa_latitudes(casa_lat_dim)
      real,intent(in) :: wrf_longitude
      real,intent(in) :: wrf_latitude
      integer,intent(inout) :: i_casa_longitude_w
      integer,intent(inout) :: j_casa_latitude_s

      ! Local variables
      integer :: ii,jj

      i_casa_longitude_w = casa_lon_dim
      do ii = 1, casa_lon_dim-1
         if (casa_longitudes(ii  ) .le. wrf_longitude .and. &
             casa_longitudes(ii+1) .gt. wrf_longitude) then
            i_casa_longitude_w = ii
         end if
      end do

      j_casa_latitude_s = casa_lat_dim
      do jj = 1, casa_lat_dim-1
         if (casa_latitudes(jj  ) .le. wrf_latitude .and. &
             casa_latitudes(jj+1) .gt. wrf_latitude) then
            j_casa_latitude_s = jj
         end if
      end do
      
   end subroutine find_casa_latlon_indices

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
   ! ROUTINE:  emission_interp
   !
   ! Emission is interpolated as the nearest grid point value. 
   ! ---------------------------------------------------------------------------

   function emission_interp(i,j,nn1,nn2,rlon,rlat, &
         nlon, nlat, elon, elat,emiss) result(TX)
!   
   implicit none
   integer i,j ! WRF grid point emission is interpolated to
   integer ic,jc ! nearest Emission data point.
   integer nn1,nn2 !WRF grid number of grid points
   real, dimension(nn1,nn2) :: rlat,rlon !wrf lat & lon
   integer nlat, nlon!Emission grid number of grid points
   real, dimension(nlon) :: elon !emission data grid lon
   real, dimension(nlat) :: elat !emission data grid  lat
   real ilatn,ilonn ! emission data resolution 0.5x0.5 deg etc
   real, dimension(nlon,nlat) :: emiss
   real :: TX !New interpolated emission

   !-local var
   real dlonr,dlatr !WRF grid dx and dy for the (i,j) grid
   real :: rrlat, rrlon,difflon
   integer qi1,qi2,qj1,qj2,ncount,ii,jj
    ilonn = elon(3)- elon(2)
    ilatn = elat(3) - elat(2)
    !print *, 'elat: ',elat
    !print *, 'ilonn, ilatn: ',ilonn, ilatn

   ! Compute the nearest emission grid point (ic, jc)
    rrlat=rlat(i,j)
    rrlon=rlon(i,j) 

    if(rrlon.lt.0)rrlon=360.+rrlon    
    difflon=(rrlon-elon(1))
    if(difflon.lt.0.)difflon=abs(difflon)
    if(difflon.ge.360.)difflon=difflon-360.
    !print *, 'difflon:',difflon
   ic =   (nint(difflon/ilonn)) + 1       
   jc =   (nint((rrlat-elat(1))/ilatn)) + 1

   dlonr=0.5*(rlon(nn1,j)-rlon(1,j))/float(nn1-1)
   dlatr=0.5*(rlat(i,nn2)-rlat(i,1))/float(nn2-1)

   qi1=int(dlonr/ilonn+0.5)
   qi2=int(dlonr/ilonn+0.5)
   qj1=int(dlatr/ilatn+0.5)
   qj2=int(dlatr/ilatn+0.5)

   ncount = 0
   TX  = 0.
!
!Emission is emiss(ic,jc) if WRF grid < emission grid
!
   do jj =max(1,jc-qj1),min(nlat,jc+qj2)
      do ii = max(1,ic-qi1),min(nlon,ic+qi2)
         ncount = ncount + 1
         TX= TX + emiss(ii,jj)
      enddo
   enddo
   TX = TX / (float(ncount) + 1.E-10) ! interpolated rate

   end function emission_interp 

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  horiz_interpolation
   !
   ! DESCRIPTION:  Oversees horizontal interpolation of CASA` data to
   ! WRF grid.
   !
   ! NOTE:  Interpolation is actually done via a lower-level bilinear
   ! interpolation scheme. 
   !
   !---------------------------------------------------------------------------

   subroutine horiz_interpolation(wrf_dim_1, wrf_dim_2, &
        wrf_latitudes, wrf_longitudes, &
        casa_lon_dim,casa_lat_dim, &
        casa_longitudes, casa_latitudes, &
        casa_lev_dim, casa_variable, &
        inverse_casa_delta_lon, inverse_casa_delta_lat, &
        casa_interp_variable)

      ! Arguments
      integer,intent(in) :: wrf_dim_1
      integer,intent(in) :: wrf_dim_2
      real,intent(in) :: wrf_latitudes(wrf_dim_1,wrf_dim_2)
      real,intent(in) :: wrf_longitudes(wrf_dim_1,wrf_dim_2)
      integer,intent(in) :: casa_lon_dim
      integer,intent(in) :: casa_lat_dim
      real,intent(in) :: casa_longitudes(casa_lon_dim)
      real,intent(in) :: casa_latitudes(casa_lat_dim)
      integer,intent(in) :: casa_lev_dim
      real,intent(in) :: casa_variable(casa_lon_dim,casa_lat_dim, &
           casa_lev_dim)
      real,intent(in) :: inverse_casa_delta_lon
      real,intent(in) :: inverse_casa_delta_lat
      real,intent(inout) :: casa_interp_variable(wrf_dim_1,wrf_dim_2, &
           casa_lev_dim)
      
      ! Local variables
      real :: wrf_latitude, wrf_longitude
      real :: casa_longitude_w, casa_latitude_s
      real :: q_sw,q_se,q_nw,q_ne,q_interp
      integer :: i_casa_longitude_w, j_casa_latitude_s
      integer :: ii_w, ii_e, jj_s, jj_n
      integer :: ii,jj,kk
      do jj = 1, wrf_dim_2
         do ii = 1, wrf_dim_1
            wrf_latitude = wrf_latitudes(ii,jj)
            wrf_longitude = wrf_longitudes(ii,jj)

            call find_casa_latlon_indices(casa_lon_dim, &
                 casa_lat_dim, casa_longitudes, casa_latitudes, &
                 wrf_longitude, wrf_latitude, &
                 i_casa_longitude_w, j_casa_latitude_s)

            ii_w = i_casa_longitude_w
            casa_longitude_w = casa_longitudes(ii_w)
            if (i_casa_longitude_w == casa_lon_dim) then
               ii_e = 1
            else
               ii_e = ii_w + 1
            end if

            jj_s = j_casa_latitude_s
            jj_n = jj_s + 1
            casa_latitude_s = casa_latitudes(jj_s)

            ! Loop through each casa level and interpolate pressure
            ! to WRF grid
            do kk = 1, casa_lev_dim
               q_sw = casa_variable(ii_w,jj_s,kk)
               q_se = casa_variable(ii_e,jj_s,kk)
               q_nw = casa_variable(ii_w,jj_n,kk)
               q_ne = casa_variable(ii_e,jj_n,kk)

               q_interp = bilinear_interp(casa_longitude_w, &
                    inverse_casa_delta_lon, &
                    wrf_longitude, &
                    casa_latitude_s, &
                    inverse_casa_delta_lat,wrf_latitude, &
                    q_sw,q_se,q_nw,q_ne)

               casa_interp_variable(ii,jj,kk) = q_interp
            end do ! kk
         end do ! jj
      end do ! ii
      
   end subroutine horiz_interpolation
   !---------------------------------------------------------------------------
   !
   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  horiz_interpolation2
   !
   ! DESCRIPTION:  Oversees horizontal interpolation of CASA` data to
   ! WRF grid. (another method for CO2 tracers) Added by Jossy Jacob
   !
   ! NOTE: Interpolation is done by averaging the grid points around, 
   ! effort to conserving mass of tracer
   !---------------------------------------------------------------------------

   subroutine horiz_interpolation2(wrf_dim_1, wrf_dim_2, &
        wrf_latitudes, wrf_longitudes, &
        casa_lon_dim,casa_lat_dim, &
        casa_longitudes, casa_latitudes, &
        casa_lev_dim, casa_variable, &
        casa_interp_variable)

      ! Arguments
      integer,intent(in) :: wrf_dim_1
      integer,intent(in) :: wrf_dim_2
      real,intent(in) :: wrf_latitudes(wrf_dim_1,wrf_dim_2)
      real,intent(in) :: wrf_longitudes(wrf_dim_1,wrf_dim_2)
      integer,intent(in) :: casa_lon_dim
      integer,intent(in) :: casa_lat_dim
      real,intent(in) :: casa_longitudes(casa_lon_dim)
      real,intent(in) :: casa_latitudes(casa_lat_dim)
      integer,intent(in) :: casa_lev_dim
      real,intent(in) :: casa_variable(casa_lon_dim,casa_lat_dim, &
           casa_lev_dim)
      real,intent(inout) :: casa_interp_variable(wrf_dim_1,wrf_dim_2, &
           casa_lev_dim)
      
      ! Local variables
      real ::q_interp
      integer :: ii,jj,kk
      do jj = 1, wrf_dim_2
         do ii = 1, wrf_dim_1          
            ! Loop through each casa level and interpolate emission 
            ! to WRF grid
            do kk = 1, casa_lev_dim

               q_interp = emission_interp(ii,jj,wrf_dim_1, wrf_dim_2, &
                  wrf_longitudes, wrf_latitudes, &
                 casa_lon_dim,casa_lat_dim,casa_longitudes, casa_latitudes, &
                 casa_variable(:,:,kk)) 

               casa_interp_variable(ii,jj,kk) = q_interp
            end do ! kk
         end do ! jj
      end do ! ii
      
   end subroutine horiz_interpolation2
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
        casa_lon_dim, casa_lat_dim, &  
        casa_longitudes, casa_latitudes, &
        casa_lev_dim, casa_variable, &
        inverse_casa_delta_lon, inverse_casa_delta_lat, &
        casa_interp_variable)

      ! Arguments
      integer,intent(in) :: wrfbdy_dim_1
      integer,intent(in) :: wrfbdy_w_dim
      real,intent(in) :: wrfbdy_latitudes(wrfbdy_dim_1,wrfbdy_w_dim)
      real,intent(in) :: wrfbdy_longitudes(wrfbdy_dim_1,wrfbdy_w_dim)
      integer,intent(in) :: casa_lon_dim
      integer,intent(in) :: casa_lat_dim
      real,intent(in) :: casa_longitudes(casa_lon_dim)
      real,intent(in) :: casa_latitudes(casa_lat_dim)
      integer,intent(in) :: casa_lev_dim
      real,intent(in) :: casa_variable(casa_lon_dim,casa_lat_dim, &
           casa_lev_dim)
      real,intent(in) :: inverse_casa_delta_lon
      real,intent(in) :: inverse_casa_delta_lat
      real,intent(inout) :: casa_interp_variable(wrfbdy_dim_1, &
           casa_lev_dim, wrfbdy_w_dim)

      ! Local variables
      real :: wrfbdy_latitude, wrfbdy_longitude
      real :: casa_longitude_w, casa_latitude_s
      real :: q_sw,q_se,q_nw,q_ne,q_interp
      integer :: i_casa_longitude_w, j_casa_latitude_s
      integer :: ii_w, ii_e, jj_s, jj_n
      integer :: ii,jj,kk

      do jj = 1, wrfbdy_w_dim
         do ii = 1, wrfbdy_dim_1
            wrfbdy_latitude = wrfbdy_latitudes(ii,jj)
            wrfbdy_longitude = wrfbdy_longitudes(ii,jj)

            call find_casa_latlon_indices(casa_lon_dim, &
                 casa_lat_dim, casa_longitudes, casa_latitudes, &
                 wrfbdy_longitude, wrfbdy_latitude, &
                 i_casa_longitude_w, j_casa_latitude_s)

            ii_w = i_casa_longitude_w
            casa_longitude_w = casa_longitudes(ii_w)
            if (i_casa_longitude_w == casa_lon_dim) then
               ii_e = 1
            else
               ii_e = ii_w + 1
            end if

            jj_s = j_casa_latitude_s
            jj_n = jj_s + 1
            casa_latitude_s = casa_latitudes(jj_s)
                  
            ! Loop through each casa level and interpolate pressure
            ! to WRF grid
            do kk = 1, casa_lev_dim
               q_sw = casa_variable(ii_w,jj_s,kk)
               q_se = casa_variable(ii_e,jj_s,kk)
               q_nw = casa_variable(ii_w,jj_n,kk)
               q_ne = casa_variable(ii_e,jj_n,kk)

               q_interp = bilinear_interp(casa_longitude_w, &
                    inverse_casa_delta_lon, &
                    wrfbdy_longitude, &
                    casa_latitude_s, &
                    inverse_casa_delta_lat,wrfbdy_latitude, &
                    q_sw,q_se,q_nw,q_ne)

               casa_interp_variable(ii,kk,jj) = q_interp
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
        casa_lev_dim,wrf_pressures,casa_horiz_interp_pressures, &
        casa_horiz_interp_var,casa_final_interp_var)

      ! Arguments
      integer, intent(in) :: wrf_x_dim
      integer, intent(in) :: wrf_y_dim
      integer, intent(in) :: wrf_z_dim
      integer, intent(in) :: wrf_t_dim
      integer, intent(in) :: casa_lev_dim
      real,intent(in) :: wrf_pressures(wrf_x_dim,wrf_y_dim,wrf_z_dim,wrf_t_dim)
      real,intent(in) :: casa_horiz_interp_pressures(wrf_x_dim,wrf_y_dim, &
           casa_lev_dim,wrf_t_dim)
      real,intent(in) :: casa_horiz_interp_var(wrf_x_dim,wrf_y_dim, &
           casa_lev_dim,wrf_t_dim)
      real,intent(out) :: casa_final_interp_var(wrf_x_dim,wrf_y_dim, &
           wrf_z_dim,wrf_t_dim)

      ! Local variables
      real :: casa_p_base, casa_log_p_base
      real :: casa_p_top, casa_log_p_top
      real :: casa_var_base, casa_var_top, casa_var_interp
      real :: wrf_p, wrf_log_p
      logical :: done
      integer :: ii,jj,kk_wrf,kk_casa,nn

      do nn = 1, wrf_t_dim
         do jj = 1, wrf_y_dim
            do ii = 1, wrf_x_dim

               do kk_wrf = 1, wrf_z_dim

                  wrf_p = wrf_pressures(ii,jj,kk_wrf,nn)

                  ! Special case:  WRF level is physically below GOCART
                  ! data.  Just use lowest GOCART level.
                  casa_p_base = casa_horiz_interp_pressures(ii,jj,1,nn)

                  if (wrf_p .ge. casa_p_base) then
                     casa_final_interp_var(ii,jj,kk_wrf,nn) = &
                          casa_horiz_interp_var(ii,jj,1,nn)
                     cycle
                  end if

                  ! Special case:  WRF level is physically above GOCART
                  ! data.  Just use highest GOCART level.
                  casa_p_top = &
                       casa_horiz_interp_pressures(ii,jj,casa_lev_dim,nn)

                  if (wrf_p .le. casa_p_top) then
                     casa_final_interp_var(ii,jj,kk_wrf,nn) = &
                          casa_horiz_interp_var(ii,jj,casa_lev_dim,&
                          nn)
                     cycle
                  end if

                  ! General case:  Linearly interpolate in log(p).
                  done = .false.
                  do kk_casa = 1, casa_lev_dim-1
                     casa_p_base = &
                          casa_horiz_interp_pressures(ii,jj,kk_casa,nn)
                     casa_p_top = &
                          casa_horiz_interp_pressures(ii,jj,kk_casa+1,nn)

                     if (casa_p_base .ge. wrf_p .and. &
                          casa_p_top .lt. wrf_p) then
                        done = .true.
                        exit
                     end if
                  end do
                  if (.not. done) then
                     print*,'ERROR, cannot find wrf level in casa!'
                     print*,'casa_p_base = ',casa_p_base
                     print*,'casa_p_top = ',casa_p_top
                     print*,'wrf_p = ',wrf_p
                     print*,'ii,jj,kk_wrf,nn = ',ii,jj,kk_wrf,nn
                     stop
                  end if
                  
                  casa_log_p_base = log(casa_p_base)
                  casa_log_p_top = log(casa_p_top)
                  wrf_log_p = log(wrf_p)

                  casa_var_base = &
                       casa_horiz_interp_var(ii,jj,kk_casa,nn)
                  casa_var_top = &
                       casa_horiz_interp_var(ii,jj,kk_casa+1,nn)
                  casa_var_interp = &
                       casa_var_base + &
                       ( (casa_var_top - casa_var_base) * &
                         (wrf_log_p - casa_log_p_base) / &
                         (casa_log_p_top - casa_log_p_base) )

                  casa_final_interp_var(ii,jj,kk_wrf,nn) = casa_var_interp

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
        wrfbdy_w_dim,wrfbdy_t_dim,casa_lev_dim,wrf_pressures, &
        casa_horiz_interp_pressures, &
        casa_horiz_interp_var,casa_final_interp_var)

      ! Arguments
      integer, intent(in) :: wrfbdy_dim_1
      integer, intent(in) :: wrfbdy_z_dim
      integer, intent(in) :: wrfbdy_w_dim
      integer, intent(in) :: wrfbdy_t_dim
      integer, intent(in) :: casa_lev_dim
      real,intent(in) :: wrf_pressures(wrfbdy_dim_1, &
           wrfbdy_z_dim,wrfbdy_w_dim,wrfbdy_t_dim)
      real,intent(in) :: casa_horiz_interp_pressures(wrfbdy_dim_1, &
           casa_lev_dim,wrfbdy_w_dim,wrfbdy_t_dim)
      real,intent(in) :: casa_horiz_interp_var(wrfbdy_dim_1, &
           casa_lev_dim,wrfbdy_w_dim,wrfbdy_t_dim)
      real,intent(out) :: casa_final_interp_var(wrfbdy_dim_1, &
           wrfbdy_z_dim, wrfbdy_w_dim,wrfbdy_t_dim)

      ! Local variables
      real :: casa_p_base, casa_log_p_base
      real :: casa_p_top, casa_log_p_top
      real :: casa_var_base, casa_var_top, casa_var_interp
      real :: wrf_p, wrf_log_p
      logical :: done
      integer :: ii,jj_wrf,jj_casa,kk,nn

      do nn = 1, wrfbdy_t_dim
         do kk = 1, wrfbdy_w_dim
            do ii = 1, wrfbdy_dim_1

               do jj_wrf = 1, wrfbdy_z_dim

                  wrf_p = wrf_pressures(ii,jj_wrf,kk,nn)

                  ! Special case:  WRF level is physically below GOCART
                  ! data.  Just use lowest GOCART level.
                  casa_p_base = casa_horiz_interp_pressures(ii,1,kk,nn)
                  if (wrf_p .ge. casa_p_base) then
                     casa_final_interp_var(ii,jj_wrf,kk,nn) = &
                          casa_horiz_interp_var(ii,1,kk,nn)
                     cycle
                  end if

                  ! Special case:  WRF level is physically above GOCART
                  ! data.  Just use highest GOCART level.
                  casa_p_top = &
                       casa_horiz_interp_pressures(ii,casa_lev_dim,kk,nn)
                  if (wrf_p .le. casa_p_top) then
                     casa_final_interp_var(ii,jj_wrf,kk,nn) = &
                          casa_horiz_interp_var(ii,casa_lev_dim,kk,&
                          nn)
                     cycle
                  end if

                  ! General case:  Linearly interpolate in log(p).
                  done = .false.
                  do jj_casa = 1, casa_lev_dim-1
                     casa_p_base = &
                          casa_horiz_interp_pressures(ii,jj_casa  ,kk,nn)
                     casa_p_top = &
                          casa_horiz_interp_pressures(ii,jj_casa+1,kk,nn)

                     if (casa_p_base .ge. wrf_p .and. &
                          casa_p_top .lt. wrf_p) then
                        done = .true.
                        exit
                     end if
                  end do
                  if (.not. done) then
                     print*,'ERROR, cannot find wrf level in casa!'
                     print*,'casa_p_base = ',casa_p_base
                     print*,'casa_p_top = ',casa_p_top
                     print*,'wrf_p = ',wrf_p
                     print*,'ii,jj_wrf,kk,nn = ',ii,jj_wrf,kk,nn
                     stop
                  end if
                  
                  casa_log_p_base = log(casa_p_base)
                  casa_log_p_top = log(casa_p_top)
                  wrf_log_p = log(wrf_p)

                  casa_var_base = &
                       casa_horiz_interp_var(ii,jj_casa  ,kk,nn)
                  casa_var_top = &
                       casa_horiz_interp_var(ii,jj_casa+1,kk,nn)
                  casa_var_interp = &
                       casa_var_base + &
                       ( (casa_var_top - casa_var_base) * &
                         (wrf_log_p - casa_log_p_base) / &
                         (casa_log_p_top - casa_log_p_base) )

                  casa_final_interp_var(ii,jj_wrf,kk,nn) = casa_var_interp

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
