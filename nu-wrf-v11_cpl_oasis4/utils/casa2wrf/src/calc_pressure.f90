
program calc_pressure

   implicit none

   include "netcdf.inc"

   ! wrfinput dimensions

   ! wrfinput arrays
   ! Perturbation geopotential (m^2 s^-2)
   real,allocatable :: PH(:,:,:,:)
   ! Base state geopotential (m^2 s^-2)
   real,allocatable :: PHB(:,:,:,:)
   ! Perturbation potential temperature (K)
   real,allocatable :: T(:,:,:,:) 
   ! Perturbation dry air mass in column (Pa)
   real,allocatable :: MU(:,:,:) 
   ! Base-state dry air mass in column (Pa)
   real,allocatable :: MUB(:,:,:)
   ! Inverse d(eta) values between full (w) levels
   real,allocatable :: RDNW(:,:)
   ! Water vapor mixing ratio (kg/kg)
   real,allocatable :: QVAPOR(:,:,:,:)

   ! Derived variables
   real,allocatable :: geopotential(:,:,:,:)
   real,allocatable :: theta(:,:,:,:)
   real,allocatable :: theta_moist(:,:,:,:)
   real,allocatable :: mu_dry(:,:,:)
   real,allocatable :: alpha_dry(:,:,:,:)
   real,allocatable :: pressure(:,:,:,:)

   ! Dimensions
   integer :: Time, bottom_top, bottom_top_stag, south_north, west_east

   ! NetCDF variables
   integer :: ncid
   integer :: dimid
   integer :: varid
   integer, dimension(2) :: start2d
   integer, dimension(3) :: start3d
   integer, dimension(4) :: start4d
   integer, dimension(2) :: count2d
   integer, dimension(3) :: count3d
   integer, dimension(4) :: count4d
   data start2d /1, 1/
   data start3d /1, 1, 1/
   data start4d /1, 1, 1, 1/

   ! Other variables
   character(len=132) :: inputfile
   integer :: status
   integer :: i,j,k,n

   ! Constants
   real,parameter :: P_0 = 100000. ! Reference pressure (Pa)
   real,parameter :: R_D = 287.    ! Dry air gas constant
   real,parameter :: C_P = 7.*R_D*0.5
   real,parameter :: C_V = C_P - R_D
   real,parameter :: GAMMA = C_P/C_V
   real,parameter :: T0 = 300.
   real,parameter :: R_V = 461.6
   real,parameter :: RV_OVR_RD = R_V/R_D

   ! Beginning of executable code
   inputfile = trim("wrfinput_d01")

   ! Open file
   status = nf_open(trim(inputfile),NF_NOWRITE,ncid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if

   ! Get Time dimension
   status = nf_inq_dimid(ncid,"Time",dimid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if
   status = nf_inq_dimlen(ncid,dimid,Time)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if

   ! Get west_east dimension
   status = nf_inq_dimid(ncid,"west_east",dimid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if
   status = nf_inq_dimlen(ncid,dimid,west_east)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if

   ! Get south_north dimension
   status = nf_inq_dimid(ncid,"south_north",dimid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if
   status = nf_inq_dimlen(ncid,dimid,south_north)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if

   ! Get bottom_top dimension
   status = nf_inq_dimid(ncid,"bottom_top",dimid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if
   status = nf_inq_dimlen(ncid,dimid,bottom_top)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if

   ! Get bottom_top_stag dimension
   status = nf_inq_dimid(ncid,"bottom_top_stag",dimid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if
   status = nf_inq_dimlen(ncid,dimid,bottom_top_stag)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if

   ! Read ph
   status = nf_inq_varid(ncid,"PH",varid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if
   allocate(ph(west_east,south_north,bottom_top_stag,Time))
   count4d = (/ west_east, south_north, bottom_top_stag, Time /)
   status = nf_get_vara_real(ncid,varid,start4d,count4d,ph)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if
   
   ! Read phb
   status = nf_inq_varid(ncid,"PHB",varid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if
   allocate(phb(west_east,south_north,bottom_top_stag,Time))
   count4d = (/ west_east, south_north, bottom_top_stag, Time /)
   status = nf_get_vara_real(ncid,varid,start4d,count4d,phb)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if

   ! Calculate geopotential
   allocate(geopotential(west_east,south_north,bottom_top_stag,Time))
   do n = 1, Time
      do k = 1, bottom_top_stag
         do j = 1, south_north
            do i = 1, west_east
               geopotential(i,j,k,n) = ph(i,j,k,n) + phb(i,j,k,n)
            end do
         end do
      end do
   end do
   deallocate(ph)
   deallocate(phb)

   ! Read t
   status = nf_inq_varid(ncid,"T",varid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if
   allocate(t(west_east,south_north,bottom_top,Time))
   count4d = (/ west_east, south_north, bottom_top, Time /)
   status = nf_get_vara_real(ncid,varid,start4d,count4d,t)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if

   ! Calculate theta
   allocate(theta(west_east,south_north,bottom_top,Time))
   do n = 1, Time
      do k = 1, bottom_top
         do j = 1, south_north
            do i = 1, west_east
               theta(i,j,k,n) = t(i,j,k,n) + T0
            end do
         end do
      end do
   end do
   deallocate(t)

   ! Read qvapor
   status = nf_inq_varid(ncid,"QVAPOR",varid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if
   allocate(qvapor(west_east,south_north,bottom_top,Time))
   count4d = (/ west_east, south_north, bottom_top, Time /)
   status = nf_get_vara_real(ncid,varid,start4d,count4d,qvapor)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if

   ! Calculate moist potential temperature
   allocate(theta_moist(west_east,south_north,bottom_top,Time))
   do n = 1, Time
      do k = 1, bottom_top
         do j = 1, south_north
            do i = 1, west_east
               theta_moist(i,j,k,n) = &
                    theta(i,j,k,n)*(1. + RV_OVR_RD*qvapor(i,j,k,n))
            end do
         end do
      end do
   end do
   deallocate(theta)
   deallocate(qvapor)

   ! Read mu
   status = nf_inq_varid(ncid,"MU",varid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if
   allocate(mu(west_east,south_north,Time))
   count3d = (/ west_east, south_north, Time /)
   status = nf_get_vara_real(ncid,varid,start3d,count3d,mu)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if

   ! Read mub
   status = nf_inq_varid(ncid,"MUB",varid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if
   allocate(mub(west_east,south_north,Time))
   count3d = (/ west_east, south_north, Time /)
   status = nf_get_vara_real(ncid,varid,start3d,count3d,mub)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if

   ! Calculate total dry air mass
   allocate(mu_dry(west_east, south_north, Time))
   do n = 1, Time
      do j = 1, south_north
         do i = 1, west_east
            mu_dry(i,j,n) = mu(i,j,n) + mub(i,j,n)
         end do
      end do
   end do
   deallocate(mu)
   deallocate(mub)

   ! Read rdnw
   status = nf_inq_varid(ncid,"RDNW",varid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if
   allocate(rdnw(bottom_top,Time))
   count2d = (/ bottom_top, Time /)
   status = nf_get_vara_real(ncid,varid,start2d,count2d,rdnw)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
   end if

   ! Close file
   status = nf_close(ncid)
   if (status .ne. NF_NOERR) then
      print*, trim(nf_strerror(status))
      stop
   end if

   ! Calculate inverse of dry air density
   allocate(alpha_dry(west_east,south_north,bottom_top,Time))
   do n = 1, Time
      do k = 1, bottom_top
         do j = 1, south_north
            do i = 1, west_east
               alpha_dry(i,j,k,n) = &
                    -1 * &
                    (geopotential(i,j,k+1,n) - geopotential(i,j,k,n)) * &
                        rdnw(k,n) / mu_dry(i,j,n)
            end do
         end do
      end do
   end do
   deallocate(geopotential)
   deallocate(rdnw)
   deallocate(mu_dry)

   ! Calculate full pressure
   allocate(pressure(west_east,south_north,bottom_top,Time))
   do n = 1, Time
      do k = 1, bottom_top
         do j = 1, south_north
            do i = 1, west_east
               pressure(i,j,k,n) = &
                    R_D*theta_moist(i,j,k,n)/P_0/alpha_dry(i,j,k,n)
               pressure(i,j,k,n) = (pressure(i,j,k,n))**GAMMA
               pressure(i,j,k,n) = P_0*pressure(i,j,k,n)
            end do
         end do
      end do
   end do
   deallocate(theta_moist)
   deallocate(alpha_dry)

   print*,'pressure(1,1,1,1) = ',pressure(1,1,1,1)


   deallocate(pressure)
end program calc_pressure
