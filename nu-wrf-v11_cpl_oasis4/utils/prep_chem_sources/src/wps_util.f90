!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE: wps_util
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! Contains utility routines to call WPS map projection from prep_chem_sources.
!
!------------------------------------------------------------------------------

module wps_util

   ! Change defaults
   implicit none
   private

   ! Public variables
   public :: wps_latlon
   type wps_latlon
      real,allocatable :: lat_m(:,:) ! Latitude at mass points
      real,allocatable :: lon_m(:,:) ! Longitude at mass points
      real :: swlat
      real :: swlon
      real :: nelat
      real :: nelon
   end type wps_latlon

   public :: wps_grids
   type(wps_latlon),allocatable :: wps_grids(:)

   ! Public routines
   public :: wps_calc_latlon

contains

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  wps_calc_latlon
   !
   ! DESCRIPTION:  Calculates lat/lon of grid using WPS map projection code.
   ! Based on logic in geogrid.exe.
   !
   !--------------------------------------------------------------------------

   subroutine wps_calc_latlon()

      ! Modules
      use mem_grid  ! from RAMS
      use llxy_module ! from WPS 

      ! Local variables
      integer :: start_mem_i, end_mem_i, start_mem_j, end_mem_j
      real :: rx,ry
      logical :: ew_extra_col, sn_extra_row
      integer :: igrid
      integer :: i,j

      ! WPS code is designed to optionally support parallel runs.  We
      ! don't actually use that here, but it's good to initialize the
      ! code, as is done in geogrid.exe.
      call parallel_start()

      ! Assemble the map projection input.  Variables on the left-hand-side
      ! are generally from llxy_module, while those on the right are generally
      ! from mem_grid.
      known_lat = centlat(1)
      known_lon = centlon(1)
      known_x = real(nnxp(1)) / 2.
      known_y = real(nnyp(1)) / 2.
      dxkm = deltaxn(1) ! Despite the name, dxkm expects distance in meters.
      dykm = deltayn(1) ! Despite the name, dykm expects distance in meters.
      truelat1 = stdlat1
      truelat2 = stdlat2
      stand_lon = stdlon

      ! Convert RAMS integer code to WPS integer code.
      iproj_type = get_wps_iproj(ihtran) 
      
      ! Loop through grids.
      grid_is_active(:) = .false.
      n_domains = ngrids
      do igrid = 1,ngrids
         ixdim(igrid) = nnxp(igrid)
         jydim(igrid) = nnyp(igrid)
         parent_grid_ratio(igrid) = nstratx(igrid)
         parent_id(igrid) = nxtnest(igrid)
         subgrid_ratio_x(igrid) = 1
         subgrid_ratio_y(igrid) = 1
         parent_ll_x(igrid) = real(ninest(igrid))
         parent_ll_y(igrid) = real(njnest(igrid))
         parent_ur_x(igrid) = &
              real(ninest(igrid))+real(ixdim(igrid)) / &
              real(nstratx(igrid))-1.
         parent_ur_y(igrid) = &
              real(njnest(igrid))+real(jydim(igrid)) / &
              real(nstraty(igrid))-1.
         grid_is_active(igrid) = .true.
      end do
      
      ! Compute parameters necessary for transformations to any nest.
      call compute_nest_locations()

      ! Process all requested domains.
      allocate(wps_grids(ngrids))
      do igrid = 1,ngrids
         if (.not. grid_is_active(igrid)) cycle

         ! Set transformations in llxy module to be with respect to current
         ! nest.
         call select_domain(igrid)

         ! We won't bother with parallel support, so just use the full
         ! dimensions.
         start_mem_i = 1
         start_mem_j = 1
         end_mem_i = ixdim(igrid)
         end_mem_j = jydim(igrid)

         allocate(wps_grids(igrid)%lat_m(end_mem_i,end_mem_j))
         allocate(wps_grids(igrid)%lon_m(end_mem_i,end_mem_j))

         ! Calculate lat/lon
         rx = 1.0
         ry = 1.0
         do i=start_mem_i, end_mem_i
            do j=start_mem_j, end_mem_j
               call xytoll(real(i-1)/rx+1.0, real(j-1)/ry+1.0, &
                    wps_grids(igrid)%lat_m(i,j), &
                    wps_grids(igrid)%lon_m(i,j), M)
           end do
         end do

         ! Calculate lat/lon of southwest corner
         i = start_mem_i
         j = start_mem_j
         call xytoll(real(i-1)/rx+0.5, real(j-1)/ry+0.5, &
              wps_grids(igrid)%swlat, &
              wps_grids(igrid)%swlon, M)

         ! Calculate lat/lon of northeast corner
         i = end_mem_i
         j = end_mem_j
         call xytoll(real(i-1)/rx+0.5, real(j-1)/ry+0.5, &
              wps_grids(igrid)%nelat, &
              wps_grids(igrid)%nelon, M)

      end do

      return
   end subroutine wps_calc_latlon

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  get_wps_iproj
   !
   ! DESCRIPTION:  Given RAMS map projection flag, return corresponding
   ! WPS flag.
   !
   !--------------------------------------------------------------------------

   function get_wps_iproj(iproj_rams) result(iproj_wps)
      
      ! Modules
      use map_utils

      ! Arguments
      integer,intent(in) :: iproj_rams

      ! Return variable
      integer :: iproj_wps

      ! Convert RAMS map projection flag to WPS
      if (iproj_rams == 1) then ! Polar stereographic
         iproj_wps = PROJ_PS
      else if (iproj_rams == 2) then ! Lambert conformal
         iproj_wps = PROJ_LC 
      else if (iproj_rams == 3) then ! Mercator
         iproj_wps = PROJ_MERC
      else
         print*,'ERROR, unknown map projection for WPS!'
         print*,'iproj_rams = ',iproj_rams
         stop
      end if

      return
   end function get_wps_iproj
   
end module wps_util
