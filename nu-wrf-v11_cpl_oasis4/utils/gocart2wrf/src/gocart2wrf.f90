!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  gocart2wrf
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Complete rewrite of original GOCART2WRF.  Reads in namelist.gocart2wrf,
! netCDF files containing GOCART aerosol data, and netCDF wrfbdy and wrfinput
! files output from WRF REAL program.  Interpolates user-selected GOCART
! variables to WRF grid with optional unit conversions, and adds data to
! wrfinput and wrfbdy files.
!
!------------------------------------------------------------------------------

program gocart2wrf

   ! Modules
   use gocart_variables_mod
   use interpolator_mod
   use sorted_datetimes_mod
   use nml_mod

   ! Force declarations
   implicit none

   ! Data types
   type(interpolator) :: interp_object

   ! Other variables
   character(len=132) :: wrf_dir
   character(len=132) :: gocart_dir
   character(len=132) :: gocart_prefix
   character(len=132) :: gocart_source
   integer :: num_wrf_domains
   real :: conversion_factor
   integer :: ii,jj
   
   ! Local constants                                                         
   character(len=*), parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
   character(len=*), parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

   type(nml) nml_g2w
   character(len=256) :: filename

   ! Read in namelist.gocart2wrf
   filename = 'namelist.gocart2wrf'
   nml_g2w = nml_new(filename)

   ! Get directory with wrf input and bdy files
   call nml_g2w%read("wrf", "wrf_dir", wrf_dir) 

   ! Get number of domains to process
   call nml_g2w%read("wrf", "max_dom", num_wrf_domains) 

   ! Get directory and prefix of GOCART files
   call nml_g2w%read("gocart_shared", "gocart_dir", gocart_dir)
   call nml_g2w%read("gocart_shared", "gocart_prefix", gocart_prefix)
   call nml_g2w%read("gocart_shared", "gocart_source", gocart_source, &
        init=.true.)

   ! Sanity check gocart_source
   do jj = 1,len(gocart_source)
      ii = index(LOWER_CASE,gocart_source(jj:jj))
      if (ii .ne. 0) then
         gocart_source(jj:jj) = UPPER_CASE(ii:ii)
      end if
   end do
   if (trim(gocart_source) .eq. "MISSING") then
      print*,'Resetting gocart_source to GEOS5'
      gocart_source = GEOS5
   end if

   ! Create interpolator object
   interp_object = create_interpolator(wrf_dir,num_wrf_domains,gocart_dir, &
        gocart_prefix,NUM_GOCART_VARIABLES,gocart_wrf_variable_names, &
        gocart_wrf_variable_units,gocart_wrf_variable_descriptions, &
        gocart_source)

   ! Loop through each species
   do ii = 1, NUM_GOCART_VARIABLES
      conversion_factor = gocart_geos5_conversion_factors(ii)
      if (trim(gocart_source) .eq. MERRA2) then
         print*,'Processing ',trim(gocart_merra2_variable_names(ii))

         ! Use same conversion factors as GEOS
         call run_interpolator(interp_object, &
              trim(gocart_merra2_variable_names(ii)), &
              trim(gocart_wrf_variable_names(ii)), &
              gocart_geos5_conversion_factors(ii))         
      else if (trim(gocart_source) .eq. OFFLINE) then
         print*,'Processing ',trim(gocart_offline_variable_names(ii))

         ! Use unique conversion factors for offline GOCART data
         call run_interpolator(interp_object, &
              trim(gocart_offline_variable_names(ii)), &
              trim(gocart_wrf_variable_names(ii)), &
              gocart_offline_conversion_factors(ii))         

      else if (trim(gocart_source) .eq. MERRAERO .or. &
           trim(gocart_source) .eq. GEOS5) then
         print*,'Processing ',trim(gocart_geos5_variable_names(ii))

         ! Use GEOS5 conversion factors
         call run_interpolator(interp_object, &
              trim(gocart_geos5_variable_names(ii)), &
              trim(gocart_wrf_variable_names(ii)), &
              gocart_geos5_conversion_factors(ii))
      else
         write(6,*),"ERROR, invalid selection for gocart_source!"
         write(6,*),'gocart_source = ',trim(gocart_source)
         stop
      endif
   end do

   ! Destroy interpolator object
   call destroy_interpolator(interp_object)
   
   ! The end
   print*,'Successful exit from gocart2wrf!'
   
end program gocart2wrf
