!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  geos2wps
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! 2.0 version of geos2wrf.  Reads HDF4 or netCDF files, extracts 2D slices of
! data as requested by user, and outputs slices in WPS intermediate format.
! No derivation of missing variables is performed--this task is deferred to
! other tools that will operate on the WPS files themselves.
!
! REVISION:
! 17 Apr 2012 - First version.
!  6 Jun 2012 - Added HDFEOS option.
!------------------------------------------------------------------------------

program geos2wps

   ! Import modules
   use GEOS2WPS_mod
   use NamelistGEOS2WPS_mod

   ! Change defaults
   implicit none

   ! Local variables
   type(NamelistGEOS2WPS) :: namelist

   ! Process namelist input
   namelist = createNamelistGEOS2WPS('namelist.geos2wps')

   ! Process the GEOS file
   call driver(namelist)

end program geos2wps
