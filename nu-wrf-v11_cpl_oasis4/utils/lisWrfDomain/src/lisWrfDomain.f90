!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  lisWrfDomain
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Reads in namelist.wps file and netCDF output files from geogrid.exe, and 
! outputs table of grid settings for insertion into a lis.config file.
! Note that the actual insertion of values will be done either by hand
! or by an external script.
!
! REVISION HISTORY:
! 09 May 2010 - Initial version (processed metgrid.exe output)
! 16 Jul 2012 - 2.0 version for geogrid.exe files.
!------------------------------------------------------------------------------

program lisWrfDomain

   ! Import modules
   use DomainData_mod

   ! Force explicit declarations
   implicit none

   ! Data structure with grid data from namelist.wps
   type(DomainData) :: wpsData

   ! Fetch grid data from namelist.wps
   wpsData = createDomainData()

   ! Write data to text file
   call writeDomainData(wpsData)

   ! Clean up
   call destroyDomainData(wpsData)

   stop
end program lisWrfDomain
