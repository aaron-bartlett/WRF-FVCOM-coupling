!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  casa_flux_filenames_mod
!
! AUTHOR:
! 
! Jossy Jacob, NASA !
! DESCRIPTION:
! Defines look-up tables with CO2 flux related files, units to be used for flux intrpolation
!------------------------------------------------------------------------------

module casa_flux_filenames_mod

   ! Change defaults
   implicit none
   private

   integer,parameter,public :: NUM_FLUX_VARIABLES = 6 ! number of FLUX related variables
!   character(len=30),parameter,public :: &
!        NUWRFname1 = 'wrfout_d01_2010-01-02_00:00:00' ! First NUWRF file name to process 
!   character(len=30),parameter,public :: &
!        NUWRFname2 = 'wrfout_d01_2010-01-09_00:00:00' ! Last NUWRF file name to process
   character(len=17),parameter,public :: &
        NPPname = 'NPP0_2010_mon0.nc' 
   character(len=21),parameter,public :: &
        FIREname = 'FIRE_2010_daily_01.nc' 
   character(len=17),parameter,public :: &
        RESPname = 'RESP_2010_mon0.nc' 
   character(len=17),parameter,public :: &
        BFUELname = 'BFUE_2010_mon0.nc' 
   character(len=17),parameter,public :: &
        OceanCO2name = 'OCO2_2010_mon0.nc' 
   character(len=17),parameter,public :: &
        FossilFuelname = 'FFUE_2010_year.nc' 


      real,dimension(NUM_FLUX_VARIABLES),parameter,public :: &
        flux_econversion_factors = &
        (/1.0, 1.0 , 1.0, 1.0, 1.0, 1.0 /)
      real,dimension(4), parameter, public :: & 
        qclimit = (/20.0,60.0, -95.0, -60.0/)

end module casa_flux_filenames_mod
