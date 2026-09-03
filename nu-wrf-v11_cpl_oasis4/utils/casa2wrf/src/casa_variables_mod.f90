!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  casa_variables_mod
!
! AUTHOR:
! Eric Kemp (Part of gocart2wrf utility)
! Modified by Jossy Jacob, NASA SSSO Oct 2013
!
! DESCRIPTION:
! Defines look-up tables summarizing CASA CO2 variables, units, etc. to be
! interpolated to WRF grids.
!------------------------------------------------------------------------------

module casa_variables_mod

   ! Change defaults
   implicit none
   private

   integer,parameter,public :: NUM_CASA_VARIABLES = 1
   character(len=9),dimension(NUM_CASA_VARIABLES),parameter,public :: &
        casa_data_variable_names = &
        (/'CO2c  '/)
   character(len=9),dimension(NUM_CASA_VARIABLES),parameter,public :: &
        casa_wrf_variable_names = &
        (/'CASACO2  '/) 
   character(len=5),dimension(NUM_CASA_VARIABLES),parameter,public :: &
        casa_data_variable_units = &
        (/'ppmv '/)
   character(len=15),dimension(NUM_CASA_VARIABLES),parameter,public :: &
        casa_wrf_variable_units = &
        (/'ppmv           '/)
   character(len=30),dimension(NUM_CASA_VARIABLES),parameter,public :: &
        casa_wrf_variable_descriptions = &
        (/"Total CASACO2 " /)
   real,dimension(NUM_CASA_VARIABLES),parameter,public :: &
        casa_conversion_factors = &
        (/1.0E+6/)

end module casa_variables_mod
