!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  casa_emission_variables_mod
!
! AUTHOR:
! Eric Kemp (Part of gocart2wrf utility)
! Modified by Jossy Jacob, NASA SSSO Oct 2013
!
! DESCRIPTION:
! Defines look-up tables summarizing CASA CO2 variables, units, etc. to be
! interpolated to WRF grids.
!------------------------------------------------------------------------------

module casa_emission_variables_mod

   ! Change defaults
   implicit none
   private

   integer,parameter,public :: NUM_CASA_eVARIABLES = 1
   character(len=9),dimension(NUM_CASA_eVARIABLES),parameter,public :: &
        casa_data_Evariable_names = &
        (/'CO2'/) 
   character(len=9),dimension(NUM_CASA_eVARIABLES),parameter,public :: &
        casa_wrf_Evariable_names = &
        (/'E_CASACO2'/) 
   character(len=5),dimension(NUM_CASA_eVARIABLES),parameter,public :: &
        casa_data_Evariable_units = &
        (/'Kg/m^2/s '/)
   character(len=15),dimension(NUM_CASA_eVARIABLES),parameter,public :: &
        casa_wrf_Evariable_units = &
	  (/'mol km^-2 hr^-1'/)
   character(len=30),dimension(NUM_CASA_eVARIABLES),parameter,public :: &
        casa_wrf_Evariable_descriptions = &
	 (/ "EMISSIONS of CASACO2 fluxes   " /)
   real,dimension(NUM_CASA_eVARIABLES),parameter,public :: &
        casa_econversion_factors = &
        (/3.0E+11/)

end module casa_emission_variables_mod
