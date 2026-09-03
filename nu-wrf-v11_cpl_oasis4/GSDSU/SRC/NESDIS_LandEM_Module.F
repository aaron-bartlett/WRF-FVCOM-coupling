!--------------------------------------------------------------------------------
!M+
! NAME:
!       NESDIS_LandEM_Module
!
! PURPOSE:
!       Module containing the microwave land emissivity model
!       Module containing the parameters related to microwave snow emissivity model
!       Module containing the microwave snow emissivity model from SSMIS
!
! CATEGORY:
!       Surface : MW Surface Land Emissivity
!       Surface : MW Surface Snow Emissivity Model Parameters
!       Surface : MW Surface Snow  Emissivity from SSMIS
!
! REFERENCES:
!
!  (1) Yan, B., F. Weng and K.Okamoto,2004: "A microwave snow emissivity model, 8th Specialist Meeting on
!
!  Microwave Radiometry and Remote Sension Applications, 24-27 February, 2004, Rome, Italy.
!
!  (2) F. Weng, B. Yan, and N. Sun, Correction of SSMIS Radiance Anomalies,NRL/JCSDA Mini-Workshop on Preparation
!
!  for SSMIS Radiance Assimilation, 26-27 October, 2005, Monterey, CA.
!
!  (3) B. Yan, F. Weng, and Tsan Mo, Calibration of DMSP F-16 Special Sensor Microwave Imager and Sounder, submitted
!
!  to the 14th Conference on Satellite Meteorology and Oceanography, 29 JAN ~ 2 FEB, 2006, ATLANTA, GA
!
!  (4) B. Yan and F. Weng, Assessments of Special Sensor Microwave Imager and Sounder (SSMIS) Data for NOAA Operational
!
!  Applications,submitted to the 9th Specialist Meeting on Microwave Radiometry and Remote Sensing Applications,
!
!  San Juan, Puerto Rico, from 28 February to 03 March 2006.
!
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       USE LandEM_Module
!
! MODULES:
!       Type_Kinds:          Module containing definitions for kinds of variable types.
!
! CONTAINS:
!   PUBLIC SUBPROGRAM:
!
!       NESDIS_LandEM     : Subroutine to calculate the microwave land emissivity
!
!       NESDIS_SSMIS_SNOW : Subroutine to calculate snow emissivity from SSMIS
!
!
!   PRIVATE SUBPROGRAM:
!
!       Canopy_Diel        : Subroutine to calculate the dielectric constants of vegetation canopy
!
!       Soil_Diel          : Subroutine to calculate the dielectric properties of soil
!
!       Snow_Diel          : Subroutine to calculate the dielectric properties of snow
!
!       Reflectance        : Subroutine to compute the surface reflectivety using fresnel equations
!
!       Transmittance      : Subroutine to compute the surface transmitance using fresnel equations
!
!       Roughness_Reflectance : Subroutine to compute the surface reflectivety using fresnel equations
!                                      for a rough surface having a standard deviation of height of sigma
!       Canopy_Optic       : Subroutine to compute optic parameters for canopy
!
!       Snow_Optic         : Subroutine to compute optic parameters for snow
!
!       Two_Stream_Solution: Subroutine to simulate microwave emissivity over land conditions using
!                                   two stream approximation RTEs
!
!       SSMIS_SNOW_TB   : Subroutine to calculate the snow microwave emissivity from SSMIS TB
!
!       SSMIS_SNOW_TBTS : Subroutine to calculate the snow microwave emissivity from SSMIS TB & TS
!
! INCLUDE FILES:
!       None.
!
! EXTERNALS:
!       None.
!
! COMMON BLOCKS:
!       None.
!
! FILES ACCESSED:
!       None.
!
! CREATION HISTORY:
!       Written by:     Banghua Yan, QSS Group Inc., Banghua.Yan@noaa.gov (06-01-2005) & (06-03-2005) & (10-Jan-2006)
!
!
!       and             Fuzhong Weng, NOAA/NESDIS/ORA, Fuzhong.Weng@noaa.gov
!
!     Fixed bugs and added quality control: Banghua Yan, Quanhua Liu, Hong Yan (09-10-2005)
!
!
!  Copyright (C) 2005 Fuzhong Weng and Banghua Yan
!
!  This program is free software; you can redistribute it and/or modify it under the terms of the GNU
!  General Public License as published by the Free Software Foundation; either version 2 of the License,
!  or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
!  License for more details.
!
!  You should have received a copy of the GNU General Public License along with this program; if not, write
!  to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!M-
!--------------------------------------------------------------------------------

MODULE NESDIS_LandEM_Module

  ! ----------
  ! Module use
  ! ----------

  USE Type_Kinds

  ! -----------------------
  ! Disable implicit typing
  ! -----------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE


  PUBLIC  :: NESDIS_LandEM
  PUBLIC  :: NESDIS_SSMIS_SNOW


  real(fp_kind),public, parameter:: zero = 0.0_fp_kind
  real(fp_kind),public, parameter:: one_tenth = 0.1_fp_kind
  real(fp_kind),public, parameter:: half = 0.5_fp_kind
  real(fp_kind),public, parameter:: one = 1.0_fp_kind
  real(fp_kind),public, parameter:: two = 2.0_fp_kind
  real(fp_kind),public, parameter:: three = 3.0_fp_kind
  real(fp_kind),public, parameter:: four = 4.0_fp_kind
  real(fp_kind),public, parameter:: pi = 3.14159_fp_kind
  real(fp_kind),public, parameter:: emissh_default = 0.25_fp_kind
  real(fp_kind),public, parameter:: emissv_default = 0.30_fp_kind

!
! Snow EM Parameters
!

  INTEGER, PUBLIC, PARAMETER :: INVALID_SNOW_TYPE   = -999


  ! ----------------
  ! Valid snow types
  ! ----------------

  INTEGER, PUBLIC, PARAMETER :: WET_SNOW            =  1
  INTEGER, PUBLIC, PARAMETER :: GRASS_AFTER_SNOW    =  2
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_A           =  3
  INTEGER, PUBLIC, PARAMETER :: POWDER_SNOW         =  4
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_B           =  5
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_C           =  6
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_D           =  7
  INTEGER, PUBLIC, PARAMETER :: THIN_CRUST_SNOW     =  8
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_E           =  9
  INTEGER, PUBLIC, PARAMETER :: BOTTOM_CRUST_SNOW_A = 10
  INTEGER, PUBLIC, PARAMETER :: SHALLOW_SNOW        = 11
  INTEGER, PUBLIC, PARAMETER :: DEEP_SNOW           = 12
  INTEGER, PUBLIC, PARAMETER :: CRUST_SNOW          = 13
  INTEGER, PUBLIC, PARAMETER :: MEDIUM_SNOW         = 14
  INTEGER, PUBLIC, PARAMETER :: BOTTOM_CRUST_SNOW_B = 15
  INTEGER, PUBLIC, PARAMETER :: THICK_CRUST_SNOW    = 16

  INTEGER(ip_kind), PUBLIC, PARAMETER                          :: N_FREQUENCY  = 10
  INTEGER(ip_kind), PUBLIC, PARAMETER                          :: N_FREQ_AMSRE = 7


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION (N_FREQUENCY)    ::                                &
  FREQUENCY_DEFAULT = (/ 4.9_fp_kind,  6.93_fp_kind, 10.65_fp_kind, 18.7_fp_kind,23.8_fp_kind,   &
                        31.4_fp_kind, 50.3_fp_kind,  52.5_fp_kind, 89.0_fp_kind,150._fp_kind/)



  ! Define sixteen MW weighted emissivity spectra for AMSU ALGORITHMS

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  WET_SNOW_EMISS = (/0.87_fp_kind,0.89_fp_kind,0.91_fp_kind,0.93_fp_kind,0.94_fp_kind,           &
                          0.94_fp_kind,0.94_fp_kind,0.93_fp_kind,0.92_fp_kind,0.90_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  GRASS_AFTER_SNOW_EMISS = (/0.91_fp_kind,0.91_fp_kind,0.92_fp_kind,0.91_fp_kind,                &
                                  0.90_fp_kind,0.90_fp_kind,0.91_fp_kind,0.91_fp_kind,           &
                                  0.91_fp_kind,0.86_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  RS_SNOW_A_EMISS = (/0.90_fp_kind,0.89_fp_kind,0.88_fp_kind,0.87_fp_kind, 0.86_fp_kind,         &
                           0.86_fp_kind,0.85_fp_kind,0.85_fp_kind,0.82_fp_kind,0.82_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  POWDER_SNOW_EMISS = (/0.91_fp_kind,0.91_fp_kind,0.93_fp_kind,0.93_fp_kind,0.93_fp_kind,        &
                             0.93_fp_kind,0.89_fp_kind,0.88_fp_kind,0.79_fp_kind,0.79_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  RS_SNOW_B_EMISS = (/0.90_fp_kind,0.89_fp_kind,0.88_fp_kind,0.85_fp_kind,0.84_fp_kind,          &
                           0.83_fp_kind,0.83_fp_kind,0.82_fp_kind,0.79_fp_kind,0.73_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  RS_SNOW_C_EMISS = (/0.90_fp_kind,0.89_fp_kind,0.86_fp_kind,0.82_fp_kind,0.80_fp_kind,          &
                           0.79_fp_kind,0.78_fp_kind,0.78_fp_kind,0.77_fp_kind,0.77_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  RS_SNOW_D_EMISS = (/0.88_fp_kind,0.86_fp_kind,0.85_fp_kind,0.80_fp_kind,0.78_fp_kind,          &
                           0.77_fp_kind,0.77_fp_kind,0.76_fp_kind,0.72_fp_kind,0.72_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  THIN_CRUST_SNOW_EMISS = (/0.93_fp_kind,0.94_fp_kind,0.96_fp_kind,0.96_fp_kind,0.95_fp_kind,    &
                                 0.93_fp_kind,0.87_fp_kind,0.86_fp_kind,0.74_fp_kind,0.65_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  RS_SNOW_E_EMISS = (/0.87_fp_kind,0.86_fp_kind,0.84_fp_kind,0.80_fp_kind,0.76_fp_kind,          &
                           0.76_fp_kind,0.75_fp_kind,0.75_fp_kind,0.70_fp_kind,0.69_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  BOTTOM_CRUST_SNOW_A_EMISS = (/0.87_fp_kind,0.86_fp_kind,0.83_fp_kind,0.77_fp_kind,             &
                                     0.73_fp_kind,0.68_fp_kind,0.66_fp_kind,0.66_fp_kind,        &
                                     0.68_fp_kind,0.67_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  SHALLOW_SNOW_EMISS = (/0.89_fp_kind,0.89_fp_kind,0.88_fp_kind,0.87_fp_kind, 0.86_fp_kind,      &
                              0.82_fp_kind,0.77_fp_kind,0.76_fp_kind,0.69_fp_kind,0.64_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  DEEP_SNOW_EMISS = (/0.88_fp_kind,0.87_fp_kind,0.86_fp_kind,0.83_fp_kind,0.81_fp_kind,          &
                           0.77_fp_kind,0.74_fp_kind,0.73_fp_kind,0.69_fp_kind,0.64_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  CRUST_SNOW_EMISS = (/0.86_fp_kind,0.86_fp_kind,0.86_fp_kind,0.85_fp_kind,0.82_fp_kind,         &
                            0.78_fp_kind,0.69_fp_kind,0.68_fp_kind,0.51_fp_kind,0.47_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  MEDIUM_SNOW_EMISS = (/0.89_fp_kind,0.88_fp_kind,0.87_fp_kind,0.83_fp_kind,0.80_fp_kind,        &
                             0.75_fp_kind,0.70_fp_kind,0.70_fp_kind,0.64_fp_kind,0.60_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  BOTTOM_CRUST_SNOW_B_EMISS = (/0.91_fp_kind,0.92_fp_kind,0.93_fp_kind,0.88_fp_kind,             &
                                     0.84_fp_kind,0.76_fp_kind,0.66_fp_kind,0.64_fp_kind,        &
                                     0.48_fp_kind,0.44_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQUENCY)  ::                                   &
  THICK_CRUST_SNOW_EMISS = (/0.94_fp_kind,0.95_fp_kind,0.97_fp_kind,0.91_fp_kind,                &
                                  0.86_fp_kind,0.74_fp_kind,0.63_fp_kind,0.63_fp_kind,           &
                                  0.50_fp_kind,0.45_fp_kind/)


!
  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION (N_FREQ_AMSRE)    ::                                &
  FREQUENCY_AMSRE = (/ 6.925_fp_kind, 10.65_fp_kind, 18.7_fp_kind,23.8_fp_kind,                   &
                        36.5_fp_kind, 89.0_fp_kind,150._fp_kind/)


 ! Define sixteen MW weighted emissivity spectra for AMSRE ALGORITHMS

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  WET_SNOW_EM_AMSRE = (/0.91_fp_kind, 0.93_fp_kind, 0.94_fp_kind, 0.95_fp_kind, 0.95_fp_kind,     &
                        0.93_fp_kind, 0.93_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  GRASS_AFTER_SNOW_EM_AMSRE = (/0.91_fp_kind, 0.92_fp_kind, 0.91_fp_kind, 0.90_fp_kind,           &
                                0.91_fp_kind, 0.91_fp_kind, 0.91_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_A_EM_AMSRE = (/0.90_fp_kind, 0.89_fp_kind, 0.88_fp_kind, 0.87_fp_kind, 0.86_fp_kind,    &
                         0.82_fp_kind, 0.82_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  POWDER_SNOW_EM_AMSRE = (/0.92_fp_kind, 0.93_fp_kind, 0.94_fp_kind, 0.94_fp_kind, 0.92_fp_kind,  &
                           0.80_fp_kind, 0.80_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_B_EM_AMSRE = (/0.87_fp_kind, 0.86_fp_kind, 0.83_fp_kind, 0.80_fp_kind, 0.79_fp_kind,    &
                         0.77_fp_kind, 0.77_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_C_EM_AMSRE = (/0.89_fp_kind, 0.88_fp_kind, 0.85_fp_kind, 0.84_fp_kind, 0.83_fp_kind,    &
                         0.79_fp_kind, 0.79_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_D_EM_AMSRE = (/0.84_fp_kind, 0.83_fp_kind, 0.82_fp_kind, 0.80_fp_kind,            &
                               0.78_fp_kind, 0.72_fp_kind, 0.72_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  THIN_CRUST_SNOW_EM_AMSRE = (/0.95_fp_kind, 0.96_fp_kind, 0.96_fp_kind, 0.95_fp_kind,            &
                               0.91_fp_kind, 0.75_fp_kind, 0.75_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_E_EM_AMSRE = (/0.80_fp_kind, 0.80_fp_kind, 0.80_fp_kind, 0.79_fp_kind,                  &
                         0.75_fp_kind, 0.70_fp_kind, 0.70_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  BOTTOM_CRUST_SNOW_A_EM_AMSRE = (/0.91_fp_kind, 0.90_fp_kind, 0.89_fp_kind, 0.87_fp_kind,        &
                                   0.82_fp_kind, 0.69_fp_kind, 0.69_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  SHALLOW_SNOW_EM_AMSRE = (/0.90_fp_kind, 0.89_fp_kind, 0.85_fp_kind,0.82_fp_kind, 0.76_fp_kind, &
                         0.65_fp_kind, 0.65_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  DEEP_SNOW_EM_AMSRE = (/0.89_fp_kind, 0.88_fp_kind, 0.86_fp_kind, 0.83_fp_kind, 0.78_fp_kind,    &
                          0.70_fp_kind, 0.70_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  CRUST_SNOW_EM_AMSRE = (/0.88_fp_kind, 0.86_fp_kind, 0.80_fp_kind, 0.75_fp_kind, 0.69_fp_kind,   &
                           0.67_fp_kind, 0.67_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  MEDIUM_SNOW_EM_AMSRE = (/0.96_fp_kind, 0.97_fp_kind, 0.92_fp_kind, 0.87_fp_kind, 0.72_fp_kind,  &
                           0.50_fp_kind, 0.50_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  BOTTOM_CRUST_SNOW_B_EM_AMSRE = (/0.93_fp_kind, 0.94_fp_kind, 0.89_fp_kind, 0.85_fp_kind,        &
                                   0.74_fp_kind, 0.48_fp_kind, 0.48_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  THICK_CRUST_SNOW_EM_AMSRE = (/0.88_fp_kind, 0.88_fp_kind, 0.87_fp_kind, 0.85_fp_kind,           &
                                0.77_fp_kind, 0.52_fp_kind, 0.52_fp_kind/)


 ! Define sixteen MW H-POL emissivity spectra for AMSRE ALGORITHMS


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  WET_SNOW_EH_AMSRE = (/0.93_fp_kind, 0.92_fp_kind, 0.93_fp_kind, 0.94_fp_kind, 0.93_fp_kind,     &
                        0.93_fp_kind, 0.90_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  GRASS_AFTER_SNOW_EH_AMSRE = (/0.91_fp_kind, 0.90_fp_kind, 0.90_fp_kind, 0.90_fp_kind,           &
                                0.91_fp_kind, 0.90_fp_kind, 0.85_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_A_EH_AMSRE = (/0.85_fp_kind, 0.85_fp_kind, 0.84_fp_kind, 0.84_fp_kind, 0.82_fp_kind,    &
                         0.80_fp_kind, 0.80_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  POWDER_SNOW_EH_AMSRE = (/0.90_fp_kind, 0.90_fp_kind, 0.92_fp_kind, 0.92_fp_kind, 0.90_fp_kind,  &
                           0.80_fp_kind, 0.79_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_B_EH_AMSRE = (/0.82_fp_kind, 0.81_fp_kind, 0.77_fp_kind, 0.76_fp_kind, 0.74_fp_kind,    &
                         0.74_fp_kind, 0.74_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_C_EH_AMSRE = (/0.84_fp_kind, 0.83_fp_kind, 0.80_fp_kind, 0.78_fp_kind, 0.77_fp_kind,   &
                         0.75_fp_kind, 0.69_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_D_EH_AMSRE = (/0.77_fp_kind, 0.77_fp_kind, 0.76_fp_kind, 0.75_fp_kind, 0.73_fp_kind,    &
                         0.71_fp_kind, 0.71_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  THIN_CRUST_SNOW_EH_AMSRE = (/0.95_fp_kind, 0.94_fp_kind, 0.95_fp_kind, 0.94_fp_kind,            &
                               0.89_fp_kind, 0.75_fp_kind, 0.65_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_E_EH_AMSRE = (/0.73_fp_kind, 0.73_fp_kind, 0.74_fp_kind, 0.72_fp_kind, 0.71_fp_kind,    &
                         0.68_fp_kind, 0.67_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  BOTTOM_CRUST_SNOW_A_EH_AMSRE = (/0.88_fp_kind, 0.87_fp_kind, 0.86_fp_kind, 0.85_fp_kind,        &
                                   0.80_fp_kind, 0.68_fp_kind, 0.63_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  SHALLOW_SNOW_EH_AMSRE = (/0.86_fp_kind, 0.84_fp_kind, 0.80_fp_kind, 0.78_fp_kind,               &
                            0.72_fp_kind, 0.62_fp_kind, 0.57_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  DEEP_SNOW_EH_AMSRE = (/0.87_fp_kind, 0.85_fp_kind, 0.83_fp_kind, 0.80_fp_kind, 0.77_fp_kind,    &
                         0.68_fp_kind, 0.62_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  CRUST_SNOW_EH_AMSRE = (/0.82_fp_kind, 0.78_fp_kind, 0.74_fp_kind, 0.71_fp_kind, 0.67_fp_kind,   &
                          0.64_fp_kind, 0.64_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  MEDIUM_SNOW_EH_AMSRE = (/0.90_fp_kind, 0.90_fp_kind, 0.89_fp_kind, 0.88_fp_kind, 0.83_fp_kind,  &
                           0.53_fp_kind, 0.48_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  BOTTOM_CRUST_SNOW_B_EH_AMSRE = (/0.87_fp_kind, 0.85_fp_kind, 0.84_fp_kind, 0.82_fp_kind,        &
                                   0.74_fp_kind, 0.53_fp_kind, 0.49_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  THICK_CRUST_SNOW_EH_AMSRE = (/0.85_fp_kind, 0.84_fp_kind, 0.83_fp_kind, 0.81_fp_kind,           &
                                0.79_fp_kind, 0.51_fp_kind, 0.46_fp_kind/)



 ! Define sixteen MW V-POL emissivity spectra for AMSRE ALGORITHMS

 REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  WET_SNOW_EV_AMSRE = (/0.96_fp_kind, 0.94_fp_kind, 0.96_fp_kind, 0.95_fp_kind, 0.94_fp_kind,     &
                        0.94_fp_kind, 0.91_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  GRASS_AFTER_SNOW_EV_AMSRE = (/0.96_fp_kind, 0.94_fp_kind, 0.95_fp_kind, 0.96_fp_kind,           &
                                0.96_fp_kind, 0.92_fp_kind, 0.87_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_A_EV_AMSRE = (/0.99_fp_kind, 0.97_fp_kind, 0.96_fp_kind, 0.96_fp_kind, 0.93_fp_kind,    &
                         0.87_fp_kind, 0.87_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  POWDER_SNOW_EV_AMSRE = (/0.98_fp_kind, 0.97_fp_kind, 0.99_fp_kind, 0.98_fp_kind, 0.96_fp_kind,  &
                           0.84_fp_kind, 0.83_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_B_EV_AMSRE = (/0.97_fp_kind, 0.95_fp_kind, 0.93_fp_kind, 0.92_fp_kind, 0.89_fp_kind,    &
                         0.84_fp_kind,0.84_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_C_EV_AMSRE = (/1.00_fp_kind, 0.97_fp_kind, 0.96_fp_kind, 0.94_fp_kind, 0.91_fp_kind,    &
                         0.84_fp_kind, 0.78_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_D_EV_AMSRE = (/0.99_fp_kind, 0.96_fp_kind, 0.93_fp_kind, 0.90_fp_kind, 0.86_fp_kind,    &
                         0.80_fp_kind, 0.80_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  THIN_CRUST_SNOW_EV_AMSRE = (/0.98_fp_kind, 0.97_fp_kind, 0.98_fp_kind, 0.97_fp_kind,            &
                               0.92_fp_kind, 0.77_fp_kind, 0.67_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  RS_SNOW_E_EV_AMSRE = (/0.98_fp_kind, 0.95_fp_kind, 0.90_fp_kind, 0.86_fp_kind, 0.82_fp_kind,    &
                         0.74_fp_kind, 0.73_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  BOTTOM_CRUST_SNOW_A_EV_AMSRE = (/0.96_fp_kind, 0.95_fp_kind, 0.95_fp_kind, 0.93_fp_kind,        &
                                   0.87_fp_kind, 0.71_fp_kind, 0.66_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  SHALLOW_SNOW_EV_AMSRE = (/0.97_fp_kind, 0.95_fp_kind, 0.94_fp_kind, 0.90_fp_kind, 0.84_fp_kind, &
                            0.68_fp_kind, 0.63_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  DEEP_SNOW_EV_AMSRE = (/0.96_fp_kind, 0.94_fp_kind, 0.92_fp_kind, 0.90_fp_kind, 0.85_fp_kind,    &
                         0.77_fp_kind, 0.71_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  CRUST_SNOW_EV_AMSRE = (/0.98_fp_kind, 0.96_fp_kind, 0.93_fp_kind, 0.90_fp_kind, 0.81_fp_kind,   &
                          0.71_fp_kind, 0.71_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  MEDIUM_SNOW_EV_AMSRE = (/0.99_fp_kind, 0.97_fp_kind, 0.98_fp_kind, 0.96_fp_kind, 0.92_fp_kind,  &
                           0.57_fp_kind, 0.52_fp_kind/)


  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  BOTTOM_CRUST_SNOW_B_EV_AMSRE = (/1.00_fp_kind, 0.97_fp_kind, 0.97_fp_kind, 0.95_fp_kind,        &
                                   0.86_fp_kind, 0.58_fp_kind, 0.54_fp_kind/)

  REAL(fp_kind), PUBLIC, PARAMETER, DIMENSION(N_FREQ_AMSRE)  ::                                   &
  THICK_CRUST_SNOW_EV_AMSRE = (/0.98_fp_kind, 0.96_fp_kind, 0.96_fp_kind, 0.94_fp_kind,           &
                                0.89_fp_kind, 0.56_fp_kind, 0.51_fp_kind/)



CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!-------------------------------------------------------------------------------------------------------------
!
! NAME:
!       NESDIS_LandEM
!
! PURPOSE:
!       Subroutine to simulate microwave emissivity over land conditions.
!
! REFERENCES:
!       Weng, F., B. Yan, and N. Grody, 2001: "A microwave land emissivity model", J. Geophys. Res., 106,
!                                             20, 115-20, 123
!
! CATEGORY:
!       CRTM : Surface : MW LandEM
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL NESDIS_LandEM
!
! INPUT ARGUMENTS:
!
!         Frequency                Frequency User defines
!                                  This is the "I" dimension
!                                  UNITS:      GHz
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!
!         Angle:                   The angle values in degree.
!                                  ** NOTE: THIS IS A MANDATORY MEMBER **
!                                  **       OF THIS STRUCTURE          **
!                                  UNITS:      Degrees
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Rank-1, (I)
!         t_skin = Land_Temperature:        The land surface temperature.
!                                  UNITS:      Kelvin, K
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!         mv = Soil_Moisture_Content:   The volumetric water content of the soil(0~1).
!                                  UNITS:      demensionless
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!         mg = Canopy_Water_Content:The gravimetric water content of the canopy (0~1)
!                                  UNITS:      demensionless
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!         veg_frac = Vegetation_Fraction:     The vegetation fraction of the surface(0:1).
!                                  UNITS:      N/A
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!         t_soil = Soil_Temperature:        The soil temperature.
!                                  UNITS:      Kelvin, K
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!         Snow_Depth:              The snow depth.
!                                  UNITS:      mm
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!
! OUTPUT ARGUMENTS:
!
!         Emissivity_H:            The surface emissivity at a horizontal polarization.
!                                  ** NOTE: THIS IS A MANDATORY MEMBER **
!                                  **       OF THIS STRUCTURE          **
!                                  UNITS:      N/A
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!         Emissivity_V:            The surface emissivity at a vertical polarization.
!                                  ** NOTE: THIS IS A MANDATORY MEMBER **
!                                  **       OF THIS STRUCTURE          **
!                                  UNITS:      N/A
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!
! INTERNAL ARGUMENTS:
!       theta       -  local zenith angle in radian
!       rhob        -  bulk volume density of the soil (1.18-1.12)
!       rhos        -  density of the solids (2.65 g.cm^3 for solid soil material)
!       sand        -  sand fraction (sand + clay = 1.0)
!       clay        -  clay fraction
!       lai         -  leaf area index (eg. lai = 4.0 for corn leaves)
!       sigma       -  surface roughness formed between medium 1 and 2,
!                      expressed as the standard deviation of roughtness height (mm)
!       leaf_thick  --  leaf thickness (mm)
!       rad         -  radius of dense medium scatterers (mm)
!       va          -  fraction volume of dense medium scatterers(0.0 - 1.0)
!       ep          -  dielectric constant of ice or sand particles, complex value
!                               (e.g, 3.0+i0.0)
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       Canopy_Diel        : Subroutine to calculate the dielectric constants of vegetation canopy
!
!       Soil_Diel          : Subroutine to calculate the dielectric properties of soil
!
!       Snow_Diel          : Subroutine to calculate the dielectric properties of snow
!
!       Reflectance        : Subroutine to compute the surface reflectivety using fresnel equations
!
!       Transmittance      : Subroutine to compute the surface transmitance using fresnel equations
!
!       Roughness_Reflectance : Subroutine to compute the surface reflectivety using fresnel equations
!                                      for a rough surface having a standard deviation of height of sigma
!       Canopy_Optic       : Subroutine to compute optic parameters for canopy
!
!       Snow_Optic         : Subroutine to compute optic parameters for snow
!
!       Two_Stream_Solution: Subroutine to simulate microwave emissivity over land conditions using
!                                   two stream approximation RTEs
!
!
! SIDE EFFECTS:
!       None.
!
! RESTRICTIONS:
!       None.
!
! COMMENTS:
!       Note the INTENT on the output SensorData argument is IN OUT rather than
!       just OUT. This is necessary because the argument may be defined upon
!       input. To prevent memory leaks, the IN OUT INTENT is a must.
!
! CREATION HISTORY:
!       Written by:     Banghua Yan, QSS Group Inc., Banghua.Yan@noaa.gov (16-May-2005)
!
!
!       and             Fuzhong Weng, NOAA/NESDIS/ORA, Fuzhong.Weng@noaa.gov
!
!      Fixed bugs and added quality controls:
!
!                        Banghua Yan, Quanhua Liu and Hong Han  (10-September-2005)
!
!      Fixed bugs and added the parameteried value of tsoil under snows suggested by Paul (11/03/06)
!
!  Copyright (C) 2005 Fuzhong Weng and Banghua Yan
!
!  This program is free software; you can redistribute it and/or modify it under the terms of the GNU
!  General Public License as published by the Free Software Foundation; either version 2 of the License,
!  or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
!  License for more details.
!
!  You should have received a copy of the GNU General Public License along with this program; if not, write
!  to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!
!------------------------------------------------------------------------------------------------------------

  SUBROUTINE NESDIS_LandEM(Angle,                                     &   ! Input
                           Frequency,                                 &   ! Input
                           mv,                                        &   ! Input
                           veg_frac,                                  &   ! Input
                           t_soil,                                    &   ! Input
                           t_skin,                                    &   ! Input
                           Snow_Depth,                                &   ! Input
                           Emissivity_H,                              &   ! Output
                           Emissivity_V)                                  ! Output

  use type_kinds, only: fp_kind, ip_kind
  implicit none

! Declare passed variables
  real(fp_kind),intent(in) :: Angle,Frequency,t_skin
  real(fp_kind),intent(out):: Emissivity_V,Emissivity_H
  real(fp_kind):: Emiss_temp_V,Emiss_temp_H
  real(fp_kind)            :: t_soil

! Declare local parameters

  integer(ip_kind),parameter :: PHYSICAL_MODEL = 1, EMPRIRICAL_METHOD = 2 
  real(fp_kind),parameter :: snow_depth_c = 10.0, tsoilc_undersnow = 280.0
  real(fp_kind),parameter:: rhob = 1.18_fp_kind  !toshii hardwire with LSM
  real(fp_kind),parameter:: rhos = 2.65_fp_kind  !toshii hardwire with LSM
  real(fp_kind),parameter:: sand = 0.8_fp_kind   !toshii hardwire with LSM  (very small impact)
  real(fp_kind),parameter:: clay = 0.2_fp_kind   !toshii hardwire with LSM

! Declare local variables

  integer(ip_kind) SNOWEM_APPROACH

  real(fp_kind) mv,veg_frac,Snow_Depth,theta,theta_i,theta_t,mu,r21_h,r21_v,r23_h,r23_v,  &
                t21_v,t21_h,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v,mge, &
                lai,leaf_thick,rad,sigma,va,ep_real,ep_imag
!toshii(added)
  real(fp_kind) :: esoil_V, esoil_H, eveg_V, eveg_H, wgt

  complex esoil, eveg, esnow, eair

  eair = cmplx(one,-zero)

  theta = Angle*pi/180.0_fp_kind

! Quality Control

  if ( (t_soil .le. 100.0_fp_kind .or. t_soil .ge. 350.0_fp_kind) .AND.    &

       (t_skin .ge. 100.0_fp_kind .and. t_skin .le. 350.0_fp_kind) ) t_soil = t_skin


  if (mv > one)  mv = one         !domensional (zero ~ one)

  if (mv < zero) mv = zero        !domensional (zero ~ one)

!
! toshii--> hypotheis
! since microwave penetrate into depth until getting saturated.
!
!  mv = 0.3 

  if(Snow_Depth <= 10.e0) Snow_Depth = 0.e0  !zero out (toshii temporal)

  if (Snow_Depth .gt.one_tenth) then

     SNOWEM_APPROACH = PHYSICAL_MODEL

     if (Snow_Depth .gt. snow_depth_c)  SNOWEM_APPROACH = EMPRIRICAL_METHOD

     GET_SNOWEM: SELECT CASE (SNOWEM_APPROACH)

     CASE (PHYSICAL_MODEL)

      ep_real = 3.2_fp_kind      !toshii this should be temperature dependent??

      ep_imag = -0.0005_fp_kind  !toshii this is too??

      sigma = one

      if (Snow_Depth > 1000.0_fp_kind) Snow_Depth = 1000.0_fp_kind   ! For deep snow, the performance of the model is poor

      va = 0.4_fp_kind + 0.0004_fp_kind*Snow_Depth

      if( va > one ) va = one      ! the fraction volume of dense medium scatterers must be less than one.

      if( va < zero ) va = zero    ! the fraction volume of dense medium scatterers must be greater than/equal to zero.

      rad = half + 0.005_fp_kind*Snow_Depth !toshii diagnostic snow grain size??

      if( rad > one ) rad = one  ! Limit for snow grain size

      if (t_soil .ge. tsoilc_undersnow) t_soil = tsoilc_undersnow

      call Snow_Diel(Frequency, ep_real, ep_imag, rad, va, esnow)

      call Soil_Diel(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)

      theta_i = asin(real(sin(theta)*csqrt(eair)/csqrt(esnow)))

      call Reflectance(esnow, eair, theta_i,  theta, r21_v, r21_h)

      call Transmittance(esnow, eair, theta_i, theta, t21_v, t21_h)

      mu  = cos(theta_i)

      theta_t = asin(real(sin(theta_i)*csqrt(esnow)/csqrt(esoil)))

      call Reflectance(esnow, esoil, theta_i, theta_t, r23_v, r23_h)

      call Roughness_Reflectance(Frequency, sigma, r23_v, r23_h)

      call Snow_Optic(Frequency,rad,Snow_Depth,va,ep_real, ep_imag,gv,gh,&

                      ssalb_v,ssalb_h,tau_v,tau_h)


      call Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &

                               r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,Emissivity_V,Emissivity_H)


    CASE (EMPRIRICAL_METHOD)


      CALL SnowEM_Default(Frequency,t_skin, Snow_Depth,Emissivity_V,Emissivity_H)


    END SELECT GET_SNOWEM


  else  !Non-snow land surface (Start from here)

     sigma = half !tunable

     if (veg_frac > one) veg_frac = one

     if (veg_frac < zero) veg_frac = zero

!     lai = 3. !(three*veg_frac + half)   !toshii change this 
     lai  = -1.*log( 1. - veg_frac )  !beer's law 

     mge = half*veg_frac           !toshii change this

     leaf_thick = 0.07_fp_kind     !tunable

     mu  = cos(theta)

     r21_h    = zero

     r21_v    = zero

     t21_h    = one

     t21_v    = one

     call Soil_Diel(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)

     theta_t = asin(real(sin(theta)*csqrt(eair)/csqrt(esoil)))

     call Reflectance(eair, esoil, theta, theta_t, r23_v, r23_h)

     call Roughness_Reflectance(Frequency, sigma, r23_v, r23_h)

     !toshii artificially double horizontal reflectance
     !r23_h = 2.*r23_h
     call Canopy_Diel(Frequency, mge, eveg)

     call Canopy_Optic(lai,Frequency,theta,eveg,leaf_thick,gv,gh,ssalb_v,ssalb_h,tau_v,tau_h)

     call Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
                              r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,Emissivity_V,Emissivity_H)

! toshi--> modification  (new hybrid scheme)
!     call Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
!                              r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,eveg_V,eveg_H)
!      esoil_V = 1.-r23_v
!      esoil_H = 1.-r23_h
!      wgt = max(0., min(1., veg_frac**0.5 ))
!      Emissivity_V = wgt * eveg_V + (1.-wgt )*esoil_V
!      Emissivity_H = wgt * eveg_H + (1.-wgt )*esoil_H



  endif

  return


END SUBROUTINE NESDIS_LandEM



!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################



subroutine SnowEM_Default(frequency,ts, depth,Emissivity_V,Emissivity_H)

!----------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .
!   prgmmr:  Banghua Yan and Fuzhong Weng               org: nesdis              date: 2005-12-01
!
! abstract: preliminary estimate of snow emissivity using  surface temperature and snow depth
!
! input argument list:
!
!      ts         -  surface temperature
!      frequency   -  frequency (ghz)
!
! output argument list:
!
!      Emissivity_V         -  snow emissivty at V-POL
!      Emissivity_H         -  snow emissivty at H-POL
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!------------------------------------------------------------------------------------------------------------

  use type_kinds, only: fp_kind, ip_kind


  implicit none

  integer(ip_kind), parameter :: new = 7

  integer(ip_kind), parameter :: NFRESH_SHALLOW_SNOW = 1, NPOWDER_SNOW = 2, NWET_SNOW = 3, NDEEP_SNOW = 4

  real(fp_kind), parameter :: twet = 270.0, tcrust = 235.0, depth_s = 50.0, depth_c = 100.0

  integer(ip_kind) :: ich,basic_snow_type

  real(fp_kind)    :: frequency,ts, depth,Emissivity_V,Emissivity_H

  real(fp_kind), dimension(new) :: ev, eh, freq


  freq(1:new) = FREQUENCY_AMSRE(1:new)

  basic_snow_type = NFRESH_SHALLOW_SNOW

  if (ts .ge. twet .and. depth .le. depth_s) then

     basic_snow_type = NWET_SNOW

  else

     if (depth .le. depth_s) then

         basic_snow_type = NFRESH_SHALLOW_SNOW

     else

        basic_snow_type = NPOWDER_SNOW

     endif

  endif

  if (ts .le. tcrust .and. depth .ge. depth_c) basic_snow_type = NDEEP_SNOW

! INITIALIZATION

  ev(1:new) = GRASS_AFTER_SNOW_EV_AMSRE(1:new)

  eh(1:new) = GRASS_AFTER_SNOW_EH_AMSRE(1:new)

  GET_SNOWTYPE: SELECT CASE (basic_snow_type)

     CASE (NFRESH_SHALLOW_SNOW)

          ev(1:new) = GRASS_AFTER_SNOW_EV_AMSRE(1:new)

          eh(1:new) = GRASS_AFTER_SNOW_EH_AMSRE(1:new)

     CASE (NPOWDER_SNOW)

          ev(1:new) = POWDER_SNOW_EV_AMSRE(1:new)

          eh(1:new) = POWDER_SNOW_EH_AMSRE(1:new)

     CASE (NWET_SNOW)

         ev(1:new) = WET_SNOW_EV_AMSRE(1:new)

         eh(1:new) = WET_SNOW_EH_AMSRE(1:new)

     CASE (NDEEP_SNOW)

         ev(1:new) = DEEP_SNOW_EV_AMSRE(1:new)

         eh(1:new) = DEEP_SNOW_EH_AMSRE(1:new)

   END SELECT GET_SNOWTYPE




! Interpolate emissivity at a certain frequency

  do ich=1,new

     if (frequency .le. freq(1)) then

         Emissivity_H = eh(1)

         Emissivity_V = ev(1)

         exit

      endif

     if (frequency .ge. freq(new)) then

         Emissivity_H = eh(new)

         Emissivity_V = ev(new)

         exit

      endif


      if (frequency .le. freq(ich)) then

           Emissivity_H = eh(ich-1) + &

                          (frequency-freq(ich-1))*(eh(ich) - eh(ich-1))/(freq(ich)-freq(ich-1))

           Emissivity_V = ev(ich-1) + &

                          (frequency-freq(ich-1))*(ev(ich) - ev(ich-1))/(freq(ich)-freq(ich-1))

          exit

      endif

  enddo

end subroutine SnowEM_Default


subroutine Canopy_Optic(lai,frequency,theta,esv,d,gv,gh,&
     ssalb_v,ssalb_h,tau_v, tau_h)

!----------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    canopy_optic compute optic parameters for canopy
!
!   prgmmr:  Fuzhong Weng and Banghua Yan                org: nesdis              date: 2000-11-28
!
! abstract: compute optic parameters for canopy
!
! program history log:
!
! input argument list:
!
!      lai         -  leaf area index
!      frequency   -  frequency (ghz)
!      theta       -  incident angle
!      esv         -  leaf dielectric constant
!      d           -  leaf thickness (mm)
!
! output argument list:
!
!      gv           -  asymmetry factor for v pol
!      gh           -  asymmetry factor for h pol
!      ssalb_v      -  single scattering albedo at v. polarization
!      ssalb_h      -  single scattering albedo at h. polarization
!      tau_v        -  optical depth at v. polarization
!      tau_h        -  optical depth at h. polarization
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!------------------------------------------------------------------------------------------------------------

  use type_kinds, only: fp_kind

  implicit none

  real(fp_kind) threshold

  real(fp_kind) frequency,theta,d,lai,ssalb_v,ssalb_h,tau_v,tau_h,gv, gh, mu

  complex  ix,k0,kz0,kz1,rhc,rvc,esv,expval1,factt,factrvc,factrhc

  real(fp_kind) rh,rvert,th,tv


  threshold=0.999_fp_kind

  mu = cos(theta)

  ix  = cmplx(zero,one)

  k0  = cmplx(two*pi*frequency/300.0_fp_kind, zero)   ! 1/mm

  kz0 = k0*mu

  kz1 = k0*sqrt(esv - sin(theta)**2)

  rhc = (kz0 - kz1)/(kz0 + kz1)

  rvc = (esv*kz0 - kz1)/(esv*kz0 + kz1)

  expval1=exp(-two*ix*kz1*d)

  factrvc=one-rvc**2*expval1

  factrhc=one-rhc**2*expval1

  factt=four*kz0*kz1*exp(ix*(kz0-kz1)*d)

  rvert = abs(rvc*(one - expval1)/factrvc)**2

  rh = abs(rhc*(one - expval1)/factrhc)**2

  th = abs(factt/((kz1+kz0)**2*factrhc))**2

  tv = abs(esv*factt/((kz1+esv*kz0)**2*factrvc))**2

  gv = half

  gh = half

  tau_v = half*lai*(two-tv-th)

  tau_h = tau_v

  ssalb_v = min((rvert+rh)/ (two-tv-th),threshold)

  ssalb_h = ssalb_v

  return

end subroutine Canopy_Optic


subroutine Snow_Optic(frequency,a,h,f,ep_real,ep_imag,gv,gh, ssalb_v,ssalb_h,tau_v,tau_h)

!-------------------------------------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    landem      comput optic parameters for snow
!
!   prgmmr: Fuzhong Weng and Banghua Yan                 org: nesdis              date: 2000-11-28
!
! abstract: compute optic parameters for snow
!
! program history log:
!
! input argument list:
!
!      theta        -  local zenith angle (degree)
!      frequency    -  frequency (ghz)
!      ep_real      -  real part of dielectric constant of particles
!      ep_imag      -  imaginary part of dielectric constant of particles
!      a            -  particle radiu (mm)
!      h            -  snow depth(mm)
!      f            -  fraction volume of snow (0.0 - 1.0)
!
! output argument list:
!
!       ssalb       -  single scattering albedo
!       tau         -  optical depth
!       g           -  asymmetry factor
!
!   important internal variables:
!
!       ks          -  scattering coeffcient (/mm)
!       ka          -  absorption coeffient (/mm)
!       kp          -  eigenvalue of two-stream approximation
!       y           -  = yr+iyi
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!----------------------------------------------------------------------------------

  use type_kinds, only: fp_kind

  implicit none

  real(fp_kind) yr,yi,ep_real,ep_imag

  real(fp_kind) frequency,a,h,f,ssalb_v,ssalb_h,tau_v,tau_h,gv,gh,k

  real(fp_kind) ks1,ks2,ks3,ks,kr1,kr2,kr3,kr,ki1,ki2,ki3,ki

  real(fp_kind) fact1,fact2,fact3,fact4,fact5


  k = two*pi/(300._fp_kind/frequency)

  yr  = (ep_real - one)/(ep_real + two)

  yi =  - ep_imag/(ep_real + two)

  fact1 = (one+two*f)**2

  fact2 = one-f*yr

  fact3 = (one-f)**4

  fact4 = f*(k*a)**3

  fact5 = one+two*f*yr

  ks1 = k*sqrt(fact2/fact5)

  ks2 = fact4*fact3/fact1

  ks3 = (yr/fact2)**2

  ks = ks1*ks2*ks3

  kr1 = fact5/fact2

  kr2 = two*ks2

  kr3 = two*yi*yr/(fact2**3)

  kr = k*sqrt(kr1+kr2*kr3)

  ki1 = three*f*yi/fact2**2

  ki2 = kr2

  ki3 = ks3

  ki  = k**2/(two*kr)*(ki1+ki2*ki3)

  gv = half

  gh = half

  ssalb_v = ks / ki

  if(ssalb_v .gt. 0.999_fp_kind) ssalb_v = 0.999_fp_kind

  ssalb_h = ssalb_v

  tau_v = two*ki*h

  tau_h = tau_v

  return

end subroutine Snow_Optic


subroutine Soil_Diel(freq,t_soil,vmc,rhob,rhos,sand,clay,esm)

!----------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    Soil_Diel   calculate the dielectric properties of soil
!
!   prgmmr: Fuzhong Weng and Banghua Yan                 org: nesdis              date: 2000-11-28
!
! abstract: compute the dilectric constant of the bare soil
!
! program history log:
!
! input argument list:
!
!      theta        -  local zenith angle (degree)
!      frequency    -  frequency (ghz)
!      t_soil       -  soil temperature
!      vmc          -  volumetric moisture content (demensionless)
!      rhob         -  bulk volume density of the soil (1.18-1.12)
!      rhos         -  density of the solids (2.65 g.cm^3 for
!                       solid soil material)
!      sand         -  sand fraction (sand + clay = 1.0)
!      clay         -  clay fraction
!
! output argument list:
!
!      esm          -  dielectric constant for bare soil
!
! important internal variables:
!
!      esof         -  the permittivity of free space
!      eswo         -  static dieletric constant
!      tauw         -  relaxation time of water
!      s            -  salinity
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!----------------------------------------------------------------------------------

  use type_kinds, only: fp_kind

  implicit none

  real(fp_kind) esof

  real(fp_kind)    f,tauw,freq,t_soil,vmc,rhob,rhos,sand,clay

  real(fp_kind)    alpha,beta,ess,rhoef,t,eswi,eswo

  complex esm,esw,es1,es2


  alpha = 0.65_fp_kind

  beta  = 1.09_fp_kind - 0.11_fp_kind*sand + 0.18_fp_kind*clay

  ess = (1.01_fp_kind + 0.44_fp_kind*rhos)**2 - 0.062_fp_kind

  rhoef = -1.645_fp_kind + 1.939_fp_kind*rhob - 0.020213_fp_kind*sand + 0.01594_fp_kind*clay

  t = t_soil - 273.0_fp_kind

  f = freq*1.0e9_fp_kind

! the permittivity at the high frequency limit

  eswi = 5.5_fp_kind

! the permittivity of free space (esof)

  esof = 8.854e-12_fp_kind

! static dieletric constant (eswo)

  eswo = 87.134_fp_kind+(-1.949e-1_fp_kind+(-1.276e-2_fp_kind+2.491e-4_fp_kind*t)*t)*t

  tauw = 1.1109e-10_fp_kind+(-3.824e-12_fp_kind+(6.938e-14_fp_kind-5.096e-16_fp_kind*t)*t)*t

  if (vmc .gt. zero) then

     es1 = cmplx(eswi, - rhoef *(rhos - rhob)/(two*pi*f*esof*rhos*vmc))

  else

     es1 = cmplx(eswi, zero)

  endif

  es2 = cmplx(eswo - eswi,zero)/cmplx(one, f*tauw)

  esw = es1 + es2

  esm = one + (ess**alpha - one)*rhob/rhos + vmc**beta*esw**alpha - vmc

  esm = esm**(one/alpha)

  if(aimag(esm) .ge. zero) esm = cmplx(real(esm), -0.0001_fp_kind)

  return

end subroutine Soil_Diel


subroutine Snow_Diel(frequency,ep_real,ep_imag,rad,frac,ep_eff)

!----------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    Snow_Diel   compute dielectric constant of snow
!
!   prgmmr: Fuzhong Weng and Banghua Yan                 org: nesdis              date: 2000-11-28
!
! abstract: compute dielectric constant of snow
!
!
! program history log:
!
! input argument list:
!
!       frequency   -  frequency (ghz)
!       ep_real     -  real part of dielectric constant of particle
!       ep_imag     -  imaginary part of dielectric constant of particle
!       rad         -  particle radiu (mm)
!       frac        -  fraction volume of snow (0.0 - 1.0)
!
! output argument list:
!
!       ep_eff      -  dielectric constant of the dense medium
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!  Copyright (C) 2005 Fuzhong Weng and Banghua Yan
!
!  This program is free software; you can redistribute it and/or modify it under the terms of the GNU
!  General Public License as published by the Free Software Foundation; either version 2 of the License,
!  or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
!  License for more details.
!
!  You should have received a copy of the GNU General Public License along with this program; if not, write
!  to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!----------------------------------------------------------------------------------

  use type_kinds, only: fp_kind

  implicit none

  real(fp_kind) ep_imag,ep_real

  real(fp_kind) frequency,rad,frac,k0,yr,yi

  complex  y,ep_r,ep_i,ep_eff,fracy

  k0 = two*pi/(300.0_fp_kind/frequency)

  yr  = (ep_real - one)/(ep_real + two)

  yi =   ep_imag/(ep_real + two)

  y = cmplx(yr, yi)

  fracy=frac*y

  ep_r = (one + two*fracy)/(one - fracy)

  ep_i = two*fracy*y*(k0*rad)**3*(one-frac)**4/((one-fracy)**2*(one+two*frac)**2)

  ep_eff = ep_r - cmplx(zero,one)*ep_i

  if (aimag(ep_eff).ge.zero) ep_eff = cmplx(real(ep_eff), -0.0001_fp_kind)

  return

end subroutine Snow_Diel


subroutine Canopy_Diel(frequency,mg,esv)

!----------------------------------------------------------------------------------
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:   canopy_diel compute the dielectric constant of the vegetation canopy
!
!   prgmmr:  Fuzhong Weng and Banghua Yan                org: nesdis              date: 2000-11-28
!
! abstract: compute the dielectric constant of the vegetation canopy geomatrical optics approximation
!
!           for vegetation canopy work for horizontal leaves
!
! program history log:
!
! input argument list:
!
!      frequency    -  frequency (ghz)
!      mg           -  gravimetric water content
!
! output argument list:
!
!      esv          -  dielectric constant of leaves
!
! remarks:
!
! references:
!
!     ulaby and el-rayer, 1987: microwave dielectric spectrum of vegetation part ii,
!           dual-dispersion model, ieee trans geosci. remote sensing, 25, 550-557
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!----------------------------------------------------------------------------------

  use type_kinds, only: fp_kind

  implicit none

  real(fp_kind)  frequency,  mg, en, vf, vb

  complex  esv, xx

  en = 1.7_fp_kind - (0.74_fp_kind - 6.16_fp_kind*mg)*mg

  vf = mg*(0.55_fp_kind*mg - 0.076_fp_kind)

  vb = 4.64_fp_kind*mg*mg/( one + 7.36_fp_kind*mg*mg)

  xx = cmplx(zero,one)

  esv = en + vf*(4.9_fp_kind + 75.0_fp_kind/(one + xx*frequency/18.0_fp_kind)-xx*(18.0_fp_kind/frequency)) + &

       vb*(2.9_fp_kind + 55.0_fp_kind/(one + sqrt(xx*frequency/0.18_fp_kind)))

  if (aimag(esv).ge.zero) esv = cmplx(real(esv), -0.0001_fp_kind)

  return

end subroutine Canopy_Diel


subroutine Reflectance(em1, em2, theta_i, theta_t, rvert, rh)

!----------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    Reflectance compute the surface reflectivity
!
!   prgmmr:                  org: nesdis              date: 2000-11-28
!
! abstract: compute the surface reflectivety using fresnel equations
!    for a rough surface having a standard deviation of height of sigma
!
! program history log:
!
! input argument list:
!      theta_i      -  incident angle (degree)
!      theta_t      -  transmitted angle (degree)
!      em1          -  dielectric constant of the medium 1
!      em2          -  dielectric constant of the medium 2
!
! output argument list:
!
!      rvert        -  reflectivity at vertical polarization
!      rh           -  reflectivity at horizontal polarization
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!----------------------------------------------------------------------------------

  use type_kinds, only: fp_kind

  implicit none

  real(fp_kind) theta_i, theta_t

  real(fp_kind) rh, rvert,cos_i,cos_t

  complex em1, em2, m1, m2, angle_i, angle_t

! compute the refractive index ratio between medium 2 and 1 using dielectric constant (n = sqrt(e))

  cos_i=cos(theta_i)

  cos_t=cos(theta_t)

  angle_i = cmplx(cos_i, zero)

  angle_t = cmplx(cos_t, zero)

  m1 = csqrt(em1)

  m2 = csqrt(em2)

  rvert = (cabs((m1*angle_t-m2*angle_i)/(m1*angle_t+m2*angle_i)))**2

  rh =    (cabs((m1*angle_i-m2*angle_t)/(m1*angle_i+m2*angle_t)))**2
!toshii change empirically
  

  return

end subroutine Reflectance

subroutine Transmittance(em1,em2,theta_i,theta_t,tv,th)

!----------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    Transmittance    calculate Transmittance
!
!   prgmmr:  Banghua Yan and Fuzhong Weng               org: nesdis              date: 2000-11-28
!
! abstract: compute Transmittance
!
! program history log:
!
! input argument list:
!
!      theta        -  local zenith angle (degree)
!      theta_i      -  incident angle (degree)
!      theta_t      -  transmitted angle (degree)
!      em1          -  dielectric constant of the medium 1
!      em2          -  dielectric constant of the medium 2
!
! output argument list:
!
!      tv           -  transmisivity at vertical polarization
!      th           -  transmisivity at horizontal polarization
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!----------------------------------------------------------------------------------

  use type_kinds, only: fp_kind

  implicit none

  real(fp_kind) theta_i, theta_t

  real(fp_kind) th, tv, rr, cos_i,cos_t

  complex em1, em2, m1, m2, angle_i, angle_t

! compute the refractive index ratio between medium 2 and 1 using dielectric constant (n = sqrt(e))

  cos_i=cos(theta_i)

  cos_t=cos(theta_t)

  angle_i = cmplx(cos_i, zero)

  angle_t = cmplx(cos_t, zero)

  m1 = csqrt(em1)

  m2 = csqrt(em2)

  rr = cabs(m2/m1)*cos_t/cos_i

  tv =  rr*(abs(two*m1*angle_i/(m1*angle_t + m2*angle_i)))**2

  th =  rr*(abs(two*m1*angle_i/(m1*angle_i + m2*angle_t)))**2

  return

end subroutine Transmittance


subroutine Roughness_Reflectance(frequency,sigma,rvert,rh)

!-------------------------------------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    rought_reflectance    calculate surface relectivity
!
!   prgmmr: Banghua Yan and Fuzhong Weng                 org: nesdis              date: 2000-11-28
!
! abstract: compute the surface reflectivety for a rough surface having a standard devoation of height of sigma
!
!
! program history log:
!
! input argument list:
!
!      frequency    -  frequency (ghz)
!
!      theta        -  local zenith angle (degree) (currently, not used here)
!
!      sigma        -  standard deviation of rough surface height
!
!                      smooth surface:0.38, medium: 1.10, rough:2.15 cm
!
!    internal variables
!
!
! output argument list:
!
!      rvert         -  reflectivity at vertical polarization
!
!      rh            -  reflectivity at horizontal polarization
!
!
!   important internal variables:
!
!      k0           -  a propagation constant or wavenumber in a free space
!
! remarks:
!
! references:
!
!   wang, j. and b. j. choudhury, 1992: passive microwave radiation from soil: examples...
!    passive microwave remote sensing of .. ed. b. j. choudhury, etal vsp.
!    also wang and choudhury (1982)
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!-------------------------------------------------------------------------------------------------------------

  use type_kinds, only: fp_kind

  implicit none

  real(fp_kind) frequency

  real(fp_kind) q, rh, rvert, rh_s, rv_s, sigma

  rh_s = 0.3_fp_kind*rh

  rv_s = 0.3_fp_kind*rvert

  q = 0.35_fp_kind*(one - exp(-0.60_fp_kind*frequency*sigma**two))

  rh = rh_s + q*(rv_s-rh_s)

  rvert = rv_s + q*(rh_s-rv_s)

  return

end subroutine Roughness_Reflectance


subroutine Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &

                                      r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,esv,esh)

!-------------------------------------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    two_stream_solution
!
!   prgmmr: Banghua Yan and Fuzhong Weng                 org: nesdis              date: 2000-11-28
!
! abstract: two stream solution
!
! REFERENCES:
!       Weng, F., B. Yan, and N. Grody, 2001: "A microwave land emissivity model", J. Geophys. Res., 106,
!                                             20, 115-20, 123
!   version: beta
!
! program history log:
!
! input argument list:
!
!      b            -  scattering layer temperature (k)         (gdas)   (not used here)
!      mu           -  cos(theta)
!      gv           -  asymmetry factor for v pol
!      gh           -  asymmetry factor for h pol
!      ssalb_v      -  single scattering albedo at v. polarization
!      ssalb_h      -  single scattering albedo at h. polarization
!      tau_v        -  optical depth at v. polarization
!      tau_h        -  optical depth at h. polarization
!      r12_v        -  reflectivity at vertical polarization   (not used here)
!      r12_h        -  reflectivity at horizontal polarization (not used here)
!      r21_v        -  reflectivity at vertical polarization
!      r21_h        -  reflectivity at horizontal polarization
!      r23_v        -  reflectivity at vertical polarization
!      r23_h        -  reflectivity at horizontal polarization
!      t21_v        -  transmisivity at vertical polarization
!      t21_h        -  transmisivity at horizontal polarization
!      t12_v        -  transmisivity at vertical polarization   (not used here)
!      t12_h        -  transmisivity at horizontal polarization (not used here)
!
! output argument list:
!
!       esh         -  emissivity for horizontal polarization
!       esv         -  emissivity for vertical polarization
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!-------------------------------------------------------------------------------------------------------------

  use type_kinds, only: fp_kind

  implicit none

  real(fp_kind) mu, gv, gh, ssalb_h, ssalb_v, tau_h,tau_v,                 &

                r21_h, r21_v, r23_h, r23_v, t21_v, t21_h, esv, esh

  real(fp_kind) alfa_v, alfa_h, kk_h, kk_v, gamma_h, gamma_v, beta_v, beta_h

  real(fp_kind) fact1,fact2

  alfa_h = sqrt((one - ssalb_h)/(one - gh*ssalb_h))

  kk_h = sqrt ((one - ssalb_h)*(one -  gh*ssalb_h))/mu

  beta_h = (one - alfa_h)/(one + alfa_h)

  gamma_h = (beta_h -r23_h)/(one-beta_h*r23_h)

  alfa_v = sqrt((one-ssalb_v)/(one - gv*ssalb_v))

  kk_v = sqrt ((one-ssalb_v)*(one - gv*ssalb_v))/mu

  beta_v = (one - alfa_v)/(one + alfa_v)

  gamma_v = (beta_v -r23_v)/(one-beta_v*r23_v)

  fact1=gamma_h*exp(-two*kk_h*tau_h)

  fact2=gamma_v*exp(-two*kk_v*tau_v)

  esh  = t21_h*(one - beta_h)*(one + fact1) /(one-beta_h*r21_h-(beta_h-r21_h)*fact1)

  esv  = t21_v*(one - beta_v)*(one + fact2) /(one-beta_v*r21_v-(beta_v-r21_v)*fact2)

  if (esh .lt. emissh_default) esh = emissh_default

  if (esv .lt. emissv_default) esv = emissv_default

  if (esh .gt. one) esh = one

  if (esv .gt. one) esv = one

  return

end subroutine Two_Stream_Solution


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!-------------------------------------------------------------------------------------------------------------
!
! NAME:
!       NESDIS_SSMIS_SNOW
!
! PURPOSE:
!       Subroutine to simulate microwave emissivity over snow conditions from SSMIS measurements
!
! REFERENCES:
!
!  (1) Yan, B., F. Weng and K.Okamoto,2004: "A microwave snow emissivity model, 8th Specialist Meeting on
!
!  Microwave Radiometry and Remote Sension Applications, 24-27 February, 2004, Rome, Italy.
!
!  (2) F. Weng, B. Yan, and N. Sun, Correction of SSMIS Radiance Anomalies,NRL/JCSDA Mini-Workshop on Preparation
!
!  for SSMIS Radiance Assimilation, 26-27 October, 2005, Monterey, CA.
!
!  (3) B. Yan, F. Weng, and Tsan Mo, Calibration of DMSP F-16 Special Sensor Microwave Imager and Sounder, submitted
!
!  to the 14th Conference on Satellite Meteorology and Oceanography, 29 JAN ~ 2 FEB, 2006, ATLANTA, GA
!
!  (4) B. Yan and F. Weng, Assessments of Special Sensor Microwave Imager and Sounder (SSMIS) Data for NOAA Operational
!
!  Applications,submitted to the 9th Specialist Meeting on Microwave Radiometry and Remote Sensing Applications,
!
!  San Juan, Puerto Rico, from 28 February to 03 March 2006.
!
! CATEGORY:
!       CRTM : Surface : MW SNOW EM
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       CALL NESDIS_SSMIS_SNOW
!
! INPUT ARGUMENTS:
!
!         Frequency                Frequency User defines
!                                  This is the "I" dimension
!                                  UNITS:      GHz
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!         User_Angle               The angle value user defines (in degree).
!                                  ** NOTE: THIS IS A MANDATORY MEMBER **
!                                  **       OF THIS STRUCTURE          **
!                                  UNITS:      Degrees
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Rank-1, (I)
!
!         TV[1:4]                  SSMIS V-POL Brightness temperatures at four frequencies.
!
!         tv(1): Vertically polarized AMSR-E brighness temperature at 19.35 GHz
!         tv(2):                                                      22.235GHz
!         tv(3):                                                      37    GHz
!         tv(4):                                                      91.655GHz
!
!         TH[1:4]                  SSMIS H-POL Brightness temperatures at four frequencies.
!
!         th(1): Horizontally polarized AMSR-E brighness temperature at 19.35 GHz
!         th(2):                                                        37    GHz
!         th(3):                                                        91.655GHz
!         th(4):                                                        150   GHz
!
!         Ts                       The surface temperature.
!                                  UNITS:      Kelvin, K
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!         Tsnow                    The snow temperature.
!                                  UNITS:      Kelvin, K
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!
! OUTPUT ARGUMENTS:
!
!         Emissivity_H:            The surface emissivity at a horizontal polarization.
!                                  ** NOTE: THIS IS A MANDATORY MEMBER **
!                                  **       OF THIS STRUCTURE          **
!                                  UNITS:      N/A
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!         Emissivity_V:            The surface emissivity at a vertical polarization.
!                                  ** NOTE: THIS IS A MANDATORY MEMBER **
!                                  **       OF THIS STRUCTURE          **
!                                  UNITS:      N/A
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Scalar
!
!
! INTERNAL ARGUMENTS:
!
!         Satellite_Angle          The angle values of SSMIS measurements (in degree).
!                                  ** NOTE: THIS IS A MANDATORY MEMBER **
!                                  **       OF THIS STRUCTURE          **
!                                  UNITS:      Degrees
!                                  TYPE:       REAL( fp_kind )
!                                  DIMENSION:  Rank-1, (I)
!
!
! CALLS:
!
!       SSMIS_Snow_TB   : Subroutine to calculate the snow microwave emissivity from SSMIS TB
!
!       SSMIS_Snow_TBTS : Subroutine to calculate the snow microwave emissivity from SSMIS TB & TS
!
!
! PROGRAM HISTORY LOG:
!   2006-01-10  yan,b -  implement the algorithm for snow emissivity
!
! SIDE EFFECTS:
!       None.
!
! RESTRICTIONS:
!       None.
!
!
! CREATION HISTORY:
!       Written by:     Banghua Yan, QSS Group Inc., Banghua.Yan@noaa.gov (10-January-2006)
!
!
!       and             Fuzhong Weng, NOAA/NESDIS/ORA, Fuzhong.Weng@noaa.gov
!
!  Copyright (C) 2005 Fuzhong Weng and Banghua Yan
!
!  This program is free software; you can redistribute it and/or modify it under the terms of the GNU
!  General Public License as published by the Free Software Foundation; either version 2 of the License,
!  or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
!  License for more details.
!
!  You should have received a copy of the GNU General Public License along with this program; if not, write
!  to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!
!------------------------------------------------------------------------------------------------------------

 subroutine NESDIS_SSMIS_SNOW(Frequency,                              & ! INPUT
                              User_Angle,                             & ! INPUT
                              tv,                                     & ! INPUT
                              th,                                     & ! INPUT
                              Ts,                                     & ! INPUT
                              Tsnow,                                  & ! INPUT
                              Emissivity_H,                           & ! OUTPUT
                              Emissivity_V)                             ! OUTPUT


use type_kinds, only: ip_kind, fp_kind


implicit none

real(fp_kind),parameter    :: Satellite_Angle = 53.2_fp_kind
integer(ip_kind),parameter :: nch = 4
integer(ip_kind)           :: ich
real(fp_kind)              :: Ts,Tsnow,Frequency,User_Angle,em_vector(2),tv(nch),th(nch)
real(fp_kind)              :: esh1,esv1,esh2,esv2,desh,desv,dem
real(fp_kind), intent(out) :: Emissivity_V,Emissivity_H

!  Initialization

   Emissivity_H = 0.82_fp_kind

   Emissivity_V = 0.85_fp_kind


do ich =1, nch

   if ( tv(ich) .le. 100.0_fp_kind .or. tv(ich) .ge. 330.0_fp_kind) return

   if ( th(ich) .le. 50.0_fp_kind .or. th(ich) .ge. 330.0_fp_kind) return

enddo


! EMISSIVITY AT SATELLITE'S MEASUREMENT ANGLE

if (Tsnow .le. 100.0_fp_kind .or. Tsnow .ge. 280.0_fp_kind) Tsnow = Ts

IF( Ts .le. 100.0_fp_kind .or. Ts .ge. 280.0_fp_kind) THEN

   call SSMIS_Snow_TB(Frequency,Satellite_Angle,tv,th,em_vector)

ELSE

  call SSMIS_Snow_TBTS(Frequency,Satellite_Angle,tv,th,Ts,em_vector)

ENDIF


! Get the emissivity angle dependence

  call NESDIS_LandEM(Satellite_Angle,Frequency,0.0_fp_kind,0.0_fp_kind,Ts,Tsnow,10.0_fp_kind,esh1,esv1)

  call NESDIS_LandEM(User_Angle,Frequency,0.0_fp_kind,0.0_fp_kind,Ts,Tsnow,10.0_fp_kind,esh2,esv2)

  desh = esh1 - esh2

  desv = esv1 - esv2

  dem = ( desh + desv ) * 0.5_fp_kind

! Emissivity at User's Angle

  Emissivity_H = em_vector(1) - dem;  Emissivity_V = em_vector(2)- dem

  if (Emissivity_H > one)         Emissivity_H = one

  if (Emissivity_V > one)         Emissivity_V = one

  if (Emissivity_H < 0.3_fp_kind) Emissivity_H = 0.3_fp_kind

  if (Emissivity_V < 0.3_fp_kind) Emissivity_V = 0.3_fp_kind


 end subroutine NESDIS_SSMIS_SNOW


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PRIVATE MODULE ROUTINES ##                      ##
!##                                                                            ##
!################################################################################
!################################################################################

 subroutine SSMIS_Snow_TB(frequency,theta,tv,th,em_vector)

!**********************************************************************************************
! Programmer:
!
!     Banghua Yan and Fuzhong Weng   ORG: NESDIS              Date: 2005-01-10
!
! Abstract:
!
!     Simulate emissivity between 5.0 and 150 GHz from AMSR-E Measurements over snow conditions
!
! Input argument list:
!
!         tv(1): Vertically polarized SSMIS brighness temperature at 19.35 GHz
!         tv(2):                                                      22.235GHz
!         tv(3):                                                      37    GHz
!         tv(4):                                                      91.655GHz
!
!         th(1): Horizontally polarized SSMIS brighness temperature at 19.35 GHz
!         th(2):                                                        37    GHz
!         th(3):                                                        91.655GHz
!         th(4):                                                        150   GHz
!
!    frequency: frequency in GHz
!
!    theta  : local zenith angle in degree  (53.2 for SSMIS)
!
! Internal argument list:
!
!   ev(1:4)   : retrieved V-POL emissivity at 19.35, 22.235, 37 and 91.655 GHz
!   eh(1:4)   : ..        H-POL               19.35, 37, 91.655 and 150 GHz
!   evin(1:5) : extended V-POL emissivity at 18.7, 23.8, 36.5, 89 and 150 GHz
!   ehin(1:5) : ..       H-POL emissivity ..
!   evext(1:7): retrieved V-POL emissivity at 6.925, 10.65, 18.7, 23.8, 36.5, 89 and 150 GHz
!   ehext(1:7): retrieved H-POL emissivity at 6.925, 10.65, 18.7, 23.8, 36.5, 89 and 150 GHz
!
! Output argument lists
!
!    em_vector(1) : horizontally polarization emissivity
!    em_vector(2) : vertically polarization emissivity
!
! Optional output argument lists:
!
!   ntype         : snow type
!
! Remarks:
!
!  Questions/comments: Please send to Fuzhong.Weng@noaa.gov and Banghua.Yan@noaa.gov
!
! Attributes:
!
!   language: f90
!
!   machine:  ibm rs/6000 sp
!
!*********************************************************************************************

  use type_kinds, only: ip_kind, fp_kind
  implicit none


  integer(ip_kind), parameter :: nch = 4, nch_ext=7, ncoe = 5
  real(fp_kind)   :: frequency,theta,em_vector(*),tv(*),th(*)
  real(fp_kind)   :: freqv(nch),freqh(nch),freq_ext(nch_ext),freq_tab(nch_ext)
  real(fp_kind)   :: ev(nch),eh(nch),evin(5),ehin(5),evext(nch_ext),ehext(nch_ext)
  real(fp_kind)   :: coev(nch*10),coeh(nch*10),coefv(nch*10),coefh(nch*10)

  integer :: ich,i,k,ntype

 data freqv/19.35, 22.235, 37.0, 91.655/
 data freqh/19.35, 37.0, 91.655, 150.0/
 data freq_ext/6.925,10.65, 19.35, 23.8, 37.0, 91.655, 150.0/

 data (coev(k),k=1,5)  /7.734714e-001, 1.655787e-003,-2.105508e-003, 1.119559e-003, 1.511662e-004/
 data (coev(k),k=11,15)/7.435058e-001,-1.869867e-003, 1.404605e-003, 1.351065e-003, 4.505454e-005/
 data (coev(k),k=21,25)/6.267738e-001,-3.987848e-004,-4.065957e-003, 6.088914e-003,-2.458872e-004/
 data (coev(k),k=31,35)/3.663826e-001, 7.859076e-003,-1.181864e-002, 2.877921e-003, 3.392604e-003/

 data (coeh(k),k=1,5)  /4.280762e-001, 2.075501e-003,-4.152281e-004, 1.788226e-003,-1.391088e-003/
 data (coeh(k),k=11,15)/3.329525e-001,-1.828298e-003, 4.149635e-003, 1.778996e-003,-1.685677e-003/
 data (coeh(k),k=21,25)/2.069542e-001,-4.345058e-004,-6.713779e-004, 7.217262e-003,-3.273140e-003/
 data (coeh(k),k=31,35)/2.593676e-001,-5.911642e-004,-1.164675e-003, 4.947175e-003,-6.785180e-004/

 freq_tab(1:nch_ext) = FREQUENCY_AMSRE(1:nch_ext)

! Initialization

  ev = 0.9

  eh = 0.9


   DO ich = 1, nch

      ev(ich) = coev(1+(ich-1)*10)

      eh(ich) = coeh(1+(ich-1)*10)

      DO i=2,ncoe

         ev(ich) = ev(ich) + coev((ich-1)*10 + i)*tv(i-1)

         eh(ich) = eh(ich) + coeh((ich-1)*10 + i)*th(i-1)

      ENDDO

    ENDDO

! Extrapolate the emissivity at references of frequencies

  evin(1) = ev(1)
  evin(2) = ev(2) + (23.8-22.235)*(ev(3)-ev(2))/(37.0-22.235)
  evin(3) = ev(3)
  evin(4) = ev(3) + (89.0-37.0)*(ev(4)-ev(3))/(91.655-37.0)
  evin(5) = ev(4)

  ehin(1) = eh(1)
  ehin(2) = eh(1) + (23.8-19.35)*(eh(2)-eh(1))/(37.0-19.35)
  ehin(3) = eh(2)
  ehin(4) = eh(2) + (89.0-37.0)*(eh(3)-eh(2))/(91.655-37.0)
  ehin(5) = eh(4)

! Extrapolate emissivity at 6.925, 10.65, and 150 GHz based upon various spectrum table

  call snowemiss_extrapolate(evin,ehin,theta,ntype,evext,ehext)

  evext(6) = ev(4)    ! 91.655 GHz
  ehext(6) = eh(4)
  if (evext(7) .gt. evext(6)) evext(7) = evext(6)

! Interpolate emissivity at a certain frequency

  do ich=1,nch_ext

     if (frequency .le. freq_ext(1)) then

         em_vector(1) = ehext(1)

         em_vector(2) = evext(1)

         exit

      endif

     if (frequency .ge. freq_ext(nch_ext)) then

         em_vector(1) = ehext(nch_ext)

         em_vector(2) = evext(nch_ext)

         exit

      endif


      if (frequency .le. freq_ext(ich)) then

           em_vector(1) = ehext(ich-1) + &

                          (frequency-freq_ext(ich-1))*(ehext(ich) - ehext(ich-1))/(freq_ext(ich)-freq_ext(ich-1))

           em_vector(2) = evext(ich-1) + &

                          (frequency-freq_ext(ich-1))*(evext(ich) - evext(ich-1))/(freq_ext(ich)-freq_ext(ich-1))

          exit

      endif

  enddo

 end subroutine  SSMIS_Snow_TB


 subroutine SSMIS_Snow_TBTS(frequency,theta,tv,th,ts,em_vector)

!**********************************************************************************************
! Programmer:
!
!     Banghua Yan and Fuzhong Weng   ORG: NESDIS              Date: 2005-01-10
!
! Abstract:
!
!     Simulate emissivity between 5.0 and 150 GHz from AMSR-E Measurements over snow conditions
!
! Input argument list:
!
!         ts   : snow surface temperature
!
!         tv(1): Vertically polarized SSMIS brighness temperature at 19.35 GHz
!         tv(2):                                                      22.235GHz
!         tv(3):                                                      37    GHz
!         tv(4):                                                      91.655GHz
!
!         th(1): Horizontally polarized SSMIS brighness temperature at 19.35 GHz
!         th(2):                                                        37    GHz
!         th(3):                                                        91.655GHz
!         th(4):                                                        150   GHz
!
!    frequency: frequency in GHz
!
!    theta  : local zenith angle in degree  (53.2 for SSMIS)
!
! Internal argument list:
!
!   ev(1:4)   : retrieved V-POL emissivity at 19.35, 22.235, 37 and 91.655 GHz
!   eh(1:4)   : ..        H-POL               19.35, 37, 91.655 and 150 GHz
!   evin(1:5) : extended V-POL emissivity at 18.7, 23.8, 36.5, 89 and 150 GHz
!   ehin(1:5) : ..       H-POL emissivity ..
!   evext(1:7): retrieved V-POL emissivity at 6.925, 10.65, 18.7, 23.8, 36.5, 89 and 150 GHz
!   ehext(1:7): retrieved H-POL emissivity at 6.925, 10.65, 18.7, 23.8, 36.5, 89 and 150 GHz
!
! Output argument lists
!
!    em_vector(1) : horizontally polarization emissivity
!    em_vector(2) : vertically polarization emissivity
!
! Optional output argument lists:
!
!   ntype         : snow type
!
! Remarks:
!
!  Questions/comments: Please send to Fuzhong.Weng@noaa.gov and Banghua.Yan@noaa.gov
!
! Attributes:
!
!   language: f90
!
!   machine:  ibm rs/6000 sp
!
!*********************************************************************************************

  use type_kinds, only: ip_kind, fp_kind
  implicit none


  integer(ip_kind), parameter :: nch = 4, nch_ext=7, ncoe = 6
  real(fp_kind)   :: frequency,ts,theta,em_vector(*),tv(*),th(*)
  real(fp_kind)   :: freqv(nch),freqh(nch),freq_ext(nch_ext),freq_tab(nch_ext)
  real(fp_kind)   :: ev(nch),eh(nch),evin(5),ehin(5),evext(nch_ext),ehext(nch_ext)
  real(fp_kind)   :: coev(nch*10),coeh(nch*10),coefv(nch*10),coefh(nch*10)

  integer :: ich,i,k,ntype

 data freqv/19.35, 22.235, 37.0, 91.655/
 data freqh/19.35, 37.0, 91.655, 150.0/
 data freq_ext/6.925,10.65, 19.35, 23.8, 37.0, 91.655, 150.0/

 data (coev(k),k=1,6)  /1.034365e+000,-1.022243e-003, 1.948826e-003, 1.296002e-004, 8.543700e-004,-2.023007e-003/
 data (coev(k),k=11,16)/1.036103e+000,-4.873319e-003, 5.951604e-003, 2.408112e-004, 8.337080e-004,-2.268835e-003/
 data (coev(k),k=21,26)/9.555490e-001,-3.773609e-003, 1.043274e-003, 4.841364e-003, 6.402949e-004,-2.549367e-003/
 data (coev(k),k=31,36)/9.260455e-001, 2.114234e-003,-3.121380e-003, 7.542429e-004, 4.901134e-003,-4.339693e-003/

 data (coeh(k),k=1,6)  /7.650207e-001, 2.147383e-003,-1.341971e-004, 1.229319e-003,-2.745743e-004,-2.094220e-003/
 data (coeh(k),k=11,16)/7.538105e-001,-1.738515e-003, 4.500655e-003, 1.080897e-003,-2.911021e-004,-2.615773e-003/
 data (coeh(k),k=21,26)/8.634035e-001,-2.944621e-004,-1.238606e-004, 6.128376e-003,-1.097898e-003,-4.080052e-003/
 data (coeh(k),k=31,36)/1.179267e+000,-3.949168e-004,-3.974248e-004, 3.421291e-003, 2.369703e-003,-5.717482e-003/



 freq_tab(1:nch_ext) = FREQUENCY_AMSRE(1:nch_ext)


! Initialization

  ev = 0.9

  eh = 0.9


   DO ich = 1, nch

      ev(ich) = coev(1+(ich-1)*10) + ts*coev(ncoe+(ich-1)*10)

      eh(ich) = coeh(1+(ich-1)*10) + ts*coeh(ncoe+(ich-1)*10)

      DO i=2,ncoe-1

         ev(ich) = ev(ich) + coev((ich-1)*10 + i)*tv(i-1)

         eh(ich) = eh(ich) + coeh((ich-1)*10 + i)*th(i-1)

      ENDDO

    ENDDO

! Extrapolate the emissivity at references of frequencies

  evin(1) = ev(1)
  evin(2) = ev(2) + (23.8-22.235)*(ev(3)-ev(2))/(37.0-22.235)
  evin(3) = ev(3)
  evin(4) = ev(3) + (89.0-37.0)*(ev(4)-ev(3))/(91.655-37.0)
  evin(5) = ev(4)

  ehin(1) = eh(1)
  ehin(2) = eh(1) + (23.8-19.35)*(eh(2)-eh(1))/(37.0-19.35)
  ehin(3) = eh(2)
  ehin(4) = eh(2) + (89.0-37.0)*(eh(3)-eh(2))/(91.655-37.0)
  ehin(5) = eh(4)

! Extrapolate emissivity at 6.925, 10.65, and 150 GHz based upon various spectrum table

  call snowemiss_extrapolate(evin,ehin,theta,ntype,evext,ehext)

  evext(6) = ev(4)    ! 91.655 GHz
  ehext(6) = eh(4)
  if (evext(7) .gt. evext(6)) evext(7) = evext(6)

! Interpolate emissivity at a certain frequency

  do ich=1,nch_ext

     if (frequency .le. freq_ext(1)) then

         em_vector(1) = ehext(1)

         em_vector(2) = evext(1)

         exit

      endif

     if (frequency .ge. freq_ext(nch_ext)) then

         em_vector(1) = ehext(nch_ext)

         em_vector(2) = evext(nch_ext)

         exit

      endif


      if (frequency .le. freq_ext(ich)) then

           em_vector(1) = ehext(ich-1) + &

                          (frequency-freq_ext(ich-1))*(ehext(ich) - ehext(ich-1))/(freq_ext(ich)-freq_ext(ich-1))

           em_vector(2) = evext(ich-1) + &

                          (frequency-freq_ext(ich-1))*(evext(ich) - evext(ich-1))/(freq_ext(ich)-freq_ext(ich-1))

          exit

      endif

  enddo



 end subroutine  SSMIS_Snow_TBTS



subroutine snowemiss_extrapolate(ev,eh,theta,ntype,evext,ehext)

!**********************************************************************************************
! Programmer:
!
!     Banghua Yan and Fuzhong Weng   ORG: NESDIS              Date: 2004-09-20
!
! Abstract:
!
!     Simulate emissivity at a given frequency based upon various sea ice emissivity look-up tables
!
! Input argument list:
!
!    ev(1): V-POL emissivity  at 18.7  GHz
!    ev(2):                      23.8  GHz
!    ev(3):                      36.5  GHz
!    ev(4):                      89    GHz

!    eh(1): H-POL emissivity at 18.7  GHz
!    eh(2):                     23.8  GHz
!    eh(3):                     36.5  GHz
!    eh(4):                     89    GHz
!
!    frequency: frequency in GHz
!
! Output argument lists
!
!    ehext[1:7] : horizontally polarization emissivity between 6.925 and 150 GHz
!    evext[1:7] : vertically polarization emissivity between 6.925 and 150 GHz
!
! Optional Output argument lists:
!
!    ntype        : snow types
!
! Remarks:
!
!  Questions/comments: Please send to Fuzhong.Weng@noaa.gov and Banghua.Yan@noaa.gov
!
! Attributes:
!
!   language: f90
!
!   machine:  ibm rs/6000 sp
!
!*********************************************************************************************


use type_kinds, only: ip_kind, fp_kind
implicit none

integer(ip_kind), parameter :: nchin=4, nch = 7,nt=16
real(fp_kind)     :: ev(*), eh(*),evext(*), ehext(*)
real(fp_kind)     :: ew_tab(nt,nch),ev_tab(nt,nch),eh_tab(nt,nch),freq(nch)
real(fp_kind)     :: emiss(nch-1),theta,angle,cons,sins
real(fp_kind)     :: delt0,delt_l,delt_h,delt_all,dmin
integer(ip_kind)  :: ich,ip,ntype

! Sixteen candidate snow emissivity spectra

  freq = FREQUENCY_AMSRE

  ew_tab(1, 1: N_FREQ_AMSRE) = WET_SNOW_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(2, 1: N_FREQ_AMSRE) = GRASS_AFTER_SNOW_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(3, 1: N_FREQ_AMSRE) = RS_SNOW_A_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(4, 1: N_FREQ_AMSRE) = POWDER_SNOW_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(5, 1: N_FREQ_AMSRE) = RS_SNOW_B_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(6, 1: N_FREQ_AMSRE) = RS_SNOW_C_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(7, 1: N_FREQ_AMSRE) = RS_SNOW_D_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(8, 1: N_FREQ_AMSRE) = THIN_CRUST_SNOW_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(9, 1: N_FREQ_AMSRE) = RS_SNOW_E_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(10, 1: N_FREQ_AMSRE) = BOTTOM_CRUST_SNOW_A_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(11, 1: N_FREQ_AMSRE) = SHALLOW_SNOW_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(12, 1: N_FREQ_AMSRE) = DEEP_SNOW_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(13, 1: N_FREQ_AMSRE) = CRUST_SNOW_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(14, 1: N_FREQ_AMSRE) = MEDIUM_SNOW_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(15, 1: N_FREQ_AMSRE) = BOTTOM_CRUST_SNOW_B_EM_AMSRE(1:N_FREQ_AMSRE)
  ew_tab(16, 1: N_FREQ_AMSRE) = THICK_CRUST_SNOW_EM_AMSRE(1:N_FREQ_AMSRE)

  ev_tab(1, 1: N_FREQ_AMSRE) = WET_SNOW_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(2, 1: N_FREQ_AMSRE) = GRASS_AFTER_SNOW_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(3, 1: N_FREQ_AMSRE) = RS_SNOW_A_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(4, 1: N_FREQ_AMSRE) = POWDER_SNOW_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(5, 1: N_FREQ_AMSRE) = RS_SNOW_B_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(6, 1: N_FREQ_AMSRE) = RS_SNOW_C_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(7, 1: N_FREQ_AMSRE) = RS_SNOW_D_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(8, 1: N_FREQ_AMSRE) = THIN_CRUST_SNOW_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(9, 1: N_FREQ_AMSRE) = RS_SNOW_E_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(10, 1: N_FREQ_AMSRE) = BOTTOM_CRUST_SNOW_A_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(11, 1: N_FREQ_AMSRE) = SHALLOW_SNOW_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(12, 1: N_FREQ_AMSRE) = DEEP_SNOW_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(13, 1: N_FREQ_AMSRE) = CRUST_SNOW_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(14, 1: N_FREQ_AMSRE) = MEDIUM_SNOW_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(15, 1: N_FREQ_AMSRE) = BOTTOM_CRUST_SNOW_B_EV_AMSRE(1:N_FREQ_AMSRE)
  ev_tab(16, 1: N_FREQ_AMSRE) = THICK_CRUST_SNOW_EV_AMSRE(1:N_FREQ_AMSRE)

  eh_tab(1, 1: N_FREQ_AMSRE) = WET_SNOW_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(2, 1: N_FREQ_AMSRE) = GRASS_AFTER_SNOW_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(3, 1: N_FREQ_AMSRE) = RS_SNOW_A_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(4, 1: N_FREQ_AMSRE) = POWDER_SNOW_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(5, 1: N_FREQ_AMSRE) = RS_SNOW_B_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(6, 1: N_FREQ_AMSRE) = RS_SNOW_C_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(7, 1: N_FREQ_AMSRE) = RS_SNOW_D_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(8, 1: N_FREQ_AMSRE) = THIN_CRUST_SNOW_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(9, 1: N_FREQ_AMSRE) = RS_SNOW_E_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(10, 1: N_FREQ_AMSRE) = BOTTOM_CRUST_SNOW_A_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(11, 1: N_FREQ_AMSRE) = SHALLOW_SNOW_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(12, 1: N_FREQ_AMSRE) = DEEP_SNOW_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(13, 1: N_FREQ_AMSRE) = CRUST_SNOW_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(14, 1: N_FREQ_AMSRE) = MEDIUM_SNOW_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(15, 1: N_FREQ_AMSRE) = BOTTOM_CRUST_SNOW_B_EH_AMSRE(1:N_FREQ_AMSRE)
  eh_tab(16, 1: N_FREQ_AMSRE) = THICK_CRUST_SNOW_EH_AMSRE(1:N_FREQ_AMSRE)


angle = theta*3.14159/180.0

cons = cos(angle)*cos(angle)

sins = sin(angle)*sin(angle)

do ich = 1, nchin

   emiss(ich) = ev(ich)*cons + eh(ich)*sins

enddo

! Find a spectrum

! INitialization

delt_l   = 0.0

delt_h   = 0.0

dmin = 0.05

delt0 = 10.0

! Initialization of ntype

  ntype = 3

  if (eh(5) .le. 0.6) ntype = 15

  if (eh(5) .le. 0.75 .and. eh(5) .gt. 0.6) ntype =  11

  if (eh(5) .le. 0.82 .and. eh(5) .gt. 0.75) ntype =  6

do ip = 1,nt

    delt_l = abs(emiss(1)-ew_tab(ip,3))

    delt_h = abs(emiss(nchin)-ew_tab(ip,nchin+2))

    delt_all = 0.0

    do ich=1,nchin

       delt_all = delt_all + abs(emiss(ich)-ew_tab(ip,ich+2))

    enddo

   if ( (delt_l .le. dmin) .and. (delt_h .le. dmin) .and. (delt_all .le. delt0) .and. &

        (abs(eh(5)-eh_tab(ip,7)) .le. dmin)) then

        ntype = ip

        delt0 = delt_all

   endif

enddo

evext(1) = ev(1) - (ev_tab(ntype,3) - ev_tab(ntype,1))
evext(2) = ev(1) - (ev_tab(ntype,3) - ev_tab(ntype,2))
evext(3:6) = ev(1:4)
evext(nch) = eh(nchin+1) + (ev(nchin) - eh(nchin))

if (evext(nch) .gt. evext(nch-1)) evext(nch) = evext(nch-1)
ev(5) = evext(nch)
ehext(1) = eh(1) - (eh_tab(ntype,3) - eh_tab(ntype,1))
ehext(2) = eh(1) - (eh_tab(ntype,3) - eh_tab(ntype,2))
ehext(3:6) = eh(1:4)
ehext(nch) = eh(nchin+1)

! quality control

  do ich =1, nch

     if (evext(ich) .gt. 1.0) evext(ich) = 1.0

     if (evext(ich) .lt. 0.3) evext(ich) = 0.3

     if (ehext(ich) .gt. 1.0) ehext(ich) = 1.0

     if (ehext(ich) .lt. 0.3) ehext(ich) = 0.3

  enddo

end subroutine snowemiss_extrapolate


END MODULE NESDIS_LandEM_Module

