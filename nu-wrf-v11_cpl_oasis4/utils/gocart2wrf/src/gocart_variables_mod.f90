!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  gocart_variables_mod
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Defines look-up tables summarizing GOCART variables, units, etc. to be
! interpolated to WRF grids.
!------------------------------------------------------------------------------

module gocart_variables_mod

   ! Change defaults
   implicit none
   private

   integer,parameter,public :: NUM_GOCART_VARIABLES = 17
   character(len=8),dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_geos5_variable_names = &
        (/'BCphobic','BCphilic','OCphobic','OCphilic', &
          'du001   ','du002   ','du003   ','du004   ', &
          'du005   ','ss002   ','ss003   ','ss004   ', &
          'ss005   ','DMS     ','MSA     ','SO2     ', &
          'SO4     '/)
   character(len=8),dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_MERRA2_variable_names = &
        (/'BCPHOBIC','BCPHILIC','OCPHOBIC','OCPHILIC', &
          'DU001   ','DU002   ','DU003   ','DU004   ', &
          'DU005   ','SS002   ','SS003   ','SS004   ', &
          'SS005   ','DMS     ','MSA     ','SO2     ', &
          'SO4     '/)
   character(len=8),dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_offline_variable_names = &
        (/'bc1     ','bc2     ','oc1     ','oc2     ', &
          'du1     ','du2     ','du3     ','du4     ', &
          'du5     ','ss1     ','ss2     ','ss3     ', &
          'ss4     ','dms     ','msa     ','so2     ', &
          'so4     '/)
   character(len=6),dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_wrf_variable_names = &
        (/'BC1   ','BC2   ','OC1   ','OC2   ', &
          'DUST_1','DUST_2','DUST_3','DUST_4', &
          'DUST_5','SEAS_1','SEAS_2','SEAS_3', &
          'SEAS_4','dms   ','msa   ','so2   ', &
          'sulf  '/)
   character(len=5),dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_geos5_variable_units = &
        (/'kg/kg','kg/kg','kg/kg','kg/kg', &
          'kg/kg','kg/kg','kg/kg','kg/kg','kg/kg', &
          'kg/kg','kg/kg','kg/kg','kg/kg', &
          'kg/kg','kg/kg','kg/kg','kg/kg' /)
   character(len=7),dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_offline_variable_units = &
        (/'mol/mol','mol/mol','mol/mol','mol/mol', &
          'kg/kg  ','kg/kg  ','kg/kg  ','kg/kg  ', &
          'kg/kg  ','kg/kg  ','kg/kg  ','kg/kg  ', &
          'kg/kg  ','mol/mol','mol/mol','mol/mol', &
          'mol/mol'/)
   character(len=12),dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_wrf_variable_units = &
        (/'ug/kg-dryair','ug/kg-dryair','ug/kg-dryair','ug/kg-dryair', &
          'ug/kg-dryair','ug/kg-dryair','ug/kg-dryair','ug/kg-dryair', &
          'ug/kg-dryair', &
          'ug/kg-dryair','ug/kg-dryair','ug/kg-dryair','ug/kg-dryair', &
          'ppmv        ','ppmv        ','ppmv        ','ppmv        '/)
   character(len=80),dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_wrf_variable_descriptions = &
        (/"Hydrophobic Black Carbon                               ", &
          "Hydrophilic Black Carbon                               ", &
          "Hydrophobic Organic Carbon                             ", &
          "Hydrophilic Organic Carbon                             ", &
          "dust size bin 1: 0.5um effective radius                ", &
          "dust size bin 2: 1.4um effective radius                ", &
          "dust size bin 3: 2.4um effective radius                ", &
          "dust size bin 4: 4.5um effective radius                ", &
          "dust size bin 5: 8.0um effective radius                ", &
          "sea-salt size bin 1: 0.3um effective radius            ", &
          "sea-salt size bin 2: 1.0um effective radius            ", &
          "sea-salt size bin 3: 3.2um effective radius            ", &
          "sea-salt size bin 4: 7.5um effective radius            ", &
          "DMS concentration                                      ", &
          "MSA concentration                                      ", &
          "SO2 concentration                                      ", &
          "SO4 concentration                                      "/)
   real,dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_geos5_conversion_factors = &
        (/1.0E+9,        1.0E+9,        1.0E+9,        1.0E+9, &
          1.0E+9,        1.0E+9,        1.0E+9,        1.0E+9, &
          1.0E+9, &
          1.0E+9,        1.0E+9,        1.0E+9,        1.0E+9, &
          0.46537842E+6, 0.30072840E+6, 0.45156299E+6, 0.30104166E+6/)
   real,dimension(NUM_GOCART_VARIABLES),parameter,public :: &
        gocart_offline_conversion_factors = &
        (/0.41522491e+9, 0.41522491e+9, 0.41522491e+9, 0.41522491e+9, &
          1.0e+9,        1.0e+9,        1.0e+9,        1.0e+9, &
          1.0e+9, &
          1.0e+9,        1.0e+9,        1.0e+9,        1.0e+9, &
          1.0e+6,        1.0e+6,        1.0e+6,        1.0e+6/)

end module gocart_variables_mod
