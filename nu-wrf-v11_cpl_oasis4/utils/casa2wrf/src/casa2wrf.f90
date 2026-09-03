!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  casa2wrf
!
! DESCRIPTION:
! Process CASA global climatological CO2 concentrations and flux emissions
! and modifies wrfbdy and wrfinput files to be used by WRF-Chem.
!
! AUTHOR:
! Jossy Jacob 
! 12/16 : Major mods (C. Cruz)
!
! Notes on some settings:
!
! flux_only = 0 : for processing CO2 concentration and flux data
!             1 : for processing CO2 flux emission data only.
! flux_interpolate = 0 : use temporal interpolation with NU-WRF state.
! flux_interpolate = 1 : use temporal interpolation with NU-WRF state
!
!-------------------------------------------------------------------------------

program casa2wrf

   ! Modules
   use casa_variables_mod
   use casa_emission_variables_mod
   use interpolator_mod
   use flux_interpolator_mod
   use sorted_datetimes_mod
   use nml_mod

   implicit none

   type (interpolator)      :: conc_interp
   type (interpolator)      :: interp2
   type (flux_interpolator) :: flux_interp

   integer :: flux_only
   integer :: flux_interpolate
   integer :: ii
   type(nml) nml_c2w
   real :: flux_dt
   
   write(6,'(1x,a)') '[MAIN] START'

   nml_c2w = nml_new('namelist.casa2wrf')   
   call nml_c2w%read("wrf","flux_only",flux_only)
   call nml_c2w%read("wrf","flux_interpolate",flux_interpolate)
   call nml_c2w%read("wrf","fluxdt",flux_dt)

   if (flux_only == 0) then
      
      write(6,'(1x,a)') '[MAIN] read/interp CO2 conc'
      conc_interp = new_interpolator(nml_c2w)
      do ii = 1, NUM_CASA_VARIABLES
          call run_interpolator(conc_interp, &
              trim(casa_data_variable_names(ii)), &
              trim(casa_wrf_variable_names(ii)), &
              casa_conversion_factors(ii))
      end do
      write(6,'(1x,a)') '[MAIN] read/interp CO2 conc DONE'
      call destroy_interpolator(conc_interp)

   end if
   
   write(6,'(1x,a)') '[MAIN] read/interp CO2 fluxes'
   if (flux_interpolate == 0) then 

      interp2 = new_interpolator2(nml_c2w)
      do ii = 1, NUM_CASA_eVARIABLES
         write(6,'(1x,a)') trim(casa_data_evariable_names(ii))
         call run_interpolator2(interp2, &
              trim(casa_data_evariable_names(ii)), &
              trim(casa_wrf_evariable_names(ii)), &
              casa_econversion_factors(ii), &
              ii, flux_dt)
      end do
      call destroy_interpolator(interp2)
      
   elseif (flux_interpolate == 1) then 

      flux_interp = new_flux_interpolator(nml_c2w) 
      do ii = 1, NUM_CASA_eVARIABLES
         call run_flux_interpolator(flux_interp, &
              trim(casa_data_evariable_names(ii)), &
              trim(casa_wrf_evariable_names(ii)), &
              flux_dt)
      end do
      call destroy_flux_interpolator(flux_interp)
      
   end if 
   write(6,'(1x,a)') '[MAIN] read/interp CO2 fluxes DONE'
   
   write(6,'(1x,a)') '[MAIN] END'
   
end program casa2wrf
