!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration and Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  derived_variables_mod
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Contains function for deriving variables based on data in wrfinput and
! wrfbdy files, including pressure.
!------------------------------------------------------------------------------

module derived_variables_mod

   ! Reset defaults
   implicit none
   private

   ! Private constants
   real, parameter :: P_0 = 1.e5        ! Reference pressure (Pa)
   real, parameter :: R_D = 287.        ! Dry air gas constant
   real, parameter :: C_P = 7.*R_D*0.5
   real, parameter :: C_V = C_P - R_D
   real, parameter :: GAMMA = C_P/C_V
   real, parameter :: T0 = 300.
   real, parameter :: R_V = 461.6
   real, parameter :: RV_OVR_RD = R_V/R_D

   ! Public routines
   public :: calc_geopotential
   public :: calc_inverse_dryair_density
   public :: calc_moist_potential_temperature
   public :: calc_full_pressure

contains
   
   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calc_full_pressure
   !
   ! DESCRIPTION:  Calculate full (dry air plus vapor) pressure.  Based on
   ! equation 2.21 of Skamarock et al. (2008).
   !
   !---------------------------------------------------------------------------

   function calc_full_pressure(theta_moist,alpha_dry) result (pressure)

      ! Arguments
      ! Moist potential temperature (K)
      real, intent(in) :: theta_moist 
      ! Inverse dry air density (m^3 kg^-1)
      real, intent(in) :: alpha_dry

      ! Return variable
      real :: pressure ! Pa

      pressure = R_D*theta_moist/P_0/alpha_dry
      pressure = pressure**GAMMA
      pressure = P_0*pressure
      
      return
   end function calc_full_pressure

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calc_moist_potential_temperature
   !
   ! DESCRIPTION:  Calculates moist potential temperature.
   !
   !---------------------------------------------------------------------------

   function calc_moist_potential_temperature(t,qvapor) result (theta_moist)

      ! Arguments
      ! Perturbation potential temperature (K)
      real,intent(in) :: t

      ! Water vapor mixing ratio (kg/kg)
      real,intent(in) :: qvapor

      ! Return variable
      ! Moist potential temperature (K)
      real :: theta_moist

      theta_moist = (t + T0)*(1. + RV_OVR_RD*qvapor)

      return
   end function calc_moist_potential_temperature

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calc_inverse_dryair_density
   !
   ! DESCRIPTION:  Calculates inverse density of dry air.  Based on equation
   ! 2.20 of Skamarock et al. (2008)
   !
   !---------------------------------------------------------------------------

   function calc_inverse_dryair_density(geopotential_above, &
        geopotential_below,rdnw,mu,mub) &
        result (alpha_dry)

      ! Arguments
      ! Total geopotential (m^2 s^2)...On w points
      real,intent(in) :: geopotential_above
      real,intent(in) :: geopotential_below

      ! Inverse d(eta) values between full (w) levels
      real,intent(in) :: rdnw
      
      ! Perturbation dry air mass of column (Pa)
      real,intent(in) :: mu
      
      ! Base-state dry air mass in column (Pa)
      real,intent(in) :: mub

      ! Return variable
      ! Inverse dry air density (m^3 kg^-1)
      real :: alpha_dry

      alpha_dry = -1 * (geopotential_above - geopotential_below) * rdnw / &
           (mu + mub)
      
      return
   end function calc_inverse_dryair_density

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calc_geopotential
   !
   ! DESCRIPTION:  Calculates total geopotential.
   !
   !---------------------------------------------------------------------------

   function calc_geopotential(ph,phb) result (geopotential)

      ! Arguments
      ! Perturbation geopotential (m^2 s^-2)
      real,intent(in) :: ph

      ! Base-state geopotential (m^2 s^-2)
      real,intent(in) :: phb

      ! Return variable
      ! Geopotential (m^2 s^-2)
      real :: geopotential

      geopotential = ph + phb

      return
   end function calc_geopotential
end module derived_variables_mod
