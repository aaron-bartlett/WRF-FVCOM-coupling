#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_lidar
 use module_simulator
 use module_mpi
 implicit none
!----------------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
! Comments:
!  Lidar simulators simulates lidar-observable attenuating backscattering coef. 
!  Optical properties are computed from RSTAR6b (opt_visir). 
!
! History:
! 02/2009  Toshi Matsui@NASA GSFC ; Initial.
!
! References: 
! Matsui, T., J. Santanello, J. J. Shi, W.-K. Tao, D. Wu, C. Peters-Lidard, E. Kemp, M. Chin, 
!   D. Starr, M. Sekiguchi, and F. Aires, (2014): Introducing multisensor satellite radiance-based 
!   evaluation for regional Earth System modeling, Journal of Geophysical Research, 119, 8450-8475, 
!   doi:10.1002/2013JD021424.
!----------------------------------------------------------------------------------------------------

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 private   ! encapsulate all variables and subourtines (non accessible)

!
! public (accessible) subroutines
!
 public :: simulator_lidar ! simulating attenuating lidar signals

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  simulator_lidar(kmax,dzm,kext,sbacktot, bcat)
 implicit none
!----------------------------------------------------------------------------------------------------
! Comments:
!  Lidar simulators simulates lidar-observable attenuating backscattering coef. 
!  Optical properties are computed from RSTAR6b (opt_visir). 
!
! History:
! 02/2009  Toshi Matsui@NASA GSFC ; Initial.
!
! References: 
!----------------------------------------------------------------------------------------------------
 integer,intent(in)  :: kmax                 !maximum vertical layers 
 real(sdsu_fpd),intent(in)  :: dzm(kmax)      !layer depths [m]
 real(sdsu_fpd),intent(in)  :: kext(kmax)     !total extinction [m-1]
 real(sdsu_fpd),intent(in)  :: sbacktot(kmax) !total backscattering [m-1 str-1]
 real(sdsu_fpd),intent(out) :: bcat(kmax)    !total attenuating backscatter coef [m-1 str-1]
 real(sdsu_fpd),parameter   :: opt_strat = 0.058 ! molecular optical depth above 20.2 km (for wavelength 550nm)
 real(sdsu_fpd) :: ms_param  !multiple scattering correction
 character(len=4) :: direction 

!
! simple multiple scattering treatment [Platt 1978]
!
  ms_param = dble(MS_Correct)  !multiple scattering correction

!
! determine direction
!
  if(ground_lidar) then
     direction = 'up  '
  else !space born
     direction = 'down'
  endif

  call single_scattering(ms_param, direction, kmax, dzm, kext , sbacktot , opt_strat,  bcat  )

  return
 end subroutine simulator_lidar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine single_scattering(ms_param, direction, kmax, dzm, kext, sbacktot, opt_strat, bcat)
 use module_simulator
 implicit none
!----------------------------------------------------------------------------------------------------
! Comments:
!  Single scattering approximation with simple multiple-scattering parameter.
!
! History:
! 12/2012  Toshi Matsui@NASA GSFC ; Debugg for multiple scattering.
! 02/2009  Toshi Matsui@NASA GSFC ; Initial.
!
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fpd),intent(in) :: ms_param
 character(len=4),intent(in) :: direction !up or down
 integer,intent(in) :: kmax           ! maximum vertical layers 
 real(sdsu_fpd),intent(in) :: dzm(kmax)      !layer depths [m]
 real(sdsu_fpd),intent(in) :: kext(kmax)     !total extinction [m-1]
 real(sdsu_fpd),intent(in) :: sbacktot(kmax) !total backscattering [m-1 str-1]
 real(sdsu_fpd),intent(in) :: opt_strat
 real(sdsu_fpd),intent(out):: bcat(kmax)  !total attenuating backscatter coef [m-1 str-1]

 integer :: k,ks,ke !vertical looping
 real(sdsu_fpd) :: tau_tot  !vertically integrated total optical depth

! calculate attenuated lidar backscattering coefficient 
  tau_tot = 0.
  bcat = DBLE(undefined) 


 do k = 1,kmax
    if(kext(k) /= DBLE(undefined) ) then
       ks = k ; exit
    endif
 enddo

 do k = kmax, 1, -1
    if(kext(k) /= DBLE(undefined) ) then
       ke = k ; exit
    endif
 enddo


direction_select: select case(direction)
case('down')  !space born

  tau_tot = opt_strat !molecular optical depth above 20.2 km
  do k = ks, ke, 1
     tau_tot = tau_tot + kext(k) * dzm(k)
     bcat(k) = sbacktot(k) * exp( -2.d0 * tau_tot * ms_param )  ![m-1 str-1]
  enddo

case('up  ')  !ground sensor

  do k = ke, ks, -1
     tau_tot = tau_tot +  kext(k) * dzm(k)
     bcat(k) = sbacktot(k) * exp( -2.d0 * tau_tot * ms_param )
  enddo

case default
     call stop_sdsu(myrank, 'MSG single_scattering: There is no such direction')
end select direction_select

 return
 end subroutine single_scattering

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end  module module_lidar

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
