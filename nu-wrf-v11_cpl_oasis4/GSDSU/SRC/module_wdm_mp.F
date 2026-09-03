!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_wdm_mp
 use module_const
 implicit none

!---------------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!
!
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
! Comments:  
!   This module deal with size, density, and effective radius of WDM  microphysics scheme. 
! 
! History: 
!  07/2016  Toshi Matsui @NASA GSFC: Initial development. 
!           
!
! References: 
!  Kyo-Sun Sunny Lim and Song-You Hong, 2010: Development of an Effective Double-Moment Cloud 
!          Microphysics Scheme with Prognostic Cloud Condensation Nuclei (CCN) for Weather and 
!          Climate Models. Mon. Wea. Rev. 138, 1587-1612, doi: 10.1175/2009MWR2968.1.
!
!  Soo Ya Bae,Song-You Hong,and Kyo-Sun Sunny Lim, Coupling WRF Double-Moment 6-Class Microphysics 
!           Schemes to RRTMG Radiation Scheme in Weather Research Forecasting Model, Advances in Meteorology
!           Volume 2016 (2016), Article ID 5070154, 11 pages, http://dx.doi.org/10.1155/2016/5070154
!-----------------------------------------------------------------------

 private  !everything is unless declare it as public. 

!
! public subroutines
!
 public  :: effectRad_wdm6_0d , &
            wdm_derive_dsd

 save

!
! define module global parameters and constants in below. 
!
!..Minimum microphys values
 real, parameter :: R1 = 1.E-12  ! mass density [kg/m3]
 real, parameter :: R2 = 1.E-6   ! number density [#/m3]
 real, parameter :: pi = 3.1415926536
 real, parameter :: bm_r = 3.0
 real, parameter :: obmr = 1.0/bm_r
 real, parameter :: cdm  = 5./3.
 real, parameter :: denr = 1000.      ! water density [kg/m3]
 real, parameter :: alpha = .12       ! .122 exponen factor for n0s
 real, parameter :: n0smax = 1.e11    ! maximum n0s (t=-90C unlimited)
 real, parameter :: n0s = 2.e6        ! temperature dependent intercept parameter snow
 real, parameter :: dens  =  100.0     ! effective density of snow aggregate [kg/m3]
 real, parameter :: qmin=1.e-07

 real,parameter :: wdm_vc = 1.,  wdm_alphac = 3.  !DSD dispersion parameter [-] for cloud
 real,parameter :: wdm_vr = 2.,  wdm_alphar = 1.  !DSD dispersion parameter [-] for rain

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine wdm_derive_dsd( class, q, re_m , d_m, dsd )
 implicit none
!--------------------------------------------------------------------------------
! Purpose : 
!
! History: 
!  07/2016  Toshi Matsui@NASA GSFC ; Initial. 
!
! References: 
!  Kyo-Sun Sunny Lim and Song-You Hong, 2010: Development of an Effective Double-Moment Cloud 
!          Microphysics Scheme with Prognostic Cloud Condensation Nuclei (CCN) for Weather and 
!          Climate Models. Mon. Wea. Rev. 138, 1587-1612, doi: 10.1175/2009MWR2968.1.
!  Soo Ya Bae,Song-You Hong,and Kyo-Sun Sunny Lim, Coupling WRF Double-Moment 6-Class Microphysics 
!           Schemes to RRTMG Radiation Scheme in Weather Research Forecasting Model, Advances in Meteorology
!           Volume 2016 (2016), Article ID 5070154, 11 pages, http://dx.doi.org/10.1155/2016/5070154
!--------------------------------------------------------------------------------
 character(len=*),intent(in) :: class  !either cloud or rain 
 real,intent(in) :: q !mass mixing ratio [kg/m3]
 real,intent(in) :: re_m ! effective radius [m]
 real,intent(in) :: d_m ! particle diameter [m]
 real,intent(out) :: dsd ! number concentration density [1/m4]

 real, parameter :: const = 6./(pi*denr)
 real :: lam !slope parameters of DSD [1/m]
 real :: Nx ! number concentration [1/m3]
 real :: v, a !dispersion parameters
 real :: gam_fac, lam_fac, d_fac, exp_fac
 real :: gfac1, gfac2


 select case(trim(class))
 case('cloud')
    v = wdm_vc
    a = wdm_alphac
    call gamma_reff(3./5., gfac1)
    lam = 1. / ( 2.* gfac1 *re_m )  ! [1/m]
 case('rain')
    v = wdm_vr
    a = wdm_alphar
    lam = 2. / re_m ! [1/m]
 case default
     stop 'MSG wdm_derive_dsd: no such class'
 end select

!
! total number concentration [#/m3]
!
 call gamma_reff(v, gfac1)
 call gamma_reff(v + 3./a, gfac2)
 gam_fac = gfac1 / gfac2
 Nx = q * gam_fac * const * lam * lam * lam  ! [#/m3]

!
! number concentration for diameter [1/m4]
!
 call gamma_reff(v, gfac1)
 gam_fac = a / gfac1
 lam_fac = lam ** (a*v)
 d_fac = d_m ** (a*v-1.)
 exp_fac = exp( -1.*(lam*d_m)**a )

 dsd = Nx * gam_fac * lam_fac * d_fac * exp_fac  ! [1/m4]

 return
 end subroutine wdm_derive_dsd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine effectRad_wdm6_0d (rec_exist,rei_exist,res_exist, &
                                t, rho, qc, nc, qi, qs,  &
                                 re_qc, re_qi, re_qs)
 implicit none
!-----------------------------------------------------------------------
! Comments:  
!   Compute radiation effective radii of cloud water, ice, and snow for double-moment microphysics.
!   These are entirely consistent with microphysics assumptions, not
!   constant or otherwise ad hoc as is internal to most radiation scheme.
!
! History: 
!  07/2016: Toshi Matsui@NASA GSFC: Applied to G-SDSU, zero-dimensionalized, and changed IO with r
!                                  e-organization.  original code is written by Soo Ya Bae, KIAPS, January 2015.
!
! References: 
!  Soo Ya Bae,Song-You Hong,and Kyo-Sun Sunny Lim, Coupling WRF Double-Moment 6-Class Microphysics 
!           Schemes to RRTMG Radiation Scheme in Weather Research Forecasting Model, Advances in Meteorology
!           Volume 2016 (2016), Article ID 5070154, 11 pages, http://dx.doi.org/10.1155/2016/5070154
!-----------------------------------------------------------------------

!..Sub arguments
 logical, intent(in) :: rec_exist,rei_exist,res_exist  !io logic
 real, intent(in)::  t      !air temperature [K]
 real, intent(in)::  rho    !dry air density [kg/m3]
 real, intent(in)::  qc     !cloud mixing ratio [kg/kg]
 real, intent(in)::  nc     !cloud number concentration [#/kg]
 real, intent(in)::  qi     !ice mixing ratio [kg/kg]
 real, intent(in)::  qs     !snow mixing ratio [kg/kg]
 real, intent(inout):: re_qc  !cloud effective radius [um]
 real, intent(inout):: re_qi  !ice effective radius [um]
 real, intent(inout):: re_qs  !snow effective radius [um]

 real, parameter :: t0c = const_Kel2Cel ! = 273.15e0 [K]

!..Local variables
 integer:: i,k
 integer :: inu_c
 real :: ni  !diagnostic ice number conc [#/kg]
 real :: rqc
 real :: rnc
 real :: rqi
 real :: rni
 real :: rqs
 real :: cdm2
 real :: temp
 real :: supcol, n0sfac, lamdas, pidnc, pidn0s
 real :: diai      ! diameter of ice in m
 double precision :: lamc

!-----------------------------------------------------------------------

 cdm2 = rgmma(cdm)

!
! cloud effective radius 
!
 if( .not. rec_exist ) then

   pidnc = const_pi * denr/6.  
   rqc = max(R1, qc*rho)  ! [kg/m3]
   rnc = max(R2, nc*rho)  ! [#/m3]
   if (rqc.gt.R1 .and. rnc.gt.R2) then
     lamc = 2.*cdm2*(pidnc*nc/rqc)**obmr  !slope
     re_qc =  max(2.51, min( sngl(1.0d0/lamc)*.1e+6 , 50.) )
   else
     re_qc = 0.
   endif 

 endif

!
! ice effective radius
!
 if( .not. rei_exist ) then

   rqi = max(R1, qi*rho)       ! [kg/m3] 
   temp = (rho*max(qi,qmin))   ! ! [kg/m3]
   temp = sqrt(sqrt(temp*temp*temp))
   ni = min(max(5.38e7*temp,1.e3),1.e6)  ! diagnostic ice number conc [#/kg]
   rni= max(R2, ni*rho)                  ! number density [#/m3]
   if (rqi.gt.R1 .and. rni.gt.R2) then
     diai = 11.9*sqrt(rqi/ni)            ! dianostic ice diameter [m]
     re_qi = max(10.01e-6,min(0.75*0.163*diai,125.e-6)) * 1.e+6 ! [um]
   else
     re_qi = 0.
   endif

 endif

!
! snow effective radius
!
 if( .not. res_exist ) then

   pidn0s =  const_pi *dens*n0s
   rqs = max(R1, qs*rho)  ! [kg/m3]
   if (rqs.gt.R1) then
     supcol = t0c-t       ! [K]
     n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),1.)
     lamdas = sqrt(sqrt(pidn0s*n0sfac/rqs))
     re_qs = max(25.e-6,min(0.5*(1./lamdas),999.e-6)) * 1.e+6 ! [um]
   else
     re_qs = 0.
   endif
 endif

 return
 end subroutine effectRad_wdm6_0d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 REAL FUNCTION rgmma(x)
 IMPLICIT NONE
!-------------------------------------------------------------------
!     rgmma function:  use infinite product form
!
  REAL :: euler
 PARAMETER (euler=0.577215664901532)
 REAL :: x, y
 INTEGER :: i
 if(x.eq.1.)then
   rgmma=0.
 else
   rgmma=x*exp(euler*x)
   do i=1,10000
       y=float(i)
       rgmma=rgmma*(1.000+x/y)*exp(-x/y)
   enddo
   rgmma=1./rgmma
 endif
 END FUNCTION rgmma

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_wdm_mp

