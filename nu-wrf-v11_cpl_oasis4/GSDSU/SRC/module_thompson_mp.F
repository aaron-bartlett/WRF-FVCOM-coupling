!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_thompson_mp
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
!   This module deal with size, density, and effective radius of Thompson microphysics scheme. 
! 
! History: 
!  06/2016  Mei Han @NASA GSFC: Add snow PSD number desity: num.
!  09/2015  Toshi Matsui@NASA GSFC: Modified for G-SDSU.
!  09/2015  Mei Han @NASA GSFC ; Initial ada
!           
! References: 
!  Han, M., S. A. Braun, W. S. Olson, P. O. G. Persson, and J.-W. Bao (2010),
!  Application of TRMM PR and TMI measurements to assess cloud microphysical
!  schemes in theMM5 for a winter storm, J. Appl.Meteor. Climatol.,
!  49, 1129-1148.
!
!  Field, P. R., R. J. Hogan, P. R. A. Brown, A. J. Illingworth, T.W. Choularton,
!  and R. J. Cotton (2005), Parameterization of ice-particle size distribution
!  from mid-latitude stratiform cloud. Q. J. R. Meteorol. Soc., 131, 1997-2017
!  1997
!
!  Thompson G., P. R. Field, R. M. Rasmussen, and W. D. Hall (2008),
!  Explicit forecasts of winter precipitation using an improved bulk microphysics
!  scheme. Part II: Implementation of a new snow parameterization,
!  Mon. Wea. Rev., 136, 5095-5115.
!---------------------------------------------------------------------------------------------------

 private  !everything is unless declare it as public. 

!
! public subroutines
!
 public  ::  thom_re, thom_snow_n, thom_snow_n_part2

!
! define module global parameters and constants in below. 
!
real(sdsu_fps),parameter :: rhow=1000.  !density of water [kg/m3] 
real(sdsu_fps),parameter :: rhog=400. 
real(sdsu_fps),parameter :: rd=287.            !NOT USE YET
real(sdsu_fps),parameter :: eps=.622           !NOT USE YET
real(sdsu_fps),parameter :: cmin=1.e-10        !NOT USE YET
real(sdsu_fps),parameter :: pi=3.1415926536
 
real(sdsu_fps),parameter :: mu_g=0.0
real(sdsu_fps),parameter :: mu_r=0.0
real(sdsu_fps),parameter :: mu_s = 0.6357
real(sdsu_fps),parameter :: Kap0 = 490.6       
real(sdsu_fps),parameter :: Kap1 = 17.46      
real(sdsu_fps),parameter :: Lam0 = 20.78    
real(sdsu_fps),parameter :: Lam1 = 3.29   

real(sdsu_fps),parameter :: am_r=pi*rhow/6.  !toshi unit correction
real(sdsu_fps),parameter :: bm_r=3.0
real(sdsu_fps),parameter :: am_s=0.069
real(sdsu_fps),parameter :: bm_s=2.0
real(sdsu_fps),parameter :: am_g=pi*rhog/6. 
real(sdsu_fps),parameter :: bm_g=3.0

real(sdsu_fps),parameter :: D0s=200.E-6  !minimum diameter of snow (m)

real(sdsu_fps),parameter :: R1=q_min_condensate

!Constants for rain:
real(sdsu_fps),parameter :: cre2=mu_r+1.
real(sdsu_fps),parameter :: crg2=gamma(cre2)
real(sdsu_fps),parameter :: org2=1./crg2
real(sdsu_fps),parameter :: cre3=bm_r+mu_r+1.
real(sdsu_fps),parameter :: crg3=gamma(cre3)   
real(sdsu_fps),parameter :: obmr=1./bm_r

!Constants for snow:
real(sdsu_fps),parameter :: cse1=bm_s+1.
real(sdsu_fps),parameter :: oams=1./am_s

!Constants for graupel:
real(sdsu_fps),parameter :: gonv_min = 1.E4
real(sdsu_fps),parameter :: gonv_max = 5.E6
real(sdsu_fps),parameter :: oge1=1./(bm_g+1.)
real(sdsu_fps),parameter :: cgg1=gamma(bm_g+1.)
real(sdsu_fps),parameter :: ogg1=1./cgg1
real(sdsu_fps),parameter :: cge2=mu_g+1.
real(sdsu_fps),parameter :: cgg2=gamma(cge2)
real(sdsu_fps),parameter :: ogg2=1./cgg2
real(sdsu_fps),parameter :: cgg3=gamma(bm_g+mu_g+1.)
real(sdsu_fps),parameter :: obmg=1./bm_g

 real(sdsu_fps),public :: sa(10), sb(10)
 integer :: ia, ib

!..For snow moments conversions (from Field et al. 2005)
 data (sa(ia),ia=1,10)/5.065339, -0.062659, -3.032362, 0.029469, -0.000285, 0.31255, 0.000204, &
 0.003199, 0.0, -0.015952/

 data (sb(ib),ib=1,10)/0.476221, -0.015896,  0.165977, 0.007468, -0.000141, 0.060366, 0.000079, &
 0.000594, 0.0, -0.003577/

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

   subroutine thom_re (rr,nr,rs,rg,temp,re_r,re_s,re_g )
!---------------------------------------------------------------------------------------------
!
! Comments:
!  Generalized gamma distributions for rain, graupel and cloud ice.
!  N(D) = N_0 * D**mu * exp(-lamda*D);  mu=0 is exponential assumed for rain and grapuel.
!  Sum of two gamma distrib for snow (Field et al. 2005).
!  N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
!     + Kap1*(M2/M3)**mu_s * D**mu_s * exp(-M2*Lam1*D/M3)]
!  M2 and M3 are the (bm_s)th and (bm_s+1)th moments respectively
!  calculated as function of ice water content and temperature.
!  Mass power law relations: mass=am*D**bm
!  Snow density is not explicitly assumed, because it uses Field's m(D)=0.069*D^2
!  Since the calculation does not require N(D), parameters like
!  Kap0, Kap1, Lam0, and Lam1 are not used here.
!
! History: 
!  09/2015  Toshi Matsui@NASA GSFC: Modified for G-SDSU.
!  09/2015  Mei Han @NASA GSFC ; Initial code adapated from Thompson scheme in WRF_V3.1.
!           
! References: 
!
! References:
!
! Han et al. (2013, 2015 to be submitted)
!  Field, P. R., R. J. Hogan, P. R. A. Brown, A. J. Illingworth, T.W. Choularton,
!  and R. J. Cotton (2005), Parameterization of ice-particle size distribution
!  from mid-latitude stratiform cloud. Q. J. R. Meteorol. Soc., 131, 1997-2017
!  1997
!
!  Thompson G., P. R. Field, R. M. Rasmussen, and W. D. Hall (2008),
!  Explicit forecasts of winter precipitation using an improved bulk microphysics
!  scheme. Part II: Implementation of a new snow parameterization,
!  Mon. Wea. Rev., 136, 5095-5115.
!----------------------------------------------------------------------------------------------

real(sdsu_fps),intent(inout) :: rr   !rain water content  [g/m3] 
real(sdsu_fps),intent(inout) :: nr  !rain number concentration per unit volume [#/m3]
real(sdsu_fps),intent(inout) :: rs   !snow water content  [g/m3]
real(sdsu_fps),intent(inout) :: rg   !graupel water content [g/kg]
real(sdsu_fps),intent(inout) :: temp ! temperature [K]

real(sdsu_fps),intent(out) :: re_r  !effective radius [micron]
real(sdsu_fps),intent(out) :: re_s  !effective radius [micron]
real(sdsu_fps),intent(out) :: re_g  !effective radius [micron]

real(sdsu_fps) ::  rhos  !snow bulk density [ ]

!local param
real(sdsu_fps) :: N0_min, N0_exp, N0_r, N0_s, N0_g

real(sdsu_fps) :: lam_exp, lamr, lams, lamg  !slope parameters

real(sdsu_fps) :: tc0, & !
                  smob , & ! The bm_s (th) moment, i.e., the 2nd moment
                  smoc , & ! 3rd moment
                  loga , & ! 
                  a    , & !
                  b    , & !
                  xDs      ! effective diameter [m]

!
!
!
 N0_min = gonv_max

!
!Rain:
!
 if (rr <= R1 .or. nr <= 1.e-7) then !The purpose to set the min for nr is to limit
      !unreasonably large value of dbz (say, > 80. dBZ)
    rr=R1
    nr=1.0 !at this situation, the contribution from rain should not be counted.
    re_r=0.
 else
    lamr = (am_r*crg3*org2*nr/(rr*1.e-3))**obmr  ! slope [m-1]

    re_r = 1.5/lamr * 1.e+6  ![micron] <- [m]

 !!!!!!Toshi: N0_r is the intercept value, if you don't need it, just delete it:
    N0_r = nr*org2*lamr**cre2
 endif

!
!Snow:
!
 if (rs <= R1) then
     rs   = R1
     re_s = 0.
 else
     tc0 = min (-0.1, temp-273.15)
     smob = (rs*1.e-3) * oams     !It is the bm_s (th) moment, i.e., the 2nd moment  

     ! Field et al. 2005 Table 2:
     ! Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
      loga = sa(1) + sa(2)*tc0 + sa(3)*cse1 &
                   + sa(4)*tc0*cse1 + sa(5)*tc0*tc0 &
                   + sa(6)*cse1*cse1 + sa(7)*tc0*tc0*cse1 &
                   + sa(8)*tc0*cse1*cse1 + sa(9)*tc0*tc0*tc0 &
                   + sa(10)*cse1*cse1*cse1
      a = 10.0**loga
      b = sb(1)+ sb(2)*tc0 + sb(3)*cse1 + sb(4)*tc0*cse1 &
                + sb(5)*tc0*tc0 + sb(6)*cse1*cse1 &
                + sb(7)*tc0*tc0*cse1 + sb(8)*tc0*cse1*cse1 &
                + sb(9)*tc0*tc0*tc0 + sb(10)*cse1*cse1*cse1
      smoc = a * smob**b            ! smoc is 3rd moment
  
      !effective diameter, xDs= 3rd moment / 2nd moment
      !effective radius: re_s= 3rd moment / (2* 2nd moment) 
      xDs=smoc/smob                 ! In the unit of m
      re_s=smoc/(2.*smob)*1.e6      ! change unit from m to micron

      !Effective Density of snow:  (important for radar and microwave)
      if (xDs > D0s) then
          rhos=0.13/xDs    ! xDs is m, rhos is Kg/m^3 (eg, for xDs=1 mm, rhos=130 Kg/m^3)
      else
          rhos=0.1    ! when xDs is < 200 micron, it is supposed to be cloud ice, instead of snow
                             ! so, I set the snow density to be 0.0. 
                             !May need to adjust it when considering cloud ice ......
      endif                           

 endif
                                          
 !
 !Graupel:
 !
 if (rg <= R1) then
     rg = R1
     re_g = 0.
 else
    N0_exp = (gonv_max-gonv_min)*0.5E0                   &
             * tanh((0.15E-3-rg*1.e-3)/0.15E-3)                &
             + (gonv_max+gonv_min)*0.5E0
    N0_min = MIN(N0_exp, N0_min)
    N0_exp = N0_min
    lam_exp = (N0_exp*am_g*cgg1/(rg*1.e-3))**oge1
    lamg = lam_exp*(cgg3*ogg2*ogg1)**obmg

    re_g = 1.5/lamg * 1.e+6  ![micron] <- [m]

 !!!!!!Toshi: N0_g is the intercept value, if you don't need it, just delete it:
    N0_g = N0_exp/(cgg2*lam_exp)*lamg**cge2
 endif


 return
 end subroutine thom_re

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine thom_snow_n (rs,temp,d_increment,d_min,d_max,imax, thom_snow_num)
!---------------------------------------------------------------------------------------------
!
! Comments:
!  Sum of two gamma distrib for snow (Thompson et al. 2008, Field et al. 2005).
!  N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
!     + Kap1*(M2/M3)**mu_s * D**mu_s * exp(-M2*Lam1*D/M3)]
!  M2 and M3 are the (bm_s)th and (bm_s+1)th moments respectively
!  calculated as function of ice water content and temperature.
!  Mass power law relations: mass=am*D**bm
!  Snow density is not explicitly assumed, because it uses Field's m(D)=0.069*D^2
!
! History:
!  06/2016  Mei Han @NASA GSFC 
!
! References:
!  Field, P. R., R. J. Hogan, P. R. A. Brown, A. J. Illingworth, T.W. Choularton,
!  and R. J. Cotton (2005), Parameterization of ice-particle size distribution
!  from mid-latitude stratiform cloud. Q. J. R. Meteorol. Soc., 131, 1997-2017
!  1997
!
!  Thompson G., P. R. Field, R. M. Rasmussen, and W. D. Hall (2008),
!  Explicit forecasts of winter precipitation using an improved bulk microphysics
!  scheme. Part II: Implementation of a new snow parameterization,
!  Mon. Wea. Rev., 136, 5095-5115.
!
!----------------------------------------------------------------------------------------------
real(sdsu_fps),intent(in) :: rs   !snow water content  [g/m3]
real(sdsu_fps),intent(in) :: temp ! temperature [K]
real(sdsu_fps),intent(in) :: d_increment  ! diameter increment of do loop [mm] --> 10micron
real(sdsu_fps),intent(in) :: d_min  ! minimum diameter [mm]  --> 2micron
real(sdsu_fps),intent(in) :: d_max   ! maximum diameter [mm]  --> 2cm
integer,intent(in) :: imax ! diameter bin (size of loop)
real(sdsu_fps),intent(out) :: thom_snow_num(0:imax) ! number density function of n(D) for snow [1/m4]

!local param
real(sdsu_fps) :: tc0, & !
                  smob , & ! The bm_s (th) moment, i.e., the 2nd moment
                  smoc , & ! 3rd moment
                  loga , & ! 
                  a    , & !
                  b    , & !
                  mmt23, & !
                  d_m  , & !    
                  d_mm    
integer :: ibin    !looping


!
!num for snow in each bin 
!
   if (rs <= R1) then
     thom_snow_num(:)=0.
   else
     tc0 = min (-0.1, temp-273.15)
     smob = (rs*1.e-3) * oams     !It is the bm_s (th) moment, i.e., the 2nd moment  

     ! Field et al. 2005 Table 2:
     ! Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
     loga = sa(1) + sa(2)*tc0 + sa(3)*cse1 &
                  + sa(4)*tc0*cse1 + sa(5)*tc0*tc0 &
                  + sa(6)*cse1*cse1 + sa(7)*tc0*tc0*cse1 &
                  + sa(8)*tc0*cse1*cse1 + sa(9)*tc0*tc0*tc0 &
                  + sa(10)*cse1*cse1*cse1
     a = 10.0**loga
     b = sb(1)+ sb(2)*tc0 + sb(3)*cse1 + sb(4)*tc0*cse1 &
              + sb(5)*tc0*tc0 + sb(6)*cse1*cse1 &
              + sb(7)*tc0*tc0*cse1 + sb(8)*tc0*cse1*cse1 &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse1*cse1*cse1
     smoc = a * smob**b            ! smoc is 3rd moment
     mmt23=smob/smoc  !ratio of 2nd over 3rd moment

!
!size loop to get thom_snow_num as an array:
!
     do ibin=0,imax
       d_mm = d_min + d_increment * FLOAT(ibin)  !diameter [mm]   (new in V3.5)
       d_m  = d_mm * 1e-3                                 ! diameter [m]
       thom_snow_num(ibin)=mmt23**3*smob*(Kap0*exp(-mmt23*Lam0*d_m)+Kap1*(mmt23*d_m)**mu_s*exp(-mmt23*Lam1*d_m))
       !print*,thom_snow_num(ibin)
     enddo

   endif

   !stop

   return
   end subroutine thom_snow_n

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine thom_snow_n_part2 (rs,temp,imax, rad_cm, thom_snow_num)
!---------------------------------------------------------------------------------------------
!
! Comments:
!  Sum of two gamma distrib for snow (Thompson et al. 2008, Field et al. 2005).
!  N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
!     + Kap1*(M2/M3)**mu_s * D**mu_s * exp(-M2*Lam1*D/M3)]
!  M2 and M3 are the (bm_s)th and (bm_s+1)th moments respectively
!  calculated as function of ice water content and temperature.
!  Mass power law relations: mass=am*D**bm
!  Snow density is not explicitly assumed, because it uses Field's m(D)=0.069*D^2
!
! History:
!  06/2016  Mei Han @NASA GSFC 
!
! References:
!  Field, P. R., R. J. Hogan, P. R. A. Brown, A. J. Illingworth, T.W. Choularton,
!  and R. J. Cotton (2005), Parameterization of ice-particle size distribution
!  from mid-latitude stratiform cloud. Q. J. R. Meteorol. Soc., 131, 1997-2017
!  1997
!
!  Thompson G., P. R. Field, R. M. Rasmussen, and W. D. Hall (2008),
!  Explicit forecasts of winter precipitation using an improved bulk microphysics
!  scheme. Part II: Implementation of a new snow parameterization,
!  Mon. Wea. Rev., 136, 5095-5115.
!
!----------------------------------------------------------------------------------------------
real(sdsu_fps),intent(in) :: rs     ! snow water content  [g/m3]
real(sdsu_fps),intent(in) :: temp   ! temperature [K]
integer,intent(in) :: imax          ! maximum bin 
real(sdsu_fps),intent(in) :: rad_cm(imax)  ! radius bin  [cm]  -
real(sdsu_fps),intent(out) :: thom_snow_num(imax) ! number density function of n(D) for snow [1/m4]

!local param
real(sdsu_fps) :: tc0, & !
                  smob , & ! The bm_s (th) moment, i.e., the 2nd moment
                  smoc , & ! 3rd moment
                  loga , & ! 
                  a    , & !
                  b    , & !
                  mmt23, & !
                  d_m  , & !    
                  d_mm    
integer :: ibin    !looping


!
!num for snow in each bin 
!
   if (rs <= R1) then
     thom_snow_num(:)=0.
   else
     tc0 = min (-0.1, temp-273.15)
     smob = (rs*1.e-3) * oams     !It is the bm_s (th) moment, i.e., the 2nd moment  

     ! Field et al. 2005 Table 2:
     ! Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
     loga = sa(1) + sa(2)*tc0 + sa(3)*cse1 &
                  + sa(4)*tc0*cse1 + sa(5)*tc0*tc0 &
                  + sa(6)*cse1*cse1 + sa(7)*tc0*tc0*cse1 &
                  + sa(8)*tc0*cse1*cse1 + sa(9)*tc0*tc0*tc0 &
                  + sa(10)*cse1*cse1*cse1
     a = 10.0**loga
     b = sb(1)+ sb(2)*tc0 + sb(3)*cse1 + sb(4)*tc0*cse1 &
              + sb(5)*tc0*tc0 + sb(6)*cse1*cse1 &
              + sb(7)*tc0*tc0*cse1 + sb(8)*tc0*cse1*cse1 &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse1*cse1*cse1
     smoc = a * smob**b            ! smoc is 3rd moment
     mmt23=smob/smoc  !ratio of 2nd over 3rd moment

!
!size loop to get thom_snow_num as an array:
!
     do ibin=1,imax
       d_m  = rad_cm(ibin) * 2. * 1.e-2    ! diameter [m]
       thom_snow_num(ibin)=mmt23**3*smob*(Kap0*exp(-mmt23*Lam0*d_m)+Kap1*(mmt23*d_m)**mu_s*exp(-mmt23*Lam1*d_m))
       !print*,thom_snow_num(ibin)
     enddo

   endif

   return

 end subroutine thom_snow_n_part2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_thompson_mp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
