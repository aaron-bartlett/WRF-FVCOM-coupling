!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

module module_p3_param
implicit none

private   !all subroutines are encapsulated.

public  :: p3_psd_init, p3_psd_param, re_cloud_p3, re_rain_p3, &
           re_ice_p3, makelut_p3_re_ice, p3_vt_liq  !accessible routine


save

!
! floating point precision
!
 integer, parameter :: fps = selected_real_kind(6, 37)      ! single precision
 integer, parameter :: fp = selected_real_kind(15, 307)    ! double precision
! integer, parameter :: fp = selected_real_kind(33, 4931)   ! quadratic precision


 real(fps),allocatable :: p3_re_ice_lut(:,:,:)  !3D look-up table of ice ire [micron]


!
! public parameter
!
 integer, parameter,public :: n_rhor    =  5  ! number of indices for i_rhor  loop           (1st; outer)
 integer, parameter,public :: n_Fr      =  4  ! number of indices for i_Fr    loop           (2nd)
 integer, parameter,public :: n_Qnorm   = 50  ! number of indices for i_Qnorm loop           (3rd)
 real,public :: Qn(n_Qnorm)  !  normalized mass (Qtot/Ntot) [kg]
 real,parameter,public :: Fr(n_Fr) = &  !rime mass fraction (Qrim/Qtot=Fr) [-]
                   (/ 0., 0.333, 0.667, 1. /)
 real,parameter,public :: rhor(n_rhor) = &  !rime density (Qrim/Brim=RHOr) [kg/m3]
                   (/ 50., 250., 450., 650., 900. /)

contains

!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

 subroutine p3_psd_init( Qnorm_in, Fr_in, RHOr_in, Dth_out, Dgr_out, Dcr_out, &
                          alpha_md_out, beta_md_out, alpha_ad_out, beta_ad_out, &
                          n0_out, lambda_out, mu_out, re_out )
 implicit none

!______________________________________________________________________________________
!
! This program creates the lookup tables for ice microphysical processes used by
! the P3 microphysics scheme.
!
!  Note:  compile with double-precision (pgf90 -r8 create_p3_lookupTable_1.f90)
!
! The lookup table generation code for WRFV4.1 is attached, written in fortran 90. 
! This generates the P3 lookup table 1, and the calculations are done over 
! 3 dimensions: Qtot/Ntot, Qrim/Qtot, and Qrim/Brim, where Qtot, Ntot, Qrim, 
! and Brim are the total ice mass mixing ratio, ice number mixing ratio, 
! rime ice mass mixing ratio, and rime volume mixing ratio. 
! Note there's actually an additional dimension in the lookup table over a scaled rain size, 
! but this doesn't matter for your purposes.
!
! The calculations are done in 3 loops corresponding to the 3 ice dimensions above, 
!  and the integers associated with these loops are "i_Qnorm", "i_Fr", and "i_rhor". 
! This is described in comments at the top of the code. Within these loops, values 
! for the normalized variables associated with dimensions are assigned, and then all 
!  the calculations are done. First, mass-dimensional and areal-dimensional relations 
!  and the critical sizes like Dth and Dgr are calculated. Because there is no closed 
!  form solution, this is done iteratively. The iterative calculation is done around lines 270-295. 
!  Once the mass-size, area-size, and critical diameters are calculated, then the size 
!  distribution parameters are calculated. These parameters are the normalized N0 (intercept), 
!  lambda (slope), and mu (shape), for gamma distributions. The parameters are obtained 
!  with a very simple root-finding technique that loops around the range of all possible 
!  values of lambda, starting around line 525. All the size distribution parameters are 
!  obtained by line 683, and you won't need any of the code below this. 
!  Keep in mind that because the calculations are based on *normalized* variables, 
!  the N0 is a normalized N0 and thus needs to be multiplied by Ntot to give the actual N0. 
!  The lambda and mu do not need to be modified.
!
!
! P3 lookup table version: v4.1
! Last modified     : 2019-03-13
! moduralized       : 2019-05-08
!______________________________________________________________________________________

 real(fp),intent(in) :: &
      Qnorm_in   , & ! Normalized mass (Qtot/Ntot)      [kg]
      Fr_in      , & ! Rime mass fraction (Qrim/Qtot)   [-]
      RHOr_in        ! Density of rimed ice (Qrim/Brim) [kg/m3 ]

! global parameter used for calculation
 real(fp),intent(out) :: &
  Dth_out, &     ! critical size separating spherical/nonspherical ice [m]
  Dgr_out , &    ! size of equal mass for graupel and unrimed ice [m]
  Dcr_out , &    ! size of equal mass for graupel and partially rimed ice [m]
  alpha_md_out(4), &! alpha for m-D relationship ( mass = alpha * D ** beta )
  beta_md_out(4),  &! beta for m-D relationship ( mass = alpha * D ** beta )
  alpha_ad_out(4), &! alpha for a-D relationshi (aera = alpha * D ** beta)
  beta_ad_out(4),  &! beta for a-D relationshi (aera = alpha * D ** beta)
  n0_out     ,  &! intercept of gamma PSD
  lambda_out ,  &! slope of gamma PSD
  mu_out         ! shape param of gamma PSD

 real(fp),intent(out) :: &
    re_out  ! effective radius [micron]


!--
! integer            :: i_rhor          ! index for rho_rime loop                [1 .. n_rhor]
! integer            :: i_Fr            ! index for rime-mass-fraction loop      [1 .. n_Fr]
! integer            :: i_Qnorm         ! index for normalized qitot loop        [1 .. n_Qnorm]
 integer            :: i_Drscale       ! index for scaled mean rain size loop   [1 .. n_Drscale]

! integer, parameter :: n_rhor    =  5  ! number of indices for i_rhor  loop           (1st; outer)
! integer, parameter :: n_Fr      =  4  ! number of indices for i_Fr    loop           (2nd)
! integer, parameter :: n_Qnorm   = 50  ! number of indices for i_Qnorm loop           (3rd)
 integer, parameter :: n_Drscale = 30  ! number of indices for scaled mean rain size  (4th; inner)
!==

 integer :: i,ii,iii,jj,kk,kkk,dumii

 real(fp) :: N,q,qdum,dum1,dum2,cs1,ds1,lam,n0,lamf,qerror,del0,c0,c1,c2,dd,sum1,sum2,       &
         sum3,sum4,xx,a0,b0,a1,b1,dum,bas1,aas1,aas2,bas2,gamma,d1,d2,delu,lamold, &
         cap,lamr,dia,amg,dv,n0dum,sum5,sum6,sum7,sum8,dg,cg,bag,aag,dcritg,      &
         Fr,csr,dsr,duml,dum3,dum4,rhodep,cgpold,m1,m2,m3,dt,mu_r,initlamr,lamv,  &
         rdumii,lammin,lammax,pi,g,p,t,rho,mu,mu_i,ds,cs,bas,aas,crp,cgp

 real(fp) :: dcrit  ! Dth: critical size separating spherical/nonspherical ice [m]
 real(fp) :: dcrits ! Dgr: size of equal mass for graupel and unrimed ice [m]
 real(fp) :: dcritr ! Dcr: size of equal mass for graupel and partially rimed ice [m]

! real :: gammq

!-- outputs from lookup table (i.e. "f1prxx" read by access_lookup_table in s/r p3_main)
 real(fp), save :: uns,ums,refl,dmm,rhomm,nagg,nrwat,qsave,nsave,vdep,    &
        eff,lsave,a_100,n_100,vdep1,nsmall,nlarge,nrwats

! real(fp), dimension(n_Qnorm,n_Drscale,n_Fr) :: qrrain,nrrain,nsrain,qsrain,ngrain
!==

 real(fp), dimension(n_Drscale) :: lamrs
 real(fp), dimension(10000)     :: fall1
 real(fp), dimension(1000)      :: fall2,fallr,num,numi
 
! real(fp), dimension(n_rhor)    :: cgp
 real(fp), dimension(150)       :: mu_r_table

 real(fp), parameter            :: Dm_max = 2000.d-6   ! max. mean ice size for lambda limiter
 real(fp), parameter            :: Dm_min =    2.d-6   ! min. mean ice size for lambda limiter

! modify for constant mu_r, set value here
 real(fp), parameter            :: mu_r_constant = 0.

 real(fp), parameter            :: thrd = 1./3.
 real(fp), parameter            :: sxth = 1./6.

 real(fp) :: smallest_value


! end of variable declaration
!-------------------------------------------------------------------------------

! set constants and parameters
 smallest_value = tiny(pi)

! assume 600 hPa, 253 K for p and T for fallspeed calcs (for reference air density)
 pi  = 3.14159        !=acos(-1.)
!pi  = 3.14159265359  !=acos(-1.)        <-- to be replaced with this one (for double precision)
 g   = 9.861                           ! gravity
 p   = 60000.                          ! air pressure (pa)
 t   = 253.15                          ! temp (K)
 rho = p/(287.15*t)                    ! air density (kg m-3)
 mu  = 1.496E-6*t**1.5/(t+120.)/rho    ! viscosity of air
 dv  = 8.794E-5*t**1.81/p              ! diffusivity of water vapor in air
 dt  = 10.                             ! time step for collection (s)

! parameters for surface roughness of ice particle used for fallspeed
! see mitchell and heymsfield 2005
 del0 = 5.83
 c0   = 0.6
 c1   = 4./(del0**2*c0**0.5)
 c2   = del0**2/4.

!--- specified mass-dimension relationship (cgs units) for unrimed crystals:

! ms = cs*D^ds
!
! for graupel:
! mg = cg*D^dg     no longer used, since bulk volume is predicted
!===

!---- Choice of m-D parameters for large unrimed ice:

! Heymsfield et al. 2006
!      ds=1.75
!      cs=0.0040157+6.06e-5*(-20.)

! sector-like branches (P1b)
!      ds=2.02
!      cs=0.00142

! bullet-rosette
!     ds=2.45
!      cs=0.00739

! side planes
!      ds=2.3
!      cs=0.00419

! radiating assemblages of plates (mitchell et al. 1990)
!      ds=2.1
!      cs=0.00239

! aggreagtes of side planes, bullets, etc. (Mitchell 1996)
!      ds=2.1
!      cs=0.0028

!-- ** note: if using one of the above (.i.e. not brown and francis, which is already in mks units),
!           then uncomment line below to convert from cgs units to mks
!      cs=cs*100.**ds/1000.
!==

! Brown and Francis (1995)
 ds = 1.9    ! m-D of beta
! cs = 0.01855 ! original (pre v2.3), based on assumption of Dmax
 cs = 0.0121 ! m-D of alpha scaled value based on assumtion of Dmean from Hogan et al. 2012, JAMC

!====

! specify m-D parameter (beta) for fully rimed ice
!  note:  cg is not constant, due to variable density
 dg = 3.


!--- projected area-diam relationship (mks units) for unrimed crystals:
!     note: projected area = aas*D^bas

! sector-like branches (P1b)
!      bas = 1.97
!      aas = 0.55*100.**bas/(100.**2)

! bullet-rosettes
!      bas = 1.57
!      aas = 0.0869*100.**bas/(100.**2)

! aggreagtes of side planes, bullets, etc.
 bas = 1.88
 aas = 0.2285*100.**bas/(100.**2)

!===

!--- projected area-diam relationship (mks units) for fully rimed ice:
!    note: projected area = aag*D^bag

! assumed non-spherical
! bag = 2.0
! aag = 0.625*100.**bag/(100.**2)

! assumed spherical:
 bag = 2.
 aag = pi*0.25
!===

 dcrit = (pi/(6.*cs)*900.)**(1./(ds-3.))   !Dth

 Dth_out = dcrit  !for memory parameter


! check to make sure projected area at dcrit not greater than than of solid sphere
! stop and give warning message if that is the case

 if (pi/4.*dcrit**2.lt.aas*dcrit**bas) then
    print*,'STOP, area > area of solid ice sphere, unrimed'
    stop
 endif

!
! open output file
!
! open(unit=1,file='./p3_lookup_table_1.dat-v4.1',status='unknown')

!.........................................................
! main loop over rime density

! alpha parameter of m-D for rimed ice
 crp =  RHOr_in *pi*sxth

! do i_rhor = 1,n_rhor  ! loop for rime density (rho_rime) n_rhor=5 currently  (crp = Qrim/Brim?)

    !print*,(pi/(4.*aas))**(1./(bas-2.))  !toshii what are these?
    !print*,(aas/aag)**(1./(bag-bas))

!------------------------------------------------------------------------

! find threshold with rimed mass added

! loop over rimed mass fraction (4 points)
! Fr below are values of rime mass fraction for the lookup table
! specific values in model are interpolated between these four points


!    do i_Fr = 1,n_Fr   ! loop for rime mass fraction, Fr=4 currently (Fr = Qrim/Qtot)

       Fr = Fr_in

!-- these lines to be replaced by Fr(i_Fr) initialization outside of loops
!       if (i_Fr.eq.1) Fr = 0.
!       if (i_Fr.eq.2) Fr = 0.333
!       if (i_Fr.eq.3) Fr = 0.667
!       if (i_Fr.eq.4) Fr = 1.
!==

!
! calculate mass-dimension relationship for partially-rimed crystals
! msr = csr*D^dsr
! and boundary of critical diameters (Dgr and Dcr)
! formula from P3 Part 1 (JAS)

! dcritr is critical size separating fully-rimed from partially-rime ice

! first guess, set cgp=crp  (alphar parameterr)
       cgp = crp

! case of no riming (Fr = 0), then we need to set dcrits and dcritr to arbitrary large values

       if (Fr == 0.) then  ! Fr == 0. 
          dcrits = 1.e+6  !Dgr
          dcritr = dcrits !Dcr
          Dgr_out = dcrits
          Dcr_out = dcritr
          csr    = cs     ! m-D of alpha
          dsr    = ds     ! m-D of beta
! case of partial riming (Fr between 0 and 1)
       elseif (  Fr > 0. .and. Fr < 1. ) then  ! Fr > 0. .and. Fr < 1.

          !
          ! iterative calculation for Dgr Dcr
          !
          loop_dth_dgr: do
             dcrits = (cs/cgp)**(1./(dg-ds))   !Dgr
             dcritr = ((1.+Fr/(1.-Fr))*cs/cgp)**(1./(dg-ds)) !Dcr
             csr    = cs*(1.+Fr/(1.-Fr))  ! m-D of alpha 
             dsr    = ds                  ! m-D of beta

             ! get mean density of vapor deposition/aggregation grown ice
             rhodep = 1./(dcritr-dcrits)*6.*cs/(pi*(ds-2.))*(dcritr**(ds-2.)-dcrits**(ds-2.))

! get density of fully-rimed ice as rime mass fraction weighted rime density plus
! density of vapor deposition/aggregation grown ice
             cgpold      = cgp
             cgp = crp*Fr+rhodep*(1.-Fr)*pi*sxth

             if (abs((cgp-cgpold)/cgp).lt.0.01) exit loop_dth_dgr

          enddo loop_dth_dgr
          Dgr_out = dcrits
          Dcr_out = dcritr

       ! case of complete riming (Fr=1.0)
       elseif(  Fr == 1.0 ) then

         ! set threshold size between partially-rimed and fully-rimed ice as arbitrary large
          dcrits = (cs/cgp)**(1./(dg-ds))  !Dgr
          dcritr = 1.e+6       ! Dcr  (here is the "arbitrary large")
          Dgr_out = dcrits
          Dcr_out = dcritr
          csr    = cgp         ! m-D of alpha 
          dsr    = dg          ! m-D of beta

       endif

!    print*,'dcrit[Dth] and dcrits[Dgr]',dcrit,dcrits



!---------------------------------------------------------------------------------------
! set up particle fallspeed arrays
! fallspeed is a function of mass dimension and projected area dimension relationships
! following mitchell and heymsfield (2005), jas

! set up array of particle fallspeed to make computationally efficient

! for high-resolution (in diameter space), ice fallspeed is stored in 'fall1' array (m/s)
! for lower-resolution (in diameter space), ice fallspeed is stored in 'fall2' array (m/s)
! rain fallspeed is stored in 'fallr' (m/s)

       dd = 2.d-6  ! m

! loop over particle size

       do jj = 1,10000

! particle size (m)
!         d1 = real(jj)*dd - 0.5*dd
          d1 = real(jj)*2.*1.d-6 - 1.d-6   ! kept like this to preserve bit-matching (increment at every 2micron)


!----- get mass-size and projected area-size relationships for given size (d1)
!          call get_mass_size
!
!toshii this is m-D relationshi requires for single-scattering calculation
!

          if (d1.le.dcrit) then  !ice crystal
             cs1  = pi*sxth*900.   !alpha (m-D)
             ds1  = 3.             !beta  (m-D)
             bas1 = 2.             !beta (a-D)
             aas1 = pi/4.          !alpha (a-D)
          else if (d1.gt.dcrit.and.d1.le.dcrits) then !dense nonspherical ice
             cs1  = cs
             ds1  = ds
             bas1 = bas
             aas1 = aas
          else if (d1.gt.dcrits.and.d1.le.dcritr) then  ! dense spherical ice (a.k.a. graupel)
             cs1  = cgp
             ds1  = dg
             bas1 = bag
             aas1 = aag
          else if (d1.gt.dcritr) then  !complete rimed particle (a.k.a. hail)
             cs1  = csr
             ds1  = dsr
             if (Fr == 0. ) then
                aas1 = aas
                bas1 = bas
             else

! for projected area, keep bas1 constant, but modify aas1 according to rimed fraction
                bas1 = bas
                dum1 = aas*d1**bas
                dum2 = aag*d1**bag
                m1   = cs1*d1**ds1
                m2   = cs*d1**ds
                m3   = cgp*d1**dg
! linearly interpolate based on particle mass
                dum3 = dum1+(m1-m2)*(dum2-dum1)/(m3-m2)
!               dum3 = (1.-Fr)*dum1+Fr*dum2
                aas1 = dum3/(d1**bas)

             endif
          endif
!=====
         !
         ! Assign cs1, ds1, bas1, aas1 x 4 regime to module parameters 
         !

!-----------------------------------------------------------
! neglect turbulent correction for aggregates for now
! correction for turbulence
!            if (d1.lt.500.d-6) then
          a0 = 0.
          b0 = 0.
!            else
!               a0=1.7e-3
!               b0=0.8
!            end if

! fall speed for particle
! best number
          xx = 2.*cs1*g*rho*d1**(ds1+2.-bas1)/(aas1*(mu*rho)**2)

! drag terms
          b1 = c1*xx**0.5/(2.*((1.+c1*xx**0.5)**0.5-1.)*(1.+c1*xx**0.5)**0.5)-a0*b0*xx** &
             b0/(c2*((1.+c1*xx**0.5)**0.5-1.)**2)

          a1 = (c2*((1.+c1*xx**0.5)**0.5-1.)**2-a0*xx**b0)/xx**b1

! velocity in terms of drag terms (toshii units? m/s?)
          fall1(jj) = a1*mu**(1.-2.*b1)*(2.*cs1*g/(rho*aas1))**b1*d1**(b1*(ds1-bas1+2.)-1.)
!print*,jj, d1*1.e6, fall1(jj)
       enddo !jj-loop

!---------------------------------------------------------------------------------
! loop around normalized Q (Qnorm)
!       do i_Qnorm = 1,n_Qnorm

! lookup table values of normalized Qnorm
! (range of mean mass diameter from ~ 1 micron to 1 cm)

          q = Qnorm_in !261.7**((i_Qnorm+10)*0.1)*1.d-18

! test values
!            N=5.e3
!            q=0.01e-3

             !print*,'&&&&&&&&&&&i_rhor',i_rhor
             !print*,'***************',i_Qnorm
             !print*,'Fr',Fr
             !print*,'q/N',q !,N     !actually N is not specifid so q/N

! initialize qerror to arbitrarily large value
             qerror = 1.e+20

!.....................................................................................
! find parameters for gamma distribution by iterative process

! size distribution for ice is assumed to be
! N(D) = n0 * D^mu_i * exp(-lam*D)

! for the given q and N, we need to find n0, mu_i, and lam

! approach for finding lambda:
! cycle through a range of lambda, find closest lambda to produce correct q

! start with lam, range of lam from 100 to 1 x 10^7 is large enough to
! cover full range over mean size from approximately 1 micron to 1 cm

             do ii = 1,9000

                lam = 1.0013**ii*100.  ! [1/m]



! get 'mu' parameter (Heymsfield 2003)
! division by 100 is to convert m-1 to cm-1
                mu_i = 0.076*(lam/100.)**0.8-2.
! make sure mu_i >= 0, otherwise size dist is infinity at D = 0
                mu_i = max(mu_i,0.d0)
! set upper limit at 6
                mu_i = min(mu_i,6.d0)

!-- for lambda limiter:

! set min lam corresponding to 2000 micron for mean size
!               dum = Dm_max+Fr*(3000.e-6)
                dum = Dm_max
                lam = max(lam,(mu_i+1.)/dum)
! set max lam corresponding to 2 micron mean size
                lam = min(lam,(mu_i+1.)/Dm_min)
! this range corresponds to range of lam of 500 to 5000000
!==

! get normalized n0
                n0 = lam**(mu_i+1.)/(gamma(mu_i+1.))

! calculate integral for each of the 4 parts of the size distribution
! check difference with respect to Qnorm (toshii- I guess this is analytic formulae of patial gamma integration)

! dum1 is integral from 0 to dcrit       (solid ice)
! dum2 is integral from dcrit to dcrits  (unrimed large ice)
! dum3 is integral from dcrits to dcritr (fully rimed ice)
! dum4 is integral from dcritr to inf    (partially rimed)

! set up m-D relationship for solid ice with D < Dcrit
                cs1  = pi*sxth*900.
                ds1  = 3.
                dum1 = lam**(-ds1-mu_i-1.)*gamma(mu_i+ds1+1.)*(1.-gammq(mu_i+ds1+1.,dcrit*lam))

                dum2 = lam**(-ds-mu_i-1.)*gamma(mu_i+ds+1.)*(gammq(mu_i+ds+1.,dcrit*lam))
                dum  = lam**(-ds-mu_i-1.)*gamma(mu_i+ds+1.)*(gammq(mu_i+ds+1.,dcrits*lam))
                dum2 = dum2-dum
                dum2 = max(dum2,0.d0)

                dum3 = lam**(-dg-mu_i-1.)*gamma(mu_i+dg+1.)*(gammq(mu_i+dg+1.,dcrits*lam))
                dum  = lam**(-dg-mu_i-1.)*gamma(mu_i+dg+1.)*(gammq(mu_i+dg+1.,dcritr*lam))
                dum3 = dum3-dum
                dum3 = max(dum3,0.d0)

                dum4 = lam**(-dsr-mu_i-1.)*gamma(mu_i+dsr+1.)*(gammq(mu_i+dsr+1.,dcritr*lam))

! sum of the integrals from the 4 regions of the size distribution
                qdum = n0*(cs1*dum1+cs*dum2+cgp*dum3+csr*dum4)

! numerical integration for test to make sure incomplete gamma function is working
!               sum1 = 0.
!               dd = 1.e-6
!               do iii=1,50000
!                  dum=real(iii)*dd
!                  if (dum.lt.dcrit) then
!                  sum1 = sum1+n0*dum**mu_i*cs1*dum**ds1*
!     1                      exp(-lam*dum)*dd
!                  else
!                  sum1 = sum1+n0*dum**mu_i*cs*dum**ds*
!     1                      exp(-lam*dum)*dd
!                  end if
!               end do
!               print*,'sum1=',sum1
!               stop

                if (ii.eq.1) then
                   qerror = abs(q-qdum)
                   lamf   = lam
                endif

! find lam with smallest difference between Qnorm and estimate of Qnorm, assign to lamf
                if (abs(q-qdum).lt.qerror) then
                   lamf   = lam
                   qerror = abs(q-qdum)
                endif

!print*,'re_rain', 1./lam * 1.d+6

             enddo !ii-loop

!stop


! check and print relative error in q to make sure it is not too large
! note: large error is possible if size bounds are exceeded!!!!!!!!!!

             !print*,'qerror (%)',qerror/q*100.

! find n0 based on final lam value
! set final lamf to 'lam' variable
! this is the value of lam with the smallest qerror
             lam = lamf    ! units [1/m]?

             lambda_out = lam  !module parameter


! recalculate mu_i based on final lam
! get 'mu' parameter (Heymsfield 2003)
! division by 100 is to convert m-1 to cm-1
             mu_i = 0.076*(lam/100.)**0.8-2.
! make sure mu_i >= 0, otherwise size dist is infinity at D = 0
             mu_i = max(mu_i,0.d0)
! set upper limit at 6
             mu_i = min(mu_i,6.d0)

             mu_out = mu_i  !module pamaeter


!            n0 = N*lam**(mu_i+1.)/(gamma(mu_i+1.))

! find n0 from lam and Qnorm
! this is done instead of finding n0 from lam and N, since N
! may need to be adjusted to constrain mean size within reasonable bounds

             dum1 = lam**(-ds1-mu_i-1.)*gamma(mu_i+ds1+1.)*(1.-gammq(mu_i+ds1+1.,dcrit*lam))

             dum2 = lam**(-ds-mu_i-1.)*gamma(mu_i+ds+1.)*(gammq(mu_i+ds+1.,dcrit*lam))
             dum  = lam**(-ds-mu_i-1.)*gamma(mu_i+ds+1.)*(gammq(mu_i+ds+1.,dcrits*lam))
             dum2 = dum2-dum

             dum3 = lam**(-dg-mu_i-1.)*gamma(mu_i+dg+1.)*(gammq(mu_i+dg+1.,dcrits*lam))
             dum  = lam**(-dg-mu_i-1.)*gamma(mu_i+dg+1.)*(gammq(mu_i+dg+1.,dcritr*lam))
             dum3 = dum3-dum

             dum4 = lam**(-dsr-mu_i-1.)*gamma(mu_i+dsr+1.)*(gammq(mu_i+dsr+1.,dcritr*lam))

             n0   = q/(cs1*dum1+cs*dum2+cgp*dum3+csr*dum4)

             n0_out = n0  !module parameter


             !print*,'lam,N0:',lam,n0
             !print*,'mu_i:',mu_i


! test final lam, N0 values
!            sum1 = 0.
!            dd = 1.e-6
!               do iii=1,50000
!                  dum=real(iii)*dd
!                  if (dum.lt.dcrit) then
!                     sum1 = sum1+n0*dum**mu_i*cs1*dum**ds1*exp(-lam*dum)*dd
!                  elseif (dum.ge.dcrit.and.dum.lt.dcrits) then
!                     sum1 = sum1+n0*dum**mu_i*cs*dum**ds*exp(-lam*dum)*dd
!                  elseif (dum.ge.dcrits.and.dum.lt.dcritr) then
!                     sum1 = sum1+n0*dum**mu_i*cg*dum**dg*exp(-lam*dum)*dd
!                  elseif (dum.ge.dcritr) then
!                     sum1 = sum1+n0*dum**mu_i*csr*dum**dsr*exp(-lam*dum)*dd
!                  endif
!               enddo
!               print*,'sum1=',sum1
!               stop

!
! toshii At this point, we have solved for all of the size distribution parameters
! we can do real size loop to integrate single scattering parameters here. 
!




!--------------------------------------------------------------------
! find max/min N to constrain mean size (i.e. lambda limiter), this is stored
! and passed to lookup table, so that N can be adjusted during
! the simulation to constrain mean size


!toshii remove below section

! limit min size to 2 micron (diameter)

             duml = (mu_i+1.)/Dm_min
             dum1 = duml**(-ds1-mu_i-1.)*gamma(mu_i+ds1+1.)*(1.-gammq(mu_i+ds1+1.,dcrit*duml))
             dum2 = duml**(-ds-mu_i-1.)*gamma(mu_i+ds+1.)*(gammq(mu_i+ds+1.,dcrit*duml))
             dum  = duml**(-ds-mu_i-1.)*gamma(mu_i+ds+1.)*(gammq(mu_i+ds+1.,dcrits*duml))
             dum2 = dum2-dum
             dum3 = duml**(-dg-mu_i-1.)*gamma(mu_i+dg+1.)*(gammq(mu_i+dg+1.,dcrits*duml))
             dum  = duml**(-dg-mu_i-1.)*gamma(mu_i+dg+1.)*(gammq(mu_i+dg+1.,dcritr*duml))
             dum3 = dum3-dum
             dum4 = duml**(-dsr-mu_i-1.)*gamma(mu_i+dsr+1.)*(gammq(mu_i+dsr+1.,dcritr*duml))

             n0dum = q/(cs1*dum1+cs*dum2+cgp*dum3+csr*dum4)

! nlarge is used by P3_MAIN to apply the lambda limiter (upper limit)
             nlarge = n0dum/(duml**(mu_i+1.)/(gamma(mu_i+1.)))

! limit max size (number-weighted) to 2 mm
!            dum  = Dm_max+Fr*(3000.e-6)
             dum  = Dm_max
             duml = (mu_i+1.)/(dum)
             dum1 = duml**(-ds1-mu_i-1.)*gamma(mu_i+ds1+1.)*(1.-gammq(mu_i+ds1+1.,dcrit*duml))
             dum2 = duml**(-ds-mu_i-1.)*gamma(mu_i+ds+1.)*(gammq(mu_i+ds+1.,dcrit*duml))
             dum  = duml**(-ds-mu_i-1.)*gamma(mu_i+ds+1.)*(gammq(mu_i+ds+1.,dcrits*duml))
             dum2 = dum2-dum
             dum3 = duml**(-dg-mu_i-1.)*gamma(mu_i+dg+1.)*(gammq(mu_i+dg+1.,dcrits*duml))
             dum  = duml**(-dg-mu_i-1.)*gamma(mu_i+dg+1.)*(gammq(mu_i+dg+1.,dcritr*duml))
             dum3 = dum3-dum
             dum4 = duml**(-dsr-mu_i-1.)*gamma(mu_i+dsr+1.)*(gammq(mu_i+dsr+1.,dcritr*duml))

             n0dum = q/(cs1*dum1+cs*dum2+cgp*dum3+csr*dum4)

! nsmall is used by P3_MAIN to apply the lambda limiter (lower limit) 
             nsmall = n0dum/(duml**(mu_i+1.)/(gamma(mu_i+1.)))

             !print*,'nmax,nmin',nsmall,nlarge

!.....................................................................................
! mass- and number-weighted mean fallspeed (m/s)
! add reflectivity
!.....................................................................................

! assume conditions for t and p as assumed above (giving rhos), then in microphysics scheme
! multiply by density correction factor (rhos/rho)^0.54, from Heymsfield et al. 2006

! fallspeed formulation from Mitchell and Heymsfield 2005

             dd = 2.d-6 ! dD for numerical integration, units of m

! initialize for numerical integration
             sum1 = 0.
             sum2 = 0.
             sum3 = 0.
             sum4 = 0.
! reflectivity
             sum5 = 0.
! mass mean size
             sum6 = 0.
! mass-weighted mean density
             sum7 = 0.

! numerically integrate over size distribution
             do ii = 1,10000

! dum = particle size
!               dum = real(ii)*dd - 0.5*dd
                dum = real(ii)*dd - 1.d-6  ! units of m

          !assign mass-size parameters (depending on size at ii)
                if (dum.le.dcrit) then
                   ds1 = 3.
                   cs1 = pi*sxth*900.
                else if (dum.gt.dcrit.and.dum.le.dcrits) then
                   ds1 = ds
                   cs1 = cs
                elseif (dum.gt.dcrits.and.dum.le.dcritr) then
                   ds1 = dg
                   cs1 = cgp
                elseif (dum.gt.dcritr) then
                   ds1 = dsr
                   cs1 = csr
                endif

! numerator of number-weighted velocity - sum1
                sum1 = sum1+fall1(ii)*dum**mu_i*exp(-lam*dum)*dd

! numerator of mass-weighted velocity - sum2
                sum2 = sum2+fall1(ii)*cs1*dum**(ds1+mu_i)*exp(-lam*dum)*dd

! total number and mass for weighting above fallspeeds
! don't need to include n0 and cs since these parameters are
! in both numerator and denominator
      !denominator of number-weighted V
                sum3 = sum3+dum**mu_i*exp(-lam*dum)*dd
      !denominator of mass-weighted V
                sum4 = sum4+cs1*dum**(ds1+mu_i)*exp(-lam*dum)*dd

! reflectivity (integral of mass moment squared)
                sum5 = sum5+n0*(cs1/917.)**2*(6./pi)**2*dum**(2.*ds1+mu_i)*exp(-lam*dum)*dd

! numerator of mass-weighted mean size
                sum6 = sum6+cs1*dum**(ds1+mu_i+1.)*exp(-lam*dum)*dd

! numerator of mass-weighted density
! particle density is defined as mass divided by volume of sphere with same D
                sum7 = sum7+(cs1*dum**ds1)**2/(pi*sxth*dum**3)*dum**mu_i*exp(-lam*dum)*dd

             enddo !ii-loop

! save mean fallspeeds for lookup table
             uns   = sum1/sum3
             ums   = sum2/sum4
             refl  = sum5
             dmm   = sum6/sum4
             rhomm = sum7/sum4




!.....................................................................................
! ice effective radius
! use definition of Francis et al. (1994), e.g., Eq. 3.11 in Fu (1996) J. Climate
!.....................................................................................
             sum1 = 0.
             sum2 = 0.
             dd   = 2.d-6

! loop over exponential size distribution
             do jj = 1,10000

! particle size 
                d1 = real(jj)*2.*1.d-6 - 1.d-6

! assinge alpha and beta values for m-D and a-D relationship for module parameters
                if (d1.le.dcrit) then
                   cs1  = pi*sxth*900.  ; alpha_md_out(1) = cs1
                   ds1  = 3.            ; beta_md_out(1)  = ds1
                   aas1 = pi*0.25       ; alpha_ad_out(1) = aas1
                   bas1 = 2.            ; beta_ad_out(1)  = bas1
                elseif (d1.gt.dcrit.and.d1.le.dcrits) then
                   cs1  = cs            ; alpha_md_out(2) = cs1
                   ds1  = ds            ; beta_md_out(2)  = ds1
                   aas1 = aas           ; alpha_ad_out(2) = aas1
                   bas1 = bas           ; beta_ad_out(2)  = bas1
                elseif (d1.gt.dcrits.and.d1.le.dcritr) then
                   cs1  = cgp           ; alpha_md_out(3) = cs1
                   ds1  = dg            ; beta_md_out(3)  = ds1
                   aas1 = aag           ; alpha_ad_out(3) = aas1
                   bas1 = bag           ; beta_ad_out(3)  = bas1
                elseif (d1.gt.dcritr) then
                   cs1  = csr           ; alpha_md_out(4) = cs1
                   ds1  = dsr           ; beta_md_out(4)  = ds1
                   if (Fr == 0.) then
                      aas1 = aas        ; alpha_ad_out(4) = aas1
                      bas1 = bas        ; beta_ad_out(4)  = bas1
                   else
! for area, keep bas1 constant, but modify aas1 according to rimed fraction
                      bas1  = bas         ; beta_ad_out(4)  = bas1
                      dum1  = aas*d1**bas
                      dum2  = aag*d1**bag
                      m1    = cs1*d1**ds1
                      m2    = cs*d1**ds
                      m3    = cgp*d1**dg
! linearly interpolate based on particle mass
                      dum3  = dum1+(m1-m2)*(dum2-dum1)/(m3-m2)
                      aas1  = dum3/(d1**bas)  ; alpha_ad_out(4) = aas1
                   endif
                endif

! n0 not included below becuase it is in both numerator and denominator
!               cs1=pi*sxth*917.
!               ds1=3.
!               aas1=pi/4.*2.
!               bas1=2.

                sum1 = sum1 + cs1*d1**ds1   *d1**mu_i* exp(-lam*d1)*dd
                sum2 = sum2 + aas1*d1**bas1 *d1**mu_i* exp(-lam*d1)*dd

!               if (d1.ge.100.e-6) then
!                  sum3 = sum3+n0*d1**mu_i*exp(-lam*d1)*dd
!                  sum4 = sum4+n0*aas1*d1**bas1*d1**mu_i*exp(-lam*d1)*dd
!               endif

             enddo !jj-loop

! calculate eff radius

! hm, calculate for eff rad for twp ice
             eff = 3.*sum1/(4.*sum2*916.7)  ! make output

             re_out = eff*1.d6  !effective radius [micron]

             !print*,'eff rad',eff, re_out
!.....................................................................................

!522  continue

         !-- this column is not used (for v2.2 and after)
         !   (kept temporarly in order to preserve order of columns in lookup_table_1)
             nsave = 1.  ! HM, set to 1 for consistency w/ old lookup table file (to verify bit matching)
         !==
             qsave = q
             lsave = lam

!toshii let's remove this lator

!       enddo !i_Qnorm-loop

!    enddo !i_Fr-loop

! enddo !i_rhor-loop (main loop over variable rime density)



 end subroutine p3_psd_init

!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

! Incomplete gamma function
! from Numerical Recipes in Fortran 77: The Art of
! Scientific Computing

 function gammq(a,x)

 real(fp) :: gammq
 real(fp),intent(in) :: a,x
 real(fp) :: gamser

! USES gcf,gser
! Returns the incomplete gamma function Q(a,x) = 1-P(a,x)

 real(fp) :: gammcf,gammser,gln

 if (x.lt.0..or.a.le.0) stop 'bad argument in gammq'
 if (x.lt.a+1.) then
    call gser(gamser,a,x,gln)
    gammq=1.-gamser
 else
    call gcf(gammcf,a,x,gln)
    gammq=gammcf
 end if

 return
 end function gammq


!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

 subroutine gser(gamser,a,x,gln)

 integer,parameter :: itmax = 100
 real(fp) :: a,gamser,gln,x

 real(fp), parameter :: eps=3.d-7
 integer :: n
 real(fp) ::  ap,del,sum,gamma

 gln = log(gamma(a))
 if (x.le.0.) then
   if (x.lt.0.) stop 'MSG subroutine gser: x < 0 in gser'
   gamser = 0.
   return
 end if

 ap=a
 sum=1./a
 del=sum
 do n=1,itmax
    ap=ap+1.
    del=del*x/ap
    sum=sum+del
    if (abs(del).lt.abs(sum)*eps) goto 1
 end do

 stop 'MSG subroutine gser:a too large, itmax too small in gser'

 1    gamser=sum*exp(-x+a*log(x)-gln)
 return
 end subroutine gser


!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

 subroutine gcf(gammcf,a,x,gln)
 real(fp) ::  a,gammcf,gln,x
 integer,parameter :: itmax=100
 real(fp),parameter :: eps=3.d-7,fpmin=1.d-30
 integer :: i
 real(fp) ::  an,b,c,d,del,h,gamma

 gln=log(gamma(a))
 b=x+1.-a
 c=1./fpmin
 d=1./b
 h=d
 do i=1,itmax
    an=-i*(i-a)
    b=b+2.
    d=an*d+b
    if(abs(d).lt.fpmin) d=fpmin
    c=b+an/c
    if(abs(c).lt.fpmin) c=fpmin
    d=1./d
    del=d*c
    h = h*del
    if(abs(del-1.).lt.eps)goto 1
 end do
 stop 'MSG subroutine gcf: a too large, itmax too small in gcf'
 1    gammcf=exp(-x+a*log(x)-gln)*h
 return
 end subroutine gcf


!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

  REAL FUNCTION gamma(X)


!D    DOUBLE PRECISION FUNCTION Dgamma(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE gamma FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE gamma
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!----------------------------------------------------------------------
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH gamma(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  gamma(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!----------------------------------------------------------------------
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
      REAL(fp) ::                  &
!D    DOUBLE PRECISION
          C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,  &
           SQRTPI/0.9189385332046727417803297E0/,                     &
           PI/3.1415926535897932384626434E0/
!D    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
!D   1     SQRTPI/0.9189385332046727417803297D0/,
!D   2     PI/3.1415926535897932384626434D0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,         &
           XINF/3.4E38/
!D    DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
!D   1     XINF/1.79D308/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
            -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,  &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!D    DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
!D   1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
!D   2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
!D   3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
!D    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
!D   1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
!D   2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
!D   3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,                     &
           -5.952379913043012E-04,7.93650793500350248E-04,                &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,   &
            5.7083835261E-03/
!D    DATA C/-1.910444077728D-03,8.4171387781295D-04,
!D   1     -5.952379913043012D-04,7.93650793500350248D-04,
!D   2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
!D   3      5.7083835261D-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
!D    CONV(I) = DBLE(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO 260 I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
  260   CONTINUE
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO 290 I=1,N
            RES=RES*Y
            Y=Y+ONE
  290     CONTINUE
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO 350 I=1,6
            SUM=SUM/YSQ+C(I)
  350     CONTINUE
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 gamma=RES
!D900 Dgamma = RES
      RETURN
! ---------- LAST LINE OF gamma ----------
      END


!--------------------------------------------------------------------------
! subroutine get_mass_size
!
! !----- get mass-size and projected area-size relationships for given size (d1)
!           if (d1.le.dcrit) then
!              cs1 = pi*sxth*900.
!              ds1 = 3.
!              bas1 = 2.
!              aas1 = pi/4.
!           else if (d1.gt.dcrit.and.d1.le.dcrits) then
!              cs1  = cs
!              ds1  = ds
!              bas1 = bas
!              aas1 = aas
!           else if (d1.gt.dcrits.and.d1.le.dcritr) then
!               cs1  = cgp
!               ds1  = dg
!               bas1 = bag
!               aas1 = aag
!           else if (d1.gt.dcritr) then
!              cs1 = csr
!              ds1 = dsr
!              if (i_Fr.eq.1) then
!                 aas1 = aas
!                 bas1 = bas
!              else
!
! ! for projected area, keep bas1 constant, but modify aas1 according to rimed fraction
!                 bas1 = bas
!                 dum1 = aas*d1**bas
!                 dum2 = aag*d1**bag
!                 m1   = cs1*d1**ds1
!                 m2   = cs*d1**ds
!                 m3   = cgp*d1**dg
! ! linearly interpolate based on particle mass
!                 dum3 = dum1+(m1-m2)*(dum2-dum1)/(m3-m2)
! !               dum3 = (1.-Fr)*dum1+Fr*dum2
!                 aas1 = dum3/(d1**bas)
!              endif
!           endif
! !=====
!
! end subroutine get_mass_size


!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

 subroutine p3_psd_param( d_m, Dth, Dgr, Dcr, alpha_md, beta_md, alpha_ad, beta_ad, n0, lambda, mu ,&
                          num_out, rhoe_out, parea_out, vt_out, p3_ice_regime )
 implicit none

 real(fp),intent(in) :: &
  d_m, &         ! particle diameter [m]
  Dth, &         ! critical size separating spherical/nonspherical ice [m]
  Dgr , &        ! size of equal mass for graupel and unrimed ice [m]
  Dcr , &        ! size of equal mass for graupel and partially rimed ice [m]
  alpha_md(4), & ! alpha for m-D relationship ( mass = alpha * D ** beta )
  beta_md(4),  & ! beta for m-D relationship ( mass = alpha * D ** beta )
  alpha_ad(4), & ! alpha for a-D relationshi (aera = alpha * D ** beta)
  beta_ad(4),  & ! beta for a-D relationshi (aera = alpha * D ** beta)
  n0       ,  &  ! intercept of gamma PSD [1/m**(4+mu)]
  lambda ,  &    ! slope of gamma PSD [1/m]
  mu             ! shape param of gamma PSD [-]

 real(fps),intent(out) :: &
   num_out,  &   !number density [m4]
   rhoe_out, &   !particle density [kg/m3]
   parea_out, &  !particle projected area [m2]
   vt_out        !terminal velocity [m/s]

 integer,intent(out) :: p3_ice_regime !ice regime in P3

 real(fp) :: &
   num,  &  !number density [m4]
   rhoe, &  !particle density [kg/m3]
   parea, & !particle projected area [m2]
   vt       !terminal velocity [m/s]


 real(fp) :: cs1, ds1, aas1, bas1, mass, volume, a1, b1, xx, axis

! assume 600 hPa, 253 K for p and T for fallspeed calcs (for reference air density)
 real(fp),parameter :: pi  = 3.141592653589793        !=acos(-1.)
 real(fp),parameter :: g   = 9.861                           ! gravity
 real(fp),parameter :: p   = 60000.                          ! air pressure (pa)
 real(fp),parameter :: t   = 253.15                          ! temp (K)
 real(fp),parameter :: rho = p/(287.15*t)                    ! air density (kg m-3)
 real(fp),parameter :: amu  = 1.496E-6*t**1.5/(t+120.)/rho    ! viscosity of air
 real(fp),parameter :: dv  = 8.794E-5*t**1.81/p              ! diffusivity of water vapor in air

! parameters for surface roughness of ice particle used for fallspeed
! see mitchell and heymsfield 2005
 real(fp),parameter :: del0 = 5.83
 real(fp),parameter :: c0   = 0.6
 real(fp),parameter :: c1   = 4./(del0**2*c0**0.5)
 real(fp),parameter :: c2   = del0**2/4.

 real(fp),parameter :: a0 = 0. ! turbulent correction term for aggregates (a0=1.7e-3)
 real(fp),parameter :: b0 = 0. ! turbulent correction term for aggregates (b0=0.8)

 
!
! check input values
!
 if( d_m <= 0.) then
   print*,'MSG p3_psd_param: strange diameter values, d_m=',d_m
   stop
 endif


!
! derive number density   
!     [1/m**(4+mu)]  [ m**mu ]         [1/m]   [m]     = [1/m4]
 num = n0          * (d_m**mu) * exp(-lambda * d_m )   ! number density [1/m4]

 num_out = real(num)  ! for output in single precision

!
! derive particle density and projected area
!
 if ( d_m <= Dth ) then  !spherical ice crystal

   cs1  = alpha_md(1) 
   ds1  = beta_md(1)
   aas1 = alpha_ad(1)
   bas1 = beta_ad(1) 
   p3_ice_regime = 1

 elseif ( d_m > Dth .and. d_m <= Dgr ) then !dedrite

   cs1  = alpha_md(2)
   ds1  = beta_md(2)
   aas1 = alpha_ad(2)
   bas1 = beta_ad(2)
   p3_ice_regime = 2

 elseif ( d_m > Dgr .and. d_m <= Dcr ) then  ! dense spherical ice (a.k.a. graupel)

   cs1  = alpha_md(3)
   ds1  = beta_md(3)
   aas1 = alpha_ad(3)
   bas1 = beta_ad(3)
   p3_ice_regime = 3

 elseif ( d_m > Dcr ) then  !partially rimed particle (snow-graupel)

   cs1  = alpha_md(4)
   ds1  = beta_md(4)
   aas1 = alpha_ad(4)
   bas1 = beta_ad(4)
   p3_ice_regime = 4

 endif

 mass = cs1*(d_m**ds1)     ! particle mass [kg]

 volume = pi*(d_m**3.)/6.  ! particle volume [m3]

 rhoe = mass / volume      ! particle density [kg/m3]

 rhoe_out = real( rhoe )   ! for output in single precision

 parea = aas1*(d_m**bas1)  ! particle projected area [m2]

 parea_out = real( parea ) ! for output in single precision


!
! fall velocity [Mitchell and Heymsfield 2005]
!

! best number
  xx = 2.*cs1*g*rho*d_m**(ds1+2.-bas1)/(aas1*(amu*rho)**2.)

! drag terms
  b1 = c1*xx**0.5/(2.*((1.+c1*xx**0.5)**0.5-1.)*(1.+c1*xx**0.5)**0.5)-a0*b0*xx** &
       b0/(c2*((1.+c1*xx**0.5)**0.5-1.)**2)

  a1 = (c2*((1.+c1*xx**0.5)**0.5-1.)**2-a0*xx**b0)/xx**b1

! fall velocity in terms of drag terms (m/s)
  vt = a1*amu**(1.-2.*b1)*(2.*cs1*g/(rho*aas1))**b1*d_m**(b1*(ds1-bas1+2.)-1.)

  vt_out = real(vt)  ! for output in single precision  

 return
 end subroutine p3_psd_param

!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
 subroutine p3_vt_liq(spc, d_m, vt)
 implicit none
 character(len=*),intent(in) :: spc  !rain or cloud
 real(fps),intent(in) :: d_m
 real(fps),intent(out) :: vt  ![m/s]

 real(fps) :: dia, amg, r_cm
 real(fps) :: pi, thrd, sxth, piov3, piov6

 pi    = 3.14159265
 thrd  = 1./3.
 sxth  = 1./6.
 piov3 = pi*thrd
 piov6 = pi*sxth

 dia = d_m    ! size bin [m]
 amg = piov6*997.*dia**3           ! mass [kg]
 amg = amg*1000.                   ! convert [kg] to [g]

 select case(trim(spc)) 
 case('rain')

 !get fallspeed as a function of size [m s-1]
 if (dia*1.e+6.le.134.43)      then
     vt = 4.5795e+3*amg**(2.*thrd)
 elseif (dia*1.e+6.lt.1511.64) then
     vt = 4.962e+1*amg**thrd
 elseif (dia*1.e+6.lt.3477.84) then
     vt = 1.732e+1*amg**sxth
 else
     vt = 9.17         
 endif

 case('cloud')

!
! cloud droplet Stoke's Law (usually works well r<30um)
!       dum          = 1.496e-6*t(i,k)**1.5/(t(i,k)+120.)  ! this is mu
!
   r_cm = d_m*100./2.  !droplet radius [cm]
   vt = 1.19*1.e+4*r_cm*r_cm  ! cloud vt [m/s]

! P3's original formulae
!       acn(i,k)     = g*rhow/(18.*dum)  ! 'a' parameter for droplet fallspeed (Stokes' law)
!                Vt_qc(i,k) = acn(i,k)*gamma(4.+bcn+mu_c(i,k))/(lamc(i,k)**bcn*gamma(mu_c(i,k)+4.))

 case default
    stop 'MSG p3_vt_liq: no such spc'
 end select

 return
 end subroutine p3_vt_liq

!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

 subroutine re_cloud_p3(qc_in,nc_in,rho_dair, mu_c, reff)
 implicit none
!
! Derive effective radius based on P3 subroutine
! mu_c is a function of total number, but lambda or effect radius. 
!
 real, intent(in)   :: qc_in ![g/m3]
 real, intent(in)   :: nc_in ![#/m3]
 real, intent(in)   :: rho_dair !air density [g/m3]
 real, intent(out)  :: mu_c  ! dispersion parameter of cloud [-]
 real, intent(out)  :: reff  ! cloud effective radius [micron]

!local variables
 real :: qc  ! mass mixing ratio [kg/kg]
 real :: nc  ! number mixing ratio [#/kg]
 real :: rho ! air density [kg/m3]
 real :: lamc           ! slope perameter [1/m]
 real :: lammin,lammax  ! minimum and maximum lambda [1/m]
 real,parameter :: pi = 3.14159
 real,parameter :: rhow = 997.   !water density [kg/m3]

!--------------------------------------------------------------------------

 qc =  qc_in / rho_dair  ![g/g] -> [kg/kg]
 nc =  nc_in / (rho_dair*1.e-3) ! [#/kg]
 rho = rho_dair*1.e-3  ! [kg/m3]

 if (qc.ge.1.e-14) then

   ! set minimum nc to prevent floating point error
    nc   = max(nc,1.e-16)
    mu_c = 0.0005714*(nc*1.e-6*rho)+0.2714
    mu_c = 1./(mu_c**2)-1.
    mu_c = max(mu_c,2.)
    mu_c = min(mu_c,15.)

    ! calculate lamc
    lamc = ( pi*rhow/6. *nc*(mu_c+3.)*(mu_c+2.)*(mu_c+1.)/qc)**(1./3.)  ![1/m]

    ! apply lambda limiters
    lammin = (mu_c+1.)*2.5e+4   ! min: 40 micron mean diameter
    lammax = (mu_c+1.)*1.e+6    ! max:  1 micron mean diameter

    !
    ! update nc according to lambda limitter
    !
    if (lamc.lt.lammin) then
        lamc = lammin
    elseif (lamc.gt.lammax) then
        lamc = lammax
    endif

    reff = 0.5*(mu_c+3.)/lamc *1.e6  ! effective radius [um]

! write( 10,FMT='(F9.3, 2x, F9.3, 2x, E9.5, 2x, E9.5)' ) reff, mu_c, qc_in, nc_in
! write(10, *)  reff, mu_c, nc*1.e-6*rho
 else

    mu_c   = 0. 
    lamc   = 0.
    reff   = 0. 
 
 endif

 end subroutine re_cloud_p3

!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

 subroutine re_rain_p3 ( qr_in, nr_in, rho_dair, reff ) 
 implicit none
! Computes and returns rain size distribution parameters

!arguments:
 real, intent(in)  :: qr_in  ! [g/m3]
 real, intent(in)  :: nr_in  ![#/m3]
 real, intent(in)  :: rho_dair !air density [g/m3]
 real, intent(out) :: reff   !air density [g/m3]

!local variables:
 real :: mu_r, lamr
 real :: lammax,lammin
 real :: qr  ! mass mixing ratio [kg/kg]
 real :: nr  ! number mixing ratio [#/kg]
 real,parameter :: pi = 3.14159
 real,parameter :: rhow = 997.   !water density [kg/m3]
!--------------------------------------------------------------------------

 qr =  qr_in / rho_dair  ! [g/g] -> [kg/kg]
 nr =  nr_in / (rho_dair*1.e-3) ! [#/kg]

 if (qr.ge.1.e-14) then

   ! use lookup table to get mu
   ! mu-lambda relationship is from Cao et al. (2008), eq. (7)

   ! find spot in lookup table
   ! (scaled N/q for lookup table parameter space_
   nr      = max(nr,1.e-16)

   mu_r = 0.  !currently mu_r is set to zero

   lamr = (pi*rhow/6.*nr*(mu_r+3.)*(mu_r+2)*(mu_r+1.)/(qr))**(1./3.)  ! recalculate slope based on mu_r
   lammax = (mu_r+1.)*1.e+5   ! check for slope
   lammin = (mu_r+1.)*1250.   ! set to small value since breakup is explicitly included (mean size 0.8 mm)

   ! apply lambda limiters for rain
   if (lamr.lt.lammin) then
       lamr = lammin
   elseif (lamr.gt.lammax) then
       lamr = lammax
   endif

   reff = 0.5*(mu_r+3.)/lamr *1.e6  ! effective radius [um]

 else

   lamr = 0.
   reff = 0. 

 endif

 return
 end subroutine re_rain_p3

!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
 
 subroutine re_ice_p3( qtot, ntot, frime_in, rhorime_in, reff )
 implicit none

 real(fps),intent(in) :: qtot   !mass conc [g/m3]
 real(fps),intent(in) :: ntot   !number con [#/m3]
 real(fps),intent(in) :: frime_in  !fraction of rimed ice [-]
 real(fps),intent(in) :: rhorime_in !density of rimed ice [kg/m3]
 real(fps),intent(out) :: reff   ! particle effective radius [micron]

!
! p3 ice output
!
 real(fps) :: &
  qnorm, &      ! normalized mass (Qtot/Ntot) [kg]
  frime  , & ! fraction of rimed ice [-]
  rhorime    ! density of rimed ice [kg/m3]


 integer :: irho_found
 integer :: ifr_found
 integer :: iqn_found

 real(fps) :: wgt1_rhor, wgt2_rhor
 real(fps) :: wgt1_fr  , wgt2_fr
 real(fps) :: wgt1_qn  , wgt2_qn

 integer            :: i_rhor          ! index for rho_rime loop                [1 .. n_rhor]
 integer            :: i_Fr            ! index for rime-mass-fraction loop      [1 .. n_Fr]
 integer            :: i_Qnorm         ! index for normalized qitot loop        [1 .. n_Qnorm]



 if(ntot == 0. .or. qtot == 0. ) then
   reff = 0. ; return
 endif

 !
 ! derive P3 PSD parameters and m-D and a-D relationship
 !
 qnorm  = qtot*1.e-3/ntot  ! [kg]

 if(qnorm > 1.e-14) then

   !
   ! density weight
   !
   rhorime = min( max(rhorime_in, 50.), 900.)
   if( rhorime == 50. ) then
      wgt1_rhor = 1. ; wgt2_rhor = 0. ; irho_found = 1
   else   
      rhor_loop: do i_rhor = 1, n_rhor-1
         wgt1_rhor = 1. ; wgt2_rhor = 0. ; irho_found = 1
         if( rhorime > rhor( i_rhor ) .and. rhorime <= rhor( i_rhor+1 ) ) then
             wgt1_rhor = ( rhor( i_rhor+1 ) - rhorime ) / ( rhor( i_rhor+1 ) - rhor( i_rhor ) )
             wgt2_rhor = 1. - wgt1_rhor 
             irho_found = i_rhor
             if( wgt1_rhor > 1.0 .or. wgt1_rhor < 0. ) stop 'wrong wgt1_rhor'
             exit rhor_loop
         endif
      enddo rhor_loop
   endif

   !
   ! fraction weight
   !
   frime = min( max(frime_in, 0.), 1.)
   if( frime == 0. ) then
      wgt1_fr = 1. ; wgt2_fr = 0. ; ifr_found = 1
   else
      fr_loop: do i_Fr = 1, n_Fr-1
         wgt1_fr = 1. ; wgt2_fr = 0. ; ifr_found = 1
         if( frime > Fr( i_Fr ) .and. frime <= Fr( i_Fr+1 ) ) then
             wgt1_fr = ( Fr( i_Fr+1 ) - frime ) / ( Fr( i_Fr+1 ) - Fr( i_Fr ) )
             wgt2_fr = 1. - wgt1_fr
             ifr_found = i_Fr
             if( wgt1_fr > 1.0 .or. wgt1_fr < 0. ) stop 'wrong wgt1_fr'
             exit fr_loop
         endif
      enddo fr_loop
   endif

   !
   ! qnorm weight
   !
   if( qnorm <=  Qn(1) ) then ! Qtot/Ntot [kg]
      wgt1_qn = 1. ; wgt2_qn = 0. ; iqn_found = 1
   elseif( qnorm > Qn(n_Qnorm) )   then
      wgt1_qn = 0. ; wgt2_qn = 1. ; iqn_found = n_Qnorm - 1
   else
      qn_loop: do i_Qnorm = 1, n_Qnorm-1
         wgt1_qn = 1. ; wgt2_qn = 0. ; iqn_found = 1
         if( qnorm > Qn( i_Qnorm ) .and. qnorm <= Qn( i_Qnorm+1 ) ) then
             wgt1_qn = ( Qn( i_Qnorm+1 ) - qnorm ) / ( Qn( i_Qnorm+1 ) - Qn( i_Qnorm ) )
             wgt2_qn = 1. - wgt1_qn
             iqn_found = i_Qnorm
             if( wgt1_qn > 1.0 .or. wgt1_qn < 0. ) stop 'wrong wgt1_qn'
             exit qn_loop
         endif
      enddo qn_loop
   endif


! print*,rhorime_in, wgt1_rhor , wgt2_rhor , irho_found
! print*,frime_in  , wgt1_fr   , wgt2_fr   , ifr_found 
! print*,qnorm     , Qn(n_Qnorm) , wgt1_qn   , wgt2_qn   , iqn_found

  !
  ! interpolate LUT 
  !
  reff = &
    ( (    p3_re_ice_lut(irho_found  , ifr_found, iqn_found) * wgt1_rhor  &
         + p3_re_ice_lut(irho_found+1, ifr_found, iqn_found) * wgt2_rhor ) * wgt1_fr &
      + &
      (    p3_re_ice_lut(irho_found  , ifr_found+1, iqn_found) * wgt1_rhor &
         + p3_re_ice_lut(irho_found+1, ifr_found+1, iqn_found) * wgt2_rhor ) * wgt2_fr ) * wgt1_qn &
      + &
    ( (    p3_re_ice_lut(irho_found  , ifr_found, iqn_found+1) * wgt1_rhor  &
         + p3_re_ice_lut(irho_found+1, ifr_found, iqn_found+1) * wgt2_rhor ) * wgt1_fr &
      + &
      (    p3_re_ice_lut(irho_found  , ifr_found+1, iqn_found+1) * wgt1_rhor &
         + p3_re_ice_lut(irho_found+1, ifr_found+1, iqn_found+1) * wgt2_rhor ) * wgt2_fr ) * wgt2_qn 

 else

   reff = 0. 

 endif



 return
 end subroutine re_ice_p3

!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
 
 subroutine makelut_p3_re_ice(masterproc)
 implicit none
 logical,intent(in) :: masterproc 
!
! 
!--
 integer            :: i_rhor          ! index for rho_rime loop                [1 .. n_rhor]
 integer            :: i_Fr            ! index for rime-mass-fraction loop      [1 .. n_Fr]
 integer            :: i_Qnorm         ! index for normalized qitot loop        [1 .. n_Qnorm]


 integer, parameter :: fp = selected_real_kind(15, 307)    ! double precision

 
!
! p3 ice output
!
 real(fp) :: &
  Dth, &         ! critical size separating spherical/nonspherical ice [m]
  Dgr , &        ! size of equal mass for graupel and unrimed ice [m]
  Dcr , &        ! size of equal mass for graupel and partially rimed ice [m]
  alpha_md(4), & ! alpha for m-D relationship ( mass = alpha * D ** beta )
  beta_md(4),  & ! beta for m-D relationship ( mass = alpha * D ** beta )
  alpha_ad(4), & ! alpha for a-D relationshi (aera = alpha * D ** beta)
  beta_ad(4),  & ! beta for a-D relationshi (aera = alpha * D ** beta)
  n0norm   ,  &  ! normalized intercept of gamma PSD [1/m**(4+mu)]
  lambda ,  &    ! slope of gamma PSD [1/m]
  mu     ,  &    ! shape param of gamma PSD [-]
  re             ! effective radius of p3 ice [micron]


 if(masterproc) print*,'MSG makelut_p3_re_ice: creating LUT for P3 ice effective radius'

!
! allocate the re_lut array
!
 if( .not. allocated( p3_re_ice_lut ) ) allocate( p3_re_ice_lut(n_rhor, n_Fr, n_Qnorm) ) 


!
! generate LUTs for P3 microphysics scheme
!
 do i_rhor = 1, n_rhor

 do i_Fr = 1, n_Fr

 do i_Qnorm = 1, n_Qnorm

    Qn( i_Qnorm )  = 261.7**(real(i_Qnorm+10)*0.1)*1.e-18 ! Qtot/Ntot [kg]

    !
    ! derive P3 PSD parameters and m-D and a-D relationship
    !
    call p3_psd_init( dble(Qn(i_Qnorm)), dble(Fr(i_Fr)), dble(rhor(i_rhor)), Dth, Dgr, Dcr, &
                       alpha_md, beta_md, alpha_ad, beta_ad, n0norm, lambda, mu , re)

    p3_re_ice_lut(i_rhor, i_Fr, i_Qnorm) = real( re ) 


 enddo ; enddo ; enddo

 end subroutine makelut_p3_re_ice

!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 

 end module module_p3_param

!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
!P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 P3 
