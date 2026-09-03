! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 module module_mp_goddard_unified

!--------------------------------------------------------------------------------------
!              = Unified Goddard Microphysics Package =
!
! NASA GSFC makes no representations about the suitability of software for any purpose.
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US
! government) nor Principal Developers (their organizations) shall be liable for any
! damages suffered by the user of this software. In addition, please do not distribute
! the software to third party.
!
!
! Comments:
!  This is an interface to drive unified Goddard microphysics packages.
!
! subroutines tree
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC : Re-unified Code (initial step)
!
! References:
!
!  Goddard single-moment 4ICE scheme
!
!  Goddard double-moment scheme
!
!----------------------------------------------------------------------------------------

!
! Previous history during 2011 - 2014
!added 
!wrf:model_layer:physics
!based on module_mp_gsfcgce_311_new_20101115_newdbz.f
!implemented conv/stra separation scheme 
!--->cal_sepa
!on 06/11/2011
!
!implemented conv/stra diagnostic variable rainncv_sepa and rainnc_sepa
!--->var_sepa
!on 11/18/2011  
!di
!
!added small fix to prevent underflow and overflow issues
!on 12/02/2011
!emk
!
!implemented 4ice writen by steve lang
!added tervrh.f_v3.0
!on 08/06/2012
!
!adding saticerh.f_v3.0
!on 08/07/2012  
!
!from gce saticerh.f_v3.08j, sgmap.f_v3.08t.clean
!consatrh.f_v3.08t, tervrh.f_v3.08t, and vqrqi.f_v3.05
!on 12/09/2014
!di
!
!
!current code
!consatrh.F_v3.08t     saticerh.F_v3.08j_TTs sgmap.F_v3.08t        tervrh.F_v3.08t       vqrqi.F_v3.04
!
!--------------------------------------------------------------------------------
 use module_wrf_error
 use module_dm

#if ( WRF_CHEM == 1)
 use     module_gocart_coupling    ! inline gocart coupling
#endif

!   use module_utility, only: wrfu_clock, wrfu_alarm
!   use module_domain, only : history_alarm, is_alarm_tstep
 use module_mp_radar

 implicit none

 private  !subroutine and module parameters within the module are not accessible in general,
          !otherwise declare it as public.

!
! public subroutines (accessible from other subroutine)
!
 public  ::                    &
    goddard_microphysics_wrf , &  ! driver for WRF
    goddard_microphysics_gce , &  ! driver for GCE
    goddard_microphysics_init     ! microphysics initialization

!
! Convert old common blocks into module parameter here (related to initialization of constant)
!

!
! run-time options
!
 character :: choice_moment             ! '1moment' or '2moment'
 logical   :: turn_on_electrification   ! if true, it will account for electric charing mechanism.
 logical   :: turn_on_aerosol           ! if true, it track in-cloud aerosol tracing process




!
! module parameters
! 

   integer, parameter :: chunk = 16

   logical, external :: wrf_dm_on_monitor

!jjs 20140117 vvvvv
!   private  ! privatize all variables/subroutines in this module excepting public parameter below
!   public :: gsfcgce_4ice_nuwrf
!jjs 20140117 ^^^^^

!jjs 1/3/2008     vvvvv

!  common /bt/
   real ::          rd1,  rd2,   al,   cp

!  common /cont/
   real ::          c38, c358, c610, c149,             &
                               c879, c172, c409,  c76,             &
                               c218, c580, c141
!  common /b3cs/
   real ::           ag,   bg,   as,   bs,             &
                                 aw,   bw,  bgh,  bgq,             &
                                bsh,  bsq,  bwh,  bwq,             &
                                 ah,   bh,  bh3,  bhh,             &
                                bhh5, bhq,  bh3_2

!  common /size/
   real ::          tnw,  tns,  tng, tnh,              &
                               roqs, roqg, roqr, roqh

!  common /rterv/
   real ::           zrc,  zgc,  zsc, zhc, vrc,        &  
                                vrc0, vrc1, vrc2, vrc3,            &
                             	 vgc,  vsc, vhc
!  common /rterv_2/ 
    real ::         zgc2, vgc2

!  common /bsnw/
   real ::              rn11a

   real ::          rn17, rn19b,                       &
                               	bnd3, rn23a,                   &
                               rn23b, rn30b,                       &
                               rn30c 

!  common /rsnw/
    real ::          alv,   alf,   als,    t0,   t00,     &
                                 avc,   afc,   asc,   esi,   rn1, rn2,     &
                                bnd2,   rn3,   rn4,   rn5,  rn50,     &
                                rn51,  rn52,  rn53,   rn6,  rn60,     &
                                rn61,  rn62,  rn63,   rn7,   rn8,     &
                                 rn9,  rn10, rn101, rn102, rn10a,     &
                               rn10b, rn10c,  rn11,  rn12,  rn14,     &
                                rn15, rn15a,  rn16, rn171, rn172,     &
                               rn17a, rn17b, rn17c,  rn18, rn18a,     &
                                rn19, rn191, rn192, rn19a,  rn20,     &
                               rn20a, rn20b,  rn30, rn30a,  rn21,     &
                               bnd21,  rn22,  rn23, rn231, rn232,     &
                                rn25,  rn31,  beta,  rn32,  rn33,     &
                               rn331, rn332,  rn34,  rn35,rnn30a,     &
                               rnn191,rnn192
   real, dimension( 31 ) ::    rn12a, rn12b, rn13, rn25a

!  common /rsnw2h/
    real ::         hn9,  hn10,  hn10a, hn14,  hn15a, hn16,  &
                                hn17, hn17a, hn19,  hn19a, hn20,  hn20b  
    real ::         gn17,gn17a,gn17a2                                     !4ice revised

!  common /rainmap/ 
    real ::         draimax

!  common /snomap/ 
  real ::  xs,sno11,sno00,dsno11,dsno00,sexp11,sexp00,stt,   &
                stexp,sbase,tslopes,dsnomin,dsnomin4,slim
!  common /grpmap/
  real ::  xg,grp11,grp00,dgrp11,dgrp00,gexp11,gexp00,gtt,   &
                gtexp,gbase,tslopeg,dgrpmin,dgrpmin4,glim
!  common /haimap/
  real :: hai00,hai11,htt0,htt1,haixp 

!  common /b3cg_2/ 
  real ::         ag2, bg2, bgh2, bgq2, roqg2, qrog2

!  common /rsnw2/ 
    real ::         rn142, rn152, rn15a2, rn17a2, rn192_2

!  common /icemass/
    real ::          ami50, ami40, ami100

!  common /bergcon/
   real, dimension( 31 ) ::    bergcon1,  bergcon2,       &
                                           bergcon3,  bergcon4
  real    :: cmin
  real    :: cpi
!
   real, dimension( 31 )  ::      aa1,  aa2
   data aa1/.7939e-7, .7841e-6, .3369e-5, .4336e-5, .5285e-5,         &
           .3728e-5, .1852e-5, .2991e-6, .4248e-6, .7434e-6,          &
           .1812e-5, .4394e-5, .9145e-5, .1725e-4, .3348e-4,          &
           .1725e-4, .9175e-5, .4412e-5, .2252e-5, .9115e-6,          &
           .4876e-6, .3473e-6, .4758e-6, .6306e-6, .8573e-6,          &
           .7868e-6, .7192e-6, .6513e-6, .5956e-6, .5333e-6,          &
           .4834e-6/
   data aa2/.4006, .4831, .5320, .5307, .5319,                        &
           .5249, .4888, .3894, .4047, .4318,                         &
           .4771, .5183, .5463, .5651, .5813,                         &
           .5655, .5478, .5203, .4906, .4447,                         &
           .4126, .3960, .4149, .4320, .4506,                         &
           .4483, .4460, .4433, .4413, .4382,                         &
           .4361/

!+---+-----------------------------------------------------------------+
!..the following 6 variables moved here to facilitate reflectivity
!.. calculation similar to other mp schemes, because when they get
!.. declared later in the code (now commented out), it makes things
!.. more difficult to integreate with the radar code.
      real ::     xnor, xnos, xnoh, xnog
      real ::     rhohail, rhograul
!.. values will be defined in subroutine fall_flux --- jjs 20150731
!      real    , parameter ::     xnor = 8.0e6
!      real    , parameter ::     xnos = 1.6e7
!      real    , parameter ::     xnoh = 2.0e5
!      real    , parameter ::     xnog = 4.0e6
!      real    , parameter ::     rhohail = 917.
!      real    , parameter ::     rhograul = 400.
!+---+-----------------------------------------------------------------+

 real,parameter :: rhowater_mks = 1.0e+3   !  water density      (kg/m^3)
 real,parameter :: rhosnow_mks  = 0.917e+3 !  snow density       (kg/m^3)

contains

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP


  subroutine goddard_microphysics_wrf(   th                               &
                       ,qv, qc                                      &
                       ,qr, qi                                      &
                       ,qs, qh                                      & ! 4ice
                       ,qg                                          &
                       ,rho, pii, p, dt_in, z                       &
                       ,ht, dz8w, w                           &
                       ,itimestep, xland, dx                        &
                       ,ids,ide, jds,jde, kds,kde                   & ! domain dims
                       ,ims,ime, jms,jme, kms,kme                   & ! memory dims
                       ,its,ite, jts,jte, kts,kte                   & ! tile   dims
                       ,rainnc, rainncv                             &
                       ,snownc, snowncv, sr                         &
                       ,graupelnc, graupelncv                       &
                       ,hailnc, hailncv                             & !hail
                       ,refl_10cm, diagflag, do_radar_ref           &
                       ,physc, physe, physd, physs, physm, physf    &
                       ,acphysc, acphyse, acphysd, acphyss, acphysm, acphysf &
                       ,refc, refr, refi, refs, refg, refh          & ! cloud effective radius
                       ,preci3d, precs3d, precg3d, prech3d, precr3d &
#if ( WRF_CHEM == 1)
                       ,aero, icn_diag, nc_diag, gid               &
                       ,chem_opt                                   &
                       ,gsfcgce_gocart_coupling                    &
#endif
                                                                   )
 implicit none
!------------------------------------------------------------------------------
! Comments: 
!  This WRF interface that drive Goddard Microphysics Package. 
!  All 3D variable (i,k,j) is chuncked to 2D vector (chunk,k).
!  3D i-k-j driver. 
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; Initial   
!        
! References: 
!
!------------------------------------------------------------------------------
 logical, optional, intent(in) :: &
       diagflag

 integer, intent(in)    ::   &
       ids,ide, jds,jde, kds,kde , &
       ims,ime, jms,jme, kms,kme , &
      its,ite, jts,jte, kts,kte 

 integer, optional, intent(in) :: do_radar_ref

 integer, intent(in)    ::   itimestep
  
 real, dimension( ims:ime , kms:kme , jms:jme ), intent(in) :: &
          rho, &
          pii, &
            p, &
         dz8w, &
            z, &
            w

 real, dimension( ims:ime , kms:kme , jms:jme ), intent(inout) :: &
          th, &
          qv, &
          qc, &  !qc? cloud?
          qr, &
          qi, &
          qs, &
          qg, &
          qh

 real, dimension( ims:ime , kms:kme , jms:jme ),                 &
       intent(inout) ::                            &
       refl_10cm, &  !
       refc,      &
       refr,      &
       refi,      &
       refs,      &
       refg,      &
       refh

 real, dimension( ims:ime , kms:kme , jms:jme ), intent(inout) ::  &
        physc, &
        physe, &
        physd, &
        physs, &
        physm, &
        physf, &
      acphysc, &
      acphyse, &
      acphysd, &
      acphyss, &
      acphysm, &
      acphysf, &
      preci3d, &
      precs3d, &
      precg3d, &
      prech3d, &
      precr3d

! 2D output
 real, dimension( ims:ime , jms:jme ), intent(inout) ::  &
        rainnc,    &
        rainncv,   &
        snownc,    &   
        snowncv,   &
        sr,        &
        graupelnc, &
        graupelncv,&
           hailnc, &
          hailncv

 real, dimension( ims:ime , jms:jme ), intent(in)  :: &
       xland, &   !land-ocean mask
          ht      !terrain height

 real, intent(in   ) :: &
        dt_in, & ! model time step [sec]
           dx    ! horizontal grid spacing [m]

#if ( WRF_CHEM == 1)
!
! for inline gocart coupling
!
 integer, parameter :: num_go = 14  ! number of the gocart aerosol species
 real, dimension( ims:ime, kms:kme, jms:jme, num_go), intent(in) :: aero
 real, dimension( ims:ime, kms:kme, jms:jme), intent(out) :: icn_diag, nc_diag
 integer,      intent(in   )    ::   gid
 integer,intent(in) :: chem_opt ! chem option
 integer,intent(in) :: gsfcgce_gocart_coupling ! coupling unit
#endif

!
! local variables
!
 integer :: kmax
 integer :: istep
 integer :: improve

 integer :: i, j, k, ip, ii, ic
 integer :: iskip, ih, icount, i24h 
 real    :: hour
 real    :: dth, dqrest, dqall, dqall1, rhotot, a1, a2 
 
 real, dimension(chunk, kms:kme):: th2d, qv2d, qc2d, qr2d
 real, dimension(chunk, kms:kme):: qi2d, qs2d, qg2d, qh2d
 real, dimension(chunk, kms:kme):: rho2d, pii2d, p2d, w2d
 real, dimension(chunk, kms:kme):: refc2d, refr2d, refi2d
 real, dimension(chunk, kms:kme):: refs2d, refg2d, refh2d
 real, dimension(chunk, kms:kme):: physc2d, physe2d, physd2d
 real, dimension(chunk, kms:kme):: physs2d, physm2d, physf2d
 real, dimension(chunk, kms:kme):: acphysc2d, acphyse2d, acphysd2d
 real, dimension(chunk, kms:kme):: acphyss2d, acphysm2d, acphysf2d
 real, dimension(chunk) :: xland1d
 real, dimension(chunk, kms:kme):: refl_10cm2d
#if ( WRF_CHEM == 1)
 real, dimension(chunk, kms:kme, num_go):: aero3d
 real, dimension(chunk, kms:kme):: icn_diag2d, nc_diag2d
#endif

#if ( WRF_CHEM == 1)
 logical,save :: first_makelut = .true. ! initial value, changed later.
#endif

!wrf radar reflectivity initialization only need to run once
 integer:: ncall=0

!
! gce options: use steve's new improvement   9/18/2009
!
! improve = 8
! new_ice_sat = 0, 1, 2, or 3 , ~ 9 
! new_ice_sat = 9
! istatmin = 180
! id = 0  without in-line staticstics
! id = 1  with in-line staticstics
! ibud = 0 no calculation of dth, dqv, dqrest and dqall
! ibud = 1 yes

#if (WRF_CHEM == 1)
!
! create ccn and in candidate look-up table. only need to run once at the first timestep of each run
!    if (itimestep .eq. 1) call makelut_ccn_icn
!
 if (first_makelut) then
   ! only execute when gocart and coupling turned on.
   if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
         chem_opt == 302 .or. chem_opt == 303) .and. &
        (gsfcgce_gocart_coupling == 1)) then
       if ( wrf_dm_on_monitor() ) print *,'goddard mp coupled with gocart aerosol'
           call makelut_ccn_icn
           first_makelut = .false. ! this value is saved for subsequent calls.
       end if
   end if
#endif

!
! set up constants used internally in goddard microphysics (should we call every time step?)
!
! call consat_s 


!
! calculte fallflux and precipiation in mks system (now in 1D driver) input are mks units in general. 
!
 kmax = kte-kts+1  !toshii change later? to (kme-kms+1)

 do j = jms, jme
    do i = ims , ime 

 call fall_flux( kmax, dt_in, vgc, vgc2, vhc, bhq, &
                 qc(i,:,j), qr(i,:,j), qi(i,:,j), qs(i,:,j), qg(i,:,j), qh(i,:,j), &
                 p(i,:,j), rho(i,:,j), th(i,:,j), pii(i,:,j), z(i,:,j), dz8w(i,:,j), ht(i,j), &
                 preci3d(i,:,j), precs3d(i,:,j), precg3d(i,:,j), prech3d(i,:,j), precr3d(i,:,j), &
                 rainnc(i,j), rainncv(i,j),              &
                 snownc(i,j), snowncv(i,j), sr(i,j),     &
                 graupelnc(i,j), graupelncv(i,j),        &
                 hailnc(i,j), hailncv(i,j)          )    
              

     enddo
 enddo

!
! emk nuwrf...moved this wrf radar reflectivity initialization to after
! fall_flux, as the rhohail and rhograul variables are set in that
! subroutine.
!

 if (ncall .eq. 0) then
!..set these variables needed for computing radar reflectivity.  these
!.. get used within radar_init to create other variables used in the
!.. radar module.
    xam_r = 3.14159*rhowater_mks/6.
    xbm_r = 3.
    xmu_r = 0.
    xam_s = 3.14159*rhosnow_mks/6.
    xbm_s = 3.
    xmu_s = 0.
    xam_g = 3.14159*rhograul/6.
    xbm_g = 3.
    xmu_g = 0.


    call radar_init
    ncall = 1
endif
!+---+-----------------------------------------------------------------+


!c negative values correction
!   iskip = 1
! 
!   if (iskip.eq.0) then
!      call negcor(qv,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,1,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(qc,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,2,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(qr,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,3,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(qi,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,4,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(qs,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,5,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(qg,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,6,             &
!                           its,ite,jts,jte,kts,kte)
!!   else if (mod(itimestep,i24h).eq.1) then
!!      print *,'no neg correction in mp at timestep=',itimestep
!   endif ! iskip

!c microphysics in gce

!$omp parallel do &
!$omp private ( ic, j, ii, i, k ) &
!$omp private ( th2d, qv2d, qc2d, qr2d ) &
!$omp private ( qi2d, qs2d, qg2d, qh2d ) &
!$omp private ( rho2d, pii2d, p2d, w2d ) &
!$omp private ( refc2d, refr2d, refi2d, refs2d, refg2d, refh2d ) &
!$omp private ( physc2d, physe2d, physd2d, physs2d, physm2d, physf2d ) &
!$omp private ( acphysc2d, acphyse2d, acphysd2d, acphyss2d, acphysm2d, acphysf2d ) &
!$omp private ( xland1d, refl_10cm2d ) &
#if ( WRF_CHEM == 1)
!$omp private ( aero3d, icn_diag2d, nc_diag2d) &
#endif
!$omp schedule(dynamic,1)


!
! vectorize 3D variables into 2D variables based on size 'chunk'.  
! (i-dim contains '(1+(ite-its+1)/chunk)' blocks of size 'chunk')
!
 istep = 1

 vector_loop: do ip = 1,((1+(ite-its+1)/chunk)*chunk)*(jte-jts+1),chunk 

    j  = jts+(ip-1)/((1+(ite-its+1)/chunk)*chunk)

    vector_check: if ( j .ge. jts .and. j .le. jte ) then ! j: [jts, jte]

      ii = its+mod((ip-1),((1+(ite-its+1)/chunk)*chunk)) ! ii: [its, ((1+(ite-its+1)/chunk)*chunk)]

      do ic=1,min(chunk,ite-ii+1)
         i = ii+ic -1
         xland1d(ic) = xland(i,j)
      enddo

      do k = kts, kte
         do ic=1,min(chunk,ite-ii+1) 
            i = ii+ic -1
            th2d(ic,k) = th(i,k,j)  ![K]
            qv2d(ic,k) = qv(i,k,j)  ![kg/kg] or [g/g]
            qc2d(ic,k) = qc(i,k,j)
            qr2d(ic,k) = qr(i,k,j)
            qi2d(ic,k) = qi(i,k,j)
            qs2d(ic,k) = qs(i,k,j)
            qg2d(ic,k) = qg(i,k,j)
            qh2d(ic,k) = qh(i,k,j)
            pii2d(ic,k) = pii(i,k,j)
            rho2d(ic,k) = rho(i,k,j) * 0.001 !mks [kg/m3] --> cgs [g/cm3]
            p2d(ic,k) = p(i,k,j) * 10.       !mks [Pa] -> cgs [Ba]
            w2d(ic,k) = w(i,k,j) * 100.      !mks [m] -> cgs [cm]
            refl_10cm2d(ic,k)=refl_10cm(i,k,j)
            refc2d(ic,k)=refc(i,k,j)  ![micron]
            refr2d(ic,k)=refr(i,k,j)
            refi2d(ic,k)=refi(i,k,j)
            refs2d(ic,k)=refs(i,k,j)
            refg2d(ic,k)=refg(i,k,j)
            refh2d(ic,k)=refh(i,k,j)
            physc2d(ic,k)=physc(i,k,j)  ![K/sec]
            physe2d(ic,k)=physe(i,k,j)
            physd2d(ic,k)=physd(i,k,j)
            physs2d(ic,k)=physs(i,k,j)
            physm2d(ic,k)=physm(i,k,j)
            physf2d(ic,k)=physf(i,k,j)
            acphysc2d(ic,k)=acphysc(i,k,j) ![K]
            acphyse2d(ic,k)=acphyse(i,k,j)
            acphysd2d(ic,k)=acphysd(i,k,j)
            acphyss2d(ic,k)=acphyss(i,k,j)
            acphysm2d(ic,k)=acphysm(i,k,j)
            acphysf2d(ic,k)=acphysf(i,k,j)
#if ( WRF_CHEM == 1)
            aero3d(ic,k,:)=aero(i,k,j,:)
#endif
         enddo
      enddo

      min_chunk: if ( min(chunk,ite-ii+1) .gt. 0 ) then

        !
        ! driver of goddard microphysics (inputs are generally cgs units)
        !
         call saticel_s( kmax, istep, dt_in, dx,     &
                         itimestep, xland1d,                           &
                        th2d, qv2d, qc2d, qr2d,                       &
                        qi2d, qs2d, qg2d, qh2d,                       &
                        rho2d, pii2d, p2d, w2d,                       & 
                        refl_10cm2d, diagflag, do_radar_ref,         & ! gt added for reflectivity calcs
                        refc2d, refr2d, refi2d, refs2d, refg2d, refh2d,  & ! cloud effective radius
                        physc2d, physe2d, physd2d, physs2d, physm2d, physf2d,  &
                        acphysc2d, acphyse2d, acphysd2d, acphyss2d, acphysm2d, acphysf2d,  &
#if ( WRF_CHEM == 1)
                        aero3d, icn_diag2d, nc_diag2d, gid,          &
                        chem_opt,                                     &
                        gsfcgce_gocart_coupling,                       &
#endif   
                        ii, j, min(chunk,ite-ii+1))       !vector chunk index

      endif min_chunk


    !
    ! Convert 2D vector back to 3D variables
    !
    do k = kts, kte
       do ic=1,min(chunk,ite-ii+1)
            i = ii+ic -1

            th(i,k,j) = th2d(ic,k)
            qv(i,k,j) = qv2d(ic,k)
            qc(i,k,j) = qc2d(ic,k)
            qr(i,k,j) = qr2d(ic,k)
            qi(i,k,j) = qi2d(ic,k)
            qs(i,k,j) = qs2d(ic,k)
            qg(i,k,j) = qg2d(ic,k)
            qh(i,k,j) = qh2d(ic,k)
            refc(i,k,j)=refc2d(ic,k)
            refr(i,k,j)=refr2d(ic,k)
            refi(i,k,j)=refi2d(ic,k)
            refs(i,k,j)=refs2d(ic,k)
            refg(i,k,j)=refg2d(ic,k)
            refh(i,k,j)=refh2d(ic,k)
            refl_10cm(i,k,j)=refl_10cm2d(ic,k)
            physc(i,k,j)=physc2d(ic,k)
            physe(i,k,j)=physe2d(ic,k)
            physd(i,k,j)=physd2d(ic,k)
            physs(i,k,j)=physs2d(ic,k)
            physm(i,k,j)=physm2d(ic,k)
            physf(i,k,j)=physf2d(ic,k)
            acphysc(i,k,j)=acphysc2d(ic,k)
            acphyse(i,k,j)=acphyse2d(ic,k)
            acphysd(i,k,j)=acphysd2d(ic,k)
            acphyss(i,k,j)=acphyss2d(ic,k)
            acphysm(i,k,j)=acphysm2d(ic,k)
            acphysf(i,k,j)=acphysf2d(ic,k)
#if ( WRF_CHEM == 1)
            icn_diag(i,k,j)=icn_diag2d(ic,k)
            nc_diag(i,k,j)=nc_diag2d(ic,k)
#endif
       enddo
    enddo

    endif vector_check

 enddo vector_loop

 return
 end subroutine goddard_microphysics_wrf

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine fall_flux ( kmax, dt, vgc, vgc2, vhc, bhq, &
                        qc1d, qr1d, qi1d, qs1d, qg1d, qh1d,                 &
                        p1d,rho1d, th1d, pi_mks1d, z1d, dz1d, topo,      &
                        preci1d, precs1d, precg1d, prech1d, precr1d,       &
                        rainnc, rainncv,              &
                        snownc, snowncv, sr,          &
                        graupelnc, graupelncv,        &
                        hailnc, hailncv ) 
                  
 implicit none
!------------------------------------------------------------------------------
! Comments: 
!  this program calculates precipitation fluxes due to terminal velocities.
!  This subroutine calculate fall flux of hydrometeor  
!  adopted from jiun-dar chern's codes of purdue regional model
!  k=1 is surface, k=kmax is top of atmosphere.
!  Input are mostly mks units. 
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; convert to 1d version & clean up + modernization. 
! adopted by jainn j. shi, 6/10/2005
! modified by goddard 7/24/2010
! modified by tao 11/12/2010
!        
! References: 
!
!------------------------------------------------------------------------------
 integer, intent(in)  :: kmax  !maximum vertical layers

 real,    intent(in)  :: &
       dt,  & ! model time step [sec]
       vgc, & ! 
       vgc2,& !
       vhc, & !
       bhq, & !
       topo   ! terrain height [m]
 
 real, dimension( kmax ), intent(inout) ::   & 
       qc1d, & ! cloud mixing ratio [kg/kg] or [g/g]
       qr1d, & ! rain mixing ratio [kg/kg] or [g/g]
       qi1d, & ! cloud ice mixing ratio [kg/kg] or [g/g]
       qs1d, & ! snow mixing ratio [kg/kg] or [g/g]
       qg1d, & ! graupel mixing ratio [kg/kg] or [g/g]
       qh1d    ! hail mixing ratio [kg/kg] or [g/g]
  
 real, dimension( kmax ), intent(in) ::   &  !1d
       th1d, &   ! potential temperature [K]
       rho1d,&   ! air density [kg/m3]
       z1d,  &   ! height [m]
       dz1d, &   ! layer thickness [m]
       p1d,  &   ! pressure [Pa]
       pi_mks1d  ! exner function [-]  

 real, intent(inout) ::      &  
       rainnc, rainncv,      &  ! accumulated [mm] and one-time step [mm/step] surface precipitation
       snownc, snowncv,      &  ! accumulated [mm] and one-time step [mm/step] surface snow []
       sr,                   &  ! snow / totao precip ratio [-]
       graupelnc, graupelncv,&  ! accumulated [mm] and one-time step [mm/step] surface graupel flux []
       hailnc, hailncv          ! accumulated [mm] and one-time step [mm/step] surface hail flux []

 real, dimension( kmax ), intent(out) ::                 & !1d
       preci1d, &   ! ice fall flux profile [kg/m2/s]
       precs1d, &   ! snow fall flux profile [kg/m2/s]
       precg1d, &   ! graupel fall flux profile [kg/m2/s]
       prech1d, &   ! hail fall flux profile [kg/m2/s]
       precr1d      ! rain fall flux profile [kg/m2/s]

!
! local variables vars
! 
 real, dimension( kmax )           :: fv
 real :: tmp1, term0
 real :: pptrain, pptsnow,        &
         pptgraul, pptice, ppthail
 real :: qrz, qiz
 real, dimension( kmax )       :: qcz, qsz, qgz, qhz,  &
                                  zz, dzw, prez, rhoz,      &
                                  orhoz,rho_cgs
 real, dimension( kmax )       :: csed, rsed, ised, ssed, gsed, hsed
 real, dimension( kmax )       :: thz, piz
 integer  :: k, i, j
 real, dimension( kmax )    :: vtr, vts, vtg, vth, vti
 real :: dtb, pi, consta, constc, gambp4,    &
         gamdp4, gam4pt5, gam4bbar
 real :: y1, y2, vr, vs, vg
 real :: vgcr, vgcr2, vscf,   &
         vhcr, vhcr2
 real :: tair, tairc, fexp
 real :: ftns, ftnsq, ftng, ftngq
 real :: const_vt, const_d, const_m, bb1, bb2
 real :: ftns0, ftng0
 real :: fros, fros0
 real :: qhz2,qgz2


! will be defined later using consat values 
 real ::     rhowater 
 real ::     rhosnow 
 real, parameter ::                              &
       constb = 0.8, constd = 0.11, o6 = 1./6.,            &
       cdrag = 0.6
 real, parameter ::     abar = 19.3, bbar = 0.37,           &
                            p0 = 1.0e5
 real, parameter ::     rhoe_s = 1.29
! for terminal velocity flux
 integer  :: min_q, max_q
 real     :: t_del_tv, del_tv, flux, fluxin, fluxout
 logical  :: notlast


!
! initialize constant
!
 dtb=dt
 pi=acos(-1.)

 xnor = tnw*1.0e8
 rhowater = roqr*1000.
 xnos = tns*1.0e8             ! consistent with consat
 rhosnow = roqs*1000.

!
! tng and roqg are assigned with hail numbers in consat if (ihail.eq.1)
!
 xnog = tng*1.0e8
 rhograul = roqg*1000.   ! mks kg/m3
 xnoh = tnh*1.0e8      
 rhohail = roqh*1000.   ! mks kg/m3

 consta=2115.0*0.01**(1-constb)
 constc=78.63*0.01**(1-constd)

!
! gamma function
!
   gambp4=gammagce(constb+4.)
   gamdp4=gammagce(constd+4.)
   gam4pt5=gammagce(4.5)
   gam4bbar=gammagce(4.+bbar)

!
! initialization (zero)
!
 preci1d(:)=0.
 precs1d(:)=0.
 precg1d(:)=0.
 prech1d(:)=0.
 precr1d(:)=0.
 ised(:)=0.
 ssed(:)=0.
 gsed(:)=0.
 hsed(:)=0.
 rsed(:)=0.

 pptrain = 0.
 pptsnow = 0.
 pptgraul = 0.
 ppthail = 0.
 pptice  = 0.

 do k = 1, kmax 
    rho_cgs(k)=rho1d(k)*0.001         !rho in cgs
    fv(k)=sqrt(rhoe_s/rho1d(k))
 enddo !k


!
!-- rain
!
 t_del_tv=0.
 del_tv=dtb
 notlast=.true.

 rain_loop: do while (notlast)
!
    min_q= kmax 
    max_q= 0 
!

    do k = 1, kmax-1 

       vtr(k)=0.
       qrz=qr1d(k)
       if (qrz .gt. cmin) then
          min_q=min0(min_q,k)
          max_q=max0(max_q,k)

          ! new codes from steve's in cgs
          tair=th1d(k)*pi_mks1d(k)
          tairc=tair-t0
          y1=qrz  
          y2=qc1d(k)  
          call vqrqi(1,rho_cgs(k),fv(k),y1,y2,tair,vtr(k))     
          vtr(k)=vtr(k) * 0.01 ! terminal velocity convert back to mks   

          if (.not. vtr(k) .gt. 0.0) cycle 

          if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(z1d(k)-topo)/vtr(k))
          else
               del_tv=amin1(del_tv,0.9*(z1d(k)-z1d(k-1))/vtr(k))

          endif
       endif
    enddo !do k

    if (max_q .ge. min_q) then
!
!- check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

       t_del_tv=t_del_tv+del_tv
!
       if ( t_del_tv .ge. dtb ) then
            notlast=.false.
            del_tv=dtb+del_tv-t_del_tv
       endif

! use small delta t to calculate the qrz flux
! termi is the qrz flux pass in the grid box through the upper boundary
! termo is the qrz flux pass out the grid box through the lower boundary
!
       fluxin=0.
       do k=max_q,min_q,-1
          qrz=qr1d(k)
          fluxout=rho1d(k)*vtr(k)*qrz
          flux=(fluxin-fluxout)/rho1d(k)/dz1d(k)
          qrz=qrz+del_tv*flux
          qrz=amax1(0.,qrz)
          qr1d(k)=qrz
          fluxin=fluxout
          rsed(k)=rsed(k)+fluxin
         enddo
         if (min_q .eq. 1) then
            pptrain=pptrain+fluxin*del_tv
         else
            qrz=qr1d(min_q-1)
            qrz=amax1(0., qrz+del_tv*  &
                          fluxin/rho1d(min_q-1)/dz1d(min_q-1) )
            qr1d(min_q-1)=qrz
         endif
!
    else
         notlast=.false.
    endif

 enddo rain_loop

!
!-- snow
!
 t_del_tv=0.
 del_tv=dtb
 notlast=.true.

 snow_loop: do while (notlast)
!
    min_q=kmax 
    max_q=0    

!
    do k = 1, kmax-1 
       vts(k)=0.

       if (qs1d(k) .gt. cmin) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)

            ! new codes from steve's in cgs
            y1 = qs1d(k)             
            vscf=vsc*fv(k)

            ftns=1.
            ftns0=1.
            tair=th1d(k)*pi_mks1d(k)
            tairc=tair-t0  !air temperature [degC]

            qhz2=qh1d(k)  !toshii why need hail here?
            qgz2=qg1d(k)  !toshii why need graupel here?
            if (k .lt. kmax-2 .and. tairc .ge. -5) then
              qhz2=qh1d(k+1) 
              qgz2=qg1d(k+1)
            endif
            call sgmap(1,qs1d(k),qg1d(k),qgz2,qh1d(k),qhz2,rho_cgs(k),tairc,ftns0)  !snow intercept
            ftns=ftns0**bsq
            fros=1.
            call sgmap(3,qs1d(k),qg1d(k),qgz2,qh1d(k),qhz2,rho_cgs(k),tairc,fros0)  !snow density
            fros=fros0**bsq
            vts(k)=max(vscf*(rho_cgs(k)*y1)**bsq/ftns/fros,0.e0)
            ! bs, bsq, vscf are defined in new consat_s
            vts(k)=vts(k) * 0.01  ! convert back to mks  

            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(z1d(k)-topo)/vts(k))
            else
               del_tv=amin1(del_tv,0.9*(z1d(k)-z1d(k-1))/vts(k))
            endif !k
      endif !cmin
   enddo  ! do k

   if (max_q .ge. min_q) then

!
!- check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

      t_del_tv=t_del_tv+del_tv

      if ( t_del_tv .ge. dtb ) then
           notlast=.false.
           del_tv=dtb+del_tv-t_del_tv
      endif

! use small delta t to calculate the qs1d flux
! termi is the qs1d flux pass in the grid box through the upper boundary
! termo is the qs1d flux pass out the grid box through the lower boundary
!
       fluxin=0.
       do k=max_q,min_q,-1
            fluxout=rho1d(k)*vts(k)*qs1d(k)
            flux=(fluxin-fluxout)/rho1d(k)/dz1d(k)
            qs1d(k)=qs1d(k)+del_tv*flux
            qs1d(k)=amax1(0.,qs1d(k))
            fluxin=fluxout
            ssed(k)=ssed(k)+fluxin
       enddo

       if (min_q .eq. 1) then
            pptsnow=pptsnow+fluxin*del_tv
       else
            qs1d(min_q-1)=amax1(0., qs1d(min_q-1)+del_tv*  &
                         fluxin/rho1d(min_q-1)/dz1d(min_q-1) )
       endif
!
    else
         notlast=.false.
    endif

 enddo snow_loop

!
!--- graupel
!
 t_del_tv=0.
 del_tv=dtb
 notlast=.true.

 graupel_loop: do while (notlast)

    min_q= kmax 
    max_q= 0    
!
    do k= 1, kmax-1 
       vtg(k)=0.

       if (qg1d(k) .gt. cmin) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)

! for graupel, based on rh (1984)
! new codes from steve's in cgs
            y1 = qg1d(k)             
            vgcr=vgc*fv(k)
            vgcr2=vgc2*fv(k)
            ftng=1.
            ftng0=1.
            tair=th1d(k)*pi_mks1d(k)
            tairc=tair-t0
            qhz2=qh1d(k)
            qgz2=qg1d(k)
            if (k .lt. kmax-2 .and. tairc .ge. -5) then
                qhz2=qh1d(k+1)
                 qgz2=qg1d(k+1)
            endif
            !ftng0 is graupel intercept
            call sgmap(2,qs1d(k),qg1d(k),qgz2,qh1d(k),qhz2,rho_cgs(k),tairc,ftng0) 
            ftng=ftng0**bgq
            vtg(k)=amax1(vgcr*(rho_cgs(k)*y1)**bgq/ftng, 0.0)
                                       ! bg, vgcr, bgq are defined in new consat_s
            vtg(k)=vtg(k) * 0.01  ! convert back to mks
            if (y1.gt.qrog2)then                              !di
                ftng=ftng0**bgq2                                                !di
                vtg(k)=amax1(vgcr2*(rho_cgs(k)*y1)**bgq2/ftng, 0.0)                    !di
                                  ! bg, vgcr, bgq are defined in new consat_s  !di
                vtg(k)=vtg(k) * 0.01  ! convert back to mks
            endif 

            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(z1d(k)-topo)/vtg(k))
            else
               del_tv=amin1(del_tv,0.9*(z1d(k)-z1d(k-1))/vtg(k))
            endif 
!
      endif !qg1d
   enddo !k

   if (max_q .ge. min_q) then

!
!- check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

       t_del_tv=t_del_tv+del_tv

       if ( t_del_tv .ge. dtb ) then
            notlast=.false.
            del_tv=dtb+del_tv-t_del_tv
       endif

!
! use small delta t to calculate the qg1d flux
! termi is the qg1d flux pass in the grid box through the upper boundary
! termo is the qg1d flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rho1d(k)*vtg(k)*qg1d(k)
            flux=(fluxin-fluxout)/rho1d(k)/dz1d(k)
            qg1d(k)=qg1d(k)+del_tv*flux
            qg1d(k)=amax1(0.,qg1d(k))
            fluxin=fluxout
            gsed(k)=gsed(k)+fluxin
         enddo
         if (min_q .eq. 1) then
            pptgraul=pptgraul+fluxin*del_tv
         else
            qg1d(min_q-1)=amax1(0., qg1d(min_q-1)+del_tv*  &
                         fluxin/rho1d(min_q-1)/dz1d(min_q-1) )
         endif

   else
      notlast=.false.
   endif

 enddo graupel_loop

!
! hail
!
 t_del_tv=0.
 del_tv=dtb
 notlast=.true.

 hail_loop: do while (notlast)

    min_q= kmax 
    max_q= 0    

    do k = 1, kmax-1 
       vth(k)=0.

       if (qh1d(k) .gt. cmin) then
           min_q=min0(min_q,k)
           max_q=max0(max_q,k)

           ! new codes from steve's in cgs
           y1 = qh1d(k)
           vhcr=vhc/sqrt(rho_cgs(k))                       
           vth(k)=amax1(vhcr*(y1*rho_cgs(k))**bhq, 0.e0)             !di
           vth(k)=vth(k) * 0.01  ! convert back to mks

           if (k .eq. 1) then
              del_tv=amin1(del_tv,0.9*(z1d(k)-topo)/vth(k))
           else
              del_tv=amin1(del_tv,0.9*(z1d(k)-z1d(k-1))/vth(k))
           endif

       endif !qh1d
   enddo !k

   if (max_q .ge. min_q) then

!
!- check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)
       t_del_tv=t_del_tv+del_tv

       if ( t_del_tv .ge. dtb ) then
            notlast=.false.
            del_tv=dtb+del_tv-t_del_tv
       endif

! use small delta t to calculate the qh1d flux
! termi is the qh1d flux pass in the grid box through the upper boundary
! termo is the qh1d flux pass out the grid box through the lower boundary
!
       fluxin=0.
       do k=max_q,min_q,-1
          fluxout=rho1d(k)*vth(k)*qh1d(k)
          flux=(fluxin-fluxout)/rho1d(k)/dz1d(k)
          qh1d(k)=qh1d(k)+del_tv*flux
          qh1d(k)=amax1(0.,qh1d(k))
          fluxin=fluxout
          hsed(k)=hsed(k)+fluxin
       enddo
       if (min_q .eq. 1) then
           ppthail=ppthail+fluxin*del_tv
       else
          qh1d(min_q-1)= amax1(0., qh1d(min_q-1)+del_tv*  &
                       fluxin/rho1d(min_q-1)/dz1d(min_q-1) )
       endif

    else
       notlast=.false.
    endif
!
 enddo hail_loop

!
! cloud ice  (03/21/02) follow vaughan t.j. phillips at gfdl
!
 t_del_tv=0.
 del_tv=dtb
 notlast=.true.
!
 ice_loop: do while (notlast)

    min_q= kmax 
    max_q= 0    

    do k = 1, kmax-1 
       qiz=qi1d(k)
       vti(k)=0.
       if (qiz .gt. cmin) then
           min_q=min0(min_q,k)
           max_q=max0(max_q,k)

           ! new codes from steve's in cgs.
           vti(k)=0.
           y1=rho1d(k) * 1000. * qiz    ! y1 in g/m**3
           if (y1 .ge. 1.e-6) then
              y1=qiz
              y2=qc1d(k)
              tair=th1d(k)*pi_mks1d(k)
              tairc=tair-t0
              call vqrqi(2,rho_cgs(k),fv(k),y1,y2,tair,vti(k))      ! di
           endif  !y1
           vti(k)=vti(k) * 0.01                       ! convert back to mks

          !  prevent divsion by zero and underflow value
           if (vti(k) .gt. 1.e-20) then
             if (k .eq. 1) then
                del_tv=amin1(del_tv,0.9*(z1d(k)-topo)/vti(k))
             else
                del_tv=amin1(del_tv,0.9*(z1d(k)-z1d(k-1))/vti(k))
             endif
           endif
        endif
   enddo

   if (max_q .ge. min_q) then

!
!- check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

       t_del_tv=t_del_tv+del_tv

       if ( t_del_tv .ge. dtb ) then
            notlast=.false.
            del_tv=dtb+del_tv-t_del_tv
       endif

! use small delta t to calculate the qiz flux
! termi is the qiz flux pass in the grid box through the upper boundary
! termo is the qiz flux pass out the grid box through the lower boundary
!

       fluxin=0.
       do k=max_q,min_q,-1
          qiz=qi1d(k)
          fluxout=rho1d(k)*vti(k)*qiz
          flux=(fluxin-fluxout)/rho1d(k)/dz1d(k)
          qiz=qiz+del_tv*flux
          qiz=amax1(0.,qiz)
          qi1d(k)=qiz
          fluxin=fluxout
          ised(k)=ised(k)+fluxin
       enddo
       if (min_q .eq. 1) then
           pptice=pptice+fluxin*del_tv
       else
           qiz=qi1d(min_q-1)
           qiz=amax1(0., qiz+del_tv*  &
               fluxin/rho1d(min_q-1)/dz1d(min_q-1) )
           qi1d(min_q-1)=qiz
       endif

    else
       notlast=.false.
    endif

 enddo ice_loop


!
! precipitation profile for each hydrometeor
!
 do k = 1, kmax 
    preci1d(k)=ised(k)
    precs1d(k)=ssed(k)
    precg1d(k)=gsed(k)
    prech1d(k)=hsed(k)
    precr1d(k)=rsed(k)
 enddo

!
! surface precipitation 
!
 snowncv = pptsnow
 snownc = snownc + pptsnow
 graupelncv = pptgraul
 graupelnc = graupelnc + pptgraul 
 hailncv = ppthail
 hailnc = hailnc + ppthail
 rainncv = pptrain + pptsnow + pptgraul + pptice + ppthail
 rainnc  = rainnc + pptrain + pptsnow + pptgraul + pptice + ppthail
 sr = 0.
 if (rainncv .gt. 0.) sr = (pptsnow + pptgraul + pptice + ppthail) / rainncv

 
 return
 end subroutine fall_flux

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine negcor ( x, rho, dz8w,                         &
                    ims,ime, jms,jme, kms,kme,              & ! memory dims
                    itimestep, ics,                         &
                    its,ite, jts,jte, kts,kte               ) ! tile   dims

!------------------------------------------------------------------------------
! Comments: 
!  This subroutine correct negative values of hydrometeor.    
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; Initial   
! adopted by jainn j. shi, 6/10/2005
! modified by goddard 7/24/2010
! modified by tao 11/12/2010
!        
! References: 
!
!------------------------------------------------------------------------------

  integer,intent(in) :: ims,ime, jms,jme, kms,kme, &
                        its,ite, jts,jte, kts,kte

  real, dimension( ims:ime , kms:kme , jms:jme ),                 &
        intent(inout) ::                                     x   
  real, dimension( ims:ime , kms:kme , jms:jme ),                 &
        intent(in   ) ::                              rho, dz8w  
  integer, intent(in   ) ::                           itimestep, ics 

!c local variables
!  real, dimension( kts:kte ) ::  y1, y2
  real   ::   a0, a1, a2
  integer :: i,j,k

  a1=0.
  a2=0.
  do k=kts,kte
     do j=jts,jte
        do i=its,ite
        a1=a1+max(x(i,k,j), 0.)*rho(i,k,j)*dz8w(i,k,j)
        a2=a2+max(-x(i,k,j), 0.)*rho(i,k,j)*dz8w(i,k,j)
        enddo
     enddo
  enddo

!  a1=0.0
!  a2=0.0
!  do k=kts,kte
!     a1=a1+y1(k)
!     a2=a2+y2(k)
!  enddo

  a0=0.0

  if (a1.ne.0.0.and.a1.gt.a2) then 
     a0=(a1-a2)/a1

  if (mod(itimestep,540).eq.0) then
     if (ics.eq.1) then
        write(61,*) 'kms=',kms,'  kme=',kme,'  kts=',kts,'  kte=',kte
        write(61,*) 'jms=',jms,'  jme=',jme,'  jts=',jts,'  jte=',jte 
        write(61,*) 'ims=',ims,'  ime=',ime,'  its=',its,'  ite=',ite 
     endif 
     if (ics.eq.1) then
         write(61,*) 'qv timestep=',itimestep
         write(61,*) '  a1=',a1,'   a2=',a2,'   a0=',a0
     else if (ics.eq.2) then
             write(61,*) 'qc timestep=',itimestep
             write(61,*) '  a1=',a1,'   a2=',a2,'   a0=',a0
     else if (ics.eq.3) then
             write(61,*) 'qr timestep=',itimestep
             write(61,*) '  a1=',a1,'   a2=',a2,'   a0=',a0
     else if (ics.eq.4) then
             write(61,*) 'qi timestep=',itimestep
             write(61,*) '  a1=',a1,'   a2=',a2,'   a0=',a0
     else if (ics.eq.5) then
             write(61,*) 'qs timestep=',itimestep
             write(61,*) '  a1=',a1,'   a2=',a2,'   a0=',a0
     else if (ics.eq.6) then
             write(61,*) 'qg timestep=',itimestep
             write(61,*) '  a1=',a1,'   a2=',a2,'   a0=',a0
     else
             write(61,*) 'wrong cloud specieis number'
     endif 
  endif 

     do k=kts,kte
        do j=jts,jte
           do i=its,ite
           x(i,k,j)=a0*amax1(x(i,k,j), 0.0)
           enddo
        enddo
     enddo
  endif

  end subroutine negcor

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine consat_s 
 implicit none
!------------------------------------------------------------------------------
! Comments: 
!  This subroutine define constant used for goddard microphysics.     
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; modernization / clean up. 
!
!   implemented into wrf  by roger shi 2006/2007                       
!   additional modifications by tao, roger and steve 2009              
!   july 25 2010                                                       
!   tao november 12 2010                                               
!   coded by tao (1989-2003), modified by s. lang (2006/07)            
!        
! References: 
!
!   tao, w.-k., and j. simpson, 1989: modeling study of a tropical     
!   squall-type convective line. j. atmos. sci., 46, 177-202.          
!                                                                      
!   tao, w.-k., j. simpson and m. mccumber, 1989: an ice-water        
!   saturation adjustment. mon. wea. rev., 117, 231-235.               
!                                                                      
!   tao, w.-k., and j. simpson, 1993: the goddard cumulus ensemble     
!   model. part i: model description. terrestrial, atmospheric and     
!   oceanic sciences, 4, 35-72.                                        
!                                                                      
!   tao, w.-k., j. simpson, d. baker, s. braun, m.-d. chou, b.         
!   ferrier,d. johnson, a. khain, s. lang,  b. lynn, c.-l. shie,       
!   d. starr, c.-h. sui, y. wang and p. wetzel, 2003: microphysics,    
!   radiation and surface processes in the goddard cumulus ensemble    
!   (gce) model, a special issue on non-hydrostatic mesoscale          
!   modeling, meteorology and atmospheric physics, 82, 97-137.         
!                                                                      
!   lang, s., w.-k. tao, r. cifelli, w. olson, j. halverson, s.        
!   rutledge, and j. simpson, 2007: improving simulations of           
!   convective system from trmm lba: easterly and westerly regimes.    
!   j. atmos. sci., 64, 1141-1164.                                     
!                                                                      
!------------------------------------------------------------------------------
 real    :: cn0

 integer k
 real :: ga3, ga4, ga5, ga7, ga8, ga9, ga3g2, ga4g2, ga5g2, ga6d      
 real :: ga3h, ga4h, ga5hh, bc1, dc1, esc, egs, erc, amc, ehs, ehg
 real :: ehw, ehi, ehr, sc13, ga6, ga5gh2, egc, cpi2, grvt, tca, dwv
 real :: dva, amw, ars, rw, cw, ci, cd1, cd2, ga3b, ga4b, ga6b, ga5bh
 real :: ga3g, ga4g, ga5gh, ga3d, ga4d, ga5dh, ac1, ac2, ac3, cc1, eri
 real :: ami, esr, eiw, ui50, ri50, cmn, y1, egi, egr, apri, bpri

 logical,parameter :: old_coeff = .false.


 print*, 'MSG consat_s: initialize microphysics constant'

!
!
!
  cmin=1.e-20
      al = 2.5e10
      cp = 1.004e7
      rd1 = 1.e-3
      rd2 = 2.2
!jjs
      cpi=4.*atan(1.)
      cpi2=cpi*cpi
      grvt=980.               !hail                 
!
      c38=3.799052e3
      c358=35.86
      c610=6.1078e3
      c149=1.496286e-5
      c879=8.794142
      c172=17.26939
      c409=4098.026
      c76=7.66
      c218=21.87456
      c580=5807.695
      c141=1.414435e7

!
      tca=2.43e3
      dwv=.226
      dva=1.718e-4
      amw=18.016
      ars=8.314e7
!      scv=2.2904487
!
      t0=273.16
      t00=238.16
      alv=2.5e10    !latent heat vaporization
      alf=3.336e9
      als=2.8336e10
      avc=alv/cp
      afc=alf/cp
      asc=als/cp
      rw=4.615e6
      cw=4.187e7
      ci=2.093e7

!
! define the coefficients used in terminal velocity
!  define the density and size distribution of precipitation
!
! tnw                  !  rain intercept (1/cm**4)
! roqr                 !  rain density (g/cm**3)
! tns                  !  snow intercept (1/cm**4)
! roqs                 !  snow density (g/cm**3)
! tng                  !  graupel intercept (1/cm**4)
! roqg                 !  graupel density (g/cm**3)
!
! hail parameters    
!      tnh = 0.002                 !hail  medium size
       tnh = 0.01
!      tnh = 0.0002                !hail large size
!      tnh = 0.02                  !hail small size
      roqh = 0.9                  !hail
      cd1=6.e-1                   !hail
      cd2=4.*grvt/(3.*cd1)        !hail
      ah=sqrt(cd2*roqh)           !hail
      bh=0.5                      !hail


!!!!!!!!check /08/08/12!!!!!!!!!!!!
         roqg=.3      !graupel bulk density
!         ag=341.84    !for bulk density of 0.4
         ag=330.22    !for bulk density of 0.3
!         ag=314.54    !for bulk density of 0.2         
         bg=.36
         ag2=544.83   !for bulk density of 0.5
         bg2=.54
         tng=.04

       roqg2=0.5
       qrog2=2.0            !g/m**3
!      qrog2=9.0            !g/m**3    !fixed size test
       qrog2=qrog2*1.e-6
       ag=330.22
       ag2=544.83
       bg=0.36
       bg2=0.54

!
! snow parameters 
!                             6/15/02
!      tns=1.
!      tns=.08             ! if ice913=1, tao's
      tns=.16              ! if ice913=0, tao's snow intercept
      roqs=.1              ! snow density g/cm
      as=78.63154          ! coefficient a in snow terminal velocity
      bs=.11               ! coefficient b in snow terminal velocity
        roqs=.05        !  snow density (g/cm**3)
        tns=0.1         !  snow intercept (1/cm**4)
        as=151.01
        bs=0.24
!
! rain parameters       
!
      aw=2115.
      bw=.8
      roqr=1.  !not defined in steve's
      tnw=.08  !not defined in steve's
      bgh=.5*bg
      bgh2=.5*bg2
      bsh=.5*bs
      bwh=.5*bw
      bhh=.5*bh
      bgq=.25*bg
      bgq2=.25*bg2
      bsq=.25*bs
      bwq=.25*bw
      bhq=.25*bh          !hail

!
!  define the snow and graupel size mapping parameters for improve=3
!
       tslopes=0.11177209   ! increase tns by  50 from 0 to -35c !sgmap000
       tslopeg=0.05756866   ! increase tng by 7.5 from 0 to -35c

       draimax=0.0500 !maximum rain diameter (cm)
       draimax=draimax**4.*roqr*cpi
       dsnomin=0.0100 !minimum snow diameter (cm) !fin16
       dsnomin4=dsnomin**4.*0.900*cpi  !density at 100 microns

       dgrpmin=0.0135 !minimum graupel diameter (cm)
       dgrpmin4=dgrpmin**4.*roqg*cpi

if(old_coeff) then
!orignal old 4ICE
       xs=0.97
       sno11=0.70   !fin15
       sno00=-0.24   !fin17
       dsno11=3.25   !fin8
       dsno00=0.40   !fin12
       sexp11=1.5    !sgmapfin17
       sexp00=0.9    !fin15
       stt=-35.      !iterate
       stexp=0.42    !fin18
       slim=1.00     !fin5
       sbase=0.04000 !sgmap00s

       xg=0.98
       grp11=0.30    !sgmap00x
       grp00=0.30    !sgmap00x
       dgrp11=2.70   !sgmap001
       dgrp00=2.60
       gexp11=0.20   !sgmap004
       gexp00=0.20   !sgmap004
       gtt=-35.      !sgmap006 made consistent
       gtexp=0.40
       glim=0.90
       gbase=0.0130  !sgmap00u

       hai00=3.25     !g/m**3        fin20
       hai11=0.50     !g/m**3        fin20
       htt0=-5.00     !deg C         fin20
       htt1=-50.00    !deg C         fin20
       haixp=2.7      !exponent      fin20

else

!
! sgmap snow
!
 xs=0.97
 sno11=0.70    !fin15
 sno00=-0.20   !post1
 dsno11=3.35   !post1
 dsno00=0.60   !post1
 sexp11=1.6    !post1
 sexp00=1.0    !post1
 stt=-35.      !iterate
 stexp=0.36    !post1
 slim=1.00     !fin5
 sbase=0.04000 !sgmap00s

!
! sgmap graupel
!
 xg=0.98
 grp11=0.30    !sgmap00x
 grp00=0.30    !sgmap00x
 dgrp11=2.70   !sgmap001
 dgrp00=2.60
 gexp11=0.20   !sgmap004
 gexp00=0.20   !sgmap004
 gtt=-35.      !sgmap006 made consistent
 gtexp=0.40
 glim=0.90
 gbase=0.0130  !sgmap00u

!
! sgmap hail
!
 hai00=3.50     !g/m**3        fin20
 hai11=0.50     !g/m**3        fin20
 htt0=-5.00     !deg c         fin20
 htt1=-50.00    !deg c         fin20
 haixp=2.5      !exponent      fin21

endif

!
! gamma function calculations
!
      ga3=2.
      ga4=6.
      ga5=24.
      ga6=120.
      ga7=720.
      ga8=5040.
      ga9=40320.
!
      ga3b  = gammagce(3.+bw)
      ga4b  = gammagce(4.+bw)
      ga6b  = gammagce(6.+bw)
      ga5bh = gammagce((5.+bw)/2.)
      ga3g  = gammagce(3.+bg)
      ga4g  = gammagce(4.+bg)
      ga5gh = gammagce((5.+bg)/2.)
      ga3d  = gammagce(3.+bs)
      ga4d  = gammagce(4.+bs)
      ga5dh = gammagce((5.+bs)/2.)

        ga4g=11.63177
        ga3g=3.3233625          
        ga5gh=1.608355         
        if(bg.eq.0.37) ga4g=9.730877 
        if(bg.eq.0.37) ga3g=2.887512
        if(bg.eq.0.37) ga5gh=1.526425
        if(bg.eq.0.36) ga4g=9.599978
        if(bg.eq.0.36) ga3g=2.857136
        if(bg.eq.0.36) ga5gh=1.520402 
        if(bg2.eq.0.54) ga4g2=12.298653
        if(bg2.eq.0.54) ga3g2=3.474196
        if(bg2.eq.0.54) ga5gh2=1.635061 
          ga3d=2.54925           
          ga4d=8.285063         
          ga5dh=1.456943
          if(bs.eq.0.57) ga3d=3.59304
          if(bs.eq.0.57) ga4d=12.82715
          if(bs.eq.0.57) ga5dh=1.655588
          if(bs.eq.0.24) ga3d=2.523508
          if(bs.eq.0.24) ga4d=8.176166
          if(bs.eq.0.24) ga5dh=1.451396
          if(bs.eq.0.11) ga3d=2.218906
          if(bs.eq.0.11) ga4d=6.900796
          if(bs.eq.0.11) ga5dh=1.382792 

      ga6d=144.93124
        if(bs.eq.0.24) ga6d=181.654791
      ga3h=gammagce(3.+bh)                     !hail
      ga4h=gammagce(4.+bh)                     !hail
      ga5hh=gammagce((5.+bh)/2.)               !hail
!
! lin et al., 1983 or lord et al., 1984 
!
      ac1=aw
      ac2=ag            ! steve only defines ac1 and ac2
      ac3=as            ! need to talk about these 3 parameters.

      bc1=bw
      cc1=as
      dc1=bs

      zrc=(cpi*roqr*tnw)**0.25
      zsc=(cpi*roqs*tns)**0.25
      zgc=(cpi*roqg*tng)**0.25
      zgc2=(cpi*roqg2*tng)**0.25
      zhc=(cpi*roqh*tnh)**0.25               !hail

      vrc=aw*ga4b/(6.*zrc**bw)

       vrc0=-26.7
       vrc1=20600./zrc
       vrc2=-204500./(zrc*zrc)
       vrc3=906000./(zrc*zrc*zrc)

      vsc=as*ga4d/(6.*zsc**bs)
      vgc=ag*ga4g/(6.*zgc**bg)
      vgc2=ag2*ga4g2/(6.*zgc2**bg2)
      vhc=ah*ga4h/(6.*zhc**bh)          !hail

      rn1=9.4e-15
      rn2=1.e-3
      bnd2=2.0e-3                    ! if ice913=0 6/15/02 tao's

      esi=0.70
      rn3=.25*cpi*tns*as*esi*ga3d

      esc=0.45
      rn4=.25*cpi*esc*tns*as*ga3d

      eri=.1                        ! 6/17/02 tao's ice913=0 (not 1)
      rn5=.25*cpi*eri*tnw
      rn50=-.267e2*ga3
      rn51=5.15e3*ga4
      rn52=-1.0225e4*ga5
      rn53=7.55e3*ga6

      ami=1./(24.*6.e-9)            ! 6/15/02 tao's
      rn6=cpi2*eri*tnw*roqr*ami
      rn60=-.267e2*ga6
      rn61=5.15e3*ga7
      rn62=-1.0225e4*ga8
      rn63=7.55e3*ga9

      esr=1.                       ! also if ice913=1 for tao's
      rn7=cpi2*esr*tnw*tns*roqs
      rn8=cpi2*esr*tnw*tns*roqr

      egs=.1
      rn9=cpi2*egs*tns*tng*roqs

      rn10=4.*tns
      rn101=.65
      rn102=.44*sqrt(as/dva)*ga5dh
      rn10a=alv*als*amw/(tca*ars)
      rn10b=alv/tca
      rn10c=ars/(dwv*amw)

      rn11=2.*cpi*tns*tca/alf
      rn11a=cw/alf

      ami50=4.8e-7*(dsnomin*1.e4/2./50.)**3
      ami100=1.51e-7    ! improve < 3
      ami40=2.46e-7*.875**3     !fin24    35 microns

       eiw=1.
       ui50=100. ! 6/15/02 tao's
       ri50=dsnomin/2.

      cmn=1.05e-15
      rn12=cpi*eiw*ui50*ri50*ri50
      do k=1,31
         y1=1.-aa2(k)
         rn13(k)=aa1(k)*y1/(ami50**y1-ami40**y1)
         rn12a(k)=rn13(k)/ami50
         rn12b(k)=aa1(k)*ami50**aa2(k)
         rn25a(k)=aa1(k)*cmn**aa2(k)

         bergcon1(k)=6.*aa1(k)*ami50**(aa2(k)-1.)
         bergcon2(k)=-2.*aa1(k)*ami50**aa2(k)*1.2

         bergcon3(k)=6.*aa2(k)/((aa2(k)+1.)*(aa2(k)+2.))        &
                   *aa1(k)*ami50**(aa2(k)-1.)
         bergcon4(k)=2.*(1.-aa2(k))/((aa2(k)+1.)*(aa2(k)+2.))   &
                   *aa1(k)*ami50**aa2(k)*1.2
      enddo
!
      egc=0.65
      rn14=.25*cpi*egc*ag*tng*ga3g
      rn142=.25*cpi*egc*ag2*tng*ga3g2     !high dens graupel

      egi=.1
      rn15=.25*cpi*egi*tng*ag*ga3g
      rn15a=.25*cpi*egi*tng*ag*ga3g
      rn15a2=.25*cpi*egi*tng*ag2*ga3g2    !high dens graupel
      rn152=.25*cpi*egi*tng*ag2*ga3g2     !high dens graupel

      egr=1.
      rn16=cpi2*egr*tng*tnw*roqr
      rn17=2.*cpi*tng
      rn17a2=.31*ga5gh2*sqrt(ag2/dva)     !high dens graupel
      rn17b=cw-ci
      rn17c=cw
      rn171=2.*cpi*tng*alv*dwv
      rn172=2.*cpi*tng*tca
      rn17a=.31*ga5gh*sqrt(ag/dva)
!
      apri=.66
      bpri=1.e-4
      bpri=0.5*bpri                        ! 6/17/02 tao's
      rn18=20.*cpi2*bpri*tnw*roqr
      rn18a=apri
      rn191=.78                     !are this same as rnn191 (listed below)?
      rn192=.31*ga5gh*sqrt(ag/dva)  !are this same as rnn192 (listed below)?
      rn192_2=.31*ga5gh2*sqrt(ag2/dva)    !high dens graupel
      rn19b=cw/alf 
      rn19=2.*cpi*tng*tca/alf
      rn19a=cw/alf
!
      rn20=2.*cpi*tng
      rn20a=als*als*amw/(tca*ars)
      rn20b=als/tca
      rn30a=alv*alv*amw/(tca*ars)
      rn30=2.*cpi*tng

      bnd3=2.e-3
      rn21=1.e-3
      bnd21=1.5e-3
!
      erc=1.
      rn22=.25*cpi*erc*tnw
!
      rn23=2.*cpi*tnw
      rn23a=.31*ga5bh*sqrt(ac1)
      rn23b=alv*alv/rw
      rn231=.78
      rn232=.31*ga3*sqrt(3.e3/dva)
!
!
!        "c0" in routine      "consat" (2d), "consatrh" (3d)
!
         cn0=1.e-8
         beta=-.6
!
      rn25=cn0
      rn30b=alv/tca
      rn30c=ars/(dwv*amw)
      rn31=1.e-17

      rn32=4.*51.545e-4
!!
      rn33=4.*tns
      rn331=.65
      rn332=.44*sqrt(as/dva)*ga5dh 
!
      amc=1./(24.*4.e-9)
      rn34=cpi2*esc*amc*as*roqs*tns*ga6d  
      rn35=alv*alv/(cp*rw)

      gn17=2.*cpi*tng                                                      !4ice revised
      gn17a=.31*ga5gh*sqrt(ag)                                             !4ice revised
      gn17a2=.31*ga5gh2*sqrt(ag2)                                          !4ice revised

      ehs=1.
      hn9=cpi2*ehs*tns*tnh*roqs
      ehg=0.3
      hn10=cpi2*ehg*tng*tnh*roqg
      ehw=1.
      hn14=.25*cpi*ehw*tnh*ga3h*ah
      ehi=1.
      hn15a=.25*cpi*ehi*tnh*ga3h*ah
      ehr=1.
      hn16=cpi2*ehr*tnh*tnw*roqr
      hn17=2.*cpi*tnh
      hn17a=.31*ga5hh*sqrt(ah)
      hn19=2.*cpi*tnh/alf
      hn19a=.31*ga5hh*sqrt(ah)
      hn10a=als*als/rw
      hn20=2.*cpi*tnh
      hn20b=.31*ga5hh*sqrt(ah)

!
! jdc add schmidt number term
!
      sc13    = 0.8420526
      rn102   = rn102 * sc13
      rn17a   = rn17a * sc13
      rn17a2  = rn17a2 * sc13
      rn192   = rn192 * sc13
      rn192_2 = rn192_2 * sc13
      rn232   = rn232 * sc13
      rn332   = rn332 * sc13


  return
 end subroutine consat_s 

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

  real function gammagce (xx)
  implicit none
!------------------------------------------------------------------------------
! Comments: 
!  Gamma function used in GCE.     
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; modernization / clean up. 
!        
! References: 
!                                                                      
!------------------------------------------------------------------------------
  real*8 cof(6),stp,half,one,fpf,x,tmp,ser
  data cof,stp /  76.18009173,-86.50532033,24.01409822, &
     -1.231739516,.120858003e-2,-.536382e-5, 2.50662827465 /
  data half,one,fpf / .5, 1., 5.5 /
!
  real xx
  real gammln
  integer j

      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do  j=1,6
         x=x+one
        ser=ser+cof(j)/x
      enddo !j
      gammln=tmp+log(stp*ser)
!jjs
      gammagce=exp(gammln)
!jjs
      return
 end function gammagce

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine sgmap(isg,qcs,qcg,qcg2,qch,qch2,r00,tairc,ftnsg)
 implicit none
!------------------------------------------------------------------------------
! Comments: 
!  Parameterize base snow/graupel intercept scaling factor - ftnsg
!  Snow requires graupel and hail for break up estimation. 
!
! isg -- flag 
!        1 ftnsg=snow intercept scaling factor
!        2 ftnsg=graupel intercept scaling factor
!        3 ftnsg=snow density scaling factor
! tns is the actual intercept
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; modernization / clean up. 
!
! References: 
!                                                                      
!------------------------------------------------------------------------------
 integer, intent(in)  :: isg
 real,    intent(in)  :: r00, tairc
 real,    intent(in)  :: qcs,qcg,qcg2,qch,qch2
 real,    intent(out) :: ftnsg

! local variables
 real :: taird, qsg, qsg1, xx, fexp  !, cpi, cmin Di
 real :: ftnsT, sno1, dsno1, sexp1, ftnsQ, tnsmax, densno
 real :: ftngT, grp1, dgrp1, gexp1, ftngQ, tngmax
 real :: hx, gx, hgx
 real :: kk
 real :: hai2, dhai1


 ftnsg=1.

 if (isg == 1 .or. isg == 3) qsg=qcs
 if (isg == 2) qsg = qcg
 if (isg == 4) qsg = qch

 if (qsg > cmin) then

   qsg1=qsg*r00*1.e6

   !
   ! snow mapping
   !
   if(isg == 1.or.isg == 3)then  !snow 4ice
      taird=min(0.,max(stt,tairc)+0.0)
      ftnsT=exp(-1.*tslopes*taird)
      sno1=sno11
      dsno1=dsno11
      sexp1=sexp11

      if (taird > stt) then
         sno1=sno00-(sno00-sno11)*(taird/stt)**stexp
         dsno1=dsno00-(dsno00-dsno11)*(taird/stt)**stexp

         sexp1=sexp00-(sexp00-sexp11)*(taird/stt)**stexp
      endif !taird

      hx=0.
      gx=0.
      hgx=1.

      if(qch2*1.e6*r00 > 0.008)    &
         hx=max(1.,qch2*1.e6*r00*125.)
      if(qcg2*1.e6*r00 > 0.040)    &
         gx=max(1.,qcg2*1.e6*r00*25.)
      hgx=hx+gx
      hgx=max(1.,hgx)

      xx=xs-xs*min(slim,max(0.,(qsg1-sno1)/dsno1)**sexp1)
      ftnsT=ftnsT**xx
      fexp=xx
      ftnsQ=1.0
      ftnsQ=(qsg1/sbase)**fexp
      ftnsg=ftnsT*ftnsQ*hgx
      tnsmax=r00*qsg/dsnomin4

      if (ftnsg*tns > tnsmax) ftnsg=tnsmax/tns               !4ice
      densno=0.001996*(ftnsg*tns/qsg/r00)**0.2995   !Brandes et al. 2007
      if(isg == 3) ftnsg=min(densno,0.9)/roqs            !4ice

    !
    ! graupel mapping
    !
    else if (isg == 2) then                      !graupel
      taird=min(0.,max(gtt,tairc)+0.0)
      ftngT=exp(-1.*tslopeg*taird)
      grp1=grp11
      dgrp1=dgrp11
      gexp1=gexp11

      if (taird > gtt) then
          grp1=grp00-(grp00-grp11)*(taird/gtt)**gtexp
          dgrp1=dgrp00-(dgrp00-dgrp11)*(taird/gtt)**gtexp
          gexp1=gexp00-(gexp00-gexp11)*(taird/gtt)**gtexp
      endif !taird

      xx=xg-xg*min(glim,max(0.0,(qsg1-grp1)/dgrp1)**gexp1)
      ftngT=ftngT**xx
      fexp=xx
      ftngQ=1.0
      ftngQ=(qsg1/gbase)**fexp
      ftnsg=ftngT*ftngQ
      tngmax=r00*qsg/dgrpmin4
      if (ftnsg*tng > tngmax) ftnsg=tngmax/tng

    !
    ! hail mapping
    !
    elseif(isg == 4)then  !hail
      hai2=hai00
      if(tairc <= htt0 .and. tairc >= htt1) &  !fin18
         hai2=hai11+(hai00-hai11)*((tairc-htt1)/(htt0-htt1))**haixp   !fin18
      if(tairc <= htt1)  hai2=hai11     !fin18
         dhai1=hai2          !x2
      if(qsg1 >= hai2)    &
         ftnsg=1.0-0.80*min(max((qsg1-hai2)/dhai1,0.),1.) !  fin16

     endif !isg

  endif !qsg

 return
 end subroutine sgmap

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine vqrqi(isg,rho_cgs,fv,qri,qc,tair,ww1)
 implicit none
!------------------------------------------------------------------------------
! Comments: 
!  compute fall speed of cloud rain and ice.
!  Scalar subroutine. 
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; modernization / clean up. (physically identical to vqrqi.F_v3.04 )
!
! References: 
!                                                                      
!------------------------------------------------------------------------------
 integer, intent(in) :: isg  ! 1 for rain 2 for ice
 real, intent(in)    :: rho_cgs, fv,qri,qc,tair
 real, intent(inout) :: ww1 

! local variables
 integer :: ic 
 real  :: y1,vr,vs,vg
 real  :: const_vt, const_d, const_m !cpi, cmin di
 real  :: bb1, bb2,ice_fall
 real  :: bin_factor, ftnw, ftnwmin
 real, dimension(7) ::  aice, vice
 data aice/1.e-6, 1.e-5, 1.e-4, 1.e-3, 0.01, 0.1, 1./
 data vice/5,15,30,35,40,45,50/
 real :: ww

 ice_fall=0.
 const_vt=1.49e4
 const_d=11.9
 const_m=1./5.38e7

 y1=rho_cgs*qri
 ww1=0.

 if (y1 > cmin) then

   if (isg == 1) then   !  rain

     ftnw=1.                                                       
     if(qc < cmin .and. tair > t0)then                      
        bin_factor = min( 0.11*(1000.*qri)**-1.27 + 0.98 , 1.30)
        ftnw=1./bin_factor**3.35                             
        ftnwmin=rho_cgs*qri/draimax                        
        if(qri < 0.001) ftnw=max(ftnw,ftnwmin/tnw)   
     endif                                               

     vs = sqrt( y1 )
     vg = sqrt( vs )
     vr = vrc0 + vrc1*vg/ftnw**0.25 + vrc2*vs/ftnw**0.50 + vrc3*vg*vs/ftnw**0.75
     ww1 =max(fv*vr, 0.e0)  ! [cm/s]

   else if (isg == 2) then  ! cloud ice

       y1=1.e6*rho_cgs*qri   ! to g/m**3

       if (y1 > 1.e-6) then
                  y1=y1*1.e-3
                  bb1=const_m*y1**0.25
                  bb2=const_d*bb1**0.5
                  ww=max(const_vt*bb2**1.31, 0.0)
                  ww1 = min( 50. , ww*100. ) ! [cm/s]
                  !if (ww1 .gt. 50.) ww1=50.   ! slang
       endif  !y1
   endif  !isg
 endif !y1

 return
 end subroutine vqrqi

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine saticel_s (kmax, istep, dt, dx,                   &
                       itimestep, xland,                               &
                       pt2d, qv2d, qc2d, qr2d,                           &
                       qi2d, qs2d, qg2d, qh2d,                         &
                       rho2d, pi2d, p2d, w2d,                          &
                       refl_10cm, diagflag, do_radar_ref,              & ! gt added for reflectivity calcs
                       refc2d, refr2d, refi2d, refs2d, refg2d, refh2d, & ! cloud effective radius
                       physc2d, physe2d, physd2d, physs2d, physm2d, physf2d,       &
                       acphysc2d, acphyse2d, acphysd2d, acphyss2d, acphysm2d, acphysf2d &
#if ( WRF_CHEM == 1)
                       ,aero, icn_diag, nc_diag, gid,                   &
                       chem_opt,                                       &
                       gsfcgce_gocart_coupling                        &
#endif
                       ,ii, j, irestrict  )  !vector parameter
 implicit none
!------------------------------------------------------------------------------
! Comments: 
!  compute ice phase microphysics and saturation processes
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; modernization / clean up. 
!
!   coded by tao (1989-2003), modified by s. lang (2006/07)               
!   implemented into wrf  by jainn shi 2006/2007                          
!   improved by s. lang (2008-2009)                                       
!   implemented by tao, jainn shi and tested by jainn shi                 
!   modified by tao, jul. 2010                                            
!   modified by tao, aug. 2010                                            
!   added 4ice code (from s. lang) into nu-wrf by di wu in 2014           
!   added aerosol coupling by jainn shi, mar. 2014                        
!   added cloud droplet eff. radius code by jainn shi, dec. 2014          
!   code optimization by j. mielikainen (2016)   
!
! References: 
!   tao, w.-k., and j. simpson, 1989: modeling study of a tropical        
!   squall-type convective line. j. atmos. sci., 46, 177-202.             
!                                                                         
!   tao, w.-k., j. simpson and m. mccumber, 1989: an ice-water            
!   saturation adjustment. mon. wea. rev., 117, 231-235.                  
!                                                                         
!                                                                         
!   tao, w.-k., and j. simpson, 1993: the goddard cumulus ensemble        
!   model. part i: model description. terrestrial, atmospheric and        
!   oceanic sciences, 4, 35-72.                                           
!                                                                         
!   tao, w.-k., j. simpson, d. baker, s. braun, m.-d. chou, b.            
!   ferrier,d. johnson, a. khain, s. lang,  b. lynn, c.-l. shie,          
!   d. starr, c.-h. sui, y. wang and p. wetzel, 2003: microphysics,       
!   radiation and surface processes in the goddard cumulus ensemble       
!   (gce) model, a special issue on non-hydrostatic mesoscale             
!   modeling, meteorology and atmospheric physics, 82, 97-137.            
!                                                                         
!   lang, s., w.-k. tao, r. cifelli, w. olson, j. halverson, s.           
!   rutledge, and j. simpson, 2007: improving simulations of              
!   convective system from trmm lba: easterly and westerly regimes.       
!   j. atmos. sci., 64, 1141-1164.                                        
!                                                                         
!   tao, w.-k., j. j. shi,  s. lang, c. peters-lidard, a. hou, s.         
!   braun, and j. simpson, 2007: new, improved bulk-microphysical         
!   schemes for studying precipitation processes in wrf. part i:          
!   comparisons with other schemes.                                       
!                                                                         
!   lang, s., w.-k. tao, x. zeng, and y. li, 2011: reducing the biases in 
!   simulated radar reflectivities from a bulk microphysics scheme:       
!   tropical convective systems. j. atmos. sci., 68, 2306-2320.           
!                                                                         
!   shi, j. j., t. matsui, w.-k. tao, c. peters-lidard, m. chin, q. tan,  
!   k. pickering, n. guy, s. lang, and e. kemp., 2014: implementation of  
!   an aerosol-cloud microphysics-radiation coupling into the nasa        
!   unified wrf:  simulation results for the 6-7 august 2006 amma special 
!   observing period. quart. j. roy. meteor. soc., 140, 2158-2175,         
!   doi:10.1002/qj.2286. 
!                                                                      
!------------------------------------------------------------------------------
 integer,intent(inout) :: istep  ! time step (from re-start)
 integer,intent(in) :: &
   kmax      , & ! maximum vertical layers
   itimestep , & ! total time step (from model initialization)
   ii        , & !
    j        , & !
   irestrict     ! global vector i-index inside local i-loops: ii+i-1

 real, intent(in) :: &
     dt ,& ! timestep [sec]
     dx    ! horizonltal grid spacing [m]

 real , dimension( chunk ) , intent(in) :: xland

 real, dimension (chunk, kmax), intent(inout) ::  &
      pt2d,& ! potential temperature [K]
      qv2d,& ! water vapor mixing ratio [kg/kg] or [g/g]
      qc2d,& ! cloud mixing ratio [kg/kg] or [g/g]
      qr2d,& ! rain mixing ratio [kg/kg] or [g/g]
      qi2d,& ! ice mixing ratio [kg/kg] or [g/g]
      qs2d,& ! snow mixing ratio [kg/kg] or [g/g]
      qg2d,& ! graupel mixing ratio [kg/kg] or [g/g]
      qh2d   ! hail mixing ratio [kg/kg] or [g/g]

 real, dimension (chunk, kmax), intent(in) :: &
      rho2d,   &  ! air density [g/cm3] cgs
      pi2d,    &  ! exner function [-]
      p2d,     &  ! pressure [Ba] or [dyn/cm2]
      w2d         ! vertical velocity [cm] 

!  variables for the calculation of effective radius of cloud species
 real, dimension (chunk, kmax) , intent(inout) ::  &
      refc2d, & ! cloud effective radius [um]
      refr2d, & ! rain effective radius [um]
      refi2d, & ! ice effective radius [um]
      refs2d, & ! snow effective radius [um]
      refg2d, & ! graupel effective radius [um]
      refh2d    ! hail effective radius [um]

 real, dimension(chunk, kmax), intent(inout):: &
      refl_10cm  ! radar reflectivity at 10cm wavelenth [dBZ]

  real, dimension (chunk, kmax), intent(inout)  ::  &
       physc2d, & ! diabatic heating rate due to condensation [K/s]
       physe2d, & ! diabatic heating rate due to evaporation [K/s]
       physd2d, & ! diabatic heating rate due to deposition [K/s]
       physs2d, & ! diabatic heating rate due to sublimatiopn [K/s]
       physm2d, & ! diabatic heating rate due to melting [K/s]
       physf2d, & ! diabatic heating rate due to freezing [K/s]
       acphysc2d, & ! accumulated diabatic heating rate due to condensation [K]
       acphyse2d, & ! accumulated diabatic heating rate due to evaporation [K]
       acphysd2d, & ! accumulated diabatic heating rate due to deposition [K]
       acphyss2d, & ! accumulated diabatic heating rate due to sublimatiopn [K]
       acphysm2d, & ! accumulated diabatic heating rate due to melting [K]
       acphysf2d    ! accumulated diabatic heating rate due to freezing [K]


!
! local variables
!

 integer i, k, kp     

 real, dimension(chunk) :: &
    afcp, alvr, ascp, avcp, rp0, pi0, pir,  &
    pr0, r00, rrs, rrq, fv0, fvs, cp409, &
    rr0, zrr, zsr, zgr, zhr, &
    cp580, cs580, cv409, vscf, vgcf, vgcf2, &
    vhcr, dwvp, r3f, r4f, r5f, r6f, &
   r12r, r14f ,r14f2, r15af, r15af2, r15f, r18r, &
    r22f, r25rt, r32rt, r331r, r332rf, &
                            r34f

 real :: bg3, bg3_2, bgh5, bgh5_2, bs3 ,bs6, bsh5, bw3 ,bw6 ,bwh5, &
         cmin, cmin1, cmin2, d2t, del, f2 ,f3, ft, qb0, r25a, r_nci, &
         sccc, sddd, seee, sfff, smmm, ssss, tb0, temp, ucog ,ucog2, &
         ucor ,ucos, ucoh, uwet, rdt, bnd1, c_nci, &
         r10t, r20t, r23t

 real :: a_1, a_2, a_3, a_4
 real :: a_11, a_22, a_33, a_44
 real :: zdry, zwet, zwet0
 real :: vap_frac
 real ::   thresh_evap !threshold of evaporation


 real, dimension (chunk, kmax) ::  fv

 real, dimension (chunk, kmax) ::  qsz, qgz,qhz

!jjs      common /badv/
  real, dimension (chunk) ::                                  &              
           vg,      zg,       &
           ps,      pg,       &
          prn,     psn,       &
        pwacs,   wgacr,       &
        pidep,    pint,       &
          qsi,     ssi,       &
          esi,     esw,       &
          qsw,      pr,       &
          ssw,   pihom,       &
         pidw,   pimlt,       &
        psaut,   qracs,       &
        psaci,   psacw,       &
        qsacw,   praci,       &
        pmlts,   pmltg,       &
        asss

!jjs      common/bsat/
  real, dimension (chunk) ::        &
        praut,   pracw,       &
         psfw,    psfi,       &
        dgacs,   dgacw,       &
        dgaci,   dgacr,       &
        pgacs,   wgacs,       &
        qgacw,   wgaci,       &
        qgacr,   pgwet,       &
        pgaut,   pracs,       &
        psacr,   qsacr,       &
         pgfr,   psmlt,       &
        pgmlt,   psdep,       &
        pgdep,   piacr,       &
          egs

!jjs      common/bsat1/
  real, dimension (chunk) ::        &
           pt,      qv,       &
           qc,      qr,       &
           qi,      qs,       &
           qg,      qh,       &
                  tair,       &
        tairc,   rtair,       &
          dep,      dd,       &
          dd1,     qvs,       &
           dm,      rq,       &
        rsub1,     col,       &
          cnd,     ern,       &
         dlt1,    dlt2,       &
         dlt3,    dlt4,       &
           zr,      vr,       &
           zs,      vs
       

!jjs      common/bsat2h/
  real, dimension (chunk) ::        &
         phfr,phmlt,                          & !4ice
         dhacw,qhacw,dhacr,qhacr,whacr,       & !4ice
         dhaci,whaci,dhacs,phacs,whacs,       & !4ice
         dhacg,whacg,phwet,phdep,phsub,       & !4ice
         pvaph,primh,scv,dwv,tca                !4ice

!jjs      common/b5/ 
!  real, dimension (chunk, kmax) ::  rho !only in satice in cgs

!jjs      common/b6/
  real, dimension (chunk, kmax) ::  p0, pi, f0, ww1
  real, dimension (chunk, kmax) ::    & 
           fd,      fe,           &
           st,      sv,           &
           sq,      sc,           &
           se,     sqa

!jjs      common/bi/ it(its:ite,jts:jte), ics(its:ite,jts:jte,4)
  integer, dimension (chunk) ::        it  
  integer, dimension (chunk, 4) ::    ics 

  integer :: i24h
  real :: r2is, r2ig, r2ih
  


! emk nuwrf
  real, dimension(chunk, kmax) :: dbz

!jjs  9/30/2009 for steve's new improvement

  integer  ::  ihalmos
  real     ::  xnsplnt, xmsplnt
  real     ::  hmtemp1, hmtemp2, hmtemp3, hmtemp4
  real     ::  ftnw, ftnwmin
  real     ::  xssi, fssi, rssi, xsubi, wssi
  real     ::  dmicrons, dmicrong, dvair, alpha
  real,   dimension (chunk) :: tairn, tairi,    &
                               ftns,  ftng,    &
                               ftns0, ftng0,    &
                               ftnh,  ftnh0,    &
                               pihms, pihmg,    &
                               pimm,  pcfr,    &
                               pssub, pgsub,    &
                               fros, fros0,    &
                               vi,    zh,    &
                               vh, pihmh,    &
                               dda0, pvapg,    &
                               pracg,qracg,    &
                               qrimh,pg2h      !4ice revised

  real,   dimension (chunk) :: y1, y2, y3, y4,  &
                               y5, y6, y7, y8 
! for xiping's new dbz code
  real     ::  hfact, sfact, yy1
  real     ::  xncld, esat, rv, rlapse_m
  real     ::  delt, bhi  !di deleted cpi
  real     ::  rc, ra, cna
  real     ::  xccld, xknud, cunnf, diffar
  real     ::  qgz2, qhz2

  real :: r11t, r19t, r19at, r30t, r33t
  real, dimension(chunk) :: r7rf, r8rf, r9rf, r16rf
  real, dimension(chunk) :: r101r, r102rf, r191r, r192rf, r192rf2
  real, dimension(chunk) :: r231r, r232rf
  real, dimension(chunk) :: h9r, h10r, h14r, h15ar, h16r, h17r, h17aq,   &
              h19aq, h19rt, h10ar, h20t, h20bq
  real, dimension(chunk) :: bin_factor, rim_frac
  real     :: term1, term2, fdwv, dwv0  !jdc water vapor diffusivity correction term
  integer  :: iter

#if ( WRF_CHEM == 1)
! for inline gocart coupling
  integer, parameter :: num_go = 14  ! number of the gocart aerosol species
  integer,      intent(in   )    ::   gid
  integer,intent(in) :: chem_opt ! emk
  integer,intent(in) :: gsfcgce_gocart_coupling ! emk
  real, dimension( chunk, kmax, num_go), intent(in) :: aero
  real, dimension( chunk, kmax ), intent(out) :: icn_diag !in concentration [#/litre]
  real, dimension( chunk, kmax ), intent(out) :: nc_diag !cloud concentration [#/cm3]

! local variables
  real :: e_sat, e_dry  !saturated and dry air water vapor [hpa, mb]
  real :: rh_rad     ! relative humidity [%]
  real :: super_sat  !super saturation [%]
  real :: ccn_out(chunk)  ! ccn conc [#/cm3] ! emk test
  real :: icn_out(chunk)  ! in conc [#/litter] ! emk test
  real :: p_liu_daum  ! autoconversion rate [g/cm3 s-1]    !
  real :: re_liu_daum ! effective radius of cloud [micron]   !
  real,parameter :: min_icn = 0.01 !minimum # conc of in [#/litre]
#endif


  real, parameter :: roqi = 0.9179    ! ice density
  real, parameter :: ccn_over_land = 1500  ! [#/cm3] climatological value
  real, parameter :: ccn_over_water = 150  ! [#/cm3] climatological value
  real :: l_cloud    ! cloud water [g/cm3] !
  real :: i_cloud    ! cloud water [g/cm3] !
  real :: mu, ccn_ref, lambda
  real :: gamfac1, gamfac3


  logical, optional, intent(in) :: diagflag
  integer, optional, intent(in) :: do_radar_ref


!
! program start 
!


! dwv0 is water vapor diffusivity at stp
 dwv0 =  0.226

!
! initialize output latent heat 
!
 if (itimestep.eq.1) then
   do k = 1, kmax !kts, kte
!dir$ vector aligned
      do i=1,irestrict
         physc2d(i,k)=0.
         physe2d(i,k)=0.
         physd2d(i,k)=0.
         physs2d(i,k)=0.
         physf2d(i,k)=0.
         physm2d(i,k)=0.
         acphysc2d(i,k)=0.
         acphyse2d(i,k)=0.
         acphysd2d(i,k)=0.
         acphyss2d(i,k)=0.
         acphysf2d(i,k)=0.
         acphysm2d(i,k)=0.
    enddo
 enddo !k

    !if ( wrf_dm_on_monitor() .and. itimestep.eq.1 .and. istep == 1) then
    !   write(6, *) ' msg saticel_s: latent heating array is zeroed at timestep = ', itimestep
    !endif

 endif


 do k= 1, kmax !kts,kte
!dir$ vector aligned
    do i=1,irestrict
       fv(i,k)=sqrt(rho2d(i,1)/rho2d(i,k))
    enddo
 enddo !k
!jjs

!
!  three classes of ice-phase   (lin et al, 1983)  
!
 d2t=dt

 r2ig=1.
 r2is=1.
 r2ih=1.

!c  tao 2007 end

 cmin=1.e-20
!jjs  10/7/2008
 cmin1=1.e-20
 cmin2=1.e-35

 xssi=0.16  ! maximum allowable ice supersaturation from sel
 rssi=0.05  ! maximum allowable residual ice supersaturation sel
 xsubi=0.70  ! ice rh below which qi must sublimate

!  hallet-mossop rime splintering parameters
 ihalmos=1
 xnsplnt=370.     ! peak # splinters per milligram of rime
 xmsplnt=4.4e-8   ! mass of a splinter (from ferrier 1994)
 hmtemp1=-2.
 hmtemp2=-4.
 hmtemp3=-6.
 hmtemp4=-8.

 it(:)=1
!!!!!!!!08/07/12
 f2=rd1*d2t     
 f3=rd2*d2t

 ft=dt/d2t

 bw3=bw+3.
 bs3=bs+3.
 bg3=bg+3.
 bh3=bh+3
 bsh5=2.5+bsh
 bgh5=2.5+bgh
 bhh5=2.5+bhh
 bwh5=2.5+bwh
 bw6=bw+6.
 bs6=bs+6.
! betah=.5*beta

 r10t=rn10*d2t
 r25a=rn25

 rdt=1./d2t
 r11t=rn11*d2t
 r19t=rn19*d2t
 r19at=rn19a*d2t
 r20t=rn20*d2t
 r23t=rn23*d2t
 r30t=rn30*d2t
 r33t=rn33*d2t
 bg3_2=bg2+3
 bgh5_2=2.5+bgh2

!jjs 20150831         
! calculate the threshold for reducing spurious evaporation
! a function of dx (grid spacing) 

 thresh_evap = -39.974 * exp(-1.194 * dx/1000.)


 if ( wrf_dm_on_monitor() .and. itimestep.eq.1 .and. istep == 1) then
!  print *,'msg gsfc 4ice saticel_s: improve=',improve
   print *,'dx, thresh_evap = ', dx, thresh_evap  
   print *,'no reduce suprious evaporation adjustment'
 endif

 rc=1.e-3               ! cloud droplet radius 10 microns
 ra=1.e-5               ! aerosol radius 0.1 microns
 cna=500.               ! contact nuclei conc per cc
 bhi=1.01e-2            ! pollen (deihl et al. 2006)

!
! ================ k-loop start =========================
!
 giant_kloop: do k=1, kmax !kts,kte
    kp=k+1
    tb0=0.
    qb0=0.

!
! vectorization for i loop for zero out all variables. 
!

!dir$ vector aligned
    do i=1,irestrict


! emk bug fix...initialize variable
#if ( WRF_CHEM == 1)
       ccn_out(i) = 0.e0
       icn_out(i) = 0.e0
       if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
             chem_opt == 302 .or. chem_opt == 303) .and. &
             (gsfcgce_gocart_coupling == 1) ) then

          icn_diag(i,k) = icn_out(i)  ! #/litre
          nc_diag(i,k) = ccn_out(i)  ! #/cm3

       else
          icn_diag(i,k) = 0.
          nc_diag(i,k) = 0.
       end if
#endif
! emk end

       rp0(i)=3.799052e3/p2d(i,k)
       pi0(i)=pi2d(i,k)
       pir(i)=1./(pi2d(i,k))
       pr0(i)=1./p2d(i,k)
       r00(i)=rho2d(i,k)
       rr0(i)=1./rho2d(i,k)
       rrs(i)=sqrt(rr0(i))
       rrq(i)=sqrt(rrs(i))
       f0(i,k)=al/cp/pi2d(i,k)
       fv0(i)=fv(i,k)
       fvs(i)=sqrt(fv(i,k))
       zrr(i)=1.e5*zrc*rrq(i)
       zsr(i)=1.e5*zsc*rrq(i)
       zgr(i)=1.e5*zgc*rrq(i)
       zhr(i)=1.e5*zhc*rrq(i)
       cp409(i)=c409*pi0(i)
       cv409(i)=c409*avc !
       cp580(i)=c580*pi0(i)
       cs580(i)=c580*asc !
       alvr(i)=r00(i)*alv
       afcp(i)=afc*pir(i)
       avcp(i)=avc*pir(i)
       ascp(i)=asc*pir(i)
       vscf(i)=vsc*fv0(i)
       vgcf(i)=vgc*fv0(i)
       vgcf2(i)=vgc2*fv0(i)
       vhcr(i)=vhc*rrs(i)
       dwvp(i)=c879*pr0(i)

       r3f(i)=rn3*fv0(i)
       r4f(i)=rn4*fv0(i)
       r5f(i)=rn5*fv0(i)
       r6f(i)=rn6*fv0(i)

       r12r(i)=rn12*r00(i)
       r14f(i)=rn14*fv0(i)
       r14f2(i)=rn142*fv0(i)
       r15f(i)=rn15*fv0(i)
       r15af(i)=rn15a*fv0(i) !4ice revised
       r15af2(i)=rn15a2*fv0(i) !4ice revised
       r18r(i)=rn18*rr0(i)
       r22f(i)=rn22*fv0(i)
       r25rt(i)=rn25*rr0(i)*d2t
       r32rt(i)=rn32*d2t*rrs(i)

       r7rf(i)=rn7*rr0(i)*fv0(i)
       r8rf(i)=rn8*rr0(i)*fv0(i)
       r9rf(i)=rn9*rr0(i)*fv0(i)
       r16rf(i)=rn16*rr0(i)*fv0(i)
       r101r(i)=rn101*rr0(i)
       r102rf(i)=rn102*rrs(i)*fvs(i)
       r191r(i)=rn191*rr0(i)
       r192rf(i)=rn192*rrs(i)*fvs(i)
       r192rf2(i)=rn192_2*rrs(i)*fvs(i)
       r331r(i)=rn331*rr0(i)
       r332rf(i)=rn332*rrs(i)*fvs(i)
       r34f(i)=rn34*fv0(i)

       r231r(i)=rn231*rr0(i)
       r232rf(i)=rn232*rrs(i)*fvs(i)

       h9r(i)=hn9*rr0(i)
       h10r(i)=hn10*rr0(i)
       h14r(i)=hn14*rrs(i)
       h15ar(i)=hn15a*rrs(i)
       h16r(i)=hn16*rr0(i)
       h17r(i)=hn17*rr0(i)
       h17aq(i)=hn17a*rrq(i)
       h19aq(i)=hn19a*rrq(i)
       h19rt(i)=hn19*rr0(i)*d2t
       h10ar(i)=hn10a*r00(i)
       h20t(i)=hn20*d2t !
       h20bq(i)=hn20b*rrq(i)

       pt(i)=pt2d(i,k)
       qv(i)=qv2d(i,k)
       qc(i)=qc2d(i,k)
       qr(i)=qr2d(i,k)
       qi(i)=qi2d(i,k)
       qs(i)=qs2d(i,k)
       qg(i)=qg2d(i,k)
       qh(i)=qh2d(i,k)
       if (qc(i) .le.  cmin) qc(i)=0.0
       if (qr(i) .le.  cmin) qr(i)=0.0
       if (qi(i) .le.  cmin) qi(i)=0.0
       if (qs(i) .le.  cmin) qs(i)=0.0
       if (qg(i) .le.  cmin) qg(i)=0.0
       if (qh(i) .le.  cmin) qh(i)=0.0
       tair(i)=(pt(i)+tb0)*pi0(i)
       tairc(i)=tair(i)-t0
       zr(i)=zrr(i)
       zs(i)=zsr(i)
       zg(i)=zgr(i)
       zh(i)=zhr(i)
       vr(i)=0.0 
       vs(i)=0.0
       vg(i)=0.0
       vi(i)=0.0
       vh(i)=0.0

       ftns(i)=1.
       ftng(i)=1.
       ftns0(i)=1.
       ftng0(i)=1.
       fros(i)=1.
       ftnh(i)=1.
       ftnh0(i)=1.

       cnd(i)=0.0
       dep(i)=0.
       ern(i)=0.0
       pint(i)=0.0
       pidep(i)=0.0

       psdep(i)=0.
       pgdep(i)=0.
       dd1(i)=0.
       dd(i)=0.
       pgsub(i)=0.
       psmlt(i)=0.
       pgmlt(i)=0.
       pimlt(i)=0.
       psacw(i)=0.
       piacr(i)=0.

       pssub(i)=0.0
       pgsub(i)=0.0

       psfw(i)=0.0
       psfi(i)=0.0
       pidep(i)=0.0

       pgfr(i)=0.
       psacr(i)=0.
       wgacr(i)=0.
       pihom(i)=0.
       pidw(i)=0.0
          
       psaut(i)=0.0
       psaci(i)=0.0
       praci(i)=0.0
       pwacs(i)=0.0
       qsacw(i)=0.0
          
       pracs(i)=0.0
       qracs(i)=0.0
       qsacr(i)=0.0
       pgaut(i)=0.0
 
       praut(i)=0.0
       pracw(i)=0.0
       pgfr(i)=0.0

       qracs(i)=0.0

       pgacs(i)=0.0
       qgacw(i)=0.0
       dgaci(i)=0.0
       dgacs(i)=0.0
       wgacs(i)=0.0
       wgaci(i)=0.0
       dgacw(i)=0.0
       dgacr(i)=0.
       pgwet(i)=0.0

       qgacr(i)=0.0

       pihom(i)=0.0
       pimlt(i)=0.0
       pidw(i)=0.0
       pimm(i)=0.0
       pcfr(i)=0.0

       pihms(i)=0.0 
       pihmg(i)=0.0
       ftns(i)=1.
       ftng(i)=1.
       pmlts(i)=0.0
       pmltg(i)=0.0

       phfr(i)=0.0
       phmlt(i)=0.0
       dhacw(i)=0.0
       qhacw(i)=0.0
       dhacr(i)=0.0
       qhacr(i)=0.0
       whacr(i)=0.0
       dhaci(i)=0.0
       whaci(i)=0.0
       dhacs(i)=0.0
       phacs(i)=0.0
       whacs(i)=0.0
       dhacg(i)=0.0
       whacg(i)=0.0
       phwet(i)=0.0
       phdep(i)=0.0
       phsub(i)=0.0
       pvapg(i)=0.0
       pvaph(i)=0.0
       primh(i)=0.0

       dlt4(i)=0.0
       dlt3(i)=0.0
       dlt2(i)=0.0

!
!     ***   y1 : dynamic viscosity of air (u)
!     ***   dwv : diffusivity of water vapor in air (pi)
!     ***   tca : thermal conductivity of air (ka)
!     ***   y2 : kinetic viscosity (v)

       y1(i)=c149*tair(i)**1.5/(tair(i)+120.)
       dwv(i)=dwvp(i)*tair(i)**1.81
       tca(i)=c141*y1(i)
       scv(i)=1./((rr0(i)*y1(i))**.1666667*dwv(i)**.3333333)

   enddo  !i-loop


!dir$ vector aligned
   do i=1,irestrict

!
!  calculating processes related to both ice and warm rain
!   compute zr,zs,zg,vr,vs,vg      
!
      if (qr(i) .gt. cmin) then
          dd(i)=r00(i)*qr(i)
          y1(i)=sqrt(dd(i))
          y2(i)=sqrt(y1(i))
          zr(i)=zrc/y2(i)
      endif

      !
      ! ice and rain terminal velocity
      !
      call vqrqi(1,r00(i),fv0(i),qr(i),qc(i),tair(i),vr(i))
      call vqrqi(2,r00(i),fv0(i),qi(i),qc(i),tair(i),vi(i))

      ftns(i)=1.
      ftns0(i)=1.

      qhz2=qh2d(i,k) 
      qgz2=qg2d(i,k) 
      if (k .lt. kmax-2 .and. tairc(i) .ge. -5) then
          qhz2=qh2d(i,k+1) 
          qgz2=qg2d(i,k+1)
      endif
     
      !
      ! PSD mapping
      !
      call sgmap(1,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftns0(i))
      call sgmap(3,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),fros0(i))

      if (qs(i) .gt. cmin) then
          dd(i)=r00(i)*qs(i)
          y1(i)=dd(i)**.25
          ftns(i)=1.
          ftns(i)=ftns0(i)**0.25
          fros(i)=1                     !improve4
          fros(i)=fros0(i)**0.25        !improve4
          zs(i)=zsc/y1(i)*ftns(i)*fros(i)  !improve4
          ftns(i)=ftns0(i)**bsq
          fros(i)=fros0(i)**bsq         !improve4
          vs(i)=max(vscf(i)*dd(i)**bsq/ftns(i)/fros(i), 0.)
      endif

      ftng(i)=1.
      ftng0(i)=1.

      !
      ! PSD mapping
      !
      call sgmap(2,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftng0(i))

      if (qg(i) .gt. cmin) then
          dd(i)=r00(i)*qg(i)
          y1(i)=dd(i)**.25
          ftng(i)=1.
          ftng(i)=ftng0(i)**0.25

          zg(i)=zgc/y1(i)*ftng(i)
          if(dd(i).gt.qrog2) zg(i)=zgc2/y1(i)*ftng(i)

          ftng(i)=ftng0(i)**bgq
          vg(i)=max(vgcf(i)*dd(i)**bgq/ftng(i), 0.0)
          if(dd(i).gt.qrog2)then
              ftng(i)=ftng0(i)**bgq2
              vg(i)=max(vgcf2(i)*dd(i)**bgq2/ftng(i), 0.e0)
          endif !improve4
      endif !qg

      !
      ! PSD mapping
      !
      call sgmap(4,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftnh0(i))

      if (qh(i) .gt. cmin) then
          dd(i)=r00(i)*qh(i)
          y1(i)=dd(i)**.25
          ftnh(i)=ftnh0(i)**0.25
          zh(i)=zhc/y1(i)*ftnh(i)
          ftnh(i)=ftnh0(i)**bhq
          vh(i)=max(vhcr(i)*dd(i)**bhq/ftnh(i), 0.0)
      endif

      if (qr(i) .le. cmin1) vr(i)=0.0
      if (qs(i) .le. cmin1) vs(i)=0.0
      if (qg(i) .le. cmin1) vg(i)=0.0
      if (qi(i) .le. cmin1) vi(i)=0.0
      if (qh(i) .le. cmin1) vh(i)=0.0

!
!     ***   y1 : dynamic viscosity of air (u)
!     ***   dwv : diffusivity of water vapor in air (pi)
!     ***   tca : thermal conductivity of air (ka)
!     ***   y2 : kinetic viscosity (v)
!

!
!*  1 * psaut : autoconversion of qi to qs                        ***1**
!*  3 * psaci : accretion of qi to qs                             ***3**
!*  4 * psacw : accretion of qc by qs (riming) (qsacw for psmlt)  ***4**
!*  5 * praci : accretion of qi by qr                             ***5**
!*  6 * piacr : accretion of qr or qg by qi                       ***6**
!* 34 * pwacs : collection of qs by qc                            **34**
!

      pihms(i)=0.0 
      pihmg(i)=0.0 
      psaut(i)=0.0
      psaci(i)=0.0
      praci(i)=0.0
      piacr(i)=0.0
      psacw(i)=0.0
      pwacs(i)=0.0
      qsacw(i)=0.0
      ftns(i)=1.
      ftng(i)=1.
      ftnh(i)=1.
      ftns(i)=ftns0(i)
      ftng(i)=ftng0(i)

      if (tair(i).lt.t0) then

          rn1=1./300.
          bnd1=6.e-5
          esi(i)=0.25
          psaut(i)=r2is*max(rn1*esi(i)*(qi(i)-bnd1*fv0(i)*fv0(i)) ,0.0) 
          ftns(i)=ftns0(i)
          ftng(i)=ftng0(i)
          fros(i)=fros0(i)
          ftnh(i)=ftnh0(i)
          esi(i)=1.0 !esi constant in consatrh via r3f/rn3; use esi( ) to make f(t)/f(q)
          dmicrons=(r00(i)*qs(i)/(roqs*fros(i))/cpi/(tns*ftns(i)))**.25*1.e4
          esi(i)=min(1.,(dmicrons/375.)**3.) ! f(dmicrons)
 
          y1(i)=1.0
          if (vs(i).gt.0.) y1(i)=abs((vs(i)-vi(i)) /vs(i))
          psaci(i)=r2is*y1(i)*r3f(i)*qi(i)/zs(i)**bs3*ftns(i)*esi(i)
          if(qs(i).le.cmin) psaci(i)=0.
          psacw(i)=r2is*r4f(i)*qc(i)/zs(i)**bs3*ftns(i)
          if(qs(i).le.cmin) psacw(i)=0.
          if (ihalmos.eq.1)then
              y2(i)=0.
              if((tairc(i).le.hmtemp1).and.(tairc(i).ge.hmtemp4)) y2(i)=0.5
              if((tairc(i).le.hmtemp2).and.(tairc(i).ge.hmtemp3)) y2(i)=1.
              pihms(i)=r2ih*psacw(i)*y2(i)*xnsplnt*1000.*xmsplnt
              psacw(i)=psacw(i)-pihms(i)
          endif
          pwacs(i)=r2is*r34f(i)*qc(i)/zs(i)**bs6*ftns(i)*fros(i)      !improve4

          if(qs(i).le.cmin) pwacs(i)=0.
             y1(i)=1./zr(i)
             y2(i)=y1(i)*y1(i)
             y3(i)=y1(i)*y2(i)
             y5(i)=1.0
          if (vr(i).gt.0.) y5(i)=abs((vr(i)-vi(i)) /vr(i))
          dd(i)=y5(i)*r5f(i)*qi(i)*y3(i)*(rn50+rn51*y1(i)+rn52*y2(i)+rn53*y3(i))
          praci(i)=max(dd(i),0.0)
          if (qr(i) .le. cmin) praci(i)=0.
          y4(i)=y3(i)*y3(i)
          dd1(i)=y5(i)*r6f(i)*qi(i)*y4(i)*(rn60+rn61*y1(i) +rn62*y2(i)+rn63*y3(i))
          piacr(i)=max(dd1(i),0.0)
          if (qr(i) .le. cmin) piacr(i)=0.
      else
          qsacw(i)=r2is*r4f(i)*qc(i)/zs(i)**bs3*ftns(i)
          if (qs(i) .le. cmin) psacw(i)=0.
      endif !tairc 

!
!* 21 * praut   autoconversion of qc to qr                        **21**
!* 22 * pracw : accretion of qc by qr                             **22**
!


#if ( WRF_CHEM == 1)
      !
      ! WRF-CHEM GOCART coupling option
      !
      ! emk...only execute when gocart and coupling are selected
      if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
            chem_opt == 302 .or. chem_opt == 303) .and. &
            (gsfcgce_gocart_coupling == 1) ) then
          !      sat vapor pressure [hpa,mb]
          e_sat = 6.11 * exp( 5423. *( 1.0/273.15 - 1./tair(i) ) )
          !      dry air vapor pressure [hpa, mb]
          e_dry = qv(i) / ( qv(i) + 0.622 ) * p2d(i,k) * 1.e-3  ! p0 in [g*cm/s2/cm2]
          rh_rad = max(1.e-6, e_dry/e_sat*100.)  ! relative humidity [%]
          super_sat = max(0.001, rh_rad - 100.e0) !super saturation [%]
          !
          ! convert gocart aerosol mass conc to cn
          !
          call mass2ccn(tair(i),super_sat,aero(i,k,:),ccn_out(i) )
          !             nc_cgs(i,k)  = max(100., ccn_out)  !diagnostic cloud droplet conc  [#/cm3]
          ccn_out(i) = max(100., ccn_out(i))

          ! rho_dryair = p0(i,k) / ( tair(i) * 2.87 * 1.0e6) ! dry air density [g/cm3]
          l_cloud = qc(i) * rho2d(i,k)             ! cloud water [g/cm3]
          !  g/g        g/cm3
          !             call auto_conversion_liu_daum( l_cloud, nc_cgs(i,k), p_liu_daum, re_liu_daum )
          call auto_conversion_liu_daum( l_cloud, ccn_out(i), p_liu_daum, re_liu_daum )

          praut(i) = p_liu_daum / rho2d(i,k)  !qc2qr autoconversion rate [g/g s-1]
      else
          praut(i)=max(rn21*(qc(i)-bnd21),0.0)  !qc2qr autoconversion rate [g/g s-1]
      end if ! if (gsfcgce_gocart_coupling == 1)
#else
      praut(i)=max(rn21*(qc(i)-bnd21),0.0)  !qc2qr autoconversion rate [g/g s-1]
#endif

      y1(i)=1./zr(i)
      y2(i)=y1(i)*y1(i)
      y3(i)=y1(i)*y2(i)
      y4(i)=r22f(i)*qc(i)*y3(i)*(rn50+rn51*y1(i)+rn52*y2(i)+rn53*y3(i))
      pracw(i)=max(y4(i), 0.0) 
      if(qr(i) .le. cmin) pracw(i)=0.

!* 12 * psfw : bergeron processes for qs (koening, 1971)          **12**
!* 13 * psfi : bergeron processes for qs                          **13**

      pidep(i)=0.0
      psfw(i)=0.0
      psfi(i)=0.0

   enddo

!dir$ vector aligned
   do i=1,irestrict

      if (tair(i) .lt. t0) then
          y1(i)=max( min(tairc(i), -1.), -31.)
          it(i)=int(abs(y1(i)))
          y1(i)=rn12a(it(i))
          y2(i)=rn12b(it(i))
          y3(i)=rn13(it(i))
          psfw(i)=r2is*max(d2t*y1(i)*(y2(i)+r12r(i)*qc(i))* qi(i), 0.0)
          psfi(i)=r2is*y3(i)*qi(i)
!
          y4(i)=1./(tair(i)-c358)
          y5(i)=1./(tair(i)-c76)
          qsw(i)=rp0(i)*exp(c172-c409*y4(i))
          qsi(i)=rp0(i)*exp(c218-c580*y5(i))
          hfact=(qv(i)+qb0-qsi(i))/(qsw(i)-qsi(i)+cmin1)
          ! add cmin1 to prevent overflow of hfact
          if(hfact.gt.1.) hfact=1.
          sfact=1
          ssi(i)=(qv(i)+qb0)/qsi(i)-1.

          fssi=min(xssi,max(rssi,xssi*(tairc(i)+44.)/(44.0-38.0)))  !improve4 max ssi f(tair)
          fssi=rssi
          wssi=w2d(i,k)/100.-2.0
          if(wssi.gt.0.) fssi=rssi+min(xssi,wssi*0.01)
          fssi=min(ssi(i),fssi)

!  steve : please check
          !
          ! heterogeneous ice nucleation [Meyers et al. ]
          !
          if (tairc(i).le.-5.) then                         !meyers
              r_nci=max(1.e-3*exp(-.639+12.96*fssi),0.528e-3)  !meyers et al. 
          else
              r_nci=min(1.e-3*exp(-.639+12.96*fssi),0.528e-3)  !meyers et al. 
          endif  !tairc

          if (r_nci.gt.15.) r_nci=15.                  !cap at 15000/liter

          !
          ! heterogeneous ice nucleation [Cooper 
          !
          if( tairc(i) .lt. -40.0 ) then
              c_nci = 5.0e-6*exp(0.304*40.0)
          else
              c_nci = 5.0e-6*exp(0.304*abs(tairc(i)))
          end if
          r_nci = c_nci
          ! jdc cap r_nci to the amount corresponding to crystal size of ami50
          r_nci = max(r_nci,r00(i)*qi(i)/ami50)

          !
          ! heterogeneous ice nucleation [DeMott et al. ]
          ! only execute when gocart and coupling is turned on.
#if ( WRF_CHEM == 1)
          if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
                chem_opt == 302 .or. chem_opt == 303) .and. &
                (gsfcgce_gocart_coupling == 1) ) then
                ! conversion rate of cloud water to ice in the bergeron porcess based on meyer + demott formulae
                !
                ! convert gocart aerosol mass conc to in
                !      p0 need to be converted from g*cm/s2/cm2 to mb (hpa)
                !                call mass2icn(p0(i,k)*0.001,tair(i),aero(i,k,:), icn_out)
!                  call mass2icn(p0(i,k)*0.001,tair(i),aero(i,k,:), icn_out,i,k)
                ! emk...mass2icn requires i,j,k indices, but only prints
                ! the values when an error occurs.  in no case are the i,j,k
                ! used to access a value in an array.  since mass2icn is
                ! shared with the older 3ice scheme which still uses j,
                ! we will simply pass a bogus value of j=1 in the call
                ! below.
                call mass2icn(p2d(i,k)*0.001,tair(i),aero(i,k,:), icn_out(i), i,1,k)

                icn_out(i) = min(1.e3, max(0.01e0 ,  icn_out(i)) )
                r_nci = icn_out(i) * 1.e-3  !demotto's formuale
          end if ! if (gsfcgce_gocart_coupling == 1)
#endif

          dd(i)=min((r00(i)*qi(i)/r_nci), ami40)   !mean cloud ice mass
          yy1=1.-aa2(it(i))                   
          sfact=(ami50**yy1-ami40**yy1)/(ami50**yy1-dd(i)**yy1)

          ! water vapor diffusivity correction term
          esi(i)  = qsi(i)/rp0(i)*c610
          y3(i)   = 1./tair(i)
          term1     = y3(i)*(rn20a*y3(i)-rn20b)
          term2     = rn10c*tair(i)/esi(i)
          dd(i)   = term1+term2
          fdwv      = dd(i)/(term1+term2*dwv0/dwv(i))
          psfw(i) = max( d2t*y1(i)*fdwv*(y2(i)*fdwv + r12r(i)*qc(i))*qi(i),0.0 )

          psfw(i)=0.0

          if (hfact.gt.0.) then
              psfi(i)=r2is*psfi(i)*hfact*sfact*fdwv
          else
              psfi(i)=0.
          endif !hfact

          if(qi(i).le.1.e-5*fv0(i)*fv0(i)) psfi(i)=0.0

      endif !tair

   enddo !i-loop


!dir$ vector aligned
   do i=1,irestrict

!
!ttt***** qg=qg+min(pgdry,pgwet)
!*  9 * pgacs : accretion of qs by qg (dgacs,wgacs: dry and wet)  ***9**
!* 14 * dgacw : accretion of qc by qg (qgacw for pgmlt)           **14**
!* 16 * dgacr : accretion of qr to qg (qgacr for pgmlt)           **16**
!*******pgdry : dgacw+dgaci+dgacr+dgacs                           ******
!* 15 * dgaci : accretion of qi by qg (wgaci for wet growth)      **15**
!* 17 * pgwet : wet growth of qg                                  **17**
!*  steve turned off pgwet, set pgwet = 0.
!*  steve turned off wgaci, set wgaci = 0.                           
!

      y1(i)=abs( vg(i)-vs(i) )
      y2(i)=zs(i)*zg(i)
      y3(i)=5./y2(i)
      y4(i)=.08*y3(i)*y3(i)
      y5(i)=.05*y3(i)*y4(i)
      y2(i)=y1(i)*(y3(i)/zs(i)**5+y4(i)/zs(i)**3 +y5(i)/zs(i))

      pgacs(i)=r2ig*r2is*r9rf(i)*y2(i)*ftns(i)*ftng(i)*fros(i) !improve4
      if(qs(i).le.cmin) pgacs(i)=0.
      if(qg(i).le.cmin) pgacs(i)=0.
      dgacs(i)=pgacs(i)
      dgacs(i)=0.0             !lang et al. 2007
!crh     wgacs(i)=10.*r9rf(i)*y2(i)
      wgacs(i)=0.0
      wgacs(i)=10.*r9rf(i)*y2(i)*ftns(i)*ftng(i)*fros(i)         !4ice revised
!        if(r00(i)*qg(i).gt.qrog2)then                                      !4ice revised
!         wgacs(i)=10.*r9rf(i)*y2(i)*ftns(i)*ftng(i)*fros(i)        !4ice revised
!        endif                                                             !4ice revised
      if(qg(i).le.cmin) wgacs(i)=0.                !4ice revised2
      if(qs(i).le.cmin) wgacs(i)=0.                !4ice revised2

      dhacs(i)=0.0
      whacs(i)=0.0
      y1(i)=abs( vh(i)-vs(i) )
      y2(i)=zs(i)*zh(i)
      y3(i)=5./y2(i)
      y4(i)=.08*y3(i)*y3(i)
      y5(i)=.05*y3(i)*y4(i)
      dd(i)=y1(i)*(y3(i)/zs(i)**5+y4(i)/zs(i)**3+y5(i)/zs(i))
      whacs(i)=r2ih*r2is*min(h9r(i)*dd(i)*ftnh(i)*ftns(i)*fros(i),qs(i)/d2t)
      if(qs(i).le.cmin) whacs(i)=0.
      if(qh(i).le.cmin) whacs(i)=0.

      dhacg(i)=0.0
      whacg(i)=0.0
      y1(i)=abs( vh(i)-vg(i) )
      y2(i)=zg(i)*zh(i)
      y3(i)=5./y2(i)
      y4(i)=.08*y3(i)*y3(i)
      y5(i)=.05*y3(i)*y4(i)
      dd(i)=y1(i)*(y3(i)/zg(i)**5+y4(i)/zg(i)**3+y5(i)/zg(i))
      whacg(i)=r2ih*r2ig*min(h10r(i)*dd(i)*ftnh(i)*ftng(i),qg(i)/d2t)
      if(r00(i)*qg(i).gt.qrog2) whacg(i)=whacg(i)/roqg*roqg2*0.5        !reduce ehg for high dens grp
      if(qg(i).le.cmin) whacg(i)=0.
      if(qh(i).le.cmin) whacg(i)=0.

      y1(i)=1./zg(i)**bg3
      esi(i)=1.0 !egc constant in consatrh via r14f/rn14; use esi( ) to make f(t)/f(q)
      dmicrong=(r00(i)*qg(i)/roqg/cpi/(tng*ftng(i)))**.25*1.e4
      if(r00(i)*qg(i).gt.qrog2)then
         y1(i)=1./zg(i)**bg3_2
         dmicrong=(r00(i)*qg(i)/roqg2/cpi/(tng*ftng(i)))**.25*1.e4
      endif
      esi(i)=min(1.,(dmicrong/500.)**1.1)       ! f(dmicrons)

      dgacw(i)= r2ig*esi(i)*r14f(i)*qc(i)*y1(i)*ftng(i)
      if(r00(i)*qg(i).gt.qrog2) dgacw(i)=r2ig*esi(i)*r14f2(i)*qc(i)*y1(i)*ftng(i)
      if(qg(i).le.cmin) dgacw(i)=0.

      dhacw(i)=0.0                                                  !4ice
      y2(i)=1./zh(i)**bh3                                           !4ice
      dhacw(i)=r2ih*max(h14r(i)*qc(i)*y2(i)*ftnh(i), 0.0)           !4ice
      if(qh(i).le.cmin) dhacw(i)=0.
        
      if(ihalmos.eq.1)then
         y2(i)=0.
         if ((tairc(i).le.hmtemp1).and.(tairc(i).ge.hmtemp4)) y2(i)=0.5
         if ((tairc(i).le.hmtemp2).and.(tairc(i).ge.hmtemp3)) y2(i)=1.
         pihmg(i)=r2ig*dgacw(i)*y2(i)*xnsplnt*1000.*xmsplnt
         dgacw(i)=r2ig*(dgacw(i)-pihmg(i))
         pihmh(i)=0.0                                                   !4ice
         pihmh(i)=r2ih*dhacw(i)*y2(i)*xnsplnt*1000.*xmsplnt             !4ice
         dhacw(i)=r2ih*(dhacw(i)-pihmh(i))                              !4ice
      endif !ihalmos

      qgacw(i)=r2ig*dgacw(i)
      qhacw(i)=r2ih*dhacw(i)                                          !4ice
      y1(i)=1./zg(i)**bg3
      y5(i)=1.0
      if (vg(i).gt.0.) y5(i)=abs((vg(i)-vi(i))/vg(i))
      dgaci(i)= r2ig*y5(i)*r15f(i)*qi(i)*y1(i)*ftng(i)
      if(qg(i).le.cmin) dgaci(i)=0.
      dgaci(i)=0.0
      wgaci(i)=0.0
      wgaci(i)=y5(i)*r15af(i)*qi(i)*y1(i)*ftng(i)                !4ice revised
      if(r00(i)*qg(i).gt.qrog2)then                                      !4ice revised
         y1(i)=1./zg(i)**bg3_2
         wgaci(i)=y5(i)*r15af2(i)*qi(i)*y1(i)*ftng(i)              !4ice revised
      endif                                                             !4ice revised
      if(qg(i).le.cmin) wgaci(i)=0.                !4ice revised2

      dhaci(i)=0.0                                                     !4ice
      whaci(i)=0.0                                                     !4ice
      y5(i)=1.0                                                      !4ice
      if(vh(i).gt.0.) y5(i)=abs((vh(i)-vi(i))/vh(i))         !4ice
      y2(i)=1./zh(i)**bh3                                          !4ice
      whaci(i)=r2ih*min(y5(i)*h15ar(i)*qi(i)*y2(i)*ftnh(i),qi(i)/d2t)
      if(qh(i).le.cmin) whaci(i)=0.

      y1(i)=abs( vg(i)-vr(i) )
      y2(i)=zr(i)*zg(i)
      y3(i)=5./y2(i)
      y4(i)=.08*y3(i)*y3(i)
      y5(i)=.05*y3(i)*y4(i)
      dd(i)=r16rf(i)*y1(i)*(y3(i)/zr(i)**5+y4(i)/zr(i)**3 +y5(i)/zr(i))*ftng(i)
      dgacr(i)=r2ig*max(dd(i),0.0)
      if(qg(i).le.cmin) dgacr(i)=0.
      if(qr(i).le.cmin) dgacr(i)=0.
      qgacr(i)=dgacr(i)

      pracg(i)=0.
      qracg(i)=0.
      y2(i)=zr(i)*zg(i)
      y3(i)=5./y2(i)
      y4(i)=.08*y3(i)*y3(i)
      y5(i)=.05*y3(i)*y4(i)
      pracg(i)=r16rf(i)/roqr*roqg*y1(i)*(y3(i)/zg(i)**5+y4(i) &
                   /zg(i)**3+y5(i)/zg(i))*ftng(i)
      if(r00(i)*qg(i).gt.qrog2) pracg(i)=pracg(i)/roqg*roqg2
      if(qg(i).le.cmin) pracg(i)=0.
      if(qr(i).le.cmin) pracg(i)=0.
      qracg(i)=min(d2t*pracg(i), qg(i))

      y1(i)=abs( vh(i)-vr(i) )
      y2(i)=zr(i)*zh(i)
      y3(i)=5./y2(i)
      y4(i)=.08*y3(i)*y3(i)
      y5(i)=.05*y3(i)*y4(i)
      dd(i)=h16r(i)*y1(i)*ftnh(i)*(y3(i)/zr(i)**5   &
		+y4(i)/zr(i)**3+y5(i)/zr(i))
      dhacr(i)=r2ih*max(dd(i), 0.0)
      if(qh(i).le.cmin) dhacr(i)=0.
      if(qr(i).le.cmin) dhacr(i)=0.
      qhacr(i)=dhacr(i)

      if (tair(i) .ge. t0) then
          dgacs(i)=0.0
          wgacs(i)=0.0   !4ice revised
          whacs(i)=0.0
          whacg(i)=0.0
          dgacw(i)=0.0
          dhacw(i)=0.0                                                   !ice4
          dgaci(i)=0.0
          wgaci(i)=0.0   !4ice revised
          whaci(i)=0.0                                                   !ice4
          dgacr(i)=0.0
          pracg(i)=0.0
          dhacr(i)=0.0
      else
          pgacs(i)=0.0
          qgacw(i)=0.0
          qhacw(i)=0.0                                                   !ice4
          qgacr(i)=0.0
          qracg(i)=0.0
          qhacr(i)=0.0
      endif
 
      pgwet(i)=0.0


      if(tair(i) .lt. t0) then   !4ice revised

         ! nuwrf corrected by steve and roger to prevent division by 0
         y1(i)=1./(alf+rn17c*max(tairc(i),-75.0)) 
         y2(i)=r191r(i)/zg(i)**2+r192rf(i)/zg(i)**bgh5                   !4ice revised chern
         if(r00(i)*qg(i).gt.qrog2)      &                                 !4ice revised chern
            y2(i)=(r191r(i)/zg(i)**2+r192rf2(i)/zg(i)**bgh5_2)              !4ice revised chern
         y4(i)=alvr(i)*dwv(i)*(rp0(i)-(qv(i)+qb0))-tca(i)*tairc(i)   !4ice revised
         dd(i)=y1(i)*(rn20*ftng(i)*y4(i)*y2(i)+(wgaci(i) &    !4ice revised chern
               +wgacs(i))*(alf+rn17b*tairc(i)))                       !4ice revised
         pgwet(i)=max(dd(i), 0.0)                                    !4ice revised
         if(qg(i).le.cmin) pgwet(i)=0.                               !4ice revised
      endif                                                              !4ice revised

      phwet(i)=0.0
      if (tair(i) .lt. t0) then
          y1(i)=1./(alf+rn17c*tairc(i))
          y3(i)=.78/zh(i)**2+h17aq(i)*scv(i)/zh(i)**bhh5
          y4(i)=alvr(i)*dwv(i)*(rp0(i)-(qv(i)+qb0))-tca(i)*tairc(i)
          dd(i)=y1(i)*(h17r(i)*ftnh(i)*y4(i)*y3(i)+(whaci(i)+whacs(i)   &
                +whacg(i))*(alf+rn17b*tairc(i)))
          phwet(i)=r2ih*max(dd(i), 0.0)
          if(qh(i).le.cmin) phwet(i)=0.
      endif !tair

  enddo


!dir$ vector aligned
  do i=1,irestrict

     !
     !  handling the negative cloud water (qc) 
     !
     y1(i)=qc(i)/d2t
     psacw(i)=min(y1(i), psacw(i))
     pihms(i)=min(y1(i), pihms(i))
     praut(i)=min(y1(i), praut(i))
     pracw(i)=min(y1(i), pracw(i))
     psfw(i)= min(y1(i), psfw(i))
     dgacw(i)=min(y1(i), dgacw(i))
     pihmg(i)=min(y1(i), pihmg(i))
     dhacw(i)=min(y1(i), dhacw(i))                              !4ice
     pihmh(i)=min(y1(i), pihmh(i))                              !4ice
     qsacw(i)=min(y1(i), qsacw(i))
     qgacw(i)=min(y1(i), qgacw(i))
     qhacw(i)=min(y1(i), qhacw(i))                              !4ice

     y1(i)=d2t*(psacw(i)+praut(i)+pracw(i)+psfw(i)            &
           +dgacw(i)+qsacw(i)+qgacw(i)+pihms(i)+pihmg(i)     &
           +dhacw(i)+qhacw(i)+pihmh(i))                          !4ice

     qc(i)=qc(i)-y1(i)

     if (qc(i) .lt. 0.0) then
         y2(i)=1.
         if (y1(i) .ne. 0.) y2(i)=qc(i)/y1(i)+1.
         psacw(i)=psacw(i)*y2(i)
         praut(i)=praut(i)*y2(i)
         pracw(i)=pracw(i)*y2(i)
         psfw(i)=psfw(i)*y2(i)
         dgacw(i)=dgacw(i)*y2(i)
         dhacw(i)=dhacw(i)*y2(i)
         qsacw(i)=qsacw(i)*y2(i)
         qgacw(i)=qgacw(i)*y2(i)
         qhacw(i)=qhacw(i)*y2(i)
         pihms(i)=pihms(i)*y2(i)
         pihmg(i)=pihmg(i)*y2(i)
         pihmh(i)=pihmh(i)*y2(i)                                   !4ice
         qc(i)=0.0
     endif

     !                                                                          !4ice revised
     ! convert wet graupel to hail                                       !4ice revised
     !                                                                          !4ice revised
     y1(i)=dgacw(i)+dgacr(i)                                     !4ice revised
     if(y1(i).lt..95*pgwet(i).or.y1(i).eq.0.)then                !4ice revised
        wgaci(i)=0.0                                                  !4ice revised
        wgacs(i)=0.0                                                  !4ice revised
     endif                                                             !4ice revised
     pg2h(i)=0.0                                                     !4ice revised

     if(y1(i).gt.0..and.pgwet(i).gt.0..and. y1(i).gt.1.0*pgwet(i))then
!test     pg2h(i)=2.0*y1(i)                                            !4ice revised
        pg2h(i)=qg(i)/d2t                                            !4ice revised
!         pg2h(i)=qg(i)*0.647232/d2t                                   !4ice revised deff
        pg2h(i)=min(pg2h(i), qg(i)/d2t)                            !4ice revised
     endif                                                             !4ice revised

     whacr(i)=phwet(i)-dhacw(i)-whaci(i)-whacs(i)-whacg(i)    !4ice
     y2(i)=dhacw(i)+dhacr(i)                                        !4ice
     if(y2(i).lt.0.95*phwet(i).or.y2(i).eq.0.) then
        whacr(i)=0.0                                                    !4ice
        whaci(i)=0.0                                                    !4ice 
        whacg(i)=0.0
        whacs(i)=0.0
     endif
          
     primh(i)=0.0
     if(y2(i).gt.0..and.phwet(i).gt.0..and. y2(i).lt..95*phwet(i)) then

        rim_frac(i)=2.0*min((1.-y2(i)/phwet(i))**2,1.0)
        rim_frac(i)=rim_frac(i)*min((tairc(i)/(t00-t0))**2,1.0)
!       primh(i)=rim_frac(i)*y2(i)
        primh(i)=rim_frac(i)*dhacw(i)
        primh(i)=min(primh(i), qh(i)/d2t)
     endif

   enddo  !i-loop


!dir$ vector aligned
   do i=1,irestrict

      !
      ! handling the negative cloud ice (qi)   
      !
      y1(i)=qi(i)/d2t
      psaut(i)=min(y1(i), psaut(i))
      psaci(i)=min(y1(i), psaci(i))
      praci(i)=min(y1(i), praci(i))
      psfi(i)= min(y1(i), psfi(i))
      dgaci(i)=min(y1(i), dgaci(i))
      wgaci(i)=min(y1(i), wgaci(i))
      whaci(i)=min(y1(i), whaci(i))                             !4ice

      qi(i)=qi(i)+d2t*(pihms(i)+pihmg(i)+pihmh(i))             !fix sel

      y1(i)=d2t*(psaut(i)+psaci(i)+praci(i)+psfi(i)       &
                    +dgaci(i)+wgaci(i)+whaci(i))      

      qi(i)=qi(i)-y1(i)

      if (qi(i) .lt. 0.0) then
          y2(i)=1.
          if (y1(i) .ne. 0.0) y2(i)=qi(i)/y1(i)+1.
          psaut(i)=psaut(i)*y2(i)
          psaci(i)=psaci(i)*y2(i)
          praci(i)=praci(i)*y2(i)
          psfi(i)=psfi(i)*y2(i)
          dgaci(i)=dgaci(i)*y2(i)
          wgaci(i)=wgaci(i)*y2(i)
          whaci(i)=whaci(i)*y2(i)                                   !4ice
          qi(i)=0.0
      endif

      wgacr(i)=qgacr(i)+qgacw(i)
      dlt3(i)=0.0
      if (qr(i) .lt. 1.e-4) dlt3(i)=1.
      dlt4(i)=1.

      if (qc(i) .gt. 5.e-4) dlt4(i)=0.0

      if (qs(i) .le. 1.e-4) dlt4(i)=1.

      if (tair(i) .ge. t0) then
          dlt3(i)=0.0
          dlt4(i)=0.0
      endif
         
      pr(i)=d2t*(qsacw(i)+praut(i)+pracw(i)+wgacr(i) +qhacw(i)-qgacr(i))
      ps(i)=d2t*(psaut(i)+psaci(i)+dlt4(i)*psacw(i) +psfw(i)+psfi(i))
      pg(i)=d2t*(dlt3(i)*praci(i)+dgaci(i)               &
                +wgaci(i)+dgacw(i)+(1.-dlt4(i))*psacw(i)+primh(i)) !4ice revised

!*  7 * pracs : accretion of qs by qr                             ***7**
!*  8 * psacr : accretion of qr by qs (qsacr for psmlt)           ***8**
!*  2 * pgaut : autoconversion of qs to qg                        ***2**
!* 18 * pgfr : freezing of qr to qg                               **18**

      qracs(i)=0.0
      y1(i)=abs(vr(i)-vs(i))
      y2(i)=zr(i)*zs(i)
      y3(i)=5./y2(i)
      y4(i)=.08*y3(i)*y3(i)
      y5(i)=.05*y3(i)*y4(i)

      pracs(i)=r2ig*r2is*(r7rf(i)*y1(i)*(y3(i)/zs(i)**5 &            !improve4
               +y4(i)/zs(i)**3+y5(i)/zs(i))*ftns(i)*fros(i))       !improve4
      if(qs(i).le.cmin) pracs(i)=0.
      if(qr(i).le.cmin) pracs(i)=0.
      qracs(i)=r2ig*r2is*min(d2t*pracs(i), qs(i))
      psacr(i)=r2is*(r8rf(i)*y1(i)*(y3(i)/zr(i)**5 &
               +y4(i)/zr(i)**3+y5(i)/zr(i))*ftns(i))
      if(qs(i).le.cmin) psacr(i)=0.
      if(qr(i).le.cmin) psacr(i)=0.
      qsacr(i)=psacr(i)

      if (tair(i) .ge. t0) then
          pracs(i)=0.0
          psacr(i)=0.0
      else
          qsacr(i)=0.0
          qracs(i)=0.0
      endif

!* 18 * pgfr : freezing of qr to qg                               **18**
      pgaut(i)=0.0
      pgfr(i)=0.0
      phfr(i)=0.0
      if (tair(i) .lt. t0) then
          y1(i)=exp(rn18a*(t0-tair(i)))
          if( qr(i).ge.0.)then
             temp = 1./zr(i)
             temp = temp*temp*temp*temp*temp*temp*temp
             phfr(i)=r2ih*max(r18r(i)*(y1(i)-1.)*temp, 0.0)
!            phfr(i)=r2ih*max(r18r(i)*(y1(i)-1.)/zr(i)**7., 0.0)
          else
             temp = 1./zr(i)
             temp = temp*temp*temp*temp*temp*temp*temp
             pgfr(i)=r2ig*max(r18r(i)*(y1(i)-1.)*temp, 0.0)
!            pgfr(i)=r2ig*max(r18r(i)*(y1(i)-1.)/zr(i)**7., 0.0)
          if(qr(i).le.cmin) pgfr(i)=0.
          endif
      endif

      !
      ! handling the negative rain water (qr)    
      ! handling the negative snow (qs)       
      !
      y1(i)=qr(i)/d2t
      y2(i)=-qh(i)/d2t
      piacr(i)=min(y1(i), piacr(i))
      dgacr(i)=min(y1(i), dgacr(i))
      dhacr(i)=min(y1(i), dhacr(i))
!     whacr(i)=min(y1(i), whacr(i))
!     whacr(i)=max(y2(i), whacr(i))
      psacr(i)=min(y1(i), psacr(i))
      pgfr(i)= min(y1(i), pgfr(i))
      phfr(i)= min(y1(i), phfr(i))
      del=0.
      if(whacr(i) .lt. 0.) del=1.
      if(del.eq.1) whacr(i)=max(whacr(i),-dhacw(i))    !fix from jdc
      dhacr(i)=min((1.-del)*whacr(i),dhacr(i))
      if(del.eq.1) dhacw(i)=dhacw(i)+del*whacr(i)      !fix from jdc
      y1(i)=(piacr(i)+dgacr(i)+dhacr(i)+psacr(i)+pgfr(i) +phfr(i))*d2t
      qr(i)=qr(i)+pr(i)+qracs(i)-del*whacr(i)*d2t+qracg(i) !fix sel
      qr(i)=qr(i)-y1(i)                                          !fix sel

      if (qr(i) .lt. 0.0) then
          y2(i)=1.
          if (y1(i) .ne. 0.0) y2(i)=qr(i)/y1(i)+1.
          piacr(i)=piacr(i)*y2(i)
          dgacr(i)=dgacr(i)*y2(i)
          dhacr(i)=dhacr(i)*y2(i)
!         if(whacr(i).gt.0.) whacr(i)=whacr(i)*y2(i)
          pgfr(i)=pgfr(i)*y2(i)
          phfr(i)=phfr(i)*y2(i)
          psacr(i)=psacr(i)*y2(i)
          qr(i)=0.0
      endif !qr

      dlt2(i)=1.
      if (qr(i) .gt. 1.e-4) dlt2(i)=0.
      if (tair(i) .ge. t0) dlt2(i)=0.
      y1(i)=qs(i)/d2t
      pgacs(i)=min(y1(i), pgacs(i))
      dgacs(i)=min(y1(i), dgacs(i))
      wgacs(i)=min(y1(i), wgacs(i))
      whacs(i)=min(y1(i), whacs(i))
      pracs(i)=min(y1(i), pracs(i))
      pwacs(i)=min(y1(i), pwacs(i))

      prn(i)=d2t*(dlt3(i)*piacr(i)+dlt2(i)*dgacr(i)+pgfr(i) +dlt2(i)*psacr(i))
      pracs(i)=(1.-dlt2(i))*pracs(i)
      pwacs(i)=(1.-dlt4(i))*pwacs(i)
      pracg(i)=(1.-dlt2(i))*pracg(i)
      psn(i)=d2t*(pgacs(i)+dgacs(i)+wgacs(i)+pracs(i)+pwacs(i)+whacs(i))
 
      qs(i)=qs(i)+ps(i)-qracs(i)-psn(i)
      if (qs(i) .lt. 0.0) then
          y2(i)=1.
          if (psn(i) .ne. 0.) y2(i)=qs(i)/psn(i)+1.
          pgacs(i)=pgacs(i)*y2(i)
          dgacs(i)=dgacs(i)*y2(i)
          wgacs(i)=wgacs(i)*y2(i)
          whacs(i)=whacs(i)*y2(i) 
          pracs(i)=pracs(i)*y2(i)
          pwacs(i)=pwacs(i)*y2(i)
          qs(i)=0.0
      endif
      psn(i)=d2t*(pgacs(i)+dgacs(i)+wgacs(i)+pracs(i)+pwacs(i))
      qg(i)=qg(i)+pg(i)+prn(i)+psn(i)-qracg(i)
      qg(i)=qg(i)-d2t*(pracs(i)+whacg(i)+pracg(i)+pg2h(i))

      if (qg(i) .lt. 0.0) then                              !fix from jdc
          y2(i)=1.                                             !fix from jdc
          if (whacg(i)+pracg(i)+pg2h(i).ne. 0.)   &                 !fix from jdc
             y2(i)=qg(i)/(d2t*(whacg(i)+pracg(i)+pg2h(i)))+1.  !fix from jdc & sel
          whacg(i)=whacg(i)*y2(i)                         !fix from jdc
          pracg(i)=pracg(i)*y2(i)                         !fix from jdc
          pg2h(i)=pg2h(i)*y2(i)                           !fix from sel
          qg(i)=0.0                                           !fix from jdc
      endif        

      qh(i)=qh(i)+d2t*(phfr(i)+(1.-dlt3(i))*piacr(i)        &
            +(1.-dlt3(i))*praci(i)+(1.-dlt2(i))*psacr(i)     &
            +pracs(i)+(1.-dlt2(i))*dgacr(i)+pracg(i)         &
            +dhacw(i)+dhacr(i)-primh(i)+whaci(i)             &
            +whacs(i)+whacg(i)+pg2h(i))

      y1(i)=d2t*(psacw(i)+psfw(i)+dgacw(i)+piacr(i)             &
                +dgacr(i)+psacr(i)+pgfr(i)+phfr(i)+pihms(i)        &
                +pihmg(i)+dhacw(i)+dhacr(i)           +pihmh(i))     &
                -qracs(i)
      pt(i)=pt(i)+afcp(i)*y1(i)
      tair(i)=(pt(i)+tb0)*pi0(i)

!* 11 * psmlt : melting of qs                                     **11**
!* 19 * pgmlt : melting of qg to qr                               **19**

      psmlt(i)=0.0
      pgmlt(i)=0.0
      phmlt(i)=0.0

      qhz2=qh2d(i,k) !qhwrf(i,k)
      qgz2=qg2d(i,k) !qgwrf(i,k)
      if (k .lt. kmax-2 .and. tairc(i) .ge. -5) then
          qhz2=qh2d(i,k+1) !qhwrf(i,k+1)
          qgz2=qg2d(i,k+1) !qgwrf(i,k+1)
      endif

      !
      ! PSD mapping
      !
      call sgmap(1,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftns0(i))
      call sgmap(2,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftng0(i))
      call sgmap(3,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),fros0(i))
      call sgmap(4,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftnh0(i))

      if (tair(i).ge.t0) then
          tairc(i)=tair(i)-t0
          ftns(i)=1.
          ftng(i)=1.
          ftnh(i)=1.
          ftns(i)=ftns0(i)
          ftng(i)=ftng0(i)
          ftnh(i)=ftnh0(i)
          dd(i)=r11t*tairc(i)*(r101r(i)/zs(i)**2+r102rf(i)/zs(i)**bsh5)*ftns(i)
          psmlt(i)=r2is*min(qs(i),max(dd(i),0.0))
          if (r00(i)*qg(i).gt.qrog2) then
             y2(i)=(r191r(i)/zg(i)**2+r192rf2(i)/zg(i)**bgh5_2)*ftng(i)
          else
             y2(i)=(r191r(i)/zg(i)**2+r192rf(i)/zg(i)**bgh5)*ftng(i)
          endif 
          dd1(i)=tairc(i)*(r19t*y2(i)+r19at*(qgacw(i) +qgacr(i)))
          pgmlt(i)=r2ig*min(qg(i),max(dd1(i),0.0))
       
          y1(i)=tca(i)*tairc(i)-alvr(i)*dwv(i)*(rp0(i)-(qv(i)+qb0))
          y3(i)=.78/zh(i)**2+h19aq(i)*scv(i)/zh(i)**bhh5
          dda0(i)=h19rt(i)*y1(i)*y3(i)*ftnh(i)+r19at*tairc(i)*(qhacw(i)+qhacr(i))
          phmlt(i)=r2ih*max(0.0, min(dda0(i), qh(i)))

          pt(i)=pt(i)-afcp(i)*(psmlt(i)+pgmlt(i)+phmlt(i))
          tair(i)=(pt(i)+tb0)*pi0(i)
          qr(i)=qr(i)+psmlt(i)+pgmlt(i)+phmlt(i)
          qs(i)=qs(i)-psmlt(i)
          qg(i)=qg(i)-pgmlt(i)
          qh(i)=qh(i)-phmlt(i)

      endif !tair  ! processes 11 & 19

! 
!* 24 * pihom : homogeneous freezing of qc to qi (t < t00)        **24**
!* 25 * pidw : deposition growth of qc to qi ( t0 < t <= t00)     **25**
!* 26 * pimlt : melting of qi to qc (t >= t0)                     **26**
!****** pimm  : immersion freezing of qc to qi (t < t0)           ******
!****** pcfr  : contact nucleation of qc to qi (t < t0)           ******

      if (qc(i).le.cmin1) qc(i)=0.0
      if (qi(i).le.cmin1) qi(i)=0.0

      ftns(i)=1.
      ftng(i)=1.
      ftns0(i)=1.
      ftng0(i)=1.

      qhz2=qh2d(i,k) 
      qgz2=qg2d(i,k) 
      if (k .lt. kmax-2 .and. tairc(i) .ge. -5) then
         qhz2=qh2d(i,k+1) 
         qgz2=qg2d(i,k+1) 
      endif
      !
      ! PSD mapping
      !
      call sgmap(1,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftns0(i))
      call sgmap(2,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftng0(i))

      if (tair(i).le.t00) then
          pihom(i)=qc(i)
      else
          pihom(i)=0.0
      endif 
      if (tair(i).ge.t0) then
          pimlt(i)=qi(i)
      else
          pimlt(i)=0.0
      endif
      pidw(i)=0.0
           
      if (tair(i) .lt. t0 .and. tair(i) .gt. t00) then
          tairc(i)=tair(i)-t0
          y1(i) = max( min(tairc(i), -1.), -31.)
          it(i) = int(abs(y1(i)))
          y2(i)=aa1(it(i))
          y3(i) = aa2(it(i))
          if (tairc(i).le.-5.)then                       !  meyers
              rtair(i)=1./(tair(i)-c76)
              y5(i)=exp(c218-c580*rtair(i))
              qsi(i)=rp0(i)*y5(i)
              ssi(i)=(qv(i)+qb0)/qsi(i)-1.
              fssi=min(xssi,max(rssi,xssi*(tairc(i)+44.)/(44.0-38.0))) !max ssi f(tair)
              fssi=rssi
              wssi=w2d(i,k)/100.-2.0
              if(wssi.gt.0.) fssi=rssi+min(xssi,wssi*0.01)
              fssi=min(ssi(i),fssi)

              r_nci=max(1.e-3*exp(-.639+12.96*fssi), 0.528e-3)      ! meyers et al.
                                                                      ! roger found the bug on 2012/04/27
              if (r_nci.gt.15.) r_nci=15.                          !cap at 15000/liter

              ! cooper curve
              if( tairc(i) .lt. -40.0 ) then
                c_nci = 5.0e-6*exp(0.304*40.0)
              else
                c_nci = 5.0e-6*exp(0.304*abs(tairc(i)))
              end if
              r_nci = c_nci
              !  cap r_nci to the amount corresponding to crystal size of ami50
              r_nci = max(r_nci,r00(i)*qi(i)/ami50)

              ! DeMott's ice nucleation (coupling with GOCART only)
#if ( WRF_CHEM == 1)
              if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
                    chem_opt == 302 .or. chem_opt == 303) .and. &
                    (gsfcgce_gocart_coupling == 1) ) then
                   ! conversion rate of cloud water to ice in the bergeron porcess based on 
                   ! meyer + demott formulae
                   !
                   ! convert gocart aerosol mass conc to in
                   ! p0 need to be converted from g*cm/s2/cm2 to mb (hpa)
                   ! call mass2icn(p0(i,k)*0.001,tair(i),aero(i,k,:), icn_out)
                   ! icn_cgs(i,k) = max(min_icn, icn_out) * 1.e-3   !in conc [#/cm3] <-- [#/litter]
                   ! icn_cgs(i,k) = icn_out * 1.e-3 !in conc [#/cm3]
                   ! r_nci = icn_cgs(i,k)  !demotto's formuale
                   r_nci = icn_out(i) * 1.e-3  !demotto's formuale
              end if
#endif

              dd(i)=(r00(i)*qi(i)/r_nci)**y3(i)                  !meyers
              pidw(i)=min(rr0(i)*d2t*y2(i)*r_nci*dd(i), qc(i)) !meyers

              ! water vapor diffusivity correction term
              esi(i)  = qsi(i)/rp0(i)*c610
              y4(i)   = 1./tair(i)
!             term1     = y4(i)*(rn10a*y4(i)-rn10b)
              term1     = y4(i)*(rn20a*y4(i)-rn20b)
              term2     = rn10c*tair(i)/esi(i)
              fdwv      = (term1+term2)/(term1+term2*dwv0/dwv(i))
              pidw(i)=min(rr0(i)*d2t*y2(i)*r_nci*dd(i)*fdwv,qc(i)) !meyers
          endif  !tairc

          pimm(i)=0.0
          pcfr(i)=0.0

          if (qc(i) .gt. 0.0) then
              y4(i) = 1./(tair(i)-c358)
              qsw(i)=rp0(i)*exp(c172-c409*y4(i))
              xncld=qc(i)/4.e-9                         !cloud number
              esat=0.6112*exp(17.67*tairc(i)/(tairc(i)+243.5))*10.
              rv=0.622*esat/(p2d(i,k)/1000.-esat)
              rlapse_m=980.616*(1.+2.5e6*rv/287./tair(i))/          &     
                       (1004.67+2.5e6*2.5e6*rv*0.622/(287.*tair(i)*tair(i)))
              delt=rlapse_m*w2d(i,k)                     !roger
              if (delt.lt.0.) delt=0.
              bhi=1.01e-2 
              pimm(i)=xncld*bhi*4.e-9*exp(-tairc(i))*delt*d2t*4.e-9
 
              xccld=xncld*r00(i)                     ! cloud number concentration
              xknud=7.37*tair(i)/(288.*ra*p2d(i,k))   ! knudsen number
              alpha=1.257+0.400*exp(-1.10/xknud)     ! cunningham correction (p&klett)
 
              cunnf=1.+alpha*xknud                   ! cunningham correction (p&klett)

              if (tairc(i).ge.0.) then
                  dvair=(1.718+0.0049*tairc(i))*1.e-4   !dynamic visc air (prupp&klett)
              else 
                  dvair=(1.718+0.0049*tairc(i)-1.2e-5*tairc(i)**2)*1.e-4  
              endif !tairc
              diffar=1.3804e-16*tair(i)/6./cpi/dvair/ra*cunnf     !aerosol diffusion via p&klett
              if (qv(i)+qb0-qsw(i).lt.0.) then                       !only when cloud evaporating
                pcfr(i)=4.e-9*4.*cpi*rc*diffar*xccld*cna*rr0(i)*d2t    !brownian part only via cotton
              endif
          else  !qc
              tairc(i)=tair(i)-t0
              y1(i)=max( min(tairc(i), -1.), -31.)
              it(i)=int(abs(y1(i)))
              y2(i)=aa1(it(i))
              y3(i)=aa2(it(i))
              y4(i)=exp(abs(.5*tairc(i)))
              dd(i)=(r00(i)*qi(i)/(r25a*y4(i)))**y3(i)
              pidw(i)=min(r25rt(i)*y2(i)*y4(i)*dd(i),qc(i))
          endif  !qc

      endif  !tair

   enddo


!dir$ vector aligned
   do i=1,irestrict

!      steve: please check
      y1(i)=pihom(i)+pidw(i)+pimm(i)+pcfr(i)-pimlt(i)

      if (y1(i).gt.qc(i)) then
          y1(i)=qc(i)
          y2(i)=1.
          y3(i)=pihom(i)+pidw(i)+pimm(i)+pcfr(i)
          if(y3(i).ne.0.) y2(i)=(qc(i)+pimlt(i))/y3(i)
          pihom(i)=pihom(i)*y2(i)
          pidw(i)=pidw(i)*y2(i)
          pimm(i)=pimm(i)*y2(i)
          pcfr(i)=pcfr(i)*y2(i)
      endif  !y1

      pt(i)=pt(i)+afcp(i)*y1(i)
      tair(i)=(pt(i)+tb0)*pi0(i)
      qc(i)=qc(i)-y1(i)
      qi(i)=qi(i)+y1(i)

!* 31 * pint  : initiation of qi                                  **31**
!* 32 * pidep : deposition of qi                                  **32**
!
!     calculation of pint uses different values of the intercept and slope for
!     the fletcher equation. also, only initiate more ice if the new number
!     concentration exceeds that already present.

      tair(i)=(pt(i)+tb0)*pi0(i)
      if (tair(i) .lt. t0) then
          if (qi(i) .le. cmin) qi(i)=0.
          tairc(i)=tair(i)-t0
          rtair(i)=1./(tair(i)-c76)
          y2(i)=exp(c218-c580*rtair(i))
          qsi(i)=rp0(i)*y2(i)
          esi(i)=c610*y2(i)
          ssi(i)=(qv(i)+qb0)/qsi(i)-1.
          y1(i)=1./tair(i)
          y3(i)=sqrt(qi(i))
          ! water vapor diffusivity correction term
          term1   = y1(i)*(rn20a*y1(i)-rn20b)
          term2   = rn10c*tair(i)/esi(i)
          dd(i) = term1+term2
          fdwv    = dd(i)/(term1+term2*dwv0/dwv(i))

          dm(i)=max(qv(i)+qb0-qsi(i),0.0)
          rsub1(i)=cs580(i)*qsi(i)*rtair(i)*rtair(i)
          dep(i)=dm(i)/(1.+rsub1(i))
          if (tairc(i).le.-5.) then     
             y4(i)=1./(tair(i)-c358)
             qsw(i)=rp0(i)*exp(c172-c409*y4(i))
             fssi=min(xssi,max(rssi,xssi*(tairc(i)+44.)/(44.-38.)))
             fssi=rssi
             wssi=w2d(i,k)/100.-2.0
             if(wssi.gt.0.) fssi=rssi+min(xssi,wssi*0.01)
             fssi=min(ssi(i),fssi)
!            r_nci=min(1.e-3*exp(-.639+12.96*fssi),1.) 
             r_nci=max(1.e-3*exp(-.639+12.96*fssi),0.528e-3)
             if (r_nci.gt.15.) r_nci=15.   

             ! cooper curve
             if( tairc(i) .lt. -40.0 ) then
                c_nci = 5.0e-6*exp(0.304*40.0)
             else
                c_nci = 5.0e-6*exp(0.304*abs(tairc(i)))
             end if
             r_nci = c_nci
             ! cap r_nci to the amount corresponding to crystal size of ami50
             r_nci = max(r_nci,r00(i)*qi(i)/ami50)

#if ( WRF_CHEM == 1)
                ! emk...only execute if gocart and coupling is on
             if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
                   chem_opt == 302 .or. chem_opt == 303) .and. &
                   (gsfcgce_gocart_coupling == 1) ) then
                   !jjs 20110602 vvvvv
                   ! conversion rate of cloud water to ice in the bergeron porcess based on meyer + demott formulae
                   !
                   ! convert gocart aerosol mass conc to in
                   !      p0 need to be converted from g*cm/s2/cm2 to mb (hpa)
                   !                call mass2icn(p0(i,k)*0.001,tair(i),aero(i,k,:), icn_out)
                   !                icn_cgs(i,k) = max(min_icn, icn_out) * 1.e-3   !in conc [#/cm3] <-- [#/litter]
                   !                call mass2icn(p0(i,k)*0.001,tair(i),aero(i,k,:), icn_out,i,k)
                   !                icn_out = min(1.e3, max(0.01e0 ,  icn_out) )
                   !                icn_cgs(i,k) = icn_out * 1.e-3 !in conc [#/cm3]
                   !                r_nci = icn_cgs(i,k)  !demotto's formuale
                   r_nci = icn_out(i) * 1.e-3  !demotto's formuale
                   !jjs 20110602 ^^^^^
             end if
#endif

             pidep(i)=max(r32rt(i)*1.e4*fssi*sqrt(r_nci)*y3(i)/     & !meyers
             dd(i)*fdwv, -qi(i))                    !fix sel
             if(qi(i).le.cmin) pidep(i)=0.
             dd(i)=max(1.e-9*r_nci/r00(i)-qi2d(i,k)*1.e-9/ami50, 0.) 
             pint(i)=max(min(dd(i),dm(i)),0.)
             pint(i)=min(pint(i)+pidep(i), dep(i))
             pt(i)=pt(i)+ascp(i)*pint(i)
             tair(i)=(pt(i)+tb0)*pi0(i)
             qv(i)=qv(i)-pint(i)
             qi(i)=qi(i)+pint(i)
         endif  !taric
      endif  !tair                  

! wrf satice has new_ice_sat option 0, 1 and 2
! steve's satice has new_ice_sat option 0, 1, 2, 3 and 9
! option 0, 1 and 2 are identical in both satice
! i added option 3 below and wrapped them with "if (improve.eq.3)"
!*****   tao et al (1989) saturation technique  ***********************
!!!!
!!!   new_ice_sat option 9 from  steve's satice
!!!   sequential, non-iterative, ssi
!!!!

      dep(i)=0.0
      cnd(i)=0.0
      tair(i)=(pt(i)+tb0)*pi0(i)
      if (tair(i).ge.t00) then
          y1(i)=1./(tair(i)-c358)
          qsw(i)=rp0(i)*exp(c172-c409*y1(i))
          dd(i)=cp409(i)*y1(i)*y1(i)
          dm(i)=qv(i)+qb0-qsw(i)
          cnd(i)=dm(i)/(1.+avcp(i)*dd(i)*qsw(i))
          !  condensation or evaporation of qc  
          cnd(i)=max(-qc(i),cnd(i))

          ! thresh_evap is a function of dx and defined in the beginning of the subroutine
          if(w2d(i,k) .ge. thresh_evap)  cnd(i)=max(0.,cnd(i)) !reduce spurious evap

          pt(i)=pt(i)+avcp(i)*cnd(i)
          tair(i)=(pt(i)+tb0)*pi0(i)
          qv(i)=qv(i)-cnd(i)
          qc(i)=qc(i)+cnd(i)
      endif

      if (tair(i).le.273.16) then
          !  deposition or sublimation of qi  
          y1(i)=1./(tair(i)-c358)
          qsw(i)=rp0(i)*exp(c172-c409*y1(i))
          y2(i)=1./(tair(i)-c76)
          qsi(i)=rp0(i)*exp(c218-c580*y2(i))
          fssi=min(xssi,max(rssi,xssi*(tair(i)-t0+44.0)/(44.0-38.0)))
          fssi=rssi
          wssi=w2d(i,k)/100.-2.0
          if(wssi.gt.0.) fssi=rssi+min(xssi,wssi*0.01)
          y3(i)=1.+min((qsw(i)-qsi(i))/qsi(i), fssi)
          y4(i)=qsi(i)*y3(i)
          if (tair(i).le.268.16.and.(qv(i)+qb0.gt.y4(i))) then
              dd1(i)=cp580(i)*y2(i)*y2(i)
              dep(i)=(qv(i)+qb0-y4(i))/(1.+ascp(i)*dd1(i)*y4(i))
          else if (qv(i)+qb0.lt.xsubi*qsi(i).and.qi(i).gt.cmin) then
              dd1(i)=cp580(i)*y2(i)*y2(i)
              dep(i)=(qv(i)+qb0-xsubi*qsi(i))/(1.+ascp(i)*dd1(i)*xsubi*qsi(i))
              dep(i)=max(-qi(i),dep(i))
          endif
          if(w2d(i,k).ge.0.) dep(i)=max(0.,dep(i))  !reduce spurious sublimation
            pt(i)=pt(i)+ascp(i)*dep(i)
            tair(i)=(pt(i)+tb0)*pi0(i)
            qv(i)=qv(i)-dep(i)
            qi(i)=qi(i)+dep(i)
      endif

!* 10 * psdep : deposition or sublimation of qs                   **10**
!* 20 * pgsub : sublimation of qg                                 **20**

      psdep(i)=0.0
      pgdep(i)=0.0
      phdep(i)=0.0
      pssub(i)=0.0
      pgsub(i)=0.0
      phsub(i)=0.0
      pvapg(i)=0.0
      pvaph(i)=0.0

      if (qc(i)+qi(i).gt.1.e-5) then
         dlt1(i)=1.
      else    
         dlt1(i)=0.
      endif

      if (tair(i) .lt. t0) then 
          rtair(i)=1./(tair(i)-c76)
          y2(i)=exp(c218-c580*rtair(i))
          qsi(i)=rp0(i)*y2(i)
          esi(i)=c610*y2(i)
          ftns(i)=1.
          ftng(i)=1. 
          ftnh(i)=1.             

          ssi(i)=(qv(i)+qb0)/qsi(i)-1.
          if(ssi(i).gt.0.) dlt1(i)=1.
          if(ssi(i).le.0.) dlt1(i)=0.
          dm(i)=qv(i)+qb0-qsi(i)
          rsub1(i)=cs580(i)*qsi(i)*rtair(i)*rtair(i)
          dd1(i)=dm(i)/(1.+rsub1(i))
          y3(i)=1./tair(i)
          ! water vapor diffusivity correction
          term1   = y3(i)*(rn20a*y3(i)-rn20b)
          term2   = rn10c*tair(i)/esi(i)
          dd(i) = term1+term2
          fdwv    = dd(i)/(term1+term2*dwv0/dwv(i))
          tairc(i)=tair(i)-t0

          ftns(i)=1.
          ftng(i)=1.
          ftnh(i)=1.
          ftns0(i)=1.
          ftng0(i)=1.
          qhz2=qh2d(i,k) !qhwrf(i,k)
          qgz2=qg2d(i,k) !qgwrf(i,k)
          if (k .lt. kmax-2 .and. tairc(i) .ge. -5) then
             qhz2=qh2d(i,k+1) !qhwrf(i,k+1)
             qgz2=qg2d(i,k+1) !qgwrf(i,k+1)
          endif

          !
          ! PSD mapping
          !
          call sgmap(1,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftns0(i))
          call sgmap(2,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftng0(i))
          call sgmap(3,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),fros0(i))
          call sgmap(4,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftnh0(i))

          ftns(i)=ftns0(i)
          ftng(i)=ftng0(i)
          ftnh(i)=ftnh0(i)

          y4(i)=r10t*ssi(i)*(r101r(i)/zs(i)**2+r102rf(i)/zs(i)**bsh5) /dd(i)*ftns(i)*fdwv
          psdep(i)=r2is*max(-qs(i), y4(i))
          if(qs(i).le.cmin) psdep(i)=0.
          dd(i)=y3(i)*(rn20a*y3(i)-rn20b)+rn10c*tair(i)/esi(i)
          y2(i)=r191r(i)/zg(i)**2+r192rf(i)/zg(i)**bgh5
          if(r00(i)*qg(i).gt.qrog2) &
              y2(i)=(r191r(i)/zg(i)**2+r192rf2(i)/zg(i)**bgh5_2)
          pgdep(i)=r2ig*max(-qg(i), r20t*ssi(i)*y2(i)/dd(i)*ftng(i)*fdwv)
          if(qg(i).le.cmin) pgdep(i)=0.
          y1(i)=h10ar(i)/(tca(i)*tair(i)**2)+1./(dwv(i)*qsi(i))
          y2(i)=.78/zh(i)**2+h20bq(i)*scv(i)/zh(i)**bhh5
          phdep(i)=r2ih*max(-qh(i), h20t(i)*ftnh(i)*ssi(i)*y2(i)/y1(i))
          if(qh(i).le.cmin) phdep(i)=0.

          pssub(i)=min(psdep(i), 0.)
          psdep(i)=max(psdep(i), 0.)
          pgsub(i)=min(pgdep(i), 0.)
          pgdep(i)=max(pgdep(i), 0.)
          phsub(i)=min(phdep(i), 0.)
          phdep(i)=max(phdep(i), 0.)

          y5(i)=min(0.,dd1(i))
          dd1(i)=max(0.,dd1(i))

          if(dlt1(i).eq.1.)then                              !di
             y1(i)=psdep(i)+pgdep(i)+phdep(i)
             if(y1(i).ge.dd1(i))then
                psdep(i)=psdep(i)/y1(i)*dd1(i)            !...
                pgdep(i)=pgdep(i)/y1(i)*dd1(i)
                phdep(i)=phdep(i)/y1(i)*dd1(i)
             endif
          endif                                                !di

          if(qc(i).le.1.e-5) then
             vap_frac=2.0*min((tairc(i)/(t00-t0))**2,1.0)
             pvapg(i)=vap_frac*pgdep(i)
             pvaph(i)=vap_frac*phdep(i)
          endif
          if(pgdep(i).gt.0.) pvapg(i)=min(pvapg(i),qg(i)+pgdep(i))
          if(phdep(i).gt.0.) pvaph(i)=min(pvaph(i),qh(i)+phdep(i))

          if(dlt1(i).eq.0.)then
             y1(i)=max(pssub(i)+pgsub(i)+phsub(i), y5(i))
             if(y5(i).gt.(pssub(i)+pgsub(i)+phsub(i)))then
             y3(i)=(pssub(i)+pgsub(i)+phsub(i))
             if(y3(i).ne.0.0)then
                pssub(i)=pssub(i)/y3(i)*y5(i)
                pgsub(i)=pgsub(i)/y3(i)*y5(i)
                phsub(i)=phsub(i)/y3(i)*y5(i)
             endif
          endif
         endif
         pssub(i)=-pssub(i)
         pgsub(i)=-pgsub(i)
         phsub(i)=-phsub(i)
         y1(i)=psdep(i)+pgdep(i)+phdep(i) -pssub(i)-pgsub(i)-phsub(i)

         pt(i)=pt(i)+ascp(i)*y1(i)
         tair(i)=(pt(i)+tb0)*pi0(i)
         qv(i)=qv(i)-y1(i)
         qs(i)=qs(i)+psdep(i)-pssub(i)+pvapg(i)+pvaph(i)
         qg(i)=qg(i)+pgdep(i)-pgsub(i)-pvapg(i)
         qh(i)=qh(i)+phdep(i)-phsub(i)-pvaph(i)
      endif   ! if (tair(i) .lt. t0)

!
!* 23 * ern : evaporation of qr (subsaturation)                   **23**
! steve did not make any improvement on this process
!
      ern(i)=0.0
      if (qr(i) .gt. 0.0) then
          tair(i)=(pt(i)+tb0)*pi0(i)
          rtair(i)=1./(tair(i)-c358)
          y2(i)=exp( c172-c409*rtair(i) )
          esw(i)=c610*y2(i)
          qsw(i)=rp0(i)*y2(i)
          ssw(i)=(qv(i)+qb0)/qsw(i)-1.
          dm(i)=qv(i)+qb0-qsw(i)
          rsub1(i)=cv409(i)*qsw(i)*rtair(i)*rtair(i)
          dd1(i)=max(-dm(i)/(1.+rsub1(i)),0.0)
          y3(i)=1./tair(i)
          ! water vapor diffusivity correction term
          term1   = y3(i)*(rn30a*y3(i)-rn10b)
          term2   = rn10c*tair(i)/esw(i)
          dd(i) = term1+term2
          fdwv    = dd(i)/(term1+term2*dwv0/dwv(i))
          ftnw=1.
          if (qr(i) .gt. cmin.and.tair(i).gt.t0) then   !no need to check qc, no qc if ssw < 0
              bin_factor(i)=0.11*(1000.*qr(i))**-1.27 + 0.98
              bin_factor(i)=min(bin_factor(i),1.30)
              ftnw=1./bin_factor(i)**3.35
              ftnwmin=r00(i)*qr(i)/draimax
          if(qr(i).le.0.001) ftnw=max(ftnw,ftnwmin/tnw)
             y4(i)=r00(i)*qr(i)
             y1(i)=sqrt(y4(i))
             y2(i)=sqrt(y1(i))
             zr(i)=zrc/y2(i)*ftnw**0.25
          endif !qr
          y1(i)=-r23t*ssw(i)*(r231r(i)/zr(i)**2+r232rf(i)/zr(i)**3)/dd(i)*ftnw*fdwv
          ern(i)=min(dd1(i),qr(i),max(y1(i),0.0))
          pt(i)=pt(i)-avcp(i)*ern(i)
          tair(i)=(pt(i)+tb0)*pi0(i)
          tairc(i)=tair(i)-t0
          qv(i)=qv(i)+ern(i)
          qr(i)=qr(i)-ern(i)
      endif !qr

!  add processes 30 & 33 fpr pmltg and pmlts
!
!* 30 * pmltg : evaporation of melting qg                         **30**
!* 33 * pmlts : evaporation of melting qs                         **33**

      pmlts(i)=0.0
      pmltg(i)=0.0

      tair(i)=(pt(i)+tb0)*pi0(i)
      tairc(i)=tair(i)-t0

      ftns0(i)=1.
      ftng0(i)=1.

      qhz2=qh2d(i,k) 
      qgz2=qg2d(i,k) 
      if (k .lt. kmax-2 .and. tairc(i) .ge. -5) then
          qhz2=qh2d(i,k+1) 
          qgz2=qg2d(i,k+1) 
      endif
      !
      ! PSD mapping
      !
      call sgmap(1,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftns0(i))
      call sgmap(2,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftng0(i))
      call sgmap(3,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),fros0(i))

      if (tair(i) .ge. t0) then
          ftns(i)=1.
          ftng(i)=1.
          ftns(i)=ftns0(i)
          ftng(i)=ftng0(i)
          rtair(i)=1./(t0-c358)
          y2(i)=exp( c172-c409*rtair(i) )
          esw(i)=c610*y2(i)
          qsw(i)=rp0(i)*y2(i)
          ssw(i)=1.-(qv(i)+qb0)/qsw(i)
          dm(i)=qsw(i)-qv(i)-qb0
          rsub1(i)=cv409(i)*qsw(i)*rtair(i)*rtair(i)
          dd1(i)=max(dm(i)/(1.+rsub1(i)),0.0)
          y3(i)=1./tair(i)
          ! water vapor diffusivity correction term
          term1   = y3(i)*(rn30a*y3(i)-rn10b)
          term2   = rn10c*tair(i)/esw(i)
          dd(i) = term1+term2
          fdwv    = dd(i)/(term1+term2*dwv0/dwv(i))
          y1(i)=ftng(i)*r30t*ssw(i)*(r191r(i)/zg(i)**2+r192rf(i)/zg(i)**bgh5)/dd(i)*fdwv
          if(r00(i)*qg(i).gt.qrog2)                   &
              y1(i)=ftng(i)*r30t*ssw(i)*(r191r(i)/zg(i)**2+r192rf2(i) &
                     /zg(i)**bgh5_2)/dd(i)*fdwv
          pmltg(i)=r2ig*min(qg(i),max(y1(i),0.0))
          y1(i)=ftns(i)*r33t*ssw(i)*(r331r(i)/zs(i)**2+r332rf(i)/zs(i)**bsh5)/dd(i)*fdwv
          pmlts(i)=r2is*min(qs(i),max(y1(i),0.0))
          y1(i)=min(pmltg(i)+pmlts(i),dd1(i))
          pmltg(i)=y1(i)-pmlts(i)
          pt(i)=pt(i)-ascp(i)*y1(i)
          tair(i)=(pt(i)+tb0)*pi0(i)
          qv(i)=qv(i)+y1(i)
          qs(i)=qs(i)-pmlts(i)
          qg(i)=qg(i)-pmltg(i)
      endif  !tair t0


  
      if (qc(i) .le. cmin) qc(i)=0.
      if (qr(i) .le. cmin) qr(i)=0.
      if (qi(i) .le. cmin) qi(i)=0.
      if (qs(i) .le. cmin) qs(i)=0.
      if (qg(i) .le. cmin) qg(i)=0.
      if (qh(i) .le. cmin) qh(i)=0.
      pt2d(i,k)=pt(i)
      qv2d(i,k)=qv(i)
      qc2d(i,k)=qc(i)
      qr2d(i,k)=qr(i)
      qi2d(i,k)=qi(i)
      qs2d(i,k)=qs(i)
      qg2d(i,k)=qg(i)
      qh2d(i,k)=qh(i)


      dd(i)=max(-cnd(i), 0.)
      cnd(i)=max(cnd(i), 0.)
      dd1(i)=max(-dep(i), 0.)  !bug fix by di
      dep(i)=max(dep(i), 0.)
      del=0.
      if(whacr(i) .lt. 0.) del=1.

      sccc=cnd(i)
      seee=dd(i)+ern(i)
      sddd=dep(i)+amax1(pint(i),0.0)+psdep(i)+pgdep(i)+phdep(i)
      ssss=dd1(i)-amin1(pint(i),0.0)+pssub(i)+pgsub(i)+phsub(i)+pmlts(i)+pmltg(i)
      smmm=psmlt(i)+pgmlt(i)+pimlt(i)+qracs(i)+phmlt(i)+qracg(i)-del*whacr(i)
      sfff=psacw(i)*d2t+piacr(i)*d2t+psfw(i)*d2t+pgfr(i)*d2t   &
           +dgacw(i)*d2t+dgacr(i)*d2t+psacr(i)*d2t+pihom(i) &
           +pidw(i)+pimm(i)+pcfr(i)+pihms(i)*d2t    &
           +pihmg(i)*d2t+phfr(i)*d2t+dhacw(i)*d2t+dhacr(i)*d2t+pihmh(i)*d2t

      !
      ! snapsot diabatic heating rate (deg k / s)
      !
      physc2d(i,k) = avc * sccc / d2t       !k/s
      physe2d(i,k) = avc * seee / d2t       !k/s
      physd2d(i,k) = asc * sddd / d2t       !k/s
      physs2d(i,k) = asc * ssss / d2t       !k/s
      physf2d(i,k) = afc * sfff / d2t       !k/s
      physm2d(i,k) = afc * smmm / d2t       !k/s

      !
      ! accmulated diabatic heating (unit: deg k, in potential temperature)
      !
      acphysc2d(i,k) = acphysc2d(i,k) + avcp(i) * sccc 
      acphyse2d(i,k) = acphyse2d(i,k) + avcp(i) * seee 
      acphysd2d(i,k) = acphysd2d(i,k) + ascp(i) * sddd 
      acphyss2d(i,k) = acphyss2d(i,k) + ascp(i) * ssss 
      acphysf2d(i,k) = acphysf2d(i,k) + afcp(i) * sfff 
      acphysm2d(i,k) = acphysm2d(i,k) + afcp(i) * smmm 

      ! radar reflectivity calculation

      !
      ! PSD mapping
      !
      call sgmap(1,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftns0(i))   !di 20160114
      call sgmap(2,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftng0(i))   !di 20160114
      call sgmap(3,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),fros0(i))   !di 20160114
      call sgmap(4,qs(i),qg(i),qgz2,qh(i),qhz2,r00(i),tairc(i),ftnh0(i))   !di 20160114

      a_1=1.e6*r00(i)*qr(i)
      a_2=1.e6*r00(i)*qs(i)
      a_3=1.e6*r00(i)*qg(i)
      a_4=1.e6*r00(i)*qh(i) 

      ucor=3071.29/tnw**0.75
      ucos=687.97*roqs**0.25/tns**0.75
      ucog=687.97*roqg**0.25/tng**0.75
      ucog2=687.97*roqg2**.25/tng**.75
      ucoh=687.97*roqh**0.25/tnh**0.75
      uwet=4.464**0.95

      ftnw=1.
      if(qr(i).gt.cmin .and. qc(i).lt.cmin)then
         bin_factor(i)=0.11*(1000.*qr(i))**-1.27 + 0.98
!        bin_factor(i)=min(bin_factor(i),1.35)
         bin_factor(i)=min(bin_factor(i),1.30)
         ftnw=1./bin_factor(i)**3.50
         ftnwmin=r00(i)*qr(i)/draimax
         if(qr(i).le.0.001) ftnw=max(ftnw,ftnwmin/tnw)
      endif
      a_11=ucor*(max(0.,a_1))**1.75/ftnw**0.75

      ftns(i)=1.
      ftng(i)=1.
      ftns(i)=ftns0(i)**0.75
      ftng(i)=ftng0(i)**0.75
      ftnh(i)=ftnh0(i)**0.75  ! di 20160114
      fros(i)=1.
      fros(i)=fros0(i)**.25
      a_22=ucos*(max(0.,a_2))**1.75/ftns(i)*fros(i)
      a_33=ucog*(max(0.,a_3))**1.75/ftng(i)
      if(a_3.ge.qrog2) a_33=ucog2*(max(0.,a_3))**1.75/ftng(i)
      a_44=ucoh*(max(0.,a_4))**1.75/ftnh(i) ! di 20160114

      if (tair(i).lt.273.16) then
          zdry = max(1.e-9,a_11+a_22+a_33+a_44) !rain,snow,graupel,hail,cloud ice,cloud water
          dbz(i,k) = 10.*alog10(zdry)
      else         
          zwet0 = a_11+uwet*(a_22+a_33+a_44)**.95
          zwet = max(1.e-9,zwet0)
          dbz(i,k) = 10.*alog10(zwet)
      endif

   enddo

!dir$ vector aligned
   do i=1,irestrict
      
#if ( WRF_CHEM == 1)
      ! emk...nuclei only calculated when coupling
      if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
            chem_opt == 302 .or. chem_opt == 303) .and. &
            (gsfcgce_gocart_coupling == 1) ) then
         icn_diag(i,k) = icn_out(i)  ! #/litre
         nc_diag(i,k) = ccn_out(i)  ! #/cm3
      else
         icn_diag(i,k) = 0
         nc_diag(i,k) = 0
      end if
#endif      

      !
      ! calculate effective radius for all cloud species
      !   eff_rad is a function of the slope parameter (lambda)

      ! rain
      if (qr2d(i,k) .lt. cmin) then
         refr2d(i,k) = 0.e0
      else
         refr2d(i,k) = eff_rad(zr(i))
      endif
      ! snow
      if (qs2d(i,k) .lt. cmin) then
         refs2d(i,k) = 0.e0
      else
         refs2d(i,k) = eff_rad(zs(i))
      endif
      ! graupel
      if (qg2d(i,k) .lt. cmin) then
         refg2d(i,k) = 0.e0
      else
         refg2d(i,k) = eff_rad(zg(i))
      endif
      ! hail
      if (qh2d(i,k) .lt. cmin) then
         refh2d(i,k) = 0.e0
      else
         refh2d(i,k) = eff_rad(zh(i))
      endif

      ! 
      ! cloud ice (temperature dependent, even if you are running WRF_CHEM)
      !
      if (qi2d(i,k) .lt. cmin) then
          refi2d(i,k) = 0.e0
      else
         ! for cloud ice effective radius depends on temperature profile, formula from gce
         refi2d(i,k) = 125.e0 +(tair(i)-243.16)*5.e0     ! [micron]
         if (tair(i) .gt. 243.16) refi2d(i,k) = 125.e0
         if (tair(i) .lt. 223.16) refi2d(i,k) = 25.e0
      endif

      !
      ! cloud water
      !
      if (qc2d(i,k) .lt. cmin) then
          refc2d(i,k) = 0.e0
      else
          l_cloud = qc2d(i,k) * rho2d(i,k)             ! cloud water [g/cm3]
#if ( WRF_CHEM == 1)
   ! when running with WRF_CHEM and using aerosol coupling in goddard mp
          ! cpi: const_pi = 4.*atan(1.)         ~ 3.1415
          ! roqr: 1.0 g/cm**3, liquid water density
          ! roqi: 0.9179, ice density
          if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
                chem_opt == 302 .or. chem_opt == 303) .and. &
                (gsfcgce_gocart_coupling == 1) ) then
               ! for cloud water, estimate lambda (slope of gamma distribution)
                 mu = min(15.e0, (1000.e0/ccn_out(i) + 2.e0))
                 gamfac3 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+3.e0) )
                 gamfac1 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+1.e0) )
                 lambda = (4.e0/3.e0*cpi*roqr*ccn_out(i)/l_cloud*gamfac1)**(1.e0/3.e0)  ! [1/cm]
                 refc2d(i,k) = 1.e0/lambda * gamfac3 * 1.e4  !effective radius [micron]
          else
   ! when running with WRF_CHEM but no aerosol coupling in goddard mp
              if (xland(i) .eq. 1.0) then
                 ccn_ref = ccn_over_land
              else if (xland(i) .eq. 2.0) then
                 ccn_ref = ccn_over_water
              else
                 print *,' xland is not 1. or 2., run stopped'
                 ! emk nuwrf
                 call wrf_error_fatal(' xland is not 1. or 2., run stopped')
!                 stop
              endif
                ! for cloud water, estimate lambda (slope of gamma distribution)
                     mu = min(15.e0, (1000.e0/ccn_ref + 2.e0))
                     gamfac3 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+3.e0) )
                     gamfac1 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+1.e0) )
                     lambda = (4.e0/3.e0*cpi*roqr*ccn_ref/l_cloud*   &
                              gamfac1)**(1.e0/3.e0)  ! [1/cm]
                     refc2d(i,k) = 1.e0/lambda * gamfac3 * 1.e4  !effective radius [micron]
          endif ! chem_opt and gsfcgce_gocart_coupling
#else
          !
          ! not running with WRF_CHEM
          ! ccn_over_land = 1500  ! [#/cm3] climatological value
          ! ccn_over_water = 150  ! [#/cm3] climatological value
          if (xland(i) .eq. 1.0) then
              ccn_ref = ccn_over_land
          else if (xland(i) .eq. 2.0) then
              ccn_ref = ccn_over_water
          endif
          ! for cloud water, estimate lambda (slope of gamma distribution)
          mu = min(15.e0, (1000.e0/ccn_ref + 2.e0))
          gamfac3 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+3.e0) )
          gamfac1 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+1.e0) )
          lambda = (4.e0/3.e0*cpi*roqr*ccn_ref/l_cloud* gamfac1)**(1.e0/3.e0)  ! [1/cm]
          refc2d(i,k) = 1.e0/lambda * gamfac3 * 1.e4  !effective radius [micron]
#endif
      endif ! qcl(i,k) < cmin test
     

   enddo !i-loop

 enddo giant_kloop

!
! replace greg thompson's dbz values with those calculated
! above.
!
 if ( present (diagflag) ) then
    if (diagflag .and. do_radar_ref == 1) then
        do k=1, kmax 
!dir$ vector aligned
           do i=1,irestrict
                refl_10cm(i,k) = max(-35.,dbz(i,k))
           enddo
        end do

   endif
 endif

!updated time step
 istep = istep + 1
 
 return
 end subroutine saticel_s
 
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
 
  subroutine auto_conversion_liu_daum( l, n, p , re)
  implicit none
!-----------------------------------------------------------------------------------------------------
! comments:
!  this subroutine compute auto conversion rate folloing Liu and Daum [2004], which account for
!  total cloud liquid water, particle number concentrations, psd lambda, broadening parameters.
!
! history:
!  08/2010  toshi matsui@nasa gsfc : initial.
!
! references:
! Liu, y. and p. h. daum, 2004: parameterization of the autoconversion process. part i: analytical
!   formulation of the kessler-type parameterizations. j. atmos. sci, 61, 1539-1548.
!-----------------------------------------------------------------------------------------------------
 real,intent(in) :: l    ! cloud liquid water [g cm-3]
 real,intent(in) :: n    ! total number concentration [# cm-3]
 real,intent(out) :: p   ! auto conversion rate [g cm-3 s-1]
 real,intent(out) :: re  ! cloud effective radius [micron]

 real :: mu   ! mu of gamma psd [-]
 real :: eta  ! eta function [cm3 g-2 s-1]
 real :: beta, beta1, beta2     ! beta function [-]
 real :: gamfac , gfac1 , gfac2 ! gamma function [-]
 real :: r6_6power  ! mean radius of the sixth moment [cm]
 real :: r6_thresh ! threshold of  mean radius of the sixth moment [cm]
 real :: r6        ! mean radius of the sixth moment [cm]
 real :: heaviside_func  ! heaviside step function (0 or 1)
 real :: lambda    ! slope of gamma size ditribution [1/cm]
! real :: no        ! intercept  [cm-4]

 real,parameter :: rc = 10.e0 * 1.e-4 ! threshold of particle radus (10 micron) [cm]
 real,parameter :: const_pi    = 3.14159e0 ! pai
 real,parameter :: const_kappa = 1.9e11    ! coefficient for water droplet collection kernel [cm-3 s-1]
                                           ! from long [1974, jas].
! real,parameter :: const_kappa = 1.9e11*10000.e0 !10000 is to adjust the order to keseller


 real,parameter :: const_rho_liq = 1.e0    ! density of liquid water [g cm-3]
 real,parameter :: eta_func = ((3.e0/(4.e0*const_pi*const_rho_liq))**2) * const_kappa  ! eta function [cm3 g-2 s-1]
                                                                                       ! a part of (eq 27b)

 logical,parameter :: no_thresh = .true.  ! logic to choose no threshold parameterization or not.

!
! when no particel, no autoconversion.
!
! emk bug fix...prevent overflow for small but non-zero values of l
! if( n <= 0.e0 .or. l <= 0.e0 ) then
 if( n <= 0.e0 .or. l <= 1.0e-32 ) then
   p = 0.e0
   return
 endif

!
! check bad values of n and l
!
! if( n < 0.e0 ) stop 'msg auto_conversion_liu_daum: n is negative, it must be positive'
! if( l < 0.e0 ) stop 'msg auto_conversion_liu_daum: l is negative, it must be positive'
! if( l > 1.e0 ) stop 'msg auto_conversion_liu_daum: l is greater than 1g/cm3.'

!
! empirical fit of mu as a function of total particle number concentrations.
! from martin et al. (1994), assign gamma shape parameter mu for cloud
! drops according to general dispersion characteristics.
! disp=~0.25 for maritime and 0.45 for continental.
! since disp=sqrt((mu+2)/(mu+1) - 1), mu varies from 15 for maritime (pristine air)
! to 2 for continental (really dirty air).  if mu = 0 --> expnential distribution (narrow dist)
!

! orig
 mu = min(15.e0, (1000.e0/n + 2.e0))


!
! gamma functions
!
 gfac1 = gamma_toshi(mu+4.e0)
 gfac2 = gamma_toshi(mu+1.e0)
 gamfac = (gfac1/gfac2)

!
! estimate lambda (slope of gamma distribution)
!
 lambda = (4.e0/3.e0*const_pi*const_rho_liq*n/l*gamfac)**(1.e0/3.e0)  ! [1/cm]


 thresh: if( no_thresh ) then !-------------------------------------------

!
! threshold of particle radius (mean radius of the sixth moment )
!
 gfac1 = gamma_toshi(mu+7.e0)
 gfac2 = gamma_toshi(mu+1.e0)
 gamfac = (gfac1/gfac2)

 r6_6power = (1.e0 / lambda)**6.e0 * gamfac   ![cm] (eq. a3)

!
! auto conversion rate (eq. 26a)
!
 p = const_kappa *   n     * r6_6power *    l      ! [g cm-3 s-1 ]
!    [cm-3 s-1]  * [#/cm3] *   [cm6]   * [g/cm3]


 else  !with threshold ---------------------------------------------------

!
! estimate eta under gamma psd
!
 beta1 = (6.e0+mu)*(5.e0+mu)*(4.e0+mu)
 beta2 = (3.e0+mu)*(2.e0+mu)*(1.e0+mu)
 beta  = beta1 / beta2

 eta = eta_func * beta  ! eta function (eq 27b) [cm3 g-2 s-1]

!
! threshold of particle radius (mean radius of the sixth moment )
!
 r6_thresh = beta * rc  ![cm] (pg 1545)

!
! mean radius of the sixth moment
!
 gfac1 = gamma_toshi(6.e0+mu+1.e0)
 gfac2 = gamma_toshi(1.e0+mu)
 gamfac = (gfac1/gfac2)**(1.e0/6.e0)

 r6 = (1.e0 / lambda) * gamfac   ![cm] (eq. a3)

!
! heaviside step function
!
 if    ( r6 - r6_thresh <= 0.e0 ) then
    heaviside_func = 0.e0
 elseif( r6 - r6_thresh > 0.e0 ) then
    heaviside_func = 1.e0
 else
    ! nuwrf emk...user wrf's library to gracefully stop mpi.
    write(wrf_err_message,*)'msg: auto_conversion_liu_daum: strange value of r6= ',r6
    call wrf_error_fatal(trim(wrf_err_message))
!   print*, 'msg: auto_conversion_liu_daum: strange value of r6= ', r6 ; stop
 endif

!
! auto conversion rate [g cm-3 s-1 ] (eq. 27a)
!
 p = eta * (1.e0/n) * (l**3) *  heaviside_func

!    [cm3 g-2 s-1] * [cm3] * [g3/cm9]

 endif thresh !------------------------------------------------------------


! optional

!
! estimate effective radius
!
 gfac1 = gamma_toshi(mu+4.e0)
 gfac2 = gamma_toshi(mu+3.e0)
 gamfac = (gfac1/gfac2)

 re = 1.e0/lambda * gamfac * 1.e4  !effective radius [micron]

!
! estimate no
!
! call gamma_function(mu+1.e0 ,gfac1)
! no = n * (lambda**(mu+1)) / gfac1

 return
 end subroutine auto_conversion_liu_daum

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

!dir$ attributes forceinline :: gamma_toshi
 real function gamma_toshi(x)

!---------------------------------------------------------------------------------------------------
! comments:
!   compute the gamma function t(x) for single precision floating point.
!       input :  x  --- argument of a(x)
!                       ( x is not equal to 0,-1,-2,... )
!
! history:
! 09/2009  toshi matsui@nasa gsfc ; adapted from sdsu.  
!
! references:
!----------------------------------------------------------------------------------------------------
 implicit double precision (a-h,o-z)
 dimension g(26)
 data g/1.0d0,0.5772156649015329d0, &
       -0.6558780715202538d0, -0.420026350340952d-1, &
        0.1665386113822915d0,-.421977345555443d-1, &
        -.96219715278770d-2, .72189432466630d-2, &
        -.11651675918591d-2, -.2152416741149d-3, &
        .1280502823882d-3, -.201348547807d-4, &
        -.12504934821d-5, .11330272320d-5, &
        -.2056338417d-6, .61160950d-8, &
         .50020075d-8, -.11812746d-8, &
        .1043427d-9, .77823d-11, &
        -.36968d-11, .51d-12, &
        -.206d-13, -.54d-14, .14d-14, .1d-15/
 real :: x

 integer :: m1, k, m

 pi=3.141592653589793d0
 if (x.eq.int(x)) then
     if (x.gt.0.0d0) then
         ga=1.0d0
         m1=int(x)-1
        do k=2,m1
           ga=ga*k
        enddo
     else
        ga=1.0d+300
     endif
  else
     if (dabs(dble(x)).gt.1.0d0) then
         z=dabs(dble(x))
         m=int(z)
         r=1.0d0
        do k=1,m
           r=r*(z-k)
        enddo
        z=z-m
     else
        z=dble(x)
     endif
     gr=g(26)
     do k=25,1,-1
        gr=gr*z+g(k)
     enddo
     ga=1.0d0/(gr*z)
     if (dabs(dble(x)).gt.1.0d0) then
         ga=ga*r
         if (x.lt.0.0d0) ga=-pi/(x*ga*dsin(pi*x))
     endif
  endif

  gamma_toshi = real(ga)

  return
  end function gamma_toshi

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine find_nan_inf_double(warning_msg, real_input, i_in,j_in,k_in)
 implicit none
!------------------------------------------------------------------------------
! Comments: 
!    Find bad values. 
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; initial
!
! References: 
!                                                                      
!------------------------------------------------------------------------------
 real,intent(inout) :: real_input  !anykind of non-dimensional input real parameters
 integer,intent(in) :: i_in, j_in, k_in
 character*(*),intent(in) :: warning_msg

!
! find infinity
!
!if( exp(-abs(real_input)) == 0.) then ! this formulae is bit slow 

 if( 1e+10/real_input == 0. ) then
    print*,'msg find_nan_inf: '//warning_msg//'infinity at',i_in,j_in,k_in
    real_input = 0.
    return
 endif

!
! find nan
! 
 if( real_input==0. .or. real_input>0. .or. real_input<0. .or. real_input>=0. .or. real_input<=0. ) then
 else
    print*,'msg find_nan_inf: '//warning_msg//'nan at',i_in,j_in,k_in
    real_input = 0.
    return
 endif

 return
 end subroutine find_nan_inf_double

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine refl10cm_gsfc (qv1d, qr1d, qs1d, qg1d,                 &
                           t1d, p1d, dbz, kts, kte, ii, jj)

 implicit none
!------------------------------------------------------------------------------
! Comments: 
!  compute radar refrlectivity (10cm - S-band) following the 4ICE mapping. 
!  1D vertical (k-loop) driver. 
!
! History:
! 09/2018  Toshi Matsui@NASA GSFC ; modernization / clean up. 
!
! References: 
!                                                                      
!------------------------------------------------------------------------------


!..sub arguments
 integer, intent(in):: kts, kte, ii, jj
 real, dimension(kts:kte), intent(in)::                            &
      qv1d, qr1d, qs1d, qg1d, t1d, p1d
 real, dimension(kts:kte), intent(inout):: dbz

!..local variables
 real, dimension(kts:kte):: temp, pres, qv, rho
 real, dimension(kts:kte):: rr, rs, rg

 double precision, dimension(kts:kte):: ilamr, ilams, ilamg
 double precision, dimension(kts:kte):: n0_r, n0_s, n0_g
 double precision:: lamr, lams, lamg
 logical, dimension(kts:kte):: l_qr, l_qs, l_qg

 real, dimension(kts:kte):: ze_rain, ze_snow, ze_graupel
 double precision:: fmelt_s, fmelt_g

 integer:: i, k, k_0, kbot, n
 logical:: melti

 double precision:: cback, x, eta, f_d
 real, parameter:: r=287.
 real, parameter:: pix=3.1415926536

!+---+

      do k = kts, kte
         dbz(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = max(1.e-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(r*temp(k)*(qv(k)+0.622))

         if (qr1d(k) .gt. 1.e-9) then
            rr(k) = qr1d(k)*rho(k)
            n0_r(k) = xnor
            lamr = (xam_r*xcrg(3)*n0_r(k)/rr(k))**(1./xcre(1))
            ilamr(k) = 1./lamr
            l_qr(k) = .true.
         else
            rr(k) = 1.e-12
            l_qr(k) = .false.
         endif

         if (qs1d(k) .gt. 1.e-9) then
            rs(k) = qs1d(k)*rho(k)
            n0_s(k) = xnos
            lams = (xam_s*xcsg(3)*n0_s(k)/rs(k))**(1./xcse(1))
            ilams(k) = 1./lams
            l_qs(k) = .true.
         else
            rs(k) = 1.e-12
            l_qs(k) = .false.
         endif

         if (qg1d(k) .gt. 1.e-9) then
            rg(k) = qg1d(k)*rho(k)
               n0_g(k) = xnog
            lamg = (xam_g*xcgg(3)*n0_g(k)/rg(k))**(1./xcge(1))
            ilamg(k) = 1./lamg
            l_qg(k) = .true.
         else
            rg(k) = 1.e-12
            l_qg(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!..locate k-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt.273.15) .and. l_qr(k)                         &
                                  .and. (l_qs(k+1).or.l_qg(k+1)) ) then
            k_0 = max(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue

!+---+-----------------------------------------------------------------+
!..assume rayleigh approximation at 10 cm wavelength. rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. integrations of m(d)*m(d)*n(d)*dd.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (l_qr(k)) ze_rain(k) = n0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (l_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/pix)*(6.0/pix)     &
                                 * (xam_s/900.0)*(xam_s/900.0)          &
                                 * n0_s(k)*xcsg(4)*ilams(k)**xcse(4)
         if (l_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/pix)*(6.0/pix)  &
                                    * (xam_g/900.0)*(xam_g/900.0)       &
                                    * n0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
      enddo


!+---+-----------------------------------------------------------------+
!..special case of melting ice (snow/graupel) particles.  assume the
!.. ice is surrounded by the liquid water.  fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. uses code from uli blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (melti .and. k_0.ge.kts+1) then
       do k = k_0-1, kts, -1

!..reflectivity contributed by melting snow
          if (l_qs(k) .and. l_qs(k_0) ) then
           fmelt_s = max(0.005d0, min(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           lams = 1./ilams(k)
           do n = 1, nrbins
              x = xam_s * xxds(n)**xbm_s
              call rayleigh_soak_wetgraupel (x,dble(xocms),dble(xobms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    cback, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = n0_s(k)*xxds(n)**xmu_s * dexp(-lams*xxds(n))
              eta = eta + f_d * cback * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = sngl(lamda4 / (pi5 * k_w) * eta)
          endif

!..reflectivity contributed by melting graupel

          if (l_qg(k) .and. l_qg(k_0) ) then
           fmelt_g = max(0.005d0, min(1.0d0-rg(k)/rg(k_0), 0.99d0))
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = xam_g * xxdg(n)**xbm_g
              call rayleigh_soak_wetgraupel (x,dble(xocmg),dble(xobmg), &
                    fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                    cback, mixingrulestring_g, matrixstring_g,          &
                    inclusionstring_g, hoststring_g,                    &
                    hostmatrixstring_g, hostinclusionstring_g)
              f_d = n0_g(k)*xxdg(n)**xmu_g * dexp(-lamg*xxdg(n))
              eta = eta + f_d * cback * simpson(n) * xdtg(n)
           enddo
           ze_graupel(k) = sngl(lamda4 / (pi5 * k_w) * eta)
          endif

       enddo
      endif

      do k = kte, kts, -1
         dbz(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo

 end subroutine refl10cm_gsfc

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 real function eff_rad(lambda)
 implicit none
!---------------------------------------------------------------------------------------------------
! comments:
! compute drop effective radius from slope parameters (lambda) of expoential size distribution.
! Scalar driver. 
!
! history:
! 02/2014  toshi matsui@nasa gsfc ; initial
!
! references:
!----------------------------------------------------------------------------------------------------
 real,intent(in)  :: lambda   ! intercept parameter [1/cm]
!
! for no particles.
!
 if ( lambda <= 0.e0  .or. isnan(lambda) ) then
      eff_rad = 0.e0
      return
 endif

!
! compute drop effective radius for exponential distribution n(d) = n0*exp(-lam*d)
!
  eff_rad = 1.5e0 / (lambda*100.) * 1.0e+6  ! [micron]

 return
 end function eff_rad

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine goddard_microphysics_gce
 implicit none
!
! This subroutine is wired with GCE microphysics driver
!

!
! Unified 1D driver of microphysics
!
 call goddard_microphysics_1d

 end subroutine goddard_microphysics_gce

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine goddard_microphysics_init( allowed_to_read, dt, &
                         ids, ide, jds, jde, kds, kde,       &
                         ims, ime, jms, jme, kms, kme,       &
                         its, ite, jts, jte, kts, kte  )

 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  One-time initialize goddard unified microphysics module.   
!  this initialization can be called either from WRF or GCE.
! 
! History: 08/2018  Toshi Matsui@NASA GSFC ;  iinit. 
!           
! References: 
!---------------------------------------------------------------------------------------------------
 logical,intent(in)  :: allowed_to_read !for MPI master processor or not
 real,intent(in)     :: dt ! time step 
 integer,intent(in)  :: ids, ide, jds, jde, kds, kde,  &  ! WRF domain bounds
                        ims, ime, jms, jme, kms, kme,  &  ! WRF memory bounds
                        its, ite, jts, jte, kts, kte      ! WRF tile bounds

!
! intialization of microphysics constant
!
 call consat_s

!
! run-time option
!
 choice_moment  = '1moment'          ! '2moment'
 turn_on_electrification = .true.    ! if true, it will account for electric charing mechanism.
 turn_on_aerosol  = .false.          ! if true, it track in-cloud aerosol tracing process


 return
 end subroutine goddard_microphysics_init

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine goddard_microphysics_init_gce ( allowed_to_read, dt )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  One-time initialize goddard unified microphysics module.   
!  this initialization can be called either from WRF or GCE.
! 
! History: 08/2018  Toshi Matsui@NASA GSFC ;  iinit. 
!           
! References: 
!---------------------------------------------------------------------------------------------------
 logical,intent(in)  :: allowed_to_read !for MPI master processor or not
 real,intent(in)     :: dt ! time step 

!
! intialization of microphysics constant
!
 call consat_s

!
! run-time option
!
 choice_moment  = '1moment'          ! '2moment'
 turn_on_electrification = .true.    ! if true, it will account for electric charing mechanism.
 turn_on_aerosol  = .false.          ! if true, it track in-cloud aerosol tracing process

 return
 end subroutine goddard_microphysics_init_gce

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine goddard_microphysics_1d
 implicit none

! this is identically called from WRF and GCE coupler


 end subroutine goddard_microphysics_1d

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 end module module_mp_goddard_unified

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

