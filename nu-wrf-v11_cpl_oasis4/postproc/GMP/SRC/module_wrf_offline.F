#include <define.h>
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
 module module_wrf_offline
 use netcdf
 use module_runtime_param 
 use module_mpi
! use module_omp_check

#if ( WRF_CHEM == 1)
 use module_gocart_coupling
#endif

!======================================================================
!<description>
!
!  WRF radiation test module
! 
!======================================================================
! grid structure in physics part of wrf
!----------------------------------------------------------------------
! the horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. all predicted
! variables are carried at half levels except w, which is at full
! levels. some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! in wrf, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  in your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!==================================================================
 implicit none

 
 save  ! all module parameters are saved. 


!
! dimesion bounds
!
  integer,save :: ids,ide,jds,jde,kds,kde     ! all (in)
  integer,save :: ims,ime,jms,jme,kms,kme     ! all (in)
  integer,save :: its,ite,jts,jte,kts,kte     ! all (in)


!
! hydrometeor index
!
 integer,parameter :: id_hmax    = 6  ! total number of  hydrometeor species id maxima
 integer,parameter :: id_ice     = 1
 integer,parameter :: id_cloud   = 2
 integer,parameter :: id_rain    = 3
 integer,parameter :: id_snow    = 4
 integer,parameter :: id_graupel = 5
 integer,parameter :: id_hail    = 6


!
! Constant physics parameters
!
 real,parameter :: degrad=3.14159/180.  ! conversion factor for degrees to radians (pi/180.) (rad/deg) 
 real,parameter :: r_d = 287.           ! gas constant for dry air (j/kg/k)
 real,parameter :: cp  = 717+r_d        ! heat capacity at constant pressure for dry air (j/kg/k)
 real,parameter :: g   = 9.8            ! acceleration due to gravity (m/s^2)
 real,parameter :: dpd = 360./365       ! degrees per day for earth's orbital position (deg/day) 
 real,parameter :: thresh = 1.0e-6      ! cloud detection threshold (tunable!!) orig = 1.0e-6
 real,parameter :: const_pi = 3.14159
 real,parameter :: const_Rd = 287.053 ! dry-air gas constant [J/K/kg] 
 real,parameter :: rhowater = 1.0e+3   !  water density      (kg/m^3)
 real,parameter :: rhosnow  = 0.917e+3 !  snow density       (kg/m^3)


 real :: ts,te !start / end time 

 real,parameter :: q_min_condensate = 1.e-6

  character :: CRM_file_list*100        !Input file names, defined in "def_file_list.f".


!
! temporal output parameter for all domain
!
 real,allocatable, dimension(:,:)   ::  out_domain2d
 real,allocatable, dimension(:,:,:) ::  out_domain3d
 real,allocatable, dimension(:,:)   ::  out_memory_id
 real,allocatable, dimension(:,:)   ::  out_tile_id



!
! WRF input parameters
!
 real,allocatable, dimension(:,:,:) ::  &
     dz8w, & ! dz between full levels (m)     (in)
     z,    & ! height above sea level (m)     (in) 
     p,    & ! pressure (pa)                  (in)
     p8w,  & ! pressure at full levels (pa)   (in)
     pi_phy, & ! exner function (-)             (in)
     t,    & ! temperature (K)                (in)
     th,   & ! potential temperature (K)                (in)
     t8w,  & ! temperature at full levels (k) (in)
     rho,  & ! density (kg/m^3)               (in) 
     w  ,  & ! vertical velocity (m/s) (in)
     qv_curr,   & ! water vapor mixing ratio (kg/kg)  (in)
     qc_curr,   & ! cloud water mixing ratio (kg/kg)  (in)
     qr_curr,   & ! rain water mixing ratio (kg/kg)   (in)
     qi_curr,   & ! cloud ice mixing ratio (kg/kg)    (in)
     qs_curr,   & ! snow mixing ratio (kg/kg)         (in)
     qg_curr,   & ! graupel mixing ration (kg/kg)     (in)
     qh_curr,   & ! hail mixing ratio (kg/kg)     (in)
     qni_curr,   & ! cloud ice number conc. [#/kg]
     pion_curr,   & ! positive ion (Nb/kg)
     nion_curr,   & ! negative ion (Nb/kg)
     cmrc_curr,   & ! cloud charge mixing ratio (C/kg)  
     cmrr_curr,   & ! rain water charge mixing ratio (C/kg)   
     cmri_curr,   & ! cloud ice charge mixing ratio (C/kg)    
     cmrs_curr,   & ! snow charge mixing ratio (C/kg)         
     cmrg_curr,   & ! graupel charge mixing ration (C/kg)     
     cmrh_curr      ! hail charge mixing ration (C/kg)  

 real,allocatable, dimension(:,:,:) ::  &
    refc,   & ! cloud effective radius [um]  (in)
    refr,   & ! rain effective radius [um]   (in)
    refi,   & ! cloud ice effective radius [um]    (in)
    refs,   & ! snow effective radius [um]  (in)
    refg,   & ! graupel effective radius [um]  (in)
    refh      ! hail effective radius [um] (in)

 real,allocatable, dimension(:,:,:) ::  &
    refl_10cm, &    ! radar reflectivty [dBZ]
    phys_tot,  &    ! [K/s]
    physc,     &
    physe,     &
    physd,     &
    physs,     &
    physm,     &
    physf,     &
    acphys_tot,  &
    acphysc,     &
    acphyse,     &
    acphysd,     &
    acphyss,     &
    acphysm,     &
    acphysf,     &
         precr3d, & ! rain precipitation rate at all levels (kg/m2/s)
         preci3d, & ! ice precipitation rate at all levels (kg/m2/s)
         precs3d, & ! snow precipitation rate at all levels (kg/m2/s)
         precg3d, & ! graupel precipitation rate at all levels (kg/m2/s)
         prech3d    ! graupel precipitation rate at all levels (kg/m2/s)


 logical,parameter :: diagflag = .true. ! 
 integer,parameter :: do_radar_ref = 1

!
! WRF Scalar variables
!
   integer :: julday    ! julian day (in)
   real :: julian ! julian day in real                             (in)
   real :: dt     !time step (sec)                (in)
   real :: gmt    ! greenwich mean time hour of model start (hour) (in)
   real :: xtime  ! time since simulation start (min)              (in)

!
! WRF 2D variables
!
   real, allocatable, dimension( :,: ) :: &
      xlat,   & ! latitude, south is negative (degree)  (in)-
      xlong,  & ! longitude, west is negative (degree)  (in)-
      xland,  & ! land-water mask (in)-
      albedo, & ! broadband albedo (between 0 and 1)    (in)-
      ht   ,  & ! terrain height [m] (in) - 
      tsk  ,  & ! surface skin temperature [K] (in) -
      rainnc, & ! grid scale precipitation (mm)
      rainncv,& ! one time step grid scale precipitation (mm/step)
      snownc, & ! grid scale snow and ice (mm)
      snowncv, & ! one time step grid scale snow and ice (mm/step)
      graupelnc, & ! grid scale graupel precip (mm)
      graupelncv,& ! one time step grid scale graupel precip (mm/step)
      hailnc    ,& ! grid scale hail precip (mm)
      hailncv   ,& ! one time step grid scale hail precip (mm/step)
      sr        ,& ! one time step mass ratio of snow to total precip
      coszen  

 integer,allocatable :: ivgtyp(:,:) !dominant vegetation type


 real,allocatable :: var2d(:,:)   ! surface pressure at interface [Pa] 
 real,allocatable, dimension(:,:,:) ::  &
     var3d     , & ! temporally 3D variables
     var3d_stag    ! temp 3D stagger


#if ( WRF_CHEM == 1)
  integer, parameter:: num_go = 14         ! # of gocart aerosol species (tgmx=14 defined in module_gocart_offline.F)
  real, allocatable :: aero( :, :, :, :)   !aerosol mass conc [g/m3]            

#endif


!
! ------ interface ------------
!

 interface dump
    module procedure dump2d
    module procedure dump3d
 end interface dump


 contains

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine allocate_all
 implicit none
 integer :: ierr

 if(masterproc) print*, 'MSG allocate_all: define dimesion and allocate global parameters'

!
! set up all dimension here. 
!

    ids=ids_wrf    ! start index for i in domain
    ide=ide_wrf    ! end index for i in domain
    jds=jds_wrf    ! start index for j in domain
    jde=jde_wrf    ! end index for j in domain
    kds=1          ! start index for k in domain
    kde=kmax       ! end index for k in domain

    ims=myi_start  ! start index for i in memory
    ime=myi_end    ! end index for i in memory
    jms=myj_start  ! start index for j in memory
    jme=myj_end    ! end index for j in memory
    kms=kds        ! start index for k in memory
    kme=kde        ! end index for k in memory

    its=ims  ! start index for i in memory
    ite=ime  ! end index for i in memory
    jts=jms  ! start index for j in memory
    jte=jme  ! end index for j in memory
    kts=kms        ! start index for k in tile
!    kte=kme-1      ! end index for k in tile (toshii- check this...)
    kte=kme

!
! allocation
!
  allocate( &
      out_domain2d( ids:ide, jds:jde ),& !2D domain output
      out_domain3d( ids:ide, jds:jde, kds:kde ),& !3D domain output
      out_memory_id( ids:ide, jds:jde ),& !
      out_tile_id  ( ids:ide, jds:jde ),& !
      xlat    ( ims:ime, jms:jme ),   & ! latitude, south is negative (degree)  (in)-
      xlong   ( ims:ime, jms:jme ),   & ! longitude, west is negative (degree)  (in)-
      xland   ( ims:ime, jms:jme ),   & ! land-water mask (in)-
      ivgtyp  ( ims:ime, jms:jme ),   & ! dominant vegation type
      albedo  ( ims:ime, jms:jme ),   & ! broadband albedo (between 0 and 1)    (in)-
      tsk     ( ims:ime, jms:jme ),   & ! surface skin temperature [K] (in) -
      rainnc  ( ims:ime, jms:jme ), & ! RAINNC &
      rainncv ( ims:ime, jms:jme ),& ! RAINNCV &
      snownc  ( ims:ime, jms:jme ), & ! SNOWNC &
      snowncv ( ims:ime, jms:jme ), & ! SNOWNCV &
      graupelnc( ims:ime, jms:jme ), & ! GRAUPELNC &
      graupelncv( ims:ime, jms:jme ),& ! GRAUPELNCV &
      hailnc    ( ims:ime, jms:jme ) ,& ! HAILNC &
      hailncv   ( ims:ime, jms:jme ) ,& ! HAILNCV
      sr        ( ims:ime, jms:jme ) ,& !
      coszen  ( ims:ime, jms:jme ),   &      !
      ht      ( ims:ime, jms:jme ),   &      ! terrain height (m) (toshii double check)
      dz8w       (ims:ime,kms:kme,jms:jme), & ! dz between full levels (m)     (in)
      z          (ims:ime,kms:kme,jms:jme), & ! height above sea level (m)     (in)
      p          (ims:ime,kms:kme,jms:jme), & ! pressure (pa)                  (in)
      p8w        (ims:ime,kms:kme,jms:jme), & ! pressure at full levels (pa)   (in)
      pi_phy     (ims:ime,kms:kme,jms:jme), & ! exner function (-)             (in)
      th         (ims:ime,kms:kme,jms:jme), & ! potential temperature (K)      (in)
      t          (ims:ime,kms:kme,jms:jme), & ! temperature (K)                (in)
      t8w        (ims:ime,kms:kme,jms:jme), & ! temperature at full levels (k) (in)
      rho        (ims:ime,kms:kme,jms:jme), & ! density (kg/m^3)               (in) 
      w          (ims:ime,kms:kme,jms:jme), & ! vertical velocity (m/s)        (in) 
      qv_curr         (ims:ime,kms:kme,jms:jme), & ! water vapor mixing ratio (kg/kg)  (in)
      qc_curr         (ims:ime,kms:kme,jms:jme), & ! cloud water mixing ratio (kg/kg)  (in)
      qr_curr         (ims:ime,kms:kme,jms:jme), & ! rain water mixing ratio (kg/kg)   (in)
      qi_curr         (ims:ime,kms:kme,jms:jme), & ! cloud ice mixing ratio (kg/kg)    (in)
      qs_curr         (ims:ime,kms:kme,jms:jme), & ! snow mixing ratio (kg/kg)         (in)
      qg_curr         (ims:ime,kms:kme,jms:jme), & ! graupel mixing ration (kg/kg)     (in)
      qh_curr         (ims:ime,kms:kme,jms:jme), & ! hail mixing ration (kg/kg)     (in)
      qni_curr         (ims:ime,kms:kme,jms:jme),& ! cloud ice nuber conc. [#/kg]
     pion_curr      (ims:ime,kms:kme,jms:jme)  , & ! positive ion (Nb/kg)
     nion_curr      (ims:ime,kms:kme,jms:jme)  ,   & ! negative ion (Nb/kg)
     cmrc_curr      (ims:ime,kms:kme,jms:jme)  ,   & ! cloud water mixing ratio (kg/kg)  (in)
     cmrr_curr      (ims:ime,kms:kme,jms:jme)  ,   & ! rain water mixing ratio (kg/kg)   (in)
     cmri_curr      (ims:ime,kms:kme,jms:jme)  ,   & ! cloud ice mixing ratio (kg/kg)    (in)
     cmrs_curr      (ims:ime,kms:kme,jms:jme)  ,   & ! snow mixing ratio (kg/kg)         (in)
     cmrg_curr      (ims:ime,kms:kme,jms:jme)  ,   & ! graupel mixing ration (kg/kg)     (in)
     cmrh_curr      (ims:ime,kms:kme,jms:jme)  ,   & ! hail mixing ratio (kg/kg)     (in)
     refc        (ims:ime,kms:kme,jms:jme), & ! cloud effective radius [um] (in)
     refr        (ims:ime,kms:kme,jms:jme), & ! rain effective radius [um]   (in)
     refi        (ims:ime,kms:kme,jms:jme), & ! cloud effective radius [um]    (in)
     refs        (ims:ime,kms:kme,jms:jme), & ! snow effective radius [um]  (in)
     refg        (ims:ime,kms:kme,jms:jme), & ! graupel effective radius [um]   (in)
     refh        (ims:ime,kms:kme,jms:jme), & ! hail effective radius [um]    (in) 
    refl_10cm    (ims:ime,kms:kme,jms:jme), &    ! radar reflectivty [dBZ]
    phys_tot     (ims:ime,kms:kme,jms:jme),  &   ! total latent heat
    physc        (ims:ime,kms:kme,jms:jme),     & ! latent heat of condensation
    physe        (ims:ime,kms:kme,jms:jme),     & ! latent heat of evaporation
    physd        (ims:ime,kms:kme,jms:jme),     & ! latent heat of deposition
    physs        (ims:ime,kms:kme,jms:jme),     & ! latent heat of 
    physm        (ims:ime,kms:kme,jms:jme),     & ! latent heat of melting
    physf        (ims:ime,kms:kme,jms:jme),     & ! 
    acphys_tot   (ims:ime,kms:kme,jms:jme),  &
    acphysc      (ims:ime,kms:kme,jms:jme),     &
    acphyse      (ims:ime,kms:kme,jms:jme),     &
    acphysd      (ims:ime,kms:kme,jms:jme),     &
    acphyss      (ims:ime,kms:kme,jms:jme),     &
    acphysm      (ims:ime,kms:kme,jms:jme),     &
    acphysf      (ims:ime,kms:kme,jms:jme),     &
     precr3d     (ims:ime,kms:kme,jms:jme), & ! rain precipitation rate at all levels (kg/m2/s)
     preci3d     (ims:ime,kms:kme,jms:jme), & ! ice precipitation rate at all levels (kg/m2/s)
     precs3d     (ims:ime,kms:kme,jms:jme), & ! snow precipitation rate at all levels (kg/m2/s)
     precg3d     (ims:ime,kms:kme,jms:jme), & ! graupel precipitation rate at all levels (kg/m2/s)
     prech3d     (ims:ime,kms:kme,jms:jme), & ! graupel precipitation rate at all levels (kg/m2/s)
      var3d_stag (ims:ime,jms:jme,kms:kme+1), & ! base geopotential
      var3d      (ims:ime,jms:jme,kms:kme),  & ! temporally 3D variables
      var2d      (ims:ime,jms:jme),    & ! surface pressure [Pa]
#if ( WRF_CHEM == 1)
      aero( ims:ime, kms:kme, jms:jme, 1:num_go) , &  ! aerosol mass conc [g/m3]            
#endif
  stat=ierr )

 if(ierr /= 0) stop 'MSG allocate_all: cannot allocate the variable'


!
! zeroed all parameters
!
      out_domain2d( ids:ide, jds:jde ) = 0.  !2D domain output
      out_domain3d( ids:ide, jds:jde, kds:kde ) = 0. !3D domain output
      out_memory_id( ids:ide, jds:jde ) = 0. !
      out_tile_id  ( ids:ide, jds:jde ) = 0. !
      xlat    ( ims:ime, jms:jme ) = 0. ! latitude, south is negative (degree)  (in)-
      xlong   ( ims:ime, jms:jme ) = 0. ! longitude, west is negative (degree)  (in)-
      xland   ( ims:ime, jms:jme ) = 0. ! land-water mask (in)-
      albedo  ( ims:ime, jms:jme ) = 0. ! broadband albedo (between 0 and 1)    (in)-
      tsk     ( ims:ime, jms:jme ) = 0. ! surface skin temperature [K] (in) -
      rainnc  ( ims:ime, jms:jme ) = 0. ! RAINNC &
      rainncv ( ims:ime, jms:jme ) = 0. ! RAINNCV &
      snownc  ( ims:ime, jms:jme ) = 0. ! SNOWNC &
      snowncv ( ims:ime, jms:jme ) = 0. ! SNOWNCV &
      graupelnc( ims:ime, jms:jme ) = 0. ! GRAUPELNC &
      graupelncv( ims:ime, jms:jme ) = 0. ! GRAUPELNCV &
      hailnc    ( ims:ime, jms:jme ) = 0. ! HAILNC &
      hailncv   ( ims:ime, jms:jme ) = 0. ! HAILNCV
      sr        ( ims:ime, jms:jme ) = 0. !
      coszen  ( ims:ime, jms:jme ) = 0.      !
      ht      ( ims:ime, jms:jme ) = 0.      ! terrain height (m) (toshii double check)
      dz8w       (ims:ime,kms:kme,jms:jme) = 0. ! dz between full levels (m)     (in)
      z          (ims:ime,kms:kme,jms:jme) = 0. ! height above sea level (m)     (in)
      p          (ims:ime,kms:kme,jms:jme) = 0. ! pressure (pa)                  (in)
      p8w        (ims:ime,kms:kme,jms:jme) = 0. ! pressure at full levels (pa)   (in)
      pi_phy     (ims:ime,kms:kme,jms:jme) = 0. ! exner function (-)             (in)
      th         (ims:ime,kms:kme,jms:jme) = 0. ! potential temperature (K)      (in)
      t          (ims:ime,kms:kme,jms:jme) = 0. ! temperature (K)                (in)
      t8w        (ims:ime,kms:kme,jms:jme) = 0. ! temperature at full levels (k) (in)
      rho        (ims:ime,kms:kme,jms:jme) = 0. ! density (kg/m^3)               (in) 
      w          (ims:ime,kms:kme,jms:jme) = 0. ! vertical velocity (m/s)        (in) 
      qv_curr         (ims:ime,kms:kme,jms:jme) = 0. ! water vapor mixing ratio (kg/kg)  (in)
      qc_curr         (ims:ime,kms:kme,jms:jme) = 0. ! cloud water mixing ratio (kg/kg)  (in)
      qr_curr         (ims:ime,kms:kme,jms:jme) = 0. ! rain water mixing ratio (kg/kg)   (in)
      qi_curr         (ims:ime,kms:kme,jms:jme) = 0. ! cloud ice mixing ratio (kg/kg)    (in)
      qs_curr         (ims:ime,kms:kme,jms:jme) = 0. ! snow mixing ratio (kg/kg)         (in)
      qg_curr         (ims:ime,kms:kme,jms:jme) = 0. ! graupel mixing ration (kg/kg)     (in)
      qh_curr         (ims:ime,kms:kme,jms:jme) = 0. ! hail mixing ration (kg/kg)     (in)
      qni_curr        (ims:ime,kms:kme,jms:jme) = 0. ! cloud ice number conc [#/kg]
     pion_curr    (ims:ime,kms:kme,jms:jme)     = 0. ! positive ion (Nb/kg)
     nion_curr    (ims:ime,kms:kme,jms:jme)     = 0. ! negative ion (Nb/kg)
     cmrc_curr      (ims:ime,kms:kme,jms:jme)   = 0. ! cloud water mixing ratio (kg/kg)  (in)
     cmrr_curr      (ims:ime,kms:kme,jms:jme)   = 0. ! rain water mixing ratio (kg/kg)   (in)
     cmri_curr      (ims:ime,kms:kme,jms:jme)   = 0. ! cloud ice mixing ratio (kg/kg)    (in)
     cmrs_curr      (ims:ime,kms:kme,jms:jme)   = 0. ! snow mixing ratio (kg/kg)         (in)
     cmrg_curr      (ims:ime,kms:kme,jms:jme)   = 0. ! graupel mixing ration (kg/kg)     (in)
     cmrh_curr      (ims:ime,kms:kme,jms:jme)   = 0. ! hail mixing ratio (kg/kg)     (in)
     refc        (ims:ime,kms:kme,jms:jme) = 0. ! cloud effective radius [um] (in)
     refr        (ims:ime,kms:kme,jms:jme) = 0. ! rain effective radius [um]   (in)
     refi        (ims:ime,kms:kme,jms:jme) = 0. ! cloud effective radius [um]    (in)
     refs        (ims:ime,kms:kme,jms:jme) = 0. ! snow effective radius [um]  (in)
     refg        (ims:ime,kms:kme,jms:jme) = 0. ! graupel effective radius [um]   (in)
     refh        (ims:ime,kms:kme,jms:jme) = 0. ! hail effective radius [um]    (in) 
    refl_10cm    (ims:ime,kms:kme,jms:jme) = 0.    ! radar reflectivty [dBZ]
    phys_tot     (ims:ime,kms:kme,jms:jme) = 0.   ! total latent heat
    physc        (ims:ime,kms:kme,jms:jme) = 0. ! latent heat of condensation
    physe        (ims:ime,kms:kme,jms:jme) = 0. ! latent heat of evaporation
    physd        (ims:ime,kms:kme,jms:jme) = 0. ! latent heat of deposition
    physs        (ims:ime,kms:kme,jms:jme) = 0. ! latent heat of 
    physm        (ims:ime,kms:kme,jms:jme) = 0. ! latent heat of melting
    physf        (ims:ime,kms:kme,jms:jme) = 0. ! 
    acphys_tot   (ims:ime,kms:kme,jms:jme) = 0.
    acphysc      (ims:ime,kms:kme,jms:jme) = 0.
    acphyse      (ims:ime,kms:kme,jms:jme) = 0.
    acphysd      (ims:ime,kms:kme,jms:jme) = 0.
    acphyss      (ims:ime,kms:kme,jms:jme) = 0.
    acphysm      (ims:ime,kms:kme,jms:jme) = 0.
    acphysf      (ims:ime,kms:kme,jms:jme) = 0.
     precr3d     (ims:ime,kms:kme,jms:jme) = 0. ! rain precipitation rate at all levels (kg/m2/s)
     preci3d     (ims:ime,kms:kme,jms:jme) = 0. ! ice precipitation rate at all levels (kg/m2/s)
     precs3d     (ims:ime,kms:kme,jms:jme) = 0. ! snow precipitation rate at all levels (kg/m2/s)
     precg3d     (ims:ime,kms:kme,jms:jme) = 0. ! graupel precipitation rate at all levels (kg/m2/s)
     prech3d     (ims:ime,kms:kme,jms:jme) = 0. ! graupel precipitation rate at all levels (kg/m2/s)
      var3d_stag (ims:ime,jms:jme,kms:kme+1) = 0. ! base geopotential
      var3d      (ims:ime,jms:jme,kms:kme) = 0. ! temporally 3D variables
      var2d      (ims:ime,jms:jme) = 0. ! surface pressure [Pa]
#if ( WRF_CHEM == 1)
      aero( ims:ime, kms:kme, jms:jme, 1:num_go) = 0.  ! aerosol mass conc [g/m3]            
#endif



 return
 end subroutine allocate_all

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine read_netcdf
 implicit none
 

 character(len=200) :: fname_in
 integer :: i,j,k  !loop
 integer :: ncid                ! # of variables to inquire  
 integer :: varid               ! # of variable id
 integer :: ierr                ! allocation stat

 real    :: press_top             ! model top pressure at interface [Pa] 
 integer :: is,ie,js,je,ks,ke,di,dj,dk


 real,  dimension( kms:kme, id_hmax ) :: &
         q1d ,& ! condensate mass concentration [g/m3]
         re1d   ! particle effective radius [um]

 real,  dimension( kms:kme ) :: t1d



 is=ims ; ie=ime ; js=jms ; je=jme ; ks=kms ; ke=kme
 di=ie-is+1 ; dj = je-js+1 ; dk=ke-ks+1

!------------------program start---------------------------

 fname_in = trim(CRMdir)//trim(CRM_file_list) 
 if(masterproc) print*,'MSG read_netcdf: Read ',fname_in

 call check( nf90_open(trim(fname_in), nf90_nowrite, ncid) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!  Scalar variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 julday = 218
 julian = real(julday)
 dt  = 100.
 gmt = 12.
 
 call check( nf90_inq_varid(ncid, 'XTIME' , varid ) )  ! minutes since simulation start
 call check( nf90_get_var(ncid, varid, xtime ) )
 xtime = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2D variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---- surface skin temperature 
 call check( nf90_inq_varid(ncid, 'TSK' , varid ) )  ! surface skin temperature [K]
 call check( nf90_get_var(ncid, varid, tsk(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))

!---- lat and lon of grid position
 call check( nf90_inq_varid(ncid, 'XLAT' , varid ) )  ! LATITUDE, SOUTH IS NEGATIVE [deg]
 call check( nf90_get_var(ncid, varid, xlat(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'XLONG' , varid ) )  ! LONGITUDE, WEST IS NEGATIVE [deg]
 call check( nf90_get_var(ncid, varid, xlong(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'XLAND' , varid ) )  ! land-water mask 
 call check( nf90_get_var(ncid, varid, xland(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))

 call check( nf90_inq_varid(ncid, 'IVGTYP' , varid ) )  ! dominant vegetation type
 call check( nf90_get_var(ncid, varid, ivgtyp(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))

!----albedo
 call check( nf90_inq_varid(ncid, 'ALBEDO' , varid ) )  ! BROADBAND ALBEDO  
 call check( nf90_get_var(ncid, varid, albedo(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))

!
 call check( nf90_inq_varid(ncid, 'HGT' , varid ) )  ! surface elevation [m]  
 call check( nf90_get_var(ncid, varid, ht(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))


!!!!!!!!!  3D variable -> first read (i,j,k) then transform (i,k,j)  !!!!!!!!!!!!!!!!
!----height and all sorts
 call check( nf90_inq_varid(ncid, 'PHB' , varid ) ) ! geopotential [m/s2] 
 call check( nf90_get_var(ncid, varid, var3d_stag(is:ie,js:je,ks:ke+1), start=(/is,js,ks/), count=(/di,dj,dk+1/)))

 do k = kms, kme
    var3d(:,:,k) = (var3d_stag(:,:,k+1))/9.8 - (var3d_stag(:,:,k))/9.8  ! thickness [m]
 enddo
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     dz8w(i,k,j) = var3d(i,j,k) ! dz between full levels (m)     (in)
 enddo ; enddo ; enddo


! mean layer height 
 do k = kms, kme
    var3d(:,:,k) = 0.5 * ( (var3d_stag(:,:,k+1))/9.8 + (var3d_stag(:,:,k))/9.8 )  !height [m]
 enddo
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     z(i,k,j) = var3d(i,j,k) ! height above sea level (m)     (in)
 enddo ; enddo ; enddo

!-----vertical velocity
 call check( nf90_inq_varid(ncid, 'W' , varid ) )   !vertical wind [m/s] (staggered height)
 call check( nf90_get_var(ncid, varid, var3d_stag(is:ie,js:je,ks:ke+1), start=(/is,js,ks/), count=(/di,dj,dk+1/)))
 do k = kms, kme
    var3d(:,:,k) = 0.5 * ( (var3d_stag(:,:,k+1))/9.8 + (var3d_stag(:,:,k))/9.8 )  !mean vertical velocity
 enddo
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    w(i,k,j) = var3d(i,j,k) ! vertical velocity at half levels (m/s) (in)
 enddo ; enddo ; enddo


!------pressure and all sorts
 call check( nf90_inq_varid(ncid, 'PSFC' , varid ) ) !surface pressure [Pa]
 call check( nf90_get_var(ncid, varid, var2d(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))
 p8w(:,1,:) = var2d(:,:) ![Pa]
 call check( nf90_inq_varid(ncid, 'P_TOP' , varid ) ) !model top pressure [Pa]
 call check( nf90_get_var(ncid, varid, press_top ) )
 p8w(:,kme,:) = press_top
 call check( nf90_inq_varid(ncid, 'PB' , varid ) )   !Base pressure [Pa]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))

 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    p(i,k,j) = var3d(i,j,k)   ! base pressure (pa)  (in)
 enddo ; enddo ; enddo

!check lator
 do i = ims, ime ; do j = jms, jme ; do k = kms+1, kme
   p8w(i,k,j) = 0.5*var3d(i,j,k) + 0.5*var3d(i,j,k-1)  !ful level (interafce) <- half (layer) [Pa]
 enddo ; enddo ; enddo


!------temperature and all sorts
! call check( nf90_inq_varid(ncid, 'TSK' , varid ) )  ! surface skin temperature
! call check( nf90_get_var(ncid, varid, t8w_in(:,:,1) ) )
 t8w(:,1,:) = tsk(:,:)

 call check( nf90_inq_varid(ncid, 'T' , varid ) )    !purturvation potential temperature [K]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))

 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    t(i,k,j) = (var3d(i,j,k)+300.) / (100000./p(i,k,j))**0.286  !layer temperature [K]
 enddo ; enddo ; enddo

 do k = kms+1, kme
    t8w(:,k,:) = 0.5*t(:,k,:) + 0.5*t(:,k-1,:) ! Full level (interafce) <- half (layer) [K]
 enddo

 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     rho(i,k,j) = p(i,k,j) / (287.*t(i,k,j))  !dry air density [kg/m3]  (toshii: check this should be moist air density?)
     pi_phy(i,k,j)  = (p(i,k,j)/100000.)**(r_d/cp)  !exner function [-] 
 enddo ; enddo ; enddo

 th = t / pi_phy  ! potential temperature [K]




!-------- humidity and hydrometeors
 call check( nf90_inq_varid(ncid, 'QVAPOR' , varid ) ) !water vapor mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))

 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qv_curr  (i,k,j) = var3d(i,j,k) ! water vapor mixing ratio (kg/kg)  (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QCLOUD' , varid ) ) !cloud water mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qc_curr  (i,k,j) = var3d(i,j,k) ! cloud water mixing ratio (kg/kg)  (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QRAIN' , varid ) ) !rain mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qr_curr  (i,k,j) = var3d(i,j,k) ! rain water mixing ratio (kg/kg)   (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QICE' , varid ) ) !cloud ice mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qi_curr  (i,k,j) = var3d(i,j,k) ! cloud ice mixing ratio (kg/kg)    (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QSNOW' , varid ) ) !snow mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qs_curr  (i,k,j) = var3d  (i,j,k)  ! snow mixing ratio (kg/kg)         (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QGRAUP' , varid ) ) !graupel mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qg_curr  (i,k,j) = var3d  (i,j,k)  ! graupel mixing ration (kg/kg)     (in)
 enddo ; enddo ; enddo

! call check( nf90_inq_varid(ncid, 'QHAIL' , varid ) ) !hail mixing ratio[kg/kg]
! call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qh_curr  (i,k,j) = 0.  ! hail mixing ration (kg/kg)     (in)
 enddo ; enddo ; enddo


!
! derive effective radius here
!
 do i = ims, ime ; do j = jms, jme 
   do k = kms, kme
      q1d(k,id_cloud) = max(0.e0, qc_curr(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_ice  ) = max(0.e0, qi_curr(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_rain ) = max(0.e0, qr_curr(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_snow ) = max(0.e0, qs_curr(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_graupel)=max(0.e0, qg_curr(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_hail ) = max(0.e0, qh_curr(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      t1d(k) = t(i,k,j)  !air temp
   enddo

   call effective_radius( kme-kms+1, q1d, t1d, re1d )

   do k = kms, kme !convert into global parameters
      refc(i,k,j) = re1d(k,id_cloud  ) 
      refi(i,k,j) = re1d(k,id_ice    ) 
      refr(i,k,j) = re1d(k,id_rain   ) 
      refs(i,k,j) = re1d(k,id_snow   ) 
      refg(i,k,j) = re1d(k,id_graupel)
      refh(i,k,j) = re1d(k,id_hail   )
   enddo

   do k = kms, kme
      if( refi(i,k,j) > 0. ) qni_curr(i,k,j) = 1000.   !ice number conc. [#/kg] 
   enddo 

 enddo ; enddo
 
#if ( WRF_CHEM == 1)

 call check( nf90_inq_varid(ncid, 'sulf' , varid ) ) ! SO4 [ppmv]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,1) = var3d(i,j,k) *1.0e-6*(p(i,k,j))*96.0 &
                          /(8.314*t(i,k,j) ) ! SO4 [g/m3] 
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'BC1' , varid ) ) ! Black Carbon 1 [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,2) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'BC2' , varid ) ) ! Black Carbon 2 [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,2) = aero(i,k,j,2) + var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'OC1' , varid ) ) ! Organic Carbon 1 (non-hygroscopic) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,3) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6 * 1.4  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'OC2' , varid ) ) ! Organic Carbon 2 (hygroscopic) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,4) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6 * 1.4  ![g/m3]
 enddo ; enddo ; enddo


 call check( nf90_inq_varid(ncid, 'SEAS_1' , varid ) ) ! Sea Salt (fine mode: 0.3re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,5) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SEAS_2' , varid ) ) ! Sea Salt (coarse mode: 1.0re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,6) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SEAS_3' , varid ) ) ! Sea Salt (coarse mode: 3.2re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,6) = aero(i,k,j,6) + var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'SEAS_4' , varid ) ) ! Sea Salt (coarse mode: 7.5re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,6) = aero(i,k,j,6) + var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DUST_1' , varid ) ) ! Dust mode 1~4 (0.5re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,7)  = var3d(i,j,k) * rho(i,k,j) * 1.0e-6 * dust_frac(1) ![g/m3] 
    aero(i,k,j,8)  = var3d(i,j,k) * rho(i,k,j) * 1.0e-6 * dust_frac(2) ![g/m3] 
    aero(i,k,j,9)  = var3d(i,j,k) * rho(i,k,j) * 1.0e-6 * dust_frac(3) ![g/m3] 
    aero(i,k,j,10) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6 * dust_frac(4) ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DUST_2' , varid ) ) ! Dust mode 5 (1.4re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,11) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DUST_3' , varid ) ) ! Dust mode 6 (2.4re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,12) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DUST_4' , varid ) ) ! Dust mode 7 (4.5re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,13) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'DUST_5' , varid ) ) ! Dust mode 8 (8.0re) [ug/kg-dryair]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
    aero(i,k,j,14) = var3d(i,j,k) * rho(i,k,j) * 1.0e-6  ![g/m3]
 enddo ; enddo ; enddo


#endif

end subroutine read_netcdf


! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine effective_radius( kmax, q, tair, re  )
 implicit none
 integer,intent(in) :: kmax           !vertical layers
 real,intent(inout)  :: q(kmax,id_hmax)  !particle mass concentration [g/m3]
 real,intent(in)  :: tair(kmax)       !air temperature [K]
 real,intent(out) :: re(kmax,id_hmax) !particle effective radius [micron]

 integer :: id,k
!
! PSD parameters
!
 real :: n0_all (id_hmax)  ! intercept of gamma PSD [1/m4]
 real :: rho_all(id_hmax)  ! bulk particle density [kg/m3]
 real :: mu_all (id_hmax)  ! mu of gamma PSD [-]

 real,parameter :: re_cloud = 10.e0


!
! Gamma PSD parameter
! This is temporally soluion, but eventually, all effective radius must be pulled from microphyscs routine. 
!
 n0_all(id_cloud)   = HUGE(0.e0)  ! intercept for cloud [1/m4] (not used)
 n0_all(id_ice)     = HUGE(0.e0)   ! intercept for ice [1/m4] (not used)
 n0_all(id_rain)    = 0.08e+8   ! intercept for rain [1/m4] (=0.08[1/cm4])
 n0_all(id_snow)    = 0.16e+8   ! intercept for snow [1/m4] (=0.16[1/cm4]) !this is default (actually T-q function)
 n0_all(id_graupel) = 0.04e+8   ! intercept for graupel [1/m4] (=0.04[1/cm4]) !his is default (actually T-q function)
 n0_all(id_hail)    = 0.002e+8  ! intercept for hail [1/m4] (=0.002[1/cm4]) 

 rho_all(id_cloud)   = HUGE(0.e0)   ! density for cloud [1/m4] (=0.08[1/cm4])
 rho_all(id_ice)     = HUGE(0.e0)   ! density for ice [1/m4] (=0.08[1/cm4])
 rho_all(id_rain)    = 1.0e+3  ! density of rain [kg/m3] (=1.0[g/cm3])
 rho_all(id_snow)    = 0.05e+3 ! density of snow [kg/m3] (=0.05[g/cm3]) (Note it's reduced from gsfc 2003 scheme.)
 rho_all(id_graupel) = 0.2e+3  ! density of graupel [kg/m3] (=0.2[g/cm3]) (Note it's reduced from gsfc 2003 scheme.)
 rho_all(id_hail)    = 0.9e+3  ! density of hail [kg/m3] (=0.9[g/cm3])

 mu_all(id_cloud)    = 0.  !
 mu_all(id_ice)      = 0.  !
 mu_all(id_rain)     = 0.  !
 mu_all(id_snow)     = 0.  !
 mu_all(id_graupel)  = 0.  !
 mu_all(id_hail)     = 0.  !


 do k = 1, kmax

    !
    ! Goddard Microophysics 3ICE 2010 Version
    !
    re(k,id_cloud) = re_cloud ! mono distribution (must be modified)
    re(k,id_ice) = max( min( 125. + (tair(k)-243.16e0)*5.e0  ,150.e0), 25.e0)  ! ice effective radius

    call re_gamma_1mmt(rho_all(id_rain), n0_all(id_rain)  , q(k,id_rain) , re(k,id_rain) )

    call re_bulk_god10_snow ( tair(k), rho_all(id_snow)   , q(k,id_snow), re(k,id_snow) )

    call re_bulk_god10_graupel( tair(k), rho_all(id_graupel), q(k,id_graupel), re(k,id_graupel) )

    call re_gamma_1mmt(rho_all(id_hail), n0_all(id_hail)  , q(k,id_hail) , re(k,id_hail) )

    !
    ! Goddard Microophysics 4ICE 2014 Version
    !


 enddo

!
!
!
 do id = 1, id_hmax ; do k = 1, kmax

    if( re(k,id) == 0.e0 ) q(k,id) = 0.e0

 enddo  ; enddo


 return
 end subroutine effective_radius

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine re_gamma_1mmt(den, n0 , wc, re )
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! Compute drop effective radius for one-moment bulk scheme using analytic solution (fast). 
!  
! History:
! 08/2007  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 real,intent(in) :: den  !density [kg/m3]
 real,intent(in) :: n0   !intercept of exponential DSD [1/m4]
 real,intent(inout) :: wc   !liquid/ice water content [g/m3] 
 real,intent(out) :: re  !effective radius [micron] 
 real :: lam   !intercept [1/m]

!
! for no particles. 
!
 
 if( wc <= 1.e-5) then
    re = 0.e0 ; wc = 0.e0
 else

!
! compute drop effective radius for exponential distribution N(D) = N0*exp(-lam*D)
!
 lam = (n0*const_pi*den/(max( wc*1.e-3 ,tiny(1.0))   ))**(0.25e0)  ![1/m]
 re = 3.e0 / 2.e0 / lam *1.0e+6  ![micron] <- [m]

 endif

 return
 end subroutine re_gamma_1mmt

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine re_bulk_god10_snow(tair, rho, q_snow, re)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!   Compute slope of exponential PSD as a function of temperature and mixing ratio for snow.
!  
! History:
! 03/2009  Toshi Matsui@NASA GSFC ; Adapted from GCE. 
!           
! References: 
!
!----------------------------------------------------------------------------------------------------
 real,intent(in) :: tair    ! air temperature [k]
 real,intent(in) :: rho     ! density [kg/m3]
 real,intent(in) :: q_snow  ! snow mass mixing ratio [g/m3]
 real,intent(out) :: re     ! effective radius [micron]

!local parameters
 real :: slope  ! slope of exponential PSD [1/cm]
 real :: roqs
 real :: tslopes, cpi, dsnomin
 real :: xs, sno11, sno00, dsno11, dsno00, sexp11, sexp00
 real :: slim11, slim00, slit1, slit0, slexp, stt, stexp, sbase
 real :: a2, tairc, tairc1 ,ftnsT, sno1, dsno1, sexp1, slim1, tairc11
 real :: sdl,  xx
 real :: fexp, fftnsT, ftnsQ, fftnsTQ, tnsmax
 real :: lambda   ! lambda [1/cm]

 real,parameter :: tns = 0.10e0  ! default snow intercept [cm-4]
!
! too small mixing ratio -> return
!
 if( q_snow == 0.e0 ) then
    re = 0.e0 ; return
 endif

 roqs= rho * 1.e-3  !Snow density [g/cm3] <- [kg/m3]

! roqs=.05    !SNOW DENSITY (G/CM^3)

 tslopes=0.03218876   ! increase tns by  5 from  0 to -50C
 cpi=4.*atan(1.)
 dsnomin=0.0185 !minimum snow diameter (cm)
 dsnomin=dsnomin**4.*roqs*cpi


!
! tunable parameters
!
 xs=0.97
 sno11=0.40
 sno00=0.10
 dsno11=0.95
 dsno00=0.50
 sexp11=0.9
 sexp00=0.9
 slim11=1.0
 slim00=1.0
 slit1=-60.
 slit0=-20.
 slexp=3.3
 stt=-20.
 stexp=0.5
 sbase=0.00167


 a2= q_snow  ! [g/m3]

 fftnstq=1.0

 tairc=min(0.,max(-50.,tair-273.16)+0.0)
 tairc1=tair-273.16

 ftnsT=exp(-1.*tslopes*tairc)
 sno1=sno11
 dsno1=dsno11
 sexp1=sexp11
 slim1=slim11
 if(tairc1.gt.slit1.and.tairc1.lt.slit0)then
    tairc11=tairc1-slit0
    sdl=slit1-slit0
    slim1=slim00-(slim00-slim11)*(tairc11/sdl)**slexp
 endif
 if(tairc1.ge.slit0) slim1=slim00
 if(tairc.gt.stt)then
     sno1=sno00-(sno00-sno11)*(tairc/(stt))**stexp
     dsno1=dsno00-(dsno00-dsno11)*(tairc/(stt))**stexp
     sexp1=sexp00-(sexp00-sexp11)*(tairc/(stt))**stexp
 endif
 xx=xs-xs*min(slim1,max(0.0,(a2-sno1)/dsno1)**sexp1)

 ftnsT=ftnsT**xx
 fftnsT=ftnsT*tns
 fexp=xx
 ftnsQ=1.0
 if(a2.gt.0.) ftnsQ=(a2/sbase)**fexp
     fftnsTQ=fftnsT*ftnsQ
     tnsmax=(q_snow*1.e-6)/dsnomin   ! intercept [1/cm4]
!             [g/cm3]     / [g cm]

 if(fftnsTQ.gt.tnsmax.and.a2.gt.1.e-20) fftnsTQ=tnsmax  ! intercept [1/cm4]

 lambda = ( fftnsTQ * const_pi * roqs / (q_snow*1.e-6) ) ** 0.25e0  ! lambda [1/cm]

 re = 3.e0 / 2.e0 / lambda * 1.e+4   !effective radius [micron]

 return
 end subroutine re_bulk_god10_snow

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine re_bulk_god10_graupel(tair, rho, q_graupel,re)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!   Compute slope of exponential PSD as a function of temperature and mixing ratio for graupel.
!  
! History:
! 03/2009  Toshi Matsui@NASA GSFC ; Adapted from GCE. 
!           
! References: 
!
!----------------------------------------------------------------------------------------------------
 real,intent(in) :: tair      ! air temperature [k]
 real,intent(in) :: rho       ! density [kg/m3]
 real,intent(in) :: q_graupel ! graupel mass mixing ratio [g/m3]
 real,intent(out) :: re       ! effective radius [micron]

!local parameters
 real :: slope    ! slope of exponential PSD [1/cm]
 real :: roqg
 real :: tslopeg, cpi, dgrpmin, xg, grp11, grp00, dgrp11, dgrp00, gexp11, gexp00
 real :: gtt, gtexp, gbase, a3, fftngtq
 real :: tairc, tairc1, ftngT, grp1, dgrp1, gexp1, xx
 real :: fexp, fftngT, ftngQ, tngmax
 real,parameter :: tng = 0.04e0  ! default snow intercept [cm-4]
 real :: lambda   ! lambda [1/cm]


!
! too small mixing ratio -> return
!
 if( q_graupel == 0.e0 ) then
    re = 0.e0 ; return
 endif

 roqg= rho * 1.e-3  !Graupel density [g/cm3] <- [kg/m3]
! roqg=.2     !graupel DENSITY (G/CM^3)

 tslopeg=0.02197225   ! increase tng by  3 from  0 to -50C

 cpi=4.*atan(1.)
 dgrpmin=0.0225 !minimum graupel diameter (cm)
 dgrpmin=dgrpmin**4.*roqg*cpi

!
! tunable parameters
!
 xg=0.93
 grp11=0.55
 grp00=0.45
 dgrp11=3.90
 dgrp00=0.25
 gexp11=0.5
 gexp00=0.7
 gtt=-20.
 gtexp=0.50
 gbase=0.0095

 a3=q_graupel  ! [g/m3]
 fftngtq=1.0

 tairc=min(0.,max(-50.,tair-273.16)+0.0)
 tairc1=tair-273.16

 tairc=min(0.,max(-50.,tair-273.16)+0.0)
 ftngT=exp(-1.*tslopeg*tairc)
 grp1=grp11
 dgrp1=dgrp11
 gexp1=gexp11
 if(tairc.gt.gtt)then
    grp1=grp00-(grp00-grp11)*(tairc/(gtt))**gtexp
    dgrp1=dgrp00-(dgrp00-dgrp11)*(tairc/(gtt))**gtexp
    gexp1=gexp00-(gexp00-gexp11)*(tairc/(gtt))**gtexp
 endif
 xx=xg-xg*min(1.,max(0.0,(a3-grp1)/dgrp1)**gexp1)
 ftngT=ftngT**xx
 fftngT=ftngT*tng
 fexp=xx
 ftngQ=1.0
 if(a3.gt.0.) ftngQ=(a3/gbase)**fexp
 fftngTQ=fftngT*ftngQ

 tngmax=(q_graupel*1.e-6)/dgrpmin  ! intercept [1/cm4]

 if(fftngTQ.gt.tngmax.and.a3.gt.1.e-20) fftngTQ=tngmax ! intercept [1/cm4]

 lambda = ( fftngTQ * const_pi * roqg / (q_graupel*1.e-6) ) ** 0.25e0  ! lambda [1/cm]

 re = 3.e0 / 2.e0 / lambda * 1.e+4   !effective radius [micron]

 return
 end subroutine re_bulk_god10_graupel

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine write_out_gmp(proc)
 implicit none
 character(len=*) :: proc

 integer :: irec
 integer :: i,j,k,n
 real,parameter :: gridsize=6.0
 character*100 :: outfile
 integer,parameter :: io=11

!
! GrADS binary output file for eneryg budget and cloud/aerosol optical depth
!
  outfile = trim(outdir)//trim(CRM_file_list)//trim(out_appendix)//'_'//trim(proc)//'.bin'

  if(masterproc) print*,'MSG write output -->',outfile

  if(masterproc) open(unit= io, file = trim(outfile), access='direct', status='replace' &
                ,recl=(ide-ids+1)*(jde-jds+1)*(kde-kds+1)*4 )
  irec = 1 


 call dump( z(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( p(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( t(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( w(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( qv_curr(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( qc_curr(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( qr_curr(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( qi_curr(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( qs_curr(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( qg_curr(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( qh_curr(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( refc(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( refr(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( refi(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( refs(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( refg(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( refh(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( phys_tot(ims:ime,kms:kme,jms:jme), io, irec )
 call dump( refl_10cm(ims:ime,kms:kme,jms:jme), io, irec )


  if(masterproc) close(io)

!
!  write out grads control file
!
 if(masterproc) then

!
! GrADS Control files
!
  outfile = trim(outdir)//trim(CRM_file_list)//trim(out_appendix)//'_'//trim(proc)//'.ctl'
  if(masterproc) print*, trim(outfile)
  open(io,file=trim(outfile))

  outfile = trim(CRM_file_list)//trim(out_appendix)//'_'//trim(proc)//'.bin'
  write(io,"(a)") 'dset  '//trim(outfile)
  write(io,*) 'title  Energy budgets '
  write(io,*) 'undef  ',-999.0
  write(io,*) 'xdef   ',(ide-ids+1),' linear 1. ',gridsize
  write(io,*) 'ydef   ',(jde-jds+1),' linear 1. ',gridsize
  write(io,*) 'zdef   ',(kde-kds+1),' linear 0. 1. '
  write(io,*) 'tdef    1 linear 00Z20jan2007  1hr'
  write(io,*) 'vars ',irec-1
  write(io,*) 'z    ',(kde-kds+1),' 0  height [m]'
  write(io,*) 'p    ',(kde-kds+1),' 0  pressure [hPa]'
  write(io,*) 't    ',(kde-kds+1),' 0  temperature [K]'
  write(io,*) 'w    ',(kde-kds+1),' 0  vertical velocity [m/s]'
  write(io,*) 'qv   ',(kde-kds+1),' 0  water vapor [kg/kg]'
  write(io,*) 'qc   ',(kde-kds+1),' 0  cloud [kg/kg]'
  write(io,*) 'qr   ',(kde-kds+1),' 0  rain [kg/kg]'
  write(io,*) 'qi   ',(kde-kds+1),' 0  cloud ice [kg/kg]'
  write(io,*) 'qs   ',(kde-kds+1),' 0  snow [kg/kg]'
  write(io,*) 'qg   ',(kde-kds+1),' 0  graupel [kg/kg]'
  write(io,*) 'qh   ',(kde-kds+1),' 0  hail [kg/kg]'
  write(io,*) 'refc ',(kde-kds+1),' 0  cloud re [um]'
  write(io,*) 'refr ',(kde-kds+1),' 0  rain re [um]'
  write(io,*) 'refi ',(kde-kds+1),' 0  ice re [um]'
  write(io,*) 'refs ',(kde-kds+1),' 0  snow re [um]'
  write(io,*) 'refg ',(kde-kds+1),' 0  graupel re [um]'
  write(io,*) 'refh ',(kde-kds+1),' 0  hail re [um]'
  write(io,*) 'lh   ',(kde-kds+1),' 0  total latent heat [K/s]'
  write(io,*) 'dbz  ',(kde-kds+1),' 0  10cm reflectivity [dBZ]'

 write(io,*) 'endvars'
 close(io)


 endif


 end subroutine write_out_gmp

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine dump2d( var2d, io, irec)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out 2D real array in GrADS (binary) format. 
! If size of var2d is smaller than domain size, var2d will be gathered into the domain2d.
!
! History:
! 11/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------

 real(sdsu_fps),dimension(ims:ime, jms:jme), intent(in)  :: var2d  ! 
 integer,intent(in)    :: io   ! io index
 integer,intent(inout) :: irec ! binary record #
 integer :: i,j


#if MPI == 2

  out_domain2d = -999.0

  do j=jms,jme ; do i=ims,ime
     out_domain2d(i,j) = var2d(i,j)
  enddo ; enddo

  call mpi_sdsu_communicate( 'TO_MASTER', out_domain2d )

 if(masterproc) write(io,rec=irec) ((out_domain2d(i,j),i=ids,ide),j=jds,jde)
 irec=irec+1

#else

  write(io,rec=irec) ((var2d(i,j),i=ids,ide),j=jds,jde)
  irec=irec+1

#endif


 return
 end subroutine dump2d

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine dump3d( var3d_ikj, io, irec)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Write out 3D real array in GrADS (binary) format. 
! If size of var3d is smaller than domain size, var3d will be gathered into the domain3d via mpi library. 
! Vertical level can be either CRM (1~mxlyr) or satellite (1~variable) levels. 
!
! History:
! 11/2009  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real(sdsu_fps),dimension(ims:ime,kms:kme,jms:jme), intent(in)  :: var3d_ikj  ! assumed 3D array input 
                                                                           !  (specify starting index)
 integer,intent(in)    :: io   ! io index
 integer,intent(inout) :: irec ! binary record #
 integer :: i,j,k
 integer :: bnd(3)      ! upper memory bound for 3D array
 integer :: ks, ke      ! starting and ending k index
 real(sdsu_fps),dimension(ims:ime,jms:jme,kms:kme) :: var3d_ijk  ! assumed 3D array input 

!
!
!
 do k = kms, kme ; do j = jms, jme ; do i = ims, ime
    var3d_ijk(i,j,k) = var3d_ikj(i,k,j)
 enddo ; enddo ; enddo
#if MPI < 2  

!    
! In this case var3d can be CRM or Satellite vertical levels 
!
    write(io,rec=irec) (((var3d_ijk(i,j,k),i=ids,ide),j=jds,jde),k=kds,kde)


#else


   !
   ! gather slave memory into master memory
   !
   if( .not. masterproc) then

   call mpi_sdsu_collect_tile( ke,    var3d_ijk ( ims:ime, jms:jme, kms:kme ),&
                               out_domain3d ( ids:ide, jds:jde, kds:kde )  )
   
   else
      out_domain3d ( ims:ime, jms:jme, kms:kme ) = var3d_ijk ( ims:ime, jms:jme, kms:kme )
   endif
     

   !
   ! write CRM-dimension output
   !
   if(masterproc) write(io,rec=irec) (((out_domain3d(i,j,k),i=ids,ide),j=jds,jde),k=kds,kde)

#endif

 !
 ! add record # for next output
 !
 irec=irec+1

 return
 end subroutine dump3d

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

end module module_wrf_offline

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
