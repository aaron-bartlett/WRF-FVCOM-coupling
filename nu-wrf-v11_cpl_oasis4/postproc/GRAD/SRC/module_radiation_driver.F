#include <define.h>
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
 module module_radiation_driver
 use netcdf
 use module_runtime_param 
 use module_ra_goddard_2014
 use module_ra_goddard_2011
 use module_mpi
! use module_omp_check

#ifdef WRF_CHEM
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
     p,    & ! pressure (pa)                  (in)
     p8w,  & ! pressure at full levels (pa)   (in)
     pi,   & ! exner function (-)             (in)
     t,    & ! temperature (K)                (in)
     t8w,  & ! temperature at full levels (k) (in)
     rho,  & ! density (kg/m^3)               (in) 
     qv,   & ! water vapor mixing ratio (kg/kg)  (in)
     qc,   & ! cloud water mixing ratio (kg/kg)  (in)
     qr,   & ! rain water mixing ratio (kg/kg)   (in)
     qi,   & ! cloud ice mixing ratio (kg/kg)    (in)
     qs,   & ! snow mixing ratio (kg/kg)         (in)
     qg,   & ! graupel mixing ration (kg/kg)     (in)
     qh      ! hail mixing ratio (kg/kg)     (in)

 real,allocatable, dimension(:,:,:) ::  &
    refc,   & ! cloud effective radius [um]  (in)
    refr,   & ! rain effective radius [um]   (in)
    refi,   & ! cloud ice effective radius [um]    (in)
    refs,   & ! snow effective radius [um]  (in)
    refg,   & ! graupel effective radius [um]  (in)
    refh      ! hail effective radius [um] (in)


!
! WRF 1D variables
!
   integer :: julday    ! julian day (in)
   real :: julian ! julian day in real                             (in)
   real :: radt   !time for calling radiation (min)                (in)
   real :: gmt    ! greenwich mean time hour of model start (hour) (in)
   real :: xtime  ! time since simulation start (min)              (in)


!
! WRF 3D output variables
!
 real, allocatable, dimension(:,:,:) :: & 
     rthraten,   & ! theta tendency due to radiation (k/s)     (inout) 
     rthratenlw, & ! theta tendency due to lw radiation (k/s)  (inout)
     rthratensw, & ! theta tendency due to sw radiation (k/s)  (inout)
     cldfra,     & ! cloud fraction (0-1)            (inout) 
     taucldi,    & ! ice cloud optical thickness     (inout)
     taucldc       ! liquid cloud optical thickness  (inout)


!
! WRF 3D variables
!
 real, allocatable, dimension( :, :, : ) :: sflxd    !SW only
       ! 1-surface downward UV+VIS beam radiation [W/m2]
       ! 2-surface downward UV+VIS diffuse radiation [W/m2]
       ! 3-surface downward NIR beam radiation [W/m2]
       ! 4-surface downward NIR diffuse radiation [W/m2]

! 
! erbe_out: Energy Budget (last dimension 1-TOA LW down, 2-TOA LW up, 3-surface LW down, 4-surface LW up)
!                                         5-TOA SW down, 6-TOA SW up, 7-surface SW down, 8-surface SW up)
!
  real,allocatable :: erbe_out(:,:,:)  ! dimension(ims:ime, jms:jme, 1:8) !earth radiation budget (W/m2)

!
! WRF 2D variables
!
   real, allocatable, dimension( :,: ) :: &
      xlat,   & ! latitude, south is negative (degree)  (in)-
      xlong,  & ! longitude, west is negative (degree)  (in)-
      albedo, & ! broadband albedo (between 0 and 1)    (in)-
      emiss,  & ! surface emissivity (between 0 and 1)  (in)-
      tsk  ,  & ! surface skin temperature [K] (in) -
      cod2d_out , & ! cloud optical depth [-]
      ctop2d_out, & ! cloud-top pressure [mb]
      swddir  , & !
      swddni  , & !
      swddif  , & ! 
      coszen  , & !
      gsw,    &   !net short wave flux at ground surface (w/m^2)     (inout)
      glw,    &   !downward long wave flux at ground surface (w/m^2) (inout)
      swdown      !downward short wave flux at ground surface (w/m^2)  (local)


 real,allocatable :: var2d(:,:)   ! surface pressure at interface [Pa] 
 real,allocatable, dimension(:,:,:) ::  &
     var3d     , & ! temporally 3D variables
     var3d_stag    ! temp 3D stagger


#ifdef WRF_CHEM
  integer, parameter:: num_go = 14         ! # of gocart aerosol species (tgmx=14 defined in module_gocart_offline.F)
  real, allocatable :: aero( :, :, :, :)   !aerosol mass conc [g/m3]            
  real, allocatable :: aod2d_out(:,:)      ! column Aeorosol Optical Depth (-)
  real, allocatable :: atop2d_out(:,:)       ! aerosol top pressure [mb]

  real, allocatable :: aod3d_out( :, :, :)      !3D Aeorosol Optical Depth (-)

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
    jms=myj_start        ! start index for j in memory
    jme=myj_end        ! end index for j in memory
    kms=kds        ! start index for k in memory
    kme=kde        ! end index for k in memory

    kts=kds        ! start index for k in tile
    kte=kde-1      ! end index for k in tile



!
! allocation
!
  allocate( &
      out_domain2d( ids:ide, jds:jde ),& !2D domain output
      out_domain3d( ids:ide, jds:jde, kds:kde ),& !3D domain output
      out_memory_id( ids:ide, jds:jde ),& !
      out_tile_id  ( ids:ide, jds:jde ),& !
      sflxd   ( ims:ime, jms:jme, 1:4 ),&  !SW only
      erbe_out( ims:ime, jms:jme, 1:8) ,&  !earth radiation budget (W/m2)
      xlat    ( ims:ime, jms:jme ),   & ! latitude, south is negative (degree)  (in)-
      xlong   ( ims:ime, jms:jme ),   & ! longitude, west is negative (degree)  (in)-
      albedo  ( ims:ime, jms:jme ),   & ! broadband albedo (between 0 and 1)    (in)-
      emiss   ( ims:ime, jms:jme ),   & ! surface emissivity (between 0 and 1)  (in)-
      tsk     ( ims:ime, jms:jme ),   & ! surface skin temperature [K] (in) -
      cod2d_out ( ims:ime, jms:jme ),   & ! column cloud optical depth
      ctop2d_out( ims:ime, jms:jme ),   & ! cloud top pressure [mb]
      gsw     ( ims:ime, jms:jme ),   & ! net short wave flux at ground surface (w/m^2)     (inout)
      glw     ( ims:ime, jms:jme ),   & ! downward long wave flux at ground surface (w/m^2) (inout)
      swdown  ( ims:ime, jms:jme ),   &      !downward short wave flux at ground surface (w/m^2)  (local)
      swddir  ( ims:ime, jms:jme ),   &      !
      swddni  ( ims:ime, jms:jme ),   &      !
      swddif  ( ims:ime, jms:jme ),   &      ! 
      coszen  ( ims:ime, jms:jme ),   &      !
      rthraten   ( ims:ime, kms:kme, jms:jme ),   & ! theta tendency due to radiation (k/s)     (inout) 
      rthratenlw ( ims:ime, kms:kme, jms:jme ),   & ! theta tendency due to lw radiation (k/s)  (inout)
      rthratensw ( ims:ime, kms:kme, jms:jme ),   & ! theta tendency due to sw radiation (k/s)  (inout)
      cldfra     ( ims:ime, kms:kme, jms:jme ),   & ! cloud fraction (0-1)            (inout) 
      taucldi    ( ims:ime, kms:kme, jms:jme ),   & ! ice cloud optical thickness     (inout)
      taucldc    ( ims:ime, kms:kme, jms:jme ),   & ! liquid cloud optical thickness  (inout)
      dz8w       (ims:ime,kms:kme,jms:jme), & ! dz between full levels (m)     (in)
      p          (ims:ime,kms:kme,jms:jme), & ! pressure (pa)                  (in)
      p8w        (ims:ime,kms:kme,jms:jme), & ! pressure at full levels (pa)   (in)
      pi         (ims:ime,kms:kme,jms:jme), & ! exner function (-)             (in)
      t          (ims:ime,kms:kme,jms:jme), & ! temperature (K)                (in)
      t8w        (ims:ime,kms:kme,jms:jme), & ! temperature at full levels (k) (in)
      rho        (ims:ime,kms:kme,jms:jme), & ! density (kg/m^3)               (in) 
      qv         (ims:ime,kms:kme,jms:jme), & ! water vapor mixing ratio (kg/kg)  (in)
      qc         (ims:ime,kms:kme,jms:jme), & ! cloud water mixing ratio (kg/kg)  (in)
      qr         (ims:ime,kms:kme,jms:jme), & ! rain water mixing ratio (kg/kg)   (in)
      qi         (ims:ime,kms:kme,jms:jme), & ! cloud ice mixing ratio (kg/kg)    (in)
      qs         (ims:ime,kms:kme,jms:jme), & ! snow mixing ratio (kg/kg)         (in)
      qg         (ims:ime,kms:kme,jms:jme), & ! graupel mixing ration (kg/kg)     (in)
      qh         (ims:ime,kms:kme,jms:jme), & ! hail mixing ration (kg/kg)     (in)
     refc        (ims:ime,kms:kme,jms:jme), & ! cloud effective radius [um] (in)
     refr        (ims:ime,kms:kme,jms:jme), & ! rain effective radius [um]   (in)
     refi        (ims:ime,kms:kme,jms:jme), & ! cloud effective radius [um]    (in)
     refs        (ims:ime,kms:kme,jms:jme), & ! snow effective radius [um]  (in)
     refg        (ims:ime,kms:kme,jms:jme), & ! graupel effective radius [um]   (in)
     refh        (ims:ime,kms:kme,jms:jme), & ! hail effective radius [um]    (in)
      var3d_stag (ims:ime,jms:jme,kms-1:kme), & ! base geopotential
      var3d      (ims:ime,jms:jme,kms:kme),  & ! temporally 3D variables
      var2d      (ims:ime,jms:jme),    & ! surface pressure [Pa]
#ifdef WRF_CHEM
      aero( ims:ime, kms:kme, jms:jme, 1:num_go) , &  ! aerosol mass conc [g/m3]            
      aod3d_out(ims:ime, kms:kme, jms:jme)       , &  ! 3D Aeorosol Optical Depth (-)
      atop2d_out( ims:ime, jms:jme ),   & ! aerosol top pressure [mb]
      aod2d_out(ims:ime, jms:jme)                  , &  ! column Aeorosol Optical Depth (-)
#endif
  stat=ierr )

 if(ierr /= 0) stop 'MSG allocate_all: cannot allocate the variable'


 return
 end subroutine allocate_all

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

 subroutine radiation_driver
 implicit none
! iroutine: radiation_driver - interface to radiation physics options

!----------------------------------------------------------------------
! local  var


 real    :: declin ! solar declination angle (rad) (in)
 real    :: solcon ! solar constant (w/m^2) (in)
 integer :: i,j,k,ij !its,ite,jts,jte,ij

 integer,save ::  num_tiles    ! number of tiles

 integer, allocatable, save,dimension(:) ::  &  !all (local)
    i_start ,& !start indices for i in tile (>ids)
    i_end   ,& !end indices for i in tile  (<ide)
    j_start ,& !start indices for j in tile (>jds)
    j_end      !end indices for j in tile (<jde)

 real, allocatable, save,dimension(:) :: &
   clocktime_sw, & !clocktime of SW scheme
   clocktime_lw    !clocktime of LW scheme


!
! logical flags
!
   logical,parameter :: warm_rain = .false. ! warm rain only            (in)
   logical,parameter :: f_qv = .true.   !flag for simulating water vapor (in)
   logical,parameter :: f_qc = .true.  !flag for simulating cloud water (in)
   logical,parameter :: f_qr = .true.  !flag for simulating rain        (in)
   logical,parameter :: f_qi = .true.  !flag for simulating cloud ice   (in)
   logical,parameter :: f_qs = .true.  !flag for simulating snow        (in)
   logical,parameter :: f_qg = .true.  !flag for simulating graupel     (in)
   logical,parameter :: f_qh = .true.  !flag for simulating graupel     (in)

   integer,parameter :: icloud =1      !0: dry simulation

! urban related variables are added to declaration
   real :: declin_urb  !urban
   real, dimension( ims:ime, jms:jme) :: cosz_urb2d  !urban
   real, dimension( ims:ime, jms:jme) :: omg_urb2d   !urban

! misc parameter for calculations
 real  :: fac,x,xt24,tloctm,hrang,xxlat
 integer :: itimestep

#ifdef WRF_CHEM
  integer,parameter :: chem_opt = 300 ! EMK
  integer,parameter :: gsfcrad_gocart_coupling = 1! EMK
#endif

#ifdef OMP
   !     Obtain thread number
  integer ::  id_thread 
#endif

 
!
! -------------------------- program start -------------------------------------
!


!
! set up OMP/MIC thread tile boundary indices
!
  if( .not. allocated(i_start) ) then

    num_tiles = num_tiles_omp

    allocate(  i_start(num_tiles)  ,& !start indices for i in tile (>ims)
               i_end  (num_tiles)  ,& !end indices for i in tile  (<ime)
               j_start(num_tiles)  ,& !start indices for j in tile (>jms)
               j_end  (num_tiles) )     !end indices for j in tile (<jme)

    do ij = 1, num_tiles
       i_start(ij) = omp_tile(ij)%is    !start indices for i in tile (>ids)
       i_end  (ij) = omp_tile(ij)%ie    !end indices for i in tile  (<ide)
       j_start(ij) = omp_tile(ij)%js  !start indices for j in tile (>jds)
       j_end  (ij) = omp_tile(ij)%je  !end indices for j in tile (<jde)
    enddo

    allocate( clocktime_sw( num_tiles ),&
              clocktime_lw( num_tiles ) )

  endif

!
! initialze time counter
!
   clocktime_sw = 0.0 ; clocktime_lw = 0.0


   itimestep = 1 ! always 1 for this driver. 

   if (lw_physics .eq. 0 .and. sw_physics .eq. 0)         return


   out_memory_id(:,:) = real(myrank)   !memory ID


!   call getinfo_omp  !check omp parameter

   !$omp parallel do   &
   !$omp private ( ij ,i,j,k,its,ite,jts,jte)


!#ifdef OMP
!   id_thread = omp_get_thread_num()  
!#endif


!
! toshi--> comeback here for thread tiling...
!
   do ij = 1 , num_tiles
     its = i_start(ij)
     ite = i_end(ij)
     jts = j_start(ij)
     jte = j_end(ij)

! initialize data

     do j=jts,jte
     do i=its,ite
        out_tile_id(i,j) = real(ij)
        gsw(i,j)=0.
        glw(i,j)=0.
        swdown(i,j)=0.
     enddo
     enddo

     do j=jts,jte
!     do k=kts,kte+1
     do i=its,ite
        rthraten(i,:,j)=0.
     enddo
!     enddo
     enddo


#ifdef WRF_CHEM
!------Read Mie Tables (read_mie)-----------------------------------------------
! Call this just once at the beggining of simulation.
! Mie tables will be read & saved in aerrosol_mod
!       call read_mie_lw
!       call write_mie_lw
!------Derive & Interporate GOCART aerosol data (init_gocart & aero_interp)------
! Call these subroutines whenever you need to update GOCART data.
! For init_gocart you need the time input, "yyyymm".
!       call init_gocart('200405')

#endif


!------- calculate constant for short wave radiation -------------

     call radconst(xtime,declin,solcon,julian,degrad,dpd) 
     declin_urb=declin  ! urban


    !
    ! estimate cosin of solar zenith angle
    !
    do j=jts,jte
    do i=its,ite
       xt24 = mod(xtime + radt * 0.5e0, 1440.e0)
       tloctm = gmt + xt24 / 60.e0 + xlong(i,j) / 15.e0
       hrang = 15.e0 * (tloctm - 12.e0) * degrad
       xxlat = xlat(i,j) * degrad
       coszen(i,j) = sin(xxlat) * sin(declin) + &
                     cos(xxlat) * cos(declin) * cos(hrang)  !cosine of solar zenith angle [-]
    enddo
    enddo

!----------cloud fraction-----------------------------------

      call cal_cldfra(cldfra,qc,qi,f_qc,f_qi,               &
                      ids,ide, jds,jde, kds,kde,            &
                      ims,ime, jms,jme, kms,kme,            &
                      its,ite, jts,jte, kts,kte             )


!--------------------LW schem selection------------------------------------

     call cpu_time(ts) 

     lwrad_select: select case(lw_physics)


        case (goddardlwscheme_2011)

             CALL wrf_debug (100, 'CALL NUWRF goddard longwave radiation scheme 2011')

             CALL goddardrad_2011(sw_or_lw='lw'                             &
                    ,rthraten=rthraten,gsf=glw,xlat=xlat,xlong=xlong   &
                    ,alb=albedo,t3d=t,p3d=p,p8w3d=p8w,pi3d=pi          &
                    ,dz8w=dz8w,rho_phy=rho,emiss=emiss,tsk=tsk         &
                    ,cldfra3d=cldfra                                   &
                    ,gmt=gmt,cp=cp,g=g,t8w=t8w                         &
                    ,julday=julday,xtime=xtime                         &
                    ,declin=declin,solcon=solcon                       &
                    ,radfrq=radt,degrad=degrad                         &
                    ,taucldi=taucldi,taucldc=taucldc                   &
                    ,warm_rain=warm_rain                               &
                    ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde &
                    ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme &
                    ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte &
                    ,qv=qv                                           &
                    ,qc=qc                                           &
                    ,qr=qr                                           &
                    ,qi=qi                                           &
                    ,qs=qs                                           &
                    ,qg=qg                                           &
                    ,f_qv=f_qv,f_qc=f_qc,f_qr=f_qr                     &
                    ,f_qi=f_qi,f_qs=f_qs,f_qg=f_qg ,f_qh=f_qh           &
                    ,erbe_out=erbe_out                                 & 
!JJS 20101021 vvvvv
#ifdef WRF_CHEM
                    ,aero=aero                                      &
                    ,CHEM_OPT=chem_opt                                 &
                    ,GSFCRAD_GOCART_COUPLING=gsfcrad_gocart_coupling   &
#endif
!JJS 20101021 ^^^^^
                                                                       )


        case (goddardlwscheme_2014)

             CALL wrf_debug (100, 'CALL NUWRF goddard longwave radiation scheme 2014')

             CALL goddardrad_2014(sw_or_lw='lw',dx=dx,iskip=goddardrad_2014_skip &
                    ,rthraten=rthraten,gsf=glw,xlat=xlat,xlong=xlong   &
                    ,alb=albedo,t3d=t,p3d=p,p8w3d=p8w,pi3d=pi          &
                    ,dz8w=dz8w,rho_phy=rho,emiss=emiss,tsk=tsk         &
                    ,cldfra3d=cldfra                                   &
                    ,gmt=gmt,cp=cp,g=g,t8w=t8w                         &
                    ,julday=julday,xtime=xtime                         &
                    ,declin=declin,solcon=solcon                       &
                    ,radfrq=radt,degrad=degrad                         &
                    ,warm_rain=warm_rain                               &
                    ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde &
                    ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme &
                    ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte &
                    ,qv=qv,qc=qc,qi=qi,qr=qr,qs=qs,qg=qg,qh=qh         &
                    ,f_qv=f_qv,f_qc=f_qc,f_qr=f_qr                     &
                    ,f_qi=f_qi,f_qs=f_qs,f_qg=f_qg ,f_qh=f_qh          &
                    ,rec3d=refc,rei3d=refi,rer3d=refr,res3d=refs       & !optional
                    ,reg3d=refg,reh3d=refh                             & !optional 
                    ,erbe_out=erbe_out                                 & 
#ifdef WRF_CHEM
                    ,AERO=aero                                         &
                    ,CHEM_OPT=chem_opt                                 &
                    ,GSFCRAD_GOCART_COUPLING=gsfcrad_gocart_coupling   &
#endif
                                                                       )

        case default
  
             if(masterproc)  print*, 'the longwave option does not exist: lw_physics = ', lw_physics
           
     end select lwrad_select    


     call cpu_time(te) 
         if(masterproc) then
            clocktime_lw(ij) = te-ts !sec
            !print*,'total CPU time for LW scheme is',te-ts
         endif


     if (lw_physics .gt. 0 ) then
        do j=jts,jte
        do k=kts,kte
        do i=its,ite
           rthratenlw(i,k,j)=rthraten(i,k,j)
        enddo
        enddo
        enddo


     endif


     call cpu_time(ts)

! ----------------SW scheme selection------------------------------------------
     swrad_select: select case(sw_physics)


        case (goddardlwscheme_2011)

             CALL wrf_debug(100, 'CALL NUWRF goddard shortwave radiation scheme 2011')

             CALL goddardrad_2011(sw_or_lw='sw'                             &
                    ,rthraten=rthraten,gsf=gsw,xlat=xlat,xlong=xlong   &
                    ,alb=albedo,t3d=t,p3d=p,p8w3d=p8w,pi3d=pi          &
                    ,dz8w=dz8w,rho_phy=rho,emiss=emiss,tsk=tsk         &
                    ,cldfra3d=cldfra                                   &
                    ,gmt=gmt,cp=cp,g=g,t8w=t8w                         &
                    ,julday=julday,xtime=xtime                         &
                    ,declin=declin,solcon=solcon                       &
                    ,radfrq=radt,degrad=degrad                         &
                    ,taucldi=taucldi,taucldc=taucldc                   &
                    ,warm_rain=warm_rain                               &
                    ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde &
                    ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme &
                    ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte &
                    ,qv=qv,qc=qc,qi=qi,qr=qr,qs=qs,qg=qg               &
                    ,f_qv=f_qv,f_qc=f_qc,f_qr=f_qr                     &
                    ,f_qi=f_qi,f_qs=f_qs,f_qg=f_qg                     &
                    ,erbe_out=erbe_out                                 &
                    ,sflxd=sflxd                                       & !optional (for LIS)
#ifdef WRF_CHEM
                    ,aod_out=aod3d_out                                   & ! optional gocart
                    ,aero=aero                                      &
                    ,CHEM_OPT=chem_opt                                 &
                    ,GSFCRAD_GOCART_COUPLING=gsfcrad_gocart_coupling   &
#endif
                                                                       )
   

        ! 
        ! for 2011 scheme only for comparison with 2014 scheme.
        !
        do j=jts,jte 
        do i=its,ite
           cod2d_out(i,j) =  sum(taucldc(i,:,j))+sum(taucldi(i,:,j))
        enddo
        enddo

        do j=jts,jte
        do i=its,ite
           do k=kte,kts,-1
              if( taucldc(i,k,j)+taucldi(i,k,j) > 0.01) then
                 ctop2d_out(i,j) = p8w(i,k,j)/100.
                 exit 
              endif 
           enddo
        enddo
        enddo


#ifdef WRF_CHEM
        do j=jts,jte 
        do i=its,ite
           aod2d_out(i,j) = sum(aod3d_out(i,:,j))
        enddo
        enddo

        do j=jts,jte 
        do i=its,ite
           atop2d_out(i,j) = -999.
           do k=kte,kts,-1
              if( aod3d_out(i,k,j) > 0.01) then
                 atop2d_out(i,j) = p8w(i,k,j)/100.
                 exit
              endif
           enddo
        enddo
        enddo
#endif
 

        case (goddardswscheme_2014)

             CALL wrf_debug(100, 'CALL NUWRF goddard shortwave radiation scheme 2014')

             CALL goddardrad_2014(sw_or_lw='sw',dx=dx,iskip=goddardrad_2014_skip &
                    ,rthraten=rthraten,gsf=gsw,xlat=xlat,xlong=xlong   &
                    ,alb=albedo,t3d=t,p3d=p,p8w3d=p8w,pi3d=pi          &
                    ,dz8w=dz8w,rho_phy=rho,emiss=emiss,tsk=tsk         &
                    ,cldfra3d=cldfra                                   &
                    ,gmt=gmt,cp=cp,g=g,t8w=t8w                         &
                    ,julday=julday,xtime=xtime                         &
                    ,declin=declin,solcon=solcon                       &
                    ,radfrq=radt,degrad=degrad                         &
                    ,warm_rain=warm_rain                               &
                    ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde &
                    ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme &
                    ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte &
                    ,qv=qv,qc=qc,qr=qr,qi=qi,qs=qs,qg=qg,qh=qh         & !optional
                    ,f_qv=f_qv,f_qc=f_qc,f_qr=f_qr                     & !optional
                    ,f_qi=f_qi,f_qs=f_qs,f_qg=f_qg ,f_qh=f_qh          & !optional
                    ,rec3d=refc,rei3d=refi,rer3d=refr,res3d=refs       & !optional
                    ,reg3d=refg,reh3d=refh                             & !optional 
                    ,erbe_out=erbe_out                                 & 
                    ,cod2d_out=cod2d_out,ctop2d_out=ctop2d_out         & !optional
                    ,sflxd=sflxd                                       & !optional
                    ,swddir=swddir,swddni=swddni,swddif=swddif         & !optional 
                    ,coszen=coszen                                     & !optional
#ifdef WRF_CHEM
                    ,aod2d_out=aod2d_out, atop2d_out=atop2d_out        & ! optional
                    ,AERO=aero                                         &
                    ,CHEM_OPT=chem_opt                                 &
                    ,GSFCRAD_GOCART_COUPLING=gsfcrad_gocart_coupling   &
#endif
                                                                       )


        case default

             if(masterproc)  print*, 'the shortwave option does not exist: sw_physics = ', sw_physics

     end select swrad_select    


     call cpu_time(te) !; if(masterproc)  print*,'total CPU time for SW scheme is',te-ts

         if(masterproc) then
            clocktime_sw(ij) = te-ts !sec
            !print*,'total CPU time for LW scheme is',te-ts
         endif

 
     if (sw_physics .gt. 0 ) then
        do j=jts,jte
        do k=kts,kte
        do i=its,ite
           rthratensw(i,k,j)=rthraten(i,k,j)-rthratenlw(i,k,j)
        enddo
        enddo
        enddo
     endif


   enddo !tile loop
   !$omp end parallel do


  if(masterproc)  then
      print*,'total CPU time for LW scheme is',sum(clocktime_lw(:))
      print*,'total CPU time for SW scheme is',sum(clocktime_sw(:))
  endif

   end subroutine radiation_driver

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

   subroutine radconst(xtime,declin,solcon,julian,                   &
                       degrad,dpd                                    )
   implicit none

! !arguments:
   real, intent(in   )      ::       degrad,dpd,xtime,julian
   real, intent(out  )      ::       declin,solcon
   real                     ::       obecl,sinob,sxlong,arg,  &
                                     decdeg,djul,rjul,eccfac
!
! !description:
! compute terms used in radiation physics 
!eop

! for short wave radiation

   declin=0.
   solcon=0.

!-----obecl : obliquity = 23.5 degree.
        
   obecl=23.5*degrad
   sinob=sin(obecl)
        
!-----calculate longitude of the sun from vernal equinox:
        
   if(julian.ge.80.)sxlong=dpd*(julian-80.)
   if(julian.lt.80.)sxlong=dpd*(julian+285.)
   sxlong=sxlong*degrad
   arg=sinob*sin(sxlong)
   declin=asin(arg)
   decdeg=declin/degrad
!----solar constant eccentricity factor (paltridge and platt 1976)
   djul=julian*360./365.
   rjul=djul*degrad
   eccfac=1.000110+0.034221*cos(rjul)+0.001280*sin(rjul)+0.000719*  &
          cos(2*rjul)+0.000077*sin(2*rjul)
   solcon=1370.*eccfac
   

   end subroutine radconst

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF

!bop
! !iroutine: cal_cldfra - compute cloud fraction
! !interface:
   subroutine cal_cldfra(cldfra,qc,qi,f_qc,f_qi,                     &
          ids,ide, jds,jde, kds,kde,                                 &
          ims,ime, jms,jme, kms,kme,                                 &
          its,ite, jts,jte, kts,kte                                  )
   implicit none
   integer,  intent(in   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

!
   real, dimension( ims:ime, kms:kme, jms:jme ), intent(out  ) ::    &
                                                             cldfra

   real, dimension( ims:ime, kms:kme, jms:jme ), intent(in   ) ::    &
                                                                 qi, &
                                                                 qc

   logical,intent(in) :: f_qc,f_qi

   integer:: i,j,k

! !description:
! compute cloud fraction from input ice and cloud water fields
! if provided.
!
! whether qi or qc is active or not is determined from the indices of
! the fields into the 4d scalar arrays in wrf. these indices are
! p_qi and p_qc, respectively, and they are passed in to the routine
! to enable testing to see if qi and qc represent active fields in
! the moisture 4d scalar array carried by wrf.
!
! if a field is active its index will have a value greater than or
! equal to param_first_scalar, which is also an input argument to
! this routine.
!eop
!---------------------------------------------------------------------

     if ( f_qi .and. f_qc ) then
        do j = jts,jte
        do k = kts,kte
        do i = its,ite
           if ( qc(i,k,j)+qi(i,k,j) .gt. thresh) then
              cldfra(i,k,j)=1.
           else
              cldfra(i,k,j)=0.
           endif
        enddo
        enddo
        enddo
     else if ( f_qc ) then
        do j = jts,jte
        do k = kts,kte
        do i = its,ite
           if ( qc(i,k,j) .gt. thresh) then
              cldfra(i,k,j)=1.
           else
              cldfra(i,k,j)=0.
           endif
        enddo
        enddo
        enddo
     else if ( f_qi ) then
        do j = jts,jte
        do k = kts,kte
        do i = its,ite
           if ( qi(i,k,j) .gt. thresh) then
              cldfra(i,k,j)=1.
           else
              cldfra(i,k,j)=0.
           endif
        enddo
        enddo
        enddo
     else
        do j = jts,jte
        do k = kts,kte
        do i = its,ite
           cldfra(i,k,j)=0.
        enddo
        enddo
        enddo
     endif

   end subroutine cal_cldfra

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!  1D variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 julday = 218
 julian = real(julday)
 radt = 10.
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

!----albedo
 call check( nf90_inq_varid(ncid, 'ALBEDO' , varid ) )  ! BROADBAND ALBEDO  
 call check( nf90_get_var(ncid, varid, albedo(is:ie,js:je), start=(/is,js/), count=(/di,dj/)))

!---- emiss
! call check( nf90_inq_varid(ncid, 'EMISS' , varid ) )  ! BROADBAND EMISSIVITY  
! call check( nf90_get_var(ncid, varid, emiss ) )
 emiss = 0.98

!!!!!!!!!  3D variable -> first read (i,j,k) then transform (i,k,j)  !!!!!!!!!!!!!!!!
!----height and all sorts
 call check( nf90_inq_varid(ncid, 'PHB' , varid ) ) ! geopotential [m/s2] 
 call check( nf90_get_var(ncid, varid, var3d_stag(is:ie,js:je,ks:ke+1), start=(/is,js,ks/), count=(/di,dj,dk+1/)))

 do k = kms, kme
    var3d(:,:,k) = (var3d_stag(:,:,k))/9.8 - (var3d_stag(:,:,k-1))/9.8  ![m]
 enddo
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     dz8w(i,k,j) = var3d(i,j,k) ! dz between full levels (m)     (in)
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
     rho(i,k,j) = p(i,k,j) / (287.*t(i,k,j))  !dry air density [kg/m3]
     pi(i,k,j)  = (p(i,k,j)/100000.)**(r_d/cp)  !exner function [-] 
 enddo ; enddo ; enddo


!-------- humidity and hydrometeors
 call check( nf90_inq_varid(ncid, 'QVAPOR' , varid ) ) !water vapor mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))

 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qv  (i,k,j) = var3d(i,j,k) ! water vapor mixing ratio (kg/kg)  (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QCLOUD' , varid ) ) !cloud water mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qc  (i,k,j) = var3d(i,j,k) ! cloud water mixing ratio (kg/kg)  (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QRAIN' , varid ) ) !rain mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qr  (i,k,j) = var3d(i,j,k) ! rain water mixing ratio (kg/kg)   (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QICE' , varid ) ) !cloud ice mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qi  (i,k,j) = var3d(i,j,k) ! cloud ice mixing ratio (kg/kg)    (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QSNOW' , varid ) ) !snow mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qs  (i,k,j) = var3d  (i,j,k)  ! snow mixing ratio (kg/kg)         (in)
 enddo ; enddo ; enddo

 call check( nf90_inq_varid(ncid, 'QGRAUP' , varid ) ) !graupel mixing ratio[kg/kg]
 call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qg  (i,k,j) = var3d  (i,j,k)  ! graupel mixing ration (kg/kg)     (in)
 enddo ; enddo ; enddo

! call check( nf90_inq_varid(ncid, 'QHAIL' , varid ) ) !hail mixing ratio[kg/kg]
! call check( nf90_get_var(ncid, varid, var3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
 do i = ims, ime ; do j = jms, jme ; do k = kms, kme
     qh  (i,k,j) = 0.  ! hail mixing ration (kg/kg)     (in)
 enddo ; enddo ; enddo


!
! derive effective radius here
!
 do i = ims, ime ; do j = jms, jme 
   do k = kms, kme
      q1d(k,id_cloud) = max(0.e0, qc(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_ice  ) = max(0.e0, qi(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_rain ) = max(0.e0, qr(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_snow ) = max(0.e0, qs(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_graupel)=max(0.e0, qg(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
      q1d(k,id_hail ) = max(0.e0, qh(i,k,j)) * rho(i,k,j) * 1.e+3  ! mass concentration [g/m3]
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

 enddo ; enddo
 
#ifdef WRF_CHEM

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

 subroutine write_out
 implicit none
 integer :: irec
 integer :: i,j,k,n
 real,parameter :: gridsize=6.0
 character*100 :: outfile
 integer,parameter :: io=11

!
! GrADS binary output file for eneryg budget and cloud/aerosol optical depth
!

  outfile = trim(outdir)//trim(CRM_file_list)//trim(out_appendix)//'.bin'
  if(masterproc) print*,'MSG write output -->',outfile

  if(masterproc) open(unit= io, file = trim(outfile), access='direct', status='replace' &
                ,recl=(ide-ids+1)*(jde-jds+1)*4 )
  irec = 1 

  do n = 1, 8
    call dump( erbe_out(ims:ime,jms:jme,n), io, irec ) ! energy budget [W/m2]
  enddo


  call dump( cod2d_out (ims:ime,jms:jme), io, irec )
  call dump( ctop2d_out(ims:ime,jms:jme), io, irec )


#ifdef WRF_CHEM
  call dump( aod2d_out(ims:ime,jms:jme), io, irec )
  call dump( atop2d_out(ims:ime,jms:jme), io, irec )
#endif

  call dump( out_memory_id(ims:ime,jms:jme), io, irec )
  call dump( out_tile_id  (ims:ime,jms:jme), io, irec )


  if(masterproc) close(io)

!
!  write out grads control file
!
 if(masterproc) then

!
! GrADS Control files
!
  outfile = trim(outdir)//trim(CRM_file_list)//trim(out_appendix)//'.ctl'
  if(masterproc) print*, trim(outfile)
  open(io,file=trim(outfile))

  outfile = trim(CRM_file_list)//trim(out_appendix)//'.bin'
  write(io,"(a)") 'dset  '//trim(outfile)
  write(io,*) 'title  Energy budgets '
  write(io,*) 'undef  ',-999.0
  write(io,*) 'xdef   ',(ide-ids+1),' linear 1. ',gridsize
  write(io,*) 'ydef   ',(jde-jds+1),' linear 1. ',gridsize
  write(io,*) 'zdef    1 linear 0. 1. '
  write(io,*) 'tdef    1 linear 00Z20jan2007  1hr'
  write(io,*) 'vars ',irec-1
    write(io,*) 'tld   0 0  TOA LW down [W/m2]'
    write(io,*) 'tlu   0 0  TOA LW up   [W/m2]'
    write(io,*) 'bld   0 0  BOA LW down [W/m2]'
    write(io,*) 'blu   0 0  BOA LW up   [W/m2]'
    write(io,*) 'tsd   0 0  TOA SW down [W/m2]'
    write(io,*) 'tsu   0 0  TOA SW up   [W/m2]'
    write(io,*) 'bsd   0 0  BOA SW down [W/m2]'
    write(io,*) 'bsu   0 0  BOA SW up [W/m2]'
    write(io,*) 'cod   0 0  cloud optical depth (vis) [-]'
    write(io,*) 'ctop  0 0  cloud-top pressure [mb]'
#ifdef WRF_CHEM
    write(io,*) 'aod   0 0  aerosol optical depth (vis) [-]'
    write(io,*) 'atop   0 0  aerosol-top pressure [mb]'
#endif
    write(io,*) 'memory  0 0  memory (vis) [-]'
    write(io,*) 'tile    0 0  tile [mb]'
 write(io,*) 'endvars'
 close(io)


!
! Atmospheric Heating Rate
!
!  write(1,rec=1) (((rthratensw(i,k,j),i=ids,ide),j=jds,jde),k=kds,kde) !SW heating [K/sec] 
!  write(1,rec=2) (((rthratenlw(i,k,j),i=ids,ide),j=jds,jde),k=kds,kde) !LW cooling [K/sec]
!  close(1)

 endif


 end subroutine write_out

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

 subroutine dump3d( var3d, io, irec)
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
 real(sdsu_fps),dimension(ims:ime,jms:jme,kms:kme), intent(in)  :: var3d  ! assumed 3D array input 
                                                                           !  (specify starting index)
 integer,intent(in)    :: io   ! io index
 integer,intent(inout) :: irec ! binary record #
 integer :: i,j,k
 integer :: bnd(3)      ! upper memory bound for 3D array
 integer :: ks, ke      ! starting and ending k index

#if MPI < 2  

!    
! In this case var3d can be CRM or Satellite vertical levels 
!
    write(io,rec=irec) (((var3d(i,j,k),i=ids,ide),j=jds,jde),k=kds,kde)


#else


   !
   ! gather slave memory into master memory
   !
   if( .not. masterproc) then

   call mpi_sdsu_collect_tile( ke,    var3d ( ims:ime, jms:jme, kms:kme ),&
                               out_domain3d ( ids:ide, jds:jde, kds:kde )  )
   
   else
      out_domain3d ( ims:ime, jms:jme, kms:kme ) = var3d ( ims:ime, jms:jme, kms:kme )
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


end module module_radiation_driver

! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
! NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF NUWRF
