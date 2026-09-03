
! NTU multi-moment bulk microphysical scheme diagnosed version for WRFV3.8.1
! Released on May. 7, 2018 by Tzu-Chin Tsai at NTU (d97229002@ntu.edu.tw)
!==========================================================================

      MODULE module_mp_ntu
      USE    module_wrf_error
      USE    module_mp_radar
      USE    module_state_description, ONLY : P_QDCN,P_QTCN,P_QCCN,    &
                                              P_QRCN,P_QNIN

      IMPLICIT NONE
      PUBLIC  :: MP_NTU
      PRIVATE :: GAMMA,GAMLN,GAMIN,GAMMP,GSER,CFG,GUESS_RC,YEQU,DYEQU, &
                 PDF,DPDF,DLNX,POLYSVP

      INTEGER, PRIVATE, PARAMETER :: ID_NH42SO4 = 1,                   &! 1. (NH4)2SO4
                                     ID_DUST = 0,                      &! 2. DUST
                                     ID_IN = 2                          ! 3. IN
      INTEGER, PRIVATE, PARAMETER :: INSPEC = 2                         ! 1. SOOT; 2. SD; 3. AD; 4. BT
      INTEGER, PRIVATE, PARAMETER :: ICE_SHAPE = 1                      ! 0. SPHERICAL; 1. ASPHERICAL; PRISTINE ICE (HABIT)
      INTEGER, PRIVATE, PARAMETER :: AGG_SHAPE = 1                      ! 0. SPHERICAL; 1. CL94; ASPHERICAL AGGREGATES
      INTEGER, PRIVATE, PARAMETER :: ICE_RHOI = 1                       ! 0. FIXED; 1. PREDICT; 2. FIXED; PRISTINE ICE DENSITY
      INTEGER, PRIVATE, PARAMETER :: ICE_RHOS = 1                       ! 0. FIXED; 1. PREDICT; 2. MAPPING; AGGREGATES DENSITY
      INTEGER, PRIVATE, PARAMETER :: ICE_RHOG = 1                       ! 0. FIXED; 1. CL1993; 2. RH1985; PREDICT GRAUPEL DENSITY
      INTEGER, PRIVATE, PARAMETER :: ICE_VENT = 2                       ! 0. OFF; 1. HP76; 2. CL94; 3. JW99; VENTILATION EFFECT
      INTEGER, PRIVATE, PARAMETER :: HAIL_VENT = 1                      ! 0. DEFAULT; 1. CHENG ET AL. (2014); VENTILATION EFFECT
      INTEGER, PRIVATE, PARAMETER :: HWET_MODE = 1                      ! 0. OFF; 1. ON; HAILSTONE WET GROWTH
      INTEGER, PRIVATE, PARAMETER :: LIQ_VTC = 1                        ! 0. RAMS/FIXED; 1.,2. CL04; CLOUD DROPS/FLUX FALL SPEED
      INTEGER, PRIVATE, PARAMETER :: LIQ_VTR = 1                        ! 0. RAMS/FIXED; 1.,2. CL04; RAIN DROPS/FLUX FALL SPEED
      INTEGER, PRIVATE, PARAMETER :: ICE_VTI = 1                        ! 0. MORR/FIXED; 1. MH05; PRISTINE ICE FALL SPEED
      INTEGER, PRIVATE, PARAMETER :: ICE_VTS = 1                        ! 0. RAMS/FIXED; 1. MH05/EM16 SYNOPTIC; 2. SH09 CRYSTAL-FACE; AGGREGATES FALL SPEED
      INTEGER, PRIVATE, PARAMETER :: ICE_VTG = 1                        ! 0. RAMS/FIXED; 1. MH05; GRAUPEL FALL SPEED
      INTEGER, PRIVATE, PARAMETER :: ICE_VTH = 1                        ! 0. RAMS/FIXED; 1. MH05; HAIL FALL SPEED
      INTEGER, PRIVATE, PARAMETER :: AFAC_3M = 1                        ! 0. FIXED; 1. VARIABLE; AFAC PARAMETER
      INTEGER, PRIVATE, PARAMETER :: AFAR_3M = 1                        ! 0. FIXED; 1. VARIABLE; AFAR PARAMETER
      INTEGER, PRIVATE, PARAMETER :: AFAI_3M = 1                        ! 0. FIXED; 1.,2. VARIABLE AFAI PARAMETER
      INTEGER, PRIVATE, PARAMETER :: AFAS_3M = 1                        ! 0. FIXED; 1.,2. VARIABLE AFAS PARAMETER
      INTEGER, PRIVATE, PARAMETER :: AFAG_3M = 1                        ! 0. FIXED; 1.,2. VARIABLE AFAG PARAMETER
      INTEGER, PRIVATE, PARAMETER :: AFAH_3M = 1                        ! 0. FIXED; 1.,2. VARIABLE AFAH PARAMETER
      INTEGER, PRIVATE, PARAMETER :: NCCN = 3                           ! for aerosol modes
      INTEGER, PRIVATE, PARAMETER :: NAER = 2                           ! for aerosol categories of CCN and IN
      INTEGER, PRIVATE, PARAMETER :: MAER = 4, NAER1 = 4, NAER2 = 1
      INTEGER, PRIVATE, PARAMETER :: NAERT = NAER1+NAER2                ! 5
      INTEGER, PRIVATE, PARAMETER :: NTBXA = 25                         ! size of tables for the first guess of Rc
      INTEGER, PRIVATE, DIMENSION(NAER) :: IBAER
      INTEGER, PRIVATE, DIMENSION(NAER) :: NAERN(1:NAER)=(/NAER1,NAER2/)
      DOUBLE PRECISION, PRIVATE, SAVE, DIMENSION(NTBXA) :: TBLRC
      DOUBLE PRECISION, PRIVATE, SAVE, DIMENSION(NTBXA,NAER) :: TBLXA
!----------------- FOR MICROPHYSICAL PROCESSES -------------------------
      REAL, PRIVATE, PARAMETER :: DTMIN = 0.01                          ! [sec] mininum dt for cond./evap. calculation
      REAL, PRIVATE, PARAMETER :: DT20S = 20.                           ! [sec] dt for coalscence and accretion calculation
      REAL, PRIVATE, PARAMETER :: PI = 3.1415926535897932384626434
      REAL, PRIVATE, PARAMETER :: SQRTPI = 9.189385332046727417803297E-1
      REAL, PRIVATE, PARAMETER :: SQRT2 = 1.4142135623730950488016887
      REAL, PRIVATE, PARAMETER :: THRD = 1./3.,     C4PI3 = 4.*PI/3.
      REAL, PRIVATE, PARAMETER :: CP = 1.00546E3,   TK0C = 2.7315E2
      REAL, PRIVATE, PARAMETER :: R = 2.87058E2,    RV = 4.61495E2
      REAL, PRIVATE, PARAMETER :: CPI = 2.093E3,    CPW = 4.218E3
      REAL, PRIVATE, PARAMETER :: CMW = 1.8015E-2,  ISEPR = 0.6
      REAL, PRIVATE, PARAMETER :: RHOSU = 8.5E4/(2.8715E2*2.7315E2)
      REAL, PRIVATE, PARAMETER :: RHOW = 9.97E2,    RHOG1 = 4.E2
      REAL, PRIVATE, PARAMETER :: RHOI0 = 9.1E2,    iRHOI0 = 1./RHOI0
      REAL, PRIVATE, PARAMETER :: RHOS0 = RHOI0/(1.+2.*ISEPR+ISEPR**2.)
      REAL, PRIVATE, PARAMETER :: RHOG0 = (0.078+0.184*6.-0.015*36.)*1000
      REAL, PRIVATE, PARAMETER :: RHOH = 9.E2,      iRHOH = 1./RHOH
      REAL, PRIVATE, PARAMETER :: RHOI1 = 5.E2,     RHOS1 = 1.E2
      REAL, PRIVATE, PARAMETER :: RHOIMIN = 5.E1,   RHOIMAX = RHOI0
      REAL, PRIVATE, PARAMETER :: ROSA = 2.35738,   ROSB = -0.64016     ! Khain et al. (2004); Taka et al. (2012) SBM
!      REAL, PRIVATE, PARAMETER :: ROSA = 12.40442,  ROSB = -0.37404     ! Khain and Sednev (1995)
!      REAL, PRIVATE, PARAMETER :: ROSA = 0.946732,  ROSB = -0.6         ! Passarelli and Srivastava (1979)
      REAL, PRIVATE, PARAMETER :: C4PI3W = 4.*PI*RHOW/3.
      REAL, PRIVATE, PARAMETER :: iAMI0 = 6./(PI*RHOI0)
      REAL, PRIVATE, PARAMETER :: G = 9.806,        iAPW = 4./(PI*RHOW)
      REAL, PRIVATE, PARAMETER :: AAW = PI/4.,      V2M3 = 6./PI
      REAL, PRIVATE, PARAMETER :: BMW = 3.,         AMW = PI*RHOW/6.
      REAL, PRIVATE, PARAMETER :: BMI0 = 3.,        AMI0 = PI*RHOI0/6.
      REAL, PRIVATE, PARAMETER :: BMS = 3.,         AMS0 = PI*RHOS0/6.
      REAL, PRIVATE, PARAMETER :: BMG = 3.,         AMG0 = PI*RHOG0/6.
      REAL, PRIVATE, PARAMETER :: BMH = 3.,         AMH = PI*RHOH/6.
      REAL, PRIVATE, PARAMETER :: iAMW = 1./AMW,    iAMH = 1./AMH
      REAL, PRIVATE, PARAMETER :: AIMM = 6.6E-1,    BIMM = 1.E2
      REAL, PRIVATE, PARAMETER :: AVSG = 8.6E-1,    BVSG = 2.8E-1
      REAL, PRIVATE, PARAMETER :: AVRH = 7.8E-1,    BVRH = 3.08E-1
      REAL, PRIVATE, PARAMETER :: AVIS = 1.,        BVIS = 1.4E-1
      REAL, PRIVATE, PARAMETER :: VENC1 = 3.09E-2,  VENC2 = 1.447E-1
      REAL, PRIVATE, PARAMETER :: VENP1 = 1.05E-2,  VENP2 = 2.28E-2
      REAL, PRIVATE, PARAMETER :: VENH1 = 0.22385,  VENH2 = 0.00101
      REAL, PRIVATE, PARAMETER :: SCATF0 = 1.89E-1, ZXMIN1 = 1.E-28     ! 0.189 = 0.176/0.93
      REAL, PRIVATE, PARAMETER :: ZXMIN2 = 1.E-10,  ZXMIN3 = 1.E-15
      REAL, PRIVATE, PARAMETER :: DBZWL = 5.3E-2,   SCATI0 = 1.76E-1    ! S = 0.11; C = 0.053; X = 0.032m
      REAL, PRIVATE, PARAMETER :: NSMALL = 1.E-2,   NSMAL1 = 1.E2
      REAL, PRIVATE, PARAMETER :: RLIMIT = 1.E-32,  SLIMIT = 1.E-2
      REAL, PRIVATE, PARAMETER :: QSMAL1 = 1.E-9,   QSMALL = 1.E-14
      REAL, PRIVATE, PARAMETER :: QLIMIT = 1.E-6,   RSMALL = 1.E-20
      REAL, PRIVATE, PARAMETER :: ASMALL = 1.E-12,  ISMALL = 2.E-17
      REAL, PRIVATE, PARAMETER :: BOLTZ = 1.38E-23, EFMIN = 0.1
      REAL, PRIVATE, PARAMETER :: VTZ0 = 5.83,      VTC0 = 0.6          ! SURFACE ROUGHNESS PARAMETERS
      REAL, PRIVATE, PARAMETER :: VTC1 = 0.151931,  VTC2 = VTZ0**2./4.  ! VTC1=4./(VTZ0**2.*VTC0**5.E-1)
      REAL, PRIVATE, PARAMETER :: VTA0 = 1.7E-3,    VTB0 = 0.           ! VTB0 = 0.8
      REAL, PRIVATE, PARAMETER :: AVC0 = 3.E7,      BVC0 = 2.           ! Morrison scheme
      REAL, PRIVATE, PARAMETER :: AVR0 = 841.997,   BVR0 = 0.8          ! Morrison scheme
      REAL, PRIVATE, PARAMETER :: AVI0 = 700.,      BVI0 = 1.           ! Morrison scheme
      REAL, PRIVATE, PARAMETER :: AVS0 = 11.72,     BVS0 = 0.41         ! Locatelli and Hobbs (1974)
      REAL, PRIVATE, PARAMETER :: AVG0 = 19.3,      BVG0 = 0.37         ! Ferrier (1994)
      REAL, PRIVATE, PARAMETER :: AVH0 = 206.89,    BVH0 = 0.6384       ! Ferrier (1994)
!      REAL, PRIVATE, PARAMETER :: AVC0 = 3173.,     BVC0 = 2.           ! RAMS
!      REAL, PRIVATE, PARAMETER :: AVR0 = 149.,      BVR0 = 0.5          ! RAMS
!      REAL, PRIVATE, PARAMETER :: AVI0 = 71.34,     BVI0 = 0.6635       ! Ferrier (1994)
!      REAL, PRIVATE, PARAMETER :: AVS0 = 3.084,     BVS0 = 0.2          ! RAMS
!      REAL, PRIVATE, PARAMETER :: AVG0 = 93.3,      BVG0 = 0.5          ! RAMS
!      REAL, PRIVATE, PARAMETER :: AVH0 = 161.,      BVH0 = 0.5          ! RAMS
      REAL, PRIVATE, PARAMETER :: VTCMAX = 1.,      VTSMAX = 10.
      REAL, PRIVATE, PARAMETER :: VTIMAX = 10.,     VTRMAX = 15.
      REAL, PRIVATE, PARAMETER :: VTGMAX = 20.,     VTHMAX = 25.
      REAL, PRIVATE, PARAMETER :: AFAC0 = 0.,       AFAR0 = 0.
      REAL, PRIVATE, PARAMETER :: AFAI0 = 3.,       AFAS0 = 0.
      REAL, PRIVATE, PARAMETER :: AFAG0 = 0.,       AFAH0 = 0.
      REAL, PRIVATE, PARAMETER :: SASMAX = 0.95,    SASMIN = 0.05
      REAL, PRIVATE, PARAMETER :: KCIMIN = 0.556,   KCIMAX = 0.999      ! KC = 0.999, AFAMAX = 2996
      REAL, PRIVATE, PARAMETER :: KCSMIN = 0.223,   KCSMAX = 0.999      ! KC = 0.995, AFAMAX = 596
      REAL, PRIVATE, PARAMETER :: KCGMIN = 0.223,   KCGMAX = 0.999      ! KC = 0.99,  AFAMAX = 296
      REAL, PRIVATE, PARAMETER :: KCHMIN = 0.223,   KCHMAX = 0.999      ! KC = 0.95,  AFAMAX = 56
      REAL, PRIVATE, PARAMETER :: AFAMAX = 3.E4,    AFAMIN = 0.
      REAL, PRIVATE, PARAMETER :: AFU = 3.125E-1,   BFU = 1.0552E-3
      REAL, PRIVATE, PARAMETER :: CFU = -2.4023,    NR0 = 8.E6
      REAL, PRIVATE, PARAMETER :: NI0 = 1.E10,      NG0 = 4.E6
      REAL, PRIVATE, PARAMETER :: DCMIN = 1.E-7,    DCMAX = 1.E-4
      REAL, PRIVATE, PARAMETER :: DIMIN = 1.E-6,    DIMAX = 5.E-3
      REAL, PRIVATE, PARAMETER :: DRMIN = 3.E-5,    DRMAX = 6.E-3
      REAL, PRIVATE, PARAMETER :: DSMIN = 5.E-5,    DSMAX = 1.E-2
      REAL, PRIVATE, PARAMETER :: DGMIN = 5.E-5,    DGMAX = 1.E-2
      REAL, PRIVATE, PARAMETER :: DHMIN = 1.E-3,    DHMAX = 4.E-2
      REAL, PRIVATE, PARAMETER :: RCMIN = DCMIN/2., RCMAX = DCMAX/2.
      REAL, PRIVATE, PARAMETER :: RRMIN = DRMIN/2., RRMAX = DRMAX/2.
      REAL, PRIVATE, PARAMETER :: DI0 = 6.E-6,      DCR = 100.E-6
      REAL, PRIVATE, PARAMETER :: MI0 = AMI0*DI0**BMI0, NH0 = 4.E4
      REAL, PRIVATE, PARAMETER :: SIG1 = -6.2685,   SIG2 = -2.7312E-1
      REAL, PRIVATE, PARAMETER :: MNR1 = -2.2920,   SIG3 = 2.2606E-1
      REAL, PRIVATE, PARAMETER :: MNR3 = 3.4708E-1, MNR2 = -3.5158E-1
      REAL, DIMENSION(0:120) :: ITBLE
      DATA ITBLE /1.000000,0.979490,0.959401,0.939723,0.920450,        &
                  0.899498,0.879023,0.857038,0.833681,0.810961,        &
                  0.783430,0.755092,0.703072,0.537032,0.467735,        &
                  0.524807,0.630957,0.812831,1.096478,1.479108,        &
                  1.905461,2.089296,2.290868,2.398833,2.454709,        &
                  2.426610,2.371374,2.290868,2.137962,1.995262,        &
                  1.862087,1.737801,1.621810,1.513561,1.396368,        &
                  1.288250,1.188502,1.096478,1.000000,0.922571,        &
                  0.851138,0.785236,0.724436,0.668344,0.616595,        &
                  0.575440,0.537032,0.501187,0.467735,0.436516,        &
                  0.407380,0.380189,0.354813,0.331131,0.316228,        &
                  0.301995,0.291743,0.285102,0.281838,0.278612,        &
                  0.275423,0.278612,0.281838,0.285102,0.291743,        &
                  0.298538,0.309030,0.319890,0.331131,0.346737,        &
                  0.367282,0.393550,0.426580,0.457088,0.489779,        &
                  0.524807,0.562341,0.609537,0.660693,0.716143,        &
                  0.785236,0.860994,0.954993,1.047129,1.148154,        &
                  1.258925,1.380384,1.496236,1.603245,1.698244,        &
                  1.778279,1.840772,1.883649,1.905461,1.905461,        &
                  1.883649,1.862087,1.840772,1.798871,1.737801,        &
                  1.698244,1.640590,1.584893,1.548817,1.513561,        &
                  1.475707,1.452112,1.428894,1.412538,1.393157,        &
                  1.377209,1.361445,1.348963,1.336596,1.327394,        &
                  1.318257,1.309182,1.303167,1.294196,1.288250,1.279381/
      REAL, DIMENSION(0:167) :: IECC
      DATA IECC /0.00,0.03,0.29,0.25,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.10,0.42,0.50,0.47,0.21,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.15,0.49,0.56,0.55,0.46,0.10,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.20,0.52,0.61,0.62,0.59,0.49,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.25,0.62,0.72,0.75,0.74,    &
                 0.71,0.68,0.57,0.00,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.30,0.70,0.80,0.84,    &
                 0.85,0.85,0.84,0.83,0.81,0.77,0.69,0.10,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.35,0.75,0.85,    &
                 0.89,0.91,0.92,0.92,0.91,0.90,0.89,0.87,0.85,0.82,    &
                 0.71,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.40,0.80,    &
                 0.89,0.93,0.94,0.95,0.95,0.95,0.95,0.95,0.95,0.94,    &
                 0.93,0.92,0.91,0.88,0.80,0.15,0.00,0.00/
      REAL, DIMENSION(0:167) :: IEPC
      DATA IEPC /0.00,0.00,0.13,0.41,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.25,0.54,0.56,0.39,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.08,0.57,0.75,0.83,0.85,0.86,0.84,    &
                 0.78,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.10,0.65,0.80,0.87,0.90,0.91,    &
                 0.91,0.91,0.88,0.78,0.00,0.00,0.00,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.12,0.67,0.81,0.88,0.91,    &
                 0.92,0.93,0.93,0.92,0.91,0.89,0.81,0.00,0.00,0.00,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.15,0.68,0.82,0.89,    &
                 0.91,0.93,0.94,0.94,0.95,0.95,0.95,0.94,0.91,0.83,    &
                 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.18,0.69,0.82,    &
                 0.89,0.92,0.94,0.95,0.95,0.96,0.96,0.96,0.96,0.95,    &
                 0.94,0.91,0.80,0.00,0.00,0.00,0.00,0.00,0.22,0.70,    &
                 0.83,0.90,0.93,0.95,0.96,0.97,0.98,0.98,0.98,0.98,    &
                 0.97,0.96,0.94,0.92,0.84,0.00,0.00,0.00/
      REAL, DIMENSION(0:44) :: AMS1                                     ! SYNOPTIC CLOUDS
      DATA AMS1 /2.728E+01,1.258E+00,6.452E-02,6.259E-03,8.658E-04,    &
                 4.657E+01,1.193E+00,3.790E-02,2.688E-03,3.101E-04,    &
                 2.260E+01,1.000E+00,5.059E-02,4.889E-03,6.826E-04,    &
                 6.176E+00,6.002E-01,5.889E-02,8.892E-03,1.647E-03,    &
                 2.982E+00,1.795E-01,1.261E-02,1.647E-03,3.128E-04,    &
                 2.783E+00,2.179E-01,1.901E-02,2.825E-03,5.699E-04,    &
                 1.630E+00,3.267E-01,6.223E-02,1.493E-02,3.902E-03,    &
                 1.145E+01,3.346E-02,1.178E-04,7.701E-05,2.463E-05,    &
                 6.606E+01,1.799E-03,1.178E-04,7.701E-05,2.463E-05/
      REAL, DIMENSION(0:44) :: BMS1                                     ! SYNOPTIC CLOUDS
      DATA BMS1 /2.792,2.455,2.085,1.748,1.411,2.846,2.449,2.015,1.618,&
                 1.221,2.773,2.429,2.053,1.710,1.367,2.642,2.371,2.073,&
                 1.802,1.530,2.556,2.254,1.923,1.621,1.320,2.549,2.276,&
                 1.977,1.704,1.431,2.495,2.322,2.133,1.960,1.787,2.686,&
                 2.064,1.382,1.382,1.382,2.863,1.732,1.382,1.382,1.382/
      REAL, DIMENSION(0:44) :: AAS1                                     ! SYNOPTIC CLOUDS
      DATA AAS1 /1.782E+00,2.778E-01,4.794E-02,1.241E-02,4.069E-03,    &
                 2.571E+00,2.635E-01,3.180E-02,6.517E-03,1.879E-03,    &
                 1.910E+00,2.522E-01,3.759E-02,8.905E-03,2.787E-03,    &
                 6.238E-01,1.697E-01,4.686E-02,1.662E-02,6.761E-03,    &
                 5.604E-01,9.194E-02,1.726E-02,5.030E-03,1.927E-03,    &
                 2.802E-01,8.345E-02,2.602E-02,1.042E-02,4.854E-03,    &
                 1.269E-01,1.106E-01,8.453E-02,5.913E-02,3.763E-02,    &
                 3.569E-01,2.314E-02,1.664E-03,1.397E-03,8.521E-04,    &
                 6.339E-01,5.656E-03,1.664E-03,1.397E-03,8.521E-04/
      REAL, DIMENSION(0:44) :: BAS1                                     ! SYNOPTIC CLOUDS
      DATA BAS1 /2.133,1.938,1.725,1.531,1.337,2.170,1.932,1.671,1.432,&
                 1.194,2.140,1.927,1.693,1.480,1.267,2.027,1.882,1.722,&
                 1.576,1.431,2.011,1.821,1.612,1.422,1.232,1.941,1.810,&
                 1.666,1.534,1.403,1.861,1.842,1.821,1.801,1.782,1.960,&
                 1.669,1.350,1.350,1.350,2.018,1.509,1.350,1.350,1.350/
!      REAL, DIMENSION(0:44) :: AMS1                                     ! ANVIL CIRRUS CLOUDS
!      DATA AMS1 /4.980E+00,6.195E-01,7.945E-02,1.510E-02,3.493E-03,    &
!                 1.510E+01,8.733E-01,5.862E-02,7.289E-03,1.292E-03,    &
!                 7.328E+00,6.589E-01,6.422E-02,1.015E-02,2.124E-03,    &
!                 9.480E+00,6.589E-01,5.230E-02,7.241E-03,1.394E-03,    &
!                 7.690E+00,4.643E-01,3.552E-02,5.371E-03,1.259E-03,    &
!                 2.252E+00,2.551E-01,3.368E-02,7.386E-03,2.220E-03,    &
!                 8.449E-01,1.976E-01,4.887E-02,1.601E-02,6.165E-03,    &
!                 2.434E+00,4.203E-02,8.896E-04,7.212E-04,3.799E-04,    &
!                 1.028E+00,1.132E-02,1.697E-04,1.748E-04,1.170E-04/
!      REAL, DIMENSION(0:44) :: BMS1                                     ! ANVIL CIRRUS CLOUDS
!      DATA BMS1 /2.614,2.387,2.138,1.911,1.683,2.726,2.426,2.098,1.799,&
!                 1.499,2.653,2.394,2.110,1.850,1.591,2.679,2.394,2.083,&
!                 1.798,1.513,2.655,2.370,2.058,1.773,1.488,2.531,2.302,&
!                 2.051,1.822,1.593,2.432,2.273,2.100,1.941,1.782,2.533,&
!                 2.105,1.637,1.637,1.637,2.446,1.956,1.419,1.419,1.419/
!      REAL, DIMENSION(0:44) :: AAS1                                     ! ANVIL CIRRUS CLOUDS
!      DATA AAS1 /6.165E-01,1.983E-01,6.398E-02,2.527E-02,1.101E-02,    &
!                 1.431E+00,2.494E-01,4.757E-02,1.336E-02,4.687E-03,    &
!                 8.135E-01,2.001E-01,5.172E-02,1.767E-02,7.106E-03,    &
!                 9.162E-01,1.898E-01,4.212E-02,1.319E-02,5.029E-03,    &
!                 8.699E-01,1.638E-01,3.472E-02,1.099E-02,4.455E-03,    &
!                 2.869E-01,9.914E-02,3.552E-02,1.602E-02,8.159E-03,    &
!                 1.379E-01,8.313E-02,4.777E-02,2.933E-02,1.815E-02,    &
!                 4.723E-02,4.236E-02,3.724E-02,3.400E-02,3.019E-02,    &
!                 4.451E-02,7.812E-03,1.836E-03,2.590E-03,3.553E-03/
!      REAL, DIMENSION(0:44) :: BAS1                                     ! ANVIL CIRRUS CLOUDS
!      DATA BAS1 /2.023,1.899,1.763,1.639,1.515,2.108,1.925,1.724,1.541,&
!                 1.357,2.051,1.900,1.735,1.584,1.434,2.063,1.894,1.708,&
!                 1.539,1.370,2.055,1.885,1.698,1.528,1.358,1.943,1.828,&
!                 1.701,1.586,1.470,1.869,1.808,1.740,1.679,1.618,1.760,&
!                 1.753,1.746,1.746,1.746,1.754,1.561,1.350,1.350,1.350/
!----------------- FOR AEROSOL PARAMETER -------------------------------
      REAL, SAVE :: NC0
      REAL, SAVE, DIMENSION(NAER) :: SENS,ASH,CMAS,DNAS,AVAN,BETA1,RXMIN
      REAL, SAVE, DIMENSION(NCCN,NAER) :: ZCCN,CNMOD,CNSTD,WMAS,RFACT

      CONTAINS
!======================================================================
      SUBROUTINE NTU_INIT(PHB,PH,ALT,P,PB,QV,SCALAR,NSCAR,XLAND,BACKIN,&
                 CCNTY,RESTART,IDS,IDE,JDS,JDE,KDS,KDE,IMS,IME,JMS,JME,&
                 KMS,KME,ITS,ITE,JTS,JTE,KTS,KTE)
!======================================================================
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: RESTART
      INTEGER, INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE,IMS,IME,JMS,JME,  &
                             KMS,KME,ITS,ITE,JTS,JTE,KTS,KTE,CCNTY,NSCAR
      REAL, INTENT(IN) :: BACKIN
      REAL, INTENT(IN), DIMENSION(IMS:IME,JMS:JME) :: XLAND             ! 1:land ; 2:ocean
      REAL, INTENT(IN), DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: QV,PHB,  &
                        PH,ALT,P,PB
      REAL, INTENT(INOUT), DIMENSION(IMS:IME,KMS:KME,JMS:JME,NSCAR) :: &
                           SCALAR
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: DZ
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: RHO,DZ8W,P_PHY
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE,NAERT) :: QAERO
      INTEGER :: I,J,K,NK,ITF,JTF

      ITF = MIN(ITE,IDE-1)
      JTF = MIN(JTE,JDE-1)
      CALL AERO_CONST(CCNTY)
      DO J = JTS,JTE
         DO K = KTS,KTE
            DO I = ITS,ITE
               DZ(I,K,J) = (PHB(I,K,J)+PH(I,K,J))/G
            ENDDO
         ENDDO
      ENDDO

      DO J = JTS,JTF
         DO K = KTS,KTE-1
            DO I = ITS,ITF
               NK = KME-K
               DZ8W(I,K,J) = DZ(I,NK+1,J)-DZ(I,NK,J)
               RHO(I,K,J) = 1./ALT(I,NK,J)*(1.+QV(I,NK,J))
               P_PHY(I,K,J) = P(I,NK,J)+PB(I,NK,J)
               QAERO(I,K,J,1) = SCALAR(I,NK,J,P_QDCN)
               QAERO(I,K,J,2) = SCALAR(I,NK,J,P_QTCN)
               QAERO(I,K,J,3) = SCALAR(I,NK,J,P_QCCN)
               QAERO(I,K,J,4) = SCALAR(I,NK,J,P_QRCN)
               QAERO(I,K,J,5) = SCALAR(I,NK,J,P_QNIN)
            ENDDO
         ENDDO
      ENDDO
      IF (.NOT.RESTART) THEN
         CALL INIT_AEROSOL(P_PHY,RHO,DZ8W,XLAND,QAERO,BACKIN,IDS,IDE,  &
              JDS,JDE,KDS,KDE,IMS,IME,JMS,JME,KMS,KME,ITS,ITE,JTS,JTE, &
              KTS,KTE)
      ENDIF
      DO J = JTS,JTF
         DO K = KTS,KTE-1
            DO I = ITS,ITF
               NK = KME-K
               SCALAR(I,NK,J,P_QDCN) = QAERO(I,K,J,1)
               SCALAR(I,NK,J,P_QTCN) = QAERO(I,K,J,2)
               SCALAR(I,NK,J,P_QCCN) = QAERO(I,K,J,3)
               SCALAR(I,NK,J,P_QRCN) = QAERO(I,K,J,4)
               SCALAR(I,NK,J,P_QNIN) = QAERO(I,K,J,5)
            ENDDO
         ENDDO
      ENDDO
!+---+-----------------------------------------------------------------+
!..Set these variables needed for computing radar reflectivity.  These get used within radar_init to create other variables used in the radar module.
      xam_r = AMW
      xbm_r = BMW
      xmu_r = AFAR0
      xam_s = PI*RHOS1/6.
      xbm_s = BMS
      xmu_s = AFAS0
      xam_g = PI*RHOG1/6.
      xbm_g = BMG
      xmu_g = AFAG0

      CALL radar_init
!+---+-----------------------------------------------------------------+

      END SUBROUTINE NTU_INIT
!======================================================================
!
!======================================================================
      SUBROUTINE FIND_RC0(XAFRC,CMODE,CSTDV,WMAS,RC,TBLXA,TBLRC)    
!======================================================================
      IMPLICIT NONE
      INTEGER, PARAMETER :: ITERMAX = 50
      INTEGER :: IM,I
!----------------- AEROSOL TRI-MODAL SIZE DISTRIBUTION ----------------
      REAL :: RC                                                        ! cut-off radius [m]
      REAL, DIMENSION(NCCN) :: WMAS,                                   &! mass weighting for each mode
                               CMODE,                                  &! modal value (ln R) of each mode [m]
                               CSTDV                                    ! standard distribution (in logarithm) of the mode
      REAL, PARAMETER :: RC_MAX = 99., RC_MIN = 1.E-9, XTOR = 1.E-5
      REAL*8, DIMENSION(NCCN) :: DWMAS,DMODE,DSTDV
!----------------- GET THE FIRST GUESS OF RC ---------------------------
      REAL*8 :: DRC,UL,X1MAFRC,Y0,Y,DY,DX,DERF,DLOG,DEXP,DSQRT,XAFRC    ! ratio of dry/total aerosol mass
      REAL*8, DIMENSION(NTBXA) :: TBLXA,                               &! table of dry aerosol mass fraction
                                  TBLRC                                 ! table of curresponding R-cut-off

      DO IM = 1,NCCN
         DWMAS(IM) = DBLE(WMAS(IM))
         DSTDV(IM) = DBLE(CSTDV(IM))
         DMODE(IM) = DLOG(DBLE(CMODE(IM)))+DBLE(3.*CSTDV(IM)*CSTDV(IM)) ! for mass, the 3rd momentum DBLE(CMODE (IM))+
      ENDDO
      UL = (1.-DERF(5.D+0/DSQRT(2.D+0)))*0.5                            ! after 5 sigma, x1mfrc will < 2.86E-7
      X1MAFRC = 1.D+0-XAFRC
      IF (XAFRC.LT.1.E-10) THEN
         RC = RC_MIN
         RETURN
      ENDIF
      IF (X1MAFRC.LT.UL) THEN
         RC = RC_MAX
         RETURN
      ENDIF
      DRC = GUESS_RC(XAFRC,TBLXA,TBLRC)                                 ! cut-off radius [m]
      Y0  = YEQU(DRC,XAFRC,DWMAS,DMODE,DSTDV)
      DO I = 1,ITERMAX
         Y  = Y0
         DY = DYEQU(DRC,DWMAS,DMODE,DSTDV)
         IF (DY.LE.1.E-50) THEN
            PRINT *,'IN FIND_RC0.F DY IS',DY
            PRINT *,I,DX,Y,DY,DEXP(DRC),DRC,XAFRC,X1MAFRC,DMODE,DSTDV, &
                    CMODE,CSTDV
            STOP
         ENDIF
         DX = -Y/DY
         IF (DABS(DX).GT.3.D+0) DX = DX/(DABS(DX/0.4))
         IF (DABS(DX).LT.XTOR) THEN
            RC = REAL(DEXP(DRC+DX))
            RC = MIN(RC_MAX,MAX(RC,RC_MIN))
            RETURN
         ENDIF
         Y0 = YEQU(DRC+DX,XAFRC,DWMAS,DMODE,DSTDV)
         IF (Y*Y0.LT.0.) THEN
            DX = DX*5.D-1
            Y0 = YEQU(DRC+DX,XAFRC,DWMAS,DMODE,DSTDV)
         ENDIF
         DRC = DRC+DX
      ENDDO
      RC = REAL(DEXP(DRC))

      END SUBROUTINE FIND_RC0
!======================================================================
!
!======================================================================
      FUNCTION GUESS_RC(XAFRC,TBLXA,TBLRC)
!======================================================================
      REAL*8 :: GUESS_RC,XAFRC,TBLXA(NTBXA),TBLRC(NTBXA)
      INTEGER :: IBNG,IEND,I

      IF (XAFRC.GE.TBLXA(NTBXA)) THEN
         GUESS_RC = TBLRC(NTBXA)
         RETURN
      ELSEIF (XAFRC.LT.TBLXA(1)) THEN
         GUESS_RC = TBLRC(1)
         RETURN
      ENDIF
      IF (XAFRC.LT.TBLXA(NTBXA/4)) THEN
         IBNG = 1
         IEND = NTBXA/4-1
      ELSEIF (XAFRC.LT.TBLXA(NTBXA*2/4)) THEN
         IBNG = NTBXA/4
         IEND = NTBXA*2/4-1
      ELSEIF (XAFRC.LT.TBLXA(NTBXA*3/4)) THEN
         IBNG = NTBXA*2/4
         IEND = NTBXA*3/4-1
      ELSE
         IBNG = NTBXA*3/4
         IEND = NTBXA-1
      ENDIF
      DO I = IEND,IBNG,-1
         IF (XAFRC.GE.TBLXA(I)) GOTO 111
      ENDDO
111   GUESS_RC = TBLRC(I)

      END FUNCTION GUESS_RC
!======================================================================
!
!======================================================================
      FUNCTION  YEQU(X,XAFRC,WEGHT,XLBAR,SIGMA)
!======================================================================
      INTEGER :: I
      REAL*8 :: YEQU,DX,XAFRC,X,DSQRT2,DSQRT
      REAL*8, DIMENSION(NCCN) :: WEGHT,XLBAR,SIGMA

      DSQRT2 = DSQRT(2.D+0)
      YEQU   = 0.D+0
      DO I = 1,NCCN
         DX = (X-XLBAR(I))/SIGMA(I)/DSQRT2
         YEQU = YEQU+WEGHT(I)*PDF(DX)
      ENDDO
      YEQU = YEQU-XAFRC

      END FUNCTION YEQU
!----------------------------------------------------------------------
      FUNCTION DYEQU(X,WEGHT,XLBAR,SIGMA)
!======================================================================
      INTEGER :: I
      REAL*8 :: DYEQU,X,DX,DSQRT2,DSQRT
      REAL*8, DIMENSION(NCCN) :: WEGHT,XLBAR,SIGMA 

      DSQRT2 = DSQRT(2.D+0)
      DYEQU  = 0.D+0
      DO I = 1,NCCN
         DX = (X-XLBAR(I))/SIGMA(I)/DSQRT2
         DYEQU = DYEQU+WEGHT(I)*DPDF(DX)/SIGMA(I)/DSQRT2
      ENDDO

      END FUNCTION DYEQU
!======================================================================
      FUNCTION  PDF(X)
      REAL*8 :: PDF,X,DERF

      PDF = (1.D+0+DERF(X))*5.D-1

      END FUNCTION PDF
!----------------------------------------------------------------------
      FUNCTION DPDF(X)
      REAL*8 :: DPDF,X,DPI,DACOS,DEXP,DSQRT

      DPI  = DACOS(-1.D+0)
      DPDF = DEXP(-X*X)/DSQRT(DPI)

      END FUNCTION DPDF
!======================================================================
!
!======================================================================
      REAL FUNCTION GAMMA(X)
!======================================================================
      IMPLICIT NONE
      INTEGER :: I,N
      LOGICAL :: PARITY
      REAL  :: CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,TWO,X,XBIG,       &
               XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      REAL, DIMENSION(7) :: C
      REAL, DIMENSION(8) :: P
      REAL, DIMENSION(8) :: Q
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P /-1.71618513886549492533811E0,2.47656508055759199108314E1,&
              -3.79804256470945635097577E2,6.29331155312818442661052E2,&
              8.66966202790413211295064E2,-3.14512729688483675254357E4,&
              -3.61444134186911729807069E4,6.64561438202405440627855E4/
      DATA Q /-3.08402300119738975254353E1,3.15350626979604161529144E2,&
             -1.01515636749021914166146E3,-3.10777167157231109440444E3,&
              2.25381184209801510330112E4,4.75584627752788110767815E3, &
             -1.34659959864969306392456E5,-1.15132259675553483497211E5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C /-1.910444077728E-3,8.4171387781295E-4,                   &
             -5.952379913043012E-4,7.93650793500350248E-4,             &
             -2.777777777777681622553E-3,8.333333333333333331554247E-2,&
             5.7083835261E-3/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y.LE.ZERO) THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
         Y   = -X
         Y1  = AINT(Y)
         RES = Y-Y1
         IF (RES.NE.ZERO) THEN
            IF (Y1.NE.AINT(Y1*HALF)*TWO) PARITY = .TRUE.
            FACT = -PI/SIN(PI*RES)
            Y = Y+ONE
         ELSE
            RES = XINF
            GOTO 900
         ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF (Y.LT.EPS) THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
         IF (Y.GE.XMININ) THEN
            RES = ONE/Y
         ELSE
            RES = XINF
            GOTO 900
         ENDIF
      ELSEIF (Y.LT.TWELVE) THEN
         Y1 = Y
         IF (Y.LT.ONE) THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
            Z = Y
            Y = Y+ONE
         ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
            N = INT(Y)-1
            Y = Y-CONV(N)
            Z = Y-ONE
         ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
         XNUM = ZERO
         XDEN = ONE
         DO I = 1,8
            XNUM = (XNUM+P(I))*Z
            XDEN = XDEN*Z+Q(I)
         END DO
         RES = XNUM/XDEN+ONE
         IF (Y1.LT.Y) THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
            RES = RES/Y1
         ELSEIF (Y1.GT.Y) THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
            DO I = 1,N
               RES = RES*Y
               Y = Y+ONE
            END DO
         ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
         IF (Y.LE.XBIG) THEN
            YSQ = Y*Y
            SUM = C(7)
            DO I = 1,6
               SUM = SUM/YSQ+C(I)
            END DO
            SUM = SUM/Y-Y+SQRTPI
            SUM = SUM+(Y-HALF)*LOG(Y)
            RES = EXP(SUM)
         ELSE
            RES = XINF
            GOTO 900
         ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT.NE.ONE) RES = FACT/RES
  900 GAMMA = RES
      RETURN

      END FUNCTION GAMMA
!======================================================================

!======================================================================
      REAL FUNCTION POLYSVP(T,TYPE)
!======================================================================
!  COMPUTE SATURATION VAPOR PRESSURE POLYSVP RETURNED IN UNITS OF PA. T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
!  REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)
      IMPLICIT NONE
      REAL :: DUM,T,DT
      INTEGER :: TYPE
! ice
      REAL :: a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i 
      DATA a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /6.11147274,0.503160820,&
           0.188439774E-1,0.420895665E-3,0.615021634E-5,0.602588177E-7,&
           0.385852041E-9,0.146898966E-11,0.252751365E-14/
! liquid
      REAL :: a0,a1,a2,a3,a4,a5,a6,a7,a8 
      DATA a0,a1,a2,a3,a4,a5,a6,a7,a8 /6.11239921,0.443987641,         &
           0.142986287E-1,0.264847430E-3,0.302950461E-5,0.206739458E-7,&
           0.640689451E-10,-0.952447341E-13,-0.976195544E-15/
      IF (TYPE.EQ.1) THEN                                               ! ICE
         DT = MAX(-80.,T-273.16)
         POLYSVP = a0i+DT*(a1i+DT*(a2i+DT*(a3i+DT*(a4i+DT*             &
                          (a5i+DT*(a6i+DT*(a7i+a8i*DT))))))) 
         POLYSVP = POLYSVP*100.
      ENDIF
      IF (TYPE.EQ.0) THEN                                               ! LIQUID
         DT = MAX(-80.,T-273.16)
         POLYSVP = a0+DT*(a1+DT*(a2+DT*(a3+DT*(a4+DT*(a5+DT*           &
                         (a6+DT*(a7+a8*DT)))))))
         POLYSVP = POLYSVP*100.
      ENDIF

      END FUNCTION POLYSVP
!======================================================================
!
!======================================================================
      REAL FUNCTION GAMLN(XX)
!======================================================================
!  Returns value of ln(GAMMA(XX)) for XX>0 (modified from "Numerical Recipes")
      IMPLICIT NONE
      REAL, INTENT(IN) :: XX
      INTEGER :: J
      REAL*8 :: ser,stp,TMP,X,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,&
                  -1.231739572450155d0,.1208650973866179d-2,-.5395239384953d-5,2.5066282746310005d0/

      X = DBLE(XX)
      y = X
      TMP = X+5.5D0
      TMP = (X+0.5D0)*LOG(TMP)-TMP
      ser = 1.000000000190015d0
      DO J = 1,6   !original
         y = y+1.D0
         ser = ser+cof(J)/y
      ENDDO
#if (DWORDSIZE == 8 && RWORDSIZE == 8)
      GAMLN = TMP+LOG(stp*ser/X)
#elif (DWORDSIZE == 8 && RWORDSIZE == 4)
      GAMLN = SNGL(TMP+LOG(stp*ser/X))
#else
!      This is a temporary hack assuming double precision is 8 bytes.
#endif

      END FUNCTION GAMLN
!======================================================================
!
!======================================================================
      REAL FUNCTION GAMIN(P,XMAX)
!======================================================================
! The incomplete gamma function below the limit of P = AFA+1 and XMAX = DI0*LAMDA
      REAL :: P,XMAX

      GAMIN = GAMMP(P,XMAX)*EXP(GAMLN(P))

      END FUNCTION GAMIN
!=======================================================================
      REAL FUNCTION GAMMP(A,X)
!=======================================================================
! The fraction of distribution below the limit of A = AFA+1 and X = Db*LAMDA
      IMPLICIT NONE
      REAL :: A,X,GAMMCF,GAMSER,GLN

      IF (X.LT.0..OR.A.LE.0.) &
         CALL wrf_error_fatal('warning : bad arguments in gammq')
      IF (X.LT.A+1.) THEN
         CALL GSER(GAMSER,A,X,GLN)
         GAMMP = GAMSER
      ELSE
         CALL CFG(GAMMCF,A,X,GLN)
         GAMMP = 1.-GAMMCF
      ENDIF
      RETURN

      END FUNCTION GAMMP
!======================================================================
      SUBROUTINE GSER(GAMSER,A,X,GLN)
!======================================================================
! USES GAMLN, Returns the incomplete gamma function P(A,X) evaluated by its series representation as GAMSER. Also returns GAMMA(A) as GLN.
      IMPLICIT NONE
      INTEGER :: N
      REAL :: A,GAMSER,GLN,X,AP,de1,summ
      INTEGER, PARAMETER :: ITMAX = 500
      REAL, PARAMETER :: EPS = 3.E-7

      GLN = GAMLN(A)
      IF (X.LE.0.) THEN
         IF (X.LT.0.) CALL wrf_error_fatal ( 'WARNING: X <0 in GSER' )
         GAMSER = 0.
         RETURN
      ENDIF
      AP = A
      summ = 1./A
      de1 = summ
      DO N = 1,ITMAX
         AP = AP+1.
         de1 = de1*X/AP
         summ = summ+de1
         IF (ABS(de1).LT.ABS(summ)*EPS) GOTO 777
      ENDDO
      CALL wrf_error_fatal ('Warning : ITMAX too small in GSER')
 777  GAMSER = summ*EXP(-X+A*LOG(X)-GLN)

      RETURN

      END SUBROUTINE GSER
!=======================================================================
      SUBROUTINE CFG(GAMMCF,A,X,GLN)
!=======================================================================
! USES GAMLN, Returns the incomplete gamma function (Q(A,X) evaluated by tis continued fraction
! representation as GAMMCF.  Also returns ln(GAMMA(A)) as GLN.  ITMAX is the maximum allowed number of iterations;
! EPS is the relative accuracy; FPMIN is a number near the smallest representable floating-point number.
      IMPLICIT NONE
      INTEGER :: I
      REAL :: A,GAMMCF,GLN,X,AN,b,c,d,de1,h
      INTEGER, PARAMETER :: ITMAX = 500
      REAL, PARAMETER :: EPS = 3.E-7
      REAL, PARAMETER :: fpmin = 1.E-30

      GLN = GAMLN(A)
      b = X+1.-A
      c = 1./fpmin
      d = 1./b
      h = d
      DO I = 1,ITMAX
         AN = -I*(I-A)
         b = b+2.
         d = AN*d+b
         IF (ABS(d).LT.fpmin) d = fpmin
         c = b+AN/c
         IF (ABS(c).LT.fpmin) c = fpmin
         d = 1./d
         de1 = d*c
         h = h*de1
         IF (ABS(de1-1.).LT.EPS) GOTO 888
      ENDDO
      CALL wrf_error_fatal ('Warning : ITMAX too small in gcf')
 888  GAMMCF = EXP(-X+A*LOG(X)-GLN)*h
      RETURN

      END SUBROUTINE CFG
!=======================================================================
!
!======================================================================
      SUBROUTINE SOLVE_AFAC(QC1D,NC1D,LAMC,MVDC,AFAC)
!======================================================================
      IMPLICIT NONE
      REAL :: QC1D,NC1D,LAMC,GC1,MVDC,AFAC,LAMCMIN,LAMCMAX,C3M1D,SIGC, &
              MNRC

      IF (NC1D.LT.NSMAL1) THEN
         NC1D = NC0
      ENDIF
      IF (AFAC_3M.EQ.0) THEN
         AFAC = AFAC0
      ELSEIF (AFAC_3M.EQ.1) THEN
         SIGC = EXP(SIG1+SIG2*LOG(NC1D)+SIG3*LOG(QC1D))
         MNRC = EXP(MNR1+MNR2*LOG(NC1D)+MNR3*LOG(QC1D))
         AFAC = MIN(MAX(SIGC/MNRC,AFAMIN),AFAMAX)
      ENDIF
      C3M1D   = QC1D*V2M3/RHOW
      GC1     = GAMLN(AFAC+1.)
      LAMC    = (EXP(GAMLN(AFAC+4.)-GC1)*NC1D/C3M1D)**THRD
      LAMCMIN = (EXP(GAMLN(AFAC+4.)-GAMLN(AFAC+1.)))**THRD/DCMAX
      LAMCMAX = (EXP(GAMLN(AFAC+4.)-GAMLN(AFAC+1.)))**THRD/DCMIN
      IF (LAMC.LT.LAMCMIN) THEN
         LAMC = LAMCMIN
         NC1D = C3M1D*EXP(GAMLN(AFAC+1.)-GAMLN(AFAC+4.)+3.*LOG(LAMC))
      ELSEIF (LAMC.GT.LAMCMAX) THEN
         LAMC = LAMCMAX
         NC1D = C3M1D*EXP(GAMLN(AFAC+1.)-GAMLN(AFAC+4.)+3.*LOG(LAMC))
      ENDIF
      MVDC = (EXP(GAMLN(AFAC+4.)-GAMLN(AFAC+1.)))**THRD/LAMC

      END SUBROUTINE SOLVE_AFAC
!======================================================================
!
!======================================================================
      SUBROUTINE SOLVE_AFAR(QR1D,NR1D,LAMR,MVDR,AFAR)
!======================================================================
      IMPLICIT NONE
      REAL :: QR1D,NR1D,LAMR,GR1,MVDR,LAMRMIN,LAMRMAX,R3M1D,AFAR,BDR

      IF (NR1D.LT.NSMALL) THEN
         NR1D = MAX(1.,(QR1D/PI/RHOW)**0.25*NR0**0.75)
      ENDIF
      R3M1D = QR1D*V2M3/RHOW
      IF (AFAR_3M.EQ.0) THEN
         AFAR = AFAR0
      ELSEIF (AFAR_3M.EQ.1) THEN
         BDR  = MIN(MAX((R3M1D/NR1D)**THRD,DRMIN),DRMAX)
         AFAR = MAX(AFAMIN,19.*TANH(0.6*(1.E3*BDR-1.8))+17.)
      ENDIF
      GR1     = GAMLN(AFAR+1.)
      LAMR    = (EXP(GAMLN(AFAR+4.)-GR1)*NR1D/R3M1D)**THRD
      LAMRMIN = (EXP(GAMLN(AFAR+4.)-GAMLN(AFAR+1.)))**THRD/DRMAX
      LAMRMAX = (EXP(GAMLN(AFAR+4.)-GAMLN(AFAR+1.)))**THRD/DRMIN
      IF (LAMR.LT.LAMRMIN) THEN
         LAMR = LAMRMIN
         NR1D = R3M1D*EXP(GAMLN(AFAR+1.)-GAMLN(AFAR+4.)+3.*LOG(LAMR))
      ELSEIF (LAMR.GT.LAMRMAX) THEN
         LAMR = LAMRMAX
         NR1D = R3M1D*EXP(GAMLN(AFAR+1.)-GAMLN(AFAR+4.)+3.*LOG(LAMR))
      ENDIF
      MVDR = (EXP(GAMLN(AFAR+4.)-GAMLN(AFAR+1.)))**THRD/LAMR

      END SUBROUTINE SOLVE_AFAR
!======================================================================

!======================================================================
      SUBROUTINE SOLVE_AFAI(TK1D,RHO,QI1D,NI1D,VI1D,FI1D,I2M1D,I3M1D,  &
                 ADAGR,ZETA,LAMI,AFAI,MVDI,RHOI,AMI,BMI,AVI,BVI,BEST)
!======================================================================
      IMPLICIT NONE
      REAL :: TK1D,RHO,QI1D,NI1D,VI1D,FI1D,I2M1D,I3M1D,ADAGR,ZETA,RHOI,&
              LAMI,AFAI,AMI,BMI,AVI,BVI,MVDI,GI1
      REAL :: I3M0,BEST0,LAMIMIN,LAMIMAX,IPF,IPG,KDX,ZETA2,ZETA4,KINV, &
              AAI,BAI,IBA1,C1X2,VTA1,VTB1,IPH2,IPG2,AVIA0,BEST,BDI

      IF (ICE_RHOI.EQ.0) THEN
         RHOI = RHOI0
         VI1D = 0.
      ELSEIF (ICE_RHOI.EQ.1) THEN
         RHOI = QI1D/(VI1D+ISMALL)
         IF (RHOI.LT.RHOIMIN) THEN
            RHOI = RHOIMIN
         ELSEIF (RHOI.GT.RHOIMAX) THEN
            RHOI = RHOIMAX
         ENDIF
         VI1D = QI1D/RHOI
      ELSEIF (ICE_RHOI.EQ.2) THEN
         RHOI = RHOI1
         VI1D = 0.
      ENDIF
      IF (NI1D.LT.NSMAL1) THEN
         NI1D = MAX(1.,(QI1D/PI/RHOI)**0.25*NI0**0.75)
      ENDIF
      IF (I2M1D.GE.ASMALL.AND.I3M1D.GE.ISMALL) THEN
         KDX = MAX(0.,I2M1D**3./(I3M1D**2.*NI1D))
         IF (KDX.GT.KCIMAX) THEN
            KDX   = KCIMIN
            I2M1D = (KDX*NI1D*I3M1D**2.)**THRD
         ELSEIF (KDX.LT.KCIMIN) THEN
            KDX   = KCIMIN
            I2M1D = (KDX*NI1D*I3M1D**2.)**THRD
         ENDIF
         AFAI = (6.*KDX-3.+SQRT(8.*KDX+1.))/(2.-2.*KDX)
         AFAI = MIN(MAX(AFAI,AFAMIN),AFAMAX)
!         LAMI = SQRT(NI1D*(AFAI+2.)*(AFAI+1.)/I2M1D)
!         LAMI = (NI1D*(AFAI+1.)*(AFAI+2.)*(AFAI+3.)/I3M1D)**THRD
         LAMI = (AFAI+3.)*I2M1D/I3M1D
      ELSEIF (I2M1D.LT.ASMALL.AND.I3M1D.GE.ISMALL) THEN
         IF (AFAI_3M.EQ.0) THEN
            AFAI  = AFAI0
            I2M1D = 0.
         ELSEIF (AFAI_3M.EQ.1) THEN
            AFAI  = AFAI0
            I2M1D = (KCIMIN*NI1D*I3M1D**2.)**THRD
         ELSEIF (AFAI_3M.EQ.2) THEN
            BDI   = MIN(MAX((I3M1D/NI1D)**THRD,DIMIN),DIMAX)
            AFAI  = MAX(AFAMIN,12.*TANH(0.7*(1.E3*BDI-1.7))+11.)
            I2M1D = 0.
         ENDIF
         LAMI = (EXP(GAMLN(AFAI+4.)-GAMLN(AFAI+1.))*NI1D/I3M1D)**THRD
      ELSE
         IF (AFAI_3M.EQ.0) THEN
            AFAI  = AFAI0
            I2M1D = 0.
         ELSEIF (AFAI_3M.EQ.1) THEN
            AFAI  = AFAI0
            I2M1D = (KCIMIN*NI1D*(QI1D*V2M3/RHOI)**2.)**THRD
         ELSEIF (AFAI_3M.EQ.2) THEN
            BDI   = MIN(MAX((QI1D*V2M3/NI1D/RHOI)**THRD,DIMIN),DIMAX)
            AFAI  = MAX(AFAMIN,12.*TANH(0.7*(1.E3*BDI-1.7))+11.)
            I2M1D = 0.
         ENDIF
         GI1  = GAMLN(AFAI+1.)
         LAMI = (EXP(GAMLN(AFAI+4.)-GI1)*PI*RHOI*NI1D/QI1D/6.)**THRD
      ENDIF
      LAMIMIN = (EXP(GAMLN(AFAI+4.)-GAMLN(AFAI+1.)))**THRD/DIMAX
      LAMIMAX = (EXP(GAMLN(AFAI+4.)-GAMLN(AFAI+1.)))**THRD/DIMIN
      IF (LAMI.LT.LAMIMIN) THEN
         LAMI = LAMIMIN
         NI1D = I3M1D*EXP(GAMLN(AFAI+1.)-GAMLN(AFAI+4.)+3.*LOG(LAMI))
      ELSEIF (LAMI.GT.LAMIMAX) THEN
         LAMI = LAMIMAX
         NI1D = I3M1D*EXP(GAMLN(AFAI+1.)-GAMLN(AFAI+4.)+3.*LOG(LAMI))
      ENDIF
      MVDI = (EXP(GAMLN(AFAI+4.)-GAMLN(AFAI+1.)))**THRD/LAMI
      IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
         I3M0 = NI1D*DI0**3.
         IF (MVDI.GT.(DI0+1.E-7)) THEN
            IF (ICE_SHAPE.EQ.0) THEN
               ZETA = 0.; ADAGR = 1.
               FI1D = I3M1D
            ELSEIF (ICE_SHAPE.EQ.1) THEN
               ZETA = LOG(FI1D/I3M1D)/LOG(I3M1D/I3M0)
               IF (ZETA.GT.0.4) THEN
                  ZETA = 0.4
                  FI1D = (I3M1D/I3M0)**ZETA*I3M1D
               ELSEIF (ZETA.LT.(-0.4)) THEN
                  ZETA = -0.4
                  FI1D = (I3M1D/I3M0)**ZETA*I3M1D
               ENDIF
               ADAGR = (1.+2.*ZETA)/(1.-ZETA)
            ENDIF
         ELSE
            FI1D = I3M1D
            ZETA = 0.; ADAGR = 1.
         ENDIF
      ELSE
         FI1D = I3M1D 
         ZETA = 0.; ADAGR = 1.
      ENDIF
      IF ((ADAGR-1.).GE.SLIMIT) THEN
         AMI = PI*RHOI*DI0**(2.-2./ADAGR)/6.
         BMI = 2./ADAGR+1.
      ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
         AMI = PI*RHOI*DI0**(1.-ADAGR)/6.
         BMI = ADAGR+2.
      ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
         AMI = PI*RHOI/6.
         BMI = BMI0
      ENDIF
      IF (ICE_VTI.EQ.0) THEN
         AVI = AVI0
         BVI = BVI0
      ELSEIF (ICE_VTI.EQ.1) THEN
         IF (I3M1D.GE.ISMALL) THEN
            KINV  = (1.72E-5*(393./(TK1D+120.))*(TK1D/TK0C)**1.5)/RHO
            BEST0 = 2.*G*NI1D/(KINV**2.)
            IF ((ADAGR-1.).GE.SLIMIT) THEN
               AAI   = PI/4./DI0**ZETA
               BAI   = 3.*(ADAGR+1.)/(ADAGR+2.)
               IPH2  = 2.*3.*ADAGR/(ADAGR+2.)
               ZETA4 = 4.*(ADAGR-1.)/(ADAGR+2.)
               IBA1  = BMI+IPH2-BAI
               BEST  = BEST0*AMI*EXP(GAMLN(IBA1+AFAI+1.)-GAMLN(AFAI+   &
                       1.)-IBA1*LOG(LAMI))/(AAI*DI0**ZETA4)
               C1X2  = VTC1*BEST**5.E-1
               VTB1  = C1X2/(1.+C1X2)**5.E-1/((1.+C1X2)**5.E-1-1.)/2.- &
                       VTA0*VTB0*BEST**VTB0/VTC2/(SQRT(1.+C1X2)-1.)**2.
               VTA1  = MAX((VTC2*((1.+C1X2)**5.E-1-1.)**2.-VTA0*BEST** &
                       VTB0)/BEST**VTB1,0.)
               AVIA0 = VTA1*KINV**(1.-2.*VTB1)
               AVI   = AVIA0*(2.*AMI*G/(AAI*DI0**ZETA4))**VTB1
               BVI   = VTB1*(BMI+IPH2-BAI)-1.
            ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
               ZETA2 = 2.*(ADAGR-1.)/(ADAGR+2.)
               AAI   = PI*DI0**ZETA2/4.
               BAI   = 2.*3./(ADAGR+2.)
               IPG2  = 2.*3./(ADAGR+2.)
               IBA1  = BMI+IPG2-BAI
               BEST  = BEST0*AMI*DI0**ZETA2*EXP(GAMLN(IBA1+AFAI+1.)-   &
                       GAMLN(AFAI+1.)-IBA1*LOG(LAMI))/AAI
               C1X2  = VTC1*BEST**5.E-1
               VTB1  = C1X2/(1.+C1X2)**5.E-1/((1.+C1X2)**5.E-1-1.)/2.- &
                       VTA0*VTB0*BEST**VTB0/VTC2/(SQRT(1.+C1X2)-1.)**2.
               VTA1  = MAX((VTC2*((1.+C1X2)**5.E-1-1.)**2.-VTA0*BEST** &
                       VTB0)/BEST**VTB1,0.)
               AVIA0 = VTA1*KINV**(1.-2.*VTB1)
               AVI   = AVIA0*(2.*AMI*G*DI0**ZETA2/AAI)**VTB1
               BVI   = VTB1*(BMI+IPG2-BAI)-1.
            ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
               AAI   = PI/4.
               BAI   = 2.
               IBA1  = BMI
               BEST  = BEST0*AMI*EXP(GAMLN(IBA1+AFAI+1.)-GAMLN(AFAI+   &
                       1.)-IBA1*LOG(LAMI))/AAI
               C1X2  = VTC1*BEST**5.E-1
               VTB1  = C1X2/(1.+C1X2)**5.E-1/((1.+C1X2)**5.E-1-1.)/2.- &
                       VTA0*VTB0*BEST**VTB0/VTC2/(SQRT(1.+C1X2)-1.)**2.
               VTA1  = MAX((VTC2*((1.+C1X2)**5.E-1-1.)**2.-VTA0*BEST** &
                       VTB0)/BEST**VTB1,0.)
               AVIA0 = VTA1*KINV**(1.-2.*VTB1)
               AVI   = AVIA0*(2.*AMI*G/AAI)**VTB1
               BVI   = VTB1*BMI-1.
            ENDIF
         ELSE
            AVI = AVI0
            BVI = BVI0
         ENDIF
      ENDIF

      END SUBROUTINE SOLVE_AFAI
!======================================================================
!
!======================================================================
      SUBROUTINE SOLVE_AFAS(TK1D,RHO,QS1D,NS1D,VS1D,FS1D,S2M1D,AFAS,   &
                 LAMS,MVDS,RHOS,SASPR,AMS,AVS,BVS,SMLF)
!======================================================================
      IMPLICIT NONE
      INTEGER :: TBIN,DBIN
      REAL :: TK1D,RHO,QS1D,NS1D,VS1D,FS1D,S2M1D,AFAS,RHOS,SASPR,LAMS, &
              AMS,AVS,BVS,MVDS,NS0,LAMSMIN,LAMSMAX,KINV,C1X2,VTA1,VTB1,&
              BEST,S3M1D,KDX,TC1D,SBA1,AMS2,AAS2,AMS3,BMS3,AAS3,BAS3,  &
              AEAR,LAMD,DMAX,SMLF,SMLR,BDS

      NS0 = MIN(2.E8,2.E6*EXP(-0.12*MIN(-1.E-3,TK1D-TK0C)))
      IF (ICE_RHOS.EQ.0) THEN
         RHOS = RHOS1
         VS1D = 0.
      ELSEIF (ICE_RHOS.EQ.1) THEN
         RHOS = QS1D/(VS1D+ISMALL)
         IF (TK1D.GE.TK0C) THEN
            IF (RHOS.GT.RHOIMAX) THEN
               RHOS = RHOIMAX
            ENDIF
         ELSE
            IF (RHOS.GT.RHOIMAX) THEN
               RHOS = RHOIMAX
            ENDIF
         ENDIF
         IF (RHOS.LT.RHOIMIN) THEN
            RHOS = RHOIMIN
         ENDIF
         VS1D = QS1D/RHOS
      ELSEIF (ICE_RHOS.EQ.2) THEN
         RHOS = ROSA**(4./(ROSB+4.))*6.**(4.*ROSB/(3.*ROSB+12.))*(PI*  &
                NS0/QS1D)**(-1.*ROSB/(ROSB+4.))
         RHOS = MIN(MAX(RHOS,RHOIMIN),RHOIMAX)
         VS1D = 0.
      ENDIF
      IF (NS1D.LT.NSMALL) THEN
         NS1D = MAX(1.,(QS1D/PI/RHOS)**0.25*NS0**0.75)
      ENDIF
      IF (AGG_SHAPE.EQ.0) THEN
         SASPR = 1.
         FS1D  = 0.
      ELSEIF (AGG_SHAPE.EQ.1) THEN
         SASPR = EXP(-1.*FS1D/QS1D)
         IF (SASPR.LT.SASMIN) THEN
            SASPR = SASMIN
            FS1D  = -1.*LOG(SASPR)*QS1D
         ELSEIF (SASPR.GT.SASMAX) THEN
            SASPR = SASMAX
            FS1D  = -1.*LOG(SASPR)*QS1D
         ENDIF
      ENDIF
      AMS   = PI*RHOS/6.
      S3M1D = QS1D/AMS
      IF (S2M1D.GE.ASMALL.AND.S3M1D.GE.ISMALL) THEN
         KDX = MAX(0.,S2M1D**3./(S3M1D**2.*NS1D))
         IF (KDX.GT.KCSMAX) THEN
            KDX   = KCSMIN
            S2M1D = (KDX*NS1D*S3M1D**2.)**THRD
         ELSEIF (KDX.LT.KCSMIN) THEN
            KDX   = KCSMIN
            S2M1D = (KDX*NS1D*S3M1D**2.)**THRD
         ENDIF
         AFAS = (6.*KDX-3.+SQRT(8.*KDX+1.))/(2.-2.*KDX)
         AFAS = MIN(MAX(AFAS,AFAMIN),AFAMAX)
!         LAMS = (NS1D*(AFAS+1.)*(AFAS+2.)*(AFAS+3.)/S3M1D)**THRD
         LAMS = (AFAS+3.)*S2M1D/S3M1D
      ELSE
         IF (AFAS_3M.EQ.0) THEN
            AFAS  = AFAS0
            S2M1D = 0.
         ELSEIF (AFAS_3M.EQ.1) THEN
            AFAS  = AFAS0
            S2M1D = (KCSMIN*NS1D*S3M1D**2.)**THRD
         ELSEIF (AFAS_3M.EQ.2) THEN
            BDS   = MIN(MAX((S3M1D/NS1D)**THRD,DSMIN),DSMAX)
            AFAS  = MAX(AFAMIN,4.5*TANH(0.5*(1.E3*BDS-5.))+5.5)
            S2M1D = 0.
         ENDIF
         LAMS = (EXP(GAMLN(AFAS+4.)-GAMLN(AFAS+1.))*NS1D/S3M1D)**THRD
      ENDIF
      LAMSMIN = (EXP(GAMLN(AFAS+4.)-GAMLN(AFAS+1.)))**THRD/DSMAX
      LAMSMAX = (EXP(GAMLN(AFAS+4.)-GAMLN(AFAS+1.)))**THRD/DSMIN
      IF (LAMS.LT.LAMSMIN) THEN
         LAMS = LAMSMIN
         NS1D = S3M1D*EXP(GAMLN(AFAS+1.)-GAMLN(AFAS+4.)+3.*LOG(LAMS))
      ELSEIF (LAMS.GT.LAMSMAX) THEN
         LAMS = LAMSMAX
         NS1D = S3M1D*EXP(GAMLN(AFAS+1.)-GAMLN(AFAS+4.)+3.*LOG(LAMS))
      ENDIF
      MVDS = (EXP(GAMLN(AFAS+4.)-GAMLN(AFAS+1.)))**THRD/LAMS
      IF (ICE_VTS.EQ.0) THEN
         AVS = AVS0
         BVS = BVS0
      ELSEIF (ICE_VTS.EQ.1.OR.ICE_VTS.EQ.2) THEN
         KINV = (1.72E-5*(393./(TK1D+120.))*(TK1D/TK0C)**1.5)/RHO
         IF (AGG_SHAPE.EQ.0) THEN
            AMS3  = AMS
            BMS3  = BMS
            AAS3  = PI/4.
            BAS3  = 2.
            LAMD  = LAMS
         ELSEIF (AGG_SHAPE.EQ.1) THEN
            TC1D  = TK1D-TK0C
            IF (ICE_VTS.EQ.1) THEN
               TBIN = MIN(MAX(INT(ABS(TC1D)/5.)-4,0),8)
               IF (MVDS.LT.1.5E-4) DBIN = 0
               IF (MVDS.GE.1.5E-4.AND.MVDS.LT.3.E-4) DBIN = 1
               IF (MVDS.GE.3.E-4.AND.MVDS.LT.1.E-3) DBIN = 2
               IF (MVDS.GE.1.E-3.AND.MVDS.LT.3.E-3) DBIN = 3
               IF (MVDS.GE.3.E-3) DBIN = 4
               DBIN = MIN(MAX(DBIN,0),4)
               AMS3 = AMS1(TBIN*5+DBIN)
               BMS3 = BMS1(TBIN*5+DBIN)
               AAS3 = AAS1(TBIN*5+DBIN)
               BAS3 = BAS1(TBIN*5+DBIN)
            ELSEIF (ICE_VTS.EQ.2) THEN
               BMS3 = 2.4+0.0085*MAX(TC1D,-65.)
               BAS3 = 2-0.19+0.0056*MAX(TC1D,-65.)
               AMS2 = 0.0102+0.00013*MAX(TC1D,-65.)
               AAS2 = 0.29+0.0035*MAX(TC1D,-65.)
!               BMS3 = 2.1+0.0036*MAX(TC1D,-65.)
!               BAS3 = 2-0.35+0.0026*MAX(TC1D,-65.)
!               AMS2 = 0.0064+0.000095*MAX(TC1D,-65.)
!               AAS2 = 0.21+0.0024*MAX(TC1D,-65.)
               AMS3 = ((AMS2*(100.*MVDS)**BMS3)/1.E3)/(MVDS**BMS3)
               AAS3 = ((AAS2*(100.*MVDS)**BAS3)/1.E4)/(MVDS**BAS3)
            ENDIF
            LAMD = MIN((AMS3*LAMS**3./AMS*EXP(GAMLN(BMS3+AFAS+1.)-     &
                   GAMLN(AFAS+4.)))**(1./BMS3),LAMS)
         ENDIF
         SBA1 = BMS3+2.-BAS3
         BEST = 2.*G*NS1D*AMS3*EXP(GAMLN(SBA1+AFAS+1.)-GAMLN(AFAS+1.)- &
                SBA1*LOG(LAMD))/(KINV**2.*AAS3)
         C1X2 = VTC1*BEST**5.E-1
         VTB1 = C1X2/(1.+C1X2)**5.E-1/((1.+C1X2)**5.E-1-1.)/2.
         VTA1 = VTC2*((1+C1X2)**5.E-1-1.)**2./BEST**VTB1
         AVS  = VTA1*KINV**(1.-2.*VTB1)*(2.*AMS3*G/AAS3)**VTB1
         BVS  = VTB1*(BMS3+2.-BAS3)-1.
         IF (TK1D.GE.TK0C.AND.SMLF.GE.1.E-8) THEN
            SMLR = 0.01195*EXP(4.411*SMLF)
            AVS  = 88.701591*EXP(SMLF*2.495827)
            BVS  = -0.182754*SMLF**2.+0.368609*SMLF+0.608571
         ELSE
            AVS  = AVS
            BVS  = BVS
         ENDIF
      ENDIF

      END SUBROUTINE SOLVE_AFAS
!======================================================================

!======================================================================
      SUBROUTINE SOLVE_AFAG(TK1D,RHO,QG1D,NG1D,VG1D,G2M1D,LAMG,AFAG,   &
                 MVDG,RHOG,AMG,AVG,BVG,GMLF)
!======================================================================
      IMPLICIT NONE
      REAL :: TK1D,RHO,QG1D,NG1D,VG1D,G2M1D,G3M1D,RHOG,LAMG,AFAG,AMG,  &
              AVG,BVG,BEST0,LAMGMIN,LAMGMAX,KDX,KINV,GG1,C1X2,VTA1,    &
              VTB1,BEST,MVDG,GMLF,GMLR,BDG

      IF (ICE_RHOG.EQ.0) THEN
         RHOG = RHOG1
         VG1D = 0.
      ELSEIF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
         RHOG = QG1D/(VG1D+ISMALL)
         IF (TK1D.GE.TK0C) THEN
            IF (RHOG.GT.RHOH) THEN
               RHOG = RHOH
            ENDIF
         ELSE
            IF (RHOG.GT.RHOG0) THEN
               RHOG = RHOG0
            ENDIF
         ENDIF
         IF (RHOG.LT.RHOIMIN) THEN
            RHOG = RHOIMIN
         ENDIF
         VG1D = QG1D/RHOG
      ENDIF
      IF (NG1D.LT.NSMALL) THEN
         NG1D = MAX(1.,(QG1D/PI/RHOG)**0.25*NG0**0.75)
      ENDIF
      AMG   = PI*RHOG/6.
      G3M1D = QG1D*V2M3/RHOG
      IF (G2M1D.GE.ASMALL.AND.G3M1D.GE.ISMALL) THEN
         KDX = MAX(0.,G2M1D**3./(G3M1D**2.*NG1D))
         IF (KDX.GT.KCGMAX) THEN
            KDX   = KCGMIN
            G2M1D = (KDX*NG1D*G3M1D**2.)**THRD
         ELSEIF (KDX.LT.KCGMIN) THEN
            KDX   = KCGMIN
            G2M1D = (KDX*NG1D*G3M1D**2.)**THRD
         ENDIF
         AFAG = (6.*KDX-3.+SQRT(8.*KDX+1.))/(2.-2.*KDX)
         AFAG = MIN(MAX(AFAG,AFAMIN),AFAMAX)
!         LAMG = (NG1D*(AFAG+1.)*(AFAG+2.)*(AFAG+3.)/G3M1D)**THRD
         LAMG = (AFAG+3.)*G2M1D/G3M1D
      ELSE
         IF (AFAG_3M.EQ.0) THEN
            AFAG  = AFAG0
            G2M1D = 0.
         ELSEIF (AFAG_3M.EQ.1) THEN
            AFAG  = AFAG0
            G2M1D = (KCGMIN*NG1D*G3M1D**2.)**THRD
         ELSEIF (AFAG_3M.EQ.2) THEN
            BDG   = MIN(MAX((G3M1D/NG1D)**THRD,DGMIN),DGMAX)
            AFAG  = MAX(AFAMIN,5.5*TANH(0.7*(1.E3*BDG-4.5))+8.5)
            G2M1D = 0.
         ENDIF
         LAMG = (EXP(GAMLN(AFAG+4.)-GAMLN(AFAG+1.))*NG1D/G3M1D)**THRD
      ENDIF
      LAMGMIN = (EXP(GAMLN(AFAG+4.)-GAMLN(AFAG+1.)))**THRD/DGMAX
      LAMGMAX = (EXP(GAMLN(AFAG+4.)-GAMLN(AFAG+1.)))**THRD/DGMIN
      IF (LAMG.LT.LAMGMIN) THEN
         LAMG = LAMGMIN
         NG1D = G3M1D*EXP(GAMLN(AFAG+1.)-GAMLN(AFAG+4.)+3.*LOG(LAMG))
      ELSEIF (LAMG.GT.LAMGMAX) THEN
         LAMG = LAMGMAX
         NG1D = G3M1D*EXP(GAMLN(AFAG+1.)-GAMLN(AFAG+4.)+3.*LOG(LAMG))
      ENDIF
      MVDG = (EXP(GAMLN(AFAG+4.)-GAMLN(AFAG+1.)))**THRD/LAMG
      IF (ICE_VTG.EQ.0) THEN
         AVG  = AVG0
         BVG  = BVG0
      ELSEIF (ICE_VTG.EQ.1) THEN
         KINV  = (1.72E-5*(393./(TK1D+120.))*(TK1D/TK0C)**1.5)/RHO
         BEST0 = 2.*G*NG1D/(KINV**2.)
         GG1   = GAMLN(AFAG+1.)
         BEST  = BEST0*AMG*EXP(GAMLN(BMG+AFAG+1.)-GG1-BMG*LOG(LAMG))/AAW
         C1X2  = VTC1*BEST**5.E-1
         VTB1  = C1X2/(1.+C1X2)**5.E-1/((1.+C1X2)**5.E-1-1.)/2.
         VTA1  = VTC2*((1+C1X2)**5.E-1-1.)**2./BEST**VTB1
         AVG   = VTA1*KINV**(1.-2.*VTB1)*(2.*AMG*G/AAW)**VTB1
         BVG   = VTB1*BMG-1.
      ENDIF

      END SUBROUTINE SOLVE_AFAG
!======================================================================
!
!======================================================================
      SUBROUTINE SOLVE_AFAH(TK1D,RHO,QH1D,NH1D,H2M1D,LAMH,AFAH,MVDH,   &
                 AVH,BVH)
!======================================================================
      IMPLICIT NONE
      REAL :: TK1D,RHO,QH1D,NH1D,H2M1D,H3M1D,LAMH,AFAH,AVH,BVH,BEST0,  &
              GH1,KDX,MVDH,LAMHMIN,LAMHMAX,KINV,C1X2,VTA1,VTB1,BEST,BDH

      IF (NH1D.LT.NSMALL) THEN
         NH1D = MAX(1.,(QH1D/PI/RHOH)**0.25*NH0**0.75)
      ENDIF
      H3M1D = QH1D*V2M3/RHOH
      IF (H2M1D.GE.ASMALL.AND.H3M1D.GE.ISMALL) THEN
         KDX = MAX(0.,H2M1D**3./(H3M1D**2.*NH1D))
         IF (KDX.GT.KCHMAX) THEN
            KDX   = KCHMIN
            H2M1D = (KDX*NH1D*H3M1D**2.)**THRD
         ELSEIF (KDX.LT.KCHMIN) THEN
            KDX   = KCHMIN
            H2M1D = (KDX*NH1D*H3M1D**2.)**THRD
         ENDIF
         AFAH = (6.*KDX-3.+SQRT(8.*KDX+1.))/(2.-2.*KDX)
         AFAH = MIN(MAX(AFAH,AFAMIN),AFAMAX)
!         LAMH = (NH1D*(AFAH+1.)*(AFAH+2.)*(AFAH+3.)/H3M1D)**THRD
         LAMH = (AFAH+3.)*H2M1D/H3M1D
      ELSE
         IF (AFAH_3M.EQ.0) THEN
            AFAH  = AFAH0
            H2M1D = 0.
         ELSEIF (AFAH_3M.EQ.1) THEN
            AFAH  = AFAH0
            H2M1D = (KCHMIN*NH1D*H3M1D**2.)**THRD
         ELSEIF (AFAH_3M.EQ.2) THEN
            BDH   = MIN(MAX((H3M1D/NH1D)**THRD,DHMIN),DHMAX)
            AFAH  = MAX(AFAMIN,3.7*TANH(0.3*(1.E3*BDH-9.))+6.5)
            H2M1D = 0.
         ENDIF
         LAMH = (EXP(GAMLN(AFAH+4.)-GAMLN(AFAH+1.))*NH1D/H3M1D)**THRD
      ENDIF
      LAMHMIN = (EXP(GAMLN(AFAH+4.)-GAMLN(AFAH+1.)))**THRD/DHMAX
      LAMHMAX = (EXP(GAMLN(AFAH+4.)-GAMLN(AFAH+1.)))**THRD/DHMIN
      IF (LAMH.LT.LAMHMIN) THEN
         LAMH = LAMHMIN
         NH1D = H3M1D*EXP(GAMLN(AFAH+1.)-GAMLN(AFAH+4.)+3.*LOG(LAMH))
      ELSEIF (LAMH.GT.LAMHMAX) THEN
         LAMH = LAMHMAX
         NH1D = H3M1D*EXP(GAMLN(AFAH+1.)-GAMLN(AFAH+4.)+3.*LOG(LAMH))
      ENDIF
      MVDH = (EXP(GAMLN(AFAH+4.)-GAMLN(AFAH+1.)))**THRD/LAMH
      IF (ICE_VTH.EQ.0) THEN
         AVH  = AVH0
         BVH  = BVH0
      ELSEIF (ICE_VTH.EQ.1) THEN
         KINV  = (1.72E-5*(393./(TK1D+120.))*(TK1D/TK0C)**1.5)/RHO
         BEST0 = 2.*G*NH1D/(KINV**2.)
         GH1   = GAMLN(AFAH+1.)
         BEST  = BEST0*AMH*EXP(GAMLN(BMH+AFAH+1.)-GH1-BMH*LOG(LAMH))/AAW
         C1X2  = VTC1*BEST**5.E-1
         VTB1  = C1X2/(1.+C1X2)**5.E-1/((1.+C1X2)**5.E-1-1.)/2.
         VTA1  = VTC2*((1+C1X2)**5.E-1-1.)**2./BEST**VTB1
         AVH   = VTA1*KINV**(1.-2.*VTB1)*(2.*AMH*G/AAW)**VTB1
         BVH   = VTB1*BMH-1.
      ENDIF

      END SUBROUTINE SOLVE_AFAH
!======================================================================
!
!======================================================================
      SUBROUTINE REFL10CM_HM(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,t1d,   &
                 p1d,dBZ,kts,kte,ii,jj)
!======================================================================
      IMPLICIT NONE
!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, DIMENSION(kts:kte), INTENT(IN) :: qv1d,qr1d,nr1d,qs1d,ns1d,&
                                              qg1d,ng1d,t1d,p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT) :: dBZ
!..Local variables
      REAL, DIMENSION(kts:kte) :: temp,pres,qv,rho,rr,nr,rs,ns,rg,ng
      DOUBLE PRECISION, DIMENSION(kts:kte) :: ilamr,ilamg,ilams
      DOUBLE PRECISION, DIMENSION(kts:kte) :: N0_r,N0_g,N0_s
      DOUBLE PRECISION :: lamr,lamg,lams
      LOGICAL, DIMENSION(kts:kte):: L_qr,L_qs,L_qg
      REAL, DIMENSION(kts:kte) :: ze_rain,ze_snow,ze_graupel
      DOUBLE PRECISION :: fmelt_s,fmelt_g
      DOUBLE PRECISION :: cback,x,eta,f_d
      INTEGER:: i,k,k_0,kbot,n
      LOGICAL:: melti
!+---+
      do k = kts,kte
         dBZ(k) = -35.0
      enddo
!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts,kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10,qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         if (qr1d(k).gt.1.E-9) then
            rr(k) = qr1d(k)*rho(k)
            nr(k) = nr1d(k)*rho(k)
            lamr = (xam_r*xcrg(3)*xorg2*nr(k)/rr(k))**xobmr
            ilamr(k) = 1./lamr
            N0_r(k) = nr(k)*xorg2*lamr**xcre(2)
            L_qr(k) = .true.
         else
            rr(k) = 1.E-12
            nr(k) = 1.E-12
            L_qr(k) = .false.
         endif
         if (qs1d(k).gt.1.E-9) then
            rs(k) = qs1d(k)*rho(k)
            ns(k) = ns1d(k)*rho(k)
            lams = (xam_s*xcsg(3)*xosg2*ns(k)/rs(k))**xobms
            ilams(k) = 1./lams
            N0_s(k) = ns(k)*xosg2*lams**xcse(2)
            L_qs(k) = .true.
         else
            rs(k) = 1.E-12
            ns(k) = 1.E-12
            L_qs(k) = .false.
         endif
         if (qg1d(k).gt.1.E-9) then
            rg(k) = qg1d(k)*rho(k)
            ng(k) = ng1d(k)*rho(k)
            lamg = (xam_g*xcgg(3)*xogg2*ng(k)/rg(k))**xobmg
            ilamg(k) = 1./lamg
            N0_g(k) = ng(k)*xogg2*lamg**xcge(2)
            L_qg(k) = .true.
         else
            rg(k) = 1.E-12
            ng(k) = 1.E-12
            L_qg(k) = .false.
         endif
      enddo
!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      melti = .false.
      k_0 = kts
      do k = kte-1,kts,-1
         if ((temp(k).gt.273.15).and.L_qr(k).and.(L_qs(k+1).or.        &
              L_qg(k+1))) then
            k_0 = MAX(k+1,k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue
!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+
      do k = kts,kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93)*(6./PI)*(6./PI)*(xam_s/&
                                   900.)*(xam_s/900.)*N0_s(k)*xcsg(4)* &
                                   ilams(k)**xcse(4)
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93)*(6./PI)*(6./PI)*(   &
                                      xam_g/900.)*(xam_g/900.)*N0_g(k)*&
                                      xcgg(4)*ilamg(k)**xcge(4)
      enddo
!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+
      if (melti.and.k_0.ge.kts+1) then
         do k = k_0-1,kts,-1
!..Reflectivity contributed by melting snow
            if (L_qs(k).and.L_qs(k_0)) then
               fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0),0.99d0))
               eta = 0.d0
               lams = 1./ilams(k)
               do n = 1,nrbins
                  x = xam_s*xxDs(n)**xbm_s
               CALL rayleigh_soak_wetgraupel(x,DBLE(xocms),DBLE(xobms),&
                    fmelt_s,melt_outside_s,m_w_0,m_i_0,lamda_radar,    &
                    CBACK,mixingrulestring_s,matrixstring_s,           &
                    inclusionstring_s,hoststring_s,hostmatrixstring_s, &
                    hostinclusionstring_s)
                  f_d = N0_s(k)*xxDs(n)**xmu_s*DEXP(-lams*xxDs(n))
                  eta = eta+f_d*CBACK*simpson(n)*xdts(n)
               enddo
               ze_snow(k) = SNGL(lamda4/(pi5*K_w)*eta)
            endif
!..Reflectivity contributed by melting graupel
            if (L_qg(k).and.L_qg(k_0)) then
               fmelt_g = MAX(0.005d0,MIN(1.0d0-rg(k)/rg(k_0),0.99d0))
               eta = 0.d0
               lamg = 1./ilamg(k)
               do n = 1,nrbins
                  x = xam_g*xxDg(n)**xbm_g
               CALL rayleigh_soak_wetgraupel(x,DBLE(xocmg),DBLE(xobmg),&
                    fmelt_g,melt_outside_g,m_w_0,m_i_0,lamda_radar,    &
                    CBACK,mixingrulestring_g,matrixstring_g,           &
                    inclusionstring_g,hoststring_g,hostmatrixstring_g, &
                    hostinclusionstring_g)
                  f_d = N0_g(k)*xxDg(n)**xmu_g*DEXP(-lamg*xxDg(n))
                  eta = eta+f_d*CBACK*simpson(n)*xdtg(n)
               enddo
               ze_graupel(k) = SNGL(lamda4/(pi5*K_w)*eta)
            endif
         enddo
      endif
      do k = kte,kts,-1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo

      END SUBROUTINE REFL10CM_HM
!======================================================================
!
!======================================================================
      SUBROUTINE AERO_CONST(CCNTY)
!======================================================================
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: CCNTY
      INTEGER, PARAMETER :: ID_SENS = 1                                 ! sensetivity test ID
      INTEGER, PARAMETER :: N1 = NCCN, N2 = 5                           ! AEROSOL MODE AND COMPONENT
      INTEGER :: IA,I,K,J,IM,IV,IB
      INTEGER, DIMENSION(NAER) :: IDAER
      REAL, DIMENSION(NAER) :: XMASS
      REAL, DIMENSION(N2) :: ASHX,CMASX,DNASX,AVANX                     ! AEROSOL SCALE HEIGHT. 3.57km-continental, 800m-marine; MOLAR WEIGHT OF AEROSOL (kg/mol); DENSITY OF SOLID AEROSOL (kg/m^3); VAN'T HOFF FACTOR
      REAL, DIMENSION(N1,N2) :: ZCCNX,CNMODX,CNSTDX                     !
      REAL, DIMENSION(N2,3) :: SENX                                     ! SENSITIVE TEST RATIO
      REAL :: XLOG3,XMTOT,XMAS,DUST_IN0,DUST_IN,DAIR,DZ0,ZH             ! Height of grid point [m]
      REAL*8 :: DMODE,D2STDV,DLNXX,DERF,DSQRT,DLOG
      DATA ZCCNX/3.4E8,6.E7,3.E6, 1.E9,8.E8,7.2E5, 6.4E9,2.3E9,3.2E6,  &! 1. MARINE (WHITBY,1978); 2. CLEAN CONTINENTAL (WHITBY,1978); 3. AVERAGE CONTINENTAL (WHITBY,1978); 4. URBAN (WHITBY,1978); 5. ICE NUECEI
                 1.06E11,3.2E10,5.4E6, 4.02E-11,4.5E9,5.4E6 /           ! 2. CLEAN CONTINENTAL (WHITBY,1978)
      DATA CNMODX/5.E-9,3.55E-8,3.1E-7,8.E-9,3.3E-8,4.6E-7,7.5E-9,     &
                  3.8E-8,5.1E-7,7.E-9,2.7E-8,4.3E-7,4.5E-9,2.35E-8,5.5E-7/
      DATA CNSTDX/0.47,0.69,0.99,0.47,0.74,0.79,0.53,0.69,0.77,0.59,   &
                  0.77,0.79,0.53,0.75,0.74/
      DATA ASHX/8.E2,3.57E3,3.57E3,2.E3,1.E4/
      DATA CMASX/0.13214,0.13214,0.13214,0.13214,0.13214/               ! MOLAR WEIGHT OF AEROSOL (kg/mol)
      DATA DNASX/1769.,1769.,1769.,1769.,1769./                         ! DENSITY OF SOLID AEROSOL (kg/m^3)
      DATA AVANX/3.,3.,3.,3.,3./                                        ! VAN'T HOFF FACTOR
      DATA SENX/1.,1.,1.,1.,1.,10.,10.,10.,10.,10.,0.1,0.1,0.1,0.1,0.1/ ! 1st, 2nd, and 3rd set for sensetivity test

      IF (CCNTY.EQ.1) THEN
         NC0 = 5.E7
      ELSEIF (CCNTY.EQ.2) THEN
         NC0 = 1.E8
      ELSEIF (CCNTY.EQ.3) THEN
         NC0 = 2.5E8
      ELSEIF (CCNTY.EQ.4) THEN
         NC0 = 5.E8
      ENDIF

      IBAER(1) = 1
      DO I = 2,NAER
         IBAER(I) = IBAER(I-1)+NAERN(I-1)
      ENDDO
      TBLRC(1)     = DLOG(1.D-9)                                        ! 0.001 um
      TBLRC(NTBXA) = DLOG(1.D-3)                                        ! 1000. um
      DO I = 2,NTBXA-1
         TBLRC(I) = TBLRC(1)+DBLE(I-1)/DBLE(NTBXA-1)*(TBLRC(NTBXA)-    &
                    TBLRC(1))                                           ! ln(rc)
      ENDDO
      IDAER = (/CCNTY,5/)
      PRINT *, 'AEROSOL ID =',IDAER
      DO J = 1,NAER
         ASH(J)   = ASHX(IDAER(J))
         DNAS(J)  = DNASX(IDAER(J))
         CMAS(J)  = CMASX(IDAER(J))
         AVAN(J)  = AVANX(IDAER(J))
         BETA1(J) = CMW*DNAS(J)*AVAN(J)/(CMAS(J)*RHOW)
         XMTOT = 0.
         DO IM = 1,N1
            ZCCN(IM,J)  = ZCCNX(IM,IDAER(J))*1.                         ! change the factor to change the initial concentration
            CNMOD(IM,J) = CNMODX(IM,IDAER(J))
            CNSTD(IM,J) = CNSTDX(IM,IDAER(J))
            XLOG3       = CNMOD(IM,J)*CNMOD(IM,J)*CNMOD(IM,J)
            WMAS(IM,J)  = ZCCN(IM,J)*XLOG3*EXP(4.5*CNSTD(IM,J)**2.)     ! factor for mass, the 3rd momentum
            RFACT(IM,J) = 1./(C4PI3*DNAS(J)*XLOG3*EXP(4.5*             &
                          CNSTD(IM,J)**2.))
            XMTOT = XMTOT+WMAS(IM,J)
         ENDDO
         DO IM = 1,N1
            WMAS(IM,J) = WMAS(IM,J)/XMTOT
         ENDDO

         DO I = 1,NTBXA
            TBLXA(I,J) = DBLE(0.)
         ENDDO
         DO IM = 1,N1
            DMODE = DLOG(DBLE(CNMOD(IM,J)))+DBLE(3.)*DBLE(CNSTD(IM,J))*DBLE(CNSTD(IM,J))
            D2STDV = DSQRT(DBLE(2.))*DBLE(CNSTD(IM,J))
            DO I = 1,NTBXA
               DLNXX = (TBLRC(I)-DMODE)/D2STDV
               TBLXA(I,J) = TBLXA(I,J)+DBLE(WMAS(IM,J))*DBLE(0.5)*     &
                            (1.D+0+DERF(DLNXX))
            ENDDO
         ENDDO                                                          ! IM=1,N1
         DO I = 1,NTBXA
            TBLXA(I,J) = DMAX1(1.D-20,TBLXA(I,J))
         ENDDO

         CALL FIND_RC0(1.D-6,CNMOD(1,J),CNSTD(1,J),WMAS(1,J),RXMIN(J), &
                       TBLXA(1,J),TBLRC)                                ! lower limit of dry aerosol fraction is set as 1.E-6
         WRITE(*,*)'FOR DRY AEROSOL FRACTION BE > 1.E-6'
         WRITE(*,'(A30,I1,A8,F6.4,A3)')'MINIMUN R-CUTOFF FOR AEROSOL'  &
               ,J,' IS SET ',RXMIN(J)*1.E6,' um'
         WRITE(*,*)'SENSITIVITY TEST ID:',ID_SENS
!----------------- FACTOR FOR SENSETIVITY ------------------------------
         SENS(J) = SENX(IDAER(J),ID_SENS)
      ENDDO

      END SUBROUTINE AERO_CONST
!======================================================================

!======================================================================
      SUBROUTINE INIT_AEROSOL(P,RHO,DZ,XLAND,QAERO,BACKIN,IDS,IDE,JDS, &
                 JDE,KDS,KDE,IMS,IME,JMS,JME,KMS,KME,ITS,ITE,JTS,JTE,  &
                 KTS,KTE)
!======================================================================
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE,IMS,IME,JMS,JME,  &
                             KMS,KME,ITS,ITE,JTS,JTE,KTS,KTE
      INTEGER :: ITF,JTF,IA,I,K,J,IM,K_PBL,IV,IB
      INTEGER, DIMENSION(NAER) :: IDAER
      REAL, INTENT(IN) :: BACKIN
      REAL, INTENT(IN), DIMENSION(IMS:IME,JMS:JME) :: XLAND             ! 1:land ; 2:ocean
      REAL, INTENT(IN), DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: P,RHO,DZ
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE,NAERT) :: QAERO
      REAL, DIMENSION(ITS:ITE,JTS:JTE,NAER) :: RASH
      REAL, DIMENSION(NAER) :: XMASS
      REAL, PARAMETER :: P_PBL = 85000.                                 ! set PBL pressure [Pa]
      REAL :: XLOG3,XMTOT,XMAS,DUST_IN0,DUST_IN,DAIR,DZ0,ZH             ! Height of grid point [m]
      REAL*8 :: DMODE,D2STDV,DLNXX,DERF,DSQRT,DLOG

      ITF = MIN(ITE,IDE-1)
      JTF = MIN(JTE,JDE-1)
      DO IA = 1,NAERT
         DO K = KTS,KTE
         DO J = JTS,JTF
         DO I = ITS,ITF
            QAERO(I,K,J,IA) = 0.
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      DO J = 1,NAER
         XMASS(J) = 0.
         DO IM = 1,NCCN
            XMAS = ZCCN(IM,J)*CNMOD(IM,J)**3.*EXP(4.5*CNSTD(IM,J)**2.)
            XMASS(J) = XMASS(J)+XMAS
         ENDDO
         XMASS(J) = XMASS(J)*C4PI3*DNAS(J)*SENS(J)
      ENDDO
      ! CC: code commented below nver happens and it fails strict
      ! compilation
      !IF (ID_DUST*ID_IN.NE.0) THEN                                      ! IF DUST EXISTS
      !   DUST_IN0 = 0.
      !   DO IM = 1,NCCN
      !      DUST_IN0 = DUST_IN0+ZCCN(IM,ID_DUST)
      !   ENDDO
      !ELSE                                                              ! NO DUST
         DUST_IN0 = BACKIN
      !ENDIF
      DO J = JTS,JTF
         DO I = ITS,ITF
            K_PBL = 1
            DO K = KTS,KTE-1
               IV = 1
               DO IA = 1,NAER
                  QAERO(I,K,J,IV) = XMASS(IA)/RHO(I,K,J)                ! dry aerosol [kg/kg]
                  RASH (I,J,IA)   = 1./ASH(IA)                          ! increase scale height if necessary
                  DUST_IN         = DUST_IN0                            ! dust # over land
                  IF (IA.EQ.ID_IN) THEN                                 ! IN
                     IF (ID_DUST.NE.0) THEN                             ! DUST IS CONSIDERED AS IN
!                        IF (XLAND(I,J).EQ.2.) THEN                      ! over ocean
!                           DUST_IN       = DUST_IN*0.01                 ! dust #
!                           RASH (I,J,IA) = RASH (I,J,IA)*0.5            ! double scale height
!                        ENDIF
                     ELSE
                        DUST_IN       = BACKIN                          ! set to background value
                        RASH (I,J,IA) = RASH (I,J,IA)*0.5               ! double scale height
                     ENDIF
                     DAIR    = RHO(I,K,J)
                     DUST_IN = MAX(BACKIN,DUST_IN)
                     QAERO(I,K,J,IV) = DUST_IN/DAIR                     ! 0.5 #/L [#/kg]*[kPa]      ! tsai
                  ELSE                                                  ! dust, NH42SO4, chemicals....
!                     IF (XLAND(I,J).EQ.2.) THEN                         ! over ocean
!                        IF (IA.EQ.ID_NH42SO4) THEN                      !
!                           QAERO(I,K,J,IV) = QAERO(I,K,J,IV)*0.05       ! reduce it
!                           RASH(I,J,IA)    = RASH(I,J,IA)*0.5           ! double scale height
!                        ELSE                                            ! dust or IN, not sulfate (the same)
!                           QAERO(I,K,J,IV) = QAERO(I,K,J,IV)*0.01       ! reduce it
!                           RASH(I,J,IA)    = RASH (I,J,IA)*0.5          ! double scale height
!                           DUST_IN         = DUST_IN0*0.01              ! dust #
!                        ENDIF
!                     ENDIF
                  ENDIF
                  IV = IV+NAERN(IA)
               ENDDO
               IF (P(I,K,J).LT.P_PBL) K_PBL = K
            ENDDO
            K_PBL = MIN(K_PBL,KTE-3+1)                                  ! at least 3 layers are within PBL
            ZH  = 0.
            DZ0 = 0.
            DO K = K_PBL-1,1,-1
               DAIR = RHO(I,K,J)
               ZH   = ZH+(DZ0+DZ(I,K,J))*0.5
               DZ0  = DZ(I,K,J)
               IV = 1
               DO IA = 1,NAER
                  QAERO(I,K,J,IV)=QAERO(I,K,J,IV)*EXP(-RASH(I,J,IA)*ZH) ! dry aerosol [kg/kg]*[kPa]
                  IV = IV+NAERN(IA)
               ENDDO
            ENDDO
            IV = 0
            DO IA = 1,NAER
               DO IB = 2,MIN(2,NAERN(IA))                               !  if NV>2
                  DO K = KTS,KTE
                     QAERO(I,K,J,IV+IB) = QAERO(I,K,J,IV+1)             ! set NV2: total = dry aerosol
                  ENDDO
               ENDDO
               IV = IV+NAERN(IA)
            ENDDO
         ENDDO                                                          ! FOR I LOOPS
      ENDDO                                                             ! FOR J LOOPS

      PRINT *,'IN init_aerosol.F,I,J,NAER',(ITS+ITE)/2,(JTS+JTE)/2,NAER
      DO K = KTS,KTE
         WRITE(*,'(I2,X,20(E12.6,X))') K,(QAERO((ITS+ITE)/2,K,         &
         (JTS+JTE)/2,IV),IV=1,NAERT)
      ENDDO

      RETURN

      END SUBROUTINE INIT_AEROSOL
!======================================================================
!
!=======================================================================
      SUBROUTINE MP_NTU(ITIMESTEP,TH,P,DZ,W,PII,DT_MP,SR,QV,QC,QR,QI,  &
                 QS,QG,QH,NC,NR,NI,NS,NG,NH,VI,VS,VG,FI,FS,AI,AS,AG,AH,&
                 I3M,QDCN,QTCN,QCCN,QRCN,QNIN,IDEPNU,RAINNC,RAINNCV,   &
                 SNOWNC,SNOWNCV,GRAPNC,GRAPNCV,HAILNC,HAILNCV,IDS,IDE, &
                 JDS,JDE,KDS,KDE,IMS,IME,JMS,JME,KMS,KME,ITS,ITE,JTS,  &
                 JTE,KTS,KTE)
!======================================================================
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE,IMS,IME,JMS,JME,  &
            KMS,KME,ITS,ITE,JTS,JTE,KTS,KTE,IDEPNU,ITIMESTEP
      REAL, INTENT(IN) :: DT_MP
      REAL, DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: RAINNC,RAINNCV,&
            SNOWNC,SNOWNCV,GRAPNC,GRAPNCV,HAILNC,HAILNCV,SR
      REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(IN) :: PII,DZ,W,P
      REAL, DIMENSION(IMS:IME,KMS:KME,JMS:JME), INTENT(INOUT) ::       &
            TH,QV,QC,QR,QI,QS,QG,QH,NC,NR,NI,NS,NG,NH,VI,VS,VG,FI,FS,  &
            AI,AS,AG,AH,I3M,QDCN,QTCN,QCCN,QRCN,QNIN
      INTEGER :: I,K,J,NK,ITF,JTF
      REAL :: DT,DTMN
      REAL, DIMENSION(KTS:KTE) :: TK1D,P1D,W1D,S1D,DZ1D,QV1D,QC1D,QR1D,&
            QI1D,QS1D,QG1D,QH1D,NC1D,NR1D,NI1D,NS1D,NG1D,NH1D,VI1D,    &
            VS1D,VG1D,FI1D,FS1D,AI1D,AS1D,AG1D,AH1D,I3M1D,QDCN1D,      &
            QTCN1D,QCCN1D,QRCN1D,QNIN1D
      REAL, DIMENSION(ITS:ITE,JTS:JTE) :: CLODNCV,ICENCV,VINCV,VSNCV,  &
            VGNCV,FINCV,FSNCV,AINCV,ASNCV,AGNCV,AHNCV,I3MNCV
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: TK

      DT = DT_MP
      DTMN = DT_MP/60.
      ITF = MIN(ITE,IDE-1)
      JTF = MIN(JTE,JDE-1)
      DO I = ITS,ITF
      DO J = JTS,JTF
         CLODNCV(I,J) = 0.; RAINNCV(I,J) = 0.; ICENCV(I,J) = 0.
         SNOWNCV(I,J) = 0.; GRAPNCV(I,J) = 0.; HAILNCV(I,J) = 0.
         VINCV(I,J) = 0.; VSNCV(I,J) = 0.; VGNCV(I,J) = 0.
         FINCV(I,J) = 0.; FSNCV(I,J) = 0.; AINCV(I,J) = 0.
         ASNCV(I,J) = 0.; SR(I,J) = 0.; AGNCV(I,J) = 0.; AHNCV(I,J) = 0.
         I3MNCV(I,J) = 0.
!----------------- write data from 3D to 1D and assign upside down -----
         DO K = KTS,KTE
            NK = KME-K
            S1D(K) = 0.
            TK1D(K) = TH(I,NK,J)*PII(I,NK,J)
            W1D(K)  = 0.5*(W(I,NK,J)+W(I,NK+1,J)) ! W at half level
            QV1D(K) = QV(I,NK,J)
            DZ1D(K) = DZ(I,NK,J)
            P1D(K)  = P(I,NK,J)
            QC1D(K) = QC(I,NK,J)
            QR1D(K) = QR(I,NK,J)
            QI1D(K) = QI(I,NK,J)
            QS1D(K) = QS(I,NK,J)
            QG1D(K) = QG(I,NK,J)
            QH1D(K) = QH(I,NK,J)
            NC1D(K) = NC(I,NK,J)
            NR1D(K) = NR(I,NK,J)
            NI1D(K) = NI(I,NK,J)
            NS1D(K) = NS(I,NK,J)
            NG1D(K) = NG(I,NK,J)
            NH1D(K) = NH(I,NK,J)
            VI1D(K) = VI(I,NK,J)
            VS1D(K) = VS(I,NK,J)
            VG1D(K) = VG(I,NK,J)
            FI1D(K) = FI(I,NK,J)
            FS1D(K) = FS(I,NK,J)
            AI1D(K) = AI(I,NK,J)
            AS1D(K) = AS(I,NK,J)
            AG1D(K) = AG(I,NK,J)
            AH1D(K) = AH(I,NK,J)
            I3M1D(K) = I3M(I,NK,J)
            QDCN1D(K) = MAX(RLIMIT,QDCN(I,NK,J))
            QTCN1D(K) = MAX(RLIMIT,QTCN(I,NK,J))
            QCCN1D(K) = MAX(RLIMIT,QCCN(I,NK,J))
            QRCN1D(K) = MAX(RLIMIT,QRCN(I,NK,J))
            QNIN1D(K) = MAX(RLIMIT,QNIN(I,NK,J))
         ENDDO

         CALL NTU_MICRO(TK1D,QV1D,DZ1D,P1D,S1D,W1D,SR(I,J),ICENCV(I,J),&
              VINCV(I,J),VSNCV(I,J),VGNCV(I,J),FINCV(I,J),FSNCV(I,J),  &
              AINCV(I,J),ASNCV(I,J),AGNCV(I,J),AHNCV(I,J),I3MNCV(I,J), &
              CLODNCV(I,J),RAINNC(I,J),RAINNCV(I,J),SNOWNC(I,J),       &
              SNOWNCV(I,J),GRAPNC(I,J),GRAPNCV(I,J),HAILNC(I,J),       &
              HAILNCV(I,J),DT,DTMN,QC1D,QR1D,QI1D,QS1D,QG1D,QH1D,NC1D, &
              NR1D,NI1D,NS1D,NG1D,NH1D,VI1D,VS1D,VG1D,FI1D,FS1D,AI1D,  &
              AS1D,AG1D,AH1D,I3M1D,QDCN1D,QTCN1D,QCCN1D,QRCN1D,QNIN1D, &
              IDEPNU,KTS,KTE)

         DO K = KTS,KTE
            NK = KME-K
            TH(I,NK,J) = TK1D(K)/PII(I,NK,J)
            QV(I,NK,J) = MAX(0.,QV1D(K))
            QC(I,NK,J) = MAX(0.,QC1D(K))
            QR(I,NK,J) = MAX(0.,QR1D(K))
            QI(I,NK,J) = MAX(0.,QI1D(K))
            QS(I,NK,J) = MAX(0.,QS1D(K))
            QG(I,NK,J) = MAX(0.,QG1D(K))
            QH(I,NK,J) = MAX(0.,QH1D(K))
            NC(I,NK,J) = MAX(0.,NC1D(K))
            NR(I,NK,J) = MAX(0.,NR1D(K))
            NI(I,NK,J) = MAX(0.,NI1D(K))
            NS(I,NK,J) = MAX(0.,NS1D(K))
            NG(I,NK,J) = MAX(0.,NG1D(K))
            NH(I,NK,J) = MAX(0.,NH1D(K))
            VI(I,NK,J) = MAX(0.,VI1D(K))
            VS(I,NK,J) = MAX(0.,VS1D(K))
            VG(I,NK,J) = MAX(0.,VG1D(K))
            FI(I,NK,J) = MAX(0.,FI1D(K))
            FS(I,NK,J) = MAX(0.,FS1D(K))
            AI(I,NK,J) = MAX(0.,AI1D(K))
            AS(I,NK,J) = MAX(0.,AS1D(K))
            AG(I,NK,J) = MAX(0.,AG1D(K))
            AH(I,NK,J) = MAX(0.,AH1D(K))
            I3M(I,NK,J) = MAX(0.,I3M1D(K))
            QDCN(I,NK,J) = MAX(0.,QDCN1D(K))
            QTCN(I,NK,J) = MAX(0.,QTCN1D(K))
            QCCN(I,NK,J) = MAX(0.,QCCN1D(K))
            QRCN(I,NK,J) = MAX(0.,QRCN1D(K))
            QNIN(I,NK,J) = MAX(0.,QNIN1D(K))
         ENDDO
!         CALL REFL10CM_HM(QV1D,QR1D,NR1D,QS1D,NS1D,QG1D,NG1D,TK1D,P1D, &
!                          DBZ11D,KTS,KTE,I,J)
!         DO K = KTS,KTE
!            NK = KME-K
!            DBZ1(I,NK,J) = MAX(-99.,DBZ11D(K))
!         ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE MP_NTU
!======================================================================
!
!======================================================================
      SUBROUTINE NTU_MICRO(TK3D,QV3D,DZ3D,P3D,S3D,W3D,SR,ICENCV,VINCV, &
                 VSNCV,VGNCV,FINCV,FSNCV,AINCV,ASNCV,AGNCV,AHNCV,      &
                 I3MNCV,CLODNCV,RAINNC,RAINNCV,SNOWNC,SNOWNCV,GRAPNC,  &
                 GRAPNCV,HAILNC,HAILNCV,DT,DTMN,QC3D,QR3D,QI3D,QS3D,   &
                 QG3D,QH3D,NC3D,NR3D,NI3D,NS3D,NG3D,NH3D,VI3D,VS3D,    &
                 VG3D,FI3D,FS3D,AI3D,AS3D,AG3D,AH3D,I3M3D,QDCN3D,      &
                 QTCN3D,QCCN3D,QRCN3D,QNIN3D,IDEPNU,KTS,KTE)
!======================================================================
      IMPLICIT NONE 
      INTEGER, INTENT(IN) :: KTS,KTE,IDEPNU
      INTEGER :: K,A,IV,IV0,IV1,IV2,IV3,IV4,IM,IT,J
      INTEGER, DIMENSION(KTS:KTE) :: HID
      REAL, DIMENSION(KTS:KTE), INTENT(IN) :: W3D,P3D
      REAL, DIMENSION(KTS:KTE) :: TK3D,TC3D,QV3D,DZ3D,S3D,QC3D,QR3D,   &
            QI3D,QS3D,QG3D,QH3D,QDCN3D,QTCN3D,QCCN3D,QRCN3D,QNIN3D,    &
            NC3D,NR3D,NI3D,NS3D,NG3D,NH3D,VI3D,VS3D,VG3D,FI3D,FS3D,    &
            AI3D,AS3D,AG3D,AH3D,I3M3D,TAIR,PRES,QVAP,ADAGR,NS0,AMI,BMI,&
            AVI,BVI,AMS,AVS,BVS,AMG,AVG,BVG,AVH,BVH,DNIVD,INHGR,DCLDMT,&
            DCLDMC,DCLDMR,DTAIR,DPDT,DQVAP,QACac,QACar,NACac,NACar,    &
            VT_QC,VT_QR,VT_QI,VT_QS,VT_QG,VT_QH,VT_NC,VT_NR,VT_NI,     &
            VT_NS,VT_NG,VT_NH,VT_VI,VT_VS,VT_VG,VT_FI,VT_FS,VT_AI,     &
            VT_AS,VT_AG,VT_AH,VTI3M,TK0,P40,RHO,CPM,GQCTR,ESW,ESI,QVSW,&
            QVSI,SSRW,SSRI,XXLV,XXLS,XXLF,XDNC,XDNR,IPF,IPG,I3M0,IAREA,&
            LAMC,LAMR,LAMI,LAMS,LAMG,LAMH,ZETA,IMTMP,ZETA2,BEST,SMLF,  &
            GMLF,RHOAJ,GUR,FUR,FSQR,FSQI,FSQS,FSQG,FSQH,SIGC,MNRC,KDX, &
            BDR,BDI,BDS,BDG,BDH,R3M3D,S3M3D,G3M3D,H3M3D,RHOI,RHOS,RHOG,&
            IASPR,SASPR,MVDC,MVDR,MVDI,MVDS,MVDG,MVDH,AFAC,AFAR,AFAI,  &
            AFAS,AFAG,AFAH
      REAL :: SR,ICENCV,CLODNCV,RAINNC,RAINNCV,SNOWNC,SNOWNCV,GRAPNC,  &
              GRAPNCV,HAILNC,HAILNCV,DT,DTMN,SDTL,SDTS,DTL,DTS,CNDMAX, &
              TMP,XMASS,TEMPC,TEMPR,TEMPT,FINCV,FSNCV,VINCV,VSNCV,     &
              VGNCV,AINCV,ASNCV,AGNCV,AHNCV,I3MNCV
      REAL, DIMENSION(KTS:KTE,3) :: VTMEAN
      LOGICAL, DIMENSION(KTS:KTE,NAER) :: DO_FR                         ! re-construct r-cutoff from XA & XC for once only
      REAL, DIMENSION(KTS:KTE,NAER) :: RX0                              ! CUTOFF RADIUS (m)
      REAL, DIMENSION(KTS:KTE,NAERT) :: AERO                            ! AEROSOL MIXING RATIO
      REAL, DIMENSION(MAER) :: ABCD
      REAL, DIMENSION(NCCN,KTS:KTE,NAER) :: ZCCNS                       ! AEROSOL # OF EACH MODE (#/kg)
      REAL, PARAMETER :: CQR1 = 1.5943E+01, CQR2 = 1.1898E+00, CQR3 = 4.0073E+00
      REAL, PARAMETER :: AVTR = 2.1454E+00, BVTR = -2.2812E-04, CVTR = 2.9676E-09

!----- AT SUBSATURATION, REMOVE SMALL AMOUNTS OF HYDROMETEORS ----------
      DO K = KTS,KTE
         AERO(K,1) = MAX(0.,QDCN3D(K))
         AERO(K,2) = MAX(0.,QTCN3D(K))
         AERO(K,3) = MAX(0.,QCCN3D(K))
         AERO(K,4) = MAX(0.,QRCN3D(K))
         AERO(K,5) = MAX(0.,QNIN3D(K))
         TC3D(K)   = TK3D(K)-TK0C
         ESW(K)    = MIN(0.99*P3D(K),POLYSVP(TK3D(K),0))
         ESI(K)    = MIN(0.99*P3D(K),POLYSVP(TK3D(K),1))
         IF (ESI(K).GT.ESW(K)) ESI(K) = ESW(K)
         SSRW(K)   = QV3D(K)/(0.622*ESW(K)/(P3D(K)-ESW(K)))-1.
         SSRI(K)   = QV3D(K)/(0.622*ESI(K)/(P3D(K)-ESI(K)))-1.
         XXLV(K)   = 3.1484E6-2370.*TK3D(K)
         XXLS(K)   = 3.15E6-2370.*TK3D(K)+0.3337E6
         XXLF(K)   = 2836310.8-(3.1484E6-2370.*TK3D(K))
         CPM(K)    = CP*(1.+0.887*QV3D(K))
         HID(K)    = MAX(MIN(NINT(ABS(TC3D(K))/0.25),120),0)
         IF (QC3D(K).GE.QSMALL.AND.NC3D(K).LT.NSMALL) THEN
            NC3D(K) = NC0
            AFAC(K) = AFAC0
         ENDIF
         IF (QR3D(K).GE.QSMALL.AND.NR3D(K).LT.NSMALL) THEN
            NR3D(K) = MAX(1.,(QR3D(K)/PI/RHOW)**0.25*NR0**0.75)
            AFAR(K) = AFAR0
         ENDIF
         IF (QI3D(K).GE.QSMALL.AND.NI3D(K).LT.NSMALL) THEN
            IF (ICE_RHOI.EQ.0) THEN
               RHOI(K)  = RHOI0
               VI3D(K)  = 0.
            ELSEIF (ICE_RHOI.EQ.1) THEN
               QVSI(K)  = 0.622*ESI(K)/(P3D(K)-ESI(K))
               RHO(K)   = P3D(K)/(TK3D(K)*(1.+0.61*QV3D(K)))/R
               INHGR(K) = ITBLE(HID(K))
               DNIVD(K) = RHOI0*EXP(-3.*MAX((QV3D(K)-QVSI(K))-5.E-5,   &
                          0.)/INHGR(K))
               RHOI(K)  = MAX(RHOIMIN,DNIVD(K))
               VI3D(K)  = QI3D(K)/RHOI(K)
            ELSEIF (ICE_RHOI.EQ.2) THEN
               RHOI(K)  = RHOI1
               VI3D(K)  = 0.
            ENDIF
            NI3D(K)  = MAX(1.,(QI3D(K)/PI/RHOI(K))**0.25*NI0**0.75)
            I3M3D(K) = QI3D(K)*V2M3/RHOI(K)
            IF (ICE_SHAPE.EQ.0) THEN
               FI3D(K)  = I3M3D(K)
            ELSEIF (ICE_SHAPE.EQ.1) THEN
               ADAGR(K) = MAX(MIN(ITBLE(HID(K)),2.),0.4)
               ZETA(K)  = (ADAGR(K)-1.)/(ADAGR(K)+2.)
               I3M0(K)  = NI3D(K)*DI0**3.
               IF (TC3D(K).LT.-30.) ZETA(K) = MIN(ZETA(K),0.1)
               IF (TC3D(K).LT.-40.) ZETA(K) = MIN(ZETA(K),0.05)
               FI3D(K)  = (I3M3D(K)/I3M0(K))**ZETA(K)*I3M3D(K)
            ENDIF
            AFAI(K) = AFAI0
            IF (AFAI_3M.EQ.0.OR.AFAI_3M.EQ.2) THEN
               AI3D(K) = 0.
            ELSEIF (AFAI_3M.EQ.1) THEN
               AI3D(K) = (KCIMIN*NI3D(K)*I3M3D(K)**2.)**THRD
            ENDIF
         ENDIF
         IF (QS3D(K).GE.QSMALL.AND.NS3D(K).LT.NSMALL) THEN
            NS0(K) = MIN(2.E8,2.E6*EXP(-0.12*MIN(-1.E-3,TK3D(K)-TK0C)))
            IF (ICE_RHOS.EQ.0) THEN
               RHOS(K) = RHOS1
               VS3D(K) = 0.
            ELSE
               IF (ICE_RHOI.EQ.0.AND.ICE_SHAPE.EQ.0) THEN
                  RHOS(K) = RHOS0
               ELSE
                  RHOS(K) = MIN(RHOIMAX,MAX(RHOIMIN,ROSA**(4./(ROSB+   &
                            4.))*6.**(4.*ROSB/(3.*ROSB+12.))*(PI*      &
                            NS0(K)/QS3D(K))**(-1.*ROSB/(ROSB+4.))))
               ENDIF
               IF (ICE_RHOS.EQ.1) THEN
                  VS3D(K) = QS3D(K)/RHOS(K)
               ELSEIF (ICE_RHOS.EQ.2) THEN
                  VS3D(K) = 0.
               ENDIF
            ENDIF
            NS3D(K) = MAX(1.,(QS3D(K)/PI/RHOS(K))**0.25*NS0(K)**0.75)
            IF (AGG_SHAPE.EQ.0) THEN
               SASPR(K) = 1.
               FS3D(K)  = 0.
            ELSEIF (AGG_SHAPE.EQ.1) THEN
               IF (ICE_RHOI.EQ.0.AND.ICE_SHAPE.EQ.0) THEN
                  SASPR(K) = 0.8
               ELSE
                  SASPR(K) = MIN(SASMAX,0.76239415+0.00010442*RHOS(K))
               ENDIF
               FS3D(K) = -1.*LOG(SASPR(K))*QS3D(K)
            ENDIF
            AFAS(K) = AFAS0
            IF (AFAS_3M.EQ.0.OR.AFAS_3M.EQ.2) THEN
               AS3D(K) = 0.
            ELSEIF (AFAS_3M.EQ.1) THEN
               S3M3D(K) = QS3D(K)*V2M3/RHOS(K)
               AS3D(K)  = (KCSMIN*NS3D(K)*S3M3D(K)**2.)**THRD
            ENDIF
         ENDIF
         IF (QG3D(K).GE.QSMALL.AND.NG3D(K).LT.NSMALL) THEN
            RHOG(K) = RHOG1
            IF (ICE_RHOG.EQ.0) THEN
               VG3D(K) = 0.
            ELSEIF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               VG3D(K) = QG3D(K)/RHOG(K)
            ENDIF
            NG3D(K) = MAX(1.,(QG3D(K)/PI/RHOG(K))**0.25*NG0**0.75)
            AFAG(K) = AFAG0
            IF (AFAG_3M.EQ.0.OR.AFAG_3M.EQ.2) THEN
               AG3D(K) = 0.
            ELSEIF (AFAG_3M.EQ.1) THEN
               G3M3D(K) = QG3D(K)*V2M3/RHOG(K)
               AG3D(K)  = (KCGMIN*NG3D(K)*G3M3D(K)**2.)**THRD
            ENDIF
         ENDIF
         IF (QH3D(K).GE.QSMALL.AND.NH3D(K).LT.NSMALL) THEN
            NH3D(K) = MAX(1.,(QH3D(K)/PI/RHOH)**0.25*NH0**0.75)
            IF (AFAH_3M.EQ.0.OR.AFAH_3M.EQ.2) THEN
               AH3D(K) = 0.
            ELSEIF (AFAH_3M.EQ.1) THEN
               H3M3D(K) = QH3D(K)*V2M3/RHOH
               AH3D(K)  = (KCHMIN*NH3D(K)*H3M3D(K)**2.)**THRD
            ENDIF
         ENDIF
         IF (SSRW(K).LT.-1.E-1) THEN
            IF (QR3D(K).LT.QSMAL1) THEN
               QV3D(K) = QV3D(K)+QR3D(K)
               TK3D(K) = TK3D(K)-QR3D(K)*XXLV(K)/CPM(K)
               QR3D(K) = 0.; NR3D(K) = 0.; MVDR(K) = 0.; AFAR(K) = 0.
            ENDIF
            IF (QC3D(K).LT.QSMAL1) THEN
               QV3D(K) = QV3D(K)+QC3D(K)
               TK3D(K) = TK3D(K)-QC3D(K)*XXLV(K)/CPM(K)
               QC3D(K) = 0.; NC3D(K) = 0.; MVDC(K) = 0.; AFAC(K) = 0.
            ENDIF
         ENDIF
         IF (SSRI(K).LT.-1.E-1) THEN
            IF (QI3D(K).LT.QSMAL1) THEN
               QV3D(K) = QV3D(K)+QI3D(K)
               TK3D(K) = TK3D(K)-QI3D(K)*XXLS(K)/CPM(K)
               QI3D(K) = 0.; NI3D(K) = 0.; MVDI(K) = 0.
               I3M3D(K) = 0.; FI3D(K) = 0.; VI3D(K) = 0.; AFAI(K) = 0.
               IASPR(K) = 1.; ADAGR(K) = 1.; RHOI(K) = 0.; AI3D(K) = 0.
            ENDIF
            IF (QS3D(K).LT.QSMAL1) THEN
               QV3D(K) = QV3D(K)+QS3D(K)
               TK3D(K) = TK3D(K)-QS3D(K)*XXLS(K)/CPM(K)
               QS3D(K) = 0.; NS3D(K) = 0.; VS3D(K) = 0.; AS3D(K) = 0.
               MVDS(K) = 0.; RHOS(K) = 0.; AFAS(K) = 0.
               FS3D(K) = 0.; SASPR(K) = 1.
            ENDIF
            IF (QG3D(K).LT.QSMAL1) THEN
               QV3D(K) = QV3D(K)+QG3D(K)
               TK3D(K) = TK3D(K)-QG3D(K)*XXLS(K)/CPM(K)
               QG3D(K) = 0.; NG3D(K) = 0.; MVDG(K) = 0.
               VG3D(K) = 0.; AG3D(K) = 0.; AFAG(K) = 0.; RHOG(K) = 0.
            ENDIF
            IF (QH3D(K).LT.QSMAL1) THEN
               QV3D(K) = QV3D(K)+QH3D(K)
               TK3D(K) = TK3D(K)-QH3D(K)*XXLS(K)/CPM(K)
               QH3D(K) = 0.; NH3D(K) = 0.; MVDH(K) = 0.
               AH3D(K) = 0.; AFAH(K) = 0.
            ENDIF
         ENDIF
!----- MELT VERY SMALL SNOW/GRAUPEL/HAIL MIXING RATIOS, ADD TO RAIN ----
         IF (TK3D(K).GE.TK0C.AND.QS3D(K).LT.QLIMIT) THEN
            QR3D(K) = QR3D(K)+QS3D(K)
            NR3D(K) = NR3D(K)+NS3D(K)
            TK3D(K) = TK3D(K)-QS3D(K)*XXLF(K)/CPM(K)
            QS3D(K) = 0.; NS3D(K) = 0.; VS3D(K) = 0.; AS3D(K) = 0.
            MVDS(K) = 0.; RHOS(K) = 0.; AFAS(K) = 0.
            FS3D(K) = 0.; SASPR(K) = 1.
         ENDIF
         IF (TK3D(K).GE.TK0C.AND.QG3D(K).LT.QLIMIT) THEN
            QR3D(K) = QR3D(K)+QG3D(K)
            NR3D(K) = NR3D(K)+NG3D(K)
            TK3D(K) = TK3D(K)-QG3D(K)*XXLF(K)/CPM(K)
            QG3D(K) = 0.; NG3D(K) = 0.; MVDG(K) = 0.
            VG3D(K) = 0.; AG3D(K) = 0.; AFAG(K) = 0.; RHOG(K) = 0.
         ENDIF
         IF (TK3D(K).GE.TK0C.AND.QH3D(K).LT.QLIMIT) THEN
            QR3D(K) = QR3D(K)+QH3D(K)
            NR3D(K) = NR3D(K)+NH3D(K)
            TK3D(K) = TK3D(K)-QH3D(K)*XXLF(K)/CPM(K)
            QH3D(K) = 0.; NH3D(K) = 0.; MVDH(K) = 0.
            AH3D(K) = 0.; AFAH(K) = 0.
         ENDIF
      ENDDO
!------- INITIALIZE VARIABLE,ONLY T&QV ARE TREATED DIFFERENTLY ---------
      DO K = KTS,KTE
         TAIR(K)  = TK3D(K)
         QVAP(K)  = QV3D(K)
         DTAIR(K) = -W3D(K)*G/CP                                        ! Lagrangian parcel dry adiabatic lapse rate
         DQVAP(K) = 0.
         PRES(K)  = P3D(K)
         TK0(K)   = TAIR(K)-DTAIR(K)*DT                                 ! diabatically decend parcel
         P40(K)   = PRES(K)*(TK0(K)/TAIR(K))**(CP/R)
         TAIR(K)  = TK0(K)
         IV1 = 1
         DO IV = 1,NAER
            IV2 = IV1+1
            DO_FR(K,IV) = .TRUE.                                        ! FIND CUT OFF RADIUS
            RX0(K,IV) = 99.                                             ! initial R cutoff
            IF (NAERN(IV).GE.2) THEN                                    ! have info of dry and total aerosol 
               AERO(K,IV1) = MIN(AERO(K,IV1),AERO(K,IV2))               ! dry aerosol < total aerosol
               DO IM = 1,NCCN                                           ! calculate CCN # [#/kg] from CCN mass
                  ZCCNS(IM,K,IV) = AERO(K,IV2)*WMAS(IM,IV)*RFACT(IM,IV)
               ENDDO
            ELSE                                                        ! only info available is total aerosol 
               DO IM = 1,NCCN                                           ! calculate CCN # [#/kg] from CCN mass
                  ZCCNS(IM,K,IV) = AERO(K,IV1)*WMAS(IM,IV)*RFACT(IM,IV)
               ENDDO
            ENDIF
            IV1 = IV1+NAERN(IV)
         ENDDO
      ENDDO                                                             ! INITILIZATION LOOPS

      SDTL = 0.                                                         ! initial integrated largedt timestep
      DO IT = 1,INT((DT-DT20S*0.01)/DT20S)+1
         DTL = MIN(DT20S,DT-SDTL)                                       ! initial largedt timestep
         DO K = KTS,KTE
            DTS   = DTL                                                 ! initial smalldt timestep (could be changed later)
            SDTS  = 0.                                                  ! initial integrated smalldt timestep
            TMP   = W3D(K)*G/CP/TK0(K)                                  ! 1-TAIR/TK0
            XDNC(K) = 0.; XDNR(K) = 0.; GQCTR(K) = 0.; ADAGR(K) = 1.
            RHOI(K) = RHOI0; RHOS(K) = RHOS0; RHOG(K) = RHOG0
            AFAC(K) = AFAC0; AFAR(K) = AFAR0; AFAI(K) = AFAI0
            AFAS(K) = AFAS0; AFAG(K) = AFAG0; AFAH(K) = AFAH0
            AVI(K) = AVI0; BVI(K) = BVI0; AVS(K) = AVS0; BVS(K) = BVS0
            AVG(K) = AVG0; BVG(K) = BVG0; AVH(K) = AVH0; BVH(K) = BVH0
            AMI(K) = AMI0; BMI(K) = BMI0; AMS(K) = AMS0; AMG(K) = AMG0
            IASPR(K) = 1.; SASPR(K) = 1.; QACac(K) = 0.; QACar(K) = 0.
 111        CONTINUE 

            PRES(K) = P40(K)*(1.-TMP*(SDTL+SDTS))**(CP/R) 
            DPDT(K) = -G*W3D(K)*PRES(K)/(R*TK0(K)*(1.-TMP*(SDTL+SDTS)))
            ESW(K)  = MIN(0.99*PRES(K),POLYSVP(TAIR(K),0))
            ESI(K)  = MIN(0.99*PRES(K),POLYSVP(TAIR(K),1))
            IF (ESI(K).GT.ESW(K)) ESI(K) = ESW(K)
            QVSW(K) = 0.622*ESW(K)/(PRES(K)-ESW(K))
            QVSI(K) = 0.622*ESI(K)/(PRES(K)-ESI(K))
            SSRW(K) = QVAP(K)/QVSW(K)-1.
            SSRI(K) = QVAP(K)/QVSI(K)-1.
            XXLV(K) = 3.1484E6-2370.*TAIR(K)
            CPM(K)  = CP*(1.+0.887*QVAP(K))
            RHO(K)  = PRES(K)/(TAIR(K)*(1.+0.61*QVAP(K)))/R             ! [KG/M**3]
            CNDMAX  = (QVAP(K)-QVSW(K))/(1.+XXLV(K)**2.*QVAP(K)/       &
                      (CPM(K)*RV*TAIR(K)**2.))
!--------------------------- CCN ACTIVATION. ---------------------------
            IF (SSRW(K).GT.RSMALL.AND.SSRW(K).GT.S3D(K)) THEN 
               DCLDMT(K) = 0.; DCLDMC(K) = 0.; DCLDMR(K) = 0.
               IV1 = 1
               DO IV = 1,NAER
                  DO J = 1,NAERN(IV)
                     ABCD(J) = AERO(K,IV1-1+J)
                  ENDDO
                  IF (AERO(K,IV1).GT.1.E-20.AND.NAERN(IV).GE.2.AND.    &
                     W3D(K).GT.1.E-5) THEN
                     CALL ACTIVA(TAIR(K),W3D(K),NC3D(K),NR3D(K),ABCD,  &
                          QACac(K),QACar(K),RX0(K,IV),ZCCNS(1,K,IV),   &
                          SSRW(K),IV,DO_FR(K,IV))
                     DCLDMC(K) = DCLDMC(K)+QACac(K)
                     DCLDMR(K) = DCLDMR(K)+QACar(K)
                     DO J = 1,NAERN(IV)
                        AERO(K,IV1-1+J) = ABCD(J)
                     ENDDO
                  ENDIF
                  IV1 = IV1+NAERN(IV)
               ENDDO
               DCLDMT(K) = DCLDMC(K)+DCLDMR(K)
               IF (DCLDMT(K).GT.RLIMIT) THEN
                  IF (DCLDMT(K).GT.CNDMAX) THEN
                     DCLDMC(K) = DCLDMC(K)*(CNDMAX/DCLDMT(K))
                     DCLDMR(K) = DCLDMR(K)*(CNDMAX/DCLDMT(K))
                     DCLDMT(K) = DCLDMC(K)+DCLDMR(K)
                  ENDIF
                  QVAP(K) = MAX(0.,QVAP(K)-DCLDMT(K))
                  QC3D(K) = MAX(0.,QC3D(K)+DCLDMC(K))
                  QR3D(K) = MAX(0.,QR3D(K)+DCLDMR(K))
                  CPM(K)  = CP*(1.+0.887*QVAP(K))
                  XXLV(K) = 3.1484E6-2370.*TAIR(K)
                  TAIR(K) = TAIR(K)+XXLV(K)*DCLDMT(K)/CPM(K)
               ENDIF
               ESW(K)  = MIN(0.99*PRES(K),POLYSVP(TAIR(K),0))
               ESI(K)  = MIN(0.99*PRES(K),POLYSVP(TAIR(K),1))
               IF (ESI(K).GT.ESW(K)) ESI(K) = ESW(K)
               SSRW(K) = QVAP(K)/(0.622*ESW(K)/(PRES(K)-ESW(K)))-1.
               SSRI(K) = QVAP(K)/(0.622*ESI(K)/(PRES(K)-ESI(K)))-1.
               RHO(K)  = PRES(K)/(TAIR(K)*(1.+0.61*QVAP(K)))/R
            ENDIF 
            S3D(K) = MAX(S3D(K),SSRW(K))
!------------------------- ICE NUCLEATION ------------------------------
            IF (TAIR(K).LT.(TK0C-2.).AND.SSRI(K).GT.1.E-4) THEN
               CALL ICENU(TAIR(K),PRES(K),DTS,RHO(K),QVAP(K),QI3D(K),  &
                    NI3D(K),VI3D(K),FI3D(K),AI3D(K),I3M3D(K),AERO(K,5),&
                    IDEPNU)
               ESW(K)  = MIN(0.99*PRES(K),POLYSVP(TAIR(K),0))
               ESI(K)  = MIN(0.99*PRES(K),POLYSVP(TAIR(K),1))
               IF (ESI(K).GT.ESW(K)) ESI(K) = ESW(K)
               SSRW(K) = QVAP(K)/(0.622*ESW(K)/(PRES(K)-ESW(K)))-1.
               SSRI(K) = QVAP(K)/(0.622*ESI(K)/(PRES(K)-ESI(K)))-1.
               RHO(K)  = PRES(K)/(TAIR(K)*(1.+0.61*QVAP(K)))/R
            ENDIF
!-----------MICROPHYSICAL TENDENCY CALCULATION IN SMALL TIMESTEP -------
            IF (QC3D(K).GE.QSMALL.OR.QR3D(K).GE.QSMALL.OR.QI3D(K).GE.  &
               QSMALL.OR.QS3D(K).GE.QSMALL.OR.QG3D(K).GE.QSMALL.OR.    &
               QH3D(K).GE.QSMALL) THEN
               CALL SMALL_DT(DTL,DTS,SDTS,DTAIR(K),DQVAP(K),DPDT(K),   &
                    RHO(K),TAIR(K),PRES(K),QVAP(K),QC3D(K),QR3D(K),    &
                    QI3D(K),QS3D(K),QG3D(K),QH3D(K),NC3D(K),NR3D(K),   &
                    NI3D(K),NS3D(K),NG3D(K),NH3D(K),VI3D(K),VS3D(K),   &
                    VG3D(K),FI3D(K),FS3D(K),AI3D(K),AS3D(K),AG3D(K),   &
                    AH3D(K),I3M3D(K),SASPR(K),XDNC(K),XDNR(K),GQCTR(K),&
                    SMLF(K),GMLF(K))
            ENDIF
!----------------- UPDATE TIME & PROGNOSTIC VARIABLES ------------------
            TAIR(K) = TAIR(K)+DTAIR(K)*DTS
            QVAP(K) = QVAP(K)+DQVAP(K)*DTS
            SDTS    = SDTS+DTS

            IF (SDTS.GE.DTL) GOTO 333
               GOTO 111
 333        CONTINUE
!----------------- RAIN DROP DEACTIVATION. -----------------------------
            IF (XDNR(K).GT.1.) THEN
               IV0 = 0
               DO IV = 1,NAER
                  IF (NAERN(IV).EQ.4) THEN
                     IV3 = IV0+3
                     IV4 = IV0+4
                     IF (NR3D(K).GT.1.) THEN
                        XMASS = AERO(K,IV4)*MIN(1.,XDNR(K)/NR3D(K))
                     ELSE
                        XMASS = AERO(K,IV4)
                     ENDIF 
                     AERO(K,IV3) = AERO(K,IV3)+XMASS
                     AERO(K,IV4) = AERO(K,IV4)-XMASS
                  ENDIF
                  IV0 = IV0+NAERN(IV)
               ENDDO
            ENDIF
!----------------- CLOUD DROP DEACTIVATION. ----------------------------
            IF (XDNC(K).GT.1.) THEN 
               IV1 = 1
               DO IV = 1,NAER
                  IV2 = IV1+1
                  IV3 = IV2+1
                  IF (NAERN(IV).GE.3) THEN
                     IF (NC3D(K).GT.1.) THEN
                        XMASS = AERO(K,IV3)*MIN(1.,XDNC(K)/NC3D(K))
                     ELSE
                        XMASS = AERO(K,IV3)
                     ENDIF 
                     AERO(K,IV1) = AERO(K,IV1)+XMASS
                     AERO(K,IV3) = AERO(K,IV3)-XMASS
                     AERO(K,IV1) = MIN(AERO(K,IV1),AERO(K,IV2))
                  ELSEIF (NAERN(IV).EQ.2) THEN
                     CALL DEACTIVA(XDNC(K),QC3D(K),NC3D(K),AERO(K,IV1),&
                          AERO(K,IV2),AERO(K,IV3),RX0(K,IV),           &
                          ZCCNS(1,K,IV),IV,DO_FR(K,IV))
                  ELSE
                     XMASS = AERO(K,IV1)*MIN(1.,XDNC(K)/               &
                            (ZCCNS(2,K,IV)+ZCCNS(3,K,IV)))              ! COARSE AND ACCUMULATION MODE
                     AERO(K,IV1) = AERO(K,IV1)+XMASS
                  ENDIF
                  IV1 = IV1+NAERN(IV)
               ENDDO 
            ENDIF 
            PRES(K) = P40(K)*(1.-TMP*(SDTL+SDTS))**(CP/R) 
            RHO(K)  = PRES(K)/(TAIR(K)*(1.+0.61*QVAP(K)))/R
         ENDDO
!-----------MICROPHYSICAL TENDENCY CALCULATION IN LARGE TIMESTEP -------
         DO K = KTS,KTE
            IF (QC3D(K).GE.QSMALL.OR.QR3D(K).GE.QSMALL.OR.QI3D(K).GE.  &
               QSMALL.OR.QS3D(K).GE.QSMALL.OR.QG3D(K).GE.QSMALL.OR.    &
               QH3D(K).GE.QSMALL) THEN
               CALL LARGE_DT(DTL,TAIR(K),QVAP(K),PRES(K),W3D(K),RHO(K),&
                    QC3D(K),QR3D(K),QI3D(K),QS3D(K),QG3D(K),QH3D(K),   &
                    NC3D(K),NR3D(K),NI3D(K),NS3D(K),NG3D(K),NH3D(K),   &
                    VI3D(K),VS3D(K),VG3D(K),FI3D(K),FS3D(K),AI3D(K),   &
                    AS3D(K),AG3D(K),AH3D(K),I3M3D(K),SASPR(K),GQCTR(K),&
                    SMLF(K),GMLF(K))
            ENDIF
            IV0 = 0
            DO IV = 1,NAER                                              ! for all aerosol compoment
               IF (NAERN(IV).GE.4) THEN                                 ! with aerosol in cloud/precip. recored
                  IV3 = IV0+3                                           ! aerosols in clouds ll
                  IV4 = IV0+4                                           ! aerosols in rains
                  IF (GQCTR(K).GT.0.) THEN                              ! + from cloud to precip.
                      XMASS = GQCTR(K)*AERO(K,IV3)                      ! + from cloud to precip.
                  ELSE                                                  ! - from precip. to cloud
                      XMASS = GQCTR(K)*AERO(K,IV4)                      ! - from precip. to cloud
                  ENDIF
                  AERO(K,IV3) = AERO(K,IV3)-XMASS                       ! update aerosol in cloud
                  AERO(K,IV4) = AERO(K,IV4)+XMASS                       ! update aerosol in precipi.
               ENDIF
               IV0 = IV0+NAERN(IV)
            ENDDO
         ENDDO 
         SDTL = SDTL+DTL
      ENDDO                                                             ! FOR IT LOOPS
!----- AT SUBSATURATION, REMOVE SMALL AMOUNTS OF HYDROMETEORS ----------
      DO K = KTS,KTE
         ESW(K)  = MIN(0.99*PRES(K),POLYSVP(TAIR(K),0))
         ESI(K)  = MIN(0.99*PRES(K),POLYSVP(TAIR(K),1))
         IF (ESI(K).GT.ESW(K)) ESI(K) = ESW(K)
         SSRW(K) = QVAP(K)/(0.622*ESW(K)/(PRES(K)-ESW(K)))-1.
         SSRI(K) = QVAP(K)/(0.622*ESI(K)/(PRES(K)-ESI(K)))-1.
         XXLV(K) = 3.1484E6-2370.*TAIR(K)
         XXLS(K) = 3.15E6-2.37E3*TAIR(K)+0.3337E6
         CPM(K)  = CP*(1.+0.887*QVAP(K))
         IF (SSRW(K).LT.-1.E-1) THEN
            IF (QR3D(K).LT.QSMAL1) THEN
               QVAP(K) = QVAP(K)+QR3D(K)
               TAIR(K) = TAIR(K)-QR3D(K)*XXLV(K)/CPM(K)
               QR3D(K) = 0.; NR3D(K) = 0.; MVDR(K) = 0.; AFAR(K) = 0.
            ENDIF
            IF (QC3D(K).LT.QSMAL1) THEN
               QVAP(K) = QVAP(K)+QC3D(K)
               TAIR(K) = TAIR(K)-QC3D(K)*XXLV(K)/CPM(K)
               QC3D(K) = 0.; NC3D(K) = 0.; MVDC(K) = 0.; AFAC(K) = 0.
            ENDIF
         ENDIF
         IF (SSRI(K).LT.-1.E-1) THEN
            IF (QI3D(K).LT.QSMAL1) THEN
               QVAP(K) = QVAP(K)+QI3D(K)
               TAIR(K) = TAIR(K)-QI3D(K)*XXLS(K)/CPM(K)
               QI3D(K) = 0.; NI3D(K) = 0.; AI3D(K) = 0.
               MVDI(K) = 0.; I3M3D(K) = 0.; FI3D(K) = 0.; VI3D(K) = 0.
               AFAI(K) = 0.; IASPR(K) = 1.; ADAGR(K) = 1.; RHOI(K) = 0.
            ENDIF
            IF (QS3D(K).LT.QSMAL1) THEN
               QVAP(K) = QVAP(K)+QS3D(K)
               TAIR(K) = TAIR(K)-QS3D(K)*XXLS(K)/CPM(K)
               QS3D(K) = 0.; NS3D(K) = 0.; VS3D(K) = 0.; AS3D(K) = 0.
               MVDS(K) = 0.; RHOS(K) = 0.; AFAS(K) = 0.
               FS3D(K) = 0.; SASPR(K) = 1.
            ENDIF
            IF (QG3D(K).LT.QSMAL1) THEN
               QVAP(K) = QVAP(K)+QG3D(K)
               TAIR(K) = TAIR(K)-QG3D(K)*XXLS(K)/CPM(K)
               QG3D(K) = 0.; NG3D(K) = 0.; MVDG(K) = 0.
               VG3D(K) = 0.; AG3D(K) = 0.; AFAG(K) = 0.; RHOG(K) = 0.
            ENDIF
            IF (QH3D(K).LT.QSMAL1) THEN
               QVAP(K) = QVAP(K)+QH3D(K)
               TAIR(K) = TAIR(K)-QH3D(K)*XXLS(K)/CPM(K)
               QH3D(K) = 0.; NH3D(K) = 0.; MVDH(K) = 0.
               AH3D(K) = 0.; AFAH(K) = 0.
            ENDIF
         ENDIF
         IF (TAIR(K).GT.TK0C.AND.QI3D(K).GT.0.) THEN
            QR3D(K) = QR3D(K)+QI3D(K)
            TAIR(K) = TAIR(K)-QI3D(K)*XXLF(K)/CPM(K)
            NR3D(K) = NR3D(K)+NI3D(K)
            QI3D(K) = 0.; NI3D(K) = 0.; MVDI(K) = 0.
            I3M3D(K) = 0.; FI3D(K) = 0.; VI3D(K) = 0.; AFAI(K) = 0.
            IASPR(K) = 1.; ADAGR(K) = 1.; RHOI(K) = 0.; AI3D(K) = 0.
         ENDIF
         IF (QC3D(K).LT.QSMALL.OR.NC3D(K).LT.NSMALL) THEN
            QVAP(K) = QVAP(K)+QC3D(K)
            TAIR(K) = TAIR(K)-QC3D(K)*XXLV(K)/CPM(K)
            QC3D(K) = 0.; NC3D(K) = 0.; MVDC(K) = 0.; AFAC(K) = 0.
         ENDIF
         IF (QR3D(K).LT.QSMALL.OR.NR3D(K).LT.NSMALL) THEN
            QVAP(K) = QVAP(K)+QR3D(K)
            TAIR(K) = TAIR(K)-QR3D(K)*XXLV(K)/CPM(K)
            QR3D(K) = 0.; NR3D(K) = 0.; MVDR(K) = 0.; AFAR(K) = 0.
         ENDIF
         IF (QI3D(K).LT.QSMALL.OR.NI3D(K).LT.NSMALL) THEN
            QVAP(K) = QVAP(K)+QI3D(K)
            TAIR(K) = TAIR(K)-QI3D(K)*XXLS(K)/CPM(K)
            QI3D(K) = 0.; NI3D(K) = 0.; MVDI(K) = 0.
            I3M3D(K) = 0.; FI3D(K) = 0.; VI3D(K) = 0.; AFAI(K) = 0.
            IASPR(K) = 1.; ADAGR(K) = 1.; RHOI(K) = 0.; AI3D(K) = 0.
         ENDIF
         IF (AI3D(K).LT.ASMALL) THEN
            AI3D(K) = 0.
         ENDIF
         IF (AS3D(K).LT.ASMALL) THEN
            AS3D(K) = 0.
         ENDIF
         IF (AG3D(K).LT.ASMALL) THEN
            AG3D(K) = 0.
         ENDIF
         IF (AH3D(K).LT.ASMALL) THEN
            AH3D(K) = 0.
         ENDIF
         IF (I3M3D(K).LT.ISMALL.OR.FI3D(K).LT.ISMALL) THEN
            FI3D(K) = 0.; I3M3D(K) = 0.; ADAGR(K) = 1.; ZETA(K) = 0.
         ENDIF
         IF (QS3D(K).LT.QSMALL.OR.NS3D(K).LT.NSMALL) THEN
            QVAP(K) = QVAP(K)+QS3D(K)
            TAIR(K) = TAIR(K)-QS3D(K)*XXLS(K)/CPM(K)
            QS3D(K) = 0.; NS3D(K) = 0.; VS3D(K) = 0.; AS3D(K) = 0.
            MVDS(K) = 0.; RHOS(K) = 0.; AFAS(K) = 0.
            FS3D(K) = 0.; SASPR(K) = 1.
         ENDIF
         IF (QG3D(K).LT.QSMALL.OR.NG3D(K).LT.NSMALL) THEN
            QVAP(K) = QVAP(K)+QG3D(K)
            TAIR(K) = TAIR(K)-QG3D(K)*XXLS(K)/CPM(K)
            QG3D(K) = 0.; NG3D(K) = 0.; MVDG(K) = 0.
            VG3D(K) = 0.; AG3D(K) = 0.; AFAG(K) = 0.; RHOG(K) = 0.
         ENDIF
         IF (QH3D(K).LT.QSMALL.OR.NH3D(K).LT.NSMALL) THEN
            QVAP(K) = QVAP(K)+QH3D(K)
            TAIR(K) = TAIR(K)-QH3D(K)*XXLS(K)/CPM(K)
            QH3D(K) = 0.; NH3D(K) = 0.; MVDH(K) = 0.
            AH3D(K) = 0.; AFAH(K) = 0.
         ENDIF
         RHO(K) = PRES(K)/(TAIR(K)*(1.+0.61*QVAP(K)))/R
         VT_QC(K) = 0.; VT_NC(K) = 0.; VT_QR(K) = 0.; VT_NR(K) = 0.
         VT_QI(K) = 0.; VT_NI(K) = 0.; VT_QS(K) = 0.; VT_NS(K) = 0.
         VT_QG(K) = 0.; VT_NG(K) = 0.; VT_QH(K) = 0.; VT_NH(K) = 0.
         VT_VI(K) = 0.; VT_VS(K) = 0.; VT_VG(K) = 0.; VT_FI(K) = 0.
         VT_FS(K) = 0.; VT_AI(K) = 0.; VT_AS(K) = 0.; VT_AG(K) = 0.
         VT_AH(K) = 0.; VTI3M(K) = 0.
         IF (QC3D(K).GE.QSMALL.OR.QR3D(K).GE.QSMALL.OR.QI3D(K).GE.     &
            QSMALL.OR.QS3D(K).GE.QSMALL.OR.QG3D(K).GE.QSMALL.OR.       &
            QH3D(K).GE.QSMALL) THEN
            CALL SEDI_FALL(TAIR(K),QC3D(K),QR3D(K),QI3D(K),QS3D(K),    &
                 QG3D(K),QH3D(K),NC3D(K),NR3D(K),NI3D(K),NS3D(K),      &
                 NG3D(K),NH3D(K),VI3D(K),VS3D(K),VG3D(K),FI3D(K),      &
                 FS3D(K),AI3D(K),AS3D(K),AG3D(K),AH3D(K),I3M3D(K),     &
                 VT_QC(K),VT_QR(K),VT_QI(K),VT_QS(K),VT_QG(K),VT_QH(K),&
                 VT_NC(K),VT_NR(K),VT_NI(K),VT_NS(K),VT_NG(K),VT_NH(K),&
                 VT_VI(K),VT_VS(K),VT_VG(K),VT_FI(K),VT_FS(K),VT_AI(K),&
                 VT_AS(K),VT_AG(K),VT_AH(K),VTI3M(K),SASPR(K),RHO(K),  &
                 SMLF(K),GMLF(K))
         ENDIF
         IF (QC3D(K).GE.QSMALL.AND.NC3D(K).GE.NSMALL) THEN
            MVDC(K) = (QC3D(K)*iAMW/NC3D(K))**THRD
            IF (MVDC(K).GT.DCR) THEN
               QR3D(K) = QR3D(K)+QC3D(K)
               NR3D(K) = NR3D(K)+NC3D(K)
               QC3D(K) = 0.; NC3D(K) = 0.; AFAC(K) = 0.
            ENDIF
         ENDIF
         IF (QR3D(K).GE.QSMALL.AND.NR3D(K).GE.NSMALL) THEN
            MVDR(K) = (QR3D(K)*iAMW/NR3D(K))**THRD
            IF (MVDR(K).LT.DCR) THEN
               QC3D(K) = QC3D(K)+QR3D(K)
               NC3D(K) = NC3D(K)+NR3D(K)
               QR3D(K) = 0.; NR3D(K) = 0.; AFAR(K) = 0.
            ENDIF
         ENDIF
         IF (QH3D(K).GE.QSMALL.AND.NH3D(K).GE.NSMALL) THEN
            MVDH(K) = (QH3D(K)*iAMH/NH3D(K))**THRD
            IF (MVDH(K).LT.DHMIN) THEN
               QG3D(K) = QG3D(K)+QH3D(K)
               NG3D(K) = NG3D(K)+NH3D(K)
               VG3D(K) = VG3D(K)+QH3D(K)/RHOG0
               QH3D(K) = 0.; NH3D(K) = 0.
               IF (AH3D(K).GE.ASMALL.AND.AFAG_3M.EQ.1) THEN
                  AG3D(K) = AG3D(K)+AH3D(K)
                  AH3D(K) = 0.; AFAH(K) = 0.
               ENDIF
            ENDIF
         ENDIF
!----------------- FOR TENDENCIES DUE TO PRECIPITATION/SEDIMENTATION ---
         TEMPC = QC3D(K)+QI3D(K)                                        ! QC+QI
         IF (TEMPC.GT.0.) THEN
            VTMEAN(K,1) = (VT_QC(K)*QC3D(K)+VT_QI(K)*QI3D(K))/TEMPC     ! SEDIMENTATION SPEED OF WET AEROSOL IN QC+QI
         ELSE
            VTMEAN(K,1) = 0.
         ENDIF
         TEMPR = QR3D(K)+QS3D(K)+QG3D(K)+QH3D(K)                        ! QR+QS+QG+QH
         IF (TEMPR.GT.0.) THEN
            VTMEAN(K,2) = (VT_QR(K)*QR3D(K)+VT_QS(K)*QS3D(K)+VT_QG(K)* &
                          QG3D(K)+VT_QH(K)*QH3D(K))/TEMPR               ! SEDIMENTATION SPEED OF WET AEROSOL IN QR+QS+QG+QH
         ELSE
            VTMEAN(K,2) = 0.
         ENDIF
         TEMPT = TEMPC+TEMPR                                            ! QC+QR+QI+QS+QG
         IF (TEMPT.GT.0.) THEN
            VTMEAN(K,3) = (VTMEAN(K,1)*TEMPC+VTMEAN(K,2)*TEMPR)/TEMPT   ! MEAN SEDIMENTATION SPEED OF WET AEROSOL IN QC+QR+QI+QS+QG+QH
          ELSE
            VTMEAN(K,3) = 0.
         ENDIF
      ENDDO
!----------------- CALCULATED FALLING TERMS ----------------------------
      CALL PTFLUX(QC3D(1),VT_QC(1),RHO(1),DZ3D(1),KTE,DT,DTMN,CLODNCV)
      CALL PTFLUX(QR3D(1),VT_QR(1),RHO(1),DZ3D(1),KTE,DT,DTMN,RAINNCV)
      CALL PTFLUX(QI3D(1),VT_QI(1),RHO(1),DZ3D(1),KTE,DT,DTMN,ICENCV)
      CALL PTFLUX(QS3D(1),VT_QS(1),RHO(1),DZ3D(1),KTE,DT,DTMN,SNOWNCV)
      CALL PTFLUX(QG3D(1),VT_QG(1),RHO(1),DZ3D(1),KTE,DT,DTMN,GRAPNCV)
      CALL PTFLUX(QH3D(1),VT_QH(1),RHO(1),DZ3D(1),KTE,DT,DTMN,HAILNCV)
      RAINNC = RAINNC+CLODNCV+RAINNCV+ICENCV+SNOWNCV+GRAPNCV+HAILNCV
      SNOWNC = SNOWNC+SNOWNCV
      GRAPNC = GRAPNC+GRAPNCV
      HAILNC = HAILNC+HAILNCV
      SR     = (ICENCV+SNOWNCV+GRAPNCV+HAILNCV)/                       &
               (CLODNCV+RAINNCV+ICENCV+SNOWNCV+GRAPNCV+HAILNCV+1.E-12)
      CALL FLFLUX(NC3D(1),VT_NC(1),RHO(1),DZ3D(1),KTE,DT)
      CALL FLFLUX(NR3D(1),VT_NR(1),RHO(1),DZ3D(1),KTE,DT)
      CALL FLFLUX(NI3D(1),VT_NI(1),RHO(1),DZ3D(1),KTE,DT)
      CALL FLFLUX(NS3D(1),VT_NS(1),RHO(1),DZ3D(1),KTE,DT)
      CALL FLFLUX(NG3D(1),VT_NG(1),RHO(1),DZ3D(1),KTE,DT)
      CALL FLFLUX(NH3D(1),VT_NH(1),RHO(1),DZ3D(1),KTE,DT)
      CALL PTFLUX(VI3D(1),VT_VI(1),RHO(1),DZ3D(1),KTE,DT,DTMN,VINCV)
      CALL PTFLUX(VS3D(1),VT_VS(1),RHO(1),DZ3D(1),KTE,DT,DTMN,VSNCV)
      CALL PTFLUX(VG3D(1),VT_VG(1),RHO(1),DZ3D(1),KTE,DT,DTMN,VGNCV)
      CALL PTFLUX(FI3D(1),VT_FI(1),RHO(1),DZ3D(1),KTE,DT,DTMN,FINCV)
      CALL PTFLUX(FS3D(1),VT_FS(1),RHO(1),DZ3D(1),KTE,DT,DTMN,FSNCV)
      CALL PTFLUX(AI3D(1),VT_AI(1),RHO(1),DZ3D(1),KTE,DT,DTMN,AINCV)
      CALL PTFLUX(AS3D(1),VT_AS(1),RHO(1),DZ3D(1),KTE,DT,DTMN,ASNCV)
      CALL PTFLUX(AG3D(1),VT_AG(1),RHO(1),DZ3D(1),KTE,DT,DTMN,AGNCV)
      CALL PTFLUX(AH3D(1),VT_AH(1),RHO(1),DZ3D(1),KTE,DT,DTMN,AHNCV)
      CALL PTFLUX(I3M3D(1),VTI3M(1),RHO(1),DZ3D(1),KTE,DT,DTMN,I3MNCV)

      IV0 = 0
      DO IV = 1,NAER
         IF (NAERN(IV).GT.3) THEN
            CALL FLFLUX(AERO(1,IV0+3),VTMEAN(1,1),RHO(1),DZ3D(1),KTE,DT)! SEDI. OF AEROSOLS IN CLOUDS
            CALL FLFLUX(AERO(1,IV0+4),VTMEAN(1,2),RHO(1),DZ3D(1),KTE,DT)! SEDI. OF AEROSOLS IN PRECIP.
         ELSEIF (NAERN(IV).EQ.3) THEN
            CALL FLFLUX(AERO(1,IV0+3),VTMEAN(1,3),RHO(1),DZ3D(1),KTE,DT)! SEDI. OF AEROSOLS IN HYDROM.
         ENDIF
         IV0 = IV0+NAERN(IV)
      ENDDO
!----------------------------------------------------------------------
      DO K = KTS,KTE
         TK3D(K) = TAIR(K)
         QV3D(K) = QVAP(K)
         RHOAJ(K) = (RHOSU/RHO(K))**0.54
         IF (QC3D(K).GE.QSMALL.AND.NC3D(K).GE.NSMALL) THEN
            CALL SOLVE_AFAC(QC3D(K),NC3D(K),LAMC(K),MVDC(K),AFAC(K))
         ELSE
            QC3D(K) = 0.; NC3D(K) = 0.; MVDC(K) = 0.; AFAC(K) = 0.
         ENDIF
         IF (QR3D(K).GE.QSMALL.AND.NR3D(K).GE.NSMALL) THEN
            CALL SOLVE_AFAR(QR3D(K),NR3D(K),LAMR(K),MVDR(K),AFAR(K))
         ELSE
            QR3D(K) = 0.; NR3D(K) = 0.; MVDR(K) = 0.; AFAR(K) = 0.
         ENDIF
!----------------------------- ICE -------------------------------------
         IF (I3M3D(K).LT.ISMALL.OR.FI3D(K).LT.ISMALL) THEN
            I3M3D(K) = 0.; IASPR(K) = 1.; FI3D(K) = 0.
         ENDIF
         IF (AI3D(K).LT.ASMALL) THEN
            AI3D(K) = 0.
         ENDIF
         IF (QI3D(K).GE.QSMALL.AND.NI3D(K).GE.NSMALL) THEN
            CALL SOLVE_AFAI(TK3D(K),RHO(K),QI3D(K),NI3D(K),VI3D(K),    &
                 FI3D(K),AI3D(K),I3M3D(K),ADAGR(K),ZETA(K),LAMI(K),    &
                 AFAI(K),MVDI(K),RHOI(K),AMI(K),BMI(K),AVI(K),BVI(K),  &
                 BEST(K))
            IF (I3M3D(K).GE.ISMALL.AND.FI3D(K).GE.ISMALL) THEN
               I3M0(K) = NI3D(K)*DI0**3.
               IF (MVDI(K).GT.(DI0+1.E-7)) THEN
                  IF (ICE_SHAPE.EQ.0) THEN
                     ZETA(K) = 0.; ADAGR(K) = 1.; IASPR(K) = 1.
                     FI3D(K) = I3M3D(K)
                  ELSEIF (ICE_SHAPE.EQ.1) THEN
                     ZETA(K) = LOG(FI3D(K)/I3M3D(K))/LOG(I3M3D(K)/     &
                               I3M0(K))
                     IF (ZETA(K).GT.0.4) THEN
                        ZETA(K) = 0.4
                        FI3D(K) = (I3M3D(K)/I3M0(K))**ZETA(K)*I3M3D(K)
                     ELSEIF (ZETA(K).LT.(-0.4)) THEN
                        ZETA(K) = -0.4
                        FI3D(K) = (I3M3D(K)/I3M0(K))**ZETA(K)*I3M3D(K)
                     ENDIF
                     ADAGR(K) = (1.+2.*ZETA(K))/(1.-ZETA(K))
                     IASPR(K) = FI3D(K)/I3M3D(K)
                  ENDIF
               ELSE
                  ADAGR(K) = 1.; ZETA(K) = 0.; IASPR(K) = 1.
                  FI3D(K) = I3M3D(K)
               ENDIF
            ENDIF
         ELSE
            QI3D(K) = 0.; NI3D(K) = 0.; AI3D(K) = 0.; MVDI(K) = 0.
            I3M3D(K) = 0.; FI3D(K) = 0.; VI3D(K) = 0.; AFAI(K) = 0.
            IASPR(K) = 1.; ADAGR(K) = 1.; RHOI(K) = 0.
         ENDIF
!-----------------------------------------------------------------------
         IF (AS3D(K).LT.ASMALL) THEN
            AS3D(K) = 0.
         ENDIF
         IF (QS3D(K).GE.QSMALL.AND.NS3D(K).GE.NSMALL) THEN
            CALL SOLVE_AFAS(TK3D(K),RHO(K),QS3D(K),NS3D(K),VS3D(K),    &
                 FS3D(K),AS3D(K),AFAS(K),LAMS(K),MVDS(K),RHOS(K),      &
                 SASPR(K),AMS(K),AVS(K),BVS(K),SMLF(K))
         ELSE
            QS3D(K) = 0.; NS3D(K) = 0.; VS3D(K) = 0.; AS3D(K) = 0.
            MVDS(K) = 0.; RHOS(K) = 0.; AFAS(K) = 0.; FS3D(K) = 0.
            SASPR(K) = 1.
         ENDIF
         IF (AG3D(K).LT.ASMALL) THEN
            AG3D(K) = 0.
         ENDIF
         IF (QG3D(K).GE.QSMALL.AND.NG3D(K).GE.NSMALL) THEN
            CALL SOLVE_AFAG(TK3D(K),RHO(K),QG3D(K),NG3D(K),VG3D(K),    &
                 AG3D(K),LAMG(K),AFAG(K),MVDG(K),RHOG(K),AMG(K),AVG(K),&
                 BVG(K),GMLF(K))
         ELSE
            QG3D(K) = 0.; NG3D(K) = 0.; MVDG(K) = 0.; VG3D(K) = 0.
            AG3D(K) = 0.; AFAG(K) = 0.; RHOG(K) = 0.
         ENDIF
         IF (AH3D(K).LT.ASMALL) THEN
            AH3D(K) = 0.
         ENDIF
         IF (QH3D(K).GE.QSMALL.AND.NH3D(K).GE.NSMALL) THEN
            CALL SOLVE_AFAH(TK3D(K),RHO(K),QH3D(K),NH3D(K),AH3D(K),    &
                 LAMH(K),AFAH(K),MVDH(K),AVH(K),BVH(K))
         ELSE
            QH3D(K) = 0.; NH3D(K) = 0.; MVDH(K) = 0.; AH3D(K) = 0.
            AFAH(K) = 0.
         ENDIF
         IF ((QC3D(K)+QR3D(K)).LT.QSMALL) THEN
            AERO(K,1) = MAX(RLIMIT,AERO(K,1)+AERO(K,3))
            AERO(K,1) = MIN(AERO(K,1),AERO(K,2))
            AERO(K,3) = 0.
         ENDIF
         IF ((QR3D(K)+QS3D(K)+QG3D(K)+QH3D(K)).LT.QSMALL) THEN
            AERO(K,1) = MAX(RLIMIT,AERO(K,1)+AERO(K,4))
            AERO(K,1) = MIN(AERO(K,1),AERO(K,2))
            AERO(K,4) = 0.
         ENDIF
         QDCN3D(K) = AERO(K,1)
         QTCN3D(K) = AERO(K,2)
         QCCN3D(K) = AERO(K,3)
         QRCN3D(K) = AERO(K,4)
         QNIN3D(K) = AERO(K,5)
      ENDDO 

      END SUBROUTINE NTU_MICRO
!======================================================================

!======================================================================
      SUBROUTINE ACTIVA(TK1D,W1D,NC1D,NR1D,ABCD,QACac,QACar,RX0,ZCCNS, &
                        SSW,IAE,DO_FR)
!======================================================================
      IMPLICIT NONE
      INTEGER :: IAE                                                    ! (I) aerosol component
      REAL :: TK1D,NC1D,NR1D,QACac,QACar,NACac,NACar,RX0,SSW
      REAL, INTENT(IN) :: W1D
      REAL*8 :: X1,X3,DMASS,DERF
      REAL, DIMENSION(NCCN) :: ZCCNS
      REAL, DIMENSION(MAER) :: ABCD
      LOGICAL :: DO_FR
      INTEGER :: IM
      REAL :: GAMAA,ALPHA,SSRW,ES0,RX9,RX1,XMASS,RACT,ZCN9,ZCN99,      &
              ZCCN10,ZCCN15,ZCCN20,ZCCN25,RACT10,RACT15,RACT20,RACT25, &
              WR10,WR15,WR20,WR25,WR10TO15,WR15TO20,WR20TO25,WR25PLUS
      REAL, PARAMETER :: RS10 = 10.E-6, RS15 = 15.E-6, RS20 = 20.E-6,  &
                         RS25 = 25.E-6

      QACac = 0.; QACar = 0.; NACac = 0.; NACar = 0.
!------ decide the number of each mode and radius for CCN cut-off ------
      IF (DO_FR) THEN
         CALL FIND_RC0(DBLE(ABCD(1)/ABCD(2)),CNMOD(1,IAE),CNSTD(1,IAE),&
              WMAS(1,IAE),RX0,TBLXA(1,IAE),TBLRC)
         DO_FR = .FALSE.
      ENDIF
      IF (RX0.EQ.1.E-9) RETURN                                          ! RC=RC_MIN, no dry aerosol available
      RX1 = RX0
!----------------- Kohler Curve parameters -----------------------------
      SSRW  = MIN(SSW,0.03)
      ES0   = 0.0761-1.55E-4*(TK1D-TK0C)                               ! [J/M**2]
      ALPHA = 2.*ES0/(RV*TK1D*RHOW)                                    ! [J/m^2]/[J/kg/K]/[K]/[kg/m^3] = [m]
      GAMAA = 4.*ALPHA**3./(27.*BETA1(IAE))                            ! [m^3]
      RX9   = (GAMAA/(SSRW**2.))**THRD                                 ! critical dry radius [m]
      RX9   = MAX(RX9,RXMIN(IAE))                                      ! set lower limit of RX
      IF (RX1.LE.RX9) RETURN
!----------------- RAINDROP ACTIVATION ---------------------------------
      IF (RX1.GT.RS10) THEN                                             ! force all aerosol > RS10 into rain, no matter RX9=??
         DMASS = 0.D+0
         DO IM = 1,NCCN
            X3 = DLNX(RX1,CNMOD(IM,IAE),CNSTD(IM,IAE),3)
            X1 = DLNX(RS10,CNMOD(IM,IAE),CNSTD(IM,IAE),3)
            DMASS = DMASS+DBLE(ZCCNS(IM))/DBLE(RFACT(IM,IAE))*(DERF(X3)-DERF(X1))/2.D+0  ! aerosol into rain
         ENDDO
         XMASS = MAX(0.,REAL(DMASS))
         XMASS = MIN(ABCD(1),XMASS)
         IF (NAERN(IAE).GE.4) THEN
            ABCD(4) = ABCD(4)+XMASS                                     ! update aerosol in rain
         ELSEIF (NAERN(IAE).EQ.3) THEN
            ABCD(3) = ABCD(3)+XMASS                                     ! update aerosol in (cloud+rain)
         ENDIF
         ABCD(1) = ABCD(1)-XMASS                                        ! update dry aerosol
         CALL RSWHITBY(W1D,RS10,ZCCN10,RACT10,BETA1(IAE),ALPHA,ZCCNS,  &
                       CNMOD(1,IAE),CNSTD(1,IAE))
         CALL RSWHITBY(W1D,RS15,ZCCN15,RACT15,BETA1(IAE),ALPHA,ZCCNS,  &
                       CNMOD(1,IAE),CNSTD(1,IAE))
         CALL RSWHITBY(W1D,RS20,ZCCN20,RACT20,BETA1(IAE),ALPHA,ZCCNS,  &
                       CNMOD(1,IAE),CNSTD(1,IAE))
         CALL RSWHITBY(W1D,RS25,ZCCN25,RACT25,BETA1(IAE),ALPHA,ZCCNS,  &
                       CNMOD(1,IAE),CNSTD(1,IAE))
         WR10 = C4PI3W*RACT10**3.                                       ! [kg]
         WR15 = C4PI3W*RACT15**3.                                       ! [kg]
         WR20 = C4PI3W*RACT20**3.                                       ! [kg]
         WR25 = C4PI3W*RACT25**3.                                       ! [kg]
         WR10TO15 = 0.5*(WR10+WR15)*(ZCCN10-ZCCN15)                     ! [kg/kg]
         WR15TO20 = 0.5*(WR15+WR20)*(ZCCN15-ZCCN20)                     ! [kg/kg]
         WR20TO25 = 0.5*(WR20+WR25)*(ZCCN20-ZCCN25)                     ! [kg/kg]
         WR25PLUS = WR25*ZCCN25                                         ! [kg/kg]
         NACar    = ZCCN10                                              ! [1/kg]
         QACar    = WR10TO15+WR15TO20+WR20TO25+WR25PLUS                 ! [kg/kg]
         RX1      = RS10                                                ! [m]
!----------------- end sectionalize rain embryo ------------------------
!----------------- single category rain embryo -------------------------
         NR1D = NR1D+NACar                                              ! [ #/kg]
      ENDIF
!----------------- CLOUD DROP ACTIVATION ------------------------------
      IF (RX1.GT.RX9) THEN
         DMASS = 0.D+0
         DO IM = 1,NCCN
            X3 = DLNX(RX1,CNMOD(IM,IAE),CNSTD(IM,IAE),3)
            X1 = DLNX(RX9,CNMOD(IM,IAE),CNSTD(IM,IAE),3)
            DMASS = DMASS+DBLE(ZCCNS(IM))/DBLE(RFACT(IM,IAE))*(DERF(X3)-DERF(X1))/2.D+0  ! aerosol into rain
         ENDDO
         XMASS = MAX(0.,REAL(DMASS))
         XMASS = MIN(XMASS,ABCD(1))
         IF (NAERN(IAE).GE.3) THEN
            ABCD(3) = ABCD(3)+XMASS                                     ! update wet aerosol in (cloud+rain)
         ENDIF
         ABCD(1) = ABCD(1)-XMASS                                        ! update dry aerosol
         CALL CCNWHITBY(W1D,RS10,RX9,ZCN99,RACT,BETA1(IAE),ALPHA,ZCCNS,&
                        CNMOD(1,IAE),CNSTD(1,IAE))
         CALL CCNWHITBY(W1D,RX1,RX9,ZCN9,RACT,BETA1(IAE),ALPHA,ZCCNS,  &
                        CNMOD(1,IAE),CNSTD(1,IAE))
         NACac = MIN(ZCN9,MAX(ZCN99-NC1D,0.))                           ! Total # should less than ZCN99
         QACac = NACac*C4PI3W*RACT**3.                                  ! [kg/kg]
         NC1D  = NC1D+NACac                                             ! [#/kg]
         RX1   = RX9
      ENDIF
      RX0 = RX1

      END SUBROUTINE ACTIVA
!======================================================================
!
!----------------- FUNCTION FOR CALCULATE lnX --------------------------
      FUNCTION DLNX(RX,XMODE,SIGMA,N)
!======================================================================
      INTEGER :: N
      REAL :: RX,XMODE,SIGMA
      REAL*8 :: DRX,DXMODE,DSTDV,DLNX,DLOG,DSQRT

      DRX    = DBLE(RX)
      DXMODE = DBLE(XMODE)
      DSTDV  = DBLE(SIGMA)
      DLNX   = (DLOG(DRX/DXMODE)-DSTDV**2.*DBLE(N))/(DSQRT(2.D+0)*DSTDV)

      END FUNCTION DLNX
!----------------- FUNCTION FOR CALCULATE lnX --------------------------
      FUNCTION DLNX2(DRX,XMODE,SIGMA,N)
      INTEGER :: N
      REAL :: XMODE,SIGMA
      REAL*8 :: DRX,DXMODE,DSTDV,DLNX2,DLOG,DSQRT

      DXMODE = DBLE(XMODE)
      DSTDV  = DBLE(SIGMA)
      DLNX2  = (DLOG(DRX/DXMODE)-DSTDV**2.*DBLE(N))/(DSQRT(2.D+0)*DSTDV)

      END FUNCTION DLNX2
!======================================================================
!
!======================================================================
      SUBROUTINE RSWHITBY(W1D,RSX,Z,RACT,BETA1,ALPHA,ZCCN,CNMOD,CNSTD)
!======================================================================
      IMPLICIT NONE
!----------------- calculate number of CCN with dry radii > R0 ---------
      INTEGER :: IM
      REAL :: RSX,                                                     &! smallest dry radius of rain embryo (input)
              RWS0,                                                    &! wet radius at S = 1
              RACT,                                                    &! activation radius (output)
              Z,                                                       &! number of CCN that fits the criterion (output)->should be the same as ZCCN (#/kg)?
              BETA1,ALPHA,TEMP,TEMP1,TEMP2,TMPX,TEMP4,ERF
      REAL, INTENT(IN) :: W1D
      REAL, DIMENSION(NCCN) :: ZCCN,CNMOD,CNSTD
!----- integrate the number concentration of aerosols with size > RACT--
      Z = 0.
      DO IM = 1,NCCN
         TEMP1 = ZCCN(IM)/2.      
         TEMP2 = SQRT2*CNSTD(IM) 
         TEMP  = (LOG(RSX/CNMOD(IM)))/TEMP2
         Z     = Z+TEMP1*(1.-ERF(TEMP))
      END DO
      RWS0  = SQRT(BETA1*RSX**3./ALPHA)                                 ! wet radius at S = 1
      TMPX  = LOG10(RSX*1.E6)
      TEMP4 = -0.61425115+(-0.66624878-0.17367658*TMPX)*TMPX
      RACT  = 10.**TEMP4*RWS0*W1D**(-0.11782)

      END SUBROUTINE RSWHITBY
!----------------------------------------------------------------------
      SUBROUTINE CCNWHITBY(W1D,RSX,RS9,ZCN9,RACT,BETA1,ALPHA,ZCCN,     &
                 CNMOD,CNSTD)
!======================================================================
      IMPLICIT NONE
!---- derive number of cloud drop activation according to Whitby's distributions
      INTEGER :: IM
      REAL :: RSX,                                                     &! previously activated smallest dry radius
              RS9,                                                     &! current activated smallest dry radius
              ZCN9,                                                    &! newly activated number concentration (#/mol)-->should be the same as ZCCN (#/kg)?
              RACT,                                                    &! wet radius of the activated drop
              RWS0,                                                    &! wet radius at S = 1
              BETA1,ALPHA,TMPX,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,ERF
      REAL, INTENT(IN) :: W1D
      REAL, DIMENSION(NCCN) :: ZCCN,CNMOD,CNSTD

      RWS0 = SQRT(BETA1*RS9**3./ALPHA)                                  ! wet radius at S = 1
      TMPX = LOG10(RSX*1.E6)
      TMP1 = -0.61425115+(-0.66624878-0.17367658*TMPX)*TMPX
      RACT = 10.**TMP1*RWS0*W1D**(-0.11782)                             ! activation size
      ZCN9 = 0.
      DO IM = 1,NCCN
         TMP2 = ZCCN(IM)/2.
         TMP3 = SQRT2*CNSTD(IM) 
         TMP4 = (LOG(RSX/CNMOD(IM)))/TMP3
         TMP5 = (LOG(RS9/CNMOD(IM)))/TMP3
         TMP6 = TMP2*(ERF(TMP4)-ERF(TMP5))
         ZCN9 = ZCN9+TMP6
      END DO

      END SUBROUTINE CCNWHITBY
!======================================================================
!
!======================================================================
      SUBROUTINE DEACTIVA(NACcv,QC1D,NC1D,DCN,TCN,WCN,RC,ZCCNS,IAE,    &
                          DO_FR)  
!======================================================================
      IMPLICIT NONE
      INTEGER :: IAE                                                    ! (I) aerosol component
      LOGICAL :: DO_FR
      REAL :: NC1D,NACcv,QC1D,DCN,TCN,WCN,                             &! (I/O) wet aerosol mixing ratio [kg/kg]
              RC                                                        ! (I/O) R-CUTOFF [um]
      REAL, DIMENSION(NCCN) :: ZCCNS                                    ! (I) # of each mode in [#/kg]
      INTEGER :: I,IM
      REAL*8 :: DNACcv,                                                &! # of deactivated cloud droplet
                DMASS,                                                 &! mass of deactivated wet aerosol
                DSQRT2,X0,X1,X2,X3,S1,S2,DERF,DEXP,DSQRT
!----------------- deactivation of cloud drops -------------------------
      DSQRT2 = DSQRT(2.D+0)
      DNACcv = DBLE(MIN(NACcv,NC1D))
      IF (QC1D.LT.RSMALL) THEN
         DNACcv = DBLE(NC1D)                                            ! deactivate all cld drop for negligiable QC
         NC1D = RLIMIT
      ENDIF
      DMASS = 0.D+0
      IF (DNACcv.GT.1.D+1) THEN
!------ decide the number of each mode and radius for CCN cut-off ------
         IF (DO_FR) THEN
            CALL FIND_RC0(DBLE(DCN/TCN),CNMOD(1,IAE),CNSTD(1,IAE),     &
                 WMAS(1,IAE),RC,TBLXA(1,IAE),TBLRC)
            DO_FR = .FALSE.
         ENDIF
         S1 = 0.D+0
         DO IM = 1,NCCN
            X0 = DLNX(RC,CNMOD(IM,IAE),CNSTD(IM,IAE),0)
            S1 = S1+DBLE(ZCCNS(IM))*(1.D+0-DERF(X0))*5.D-1              ! # for all with r > RC
         ENDDO
         IF (DNACcv.GT.S1) THEN
            DCN = DCN+WCN
            WCN = 0.
            CALL FIND_RC0(DBLE(DCN/TCN),CNMOD(1,IAE),CNSTD(1,IAE),     &
                 WMAS(1,IAE),RC,TBLXA(1,IAE),TBLRC)
            RETURN
         ENDIF
         S1 = 0.
         DO IM = 1,NCCN
            X0 = DLNX(RC,CNMOD(IM,IAE),CNSTD(IM,IAE),0)
            S1 = S1+DBLE(ZCCNS(IM)/CNSTD(IM,IAE))/DSQRT2*DPDF(X0)       ! # between 0 & RC
         ENDDO
         X1 = DBLE(RC)
         X2 = DBLE(RC)*DEXP(DNACcv/S1)
         DO I = 1,4                                                     ! iteratively fixed dN 3 times should be enough??
            S1 = 0.D+0
            S2 = 0.D+0
            DO IM = 1,NCCN
               X0 = DLNX2(X1,CNMOD(IM,IAE),CNSTD(IM,IAE),0)
               X3 = DLNX2(X2,CNMOD(IM,IAE),CNSTD(IM,IAE),0)
               S1 = S1+DBLE(ZCCNS(IM)/CNSTD(IM,IAE))/DSQRT2*DPDF(X3)
               S2 = S2+DBLE(ZCCNS(IM))*(DERF(X3)-DERF(X0))*5.D-1        ! # between X1 & X2
            ENDDO
            DNACcv = DNACcv-S2
            X1 = X2
            X2 = X2*DEXP(DNACcv/S1)
         ENDDO
         X1 = DBLE(RC)
         X2 = MAX(X1,X2)
         DO IM = 1,NCCN
            X0 = DLNX2(X1,CNMOD(IM,IAE),CNSTD(IM,IAE),3)
            X3 = DLNX2(X2,CNMOD(IM,IAE),CNSTD(IM,IAE),3)
            DMASS = DMASS+DBLE(ZCCNS(IM))/DBLE(RFACT(IM,IAE))*(DERF(X3)-DERF(X0))*5.D-1 ! MASS between X0 & X3
         ENDDO
         DCN = DCN+MIN(WCN,REAL(DMASS))
         WCN = WCN-MIN(WCN,REAL(DMASS))
         RC  = REAL(X2)
      ENDIF

      END SUBROUTINE DEACTIVA
!======================================================================
!
!======================================================================
      SUBROUTINE ICENU(TK1D,P1D,DT,RHO,QV1D,QI1D,NI1D,VI1D,FI1D,AI1D,  &
                 I3M1D,QNIN,IDEPNU)
!======================================================================
      IMPLICIT NONE
      INTEGER :: IDEPNU
      REAL :: TK1D,P1D,DT,RHO,QV1D,QI1D,NI1D,VI1D,FI1D,AI1D,I3M1D,QNIN,&
              QVSI,XXLS,ESI,QNDvi,NNDvi,VNDvi,FNDvi,ANDvi,INDvi,NVI0,  &
              EPA,RGDEP,SFCTNV,ICED,GGDEP,SRI,ARDEP0,ARDEP,COSM1,GEOF1,&
              IJDEP0,IJDEP,DC1,DC2,DC3,TC1D,INR0,DANGLE,DACTE,QVDMAX,  &
              NNUMAX,SSRI,CPM

      QNDvi = 0.;NNDvi = 0.;VNDvi = 0.;FNDvi = 0.;ANDvi = 0.;INDvi = 0.

      TC1D   = TK1D-TK0C
      XXLS   = 3.15E6-2.37E3*TK1D+0.3337E6
      ESI    = MIN(0.99*P1D,POLYSVP(TK1D,1))
      QVSI   = 0.622*ESI/(P1D-ESI)
      SRI    = QV1D/QVSI
      SSRI   = QV1D/QVSI-1.
      QVDMAX = (QV1D-QVSI)/(1.+XXLS**2.*QV1D/(CP*RV*TK1D**2.))          ! SUPERSATURATED WATER FIXED DUE TO TEMPERATURE CHANGE
      NNUMAX = 0.98*QVDMAX/MI0                                          ! INITIATION OF CLOUD ICE CONSTRAINED by 0.98*SUPICE
      ICED   = 916.7-0.175*TC1D-5.E-4*TC1D**2.                          ! ICE DENSITY
      SRI    = QV1D/QVSI                                                ! SATURATION RATIO OVER ICE
      SFCTNV = ((76.1-0.155*TC1D)+(28.5+0.25*TC1D))*1.E-3               ! SURFACE TENSION OVER ICE/VAPOR

      IF (IDEPNU.EQ.0) THEN                                             ! NO DEPOSITION NUCLEATION
         NNDvi = 0.
      ELSEIF (IDEPNU.EQ.1) THEN                                         ! Fletcher (1962)
         NVI0  = 1.E-2*EXP(0.6*(TK0C-MAX(TK1D,2.46E2)))
         NNDvi = MIN(QNIN,NNUMAX,MAX(0.,NVI0/RHO-NI1D))
      ELSEIF (IDEPNU.EQ.2) THEN                                         ! Cooper
         NVI0  = 5.*EXP(0.304*(TK0C-MAX(TK1D,2.46E2)))
         NNDvi = MIN(QNIN,NNUMAX,MAX(0.,NVI0/RHO-NI1D))
      ELSEIF (IDEPNU.EQ.3) THEN                                         ! Huffman (1973)
         NVI0  = 1007.08*SSRI**4.5*1.E3
         NNDvi = MIN(QNIN,NNUMAX,MAX(0.,NVI0/RHO-NI1D))
      ELSEIF (IDEPNU.EQ.4) THEN                                         ! Meyers et al (1992)
         NVI0  = 1.E3*EXP(1.296E1*SSRI-6.39E-1)
         NNDvi = MIN(QNIN,NNUMAX,MAX(0.,NVI0/RHO-NI1D))
      ELSEIF (IDEPNU.EQ.5.AND.INSPEC.EQ.4) THEN                         ! Bacteria(P. syringae):Yankofsky 1981
         NVI0  = QNIN*1.1776*EXP(-89318./((LOG(SRI))**2.*TK1D**3.))
         NNDvi = MIN(QNIN,NNUMAX,MAX(0.,NVI0*DT))
      ELSEIF (IDEPNU.EQ.6) THEN                                         ! DUST (Field 2006)(Asia)
         NVI0  = QNIN*4.337E-1*EXP(-2.E6/((LOG(SRI))**2.*TK1D**3.))
         NNDvi = MIN(QNIN,NNUMAX,MAX(0.,NVI0*DT))
      ELSEIF (IDEPNU.EQ.7) THEN
         IF (INSPEC.EQ.1) THEN                                          ! SOOT
            DC1 = -0.5411955; DC2 = 1.879918; DC3 = 1.607947            ! CHEN Table3
            INR0 = 0.4E-7; DANGLE = 28.; DACTE = -20.E-20               ! HOOSE Table3
         ELSEIF (INSPEC.EQ.2) THEN                                      ! SAHARAN DUST
            DC1 = -0.3353619; DC2 = 1.990979; DC3 = 2.175539            ! CHEN Table3
            INR0 = 1.75E-7; DANGLE = 5.06; DACTE = 3.35E-20             ! CHEN Table6
         ELSEIF (INSPEC.EQ.3) THEN                                      ! ASIAN DUST
            DC1 = -0.3598818; DC2 = 1.982032; DC3 = 2.025390            ! CHEN Table3
            INR0 = 2.E-7; DANGLE = 8.1; DACTE = 1.82E-20                ! CHEN Table6
         ENDIF
         EPA    = ESI*SRI
         RGDEP  = 2.*2.99E-26*SFCTNV/(ICED*BOLTZ*TK1D*LOG(SRI))         ! GERM RADIUS FOR DEPOSITION NUCLEATION
         GGDEP  = C4PI3*SFCTNV*RGDEP**2.                                ! HOMOGENEOUS ENERGY
         ARDEP0 = SFCTNV**5.E-1*(BOLTZ*TK1D)**(-1.5)/1.E13
         ARDEP  = EPA**2.*ARDEP0/ICED
         COSM1  = LOG(1.-COS(DANGLE*PI/1.8E2))
         GEOF1  = MIN(1.,EXP(DC1+DC2*COSM1+DC3*RGDEP/INR0))
         IJDEP0 = EXP((-1.*DACTE-GEOF1*GGDEP)/(BOLTZ*TK1D))
         IJDEP  = ARDEP*INR0**2.*GEOF1**5.E-1*IJDEP0
         NVI0   = QNIN*(1.-EXP(-IJDEP*DT))
         NNDvi  = MIN(QNIN,NNUMAX,MAX(0.,NVI0))
      ELSEIF (IDEPNU.EQ.8) THEN
         NVI0  = MAX(0.,5.94E-5*((TK0C-TK1D)**3.33)*((1.E6/1.E6)**     &! set naer05=1e6/kg (derived to cm-3 at STP below)
                 (0.0264*(TK0C-TK1D)+3.3E-3))*1.E3)                     ! units of #/kg
         NNDvi = MIN(QNIN,NNUMAX,MAX(0.,NVI0/RHO-NI1D))
      ENDIF
      QNDvi = NNDvi*MI0
      VNDvi = QNDvi*iRHOI0
      FNDvi = QNDvi*1.*iAMI0                                            ! ISOMETRICa
      ANDvi = (KCIMIN*NNDvi*INDvi**2.)**THRD
      INDvi = QNDvi*iAMI0
      QNIN  = MAX(0.,QNIN-NNDvi)
      QV1D  = MAX(0.,QV1D-QNDvi)
      QI1D  = MAX(0.,QI1D+QNDvi)
      NI1D  = MAX(0.,NI1D+NNDvi)
      IF (ICE_RHOI.EQ.0.OR.ICE_RHOI.EQ.2) THEN
         VI1D = 0.
      ELSEIF (ICE_RHOI.EQ.1) THEN
         VI1D = MAX(0.,VI1D+VNDvi)
      ENDIF
      FI1D  = MAX(0.,FI1D+FNDvi)
      I3M1D = MAX(0.,I3M1D+INDvi)
      IF (AFAI_3M.EQ.0.OR.AFAI_3M.EQ.2) THEN
         AI1D = 0.
      ELSEIF (AFAI_3M.EQ.1) THEN
         AI1D = MAX(0.,AI1D+ANDvi)
      ENDIF
      CPM    = CP*(1.+0.887*QV1D)
      XXLS   = 3.15E6-2370.*TK1D+0.3337E6
      TK1D   = TK1D+XXLS*QNDvi/CPM

      END SUBROUTINE ICENU
!======================================================================
!
!======================================================================
      SUBROUTINE SEDI_FALL(TK1D,QC1D,QR1D,QI1D,QS1D,QG1D,QH1D,NC1D,    &
                 NR1D,NI1D,NS1D,NG1D,NH1D,VI1D,VS1D,VG1D,FI1D,FS1D,    &
                 AI1D,AS1D,AG1D,AH1D,I3M1D,VTQC,VTQR,VTQI,VTQS,VTQG,   &
                 VTQH,VTNC,VTNR,VTNI,VTNS,VTNG,VTNH,VTVI,VTVS,VTVG,    &
                 VTFI,VTFS,VTAI,VTAS,VTAG,VTAH,VTI3M,SASPR,RHO,SMLF,   &
                 GMLF)
!======================================================================
      IMPLICIT NONE
      REAL :: TK1D,QC1D,QR1D,QI1D,QS1D,QG1D,QH1D,NC1D,NR1D,NI1D,NS1D,  &
              NG1D,NH1D,VI1D,VS1D,VG1D,FI1D,FS1D,AI1D,AS1D,AG1D,AH1D,  &
              I3M1D,VTQC,VTQR,VTQI,VTQS,VTQG,VTQH,VTNC,VTNR,VTNI,VTNS, &
              VTNG,VTNH,VTVI,VTVS,VTVG,VTFI,VTFS,VTAI,VTAS,VTAG,VTAH,  &
              VTI3M,QRHO,RHOAJ,RHO,MVRC,MVRR,FUC,FUR,GUC,GUR,LAMC,LAMR,&
              LAMI,LAMS,LAMG,LAMH,RHOI,RHOS,RHOG,AFAC,AFAR,AFAI,AFAS,  &
              AFAG,AFAH,ADAGR,ZETA3,AMI,BMI,AMS,AMG,AVI,BVI,AVS,BVS,   &
              AVG,BVG,AVH,BVH,ZETA,FSQC,FSQR,FSQI,FSQS,FSQG,FSQH,FSNC, &
              FSNR,FSNI,FSNS,FSNG,FSNH,FSVI,FSVS,FSVG,FSFI,FSAI,FSAS,  &
              FSAG,FSAH,MVDC,MVDR,MVDI,MVDS,MVDG,MVDH,BEST,SASPR,SMLF, &
              GMLF
      REAL, PARAMETER :: AVTC = 8.8462E+02, BVTC = 9.7593E+07, CVTC = -3.4249E+11
      REAL, PARAMETER :: AVTR = 2.1454E+00, BVTR = -2.2812E-04, CVTR = 2.9676E-09
      REAL, PARAMETER :: CQC1 = 2.0901E+01, CQC2 = 9.9111E-01, CQC3 = 4.4182E+00
      REAL, PARAMETER :: CNC1 = 1.8276E+01, CNC2 = 1.0015E+00, CNC3 = 1.9838E+00
      REAL, PARAMETER :: CQR1 = 1.5943E+01, CQR2 = 1.1898E+00, CQR3 = 4.0073E+00
      REAL, PARAMETER :: CNR1 = 9.4791E+00, CNR2 = 9.7607E-01, CNR3 = 1.0858E+00

      QRHO  = SQRT(RHO)
      RHOAJ = (RHOSU/RHO)**0.54
      IF (QC1D.GE.QSMALL) THEN
         CALL SOLVE_AFAC(QC1D,NC1D,LAMC,MVDC,AFAC)
         IF (LIQ_VTC.EQ.0) THEN
            FSQC = EXP(GAMLN(BVC0+BMW+AFAC+1.)-GAMLN(BMW+AFAC+1.)-     &
                   BVC0*LOG(LAMC))
            FSNC = EXP(GAMLN(BVC0+AFAC+1.)-GAMLN(AFAC+1.)-BVC0*        &
                   LOG(LAMC))
            VTQC = RHOAJ*FSQC*AVC0
            VTNC = RHOAJ*FSNC*AVC0
         ELSEIF (LIQ_VTC.EQ.1) THEN
            MVRC = MIN(MAX((QC1D/NC1D/C4PI3W)**THRD,RCMIN),RCMAX)
            GUC  = EXP(EXP(AFU+BFU*(LOG(MVRC))**3.+CFU*QRHO**3.))
            VTQC = EXP(CQC1+CQC2*LOG(NC1D)+CQC3*LOG(MVRC))*GUC/QC1D
            VTNC = EXP(CNC1+CNC2*LOG(NC1D)+CNC3*LOG(MVRC))*GUC/NC1D
         ELSEIF (LIQ_VTC.EQ.2) THEN
            MVRC = MIN(MAX((QC1D/NC1D/C4PI3W)**THRD,RCMIN),RCMAX)
            FUC  = EXP(EXP(AFU+BFU*(LOG(MVRC))**3.+CFU*QRHO**3.))
            VTQC = MVRC*(AVTC+BVTC*MVRC+CVTC*MVRC**2.)*FUC
            VTNC = MVRC*(AVTC+BVTC*MVRC+CVTC*MVRC**2.)*FUC
         ENDIF
         VTQC = MIN(VTQC,VTCMAX)
         VTNC = MIN(VTNC,VTCMAX)
      ENDIF
      IF (QR1D.GE.QSMALL) THEN
         CALL SOLVE_AFAR(QR1D,NR1D,LAMR,MVDR,AFAR)
         IF (LIQ_VTR.EQ.0) THEN
            FSQR = EXP(GAMLN(BVR0+BMW+AFAR+1.)-GAMLN(BMW+AFAR+1.)-     &
                   BVR0*LOG(LAMR))
            FSNR = EXP(GAMLN(BVR0+AFAR+1.)-GAMLN(AFAR+1.)-BVR0*        &
                   LOG(LAMR))
            VTQR = RHOAJ*FSQR*AVR0
            VTNR = RHOAJ*FSNR*AVR0
         ELSEIF (LIQ_VTR.EQ.1) THEN
            MVRR = MIN(MAX((QR1D/NR1D/C4PI3W)**THRD,RRMIN),RRMAX)
            GUR  = EXP(EXP(AFU+BFU*(LOG(MVRR))**3.+CFU*QRHO**3.))
            VTQR = EXP(CQR1+CQR2*LOG(NR1D)+CQR3*LOG(MVRR))*GUR/QR1D
            VTNR = EXP(CNR1+CNR2*LOG(NR1D)+CNR3*LOG(MVRR))*GUR/NR1D
         ELSEIF (LIQ_VTR.EQ.2) THEN
            MVRR = MIN(MAX((QR1D/NR1D/C4PI3W)**THRD,RRMIN),RRMAX)
            FUR  = EXP(EXP(AFU+BFU*(LOG(MVRR))**3.+CFU*QRHO**3.))
            VTQR = EXP(AVTR+BVTR/MVRR+CVTR/(MVRR**2.))*FUR
            VTNR = EXP(AVTR+BVTR/MVRR+CVTR/(MVRR**2.))*FUR
         ENDIF
         VTQR = MIN(VTQR,VTRMAX)
         VTNR = MIN(VTNR,VTRMAX)
      ENDIF
      IF (QI1D.GE.QSMALL) THEN
         CALL SOLVE_AFAI(TK1D,RHO,QI1D,NI1D,VI1D,FI1D,AI1D,I3M1D,ADAGR,&
              ZETA,LAMI,AFAI,MVDI,RHOI,AMI,BMI,AVI,BVI,BEST)
         FSQI = EXP(GAMLN(BVI+BMI+AFAI+1.)-GAMLN(BMI+AFAI+1.)-BVI*     &
                LOG(LAMI))
         FSNI = EXP(GAMLN(BVI+AFAI+1.)-GAMLN(AFAI+1.)-BVI*LOG(LAMI))
         FSVI = EXP(GAMLN(BVI+AFAI+4.)-GAMLN(AFAI+4.)-BVI*LOG(LAMI))
         VTQI = MIN(RHOAJ*FSQI*AVI,VTIMAX)
         VTNI = MIN(RHOAJ*FSNI*AVI,VTIMAX)
         VTVI = MIN(RHOAJ*FSVI*AVI,VTIMAX)
        IF (AI1D.GE.ASMALL) THEN
            FSAI = EXP(GAMLN(BVI+AFAI+3.)-GAMLN(AFAI+3.)-BVI*LOG(LAMI))
            VTAI = MIN(RHOAJ*FSAI*AVI,VTIMAX)
         ENDIF
         IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
            ZETA3 = 3.*(ADAGR-1.)/(ADAGR+2.)
            FSFI  = EXP(GAMLN(BVI+ZETA3+AFAI+4.)-GAMLN(ZETA3+AFAI+4.)- &
                    BVI*LOG(LAMI))
            VTFI  = MIN(RHOAJ*FSFI*AVI,VTIMAX)
            VTI3M = MIN(RHOAJ*FSVI*AVI,VTIMAX)
         ENDIF
      ENDIF
      IF (QS1D.GE.QSMALL) THEN
         CALL SOLVE_AFAS(TK1D,RHO,QS1D,NS1D,VS1D,FS1D,AS1D,AFAS,LAMS,  &
              MVDS,RHOS,SASPR,AMS,AVS,BVS,SMLF)
         FSQS = EXP(GAMLN(BVS+BMS+AFAS+1.)-GAMLN(BMS+AFAS+1.)-BVS*     &
                LOG(LAMS))
         FSNS = EXP(GAMLN(BVS+AFAS+1.)-GAMLN(AFAS+1.)-BVS*LOG(LAMS))
         FSVS = EXP(GAMLN(BVS+AFAS+4.)-GAMLN(AFAS+4.)-BVS*LOG(LAMS))
         VTQS = MIN(RHOAJ*FSQS*AVS,VTSMAX)
         VTNS = MIN(RHOAJ*FSNS*AVS,VTSMAX)
         VTVS = MIN(RHOAJ*FSVS*AVS,VTSMAX)
         VTFS = VTQS
         IF (AS1D.GE.ASMALL) THEN
            FSAS = EXP(GAMLN(BVS+AFAS+3.)-GAMLN(AFAS+3.)-BVS*LOG(LAMS))
            VTAS = MIN(RHOAJ*FSAS*AVS,VTSMAX)
         ENDIF
      ENDIF
      IF (QG1D.GE.QSMALL) THEN
         CALL SOLVE_AFAG(TK1D,RHO,QG1D,NG1D,VG1D,AG1D,LAMG,AFAG,MVDG,  &
              RHOG,AMG,AVG,BVG,GMLF)
         FSQG = EXP(GAMLN(BVG+BMG+AFAG+1.)-GAMLN(BMG+AFAG+1.)-BVG*     &
                LOG(LAMG))
         FSNG = EXP(GAMLN(BVG+AFAG+1.)-GAMLN(AFAG+1.)-BVG*LOG(LAMG))
         FSVG = EXP(GAMLN(BVG+AFAG+4.)-GAMLN(AFAG+4.)-BVG*LOG(LAMG))
         VTQG = MIN(RHOAJ*FSQG*AVG,VTGMAX)
         VTNG = MIN(RHOAJ*FSNG*AVG,VTGMAX)
         VTVG = MIN(RHOAJ*FSVG*AVG,VTGMAX)
         IF (AG1D.GE.ASMALL) THEN
            FSAG = EXP(GAMLN(BVG+AFAG+3.)-GAMLN(AFAG+3.)-BVG*LOG(LAMG))
            VTAG = MIN(RHOAJ*FSAG*AVG,VTGMAX)
         ENDIF
      ENDIF
      IF (QH1D.GE.QSMALL) THEN
         CALL SOLVE_AFAH(TK1D,RHO,QH1D,NH1D,AH1D,LAMH,AFAH,MVDH,AVH,BVH)
         FSQH = EXP(GAMLN(BVH+BMH+AFAH+1.)-GAMLN(BMH+AFAH+1.)-BVH*     &
                LOG(LAMH))
         FSNH = EXP(GAMLN(BVH+AFAH+1.)-GAMLN(AFAH+1.)-BVH*LOG(LAMH))
         VTQH = MIN(RHOAJ*FSQH*AVH,VTHMAX)
         VTNH = MIN(RHOAJ*FSNH*AVH,VTHMAX)
         IF (AH1D.GE.ASMALL) THEN
            FSAH = EXP(GAMLN(BVH+AFAH+3.)-GAMLN(AFAH+3.)-BVH*LOG(LAMH))
            VTAH = MIN(RHOAJ*FSAH*AVH,VTHMAX)
         ENDIF
      ENDIF

      END SUBROUTINE SEDI_FALL
!======================================================================
!
!======================================================================
      SUBROUTINE PTFLUX(Q1D,VT1D,RHO,DZ,NK,DT,DTMN,PRT1D) 
!======================================================================
      IMPLICIT NONE
      INTEGER :: K,NK,NS,NSTEP
      INTEGER, PARAMETER :: MAXSTP = 1000
      REAL :: PRT1D,                                                   &! input/output ; accumulative prcipitation (cm)
              DT,                                                      &! time step for integration (=2*DTMN*60. for leapfrog scheme)
              DTMN                                                      ! time step (min)
      REAL, DIMENSION(NK) :: Q1D,VT1D,RHO,DZ,DQDT
!------------ FOR TENDENCIES DUE TO PRECIPITATION/SEDIMENTATION -------
      NSTEP = 1
      DO K = 1,NK
         NSTEP = MAX(NSTEP,INT(VT1D(K)*DT/DZ(K)+1.))
         IF (NSTEP.GT.80) THEN
             PRINT *,'IN PREFLUX',K,NSTEP,DT,VT1D(K),DZ(K),Q1D(K),RHO(K)
         ENDIF
      ENDDO
      IF (NSTEP.GT.MAXSTP) THEN
         PRINT *,'NSTEP FOR PRECIP. IS: ',NSTEP,VT1D,DT,DZ
         STOP 
      ENDIF 
      DO NS = 1,NSTEP
         DQDT(1) = -VT1D(1)*Q1D(1)/DZ(1)
         DO K = 2,NK,1
            DQDT(K) = (VT1D(K-1)*RHO(K-1)*Q1D(K-1)-                    &
                       VT1D(K)*RHO(K)*Q1D(K))/(DZ(K)*RHO(K))
         ENDDO
         PRT1D = PRT1D+VT1D(NK)*RHO(NK)*Q1D(NK)*DTMN*60./REAL(NSTEP)    ! accumulate precipitation [kg/m^2]; 1 kg/m^2 = 1 mm 
         DO K = 1,NK
            Q1D(K) = Q1D(K)+DQDT(K)*DT/REAL(NSTEP) 
         ENDDO 
      ENDDO 

      END SUBROUTINE PTFLUX
!======================================================================
!
!======================================================================
      SUBROUTINE FLFLUX(Q1D,VT1D,RHO,DZ,NK,DT) 
!======================================================================
      IMPLICIT NONE
      INTEGER :: NSTEP,K,NK,NS
      INTEGER, PARAMETER :: MAXSTP = 1000
      REAL :: DT
      REAL, DIMENSION(NK) :: Q1D,VT1D,RHO,DZ,DQDT
!-- INPUT DATA FOR STEPSIZE RK4 CALCULATION FOR TENDENCIES DUE TO PRECIPITATION/SEDIMENTATION
      NSTEP = 1
      DO K = 1,NK
         NSTEP = MAX(NSTEP,INT(VT1D(K)*DT/DZ(K)+1.))
      ENDDO
      DO NS = 1,NSTEP
         DQDT(1) = -VT1D(1)*Q1D(1)/DZ(1) 
         DO K = 2,NK,1
            DQDT(K) = (VT1D(K-1)*RHO(K-1)*Q1D(K-1)-                    &
                       VT1D(K)*RHO(K)*Q1D(K))/(DZ(K)*RHO(K))
         ENDDO
         DO K = 1,NK
            Q1D(K) = Q1D(K)+DQDT(K)*DT/REAL(NSTEP) 
         ENDDO 
      ENDDO 

      END SUBROUTINE FLFLUX
!======================================================================
!
!======================================================================
      SUBROUTINE SMALL_DT(DT,DTS,SDTS,DTKDT,DQVDT,DPDT,RHO,TK1D,P1D,   &
                 QV1D,QC1D,QR1D,QI1D,QS1D,QG1D,QH1D,NC1D,NR1D,NI1D,    &
                 NS1D,NG1D,NH1D,VI1D,VS1D,VG1D,FI1D,FS1D,AI1D,AS1D,    &
                 AG1D,AH1D,I3M1D,SASPR,XDNC,XDNR,GQCTR,SMLF,GMLF)
!======================================================================
      IMPLICIT NONE
      INTEGER :: I,HID
      REAL :: DT,DTS,SDTS,                                             &! TIME STEP(s)& SUB-TIME STEP(s)& DTS MUST<DT-SDTS
              DTKDT,DQVDT,DPDT                                          ! TEMPERATURE,WATER VAPOR,PRESSURE FORCING
      REAL :: RHO,P1D,CPM,ESW,ESI,QVSW,QVSI,SSRW,SSRI,XXLV,XXLS,XXLF,  &
              QV1D,TK1D,QC1D,QR1D,QI1D,QS1D,QG1D,QH1D,GQCTR,NC1D,NR1D, &
              NI1D,NS1D,NG1D,NH1D,VI1D,VS1D,VG1D,FI1D,FS1D,AI1D,AS1D,  &
              AG1D,AH1D,I3M1D,AFAC,AFAR,AFAI,AFAS,AFAG,AFAH,ADAGR,ZETA,&
              AMI,BMI,AMS,AMG,AVI,BVI,AVS,BVS,AVG,BVG,AVH,BVH,DV,MUA,  &
              KAP,SCN,TC1D,iDT,SASPR,ELCLD,ELDLD,ELCLC,CAPS,SMLF,GMLF
      REAL :: QVTEND,QCTEND,QRTEND,QITEND,QSTEND,QGTEND,QHTEND,RATIO,  &
              QVSOUR,QCSOUR,QRSOUR,QISOUR,QSSOUR,QGSOUR,QHSOUR,QVSINK, &
              QCSINK,QRSINK,QISINK,QSSINK,QGSINK,QHSINK,NCSOUR,NISOUR, &
              NCSINK,NISINK,NSSOUR,NSSINK,FISOUR,FISINK,FSSOUR,FSSINK, &
              VISOUR,VISINK,VSSOUR,VSSINK,VGSOUR,VGSINK,AISOUR,AISINK, &
              ASSOUR,ASSINK,AGSOUR,AGSINK,AHSOUR,AHSINK,IISOUR,IISINK, &
              QVWTEND,QVITEND,DRHIDT,DRHWDT
      REAL :: XDNC                                                      ! deactivated cloud drops (#/kg) = NACcv*DTS
      REAL :: XDNR                                                      ! raindrop evap into cloud drop
      REAL :: DTSI                                                      ! TIME TO REACH RHI=100% (NEGATIVE MEANS THAT TIME IS)
      REAL :: DTSW                                                      ! TIME TO REACH RHW=100% (NEVER HAPPEN OR ALREADY PAST)
      REAL :: QACcv,QACrc,QVDvc,QVDvr,QVDvi,QVDvs,QVDvg,QVDvh,QSBiv,   &
              QSBsv,QSBgv,QSBhv,QEVcv,QEVrv,QEVsv,QEVgv,QEVhv,NACcv,   &
              NACrc,NSBiv,NSBsv,NSBgv,NSBhv,VVDvi,VSBiv,VVDvg,VSBgv,   &
              VEVgv,FVDvi,VVDvs,VSBsv,VEVsv,FSBiv,FVDvs,FSBsv,FEVsv,   &
              AVDvi,ASBiv,AVDvs,ASBsv,AEVsv,AVDvg,ASBgv,AEVgv,AVDvh,   &
              ASBhv,AEVhv,IVDvi,ISBiv
      REAL :: MVDC,MVDR,MVDI,MVDS,MVDG,MVDH,MVRC,MVRR,BTMP,BSTMP,BGTMP,&
              BHTMP,ABW,ABI,RHOAJ,LMVRC,LMVRR,RHOI,RHOS,RHOG,VENQS,    &
              VENQG,VENQH,SUMDEP,VDMAX,EVMAX,SBMAX,SUMCND,SUMEVP,      &
              SUMSUB,H2Z,H4Z,ZETA2,ZETA3,ZETA4,ZETA5,GVHAB,IPH,IPG,GI1,&
              INHGR,BEST,DNIVD,DNSVD,DNGVD,QTMP0,QTMP1,QTMP2,QTMP3,    &
              QTMP4,QTMP5,QTMP6,QTMP7,QTMP8,QTMP9,FTMP0,FTMP1,FTMP2,   &
              FTMP3,FTMP4,FTMP5,FTMP6,FTMP7,FTMP8,FTMP9,ATMP0,ATMP1,   &
              ATMP2,ATMP3,ATMP4,ATMP5,ATMP6,ATMP7,ATMP8,ATMP9,LAMC,    &
              LAMR,LAMI,LAMS,LAMG,LAMH,VENQI,VENFI,VENAI,VENAS,VENAG,  &
              VENAH,VENIC,VENIA,VENQI0,RAT1,RAT2,LLMI,LLMS,LLMG,LLMH
      REAL, PARAMETER :: TORR = 1.E-2                                   ! for 0.01 changes
      REAL, PARAMETER :: AQ1 = 6.6793E+0, BQ1 = 1.0090E+0, CQ1 = 1.4095E+0
      REAL, PARAMETER :: AQ2 = 9.9912E+0, BQ2 =-4.7678E-1, CQ2 =-3.1388E-2
      REAL, PARAMETER :: AN9 =-1.0593E+0, BN9 = 8.9774E-1, CN9 =-2.8403E-1, DN9 = 1.6328E+0
      REAL, PARAMETER :: AN10= 8.2841E+0, BN10= 9.7219E-1, CN10=-5.0808E-1
      REAL :: ZC1,ZC2,ZC3,ZC4,ZP1,ZP2,ZP3,ZP4
      DATA ZC1,ZC2,ZC3,ZC4/0.69509913,-0.46685819,0.30490087,1.62148100/
      DATA ZP1,ZP2,ZP3,ZP4/0.36793126,1.82782890,0.63206874,-1.00164090/

      QVTEND = 0.; QVWTEND = 0.; QVITEND = 0.; QCTEND = 0.; QRTEND = 0.
      QITEND = 0.; QSTEND  = 0.; QGTEND  = 0.; QHTEND = 0.
      MVRC  = 0.;  MVRR  = 0.;   MVDC  = 0.;   MVDR  = 0.;  MVDI  = 0.
      MVDS  = 0.;  MVDG  = 0.;   MVDH  = 0.;   QACcv = 0.;  QACrc = 0.
      QVDvc = 0.;  QVDvr = 0.;   QVDvi = 0.;   QVDvs = 0.;  QVDvg = 0.
      QVDvh = 0.;  QSBiv = 0.;   QSBsv = 0.;   QSBgv = 0.;  QSBhv = 0.
      QEVcv = 0.;  QEVrv = 0.;   QEVsv = 0.;   QEVgv = 0.;  QEVhv = 0.
      NACcv = 0.;  NACrc = 0.;   NSBiv = 0.;   NSBsv = 0.;  NSBgv = 0.
      NSBhv = 0.;  VVDvi = 0.;   VSBiv = 0.;   VVDvs = 0.;  VSBsv = 0.
      VEVsv = 0.;  VVDvg = 0.;   VSBgv = 0.;   VEVgv = 0.;  FVDvi = 0.
      FSBiv = 0.;  FVDvs = 0.;   FSBsv = 0.;   FEVsv = 0.;  IVDvi = 0.
      ISBiv = 0.;  AVDvi = 0.;   ASBiv = 0.;   AVDvs = 0.;  ASBsv = 0.
      AEVsv = 0.;  AVDvg = 0.;   ASBgv = 0.;   AEVgv = 0.;  AVDvh = 0.
      ASBhv = 0.;  AEVhv = 0.

      ESW    = MIN(0.99*P1D,POLYSVP(TK1D,0))
      ESI    = MIN(0.99*P1D,POLYSVP(TK1D,1))
      IF (ESI.GT.ESW) ESI = ESW
      QVSW   = 0.622*ESW/(P1D-ESW)
      QVSI   = 0.622*ESI/(P1D-ESI)
      SSRW   = QV1D/QVSW-1.
      SSRI   = QV1D/QVSI-1.
      XXLV   = 3.1484E6-2370.*TK1D
      XXLS   = 3.15E6-2370.*TK1D+0.3337E6
      XXLF   = 2836310.8-(3.1484E6-2370.*TK1D)
      ELDLD  = (1.+XXLS*XXLS*QV1D/(CP*RV*TK1D**2.))
      ELCLD  = (1.+XXLV*XXLS*QV1D/(CP*RV*TK1D**2.))
      ELCLC  = (1.+XXLV*XXLV*QV1D/(CP*RV*TK1D**2.))
      DRHWDT = (1.+SSRW)*(DPDT/P1D+DQVDT/QV1D-XXLV*DTKDT/(RV*TK1D**2.)) ! caused by large scale QV & P T changes
      DRHIDT = (1.+SSRI)*(DPDT/P1D+DQVDT/QV1D-XXLS*DTKDT/(RV*TK1D**2.)) ! caused by large scale QV & P T changes

      IF (ABS(DRHWDT).LT.RLIMIT) THEN
         DTSW = 1.E7                                                    ! very long time to reach saturation
      ELSE
         DTSW = -SSRW/DRHWDT*1.01                                       ! time to reach saturation
      ENDIF                                                             ! negative means the right time has past
      IF (ABS(DRHIDT).LT.RLIMIT) THEN
         DTSI = 1.E7                                                    ! very long time to reach saturation
      ELSE
         DTSI = -SSRI/DRHIDT*1.01                                       ! time to reach saturation
      ENDIF                                                             ! negative means the right time has past
      IF (TK1D.GT.TK0C) THEN
         IF (DTSW.GT.0.) DTS = MIN(DT,DTSW,DT-SDTS)
      ELSE 
         IF (DTSI.GT.0.) DTS = MIN(DT,DTSI,DT-SDTS)
         IF (DTSW.GT.0.) DTS = MIN(DT,DTS,DTSW,DT-SDTS)
      ENDIF 
      DTS = MAX(DTMIN,MIN(DT,DT-SDTS))

      iDT    = 1./DTS
      TC1D   = TK1D-TK0C
      DV     = 2.11E-5*(TK1D/TK0C)**1.94*(101325./P1D)                  ! DIFFUSIVITY OF WATER VAPOR IN AIR (PRUPPACHER&KLETT:13-3)
      MUA    = 1.72E-5*(393./(TK1D+120.))*(TK1D/TK0C)**1.5              ! DYNAMIC VISCOSITY OF AIR (ROGER AND YAU,P.102)
      KAP    = 2.3971E-2+0.0078E-2*TC1D                                 ! THERMAL CONDUCTIVITY OF AIR
      SCN    = (MUA/(RHO*DV))**THRD                                     ! SCHMIDT NUMBER
      CPM    = CP*(1.+0.887*QV1D)
      ABW    = TK1D*RV/ESW/DV+XXLV*(XXLV/TK1D/RV-1.)/TK1D/KAP
      ABI    = TK1D*RV/ESI/DV+XXLS*(XXLS/TK1D/RV-1.)/TK1D/KAP
      RHOAJ  = (RHOSU/RHO)**0.54
      HID    = MAX(MIN(NINT(ABS(TC1D)/0.25),120),0)
      INHGR  = ITBLE(HID)
      GVHAB  = (INHGR-1.)/(INHGR+2.)+1.

      IF (QC1D.GE.QSMALL) THEN
         CALL SOLVE_AFAC(QC1D,NC1D,LAMC,MVDC,AFAC)
         MVRC  = MIN(MAX((QC1D/NC1D/C4PI3W)**THRD,RCMIN),RCMAX)
         LMVRC = LOG(MVRC)
      ENDIF
      IF (QR1D.GE.QSMALL) THEN
         CALL SOLVE_AFAR(QR1D,NR1D,LAMR,MVDR,AFAR)
         MVRR  = MIN(MAX((QR1D/NR1D/C4PI3W)**THRD,RRMIN),RRMAX)
         LMVRR = LOG(MVRR)
      ENDIF
      IF (QI1D.GE.QSMALL) THEN
         CALL SOLVE_AFAI(TK1D,RHO,QI1D,NI1D,VI1D,FI1D,AI1D,I3M1D,ADAGR,&
              ZETA,LAMI,AFAI,MVDI,RHOI,AMI,BMI,AVI,BVI,BEST)
         GI1  = GAMLN(AFAI+1.)
         LLMI = LOG(LAMI)
      IF (ICE_VENT.EQ.3) THEN
         IF ((ADAGR-1.).GE.SLIMIT) THEN
            BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
            IPH   = 3.*ADAGR/(ADAGR+2.)
            ZETA2 = 2.*(ADAGR-1.)/(ADAGR+2.)
            ZETA3 = 3.*(ADAGR-1.)/(ADAGR+2.)
            ZETA4 = 4.*(ADAGR-1.)/(ADAGR+2.)
            ZETA5 = 5.*(ADAGR-1.)/(ADAGR+2.)
            H2Z   = ZC2*ZETA
            H4Z   = ZC4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            QTMP2 = LLMI*(H2Z+BVI/2.+IPH/2.+1.)
            QTMP3 = LLMI*(H4Z+BVI/2.+IPH/2.+1.)
            QTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPH/2.+AFAI+2.)-GI1-QTMP2)
            QTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPH/2.+AFAI+2.)-GI1-QTMP3)
            QTMP6 = LLMI*(H2Z+BVI+IPH+1.)
            QTMP7 = LLMI*(H4Z+BVI+IPH+1.)
            QTMP8 = EXP(GAMLN(H2Z+BVI+IPH+AFAI+2.)-GI1-QTMP6)
            QTMP9 = EXP(GAMLN(H4Z+BVI+IPH+AFAI+2.)-GI1-QTMP7)
            FTMP0 = EXP(GAMLN(H2Z+ZETA3+AFAI+2.)-GI1-LLMI*(H2Z+ZETA3+  &
                    1.))
            FTMP1 = EXP(GAMLN(H4Z+ZETA3+AFAI+2.)-GI1-LLMI*(H4Z+ZETA3+  &
                    1.))
            FTMP2 = LLMI*(H2Z+BVI/2.+IPH/2.+ZETA3+1.)
            FTMP3 = LLMI*(H4Z+BVI/2.+IPH/2.+ZETA3+1.)
            FTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPH/2.+ZETA3+AFAI+2.)-GI1-    &
                    FTMP2)
            FTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPH/2.+ZETA3+AFAI+2.)-GI1-    &
                    FTMP3)
            FTMP6 = LLMI*(H2Z+BVI+IPH+ZETA3+1.)
            FTMP7 = LLMI*(H4Z+BVI+IPH+ZETA3+1.)
            FTMP8 = EXP(GAMLN(H2Z+ZETA3+BVI+IPH+AFAI+2.)-GI1-FTMP6)
            FTMP9 = EXP(GAMLN(H4Z+ZETA3+BVI+IPH+AFAI+2.)-GI1-FTMP7)
            ATMP0 = EXP(GAMLN(H2Z+AFAI+1.)-GI1-LLMI*H2Z)
            ATMP1 = EXP(GAMLN(H4Z+AFAI+1.)-GI1-LLMI*H4Z)
            ATMP2 = LLMI*(H2Z+BVI/2.+IPH/2.+ZETA3)
            ATMP3 = LLMI*(H4Z+BVI/2.+IPH/2.+ZETA3)
            ATMP4 = EXP(GAMLN(H2Z+BVI/2.+IPH/2.+ZETA3+AFAI+1.)-GI1-    &
                    ATMP2)
            ATMP5 = EXP(GAMLN(H4Z+BVI/2.+IPH/2.+ZETA3+AFAI+1.)-GI1-    &
                    ATMP3)
            ATMP6 = LLMI*(H2Z+BVI+IPH)
            ATMP7 = LLMI*(H4Z+BVI+IPH)
            ATMP8 = EXP(GAMLN(H2Z+BVI+IPH+AFAI+1.)-GI1-ATMP6)
            ATMP9 = EXP(GAMLN(H4Z+BVI+IPH+AFAI+1.)-GI1-ATMP7)
            VENQI = ZC1*QTMP0/DI0**H2Z+ZC3*QTMP1/DI0**H4Z+VENC1*ZC1*   &
                    BTMP*QTMP4/DI0**(H2Z+ZETA)+VENC1*ZC3*BTMP*QTMP5/   &
                    DI0**(H4Z+ZETA)+VENC2*ZC1*BTMP**2.*QTMP8/DI0**(H2Z+&
                    ZETA2)+VENC2*ZC3*BTMP**2.*QTMP9/DI0**(H4Z+ZETA2)
            VENFI = ZC1*FTMP0/DI0**(ZETA3+H2Z)+ZC3*FTMP1/DI0**(ZETA3+  &
                    H4Z)+VENC1*ZC1*BTMP*FTMP4/DI0**(H2Z+ZETA4)+VENC1*  &
                    ZC3*BTMP*FTMP5/DI0**(H4Z+ZETA4)+VENC2*ZC1*BTMP**2.*&
                    FTMP8/DI0**(H2Z+ZETA5)+VENC2*ZC3*BTMP**2.*FTMP9/   &
                    DI0**(H4Z+ZETA5)
            VENAI = ZC1*ATMP0/DI0**H2Z+ZC3*ATMP1/DI0**H4Z+VENC1*ZC1*   &
                    BTMP*ATMP4/DI0**(H2Z+ZETA)+VENC1*ZC3*BTMP*ATMP5/   &
                    DI0**(H4Z+ZETA)+VENC2*ZC1*BTMP**2.*ATMP8/DI0**(H2Z+&
                    ZETA2)+VENC2*ZC3*BTMP**2.*ATMP9/DI0**(H4Z+ZETA2)
         ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
            BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
            IPG   = 3./(ADAGR+2.)
            ZETA2 = 2.*(ADAGR-1.)/(ADAGR+2.)
            ZETA3 = 3.*(ADAGR-1.)/(ADAGR+2.)
            ZETA4 = 2.5*(ADAGR-1.)/(ADAGR+2.)
            H2Z   = ZP2*ZETA
            H4Z   = ZP4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            QTMP2 = LLMI*(H2Z+BVI/2.+IPG/2.+1.)
            QTMP3 = LLMI*(H4Z+BVI/2.+IPG/2.+1.)
            QTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPG/2.+AFAI+2.)-GI1-QTMP2)
            QTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPG/2.+AFAI+2.)-GI1-QTMP3)
            QTMP6 = LLMI*(H2Z+BVI+IPG+1.)
            QTMP7 = LLMI*(H4Z+BVI+IPG+1.)
            QTMP8 = EXP(GAMLN(H2Z+BVI+IPG+AFAI+2.)-GI1-QTMP6)
            QTMP9 = EXP(GAMLN(H4Z+BVI+IPG+AFAI+2.)-GI1-QTMP7)
            FTMP0 = EXP(GAMLN(H2Z+ZETA3+AFAI+2.)-GI1-LLMI*(H2Z+ZETA3+  &
                    1.))
            FTMP1 = EXP(GAMLN(H4Z+ZETA3+AFAI+2.)-GI1-LLMI*(H4Z+ZETA3+  &
                    1.))
            FTMP2 = LLMI*(H2Z+BVI/2.+IPG/2.+ZETA3+1.)
            FTMP3 = LLMI*(H4Z+BVI/2.+IPG/2.+ZETA3+1.)
            FTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPG/2.+ZETA3+AFAI+2.)-GI1-    &
                    FTMP2)
            FTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPG/2.+ZETA3+AFAI+2.)-GI1-    &
                    FTMP3)
            FTMP6 = LLMI*(H2Z+BVI+IPG+ZETA3+1.)
            FTMP7 = LLMI*(H4Z+BVI+IPG+ZETA3+1.)
            FTMP8 = EXP(GAMLN(H2Z+BVI+IPG+ZETA3+AFAI+2.)-GI1-FTMP6)
            FTMP9 = EXP(GAMLN(H4Z+BVI+IPG+ZETA3+AFAI+2.)-GI1-FTMP7)
            ATMP0 = EXP(GAMLN(H2Z+AFAI+1.)-GI1-LLMI*H2Z)
            ATMP1 = EXP(GAMLN(H4Z+AFAI+1.)-GI1-LLMI*H4Z)
            ATMP2 = LLMI*(H2Z+BVI/2.+IPG/2.)
            ATMP3 = LLMI*(H4Z+BVI/2.+IPG/2.)
            ATMP4 = EXP(GAMLN(H2Z+BVI/2.+IPG/2.+AFAI+1.)-GI1-ATMP2)
            ATMP5 = EXP(GAMLN(H4Z+BVI/2.+IPG/2.+AFAI+1.)-GI1-ATMP3)
            ATMP6 = LLMI*(H2Z+BVI+IPG)
            ATMP7 = LLMI*(H4Z+BVI+IPG)
            ATMP8 = EXP(GAMLN(H2Z+BVI+IPG+AFAI+1.)-GI1-ATMP6)
            ATMP9 = EXP(GAMLN(H4Z+BVI+IPG+AFAI+1.)-GI1-ATMP7)
            VENQI = ZP1*QTMP0/DI0**H2Z+ZP3*QTMP1/DI0**H4Z+VENP1*ZP1*   &
                    BTMP*QTMP4/DI0**(H2Z-ZETA/2.)+VENP1*ZP3*BTMP*QTMP5/&
                    DI0**(H4Z-ZETA/2.)+VENP2*ZP1*BTMP**2.*QTMP8/DI0**  &
                    (H2Z-ZETA)+VENP2*ZP3*BTMP**2.*QTMP9/DI0**(H4Z-ZETA)
            VENFI = ZP1*FTMP0/DI0**(H2Z+ZETA3)+ZP3*FTMP1/DI0**(H4Z+    &
                    ZETA3)+VENP1*ZP1*BTMP*FTMP4/DI0**(H2Z+ZETA4)+VENP1*&
                    ZP3*BTMP*FTMP5/DI0**(H4Z+ZETA4)+VENP2*ZP1*BTMP**2.*&
                    FTMP8/DI0**(H2Z+ZETA2)+VENP2*ZP3*BTMP**2.*FTMP9/   &
                    DI0**(H4Z+ZETA2)
            VENAI = ZP1*ATMP0/DI0**H2Z+ZP3*ATMP1/DI0**H4Z+VENP1*ZP1*   &
                    BTMP*ATMP4/DI0**(H2Z-ZETA/2.)+VENP1*ZP3*BTMP*ATMP5/&
                    DI0**(H4Z-ZETA/2.)+VENP2*ZP1*BTMP**2.*ATMP8/DI0**( &
                    H2Z-ZETA)+VENP2*ZP3*BTMP**2.*ATMP9/DI0**(H4Z-ZETA)
         ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
            BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
            QTMP0 = EXP(GAMLN(AFAI+2.)-GI1-LOG(LAMI))
            QTMP1 = LLMI*(1.5+BVI/2.)
            QTMP2 = EXP(GAMLN(BVI/2.+AFAI+2.5)-GI1-QTMP1)
            ATMP1 = LLMI*(0.5+BVI/2.)
            ATMP2 = EXP(GAMLN(BVI/2.+AFAI+1.5)-GI1-ATMP1)
            VENQI = AVSG*QTMP0+BVSG*BTMP*QTMP2
            VENFI = VENQI
            VENAI = AVSG+BVSG*BTMP*ATMP2
         ENDIF
      ELSEIF (ICE_VENT.EQ.1.OR.ICE_VENT.EQ.2) THEN
         IF ((ADAGR-1.).GE.SLIMIT) THEN
            IPH   = 3.*ADAGR/(ADAGR+2.)
            IPG   = 3./(ADAGR+2.)
            ZETA2 = 2.*(ADAGR-1.)/(ADAGR+2.)
            ZETA3 = 3.*(ADAGR-1.)/(ADAGR+2.)
            ZETA4 = 4.*(ADAGR-1.)/(ADAGR+2.)
            ZETA5 = 5.*(ADAGR-1.)/(ADAGR+2.)
            H2Z   = ZC2*ZETA
            H4Z   = ZC4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            FTMP0 = EXP(GAMLN(H2Z+ZETA3+AFAI+2.)-GI1-LLMI*(H2Z+ZETA3+  &
                    1.))
            FTMP1 = EXP(GAMLN(H4Z+ZETA3+AFAI+2.)-GI1-LLMI*(H4Z+ZETA3+  &
                    1.))
            ATMP0 = EXP(GAMLN(H2Z+AFAI+1.)-GI1-LLMI*H2Z)
            ATMP1 = EXP(GAMLN(H4Z+AFAI+1.)-GI1-LLMI*H4Z)
            IF (BEST.LE.1.) THEN
               BTMP  = SCN**2.*(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(H2Z+BVI+IPH+1.)
               QTMP3 = LLMI*(H4Z+BVI+IPH+1.)
               QTMP4 = EXP(GAMLN(H2Z+BVI+IPH+AFAI+2.)-GI1-QTMP2)
               QTMP5 = EXP(GAMLN(H4Z+BVI+IPH+AFAI+2.)-GI1-QTMP3)
               FTMP2 = LLMI*(H2Z+BVI+IPH+ZETA3+1.)
               FTMP3 = LLMI*(H4Z+BVI+IPH+ZETA3+1.)
               FTMP4 = EXP(GAMLN(H2Z+BVI+IPH+ZETA3+AFAI+2.)-GI1-FTMP2)
               FTMP5 = EXP(GAMLN(H4Z+BVI+IPH+ZETA3+AFAI+2.)-GI1-FTMP3)
               ATMP2 = LLMI*(H2Z+BVI+IPH)
               ATMP3 = LLMI*(H4Z+BVI+IPH)
               ATMP4 = EXP(GAMLN(H2Z+BVI+IPH+AFAI+1.)-GI1-ATMP2)
               ATMP5 = EXP(GAMLN(H4Z+BVI+IPH+AFAI+1.)-GI1-ATMP3)
               VENQI = AVIS*ZC1*QTMP0/DI0**H2Z+AVIS*ZC3*QTMP1/DI0**H4Z+&
                       BVIS*ZC1*BTMP*QTMP4/DI0**(H2Z+ZETA2)+BVIS*ZC3*  &
                       BTMP*QTMP5/DI0**(H4Z+ZETA2)
               VENFI = AVIS*ZC1*FTMP0/DI0**(ZETA3+H2Z)+AVIS*ZC3*FTMP1/ &
                       DI0**(H4Z+ZETA3)+BVIS*ZC1*BTMP*FTMP4/DI0**(H2Z+ &
                       ZETA5)+BVIS*ZC3*BTMP*FTMP5/DI0**(H4Z+ZETA5)
               VENAI = AVIS*ZC1*ATMP0/DI0**H2Z+AVIS*ZC3*ATMP1/DI0**H4Z+&
                       BVIS*ZC1*BTMP*ATMP4/DI0**(H2Z+ZETA2)+BVIS*ZC3*  &
                       BTMP*ATMP5/DI0**(H4Z+ZETA2)
               IF (ICE_VENT.EQ.2) THEN
                  QTMP6 = LLMI*(BVI+IPH)
                  QTMP7 = LLMI*(BVI+IPG)
                  QTMP8 = EXP(GAMLN(BVI+IPH+AFAI+1.)-GI1-QTMP6)
                  QTMP9 = EXP(GAMLN(BVI+IPG+AFAI+1.)-GI1-QTMP7)
                  VENIC = AVIS+BVIS*BTMP*QTMP8/DI0**ZETA2
                  VENIA = AVIS+BVIS*BTMP*QTMP9*DI0**ZETA
                  INHGR = INHGR*VENIC/VENIA
                  GVHAB = (INHGR-1.)/(INHGR+2.)+1.
               ENDIF
            ELSEIF (BEST.GT.1.) THEN
               BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(H2Z+BVI/2.+IPH/2.+1.)
               QTMP3 = LLMI*(H4Z+BVI/2.+IPH/2.+1.)
               QTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPH/2.+AFAI+2.)-GI1-QTMP2)
               QTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPH/2.+AFAI+2.)-GI1-QTMP3)
               FTMP2 = LLMI*(H2Z+BVI/2.+IPH/2.+ZETA3+1.)
               FTMP3 = LLMI*(H4Z+BVI/2.+IPH/2.+ZETA3+1.)
               FTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPH/2.+ZETA3+AFAI+2.)-GI1- &
                       FTMP2)
               FTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPH/2.+ZETA3+AFAI+2.)-GI1- &
                       FTMP3)
               ATMP2 = LLMI*(H2Z+BVI/2.+IPH/2.+ZETA3)
               ATMP3 = LLMI*(H4Z+BVI/2.+IPH/2.+ZETA3)
               ATMP4 = EXP(GAMLN(H2Z+BVI/2.+IPH/2.+ZETA3+AFAI+1.)-GI1- &
                       ATMP2)
               ATMP5 = EXP(GAMLN(H4Z+BVI/2.+IPH/2.+ZETA3+AFAI+1.)-GI1- &
                       ATMP3)
               VENQI = AVSG*ZC1*QTMP0/DI0**H2Z+AVSG*ZC3*QTMP1/DI0**H4Z+&
                       BVSG*ZC1*BTMP*QTMP4/DI0**(H2Z+ZETA)+BVSG*ZC3*   &
                       BTMP*QTMP5/DI0**(H4Z+ZETA)
               VENFI = AVSG*ZC1*FTMP0/DI0**(ZETA3+H2Z)+AVSG*ZC3*FTMP1/ &
                       DI0**(H4Z+ZETA3)+BVSG*ZC1*BTMP*FTMP4/DI0**(H2Z+ &
                       ZETA4)+BVSG*ZC3*BTMP*FTMP5/DI0**(H4Z+ZETA4)
               VENAI = AVSG*ZC1*ATMP0/DI0**H2Z+AVSG*ZC3*ATMP1/DI0**H4Z+&
                       BVSG*ZC1*BTMP*ATMP4/DI0**(H2Z+ZETA)+BVSG*ZC3*   &
                       BTMP*ATMP5/DI0**(H4Z+ZETA)
               IF (ICE_VENT.EQ.2) THEN
                  QTMP6 = LLMI*(BVI/2.+IPH/2.)
                  QTMP7 = LLMI*(BVI/2.+IPG/2.)
                  QTMP8 = EXP(GAMLN(BVI/2.+IPH/2.+AFAI+1.)-GI1-QTMP6)
                  QTMP9 = EXP(GAMLN(BVI/2.+IPG/2.+AFAI+1.)-GI1-QTMP7)
                  VENIC = AVSG+BVSG*BTMP*QTMP8/DI0**ZETA
                  VENIA = AVSG+BVSG*BTMP*QTMP9*DI0**(ZETA/2.)
                  INHGR = INHGR*VENIC/VENIA
                  GVHAB = (INHGR-1.)/(INHGR+2.)+1.
               ENDIF
            ENDIF
         ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
            IPG   = 3./(ADAGR+2.)
            IPH   = 3.*ADAGR/(ADAGR+2.)
            ZETA2 = 2.*(ADAGR-1.)/(ADAGR+2.)
            ZETA3 = 3.*(ADAGR-1.)/(ADAGR+2.)
            ZETA4 = 2.5*(ADAGR-1.)/(ADAGR+2.)
            H2Z   = ZP2*ZETA
            H4Z   = ZP4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            FTMP0 = EXP(GAMLN(H2Z+ZETA3+AFAI+2.)-GI1-LLMI*(H2Z+ZETA3+  &
                    1.))
            FTMP1 = EXP(GAMLN(H4Z+ZETA3+AFAI+2.)-GI1-LLMI*(H4Z+ZETA3+  &
                    1.))
            ATMP0 = EXP(GAMLN(H2Z+AFAI+1.)-GI1-LLMI*H2Z)
            ATMP1 = EXP(GAMLN(H4Z+AFAI+1.)-GI1-LLMI*H4Z)
            IF (BEST.LE.1.) THEN
               BTMP  = SCN**2.*(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(H2Z+BVI+IPG+1.)
               QTMP3 = LLMI*(H4Z+BVI+IPG+1.)
               QTMP4 = EXP(GAMLN(H2Z+BVI+IPG+AFAI+2.)-GI1-QTMP2)
               QTMP5 = EXP(GAMLN(H4Z+BVI+IPG+AFAI+2.)-GI1-QTMP3)
               FTMP2 = LLMI*(H2Z+BVI+IPG+ZETA3+1.)
               FTMP3 = LLMI*(H4Z+BVI+IPG+ZETA3+1.)
               FTMP4 = EXP(GAMLN(H2Z+BVI+IPG+ZETA3+AFAI+2.)-GI1-FTMP2)
               FTMP5 = EXP(GAMLN(H4Z+BVI+IPG+ZETA3+AFAI+2.)-GI1-FTMP3)
               ATMP2 = LLMI*(H2Z+BVI+IPG)
               ATMP3 = LLMI*(H4Z+BVI+IPG)
               ATMP4 = EXP(GAMLN(H2Z+BVI+IPG+AFAI+1.)-GI1-ATMP2)
               ATMP5 = EXP(GAMLN(H4Z+BVI+IPG+AFAI+1.)-GI1-ATMP3)
               VENQI = AVIS*ZP1*QTMP0/DI0**H2Z+AVIS*ZP3*QTMP1/DI0**H4Z+&
                       BVIS*ZP1*BTMP*QTMP4/DI0**(H2Z-ZETA)+BVIS*ZP3*   &
                       BTMP*QTMP5/DI0**(H4Z-ZETA)
               VENFI = AVIS*ZP1*FTMP0/DI0**(ZETA3+H2Z)+AVIS*ZP3*FTMP1/ &
                       DI0**(H4Z+ZETA3)+BVIS*ZP1*BTMP*FTMP4/DI0**(H2Z+ &
                       ZETA2)+BVIS*ZP3*BTMP*FTMP5/DI0**(H4Z+ZETA2)
               VENAI = AVIS*ZP1*ATMP0/DI0**H2Z+AVIS*ZP3*ATMP1/DI0**H4Z+&
                       BVIS*ZP1*BTMP*ATMP4/DI0**(H2Z-ZETA)+BVIS*ZP3*   &
                       BTMP*ATMP5/DI0**(H4Z-ZETA)
               IF (ICE_VENT.EQ.2) THEN
                  QTMP6 = LLMI*(BVI+IPH)
                  QTMP7 = LLMI*(BVI+IPG)
                  QTMP8 = EXP(GAMLN(BVI+IPH+AFAI+1.)-GI1-QTMP6)
                  QTMP9 = EXP(GAMLN(BVI+IPG+AFAI+1.)-GI1-QTMP7)
                  VENIC = AVIS+BVIS*BTMP*QTMP8/DI0**ZETA2
                  VENIA = AVIS+BVIS*BTMP*QTMP9*DI0**ZETA
                  INHGR = INHGR*VENIC/VENIA
                  GVHAB = (INHGR-1.)/(INHGR+2.)+1.
               ENDIF
            ELSEIF (BEST.GT.1.) THEN
               BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(H2Z+BVI/2.+IPG/2.+1.)
               QTMP3 = LLMI*(H4Z+BVI/2.+IPG/2.+1.)
               QTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPG/2.+AFAI+2.)-GI1-QTMP2)
               QTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPG/2.+AFAI+2.)-GI1-QTMP3)
               FTMP2 = LLMI*(H2Z+BVI/2.+IPG/2.+ZETA3+1.)
               FTMP3 = LLMI*(H4Z+BVI/2.+IPG/2.+ZETA3+1.)
               FTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPG/2.+ZETA3+AFAI+2.)-GI1- &
                       FTMP2)
               FTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPG/2.+ZETA3+AFAI+2.)-GI1- &
                       FTMP3)
               ATMP2 = LLMI*(H2Z+BVI/2.+IPG/2.)
               ATMP3 = LLMI*(H4Z+BVI/2.+IPG/2.)
               ATMP4 = EXP(GAMLN(H2Z+BVI/2.+IPG/2.+AFAI+1.)-GI1-ATMP2)
               ATMP5 = EXP(GAMLN(H4Z+BVI/2.+IPG/2.+AFAI+1.)-GI1-ATMP3)
               VENQI = AVSG*ZP1*QTMP0/DI0**H2Z+AVSG*ZP3*QTMP1/DI0**H4Z+&
                       BVSG*ZP1*BTMP*QTMP4/DI0**(H2Z-ZETA/2.)+BVSG*ZP3*&
                       BTMP*QTMP5/DI0**(H4Z-ZETA/2.)
               VENFI = AVSG*ZP1*FTMP0/DI0**(ZETA3+H2Z)+AVSG*ZP3*FTMP1/ &
                       DI0**(H4Z+ZETA3)+BVSG*ZP1*BTMP*FTMP4/DI0**(H2Z+ &
                       ZETA4)+BVSG*ZP3*BTMP*FTMP5/DI0**(H4Z+ZETA4)
               VENAI = AVSG*ZP1*ATMP0/DI0**H2Z+AVSG*ZP3*ATMP1/DI0**H4Z+&
                       BVSG*ZP1*BTMP*ATMP4/DI0**(H2Z-ZETA/2.)+BVSG*ZP3*&
                       BTMP*ATMP5/DI0**(H4Z-ZETA/2.)
               IF (ICE_VENT.EQ.2) THEN
                  QTMP6 = LLMI*(BVI/2.+IPH/2.)
                  QTMP7 = LLMI*(BVI/2.+IPG/2.)
                  QTMP8 = EXP(GAMLN(BVI/2.+IPH/2.+AFAI+1.)-GI1-QTMP6)
                  QTMP9 = EXP(GAMLN(BVI/2.+IPG/2.+AFAI+1.)-GI1-QTMP7)
                  VENIC = AVSG+BVSG*BTMP*QTMP8/DI0**ZETA
                  VENIA = AVSG+BVSG*BTMP*QTMP9*DI0**(ZETA/2.)
                  INHGR = INHGR*VENIC/VENIA
                  GVHAB = (INHGR-1.)/(INHGR+2.)+1.
               ENDIF
            ENDIF
         ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
            QTMP0 = EXP(GAMLN(AFAI+2.)-GI1-LOG(LAMI))
            IF (BEST.LE.1.) THEN
               BTMP  = SCN**2.*(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(2.+BVI)
               QTMP3 = EXP(GAMLN(BVI+AFAI+3.)-GI1-QTMP2)
               ATMP2 = LLMI*(1.+BVI)
               ATMP3 = EXP(GAMLN(BVI+AFAI+2.)-GI1-ATMP2)
               VENQI = AVIS*QTMP0+BVIS*BTMP*QTMP3
               VENFI = VENQI
               VENAI = AVIS+BVIS*BTMP*ATMP3
            ELSEIF (BEST.GT.1.) THEN
               BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(1.5+BVI/2.)
               QTMP3 = EXP(GAMLN(BVI/2.+AFAI+2.5)-GI1-QTMP2)
               ATMP2 = LLMI*(0.5+BVI/2.)
               ATMP3 = EXP(GAMLN(BVI/2.+AFAI+1.5)-GI1-ATMP2)
               VENQI = AVSG*QTMP0+BVSG*BTMP*QTMP3
               VENFI = VENQI
               VENAI = AVSG+BVSG*BTMP*ATMP3
            ENDIF
         ENDIF
      ELSEIF (ICE_VENT.EQ.0) THEN
         IF ((ADAGR-1.).GE.SLIMIT) THEN
            H2Z   = ZC2*ZETA
            H4Z   = ZC4*ZETA
            ZETA3 = 3.*(ADAGR-1.)/(ADAGR+2.)
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            FTMP0 = EXP(GAMLN(H2Z+ZETA3+AFAI+2.)-GI1-LLMI*(H2Z+ZETA3+  &
                    1.))
            FTMP1 = EXP(GAMLN(H4Z+ZETA3+AFAI+2.)-GI1-LLMI*(H4Z+ZETA3+  &
                    1.))
            ATMP0 = EXP(GAMLN(H2Z+AFAI+1.)-GI1-LLMI*H2Z)
            ATMP1 = EXP(GAMLN(H4Z+AFAI+1.)-GI1-LLMI*H4Z)
            VENQI = ZC1*QTMP0/DI0**H2Z+ZC3*QTMP1/DI0**H4Z
            VENFI = ZC1*FTMP0/DI0**(H2Z+ZETA3)+ZC3*FTMP1/DI0**(H4Z+    &
                    ZETA3)
            VENAI = ZC1*ATMP0/DI0**H2Z+ZC3*ATMP1/DI0**H4Z
         ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
            H2Z   = ZP2*ZETA
            H4Z   = ZP4*ZETA
            ZETA3 = 3.*(ADAGR-1.)/(ADAGR+2.)
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            FTMP0 = EXP(GAMLN(H2Z+ZETA3+AFAI+2.)-GI1-LLMI*(H2Z+ZETA3+  &
                    1.))
            FTMP1 = EXP(GAMLN(H4Z+ZETA3+AFAI+2.)-GI1-LLMI*(H4Z+ZETA3+  &
                    1.))
            ATMP0 = EXP(GAMLN(H2Z+AFAI+1.)-GI1-LLMI*H2Z)
            ATMP1 = EXP(GAMLN(H4Z+AFAI+1.)-GI1-LLMI*H4Z)
            VENQI = ZP1*QTMP0/DI0**H2Z+ZP3*QTMP1/DI0**H4Z
            VENFI = ZP1*FTMP0/DI0**(H2Z+ZETA3)+ZP3*FTMP1/DI0**(H4Z+    &
                    ZETA3)
            VENAI = ZP1*ATMP0/DI0**H2Z+ZP3*ATMP1/DI0**H4Z
         ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
            VENQI = EXP(GAMLN(AFAI+2.)-GI1-LOG(LAMI))
            VENFI = VENQI
            VENAI = 1.
         ENDIF
      ENDIF                                                             ! ICE_VENT
      IF (ICE_SHAPE.EQ.1.AND.AFAI.LE.20.) THEN
         IF (ICE_VENT.EQ.3) THEN
            BTMP   = SCN*SQRT(AVI*RHOAJ/MUA)
            QTMP0  = EXP(GAMLN(AFAI+2.)-GI1-LOG(LAMI))
            QTMP1  = LLMI*(1.5+BVI/2.)
            QTMP2  = EXP(GAMLN(BVI/2.+AFAI+2.5)-GI1-QTMP1)
            ATMP1  = LLMI*(0.5+BVI/2.)
            ATMP2  = EXP(GAMLN(BVI/2.+AFAI+1.5)-GI1-ATMP1)
            VENQI0 = AVSG*QTMP0+BVSG*BTMP*QTMP2
         ELSEIF (ICE_VENT.EQ.1.OR.ICE_VENT.EQ.2) THEN
            QTMP0 = EXP(GAMLN(AFAI+2.)-GI1-LOG(LAMI))
            IF (BEST.LE.1.) THEN
               BTMP   = SCN**2.*(AVI*RHOAJ/MUA)
               QTMP2  = LLMI*(2.+BVI)
               QTMP3  = EXP(GAMLN(BVI+AFAI+3.)-GI1-QTMP2)
               ATMP2  = LLMI*(1.+BVI)
               ATMP3  = EXP(GAMLN(BVI+AFAI+2.)-GI1-ATMP2)
               VENQI0 = AVIS*QTMP0+BVIS*BTMP*QTMP3
            ELSEIF (BEST.GT.1.) THEN
               BTMP   = SCN*SQRT(AVI*RHOAJ/MUA)
               QTMP2  = LLMI*(1.5+BVI/2.)
               QTMP3  = EXP(GAMLN(BVI/2.+AFAI+2.5)-GI1-QTMP2)
               ATMP2  = LLMI*(0.5+BVI/2.)
               ATMP3  = EXP(GAMLN(BVI/2.+AFAI+1.5)-GI1-ATMP2)
               VENQI0 = AVSG*QTMP0+BVSG*BTMP*QTMP3
            ENDIF
         ELSEIF (ICE_VENT.EQ.0) THEN
            VENQI0 = EXP(GAMLN(AFAI+2.)-GI1-LOG(LAMI))
         ENDIF
      ENDIF                                                             ! DI0_CORRECTION
      ELSE
         ADAGR = 1.
         ZETA = 0.
      ENDIF                                                             ! QI SMALL
      IF (QS1D.GE.QSMALL) THEN
         CALL SOLVE_AFAS(TK1D,RHO,QS1D,NS1D,VS1D,FS1D,AS1D,AFAS,LAMS,  &
              MVDS,RHOS,SASPR,AMS,AVS,BVS,SMLF)
         LLMS  = LOG(LAMS)
         BSTMP = SCN*SQRT(AVS*RHOAJ/MUA)
         QTMP1 = LLMS*(1.5+BVS/2.)
         QTMP2 = EXP(GAMLN(AFAS+2.)-GAMLN(AFAS+1.)-LOG(LAMS))
         QTMP3 = EXP(GAMLN(BVS/2.+AFAS+2.5)-GAMLN(AFAS+1.)-QTMP1)
         ATMP1 = LLMS*(0.5+BVS/2.)
         ATMP2 = EXP(GAMLN(BVS/2.+AFAS+1.5)-GAMLN(AFAS+1.)-ATMP1)
         CAPS  = ZP1*SASPR**(ZP2/3.)+ZP3*SASPR**(ZP4/3.)
         VENQS = AVSG*QTMP2*CAPS+BVSG*BSTMP*QTMP3*CAPS
         VENAS = AVSG*CAPS+BVSG*BSTMP*ATMP2*CAPS
      ENDIF
      IF (QG1D.GE.QSMALL) THEN
         CALL SOLVE_AFAG(TK1D,RHO,QG1D,NG1D,VG1D,AG1D,LAMG,AFAG,MVDG,  &
              RHOG,AMG,AVG,BVG,GMLF)
         LLMG  = LOG(LAMG)
         BGTMP = SCN*SQRT(AVG*RHOAJ/MUA)
         QTMP1 = LLMG*(1.5+BVG/2.)
         QTMP2 = EXP(GAMLN(AFAG+2.)-GAMLN(AFAG+1.)-LOG(LAMG))
         QTMP3 = EXP(GAMLN(BVG/2.+AFAG+2.5)-GAMLN(AFAG+1.)-QTMP1)
         ATMP1 = LLMG*(0.5+BVG/2.)
         ATMP2 = EXP(GAMLN(BVG/2.+AFAG+1.5)-GAMLN(AFAG+1.)-ATMP1)
         VENQG = AVSG*QTMP2+BVSG*BGTMP*QTMP3
         VENAG = AVSG+BVSG*BGTMP*ATMP2
      ENDIF
      IF (QH1D.GE.QSMALL) THEN
         CALL SOLVE_AFAH(TK1D,RHO,QH1D,NH1D,AH1D,LAMH,AFAH,MVDH,AVH,BVH)
         LLMH = LOG(LAMH)
         IF (HAIL_VENT.EQ.0) THEN
            QTMP1 = LLMH*(1.5+BVH/2.)
            QTMP2 = EXP(GAMLN(AFAH+2.)-GAMLN(AFAH+1.)-LOG(LAMH))
            QTMP3 = EXP(GAMLN(BVH/2.+AFAH+2.5)-GAMLN(AFAH+1.)-QTMP1)
            ATMP1 = LLMH*(0.5+BVH/2.)
            ATMP2 = EXP(GAMLN(BVH/2.+AFAH+1.5)-GAMLN(AFAH+1.)-ATMP1)
            VENQH = AVRH*QTMP2+BVRH*BHTMP*QTMP3
            VENAH = AVRH+BVRH*BHTMP*ATMP2
         ELSEIF (HAIL_VENT.EQ.1) THEN
            QTMP1 = LLMH*(1.5+BVH/2.)
            QTMP2 = EXP(GAMLN(AFAH+2.)-GAMLN(AFAH+1.)-LOG(LAMH))
            QTMP3 = EXP(GAMLN(BVH/2.+AFAH+2.5)-GAMLN(AFAH+1.)-QTMP1)
            QTMP4 = LLMH*(2.+BVH)
            QTMP5 = EXP(GAMLN(BVH+AFAH+3.)-GAMLN(AFAH+1.)-QTMP4)
            ATMP1 = LLMH*(0.5+BVH/2.)
            ATMP2 = EXP(GAMLN(BVH/2.+AFAH+1.5)-GAMLN(AFAH+1.)-ATMP1)
            ATMP3 = LLMH*(1.+BVH)
            ATMP4 = EXP(GAMLN(BVH+AFAH+2.)-GAMLN(AFAH+1.)-ATMP3)
            VENQH = QTMP2+VENH1*BHTMP*QTMP3+VENH2*BHTMP**2.*QTMP5
            VENAH = 1.+VENH1*BHTMP*ATMP2+VENH2*BHTMP**2.*ATMP4
         ENDIF
      ENDIF
!----------------- WARM RAIN PROCESSES ---------------------------------
      IF (QC1D.GE.QSMALL) THEN
         IF (SSRW.GT.RSMALL) THEN
            QVDvc = SSRW*EXP(AQ1+BQ1*LOG(NC1D)+CQ1*LMVRC)/ABW           ! diffusion growth of cloud drops
         ELSEIF (SSRW.LT.(-1.*RSMALL)) THEN
            QEVcv = MIN(SSRW*EXP(AQ1+BQ1*LOG(NC1D)+CQ1*LMVRC)/ABW,0.)   ! diffusion growth of cloud drops
            NACcv = -1.*EXP(AN9+BN9*LOG(NC1D)+CN9*LMVRC+DN9*LOG(-SSRW)) ! deactivation of cloud drop into CN
            QACcv = NACcv*4.1888E-15                                    ! deactivation of cloud drop into CN
            QACcv = MAX(QACcv,-1.*QC1D*iDT)
            NACcv = MAX(NACcv,-1.*NC1D*iDT)
         ENDIF
      ENDIF
      IF (QR1D.GE.QSMALL) THEN
         IF (SSRW.GT.RSMALL) THEN
            QVDvr = SSRW*NR1D*EXP(AQ2+(BQ2+CQ2*LMVRR)*LMVRR**2.)/ABW    ! diffusion growth of raindrops
         ELSEIF (SSRW.LT.(-1.*RSMALL)) THEN
            QEVrv = SSRW*NR1D*EXP(AQ2+(BQ2+CQ2*LMVRR)*LMVRR**2.)/ABW    ! diffusion growth of raindrops
            NACrc = SSRW*EXP(AN10+BN10*LOG(NR1D)+CN10*LMVRR)/ABW        ! rain drop evaporate into cloud drop
            QACrc = NACrc*5.236E-10                                     ! rain drop evaporate into cloud drop
            QEVrv = MIN(QEVrv,0.)
            QACrc = MAX(QACrc,-1.*QR1D*iDT)
            NACrc = MAX(NACrc,-1.*NR1D*iDT)
         ENDIF
      ENDIF
      IF ((QVDvc+QVDvr).GE.QSMALL) THEN
         VDMAX  = (QV1D-QVSW)/(1.+XXLV**2.*QV1D/(CPM*RV*TK1D**2.))*iDT
         SUMCND = QVDvc+QVDvr
         IF (SUMCND.GT.VDMAX.AND.VDMAX.GE.QSMALL) THEN
            RATIO = MIN(1.,VDMAX/(SUMCND+QSMALL))
            QVDvc = QVDvc*RATIO
            QVDvr = QVDvr*RATIO
         ENDIF
      ENDIF
!-------------- DEPOSITION/SUBLIMATION TERMS FOR QI,QS,QG,QH -----------
      IF (TK1D.LT.TK0C) THEN
         IF (QI1D.GE.QSMALL) THEN
            DNIVD = RHOI0*EXP(-3.*MAX((QV1D-QVSI)-5.E-5,0.)/INHGR)
            DNIVD = MIN(MAX(DNIVD,RHOIMIN),RHOIMAX)
         IF (ICE_SHAPE.EQ.1.AND.AFAI.LE.20.) THEN
            RAT1  = MIN(1.,MAX(0.,ABS(GAMMP(AFAI+1.,DI0*LAMI))))
            RAT2  = 1.-RAT1
            QVDvi = 2.*PI*NI1D*(RAT2*VENQI+RAT1*VENQI0)*SSRI/ABI
            IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
               FVDvi = 12.*NI1D*(RAT2*VENFI*GVHAB+RAT1*VENQI0)*SSRI/   &
                       ABI/DNIVD
               IVDvi = 12.*NI1D*(RAT2*VENQI+RAT1*VENQI0)*SSRI/ABI/DNIVD
            ENDIF
         ELSE
            QVDvi = 2.*PI*NI1D*VENQI*SSRI/ABI
            IF (AI1D.GE.ASMALL) THEN
               AVDvi = 8.*NI1D*VENAI*SSRI/ABI/DNIVD
            ENDIF
            IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
               FVDvi = 12.*NI1D*GVHAB*VENFI*SSRI/ABI/DNIVD
               IVDvi = 12.*NI1D*VENQI*SSRI/ABI/DNIVD
            ENDIF
         ENDIF
            VVDvi = QVDvi/DNIVD
            IF (QVDvi.LT.0.) THEN
               QSBiv = MIN(QVDvi,0.)
               VSBiv = MIN(VVDvi,0.)
               QVDvi = 0.
               VVDvi = 0.
               IF (AI1D.GE.ASMALL) THEN
                  ASBiv = MIN(AVDvi,0.)
                  AVDvi = 0.
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FSBiv = MIN(FVDvi,0.)
                  ISBiv = MIN(IVDvi,0.)
                  FVDvi = 0.
                  IVDvi = 0.
               ENDIF
            ENDIF
         ENDIF                                                          ! QISMALL
         IF (QS1D.GE.QSMALL) THEN
            DNSVD = RHOS*EXP(-3.*MAX((QV1D-QVSI)-5.E-5,0.)/1.)
            DNSVD = MIN(MAX(DNSVD,RHOIMIN),RHOIMAX)
            QVDvs = 2.*PI*NS1D*VENQS*SSRI/ABI
            IF (ICE_RHOS.EQ.1.AND.AGG_SHAPE.EQ.1) THEN
               VVDvs = QVDvs/DNSVD
            ELSE
               VVDvs = QVDvs/RHOS0
            ENDIF
            IF (AGG_SHAPE.EQ.1) THEN
               FVDvs = QVDvs*MAX(0.,-1.*LOG(SASPR))
            ENDIF
            IF (AS1D.GE.ASMALL) THEN
               AVDvs = 8.*NS1D*VENAS*SSRI/ABI/DNSVD
            ENDIF
            IF (QVDvs.LT.0.) THEN
               QSBsv = MIN(QVDvs,0.)
               VSBsv = MIN(VVDvs,0.)
               FSBsv = MIN(FVDvs,0.)
               QVDvs = 0.
               VVDvs = 0.
               FVDvs = 0.
               IF (AS1D.GE.ASMALL) THEN
                  ASBsv = MIN(AVDvs,0.)
                  AVDvs = 0.
               ENDIF
            ENDIF
         ENDIF
         IF (QG1D.GE.QSMALL) THEN
            DNGVD = RHOG*EXP(-3.*MAX((QV1D-QVSI)-5.E-5,0.))
            DNGVD = MIN(MAX(DNGVD,RHOIMIN),RHOG0)
            QVDvg = 2.*PI*NG1D*VENQG*SSRI/ABI
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               VVDvg = QVDvg/DNGVD
            ENDIF
            IF (AG1D.GE.ASMALL) THEN
               AVDvg = 8.*NG1D*VENAG*SSRI/ABI/DNGVD
            ENDIF
            IF (QVDvg.LT.0.) THEN
               QSBgv = MIN(QVDvg,0.)
               QVDvg = 0.
               VSBgv = MIN(VVDvg,0.)
               VVDvg = 0.
               IF (AG1D.GE.ASMALL) THEN
                  ASBgv = MIN(AVDvg,0.)
                  AVDvg = 0.
               ENDIF
            ENDIF
         ENDIF
         IF (QH1D.GE.QSMALL) THEN
            QVDvh = 2.*PI*NH1D*VENQH*SSRI/ABI
            IF (AH1D.GE.ASMALL) THEN
               AVDvh = 8.*NH1D*VENAH*SSRI/ABI/RHOH
            ENDIF
            IF (QVDvh.LT.0.) THEN
               QSBhv = MIN(QVDvh,0.)
               QVDvh = 0.
               IF (AH1D.GE.ASMALL) THEN
                  ASBhv = MIN(AVDvh,0.)
                  AVDvh = 0.
               ENDIF
            ENDIF
         ENDIF
         IF ((QVDvi+QVDvs+QVDvg+QVDvh).GE.QSMALL) THEN
            VDMAX = (QV1D-QVSI)/(1.+XXLS**2.*QV1D/(CPM*RV*TK1D**2.))*iDT
            SUMDEP = QVDvi+QVDvs+QVDvg+QVDvh
            IF (SUMDEP.GT.VDMAX.AND.VDMAX.GE.QSMALL) THEN
               RATIO = MIN(1.,VDMAX/(SUMDEP+QSMALL))
               QVDvi = QVDvi*RATIO
               QVDvs = QVDvs*RATIO
               QVDvg = QVDvg*RATIO
               QVDvh = QVDvh*RATIO
               VVDvi = VVDvi*RATIO
               VVDvs = VVDvs*RATIO
               VVDvg = VVDvg*RATIO
               FVDvi = FVDvi*RATIO
               FVDvs = FVDvs*RATIO
               IVDvi = IVDvi*RATIO
               AVDvi = AVDvi*RATIO
               AVDvs = AVDvs*RATIO
               AVDvg = AVDvg*RATIO
               AVDvh = AVDvh*RATIO
            ENDIF
         ENDIF
         IF ((QSBiv+QSBsv+QSBgv+QSBhv).LT.0.) THEN
            SBMAX = (QV1D-QVSI)/(1.+XXLS**2.*QV1D/(CPM*RV*TK1D**2.))*iDT
            SUMSUB = QSBiv+QSBsv+QSBgv+QSBhv
            IF (SBMAX.LT.0..AND.SUMSUB.LT.SBMAX*0.9999) THEN
               QSBiv = QSBiv*MIN(1.,0.9999*SBMAX/SUMSUB)
               QSBsv = QSBsv*MIN(1.,0.9999*SBMAX/SUMSUB)
               QSBgv = QSBgv*MIN(1.,0.9999*SBMAX/SUMSUB)
               QSBhv = QSBhv*MIN(1.,0.9999*SBMAX/SUMSUB)
               VSBiv = VSBiv*MIN(1.,0.9999*SBMAX/SUMSUB)
               VSBsv = VSBsv*MIN(1.,0.9999*SBMAX/SUMSUB)
               VSBgv = VSBgv*MIN(1.,0.9999*SBMAX/SUMSUB)
               FSBiv = FSBiv*MIN(1.,0.9999*SBMAX/SUMSUB)
               FSBsv = FSBsv*MIN(1.,0.9999*SBMAX/SUMSUB)
               ISBiv = ISBiv*MIN(1.,0.9999*SBMAX/SUMSUB)
               ASBiv = ASBiv*MIN(1.,0.9999*SBMAX/SUMSUB)
               ASBsv = ASBsv*MIN(1.,0.9999*SBMAX/SUMSUB)
               ASBgv = ASBgv*MIN(1.,0.9999*SBMAX/SUMSUB)
               ASBhv = ASBhv*MIN(1.,0.9999*SBMAX/SUMSUB)
            ENDIF
         ENDIF
      ELSE                                                              ! TK1D.GE.TK0C
         IF (QS1D.GE.QSMALL.AND.SSRW.LT.-1.*RSMALL) THEN
            QEVsv = 2.*PI*NS1D*VENQS*SSRW/ABW
            QEVsv = MAX(MIN(QEVsv,0.),-1.*QS1D*iDT)
            IF (ICE_RHOS.EQ.1) THEN
               VEVsv = MAX(MIN(QEVsv/RHOS,0.),-1.*VS1D*iDT)
            ENDIF
            IF (AGG_SHAPE.EQ.1) THEN
               FEVsv = QEVsv*MAX(0.,-1.*LOG(SASPR))
            ENDIF
            IF (AS1D.GE.ASMALL) THEN
               AEVsv = 8.*NS1D*VENAS*SSRW/ABW/RHOS
               AEVsv = MAX(MIN(AEVsv,0.),-1.*AS1D*iDT)
            ENDIF
         ENDIF
         IF (QG1D.GE.QSMALL.AND.SSRW.LT.-1.*RSMALL) THEN
            QEVgv = 2.*PI*NG1D*VENQG*SSRW/ABW
            QEVgv = MAX(MIN(QEVgv,0.),-1.*QG1D*iDT)
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               VEVgv = MAX(MIN(QEVgv/RHOG,0.),-1.*VG1D*iDT)
            ENDIF
            IF (AG1D.GE.ASMALL) THEN
               AEVgv = 8.*NG1D*VENAG*SSRW/ABW/RHOG
               AEVgv = MAX(MIN(AEVgv,0.),-1.*AG1D*iDT)
            ENDIF
         ENDIF
         IF (QH1D.GE.QSMALL.AND.SSRW.LT.-1.*RSMALL) THEN
            QEVhv = 2.*PI*NH1D*VENQH*SSRW/ABW
            QEVhv = MAX(MIN(QEVhv,0.),-1.*QH1D*iDT)
            IF (AH1D.GE.ASMALL) THEN
               AEVhv = 8.*NH1D*VENAH*SSRW/ABW/RHOH
               AEVhv = MAX(MIN(AEVhv,0.),-1.*AH1D*iDT)
            ENDIF
         ENDIF
         IF ((QEVcv+QEVrv+QEVsv+QEVgv+QEVhv).LT.0.) THEN
            EVMAX = (QV1D-QVSW)/(1.+XXLV**2.*QV1D/(CPM*RV*TK1D**2.))*iDT
            SUMEVP = QEVcv+QEVrv+QEVsv+QEVgv+QEVhv
            IF (EVMAX.LT.0..AND.SUMEVP.LT.EVMAX*0.9999) THEN
               QEVcv = QEVcv*MIN(1.,0.9999*EVMAX/SUMEVP)
               QEVrv = QEVrv*MIN(1.,0.9999*EVMAX/SUMEVP)
               QEVsv = QEVsv*MIN(1.,0.9999*EVMAX/SUMEVP)
               QEVgv = QEVgv*MIN(1.,0.9999*EVMAX/SUMEVP)
               QEVhv = QEVhv*MIN(1.,0.9999*EVMAX/SUMEVP)
               VEVsv = VEVsv*MIN(1.,0.9999*EVMAX/SUMEVP)
               VEVgv = VEVgv*MIN(1.,0.9999*EVMAX/SUMEVP)
               FEVsv = FEVsv*MIN(1.,0.9999*EVMAX/SUMEVP)
               AEVsv = AEVsv*MIN(1.,0.9999*EVMAX/SUMEVP)
               AEVgv = AEVgv*MIN(1.,0.9999*EVMAX/SUMEVP)
               AEVhv = AEVhv*MIN(1.,0.9999*EVMAX/SUMEVP)
            ENDIF
         ENDIF
      ENDIF                                                             ! TEMPERATURE LOOP
      IF (AI1D.LT.ASMALL) THEN
         AI1D = 0.
      ENDIF
      IF (AS1D.LT.ASMALL) THEN
         AS1D = 0.
      ENDIF
      IF (AG1D.LT.ASMALL) THEN
         AG1D = 0.
      ENDIF
      IF (AH1D.LT.ASMALL) THEN
         AH1D = 0.
      ENDIF
      IF (I3M1D.LT.ISMALL.OR.FI1D.LT.ISMALL) THEN
         FI1D = 0.; I3M1D = 0.; ADAGR = 1.; ZETA = 0.
      ENDIF

      QVWTEND = -QACcv-QVDvc-QVDvr-QEVcv-QEVrv-QEVsv-QEVgv-QEVhv
      QVITEND = -QVDvi-QVDvs-QVDvg-QVDvh-QSBiv-QSBsv-QSBgv-QSBhv
      QVTEND  = QVWTEND+QVITEND
      QCTEND  = QACcv-QACrc+QVDvc+QEVcv
      QRTEND  = QACrc+QVDvr+QEVrv
      QITEND  = QVDvi+QSBiv
      QSTEND  = QVDvs+QEVsv+QSBsv
      QGTEND  = QVDvg+QEVgv+QSBgv
      QHTEND  = QVDvh+QEVhv+QSBhv

      DTS = MIN(TORR*MAX(QV1D,1.E-7)/MAX(ABS(QVTEND),QSMALL),          &! limited by d(QV)
                TORR*MAX(QC1D,QLIMIT)/MAX(ABS(QCTEND),QSMALL),         &! limited by d(QC)
                TORR*MAX(QR1D,QLIMIT)/MAX(ABS(QRTEND),QSMALL),         &! limited by d(QR)
                TORR*MAX(QI1D,QLIMIT)/MAX(ABS(QITEND),QSMALL),         &! limited by d(QI)
                TORR*MAX(QS1D,QLIMIT)/MAX(ABS(QSTEND),QSMALL),         &! limited by d(QS)
                TORR*MAX(QG1D,QLIMIT)/MAX(ABS(QGTEND),QSMALL),         &! limited by d(QG)
                TORR*MAX(QH1D,QLIMIT)/MAX(ABS(QHTEND),QSMALL),         &! limited by d(QH)
                DT,DTS,DT-SDTS)
      DTS = MIN(DT-SDTS,MAX(DTMIN,DTS))

!--- CALCULATE TIME TO REACH RHW=100%, RHI=100% 
!--- NEGATIVE VALUE MEANS TIME HAS PAST AND WOULD NEVER REACH 100% IN THE FUTURE

      DTSI = -(QVWTEND*ELCLD+QVITEND*ELDLD)                             ! rate of saturation vapor over ice change   
      DTSW = -(QVWTEND*ELCLC+QVITEND*ELCLD)                             ! rate of saturation vapor over liquid change 
      IF (ABS(DTSI).LT.RLIMIT) THEN
         DTSI = 1.E5                                                    ! AN ARBITARY VALUE
      ELSE
         DTSI = SSRI*QV1D/((1.+SSRI)*DTSI)                              ! time to reach saturation over ice 
      ENDIF
      IF (ABS(DTSW).LT.RLIMIT) THEN
         DTSW = 1.E5                                                    ! AN ARBITARY VALUE
      ELSE
         DTSW = SSRW*QV1D/((1.+SSRW)*DTSW)                              ! time to reach saturation over liquid 
      ENDIF

      IF (TK1D.GT.TK0C) THEN 
         IF(DTSW.GT.0.) DTS = MAX(MIN(DTS,DTSW),DTMIN)
      ELSE 
         IF(DTSI.GT.0.) DTS = MAX(MIN(DTS,DTSI),DTMIN)
         IF(DTSW.GT.0.) DTS = MAX(MIN(DTS,DTSW),DTMIN)
      ENDIF

      QVSOUR = QV1D+(-QACcv-QEVcv-QEVrv-QEVsv-QEVgv-QEVhv-QSBiv-QSBsv- &
               QSBgv-QSBhv)*DTS
      QVSINK = (QVDvc+QVDvr+QVDvi+QVDvs+QVDvg+QVDvh)*DTS
      IF (QVSINK.GT.QVSOUR.AND.QVSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QVSOUR/(QVSINK+QSMALL))
         QVDvc = QVDvc*RATIO; QVDvr = QVDvr*RATIO
         QVDvi = QVDvi*RATIO; QVDvs = QVDvs*RATIO
         QVDvg = QVDvg*RATIO; QVDvh = QVDvh*RATIO
      ENDIF
      QCSOUR = QC1D+(QVDvc-QACrc)*DTS
      QCSINK = (-QACcv-QEVcv)*DTS
      IF (QCSINK.GT.QCSOUR.AND.QCSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QCSOUR/(QCSINK+QSMALL))
         QACcv = QACcv*RATIO; QEVcv = QEVcv*RATIO
      ENDIF
      QRSOUR = QR1D+QVDvr*DTS
      QRSINK = (-QACrc-QEVrv)*DTS
      IF (QRSINK.GT.QRSOUR.AND.QRSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QRSOUR/(QRSINK+QSMALL))
         QACrc = QACrc*RATIO; QEVrv = QEVrv*RATIO
      ENDIF
      QISOUR = QI1D+QVDvi*DTS
      QISINK = (-QSBiv)*DTS
      IF (QISINK.GT.QISOUR.AND.QISOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QISOUR/(QISINK+QSMALL))
         QSBiv = QSBiv*RATIO
      ENDIF
      QSSOUR = QS1D+QVDvs*DTS
      QSSINK = (-QEVsv-QSBsv)*DTS
      IF (QSSINK.GT.QSSOUR.AND.QSSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QSSOUR/(QSSINK+QSMALL))
         QSBsv = QSBsv*RATIO; QEVsv = QEVsv*RATIO
      ENDIF
      QGSOUR = QG1D+QVDvg*DTS
      QGSINK = (-QSBgv-QEVgv)*DTS
      IF (QGSINK.GT.QGSOUR.AND.QGSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QGSOUR/(QGSINK+QSMALL))
         QSBgv = QSBgv*RATIO; QEVgv = QEVgv*RATIO
      ENDIF
      QHSOUR = QH1D+QVDvh*DTS
      QHSINK = (-QSBhv-QEVhv)*DTS
      IF (QHSINK.GT.QHSOUR.AND.QHSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QHSOUR/(QHSINK+QSMALL))
         QSBhv = QSBhv*RATIO; QEVhv = QEVhv*RATIO
      ENDIF
      IF (QI1D.GE.QSMALL) THEN
         NSBiv = MIN(QSBiv*NI1D/QI1D,0.)
      ENDIF
      IF (QS1D.GE.QSMALL) THEN
         NSBsv = MIN(QSBsv*NS1D/QS1D,0.)
      ENDIF
      IF (QG1D.GE.QSMALL) THEN
         NSBgv = MIN(QSBgv*NG1D/QG1D,0.)
      ENDIF
      IF (QH1D.GE.QSMALL) THEN
         NSBhv = MIN(QSBhv*NH1D/QH1D,0.)
      ENDIF
      NCSOUR = NC1D+(-NACrc)*DTS
      NCSINK = (-NACcv)*DTS
      IF (NCSINK.GT.NCSOUR.AND.NCSOUR.GE.NSMALL) THEN
         RATIO = MIN(1.,NCSOUR/(NCSINK+NSMALL))
         NACcv = NACcv*RATIO
      ENDIF
      NISOUR = NI1D
      NISINK = (-NSBiv)*DTS
      IF (NISINK.GT.NISOUR.AND.NISOUR.GE.NSMALL) THEN
         RATIO = MIN(1.,NISOUR/(NISINK+NSMALL))
         NSBiv = NSBiv*RATIO
      ENDIF
      NSSOUR = NS1D
      NSSINK = (-NSBsv)*DTS
      IF (NSSINK.GT.NSSOUR.AND.NSSOUR.GE.NSMALL) THEN
         RATIO = MIN(1.,NSSOUR/(NSSINK+NSMALL))
         NSBsv = NSBsv*RATIO
      ENDIF
      VISOUR = VI1D+VVDvi*DTS
      VISINK = (-VSBiv)*DTS
      IF (VISINK.GT.VISOUR.AND.VISOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,VISOUR/(VISINK+ISMALL))
         VSBiv = VSBiv*RATIO
      ENDIF
      VSSOUR = VS1D+VVDvs*DTS
      VSSINK = (-VEVsv-VSBsv)*DTS
      IF (VSSINK.GT.VSSOUR.AND.VSSOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,VSSOUR/(VSSINK+ISMALL))
         VSBsv = VSBsv*RATIO; VEVsv = VEVsv*RATIO
      ENDIF
      VGSOUR = VG1D+VVDvg*DTS
      VGSINK = (-VSBgv-VEVgv)*DTS
      IF (VGSINK.GT.VGSOUR.AND.VGSOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,VGSOUR/(VGSINK+ISMALL))
         VSBgv = VSBgv*RATIO; VEVgv = VEVgv*RATIO
      ENDIF
      FISOUR = FI1D+FVDvi*DTS
      FISINK = (-FSBiv)*DTS
      IF (FISINK.GT.FISOUR.AND.FISOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,FISOUR/(FISINK+ISMALL))
         FSBiv = FSBiv*RATIO
      ENDIF
      FSSOUR = FS1D+FVDvs*DTS
      FSSINK = (-FEVsv-FSBsv)*DTS
      IF (FSSINK.GT.FSSOUR.AND.FSSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,FSSOUR/(FSSINK+QSMALL))
         FSBsv = FSBsv*RATIO; FEVsv = FEVsv*RATIO
      ENDIF
      AISOUR = AI1D+AVDvi*DTS
      AISINK = (-ASBiv)*DTS
      IF (AISINK.GT.AISOUR.AND.AISOUR.GE.ASMALL) THEN
         RATIO = MIN(1.,AISOUR/(AISINK+ASMALL))
         ASBiv = ASBiv*RATIO
      ENDIF
      ASSOUR = AS1D+AVDvs*DTS
      ASSINK = (-ASBsv-AEVsv)*DTS
      IF (ASSINK.GT.ASSOUR.AND.ASSOUR.GE.ASMALL) THEN
         RATIO = MIN(1.,ASSOUR/(ASSINK+ASMALL))
         ASBsv = ASBsv*RATIO; AEVsv = AEVsv*RATIO
      ENDIF
      AGSOUR = AG1D+AVDvg*DTS
      AGSINK = (-ASBgv-AEVgv)*DTS
      IF (AGSINK.GT.AGSOUR.AND.AGSOUR.GE.ASMALL) THEN
         RATIO = MIN(1.,AGSOUR/(AGSINK+ASMALL))
         ASBgv = ASBgv*RATIO; AEVgv = AEVgv*RATIO
      ENDIF
      AHSOUR = AH1D+AVDvh*DTS
      AHSINK = (-ASBhv-AEVhv)*DTS
      IF (AHSINK.GT.AHSOUR.AND.AHSOUR.GE.ASMALL) THEN
         RATIO = MIN(1.,AHSOUR/(AHSINK+ASMALL))
         ASBhv = ASBhv*RATIO; AEVhv = AEVhv*RATIO
      ENDIF
      IISOUR = I3M1D+IVDvi*DTS
      IISINK = (-ISBiv)*DTS
      IF (IISINK.GT.IISOUR.AND.IISOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,IISOUR/(IISINK+ISMALL))
         ISBiv = ISBiv*RATIO
      ENDIF
!----------------- update values of variables -------------------------
      QV1D = MAX(0.,QV1D+(-QACcv-QVDvc-QVDvr-QEVcv-QEVrv-QEVsv-QEVgv-  &
             QEVhv-QVDvi-QVDvs-QVDvg-QVDvh-QSBiv-QSBsv-QSBgv-QSBhv)*DTS)
      QC1D = MAX(0.,QC1D+(QACcv-QACrc+QVDvc+QEVcv)*DTS)
      QR1D = MAX(0.,QR1D+(QACrc+QVDvr+QEVrv)*DTS)
      QI1D = MAX(0.,QI1D+(QVDvi+QSBiv)*DTS)
      QS1D = MAX(0.,QS1D+(QVDvs+QSBsv+QEVsv)*DTS)
      QG1D = MAX(0.,QG1D+(QVDvg+QSBgv+QEVgv)*DTS)
      QH1D = MAX(0.,QH1D+(QVDvh+QSBhv+QEVhv)*DTS)
      NC1D = MAX(0.,NC1D+(NACcv-NACrc)*DTS)
      NR1D = MAX(0.,NR1D+(NACrc)*DTS)
      NI1D = MAX(0.,NI1D+(NSBiv)*DTS)
      NS1D = MAX(0.,NS1D+(NSBsv)*DTS)
      NG1D = MAX(0.,NG1D+(NSBgv)*DTS)
      NH1D = MAX(0.,NH1D+(NSBhv)*DTS)
      IF (ICE_RHOI.EQ.0.OR.ICE_RHOI.EQ.2) THEN
         VI1D = 0.
      ELSEIF (ICE_RHOI.EQ.1) THEN
         VI1D = MAX(0.,VI1D+(VVDvi+VSBiv)*DTS)
      ENDIF
      IF (ICE_RHOS.EQ.0.OR.ICE_RHOS.EQ.2) THEN
         VS1D = 0.
      ELSEIF (ICE_RHOS.EQ.1) THEN
         VS1D = MAX(0.,VS1D+(VVDvs+VSBsv+VEVsv)*DTS)
      ENDIF
      IF (ICE_RHOG.EQ.0) THEN
         VG1D = 0.
      ELSEIF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
         VG1D = MAX(0.,VG1D+(VVDvg+VSBgv+VEVgv)*DTS)
      ENDIF
      FI1D = MAX(0.,FI1D+(FVDvi+FSBiv)*DTS)
      IF (AGG_SHAPE.EQ.0) THEN
         FS1D = 0.
      ELSEIF (AGG_SHAPE.EQ.1) THEN
         FS1D = MAX(0.,FS1D+(FVDvs+FSBsv+FEVsv)*DTS)
      ENDIF
      I3M1D = MAX(0.,I3M1D+(IVDvi+ISBiv)*DTS)
      IF (AFAI_3M.EQ.0.OR.AFAI_3M.EQ.2) THEN
         AI1D = 0.
      ELSEIF (AFAI_3M.EQ.1) THEN
         AI1D = MAX(0.,AI1D+(AVDvi+ASBiv)*DTS)
      ENDIF
      IF (AFAS_3M.EQ.0.OR.AFAS_3M.EQ.2) THEN
         AS1D = 0.
      ELSEIF (AFAS_3M.EQ.1) THEN
         AS1D = MAX(0.,AS1D+(AVDvs+ASBsv+AEVsv)*DTS)
      ENDIF
      IF (AFAG_3M.EQ.0.OR.AFAG_3M.EQ.2) THEN
         AG1D = 0.
      ELSEIF (AFAG_3M.EQ.1) THEN
         AG1D = MAX(0.,AG1D+(AVDvg+ASBgv+AEVgv)*DTS)
      ENDIF
      IF (AFAH_3M.EQ.0.OR.AFAH_3M.EQ.2) THEN
         AH1D = 0.
      ELSEIF (AFAH_3M.EQ.1) THEN
         AH1D = MAX(0.,AH1D+(AVDvh+ASBhv+AEVhv)*DTS)
      ENDIF
      XDNC  = XDNC+NACcv*DTS
      XDNR  = XDNR+NACrc*DTS
      GQCTR = QACrc*DTS
      TK1D  = TK1D+((QACcv+QACrc+QVDvc+QVDvr+QEVcv+QEVrv+QEVsv+QEVgv+  &
              QEVhv)*XXLV+(QVDvi+QVDvs+QVDvg+QVDvh+QSBiv+QSBsv+QSBgv+  &
              QSBhv)*XXLS)/CPM*DTS

      IF (QC1D.GE.QSMALL.AND.NC1D.GE.NSMALL) THEN
         MVDC = (QC1D*iAMW/NC1D)**THRD
         IF (MVDC.GT.DCR) THEN
            QR1D = QR1D+QC1D
            NR1D = NR1D+NC1D
            QC1D = 0.; NC1D = 0.
         ENDIF
      ENDIF
      IF (QR1D.GE.QSMALL.AND.NR1D.GE.NSMALL) THEN
         MVDR = (QR1D*iAMW/NR1D)**THRD
         IF (MVDR.LT.DCR) THEN
            QC1D = QC1D+QR1D
            NC1D = NC1D+NR1D
            QR1D = 0.; NR1D = 0.
         ENDIF
      ENDIF
      IF (QH1D.GE.QSMALL.AND.NH1D.GE.NSMALL) THEN
         MVDH = (QH1D*iAMH/NH1D)**THRD
         IF (MVDH.LT.DHMIN) THEN
            QG1D = QG1D+QH1D
            NG1D = NG1D+NH1D
            VG1D = VG1D+QH1D/RHOG0
            QH1D = 0.; NH1D = 0.
            IF (AH1D.GE.ASMALL.AND.AFAG_3M.EQ.1) THEN
               AG1D = AG1D+AH1D
               AH1D = 0.
            ENDIF
         ENDIF
      ENDIF

      END SUBROUTINE SMALL_DT
!======================================================================
!
!======================================================================
      SUBROUTINE LARGE_DT(DT,TK1D,QV1D,P1D,W1D,RHO,QC1D,QR1D,QI1D,QS1D,&
                 QG1D,QH1D,NC1D,NR1D,NI1D,NS1D,NG1D,NH1D,VI1D,VS1D,    &
                 VG1D,FI1D,FS1D,AI1D,AS1D,AG1D,AH1D,I3M1D,SASPR,GQCTR, &
                 SMLF,GMLF)
!=======================================================================
      IMPLICIT NONE
      INTEGER :: I,WBIN,CBIN,PBIN,TBIN,DBIN
      REAL, INTENT(IN) :: W1D
      REAL :: TK1D,QV1D,P1D,TC1D,DT,GQCTR,QI1D,AI1D,AS1D,AG1D,AH1D,    &
              I3M1D,QNIN,VI1D,VS1D,VG1D,FI1D,FS1D,AFAC,AFAR,AFAI,AFAS, &
              AFAG,AFAH,ADAGR,QC1D,QR1D,QS1D,QG1D,QH1D,NC1D,NR1D,NI1D, &
              NS1D,NG1D,NH1D,RHOI,RHOS,RHOG,RHOIS,iRHOI,iRHOS,iRHOG,   &
              ZETA,ZETA2,ZETA3,iAMI,iAMG,DI0Z1,DI0Z2,DI0Z3,DI0Z4,IPG,  &
              IPH,IPF,SASP1,SASP2,SASP3,SASPR,DNRI,DNIRM,DNIAG,DNSAG,  &
              DNSAC,DNGAC,DNSRM,DNGRM,AMI,BMI,AMG,AMS,RHO,iRHO,AVI,BVI,&
              AVS,BVS,AVG,BVG,AVH,BVH,CAPS,SASR1,SASR2,SASR3,SASR4
      REAL :: QHOci,QHOrg,QNMci,QNMrg,QNCci,QFZci,QFZrg,QMLir,QMLic,   &
              QMLsr,QMLgr,QMLhr,QIMii,QIMcsi,QIMrsi,QIMcgi,QIMrgi,     &
              QBKrc,QINig,QINsg,QCNcr,QCNis,QCNgh,QRMci,QRMcs,QRMcg,   &
              QRMch,QRMri,QRMrs,QRMrg,QRMrh,QCLcr,QCLir,QCLis,QCLig,   &
              QCLih,QCLsr,QCLsg,QCLsh,QCLgr,QCLirg,QCLsrg,QCLgrg,QHwsh,&
              QHwml,QHdrm,QCLss
      REAL :: NHOci,NHOrg,NNMci,NNMrg,NNCci,NFZci,NFZrg,NMLir,NMLic,   &
              NMLsr,NMLgr,NMLhr,NIMii,NIMcsi,NIMrsi,NIMcgi,NIMrgi,     &
              NBKrc,NBKrr,NINig,NINsg,NCNcr,NiCNis,NgCNgh,NhCNgh,NRMci,&
              NRMcs,NRMcg,NRMch,NRMri,NRMrs,NRMrg,NRMrh,NCLcc,NCLcr,   &
              NCLrr,NCLir,NCLis,NCLig,NCLih,NCLsr,NCLss,NCLsg,NCLsh,   &
              NCLgr,NCLirg,NCLsrg,NCLgrg,NHwsh,NsCNis
      REAL :: VIMii,ViCLis,VsCLis,ViCLig,VgCLig,VCLih,ViCNis,VsCNis,   &
              ViINig,VgINig,VFZci,VFZrg,VIMcsi,VMLir,VMLic,VRMci,VCLir,&
              VHOci,VHOrg,VNMci,VNMrg,VNCci,VIMcgi,VIMrsi,VIMrgi,VCLss,&
              VsCLsg,VgCLsg,VCLsh,VsINsg,VgINsg,VMLsr,VRMcs,VCLsr,     &
              VMLgr,VRMcg,VCNgh,VCLgr,VCLirg,VCLsrg,VCLgrg,FHOci,FNMci,&
              FNCci,FRMci,FMLir,FMLic,FCLir,FiCLis,FCLig,FIMcsi,FCLih, &
              FiCNis,FINig,FFZci,FIMii,FIMcgi,FIMrsi,FIMrgi,FsCNis,    &
              FsCLis,FCLsg,FCLsh,FINsg,FMLsr,FRMcs,FCLsr,FCLss
      REAL :: IHOci,INMci,INCci,IRMci,IMLir,IMLic,ICLir,ICLis,ICLig,   &
              IIMcsi,ICLih,ICNis,IINig,IFZci,IIMii,IIMcgi,IIMrsi,      &
              IIMrgi,AHOci,ANMci,ANCci,ARMci,AMLir,AMLic,ACLir,AiCLis, &
              AiCLig,AIMcsi,AiCLih,AiCNis,AsCNis,AiINig,AgINig,AFZci,  &
              AIMii,AIMcgi,AIMrsi,AIMrgi,AMLsr,AsCLis,AsCLsg,AsCLsh,   &
              ARMcs,ACLsr,AgCLig,AsINsg,AgINsg,AMLgr,ARMcg,AgCNgh,     &
              AhCNgh,AgCLsg,ACLgr,ACLirg,ACLsrg,ACLgrg,AHOrg,ANMrg,    &
              AFZrg,AMLhr,ACLss,AhCLih,AhCLsh,ARMch,ARMrh,AHdrm,AHwml
      REAL :: QCSOUR,QCSINK,QRSOUR,QRSINK,QISOUR,QISINK,QSSOUR,QSSINK, &
              QGSOUR,QGSINK,QHSOUR,QHSINK,NCSOUR,NCSINK,NRSOUR,NRSINK, &
              NISOUR,NISINK,NSSOUR,NSSINK,NGSOUR,NGSINK,NHSOUR,NHSINK, &
              FISOUR,FISINK,FSSOUR,FSSINK,VISOUR,VISINK,VSSOUR,VSSINK, &
              VGSOUR,VGSINK,AISOUR,AISINK,ASSOUR,ASSINK,AGSOUR,AGSINK, &
              AHSOUR,AHSINK,IISOUR,IISINK
      REAL :: RHOAJ,CPM,SSRI,XXLS,XXLF,XXLV,QVSI,ESI,ABI,iDT,SSRW,KAP, &
              DV,MUA,SCN,ACTW,ESW,SFCTNW,SFCTNV,EPA,QVSW,TQCI,TQRSG,FF,&
              RFZ,ICED,VOLMC,IJHOF,AREAC,RGIMF,GGIMF,ARIMF0,ARIMF,     &
              COSM2,GEOF2,IJIMF0,IJIMF,GGCNT,NGCNT0,NGCNT1,NGCNT,TCC,  &
              CNTGG0,CNTGG,PSIA,KDIFF,KAPA,CNTKN,CNTFT0,CNTFT1,CNTF1,  &
              CNTF2,CNTFT,IJCNT1,IJCNT2,IJCNT3,DC1,DC2,DC3,INR0,FANGLE,&
              FACTE,DACTE,RATIO,WNRE,CNRE,PNRE,MRATO,MVRC,MVRR,MVDC,   &
              MVDR,MVDI,MVDS,MVDG,MVDH,LMVRC,LMVRR,TMP1,GUC,FUC,FUR,   &
              GUR,RMcsq,RMrsq,RMcgq,RMrgq,RMchq,RMrhq,RMcsa,RMrsa,     &
              RMcga,RMrga,RMcha,RMrha,VTQ0,VTQC,VTQR,VTQI,VTQS,VTQG,   &
              VTQH,VTN0,VTNC,VTNR,VTNI,VTNS,VTNG,VTNH,VTAX,VTAI,VTAS,  &
              VTAG,VTAH,VTFI,VTI3M,VTQIC,VTQSC,VTQGC,VTQHC,VTQRI,VTQRS,&
              VTQRG,VTQRH,VTQIS,VTQIG,VTQIH,VTQSG,VTQSH,VTNIC,VTNSC,   &
              VTNGC,VTNHC,VTNRI,VTNRS,VTNRG,VTNRH,VTNIS,VTNIG,VTNIH,   &
              VTNSG,VTNSH,VTAIS,VTAIG,VTAIH,VTASG,VTASH,EIS,EIG,EIH,   &
              ESG,ESH,EII,ESS,ECI,ECS,ECG,ECH,ERI,ERS,ERG,ERH,EII1,    &
              EII2,EIS1,EIS2,ESS1,ESS2,SMLTQ,GMLTQ,HMLTQ,SMLTA,GMLTA,  &
              HMLTA,VENQS,VENQG,VENQH,VENAS,VENAG,VENAH,STOKE,DSLL,    &
              NGTAL,MVDX,ABW,FSQC,FSNC,FSQR,FSNR,FSQI,FSNI,FSQS,FSNS,  &
              FSQG,FSNG,FSQH,FSNH,FSAI,FSAS,FSAG,FSAH,FSVI,FSFI,DMWDT, &
              DMIDT,HIdqv,HSdqv,HGdqv,HHdqv,HCwqv,HRwqv,HSwqv,HGwqv,   &
              HHwqv,HHdcd,HHdrm,HHdcl,HHdtt,HSwcd,HGwcd,HHwcd,HSwrm,   &
              HGwrm,HHwrm,VDMAX,SBMAX,EVMAX,SUMDEP,SUMSUB,SUMEVP,      &
              SUMCND,MLMAX,ESW0,ESI0,QVSW0,QVSI0,XXLF0,MLWM,MWT,MLWC,  &
              DH9,MCORE,ICOR1,DIT,BEST,VENQI,QTMP0,QTMP1,QTMP2,QTMP3,  &
              QTMP4,QTMP5,QTMP6,QTMP7,QTMP8,QTMP9,ATMP1,ATMP2,ATMP3,   &
              ATMP4,NCLS2,BTMP,BSTMP,BGTMP,BHTMP,H2Z,H4Z,QRMC1,NRMC1,  &
              FRMC1,IRMC1,ARMC1,QRMR1,NRMR1,ARMR1,QFZC1,NFZC1,QFZR1,   &
              NFZR1,AFZR1,QCLI1,NCLI1,QCNI1,NCNI1,ACNI1,FCNI1,QCLG1,   &
              NCLG1,ACLG1,ACLI1,FCLI1,QCLS1,NCLS1,ACLS1,LIM1,LIM2,LIM3,&
              LIM4,LIM5,LIM6,LIMA,LIMB,LIMC,LIMD,LIME,LIMF,LIMG,LIMH,  &
              LIMI,DICC,DICA,DIAC,DIAA,DIC0,DIA0,DIF0,DSS0,DSL0,DIC2,  &
              DIA2,DSLS,DGLS,SMLR,GMLR,SMLF,GMLF,AVWSG,BVWSG,CAPWS,    &
              WSAPR,RHOWS,RHOWG,VENWS,VENWG,SSRW0,SSRI0,LLMI,LLMS,LLMG,&
              LLMH,KDX,RHOIW,RHOSW,RHOGW,RHOIG,RHOIH,RHOSG,RHOSH,RHOHW,&
              KINV,BEST0,C1X2,VTB1,VTA1
      REAL :: LAMC,LAMR,LAMS,LAMG,LAMH,LAMI,GC2,GC3,GC4,GC5,GC6,GC7,   &
              GR2,GR3,GR4,GR5,GR6,GR7,GI1,GI2,GI2G1,GI2G2,GI2H1,GI3,   &
              GI3H1,GIF1,GIF2,GIG1,GIG2,GIH1,GIH2,GIM1,GIM2,GI4,GI5,   &
              GIM2H1,GI2H3,GI2HG1,GIH2G1,GIF3,GI2G3,GI3G1,GIZM2H1,GIG3,&
              GIM2G1,GIZM2G1,GIM3,GIMF1,GIZMF1,GIMG1,GIZMG1,GIMH1,     &
              GIZMH1,GIZM1,GIZ1,GIZ2G1,GIZF1,GIZG1,GIZH1,GIH3,GS2,GS3, &
              GS4,GS5,GSM1,GSM2,GSM3,GG2,GG3,GG4,GG5,GGM1,GGM2,GGM3,   &
              GH2,GH3,GH4,GH5,Z32G,Z32H,Z3BMI
      REAL, PARAMETER :: AVTC = 8.8462E+02, BVTC = 9.7593E+07, CVTC = -3.4249E+11
      REAL, PARAMETER :: AVTR = 2.1454E+00, BVTR = -2.2812E-04, CVTR = 2.9676E-09
      REAL, PARAMETER :: CQC1 = 2.0901E+01, CQC2 = 9.9111E-01, CQC3 = 4.4182E+00
      REAL, PARAMETER :: CNC1 = 1.8276E+01, CNC2 = 1.0015E+00, CNC3 = 1.9838E+00
      REAL, PARAMETER :: CQR1 = 1.5943E+01, CQR2 = 1.1898E+00, CQR3 = 4.0073E+00
      REAL, PARAMETER :: CNR1 = 9.4791E+00, CNR2 = 9.7607E-01, CNR3 = 1.0858E+00
      REAL, PARAMETER :: AQ1 = 6.6793E+0, BQ1 = 1.0090E+0, CQ1 = 1.4095E+0
      REAL, PARAMETER :: AQ2 = 9.9912E+0, BQ2 =-4.7678E-1, CQ2 =-3.1388E-2
      REAL, PARAMETER :: AN3 = -4.3561E0, BN3 = 1.9934E0,  CN3 = 1.6465E-2
      REAL, PARAMETER :: AN4 = -4.0731E1, BN4 = 5.3720E5,  CN4 = -2.0139E-5, AQ4 = -2.1370E1, BQ4 = 1.9899E9
      REAL, PARAMETER :: AN5 = 1.5519E1,  BN5 = 3.1491E-0, CN5 = 4.3989E-1,  AQ5 = 2.0090E1,  BQ5 = 2.9626E0, CQ5 = 3.2358E0
      REAL, PARAMETER :: AN6 = -1.8239E1, BN6 = 2.2956E0,  CN6 = -2.3261E-4
      REAL, PARAMETER :: AN7 = -1.7431E2, BN7 = 2.6031E5,  CN7 = -9.3613E7
      REAL, PARAMETER :: AN8 = -1.6185E2, BN8 = 2.2786E5,  CN8 = -7.6988E7, AQ8 = -2.3531E1, BQ8 = 9.8271E-1, CQ8 = -1.3202E-1
      REAL, PARAMETER :: DSHED = 1.E-3, XISP = 0.25, AX = 6., BX = 10.
      REAL :: ZC1,ZC2,ZC3,ZC4,ZP1,ZP2,ZP3,ZP4,ZZ1,ZZ2,ZZ3,ZZ4,ZZ5,ZZ6, &
              ZZ7,ZZ8
      DATA ZC1,ZC2,ZC3,ZC4/0.69509913,-0.46685819,0.30490087,1.62148100/
      DATA ZP1,ZP2,ZP3,ZP4/0.36793126,1.82782890,0.63206874,-1.00164090/
      DATA ZZ1,ZZ2,ZZ3/0.91767336,-0.37871035,1.03995157/

      IF (INSPEC.EQ.1) THEN                                             ! SOOT
         DC1 = -0.5411955; DC2 = 1.879918; DC3 = 1.607947               ! CHEN Table3
         INR0 = 0.4E-7; FANGLE = 33.2; FACTE = 13.8E-20; DACTE = -20.E-20
      ELSEIF (INSPEC.EQ.2) THEN                                         ! SAHARAN DUST
         DC1 = -0.3353619; DC2 = 1.990979; DC3 = 2.175539               ! CHEN Table3
         INR0 = 1.75E-7; FANGLE = 30.98; FACTE = 15.7E-20; DACTE = 3.35E-20
      ELSEIF (INSPEC.EQ.3) THEN                                         ! ASIAN DUST
         DC1 = -0.3598818; DC2 = 1.982032; DC3 = 2.025390               ! CHEN Table3
         INR0 = 2.E-7; FANGLE = 30.98; FACTE = 15.7E-20; DACTE = 1.82E-20
      ELSEIF (INSPEC.EQ.4) THEN                                         ! BACTERIA
         DC1 = -0.3541568; DC2 = 1.983928; DC3 = 2.019456               ! CHEN Table3
         INR0 = 5.E-7; FANGLE = 14.82; FACTE = 17.6E-20; DACTE = 1.82E-20 ! DACTE NOT SURE
      ENDIF

      MVDC   = 0.; MVDR   = 0.; MVDI   = 0.; MVDS   = 0.; MVDG   = 0.
      MVDH   = 0.; MVRC   = 0.; MVRR   = 0.; SMLF   = 0.; GMLF   = 0.
      FF     = 0.; DMWDT  = 0.; DMIDT  = 0.; NsCNis = 0.
      HIdqv  = 0.; HSdqv  = 0.; HGdqv  = 0.; HHdqv  = 0.; HCwqv  = 0.
      HRwqv  = 0.; HSwqv  = 0.; HGwqv  = 0.; HHwqv  = 0.; HHdtt  = 0.
      QHOci  = 0.; QHOrg  = 0.; QNMci  = 0.; QNMrg  = 0.; QNCci  = 0.
      QFZci  = 0.; QFZrg  = 0.; QMLir  = 0.; QMLic  = 0.; QMLsr  = 0.
      QMLgr  = 0.; QMLhr  = 0.; QIMii  = 0.; QIMcsi = 0.; QIMrsi = 0.
      QIMcgi = 0.; QIMrgi = 0.; QBKrc  = 0.; QINig  = 0.; QINsg  = 0.
      QCNcr  = 0.; QCNis  = 0.; QCNgh  = 0.; QRMci  = 0.; QRMcs  = 0.
      QRMcg  = 0.; QRMch  = 0.; QRMri  = 0.; QRMrs  = 0.; QRMrg  = 0.
      QRMrh  = 0.; QCLcr  = 0.; QCLir  = 0.; QCLis  = 0.; QCLig  = 0.
      QCLih  = 0.; QCLsr  = 0.; QCLsg  = 0.; QCLsh  = 0.; QCLgr  = 0.
      QCLirg = 0.; QCLsrg = 0.; QCLgrg = 0.; QHwsh  = 0.; QHdrm  = 0.
      QHwml  = 0.; NHOci  = 0.; NHOrg  = 0.; NNMci  = 0.; NNMrg  = 0.
      NNCci  = 0.; NFZci  = 0.; NFZrg  = 0.; NMLir  = 0.; NMLic  = 0.
      NMLsr  = 0.; NMLgr  = 0.; NMLhr  = 0.; NIMii  = 0.; NIMcsi = 0.
      NIMrsi = 0.; NIMcgi = 0.; NIMrgi = 0.; NBKrc  = 0.; NBKrr  = 0.
      NINig  = 0.; NINsg  = 0.; NCNcr  = 0.; NiCNis = 0.; NgCNgh = 0.
      NRMci  = 0.; NRMcs  = 0.; NRMcg  = 0.; NRMch  = 0.; NRMri  = 0.
      NRMrs  = 0.; NRMrg  = 0.; NRMrh  = 0.; NCLcc  = 0.; NCLcr  = 0.
      NCLrr  = 0.; NCLir  = 0.; NCLis  = 0.; NCLig  = 0.; NCLih  = 0.
      NCLsr  = 0.; NCLss  = 0.; NCLsg  = 0.; NCLsh  = 0.; NCLgr  = 0.
      NCLirg = 0.; NCLsrg = 0.; NCLgrg = 0.; NHwsh  = 0.; NhCNgh = 0.
      VCLss  = 0.; ViCNis = 0.; VIMii  = 0.; ViCLis = 0.; VsCLis = 0.
      ViCLig = 0.; VgCLig = 0.; VCLih  = 0.; VsCNis = 0.; ViINig = 0.
      VgINig = 0.; VFZci  = 0.; VFZrg  = 0.; VMLir  = 0.; VMLic  = 0.
      VRMci  = 0.; VCLir  = 0.; VHOci  = 0.; VHOrg  = 0.; VNMci  = 0.
      VNMrg  = 0.; VNCci  = 0.; VIMcsi = 0.; VIMcgi = 0.; VIMrsi = 0.
      VIMrgi = 0.; VCLsr  = 0.; VsCLsg = 0.; VgCLsg = 0.; VCLsh  = 0.
      VsINsg = 0.; VgINsg = 0.; VMLsr  = 0.; VRMcs  = 0.; VMLgr  = 0.
      VRMcg  = 0.; VCNgh  = 0.; VCLgr  = 0.; VCLirg = 0.; VCLsrg = 0.
      VCLgrg = 0.; FIMcsi = 0.; FIMcgi = 0.; FIMrsi = 0.; FIMrgi = 0.
      FIMii  = 0.; FHOci  = 0.; FNMci  = 0.; FNCci  = 0.; FFZci  = 0.
      FMLir  = 0.; FCLir  = 0.; FiCLis = 0.; FCLig  = 0.; FCLih  = 0.
      FRMci  = 0.; FiCNis = 0.; FINig  = 0.; FMLic  = 0.; FsCNis = 0.
      FsCLis = 0.; FCLsg  = 0.; FCLsh  = 0.; FINsg  = 0.; FMLsr  = 0.
      FRMcs  = 0.; FCLsr  = 0.; FCLss  = 0.; IIMcsi = 0.; IIMcgi = 0.
      IIMrsi = 0.; IIMrgi = 0.; IIMii  = 0.; IHOci  = 0.; INMci  = 0.
      INCci  = 0.; IFZci  = 0.; IMLir  = 0.; IMLic  = 0.; ICLir  = 0.
      ICLis  = 0.; ICLig  = 0.; ICLih  = 0.; IRMci  = 0.; ICNis  = 0.
      IINig  = 0.; AIMcsi = 0.; AIMcgi = 0.; AIMrsi = 0.; AIMrgi = 0.
      AIMii  = 0.; AHOci  = 0.; ANMci  = 0.; ANCci  = 0.; AFZci  = 0.
      AMLic  = 0.; AMLir  = 0.; ACLir  = 0.; AiCLis = 0.; AiCLig = 0.
      AiCLih = 0.; ARMci  = 0.; AiCNis = 0.; AiINig = 0.; AgINig = 0.
      AsINsg = 0.; AgINsg = 0.; AsCLis = 0.; AsCLsg = 0.; AsCLsh = 0.
      ARMcs  = 0.; ACLsr  = 0.; AMLsr  = 0.; AgCLig = 0.; AgCLsg = 0.
      AMLgr  = 0.; ARMcg  = 0.; AgCNgh = 0.; ACLgr  = 0.; ACLirg = 0.
      ACLsrg = 0.; ACLgrg = 0.; AHOrg  = 0.; ANMrg  = 0.; AFZrg  = 0.
      AMLhr  = 0.; AhCLih = 0.; AhCLsh = 0.; ARMch  = 0.; ARMrh  = 0.
      AHdrm  = 0.; AHwml  = 0.; AsCNis = 0.; AhCNgh = 0.; ACLss  = 0.

      iDT    = 1./DT
      iRHO   = 1./RHO
      TC1D   = TK1D-TK0C
      IF (TK1D.GE.265..AND.TK1D.LT.268.) FF = (TC1D+8.)*THRD
      IF (TK1D.GE.268..AND.TK1D.LE.270.) FF = (-3.-TC1D)*0.5
      RHOAJ  = (RHOSU/RHO)**0.54
      CPM    = CP*(1.+0.887*QV1D)
      ESW    = POLYSVP(TK1D,0)
      ESI    = POLYSVP(TK1D,1)
      QVSW   = 0.622*ESW/(P1D-ESW)
      QVSI   = 0.622*ESI/(P1D-ESI)
      SSRW   = QV1D/QVSW-1.
      SSRI   = QV1D/QVSI-1.
      ESW0   = POLYSVP(TK0C,0)
      QVSW0  = 0.622*ESW0/(P1D-ESW0)
      ESI0   = POLYSVP(TK0C,1)
      QVSI0  = 0.622*ESI0/(P1D-ESI0)
      XXLV   = 3.1484E6-2370.*TK1D
      XXLS   = 3.15E6-2370.*TK1D+0.3337E6
      XXLF   = 2836310.8-(3.1484E6-2370.*TK1D)
      XXLF0  = 2836310.8-(3.1484E6-2370.*TK0C)
      SSRW0  = QV1D/QVSW0-1.
      SSRI0  = QV1D/QVSI0-1.
      DV     = 2.11E-5*(TK1D/TK0C)**1.94*(101325./P1D)
      KAP    = 2.3971E-2+0.0078E-2*TC1D
      MUA    = 1.72E-5*(393./(TK1D+120.))*(TK1D/TK0C)**1.5
      SCN    = (MUA/(RHO*DV))**THRD
      ICED   = 916.7-0.175*TC1D-5.E-4*TC1D**2.                          ! ICE DENSITY
      ABW    = TK1D*RV/ESW/DV+XXLV*(XXLV/TK1D/RV-1.)/TK1D/KAP
      ABI    = TK1D*RV/ESI/DV+XXLS*(XXLS/TK1D/RV-1.)/TK1D/KAP
      SFCTNV = (0.095*TC1D+104.6)*1.E-3                                 ! SURFACE TENSION OVER ICE/VAPOR
      SFCTNW = (28.+0.25*TC1D)*1.E-3
      COSM2  = LOG(1.-COS(FANGLE*PI/180.))
      ACTW   = LOG(1.*ESW/ESI)                                          ! WATER ACTIVITY NEED TO BE MODIFIED
      EIH    = MIN(1.,0.01*EXP(0.1*TC1D))                               ! FERRIER ET AL., 1995
      ESH    = MIN(1.,0.01*EXP(0.1*TC1D))                               ! FERRIER ET AL., 1995

      IF (QC1D.GE.QSMALL) THEN
         CALL SOLVE_AFAC(QC1D,NC1D,LAMC,MVDC,AFAC)
         GC2   = EXP(GAMLN(AFAC+2.)-GAMLN(AFAC+1.)-LOG(LAMC))
         GC3   = EXP(GAMLN(AFAC+3.)-GAMLN(AFAC+1.)-2.*LOG(LAMC))
         GC4   = EXP(GAMLN(AFAC+4.)-GAMLN(AFAC+1.)-3.*LOG(LAMC))
         GC5   = EXP(GAMLN(AFAC+5.)-GAMLN(AFAC+1.)-4.*LOG(LAMC))
         GC6   = EXP(GAMLN(AFAC+6.)-GAMLN(AFAC+1.)-5.*LOG(LAMC))
         GC7   = EXP(GAMLN(AFAC+7.)-GAMLN(AFAC+1.)-6.*LOG(LAMC))
         MVRC  = MIN(MAX((QC1D/NC1D/C4PI3W)**THRD,RCMIN),RCMAX)
         LMVRC = LOG(MVRC)
         GUC   = EXP(EXP(AFU+BFU*LMVRC**3.+CFU*SQRT(RHO)**3.))
         IF (LIQ_VTC.EQ.0) THEN
            FSQC = EXP(GAMLN(BVC0+BMW+AFAC+1.)-GAMLN(BMW+AFAC+1.)-     &
                   BVC0*LOG(LAMC))
            FSNC = EXP(GAMLN(BVC0+AFAC+1.)-GAMLN(AFAC+1.)-BVC0*        &
                   LOG(LAMC))
            VTQC = RHOAJ*FSQC*AVC0
            VTNC = RHOAJ*FSNC*AVC0
         ELSEIF (LIQ_VTC.EQ.1) THEN
            VTQC = EXP(CQC1+CQC2*LOG(NC1D)+CQC3*LMVRC)*GUC/QC1D
            VTNC = EXP(CNC1+CNC2*LOG(NC1D)+CNC3*LMVRC)*GUC/NC1D
         ELSEIF (LIQ_VTC.EQ.2) THEN
            FUC  = EXP(EXP(AFU+BFU*(LMVRC)**3.+CFU*SQRT(RHO)**3.))
            VTQC = MVRC*(AVTC+BVTC*MVRC+CVTC*MVRC**2.)*FUC
            VTNC = MVRC*(AVTC+BVTC*MVRC+CVTC*MVRC**2.)*FUC
         ENDIF
         VTQC  = MIN(VTQC,VTCMAX)
         VTNC  = MIN(VTNC,VTCMAX)
         QRMC1 = PI*PI*RHOW*NC1D/24.
         NRMC1 = PI*NC1D/4.
         FRMC1 = PI*RHOW*NC1D
         IRMC1 = PI*RHOW*NC1D/4.
         ARMC1 = PI*NC1D/4.
         QFZC1 = PI*PI*RHOW*NC1D*GC7/36.
         NFZC1 = PI*NC1D*GC4/6.
         HCwqv = SSRW0*EXP(AQ1+BQ1*LOG(NC1D)+CQ1*LMVRC)/ABW
      ENDIF
      IF (QR1D.GE.QSMALL) THEN
         CALL SOLVE_AFAR(QR1D,NR1D,LAMR,MVDR,AFAR)
         GR2   = EXP(GAMLN(AFAR+2.)-GAMLN(AFAR+1.)-LOG(LAMR))
         GR3   = EXP(GAMLN(AFAR+3.)-GAMLN(AFAR+1.)-2.*LOG(LAMR))
         GR4   = EXP(GAMLN(AFAR+4.)-GAMLN(AFAR+1.)-3.*LOG(LAMR))
         GR5   = EXP(GAMLN(AFAR+5.)-GAMLN(AFAR+1.)-4.*LOG(LAMR))
         GR6   = EXP(GAMLN(AFAR+6.)-GAMLN(AFAR+1.)-5.*LOG(LAMR))
         GR7   = EXP(GAMLN(AFAR+6.)-GAMLN(AFAR+1.)-6.*LOG(LAMR))
         MVRR  = MIN(MAX((QR1D/NR1D/C4PI3W)**THRD,RRMIN),RRMAX)
         LMVRR = LOG(MVRR)
         GUR   = EXP(EXP(AFU+BFU*LMVRR**3.+CFU*SQRT(RHO)**3.))
         IF (LIQ_VTR.EQ.0) THEN
            FSQR = EXP(GAMLN(BVR0+BMW+AFAR+1.)-GAMLN(BMW+AFAR+1.)-     &
                   BVR0*LOG(LAMR))
            FSNR = EXP(GAMLN(BVR0+AFAR+1.)-GAMLN(AFAR+1.)-BVR0*        &
                   LOG(LAMR))
            VTQR = RHOAJ*FSQR*AVR0
            VTNR = RHOAJ*FSNR*AVR0
         ELSEIF (LIQ_VTR.EQ.1) THEN
            VTQR = EXP(CQR1+CQR2*LOG(NR1D)+CQR3*LMVRR)*GUR/QR1D
            VTNR = EXP(CNR1+CNR2*LOG(NR1D)+CNR3*LMVRR)*GUR/NR1D
         ELSEIF (LIQ_VTR.EQ.2) THEN
            FUR  = EXP(EXP(AFU+BFU*(LOG(MVRR))**3.+CFU*SQRT(RHO)**3.))
            VTQR = EXP(AVTR+BVTR/MVRR+CVTR/(MVRR**2.))*FUR
            VTNR = EXP(AVTR+BVTR/MVRR+CVTR/(MVRR**2.))*FUR
         ENDIF
         VTQR  = MIN(VTQR,VTRMAX)
         VTNR  = MIN(VTNR,VTRMAX)
         QRMR1 = PI*PI*RHOW*NR1D/24.
         NRMR1 = PI*NR1D/4.
         ARMR1 = PI*NR1D/4.
         QFZR1 = PI*PI*RHOW*NR1D*GR7/36.
         NFZR1 = PI*NR1D*GR4/6.
         AFZR1 = PI*NR1D*GR6/6.
         HRwqv = SSRW0*NR1D*EXP(AQ2+(BQ2+CQ2*LMVRR)*LMVRR**2.)/ABW
      ENDIF
      IF (QI1D.GE.QSMALL) THEN
         CALL SOLVE_AFAI(TK1D,RHO,QI1D,NI1D,VI1D,FI1D,AI1D,I3M1D,ADAGR,&
              ZETA,LAMI,AFAI,MVDI,RHOI,AMI,BMI,AVI,BVI,BEST)
         IPG    = 3./(ADAGR+2.)
         IPH    = 3.*ADAGR/(ADAGR+2.)
         IPF    = 3.*(ADAGR+1.)/(ADAGR+2.)
         ZETA2  = 2.*(ADAGR-1.)/(ADAGR+2.)
         ZETA3  = 3.*(ADAGR-1.)/(ADAGR+2.)
         DI0Z1  = DI0**ZETA
         DI0Z2  = DI0**ZETA2
         DI0Z3  = DI0**ZETA3
         DI0Z4  = DI0**(4.*ZETA)
         iRHOI  = 1./RHOI
         iAMI   = 1./AMI
         QCLI1  = PI*AMI*NI1D/4.
         QCNI1  = PI*XISP*AMI*NI1D/6.
         NCLI1  = PI*NI1D/4.
         NCNI1  = PI*XISP*NI1D/6.
         ACLI1  = PI*NI1D/4.
         ACNI1  = PI*XISP*NI1D/6.
         FCLI1  = PI*NI1D/4./DI0Z3
         FCNI1  = AMI*XISP*iRHOI*NI1D
         Z32G   = ZETA3+2.*IPG
         Z32H   = ZETA3+2.*IPH
         Z3BMI  = ZETA3+BMI
         LLMI   = LOG(LAMI)
         GI1    = GAMLN(AFAI+1.)
         GI2    = EXP(GAMLN(AFAI+2.)-GI1-LLMI)
         GI3    = EXP(GAMLN(AFAI+3.)-GI1-2.*LLMI)
         GI4    = EXP(GAMLN(AFAI+4.)-GI1-3.*LLMI)
         GI5    = EXP(GAMLN(AFAI+5.)-GI1-4.*LLMI)
         GIM1   = EXP(GAMLN(AFAI+BMI+1.)-GI1-BMI*LLMI)
         GIM2   = EXP(GAMLN(AFAI+BMI+2.)-GI1-(BMI+1.)*LLMI)
         GIM3   = EXP(GAMLN(AFAI+BMI+3.)-GI1-(BMI+2.)*LLMI)
         GIF1   = EXP(GAMLN(AFAI+IPF+1.)-GI1-IPF*LLMI)
         GIF2   = EXP(GAMLN(AFAI+IPF+2.)-GI1-(IPF+1.)*LLMI)
         GIF3   = EXP(GAMLN(AFAI+IPF+3.)-GI1-(IPF+2.)*LLMI)
         GIG1   = EXP(GAMLN(AFAI+IPG+1.)-GI1-IPG*LLMI)
         GIG2   = EXP(GAMLN(AFAI+IPG+2.)-GI1-(IPG+1.)*LLMI)
         GIG3   = EXP(GAMLN(AFAI+IPG+3.)-GI1-(IPG+2.)*LLMI)
         GIH1   = EXP(GAMLN(AFAI+IPH+1.)-GI1-IPH*LLMI)
         GIH2   = EXP(GAMLN(AFAI+IPH+2.)-GI1-(IPH+1.)*LLMI)
         GIH3   = EXP(GAMLN(AFAI+IPH+3.)-GI1-(IPH+2.)*LLMI)
         GIZ1   = EXP(GAMLN(AFAI+ZETA3+1.)-GI1-ZETA3*LLMI)
         GI2G1  = EXP(GAMLN(AFAI+2.*IPG+1.)-GI1-2.*IPG*LLMI)
         GI2G2  = EXP(GAMLN(AFAI+2.*IPG+2.)-GI1-(2.*IPG+1.)*LLMI)
         GI2G3  = EXP(GAMLN(AFAI+2.*IPG+3.)-GI1-(2.*IPG+2.)*LLMI)
         GI3G1  = EXP(GAMLN(AFAI+3.*IPG+1.)-GI1-3.*IPG*LLMI)
         GI2H1  = EXP(GAMLN(AFAI+2.*IPH+1.)-GI1-2.*IPH*LLMI)
         GI2H3  = EXP(GAMLN(AFAI+2.*IPH+3.)-GI1-(2.*IPH+2.)*LLMI)
         GI3H1  = EXP(GAMLN(AFAI+3.*IPH+1.)-GI1-3.*IPH*LLMI)
         GIMF1  = EXP(GAMLN(AFAI+BMI+IPF+1.)-GI1-(BMI+IPF)*LLMI)
         GIMG1  = EXP(GAMLN(AFAI+BMI+IPG+1.)-GI1-(BMI+IPG)*LLMI)
         GIMH1  = EXP(GAMLN(AFAI+BMI+IPH+1.)-GI1-(BMI+IPH)*LLMI)
         GIZM1  = EXP(GAMLN(AFAI+ZETA3+BMI+1.)-GI1-(ZETA3+BMI)*LLMI)
         GIZF1  = EXP(GAMLN(AFAI+ZETA3+IPF+1.)-GI1-(ZETA3+IPF)*LLMI)
         GIZG1  = EXP(GAMLN(AFAI+ZETA3+IPG+1.)-GI1-(ZETA3+IPG)*LLMI)
         GIZH1  = EXP(GAMLN(AFAI+ZETA3+IPH+1.)-GI1-(ZETA3+IPH)*LLMI)
         GI2HG1 = EXP(GAMLN(AFAI+2.*IPH+IPG+1.)-GI1-(2.*IPH+IPG)*LLMI)
         GIH2G1 = EXP(GAMLN(AFAI+IPH+2.*IPG+1.)-GI1-(IPH+2.*IPG)*LLMI)
         GIM2G1 = EXP(GAMLN(AFAI+BMI+2.*IPG+1.)-GI1-(BMI+2.*IPG)*LLMI)
         GIM2H1 = EXP(GAMLN(AFAI+BMI+2.*IPH+1.)-GI1-(BMI+2.*IPH)*LLMI)
         GIZ2G1 = EXP(GAMLN(AFAI+Z32G+1.)-GI1-Z32G*LLMI)
         GIZMF1 = EXP(GAMLN(AFAI+Z3BMI+IPF+1.)-GI1-(Z3BMI+IPF)*LLMI)
         GIZMG1 = EXP(GAMLN(AFAI+Z3BMI+IPG+1.)-GI1-(Z3BMI+IPG)*LLMI)
         GIZMH1 = EXP(GAMLN(AFAI+Z3BMI+IPH+1.)-GI1-(Z3BMI+IPH)*LLMI)
         GIZM2G1 = EXP(GAMLN(AFAI+Z32G+BMI+1.)-GI1-(Z32G+BMI)*LLMI)
         GIZM2H1 = EXP(GAMLN(AFAI+Z32H+BMI+1.)-GI1-(Z32H+BMI)*LLMI)
         FSQI  = EXP(GAMLN(BVI+BMI+AFAI+1.)-GAMLN(BMI+AFAI+1.)-BVI*LLMI)
         FSNI  = EXP(GAMLN(BVI+AFAI+1.)-GI1-BVI*LLMI)
         VTQI  = MIN(RHOAJ*FSQI*AVI,VTIMAX)
         VTNI  = MIN(RHOAJ*FSNI*AVI,VTIMAX)
         IF (AI1D.GE.ASMALL) THEN
            FSAI = EXP(GAMLN(BVI+AFAI+3.)-GAMLN(AFAI+3.)-BVI*LLMI)
            VTAI = MIN(RHOAJ*FSAI*AVI,VTIMAX)
         ENDIF
         IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
            FSFI  = EXP(GAMLN(BVI+ZETA3+AFAI+4.)-GAMLN(ZETA3+AFAI+4.)- &
                    BVI*LLMI)
            FSVI  = EXP(GAMLN(BVI+AFAI+4.)-GAMLN(AFAI+4.)-BVI*LLMI)
            VTFI  = MIN(RHOAJ*FSFI*AVI,VTIMAX)
            VTI3M = MIN(RHOAJ*FSVI*AVI,VTIMAX)
         ENDIF
      IF (ICE_VENT.EQ.3) THEN
         IF ((ADAGR-1.).GE.SLIMIT) THEN
            BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
            H2Z   = ZC2*ZETA
            H4Z   = ZC4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            QTMP2 = LLMI*(H2Z+BVI/2.+IPH/2.+1.)
            QTMP3 = LLMI*(H4Z+BVI/2.+IPH/2.+1.)
            QTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPH/2.+AFAI+2.)-GI1-QTMP2)
            QTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPH/2.+AFAI+2.)-GI1-QTMP3)
            QTMP6 = LLMI*(H2Z+BVI+IPH+1.)
            QTMP7 = LLMI*(H4Z+BVI+IPH+1.)
            QTMP8 = EXP(GAMLN(H2Z+BVI+IPH+AFAI+2.)-GI1-QTMP6)
            QTMP9 = EXP(GAMLN(H4Z+BVI+IPH+AFAI+2.)-GI1-QTMP7)
            VENQI = ZC1*QTMP0/DI0**H2Z+ZC3*QTMP1/DI0**H4Z+VENC1*ZC1*   &
                    BTMP*QTMP4/DI0**(H2Z+ZETA)+VENC1*ZC3*BTMP*QTMP5/   &
                    DI0**(H4Z+ZETA)+VENC2*ZC1*BTMP**2.*QTMP8/DI0**(H2Z+&
                    ZETA2)+VENC2*ZC3*BTMP**2.*QTMP9/DI0**(H4Z+ZETA2)
         ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
            BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
            H2Z   = ZP2*ZETA
            H4Z   = ZP4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            QTMP2 = LLMI*(H2Z+BVI/2.+IPG/2.+1.)
            QTMP3 = LLMI*(H4Z+BVI/2.+IPG/2.+1.)
            QTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPG/2.+AFAI+2.)-GI1-QTMP2)
            QTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPG/2.+AFAI+2.)-GI1-QTMP3)
            QTMP6 = LLMI*(H2Z+BVI+IPG+1.)
            QTMP7 = LLMI*(H4Z+BVI+IPG+1.)
            QTMP8 = EXP(GAMLN(H2Z+BVI+IPG+AFAI+2.)-GI1-QTMP6)
            QTMP9 = EXP(GAMLN(H4Z+BVI+IPG+AFAI+2.)-GI1-QTMP7)
            VENQI = ZP1*QTMP0/DI0**H2Z+ZP3*QTMP1/DI0**H4Z+VENP1*ZP1*   &
                    BTMP*QTMP4/DI0**(H2Z-ZETA/2.)+VENP1*ZP3*BTMP*QTMP5/&
                    DI0**(H4Z-ZETA/2.)+VENP2*ZP1*BTMP**2.*QTMP8/DI0**  &
                    (H2Z-ZETA)+VENP2*ZP3*BTMP**2.*QTMP9/DI0**(H4Z-ZETA)
         ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
            BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
            QTMP0 = EXP(GAMLN(AFAI+2.)-GI1-LOG(LAMI))
            QTMP1 = LLMI*(1.5+BVI/2.)
            QTMP2 = EXP(GAMLN(BVI/2.+AFAI+2.5)-GI1-QTMP1)
            VENQI = AVSG*QTMP0+BVSG*BTMP*QTMP2
         ENDIF
      ELSEIF (ICE_VENT.EQ.1.OR.ICE_VENT.EQ.2) THEN
         IF ((ADAGR-1.).GE.SLIMIT) THEN
            H2Z   = ZC2*ZETA
            H4Z   = ZC4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            IF (BEST.LE.1.) THEN
               BTMP  = SCN**2.*(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(H2Z+BVI+IPH+1.)
               QTMP3 = LLMI*(H4Z+BVI+IPH+1.)
               QTMP4 = EXP(GAMLN(H2Z+BVI+IPH+AFAI+2.)-GI1-QTMP2)
               QTMP5 = EXP(GAMLN(H4Z+BVI+IPH+AFAI+2.)-GI1-QTMP3)
               VENQI = AVIS*ZC1*QTMP0/DI0**H2Z+AVIS*ZC3*QTMP1/DI0**H4Z+&
                       BVIS*ZC1*BTMP*QTMP4/DI0**(H2Z+ZETA2)+BVIS*ZC3*  &
                       BTMP*QTMP5/DI0**(H4Z+ZETA2)
            ELSEIF (BEST.GT.1.) THEN
               BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(H2Z+BVI/2.+IPH/2.+1.)
               QTMP3 = LLMI*(H4Z+BVI/2.+IPH/2.+1.)
               QTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPH/2.+AFAI+2.)-GI1-QTMP2)
               QTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPH/2.+AFAI+2.)-GI1-QTMP3)
               VENQI = AVSG*ZC1*QTMP0/DI0**H2Z+AVSG*ZC3*QTMP1/DI0**H4Z+&
                       BVSG*ZC1*BTMP*QTMP4/DI0**(H2Z+ZETA)+BVSG*ZC3*   &
                       BTMP*QTMP5/DI0**(H4Z+ZETA)
            ENDIF
         ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
            H2Z   = ZP2*ZETA
            H4Z   = ZP4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            IF (BEST.LE.1.) THEN
               BTMP  = SCN**2.*(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(H2Z+BVI+IPG+1.)
               QTMP3 = LLMI*(H4Z+BVI+IPG+1.)
               QTMP4 = EXP(GAMLN(H2Z+BVI+IPG+AFAI+2.)-GI1-QTMP2)
               QTMP5 = EXP(GAMLN(H4Z+BVI+IPG+AFAI+2.)-GI1-QTMP3)
               VENQI = AVIS*ZP1*QTMP0/DI0**H2Z+AVIS*ZP3*QTMP1/DI0**H4Z+&
                       BVIS*ZP1*BTMP*QTMP4/DI0**(H2Z-ZETA)+BVIS*ZP3*   &
                       BTMP*QTMP5/DI0**(H4Z-ZETA)
            ELSEIF (BEST.GT.1.) THEN
               BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
               QTMP2 = LLMI*(H2Z+BVI/2.+IPG/2.+1.)
               QTMP3 = LLMI*(H4Z+BVI/2.+IPG/2.+1.)
               QTMP4 = EXP(GAMLN(H2Z+BVI/2.+IPG/2.+AFAI+2.)-GI1-QTMP2)
               QTMP5 = EXP(GAMLN(H4Z+BVI/2.+IPG/2.+AFAI+2.)-GI1-QTMP3)
               VENQI = AVSG*ZP1*QTMP0/DI0**H2Z+AVSG*ZP3*QTMP1/DI0**H4Z+&
                       BVSG*ZP1*BTMP*QTMP4/DI0**(H2Z-ZETA/2.)+BVSG*ZP3*&
                       BTMP*QTMP5/DI0**(H4Z-ZETA/2.)
            ENDIF
         ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
            QTMP0 = EXP(GAMLN(AFAI+2.)-GI1-LOG(LAMI))
            IF (BEST.LE.1.) THEN
               BTMP  = SCN**2.*(AVI*RHOAJ/MUA)
               QTMP1 = LLMI*(BVI+2.)
               QTMP2 = EXP(GAMLN(BVI+AFAI+3.)-GI1-QTMP1)
               VENQI = AVIS*QTMP0+BVIS*BTMP*QTMP2
            ELSEIF (BEST.GT.1.) THEN
               BTMP  = SCN*SQRT(AVI*RHOAJ/MUA)
               QTMP1 = LLMI*(1.5+BVI/2.)
               QTMP2 = EXP(GAMLN(BVI/2.+AFAI+2.5)-GI1-QTMP1)
               VENQI = AVSG*QTMP0+BVSG*BTMP*QTMP2
            ENDIF
         ENDIF
      ELSEIF (ICE_VENT.EQ.0) THEN
         IF ((ADAGR-1.).GE.SLIMIT) THEN
            H2Z   = ZC2*ZETA
            H4Z   = ZC4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            VENQI = ZC1*QTMP0/DI0**H2Z+ZC3*QTMP1/DI0**H4Z
         ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
            H2Z   = ZP2*ZETA
            H4Z   = ZP4*ZETA
            QTMP0 = EXP(GAMLN(H2Z+AFAI+2.)-GI1-LLMI*(H2Z+1.))
            QTMP1 = EXP(GAMLN(H4Z+AFAI+2.)-GI1-LLMI*(H4Z+1.))
            VENQI = ZP1*QTMP0/DI0**H2Z+ZP3*QTMP1/DI0**H4Z
         ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
            VENQI = EXP(GAMLN(AFAI+2.)-GI1-LOG(LAMI))
         ENDIF
      ENDIF                                                             ! ICE_VENT
      HIdqv = 2.*PI*NI1D*VENQI*XXLS*SSRI0/ABI
      ENDIF                                                             ! QI QSMALL
      IF (QS1D.GE.QSMALL) THEN
         CALL SOLVE_AFAS(TK1D,RHO,QS1D,NS1D,VS1D,FS1D,AS1D,AFAS,LAMS,  &
              MVDS,RHOS,SASPR,AMS,AVS,BVS,SMLF)
         LLMS  = LOG(LAMS)
         GS2   = EXP(GAMLN(AFAS+2.)-GAMLN(AFAS+1.)-LLMS)
         GS3   = EXP(GAMLN(AFAS+3.)-GAMLN(AFAS+1.)-2.*LLMS)
         GS4   = EXP(GAMLN(AFAS+4.)-GAMLN(AFAS+1.)-3.*LLMS)
         GS5   = EXP(GAMLN(AFAS+5.)-GAMLN(AFAS+1.)-4.*LLMS)
         GSM1  = EXP(GAMLN(AFAS+BMS+1.)-GAMLN(AFAS+1.)-BMS*LLMS)
         GSM2  = EXP(GAMLN(AFAS+BMS+2.)-GAMLN(AFAS+1.)-(BMS+1.)*LLMS)
         GSM3  = EXP(GAMLN(AFAS+BMS+3.)-GAMLN(AFAS+1.)-(BMS+2.)*LLMS)
         iRHOS = 1./RHOS
         FSQS  = EXP(GAMLN(BVS+BMS+AFAS+1.)-GAMLN(BMS+AFAS+1.)-BVS*LLMS)
         FSNS  = EXP(GAMLN(BVS+AFAS+1.)-GAMLN(AFAS+1.)-BVS*LLMS)
         VTQS  = MIN(RHOAJ*FSQS*AVS,VTSMAX)
         VTNS  = MIN(RHOAJ*FSNS*AVS,VTSMAX)
         IF (AS1D.GE.ASMALL) THEN
            FSAS = EXP(GAMLN(BVS+AFAS+3.)-GAMLN(AFAS+3.)-BVS*LLMS)
            VTAS = MIN(RHOAJ*FSAS*AVS,VTSMAX)
         ENDIF
         QCLS1 = PI*AMS*NS1D/4.
         NCLS1 = PI*NS1D/4.
         ACLS1 = PI*NS1D/4.
         BSTMP = SCN*SQRT(AVS*RHOAJ/MUA)
         QTMP1 = LLMS*(1.5+BVS/2.)
         QTMP2 = EXP(GAMLN(AFAS+2.)-GAMLN(AFAS+1.)-LLMS)
         QTMP3 = EXP(GAMLN(BVS/2.+AFAS+2.5)-GAMLN(AFAS+1.)-QTMP1)
         CAPS  = ZP1*SASPR**(ZP2/3.)+ZP3*SASPR**(ZP4/3.)
         VENQS = AVSG*QTMP2*CAPS+BVSG*BSTMP*QTMP3*CAPS
         SASR1 = SASPR**(-1./3.)
         SASR2 = SASPR**(-2./3.)
         SASR3 = SASPR**(1./3.)
         SASR4 = SASPR**(2./3.)
         HSdqv = 2.*PI*NS1D*VENQS*XXLS*SSRI0/ABI
         HSwqv = 2.*PI*NS1D*VENQS*XXLV*SSRW0/ABW
      ENDIF
      IF (QG1D.GE.QSMALL) THEN
         CALL SOLVE_AFAG(TK1D,RHO,QG1D,NG1D,VG1D,AG1D,LAMG,AFAG,MVDG,  &
              RHOG,AMG,AVG,BVG,GMLF)
         LLMG  = LOG(LAMG)
         GG2   = EXP(GAMLN(AFAG+2.)-GAMLN(AFAG+1.)-LLMG)
         GG3   = EXP(GAMLN(AFAG+3.)-GAMLN(AFAG+1.)-2.*LLMG)
         GG4   = EXP(GAMLN(AFAG+4.)-GAMLN(AFAG+1.)-3.*LLMG)
         GG5   = EXP(GAMLN(AFAG+5.)-GAMLN(AFAG+1.)-4.*LLMG)
         GGM1  = EXP(GAMLN(AFAG+BMG+1.)-GAMLN(AFAG+1.)-BMG*LLMG)
         GGM2  = EXP(GAMLN(AFAG+BMG+2.)-GAMLN(AFAG+1.)-(BMG+1.)*LLMG)
         GGM3  = EXP(GAMLN(AFAG+BMG+3.)-GAMLN(AFAG+1.)-(BMG+2.)*LLMG)
         iAMG  = 1./AMG
         iRHOG = 1./RHOG
         FSQG  = EXP(GAMLN(BVG+BMG+AFAG+1.)-GAMLN(BMG+AFAG+1.)-BVG*LLMG)
         FSNG  = EXP(GAMLN(BVG+AFAG+1.)-GAMLN(AFAG+1.)-BVG*LLMG)
         VTQG  = MIN(RHOAJ*FSQG*AVG,VTGMAX)
         VTNG  = MIN(RHOAJ*FSNG*AVG,VTGMAX)
         IF (AG1D.GE.ASMALL) THEN
            FSAG = EXP(GAMLN(BVG+AFAG+3.)-GAMLN(AFAG+3.)-BVG*LLMG)
            VTAG = MIN(RHOAJ*FSAG*AVG,VTGMAX)
         ENDIF
         QCLG1 = PI*AMG*NG1D/4.
         NCLG1 = PI*NG1D/4.
         ACLG1 = PI*NG1D/4.
         BGTMP = SCN*SQRT(AVG*RHOAJ/MUA)
         QTMP1 = LLMG*(1.5+BVG/2.)
         QTMP2 = EXP(GAMLN(AFAG+2.)-GAMLN(AFAG+1.)-LLMG)
         QTMP3 = EXP(GAMLN(BVG/2.+AFAG+2.5)-GAMLN(AFAG+1.)-QTMP1)
         VENQG = AVSG*QTMP2+BVSG*BGTMP*QTMP3
         HGdqv = 2.*PI*NG1D*VENQG*XXLS*SSRI0/ABI
         HGwqv = 2.*PI*NG1D*VENQG*XXLV*SSRW0/ABW
      ENDIF
      IF (QH1D.GE.QSMALL) THEN
         CALL SOLVE_AFAH(TK1D,RHO,QH1D,NH1D,AH1D,LAMH,AFAH,MVDH,AVH,BVH)
         LLMH  = LOG(LAMH)
         GH2   = EXP(GAMLN(AFAH+2.)-GAMLN(AFAH+1.)-LLMH)
         GH3   = EXP(GAMLN(AFAH+3.)-GAMLN(AFAH+1.)-2.*LLMH)
         GH4   = EXP(GAMLN(AFAH+4.)-GAMLN(AFAH+1.)-3.*LLMH)
         GH5   = EXP(GAMLN(AFAH+5.)-GAMLN(AFAH+1.)-4.*LLMH)
         FSQH  = EXP(GAMLN(BVH+BMH+AFAH+1.)-GAMLN(BMH+AFAH+1.)-BVH*LLMH)
         FSNH  = EXP(GAMLN(BVH+AFAH+1.)-GAMLN(AFAH+1.)-BVH*LLMH)
         VTQH  = MIN(RHOAJ*FSQH*AVH,VTHMAX)
         VTNH  = MIN(RHOAJ*FSNH*AVH,VTHMAX)
         IF (AH1D.GE.ASMALL) THEN
            FSAH = EXP(GAMLN(BVH+AFAH+3.)-GAMLN(AFAH+3.)-BVH*LLMH)
            VTAH = MIN(RHOAJ*FSAH*AVH,VTHMAX)
         ENDIF
         BHTMP = SCN*SQRT(AVH*RHOAJ/MUA)
         IF (HAIL_VENT.EQ.0) THEN
            QTMP1 = LLMH*(1.5+BVH/2.)
            QTMP2 = EXP(GAMLN(AFAH+2.)-GAMLN(AFAH+1.)-LOG(LAMH))
            QTMP3 = EXP(GAMLN(BVH/2.+AFAH+2.5)-GAMLN(AFAH+1.)-QTMP1)
            ATMP1 = LLMH*(0.5+BVH/2.)
            ATMP2 = EXP(GAMLN(BVH/2.+AFAH+1.5)-GAMLN(AFAH+1.)-ATMP1)
            VENQH = AVRH*QTMP2+BVRH*BHTMP*QTMP3
            VENAH = AVRH+BVRH*BHTMP*ATMP2
         ELSEIF (HAIL_VENT.EQ.1) THEN
            QTMP1 = LLMH*(1.5+BVH/2.)
            QTMP2 = EXP(GAMLN(AFAH+2.)-GAMLN(AFAH+1.)-LOG(LAMH))
            QTMP3 = EXP(GAMLN(BVH/2.+AFAH+2.5)-GAMLN(AFAH+1.)-QTMP1)
            QTMP4 = LLMH*(2.+BVH)
            QTMP5 = EXP(GAMLN(BVH+AFAH+3.)-GAMLN(AFAH+1.)-QTMP4)
            ATMP1 = LLMH*(0.5+BVH/2.)
            ATMP2 = EXP(GAMLN(BVH/2.+AFAH+1.5)-GAMLN(AFAH+1.)-ATMP1)
            ATMP3 = LLMH*(1.+BVH)
            ATMP4 = EXP(GAMLN(BVH+AFAH+2.)-GAMLN(AFAH+1.)-ATMP3)
            VENQH = QTMP2+VENH1*BHTMP*QTMP3+VENH2*BHTMP**2.*QTMP5
            VENAH = 1.+VENH1*BHTMP*ATMP2+VENH2*BHTMP**2.*ATMP4
         ENDIF
         HHdqv = 2.*PI*NH1D*VENQH*XXLS*SSRI0/ABI
         HHwqv = 2.*PI*NH1D*VENQH*XXLV*SSRW0/ABW
      ENDIF
!----------------- CHEN&LIU WARM CLOUD PROCESSES ----------------------
      IF (QC1D.GE.QSMALL) THEN
         NCLcc = EXP(AN3+BN3*LOG(NC1D)+CN3*LMVRC**3.)*GUC
         NCNcr = NC1D*NC1D*EXP(AN4+BN4*MVRC+CN4/MVRC)*GUC
         QCNcr = NCNcr*EXP(AQ4+BQ4*QC1D/NC1D)
         NCNcr = MIN(NCNcr,NC1D*iDT)
         QCNcr = MIN(QCNcr,QC1D*iDT)
      ENDIF
      IF (QR1D.GE.QSMALL) THEN
         NCLrr = MIN(EXP(AN6+BN6*LOG(NR1D)+CN6/MVRR)*GUR,NR1D*iDT)
         NBKrr = NCLrr*EXP(AN7+(BN7+CN7*MVRR)*MVRR)
         NBKrc = NCLrr*EXP(AN8+(BN8+CN8*MVRR)*MVRR)
         QBKrc = MIN(NBKrc**BQ8*EXP(AQ8+CQ8*LMVRR),QR1D*iDT)
         NBKrc = MIN(NBKrc,NR1D*iDT)
      ENDIF
      IF (QC1D.GE.QSMALL.AND.QR1D.GE.QSMALL) THEN
         NCLcr = NC1D*NR1D*GUR*EXP(AN5+BN5*LMVRR+CN5*LMVRC)
         QCLcr = NC1D*NR1D*GUR*EXP(AQ5+BQ5*LMVRR+CQ5*LMVRC)
         NCLcr = MIN(NCLcr,NC1D*iDT)
         QCLcr = MIN(QCLcr,QC1D*iDT)
      ENDIF
      IF (TK1D.LT.TK0C) THEN
         IF ((HIdqv+HSdqv+HGdqv+HHdqv).GE.QSMALL) THEN
            VDMAX  = XXLS*(QV1D-QVSI0)/(1.+XXLS**2.*QV1D/(CPM*RV*      &
                     TK1D**2.))*iDT
            SUMDEP = HIdqv+HSdqv+HGdqv+HHdqv
            IF (SUMDEP.GT.VDMAX.AND.VDMAX.GE.QSMALL) THEN
               RATIO = MIN(1.,VDMAX/(SUMDEP+QSMALL))
               HIdqv = HIdqv*RATIO
               HSdqv = HSdqv*RATIO
               HGdqv = HGdqv*RATIO
               HHdqv = HHdqv*RATIO
            ENDIF
         ENDIF
         IF ((HIdqv+HSdqv+HGdqv+HHdqv).LT.0.) THEN
            SBMAX  = XXLS*(QV1D-QVSI0)/(1.+XXLS**2.*QV1D/(CPM*RV*      &
                     TK1D**2.))*iDT
            SUMSUB = HIdqv+HSdqv+HGdqv+HHdqv
            IF (SBMAX.LT.0..AND.SUMSUB.LT.SBMAX*0.9999) THEN
               HIdqv = HIdqv*MIN(1.,0.9999*SBMAX/SUMSUB)
               HSdqv = HSdqv*MIN(1.,0.9999*SBMAX/SUMSUB)
               HGdqv = HGdqv*MIN(1.,0.9999*SBMAX/SUMSUB)
               HHdqv = HHdqv*MIN(1.,0.9999*SBMAX/SUMSUB)
            ENDIF
         ENDIF
!---------- HOMO/HETER FREEZING OF DROPLETS AND RAIN DROPS  ------------
         IF (QC1D.GE.QSMALL) THEN                                       ! DeMott et al. (1994)
            TMP1  = TC1D*TC1D
            VOLMC = C4PI3*(MVRC*1.E6)**3.
            IJHOF = (10.**MAX(-20.,(-606.3952-52.6611*TC1D-1.7439*TMP1-&
                    0.0265*TMP1*TC1D-1.536e-4*TMP1**2.)))
            RFZ   = 1.-EXP(-IJHOF*VOLMC*DT)
            IF (TK1D.GT.243.15) RFZ = 0.
            IF (TK1D.LT.223.15) RFZ = 1.
            NHOci = RFZ*NC1D*iDT
            QHOci = RFZ*QC1D*iDT
            VHOci = QHOci*iRHOI0
            FHOci = QHOci*1.*iAMI0
            IHOci = QHOci*iAMI0
            AHOci = (KCIMIN*NHOci*IHOci**2.)**THRD
         ENDIF
         IF (TK1D.LT.269.15.AND.QC1D.GE.QSMALL) THEN
            QNMci = QFZC1*BIMM*EXP(AIMM*(TK0C-TK1D)-1.)
            NNMci = NFZC1*BIMM*EXP(AIMM*(TK0C-TK1D)-1.)
            QNMci = MIN(QNMci,QC1D*iDT)
            NNMci = MIN(NNMci,NC1D*iDT)
            VNMci = QNMci*iRHOI0
            FNMci = QNMci*1.*iAMI0                                      ! ISOMETRIC
            INMci = QNMci*iAMI0
            ANMci = (KCIMIN*NNMci*INMci**2.)**THRD
         ENDIF
         IF (QC1D.GE.QSMALL.AND.TK1D.LT.271.15) THEN
            RGIMF  = (2.*2.99E-26*SFCTNW)/(ICED*BOLTZ*TK1D*ACTW)
            GGCNT  = C4PI3*SFCTNV*RGIMF**2.
            EPA    = ESI*QV1D/QVSI
            GEOF2  = EXP(DC1+DC2*COSM2+DC3*RGIMF/INR0)
            NGCNT0 = EXP((-DACTE-GEOF2*GGCNT)/(BOLTZ*TK1D))
            NGCNT1 = EPA/(1.61E-11*TK1D**(0.5))                         ! 1.61E-11 = vs*(2*PI*mw*k)**0.5
!            NGCNT  = MAX(4.*PI*INR0**2.*NGCNT0*NGCNT1,1.)*QNIN
            CNTGG0 = (XXLV**2.-XXLV*RV*TK1D)/(KAP*RV*TK1D**2.)
            CNTGG  = (RV*TK1D/(ESW*DV)+CNTGG0)**(-1.)/RHOW
            KDIFF  = 9.1018E-11*TK1D**2.+8.8197E-8*TK1D-1.0654E-5       ! AIR THERMAL DIFFUSIVITY (COTTON ET AL., 1986)
            TCC    = TC1D+CNTGG*SSRW*XXLV/KDIFF                         ! DROPLETS TEMPERATURE
            NGCNT  = EXP(4.11-0.262*TCC)
            CNTKN  = 6.6E-8*TK1D*101325./(293.15*P1D*INR0)              ! KNUDSEN NUMBER
            PSIA   = -1.*BOLTZ*TCC*(1.+CNTKN)/(6.*PI*INR0*MUA)          ! AEROSOL DIFFUSIVITY
            KAPA   = 5.39E5                                             ! AEROSOL THERMAL CONDUCTIVITY
            CNTFT0 = 0.4*(1.+1.45*CNTKN+0.4*CNTKN*EXP(-1./CNTKN))
            CNTFT1 = (1.+3.*CNTKN)*(2.*KAP+5.*KAPA*CNTKN+KAPA)
            CNTFT  = CNTFT0*(KAP+2.5*CNTKN*KAPA)/CNTFT1
            CNTF1  = 4.*PI*NGCNT*NC1D*MVRC
            CNTF2  = KAP*(TC1D-TCC)/P1D
            IJCNT1 = -1.*CNTF1*CNTF2*RV*TK1D/(XXLV*RHO)                 ! DIFFUSIOPHORESIS
            IJCNT2 = CNTF1*CNTF2*CNTFT*iRHO                             ! THERMEOPHORESIS
            IJCNT3 = CNTF1*PSIA                                         ! BROWNIAN DIFFUSION
!            NNCci  = MIN(MAX(0.,IJCNT1+IJCNT2+IJCNT3),NC1D*iDT)
!            QNCci  = NNCci*QC1D/NC1D
!            VNCci  = QNCci*iRHOI0
!            FNCci  = QNCci*1.*iAMI0                                     ! ISOMETRIC
!            INCci  = QNCci*iAMI0
!            ANCci  = (KCIMIN*NNCci*INCci**2.)**THRD
         ENDIF
         IF (TK1D.LE.233.15.AND.QR1D.GE.QSMALL) THEN
            QHOrg = QR1D*iDT
            NHOrg = NR1D*iDT
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               VHOrg = QHOrg/RHOG0
            ENDIF
            AHOrg = (KCGMAX*NHOrg*(QHOrg*V2M3/RHOG0)**2.)**THRD
         ENDIF
         IF (TK1D.LT.269.15.AND.QR1D.GE.QSMALL) THEN
            QNMrg = QFZR1*BIMM*EXP(AIMM*(TK0C-TK1D)-1.)
            NNMrg = NFZR1*BIMM*EXP(AIMM*(TK0C-TK1D)-1.)
            QNMrg = MIN(QNMrg,QR1D*iDT)
            NNMrg = MIN(NNMrg,NR1D*iDT)
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               VNMrg = QNMrg/RHOG0
            ENDIF
            ANMrg = AFZR1*BIMM*EXP(AIMM*(TK0C-TK1D)-1.)
         ENDIF
!------------------ AUTO-CONVERSION PROCESSES -------------------------
         IF (QI1D.GE.QSMALL.AND.MVDI.GE.5.E-5) THEN
!            EII1 = 0.1*EXP(0.1*TC1D)                                    ! STRAKA AND MANSELL (2005)
            EII1 = 10.**(3.5E-2*TC1D-0.7)                               ! BASEN ON COTTON ET AL. (1986)
            IF (ICE_RHOI.EQ.1) THEN
               EII2 = 1.-RHOI/RHOI0
            ELSE
               EII2 = 0.
            ENDIF
            EII = MIN(MAX(EII1,EII2,0.),1.)
            IF ((ADAGR-1.).GE.SLIMIT) THEN
               QCNis = QCNI1*VTQI*EII*NI1D*(GIM2H1/DI0Z4+2.*GIMF1/     &
                       DI0Z1+DI0Z2*GIM2G1)
               NiCNis = NCNI1*VTNI*EII*NI1D*(GI2H1/DI0Z4+2.*GIF1/DI0Z1+&
                        DI0Z2*GI2G1)
               IF (AI1D.GE.ASMALL) THEN
                  AiCNis = ACNI1*VTAI*EII*NI1D*(GI2H3/DI0Z4+2.*GIF3/   &
                           DI0Z1+DI0Z2*GI2G3)
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FiCNis = FCNI1*VTFI*EII*NI1D/DI0Z3*(GIZM2H1/DI0Z4+2.*&
                           GIZMF1/DI0Z1+DI0Z2*GIZM2G1)
                  ICNis = FCNI1*VTI3M*EII*NI1D*(GIM2H1/DI0Z4+2.*GIMF1/ &
                          DI0Z1+DI0Z2*GIM2G1)
               ENDIF
            ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
               QCNis = QCNI1*VTQI*EII*NI1D*4.*DI0Z2*GIM2G1
               NiCNis = NCNI1*VTNI*EII*NI1D*4.*DI0Z2*GI2G1
               IF (AI1D.GE.ASMALL) THEN
                  AiCNis = ACNI1*VTAI*EII*NI1D*4.*DI0Z2*GI2G3
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FiCNis = FCNI1*VTFI*EII*NI1D*4.*DI0Z2*GIZM2G1/DI0Z3
                  ICNis = FCNI1*VTI3M*EII*NI1D*4.*DI0Z2*GIM2G1
               ENDIF
            ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
               QCNis = QCNI1*VTQI*EII*NI1D*4.*GIM3
               NiCNis = NCNI1*VTNI*EII*NI1D*4.*GI3
               IF (AI1D.GE.ASMALL) THEN
                 AiCNis = ACNI1*VTAI*EII*NI1D*4.*GI5
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FiCNis = FCNI1*VTFI*EII*NI1D*4.*GIM3
                  ICNis = FCNI1*VTI3M*EII*NI1D*4.*GIM3
               ENDIF
            ENDIF
            QCNis  = MIN(QCNis,QI1D*iDT)
            NiCNis = MIN(NiCNis,NI1D*iDT)
            AiCNis = MIN(AiCNis,AI1D*iDT)
            FiCNis = MIN(FiCNis,FI1D*iDT)
            ICNis  = MIN(ICNis,I3M1D*iDT)
            ViCNis = MIN(QCNis*iRHOI,VI1D*iDT)
            IF (ICE_RHOS.EQ.1.OR.AGG_SHAPE.EQ.1) THEN
               IF (ICE_RHOI.EQ.0.AND.ICE_SHAPE.EQ.0) THEN
                  VsCNis = QCNis/RHOS0
                  FsCNis = QCNis*MAX(0.,-1.*LOG(0.8))
                  RATIO  = (RHOI/RHOS0)**THRD
                  NsCNis = QCNis*V2M3/(RHOS0*(RATIO*1.2599*MVDI)**3.)
               ELSE
                  LIM1 = 0.7*(1.-RHOI/RHOI0)                            ! CONTACT ANGLE = 45
                  LIM2 = SQRT(1.-LIM1*LIM1)
                  DICC = DI0**(-6.*ZETA)*GI3H1
                  DICA = GI2HG1/DI0Z3
                  DIAC = GIH2G1
                  DIAA = DI0Z3*GI3G1
                  DIA0 = DI0Z1*GIG1
                  DIC0 = GIH1/DI0Z2
                  IF ((ADAGR-1.).GE.SLIMIT) THEN
                     LIM3  = LIM1**2.*(1.+ISEPR*LIM2)
                     LIM4  = LIM1*(2.+2.*LIM2+2.*ISEPR*LIM2+LIM1**2.*  &
                             ISEPR+2.*LIM2**2.*ISEPR)
                     LIM5  = 1.+2.*LIM2+LIM2*ISEPR+LIM2**2.+2.*LIM1**  &
                             2.*ISEPR+2.*LIM2**2.*ISEPR+2.*LIM1**2.*   &
                             LIM2*ISEPR+LIM2**3.*ISEPR
                     LIM6  = LIM1*ISEPR*(1.+2*LIM2+LIM2**2.)
                     DNIAG = 2.*GI4*RHOI/(LIM3*DICC+LIM4*DICA+LIM5*    &
                             DIAC+LIM6*DIAA)
                     SASP1 = MIN((1.+LIM2)*DIA0+LIM1*DIC0,(1.+LIM2*    &
                             ISEPR)*DIC0+LIM1*ISEPR*DIA0)/MAX((1.+     &
                             LIM2)*DIA0+LIM1*DIC0,(1.+LIM2*ISEPR)*     &
                             DIC0+LIM1*ISEPR*DIA0)
                  ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
                     LIM3  = LIM1*(1.+2.*ISEPR*LIM2+LIM2**2.*ISEPR**2.)
                     LIM4  = 1.+LIM2+2.*LIM2*ISEPR+2.*LIM1**2.*ISEPR+  &
                             2.*LIM2**2.*ISEPR+LIM2**2.*ISEPR**2.+     &
                             LIM2**3.*ISEPR**2.+2.*LIM1**2.*LIM2*      &
                             ISEPR**2.
                     LIM5  = LIM1*ISEPR*(2.+2.*LIM2+2.*LIM2*ISEPR+     &
                             LIM1**2.*ISEPR+2.*LIM2**2.*ISEPR)
                     LIM6  = LIM1**2.*ISEPR**2.*(1.+LIM2)
                     DNIAG = 2.*GI4*RHOI/(LIM3*DIAA+LIM4*DIAC+LIM5*    &
                             DICA+LIM6*DICC)
                     SASP1 = MIN((1.+LIM2)*DIC0+LIM1*DIA0,(1.+LIM2*    &
                             ISEPR)*DIA0+LIM1*ISEPR*DIC0)/MAX((1.+     &
                             LIM2)*DIC0+LIM1*DIA0,(1.+LIM2*ISEPR)*     &
                             DIA0+LIM1*ISEPR*DIC0)
                  ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
                     LIM3  = 1.+LIM1+LIM2+2.*LIM1*ISEPR+2.*LIM2*ISEPR+ &
                             2.*LIM1**2.*ISEPR+2.*LIM2**2.*ISEPR
                     LIM4  = 4.*LIM1*LIM2*ISEPR+LIM1**2.*ISEPR**2.+    &
                             LIM2**2.*ISEPR**2.+2.*LIM1*LIM2*ISEPR**2.
                     LIM5  = LIM1**3.*ISEPR**2.+LIM2**3.*ISEPR**2.
                     LIM6  = 3.*LIM1*LIM2**2.*ISEPR**2.+3.*LIM1**2.*   &
                             LIM2*ISEPR**2.
                     DNIAG = 2.*RHOI/(LIM3+LIM4+LIM5+LIM6)
                     SASP1 = MIN(1.+LIM2*ISEPR+LIM1*ISEPR,1.+LIM2+LIM1)&
                             /MAX(1.+LIM2*ISEPR+LIM1*ISEPR,1.+LIM2+LIM1)
                  ENDIF
                  DNIAG  = MIN(MAX(DNIAG,RHOIMIN),RHOS0)
                  SASP1  = MIN(MAX(SASP1,SASMIN),SASMAX)
                  VsCNis = QCNis/DNIAG
                  FsCNis = QCNis*MAX(0.,-1.*LOG(SASP1))
                  RATIO  = (RHOI/DNIAG)**THRD
                  NsCNis = QCNis*V2M3/(DNIAG*(RATIO*1.2599*MVDI)**3.)
               ENDIF
            ELSE
               RATIO  = (RHOI/RHOS)**THRD
               NsCNis = QCNis*V2M3/(RHOS*(RATIO*1.2599*MVDI)**3.)
            ENDIF
            IF ((ADAGR-1.).GE.SLIMIT) THEN
               AsCNis = ACNI1*VTAI*EII*NI1D*(GI2H3/DI0Z4+2.*GIF3/DI0Z1+&
                        DI0Z2*GI2G3)*(1.5874*RATIO**2.-1.)
            ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
               AsCNis = ACNI1*VTAI*EII*NI1D*4.*DI0Z2*GI2G3*(1.5874*    &
                        RATIO**2.-1.)
            ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
               AsCNis = ACNI1*VTAI*EII*NI1D*4.*GI5*(1.5874*RATIO**2.-1.)
            ENDIF
            AsCNis = MIN(AsCNis,AI1D*iDT)
         ENDIF
         IF (QG1D.GE.QSMALL.AND.W1D.GT.3.) THEN
            DSLL = 2.*1.E-2*(EXP(MIN(20.,-TC1D/(1.1E4*RHO*(QC1D+QR1D)- &
                   1.3E3*RHO*QI1D+1.)))-1.)
            DSLL = MIN(1.,MAX(1.E-4,DSLL))
            IF (AFAG.LE.20.) THEN
               RATIO = MIN(1.,MAX(0.,ABS(GAMMP(AFAG+1.,DSLL*LAMG))))
               NGTAL = NG1D*(1.-RATIO)
               IF (NGTAL.GE.NSMALL) THEN
                  QCNgh  = MIN((1.-RATIO)*QG1D*iDT,QG1D*iDT)
                  NgCNgh = MIN((1.-RATIO)*NG1D*iDT,NG1D*iDT)
                  IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
                     VCNgh = MIN((1.-RATIO)*VG1D*iDT,VG1D*iDT)
                  ENDIF
                  IF (AG1D.GE.ASMALL) THEN
                     AgCNgh = MIN((1.-RATIO)*AG1D*iDT,AG1D*iDT)
                  ENDIF
               ENDIF
            ELSE
               IF (MVDG.GE.DSLL) THEN
                  QCNgh  = QG1D*iDT
                  NgCNgh = NG1D*iDT
                  IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
                     VCNgh = VG1D*iDT
                  ENDIF
                  IF (AG1D.GE.ASMALL) THEN
                     AgCNgh = AG1D*iDT
                  ENDIF
               ENDIF
            ENDIF
            NhCNgh = MIN(QCNgh/(AMH*DSLL**BMH),NG1D*iDT)
            IF (AG1D.GE.ASMALL) THEN
               AhCNgh = (KCHMIN*NhCNgh*(QCNgh*iAMH)**2.)**THRD
            ENDIF
         ENDIF
!----------------- RIMING OF CLOUD DROPLETS ----------------------------
         IF (QI1D.GE.QSMALL.AND.QC1D.GE.QSMALL) THEN
            VTQ0  = VTQI-VTQC
            VTN0  = VTNI-VTNC
            VTQIC = SQRT(VTQ0*VTQ0+0.04*VTQI*VTQC)
            VTNIC = SQRT(VTN0*VTN0+0.04*VTNI*VTNC)
            IF (ICE_RHOG.EQ.0) THEN
               DNIRM = RHOG1
            ELSEIF (ICE_RHOG.EQ.1) THEN
               DNRI  = MIN(MAX(-5.E5*MVDC*VTQIC/TC1D,0.),6.)
               DNIRM = 1.E3*(0.078+0.184*DNRI-0.015*DNRI**2.)
            ELSEIF (ICE_RHOG.EQ.2) THEN
               DNIRM = 3.E2*(-5.E5*MVDC*0.6*VTQIC/TC1D)**0.44
            ENDIF
            DNIRM = MIN(MAX(DNIRM,RHOIMIN),RHOG0)
            RHOIW = (QI1D+QC1D)/(QI1D/RHOI+QC1D/RHOW+ISMALL)
            RATIO = (RHOIW/RHOI)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            IF ((ADAGR-1.).GE.SLIMIT.AND.MVDI.GE.1.E-4) THEN
               WBIN = MIN(MAX(NINT(MVDC*1.E6/10.),0),20)
               CNRE = VTQI*DI0**(-2.*ZETA)*MVDI**IPH*RHO/MUA
               IF (CNRE.LT.0.4) CBIN = 0
               IF (CNRE.GE.0.4.AND.CNRE.LT.0.6) CBIN = 1
               IF (CNRE.GE.0.6.AND.CNRE.LT.0.9) CBIN = 2
               IF (CNRE.GE.0.9.AND.CNRE.LT.1.5) CBIN = 3
               IF (CNRE.GE.1.5.AND.CNRE.LT.3.5) CBIN = 4
               IF (CNRE.GE.3.5.AND.CNRE.LT.7.5) CBIN = 5
               IF (CNRE.GE.7.5.AND.CNRE.LT.15.) CBIN = 6
               IF (CNRE.GE.15.) CBIN = 7
               CBIN  = MIN(MAX(CBIN,0),7)
               ECI   = MAX(0.,IECC(CBIN*21+WBIN))
               QRMci = QRMC1*ECI*VTQIC*NI1D*(GIF1*GC4/DI0Z1+GIG1*GC5*  &
                       DI0Z1+GIH1*GC5/DI0Z2+GC6)
               NRMci = NRMC1*ECI*VTNIC*NI1D*(GIF1/DI0Z1+GIG1*GC2*DI0Z1+&
                       GIH1*GC2/DI0Z2+GC3)
               IF (AI1D.GE.ASMALL.AND.(AX*MVDC).GE.MVDI) THEN
                  ARMci = ARMC1*ECI*VTQIC*NI1D*(ZZ4*GIF1*GC3/DI0Z1+ZZ4*&
                          DI0Z1*GIG1*GC4+ZZ4*GIH1*GC4/DI0Z2+ZZ4*GC5+   &
                          ZZ5*GIF2*GC2/DI0Z1+ZZ5*DI0Z1*GIG2*GC3+ZZ5*   &
                          GIH2*GC3/DI0Z2+ZZ5*GI2*GC4+ZZ6*GIF3/DI0Z1+   &
                          ZZ6*DI0Z1*GIG3*GC2+ZZ6*GIH3*GC2/DI0Z2+ZZ6*   &
                          GI3*GC3)
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FRMci = FRMC1*ECI*VTQIC*NI1D*3./DNIRM/8./DI0Z3*(     &
                          GIZF1*GC4/DI0Z1+GIZG1*GC5*DI0Z1+GIZH1*GC5/   &
                          DI0Z2+GIZ1*GC6)
                  IRMci = IRMC1*ECI*VTQIC*NI1D/DNIRM*(GIF1*GC4/DI0Z1+  &
                          GIG1*GC5*DI0Z1+GIH1*GC5/DI0Z2+GC6)
               ENDIF
            ELSEIF ((1.-ADAGR).GE.SLIMIT.AND.MVDI.GE.1.5E-4) THEN
               WBIN  = MIN(MAX(NINT(MVDC*1.E6/10.),0),20)
               PNRE  = VTQI*DI0**ZETA*MVDI**IPG*RHO/MUA
               IF (PNRE.LT.1.5) PBIN = 0
               IF (PNRE.GE.1.5.AND.PNRE.LT.6.0) PBIN = 1
               IF (PNRE.GE.6.0.AND.PNRE.LT.15.) PBIN = 2
               IF (PNRE.GE.15..AND.PNRE.LT.28.) PBIN = 3
               IF (PNRE.GE.28..AND.PNRE.LT.43.) PBIN = 4
               IF (PNRE.GE.43..AND.PNRE.LT.70.) PBIN = 5
               IF (PNRE.GE.70..AND.PNRE.LT.105.) PBIN = 6
               IF (PNRE.GE.105.) PBIN = 7
               PBIN  = MIN(MAX(PBIN,0),7)
               ECI   = MAX(0.,IEPC(PBIN*21+WBIN))
               QRMci = QRMC1*ECI*VTQIC*NI1D*(GI2G1*GC4*DI0Z2+2.*GIG1*  &
                       GC5*DI0Z1+GC6)
               NRMci = NRMC1*ECI*VTNIC*NI1D*(GI2G1*DI0Z2+2.*GIG1*GC2*  &
                       DI0Z1+GC3)
               IF (AI1D.GE.ASMALL.AND.(AX*MVDC).GE.MVDI) THEN
                  ARMci = ARMC1*ECI*VTQIC*NI1D*(ZZ4*DI0Z2*GI2G1*GC3+2.*&
                          ZZ4*DI0Z1*GIG1*GC4+ZZ4*GC5+ZZ5*DI0Z2*GI2G2*  &
                          GC2+2.*ZZ5*DI0Z1*GIG2*GC3+ZZ5*GI2*GC4+ZZ6*   &
                          DI0Z2*GI2G3+2.*ZZ6*DI0Z1*GIG3*GC2+ZZ6*GI3*GC3)
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FRMci = FRMC1*ECI*VTQIC*NI1D/DNIRM/2./DI0Z3*(GIZ2G1* &
                          GC4*DI0Z2+2.*GIZG1*GC5*DI0Z1+GIZ1*GC6)
                  IRMci = IRMC1*ECI*VTQIC*NI1D/DNIRM*(GI2G1*GC4*DI0Z2+ &
                          2.*GIG1*GC5*DI0Z1+GC6)
               ENDIF
            ELSEIF (ABS(ADAGR-1.).LT.SLIMIT.AND.MVDI.GE.2.E-4) THEN
               ECI   = 0.7
               QRMci = QRMC1*ECI*VTQIC*NI1D*(GI3*GC4+2.*GI2*GC5+GC6)
               NRMci = NRMC1*ECI*VTNIC*NI1D*(GI3+2.*GI2*GC2+GC3)
               IF (AI1D.GE.ASMALL.AND.(AX*MVDC).GE.MVDI) THEN
                  ARMci = ARMC1*ECI*VTQIC*NI1D*(ZZ4*GI3*GC3+2.*ZZ4*GI2*&
                          GC4+ZZ4*GC5+ZZ5*GI4*GC2+2.*ZZ5*GI3*GC3+ZZ5*  &
                          GI2*GC4+ZZ6*GI5+2.*ZZ6*GI4*GC2+ZZ6*GI3*GC3)
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FRMci = FRMC1*ECI*VTQIC*NI1D/DNIRM/4.*(GI3*GC4+2.*   &
                          GI2*GC5+GC6)
                  IRMci = IRMC1*ECI*VTQIC*NI1D/DNIRM*(GI3*GC4+2.*GI2*  &
                          GC5+GC6)
               ENDIF
            ENDIF
            QRMci = MIN(QRMci,QC1D*iDT)
            NRMci = MIN(NRMci,NC1D*iDT)
            VRMci = QRMci/DNIRM
            FRMci = MIN(FRMci,QC1D*iDT/AMW)
            IRMci = MIN(IRMci,QC1D*iDT/AMW)
            ARMci = MIN(ARMci,QC1D*iDT*iAPW/MVDC)
            IF (QRMci.GT.0.) THEN
               QINig  = MIN(2.*QRMci,QI1D*iDT)
               NINig  = MIN(QINig*NI1D/QI1D,NI1D*iDT)
               ViINig = MIN(VRMci+QRMci*VI1D/QI1D,VI1D*iDT)
               IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
                  VgINig = 2.*VRMci
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FINig = MIN(FRMci+QRMci*FI1D/QI1D,FI1D*iDT)
                  IINig = MIN(IRMci+QRMci*I3M1D/QI1D,I3M1D*iDT)
               ENDIF
               IF (AI1D.GE.ASMALL) THEN
                  AiINig = MIN(ARMci+QRMci*AI1D/QI1D,AI1D*iDT)
                  AgINig = (KCGMAX*NINig*(QINig*V2M3/DNIRM)**2.)**THRD
               ENDIF
            ENDIF
         ENDIF
         IF (QS1D.GE.QSMALL.AND.QC1D.GE.QSMALL.AND.MVDS.GE.1.5E-4) THEN
            ECS   = MIN(MVDC,3.E-5)*3.333E4*SQRT(MIN(MVDS,1.E-3)*1.E3)  ! Approximation based on Pruppacher&Klett(1997) Fig.14-11
            VTQ0  = VTQS-VTQC
            VTN0  = VTNS-VTNC
            VTQSC = SQRT(VTQ0*VTQ0+0.04*VTQS*VTQC)
            VTNSC = SQRT(VTN0*VTN0+0.04*VTNS*VTNC)
            IF (ICE_RHOG.EQ.0) THEN
               DNSRM = RHOG1
            ELSEIF (ICE_RHOG.EQ.1) THEN
               DNRI  = MIN(MAX(-5.E5*MVDC*VTQSC/TC1D,0.),6.)
               DNSRM = 1.E3*(0.078+0.184*DNRI-0.015*DNRI**2.)
            ELSEIF (ICE_RHOG.EQ.2) THEN
               DNSRM = 3.E2*(-5.E5*MVDC*0.6*VTQSC/TC1D)**0.44
            ENDIF
            DNSRM = MIN(MAX(DNSRM,RHOIMIN),RHOG0)
            RHOSW = (QS1D+QC1D)/(QS1D/RHOS+QC1D/RHOW+ISMALL)
            RATIO = (RHOSW/RHOS)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            QRMcs = QRMC1*ECS*VTQSC*NS1D*(SASR2*GS3*GC4+SASR1*2.*GS2*  &
                    GC5+GC6)
            NRMcs = NRMC1*ECS*VTNSC*NS1D*(SASR2*GS3+SASR1*2.*GS2*GC2+  &
                    GC3)
            QRMcs = MIN(QRMcs,QC1D*iDT)
            NRMcs = MIN(NRMcs,NC1D*iDT)
            IF (ICE_RHOS.EQ.1.OR.AGG_SHAPE.EQ.1) THEN
               IF (ICE_RHOI.EQ.0.AND.ICE_SHAPE.EQ.0) THEN
                  VRMcs = QRMcs/RHOS0
                  FRMcs = QRMcs*MAX(0.,-1.*LOG(0.8))
               ELSE
                  VRMcs = QRMcs/DNSRM
                  FRMcs = QRMcs*MAX(0.,-1.*LOG(SASPR))
               ENDIF
            ENDIF
            IF (AS1D.GE.ASMALL.AND.(AX*MVDC).GE.MVDS) THEN
               ARMcs = ARMC1*ECS*VTQSC*NS1D*(ZZ4*SASR2*GS3*GC3+2.*ZZ4* &
                       SASR1*GS2*GC4+ZZ4*GC5+ZZ5*SASR2*GS4*GC2+2.*ZZ5* &
                       SASR1*GS3*GC3+ZZ5*GS2*GC4+ZZ6*SASR2*GS5+2.*ZZ6* &
                       SASR1*GS4*GC2+ZZ6*GS3*GC3)
               ARMcs = MIN(ARMcs,QC1D*iDT*iAPW/MVDC)
            ENDIF
            IF (QRMcs.GT.0.) THEN
               NIMcsi = 3.5E8*QRMcs*FF
               QIMcsi = MIN(NIMcsi*MI0,QRMcs)
               VIMcsi = QIMcsi*iRHOI0
               FIMcsi = QIMcsi*1.*iAMI0
               IIMcsi = QIMcsi*iAMI0
               AIMcsi = (KCIMIN*NIMcsi*IIMcsi**2.)**THRD
               QINsg  = MIN(2.*QRMcs,QS1D*iDT)
               NINsg  = MIN(QINsg*NS1D/QS1D,NS1D*iDT)
               IF (ICE_RHOS.EQ.1) THEN
                  VsINsg = MIN(VRMcs+QRMcs*VS1D/QS1D,VS1D*iDT)
               ENDIF
               IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
                  VgINsg = 2.*VRMcs
               ENDIF
               IF (AGG_SHAPE.EQ.1) THEN
                  FINsg = MIN(QINsg*MAX(0.,-1.*LOG(SASPR)),FS1D*iDT)
               ENDIF
               IF (AS1D.GE.ASMALL) THEN
                  AsINsg = MIN(ARMcs+QRMcs*AS1D/QS1D,AS1D*iDT)
                  AgINsg = (KCGMAX*NINsg*(QINsg*V2M3/DNSRM)**2.)**THRD
               ENDIF
            ENDIF
         ENDIF
         IF (QG1D.GE.QSMALL.AND.QC1D.GE.QSMALL.AND.MVDG.GE.2.E-4) THEN
            STOKE = (RHOW*VTQG*MVDC**2.)/(9.*MUA*MVDG)
            STOKE = MAX(1.5,MIN(10.,STOKE))
            ECG   = 5.5E-1*LOG10(2.51*STOKE)                            ! parameterization based on Cober and List,1993 [JAS]
            VTQ0  = VTQG-VTQC
            VTN0  = VTNG-VTNC
            VTQGC = SQRT(VTQ0*VTQ0+0.04*VTQG*VTQC)
            VTNGC = SQRT(VTN0*VTN0+0.04*VTNG*VTNC)
            IF (ICE_RHOG.EQ.0) THEN
               DNGRM = RHOG1
            ELSEIF (ICE_RHOG.EQ.1) THEN
               DNRI  = MIN(MAX(-5.E5*MVDC*VTQGC/TC1D,0.),6.)
               DNGRM = 1.E3*(0.078+0.184*DNRI-0.015*DNRI**2.)
            ELSEIF (ICE_RHOG.EQ.2) THEN
               DNGRM = 3.E2*(-5.E5*MVDC*0.6*VTQGC/TC1D)**0.44
            ENDIF
            DNGRM = MIN(MAX(DNGRM,RHOIMIN),RHOG0)
            RHOGW = (QG1D+QC1D)/(QG1D/RHOG+QC1D/RHOW+ISMALL)
            RATIO = (RHOGW/RHOG)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            QRMcg = QRMC1*ECG*VTQGC*NG1D*(GG3*GC4+2.*GG2*GC5+GC6)
            NRMcg = NRMC1*ECG*VTNGC*NG1D*(GG3+2.*GG2*GC2+GC3)
            QRMcg = MIN(QRMcg,QC1D*iDT)
            NRMcg = MIN(NRMcg,NC1D*iDT)
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               VRMcg = QRMcg/DNGRM
            ENDIF
            IF (AG1D.GE.ASMALL.AND.(AX*MVDC).GE.MVDG) THEN
               ARMcg = ARMC1*ECG*VTQGC*NG1D*(ZZ4*GG3*GC3+2.*ZZ4*GG2*   &
                       GC4+ZZ4*GC5+ZZ5*GG4*GC2+2.*ZZ5*GG3*GC3+ZZ5*GG2* &
                       GC4+ZZ6*GG5+2.*ZZ6*GG4*GC2+ZZ6*GG3*GC3)
               ARMcg = MIN(ARMcg,QC1D*iDT*iAPW/MVDC)
            ENDIF
            IF (QRMcg.GT.0.) THEN
               NIMcgi = 3.5E8*QRMcg*FF
               QIMcgi = MIN(NIMcgi*MI0,QRMcg)
               VIMcgi = QIMcgi*iRHOI0
               FIMcgi = QIMcgi*1.*iAMI0
               IIMcgi = QIMcgi*iAMI0
               AIMcgi = (KCIMIN*NIMcgi*IIMcgi**2.)**THRD
            ENDIF
         ENDIF
         IF (QH1D.GE.QSMALL.AND.QC1D.GE.QSMALL) THEN
            ECH   = EXP(-8.68E-7*MVDC**(-1.6)*MVDH)                     ! Ziegler (1985) A24
            VTQ0  = VTQH-VTQC
            VTN0  = VTNH-VTNC
            VTQHC = SQRT(VTQ0*VTQ0+0.04*VTQH*VTQC)
            VTNHC = SQRT(VTN0*VTN0+0.04*VTNH*VTNC)
            RHOHW = (QH1D+QC1D)/(QH1D/RHOH+QC1D/RHOW+ISMALL)
            RATIO = (RHOHW/RHOH)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            QRMch = QRMC1*ECH*VTQHC*NH1D*(GH3*GC4+2.*GH2*GC5+GC6)
            NRMch = NRMC1*ECH*VTNHC*NH1D*(GH3+2.*GH2*GC2+GC3)
            QRMch = MIN(QRMch,QC1D*iDT)
            NRMch = MIN(NRMch,NC1D*iDT)
            IF (AH1D.GE.ASMALL.AND.(AX*MVDC).GE.MVDH) THEN
               ARMch = ARMC1*ECH*VTQHC*NH1D*(ZZ4*GH3*GC3+2.*ZZ4*GH2*   &
                       GC4+ZZ4*GC5+ZZ5*GH4*GC2+2.*ZZ5*GH3*GC3+ZZ5*GH2* &
                       GC4+ZZ6*GH5+2.*ZZ6*GH4*GC2+ZZ6*GH3*GC3)
               ARMch = MIN(ARMch,QC1D*iDT*iAPW/MVDC)
            ENDIF
         ENDIF
!----------------- RIMING/COLLECTION OF RAIN DROPS ---------------------
         IF (QR1D.GE.QSMALL.AND.QI1D.GE.QSMALL) THEN
            ERI   = 1.
            VTQ0  = VTQR-VTQI
            VTN0  = VTNR-VTNI
            VTQRI = SQRT(VTQ0*VTQ0+0.04*VTQR*VTQI)
            VTNRI = SQRT(VTN0*VTN0+0.04*VTNR*VTNI)
            IF ((ADAGR-1.).GE.SLIMIT.AND.MVDR.GE.MVDI) THEN
               QCLir = QCLI1*ERI*VTQRI*NR1D*(GIMF1/DI0Z1+GIMG1*GR2*    &
                       DI0Z1+GIMH1*GR2/DI0Z2+GIM1*GR3)
               NCLir = NCLI1*ERI*VTNRI*NR1D*(GIF1/DI0Z1+GIG1*GR2*DI0Z1+&
                       GIH1*GR2/DI0Z2+GR3)
               IF (AI1D.GE.ASMALL) THEN
                  ACLir = ACLI1*ERI*VTQRI*NR1D*(GIF3/DI0Z1+GIG3*GR2*  &
                          DI0Z1+GIH3*GR2/DI0Z2+GI3*GR3)
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FCLir = FCLI1*ERI*VTQRI*NR1D*(GIZMF1/DI0Z1+GIZMG1*   &
                          GR2*DI0Z1+GIZMH1*GR2/DI0Z2+GIZM1*GR3)
                  ICLir = NCLI1*ERI*VTQRI*NR1D*(GIMF1/DI0Z1+GIMG1*GR2* &
                          DI0Z1+GIMH1*GR2/DI0Z2+GIM1*GR3)
               ENDIF
            ELSEIF ((1.-ADAGR).GE.SLIMIT.AND.MVDR.GE.MVDI) THEN
               QCLir = QCLI1*ERI*VTQRI*NR1D*(GIM2G1*DI0Z2+2.*GIMG1*GR2*&
                       DI0Z1+GIM1*GR3)
               NCLir = NCLI1*ERI*VTNRI*NR1D*(GI2G1*DI0Z2+2.*GIG1*GR2*  &
                       DI0Z1+GR3)
               IF (AI1D.GE.ASMALL) THEN
                  ACLir = ACLI1*ERI*VTQRI*NR1D*(GI2G3*DI0Z2+2.*GIG3*  &
                          GR2*DI0Z1+GI3*GR3)
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FCLir = FCLI1*ERI*VTQRI*NR1D*(GIZM2G1*DI0Z2+2.*      &
                          GIZMG1*GR2*DI0Z1+GIZM1*GR3)
                  ICLir = NCLI1*ERI*VTQRI*NR1D*(GIM2G1*DI0Z2+2.*GIMG1* &
                          GR2*DI0Z1+GIM1*GR3)
               ENDIF
            ELSEIF (ABS(ADAGR-1.).LT.SLIMIT.AND.MVDR.GE.MVDI) THEN
               QCLir = QCLI1*ERI*VTQRI*NR1D*(GIM3+2.*GIM2*GR2+GIM1*GR3)
               NCLir = NCLI1*ERI*VTNRI*NR1D*(GI3+2.*GI2*GR2+GR3)
               IF (AI1D.GE.ASMALL) THEN
                  ACLir = ACLI1*ERI*VTQRI*NR1D*(GI5+2.*GI4*GR2+GI3*GR3)
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FCLir = NCLI1*ERI*VTQRI*NR1D*(GIM3+2.*GIM2*GR2+GIM1* &
                          GR3)
                  ICLir = NCLI1*ERI*VTQRI*NR1D*(GIM3+2.*GIM2*GR2+GIM1* &
                          GR3)
               ENDIF
            ENDIF
            QCLir = MIN(QCLir,QI1D*iDT)
            NCLir = MIN(NCLir,NI1D*iDT)
            VCLir = MIN(QCLir*iRHOI,VI1D*iDT)
            ACLir = MIN(ACLir,AI1D*iDT)
            FCLir = MIN(FCLir,FI1D*iDT)
            ICLir = MIN(ICLir,I3M1D*iDT)
            IF (ICE_RHOG.EQ.0) THEN
               DNIRM = RHOG1
            ELSEIF (ICE_RHOG.EQ.1) THEN
               DNRI  = MIN(MAX(-5.E5*MVDR*VTQRI/TC1D,0.),6.)
               DNIRM = 1.E3*(0.078+0.184*DNRI-0.015*DNRI**2.)
            ELSEIF (ICE_RHOG.EQ.2) THEN
               DNIRM = 3.E2*(-5.E5*MVDR*0.6*VTQRI/TC1D)**0.44
            ENDIF
            DNIRM = MIN(MAX(DNIRM,RHOIMIN),RHOG0)
            IF (MVDI.GE.MVDR) THEN
            IF ((ADAGR-1.).GE.SLIMIT.AND.MVDI.GE.1.E-4) THEN
               QRMri = QRMR1*ERI*VTQRI*NI1D*(GIF1*GR4/DI0Z1+GIG1*GR5*  &
                       DI0Z1+GIH1*GR5/DI0Z2+GR6)
               NRMri = NRMR1*ERI*VTNRI*NI1D*(GIF1/DI0Z1+GIG1*GR2*DI0Z1+&
                       GIH1*GR2/DI0Z2+GR3)
            ELSEIF ((1.-ADAGR).GE.SLIMIT.AND.MVDI.GE.1.5E-4) THEN
               QRMri = QRMR1*ERI*VTQRI*NI1D*(GI2G1*GR4*DI0Z2+2.*GIG1*  &
                       GR5*DI0Z1+GR6)
               NRMri = NRMR1*ERI*VTNRI*NI1D*(GI2G1*DI0Z2+2.*GIG1*GR2*  &
                       DI0Z1+GR3)
            ELSEIF (ABS(ADAGR-1.).LT.SLIMIT.AND.MVDI.GE.2.E-4) THEN
               QRMri = QRMR1*ERI*VTQRI*NI1D*(GI3*GR4+2.*GI2*GR5+GR6)
               NRMri = NRMR1*ERI*VTNRI*NI1D*(GI3+2.*GI2*GR2+GR3)
            ENDIF
            QRMri = MIN(QRMri,QR1D*iDT)
            NRMri = MIN(NRMri,NR1D*iDT)
            ENDIF
            IF (MVDI.GE.MVDR.AND.(BX*MVDR).GE.MVDI) THEN
               MVDX = SQRT(ZZ1*MVDR**2.+ZZ2*MVDR*MVDI+ZZ3*MVDI**2.)
            ELSEIF (MVDR.GT.MVDI.AND.(BX*MVDI).GE.MVDR) THEN
               MVDX = SQRT(ZZ1*MVDI**2.+ZZ2*MVDR*MVDI+ZZ3*MVDR**2.)
            ELSE
               MVDX = MAX(MVDI,MVDR)
            ENDIF
            QCLirg = QRMri+QCLir
            RHOIW  = (QI1D+QR1D)/(QI1D/RHOI+QR1D/RHOW+ISMALL)
            RATIO  = (RHOIW/DNIRM)**THRD
            NCLirg = QCLirg*V2M3/(DNIRM*(RATIO*MVDX)**3.)
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               VCLirg = QCLirg/DNIRM
            ENDIF
            ACLirg = (KCIMIN*NCLirg*(QCLirg*V2M3/DNIRM)**2.)**THRD
         ENDIF
         IF (QR1D.GE.QSMALL.AND.QS1D.GE.QSMALL) THEN
            ERS   = 1.
            VTQ0  = VTQR-VTQS
            VTN0  = VTNR-VTNS
            VTQRS = SQRT(VTQ0*VTQ0+0.04*VTQR*VTQS)
            VTNRS = SQRT(VTN0*VTN0+0.04*VTNR*VTNS)
            IF (ICE_RHOG.EQ.0) THEN
               DNSRM = RHOG1
            ELSEIF (ICE_RHOG.EQ.1) THEN
               DNRI  = MIN(MAX(-5.E5*MVDR*VTQRS/TC1D,0.),6.)
               DNSRM = 1.E3*(0.078+0.184*DNRI-0.015*DNRI**2.)
            ELSEIF (ICE_RHOG.EQ.2) THEN
               DNSRM = 3.E2*(-5.E5*MVDR*0.6*VTQRS/TC1D)**0.44
            ENDIF
            DNSRM = MIN(MAX(DNSRM,RHOIMIN),RHOG0)
            IF (MVDS.GE.2.E-4.AND.MVDS.GE.MVDR) THEN
               QRMrs = QRMR1*ERS*VTQRS*NS1D*(SASR2*GS3*GR4+SASR1*2.*   &
                       GS2*GR5+GR6)
               NRMrs = NRMR1*ERS*VTNRS*NS1D*(SASR2*GS3+SASR1*2.*GS2*   &
                       GR2+GR3)
               QRMrs = MIN(QRMrs,QR1D*iDT)
               NRMrs = MIN(NRMrs,NR1D*iDT)
            ENDIF
            IF (MVDR.GE.MVDS) THEN
               QCLsr = QCLS1*ERS*VTQRS*NR1D*(SASR2*GSM3+SASR1*2.*GR2*  &
                       GSM2+GR3*GSM1)
               NCLsr = NCLS1*ERS*VTNRS*NR1D*(SASR2*GS3+SASR1*2.*GR2*   &
                       GS2+GR3)
               QCLsr = MIN(QCLsr,QS1D*iDT)
               NCLsr = MIN(NCLsr,NS1D*iDT)
               IF (ICE_RHOS.EQ.1) THEN
                  VCLsr = MIN(QCLsr*iRHOS,VS1D*iDT)
               ENDIF
               IF (AGG_SHAPE.EQ.1) THEN
                  FCLsr = QCLsr*MAX(0.,-1.*LOG(SASPR))
               ENDIF
               IF (AS1D.GE.ASMALL) THEN
                  ACLsr = ACLS1*ERS*VTQRS*NR1D*(SASR2*GS5+SASR1*2.*GR2*&
                          GS4+GR3*GS3)
                  ACLsr = MIN(ACLsr,AS1D*iDT)
               ENDIF
            ENDIF
            IF (QRMrs.GT.0.) THEN
               NIMrsi = 3.5E8*QRMrs*FF
               QIMrsi = MIN(NIMrsi*MI0,QRMrs)
               VIMrsi = QIMrsi*iRHOI0
               FIMrsi = QIMrsi*1.*iAMI0
               IIMrsi = QIMrsi*iAMI0
               AIMrsi = (KCIMIN*NIMrsi*IIMrsi**2.)**THRD
            ENDIF
            IF (MVDS.GE.MVDR.AND.(BX*MVDR).GE.MVDS) THEN
               MVDX = SQRT(ZZ1*MVDR**2.+ZZ2*MVDR*MVDS+ZZ3*MVDS**2.)
            ELSEIF (MVDR.GT.MVDS.AND.(BX*MVDS).GE.MVDR) THEN
               MVDX = SQRT(ZZ1*MVDS**2.+ZZ2*MVDR*MVDS+ZZ3*MVDR**2.)
            ELSE
               MVDX = MAX(MVDS,MVDR)
            ENDIF
            QCLsrg = QRMrs+QCLsr
            RHOSW  = (QS1D+QR1D)/(QS1D/RHOS+QR1D/RHOW+ISMALL)
            RATIO  = (RHOSW/DNSRM)**THRD
            NCLsrg = QCLsrg*V2M3/(DNSRM*(RATIO*MVDX)**3.)
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               VCLsrg = QCLsrg/DNSRM
            ENDIF
            ACLsrg = (KCIMIN*NCLsrg*(QCLsrg*V2M3/DNSRM)**2.)**THRD
         ENDIF
         IF (QR1D.GE.QSMALL.AND.QG1D.GE.QSMALL) THEN
            VTQ0  = VTQR-VTQG
            VTN0  = VTNR-VTNG
            VTQRG = SQRT(VTQ0*VTQ0+0.04*VTQR*VTQG)
            VTNRG = SQRT(VTN0*VTN0+0.04*VTNR*VTNG)
            STOKE = (RHOW*ABS(VTQ0)*MVDR**2.)/(9.*MUA*MVDG)
            STOKE = MAX(1.5,MIN(10.,STOKE))
            ERG   = 5.5E-1*LOG10(2.51*STOKE)                            ! parameterization based on Cober and List, 1993 [JAS]
            IF (ICE_RHOG.EQ.0) THEN
               DNGRM = RHOG1
            ELSEIF (ICE_RHOG.EQ.1) THEN
               DNRI  = MIN(MAX(-5.E5*MVDR*VTQRG/TC1D,0.),6.)
               DNGRM = 1.E3*(0.078+0.184*DNRI-0.015*DNRI**2.)
            ELSEIF (ICE_RHOG.EQ.2) THEN
               DNGRM = 3.E2*(-5.E5*MVDR*0.6*VTQRG/TC1D)**0.44
            ENDIF
            DNGRM = MIN(MAX(DNGRM,RHOIMIN),RHOG0)
            IF (MVDG.GE.2.E-4.AND.MVDG.GE.MVDR) THEN
               QRMrg = QRMR1*ERG*VTQRG*NG1D*(GG3*GR4+2.*GG2*GR5+GR6)
               NRMrg = NRMR1*ERG*VTNRG*NG1D*(GG3+2.*GG2*GR2+GR3)
               QRMrg = MIN(QRMrg,QR1D*iDT)
               NRMrg = MIN(NRMrg,NR1D*iDT)
            ENDIF
            IF (MVDR.GE.MVDG) THEN
               QCLgr = QCLG1*ERG*VTQRG*NR1D*(GGM3+2.*GR2*GGM2+GR3*GGM1)
               NCLgr = NCLG1*ERG*VTNRG*NR1D*(GG3+2.*GR2*GG2+GR3)
               QCLgr = MIN(QCLgr,QG1D*iDT)
               NCLgr = MIN(NCLgr,NG1D*iDT)
               IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
                  VCLgr = MIN(QCLgr*iRHOG,VG1D*iDT)
               ENDIF
               IF (AG1D.GE.ASMALL) THEN
                  ACLgr = ACLG1*ERG*VTQRG*NR1D*(GG5+2.*GR2*GG4+GR3*GG3)
                  ACLgr = MIN(ACLgr,AG1D*iDT)
               ENDIF
            ENDIF
            IF (QRMrg.GT.0.) THEN
               NIMrgi = 3.5E8*QRMrg*FF
               QIMrgi = MIN(NIMrgi*MI0,QRMrg)
               VIMrgi = QIMrgi*iRHOI0
               FIMrgi = QIMrgi*1.*iAMI0
               IIMrgi = QIMrgi*iAMI0
               AIMrgi = (KCIMIN*NIMrgi*IIMrgi**2.)**THRD
            ENDIF
            IF (MVDG.GE.MVDR.AND.(BX*MVDR).GE.MVDG) THEN
               MVDX = SQRT(ZZ1*MVDR**2.+ZZ2*MVDR*MVDG+ZZ3*MVDG**2.)
            ELSEIF (MVDR.GT.MVDG.AND.(BX*MVDG).GE.MVDR) THEN
               MVDX = SQRT(ZZ1*MVDG**2.+ZZ2*MVDR*MVDG+ZZ3*MVDR**2.)
            ELSE
               MVDX = MAX(MVDG,MVDR)
            ENDIF
            QCLgrg = QRMrg+QCLgr
            RHOGW  = (QG1D+QR1D)/(QG1D/RHOG+QR1D/RHOW+ISMALL)
            RATIO  = (RHOGW/DNGRM)**THRD
            NCLgrg = QCLgrg*V2M3/(DNGRM*(RATIO*MVDX)**3.)
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               VCLgrg = QCLgrg/DNGRM
            ENDIF
            ACLgrg = (KCIMIN*NCLgrg*(QCLgrg*V2M3/DNGRM)**2.)**THRD
         ENDIF
         IF (QH1D.GE.QSMALL.AND.QR1D.GE.QSMALL.AND.MVDH.GE.MVDR) THEN
            ERH   = 1.
            VTQ0  = VTQR-VTQH
            VTN0  = VTNR-VTNH
            VTQRH = SQRT(VTQ0*VTQ0+0.04*VTQR*VTQH)
            VTNRH = SQRT(VTN0*VTN0+0.04*VTNR*VTNH)
            RHOHW = (QH1D+QR1D)/(QH1D/RHOH+QR1D/RHOW+ISMALL)
            RATIO = (RHOHW/RHOH)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            QRMrh = QRMR1*ERH*VTQRH*NH1D*(GH3*GR4+2.*GH2*GR5+GR6)
            NRMrh = NRMR1*ERH*VTNRH*NH1D*(GH3+2.*GH2*GR2+GR3)
            QRMrh = MIN(QRMrh,QR1D*iDT)
            NRMrh = MIN(NRMrh,NR1D*iDT)
            IF (AH1D.GE.ASMALL.AND.(AX*MVDR).GE.MVDH) THEN
               ARMrh = ARMR1*ERH*VTQRH*NH1D*(ZZ4*GH3*GR3+2.*ZZ4*GH2*   &
                       GR4+ZZ4*GR5+ZZ5*GH4*GR2+2.*ZZ5*GH3*GR3+ZZ5*GH2* &
                       GR4+ZZ6*GH5+2.*ZZ6*GH4*GR2+ZZ6*GH3*GR3)
               ARMrh = MIN(ARMrh,QR1D*iDT*iAPW/MVDR)
            ENDIF
         ENDIF
!----------------- ICE-PHASE COLLECTION&CONVERSION*INITILIZATION -------
         IF (QS1D.GE.QSMALL.AND.QI1D.GE.QSMALL.AND.MVDS.GE.MVDI) THEN
!            EIS1  = MIN(1.,0.05*EXP(0.1*TC1D))                          ! FERRIER ET AL., 1995
            EIS1  = 10.**(3.5E-2*TC1D-0.7)
            RHOIS = (QI1D+QS1D)/(QI1D/RHOI+QS1D/RHOS+ISMALL)
            IF (ICE_RHOI.NE.2.AND.ICE_RHOS.NE.0) THEN
               EIS2 = 1.-RHOIS/RHOI0
            ELSE
               EIS2 = 0.
            ENDIF
            EIS   = MIN(MAX(EIS1,EIS2,0.),1.)
            RATIO = (RHOIS/RHOS)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            ZZ7   = ZZ1*RATIO**2.-1.
            ZZ8   = ZZ3*RATIO**2.
            VTQ0  = VTQS-VTQI
            VTN0  = VTNS-VTNI
            VTQIS = SQRT(VTQ0*VTQ0+0.04*VTQI*VTQS)
            VTNIS = SQRT(VTN0*VTN0+0.04*VTNI*VTNS)
            IF (AS1D.GE.ASMALL.OR.AI1D.GE.ASMALL) THEN
               VTAX  = VTAS-VTAI
               VTAIS = SQRT(VTAX*VTAX+0.04*VTAI*VTAS)
            ENDIF
            IF ((ADAGR-1.).GE.SLIMIT) THEN
               QCLis = QCLI1*EIS*VTQIS*NS1D*(GIMF1/DI0Z1+SASR1*GIMG1*  &
                       GS2*DI0Z1+SASR1*GIMH1*GS2/DI0Z2+SASR2*GIM1*GS3)
               NCLis = NCLI1*EIS*VTNIS*NS1D*(GIF1/DI0Z1+SASR1*GIG1*GS2*&
                       DI0Z1+SASR1*GIH1*GS2/DI0Z2+SASR2*GS3)
               IF (AI1D.GE.ASMALL) THEN
                  AiCLis = ACLI1*EIS*VTAIS*NS1D*(GIF3/DI0Z1+SASR1*GIG3*&
                           GS2*DI0Z1+SASR1*GIH3*GS2/DI0Z2+SASR2*GI3*GS3)
               ENDIF
               IF (AS1D.GE.ASMALL) THEN
                  IF (MVDS.GE.MVDI.AND.(AX*MVDI).GE.MVDS) THEN
                     AsCLis = ACLI1*EIS*VTAIS*NS1D*(ZZ4*GIF3/DI0Z1+ZZ4*&
                              SASR1*GIG3*GS2*DI0Z1+ZZ4*SASR1*GIH3*GS2/ &
                              DI0Z2+ZZ4*SASR2*GI3*GS3+ZZ5*GIF2*GS2/    &
                              DI0Z1+ZZ5*SASR1*GIG2*GS3*DI0Z1+ZZ5*SASR1*&
                              GIH2*GS3/DI0Z2+ZZ5*SASR2*GI2*GS4+ZZ6*    &
                              GIF1*GS3/DI0Z1+ZZ6*SASR1*GIG1*GS4*DI0Z1+ &
                              ZZ6*SASR1*GIH1*GS4/DI0Z2+ZZ6*SASR2*GS5)
                  ELSEIF (MVDI.GT.MVDS.AND.(AX*MVDS).GE.MVDI) THEN
                     AsCLis = ACLI1*EIS*VTAIS*NS1D*(ZZ7*GIF1*GS3/DI0Z1+&
                              ZZ7*SASR1*GIG1*GS4*DI0Z1+ZZ7*SASR1*GIH1* &
                              GS4/DI0Z2+ZZ7*SASR2*GS5+ZZ5*GIF2*GS2/    &
                              DI0Z1+ZZ5*SASR1*GIG2*GS3*DI0Z1+ZZ5*SASR1*&
                              GIH2*GS3/DI0Z2+ZZ5*SASR2*GI2*GS4+ZZ8*    &
                              GIF3/DI0Z1+ZZ8*SASR1*GIG3*GS2*DI0Z1+ZZ8* &
                              SASR1*GIH3*GS2/DI0Z2+ZZ8*SASR2*GI3*GS3)
                  ENDIF
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FiCLis = FCLI1*EIS*VTQIS*NS1D*(GIZMF1/DI0Z1+SASR1*   &
                           GIZMG1*GS2*DI0Z1+SASR1*GIZMH1*GS2/DI0Z2+    &
                           SASR2*GIZM1*GS3)
                  ICLis = NCLI1*EIS*VTQIS*NS1D*(GIMF1/DI0Z1+SASR1*     &
                          GIMG1*GS2*DI0Z1+SASR1*GIMH1*GS2/DI0Z2+SASR2* &
                          GIM1*GS3)
               ENDIF
            ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
               QCLis = QCLI1*EIS*VTQIS*NS1D*(GIM2G1*DI0Z2+SASR1*2.*    &
                       GIMG1*GS2*DI0Z1+SASR2*GIM1*GS3)
               NCLis = NCLI1*EIS*VTNIS*NS1D*(GI2G1*DI0Z2+SASR1*2.*GIG1*&
                       GS2*DI0Z1+SASR2*GS3)
               IF (AI1D.GE.ASMALL) THEN
                  AiCLis = ACLI1*EIS*VTAIS*NS1D*(GI2G3*DI0Z2+SASR1*2.* &
                           GIG3*GS2*DI0Z1+SASR2*GI3*GS3)
               ENDIF
               IF (AS1D.GE.ASMALL) THEN
                  IF (MVDS.GE.MVDI.AND.(AX*MVDI).GE.MVDS) THEN
                     AsCLis = ACLI1*EIS*VTAIS*NS1D*(ZZ4*GI2G3*DI0Z2+2.*&
                              ZZ4*SASR1*GIG3*GS2*DI0Z1+ZZ4*SASR2*GI3*  &
                              GS3+ZZ5*GI2G2*GS2*DI0Z2+2.*ZZ5*SASR1*    &
                              GIG2*GS3*DI0Z1+ZZ5*SASR2*GI2*GS4+ZZ6*    &
                              GI2G1*GS3*DI0Z2+2.*ZZ6*SASR1*GIG1*GS4*   &
                              DI0Z1+ZZ6*SASR2*GS5)
                  ELSEIF (MVDI.GT.MVDS.AND.(AX*MVDS).GE.MVDI) THEN
                     AsCLis = ACLI1*EIS*VTAIS*NS1D*(ZZ7*GI2G1*GS3*     &
                              DI0Z2+2.*ZZ7*SASR1*GIG1*GS4*DI0Z1+ZZ7*   &
                              SASR2*GS5+ZZ5*GI2G2*GS2*DI0Z2+2.*ZZ5*    &
                              SASR1*GIG2*GS3*DI0Z1+ZZ5*SASR2*GI2*GS4+  &
                              ZZ8*GI2G3*DI0Z2+2.*ZZ8*SASR1*GIG3*GS2*   &
                              DI0Z1+ZZ8*SASR2*GI3*GS3)
                  ENDIF
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FiCLis = FCLI1*EIS*VTQIS*NS1D*(GIZM2G1*DI0Z2+SASR1*  &
                           2.*GIZMG1*GS2*DI0Z1+SASR2*GIZM1*GS3)
                  ICLis = NCLI1*EIS*VTQIS*NS1D*(GIM2G1*DI0Z2+SASR1*2.* &
                          GIMG1*GS2*DI0Z1+SASR2*GIM1*GS3)
               ENDIF
            ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
               QCLis = QCLI1*EIS*VTQIS*NS1D*(GIM3+SASR1*2.*GIM2*GS2+   &
                       SASR2*GIM1*GS3)
               NCLis = NCLI1*EIS*VTNIS*NS1D*(GI3+SASR1*2.*GI2*GS2+     &
                       SASR2*GS3)
               IF (AI1D.GE.ASMALL) THEN
                  AiCLis = ACLI1*EIS*VTAIS*NS1D*(GI5+SASR1*2.*GI4*GS2+ &
                           SASR2*GI3*GS3)
               ENDIF
               IF (AS1D.GE.ASMALL) THEN
                  IF (MVDS.GE.MVDI.AND.(AX*MVDI).GE.MVDS) THEN
                     AsCLis = ACLI1*EIS*VTAIS*NS1D*(ZZ4*GI5+2.*ZZ4*    &
                              SASR1*GI4*GS2+ZZ4*SASR2*GI3*GS3+ZZ5*GI4* &
                              GS2+2.*ZZ5*SASR1*GI3*GS3+ZZ5*SASR2*GI2*  &
                              GS4+ZZ6*GI3*GS3+2.*ZZ6*SASR1*GI2*GS4+ZZ6*&
                              SASR2*GS5)
                  ELSEIF (MVDI.GT.MVDS.AND.(AX*MVDS).GE.MVDI) THEN
                     AsCLis = ACLI1*EIS*VTAIS*NS1D*(ZZ7*GI3*GS3+2.*ZZ7*&
                              SASR1*GI2*GS4+ZZ7*SASR2*GS5+ZZ5*GI4*GS2+ &
                              2.*ZZ5*SASR1*GI3*GS3+ZZ5*SASR2*GI2*GS4+  &
                              ZZ8*GI5+2.*ZZ8*SASR1*GI4*GS2+ZZ8*SASR2*  &
                              GI3*GS3)
                  ENDIF
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FiCLis = NCLI1*EIS*VTQIS*NS1D*(GIM3+SASR1*2.*GIM2*   &
                           GS2+SASR2*GIM1*GS3)
                  ICLis = NCLI1*EIS*VTQIS*NS1D*(GIM3+SASR1*2.*GIM2*GS2+&
                          SASR2*GIM1*GS3)
               ENDIF
            ENDIF
            QCLis  = MIN(QCLis,QI1D*iDT)
            NCLis  = MIN(NCLis,NI1D*iDT)
            AiCLis = MIN(AiCLis,AI1D*iDT)
            AsCLis = MIN(AsCLis,AI1D*iDT)
            ViCLis = MIN(QCLis*iRHOI,VI1D*iDT)
            FiCLis = MIN(FiCLis,FI1D*iDT)
            ICLis  = MIN(ICLis,I3M1D*iDT)
            IF (ICE_RHOS.EQ.1.OR.AGG_SHAPE.EQ.1) THEN
               IF (ICE_RHOI.EQ.0.AND.ICE_SHAPE.EQ.0) THEN
                  VsCLis = QCLis/RHOS0
                  FsCLis = QCLis*MAX(0.,-1.*LOG(0.8))
               ELSE
                  LIM1 = 0.7*(1.-RHOIS/RHOI0)
                  LIM2 = SQRT(1.-LIM1*LIM1)
                  LIMA = SASR2*LIM2*ISEPR+2.*SASR3*LIM1
                  LIMB = SASR2*LIM1*ISEPR+2.*SASR3*LIM2
                  LIMC = LIM1*(2.*SASR1*LIM2*ISEPR+SASR4*LIM1)
                  LIMD = LIM2*(2.*SASR1*LIM1*ISEPR+SASR4*LIM2)
                  LIME = 2.*(SASR1*LIM1**2.*ISEPR+SASR1*LIM2**2.*ISEPR+&
                         SASR4*LIM1*LIM2)
                  LIMF = LIM1**2.*LIM2*ISEPR
                  LIMG = LIM1*ISEPR*(LIM1**2.+2.*LIM2**2.)
                  LIMH = LIM2*ISEPR*(2.*LIM1**2.+LIM2**2.)
                  LIMI = LIM1*LIM2**2.*ISEPR
                  IF (AGG_SHAPE.EQ.1) THEN
                     DSS0 = SASR4*GS2
                     DSL0 = SASR1*GS2
                     DICC = DI0**(-6.*ZETA)*GI3H1
                     DICA = GI2HG1/DI0Z3
                     DIAC = GIH2G1
                     DIAA = DI0Z3*GI3G1
                     DIA0 = DI0Z1*GIG1
                     DIC0 = GIH1/DI0Z2
                     DIF0 = GIF1/DI0Z1
                     DIC2 = GI2H1/DI0Z4
                     DIA2 = DI0Z2*GI2G1
                  ENDIF
                  IF ((ADAGR-1.).GE.SLIMIT) THEN
                     DSLS  = GS4+LIMA*GS3*DIC0+LIMB*GS3*DIA0+LIMC*GS2* &
                             DIC2+LIMD*GS2*DIA2+LIME*GS2*DIF0+LIMF*    &
                             DICC+LIMG*DICA+LIMH*DIAC+LIMI*DIAA
                     SASP2 = (DSS0+LIM1*ISEPR*DIA0+LIM2*ISEPR*DIC0)/   &
                             (DSL0+LIM1*DIC0+LIM2*DIA0)
                  ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
                     DSLS  = GS4+LIMA*GS3*DIA0+LIMB*GS3*DIC0+LIMC*GS2* &
                             DIA2+LIMD*GS2*DIC2+LIME*GS2*DIF0+LIMF*    &
                             DIAA+LIMG*DIAC+LIMH*DICA+LIMI*DICC
                     SASP2 = (DSS0+LIM1*ISEPR*DIC0+LIM2*ISEPR*DIA0)/   &
                             (DSL0+LIM1*DIA0+LIM2*DIC0)
                  ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
                     DSLS  = GS4+(LIMA+LIMB)*GS3*GI2+(LIMC+LIMD+LIME)* &
                             GS2*GI3+(LIMF+LIMG+LIMH+LIMI)*GI4
                     SASP2 = (DSS0+LIM1*ISEPR*GI2+LIM2*ISEPR*GI2)/     &
                             (DSL0+LIM1*GI2+LIM2*GI2)
                  ENDIF
                  DNSAC  = (RHOI*GI4+RHOS*GS4)/DSLS
                  DNSAC  = MIN(MAX(DNSAC,RHOIMIN),RHOIMAX)
                  SASP2  = MIN(MAX(SASP2,SASMIN),SASMAX)
                  VsCLis = QCLis/DNSAC
                  FsCLis = QCLis*MAX(0.,-1.*LOG(SASP2))
               ENDIF
            ENDIF
         ENDIF
         IF (QG1D.GE.QSMALL.AND.QI1D.GE.QSMALL) THEN
!            EIG   = MIN(1.,0.01*EXP(0.1*TC1D))                          ! FERRIER ET AL., 1995
            EIG   = 0.
            RHOIG = (QI1D+QG1D)/(QI1D/RHOI+QG1D/RHOG+ISMALL)
            RATIO = (RHOIG/RHOG)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            ZZ7   = ZZ1*RATIO**2.-1.
            ZZ8   = ZZ3*RATIO**2.
            VTQ0  = VTQG-VTQI
            VTN0  = VTNG-VTNI
            VTQIG = SQRT(VTQ0*VTQ0+0.04*VTQI*VTQG)
            VTNIG = SQRT(VTN0*VTN0+0.04*VTNI*VTNG)
            IF (AG1D.GE.ASMALL.OR.AI1D.GE.ASMALL) THEN
               VTAX  = VTAG-VTAI
               VTAIG = SQRT(VTAX*VTAX+0.04*VTAI*VTAG)
            ENDIF
            IF ((ADAGR-1.).GE.SLIMIT) THEN
               QCLig = QCLI1*EIG*VTQIG*NG1D*(GIMF1/DI0Z1+GIMG1*GG2*    &
                       DI0Z1+GIMH1*GG2/DI0Z2+GIM1*GG3)
               NCLig = NCLI1*EIG*VTNIG*NG1D*(GIF1/DI0Z1+GIG1*GG2*DI0Z1+&
                       GIH1*GG2/DI0Z2+GG3)
               IF (AI1D.GE.ASMALL) THEN
                  AiCLig = ACLI1*EIG*VTAIG*NG1D*(GIF3/DI0Z1+GIG3*GG2*  &
                           DI0Z1+GIH3*GG2/DI0Z2+GI3*GG3)
               ENDIF
               IF (AG1D.GE.ASMALL) THEN
                  IF (MVDG.GE.MVDI.AND.(AX*MVDI).GE.MVDG) THEN
                     AgCLig = ACLI1*EIG*VTAIG*NG1D*(ZZ4*GIF3/DI0Z1+ZZ4*&
                              GIG3*GG2*DI0Z1+ZZ4*GIH3*GG2/DI0Z2+ZZ4*   &
                              GI3*GG3+ZZ5*GIF2*GG2/DI0Z1+ZZ5*GIG2*GG3* &
                              DI0Z1+ZZ5*GIH2*GG3/DI0Z2+ZZ5*GI2*GG4+ZZ6*&
                              GIF1*GG3/DI0Z1+ZZ6*GIG1*GG4*DI0Z1+ZZ6*   &
                              GIH1*GG4/DI0Z2+ZZ6*GG5)
                  ELSEIF (MVDI.GT.MVDG.AND.(AX*MVDG).GE.MVDI) THEN
                     AgCLig = ACLI1*EIG*VTAIG*NG1D*(ZZ7*GIF1*GG3/DI0Z1+&
                              ZZ7*GIG1*GG4*DI0Z1+ZZ7*GIH1*GG4/DI0Z2+   &
                              ZZ7*GG5+ZZ5*GIF2*GG2/DI0Z1+ZZ5*GIG2*GG3* &
                              DI0Z1+ZZ5*GIH2*GG3/DI0Z2+ZZ5*GI2*GG4+ZZ8*&
                              GIF3/DI0Z1+ZZ8*GIG3*GG2*DI0Z1+ZZ8*GIH3*  &
                              GG2/DI0Z2+ZZ8*GI3*GG3)
                  ENDIF
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FCLig = FCLI1*EIG*VTQIG*NG1D*(GIZMF1/DI0Z1+GIZMG1*   &
                          GG2*DI0Z1+GIZMH1*GG2/DI0Z2+GIZM1*GG3)
                  ICLig = NCLI1*EIG*VTQIG*NG1D*(GIMF1/DI0Z1+GIMG1*GG2* &
                          DI0Z1+GIMH1*GG2/DI0Z2+GIM1*GG3)
               ENDIF
            ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
               QCLig = QCLI1*EIG*VTQIG*NG1D*(GIM2G1*DI0Z2+2.*GIMG1*GG2*&
                       DI0Z1+GIM1*GG3)
               NCLig = NCLI1*EIG*VTNIG*NG1D*(GI2G1*DI0Z2+2.*GIG1*GG2*  &
                       DI0Z1+GG3)
               IF (AI1D.GE.ASMALL) THEN
                  AiCLig = ACLI1*EIG*VTAIG*NG1D*(GI2G3*DI0Z2+2.*GIG3*  &
                           GG2*DI0Z1+GI3*GG3)
               ENDIF
               IF (AG1D.GE.ASMALL) THEN
                  IF (MVDG.GE.MVDI.AND.(AX*MVDI).GE.MVDG) THEN
                     AgCLig = ACLI1*EIG*VTAIG*NG1D*(ZZ4*GI2G3*DI0Z2+2.*&
                              ZZ4*GIG3*GG2*DI0Z1+ZZ4*GI3*GG3+ZZ5*GI2G2*&
                              GG2*DI0Z2+2.*ZZ5*GIG2*GG3*DI0Z1+ZZ5*GI2* &
                              GG4+ZZ6*GI2G1*GG3*DI0Z2+2.*ZZ6*GIG1*GG4* &
                              DI0Z1+ZZ6*GG5)
                  ELSEIF (MVDI.GT.MVDG.AND.(AX*MVDG).GE.MVDI) THEN
                     AgCLig = ACLI1*EIG*VTAIG*NG1D*(ZZ7*GI2G1*GG3*     &
                              DI0Z2+2.*ZZ7*GIG1*GG4*DI0Z1+ZZ7*GG5+ZZ5* &
                              GI2G2*GG2*DI0Z2+2.*ZZ5*GIG2*GG3*DI0Z1+   &
                              ZZ5*GI2*GG4+ZZ8*GI2G3*DI0Z2+2.*ZZ8*GIG3* &
                              GG2*DI0Z1+ZZ8*GI3*GG3)
                  ENDIF
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FCLig = FCLI1*EIG*VTQIG*NG1D*(GIZM2G1*DI0Z2+2.*      &
                          GIZMG1*GG2*DI0Z1+GIZM1*GG3)
                  ICLig = NCLI1*EIG*VTQIG*NG1D*(GIM2G1*DI0Z2+2.*GIMG1* &
                          GG2*DI0Z1+GIM1*GG3)
               ENDIF
            ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
               QCLig = QCLI1*EIG*VTQIG*NG1D*(GIM3+2.*GIM2*GG2+GIM1*GG3)
               NCLig = NCLI1*EIG*VTNIG*NG1D*(GI3+2.*GI2*GG2+GG3)
               IF (AI1D.GE.ASMALL) THEN
                  AiCLig = ACLI1*EIG*VTAIG*NG1D*(GI5+2.*GI4*GG2+GI3*GG3)
               ENDIF
               IF (AG1D.GE.ASMALL) THEN
                  IF (MVDG.GE.MVDI.AND.(AX*MVDI).GE.MVDG) THEN
                     AgCLig = ACLI1*EIG*VTAIG*NG1D*(ZZ4*GI5+2.*ZZ4*GI4*&
                              GG2+ZZ4*GI3*GG3+ZZ5*GI4*GG2+2.*ZZ5*GI3*  &
                              GG3+ZZ5*GI2*GG4+ZZ6*GI3*GG3+2.*ZZ6*GI2*  &
                              GG4+ZZ6*GG5)
                  ELSEIF (MVDI.GT.MVDG.AND.(AX*MVDG).GE.MVDI) THEN
                     AgCLig = ACLI1*EIG*VTAIG*NG1D*(ZZ7*GI3*GG3+2.*ZZ7*&
                              GI2*GG4+ZZ7*GG5+ZZ5*GI4*GG2+2.*ZZ5*GI3*  &
                              GG3+ZZ5*GI2*GG4+ZZ8*GI5+2.*ZZ8*GI4*GG2+  &
                              ZZ8*GI3*GG3)
                  ENDIF
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FCLig = NCLI1*EIG*VTQIG*NG1D*(GIM3+2.*GIM2*GG2+GIM1* &
                          GG3)
                  ICLig = NCLI1*EIG*VTQIG*NG1D*(GIM3+2.*GIM2*GG2+GIM1* &
                          GG3)
               ENDIF
            ENDIF
            QCLig  = MIN(QCLig,QI1D*iDT)
            NCLig  = MIN(NCLig,NI1D*iDT)
            AiCLig = MIN(AiCLig,AI1D*iDT)
            AgCLig = MIN(AgCLig,AI1D*iDT)
            ViCLig = MIN(QCLig*iRHOI,VI1D*iDT)
            FCLig  = MIN(FCLig,FI1D*iDT)
            ICLig  = MIN(ICLig,I3M1D*iDT)
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               LIM1 = 0.7*(1.-RHOIG/RHOI0)
               LIM2 = SQRT(1.-LIM1*LIM1)
               LIMA = LIM2*ISEPR+2.*LIM1
               LIMB = LIM1*ISEPR+2.*LIM2
               LIMC = LIM1*(2.*LIM2*ISEPR+LIM1)
               LIMD = LIM2*(2.*LIM1*ISEPR+LIM2)
               LIME = 2.*(LIM1**2.*ISEPR+LIM2**2.*ISEPR+LIM1*LIM2)
               LIMF = LIM1**2.*LIM2*ISEPR
               LIMG = LIM1*ISEPR*(LIM1**2.+2.*LIM2**2.)
               LIMH = LIM2*ISEPR*(2.*LIM1**2.+LIM2**2.)
               LIMI = LIM1*LIM2**2.*ISEPR
               DICC = DI0**(-6.*ZETA)*GI3H1
               DICA = GI2HG1/DI0Z3
               DIAC = GIH2G1
               DIAA = DI0Z3*GI3G1
               DIA0 = DI0Z1*GIG1
               DIC0 = GIH1/DI0Z2
               DIF0 = GIF1/DI0Z1
               DIC2 = GI2H1/DI0Z4
               DIA2 = DI0Z2*GI2G1
               IF ((ADAGR-1.).GE.SLIMIT) THEN
                  DGLS = GG4+LIMA*GG3*DIC0+LIMB*GG3*DIA0+LIMC*GG2*DIC2+&
                         LIMD*GG2*DIA2+LIME*GG2*DIF0+LIMF*DICC+LIMG*   &
                         DICA+LIMH*DIAC+LIMI*DIAA
               ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
                  DGLS = GG4+LIMA*GG3*DIA0+LIMB*GG3*DIC0+LIMC*GG2*DIA2+&
                         LIMD*GG2*DIC2+LIME*GG2*DIF0+LIMF*DIAA+LIMG*   &
                         DIAC+LIMH*DICA+LIMI*DICC
               ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
                  DGLS = GG4+(LIMA+LIMB)*GG3*GI2+(LIMC+LIMD+LIME)*GG2* &
                         GI3+(LIMF+LIMG+LIMH+LIMI)*GI4
               ENDIF
               DNGAC  = (RHOI*GI4+RHOG*GG4)/DGLS
               DNGAC  = MIN(MAX(DNGAC,RHOIMIN),RHOG0)
               VgCLig = QCLig/DNGAC
            ENDIF
         ENDIF
         IF (QH1D.GE.QSMALL.AND.QI1D.GE.QSMALL) THEN
            RHOIH = (QI1D+QH1D)/(QI1D/RHOI+QH1D/RHOH+ISMALL)
            RATIO = (RHOIH/RHOH)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            ZZ7   = ZZ1*RATIO**2.-1.
            ZZ8   = ZZ3*RATIO**2.
            VTQ0  = VTQH-VTQI
            VTN0  = VTNH-VTNI
            VTQIH = SQRT(VTQ0*VTQ0+0.04*VTQI*VTQH)
            VTNIH = SQRT(VTN0*VTN0+0.04*VTNI*VTNH)
            IF (AH1D.GE.ASMALL.OR.AI1D.GE.ASMALL) THEN
               VTAX  = VTAH-VTAI
               VTAIH = SQRT(VTAX*VTAX+0.04*VTAI*VTAH)
            ENDIF
            IF ((ADAGR-1.).GE.SLIMIT) THEN
               QCLih = QCLI1*EIH*VTQIH*NH1D*(GIMF1/DI0Z1+GIMG1*GH2*    &
                       DI0Z1+GIMH1*GH2/DI0Z2+GIM1*GH3)
               NCLih = NCLI1*EIH*VTNIH*NH1D*(GIF1/DI0Z1+GIG1*GH2*DI0Z1+&
                       GIH1*GH2/DI0Z2+GH3)
               IF (AI1D.GE.ASMALL) THEN
                  AiCLih = ACLI1*EIH*VTAIH*NH1D*(GIF3/DI0Z1+GIG3*GH2*  &
                           DI0Z1+GIH3*GH2/DI0Z2+GI3*GH3)
               ENDIF
               IF (AH1D.GE.ASMALL) THEN
                  IF (MVDH.GE.MVDI.AND.(AX*MVDI).GE.MVDH) THEN
                     AhCLih = ACLI1*EIH*VTAIH*NH1D*(ZZ4*GIF3/DI0Z1+ZZ4*&
                              GIG3*GH2*DI0Z1+ZZ4*GIH3*GH2/DI0Z2+ZZ4*   &
                              GI3*GH3+ZZ5*GIF2*GH2/DI0Z1+ZZ5*GIG2*GH3* &
                              DI0Z1+ZZ5*GIH2*GH3/DI0Z2+ZZ5*GI2*GH4+ZZ6*&
                              GIF1*GH3/DI0Z1+ZZ6*GIG1*GH4*DI0Z1+ZZ6*   &
                              GIH1*GH4/DI0Z2+ZZ6*GH5)
                  ELSEIF (MVDI.GT.MVDH.AND.(AX*MVDH).GE.MVDI) THEN
                     AhCLih = ACLI1*EIH*VTAIH*NH1D*(ZZ7*GIF1*GH3/DI0Z1+&
                              ZZ7*GIG1*GH4*DI0Z1+ZZ7*GIH1*GH4/DI0Z2+   &
                              ZZ7*GH5+ZZ5*GIF2*GH2/DI0Z1+ZZ5*GIG2*GH3* &
                              DI0Z1+ZZ5*GIH2*GH3/DI0Z2+ZZ5*GI2*GH4+ZZ8*&
                              GIF3/DI0Z1+ZZ8*GIG3*GH2*DI0Z1+ZZ8*GIH3*  &
                              GH2/DI0Z2+ZZ8*GI3*GH3)
                  ENDIF
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FCLih = FCLI1*EIH*VTQIH*NH1D*(GIZMF1/DI0Z1+GIZMG1*   &
                          GH2*DI0Z1+GIZMH1*GH2/DI0Z2+GIZM1*GH3)
                  ICLih = NCLI1*EIH*VTQIH*NH1D*(GIMF1/DI0Z1+GIMG1*GH2* &
                          DI0Z1+GIMH1*GH2/DI0Z2+GIM1*GH3)
               ENDIF
            ELSEIF ((1.-ADAGR).GE.SLIMIT) THEN
               QCLih = QCLI1*EIH*VTQIH*NH1D*(GIM2G1*DI0Z2+2.*GIMG1*GH2*&
                       DI0Z1+GIM1*GH3)
               NCLih = NCLI1*EIH*VTNIH*NH1D*(GI2G1*DI0Z2+2.*GIG1*GH2*  &
                       DI0Z1+GH3)
               IF (AI1D.GE.ASMALL) THEN
                  AiCLih = ACLI1*EIH*VTAIH*NH1D*(GI2G3*DI0Z2+2.*GIG3*  &
                           GH2*DI0Z1+GI3*GH3)
               ENDIF
               IF (AH1D.GE.ASMALL) THEN
                  IF (MVDH.GE.MVDI.AND.(AX*MVDI).GE.MVDH) THEN
                     AhCLih = ACLI1*EIH*VTAIH*NH1D*(ZZ4*GI2G3*DI0Z2+2.*&
                              ZZ4*GIG3*GH2*DI0Z1+ZZ4*GI3*GH3+ZZ5*GI2G2*&
                              GH2*DI0Z2+2.*ZZ5*GIG2*GH3*DI0Z1+ZZ5*GI2* &
                              GH4+ZZ6*GI2G1*GH3*DI0Z2+2.*ZZ6*GIG1*GH4* &
                              DI0Z1+ZZ6*GH5)
                  ELSEIF (MVDI.GT.MVDH.AND.(AX*MVDH).GE.MVDI) THEN
                     AhCLih = ACLI1*EIH*VTAIH*NH1D*(ZZ7*GI2G1*GH3*     &
                              DI0Z2+2.*ZZ7*GIG1*GH4*DI0Z1+ZZ7*GH5+ZZ5* &
                              GI2G2*GH2*DI0Z2+2.*ZZ5*GIG2*GH3*DI0Z1+   &
                              ZZ5*GI2*GH4+ZZ8*GI2G3*DI0Z2+2.*ZZ8*GIG3* &
                              GH2*DI0Z1+ZZ8*GI3*GH3)
                  ENDIF
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FCLih = FCLI1*EIH*VTQIH*NH1D*(GIZM2G1*DI0Z2+2.*      &
                          GIZMG1*GH2*DI0Z1+GIZM1*GH3)
                  ICLih = NCLI1*EIH*VTQIH*NH1D*(GIM2G1*DI0Z2+2.*GIMG1* &
                          GH2*DI0Z1+GIM1*GH3)
               ENDIF
            ELSEIF (ABS(ADAGR-1.).LT.SLIMIT) THEN
               QCLih = QCLI1*EIH*VTQIH*NH1D*(GIM3+2.*GIM2*GH2+GIM1*GH3)
               NCLih = NCLI1*EIH*VTNIH*NH1D*(GI3+2.*GI2*GH2+GH3)
               IF (AI1D.GE.ASMALL) THEN
                  AiCLih = ACLI1*EIH*VTAIH*NH1D*(GI5+2.*GI4*GH2+GI3*GH3)
               ENDIF
               IF (AH1D.GE.ASMALL) THEN
                  IF (MVDH.GE.MVDI.AND.(AX*MVDI).GE.MVDH) THEN
                     AhCLih = ACLI1*EIH*VTAIH*NH1D*(ZZ4*GI5+2.*ZZ4*GI4*&
                              GH2+ZZ4*GI3*GH3+ZZ5*GI4*GH2+2.*ZZ5*GI3*  &
                              GH3+ZZ5*GI2*GH4+ZZ6*GI3*GH3+2.*ZZ6*GI2*  &
                              GH4+ZZ6*GH5)
                  ELSEIF (MVDI.GT.MVDH.AND.(AX*MVDH).GE.MVDI) THEN
                     AhCLih = ACLI1*EIH*VTAIH*NH1D*(ZZ7*GI3*GH3+2.*ZZ7*&
                              GI2*GH4+ZZ7*GH5+ZZ5*GI4*GH2+2.*ZZ5*GI3*  &
                              GH3+ZZ5*GI2*GH4+ZZ8*GI5+2.*ZZ8*GI4*GH2+  &
                              ZZ8*GI3*GH3)
                  ENDIF
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FCLih = NCLI1*EIH*VTQIH*NH1D*(GIM3+2.*GIM2*GH2+GIM1* &
                          GH3)
                  ICLih = NCLI1*EIH*VTQIH*NH1D*(GIM3+2.*GIM2*GH2+GIM1* &
                          GH3)
               ENDIF
            ENDIF
            QCLih  = MIN(QCLih,QI1D*iDT)
            NCLih  = MIN(NCLih,NI1D*iDT)
            VCLih  = MIN(QCLih*iRHOI,VI1D*iDT)
            AiCLih = MIN(AiCLih,AI1D*iDT)
            AhCLih = MIN(AhCLih,AI1D*iDT)
            FCLih  = MIN(FCLih,FI1D*iDT)
            ICLih  = MIN(ICLih,I3M1D*iDT)
         ENDIF
         IF (QS1D.GE.QSMALL.AND.QG1D.GE.QSMALL) THEN
            ESG   = MIN(1.,0.01*EXP(0.1*TC1D))                         ! FERRIER ET AL., 1995
            RHOSG = (QS1D+QG1D)/(QS1D/RHOS+QG1D/RHOG+ISMALL)
            RATIO = (RHOSG/RHOG)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            ZZ7   = ZZ1*RATIO**2.-1.
            ZZ8   = ZZ3*RATIO**2.
            VTQ0  = VTQG-VTQS
            VTN0  = VTNG-VTNS
            VTQSG = SQRT(VTQ0*VTQ0+0.04*VTQS*VTQG)
            VTNSG = SQRT(VTN0*VTN0+0.04*VTNS*VTNG)
            IF (AS1D.GE.ASMALL.OR.AG1D.GE.ASMALL) THEN
               VTAX  = VTAG-VTAS
               VTASG = SQRT(VTAX*VTAX+0.04*VTAS*VTAG)
            ENDIF
            QCLsg = QCLS1*ESG*VTQSG*NG1D*(SASR2*GSM3+SASR1*2.*GG2*GSM2+&
                    GSM1*GG3)
            NCLsg = NCLS1*ESG*VTNSG*NG1D*(SASR2*GS3+SASR1*2.*GS2*GG2+  &
                    GG3)
            QCLsg = MIN(QCLsg,QS1D*iDT)
            NCLsg = MIN(NCLsg,NS1D*iDT)
            IF (ICE_RHOS.EQ.1) THEN
               VsCLsg = MIN(QCLsg*iRHOS,VS1D*iDT)
            ENDIF
            IF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
               LIM1   = 0.7*(1.-RHOSG/RHOI0)
               LIM2   = SQRT(1.-LIM1*LIM1)
               LIMA   = SASR4*(LIM1*ISEPR+2.*LIM2)+LIM2*ISEPR*SASR1+   &
                        LIM1*(SASR1+SASR3)
               LIMB   = LIM1**2.+2.*ISEPR*LIM1*LIM2*(SASR2+SASR4**2.)
               LIMC   = LIM2**2.*SASR4**2.+2.*ISEPR*SASR3*(LIM1**2.+   &
                        LIM2**2.)+LIM1*LIM2*(SASPR+SASR3)
               LIMD   = ISEPR*(LIM1**3.+LIM2**3.*SASPR+2.*LIM1**2.*    &
                        LIM2*SASPR+LIM1**2.*LIM2/SASPR+2.*LIM1*        &
                        LIM2**2.+LIM1*LIM2**2.*SASPR**2.)
               DGLS   = GG4+LIMA*GG3*GS2+(LIMB+LIMC)*GG2*GS3+LIMD*GS4
               DNGAC  = MAX(MIN((RHOS*GS4+RHOG*GG4)/DGLS,RHOG0),RHOIMIN)
               VgCLsg = MIN(QCLsg/DNGAC,VS1D*iDT)
            ENDIF
            IF (AGG_SHAPE.EQ.1) THEN
               FCLsg = MIN(QCLsg*MAX(0.,-1.*LOG(SASPR)),FS1D*iDT)
            ENDIF
            IF (AS1D.GE.ASMALL) THEN
               AsCLsg = ACLS1*ESG*VTASG*NG1D*(SASR2*GS5+SASR1*2.*GG2*  &
                        GS4+GG3*GS3)
               AsCLsg = MIN(AsCLsg,AS1D*iDT)
            ENDIF
            IF (AG1D.GE.ASMALL) THEN
               IF (MVDG.GE.MVDS.AND.(AX*MVDS).GE.MVDG) THEN
                  AgCLsg = ACLS1*ESG*VTASG*NG1D*(ZZ4*SASR2*GS5+2.*ZZ4* &
                           SASR1*GS4*GG2+ZZ4*GG3*GS3+ZZ5*SASR2*GS4*GG2+&
                           2.*ZZ5*SASR1*GS3*GG3+ZZ5*GS2*GG4+ZZ6*SASR2* &
                           GS3*GG3+2.*ZZ6*SASR1*GS2*GG4+ZZ6*GG5)
               ELSEIF (MVDS.GT.MVDG.AND.(AX*MVDG).GE.MVDS) THEN
                  AgCLsg = ACLS1*ESG*VTASG*NG1D*(ZZ7*SASR2*GS3*GG3+2.* &
                           ZZ7*SASR1*GS2*GG4+ZZ7*GG5+ZZ5*SASR2*GS4*GG2+&
                           2.*ZZ5*SASR1*GS3*GG3+ZZ5*GS2*GG4+ZZ8*SASR2* &
                           GS5+2.*ZZ8*SASR1*GS4*GG2+ZZ8*GS3*GG3)
               ENDIF
               AgCLsg = MIN(AgCLsg,AS1D*iDT)
            ENDIF
         ENDIF
         IF (QH1D.GE.QSMALL.AND.QS1D.GE.QSMALL) THEN
            RHOSH = (QS1D+QH1D)/(QS1D/RHOS+QH1D/RHOH+ISMALL)
            RATIO = (RHOSH/RHOH)**THRD
            ZZ4   = ZZ1*RATIO**2.
            ZZ5   = ZZ2*RATIO**2.
            ZZ6   = ZZ3*RATIO**2.-1.
            ZZ7   = ZZ1*RATIO**2.-1.
            ZZ8   = ZZ3*RATIO**2.
            VTQ0  = VTQH-VTQS
            VTN0  = VTNH-VTNS
            VTQSH = SQRT(VTQ0*VTQ0+0.04*VTQS*VTQH)
            VTNSH = SQRT(VTN0*VTN0+0.04*VTNS*VTNH)
            IF (AH1D.GE.ASMALL.OR.AS1D.GE.ASMALL) THEN
               VTAX  = VTAH-VTAS
               VTASH = SQRT(VTAX*VTAX+0.04*VTAS*VTAH)
            ENDIF
            QCLsh = QCLS1*ESH*VTQSH*NH1D*(SASR2*GSM3+SASR1*2.*GH2*GSM2+&
                    GH3*GSM1)
            NCLsh = NCLS1*ESH*VTNSH*NH1D*(SASR2*GS3*+SASR1*2.*GH2*GS2+ &
                    GH3)
            QCLsh = MIN(QCLsh,QS1D*iDT)
            NCLsh = MIN(NCLsh,NS1D*iDT)
            IF (ICE_RHOS.EQ.1) THEN
               VCLsh = MIN(QCLsh*iRHOS,VS1D*iDT)
            ENDIF
            IF (AGG_SHAPE.EQ.1) THEN
               FCLsh = QCLsh*MAX(0.,-1.*LOG(SASPR))
            ENDIF
            IF (AS1D.GE.ASMALL) THEN
               AsCLsh = ACLS1*ESH*VTASH*NH1D*(SASR2*GS5+SASR1*2.*GH2*  &
                        GS4+GH3*GS3)
               AsCLsh = MIN(AsCLsh,AS1D*iDT)
            ENDIF
            IF (AH1D.GE.ASMALL) THEN
               IF (MVDH.GE.MVDS.AND.(AX*MVDS).GE.MVDH) THEN
                  AhCLsh = ACLS1*ESH*VTASH*NH1D*(ZZ4*SASR2*GS5+2.*ZZ4* &
                           SASR1*GS4*GH2+ZZ4*GH3*GS3+ZZ5*SASR2*GS4*GH2+&
                           2.*ZZ5*SASR1*GS3*GH3+ZZ5*GS2*GH4+ZZ6*SASR2* &
                           GS3*GH3+2.*ZZ6*SASR1*GS2*GH4+ZZ6*GH5)
               ELSEIF (MVDS.GT.MVDH.AND.(AX*MVDH).GE.MVDS) THEN
                  AhCLsh = ACLS1*ESH*VTASH*NH1D*(ZZ7*SASR2*GS3*GH3+2.* &
                           ZZ7*SASR1*GS2*GH4+ZZ7*GH5+ZZ5*SASR2*GS4*GH2+&
                           2.*ZZ5*SASR1*GS3*GH3+ZZ5*GS2*GH4+ZZ8*SASR2* &
                           GS5+2.*ZZ8*SASR1*GS4*GH2+ZZ8*GS3*GH3)
               ENDIF
               AhCLsh = MIN(AhCLsh,AS1D*iDT)
            ENDIF
         ENDIF
         IF (QS1D.GE.QSMALL.AND.MVDS.GE.2.E-4) THEN
            ESS1 = 0.1*MIN(1.,0.05*EXP(0.1*TC1D))                      ! FERRIER ET AL., 1995
            IF (ICE_RHOS.NE.0) THEN
               ESS2 = 0.001*(1.-RHOS/RHOI0)
            ELSE
               ESS2 = 0.
            ENDIF
!            ESS   = MIN(MAX(ESS1,ESS2,0.),1.)
            ESS   = ESS1
            QCLss = 2.*PI*XISP*VTQS*ESS*SASR2/3.*NS1D**2.*GSM3*AMS
            QCLss = MIN(QCLss,0.5*QS1D*iDT)
            NCLss = 2.*PI*XISP*VTNS*ESS*SASR2/3.*NS1D**2.*GS3
            IF (ICE_RHOS.EQ.1) THEN
               IF (ICE_RHOI.EQ.0.AND.ICE_SHAPE.EQ.0) THEN
                  VCLss = MIN(QCLss/RHOS0,0.5*VS1D*iDT)
                  NCLss = NCLss-QCLss*V2M3/(RHOS0*(1.26*MVDS)**3.)
               ELSE
                  LIM1  = 0.7*(1.-RHOS/RHOI0)
                  LIM2  = SQRT(1.-LIM1*LIM1)
                  LIMA  = 1.+LIM1/SASPR+LIM2+2.*LIM1*ISEPR*SASPR+2.*   &
                          LIM2*ISEPR
                  LIMB  = 2.*ISEPR*(LIM1**2.+LIM2**2.+LIM1*LIM2*SASPR+ &
                          LIM1*LIM2/SASPR)
                  LIMC  = ISEPR**2.*(LIM1**2.*SASPR**2.+LIM2**2.+2.*   &
                          LIM1*LIM2*SASPR+LIM1**3.*SASPR+LIM2**3.)
                  LIMD  = LIM1*LIM2*ISEPR**2.*(2.*LIM1+LIM1*SASPR**2.+ &
                          LIM2/SASPR+2.*LIM2*SASPR)
                  DSLS  = LIMA+LIMB+LIMC+LIMD
                  DNSAG = MAX(MIN(2.*RHOS/DSLS,RHOS0),RHOIMIN)
                  VCLss = MIN(QCLss/DNSAG,0.5*VS1D*iDT)
                  RATIO = (RHOS/DNSAG)**THRD
                  NCLss = NCLss-QCLss*V2M3/(DNSAG*(RATIO*1.26*MVDS)**3.)
               ENDIF
            ELSE
               NCLss = NCLss-QCLss*V2M3/(RHOS*(1.26*MVDS)**3.)
            ENDIF
            NCLss = MAX(0.,MIN(NCLss,0.5*NS1D*iDT))
            IF (AGG_SHAPE.EQ.1) THEN
               IF (ICE_RHOI.EQ.0.AND.ICE_SHAPE.EQ.0) THEN
                  FCLss = MIN(QCLss*MAX(0.,-1.*LOG(0.8)),0.5*FS1D*iDT)
               ELSE
                  LIM1  = 0.7*(1.-RHOS/RHOS0)
                  LIM2  = SQRT(1.-LIM1*LIM1)
                  SASP3 = MIN(SASR4*(1.+LIM2)+LIM1*SASR1,SASR1*(1.+    &
                          LIM2*ISEPR)+LIM1*ISEPR*SASR4)/MAX(SASR4*(1.+ &
                          LIM2)+LIM1*SASR1,SASR1*(1.+LIM2*ISEPR)+LIM1* &
                          ISEPR*SASR4)
                  SASP3 = MIN(MAX(SASP3,SASMIN),SASMAX)
                  FCLss = MIN(QCLss*MAX(0.,-1.*LOG(SASP3)),0.5*FS1D*iDT)
               ENDIF
            ENDIF
            IF (AS1D.GE.ASMALL) THEN
               ACLss = 2.*PI*XISP*VTAS*ESS*SASR2/3.*NS1D**2.*GS5*0.4126
               ACLss = MIN(ACLss,0.5*AS1D*iDT)
            ENDIF
         ENDIF
      ELSE                                                              ! WARM CLOUD START
         IF ((HCwqv+HRwqv+HSwqv+HGwqv+HHwqv).GE.QSMALL) THEN
            VDMAX = XXLV*(QV1D-QVSW0)/(1.+XXLV**2.*QV1D/(CPM*RV*TK1D** &
                    2.))*iDT
            SUMCND = HCwqv+HRwqv+HSwqv+HGwqv+HHwqv
            IF (SUMCND.GT.VDMAX.AND.VDMAX.GE.QSMALL) THEN
               RATIO = MIN(1.,VDMAX/(SUMCND+QSMALL))
               HCwqv = HCwqv*RATIO
               HRwqv = HRwqv*RATIO
               HSwqv = HSwqv*RATIO
               HGwqv = HGwqv*RATIO
               HHwqv = HHwqv*RATIO
            ENDIF
         ENDIF
         IF ((HCwqv+HRwqv+HSwqv+HGwqv+HHwqv).LT.0.) THEN
            EVMAX = XXLV*(QV1D-QVSW0)/(1.+XXLV**2.*QV1D/(CPM*RV*TK1D** &
                    2.))*iDT
            SUMEVP = HCwqv+HRwqv+HSwqv+HGwqv+HHwqv
            IF (EVMAX.LT.0..AND.SUMEVP.LT.EVMAX*0.9999) THEN
               HCwqv = HCwqv*MIN(1.,0.9999*EVMAX/SUMEVP)
               HRwqv = HRwqv*MIN(1.,0.9999*EVMAX/SUMEVP)
               HSwqv = HSwqv*MIN(1.,0.9999*EVMAX/SUMEVP)
               HGwqv = HGwqv*MIN(1.,0.9999*EVMAX/SUMEVP)
               HHwqv = HHwqv*MIN(1.,0.9999*EVMAX/SUMEVP)
            ENDIF
         ENDIF
!----------------- MELTING PROXESSES -----------------------------------
         IF (QI1D.GE.QSMALL) THEN
            IF (MVDI.GT.DCR) THEN
               QMLir = -1.*QI1D*iDT
               VMLir = -1.*VI1D*iDT
               IF (AI1D.GE.ASMALL) THEN
                  AMLir = -1.*AI1D*iDT
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FMLir = -1.*FI1D*iDT
                  IMLir = -1.*I3M1D*iDT
               ENDIF
            ELSE
               QMLic = -1.*QI1D*iDT
               VMLic = -1.*VI1D*iDT
               IF (AI1D.GE.ASMALL) THEN
                  AMLic = -1.*AI1D*iDT
               ENDIF
               IF (I3M1D.GE.ISMALL.AND.FI1D.GE.ISMALL) THEN
                  FMLic = -1.*FI1D*iDT
                  IMLic = -1.*I3M1D*iDT
               ENDIF
            ENDIF
            NMLir = MIN(QMLir*NI1D/QI1D,0.)
            NMLic = MIN(QMLic*NI1D/QI1D,0.)
         ENDIF
         IF (QS1D.GE.QSMALL) THEN
            RMcsq = 0.
            RMrsq = 0.
            RMcsa = 0.
            RMrsa = 0.
            ZZ6   = ZZ3-1.
            IF (QC1D.GE.QSMALL) THEN
               ECS   = 1.
               VTQ0  = VTQS-VTQC
               VTQSC = SQRT(VTQ0*VTQ0+4.E-2*VTQS*VTQC)
               RMcsq = QRMC1*ECS*VTQSC*NS1D*(GS3*GC4+2.*GS2*GC5+GC6)
               RMcsq = MIN(RMcsq,QC1D*iDT)
               IF (AS1D.GE.ASMALL.AND.(AX*MVDC).GE.MVDS) THEN
                  RMcsa = ARMC1*ECS*VTQSC*NS1D*(ZZ1*GS3*GC3+2.*ZZ1*GS2*&
                          GC4+ZZ1*GC5+ZZ2*GS4*GC2+2.*ZZ2*GS3*GC3+ZZ2*  &
                          GS2*GC4+ZZ6*GS5+2.*ZZ6*GS4*GC2+ZZ6*GS3*GC3)
                  RMcsa = MIN(RMcsa,QC1D*iDT*iAPW/MVDC)
               ENDIF
            ENDIF
            IF (QR1D.GE.QSMALL) THEN
               ERS   = 1.
               VTQ0  = VTQR-VTQS
               VTQRS = SQRT(VTQ0*VTQ0+4.E-2*VTQR*VTQS)
               RMrsq = QRMR1*ERS*VTQRS*NS1D*(GS3*GR4+2.*GS2*GR5+GR6)
               RMrsq = MIN(RMrsq,QR1D*iDT)
               IF (AS1D.GE.ASMALL.AND.(AX*MVDR).GE.MVDS) THEN
                  RMrsa = ARMR1*ERS*VTQRS*NS1D*(ZZ1*GS3*GR3+2.*ZZ1*GS2*&
                          GR4+ZZ1*GR5+ZZ2*GS4*GR2+2.*ZZ2*GS3*GR3+ZZ2*  &
                          GS2*GR4+ZZ6*GS5+2.*ZZ6*GS4*GR2+ZZ6*GS3*GR3)
                  RMrsa = MIN(RMrsa,QR1D*iDT*iAPW/MVDR)
               ENDIF
            ENDIF
            MLMAX = 1.-MIN(1.,RHOS/RHOIMAX)
            HSwcd = 2.*PI*NS1D*KAP*VENQS*(TK0C-TK1D)
            HSwrm = CPW*(TK0C-TK1D)*(RMrsq+RMcsq)
            SMLF  = MAX(0.,-1.*(HSwcd+HSwqv+HSwrm)*DT/XXLF0)/QS1D
            SMLF  = MIN(MLMAX,MAX(0.,SMLF))
            SMLR  = 0.01195*EXP(4.411*SMLF)
            AVWSG = SMLR*AVRH+(1.-SMLR)*AVSG
            BVWSG = SMLR*BVRH+(1.-SMLR)*BVSG
            RHOWS = SMLF*RHOW+(1.-SMLF)*RHOS
            IF (ICE_VTS.EQ.1.OR.ICE_VTS.EQ.2) THEN
               IF (SMLF.GE.1.E-8) THEN
                  AVS = 88.701591*EXP(SMLF*2.495827)
                  BVS = -0.182754*SMLF**2.+0.368609*SMLF+0.608571
               ELSE
                  AVS = AVS
                  BVS = BVS
               ENDIF
            ENDIF
            BSTMP = SCN*SQRT(AVS*RHOAJ/MUA)
            QTMP1 = LLMS*(1.5+BVS/2.)
            QTMP3 = EXP(GAMLN(BVS/2.+AFAS+2.5)-GAMLN(AFAS+1.)-QTMP1)
            IF (AGG_SHAPE.EQ.0) THEN
               WSAPR = 1.
            ELSEIF (AGG_SHAPE.EQ.1) THEN
               WSAPR = EXP((1.-SMLF)*LOG(SASPR))
            ENDIF
            CAPWS = ZP1*WSAPR**(ZP2/3.)+ZP3*WSAPR**(ZP4/3.)
            VENWS = AVWSG*QTMP2*CAPWS+BVWSG*BSTMP*QTMP3*CAPWS
            SMLTQ = CPW/XXLF*(TK0C-TK1D)*(RMcsq+RMrsq)
            QMLsr = 2.*PI*NS1D*KAP*(TK0C-TK1D)/XXLF*VENWS+SMLTQ     ! NEGATIVE
            QMLsr = MAX(MIN(QMLsr,0.),-1.*QS1D*iDT)
            NMLsr = MIN(QMLsr*NS1D/QS1D,0.)
            IF (ICE_RHOS.EQ.1) THEN
               VMLsr = ((QS1D+QMLsr*DT)/RHOWS-VS1D)*iDT
               VMLsr = MAX(VMLsr,-1.*VS1D*iDT)
            ENDIF
            IF (AGG_SHAPE.EQ.1) THEN
               FMLsr = (-1.*FS1D-LOG(WSAPR)*(QS1D+QMLsr*DT))*iDT
               FMLsr = MAX(FMLsr,-1.*FS1D*iDT)
            ENDIF
            IF (AS1D.GE.ASMALL) THEN
               ATMP1 = LOG(LAMS)*(0.5+BVS/2.)
               ATMP2 = EXP(GAMLN(BVS/2.+AFAS+1.5)-GAMLN(AFAS+1.)-ATMP1)
               VENAS = AVWSG*CAPWS+BVWSG*BSTMP*ATMP2*CAPWS
               SMLTA = CPW/XXLF*(TK0C-TK1D)*(RMcsa+RMrsa)
               AMLsr = 8.*NS1D*KAP*(TK0C-TK1D)/XXLF/RHOWS*VENAS+SMLTA
               AMLsr = MAX(MIN(AMLsr,0.),-1.*AS1D*iDT)
            ENDIF
         ENDIF
         IF (QG1D.GE.QSMALL) THEN
            RMcgq = 0.
            RMrgq = 0.
            RMcga = 0.
            RMrga = 0.
            ZZ6   = ZZ3-1.
            IF (QC1D.GE.QSMALL) THEN
               ECG   = 1.
               VTQ0  = VTQG-VTQC
               VTQGC = SQRT(VTQ0*VTQ0+0.04*VTQG*VTQC)
               RMcgq = QRMC1*ECG*VTQGC*NG1D*(GG3*GC4+2.*GG2*GC5+GC6)
               RMcgq = MIN(RMcgq,QC1D*iDT)
               IF (AG1D.GE.ASMALL.AND.(AX*MVDC).GE.MVDG) THEN
                  RMcga = ARMC1*ECG*VTQGC*NG1D*(ZZ1*GG3*GC3+2.*ZZ1*GG2*&
                          GC4+ZZ1*GC5+ZZ2*GG4*GC2+2.*ZZ2*GG3*GC3+ZZ2*  &
                          GG2*GC4+ZZ6*GG5+2.*ZZ6*GG4*GC2+ZZ6*GG3*GC3)
                  RMcga = MIN(RMcga,QC1D*iDT*iAPW/MVDC)
               ENDIF
            ENDIF
            IF (QR1D.GE.QSMALL) THEN
               ERG  = 1.
               VTQ0  = VTQR-VTQG
               VTQRG = SQRT(VTQ0*VTQ0+0.04*VTQR*VTQG)
               RMrgq = QRMR1*ERG*VTQRG*NG1D*(GG3*GR4+2.*GG2*GR5+GR6)
               RMrgq = MIN(RMrgq,QR1D*iDT)
               IF (AG1D.GE.ASMALL.AND.(AX*MVDR).GE.MVDG) THEN
                  RMrga = ARMR1*ERG*VTQRG*NG1D*(ZZ1*GG3*GR3+2.*ZZ1*GG2*&
                          GR4+ZZ1*GR5+ZZ2*GG4*GR2+2.*ZZ2*GG3*GR3+ZZ2*  &
                          GG2*GR4+ZZ6*GG5+2.*ZZ6*GG4*GR2+ZZ6*GG3*GR3)
                  RMrga = MIN(RMrga,QR1D*iDT*iAPW/MVDR)
               ENDIF
            ENDIF
            MLMAX = 1.-MIN(1.,RHOG/RHOIMAX)
            HGwcd = 2.*PI*NG1D*KAP*VENQG*(TK0C-TK1D)
            HGwrm = CPW*(TK0C-TK1D)*(RMrgq+RMcgq)
            GMLF  = MAX(0.,-1.*(HGwcd+HGwqv+HGwrm)*DT/XXLF0)/QG1D
            GMLF  = MIN(MLMAX,MAX(0.,GMLF))
            GMLR  = 0.01195*EXP(4.411*GMLF)
            AVWSG = GMLR*AVRH+(1.-GMLR)*AVSG
            BVWSG = GMLR*BVRH+(1.-GMLR)*BVSG
            RHOWG = GMLF*RHOW+(1.-GMLF)*RHOG
            IF (ICE_VTG.EQ.1) THEN
               KINV  = (1.72E-5*(393./(TK1D+120.))*(TK1D/TK0C)**1.5)/RHO
               BEST0 = 2.*G*NG1D/(KINV**2.)
               BEST  = 2.*BEST0*RHOWG*EXP(GAMLN(BMG+AFAG+1.)-GAMLN(    &
                       AFAG+1.)-BMG*LOG(LAMG))/3.
               C1X2  = VTC1*BEST**5.E-1
               VTB1  = C1X2/(1.+C1X2)**5.E-1/((1.+C1X2)**5.E-1-1.)/2.
               VTA1  = VTC2*((1+C1X2)**5.E-1-1.)**2./BEST**VTB1
               AVG   = VTA1*KINV**(1.-2.*VTB1)*(4.*RHOWG*G/3.)**VTB1
               BVG   = VTB1*BMG-1.
            ENDIF
            BGTMP = SCN*SQRT(AVG*RHOAJ/MUA)
            QTMP1 = LLMG*(1.5+BVG/2.)
            QTMP3 = EXP(GAMLN(BVG/2.+AFAG+2.5)-GAMLN(AFAG+1.)-QTMP1)
            VENWG = AVWSG*QTMP2+BVWSG*BGTMP*QTMP3
            GMLTQ = CPW/XXLF*(TK0C-TK1D)*(RMcgq+RMrgq)
            QMLgr = 2.*PI*NG1D*KAP*(TK0C-TK1D)/XXLF*VENWG+GMLTQ     ! NEGATIVE
            QMLgr = MAX(MIN(QMLgr,0.),-1.*QG1D*iDT)
            NMLgr = MIN(QMLgr*NG1D/QG1D,0.)
            VMLgr = MAX(((QG1D+QMLgr*DT)/RHOWG-VG1D)*iDT,-1.*VG1D*iDT)
            IF (AG1D.GE.ASMALL) THEN
               ATMP1 = LOG(LAMG)*(0.5+BVG/2.)
               ATMP2 = EXP(GAMLN(BVG/2.+AFAG+1.5)-GAMLN(AFAG+1.)-ATMP1)
               VENAG = AVWSG+BVWSG*BGTMP*ATMP2
               GMLTA = CPW/XXLF*(TK0C-TK1D)*(RMcga+RMrga)
               AMLgr = 8.*NG1D*KAP*(TK0C-TK1D)/XXLF/RHOWG*VENAG+GMLTA ! NEGATIVE
               AMLgr = MAX(MIN(AMLgr,0.),-1.*AG1D*iDT)
            ENDIF
         ENDIF
         IF (QH1D.GE.QSMALL) THEN
            RMchq = 0.
            RMrhq = 0.
            RMcha = 0.
            RMrha = 0.
            ZZ6   = ZZ3-1.
            IF (QC1D.GE.QSMALL) THEN
               ECH   = 1.
               VTQ0  = VTQH-VTQC
               VTQHC = SQRT(VTQ0*VTQ0+0.04*VTQH*VTQC)
               RMchq = QRMC1*ECH*VTQHC*NH1D*(GH3*GC4+2.*GH2*GC5+GC6)
               RMchq = MIN(RMchq,QC1D*iDT)
               IF (AH1D.GE.ASMALL.AND.(AX*MVDC).GE.MVDH) THEN
                  RMcha = ARMC1*ECH*VTQHC*NH1D*(ZZ1*GH3*GC3+2.*ZZ1*GH2*&
                          GC4+ZZ1*GC5+ZZ2*GH4*GC2+2.*ZZ2*GH3*GC3+ZZ2*  &
                          GH2*GC4+ZZ6*GH5+2.*ZZ6*GH4*GC2+ZZ6*GH3*GC3)
                  RMcha = MIN(RMcha,QC1D*iDT*iAPW/MVDC)
               ENDIF
            ENDIF
            IF (QR1D.GE.QSMALL) THEN
               ERH   = 1.
               VTQ0  = VTQH-VTQR
               VTQRH = SQRT(VTQ0*VTQ0+0.04*VTQH*VTQR)
               RMrhq = QRMR1*ERH*VTQRH*NH1D*(GH3*GR4+2.*GH2*GR5+GR6)
               RMrhq = MIN(RMrhq,QR1D*iDT)
               IF (AH1D.GE.ASMALL.AND.(AX*MVDR).GE.MVDH) THEN
                  RMrha = ARMR1*ERH*VTQRH*NH1D*(ZZ1*GH3*GR3+2.*ZZ1*GH2*&
                          GR4+ZZ1*GR5+ZZ2*GH4*GR2+2.*ZZ2*GH3*GR3+ZZ2*  &
                          GH2*GR4+ZZ6*GH5+2.*ZZ6*GH4*GR2+ZZ6*GH3*GR3)
                  RMrha = MIN(RMrha,QR1D*iDT*iAPW/MVDR)
               ENDIF
            ENDIF
            HMLTQ = CPW/XXLF*(TK0C-TK1D)*(RMrhq+RMchq)
            QMLhr = 2.*PI*NH1D*KAP*(TK0C-TK1D)/XXLF*VENQH+HMLTQ    ! NEGATIVE
            QMLhr = MAX(MIN(QMLhr,0.),-1.*QH1D*iDT)
            NMLhr = MIN(QMLhr*NH1D/QH1D,0.)
            IF (AH1D.GE.ASMALL) THEN
               HMLTA = CPW/XXLF*(TK0C-TK1D)*(RMrha+RMcha)
               AMLhr = 8.*NH1D*KAP*(TK0C-TK1D)/XXLF/RHOH*VENAH+HMLTA ! NEGATIVE
               AMLhr = MAX(MIN(AMLhr,0.),-1.*AH1D*iDT)
            ENDIF
         ENDIF
      ENDIF                                                             ! TEMPERATURE LOOPS
!----------------- HAILSTONE GROWTH HEATBALANCE ------------------------
      IF (HWET_MODE.EQ.1) THEN
         DMWDT = QRMch+QRMrh
         DMIDT = QCLih+QCLsh
         IF (QH1D.GE.QLIMIT.AND.DMWDT.GE.QSMALL) THEN
            HHdcd = 2.*PI*NH1D*KAP*VENQH*(TK1D-TK0C)
            HHdrm = DMWDT*(XXLF+CPW*(TK1D-TK0C))
            HHdcl = DMIDT*CPI*(TK1D-TK0C)
            HHdtt = HHdcd+HHdqv+HHdrm+HHdcl
            IF (HHdtt.GT.0.) THEN                                       ! WET-GROWTH MODE
               QCLih = QCLih/EIH; NCLih = NCLih/EIH; VCLih = VCLih/EIH
               FCLih = FCLih/EIH; ICLih = ICLih/EIH
               QCLsh = QCLsh/ESH; NCLsh = NCLsh/ESH; VCLsh = VCLsh/ESH
               AiCLih = AiCLih/EIH; AhCLih = AhCLih/EIH
               AsCLsh = AsCLsh/ESH; AhCLsh = AhCLsh/ESH
               MCORE = (QH1D*iDT+QCLih+QCLsh)
               HHwcd = 2.*PI*NH1D*KAP*VENQH*(TK0C-TK1D)
               HHwrm = CPW*(TK0C-TK1D)*DMWDT
               MLWM  = MAX(0.,-1.*(HHwcd+HHwqv+HHwrm)/XXLF0)
               IF (MLWM.GT.DMWDT) THEN
                  QHwml = MAX(-1.*MCORE,DMWDT-MLWM)
                  AHwml = QHwml*4.*iRHOH/PI/MVDH
               ENDIF
               MLWC = 0.268+0.1389*MCORE*1.E3
               IF ((MLWM*1.E3).GT.MLWC) THEN
                  DH9   = MVDH*1.E2-9.E-1
                  MWT   = (QH1D*iDT+QCLih+QCLsh+DMWDT)*1.E3
                  ICOR1 = MWT/(1.+10.67*DH9-10.81*DH9**2.+10.26*DH9**3.)
                  QHwsh = MAX(0.,MLWM-MAX(0.,ICOR1*1.E-3))
                  NHwsh = QHwsh/(AMW*DSHED**BMW)
               ENDIF
            ELSE                                                        ! DRY-GROWTH MODE
               QHdrm = DMWDT
               AHdrm = ARMch+ARMrh
            ENDIF
         ENDIF
      ELSEIF (HWET_MODE.EQ.0) THEN
         QHdrm = QRMch+QRMrh
         AHdrm = ARMch+ARMrh
      ENDIF                                                             ! HAILSTONE WET GROWTH MODE OPTION
!----------------- CONSERVATION OF HYDROMETEORS ------------------------
      QCSOUR = QC1D+(QBKrc-QMLic)*DT
      QCSINK = (QHOci+QNMci+QNCci+QRMci+QRMcs+QRMcg+QRMch+QCLcr+QCNcr+ &
               QIMcsi+QIMcgi)*DT
      IF (QCSINK.GT.QCSOUR.AND.QCSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QCSOUR/(QCSINK+QSMALL))
         QHOci = QHOci*RATIO; QNMci = QNMci*RATIO; QNCci = QNCci*RATIO
         QRMci = QRMci*RATIO; QRMcs = QRMcs*RATIO; QRMcg = QRMcg*RATIO
         QRMch = QRMch*RATIO; QCLcr = QCLcr*RATIO; QCNcr = QCNcr*RATIO
         QIMcsi = QIMcsi*RATIO; QIMcgi = QIMcgi*RATIO
      ENDIF
      NCSOUR = NC1D+(NBKrc-NMLic)*DT
      NCSINK = (NHOci+NNMci+NNCci+NRMci+NRMcs+NRMcg+NRMch+NCLcr+NCLcc+ &
               NCNcr+NIMcsi+NIMcgi)*DT
      IF (NCSINK.GT.NCSOUR.AND.NCSOUR.GE.NSMALL) THEN
         RATIO = MIN(1.,NCSOUR/(NCSINK+NSMALL))
         NHOci = NHOci*RATIO; NNMci = NNMci*RATIO; NNCci = NNCci*RATIO
         NRMci = NRMci*RATIO; NRMcs = NRMcs*RATIO; NRMcg = NRMcg*RATIO
         NRMch = NRMch*RATIO; NCLcr = NCLcr*RATIO; NCLcc = NCLcc*RATIO
         NCNcr = NCNcr*RATIO; NIMcsi = NIMcsi*RATIO
         NIMcgi = NIMcgi*RATIO
      ENDIF
      QRSOUR = QR1D+(QCLcr+QCNcr-QMLir-QMLsr-QMLgr-QMLhr+QHwsh-QHwml)*DT
      QRSINK = (QHOrg+QNMrg+QRMri+QRMrs+QRMrg+QRMrh+QBKrc+QIMrsi+      &
               QIMrgi)*DT
      IF (QRSINK.GT.QRSOUR.AND.QRSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QRSOUR/(QRSINK+QSMALL))
         QHOrg = QHOrg*RATIO; QNMrg = QNMrg*RATIO; QRMri = QRMri*RATIO
         QRMrs = QRMrs*RATIO; QRMrg = QRMrg*RATIO; QRMrh = QRMrh*RATIO
         QBKrc = QBKrc*RATIO; QIMrsi = QIMrsi*RATIO
         QIMrgi = QIMrgi*RATIO
      ENDIF
      NRSOUR = NR1D+(NCNcr+NBKrr-NMLir-NMLsr-NMLgr-NMLhr+NCLcr+NHwsh)*DT
      NRSINK = (NHOrg+NNMrg+NRMri+NRMrs+NRMrg+NRMrh+NCLrr+NBKrc+NIMrsi+&
                NIMrgi)*DT
      IF (NRSINK.GT.NRSOUR.AND.NRSOUR.GE.NSMALL) THEN
         RATIO = MIN(1.,NRSOUR/(NRSINK+NSMALL))
         NHOrg = NHOrg*RATIO; NNMrg = NNMrg*RATIO; NRMri = NRMri*RATIO
         NRMrs = NRMrs*RATIO; NRMrg = NRMrg*RATIO; NRMrh = NRMrh*RATIO
         NCLrr = NCLrr*RATIO; NBKrc = NBKrc*RATIO; NIMrsi = NIMrsi*RATIO
         NIMrgi = NIMrgi*RATIO
      ENDIF
      QISOUR = QI1D+(QIMcsi+QIMcgi+QIMrsi+QIMrgi+QHOci+QNMci+QNCci+    &
                     QRMci)*DT
      QISINK = (-QMLir-QMLic+QCLir+QCLis+QCLig+QCLih+QINig+QCNis)*DT
      IF (QISINK.GT.QISOUR.AND.QISOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QISOUR/(QISINK+QSMALL))
         QMLir = QMLir*RATIO; QMLic = QMLic*RATIO; QCLir = QCLir*RATIO
         QCLis = QCLis*RATIO; QCLig = QCLig*RATIO; QCLih = QCLih*RATIO
         QINig = QINig*RATIO; QCNis = QCNis*RATIO
      ENDIF
      NISOUR = NI1D+(NIMcsi+NIMcgi+NIMrsi+NIMrgi+NHOci+NNMci+NNCci)*DT
      NISINK = (-NMLir-NMLic+NCLir+NCLis+NCLig+NCLih+NINig+NiCNis)*DT
      IF (NISINK.GT.NISOUR.AND.NISOUR.GE.NSMALL) THEN
         RATIO = MIN(1.,NISOUR/(NISINK+NSMALL))
         NCLis = NCLis*RATIO; NCLig = NCLig*RATIO; NCLih = NCLih*RATIO
         NINig = NINig*RATIO; NMLir = NMLir*RATIO; NMLic = NMLic*RATIO
         NCLir = NCLir*RATIO; NiCNis = NiCNis*RATIO
      ENDIF
      VISOUR = VI1D+(VIMcsi+VIMcgi+VIMrsi+VIMrgi+VHOci+VNMci+VNCci+    &
                     VRMci)*DT
      VISINK = (-VMLir-VMLic+VCLir+ViCLis+ViCLig+VCLih+ViINig+ViCNis)*DT
      IF (VISINK.GT.VISOUR.AND.VISOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,VISOUR/(VISINK+ISMALL))
         VMLir = VMLir*RATIO; VMLic = VMLic*RATIO; VCLir = VCLir*RATIO
         ViCLis = ViCLis*RATIO; ViCLig = ViCLig*RATIO
         VCLih = VCLih*RATIO; ViINig = ViINig*RATIO
         ViCNis = ViCNis*RATIO
      ENDIF
      FISOUR = FI1D+(FIMcsi+FIMcgi+FIMrsi+FIMrgi+FHOci+FNMci+FNCci+    &
                     FRMci)*DT
      FISINK = (-FMLir-FMLic+FCLir+FiCLis+FCLig+FCLih+FINig+FiCNis)*DT
      IF (FISINK.GT.FISOUR.AND.FISOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,FISOUR/(FISINK+ISMALL))
         FMLir = FMLir*RATIO; FMLic = FMLic*RATIO; FCLir = FCLir*RATIO
         FiCLis = FiCLis*RATIO; FCLig = FCLig*RATIO; FCLih = FCLih*RATIO
         FINig = FINig*RATIO; FiCNis = FiCNis*RATIO
      ENDIF
      AISOUR = AI1D+(AIMcsi+AIMcgi+AIMrsi+AIMrgi+AHOci+ANMci+ANCci+    &
                     ARMci)*DT
      AISINK = (ACLir-AMLir-AMLic+AiCLis+AiCLig+AiCLih+AiINig+AiCNis)*DT
      IF (AISINK.GT.AISOUR.AND.AISOUR.GE.ASMALL) THEN
         RATIO = MIN(1.,AISOUR/(AISINK+ASMALL))
         AMLir = AMLir*RATIO; AMLic = AMLic*RATIO
         ACLir = ACLir*RATIO; AiCLis = AiCLis*RATIO
         AiCLig = AiCLig*RATIO; AiCLih = AiCLih*RATIO
         AiINig = AiINig*RATIO; AiCNis = AiCNis*RATIO
      ENDIF
      IISOUR = I3M1D+(IIMcsi+IIMcgi+IIMrsi+IIMrgi+IHOci+INMci+INCci+   &
                     IRMci)*DT
      IISINK = (-IMLir-IMLic+ICLir+ICLis+ICLig+ICLih+IINig+ICNis)*DT
      IF (IISINK.GT.IISOUR.AND.IISOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,IISOUR/(IISINK+ISMALL))
         IMLir = IMLir*RATIO; IMLic = IMLic*RATIO; ICLir = ICLir*RATIO
         ICLis = ICLis*RATIO; ICLig = ICLig*RATIO; ICLih = ICLih*RATIO
         IINig = IINig*RATIO; ICNis = ICNis*RATIO
      ENDIF
      QSSOUR = QS1D+(QCNis+QCLis+QRMcs)*DT
      QSSINK = (-QMLsr+QCLsr+QCLsg+QCLsh+QINsg)*DT
      IF (QSSINK.GT.QSSOUR.AND.QSSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QSSOUR/(QSSINK+QSMALL))
         QMLsr = QMLsr*RATIO; QCLsg = QCLsg*RATIO; QCLsh = QCLsh*RATIO
         QINsg = QINsg*RATIO; QCLsr = QCLsr*RATIO
      ENDIF
      NSSOUR = NS1D+(NsCNis)*DT
      NSSINK = (-NMLsr+NCLsr+NCLsg+NCLsh+NCLss+NINsg)*DT
      IF (NSSINK.GT.NSSOUR.AND.NSSOUR.GE.NSMALL) THEN
         RATIO = MIN(1.,NSSOUR/(NSSINK+NSMALL))
         NMLsr = NMLsr*RATIO; NCLsg = NCLsg*RATIO; NCLsh = NCLsh*RATIO
         NCLss = NCLss*RATIO; NINsg = NINsg*RATIO; NCLsr = NCLsr*RATIO
      ENDIF
      VSSOUR = VS1D+(VsCNis+VsCLis+VRMcs+VCLss)*DT
      VSSINK = (-VMLsr+VCLsr+VsCLsg+VCLsh+VsINsg)*DT
      IF (VSSINK.GT.VSSOUR.AND.VSSOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,VSSOUR/(VSSINK+ISMALL))
         VMLsr = VMLsr*RATIO; VsCLsg = VsCLsg*RATIO; VCLsh = VCLsh*RATIO
         VsINsg = VsINsg*RATIO; VCLsr = VCLsr*RATIO
      ENDIF
      FSSOUR = FS1D+(FsCNis+FsCLis+FRMcs+FCLss)*DT
      FSSINK = (-FMLsr+FCLsr+FCLsg+FCLsh+FINsg)*DT
      IF (FSSINK.GT.FSSOUR.AND.FSSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,FSSOUR/(FSSINK+QSMALL))
         FMLsr = FMLsr*RATIO; FCLsg = FCLsg*RATIO; FCLsh = FCLsh*RATIO
         FINsg = FINsg*RATIO; FCLsr = FCLsr*RATIO
      ENDIF
      ASSOUR = AS1D+(AsCNis+AsCLis+ARMcs)*DT
      ASSINK = (-AMLsr+ACLsr+AsCLsg+AsCLsh+AsINsg+ACLss)*DT
      IF (ASSINK.GT.ASSOUR.AND.ASSOUR.GE.ASMALL) THEN
         RATIO = MIN(1.,ASSOUR/(ASSINK+ASMALL))
         AMLsr = AMLsr*RATIO; AsINsg = AsINsg*RATIO; ACLsr = ACLsr*RATIO
         AsCLsg = AsCLsg*RATIO; AsCLsh = AsCLsh*RATIO
         ACLss = ACLss*RATIO
      ENDIF
      QGSOUR = QG1D+(QHOrg+QNMrg+QCLig+QCLsg+QINig+QINsg+QRMcg+QCLirg+ &
               QCLsrg+QCLgrg)*DT
      QGSINK = (-QMLgr+QCNgh+QCLgr)*DT
      IF (QGSINK.GT.QGSOUR.AND.QGSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QGSOUR/(QGSINK+QSMALL))
         QMLgr = QMLgr*RATIO; QCNgh = QCNgh*RATIO; QCLgr = QCLgr*RATIO
      ENDIF
      NGSOUR = NG1D+(NHOrg+NNMrg+NINig+NINsg+NCLirg+NCLsrg+NCLgrg)*DT
      NGSINK = (-NMLgr+NgCNgh+NCLgr)*DT
      IF (NGSINK.GT.NGSOUR.AND.NGSOUR.GE.NSMALL) THEN
         RATIO = MIN(1.,NGSOUR/(NGSINK+NSMALL))
         NMLgr = NMLgr*RATIO; NgCNgh = NgCNgh*RATIO; NCLgr = NCLgr*RATIO
      ENDIF
      VGSOUR = VG1D+(VHOrg+VNMrg+VgCLig+VgCLsg+VgINig+VgINsg+VRMcg+    &
               VCLirg+VCLsrg+VCLgrg)*DT
      VGSINK = (-VMLgr+VCNgh+VCLgr)*DT
      IF (VGSINK.GT.VGSOUR.AND.VGSOUR.GE.ISMALL) THEN
         RATIO = MIN(1.,VGSOUR/(VGSINK+ISMALL))
         VMLgr = VMLgr*RATIO; VCNgh = VCNgh*RATIO; VCLgr = VCLgr*RATIO
      ENDIF
      AGSOUR = AG1D+(AHOrg+ANMrg+AgCLig+AgCLsg+AgINig+AgINsg+ARMcg+    &
               ACLirg+ACLsrg+ACLgrg)*DT
      AGSINK = (-AMLgr+AgCNgh+ACLgr)*DT
      IF (AGSINK.GT.AGSOUR.AND.AGSOUR.GE.ASMALL) THEN
         RATIO = MIN(1.,AGSOUR/(AGSINK+ASMALL))
         AMLgr = AMLgr*RATIO; AgCNgh = AgCNgh*RATIO; ACLgr = ACLgr*RATIO
      ENDIF
      QHSOUR = QH1D+(QCNgh+QCLih+QCLsh+QHdrm)*DT
      QHSINK = (-QMLhr-QHwml)*DT
      IF (QHSINK.GT.QHSOUR.AND.QHSOUR.GE.QSMALL) THEN
         RATIO = MIN(1.,QHSOUR/(QHSINK+QSMALL))
         QMLhr = QMLhr*RATIO; QHwml = QHwml*RATIO
      ENDIF
      NHSOUR = NH1D+(NhCNgh)*DT
      NHSINK = (-NMLhr)*DT
      IF (NHSINK.GT.NHSOUR.AND.NHSOUR.GE.NSMALL) THEN
         RATIO = MIN(1.,NHSOUR/(NHSINK+NSMALL))
         NMLhr = NMLhr*RATIO
      ENDIF
      AHSOUR = AH1D+(AhCNgh+AhCLih+AhCLsh+AHdrm)*DT
      AHSINK = (-AMLhr-AHwml)*DT
      IF (AHSINK.GT.AHSOUR.AND.AHSOUR.GE.ASMALL) THEN
         RATIO = MIN(1.,AHSOUR/(AHSINK+ASMALL))
         AMLhr = AMLhr*RATIO; AHwml = AHwml*RATIO
      ENDIF

      QFZci = QHOci+QNMci+QNCci
      NFZci = NHOci+NNMci+NNCci
      VFZci = VHOci+VNMci+VNCci
      FFZci = FHOci+FNMci+FNCci
      AFZci = AHOci+ANMci+ANCci
      IFZci = IHOci+INMci+INCci
      QFZrg = QHOrg+QNMrg
      NFZrg = NHOrg+NNMrg
      VFZrg = VHOrg+VNMrg
      AFZrg = AHOrg+ANMrg
      QIMii = QIMcsi+QIMcgi+QIMrsi+QIMrgi
      NIMii = NIMcsi+NIMcgi+NIMrsi+NIMrgi
      VIMii = VIMcsi+VIMcgi+VIMrsi+VIMrgi
      AIMii = AIMcsi+AIMcgi+AIMrsi+AIMrgi
      FIMii = FIMcsi+FIMcgi+FIMrsi+FIMrgi
      IIMii = IIMcsi+IIMcgi+IIMrsi+IIMrgi
      GQCTR = (QRMcs+QRMcg+QRMch+QCLir+QCLcr+QCLis+QCLig+QCLih+QCNcr+  &
              QCNis+QINig-QBKrc-QIMcsi-QIMcgi-QMLir)*DT
!----------------- UPDATE VARIABLES ------------------------------------
      QC1D = MAX(0.,QC1D+(QBKrc-QCLcr-QCNcr-QFZci-QMLic-QRMci-QRMcs-   &
             QRMcg-QRMch-QIMcsi-QIMcgi)*DT)
      QR1D = MAX(0.,QR1D+(QCLcr+QCNcr-QBKrc-QFZrg-QIMrsi-QIMrgi-QMLir- &
             QMLsr-QMLgr-QMLhr-QRMri-QRMrs-QRMrg-QRMrh+QHwsh-QHwml)*DT)
      QI1D = MAX(0.,QI1D+(QIMii-QCNis-QCLis-QCLig-QCLih-QINig+QFZci+   &
             QMLir+QMLic+QRMci-QCLir)*DT)
      QS1D = MAX(0.,QS1D+(QCNis+QCLis-QCLsg-QCLsh-QINsg+QMLsr+QRMcs-   &
             QCLsr)*DT)
      QG1D = MAX(0.,QG1D+(QFZrg+QCLig+QCLsg+QINig+QINsg+QMLgr+QRMcg-   &
             QCNgh-QCLgr+QCLirg+QCLsrg+QCLgrg)*DT)
      QH1D = MAX(0.,QH1D+(QMLhr+QCNgh+QHdrm+QCLih+QCLsh+QHwml)*DT)
      NC1D = MAX(0.,NC1D+(NBKrc-NCLcr-NCLcc-NCNcr-NFZci-NRMci-NRMcs-   &
             NRMcg-NRMch-NMLic-NIMcsi-NIMcgi)*DT)
      NR1D = MAX(0.,NR1D+(NCNcr+NBKrr-NBKrc-NCLrr-NIMrsi-NIMrgi-NFZrg- &
             NMLir-NMLsr-NMLgr-NMLhr-NRMri-NRMrs-NRMrg-NRMrh+NHwsh)*DT)
      NI1D = MAX(0.,NI1D+(NIMii-NiCNis-NCLis-NCLig-NCLih-NINig+NFZci+  &
             NMLir+NMLic-NCLir)*DT)
      NS1D = MAX(0.,NS1D+(NsCNis-NCLsg-NCLsh-NCLss-NINsg+NMLsr-NCLsr)* &
             DT)
      NG1D = MAX(0.,NG1D+(NFZrg+NINig+NINsg+NMLgr-NgCNgh-NCLgr+NCLirg+ &
             NCLsrg+NCLgrg)*DT)
      NH1D = MAX(0.,NH1D+(NMLhr+NhCNgh)*DT)
      IF (ICE_RHOI.EQ.0.OR.ICE_RHOI.EQ.2) THEN
         VI1D = 0.
      ELSEIF (ICE_RHOI.EQ.1) THEN
         VI1D = MAX(0.,VI1D+(VIMii-ViCNis-ViCLis-ViCLig-VCLih-ViINig+  &
                VFZci+VMLir+VMLic+VRMci-VCLir)*DT)
      ENDIF
      IF (ICE_RHOS.EQ.0.OR.ICE_RHOS.EQ.2) THEN
         VS1D = 0.
      ELSEIF (ICE_RHOS.EQ.1) THEN
         VS1D = MAX(0.,VS1D+(VsCNis+VsCLis-VsCLsg-VCLsh-VsINsg+VMLsr+  &
                VRMcs-VCLsr+VCLss)*DT)
      ENDIF
      IF (ICE_RHOG.EQ.0) THEN
         VG1D = 0.
      ELSEIF (ICE_RHOG.EQ.1.OR.ICE_RHOG.EQ.2) THEN
         VG1D = MAX(0.,VG1D+(VFZrg+VgCLig+VgCLsg+VgINig+VgINsg+VMLgr+  &
                VRMcg-VCNgh-VCLgr+VCLirg+VCLsrg+VCLgrg)*DT)
      ENDIF
      FI1D = MAX(0.,FI1D+(FIMii-FiCNis-FiCLis-FCLig-FCLih-FINig+FFZci+ &
             FMLir+FMLic+FRMci-FCLir)*DT)
      IF (AGG_SHAPE.EQ.0) THEN
         FS1D = 0.
      ELSEIF (AGG_SHAPE.EQ.1) THEN
         FS1D = MAX(0.,FS1D+(FsCNis+FsCLis-FCLsg-FCLsh-FINsg+FMLsr+    &
                FRMcs-FCLsr+FCLss)*DT)
      ENDIF
      I3M1D = MAX(0.,I3M1D+(IIMii-ICNis-ICLis-ICLig-ICLih-IINig+IFZci+ &
              IMLir+IMLic+IRMci-ICLir)*DT)
      IF (AFAI_3M.EQ.0.OR.AFAI_3M.EQ.2) THEN
         AI1D = 0.
      ELSEIF (AFAI_3M.EQ.1) THEN
         AI1D = MAX(0.,AI1D+(AIMii-AiCNis-AiCLis-AiCLig-AiCLih-AiINig+ &
                AFZci+AMLir+AMLic+ARMci-ACLir)*DT)
      ENDIF
      IF (AFAS_3M.EQ.0.OR.AFAS_3M.EQ.2) THEN
         AS1D = 0.
      ELSEIF (AFAS_3M.EQ.1) THEN
         AS1D = MAX(0.,AS1D+(AsCNis+AsCLis-AsCLsg-AsCLsh-AsINsg+AMLsr+ &
                ARMcs-ACLsr-ACLss)*DT)
      ENDIF
      IF (AFAG_3M.EQ.0.OR.AFAG_3M.EQ.2) THEN
         AG1D = 0.
      ELSEIF (AFAG_3M.EQ.1) THEN
         AG1D = MAX(0.,AG1D+(AFZrg+AgCLig+AgCLsg+AgINig+AgINsg+AMLgr+  &
                ARMcg-AgCNgh-ACLgr+ACLirg+ACLsrg+ACLgrg)*DT)
      ENDIF
      IF (AFAH_3M.EQ.0.OR.AFAH_3M.EQ.2) THEN
         AH1D = 0.
      ELSEIF (AFAH_3M.EQ.1) THEN
         AH1D = MAX(0.,AH1D+(AMLhr+AhCNgh+AHdrm+AhCLih+AhCLsh+AHwml)*DT)
      ENDIF
!      QNIN  = MAX(0.,QNIN-(NNCci+NNMci)*DT)
      TK1D = TK1D+((QFZci+QFZrg+QMLir+QMLic+QMLsr+QMLgr+QMLhr+QRMci+   &
             QRMcs+QRMcg+QRMri+QRMrs+QRMrg+QIMii+QHdrm+QHwml)*XXLF)/   &
             CPM*DT

      IF (QC1D.GE.QSMALL.AND.NC1D.GE.NSMALL) THEN
         MVDC = (QC1D*iAMW/NC1D)**THRD
         IF (MVDC.GT.DCR) THEN
            QR1D = QR1D+QC1D
            NR1D = NR1D+NC1D
            QC1D = 0.; NC1D = 0.
         ENDIF
      ENDIF
      IF (QR1D.GE.QSMALL.AND.NR1D.GE.NSMALL) THEN
         MVDR = (QR1D*iAMW/NR1D)**THRD
         IF (MVDR.LT.DCR) THEN
            QC1D = QC1D+QR1D
            NC1D = NC1D+NR1D
            QR1D = 0.; NR1D = 0.
         ENDIF
      ENDIF
      IF (QH1D.GE.QSMALL.AND.NH1D.GE.NSMALL) THEN
         MVDH = (QH1D*iAMH/NH1D)**THRD
         IF (MVDH.LT.DHMIN) THEN
            QG1D = QG1D+QH1D
            NG1D = NG1D+NH1D
            VG1D = VG1D+QH1D/RHOG0
            QH1D = 0.; NH1D = 0.
            IF (AH1D.GE.ASMALL.AND.AFAG_3M.EQ.1) THEN
               AG1D = AG1D+AH1D
               AH1D = 0.
            ENDIF
         ENDIF
      ENDIF
      IF (GQCTR.GT.0.) THEN
         TQCI = QC1D+QI1D 
         IF (TQCI.GT.RLIMIT*2.) THEN
            GQCTR = MIN(1.,GQCTR/TQCI)
         ELSE
            GQCTR = 0.
         ENDIF
      ELSE
         TQRSG = QR1D+QS1D+QG1D+QH1D
         IF (TQRSG.GT.RLIMIT*3.) THEN 
            GQCTR = MAX(-1.,GQCTR/TQRSG) 
         ELSE 
            GQCTR = 0.
         ENDIF 
      ENDIF

      END SUBROUTINE LARGE_DT
!======================================================================
      END MODULE MODULE_MP_NTU
!======================================================================
