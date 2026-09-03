!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 
!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 

 module module_nssl_2mmt
 use module_const
 implicit none


 private  !everything is unless declare it as public. 

!
! public subroutines
!
 public  :: re_nssl_volume


  integer, private :: eqtset = 1 ! Flag for use with cm1 to use alternate equation set (changes latent heating rates)
   double precision, parameter, public :: zscale = 1.0d0 ! 1.000e-10
   double precision, parameter, public :: zscaleinv = 1.0d0/zscale ! 1.000e-10


  real, parameter :: warmonly = 0.0 ! testing parameter, set to 1.0 to reduce to warm-rain physics (ice variables stay zero)

  logical, parameter :: lwsm6 = .false. ! act like wsm6 for some single moment interactions

! some constants from WSM6
  real, parameter  :: dimax = 500.e-6    ! limited maximum value for the cloud-ice diamter
  real, parameter  :: roqimax = 2.08e22*dimax**8

! Params for dbz:
  integer  :: iuseferrier = 1  ! =1: use dry graupel only from Ferrier 1994; = 0: Use Smith (wet graupel)
  integer  :: idbzci      = 0
  integer  :: iusewetgraupel = 1 ! =1 to turn on use of QHW for graupel reflectivity (only for ZVDM -- mixedphase)
                                 ! =2 turn on for graupel density less than 300. only 
  integer  :: iusewethail = 0 ! =1 to turn on use of QHW for graupel reflectivity (only for ZVDM -- mixedphase)
  integer  :: iusewetsnow = 1 ! =1 to turn on diagnosed bright band

! microphysics

  real, private :: rho_qr = 1000., cnor = 8.0e5  ! cnor is set in namelist!!  rain params
  real, private :: rho_qs =  100., cnos = 3.0e6  ! set in namelist!!  snow params
  real, private :: rho_qh =  500., cnoh = 4.0e5  ! set in namelist!!  graupel params
  real, private :: rho_qhl=  900., cnohl = 4.0e4 ! set in namelist!!  hail params

  real, private :: hdnmn  = 170.0  ! minimum graupel density (for variable density graupel)
  real, private :: hldnmn = 500.0  ! minimum hail density (for variable density hail)

  real :: cnohmn  = 1.e-2 ! minimum intercept for 2-moment graupel (alphah < 0.5)
  real :: cnohlmn = 1.e-2 ! minimum intercept for 2-moment hail (alphahl < 0.5)

! Autoconversion parameters

  real   , private :: qcmincwrn      = 2.0e-3    ! qc threshold for autonconversion (LFO; for 10ICE use qminrncw for ircnw != 5)
  real   , private :: cwdiap         = 20.0e-6   ! threshold diameter of cloud drops (Ferrier 1994 autoconversion)
  real   , private :: cwdisp         = 0.15      ! assume droplet dispersion parameter (can be 0.3 for maritime)
  real   , private  :: ccn            = 0.6e+09   ! set in namelist!! Central plains CCN value
  real   , private :: qccn             ! ccn "mixing ratio"
  integer, private :: iauttim        = 1         ! 10-ice rain delay flag
  real   , private :: auttim         = 300.      ! 10-ice rain delay time
  real   , private :: qcwmntim       = 1.0e-5    ! 10-ice rain delay min qc for time accrual

  logical, parameter :: invertccn = .false. ! =true for base state of ccn=0, =false for ccn initialized in the base state
  logical :: restoreccn = .false. ! whether or not to nudge CCN back to base state (qccn) (only applies if CCNA is NOT predicted)
  real    :: ccntimeconst = 600.  ! time constant for CCN restore (either for CCNA or when restoreccn = true)


! sedimentation flags
! itfall -> 0 = 1st order fallout (other options removed)
! iscfall, infall -> fallout options for charge and number concentration, respectively
!                    1 = mass-weighted fall speed; 2 = number-weighted fallspeed.
  integer, private :: itfall = 0
  integer, private :: iscfall = 1
  integer, private :: irfall = -1
  logical, private :: do_accurate_sedimentation = .false. ! if true, recalculate fall speeds on sub time steps; (more expensive)
                                                         ! if false, reuse fall speeds on multiple steps (can have a noticeable speedup)
                                                         ! Mainly is an issue for small dz near the surface. 
  integer, private :: infall = 4   ! 0 -> uses number-wgt for N; NO correction applied (results in excessive size sorting)
                          ! 1 -> uses mass-weighted fallspeed for N ALWAYS
                          ! 2 -> uses number-wgt for N and mass-weighted correction for N (Method II in Mansell, 2010 JAS)
                          ! 3 -> uses number-wgt for N and Z-weighted correction for N (Method I in Mansell, 2010 JAS)
                          ! 4 -> Hybrid of 2 and 3: Uses minimum N from each method (z-wgt and m-wgt corrections) (Method I+II in Mansell, 2010 JAS)
                          ! 5 -> uses number-wgt for N and uses average of N-wgt and q-wgt instead of Max.
  real    :: rainfallfac = 1.0 ! factor to adjust rain fall speed (single moment only)
  integer, private :: icdx = 3 ! (graupel) 0=Ferrier; 1=leave drag coef. cd fixed; 2=vary by density, 4=set by user with cdxmin,cdxmax,etc.
  integer, private :: icdxhl = 3 ! (hail) 0=Ferrier; 1=leave drag coef. cd fixed; 2=vary by density, 4=set by user with cdxmin,cdxmax,etc.
  real   , private :: cdhmin = 0.45, cdhmax = 0.8        ! defaults for graupel (icdx=4)
  real   , private :: cdhdnmin = 500., cdhdnmax = 800.0  ! defaults for graupel (icdx=4)
  real   , private :: cdhlmin = 0.45, cdhlmax = 0.6      ! defaults for hail (icdx=4)
  real   , private :: cdhldnmin = 500., cdhldnmax = 800.0  ! defaults for hail (icdx=4)
  real   , private :: vtmaxsed = 70. ! Limit on fall speed (m/s, all moments) for sedimentation calculations. Not applied to fall speeds for microphysical rates

  integer :: rssflg = 1   ! Rain size-sorting allowed (1, default), or disallowed (0).  If 0, sets N and Z-weighted fall speeds to q-weighted value
  integer :: sssflg = 1   ! As above but for snow
  integer :: hssflg = 1   ! As above but for graupel
  integer :: hlssflg = 1  ! As above but for hail

! input flags

  integer, private :: ndebug = -1, ncdebug = 0
  integer, private :: ipconc = 5
  integer, private :: ichaff = 0
  integer, private :: ilimit = 0

  real, private :: constccw = -1.

  real, private :: cimn = 1.0e3, cimx = 1.0e6


  real   , private :: ifrzg = 1.0 ! fraction of frozen drops going to graupel. 1=freeze all rain to graupel, 0=freeze all to hail
  real   , private :: ifrzs = 1.0 ! fraction of small frozen drops going to snow. 1=freeze rain to snow, 0=freeze to cloud ice
  real   , private :: ffrzs = 0.0 ! fraction of other initiated cloud ice going to snow. 1=freeze rain to snow, 0=freeze to cloud ice
  integer, private :: irwfrz = 1 ! compute total rain that can freeze (checks heat budget)
  integer, private :: irimtim = 0 ! future use
!  integer, private :: infdo = 1   ! 1 = calculate number-weighted fall speeds

  integer, private :: irimdenopt = 1 ! = 1 for default Maklin; = 2 for experimental Cober and List (1993)
  real   , private :: rimc1 = 300.0, rimc2 = 0.44  ! rime density coeff. and power (Default Heymsfield and Pflaum, 1985)
  real   , private :: rimc3 = 170.0                ! minimum rime density
  real    :: rimc4 = 900.0                ! maximum rime density
  real   , private :: rimtim = 120.0               ! cut-off rime time (10ICE)
  real   , private :: eqtot = 1.0e-9               ! threshold for mass budget reporting

  integer, private :: ireadmic = 0

  integer, private :: iccwflg = 1     ! sets max size of first droplets in parcel to 4 micron radius (in two-moment liquid)
                             ! (first nucleation is done with a KW sat. adj. step)
  integer, private :: issfilt = 0     ! flag to turn on filtering of supersaturation field
  integer, private :: irenuc = 2      ! =1 to always allow renucleation of droplets within the cloud
                                      ! =2 renucleation following Twomey/Cohard&Pinty
                                      ! =7 New renucleation that requires prediction of the number of activated nuclei
                             ! i.e., not only at cloud base
  integer, private :: irenuc3d = 0      ! =1 to include horizontal gradient in renucleation of droplets within the cloud
  real    :: renucfrac = 0.0 ! = 0 : cnuc = cwccn
                             ! = 1 : cnuc = actual available CCN
                             ! otherwise cnuc = cwccn*(1. - renufrac) + ccnc(1:ngscnt)*renucfrac
  real   , private :: cck = 0.6       ! exponent in Twomey expression
  real   , private :: ciintmx = 1.0e6 ! limit on ice concentration from primary nucleation

  real   , private :: cwccn ! , cwmasn,cwmasx
  real   , private :: ccwmx

  integer, private :: idocw = 1, idorw = 1, idoci = 1, idoir = 1, idoip = 1, idosw = 1
  integer, private :: idogl = 1, idogm = 1, idogh = 1, idofw = 1, idohw = 1, idohl = 1
!  integer, private :: ido(3:14) = / 12*1 /


! 0,2, 5.00e-10, 1, 0, 0, 0      : itype1,itype2,cimas0,icfn,ihrn,ibfc,iacr
  integer, private :: itype1 = 0, itype2 = 2  ! controls Hallett-Mossop process
  integer, private :: icenucopt = 1       ! =1 Meyers/Ferrier primary ice nucleation; =2 Thompson/Cooper, =3 Phillips (Meyers/Demott)
  integer, private :: icfn = 2                ! contact freezing: 0 = off; 1 = hack (ok for single moment); 2 = full Cotton/Meyers version
  integer, private :: ihrn = 0            ! Hobbs-Rangno ice multiplication (Ferrier, 1994; use in 10-ice only)
  integer, private :: ibfc = 1            ! Flag to use Bigg freezing on droplets (0 = off (uses alternate freezing), 1 = on)
  integer, private :: iremoveqwfrz = 1    ! Whether to remove (=1) or not (=0) the newly-frozen cloud droplets (ibfc=1) from the CWC used for charge separation
  integer, private :: iacr = 2            ! Flag for drop contact freezing with crytals
                                 ! (0=off; 1=drops > 500micron diameter; 2 = > 300micron)
  integer, private :: ibfr = 2            ! Flag for Bigg freezing conversion of freezing drops to graupel
                                 ! (1=min graupel size is vr1mm; 2=use min size of dfrz, 5= as for 2 and apply dbz conservation)
  integer, private :: ibiggopt = 1        ! 1 = old Bigg; 2 = experimental Bigg (only for imurain = 1, however)
  integer :: ibiggsmallrain = 0  ! 1 = When rain is too small, freeze none to graupel and send all to snow (experimental)
  integer, private :: iacrsize = 5        ! assumed min size of drops freezing by capture
                                 !  1: > 500 micron diam
                                 !  2: > 300 micron
                                 !  3: > 40 micron
                                 !  4: all sizes
                                 !  5: > 150 micron (only for imurain = 1)
  real   , private :: cimas0 = 6.62e-11   ! default mass of Hallett-Mossop crystals
                                 ! 6.62e-11kg results in half the diam. (60 microns) of old default value of 5.0e-10
  real   , private :: cimas1 = 6.88e-13   ! default mass of new ice crystals
  real   , private :: splintermass = 6.88e-13
  real   , private :: cfnfac = 0.1        ! Hack factor that goes with icfn=1
  integer, private :: iscni = 4           ! default option for ice crystal aggregation/conversion to snow
  real   , private :: fscni = 1.0         ! factor for calculating cscni
  logical, private :: imeyers5 = .false.  ! .false.=off, true=on for Meyers ice nucleation for temp > -5 C
  real   , private :: dmincw = 15.0e-6    ! minimum droplet diameter for collection for iehw=3
  integer, private :: iehw = 1            ! 0 -> ehw=ehw0; 1 -> old ehw; 2 -> test ehw with Mason table data
  integer, private :: iehlw = 1           ! 0 -> ehlw=ehlw0; 1 -> old ehlw; 2 -> test ehlw with Mason table data
                                 ! For ehw/ehlw = 1, ehw0/ehlw0 act as maximum limit on collection efficiency (defaults are 1.0)
  integer, private :: ierw = 1            ! for single-moment rain (LFO/Z)
  integer, private :: iehr0c = 0          ! 0 -> no collection for T > 0C;  1 -> turn on collection/shedding for T > 0C
  integer, private :: iehlr0c = 0         ! 0 -> no collection for T > 0C;  1 -> turn on collection/shedding for T > 0C
  real   , private :: ehw0 = 0.5          ! constant or max assumed graupel-droplet collection efficiency
  real   , private :: erw0 = 1.0          ! constant assumed rain-droplet collection efficiency
  real   , private :: ehlw0 = 0.75        ! constant or max assumed hail-droplet collection efficiency
  real    :: ehr0 = 1.0          ! constant or max assumed graupel-rain collection efficiency
  real    :: ehlr0 = 1.0         ! constant or max assumed hail-rain collection efficiency
  real   , private :: exwmindiam = 0.0    ! minimum diameter of droplets for riming. If set > 0, will exclude that fraction of mass/number from accretion (idea from Furtado and Field 2017 JAS but also Fierro and Mansell 2017)


  real   , private :: esilfo0 = 1.0       ! factor for LFO collection efficiency of snow for cloud ice.
  real   , private :: ehslfo0 = 1.0       ! factor for LFO collection efficiency of hail/graupel for snow.

  integer, private :: ircnw    = 5        ! single-moment warm-rain autoconversion option.  5= Ferrier 1994.
  real   , private :: qminrncw = 2.0e-3   ! qc threshold for rain autoconversion (NA for ircnw=5)

  integer, private :: iqcinit = 2         ! For ZVDxx schemes, flag to choose which way to initialize droplets
                                 ! 1 = Soong-Ogura adjustment
                                 ! 2 = Saturation adjustment to value of ssmxinit
                                 ! 3 = KW adjustment

  real   , private :: ssmxinit = 0.4      ! saturation percentage to adjust down to for initial cloud
                                 ! formation (ZVDxx scheme only)

  real   , private :: ewfac = 1.0         ! hack factor applied to graupel and hail collection eff. for droplets
  real   , private :: eii0 = 0.1 ,eii1 = 0.1  ! graupel-crystal coll. eff. parameters: eii0*exp(eii1*min(temcg(mgs),0.0))
                                     ! set eii1 = 0 to get a constant value of eii0
  real   , private :: eii0hl = 0.2 ,eii1hl = 0.0  ! hail-crystal coll. eff. parameters: eii0hl*exp(eii1hl*min(temcg(mgs),0.0))
                                     ! set eii1hl = 0 to get a constant value of eii0hl
  real   , private :: eri0 = 0.1   ! rain efficiency to collect ice crystals
  real   , private :: eri_cimin = 10.e-6      ! minimum ice crystal diameter for collection by rain
  real   , private :: esi0 = 0.1              ! linear factor in snow-ice collection efficiency
  real   , private :: ehs0 = 0.1, ehs1 = 0.1  ! graupel-snow coll. eff. parameters: ehs0*exp(ehs1*min(temcg(mgs),0.0))
                                     ! set ehs1 = 0 to get a constant value of ehs0
  real   , private :: ess0 = 1.0, ess1 = 0.05 ! snow aggregation coefficients: ess0*exp(ess1*min(temcg(mgs),0.0))
                                     ! set ess1 = 0 to get a constant value of ess0
  real   , private :: esstem1 = -25.  ! lower temperature where snow aggregation turns on
  real   , private :: esstem2 = -20.  ! higher temperature for linear ramp of ess from zero at esstem1 to formula value at esstem2
  real   , private :: ehsfrac = 1.0           ! multiplier for graupel collection efficiency in wet growth
  real   , private :: ehimin = 0.0 ! Minimum collection efficiency (graupel - ice crystal)
  real   , private :: ehimax = 1.0 ! Maximum collection efficiency (graupel - ice crystal)
  real   , private :: ehsmax = 0.5 ! Maximum collection efficiency (graupel - snow)
  real   , private :: ecollmx = 0.5 ! Maximum collision efficiency for graup/hail with ice; used only for charging rates
  integer, private :: iglcnvi = 1  ! flag for riming conversion from cloud ice to rimed ice/graupel
  integer, private :: iglcnvs = 2  ! flag for conversion from snow to rimed ice/graupel

  real   , private :: rz          ! reflectivity conservation factor for graupel/rain
                         ! now calculated in icezvd_dr.F from alphah and rnu
                         ! currently only used for graupel melting to rain
  real   , private :: rzhl        ! reflectivity conservation factor for hail/rain
                         ! now calculated in icezvd_dr.F from alphahl and rnu

  real   , private :: rzs     ! reflectivity conservation factor for snow(imusnow=3) with rain (imurain=1)

  real   , private :: alphahacx = 0.0 ! assumed minimum shape parameter for zhacw and zhacr

  real   , private :: fconv = 1.0  ! factor to boost max graupel depletion by riming conversions in 10ICE

  real   , private :: rg0 = 400.0  ! reference graupel density for graupel fall speed

  integer, private :: rcond = 2    ! (Z only) rcond = 2 includes rain condensation in loop with droplet condensation
                                   ! 0 = no condensation on rain; 1 = bulk condensation on rain
  integer, parameter, private :: icond = 1    ! (Z only) icond = 1 calculates ice deposition (crystals and snow) BEFORE droplet condensation
                          ! icond = 2 does not work (intended to calc. dep in loop with droplet cond.)

  real   , private :: dfrz = 0.15e-3 ! 0.25e-3  ! minimum diameter of frozen drops from Bigg freezing (used for vfrz) for iacr > 1
                            ! and for ciacrf for iacr=4
  real   , private :: dmlt = 3.0e-3  ! maximum diameter for rain melting from graupel and hail
  real   , private :: dshd = 1.0e-3  ! nominal diameter for drops shed from graupel/hail

  integer, private :: ihmlt = 2      ! 1=old melting with vmlt; 2=new melting using mean volume diam of graupel/hail
  integer, private :: imltshddmr = 2 ! 0 (default)=mean diameter of drops produced during melting+shedding as before (using mean diameter of graupel/hail
                            ! and max mean diameter of rain)
                            ! 1=new method where mean diameter of rain during melting is adjusted linearly downward 
                            ! toward 3 mm for large (> sheddiam) graupel and hail, to take into account shedding of 
                            ! smaller drops.  sheddiam0 controls the size of graupel/hail above which the assumed 
                            ! mean diameter of rain is set to 3 mm
                            ! Only valid for ihmlt = 2 for ZVD(H) but also applies to ZVD(H)M
                            ! 2 = method that sets the resulting rain size ( vshdgs ) according to the mass-weighted diameter of the ice


  integer, private :: nsplinter = 0  ! number of ice splinters per freezing drop, if negative, then per resulting graupel particle
  real,    private :: lawson_splinter_fac = 2.5e-11  ! constant in Lawson et al. (2015, JAS) for ice particle production from freezing drops
  integer, private :: isnwfrac = 0   ! 0= no snow fragmentation; 1 = turn on snow fragmentation (Schuur, 2000)

!  integer, private :: denscale = 1  ! 1=scale num. conc. and charge by air density for advection, 0=turn off for comparison

  real, private  :: qhdpvdn = -1.
  real, private  :: qhacidn = -1.

  logical, private :: mixedphase = .false.   ! .false.=off, true=on to include mixed phase graupel
  integer, private :: imixedphase = 0
  logical, private :: qsdenmod = .false.     ! true = modify snow density by linear interpolation of snow and rain density
  logical, private :: qhdenmod = .false.     ! true = modify graupel density by linear interpolation of graupel and rain density
  logical, private :: qsvtmod = .false.      ! true = modify snow fall speed by linear interpolation of snow and rain vt
  real   , private :: sheddiam   = 8.0e-03  ! minimum diameter of graupel before shedding occurs
  real    :: sheddiamlg = 10.0e-03  ! diameter of hail to use fwmlarge
  real    :: sheddiam0  = 20.0e-03  ! diameter of hail at which all water is shed

  integer :: ifwmhopt = 1 ! option for calculating maximum liquid fraction when fwmh and/or fwmhl is set to -1
                          ! 1 = maximum based on size of maximum mass diameter
                          ! 2 = integrate over spectrum for maximum liquid (experimental)

  real   , private :: fwms = 0.5 ! maximum liquid water fraction on snow
  real   , private :: fwmh = 0.5 ! maximum liquid water fraction on graupel
  real   , private :: fwmhl = 0.5 ! maximum liquid water fraction on hail
  real    :: fwmlarge = 0.2 ! maximum liquid water fraction on hail larger than sheddiam
  integer :: ifwmfall = 0   ! whether to interpolate toward rain fall speed for graupel and hail
                            ! when diam < sheddiam and liquid fraction is predicted (0=no, 1=yes)

  logical :: rescale_high_alpha = .false.  ! whether to rescale number. conc. when alpha = alphamax (3-moment only)
  logical :: rescale_low_alpha = .true.    ! whether to rescale Z (graupel/hail) when alpha = alphamin (3-moment only)
  logical :: rescale_low_alphar = .true.    ! whether to rescale Z for rain when alpha = alphamin (3-moment only)

  real, parameter :: alpharmax = 8. ! limited for rwvent calculation

  integer, private ::  ihlcnh = 1  ! which graupel -> hail conversion to use
                          ! 1 = Milbrandt and Yau (2005) using Ziegler 1985 wet growth diameter
                          ! 2 = Straka and Mansell (2005) conversion using size threshold
  real, private :: hlcnhdia = 1.e-3 ! threshold diameter for graupel -> hail conversion for ihlcnh = 1 option.
  real, private :: hlcnhqmin = 0.1e-3 ! minimum graupel mass content for graupel -> hail conversion (ihlcnh = 1)
  real   , private :: hldia1 = 20.0e-3  ! threshold diameter for graupel -> hail conversion for ihlcnh = 2 option.
  integer :: icvhl2h = 0   ! allow conversion of hail back to graupel when hail density gets close to minimum allowed

  integer, private :: imurain = 1 ! 3 for gamma-volume, 1 for gamma-diameter DSD for rain.
  integer, private :: imusnow = 3 ! 3 for gamma-volume, 1 for gamma-diameter DSD for snow (=1 NOT IMPLEMENTED!!).
  integer, private :: iturbenhance = 0 ! warm-rain collision enhancement
                              ! 1 = enhance autoconversion only
                              ! 2 = add rain collection of cloud
                              ! 3 = add rain self-collection
  integer, private :: isedonly = 0 ! 1= only do sedimentation and skip other microphysics
  integer, private :: iferwisventr = 2 ! =1 for Ferrier rwvent, =2 for Wisner rwvent (imurain=1)
  integer, private :: izwisventr   = 2 ! =1 for old Ziegler rwvent, =2 for Wisner-style rwvent (imurain=3)
  integer :: iresetmoments = 0 ! if >0, then set all moments to zero when one of them is zero (3-moment only)
  integer, private :: imaxdiaopt    = 3 ! = 1 use mean diameter for breakup
                               ! = 2 use maximum mass diameter for breakup
                               ! = 3 use mass-weighted diameter for breakup
  integer, private :: dmrauto       = 0 ! = -1 no limiter on crcnw
                              ! =  0 limit crcnw when qr > 1.2*L (Cohard-Pinty 2002)
                              ! =  1 DTD version based on MY code
                              ! =  2 DTD mass-weighted version based on MY code
                              ! =  3 Milbrandt version (from Cohard and Pinty code
  real, parameter :: alpharaut = 0.0 ! MY2005 for autoconversion
  real    :: cxmin = 1.e-4  ! threshold cutoff for number concentration
  real    :: zxmin = 1.e-28 ! threshold cutoff for reflectivity moment

  integer :: ithompsoncnoh = 0 ! For single moment graupel only
                           ! 0 = fixed intercept
                           ! 1 = intercept based on graupel mass

  integer :: ivhmltsoak = 1   ! 0=off, 1=on : flag to simulate soaking (graupel/hail) during melting 
                         ! when liquid fraction is not predicted
  integer, private :: ioldlimiter = 0 ! test switch for new(=0) or old(=1) size limiter at the end of GS for 3-moment categories
  integer, private :: isnowfall = 2   ! Option for choosing between snow fall speed parameters
                         ! 1 = original Zrnic et al. (Mansell et al. 2010)
                         ! 2 = Ferrier 1994 (results in slower fall speeds)

  integer, private :: isnowdens = 1   ! Option for choosing between snow density options
                             ! 1 = constant of 100 kg m^-3
                             ! 2 = Option based on Cox 

  integer, private  :: ibiggsnow   = 3 ! 1 = switch conversion over to snow for small frozen drops from Bigg freezing
                                       ! 2 = switch conversion over to snow for small frozen drops from rain-ice interaction
                                       ! 3 = switch conversion over to snow for small frozen drops from both

  integer, private  :: ixtaltype = 1 ! =1 column, =2 disk (similar to Takahashi)

  real, private  :: takshedsize1 = 0.15 ! diameter (cm) of drop shed from ice with D > 1.9 cm
  real, private  :: takshedsize2 = 0.3 ! diameter (cm) of drop shed from ice with D < 1.9 cm and D > 0.8 cm

  real, private     :: evapfac     = 1.0 ! Multiplier on rain evaporation rate
  real, private     :: depfac      = 1.0 ! Multiplier on graupel/hail deposition/sublimation rate

  integer, private :: ibinhmlr = 0  ! =1 use incomplete gammas to determine melting from larger and smaller sizes of graupel, and appropriate shed drop sizes 
                           ! =2 to test melting by temporary bins
  integer, private :: ibinhlmlr = 0  ! =1 use incomplete gammas to determine melting from larger and smaller sizes of hail, and appropriate shed drop sizes 
                           ! =2 to test melting by temporary bins
  real, private :: snowmeltdia = 0 ! If nonzero, sets the size of rain drops from melting snow.
  real, private :: delta_alphamlr = 0.5 ! offset from alphamax at which melting does not further collapse the shape parameter

  integer :: iqvsopt = 0 ! =0 use old default for tabqvs; =1 use Bolton formulation (Rogers and Yau)

  real    :: maxsupersat = 1.9 ! maximum supersaturation ratio, above which a saturation adustment is done

  integer :: ihailrain = 0 ! =0 default turns off graupel/hail rain collection for T < tfr
                           ! =1  turns on graupel/hail rain collection for T < tfr (old behavior)


  integer, parameter :: icespheres = 0 ! turn ice spheres (frozen droplets) on (1) or off (0). NOT COMPLETE IN WRF/ARPS/CM1 CODE!
  integer, parameter :: lqmx = 30
  integer, parameter :: lt = 1
  integer, parameter :: lv = 2
  integer, parameter :: lc = 3
  integer, parameter :: lr = 4
  integer, parameter :: li = 5
  integer, private :: lis = 0
  integer, private :: ls = 6
  integer, private :: lh = 7
  integer, private :: lhl = 0

  integer, private  :: lccn = 9 ! 0 or 9, other indices adjusted accordingly
  integer, private :: lccna = 0
  integer, private :: lcina = 0
  integer, private :: lcin = 0
  integer, private :: lnc = 9
  integer, private :: lnr = 10
  integer, private :: lni = 11
  integer, private :: lnis = 0
  integer, private :: lns = 12
  integer, private :: lnh = 13
  integer, private :: lnhl = 0
  integer, private :: lss = 0
  integer :: lvh = 15

  integer, private :: lhab = 8
  integer, private :: lg = 7

! Particle volume

  integer :: lvi = 0
  integer :: lvs = 0
  integer :: lvgl = 0
  integer :: lvgm = 0
  integer :: lvgh = 0
  integer :: lvf = 0
!  integer :: lvh = 16
  integer :: lvhl = 0

! liquid water fraction (not predicted here but tested for)
  integer :: lhw = 0
  integer :: lsw = 0
  integer :: lhlw = 0

! reflectivity (6th moment) ! not predicted here but may be tested against

  integer :: lzr = 0
  integer :: lzi = 0
  integer :: lzs = 0
  integer :: lzgl = 0
  integer :: lzgm = 0
  integer :: lzgh = 0
  integer :: lzf = 0
  integer :: lzh = 0
  integer :: lzhl = 0

! Space charge

  integer :: lscw = 0
  integer :: lscr = 0
  integer :: lsci = 0
  integer :: lscis = 0
  integer :: lscs = 0
  integer :: lsch = 0
  integer :: lschl = 0
  integer :: lscwi = 0
  integer :: lscpi = 0
  integer :: lscni = 0
  integer :: lscpli = 0
  integer :: lscnli = 0
  integer :: lschab = 0

  integer :: lscb = 0
  integer :: lsce = 0
  integer :: lsceq = 0

!  integer, parameter :: lscmx = 100

  integer :: lne = 0 ! last varible for transforming

  real :: cnoh0 = 4.0e+5
  real :: hwdn1 = 700.0

  real    :: alphai  = 0.0 ! shape parameter for ZIEG ice crystals ! not currently used
  real    :: alphas  = 0.0 ! shape parameter for ZIEG snow         ! used only for single moment
  real    :: alphar  = 0.0 ! shape parameter for rain (imurain=1 only)
  real, private    :: alphah  = 0.0 ! set in namelist!! shape parameter for ZIEG graupel
  real, private    :: alphahl = 1.0 ! set in namelist!! shape parameter for ZIEG hail

  real    :: dmuh    = 1.0  ! power in exponential part (graupel)
  real    :: dmuhl   = 1.0  ! power in exponential part (hail)

  real, parameter :: alphamax = 15.
  real, parameter :: alphamin = 0.
  real, parameter :: rnumin = -0.8
  real, parameter :: rnumax = 15.0


  real            :: cnu = 0.0
  real, parameter :: rnu = -0.8, snu = -0.8, cinu = 0.0
!      parameter ( cnu = 0.0, rnu = -0.8, snu = -0.8, cinu = 0.0 )

  real xnu(lc:lqmx) ! 1st shape parameter (mass)
  real xmu(lc:lqmx) ! 2nd shape parameter (mass)
  real dnu(lc:lqmx) ! 1st shape parameter (diameter)
  real dmu(lc:lqmx) ! 2nd shape parameter (diameter)

  real ax(lc:lqmx)
  real bx(lc:lqmx)
  real fx(lc:lqmx)

      real da0 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real dab0(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real dab1(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real da1 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real bb  (lc:lqmx)

! put ipelec here for now....
  integer :: ipelec = 0
  integer :: isaund = 0
  logical :: idonic = .false.
  integer, private :: elec_on_time = -1     ! time (seconds) to turn on charge separation.
  integer, private :: elec_ramp_time = 0   ! time (interval) for linear ramp after elec_on_time 
                                   ! (i.e., linear factor on chg sep to smoothly turn on elec)
                                   ! full charging rate is achieved at time = elec_on_time + elec_ramp_time
  integer :: jchgs = 3  ! number of points near boundary where charging is turned off (to keep lightning from getting wonky)
  integer :: jchgn = 2
  integer :: ichge = 3
  integer :: ichgw = 2
  real    :: charging_border = 4000. ! width of no-charging zone from boundary
  integer :: nic_noliq = 0  ! switch for NIC without riming (MST, 2006)
  real    :: qc_noliq = 0.1 ! charge per collision (in fC) for non-riming NIC.
  logical :: largeion = .false.   ! true: use large ion category 
  integer :: nonigrd = 0
  integer :: ftauopt = 1
  real    :: rgard1  = 1.0
  real    :: rarfac  = 1.0
  real    :: trever = -15.0
  real    :: cidiamin = 0.0  ! minimum crystal size for charging

  real    :: erbnd  = 0.01
  real    :: fdgt   = 1.0
  real    :: costhe = 0.35
  real    :: delqnsa = -50.0e-14
  real    :: delqxsa =  50.0e-14
  real    :: delqnsb = -2.0e-14
  real    :: delqxsb =  2.0e-14
  real    :: delqnia = -50.0e-14
  real    :: delqxia =  50.0e-14
  real    :: scxacymax = 3000.e-12 ! Maximum charging rate magnitude (C s-1 m-3)

!
! Takahashi lookup table
!
      integer, parameter :: nlwc=30, ntem=31
      real takalu(0:ntem,0:nlwc)

      real, private    :: delqnw = -1.0e-10!-1.0e-12 !
      real, private    :: delqxw =  1.0e-10! 1.0e-12 !
      real :: tindmn = 233, tindmx = 298.0  ! min and max temperatures where inductive charging is allowed
!
!  gamma function lookup table
!
      integer ngm0,ngm1,ngm2
      parameter (ngm0=3001,ngm1=500,ngm2=500)
      double precision, parameter :: dgam = 0.01, dgami = 100.
      double precision gmoi(0:ngm0) ! ,gmod(0:ngm1,0:ngm2),gmdi(0:ngm1,0:ngm2)

      integer, parameter :: nqiacralpha =  240 !480 ! 240 ! 120 ! 15
      integer, parameter :: nqiacrratio =  100 ! 500 !50  ! 25
      real,    parameter :: maxratiolu = 25.
      real,    parameter :: maxalphalu = 15.
      real,    parameter :: dqiacralpha = maxalphalu/Float(nqiacralpha), dqiacrratio = maxratiolu/Float(nqiacrratio)
      real,    parameter :: dqiacrratioinv = 1./dqiacrratio, dqiacralphainv = 1./dqiacralpha
      real :: ciacrratio(0:nqiacrratio,0:nqiacralpha)
      real :: qiacrratio(0:nqiacrratio,0:nqiacralpha)
      real :: ziacrratio(0:nqiacrratio,0:nqiacralpha)
      double precision :: gamxinflu(0:nqiacrratio,0:nqiacralpha,10,2) ! last index for graupel (1) or hail (2)

    integer, parameter :: ngdnmm = 9
    real :: mmgraupvt(ngdnmm,3)  ! Milbrandt and Morrison (2013) fall speed coefficients for graupel/hail

    DATA mmgraupvt(:,1) / 50., 150., 250., 350., 450., 550., 650., 750., 850./
    DATA mmgraupvt(:,2) / 62.923, 94.122, 114.74, 131.21, 145.26, 157.71, 168.98, 179.36, 189.02 /
    DATA mmgraupvt(:,3) / 0.67819, 0.63789, 0.62197, 0.61240, 0.60572, 0.60066, 0.59663, 0.59330, 0.59048 /

      integer lsc(lc:lqmx)
      integer ln(lc:lqmx)
      integer ipc(lc:lqmx)
      integer lvol(lc:lqmx)
      integer lz(lc:lqmx)
      integer lliq(li:lqmx)
      integer denscale(lc:lqmx) ! flag for density scaling (mixing ratio conversion)

      integer ido(lc:lqmx)
      logical ldovol

      real xdn0(lc:lqmx)     !bulk density for different hydrometeor [kg/m3]
      real xdnmx(lc:lqmx), xdnmn(lc:lqmx)
      real cdx(lc:lqmx)
      real cno(lc:lqmx)
      real xvmn(lc:lqmx), xvmx(lc:lqmx)
      real qxmin(lc:lqmx)

      integer nqsat
      parameter (nqsat=1000001) ! (nqsat=20001)
      real fqsat,fqsati
      parameter (fqsat=0.002,fqsati=1./fqsat)
      real tabqvs(nqsat),tabqis(nqsat),dtabqvs(nqsat),dtabqis(nqsat)

!
!  constants
!
      real, parameter :: cp608 = 0.608          ! constant used in conversion of T to Tv
      real, parameter :: ar = 841.99666         ! rain terminal velocity power law coefficient (LFO)
      real, parameter :: br = 0.8               ! rain terminal velocity power law coefficient (LFO)
      real, parameter :: aradcw = -0.27544      !
      real, parameter :: bradcw = 0.26249e+06   !
      real, parameter :: cradcw = -1.8896e+10   !
      real, parameter :: dradcw = 4.4626e+14    !
      real, parameter :: bta1 = 0.6             ! beta-1 constant used for ice nucleation by deposition (Ferrier 94, among others)
      real, parameter :: cnit = 1.0e-02         ! No for ice nucleation by deposition (Cotton et al. 86)
      real, parameter :: dragh = 0.60           ! coefficient used to adjust fall speed for hail versus graupel (Pruppacher and Klett 78)
      real, parameter :: dnz00 = 1.225          ! reference/MSL air density
      real, parameter :: rho00 = 1.225          ! reference/MSL air density
!      cs = 4.83607122       ! snow terminal velocity power law coefficient (LFO)
!      ds = 0.25             ! snow terminal velocity power law coefficient (LFO)
!  new values for  cs and ds
      real, parameter :: cs = 12.42             ! snow terminal velocity power law coefficient 
      real, parameter :: ds = 0.42              ! snow terminal velocity power law coefficient 
      real, parameter :: pi = 3.141592653589793
      real, parameter :: piinv = 1./pi
      real, parameter :: pid4 = pi/4.0

      real, parameter :: gr = 9.8

!
! max and min mean volumes
!
      real xvrmn, xvrmx0  ! min, max rain volumes
      real xvsmn, xvsmx  ! min, max snow volumes
      real xvfmn, xvfmx  ! min, max frozen drop volumes
      real xvgmn, xvgmx  ! min, max graupel volumes
      real xvhmn, xvhmn0, xvhmx, xvhmx0  ! min, max hail volumes
      real xvhlmn, xvhlmx  ! min, max lg hail volumes

      real, parameter :: dhlmn = 0.3e-3, dhlmx = 40.e-3
      real, parameter :: dhmn0 = 0.3e-3
      real, private :: dhmn = dhmn0, dhmx = -1.

      real, parameter :: cwradn = 2.5e-6, xcradmn = cwradn    ! minimum radius
      real, parameter :: cwradx = 60.e-6, xcradmx = cwradx    ! maximum radius
      real, parameter :: cwc1 = 6.0/(pi*1000.)

!      parameter( xvcmn=4.188e-18 )   ! mks  min volume = 3 micron radius
      real, parameter :: xvcmn=0.523599*(2.*cwradn)**3    ! mks  min volume = 2.5 micron radius
      real, parameter :: xvcmx=0.523599*(2.*xcradmx)**3    ! mks  min volume = 2.5 micron radius
      real, parameter :: cwmasn = 1000.*xvcmn   ! minimum mass, defined by radius of 5.0e-6
      real, parameter :: cwmasx = 1000.*xvcmx   ! maximum mass, defined by radius of 50.0e-6
      real, parameter :: cwmasn5 = 1000.*0.523599*(2.*5.0e-6)**3 !  5.23e-13

      real, parameter :: xvimn=0.523599*(2.*5.e-6)**3    ! mks  min volume = 5 micron radius
      real, parameter :: xvimx=0.523599*(2.*1.e-3)**3    ! mks  max volume = 1 mm radius (solid sphere approx)

      real, private   :: xvdmx = -1.0 ! 3.0e-3
      real     :: xvrmx
      parameter( xvrmn=0.523599*(80.e-6)**3, xvrmx0=0.523599*(6.e-3)**3 ) !( was 4.1887e-9 )  ! mks
      parameter( xvsmn=0.523599*(0.01e-3)**3, xvsmx=0.523599*(10.e-3)**3 ) !( was 4.1887e-9 )  ! mks
      parameter( xvfmn=0.523599*(0.1e-3)**3, xvfmx=0.523599*(10.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      parameter( xvgmn=0.523599*(0.1e-3)**3, xvgmx=0.523599*(10.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      parameter( xvhmn0=0.523599*(0.3e-3)**3, xvhmx0=0.523599*(20.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      parameter( xvhlmn=0.523599*(dhlmn)**3, xvhlmx=0.523599*(dhlmx)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3

!
!  electrical permitivity of air C / (N m**2) -  check the units
!
      real eperao
      parameter (eperao  = 8.8592e-12 )

      real ec,eci  ! fundamental unit of charge
      parameter (ec = 1.602e-19)
      parameter (eci = 1.0/ec)

      real    :: scwppmx = 20.0e-12
      real    :: scippmx = 20.0e-12
!
!  constants
!
      real, parameter :: c1f3 = 1.0/3.0

      real, parameter :: cai = 21.87455
      real, parameter :: caw = 17.2693882
      real, parameter :: cbi = 7.66
      real, parameter :: cbw = 35.86

      real, parameter :: cbwbolton = 29.65 ! constants for Bolton formulation
      real, parameter :: cawbolton = 17.67

      real, parameter :: tfr = 273.15, tfrh = 233.15

      real, parameter :: cp = 1004.0, rd = 287.04
      real, parameter :: cpi = 1./cp
      real, parameter :: cap = rd/cp, poo = 1.0e+05

      real, parameter :: rw = 461.5              ! gas const. for water vapor
      real, parameter :: advisc0 = 1.832e-05     ! reference dynamic viscosity (SMT; see Beard & Pruppacher 71)
      real, parameter :: advisc1 = 1.718e-05     ! dynamic viscosity constant used in thermal conductivity calc
      real, parameter :: tka0 = 2.43e-02         ! reference thermal conductivity
      real, parameter :: tfrcbw = tfr - cbw
      real, parameter :: tfrcbi = tfr - cbi

     ! GHB: Needed for eqtset=2 in cm1
!     REAL, PRIVATE ::      cv = cp - rd
     real, private, parameter ::      cv = 717.0             ! specific heat at constant volume - air
     REAL, PRIVATE, parameter ::      cvv = 1408.5
     REAL, PRIVATE, parameter ::      cpl = 4190.0
     REAL, PRIVATE, parameter ::      cpigb = 2106.0
     ! GHB

      real, parameter ::  bfnu0 = (rnu + 2.0)/(rnu + 1.0)
      real :: ventr, ventrn, ventc, c1sw

      real :: cckm,ccne,ccnefac,cnexp,CCNE0

      integer :: na = 9
      integer :: nxtra = 1
      real gf4p5, gf4ds, gf4br
      real gsnow1, gsnow53, gsnow73
      real gfcinu1, gfcinu1p47, gfcinu2p47

      real :: cwchtmp0 = 1.0
      real :: cwchltmp0 = 1.0

      real    :: esctot = 1.0e-13

      integer iexy(lc:lqmx,lc:lqmx)
      integer :: ieswi = 1,  ieswc = 1, ieswr = 0
      integer :: iehlsw = 1, iehli = 1,  iehlc = 1, iehlr = 0
      integer :: iehwsw = 1, iehwi = 1,  iehwc = 1, iehwr = 0

      logical, parameter :: do_satadj_for_wrfchem = .true.


 contains

!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 
!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 

 subroutine re_nssl_volume( q_c,q_i,q_s, n_c,n_i,n_s, mu_c,mu_i,mu_s , &
                            rho_c,rho_i,rho_s, re_c,re_i,re_s )

!c: cloud droplet, i: cloud ice, s: snow aggregate
 real,intent(inout) :: q_c, q_i, q_s    ! mass [g/m3]
 real,intent(inout) :: n_c, n_i, n_s    ! number conc. [g/m3] and [#/m3]
 real,intent(in) :: mu_c, mu_i, mu_s ! shape parameters [-]
 real,intent(in) :: rho_c, rho_i, rho_s ! effective density [kg/m3]

 real,intent(out) :: re_c, re_i, re_s ! effective radius [um]

! local parameters
 real :: gamc1,gamc2,gami1, gami2,gams1,gams2, factor_c,factor_i,factor_s  !gamma parameters
 real :: lam_c, lam_i, lam_s  !slope parameters [1/m]

 real, parameter :: q_min_condensate = 1.0e-5  ! siginificant minimum condensate mixing ratio [g/m3]

!
! preparation 
!
 gamc1 = Gamma_sp(2. + mu_c)    ! cloud droplet
 gamc2 = Gamma_sp(1  + mu_c)
 gami1 = Gamma_sp(2. + mu_i)    ! cloud ice
 gami2 = Gamma_sp(1  + mu_i)
 gams1 = Gamma_sp(2. + mu_s)    ! snow
 gams2 = Gamma_sp(1. + mu_s)     

 factor_c = (1. + mu_c)*Gamma_sp(1. + mu_c)/Gamma_sp(5./3. + mu_c)
 factor_i = (1. + mu_i)*Gamma_sp(1. + mu_i)/Gamma_sp(5./3. + mu_i)
 factor_s = (1. + mu_s)*Gamma_sp(1. + mu_s)/Gamma_sp(5./3. + mu_s)


!
! cloud droplets 
!
 rec_if: if( q_c > q_min_condensate ) then
    lam_c = ((n_c*(Pi/6.)*rho_c*gamc1)/(0.001*q_c*gamc2))**(1./3.)  !slope parameter [1/m]
    if(lam_c <= const_smallest_real ) then
       re_c = 0. ; q_c  = 0. ; n_c  = 0.
    else
       re_c =   1.e6 * 0.5*factor_c/lam_c ! effective radius [micron] 
    endif
!print*, lam_c, q_c, n_c, const_smallest_real
!    if( 1.e6*0.5*factor_c/max(const_smallest_real,lam_c) > const_largest_real ) print*, 'overflow',factor_c, lam_c
 else
    re_c = 0.
    q_c  = 0.
    n_c  = 0.
 endif rec_if


!
! cloud ice
!
 rei_if: if( q_i > q_min_condensate ) then
    lam_i = ((n_i*(Pi/6.)*rho_i*gamc1)/(0.001*q_i*gamc2))**(1./3.)  !slope parameter [1/m]
    if( lam_i <= const_smallest_real ) then
       re_i = 0. ; q_i  = 0. ; n_i  = 0.
    else
       re_i =  1.e6 *  0.5*factor_i/lam_i  ! effective radius [micron] 
    endif
 else
    re_i = 0.
    q_i  = 0.
    n_i  = 0.
 endif rei_if


!
! snow aggregate
!
 res_if: if( q_s > q_min_condensate ) then
    lam_s = ((n_s*(Pi/6.)*rho_s*gamc1)/(0.001*q_s*gamc2))**(1./3.)  !slope parameter [1/m]
    if( lam_s <= const_smallest_real ) then
       re_s = 0. ; q_s  = 0. ; n_s  = 0.
    else
       re_s =   1.e6 * 0.5*factor_s/lam_s   ! effective radius [micron] 
    endif
 else
    re_s = 0.
    q_s  = 0.
    n_s  = 0.
 endif res_if


 return

 end subroutine re_nssl_volume

!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 
!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 

 SUBROUTINE setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,cnostmp, &
                      xmas,vtxbar,xdn,xvmn0,xvmx0,xv,cdx,            &
                      ipconc1,ndebug1,ngs,nz,kgs,fadvisc,   &
                      cwmasn,cwmasx,cwradn,cnina,cimna,cimxa,      &
                      itype1a,itype2a,temcg,infdo,alpha,ildo,axh,bxh,axhl,bxhl)


 implicit none

  integer ngscnt
  integer :: ngs0,ngs,nz
!      integer infall    ! whether to calculate number-weighted fall speeds

      real xv(ngs,lc:lhab)
      real qx(ngs,lv:lhab)       ! mass conc. [kg/m3]
      real qxw(ngs,ls:lhab)
      real cx(ngs,lc:lhab)       ! number conc. [#/m3]
      real vtxbar(ngs,lc:lhab,3)
      real xmas(ngs,lc:lhab)
      real xdn(ngs,lc:lhab)
      real xdia(ngs,lc:lhab,3)
      real xvmn0(lc:lhab), xvmx0(lc:lhab)
      real qxmin(lc:lhab)
      real cdx(lc:lhab)
      real alpha(ngs,lc:lhab)  !shape parameter of gamma DSD (using diameter)

      real rho0(ngs),rhovt(ngs),temcg(ngs)
      real cno(lc:lhab)
      real cnostmp(ngs)

      real cwc1, cimna, cimxa
      real cnina(ngs)
      integer kgs(ngs)
      real fadvisc(ngs)
      real fsw

      integer ipconc1
      integer ndebug1

      integer, intent (in) :: itype1a,itype2a,infdo
      integer, intent (in) :: ildo ! which species to do, or all if ildo=0

      real :: axh(ngs),bxh(ngs)
      real :: axhl(ngs),bxhl(ngs)

! Local vars



      real swmasmx, dtmp
      real cd
      real cwc0 ! ,cwc1
      real :: cwch(ngscnt), cwchl(ngscnt)
      real :: cwchtmp,cwchltmp,xnutmp
      real pii
      real cimasx,cimasn
      real cwmasn,cwmasx,cwradn
      real cwrad
      real vr,rnux
      real alp

      real ccimx

      integer mgs

      real arx,frx,vtrain,fw
      real fwlo,fwhi,rfwdiff
      real ar,br,cs,ds
!      real gf4p5, gf4ds, gf4br, ifirst, gf1ds
!      real gfcinu1, gfcinu1p47, gfcinu2p47
      real gr
      real rwrad,rwdia
      real mwfac
      integer il

!      save gf4p5, gf4ds, gf4br, ifirst, gf1ds
!      save gfcinu1, gfcinu1p47, gfcinu2p47
!      data ifirst /0/

      real bta1,cnit
      parameter ( bta1 = 0.6, cnit = 1.0e-02 )
      real x,y,tmp,del
      real aax,bbx,delrho
      integer :: indxr
      real mwt, nwt
      real, parameter :: rho00 = 1.225
      integer i
      real xvbarmax

      integer l1, l2


!
! set values
!
!      cwmasn = 5.23e-13  ! radius of 5.0e-6
!      cwradn = 5.0e-6
!      cwmasx = 5.25e-10  ! radius of 50.0e-6

      fwlo = 0.2                ! water fraction to start weighting toward rain fall speed
      fwhi = 0.4                ! water fraction at which rain fall speed only is used
      rfwdiff = 1./(fwhi - fwlo)

!      pi = 4.0*atan(1.0)
      pii = piinv ! 1.0/pi

      arx = 10.
      frx = 516.575 ! raind fit parameters for arx*(1 - Exp(-fx*d)), where d is rain diameter in meters.

      ar = 841.99666
      br = 0.8
      gr = 9.8
!  new values for  cs and ds
      cs = 12.42
      ds = 0.42

      IF ( ildo == 0 ) THEN
        l1 = lc
        l2 = lhab
      ELSE
        l1 = ildo
        l2 = ildo
      ENDIF


        IF ( lh  .gt. 1 ) THEN
          IF ( dmuh == 1.0 ) THEN
            cwchtmp = ((3. + dnu(lh))*(2. + dnu(lh))*(1.0 + dnu(lh)))**(-1./3.)
          ELSE
            cwchtmp = 6.0*pii*gamma_sp( (xnu(lh) + 1.)/xmu(lh) )/gamma_sp( (xnu(lh) + 2.)/xmu(lh) )
          ENDIF
        ENDIF
        IF ( lhl .gt. 1 ) THEN
          IF ( dmuhl == 1.0 ) THEN
            cwchltmp = ((3. + dnu(lhl))*(2. + dnu(lhl))*(1.0 + dnu(lhl)))**(-1./3.)
          ELSE
            cwchltmp = 6.0*pii*gamma_sp( (xnu(lhl) + 1)/xmu(lhl) )/gamma_sp( (xnu(lhl) + 2)/xmu(lhl) )
          ENDIF
        ENDIF

        IF ( ipconc .le. 5 ) THEN  ! ipconc = 5

          IF ( lh  .gt. 1 ) cwch(:) =  cwchtmp
          IF ( lhl .gt. 1 ) cwchl(:) = cwchltmp

        ELSE

          DO mgs = 1,ngscnt

          IF ( lh  .gt. 1 .and. ( ildo == 0 .or. ildo == lh ) ) THEN
           IF ( qx(mgs,lh) .gt. qxmin(lh) ) THEN
            IF ( dmuh == 1.0 ) THEN
              cwch(mgs) = ((3. + alpha(mgs,lh))*(2. + alpha(mgs,lh))*(1.0 + alpha(mgs,lh)))**(-1./3.)
             ELSE
             xnutmp = (alpha(mgs,lh) - 2.0)/3.0
             cwch(mgs) =  6.0*pii*gamma_sp( (xnutmp + 1.)/xmu(lh) )/gamma_sp( (xnutmp + 2.)/xmu(lh) )
            ENDIF
           ELSE
             cwch(mgs) = cwchtmp
           ENDIF
          ENDIF
          IF ( lhl .gt. 1 .and. ( ildo == 0 .or. ildo == lhl ) ) THEN
           IF ( qx(mgs,lhl) .gt. qxmin(lhl) ) THEN
            IF ( dmuhl == 1.0 ) THEN
              cwchl(mgs) = ((3. + alpha(mgs,lhl))*(2. + alpha(mgs,lhl))*(1.0 + alpha(mgs,lhl)))**(-1./3.)
             ELSE
             xnutmp = (alpha(mgs,lhl) - 2.0)/3.0
             cwchl(mgs) = 6.0*pii*gamma_sp( (xnutmp + 1)/xmu(lhl) )/gamma_sp( (xnutmp + 2)/xmu(lhl) )
            ENDIF
           ELSE
             cwchl(mgs) = cwchltmp
           ENDIF
          ENDIF

          ENDDO

        ENDIF


      cimasn = Min( cimas0, 6.88e-13)
      cimasx = 1.0e-8
      ccimx = 5000.0e3   ! max of 5000 per liter

      cwc1 = 6.0/(pi*1000.)
      cwc0 = pii ! 6.0*pii
      mwfac = 6.0**(1./3.)


      if (ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set scale diameter'
!


!
!  cloud water variables
! ################################################################
!
!  DROPLETS
!
!
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set cloud water variables'

 IF ( ildo == 0 .or. ildo == lc ) THEN

   do mgs = 1,ngscnt   !mgs is location index  (ngscnt is maximum location idex)
      xv(mgs,lc) = 0.0
      xdn(mgs,lc) = 1000.0 ! liquid density [kg/m3] toshii

    IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN 

        xmas(mgs,lc) = min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )  ! mean particle mass [kg]
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)                                     ! mean particle volume [m3] (assuming xdn is mean density)

 
      xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**(1./3.)  ! d0: cwc1 = 6.0/(pi*1000.)  [m]
      xdia(mgs,lc,2) = xdia(mgs,lc,1)**2             ! d0^2                       [m2]
      xdia(mgs,lc,3) = xdia(mgs,lc,1)                ! dbar = 1./lambda           [m]  DSD params
      cwrad = 0.5*xdia(mgs,lc,1)

      IF ( fadvisc(mgs) > 0.0 ) THEN
        vtxbar(mgs,lc,1) =  (2.0*gr*xdn(mgs,lc) *(cwrad**2)) /(9.0*fadvisc(mgs))  ! terminal velocity [m/s] ??
      ELSE
        vtxbar(mgs,lc,1) = 0.0
      ENDIF


    ELSE
       xmas(mgs,lc) = cwmasn
       IF ( ipconc .le. 1 ) cx(mgs,lc) = 0.01
       xdia(mgs,lc,1) = 2.*cwradn
       xdia(mgs,lc,2) = 4.*cwradn**2
       xdia(mgs,lc,3) = xdia(mgs,lc,1)
       vtxbar(mgs,lc,1) = 0.0

   ENDIF !} qcw .gt. qxmin(lc)

  end do

 ENDIF



!
! cloud ice variables
! columns
!
! ################################################################
!
!  CLOUD ICE
!
 if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set cip'

 IF ( li .gt. 1 .and. ( ildo == 0 .or. ildo == li ) ) THEN

   do mgs = 1,ngscnt
      xdn(mgs,li)  = 900.0  !ice deinsity [kg/m3]


    IF ( qx(mgs,li) .gt. qxmin(li) ) THEN

      cx(mgs,li) = Max(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasx)
      cx(mgs,li) = Min(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasn)

      xmas(mgs,li) = max( qx(mgs,li)*rho0(mgs)/cx(mgs,li), cimasn )  !  mean particle mass [kg]

      if ( xmas(mgs,li) .gt. 0.0 ) THEN !
       IF ( ixtaltype == 1 ) THEN ! column  ixtaltype = 1 now
        xdia(mgs,li,1) = 0.1871*(xmas(mgs,li)**(0.3429))     ! d0   [m]  
        xdia(mgs,li,3) = 0.1871*(xmas(mgs,li)**(0.3429))     ! dbar = 1/lambda [m]  --> DSD param
       ELSEIF  ( ixtaltype == 2 ) THEN ! disk
        xdia(mgs,li,1) = 0.277823*xmas(mgs,li)**0.359971
        xdia(mgs,li,3) = 0.277823*xmas(mgs,li)**0.359971
       ENDIF
      end if

      xv(mgs,li) = xmas(mgs,li)/xdn(mgs,li)   ! mean particle volume [m3]

      IF ( ixtaltype == 1 ) THEN ! column
        tmp = (67056.6300748612*rhovt(mgs))/  &
        (((1.0 + cinu)/xv(mgs,li))**0.4716666666666667*gfcinu1)
        vtxbar(mgs,li,2) = tmp*gfcinu1p47
        vtxbar(mgs,li,1) = tmp*gfcinu2p47/(1. + cinu)
        vtxbar(mgs,li,3) = vtxbar(mgs,li,1)
        ELSEIF  ( ixtaltype == 2 ) THEN ! disk -- but just use column fall speed for now
        tmp = (67056.6300748612*rhovt(mgs))/  &
        (((1.0 + cinu)/xv(mgs,li))**0.4716666666666667*gfcinu1)
        vtxbar(mgs,li,2) = tmp*gfcinu1p47
        vtxbar(mgs,li,1) = tmp*gfcinu2p47/(1. + cinu)
        vtxbar(mgs,li,3) = vtxbar(mgs,li,1)

      ENDIF

      xdia(mgs,li,2) = xdia(mgs,li,1)**2   !d0^2 [m2]

    ELSE
       xmas(mgs,li) = 1.e-13
       xdn(mgs,li)  = 900.0
       xdia(mgs,li,1) = 1.e-7
       xdia(mgs,li,2) = (1.e-14)
       xdia(mgs,li,3) = 1.e-7
       vtxbar(mgs,li,1) = 0.0
    ENDIF

  end do

 ENDIF ! li .gt. 1


! ################################################################
!
!  RAIN:  alpha is shape parameter (alpha = alphar = 0)
!

!
 IF ( ildo == 0 .or. ildo == lr ) THEN

  do mgs = 1,ngscnt
     xdn(mgs,lr) = 1000.0 !liquid density [kg/m3]

     if ( qx(mgs,lr) .gt. qxmin(lr) ) then

        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))  ! mean particle volume [m3]??
        xvbarmax = xvmx(lr)/((4. + alpha(mgs,lr))**3/((3. + alpha(mgs,lr))*(2. + alpha(mgs,lr))*(1. + alpha(mgs,lr))))  ! maximum xv??

        IF ( xv(mgs,lr) .gt. xvbarmax ) THEN  !set up mim/max values
          xv(mgs,lr) = xvbarmax
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvbarmax*xdn(mgs,lr))  !why cx (number conc.) is calculated here? already input??
        ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        ENDIF

        xmas(mgs,lr) = xv(mgs,lr)*xdn(mgs,lr)         ! mean particle mass [kg]
        xdia(mgs,lr,3) = (xmas(mgs,lr)*cwc1)**(1./3.) ! dbar = 1/lambda [m]  
        xdia(mgs,lr,1) = (6.*pii*xv(mgs,lr)/((alpha(mgs,lr)+3.)*(alpha(mgs,lr)+2.)*(alpha(mgs,lr)+1.)))**(1./3.)  ! d0 [m] DSD params

    else

        xdia(mgs,lr,1) = 1.e-9
        xdia(mgs,lr,3) = 1.e-9
        xmas(mgs,lr) = xdn(mgs,lr)*(pi/6.)*xdia(mgs,lr,1)**3

    end if !qxmin
    xdia(mgs,lr,2) = xdia(mgs,lr,1)**2

  end do

 ENDIF
! ################################################################
!
!  SNOW
!

 IF ( ls .gt. 1 .and. ( ildo == 0 .or. ildo == ls ) ) THEN

  do mgs = 1,ngscnt

     if ( qx(mgs,ls) .gt. qxmin(ls) ) then

        xdn(mgs,ls) = 100. ! snow density [kg/m3]

        xmas(mgs,ls) =  rho0(mgs)*qx(mgs,ls)/(Max(1.0e-9,cx(mgs,ls)))

        xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*Max(1.0e-9,cx(mgs,ls)))
        xdia(mgs,ls,1) = (xv(mgs,ls)*cwc0*6.0)**(1./3.)  ! cwc0 = pii

        IF ( xv(mgs,ls) .lt. xvmn(ls) .and. isnowdens == 1) THEN
          xv(mgs,ls) = Max( xvmn(ls),xv(mgs,ls) )
          xmas(mgs,ls) = xv(mgs,ls)*xdn(mgs,ls)
          cx(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xmas(mgs,ls))
          xdia(mgs,ls,1) = (xv(mgs,ls)*cwc0*6.0)**(1./3.)
        ENDIF

        IF ( xv(mgs,ls) .gt. xvmx(ls)*Max(1.,100./Min(100.,xdn(mgs,ls))) ) THEN
          xv(mgs,ls) = Min( xvmx(ls), Max( xvmn(ls),xv(mgs,ls) ) )
          xmas(mgs,ls) = 0.106214*xv(mgs,ls)**(2./3.)
          cx(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xmas(mgs,ls))
          xdn(mgs,ls) = 0.0346159*Sqrt(cx(mgs,ls)/(qx(mgs,ls)*rho0(mgs)) )
          xdia(mgs,ls,1) = Sqrt( xmas(mgs,ls)/0.069 )
        ENDIF
        
        ! re-calculate after max/min setting.
        xdia(mgs,ls,1) = (xv(mgs,ls)*cwc0*6.0)**(1./3.)  ! d0 [m]
        xdia(mgs,ls,3) = xdia(mgs,ls,1)                  ! dbar = 1/lambda [m]  --> DSD param


      else

        xdia(mgs,ls,1) = 1.e-9
        xdia(mgs,ls,3) = 1.e-9
        cx(mgs,ls) = 0.0

        IF ( isnowdens == 2 ) xdn(mgs,ls) = 90. ! Set values according to Cox relationship

      end if ! qxmin

      xdia(mgs,ls,2) = xdia(mgs,ls,1)**2

  end do

 ENDIF ! ls .gt 1

!
!
! ################################################################
!
!  GRAUPEL
!

 IF ( lh .gt. 1 .and. ( ildo == 0 .or. ildo == lh ) ) THEN

   do mgs = 1,ngscnt
      !xdn will come from qg / qvg

      if ( qx(mgs,lh) .gt. qxmin(lh) ) then

        xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*Max(1.0e-9,cx(mgs,lh)))
        xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)          !mean particle mass

        IF ( xv(mgs,lh) .lt. xvmn(lh) .or. xv(mgs,lh) .gt. xvmx(lh) ) THEN
          xv(mgs,lh) = Min( xvmx(lh), Max( xvmn(lh),xv(mgs,lh) ) )
          xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)
          cx(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xmas(mgs,lh))
        ENDIF

         xdia(mgs,lh,3) = (xv(mgs,lh)*6.*pii)**(1./3.) ! dbar = 1/lambda [m]

         IF ( dmuh == 1.0 ) THEN  !choice here
           xdia(mgs,lh,1) = cwch(mgs)*xdia(mgs,lh,3)    ! d0 [m] DSD params
         ELSE
           xdia(mgs,lh,1) = (xv(mgs,lh)*cwch(mgs))**(1./3.)
         ENDIF

      else
         xdia(mgs,lh,1) = 1.e-9
         xdia(mgs,lh,3) = 1.e-9
      end if

      xdia(mgs,lh,2) = xdia(mgs,lh,1)**2

   end do

 ENDIF

!
! ################################################################
!
!  HAIL
!

      IF ( lhl .gt. 1 .and. ( ildo == 0 .or. ildo == lhl ) ) THEN

      do mgs = 1,ngscnt
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*Max(1.0e-9,cx(mgs,lhl)))
        xmas(mgs,lhl) = xv(mgs,lhl)*xdn(mgs,lhl)
!        write(0,*) 'setvt: xv = ',xv(mgs,lhl),xdn(mgs,lhl),cx(mgs,lhl),xmas(mgs,lhl),qx(mgs,lhl)

        IF ( xv(mgs,lhl) .lt. xvmn(lhl) .or. xv(mgs,lhl) .gt. xvmx(lhl) ) THEN
          xv(mgs,lhl) = Min( xvmx(lhl), Max( xvmn(lhl),xv(mgs,lhl) ) )
          xmas(mgs,lhl) = xv(mgs,lhl)*xdn(mgs,lhl)
          cx(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xmas(mgs,lhl))
        ENDIF

        xdia(mgs,lhl,3) = (xv(mgs,lhl)*6./pi)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
         IF ( dmuhl == 1.0 ) THEN
           xdia(mgs,lhl,1) = cwchl(mgs)*xdia(mgs,lhl,3)
         ELSE
           xdia(mgs,lhl,1) = (xv(mgs,lhl)*cwchl(mgs))**(1./3.)
         ENDIF

!        write(0,*) 'setvt: xv = ',xv(mgs,lhl),xdn(mgs,lhl),cx(mgs,lhl),xdia(mgs,lhl,3)
      ELSE
      xdia(mgs,lhl,1) = &
       (qx(mgs,lhl)*rho0(mgs)/(pi*xdn(mgs,lhl)*cno(lhl)))**(0.25)
      cx(mgs,lhl) = cno(lhl)*xdia(mgs,lhl,1)
      xv(mgs,lhl) = Max(xvmn(lhl), rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl)) )
      xdia(mgs,lhl,3) = (xv(mgs,lhl)*6./pi)**(1./3.)
      end if
      else
      xdia(mgs,lhl,1) = 1.e-9
      xdia(mgs,lhl,3) = 1.e-9
      end if
      xdia(mgs,lhl,2) = xdia(mgs,lhl,1)**2
      end do

      ENDIF
     

 
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set terminal velocities'
!
!
! ################################################################
!
!  RAIN
!
      IF ( ildo == 0 .or. ildo == lr ) THEN
      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      IF ( ipconc .lt. 3 ) THEN
        vtxbar(mgs,lr,1) = rainfallfac*(ar*gf4br/6.0)*(xdia(mgs,lr,1)**br)*rhovt(mgs)
!        write(91,*) 'vtxbar: ',vtxbar(mgs,lr,1),mgs,gf4br,xdia(mgs,lr,1),rhovt(mgs)
      ELSE

        IF ( imurain == 1 ) THEN ! DSD of Diameter

        ! using functional form of  arx*(1 - Exp(-frx*diameter) ), with arx =       arx = 10.
        !  and frx = 516.575 ! raind fit parameters for arx*(1 - Exp(-fx*d)), where d is rain diameter in meters.
        ! Similar form as in Atlas et al. (1973), who had 9.65 - 10.3*Exp[-600 * d]


          alp = alpha(mgs,lr)

          vtxbar(mgs,lr,1) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 4.0) ) ! mass weighted

          IF ( infdo .ge. 1 .and. rssflg == 1 ) THEN
            vtxbar(mgs,lr,2) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 1.0) ) ! number weighted
          ELSE
            vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)
          ENDIF

          IF ( infdo .ge. 2 .and. rssflg == 1 ) THEN
            vtxbar(mgs,lr,3) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 7.0) ) ! z-weighted
          ELSE
            vtxbar(mgs,lr,3) = vtxbar(mgs,lr,1)
          ENDIF

!          write(91,*) 'setvt: alp,vn,vm,vz = ',alp,vtxbar(mgs,lr,2), vtxbar(mgs,lr,1), vtxbar(mgs,lr,3),alpha(mgs,lr)

        ELSEIF ( imurain == 3 ) THEN ! DSD of Volume

        IF ( lzr < 1 ) THEN ! not 3-moment rain
        rwdia = Min( xdia(mgs,lr,1), 8.0e-3 )

         vtxbar(mgs,lr,1) = rhovt(mgs)*6.0*pii*( 0.04771 + 3788.0*rwdia -  &
             1.105e6*rwdia**2 + 1.412e8*rwdia**3 - 6.527e9*rwdia**4)

        IF ( infdo .ge. 1 ) THEN
         vtxbar(mgs,lr,2) = (0.09112 + 2714.0*rwdia - 4.872e5*rwdia**2 +  &
                 4.495e7*rwdia**3 - 1.626e9*rwdia**4)*rhovt(mgs)
        ENDIF

        IF ( infdo .ge. 2 ) THEN ! Z-weighted fall speed
        vtxbar(mgs,lr,3)  = rhovt(mgs)*(  &
            0.0911229 +                  &
       9246.494*(rwdia) -               &
       3.2839926e6*(rwdia**2) +          &
       4.944093e8*(rwdia**3) -          &
       2.631718e10*(rwdia**4) )
       ENDIF

        ELSE ! 3-moment rain, gamma-volume

        vr = xv(mgs,lr)
        rnux = alpha(mgs,lr)

        IF ( infdo .ge. 1 .and. rssflg == 1) THEN ! number-weighted; DTD: added size-sorting flag
        vtxbar(mgs,lr,2) = rhovt(mgs)*                             &
          (((1. + rnux)/vr)**(-1.333333)*                         &
         (0.0911229*((1. + rnux)/vr)**1.333333*Gamma_sp(1. + rnux) + &
           (5430.3131*(1. + rnux)*Gamma_sp(4./3. + rnux))/           &
            vr - 1.0732802e6*((1. + rnux)/vr)**0.6666667*         &
            Gamma_sp(1.666667 + rnux) +                              &
           8.584110982429507e7*((1. + rnux)/vr)**(1./3.)*         &
            Gamma_sp(2. + rnux) -                                    &
           2.3303765697228556e9*Gamma_sp(7./3. + rnux)))/            &
       Gamma_sp(1. + rnux)
        ENDIF

!  mass-weighted
       vtxbar(mgs,lr,1)  = rhovt(mgs)*                                                 &
        (0.0911229*(1 + rnux)**1.3333333333333333*Gamma_sp(2. + rnux) +                  &
         5430.313059683277*(1 + rnux)*vr**0.3333333333333333*                         &
          Gamma_sp(2.333333333333333 + rnux) -                                           &
         1.0732802065650471e6*(1 + rnux)**0.6666666666666666*vr**0.6666666666666666*  &
          Gamma_sp(2.6666666666666667 + rnux) +                                          &
         8.584110982429507e7*(1 + rnux)**0.3333333333333333*vr*Gamma_sp(3 + rnux) -      &
         2.3303765697228556e9*vr**1.3333333333333333*                                 &
          Gamma_sp(3.333333333333333 + rnux))/                                           &
       ((1 + rnux)**2.333333333333333*Gamma_sp(1 + rnux))

        IF(infdo .ge. 1 .and. rssflg == 0) THEN ! No size-sorting, set N-weighted fall speed to mass-weighted
          vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)
        ENDIF

        IF ( infdo .ge. 2 .and. rssflg == 1) THEN ! Z-weighted fall speed
        vtxbar(mgs,lr,3)  =   rhovt(mgs)*                                          &
       ((1. + rnux)*(0.0911229*(1 + rnux)**1.3333333333333333*Gamma_sp(3. + rnux) +  &
           5430.313059683277*(1 + rnux)*vr**0.3333333333333333*                   &
            Gamma_sp(3.3333333333333335 + rnux) -                                    &
           1.0732802065650471e6*(1 + rnux)**0.6666666666666666*                   &
            vr**0.6666666666666666*Gamma_sp(3.6666666666666665 + rnux) +             &
           8.5841109824295e7*(1 + rnux)**0.3333333333333333*vr*Gamma_sp(4. + rnux) - &
           2.3303765697228556e9*vr**1.3333333333333333*                           &
            Gamma_sp(4.333333333333333 + rnux)))/                                    &
       ((1 + rnux)**3.3333333333333335*(2 + rnux)*Gamma_sp(1 + rnux))

!         write(0,*) 'setvt: mgs,lzr,infdo = ',mgs,lzr,infdo
!         write(0,*) 'vt1,2,3 = ',vtxbar(mgs,lr,1),vtxbar(mgs,lr,2),vtxbar(mgs,lr,3)
        ELSEIF (infdo .ge. 2) THEN ! No size-sorting, set Z-weighted fall speed to mass-weighted
          vtxbar(mgs,lr,3) = vtxbar(mgs,lr,1)
        ENDIF


        ENDIF
       ENDIF ! imurain

!        IF ( rwrad*mwfac .gt. 6.0e-4  ) THEN
!          vtxbar(mgs,lr,1) = 20.1*Sqrt(100.*rwrad*mwfac)*rhovt(mgs)
!        ELSE
!          vtxbar(mgs,lr,1) = 80.0e2*rwrad*rhovt(mgs)*mwfac
!        ENDIF
!        IF ( rwrad .gt. 6.0e-4  ) THEN
!          vtxbar(mgs,lr,2) = 20.1*Sqrt(100.*rwrad)*rhovt(mgs)
!        ELSE
!          vtxbar(mgs,lr,2) = 80.0e2*rwrad*rhovt(mgs)
!        ENDIF
      ENDIF ! ipconc
      else  ! qr < qrmin
      vtxbar(mgs,lr,1) = 0.0
      vtxbar(mgs,lr,2) = 0.0
      end if
      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set rain vt'

      ENDIF
!
! ################################################################
!
!  SNOW !Zrnic et al. (1993)
!
      IF ( ls .gt. 1 .and. ( ildo == 0 .or. ildo == ls ) ) THEN
      do mgs = 1,ngscnt
      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
        IF ( ipconc .ge. 4 ) THEN
         if ( mixedphase .and. qsvtmod ) then
         else
          IF ( isnowfall == 1 ) THEN
           ! original (Zrnic et al. 1993)
           vtxbar(mgs,ls,1) = 5.72462*rhovt(mgs)*(xv(mgs,ls))**(1./12.)
          ELSEIF ( isnowfall == 2 ) THEN
          ! Ferrier:
            IF ( isnowdens == 1 ) THEN
              vtxbar(mgs,ls,1) = 11.9495*rhovt(mgs)*(xv(mgs,ls))**(0.14)
            ELSE
              vtxbar(mgs,ls,1) = 11.9495*rhovt(mgs)*(xv(mgs,ls)*xdn(mgs,ls)/100.)**(0.14)
            ENDIF
          ENDIF

          IF(sssflg == 1) THEN
            IF ( isnowfall == 1 ) THEN
              vtxbar(mgs,ls,2) = 4.04091*rhovt(mgs)*(xv(mgs,ls))**(1./12.)
            ELSEIF ( isnowfall == 2 ) THEN
            ! Ferrier:
              IF ( isnowdens == 1 ) THEN
                vtxbar(mgs,ls,2) = 7.02909*rhovt(mgs)*(xv(mgs,ls))**(0.14)  ! bug fix 11/15/2015: was rewriting to mass fall speed vtxbar(mgs,ls,1)
              ELSE
                vtxbar(mgs,ls,2) = 7.02909*rhovt(mgs)*(xv(mgs,ls)*xdn(mgs,ls)/100.)**(0.14)  ! bug fix 11/15/2015: was rewriting to mass fall speed vtxbar(mgs,ls,1)
              ENDIF
            ENDIF
          ELSE
            vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
          ENDIF
          vtxbar(mgs,ls,3) = vtxbar(mgs,ls,1)
         endif
        ELSE
         vtxbar(mgs,ls,1) = (cs*gf4ds/6.0)*(xdia(mgs,ls,1)**ds)*rhovt(mgs)
         vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
        ENDIF
      else
      vtxbar(mgs,ls,1) = 0.0
      end if
      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set snow vt'

      ENDIF ! ls .gt. 1
!
!
! ################################################################
!
!  GRAUPEL !Wisner et al. (1972)
!
      IF ( lh .gt. 1 .and. ( ildo == 0 .or. ildo == lh ) ) THEN

      do mgs = 1,ngscnt
      vtxbar(mgs,lh,1) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
       IF ( icdx .eq. 1 ) THEN
         cd = cdx(lh)
       ELSEIF ( icdx .eq. 2 ) THEN
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(xdnmx(lh) - xdn(mgs,lh))/(xdnmx(lh)-xdnmn(lh)) ) )
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lh))/(900. - 300.) ) )
         cd = Max(0.45, Min(1.0, 0.45 + 0.35*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
!         cd = Max(0.55, Min(1.0, 0.55 + 0.25*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
       ELSEIF ( icdx .eq. 3 ) THEN
!         cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 300.) ) )
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
       ELSEIF ( icdx .eq. 4 ) THEN
         cd = Max(cdhmin, Min(cdhmax, cdhmin + (cdhmax-cdhmin)* &
             (cdhdnmax - Max( cdhdnmin, Min( cdhdnmax, xdn(mgs,lh) ) ) )/(cdhdnmax - cdhdnmin) ) )
       ELSEIF ( icdx .eq. 5 ) THEN
         cd = cdx(lh)*(xdn(mgs,lh)/rho_qh)**(2./3.)
       ELSEIF ( icdx .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
         indxr = Int( (xdn(mgs,lh)-50.)/100. ) + 1
         indxr = Min( ngdnmm, Max(1,indxr) )


         delrho = Max( 0.0, 0.01*(xdn(mgs,lh) - mmgraupvt(indxr,1)) )
         IF ( indxr < ngdnmm ) THEN

          axh(mgs) = mmgraupvt(indxr,2) + delrho*(mmgraupvt(indxr+1,2) - mmgraupvt(indxr,2) )
          bxh(mgs) = mmgraupvt(indxr,3) + delrho*(mmgraupvt(indxr+1,3) - mmgraupvt(indxr,3) )


         ELSE
           axh(mgs) = mmgraupvt(indxr,2)
          bxh(mgs) = mmgraupvt(indxr,3)
         ENDIF

         aax = axh(mgs)
         bbx = bxh(mgs)

       ENDIF

      IF ( alpha(mgs,lh) .eq. 0.0 .and. icdx > 0 .and. icdx /= 6 ) THEN
      vtxbar(mgs,lh,1) = (gf4p5/6.0)*  &
       Sqrt( (xdn(mgs,lh)*xdia(mgs,lh,1)*4.0*gr) /  &
         (3.0*cd*rho0(mgs)) )
      ELSE
        IF ( icdx /= 6 ) bbx = bx(lh)
        tmp = 4. + alpha(mgs,lh) + bbx
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!        aax = Max( 1.0, Min(2.0, (xdn(mgs,lh)/400.) ) )
!        vtxbar(mgs,lh,1) =  rhovt(mgs)*aax*ax(lh)*(xdia(mgs,lh,1)**bx(lh)*x)/y

        IF ( icdx > 0 .and. icdx /= 6) THEN
          aax = Sqrt(4.0*xdn(mgs,lh)*gr/(3.0*cd*rho00))
          vtxbar(mgs,lh,1) =  rhovt(mgs)*aax* Sqrt(xdia(mgs,lh,1)) * x/y
          axh(mgs) = aax
          bxh(mgs) = bbx
        ELSEIF (icdx == 6 ) THEN
          vtxbar(mgs,lh,1) =  rhovt(mgs)*aax* xdia(mgs,lh,1)**bbx * x/y
        ELSE
          axh(mgs) = ax(lh)
          bxh(mgs) = bx(lh)
          vtxbar(mgs,lh,1) =  rhovt(mgs)*ax(lh)*(xdia(mgs,lh,1)**bx(lh)*x)/y
        ENDIF

      ENDIF

      IF ( lwsm6 .and. ipconc == 0 ) THEN
         vtxbar(mgs,lh,1) = (330.*gf4br/6.0)*(xdia(mgs,lh,1)**br)*rhovt(mgs)
      ENDIF

      end if
      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt'

      ENDIF ! lh .gt. 1
!
!
!
! ################################################################
!
!  HAIL
!
      IF ( lhl .gt. 1 .and. ( ildo == 0 .or. ildo == lhl ) ) THEN

      do mgs = 1,ngscnt
      vtxbar(mgs,lhl,1) = 0.0
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then

       IF ( icdxhl .eq. 1 ) THEN
         cd = cdx(lhl)
       ELSEIF ( icdxhl .eq. 3 ) THEN
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
       ELSEIF ( icdxhl .eq. 4 ) THEN
         cd = Max(cdhlmin, Min(cdhlmax, cdhlmin + (cdhlmax-cdhlmin)*  &
            (cdhldnmax - Max( cdhldnmin, Min( cdhldnmax, xdn(mgs,lhl) ) ) )/(cdhldnmax - cdhldnmin) ) )
       ELSEIF ( icdxhl .eq. 5 ) THEN
         cd = cdx(lh)*(xdn(mgs,lhl)/rho_qh)**(2./3.)
       ELSEIF ( icdxhl .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
         indxr = Int( (xdn(mgs,lhl)-50.)/100. ) + 1
         indxr = Min( ngdnmm, Max(1,indxr) )


         delrho = Max( 0.0, 0.01*(xdn(mgs,lhl) - mmgraupvt(indxr,1)) )
         IF ( indxr < ngdnmm ) THEN

          axhl(mgs) = mmgraupvt(indxr,2) + delrho*(mmgraupvt(indxr+1,2) - mmgraupvt(indxr,2) )
          bxhl(mgs) = mmgraupvt(indxr,3) + delrho*(mmgraupvt(indxr+1,3) - mmgraupvt(indxr,3) )


         ELSE
          axhl(mgs) = mmgraupvt(indxr,2)
          bxhl(mgs) = mmgraupvt(indxr,3)
         ENDIF

         aax = axhl(mgs)
         bbx = bxhl(mgs)

       ELSE
         cd = Max(0.45, Min(0.6, 0.45 + 0.15*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 500.) ) )
       ENDIF

      IF ( alpha(mgs,lhl) .eq. 0.0 .and. icdxhl > 0 .and. icdxhl /= 6) THEN
      vtxbar(mgs,lhl,1) = (gf4p5/6.0)* &
       Sqrt( (xdn(mgs,lhl)*xdia(mgs,lhl,1)*4.0*gr) / &
         (3.0*cd*rho0(mgs)) )
      ELSE
        IF ( icdxhl /= 6 ) bbx = bx(lhl)
        tmp = 4. + alpha(mgs,lhl) + bbx
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
        tmp = 4. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        IF ( icdxhl > 0 .and. icdxhl /= 6) THEN
          aax = Sqrt(4.0*xdn(mgs,lhl)*gr/(3.0*cd*rho00))
          vtxbar(mgs,lhl,1) =  rhovt(mgs)*aax* Sqrt(xdia(mgs,lhl,1)) * x/y
          axhl(mgs) = aax
          bxhl(mgs) = bbx
        ELSEIF ( icdxhl == 6 ) THEN
          vtxbar(mgs,lhl,1) =  rhovt(mgs)*aax* (xdia(mgs,lhl,1))**bbx * x/y
        ELSE
          axhl(mgs) = ax(lhl)
          bxhl(mgs) = bx(lhl)
         vtxbar(mgs,lhl,1) =  rhovt(mgs)*(ax(lhl)*xdia(mgs,lhl,1)**bx(lhl)*x)/y
        ENDIF

      ENDIF


      end if
      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt'

      ENDIF ! lhl .gt. 1


      IF ( infdo .ge. 1 ) THEN

!      DO il = lc,lhab
!      IF ( il .ne. lr ) THEN
        DO mgs = 1,ngscnt
          vtxbar(mgs,lc,2) = vtxbar(mgs,lc,1)
        IF ( li .gt. 1 ) THEN
!          vtxbar(mgs,li,2) = rhovt(mgs)*49420.*1.25447*xdia(mgs,li,1)**(1.415) ! n-wgt (Ferrier 94)
!          vtxbar(mgs,li,2) = vtxbar(mgs,li,1)

! test print stuff...
!          IF ( xdia(mgs,li,1) .gt. 200.e-6 ) THEN
!            tmp = (xv(mgs,li)*cwc0)**(1./3.)
!            x = rhovt(mgs)*49420.*40.0005/5.40662*tmp**(1.415)
!            y = rhovt(mgs)*49420.*1.25447*tmp**(1.415)
!            write(6,*) 'Ice fall: ',vtxbar(mgs,li,1),x,y,tmp,xdia(mgs,li,1)
!          ENDIF
        ENDIF
!          vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
        ENDDO

        IF ( lg .gt. lr ) THEN

        DO il = lg,lhab
         IF ( ildo == 0 .or. ildo == il ) THEN

            DO mgs = 1,ngscnt
             IF ( qx(mgs,il) .gt. qxmin(il) ) THEN
              IF ( (il .eq. lh .and. hssflg == 1) .or. ( lhl .gt. 1 .and. il .eq. lhl .and. hlssflg == 1) ) THEN ! DTD: added flag for size-sorting
              ! DTD: allow for setting of number-weighted and z-weighted fall speeds to the mass-weighted value,
              ! effectively turning off size-sorting

              IF ( il .eq. lh ) THEN ! {

               IF ( icdx .eq. 1 ) THEN
                 cd = cdx(lh)
               ELSEIF ( icdx .eq. 2 ) THEN
!                 cd = Max(0.6, Min(1.0, 0.6 + 0.4*(xdnmx(lh) - xdn(mgs,lh))/(xdnmx(lh)-xdnmn(lh)) ) )
!                 cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lh))/(900. - 300.) ) )
                 cd = Max(0.45, Min(1.0, 0.45 + 0.35*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
!                 cd = Max(0.55, Min(1.0, 0.55 + 0.25*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
               ELSEIF ( icdx .eq. 3 ) THEN
!                 cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 170.0, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
                 cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
               ELSEIF ( icdx .eq. 4 ) THEN
                 cd = Max(cdhmin, Min(cdhmax, cdhmin + (cdhmax-cdhmin)* &
                 (cdhdnmax - Max( cdhdnmin, Min( cdhdnmax, xdn(mgs,lh) ) ) )/(cdhdnmax - cdhdnmin) ) )
               ELSEIF ( icdx .eq. 5 ) THEN
                 cd = cdx(lh)*(xdn(mgs,lh)/rho_qh)**(2./3.)
               ELSEIF ( icdx .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
                  aax = axh(mgs)
                  bbx = bxh(mgs)
               ENDIF

              ELSEIF ( lhl .gt. 1 .and. il .eq. lhl ) THEN

               IF ( icdxhl .eq. 1 ) THEN
                 cd = cdx(lhl)
               ELSEIF ( icdxhl .eq. 3 ) THEN
!               cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 300.) ) )
                cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
               ELSEIF ( icdxhl .eq. 4 ) THEN
                cd = Max(cdhlmin, Min(cdhlmax, cdhlmin + (cdhlmax-cdhlmin)*  &
                    (cdhldnmax - Max( cdhldnmin, Min( cdhldnmax, xdn(mgs,lhl) ) ) )/(cdhldnmax - cdhldnmin) ) )
               ELSEIF ( icdxhl == 5 ) THEN
!                cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lhl))/(900. - 300.) ) )
!                cd = Max(0.5, Min(0.8, 0.5 + 0.3*(xdnmx(lhl) - xdn(mgs,lhl))/(xdnmx(lhl)-xdnmn(lhl)) ) )
                 cd = Max(0.45, Min(0.6, 0.45 + 0.15*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 500.) ) )
               ELSEIF ( icdxhl .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
                  aax = axhl(mgs)
                  bbx = bxhl(mgs)
               ENDIF

              ENDIF ! }

               IF ( alpha(mgs,il) .eq. 0. .and. infdo .lt. 2 .and.   &
               ( ( il==lh .and. icdx > 0 .and. icdx /= 6) .or. ( il==lhl .and. icdxhl > 0 .and. icdxhl /= 6 ) ) ) THEN ! {
                 vtxbar(mgs,il,2) =   &
                   Sqrt( (xdn(mgs,il)*xdia(mgs,il,1)*pi*gr) / &
                     (3.0*cd*rho0(mgs)) )

               ELSE
               IF ( il == lh  .and. icdx   /= 6 ) bbx = bx(il)
               IF ( il == lhl .and. icdxhl /= 6 ) bbx = bx(il)
               tmp = 1. + alpha(mgs,il) + bbx
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
               tmp = 1. + alpha(mgs,il)
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

                 IF ( il .eq. lh  .or. il .eq. lhl) THEN ! {
                   IF ( ( il==lh .and. icdx > 0 ) ) THEN
                     IF ( icdx /= 6 ) THEN
                      aax = Sqrt(4.0*xdn(mgs,il)*gr/(3.0*cd*rho00))
                      vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**0.5 * x/y
                     ELSE !  (icdx == 6 ) THEN
                       vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bbx * x/y
                     ENDIF
!                   ELSE
!                     aax = ax(il)
!                     vtxbar(mgs,il,2) =  rhovt(mgs)*ax(il)*(xdia(mgs,il,1)**bx(il)*x)/y
!                   ENDIF

                   ELSEIF ( ( il==lhl .and. icdxhl > 0 ) ) THEN
                     IF ( icdxhl /= 6 ) THEN
                       aax = Sqrt(4.0*xdn(mgs,il)*gr/(3.0*cd*rho00))
                       vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**0.5 * x/y
                     ELSE ! ( icdxhl == 6 )
                       vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bbx * x/y
                     ENDIF
                   ELSE
                     aax = ax(il)
                     vtxbar(mgs,il,2) =  rhovt(mgs)*ax(il)*(xdia(mgs,il,1)**bx(il)*x)/y
                   ENDIF

!                  vtxbar(mgs,il,2) =  &
!     &               rhovt(mgs)*(xdn(mgs,il)/400.)*(75.715*xdia(mgs,il,1)**0.6* &
!     &               x)/y
!                  vtxbar(mgs,il,2) =  &
!     &               rhovt(mgs)*(xdn(mgs,il)/400.)*(ax(il)*xdia(mgs,il,1)**bx(il)* &
!     &               x)/y
                  IF ( infdo .ge. 2 ) THEN ! Z-weighted
                   vtxbar(mgs,il,3) = rhovt(mgs)*                 &
                     (aax*(1.0/xdia(mgs,il,1) )**(- bbx)*  &
                      Gamma_sp(7.0 + alpha(mgs,il) + bbx))/Gamma_sp(7. + alpha(mgs,il))
!     &                (aax*(1.0/xdia(mgs,il,1) )**(- bx(il))*  &
!     &                 Gamma_sp(7.0 + alpha(mgs,il) + bx(il)))/Gamma_sp(7. + alpha(mgs,il))
                  ENDIF

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt3'

                 ELSE ! hail
                  vtxbar(mgs,il,2) =  &
                    rhovt(mgs)*(ax(il)*xdia(mgs,il,1)**bx(il)* &
                    x)/y

                 IF ( infdo .ge. 2 ) THEN ! Z-weighted
                  vtxbar(mgs,il,3) = rhovt(mgs)*                 &
                   (aax*(1.0/xdia(mgs,il,1) )**(- bbx)*  &
                    Gamma_sp(7.0 + alpha(mgs,il) + bbx))/Gamma_sp(7. + alpha(mgs,il))
!     &              (ax(il)*(1.0/xdia(mgs,il,1) )**(- bx(il))*  &
!     &               Gamma_sp(7.0 + alpha(mgs,il) + bx(il)))/Gamma_sp(7. + alpha(mgs,il))
                  ENDIF
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt4'

                 ENDIF ! }
!     &             Gamma_sp(1.0 + dnu(il) + 0.6)/Gamma_sp(1. + dnu(il))
               ENDIF ! }

!              IF ( infdo .ge. 2 ) THEN ! Z-weighted
!               vtxbar(mgs,il,3) = rhovt(mgs)*                 &
!     &            (ax(il)*(1.0/xdia(mgs,il,1) )**(- bx(il))*  &
!     &             Gamma_sp(7.0 + alpha(mgs,il) + bx(il)))/Gamma_sp(7. + alpha(mgs,il))
!              ENDIF

!               IF ( lhl .gt. 1 .and. il .eq. lhl ) THEN
!                write(0,*) 'setvt: ',qx(mgs,il),xdia(mgs,il,1),xdia(mgs,il,3),dnu(il),ax(il),bx(il)
!               ENDIF
             ELSEIF ( (il .eq. lh .and. hssflg == 0) .or. ( lhl .gt. 1 .and. il .eq. lhl .and. hlssflg == 0) ) THEN ! no size-sorting for graupel or hail
              vtxbar(mgs,il,2) = vtxbar(mgs,il,1)
              vtxbar(mgs,il,3) = vtxbar(mgs,il,1)
             ELSE ! not lh or lhl
              vtxbar(mgs,il,2) = &
                 Sqrt( (xdn(mgs,il)*xdia(mgs,il,1)*pi*gr) /  &
                   (3.0*cdx(il)*rho0(mgs)) )
              vtxbar(mgs,il,3) = vtxbar(mgs,il,1)

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt5'


              ENDIF
             ELSE ! qx < qxmin
              vtxbar(mgs,il,2) = 0.0

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt6'

             ENDIF
           ENDDO ! mgs

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt7'

        ENDIF
        ENDDO ! il

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt8'

        ENDIF ! lg .gt. 1 

!      ENDIF
!      ENDDO

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt9'

!       DO mgs = 1,ngscnt
!        IF ( qx(mgs,lr) > qxmin(lr) ) THEN
!         write(0,*) 'setvt2: mgs,lzr,infdo = ',mgs,lzr,infdo
!         write(0,*) 'vt1,2,3 = ',vtxbar(mgs,lr,1),vtxbar(mgs,lr,2),vtxbar(mgs,lr,3)
!        ENDIF
!       ENDDO

      ENDIF ! infdo .ge. 1 


!############ SETVTZ ############################

      RETURN
      END SUBROUTINE setvtz

!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 
!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 

      REAL FUNCTION GAMMA_SP(xx)

      implicit none
      real xx
      integer j

! Double precision ser,stp,tmp,x,y,cof(6)

      real*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d+0,  &
     &            -86.50532032941677d0,   &
     &             24.01409824083091d0,   &
     &             -1.231739572450155d0,  &
     &              0.1208650973866179d-2,&
     &             -0.5395239384953d-5,   &
     &              2.5066282746310005d0/

      IF ( xx <= 0.0 ) THEN
        write(0,*) 'Argument to gamma must be > 0!! xx = ',xx
        STOP
      ENDIF

      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = (x + 0.5d0)*Log(tmp) - tmp
      ser = 1.000000000190015d0
      DO j=1,6
        y = y + 1.0d0
        ser = ser + cof(j)/y
      END DO
      gamma_sp = Exp(tmp + log(stp*ser/x))

      RETURN
      END FUNCTION GAMMA_SP

!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 
!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 
 
 end module module_nssl_2mmt

!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 
!NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL NSSL 

