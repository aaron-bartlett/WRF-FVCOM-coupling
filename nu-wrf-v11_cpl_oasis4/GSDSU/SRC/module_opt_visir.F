#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_opt_visir
 use module_simulator  
 use module_surface
 use module_lidar
 use module_gmp_4ice
 use module_mpi
 use module_icescat_visir
 use module_icescat_yang_visir
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
!
! Comments:  
!   This module deal with visible-IR particle single scattering and gaseous absorption/scattering for lidar and
!   visir simulators. Single-scattering LUTs are pre-computed. Cloud and aerosol PSDs assumptions and density 
!   are maded consistent to broadband modules. Note that visir module define double-precision floating point. 
!   , while Goddard SDSU module define single-precision floating point. 
!
! 
! History:
!  07/2010  Toshi Matsui@NASA GSFC : Clearn up / re-organization of subroutine
!  04/2009  Toshi Matsui@NASA GSFC : Add lidar component. 
!  03/2009  Toshi Matsui@NASA GSFC : RSTAR6bF90 was modulated to Goddard SDSU framework
! ~03/2009  Miho Sekiguchi@        : Upgrade RSTAR from f77 to f90. 
!           
!---------------------------------------------------------------------------------------------------

  private  !encapselate all parameters 
           !(You need add public statement for parameter to interact outside of this module)

!
! public subourtines (accesible from other routines/modules)
!
  public :: opt_visir, &
            makelut_visir, &
            makelut_visir_grnd,&
            read_rfi,&
            mie_rams_visir,&
            mie_gen_visir,&
            mie_gocart,&
            mie_ucd,&
            makelut_visir_aerosol


  integer,parameter,public :: &  ! Array Index
      kna0 = 1   ,&  ! max number of solar incidences (always 1 for SDSU)
      kna1u= 1   ,&  ! max number of emergent nadir angles (always 1 for SDSU)
      knfi = 1   ,&  ! max number of azimuthal angles  (always 1 for SDSU)
      knw0 = 50  ,&  ! max number of sub-wavelength grids
      knw  = max_chan ,& ! max number of wavelength grid (=mxwavel) 
      kntau= 501 ,&!mxlyr+1 ,& ! max number of optical thickness
      knln = 500 ,&!mxlyr   ,& ! max number of layers 
      kndm = N_stream,& ! max number of matrix size (=quadrature number)
      knl  = 100 ,&  ! max number of sounding 
      knm0 = 28  ,&  ! max number of gas species
      kpoly = 4  ,&  ! max number of input cloud/aerosol species (not used for SDSU)
      kaw  = 17  ,&  ! max number of aerosol growth
      kwlv = 210 ,&  ! max number of refractive index
      knv  = 9   ,&  ! max number of fundamental particle species
      kptc = 20  ,&  ! max number of cloud/aerosol species model
      klgn1= N_stream*2+1 ,& ! max number of phase function 
      kp   = 26     ,& ! number of pressure grids
      kt   = 3      ,& ! number of temperature grids
      kmol = 7      ,& ! number of major gases
      kmol1= kmol+1 ,& ! number of major gases + self
      nch_k=2       ,& ! for k-distribution table
      dintvl=1.d3   ,& ! for k-distirubition
      kda  = 2      ,& ! not used
      knang= 80     ,& ! max number of angular point in phase function  
      kplk1= 4      ,& ! 2 SDSU  !! 
      kww  = 200    ,& ! for k-distribution table
      kch  = 2         ! for k-distribution table

! for rfidxb
  integer,parameter,public:: & ! for mie routine 
      kintvl = 75    ,& !
      kpol = 1       ,& ! 
      knm  = knm0    ,& ! 
      kgs  = 53      ,& ! 
      kna1 = 4       ,& !
      knln1 = knln + 1 ,& ! 
      knw01 = knw0 + 1 ,& ! 
      kndm2 = 2 * kndm ,& ! 
      knlnm = knln1    ,& !
      knlnm1 = knlnm + 1 ,& !
      krowij = kndm * (kndm + 1) / 2  ,& !
      knltm=knln1       ,& !
      knsb=knln1        ,& !
      knlt=knln1

  real(8),parameter,public :: pi= 3.14159265358979d0 ,rad=1.745329251994330d-002
  real(8),parameter,public :: pstd=1013.25d0, tstd=273.15d0

! rfi stuff from refractive index table
  integer,save:: nv              ! number of fundamental species (1-8)
  integer,save:: nwlv            ! number of wavelenggths for refractive index(rfi) tables
  real(8),save:: wlv(kwlv)       ! wavelengths [micron] for rfi tables
  real(8),save:: rfi(kwlv,knv,2) ! refractive index of fundamental materials (mr,mi) = mr - i mi
                                 ! with log-regular wavelength interval

! for abgask
 real(8), allocatable :: taug(:,:,:)    !! optical thickness of gas absorption: dimension(1:knln,1:kch,1:kww)

 real(8),parameter :: q_min_conden = q_min_condensate*1.0d-6  ! minimum detectable condensate amount from CRM [g/cm3] 
                                                              ! consistent criteria to other modules
 real(8),parameter :: q_min_aero   = q_min_aerosol*1.0d-6     ! minimum detectable condensate amount from Aerosol Model [g/cm3]


!
! LUT for surface reflection/source matrix
!
  integer, parameter,public :: mxpts_u10   = 25   ! LUT dimension for u10 
  integer, parameter,public :: mxpts_cosz  = 100  ! LUT dimension for u10 

  real(8),save,allocatable,public :: orf_lut(:,:,:,:,:) !Look-up table of surface reflection matrix
  real(8),save,allocatable,public :: srr_lut(:,:,:,:,:) !Look-up table of surface source matrix

 integer,public :: iw_global ! channel order 


!
! LUT
!
 real(sdsu_fpd),allocatable :: &
     cextp_lut_array(:,:,:),  & ! extinction cross section array [cm2/cm3] or [cm-1]
     cabsp_lut_array(:,:,:),  & ! absorption cross section array [cm2/cm3] or [cm-1]
     php_lut_array(:,:,:,:)     ! volume scattering phase function array [cm-1 str-1]

 real(sdsu_fpd),allocatable :: &
     cextp_aero_lut_array(:,:,:),  & ! extinction cross section array [cm2/cm3] or [cm-1]
     cabsp_aero_lut_array(:,:,:),  & ! absorption cross section array [cm2/cm3] or [cm-1]
     php_aero_lut_array(:,:,:,:)     ! volume scattering phase function array [cm-1 str-1]

!
! Misc
!

 character :: err*64  !error messages

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_visir(sensor,umu_grid )
 implicit none

!----------------------------------------------------------------------------------------------------
!
! Comments:
!  This subroutine perpares the visible IR (0.17 ~ 1000 microns) optical properties
!  for a given atmosphere and condensates amounts. It assumes a plane parallel 
!  atmosphere divided into several homogeneous sublayers with underlying ground or ocean surface.
!  Microphysical assumptions such as drop-size distribution and ice-particle
!  characteristics are built in the Mie calculation subroutines (mie_sbm or mie_gen_visir).
!  for details on prescribed microphysical assumptions, which users
!  might want to modify. Dielectric properties of fundamental particles are
!  prepared from Look-up table (). Gaseous absorptions are accounted through the non-linear 
!  fitting k-distribution method of Sekiguchi and Nakajima (2008).
!  It contains a k-distribution table with HITRAN2004,
!  but this version is now only treated major seven gases.
!  (H2O, CO2, O3, N2O, CO, CH4, O2).
!  N2 is needed for sattelite analyses, but not included.
!
!  Pre-computed lookup tables containing extinction and scattering cross
!  sections are optionally available to save time to calculate them online.
!  The Mie calculation subroutines mentioned above are shared by the program
!  to construct the lookup table. Any changes made to the Mie routines
!  therefore will be automatically reflected in the code to creat the
!  lookup tables when it is compiled next. 
!
!
! Methods: Calling sequence 
!
! opt_visir (subroutine to compute optical properties at visible and IR spectrum)
!   |-- read_sounding (initialization of atmospheric profile from 'MLATMD')
!   |-- read_rfi (initialize cloud/aerosol optical parameter from 'rhi_table.txt' )
!   |-- wvcal (conversion of units for water vapor)
!   |-- sdp2 (derive absorption coefficient from pt table and fitting)
!   |-- fsoli (wavelength integrated solar incidence)
!   !-- mie_sbm or mie_gen_visir (setup particle parameters)
!   !     |-- rfidxb_sdsu (estimate total refractive index & get parameters for wet particle polydispersions)
!   !           |-- gtph4b (get phase function and cross section)
!   !           |     |-- smllop (optical cross sections for small particles)
!   !           |     |-- nons (empirical correction for nonspherical particle)
!   |-- abgask (calculate absorption coefficients)
!   |     |-- tdok2 ()
!   |-- qgausn (compute weights and abscissae for ordinary gaussian quadrature)
!   |-- simulator_visir(solve the radiative transfer in atmosphere system)
!   |      or
!   |-- simulator_lidar(compute lidar backscattering)
!
! History:
!  03/2009  Toshi Matsui@NASA GSFC : RSTAR6bF90 was modulated to SDSU framework
!  ~11/2008  Miho Sekiguchi        : F90 version of RSTAR 6b (See origional LONG history below)
!
!
! References: 
! Nakajima, T., and M. Tanaka, 1986: Matrix formulations for the transfer
!   of solar radiation in a plane-parallel scattering atmosphere. J.
!   Quant. Spectrosc. Radiat. Transfer, 35, 13-21.
!
! Nakajima, T., and M. Tanaka, 1988: Algorithms for radiative intensity
!   calculations in moderately thick atmospheres using a truncation
!   approximation. J. Quant. Spectrosc. Radiat. Transfer , 40, 51-69.
!
! Sekiguchi, M., and T. Nakaima, 2008: A study of the absorption process and
!   its computational optimization in an atmospheric general circulation model.
!   J. Quant. Spectrosc. Radiat. Transfer, submitted.
!
!
!==================================================================================================
!--- Original history (rstar main)
! 95. 2.20 aerosol growing with relative humidity (I. Lensky)
!     8.11 introduction of taumax
!          modified vlspcv
! 95. 9.15 re-organized by t. nakajima, a. higurashi
! 95.12.30 updated humidity, icn definitions in rstr5
! 96. 2. 7 registered aerd in lbr2
! 96. 3.11 allow non-spherical parameters
!     5. 5 add ncomp, change ints5b arguments
! 98. 1. 7
!   u10 > 0.01 for rstr and pstr.  u10 < 0.01 will have significant
!   error in radiances due to angular integration errors for surface
!   reflection.  recommend that u10 is set to be 0.01 for flat surface
!   case
! 07.11.30 rstar6 with new gas absorption k-distribution (M. Sekiguchi)
! 08. 2.19 added flxd00 (T. Nakajima)
! 08. 2.19 rstr6b without ifrh and trh (T. Nakajima)
! 08. 8.14 add dintvl (M.Sekiguchi)
! 08. 9. 5 
!--- history (rstar's rstar 6b)
! 95. 9.15 created
! 95.12.28 corrected the comments for icn=1, and added icn=3 case.
! 95.12.30 modified the relative humidity profile calculation
! 96. 3.11 added non-spherical parameters
!     5. 5 introduce ncomp and some bug-fix for cnprf1
!     5. 9 cextp() -> cextp2, cabsp() -> cabsp2, voldp()->voldp2, volp()->volp2
!     5.24 change the structure of do 18
! 96.12.27 bug in initialization with do 17
! 96.12.28 change the meaning of rx
! 97. 1.14  20.0->200.0
! 97. 3.17  get scr and sci from subroutine rfidxb.
!           call rtrn21 with scr and sci.
! 98.1.7
!   u10 > 0.01 for rstr and pstr.  u10 < 0.01 will have significant
!   error in radiances due to angular integration errors for surface
!   reflection.  recommend that u10 is set to be 0.01 for flat surface
!   case
! 06.06.25 version up to 6. (miho)
! 07.03.28 fixed some bugs. (miho)
! 07.11.30 dynamical => statical dimension (nln=>knln, nch=>kch) (miho) 
! 08. 2.19 flxd00 for direct solar flux
! 08. 2.19 added notes for thk0
! 08. 2.19 deleted ifrh and trh
! 08. 8.14 rtrn2 --> rtrn21, add dintvl
! 08. 9.19 modified free-style form
!
!
!--- input
! init       i        if 1 initilize the routine
!                     0 skip wavelength-independent part
! na0        i        number of solar incidences.
! th0     r(na0)      solar zenith angles 0-90 degrees
! na1        i        number of emergent nadir angles in sphere.
! th1    r(kna1u)     emergent nadir angles.
!                      0-90 for downward; 90-180 for upward
! nfi        i        number of azimuthal angles.
! fi     r(knfi)      azimuthal angles in degrees.
! wl         r        center wavelengths (micron).
! dw         r        scaling factor to define scaled wavelengths rx.
!                     give 0 for monochromati! calculations.
! nw0        r        number of sub-wavelength grids for calculating response
!                      function.
! rx     r(knw0+1)    dw*rx(i) = boundaries in micron for sub-wavelenth integration
!                     of the response function rf.  -0.5 <= rx <= 0.5
! rf      r(knw0)     integrated response function of radiometer over
!                       wavelength interval of (rx(i)*dw, rx(i+1)*dw)
!                     averaged radiance is given as
!                       l = sum (i=1,n; radiance(i)*rf(i))/sum(i=1,n; rf(i))
!                     give nw0=1, rx(1)=rx(2)=0, rf(1)=1 for monochromatic
!                     calculations.
! galb       r        ground albedo if indg=0
!                     u10 (m/sec)   if indg>0; u10>0.01
! wlcn       r        scaling wavelength for icn=2,3
! npoly      i        number of external mixtures
! ncomp   i(npoly)    number of aerosol types for the external mixture
! mptc    i(kpoly,4)  aerosol types for the external mixture
! cnpt    r(npoly)    total volume [cm3/cm2] (icn=0,1) or optical thickness (icn=2,3)
!                      for each polydispersion.
! nl         i        number of atmospheric layers to define model atmosphere
! alt     r(knl)      height at interfaces of layers (km), top to bottom
! prs     r(knl)      pressure    at interfaces of layers (mb)
! tmp     r(knl)      temperature at interfaces of layers (k)
!! nmol       i        number of gases
! cng    r(knl,knm0)  gas concentration (ppmv)
! nptc       i        number of defined aerosol types
! cnp    r(knl,kptc)  dry volume concentration profile (relative unit)
! ispcvp i(3,kptc)    fundamental materials for 3 comp. internal mixture (1-8)
! rfracp r(3,kptc)    dry volume fraction of the dry mixture
! asphr  r(3,kptc)    non-spherical parameters (x0, g, r)
! rop    r(kptc)      paricle density relative to water
! nv         i        number of fundamental species (1-8)
! nwlv       i        number of wavelengths for refractive index tables
! wlv    r(kwl)       wavelengths (micron)
! rfi   r(kwl,knv,2)  refractive index of fundamental materials (mr, mi)
!                     =mr - i mi
!                     with log-regular wavelength interval
! nln        i        number of atmospheric layers for radiative transfer.
! ipbf    i(nln1)     interface number to define sublayers to construct
!                      the transfer atmosphere in mlatm.
!                     ipbf(1)=50, ipbf(nln+1)=1
!                     if >0 then field calculations
!                        <0 then no field calculations
!--- output
! thk0    r(nln,10)   (l,1): total optical thickness for extinction (tau)
!                     (l,2): particle optical thickness for extinction (taua)
!                     (l,3): optical thickness for rayleigh scattering (taur)
!                     (l,4): single scattering albedo (w)
!                     (l,5): g1 (1st legendre polynomial moment)
!                     (l,6): g2 (2nd legendre polynomial moment) 
!   these values are not for monochromatic but for wavelength-mean values
!   other important parameters can be obtained as follows:
!   optical thickness for gas absorption: taug= tau - taua - taur
!   single scattering albedo of particles: wa = (w*tau-taur)/taua
!   legendre poly. moments of particles: gan= (w*tau*gn - grn*taur)/(wa*taua)
!    because  tau= taua + taur + taug
!             w*tau = wa*taua + taur
!             w*tau*gn = wa*taua*gan + grn*taur
!    where gr is the phase function of rayleigh scattering phase function
!   good approximation is: gr1= 0 and gr2= 0.1
! sol        r        solar incident irradiance (w/m2/micron)
! flxd0  r(kna0,kntau) downward flux (w/m2/micron)
! flxd00 r(kna0,kntau) downward direct flux (horizontal plane flux)
! flxu0  r(kna0,kntau) upward flux (w/m2/micron)
! ai0    r(kna1u,kna0,knfi,kntau) radiance (w/m2/micron/str)
!-----------------------------------------------------------------------------------
 character(len=5),intent(in) :: sensor  !'visir' or 'lidar'
 real(sdsu_fps),intent(in),optional :: umu_grid(:,:) ! cosine of local incident angle [-]

! for main
  character:: fnm*100     !! file name
! for data
  integer:: isol          !! solar incident flag (0:night, 1:day)
  integer:: inda          !! radiance flag (0: flux only, 1: flux & radiance)
  integer:: indg          !! surface flag (0:Lambert, 1: Ocean)
                          ! -1: no ground surface
                          ! 0: lambert surface
                          ! 1: ocean surface initialization
                          ! 2: ocean surface with no initialization !Do not use for SDSU since u10 is grid dependent (toshi)
                          ! when indg>0 and imthd>0 then cabspsingle scattering correction
                          ! for ocean surface reflection
  integer:: imthd         !! flux/intensity calculation modes
                          !-1: nt,  0: dms-method  for intensity/flux
                          !1: ms,  2:tms,  3:ims-method for intensity.
                          !when indg>0 and imthd>0 then single scattering correction
                          !for ocean surface reflection
  integer:: nda           !! number of (stream) nunber of nadir-quadrature angles in the hemisphere.

  integer:: na0           !! number of solar zenith angles
  real(8):: th0(kna0)     !! solar zenith angles (degree)

  integer:: na1           !! number of zenith angles of emerging radiance
  real(8):: th1(kna1u)    !! zenith angles of emerging radiances (degree)

  integer:: nfi           !! number of azimuthal angles
  real(8):: fi(knfi)      !! azimuthal angles (degree)

  integer:: nw0           !! number of sub-wavelength grids
  real(8):: rx(knw01)     !! boundaries for sub-wavelength integration
  real(8):: rf(knw0)      !! integrated response function

  integer:: nw            !! number of center wavelength
  real(8):: wl(knw)       !! center wavelength (micron)
  real(8):: dw(knw)       !! scaling factor to define scaled wavelengths
  real(8):: galb(knw)     !! ground albedo (indg=0)/ U10 [m/sec] (indg=1)
  real(8):: galb_gwind 
  integer:: matm          !! atmospheric sounding index [1~6]
  integer:: nln           !! number of sublayers

  integer:: ipbf(knln1)   !! interface number to define sublayers
  integer:: ifrh          !! humidity flag (0: model, 1: use trh)
  real(8):: trh           !! total relative humiditiy

  real(8):: wlcn          !! scaling wavelength [cm] (icn=2 or 3)

! for INTS5B
  integer:: nl            !! number of atmospheric layer 
  real(8):: alt(knl)      !! altitude [km]
  real(8):: prs(knl)      !! pressure [mb]
  real(8):: tmp(knl)      !! temperature [K]
  real(8):: cng(knl,knm0) !! gas concentration [ppmv]

! for RSTR6B
  integer:: init          !! initialization flag for RSTR6B

! for output
  real(8):: thk0(knln,10) !! calculated optical parameters
  real(8):: sol           !! solar incident irradiance [W/m2/micron]
  real(8):: flxd0(kna0,kntau) !! downward flux [W/m2/micron]
  real(8):: flxu0(kna0,kntau) !! upward flux [W/m2/micron]
  real(8):: flxd00(kna0,kntau) !! downward direct flux [W/m2/micron]
  real(8):: ai0(kna1u,kna0,knfi,kntau) !! radiance [W/m2/micron/str]

! Lidar
  integer :: n_diag, k_rev, k_inst
  real(8) :: thk0_inst(mxlyr_lidar,5) !! calculated optical parameters
  real(8) :: kext_inst(mxlyr_lidar) , kext_clear_inst(mxlyr_lidar)
  real(8) :: sbacktot_inst(mxlyr_lidar), sbacktot_p22_inst(mxlyr_lidar), sbacktot_clear_inst(mxlyr_lidar)
  real(8) :: bcat_inst(mxlyr_lidar), bcat_p22_inst(mxlyr_lidar),  bcat_clear_inst(mxlyr_lidar)
  real(8) :: bcat0_inst(mxlyr_lidar),bcat0_p22_inst(mxlyr_lidar), bcat0_clear_inst(mxlyr_lidar)
  real(8) :: dzm_inst(mxlyr_lidar)   ! [m]

! for work
  integer:: i,l,iw   !! do-loop indicator
  integer:: ierr          !! error indicator
  integer:: lf            !! layer counter

! for wvcal
  integer:: m,m1,ind
  real(8):: p1,ppmv_h2o,gm3,e,gm3s,es

! for rfidxb
  real(8):: scr                  !! relative water reflactive index (real part)
  real(8):: sci                  !! relative water reflactive index (imaginary part)
  integer,save:: nang                 !! number of scaattering angles
  real(8),save:: ang(knang)           !! scattering angles [degree]
  real(8):: cextp_fake  !extinction cross section [cm2/cm3] or [cm-1]
  real(8):: cabsp_fake !absorption cross section [cm2/cm3] or [cm-1]
  real(8),dimension(knang) :: php_fake   !! volume scattering phase function


  type ( particle_gen_r8 ) :: & !particle_gen is defined in module_simulater
     cextp_gen,  & ! !extinction cross section [cm2/cm3] or [cm-1]
     cabsp_gen,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_gen(knang)! volume scattering phase function [-]

  type ( particle_rams_r8 ) :: & !particle_rams is defined in module_simulater
     cextp_rams,  & ! !extinction cross section [cm2/cm3] or [cm-1]
     cabsp_rams,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_rams(knang)! volume scattering phase function [-]

  type ( particle_sbm_r8 ) :: & !particle_sbm is defined in module_simulater
     cextp_sbm,  & ! !extinction cross section [cm2/cm3] or [cm-1]
     cabsp_sbm,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_sbm(knang)! volume scattering phase function [-]

  type ( particle_gocart_r8 ) :: & !particle_gocart is defined in module_simulater
     cextp_gocart,  & ! !extinction cross section [cm2/cm3] or [cm-1]
     cabsp_gocart,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_gocart(knang)! volume scattering phase function [-]

  type ( particle_ucd_r8 ) :: & !particle_ucd is defined in module_simulater
     cextp_ucd,  & ! !extinction cross section [cm2/cm3] or [cm-1]
     cabsp_ucd,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_ucd(knang)! volume scattering phase function [-]

 integer :: nt,ns !temp and rh index
 integer :: ispc  !rams species index
 real(8) :: cfmas, pwmas !alpha_m, beta_m for mass-diameter index


! for sdp2
  real(8),save:: dp(knln)             !! weights in pressure grids
  integer,save:: np(knln)             !! numbers in pressure grids

! for abgask
  integer:: nww                  !! number of integration
  real(8):: wlb                  !! longer wavelength boundary
  real(8),save:: amtpb(kmol1,knln)    !! gas concentration
  real(8),save:: tl(knln)             !! temperature in each layer

! for qgausn
  real(8):: gwt(kch)             !! Gaussian weights 
  real(8):: gmu(kch)             !! Gaussian points

! for rtrn21
  integer,save:: indt                 !! flag of optical depth (1:utau, 2:dpt)
  integer,save:: indp                 !! flag of phase function
  real(8):: am1u(kna1u)          !! cosine of the emergent nadir angles
  real(8):: am0(kna0)            !! cosine of solar zenith angles
  real(8):: tau(knln)            !! total optical thickness of sublayers
  real(8),dimension(knln) :: omg,omga,omgc,omgr    !! single scattering albedo (total,particle,rayleigh)
  integer:: nlgn1(knln)          !! maximum order of moments + 1 (indp<=0;imthd=3) (N_stream x 2 + 1)
  real(8):: g(klgn1,knln)        !! Legendre moments of phase function (indp<=0)
  real(8),save:: epsp            !! truncation criterion of phase function moments
  real(8),save:: epsu            !! convergence criterion of intensity (inda>0)
  real(8):: fsol                 !! solar irradiance at the system top
  integer:: nplk1                !! number of order to approximate plank + 1
  real(8):: cplk(kplk1,knln)     !! coefficient of planck function 
  real(8):: bgnd                 !! planck function at ground 
  integer,save:: ntau            !! number of user defined optical depth
  real(8):: utau(kntau)          !! optical depths where the field is calculated
  real(8):: flxd(kna0,kntau)     !! downward flux
  real(8):: flxu(kna0,kntau)     !! upward flux
!  real(8):: ai(kna1u,kna0,knfi,kntau) !! intensity (not used for SDSU)

  real(8):: ai_sdsu(kntau)  !intensity output

! for initialization
  integer:: k,k1
  real(8):: rhl
  real(8):: dz,dzcm(knln),dzm(knln)

! for work
  integer:: iws,iww,ik

  integer,save:: nln1                 !! nln + 1
  integer,save:: ipb(knln1)           !! abs(ipbf)
  real(8),save:: gtmp                 !! ground temperature
  real(8):: pl(knln)             !! pressure in each layer
  real(8),save:: rh(knln)             !! relative humidity in each layer
  real(8),save:: phsm(knang)          !! Rayleigh phase function
  real(8):: wgt_res  !weigting function for response function

  real(8):: wla,wl1,wlc,wn       !! wavelenghs
  real(8):: sunir                !! function of solar irradiance

  real(8):: tray,taur(knln)      !! Rayleigh
  real(8):: tauc(knln)    ! condensate optical depth
  real(8):: taua(knln)    ! aerosol optical depth
  real(8),allocatable,dimension(:,:) :: phsf,phsfa,phsfc,phsfr ! phase function (indp=1) total,particle,rayleigh [str-1]
  real(8):: wgt_k(kch)  !weighting function of k-distribution

  real(8):: taumax
  real(8):: dpt
  real(8):: bplnk,b1,b2
  
  real(8):: wgt_k_res  ! joint weigting function = wgt_k * wgt_res
  real(8):: wgt_tot  !total weight = sum ( wgt_k_res)
  real(8):: sca2(knln)
  real(8):: thks

! for SDSU
  integer :: nx,ny,sl  ! horizontal model loop
  real :: w1, w2 ! linear weights
  real(8) :: cng_crm(knl,knm0) ! temporal gas concentrations interpolated for CRM pressure [ppmv]
  real(8) :: n0_fake(nbin) !fake n0 for initialization

  real(8) :: alpha(knln)  !lidar ratio
  real(8) :: sbacktot(knln) , sbacktot_p22(knln), sbacktot_clear(knln) !total backscattering coefficient [m-1 str-1]
  logical :: logic_opened  
  real(8) :: bcat(knln),bcat0(knln)  !lidar attenuating backscattering coef [m-1 str-1]
  real(8) :: bcat_p22(knln),bcat0_p22(knln) !perpendicular polarization of lidar attenuating backscattering coef [m-1 str-1]
  real(8) :: bcat_clear(knln),bcat0_clear(knln)  !clear-sky lidar attenuating backscattering coef [m-1 str-1]
  real(8) :: kext(knln), kext_clear(knln)  !total extinction [m-1]  
  real(8) :: solzen !solar zenith angle [deg]
  real(8) :: psum,phs
  integer :: ll,iw_r
  real(8) :: dhgt(knln),hgt_stag(knln)

  real(8) :: start_wave_number  !starting wavenumber [1/cm]
  real(8) :: end_wave_number    !ending wavenumber [1/cm]
  real(8) :: wave_number !wavenumber [1/cm]
 character(len=20) :: sensor_name

!----------------------------------------------------------------------------------------

!
! init
!
  nlgn1 = 0
  galb  = 0.d0
  cplk  = 0.d0
  utau  = 0.d0


 if (trim(sensor) == 'lidar') then
  att_B       = undefined 
  att_p22_B   = undefined
  lidar_sback = undefined 
  lidar_ratio = undefined 
  optical_depth = undefined
 endif

!--exec

!
! allocate memory
!
  if( .not. allocated(taug) ) then
      allocate( taug(1:knln,1:kch,1:kww) , stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank, 'MSG opt_visir: allocation error -> Terminate program. ')
  endif

  if( .not. allocated(phsf) ) then
      allocate( phsf(1:knang,1:knln)  , &
                phsfa(1:knang,1:knln) , &
                phsfc(1:knang,1:knln) , &
                phsfr(1:knang,1:knln) , &
                stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank, 'MSG opt_visir: allocation error -> Terminate program. ')
  endif



! ------------------- SDSU constant parameter --------------------
 nda   = N_stream ! # of stream (N_stream = 1 -> 1x2= 2 stream: for radiance at least 4 stream)
 imthd = 2        ! flux/intensity calculation modes (0-DM,1-MS,2-TMS,3-IMS) 
                  ! 2-TMS is good for reflected radiance.
 inda  = 1        ! radiance flag (0: flux only, 1: flux & radiance)
 na0 = 1            ! solar zenith angle 
 na1 = 1            ! satellite incident angle
 th1(1) = znth_obs  ! veiwing zenith angle [deg] from module_simulator (=sensor incident angle)
 nfi = 1            ! azimuthal angle
 fi(1) = azmth      ! azimuth angle between the sun and sensor [deg]
 nw0=1 ; rx=0 ; rf=1     ! monocromatic (line) radiation (modify here lator )

sensor_select0: select case(sensor)
case('visir')
     nw = mxwavel            ! # of wavelengths (channels)
     wl(1:nw) = wavel(1:nw)  ! center wavelength [micron] 
     sensor_name = visir_sensor
     if(masterproc) print*,'MSG opt_visir: Start VisibleIR Simulator '
     if(masterproc .and. verbose_SDSU) print*,'... Be Patient. Very slow, but running ...'

case('lidar')
     nw = mxwavel_lidar           ! # of wavelengths (channels)
     wl(1:nw) = wavel_lidar(1:nw) ! center wavelength [micron]
     sensor_name = lidar_sensor
     if(masterproc) print*,'MSG opt_visir: Start Lidar Simulator '
case default ; call stop_sdsu(myrank, 'MSG opt_visir: There is no such sensor')
end select sensor_select0

 dw(1:nw) = 0.           ! scaling factor to define scaled wavelengths (0.-monocromatic )   

 nln  = mxlyr  ! # of vertical layers (CRM layers) 
 nln1 = nln+1  ! nln+1
 do i=1,nln1
    ipbf(i)=nln1-i+1 ! profile index: +1 is removed since SDSU level is 0~mxlyr
 end do

 matm = 1      ! Sounding climatology index [1~6] 
               ! 1-tropics (toshii modify here lator for interpolation)

 wlcn  = 0. !scaling wavelngth [cm] (not used for SDSU) 
! ---------------- end of SDSU constant parameter ------------------

!
! allocate LUT parameter
!
 if( lut_visir  ) then

   !
   ! reading netcdf LUT
   !
   call read_lut_visir_nc(sensor_name,nw)

   if(account_aerosol) then
     call read_lut_visir_aero_nc(trim(sensor_name),nw)
   endif 

  endif !lut

!
! initialization of atmosphere profiles from MLATMD (each values are level (staggered))
!
  alt=0.d0 ; prs=0.d0 ; tmp=0.d0 ; cng=0.d0
  call read_sounding(matm,nl,alt,prs,tmp,cng)

!
! initialization of fundamental particle refractive index (rfi)
! outputs are stored in the saved global parameter in module paras
!
  call read_rfi

!
! prepare LUT of water surface reflectance matrix
!
  if( visir ) &
  call readlut_visir_grnd

!
! start horizontal loop 
!

 domain_loop: do ny = myj_start, myj_end  ; do nx = myi_start, myi_end



!
! if umu_grid present (<-- from scan simulator), use local incident angle
!
   if( present( umu_grid ) ) then
       if( umu_grid(nx,ny) == undefined ) cycle

       th1(1) = acos( umu_grid(nx,ny) )/const_degrad !veiwing zenith angle [deg] from module_scan (=sensor incident angle) 
   endif

!
! solar zenith angle for daytime simulation (2D loop lator)
!
     if(znth_slr /= -999.) then
          th0(1) = znth_slr ! solar zenith angle [deg] from module_simulator
     else
          th0(1) = acos( surface(nx,ny)%cosz )/ (pi/180.d0)
     endif

     if(th0(1) >= 0.d0 .and. th0(1) < 90.d0 ) then    !- daytime
        isol=1            ! (0:night, 1:day)
     else                             !- nighttime
        isol=0            ! (0:night, 1:day)
        th0(1)=0.
     endif

! remove this lator
     ifrh = 0         ! humidity flag (0: model, 1: use trh) -> not use
     trh = 0.98       ! relative humidity -> not used

!
! interpolate CO2,O3,N2O,CO,CH4,O2 conc [mmpv] by CRM pressure level [mb]
!
     cng_crm = 0. !initialize 

     do l=1,nln1 ! CRM level loop
        p1 = DBLE(atmos_stag(nx,ny,l-1)%press)       ! pressure [mb]

        do sl = 1,nl-1 ! sounding level loop
           if    ( p1 > prs(1) ) then               ! use bottom level of sounding
               cng_crm(l,2:kmol) =  cng(1 ,2:kmol)
               exit
           elseif( p1 < prs(nl) ) then              ! use top level of sounding
               cng_crm(l,2:kmol) =  cng(nl,2:kmol) 
               exit
           elseif( p1 <= prs(sl) .and. p1 >= prs(sl+1) ) then  ! interpolate
               w1 = ( p1-prs(sl+1) ) / ( prs(sl) - prs(sl+1) )       
               w2 = 1.- w1
               cng_crm(l,2:kmol) = w1*cng(sl,2:kmol) + w2*cng(sl+1,2:kmol) !interpolated gas conc at CRM level
               exit 
           endif
        enddo !sl

     enddo !nln1


!
! Replacing Sounding level parameters (alt,prs,tmp,cng) by CRM level parameters
!
     do l=1,nln1 ! CRM level profile loop
        alt(l) = DBLE(atmos_stag(nx,ny,l-1)%hgt)       ! altitude (interface) [km]
        prs(l) = DBLE(atmos_stag(nx,ny,l-1)%press)     ! pressure [mb]
        tmp(l) = DBLE(atmos_stag(nx,ny,l-1)%t_air)     ! temperature [K]
        cng(l,2:kmol) = cng_crm(l,2:kmol) ! gas concenration [ppmv]  (excepting H2O)
     enddo !nln1

!
! loop for wavelengths
!
     init=1
     iw_loop: do iw=1,nw

     iw_global = iw  ! module parameter (used in simulator_visir)


sensor_select1: select case(sensor)
case('visir')

! 
! surface flag (water or land)  & boundary conditions
!
     if(surface(nx,ny)%iland == 1 ) then !land

       indg = 0                                    !surface flag (0:Lambert surface (land) )

       if( wavel(nw) < DBLE(tb_or_rad) ) then  !solar radiation spectrum
           galb(iw) =  DBLE(albedo_visir(nx,ny,iw))  !modis spectrum albedo from module_surface
       else  ! IR spectrum

           if( emiss_visir(nx,ny,iw) == undefined ) then !treat as water
             indg  = 1                               !surface flag (1: Specular surface (Ocean) )
             galb(1:nw) = DBLE(surface(nx,ny)%u10m)  !U10 [m/sec] (indg=1)  -> BRDF albeo will be computed in routine
           else
              galb(iw) =  DBLE(1.e0-emiss_visir(nx,ny,iw))  !MODIS-IR BF emissivity database 
           endif

       endif

     elseif(surface(nx,ny)%iland == 2 ) then !water

       indg  = 1                                   !surface flag (1: Specular surface (Ocean) )

       galb(1:nw) = DBLE(surface(nx,ny)%u10m)  !U10 [m/sec] (indg=1)  -> BRDF albeo will be computed in routine

     else
       call stop_sdsu(myrank, 'MSG: opt_visir: No such ground type' )
     endif


case('lidar')
case default ; call stop_sdsu(myrank, 'MSG opt_visir: There is no such sensor')
end select sensor_select1


!
! init
!
       if(init==1) then 
          init=0
          if(inda.ge.1) inda = 1
          indp = 1     ! flag of phase function
          epsp = 0.d0 ; epsu = 0.d0  ! truncation criterion of phase function moments
                                     ! convergence criterion of intensity (inda>0)
       endif  !(init=1)

!
! level checking
!
        ipb(1:nln1)=abs(ipbf(1:nln1))

        ntau=count(ipbf(1:nln1)>0)

        if(ntau<=0) call stop_sdsu(myrank, 'no interface for calculating radiation field')
        if(ntau>kntau) call stop_sdsu(myrank, 'number of field calculation layers larger than kntau')


!toshii (I don't understand this parameter well...)
!        if(ntau==nln1) then
!           indt=2     ! flag of optical depth (1:utau, 2:dpt)
!        else
           indt=1     ! flag of optical depth (1:utau, 2:dpt)
!        endif

        gtmp= DBLE(surface(nx,ny)%t_skin) ! skin temperature [K]

        amtpb(1:kmol,1:nln)=0.d0  ! initialize gas conc.

        do l = 1 ,nln  !CRM layer loop

           k  = ipb(l)-1
           k1 = k+1                
!
! pressure, temperature, and relative humidity layer value
!
           pl(l) = DBLE(atmos(nx,ny,k)%press)       ! pressure [mb] (k=1 is the bottom layer) for sdp2
           tl(l) = DBLE(atmos(nx,ny,k)%t_air)       ! temperature [K] for abgask
           rh(l) = DBLE(atmos(nx,ny,k)%rh)/100.d0  ! RH [-] (0.~1.) for rfidxb

           dhgt(l) = DBLE(atmos(nx,ny,k)%dhgt)  !height [km]

           ind=1 !compute ppmv for a given tl,pl,rh  (relative humidity, 0-1)
           call wvcal(ind,pl(l),tl(l),rh(l),ppmv_h2o,gm3,e,gm3s,es)

           dz=abs(alt(k)-alt(k1)) ; dzcm(l)=dz*1.0d5      !dz [km] , dzcm [cm]
! Modifying dz for computing gas concentration via hydrostatic-ideal gas law
           dz=8.31451d0*273.15d0*(prs(k)-prs(k1))/(1013.25d0*9.8d0*28.8d0)

! gas amount (within layer loop) (level -> layer value)

           amtpb(1,l)=amtpb(1,l)+ppmv_h2o*dz*1.d-1  !for H2O

           do m=2,kmol  !CO2,O3,N2O,CO,CH4,O2 
              amtpb(m,l)=amtpb(m,l)+0.5d0*(cng(k,m)+cng(k1,m))*dz*1.d-1   ! cng is gas concentration [ppmv]
           enddo                                                            !1.d-1=1.d-6*1.d5

           amtpb(kmol1,l)=amtpb(1,l)**2/(amtpb(1,l)+dz*1.d5)  ! for kmol1  dummy

        enddo  ! l CRM layer loop

        do l = 1 ,nln+1  !CRM stagger layer loop
           k  = ipb(l)-1
           hgt_stag(l) = DBLE(atmos_stag(nx,ny,k)%hgt)  ! height at interface [km]
        enddo
!
! pressure index in each layer, required for abgask.
!
        call sdp2(nln,pl,dp,np)

!
! initialization
!

        am0(1:na0)= cos(th0(1:na0)*rad) !toshii
        am1u(1:na1)=-cos(th1(1:na1)*rad) !toshii - this negative depends on how user define th1

        sol=0.d0 ; wgt_tot=0.d0 ; wgt_k_res=0.d0 ; thk0(1:nln,1:10)=0.d0 ; thk0_inst(:,:) = 0.d0
        flxd0(1:na0,1:ntau)=0.d0 ; flxd00(1:na0,1:ntau)=0.d0
        flxu0(1:na0,1:ntau)=0.d0 ; ai0(1:na1,1:na0,1:nfi,1:ntau)=0.d0

        bcat0 = 0.d0  ; bcat0_clear = 0.d0 ; bcat0_inst = 0.d0 ; bcat0_clear_inst = 0.d0
        bcat0_p22 = 0.d0  ; bcat0_p22_inst = 0.d0 
       

!
! sub-wavelength loop start
!
        iws_loop: do iws=1,nw0

!
! compute wavelength
!
          wla=wl(iw)+dw(iw)*rx(iws) ![micron]
          if(wla<0.d0) then
             call stop_sdsu(myrank, 'wl+dw*rx hits a negative value')
          endif
          wlb=wl(iw)+dw(iw)*rx(iws+1) ![micron]
          wl1=0.5d0*(wla+wlb)  !wavelength [micron]
          wgt_res =rf(iws)   ! weighting function
          wlc=wl1*1.0d-4  ! wavelength [micron] -> [cm]

!
! get solar constant 
!
          if(wl1<5.d0) then  !solar radiation spectrum range...
             if(wlb==wla) then
                wn=1.d4/wl1
                fsol=sunir(wn)  !function 
             else
                call fsoli(wla,wlb,fsol)
                fsol=fsol/(wlb-wla)
             endif
          else ! IR spectrum
             fsol=0.d0
          endif

          if(isol==0) fsol=0 ! make sure

          if(wl1<=2.d0) then
             nplk1=0  ! # of order to approximate plank + 1
          else
             nplk1=2  ! # of order to approximate plank + 1
          endif

! initialize
          taua = 0.d0  ; omga = 0.d0
          tauc = 0.d0  ; omgc = 0.d0
          phsf=0.d0 ; phsfa=0.d0 ; phsfr=0.d0 ; phsfc=0.d0 

!
! Initialize mie routine (Get nang, ang, scr, sci)
!

          mic_select1: select case(trim(type_microphysics))
          case('GEN')

             call mie_gen_visir('ini',  1.0d0, 15.0d-0, 0.d0, 1.0d+3, wlc,rhl, &            !get nang
                                 nang, ang, scr, sci, cextp_fake,cabsp_fake,php_fake )

          case('SBM')

             n0_fake(1:nbin) = 100. ![1/cm4]
             call mie_sbm_visir('ini',dble(rad_sbm%liq),dble(drad_sbm%liq),&
                               dble(den_sbm%liq),1.0d-6,n0_fake,wlc,rhl, & !get nang
                       nang, ang, scr, sci, cextp_fake,cabsp_fake,php_fake )


          case('RAMS')

             ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  !get nang 
             call mie_rams_visir('ini',cfmas,pwmas,  1.0d-6,              &
                                  15.0d-0, 2.0d-0, wlc, rhl,              &
                                  nang, ang, scr, sci, cextp_fake,cabsp_fake,php_fake )

          case default 

             call mie_gen_visir('ini',  1.0d0, 15.0d-0, 0.d0, 1.0d+3, wlc,rhl, &            !get nang
                                 nang, ang, scr, sci, cextp_fake,cabsp_fake,php_fake )                    

          end select mic_select1


          if(account_aerosol) then

            select case(trim(aerosol_microphysics) )
            case('GOCART')

               call mie_gocart('ini',                              1.0d-6, wlc,rhl, &   !get nang
                                nang, ang, scr, sci, cextp_fake,cabsp_fake,php_fake )
            case('UCD')

               call mie_ucd('ini',                              1.0d-6, wlc,rhl, &   !get nang
                             nang, ang, scr, sci, cextp_fake,cabsp_fake,php_fake )

            end select

          endif


!
! Get particle (aerosols, clouds, precipitation) optical properties
!

          do l=1,nln  !vertical loop

             rhl=rh(l)  !relative humidity

             mic_select2: select case(trim(type_microphysics))
             case('GEN')

               !
               ! compute optical properties for different species in GCE microphysics
               !
               call opt_gen_visir( iw, wlc, rhl, q_gen(nx,ny,ipb(l+1)), re_gen(nx,ny,ipb(l+1)), &
                                   cextp_gen, cabsp_gen, php_gen)

               !
               ! get total particle optical properties
               !
               call total_opt_gen_visir( cextp_gen,cabsp_gen,php_gen,dzcm(l),nang,&
                                         tauc(l),omgc(l),phsfc(1:knang,l) )


             case('RAMS')

               nt = max( 1,min( 31,-nint(  (atmos(nx,ny,ipb(l+1))%t_air-273.15)  )))  ! Temperature index--1 to 31
               ns = max(79,min(100, nint( atmos(nx,ny,ipb(l+1))%rh          )))       ! RH index--79 to 100
               !
               ! compute optical properties for different species in RAMS microphysics
               !
               call opt_rams_visir(iw, nt, ns, wlc, rhl, q_rams(nx,ny,ipb(l+1)), re_rams(nx,ny,ipb(l+1)), &
                                   cextp_rams, cabsp_rams, php_rams)

               !
               ! get total particle optical properties
               !
               call total_opt_rams_visir( cextp_rams,cabsp_rams,php_rams,dzcm(l),nang,& 
                                          tauc(l),omgc(l),phsfc(1:knang,l) ) 


             case('SBM')

               !
               ! compute optical properties for different species in HUCM microphysics
               !
               call opt_sbm_visir(iw, wlc, rhl, q_sbm(nx,ny,ipb(l+1)), n_sbm(nx,ny,ipb(l+1),1:nbin), &
                                  cextp_sbm, cabsp_sbm, php_sbm)

               !
               ! get total particle optical properties
               !
               call total_opt_sbm_visir( cextp_sbm,cabsp_sbm,php_sbm,dzcm(l),nang,&
                                          tauc(l),omgc(l),phsfc(1:knang,l) )


             case default 

                 tauc(l) = 0.d0  ; omgc(l) = 0.d0  ;  phsfc(1:nang,l) = 0.d0

             end select mic_select2

            if(account_aerosol) then


              select case(trim(aerosol_microphysics) )
              case('GOCART')

               !
               ! compute optical properties for different species in HUCM microphysics
               !
               call opt_gocart_visir( iw, wlc, rhl, q_gocart(nx,ny,ipb(l+1)), &
                                      cextp_gocart, cabsp_gocart, php_gocart )

               !
               ! get total particle optical properties
               !
               call total_opt_gocart_visir( cextp_gocart,cabsp_gocart,php_gocart,dzcm(l),nang,&
                                            taua(l),omga(l),phsfa(1:knang,l) )


              case('UCD')

               !
               ! compute optical properties for different species in HUCM microphysics
               !
               call opt_ucd_visir( iw, wlc, rhl, q_ucd(nx,ny,ipb(l+1)), &
                                   cextp_ucd, cabsp_ucd, php_ucd )

               !
               ! get total particle optical properties
               !
               call total_opt_ucd_visir( cextp_ucd,cabsp_ucd,php_ucd,dzcm(l),nang,&
                                         taua(l),omga(l),phsfa(1:knang,l) )


              end select

            endif !aerosol

          enddo !nln


!
! Rayleigh Scattering for given ang & nang
! If you want to check pure phase function (not normalized by solid angle)
! calculate 4*pi*phsm (this is non-dimensional purely probability density)
!
          phsm(1:nang) = 3.d0/16.d0/pi*(1.d0+cos(ang(1:nang)*rad)**2)  !rayleigh phase function [str-1]

          tray=0.00864d0/wl1**(3.916d0+0.074d0*wl1+0.05d0/wl1)  !wl1 is wavlength in um

          do l=1,nln
             taur(l)=abs(prs(ipb(l))-prs(ipb(l+1)))/pstd*tray ! optical depth due to rayleigh scattering [-]
             omgr(l)=1.        ! single scattering albedo (zero absorption from Rayleigh scatter of molecule) [-]
             phsfr(1:nang,l)=phsm(1:nang)*taur(l)*omgr(l)    ! phase function of rayleigh scattering [str-1]
          enddo !l loop

!
! total phase function  (here defined as phsf## = ph*omega*tau)
!
          do l=1,nln

             !total            aerosol         cloud           rayleigh
             phsf(1:nang,l)= phsfa(1:nang,l)+phsfc(1:nang,l)+phsfr(1:nang,l)   ! phase function [str-1]

             if(sensor == 'lidar') then
                 dzm(l)      = dzcm(l)/100.                  ! layer depths [m]
                 kext(l)     = (taua(l) + tauc(l) + taur(l)) / dzm(l)  ! total extinction [m-1]
                 kext_clear(l)  = taur(l) / dzm(l)         ! clear-sky total extinction [m-1]
                 sbacktot(l) = phsf(nang,l)        / dzm(l)  ! total backscattering [m-1 str-1]
                 sbacktot_clear(l) = phsfr(nang,l) / dzm(l)  ! rayleigh total backscattering [m-1 str-1]
                 alpha(l)    = kext(l) / sbacktot(l)         ! lidar ratio [str]

                 ! note does not account polarization in Rayleigh scattering yet...(toshii)
                 sbacktot_p22(l) =  (phsfa(nang+1,l)+phsfc(nang+1,l)+phsfr(nang,l)) / dzm(l)  !perpendicular polarization
  
             endif

          enddo

! 
! get gas absorption : abgask (the k-distribution using hitran 2004)
!
          nww=int(log10(1.d4/wla)*dintvl+1.d0)-int(log10(1.d4/wlb)*dintvl+1.d0)+1

          call abgask(nln,nww,wlb,dp,np,amtpb,tl)
          call qgausn(gwt,gmu,nch_k)

          do ik=1,nch_k
             wgt_k(ik)=2.d0*gwt(ik)*gmu(ik)
          enddo

!
! abgask integration loop
!
           iww_loop: do iww=1,nww  

!
! k-dist file loop channel
!
              ik_loop: do ik=1,nch_k
!
! total optical properties
!
                 tau(1:nln)=taug(1:nln,ik,iww)+taua(1:nln)+tauc(1:nln)+taur(1:nln)  ! total optical depths
                 omg(1:nln)=( omga(1:nln)*taua(1:nln) + omgc(1:nln)*tauc(1:nln) + omgr(1:nln)*taur(1:nln) ) &
                             / tau(1:nln)                                ! single scattering albedo 


!
! filter optical depths for numerical stability 
! 
                 do l=1,nln
                    taumax=200.d0/(1.04d0-min(1.d0,omg(l)))  ! EMK omg become > 1.0 (check in details)
                    tau(l)=min(max(tau(l),1.0d-7),taumax)
                 enddo
! indt=1 case
                 if(indt==1) then
                    lf=0; dpt=0.d0
                    do l=1,nln1
                       if(ipbf(l)>0) then
                          lf=lf+1
                          utau(lf)=dpt
                       endif
                       dpt=dpt+tau(l)
                    enddo
                 endif

! plank function fitting (nplk1=3)
                 bgnd=0.d0
                 if(nplk1>0) then
                    b2=bplnk(wl1,tmp(ipb(1)))
                    do l=1,nln
                       b1=b2
                       b2=bplnk(wl1,tmp(ipb(l+1)))
                       cplk(1,l)=b1
                       cplk(2,l)=(b2-b1)/tau(l)
                    enddo
                    if(indg==0) then     ! lambert surfacce (land)
                       bgnd=(1.-galb(iw))*bplnk(wl1,gtmp)
                    else if(indg>0) then ! specular surface (water)
                       bgnd=bplnk(wl1,gtmp)
                    endif
                 endif

!
! Weighting functions for output
!
                 wgt_k_res = wgt_k(ik)*wgt_res
                 wgt_tot  = wgt_tot + wgt_k_res


!
! Sensor dependent
!
sensor_select2: select case(sensor)
 case('visir')

        galb_gwind = galb(iw) 
!
! Drive radiative transfer  (this was rtrn21)
!
!am1u = -1.0d0
                call simulator_visir(indg,inda,indt,indp,imthd,nda,na1,am1u,na0,am0, &
                     nfi,fi,nln,tau,omg,nlgn1,g,nang,ang,phsf,epsp,epsu,galb_gwind,&
                     fsol,nplk1,cplk,bgnd,ntau,utau,scr,sci,flxd,flxu,ai_sdsu)

!
! Output (06.03.29 miho)
!
                     sol  = sol + fsol* wgt_k_res
                     thk0(1:nln,1)=thk0(1:nln,1)+tau (1:nln)*wgt_k_res  !total tau
                     thk0(1:nln,2)=thk0(1:nln,2)+(taua(1:nln)+tauc(1:nln))*wgt_k_res  !particle tau
                     thk0(1:nln,3)=thk0(1:nln,3)+taur(1:nln)*wgt_k_res  !rayleigh tau
                     sca2(1:nln)=omg(1:nln)*tau(1:nln)*wgt_k_res     
                     thk0(1:nln,4)=thk0(1:nln,4)+sca2(1:nln)            !single scattering albedo
                     thk0(1:nln,5)=thk0(1:nln,5)+g(2,1:nln)*sca2(1:nln) !g1 (1st legendre polynomial moment)
                     thk0(1:nln,6)=thk0(1:nln,6)+g(3,1:nln)*sca2(1:nln) !g2 (2nd legendre polynomial moment)

                     lf=0
                     thks=0.d0
                     do l=1,nln1
                        if(l>=2) thks=thks+tau(l-1)
                        if(ipbf(l)>0) then
                           lf=lf+1
                           flxd0(1:na0,lf)=flxd0(1:na0,lf)+flxd(1:na0,lf)*wgt_k_res
                           flxd00(1:na0,lf)=flxd00(1:na0,lf) &
                                            +fsol*am0(1:na0)*exp(-thks/am0(1:na0))*wgt_k_res
                           flxu0(1:na0,lf)=flxu0(1:na0,lf)+flxu(1:na0,lf)*wgt_k_res
                           if(inda/=0) then
                             ai0(1:na1,1:na0,1:nfi,lf) = ai0(1:na1,1:na0,1:nfi,lf) &
!                                                      + ai(1:na1,1:na0,1:nfi,lf)*wgt_k_res
                                                       + ai_sdsu(lf)*wgt_k_res
                           endif
                        endif
                     enddo


case('lidar')

                !
                ! for diagnostic
                !
                thk0(1:nln,1)=thk0(1:nln,1)+ tau (1:nln)*wgt_k_res     !total tau
                thk0(1:nln,2)=thk0(1:nln,2)+ (taua(1:nln)+tauc(1:nln))*wgt_k_res  !particle tau
                thk0(1:nln,3)=thk0(1:nln,3)+ taur(1:nln)*wgt_k_res     ! rayleigh tau
                thk0(1:nln,4)=thk0(1:nln,4)+sbacktot(1:nln)*wgt_k_res ! total backscattering  [m-1 str-1]
                thk0(1:nln,5)=thk0(1:nln,5)+ alpha(1:nln)*wgt_k_res    ! lidar ratio

             if( inst_profile_lidar ) then  !instrumental vertical cordinate

               !
               ! interpolate/extrapolate vertical cordinate into instrumental level
               !   
               call refine_range_bin( nln, hgt_stag(1:nln+1), dhgt(1:nln), &
                                       mxlyr_lidar,  hgt_stag_lidar, &
                                       kext(1:nln), kext_inst(1:mxlyr_lidar)  &
                                      )

               call refine_range_bin( nln, hgt_stag(1:nln+1), dhgt(1:nln), &
                                       mxlyr_lidar,  hgt_stag_lidar, &
                                       kext_clear(1:nln), kext_clear_inst(1:mxlyr_lidar)  &
                                      )

               call refine_range_bin( nln, hgt_stag(1:nln+1), dhgt(1:nln), &
                                       mxlyr_lidar,  hgt_stag_lidar, &
                                       sbacktot(1:nln), sbacktot_inst(1:mxlyr_lidar)  &
                                      )

               call refine_range_bin( nln, hgt_stag(1:nln+1), dhgt(1:nln), &
                                       mxlyr_lidar,  hgt_stag_lidar, &
                                       sbacktot_p22(1:nln), sbacktot_p22_inst(1:mxlyr_lidar)  &
                                      )

               call refine_range_bin( nln, hgt_stag(1:nln+1), dhgt(1:nln), &
                                       mxlyr_lidar,  hgt_stag_lidar, &
                                       sbacktot_clear(1:nln), sbacktot_clear_inst(1:mxlyr_lidar)  &
                                      )

                do k_inst = 1, mxlyr_lidar
                   k_rev = mxlyr_lidar - k_inst + 1
                     dzm_inst(k_rev) = dhgt_lidar(k_inst) * 1000.  !layer depth [m]
                enddo 

               !
               ! compute attenuating back scatter profiles
               !
                call simulator_lidar(mxlyr_lidar,dzm_inst,kext_inst,sbacktot_inst,bcat_inst) 

                call simulator_lidar(mxlyr_lidar,dzm_inst,kext_inst,sbacktot_p22_inst,bcat_p22_inst)

                call simulator_lidar(mxlyr_lidar,dzm_inst,kext_clear_inst,sbacktot_clear_inst,&
                                     bcat_clear_inst)  !clear sky only (Rayleigh atmosphere)

               !
               ! attenuating backscattering [m-1 str-1]
               !
               do k_inst = 1,mxlyr_lidar 
                   if( bcat_inst(k_inst) /= DBLE(undefined) ) then
                       bcat0_inst(k_inst) =bcat0_inst(k_inst) + bcat_inst(k_inst)*wgt_k_res    
                       bcat0_p22_inst(k_inst) =bcat0_p22_inst(k_inst) + bcat_p22_inst(k_inst)*wgt_k_res
                   else
                       bcat0_inst(k_inst) = DBLE(undefined)
                       bcat0_p22_inst(k_inst) = DBLE(undefined)
                   endif

                   if( bcat_clear_inst(k_inst) /= DBLE(undefined) ) then
                       bcat0_clear_inst(k_inst) = bcat0_clear_inst(k_inst) + bcat_clear_inst(k_inst)*wgt_k_res
                   else
                       bcat0_clear_inst(k_inst) = DBLE(undefined)
                   endif

               enddo 


           
             else !CRM level

               !
               ! compute attenuating back scatter profiles
               !
                call simulator_lidar(nln,dzm,kext,sbacktot,bcat)

                call simulator_lidar(nln,dzm,kext,sbacktot_p22,bcat_p22)

                call simulator_lidar(nln,dzm,kext,sbacktot_clear,bcat_clear)

               !
               ! output
               !
                bcat0(1:nln) =bcat0(1:nln) + bcat(1:nln)*wgt_k_res     ! attenuating backscattering [m-1 str-1]
                bcat0_p22(1:nln) =bcat0_p22(1:nln) + bcat_p22(1:nln)*wgt_k_res !perpendicular polariztion component [m-1 str-1]
                bcat0_clear(1:nln) =bcat0_clear(1:nln) + bcat_clear(1:nln)*wgt_k_res     ! attenuating backscattering [m-1 str-1]

             endif

case default ; call stop_sdsu(myrank, 'MSG opt_visir: There is no such sensor')
end select sensor_select2

              enddo ik_loop !end of k-dist file loop channel

           enddo iww_loop !end of abgask integration loop

       enddo iws_loop  !end of sub-wavelength loop

!
! Sensor Dependent
!
sensor_select3: select case(sensor)
case('visir')

!
! output (total weighting)
!
           sol=sol/wgt_tot
           thk0(1:nln,1)=thk0(1:nln,1)/wgt_tot
           thk0(1:nln,2)=thk0(1:nln,2)/wgt_tot
           thk0(1:nln,3)=thk0(1:nln,3)/wgt_tot
           thk0(1:nln,4)=thk0(1:nln,4)/wgt_tot/thk0(1:nln,1)
           thk0(1:nln,5)=thk0(1:nln,5)/wgt_tot/thk0(1:nln,1)
           thk0(1:nln,6)=thk0(1:nln,6)/wgt_tot/thk0(1:nln,1)

           lf=0
           do l=1,nln1
              if(ipbf(l)>0) then
                 lf=lf+1
                 flxd0(1:na0,lf)=flxd0(1:na0,lf)/wgt_tot
                 flxd00(1:na0,lf)=flxd00(1:na0,lf)/wgt_tot
                 flxu0(1:na0,lf)=flxu0(1:na0,lf)/wgt_tot
                 if(inda/=0) then
                   ai0(1:na1,1:na0,1:nfi,lf)=ai0(1:na1,1:na0,1:nfi,lf)/wgt_tot
                 endif
              endif
           enddo

!
! for SDSU global parameter 
!
           if(wl1 < DBLE(tb_or_rad)) then ! radiance output
                radiance_out(nx,ny,iw) = ai0(1,1,1,1)    ! radiance [w/m2/micron/str]
           else !Tb output
                radiance_out(nx,ny,iw) =  1.438786d+4/wl1/log(1.1911d+8/wl1**5/ai0(1,1,1,1)+1.) ! Tb [K]
           endif
!           print*,'radiance',nx,iw, radiance_out(nx,ny,iw) !,qcol_rams(nx,ny)%ice1,qcol_rams(nx,ny)%ice2,qcol_rams(nx,ny)%snow


case('lidar')

           !
           ! output onto the global module parameter
           !
           thk0(1:nln,1)=thk0(1:nln,1)/wgt_tot  !total tau
           thk0(1:nln,2)=thk0(1:nln,2)/wgt_tot  !particle tau
           thk0(1:nln,3)=thk0(1:nln,3)/wgt_tot  !rayleigh tau
           thk0(1:nln,4)=thk0(1:nln,4)/wgt_tot  !total backscattering  [m-1 str-1]
           thk0(1:nln,5)=thk0(1:nln,5)/wgt_tot  !lidar ratio

           if( inst_profile_lidar ) then  !instrumental vertical cordinate

             do n_diag = 1, 5
                call refine_range_bin( nln, hgt_stag(1:nln+1), dhgt(1:nln), &
                                       mxlyr_lidar,  hgt_stag_lidar, &
                                       thk0(1:nln,n_diag), thk0_inst(1:mxlyr_lidar,n_diag)  &
                                       )
             enddo

             do k_inst = 1,mxlyr_lidar
                if( bcat_inst(k_inst) /= DBLE(undefined) ) then
                    bcat0_inst(k_inst) =bcat0_inst(k_inst)/wgt_tot
                    bcat0_p22_inst(k_inst) =bcat0_p22_inst(k_inst)/wgt_tot
                else
                    bcat0_inst(k_inst) = DBLE(undefined)
                    bcat0_p22_inst(k_inst) = DBLE(undefined) 
                endif

                if( bcat_clear_inst(k_inst) /= DBLE(undefined) ) then
                    bcat0_clear_inst(k_inst) =bcat0_clear_inst(k_inst)/wgt_tot
                else
                    bcat0_clear_inst(k_inst) = DBLE(undefined)
                endif

            enddo

            do k_inst = 1, mxlyr_lidar
               l = mxlyr_lidar - k_inst + 1
               if( bcat0_inst( k_inst  ) < 0. ) then
                  att_B      (nx,ny,l,iw) = undefined  ! total attenuating backscattering [km-1 str-1]
                  att_p22_B  (nx,ny,l,iw) = undefined  ! perpendicular polarization component [km-1 str-1]
                  att_clear_B(nx,ny,l,iw) = undefined  ! clear-sky attenuating backscattering [km-1 str-1]
                  lidar_sback(nx,ny,l,iw) = undefined  ! total backscattering  [km-1 str-1]
                  lidar_ratio(nx,ny,l,iw) = undefined  ! lidar ratio
                  optical_depth(nx,ny,l,iw)= undefined
               else
                  att_B      (nx,ny,l,iw) = REAL(bcat0_inst( k_inst  )*1000.d0) ! total attenuating backscattering [km-1 str-1]
                  att_p22_B  (nx,ny,l,iw) = REAL(bcat0_p22_inst( k_inst  )*1000.d0) ! perpendicular polarization component 
                                                                                    ! of total attentuaing backscattering [km-1 str-1]
                  att_clear_B(nx,ny,l,iw) = REAL(bcat0_clear_inst( k_inst  )*1000.d0) ! clear-sky attenuating backscattering [km-1 str-1]
                  lidar_sback(nx,ny,l,iw)  = REAL(thk0_inst ( k_inst,4)*1000.d0) ! total backscattering  [km-1 str-1]
                  lidar_ratio(nx,ny,l,iw)  = REAL(thk0_inst ( k_inst,5)        ) ! lidar ratio
                  optical_depth(nx,ny,l,iw)= REAL(thk0_inst ( k_inst,2)        ) ! particle optical depth [-]
               endif

            enddo


           else 
             ! 
             ! normalized by total weight. 
             !

             bcat0(1:nln) =bcat0(1:nln)/wgt_tot    !attenuating backscattering [m-1 str-1]  
             bcat0_p22(1:nln) =bcat0_p22(1:nln)/wgt_tot  !perpendicular polarization attenuating backscattering [m-1 str-1] 
             bcat0_clear(1:nln) =bcat0_clear(1:nln)/wgt_tot    !attenuating backscattering [m-1 str-1]  

             !
             ! for SDSU global parameter (reverse profiles)
             !           
             do l=1,nln
                if( bcat0( ipb(l+1) )  < 0.) then
                  att_B      (nx,ny,l,iw) = undefined ! total attenuating backscattering [km-1 str-1]
                  att_p22_B  (nx,ny,l,iw) = undefined  ! perpendicular polarization component [km-1 str-1]
                  att_clear_B      (nx,ny,l,iw) = undefined ! total attenuating backscattering [km-1 str-1]
                  lidar_sback(nx,ny,l,iw) = undefined ! total backscattering  [km-1 str-1]
                  lidar_ratio(nx,ny,l,iw) = undefined ! lidar ratio
                  optical_depth(nx,ny,l,iw)= undefined! optical depth [-]
                else
                  att_B      (nx,ny,l,iw) = bcat0( ipb(l+1) )   *1000.d0 ! total attenuating backscattering [km-1 str-1]
                  att_p22_B  (nx,ny,l,iw) = bcat0_p22( ipb(l+1) )   *1000.d0 !perpendicular polarization component [km-1 str-1]
                  att_clear_B(nx,ny,l,iw) = bcat0_clear( ipb(l+1) )*1000.d0 ! clear-sky attenuating backscattering [km-1 str-1]
                  lidar_sback(nx,ny,l,iw) = thk0 ( ipb(l+1) ,4) *1000.d0 ! total backscattering  [km-1 str-1]
                  lidar_ratio(nx,ny,l,iw) = thk0 ( ipb(l+1) ,5)          ! lidar ratio
                  optical_depth(nx,ny,l,iw)= thk0( ipb(l+1), 2)          ! optical depth [-]
                endif
             enddo
           endif

case default ; call stop_sdsu(myrank, 'MSG opt_visir: There is no such sensor')
end select sensor_select3


     enddo iw_loop !end of channels wavelength loop

 enddo ; enddo domain_loop !end of domain nx ny loop


#if MPI == 2
if(numproc_tot>1) then
!
! let's send each tile into the master processor. 
!
 sensor_select4: select case(sensor)
 case('visir')
   if(numproc_tot>1) call mpi_sdsu_communicate( 'TO_MASTER',radiance_out )
 case('lidar')
   if(numproc_tot>1) call mpi_sdsu_communicate( 'TO_MASTER',att_B )
   if(numproc_tot>1) call mpi_sdsu_communicate( 'TO_MASTER',att_p22_B )
   if(numproc_tot>1) call mpi_sdsu_communicate( 'TO_MASTER',lidar_sback )
   if(numproc_tot>1) call mpi_sdsu_communicate( 'TO_MASTER',lidar_ratio )
   if(numproc_tot>1) call mpi_sdsu_communicate( 'TO_MASTER',optical_depth )

 case default ; call stop_sdsu(myrank, 'MSG opt_visir: There is no such sensor')
 end select sensor_select4

endif
#endif


 return
 end subroutine opt_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine read_rfi
! get parameters of refractive index of fundamental particle spceies
!--- history
! 95. 1.25  created from getph3
! 95. 1.24  dintpm -> aerd (now includs also aerosol water uptake data)
! 95  9.13  with modified aerd file format
! 96. 3.11  added non-spherical parameters
!     5. 5  awcrp(*,*,2)
! 08.10.31  modified fortran 90 version.
!--- input
! ius       i         read unit number of particle parameter file
!--- output
! nv        i         number of fundamental species (1-8)
!                     1 water
!                     2 ice
!                     3 dust-like
!                     4 water-soluble
!                     5 sea salt
!                     6 soot
!                     7 volcanic-ash
!                     8 75%H2SO4
!                     9 Yellow sand
!
! nwlv      i         number of wavelengths for refractive index tables
! wlv    r(kwl)       wavelengths (micron)
! rfi   r(kwl,knv,2)  refractive index of fundamental materials (mr, mi)
!                     =mr - i mi
!                     with log-regular wavelength interval
! 
!---
 ! use paras,only: kptc,knl,kaw,kwlv,knv,nv,nwlv,wlv,rfi

! work
 character(1)::ch
 integer::i,j,k,l,m,iv
 real(8)::rmin,rmax,nmode
 integer,parameter :: ius = 15
 character(len=200) :: path_file
 logical :: logic_opened
!---

 path_file=trim(sdsu_dir_data)//'rfi_table.txt'
 inquire(file=trim(path_file),opened=logic_opened)
 if( .not. logic_opened) then

    if(masterproc .and. verbose_SDSU) print*, 'MSG read_rfi: initialize refractive index of fundamental species'
!    open(ius,file=path_file,status='old')
    call open_oldfile(ius, trim(path_file) )

    rewind ius
!  number concentration profiles (relative unit)
!  number is conservative in growing aerosols
    read(ius,'(a1)') ch   !comments
    read(ius,*) rmin,rmax  !minimum & maximum radius

! get refractive index table (This is necessary).
    read(ius,'(a1)') ch     !comments
    read(ius,*) nv,nwlv     !nv is # of fundamental species, nwlv is # of wavelength for deining refractive index.
    read(ius,*) wlv(1:nwlv) !wavelength
    do iv=1,nv
       read(ius,'(a1)') ch          !comments
       read(ius,*) rfi(1:nwlv,iv,1) !refractive index (mr)
       read(ius,*) rfi(1:nwlv,iv,2) !refractive index (mi)  =mr - i mi
    enddo

 endif 

!  close(ius)

  return
end subroutine read_rfi

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine read_sounding(matm,nl_mlat,alt,prs,tmp,cng)
                 
! initilization of vertical profiles of gaseous conc from MLATMD
!--- history
! 95. 9.20  created
! 96. 3.11  kpcl=20 for more aerosol models and non-spherical parameters
!           getpr1 included
!     5. 5  drop npoly and get all the polydispersion information
! 01.10.14  changed for fortran90 version
! 08.09.20  modified free-style form
!--- input
! iua        i        device number for atmospheric model file
! matm       i        model number of atmosphere (afgl)
!                      1: tropical
!                      2: mid-latitude summer    3: mid-latitude winter
!                      4: high-lat summer        5: high-lat winter
!                      6: us standard
!--- output
! nl         i        number of atmospheric layers
! alt     r(knl)      height at interfaces of layers (km), top to bottom
! prs     r(knl)      pressure    at interfaces of layers (mb)
! tmp     r(knl)      temperature at interfaces of layers (k)
!---Declarations
 ! use paras,only: kptc,knl,knm0,kaw,kwlv,knv
  implicit none

! input
  integer,intent(in):: matm      !! atmospheric number

! output
  integer,intent(out):: nl_mlat  !! number of atmospheric layer
  real(8),intent(out):: alt(knl) !! altitude [km]
  real(8),intent(out):: prs(knl) !! pressure [mb]
  real(8),intent(out):: tmp(knl) !! temperature [K]
  real(8),intent(out):: cng(knl,knm0) !! gas concentration [ppmv]

! for mlatm
  integer,parameter:: knm1=7     !! number of major gas
  integer,parameter:: knm2=21    !! number of minor gas
  integer,parameter:: katm=6     !! number of standard atmospheric condition
  integer,parameter:: knm=knm1+knm2
  character(8)::idm(knm)
  integer,save::idms(knm,10),nl
  real(8),save::wmol(knm,10),rams(knm,10)
  real(8),save::pmatm(knl,katm),tmatm(knl,katm),dnsty(knl,katm)
  real(8),save::amol(knl,knm1,katm),trac(knl,knm2)

! input
  integer,parameter ::iu=13

!--- work area
  character ch*1
  integer::nm1,nm2,natm,nm,i,j,k,l,m, ns
  real(8)::airm
  logical :: logic_opened
  character:: fnm*100     !! file name

!
! initilization start
!

  fnm=trim(sdsu_dir_data)//'MLATMD'   
  inquire(file=trim(fnm),opened=logic_opened)

  if( .not. logic_opened)  then
    if(masterproc .and. verbose_SDSU) print*,'MSG mlatm: open/read MLATMD'

!    open(iu,file=fnm,status='old')
    call open_oldfile(iu, trim(fnm) )

! read numbers
    read(iu,'(a1)') ch
    read(iu,*) nm1,nm2,natm,nl
    nm=nm1+nm2

! read molecular information
    read(iu,'(a1)') ch
    read(iu,*) airm
    read(iu,'(a1)') ch
    idms(:,:)=0
    rams(:,:)=0.0d0
    do i=1,nm
       read(iu,'(a8)') idm(i)
       read(iu,*) ns
       read(iu,*) idms(i,1:ns)
       read(iu,*) wmol(i,1:ns)
       read(iu,*) rams(i,1:ns)
    enddo

! read altitude
    read(iu,'(a1)') ch
    read(iu,*) alt(1:nl)

! read pressure
    do j=1,natm
       read(iu,'(a1)') ch
       read(iu,*) pmatm(1:nl,j)
    enddo

! read temperature
    do j=1,natm
       read(iu,'(a1)') ch
       read(iu,*) tmatm(1:nl,j)
    enddo

! read molecular profiles
    do k=1,nm1
       do j=1,natm
          read(iu,'(a1)') ch
          read(iu,*) amol(1:nl,k,j)
       enddo
    enddo

    do j=1,natm
       read(iu,'(a1)') ch
       read(iu,*) dnsty(1:nl,j)
    enddo

! read trace gase profiles
    do  k=1,nm2
       read(iu,'(a1)') ch
       read(iu,*) trac(1:nl,k)
    enddo

 endif !logic_opened

!
! initialization end
!


!
! transfer to out parameter. (toshii: Need interpolation of matm depends on seasons & latitude)
!
  nl_mlat = nl
  prs(1:knl)=pmatm(1:knl,matm)               !pressure
  tmp(1:knl)=tmatm(1:knl,matm)               !temperature
  cng(1:knl,1:knm1)=amol(1:knl,1:knm1,matm)  !major gas only 1~7
                                             !trac (trace gasues) is not used...

  !nmol=knm1 !number of gaseus species
  return
end subroutine read_sounding

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine wvcal(ind,p,t,rh,ppmv,gm3,e,gm3s,es)
! conversion of units for water vapor
! ppmv =n/na =e/(p-e)*1.0e6, rh =e/es =gm3/gm3s
! ev = nrt, (p-e)v = na rt
!--- history
! 95. 9.30  created
! 08. 9.20  modified Fortran 90 free-style form
!--- input
! ind   i    1: give rh  (relative humidity, 0-1)
!            2: give ppmv (volume mixing ratio in ppmv)
!            3: give gm3 (mass mixing ratio in g/m3)
!            4: give e    (water vapor pressure in hpa)
! p     r    atmospheric pressure (hPa)
! t     r    temperature (k)
! one of rh, ppmv, gm3, or e
!--- output
! other than input variable among rh, ppmv and gm3.
! gm3s  r    saturation water vapor content (g/m3)
! es    r    saturation water vapor pressure (hPa)
!--- gas constant (cgs) and molecular weights of water vapor
!           and air (g/mol)
!--
  implicit none
! for input & output
  integer,intent(in)::ind        !! indicator
  real(8),intent(in)::p,t        !! pressure and temperature
  real(8),intent(inout)::rh,ppmv,gm3,e,gm3s,es !! quantities of watere vaour

! work
  real(8):: gcm3,gcm3s
  real(8):: wvsat                !! function for saturated
  real(8),parameter::r=8.314d7, w=18.02d0, air=28.964d0
!--exec
! satuaration water vapor (g/cm3)
  gcm3s=wvsat(t)*1.0d-6
  es=gcm3s/w*r*t/1.0d3
  if(ind==1) then
     gcm3=rh*gcm3s
     e=rh*es
     ppmv=e/(p-e)*1.0d6
  elseif(ind==2) then
     e=p*ppmv/(1.0d6+ppmv)
     rh=e/es
     gcm3=rh*gcm3s
  elseif(ind==3) then
     gcm3=gm3*1.0d-6
     rh=gcm3/gcm3s
     e=rh*es
     ppmv=e/(p-e)*1.0d6
  elseif(ind==4) then
     rh=e/es
     gcm3=rh*gcm3s
     ppmv=e/(p-e)*1.0d6
  endif
  gm3=gcm3*1.0d6
  gm3s=gcm3s*1.0d6
  if(rh>1.d0) rh=1.d0
  return
end subroutine wvcal


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine abgask(nln,nww,wlb,dp,np,amtpb,t)
! calculate absorption coefficients.
!--- history
! 06.03.17 created (miho)
! 07.11.30 fix dimension (nln => knln, nch => kch) 
! 08.08.14 add dintvl
! 08. 9.26 modified Fortran 90 free-style form
!--- input
! iug      i     read unit number of gas line absorption para file
! nln      i     number of layers
! wla      r     the lower boundary of wavelength
! wlb      r     the upper boundary of wavelength
! dp     r(nln)  residual of the pressure   
! np     i(nln)  grid number of the pressure
! t      r(nln)  temperature of the layer
!--- output
! taug  r(nln,kch,kww) optical thickness (global parameter in module paras)
!---
  !use paras, only: kmol,kp,kt,nch_k,kmol1,kww,knln,kch, taug, dintvl
  implicit none
! input
  integer,intent(in):: nln       !! device number
  integer,intent(in):: nww       !! device number
  real(8),intent(in):: wlb       !! device number
  real(8),intent(in):: dp(knln)  !! device number
  integer,intent(in):: np(knln)  !! device number
  real(8),intent(in):: amtpb(kmol1,knln) !! gas concentration
  real(8),intent(in):: t(knln)   !! temperature in layer


! work
  integer:: irec,iww,l,ich,imol,it
  real(4):: ak(kch,kp,kt,kmol1)
  real(8):: akt(kt,knln),knu(knln)

  integer,parameter:: iug=16   !! file number
  character(len=200) :: fnm
  logical :: logic_opened 

!--- exec
  taug(1:nln,1:nch_k,1:nww)=0.d0

  fnm=trim(sdsu_dir_data)//'ckd.g.ch_2'   
  inquire(file=trim(fnm),opened=logic_opened)
  if( .not. logic_opened) then 
     open(iug,file=fnm,status='old',access='direct',form='unformatted', &
          recl=4*kp*kt*nch_k*kmol1)
     if(masterproc .and. verbose_SDSU) print*,'MSG abgask: k-ditribution table (',trim(fnm),') is opened.'
  endif

! find the position of the beginning of the data.
! wlb = 10+10**(irec/dintvl)
  irec=int((log10(1.d4/wlb)-1.d0)*dintvl)
  do iww=1,nww
     irec=irec+1
     read(iug,rec=irec) ak(1:nch_k,1:kp,1:kt,1:kmol1)
     if(irec.lt.1) then
        call stop_sdsu(myrank, 'wavenumber out of the range [0.2-1000 micron]')
        return
     endif

     do ich=1,nch_k
        do imol=1,kmol1
           do l=1,nln
              do it=1,kt
                 akt(it,l)=ak(ich,np(l)-1,it,imol) &
                      +(ak(ich,np(l),it,imol)-ak(ich,np(l)-1,it,imol))*dp(l)
              enddo
           enddo
           call tdok2(nln,akt,t,knu)
           taug(1:nln,ich,iww)=taug(1:nln,ich,iww)+knu(1:nln)*amtpb(imol,1:nln) !this is global para in module paras
        enddo
     enddo
  enddo
  return
end subroutine abgask

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine sdp2(nln,p,dp,np)
! derive absorption coefficient from pt table and fitting.
!--history
! 02.02.27  Modified knupt in ckdng.f (Zhang Hua)
! 02.11.25  Modified from knpt1.      (Miho Sekiguchi)
! 03.09.09  Modified from sdp for free source form.(Miho Sekiguchi)
! 07.03.28  Modified about dimension nln => knln. (Miho Sekiguchi)
! 08. 9.25 modified Fortran 90 free-style form
!--
  !use paras, only: kp,knln
  implicit none

! input & output
  integer,intent(in):: nln       !! number of layer
  real(8),intent(in):: p(knln)   !! pressure in layer
  real(8),intent(out):: dp(knln) !! ratio of pressure grid
  integer,intent(out):: np(knln) !! number of pressure grid

! works
  integer:: l,ip
  real(8):: prs(1:26)=(/ &
       1.0130d+3,  0.63096d+3, 0.39811d+3,  0.25119d+3,  0.15849d+3, &
       1.0000d+2,  0.63096d+2, 0.39811d+2,  0.25119d+2,  0.15849d+2, &
       1.0000d+1,  0.63096d+1, 0.39811d+1,  0.25119d+1,  0.15849d+1, &
       1.0000d+0,  0.63096d+0, 0.39811d+0,  0.25119d+0,  0.15849d+0, &
       1.0000d-1,  0.63096d-1, 0.39811d-1,  0.25119d-1,  0.15849d-1, &
       1.0000d-2/)               !! pressure grid
!--exec
  do l=1,nln
     if(p(l) < prs(kp)) then
        np(l)=kp
        dp(l)=1.d0
        cycle
     endif
     do ip = 1, kp
        if (p(l) >= prs(ip)) exit
     enddo
     np(l)=ip
     if(np(l) > kp) np(l)=kp
     if(np(l) <= 1) np(l)=2
     dp(l)=log10(p(l)/prs(np(l)-1))/log10(prs(np(l))/prs(np(l)-1))
  enddo
end subroutine sdp2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine tdok2(nln,akt,t,knu)
!----------------------------------------------------------------------c
!                using Shi's formula: (t/to)**(a+bt)
!----------------------------------------------------------------------c
  !use paras, only: kt,knln
  implicit none
  integer,intent(in):: nln       !! number of layer

! 09.05.18 Bug fixed (Miho, reported by Hiro Masunaga)
  real(8),intent(inout):: akt(kt,knln)   !! absorption coefficients at T grid
  real(8),intent(in):: t(knln)   !! temperature in layer
  real(8),intent(out):: knu(knln) !! absorption coefficients
! works
  integer:: l
  real(8):: a,b,at1,at2,bt1,bt2
!--exec
  at1 = log10(200.d0/260.d0); at2 = log10(320.d0/260.d0)
  
  do l=1,nln
     if(akt(2,l) <= -14.d0.and.akt(1,l) > -14.d0.and. &
          akt(3,l) > -14.d0) akt(2,l)=(akt(1,l)+akt(3,l))*0.5d0
     bt1 = akt(1,l)-akt(2,l)
     bt2 = akt(3,l)-akt(2,l)
     if(bt1 == 0.d0 .and. bt2 == 0.d0) then
        knu(l)=10**akt(2,l); cycle
     endif
     b   = (bt2/at2-bt1/at1)/120.d0
     a   = bt2/at2 - b*320.d0
     
     knu(l)=10.d0**(akt(2,l))*(t(l)/260.d0)**(a+b*t(l))
  enddo
  return
 end subroutine tdok2


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_gen_visir(iw, wlc, rhl, q, re, cextp_gen, cabsp_gen, php_gen)
 implicit none

 integer,intent(in) :: iw  !index of wavelength loop
 real(8),intent(in) :: wlc ! wavelength [cm]
 real(8),intent(in) :: rhl ! relativev humidity  (0~1)
 type ( particle_gen ),intent(in) :: & ! particle_gen is defined in module_simulater
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 type ( particle_gen_r8 ),intent(out) :: & !particle_gen is defined in module_simulater
      cextp_gen,  & ! !extinction cross section [cm2/cm3] or [cm-1]
      cabsp_gen,  & ! absorption cross section [cm2/cm3] or [cm-1]
      php_gen(knang)! volume scattering phase function [cm-1  str-1]


! local parameters
 integer :: ispc       ! RAMS spicies index
 integer :: nang       ! number of scaattering angles
 real(8) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(8) :: ang(knang) ! scattering angles [degree]
 real(8):: scr         ! relative water reflactive index (real part)
 real(8):: sci         ! relative water reflactive index (imaginary part)
 integer :: n

 real(8) :: cextp
 real(8) :: cabsp
 real(8) :: php(knang)

 integer :: j
 integer :: idx_spc
 real(sdsu_fpd) :: wgt1, wgt2
 real(sdsu_fps) :: re_in, q_in

!
! initialization
!
 call init_microphysics( cextp_gen,'zero')
 call init_microphysics( cabsp_gen,'zero')
 do n = 1, knang
 call init_microphysics( php_gen(n),'zero')
 enddo


 LUT_GCE: if(lut_visir) then !using LUT approach

!
! cloud species
!
 idx_spc = 1  ; q_in = q%cloud ; re_in = re%cloud

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_gen%cloud = cextp * dble(q_in) 
 cabsp_gen%cloud = cabsp * dble(q_in) 
 php_gen%cloud   = php   * dble(q_in)  

!
! rain species
!
 idx_spc = 2  ; q_in = q%rain ; re_in = re%rain

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_gen%rain = cextp * dble(q_in)
 cabsp_gen%rain = cabsp * dble(q_in)
 php_gen%rain   = php   * dble(q_in)


!
! ice species
!
 idx_spc = 3  ; q_in = q%ice ; re_in = re%ice

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_gen%ice = cextp * dble(q_in)
 cabsp_gen%ice = cabsp * dble(q_in)
 php_gen%ice   = php   * dble(q_in)

!
! snow species
!
 idx_spc = 4  ; q_in = q%snow ; re_in = re%snow

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_gen%snow = cextp * dble(q_in)
 cabsp_gen%snow = cabsp * dble(q_in)
 php_gen%snow   = php   * dble(q_in)

!
! graupel species
!
 idx_spc = 5  ; q_in = q%graupel ; re_in = re%graupel

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_gen%graupel = cextp * dble(q_in)
 cabsp_gen%graupel = cabsp * dble(q_in)
 php_gen%graupel   = php   * dble(q_in)

!
! hail species
!
 idx_spc = 6  ; q_in = q%hail ; re_in = re%hail

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_gen%hail = cextp * dble(q_in)
 cabsp_gen%hail = cabsp * dble(q_in)
 php_gen%hail   = php   * dble(q_in)



 else  !explicit mie calculaiton

    ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  !cloud use RAMS cloud1 PSD
    call mie_rams_visir('qc1',cfmas,pwmas,q%cloud*1.0d-6  ,             & !cloud
                        dble(re%cloud), dble(gnu_rams%cloud1), wlc, rhl  ,             &
                        nang, ang, scr, sci, cextp_gen%cloud, cabsp_gen%cloud, php_gen%cloud )

    call mie_gen_visir('qr ', dble(q%rain),                             & !rain
                       dble(re%rain), dble(mu_gen%rain), dble(rho_gen%rain), wlc,rhl,      &
                       nang, ang, scr, sci, cextp_gen%rain, cabsp_gen%rain, php_gen%rain )

    ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice use RAMS ice1 col PSD
    call mie_rams_visir('qi1', cfmas, pwmas,q%ice*1.0d-6,               & !ice
                        dble(re%ice),dble(gnu_rams%ice1), wlc, rhl,                    &
                        nang, ang, scr, sci, cextp_gen%ice, cabsp_gen%ice, php_gen%ice  )

    call mie_gen_visir('qs ', dble(q%snow),                             & !snow
                       dble(re%snow),dble(mu_gen%snow), dble(rho_gen%snow), wlc,rhl,       &
                       nang, ang, scr, sci, cextp_gen%snow, cabsp_gen%snow, php_gen%snow )

    call mie_gen_visir('qg ', dble(q%graupel),                          & !graupel
                       dble(re%graupel),dble(mu_gen%graupel), dble(rho_gen%graupel), wlc,rhl, &
                       nang, ang, scr, sci, cextp_gen%graupel, cabsp_gen%graupel, php_gen%graupel )

    call mie_gen_visir('qh ', dble(q%hail),                             & !hail
                        dble(re%hail),dble(mu_gen%hail), dble(rho_gen%hail), wlc,rhl,       &
                        nang, ang, scr, sci, cextp_gen%hail, cabsp_gen%hail, php_gen%hail  )

 endif LUT_GCE


 return
 end subroutine opt_gen_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_rams_visir(iw, nt, ns, wlc, rhl, q, re, cextp_rams, cabsp_rams, php_rams)
 implicit none

 integer,intent(in) :: iw  !index of wavelength loop
 integer,intent(in) :: nt,ns !
 real(8),intent(in) :: wlc ! wavelength [cm]
 real(8),intent(in) :: rhl ! relativev humidity  (0~1)
 type ( particle_rams ),intent(in) :: & ! particle_rams is defined in module_simulater
     q,  & ! mixing ratio [g/m3] 
     re    ! effective radius  [micron]

 type ( particle_rams_r8 ),intent(out) :: & !particle_rams is defined in module_simulater
      cextp_rams,  & ! !extinction cross section [cm2/cm3] or [cm-1]
      cabsp_rams,  & ! absorption cross section [cm2/cm3] or [cm-1]
      php_rams(knang)! volume scattering phase function [cm-1  str-1]


! local parameters
 integer :: ispc       ! RAMS spicies index
 integer :: nang       ! number of scaattering angles
 real(8) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(8) :: ang(knang) ! scattering angles [degree]
 real(8):: scr         ! relative water reflactive index (real part)
 real(8):: sci         ! relative water reflactive index (imaginary part)
 integer :: n

 real(8) :: cextp
 real(8) :: cabsp
 real(8) :: php(knang)

 integer :: j
 integer :: idx_spc
 real(sdsu_fpd) :: wgt1, wgt2
 real(sdsu_fps) :: re_in, q_in


!
! initialization
!
 call init_microphysics( cextp_rams,'zero')
 call init_microphysics( cabsp_rams,'zero')
 do n = 1, knang
    call init_microphysics( php_rams(n),'zero')
 enddo

 LUT_RAMS: if(lut_visir) then !using LUT approach

!
! cloud1 species
!
 idx_spc = 1  ; q_in = q%cloud1 ; re_in = re%cloud1

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_rams%cloud1 = cextp * dble(q_in)
 cabsp_rams%cloud1 = cabsp * dble(q_in)
 php_rams%cloud1   = php   * dble(q_in)


!
! cloud2 species
!
 idx_spc = 16  ; q_in = q%cloud2 ; re_in = re%cloud2
 
 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else                   
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1) 
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_rams%cloud2 = cextp * dble(q_in)
 cabsp_rams%cloud2 = cabsp * dble(q_in)
 php_rams%cloud2   = php   * dble(q_in)


!
! rain species
!
 idx_spc = 2  ; q_in = q%rain ; re_in = re%rain
 
 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_rams%rain = cextp * dble(q_in)
 cabsp_rams%rain = cabsp * dble(q_in)
 php_rams%rain   = php   * dble(q_in)


!
! ice 1 habit
!
 ispc = rams_jhabtab(nt,ns,1) !ice1 habit
 idx_spc = ispc  ; q_in = q%ice1 ; re_in = re%ice1

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_rams%ice1 = cextp * dble(q_in)
 cabsp_rams%ice1 = cabsp * dble(q_in)
 php_rams%ice1   = php   * dble(q_in)
 

!
! ice 2 habit
!
 ispc = rams_jhabtab(nt,ns,2) !ice2 habit
 idx_spc = ispc  ; q_in = q%ice2 ; re_in = re%ice2

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_rams%ice2 = cextp * dble(q_in)
 cabsp_rams%ice2 = cabsp * dble(q_in)
 php_rams%ice2   = php   * dble(q_in)


!
! snow habit
!
 idx_spc = 14  ; q_in = q%snow ; re_in = re%snow

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_rams%snow = cextp * dble(q_in)
 cabsp_rams%snow = cabsp * dble(q_in)
 php_rams%snow   = php   * dble(q_in)

!
! graupel habit
!
 idx_spc = 15  ; q_in = q%graupel ; re_in = re%graupel

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_rams%graupel = cextp * dble(q_in)
 cabsp_rams%graupel = cabsp * dble(q_in)
 php_rams%graupel   = php   * dble(q_in)

!
! hail habit
!
 idx_spc = 16  ; q_in = q%hail ; re_in = re%hail

 if(q_in <= q_min_condensate .or. re_in == 0. ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(re_in < pts_re_visir(1)             ) re_in = pts_re_visir(1)
    if(re_in > pts_re_visir(mxpts_re_visir)) re_in = pts_re_visir(mxpts_re_visir)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_re_visir-1
      if( re_in >= pts_re_visir(j) .and. re_in <= pts_re_visir(j+1) ) then
         wgt2 = dble( ( re_in - pts_re_visir(j) )    / ( pts_re_visir(j+1) - pts_re_visir(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_lut_array(j,idx_spc,iw) + wgt2*cextp_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_lut_array(j,idx_spc,iw) + wgt2*cabsp_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_rams%hail = cextp * dble(q_in)
 cabsp_rams%hail = cabsp * dble(q_in)
 php_rams%hail   = php   * dble(q_in)


 else  !explicit mie calculaiton

    ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud1
    call mie_rams_visir('qc1',cfmas,pwmas,q%cloud1*1.0d-6,&
                     dble(re%cloud1), dble(gnu_rams%cloud1), wlc, rhl,&
                     nang, ang, scr, sci, cextp_rams%cloud1, cabsp_rams%cloud1, php_rams%cloud1 )

    ispc = 16; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud2
    call mie_rams_visir('qc2',cfmas,pwmas,q%cloud2*1.0d-6,&
                     dble(re%cloud2), dble(gnu_rams%cloud2), wlc, rhl,&
                     nang, ang, scr, sci, cextp_rams%cloud2, cabsp_rams%cloud2, php_rams%cloud2 )

    ispc = 2 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! rain
    call mie_rams_visir('qr ', cfmas, pwmas,q%rain*1.0d-6,&
                     dble(re%rain),dble(gnu_rams%rain), wlc, rhl,     &
                     nang, ang, scr, sci, cextp_rams%rain, cabsp_rams%rain, php_rams%rain )

    ispc = rams_jhabtab(nt,ns,1) ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1
    call mie_rams_visir('qi1', cfmas, pwmas,q%ice1*1.0d-6,&
                     dble(re%ice1),dble(gnu_rams%ice1), wlc, rhl,     &
                     nang, ang, scr, sci, cextp_rams%ice1, cabsp_rams%ice1, php_rams%ice1 )

    ispc = rams_jhabtab(nt,ns,2) ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2
    call mie_rams_visir('qi2', cfmas, pwmas,q%ice2*1.0d-6,&
                     dble(re%ice2),dble(gnu_rams%ice2), wlc, rhl,      &
                     nang, ang, scr, sci, cextp_rams%ice2, cabsp_rams%ice2, php_rams%ice2 )


    ispc = 5 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! snow
    call mie_rams_visir('qs ', cfmas, pwmas,q%snow*1.0d-6,&
                     dble(re%snow),dble(gnu_rams%snow), wlc, rhl,     &
                     nang, ang, scr, sci, cextp_rams%snow, cabsp_rams%snow, php_rams%snow )

    ispc = 6 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! graupel
    call mie_rams_visir('qg ', cfmas, pwmas,q%graupel*1.0d-6,&
                     dble(re%graupel),dble(gnu_rams%graupel), wlc, rhl,  &
                     nang, ang, scr, sci, cextp_rams%graupel, cabsp_rams%graupel, php_rams%graupel )

    ispc = 7 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! hail
    call mie_rams_visir('qh ', cfmas, pwmas,q%hail*1.0d-6,&
                     dble(re%hail),dble(gnu_rams%hail), wlc, rhl,     &
                     nang, ang, scr, sci, cextp_rams%hail, cabsp_rams%hail, php_rams%hail )

 endif LUT_RAMS



 return
 end subroutine opt_rams_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_sbm_visir(iw, wlc, rhl, q, n, cextp_sbm, cabsp_sbm, php_sbm)
 implicit none

 integer,intent(in) :: iw  !index of wavelength loop
 real(8),intent(in) :: wlc ! wavelength [cm]
 real(8),intent(in) :: rhl ! relativev humidity  (0~1)
 type ( particle_sbm ),intent(in) :: & ! particle_sbm is defined in module_simulater
     q,  &        ! mixing ratio [g/m3] 
     n(1:nbin)    ! particle # concentration [1/m4]

 type ( particle_sbm_r8 ),intent(out) :: & !particle_sbm is defined in module_simulater
      cextp_sbm,  & ! !extinction cross section [cm2/cm3] or [cm-1]
      cabsp_sbm,  & ! absorption cross section [cm2/cm3] or [cm-1]
      php_sbm(knang)! volume scattering phase function [cm-1  str-1]


! local parameters
 integer :: ispc       ! RAMS spicies index
 integer :: nang       ! number of scaattering angles
 real(8) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(8) :: ang(knang) ! scattering angles [degree]
 real(8):: scr         ! relative water reflactive index (real part)
 real(8):: sci         ! relative water reflactive index (imaginary part)


! initialize
 cextp_sbm%liq    =0.d0  ;  cabsp_sbm%liq    =0.d0  ;  php_sbm%liq    =0.d0
 cextp_sbm%ice_col=0.d0  ;  cabsp_sbm%ice_col=0.d0  ;  php_sbm%ice_col=0.d0
 cextp_sbm%ice_pla=0.d0  ;  cabsp_sbm%ice_pla=0.d0  ;  php_sbm%ice_pla=0.d0
 cextp_sbm%ice_den=0.d0  ;  cabsp_sbm%ice_den=0.d0  ;  php_sbm%ice_den=0.d0
 cextp_sbm%snow   =0.d0  ;  cabsp_sbm%snow   =0.d0  ;  php_sbm%snow   =0.d0
 cextp_sbm%graupel=0.d0  ;  cabsp_sbm%graupel=0.d0  ;  php_sbm%graupel=0.d0
 cextp_sbm%hail   =0.d0  ;  cabsp_sbm%hail   =0.d0  ;  php_sbm%hail   =0.d0

 call mie_sbm_visir('liq    ',dble(rad_sbm%liq),dble(drad_sbm%liq),&
                     dble(den_sbm%liq), q%liq*1.0d-6, &
                     n(1:nbin)%liq*1.0d-8, wlc,rhl, &
                     nang, ang, scr, sci, cextp_sbm%liq, cabsp_sbm%liq, php_sbm%liq )

 call mie_sbm_visir('ice_col',dble(rad_sbm%ice_col),dble(drad_sbm%ice_col), &
                     dble(den_sbm%ice_col), q%ice_col*1.0d-6,&
                     n(1:nbin)%ice_col*1.0d-8, wlc,rhl, &
                     nang, ang, scr, sci, cextp_sbm%ice_col, cabsp_sbm%ice_col, php_sbm%ice_col )

 call mie_sbm_visir('ice_pla',dble(rad_sbm%ice_pla),dble(drad_sbm%ice_pla),&
                     dble(den_sbm%ice_pla), q%ice_pla*1.0d-6,&
                     n(1:nbin)%ice_pla*1.0d-8, wlc,rhl, &
                     nang, ang, scr, sci, cextp_sbm%ice_pla, cabsp_sbm%ice_pla, php_sbm%ice_pla )

 call mie_sbm_visir('ice_den',dble(rad_sbm%ice_den),dble(drad_sbm%ice_den),&
                     dble(den_sbm%ice_den), q%ice_den*1.0d-6,&
                     n(1:nbin)%ice_den*1.0d-8, wlc,rhl, &
                     nang, ang, scr, sci, cextp_sbm%ice_den, cabsp_sbm%ice_den, php_sbm%ice_den )

 call mie_sbm_visir('snow   ',dble(rad_sbm%snow),dble(drad_sbm%snow),&
                    dble(den_sbm%snow), q%snow*1.0d-6,&
                     n(1:nbin)%snow*1.0d-8, wlc,rhl, &
                     nang, ang, scr, sci, cextp_sbm%snow, cabsp_sbm%snow, php_sbm%snow )

 call mie_sbm_visir('graupel',dble(rad_sbm%graupel),dble(drad_sbm%graupel),&
                     dble(den_sbm%graupel), q%graupel*1.0d-6,&
                     n(1:nbin)%graupel*1.0d-8, wlc,rhl, &
                     nang, ang, scr, sci, cextp_sbm%graupel, cabsp_sbm%graupel, php_sbm%graupel )

 call mie_sbm_visir('hail   ',dble(rad_sbm%hail),dble(drad_sbm%hail),&
                     dble(den_sbm%hail), q%hail*1.0d-6,&
                     n(1:nbin)%hail*1.0d-8, wlc,rhl, &
                     nang, ang, scr, sci, cextp_sbm%hail, cabsp_sbm%hail, php_sbm%hail )

 return
 end subroutine opt_sbm_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_gocart_visir(iw, wlc, rhl, q, cextp_gocart, cabsp_gocart, php_gocart)
 implicit none

 integer,intent(in) :: iw  !index of wavelength loop
 real(8),intent(in) :: wlc ! wavelength [cm]
 real(8),intent(in) :: rhl ! relativev humidity  (0~1)

 type ( particle_gocart ),intent(in) :: & ! particle_gocart is defined in module_simulater
     q     ! mixing ratio [g/m3] 


 type ( particle_gocart_r8 ),intent(out) :: & !particle_gocart is defined in module_simulater
      cextp_gocart,  & ! !extinction cross section [cm2/cm3] or [cm-1]
      cabsp_gocart,  & ! absorption cross section [cm2/cm3] or [cm-1]
      php_gocart(knang)! volume scattering phase function [cm-1  str-1]


! local parameters
 integer :: ispc       ! RAMS spicies index
 integer :: nang       ! number of scaattering angles
 real(8) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(8) :: ang(knang) ! scattering angles [degree]
 real(8):: scr         ! relative water reflactive index (real part)
 real(8):: sci         ! relative water reflactive index (imaginary part)


 real(8) :: cextp
 real(8) :: cabsp
 real(8) :: php(knang)

 integer :: j,n
 integer :: idx_spc
 real(sdsu_fpd) :: wgt1, wgt2
 real(sdsu_fps) :: rh_in, q_in

!
! initialization
!
 call init_aerosol( cextp_gocart,'zero')
 call init_aerosol( cabsp_gocart,'zero')
 do n = 1, knang
 call init_aerosol( php_gocart(n),'zero')
 enddo
 

 LUT_GOCART: if(lut_visir) then !using LUT approach

!
! sulfate
!
 idx_spc = 1  ; q_in = q%so4 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
    if(rh_in < pts_rh(1)             ) rh_in = pts_rh(1)
    if(rh_in > pts_rh(mxpts_rh)) rh_in = pts_rh(mxpts_rh)

    !   
    ! Linear interpolation
    !
    do j = 1, mxpts_rh-1
      if( rh_in >= pts_rh(j) .and. rh_in <= pts_rh(j+1) ) then
         wgt2 = dble( ( rh_in - pts_rh(j) )    / ( pts_rh(j+1) - pts_rh(j) ) )
         wgt1 = dble( 1.0 - wgt2 )

         cextp = wgt1*cextp_aero_lut_array(j,idx_spc,iw) + wgt2*cextp_aero_lut_array(j+1,idx_spc,iw)
         cabsp = wgt1*cabsp_aero_lut_array(j,idx_spc,iw) + wgt2*cabsp_aero_lut_array(j+1,idx_spc,iw)
         php(1:knang) = wgt1*php_aero_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_aero_lut_array(1:knang,j+1,idx_spc,iw)
         exit
      endif
    enddo
 endif

! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_gocart%so4 = cextp * dble(q_in)
 cabsp_gocart%so4 = cabsp * dble(q_in)
 php_gocart%so4   = php   * dble(q_in)

!
! black carbon 
!
 idx_spc = 2  ; q_in = q%blc ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = wgt1*cextp_aero_lut_array(j,idx_spc,iw) + wgt2*cextp_aero_lut_array(j+1,idx_spc,iw)
     cabsp = wgt1*cabsp_aero_lut_array(j,idx_spc,iw) + wgt2*cabsp_aero_lut_array(j+1,idx_spc,iw)
     php(1:knang) = wgt1*php_aero_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_aero_lut_array(1:knang,j+1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_gocart%blc = cextp * dble(q_in)
 cabsp_gocart%blc = cabsp * dble(q_in)
 php_gocart%blc   = php   * dble(q_in)

!
! non hygroscopic OC 
!
 idx_spc = 3  ; q_in = q%ocn ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw) 
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw) 
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])
 cextp_gocart%ocn = cextp * dble(q_in)
 cabsp_gocart%ocn = cabsp * dble(q_in)
 php_gocart%ocn   = php   * dble(q_in)

!
! hygroscopic OC 
!
 idx_spc = 4  ; q_in = q%och ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = wgt1*cextp_aero_lut_array(j,idx_spc,iw) + wgt2*cextp_aero_lut_array(j+1,idx_spc,iw)
     cabsp = wgt1*cabsp_aero_lut_array(j,idx_spc,iw) + wgt2*cabsp_aero_lut_array(j+1,idx_spc,iw)
     php(1:knang) = wgt1*php_aero_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_aero_lut_array(1:knang,j+1,idx_spc,iw)
 endif                          
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%och = cextp * dble(q_in)
 cabsp_gocart%och = cabsp * dble(q_in)
 php_gocart%och   = php   * dble(q_in)

!
! fine-mode sea salt
!
 idx_spc = 5  ; q_in = q%ssa ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = wgt1*cextp_aero_lut_array(j,idx_spc,iw) + wgt2*cextp_aero_lut_array(j+1,idx_spc,iw)
     cabsp = wgt1*cabsp_aero_lut_array(j,idx_spc,iw) + wgt2*cabsp_aero_lut_array(j+1,idx_spc,iw)
     php(1:knang) = wgt1*php_aero_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_aero_lut_array(1:knang,j+1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%ssa = cextp * dble(q_in)
 cabsp_gocart%ssa = cabsp * dble(q_in)
 php_gocart%ssa   = php   * dble(q_in)

!
! coarse-mode sea salt
!
 idx_spc = 6  ; q_in = q%ssc ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = wgt1*cextp_aero_lut_array(j,idx_spc,iw) + wgt2*cextp_aero_lut_array(j+1,idx_spc,iw)
     cabsp = wgt1*cabsp_aero_lut_array(j,idx_spc,iw) + wgt2*cabsp_aero_lut_array(j+1,idx_spc,iw)
     php(1:knang) = wgt1*php_aero_lut_array(1:knang,j,idx_spc,iw) + wgt2*php_aero_lut_array(1:knang,j+1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%ssc = cextp * dble(q_in)
 cabsp_gocart%ssc = cabsp * dble(q_in)
 php_gocart%ssc   = php   * dble(q_in)


!
! dust1
!
 idx_spc = 7  ; q_in = q%du1 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw) 
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%du1 = cextp * dble(q_in)
 cabsp_gocart%du1 = cabsp * dble(q_in)
 php_gocart%du1   = php   * dble(q_in)

!
! dust2
!
 idx_spc = 8  ; q_in = q%du2 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw) 
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw) 
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%du2 = cextp * dble(q_in)
 cabsp_gocart%du2 = cabsp * dble(q_in)
 php_gocart%du2   = php   * dble(q_in)

!
! dust3
!
 idx_spc = 9  ; q_in = q%du3 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw) 
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw) 
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%du3 = cextp * dble(q_in)
 cabsp_gocart%du3 = cabsp * dble(q_in)
 php_gocart%du3   = php   * dble(q_in)

!
! dust4
!
 idx_spc = 10  ; q_in = q%du4 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw) 
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw) 
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%du4 = cextp * dble(q_in)
 cabsp_gocart%du4 = cabsp * dble(q_in)
 php_gocart%du4   = php   * dble(q_in)

!
! dust5
!
 idx_spc = 11  ; q_in = q%du5 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw) 
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw) 
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%du5 = cextp * dble(q_in)
 cabsp_gocart%du5 = cabsp * dble(q_in)
 php_gocart%du5   = php   * dble(q_in)

!
! dust6
!
 idx_spc = 12  ; q_in = q%du6 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw) 
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw) 
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%du6 = cextp * dble(q_in)
 cabsp_gocart%du6 = cabsp * dble(q_in)
 php_gocart%du6   = php   * dble(q_in)

!
! dust7
!
 idx_spc = 13  ; q_in = q%du7 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw) 
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw) 
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%du7 = cextp * dble(q_in)
 cabsp_gocart%du7 = cabsp * dble(q_in)
 php_gocart%du7   = php   * dble(q_in)

!
! dust8
!
 idx_spc = 14  ; q_in = q%du8 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw)
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_gocart%du8 = cextp * dble(q_in)
 cabsp_gocart%du8 = cabsp * dble(q_in)
 php_gocart%du8   = php   * dble(q_in)


 else

    call mie_gocart('so4', q%so4*1.0d-6, wlc,rhl,         &   !sulfur and its precure
                     nang, ang, scr, sci, cextp_gocart%so4 , cabsp_gocart%so4 , php_gocart%so4 )
    call mie_gocart('blc', q%blc*1.0d-6, wlc,rhl,         &   !black carbon (soot)
                     nang, ang, scr, sci, cextp_gocart%blc , cabsp_gocart%blc , php_gocart%blc )
    call mie_gocart('ocn', q%ocn*1.0d-6, wlc,rhl,         &   !non hygroscopic OC
                     nang, ang, scr, sci, cextp_gocart%ocn , cabsp_gocart%ocn , php_gocart%ocn )
    call mie_gocart('och', q%och*1.0d-6, wlc,rhl,         &   !hygroscopic OC
                     nang, ang, scr, sci, cextp_gocart%och , cabsp_gocart%och , php_gocart%och )
    call mie_gocart('ssa', q%ssa*1.0d-6, wlc,rhl,         &   !sea-salt accumulation mode
                     nang, ang, scr, sci, cextp_gocart%ssa , cabsp_gocart%ssa , php_gocart%ssa )
    call mie_gocart('ssc', q%ssc*1.0d-6, wlc,rhl,         &   !sea-salt coarse mode
                     nang, ang, scr, sci, cextp_gocart%ssc , cabsp_gocart%ssc , php_gocart%ssc )
    call mie_gocart('du1', q%du1*1.0d-6, wlc,rhl,         &   !dust mode 1
                     nang, ang, scr, sci, cextp_gocart%du1 , cabsp_gocart%du1 , php_gocart%du1 )
    call mie_gocart('du2', q%du2*1.0d-6, wlc,rhl,         &   !dust mode 2
                     nang, ang, scr, sci, cextp_gocart%du2 , cabsp_gocart%du2 , php_gocart%du2 )
    call mie_gocart('du3', q%du3*1.0d-6, wlc,rhl,         &   !dust mode 3
                     nang, ang, scr, sci, cextp_gocart%du3 , cabsp_gocart%du3 , php_gocart%du3 )
    call mie_gocart('du4', q%du4*1.0d-6, wlc,rhl,         &   !dust mode 4
                     nang, ang, scr, sci, cextp_gocart%du4 , cabsp_gocart%du4 , php_gocart%du4 )
    call mie_gocart('du5', q%du5*1.0d-6, wlc,rhl,         &   !dust mode 5
                     nang, ang, scr, sci, cextp_gocart%du5 , cabsp_gocart%du5 , php_gocart%du5 )
    call mie_gocart('du6', q%du6*1.0d-6, wlc,rhl,         &   !dust mode 6
                     nang, ang, scr, sci, cextp_gocart%du6 , cabsp_gocart%du6 , php_gocart%du6 )
    call mie_gocart('du7', q%du7*1.0d-6, wlc,rhl,         &   !dust mode 7
                     nang, ang, scr, sci, cextp_gocart%du7 , cabsp_gocart%du7 , php_gocart%du7 )
    call mie_gocart('du8', q%du8*1.0d-6, wlc,rhl,         &   !dust mode 8
                     nang, ang, scr, sci, cextp_gocart%du8 , cabsp_gocart%du8 , php_gocart%du8 )

 endif LUT_GOCART


 return
 end subroutine opt_gocart_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine opt_ucd_visir(iw, wlc, rhl, q, cextp_ucd, cabsp_ucd, php_ucd)
 implicit none

 integer,intent(in) :: iw  !index of wavelength loop
 real(8),intent(in) :: wlc ! wavelength [cm]
 real(8),intent(in) :: rhl ! relativev humidity  (0~1)

 type ( particle_ucd ),intent(in) :: & ! particle_ucd is defined in module_simulater
     q     ! mixing ratio [g/m3] 


 type ( particle_ucd_r8 ),intent(out) :: & !particle_ucd is defined in module_simulater
      cextp_ucd,  & ! !extinction cross section [cm2/cm3] or [cm-1]
      cabsp_ucd,  & ! absorption cross section [cm2/cm3] or [cm-1]
      php_ucd(knang)! volume scattering phase function [cm-1  str-1]


! local parameters
 integer :: ispc       ! spicies index
 integer :: nang       ! number of scaattering angles
 real(8) :: cfmas, pwmas !alpha_mass and beta_mass for RAMS mass-dimater equation
 real(8) :: ang(knang) ! scattering angles [degree]
 real(8):: scr         ! relative water reflactive index (real part)
 real(8):: sci         ! relative water reflactive index (imaginary part)

 real(8) :: cextp
 real(8) :: cabsp
 real(8) :: php(knang)

 integer :: j,n
 integer :: idx_spc
 real(sdsu_fpd) :: wgt1, wgt2
 real(sdsu_fps) :: rh_in, q_in

!
! initialization
!
 call init_aerosol( cextp_ucd,'zero')
 call init_aerosol( cabsp_ucd,'zero')
 do n = 1, knang
 call init_aerosol( php_ucd(n),'zero')
 enddo


 LUT_UCD: if(lut_visir) then !using LUT approach

!
! dust1
!
 idx_spc = 1  ; q_in = q%du1 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw)
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_ucd%du1 = cextp * dble(q_in)
 cabsp_ucd%du1 = cabsp * dble(q_in)
 php_ucd%du1   = php   * dble(q_in)

!
! dust2
!
 idx_spc = 2  ; q_in = q%du2 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw)
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_ucd%du2 = cextp * dble(q_in)
 cabsp_ucd%du2 = cabsp * dble(q_in)
 php_ucd%du2   = php   * dble(q_in)

!
! dust3
!
 idx_spc = 3  ; q_in = q%du3 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw)
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_ucd%du3 = cextp * dble(q_in)
 cabsp_ucd%du3 = cabsp * dble(q_in)
 php_ucd%du3   = php   * dble(q_in)

!
! dust4
!
 idx_spc = 4  ; q_in = q%du4 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw)
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_ucd%du4 = cextp * dble(q_in)
 cabsp_ucd%du4 = cabsp * dble(q_in)
 php_ucd%du4   = php   * dble(q_in)

!
! dust5
!
 idx_spc = 5  ; q_in = q%du5 ; rh_in = rhl
 if(q_in <= q_min_aerosol ) then
     cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 else
     cextp = cextp_aero_lut_array(1,idx_spc,iw)
     cabsp = cabsp_aero_lut_array(1,idx_spc,iw)
     php(1:knang) = php_aero_lut_array(1:knang,1,idx_spc,iw)
 endif
! Re-scale from unit mass mixing ratio (1[g/m3]) to actual mass mixing ratio (q [g/m3])   
 cextp_ucd%du5 = cextp * dble(q_in)
 cabsp_ucd%du5 = cabsp * dble(q_in)
 php_ucd%du5   = php   * dble(q_in)


 else

     call mie_ucd('du1', q%du1*1.0d-6, wlc,rhl,         &   !dust mode 1
                  nang, ang, scr, sci, cextp_ucd%du1 , cabsp_ucd%du1 , php_ucd%du1 )
     call mie_ucd('du2', q%du2*1.0d-6, wlc,rhl,         &   !dust mode 2
                  nang, ang, scr, sci, cextp_ucd%du2 , cabsp_ucd%du2 , php_ucd%du2 )
     call mie_ucd('du3', q%du3*1.0d-6, wlc,rhl,         &   !dust mode 3
                  nang, ang, scr, sci, cextp_ucd%du3 , cabsp_ucd%du3 , php_ucd%du3 )
     call mie_ucd('du4', q%du4*1.0d-6, wlc,rhl,         &   !dust mode 4
                  nang, ang, scr, sci, cextp_ucd%du4 , cabsp_ucd%du4 , php_ucd%du4 )
     call mie_ucd('du5', q%du5*1.0d-6, wlc,rhl,         &   !dust mode 5
                  nang, ang, scr, sci, cextp_ucd%du5 , cabsp_ucd%du5 , php_ucd%du5 )

 endif LUT_UCD


 return
 end subroutine opt_ucd_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mie_gen_visir(spc, q, re, mu, rho, wlc, rh, nang, ang, crw1, ciw1, cextp, cabsp, php  )
 implicit none
 
 character(len=3),intent(in) :: spc  ! particle species character index
 real(8),intent(in):: q          ! mass mixing ratio of particle [g/m3]
 real(8),intent(in):: re         ! effective particle radius [micron]
 real(8),intent(in):: rho        ! particle density [kg/m3]
 real(8),intent(in) :: mu        ! mu of gamma size distribution [-]
 real(8),intent(in):: wlc        ! wavelength [cm]
 real(8),intent(in):: rh         ! relativev humidity  (0~1)
 integer,intent(out):: nang       ! number of scaattering angles
 real(8),intent(out):: ang(knang) ! scattering angles [degree]
 real(8),intent(out):: crw1       ! relative reflactive index (real part)
 real(8),intent(out):: ciw1       ! relative reflactive index (imaginary part)
 real(8),intent(out):: cextp      ! extinction cross section [cm-1]
 real(8),intent(out):: cabsp      ! absorption cross section [cm-1]
 real(8),intent(out):: php(knang) ! volume scattering phase function [cm-1  str-1]

 integer :: ispcvp(3)     ! internal mix of fundamental materials
 real(8) :: rfracp(3)     ! dry component volume fractions
 real(8) :: asphr(3)      ! non-spherical parameters
 real(8) :: rop           ! particle density of dry mixture [g/cm3]
 real(8) :: psd_params(6,33)   ! params to define volume size dist.
 real(8) :: rm       ! mode and effective radius for aerosols [um]

 real(8) :: lambda !lambda of Marshal-Palmer PSD [1/m]
 real(8) :: n0     !intercept of Marshal-Palmer PSD [1/m4] 

 real(8) :: wgt1,wgt2
 integer :: irh ,ispc 
 real(8) :: wfrac ! water fraction
 real(8) :: rmdry ! mode radius for dry particle
 logical :: success
!
! 1moment case
! PSD parameter packet psd_params(6,33)=(i,j)
! 
!   i   1       2       3     4     5    6  
! j -----------------------------------------
! 1 |   r   psd_type   PSD constant parameter
!   |
! 2 | nmode   
!   |
! 3 | rmin 
!   |
! 4 | rmax
!   |
!

!--exec

!initialise
 cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
 psd_params = 0. 

 if( trim(spc) /= 'ini' .and. q*1.0d-6 <= q_min_conden) return


 if( ss_opt_visir == 1 ) then

!  call interp_icescat_visir(q, re, wlc, cextp, cabsp, php, success )

 call ss_visir_gen_micro  (trim(spc), wlc*1.d+4, q, re, mu, rho, &
                            cextp, cabsp, php, success)
  if(success)  return  ! if not sucess go to old routine below....

 endif



!
! compute lambda and n0 for given re and q  
!
 lambda = 3.d0 / 2.d0 / re  * 1.0d+6   ! lambda of Marshal-Palmer PSD [1/m]  (corrected by Hashino)

 n0 = ( (2.0d0 / (re* 1.0d-6))**4  )  * q * 1.0d-3 / (pi * rho)  !intercept of Marshal-Palmer PSD [1/m4]



!
! species-independent constant optical parameter
!
  psd_params(1,2) = 1.d0      !# of modes is always 1 for SDSU

 select case(trim(spc))
 case('qc','qi','ini')  ! non-precip
   psd_params(1,3) =     1.0d-4  ! 1micron        (must start from this for stability purpose)
   psd_params(1,4) = 100.0d-4    ! 100micron  (consistent to module_opt_micro)
 case('qr','qs','qg','qh') ! precip
   psd_params(1,3) =   1.0d-4  ! 1micron        (must start from this for stability purpose)
   psd_params(1,4) = 10000.0d-4    ! 10000micron = 1cm (consistent to module_opt_micro)
 case default
    print*, 'MSG mie_rams_visir: There is no such spc=',spc
   call stop_sdsu(myrank,'')
 end select

  psd_params(2,1) = 4                 ! Marshal-Palmer PSD: n0 * exp(-lambda*D) , D=2*r  !(dN/dr is [1/cm4])
  psd_params(3,1) = n0 * 1.0d-8       ! intercept [1/cm4] <- [1/m4]
  psd_params(4,1) = lambda * 1.0d-2   ! lambda [1/cm] <- [1/m]

  rop = 1.0             ! suedo density

!
! species-dependent constant optical parameter (only 1mode)
!
  spc_select: select case(trim(spc))

  case('qr','ini') !rain
     ispcvp(1)=1     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case('qs') !snow
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case('qg') !graupel
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case('qh') !hail
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case default
      print*, 'MSG mie_gen_visir: There is no such spc=',spc
      call stop_sdsu(myrank,'')  
  end select spc_select

!
! compute unified refractive index & optical properties
!
  call rfidxb_sdsu( wlc,rh,ispcvp,rfracp,asphr,rop,psd_params, &
                   nang,ang,crw1,ciw1,cextp,cabsp,php )

!
! spherical particle assume no polarization --> P11(180deg) = P22(180deg)  
!
 php(nang+1) = php(nang)  !P22 phase function at 180deg [cm-1 str-1]


 return
 end subroutine mie_gen_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mie_gocart(spc, q, wlc, rh, nang, ang, crw1, ciw1, cextp, cabsp, php  )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute mie scattering from GOCART aerosol odel. 
!  Currently it assumes log-normal PSD using the similar re/rm identical to GOCART.
! 
! History:
! 10/2017  Toshi Matsui@NASA GSFC ; Revised log-normal
! 10/2010  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 
 character(len=3),intent(in) :: spc  ! particle species character index
 real(8),intent(in):: q         ! mass mixing ratio of particle [g/cm3]
 real(8),intent(in):: wlc       ! wavenumber [cm]
 real(8),intent(in):: rh        ! relativev humidity  (0~1)
 integer,intent(out):: nang       !! number of scaattering angles
 real(8),intent(out):: ang(knang) !! scattering angles [degree]
 real(8),intent(out):: crw1     !! relative reflactive index (real part)
 real(8),intent(out):: ciw1     !! relative reflactive index (imaginary part)
 real(8),intent(out):: cextp    ! extinction cross section [cm-1]
 real(8),intent(out):: cabsp    ! absorption cross section [cm-1]
 real(8),intent(out):: php(knang) ! volume scattering phase function [cm-1  str-1]

 integer :: ispcvp(3)     ! internal mix of fundamental materials
 real(8) :: rfracp(3)     ! dry component volume fractions
 real(8) :: asphr(3)      ! non-spherical parameters
 real(8) :: rop            ! particle density of dry mixture [g/cm3]
 real(8) :: psd_params(6,33)   ! params to define volume size dist.
 real(8) :: rm,re       !mode and effective radius for aerosols [um]
 real(8) :: sigma_aero

 real(8) :: wgt1,wgt2
 integer :: irh ,ispc 
 real(8) :: wfrac !water fraction
 real(8) :: rmdry !mode radius for dry particle


 real(8),parameter ::  dust_nonspherical_p22 = 0.8   !empirical P11 to P22 ratio 
 real(8),parameter ::  dust_nospherical_ratio = 0.5

 real(8),parameter ::  dust_depol_ratio = 0.8
 logical,parameter :: new_lognormal_psd = .true.

 real(8) :: rm_cm  ! mode radius [cm]
 real(8) :: sigma ! log of standard deviation
 real(8) :: Ntot  !total number concentration

!
! 1moment case
! PSD parameter packet psd_params(6,33)=(i,j)
! 
!   i   1       2       3     4     5    6  
! j -----------------------------------------
! 1 |   r   psd_type   PSD constant parameter
!   |
! 2 | nmode   
!   |
! 3 | rmin 
!   |
! 4 | rmax
!   |
!

!--exec

!
!initialise
!
  cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
  psd_params = 0. 

  if( trim(spc) /= 'ini' .and. q <= q_min_aero) return
  

!
! find RH index
!
  if( rh > 0.99d0 ) then
      irh = 98 ; wgt1 = 0.d0 ; wgt2 = 1.0d0
  elseif( rh < 0.d0 ) then
      irh =  0 ; wgt1 = 1.d0 ; wgt2 = 0.0d0
  else !
      irh = INT(rh*100.d0) 
      wgt2 = rh*100.d0 - REAL(irh)
      wgt1 = 1.0d0 - wgt2
  endif
 if(irh > 98) irh = 98  !make sure to avoid bounds problem. 

!
! species-independent constant optical parameter
!
 psd_params(1,2) = 1.d0    !# of modes is always 1 for SDSU
 psd_params(1,3) = 0.01d-4 ! minimum particle radius [cm]  0.01 micron
 psd_params(1,4) = 30.0d-4 ! maximum particle radius [cm] 30.0  micron

!
! species-dependent constant optical parameter (only 1mode)
!
  spc_select: select case(trim(spc))

! fundamental material used for 3 component mixture
!                     1 water
!                     2 ice
!                     3 dust-like
!                     4 water-soluble
!                     5 sea salt
!                     6 soot
!                     7 volcanic-ash
!                     8 75%H2SO4
!                     9 Yellow sand

  case('so4','ini') !sulfate and its precursor
     ispcvp(1)=8     ; ispcvp(2)=1    ; ispcvp(3)=0     ! 3 component mixture
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere
     rop=1.0                      ! density [g/cm]
     ispc = 1
     rmdry = rm_gocart(ispc,0)   !dry mode radius [um]
     rm=wgt1*rm_gocart(ispc,int(irh))+wgt2*rm_gocart(ispc,int(irh+1)) !mode radius [um]
     re=wgt1*re_gocart(ispc,int(irh))+wgt2*re_gocart(ispc,int(irh+1)) !effective radius [um]
     wfrac = 1.d0 - (rmdry/rm)**3  !water fraction
     rfracp(1)=1.d0-wfrac ; rfracp(2)=wfrac ; rfracp(3)=0.00  ! dry component volume fractions (sphere)

  case('blc') !black carbon (soot)
     ispcvp(1)=6     ; ispcvp(2)=1    ; ispcvp(3)=0     ! 3 component mixture
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere
     rop=1.0                      ! density [g/cm]
     ispc = 3
     rmdry = rm_gocart(ispc,0)   !dry mode radius [um]
     rm=wgt1*rm_gocart(ispc,int(irh))+wgt2*rm_gocart(ispc,int(irh+1)) !mode radius [um]
     re=wgt1*re_gocart(ispc,int(irh))+wgt2*re_gocart(ispc,int(irh+1)) !effective radius [um]
     wfrac = 1.d0 - (rmdry/rm)**3  !water fraction
     rfracp(1)=1.d0-wfrac ; rfracp(2)=wfrac ; rfracp(3)=0.00  ! dry component volume fractions (sphere)

  case('ocn','och') !non-hygroscopic or hygroscopic organic carbon
     ispcvp(1)=4     ; ispcvp(2)=1    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere
     rop=1.0                      ! density [g/cm]
     if( spc == 'ocn') then
        ispc = 2 
        rmdry = rm_gocart(ispc,0)  !dry mode radius
        rm=     rm_gocart(ispc,0)  !mode radius [um]
        re=     re_gocart(ispc,0)  !effective radius [um]
     elseif( spc == 'och' ) then
        ispc = 2 
        rmdry = rm_gocart(ispc,0)  !dry mode radius
        rm=wgt1*rm_gocart(ispc,int(irh))+wgt2*rm_gocart(ispc,int(irh+1)) !mode radius [um]
        re=wgt1*re_gocart(ispc,int(irh))+wgt2*re_gocart(ispc,int(irh+1)) !effective radius [um]
     endif
     wfrac = 1.d0 - (rmdry/rm)**3  !water fraction
     rfracp(1)=1.d0-wfrac ; rfracp(2)=wfrac ; rfracp(3)=0.00  ! dry component volume fractions (sphere)

  case('ssa','ssc') !non-hygroscopic or hygroscopic organic carbon
     ispcvp(1)=5     ; ispcvp(2)=1    ; ispcvp(3)=0     ! 3 component mixture
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere
     rop=1.0                      ! density [g/cm]
     if( spc == 'ssa') then
        ispc = 4
        rmdry = rm_gocart(ispc,0)  !dry mode radius
        rm=wgt1*rm_gocart(ispc,int(irh))+wgt2*rm_gocart(ispc,int(irh+1)) !mode radius [um]
        re=wgt1*re_gocart(ispc,int(irh))+wgt2*re_gocart(ispc,int(irh+1)) !effective radius [um]
     elseif( spc == 'ssc' ) then
        ispc = 5
        rmdry = rm_gocart(ispc,0)
        rm=wgt1*rm_gocart(ispc,int(irh))+wgt2*rm_gocart(ispc,int(irh+1)) !mode radius [um]
        re=wgt1*re_gocart(ispc,int(irh))+wgt2*re_gocart(ispc,int(irh+1)) !effective radius [um]
     endif
     wfrac = 1.d0 - (rmdry/rm)**3   !water fraction
     rfracp(1)=1.d0-wfrac ; rfracp(2)=wfrac ; rfracp(3)=0.00  ! dry component volume fractions (sphere)

  case('du1','du2','du3','du4','du5','du6','du7','du8') !dust mode 1 ~ 8
     ispcvp(1)=3     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere
!     asphr (1)=7.0 ; asphr (2)=10.0 ; asphr (3)=1.1  ![Nakajima et al. ]
     rop= 2.56                     ! density [g/cm]

     psd_params(2,1)=2        ! log-normal: vj = c * exp((ln(r/rm)/ln(s))**2 / 2)

     if( spc == 'du1') then !these values are consistent to GOCART spectrum AOD calculations
       rm=rm_gocart(6,1) ; re=re_gocart(6,1)  !mode & effective radius [um]
     elseif( spc == 'du2') then
       rm=rm_gocart(6,2) ; re=re_gocart(6,2)  !mode & effective radius [um]
     elseif( spc == 'du3') then
       rm=rm_gocart(6,3) ; re=re_gocart(6,3)  !mode & effective radius [um]
     elseif( spc == 'du4') then
       rm=rm_gocart(6,4) ; re=re_gocart(6,4)  !mode & effective radius [um]
     elseif( spc == 'du5') then
       rm=rm_gocart(6,5) ; re=re_gocart(6,5)  !mode & effective radius [um]
     elseif( spc == 'du6') then
       rm=rm_gocart(6,6) ; re=re_gocart(6,6)  !mode & effective radius [um]
     elseif( spc == 'du7') then
       rm=rm_gocart(6,7) ; re=re_gocart(6,7)  !mode & effective radius [um]
     elseif( spc == 'du8') then
       rm=rm_gocart(6,8) ; re=re_gocart(6,8)  !mode & effective radius [um]
     endif

  case default
      print*, 'MSG mie_gocart: There is no such spc=',spc
      call stop_sdsu(myrank,'')  
  end select spc_select

 if( new_lognormal_psd ) then

    psd_params(2,1)=8        ! new log normal
    sigma = 2. * log( re/rm ) / 5.  ! sigma = ln2(S)  (0.48 ~ 1.5)
    rm_cm = rm * 1.0d-4  !  [cm]
    !total number concentration 
    Ntot = 3. * q / ( 4.*pi*rop*(rm_cm**3) * exp( 9.*(sigma**2)/2. ) )
    psd_params(3,1) = Ntot          ! total number concentration [1/cm3]
    psd_params(4,1) = sigma         ! s: log standard deviation [-]
    psd_params(5,1) = log( rm_cm )  !mu: mean --> log(rm) [cm] <- [um] 

 else

     psd_params(2,1)=2        ! log-normal: vj = c * exp((ln(r/rm)/ln(s))**2 / 2)
     sigma_aero = sqrt( 0.4d0*log(re/rm) )   !~0.45
     psd_params(3,1)=q / (sqrt(2.d0*pi) * sigma_aero)   ! c:
     psd_params(4,1)=sigma_aero  ! s: log standard deviation
     psd_params(5,1)=rm*1.0d-4   ! rm: mode radius [cm] <- [um] 

 endif

!
! compute unified refractive index & optical properties
!
  call rfidxb_sdsu( wlc,rh,ispcvp,rfracp,asphr,rop,psd_params, &
                   nang,ang,crw1,ciw1,cextp,cabsp,php )



 select case(trim(spc) )
 case('du1','du2','du3','du4','du5','du6','du7','du8') !dust mode 1 ~ 8

   !
   ! assume non-spherical shape impact on P22 backscatter 
   !
   php(nang+1) = php(nang) * dust_nonspherical_p22  ! P22 phase function at 180deg [cm-1 str-1]

   !
   ! dust non-spherical backscatter correction 
   ! This suppress backscatter peak from Mie scattering about 50%. 
   !
   php(nang)   = php(nang)   * dust_nospherical_ratio
   php(nang+1) = php(nang+1) * dust_nospherical_ratio

 case default
   !
   ! spherical particle assume no polarization --> P11(180deg) = P22(180deg)  
   !
   php(nang+1) = php(nang)                     ! P22 phase function at 180deg [cm-1 str-1]

 end select


 return
 end subroutine mie_gocart

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mie_ucd(spc, q, wlc, rh, nang, ang, crw1, ciw1, cextp, cabsp, php  )
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
!  Compute mie scattering from UCD dust model. 
!  Currently it assumes log-normal PSD using the similar re/rm relationship to GOCART,
!  although it assume fixed size (bin) approach. 
! 
! History:
! 10/2017  Toshi Matsui@NASA GSFC ; Initial 
!           
! References: 
!----------------------------------------------------------------------------------------------------
 
 character(len=3),intent(in) :: spc  ! particle species character index
 real(8),intent(in):: q         ! mass mixing ratio of particle [g/cm3]
 real(8),intent(in):: wlc       ! wavenumber [cm]
 real(8),intent(in):: rh        ! relativev humidity  (0~1)
 integer,intent(out):: nang       !! number of scaattering angles
 real(8),intent(out):: ang(knang) !! scattering angles [degree]
 real(8),intent(out):: crw1     !! relative reflactive index (real part)
 real(8),intent(out):: ciw1     !! relative reflactive index (imaginary part)
 real(8),intent(out):: cextp    ! extinction cross section [cm-1]
 real(8),intent(out):: cabsp    ! absorption cross section [cm-1]
 real(8),intent(out):: php(knang) ! volume scattering phase function [cm-1  str-1]

 integer :: ispcvp(3)     ! internal mix of fundamental materials
 real(8) :: rfracp(3)     ! dry component volume fractions
 real(8) :: asphr(3)      ! non-spherical parameters
 real(8) :: rop            ! particle density of dry mixture [g/cm3]
 real(8) :: psd_params(6,33)   ! params to define volume size dist.
 real(8) :: rm, re, rbin, dr    ! radius and dr for aerosols [um]
 real(8) :: rm_cm  ! mode radius [cm]
 real(8) :: sigma_aero
 real(8) :: sigma ! log of standard deviation
 real(8) :: Ntot  !total number concentration

 real(8) :: wgt1,wgt2
 integer :: irh ,ispc 
 real(8) :: wfrac !water fraction
 real(8) :: rmdry !mode radius for dry particle

 real(8),parameter ::  dust_nonspherical_p22 = 0.8   !empirical P11 to P22 ratio 
 real(8),parameter ::  dust_nospherical_ratio = 0.5

 real(8),parameter :: stdev_const = 0.7  !similar to GOCART assumption 
 logical,parameter :: lognormal_psd_assum = .true.
 logical,parameter :: new_lognormal_psd = .true.

!
!initialise
!
  cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
  psd_params = 0. 

  if( trim(spc) /= 'ini' .and. q <= q_min_aero) return
!  if( trim(spc) /= 'ini' .and. q <= 0. ) return !q_min_aero) return
  

!
! find RH index
!
  if( rh > 0.99d0 ) then
      irh = 98 ; wgt1 = 0.d0 ; wgt2 = 1.0d0
  elseif( rh < 0.d0 ) then
      irh =  0 ; wgt1 = 1.d0 ; wgt2 = 0.0d0
  else !
      irh = INT(rh*100.d0) 
      wgt2 = rh*100.d0 - REAL(irh)
      wgt1 = 1.0d0 - wgt2
  endif
 if(irh > 98) irh = 98  !make sure to avoid bounds problem. 

!
! species-independent constant optical parameter
!
 psd_params(1,2) = 1.d0    !# of modes is always 1 for SDSU
 psd_params(1,3) = 0.01d-4 ! minimum particle radius [cm]  0.01 micron
 psd_params(1,4) = 30.0d-4 ! maximum particle radius [cm] 30.0  micron

!
! species-dependent constant optical parameter (only 1mode)
!

! fundamental material used for 3 component mixture
!                     1 water
!                     2 ice
!                     3 dust-like
!                     4 water-soluble
!                     5 sea salt
!                     6 soot
!                     7 volcanic-ash
!                     8 75%H2SO4
!                     9 Yellow sand

 ispcvp(1)=3     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
 rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
 asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere
!     asphr (1)=7.0 ; asphr (2)=10.0 ; asphr (3)=1.1  ! non spherical (does not work well) [Nakajima et al. ]

 rop=2.56                 ! density [g/cm3]

 if( lognormal_psd_assum ) then

!
! 1moment case
! PSD parameter packet psd_params(6,33)=(i,j)
! 
!   i   1       2       3     4     5    6  
! j -----------------------------------------
! 1 |   r   psd_type   PSD constant parameter
!   |
! 2 | nmode   
!   |
! 3 | rmin 
!   |
! 4 | rmax
!   |
!
     psd_params(2,1)=8        ! new log normal

     if( spc == 'du1') then !these values are consistent to GOCART spectrum AOD calculations (re = rm*3)
       rm=0.25  ; re = 0.75 !mode & stdev [um]
     elseif( spc == 'du2') then
       rm=0.5   ; re = 1.5 !mode & stdev [um]
     elseif( spc == 'du3') then
       rm=1.0   ; re = 3.0 !mode & stdev [um]
     elseif( spc == 'du4') then
       rm=2.0   ; re = 6.0 !mode & stdev [um]
     elseif( spc == 'du5') then
       rm=4.0   ; re = 12.0 !mode & stdev [um]
     elseif( spc == 'ini') then
       rm=4.0   ; re = 12.0 !mode & stdev [um]
     endif

    sigma = 2. * log( re/rm ) / 5.  ! sigma = ln2(S)  (0.48 ~ 1.5)

    rm_cm = rm * 1.0d-4  !  [cm]

    !total number concentration 
    Ntot = 3. * q / ( 4.*pi*rop*(rm_cm**3) * exp( 9.*(sigma**2)/2. ) )

    psd_params(3,1) = Ntot          ! total number concentration [1/cm3]
    psd_params(4,1) = sigma         ! s: log standard deviation [-]
    psd_params(5,1) = log( rm_cm )  !mu: mean --> log(rm) [cm] <- [um] 

!print*, 'Ntot= ',Ntot



 else  !bin approach (single size) DOES NOT WORK WELL

!
! UCD aerosol case (single bin)
! PSD parameter packet psd_params(6,33)=(i,j)  
! n0: # conc per SBM bin
! 
!   i   1       2         3         4       5  
! j -----------------------------------------
! 1 |   r   psd_type      n0      rad   drad(1)
!   |
!
    psd_params(2,1)=7        ! single bin

    rbin = 0.25
    if( spc == 'du1') then !these values are consistent to GOCART spectrum AOD calculations
     rbin = 0.25   ! radius [um]
    elseif( spc == 'du2') then
     rbin = 0.5    ! radius [um]
    elseif( spc == 'du3') then
     rbin = 1.0    ! radius [um] (toshii)
    elseif( spc == 'du4') then
     rbin = 2.0    ! radius [um]
    elseif( spc == 'du5') then
     rbin = 4.0    ! radius [um]
    endif

    dr = rbin * 1.0d-4  ! bin width (artificial) [cm]
    psd_params(3,1) = 3.* q / (4. * const_dpi * (rbin**3.) * rop ) / dr  ! dN/dr [1/cm4]
    psd_params(4,1) = rbin * 1.0d-4   ! radius [cm]
    psd_params(5,1) = dr              ! dr  [cm]

!print*, 'dN/dr=', psd_params(3,1), 'Ntot', psd_params(3,1)*dr
!stop toshii


  endif

    !
! compute unified refractive index & optical properties
!
  call rfidxb_sdsu( wlc,rh,ispcvp,rfracp,asphr,rop,psd_params, &
                   nang,ang,crw1,ciw1,cextp,cabsp,php )

!
! assume non-spherical shape from Mie routine by using depolarization ratio
!
 php(nang+1) = php(nang) * dust_nonspherical_p22 ! P22 phase function at 180deg [cm-1 str-1]

!
! dust non-spherical backscatter correction 
! This suppress backscatter peak from Mie scattering about 50%. 
!
 php(nang)   = php(nang)   * dust_nospherical_ratio
 php(nang+1) = php(nang+1) * dust_nospherical_ratio



 return
 end subroutine mie_ucd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mie_sbm_visir(spc,rad,drad, rho, q, n0, wlc, rh, nang, ang, crw1, ciw1, cextp, cabsp, php  )
 implicit none
 
 character(len=*),intent(in) :: spc ! particle species character index
 real(8),intent(in):: rad(nbin)     ! radius of size bin [cm]
 real(8),intent(in):: drad(nbin)    ! width of size bin [cm]
 real(8),intent(in):: rho(nbin)     ! effective density of size bin [g/cm3]
 real(8),intent(in):: q             !total mass mixing ratio of particle for each nbin [g/cm3]
 real(8),intent(in):: n0(nbin)      ! dN/dr: # conc per radius [cm-4] for each nbin
 real(8),intent(in):: wlc           ! wave length [cm]
 real(8),intent(in):: rh            ! relativev humidity  (0~1)
 integer,intent(out):: nang         ! number of scaattering angles
 real(8),intent(out):: ang(knang)   ! scattering angles [degree]
 real(8),intent(out):: crw1         ! relative reflactive index (real part)
 real(8),intent(out):: ciw1         ! relative reflactive index (imaginary part)
 real(8),intent(out):: cextp        ! extinction cross section [cm-1]
 real(8),intent(out):: cabsp        ! absorption cross section [cm-1]
 real(8),intent(out):: php(knang)   ! volume scattering phase function [cm-1  str-1]


 integer :: ispcvp(3)     ! internal mix of fundamental materials
 real(8) :: rfracp(3)     ! dry component volume fractions
 real(8) :: asphr(3)      ! non-spherical parameters
 real(8) :: rop            ! particle density of dry mixture [g/cm3]
 real(8) :: psd_params(6,33)   ! params to define volume size dist.

 logical :: success

!
! SBM case
! PSD parameter packet psd_params(6,33)=(i,j)  
! n0: # conc per SBM bin
! 
!   i   1       2         3         4       5  
! j -----------------------------------------
! 1 |   r   psd_type   n0(1)   rad(1)   drad(1)
!   |
! 2 | nmode            n0(2)   rad(2)   drad(2)
!   |
! 3 | rmin             n0(3)   rad(3)   drad(3)
!   |
! 4 | rmax             n0(4)   rad(4)   drad(4)
! . |                     .       .         . 
! . |                     .       .         .
!33 |                  n0(33)  rad(33)  drad(33)
!

!--exec

!initialise
  cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.
  psd_params = 0. 

  if( trim(spc) /= 'ini' .and. q <= q_min_conden) return

 if( ss_opt_visir == 1 ) then

     call ss_visir_sbm_micro  (trim(spc), wlc*1.d+4, rad, drad, q, n0, rho, &
                               cextp, cabsp, php, success)
     if(success)  return  ! if not sucess go to old routine below....

 endif



!
! species-independent constant optical parameter
!
 psd_params(1,2) = 1.      !# of modes is always 1 for SDSU
 psd_params(1,3) =     0.01d-4  ! 0.01micron        (must start from this for stability purpose)
 psd_params(1,4) = 10000.0d-4   ! 10000micron = 1cm (consistent to module_opt_micro)


 psd_params(2,1)=5        ! bin-by-bin
 psd_params(3,1:nbin)= n0 (1:nbin)  ! # conc per bin [1/cm4]
 psd_params(4,1:nbin)=rad (1:nbin)  ! center radius of size bin [cm] 
 psd_params(5,1:nbin)=drad(1:nbin)  ! width of size bin [cm] 

 rop=1.0             ! density

!
! species-dependent constant optical parameter
!
  spc_select: select case(trim(spc))

  case('liq','ini')  !cloud water and rain
     ispcvp(1)=1     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions 
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere 
 
  case('ice_col') !cloud ice
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case('ice_pla') !cloud ice
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case('ice_den') !cloud ice
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere


  case('snow') !snow
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case('graupel') !graupel
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case('hail') !hail
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case default
      print*,'spc=',spc
      call stop_sdsu(myrank, 'MSG mie_sbm_visir: There is no such spc')
  end select spc_select

!
! compute unified refractive index & optical properties
!
  call rfidxb_sdsu( wlc,rh,ispcvp,rfracp,asphr,rop,psd_params, &
                   nang,ang,crw1,ciw1,cextp,cabsp,php )

!
! spherical particle assume no polarization --> P11(180deg) = P22(180deg)  
!
 php(nang+1) = php(nang)  !P22 phase function at 180deg [cm-1 str-1]


 return
 end subroutine mie_sbm_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mie_rams_visir(spc, cfmas, pwmas, q, re, gnu, wlc, rh, nang, ang, crw1, ciw1, cextp, cabsp, php  )
 implicit none
 
 character(len=3),intent(in) :: spc  ! particle species character index
 real(8),intent(in):: cfmas,pwmas !alpha_m and beta_m for mass-diameter relationship
 real(8),intent(in):: q         ! mass mixing ratio of particle [g/cm3]
 real(8),intent(in):: re        ! particle effective radius [micron]
 real(8),intent(in):: gnu       ! shape parameter for generalized gamma disribution
 real(8),intent(in):: wlc       ! wavelength [cm]
 real(8),intent(in):: rh        ! relativev humidity  (0~1)
 integer,intent(out):: nang       !! number of scaattering angles
 real(8),intent(out):: ang(knang) !! scattering angles [degree]
 real(8),intent(out):: crw1     !! relative reflactive index (real part)
 real(8),intent(out):: ciw1     !! relative reflactive index (imaginary part)
 real(8),intent(out):: cextp    ! extinction cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: cabsp    ! absorption cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: php(knang) ! volume scattering phase function [cm-1  str-1]

 integer :: ispcvp(3)     ! internal mix of fundamental materials
 real(8) :: rfracp(3)     ! dry component volume fractions
 real(8) :: asphr(3)      ! non-spherical parameters
 real(8) :: rop            ! particle density of dry mixture [g/cm3]
 real(8) :: psd_params(6,33)   ! params to define volume size dist.

 integer :: ispc

 real(8) :: gamfac , gfac1, gfac2 ! gamma parameter
 real(8) :: dn  ! characteristic diameter [m]
 real(8) :: mean_mass ! mean particle mass [kg]
 real(8) :: ntot  ! total particle number conc [#/m3]
 logical :: success

!
! 1moment case
! PSD parameter packet psd_params(6,33)=(i,j)
! 
!   i   1       2       3     4     5    6  
! j -----------------------------------------
! 1 |   r   psd_type   PSD constant parameter
!   |
! 2 | nmode   
!   |
! 3 | rmin 
!   |
! 4 | rmax
!   |
!


!--exec
!initialise
  cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
  psd_params = 0. 
 
  if( trim(spc) /= 'ini' .and. q <= q_min_conden) return  

!
! non -
!

 if( ss_opt_visir == 1 ) then

!  call interp_icescat_visir(q, re, wlc, cextp, cabsp, php, success )

  call ss_visir_rams_micro(trim(spc), wlc*1.d+4, q*1.d6, re, gnu, cfmas, pwmas, &
                            cextp, cabsp, php, success )
  if(success)  return  ! if not sucess go to old routine below....

 endif


!
! Get characteristic diameter from effective radius (re) and PSD shape parameter (gnu)
!
  call gamma_reff_r8(gnu+3.0,gfac1)
  call gamma_reff_r8(gnu+2.0,gfac2)
  gamfac = gfac1/gfac2
  dn = 2.0 * re / gamfac * 1e-6  !charactristic diameter [m]

!
! Get total number concentration for a given characteristic diameter, mass mixing ratio, PSD shape parameter
! , and gamma factors.
!
  call gamma_reff_r8(gnu+pwmas,gfac1)
  call gamma_reff_r8(gnu      ,gfac2)
  gamfac = gfac1/gfac2
  mean_mass = cfmas * ( dn ** pwmas ) * gamfac  ! mean particle mass [kg]
  ntot = q * 1.0d+3 / mean_mass                 ! total particle number concentration [#/m3]

!
! species-independent constant optical parameter
!
 psd_params(1,2) = 1.d0      !# of modes is always 1 for SDSU

 select case(trim(spc))
 case('qc1','qc2','qi1','qi2','ini')  ! non-precip
   psd_params(1,3) =     1.0d-4  ! 1micron        (must start from this for stability purpose)
   psd_params(1,4) = 100.0d-4    ! 100micron  (consistent to module_opt_micro)
 case('qr','qs','qg','qh') ! precip
   psd_params(1,3) =     1.0d-4  ! 1micron        (must start from this for stability purpose)
   psd_params(1,4) = 10000.0d-4    ! 10000micron = 1cm (consistent to module_opt_micro)
 case default
    print*, 'MSG mie_rams_visir: There is no such spc=',spc
   call stop_sdsu(myrank,'')
 end select 


 psd_params(2,1) = 6    !index for generlized gamma distribution
 psd_params(3,1) = ntot !total particle number concentration [#/m3]
 psd_params(4,1) = dn   !charactristic diameter [m]
 psd_params(5,1) = gnu  !psd shape parameter for generlized gamma distribution
 rop=1.0                ! density


!
! species-dependent constant optical parameter (only 1mode)
!
  spc_select: select case(trim(spc))

  case('qc1','qc2','qr','ini')  ! liquid
     ispcvp(1)=1     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions 
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere 

  case('qi1','qi2')    !cloud ice mode1 and mode2
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

  case('qs','qg','qh') ! riming/aggregation
     ispcvp(1)=2     ; ispcvp(2)=0    ; ispcvp(3)=0     ! 3 component mixture
     rfracp(1)=1.00  ; rfracp(2)=0.00 ; rfracp(3)=0.00  ! dry component volume fractions (sphere)
     asphr (1)=1.0d9 ; asphr (2)=10.0 ; asphr (3)=1.1   ! nonspherical parameter (x0, G, r): sphere

!    asphr (1)=7.0 ; asphr (2)=10.0 ; asphr (3)=1.1      !toshii let's give more isotropic phase function 
                                                          !but it creats NaN for very small re

  case default
      print*, 'MSG mie_rams_visir: There is no such spc=',spc
      call stop_sdsu(myrank,'') 
  end select spc_select

!
! compute unified refractive index & optical properties
!
  call rfidxb_sdsu( wlc,rh,ispcvp,rfracp,asphr,rop,psd_params, &
                   nang,ang,crw1,ciw1,cextp,cabsp,php )

!
! spherical particle assume no polarization --> P11(180deg) = P22(180deg)  
!
 php(nang+1) = php(nang)  !P22 phase function at 180deg [cm-1 str-1]

 return
 end subroutine mie_rams_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine rfidxb_sdsu( wlc,rh, ispcvp,rfracp,asphr,rop,psd_params, &
                         nang,ang,crw1,ciw1,cextp,cabsp,php)

! Get parameters for wet particle polydispersions
! Compute refractive index for internal mixture, including
! hygroscopic growth.
!--- history
! 95. 2.22 created from getpar
! 95. 9.16 modified
! 96. 2.29  2nd order polynomial interpolation (t.y.nakajima)
! 96. 3.11 allow non-spherical parameters
! 96. 7.19 debugged (takayabu, t. y. nakajima)
! 96. 7.19 log for ci and cr interpolation.(t.y.nakajima)
!     5. 5 awcrp(*,*,2)
! 97. 3.17 crw1,ciw1 on subroutine rfidxb by t.y.nakajima
! 97. 5.6  same debug as
! 08. 9.19 modified Fortran 90 free-style form
!--- input
! wlc      r          wavelength (cm)
! rh       r          relative humidity

! mp       i          aerosol number
! ispcvp i(3,kptc)    fundamental materials for 3 comp. internal mixture (1-8)
! rfracp r(3,kptc)    dry volume fraction of the dry mixture
! asphr  r(3,kptc)    non-spherical parameters (x0, g, r)
! rop    r(kptc)      paricle density relative to water
! pr  r(6,4)  dv/dlnr parameters for the dry mixture
!                     see volume_spectrum (old name vlspc2)
!                     c-values (coefficients of volume spectrum) are ralative
! nv        i         number of fundamental species (1-8)
! nwlv      i         number of wavelengths for refractive index tables
! wlv    r(kwl)       wavelengths (micron)
! rfi   r(kwl,knv,2)  refractive index of fundamental materials (mr, mi)
!                     =mr - i mi
!                     with log-regular wavelength interval
!--- output
! cextp    r          extinction cross section [cm-1]
! cabsp    r          absorption cross section [cm-1]
! nang     i          number of scattering angles
! ang    r(nang)      scattering angle in degrees
! php    r(knang)     volume scattering phase function [cm-1  str-1]
!---
!  use paras, only: kptc,kaw,kwlv,knv,knang,kintvl,kpol,nv,nwlv,wlv,rfi
  implicit none
  save

! input (nv,nwlv,wlv,rfi are global paramter from module pras)
  real(8),intent(in):: wlc      !! wavelength [cm]
  real(8),intent(in):: rh        !! relativev humidity
  integer,intent(in):: ispcvp(3)    ! internal mix of fundamental materials
  real(8),intent(in):: rfracp(3)    ! dry component volume fractions
  real(8),intent(in):: asphr(3)     ! non-spherical parameters
  real(8),intent(in):: rop          ! particle density of dry mixture [g/cm3]
  real(8),intent(inout):: psd_params(6,33)     ! parameter packet for particle size distribution (PSD) 

! output
  integer,intent(out):: nang       !! number of scaattering angles
  real(8),intent(out):: ang(knang) !! scattering angles [degree]
  real(8),intent(out):: crw1     !! relative reflactive index (real part)
  real(8),intent(out):: ciw1     !! relative reflactive index (imaginary part)
  real(8),intent(out):: cextp    !! extinction cross section [cm-1]
  real(8),intent(out):: cabsp    !! absorption cross section [cm-1]
  real(8),intent(out):: php(knang)    !! volume scattering phase function [cm-1  str-1]

! for rawb and wet grouth of aerosol
  real(8):: rmd                  !! dry aerosol radius relative to water
  real(8):: rmw                  !! wet aerosol radius relative to water
  real(8):: aw(kaw), rmmd(kaw)   !! water activity
  real(8):: rawb                 !! mode radius of wet aerosols [cm]

! for wavelength
  real(8):: wl1,dwl,pint
  integer:: iwl,iwlt

! for refractive index
  real(8):: cr0,ci0,cr01,ci01,cr1,ci1
  integer:: nintp
  real(8):: x(3),yr(3),yi(3)
  real(8):: bintp

! for gtph4b
  real(8):: x0                   !! critical size parameter for mie scaattering
  real(8):: gg                   !! asymmetry parameter for transmittetd ray
  real(8):: rp                   !! surface area fraction to that of sphere

  integer:: ipol                 !! number of polarization component
  real(8):: ph(knang,kpol)       !! volume phase function
  real(8):: cext                 !! extinction cross section [cm-1]
  real(8):: cabst                 !! absorption cross section [cm-1]
  real(8):: cg                   !! geometrical cross section [cm-1]
  real(8):: vl                   !! volume [cm3/cm3]
  integer:: intvl                !! number of size interval in kernel file
  real(8):: szp(kintvl)          !! size parameter

! for work
  integer:: m,isp,is
!--exec


! non-spherical parameters
  x0=asphr(1)
  gg=asphr(2)
  rp=asphr(3)

! initialization
  cextp=0.d0 ; cabsp=0.d0 ; php(1:knang)=0.d0

!
! compute internal mixture
!

! search wavelength in refractive index table
     wl1=wlc*1.0e4   ! wavelength [cm] -> [micron] 
     dwl=log(wl1/wlv(1))/log(wlv(nwlv)/wlv(1))*(nwlv-1)+1  !
     iwl=int(dwl) !wavenumber index

! 3 dry aerosol components refractive index:
     cr0=0.d0; ci0=0.d0

     do isp=1,3
        is=ispcvp(isp)    ! internal mix of fundamental materials
        if(is==0) exit

! interpolation to get corresponding cr1 and ci1 to wl1

        nintp=3
        if(nintp>nwlv) nintp=nwlv
        iwlt=iwl

        if(iwl>=nwlv-1) iwlt=nwlv-2

        !
        ! Filter for wavelength out of range indicated by identifying negative iwlt
        !
        if( iwlt <= 0 ) then 
!           print*,'MSG rfidxb_sdsu: Ooops iwlt is negative --> wavelength may be too small, return' 
           cextp        = dble(undefined)  ! extinction cross section [cm-1]
           cabsp        = dble(undefined)  ! absorption cross section [cm-1]
           php(1:knang) = dble(undefined)  ! volume scattering phase function
           return
        endif

        x (1:nintp)=log(wlv(iwlt:iwlt+nintp-1))
        yr(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,is,1) )
        yi(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,is,2) )
        pint=log(wl1)

        cr01=exp( bintp(pint,nintp,x,yr) )
        ci01=exp( bintp(pint,nintp,x,yi) )

! dry aerosol refractive index: (considering homogeneous mixture)
        cr0=cr0+cr01*rfracp(isp) ! just multiplying dry component volume fractions 
        ci0=ci0+ci01*rfracp(isp)
     enddo !isp loop

     if(cr0<=0.d0) then
        call stop_sdsu(myrank, 'MSG rfidxb_sdsu: cr0<= 0. particle species does not exist')
     endif

! water refractive index
     yr(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,1,1) )
     yi(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,1,2) )
     pint=log(wl1)

     crw1=exp( bintp(pint,nintp,x,yr) )  !output needed for water surface
     ciw1=exp( bintp(pint,nintp,x,yi) )  !output needed for water surface

! wet aerosol refractive index: (this values accounts for internal & water mixing)
     cr1=cr0
     ci1=ci0

!
! up to here, the code was designed to get total (internal mixture) refractive index
!

! get optical parameters toshii
     call gtph4b_org(wlc,cr1,ci1,psd_params,x0,gg,rp,nang,ipol,ang,ph, &
                 cext,cabst,cg,vl,intvl,szp)
     cextp=cext  ! extinction cross section [cm-1]
     cabsp=cabst  ! absorption cross section [cm-1]
     php(1:knang)=ph(1:knang,1) ! volume scattering phase function  [cm-1  str-1]

  return
 end subroutine rfidxb_sdsu

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gtph4b_org(wlc,cr3,ci3,psd_params,x0,gg,rp,nang,ipol,ang,ph, &
                   cext,cabst,cg,vl,intvl,szp)
! get phase function and cross section 
! from table of interpolated refractive index
! with log-regular grid wavelengths.
! cr3,ci3 assigned; assume ipol=1
! pollack and cuzzi method implemented
!--- history
! 95. 9.14  created from getph3
! 96. 2.29  2nd order polynomial interpolation (t.y.nakajima)
! 96. 3.11  pollack and cuzzi method implemented (noguchi, tokai)
! 96. 4. 4  dcg(ir)=dcg(ir)+cg1 -> f(indr.le.1) then ... else...endif
! 08. 9.19 modified Fortran 90 free-style form
!--- input
! wlc     r      wavelength in cm
! cr3     r      real part of refractive index
! ci3     r      imaginary part of refractive index (positive value)
! pr   r(6,4)    parameter packet of size distribution (see volume_spectrum)
! x0      r      critical size parameter for mie scattering (pollack-cuzzi)
!                give 1.0e9 for mie theory only calculations.
! gg      r      asymmetry parameter for transmitted ray (pollcack-cuzzi)
!                any value for mie theory only calculations.
! rp      r      surface area fraction to that of sphere (about 1.1-1.3)
!                any value for mie theory only calculations.
!--- output
! nang    i      nbr of scattering angle
! intvl   i      nbr of size parameter
! ipol    i      nbr of polarization component
! ang  r(nang)   scattering angle in degree
! ph   r(knang,ipol) volume phase function
! cext   r       extinction cross section (cm-1)
! cabst   r       absorption cross section (cm-1)
! cg     r       geometrical cross section (cm-1)
! vl     r       volume (cm3/cm3)
! intvl  i       number of size interval in the kernel file
! szp  r(intvl)
!---
  !use paras, only: knang, kintvl,kpol,pi,rad
  implicit none
  save

! for input
  integer,parameter :: iuk=14       !! device number
  real(8),intent(in):: wlc       !! wavelength [cm]
  real(8),intent(in):: cr3       !! refractive index (real part)
  real(8),intent(in):: ci3       !! refractive index (imaginary part)
  real(8),intent(inout):: psd_params(6,33)   !! PSD parameter packet
  real(8),intent(in):: x0        !! critical size parameter for mie scaattering
  real(8),intent(in):: gg        !! asymmetry parameter for transmittetd ray
  real(8),intent(in):: rp        !! surface area fraction to that of sphere

! for output
  integer,intent(inout):: nang     !! number of scattering angles
  integer,intent(inout):: ipol     !! number of polarization component
  real(8),intent(inout):: ang(knang)   !! scattering angles [degree]
  real(8),intent(out):: ph(knang,kpol) !! volume phase function
  real(8),intent(out):: cext         !! extinction cross section
  real(8),intent(out):: cabst         !! absorption cross section
  real(8),intent(out):: cg           !! geometrical cross section
  real(8),intent(out):: vl           !! volume
  integer,intent(inout):: intvl        !! number of size interval in kernel file
  real(8),intent(inout):: szp(kintvl)  !! size parameter
! v
  integer :: nang_save
  real(8) :: ang_save(knang)

! initialization
  integer:: init=1
  real(8),save :: del,rszp2,phm(knang)

! kernel parameters
  integer,parameter:: krf=200        !! number of refractive index
  integer,parameter:: kr=20          !! number of real part in each section
  integer,parameter:: ki=20          !! number of imaginary part in each sec
  integer,parameter:: ksec=3         !! number of section
  integer,parameter:: knang2=knang+2

  integer:: k,i,j,l,ip,jp,kp,ipt,jpt,ir,ix,iy
  integer:: ierr
  integer:: nrf
  character:: ch*1
  real(8):: cr(krf),ci(krf)
  real(8),allocatable :: q(:,:,:,:) !dimension (kintvl,kpol,krf,knang2)
  real(8):: dcr,crt(kr,ksec),cit(ki,ksec)
  integer:: ict(kr,ki,ksec)
  integer:: ns,nt
  integer:: nav

  integer:: nintpr,nintpi,irf(3,3)
  real(8):: xr(3),xi(3)

  integer:: intvl_rayleigh,intvl_tot,ir8
  integer:: indr
  real(8):: r1,r2,r_center,dr,dlr,v1,dv_dlnr,cg1,rmin_mie
  real(8):: volume_spectrum

  real(8):: cof
  real(8):: xi2(3),y(3),zi(3),z,dci
  integer:: ic,irfi
  real(8):: bintp,pint

! for smllop 
  real(8):: x,qext_rayleigh,qsca_rayleigh,qabs_rayleigh

! for nons
  integer:: ix0
  real(8):: cr1,ci1
  real(8):: dcgl,dcgs
  real(8):: rmin,rmax
  real(8):: wvn

  real(8):: dabs(kintvl)
  real(8):: dcg(kintvl)
  real(8):: dext(kintvl)
  real(8):: phs(knang)


! kernel data parameters
  integer:: nr(1:ksec)  =(/    3,     9,    18/)
  integer:: ni(1:ksec)  =(/   16,    10,     4/)
  real(8):: crmn(1:ksec)=(/1.3d0, 1.0d0, 1.0d0/)
  real(8):: crmx(1:ksec)=(/1.5d0, 1.8d0, 2.7d0/)
  real(8):: cimn(1:ksec)=(/1.d-9, 1.d-4, 1.d-1/)
  real(8):: cimx(1:ksec)=(/1.d-4, 1.d-1, 1.0d0/)

! for kernal file
  character (len=200) :: path_file
  logical :: logic_opened

!--exec
!  if(init>0) then
!     init=0

 path_file=trim(sdsu_dir_data)//'KRNL.OUT'
 inquire(file=trim(path_file),opened=logic_opened)

 if( .not. logic_opened) then

    if(masterproc .and. verbose_SDSU) print*,'MSG gtph4b: open/reading KRNL.OUT'
!    open(iuk,file=trim(path_file),status='old')
    call open_oldfile(iuk, trim(path_file) )

!
! allocate memory
!
    if( .not. allocated(q) ) then
      allocate( q(1:kintvl,1:kpol,1:krf,1:knang2) , stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank, 'MSG gtph4b: allocation error -> Terminate program. ')
    endif

     rewind iuk
! get kernel ang, szp, ak, cext, cabs
     read(iuk,*) intvl, nang, ipol  ! 73, 74, 1  
     if (nang > knang) then
        call stop_sdsu(myrank, 'nang > knang')
     end if
     if (ipol > kpol) then
        call stop_sdsu(myrank, 'ipol > kpol' )
     end if
     if(ipol/=1) then
        call stop_sdsu(myrank, 'assume only ipol=1 case in getph4')
     endif

     read(iuk,*) szp(1:intvl)  ! size parameter (0.1371 ~ 938.0: Mie regime)
     read(iuk,*) ang(1:nang)   ! scattering angle [deg] (0-180deg) 0: forward 180: backscatter


     del=log(szp(intvl)/szp(1))/dble(intvl-1)
     rszp2=exp(del*0.5d0)
     phm(1:nang)=3.0/16.0/pi*(1.d0+cos(ang(1:nang)*rad)**2)  !rayleigh phase function

!! read kernel
     do ir=1,krf
        read(iuk,'(a1)',iostat=ierr) ch !comments
        if(ierr/=0) then
           nrf=ir; exit
        endif
        read(iuk,*) cr(ir),ci(ir)
        ci(ir)=abs(ci(ir))
        do i = 1, intvl   ! 1~73
           do j = 1, ipol  ! 1~1
              read(iuk,*) q(i,j,ir,1:nang) 
           enddo
        enddo
        read(iuk,*) q(1:intvl,1,ir,nang+1:nang+2)
     enddo



!! get refractive index tables (3 tables)
     do k=1,ksec
        dcr=(crmx(k)-crmn(k))/dble(nr(k)-1)
        do i=1,nr(k)
           cr1=crmn(k)+dcr*dble(i-1)
           if(abs(cr1-1.0) <= 0.d0) cr1=1.01
           crt(i,k)=cr1
        enddo

        dci=log(cimx(k)/cimn(k))/dble(ni(k)-1)
        do i=1,ni(k)
           ci1=dci*dble(i-1)
           cit(i,k)=cimn(k)*exp(ci1)
        enddo
     enddo


!! find grids for kernel refractive indices
     do k=1,ksec
        ict(1:nr(k),1:ni(k),k)=0
     enddo

     do ir=1,nrf
        cr1=cr(ir); ci1=ci(ir)
        do k=1,ksec
           do i=1,nr(k)
              if(abs(crt(i,k)/max(cr1,1.d-18)-1.d0) <= 0.01d0) exit
           enddo
           do j=1,ni(k)
              if(abs(cit(j,k)/max(ci1,1.d-18)-1.d0) <= 0.01d0) exit
           enddo
           ict(i,j,k)=ir
        enddo
     enddo

     do k=1,ksec; do i=1,nr(k); do j=1,ni(k)
        if(ict(i,j,k)==0) then
           call stop_sdsu(myrank, 'no matching grid')
        endif
     enddo; enddo; enddo

! number of averaging

     nav=10

  nang_save = nang
  ang_save(1:nang) = ang(1:nang)

  endif  !logic_opened

!
! initialization end 
!

!
! Make sure you have nang & ang without initialization
!
  nang        = nang_save
  ang(1:nang) = ang_save(1:nang)

  wvn=2.d0*pi/wlc  ! wavenumber [1/cm]

! find the domain
  do kp=1,ksec
     if(ci3 <= cimx(kp).or.kp==ksec) then
        ns=nr(kp); nt=ni(kp)
        exit
     endif
  enddo


  do ip=1,ns-1
     if(cr3 <= crt(ip,kp)) exit
  enddo
  do jp=1,nt-1
     if(ci3.le.cit(jp,kp)) exit
  enddo

  ip=max(1, ip-1)
  jp=max(1, jp-1)

! get grid points
!! 96.2.29
  nintpr=3
  nintpi=3
  if(nintpr > nr(kp)) nintpr=nr(kp)
  if(nintpi > ni(kp)) nintpi=ni(kp)

  ipt=ip
  jpt=jp
  if(ip >= nr(kp)-1) ipt=nr(kp)-2
  if(jp >= ni(kp)-1) jpt=ni(kp)-2
  xr(1:nintpr)=cr(ict(ipt:ipt+nintpr-1,jpt             ,kp))
  xi(1:nintpi)=ci(ict(ipt             ,jpt:jpt+nintpi-1,kp))

  irf(1:nintpr,1:nintpi)=ict(ipt:ipt+nintpr-1,jpt:jpt+nintpi-1,kp)

! get interpolated optical constants
  cg=0.d0; vl=0.d0
  cext=0.d0; cabst=0.d0
  dcgl=0.d0; dcgs=0.d0

  do ix0=1,intvl
     if(szp(ix0) > x0) exit
  enddo
  phs(1:nang)=0.d0; ph(1:nang,1)=0.d0

  rmin=psd_params(1,3)   !minimum radius [cm] 
  rmax=psd_params(1,4)   !maximum radius  [cm]

  rmin_mie=szp(1)/wvn/rszp2  ! minimum raidus from size parameter LUT [cm]

  if(rmin_mie > rmin) then
     intvl_rayleigh=int(log(rmin_mie/rmin)/del)+1 ! define additional interval for rayleigh regime
     r2=rmin_mie*exp(-intvl_rayleigh*del)  ! define minimum size to start size loop
  else
     intvl_rayleigh=0
  endif
  intvl_tot=intvl+intvl_rayleigh ! entire size bin for loop


  size_loop: do ir8=1, intvl_tot  !size loop
     ir=ir8-intvl_rayleigh

     scat_regime: if(ir <= 0) then
        indr=1  ! rayleigh regime
        r1=r2
        r2=r1*rszp2**2
     else
        dcg(ir)=0.d0
        indr=2   ! mie regimea
        ! compute bounds from size parameter (szp)
        ! so this interval change for different wavelength...
        r1=szp(ir)/wvn/rszp2  ![cm]
        r2=szp(ir)/wvn*rszp2  ![cm]
     endif scat_regime

     dr=(r2-r1)/dble(nav) ! size interval (bin width) for sub-interval [cm]

     r_center=r1-dr*0.5d0       ! center point of radius

     v1=0.d0

     !
     ! this sub-size loop is to more precisely calculate total volume concentrations
     ! between size bins. 
     !
     sub_loop: do j=1,nav  !1~10 averaging loop

        r_center=r_center+dr  !increment of every dr

        dlr=log(r_center/(r_center-dr))    ! d ln r [-]
        psd_params(1,1) = r_center         ! radius of particle [cm]

        !
        ! compute volume spectrum of partile polydisperison: v(x) = dv / d ln r
        ! This contains all size parameter packet PSD type and size bounds. 
        !
        dv_dlnr = volume_spectrum(psd_params)      !  volume spectrum of partile polydisperison: v(x) = dv / d ln r

        cg1=dv_dlnr*3.d0/4.d0/r_center*dlr
        cg=cg+cg1

        vl = vl + dv_dlnr * dlr ! total volume (output) = dv = dv / d ln r * (d ln r) [cm3/cm3]

        v1 = v1 + dv_dlnr     ! integrating volume spectrum  [cm3/cm3]

!print*, 'ir8=',ir8,'j=',j,'r_center=',r_center,'dv_dlnr=',dv_dlnr

        if(ir <= ix0)then
           dcgs=dcgs+cg1
        else
           dcgl=dcgl+cg1
        endif

        if(indr <= 1) then  ! rayleigh regime
           dcg(1)=dcg(1)+cg1
        else
           dcg(ir)=dcg(ir)+cg1
        endif
     enddo sub_loop

     v1=v1/dble(nav)  ! derive mean value of volume spectrum [cm3/cm3]

!toshii
!write(*,FMT='(A3,1x,I3,1x,A3,1x,F9.4,1x,A3,1x,E9.4,1x,A3,1x,E9.4)' ) 'ir=',ir8,'r_center=',r_center*1.d+4,'v1=',v1,'dv_dlnr=', dv_dlnr

!write(*,FMT='(I3,1x,F9.3,1x, F9.3, 2x, E9.3, 2x, E9.3, 2x, E9.3, 2x, E9.3, 2x, E9.3)') &


    !
    ! interpolation for output
    !
    rayleigh: if(indr==1) then !rayleigh regime

        x=(r1+r2)*0.5d0*wvn  ! size parameter [-]
        call smllop(x,cr3,ci3,qext_rayleigh,qsca_rayleigh) ! for small opitcal depths
        qabs_rayleigh = qext_rayleigh - qsca_rayleigh    !absorption efficiency factor
        cof = del*3.d0/4.d0/x*v1*wvn

        !
        ! integrate rayleigh single-scattering parameter over size (output)
        !
        cext  = cext  + qext_rayleigh * cof
        cabst = cabst + qabs_rayleigh * cof
        ph(1:nang,1) = ph(1:nang,1) + phm(1:nang) * cof * qsca_rayleigh

     else  !Mie regime

        do l=1,nang+2
          ! interpolation to get corresponding q to cr3 and ci3
           xi2(1:nintpi)=log(xi(1:nintpi))
           do ix=1,nintpr
              ic=0
              do iy=1,nintpi
                 ic=ic+1
                 irfi=irf(ix,iy)
                 y(ic)= log(q(ir, 1, irfi, l))
              enddo
              pint=log(ci3)
              zi(ix)= bintp(pint,nintpi,xi2,y)
           enddo
           z=bintp(cr3,nintpr,xr,zi)
           z=exp(z)*v1*wvn

           !
           ! integrate single-scattering parameter over size (output)
           !
           if(l <= nang) then
              ph(l,1) = ph(l,1) + z
              if(ir <= ix0) phs(l) = phs(l) + z
           else if(l == nang+1) then
              cext = cext + z
              dext(ir) = z  !for nons
           else if(l==nang+2) then
              cabst=cabst+z  ! volume absorption cross section (output) [cm2/cm3] or [cm-1]
              dabs(ir)=z    ! for nons
           endif
        enddo

     endif rayleigh

  enddo size_loop

!
! semi-empirical parameterization for non-spherical particles
!
  if(indr == 1) return ! rayleigh regime
  call nons(intvl,ix0,nang,nav,cg,cr3,ci3,dcgs,dcgl,gg,rmax,rmin, &
       rp,wvn,x0,ang,dabs,dcg,dext,phs,szp,ph,cext,cabst)
  return
 end subroutine gtph4b_org

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gtph4b(wlc,cr3,ci3,psd_params,x0,gg,rp,nang,ipol,ang,ph, &
                   cext,cabst,cg,vl,intvl,szp)
! get phase function and cross section 
! from table of interpolated refractive index
! with log-regular grid wavelengths.
! cr3,ci3 assigned; assume ipol=1
! pollack and cuzzi method implemented
!--- history
! 95. 9.14  created from getph3
! 96. 2.29  2nd order polynomial interpolation (t.y.nakajima)
! 96. 3.11  pollack and cuzzi method implemented (noguchi, tokai)
! 96. 4. 4  dcg(ir)=dcg(ir)+cg1 -> f(indr.le.1) then ... else...endif
! 08. 9.19 modified Fortran 90 free-style form
!
! 07/2015 Toshi Matsui@NASA GSFC: Adding geometric regime and
!                                 associated significant changes in coding. 
!
!--- input
! wlc     r      wavelength in cm
! cr3     r      real part of refractive index
! ci3     r      imaginary part of refractive index (positive value)
! pr   r(6,4)    parameter packet of size distribution (see volume_spectrum)
! x0      r      critical size parameter for mie scattering (pollack-cuzzi)
!                give 1.0e9 for mie theory only calculations.
! gg      r      asymmetry parameter for transmitted ray (pollcack-cuzzi)
!                any value for mie theory only calculations.
! rp      r      surface area fraction to that of sphere (about 1.1-1.3)
!                any value for mie theory only calculations.
!--- output
! nang    i      nbr of scattering angle
! intvl   i      nbr of size parameter
! ipol    i      nbr of polarization component
! ang  r(nang)   scattering angle in degree
! ph   r(knang,ipol) volume phase function [cm-1  str-1]
! cext   r       extinction cross section (cm-1)  [cm2/cm3] or [cm-1]
! cabst   r       absorption cross section (cm-1) [cm2/cm3] or [cm-1]
! cg     r       geometrical cross section (cm-1) [cm2/cm3] or [cm-1]
! vl     r       volume (cm3/cm3)
! intvl  i       number of size interval in the kernel file
! szp  r(intvl)
!---
  !use paras, only: knang, kintvl,kpol,pi,rad
  implicit none
  save

! for input
  integer,parameter :: iuk=14       !! device number
  real(8),intent(in):: wlc       !! wavelength [cm]
  real(8),intent(in):: cr3       !! refractive index (real part)
  real(8),intent(in):: ci3       !! refractive index (imaginary part)
  real(8),intent(inout):: psd_params(6,33)   !! PSD parameter packet
  real(8),intent(in):: x0        !! critical size parameter for mie scaattering
  real(8),intent(in):: gg        !! asymmetry parameter for transmittetd ray
  real(8),intent(in):: rp        !! surface area fraction to that of sphere

! for output
  integer,intent(inout):: nang     !! number of scattering angles
  integer,intent(inout):: ipol     !! number of polarization component
  real(8),intent(inout):: ang(knang)   !! scattering angles [degree]
  real(8),intent(out):: ph(knang,kpol) !! volume phase function
  real(8),intent(out):: cext         !! extinction cross section
  real(8),intent(out):: cabst         !! absorption cross section
  real(8),intent(out):: cg           !! geometrical cross section
  real(8),intent(out):: vl           !! volume (cm3/cm3)
  integer,intent(inout):: intvl        !! number of size interval in kernel file
  real(8),intent(inout):: szp(kintvl)  !! size parameter
! v
  integer :: nang_save
  real(8) :: ang_save(knang)

! initialization
  integer:: init=1
  real(8),save :: del,rszp2,phm(knang)

! kernel parameters
  integer,parameter:: krf=200        !! number of refractive index
  integer,parameter:: kr=20          !! number of real part in each section
  integer,parameter:: ki=20          !! number of imaginary part in each sec
  integer,parameter:: ksec=3         !! number of section
  integer,parameter:: knang2=knang+2

  integer:: k,i,j,l,ip,jp,kp,ipt,jpt,ir,ix,iy
  integer:: ierr
  integer:: nrf
  character:: ch*1
  real(8):: cr(krf),ci(krf)
  real(8),allocatable :: q(:,:,:,:) !dimension (kintvl,kpol,krf,knang2)
  real(8):: dcr,crt(kr,ksec),cit(ki,ksec)
  integer:: ict(kr,ki,ksec)
  integer:: ns,nt
  integer:: nav

  integer:: nintpr,nintpi,irf(3,3)
  real(8):: xr(3),xi(3)

  integer:: intvl_rayleigh, intvl_geom, intvl_tot, ir8
  integer:: indr
  real(8):: r1,r2,r_center,dr,dlr,v1,dv_dlnr,cg1,rmin_mie, rmax_mie
  real(8):: volume_spectrum

  real(8):: cof
  real(8):: xi2(3),y(3),zi(3),z,z_bin,dci
  integer:: ic,irfi
  real(8):: bintp,pint

! for smllop 
  real(8):: x              !shape parameter
  real(8):: qext_rayleigh  !extinction efficiency [-] 
  real(8):: qsca_rayleigh  !scattering efficiency [-] 
  real(8):: qabs_rayleigh  !absorption efficiency [-] 


! for nons
  integer:: ix0
  real(8):: cr1,ci1
  real(8):: dcgl,dcgs
  real(8):: rmin,rmax
  real(8):: wvn

  real(8):: dabss(kintvl)
  real(8):: dcg(kintvl)
  real(8):: dext(kintvl)
  real(8):: phs(knang)

! kernel data parameters
  integer:: nr(1:ksec)  =(/    3,     9,    18/)
  integer:: ni(1:ksec)  =(/   16,    10,     4/)
  real(8):: crmn(1:ksec)=(/1.3d0, 1.0d0, 1.0d0/)
  real(8):: crmx(1:ksec)=(/1.5d0, 1.8d0, 2.7d0/)
  real(8):: cimn(1:ksec)=(/1.d-9, 1.d-4, 1.d-1/)
  real(8):: cimx(1:ksec)=(/1.d-4, 1.d-1, 1.0d0/)

! for kernal file
  character (len=200) :: path_file
  logical :: logic_opened

! for geometric opt
  integer :: ir_mie

! for non-pherical ice particle
  real(8) :: num   !number concentrations [#/cm3]
  real(8) :: across ! total cross sectional area [cm2/cm3]

!--exec
!  if(init>0) then
!     init=0

 path_file=trim(sdsu_dir_data)//'KRNL.OUT'
 inquire(file=trim(path_file),opened=logic_opened)

 if( .not. logic_opened) then

    if(masterproc .and. verbose_SDSU) print*,'MSG gtph4b: open/reading KRNL.OUT'
!    open(iuk,file=trim(path_file),status='old')
    call open_oldfile(iuk, trim(path_file) )

!
! allocate memory
!
    if( .not. allocated(q) ) then
      allocate( q(1:kintvl,1:kpol,1:krf,1:knang2) , stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank, 'MSG gtph4b: allocation error -> Terminate program. ')
    endif

     rewind iuk
! get kernel ang, szp, ak, cext, cabs
     read(iuk,*) intvl, nang, ipol  ! 73, 74, 1  
     if (nang > knang) then
        call stop_sdsu(myrank, 'nang > knang')
     end if
     if (ipol > kpol) then
        call stop_sdsu(myrank, 'ipol > kpol' )
     end if
     if(ipol/=1) then
        call stop_sdsu(myrank, 'assume only ipol=1 case in getph4')
     endif

     read(iuk,*) szp(1:intvl)  ! size parameter (0.1371 ~ 938.0: Mie regime)
     read(iuk,*) ang(1:nang)   ! scattering angle [deg] (0-180deg) 0: forward 180: backscatter

     del=log(szp(intvl)/szp(1))/dble(intvl-1)
     rszp2=exp(del*0.5d0)
     phm(1:nang)=3.0/16.0/pi*(1.d0+cos(ang(1:nang)*rad)**2)  !rayleigh phase function

!! read kernel
     do ir=1,krf
        read(iuk,'(a1)',iostat=ierr) ch !comments
        if(ierr/=0) then
           nrf=ir; exit
        endif
        read(iuk,*) cr(ir),ci(ir)
        ci(ir)=abs(ci(ir))
        do i = 1, intvl   ! 1~73
           do j = 1, ipol  ! 1~1
              read(iuk,*) q(i,j,ir,1:nang) 
           enddo
        enddo
        read(iuk,*) q(1:intvl,1,ir,nang+1:nang+2)
     enddo



!! get refractive index tables (3 tables)
     do k=1,ksec
        dcr=(crmx(k)-crmn(k))/dble(nr(k)-1)
        do i=1,nr(k)
           cr1=crmn(k)+dcr*dble(i-1)
           if(abs(cr1-1.0) <= 0.d0) cr1=1.01
           crt(i,k)=cr1
        enddo

        dci=log(cimx(k)/cimn(k))/dble(ni(k)-1)
        do i=1,ni(k)
           ci1=dci*dble(i-1)
           cit(i,k)=cimn(k)*exp(ci1)
        enddo
     enddo


!! find grids for kernel refractive indices
     do k=1,ksec
        ict(1:nr(k),1:ni(k),k)=0
     enddo

     do ir=1,nrf
        cr1=cr(ir); ci1=ci(ir)
        do k=1,ksec
           do i=1,nr(k)
              if(abs(crt(i,k)/max(cr1,1.d-18)-1.d0) <= 0.01d0) exit
           enddo
           do j=1,ni(k)
              if(abs(cit(j,k)/max(ci1,1.d-18)-1.d0) <= 0.01d0) exit
           enddo
           ict(i,j,k)=ir
        enddo
     enddo

     do k=1,ksec; do i=1,nr(k); do j=1,ni(k)
        if(ict(i,j,k)==0) then
           call stop_sdsu(myrank, 'no matching grid')
        endif
     enddo; enddo; enddo

! number of averaging
! toshii
     nav=10   

  nang_save = nang
  ang_save(1:nang) = ang(1:nang)

  endif  !logic_opened

!
! initialization end 
!

!
! Make sure you have nang & ang without initialization
!
  nang        = nang_save 
  ang(1:nang) = ang_save(1:nang) 

  wvn=2.d0*pi/wlc  ! wavenumber [1/cm]

! find the domain
  do kp=1,ksec
     if(ci3 <= cimx(kp).or.kp==ksec) then
        ns=nr(kp); nt=ni(kp)
        exit
     endif
  enddo

  do ip=1,ns-1
     if(cr3 <= crt(ip,kp)) exit
  enddo
  do jp=1,nt-1
     if(ci3.le.cit(jp,kp)) exit
  enddo

  ip=max(1, ip-1)
  jp=max(1, jp-1)

! get grid points
!! 96.2.29
  nintpr=3
  nintpi=3
  if(nintpr > nr(kp)) nintpr=nr(kp)
  if(nintpi > ni(kp)) nintpi=ni(kp)

  ipt=ip
  jpt=jp
  if(ip >= nr(kp)-1) ipt=nr(kp)-2
  if(jp >= ni(kp)-1) jpt=ni(kp)-2
  xr(1:nintpr)=cr(ict(ipt:ipt+nintpr-1,jpt             ,kp))
  xi(1:nintpi)=ci(ict(ipt             ,jpt:jpt+nintpi-1,kp))

  irf(1:nintpr,1:nintpi)=ict(ipt:ipt+nintpr-1,jpt:jpt+nintpi-1,kp)

!
! initialie optical properties
! get interpolated optical constants
!
  cg=0.d0; vl=0.d0
  cext=0.d0; cabst=0.d0
  dcgl=0.d0; dcgs=0.d0
  phs(1:nang)=0.d0; ph(1:nang,1)=0.d0

  do ix0=1,intvl
!     print*, ix0, 'ranges of size parameter',szp(ix0)  !0.1 ~ 938
     if(szp(ix0) > x0) exit
  enddo

  rmin = psd_params(1,3)   !minimum radius [cm] 
  rmax = psd_params(1,4)   !maximum radius  [cm]

  rmin_mie = szp(1)/wvn/rszp2     ! minimum radius from size parameters LUT.  [cm]
  rmax_mie = szp(intvl)/wvn/rszp2 ! maximum radius from size parameters LUT. [cm]

  if(rmin_mie > rmin) then
     intvl_rayleigh = int(log(rmin_mie/rmin)/del) + 1 ! define additional interval for rayleigh regime
     r2 = rmin_mie*exp(-intvl_rayleigh*del)  ! define minimum size to start size loop
  else
     intvl_rayleigh = 0
  endif

 if( rmax > rmax_mie ) then
     intvl_geom = int(log(rmax/rmax_mie)/del) - 1  ! define additional interval for geometric regime
 endif

  intvl_tot = intvl + intvl_rayleigh  + intvl_geom  ! entire size bin for loop

!toshii
! print*, 'intvl',  intvl_tot , intvl_rayleigh , intvl_geom  ! entire size bin for loop

!
! let's make algorithm to extend the size loop up to the rmax (1cm for precip species). 
!  
  size_loop: do ir8 = 1, intvl_tot   !size loop of original LUT database (73 bin)

     ir = ir8 - intvl_rayleigh

     scat_regime: if(ir <= 0) then  ! rayleigh scattering regime

        indr = 1    ! rayleigh regime
        r1   = r2
        r2   = r1*rszp2**2

     elseif( ir8 > intvl + intvl_rayleigh ) then ! geometric regime

       indr = 3    ! geometric regime
       r1 = r2
       r2 = r1*rszp2**2 

     else  ! mie scattering regime 

        dcg(ir)=0.d0
        indr=2    ! mie regime
        ! compute bounds from size parameter (szp)
        ! so this interval change for different wavelength...
        r1 = szp(ir) / wvn / rszp2   ![cm]
        r2 = szp(ir) / wvn * rszp2   ![cm]

     endif scat_regime

     dr = (r2-r1)/dble(nav) ! size interval (bin width) for sub-interval [cm]

     r_center = r1-dr*0.5d0       ! center point radius [cm]

     v1=0.d0

     !
     ! this sub-size loop is to more precisely calculate total volume concentrations
     ! between size bins. 
     !
     sub_loop: do j=1,nav  !nav = 10  

        r_center=r_center+dr  !increment of every dr 

        dlr = log(r_center/(r_center-dr))  ! d ln r
        psd_params(1,1) = r_center   ! radius of particle [cm]

        !
        ! compute volume spectrum of partile polydisperison: v(x) = dv / d ln r
        ! This contains all size parameter packet PSD type and size bounds. 
        !
        dv_dlnr = volume_spectrum(psd_params)      ! dv_dlnr= v(x) = dv / d ln r  [cm3/cm3]

        !
        ! Note that dlr is dimentionless: dlr = d ln( r ) = ln ( r/r_ref )  ! where r_reff = 1cm
        !
        vl = vl + dv_dlnr * dlr  ! total volume (output) = dv = dv / d ln r * (d ln r) [cm3/cm3]

        v1=v1+dv_dlnr  !integrating volume spectrum  [cm3/cm3]

        ! 
        ! for nons ( empirical correction for non-spherical scattering )
        !
!        cg1=dv_dlnr*3.d0/4.d0/r_center*dlr
!        cg=cg+cg1
!        if(ir <= ix0)then
!           dcgs=dcgs+cg1  !for nons
!        else
!           dcgl=dcgl+cg1  !for nons
!        endif
!        if(indr <= 1) then
!           dcg(1)=dcg(1)+cg1  !for nons
!        else
!           dcg(ir)=dcg(ir)+cg1 !for nons
!        endif

     enddo sub_loop

     v1 = v1/dble(nav)  ! derive mean value of volume spectrum [cm3/cm3]

     num = v1 / ( 4. * 3.14159 * ( r_center**3. ) / 3. )   ! number concentrations [#/cm3]

!total cross sectional area
     across = num * 3.14159 * (r_center**2.)  ! [cm2/cm3]

    !
    ! interpolation for output
    !
     rayleigh: if(indr==1) then  ! Rayleigh regime

        x=(r1+r2)*0.5d0*wvn  ! size parameter [-]
        call smllop(x,cr3,ci3,qext_rayleigh,qsca_rayleigh)  !for small opitcal depths
        qabs_rayleigh = qext_rayleigh - qsca_rayleigh !absorption efficiency factor
        cof = del*3.d0/4.d0/x * v1 * wvn

        !
        ! integrate rayleigh single-scattering parameter over size (output)
        !
        cext  = cext+qext_rayleigh*cof        ! extinction cross section (output) [cm2/cm3] or [cm-1]
        cabst = cabst+qabs_rayleigh*cof       ! absorption cross section (output) [cm2/cm3] or [cm-1]
        ph(1:nang,1) = ph(1:nang,1)+phm(1:nang)*cof*qsca_rayleigh  ! rayleigh phase function (output)  [cm-1  str-1]

!print*,ir8,'Cext [cm2/cm3] ', qext_rayleigh*cof , 'Ray Qext [-]=', qext_rayleigh*cof/across

     elseif( indr > 1) then  ! Mie or Geometric regime 

        if( indr == 2 ) ir_mie = ir    !mie
        if( indr == 3 ) ir_mie = intvl !geometric assumption using maximum value of Mie

        do l=1,nang+2

           ! interpolation to get corresponding q to cr3 and ci3
           xi2(1:nintpi)=log(xi(1:nintpi))
           do ix=1,nintpr
              ic=0
              do iy=1,nintpi
                 ic=ic+1
                 irfi=irf(ix,iy)
                 y(ic)= log(q(ir_mie, 1, irfi, l))
              enddo
              pint=log(ci3)
              zi(ix)= bintp(pint,nintpi,xi2,y)
           enddo
           z = bintp(cr3,nintpr,xr,zi)
           z = exp(z) * v1 * wvn  !   (wvn = 2.d0*pi/wlc)

           !
           ! integrate single-scattering parameter over size (output)
           !
           if(l <= nang) then
              ph(l,1) = ph(l,1) + z              ! volume phase function (output) [cm-1  str-1]
!              if(ir <= ix0) phs(l) = phs(l) + z  ! for nons
           else if(l == nang+1) then
              cext = cext + z     ! volume extinction cross section (output) [cm2/cm3] or [cm-1]

!toshii
!print*,ir8, 'cext_tot',cext, 'Cext [cm2/cm3]', z, 'Qext [-]=', z/across, across


!              dext(ir) = z        ! for nons
           else if(l==nang+2) then
              cabst = cabst+z   ! volume absorption cross section (output) [cm2/cm3] or [cm-1]
!              dabss(ir) = z     ! for nons
           endif
        enddo

     endif rayleigh

!check output
!write(*,FMT='(I3,1x,F9.3,1x, F9.3, 2x, E9.3, 2x, E9.3, 2x, E9.3, 2x, E9.3, 2x, E9.3)') &
!ir8,real(r_center*1.d4), real(r_center*0.5d0*wvn), cext, cabst, ph(1,1), ph(nang,1), v1

  enddo size_loop

!
! semi-empirical parameterization for non-spherical particles
!
!  if(indr == 1) return   !rayleigh regime
!  call nons(intvl,ix0,nang,nav,cg,cr3,ci3,dcgs,dcgl,gg,rmax,rmin, &
!       rp,wvn,x0,ang,dabss,dcg,dext,phs,szp,ph,cext,cabst)

!toshii
!stop

  return
 end subroutine gtph4b

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine nons(intvl,ix0,nang,nav,cg,cr1,ci1,dcgs,dcgl,gg,rmax,rmin, &
     rp,wvn,x0,ang,dabss,dcg,dext,phs,szp,ph,cext,cabst)
! semi empirical theory to estimate optical properties of non-spherical particles
  !use paras, only: knang,kintvl,kpol,pi,rad
  implicit none
! input
  integer,intent(in):: intvl
  integer,intent(in):: ix0
  integer,intent(in):: nang
  integer,intent(in):: nav
  real(8),intent(in):: cg
  real(8),intent(in):: cr1
  real(8),intent(in):: ci1
  real(8),intent(in):: dcgs
  real(8),intent(in):: dcgl
  real(8),intent(in):: gg
  real(8),intent(in):: rmax
  real(8),intent(in):: rmin
  real(8),intent(in):: rp
  real(8),intent(in):: wvn
  real(8),intent(in):: x0
  real(8),intent(in):: ang(knang)
  real(8),intent(in):: dabss(kintvl)
  real(8),intent(in):: dcg(kintvl)
  real(8),intent(in):: dext(kintvl)
  real(8),intent(inout):: phs(knang)
  real(8),intent(in):: szp(kintvl)
  real(8),intent(out):: ph(knang,kpol)
  real(8),intent(out):: cext
  real(8),intent(out):: cabst

! work
  integer:: imn,imx,ia,ir,jj,i
  real(8):: dcb
  real(8):: xmin,xmax

  real(8):: dm2
  real(8):: didd,az,sn,sn2,cs,cs2
  real(8):: xb,xbz
  real(8):: dx,dx1,dy,dy1,dj,dd,sq
  real(8):: did(knang),dir(knang),dit(knang)

  real(8):: phintg
  real(8):: zs,zd,zr,zt
  real(8):: qss,qsa,qls,qla
  real(8):: fls,qs,qa,qe,qd,qr,qt
!--exec
! semi empirical theory
  if(ix0 > intvl) return
  dcb=0.d0
  xmin=rmin*wvn
  xmax=rmax*wvn
  do imn=1,kintvl
     if(szp(imn) >= xmin) exit
  enddo

! 09.06.15 Bug fixed (Miho, Reported by Toshi Matsui)
!  do imx=intvl+1,1,-1
  do imx=intvl,1,-1
     if(szp(imx) <= xmax) exit
  enddo

  if(imn/=1) imn=imn-1
  if(imx <= ix0) return

  dm2=cr1**2+ci1**2
  do ia=1,nang
     didd=0.d0
     az=ang(ia)
     if(az < 0.018d0) az=0.018d0
     if(az > 179.982d0) az=179.982d0
     sn=sin(az*rad)
     sn2=sin(az*rad*0.5d0)
     cs=cos(az*rad)
     cs2=cos(az*rad*0.5d0)
     do ir=ix0,imx
        xb=sqrt(rp)*szp(ir)
        xbz=xb*sn
        dj=0.d0
        do jj=1,nang-1
           dx1=ang(jj+1)*rad
           dx=ang(jj)*rad
           dy1=cos(dx1-xbz*sin(dx1))
           dy=cos(dx-xbz*sin(dx))
           dj=dj+(dy1+dy)*(dx1-dx)*0.5d0
        enddo
        dd=rp*xb*xb*dj*dj*0.5d0*(1.d0+cs**2)/(pi*xbz*xbz)*dcg(ir)
        didd=didd+dd
     enddo
     if(ang(ia) <= 30.d0) then
        did(ia)=didd
        dcb=didd
     else
        did(ia)=0.d0
     endif
     sq=sqrt(dm2-1.d0+sn2**2)
     dir(ia)=1.d0-2.d0*sn2*sq*(1.d0/(sn2+sq)**2+dm2/(dm2*sn2+sq)**2)
     dit(ia)=exp(1.d0-log(gg)/log(1.53171d0*pi)*ang(ia)*rad)
  enddo

  do i=1,nang
     if(did(i) > dcb) then
        did(i)=did(i)-dcb
     else
        did(i)=0.d0
     endif
  enddo

  zs=phintg(nang,ang,phs)
  zd=phintg(nang,ang,did) !/4.0/pi
  zr=phintg(nang,ang,dir) !/4.0/pi
  zt=phintg(nang,ang,dit) !/4.0/pi

  phs(1:nang)=phs(1:nang)/zs

!  if(zd == 0.d0) then
!    did(1:nang)= 0.d0
!  else
  did(1:nang)=did(1:nang)/zd
!  endif
 
  dir(1:nang)=dir(1:nang)/zr
  dit(1:nang)=dit(1:nang)/zt

  qss=0.d0
  qsa=0.d0
  qls=0.d0
  qla=0.d0
  do ir=imn,imx
     if(ir <= ix0)then
        qss=qss+dext(ir)-dabss(ir)
        qsa=qsa+dabss(ir)
     else
        qls=qls+dext(ir)-dabss(ir)
        qla=qla+dabss(ir)
     endif
  enddo
  qss=qss/dcgs
  qsa=qsa/dcgs
  qls=qls/dcgl
  qla=qla/dcgl
  fls=dcgl/cg
  qs=qls*fls*rp+(1-fls)*qss
  qa=qla*fls+(1-fls)*qsa
  qe=qs+qa
  qd=1.d0
  qr=zr /4.0/pi
  qt=qls-qd-qr
  ph(1:nang,1)=phs(1:nang)*(1-fls)*qss+rp*fls &
       *(did(1:nang)*qd+dir(1:nang)*qr+dit(1:nang)*qt)
  ph(1:nang,1)=ph(1:nang,1)*cg
  cext=qe*cg
  cabst=qa*cg

  return
end subroutine nons

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine smllop(x,cr,ci,qext,qsca)
! optical cross sections for small particles (x<0.1)
! in this approximation
!     q(ang,1)=3.0/8.0/pi*qsca
!     q(ang,2)=q(ang,1)*cos(ang)**2
!     q(ang,3)=q(ang,1)*cos(ang)
!     q(ang,4)=0
!--- history
! 95. 9.19 created
! 08. 9.25 modified Fortran 90 free-style form
!--- input
! x       r       size parameter (<1.5)
! cr      r       m= cr - i ci
! ci      r       ci>0
!--- output
! qext    r       extinction efficiency factor   = Cext / (pai*r*r)
! qsca    r       scattering efficiency factor   = Csca / (pai*r*r)
!---
  implicit none

! input & output
  real(8),intent(in):: x         !! size parameter
  real(8),intent(in):: cr        !! refractive index (real part)
  real(8),intent(in):: ci        !! refractive index (imaginary part)
  real(8),intent(out):: qext     !! extinction efficiency factor
  real(8),intent(out):: qsca     !! scattering efficiency factor

! work
  real(8):: a,b,g,z1,z2,e1,e3,e4,s4
!--
  a=(cr**2+ci**2)**2
  b=cr**2-ci**2
  g=cr*ci
  z1=a+4.d0*b+4.d0
  z2=4.d0*a+12.d0*b+9.d0
  e1=24.d0*g/z1
  e3=g*(4.d0/15.d0+20.d0/3.d0/z2+24.d0/5.d0/z1**2*(7.d0*a+4.d0*b-20.d0))
  e4=8.d0/3.d0/z1**2*((a+b-2.d0)**2-36.d0*g**2)
  s4=8.d0/3.d0/z1**2*((a+b-2.d0)**2+36.d0*g**2)
  qext=e1*x+e3*x**3+e4*x**4
  qsca=s4*x**4
  return
end subroutine smllop



!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine fsoli(wl1,wl2,fsol)
! wavelength integrated solar incidence.
!---history
! 94. 5.30 created by M.Tsukamoto
! 97. 4.19 tuned parameters by T. Nakajima
! 03. 9. 9 Changed for free source form (Fortran 90) by M. Sekiguchi 
!---input
! wl1: min. wavelength of the band (micron)
! wl2: max. wavelength of the band (micron)
!---output
! fsol: integrated solar irradiance (w/m2)
!---
  implicit none
  real(8),intent(in)::wl1,wl2
  real(8),intent(out)::fsol
  integer,parameter::knx=3000,knm=10
  real(8),parameter::dwwlc=1.d6
  integer::n,i
  real(8)::wll1,wll2,wl,dwll,wn,sunir
!--exec
  wll1=log(wl1); wll2=log(wl2)
  n=min(int((wll2-wll1)*dwwlc),knx); n=max(n,knm)
  dwll=(wll2-wll1)/dble(n)
  fsol=0.d0
  do i=1,n
     wl=exp(wll1+(dble(i)-0.5d0)*dwll); wn=1.d4/wl
     fsol=fsol+wl*sunir(wn)*dwll
  enddo
  return
end subroutine fsoli


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine total_opt_sbm_visir(cextp_sbm,cabsp_sbm,php_sbm,dzcm,nang, &
                                 taua,omga,phsfa ) 
 implicit none

 type (particle_sbm_r8), intent(in) :: & 
     cextp_sbm,  & ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_sbm,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_sbm(knang)! volume scattering phase function [cm-1  str-1]

 real(8),intent(in) :: dzcm  !layer height [cm]
 integer,intent(in) :: nang

 real(8),intent(out) :: taua !total particle optical depths [-]
 real(8),intent(out) :: omga !total particle single scattering albedo
 real(8),intent(out) :: phsfa(knang) !total phase function ( phsf = ph * tau * omg ) [str-1]
  
!
! total particle optical depths from extinction cross section [cm2/cm3] or [cm-1]
!
             taua    = (&
                          cextp_sbm%liq    & ! liquid (cloud+rain) 
                        + cextp_sbm%ice_col& ! ice column
                        + cextp_sbm%ice_pla& ! ice plate
                        + cextp_sbm%ice_den& ! ice dendride
                        + cextp_sbm%snow   & ! snow 
                        + cextp_sbm%graupel& ! graupel
                        + cextp_sbm%hail   & ! hail
                       ) * dzcm
!
! total particle extinction for scattering components
! Exact formulae is  omga = sum [ ( (cextp-cabsp) / cextp ) * (cextp*dzcm) ]  / taua
!                                    ---------------------     -----------
!                                          omega                  tau
!
             if(taua <= 0.d0) then
                omga = 0.d0
             else
                omga = ( &
                          (cextp_sbm%liq    - cabsp_sbm%liq    ) & ! liquid (cloud+rain)
                        + (cextp_sbm%ice_col- cabsp_sbm%ice_col) & ! ice column
                        + (cextp_sbm%ice_pla- cabsp_sbm%ice_pla) & ! ice plate
                        + (cextp_sbm%ice_den- cabsp_sbm%ice_den) & ! ice dendride
                        + (cextp_sbm%snow   - cabsp_sbm%snow   ) & ! snow 
                        + (cextp_sbm%graupel- cabsp_sbm%graupel) & ! graupel
                        + (cextp_sbm%hail   - cabsp_sbm%hail   ) & ! hail
                        ) * dzcm / taua
             endif

!
! total particle scattering phase function (here defined as Ph = ph * tau * omg) [str-1]
! Exact formulae is phsfa = sum [ ( php/(cextp-cabsp) ) * ( (cextp-cabsp)/cextp ) * (cextp*dzcm)  ] 
!                                   ------------------     ---------------------     ----------    
!                                     phase function                omega               tau       
!
             if(omga <= 0. ) then
                phsfa(1:knang) = 0.d0 
             else
                phsfa(1:knang)   = ( &
                                   php_sbm(1:knang)%liq    & ! liquid (cloud+rain)
                                 + php_sbm(1:knang)%ice_col& ! ice column
                                 + php_sbm(1:knang)%ice_pla& ! ice plate
                                 + php_sbm(1:knang)%ice_den& ! ice dendride
                                 + php_sbm(1:knang)%snow   & ! snow 
                                 + php_sbm(1:knang)%graupel& ! graupel
                                 + php_sbm(1:knang)%hail   & ! hail
                                  ) * dzcm

             endif


 return
 end subroutine total_opt_sbm_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine total_opt_gen_visir(cextp_gen,cabsp_gen,php_gen,dzcm,nang, &
                                 taua,omga,phsfa ) 
 implicit none

 type (particle_gen_r8), intent(in) :: & 
     cextp_gen,  & ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_gen,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_gen(knang)! volume scattering phase function [cm-1  str-1]

 real(8),intent(in) :: dzcm  !layer height [cm]
 integer,intent(in) :: nang

 real(8),intent(out) :: taua !total particle optical depths [-]
 real(8),intent(out) :: omga !total particle single scattering albedo
 real(8),intent(out) :: phsfa(knang) !total phase function ( phsf = ph * tau * omg ) [str-1]
  
!
! total particle optical depths from extinction cross section [cm2/cm3] or [cm-1]
!
             taua    = (&
                          cextp_gen%cloud  & ! cloud water
                        + cextp_gen%rain   & ! rain 
                        + cextp_gen%ice    & ! cloud ice 
                        + cextp_gen%snow   & ! snow 
                        + cextp_gen%graupel& ! graupel
                        + cextp_gen%hail   & ! hail
                       ) * dzcm
!
! total particle extinction for scattering components
! Exact formulae is  omga = sum [ ( (cextp-cabsp) / cextp ) * (cextp*dzcm) ]  / taua
!                                    ---------------------     -----------
!                                          omega                  tau
!
             if(taua <= 0.d0) then
                omga = 0.d0
             else
                omga = ( &
                          (cextp_gen%cloud  - cabsp_gen%cloud  ) & ! cloud water
                        + (cextp_gen%rain   - cabsp_gen%rain   ) & ! rain 
                        + (cextp_gen%ice    - cabsp_gen%ice    ) & ! cloud ice 
                        + (cextp_gen%snow   - cabsp_gen%snow   ) & ! snow 
                        + (cextp_gen%graupel- cabsp_gen%graupel) & ! graupel
                        + (cextp_gen%hail   - cabsp_gen%hail   ) & ! hail
                        ) * dzcm / taua
             endif

!
! total particle scattering phase function (here defined as Ph = ph * tau * omg ) [str-1]
! Exact formulae is phsfa = sum [ ( php/(cextp-cabsp) ) * ( (cextp-cabsp)/cextp ) * (cextp*dzcm)  ] 
!                                   ------------------     ---------------------     ----------    
!                                     phase function                omega               tau       
!
             if(omga <= 0. ) then
                phsfa(1:knang) = 0.d0 
             else
                phsfa(1:knang)   = ( &
                                   php_gen(1:knang)%cloud  & ! cloud water
                                 + php_gen(1:knang)%rain   & ! rain
                                 + php_gen(1:knang)%ice    & ! cloud ice
                                 + php_gen(1:knang)%snow   & ! snow
                                 + php_gen(1:knang)%graupel& ! graupel
                                 + php_gen(1:knang)%hail   & ! hail
                                  ) * dzcm

             endif


 return
 end subroutine total_opt_gen_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine total_opt_rams_visir(cextp_rams,cabsp_rams,php_rams,dzcm,nang, &
                                 taua,omga,phsfa ) 
 implicit none

 type (particle_rams_r8), intent(in) :: & 
     cextp_rams,  & ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_rams,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_rams(knang)! volume scattering phase function [cm-1  str-1]

 real(8),intent(in) :: dzcm  !layer height [cm]
 integer,intent(in) :: nang

 real(8),intent(out) :: taua !total particle optical depths [-]
 real(8),intent(out) :: omga !total particle single scattering albedo
 real(8),intent(out) :: phsfa(knang) !total phase function ( phsf = ph * tau * omg ) [str-1]
  
!
! total particle optical depths from extinction cross section [cm2/cm3] or [cm-1]
!
             taua    = (&
                          cextp_rams%cloud1 & ! cloud mode1 (small)
                        + cextp_rams%cloud2 & ! cloud mode2 (large)
                        + cextp_rams%rain   & ! rain 
                        + cextp_rams%ice1   & ! ice mode1 (small)
                        + cextp_rams%ice2   & ! ice mode2 (large)
                        + cextp_rams%snow   & ! snow 
                        + cextp_rams%graupel& ! graupel
                        + cextp_rams%hail   & ! hail
                       ) * dzcm
!
! total particle extinction for scattering components
! Exact formulae is  omga = sum [ ( (cextp-cabsp) / cextp ) * (cextp*dzcm) ]  / taua
!                                    ---------------------     -----------
!                                          omega                  tau
!
             if(taua <= 0.d0) then
                omga = 0.d0
             else
                omga = ( &
                          (cextp_rams%cloud1 - cabsp_rams%cloud1 ) & ! cloud mode1 (small)
                        + (cextp_rams%cloud2 - cabsp_rams%cloud2 ) & ! cloud mode2 (large)
                        + (cextp_rams%rain   - cabsp_rams%rain   ) & ! rain 
                        + (cextp_rams%ice1   - cabsp_rams%ice1   ) & ! ice mode1 (small)
                        + (cextp_rams%ice2   - cabsp_rams%ice2   ) & ! ice mode2 (large)
                        + (cextp_rams%snow   - cabsp_rams%snow   ) & ! snow 
                        + (cextp_rams%graupel- cabsp_rams%graupel) & ! graupel
                        + (cextp_rams%hail   - cabsp_rams%hail   ) & ! hail
                        ) * dzcm / taua
             endif

!
! total particle scattering phase function  (here defined as Ph = ph * tau * omg ) [str-1]
! Exact formulae is phsfa = sum [ ( php/(cextp-cabsp) ) * ( (cextp-cabsp)/cextp ) * (cextp*dzcm)  ] 
!                                   ------------------     ---------------------     ----------    
!                                     phase function                omega               tau       
!
             if(omga <= 0. ) then
                phsfa(1:knang) = 0.d0 
             else
                phsfa(1:knang)   = ( &
                                   php_rams(1:knang)%cloud1 & ! cloud mode1 (small)
                                 + php_rams(1:knang)%cloud2 & ! cloud mode2 (large)
                                 + php_rams(1:knang)%rain   & ! rain 
                                 + php_rams(1:knang)%ice1   & ! ice mode1 (small)
                                 + php_rams(1:knang)%ice2   & ! ice mode2 (large)
                                 + php_rams(1:knang)%snow   & ! snow 
                                 + php_rams(1:knang)%graupel& ! graupel
                                 + php_rams(1:knang)%hail   & ! hail
                                  ) * dzcm

             endif


 return
 end subroutine total_opt_rams_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine total_opt_gocart_visir(cextp_gocart,cabsp_gocart,php_gocart,dzcm,nang, &
                                 taua,omga,phsfa )
 implicit none

 type (particle_gocart_r8), intent(in) :: &
     cextp_gocart,  & ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_gocart,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_gocart(knang)! volume scattering phase function [cm-1]

 real(8),intent(in) :: dzcm  !layer height [cm]
 integer,intent(in) :: nang

 real(8),intent(out) :: taua !total particle optical depths [-]
 real(8),intent(out) :: omga !total particle single scattering albedo
 real(8),intent(out) :: phsfa(knang) !total phase function ( phsf = ph * tau * omg ) [str-1]


!
! total particle optical depths from extinction cross section [cm2/cm3] or [cm-1]
!
             taua   = (   &
                          cextp_gocart%so4 & !sulfate and its precurser
                        + cextp_gocart%blc & !black carbon (soot)
                        + cextp_gocart%ocn & !non hygroscopic OC
                        + cextp_gocart%och & !hygroscopic OC
                        + cextp_gocart%ssa & !sea-salt accumulation mode
                        + cextp_gocart%ssc & !sea-salt coarse mode
                        + cextp_gocart%du1 & !dust mode 1
                        + cextp_gocart%du2 & !dust mode 2
                        + cextp_gocart%du3 & !dust mode 3
                        + cextp_gocart%du4 & !dust mode 4
                        + cextp_gocart%du5 & !dust mode 5
                        + cextp_gocart%du6 & !dust mode 6
                        + cextp_gocart%du7 & !dust mode 7
                        + cextp_gocart%du8 & !dust mode 8
                       ) * dzcm
!
! total particle extinction for scattering components
! Exact formulae is  omga = sum [ ( (cextp-cabsp) / cextp ) * (cextp*dzcm) ]  / taua
!                                    ---------------------     -----------
!                                          omega                  tau
!
             if(taua <= 0.) then
                omga = 0.d0
             else
                omga = (  & ! 
                             (cextp_gocart%so4-cabsp_gocart%so4) & !sulfate and its precurser
                           + (cextp_gocart%blc-cabsp_gocart%blc) & !black carbon (soot)
                           + (cextp_gocart%ocn-cabsp_gocart%ocn) & !non hygroscopic OC
                           + (cextp_gocart%och-cabsp_gocart%och) & !hygroscopic OC
                           + (cextp_gocart%ssa-cabsp_gocart%ssa) & !sea-salt accumulation mode
                           + (cextp_gocart%ssc-cabsp_gocart%ssc) & !sea-salt coarse mode
                           + (cextp_gocart%du1-cabsp_gocart%du1) & !dust mode 1
                           + (cextp_gocart%du2-cabsp_gocart%du2) & !dust mode 2
                           + (cextp_gocart%du3-cabsp_gocart%du3) & !dust mode 3
                           + (cextp_gocart%du4-cabsp_gocart%du4) & !dust mode 4
                           + (cextp_gocart%du5-cabsp_gocart%du5) & !dust mode 5
                           + (cextp_gocart%du6-cabsp_gocart%du6) & !dust mode 6
                           + (cextp_gocart%du7-cabsp_gocart%du7) & !dust mode 7
                           + (cextp_gocart%du8-cabsp_gocart%du8) & !dust mode 8
                          ) * dzcm / taua

             endif
!
! total particle scattering phase function  (here defined as Ph = ph * tau * omg ) [str-1]
! Exact formulae is phsfa = sum [ ( php/(cextp-cabsp) ) * ( (cextp-cabsp)/cextp ) * (cextp*dzcm)  ] 
!                                   ------------------     ---------------------     ----------    
!                                     phase function                omega               tau       
!
             if(omga <= 0. ) then
                phsfa(1:knang) = 0.d0 
             else
                phsfa(1:knang) = (  php_gocart(1:knang)%so4 & !sulfate and its precurser
                                 +  php_gocart(1:knang)%blc & !black carbon (soot)
                                 +  php_gocart(1:knang)%ocn & !non hygroscopic OC
                                 +  php_gocart(1:knang)%och & !hygroscopic OC
                                 +  php_gocart(1:knang)%ssa & !sea-salt accumulation mode
                                 +  php_gocart(1:knang)%ssc & !sea-salt coarse mode
                                 +  php_gocart(1:knang)%du1 & !dust mode 1
                                 +  php_gocart(1:knang)%du2 & !dust mode 2
                                 +  php_gocart(1:knang)%du3 & !dust mode 3
                                 +  php_gocart(1:knang)%du4 & !dust mode 4
                                 +  php_gocart(1:knang)%du5 & !dust mode 5
                                 +  php_gocart(1:knang)%du6 & !dust mode 6
                                 +  php_gocart(1:knang)%du7 & !dust mode 7
                                 +  php_gocart(1:knang)%du8 & !dust mode 8
                                ) * dzcm
             endif


 return
 end subroutine total_opt_gocart_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine total_opt_ucd_visir(cextp_ucd,cabsp_ucd,php_ucd,dzcm,nang, &
                                 taua,omga,phsfa )
 implicit none

 type (particle_ucd_r8), intent(in) :: &
     cextp_ucd,  & ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_ucd,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_ucd(knang)! volume scattering phase function [cm-1]

 real(8),intent(in) :: dzcm  !layer height [cm]
 integer,intent(in) :: nang

 real(8),intent(out) :: taua !total particle optical depths [-]
 real(8),intent(out) :: omga !total particle single scattering albedo
 real(8),intent(out) :: phsfa(knang) !total phase function ( phsf = ph * tau * omg ) [str-1]


!
! total particle optical depths from extinction cross section [cm2/cm3] or [cm-1]
!
             taua   = (   &
                          cextp_ucd%du1 & !dust 1
                        + cextp_ucd%du2 & !dust 2
                        + cextp_ucd%du3 & !dust 3
                        + cextp_ucd%du4 & !dust 4
                        + cextp_ucd%du5 & !dust 5
                       ) * dzcm
!
! total particle extinction for scattering components
! Exact formulae is  omga = sum [ ( (cextp-cabsp) / cextp ) * (cextp*dzcm) ]  / taua
!                                    ---------------------     -----------
!                                          omega                  tau
!
             if(taua <= 0.) then
                omga = 0.d0
             else
                omga = (  & ! 
                             (cextp_ucd%du1-cabsp_ucd%du1) & !dust 1
                           + (cextp_ucd%du2-cabsp_ucd%du2) & !dust 2
                           + (cextp_ucd%du3-cabsp_ucd%du3) & !dust 3
                           + (cextp_ucd%du4-cabsp_ucd%du4) & !dust 4
                           + (cextp_ucd%du5-cabsp_ucd%du5) & !dust 5
                          ) * dzcm / taua

             endif
!
! total particle scattering phase function  (here defined as Ph = ph * tau * omg ) [str-1]
! Exact formulae is phsfa = sum [ ( php/(cextp-cabsp) ) * ( (cextp-cabsp)/cextp ) * (cextp*dzcm)  ] 
!                                   ------------------     ---------------------     ----------    
!                                     phase function                omega               tau       
!
             if(omga <= 0. ) then
                phsfa(1:knang) = 0.d0 
             else
                phsfa(1:knang) = ( &
                                    php_ucd(1:knang)%du1 & !dust 1
                                 +  php_ucd(1:knang)%du2 & !dust 2
                                 +  php_ucd(1:knang)%du3 & !dust 3
                                 +  php_ucd(1:knang)%du4 & !dust 4
                                 +  php_ucd(1:knang)%du5 & !dust 5
                                ) * dzcm

             endif


 return
 end subroutine total_opt_ucd_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine makelut_visir(sensor)
 implicit none

!--------------------------------------------------------------------------------------------------
! Comments: 
! Construct the look-up tables of Mie coefficients
! for visir radiometer and lidar simulations.
!
! History:
! 05/2009  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 character(len=5),intent(in) :: sensor  !'visir' or 'lidar'
 character(len=30) :: fname_out(100)
 character(len=100):: char_wav
 integer :: max_spc ! maximum particle species
 integer :: ispc  !RAMS particle habit index
 real*4 :: pbck   ! backscatter phase function/(4*pi)
 real*4 :: frequency
 integer :: nspec,i,j,nwav
 integer :: iounit 
 logical :: fexist
 integer :: mxwav ! maximum # of wavelength channels 
 real*4 :: cfmas, pwmas ! alpha_m & beta_m 
 integer :: ns,nt !humidity and temperature index

 real(sdsu_fps) :: rhoe  !snow density [kg/m3]

 integer,save:: nang                 !! number of scaattering angles
 real(8),save:: ang(knang)           !! scattering angles [degree]
 real(8):: cextp_fake  !extinction cross section [cm2/cm3] or [cm-1]
 real(8):: cabsp_fake !absorption cross section [cm2/cm3] or [cm-1]
 real(8),dimension(knang) :: php_fake   !! volume scattering phase function

 real(8):: scr                  !! relative water reflactive index (real part)
 real(8):: sci                  !! relative water reflactive index (imaginary part)

 real(8),parameter :: q_unit = 1.0d0  !unit mass [g/m3]
 real(8) :: &
     wlc      ,  & !wavelength [cm]
     rhl      ,  & ! RH [0~1]
     cextp_lut,  & ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_lut,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_lut(knang)! volume scattering phase function [cm-1 str-1]

 integer,parameter :: iuk = 120
 integer :: intvl, ipol
  real(8):: szp(kintvl)          !! size parameter
 logical :: logic_opened
 integer :: ierr
 character(len=20) :: sensor_name

!
! NetCDF parameters
!
 integer :: ncid, id
 integer :: wav_dimid, ref_dimid, spc_dimid, nang_dimid
 integer, parameter :: ndims3d = 3, ndims4d = 4
 integer :: dimids3d(ndims3d), dimids4d(ndims4d)
 integer :: id_cextp, id_cabsp, id_php

!
! initialization of fundamental particle refractive index (rfi)
! outputs are stored in the saved global parameter in module paras
!
  call read_rfi


 if(sensor == 'lidar') then
    mxwav = mxwavel_lidar
    sensor_name = trim(lidar_sensor)
   write(*,*) '- MSG makelut_visir: LUT for lidar simulation - '
 elseif(sensor == 'visir') then
    mxwav = mxwavel
    sensor_name = trim(visir_sensor)
    write(*,*) '- MSG makelut_visir: LUT for visir simulation - '
 else
    call stop_sdsu(myrank, 'MSG makelut_visir: There is no such sensor')
 endif


 select case(trim(type_microphysics))
   case('GEN')    ! GCE (LIN_type) 1moment bulk microphsics 
    max_spc = 6
   case('RAMS')  !GCE RAMS 1- and 2-moment microphysics 
    max_spc = 16
   case('SBM')
    call stop_sdsu(myrank,'Not yet supporoting SBM')
   case default
 end select


!
! allocate array parameters
!
 if( .not. allocated(cextp_lut_array) ) then

   allocate( cextp_lut_array (          mxpts_re_visir, max_spc, max_chan ), &
             cabsp_lut_array (          mxpts_re_visir, max_spc, max_chan ), &
             php_lut_array   ( 1:knang, mxpts_re_visir, max_spc, max_chan ), &
             stat=ierr )
   cextp_lut_array = 0.  ; cabsp_lut_array = 0. ; php_lut_array = 0.

   if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_visir: allocation error -> Terminate program.')

 endif



 !
 ! open netcdf file
 !
 sdsu_io_file = trim(sdsu_dir_sslut)//'VISIR_CLOUD_'//trim(cloud_microphysics)//'-'//trim(sensor_name)//'.nc'

 inquire(file=trim(sdsu_io_file),exist=fexist)
 if( lut_replace ) fexist = .false.
 if (fexist) then !skip
     print*,'MSG makelut_visir; ',trim(sdsu_io_file),' exist.'
     return
 endif

 print*,'MSG makelut_visir; Generating VisibleIR LUT ',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )               !32bit output

 !
 ! Define netcdf dimensions
 !
 call check( nf90_def_dim(ncid, 'wav_max'  , max_chan         , wav_dimid) )
 call check( nf90_def_dim(ncid, 're_max'   , mxpts_re_visir   , ref_dimid) )
 call check( nf90_def_dim(ncid, 'spc_max'  , max_spc          , spc_dimid) )
 call check( nf90_def_dim(ncid, 'nang_max' , knang            , nang_dimid) )

 dimids3d  = (/             ref_dimid, spc_dimid, wav_dimid /)
 dimids4d  = (/ nang_dimid, ref_dimid, spc_dimid, wav_dimid /)

!
! global parameter
!
 if(sensor == 'lidar') then
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'wavelength[micron]', wavel_lidar(1:mxwav)  ) )
 elseif(sensor == 'visir') then
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'wavelength[micron]', wavel(1:mxwav)  ) )
 endif


 !
 ! output parameter definition
 !
 call check( nf90_def_var(ncid, 'cextp', NF90_DOUBLE, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'description',' extinction cross section ') )
 call check( nf90_put_att(ncid, id, 'units', '[cm2/cm3] or [cm-1]') )
 id_cextp = id

 call check( nf90_def_var(ncid, 'cabsp', NF90_DOUBLE, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'description',' absorption cross section ') )
 call check( nf90_put_att(ncid, id, 'units', '[cm2/cm3] or [cm-1]') )
 id_cabsp = id

 call check( nf90_def_var(ncid, 'php', NF90_DOUBLE, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' volume scattering phase function ') )
 call check( nf90_put_att(ncid, id, 'units', '[cm-1 str-1]') )
 id_php = id

 call check( nf90_enddef(ncid) )  ! End define mode.

 print*,''

  WAV_LOOP : do nwav = 1 , mxwav

    if(sensor == 'visir') then
       wlc = DBLE(wavel(nwav))*1.0d-4  ![cm]
       write(*,*) 'Wavelength = ',wavel(nwav),' micron'
       char_wav = nch_wavel(nwav)
    elseif(sensor == 'lidar') then
       wlc = DBLE(wavel_lidar(nwav))*1.0d-4 ![cm]
       write(*,*) 'Wavelength = ',wavel_lidar(nwav),' micron'
       char_wav = nch_wavel_lidar(nwav)
    else
       call stop_sdsu(myrank,'MSG makelut_visir: There is no such sensor')
    endif

!
! species loop 
!
 SPC_LOOP: do nspec = 1, max_spc

!
!       Calculate scattering parameters per unit mass mixing ratio 1 [g/cm3]
!
    RE_LOOP: do i = 1, mxpts_re_visir
       cextp_lut=0.d0 ;cabsp_lut=0.d0 ;php_lut=0.d0
       rhl = 0.1d0 ! RH (0~1) not used for condensates (but aerosols)

       mic_select2: select case(trim(type_microphysics))
       case('GEN')      ! general bulk microphsics 

           NSPEC_GCE: if(nspec==1) then  !cloud use gamma PSD in RAMS (cloud1 class)

!             if( trim(cloud_microphysics) == 'THOM') then

!                call mie_gen_visir('qc ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%cloud), &
!                                   dble(rho_gen%cloud), wlc,rhl, & !rain
!                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )

!             else

              ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud
              call mie_rams_visir('qc1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%cloud1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

!             endif

            elseif(nspec==2) then           
                call mie_gen_visir('qr ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%rain), &
                                   dble(rho_gen%rain), wlc,rhl, & !rain
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )

            elseif(nspec==3) then !ice use gamma PSD in RAMS (ice column)

              if( trim(cloud_microphysics) == 'MORR' .or. trim(cloud_microphysics) == 'MORRH' .or.&
                  trim(cloud_microphysics) == 'THOM' ) then !  twomoment

                call mie_gen_visir('qi ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%ice), &
                                   dble(rho_gen%ice), wlc,rhl, & !rain
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )


              else

               ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice
               call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , & 
                    dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl        , &
                    nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

              endif


            elseif(nspec==4) then


              if( trim(cloud_microphysics) == 'GMP4ICE' ) then
               rhoe = 1.e+3 * max(0.05, min(0.9e0, 142.36 * ( pts_re_visir(i)**(-0.95) ) ) ) ! re-dependent bulk snow density [kg/m3]
              elseif( trim(cloud_microphysics) == 'THOM' ) then
               rhoe = min(900. , 0.13 / (pts_re(j)*2.*1.e-6) ) ! Kg/m^3
              else
                rhoe = rho_gen%snow
              endif

                call mie_gen_visir('qs ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%snow), &
                                   dble(rhoe), wlc,rhl, & !snow
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )


            elseif(nspec==5) then
                call mie_gen_visir('qg ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%graupel), &
                                    dble(rho_gen%graupel), wlc,rhl, & !graupel
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )


            elseif(nspec==6) then
                call mie_gen_visir('qh ', q_unit, dble(pts_re_visir(i)), dble(mu_gen%hail), &
                                   dble(rho_gen%hail), wlc,rhl, & !hail
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )
!toshii: check
!print*, nspec, pts_re_visir(i), cextp_lut, cabsp_lut 
!print*, nspec, pts_re_visir(i), cextp_lut / php_lut(nang)

            else
               call stop_sdsu(myrank, 'MSG makelut_visir: there is no such nspec')
            endif NSPEC_GCE
          

       case('RAMS')  !GCE RAMS 1- & 2-moment microphysics 

            NSPEC_RAMS: if(nspec==1) then
              ispc = 1 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud1
              call mie_rams_visir('qc1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%cloud1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==2) then
              ispc = 16; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! cloud2
              call mie_rams_visir('qc2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%cloud2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==3) then
              ispc = 2 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! rain 
              call mie_rams_visir('qr ',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%rain), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==4) then
              ispc = 3 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 col
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==5) then
              ispc = 8 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 hex
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==6) then
              ispc = 9 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 den
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==7) then
              ispc = 10 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 ndl
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==8) then
              ispc = 11 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice1 ros
              call mie_rams_visir('qi1',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice1), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==9) then
              ispc = 4 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 col
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==10) then
              ispc = 12 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 hex
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==11) then
              ispc = 13 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 den
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==12) then
              ispc = 14 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 ndl
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==13) then
              ispc = 15 ;cfmas = rams_dstprms(2,ispc) ;pwmas = rams_dstprms(3,ispc)  ! ice2 ros
              call mie_rams_visir('qi2',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%ice2), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==14) then
              ispc = 5 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! snow aggregate
              call mie_rams_visir('qs ',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%snow), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==15) then
              ispc = 6 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! graupel
              call mie_rams_visir('qg ',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%graupel), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            elseif(nspec==16) then
               ispc = 7 ; cfmas = rams_dstprms(2,ispc)  ; pwmas = rams_dstprms(3,ispc)  ! hail
              call mie_rams_visir('qh ',dble(cfmas),dble(pwmas), q_unit*1.0d-6 , &
                   dble(pts_re_visir(i)), dble(gnu_rams%hail), wlc, rhl      , &
                   nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut           )

            else
               call stop_sdsu(myrank, 'MSG makelut_visir: there is no such nspec')
            endif NSPEC_RAMS

       case('SBM')

         call stop_sdsu(myrank, 'MSG makelut_visir: do not use lut option for SBM')

       case default

         ! do nothing 

       end select mic_select2

     
      ! if nan happen, usually, wavelength is out of ranges.
      if( isnan(cextp_lut) ) then
         call stop_sdsu(myrank,'MSG makelut_visir: single scattering properties are NaN')
      endif 

!       write(iounit,90) cextp_lut, cabsp_lut, php_lut(1:knang)

      ! transfering to array data
       cextp_lut_array(       i,nspec,nwav) = cextp_lut
       cabsp_lut_array(       i,nspec,nwav) = cabsp_lut
       php_lut_array(1:knang, i,nspec,nwav) = php_lut(1:knang)

!toshii: check
!print*, nspec, pts_re_visir(i), cextp_lut, cabsp_lut

     enddo RE_LOOP
  enddo SPC_LOOP

  enddo WAV_LOOP

!
! write output in NetCDF
!
 call check( nf90_put_var(ncid, id_cextp , cextp_lut_array ) )
 call check( nf90_put_var(ncid, id_cabsp , cabsp_lut_array ) )
 call check( nf90_put_var(ncid, id_php   , php_lut_array   ) )

 call check( nf90_close(ncid) )  ! close nc file


  90  Format ( 1P82E14.6 )

 print*,''

 return
 end subroutine makelut_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine makelut_visir_aerosol(sensor)
 implicit none

!--------------------------------------------------------------------------------------------------
! Comments: 
! Construct the look-up tables of Mie coefficients
! for visir radiometer and lidar simulations.
!
! History:
! 05/2009  Toshi Matsui@NASA GSFC ; Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 character(len=5),intent(in) :: sensor  !'visir' or 'lidar'
 character(len=30) :: fname_out(100)
 character(len=100):: char_wav
 integer :: max_spc ! maximum particle species
 integer :: ispc  !RAMS particle habit index
 real*4 :: pbck   ! backscatter phase function/(4*pi)
 real*4 :: frequency
 integer :: nspec,i,j,nwav
 integer :: iounit 
 logical :: fexist
 integer :: mxwav ! maximum # of wavelength channels 

 integer,save:: nang                 !! number of scaattering angles
 real(8),save:: ang(knang)           !! scattering angles [degree]
 real(8):: scr                  !! relative water reflactive index (real part)
 real(8):: sci                  !! relative water reflactive index (imaginary part)

 real(8),parameter :: q_unit = 1.0d0  !unit mass [g/m3]
 real(8) :: &
     wlc      ,  & !wavelength [cm]
     rhl      ,  & ! RH [0~1]
     cextp_lut,  & ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_lut,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_lut(knang),&! volume scattering phase function [-]
     cextp_lut_previous,  & ! extinction cross section [cm2/cm3] or [cm-1]
     cabsp_lut_previous,  & ! absorption cross section [cm2/cm3] or [cm-1]
     php_lut_previous(knang)! volume scattering phase function [-]

 integer :: ierr
 character(len=20) :: sensor_name

!
! NetCDF parameters
!
 integer :: ncid, id
 integer :: wav_dimid, rh_dimid, spc_dimid, nang_dimid
 integer, parameter :: ndims3d = 3, ndims4d = 4
 integer :: dimids3d(ndims3d), dimids4d(ndims4d)
 integer :: id_cextp, id_cabsp, id_php

!
! initialization of fundamental particle refractive index (rfi)
! outputs are stored in the saved global parameter in module paras
!
 call read_rfi


 if(sensor == 'lidar') then
    mxwav = mxwavel_lidar
    sensor_name = trim(lidar_sensor)
   write(*,*) '- MSG makelut_visir_aerosol: LUT for lidar simulation - '
 elseif(sensor == 'visir') then
    mxwav = mxwavel
    sensor_name = trim(visir_sensor)
    write(*,*) '- MSG makelut_visir_aerosol: LUT for visir simulation - '
 else
    call stop_sdsu(myrank, 'MSG makelut_visir_aerosol: There is no such sensor')
 endif


 select case(trim(aerosol_microphysics) )
 case('GOCART')
    max_spc = 14
 case('UCD')
    max_spc = 5
 end select

!
! allocate array parameters
!
 if( .not. allocated(cextp_aero_lut_array) ) then

   allocate( cextp_aero_lut_array (          mxpts_rh, max_spc , max_chan ), &
             cabsp_aero_lut_array (          mxpts_rh, max_spc , max_chan ), &
             php_aero_lut_array   ( 1:knang, mxpts_rh, max_spc , max_chan ), &
             stat=ierr )
   cextp_aero_lut_array = 0. ; cabsp_aero_lut_array = 0. ; php_aero_lut_array = 0. 
   if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_visir_aero: allocation error -> Terminate program.')

 endif

!
! open netcdf file
!
 sdsu_io_file = trim(sdsu_dir_sslut)//'VISIR_AEROSOL_'//trim(aerosol_microphysics)//'-'//trim(sensor_name)//'.nc'

 inquire(file=trim(sdsu_io_file),exist=fexist)
 if( lut_replace ) fexist = .false.
 if (fexist) then !skip
    print*,'MSG makelut_visir_aero; ',trim(sdsu_io_file),' exist.'
    return
 endif

 print*,'MSG makelut_visir_aero; Generating VisibleIR LUT ',trim(sdsu_io_file)

 call check( nf90_create(trim(sdsu_io_file), nf90_clobber, ncid) )               !32bit output

!
! Define netcdf dimensions
!
 call check( nf90_def_dim(ncid, 'wav_max'  , max_chan   , wav_dimid ) )
 call check( nf90_def_dim(ncid, 'rh_max'   , mxpts_rh   , rh_dimid  ) )
 call check( nf90_def_dim(ncid, 'spc_max'  , max_spc    , spc_dimid ) )
 call check( nf90_def_dim(ncid, 'nang_max' , knang      , nang_dimid) )

 dimids3d  = (/             rh_dimid, spc_dimid, wav_dimid /)
 dimids4d  = (/ nang_dimid, rh_dimid, spc_dimid, wav_dimid /)

!
! global parameter
!
 if(sensor == 'lidar') then
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'wavelength[micron]', wavel_lidar(1:mxwav)  ) )
 elseif(sensor == 'visir') then
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'wavelength[micron]', wavel(1:mxwav)  ) )
 endif

 !
 ! output parameter definition
 !
 call check( nf90_def_var(ncid, 'cextp', NF90_DOUBLE, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'description',' extinction cross section ') )
 call check( nf90_put_att(ncid, id, 'units', '[cm2/cm3] or [cm-1]') )
 id_cextp = id

 call check( nf90_def_var(ncid, 'cabsp', NF90_DOUBLE, dimids3d, id) )
 call check( nf90_put_att(ncid, id, 'description',' absorption cross section ') )
 call check( nf90_put_att(ncid, id, 'units', '[cm2/cm3] or [cm-1]') )
 id_cabsp = id

 call check( nf90_def_var(ncid, 'php', NF90_DOUBLE, dimids4d, id) )
 call check( nf90_put_att(ncid, id, 'description',' volume scattering phase function ') )
 call check( nf90_put_att(ncid, id, 'units', '[cm-1 str-1]') )
 id_php = id

 call check( nf90_enddef(ncid) )  ! End define mode.

  
 WAV_LOOP : do nwav = 1 , mxwav

      if(sensor == 'visir') then
          wlc = DBLE(wavel(nwav))*1.0d-4  ![cm]
          write(*,*) 'Wavelength = ',wavel(nwav),' micron'
          char_wav = nch_wavel(nwav)
      elseif(sensor == 'lidar') then
          wlc = DBLE(wavel_lidar(nwav))*1.0d-4 ![cm]
          write(*,*) 'Wavelength = ',wavel_lidar(nwav),' micron'
          char_wav = nch_wavel_lidar(nwav)
      else
         call stop_sdsu(myrank, 'MSG makelut_visir_aerosol: There is no such sensor')
      endif
   print*,''


!
! species loop 
!
 SPC_LOOP: do nspec = 1, max_spc

   sdsu_io_file = trim(sdsu_dir_sslut)//fname_out(nspec)

!
! inquire Mie LUTs exist or not (if not create LUT)
!

!  inquire(file=sdsu_io_file,exist=fexist)
!  if( lut_replace ) fexist = .false.
!  INQ: if (fexist) then !skip
!      print*,'MSG makelut_visir_aerosol; ',trim(fname_out(nspec)),' exists.'
!  else !create lut
!      print*,'MSG makelut_visir_aerosol; Create LUT ',trim(fname_out(nspec))
!      iounit=200
!      open(iounit,file=sdsu_io_file,STATUS='REPLACE')

!
!       Calculate scattering parameters per unit mass mixing ratio 1 [g/cm3]
!
    RH_LOOP: do i = 1, mxpts_rh
       cextp_lut=0.d0 ;cabsp_lut=0.d0 ;php_lut=0.d0

       rhl = DBLE( pts_rh(i) )   ! RH (0~1) for condensates (aerosols)


       select case(trim(aerosol_microphysics) )
       case('GOCART')

            NSPEC_GOCART: if(nspec==1) then
                call mie_gocart('so4', q_unit*1.0d-6, wlc,rhl,         &   !sulfur and its precure
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut )

            elseif(nspec==2) then
                call mie_gocart('blc', q_unit*1.0d-6, wlc,rhl,         &   !black carbon (soot)
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==3) then
                call mie_gocart('ocn', q_unit*1.0d-6, wlc, rhl,         &   !non hygroscopic OC
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==4) then
                call mie_gocart('och', q_unit*1.0d-6, wlc,rhl,         &   !hygroscopic OC
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==5) then
                call mie_gocart('ssa', q_unit*1.0d-6, wlc,rhl,         &   !sea-salt accumulation mode
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==6) then
                call mie_gocart('ssc', q_unit*1.0d-6, wlc,rhl,         &   !sea-salt coarse mode
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==7) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du1', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 1
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==8) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du2', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 2 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==9) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du3', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 3 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==10) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du4', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 4 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==11) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du5', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 5
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==12) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du6', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 6
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==13) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du7', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 7
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==14) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_gocart('du8', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 8
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            else
               call stop_sdsu(myrank, 'MSG makelut_visir_aerosol: there is no such nspec')
            endif NSPEC_GOCART

            !
            ! For RH index greather than 1 and dust species, SS parameters are constatnt (not RH function).
            !

            if(i==1 .and. nspec>=7) then
               cextp_lut_previous=cextp_lut 
               cabsp_lut_previous= cabsp_lut 
               php_lut_previous(1:knang)=php_lut(1:knang) 
            endif

            if(i>1 .and. nspec>=7) then
               cextp_lut=cextp_lut_previous 
               cabsp_lut= cabsp_lut_previous 
               php_lut(1:knang)=php_lut_previous(1:knang)  
            endif

         case('UCD')

            NSPEC_UCD: if(nspec==1) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du1', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 1
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==2) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du2', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 2 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==3) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du3', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 3 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==4) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du4', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 4 
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )

            elseif(nspec==5) then
                if(i==1) &  !there is no RH function in dust, so calculation is just done for RH index = 1
                call mie_ucd('du5', q_unit*1.0d-6, wlc,rhl,         &   !dust mode 5
                                  nang, ang, scr, sci, cextp_lut, cabsp_lut, php_lut  )
            else
               call stop_sdsu(myrank, 'MSG makelut_visir_aerosol: there is no such nspec')
            endif NSPEC_UCD

            !
            ! For RH index greather than 1 and dust species, SS parameters are constatnt (not RH function).
            !

            if( i==1 ) then
               cextp_lut_previous=cextp_lut
               cabsp_lut_previous= cabsp_lut
               php_lut_previous(1:knang)=php_lut(1:knang)

! toshii
! print*,'nspec=',nspec, 'ext=',cextp_lut, 'P11=',php_lut(nang), 'P22=',php_lut(nang+1)

            endif

            if( i>1 ) then
               cextp_lut=cextp_lut_previous
               cabsp_lut= cabsp_lut_previous
               php_lut(1:knang)=php_lut_previous(1:knang)
            endif
         
         end select

      ! if nan happen, usually, wavelength is out of ranges.
      if( isnan(cextp_lut) ) then
         call stop_sdsu(myrank,'MSG makelut_visir: single scattering properties are NaN')
      endif 

       ! transfering to array data
       cextp_aero_lut_array(      i,nspec,nwav) = cextp_lut
       cabsp_aero_lut_array(      i,nspec,nwav) = cabsp_lut
       php_aero_lut_array(1:knang,i,nspec,nwav) = php_lut(1:knang)

!       write(iounit,90) cextp_lut, cabsp_lut, php_lut(1:knang)
     enddo RH_LOOP
  enddo SPC_LOOP

  enddo WAV_LOOP

!
! write output in NetCDF
!
 call check( nf90_put_var(ncid, id_cextp , cextp_aero_lut_array ) )
 call check( nf90_put_var(ncid, id_cabsp , cabsp_aero_lut_array ) )
 call check( nf90_put_var(ncid, id_php   , php_aero_lut_array   ) )

 call check( nf90_close(ncid) )  ! close nc file

  90  Format ( 1P82E14.6 )

 print*,''

 return
 end subroutine makelut_visir_aerosol

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_lut_visir_nc(sensor_name,mxwav)
 implicit none
 character(len=20),intent(in) :: sensor_name
 integer,intent(in) :: mxwav
 integer :: max_spc
 integer :: ierr, ncid, varid

!
! 
!
 select case(trim(type_microphysics))
 case('GEN')    
   max_spc = 6
 case('RAMS')   
   max_spc = 16
 case('SBM')
   call stop_sdsu(myrank,'Not yet supporoting SBM')
 case default
 end select

!
! allocate array parameters
!
 if( .not. allocated(cextp_lut_array) ) then

   allocate( cextp_lut_array ( mxpts_re_visir, max_spc , max_chan ), &
             cabsp_lut_array ( mxpts_re_visir, max_spc , max_chan ), &
             php_lut_array   ( 1:knang, mxpts_re_visir, max_spc, max_chan ), &
             stat=ierr )
   cextp_lut_array = 0.  ; cabsp_lut_array = 0. ; php_lut_array = 0. 
   if (ierr /= 0) call stop_sdsu(myrank,'MSG makelut_visir: allocation error -> Terminate program.')

 endif

!
! open netcdf file
!
 sdsu_io_file = trim(sdsu_dir_sslut)//'VISIR_CLOUD_'//trim(cloud_microphysics)//'-'//trim(sensor_name)//'.nc'
 if(masterproc) print*, 'MSG read_lut_visir_nc: reading visible-IR LUT:',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! extinction cross section [cm2/cm3] or [cm-1]
!
 call check( nf90_inq_varid(ncid, 'cextp' , varid ))
 call check( nf90_get_var(ncid, varid, cextp_lut_array ))

!
! absorption cross section [cm2/cm3] or [cm-1]
!
 call check( nf90_inq_varid(ncid, 'cabsp' , varid ))
 call check( nf90_get_var(ncid, varid, cabsp_lut_array ))

!
! volume scattering phase function [cm-1 str-1]
!
 call check( nf90_inq_varid(ncid, 'php' , varid ))
 call check( nf90_get_var(ncid, varid, php_lut_array ))

!
! close netCDF file
!
  call check( nf90_close(ncid) )

 return
 end subroutine read_lut_visir_nc
 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_lut_visir_aero_nc(sensor_name,mxwav)
 implicit none
 character(len=20),intent(in) :: sensor_name
 integer,intent(in) :: mxwav
 integer :: max_spc
 integer :: ierr, ncid, varid

!
!
!
 select case(trim(aerosol_microphysics) )
 case('GOCART')
   max_spc = 14
 case('UCD')
   max_spc = 5
 end select

!
! allocate array parameters
!
 if( .not. allocated(cextp_aero_lut_array) ) then

   allocate( cextp_aero_lut_array (          mxpts_rh, max_spc, max_chan ), &
             cabsp_aero_lut_array (          mxpts_rh, max_spc, max_chan ), &
             php_aero_lut_array   ( 1:knang, mxpts_rh, max_spc, max_chan ), &
             stat=ierr )
   cextp_aero_lut_array = 0. ; cabsp_aero_lut_array = 0. ; php_aero_lut_array = 0. 

   if (ierr /= 0) call stop_sdsu(myrank,'MSG read_lut_visir_aero_nc: allocation error -> Terminate program.')

 endif

!
! open netcdf file
!
 sdsu_io_file = trim(sdsu_dir_sslut)//'VISIR_AEROSOL_'//trim(aerosol_microphysics)//'-'//trim(sensor_name)//'.nc'
 if(masterproc) print*, 'MSG read_lut_visir_aero_nc: reading visible-IR LUT:',trim(sdsu_io_file)
 call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncid) )

!
! extinction cross section [cm2/cm3] or [cm-1]
!
 call check( nf90_inq_varid(ncid, 'cextp' , varid ))
 call check( nf90_get_var(ncid, varid, cextp_aero_lut_array ))

!
! absorption cross section [cm2/cm3] or [cm-1]
!
 call check( nf90_inq_varid(ncid, 'cabsp' , varid ))
 call check( nf90_get_var(ncid, varid, cabsp_aero_lut_array ))

!
! volume scattering phase function [cm-1 str-1]
!
 call check( nf90_inq_varid(ncid, 'php' , varid ))
 call check( nf90_get_var(ncid, varid, php_aero_lut_array ))

!
! close netCDF file
!
  call check( nf90_close(ncid) )

 return
 end subroutine read_lut_visir_aero_nc
 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mie_lut_visir( q , re, cextp_lut, cabsp_lut, php_lut, pts_ref, mxpts_ref, &
                          cextp_r8, cabsp_r8, php_r8 )
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the extinction/scattering coefficients for
! given Re using Mie LUTs (per unit q), 
! and optal values are re-scaled to actual q (hydrometeor content).
! 
! History:
!   05/2009  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------

 real(4),intent(in) :: q  ! particle mass mixing ratio [g/m3]
 real(4) :: re            ! particle effective radius [micron]
 real(4),intent(in) :: cextp_lut(mxpts_ref)  !extinction cross section [cm2/cm3] or [cm-1]
 real(4),intent(in) :: cabsp_lut(mxpts_ref)  !apsorption cross section [cm2/cm3] or [cm-1]
 real(4),intent(in) :: php_lut(knang,mxpts_ref) !phase function
 real(4),intent(in) :: pts_ref(mxpts_ref)  !LUT of re index [micron]
 integer,intent(in) :: mxpts_ref

 real(8),intent(out) :: cextp_r8      !extinction cross section [cm2/cm3] or [cm-1]
 real(8),intent(out) :: cabsp_r8      !absorption cross section [cm2/cm3] or [cm-1]
 real(8),intent(out) :: php_r8(knang) !phase function

! local parameter
 integer :: j !loop index
 real(4) :: cextp
 real(4) :: cabsp
 real(4) :: php(knang)
 real(4) :: wgt1,wgt2

!
!initialise
!
 cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0

 if(q <= q_min_conden) return


 if(re < pts_ref(1)        ) re = pts_ref(1)
 if(re > pts_ref(mxpts_ref)) re = pts_ref(mxpts_ref)

!
! Linear interpolation
!
 REF_LOOP: do j = 1,mxpts_ref-1

  if( re >= pts_ref(j) .and. re <= pts_ref(j+1) ) then
      wgt2 =  ( re - pts_ref(j) )    / ( pts_ref(j+1) - pts_ref(j) )
      wgt1 = 1.0 - wgt2

      cextp = wgt1*cextp_lut(j) + wgt2*cextp_lut(j+1) 
      cabsp = wgt1*cabsp_lut(j) + wgt2*cabsp_lut(j+1) 
      php(1:knang) = wgt1*php_lut(1:knang,j) + wgt2*php_lut(1:knang,j+1) 
 
      exit REF_LOOP
  endif 

 enddo REF_LOOP


!
! Re-scale from unit mass mixing ratio (1[g/cm3]) to actual mass mixing ratio (q [g/cm3])
! 
  cextp_r8 = dble( cextp * q ) 
  cabsp_r8 = dble( cabsp * q ) 
  php_r8(1:knang) = dble( php(1:knang) * q ) 

 return
 end subroutine mie_lut_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mie_lut_visir_aerosol( q , rh, cextp_lut, cabsp_lut, php_lut, pts_rh, mxpts_rh, &
                                   cextp_r8, cabsp_r8, php_r8 )
 implicit none

!---------------------------------------------------------------------------------------------------
! Comments:  
! This subroutine find the extinction/scattering coefficients for
! given Re using Mie LUTs (per unit q), 
! and optal values are re-scaled to actual q (hydrometeor content).
! 
! History:
!   05/2009  Toshi Matsui@NASA GSFC ; Initial     
!           
! References: 
!---------------------------------------------------------------------------------------------------

 real(4),intent(in) :: q  ! particle mass mixing ratio [g/m3]
 real(4) :: rh            ! relative humidity (0~1)
 real(4),intent(in) :: cextp_lut(mxpts_rh)  !extinction cross section [cm2/cm3] or [cm-1]
 real(4),intent(in) :: cabsp_lut(mxpts_rh)  !apsorption cross section [cm2/cm3] or [cm-1]
 real(4),intent(in) :: php_lut(knang,mxpts_rh) !phase function
 real(4),intent(in) :: pts_rh(mxpts_rh)  !LUT of RH index [micron]
 integer,intent(in) :: mxpts_rh

 real(8),intent(out) :: cextp_r8      !extinction cross section [cm2/cm3] or [cm-1]
 real(8),intent(out) :: cabsp_r8      !absorption cross section [cm2/cm3] or [cm-1]
 real(8),intent(out) :: php_r8(knang) !phase function

! local parameter
 integer :: j !loop index
 real(4) :: cextp
 real(4) :: cabsp
 real(4) :: php(knang)
 real(4) :: wgt1,wgt2

!
!initialise
!
 cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0

 if(q <= q_min_aero) return

 if(rh < pts_rh(1)        ) rh = pts_rh(1)
 if(rh > pts_rh(mxpts_rh) ) rh = pts_rh(mxpts_rh)

!
! Linear interpolation
!
 REF_LOOP: do j = 1,mxpts_rh-1

  if( rh >= pts_rh(j) .and. rh <= pts_rh(j+1) ) then
      wgt2 =  ( rh - pts_rh(j) )    / ( pts_rh(j+1) - pts_rh(j) )
      wgt1 = 1.0 - wgt2

      cextp = wgt1*cextp_lut(j) + wgt2*cextp_lut(j+1) 
      cabsp = wgt1*cabsp_lut(j) + wgt2*cabsp_lut(j+1) 
      php(1:knang) = wgt1*php_lut(1:knang,j) + wgt2*php_lut(1:knang,j+1) 
 
      exit REF_LOOP
  endif 

 enddo REF_LOOP


!
! Re-scale from unit mass mixing ratio (1[g/cm3]) to actual mass mixing ratio (q [g/cm3])
! 
  cextp_r8 = dble( cextp * q ) 
  cabsp_r8 = dble( cabsp * q ) 
  php_r8(1:knang) = dble( php(1:knang) * q ) 

 return
 end subroutine mie_lut_visir_aerosol

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine makelut_visir_grnd
 implicit none

 integer :: nda  !# of stream
 integer :: nw   !# of wavelength
 real(8) :: wl(knw)   !wavelength [micron]
 real(8) :: wlc          ! wavelength [cm]
 real(8) :: amua(kndm)   ! cosine of the zenith-quadrature angles. 1 to 0.
 real(8) :: wa(kndm)     ! corresponding quadrature weights.
 real(8):: cr            ! relative water reflactive index (real part)
 real(8):: ci            ! relative water reflactive index (imaginary part)
 integer :: m1,m         ! Fourier order
 integer :: iw,iu,ic,n1,n2   ! loop indice
 integer :: ierr         ! allocation stat
 real(8) :: r(kndm,kndm)   ! reflection matrix
 real(8) :: sr(kndm)      ! source matrix
 real(8) :: u10          ! 10m wind speed [m/s] (0m/s cause numerical instability)
 real(8) :: cosz         ! cosine of solar zenith angle
 integer,parameter :: iu_lut = 2
 
 real(8) :: du10         ! increment of u10 [m/s]
 character(len=50) :: fname_out(2)
 logical :: fexist


!
! initialization of fundamental particle refractive index (rfi)
! outputs are stored in the saved global parameter in module paras
!
  call read_rfi


 nda = N_stream  !N_stream is module par
 nw = mxwavel    !mxwavel is module par
 wl(1:nw) = DBLE(wavel(1:nw))  !wavel is module para

 du10 = 25.d0 / REAL(mxpts_u10) !mxpts_u10 is module para

!
! get gausian cosz and weights
!
 call  qgausn( wa, amua, nda ) ! nda is input
 amua(1:nda)=amua(nda:1:-1)    ! re-order
 wa(1:nda)=wa(nda:1:-1)        ! re-order

 do iw=1,nw  ! wavelength loop

    wlc = wl(iw) * 1.0d-4              ! wavelength [micron] -> [cm]

    fname_out(1) = 'ORF_SURFACE_'//trim(nch_wavel(iw))//'.dat'  !nch_wavel is module para
    fname_out(2) = 'SRR_SURFACE_'//trim(nch_wavel(iw))//'.dat'

!
! get water reflative index
!
    call water_rhi( wlc, cr, ci )

!
! inquire ORF LUT exist or not (if not create LUT)
!
    sdsu_io_file = trim(sdsu_dir_sslut)//fname_out(1)
    inquire(file=trim(sdsu_io_file),exist=fexist)

    if(lut_grnd_replace) fexist = .false. !for debug  lut_grnd_replace is module para

    if (fexist) then !read lut

      print*,'MSG makelut_visir_grnd; ',trim(fname_out(1)),' exists.'

    else !create lut

      print*,'MSG makelut_visir_grnd; Create ',trim(fname_out(1))
      open(iu_lut,file=trim(sdsu_io_file),STATUS='REPLACE')
      do m = 0, 2*nda-1  ! Fourier order loop
         m1=m+1
         do iu = 1,mxpts_u10  !wind loop
            u10 = REAL(iu)*du10
            call ocnr11(m,nda,amua,wa,cr,ci,u10,r)
            do n1 = 1,nda ; do n2 = 1,nda
               write(iu_lut,*) r(n1,n2)
            enddo ; enddo
         enddo ! iu
      enddo ! m
      close(iu_lut)

    endif

 90 Format ( 1P1E24.15 )

!
! inquire SRR LUT exist or not (if not create LUT)
!
    sdsu_io_file = trim(sdsu_dir_sslut)//fname_out(2)
    inquire(file=trim(sdsu_io_file),exist=fexist)
    if(lut_grnd_replace) fexist = .false. !for debug

    if (fexist) then !read lut

      print*,'MSG makelut_visir_grnd; ',trim(fname_out(2)),' exists.'

    else !create lut

      print*,'MSG makelut_visir_grnd; Create ',trim(fname_out(2))
      open(iu_lut,file=trim(sdsu_io_file),STATUS='REPLACE')

      do m = 0, 2*nda-1  ! Fourier order loop
         m1=m+1
         do iu = 1,mxpts_u10  !wind loop
            u10 = REAL(iu)*du10     
            do ic = 0, mxpts_cosz
               cosz = REAL(ic)*0.01d0
               call ocnr31(m,nda,amua,wa,cosz,cr,ci,u10,sr)
               do n1 = 1,nda
                  write(iu_lut,*) sr(n1)
               enddo
            enddo !ic
         enddo ! iu
      enddo ! m
      close(iu_lut)
    
    endif
 
 enddo ! iw


 return
 end subroutine makelut_visir_grnd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine readlut_visir_grnd
 implicit none

 integer :: nda  !# of stream
 integer :: nw   !# of wavelength
 real(8) :: wl(knw)   !wavelength [micron]
 real(8) :: wlc          ! wavelength [cm]
 integer :: m1,m         ! Fourier order
 integer :: iw,iu,ic,n1,n2   ! loop indice
 integer :: ierr         ! allocation stat
 real(8) :: r(kndm,kndm)   ! reflection matrix
 real(8) :: sr(kndm)      ! source matrix
 real(8) :: u10          ! 10m wind speed [m/s] (0m/s cause numerical instability)
 real(8) :: cosz         ! cosine of solar zenith angle
 integer,parameter :: iu_lut = 2
 
 real(8) :: du10         ! increment of u10 [m/s]
 character(len=50) :: fname_out(2)
 logical :: fexist

 nda = N_stream    !from module_para
 nw = mxwavel      !from module_para
 wl(1:nw) = DBLE(wavel(1:nw))  !from module_para

 du10 = 25.d0 / REAL(mxpts_u10) !from module_para

!
! allocate orf_lut & srr_lut
!
  if( .not. allocated(orf_lut) ) then
      allocate( orf_lut(1:nda, 1:nda, 1:mxpts_u10, 1:nda*2,nw) , stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank, 'MSG grndo1: allocation error -> Terminate program. ')
  endif
  if( .not. allocated(srr_lut) ) then
      allocate( srr_lut(1:nda, 0:mxpts_cosz, 1:mxpts_u10, 1:nda*2,nw) , stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG grndo1: allocation error -> Terminate program. ')
  endif
 
!
! get gausian cosz and weights
!

 do iw=1,nw  ! wavelength loop

    wlc = wl(iw) * 1.0d-4              ! wavelength [micron] -> [cm]

    fname_out(1) = 'ORF_SURFACE_'//trim(nch_wavel(iw))//'.dat'
    fname_out(2) = 'SRR_SURFACE_'//trim(nch_wavel(iw))//'.dat'

!
! inquire ORF LUT exist or not (if not create LUT)
!
    sdsu_io_file = trim(sdsu_dir_sslut)//fname_out(1)    !from module_para
    inquire(file=trim(sdsu_io_file),exist=fexist)

    if (fexist) then !read lut

      if(masterproc .and. verbose_SDSU) print*,'MSG readlut_visir_grnd; Reading ',trim(fname_out(1))

!      open(iu_lut,file=trim(sdsu_io_file),STATUS='OLD')
      call open_oldfile(iu_lut, trim(sdsu_io_file) )
      do m = 0, 2*nda-1  ! Fourier order loop
         m1=m+1
         do iu = 1,mxpts_u10  !wind loop
            do n1 = 1,nda ; do n2 = 1,nda
               read(iu_lut,*) orf_lut(n1,n2,iu,m1,iw)
            enddo ; enddo
         enddo ! iu
      enddo ! m
      close(iu_lut)

    else !stop

      print*,'MSG readlut_visir_grnd',trim(fname_out(1)),'does not exisit.'
      call stop_sdsu(myrank,'') 

    endif

 90 Format ( 1P1E24.15 )

!
! inquire SRR LUT exist or not (if not create LUT)
!
    sdsu_io_file = trim(sdsu_dir_sslut)//fname_out(2)
    inquire(file=trim(sdsu_io_file),exist=fexist)

    if (fexist) then !read lut
      if(masterproc .and. verbose_SDSU) print*,'MSG readlut_visir_grnd; Reading ',trim(fname_out(2))
!      open(iu_lut,file=trim(sdsu_io_file),STATUS='OLD')
      call open_oldfile(iu_lut, trim(sdsu_io_file) )

      do m = 0, 2*nda-1  ! Fourier order loop
         m1=m+1
         do iu = 1,mxpts_u10  !wind loop
            do ic = 0, mxpts_cosz
               do n1 = 1,nda
                  read(iu_lut,*) srr_lut(n1,ic,iu,m1,iw)
               enddo
            enddo !ic
         enddo ! iu
      enddo ! m
      close(iu_lut)

    else !stop

      print*, 'MSG readlut_visir_grnd; ',trim(fname_out(2)),'does not exisit.'
       call stop_sdsu(myrank,'') 
    endif
 
 enddo ! iw

 return
 end subroutine readlut_visir_grnd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine water_rhi(wlc,crw,ciw)
 implicit none

! (nv,nwlv,wlv,rfi are global paramter from module pras)
 real(8),intent(in):: wlc      !! wavenumber [cm]
 real(8),intent(out):: crw     !! relative reflactive index (real part) of water
 real(8),intent(out):: ciw     !! relative reflactive index (imaginary part) of water

! for wavelength
 real(8):: wl1,dwl,pint
 integer:: iwl,iwlt

! for refractive index
  integer:: nintp
  real(8):: x(3),yr(3),yi(3)
  real(8):: bintp

! search wavelength in refractive index table
  wl1=wlc*1.0e4   ! wavenumber [cm] -> [micron] 
  dwl=log(wl1/wlv(1))/log(wlv(nwlv)/wlv(1))*(nwlv-1)+1  !
  iwl=int(dwl) !wavenumber index

  nintp=3
  if(nintp>nwlv) nintp=nwlv
  iwlt=iwl
  if(iwl>=nwlv-1) iwlt=nwlv-2

! water refractive index:

  x (1:nintp)=log(wlv(iwlt:iwlt+nintp-1))
  yr(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,1,1) )
  yi(1:nintp)=log(rfi(iwlt:iwlt+nintp-1,1,2) )
  pint=log(wl1)

  crw=exp( bintp(pint,nintp,x,yr) )  !output needed for water surface
  ciw=exp( bintp(pint,nintp,x,yi) )  !output needed for water surface

 return
 end subroutine water_rhi

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine ss_visir_gen_micro(spc, wave_um, wc, re, mu, dens, &
                                cextp, cabsp, php, success)
 implicit none
!--------------------------------------------------------------------------------
! Purpose : Compute the extinction, absorption, asymmetry parameter and
!           phase function for a given water content and effective radius of condensates , and
!           a particle size distribution parameter for generalized gamma distribution
!
! 07/2015   Toshi Matsui @ NASA GSFC: Initial
!--------------------------------------------------------------------------------
 character*(*),intent(in) :: spc !species character
 real(8),intent(in) :: wave_um ! wavelength [um]
 real(8),intent(in) :: wc    ! water content of condensate [g/m3]
 real(8),intent(in) :: re    ! drop effective radius [micron] 
 real(8),intent(in) :: mu    ! PSD shape parameter for generalized gamma distribution
 real(8),intent(in) :: dens  ! effective density of hydrometero [kg/m**3]

 real(8),intent(out):: cextp    ! extinction cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: cabsp    ! absorption cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: php(knang) ! scattering phase function [cm-1  str-1]
 logical,intent(out) :: success ! success criteria of interpolation


!--------Local Parameters
 integer :: i    !looping 
 integer :: imax !max for looping

 real(8) ::&
   rad   ,& ! radius of particle [mm]
   d     ,& ! particle diameter 
   num   ,& ! a particle number density per radius increment [1/m4]
   gamfac,& ! gamma PSD factors [-]
   gfac1 ,& !
   gfac2 ,& !
   lam   ,& ! slope lambda [1/m]
   n0    ,& ! intercept  [1/m**(4+mu)]
   d_increment,& ! diameter increment over size loops [mm]
   d_increment_m,&! diameter increment over size loops [m]
   d_min      ,& ! minium diameter of particle [mm]
   d_max      ,& ! maximum diameter of particle [mm]
   d_m        ,& ! particle diameter [mm]
   d_mm       ,& ! particle diameter [mm]
   mean_mass  ,& ! mean mass [kg]
   fgamma     ,& ! gamma function
   ntot       ,& ! total particle number concentrations [#/m3]
   factor        ! factor = num * pi * 0.25 * d_mm * d_mm * d_increment * 1.e-6


 real(8),parameter :: dr = 0.05          ! increment bin of radius [mm]
 real(8),parameter :: densice = 0.917d+3 ! solid ice density [kg/m3]
 real(8),parameter :: densliq = 1.0d+3   ! liquid water density [kg/m3]
 real(8) :: d_um_solid  !solid ice diameter [um]
 real(8) :: rad_solid_micron  !solid ice radius [micron]
 real(8) :: area   !area of fully particle [mm2]
 character(len=100) :: icescat_spc

 real :: qext, w0, p11(nScatAng), p22(nScatAng)

 real(8) :: cscat     !scattering cross section [cm-1]
 real(8) :: cext_sum  !sum of extinction cross section [cm-1]
 real(8) :: cabsp_sum !sum of  absorption cross section [cm-1]
 real(8) :: p11_sum(1:nScatAng)  !sum of P11 phase function [cm-1 str-1]
 real(8) :: p22_sum(1:nScatAng)  !sum of P22 phase function [cm-1 str-1]

!
! Begin by checking if hydrometeors of this species have enough water content.
! If not, set scattering parameters to zero and return.
!
  cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
  success = .true.
  if(wc .lt. q_min_conden) return

!
! prepare ICESCAT  (options are tunable)
!
 select case( trim(spc) ) 
 case('qi')
   icescat_spc = 'solid_bullet_rosettes'
 case('qs')
   icescat_spc = 'large_aggregate_plates'
 case('qg')
   icescat_spc = 'droxtals'
 case('qh')
   icescat_spc = 'droxtals'

 case default
   !false for spherical liquid species...
   success = .false.
   return
 end select

 call interp_wave_icescat_yang_visir(real(wave_um),icescat_spc,success)
 if(.not. success) return


!
! Get characteristic diameter from effective radius (re) and PSD shape parameter (mu)
!
  call gamma_reff_r8(mu+4.0,gfac1)
  call gamma_reff_r8(mu+3.0,gfac2)
  gamfac = gfac1/gfac2


! derive lambda
  lam = 1.d0 / (2.d0*re*1.d-6) * gamfac  !slope [1/m]

! derive intercept
!       [g/m3]                     [m3/g]        [1/m]
  n0  = wc * 6.d0 / pi / (dens*1000.d0) * (lam**(4.+mu))  / gfac1  ! [1/m**(4+mu)]


!
! loop bounds for size integration
!
 if( trim(spc) == 'qi' )  then  !non-precip

  d_increment = 0.001  ! [mm] 0.01mm=1micron
  d_min =  0.002  ![mm]  --> 0.002mm=2micron
  d_max = 0.2  ! [mm]  --> 0.2mm=200micron

 else  !precip species

  d_increment = 0.10  ! [mm]  0.1mm=100micron
  d_min =  0.2  ![mm]  --> 0.2mm=200micron
  d_max = 20.  ! [mm]  --> 20mm=2cm

 endif


!
! special case (Thompson snow class starts from ice crystal size)
!
  if( trim(cloud_microphysics) == 'THOM' .and. &
      trim(spc) == 'qs'                  ) then  !thompson spectrum density of snow
     d_increment = 0.01  ! [mm] --> 10micron
     d_min =  0.002  ![mm]  --> 2micron
     d_max = 20.  ! [mm]  --> 2cm
  endif

 d_increment_m = d_increment * 1.d-3  ! [m]

!
! loop bound
!
 imax = nint( (d_max-d_min) /d_increment)

!
! At this point, special thompson DSD is not applied, since it requires temperature information, 
! which is not dealt in the visible IR simulator optical routine...(toshi)
!
!  if( trim(cloud_microphysics) == 'THOM' .and. &
!      trim(spc) == 'qs'                  ) then  !thompson spectrum density of snow
!      if( .not. allocated(thom_snow_num) ) allocate(  thom_snow_num(0:imax)   )
!      call thom_snow_n (wc,temp,d_increment,d_min,d_max,imax, thom_snow_num )
!  endif


!
! initialize summing parameters
!
 cext_sum  = 0.   
 cabsp_sum = 0.   
 p11_sum(1:nScatAng) = 0.
 p22_sum(1:nScatAng) = 0.

!
! loop over particle size
!
  SIZE_LOOP: do i=0,imax

!
! Compute diameter, size parameter particle number density, and particle density
!
     d_mm = d_increment * 0.5 + d_increment * FLOAT(i)  !diameter [mm]
     d_m  = d_mm * 1e-3                        ! diameter but using different unit [m]

    !
    ! Get N(D)
    !    [1/m**(4+mu)]  [ m**mu ]      [1/m]  [m]   = [1/m4]
    num = n0           * d_m**mu * exp(-lam * d_m )   ! particle number density [1/m4]

!
! Thompson snow number density calculated in module_thompson_mp.F:
!
!     if( trim(cloud_microphysics) == 'THOM' .and. &
!         trim(spc) == 'qs'                  ) then
!       num = thom_snow_num(i)
!     endif


    !solid-sphere radius for interpolation
    d_um_solid = 1.d+3 * ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![um]



    !
    ! ICESCAT ption for ice crystal
    !
    call interp_size_icescat_yang_visir(real(d_um_solid) ,qext, w0, p11, p22 )

    !fluffy particle cross section for estimating efficiency  <-> cross section
    area = (d_m*0.5)*(d_m*0.5)*pi  !area of fully particle [m2]


    !
    ! intgrate for size distributions 
    !
    factor = num * d_increment_m * 1.d-2  ! 1.d-2 is unit convesion from [m-1] to [cm-1]

    cext_sum  = cext_sum  + area * dble(qext) * factor                      !sum of extinction cross section [cm-1]
!    print*,i,  num, factor, n0

    cabsp_sum = cabsp_sum + area * dble(qext) * factor * (1.d0 - dble(w0))  !sum of absorption cross section [cm-1]

    cscat = area * dble(qext) * factor * dble(w0)  !scattering cross section [cm-1]

    p11_sum(1:nScatAng) = p11_sum(1:nScatAng) + cscat*dble(p11(1:nScatAng))/( 4.d0 * pi ) 
    !sum of P11 phase function [cm-1 str-1]

    p22_sum(1:nScatAng) = p22_sum(1:nScatAng) + cscat*dble(p22(1:nScatAng))/( 4.d0 * pi ) 
    !sum of P22 phase function [cm-1 str-1]


  enddo SIZE_LOOP

!
! output
!
 cextp    = cext_sum  ! extinction cross section [cm-1] = [cm2/cm3]
 cabsp    = cabsp_sum ! absorption cross section [cm-1] = [cm2/cm3]
 php(1:nScatAng) =  p11_sum(1:nScatAng)   ! scattering phase function [cm-1  str-1]
 php(nScatAng+1) =  p22_sum(nScatAng)     ! P22 scattering phase function at 180deg [cm-1  str-1]
 
!print*, cextp,cabsp, php(nScatAng)
!stop

 return
 end subroutine ss_visir_gen_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine ss_visir_rams_micro(spc, wave_um, wc, re, gnu, cfmas, pwmas, &
                                cextp, cabsp, php, success)
 implicit none
!--------------------------------------------------------------------------------
! Purpose : Compute the extinction, absorption, asymmetry parameter and
!           phase function for a given water content and effective radius of condensates , and
!           a particle size distribution parameter for RAMS 1moment & 2moment microphyiscs
!
! 07/2015   Toshi Matsui @ NASA GSFC: Initial
!--------------------------------------------------------------------------------
 character*(*),intent(in) :: spc !species character
 real(8),intent(in) :: wave_um ! wavelength [um]
 real(8),intent(in) :: wc    ! water content of condensate [g/m3]
 real(8),intent(in) :: re    ! drop effective radius [micron] 
 real(8),intent(in) :: gnu   ! PSD shape parameter for generalized gamma distribution
 real(8),intent(in) :: cfmas  ! alpha_m [kg /m**beta_m]
 real(8),intent(in) :: pwmas  ! beta_m  [-]

 real(8),intent(out):: cextp    ! extinction cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: cabsp    ! absorption cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: php(knang) ! scattering phase function [cm-1  str-1]
 logical,intent(out) :: success ! success criteria of interpolation


!--------Local Parameters
 integer :: i    !looping 
 integer :: imax !max for looping

 real(8) ::&
   rad   ,& ! radius of particle [mm]
   d     ,& ! particle diameter 
   dens  ,& ! particle density [kg/m3]
   num   ,& ! a particle number density per radius increment [1/m4]
   gamfac,& ! gamma PSD factors [-]
   gfac1 ,& !
   gfac2 ,& !
   dn         ,& ! charactristic diameter [m]
   d_increment,& ! diameter increment over size loops [mm]
   d_increment_m,&! diameter increment over size loops [m]
   d_min      ,& ! minium diameter of particle [mm]
   d_max      ,& ! maximum diameter of particle [mm]
   d_m        ,& ! particle diameter [mm]
   d_mm       ,& ! particle diameter [mm]
   mean_mass  ,& ! mean mass [kg]
   fgamma     ,& ! gamma function
   ntot       ,& ! total particle number concentrations [#/m3]
   factor        ! factor = num * pi * 0.25 * d_mm * d_mm * d_increment * 1.e-6


 real(8),parameter :: dr = 0.05          ! increment bin of radius [mm]
 real(8),parameter :: densice = 0.917d+3 ! solid ice density [kg/m3]
 real(8),parameter :: densliq = 1.0d+3   ! liquid water density [kg/m3]
 real(8) :: d_um_solid  !solid ice diameter [um]
 real(8) :: rad_solid_micron  !solid ice radius [micron]
 real(8) :: area   !area of fully particle [mm2]
 character(len=100) :: icescat_spc

 real :: qext, w0, p11(nScatAng), p22(nScatAng)

 real(8) :: cscat     !scattering cross section [cm-1]
 real(8) :: cext_sum  !sum of extinction cross section [cm-1]
 real(8) :: cabsp_sum !sum of  absorption cross section [cm-1]
 real(8) :: p11_sum(1:nScatAng)  !sum of P11 phase function [cm-1 str-1]
 real(8) :: p22_sum(1:nScatAng)  !sum of P22 phase function [cm-1 str-1]

!
! Begin by checking if hydrometeors of this species have enough water content.
! If not, set scattering parameters to zero and return.
!
  cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
  success = .true.
  if(wc .lt. q_min_conden) return

!
! prepare ICESCAT  (options are tunable)
!
 select case( trim(spc) ) 
 case('qi1')
   icescat_spc = 'solid_bullet_rosettes'
 case('qi2')
   icescat_spc = 'solid_bullet_rosettes'
 case('qs')
   icescat_spc = 'large_aggregate_plates'
 case('qg')
   icescat_spc = 'droxtals'
 case('qh')
   icescat_spc = 'droxtals'

 case default
   !false for spherical liquid species...
   success = .false.
   return
 end select

 call interp_wave_icescat_yang_visir(real(wave_um),icescat_spc,success)
 if(.not. success) return


! wave=wave_um*1.e-3  !wavelgnth [mm] 

!
! Get characteristic diameter from effective radius (re) and PSD shape parameter (gnu)
!
  call gamma_reff_r8(gnu+3.0,gfac1)
  call gamma_reff_r8(gnu+2.0,gfac2)
  gamfac = gfac1/gfac2
  dn = 2.0 * re / gamfac * 1d-6  !charactristic diameter [m]

!
! get total number concentration for a given characteristic diameter, mass mixing ratio, PSD shape parameter
! , and gamma factors.
!

  call gamma_reff_r8(gnu+pwmas,gfac1)
  call gamma_reff_r8(gnu      ,gfac2)
  gamfac = gfac1/gfac2
  mean_mass = cfmas * ( dn ** pwmas ) * gamfac  !mean mass [kg]
  ntot = wc * 1e-3 / mean_mass                   !total particle number concentration [#/m3]


!
! loop bounds for size integration
!
 if( trim(spc) == 'qi1' )  then  !non-precip

  d_increment = 0.001  ! [mm] 0.01mm=1micron
  d_min =  0.002  ![mm]  --> 0.002mm=2micron
  d_max = 0.2  ! [mm]  --> 0.2mm=200micron

 else  !precip species

  d_increment = 0.10  ! [mm]  0.1mm=100micron
  d_min =  0.2  ![mm]  --> 0.2mm=200micron
  d_max = 20.  ! [mm]  --> 20mm=2cm

 endif

 d_increment_m = d_increment * 1.d-3  ! [m]

!
! loop bound
!
 imax = nint( (d_max-d_min) /d_increment)


!
! initialize summing parameters
!
 cext_sum  = 0.   
 cabsp_sum = 0.   
 p11_sum(1:nScatAng) = 0.
 p22_sum(1:nScatAng) = 0.

!
! loop over particle size
!
  SIZE_LOOP: do i=0,imax

!
! Compute diameter, size parameter particle number density, and particle density
!
     d_mm = d_increment * 0.5 + d_increment * FLOAT(i)  !diameter [mm]
     d_m  = d_mm * 1e-3                        ! diameter but using different unit [m]
     dens = (6.0 * cfmas / pi) * (d_m ** (pwmas-3.e0)) ! particle density [kg/m3]


!
! Gamma function
!
     call gamma_function_r8(gnu,d_m,dn, fgamma)   ! get gamma function
     num = ntot * fgamma                       ! particle number density [1/m4]


    !solid-sphere radius for interpolation
    d_um_solid = 1.d+3 * ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![um]

    !
    ! ICESCAT ption for ice crystal
    !
    call interp_size_icescat_yang_visir(real(d_um_solid) ,qext, w0, p11, p22 )

    !fluffy particle cross section for estimating efficiency  <-> cross section
    area = (d_m*0.5)*(d_m*0.5)*pi  !area of fully particle [m2]

    !
    ! intgrate for size distributions 
    !
    factor = num * d_increment_m * 1.d-2  ! 1.d-2 is unit convesion from [m-1] to [cm-1]

    cext_sum  = cext_sum  + area * dble(qext) * factor                      !sum of extinction cross section [cm-1]

    cabsp_sum = cabsp_sum + area * dble(qext) * factor * (1.d0 - dble(w0))  !sum of absorption cross section [cm-1]

    cscat = area * dble(qext) * factor * dble(w0)  !scattering cross section [cm-1]

    p11_sum(1:nScatAng) = p11_sum(1:nScatAng) + cscat*dble(p11(1:nScatAng))/( 4.d0 * pi ) 
    !sum of P11 phase function [cm-1 str-1]

    p22_sum(1:nScatAng) = p22_sum(1:nScatAng) + cscat*dble(p22(1:nScatAng))/( 4.d0 * pi ) 
    !sum of P22 phase function [cm-1 str-1]


  enddo SIZE_LOOP

!
! output
!
 cextp    = cext_sum  ! extinction cross section [cm-1] = [cm2/cm3]
 cabsp    = cabsp_sum ! absorption cross section [cm-1] = [cm2/cm3]
 php(1:nScatAng) =  p11_sum(1:nScatAng)   ! scattering phase function [cm-1  str-1]
 php(nScatAng+1) =  p22_sum(nScatAng)     ! P22 scattering phase function at 180deg [cm-1  str-1]
 

 return
 end subroutine ss_visir_rams_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine ss_visir_sbm_micro(spc, wave_um, rad_bin, drad_bin, wc, n0_bin, rho_bin, &
                                cextp, cabsp, php, success)
 implicit none
!--------------------------------------------------------------------------------
! Purpose : Compute the extinction, absorption, asymmetry parameter and
!           phase function for a given water content and effective radius of condensates , and
!           a particle size distribution parameter for HUCM SBM 
!
! 07/2015   Toshi Matsui @ NASA GSFC: Initial
!--------------------------------------------------------------------------------
 character*(*),intent(in) :: spc !species character
 real(8),intent(in) :: wave_um ! wavelength [um]
 real(8),intent(in):: rad_bin(nbin)     ! radius of size bin [cm]
 real(8),intent(in):: drad_bin(nbin)    ! width of size bin [cm]
 real(8),intent(in) :: wc               ! water content of condensate [g/m3]
 real(8),intent(in):: n0_bin(nbin)      ! dN/dr: # conc per radius [cm-4] for each nbin
 real(8),intent(in):: rho_bin(nbin)     ! effective density of size bin [g/cm3]

 real(8),intent(out):: cextp    ! extinction cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: cabsp    ! absorption cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: php(knang) ! scattering phase function [cm-1  str-1]
 logical,intent(out) :: success ! success criteria of interpolation


!--------Local Parameters
 integer :: i    !looping 
 integer :: imax !max for looping

 real(8) ::&
   rad   ,& ! radius of particle [mm]
   d     ,& ! particle diameter 
   dens  ,& ! particle density [kg/m3]
   num   ,& ! a particle number density per radius increment [1/m4]
   gamfac,& ! gamma PSD factors [-]
   gfac1 ,& !
   gfac2 ,& !
   dn         ,& ! charactristic diameter [m]
   d_increment_m,&! diameter increment over size loops [m]
   d_min      ,& ! minium diameter of particle [mm]
   d_max      ,& ! maximum diameter of particle [mm]
   d_m        ,& ! particle diameter [mm]
   d_mm       ,& ! particle diameter [mm]
   mean_mass  ,& ! mean mass [kg]
   fgamma     ,& ! gamma function
   ntot       ,& ! total particle number concentrations [#/m3]
   factor        ! factor = num * pi * 0.25 * d_mm * d_mm * d_increment * 1.e-6


 real(8),parameter :: dr = 0.05          ! increment bin of radius [mm]
 real(8),parameter :: densice = 0.917d+3 ! solid ice density [kg/m3]
 real(8),parameter :: densliq = 1.0d+3   ! liquid water density [kg/m3]
 real(8) :: d_um_solid  !solid ice diameter [um]
 real(8) :: rad_solid_micron  !solid ice radius [micron]
 real(8) :: area   !area of fully particle [mm2]
 character(len=100) :: icescat_spc

 real :: qext, w0, p11(nScatAng), p22(nScatAng)

 real(8) :: cscat     !scattering cross section [cm-1]
 real(8) :: cext_sum  !sum of extinction cross section [cm-1]
 real(8) :: cabsp_sum !sum of  absorption cross section [cm-1]
 real(8) :: p11_sum(1:nScatAng)  !sum of P11 phase function [cm-1 str-1]
 real(8) :: p22_sum(1:nScatAng)  !sum of P22 phase function [cm-1 str-1]

!
! Begin by checking if hydrometeors of this species have enough water content.
! If not, set scattering parameters to zero and return.
!
  cextp = 0.d0 ; cabsp = 0.d0 ; php = 0.d0
  success = .true.
  if(wc .lt. q_min_conden) return

!
! prepare ICESCAT  (options are tunable)
!
 select case( trim(spc) ) 
 case('ice_col')
   icescat_spc = 'solid_columns'
 case('ice_pla')
   icescat_spc = 'plates'
 case('ice_den')
   icescat_spc = 'plates'
 case('snow')
   icescat_spc = 'large_aggregate_plates'
 case('graupel')
   icescat_spc = 'droxtals'
 case('hail')
   icescat_spc = 'droxtals'
 case default
   !false for spherical liquid species...
   success = .false.
   return
 end select

 call interp_wave_icescat_yang_visir(real(wave_um),icescat_spc,success)
 if(.not. success) return

!
! initialize summing parameters
!
 cext_sum  = 0.   
 cabsp_sum = 0.   
 p11_sum(1:nScatAng) = 0.
 p22_sum(1:nScatAng) = 0.

!
! loop over particle size
!
  SIZE_LOOP: do i=1, nbin

!
! Compute diameter, size parameter particle number density, and particle density
!
     d_mm = rad_bin(i)*20.  !diameter [mm]
     d_m  = d_mm / 1000.       ! diameter [m]
     dens = rho_bin(i) * 1.d+3 ! particle density [kg/m3]
     d_increment_m = drad_bin(i) / 100.


!
! Gamma function
!
     num = n0_bin(i) * 1.d+8     ! particle number density [1/m4]


    !solid-sphere radius for interpolation
    d_um_solid = 1.d+3 * ( d_mm * d_mm * d_mm * (dens/densice) ) ** (1./3.)  ![um]

    !
    ! ICESCAT ption for ice crystal
    !
    call interp_size_icescat_yang_visir(real(d_um_solid) ,qext, w0, p11, p22 )

    !fluffy particle cross section for estimating efficiency  <-> cross section
    area = (d_m*0.5)*(d_m*0.5)*pi  !area of fully particle [m2]

    !
    ! intgrate for size distributions 
    !
    factor = num * d_increment_m * 1.d-2  ! 1.d-2 is unit convesion from [m-1] to [cm-1]

    cext_sum  = cext_sum  + area * dble(qext) * factor                      !sum of extinction cross section [cm-1]

    cabsp_sum = cabsp_sum + area * dble(qext) * factor * (1.d0 - dble(w0))  !sum of absorption cross section [cm-1]

    cscat = area * dble(qext) * factor * dble(w0)  !scattering cross section [cm-1]

    p11_sum(1:nScatAng) = p11_sum(1:nScatAng) + cscat*dble(p11(1:nScatAng))/( 4.d0 * pi ) 
    !sum of P11 phase function [cm-1 str-1]

    p22_sum(1:nScatAng) = p22_sum(1:nScatAng) + cscat*dble(p22(1:nScatAng))/( 4.d0 * pi ) 
    !sum of P22 phase function [cm-1 str-1]


  enddo SIZE_LOOP

!
! output
!
 cextp    = cext_sum  ! extinction cross section [cm-1] = [cm2/cm3]
 cabsp    = cabsp_sum ! absorption cross section [cm-1] = [cm2/cm3]
 php(1:nScatAng) =  p11_sum(1:nScatAng)   ! scattering phase function [cm-1  str-1]
 php(nScatAng+1) =  p22_sum(nScatAng)     ! P22 scattering phase function at 180deg [cm-1  str-1]
 

 return
 end subroutine ss_visir_sbm_micro

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_opt_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function bintp(xx,n,x,y)
! 2nd order polynomial interpolation
!--- history
! 88. 1.14 created
!     3. 7 linear extrapolation for out-ou-boundary data
! 08. 9.19 modified Fortran 90 free-style form
!--- input
! xx    r     interpolation point
! n     r     nbr of data
! x   r(nn)   independent variable,  x(i) .lt. x(i+1)
! y   r(nn)   dependent   variable
!--- output
! bintp r     output
!--
  implicit none

! for input & output
  integer,intent(in)::n
  real(8),intent(in)::xx,x(n),y(n)
  real(8):: bintp
! for work
  integer::i,i1,i2,i3,k1,k2
!--exec
! 1 point
  if(n <= 1) then
     bintp=y(1)
     return
  endif

! 2 points or xx.le.x(1)
  if(xx <= x(1) .or. n==2) then
     bintp=y(1)+(y(2)-y(1))*(xx-x(1))/(x(2)-x(1))
     return
  endif
! xx.ge.xx(n)
  if(xx > x(n)) then
     bintp=y(n-1)+(y(n)-y(n-1))*(xx-x(n-1))/(x(n)-x(n-1))
     return
  endif

! 3 points
  do i=1,n
     if(xx < x(i)) goto 10
  enddo
  i=n
10 k1=max(1,i-1); k2=k1+2
  if(k2 > n) k2=n; k1=k2-2

  bintp=0.0d0; i2=k1; i3=k1+1
  do i=1,3
     i1=i2; i2=i3; i3=mod(i+1,3)+k1
     bintp=bintp+(xx-x(i2))*(xx-x(i3))/(x(i1)-x(i2))/(x(i1)-x(i3))*y(i1)
  enddo
  return
end function bintp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function volume_spectrum(pr)
use module_simulator
! volume spectrum of partile polydisperison: v(x) = dv / d ln r  [cm3/cm3]
!--- history
! 95. 9.14 created with parameter packet only
!--- input
! pr    r(6,33)      PSD parameter packet
!
!    pr(1,1)=r       particle radius in cm
!    pr(1,2)=nmode   number of mode radius
!    pr(1,3)=rmin    minimum particle radius in cm
!    pr(1,4)=rmax    maximum particle radius in cm
!
!
!  pr(2,1): type of function (itp) for the mode.
!   itp=1: power law
!     pr(3,1)=c, pr(4,1)=r0,  pr(5,1)=p
!     vj = c * (r/r0)**(4-p) if r>r0; = c * (r/r0)**4 if r<r0
!   itp=2: log-normal
!     pr(3,1)=c, pr(4,1)=s,   pr(5,1)=rm
!     vj = c * exp((ln(r/rm)/ln(s))**2 / 2)
!
!   itp=3: modified gamma
!     pr(3,1)=c, pr(4,1)=alfa, pr(5,1)=beta, pr(6,1)=gamma
!     vj = c * (r1)**(alfa+4) exp (-beta*r1**gamma) where r1=r*1.0e4
!
!   itp=4: MARSHALL-PALMER (MODIFIED by H.MASUNAGA (01/06/21))
!     pr(3,1) = n0 [cm^-4], PR(4,j) = lambda [cm^-1]
!     dN/dr = n0 * exp(-lambda*D) , D=2*r  !(dN/dr is [1/cm4])
!     vj = 4 * pi * r**4 * dN/dr / 3
!
!   itp=5: Explicite spectra-bin microphysics (modified by T.Matsui (07/2007))
!     dN/dr = n0(r)
!     pr(3,1:33) = n0(1:33) [cm^-4]
!     pr(4,1:33) = rad(1:33) [cm] 
!     pr(5,1:33) = drad(1:33) [cm] 
!     vj = 4 * pi * r**4 * dN/dr / 3
!   
!   itp=6: generalized gamma distributin  (MODIFIED by T.Matsui (04/2009))
!     pr(3,1) = ntot [m^-3], pr(4,j) = dn [m] , pr(5,j) = gnu
!     dN/dr = ntot * fgamma * 1.0d-8       ! # conc per radius  [1/cm4] 
!     vj = 4 * pi * r**4 * dN/dr / 3
!
!--- output
! volume_spectrum   rf       dv/d ln r  [cm3/cm3]
!--

  implicit none
  real(8),intent(in)::pr(6,33)    !! parameter packet
  integer::itp,m
  real(8)::r,rmin,rmax,e1,r1
  real(8)::c,rc,pdndr,pn,s,rm,alf,bet,gam,diamtr,dndr,n0,lambda
  real(8)::volume_spectrum
  real(8),parameter :: pi=3.14159265358979323846d0
  real(8):: left, right  ![cm]
  integer :: ibin !for SBM size bin loop

  real(8) :: d_m, ntot, dn, gnu,fgamma  !for RAMS generlized gamma distribution

  real(8) :: Nt    ! total conc [1/cm3]
  real(8) :: sigma ! s: log standard deviation [-]
  real(8) :: mu    ! mean --> log(rm) [cm] <- [um] 
  real(8) :: exp_func, const

!--exec
  r=pr(1,1)  !radius [cm] 
  rmin =pr(1,3)  !minimum radius [cm]
  rmax =pr(1,4)  !maximum radisu [cm]

  volume_spectrum=0.0d0   !initialize total particle volume as zero.

!  print*,'vlspc', rmin, r, rmax
  if(r < rmin .or. r > rmax) return   !size bounds here

     itp=int(pr(2,1))  !PSD type

     select_psd: select case(itp)
     case(1)  !power law
        c=pr(3,1); rc=pr(4,1); pdndr=pr(5,1)
        if(r <= rc) then
           pn=4.d0
        else
           pn=4.d0-pdndr
        endif
        e1=pn*log(r/rc)
        if(e1 > -100.d0) volume_spectrum = c * exp(e1)

     case(2) ! log-normal
        c =pr(3,1); s =pr(4,1); rm=pr(5,1)
        e1=-0.5d0*(log(r/rm)/log(s))**2
        if(e1 > -100.d0) volume_spectrum = c * exp(e1)

     case(3) ! modified gamma
        r1=r*1.0d4; c=pr(3,1); alf=pr(4,1); bet=pr(5,1); gam=pr(6,1)
        e1=(alf+4.d0)*log(r1)-bet*r1**gam
        if(e1 > -100.d0) volume_spectrum = c * exp(e1)

     case(4) ! marshall-palmer --> gamma distribution
         diamtr = r * 2. ; n0=pr(3,1) ; lambda=pr(4,1) !; gnu = pr(5,1)
         dndr = n0 * exp(-lambda * diamtr)  ![1/cm4]
!        dndr = n0 * (dimtr**(-gnu)) * exp(-lambda * diamtr)  ![1/cm4] gamma formulae check later
         volume_spectrum = 4.d0 * pi * dndr * r**4 / 3.d0     ! [cm3/cm3] 
     
     case(5)  ! spectra-bin microphysics
         do ibin = 1, 33 
           left  = pr(4,ibin) - 0.5d0 * pr(5,ibin)   !left boundary of size bin
           right = pr(4,ibin) + 0.5d0 * pr(5,ibin)  !right boundary of size bin 
           
           if(left < 0.d0 ) call stop_sdsu(myrank, 'MSG volume_spectrum: left becomes negative')
           if( r >= left .and. r <= right) then  
               dndr = pr(3,ibin) !# conc per radius [1/cm4]
               volume_spectrum = 4.d0 * pi * dndr * r**4 / 3.d0  ![cm3/cm3]
               exit
           endif
         enddo 

     case(6) ! generalized gamma distribution

         d_m= 2.d0*r*1.d-2  ; ntot=pr(3,1) ; dn = pr(4,1) ; gnu = pr(5,1)      
!     pr(3,1) = ntot [m^-3], pr(4,j) = dn [m] , pr(5,j) = gnu
        ! d_m = 1d-3

         call gamma_function_r8(gnu,d_m,dn, fgamma)   ! get gamma function ( d_n [m], dn [m] )

         dndr =  ntot * fgamma * 1.0d-8 * 2.  

         volume_spectrum = 4.d0 * pi * dndr * r**4 / 3.d0   ![cm3/cm3]

     case(7)  ! mono (single size) distribution

! psd_params(3,1) = 3.* q / (4. * const_dpi * (rbin**3.) * rop ) / dr  ! dN/dr [1/cm4]
! psd_params(4,1) = rbin * 1.0d-4   ! radius [cm]
! psd_params(5,1) = dr              ! dr  [cm]

         left  = pr(4,1) - 0.5d0 * pr(5,1)  !left boundary of size bin [cm]
         right = pr(4,1) + 0.5d0 * pr(5,1)  !right boundary of size bin  [cm]
         if(left < 0.d0 ) call stop_sdsu(myrank, 'MSG volume_spectrum: left becomes negative')

         if( r >= left .and. r <= right) then
             dndr = pr(3,1)                        ! # conc per radius [1/cm4]
!             volume_spectrum = 4.d0 * pi * dndr * r**4 / 3.d0  ! [cm3/cm3]
              volume_spectrum = 4.d0 * pi * dndr * r**3 / 3.d0 * pr(5,1)  ! [cm3/cm3]

             !print*,'volume_spectrum', r, dndr, volume_spectrum
         endif

     case(8) !new log-normal

!    psd_params(3,1) = Nt ! total number concentration [1/cm3]
!    psd_params(4,1) = sigma  ! s: log standard deviation [-]
!    psd_params(5,1) = log( rm*1.0d-4 )  !mu: mean --> log(rm) [cm] <- [um] 

        Nt    = pr(3,1)
        sigma = pr(4,1)
        mu    = pr(5,1)

        exp_func = exp( -( (log(r) - mu)**2 ) / (2.*sigma*sigma) )

!        const = Nt / ( sqrt(2.*pi) * sigma * r )   !found bug this is dN/dr
        const = Nt / ( sqrt(2.*pi) * sigma  ) ! -->  dN/dln(r)

        volume_spectrum = const * exp_func * 4. * pi * (r **3.) / 3.  ! dV/dln(r)

     case default ; call stop_sdsu(myrank, 'MSG function volume_spectrum: There is no such itp')
     end select select_psd
      

  return
end function volume_spectrum

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function bplnk(w,t)
! plank function with respect to wavelength
!--- history
! 89. 8. 1  created
! 08. 9.25 modified Fortran 90 free-style form
!--- input 
! w      r      wavelength (micron)
! t      r      absolute temperature (K)
!--- output
! bplnk  rf     plank function (W/m2/str/micron)
!--
  implicit none

! input 
  real(8),intent(in)::w,t

  real(8)::x,expfn,bplnk
!--exec
  if(w*t <= 0.d0) then
     bplnk=0.d0
  else
     x=1.438786d4/w/t
     bplnk=1.1911d8/w**5/(expfn(x)-1.d0)
  endif
  return
end function bplnk

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function sunir(v)
!--- history
! 92. 4.24  registered from lowtran-7
! 03. 9. 9 Changed for free source form (Fortran 90) by M. Sekiguchi 
!--- input
! v        r        wavenumber (cm-1) 0 to 57490 cm-1
!--- output
! sunir    rf       the extra-terrestrial solar irradiance (watt/m2/micron)  
!                   if 0 then v is out of range
!--  uses  block data solar   which contains the values for solara +
!
!  use sundat, only: solara,solarb,solar
  implicit none
  real(8),intent(in)::v
  integer::i
  real(8)::p,wm,w0,w1,w2
  real(8),parameter::a=3.50187d-13,b=3.93281
  real(8)::sunir
  real(8),save::solara(1440),solarb(2910)

!
!--exec
!     wm, w0, w1, w2  are statement functions used by
!     the 4 point lagrange interpolation
  wm(p) = p*(p - 1)*(p - 2)
  w0(p) = 3*(p**2 - 1)*(p - 2)
  w1(p) = 3*p*(p + 1)*(p - 2)
  w2(p) = p*(p**2 - 1)

!  call solar(solara,solarb)
!     solar spectrum from      0 to    800 cm-1,  in steps of 20 cm-
    solara(1:41)=(/ &
         0.0000d+00, 4.5756d-08, 7.0100d-07, 3.4580d-06, 1.0728d-05, &
         2.5700d-05, 5.2496d-05, 9.6003d-05, 1.6193d-04, 2.5766d-04, &
         3.9100d-04, 5.6923d-04, 8.0203d-04, 1.1006d-03, 1.4768d-03, &
         1.9460d-03, 2.5213d-03, 3.2155d-03, 4.0438d-03, 5.0229d-03, &
         6.1700d-03, 7.5145d-03, 9.0684d-03, 1.0853d-02, 1.2889d-02, &
         1.5213d-02, 1.7762d-02, 2.0636d-02, 2.3888d-02, 2.7524d-02, &
         3.1539d-02, 3.5963d-02, 4.0852d-02, 4.6236d-02, 5.2126d-02, &
         5.8537d-02, 6.5490d-02, 7.3017d-02, 8.1169d-02, 9.0001d-02, &
         9.9540d-02 /)
!     solar spectrum from    820 to   3680 cm-1,  in steps of 20 cm-
    solara(42:185)=(/ &
         .10980, .12080, .13260, .14520, .15860, .17310, .18850, .20490, &
         .22240, .24110, .26090, .28200, .30430, .32790, .35270, .37890, &
         .40650, .43550, .46600, .49800, .53160, .56690, .60390, .64260, &
         .68320, .72560, .76990, .81620, .86440, .91470, .96710, 1.0220, &
         1.0780, 1.1370, 1.1990, 1.2630, 1.3290, 1.3990, 1.4710, 1.5460, &
         1.6250, 1.7060, 1.7910, 1.8800, 1.9710, 2.0670, 2.1660, 2.2680, &
         2.3740, 2.4840, 2.5970, 2.7140, 2.8350, 2.9600, 3.0890, 3.2210, &
         3.3570, 3.4980, 3.6420, 3.7900, 3.9440, 4.1040, 4.2730, 4.4450, &
         4.6150, 4.7910, 4.9830, 5.1950, 5.4210, 5.6560, 5.8930, 6.1270, &
         6.3560, 6.5820, 6.8080, 7.0360, 7.2700, 7.5170, 7.7890, 8.0910, &
         8.4070, 8.7120, 8.9900, 9.2490, 9.5000, 9.7550, 10.010, 10.250, &
         10.480, 10.700, 10.950, 11.230, 11.550, 11.900, 12.250, 12.600, &
         12.930, 13.250, 13.530, 13.780, 14.040, 14.320, 14.660, 15.070, &
         15.530, 16.011, 16.433, 16.771, 17.077, 17.473, 17.964, 18.428, &
         18.726, 18.906, 19.141, 19.485, 19.837, 20.160, 20.509, 21.024, &
         21.766, 22.568, 23.190, 23.577, 23.904, 24.335, 24.826, 25.236, &
         25.650, 26.312, 27.208, 27.980, 28.418, 28.818, 29.565, 30.533, &
         31.247, 31.667, 32.221, 33.089, 33.975, 34.597, 35.004, 35.395/)
!     solar spectrum from   3700 to   6560 cm-1,  in steps of 20 cm-
    solara(186:329)=(/ &
         36.026, 36.985, 37.890, 38.401, 38.894, 39.857, 40.926, 41.570, &
         42.135, 43.083, 44.352, 45.520, 45.982, 46.281, 48.335, 51.987, &
         54.367, 54.076, 52.174, 50.708, 52.153, 55.707, 56.549, 54.406, &
         53.267, 56.084, 61.974, 64.406, 60.648, 55.146, 53.067, 57.476, &
         64.645, 68.348, 69.055, 69.869, 70.943, 71.662, 72.769, 74.326, &
         75.257, 74.883, 73.610, 73.210, 74.886, 78.042, 80.204, 80.876, &
         82.668, 84.978, 86.244, 88.361, 91.998, 95.383, 98.121, 100.29, &
         100.64, 99.997, 101.82, 105.06, 107.50, 109.99, 112.45, 113.90, &
         113.79, 119.23, 121.96, 124.58, 127.14, 125.19, 124.37, 125.00, &
         127.88, 130.67, 131.98, 133.74, 136.69, 136.18, 135.02, 137.44, &
         138.44, 137.25, 136.35, 142.60, 144.54, 148.37, 151.90, 151.55, &
         155.35, 157.59, 159.70, 162.28, 168.44, 171.43, 169.82, 170.33, &
         172.28, 176.68, 181.92, 186.06, 187.85, 186.00, 189.82, 189.35, &
         192.86, 202.00, 209.63, 205.76, 212.88, 215.63, 216.51, 219.20, &
         220.29, 221.12, 227.12, 229.97, 233.23, 233.95, 234.52, 234.45, &
         235.77, 239.80, 243.11, 241.19, 242.34, 243.69, 242.84, 246.19, &
         246.11, 246.76, 251.75, 255.38, 258.74, 260.26, 263.40, 268.68, &
         271.81, 272.95, 273.93, 274.74, 274.43, 279.69, 287.76, 287.72/)
!     solar spectrum from   6580 to   9440 cm-1,  in steps of 20 cm-
    solara(330:473)=(/ &
         287.96, 290.01, 291.92, 295.28, 296.78, 300.46, 302.19, 299.14, &
         301.43, 305.68, 309.29, 310.63, 313.24, 314.61, 309.58, 318.81, &
         320.54, 321.62, 328.58, 331.66, 337.20, 345.62, 345.54, 342.96, &
         344.38, 346.23, 349.17, 351.79, 354.71, 356.97, 358.29, 362.29, &
         364.15, 364.97, 367.81, 368.98, 369.07, 372.17, 377.79, 381.25, &
         384.22, 388.66, 393.58, 396.98, 398.72, 400.61, 404.06, 408.23, &
         412.47, 415.58, 416.17, 416.53, 419.55, 425.88, 433.30, 437.73, &
         438.13, 439.79, 441.51, 438.71, 434.25, 437.54, 448.95, 448.86, &
         439.46, 437.10, 439.34, 444.33, 455.00, 467.05, 473.04, 469.64, &
         467.53, 473.78, 477.50, 477.50, 480.96, 483.94, 482.19, 479.08, &
         482.09, 493.43, 498.40, 492.05, 489.53, 493.34, 495.51, 496.52, & 
         499.57, 504.65, 509.68, 512.00, 512.05, 512.31, 515.00, 520.70, &
         527.30, 531.88, 532.16, 530.48, 532.33, 539.26, 548.57, 553.00, &
         548.96, 546.05, 551.00, 556.41, 557.21, 557.85, 560.95, 564.02, &
         565.57, 566.38, 567.88, 571.48, 576.68, 581.54, 586.51, 593.62, &
         600.70, 602.79, 601.39, 603.00, 606.88, 605.95, 600.97, 600.79, &
         607.21, 612.87, 614.13, 614.39, 616.61, 620.53, 625.19, 629.78, &
         633.79, 637.31, 640.47, 642.53, 642.62, 641.93, 643.11, 646.68/)
!     solar spectrum from   9460 to  12320 cm-1,  in steps of 20 cm-
    solara(474:617)=(/ &
         650.57, 654.30, 660.95, 672.10, 682.31, 684.89, 682.20, 682.53, &
         687.79, 691.42, 689.62, 688.14, 693.71, 703.25, 708.07, 706.22, &
         704.64, 708.97, 717.35, 725.43, 731.08, 734.17, 735.41, 736.60, &
         739.34, 742.90, 745.04, 744.29, 742.44, 749.53, 755.70, 758.82, &
         766.31, 761.53, 762.09, 769.68, 764.18, 763.75, 768.88, 762.69, & 
         753.93, 762.38, 765.79, 772.19, 760.67, 762.10, 766.76, 766.98, &
         769.35, 773.50, 766.84, 763.60, 773.82, 777.18, 779.61, 792.48, &
         797.54, 787.81, 793.75, 805.96, 804.77, 806.62, 821.72, 830.28, &
         827.54, 831.06, 830.20, 826.22, 823.28, 822.18, 833.92, 854.58, &
         859.80, 862.56, 871.16, 875.16, 867.67, 863.87, 883.30, 893.40, &
         897.74, 905.24, 905.38, 911.07, 930.21, 939.24, 934.74, 935.15, &
         942.38, 948.13, 947.00, 951.88, 960.12, 951.88, 954.22, 959.07, &
         963.36, 980.16, 983.66, 978.76, 979.38, 985.24, 977.08, 919.94, &
         899.68, 962.91, 997.17, 999.93, 995.65, 999.93, 1014.9, 951.57, &
         893.52, 955.14, 1003.1, 990.13, 978.79, 1011.2, 1034.7, 1031.9, &
         1029.9, 1039.7, 1045.5, 1044.1, 1049.6, 1056.1, 1049.8, 1038.0, &
         1051.9, 1072.2, 1075.5, 1077.0, 1079.3, 1078.0, 1075.7, 1079.7, &
         1081.0, 1069.8, 1078.4, 1104.3, 1111.4, 1111.7, 1117.6, 1119.6/)

!     solar spectrum from  12340 to  15200 cm-1,  in steps of 20 cm-
    solara(618:761)=(/ &
         1109.3, 1100.6, 1112.9, 1122.7, 1119.5, 1123.9, 1136.1, 1143.7, &
         1140.5, 1141.2, 1151.5, 1148.7, 1138.3, 1141.0, 1150.6, 1160.1, &
         1170.6, 1177.7, 1179.8, 1181.7, 1182.4, 1179.8, 1181.8, 1188.3, &
         1190.0, 1191.4, 1197.0, 1196.0, 1192.2, 1200.6, 1210.4, 1209.1, &
         1207.5, 1205.3, 1193.3, 1192.9, 1220.0, 1243.3, 1245.4, 1241.5, &
         1240.2, 1241.1, 1244.0, 1248.5, 1253.2, 1257.1, 1259.9, 1261.9, &
         1263.6, 1265.7, 1269.6, 1277.0, 1284.2, 1284.4, 1282.7, 1287.2, &
         1286.8, 1272.3, 1262.2, 1270.7, 1288.8, 1304.8, 1311.8, 1312.2, &
         1314.4, 1320.2, 1326.2, 1328.4, 1325.3, 1322.5, 1325.4, 1334.6, & 
         1346.4, 1354.0, 1353.7, 1347.3, 1338.3, 1331.0, 1329.7, 1338.0, &
         1351.9, 1363.0, 1368.8, 1372.0, 1375.9, 1382.1, 1387.8, 1388.8, &
         1388.2, 1392.2, 1401.7, 1412.9, 1418.2, 1410.7, 1395.9, 1385.7, &
         1388.1, 1405.0, 1424.0, 1428.1, 1422.2, 1423.6, 1434.5, 1445.2, &
         1450.7, 1451.8, 1451.5, 1453.9, 1459.9, 1466.9, 1471.3, 1469.4, &
         1462.5, 1460.4, 1468.9, 1481.8, 1490.8, 1495.3, 1497.9, 1500.7, &
         1505.2, 1510.0, 1512.3, 1512.7, 1515.6, 1521.6, 1524.2, 1520.7, &
         1520.3, 1531.6, 1545.7, 1548.2, 1541.7, 1542.2, 1553.6, 1563.6, &
         1563.6, 1559.9, 1561.3, 1569.9, 1581.6, 1577.6, 1529.7, 1447.0/)

!    solar spectrum from  15220 to  18080 cm-1,  in steps of 20 cm-
    solara(762:905)=(/ &
         1396.9, 1428.7, 1506.4, 1567.1, 1594.0, 1606.1, 1613.5, 1609.0, &
         1588.6, 1567.8, 1567.3, 1587.2, 1610.2, 1624.4, 1630.2, 1630.9, &
         1628.1, 1622.3, 1616.9, 1618.9, 1631.6, 1648.1, 1658.2, 1659.7, &
         1658.1, 1658.0, 1659.4, 1660.4, 1659.2, 1653.7, 1645.3, 1642.1, &
         1652.7, 1674.2, 1694.1, 1700.6, 1703.4, 1697.6, 1654.5, 1644.4, &
         1661.6, 1676.3, 1707.7, 1703.1, 1710.8, 1732.3, 1716.5, 1719.6, &
         1729.6, 1683.1, 1628.5, 1683.5, 1727.0, 1707.8, 1689.4, 1698.4, &
         1733.1, 1737.8, 1714.1, 1734.6, 1750.1, 1750.1, 1760.3, 1764.3, &
         1765.3, 1769.4, 1779.9, 1793.0, 1765.1, 1729.4, 1745.9, 1753.4, &
         1758.1, 1775.0, 1768.4, 1767.9, 1789.5, 1806.6, 1799.3, 1782.6, &
         1779.3, 1792.1, 1809.7, 1808.0, 1794.4, 1818.6, 1774.2, 1648.5, &
         1674.3, 1789.3, 1847.2, 1848.3, 1812.9, 1796.4, 1840.3, 1868.3, &
         1864.6, 1873.2, 1872.2, 1856.0, 1845.0, 1842.4, 1823.9, 1795.1, &
         1819.6, 1861.5, 1857.7, 1838.6, 1840.5, 1863.5, 1876.8, 1884.4, &
         1894.9, 1875.2, 1821.2, 1779.4, 1810.2, 1855.3, 1831.8, 1837.3, &
         1882.3, 1866.4, 1819.6, 1804.8, 1831.4, 1861.6, 1867.1, 1862.9, &
         1851.9, 1834.7, 1835.2, 1845.1, 1831.9, 1803.6, 1792.5, 1821.8, &
         1845.8, 1832.3, 1847.6, 1894.2, 1909.2, 1901.0, 1891.2, 1869.9/)

!    solar spectrum from  18100 to  20960 cm-1,  in steps of 20 cm-
    solara(906:1049)=(/ &
         1854.4, 1865.8, 1873.7, 1868.8, 1881.7, 1897.1, 1884.2, 1856.2, &
         1840.6, 1855.1, 1885.3, 1903.6, 1900.1, 1887.4, 1887.7, 1879.0, &
         1844.5, 1844.1, 1877.1, 1847.3, 1785.1, 1792.6, 1848.7, 1894.4, & 
         1908.8, 1892.8, 1867.4, 1885.6, 1959.9, 1971.9, 1895.8, 1883.5, &
         1917.6, 1853.8, 1793.0, 1875.6, 1974.0, 1975.7, 1943.9, 1926.4, &
         1914.4, 1902.7, 1882.5, 1813.3, 1710.8, 1717.9, 1859.7, 1965.1, &
         1970.1, 1941.4, 1902.5, 1852.0, 1836.3, 1879.3, 1901.6, 1862.9, &
         1839.1, 1840.9, 1780.0, 1684.9, 1677.3, 1718.7, 1697.3, 1684.3, &
         1784.5, 1898.0, 1910.3, 1877.2, 1866.6, 1862.6, 1860.3, 1899.7, &
         1971.0, 1999.9, 1970.9, 1936.5, 1922.8, 1922.8, 1924.0, 1917.2, &
         1912.0, 1926.2, 1959.7, 1995.4, 1995.9, 1938.8, 1883.5, 1894.7, &
         1933.3, 1935.1, 1899.3, 1852.7, 1820.2, 1821.5, 1865.2, 1935.5, &
         1966.1, 1919.6, 1881.2, 1931.5, 2015.6, 2050.0, 2021.4, 1960.8, &
         1938.2, 1997.0, 2051.0, 2003.4, 1912.1, 1880.2, 1895.2, 1898.0, &
         1898.8, 1938.3, 1994.2, 2010.0, 1982.4, 1948.8, 1927.3, 1911.6, &
         1877.7, 1791.6, 1679.8, 1645.0, 1727.3, 1845.2, 1926.2, 1973.4, &
         2005.2, 2021.6, 2021.8, 2025.7, 2054.3, 2086.5, 2082.6, 2052.9, &
         2047.1, 2070.2, 2072.4, 2038.1, 2020.2, 2049.9, 2074.0, 2038.1/)

!    solar spectrum from  20980 to  23840 cm-1,  in steps of 20 cm-
    solara(1050:1193)=(/ &
         1978.6, 1963.5, 1996.8, 2037.5, 2057.5, 2048.2, 2018.4, 1999.2, &
         2011.4, 2039.5, 2056.0, 2040.2, 1981.8, 1911.4, 1891.8, 1938.3, & 
         1991.7, 2005.5, 2000.8, 2011.3, 2022.7, 1997.5, 1947.7, 1936.3, &
         1986.6, 2037.9, 2032.8, 1995.7, 1984.0, 2012.0, 2055.5, 2091.6, &
         2106.5, 2094.9, 2070.4, 2052.8, 2046.7, 2043.8, 2035.5, 2016.6, &
         1988.4, 1973.3, 1999.0, 2057.4, 2103.8, 2109.4, 2089.4, 2068.5, &
         2051.8, 2031.2, 2005.9, 1986.7, 1981.5, 1979.4, 1964.1, 1943.6, &
         1951.8, 2007.3, 2083.2, 2139.1, 2158.0, 2143.3, 2103.2, 2050.9, &
         2001.9, 1974.5, 1988.0, 2037.8, 2075.1, 2050.6, 1971.5, 1884.5, &
         1828.5, 1820.9, 1866.4, 1935.3, 1974.2, 1958.7, 1925.1, 1920.2, &
         1949.7, 1984.6, 1996.4, 1966.4, 1884.8, 1781.9, 1726.8, 1759.4, &
         1817.4, 1800.4, 1692.6, 1593.2, 1598.6, 1700.3, 1823.8, 1909.7, &
         1937.7, 1902.5, 1822.4, 1737.8, 1683.2, 1666.8, 1682.7, 1715.3, &
         1734.1, 1712.4, 1668.2, 1655.0, 1698.1, 1727.2, 1636.9, 1415.7, &
         1204.2, 1155.8, 1278.4, 1450.0, 1560.5, 1595.1, 1587.8, 1570.6, &
         1565.8, 1590.3, 1640.5, 1688.4, 1708.1, 1703.6, 1700.7, 1718.5, &
         1749.0, 1772.2, 1772.5, 1745.2, 1690.2, 1624.9, 1589.0, 1618.5, &
         1701.3, 1783.2, 1816.4, 1800.7, 1765.0, 1734.1, 1714.6, 1705.0/)

!    solar spectrum from  23860 to  26720 cm-1,  in steps of 20 cm-
    solara(1194:1337)=(/ &
         1701.6, 1696.6, 1682.0, 1661.4, 1657.2, 1693.0, 1763.2, 1826.5, &
         1841.6, 1806.1, 1755.6, 1725.8, 1724.2, 1736.8, 1749.0, 1756.1, &
         1759.5, 1762.1, 1770.2, 1791.7, 1826.8, 1848.9, 1819.6, 1720.7, &
         1595.5, 1513.9, 1522.5, 1602.0, 1706.2, 1793.4, 1837.9, 1820.3, &
         1738.3, 1631.1, 1553.1, 1539.2, 1574.3, 1623.9, 1660.6, 1676.8, &
         1673.1, 1652.9, 1626.4, 1606.7, 1604.2, 1620.9, 1654.5, 1701.2, &
         1752.2, 1796.2, 1822.8, 1827.4, 1808.5, 1767.0, 1713.9, 1667.3, &
         1643.7, 1643.5, 1652.5, 1655.3, 1638.7, 1592.2, 1506.4, 1377.3, &
         1209.5, 1010.5, 807.59, 666.84, 664.53, 835.23, 1099.6, 1330.7, &
         1423.2, 1363.7, 1194.1, 961.77, 725.04, 551.29, 504.01, 596.30, &
         775.15, 975.62, 1150.2, 1287.2, 1386.1, 1447.5, 1473.7, 1468.5, &
         1435.2, 1376.9, 1296.0, 1195.5, 1085.3, 985.40, 917.25, 894.59, &
         910.86, 951.53, 1001.7, 1046.4, 1070.7, 1061.2, 1021.2, 977.16, &
         959.15, 982.06, 1020.5, 1032.6, 983.44, 879.83, 762.66, 675.28, &
         643.33, 662.65, 721.49, 808.35, 913.24, 1027.0, 1139.9, 1236.2, &
         1293.2, 1287.1, 1210.4, 1102.1, 1021.6, 1022.8, 1109.3, 1232.6, &
         1337.0, 1383.1, 1372.8, 1324.7, 1257.7, 1188.8, 1133.5, 1106.5, &
         1113.7, 1136.8, 1147.9, 1121.4, 1054.1, 968.10, 889.19, 837.87/)

!    solar spectrum from  26740 to  28780 cm-1,  in steps of 20 cm-
    solara(1338:1440)=(/ &
         817.64, 823.72, 851.04, 896.53, 959.85, 1041.2, 1137.6, 1231.2, &
         1294.4, 1299.9, 1241.2, 1155.0, 1092.0, 1097.1, 1170.2, 1263.5, &
         1322.4, 1307.4, 1233.6, 1146.1, 1090.8, 1092.5, 1134.6, 1188.9, &
         1228.9, 1245.5, 1248.5, 1250.3, 1260.5, 1274.6, 1279.5, 1261.8, &
         1214.3, 1145.4, 1069.6, 1001.4, 952.52, 930.48, 941.68, 990.34, &
         1064.4, 1135.2, 1171.5, 1149.1, 1076.3, 984.35, 906.25, 868.17, &
         873.75, 915.33, 984.41, 1067.2, 1137.1, 1163.1, 1115.5, 990.55, &
         830.93, 692.29, 627.44, 654.10, 739.24, 838.88, 911.69, 941.90, &
         944.42, 939.58, 946.10, 970.23, 1005.2, 1042.4, 1073.8, 1097.0, &
         1114.3, 1128.8, 1142.9, 1153.4, 1152.4, 1131.5, 1084.2, 1016.7, &
         945.95, 890.37, 866.15, 876.54, 913.13, 966.10, 1025.4, 1080.2, &
         1119.0, 1102.7, 1243.5, 1209.9, 1079.2, 852.20, 956.80, 842.31, &
         897.44, 1081.8, 914.23, 993.09, 1049.8, 844.95, 839.16/)
     
!    solar spectrum from  28400 to  29830 cm-1,  in steps of 10 cm-
    solarb(1:144)=(/ &
         876.54, 892.17, 913.13, 938.18, 966.10, 995.62, 1025.4, 1054.1, &
         1080.2, 1102.1, 1119.0, 1132.2, 1102.7, 1159.3, 1243.5, 1238.3, &
         1209.9, 1196.2, 1079.2, 895.60, 852.20, 935.59, 956.80, 897.09, &
         842.31, 821.15, 897.44, 1042.7, 1081.8, 988.79, 914.23, 929.38, &
         993.09, 1041.9, 1049.8, 984.33, 844.95, 770.76, 839.16, 939.65, &
         1026.1, 1121.1, 1162.6, 1142.6, 1077.9, 1027.3, 1078.2, 1094.3, &
         969.83, 853.72, 849.91, 909.12, 995.68, 1095.0, 1146.9, 1086.3, &
         1010.4, 1065.4, 1128.9, 1080.6, 987.93, 898.18, 835.20, 771.63, &
         687.12, 614.52, 606.14, 737.09, 908.13, 997.64, 1080.6, 1126.3, &
         1056.7, 1028.4, 1141.7, 1252.6, 1225.3, 1103.2, 1038.6, 1043.4, &
         1002.9, 965.51, 1035.0, 1150.7, 1200.9, 1152.0, 1068.5, 995.84, &
         889.52, 818.48, 907.01, 1042.2, 1055.6, 1000.6, 972.00, 985.72, &
         1027.2, 1054.8, 1078.0, 1126.6, 1205.3, 1245.7, 1201.0, 1144.7, &
         1097.5, 1030.1, 926.85, 836.71, 864.11, 993.50, 1075.3, 1032.6, &
         1008.9, 1066.1, 1067.4, 1004.8, 971.54, 923.18, 815.71, 799.70, &
         946.19, 1100.1, 1126.4, 1032.2, 895.14, 784.30, 734.77, 726.53, &
         726.88, 765.54, 863.90, 992.24, 1070.9, 1028.1, 858.78, 647.15, &
         563.18, 679.98, 906.40, 1094.3, 1155.3, 1124.3, 1098.4, 1109.5/)

!    solar spectrum from  29840 to  31270 cm-1,  in steps of 10 cm-
    solarb(145:288)=(/ &
         1076.2, 944.17, 849.20, 928.54, 1062.0, 1118.9, 1119.2, 1074.6, &
         1005.8, 980.02, 999.11, 1002.4, 939.78, 838.12, 816.13, 908.73, &
         1014.9, 1058.3, 1043.7, 987.54, 946.35, 981.40, 1055.8, 1094.3, &
         1028.3, 916.41, 908.99, 991.83, 1049.6, 1076.2, 1093.5, 1076.3, &
         1014.5, 949.61, 947.26, 1001.2, 1051.5, 1072.8, 1068.0, 1012.5, &
         907.81, 866.30, 950.89, 1037.5, 1079.5, 1183.9, 1291.3, 1268.6, &
         1199.3, 1188.6, 1188.0, 1186.6, 1198.2, 1171.3, 1132.6, 1131.6, &
         1096.0, 971.10, 847.07, 836.62, 922.78, 990.99, 987.51, 969.24, &
         981.46, 981.36, 971.95, 985.34, 1003.0, 1037.2, 1071.2, 1065.7, &
         1026.7, 984.84, 1002.7, 1070.3, 1117.5, 1116.0, 1048.9, 965.34, &
         972.27, 1045.7, 1096.6, 1127.5, 1133.5, 1099.6, 1079.3, 1082.9, &
         1026.8, 927.50, 879.08, 858.83, 831.01, 807.82, 789.56, 813.75, &
         893.46, 937.62, 901.56, 864.46, 873.35, 891.03, 862.46, 810.30, &
         787.36, 752.93, 715.34, 708.07, 728.93, 786.79, 807.73, 736.28, &
         645.08, 616.90, 649.17, 691.77, 749.18, 820.21, 820.68, 791.26, &
         854.27, 940.56, 956.38, 909.42, 824.18, 767.17, 722.06, 653.42, &
         624.67, 633.73, 655.14, 707.93, 784.94, 880.79, 961.15, 985.60, &
         986.18, 966.53, 921.47, 888.89, 855.85, 851.66, 886.78, 850.97/)

!    solar spectrum from  31280 to  32710 cm-1,  in steps of 10 cm-
    solarb(289:432)=(/ &
         766.97, 738.95, 724.53, 657.61, 587.77, 616.86, 760.61, 903.23, &
         917.27, 838.49, 784.80, 759.41, 719.61, 671.48, 624.63, 588.57, &
         574.70, 596.68, 698.02, 866.39, 974.82, 960.37, 930.10, 962.65, &
         1007.1, 1001.9, 926.29, 816.64, 763.25, 772.93, 762.66, 729.39, &
         725.01, 727.16, 672.73, 581.42, 520.97, 488.80, 478.60, 542.08, &
         663.71, 749.48, 785.87, 811.05, 818.19, 813.80, 824.54, 836.62, &
         799.66, 728.00, 660.36, 559.28, 473.28, 550.16, 752.04, 885.84, &
         906.80, 912.21, 929.32, 899.72, 830.20, 774.56, 736.42, 724.09, &
         740.12, 754.11, 764.96, 780.76, 788.94, 784.87, 758.80, 725.91, &
         751.84, 804.24, 777.73, 703.36, 665.27, 663.99, 679.36, 706.09, &
         757.57, 836.09, 880.02, 881.18, 907.91, 929.26, 894.32, 874.01, &
         918.56, 953.50, 922.32, 866.61, 836.54, 825.28, 752.54, 586.02, &
         427.46, 374.05, 437.23, 534.32, 556.74, 563.11, 629.31, 631.26, &
         518.76, 438.31, 460.31, 530.45, 608.50, 657.99, 662.08, 686.17, &
         775.18, 843.11, 797.46, 685.33, 611.33, 628.74, 711.36, 754.94, &
         728.80, 722.79, 726.38, 679.68, 665.83, 710.48, 723.10, 724.09, &
         760.18, 784.01, 742.78, 634.33, 546.55, 563.54, 611.03, 623.16, &
         665.36, 743.55, 764.46, 671.14, 513.18, 401.86, 405.77, 515.72/)

!    solar spectrum from  32720 to  34150 cm-1,  in steps of 10 cm-
    solarb(433:576)=(/ &
         639.90, 677.85, 679.55, 759.33, 848.11, 819.89, 751.75, 710.50, & 
         615.33, 525.09, 583.35, 715.23, 767.53, 739.10, 664.05, 580.57, &
         572.85, 634.13, 648.77, 561.27, 497.72, 591.71, 737.83, 794.19, &
         802.51, 799.33, 735.79, 658.41, 659.47, 718.18, 761.67, 697.24, &
         545.14, 474.47, 526.96, 597.65, 584.74, 447.28, 291.35, 261.28, &
         330.26, 401.96, 466.32, 531.26, 572.34, 584.86, 585.17, 569.46, &
         558.27, 559.41, 512.02, 426.37, 378.14, 398.26, 473.49, 542.18, &
         531.76, 437.48, 341.85, 305.82, 299.88, 328.12, 440.04, 586.46, &
         660.32, 625.22, 510.26, 418.85, 447.36, 534.89, 605.86, 667.07, &
         687.31, 636.79, 549.63, 472.88, 419.53, 370.06, 327.98, 320.49, &
         354.00, 399.17, 450.98, 528.34, 608.25, 696.07, 774.28, 760.75, &
         690.58, 648.20, 580.63, 477.96, 453.91, 488.74, 464.02, 421.59, &
         444.32, 446.59, 375.95, 342.13, 397.49, 510.97, 646.38, 725.14, &
         703.06, 639.06, 619.10, 654.66, 665.99, 611.40, 580.22, 607.29, &
         591.05, 542.30, 583.82, 673.02, 673.21, 582.44, 465.73, 377.25, &
         377.04, 487.27, 607.93, 617.52, 583.46, 601.68, 615.94, 575.47, &
         541.63, 542.06, 522.28, 472.49, 423.29, 438.09, 556.72, 664.34, &
         669.88, 657.45, 684.71, 705.70, 683.11, 600.81, 509.90, 497.64/)

!    solar spectrum from  34160 to  35590 cm-1,  in steps of 10 cm-
    solarb(577:720)=(/ &
         511.07, 496.07, 500.32, 518.70, 529.91, 563.00, 609.20, 626.49, &
         622.11, 615.72, 600.44, 591.26, 598.12, 593.07, 590.94, 631.58, &
         696.48, 718.48, 676.11, 631.56, 619.64, 620.53, 624.10, 636.56, &
         658.02, 688.78, 724.81, 742.60, 722.31, 675.86, 665.96, 704.73, &
         703.70, 645.00, 598.26, 587.77, 590.94, 575.93, 528.03, 477.92, &
         457.52, 456.80, 454.91, 448.65, 445.47, 445.38, 444.43, 446.04, &
         455.91, 468.02, 454.34, 393.32, 301.22, 211.44, 167.11, 193.99, &
         254.01, 305.35, 353.03, 385.08, 387.03, 391.60, 406.20, 415.34, &
         435.34, 469.77, 492.15, 472.73, 409.86, 353.25, 340.68, 355.27, &
         379.77, 401.81, 409.67, 406.89, 393.16, 378.89, 375.20, 373.52, &
         360.19, 322.79, 273.55, 237.76, 212.33, 184.80, 156.20, 127.75, &
         96.269, 68.806, 62.047, 77.143, 100.47, 127.56, 159.88, 194.05, &
         225.20, 254.64, 285.75, 300.14, 294.40, 308.92, 340.83, 346.26, &
         336.29, 347.54, 373.81, 388.78, 372.68, 325.29, 294.40, 317.56, &
         360.30, 378.08, 374.22, 374.03, 383.34, 387.88, 377.55, 356.96, &
         340.67, 328.71, 314.00, 316.91, 344.51, 355.54, 335.66, 318.68, &
         318.65, 322.43, 318.61, 304.92, 284.84, 268.13, 265.80, 273.55, &
         274.18, 252.38, 215.04, 188.60, 181.31, 181.31, 180.78, 175.24 /)

!    solar spectrum from  35600 to  37030 cm-1,  in steps of 10 cm-
    solarb(721:864)=(/ &
         162.06, 145.08, 128.76, 113.76, 98.078, 83.072, 76.222, 78.359, &
         78.434, 74.235, 75.843, 80.321, 77.859, 70.298, 64.651, 67.049, &
         77.810, 83.167, 75.286, 71.202, 80.549, 92.008, 100.17, 108.63, &
         119.44, 130.78, 142.31, 158.94, 177.12, 186.40, 186.60, 181.47, &
         175.30, 175.54, 179.00, 177.04, 172.60, 172.67, 178.98, 193.77, &
         215.13, 233.62, 252.05, 277.68, 298.91, 298.40, 280.81, 274.21, &
         286.52, 285.46, 259.71, 241.39, 246.98, 259.87, 274.27, 298.47, &
         316.85, 303.19, 263.69, 229.31, 227.90, 256.12, 281.58, 300.19, &
         310.56, 279.54, 211.93, 152.18, 129.94, 147.47, 181.62, 215.37, &
         239.50, 233.12, 191.55, 139.41, 110.51, 118.93, 134.79, 129.05, &
         124.39, 143.53, 158.29, 141.84, 116.32, 111.59, 128.93, 149.17, &
         153.44, 145.63, 148.52, 159.25, 155.84, 154.17, 177.28, 203.40, &
         207.35, 205.27, 222.85, 253.18, 271.28, 279.27, 302.17, 321.47, &
         288.83, 230.14, 206.40, 213.22, 216.49, 207.46, 196.20, 195.21, &
         202.03, 194.33, 164.86, 136.65, 123.87, 128.14, 161.89, 216.99, &
         253.68, 249.26, 222.89, 213.11, 243.64, 293.10, 309.42, 286.40, &
         269.61, 272.23, 271.67, 265.84, 265.61, 264.77, 266.03, 289.51, &
         325.67, 337.34, 321.17, 300.30, 282.60, 287.14, 322.06, 335.79 /)

!    solar spectrum from  37040 to  38470 cm-1,  in steps of 10 cm-
    solarb(865:1008)=(/ &
         297.22, 254.10, 243.47, 239.49, 219.32, 211.94, 239.28, 271.43,  &
         279.37, 272.26, 264.77, 250.52, 229.93, 222.15, 235.30, 256.79,  &
         275.28, 286.92, 284.85, 269.52, 255.05, 253.46, 263.22, 274.78,  &
         279.19, 270.17, 249.41, 229.04, 221.64, 231.38, 252.70, 280.64,  &
         310.06, 328.33, 325.01, 290.26, 238.97, 223.38, 257.24, 282.60,  &
         264.32, 243.34, 253.18, 272.89, 271.32, 256.12, 260.24, 271.35,  &
         257.11, 236.61, 238.72, 248.92, 255.90, 272.04, 291.78, 297.40,  &
         288.09, 283.28, 292.92, 301.74, 309.07, 322.05, 320.42, 295.43,  &
         269.65, 254.41, 240.88, 228.18, 221.23, 213.72, 201.23, 197.17,  &
         212.29, 233.39, 247.65, 261.74, 286.17, 322.49, 349.47, 338.28,  &
         297.06, 261.55, 252.28, 264.65, 286.92, 298.94, 280.45, 244.37,  &
         213.47, 193.03, 182.07, 168.54, 143.12, 114.10, 89.615, 73.589,  &
         73.990, 87.912, 96.265, 94.813, 96.604, 102.30, 102.15, 103.07,  &
         117.81, 137.41, 146.09, 144.28, 137.89, 128.11, 122.82, 128.19,  &
         130.66, 117.31, 98.912, 93.397, 105.63, 122.73, 126.39, 113.05,  &
         92.317, 76.340, 69.032, 66.324, 71.280, 87.431, 105.94, 114.02,  &
         107.91, 91.872, 75.208, 69.123, 75.930, 90.928, 109.71, 125.70,  &
         135.79, 141.14, 138.14, 121.33, 91.806, 63.497, 52.106, 59.555  /)

!    solar spectrum from  38480 to  39910 cm-1,  in steps of 10 cm-
    solarb(1009:1152)=(/ &
         81.015, 106.67, 118.97, 116.36, 110.82, 100.88, 89.056, 90.431, &
         104.41, 114.95, 124.85, 148.87, 171.72, 167.22, 142.25, 118.42, &
         98.653, 78.908, 68.133, 77.286, 100.93, 120.08, 125.49, 131.79, &
         155.69, 180.75, 181.81, 166.77, 150.06, 133.24, 116.14, 97.728, &
         81.629, 76.695, 87.607, 110.23, 134.88, 149.13, 147.64, 139.88, &
         135.19, 135.07, 138.00, 136.73, 128.84, 122.22, 120.48, 121.98, &
         123.08, 116.30, 101.43, 86.303, 74.719, 68.800, 71.327, 80.626, &
         90.485, 96.739, 100.69, 100.81, 93.677, 84.740, 81.532, 82.893, &
         84.564, 87.584, 91.780, 91.272, 87.014, 87.386, 90.149, 84.917, &
         71.266, 57.873, 51.863, 53.876, 57.909, 58.508, 57.020, 57.432, &
         60.671, 64.667, 67.362, 67.511, 64.233, 59.035, 55.697, 56.636, &
         59.400, 59.070, 56.522, 55.834, 55.860, 54.039, 51.976, 52.344, &
         54.667, 56.450, 56.751, 56.769, 58.002, 60.029, 59.602, 53.134, &
         42.926, 35.588, 33.447, 35.171, 39.379, 44.371, 47.745, 46.933, &
         42.441, 37.879, 35.595, 36.458, 41.048, 47.300, 51.098, 50.024, &
         45.331, 41.282, 40.082, 40.000, 39.104, 37.329, 36.632, 37.792, &
         39.189, 41.058, 45.214, 50.737, 54.281, 55.015, 56.138, 60.931, &
         67.383, 69.534, 65.159, 56.372, 47.326, 44.322, 49.944, 59.696 /)

!    solar spectrum from  39920 to  41350 cm-1,  in steps of 10 cm-
    solarb(1153:1296)=(/ &
         67.929, 71.334, 69.905, 65.620, 59.303, 54.016, 55.880, 65.155, &
         74.065, 76.217, 73.506, 71.406, 70.849, 69.749, 69.268, 71.380, &
         72.721, 68.929, 61.665, 54.896, 47.420, 38.325, 32.219, 31.243, &
         33.310, 35.358, 35.623, 36.840, 41.551, 47.499, 51.176, 50.344, &
         45.362, 38.341, 33.130, 33.801, 40.140, 49.121, 55.385, 55.174, &
         50.450, 46.511, 47.495, 51.883, 56.354, 59.603, 61.584, 63.215, &
         64.603, 64.101, 59.027, 50.956, 47.633, 52.543, 58.883, 59.829, &
         57.617, 56.727, 57.371, 57.898, 57.177, 55.129, 52.952, 52.018, &
         52.186, 52.044, 50.269, 46.592, 42.515, 40.755, 41.887, 44.119, &
         46.536, 48.858, 50.490, 51.919, 54.085, 54.707, 51.927, 49.449, &
         49.865, 50.933, 50.496, 48.616, 46.717, 46.070, 46.263, 46.733, &
         48.009, 50.187, 52.420, 53.536, 52.507, 51.380, 53.214, 56.985, &
         60.614, 63.139, 63.999, 63.869, 65.100, 69.385, 74.743, 78.184, &
         78.103, 74.113, 67.371, 60.849, 58.924, 62.682, 68.032, 69.117, &
         64.604, 59.110, 55.998, 56.838, 61.778, 65.874, 65.079, 63.038, &
         64.809, 69.911, 74.841, 76.439, 73.587, 68.853, 67.497, 72.675, &
         80.602, 83.422, 78.957, 72.228, 66.737, 62.842, 61.535, 63.574, &
         69.248, 76.577, 79.922, 77.755, 73.938, 70.518, 68.003, 66.339 /)

!    solar spectrum from  41360 to  42790 cm-1,  in steps of 10 cm-
    solarb(1297:1440)=(/ &
         63.979, 61.098, 59.421, 58.103, 55.741, 52.549, 48.079, 42.578, &
         38.373, 37.297, 37.455, 34.861, 30.483, 29.634, 34.734, 42.460, &
         47.066, 45.848, 40.157, 34.290, 31.584, 30.650, 29.054, 27.788, &
         30.427, 37.570, 44.196, 46.880, 47.848, 49.166, 49.180, 45.002, &
         38.135, 35.055, 38.095, 41.750, 40.899, 35.722, 28.884, 24.835, &
         28.670, 39.646, 50.310, 55.725, 57.401, 58.110, 59.406, 59.360, &
         53.420, 43.004, 34.787, 33.697, 39.682, 47.554, 52.605, 53.632, &
         51.001, 45.266, 37.844, 31.030, 25.936, 22.799, 21.882, 23.484, &
         27.857, 33.447, 37.319, 39.195, 42.826, 50.398, 58.752, 63.301, &
         61.094, 53.532, 46.046, 41.118, 37.646, 36.304, 40.426, 50.893, &
         61.553, 65.395, 62.680, 58.087, 54.622, 51.330, 46.874, 42.870, &
         40.547, 39.760, 40.217, 40.359, 39.559, 40.667, 46.260, 53.413, &
         56.041, 52.566, 46.674, 41.073, 35.511, 31.231, 31.082, 35.955, &
         45.199, 55.464, 61.802, 63.505, 61.850, 56.412, 49.388, 46.369, &
         50.058, 56.694, 60.884, 61.030, 58.107, 54.303, 51.940, 50.508, &
         46.749, 39.155, 31.535, 28.959, 30.973, 32.670, 31.567, 29.340, &
         27.275, 25.184, 24.264, 27.068, 34.296, 42.475, 47.230, 47.425, &
         44.435, 40.538, 36.868, 33.020, 29.405, 28.753, 34.079, 44.246 /)

!    solar spectrum from  42800 to  44230 cm-1,  in steps of 10 cm-
    solarb(1441:1584)=(/ &
         53.780, 57.974, 56.376, 51.200, 45.308, 40.273, 35.900, 33.344, &
         34.011, 36.858, 41.283, 47.374, 53.088, 56.201, 55.633, 50.843, &
         43.997, 38.767, 36.248, 36.380, 40.762, 50.700, 63.371, 73.432, &
         76.418, 70.373, 58.741, 47.034, 38.598, 34.664, 35.794, 42.084, &
         49.973, 54.338, 53.956, 52.287, 52.778, 55.571, 59.034, 60.268, &
         56.247, 47.362, 38.056, 32.889, 31.739, 31.734, 32.476, 35.060, &
         39.091, 43.398, 48.131, 53.574, 58.749, 63.599, 68.971, 73.421, &
         73.861, 69.003, 60.557, 51.865, 44.879, 42.060, 44.802, 47.950, &
         46.882, 42.973, 39.293, 37.711, 37.137, 35.222, 32.243, 30.488, &
         32.605, 40.429, 51.099, 57.710, 57.150, 52.992, 50.275, 49.986, &
         49.778, 48.371, 46.421, 44.604, 42.730, 41.244, 41.565, 43.805, &
         47.013, 48.992, 46.428, 40.595, 37.840, 42.353, 52.248, 60.529, &
         61.566, 56.800, 52.041, 52.260, 57.077, 61.019, 60.712, 57.048, &
         51.481, 46.352, 44.366, 44.947, 45.478, 44.944, 43.825, 42.105, &
         39.466, 36.826, 35.907, 36.357, 35.661, 33.947, 33.690, 34.429, &
         34.000, 32.645, 31.410, 30.281, 29.409, 29.127, 29.326, 29.869, &
         30.601, 31.311, 32.099, 32.779, 32.757, 32.098, 31.975, 33.484, &
         36.048, 39.169, 43.365, 47.244, 48.214, 45.786, 41.586, 38.775 /)

!    solar spectrum from  44240 to  45670 cm-1,  in steps of 10 cm-
    solarb(1585:1728)=(/ &
         40.753, 46.752, 51.684, 52.597, 51.449, 50.684, 49.450, 46.747, &
         45.369, 47.685, 50.240, 48.961, 46.693, 48.600, 53.694, 56.465, &
         54.341, 50.722, 49.877, 51.246, 52.088, 52.765, 56.254, 63.326, &
         69.744, 71.066, 68.349, 65.123, 62.551, 59.195, 53.705, 48.161, &
         46.236, 47.710, 49.660, 50.799, 51.836, 54.537, 59.647, 64.707, &
         65.844, 61.634, 55.570, 54.083, 58.781, 64.888, 69.777, 74.008, &
         76.492, 76.226, 74.746, 74.941, 77.801, 79.619, 76.190, 67.190, &
         55.231, 45.813, 43.141, 45.647, 49.466, 52.231, 52.221, 48.886, &
         44.716, 42.613, 43.385, 45.968, 48.121, 48.998, 49.885, 50.707, &
         49.893, 48.319, 48.198, 50.280, 53.830, 55.914, 54.822, 52.939, &
         51.944, 49.438, 42.956, 34.614, 28.100, 24.503, 24.203, 27.839, &
         34.604, 41.615, 45.324, 45.444, 45.527, 47.179, 45.756, 36.862, &
         26.037, 20.569, 20.329, 24.263, 30.863, 35.939, 36.711, 35.693, &
         37.256, 40.862, 44.416, 48.800, 54.182, 57.655, 58.427, 59.965, &
         63.940, 66.820, 65.465, 59.482, 49.396, 39.422, 34.182, 35.388, &
         42.875, 52.034, 57.595, 59.093, 57.272, 52.172, 45.493, 39.419, &
         35.581, 35.902, 40.354, 46.732, 53.309, 58.781, 61.785, 59.255, &
         50.030, 41.567, 40.523, 43.584, 44.875, 42.754, 40.077, 39.941 /)

!    solar spectrum from  45680 to  47110 cm-1,  in steps of 10 cm-
    solarb(1729:1872)=(/ &
         40.977, 39.567, 34.955, 30.424, 31.039, 38.687, 47.480, 49.830, &
         46.790, 44.829, 46.546, 50.415, 54.602, 57.656, 58.463, 57.276, &
         55.621, 54.514, 53.338, 50.026, 42.817, 33.636, 27.134, 25.516, &
         27.897, 31.392, 32.125, 29.463, 26.581, 25.956, 27.737, 31.175, &
         34.959, 37.671, 38.641, 37.958, 36.733, 35.681, 33.877, 30.849, &
         28.059, 27.615, 29.319, 29.375, 25.390, 20.659, 19.484, 22.297, &
         27.282, 32.467, 35.906, 37.137, 37.895, 39.130, 39.777, 39.872, &
         40.778, 42.317, 42.934, 40.430, 34.227, 27.701, 23.880, 22.174, &
         21.639, 22.589, 25.184, 29.017, 32.981, 36.110, 38.580, 41.239, &
         44.426, 46.939, 47.010, 44.165, 39.659, 35.556, 32.838, 31.546, &
         32.676, 36.963, 42.333, 44.931, 43.704, 40.943, 37.973, 35.199, &
         33.574, 33.339, 34.185, 36.347, 39.963, 43.964, 47.162, 48.987, &
         48.976, 47.948, 48.004, 49.892, 51.065, 47.834, 40.489, 32.665, &
         26.795, 24.461, 26.655, 31.928, 37.634, 41.345, 40.956, 36.827, &
         32.110, 28.612, 26.482, 26.602, 28.831, 30.877, 30.976, 30.063, &
         29.887, 30.305, 29.974, 28.265, 26.517, 27.066, 30.403, 34.539, &
         37.104, 37.598, 37.252, 37.060, 36.498, 34.167, 29.814, 24.192, &
         18.515, 15.086, 15.040, 17.158, 20.807, 25.682, 30.352, 34.203 /)

!    solar spectrum from  47120 to  48550 cm-1,  in steps of 10 cm-
    solarb(1873:2016)=(/ &
         37.902, 42.531, 47.832, 50.509, 48.019, 42.616, 38.321, 37.370, &
         40.172, 44.395, 46.132, 43.911, 38.396, 31.379, 26.275, 25.075, &
         26.652, 28.963, 31.168, 34.168, 38.050, 40.231, 38.347, 32.741, &
         26.199, 21.863, 20.249, 20.185, 21.726, 25.562, 30.318, 33.431, &
         34.453, 34.959, 36.374, 37.870, 36.655, 31.966, 25.920, 21.264, &
         20.663, 24.658, 30.263, 34.021, 34.336, 31.356, 26.926, 23.109, &
         20.867, 20.684, 22.416, 24.878, 26.779, 27.334, 26.537, 25.210, &
         24.013, 22.944, 21.800, 20.449, 19.290, 19.528, 21.742, 24.125, &
         23.994, 21.559, 19.555, 18.915, 18.342, 17.335, 16.549, 16.479, &
         17.211, 18.445, 19.294, 18.980, 17.912, 17.156, 17.103, 17.256, &
         16.925, 15.842, 14.485, 13.683, 13.647, 13.914, 14.009, 13.770, &
         13.456, 13.399, 13.547, 13.760, 14.060, 14.427, 14.644, 14.438, &
         13.986, 13.749, 13.927, 14.390, 14.759, 14.822, 14.679, 14.448, &
         14.186, 13.937, 13.754, 13.657, 13.540, 13.308, 13.053, 12.841, &
         12.704, 12.742, 12.811, 12.662, 12.355, 12.100, 12.003, 12.014, &
         12.067, 12.223, 12.444, 12.472, 12.164, 11.732, 11.515, 11.619, &
         11.873, 12.028, 11.947, 11.722, 11.399, 10.930, 10.473, 10.205, &
         10.224, 10.694, 11.468, 12.007, 12.083, 11.905, 11.498, 10.891 /)

!   solar spectrum from  48560 to  49990 cm-1,  in steps of 10 cm-
    solarb(2017:2160)=(/ &
         10.575, 10.846, 11.353, 11.612, 11.411, 10.876, 10.383, 10.305, &
         10.695, 11.245, 11.636, 11.828, 11.918, 11.865, 11.674, 11.510, &
         11.407, 11.303, 11.216, 11.143, 11.039, 10.983, 11.004, 10.900, &
         10.653, 10.562, 10.781, 11.186, 11.605, 11.806, 11.582, 11.056, &
         10.567, 10.335, 10.408, 10.729, 11.165, 11.540, 11.646, 11.372, &
         10.933, 10.524, 9.9973, 9.3783, 8.9883, 9.0163, 9.4125, 9.9179, &
         10.278, 10.472, 10.553, 10.575, 10.519, 10.216, 9.6821, 9.1499, &
         8.7057, 8.3894, 8.3442, 8.6241, 9.1371, 9.7184, 10.191, 10.443, &
         10.458, 10.289, 9.9772, 9.5829, 9.3097, 9.3195, 9.4694, 9.5182, &
         9.4326, 9.2478, 8.8197, 7.9809, 6.9996, 6.4856, 6.7462, 7.5406, &
         8.2813, 8.7258, 9.0682, 9.1665, 8.8637, 8.4638, 8.2393, 8.1656, &
         8.1880, 8.3578, 8.6488, 8.8980, 9.0117, 9.0659, 9.1955, 9.4207, &
         9.5526, 9.4237, 9.1290, 8.8441, 8.6138, 8.4237, 8.2979, 8.2598, &
         8.2859, 8.3475, 8.4533, 8.6285, 8.8310, 8.8866, 8.6750, 8.3312, &
         8.0091, 7.7296, 7.6239, 7.8692, 8.2725, 8.4086, 8.2515, 8.0914, &
         8.0003, 7.9367, 7.9266, 7.9580, 8.0492, 8.2376, 8.4263, 8.4811, &
         8.3309, 8.0263, 7.7632, 7.6987, 7.8124, 7.9390, 8.0183, 8.0816, &
         8.0428, 7.8923, 7.6963, 7.4969, 7.4013, 7.4289, 7.4489, 7.4059 /)

!    solar spectrum from  50000 to  51430 cm-1,  in steps of 10 cm-
    solarb(2161:2304)=(/ &
         7.4198, 7.5261, 7.5252, 7.3239, 7.1263, 7.1423, 7.3340, 7.5049, &
         7.5484, 7.5319, 7.5163, 7.4995, 7.5728, 7.8104, 8.0588, 8.0948, &
         7.9140, 7.6978, 7.5116, 7.2138, 6.8063, 6.5430, 6.5232, 6.5869, &
         6.5610, 6.3984, 6.1889, 6.0587, 6.0676, 6.1988, 6.3140, 6.2527, &
         6.0929, 6.0277, 6.0941, 6.3031, 6.6594, 6.9398, 6.9566, 6.8310, &
         6.7374, 6.6812, 6.6558, 6.8336, 7.2020, 7.4012, 7.2950, 7.0488, &
         6.7966, 6.6293, 6.5868, 6.5980, 6.6007, 6.6501, 6.7627, 6.7853, &
         6.6321, 6.4856, 6.5198, 6.6486, 6.7271, 6.7227, 6.6696, 6.6189, &
         6.5979, 6.6188, 6.7110, 6.8343, 6.8750, 6.8250, 6.7885, 6.8266, &
         6.8556, 6.8068, 6.8377, 7.0467, 7.2779, 7.4139, 7.4712, 7.4621, &
         7.4071, 7.3592, 7.3372, 7.3220, 7.2938, 7.2531, 7.2052, 7.1335, &
         7.0298, 6.8533, 6.5535, 6.2227, 6.0139, 5.9384, 5.9038, 5.8568, &
         5.7909, 5.7326, 5.7745, 5.9608, 6.1865, 6.3681, 6.4997, 6.5437, &
         6.4637, 6.2708, 6.0451, 5.9557, 6.0855, 6.2542, 6.2454, 6.0795, &
         5.9102, 5.8447, 5.9218, 6.1063, 6.2895, 6.3271, 6.1097, 5.7421, &
         5.4452, 5.2981, 5.3256, 5.4935, 5.6819, 5.8245, 5.8933, 5.9630, &
         6.1703, 6.4525, 6.6325, 6.6965, 6.7185, 6.6238, 6.3107, 5.9241, &
         5.6987, 5.6651, 5.7428, 5.8790, 5.9715, 5.9618, 5.9674, 6.0754 /)

!    solar spectrum from  51440 to  52870 cm-1,  in steps of 10 cm-
    solarb(2305:2448)=(/ &
         6.2541, 6.4300, 6.4968, 6.4564, 6.4082, 6.3024, 6.0135, 5.6431, &
         5.3963, 5.2989, 5.2635, 5.2227, 5.1279, 4.9315, 4.6348, 4.3168, &
         4.0151, 3.6625, 3.2906, 3.1028, 3.1349, 3.1994, 3.2596, 3.4144, &
         3.5949, 3.6534, 3.6296, 3.6281, 3.5876, 3.4292, 3.2659, 3.2284, &
         3.2576, 3.3002, 3.4535, 3.7372, 4.0573, 4.3558, 4.5999, 4.7781, &
         4.8855, 4.8999, 4.8392, 4.7624, 4.7059, 4.6981, 4.7666, 4.8453, &
         4.8236, 4.7293, 4.6861, 4.7132, 4.7725, 4.8713, 4.9596, 4.9527, &
         4.8957, 4.9252, 5.0736, 5.2229, 5.2505, 5.1537, 5.0156, 4.8880, &
         4.7686, 4.6549, 4.5534, 4.4828, 4.4661, 4.5040, 4.5905, 4.7033, &
         4.7858, 4.8334, 4.9283, 5.0377, 5.0065, 4.8471, 4.6828, 4.5586, &
         4.4812, 4.4314, 4.3903, 4.3830, 4.4066, 4.3900, 4.2973, 4.1978, &
         4.1462, 4.1084, 4.1495, 4.3897, 4.6859, 4.8206, 4.7938, 4.6781, &
         4.5222, 4.3959, 4.3358, 4.2947, 4.2259, 4.1452, 4.1060, 4.1462, &
         4.2149, 4.2549, 4.3061, 4.3742, 4.3738, 4.2718, 4.1389, 4.0405, &
         3.9457, 3.8127, 3.7099, 3.7344, 3.8589, 3.9598, 3.9525, 3.8377, &
         3.6708, 3.5357, 3.4929, 3.5375, 3.6381, 3.7890, 3.9671, 4.0995, &
         4.1421, 4.1302, 4.1235, 4.1623, 4.2506, 4.2948, 4.2231, 4.0993, &
         3.9680, 3.9475, 4.1958, 4.5131, 4.6101, 4.5130, 4.3474, 4.1749 /)

!    solar spectrum from  52880 to  54310 cm-1,  in steps of 10 cm-
    solarb(2449:2592)=(/ &
         4.0467, 3.9956, 4.0078, 4.0374, 4.0255, 3.9379, 3.8192, 3.7529, &
         3.7675, 3.8260, 3.8654, 3.8518, 3.8148, 3.8028, 3.8098, 3.7934, &
         3.7660, 3.7944, 3.8689, 3.8978, 3.8856, 3.8923, 3.8570, 3.6940, &
         3.4693, 3.3222, 3.2824, 3.2887, 3.3039, 3.3222, 3.3313, 3.3326, &
         3.3482, 3.3807, 3.4188, 3.4602, 3.4972, 3.5151, 3.5155, 3.5165, &
         3.5258, 3.5406, 3.5478, 3.5345, 3.5339, 3.5820, 3.6396, 3.6448, &
         3.5872, 3.5112, 3.4804, 3.5257, 3.6238, 3.7290, 3.8023, 3.8024, &
         3.7268, 3.6578, 3.6439, 3.6422, 3.6373, 3.6397, 3.6410, 3.6494, &
         3.6608, 3.6251, 3.5212, 3.4020, 3.2845, 3.1230, 2.9483, 2.8515, &
         2.8432, 2.8638, 2.8967, 2.9505, 3.0025, 3.0552, 3.1106, 3.1178, &
         3.0596, 2.9854, 2.9316, 2.8903, 2.8590, 2.8500, 2.8450, 2.8121, &
         2.7626, 2.7424, 2.7667, 2.8024, 2.8165, 2.8111, 2.8128, 2.8569, &
         2.9659, 3.1062, 3.1990, 3.2128, 3.2088, 3.2391, 3.2661, 3.2364, &
         3.1173, 2.9094, 2.6952, 2.5324, 2.3959, 2.2953, 2.2510, 2.2245, &
         2.1811, 2.1301, 2.1482, 2.3257, 2.5856, 2.7226, 2.6495, 2.4508, &
         2.2444, 2.0850, 1.9891, 1.9843, 2.0816, 2.2233, 2.3248, 2.3551, &
         2.3479, 2.3606, 2.4296, 2.5361, 2.6128, 2.6216, 2.6069, 2.6196, &
         2.6464, 2.6427, 2.5823, 2.4682, 2.3320, 2.2405, 2.2637, 2.3973 /)

!    solar spectrum from  54320 to  55750 cm-1,  in steps of 10 cm-
    solarb(2593:2736)=(/ &
         2.5524, 2.6891, 2.8508, 3.0103, 3.0681, 3.0064, 2.9114, 2.8609, &
         2.8517, 2.8374, 2.7894, 2.7288, 2.7138, 2.7729, 2.8707, 2.9536, &
         2.9953, 2.9911, 2.9398, 2.8550, 2.7732, 2.7303, 2.7366, 2.7650, &
         2.7705, 2.7374, 2.6830, 2.6218, 2.5663, 2.5341, 2.5351, 2.5681, &
         2.6124, 2.6305, 2.6024, 2.5431, 2.4840, 2.4546, 2.4684, 2.5100, &
         2.5445, 2.5532, 2.5564, 2.5889, 2.6616, 2.7553, 2.8466, 2.9290, &
         2.9958, 3.0175, 2.9774, 2.8990, 2.8001, 2.6927, 2.6171, 2.5931, &
         2.5809, 2.5276, 2.4284, 2.3365, 2.3162, 2.3855, 2.4872, 2.5455, &
         2.5773, 2.6809, 2.9720, 3.5757, 4.4006, 5.0044, 5.0295, 4.5135, &
         3.7071, 2.9059, 2.3600, 2.1418, 2.1119, 2.0871, 2.0301, 2.0043, &
         2.0361, 2.0963, 2.1520, 2.1878, 2.1955, 2.1864, 2.1899, 2.2170, &
         2.2574, 2.2895, 2.2783, 2.2148, 2.1641, 2.2343, 2.4726, 2.8119, &
         3.1288, 3.2984, 3.2206, 2.8859, 2.4473, 2.1436, 2.0729, 2.1391, &
         2.2171, 2.2580, 2.2654, 2.2481, 2.2103, 2.1657, 2.1356, 2.1321, &
         2.1438, 2.1461, 2.1396, 2.1460, 2.1588, 2.1581, 2.1481, 2.1343, &
         2.1101, 2.0754, 2.0400, 2.0121, 1.9930, 1.9799, 1.9699, 1.9613, &
         1.9537, 1.9454, 1.9312, 1.9058, 1.8726, 1.8470, 1.8465, 1.8693, &
         1.8844, 1.8635, 1.8143, 1.7618, 1.7188, 1.6853, 1.6656, 1.6708 /)

!    solar spectrum from  55760 to  57190 cm-1,  in steps of 10 cm-
    solarb(2737:2880)=(/ &
         1.7036, 1.7519, 1.8120, 1.9015, 2.0124, 2.0980, 2.1385, 2.1481, &
         2.1347, 2.1086, 2.0953, 2.1062, 2.1095, 2.0685, 2.0001, 1.9461, &
         1.9194, 1.9088, 1.9023, 1.8977, 1.9049, 1.9300, 1.9588, 1.9635, &
         1.9357, 1.9019, 1.8887, 1.8939, 1.9018, 1.9038, 1.8975, 1.8747, &
         1.8289, 1.7716, 1.7303, 1.7330, 1.7900, 1.8782, 1.9548, 1.9907, &
         1.9807, 1.9430, 1.9173, 1.9218, 1.9203, 1.8717, 1.7832, 1.6965, &
         1.6389, 1.6077, 1.5924, 1.5818, 1.5583, 1.5142, 1.4616, 1.4237, &
         1.4252, 1.4834, 1.5970, 1.7410, 1.8771, 1.9784, 2.0451, 2.0872, &
         2.0909, 2.0384, 1.9573, 1.9002, 1.8824, 1.8663, 1.8193, 1.7540, &
         1.6874, 1.6222, 1.5726, 1.5450, 1.5290, 1.5312, 1.5699, 1.6411, &
         1.7186, 1.7678, 1.7546, 1.6623, 1.5115, 1.3588, 1.2605, 1.2348, &
         1.2611, 1.3091, 1.3588, 1.3884, 1.3800, 1.3482, 1.3224, 1.3159, &
         1.3437, 1.4142, 1.4950, 1.5443, 1.5521, 1.5282, 1.4902, 1.4606, &
         1.4465, 1.4398, 1.4399, 1.4544, 1.4760, 1.4781, 1.4506, 1.4229, &
         1.4185, 1.4221, 1.4119, 1.3908, 1.3779, 1.3813, 1.3933, 1.4087, &
         1.4268, 1.4417, 1.4408, 1.4188, 1.3861, 1.3548, 1.3261, 1.2980, &
         1.2769, 1.2731, 1.2856, 1.3002, 1.3056, 1.2987, 1.2817, 1.2590, &
         1.2291, 1.1868, 1.1428, 1.1183, 1.1141, 1.1120, 1.1009, 1.0797 /)

!    solar spectrum from  57200 to  57490 cm-1,  in steps of 10 cm-
    solarb(2881:2910)=(/ &
         1.0523, 1.0284, 1.0251, 1.0577, 1.1195, 1.1791, 1.2061, 1.2013, &
         1.1936, 1.2000, 1.2040, 1.1824, 1.1489, 1.1400, 1.1539, 1.1629, &
         1.1617, 1.1586, 1.1564, 1.1572, 1.1565, 1.1399, 1.1037, 1.0627, &
         1.0341, 1.0223, 1.0199, 1.0188, 1.0174, 1.0163  /)

!     if  v  is too small,  write warning  +  return sunir = 0
  if(v < 0.0) then
     sunir = 0.0
     write(6, 900) v
     return

  elseif( v >= 0.0  .and.  v < 100.0 ) then
!     for low frequencies use a power law approximation
     sunir = a*v**b
     return
     
  elseif( v >= 100.0  .and.  v <= 28420.0 ) then
!     use  4 point interpolation  on  array  solara
!     which is at  20 cm-1  spacing  from 0 to 28720 cm-1
     i = 1 + int(v/20.0)
     p = (mod(int(v),20)+(v-int(v)))/20.0
     sunir = ( w2(p)*solara(i+2) - w1(p)*solara(i+1) + &
          w0(p)*solara(i) - wm(p)*solara(i-1) ) / 6.
     return

  elseif( v >= 28420.0  .and.  v <= 57470.0 ) then
!     use  4 point interpolation  on  array  solarb
!     which is at  10 cm-1  spacing  from 28400 to 57490 cm-1
     i = int(v/10.0) - 2839
     p = (mod(int(v),10)+(v-int(v)))/10.0
     sunir = ( w2(p)*solarb(i+2) - w1(p)*solarb(i+1) + &
          w0(p)*solarb(i) - wm(p)*solarb(i-1) ) / 6.
     return
!
  elseif( v > 57470.0 ) then
!     if  v  is too large,  write warning  +  return sunir = 0
     sunir = 0.0
     write(6, 900) v
     return
  endif

  return
900 format('0 *****  warning - input frequency = ', 1pg12.5, 'cm-1', &
         /, '   outside valid range of 0 to 57470 cm-1    *******', / )
end function sunir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function wvsat(t)
! saturation water vapor amount (LOWTRAN)
!  ws=a*exp(18.9766-14.9595a-2.4388a**2) g/m3
!  a=t0/t, t0=273.15, t in k (t0-50,t0+50)
!--- history
! 93. 3. 3 created
! 08. 9.20  modified Fortran 90 free-style form
!--- input
! t       r      absolute temperature (K)
!                t0-50 to t0+50
!--- output
! wvsat  rf      water vapor amount (g/m3)
!---
  implicit none
  real(8),intent(in)::t
  real(8)::a,wvsat
!--exec
  a=273.15d0/t
  wvsat=a*dexp(18.9766d0-14.9595d0*a-2.4388d0*a**2)
  return
end function wvsat


