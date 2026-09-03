#include <define_CPP.h>
!############################################################################
! define AIRS parameters here before compiling. 
! 1 and 2 is for AIRS post-Nov2003, while 3 & 4 are for AIRS pre-Nov2003.
!
! 1 - airs_PGEv6_postNov2003:
! 2 - airs_PGEv6_postNov2003_wcon_nte:
! 3 - airs_PGEv6_preNov2003:
! 4 - airs_PGEv6_preNov2003_wcon_nte:
!
# define AIRS_PARAM 1
!############################################################################
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_sarta
 use module_simulator
 use module_mpi
 use module_surface
 use module_bf_emiss
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
!   This is wrapper module for Stand Alone RTA (SARTA) code developed by L.L.Strow/S.Hannon
!   at University of Maryland Baltimore County [UMBC].
!   Essentially, original SARTA code was F90ed and coupled with G-SDSU core module. 
!   Other than that, there are no significant change in the code.
! 
! History: 
!  11/2013  Toshi Matsui@NASA GSFC ; Initial coupling is done.
!           
! References: 
!---------------------------------------------------------------------------------------------------

 private

 public :: sarta_driver

 character(len=200),parameter :: sarta_dir_data = './../DATAFILES/' ! Directory path for the various datafiles 

 save  !all parameter is saved.

 logical,parameter :: LRHOT=.FALSE. ! use input rho for reflected thermal


#if AIRS_PARAM == 1

!=======================================================================
!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS (Atmospheric Infra-Red Sounder)
!
!    incFTC
!
!F77====================================================================


!ROUTINE NAME:
!    incFTC (include file)


!ABSTRACT:
!    Include file consisting of parameter statements to size various
!    arrays in the USEFAST related routines source code.


!CALL PROTOCOL:
!    none (include file)


!INPUT PARAMETERS:
!    none


!OUTPUT PARAMETERS:
!    none


!INPUT/OUTPUT PARAMETERS:
!    none


!RETURN VALUES:
!    none


!PARENT(S):
!    CALOKW
!    CALOWP
!    CALPAR
!    CALRAD
!    CALT1
!    CALT2
!    CALT3
!    CALT4
!    CALT5
!    CALT6
!    CALT7
!    FAKETZ
!    RDCOEF
!    RDLIST
!    RDPROF
!    RDSUN
!    SUNPAR
!    SARTA


!ROUTINES CALLED:
!    none


!FILES ACCESSED:
!    none


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    Include file for the October 2011 100 layer AIRS
!    Stand Alone RTA (SARTA) code by L.L.Strow/S.Hannon.
!
!    Parameter statements for the FTC routines.


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
! Date        Programmer     Comments
! ----------- -------------- -------------------------------------------
!  1 Dec 1994 Scott Hannon   Created
! 31 Jan 1997 Scott Hannon   Re-wrote for FWO+FOW+FMW+FCOW=Feb97 FTC
!  3 Sep 1997 Scott Hannon   Re-wrote for sets 1 - 7
! 30 Sep 1997 Scott Hannon   Added NCO2 and MXCHNC
! 26 Feb 1998 Scott Hannon   Added OPTRAN variables for water, and
!                            changed both N1H2O & N3H2O from 13 to 11
! 23 Sep 1999 Scott Hannon   Change number of channel dimensions for
!                            new Sep99 version of fast model.
!  5 Apr 2000 Scott Hannon   Added MXEMIS
!  4 Aug 2000 Scott Hannon   Changes values for use with testfast
! 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
! 23 Jan 2001 Scott Hannon   Update values of C1 & C2
! 15 Feb 2001 Scott Hannon   Add MAXPRO, CO2STD, IOERR, IOINFO,
!                            MXGAS, CSARTA, and all filenames
! 24 Apr 2001 Scott Hannon   Add MXMIEA and FNMIEA
! 14 Aug 2001 Scott Hannon   Add FNMIEE and FNMIEG
! 21 Nov 2001 Scott Hannon   Add VSARTA, VSCOEF, & VCLOUD; remove CSARTA
! 12 Sep 2002 Scott Hannon   Updated for m135f (-13.5 um with fringes)
! 17 Dec 2002 Scott Hannon   Updated for revised(Dec02) m135f
!  3 Jan 2003 Scott Hannon   Updated VSARTA for version 1.04
! 06 Feb 2004 Scott Hannon   Add FNTMLT & update VSARTA for v1.05
! 07 Apr 2005 Scott Hannon   NFCOEF increased from 5 to 6 for v1.06
! 18 May 2005 Scott Hannon   update for HNO3 version
! 29 Jun 2005 Scott Hannon   "trace" version v1.07 with CO2,SO2,HNO3,N2O
! 13 Oct 2005 Scott Hannon   Add variables for non-LTE
! 22 Nov 2005 Scott Hannon   Replace set1,set2,CO2 coefs for new M12
! 02 May 2007 Scott Hannon   Added XSALT
! 08 May 2008 Scott Hannon   Updated for v1.08; change most filenames
!                            from FN* to double filenames FA* and FB*;
!                            add YOFFA,YOFFB,YOFMIN,YOFMAX,YOFDEF
! 13 May 2008 Scott Hannon   Add CO2NTE and NTEBOT for non-LTE calc and
!                            increase NNCOEF from 6 to 7
! 12 May 2009 Scott Hannon   Add VTUNNG string; delete VCLOUD
! 11 Oct 2011 Scott Hannon   Add MAXABW; replace FNSCAL with ABFILE

!END====================================================================
!
!-----------------------------------------------------------------------
!      IMPLICIT NONE
!-----------------------------------------------------------------------
! Note: having an "implicit none" in both the include file & the main
! source code will cause some compilers to complain.
!       IMPLICIT NONE


!-----------------------------------------------------------------------
!      INCLUDE FILES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXECUTABLE CODE
!-----------------------------------------------------------------------
!      none

!      -----------------------------------------------------------------
!      Assign SARTA version strings
!      -----------------------------------------------------------------
!      The version strings consists of 3 parts: version number, date,
!      and comment.  The version date should be updated to the
!      current date whenever any portion of the code is updated.  The
!      version number consists of two parts; a major version to the
!      left of the decimal point, and a minor version to the right.
!      The major number should be incremented only when major changes
!      have been made to the overall SARTA code.  The minor number
!      should be incremented only when minor but non-trivial changes
!      are made to the code.  Bug fixes should generally be handled
!      with the version date, but a fix for a serious bug may warrant
!      a change to the minor version number.
!      See the "Doc/last_update.txt" file for a description of the
!      changes associated with every change of VSARTA.
!
       CHARACTER(len=40) :: VSARTA  ! SARTA source code version
       CHARACTER(len=40) :: VSCOEF  ! SARTA coefficient version
       CHARACTER(len=40) :: VTUNNG  ! optical depth tuning version
!      version template    '#.## YYYY-MM-DD <--------comment------->'
       PARAMETER( VSARTA = '1.08 2011-11-21 PGEv6 AB20120210' )
       PARAMETER( VSCOEF = 'AIRS 2008-05-07 m130 m150 CO2=385' )
       PARAMETER( VTUNNG = 'PGEv6; refprof N2O x1/1.04' )

!      *********
!      VARIABLES
!      *********
!      Note: these should not be changed by the user
!
!      ------------------------
!      Constants and other data
!      ------------------------
       REAL     PI ! pi, circle circumference/diameter (3.1415926)
       REAL RADSUN ! radius of the sun (6.956E+8 m)
       REAL     C1 ! radiation constant c1 (1.1911E-8  W/(m2.st.(cm-1)4)
       REAL     C2 ! radiation constant c2 (1.4387863 K/cm-1)
       PARAMETER(    PI = 3.1415926)
       PARAMETER(RADSUN = 6.956E+8)
!
!cc    Previously used values; agrees w/JPL pre-Dec2000
!cc    PARAMETER(    C1 = 1.1910439E-8)  ! JPL value is 1E+3 bigger
!cc    PARAMETER(    C2 = 1.4387687)
!
!      Current values (CODATA98 from NIST); agrees w/JPL Dec2000
       PARAMETER(  C1 = 1.191042722E-8)  ! JPL value is 1E+3 bigger
       PARAMETER(  C2 = 1.4387752)
!
       REAL CO2STD ! standard CO2 PPMV mixing ratio (385)
!      Note: mean global CO2 approximately 385 in 2009
       PARAMETER( CO2STD = 385.0 )
!
       REAL  XSALT ! expected nominal satellite altitude (km)
       PARAMETER( XSALT = 705.0 )
!
!
!      -----------------------------------
!      Channels and layers other variables
!      -----------------------------------
       INTEGER MAXLAY ! # of layers (100)
       INTEGER   NSET ! # of coefficient data sets (7)
       INTEGER MXCHAN ! max total # of channels (2378)
       INTEGER NFCOEF ! # of downwelling thermal "F" factor coefs
       INTEGER MXEMIS ! max # of input emis/rho data points
       INTEGER MAXPRO ! max # of user specified profiles
       INTEGER  MXGAS ! max # of gases in user profile
       INTEGER MXMIEA ! max # of mie particle sizes
       PARAMETER(MAXLAY = 100)
       PARAMETER(  NSET = 7)
       PARAMETER(MXCHAN = 2834)
       PARAMETER(NFCOEF = 6)
       PARAMETER(MXEMIS = 100)
       PARAMETER(MAXPRO = 25)
       PARAMETER( MXGAS = 44)
       PARAMETER(MXMIEA = 10)
!
!***********************************************************************
!      Variables for the coefficient sets
!***********************************************************************
!
!      --------------
!      For set1 = FWO
!      -------------
!      Used in part by modules: 12, 11, 10, 9, 8, 7, 6, 5, 3, 4b, 4a
       INTEGER MXCHN1 ! max # of channels for set1 = FWO (1461)
       INTEGER  N1CON ! # of water con predictors/coefs for set1 (5)
       INTEGER  N1FIX ! # of "fixed" predictors/coefs for set1 (8)
       INTEGER  N1H2O ! # of water predictors/coefs for set1 (13)
       INTEGER   N1O3 ! # of ozone predictors/coefs for set1 (5)
       INTEGER N1COEF ! total # of coefs for set1
       PARAMETER(MXCHN1 = 1461)
       PARAMETER( N1CON = 7)
       PARAMETER( N1FIX = 8)
       PARAMETER( N1H2O = 11)
       PARAMETER(  N1O3 = 5)
       PARAMETER(N1COEF = N1CON + N1FIX + N1H2O + N1O3 )
!
!
!      --------------
!      For set2 = FOW
!      --------------
!      Used in part by modules: 6, 5
       INTEGER MXCHN2 ! max # of channels for set2 = FOW  (325)
       INTEGER  N2CON ! # of water con predictors/coefs for set2 (5)
       INTEGER  N2FIX ! # of "fixed" predictors/coefs for set2 (8)
       INTEGER   N2O3 ! # of ozone predictors/coefs for set2 (10)
       INTEGER  N2H2O ! # of water predictors/coefs for set2 (11)
       INTEGER N2COEF ! total # of coefs for set2
       PARAMETER(MXCHN2 = 325)
       PARAMETER( N2CON = 7)
       PARAMETER( N2FIX = 8)
       PARAMETER(  N2O3 = 10)
       PARAMETER( N2H2O = 11)
       PARAMETER(N2COEF = N2CON + N2FIX + N2O3 + N2H2O )
!
!
!      --------------
!      For set3 = FMW
!      --------------
!      Used in part by modules: 4d, 4c, 3
       INTEGER MXCHN3 ! max # of channels for set3 = FMW  (396)
       INTEGER  N3CON ! # of water con predictors/coefs for set3 (5)
       INTEGER  N3FIX ! # of "fixed" predictors/coefs for set3 (8)
       INTEGER  N3CH4 ! # of methane predictors/coefs for set3 (9)
       INTEGER  N3H2O ! # of water predictors/coefs for set3 (13)
       INTEGER N3COEF ! total # of coefs for set3
       PARAMETER(MXCHN3 = 396)
       PARAMETER( N3CON = 7)
       PARAMETER( N3FIX = 8)
       PARAMETER( N3CH4 = 9)
       PARAMETER( N3H2O = 11)
       PARAMETER(N3COEF = N3CON + N3FIX + N3CH4 + N3H2O )
!
!
!      ---------------
!      For set4 = sun FCOW
!      ---------------
!      Used in part by modules: 2b
       INTEGER MXCHN4 ! max # of channels for set4 = FCOW (85)
       INTEGER  N4CON ! # of water con predictors/coefs for set4 (5)
       INTEGER  N4FIX ! # of "fixed" predictors/coefs for set4 (11)
       INTEGER   N4CO ! # of CO predictors/coefs for set4 (11)
       INTEGER   N4O3 ! # of ozone predictors/coefs for set4 (3)
       INTEGER  N4H2O ! # of water predictors/coefs for set4 (13)
       INTEGER N4COEF ! total # of coefs for set4
       PARAMETER(MXCHN4 = 85)
       PARAMETER( N4CON = 7)
       PARAMETER( N4FIX = 11)
       PARAMETER(  N4CO = 11)
       PARAMETER(  N4O3 = 3)
       PARAMETER( N4H2O = 13)
       PARAMETER(N4COEF = N4CON + N4FIX + N4CO + N4O3 + N4H2O )
!
!
!      -----------------------
!      For set5 = sun BFSW
!      -----------------------
!      Used in part by modules: 2b, 1b
       INTEGER MXCHN5 ! max # of channels for set5 = BFSW (210)
       INTEGER  N5CON ! # of water con predictors/coefs for set5 (5)
       INTEGER  N5FIX ! # of "fixed" predictors/coefs for set5 (11)
       INTEGER  N5H2O ! # of water predictors/coefs for set5 (3)
       INTEGER   N5O3 ! # of ozone predictors/coefs for set5 (1)
       INTEGER N5COEF ! total # of coefs for set5
       PARAMETER(MXCHN5 = 210)
       PARAMETER( N5CON = 7)
       PARAMETER( N5FIX = 11)
       PARAMETER( N5H2O = 3)
       PARAMETER(  N5O3 = 1)
       PARAMETER(N5COEF = N5CON + N5FIX + N5H2O + N5O3 )
!
!
!      -----------------------
!      For set6 = sun MFMW
!      -----------------------
!      Used in part by modules: 1b, 2a
       INTEGER MXCHN6 ! max # of channels for set6 = MFMW (217)
       INTEGER  N6CON ! # of water con predictors/coefs for set6 (5)
       INTEGER  N6FIX ! # of "fixed" predictors/coefs for set6 (8)
       INTEGER  N6H2O ! # of water predictors/coefs for set6 (7)
       INTEGER   N6O3 ! # of ozone predictors/coefs for set6 (1)
       INTEGER N6COEF ! total # of coefs for set6
       PARAMETER(MXCHN6 = 217)
       PARAMETER( N6CON = 7 )
       PARAMETER( N6FIX = 8 )
       PARAMETER( N6H2O = 7 )
       PARAMETER(  N6O3 = 1 )
       PARAMETER(N6COEF = N6CON + N6FIX + N6H2O + N6O3 )
!
!
!      -----------------------
!      For set7 = sun MFBW
!      -----------------------
!      Used in part by modules: 2a, 1a
       INTEGER MXCHN7 ! max # of channels for set7 = MFBW (140)
       INTEGER  N7CON ! # of water con predictors/coefs for set7 (5)
       INTEGER  N7FIX ! # of "fixed" predictors/coefs for set7 (8)
       INTEGER  N7H2O ! # of water predictors/coefs for set7 (13)
       INTEGER   N7O3 ! # of ozone predictors/coefs for set7 (1)
       INTEGER N7COEF ! total # of coefs for set7
       PARAMETER(MXCHN7 = 140)
       PARAMETER( N7CON = 7)
       PARAMETER( N7FIX = 8)
       PARAMETER( N7H2O = 13)
       PARAMETER(  N7O3 = 1)
       PARAMETER(N7COEF = N7CON + N7FIX + N7H2O + N7O3 )
!
!
!      ---------------
!      For trace gases predictors
!      ---------------
       INTEGER NTRACE ! number of trace gas perturbation predictors (7)
       PARAMETER(NTRACE = 7)
!
!
!      ----------------
!      For variable CO2
!      ----------------
!      Used in part by modules: 12, 11, 10, 9, 7, 6, 5, 2b, 1b, 2a
       INTEGER MXCHNC ! max # of channels with CO2 pert coefs (1082)
       INTEGER   NCO2 ! number of CO2 coefficients
       PARAMETER(MXCHNC = 1082)
       PARAMETER(  NCO2 = 5)
!
!
!      ----------------
!      For variable SO2
!      ----------------
       INTEGER MXCHNS ! max # of channels with SO2 pert coefs (602)
       INTEGER   NSO2 ! number of SO2 coefficients
       PARAMETER(MXCHNS = 602)
       PARAMETER(  NSO2 = 4)
!
!
!      -----------------
!      For variable HNO3
!      -----------------
       INTEGER MXCHNH ! max # of channels with HNO3 pert coefs (383)
       INTEGER  NHNO3 ! number of HNO3 coefficients
       PARAMETER(MXCHNH = 383)
       PARAMETER( NHNO3 = 4)
!
!
!      -----------------
!      For variable N2O
!      -----------------
       INTEGER MXCHNN ! max # of channels with N2O pert coefs (586)
       INTEGER   NN2O ! number of N2O coefficients
       PARAMETER(MXCHNN = 586)
       PARAMETER(  NN2O = 7)
!
!
!      ----------------------
!      For OPTRAN water coefs
!      ----------------------
!      Used in part by modules:
       INTEGER MXCHNW ! max # of channelss with OPTRAN H2O coefs (754)
       INTEGER MXOWLY ! number of OPTRAN water layers
       INTEGER NOWAVG ! # of OPTRAN water average profile values (4)
       INTEGER NH2O   ! number of OPTRAN H2O predictors/coefs (9)
       PARAMETER(MXCHNW = 754)
       PARAMETER(MXOWLY = 300)
       PARAMETER(NOWAVG = 4)
       PARAMETER(  NH2O = 9)
!
!      -----------
!      For non-LTE
!      -----------
       INTEGER MXCNTE ! max # of channels for non-LTE (203)
       INTEGER NNCOEF ! # of coefs for non-LTE
       INTEGER NTEBOT ! bottom layer for CO2TOP calc
       REAL CO2NTE ! ref CO2 mixing ratio for non-LTE coefs (ppmv)
       PARAMETER(MXCNTE = 203)
       PARAMETER(NNCOEF = 7)
       PARAMETER(NTEBOT = 10)
       PARAMETER(CO2NTE = 370.0)
!
!      ----------------
!      Single filenames
!      ----------------
       CHARACTER(len=80 ) :: FNFX   ! coef fx
       CHARACTER(len=80 ) :: FNPREF ! reference profile
       PARAMETER(FNFX  = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/fx.txt')
       PARAMETER(FNPREF= &
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/profref_trace385'
!      refprof N2O scaled by 1/1.04 for tuning purposes
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'profref_trace385tuned')
!
!      ----------------
!      Double filenames
!      ----------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       REAL YOFFA          ! yoffset for database "A"
       CHARACTER*80 FACOF1 ! coef set1
       CHARACTER*80 FACOF2 ! coef set2
       CHARACTER*80 FACOF3 ! coef set3
       CHARACTER*80 FACOF4 ! coef set4
       CHARACTER*80 FACOF5 ! coef set5
       CHARACTER*80 FACOF6 ! coef set6
       CHARACTER*80 FACOF7 ! coef set7
       CHARACTER*80 FACO2  ! coef CO2
       CHARACTER*80 FAN2O  ! coef N2O
       CHARACTER*80 FASO2  ! coef SO2
       CHARACTER*80 FAHNO3 ! coef HNO3
       CHARACTER*80 FAOPTR ! coef optran
       CHARACTER*80 FATHER ! coef therm
       CHARACTER*80 FACOFN ! coef non-LTE
       CHARACTER*80 FASUN  ! solar data
       CHARACTER*80 FATMLT ! tuning multiplier
!
       PARAMETER(YOFFA=-13.0)
       PARAMETER(FACOF1= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set1_m130.dat')
       PARAMETER(FACOF2= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set2_m130.dat')
       PARAMETER(FACOF3= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set3_m130.dat')
       PARAMETER(FACOF4= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set4_m130.dat')
       PARAMETER(FACOF5= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set5_m130.dat')
       PARAMETER(FACOF6= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set6_m130.dat')
       PARAMETER(FACOF7= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set7_m130.dat')
       PARAMETER(FACO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'CO2_5term_m130.dat')
       PARAMETER(FAN2O = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'N2O_m130.dat')
       PARAMETER(FASO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'SO2_m130.dat')
       PARAMETER(FAHNO3 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'HNO3_m130.dat')
       PARAMETER(FAOPTR= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'optran_m130.dat')
       PARAMETER(FATHER= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'therm_m130.dat')
       PARAMETER(FACOFN= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'nonLTE7_m130.dat')
       PARAMETER(FASUN = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Solar/solar_m130.txt')
       PARAMETER(FATMLT= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'tunmlt_PGEv6.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_ones.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_wcon_nte.txt')
!cc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
       REAL YOFFB          ! yoffset for database "B"
       CHARACTER*80 FBCOF1 ! coef set1
       CHARACTER*80 FBCOF2 ! coef set2
       CHARACTER*80 FBCOF3 ! coef set3
       CHARACTER*80 FBCOF4 ! coef set4
       CHARACTER*80 FBCOF5 ! coef set5
       CHARACTER*80 FBCOF6 ! coef set6
       CHARACTER*80 FBCOF7 ! coef set7
       CHARACTER*80 FBCO2  ! coef CO2
       CHARACTER*80 FBN2O  ! coef N2O
       CHARACTER*80 FBSO2  ! coef SO2
       CHARACTER*80 FBHNO3 ! coef HNO3
       CHARACTER*80 FBOPTR ! coef optran
       CHARACTER*80 FBTHER ! coef therm
       CHARACTER*80 FBCOFN ! coef non-LTE
       CHARACTER*80 FBSUN  ! solar data
       CHARACTER*80 FBTMLT ! tuning multiplier
!
       PARAMETER(YOFFB=-15.0)
       PARAMETER(FBCOF1= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set1_m150.dat')
       PARAMETER(FBCOF2= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set2_m150.dat')
       PARAMETER(FBCOF3= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set3_m150.dat')
       PARAMETER(FBCOF4= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set4_m150.dat')
       PARAMETER(FBCOF5= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set5_m150.dat')
       PARAMETER(FBCOF6= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set6_m150.dat')
       PARAMETER(FBCOF7= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set7_m150.dat')
       PARAMETER(FBCO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'CO2_5term_m150.dat')
       PARAMETER(FBN2O = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'N2O_m150.dat')
       PARAMETER(FBSO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'SO2_m150.dat')
       PARAMETER(FBHNO3 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'HNO3_m150.dat')
       PARAMETER(FBOPTR= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'optran_m150.dat')
       PARAMETER(FBTHER= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'therm_m150.dat')
       PARAMETER(FBCOFN= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'nonLTE7_m150.dat')
       PARAMETER(FBSUN = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Solar/solar_m150.txt')
       PARAMETER(FBTMLT= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'tunmlt_PGEv6.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_ones.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_wcon_nte.txt')
!cc
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       CHARACTER*80 ABFILE ! A/B weight filename
       PARAMETER(ABFILE= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'ABweight_table_20120210.txt')
!
!      Number of A/B weight epochs
!      Note: must match number in ABFILE
       INTEGER MAXABW
       PARAMETER(MAXABW=6)

!      ----------------------
!      Default & min/max YOFF
!      ----------------------
!      Note: these should be based YOFFA and YOFFB as well as the
!      expected range of variation of YOFF
       REAL YOFDEF   ! default YOFF
       REAL YOFMIN   ! min allowed YOFF
       REAL YOFMAX   ! max allowed YOFF
       PARAMETER(YOFDEF=-14.0)
       PARAMETER(YOFMIN=-16.0)
       PARAMETER(YOFMAX=-11.0)
!
!      ----------------
!      I/O unit numbers
!      ----------------
!      Note: these units are not explicitly openned by the sarta code,
!      they should be set to standard I/O units for your compiler
       INTEGER IOINFO  ! unit number for non-error info messages (6)
       INTEGER IOERR   ! unit number for error messages (2 or 6)
       PARAMETER( IOINFO = 6 )
       PARAMETER( IOERR = 0 )
!
!
!      -----------------
!      Allowed input GUC (Gas Units Code number)
!      -----------------
       INTEGER GUCIN  ! The one & only allowed input GUC number
!      Note: GUCIN must be 1 or 2.  All gases in the input RTP
!      must be of this type.
!       PARAMETER( GUCIN = 1 ) ! GUC number for:  molecules/cm^2
!       PARAMETER( GUCIN = 2 ) ! GUC number for:  kilomoles/cm^2
        PARAMETER( GUCIN = 3 ) ! GUC number for:  ppmv
 

!      End of include file
!=======================================================================

#elif AIRS_PARAM == 2

!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS (Atmospheric Infra-Red Sounder)
!
!    incFTC
!
!F77====================================================================


!ROUTINE NAME:
!    incFTC (include file)


!ABSTRACT:
!    Include file consisting of parameter statements to size various
!    arrays in the USEFAST related routines source code.


!CALL PROTOCOL:
!    none (include file)


!INPUT PARAMETERS:
!    none


!OUTPUT PARAMETERS:
!    none


!INPUT/OUTPUT PARAMETERS:
!    none


!RETURN VALUES:
!    none


!PARENT(S):
!    CALOKW
!    CALOWP
!    CALPAR
!    CALRAD
!    CALT1
!    CALT2
!    CALT3
!    CALT4
!    CALT5
!    CALT6
!    CALT7
!    FAKETZ
!    RDCOEF
!    RDLIST
!    RDPROF
!    RDSUN
!    SUNPAR
!    SARTA


!ROUTINES CALLED:
!    none


!FILES ACCESSED:
!    none


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    Include file for the October 2011 100 layer AIRS
!    Stand Alone RTA (SARTA) code by L.L.Strow/S.Hannon.
!
!    Parameter statements for the FTC routines.


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
! Date        Programmer     Comments
! ----------- -------------- -------------------------------------------
!  1 Dec 1994 Scott Hannon   Created
! 31 Jan 1997 Scott Hannon   Re-wrote for FWO+FOW+FMW+FCOW=Feb97 FTC
!  3 Sep 1997 Scott Hannon   Re-wrote for sets 1 - 7
! 30 Sep 1997 Scott Hannon   Added NCO2 and MXCHNC
! 26 Feb 1998 Scott Hannon   Added OPTRAN variables for water, and
!                            changed both N1H2O & N3H2O from 13 to 11
! 23 Sep 1999 Scott Hannon   Change number of channel dimensions for
!                            new Sep99 version of fast model.
!  5 Apr 2000 Scott Hannon   Added MXEMIS
!  4 Aug 2000 Scott Hannon   Changes values for use with testfast
! 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
! 23 Jan 2001 Scott Hannon   Update values of C1 & C2
! 15 Feb 2001 Scott Hannon   Add MAXPRO, CO2STD, IOERR, IOINFO,
!                            MXGAS, CSARTA, and all filenames
! 24 Apr 2001 Scott Hannon   Add MXMIEA and FNMIEA
! 14 Aug 2001 Scott Hannon   Add FNMIEE and FNMIEG
! 21 Nov 2001 Scott Hannon   Add VSARTA, VSCOEF, & VCLOUD; remove CSARTA
! 12 Sep 2002 Scott Hannon   Updated for m135f (-13.5 um with fringes)
! 17 Dec 2002 Scott Hannon   Updated for revised(Dec02) m135f
!  3 Jan 2003 Scott Hannon   Updated VSARTA for version 1.04
! 06 Feb 2004 Scott Hannon   Add FNTMLT & update VSARTA for v1.05
! 07 Apr 2005 Scott Hannon   NFCOEF increased from 5 to 6 for v1.06
! 18 May 2005 Scott Hannon   update for HNO3 version
! 29 Jun 2005 Scott Hannon   "trace" version v1.07 with CO2,SO2,HNO3,N2O
! 13 Oct 2005 Scott Hannon   Add variables for non-LTE
! 22 Nov 2005 Scott Hannon   Replace set1,set2,CO2 coefs for new M12
! 02 May 2007 Scott Hannon   Added XSALT
! 08 May 2008 Scott Hannon   Updated for v1.08; change most filenames
!                            from FN* to double filenames FA* and FB*;
!                            add YOFFA,YOFFB,YOFMIN,YOFMAX,YOFDEF
! 13 May 2008 Scott Hannon   Add CO2NTE and NTEBOT for non-LTE calc and
!                            increase NNCOEF from 6 to 7
! 12 May 2009 Scott Hannon   Add VTUNNG string; delete VCLOUD
! 11 Oct 2011 Scott Hannon   Add MAXABW; replace FNSCAL with ABFILE

!END====================================================================
!
!-----------------------------------------------------------------------
!      IMPLICIT NONE
!-----------------------------------------------------------------------
! Note: having an "implicit none" in both the include file & the main
! source code will cause some compilers to complain.
!       IMPLICIT NONE


!-----------------------------------------------------------------------
!      INCLUDE FILES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXECUTABLE CODE
!-----------------------------------------------------------------------
!      none

!      -----------------------------------------------------------------
!      Assign SARTA version strings
!      -----------------------------------------------------------------
!      The version strings consists of 3 parts: version number, date,
!      and comment.  The version date should be updated to the
!      current date whenever any portion of the code is updated.  The
!      version number consists of two parts; a major version to the
!      left of the decimal point, and a minor version to the right.
!      The major number should be incremented only when major changes
!      have been made to the overall SARTA code.  The minor number
!      should be incremented only when minor but non-trivial changes
!      are made to the code.  Bug fixes should generally be handled
!      with the version date, but a fix for a serious bug may warrant
!      a change to the minor version number.
!      See the "Doc/last_update.txt" file for a description of the
!      changes associated with every change of VSARTA.
!
       CHARACTER*40 VSARTA  ! SARTA source code version
       CHARACTER*40 VSCOEF  ! SARTA coefficient version
       CHARACTER*40 VTUNNG  ! optical depth tuning version
!      version template    '#.## YYYY-MM-DD <--------comment------->'
       PARAMETER( VSARTA = '1.08 2011-11-21 PGEv6 AB20120210' )
       PARAMETER( VSCOEF = 'AIRS 2008-05-07 m130 m150 CO2=385' )
       PARAMETER( VTUNNG = 'wcon nte' )

!      *********
!      VARIABLES
!      *********
!      Note: these should not be changed by the user
!
!      ------------------------
!      Constants and other data
!      ------------------------
       REAL     PI ! pi, circle circumference/diameter (3.1415926)
       REAL RADSUN ! radius of the sun (6.956E+8 m)
       REAL     C1 ! radiation constant c1 (1.1911E-8  W/(m2.st.(cm-1)4)
       REAL     C2 ! radiation constant c2 (1.4387863 K/cm-1)
       PARAMETER(    PI = 3.1415926)
       PARAMETER(RADSUN = 6.956E+8)
!
!cc    Previously used values; agrees w/JPL pre-Dec2000
!cc    PARAMETER(    C1 = 1.1910439E-8)  ! JPL value is 1E+3 bigger
!cc    PARAMETER(    C2 = 1.4387687)
!
!      Current values (CODATA98 from NIST); agrees w/JPL Dec2000
       PARAMETER(  C1 = 1.191042722E-8)  ! JPL value is 1E+3 bigger
       PARAMETER(  C2 = 1.4387752)
!
       REAL CO2STD ! standard CO2 PPMV mixing ratio (385)
!      Note: mean global CO2 approximately 385 in 2009
       PARAMETER( CO2STD = 385.0 )
!
       REAL  XSALT ! expected nominal satellite altitude (km)
       PARAMETER( XSALT = 705.0 )
!
!
!      -----------------------------------
!      Channels and layers other variables
!      -----------------------------------
       INTEGER MAXLAY ! # of layers (100)
       INTEGER   NSET ! # of coefficient data sets (7)
       INTEGER MXCHAN ! max total # of channels (2378)
       INTEGER NFCOEF ! # of downwelling thermal "F" factor coefs
       INTEGER MXEMIS ! max # of input emis/rho data points
       INTEGER MAXPRO ! max # of user specified profiles
       INTEGER  MXGAS ! max # of gases in user profile
       INTEGER MXMIEA ! max # of mie particle sizes
       PARAMETER(MAXLAY = 100)
       PARAMETER(  NSET = 7)
       PARAMETER(MXCHAN = 2834)
       PARAMETER(NFCOEF = 6)
       PARAMETER(MXEMIS = 100)
       PARAMETER(MAXPRO = 25)
       PARAMETER( MXGAS = 44)
       PARAMETER(MXMIEA = 10)
!
!***********************************************************************
!      Variables for the coefficient sets
!***********************************************************************
!
!      --------------
!      For set1 = FWO
!      -------------
!      Used in part by modules: 12, 11, 10, 9, 8, 7, 6, 5, 3, 4b, 4a
       INTEGER MXCHN1 ! max # of channels for set1 = FWO (1461)
       INTEGER  N1CON ! # of water con predictors/coefs for set1 (5)
       INTEGER  N1FIX ! # of "fixed" predictors/coefs for set1 (8)
       INTEGER  N1H2O ! # of water predictors/coefs for set1 (13)
       INTEGER   N1O3 ! # of ozone predictors/coefs for set1 (5)
       INTEGER N1COEF ! total # of coefs for set1
       PARAMETER(MXCHN1 = 1461)
       PARAMETER( N1CON = 7)
       PARAMETER( N1FIX = 8)
       PARAMETER( N1H2O = 11)
       PARAMETER(  N1O3 = 5)
       PARAMETER(N1COEF = N1CON + N1FIX + N1H2O + N1O3 )
!
!
!      --------------
!      For set2 = FOW
!      --------------
!      Used in part by modules: 6, 5
       INTEGER MXCHN2 ! max # of channels for set2 = FOW  (325)
       INTEGER  N2CON ! # of water con predictors/coefs for set2 (5)
       INTEGER  N2FIX ! # of "fixed" predictors/coefs for set2 (8)
       INTEGER   N2O3 ! # of ozone predictors/coefs for set2 (10)
       INTEGER  N2H2O ! # of water predictors/coefs for set2 (11)
       INTEGER N2COEF ! total # of coefs for set2
       PARAMETER(MXCHN2 = 325)
       PARAMETER( N2CON = 7)
       PARAMETER( N2FIX = 8)
       PARAMETER(  N2O3 = 10)
       PARAMETER( N2H2O = 11)
       PARAMETER(N2COEF = N2CON + N2FIX + N2O3 + N2H2O )
!
!
!      --------------
!      For set3 = FMW
!      --------------
!      Used in part by modules: 4d, 4c, 3
       INTEGER MXCHN3 ! max # of channels for set3 = FMW  (396)
       INTEGER  N3CON ! # of water con predictors/coefs for set3 (5)
       INTEGER  N3FIX ! # of "fixed" predictors/coefs for set3 (8)
       INTEGER  N3CH4 ! # of methane predictors/coefs for set3 (9)
       INTEGER  N3H2O ! # of water predictors/coefs for set3 (13)
       INTEGER N3COEF ! total # of coefs for set3
       PARAMETER(MXCHN3 = 396)
       PARAMETER( N3CON = 7)
       PARAMETER( N3FIX = 8)
       PARAMETER( N3CH4 = 9)
       PARAMETER( N3H2O = 11)
       PARAMETER(N3COEF = N3CON + N3FIX + N3CH4 + N3H2O )
!
!
!      ---------------
!      For set4 = sun FCOW
!      ---------------
!      Used in part by modules: 2b
       INTEGER MXCHN4 ! max # of channels for set4 = FCOW (85)
       INTEGER  N4CON ! # of water con predictors/coefs for set4 (5)
       INTEGER  N4FIX ! # of "fixed" predictors/coefs for set4 (11)
       INTEGER   N4CO ! # of CO predictors/coefs for set4 (11)
       INTEGER   N4O3 ! # of ozone predictors/coefs for set4 (3)
       INTEGER  N4H2O ! # of water predictors/coefs for set4 (13)
       INTEGER N4COEF ! total # of coefs for set4
       PARAMETER(MXCHN4 = 85)
       PARAMETER( N4CON = 7)
       PARAMETER( N4FIX = 11)
       PARAMETER(  N4CO = 11)
       PARAMETER(  N4O3 = 3)
       PARAMETER( N4H2O = 13)
       PARAMETER(N4COEF = N4CON + N4FIX + N4CO + N4O3 + N4H2O )
!
!
!      -----------------------
!      For set5 = sun BFSW
!      -----------------------
!      Used in part by modules: 2b, 1b
       INTEGER MXCHN5 ! max # of channels for set5 = BFSW (210)
       INTEGER  N5CON ! # of water con predictors/coefs for set5 (5)
       INTEGER  N5FIX ! # of "fixed" predictors/coefs for set5 (11)
       INTEGER  N5H2O ! # of water predictors/coefs for set5 (3)
       INTEGER   N5O3 ! # of ozone predictors/coefs for set5 (1)
       INTEGER N5COEF ! total # of coefs for set5
       PARAMETER(MXCHN5 = 210)
       PARAMETER( N5CON = 7)
       PARAMETER( N5FIX = 11)
       PARAMETER( N5H2O = 3)
       PARAMETER(  N5O3 = 1)
       PARAMETER(N5COEF = N5CON + N5FIX + N5H2O + N5O3 )
!
!
!      -----------------------
!      For set6 = sun MFMW
!      -----------------------
!      Used in part by modules: 1b, 2a
       INTEGER MXCHN6 ! max # of channels for set6 = MFMW (217)
       INTEGER  N6CON ! # of water con predictors/coefs for set6 (5)
       INTEGER  N6FIX ! # of "fixed" predictors/coefs for set6 (8)
       INTEGER  N6H2O ! # of water predictors/coefs for set6 (7)
       INTEGER   N6O3 ! # of ozone predictors/coefs for set6 (1)
       INTEGER N6COEF ! total # of coefs for set6
       PARAMETER(MXCHN6 = 217)
       PARAMETER( N6CON = 7 )
       PARAMETER( N6FIX = 8 )
       PARAMETER( N6H2O = 7 )
       PARAMETER(  N6O3 = 1 )
       PARAMETER(N6COEF = N6CON + N6FIX + N6H2O + N6O3 )
!
!
!      -----------------------
!      For set7 = sun MFBW
!      -----------------------
!      Used in part by modules: 2a, 1a
       INTEGER MXCHN7 ! max # of channels for set7 = MFBW (140)
       INTEGER  N7CON ! # of water con predictors/coefs for set7 (5)
       INTEGER  N7FIX ! # of "fixed" predictors/coefs for set7 (8)
       INTEGER  N7H2O ! # of water predictors/coefs for set7 (13)
       INTEGER   N7O3 ! # of ozone predictors/coefs for set7 (1)
       INTEGER N7COEF ! total # of coefs for set7
       PARAMETER(MXCHN7 = 140)
       PARAMETER( N7CON = 7)
       PARAMETER( N7FIX = 8)
       PARAMETER( N7H2O = 13)
       PARAMETER(  N7O3 = 1)
       PARAMETER(N7COEF = N7CON + N7FIX + N7H2O + N7O3 )
!
!
!      ---------------
!      For trace gases predictors
!      ---------------
       INTEGER NTRACE ! number of trace gas perturbation predictors (7)
       PARAMETER(NTRACE = 7)
!
!
!      ----------------
!      For variable CO2
!      ----------------
!      Used in part by modules: 12, 11, 10, 9, 7, 6, 5, 2b, 1b, 2a
       INTEGER MXCHNC ! max # of channels with CO2 pert coefs (1082)
       INTEGER   NCO2 ! number of CO2 coefficients
       PARAMETER(MXCHNC = 1082)
       PARAMETER(  NCO2 = 5)
!
!
!      ----------------
!      For variable SO2
!      ----------------
       INTEGER MXCHNS ! max # of channels with SO2 pert coefs (602)
       INTEGER   NSO2 ! number of SO2 coefficients
       PARAMETER(MXCHNS = 602)
       PARAMETER(  NSO2 = 4)
!
!
!      -----------------
!      For variable HNO3
!      -----------------
       INTEGER MXCHNH ! max # of channels with HNO3 pert coefs (383)
       INTEGER  NHNO3 ! number of HNO3 coefficients
       PARAMETER(MXCHNH = 383)
       PARAMETER( NHNO3 = 4)
!
!
!      -----------------
!      For variable N2O
!      -----------------
       INTEGER MXCHNN ! max # of channels with N2O pert coefs (586)
       INTEGER   NN2O ! number of N2O coefficients
       PARAMETER(MXCHNN = 586)
       PARAMETER(  NN2O = 7)
!
!
!      ----------------------
!      For OPTRAN water coefs
!      ----------------------
!      Used in part by modules:
       INTEGER MXCHNW ! max # of channelss with OPTRAN H2O coefs (754)
       INTEGER MXOWLY ! number of OPTRAN water layers
       INTEGER NOWAVG ! # of OPTRAN water average profile values (4)
       INTEGER NH2O   ! number of OPTRAN H2O predictors/coefs (9)
       PARAMETER(MXCHNW = 754)
       PARAMETER(MXOWLY = 300)
       PARAMETER(NOWAVG = 4)
       PARAMETER(  NH2O = 9)
!
!      -----------
!      For non-LTE
!      -----------
       INTEGER MXCNTE ! max # of channels for non-LTE (203)
       INTEGER NNCOEF ! # of coefs for non-LTE
       INTEGER NTEBOT ! bottom layer for CO2TOP calc
       REAL CO2NTE ! ref CO2 mixing ratio for non-LTE coefs (ppmv)
       PARAMETER(MXCNTE = 203)
       PARAMETER(NNCOEF = 7)
       PARAMETER(NTEBOT = 10)
       PARAMETER(CO2NTE = 370.0)
!
!      ----------------
!      Single filenames
!      ----------------
       CHARACTER*80 FNFX   ! coef fx
       CHARACTER*80 FNPREF ! reference profile
       PARAMETER(FNFX  = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/fx.txt')
       PARAMETER(FNPREF= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/profref_trace385')
!
!      ----------------
!      Double filenames
!      ----------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       REAL YOFFA          ! yoffset for database "A"
       CHARACTER*80 FACOF1 ! coef set1
       CHARACTER*80 FACOF2 ! coef set2
       CHARACTER*80 FACOF3 ! coef set3
       CHARACTER*80 FACOF4 ! coef set4
       CHARACTER*80 FACOF5 ! coef set5
       CHARACTER*80 FACOF6 ! coef set6
       CHARACTER*80 FACOF7 ! coef set7
       CHARACTER*80 FACO2  ! coef CO2
       CHARACTER*80 FAN2O  ! coef N2O
       CHARACTER*80 FASO2  ! coef SO2
       CHARACTER*80 FAHNO3 ! coef HNO3
       CHARACTER*80 FAOPTR ! coef optran
       CHARACTER*80 FATHER ! coef therm
       CHARACTER*80 FACOFN ! coef non-LTE
       CHARACTER*80 FASUN  ! solar data
       CHARACTER*80 FATMLT ! tuning multiplier
!
       PARAMETER(YOFFA=-13.0)
       PARAMETER(FACOF1= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set1_m130.dat')
       PARAMETER(FACOF2= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set2_m130.dat')
       PARAMETER(FACOF3= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set3_m130.dat')
       PARAMETER(FACOF4= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set4_m130.dat')
       PARAMETER(FACOF5= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set5_m130.dat')
       PARAMETER(FACOF6= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set6_m130.dat')
       PARAMETER(FACOF7= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set7_m130.dat')
       PARAMETER(FACO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'CO2_5term_m130.dat')
       PARAMETER(FAN2O = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'N2O_m130.dat')
       PARAMETER(FASO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'SO2_m130.dat')
       PARAMETER(FAHNO3 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'HNO3_m130.dat')
       PARAMETER(FAOPTR= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'optran_m130.dat')
       PARAMETER(FATHER= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'therm_m130.dat')
       PARAMETER(FACOFN= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'nonLTE7_m130.dat')
       PARAMETER(FASUN = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Solar/solar_m130.txt')
       PARAMETER(FATMLT= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'tunmlt_wcon_nte.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_ones.txt')
!cc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
       REAL YOFFB          ! yoffset for database "B"
       CHARACTER*80 FBCOF1 ! coef set1
       CHARACTER*80 FBCOF2 ! coef set2
       CHARACTER*80 FBCOF3 ! coef set3
       CHARACTER*80 FBCOF4 ! coef set4
       CHARACTER*80 FBCOF5 ! coef set5
       CHARACTER*80 FBCOF6 ! coef set6
       CHARACTER*80 FBCOF7 ! coef set7
       CHARACTER*80 FBCO2  ! coef CO2
       CHARACTER*80 FBN2O  ! coef N2O
       CHARACTER*80 FBSO2  ! coef SO2
       CHARACTER*80 FBHNO3 ! coef HNO3
       CHARACTER*80 FBOPTR ! coef optran
       CHARACTER*80 FBTHER ! coef therm
       CHARACTER*80 FBCOFN ! coef non-LTE
       CHARACTER*80 FBSUN  ! solar data
       CHARACTER*80 FBTMLT ! tuning multiplier
!
       PARAMETER(YOFFB=-15.0)
       PARAMETER(FBCOF1= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set1_m150.dat')
       PARAMETER(FBCOF2= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set2_m150.dat')
       PARAMETER(FBCOF3= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set3_m150.dat')
       PARAMETER(FBCOF4= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set4_m150.dat')
       PARAMETER(FBCOF5= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set5_m150.dat')
       PARAMETER(FBCOF6= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set6_m150.dat')
       PARAMETER(FBCOF7= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set7_m150.dat')
       PARAMETER(FBCO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'CO2_5term_m150.dat')
       PARAMETER(FBN2O = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'N2O_m150.dat')
       PARAMETER(FBSO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'SO2_m150.dat')
       PARAMETER(FBHNO3 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'HNO3_m150.dat')
       PARAMETER(FBOPTR= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'optran_m150.dat')
       PARAMETER(FBTHER= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'therm_m150.dat')
       PARAMETER(FBCOFN= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'nonLTE7_m150.dat')
       PARAMETER(FBSUN = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Solar/solar_m150.txt')
       PARAMETER(FBTMLT= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'tunmlt_wcon_nte.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_ones.txt')
!cc
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       CHARACTER*80 ABFILE ! A/B weight filename
       PARAMETER(ABFILE= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'ABweight_table_20120210.txt')
!
!      Number of A/B weight epochs
!      Note: must match number in ABFILE
       INTEGER MAXABW
       PARAMETER(MAXABW=6)

!      ----------------------
!      Default & min/max YOFF
!      ----------------------
!      Note: these should be based YOFFA and YOFFB as well as the
!      expected range of variation of YOFF
       REAL YOFDEF   ! default YOFF
       REAL YOFMIN   ! min allowed YOFF
       REAL YOFMAX   ! max allowed YOFF
       PARAMETER(YOFDEF=-14.0)
       PARAMETER(YOFMIN=-16.0)
       PARAMETER(YOFMAX=-11.0)
!
!      ----------------
!      I/O unit numbers
!      ----------------
!      Note: these units are not explicitly openned by the sarta code,
!      they should be set to standard I/O units for your compiler
       INTEGER IOINFO  ! unit number for non-error info messages (6)
       INTEGER IOERR   ! unit number for error messages (2 or 6)
       PARAMETER( IOINFO = 6 )
       PARAMETER( IOERR = 0 )
!
!
!      -----------------
!      Allowed input GUC (Gas Units Code number)
!      -----------------
       INTEGER GUCIN  ! The one & only allowed input GUC number
!      Note: GUCIN must be 1 or 2.  All gases in the input RTP
!      must be of this type.
!       PARAMETER( GUCIN = 1 ) ! GUC number for:  molecules/cm^2
!       PARAMETER( GUCIN = 2 ) ! GUC number for:  kilomoles/cm^2
        PARAMETER( GUCIN = 3 ) ! GUC number for:  ppmv

!      End of include file
!=======================================================================

#elif AIRS_PARAM == 3

!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS (Atmospheric Infra-Red Sounder)
!
!    incFTC
!
!F77====================================================================


!ROUTINE NAME:
!    incFTC (include file)


!ABSTRACT:
!    Include file consisting of parameter statements to size various
!    arrays in the USEFAST related routines source code.


!CALL PROTOCOL:
!    none (include file)


!INPUT PARAMETERS:
!    none


!OUTPUT PARAMETERS:
!    none


!INPUT/OUTPUT PARAMETERS:
!    none


!RETURN VALUES:
!    none


!PARENT(S):
!    CALOKW
!    CALOWP
!    CALPAR
!    CALRAD
!    CALT1
!    CALT2
!    CALT3
!    CALT4
!    CALT5
!    CALT6
!    CALT7
!    FAKETZ
!    RDCOEF
!    RDLIST
!    RDPROF
!    RDSUN
!    SUNPAR
!    SARTA


!ROUTINES CALLED:
!    none


!FILES ACCESSED:
!    none


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    Include file for the May 2008 100 layer AIRS fast
!    Stand Alone RTA (SARTA) code by L.L.Strow/S.Hannon.
!
!    Parameter statements for the FTC routines.


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
! Date        Programmer     Comments
! ----------- -------------- -------------------------------------------
!  1 Dec 1994 Scott Hannon   Created
! 31 Jan 1997 Scott Hannon   Re-wrote for FWO+FOW+FMW+FCOW=Feb97 FTC
!  3 Sep 1997 Scott Hannon   Re-wrote for sets 1 - 7
! 30 Sep 1997 Scott Hannon   Added NCO2 and MXCHNC
! 26 Feb 1998 Scott Hannon   Added OPTRAN variables for water, and
!                            changed both N1H2O & N3H2O from 13 to 11
! 23 Sep 1999 Scott Hannon   Change number of channel dimensions for
!                            new Sep99 version of fast model.
!  5 Apr 2000 Scott Hannon   Added MXEMIS
!  4 Aug 2000 Scott Hannon   Changes values for use with testfast
! 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
! 23 Jan 2001 Scott Hannon   Update values of C1 & C2
! 15 Feb 2001 Scott Hannon   Add MAXPRO, CO2STD, IOERR, IOINFO,
!                            MXGAS, CSARTA, and all filenames
! 24 Apr 2001 Scott Hannon   Add MXMIEA and FNMIEA
! 14 Aug 2001 Scott Hannon   Add FNMIEE and FNMIEG
! 21 Nov 2001 Scott Hannon   Add VSARTA, VSCOEF, & VCLOUD; remove CSARTA
! 12 Sep 2002 Scott Hannon   Updated for m135f (-13.5 um with fringes)
! 17 Dec 2002 Scott Hannon   Updated for revised(Dec02) m135f
!  3 Jan 2003 Scott Hannon   Updated VSARTA for version 1.04
! 06 Feb 2004 Scott Hannon   Add FNTMLT & update VSARTA for v1.05
! 07 Apr 2005 Scott Hannon   NFCOEF increased from 5 to 6 for v1.06
! 18 May 2005 Scott Hannon   update for HNO3 version
! 29 Jun 2005 Scott Hannon   "trace" version v1.07 with CO2,SO2,HNO3,N2O
! 13 Oct 2005 Scott Hannon   Add variables for non-LTE
! 22 Nov 2005 Scott Hannon   Replace set1,set2,CO2 coefs for new M12
! 02 May 2007 Scott Hannon   Added XSALT
! 08 May 2008 Scott Hannon   Updated for v1.08; change most filenames
!                            from FN* to double filenames FA* and FB*;
!                            add YOFFA,YOFFB,YOFMIN,YOFMAX,YOFDEF
! 13 May 2008 Scott Hannon   Add CO2NTE and NTEBOT for non-LTE calc and
!                            increase NNCOEF from 6 to 7
! 12 May 2009 Scott Hannon   Add VTUNNG string; delete VCLOUD
! 04 Aug 2009 Scott Hannon   Add FNSCAL for speccal text file
! 11 Oct 2011 Scott Hannon   Add MAXABW; replace FNSCAL with ABFILE

!END====================================================================
!
!-----------------------------------------------------------------------
!      IMPLICIT NONE
!-----------------------------------------------------------------------
! Note: having an "implicit none" in both the include file & the main
! source code will cause some compilers to complain.
!       IMPLICIT NONE


!-----------------------------------------------------------------------
!      INCLUDE FILES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXECUTABLE CODE
!-----------------------------------------------------------------------
!      none

!      -----------------------------------------------------------------
!      Assign SARTA version strings
!      -----------------------------------------------------------------
!      The version strings consists of 3 parts: version number, date,
!      and comment.  The version date should be updated to the
!      current date whenever any portion of the code is updated.  The
!      version number consists of two parts; a major version to the
!      left of the decimal point, and a minor version to the right.
!      The major number should be incremented only when major changes
!      have been made to the overall SARTA code.  The minor number
!      should be incremented only when minor but non-trivial changes
!      are made to the code.  Bug fixes should generally be handled
!      with the version date, but a fix for a serious bug may warrant
!      a change to the minor version number.
!      See the "Doc/last_update.txt" file for a description of the
!      changes associated with every change of VSARTA.
!
       CHARACTER*40 VSARTA  ! SARTA source code version
       CHARACTER*40 VSCOEF  ! SARTA coefficient version
       CHARACTER*40 VTUNNG  ! optical depth tuning version
!      version template    '#.## YYYY-MM-DD <--------comment------->'
       PARAMETER( VSARTA = '1.08 2011-11-21 PGEv6 AB20120210' )
       PARAMETER( VSCOEF = 'AIRS 2008-05-07 m130x m140x CO2=370' )
       PARAMETER( VTUNNG = 'PGEv6; refprof N2O x1/1.04' )

!      *********
!      VARIABLES
!      *********
!      Note: these should not be changed by the user
!
!      ------------------------
!      Constants and other data
!      ------------------------
       REAL     PI ! pi, circle circumference/diameter (3.1415926)
       REAL RADSUN ! radius of the sun (6.956E+8 m)
       REAL     C1 ! radiation constant c1 (1.1911E-8  W/(m2.st.(cm-1)4)
       REAL     C2 ! radiation constant c2 (1.4387863 K/cm-1)
       PARAMETER(    PI = 3.1415926)
       PARAMETER(RADSUN = 6.956E+8)
!
!cc    Previously used values; agrees w/JPL pre-Dec2000
!cc    PARAMETER(    C1 = 1.1910439E-8)  ! JPL value is 1E+3 bigger
!cc    PARAMETER(    C2 = 1.4387687)
!
!      Current values (CODATA98 from NIST); agrees w/JPL Dec2000
       PARAMETER(  C1 = 1.191042722E-8)  ! JPL value is 1E+3 bigger
       PARAMETER(  C2 = 1.4387752)
!
       REAL CO2STD ! standard CO2 PPMV mixing ratio (370)
!      Note: mean global CO2 approximately 370 in 2002
       PARAMETER( CO2STD = 370.0 )
!
       REAL  XSALT ! expected nominal satellite altitude (km)
       PARAMETER( XSALT = 705.0 )
!
!
!      -----------------------------------
!      Channels and layers other variables
!      -----------------------------------
       INTEGER MAXLAY ! # of layers (100)
       INTEGER   NSET ! # of coefficient data sets (7)
       INTEGER MXCHAN ! max total # of channels (2378)
       INTEGER NFCOEF ! # of downwelling thermal "F" factor coefs
       INTEGER MXEMIS ! max # of input emis/rho data points
       INTEGER MAXPRO ! max # of user specified profiles
       INTEGER  MXGAS ! max # of gases in user profile
       INTEGER MXMIEA ! max # of mie particle sizes
       PARAMETER(MAXLAY = 100)
       PARAMETER(  NSET = 7)
       PARAMETER(MXCHAN = 2834)
       PARAMETER(NFCOEF = 6)
       PARAMETER(MXEMIS = 100)
       PARAMETER(MAXPRO = 25)
       PARAMETER( MXGAS = 44)
       PARAMETER(MXMIEA = 10)
!
!***********************************************************************
!      Variables for the coefficient sets
!***********************************************************************
!
!      --------------
!      For set1 = FWO
!      -------------
!      Used in part by modules: 12, 11, 10, 9, 8, 7, 6, 5, 3, 4b, 4a
       INTEGER MXCHN1 ! max # of channels for set1 = FWO (1461)
       INTEGER  N1CON ! # of water con predictors/coefs for set1 (5)
       INTEGER  N1FIX ! # of "fixed" predictors/coefs for set1 (8)
       INTEGER  N1H2O ! # of water predictors/coefs for set1 (13)
       INTEGER   N1O3 ! # of ozone predictors/coefs for set1 (5)
       INTEGER N1COEF ! total # of coefs for set1
       PARAMETER(MXCHN1 = 1461)
       PARAMETER( N1CON = 7)
       PARAMETER( N1FIX = 8)
       PARAMETER( N1H2O = 11)
       PARAMETER(  N1O3 = 5)
       PARAMETER(N1COEF = N1CON + N1FIX + N1H2O + N1O3 )
!
!
!      --------------
!      For set2 = FOW
!      --------------
!      Used in part by modules: 6, 5
       INTEGER MXCHN2 ! max # of channels for set2 = FOW  (325)
       INTEGER  N2CON ! # of water con predictors/coefs for set2 (5)
       INTEGER  N2FIX ! # of "fixed" predictors/coefs for set2 (8)
       INTEGER   N2O3 ! # of ozone predictors/coefs for set2 (10)
       INTEGER  N2H2O ! # of water predictors/coefs for set2 (11)
       INTEGER N2COEF ! total # of coefs for set2
       PARAMETER(MXCHN2 = 325)
       PARAMETER( N2CON = 7)
       PARAMETER( N2FIX = 8)
       PARAMETER(  N2O3 = 10)
       PARAMETER( N2H2O = 11)
       PARAMETER(N2COEF = N2CON + N2FIX + N2O3 + N2H2O )
!
!
!      --------------
!      For set3 = FMW
!      --------------
!      Used in part by modules: 4d, 4c, 3
       INTEGER MXCHN3 ! max # of channels for set3 = FMW  (396)
       INTEGER  N3CON ! # of water con predictors/coefs for set3 (5)
       INTEGER  N3FIX ! # of "fixed" predictors/coefs for set3 (8)
       INTEGER  N3CH4 ! # of methane predictors/coefs for set3 (9)
       INTEGER  N3H2O ! # of water predictors/coefs for set3 (13)
       INTEGER N3COEF ! total # of coefs for set3
       PARAMETER(MXCHN3 = 396)
       PARAMETER( N3CON = 7)
       PARAMETER( N3FIX = 8)
       PARAMETER( N3CH4 = 9)
       PARAMETER( N3H2O = 11)
       PARAMETER(N3COEF = N3CON + N3FIX + N3CH4 + N3H2O )
!
!
!      ---------------
!      For set4 = sun FCOW
!      ---------------
!      Used in part by modules: 2b
       INTEGER MXCHN4 ! max # of channels for set4 = FCOW (85)
       INTEGER  N4CON ! # of water con predictors/coefs for set4 (5)
       INTEGER  N4FIX ! # of "fixed" predictors/coefs for set4 (11)
       INTEGER   N4CO ! # of CO predictors/coefs for set4 (11)
       INTEGER   N4O3 ! # of ozone predictors/coefs for set4 (3)
       INTEGER  N4H2O ! # of water predictors/coefs for set4 (13)
       INTEGER N4COEF ! total # of coefs for set4
       PARAMETER(MXCHN4 = 85)
       PARAMETER( N4CON = 7)
       PARAMETER( N4FIX = 11)
       PARAMETER(  N4CO = 11)
       PARAMETER(  N4O3 = 3)
       PARAMETER( N4H2O = 13)
       PARAMETER(N4COEF = N4CON + N4FIX + N4CO + N4O3 + N4H2O )
!
!
!      -----------------------
!      For set5 = sun BFSW
!      -----------------------
!      Used in part by modules: 2b, 1b
       INTEGER MXCHN5 ! max # of channels for set5 = BFSW (210)
       INTEGER  N5CON ! # of water con predictors/coefs for set5 (5)
       INTEGER  N5FIX ! # of "fixed" predictors/coefs for set5 (11)
       INTEGER  N5H2O ! # of water predictors/coefs for set5 (3)
       INTEGER   N5O3 ! # of ozone predictors/coefs for set5 (1)
       INTEGER N5COEF ! total # of coefs for set5
       PARAMETER(MXCHN5 = 210)
       PARAMETER( N5CON = 7)
       PARAMETER( N5FIX = 11)
       PARAMETER( N5H2O = 3)
       PARAMETER(  N5O3 = 1)
       PARAMETER(N5COEF = N5CON + N5FIX + N5H2O + N5O3 )
!
!
!      -----------------------
!      For set6 = sun MFMW
!      -----------------------
!      Used in part by modules: 1b, 2a
       INTEGER MXCHN6 ! max # of channels for set6 = MFMW (217)
       INTEGER  N6CON ! # of water con predictors/coefs for set6 (5)
       INTEGER  N6FIX ! # of "fixed" predictors/coefs for set6 (8)
       INTEGER  N6H2O ! # of water predictors/coefs for set6 (7)
       INTEGER   N6O3 ! # of ozone predictors/coefs for set6 (1)
       INTEGER N6COEF ! total # of coefs for set6
       PARAMETER(MXCHN6 = 217)
       PARAMETER( N6CON = 7 )
       PARAMETER( N6FIX = 8 )
       PARAMETER( N6H2O = 7 )
       PARAMETER(  N6O3 = 1 )
       PARAMETER(N6COEF = N6CON + N6FIX + N6H2O + N6O3 )
!
!
!      -----------------------
!      For set7 = sun MFBW
!      -----------------------
!      Used in part by modules: 2a, 1a
       INTEGER MXCHN7 ! max # of channels for set7 = MFBW (140)
       INTEGER  N7CON ! # of water con predictors/coefs for set7 (5)
       INTEGER  N7FIX ! # of "fixed" predictors/coefs for set7 (8)
       INTEGER  N7H2O ! # of water predictors/coefs for set7 (13)
       INTEGER   N7O3 ! # of ozone predictors/coefs for set7 (1)
       INTEGER N7COEF ! total # of coefs for set7
       PARAMETER(MXCHN7 = 140)
       PARAMETER( N7CON = 7)
       PARAMETER( N7FIX = 8)
       PARAMETER( N7H2O = 13)
       PARAMETER(  N7O3 = 1)
       PARAMETER(N7COEF = N7CON + N7FIX + N7H2O + N7O3 )
!
!
!      ---------------
!      For trace gases predictors
!      ---------------
       INTEGER NTRACE ! number of trace gas perturbation predictors (7)
       PARAMETER(NTRACE = 7)
!
!
!      ----------------
!      For variable CO2
!      ----------------
!      Used in part by modules: 12, 11, 10, 9, 7, 6, 5, 2b, 1b, 2a
       INTEGER MXCHNC ! max # of channels with CO2 pert coefs (1082)
       INTEGER   NCO2 ! number of CO2 coefficients
       PARAMETER(MXCHNC = 1082)
       PARAMETER(  NCO2 = 5)
!
!
!      ----------------
!      For variable SO2
!      ----------------
       INTEGER MXCHNS ! max # of channels with SO2 pert coefs (602)
       INTEGER   NSO2 ! number of SO2 coefficients
       PARAMETER(MXCHNS = 602)
       PARAMETER(  NSO2 = 4)
!
!
!      -----------------
!      For variable HNO3
!      -----------------
       INTEGER MXCHNH ! max # of channels with HNO3 pert coefs (383)
       INTEGER  NHNO3 ! number of HNO3 coefficients
       PARAMETER(MXCHNH = 383)
       PARAMETER( NHNO3 = 4)
!
!
!      -----------------
!      For variable N2O
!      -----------------
       INTEGER MXCHNN ! max # of channels with N2O pert coefs (586)
       INTEGER   NN2O ! number of N2O coefficients
       PARAMETER(MXCHNN = 586)
       PARAMETER(  NN2O = 7)
!
!
!      ----------------------
!      For OPTRAN water coefs
!      ----------------------
!      Used in part by modules:
       INTEGER MXCHNW ! max # of channelss with OPTRAN H2O coefs (754)
       INTEGER MXOWLY ! number of OPTRAN water layers
       INTEGER NOWAVG ! # of OPTRAN water average profile values (4)
       INTEGER NH2O   ! number of OPTRAN H2O predictors/coefs (9)
       PARAMETER(MXCHNW = 754)
       PARAMETER(MXOWLY = 300)
       PARAMETER(NOWAVG = 4)
       PARAMETER(  NH2O = 9)
!
!      -----------
!      For non-LTE
!      -----------
       INTEGER MXCNTE ! max # of channels for non-LTE (203)
       INTEGER NNCOEF ! # of coefs for non-LTE
       INTEGER NTEBOT ! bottom layer for CO2TOP calc
       REAL CO2NTE ! ref CO2 mixing ratio for non-LTE coefs (ppmv)
       PARAMETER(MXCNTE = 203)
       PARAMETER(NNCOEF = 7)
       PARAMETER(NTEBOT = 10)
       PARAMETER(CO2NTE = 370.0)
!
!      ----------------
!      Single filenames
!      ----------------
       CHARACTER*80 FNFX   ! coef fx
       CHARACTER*80 FNPREF ! reference profile
       PARAMETER(FNFX  = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/fx.txt')
       PARAMETER(FNPREF= &
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/profref_trace370')
!      refprof N2O scaled by 1/1.04 for tuning purposes
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'profref_trace370tuned')
!
!      ----------------
!      Double filenames
!      ----------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       REAL YOFFA          ! yoffset for database "A"
       CHARACTER*80 FACOF1 ! coef set1
       CHARACTER*80 FACOF2 ! coef set2
       CHARACTER*80 FACOF3 ! coef set3
       CHARACTER*80 FACOF4 ! coef set4
       CHARACTER*80 FACOF5 ! coef set5
       CHARACTER*80 FACOF6 ! coef set6
       CHARACTER*80 FACOF7 ! coef set7
       CHARACTER*80 FACO2  ! coef CO2
       CHARACTER*80 FAN2O  ! coef N2O
       CHARACTER*80 FASO2  ! coef SO2
       CHARACTER*80 FAHNO3 ! coef HNO3
       CHARACTER*80 FAOPTR ! coef optran
       CHARACTER*80 FATHER ! coef therm
       CHARACTER*80 FACOFN ! coef non-LTE
       CHARACTER*80 FASUN  ! solar data
       CHARACTER*80 FATMLT ! tuning multiplier
!
       PARAMETER(YOFFA=-13.0)
       PARAMETER(FACOF1= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set1_m130x370.dat')
       PARAMETER(FACOF2= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set2_m130x370.dat')
       PARAMETER(FACOF3= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set3_m130x370.dat')
       PARAMETER(FACOF4= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set4_m130x370.dat')
       PARAMETER(FACOF5= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set5_m130x370.dat')
       PARAMETER(FACOF6= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set6_m130x370.dat')
       PARAMETER(FACOF7= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set7_m130x370.dat')
       PARAMETER(FACO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'CO2_5term_m130x370.dat')
       PARAMETER(FAN2O = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'N2O_m130x370.dat')
       PARAMETER(FASO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'SO2_m130x370.dat')
       PARAMETER(FAHNO3 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'HNO3_m130x370.dat')
       PARAMETER(FAOPTR= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'optran_m130x370.dat')
       PARAMETER(FATHER= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'therm_m130x370.dat')
       PARAMETER(FACOFN= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'nonLTE7_m130x.dat')
       PARAMETER(FASUN = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Solar/solar_m130x.txt')
!
       PARAMETER(FATMLT= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'tunmlt_PGEv6.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_ones.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_wcon_nte.txt')
!cc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
       REAL YOFFB          ! yoffset for database "B"
       CHARACTER*80 FBCOF1 ! coef set1
       CHARACTER*80 FBCOF2 ! coef set2
       CHARACTER*80 FBCOF3 ! coef set3
       CHARACTER*80 FBCOF4 ! coef set4
       CHARACTER*80 FBCOF5 ! coef set5
       CHARACTER*80 FBCOF6 ! coef set6
       CHARACTER*80 FBCOF7 ! coef set7
       CHARACTER*80 FBCO2  ! coef CO2
       CHARACTER*80 FBN2O  ! coef N2O
       CHARACTER*80 FBSO2  ! coef SO2
       CHARACTER*80 FBHNO3 ! coef HNO3
       CHARACTER*80 FBOPTR ! coef optran
       CHARACTER*80 FBTHER ! coef therm
       CHARACTER*80 FBCOFN ! coef non-LTE
       CHARACTER*80 FBSUN  ! solar data
       CHARACTER*80 FBTMLT ! tuning multiplier
!
       PARAMETER(YOFFB=-14.0)
       PARAMETER(FBCOF1= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set1_m140x370.dat')
       PARAMETER(FBCOF2= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set2_m140x370.dat')
       PARAMETER(FBCOF3= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set3_m140x370.dat')
       PARAMETER(FBCOF4= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set4_m140x370.dat')
       PARAMETER(FBCOF5= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set5_m140x370.dat')
       PARAMETER(FBCOF6= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set6_m140x370.dat')
       PARAMETER(FBCOF7= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set7_m140x370.dat')
       PARAMETER(FBCO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'CO2_5term_m140x370.dat')
       PARAMETER(FBN2O = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'N2O_m140x370.dat')
       PARAMETER(FBSO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'SO2_m140x370.dat')
       PARAMETER(FBHNO3 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'HNO3_m140x370.dat')
       PARAMETER(FBOPTR= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'optran_m140x370.dat')
       PARAMETER(FBTHER= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'therm_m140x370.dat')
       PARAMETER(FBCOFN= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'nonLTE7_m140x.dat')
       PARAMETER(FBSUN = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Solar/solar_m140x.txt')
!
       PARAMETER(FBTMLT= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'tunmlt_PGEv6.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_ones.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_wcon_nte.txt')
!cc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       CHARACTER*80 ABFILE ! spectral calibration
       PARAMETER(ABFILE= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'ABweight_table_20120210.txt')
!
!      Number of A/B weight epochs
!      Note: must match number in ABFILE
       INTEGER MAXABW
       PARAMETER(MAXABW=6)
!
!      ----------------------
!      Default & min/max YOFF
!      ----------------------
!      Note: these should be based YOFFA and YOFFB as well as the
!      expected range of variation of YOFF
       REAL YOFDEF   ! default YOFF
       REAL YOFMIN   ! min allowed YOFF
       REAL YOFMAX   ! max allowed YOFF
       PARAMETER(YOFDEF=-13.5)
       PARAMETER(YOFMIN=-16.0)
       PARAMETER(YOFMAX=-11.0)

!      ----------------
!      I/O unit numbers
!      ----------------
!      Note: these units are not explicitly openned by the sarta code,
!      they should be set to standard I/O units for your compiler
       INTEGER IOINFO  ! unit number for non-error info messages (6)
       INTEGER IOERR   ! unit number for error messages (2 or 6)
       PARAMETER( IOINFO = 6 )
       PARAMETER( IOERR = 0 )
!
!
!      -----------------
!      Allowed input GUC (Gas Units Code number)
!      -----------------
       INTEGER GUCIN  ! The one & only allowed input GUC number
!      Note: GUCIN must be 1 or 2.  All gases in the input RTP
!      must be of this type.
!       PARAMETER( GUCIN = 1 ) ! GUC number for:  molecules/cm^2
!       PARAMETER( GUCIN = 2 ) ! GUC number for:  kilomoles/cm^2
        PARAMETER( GUCIN = 3 ) ! GUC number for:  ppmv

!      End of include file
!=======================================================================

#elif AIRS_PARAM == 4

!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS (Atmospheric Infra-Red Sounder)
!
!    incFTC
!
!F77====================================================================


!ROUTINE NAME:
!    incFTC (include file)


!ABSTRACT:
!    Include file consisting of parameter statements to size various
!    arrays in the USEFAST related routines source code.


!CALL PROTOCOL:
!    none (include file)


!INPUT PARAMETERS:
!    none


!OUTPUT PARAMETERS:
!    none


!INPUT/OUTPUT PARAMETERS:
!    none


!RETURN VALUES:
!    none


!PARENT(S):
!    CALOKW
!    CALOWP
!    CALPAR
!    CALRAD
!    CALT1
!    CALT2
!    CALT3
!    CALT4
!    CALT5
!    CALT6
!    CALT7
!    FAKETZ
!    RDCOEF
!    RDLIST
!    RDPROF
!    RDSUN
!    SUNPAR
!    SARTA


!ROUTINES CALLED:
!    none


!FILES ACCESSED:
!    none


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    Include file for the May 2008 100 layer AIRS fast
!    Stand Alone RTA (SARTA) code by L.L.Strow/S.Hannon.
!
!    Parameter statements for the FTC routines.


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
! Date        Programmer     Comments
! ----------- -------------- -------------------------------------------
!  1 Dec 1994 Scott Hannon   Created
! 31 Jan 1997 Scott Hannon   Re-wrote for FWO+FOW+FMW+FCOW=Feb97 FTC
!  3 Sep 1997 Scott Hannon   Re-wrote for sets 1 - 7
! 30 Sep 1997 Scott Hannon   Added NCO2 and MXCHNC
! 26 Feb 1998 Scott Hannon   Added OPTRAN variables for water, and
!                            changed both N1H2O & N3H2O from 13 to 11
! 23 Sep 1999 Scott Hannon   Change number of channel dimensions for
!                            new Sep99 version of fast model.
!  5 Apr 2000 Scott Hannon   Added MXEMIS
!  4 Aug 2000 Scott Hannon   Changes values for use with testfast
! 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
! 23 Jan 2001 Scott Hannon   Update values of C1 & C2
! 15 Feb 2001 Scott Hannon   Add MAXPRO, CO2STD, IOERR, IOINFO,
!                            MXGAS, CSARTA, and all filenames
! 24 Apr 2001 Scott Hannon   Add MXMIEA and FNMIEA
! 14 Aug 2001 Scott Hannon   Add FNMIEE and FNMIEG
! 21 Nov 2001 Scott Hannon   Add VSARTA, VSCOEF, & VCLOUD; remove CSARTA
! 12 Sep 2002 Scott Hannon   Updated for m135f (-13.5 um with fringes)
! 17 Dec 2002 Scott Hannon   Updated for revised(Dec02) m135f
!  3 Jan 2003 Scott Hannon   Updated VSARTA for version 1.04
! 06 Feb 2004 Scott Hannon   Add FNTMLT & update VSARTA for v1.05
! 07 Apr 2005 Scott Hannon   NFCOEF increased from 5 to 6 for v1.06
! 18 May 2005 Scott Hannon   update for HNO3 version
! 29 Jun 2005 Scott Hannon   "trace" version v1.07 with CO2,SO2,HNO3,N2O
! 13 Oct 2005 Scott Hannon   Add variables for non-LTE
! 22 Nov 2005 Scott Hannon   Replace set1,set2,CO2 coefs for new M12
! 02 May 2007 Scott Hannon   Added XSALT
! 08 May 2008 Scott Hannon   Updated for v1.08; change most filenames
!                            from FN* to double filenames FA* and FB*;
!                            add YOFFA,YOFFB,YOFMIN,YOFMAX,YOFDEF
! 13 May 2008 Scott Hannon   Add CO2NTE and NTEBOT for non-LTE calc and
!                            increase NNCOEF from 6 to 7
! 12 May 2009 Scott Hannon   Add VTUNNG string; delete VCLOUD
! 04 Aug 2009 Scott Hannon   Add FNSCAL for speccal text file
! 11 Oct 2011 Scott Hannon   Add MAXABW; replace FNSCAL with ABFILE

!END====================================================================
!
!-----------------------------------------------------------------------
!      IMPLICIT NONE
!-----------------------------------------------------------------------
! Note: having an "implicit none" in both the include file & the main
! source code will cause some compilers to complain.
!       IMPLICIT NONE


!-----------------------------------------------------------------------
!      INCLUDE FILES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXECUTABLE CODE
!-----------------------------------------------------------------------
!      none

!      -----------------------------------------------------------------
!      Assign SARTA version strings
!      -----------------------------------------------------------------
!      The version strings consists of 3 parts: version number, date,
!      and comment.  The version date should be updated to the
!      current date whenever any portion of the code is updated.  The
!      version number consists of two parts; a major version to the
!      left of the decimal point, and a minor version to the right.
!      The major number should be incremented only when major changes
!      have been made to the overall SARTA code.  The minor number
!      should be incremented only when minor but non-trivial changes
!      are made to the code.  Bug fixes should generally be handled
!      with the version date, but a fix for a serious bug may warrant
!      a change to the minor version number.
!      See the "Doc/last_update.txt" file for a description of the
!      changes associated with every change of VSARTA.
!
       CHARACTER*40 VSARTA  ! SARTA source code version
       CHARACTER*40 VSCOEF  ! SARTA coefficient version
       CHARACTER*40 VTUNNG  ! optical depth tuning version
!      version template    '#.## YYYY-MM-DD <--------comment------->'
       PARAMETER( VSARTA = '1.08 2011-11-21 PGEv6 AB20120210' )
       PARAMETER( VSCOEF = 'AIRS 2008-05-07 m130x m140x CO2=370' )
       PARAMETER( VTUNNG = 'wcon nte' )

!      *********
!      VARIABLES
!      *********
!      Note: these should not be changed by the user
!
!      ------------------------
!      Constants and other data
!      ------------------------
       REAL     PI ! pi, circle circumference/diameter (3.1415926)
       REAL RADSUN ! radius of the sun (6.956E+8 m)
       REAL     C1 ! radiation constant c1 (1.1911E-8  W/(m2.st.(cm-1)4)
       REAL     C2 ! radiation constant c2 (1.4387863 K/cm-1)
       PARAMETER(    PI = 3.1415926)
       PARAMETER(RADSUN = 6.956E+8)
!
!cc    Previously used values; agrees w/JPL pre-Dec2000
!cc    PARAMETER(    C1 = 1.1910439E-8)  ! JPL value is 1E+3 bigger
!cc    PARAMETER(    C2 = 1.4387687)
!
!      Current values (CODATA98 from NIST); agrees w/JPL Dec2000
       PARAMETER(  C1 = 1.191042722E-8)  ! JPL value is 1E+3 bigger
       PARAMETER(  C2 = 1.4387752)
!
       REAL CO2STD ! standard CO2 PPMV mixing ratio (370)
!      Note: mean global CO2 approximately 370 in 2002
       PARAMETER( CO2STD = 370.0 )
!
       REAL  XSALT ! expected nominal satellite altitude (km)
       PARAMETER( XSALT = 705.0 )
!
!
!      -----------------------------------
!      Channels and layers other variables
!      -----------------------------------
       INTEGER MAXLAY ! # of layers (100)
       INTEGER   NSET ! # of coefficient data sets (7)
       INTEGER MXCHAN ! max total # of channels (2378)
       INTEGER NFCOEF ! # of downwelling thermal "F" factor coefs
       INTEGER MXEMIS ! max # of input emis/rho data points
       INTEGER MAXPRO ! max # of user specified profiles
       INTEGER  MXGAS ! max # of gases in user profile
       INTEGER MXMIEA ! max # of mie particle sizes
       PARAMETER(MAXLAY = 100)
       PARAMETER(  NSET = 7)
       PARAMETER(MXCHAN = 2834)
       PARAMETER(NFCOEF = 6)
       PARAMETER(MXEMIS = 100)
       PARAMETER(MAXPRO = 25)
       PARAMETER( MXGAS = 44)
       PARAMETER(MXMIEA = 10)
!
!***********************************************************************
!      Variables for the coefficient sets
!***********************************************************************
!
!      --------------
!      For set1 = FWO
!      -------------
!      Used in part by modules: 12, 11, 10, 9, 8, 7, 6, 5, 3, 4b, 4a
       INTEGER MXCHN1 ! max # of channels for set1 = FWO (1461)
       INTEGER  N1CON ! # of water con predictors/coefs for set1 (5)
       INTEGER  N1FIX ! # of "fixed" predictors/coefs for set1 (8)
       INTEGER  N1H2O ! # of water predictors/coefs for set1 (13)
       INTEGER   N1O3 ! # of ozone predictors/coefs for set1 (5)
       INTEGER N1COEF ! total # of coefs for set1
       PARAMETER(MXCHN1 = 1461)
       PARAMETER( N1CON = 7)
       PARAMETER( N1FIX = 8)
       PARAMETER( N1H2O = 11)
       PARAMETER(  N1O3 = 5)
       PARAMETER(N1COEF = N1CON + N1FIX + N1H2O + N1O3 )
!
!
!      --------------
!      For set2 = FOW
!      --------------
!      Used in part by modules: 6, 5
       INTEGER MXCHN2 ! max # of channels for set2 = FOW  (325)
       INTEGER  N2CON ! # of water con predictors/coefs for set2 (5)
       INTEGER  N2FIX ! # of "fixed" predictors/coefs for set2 (8)
       INTEGER   N2O3 ! # of ozone predictors/coefs for set2 (10)
       INTEGER  N2H2O ! # of water predictors/coefs for set2 (11)
       INTEGER N2COEF ! total # of coefs for set2
       PARAMETER(MXCHN2 = 325)
       PARAMETER( N2CON = 7)
       PARAMETER( N2FIX = 8)
       PARAMETER(  N2O3 = 10)
       PARAMETER( N2H2O = 11)
       PARAMETER(N2COEF = N2CON + N2FIX + N2O3 + N2H2O )
!
!
!      --------------
!      For set3 = FMW
!      --------------
!      Used in part by modules: 4d, 4c, 3
       INTEGER MXCHN3 ! max # of channels for set3 = FMW  (396)
       INTEGER  N3CON ! # of water con predictors/coefs for set3 (5)
       INTEGER  N3FIX ! # of "fixed" predictors/coefs for set3 (8)
       INTEGER  N3CH4 ! # of methane predictors/coefs for set3 (9)
       INTEGER  N3H2O ! # of water predictors/coefs for set3 (13)
       INTEGER N3COEF ! total # of coefs for set3
       PARAMETER(MXCHN3 = 396)
       PARAMETER( N3CON = 7)
       PARAMETER( N3FIX = 8)
       PARAMETER( N3CH4 = 9)
       PARAMETER( N3H2O = 11)
       PARAMETER(N3COEF = N3CON + N3FIX + N3CH4 + N3H2O )
!
!
!      ---------------
!      For set4 = sun FCOW
!      ---------------
!      Used in part by modules: 2b
       INTEGER MXCHN4 ! max # of channels for set4 = FCOW (85)
       INTEGER  N4CON ! # of water con predictors/coefs for set4 (5)
       INTEGER  N4FIX ! # of "fixed" predictors/coefs for set4 (11)
       INTEGER   N4CO ! # of CO predictors/coefs for set4 (11)
       INTEGER   N4O3 ! # of ozone predictors/coefs for set4 (3)
       INTEGER  N4H2O ! # of water predictors/coefs for set4 (13)
       INTEGER N4COEF ! total # of coefs for set4
       PARAMETER(MXCHN4 = 85)
       PARAMETER( N4CON = 7)
       PARAMETER( N4FIX = 11)
       PARAMETER(  N4CO = 11)
       PARAMETER(  N4O3 = 3)
       PARAMETER( N4H2O = 13)
       PARAMETER(N4COEF = N4CON + N4FIX + N4CO + N4O3 + N4H2O )
!
!
!      -----------------------
!      For set5 = sun BFSW
!      -----------------------
!      Used in part by modules: 2b, 1b
       INTEGER MXCHN5 ! max # of channels for set5 = BFSW (210)
       INTEGER  N5CON ! # of water con predictors/coefs for set5 (5)
       INTEGER  N5FIX ! # of "fixed" predictors/coefs for set5 (11)
       INTEGER  N5H2O ! # of water predictors/coefs for set5 (3)
       INTEGER   N5O3 ! # of ozone predictors/coefs for set5 (1)
       INTEGER N5COEF ! total # of coefs for set5
       PARAMETER(MXCHN5 = 210)
       PARAMETER( N5CON = 7)
       PARAMETER( N5FIX = 11)
       PARAMETER( N5H2O = 3)
       PARAMETER(  N5O3 = 1)
       PARAMETER(N5COEF = N5CON + N5FIX + N5H2O + N5O3 )
!
!
!      -----------------------
!      For set6 = sun MFMW
!      -----------------------
!      Used in part by modules: 1b, 2a
       INTEGER MXCHN6 ! max # of channels for set6 = MFMW (217)
       INTEGER  N6CON ! # of water con predictors/coefs for set6 (5)
       INTEGER  N6FIX ! # of "fixed" predictors/coefs for set6 (8)
       INTEGER  N6H2O ! # of water predictors/coefs for set6 (7)
       INTEGER   N6O3 ! # of ozone predictors/coefs for set6 (1)
       INTEGER N6COEF ! total # of coefs for set6
       PARAMETER(MXCHN6 = 217)
       PARAMETER( N6CON = 7 )
       PARAMETER( N6FIX = 8 )
       PARAMETER( N6H2O = 7 )
       PARAMETER(  N6O3 = 1 )
       PARAMETER(N6COEF = N6CON + N6FIX + N6H2O + N6O3 )
!
!
!      -----------------------
!      For set7 = sun MFBW
!      -----------------------
!      Used in part by modules: 2a, 1a
       INTEGER MXCHN7 ! max # of channels for set7 = MFBW (140)
       INTEGER  N7CON ! # of water con predictors/coefs for set7 (5)
       INTEGER  N7FIX ! # of "fixed" predictors/coefs for set7 (8)
       INTEGER  N7H2O ! # of water predictors/coefs for set7 (13)
       INTEGER   N7O3 ! # of ozone predictors/coefs for set7 (1)
       INTEGER N7COEF ! total # of coefs for set7
       PARAMETER(MXCHN7 = 140)
       PARAMETER( N7CON = 7)
       PARAMETER( N7FIX = 8)
       PARAMETER( N7H2O = 13)
       PARAMETER(  N7O3 = 1)
       PARAMETER(N7COEF = N7CON + N7FIX + N7H2O + N7O3 )
!
!
!      ---------------
!      For trace gases predictors
!      ---------------
       INTEGER NTRACE ! number of trace gas perturbation predictors (7)
       PARAMETER(NTRACE = 7)
!
!
!      ----------------
!      For variable CO2
!      ----------------
!      Used in part by modules: 12, 11, 10, 9, 7, 6, 5, 2b, 1b, 2a
       INTEGER MXCHNC ! max # of channels with CO2 pert coefs (1082)
       INTEGER   NCO2 ! number of CO2 coefficients
       PARAMETER(MXCHNC = 1082)
       PARAMETER(  NCO2 = 5)
!
!
!      ----------------
!      For variable SO2
!      ----------------
       INTEGER MXCHNS ! max # of channels with SO2 pert coefs (602)
       INTEGER   NSO2 ! number of SO2 coefficients
       PARAMETER(MXCHNS = 602)
       PARAMETER(  NSO2 = 4)
!
!
!      -----------------
!      For variable HNO3
!      -----------------
       INTEGER MXCHNH ! max # of channels with HNO3 pert coefs (383)
       INTEGER  NHNO3 ! number of HNO3 coefficients
       PARAMETER(MXCHNH = 383)
       PARAMETER( NHNO3 = 4)
!
!
!      -----------------
!      For variable N2O
!      -----------------
       INTEGER MXCHNN ! max # of channels with N2O pert coefs (586)
       INTEGER   NN2O ! number of N2O coefficients
       PARAMETER(MXCHNN = 586)
       PARAMETER(  NN2O = 7)
!
!
!      ----------------------
!      For OPTRAN water coefs
!      ----------------------
!      Used in part by modules:
       INTEGER MXCHNW ! max # of channelss with OPTRAN H2O coefs (754)
       INTEGER MXOWLY ! number of OPTRAN water layers
       INTEGER NOWAVG ! # of OPTRAN water average profile values (4)
       INTEGER NH2O   ! number of OPTRAN H2O predictors/coefs (9)
       PARAMETER(MXCHNW = 754)
       PARAMETER(MXOWLY = 300)
       PARAMETER(NOWAVG = 4)
       PARAMETER(  NH2O = 9)
!
!      -----------
!      For non-LTE
!      -----------
       INTEGER MXCNTE ! max # of channels for non-LTE (203)
       INTEGER NNCOEF ! # of coefs for non-LTE
       INTEGER NTEBOT ! bottom layer for CO2TOP calc
       REAL CO2NTE ! ref CO2 mixing ratio for non-LTE coefs (ppmv)
       PARAMETER(MXCNTE = 203)
       PARAMETER(NNCOEF = 7)
       PARAMETER(NTEBOT = 10)
       PARAMETER(CO2NTE = 370.0)
!
!      ----------------
!      Single filenames
!      ----------------
       CHARACTER*80 FNFX   ! coef fx
       CHARACTER*80 FNPREF ! reference profile
       PARAMETER(FNFX  = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/fx.txt')
       PARAMETER(FNPREF= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/profref_trace370')
!
!      ----------------
!      Double filenames
!      ----------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       REAL YOFFA          ! yoffset for database "A"
       CHARACTER*80 FACOF1 ! coef set1
       CHARACTER*80 FACOF2 ! coef set2
       CHARACTER*80 FACOF3 ! coef set3
       CHARACTER*80 FACOF4 ! coef set4
       CHARACTER*80 FACOF5 ! coef set5
       CHARACTER*80 FACOF6 ! coef set6
       CHARACTER*80 FACOF7 ! coef set7
       CHARACTER*80 FACO2  ! coef CO2
       CHARACTER*80 FAN2O  ! coef N2O
       CHARACTER*80 FASO2  ! coef SO2
       CHARACTER*80 FAHNO3 ! coef HNO3
       CHARACTER*80 FAOPTR ! coef optran
       CHARACTER*80 FATHER ! coef therm
       CHARACTER*80 FACOFN ! coef non-LTE
       CHARACTER*80 FASUN  ! solar data
       CHARACTER*80 FATMLT ! tuning multiplier
!
       PARAMETER(YOFFA=-13.0)
       PARAMETER(FACOF1= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set1_m130x370.dat')
       PARAMETER(FACOF2= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set2_m130x370.dat')
       PARAMETER(FACOF3= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set3_m130x370.dat')
       PARAMETER(FACOF4= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set4_m130x370.dat')
       PARAMETER(FACOF5= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set5_m130x370.dat')
       PARAMETER(FACOF6= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set6_m130x370.dat')
       PARAMETER(FACOF7= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set7_m130x370.dat')
       PARAMETER(FACO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'CO2_5term_m130x370.dat')
       PARAMETER(FAN2O = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'N2O_m130x370.dat')
       PARAMETER(FASO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'SO2_m130x370.dat')
       PARAMETER(FAHNO3 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'HNO3_m130x370.dat')
       PARAMETER(FAOPTR= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'optran_m130x370.dat')
       PARAMETER(FATHER= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'therm_m130x370.dat')
       PARAMETER(FACOFN= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'nonLTE7_m130x.dat')
       PARAMETER(FASUN = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Solar/solar_m130x.txt')
!
       PARAMETER(FATMLT= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'tunmlt_wcon_nte.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_ones.txt')
!cc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
       REAL YOFFB          ! yoffset for database "B"
       CHARACTER*80 FBCOF1 ! coef set1
       CHARACTER*80 FBCOF2 ! coef set2
       CHARACTER*80 FBCOF3 ! coef set3
       CHARACTER*80 FBCOF4 ! coef set4
       CHARACTER*80 FBCOF5 ! coef set5
       CHARACTER*80 FBCOF6 ! coef set6
       CHARACTER*80 FBCOF7 ! coef set7
       CHARACTER*80 FBCO2  ! coef CO2
       CHARACTER*80 FBN2O  ! coef N2O
       CHARACTER*80 FBSO2  ! coef SO2
       CHARACTER*80 FBHNO3 ! coef HNO3
       CHARACTER*80 FBOPTR ! coef optran
       CHARACTER*80 FBTHER ! coef therm
       CHARACTER*80 FBCOFN ! coef non-LTE
       CHARACTER*80 FBSUN  ! solar data
       CHARACTER*80 FBTMLT ! tuning multiplier
!
       PARAMETER(YOFFB=-14.0)
       PARAMETER(FBCOF1= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set1_m140x370.dat')
       PARAMETER(FBCOF2= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set2_m140x370.dat')
       PARAMETER(FBCOF3= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set3_m140x370.dat')
       PARAMETER(FBCOF4= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set4_m140x370.dat')
       PARAMETER(FBCOF5= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set5_m140x370.dat')
       PARAMETER(FBCOF6= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set6_m140x370.dat')
       PARAMETER(FBCOF7= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'set7_m140x370.dat')
       PARAMETER(FBCO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'CO2_5term_m140x370.dat')
       PARAMETER(FBN2O = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'N2O_m140x370.dat')
       PARAMETER(FBSO2 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'SO2_m140x370.dat')
       PARAMETER(FBHNO3 = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'HNO3_m140x370.dat')
       PARAMETER(FBOPTR= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'optran_m140x370.dat')
       PARAMETER(FBTHER= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'therm_m140x370.dat')
       PARAMETER(FBCOFN= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'nonLTE7_m140x.dat')
       PARAMETER(FBSUN = &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Solar/solar_m140x.txt')
!
       PARAMETER(FBTMLT= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'tunmlt_wcon_nte.txt')
!cc
!     $ '/asl/data/sarta_database/Data_AIRS_apr08/Coef/'
!     $ // 'tunmlt_ones.txt')
!cc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       CHARACTER*80 ABFILE ! spectral calibration
       PARAMETER(ABFILE= &
       trim(sarta_dir_data)//'Data_AIRS_apr08/Coef/' &
       // 'ABweight_table_20120210.txt')
!
!      Number of A/B weight epochs
!      Note: must match number in ABFILE
       INTEGER MAXABW
       PARAMETER(MAXABW=6)
!
!      ----------------------
!      Default & min/max YOFF
!      ----------------------
!      Note: these should be based YOFFA and YOFFB as well as the
!      expected range of variation of YOFF
       REAL YOFDEF   ! default YOFF
       REAL YOFMIN   ! min allowed YOFF
       REAL YOFMAX   ! max allowed YOFF
       PARAMETER(YOFDEF=-13.5)
       PARAMETER(YOFMIN=-16.0)
       PARAMETER(YOFMAX=-11.0)

!      ----------------
!      I/O unit numbers
!      ----------------
!      Note: these units are not explicitly openned by the sarta code,
!      they should be set to standard I/O units for your compiler
       INTEGER IOINFO  ! unit number for non-error info messages (6)
       INTEGER IOERR   ! unit number for error messages (2 or 6)
       PARAMETER( IOINFO = 6 )
       PARAMETER( IOERR = 0 )
!
!
!      -----------------
!      Allowed input GUC (Gas Units Code number)
!      -----------------
       INTEGER GUCIN  ! The one & only allowed input GUC number
!      Note: GUCIN must be 1 or 2.  All gases in the input RTP
!      must be of this type.
!       PARAMETER( GUCIN = 1 ) ! GUC number for:  molecules/cm^2
!       PARAMETER( GUCIN = 2 ) ! GUC number for:  kilomoles/cm^2
        PARAMETER( GUCIN = 3 ) ! GUC number for:  ppmv

!      End of include file

#endif


!=======================================================================
!
! RTP Fortran API structures and parameters
! Version 2.01
! The Fortran structures defined here, RTPHEAD, RTPPROF, and
! RTPATTR, must match the corresponding C structures rtp_head,
! rtp_prof, and rtpfatt, and the parameters set below must have
! the same values as the corresponding C #define parameters.
!
! Note: total record size for header or profile records
! may not exceed 50 kB
!
! See rtp.h and rtpspec.pdf for more information on the fields
! defined below.
!
! --------------
! RTP parameters
! --------------
        integer BAD ! value to be used if no data
        integer LEVPRO ! levels-type profile flag
        integer LAYPRO ! layers-type profile flag
        integer AIRSLAY ! AIRS layers-type profile flag
        integer PROFBIT ! Summed radiance bit flags
        integer IRCALCBIT ! bit flag for calculated radiance
        integer IROBSVBIT ! bit flag for observed radiance
        integer PFIELDSMAX ! max allowed value of PROFBIT
        integer MAXEMIS ! max num of emissivity/rho points
        integer MAXGAS ! max number of gases
        integer MAXGASID ! max gas ID value
        integer MAXLEV ! max number of levels
        integer MAXCHAN ! max number of channels
        integer MAXPNOTE ! max profile comment string
        integer MAXUDEF ! max profile udef values
        integer MAXIUDEF ! max profile and header iudef values
        integer MAXOPEN ! max num of open RTP files
        integer MAXNATTR ! max number of attributes
        integer MAXVNAME ! max num of chars in field or vdata name
       integer MAXANAME ! max num of chars in attribute name
        integer MAXATEXT ! max num of chars in attribute text
        integer MAXCALF ! max cal flag bytes, ceil(MAXCHAN/4)*4
        integer MAXPN4 ! true max for pnote, ceil(MAXPNOTE/4)*4
        ! the following parameters must match the values set in rtp.h
        !
        parameter ( BAD = -9999 )
        parameter ( LEVPRO = 0 )
        parameter ( LAYPRO = 1 )
        parameter ( AIRSLAY = 2 )
        parameter ( PROFBIT = 1 )
        parameter ( IRCALCBIT = 2 )
        parameter ( IROBSVBIT = 4 )
        parameter ( PFIELDSMAX = 7 )
        parameter ( MAXEMIS = 100 )
        parameter ( MAXGAS = 16 )
        parameter ( MAXGASID = 303 )
        parameter ( MAXLEV = 120 )
        parameter ( MAXCHAN = 4231 )
        parameter ( MAXPNOTE = 80 )
        parameter ( MAXUDEF = 20 )
        parameter ( MAXIUDEF = 10 )
        parameter ( MAXOPEN = 8 )
        parameter ( MAXNATTR = 32 )
        parameter ( MAXCALF = ((MAXCHAN-1)/4+1)*4 )
        parameter ( MAXPN4 = ((MAXPNOTE-1)/4+1)*4 )
        ! the following parameters must match the values set in pvdefs.h
        ! and also the field sizes in the RTPATTR structure, defined below
        !
        parameter ( MAXVNAME = 64 )
        parameter ( MAXANAME = 64 )
        parameter ( MAXATEXT = 1024 )
! --------------------
! RTP header structure
! --------------------
!
! STRUCTURE /RTPHEAD/  --> now F90 derived type
    type RTPHEAD
        ! profile data
          integer*4 ptype ! profile type
          integer*4 pfields ! profile field set
          real*4 pmin ! min plevs value (not used)
          real*4 pmax ! max plevs value (not used)
          integer*4 ngas ! number of gases
          integer*4 glist(MAXGAS) ! constituent gas list
          integer*4 gunit(MAXGAS) ! constituent gas units
        ! radiance data
          integer*4 pltfid ! platform ID/code number
          integer*4 instid ! instrument ID/code number
          integer*4 nchan ! number of channels
          integer*4 ichan(MAXCHAN) ! channel ID numbers
          real*4 vchan(MAXCHAN) ! channel center freq.
          real*4 vcmin ! chan set min freq, including wings (not used)
          real*4 vcmax ! chan set max freq, including wings (not used)
        ! maxes for profile fields
          ! these fields are not saved explicitly in the HDF file
          integer*4 memis ! max number of emis/rho points
          integer*4 mlevs ! max number of pressure level
        ! user defined fields
          integer*4 iudef(MAXIUDEF) ! user-defined integer array
          integer*4 itype ! user-defined integer
! END STRUCTURE
  end type RTPHEAD
! ---------------------
! RTP profile structure
! ---------------------
!
! STRUCTURE /RTPPROF/  --> now F90 derived type
      type RTPPROF
        ! profile location/time
          real*4 plat ! profile latitude
          real*4 plon ! profile longitude
          real*8 ptime ! profile time
          ! surface data
          real*4 stemp ! surface temperature
          real*4 salti ! surface altitude
          real*4 spres ! surface pressure
          real*4 landfrac ! land fraction
          integer*4 landtype ! land type code
          real*4 wspeed ! wind speed
          integer*4 nemis ! number of emis. pts
          real*4 efreq(MAXEMIS) ! emissivity freq's
          real*4 emis(MAXEMIS) ! surface emissivities
          real*4 rho(MAXEMIS) ! surface reflectance
          ! atmospheric data
          integer*4 nlevs ! number of press levels
          real*4 plevs(MAXLEV) ! pressure levels
          real*4 palts(MAXLEV) ! level altitudes
          real*4 ptemp(MAXLEV) ! temperature profile
          real*4 gamnt(MAXLEV,MAXGAS) ! gas amounts
          real*4 gtotal(MAXGAS) ! total column gas amount
          real*4 gxover(MAXGAS) ! gas crossover press
          real*4 txover ! temperature crossover press
          real*4 co2ppm ! CO2 mixing ratio
          ! clear flag/code
          integer*4 clrflag ! clear flag/code
          ! cloud1 data
          integer*4 ctype ! cloud type code
          real*4 cfrac ! cloud fraction
          real*4 cemis(MAXEMIS) ! cloud top emissivity
          real*4 crho(MAXEMIS) ! cloud top reflectivity
          real*4 cprtop ! cloud top pressure
          real*4 cprbot ! cloud bottom pressure
          real*4 cngwat ! cloud non-gas water
          real*4 cpsize ! cloud particle size
          real*4 cstemp ! cloud surface temperature
          ! cloud2 data
          integer*4 ctype2 ! cloud2 type code
          real*4 cfrac2 ! cloud2 fraction
          real*4 cemis2(MAXEMIS) ! cloud2 top emissivity
          real*4 crho2(MAXEMIS) ! cloud2 top reflectivity
          real*4 cprtop2 ! cloud2 top pressure
          real*4 cprbot2 ! cloud2 bottom pressure
          real*4 cngwat2 ! cloud2 non-gas water
          real*4 cpsize2 ! cloud2 particle size
          real*4 cstemp2 ! cloud2 surface temperature
          real*4 cfrac12 ! cloud1+2 fraction
        ! radiance orientation data
          real*4 pobs ! observation pressure
          real*4 zobs ! observation height
          integer*4 upwell ! radiation direction
          real*4 scanang ! scan angle
          real*4 satzen ! satellite zenith angle
          real*4 satazi ! satellite azimuth angle
          ! sun info
          real*4 solzen ! sun zenith angle
          real*4 solazi ! sun azimuth angle
          real*4 sundist ! Earth-Sun distance
          real*4 glint ! glint distance or flag
          ! observation location/time
          real*4 rlat ! radiance obs lat.
          real*4 rlon ! radiance obs lon.
! integer*4 rfill ! align rtime on 8 byte bndry
          real*8 rtime ! radiance obs time
          ! observation indices
          integer*4 findex ! file (granule) index
          integer*4 atrack ! along-track index
          integer*4 xtrack ! cross-track index
          integer*4 ifov ! field of view index
        ! observed radiance data
          real*4 robs1(MAXCHAN) ! obs radiance
          character*1 calflag(MAXCALF) ! obs rad per chan calib/qual flags
          integer*4 robsqual ! obs rad overall quality flag/code
          real*4 freqcal ! frequency calibration
        ! calculated radiance data
          real*4 rcalc(MAXCHAN) ! calc radiance
        ! user defined fields
          character*80 pnote ! profile annotation, size MAXPN4
          real*4 udef(MAXUDEF) ! user-defined real array
          integer*4 iudef(MAXIUDEF) ! user-defined integer array
          integer*4 itype ! user0defined integer
! END STRUCTURE
       end type RTPPROF
! -----------------------
! RTP attribute structure
! -----------------------
!
! fname is the name of the field the attribute is to be associated
! with, 'header' for a general header attribute, or 'profiles'
! for a general profile attribute. Its size declaration should
! be the same as the parameter MAXVNAME.
!
! aname is the attribute name, e.g., 'units' for a field attribute,
! or 'TITLE', for a general header attribute. Its size should
! also be the same as the parameter MAXANAME.
!
! atext is the attribute text, 'e.g., '48 Fitting Profiles' might be
! the atext of the header 'TITLE' attribute. Its size should be
! the same as the parameter MAXATEXT.
!
! STRUCTURE /RTPATTR/
        type RTPATTR
          character*64 fname ! associated field name
          character*64 aname ! attribute name
          character*1024 atext ! attribute text
        end type RTPATTR
! END STRUCTURE


!
! gas 
!
 integer,parameter :: ngas    = 8
 integer,dimension(ngas),parameter :: glist   = (/1,2,3,4,5,6,9,12/)
 integer,dimension(ngas),parameter :: gunit   = (/3,3,3,3,3,3,3,3/) !ppmv coupling with G-SDSU

!
! emissivity
!
 integer :: nemis_ap        ! number of emissivity anchor point 
 real :: efreq_ap(MAXEMIS)  ! frequency of emissivity anchor point
 real :: emis_ap (MAXEMIS)  ! emissivity at anchor point
 real :: rho_ap  (MAXEMIS)  ! albedo at anchor point



       REAL PLEV(101)
       INTEGER :: IP
!
       DATA ( PLEV(IP), IP = 1, 50 ) &
               / 0.0050, 0.0161, 0.0384, 0.0769, 0.1370, &
                    0.2244, 0.3454, 0.5064, 0.7140, 0.9753, &
                    1.2972, 1.6872, 2.1526, 2.7009, 3.3398, &
                    4.0770, 4.9204, 5.8776, 6.9567, 8.1655, &
                    9.5119, 11.0038, 12.6492, 14.4559, 16.4318, &
                   18.5847, 20.9224, 23.4526, 26.1829, 29.1210, &
                   32.2744, 35.6505, 39.2566, 43.1001, 47.1882, &
                   51.5278, 56.1260, 60.9895, 66.1253, 71.5398, &
                   77.2396, 83.2310, 89.5204, 96.1138, 103.0172, &
                  110.2366, 117.7775, 125.6456, 133.8462, 142.3848/
       DATA ( PLEV(IP), IP = 51, 101 ) &
               / 151.2664, 160.4959, 170.0784, 180.0183, 190.3203, &
                  200.9887, 212.0277, 223.4415, 235.2338, 247.4085, &
                  259.9691, 272.9191, 286.2617, 300.0000, 314.1369, &
                  328.6753, 343.6176, 358.9665, 374.7241, 390.8926, &
                  407.4738, 424.4698, 441.8819, 459.7118, 477.9607, &
                  496.6298, 515.7200, 535.2322, 555.1669, 575.5248, &
                  596.3062, 617.5112, 639.1398, 661.1920, 683.6673, &
                  706.5654, 729.8857, 753.6275, 777.7897, 802.3714, &
                  827.3713, 852.7880, 878.6201, 904.8659, 931.5236, &
                  958.5911, 986.0666, 1013.9476, 1042.2319, 1070.9170, &
                 1100.0000 /
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALNTE (for Non-local Thermodynamic Equilibrium)
!
!F77====================================================================
!ROUTINE NAME:
! CALNTE
!ABSTRACT:
! Adjust a LTE atmospheric radiance for a non-LTE upper atmosphere.
!CALL PROTOCOL:
! CALNTE( INDCHN, TEMP, SUNCOS, SCOS1, VSEC1,
! NCHNTE, CLISTN, COEFN, CO2TOP, RAD)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INT arr INDCHN channel indices none
! REAL arr TEMP temperature profile Kelvin
! REAL SUNCOS solzen cosine at surface none
! REAL SCOS1 solzen cosine at layer1 none
! REAL VSEC1 satzen secant at layer1 none
! INTEGER NCHNTE number of non-LTE channels none
! INT arr CLISTN non-LTE channel list none
! REAL arr COEFN non-LTE coefficients various
! REAL arr CO2TOP top layers CO2 mixing ratio ppmv
!OUTPUT PARAMETERS:
! none
!INPUT/OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr RAD radiance W/(m^2.str.cm^-1)
!RETURN VALUES:
! none
!PARENT(S):
! SARTA
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! May 2008 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! In the upper atmosphere where the air is thin, the strong CO2
! absoption bands in the 4 um region can absorb solar radiance
! faster than collisons with other air molecules can re-distribute
! the energy. The CO2 is no longer in thermodynamic equilibrium
! with its surroundings, which results in a change to the CO2
! vibrational band population statistics and its effective
! radiating temperature. This code applies a regression based
! adjustment to the input radiance to account for non-LTE effects.
! Coefficients and predictors are multiplied together and summed
! to calculate the change in radiance for non-LTE conditions.
!
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ------------------------------------------
! 15 Mar 2005 Scott Hannon Created
! 13 Oct 2005 Scott Hannon MXCHNN renamed MXCNTE to avoid conflict
! with MXCHNN used with N2O
! 14 May 2008 Scott Hannon Add CO2 adjustment using 7th coef; pass in
! CO2TOP
!END====================================================================
! =================================================================
       SUBROUTINE CALNTE ( INDCHN, FREQ, TEMP, SUNCOS, SCOS1, VSEC1, &
          NCHNTE, CLISTN, COEFN, CO2TOP, RAD , BT)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! QIKEXP
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER INDCHN(MXCHAN)
       REAL,INTENT(IN) :: FREQ(MXCHAN)

       REAL TEMP(MAXLAY)
       REAL SUNCOS ! solar zenith angle cosine at surface
       REAL SCOS1 ! solar zenith angle cosine at layer1
       REAL VSEC1 ! satellite view zenith angle secant at layer1
       INTEGER NCHNTE
       INTEGER CLISTN(MXCNTE)
       REAL COEFN(NNCOEF,MXCNTE)
       REAL CO2TOP ! CO2 mixing ratio in top layers (ppmv)
!
! Input/Output
       REAL ,INTENT(INOUT) :: RAD(MXCHAN)
       REAL ,INTENT(INOUT) :: BT(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER J
       REAL DRAD
       REAL PRED1
       REAL PRED2
       REAL PRED3
       REAL PRED4
       REAL PRED5
       REAL PRED6
       REAL THIGH
 REAL ::  C1V3
 REAL ::  C2V

!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Calculate the channel independent non-LTE predictors
       THIGH = (TEMP(1) + TEMP(2) + TEMP(3) + TEMP(4) + TEMP(5))/5.0
       PRED1 = 1.0
       PRED2 = SCOS1
       PRED3 = SCOS1*SCOS1
       PRED4 = SCOS1*VSEC1
       PRED5 = SCOS1*THIGH
       PRED6 = SUNCOS
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO I=1,NCHNTE
!
! Index for RAD
          J=INDCHN( CLISTN(I) )
!
          DRAD=( COEFN(1,I)*PRED1 ) + &
               ( COEFN(2,I)*PRED2 ) + &
               ( COEFN(3,I)*PRED3 ) + &
               ( COEFN(4,I)*PRED4 ) + &
               ( COEFN(5,I)*PRED5 ) + &
               ( COEFN(6,I)*PRED6 )
!
! Adjust DRAD for CO2 mixing ratio
          DRAD=DRAD*(COEFN(7,I)*(CO2TOP - CO2NTE) + 1.0)
!
! Adjust RAD for the non-LTE contribution
          RAD(J) = RAD(J) + DRAD/1000.0 ! convert DRAD to Watts

!
! Calc c1*v^3 and c2*v
          C1V3=C1*(FREQ(J)**3)
          C2V=C2*FREQ(J)
! Compute Brightness Temperature
          BT(J)=C2V/LOG( 1.0 + C1V3/RAD(J) )


!
!
       ENDDO
! End loops on channel number (frequency)
!
       RETURN
       END SUBROUTINE CALNTE
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALOKW
!
!F77====================================================================
!ROUTINE NAME:
! CALOKW
!ABSTRACT:
! Calculate the OPTRAN derived water pressure layer effective
! optical depth for a single channel.
!CALL PROTOCOL:
! CALOKW ( LBOT, ICHAN, LOPMIN, LOPMAX, LOPLOW, LOPUSE,
! H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER LBOT bottom pres layer number none
! INTEGER ICHAN OPTRAN water channel index none
! INTEGER LOPMIN min OPTRAN level to use none
! INTEGER LOPMAX max OPTRAN level to use none
! INTEGER LOPLOW low OPTRAN bracketing lev none
! LOG arr LOPUSE Need this OPTRAN level? none
! REAL arr H2OPRD OPTRAN water predictors various
! REAL arr COFH2O OPTRAN H2O fast trans coef various
! REAL arr WAOP OPTRAN layer water amounts kiloMoles/cm^2
! REAL arr DAOP OPTRAN-to-AIRS interp fact none
! REAL arr WAANG AIRS layer water amounts kiloMoles/cm^2
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr KW AIRS H2O layer eff op dep none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! CALT1
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! March 1998 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The OPTRAN predictors and fast transmittance coefficients are
! used to calculate the water effective absorption coefficient on
! on the OPTRAN level grid. Only the OPTRAN levels actually needed
! are calculated.
! Note: the COFH2O*H2OPRD result must be divided by WAOP, a scaling
! factor which was originally applied during the fast transmittance
! coefficient regression.
! The OPTRAN absorption coefficients are then interpolated onto the
! 100 AIRS layers and multiplied by the AIRS layer water amount (to
! convert absorption coefficient into optical depth).
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- --------------------------------------
! 27 Feb 1998 Scott Hannon Created
! 26 Aug 1998 Scott Hannon Add LBOT to call; loop on LBOT instead
! of MAXLAY
!END====================================================================
! =================================================================
       SUBROUTINE CALOKW ( LBOT, ICHAN, LOPMIN, LOPMAX, LOPLOW, LOPUSE, &
          H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER LBOT
       INTEGER ICHAN
       INTEGER LOPMIN
       INTEGER LOPMAX
       INTEGER LOPLOW(MAXLAY)
       LOGICAL LOPUSE(MXOWLY)
       REAL DAOP(MAXLAY)
       REAL WAANG(MAXLAY)
       REAL WAOP(MXOWLY)
       REAL H2OPRD( NH2O,MXOWLY)
       REAL COFH2O( NH2O,MXOWLY,MXCHNW)
! Output
       REAL KW(MAXLAY)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER L
       INTEGER LOP
       REAL KWOP(MXOWLY)
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ---------------------------------
! Loop over the OPTRAN water levels
! ---------------------------------
! Only do calc for OPTRAN levels that are needed
       DO LOP=LOPMIN,LOPMAX
          IF (LOPUSE(LOP)) THEN
             KWOP(LOP)= &
                COFH2O(1,LOP,ICHAN)*H2OPRD(1,LOP) + &
                COFH2O(2,LOP,ICHAN)*H2OPRD(2,LOP) + &
                COFH2O(3,LOP,ICHAN)*H2OPRD(3,LOP) + &
                COFH2O(4,LOP,ICHAN)*H2OPRD(4,LOP) + &
                COFH2O(5,LOP,ICHAN)*H2OPRD(5,LOP) + &
                COFH2O(6,LOP,ICHAN)*H2OPRD(6,LOP) + &
                COFH2O(7,LOP,ICHAN)*H2OPRD(7,LOP) + &
                COFH2O(8,LOP,ICHAN)*H2OPRD(8,LOP) + &
                COFH2O(9,LOP,ICHAN)*H2OPRD(9,LOP)
! Remove WAOP scaling factor
             KWOP(LOP)=KWOP(LOP)/WAOP(LOP)
! Check for negative value
             IF (KWOP(LOP) .LT. 0.0E+0) KWOP(LOP)=0.0E+0
          ENDIF
       ENDDO
!
! -------------------------
! Loop over the AIRS layers
! -------------------------
       DO L=1,LBOT
!
! Interpolate abs coef and convert to optical depth
          KW(L)=( DAOP(L)*( KWOP(LOPLOW(L) + 1) - &
             KWOP(LOPLOW(L)) ) + KWOP(LOPLOW(L)) )*WAANG(L)
          IF (KW(L) .LT. 0.0E+0) KW(L)=0.0E+0
!
       ENDDO
!
       RETURN
       END SUBROUTINE CALOKW
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
! This version uses amounts at center of layer, not lower boundary
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALOWP
!
!F77====================================================================
!ROUTINE NAME:
! CALOWP
!ABSTRACT:
! Calculate the OPTRAN water (H2O) predictors for a profile.
!CALL PROTOCOL:
! CALOWP ( LBOT, WAMNT, P, T, SECANG, WAZOP, WAVGOP,
! WAANG, LOPMIN, LOPMAX, LOPUSE, H2OPRD, LOPLOW, DAOP )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER LBOT bottom layer number none
! REAL arr WAMNT profile layer water kiloMoles/cm^2
! REAL arr P layer pressures atmospheres
! REAL arr T profile temperature K
! REAL arr SECANG secant of path angle none
! REAL arr WAZOP OPTRAN l-to-s water grid kiloMoles/cm^2
! REAL arr WAVGOP OPTRAN average preds various
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr WAANG water amount in layer kiloMoles/cm^2
! INTEGER LOPMIN min OPTRAN level to use none
! INTEGER LOPMAX max OPTRAN level to use none
! LOG arr LOPUSE OPTRAN level needed? none
! REAL arr H2OPRD OPTRAN predictors various
! INTEGER LOPLOW low bracketing OPTRAN lev none
! REAL arr DAOP OPTRAN-to-AIRS interp frac none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! March 1998 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! Assumes the user has supplied vaguely realistic profile amounts
! and temperatures.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- --------------------------------------
! 27 Feb 1998 Scott Hannon Created
! 26 Aug 1998 Scott Hannon Add LBOT to call; loop on LBOT instead
! of MAXLAY
!END====================================================================
! =================================================================
       SUBROUTINE CALOWP ( LBOT, WAMNT, P, T, SECANG, WAZOP, WAVGOP, &
          WAANG, LOPMIN, LOPMAX, LOPUSE, H2OPRD, LOPLOW, DAOP )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER LBOT
       REAL WAMNT(MAXLAY)
       REAL P(MAXLAY)
       REAL T(MAXLAY)
       REAL SECANG(MAXLAY)
       REAL WAZOP(MXOWLY)
       REAL WAVGOP(NOWAVG,MXOWLY)
!
! Output
       REAL WAANG(MAXLAY)
       INTEGER LOPMIN
       INTEGER LOPMAX
       REAL H2OPRD( NH2O,MXOWLY)
       LOGICAL LOPUSE(MXOWLY)
       INTEGER LOPLOW(MAXLAY)
       REAL DAOP(MAXLAY)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER L
       INTEGER LL
       INTEGER LOP
       INTEGER LOPL
       INTEGER LOPU
       INTEGER LU
       REAL WAZ(MAXLAY)
       REAL WAZSUM
       REAL WPZSUM
       REAL WTZSUM
       REAL PZ(MAXLAY)
       REAL TZ(MAXLAY)
       REAL DA
       REAL POP
       REAL TOP
       REAL PZOP
       REAL TZOP
       REAL ANGOP
       LOGICAL LAST
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Initialize amount above sums
       WAZSUM=0.0E+0
       WPZSUM=0.0E+0
       WTZSUM=0.0E+0
!
! ---------------------------------------
! Calculate raw predictors for all layers
! ---------------------------------------
       DO L=1,LBOT
!
! Layer amount*angle
          WAANG(L)=WAMNT(L)*SECANG(L)
!
! Center-of-layer-to-space amount*angle
! Note: do this before updating AZSUM
          WAZ(L)=5.0E-1*WAANG(L) + WAZSUM
!
! Bottom-of-layer-to-space amount sum
          WAZSUM=WAANG(L) + WAZSUM
!
! Pressure above sum
          WPZSUM=WAANG(L)*P(L) + WPZSUM
          PZ(L)=WPZSUM/WAZSUM
!
! Temperature above sum
          WTZSUM=WAANG(L)*T(L) + WTZSUM
          TZ(L)=WTZSUM/WAZSUM
!
       ENDDO
!
! --------------------------------------------------
! Find the max OPTRAN level that is less than WAZ(1)
! --------------------------------------------------
       LOPMIN=1
 30 IF (WAZOP(LOPMIN+1) .LT. WAZ(1)) THEN
          LOPMIN=LOPMIN + 1
          GOTO 30
       ENDIF
!
! Initialize the upper and lower (pressure) layer index
       LL=1
       LU=2
       LAST=.FALSE.
!
! ----------------------------------------
! Loop over the OPTRAN layers (while loop)
! ----------------------------------------
       LOP=LOPMIN
 10 IF (LOP .LE. MXOWLY) THEN
!
! --------------------------------------------------------
! Find the two pressure layers closest to the OPTRAN layer
! --------------------------------------------------------
 20 IF (WAZ(LU) .LT. WAZOP(LOP)) THEN
             IF (LU .LT. LBOT) THEN
                LL=LU
                LU=LU + 1
                GOTO 20
             ELSE
                LAST=.TRUE.
             ENDIF
          ENDIF
!
! Compute the interpolation fractor
          DA=(WAZOP(LOP) - WAZ(LL))/(WAZ(LU) - WAZ(LL))
!
! Do the interpolation
          POP=( DA*( P(LU) - P(LL) ) + P(LL) )/WAVGOP(1,LOP)
          TOP=( DA*( T(LU) - T(LL) ) + T(LL) )/WAVGOP(2,LOP)
          PZOP=( DA*( PZ(LU) - PZ(LL) ) + PZ(LL) )/WAVGOP(3,LOP)
          TZOP=( DA*( TZ(LU) - TZ(LL) ) + TZ(LL) )/WAVGOP(4,LOP)
          ANGOP=DA*( SECANG(LU) - SECANG(LL) ) + SECANG(LL)
!
! Assign the predictors
          H2OPRD(1,LOP)=1.0E+0
          H2OPRD(2,LOP)=POP
          H2OPRD(3,LOP)=TOP
          H2OPRD(4,LOP)=SQRT( POP )
          H2OPRD(5,LOP)=TOP**2
          H2OPRD(6,LOP)=POP*TOP
          H2OPRD(7,LOP)=ANGOP
          H2OPRD(8,LOP)=PZOP
          H2OPRD(9,LOP)=TZOP
!
! Update LOP and loop
          IF (LAST .EQV. .TRUE.) THEN
             LOPMAX=LOP
! Set LOP > MXOWLY to exit loop over LOP
             LOP=MXOWLY + 1
          ELSE
             LOP=LOP + 1
          ENDIF
          GOTO 10
!
       ENDIF
! End while loop over LOP
!
! -----------------
! Initialize LOPUSE
! -----------------
       DO LOP=1,MXOWLY
          LOPUSE(LOP)=.FALSE.
       ENDDO
!
! ---------------------------------------
! Determine what OPTRAN layers are needed
! ---------------------------------------
! Initialize LOPL and LOPU
       LOPL=LOPMIN
       LOPU=LOPMIN + 1
!
! Loop over the AIRS pressure layers
       DO L=1,LBOT
! Find the two OPTRAN levels that bracket the AIRS layer
 40 IF (WAZOP(LOPU) .LT. WAZ(L) .AND. LOPU .LT. LOPMAX) THEN
             LOPL=LOPU
             LOPU=LOPU + 1
             GOTO 40
          ENDIF
!
          LOPUSE(LOPL)=.TRUE.
          LOPUSE(LOPU)=.TRUE.
! Assign the lower OPTRAN level
          LOPLOW(L)=LOPL
! Assign the interpolation fraction
          DAOP(L)=(WAZ(L) - WAZOP(LOPL))/(WAZOP(LOPU) - WAZOP(LOPL))
       ENDDO
!
       RETURN
       END SUBROUTINE CALOWP
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALPAR
!
!F77====================================================================
!ROUTINE NAME:
! CALPAR
!ABSTRACT:
! Calculate the fast transmittance code temperature/amount/angle
! dependent variables for a profile.
!CALL PROTOCOL:
! CALPAR (LBOT, RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RNAMNT,
! $ PTEMP,PFAMNT,PWAMNT,POAMNT,PCAMNT,PMAMNT,PNAMNT,
! $ RES,SECANG, ALAT, FX, DZREF,
! $ LCO2, LN2O, LSO2, LHNO3,LCO2PM,FIXMUL,CONPRD,
! $ FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7,
! $ WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7,
! $ OPRED1,OPRED2, OPRED4,OPRED5,OPRED6,OPRED7,
! $ MPRED3,CPRED4,TRCPRD,CO2MLT,SO2MLT,HNOMLT,N2OMLT )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER LBOT bottom layer number none
! LOGICAL LCO2 CO2 profile switch none
! LOGICAL LN2O N2O profile switch none
! LOGICAL LSO2 SO2 profile switch none
! LOGICAL LHNO3 HNO3 profile switch none
! LOGICAL LCO2PM CO2 ppmv profile switch none
! REAL ALAT profile latitude degrees (-90 to +90)
! REAL arr DZREF ref prof layer thickness meters
! REAL arr FX fixed gases adjustment none
! REAL arr PTEMP profile temperature K
! REAL arr PCAMNT prof carbon monoxide amnt kiloMoles/cm^2
! REAL arr PFAMNT profile CO2 gas amount kiloMoles/cm^2
! REAL arr PHAMNT profile HNO3 gas amount kiloMoles/cm^2
! REAL arr PMAMNT profile methane amount kiloMoles/cm^2
! REAL arr PNAMNT profile N2O amount kiloMoles/cm^2
! REAL arr POAMNT profile ozone amount kiloMoles/cm^2
! REAL arr PRES layer pressures atm
! REAL arr PSAMNT profile SO2 amount kiloMoles/cm^2
! REAL arr PWAMNT profile water amount kiloMoles/cm^2
! REAL arr RTEMP reference temperature K
! REAL arr RCAMNT ref carbon monoxide amount kiloMoles/cm^2
! REAL arr RFAMNT reference CO2 amount kiloMoles/cm^2
! REAL arr RHAMNT reference HNO3 amount kiloMoles/cm^2
! REAL arr RMAMNT reference methane amount kiloMoles/cm^2
! REAL arr RNAMNT reference N2O amount kiloMoles/cm^2
! REAL arr ROAMNT reference ozone amount kiloMoles/cm^2
! REAL arr RSAMNT reference SO2 amount kiloMoles/cm^2
! REAL arr RWAMNT reference water amount kiloMoles/cm^2
! REAL arr SECANG secant of path angle none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr CO2MLT CO2 multiplier none
! REAL arr CPRED4 carbon monoxide pred set4 various
! REAL arr FIXMUL fixed amount multiplier none
! REAL arr FPRED1 fixed predictors set1 various
! REAL arr FPRED2 fixed predictors set2 various
! REAL arr FPRED3 fixed predictors set3 various
! REAL arr FPRED4 fixed predictors set4 various
! REAL arr FPRED5 fixed predictors set5 various
! REAL arr FPRED6 fixed predictors set6 various
! REAL arr FPRED7 fixed predictors set7 various
! REAL arr HNOMLT HNO3 multiplier none
! REAL arr MPRED3 methane predictors set3 various
! REAL arr N2OMLT N2O multiplier none
! REAL arr OPRED1 ozone predictors set1 various
! REAL arr OPRED2 ozone predictors set2 various
! REAL arr OPRED4 ozone predictors set4 various
! REAL arr OPRED5 ozone predictors set5 various
! REAL arr OPRED6 ozone predictors set6 variou
! REAL arr OPRED7 ozone predictors set7 various
! REAL arr SO2MLT SO2 multiplier none
! REAL arr TRCPRD trace gas pert predictors various
! REAL arr WPRED1 water predictors set1 various
! REAL arr WPRED2 water predictors set2 various
! REAL arr WPRED3 water predictors set3 various
! REAL arr WPRED4 water predictors set4 various
! REAL arr WPRED5 water predictors set5 various
! REAL arr WPRED6 water predictors set6 various
! REAL arr WPRED7 water predictors set7 various
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! May 2008 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! Rapid transmittace algorithm predictors consisting of various gas
! amount and temperature ratios and offsets relative to a reference
! profile are calculated.
!
! ===================================================================
! The FTC profile variables computed for each layer are:
!
! ---------------------------------
! CONPRD: water continuum predictors (7 terms)
! 1) a*W/Tr^2 2) a*(W/Tr^2)^2 3) a*W/Tr 4) a*W^2/Tr
! 5) a*(W/Tr)^2 6) a*W/Tr^4 7) a*Wr
!
! -------------------------------
! Fixed predictors
!
! FPRED1: FWO (8 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) a*Trz/Tr
!
! FPRED2: FOW (8 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) a*Trz/Tr
!
! FPRED3: FMW (8 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) a*Trz/Tr
!
! FPRED4: FCOW (11 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) a^2*Trz
! 9) a^2*Tr 10) a^3 11) sqrt(a)
!
! FPRED5: FWO (11 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) a*Trz/Tr
! 9) a^2*Tr 10) sqrt(a) 11) Trz
!
! FPRED6: FWO (8 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) sqrt(a)
!
! FPRED7: FWO (8 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) sqrt(a)
!
! ---------------------------------
! Water predictors
!
! WPRED1: FWO (11 terms):
! 1) W*a 2) sqrt(W*a) 3) W*a*W/Wz
! 4) W*a*dT 5) (W*a)^2 6) sqrt(W*a)*dT
! 7) root^4(W*a) 8) sqrt(W*a)*W/Wz 9) (W*a)^3
! 10) W 11) W*a*dT*|dT|
!
! WPRED2: FOW (11 terms):
! 1) W*a 2) sqrt(W*a) 3) W*a*dT
! 4) W*a*Ox*a 5) (W*a)^2 6) root^4(W*a)
! 7) sqrt(W*a)*dT 8) W*a*W/Wz 9) (W*a)^3
! 10) W*a*(Ox*a)^2 11) sqrt(W*a)*W/Wz
!
! WPRED3: FMW (11 terms):
! 1) W*a 2) sqrt(W*a) 3) W*a*W/Wz
! 4) W*a*dT 5) (W*a)^2 6) sqrt(W*a)*dT
! 7) root^4(W*a) 8) (W*a)^3 9) W
! 10) sqrt(W*a)*W/Wz 11) sqrt(W*a)*Mz*a
!
! WPRED4: FCOW (13 terms):
! 1) W*a 2) W 3) sqrt(W*a)
! 4) W*a*dT 5) (W*a)^2 6) sqrt(W*a)*dT
! 7) root^4(W*a) 8) W*a*W/Wz 9) W*a^2
! 10) (W*a)^3 11) W*a*Cz*a 12) sqrt(W*a)*W/Wz
! 13) W*a^2*dT
!
! WPRED5: FWO bfsw (3 terms):
! 1) W*a 2) (W*a)^3/2 3) W*a*dT
!
! WPRED6: FWO mfmw (7 terms):
! 1) W*a 2) (W*a)^3/2 3) W*a*dT
! 4) (W*a)^2 5) (W*a)^3/2*dT 6) (W*a)^3
! 7) W*a^2
!
! WPRED7: FWO mfbw (13 terms):
! 1) W*a 2) (W*a)^3/2 3) W*a*dT
! 4) (W*a)^2 5) (W*a)^3/2*dT 6) (W*a)^3
! 7) W*a^2 8) W*a*W/Wz 9) (W*a)^3/2*W/Wz
! 10) (W*a)^5/4 11) (W*a)^2*W/Wz 12) W^2*a
! 13) (W*a)^7/4
!
! ---------------------------
! Ozone predictors
!
! OPRED1: FWO (5 terms):
! 1) O*a 2) sqrt(O*a) 3) O*a*dT
! 4) (O*a)^2 5) sqrt(O*a)*dT
!
! OPRED2: FOW (10 terms):
! 1) O*a 2) sqrt(O*a) 3) O*a*dT
! 4) (O*a)^2 5) sqrt(O*a)*dT 6) O*a*O/Ox
! 7) sqrt(O*a)*O/Ox 8) O*a*Oz/Ox 9) O*a*sqrt(Ox*a)
! 10) O*a*TOz*a
!
! OPRED4: FCOW (3 terms):
! 1) O*a 2) sqrt(O*a) 3) O*a*dT
!
! OPRED5: FWO bfsw (1 term):
! 1) O*a
!
! OPRED6: FWO mfmw (1 term):
! 1) O*a
!
! OPRED7: FWO mfbw (1 term):
! 1) O*a
!
! ---------------------------
! CPRED4: carbon monoxide predictors (11 terms):
! 1) C*a 2) sqrt(C*a) 3) C*a*dT
! 4) (C*a)^2 5) C*a*C/Cz 6) sqrt(C*a)*dT
! 7) root^4(C*a) 8) sqrt(C*a)*C/Cz 9) C
!
! ---------------------------
! MPRED3: methane predictors (9 terms):
! 1) M*a 2) sqrt(M*a) 3) M*a*dT
! 4) (M*a)^2 5) M*a^2 6) Mz*a
! 7) M*dT 8) TMz*a 9) sqrt(Mz*a)
!
! ---------------------------
! CO2PRD: CO2 perturbation coefs (4 terms):
! 1) a 2) Tr 3) a*Tr 4) a*Tr^2
!
! -----
! where:
! "a" is the secant of the viewing angle SECANG
! "Tr" is the temperature ratio PTEMP/RTEMP
! "Trz" is the pressure weighted temperature ratio above, i.e.
! the sum i=2 to i=L of { P(i) * ( P(i) - P(i-1) )* Tr(i-1) }
! where "P" is the pressure PRES and "L" is the layer number, and
! Trz(L=1)=0
! "W" is the water amount ratio PWAMNT/RWAMNT
! "dT" is the temperature offset PTEMP-RTEMP
! "Wz" is the pressure weighted water amount above ratio, the
! sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PWAMNT(i) },
! divided by the same sum except using RWAMNT instead of PWAMNT.
! For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
! "O" is the ozone amount ratio POAMNT/ROAMNT
! "Oz" is the pressure weighted ozone amount above ratio, the
! sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * POAMNT(i) },
! divided by the same sum except using ROAMNT instead of POAMNT.
! For these sums, term P(0) is defined as P(0)=2*P(1) - P(2)
! "Ox" is the unweighted ozone amount above ratio, the
! sum i=1 to i=L of { POAMNT(i) },
! divided by the same sum except using ROAMNT instead of POAMNT.
! For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
! "TOz" is the pressure and ozone weighted temperature ratio above,
! sum i=2 to i=L of { P(i) * ( P(i)-P(i-1) )* dT(i-1) * O(i-1) }
! and TOz(L=1)=0
! "C" is the carbon monoxide amount ratio POAMNT/ROAMNT
! "Cz" is the pressure weighted CO amount above ratio, the
! sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PCAMNT(i) },
! divided by the same sum except using RCAMNT instead of PCAMNT.
! For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
! "M" is the methane amount ratio PMAMNT/RMAMNT
! "Mz" is the pressure weighted methane amount above ratio, the
! sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PMAMNT(i) },
! divided by the same sum except using RMAMNT instead of PMAMNT.
! For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
! "TMz" is the pressure and methane weighted temperature ratio above,
! sum i=2 to i=L of { P(i) * ( P(i)-P(i-1) )* Tr(i-1) * M(i-1) }
! and TMz(L=1)=0
!
! -----------------------------------------------------
! FIXMUL: the not-quite-fixed "fixed" amount multiplier. The
! value should be close to (within a few percent of) unity.
! This term adjusts for the effects of water vapor displacement
! and latitude dependent gravity. The equations used below are
! a combination of analytic adjustments for water and gravity,
! as well as trial-and-error fudge factors to make it all work
! accurately for any realistic surface pressure and altitude.
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! Assumes the user has supplied vaguely realistic profile amounts
! and temperatures.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- --------------------------------------
! 1 Dec 1994 Scott Hannon Created
! 10 Apr 1995 Scott Hannon New header comments; redefined WZ;
! changed SECANG to array
! 6 Sep 1995 Scott Hannon Correct WZ for top layer
! 3 Feb 1997 Scott Hannon Re-wrote it for FWO+FOW+FMW+FCOW
! 7 Jul 1997 Scott Hannon Re-wrote it for sets 1 thru 7
! 30 Sep 1997 Scott Hannon Added CO2PRD
! 5 Mar 1998 Scott Hannon Deleted water pred 12 & 13 of set 1 & 3
! 26 Aug 1998 Scott Hannon Add LBOT to call; loop on LBOT instead
! of MAXLAY
! 31 Mar 2000 Scott Hannon Change FIXMUL equation; add ALAT input
! var; add FX and DZREF data; add
! PWATER, PMULT, and GSCAL local vars.
! 11 Aug 2000 Scott Hannon Change from 4 to 5 term H2O continuum
! 17 Aug 2000 Scott Hannon Add FX & DZREF input vars
! 12 Sep 2002 Scott Hannon Add predictors 6 & 7 to H2O con
! 28 Jun 2005 Scott Hannon "trace" version for CO2,SO2,HNO3,N2O
! 23 Jan 2008 Scott Hannon Add LCO2,LN2O,LSO2,LHNO3 switches for
! perturbation multiplier calcs; add
! LCO2PM to allow CO2 ppmv profile
! 14 May 2008 Scott Hannon Add no prof CO2MLT calc; add CO2TOP and
! CO2PPM to call; add CO2TOP calc
!END====================================================================
! =================================================================
       SUBROUTINE CALPAR ( LBOT, &
          RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RSAMNT,RHAMNT,RNAMNT, &
          PTEMP,PFAMNT,PWAMNT,POAMNT,PCAMNT,PMAMNT,PSAMNT,PHAMNT,PNAMNT, &
           PRES,SECANG, ALAT, FX, DZREF, &
           LCO2, LN2O, LSO2, LHNO3,LCO2PM,CO2PPM,CO2TOP, &
         FIXMUL,CONPRD, &
         FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7, &
         WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7, &
         OPRED1,OPRED2, OPRED4,OPRED5,OPRED6,OPRED7, &
         MPRED3,CPRED4,TRCPRD,CO2MLT,SO2MLT,HNOMLT,N2OMLT )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER LBOT
       REAL RTEMP(MAXLAY)
       REAL RFAMNT(MAXLAY)
       REAL RWAMNT(MAXLAY)
       REAL ROAMNT(MAXLAY)
       REAL RCAMNT(MAXLAY)
       REAL RMAMNT(MAXLAY)
       REAL RSAMNT(MAXLAY)
       REAL RHAMNT(MAXLAY)
       REAL RNAMNT(MAXLAY)
       REAL PTEMP(MAXLAY)
       REAL PFAMNT(MAXLAY)
       REAL PWAMNT(MAXLAY)
       REAL POAMNT(MAXLAY)
       REAL PCAMNT(MAXLAY)
       REAL PMAMNT(MAXLAY)
       REAL PSAMNT(MAXLAY)
       REAL PHAMNT(MAXLAY)
       REAL PNAMNT(MAXLAY)
       REAL PRES(MAXLAY)
       REAL SECANG(MAXLAY)
       REAL ALAT
       REAL FX(MAXLAY)
       REAL DZREF(MAXLAY)
       LOGICAL LCO2
       LOGICAL LN2O
       LOGICAL LSO2
       LOGICAL LHNO3
       LOGICAL LCO2PM
       REAL CO2PPM
!
! Output
       REAL CO2TOP
       REAL FIXMUL(MAXLAY)
       REAL CONPRD( N1CON,MAXLAY)
       REAL FPRED1( N1FIX,MAXLAY)
       REAL FPRED2( N2FIX,MAXLAY)
       REAL FPRED3( N3FIX,MAXLAY)
       REAL FPRED4( N4FIX,MAXLAY)
       REAL FPRED5( N5FIX,MAXLAY)
       REAL FPRED6( N6FIX,MAXLAY)
       REAL FPRED7( N7FIX,MAXLAY)
       REAL WPRED1( N1H2O,MAXLAY)
       REAL WPRED2( N2H2O,MAXLAY)
       REAL WPRED3( N3H2O,MAXLAY)
       REAL WPRED4( N4H2O,MAXLAY)
       REAL WPRED5( N5H2O,MAXLAY)
       REAL WPRED6( N6H2O,MAXLAY)
       REAL WPRED7( N7H2O,MAXLAY)
       REAL OPRED1( N1O3,MAXLAY)
       REAL OPRED2( N2O3,MAXLAY)
       REAL OPRED4( N4O3,MAXLAY)
       REAL OPRED5( N5O3,MAXLAY)
       REAL OPRED6( N6O3,MAXLAY)
       REAL OPRED7( N7O3,MAXLAY)
       REAL MPRED3( N3CH4,MAXLAY)
       REAL CPRED4( N4CO,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       REAL CO2MLT(MAXLAY)
       REAL SO2MLT(MAXLAY)
       REAL HNOMLT(MAXLAY)
       REAL N2OMLT(MAXLAY)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER L
       REAL PDP
       REAL PNORM
       REAL DT
       REAL TR
       REAL TZ
       REAL TRZ
       REAL A_F
       REAL A_W
       REAL WZREF
       REAL WZ
       REAL AZ_W
       REAL A_O
       REAL XZREF
       REAL XZ
       REAL XZ_O
       REAL OZREF
       REAL OZ
       REAL AZ_O
       REAL TOZ
       REAL TAZ_O
       REAL A_C
       REAL CZ
       REAL CZREF
       REAL AZ_C
       REAL A_M
       REAL MZREF
       REAL MZ
       REAL AZ_M
       REAL TMZ
       REAL TAZ_M
       REAL TJUNKS
       REAL WJUNKA
       REAL WJUNKR
       REAL WJUNKS
       REAL WJUNKZ
       REAL WJUNK4
       REAL OJUNKA
       REAL OJUNKR
       REAL OJUNKZ
       REAL OJUNKX
       REAL CJUNKA
       REAL CJUNKR
       REAL CJUNKS
       REAL CJUNKZ
       REAL MJUNKA
       REAL MJUNKR
       REAL MJUNKZ
! Variables for fixed gases adjustment
       REAL PWATER
       REAl GSCAL
! variables with DATA assignments
       REAL PMULT
       REAL STDDEN
       REAL STDTMP
       REAL KMOLE
! Data statments
       DATA PMULT /0.58/ ! fudge factor * (0.622=M_H2O/M_AIR)
       DATA STDDEN /2.6867E+19/ ! Loschmidt aka standard density
       DATA STDTMP /273.15/ ! Standard Temperature
       DATA KMOLE /6.022045E+26/ ! 1000 * Avagadro's Number
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Calc the fixed gases gravity correction factor
       GSCAL=( 9.78050518 + 0.0518017*( COS( (ALAT - 90.0)*PI/180) &
          )**2 )/9.80683613
!
! Initialize the sum terms to zero
       PNORM=0.0E+0
       TZ=0.0E+0
       WZREF=0.0E+0
       WZ=0.0E+0
       XZREF=0.0E+0
       XZ=0.0E+0
       OZREF=0.0E+0
       OZ=0.0E+0
       TOZ=0.0E+0
       CZREF=0.0E+0
       CZ=0.0E+0
       MZREF=0.0E+0
       MZ=0.0E+0
       TMZ=0.0E+0
       CO2TOP=0.0E+0
!
! --------------------
! Loop over the layers
! --------------------
       DO L=1,LBOT
!
! ---------------------------
! Calculate the basic profile
! dependent predictors.
! ---------------------------
!
          IF (L .EQ. 1) THEN
             PDP=PRES(1)*( PRES(2) - PRES(1))
             TRZ=0.0E+0
             TAZ_O=0.0E+0
             TAZ_M=0.0E+0
          ELSE
             PDP=PRES(L)*( PRES(L) - PRES(L-1) )
             PNORM=PNORM + PDP
!
! Note: TRZ, TOZ, and TMZ use layer-above terms
             TZ=TZ + PDP*TR
             TRZ=TZ/PNORM
!
             TOZ=TOZ + PDP*DT*A_O
             TAZ_O=TOZ/PNORM
!
             TMZ=TMZ + PDP*TR*A_M
             TAZ_M=TMZ/PNORM
          ENDIF
!
! Temperature terms
          DT=PTEMP(L) - RTEMP(L)
          TR=PTEMP(L)/RTEMP(L)
!
! Calc the fixed gases correction term for this layer
          PWATER=KMOLE*PWAMNT(L)*PTEMP(L)/(STDDEN*STDTMP*100*DZREF(L))
          A_F=( 1 - PMULT*PWATER/PRES(L) )/( FX(L)*GSCAL )
!cc
! for testing
! A_F=1.0
!cc
!
! Water terms
          A_W=PWAMNT(L)/RWAMNT(L)
          WZREF=WZREF + PDP*RWAMNT(L)
          WZ=WZ + PDP*PWAMNT(L)
          AZ_W=WZ/WZREF
!
! Ozone terms
          A_O=POAMNT(L)/ROAMNT(L)
          XZREF=XZREF + ROAMNT(L)
          XZ=XZ + POAMNT(L)
          XZ_O=XZ/XZREF
          OZREF=OZREF + PDP*ROAMNT(L)
          OZ=OZ + PDP*POAMNT(L)
          AZ_O=OZ/OZREF
!
! Carbon monoxide terms
          A_C=PCAMNT(L)/RCAMNT(L)
          CZREF=CZREF + PDP*RCAMNT(L)
          CZ=CZ + PDP*PCAMNT(L)
          AZ_C=CZ/CZREF
!
! Methane terms
          A_M=PMAMNT(L)/RMAMNT(L)
          MZREF=MZREF + PDP*RMAMNT(L)
          MZ=MZ + PDP*PMAMNT(L)
          AZ_M=MZ/MZREF
!
! ----------------------
! Load up the predictors
! ----------------------
!
! -----
! Fixed (for FWO, FOW, FMW, & FCOW)
! -----
          TJUNKS=TR*TR
          FIXMUL(L)=A_F
          FPRED1(1,L)=SECANG(L)
          FPRED1(2,L)=SECANG(L)*SECANG(L)
          FPRED1(3,L)=SECANG(L)*TR
          FPRED1(4,L)=SECANG(L)*TJUNKS
          FPRED1(5,L)=TR
          FPRED1(6,L)=TJUNKS
          FPRED1(7,L)=SECANG(L)*TRZ
          FPRED1(8,L)=SECANG(L)*TRZ/TR
!
          FPRED2(1,L)=FPRED1(1,L)
          FPRED2(2,L)=FPRED1(2,L)
          FPRED2(3,L)=FPRED1(3,L)
          FPRED2(4,L)=FPRED1(4,L)
          FPRED2(5,L)=FPRED1(5,L)
          FPRED2(6,L)=FPRED1(6,L)
          FPRED2(7,L)=FPRED1(7,L)
          FPRED2(8,L)=FPRED1(8,L)
!
          FPRED3(1,L)=FPRED1(1,L)
          FPRED3(2,L)=FPRED1(2,L)
          FPRED3(3,L)=FPRED1(3,L)
          FPRED3(4,L)=FPRED1(4,L)
          FPRED3(5,L)=FPRED1(5,L)
          FPRED3(6,L)=FPRED1(6,L)
          FPRED3(7,L)=FPRED1(7,L)
          FPRED3(8,L)=FPRED1(8,L)
!
          FPRED4(1,L)=SECANG(L)
          FPRED4(2,L)=SECANG(L)*SECANG(L)
          FPRED4(3,L)=SECANG(L)*TR
          FPRED4(4,L)=SECANG(L)*TJUNKS
          FPRED4(5,L)=TR
          FPRED4(6,L)=TJUNKS
          FPRED4(7,L)=SECANG(L)*TRZ
          FPRED4(8,L)=SECANG(L)*SECANG(L)*TRZ
          FPRED4(9,L)=SECANG(L)*SECANG(L)*TR
          FPRED4(10,L)=SECANG(L)*SECANG(L)*SECANG(L)
          FPRED4(11,L)=SQRT(SECANG(L))
!
! Fixed predictors for FWO sun bfsw = set5
          FPRED5(1,L)=SECANG(L)
          FPRED5(2,L)=SECANG(L)*SECANG(L)
          FPRED5(3,L)=SECANG(L)*TR
          FPRED5(4,L)=SECANG(L)*TJUNKS
          FPRED5(5,L)=TR
          FPRED5(6,L)=TJUNKS
          FPRED5(7,L)=SECANG(L)*TRZ
          FPRED5(8,L)=SECANG(L)*TRZ/TR
          FPRED5(9,L)=SECANG(L)*SECANG(L)*TR
          FPRED5(10,L)=SQRT(SECANG(L))
          FPRED5(11,L)=TRZ
!
! Fixed predictors for FWO sun mfmw = set6
          FPRED6(1,L)=SECANG(L)
          FPRED6(2,L)=SECANG(L)*SECANG(L)
          FPRED6(3,L)=SECANG(L)*TR
          FPRED6(4,L)=SECANG(L)*TJUNKS
          FPRED6(5,L)=TR
          FPRED6(6,L)=TJUNKS
          FPRED6(7,L)=SECANG(L)*TRZ
          FPRED6(8,L)=SQRT(SECANG(L))
!
! Fixed predictors for FWO sun mfbw = set7
          FPRED7(1,L)=SECANG(L)
          FPRED7(2,L)=SECANG(L)*SECANG(L)
          FPRED7(3,L)=SECANG(L)*TR
          FPRED7(4,L)=SECANG(L)*TJUNKS
          FPRED7(5,L)=TR
          FPRED7(6,L)=TJUNKS
          FPRED7(7,L)=SECANG(L)*TRZ
          FPRED7(8,L)=SQRT(SECANG(L))
!
!
! -----
! Ozone
! -----
          OJUNKA=SECANG(L)*A_O
          OJUNKR=SQRT( OJUNKA )
          OJUNKZ=OJUNKA/XZ_O
          OJUNKX=SECANG(L)*XZ_O
!
! Ozone predictors for FWO = set1
          OPRED1(1,L)=OJUNKA
          OPRED1(2,L)=OJUNKR
          OPRED1(3,L)=OJUNKA*DT
          OPRED1(4,L)=OJUNKA*OJUNKA
          OPRED1(5,L)=OJUNKR*DT
!
! ozone predictors for FOW = set2
          OPRED2( 1,L)=OJUNKA
          OPRED2( 2,L)=OJUNKR
          OPRED2( 3,L)=OJUNKA*DT
          OPRED2( 4,L)=OJUNKA*OJUNKA
          OPRED2( 5,L)=OJUNKR*DT
          OPRED2( 6,L)=OJUNKZ*A_O
          OPRED2( 7,L)=OJUNKR*A_O/XZ_O
          OPRED2( 8,L)=OJUNKZ*AZ_O
          OPRED2( 9,L)=OJUNKA*SQRT( OJUNKX )
          OPRED2(10,L)=OJUNKA*TAZ_O*SECANG(L)
!
! There are no ozone predictors for set3 = FMW (the ozone
! absorption in the region covered by FMW is negligible).
! ozone predictors for FCOW = set4
          OPRED4(1,L)=OJUNKA
          OPRED4(2,L)=OJUNKR
          OPRED4(3,L)=OJUNKA*DT
!
! ozone predictors for FWO sun bfsw = set5
          OPRED5(1,L)=OJUNKA
!
! ozone predictors for FWO sun mfmw = set6
          OPRED6(1,L)=OJUNKA
!
! ozone predictors for FWO sun mfbw = set7
          OPRED7(1,L)=OJUNKA
!
!
! -------
! Methane for FMW = set3
! -------
          MJUNKA=SECANG(L)*A_M
          MJUNKR=SQRT(MJUNKA)
          MJUNKZ=SECANG(L)*AZ_M
          MPRED3(1,L)=MJUNKA
          MPRED3(2,L)=MJUNKR
          MPRED3(3,L)=MJUNKA*DT
          MPRED3(4,L)=MJUNKA*MJUNKA
          MPRED3(5,L)=MJUNKA*SECANG(L)
          MPRED3(6,L)=MJUNKZ
          MPRED3(7,L)=A_M*DT
          MPRED3(8,L)=TAZ_M*SECANG(L)
          MPRED3(9,L)=SQRT( MJUNKZ )
!
!
! -----
! Water
! -----
          WJUNKA=SECANG(L)*A_W
          WJUNKR=SQRT( WJUNKA )
          WJUNKS=WJUNKA*WJUNKA
          WJUNKZ=WJUNKA*A_W/AZ_W
          WJUNK4=SQRT( WJUNKR )
!
! Water predictors for FWO = set1
          WPRED1( 1,L)=WJUNKA
          WPRED1( 2,L)=WJUNKR
          WPRED1( 3,L)=WJUNKZ
          WPRED1( 4,L)=WJUNKA*DT
          WPRED1( 5,L)=WJUNKS
          WPRED1( 6,L)=WJUNKR*DT
          WPRED1( 7,L)=WJUNK4
          WPRED1( 8,L)=WJUNKZ/WJUNKR
          WPRED1( 9,L)=WJUNKS*WJUNKA
          WPRED1(10,L)=A_W
          WPRED1(11,L)=WJUNKA*DT*ABS( DT )
!
! water predictors for FOW = set2
          WPRED2( 1,L)=WJUNKA
          WPRED2( 2,L)=WJUNKR
          WPRED2( 3,L)=WJUNKA*DT
          WPRED2( 4,L)=WJUNKA*OJUNKX
          WPRED2( 5,L)=WJUNKS
          WPRED2( 6,L)=WJUNK4
          WPRED2( 7,L)=WJUNKR*DT
          WPRED2( 8,L)=WJUNKZ
          WPRED2( 9,L)=WJUNKA*WJUNKS
          WPRED2(10,L)=WJUNKA*OJUNKX*OJUNKX
          WPRED2(11,L)=WJUNKZ/WJUNKR
!
! water predictors for FMW = set3
          WPRED3( 1,L)=WJUNKA
          WPRED3( 2,L)=WJUNKR
          WPRED3( 3,L)=WJUNKZ
          WPRED3( 4,L)=WJUNKA*DT
          WPRED3( 5,L)=WJUNKS
          WPRED3( 6,L)=WJUNKR*DT
          WPRED3( 7,L)=WJUNK4
          WPRED3( 8,L)=WJUNKS*WJUNKA
          WPRED3( 9,L)=A_W
          WPRED3(10,L)=WJUNKZ/WJUNKR
          WPRED3(11,L)=WJUNKR*MJUNKZ
!
! water predictors for FCOW = set4
          WPRED4( 1,L)=WJUNKA
          WPRED4( 2,L)=A_W
          WPRED4( 3,L)=WJUNKR
          WPRED4( 4,L)=WJUNKA*DT
          WPRED4( 5,L)=WJUNKS
          WPRED4( 6,L)=WJUNKR*DT
          WPRED4( 7,L)=WJUNK4
          WPRED4( 8,L)=WJUNKZ
          WPRED4( 9,L)=WJUNKA*SECANG(L)
          WPRED4(10,L)=WJUNKS*WJUNKA
          WPRED4(11,L)=WJUNKA*AZ_C*SECANG(L)
          WPRED4(12,L)=WJUNKZ/WJUNKR
          WPRED4(13,L)=WJUNKA*DT*SECANG(L)
!
! Water predictors for FWO sun bfsw = set5
          WPRED5( 1,L)=WJUNKA
          WPRED5( 2,L)=WJUNKA*WJUNKR
          WPRED5( 3,L)=WJUNKA*DT
!
! Water predictors for FWO sun mfmw = set6
          WPRED6( 1,L)=WJUNKA
          WPRED6( 2,L)=WJUNKA*WJUNKR
          WPRED6( 3,L)=WJUNKA*DT
          WPRED6( 4,L)=WJUNKS
          WPRED6( 5,L)=WJUNKA*WJUNKR*DT
          WPRED6( 6,L)=WJUNKA*WJUNKS
          WPRED6( 7,L)=WJUNKA*SECANG(L)
!
! Water predictors for FWO sun mfbw = set7
          WPRED7( 1,L)=WJUNKA
          WPRED7( 2,L)=WJUNKA*WJUNKR
          WPRED7( 3,L)=WJUNKA*DT
          WPRED7( 4,L)=WJUNKS
          WPRED7( 5,L)=WJUNKA*WJUNKR*DT
          WPRED7( 6,L)=WJUNKA*WJUNKS
          WPRED7( 7,L)=WJUNKA*SECANG(L)
          WPRED7( 8,L)=WJUNKZ
          WPRED7( 9,L)=WJUNKZ*WJUNKR
          WPRED7(10,L)=WJUNKA*WJUNK4
          WPRED7(11,L)=WJUNKA*WJUNKZ
          WPRED7(12,L)=WJUNKA*A_W
          WPRED7(13,L)=WJUNKS/WJUNK4
!
! ---------------
! Water continuum (for FWO, FOW, FMW, FCOW)
! ---------------
          CONPRD(1,L)=WJUNKA/TJUNKS
          CONPRD(2,L)=CONPRD(1,L)*A_W/TJUNKS
          CONPRD(3,L)=WJUNKA/TR
          CONPRD(4,L)=CONPRD(3,L)*A_W
          CONPRD(5,L)=CONPRD(1,L)*A_W
          CONPRD(6,L)=CONPRD(1,L)/TJUNKS
          CONPRD(7,L)=WJUNKA
!
!
! ---------------
! Carbon monoxide for FCOW = set4
! ---------------
          CJUNKA=SECANG(L)*A_C
          CJUNKR=SQRT( CJUNKA )
          CJUNKS=CJUNKA*CJUNKA
          CJUNKZ=CJUNKA*A_C/AZ_C
          CPRED4(1,L)=CJUNKA
          CPRED4(2,L)=CJUNKR
          CPRED4(3,L)=CJUNKA*DT
          CPRED4(4,L)=CJUNKS
          CPRED4(5,L)=CJUNKZ
          CPRED4(6,L)=CJUNKR*DT
          CPRED4(7,L)=SQRT( CJUNKR )
          CPRED4(8,L)=CJUNKZ/CJUNKR
          CPRED4(9,L)=A_C
          CPRED4(10,L)=CJUNKA*SECANG(L)
          CPRED4(11,L)=CJUNKR*SECANG(L)
!
! ---------------
! trace gas perturbation predictors
! ---------------
! The first 4 trace predictors are used by all trace gases
          TRCPRD(1,L)=SECANG(L)
          TRCPRD(2,L)=TR
          TRCPRD(3,L)=SECANG(L)*TR
          TRCPRD(4,L)=SECANG(L)*TJUNKS
! The last 3 trace predictors are only used by N2O
          TRCPRD(5,L)=SECANG(L)*SECANG(L)
          TRCPRD(6,L)=1.0
          TRCPRD(7,L)=SQRT( SECANG(L) )
!
          IF (LCO2) THEN
             IF (LCO2PM) THEN
                CO2MLT(L)=100.0*(PFAMNT(L) - CO2STD)/(3.0*CO2STD)
             ELSE
! CO2 mult=1 when prof amount = 1.03 * ref amount
                CO2MLT(L)=33.3333*( PFAMNT(L) - FIXMUL(L)*RFAMNT(L) )/ &
                   RFAMNT(L)
! Ignore changes in CO2 of less than ~0.03%
                IF (ABS(CO2MLT(L)) .LT. 1E-2) CO2MLT(L)=0.0
             ENDIF
          ELSE
             CO2MLT(L)=100.0*(CO2PPM - CO2STD)/(3.0*CO2STD)
          ENDIF
          IF (L .LE. NTEBOT) THEN
             CO2TOP=CO2TOP + CO2STD*(1.0 + CO2MLT(L)*3.0E-2)
          ENDIF
!
          IF (LN2O) THEN
! N2O mult=-1 when prof amount = 0.75 * ref amount
             N2OMLT(L)=4.0*( PNAMNT(L) - FIXMUL(L)*RNAMNT(L) )/ &
                RNAMNT(L)
! Ignore changes in N2O less than ~0.3%
             IF (ABS(N2OMLT(L)) .LT. 1E-2) N2OMLT(L)=0.0
          ELSE
             N2OMLT(L)=0.0
          ENDIF
!
          IF (LSO2) THEN
! SO2 mult=1 when prof amount = 1000 * ref amount
             SO2MLT(L)=1.0010E-3*( PSAMNT(L) - FIXMUL(L)*RSAMNT(L) )/ &
                RSAMNT(L)
! Ignore changes in SO2 of less than ~10%
             IF (ABS(SO2MLT(L)) .LT. 1E-4) SO2MLT(L)=0.0
          ELSE
             SO2MLT(L)=0.0
          ENDIF
!
          IF (LHNO3) THEN
! HNO3 mult=1 when prof amount = 2 * ref amount
             HNOMLT(L)=( PHAMNT(L) - FIXMUL(L)*RHAMNT(L) )/ &
                RHAMNT(L)
! Ignore changes in HNO3 less than ~1%
             IF (ABS(HNOMLT(L)) .LT. 1E-2) HNOMLT(L)=0.0
          ELSE
             HNOMLT(L)=0.0
          ENDIF
!
!cc this block for testing
! N2OMLT(L)=0.0
! SO2MLT(L)=0.0
! HNOMLT(L)=0.0
!cc
!
       ENDDO
! End loop over layers
!
! Convert CO2TOP from sum to mean
       CO2TOP=CO2TOP/AMIN0(NTEBOT, LBOT)
!
!
       RETURN
       END SUBROUTINE CALPAR
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALRAD
!
!F77====================================================================
!ROUTINE NAME:
! CALRAD
!ABSTRACT:
! Calculate a profile's radiance.
!CALL PROTOCOL:
! CALRAD ( NCHAN, FREQ, TAU, TP, TBOT, EBOT, LBOT,
! SUNCOS, RHOSUN, DISTES, HSUN, TAUZSN,
! SEC, RHOTHR, LABOVE, COEFF, TAUZ, RAD, BT)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr COEFF thermal F factor coefs various
! REAL DISTES Earth-sun distance meters
! REAL arr EBOT bottom surface emissivity none
! REAL arr FREQ channel frequencies cm^-1
! INT arr LABOVE layer-above for thermal none
! INTEGER LBOT bottom layer none
! INTEGER NCHAN number of channels none
! REAL arr RHOSUN reflectivity for solar 1/steradian
! REAL arr RHOTHR reflectivity for thermal 1/steradian
! REAL SEC bottlom path angle secant none
! REAL SUNCOS sun angle cosine none
! REAL arr TAU effective layer trans none
! REAL arr TAUZ layer-to-space trans none
! REAL arr TAUZSN eff sun angle l-to-s trans none
! REAL TBOT bottom surface temperature Kelvin
! REAL arr TP temperature profile Kelvin
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr BT brightness temperature Kelvin
! REAL arr RAD radiance W/(m^2.str.cm^-1)
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! March 1998 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The radiance is calculated for each channels in turn. The rad
! is a sum of four components: surface, upwelling (non-reflected)
! atmospheric, reflected downwelling atmospheric thermal, and
! reflected solar. No scattering.
!
! Comment: this routine could easily be re-written to use layer-to-
! space transmittances rather than layer transittances. Currently
! the CALT# routines compute layer transmittances (since it's a bit
! faster and more accurate when using the QIKEXP function.
!
! ===================================================================
! Computes black body emissions for each layer using the Planck
! equation:
! planck = c1*v^3/( exp(c2*v/T) - 1 )
! where c1 and c2 are the radiation constants, T is the temperature
! TP, and v is the frequency FREQ.
!
! We assume the layers emit radiances of
! rad_layer = (1 - tau)*planck
! where tau is the layer transmittance TAU.
!
! The total radiance leaving the bottom surface and going upward
! is the surface emission and reflected solar & thermal. The
! reflected solar and thermal are handled as seperate terms added
! to radiance arriving at the satellite.
! rad_surface = e*planck
! where e is the bottom surface emissivity EBOT, and the surface is
! at temperature TBOT.
!
! We trace the upward radiance thru the atmosphere and determine
! the total radiance leaving the top layer (and then reaching the
! satellite) is:
! the sum L=L_bot downto 1 of { rad(L-1)*tau(L) + rad(L) }
! where rad(L_bot-1) = rad_surface, and rad(1) = RAD.
!
! The reflected solar term is based on an approximation suggested
! by J.Susskind et al. The reflected solar radiance reaching the
! satellite is given by
! Rsun = rho * omega * TAUZSN * Hsun
! where omega is the solid angle of the sun as seen from Earth,
! Hsun is the (non-reflected) solar radiance at the top of the
! atmosphere, and TAUZSN is (surface) layer-to-space transmittance
! of a path along an effective total angle defined as
! secant_eff = secant_view + secant_sun
! Note that this requires a seperate transmittance calculation
! at the effective sun angle. Hsun is passed to this routine (it
! is close to planck for 5800 K), while omega is computed using
! the distance of the Earth from the sun DISTES
! omega = pi * ( radius_sun / distance_Earth_sun )^2
!
! The reflected downwelling thermal is another approximation based
! on a method suggested by Susskind et al. It uses the viewing
! angle layer-to-space transmittance TAUZ, the radiance of a
! single layer somewhere above the surface, and a parameterized "F"
! factor (which is sort of a fudge factor determined by regression).
! Rtherm = rho * pi * F * planck * TAUZ*(1-TAUZ)
! where the planck radiance is computed for layer L = LBOT - LABOVE.
!
! For convenience we also output brightness temperature, which is
! related to the radiance by inverting the planck equation:
! BT = c2*v/ln( 1 + c1*v^3/RAD )
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! Currently this routine does not handle scattering or clouds.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 2 Sep 1997 Scott Hannon Created from an extensive re-write of
! our Feb97 CALRAD routine for Mar98 FTC.
! 07 Apr 2005 Scott Hannon Change refl therm calc
! 09 May 2005 Scott Hannon Add default F=1 to refl therm so that
! RTHERM is always calculated
! 13 Dec 2005 Scott Hannon Add limits check to "F"
!END====================================================================
! =================================================================
       SUBROUTINE CALRAD ( NCHAN, FREQ, TAU, TP, TBOT, EBOT, LBOT, &
          SUNCOS, RHOSUN, DISTES, HSUN, TAUZSN, &
          SEC, RHOTHR, LABOVE, COEFF, TAUZ, RAD_OUT, BT_OUT )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
 REAL,INTENT(OUT) :: RAD_OUT(MXCHAN) !Radiance [Watts/m^2/cm^-1]
 REAL,INTENT(OUT) :: BT_OUT(MXCHAN)  !Brightness temperature [K]

       INTEGER NCHAN
       REAL FREQ(MXCHAN)
       REAL TAU(MAXLAY,MXCHAN)
       REAL TP(MAXLAY)
       REAL TBOT
       REAL EBOT(MXCHAN)  !bottom emissivity
       INTEGER LBOT
       REAL SUNCOS
       REAL RHOSUN(MXCHAN)
       REAL DISTES
       REAL HSUN(MXCHAN)
       REAL TAUZSN(MXCHAN)
       REAL SEC
       REAL RHOTHR(MXCHAN)
       INTEGER LABOVE(MXCHAN)
       REAL COEFF(NFCOEF,MXCHAN)
       REAL TAUZ(MXCHAN)
       REAL RAD(MXCHAN)
       REAL BT(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER L
       INTEGER LTHERM
       REAL RPLNCK(MAXLAY)
       REAL C1V3
       REAL C2V
       REAL SUNFAC
       REAL RSUN
       REAL F
       REAL RTHERM
! Downwelling atmospheric thermal emission terms
       REAL TDOWNN ! "near-side" layer-to-surface trans
       REAL TDOWNF ! "far-side" layer-to-surface trans
       REAL RDOWN ! downward radiance
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Note: on average, DISTES = 1.496E+11 m. The exact value varies
! with time since the Earth's orbit is slightly elliptical.
       SUNFAC=SUNCOS*PI*(RADSUN/DISTES)**2
! Note: PI*(RADSUN/DISTES)^2 = omega = solid angle [steradians] of
! the sun as seen from Earth. The above equation is actually an
! approximation for the case DISTES >> RADSUN. The exact equation
! is omega = 2*pi*(1 - DISTES/sqrt(DISTES^2 + RADSUN^2)).
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------


  LOOP_CHAN: DO I=1,NCHAN
!
! Calc c1*v^3 and c2*v
          C1V3=C1*(FREQ(I)**3)
          C2V=C2*FREQ(I)

!
! ---------------------------
! Calc the up going radiance
! (excluding reflected terms)
! ---------------------------
! Initialize the upward radiance with
! the bottom surface emission
          RAD(I)=EBOT(I)*C1V3/( EXP( C2V/TBOT ) - 1.0 )

!
! Initialize downwelling thermal terms
          RDOWN=0.0
          TDOWNN=1.0
!
! Loop upward over the layers
          DO L=LBOT,1,-1
! Calculate the Planck function for this layer
             RPLNCK(L)=C1V3/( EXP( C2V/TP(L) ) - 1.0 )
!
! Calc the upward radiance thru and from this layer
             RAD(I)=( RAD(I)*TAU(L,I) ) + &
                ( RPLNCK(L)*(1.0E+0 - TAU(L,I)) )
! Calc the downward radiance from this layer
             TDOWNF=TDOWNN*TAU(L,I)
             RDOWN = RDOWN + ( RPLNCK(L)*(TDOWNN - TDOWNF) )
             TDOWNN=TDOWNF

!print*,I,L,RPLNCK(L),RAD(I),TDOWNF,RDOWN,TDOWNN

!cc
! if (I .EQ. 1928) then
! print*, RAD(I)
! endif
!cc
          ENDDO
!
! --------------------------
! Calc the reflected solar
! rad reaching the satellite
! --------------------------
          RSUN=RHOSUN(I)*SUNFAC*HSUN(I)*TAUZSN(I)
!
! ----------------------------------
! Calc the reflected downwelling
! thermal rad reaching the satellite
! ----------------------------------
          F=1.0
          IF (TAUZ(I) .GT. 0.0005) THEN
             F= COEFF(1,I) + &
                ( COEFF(2,I)/SEC ) + &
                ( COEFF(3,I)*TAUZ(I) ) + &
                ( COEFF(4,I)*TAUZ(I)*TAUZ(I) ) + &
                ( COEFF(5,I)*TAUZ(I)/SEC ) + &
                ( COEFF(6,I)*TAUZ(I)/RDOWN )
! Truncate F at limits as needed
             F = MAX( MIN(F,2.09), 0.696 )
          ENDIF
          RTHERM=RHOTHR(I)*PI*RDOWN*F*TAUZ(I)
!
! --------------------------------------------------
! Add on the reflected solar and downwelling thermal
! --------------------------------------------------
!cc
! for testing
! RTHERM=0.0
! RSUN=0.0
!cc

   RAD_OUT(I)=RAD(I) + RSUN + RTHERM

!
! ------------------------------------------
! Convert radiance to brightness temperature
! ------------------------------------------
!cc removed 07 Apr 2005 since never used
   BT_OUT(I)=C2V/LOG( 1.0 + C1V3/RAD(I) )

!  print*,'in CALRAD', I, RAD_OUT(I), BT_OUT(I)
!cc
!

   ENDDO LOOP_CHAN

! End loops on channel number (frequency)
!

  RETURN
  END SUBROUTINE CALRAD
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALRAD
!
!F77====================================================================
!ROUTINE NAME:
! CALRAD
!ABSTRACT:
! Calculate a profile's radiance.
!CALL PROTOCOL:
! CALRAD ( NCHAN, FREQ, TAU, TP, TBOT, EBOT, LBOT,
! SUNCOS, RHOSUN, DISTES, HSUN, TAUZSN,
! SEC, RHOTHR, LABOVE, COEFF, TAUZ, RAD, BT)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr COEFF thermal F factor coefs various
! REAL DISTES Earth-sun distance meters
! REAL arr EBOT bottom surface emissivity none
! REAL arr FREQ channel frequencies cm^-1
! INT arr LABOVE layer-above for thermal none
! INTEGER LBOT bottom layer none
! INTEGER NCHAN number of channels none
! REAL arr RHOSUN reflectivity for solar 1/steradian
! REAL arr RHOTHR reflectivity for thermal 1/steradian
! REAL SEC bottlom path angle secant none
! REAL SUNCOS sun angle cosine none
! REAL arr TAU effective layer trans none
! REAL arr TAUZ layer-to-space trans none
! REAL arr TAUZSN eff sun angle l-to-s trans none
! REAL TBOT bottom surface temperature Kelvin
! REAL arr TP temperature profile Kelvin
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr BT brightness temperature Kelvin
! REAL arr RAD radiance W/(m^2.str.cm^-1)
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! March 1998 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The radiance is calculated for each channels in turn. The rad
! is a sum of four components: surface, upwelling (non-reflected)
! atmospheric, reflected downwelling atmospheric thermal, and
! reflected solar. No scattering.
!
! Comment: this routine could easily be re-written to use layer-to-
! space transmittances rather than layer transittances. Currently
! the CALT# routines compute layer transmittances (since it's a bit
! faster and more accurate when using the QIKEXP function.
!
! ===================================================================
! Computes black body emissions for each layer using the Planck
! equation:
! planck = c1*v^3/( exp(c2*v/T) - 1 )
! where c1 and c2 are the radiation constants, T is the temperature
! TP, and v is the frequency FREQ.
!
! We assume the layers emit radiances of
! rad_layer = (1 - tau)*planck
! where tau is the layer transmittance TAU.
!
! The total radiance leaving the bottom surface and going upward
! is the surface emission and reflected solar & thermal. The
! reflected solar and thermal are handled as seperate terms added
! to radiance arriving at the satellite.
! rad_surface = e*planck
! where e is the bottom surface emissivity EBOT, and the surface is
! at temperature TBOT.
!
! We trace the upward radiance thru the atmosphere and determine
! the total radiance leaving the top layer (and then reaching the
! satellite) is:
! the sum L=L_bot downto 1 of { rad(L-1)*tau(L) + rad(L) }
! where rad(L_bot-1) = rad_surface, and rad(1) = RAD.
!
! The reflected solar term is based on an approximation suggested
! by J.Susskind et al. The reflected solar radiance reaching the
! satellite is given by
! Rsun = rho * omega * TAUZSN * Hsun
! where omega is the solid angle of the sun as seen from Earth,
! Hsun is the (non-reflected) solar radiance at the top of the
! atmosphere, and TAUZSN is (surface) layer-to-space transmittance
! of a path along an effective total angle defined as
! secant_eff = secant_view + secant_sun
! Note that this requires a seperate transmittance calculation
! at the effective sun angle. Hsun is passed to this routine (it
! is close to planck for 5800 K), while omega is computed using
! the distance of the Earth from the sun DISTES
! omega = pi * ( radius_sun / distance_Earth_sun )^2
!
! The reflected downwelling thermal is another approximation based
! on a method suggested by Susskind et al. It uses the viewing
! angle layer-to-space transmittance TAUZ, the radiance of a
! single layer somewhere above the surface, and a parameterized "F"
! factor (which is sort of a fudge factor determined by regression).
! Rtherm = rho * pi * F * planck * TAUZ*(1-TAUZ)
! where the planck radiance is computed for layer L = LBOT - LABOVE.
!
! For convenience we also output brightness temperature, which is
! related to the radiance by inverting the planck equation:
! BT = c2*v/ln( 1 + c1*v^3/RAD )
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! Currently this routine does not handle scattering or clouds.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 2 Sep 1997 Scott Hannon Created from an extensive re-write of
! our Feb97 CALRAD routine for Mar98 FTC.
! 07 Apr 2005 Scott Hannon Change refl therm calc
! 09 May 2005 Scott Hannon Add default F=1 to refl therm so that
! RTHERM is always calculated
! 13 Dec 2005 Scott Hannon Add limits check to "F"
! 03 Jul 2007 S.Hannon/L.Chen - PGE variant with added check of RSUN
! and different limits on F check
!END====================================================================
! =================================================================
       SUBROUTINE CALRAD_PGE ( NCHAN, FREQ, TAU, TP, TBOT, EBOT, LBOT, &
          SUNCOS, RHOSUN, DISTES, HSUN, TAUZSN, &
          SEC, RHOTHR, LABOVE, COEFF, TAUZ, RAD, BT)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
       INTEGER NCHAN
       REAL FREQ(MXCHAN)
       REAL TAU(MAXLAY,MXCHAN)
       REAL TP(MAXLAY)
       REAL TBOT
       REAL EBOT(MXCHAN)
       INTEGER LBOT
       REAL SUNCOS
       REAL RHOSUN(MXCHAN)
       REAL DISTES
       REAL HSUN(MXCHAN)
       REAL TAUZSN(MXCHAN)
       REAL SEC
       REAL RHOTHR(MXCHAN)
       INTEGER LABOVE(MXCHAN)
       REAL COEFF(NFCOEF,MXCHAN)
       REAL TAUZ(MXCHAN)
       REAL RAD(MXCHAN)
       REAL BT(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER L
       INTEGER LTHERM
       REAL RPLNCK(MAXLAY)
       REAL C1V3
       REAL C2V
       REAL SUNFAC
       REAL RSUN
       REAL F
       REAL RTHERM
! Downwelling atmospheric thermal emission terms
       REAL TDOWNN ! "near-side" layer-to-surface trans
       REAL TDOWNF ! "far-side" layer-to-surface trans
       REAL RDOWN ! downward radiance
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Note: on average, DISTES = 1.496E+11 m. The exact value varies
! with time since the Earth's orbit is slightly elliptical.
       SUNFAC=SUNCOS*PI*(RADSUN/DISTES)**2
! Note: PI*(RADSUN/DISTES)^2 = omega = solid angle [steradians] of
! the sun as seen from Earth. The above equation is actually an
! approximation for the case DISTES >> RADSUN. The exact equation
! is omega = 2*pi*(1 - DISTES/sqrt(DISTES^2 + RADSUN^2)).
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO 210 I=1,NCHAN
!
! Calc c1*v^3 and c2*v
          C1V3=C1*(FREQ(I)**3)
          C2V=C2*FREQ(I)
!
! ---------------------------
! Calc the up going radiance
! (excluding reflected terms)
! ---------------------------
! Initialize the upward radiance with
! the bottom surface emission
          RAD(I)=EBOT(I)*C1V3/( EXP( C2V/TBOT ) - 1.0 )
!
! Initialize downwelling thermal terms
          RDOWN=0.0
          TDOWNN=1.0
!
! Loop upward over the layers
          DO L=LBOT,1,-1
! Calculate the Planck function for this layer
             RPLNCK(L)=C1V3/( EXP( C2V/TP(L) ) - 1.0 )
!
! Calc the upward radiance thru and from this layer
             RAD(I)=( RAD(I)*TAU(L,I) ) + &
                ( RPLNCK(L)*(1.0E+0 - TAU(L,I)) )
! Calc the downward radiance from this layer
             TDOWNF=TDOWNN*TAU(L,I)
             RDOWN = RDOWN + ( RPLNCK(L)*(TDOWNN - TDOWNF) )
             TDOWNN=TDOWNF
!cc
! if (I .EQ. 1928) then
! print*, RAD(I)
! endif
!cc
          ENDDO
!
! --------------------------
! Calc the reflected solar
! rad reaching the satellite
! --------------------------
          RSUN=RHOSUN(I)*SUNFAC*HSUN(I)*TAUZSN(I)
!cc PGE variant
          IF (RSUN .LT. 0) RSUN=0.0
!cc
!
! ----------------------------------
! Calc the reflected downwelling
! thermal rad reaching the satellite
! ----------------------------------
          F=1.0
          IF (TAUZ(I) .GT. 0.0005) THEN
             F= COEFF(1,I) + &
                ( COEFF(2,I)/SEC ) + &
                ( COEFF(3,I)*TAUZ(I) ) + &
                ( COEFF(4,I)*TAUZ(I)*TAUZ(I) ) + &
                ( COEFF(5,I)*TAUZ(I)/SEC ) + &
                ( COEFF(6,I)*TAUZ(I)/RDOWN )
! Truncate F at limits as needed
!cc PGE variant
! F = MAX( MIN(F,2.09), 0.696 )
             F = MAX( MIN(F,2.5), 0.4 )
!cc
          ENDIF
          RTHERM=RHOTHR(I)*PI*RDOWN*F*TAUZ(I)
!
! --------------------------------------------------
! Add on the reflected solar and downwelling thermal
! --------------------------------------------------
!cc
! for testing
! RTHERM=0.0
! RSUN=0.0
!cc
          RAD(I)=RAD(I) + RSUN + RTHERM
!
! ------------------------------------------
! Convert radiance to brightness temperature
! ------------------------------------------
!cc removed 07 Apr 2005 since never used
! BT(I)=C2V/LOG( 1.0 + C1V3/RAD(I) )
!cc
!
 210 CONTINUE
! End loops on channel number (frequency)
!
       RETURN
       END SUBROUTINE CALRAD_PGE
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALT1 (for set1 = FWO) version with trace gases
!
!F77====================================================================
!ROUTINE NAME:
! CALT1
!ABSTRACT:
! Calculate the transmittance for set1 using the predictors
! and the fast transmittance coefficients.
!CALL PROTOCOL:
! CALT1 ( INDCHN, NLAY, BLMULT, NCHN1, CLIST1, COEF1,
! FIXMUL, CONPD1, FPRED1, WPRED1, OPRED1, TRCPRD,
! INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
! INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT,
! INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX,
! LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL BLMULT bottom layer opt depth mult none
! INT arr CLIST1 set channel list none
! REAL arr COEF1 set1 fast trans coefs various
! REAL arr CONPD1 set1 H2O continuum preds various
! REAL arr FIXMUL fixed amount mult (~1.0) none
! REAL arr FPRED1 set1 fixed gases preds various
! INT arr INDCHN channel indices none
! INTEGER NLAY Number of layers to bottom none
! INTEGER NCHN1 set1 number of channels none
! REAL arr OPRED1 set1 ozone predictors various
! REAL arr WPRED1 set1 water predictors various
! REAL arr TRCPRD trace gas pert predictors various
! INT arr INDCO2 CO2 pert chan indices none
! REAL arr COFCO2 CO2 pert coefs various
! REAL arr CO2MLT CO2 pert multiplier none
! INT arr INDSO2 SO2 pert chan indices none
! REAL arr COFSO2 SO2 pert coefs various
! REAL SO2MLT SO2 pert multiplier none
! INT arr INDHNO HNO3 pert chan indices none
! REAL arr COFHNO HNO3 pert coefs various
! REAL arr HNOMLT HNO3 pert multiplier none
! INT arr INDN2O N2O pert chan indices none
! REAL arr COFN2O N2O pert coefs various
! REAL arr N2OMLT N2O pert multiplier none
! INT arr INDH2O OPTRAN H2O chan indices none
! REAL arr H2OPRD OPTRAN H2O predictors various
! REAL arr COFH2O OPTRAN H2O coefs various
! INTEGER LOPMAX OPTRAN max level none
! INTEGER LOPLOW OPTRAN low bracketing level none
! LOG arr LOPUSE OPTRAN level needed? none
! REAL arr WAOP OPTRAN layer water amounts kilomoles/cm^2
! REAL arr DAOP OPTRAN-to-AIRS interp fac none
! REAL arr WAANG AIRS layer water amounts kilomoles/cm^2
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr TAU effective layer trans none
! REAL arr TAUZ layer-to-space trans none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2000 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The FTC coefficents and profile FTC variables are multiplied
! together and summed to calculate the effective layer
! transmittances. Fixed, water, and ozone transmittances are each
! checked individually to be sure they give 0 < trans < 1.
!
! ===================================================================
! The routine loops over the selected channels of set1. For each
! channel, it first decides if needs to do a calculation for
! variable CO2, and also if it needs to do an OPTRAN water calc (if
! so, it does so immediately). The program then loops downward over
! all the layers and computes the layer transmittances TAU.
!
! The water continuum absorption coefficient is:
! k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
!
! The layer effective fixed gas absorption coefficient is:
! k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
! The layer effective water lines absorption coefficient is:
! k_water = the sum i=1 to 11 of { COEF(5+8+i)*WPRED(i) }
!
! The layer effective ozone absorption coefficient is:
! k_ozone = the sum i=1 to 5 of { COEF(5+8+11+i)*OPRED(i) }
!
! where
! "COEF" are the fast transmittance coefficients COEF1
! "CONPRD" are the water continuum predictors CONPRD
! "FPRED" are the fixed gases predictors FPRED1
! "WPRED" are the water lines predictors WPRED1
! "OPRED" are the ozone predictors OPRED1
!
! The total layer effective transmittance TAU is:
! TAU = exp( -[ k_con + k_fixed + k_water + k_ozone ])
!
! To help speed up the exponential calculations, we use our own
! "EXP" replacement function called QIKEXP which uses just the
! first few series expansion terms for exp(x) if x is suitably small.
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! Dec 1 1994 Scott Hannon Created
! 3 Feb 1997 Scott Hannon Re-wrote (from CALTAU) for FWO
! 3 Sep 1997 Scott Hannon Added TAUZ and BLMULT
! 30 Sep 1997 Scott Hannon Added variable CO2
! 27 Feb 1998 Scott Hannon Added OPTRAN H2O
! 6 Mar 1998 Scott Hannon Deleted water preds 12 & 13 and shifted
! ozone coef indices to 24-28 (was 26-30)
! 4 May 1998 Scott Hannon Fix error: INDH2O(MXCHAN) not (MXCHNW)
! 26 Aug 1998 Scott Hannon Add NLAY to call to CALOKW
! 11 Aug 2000 Scott Hannon Change from 4 to 5 term H2O continuum
! 12 Sep 2002 Scott Hannon Add predictors 6 & 7 to H2O con
! 18 May 2005 Scott Hannon Add HNO3 based on SO2 code
! 28 Jun 2005 Scott Hannon "trace" version for CO2,SO2,HNO3,N2O
! 13 Sep 2010 Scott Hannon Add 5th CO2 coef
!END====================================================================
! =================================================================
       SUBROUTINE CALT1 ( INDCHN, NLAY, BLMULT, NCHN1, CLIST1, COEF1, &
          FIXMUL, CONPD1, FPRED1, WPRED1, OPRED1, TRCPRD, &
          INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
          INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, &
          INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE, &
            WAOP, DAOP, WAANG, TAU, TAUZ)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! QIKEXP
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER INDCHN(MXCHAN)
       INTEGER NLAY
       REAL BLMULT
       INTEGER NCHN1
       INTEGER CLIST1(MXCHN1)
       REAL COEF1(N1COEF,MAXLAY,MXCHN1)
       REAL FIXMUL(MAXLAY)
       REAL CONPD1( N1CON,MAXLAY)
       REAL FPRED1( N1FIX,MAXLAY)
       REAL WPRED1( N1H2O,MAXLAY)
       REAL OPRED1( N1O3,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2( NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT(MAXLAY)
       INTEGER INDSO2(MXCHAN)
       REAL COFSO2( NSO2,MAXLAY,MXCHNS)
       REAL SO2MLT(MAXLAY)
       INTEGER INDHNO(MXCHAN)
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH)
       REAL HNOMLT(MAXLAY)
      INTEGER INDN2O(MXCHAN)
       REAL COFN2O( NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
       INTEGER INDH2O(MXCHAN)
       REAL H2OPRD( NH2O,MXOWLY)
       REAL COFH2O( NH2O,MXOWLY,MXCHNW)
       INTEGER LOPMIN
       INTEGER LOPMAX
       INTEGER LOPLOW(MAXLAY)
       LOGICAL LOPUSE(MXOWLY)
       REAL WAOP(MXOWLY)
       REAL DAOP(MAXLAY)
       REAL WAANG(MAXLAY)
!
! Output
       REAL TAU(MAXLAY,MXCHAN)
       REAL TAUZ(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER ICO2
       INTEGER IHNO3
       INTEGER ILAY
       INTEGER IN2O
       INTEGER ISO2
       INTEGER J
       REAL DK
       REAL DKCO2
       REAL DKHNO3
       REAL DKN2O
       REAL DKSO2
       REAL KCON
       REAL KFIX
       REAL KLAYER
       REAL KOZO
       REAL KZ
       REAL KZFW
       LOGICAL LCO2
       LOGICAL LH2O
       LOGICAL LHNO3
       LOGICAL LN2O
       LOGICAL LSO2
!
! for function QIKEXP
!       REAL QIKEXP
!
! for CALOKW
       INTEGER IH2O
       REAL KW(MAXLAY)
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO I=1,NCHN1
!
! Array index of channel in TAU
          J=INDCHN( CLIST1(I) )
!
! Determine whether or not to do variable CO2 calc
          ICO2=INDCO2( CLIST1(I) )
          IF (ICO2 .GT. 0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable SO2 calc
          ISO2=INDSO2( CLIST1(I) )
          IF (ISO2 .GT. 0) THEN
             LSO2=.TRUE.
          ELSE
             LSO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable HNO3 calc
          IHNO3=INDHNO( CLIST1(I) )
          IF (IHNO3 .GT. 0) THEN
             LHNO3=.TRUE.
          ELSE
             LHNO3=.FALSE.
          ENDIF
!
! Determine whether or not to do variable N2O calc
          IN2O=INDN2O( CLIST1(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
!
! -------------------------
! Do OPTRAN water if needed
! -------------------------
          IH2O=INDH2O( CLIST1(I) )
          IF (IH2O .GT. 0) THEN
             LH2O=.FALSE.
! Calc OPTRAN water
!
             CALL CALOKW( NLAY, IH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE, &
                H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW )
!
          ELSE
             LH2O=.TRUE.
          ENDIF
!
! Initialize the layer-to-space optical depth
          KZ=0.0E+0
          KZFW=0.0E+0
!
! ------------------------------
! Loop on layers (top to bottom)
! ------------------------------
          DO ILAY=1,NLAY
! ---------------------------
! Compute the water continuum
! ---------------------------
             KCON=( COEF1(1,ILAY,I)*CONPD1(1,ILAY) ) + &
                  ( COEF1(2,ILAY,I)*CONPD1(2,ILAY) ) + &
                  ( COEF1(3,ILAY,I)*CONPD1(3,ILAY) ) + &
                  ( COEF1(4,ILAY,I)*CONPD1(4,ILAY) ) + &
                  ( COEF1(5,ILAY,I)*CONPD1(5,ILAY) ) + &
                  ( COEF1(6,ILAY,I)*CONPD1(6,ILAY) ) + &
                  ( COEF1(7,ILAY,I)*CONPD1(7,ILAY) )
!
             IF (KCON .LT. 0.0E+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 1.0E+1) THEN
                KCON=1.0E+1
             ENDIF
!
! -----------------------------
! Calc the fixed gases abs coef
! -----------------------------
             KFIX=( COEF1( 8,ILAY,I)*FPRED1(1,ILAY) ) + &
                  ( COEF1( 9,ILAY,I)*FPRED1(2,ILAY) ) + &
                  ( COEF1(10,ILAY,I)*FPRED1(3,ILAY) ) + &
                  ( COEF1(11,ILAY,I)*FPRED1(4,ILAY) ) + &
                  ( COEF1(12,ILAY,I)*FPRED1(5,ILAY) ) + &
                  ( COEF1(13,ILAY,I)*FPRED1(6,ILAY) ) + &
                  ( COEF1(14,ILAY,I)*FPRED1(7,ILAY) ) + &
                  ( COEF1(15,ILAY,I)*FPRED1(8,ILAY) )
!
             KFIX=KFIX*FIXMUL(ILAY)
!
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
!
! --------------------------
! Compute the water abs coef
! --------------------------
             IF (LH2O) THEN
! Not an OPTRAN water channel
                KW(ILAY)= &
                     ( COEF1(16,ILAY,I)*WPRED1( 1,ILAY) ) + &
                     ( COEF1(17,ILAY,I)*WPRED1( 2,ILAY) ) + &
                     ( COEF1(18,ILAY,I)*WPRED1( 3,ILAY) ) + &
                     ( COEF1(19,ILAY,I)*WPRED1( 4,ILAY) ) + &
                     ( COEF1(20,ILAY,I)*WPRED1( 5,ILAY) ) + &
                     ( COEF1(21,ILAY,I)*WPRED1( 6,ILAY) ) + &
                     ( COEF1(22,ILAY,I)*WPRED1( 7,ILAY) ) + &
                     ( COEF1(23,ILAY,I)*WPRED1( 8,ILAY) ) + &
                     ( COEF1(24,ILAY,I)*WPRED1( 9,ILAY) ) + &
                     ( COEF1(25,ILAY,I)*WPRED1(10,ILAY) ) + &
                     ( COEF1(26,ILAY,I)*WPRED1(11,ILAY) )
!
                IF (KW(ILAY) .LT. 0.0E+0) KW(ILAY)=0.0E+0
             ENDIF
!
! --------------------------
! Compute the ozone abs coef
! --------------------------
             KOZO=( COEF1(27,ILAY,I)*OPRED1(1,ILAY) ) + &
                  ( COEF1(28,ILAY,I)*OPRED1(2,ILAY) ) + &
                  ( COEF1(29,ILAY,I)*OPRED1(3,ILAY) ) + &
                  ( COEF1(30,ILAY,I)*OPRED1(4,ILAY) ) + &
                  ( COEF1(31,ILAY,I)*OPRED1(5,ILAY) )
!
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+1) THEN
                KOZO=1.0E+1
             ENDIF
!
! Update KZFW
             KZFW=KZFW + KFIX + KW(ILAY)
!
! ----------------------------------
! Calc the total layer transmittance
! ----------------------------------
!
!cccc
! This block is usually commented out and is only uncommented for
! testing purposes.
!c
! kcon=0.0E+0
! kfix=0.0E+0
! kw(ilay)=0.0E+0
! kozo=0.0E+0
!cccc
!
! ----------------------------
! Calc change in total optical
! depth due to variable CO2
! ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) + &
                      ( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) + &
                      ( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) + &
                      ( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) + &
                      ( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable SO2
! ----------------------------
             IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0) THEN
                DKSO2=( COFSO2(1,ILAY,ISO2)*TRCPRD(1,ILAY) ) + &
                      ( COFSO2(2,ILAY,ISO2)*TRCPRD(2,ILAY) ) + &
                      ( COFSO2(3,ILAY,ISO2)*TRCPRD(3,ILAY) ) + &
                      ( COFSO2(4,ILAY,ISO2)*TRCPRD(4,ILAY) )
                DKSO2=DKSO2*SO2MLT(ILAY)
             ELSE
                DKSO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable HNO3
! ----------------------------
             IF (LHNO3 .AND. HNOMLT(ILAY) .NE. 0) THEN
                DKHNO3=( COFHNO(1,ILAY,IHNO3)*TRCPRD(1,ILAY) ) + &
                       ( COFHNO(2,ILAY,IHNO3)*TRCPRD(2,ILAY) ) + &
                       ( COFHNO(3,ILAY,IHNO3)*TRCPRD(3,ILAY) ) + &
                       ( COFHNO(4,ILAY,IHNO3)*TRCPRD(4,ILAY) )
                DKHNO3=DKHNO3*HNOMLT(ILAY)
             ELSE
                DKHNO3=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable N2O
! ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) + &
                      ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) + &
                      ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) + &
                      ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) + &
                      ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) + &
                      ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) + &
                      ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
!
! ------------------------------------------
! Calc total optical depth and transmittance
! ------------------------------------------
! Calc total layer optical depth
!cc
! this block for testing
! DKHNO3=0.0
! DKSO2=0.0
! DKCO2=0.0
! DKN2O=0.0
!cc
! Limit -DK so it can never totally totally cancel KFIX
             DK = DKCO2 + DKSO2 + DKHNO3 + DKN2O
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF
             KLAYER = KCON + KFIX + KW(ILAY) + KOZO + DK
!
! Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
!
! Calc layer-to-space optical depth
             KZ=KZ + KLAYER
!
! Calc effective layer transmittance
             TAU(ILAY,J)=QIKEXP(-KLAYER)
!
          ENDDO
! End loop on levels
!
! Calc surface-to-space transmittance
          TAUZ(J)=QIKEXP(-KZ)
!
       ENDDO
! End loops on channel number (frequency)
!
       RETURN
       END SUBROUTINE CALT1
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALT2 (for set2 = FOW) version with trace gases
!
!F77====================================================================
!ROUTINE NAME:
! CALT2
!ABSTRACT:
! Calculate the transmittance for set2 using the predictors and the
! fast transmittance coefficients.
!CALL PROTOCOL:
! CALT2 ( INDCHN, NLAY, BLMULT, NCHN2, CLIST2, COEF2, FIXMUL,
! CONPD2, FPRED2, OPRED2, WPRED2, TRCPRD,
! INDCO2, COFCO2, CO2MLT, INDSO2, SOFCO2, SO2MLT,
! INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, TAU, TAUZ )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INT arr INDCHN channel indices none
! INTEGER NLAY number of layers to bottom none
! REAL BLMULT bottom layer opt depth mult none
! INTEGER NCHN2 set2 number of channels none
! INT arr CLIST2 set2 channel list none
! REAL arr COEF2 set2 fast trans coefs various
! REAL arr FIXMUL fixed amount mult (~1.0) none
! REAL arr CONPD2 set2 H2O continuum preds various
! REAL arr FPRED2 set2 fixed gases preds various
! REAL arr OPRED2 set2 ozone predictors various
! REAL arr WPRED2 set2 water predictors various
! REAL arr TRCPRD Trace gas pert predictors various
! INT arr INDCO2 CO2 pert chan indices none
! REAL arr COFCO2 CO2 pert coefs various
! REAL arr CO2MLT CO2 pert multiplier none
! INT arr INDSO2 SO2 pert chan indices none
! REAL arr COFSO2 SO2 pert coefs various
! REAL arr SO2MLT SO2 pert multiplier none
! INT arr INDHNO HNO3 pert chan indices none
! REAL arr COFHNO HNO3 pert coefs various
! REAL arr HNOMLT HNO3 pert multiplier none
! INT arr INDN2O N2O pert chan indices none
! REAL arr COFN2O N2O pert coefs various
! REAL arr N2OMLT N2O pert multiplier none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr TAU effective layer trans none
! REAL arr TAUZ layer-to-space trans none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2000 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The fast trans coefficents and predictors are multiplied
! together and summed to calculate the effective layer
! transmittances. Fixed, ozone, and water transmittances are each
! checked individually to be sure they give 0 < trans < 1.
!
! ===================================================================
! Loops downward over all the layers for each of the NCHN2 channels
! to compute the layer transmittances TAU.
!
! The water continuum absorption coefficient is:
! k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
!
! The layer effective fixed gas absorption coefficient is:
! k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
! The layer effective ozone absorption coefficient is:
! k_ozone = the sum i=1 to 10 of { COEF(5+8+i)*OPRED(i) }
!
! The layer effective water lines absorption coefficient is:
! k_water = the sum i=1 to 11 of { COEF(5+8+10+i)*WPRED(i) }
!
! where
! "COEF" are the fast transmittance coefficients COEF2
! "CONPRD" are the water continuum predictors CONPRD
! "FPRED" are the fixed gases predictors FPRED2
! "OPRED" are the ozone predictors OPRED2
! "WPRED" are the water lines predictors WPRED2
!
! The total layer effective transmittance TAU is:
! TAU = exp( -[ k_con + k_fixed + k_ozone + k_water])
!
! To help speed up the exponential calculations, we use our own
! "EXP" replacement function called QIKEXP which uses just the
! first few series expansion terms for exp(x) if x is suitably small.
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! Dec 1 1994 Scott Hannon Created
! 3 Feb 1997 Scott Hannon Re-wrote (from CALTAU) for FOW
! 3 Sep 1997 Scott Hannon Added TAUZ and BLMULT
! 30 Sep 1997 Scott Hannon Added variable CO2
! 11 Aug 2000 Scott Hannon Change from 4 to 5 term H2O continuum
! 12 Sep 2002 Scott Hannon Add predictors 6 & 7 to H2O con
! 18 May 2005 Scott Hannon Add HNO3 based on SO2 code
! 28 Jun 2005 Scott Hannon "trace" version for CO2,SO2,HNO3,N2O.
! 13 Sep 2010 Scott Hannon Add 5th CO2 coef
!END====================================================================
! =================================================================
       SUBROUTINE CALT2 ( INDCHN, NLAY, BLMULT, NCHN2, CLIST2, COEF2, &
          FIXMUL, CONPD2, FPRED2, OPRED2, WPRED2, TRCPRD, &
          INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
          INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, TAU, TAUZ )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! QIKEXP
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER INDCHN(MXCHAN)
       INTEGER NLAY
       REAL BLMULT
       INTEGER NCHN2
       INTEGER CLIST2(MXCHN2)
       REAL COEF2(N2COEF,MAXLAY,MXCHN2)
       REAL FIXMUL(MAXLAY)
       REAL CONPD2( N2CON,MAXLAY)
       REAL FPRED2( N2FIX,MAXLAY)
       REAL OPRED2( N2O3,MAXLAY)
       REAL WPRED2( N2H2O,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2( NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT(MAXLAY)
       INTEGER INDSO2(MXCHAN)
       REAL COFSO2( NSO2,MAXLAY,MXCHNS)
       REAL SO2MLT(MAXLAY)
       INTEGER INDHNO(MXCHAN)
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH)
       REAL HNOMLT(MAXLAY)
       INTEGER INDN2O(MXCHAN)
       REAL COFN2O( NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
!
! Output
       REAL TAU(MAXLAY,MXCHAN)
       REAL TAUZ(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER ICO2
       INTEGER IHNO3
       INTEGER ILAY
       INTEGER IN2O
       INTEGER ISO2
       INTEGER J
       REAL DK
       REAL DKCO2
       REAL DKHNO3
       REAL DKN2O
       REAL DKSO2
       REAL KCON
       REAL KFIX
       REAL KLAYER
       REAL KOZO
       REAL KWAT
       REAL KZ
       LOGICAL LCO2
       LOGICAL LHNO3
       LOGICAL LN2O
       LOGICAL LSO2
!
! for function QIKEXP
!       REAL QIKEXP
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO I=1,NCHN2
!
! Index for TAU
          J=INDCHN( CLIST2(I) )
!
! Determine whether or not to do variable CO2
          ICO2=INDCO2( CLIST2(I) )
          IF (ICO2 .GT. 0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable CO2
          ISO2=INDSO2( CLIST2(I) )
          IF (ISO2 .GT. 0) THEN
             LSO2=.TRUE.
          ELSE
             LSO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable HNO3
          IHNO3=INDHNO( CLIST2(I) )
          IF (IHNO3 .GT. 0) THEN
             LHNO3=.TRUE.
          ELSE
             LHNO3=.FALSE.
          ENDIF
!
! Determine whether or not to do variable N2O
          IN2O=INDN2O( CLIST2(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
!
! Initialize the layer-to-space optical depth
          KZ=0.0E+0
!
! ------------------------------
! Loop on layers (top to ground)
! ------------------------------
          DO ILAY=1,NLAY
! ---------------------------
! Compute the water continuum
! ---------------------------
             KCON=( COEF2(1,ILAY,I)*CONPD2(1,ILAY) ) + &
                  ( COEF2(2,ILAY,I)*CONPD2(2,ILAY) ) + &
                  ( COEF2(3,ILAY,I)*CONPD2(3,ILAY) ) + &
                  ( COEF2(4,ILAY,I)*CONPD2(4,ILAY) ) + &
                  ( COEF2(5,ILAY,I)*CONPD2(5,ILAY) ) + &
                  ( COEF2(6,ILAY,I)*CONPD2(6,ILAY) ) + &
                  ( COEF2(7,ILAY,I)*CONPD2(7,ILAY) )
!
             IF (KCON .LT. 0.0+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 1.0E+1) THEN
                KCON=1.0E+1
             ENDIF
!
! -----------------------------
! Calc the fixed gases abs coef
! -----------------------------
             KFIX=( COEF2( 8,ILAY,I)*FPRED2(1,ILAY) ) + &
                  ( COEF2( 9,ILAY,I)*FPRED2(2,ILAY) ) + &
                  ( COEF2(10,ILAY,I)*FPRED2(3,ILAY) ) + &
                  ( COEF2(11,ILAY,I)*FPRED2(4,ILAY) ) + &
                  ( COEF2(12,ILAY,I)*FPRED2(5,ILAY) ) + &
                  ( COEF2(13,ILAY,I)*FPRED2(6,ILAY) ) + &
                  ( COEF2(14,ILAY,I)*FPRED2(7,ILAY) ) + &
                  ( COEF2(15,ILAY,I)*FPRED2(8,ILAY) )
!
             KFIX=KFIX*FIXMUL(ILAY)
!
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
!
! --------------------------
! Compute the ozone abs coef
! --------------------------
             KOZO=( COEF2(16,ILAY,I)*OPRED2( 1,ILAY) ) + &
                  ( COEF2(17,ILAY,I)*OPRED2( 2,ILAY) ) + &
                  ( COEF2(18,ILAY,I)*OPRED2( 3,ILAY) ) + &
                  ( COEF2(19,ILAY,I)*OPRED2( 4,ILAY) ) + &
                  ( COEF2(20,ILAY,I)*OPRED2( 5,ILAY) ) + &
                  ( COEF2(21,ILAY,I)*OPRED2( 6,ILAY) ) + &
                  ( COEF2(22,ILAY,I)*OPRED2( 7,ILAY) ) + &
                  ( COEF2(23,ILAY,I)*OPRED2( 8,ILAY) ) + &
                  ( COEF2(24,ILAY,I)*OPRED2( 9,ILAY) ) + &
                  ( COEF2(25,ILAY,I)*OPRED2(10,ILAY) )
!
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+1) THEN
                KOZO=1.0E+1
             ENDIF
!
! --------------------------
! Compute the water abs coef
! --------------------------
             KWAT=( COEF2(26,ILAY,I)*WPRED2( 1,ILAY) ) + &
                  ( COEF2(27,ILAY,I)*WPRED2( 2,ILAY) ) + &
                  ( COEF2(28,ILAY,I)*WPRED2( 3,ILAY) ) + &
                  ( COEF2(29,ILAY,I)*WPRED2( 4,ILAY) ) + &
                  ( COEF2(30,ILAY,I)*WPRED2( 5,ILAY) ) + &
                  ( COEF2(31,ILAY,I)*WPRED2( 6,ILAY) ) + &
                  ( COEF2(32,ILAY,I)*WPRED2( 7,ILAY) ) + &
                  ( COEF2(33,ILAY,I)*WPRED2( 8,ILAY) ) + &
                  ( COEF2(34,ILAY,I)*WPRED2( 9,ILAY) ) + &
                  ( COEF2(35,ILAY,I)*WPRED2(10,ILAY) ) + &
                  ( COEF2(36,ILAY,I)*WPRED2(11,ILAY) )
!
             IF (KWAT .LT. 0.0E+0) THEN
                KWAT=0.0E+0
             ELSEIF( KWAT .GT. 1.0E+1) THEN
                KWAT=1.0E+1
             ENDIF
!
! ----------------------------------
! Calc the total layer transmittance
! ----------------------------------
!
!cccc
! This block is usually commented out and is only uncommented for
! testing purposes.
!
! kcon=0.0E+0
! kfix=0.0E+0
! kozo=0.0E+0
! kwat=0.0E+0
!cccc
!
! ----------------------------
! Calc change in total optical
! depth due to variable CO2
! ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) + &
                      ( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) + &
                      ( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) + &
                      ( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) + &
                      ( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable SO2
! ----------------------------
             IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0) THEN
                DKSO2=( COFSO2(1,ILAY,ISO2)*TRCPRD(1,ILAY) ) + &
                      ( COFSO2(2,ILAY,ISO2)*TRCPRD(2,ILAY) ) + &
                      ( COFSO2(3,ILAY,ISO2)*TRCPRD(3,ILAY) ) + &
                      ( COFSO2(4,ILAY,ISO2)*TRCPRD(4,ILAY) )
                DKSO2=DKSO2*SO2MLT(ILAY)
             ELSE
                DKSO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable HNO3
! ----------------------------
             IF (LHNO3 .AND. HNOMLT(ILAY) .NE. 0) THEN
                DKHNO3=( COFHNO(1,ILAY,IHNO3)*TRCPRD(1,ILAY) ) + &
                       ( COFHNO(2,ILAY,IHNO3)*TRCPRD(2,ILAY) ) + &
                       ( COFHNO(3,ILAY,IHNO3)*TRCPRD(3,ILAY) ) + &
                       ( COFHNO(4,ILAY,IHNO3)*TRCPRD(4,ILAY) )
                DKHNO3=DKHNO3*HNOMLT(ILAY)
             ELSE
                DKHNO3=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable N2O
! ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) + &
                      ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) + &
                      ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) + &
                      ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) + &
                      ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) + &
                      ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) + &
                      ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
!
! ------------------------------------------
! Calc total optical depth and transmittance
! ------------------------------------------
! Calc total layer optical depth
!cc
! this block for testing
! DKCO2=0.0
! DKSO2=0.0
! DKHNO3=0.0
! DKN2O=0.0
!cc
! Limit -DK so it can never totally totally cancel KFIX
             DK = DKCO2 + DKSO2 + DKHNO3 + DKN2O
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF
             KLAYER = KCON + KFIX + KOZO + KWAT + DK
!
! Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
!
! Calc layer-to-space optical depth
             KZ=KZ + KLAYER
!
! Calc effective layer transmittance
             TAU(ILAY,J)=QIKEXP(-KLAYER)
!
          ENDDO
! End loop on levels
!
! Convert KZ to TAUZ
          TAUZ(J)=QIKEXP(-KZ)
!
       ENDDO
! End loops on channel number (frequency)
!
       RETURN
       END SUBROUTINE CALT2
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALT3 (for set3 = FMW) version with trace gases (no CO2)
!
!F77====================================================================
!ROUTINE NAME:
! CALT3
!ABSTRACT:
! Calculate the transmittance for set3 using the prdictor and the
! fast transmittance coefficients.
!CALL PROTOCOL:
! CALT3 ( INDCHN, NLAY, BLMULT, NCHN3, CLIST3, COEF3,
! FIXMUL, CONPD3, FPRED3, MPRED3, WPRED3, TRCPRD,
! INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT,
! INDN2O, COFN2O, N2OMLT, INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX,
! LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INT arr INDCHN channel indices none
! INTEGER NLAY number of layers to bottom none
! REAL BLMULT bottom layer opt depth mult none
! INTEGER NCHN3 set3 number of channels none
! INT arr CLIST3 set3 channel list none
! REAL arr COEF3 set3 fast trans coefs various
! REAL arr FIXMUL fixed amount mult (~1.0) none
! REAL arr CONPD3 set3 H2O continuum preds various
! REAL arr FPRED3 set3 fixed gases preds various
! REAL arr MPRED3 set3 methane predictors various
! REAL arr WPRED3 set3 water predictors various
! REAL arr TRCPRD trace gas pert predictors various
! INT arr INDSO2 SO2 pert chan indices none
! REAL arr COFSO2 SO2 pert coefs various
! REAL arr SO2MLT SO2 pert multiplier none
! INT arr INDHNO HNO3 pert chan indices none
! REAL arr COFHNO HNO3 pert coefs various
! REAL arr HNOMLT HNO3 pert multiplier none
! INT arr INDN2O N2O pert chan indices none
! REAL arr COFN2O N2O pert coefs various
! REAL arr N2OMLT N2O pert multiplier none
! INT arr INDH2O OPTRAN H2O chan indices none
! REAL arr H2OPRD OPTRAN H2O predictors various
! REAL arr COFH2O OPTRAN H2O coefs various
! INTEGER LOPMAX OPTRAN max level none
! INTEGER LOPLOW OPTRAN low bracketing level none
! LOG arr LOPUSE OPTRAN level needed? none
! REAL arr WAOP OPTRAN layer water amounts kilomoles/cm^2
! REAL arr DAOP OPTRAN-to-AIRS interp fac none
! REAL arr WAANG AIRS layer water amounts kilomoles/cm^2
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr TAU effective layer trans none
! REAL arr TAUZ layer-to-space trans none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2000 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The fast trans coefficents and predictors are multiplied
! together and summed to calculate the effective layer
! transmittances. Fixed, methane, and water transmittances are each
! checked individually to be sure they give 0 < trans < 1.
!
! ===================================================================
! Loops downward over all the layers for each of the NCHN3 channels
! to compute the layer transmittances TAU.
!
! The water continuum absorption coefficient is:
! k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
!
! The layer effective fixed gas absorption coefficient is:
! k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
! The layer effective methane absorption coefficient is:
! k_methane = the sum i=1 to 9 of { COEF(5+8+i)*OPRED(i) }
!
! The layer effective water lines absorption coefficient is:
! k_water = the sum i=1 to 11 of { COEF(5+8+9+i)*WPRED(i) }
!
! where
! "COEF" are the fast transmittance coefficients COEF3
! "CONPRD" are the water continuum predictors CONPRD
! "FPRED" are the fixed gases predictors FPRED3
! "MPRED" are the methane predictors OPRED3
! "WPRED" are the water lines predictors WPRED3
!
! The total layer effective transmittance is:
! TAU = exp( -[ k_con + k_fixed + k_methane + k_water])
!
! To help speed up the exponential calculations, we use our own
! "EXP" replacement function called QIKEXP which uses just the
! first few series expansion terms for exp(x) if x is suitably small.
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! Dec 1 1994 Scott Hannon Created
! 3 Feb 1997 Scott Hannon Re-wrote (from CALTAU) for FMW
! 3 Sep 1997 Scott Hannon Added TAUZ and BLMULT
! 5 Mar 1998 Scott Hannon Added OPTRAN water and deleted water
! preds 12 & 13
! 4 May 1998 Scott Hannon Fix error: INDH2O(MXCHAN) not (MXCHNW)
! 26 Aug 1998 Scott Hannon Fix mistake: loop on NLAY not MAXLAY;
! Add NLAY to call to CALOKW
! 11 Aug 2000 Scott Hannon Change from 4 to 5 term H2O continuum
! 12 Sep 2002 Scott Hannon Add predictors 6 & 7 to H2O con
! 25 Apr 2003 Scott Hannon Add HNO3 based on SO2 code
! 28 Jun 2005 Scott Hannon "trace" version for SO2,HNO3,N2O
!END====================================================================
! =================================================================
       SUBROUTINE CALT3 ( INDCHN, NLAY, BLMULT, NCHN3, CLIST3, COEF3, &
          FIXMUL, CONPD3, FPRED3, MPRED3, WPRED3, TRCPRD, &
          INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT, &
          INDN2O, COFN2O, N2OMLT, INDH2O, H2OPRD, COFH2O, &
          LOPMIN, LOPMAX, LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! QIKEXP
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER INDCHN(MXCHAN)
       INTEGER NLAY
       REAL BLMULT
       INTEGER NCHN3
       INTEGER CLIST3(MXCHN3)
       REAL COEF3(N3COEF,MAXLAY,MXCHN3)
       REAL FIXMUL(MAXLAY)
       REAL CONPD3( N3CON,MAXLAY)
       REAL FPRED3( N3FIX,MAXLAY)
       REAL MPRED3( N3CH4,MAXLAY)
       REAL WPRED3( N3H2O,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDSO2(MXCHAN)
       REAL COFSO2( NSO2,MAXLAY,MXCHNS)
       REAL SO2MLT(MAXLAY)
       INTEGER INDHNO(MXCHAN)
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH)
       REAL HNOMLT(MAXLAY)
       INTEGER INDN2O(MXCHAN)
       REAL COFN2O( NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
       INTEGER INDH2O(MXCHAN)
       REAL H2OPRD( NH2O,MXOWLY)
       REAL COFH2O( NH2O,MXOWLY,MXCHNW)
       INTEGER LOPMIN
       INTEGER LOPMAX
       INTEGER LOPLOW(MAXLAY)
       LOGICAL LOPUSE(MXOWLY)
       REAL WAOP(MXOWLY)
       REAL DAOP(MAXLAY)
       REAL WAANG(MAXLAY)
!
! Output
       REAL TAU(MAXLAY,MXCHAN)
       REAL TAUZ(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER IHNO3
       INTEGER IN2O
       INTEGER ILAY
       INTEGER ISO2
       INTEGER J
       REAL DK
       REAL DKHNO3
       REAL DKN2O
       REAL DKSO2
       REAL KCON
       REAL KFIX
       REAL KMET
       REAL KZ
       REAL KZFMW
       REAL KLAYER
       LOGICAL LH2O
       LOGICAL LHNO3
       LOGICAL LN2O
       LOGICAL LSO2
!
! for function QIKEXP
!       REAL QIKEXP
!
! for CALOKW
       INTEGER IH2O
       REAL KW(MAXLAY)
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO I=1,NCHN3
!
! Index for TAU
          J=INDCHN( CLIST3(I) )
!
!
! Determine whether or not to do variable SO2
          ISO2=INDSO2( CLIST3(I) )
          IF (ISO2 .GT. 0) THEN
             LSO2=.TRUE.
          ELSE
             LSO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable HNO3
          IHNO3=INDHNO( CLIST3(I) )
          IF (IHNO3 .GT. 0) THEN
             LHNO3=.TRUE.
          ELSE
             LHNO3=.FALSE.
          ENDIF
!
! Determine whether or not to do variable N2O
          IN2O=INDN2O( CLIST3(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
!
! -------------------------
! Do OPTRAN water if needed
! -------------------------
          IH2O=INDH2O( CLIST3(I) )
          IF (IH2O .GT. 0) THEN
             LH2O=.FALSE.
! Calc OPTRAN water
!
             CALL CALOKW( NLAY, IH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE, &
                H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW )
!
          ELSE
             LH2O=.TRUE.
          ENDIF
!
! Initialize the layer-to-space optical depth
          KZ=0.0E+0
          KZFMW=0.0E+0
!
! ------------------------------
! Loop on layers (top to ground)
! ------------------------------
          DO ILAY=1,NLAY
!
! ---------------------------
! Compute the water continuum
! ---------------------------
             KCON=( COEF3(1,ILAY,I)*CONPD3(1,ILAY) ) + &
                  ( COEF3(2,ILAY,I)*CONPD3(2,ILAY) ) + &
                  ( COEF3(3,ILAY,I)*CONPD3(3,ILAY) ) + &
                  ( COEF3(4,ILAY,I)*CONPD3(4,ILAY) ) + &
                  ( COEF3(5,ILAY,I)*CONPD3(5,ILAY) ) + &
                  ( COEF3(6,ILAY,I)*CONPD3(6,ILAY) ) + &
                  ( COEF3(7,ILAY,I)*CONPD3(7,ILAY) )
!
             IF (KCON .LT. 0.0E+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 1.0E+1) THEN
                KCON=1.0E+1
             ENDIF
!
! -----------------------------
! Calc the fixed gases abs coef
! -----------------------------
             KFIX=( COEF3( 8,ILAY,I)*FPRED3(1,ILAY) ) + &
                  ( COEF3( 9,ILAY,I)*FPRED3(2,ILAY) ) + &
                  ( COEF3(10,ILAY,I)*FPRED3(3,ILAY) ) + &
                  ( COEF3(11,ILAY,I)*FPRED3(4,ILAY) ) + &
                  ( COEF3(12,ILAY,I)*FPRED3(5,ILAY) ) + &
                  ( COEF3(13,ILAY,I)*FPRED3(6,ILAY) ) + &
                  ( COEF3(14,ILAY,I)*FPRED3(7,ILAY) ) + &
                  ( COEF3(15,ILAY,I)*FPRED3(8,ILAY) )
!
             KFIX=KFIX*FIXMUL(ILAY)
!
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
!
! ----------------------------
! Compute the methane abs coef
! ----------------------------
             KMET=( COEF3(16,ILAY,I)*MPRED3(1,ILAY) ) + &
                  ( COEF3(17,ILAY,I)*MPRED3(2,ILAY) ) + &
                  ( COEF3(18,ILAY,I)*MPRED3(3,ILAY) ) + &
                  ( COEF3(19,ILAY,I)*MPRED3(4,ILAY) ) + &
                  ( COEF3(20,ILAY,I)*MPRED3(5,ILAY) ) + &
                  ( COEF3(21,ILAY,I)*MPRED3(6,ILAY) ) + &
                  ( COEF3(22,ILAY,I)*MPRED3(7,ILAY) ) + &
                  ( COEF3(23,ILAY,I)*MPRED3(8,ILAY) ) + &
                  ( COEF3(24,ILAY,I)*MPRED3(9,ILAY) )
!
             IF (KMET .LT. 0.0E+0) THEN
                KMET=0.0E+0
             ELSEIF (KMET .GT. 1.0E+1) THEN
                KMET=1.0E+1
             ENDIF
!
! --------------------------
! Compute the water abs coef
! --------------------------
             IF (LH2O) THEN
! Not an OPTRAN water channel
                KW(ILAY)= &
                     ( COEF3(25,ILAY,I)*WPRED3( 1,ILAY) ) + &
                     ( COEF3(26,ILAY,I)*WPRED3( 2,ILAY) ) + &
                     ( COEF3(27,ILAY,I)*WPRED3( 3,ILAY) ) + &
                     ( COEF3(28,ILAY,I)*WPRED3( 4,ILAY) ) + &
                     ( COEF3(29,ILAY,I)*WPRED3( 5,ILAY) ) + &
                     ( COEF3(30,ILAY,I)*WPRED3( 6,ILAY) ) + &
                     ( COEF3(31,ILAY,I)*WPRED3( 7,ILAY) ) + &
                     ( COEF3(32,ILAY,I)*WPRED3( 8,ILAY) ) + &
                     ( COEF3(33,ILAY,I)*WPRED3( 9,ILAY) ) + &
                     ( COEF3(34,ILAY,I)*WPRED3(10,ILAY) ) + &
                     ( COEF3(35,ILAY,I)*WPRED3(11,ILAY) )
!
                IF (KW(ILAY) .LT. 0.0E+0) KW(ILAY)=0.0E+0
             ENDIF
!
! Update KZFMW
             KZFMW=KZFMW + KFIX + KMET + KW(ILAY)
!
! ----------------------------------
! Calc the total layer transmittance
! ----------------------------------
!
!cccc
! This block is usually commented out and is only uncommented for
! testing purposes.
!
! kcon=0.0
! kfix=0.0
! kmet=0.0
! kw(ilay)=0.0
!cccc
! ----------------------------
! Calc change in total optical
! depth due to variable SO2
! ----------------------------
             IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0) THEN
                DKSO2=( COFSO2(1,ILAY,ISO2)*TRCPRD(1,ILAY) ) + &
                      ( COFSO2(2,ILAY,ISO2)*TRCPRD(2,ILAY) ) + &
                      ( COFSO2(3,ILAY,ISO2)*TRCPRD(3,ILAY) ) + &
                      ( COFSO2(4,ILAY,ISO2)*TRCPRD(4,ILAY) )
                DKSO2=DKSO2*SO2MLT(ILAY)
             ELSE
                DKSO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable HNO3
! ----------------------------
             IF (LHNO3 .AND. HNOMLT(ILAY) .NE. 0) THEN
                DKHNO3=( COFHNO(1,ILAY,IHNO3)*TRCPRD(1,ILAY) ) + &
                       ( COFHNO(2,ILAY,IHNO3)*TRCPRD(2,ILAY) ) + &
                       ( COFHNO(3,ILAY,IHNO3)*TRCPRD(3,ILAY) ) + &
                       ( COFHNO(4,ILAY,IHNO3)*TRCPRD(4,ILAY) )
                DKHNO3=DKHNO3*HNOMLT(ILAY)
             ELSE
                DKHNO3=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable N2O
! ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) + &
                      ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) + &
                      ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) + &
                      ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) + &
                      ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) + &
                      ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) + &
                      ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
!
! Calc the total layer effective optical depth
!cc
! this block for testing
! DKSO2=0.0
! DKHNO3=0.0
! DKN2O=0.0
!cc
! Limit -DK so it can never totally totally cancel KFIX
             DK = DKSO2 + DKHNO3 + DKN2O
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF
             KLAYER = KCON + KFIX + KMET + KW(ILAY) + DK
!
! Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
!
             KZ=KZ + KLAYER
             TAU(ILAY,J)=QIKEXP(-KLAYER)
!
          ENDDO
! End loop on levels
!
! Convert KZ to TAUZ
          TAUZ(J)=QIKEXP(-KZ)
!
       ENDDO
! End loops on channel number (frequency)
!
       RETURN
       END SUBROUTINE CALT3
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALT4 (for set4 = FCOW) version with trace gases (no SO2 or HNO3)
!
!F77====================================================================
!ROUTINE NAME:
! CALT4
!ABSTRACT:
! Calculate the transmittance for set4 using the predictors and the
! fast transmittance coefficients.
!CALL PROTOCOL:
! CALT4( LTAU, INDCHN, NLAY, BLMULT, NCHN4, CLIST4, COEF4,
! FIXMUL, CONPD4, FPRED4, CPRED4, OPRED4, WPRED4, TRCPRD, INDCO2,
! COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZ )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! LOGICAL LTAU Calc all layer trans? none
! INT arr INDCHN channel indices none
! INTEGER NLAY number of layers to bottom none
! REAL BLMULT bottom layer opt depth mult none
! INTEGER NCHN4 set4 number of channels none
! INT arr CLIST4 set4 channel list none
! REAL arr COEF4 set4 fast trans coefs various
! REAL arr FIXMUL fixed amount mult (~1.0) none
! REAL arr CONPD4 set4 H2O continuum preds various
! REAL arr FPRED4 set4 fixed gases preds various
! REAL arr CPRED4 set4 carbon monoxide preds various
! REAL arr OPRED4 set4 ozone predictors various
! REAL arr WPRED4 set4 water predictors various
! REAL arr TRCPRD trace gases pert predictors various
! INT arr INDCO2 CO2 pert chan indices none
! REAL arr COFCO2 CO2 pert coefs various
! REAL arr CO2MLT CO2 pert multiplier none
! INT arr INDN2O N2O pert chan indices none
! REAL arr COFN2O N2O pert coefs various
! REAL arr N2OMLT N2O pert multiplier none
! REAL XZ optical depth mult for TAUZ none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr TAU effective layer trans none
! REAL arr TAUZ layer-to-space trans none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2000 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The fast trans coefficents and predictors are multiplied
! together and summed to calculate the effective layer
! transmittances. Fixed, CO, ozone, and water transmittances are
! each checked individually to be sure they give 0 < trans < 1.
!
! ===================================================================
! Loops downward over all the layers for each of the NCHN4 channels
! to compute the layer transmittances TAU.
!
! The water continuum absorption coefficient is:
! k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
!
! The layer effective fixed gas absorption coefficient is:
! k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
! The layer effective CO absorption coefficient is:
! k_co = the sum i=1 to 9 of { COEF(5+8+i)*OPRED(i) }
!
! The layer effective ozone absorption coefficient is:
! k_ozone = the sum i=1 to 3 of { COEF(5+8+9+i)*OPRED(i) }
!
! The layer effective water lines absorption coefficient is:
! k_water = the sum i=1 to 11 of { COEF(5+8+9+3+i)*WPRED(i) }
!
! where
! "COEF" are the fast transmittance coefficients COEF4
! "CONPRD" are the water continuum predictors CONPRD
! "FPRED" are the fixed gases predictors FPRED4
! "CPRED" are the carbon monoxide predictors CPRED4
! "OPRED" are the ozone predictors OPRED4
! "WPRED" are the water lines predictors WPRED4
!
! The total layer effective transmittance is:
! TAU = exp( -[ k_con + k_fixed + k_co + k_ozone + k_water])
! TAU is only calc'ed if LTAU is TRUE.
!
! To help speed up the exponential calculations, we use our own
! "EXP" replacement function called QIKEXP which uses just the
! first few series expansion terms for exp(x) if x is suitably small.
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 1 Dec 1994 Scott Hannon Created
! 3 Feb 1997 Scott Hannon Re-wrote (from CALTAU) for FCOW
! 3 Sep 1997 Scott Hannon Re-wrote for sun and added TAUZ & BLMULT
! 30 Sep 1997 Scott Hannon Added variable CO2
! 27 Feb 1998 Scott Hannon Added LTAU
! 11 Aug 2000 Scott Hannon Change from 4 to 5 term H2O continuum
! 12 Sep 2002 Scott Hannon Add predictors 6 & 7 to H2O con
! 3 Jan 2003 Scott Hannon Add XZ
! 12 Oct 2004 Scott Hannon Change CO2MLT from scaler to vector
! 28 Jun 2005 Scott Hannon "trace" version for CO2,N2O
! 30 Apr 2008 Scott Hannon Change CO2 from 4 to 5 predictors
!END====================================================================
! =================================================================
       SUBROUTINE CALT4 ( LTAU, INDCHN, NLAY, BLMULT, NCHN4, CLIST4, &
          COEF4, FIXMUL, CONPD4, FPRED4, CPRED4, OPRED4, WPRED4, TRCPRD, &
          INDCO2, COFCO2, CO2MLT,INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZ )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! QIKEXP
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       LOGICAL LTAU
       INTEGER INDCHN(MXCHAN)
       INTEGER NLAY
       REAL BLMULT
       INTEGER NCHN4
       INTEGER CLIST4(MXCHN4)
       REAL COEF4(N4COEF,MAXLAY,MXCHN4)
       REAL FIXMUL(MAXLAY)
       REAL CONPD4( N4CON,MAXLAY)
       REAL FPRED4( N4FIX,MAXLAY)
       REAL CPRED4( N4CO,MAXLAY)
       REAL OPRED4( N4O3,MAXLAY)
       REAL WPRED4( N4H2O,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2( NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT(MAXLAY)
       INTEGER INDN2O(MXCHAN)
       REAL COFN2O( NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
       REAL XZ
!
! Output
       REAL TAU(MAXLAY,MXCHAN)
       REAL TAUZ(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER ICO2
       INTEGER ILAY
       INTEGER IN2O
       INTEGER J
       REAL DK
       REAL DKCO2
       REAL DKN2O
       REAL KCO
       REAL KCON
       REAL KFIX
       REAL KLAYER
       REAL KOZO
       REAL KWAT
       REAL KZ
       LOGICAL LCO2
       LOGICAL LN2O
!
! for function QIKEXP
!       REAL QIKEXP
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO I=1,NCHN4
!
! Index for TAU
          J=INDCHN( CLIST4(I) )
!
! Determine whether or not to do variable CO2
          ICO2=INDCO2( CLIST4(I) )
          IF (ICO2 .GT. 0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable N2O
          IN2O=INDN2O( CLIST4(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
!
! Initialize the layer-to-space optical depth
          KZ=0.0E+0
!
! ------------------------------
! Loop on layers (top to ground)
! ------------------------------
          DO ILAY=1,NLAY
!
! ---------------------------
! Compute the water continuum
! ---------------------------
             KCON=( COEF4(1,ILAY,I)*CONPD4(1,ILAY) ) + &
                  ( COEF4(2,ILAY,I)*CONPD4(2,ILAY) ) + &
                  ( COEF4(3,ILAY,I)*CONPD4(3,ILAY) ) + &
                  ( COEF4(4,ILAY,I)*CONPD4(4,ILAY) ) + &
                  ( COEF4(5,ILAY,I)*CONPD4(5,ILAY) ) + &
                  ( COEF4(6,ILAY,I)*CONPD4(6,ILAY) ) + &
                  ( COEF4(7,ILAY,I)*CONPD4(7,ILAY) )
!
             IF (KCON .LT. 0.0E+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 1.0E+1) THEN
                KCON=1.0E+1
             ENDIF
!
! -----------------------------
! Calc the fixed gases abs coef
! -----------------------------
             KFIX=( COEF4( 8,ILAY,I)*FPRED4( 1,ILAY) ) + &
                  ( COEF4( 9,ILAY,I)*FPRED4( 2,ILAY) ) + &
                  ( COEF4(10,ILAY,I)*FPRED4( 3,ILAY) ) + &
                  ( COEF4(11,ILAY,I)*FPRED4( 4,ILAY) ) + &
                  ( COEF4(12,ILAY,I)*FPRED4( 5,ILAY) ) + &
                  ( COEF4(13,ILAY,I)*FPRED4( 6,ILAY) ) + &
                  ( COEF4(14,ILAY,I)*FPRED4( 7,ILAY) ) + &
                  ( COEF4(15,ILAY,I)*FPRED4( 8,ILAY) ) + &
                  ( COEF4(16,ILAY,I)*FPRED4( 9,ILAY) ) + &
                  ( COEF4(17,ILAY,I)*FPRED4(10,ILAY) ) + &
                  ( COEF4(18,ILAY,I)*FPRED4(11,ILAY) )
!
             KFIX=KFIX*FIXMUL(ILAY)
!
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
!
! -----------------------
! Compute the CO abs coef
! -----------------------
             KCO=( COEF4(19,ILAY,I)*CPRED4( 1,ILAY) ) + &
                 ( COEF4(20,ILAY,I)*CPRED4( 2,ILAY) ) + &
                 ( COEF4(21,ILAY,I)*CPRED4( 3,ILAY) ) + &
                 ( COEF4(22,ILAY,I)*CPRED4( 4,ILAY) ) + &
                 ( COEF4(23,ILAY,I)*CPRED4( 5,ILAY) ) + &
                 ( COEF4(24,ILAY,I)*CPRED4( 6,ILAY) ) + &
                 ( COEF4(25,ILAY,I)*CPRED4( 7,ILAY) ) + &
                 ( COEF4(26,ILAY,I)*CPRED4( 8,ILAY) ) + &
                 ( COEF4(27,ILAY,I)*CPRED4( 9,ILAY) ) + &
                 ( COEF4(28,ILAY,I)*CPRED4(10,ILAY) ) + &
                 ( COEF4(29,ILAY,I)*CPRED4(11,ILAY) )
!
             IF (KCO .LT. 0.0E+0) THEN
                KCO=0.0E+0
             ELSEIF (KCO .GT. 1.0E+1) THEN
                KCO=1.0E+1
             ENDIF
!
! --------------------------
! Compute the ozone abs coef
! --------------------------
             KOZO=( COEF4(30,ILAY,I)*OPRED4(1,ILAY) ) + &
                  ( COEF4(31,ILAY,I)*OPRED4(2,ILAY) ) + &
                  ( COEF4(32,ILAY,I)*OPRED4(3,ILAY) )
!
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+1) THEN
                KOZO=1.0E+1
             ENDIF
!
! --------------------------
! Compute the water abs coef
! --------------------------
             KWAT=( COEF4(33,ILAY,I)*WPRED4( 1,ILAY) ) + &
                  ( COEF4(34,ILAY,I)*WPRED4( 2,ILAY) ) + &
                  ( COEF4(35,ILAY,I)*WPRED4( 3,ILAY) ) + &
                  ( COEF4(36,ILAY,I)*WPRED4( 4,ILAY) ) + &
                  ( COEF4(37,ILAY,I)*WPRED4( 5,ILAY) ) + &
                  ( COEF4(38,ILAY,I)*WPRED4( 6,ILAY) ) + &
                  ( COEF4(39,ILAY,I)*WPRED4( 7,ILAY) ) + &
                  ( COEF4(40,ILAY,I)*WPRED4( 8,ILAY) ) + &
                  ( COEF4(41,ILAY,I)*WPRED4( 9,ILAY) ) + &
                  ( COEF4(42,ILAY,I)*WPRED4(10,ILAY) ) + &
                  ( COEF4(43,ILAY,I)*WPRED4(11,ILAY) ) + &
                  ( COEF4(44,ILAY,I)*WPRED4(12,ILAY) ) + &
                  ( COEF4(45,ILAY,I)*WPRED4(13,ILAY) )
!
             IF (KWAT .LT. 0.0E+0) THEN
                KWAT=0.0E+0
             ELSEIF( KWAT .GT. 1.0E+1) THEN
                KWAT=1.0E+1
             ENDIF
!
! ----------------------------------
! Calc the total layer transmittance
! ----------------------------------
!
!cccc
! This block is usually commented out and is only uncommented for
! testing purposes.
!
! kcon=0.0E+0
! kfix=0.0E+0
! kco =0.0E+0
! kozo=0.0E+0
! kwat=0.0E+0
!cccc
!
! ----------------------------
! Calc change in total optical
! depth due to variable CO2
! ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) + &
                      ( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) + &
                      ( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) + &
                      ( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) + &
                      ( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable N2O
! ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) + &
                      ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) + &
                      ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) + &
                      ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) + &
                      ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) + &
                      ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) + &
                      ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
!
!cc
! this block for testing
! DKCO2=0.0
! DKN2O=0.0
!cc
! Limit -DK so it can never totally totally cancel KFIX
             DK = DKCO2 + DKN2O
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF
! Calc total layer optical depth
             KLAYER = KCON + KFIX + KCO + KOZO + KWAT + DK
!
! Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
!
! Calc layer-to-space optical depth
             KZ=KZ + KLAYER
!
! Calc effective layer transmittance
             IF (LTAU) TAU(ILAY,J)=QIKEXP(-KLAYER)
!
          ENDDO
! End loop on levels
!
! Convert KZ to TAUZ
          TAUZ(J)=QIKEXP(-KZ*XZ)
!
       ENDDO
! End loops on channel number (frequency)
!
       RETURN
       END SUBROUTINE CALT4
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALT5 (set5=FWO sun bfsw) version with trace gases (no SO2 or HNO3)
!
!F77====================================================================
!ROUTINE NAME:
! CALT5
!ABSTRACT:
! Calculate the transmittance for set5 using the predictors and the
! fast transmittance coefficients.
!CALL PROTOCOL:
! CALT5( LTAU, INDCHN, NLAY, BLMULT, NCHN5, CLIST5, COEF5,
! FIXMUL, CONPD5, FPRED5, WPRED5, OPRED5, TRCPRD, INDCO2, COFCO2,
! CO2MLT, INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZ )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! LOGICAL LTAU Calc all layer trans? none
! INT arr INDCHN channel indices none
! INTEGER NLAY number of layers to bottom none
! REAL BLMULT bottom layer opt depth mult none
! INTEGER NCHN5 set5 number of channels none
! INT arr CLIST5 set5 channel list none
! REAL arr COEF5 set5 fast trans coefs various
! REAL arr FIXMUL fixed amount mult (~1.0) none
! REAL arr CONPD5 set5 H2O continuum preds various
! REAL arr FPRED5 set5 fixed gases preds various
! REAL arr WPRED5 set5 water predictors various
! REAL arr OPRED5 set5 ozone predictors various
! REAL arr TRCPRD trace gas pert predictors various
! INT arr INDCO2 CO2 pert chan indices none
! REAL arr COFCO2 CO2 pert coefs various
! REAL arr CO2MLT CO2 pert multiplier none
! INT arr INDN2O N2O pert chan indices none
! REAL arr COFN2O N2O pert coefs various
! REAL arr N2OMLT N2O pert multiplier none
! REAL XZ optical depth mult for TAUZ none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr TAU effective layer trans none
! REAL arr TAUZ layer-to-space trans none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2000 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The fast trans coefficents and predictors are multiplied
! together and summed to calculate the effective layer
! transmittances. Fixed, water, and ozone transmittances are each
! checked individually to be sure they give 0 < trans < 1.
!
! ===================================================================
! Loops downward over all the layers for each of the NCHN5 channels
! to compute the layer transmittances TAU.
!
! The water continuum absorption coefficient is:
! k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
!
! The layer effective fixed gas absorption coefficient is:
! k_fixed = the sum i=1 to 11 of { COEF(5+i)*FPRED(i) }
!
! The layer effective water lines absorption coefficient is:
! k_water = the sum i=1 to 3 of { COEF(5+11+i)*WPRED(i) }
!
! The layer effective ozone absorption coefficient is:
! k_ozone = COEF(5+11+3+1)*OPRED(1)
!
! where
! "COEF" are the fast transmittance coefficients COEF5
! "CONPRD" are the water continuum predictors CONPRD
! "FPRED" are the fixed gases predictors FPRED5
! "WPRED" are the water lines predictors WPRED5
! "OPRED" are the ozone predictors OPRED5
!
! The total layer effective transmittance TAU is:
! TAU = exp( -[ k_con + k_fixed + k_water + k_ozone ])
! TAU is only calc'ed if LTAU is TRUE.
!
! To help speed up the exponential calculations, we use our own
! "EXP" replacement function called QIKEXP which uses just the
! first few series expansion terms for exp(x) if x is suitably small.
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 26 Jun 1997 Scott Hannon Created for set5
! 3 Sep 1997 Scott Hannon Added TAUZ and BLMULT
! 30 Sep 1997 Scott Hannon Added variable CO2
! 27 Feb 1998 Scott Hannon Added LTAU
! 11 Aug 2000 Scott Hannon Change from 4 to 5 term H2O continuum
! 12 Sep 2002 Scott Hannon Add predictors 6 & 7 to H2O con
! 3 Jan 2003 Scott Hannon Add XZ
! 12 Oct 2004 Scott Hannon Change CO2MLT from scaler to vector
! 28 Jun 2005 Scott Hannon "trace" version for CO2,N2O
! 30 Apr 2008 Scott Hannon Change CO2 from 4 to 5 predictors
!END====================================================================
! =================================================================
       SUBROUTINE CALT5 ( LTAU, INDCHN, NLAY, BLMULT, NCHN5, CLIST5, &
          COEF5, FIXMUL, CONPD5, FPRED5, WPRED5, OPRED5, TRCPRD, INDCO2, &
          COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZ )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! QIKEXP
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       LOGICAL LTAU
       INTEGER INDCHN(MXCHAN)
       INTEGER NLAY
       REAL BLMULT
       INTEGER NCHN5
       INTEGER CLIST5(MXCHN5)
       REAL COEF5(N5COEF,MAXLAY,MXCHN5)
       REAL FIXMUL(MAXLAY)
       REAL CONPD5( N5CON,MAXLAY)
       REAL FPRED5( N5FIX,MAXLAY)
       REAL WPRED5( N5H2O,MAXLAY)
       REAL OPRED5( N5O3,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2( NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT(MAXLAY)
       INTEGER INDN2O(MXCHAN)
       REAL COFN2O( NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
       REAL XZ
!
! Output
       REAL TAU(MAXLAY,MXCHAN)
       REAL TAUZ(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER ICO2
       INTEGER ILAY
       INTEGER IN2O
       INTEGER J
       REAL DK
       REAL DKCO2
       REAL DKN2O
       REAL KCON
       REAL KFIX
       REAL KLAYER
       REAL KOZO
       REAL KWAT
       REAL KZ
       LOGICAL LCO2
       LOGICAL LN2O
!
! for function QIKEXP
!       REAL QIKEXP
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO I=1,NCHN5
!
! Index for TAU
          J=INDCHN( CLIST5(I) )
!
! Determine whether or not to do variable CO2
          ICO2=INDCO2( CLIST5(I) )
          IF (ICO2 .GT. 0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable CO2
          IN2O=INDN2O( CLIST5(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
!
! Initialize the layer-to-space optical depth
          KZ=0.0E+0
!
! ------------------------------
! Loop on layers (top to ground)
! ------------------------------
          DO ILAY=1,NLAY
!
! ---------------------------
! Compute the water continuum
! ---------------------------
             KCON=( COEF5(1,ILAY,I)*CONPD5(1,ILAY) ) + &
                  ( COEF5(2,ILAY,I)*CONPD5(2,ILAY) ) + &
                  ( COEF5(3,ILAY,I)*CONPD5(3,ILAY) ) + &
                  ( COEF5(4,ILAY,I)*CONPD5(4,ILAY) ) + &
                  ( COEF5(5,ILAY,I)*CONPD5(5,ILAY) ) + &
                  ( COEF5(6,ILAY,I)*CONPD5(6,ILAY) ) + &
                  ( COEF5(7,ILAY,I)*CONPD5(7,ILAY) )
!
             IF (KCON .LT. 0.0E+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 1.0E+1) THEN
                KCON=1.0E+1
             ENDIF
!
! -----------------------------
! Calc the fixed gases abs coef
! -----------------------------
             KFIX=( COEF5( 8,ILAY,I)*FPRED5( 1,ILAY) ) + &
                  ( COEF5( 9,ILAY,I)*FPRED5( 2,ILAY) ) + &
                  ( COEF5(10,ILAY,I)*FPRED5( 3,ILAY) ) + &
                  ( COEF5(11,ILAY,I)*FPRED5( 4,ILAY) ) + &
                  ( COEF5(12,ILAY,I)*FPRED5( 5,ILAY) ) + &
                  ( COEF5(13,ILAY,I)*FPRED5( 6,ILAY) ) + &
                  ( COEF5(14,ILAY,I)*FPRED5( 7,ILAY) ) + &
                  ( COEF5(15,ILAY,I)*FPRED5( 8,ILAY) ) + &
                  ( COEF5(16,ILAY,I)*FPRED5( 9,ILAY) ) + &
                  ( COEF5(17,ILAY,I)*FPRED5(10,ILAY) ) + &
                  ( COEF5(18,ILAY,I)*FPRED5(11,ILAY) )
!
             KFIX=KFIX*FIXMUL(ILAY)
!
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
!
!
! --------------------------
! Compute the water abs coef
! --------------------------
             KWAT=( COEF5(19,ILAY,I)*WPRED5( 1,ILAY) ) + &
                  ( COEF5(20,ILAY,I)*WPRED5( 2,ILAY) ) + &
                  ( COEF5(21,ILAY,I)*WPRED5( 3,ILAY) )
!
             IF (KWAT .LT. 0.0E+0) THEN
                KWAT=0.0E+0
             ELSEIF( KWAT .GT. 1.0E+1) THEN
                KWAT=1.0E+1
             ENDIF
!
!
! --------------------------
! Compute the ozone abs coef
! --------------------------
             KOZO=( COEF5(22,ILAY,I)*OPRED5(1,ILAY) )
!
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+1) THEN
                KOZO=1.0E+1
             ENDIF
!
!
! ----------------------------------
! Calc the total layer transmittance
! ----------------------------------
!
!cccc
! This block is usually commented out and is only uncommented for
! testing purposes.
!
! kcon=0.0E+0
! kfix=0.0E+0
! kwat=0.0E+0
! kozo=0.0E+0
!cccc
!
! ----------------------------
! Calc change in total optical
! depth due to variable CO2
! ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0.0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) + &
                      ( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) + &
                      ( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) + &
                      ( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) + &
                      ( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable N2O
! ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0.0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) + &
                      ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) + &
                      ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) + &
                      ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) + &
                      ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) + &
                      ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) + &
                      ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
!
!cc
! this block for testing
! DKCO2=0.0
! DKN2O=0.0
!cc
! Limit -DK so it can never totally totally cancel KFIX
             DK = DKCO2 + DKN2O
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF
! Calc total layer optical depth
             KLAYER = KCON + KFIX + KWAT + KOZO + DK
!
! Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
!
! Calc layer-to-space optical depth
             KZ=KZ + KLAYER
!
! Calc effective layer transmittance
             IF (LTAU) TAU(ILAY,J)=QIKEXP(-KLAYER)
!
          ENDDO
! End loop on levels
!
! Convert KZ to TAUZ
          TAUZ(J)=QIKEXP(-KZ*XZ)
!
       ENDDO
! End loops on channel number (frequency)
!
       RETURN
       END SUBROUTINE CALT5
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALT6 (set6=FWO sun mfmw) version for trace gases (no HNO3)
!
!F77====================================================================
!ROUTINE NAME:
! CALT6
!ABSTRACT:
! Calculate the transmittance for set6 using the predictors and the
! fast transmittance coefficients.
!CALL PROTOCOL:
! CALT6( LTAU, INDCHN, NLAY, BMULT, NCHN6, CLIST6, COEF6,
! FIXMUL, CONPD6, FPRED6, WPRED6, OPRED6, TRCPRD,
! INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
! INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZ )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! LOGICAL LTAU Calc all layer trans? none
! INT arr INDCHN channel indices none
! INTEGER NLAY number of layers to bottom none
! REAL BLMULT bottom layer opt depth mult none
! INTEGER NCHN6 set6 number of channels none
! INT arr CLIST6 set6 channel list none
! REAL arr COEF6 set6 fast trans coefs various
! REAL arr FIXMUL fixed amount mult (~1.0) none
! REAL arr CONPD6 set6 H2O continuum preds various
! REAL arr FPRED6 set6 fixed gases preds various
! REAL arr WPRED6 set6 water predictors various
! REAL arr OPRED6 set6 ozone predictors various
! REAL arr TRCPRD trace gas pert predictors various
! INT arr INDCO2 CO2 pert chan indices none
! REAL arr COFCO2 CO2 pert coefs various
! REAL arr CO2MLT CO2 pert multiplier none
! REAL arr SO2PRD SO2 pert predictors various
! INT arr INDSO2 SO2 pert chan indices none
! REAL arr COFSO2 SO2 pert coefs various
! REAL arr SO2MLT SO2 pert multiplier none
! INT arr INDN2O N2O pert chan indices none
! REAL arr COFN2O N2O pert coefs various
! REAL arr N2OMLT N2O pert multiplier none
! REAL XZ optical depth mult for TAUZ none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr TAU effective layer trans none
! REAL arr TAUZ layer-to-space trans none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2000 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The fast trans coefficents and predictors are multiplied
! together and summed to calculate the effective layer
! transmittances. Fixed, water, and ozone transmittances are each
! checked individually to be sure they give 0 < trans < 1.
!
! ===================================================================
! Loops downward over all the layers for each of the NCHN6 channels
! to compute the layer transmittances TAU.
!
! The water continuum absorption coefficient is:
! k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
!
! The layer effective fixed gas absorption coefficient is:
! k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
! The layer effective water lines absorption coefficient is:
! k_water = the sum i=1 to 7 of { COEF(5+8+i)*WPRED(i) }
!
! The layer effective ozone absorption coefficient is:
! k_ozone = COEF(5+8+7+1)*OPRED(1)
!
! where
! "COEF" are the fast transmittance coefficients COEF5
! "CONPRD" are the water continuum predictors CONPRD
! "FPRED" are the fixed gases predictors FPRED5
! "WPRED" are the water lines predictors WPRED5
! "OPRED" are the ozone predictors OPRED5
!
! The total layer effective transmittance TAU is:
! TAU = exp( -[ k_con + k_fixed + k_water + k_ozone ])
! TAU is only calc'ed if LTAU is TRUE.
!
! To help speed up the exponential calculations, we use our own
! "EXP" replacement function called QIKEXP which uses just the
! first few series expansion terms for exp(x) if x is suitably small.
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 3 Jul 1997 Scott Hannon Created for set6
! 3 Sep 1997 Scott Hannon Added TAUZ and BLMULT
! 30 Sep 1997 Scott Hannon Added variable CO2
! 27 Feb 1998 Scott Hannon Added LTAU
! 11 Aug 2000 Scott Hannon Change from 4 to 5 term H2O continuum
! 12 Sep 2002 Scott Hannon Add predictors 6 & 7 to H2O con
! 3 Jan 2003 Scott Hannon Add XZ
! 6 Feb 2003 Scott Hannon Bug fix - ozone use coef 23, not 21
! 25 Apr 2003 Scott Hannon Add SO2
! 28 Jun 2005 Scott Hannon "trace" version with CO2,SO2,N2O
! 30 Apr 2008 Scott Hannon Change CO2 from 4 to 5 predictors
!END====================================================================
! =================================================================
       SUBROUTINE CALT6 ( LTAU, INDCHN, NLAY, BLMULT, NCHN6, CLIST6, &
          COEF6, FIXMUL, CONPD6, FPRED6, WPRED6, OPRED6, TRCPRD, &
          INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
          INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZ)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! QIKEXP
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       LOGICAL LTAU
       INTEGER INDCHN(MXCHAN)
       INTEGER NLAY
       REAL BLMULT
       INTEGER NCHN6
       INTEGER CLIST6(MXCHN6)
       REAL COEF6(N6COEF,MAXLAY,MXCHN6)
       REAL FIXMUL(MAXLAY)
       REAL CONPD6( N6CON,MAXLAY)
       REAL FPRED6( N6FIX,MAXLAY)
       REAL WPRED6( N6H2O,MAXLAY)
       REAL OPRED6( N6O3,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2( NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT(MAXLAY)
       INTEGER INDSO2(MXCHAN)
       REAL COFSO2( NSO2,MAXLAY,MXCHNS)
       REAL SO2MLT(MAXLAY)
       INTEGER INDN2O(MXCHAN)
       REAL COFN2O( NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
       REAL XZ
!
! Output
       REAL TAU(MAXLAY,MXCHAN)
       REAL TAUZ(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER ICO2
       INTEGER ILAY
       INTEGER IN2O
       INTEGER ISO2
       INTEGER J
       REAL DK
       REAL DKCO2
       REAL DKN2O
       REAL DKSO2
       REAL KCON
       REAL KFIX
       REAL KLAYER
       REAL KOZO
       REAL KWAT
       REAL KZ
       LOGICAL LCO2
       LOGICAL LN2O
       LOGICAL LSO2
!
! for function QIKEXP
!       REAL QIKEXP
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO I=1,NCHN6
!
! Index for TAU
          J=INDCHN( CLIST6(I) )
!
! Determine whether or not to do variable CO2
          ICO2=INDCO2( CLIST6(I) )
          IF (ICO2 .GT. 0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable SO2
          ISO2=INDSO2( CLIST6(I) )
          IF (ISO2 .GT. 0) THEN
             LSO2=.TRUE.
          ELSE
             LSO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable N2O
          IN2O=INDN2O( CLIST6(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
!
! Initialize the layer-to-space optical depth
          KZ=0.0E+0
!
! ------------------------------
! Loop on layers (top to ground)
! ------------------------------
          DO ILAY=1,NLAY
!
! ---------------------------
! Compute the water continuum
! ---------------------------
             KCON=( COEF6(1,ILAY,I)*CONPD6(1,ILAY) ) + &
                  ( COEF6(2,ILAY,I)*CONPD6(2,ILAY) ) + &
                  ( COEF6(3,ILAY,I)*CONPD6(3,ILAY) ) + &
                  ( COEF6(4,ILAY,I)*CONPD6(4,ILAY) ) + &
                  ( COEF6(5,ILAY,I)*CONPD6(5,ILAY) ) + &
                  ( COEF6(6,ILAY,I)*CONPD6(6,ILAY) ) + &
                  ( COEF6(7,ILAY,I)*CONPD6(7,ILAY) )
!
             IF (KCON .LT. 0.0E+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 1.0E+1) THEN
                KCON=1.0E+1
             ENDIF
!
! -----------------------------
! Calc the fixed gases abs coef
! -----------------------------
             KFIX=( COEF6( 8,ILAY,I)*FPRED6( 1,ILAY) ) + &
                  ( COEF6( 9,ILAY,I)*FPRED6( 2,ILAY) ) + &
                  ( COEF6(10,ILAY,I)*FPRED6( 3,ILAY) ) + &
                  ( COEF6(11,ILAY,I)*FPRED6( 4,ILAY) ) + &
                  ( COEF6(12,ILAY,I)*FPRED6( 5,ILAY) ) + &
                  ( COEF6(13,ILAY,I)*FPRED6( 6,ILAY) ) + &
                  ( COEF6(14,ILAY,I)*FPRED6( 7,ILAY) ) + &
                  ( COEF6(15,ILAY,I)*FPRED6( 8,ILAY) )
!
             KFIX=KFIX*FIXMUL(ILAY)
!
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
!
!
! --------------------------
! Compute the water abs coef
! --------------------------
             KWAT=( COEF6(16,ILAY,I)*WPRED6( 1,ILAY) ) + &
                  ( COEF6(17,ILAY,I)*WPRED6( 2,ILAY) ) + &
                  ( COEF6(18,ILAY,I)*WPRED6( 3,ILAY) ) + &
                  ( COEF6(19,ILAY,I)*WPRED6( 4,ILAY) ) + &
                  ( COEF6(20,ILAY,I)*WPRED6( 5,ILAY) ) + &
                  ( COEF6(21,ILAY,I)*WPRED6( 6,ILAY) ) + &
                  ( COEF6(22,ILAY,I)*WPRED6( 7,ILAY) )
!
             IF (KWAT .LT. 0.0E+0) THEN
                KWAT=0.0E+0
             ELSEIF( KWAT .GT. 1.0E+1) THEN
                KWAT=1.0E+1
             ENDIF
!
!
! --------------------------
! Compute the ozone abs coef
! --------------------------
             KOZO=( COEF6(23,ILAY,I)*OPRED6(1,ILAY) )
!
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+1) THEN
                KOZO=1.0E+1
             ENDIF
!
!
! ----------------------------------
! Calc the total layer transmittance
! ----------------------------------
!
!cccc
! This block is usually commented out and is only uncommented for
! testing purposes.
!
! kcon=0.0E+0
! kfix=0.0E+0
! kwat=0.0E+0
! kozo=0.0E+0
!cccc
!
! ----------------------------
! Calc change in total optical
! depth due to variable CO2
! ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0.0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) + &
                      ( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) + &
                      ( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) + &
                      ( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) + &
                      ( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable SO2
! ----------------------------
             IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0.0) THEN
                DKSO2=( COFSO2(1,ILAY,ISO2)*TRCPRD(1,ILAY) ) + &
                      ( COFSO2(2,ILAY,ISO2)*TRCPRD(2,ILAY) ) + &
                      ( COFSO2(3,ILAY,ISO2)*TRCPRD(3,ILAY) ) + &
                      ( COFSO2(4,ILAY,ISO2)*TRCPRD(4,ILAY) )
                DKSO2=DKSO2*SO2MLT(ILAY)
             ELSE
                DKSO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable N2O
! ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0.0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) + &
                      ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) + &
                      ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) + &
                      ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) + &
                      ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) + &
                      ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) + &
                      ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
!
!cc
! this block for testing
! DKCO2=0.0
! DKSO2=0.0
! DKN2O=0.0
!cc
! Limit -DK so it can never totally totally cancel KFIX
             DK = DKCO2 + DKSO2 + DKN2O
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF
! Calc total layer optical depth
             KLAYER = KCON + KFIX + KWAT + KOZO + DK
!
! Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
!
! Calc layer-to-space optical depth
             KZ=KZ + KLAYER
!
! Calc effective layer transmittance
             IF (LTAU) TAU(ILAY,J)=QIKEXP(-KLAYER)
!
          ENDDO
! End loop on levels
!
! Convert KZ to TAUZ
          TAUZ(J)=QIKEXP(-KZ*XZ)
!
       ENDDO
! End loops on channel number (frequency)
!
       RETURN
       END SUBROUTINE CALT6
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! CALT7 (set7=FWO sun mfbw) version for trace gases (no SO2 or HNO3)
!
!F77====================================================================
!ROUTINE NAME:
! CALT7
!ABSTRACT:
! Calculate the transmittance for set7 using the predictors and the
! fast transmittance coefficients.
!CALL PROTOCOL:
! CALT7( LTAU, INDCHN, NLAY, BLMULT, NCHN7, CLIST7, COEF7,
! FIXMUL, CONPD7, FPRED7, WPRED7, OPRED7, TRCPRD, INDCO2, COFCO2,
! CO2MLT, INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZ )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! LOGICAL LTAU calc all layer trans? none
! INT arr INDCHN channel indices none
! INTEGER NLAY number of layers to bottom none
! REAL BLMULT bottom layer opt depth mult none
! INTEGER NCHN7 set7 number of channels none
! INT arr CLIST7 set7 channel list none
! REAL arr COEF7 set7 fast trans coefs various
! REAL arr FIXMUL fixed amount mult (~1.0) none
! REAL arr CONPD7 set7 H2O continuum preds various
! REAL arr FPRED7 set7 fixed gases preds various
! REAL arr WPRED7 set7 water predictors various
! REAL arr OPRED7 set7 ozone predictors various
! REAL arr TRCPRD trace gas pert predictors various
! INT arr INDCO2 CO2 pert chan indices none
! REAL arr COFCO2 CO2 pert coefs various
! REAL arr CO2MLT CO2 pert multiplier none
! INT arr INDN2O N2O pert chan indices none
! REAL arr COFN2O N2O pert coefs various
! REAL arr N2OMLT N2O pert multiplier none
! REAL XZ optical depth mult for TAUZ none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr TAU effective layer trans none
! REAL arr TAUZ layer-to-space trans none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2000 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! The fast trans coefficents and predictors are multiplied
! together and summed to calculate the effective layer
! transmittances. Fixed, water, and ozone transmittances are each
! checked individually to be sure they give 0 < trans < 1.
!
! ===================================================================
! Loops downward over all the layers for each of the NCHN7 channels
! to compute the layer transmittances TAU.
!
! The water continuum absorption coefficient is:
! k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
!
! The layer effective fixed gas absorption coefficient is:
! k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
! The layer effective water lines absorption coefficient is:
! k_water = the sum i=1 to 13 of { COEF(5+8+i)*WPRED(i) }
!
! The layer effective ozone absorption coefficient is:
! k_ozone = COEF(5+8+13+1)*OPRED(1)
!
! where
! "COEF" are the fast transmittance coefficients COEF4
! "CONPRD" are the water continuum predictors CONPRD
! "FPRED" are the fixed gases predictors FPRED
! "WPRED" are the water lines predictors WPRED4
! "OPRED" are the ozone predictors OPRED4
!
! The total layer effective transmittance TAU is:
! TAU = exp( -[ k_con + k_fixed + k_water + k_ozone ])
! TAU is only calc'ed if LTAU is TRUE.
!
! To help speed up the exponential calculations, we use our own
! "EXP" replacement function called QIKEXP which uses just the
! first few series expansion terms for exp(x) if x is suitably small.
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 7 Jul 1997 Scott Hannon Created for set7
! 3 Sep 1997 Scott Hannon Added TAUZ and BLMULT
! 30 Sep 1997 Scott Hannon Added variable CO2
! 27 Feb 1998 Scott Hannon Added LTAU
! 11 Aug 2000 Scott Hannon Change from 4 to 5 term H2O continuum
! 12 Sep 2002 Scott Hannon Add predictors 6 & 7 to H2O con
! 3 Jan 2003 Scott Hannon Add XZ
! 12 Oct 2004 Scott Hannon Change CO2MLT from scaler to vector
! 28 Jun 2005 Scott Hannon "trace" version for CO2,N2O
! 30 Apr 2008 Scott Hannon Change CO2 from 4 to 5 predictors
!END====================================================================
! =================================================================
       SUBROUTINE CALT7 ( LTAU, INDCHN, NLAY, BLMULT, NCHN7, CLIST7, &
          COEF7, FIXMUL, CONPD7, FPRED7, WPRED7, OPRED7, TRCPRD, INDCO2, &
          COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZ )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! QIKEXP
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       LOGICAL LTAU
       INTEGER INDCHN(MXCHAN)
       INTEGER NLAY
       REAL BLMULT
       INTEGER NCHN7
       INTEGER CLIST7(MXCHN7)
       REAL COEF7(N7COEF,MAXLAY,MXCHN7)
       REAL FIXMUL(MAXLAY)
       REAL CONPD7( N7CON,MAXLAY)
       REAL FPRED7( N7FIX,MAXLAY)
       REAL WPRED7( N7H2O,MAXLAY)
       REAL OPRED7( N7O3,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2( NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT(MAXLAY)
       INTEGER INDN2O(MXCHAN)
       REAL COFN2O( NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
       REAL XZ
!
! Output
       REAL TAU(MAXLAY,MXCHAN)
       REAL TAUZ(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER ICO2
       INTEGER ILAY
       INTEGER IN2O
       INTEGER J
       REAL DK
       REAL DKCO2
       REAL DKN2O
       REAL KCON
       REAL KFIX
       REAL KOZO
       REAL KLAYER
       REAL KWAT
       REAL KZ
       LOGICAL LCO2
       LOGICAL LN2O
!
! for function QIKEXP
!       REAL QIKEXP
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO I=1,NCHN7
!
! Index for TAU
          J=INDCHN( CLIST7(I) )
!
! Determine whether or not to do variable CO2
          ICO2=INDCO2( CLIST7(I) )
          IF (ICO2 .GT. 0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
!
! Determine whether or not to do variable CO2
          IN2O=INDN2O( CLIST7(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
!
! Initialize the layer-to-space optical depth
          KZ=0.0E+0
!
! ------------------------------
! Loop on layers (top to ground)
! ------------------------------
          DO ILAY=1,NLAY
!
! ---------------------------
! Compute the water continuum
! ---------------------------
             KCON=( COEF7(1,ILAY,I)*CONPD7(1,ILAY) ) + &
                  ( COEF7(2,ILAY,I)*CONPD7(2,ILAY) ) + &
                  ( COEF7(3,ILAY,I)*CONPD7(3,ILAY) ) + &
                  ( COEF7(4,ILAY,I)*CONPD7(4,ILAY) ) + &
                  ( COEF7(5,ILAY,I)*CONPD7(5,ILAY) ) + &
                  ( COEF7(6,ILAY,I)*CONPD7(6,ILAY) ) + &
                  ( COEF7(7,ILAY,I)*CONPD7(7,ILAY) )
!
             IF (KCON .LT. 0.0E+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 1.0E+1) THEN
                KCON=1.0E+1
             ENDIF
!
! -----------------------------
! Calc the fixed gases abs coef
! -----------------------------
             KFIX=( COEF7( 8,ILAY,I)*FPRED7( 1,ILAY) ) + &
                  ( COEF7( 9,ILAY,I)*FPRED7( 2,ILAY) ) + &
                  ( COEF7(10,ILAY,I)*FPRED7( 3,ILAY) ) + &
                  ( COEF7(11,ILAY,I)*FPRED7( 4,ILAY) ) + &
                  ( COEF7(12,ILAY,I)*FPRED7( 5,ILAY) ) + &
                  ( COEF7(13,ILAY,I)*FPRED7( 6,ILAY) ) + &
                  ( COEF7(14,ILAY,I)*FPRED7( 7,ILAY) ) + &
                  ( COEF7(15,ILAY,I)*FPRED7( 8,ILAY) )
!
             KFIX=KFIX*FIXMUL(ILAY)
!
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
!
!
! --------------------------
! Compute the water abs coef
! --------------------------
             KWAT=( COEF7(16,ILAY,I)*WPRED7( 1,ILAY) ) + &
                  ( COEF7(17,ILAY,I)*WPRED7( 2,ILAY) ) + &
                  ( COEF7(18,ILAY,I)*WPRED7( 3,ILAY) ) + &
                  ( COEF7(19,ILAY,I)*WPRED7( 4,ILAY) ) + &
                  ( COEF7(20,ILAY,I)*WPRED7( 5,ILAY) ) + &
                  ( COEF7(21,ILAY,I)*WPRED7( 6,ILAY) ) + &
                  ( COEF7(22,ILAY,I)*WPRED7( 7,ILAY) ) + &
                  ( COEF7(23,ILAY,I)*WPRED7( 8,ILAY) ) + &
                  ( COEF7(24,ILAY,I)*WPRED7( 9,ILAY) ) + &
                  ( COEF7(25,ILAY,I)*WPRED7(10,ILAY) ) + &
                  ( COEF7(26,ILAY,I)*WPRED7(11,ILAY) ) + &
                  ( COEF7(27,ILAY,I)*WPRED7(12,ILAY) ) + &
                  ( COEF7(28,ILAY,I)*WPRED7(13,ILAY) )
!
             IF (KWAT .LT. 0.0E+0) THEN
                KWAT=0.0E+0
             ELSEIF( KWAT .GT. 1.0E+1) THEN
                KWAT=1.0E+1
             ENDIF
!
!
! --------------------------
! Compute the ozone abs coef
! --------------------------
             KOZO=( COEF7(29,ILAY,I)*OPRED7(1,ILAY) )
!
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+1) THEN
                KOZO=1.0E+1
             ENDIF
!
!
! ----------------------------------
! Calc the total layer transmittance
! ----------------------------------
!
!cccc
! This block is usually commented out and is only uncommented for
! testing purposes.
!
! kcon=0.0E+0
! kfix=0.0E+0
! kwat=0.0E+0
! kozo=0.0E+0
!cccc
!
! ----------------------------
! Calc change in total optical
! depth due to variable CO2
! ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0.0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) + &
                      ( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) + &
                      ( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) + &
                      ( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) + &
                      ( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
!
! ----------------------------
! Calc change in total optical
! depth due to variable N2O
! ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0.0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) + &
                      ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) + &
                      ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) + &
                      ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) + &
                      ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) + &
                      ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) + &
                      ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
!
!cc
! this block for testing
! DKCO2=0.0
! DKN2O=0.0
!cc
! Limit -DK so it can never totally totally cancel KFIX
             DK = DKCO2 + DKN2O
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF
! Calc total layer optical depth
             KLAYER=KCON + KFIX + KWAT + KOZO + DK
!
! Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
!
! Calc layer-to-space optical depth
             KZ=KZ + KLAYER
!
! Calc effective layer transmittance
             IF (LTAU) TAU(ILAY,J)=QIKEXP(-KLAYER)
!
          ENDDO
! End loop on levels
!
! Convert KZ to TAUZ
          TAUZ(J)=QIKEXP(-KZ*XZ)
!
       ENDDO
! End loops on channel number (frequency)
!
       RETURN
       END SUBROUTINE CALT7
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! FAKETZ
!
!F77====================================================================
!ROUTINE NAME:
! FAKETZ
!ABSTRACT:
! Calculate a "fake" layer-to-space transmittance.
!CALL PROTOCOL:
! FAKETZ ( NFAKE, INDFAK, TAUZ, SEC, SECFAK, TAUZFK )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER INDFAK array indices for fake none
! INTEGER NFAKE number of fake points none
! REAL SEC angle secant for TAUZ none
! REAL SECFAK angle secant for TAUZFK none
! REAL arr TAUZ layer-to-space trans none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr TAUZFK fake layer-to-space trans none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! March 1998 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! A "fake" layer-to-space transmittance is calculated for some
! arbitrary angle by scaling the optical depth by the ratio of
! the angle secants. The exact form of the calculation is:
! TAUZFK = EXP( LN(TAUZ) * SECFAK/SEC )
! This is a crude approximation of the correct value.
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! This is a crude approximation of the correct value.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! Aug 27 1997 Scott Hannon Created
! Aug 27 1998 Scott Hannon Fix bug for case when TAUZ=0
!END====================================================================
! =================================================================
       SUBROUTINE FAKETZ ( NFAKE, INDFAK, TAUZ, SEC, SECFAK, TAUZFK )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
       INTEGER NFAKE
       INTEGER INDFAK(MXCHAN)
       REAL TAUZ(MXCHAN)
       REAL SEC
       REAL SECFAK
       REAL TAUZFK(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER ICHAN
       REAL RATSEC
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Calc ratio of secants
       RATSEC=SECFAK/SEC
!
! ---------------------------
! Loop on channel (frequency)
! ---------------------------
       DO I=1,NFAKE
!
          ICHAN=INDFAK(I)
! Be careful to avoid log(0)
          IF (TAUZ(ICHAN) .GT. 1E-8) THEN
             TAUZFK(ICHAN)=EXP( RATSEC*LOG( TAUZ(ICHAN) ) )
          ELSE
             TAUZFK(ICHAN)=1E-10
          ENDIF
!
       ENDDO
!
       RETURN
       END SUBROUTINE FAKETZ
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! GETBOT
!
!F77====================================================================
!ROUTINE NAME:
! GETBOT
!ABSTRACT:
! Calculate the bottom layer number and fractional multiplier
! based on the supplied surface pressure and temperature profile
!CALL PROTOCOL:
! GETBOT( NLAY, PLEV, PSURF, LBOT, BLMULT )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER NLAY number of profile layers none
! REAL arr PLEV layer pres level boundaries mb
! REAL PSURF surface pressure mb
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER LBOT bottom layer number none
! REAL BLMULT bot layer fractional mult none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! May 2001 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon/H.Motteler.
!
! This routine starts at layer 100 and loops upward until it finds
! the layer bounding PSURF. It then computes the fraction of this
! bottom layer above PSURF. A bottom layer thinner than 5% of the
! full layer thickness is avoided; in this case the layer directly
! above is used instead with a fraction slightly larger than 1.
!
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! Assumes the user has supplied vaguely realistic profile amounts
! and temperatures.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- --------------------------------------
! 31 Mar 2000 Scott Hannon Created
! 1 May 2001 Scott Hannon Add DELPX and check bottom thickness
! 2 May 2001 Scott Hannon PLEV changed from local data to input
! 17 Dec 2004 Scott Hannon Add NLAY to call; trap LBOT>NLAY;
! add warning for excessive BLMULT
! 24 Jun 2005 Scott Hannon "10" loop changed to start on DELPX
! assignment rather than IF line below.
!END====================================================================
! =================================================================
       SUBROUTINE GETBOT ( NLAY, PLEV, PSURF, LBOT, BLMULT )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER,intent(in) :: NLAY ! number of profile layers
       REAL,intent(in) :: PLEV(MAXLAY+1) ! layer pressure level boundaries
       REAL,intent(in) :: PSURF ! surface pressure
!
! Output
       INTEGER,intent(out) :: LBOT ! bottom layer number
       REAL,intent(out) :: BLMULT ! bottom layer fractional multiplier
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER LBOTX ! unrestricted bottom layer number
       REAL DELPX ! 5% of layer thickness in pressure
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Determine LBOT from PSURF by comparing to PLEV

!  
! Still don't understand BLMULT, but at least giving 1.e0 does not
! modify any parameters (e.g., optical thickness layer-to-space at the bottome layer).
! So just leave as 1.e0 now. 
! Toshi Matsui Oct 2013
!
 LBOT   = NLAY
 BLMULT = 1.e0
 RETURN  !skipping the original process below.

       LBOTX=MAXLAY
 10 DELPX=0.05*( PLEV(LBOTX+1) - PLEV(LBOTX) )

       IF (PSURF .LT. PLEV(LBOTX)+DELPX) THEN
          LBOTX=LBOTX - 1
          GOTO 10
       ENDIF
!
       IF (LBOTX .GT. NLAY) THEN
          LBOT = NLAY
       ELSE
          LBOT = LBOTX
       ENDIF


!
! Calc bottom layer multiplier (fractional layer)
       BLMULT = (PSURF - PLEV(LBOT))/(PLEV(LBOT+1) - PLEV(LBOT))

       IF (BLMULT .GT. 1.3) THEN
          WRITE(IOINFO,1010) BLMULT, LBOTX, NLAY
 1010 FORMAT('WARNING :: excessive BLMULT=',F6.3,'; optimal LBOT=', &
          I3,' but layers end at NLAY=',I3)
       ENDIF
!
       RETURN
       END SUBROUTINE GETBOT
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! MEAN_T
!
!F77====================================================================
!ROUTINE NAME:
! MEAN_T
!ABSTRACT:
! Convert PGE L2 pseudo-level T's to layer mean T's.
!CALL PROTOCOL
! MEAN_T( LBOT, PLEV, PSURF, TPSEUD, TLAY )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER LBOT number of layers none
! REAL arr PLEV pressure levels millibar
! REAL PSURF surface pressure millibar
! REAL TPSEUD pseudo level temperatures Kelvin
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr TLAY layer mean temperature Kelvin
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! SARTA
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! none
!COMMON BLOCKS
! none
!DESCRIPTION:
! Converts the AIRS PGE Level2 pseudo-level temperature profile
! to the layer mean temperature profile as required by the RTA.
! The layer mean T's are defined in the PGE as the average of
! the two adjacent level pseudo-level T's.
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! SARTA assumes an input layer profile is full layers, but
! this routine has a built in adjustment for the bottom
! factional layer T. Thus SARTA needs to be aware of this
! and skip the usual adjustment for the bottom layer.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 05 Feb 2004 Scott Hannon Created; based on translation formula
! from the PGE's "meantemp.F"
!
!
!END =====================================================================
!
! =================================================================
       SUBROUTINE MEAN_T(LBOT, PLEV, PSURF, TPSEUD, TLAY)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER LBOT ! # of layers for this profile
       REAL PSURF ! surface pressure
       REAL PLEV(MAXLAY+1) ! pressure levels
       REAL TPSEUD(MAXLAY) ! pseudo level temperature at PLEV(L+1)
!
! Output
       REAL TLAY(MAXLAY) ! layer mean temperature
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER L ! layer index
       REAL TSURFA ! Air temperature at surface
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
! Do top layer (special case)
       TLAY(1) = TPSEUD(1)
!
! Loop down over the layers
       DO L=2,LBOT-1
          TLAY(L) = 0.5*( TPSEUD(L-1) + TPSEUD(L) )
       ENDDO
!
! Interpolate to get air temperature at the surface
       TSURFA = TPSEUD(LBOT-1) + ( TPSEUD(LBOT) - TPSEUD(LBOT-1) )* &
          ( PSURF - PLEV(LBOT) )/( PLEV( LBOT+1) - PLEV(LBOT) )
!
! Do bottom layer (special case)
       TLAY(LBOT) = 0.5*( TPSEUD(LBOT-1) + TSURFA )
!
       RETURN
       END SUBROUTINE MEAN_T
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
! version for sarta
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! OPNRTP version with trace gases
!
!F77====================================================================
!ROUTINE NAME: OPNRTP
!ABSTRACT:
! Open and check input RTP file.
!CALL PROTOCOL:
! OPNRTP(FIN, LRHOT, PTYPE, NCHAN, FCHAN, LSTCHN, INDCHN,
! IH2O, IO3, ICO, ICH4, ICO2, ISO2, IHNO3, IN2O,
! IOPCI, HEAD, HATT, PATT, LCO2PM)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! CHAR*80 FIN input RTP file name none
! LOGICAL LRHOT force refl therm rho? none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER PTYPE profile type none
! INTEGER NCHAN number of channels none
! INTEGER FCHAN channel frequencies cm^-1
! INTEGER LSTCHN list of channel numbers none (1-2378)
! INTEGER INDCHN indices of channels none
! INTEGER IH2O index of H2O in gamnt none
! INTEGER IO3 index of O3 in gamnt none
! INTEGER ICO index of CO in gamnt none
! INTEGER ICH4 index of CH4 in gamnt none
! INTEGER ICO2 index of CO2 in gamnt none
! INTEGER ISO2 index of SO2 in gamnt none
! INTEGER IHNO3 index of HNO3 in gamnt none
! INTEGER IN2O index of N2O in gamnt none
! INTEGER IOPCI input RTP file I/O unit none
! INTEGER IOPCO output RTP file I/O unit none
! STRUCT HEAD RTP header structure various
! STRUCT HATT RTP header attributes none
! STRUCT PATT RTP profile attributes none
! LOGICAL LCO2PM CO2 profile in ppmv? none
!INPUT/OUTPUT PARAMETERS: none
!RETURN VALUES: none
!PARENT(S): sarta_rtp
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! IOPCI : input RTP file I/O unit ("profile channel")
!COMMON BLOCKS: none
!DESCRIPTION:
! Opens the input RTP file and reads the header info.
! Checks the header info.
!ALGORITHM REFERENCES: see DESCRIPTION
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
!------------ ----------------- ----------------------------------------
! 13 Feb 2001 Scott Hannon Created
! 23 Feb 2001 Scott Hannon Added PROF.gunit and GUCIN check
! 28 Feb 2001 Scott Hannon Add IOUN and read in chan freq file
! 14 Mar 2001 Scott Hannon Add HEAD, HATT, and PATT to call
! paramters. Removed open of output
! RTP; now done outside this routine.
! 13 Sep 2001 Scott Hannon Added AIRSLAY to ptype check
! 21 Nov 2001 Scott Hannon Remove CSARTA; change comment string
! to use VSARTA, VSCOEF, & VCLOUD; add
! CJUNK2, CJUNK3, & COMMNT; CJUNK
! decreased from 80 to 40;
! 20 Sep 2002 Scott Hannon If exists, overwrite old "sarta" hattr
! 05 Aug 2003 Scott Hannon Correct FIN to CHAR*80 (not 70)
! 06 Feb 2004 Scott Hannon Add LRHOT & PTYPE to arguments and add
! associated code.
! 18 May 2005 Scott Hannon Add HNO3 based on from SO2 code
! 23 Jun 2005 Scott Hannon "trace" version for CO2,SO2,HNO3,N2O
! 23 Jan 2008 Scott Hannon Add LCO2PM to allow CO2 profile in ppmv
! 24 Oct 2008 Scott Hannon Minor update for rtpV201
! 12 May 2009 Scott Hannon Change VCLOUD to VTUNNG in "sarta" HATT
!END====================================================================
! =================================================================
       SUBROUTINE OPNRTP(PTYPE, NCHAN, FCHAN, LSTCHN, &
          INDCHN, IH2O, IO3, ICO, ICH4, ICO2, ISO2, IHNO3, IN2O, &
          HEAD, PATT, LCO2PM)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
! include 'rtpdefs.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! From "util.f"
! function LENNB = length of string excluding trailing blanks
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
!
! Output
       INTEGER,intent(out):: PTYPE ! profile type
       INTEGER,intent(out):: NCHAN ! number of channels
       REAL,intent(out)   :: FCHAN(MXCHAN) ! channel freqs
       INTEGER,intent(out):: LSTCHN(MXCHAN) ! channel ID numbers
       INTEGER,intent(out):: INDCHN(MXCHAN) ! indices of channels
       INTEGER,intent(out):: IH2O ! index of H2O in gamnt
       INTEGER,intent(out):: IO3 ! index of O3 in gamnt
       INTEGER,intent(out):: ICO ! index of CO in gamnt
       INTEGER,intent(out):: ICH4 ! index of CH4 in gamnt
       INTEGER,intent(out):: ICO2 ! index of CO2 in gamnt
       INTEGER,intent(out):: ISO2 ! index of SO2 in gamnt
       INTEGER,intent(out):: IHNO3 ! index of HNO3 in gamnt
       INTEGER,intent(out):: IN2O ! index of N2O in gamnt

!
! Structures (see "rtpdefs.f")
! RECORD /RTPHEAD/ HEAD ! header data
! RECORD /RTPATTR/ HATT(MAXNATTR) ! header attributes
! RECORD /RTPATTR/ PATT(MAXNATTR) ! profile attributes
       type( RTPHEAD ) :: HEAD ! header data
!       type( RTPATTR ) :: HATT(MAXNATTR) ! header attributes
       type( RTPPROF ) :: PATT(MAXNATTR) ! profile attributes
       LOGICAL,intent(out):: LCO2PM ! CO2 profile in ppmv?
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER IC
! INTEGER IG
       INTEGER NGASI ! number of gases in input file
       INTEGER GLISTI( MXGAS) ! list of gas IDs in input file
       INTEGER J
       INTEGER K
!       INTEGER LENNB ! for function LENNB
!       INTEGER NHATT ! counter for # of header attributes
! INTEGER NPATT ! counter for # of profile attributes
       INTEGER rtpopen ! function rtpopen
       INTEGER STATUS ! status of RTP file open
       CHARACTER*1 MODE ! mode for rtpopen: "c"=create, "r"=read
       CHARACTER*1 CRHOT ! LRHOT converted to character T or F
! CHARACTER*14 CUNITS ! string for gamnt units
       CHARACTER*40 CJUNK ! junk/work string
       CHARACTER*40 CJUNK2 ! another junk/work string
       CHARACTER*40 CJUNK3 ! yet another junk/work string
       CHARACTER*256 COMMNT ! comment string
!
! for N2BITS and BITS2N
       INTEGER*4 NUMBER
       LOGICAL LFLAGS(32)
!
       LOGICAL LNEED ! needed gas?
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE begins below
!***********************************************************************
!***********************************************************************
!
! -------------------
! Open RTP input file
! -------------------
       MODE='r'

!  call read_sarta_param(HEAD,PATT)

! STATUS=rtpopen(FIN, MODE, HEAD, HATT, PATT, IOPCI)
!cc
! print *, 'read open status = ', STATUS
!cc
! -------------------------
! Quick checks of input RTP
! -------------------------
       PTYPE= 1   ! 1 or 2

       IF (PTYPE .NE. LAYPRO .AND. PTYPE .NE. AIRSLAY) THEN
          WRITE(IOERR,1003)
 1003 FORMAT('Error! input RTP ptype must be LAYPRO or AIRSLAY')
          STOP
       ENDIF
! Note: if no RHO data will use (1-emis)/pi
!cc
! Removed 26 April 2001 by Scott Hannon since mlev may be less than MAXLAY+1
!cc

       NCHAN= nchan_sarta  !total number of channels to be simulated. 
       IF (NCHAN .LT. 1) THEN
          WRITE(IOERR,1007)
 1007 FORMAT('Error! input RTP has no channel info')
          STOP
       ENDIF

!
! -----------
! Check gases
! -----------
       LCO2PM=.FALSE.
       IH2O =-1
       ICO2 =-1
       IO3 =-1
       IN2O =-1
       ICO =-1
       ICH4 =-1
       ISO2 =-1
       IHNO3=-1
!
! Loop over gases
       NGASI=ngas

       DO I=1,NGASI
          GLISTI(I)=glist(I)
          LNEED=.FALSE.
!
! Determine indices of needed gases
! Note: will abort if a needed gas is not present
          IF (GLISTI(I) .EQ. 1) THEN
             IH2O=I
             LNEED=.TRUE.
          ENDIF
          IF (GLISTI(I) .EQ. 3) THEN
             IO3=I
             LNEED=.TRUE.
          ENDIF
          IF (GLISTI(I) .EQ. 5) THEN
             ICO=I
             LNEED=.TRUE.
          ENDIF
          IF (GLISTI(I) .EQ. 6) THEN
             ICH4=I
             LNEED=.TRUE.
          ENDIF
!
! Determine indices of trace gases
! Note: will use reference amount if a trace gas is not present
! Exception: CO2 will use CO2PPM
          IF (GLISTI(I) .EQ. 2) THEN
             ICO2=I
             LNEED=.TRUE.
          ENDIF
          IF (GLISTI(I) .EQ. 4) THEN
             IN2O=I
             LNEED=.TRUE.
          ENDIF
          IF (GLISTI(I) .EQ. 9) THEN
             ISO2=I
             LNEED=.TRUE.
          ENDIF
          IF (GLISTI(I) .EQ. 12) THEN
             IHNO3=I
             LNEED=.TRUE.
          ENDIF
!
! Check gas units
          IF (LNEED) THEN
             IF (gunit(I) .NE. GUCIN) THEN
                IF (I .EQ. ICO2 .AND. gunit(I) .EQ. 10) THEN
                   LCO2PM=.TRUE.
      print *,'CO2 profile in ppmv'
                ELSE
                   WRITE(IOERR,1020) GUCIN, I, gunit(I)
 1020 FORMAT('ERROR &! Wrong gas units code number. ', &
                   'Need ',I3,' but gunit(',I2,')=',I3)
                   STOP
                ENDIF
             ENDIF
          ENDIF
!
       ENDDO
!
! Abort if a needed gas is not present
       IF (IH2O .LT. 1) THEN
          WRITE(IOERR,1030) 1, 'H2O '
          STOP
       ELSEIF (IO3 .LT. 1) THEN
          WRITE(IOERR,1030) 3, 'O3  '
          STOP
       ELSEIF (ICO .LT. 1) THEN
          WRITE(IOERR,1030) 5, 'CO  '
          STOP
       ELSEIF (ICH4 .LT. 1) THEN
          WRITE(IOERR,1030) 6, 'CH4 '
          STOP
       ENDIF
 1030 FORMAT('Error! input files does not contain gas ',I2,' = ',A4)

! Print a warning if a trace gas is not present
       IF (ICO2 .LT. 1) THEN
          WRITE(IOERR,1035) 2, 'CO2 '
       ENDIF
       IF (IN2O .LT. 1) THEN
          WRITE(IOERR,1035) 4, 'N2O '
       ENDIF
       IF (ISO2 .LT. 1) THEN
          WRITE(IOERR,1035) 9, 'SO2 '
       ENDIF
       IF (IHNO3 .LT. 1) THEN
          WRITE(IOERR,1035) 12, 'HNO3'
       ENDIF
 1035 FORMAT('Warning! input files does not contain gas ',I2,' = ',A4)

! -------------------
! Create channel list
! -------------------
! Initialize channel index list

       DO I=1,MXCHAN  !MXCHAN is 2834
          INDCHN(I)=0
       ENDDO


       K=0 ! initialize counter
       DO I=1,NCHAN
          J= ichan_sarta(I)  !HEAD%ichan(I) ! channel ID
          FCHAN(I)= undefined !HEAD%vchan(I) ! channel freq (or junk if unfilled)
!
          IF ((J .LT. 1) .OR. (J .GT. MXCHAN)) THEN
             WRITE(IOERR,1042) MXCHAN, J
 1042 FORMAT('Error &! Channel number is out of range.',/, &
             'Range is 1 to ',I4,', but input RTP has ',I7)
             STOP
          ENDIF
!
          IF (INDCHN(J) .EQ. 0) THEN
! Not a repeat
             K=K + 1 ! increment counter (should be same as I)
             LSTCHN(K)=J   !Here is list channel
             INDCHN(J)=K   !this is channel idice
!
          ELSE
             WRITE(IOERR,1044) J
 1044 FORMAT('ERROR &! input RTP has repeat of ', &
             'channel ',I4)
             STOP
          ENDIF
       ENDDO

       RETURN
       END SUBROUTINE OPNRTP

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! QIKEXP
!
!F77====================================================================
!ROUTINE NAME:
! QIKEXP (real function)
!ABSTRACT:
! Quick approximate calculation of e^x.
!CALL PROTOCOL
! QIKEXP( XVAL )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL XVAL exponent none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL fun QIKEXP exp(XVAL) none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! CALT1
! CALT2
! CALT3
! CALT4
! CALT5
! CALT6
! CALT7
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! none
!COMMON BLOCKS
! none
!DESCRIPTION:
! March 1998 version of the 100 layer AIRS Fast Transmittance
! Code by L.Strow/S.Hannon.
!
! ===================================================================
! Quick exponential calculation of e^x. The function looks at x and
! if it is small, it does the exponential calculation by using just
! the first few terms of the series expansion:
! exp(x) = sum i=0 to inf of { x^n/n! }.
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! None, but keep in mind it is only quicker than EXP if "x" is small.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 3 Feb 1997 Scott Hannon Created as a stand-alone function
! 18 Jul 1997 Scott Hannon Changed from 3 to 4 xval regions; same
! speed but more accurate
!END====================================================================
! =================================================================
       REAL FUNCTION QIKEXP( XVAL )
! =================================================================
!
! QuIcK EXPonential
!
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
       REAL XVAL
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       REAL AXVAL
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Absolute value of XVAL
       AXVAL=ABS(XVAL)
!
       IF (AXVAL .LT. 4.0E-03) THEN
! Use the first two series terms only:
          QIKEXP=1.0E+0 + XVAL
!
       ELSEIF (AXVAL .LT. 3.6E-02) THEN
! Use the first three series terms only:
          QIKEXP=1.0E+0 + XVAL + XVAL*XVAL*5.0E-1
!
       ELSEIF (AXVAL .LT. 1.2E-01) THEN
! Use the first four series terms only:
          QIKEXP=1.0E+0 + ( XVAL*( 6.0E+0 + (XVAL*(3.0E+0 + XVAL)) ) &
             /6.0E+0 )
!
       ELSE
          QIKEXP=EXP(XVAL)
       ENDIF
!
       RETURN
       END FUNCTION QIKEXP
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! RDCOEF
!
!F77====================================================================
!ROUTINE NAME:
! RDCOEF
!ABSTRACT:
! Read in the AIRS fast transmittance coefficients.
!CALL PROTOCOL
! RDCOEF ( AORB, IOUN, NCHAN, INDCHN, SETCHN,
! $ NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7,
! $ CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
! $ COEF1, COEF2, COEF3, COEF4, COEF5, COEF6, COEF7,
! $ NCHCO2, INDCO2, CLICO2, COFCO2,
! $ NCHN2O, INDN2O, CLIN2O, COFN2O,
! $ NCHSO2, INDSO2, CLISO2, COFSO2,
! $ NCHHNO, INDHNO, CLIHNO, COFHNO,
! $ NCHH2O, INDH2O, CLIH2O, WAZOP, WAVGOP, COFH2O,
! $ FREQ, COEFF, NCHNTE, CLISTN, COEFN, FX)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! CHAR*1 AORB specify A or B none
! INT arr INDCHN indices of channels none
! INTEGER IOUN I/O unit number none
! INTEGER NCHAN number of channels none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INT arr CLICO2 CO2 channel list none
! INT arr CLIHNO HNO3 channel list none
! INT arr CLIH2O H2O channel list none
! INT arr CLIN2O N2O channel list none
! INT arr CLISO2 SO2 channel list none
! INT arr CLIST1 set1 channel list none
! INT arr CLIST2 set2 channel list none
! INT arr CLIST3 set3 channel list none
! INT arr CLIST4 set4 channel list none
! INT arr CLIST5 set5 channel list none
! INT arr CLIST6 set6 channel list none
! INT arr CLIST7 set7 channel list none
! INT arr CLISTN non-LTE channel list none
! REAL arr COEF1 set1 fast trans coefs various
! REAL arr COEF2 set2 fast trans coefs various
! REAL arr COEF3 set3 fast trans coefs various
! REAL arr COEF4 set4 fast trans coefs various
! REAL arr COEF5 set5 fast trans coefs various
! REAL arr COEF6 set6 fast trans coefs various
! REAL arr COEF7 set7 fast trans coefs various
! REAL arr COEFF thermal "F" factor coefs various
! REAL arr COEFN non-LTE coefficients various
! REAL arr COFCO2 CO2 perturbation coefs various
! REAL arr COFSO2 SO2 perturbation coefs various
! REAL arr COFHNO HNO3 perturbation coefs various
! REAL arr COFN2O N2O perturbation coefs various
! REAL arr COFH2O OPTRAN H2O coefs various
! REAL arr FREQ channel freqs cm-1
! REAL arr FX fixed gases adjustment none
! INT arr INDCO2 CO2 pert channel indices none
! INT arr INDSO2 SO2 pert channel indices none
! INT arr INDHNO HNO3 pert channel indices none
! INT arr INDN2O N2O pert channel indices none
! INT arr INDH2O OPTRAN H2O channel indices none
! INTEGER NCHN1 set1 number of channels none
! INTEGER NCHN2 set2 number of channels none
! INTEGER NCHN3 set3 number of channels none
! INTEGER NCHN4 set4 number of channels none
! INTEGER NCHN5 set5 number of channels none
! INTEGER NCHN6 set6 number of channels none
! INTEGER NCHN7 set7 number of channels none
! INTEGER NCHCO2 number of CO2 pert chans none
! INTEGER NCHH2O number of OPTRAN H2O chans none
! INTEGER NCHHNO number of HNO3 pert chans none
! INTEGER NCHN2O number of N2O pert chans none
! INTEGER NCHNTE non-LTE number of channels none
! INTEGER NCHSO2 number of SO2 pert chans none
! REAL arr WAZOP OPTRAN water grid kiloMoles/cm^2
! REAL arr WAVGOP OPTRAN water pred averges various
! INT arr SETCHN set# (1-7) chan belongs to none (integer, 1 - 7)
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
! unit IOUN : input file, binary FORTRAN data file. The file is
! opened, read, and closed. This is done 10 times, once per
! each of the 7 coef sets, and once each for the variable CO2,
! OPTRAN water, and thermal F factor coefs.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2009 version of the 100 layer AIRS Fast Transmittance
! Code by L.Strow/S.Hannon.
!
! Seven sets of binary data files containing the main fast
! transmittance coefficients are opened and read one channel at
! a time. The seven sets of coefs are each stored in their own
! arrays. Next, preturbation coefficients for four trace gases,
! (CO2, SO2, HNO3, & N2O) are read in from four binary files.
! Next, OPTRAN water fast trans coefs for some channels are read
! in from a binary file file. The header of the OPTRAN file
! specifies 300 OPTRAN water levels, and also the mean value of
! 4 predictor terms for each of the levels. Next, comes the
! read of the binary file with the reflected downwelling thermal
! radiance "F factor" coefficients. Next is a read of the "FX"
! fixed gases adjustment term from an ASCII text file. Lastly
! comes the read of the non-LTE coefficients from a binary file.
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- -------------------------------------------
! 01 Dec 1994 Scott Hannon Created
! 21 Dec 1994 Scott Hannon Fixed error with IOPF (now assigned)
! 05 Feb 1997 Scott Hannon Re-wrote for FWO+FOW+FMW+FCOW.
! 28 Aug 1997 Scott Hannon Re-wrote for sets 1 - 7 and thermal
! 30 Sep 1997 Scott Hannon Added COFCO2 and INDCO2
! 27 Feb 1998 Scott Hannon Added COFH2O, INDH2O, WAZOP, & WAVGOP
! 17 Aug 2000 Scott Hannon Add FX
! 12 Feb 2001 Scott Hannon hardcoded filenames instead of prompts
! 18 May 2005 Scott Hannon Add HNO3 based on SO2 code
! 28 Jun 2005 Scott Hannon "trace" version for CO2,SO2,HNO3,N2O
! 13 Oct 2005 Scott Hannon Add non-LTE variables
! 09 May 2008 Scott Hannon Remove LABOVE; add NCHCO2, NCHN2O,
! NCHSO2, NCHHNO, NCHH2O; reorder call;
! comment out LACHAN
! 03 Aug 2009 Scott Hannon Add CLI<gas> for CO2,N2O,HNO,H2O,SO2
!END====================================================================
! =================================================================
       SUBROUTINE RDCOEF ( AORB, IOUN, NCHAN, INDCHN, SETCHN, &
           NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7, &
          CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
           COEF1, COEF2, COEF3, COEF4, COEF5, COEF6, COEF7, &
          NCHCO2, INDCO2, CLICO2, COFCO2, &
          NCHN2O, INDN2O, CLIN2O, COFN2O, &
          NCHSO2, INDSO2, CLISO2, COFSO2, &
          NCHHNO, INDHNO, CLIHNO, COFHNO, &
          NCHH2O, INDH2O, CLIH2O, WAZOP, WAVGOP, COFH2O, &
            FREQ, COEFF, NCHNTE, CLISTN, COEFN, FX)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       CHARACTER*1 AORB
       INTEGER IOUN
       INTEGER NCHAN
       INTEGER INDCHN(MXCHAN)
!
! Output
       INTEGER SETCHN(MXCHAN)
       INTEGER NCHN1
       INTEGER NCHN2
       INTEGER NCHN3
       INTEGER NCHN4
       INTEGER NCHN5
       INTEGER NCHN6
       INTEGER NCHN7
       INTEGER CLIST1(MXCHN1)
       INTEGER CLIST2(MXCHN2)
       INTEGER CLIST3(MXCHN3)
       INTEGER CLIST4(MXCHN4)
       INTEGER CLIST5(MXCHN5)
       INTEGER CLIST6(MXCHN6)
       INTEGER CLIST7(MXCHN7)
       REAL COEF1(N1COEF,MAXLAY,MXCHN1)
       REAL COEF2(N2COEF,MAXLAY,MXCHN2)
       REAL COEF3(N3COEF,MAXLAY,MXCHN3)
       REAL COEF4(N4COEF,MAXLAY,MXCHN4)
       REAL COEF5(N5COEF,MAXLAY,MXCHN5)
       REAL COEF6(N6COEF,MAXLAY,MXCHN6)
       REAL COEF7(N7COEF,MAXLAY,MXCHN7)
       INTEGER NCHCO2
       INTEGER INDCO2(MXCHAN)
       INTEGER CLICO2(MXCHNC)
       REAL COFCO2( NCO2,MAXLAY,MXCHNC)
       INTEGER NCHN2O
       INTEGER INDN2O(MXCHAN)
       INTEGER CLIN2O(MXCHNN)
       REAL COFN2O( NN2O,MAXLAY,MXCHNN)
       INTEGER NCHSO2
       INTEGER INDSO2(MXCHAN)
       INTEGER CLISO2(MXCHNS)
       REAL COFSO2( NSO2,MAXLAY,MXCHNS)
       INTEGER NCHHNO
       INTEGER INDHNO(MXCHAN)
       INTEGER CLIHNO(MXCHNH)
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH)
       INTEGER NCHH2O
       INTEGER INDH2O(MXCHAN)
       INTEGER CLIH2O(MXCHNW)
       REAL WAZOP(MXOWLY)
       REAL WAVGOP(NOWAVG,MXOWLY)
       REAL COFH2O( NH2O,MXOWLY,MXCHNW)
       REAL FREQ(MXCHAN)
       REAL COEFF(NFCOEF,MXCHAN)
       INTEGER NCHNTE
       INTEGER CLISTN(MXCNTE)
       REAL COEFN(NNCOEF,MXCNTE)
       REAL FX(MAXLAY)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       CHARACTER*80 FNCOF1 ! coef set1
       CHARACTER*80 FNCOF2 ! coef set2
       CHARACTER*80 FNCOF3 ! coef set3
       CHARACTER*80 FNCOF4 ! coef set4
       CHARACTER*80 FNCOF5 ! coef set5
       CHARACTER*80 FNCOF6 ! coef set6
       CHARACTER*80 FNCOF7 ! coef set7
       CHARACTER*80 FNCO2 ! coef CO2
       CHARACTER*80 FNSO2 ! coef SO2
       CHARACTER*80 FNHNO3 ! coef HNO3
       CHARACTER*80 FNN2O ! coef N2O
       CHARACTER*80 FNOPTR ! coef optran
       CHARACTER*80 FNTHER ! coef therm
       CHARACTER*80 FNCOFN ! coef non-LTE
       CHARACTER*80 CLINE
       REAL FRQCHN
       REAL FCHAN(NFCOEF)
       REAL RJUNK
       INTEGER I
       INTEGER IC
       INTEGER ICHAN
       INTEGER IERR
       INTEGER IL
       INTEGER J
       INTEGER ICOUNT
! INTEGER LACHAN
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Coefficient filenames
       IF (AORB .EQ. 'A') THEN
           FNCOF1=FACOF1
           FNCOF2=FACOF2
           FNCOF3=FACOF3
           FNCOF4=FACOF4
           FNCOF5=FACOF5
           FNCOF6=FACOF6
           FNCOF7=FACOF7
           FNCO2 =FACO2
           FNSO2 =FASO2
           FNHNO3=FAHNO3
           FNN2O =FAN2O
           FNOPTR=FAOPTR
           FNTHER=FATHER
           FNCOFN=FACOFN
       ELSE
           FNCOF1=FBCOF1
           FNCOF2=FBCOF2
           FNCOF3=FBCOF3
           FNCOF4=FBCOF4
           FNCOF5=FBCOF5
           FNCOF6=FBCOF6
           FNCOF7=FBCOF7
           FNCO2 =FBCO2
           FNSO2 =FBSO2
           FNHNO3=FBHNO3
           FNN2O =FBN2O
           FNOPTR=FBOPTR
           FNTHER=FBTHER
           FNCOFN=FBCOFN
       ENDIF
! Initialize "set"-independent index arrays
       DO I=1,MXCHAN
! Trace gases
          INDCO2(I)=0
          INDSO2(I)=0
          INDHNO(I)=0
          INDN2O(I)=0
! OPTRAN water
          INDH2O(I)=0
       ENDDO
!
! ----------
! Read set 1
! ----------
       OPEN(UNIT=IOUN,FILE=FNCOF1,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCOF1
 1020 FORMAT('Error ',I5,' openning file:',/,A80)
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHN1
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COEF1(IC,IL,J),IC=1,N1COEF), &
             IL=1,MAXLAY)


!
          SETCHN(ICHAN)=1
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIST1(J)=ICHAN
             FREQ( INDCHN(ICHAN) )=FRQCHN
             J=J + 1
          ENDIF
       ENDDO
       NCHN1=J - 1
!
       CLOSE(IOUN)

!
!
! ----------
! Read set 2
! ----------
       OPEN(UNIT=IOUN,FILE=FNCOF2,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCOF2
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHN2
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COEF2(IC,IL,J),IC=1,N2COEF), &
             IL=1,MAXLAY)
!

          SETCHN(ICHAN)=2
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIST2(J)=ICHAN
             FREQ( INDCHN(ICHAN) )=FRQCHN
             J=J + 1
          ENDIF
       ENDDO
       NCHN2=J - 1
!
       CLOSE(IOUN)

!
!
! ----------
! Read set 3
! ----------
       OPEN(UNIT=IOUN,FILE=FNCOF3,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCOF3
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHN3
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COEF3(IC,IL,J),IC=1,N3COEF), &
             IL=1,MAXLAY)
!
          SETCHN(ICHAN)=3
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIST3(J)=ICHAN
             FREQ( INDCHN(ICHAN) )=FRQCHN
             J=J + 1
          ENDIF
       ENDDO
       NCHN3=J - 1
!
       CLOSE(IOUN)
!
!
! ----------
! Read set 4
! ----------
       OPEN(UNIT=IOUN,FILE=FNCOF4,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCOF4
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHN4
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COEF4(IC,IL,J),IC=1,N4COEF), &
             IL=1,MAXLAY)
!
          SETCHN(ICHAN)=4
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIST4(J)=ICHAN
             FREQ( INDCHN(ICHAN) )=FRQCHN
             J=J + 1
          ENDIF
       ENDDO
       NCHN4=J - 1
!
       CLOSE(IOUN)
!
!
! ----------
! Read set 5
! ----------
       OPEN(UNIT=IOUN,FILE=FNCOF5,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCOF5
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHN5
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COEF5(IC,IL,J),IC=1,N5COEF), &
             IL=1,MAXLAY)
!
          SETCHN(ICHAN)=5
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIST5(J)=ICHAN
             FREQ( INDCHN(ICHAN) )=FRQCHN
             J=J + 1
          ENDIF
       ENDDO
       NCHN5=J - 1
!
       CLOSE(IOUN)
!
!
! ----------
! Read set 6
! ----------
       OPEN(UNIT=IOUN,FILE=FNCOF6,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCOF6
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHN6
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COEF6(IC,IL,J),IC=1,N6COEF), &
             IL=1,MAXLAY)
!
          SETCHN(ICHAN)=6
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIST6(J)=ICHAN
             FREQ( INDCHN(ICHAN) )=FRQCHN
             J=J + 1
          ENDIF
       ENDDO
       NCHN6=J - 1
!
       CLOSE(IOUN)
!
!
! ----------
! Read set 7
! ----------
       OPEN(UNIT=IOUN,FILE=FNCOF7,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCOF7
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHN7
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COEF7(IC,IL,J),IC=1,N7COEF), &
             IL=1,MAXLAY)
!
          SETCHN(ICHAN)=7
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIST7(J)=ICHAN
             FREQ( INDCHN(ICHAN) )=FRQCHN
             J=J + 1
          ENDIF
       ENDDO
       NCHN7=J - 1
!
       CLOSE(IOUN)
!
!
! ---------------------------
! Read CO2 perturbation coefs
! ---------------------------
       OPEN(UNIT=IOUN,FILE=FNCO2,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCO2
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHNC
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COFCO2(IC,IL,J),IC=1,NCO2), &
             IL=1,MAXLAY)
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLICO2(J)=ICHAN
             INDCO2(ICHAN)=J
             J=J + 1
          ENDIF
       ENDDO
       NCHCO2=J - 1
!
       CLOSE(IOUN)
!
!
! ---------------------------
! Read SO2 perturbation coefs
! ---------------------------
       OPEN(UNIT=IOUN,FILE=FNSO2,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNSO2
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHNS
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COFSO2(IC,IL,J),IC=1,NSO2), &
             IL=1,MAXLAY)
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLISO2(J)=ICHAN
             INDSO2(ICHAN)=J
             J=J + 1
          ENDIF
       ENDDO
       NCHSO2=J - 1
!
       CLOSE(IOUN)
!
!
! ---------------------------
! Read HNO3 perturbation coefs
! ---------------------------
       OPEN(UNIT=IOUN,FILE=FNHNO3,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNHNO3
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHNH
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COFHNO(IC,IL,J),IC=1,NHNO3), &
             IL=1,MAXLAY)
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIHNO(J)=ICHAN
             INDHNO(ICHAN)=J
             J=J + 1
          ENDIF
       ENDDO
       NCHHNO=J - 1
!
       CLOSE(IOUN)
!
! ---------------------------
! Read N2O perturbation coefs
! ---------------------------
       OPEN(UNIT=IOUN,FILE=FNN2O,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNN2O
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCHNN
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COFN2O(IC,IL,J),IC=1,NN2O), &
             IL=1,MAXLAY)
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIN2O(J)=ICHAN
             INDN2O(ICHAN)=J
             J=J + 1
          ENDIF
       ENDDO
       NCHN2O=J - 1
!
       CLOSE(IOUN)
!
! ---------------------
! Read OPTRAN H2O coefs
! ---------------------
       OPEN(UNIT=IOUN,FILE=FNOPTR,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNOPTR
          STOP
       ENDIF
!
       READ(IOUN) (WAZOP(IL),IL=1,MXOWLY)
       DO IC=1,NOWAVG
! Read the header section
          READ(IOUN) (WAVGOP(IC,IL),IL=1,MXOWLY)
       ENDDO
!
       J=1
       DO I=1,MXCHNW
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, ((COFH2O(IC,IL,J),IC=1,NH2O), &
             IL=1,MXOWLY)
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLIH2O(J)=ICHAN
             INDH2O(ICHAN)=J
             J=J + 1
          ENDIF
       ENDDO
       NCHH2O=J - 1
!
       CLOSE(IOUN)
!
!
! -----------------------------------------------
! Read the downward thermal F factor coefficients
! -----------------------------------------------
       OPEN(UNIT=IOUN,FILE=FNTHER,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNTHER
          STOP
       ENDIF
!
       DO I=1,MXCHAN
! Read data for this frequency/channel
!cc changed 18 May 2005
! READ(IOUN) ICHAN, FRQCHN, LACHAN, (FCHAN(IC),IC=1,NFCOEF)
          READ(IOUN) ICHAN, FRQCHN, (FCHAN(IC),IC=1,NFCOEF)
! LACHAN=-1 ! assign dummy value
!cc
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             DO IC=1,NFCOEF
                COEFF(IC,INDCHN(ICHAN))=FCHAN(IC)
             ENDDO
          ENDIF
       ENDDO
!
       CLOSE(IOUN)
!
!
! -------
! Read FX
! -------
       OPEN(UNIT=IOUN,FILE=FNFX,FORM='FORMATTED',STATUS='OLD', &
             CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNFX
          STOP
       ENDIF
!
! Read the file
       ICOUNT=0
 10 READ(IOUN,9000,END=99) CLINE
 9000 FORMAT(A80)
       IF (CLINE(1:1) .NE. '!') THEN
! Note: fx file format is: layer_number fx_value
          READ(CLINE,*) IC, RJUNK
          ICOUNT=ICOUNT + 1
          FX(IC)=RJUNK
       ENDIF
       GOTO 10
!
 99 CLOSE(IOUN)
!
       IF (ICOUNT .NE. MAXLAY) THEN
          WRITE(6,1047) MAXLAY, ICOUNT
 1047 FORMAT('Error &! Unexpected number of layers in fx file.',/, &
          'Expected fx to have ',I4,' layers, but found ',I4)
       ENDIF
!
! ------------
! Read non-LTE
! ------------
       OPEN(UNIT=IOUN,FILE=FNCOFN,FORM='UNFORMATTED',STATUS='OLD', &
            CONVERT='BIG_ENDIAN',IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNCOFN
          STOP
       ENDIF
!
       J=1
       DO I=1,MXCNTE
! Read data for this frequency/channel
          READ(IOUN) ICHAN, FRQCHN, (COEFN(IC,J),IC=1,NNCOEF)
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             CLISTN(J)=ICHAN
             J=J + 1
          ENDIF
       ENDDO
       NCHNTE=J - 1
!
       CLOSE(IOUN)
! ---------------------------------------------
! Make sure all channels on the list were found
! ---------------------------------------------
       ICOUNT=NCHN1 + NCHN2 + NCHN3 + NCHN4 + NCHN5 + NCHN6 + NCHN7
       IF (ICOUNT .NE. NCHAN) THEN
          WRITE(6,1050) NCHAN, ICOUNT
 1050 FORMAT('Error &! Unexpected number of channels found.',/, &
          'The channel list had ',I4,' channels, but found ',I4)
       ENDIF
!
! ----------------------------
! Show summary of channel sets
! ----------------------------
!cc
! WRITE(6,1060) 1, NCHN1
! 1060 FORMAT('Number of channels for set',I1,' = ',I4)
! WRITE(6,1060) 2, NCHN2
! WRITE(6,1060) 3, NCHN3
! WRITE(6,1060) 4, NCHN4
! WRITE(6,1060) 5, NCHN5
! WRITE(6,1060) 6, NCHN6
! WRITE(6,1060) 7, NCHN7
!cc
!
       RETURN
       END SUBROUTINE RDCOEF
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! RDDYAB
!
!F77====================================================================
!ROUTINE NAME:
! RDDYAB
!ABSTRACT:
! Read in the A/B weight table
!CALL PROTOCOL:
! RDDYAB ( IOUN, ETIME, FREQAB, WIDTH, INDMOD, IABWGT, DELYAB )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER IOUN I/O unit number
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! DOUB arr ETIME A/B weight epoch start time TAI1993
! REAL arr FREQAB channel frequency cm^-1
! REAL arr WIDTH channel full width half max cm^-1
! INT arr INDMOD module index none
! INT arr IABWGT A/B weight index none
! REAL arr DELYAB delta Yoffset for A,B,AB um
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! SARTA
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
! unit IOUN : input file, ASCII text file
!COMMON BLOCKS
! none
!DESCRIPTION:
! October 2011 version of the 100 layer AIRS Fast Transmittance Code
! by L.Strow/S.Hannon.
!
! Read an ASCII text file containing delta Yoffset and A/B weight
! information.
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- -----------------------------------------
! 11 Oct 2011 Scott Hannon Created
! 21 Oct 2011 Scott Hannon bug fix: ETIME needs to be double
!END====================================================================
! =================================================================
       SUBROUTINE RDDYAB(IOUN, ETIME, FREQAB, WIDTH, INDMOD, IABWGT, &
          DELYAB)
! =================================================================
!
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
!
! Input
       INTEGER IOUN ! I/O unit number
! Output
       DOUBLE PRECISION ETIME(MAXABW) ! A/B weight epoch start time
       REAL FREQAB(MXCHAN) ! channel frequency
       REAL WIDTH(MXCHAN) ! channel width
       INTEGER INDMOD(MXCHAN) ! module index
       INTEGER IABWGT(MXCHAN,MAXABW) ! A/B weight index (for DELYAB)
       REAL DELYAB(MXCHAN,3) ! Delta Yoffset for A,B,AB
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I ! generic looping variable
       INTEGER J ! generic looping variable
       INTEGER IERR
       INTEGER IJUNK
       INTEGER NEPOCH
       CHARACTER*80 CLINE
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE follows...
!***********************************************************************
!***********************************************************************
!
! ---------------------
! Open the profile file
! ---------------------
       OPEN(UNIT=IOUN,FILE=ABFILE,STATUS='OLD',FORM='FORMATTED', &
          IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(IOERR,1010) IERR, ABFILE
 1010 FORMAT('Error ',I5,' openning A/B weight table:', &
             /,A80)
          STOP
       ENDIF
!
! ----------------------------------------
! Skip any comments at the top of the file
! ----------------------------------------
 10 READ(IOUN,9000) CLINE
 9000 FORMAT(A80)
       IF (CLINE(1:1) .EQ. '!') THEN
          GOTO 10
       ELSE
          BACKSPACE(IOUN)
       ENDIF
!
! ------------------------------------
! Read the number of A/B weight epochs
! ------------------------------------
       READ(IOUN,*) NEPOCH
       IF (NEPOCH .NE. MAXABW) THEN
          WRITE(IOERR,1020) NEPOCH, MAXABW
 1020 FORMAT('Error: A/B weight table NEPOCH=',I2, &
             '<> MAXABW=',I2)
       ENDIF
!
! --------------------------
! Read the epoch start times
! --------------------------
       DO I=1,NEPOCH
          READ(IOUN,*) ETIME(I)
       ENDDO
!
! -----------------
! Read in the table
! -----------------
       DO I=1,MXCHAN
          READ(IOUN,*) IJUNK, INDMOD(I), FREQAB(I), WIDTH(I), &
             DELYAB(I,1), DELYAB(I,2), DELYAB(I,3), &
             (IABWGT(I,J),J=1,NEPOCH)
       ENDDO
!
!cc uncomment for testing
! DO I=1,MXCHAN
! print *, I, INDMOD(I), FREQAB(I), WIDTH(I),
! $ DELYAB(I,1), DELYAB(I,2), DELYAB(I,3),
! $ (IABWGT(I,J),J=1,NEPOCH)
! ENDDO
!cc
! ----------------------
! Close the profile file
! ----------------------
       CLOSE(IOUN)
!
       RETURN
       END SUBROUTINE RDDYAB
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
! rdinfo processes command line arguments
!
! sarta fin=input.rtp fout=output.rtp listp=1,2,3
!
!
! to compile
! Absoft/Linux: f77 -N109 -o klayers $(SRC) -lU77
! SGI Irix: no special compiler options are needed
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! RDINFO_sarta
!
!F77====================================================================
!ROUTINE NAME: RDINFO
!ABSTRACT:
! Get info about the sarta run: the names of input & output
! files, the channel list, and list of profile numbers.
!CALL PROTOCOL:
! RDINFO(FIN, FOUT, LRHOT, NWANTP, LISTP)
!INPUT PARAMETERS:
! none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! CHAR*80 FIN input filename none
! CHAR*80 FOUT output filename none
! LOGICAL LRHOT force RHO for refl thermal? none
! INTEGER NWANTP Number of desired profiles none
! INT arr LISTP List of desired prof nums none
!INPUT/OUTPUT PARAMETERS: none
!RETURN VALUES: none
!PARENT(S): SARTA_rtp
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! none
!COMMON BLOCKS:
! none
!DESCRIPTION:
! Gets various info about sarta run.
!
! Each command line argument is of the form <variable>=<value>
! Each <variable>=<value> string must be 80 char or less. It is
! not necessary to enclose values in quotes unless they contain
! blanks. The recognized command-line variables are:
!
! fin : name of input file
!
! fout : name of output file
!
! lrhot : force reflected thermal rho?; logical. true/false of T/F
! If true, the refl therm will use rho=(1-emis)/pi rather than
! the rho (if any) from the input file.
!
! listp : list of desired profile numbers (all other profiles will
! be ignored). If "listp" is not specified, SARTA will process
! all profiles.
!
! The listp profile numbers may be specified either as a
! sequence of integers separated by a comma, or alternately as
! a quoted string containing integers separated by a blank space.
! Examples:
! listp=1,2,3,4,5
! listp='1 2 3 4 5'
! Due to the 80 char limit, the maximum number of entries
! in listp is limited. (Eg 15 four digit numbers, or
! 25 two digit numbers. MAXPRO is the hardcoded limit.)
!
!ALGORITHM REFERENCES: see DESCRIPTION
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
!------------ ----------------- ----------------------------------------
! 13 Feb 2001 H.Motteler/S.Hannon Re-write of KLAYERS version
! 28 Nov 2001 Scott Hannon Remove command-line argument "nwantp"
! 5 Dec 2001 Scott Hannon Remove unused local var LENNB
! 05 Aug 2003 Scott Hannon Correct FIN & FOUT to CHAR*80 (not 70)
! 06 Feb 2004 Scott Hannon Add LRHOT argument and associated code
!END====================================================================
! =================================================================
       SUBROUTINE RDINFO(FIN, FOUT, LRHOT, NWANTP, LISTP)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! From "util.f"
! subroutine UPCASE = converts a string to upper case
! function STR2BO = converts true/false string to boolean (LOGICAL)
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input:
! none
!
! Output:
       CHARACTER*80 FIN
       CHARACTER*80 FOUT
       LOGICAL LRHOT
       INTEGER NWANTP
       INTEGER LISTP(MAXPRO)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER IARGC
       INTEGER IP
       INTEGER IPJUNK(MXGAS+1) ! junk gas id work array
       INTEGER J
       INTEGER K
       INTEGER NARGS ! number of arguments
       INTEGER SORTED ! flag for sorting
       CHARACTER*80 BUF
       CHARACTER*80 VAL
       CHARACTER*80 VAR
       LOGICAL LLISTP
!       LOGICAL STR2BO
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE begins below
!***********************************************************************
!***********************************************************************

! print*,'MSG RDINFO: start '
! ------------
! Set defaults
! ------------
       FIN='sarta_in.rtp' ! input filename
       FOUT='sarta_out.rtp' ! output filename
       NWANTP=-1 ! do sarta for all profiles found in input file
       LRHOT=.FALSE. ! use input rho for reflected thermal
!
! -----------------------------------------------------------------
! Loop on program parameters
! --------------------------
! Determine the number of command-line arguments
! Toshi: For example, if you put some argument in excecutable
!>SARTA.x option1  option2 option3
! NARGS would be 3. 
       NARGS=IARGC()


!
! Loop over the command-line arguments
       LLISTP=.FALSE.
       DO I = 1, NARGS
!
! Pull out the ith argument
          CALL GETARG(I, BUF)
!
! Find the "=" character in the command-line argument string
          J=INDEX(BUF, '=')
!
          IF (J .NE. 0) THEN
!
! Name of variable
             VAR = BUF(1:J-1)
             CALL UPCASE(VAR)
!
! Specified value
             VAL = BUF(J+1:LEN(BUF))
!
! Big "IF" to set parameters
! ----------------------------
             IF (VAR(1:3) .EQ. 'FIN') THEN
                FIN=VAL
             ELSEIF (VAR(1:4) .EQ. 'FOUT') THEN
                FOUT=VAL
             ELSEIF (VAR(1:5) .EQ. 'LRHOT') THEN
                LRHOT=STR2BO(VAL)
             ELSEIF (VAR(1:5) .EQ. 'LISTP') THEN
                LLISTP=.TRUE.
!
! Read the indices of the desired profiles
                K=1
 10 IF (K .GT. MAXPRO) THEN
                   WRITE(6,1017)
 1017 FORMAT('ERROR &! bad LISTP, ', &
                   'either an unrecognized value or too many entries')
                   STOP
                ENDIF
                READ(VAL,*,END=19) (IPJUNK(IP),IP=1,K)
                K=K + 1 ! increment count of profiles
                GOTO 10 ! loop to next entry
 19 CONTINUE
                K=K - 1 ! number of profiles
                DO IP=1,K
                   LISTP(IP)=IPJUNK(IP)
                ENDDO
                NWANTP=K
!
             ELSE
                WRITE(6,1020) VAR(1:6)
 1020 FORMAT('Unknown command-line argument: ',A6)
                STOP
             ENDIF
          ENDIF
       ENDDO ! end of loop over command-line arguments
! -----------------------------------------------------------------
! -------------------------------------
! Sort prof numbers & check for repeats
! -------------------------------------
       IF (NWANTP .GT. 0) THEN
!
! Sort in ascending order
          SORTED=1 ! initialize flag for first pass
 30 IF (SORTED .EQ. 1) THEN
             SORTED=0 ! initialize flag for this loop
             DO K=1,NWANTP-1
                IF (LISTP(K) .GT. LISTP(K+1)) THEN
                   IP=LISTP(K)
                   LISTP(K)=LISTP(K+1)
                   LISTP(K+1)=IP
                   SORTED=1 ! set flag to indicate ordering was altered
                ENDIF
             ENDDO
             GOTO 30
          ENDIF
!
! Check for repeats
          DO K=1,NWANTP-1
             IF (LISTP(K) .EQ. LISTP(K+1)) THEN
                WRITE(6,1045) LISTP(K)
 1045 FORMAT('ERROR &! profile ',I2, &
                ' appears more than once in LISTP')
                STOP
             ENDIF
          ENDDO
       ENDIF
!
       RETURN
       END SUBROUTINE RDINFO
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! RDPROF version with SO2, HNO3, & N2O trace gases
!
!F77====================================================================
!ROUTINE NAME:
! RDPROF
!ABSTRACT:
! Read in an AIRS FTC formatted profile. Temperature, amounts, etc.
!CALL PROTOCOL:
! RDPROF ( IOUN, PFILE, PNAM, ALT, PRES, TEMP, FAMNT, WAMNT, OAMNT,
! CAMNT, MAMNT, SAMNT, HAMNT, NAMNT )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER IOUN I/O unit number none
! CHAR*80 PFILE filename for desired prof none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr ALT layer altitudes m
! REAL arr CAMNT carbon monoxide amount k.mol/cm2
! REAL arr DZ layer thickness m
! REAL arr FAMNT fixed gases amount k.mol/cm2
! REAL arr HAMNT HNO3 amount k.mol/cm2
! REAL arr MAMNT CH4 amount k.mol/cm2
! REAL arr NAMNT N2O amount k.mol/cm2
! REAL arr OAMNT O3 amount k.mol/cm2
! CHAR*40 PNAM profile name/comment none
! REAL arr SAMNT SO2 amount k.mol/cm2
! REAL arr TEMP temperature K
! REAL arr WAMNT water amount k.mol/cm2
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
! unit IOUN : input file, ASCII profile file
!COMMON BLOCKS
! none
!DESCRIPTION:
! March 1998 version of the 100 layer AIRS Fast Transmittance Code
! by L.Strow/S.Hannon.
!
! An ASCII file containing a profile (in the expected format) is
! read in, and the relevant temperature and amount profiles are
! passed back to the calling program.
!
! ===================================================================
! Reads profile data from a text file named PFILE. The data consists
! of a profile name/description PNAM, followed by rows of data for
! all 100 layers (in lowest to highest altitude order) consisting of
! columns:
! i, Z, dZ, P, T, F, W, O, C, M, S, H, N
! where:
! "i" is a layer number counter (value ignored)
! "Z" is the layer average altitude
! "dZ" is the layer thickness
! "P" is the layer slab average pressure PRES
! "T" is the layer slab average temperature TEMP (for "fixed")
! "F" is the "fixed" (CO2) gases amount
! "W" is the water (H2O) amount
! "O" is the ozone (O3) amount
! "C" is the carbon monoxide (CO) amount
! "M" is the methane (CH4) amount
! "S" is the sulfur dioxide (SO2) amount
! "H" is the nitric acid (HNO3) amount
! "N" is the nitrous oxide (N2O) amount
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- --------------------------------------
! Dec 1 1994 Scott Hannon Created; reading of FAMNT disabled
! Apr 10 1995 Scott Hannon New header comments; FAMNT enabled;
! added ALT
! Jun 23 1995 Scott Hannon Correct some comments
! Jul 3 1995 Scott Hannon Added parameter DZ for layer thickness
! 3 Feb 1997 Scott Hannon Add IOUN, CAMNT & MAMNT
! 18 May 2005 Scott Hannon Add HNO3 & N2O based on SO2 code
!END====================================================================
! =================================================================
       SUBROUTINE RDPROF (IOUN, PFILE, PNAM, ALT, DZ, PRES, TEMP, &
          FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT)
! =================================================================
!
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER IOUN
       CHARACTER*80 PFILE
!
! Output
! REAL arr ALT layer altitudes m
! REAL arr CAMNT carbon monoxide amount k.mol/cm2
! REAL arr DZ layer thickness m
! REAL arr FAMNT fixed gases amount k.mol/cm2
! REAL arr HAMNT HNO3 amount k.mol/cm2
! REAL arr MAMNT CH4 amount k.mol/cm2
! REAL arr NAMNT N2O amount k.mol/cm2
! REAL arr OAMNT O3 amount k.mol/cm2
! CHAR*40 PNAM profile name/comment none
! REAL arr SAMNT SO2 amount k.mol/cm2
! REAL arr TEMP temperature K
! REAL arr WAMNT water amount k.mol/cm2

       CHARACTER*40 PNAM
       REAL ALT(MAXLAY)   !ALT layer altitudes m
       REAL DZ(MAXLAY)    !layer thickness m
       REAL PRES(MAXLAY)  !
       REAL TEMP(MAXLAY)
       REAL FAMNT(MAXLAY) !fixed gases amount k.mol/cm2
       REAL WAMNT(MAXLAY) !water amount k.mol/cm2
       REAL OAMNT(MAXLAY) !
       REAL CAMNT(MAXLAY)
       REAL MAMNT(MAXLAY)
       REAL SAMNT(MAXLAY)
       REAL HAMNT(MAXLAY)
       REAL NAMNT(MAXLAY)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER IERR
       INTEGER IJUNK
       INTEGER L
       CHARACTER*80 CLINE
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE follows...
!***********************************************************************
!***********************************************************************
!
! ---------------------
! Open the profile file
! ---------------------
! print*,'MSG RDPROF: Opening and reading ',trim(PFILE)

       OPEN(UNIT=IOUN,FILE=PFILE,STATUS='OLD',FORM='FORMATTED', &
          IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1010) IERR, PFILE
 1010 FORMAT('Error ',I5,' openning profile file:',/,A80)
          STOP
       ENDIF
!
! ----------------------------------------
! Skip any comments at the top of the file
! ----------------------------------------
 10 READ(IOUN,9000) CLINE
 9000 FORMAT(A80)
       IF (CLINE(1:1) .EQ. '!') THEN
          GOTO 10
       ELSE
          BACKSPACE(IOUN)
       ENDIF
!
! -------------------------------
! Read the profile's name/comment
! -------------------------------
       READ(IOUN,9010) PNAM
 9010 FORMAT(A40)
!
! --------------------------------------------------
! Read in the temperature and amounts for each layer
! --------------------------------------------------
! Note: read the layers in reverse order.
       DO L=MAXLAY,1,-1
! Layer number, altitude, thickness, pressure, temperature,
! fixed, H2O, O3, CO, and CH4 amounts
          READ(IOUN,*) IJUNK, ALT(L), DZ(L), PRES(L), TEMP(L), &
             FAMNT(L), WAMNT(L), OAMNT(L), CAMNT(L), MAMNT(L), &
             SAMNT(L), HAMNT(L), NAMNT(L)


! toshi: FAMNT is k.mol/cm2, which is integrated path along thickness of each layer.
! Thus, for converting this to usual volume mixing ratio [ppmv], you need following eqs.
! print*, FAMNT(L) / (DZ(L)*100.)   / ( (PRES(L)*1.e5)/(8.314* TEMP(L))*1.e-9 ) * 1.e6    ! [ppmv]

       ENDDO


!
! ----------------------
! Close the profile file
! ----------------------
       CLOSE(IOUN)
!
       RETURN
       END SUBROUTINE RDPROF
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
! version for sarta
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! RDRTP version with trace gases
!
!F77====================================================================
!ROUTINE NAME: RDRTP
!ABSTRACT:
! Read a profile from a previously openned RTP file
!CALL PROTOCOL:
! RDRTP( LWANT, IPROF, IOPCI,
! IH2O, IO3, ICO, ICH4, ICO2, ISO2, IHNO3, IN2O,
! PTYPE, RALT,LCO2PM,NLAY, NEMIS, LAT, LON, SATANG, SATZEN,
! ZSAT, SUNANG, PSURF, TSURF, CO2PPM, FEMIS, EMIS, RHO,
! TEMP, WAMNT, OAMNT, CAMNT, MAMNT, FAMNT, SAMNT, HAMNT, NAMNT,
! ALT, PROF, ISTAT, YOFF, YMODP1, DOPPPM, OPHASE )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! LOGICAL LWANT do you want this profile? none
! INTEGER IPROF profile number none
! INTEGER IOPCI input RTP file I/O number none
! INTEGER IH2O index of H2O in gamnt none
! INTEGER IO3 index of O3 in gamnt none
! INTEGER ICO index of CO in gamnt none
! INTEGER ICH4 index of CH4 in gamnt none
! INTEGER ICO2 index of CO2 in gamnt none
! INTEGER ISO2 index of SO2 in gamnt none
! INTEGER IHNO3 index of HNO3 in gamnt none
! INTEGER IN2O index of N2O in gamnt none
! INTEGER PTYPE profile type code number none
! REAL arr RALT ref prof layer altitudes meters
! LOGICAL LCO2PM CO2 profile in ppmv? none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER NLAY number of used layers none
! INTEGER NEMIS number of used emis pts none
! REAL LAT latitude degrees
! REAL LON longitude degrees
! REAL SATANG satellite scan angle degrees
! REAL SATZEN satellite zenith angle degrees
! REAL ZSAT satellite altitude kilometer
! REAL SUNANG sun zenith angle degrees
! REAL PSURF surface pressure millibars
! REAL TSURF surface skin temperature Kelvin
! REAL CO2PPM mean trop-strat CO2 mix rat PPMV
! REAL arr FEMIS emis freq points cm^-1
! REAL arr EMIS emis points none (0 to 1)
! REAL arr RHO rho points none (0 to 1/pi)
! REAL arr TEMP layer temperature Kelvin
! REAL arr WAMNT layer Water vapor amount */cm^2
! REAL arr OAMNT layer Ozone amount */cm^2
! REAL arr CAMNT layer CO amount */cm^2
! REAL arr MAMNT layer Methane amount */cm^2
! REAL arr FAMNT layer CO2 amount */cm^2
! REAL arr SAMNT layer SO2 amount */cm^2
! REAL arr HAMNT layer HNO3 amount */cm^2
! REAL arr NAMNT layer N2O amount */cm^2
! REAL arr ALT layer average altitude meters
! STRUCT PROF RTP profile structure various
! INTEGER ISTAT I/O status none
! REAL YOFF overall Yoffset um
! REAL arr YMODP1 Yoffset of modules+1 um
! REAL DOPPPM Doppler shift in ppm ppm
! REAL OPHASE orbit phase degrees
! note: units "*/cm^2" can be either kilomoles/cm^2 or molecules/cm^2
!INPUT/OUTPUT PARAMETERS: none
!RETURN VALUES: none
!PARENT(S): SARTA
!ROUTINES CALLED: none
!FILES ACCESSED:
! Input RTP file withI/O number IOPCI
! unit IOERR: error messages
! unit IOINFO: info/warning messages
!COMMON BLOCKS: none
!DESCRIPTION:
! Reads a single profile from a previously openned RTP file.
! The routine expects to find the data specified in the header
! of the input RTP file.
!ALGORITHM REFERENCES: see DESCRIPTION
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
!------------ ----------------- ----------------------------------------
! 14 Feb 2001 Scott Hannon created based on klayers version
! 13 Sep 2001 Scott Hannon Added checks of PSURF & TSURF
! 31 Oct 2002 Scott Hannon Add output vars SATZEN and ZSAT
! 20 Dec 2004 Scott Hannon Add PTYPE to call; fix error in NLAY
! when PTYPE=AIRSLAY; add error trap
! for LAT
! 18 May 2005 Scott Hannon Add HNO3 based on SO2 code
! 23 Jun 2005 Scott Hannon "trace" version for CO2,SO2,HNO3,N2O
! 23 Jan 2008 Scott Hannon Add LCO2PM to allow CO2 profile in ppmv;
! fix bug in CO2PPM check
! 09 May 2008 Scott Hannon Add YOFF
! 24 Oct 2008 Scott Hannon Update for RTP v2.01; emis & rho now
! use the same freq points; set RHO
! to (1-e)/pi if input < 0
! 11 Oct 2011 Scott Hannon Add output arguments YMODP1, DOPPPM,
! and OPHASE
! 21 Nov 2011 Scott Hannon Set ymodp1(18) to yofdef
!END====================================================================
! =================================================================

!
! Now this routine is a coupler between G-SDSU-defined model parameters
! to the SARTA module parameters (Toshi Matsui OCT 2013).
!

       SUBROUTINE RDRTP(atmos1p, atmos_stag1p,surf1p, emiss_sarta, HEAD, PATT,  &
          IH2O, IO3, ICO, ICH4, ICO2, ISO2, IHNO3, IN2O, PTYPE, RALT, &
          LCO2PM, &
          NLAY, NEMIS, LAT, LON, SATANG, SATZEN, ZSAT, SUNANG, &
          PSURF, TSURF, CO2PPM, FEMIS, EMIS, RHO, &
          TEMP, WAMNT, OAMNT, CAMNT, MAMNT, FAMNT, SAMNT, HAMNT, NAMNT, &
          ALT, PROF, ISTAT, YOFF, YMODP1, DOPPPM, OPHASE )
! =================================================================
 use module_gas_profile
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
! include 'rtpdefs.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
       !REAL C4TOR4
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input parameters:
 type ( atmos_parameter )        ,intent(in) :: atmos1p(:)      !G-SDSU atmospheric layer profile (assumed length)
 type ( atmos_stagger_parameter ),intent(in) :: atmos_stag1p(:) !G-SDSU atmospheric level profile (assumed length)
 type ( surface_parameter )      ,intent(in) :: surf1p          !G-SDSU surface point 
 real,intent(in) :: emiss_sarta(anchor_bs_emiss)  !sepctral emissivity 
 type( RTPHEAD ),intent(inout) :: HEAD ! sarta header data
 type( RTPPROF ),intent(inout) :: PATT ! profile attributes

!       LOGICAL LWANT ! do we want this profile?
       INTEGER IPROF ! number of current profile
       INTEGER IH2O ! index of H2O in gamnt
       INTEGER IO3 ! index of O3 in gamnt
       INTEGER ICO ! index of CO in gamnt
       INTEGER ICH4 ! index of CH4 in gamnt
       INTEGER ICO2 ! index of CO2 in gamnt
       INTEGER ISO2 ! index of SO2 in gamnt
       INTEGER IHNO3 ! index of HNO3 in gamnt
       INTEGER IN2O ! index of N2O in gamnt
       INTEGER PTYPE ! profile type code number
       REAL RALT(MAXLAY) ! ref prof layer average altitudes
       LOGICAL LCO2PM ! CO2 profile in ppmv?
! Output parameters:
       INTEGER NLAY
       INTEGER NEMIS
       REAL LAT
       REAL LON
       REAL SATANG
       REAL SATZEN
       REAL ZSAT
       REAL SUNANG     !solar zenith angle (at 0 altitude)
       REAL PSURF      !surface pressure [mb]
       REAL TSURF      !surface skintemperature [K]
       REAL CO2PPM
       REAL FEMIS(MXEMIS)
       REAL EMIS(MXEMIS)
       REAL RHO(MXEMIS)
       REAL GAS(MAXLAY)
       REAL TEMP(MAXLAY)
       REAL PRES(MAXLAY)
       REAL THICK(MAXLAY)
       REAL WAMNT(MAXLAY)
       REAL OAMNT(MAXLAY)
       REAL CAMNT(MAXLAY)
       REAL MAMNT(MAXLAY)
       REAL FAMNT(MAXLAY)
       REAL SAMNT(MAXLAY)
       REAL HAMNT(MAXLAY)
       REAL NAMNT(MAXLAY)
       REAL ALT(MAXLAY)
!
! Profile data structure
! RECORD /RTPPROF/ PROF
       type (RTPPROF) :: PROF
!
       INTEGER ISTAT
       REAL YOFF
       REAL YMODP1(18)
       REAL DOPPPM
       REAL OPHASE
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I, I_REV ! generic
       INTEGER ICO2X ! index for reading CO2 in gamnt
       INTEGER ISO2X ! index for reading SO2 in gamnt
       INTEGER IHNOX ! index for reading HNO3 in gamnt
       INTEGER IN2OX ! index for reading N2O in gamnt
       INTEGER INDEND(4)! indices for big or little endian REAL*4
       INTEGER L ! layer looping
       INTEGER LR ! reversed layer looping
       INTEGER NLEV ! number of levels
       INTEGER rtpread ! for calling read rtp interface routine
       REAL ZSURF ! surface altitude (read but ignored for now)
       REAL RJUNK1 ! generic junk/work
       REAL RJUNK2 ! generic junk/work
       CHARACTER(len=1) :: C4(4)

 real :: esat     ! saturated vapor pressure [hPa]
 real :: ewat     ! water vapor pressure [hPa]
 real :: h2o_ppmv ! water vapor mixing [ppmv] 
 real :: r_earth  ! earth radius [km]
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE begins below
!***********************************************************************
!***********************************************************************

      PROF = PATT  !exchange from PATT to PROF

! Don't need below for G-SDSU. Toshi Matsui, Oct 2013.
!       IF (ISTAT .EQ. -1) GOTO 9999 ! reached end of file
!       IF (.NOT. LWANT) GOTO 9999 ! skip prof if not wanted

!
! --------------------
! Pull out needed data
! --------------------
! Latitude & longitude
       LAT=surf1p%lat   !latitude [deg] PROF%plat
       LON=surf1p%lon   !longitude [deg] PROF%plon

!
! Can not correctly compute gravity if LAT is bad
       IF (ABS(LAT) .GT. 90.01) THEN
          WRITE(IOERR,1005) IPROF, LAT
 1005 FORMAT('ERROR &! input profile PROF(',I4,').plat=',1PE11.4, &
          ' is out of range -90 to +90')
          STOP
       ENDIF

! Note: LON is not currently used and thus need not be checked
!
! Number of levels
       NLEV=size(atmos_stag1p(:)%press)  !          PROF%nlevs

       IF (NLEV .LT. 2 .OR. NLEV .GT. MAXLAY+1) THEN
          WRITE(IOERR,1010) IPROF, NLEV, MAXLAY+1
 1010 FORMAT('ERROR &! input profile PROF(',I4,').nlevs=',I4, &
          ' is out of range 2 to ',I3)
          STOP
       ENDIF
!
! Number of layers
       IF (PTYPE .EQ. AIRSLAY) THEN
! Special case for AIRS pseudo-levels
          NLAY=NLEV
       ELSE
          NLAY=NLEV - 1
       ENDIF

!
! Assign read indices for trace gas CO2
       IF (ICO2 .LT. 1) THEN
          CO2PPM=PROF%co2ppm
          IF (CO2PPM .LT. -998) CO2PPM=CO2STD
          RJUNK1=0.8*CO2STD
          RJUNK2=1.2*CO2STD
          IF (CO2PPM .LT. RJUNK1 .OR. CO2PPM .GT. RJUNK2) THEN
             WRITE(IOERR,1015) IPROF, CO2PPM, RJUNK1, RJUNK2
 1015 FORMAT('Warning &! PROF(',I4,').co2ppm=',1PE10.3, &
             ' is outside allowed range ',0PF5.1,' to ',F5.1)
          ENDIF
! Set ICO2X to any valid gas index; data will be read but ignored
          ICO2X=IH2O
       ELSE
          CO2PPM=-9999
          ICO2X=ICO2
       ENDIF


! Assign read indices for trace gas SO2
       IF (ISO2 .LT. 1) THEN
          ISO2X=IH2O
       ELSE
          ISO2X=ISO2
       ENDIF
! Assign read indices for trace gas HNO3
       IF (IHNO3 .LT. 1) THEN
          IHNOX=IH2O
       ELSE
          IHNOX=IHNO3
       ENDIF
! Assign read indices for trace gas N2O
       IF (IN2O .LT. 1) THEN
          IN2OX=IH2O
       ELSE
          IN2OX=IN2O
       ENDIF


!
! Angles 
!
     if( sol_zen_airs == undefined ) then
       SUNANG = acos( surf1p%cosz ) / const_degrad ! solar zenith angle from model [deg]
     else
       SUNANG = sol_zen_airs  ! solar zenith angle from configuration [deg]
     endif
     SATANG = view_ang_airs ! sensor viewing angle [deg]
     SATZEN = undefined     ! satellite zenith angle 

!
! Satellite altitude above ellipsoid surface (convert m to km)
    call earth_radius_single( LAT*const_degrad, r_earth )

    ZSAT=  7077.75 - r_earth ! 690. altitude of satellite [km]

!
! Surface
!
       PSURF = atmos_stag1p(1)%press ! surface pressure [mb]   PROF%spres
       TSURF = surf1p%t_skin         ! skin temperature [K]    PROF%stemp
       ZSURF = surf1p%elev           ! surface elevation [m]  PROF%salti ! note: ZSURF is currently ignored
       IF (PSURF .LE. 0) THEN
          WRITE(IOERR,1017) IPROF
 1017 FORMAT('ERROR! Prof(',I4,') has no surface pressure')
          STOP
       ENDIF
       IF (TSURF .LE. 0) THEN
          WRITE(IOERR,1018) IPROF
 1018 FORMAT('ERROR! Prof(',I4,') has no surface temperature')
          STOP
       ENDIF



!
! Emissivity (range 0 to 1) and Reflectance (range 0 to 1/pi)
!

  if( emiss_sarta(1) /= undefined) then
       NEMIS=anchor_bs_emiss

       DO I=1,NEMIS
          I_REV = NEMIS - I + 1 
          FEMIS(I)=1. / wave_bs_emiss(I_REV) * 1.e+4  !  wave --> wave number [1/cm]
          EMIS(I) = emiss_sarta(I_REV)    ! emissivity at anchor points
          RHO(I)  = (1. - EMIS(I))/PI
!print*,I,  FEMIS(I), EMIS(I), RHO(I)
       ENDDO
!stop


  else   !default value from table

       NEMIS=nemis_ap

       IF (NEMIS .EQ. 0) THEN
          WRITE(IOERR,1020) IPROF
 1020 FORMAT('ERROR! PROF(',I4,') has no emissivity')
          STOP
       ENDIF
       DO I=1,NEMIS
          FEMIS(I)=efreq_ap(I)   ! frequency of emissivity anchor ponts. 
          EMIS(I)= emis_ap(I)    ! emissivity at anchor points
          IF (rho_ap(I) .LT. 0.0) THEN
             RHO(I)=(1 - EMIS(I))/PI
          ELSE
             RHO(I)=rho_ap(I)
          ENDIF
       ENDDO

 endif

!
! ----------------------------------
! Get layer temperature & gas amount
! ----------------------------------
       IF (GUCIN .EQ. 1) THEN  ! GUC number for: molecules/cm^2
          ! Input gas units are molecules/cm^2; convert to kilomoles/cm^2
          IF (PROF%plevs(1) .LT. PROF%plevs(NLEV)) THEN
          ! Prof is in top-down order(L=1 is TOA;  L=NLAY is BOA)
             DO L=1,NLAY
                TEMP(L)=PROF%ptemp(L)
                WAMNT(L)=PROF%gamnt(L,IH2O )/6.02214199E+26
                OAMNT(L)=PROF%gamnt(L,IO3 )/6.02214199E+26
                CAMNT(L)=PROF%gamnt(L,ICO )/6.02214199E+26
                MAMNT(L)=PROF%gamnt(L,ICH4 )/6.02214199E+26
                FAMNT(L)=PROF%gamnt(L,ICO2X)/6.02214199E+26
                SAMNT(L)=PROF%gamnt(L,ISO2X)/6.02214199E+26
                HAMNT(L)=PROF%gamnt(L,IHNOX)/6.02214199E+26
                NAMNT(L)=PROF%gamnt(L,IN2OX)/6.02214199E+26
                ALT(L)=0.5*( PROF%palts(L) + PROF%palts(L+1) )
             ENDDO

             IF (LCO2PM) THEN
                DO L=1,NLAY
                   FAMNT(L)=PROF%gamnt(L,ICO2X)
                ENDDO
             ENDIF
          ELSE
! Prof is in bottom-up order
             DO L=1,NLAY
                LR=1 + NLAY - L ! reversed layer index
                TEMP(L)=PROF%ptemp(LR)
                WAMNT(L)=PROF%gamnt(LR,IH2O )/6.02214199E+26
                OAMNT(L)=PROF%gamnt(LR,IO3 )/6.02214199E+26
                CAMNT(L)=PROF%gamnt(LR,ICO )/6.02214199E+26
                MAMNT(L)=PROF%gamnt(LR,ICH4 )/6.02214199E+26
                FAMNT(L)=PROF%gamnt(LR,ICO2X)/6.02214199E+26
                SAMNT(L)=PROF%gamnt(LR,ISO2X)/6.02214199E+26
                HAMNT(L)=PROF%gamnt(LR,IHNOX)/6.02214199E+26
                NAMNT(L)=PROF%gamnt(LR,IN2OX)/6.02214199E+26
                ALT(L)=0.5*( PROF%palts(LR) + PROF%palts(LR+1) )
             ENDDO
             IF (LCO2PM) THEN
                DO L=1,NLAY
                   LR=1 + NLAY - L ! reversed layer index
                   FAMNT(L)=PROF%gamnt(LR,ICO2X)
                ENDDO
             ENDIF
          ENDIF
!
       ELSEIF (GUCIN .EQ. 2) THEN
! Input gas units are kilomoles/cm^2
          IF (PROF%plevs(1) .LT. PROF%plevs(NLEV)) THEN
! Prof is in top-down order
             DO L=1,NLAY
                TEMP(L)=PROF%ptemp(L)
                WAMNT(L)=PROF%gamnt(L,IH2O)
                OAMNT(L)=PROF%gamnt(L,IO3)
                CAMNT(L)=PROF%gamnt(L,ICO)
                MAMNT(L)=PROF%gamnt(L,ICH4)
                FAMNT(L)=PROF%gamnt(L,ICO2X)
                SAMNT(L)=PROF%gamnt(L,ISO2X)
                HAMNT(L)=PROF%gamnt(L,IHNOX)
                NAMNT(L)=PROF%gamnt(L,IN2OX)
                ALT(L)=0.5*( PROF%palts(L) + PROF%palts(L+1) )
             ENDDO
          ELSE
! Prof is in bottom-up order
             DO L=1,NLAY
                LR=1 + NLAY - L ! reversed layer index
                TEMP(L)=PROF%ptemp(LR)
                WAMNT(L)=PROF%gamnt(LR,IH2O)
                OAMNT(L)=PROF%gamnt(LR,IO3)
                CAMNT(L)=PROF%gamnt(LR,ICO)
                MAMNT(L)=PROF%gamnt(LR,ICH4)
                FAMNT(L)=PROF%gamnt(LR,ICO2X)
                SAMNT(L)=PROF%gamnt(LR,ISO2X)
                HAMNT(L)=PROF%gamnt(LR,IHNOX)
                NAMNT(L)=PROF%gamnt(LR,IN2OX)
                ALT(L)=0.5*( PROF%palts(LR) + PROF%palts(LR+1) )
             ENDDO
          ENDIF

       ELSEIF( GUCIN .EQ. 3) THEN  ! GUC number for: ppmv  (coupling with G-SDSU)

          DO L=1,NLAY
             LR=1 + NLAY - L ! reversed layer index
                PRES (L)=atmos1p(LR)%press                ! air pressure [mb] = [hPa]
                TEMP (L)=atmos1p(LR)%t_air                ! air temperature [K]      PROF%ptemp(LR)
                ALT  (L)=atmos1p(LR)%hgt *1000.e0         ! height [m] 
                THICK(L)=atmos1p(LR)%dhgt*1000.e0*100.e0  ! thickness [cm]
               
                !
                ! convert water vapor units
                !
                esat = 0.611e0 * exp (const_Lv_Rv* ( 1.e0/const_Kel2Cel - 1.e0/atmos1p(LR)%t_air ) ) *10.e0  ! [hPa]
                ewat = esat * atmos1p(LR)%rh * 1.e-2                      ! water vapor pressure [hPa]
                h2o_ppmv = 1.e6 * ewat / atmos1p(LR)%press                ! water vapor mixing [ppmv] 
                WAMNT(L) = h2o_ppmv*1.e-6 *(atmos1p(LR)%press * 0.01 )/ &
                           (8.315*atmos1p(LR)%t_air )*1.e-9 * THICK(L) *1.e4   ! H2O [kilo mole/cm2]

          ENDDO

          call interp_gas_profile( INT(sdsu_mm), LAT, PRES(1:NLAY), 'O3', 'kmol/cm3', GAS(1:NLAY) )
          DO L=1,NLAY
             OAMNT(L)= GAS(L) * THICK(L)   ! O3 [kmol/cm2]    
          ENDDO

          call interp_gas_profile( INT(sdsu_mm), LAT, PRES(1:NLAY), 'CO', 'kmol/cm3', GAS(1:NLAY) )
          DO L=1,NLAY
             CAMNT(L)= GAS(L) * THICK(L)   ! CO [kmol/cm2]   
          ENDDO

          call interp_gas_profile( INT(sdsu_mm), LAT, PRES(1:NLAY), 'CH4', 'kmol/cm3', GAS(1:NLAY) )
          DO L=1,NLAY
             MAMNT(L)= GAS(L) * THICK(L)   ! CH4 [kmol/cm2]  
          ENDDO

          call interp_gas_profile( INT(sdsu_mm), LAT, PRES(1:NLAY), 'CO2', 'kmol/cm3', GAS(1:NLAY) )
          DO L=1,NLAY
             FAMNT(L)= GAS(L) * THICK(L)   ! CO2 [kmol/cm2]  
          ENDDO

          call interp_gas_profile( INT(sdsu_mm), LAT, PRES(1:NLAY), 'SO2', 'kmol/cm3', GAS(1:NLAY) )
          DO L=1,NLAY
             SAMNT(L)= GAS(L) * THICK(L)   ! SO2 [kmol/cm2]  
          ENDDO

          call interp_gas_profile( INT(sdsu_mm), LAT, PRES(1:NLAY), 'HNO3', 'kmol/cm3', GAS(1:NLAY) )
          DO L=1,NLAY
             HAMNT(L)= GAS(L) * THICK(L)   ! HNO [kmol/cm2]  
          ENDDO

          call interp_gas_profile( INT(sdsu_mm), LAT, PRES(1:NLAY), 'N2O', 'kmol/cm3', GAS(1:NLAY) )
          DO L=1,NLAY
             NAMNT(L)= GAS(L) * THICK(L)   ! N2O [kmol/cm2]  
          ENDDO


!             DO L=1,NLAY
!print*, L, TEMP(L), WAMNT(L)/THICK(L), OAMNT(L)/THICK(L), FAMNT(L)/THICK(L), THICK(L)*1.e-2
!             ENDDO
!stop
!BUG FOUND all gas amount is 1.e4 small, so it should be unit problem


       ELSE ! empty else

         print*,'No such GUCIN =',GUCIN
         stop

       ENDIF
!
! -----------------
! Default altitudes
! -----------------
      PROF%palts(1) = 0.
       IF (PROF%palts(1) .LT. -998 .OR. PTYPE .EQ. AIRSLAY) THEN
! No altitudes, use reference
          DO L=1,NLAY
             ALT(L)=RALT(L)  !this is from default 100-level profile
          ENDDO
       ENDIF

!
! for now ignore clouds
!

       PROF%freqcal = 0
! Read YOFF (grating model y-offset)
       IF (PROF%freqcal .LT. YOFMIN .OR. PROF%freqcal .GT. YOFMAX) THEN
          YOFF=YOFDEF  !use this ... , since PROF%freqcal is undefined.
!cc maybe a bad idea to change the value
! PROF.freqcal=YOFDEF
!cc
       ELSE
          YOFF=PROF%freqcal
       ENDIF

!
! Read module Yoffsets, Doppler, and orbit phase from prof.pnote
! Initialize data to safe values
       DO L=1,18
          YMODP1(L)=0
       ENDDO
       DOPPPM=0
       OPHASE=-9999

    YOFF_IF:  IF (PROF%pnote(1:4) .EQ. 'YOFF') THEN  !not used 
! Try to guess the byte order of REAL*4 data in PROF%pnote
          INDEND(1)=5
          INDEND(2)=6
          INDEND(3)=7
          INDEND(4)=8
          C4(1) = PROF%pnote(INDEND(1):INDEND(1))
          C4(2) = PROF%pnote(INDEND(2):INDEND(2))
          C4(3) = PROF%pnote(INDEND(3):INDEND(3))
          C4(4) = PROF%pnote(INDEND(4):INDEND(4))
          RJUNK1 = C4TOR4(C4)
          IF (RJUNK1 .GT. 1 .OR. RJUNK1 .LT. -9999) THEN
             INDEND(1)=8
             INDEND(2)=7
             INDEND(3)=6
             INDEND(4)=5
             C4(1) = PROF%pnote(INDEND(1):INDEND(1))
             C4(2) = PROF%pnote(INDEND(2):INDEND(2))
             C4(3) = PROF%pnote(INDEND(3):INDEND(3))
             C4(4) = PROF%pnote(INDEND(4):INDEND(4))
             RJUNK2 = C4TOR4(C4)
             IF (RJUNK2 .GT. 1 .OR. RJUNK2 .LT. -9999) THEN
                WRITE(IOERR,1030) RJUNK1, RJUNK2
 1030 FORMAT('Error: pnote data has unknown endian: try1=', &
                   1PE10.3,', try2=',E10.3)
             ENDIF
          ENDIF
!
! Read the data
! Module yoffsets (um)
          DO I=1,17
             C4(1) = PROF%pnote(INDEND(1):INDEND(1))
             C4(2) = PROF%pnote(INDEND(2):INDEND(2))
             C4(3) = PROF%pnote(INDEND(3):INDEND(3))
             C4(4) = PROF%pnote(INDEND(4):INDEND(4))
             YMODP1(I) = C4TOR4(C4)
             DO L=1,4
                INDEND(L)=INDEND(L) + 4
             ENDDO
          ENDDO
! Set fake channels to YOFDEF
          YMODP1(18) = YOFDEF
!
! Doppler (ppm)
          C4(1) = PROF%pnote(INDEND(1):INDEND(1))
          C4(2) = PROF%pnote(INDEND(2):INDEND(2))
          C4(3) = PROF%pnote(INDEND(3):INDEND(3))
          C4(4) = PROF%pnote(INDEND(4):INDEND(4))
          DOPPPM = C4TOR4(C4)
!
! Orbit phase (degrees)
          DO L=1,4
             INDEND(L)=INDEND(L) + 4
          ENDDO
          C4(1) = PROF%pnote(INDEND(1):INDEND(1))
          C4(2) = PROF%pnote(INDEND(2):INDEND(2))
          C4(3) = PROF%pnote(INDEND(3):INDEND(3))
          C4(4) = PROF%pnote(INDEND(4):INDEND(4))
          OPHASE = C4TOR4(C4)
    ENDIF YOFF_IF


!
 9999 RETURN
       END SUBROUTINE RDRTP
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! RDSUN
!
!F77====================================================================
!ROUTINE NAME:
! RDSUN
!ABSTRACT:
! Read in the AIRS solar radiance data
!CALL PROTOCOL
! RDSUN ( AORB, IOUN, INDCHN, HSUN )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! CHAR*1 AORB specify A or B none
! INT arr INDCHN indices of channels none
! INTEGER IOUN I/O unit number none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr HSUN solar radiance W/(m2.str.cm-1)
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
! unit IOUN : input file, ASCII text file. The file is opened,
! read, and closed.
!COMMON BLOCKS
! none
!DESCRIPTION:
! May 2008 version of the 100 layer AIRS Fast Transmittance
! Code by L.Strow/S.Hannon.
!
! Reads in a text file with solar radiance data for each AIRS
! channel. This is the solar rad direct from the sun at the top
! of Earth's atmosphere.
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 12 Sep 1997 Scott Hannon Created
! 12 Feb 2001 Scott Hannon hardcoded filename instead of prompt
! 09 May 2008 Scott Hannon Add AORB
!END====================================================================
! =================================================================
       SUBROUTINE RDSUN ( AORB, IOUN, INDCHN, HSUN )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       CHARACTER*1 AORB
       INTEGER IOUN
       INTEGER INDCHN(MXCHAN)
!
! Output
       REAL HSUN(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       CHARACTER*80 FNSUN
       CHARACTER*80 CLINE
       INTEGER I
       INTEGER IERR
       INTEGER ICHAN
       REAL FRQCHN
       REAL SUNCHN
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Solar filename
       IF (AORB .EQ. 'A') THEN
          FNSUN=FASUN
       ELSE
          FNSUN=FBSUN
       ENDIF
! ----------------------------
! Open the solar radiance file
! ----------------------------
       OPEN(UNIT=IOUN,FILE=FNSUN,FORM='FORMATTED',STATUS='OLD', &
          IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNSUN
 1020 FORMAT('Error ',I5,' openning file:',/,A80)
          STOP
       ENDIF
!
! Initialize the channel counter
       I=0
! -----------------------
! Read the solar rad file
! -----------------------
! Read a line of text from the file
 10 READ(IOUN,9000,END=910) CLINE
 9000 FORMAT(A80)
!
! Determine if the text line is data or a comment
       IF (CLINE(1:1) .NE. '!') THEN
!
! It's data, so increment the channel counter
          I=I+1
!
! Read the data from the text line
          READ(CLINE,*) ICHAN, FRQCHN, SUNCHN
!
! Check to be sure the channel value is OK
          IF ((ICHAN .LT. 1) .OR. (ICHAN .GT. MXCHAN)) THEN
             WRITE(6,1040) MXCHAN, ICHAN
 1040 FORMAT('Error &! Channel number is out of range.',/, &
             'Range is 1 to ',I4,', but channel list has ',I7,'.')
             STOP
          ENDIF
!
! Keep the data if the current channel is on the list
          IF (INDCHN(ICHAN) .NE. 0) THEN
             HSUN( INDCHN(ICHAN) )=SUNCHN
          ENDIF
!
       ENDIF
!
! Read the next line
       IF (I .LT. MXCHAN) GOTO 10
! Note: this routine expects data for every channel
!
 910 CLOSE(IOUN)
!
       RETURN
       END SUBROUTINE RDSUN
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! SACONV
!
!F77====================================================================
!ROUTINE NAME:
! SACONV (real function)
!ABSTRACT:
! Function to convert the surface solar zenith angle SZA into the
! local solar angle at altitude ALT.
!CALL PROTOCOL
! SACONV( SZA, ALT )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL ALT Average layer altitude meters
! REAL SZA Solar Zenith Angle degrees
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL fun SACONV local solar zenith angle radians
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! none
!COMMON BLOCKS
! none
!DESCRIPTION:
! March 1998 version of the 100 layer AIRS Fast Transmittance
! Code by L.Strow/S.Hannon.
!
! ===================================================================
! Function to convert the Solar Zenith Angle SZA at the Earth's
! surface into a the local solar angle at altitude ALT.
! The local solar angle generally varies slightly with altitude
! due to the curvature of the Earth and its atmosphere.
! The effect is largest at the maximum solar zenith angle, and
! disappears as the solar zenith angle approaches 0 degrees.
!
! Currently this function only considers the geometry of the
! situation, and no refractive effects are included.
!
! The layers of the atmosphere may be considered as concentric
! rings with some average altitude. A ray traced thru these rings
! at any viewing angle other than nadir will have a slightly
! different angle (relative to the outward radial at the point
! of intersection) in each ring.
!
! The local solar angle may be calculated (using The Law of
! Sines) if we know:
! The solar zenith angle, SZA, at the Earth's surface (ALT=0)
! The layer altitude, ALT.
! The radius of the Earth, RE.
!
! The solution uses the law of sines and sin(180 - x) = sin(x)
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! No refractive effects have been included.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 25 Feb 1998 Scott Hannon Created
!END====================================================================
! =================================================================
       REAL FUNCTION SACONV( SZA, ALT )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
       REAL SZA
       REAL ALT
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       REAL CONV
       REAL RE
       REAL RA
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ------------------
! Assign some values
! ------------------
! CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02
!
! RE = radius of the Earth (in km)
       RE=6.37E+03
!
! RA = radius of the point to calc the angle at (in km)
! Note: need to convert altitude in meters to kilometers
       RA=RE + (ALT/1000.0)
!
! -----------------
! Do the conversion
! -----------------
!
       SACONV=ASIN( (RE/RA) * SIN(CONV*SZA) )
!
       RETURN
       END FUNCTION SACONV
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! SETEMS
!
!F77====================================================================
!ROUTINE NAME:
! SETEMS
!ABSTRACT:
! Assign the emissivity and reflectivity data for every channel
! being used.
!CALL PROTOCOL
! RDEMIS ( NCHAN, NEMIS, FREQ, FEMIS, XEMIS, XRHO,
! LRHOT, EMIS, RHOSUN, RHOTHR )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER NCHAN number of channels none
! INTEGER NEMIS number of emis pts none
! REAL arr FREQ channel frequencies cm^-1
! REAL arr FEMIS raw emis freq points cm^-1
! REAL arr XEMIS raw emis points none (0 to 1)
! REAL arr XRHO raw reflec points 1/steradian
! LOGICAL LRHOT force refl therm rho? none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr EMIS surface emissivity none
! REAL arr RHOSUN reflectivity for solar 1/steradian
! REAL arr RHOTHR reflectivity for thermal 1/steradian
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2000 version of the 100 layer AIRS Fast Transmittance
! Code by L.Strow/S.Hannon.
!
! Loops over the passed in EMIS & RHO points and interpolates
! or extrapolates them onto the channel freqs being used.
! Interpolations are linear in wavelength (not frequency).
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- -------------------------------------------
! 14 Feb 2001 Scott Hannon Created based on rdemis.f
! 06 Feb 2004 Scott Hannon Add LRHOT argument and associated code
! 24 Oct 2008 Scott Hannon Update for rtpV201: input emis & rho is
! now specified with same freq points
!END====================================================================
! =================================================================
       SUBROUTINE SETEMS ( NCHAN, NEMIS, FREQ, FEMIS, XEMIS, &
          XRHO, LRHOT, EMIS, RHOSUN, RHOTHR )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER NCHAN
       INTEGER NEMIS
       REAL FREQ(MXCHAN)
       REAL FEMIS(MXEMIS)
       REAL XEMIS(MXEMIS)
       REAL XRHO(MXEMIS)
       LOGICAL LRHOT
!
! Output
       REAL EMIS(MXCHAN)
       REAL RHOSUN(MXCHAN)
       REAL RHOTHR(MXCHAN)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I ! generic looping variable
       INTEGER J ! generic looping variable
       INTEGER SORTED ! conrol flag for sorts
       REAL DELELM(MXEMIS) ! delta for emis wavelength
       REAL DELEMS(MXEMIS) ! delta for emis
       REAL DELRHO(MXEMIS) ! delta for rho
       REAL DX ! delta X
       REAL ELAM(MXEMIS) ! emis wavelength
       REAL RJUNK1 ! generic junk/work variable
       REAL RJUNK2 ! generic junk/work variable
       REAL RJUNK3 ! generic junk/work variable
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ---------------------------
! Sort the emis & rho by freq
! ---------------------------
       SORTED=1
 10 IF (SORTED .EQ. 1) THEN
          SORTED=0
          DO I=1,NEMIS-1
             IF (FEMIS(I) .GT. FEMIS(I+1)) THEN
                RJUNK1=FEMIS(I)
                RJUNK2=XEMIS(I)
                RJUNK3=XRHO(I)
                FEMIS(I)=FEMIS(I+1)
                XEMIS(I)=XEMIS(I+1)
                XRHO(I)=XRHO(I+1)
                FEMIS(I+1)=RJUNK1
                XEMIS(I+1)=RJUNK2
                XRHO(I+1)=RJUNK3
                SORTED=1
             ENDIF
          ENDDO
          GOTO 10
       ENDIF
!
! ----------------------------------------
! Calc wavelength and interpolation deltas
! ----------------------------------------
       ELAM(1)=1.0E+4/FEMIS(1) ! wavelength in microns
       DO J=2,NEMIS
          ELAM(J)=1.0E+4/FEMIS(J)
          DELELM(J-1)=ELAM(J) - ELAM(J-1)
          DELEMS(J-1)=XEMIS(J) - XEMIS(J-1)
          DELRHO(J-1)=XRHO(J) - XRHO(J-1)
       ENDDO
!
! ------------------------------------
! Interpolate emis & rho onto channels
! ------------------------------------
! Loop over channels (note: chan freqs may be in any order)
       DO I=1,NCHAN
!
! Emissivity
          IF ( FREQ(I) .LE. FEMIS(1) ) THEN
             EMIS(I)=XEMIS(1)
             RHOSUN(I)=XRHO(1)
          ELSEIF ( FREQ(I) .GE. FEMIS(NEMIS) ) THEN
             EMIS(I)=XEMIS(NEMIS)
             RHOSUN(I)=XRHO(NEMIS)
          ELSE
! Determine the index of the upper bounding FEMIS
             J=2
 20 IF ( FEMIS(J) .LT. FREQ(I) ) THEN
                J=J + 1
                GOTO 20
             ENDIF
! Convert to lower boundary index
             J=J - 1
!
             DX=(( 1.0E+4/FREQ(I) - ELAM(J) )/DELELM(J))
             EMIS(I) = XEMIS(J) + DX*DELEMS(J)
             RHOSUN(I) = XRHO(J) + DX*DELRHO(J)
          ENDIF
!
          IF (LRHOT) THEN
! Force reflected thermal rho = (1-emis)/pi
             RHOTHR(I)=(1.0 - EMIS(I))/PI
          ELSE
! Reflected thermal uses the same rho as solar
             RHOTHR(I)=RHOSUN(I)
          ENDIF
!

       ENDDO

  RETURN
  END SUBROUTINE SETEMS
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! SETYOF
!
!F77====================================================================
!ROUTINE NAME:
! SETYOF
!ABSTRACT:
! Set the Yoffset for each AIRS channel
!CALL PROTOCOL
! SETYOF( NCHAN, LSTCHN, RTIME, YOFF, YMODP1, DOPPPM, ETIME,
! FREQAB, WIDTH, INDMOD, IABWGT, DELYAB, YOFCHN )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER NCHAN number of channels none
! INT arr LSTCHN channel list channel ID
! DOUBLE RTIME observation time TAI1993
! REAL arr YMODP1 Yoffset of modules+1 um
! REAL DOPPPM Doppler shift ppm in frequency
! DOUBLE ar ETIME A/B weight epoch time TAI1993
! INT arr INDMOD module index none
! INT arr IABWGT A/B weight index none
! REAL arr DELYAB Delta Yoffset A,B,AB um
!INPUT/OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL YOFF overall Yoffset um
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr YOFCHN total channel Yoffset um
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
! unit IOUN : input text file of delta Yoffsets and A/B weights
!COMMON BLOCKS
! none
!DESCRIPTION:
! October 2011 version of SARTA v1.08 "df" code by L.Strow/S.Hannon.
!
! Read in a table of delta Yoffset and A/B weights and apply to
! the user input Yoffset and Doppler shift to generate a total
! Yoffset for each channel.
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- -------------------------------------------
! 11 Oct 2011 Scott Hannon Created
!END====================================================================
! =================================================================
       SUBROUTINE SETYOF(NCHAN, LSTCHN, RTIME, YOFF, YMODP1, DOPPPM, &
          ETIME, FREQAB, WIDTH, INDMOD, IABWGT, DELYAB, YOFCHN)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! INPUT
       INTEGER NCHAN ! number of channels in use
       INTEGER LSTCHN(MXCHAN) ! channel list
       DOUBLE PRECISION RTIME ! observation time
       REAL YMODP1(18) ! grating model module Yoffset
       REAL DOPPPM ! Doppler shift in ppm
       DOUBLE PRECISION ETIME(MAXABW) ! A/B weight epoch start time
       REAL FREQAB(MXCHAN) ! frequency
       REAL WIDTH(MXCHAN) ! width
       INTEGER INDMOD(MXCHAN) ! module index
       INTEGER IABWGT(MXCHAN,MAXABW) ! A/B weight index
       REAL DELYAB(MXCHAN,3) ! Delta Yoffset for A,B,AB
! INPUT/OUTPUT
       REAL YOFF ! overall Yoffset
! OUTPUT
       REAL YOFCHN(MXCHAN) ! total channel Yoffset
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I ! generic looping variable
       INTEGER ICHAN ! channel ID
       INTEGER IEPOCH ! A/B weight epoch index
       REAL DOPPLR ! Doppler shift in um
       REAL YSUM ! sum of Yoffsets
!
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Determine A/B weight epoch
       IEPOCH = 1
       DO I=2,MAXABW
          IF (RTIME .GE. ETIME(I)) IEPOCH = I
!cc uncomment for testing
! write(6,1233) RTIME, ETIME(I)
! 1233 format(X,F12.1,X,F12.1)
!cc
       ENDDO
!
! Set the Yoffsets for each channel
       IF (YMODP1(1) .LT. -10 .AND. YMODP1(1) .GT. -18) THEN
          YSUM=0
          DO I=1,NCHAN
             ICHAN = LSTCHN(I)
!
! Convert Doppler in ppm to um.
! Note: the WIDTH is equivalent to 100 um of delta Yoffset,
! and WIDTH = FREQAB/resolving_power. Freqs get larger with
! decreasing Yoffset, so a positive ppm shift is a negative
! um shift.
!CC Use nominal resolving power = 1200
             DOPPLR = -0.12*DOPPPM ! 0.12 = 1E-6*100*FREQ/(FREQ/1200)
!CC Use the actual resolving power
! DOPPLR = -1E-4*DOPPPM*FREQAB(ICHAN)/WIDTH(ICHAN)
!CC
             YOFCHN(I) = YMODP1(INDMOD(ICHAN)) + DOPPLR &
                + DELYAB(ICHAN,IABWGT(ICHAN,IEPOCH))
             YSUM = YSUM + YOFCHN(I)
!cc uncomment for testing
! if (ICHAN .EQ. 22) then
! write(6,1234) ICHAN, INDMOD(ICHAN), IEPOCH,
! $ IABWGT(ICHAN,IEPOCH), YMODP1(INDMOD(ICHAN)),
! $ DOPPLR, DELYAB(ICHAN,IABWGT(ICHAN,IEPOCH))
! 1234 format('ID, iModule, iEpoch, iABwgt, YmodP1, Doppler, dYab=',
! $ X,I4,X,I2,X,I1,X,I1,X,F7.3,X,F6.3,X,F6.3)
! endif
!cc
          ENDDO
          YOFF=YSUM/FLOAT(NCHAN)
       ELSE

          IF (YOFF .LT. -10 .AND. YOFF .GT. -18) THEN
             DO I=1,NCHAN
                YOFCHN(I)=YOFF  !this is default
             ENDDO
          ELSE
             YOFF=YOFDEF
             DO I=1,NCHAN
                YOFCHN(I)=YOFDEF
             ENDDO
          ENDIF
       ENDIF
!
!
       RETURN
       END SUBROUTINE SETYOF
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! SUMCOF
!
!F77====================================================================
!ROUTINE NAME:
! SUMCOF
!ABSTRACT:
! Do a weighted sum of two coef/etc databases ("A" and "B")
!CALL PROTOCOL
! SUMCOF( NCHAN, NCHNTE, NCHN1, NCHN2, NCHN3, NCHN4, NCHN5,
! $ NCHN6, NCHN7, NCHCO2, NCHN2O, NCHSO2, NCHHNO, NCHH2O,
! $ LSTCHN, CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
! $ CLISTN, CLICO2, CLIN2O, CLISO2, CLIHNO, CLIH2O, YOFCHN,
! $ FREQA, HSUNA, COEFFA, COEFNA, COEF1A, COEF2A, COEF3A, COEF4A,
! $ COEF5A, COEF6A, COEF7A, COFCOA, COFN2A, COFSOA, COFHNA, COFH2A,
! $ FREQB, HSUNB, COEFFB, COEFNB, COEF1B, COEF2B, COEF3B, COEF4B,
! $ COEF5B, COEF6B, COEF7B, COFCOB, COFN2B, COFSOB, COFHNB, COFH2B,
! $ FREQ, HSUN, COEFF, COEFN, COEF1, COEF2, COEF3, COEF4,
! $ COEF5, COEF6, COEF7, COFCO2, COFN2O, COFSO2, COFHNO, COFH2O)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER NCHAN number of channels none
! INTEGER NCHNTE non-LTE number of channels none
! INTEGER NCHN1 set1 number of channels none
! INTEGER NCHN2 set2 number of channels none
! INTEGER NCHN3 set3 number of channels none
! INTEGER NCHN4 set4 number of channels none
! INTEGER NCHN5 set5 number of channels none
! INTEGER NCHN6 set6 number of channels none
! INTEGER NCHN7 set7 number of channels none
! INTEGER NCHCO2 number of CO2 channels none
! INTEGER NCHN2O number of N2O channels none
! INTEGER NCHSO2 number of SO2 channels none
! INTEGER NCHHNO number of HNO3 channels none
! INTEGER NCHH2O number of OPTRAN H2O chans none
! INT arr LSTCHN channel list channel ID
! INT arr CLIST1 set1 channel list channel ID
! INT arr CLIST2 set2 channel list channel ID
! INT arr CLIST3 set3 channel list channel ID
! INT arr CLIST4 set4 channel list channel ID
! INT arr CLIST5 set5 channel list channel ID
! INT arr CLIST6 set6 channel list channel ID
! INT arr CLIST7 set7 channel list channel ID
! INT arr CLISTN non-LTE channel list channel ID
! INT arr CLICO2 CO2 channel list channel ID
! INT arr CLIN2O N2O channel list channel ID
! INT arr CLISO2 SO2 channel list channel ID
! INT arr CLIHNO HNO3 channel list channel ID
! INT arr CLIH2O H2O channel list channel ID
! REAL arr YOFCHN channel total Yoffset um
!
! REAL arr FREQA channel frequencies cm^-1
! REAL arr HSUNA solar irradiance mW/m^2/cm^-1 ?
! REAL arr COEFFA refl thermal "F" coefs various
! REAL arr COEFNA nonLTE coefs various
! REAL arr COEF1A set1 fast trans coefs various
! REAL arr COEF2A set2 fast trans coefs various
! REAL arr COEF3A set3 fast trans coefs various
! REAL arr COEF4A set4 fast trans coefs various
! REAL arr COEF5A set5 fast trans coefs various
! REAL arr COEF6A set6 fast trans coefs various
! REAL arr COEF7A set7 fast trans coefs various
! REAL arr COFCOA CO2 perturbation coefs various
! REAL arr COFN2A N2O perturbation coefs various
! REAL arr COFSOA SO2 perturbation coefs various
! REAL arr COFHNA HNO3 perturbation coefs various
! REAL arr COFH2A OPTRAN H2O coefs various
!
! REAL arr FREQB channel frequencies cm^-1
! REAL arr HSUNB solar irradiance mW/m^2/cm^-1 ?
! REAL arr COEFFB refl thermal "F" coefs various
! REAL arr COEFNB nonLTE coefs various
! REAL arr COEF1B set1 fast trans coefs various
! REAL arr COEF2B set2 fast trans coefs various
! REAL arr COEF3B set3 fast trans coefs various
! REAL arr COEF4B set4 fast trans coefs various
! REAL arr COEF5B set5 fast trans coefs various
! REAL arr COEF6B set6 fast trans coefs various
! REAL arr COEF7B set7 fast trans coefs various
! REAL arr COFCOB CO2 perturbation coefs various
! REAL arr COFN2B N2O perturbation coefs various
! REAL arr COFSOB SO2 perturbation coefs various
! REAL arr COFHNB HNO3 perturbation coefs various
! REAL arr COFH2B OPTRAN H2O coefs various
!INPUT/OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr FREQ channel frequencies cm^-1
! REAL arr HSUN solar irradiance mW/m^2/cm^-1 ?
! REAL arr COEFF refl thermal "F" coefs various
! REAL arr COEFN nonLTE coefs various
! REAL arr COEF1 set1 fast trans coefs various
! REAL arr COEF2 set2 fast trans coefs various
! REAL arr COEF3 set3 fast trans coefs various
! REAL arr COEF4 set4 fast trans coefs various
! REAL arr COEF5 set5 fast trans coefs various
! REAL arr COEF6 set6 fast trans coefs various
! REAL arr COEF7 set7 fast trans coefs various
! REAL arr COFCO2 CO2 perturbation coefs various
! REAL arr COFN2O N2O perturbation coefs various
! REAL arr COFSO2 SO2 perturbation coefs various
! REAL arr COFHNO HNO3 perturbation coefs various
! REAL arr COFH2O OPTRAN H2O coefs various
!OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! SARTA
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
! unit IOUN : input text file of tuning adjustments.
!COMMON BLOCKS
! none
!DESCRIPTION:
! October 2011 version of SARTA v1.08 "df" code by L.Strow/S.Hannon.
!
! The routine takes a pair of fast transmittace coefficients
! and sums them.
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- -------------------------------------------
! 09 May 2008 Scott Hannon Created
! 03 Aug 2009 Scott Hannon Added DFCAL; add LSTCHN & CLI<*> to call
! 11 Oct 2011 Scott Hannon Remove DFCAL and replace YOFF with YOFCHN
!END====================================================================
! =================================================================
       SUBROUTINE SUMCOF(NCHAN, NCHNTE, NCHN1, NCHN2, NCHN3,NCHN4,NCHN5, &
          NCHN6, NCHN7, NCHCO2, NCHN2O, NCHSO2, NCHHNO, NCHH2O, &
         LSTCHN, CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
         CLISTN, CLICO2, CLIN2O, CLISO2, CLIHNO, CLIH2O, YOFCHN, &
          FREQA, HSUNA, COEFFA, COEFNA, COEF1A, COEF2A, COEF3A, COEF4A, &
         COEF5A, COEF6A, COEF7A, COFCOA, COFN2A, COFSOA, COFHNA, COFH2A, &
          FREQB, HSUNB, COEFFB, COEFNB, COEF1B, COEF2B, COEF3B, COEF4B, &
         COEF5B, COEF6B, COEF7B, COFCOB, COFN2B, COFSOB, COFHNB, COFH2B, &
           FREQ, HSUN, COEFF, COEFN, COEF1, COEF2, COEF3, COEF4, &
          COEF5, COEF6, COEF7, COFCO2, COFN2O, COFSO2, COFHNO, COFH2O)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! INPUT
       INTEGER NCHAN ! number of channels in use
       INTEGER NCHNTE ! number of nonLTE channels in use
       INTEGER NCHN1 ! number of set1 channels in use
       INTEGER NCHN2 ! number of set2 channels in use
       INTEGER NCHN3 ! number of set3 channels in use
       INTEGER NCHN4 ! number of set4 channels in use
       INTEGER NCHN5 ! number of set5 channels in use
       INTEGER NCHN6 ! number of set6 channels in use
       INTEGER NCHN7 ! number of set7 channels in use
       INTEGER NCHCO2 ! number of CO2 perturbation channels in use
       INTEGER NCHN2O ! number of N2O perturbation channels in use
       INTEGER NCHSO2 ! number of SO2 perturbation channels in use
       INTEGER NCHHNO ! number of HNO3 perturbation channels in use
       INTEGER NCHH2O ! number of OPTRAN H2O channels in use
       INTEGER LSTCHN(MXCHAN) ! channel list
       INTEGER CLIST1(MXCHN1) ! set1 channel list
       INTEGER CLIST2(MXCHN2) ! set2 channel list
       INTEGER CLIST3(MXCHN3) ! set3 channel list
       INTEGER CLIST4(MXCHN4) ! set4 channel list
       INTEGER CLIST5(MXCHN5) ! set5 channel list
       INTEGER CLIST6(MXCHN6) ! set6 channel list
       INTEGER CLIST7(MXCHN7) ! set7 channel list
       INTEGER CLISTN(MXCNTE) ! non-LTE channel list
       INTEGER CLICO2(MXCHNC) ! CO2 channel list
       INTEGER CLIN2O(MXCHNN) ! N2O channel list
       INTEGER CLISO2(MXCHNS) ! SO2 channel list
       INTEGER CLIHNO(MXCHNH) ! HNO3 channel list
       INTEGER CLIH2O(MXCHNW) ! H2O channel list
!
       REAL YOFCHN(MXCHAN) ! channel total Yoffset
!
! Database "A"
       REAL FREQA(MXCHAN) ! channel frequencies
       REAL HSUNA(MXCHAN) ! solar irradiance
       REAL COEFFA(NFCOEF,MXCHAN) ! reflected thermal "F" coefs
       REAL COEFNA(NNCOEF,MXCNTE) ! nonLTE coefs
       REAL COEF1A(N1COEF,MAXLAY,MXCHN1) ! set1 coefs
       REAL COEF2A(N2COEF,MAXLAY,MXCHN2) ! set2 coefs
       REAL COEF3A(N3COEF,MAXLAY,MXCHN3) ! set3 coefs
       REAL COEF4A(N4COEF,MAXLAY,MXCHN4) ! set4 coefs
       REAL COEF5A(N5COEF,MAXLAY,MXCHN5) ! set5 coefs
       REAL COEF6A(N6COEF,MAXLAY,MXCHN6) ! set6 coefs
       REAL COEF7A(N7COEF,MAXLAY,MXCHN7) ! set7 coefs
       REAL COFCOA( NCO2,MAXLAY,MXCHNC) ! CO2 pert coefs
       REAL COFN2A( NN2O,MAXLAY,MXCHNN) ! N2O pert coefs
       REAL COFSOA( NSO2,MAXLAY,MXCHNS) ! SO2 pert coefs
       REAL COFHNA( NHNO3,MAXLAY,MXCHNH) ! HNO3 pert coefs
       REAL COFH2A( NH2O,MXOWLY,MXCHNW) ! H2O OPTRAN coefs
!
! Database "B"
       REAL FREQB(MXCHAN) ! channel frequencies
       REAL HSUNB(MXCHAN) ! solar irradiance
       REAL COEFFB(NFCOEF,MXCHAN) ! reflected thermal "F" coefs
       REAL COEFNB(NNCOEF,MXCNTE) ! nonLTE coefs
       REAL COEF1B(N1COEF,MAXLAY,MXCHN1) ! set1 coefs
       REAL COEF2B(N2COEF,MAXLAY,MXCHN2) ! set2 coefs
       REAL COEF3B(N3COEF,MAXLAY,MXCHN3) ! set3 coefs
       REAL COEF4B(N4COEF,MAXLAY,MXCHN4) ! set4 coefs
       REAL COEF5B(N5COEF,MAXLAY,MXCHN5) ! set5 coefs
       REAL COEF6B(N6COEF,MAXLAY,MXCHN6) ! set6 coefs
       REAL COEF7B(N7COEF,MAXLAY,MXCHN7) ! set7 coefs
       REAL COFCOB( NCO2,MAXLAY,MXCHNC) ! CO2 pert coefs
       REAL COFN2B( NN2O,MAXLAY,MXCHNN) ! N2O pert coefs
       REAL COFSOB( NSO2,MAXLAY,MXCHNS) ! SO2 pert coefs
       REAL COFHNB( NHNO3,MAXLAY,MXCHNH) ! HNO3 pert coefs
       REAL COFH2B( NH2O,MXOWLY,MXCHNW) ! H2O OPTRAN coefs
! INPUT/OUTPUT
! Weighted sum database to be use for calcs
       REAL FREQ(MXCHAN) ! channel frequencies
       REAL HSUN(MXCHAN) ! solar irradiance
       REAL COEFF(NFCOEF,MXCHAN) ! reflected thermal "F" coefs
       REAL COEFN(NNCOEF,MXCNTE) ! nonLTE coefs
       REAL COEF1(N1COEF,MAXLAY,MXCHN1) ! set1 coefs
       REAL COEF2(N2COEF,MAXLAY,MXCHN2) ! set2 coefs
       REAL COEF3(N3COEF,MAXLAY,MXCHN3) ! set3 coefs
       REAL COEF4(N4COEF,MAXLAY,MXCHN4) ! set4 coefs
       REAL COEF5(N5COEF,MAXLAY,MXCHN5) ! set5 coefs
       REAL COEF6(N6COEF,MAXLAY,MXCHN6) ! set6 coefs
       REAL COEF7(N7COEF,MAXLAY,MXCHN7) ! set7 coefs
       REAL COFCO2( NCO2,MAXLAY,MXCHNC) ! CO2 pert coefs
       REAL COFN2O( NN2O,MAXLAY,MXCHNN) ! N2O pert coefs
       REAL COFSO2( NSO2,MAXLAY,MXCHNS) ! SO2 pert coefs
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH) ! HNO3 pert coefs
       REAL COFH2O( NH2O,MXOWLY,MXCHNW) ! H2O OPTRAN coefs
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I ! generic looping variable
       INTEGER J ! generic looping variable
       INTEGER L ! generic looping variable
       INTEGER ICHAN ! channel ID
       REAL WGTA(MXCHAN) ! fractional weight of database "A"
       REAL WGTB(MXCHAN) ! fractional weight of database "B"
       REAL RJUNK ! generic work variable
!
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Calc summing weights
       DO I=1,NCHAN
          ICHAN = LSTCHN(I)
          RJUNK = YOFCHN(I)
          WGTA(ICHAN)=(RJUNK - YOFFB)/(YOFFA - YOFFB)
          WGTB(ICHAN)=(RJUNK - YOFFA)/(YOFFB - YOFFA)
!
! Adjust freq, sun, and reflected thermal "F" coefs
          FREQ(I)=FREQA(I)*WGTA(ICHAN) + FREQB(I)*WGTB(ICHAN)
          HSUN(I)=HSUNA(I)*WGTA(ICHAN) + HSUNB(I)*WGTB(ICHAN)
          DO J=1,NFCOEF
             COEFF(J,I)=COEFFA(J,I)*WGTA(ICHAN) + &
                        COEFFB(J,I)*WGTB(ICHAN)
          ENDDO
       ENDDO
! Set 1
       DO I=1,NCHN1
          ICHAN = CLIST1(I)
          DO L=1,MAXLAY
             DO J=1,N1COEF
                COEF1(J,L,I)=COEF1A(J,L,I)*WGTA(ICHAN) + &
                             COEF1B(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! Set 2
       DO I=1,NCHN2
          ICHAN = CLIST2(I)
          DO L=1,MAXLAY
             DO J=1,N2COEF
                COEF2(J,L,I)=COEF2A(J,L,I)*WGTA(ICHAN)+ &
                             COEF2B(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! Set 3
       DO I=1,NCHN3
          ICHAN = CLIST3(I)
          DO L=1,MAXLAY
             DO J=1,N3COEF
                COEF3(J,L,I)=COEF3A(J,L,I)*WGTA(ICHAN) + &
                             COEF3B(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! Set 4
       DO I=1,NCHN4
          ICHAN = CLIST4(I)
          DO L=1,MAXLAY
             DO J=1,N4COEF
                COEF4(J,L,I)=COEF4A(J,L,I)*WGTA(ICHAN) + &
                             COEF4B(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! Set 5
       DO I=1,NCHN5
          ICHAN = CLIST5(I)
          DO L=1,MAXLAY
             DO J=1,N5COEF
                COEF5(J,L,I)=COEF5A(J,L,I)*WGTA(ICHAN) + &
                             COEF5B(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! Set 6
       DO I=1,NCHN6
          ICHAN = CLIST6(I)
          DO L=1,MAXLAY
             DO J=1,N6COEF
                COEF6(J,L,I)=COEF6A(J,L,I)*WGTA(ICHAN) + &
                             COEF6B(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! Set 7
       DO I=1,NCHN7
          ICHAN = CLIST7(I)
          DO L=1,MAXLAY
             DO J=1,N7COEF
                COEF7(J,L,I)=COEF7A(J,L,I)*WGTA(ICHAN) + &
                             COEF7B(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! OPTRAN
       DO I=1,NCHH2O
          ICHAN = CLIH2O(I)
          DO L=1,MXOWLY
             DO J=1,NH2O
                COFH2O(J,L,I)=COFH2A(J,L,I)*WGTA(ICHAN) + &
                              COFH2B(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! CO2 perturbation coefficients
       DO I=1,NCHCO2
          ICHAN = CLICO2(I)
          DO L=1,MAXLAY
             DO J=1,NCO2
                COFCO2(J,L,I)=COFCOA(J,L,I)*WGTA(ICHAN) + &
                              COFCOB(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! N2O perturbation coefficients
       DO I=1,NCHN2O
          ICHAN = CLIN2O(I)
          DO L=1,MAXLAY
             DO J=1,NN2O
                COFN2O(J,L,I)=COFN2A(J,L,I)*WGTA(ICHAN) + &
                              COFN2B(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! SO2 perturbation coefficients
       DO I=1,NCHSO2
          ICHAN = CLISO2(I)
          DO L=1,MAXLAY
             DO J=1,NSO2
                COFSO2(J,L,I)=COFSOA(J,L,I)*WGTA(ICHAN) + &
                              COFSOB(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! HNO3 perturbation coefficients
       DO I=1,NCHHNO
          ICHAN = CLIHNO(I)
          DO L=1,MAXLAY
             DO J=1,NHNO3
                COFHNO(J,L,I)=COFHNA(J,L,I)*WGTA(ICHAN) + &
                              COFHNB(J,L,I)*WGTB(ICHAN)
             ENDDO
          ENDDO
       ENDDO
! Adjust non-LTE coefs
       DO I=1,NCHNTE
          ICHAN = CLISTN(I)
          DO J=1,NNCOEF
             COEFN(J,I)=COEFNA(J,I)*WGTA(ICHAN) + &
                        COEFNB(J,I)*WGTB(ICHAN)
          ENDDO
       ENDDO
!
       RETURN
       END SUBROUTINE SUMCOF
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! SUNPAR version with trace gases
!
!F77====================================================================
!ROUTINE NAME:
! SUNPAR
!ABSTRACT:
! Calculate the fast transmittance code temperature/amount/angle
! dependent predictors for a profile at the effective sun angle.
!CALL PROTOCOL:
! SUNPAR ( LBOT, RTEMP, RWAMNT, ROAMNT, RCAMNT,
! $ PTEMP, PWAMNT, POAMNT, PCAMNT,
! $ PRES, SECANG, CONPRD,
! $ FPRED4, FPRED5, FPRED6, FPRED7,
! $ WPRED4, WPRED5, WPRED6, WPRED7,
! $ OPRED4, OPRED5, OPRED6, OPRED7,
! $ CPRED4, TRCPRD )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER LBOT bottom layer number none
! REAL arr PTEMP profile temperature Kelvin
! REAL arr PCAMNT prof carbon monoxide amnt kiloMoles/cm^2
! REAL arr POAMNT profile ozone amount kiloMoles/cm^2
! REAL arr PRES layer pressures atmospheres
! REAL arr PWAMNT profile water amount kiloMoles/cm^2
! REAL arr RTEMP reference temperature Kelvin
! REAL arr RCAMNT ref carbon monoxide amount kiloMoles/cm^2
! REAL arr ROAMNT reference ozone amount kiloMoles/cm^2
! REAL arr RWAMNT reference water amount kiloMoles/cm^2
! REAL arr SECANG secant of path angle none
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr CPRED4 carbon monoxide pred set4 various
! REAL arr FPRED4 fixed predictors set4 various
! REAL arr FPRED5 fixed predictors set5 various
! REAL arr FPRED6 fixed predictors set6 various
! REAL arr FPRED7 fixed predictors set7 various
! REAL arr OPRED4 ozone predictors set4 various
! REAL arr OPRED5 ozone predictors set5 various
! REAL arr OPRED6 ozone predictors set6 variou
! REAL arr OPRED7 ozone predictors set7 various
! REAL arr TRCPRD trace gas pert predictors various
! REAL arr WPRED4 water predictors set4 various
! REAL arr WPRED5 water predictors set5 various
! REAL arr WPRED6 water predictors set6 various
! REAL arr WPRED7 water predictors set7 various
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2000 version of the 100 layer AIRS Fast Transmittance
! Code by L.L.Strow/S.Hannon.
!
! Rapid transmittace algorithm predictors consisting of various gas
! amount and temperature ratios and offsets relative to a reference
! profile are calculated. Only sets 4 - 7 are calculated, as these
! are the only sets fit for extreme sun angles (up to secant 9).
!
! ===================================================================
! The FTC profile variables computed for each layer are:
!
! ---------------------------------
! CONPRD: water continuum predictors (7 terms)
! 1) a*W/Tr^2 2) a*(W/Tr^2)^2 3) a*W/Tr 4) a*W^2/Tr
! 5) a*(W/Tr)^2 6) a*W/Tr^4 7) a*Wr
!
! -------------------------------
! Fixed predictors
!
! FPRED4: FCOW (11 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) a^2*Trz
! 9) a^2*Tr 10) a^3 11) sqrt(a)
!
! FPRED5: FWO (11 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) a*Trz/Tr
! 9) a^2*Tr 10) sqrt(a) 11) Trz
!
! FPRED6: FWO (8 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) sqrt(a)
!
! FPRED7: FWO (8 terms):
! 1) a 2) a^2 3) a*Tr 4) a*Tr^2
! 5) Tr 6) Tr^2 7) a*Trz 8) sqrt(a)
!
! ---------------------------------
! Water predictors
!
! WPRED4: FCOW (13 terms):
! 1) W*a 2) W 3) sqrt(W*a)
! 4) W*a*dT 5) (W*a)^2 6) sqrt(W*a)*dT
! 7) root^4(W*a) 8) W*a*W/Wz 9) W*a^2
! 10) (W*a)^3 11) W*a*Cz*a 12) sqrt(W*a)*W/Wz
! 13) W*a^2*dT
!
! WPRED5: FWO bfsw (3 terms):
! 1) W*a 2) (W*a)^3/2 3) W*a*dT
!
! WPRED6: FWO mfmw (7 terms):
! 1) W*a 2) (W*a)^3/2 3) W*a*dT
! 4) (W*a)^2 5) (W*a)^3/2*dT 6) (W*a)^3
! 7) W*a^2
!
! WPRED7: FWO mfbw (13 terms):
! 1) W*a 2) (W*a)^3/2 3) W*a*dT
! 4) (W*a)^2 5) (W*a)^3/2*dT 6) (W*a)^3
! 7) W*a^2 8) W*a*W/Wz 9) (W*a)^3/2*W/Wz
! 10) (W*a)^5/4 11) (W*a)^2*W/Wz 12) W^2*a
! 13) (W*a)^7/4
!
! ---------------------------
! Ozone predictors
!
! OPRED4: FCOW (3 terms):
! 1) O*a 2) sqrt(O*a) 3) O*a*dT
!
! OPRED5: FWO bfsw (1 term):
! 1) O*a
!
! OPRED6: FWO mfmw (1 term):
! 1) O*a
!
! OPRED7: FWO mfbw (1 term):
! 1) O*a
!
! ---------------------------
! CPRED4: carbon monoxide predictors (11 terms):
! 1) C*a 2) sqrt(C*a) 3) C*a*dT
! 4) (C*a)^2 5) C*a*C/Cz 6) sqrt(C*a)*dT
! 7) root^4(C*a) 8) sqrt(C*a)*C/Cz 9) C
!
! ---------------------------
! CO2PRD: CO2 perturbation coefs (4 terms):
! 1) a 2) Tr 3) a*Tr 4) a*Tr^2
!
! -----
! where:
! "a" is the secant of the viewing angle SECANG
! "Tr" is the temperature ratio PTEMP/RTEMP
! "Trz" is the pressure weighted temperature ratio above, i.e.
! the sum i=2 to i=L of { P(i) * ( P(i) - P(i-1) )* Tr(i-1) }
! where "P" is the pressure PRES and "L" is the layer number, and
! Trz(L=1)=0
! "W" is the water amount ratio PWAMNT/RWAMNT
! "dT" is the temperature offset PTEMP-RTEMP
! "Wz" is the pressure weighted water amount above ratio, the
! sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PWAMNT(i) },
! divided by the same sum except using RWAMNT instead of PWAMNT.
! For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
! "O" is the ozone amount ratio POAMNT/ROAMNT
! "C" is the carbon monoxide amount ratio POAMNT/ROAMNT
! "Cz" is the pressure weighted CO amount above ratio, the
! sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PCAMNT(i) },
! divided by the same sum except using RCAMNT instead of PCAMNT.
! For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! Assumes vaguely realistic profile amounts and temperatures, else
! there might be divide by zero problems, etc.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 27 Aug 1997 Scott Hannon Created from calpar
! 30 Sep 1997 Scott Hannon Added variable CO2
! 26 Aug 1998 Scott Hannon Add LBOT to call; loop on LBOT instead
! of MAXLAY
! 24 Aug 2000 Scott Hannon Remove FIXMUL (calc'ed in CALPAR)
! 18 Sep 2002 Scott Hannon Add predictors 6 & 7 to H2O con
! 25 Apr 2003 Scott Hannon Add SO2
! 23 Jun 2005 Scott Hannon "trace" version for CO2, SO2, & HNO3,
! with all using the same predictors.
! 13 Oct 2005 S.Hannon/C.Barnet bug fix: assign TRCPRD 1-7 (was 1-4)
!END====================================================================
! =================================================================
       SUBROUTINE SUNPAR ( LBOT, &
          RTEMP, RWAMNT, ROAMNT, RCAMNT, &
          PTEMP, PWAMNT, POAMNT, PCAMNT, &
          PRES, SECANG, CONPRD, &
          FPRED4, FPRED5, FPRED6, FPRED7, &
          WPRED4, WPRED5, WPRED6, WPRED7, &
          OPRED4, OPRED5, OPRED6, OPRED7, &
          CPRED4, TRCPRD )
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input
       INTEGER LBOT
       REAL RTEMP(MAXLAY)
       REAL RWAMNT(MAXLAY)
       REAL ROAMNT(MAXLAY)
       REAL RCAMNT(MAXLAY)
       REAL PTEMP(MAXLAY)
       REAL PWAMNT(MAXLAY)
       REAL POAMNT(MAXLAY)
       REAL PCAMNT(MAXLAY)
       REAL PRES(MAXLAY)
       REAL SECANG(MAXLAY)
!
! Output
       REAL CONPRD( N1CON,MAXLAY)
       REAL FPRED4( N4FIX,MAXLAY)
       REAL FPRED5( N5FIX,MAXLAY)
       REAL FPRED6( N6FIX,MAXLAY)
       REAL FPRED7( N7FIX,MAXLAY)
       REAL WPRED4( N4H2O,MAXLAY)
       REAL WPRED5( N5H2O,MAXLAY)
       REAL WPRED6( N6H2O,MAXLAY)
       REAL WPRED7( N7H2O,MAXLAY)
       REAL OPRED4( N4O3,MAXLAY)
       REAL OPRED5( N5O3,MAXLAY)
       REAL OPRED6( N6O3,MAXLAY)
       REAL OPRED7( N7O3,MAXLAY)
       REAL CPRED4( N4CO,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER L
       REAL PDP
       REAL PNORM
       REAL DT
       REAL TR
       REAL TZ
       REAL TRZ
! REAL A_F ! unused so removed 14 Feb 2001
       REAL A_W
       REAL WZREF
       REAL WZ
       REAL AZ_W
       REAL A_O
       REAL A_C
       REAL CZ
       REAL CZREF
       REAL AZ_C
       REAL TJUNKS
       REAL WJUNKA
       REAL WJUNKR
       REAL WJUNKS
       REAL WJUNKZ
       REAL WJUNK4
       REAL OJUNKA
       REAL CJUNKA
       REAL CJUNKR
       REAL CJUNKS
       REAL CJUNKZ
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
!
! Initialize the sum terms to zero
       PNORM=0.0E+0
       TZ=0.0E+0
       WZREF=0.0E+0
       WZ=0.0E+0
       CZREF=0.0E+0
       CZ=0.0E+0
!
! --------------------
! Loop over the layers
! --------------------
       DO L=1,LBOT
!
! ---------------------------
! Calculate the basic profile
! dependent predictors.
! ---------------------------
!
          IF (L .EQ. 1) THEN
             PDP=PRES(1)*( PRES(2) - PRES(1))
             TRZ=0.0E+0
          ELSE
             PDP=PRES(L)*( PRES(L) - PRES(L-1) )
             PNORM=PNORM + PDP
!
! Note: TRZ, TOZ, and TMZ use layer-above terms
             TZ=TZ + PDP*TR
             TRZ=TZ/PNORM
          ENDIF
!
! Temperature terms
          DT=PTEMP(L) - RTEMP(L)
          TR=PTEMP(L)/RTEMP(L)
!
! Water terms
          A_W=PWAMNT(L)/RWAMNT(L)
          WZREF=WZREF + PDP*RWAMNT(L)
          WZ=WZ + PDP*PWAMNT(L)
          AZ_W=WZ/WZREF
!
! Ozone terms
          A_O=POAMNT(L)/ROAMNT(L)
!
! Carbon monoxide terms
          A_C=PCAMNT(L)/RCAMNT(L)
          CZREF=CZREF + PDP*RCAMNT(L)
          CZ=CZ + PDP*PCAMNT(L)
          AZ_C=CZ/CZREF
!
! ----------------------
! Load up the predictors
! ----------------------
!
! -----
! Fixed
! -----
          TJUNKS=TR*TR
!
          FPRED4(1,L)=SECANG(L)
          FPRED4(2,L)=SECANG(L)*SECANG(L)
          FPRED4(3,L)=SECANG(L)*TR
          FPRED4(4,L)=SECANG(L)*TJUNKS
          FPRED4(5,L)=TR
          FPRED4(6,L)=TJUNKS
          FPRED4(7,L)=SECANG(L)*TRZ
          FPRED4(8,L)=SECANG(L)*SECANG(L)*TRZ
          FPRED4(9,L)=SECANG(L)*SECANG(L)*TR
          FPRED4(10,L)=SECANG(L)*SECANG(L)*SECANG(L)
          FPRED4(11,L)=SQRT(SECANG(L))
!
! Fixed predictors for FWO sun bfsw = set5
          FPRED5(1,L)=SECANG(L)
          FPRED5(2,L)=SECANG(L)*SECANG(L)
          FPRED5(3,L)=SECANG(L)*TR
          FPRED5(4,L)=SECANG(L)*TJUNKS
          FPRED5(5,L)=TR
          FPRED5(6,L)=TJUNKS
          FPRED5(7,L)=SECANG(L)*TRZ
          FPRED5(8,L)=SECANG(L)*TRZ/TR
          FPRED5(9,L)=SECANG(L)*SECANG(L)*TR
          FPRED5(10,L)=SQRT(SECANG(L))
          FPRED5(11,L)=TRZ
!
! Fixed predictors for FWO sun mfmw = set6
          FPRED6(1,L)=SECANG(L)
          FPRED6(2,L)=SECANG(L)*SECANG(L)
          FPRED6(3,L)=SECANG(L)*TR
          FPRED6(4,L)=SECANG(L)*TJUNKS
          FPRED6(5,L)=TR
          FPRED6(6,L)=TJUNKS
          FPRED6(7,L)=SECANG(L)*TRZ
          FPRED6(8,L)=SQRT(SECANG(L))
!
! Fixed predictors for FWO sun mfbw = set7
          FPRED7(1,L)=SECANG(L)
          FPRED7(2,L)=SECANG(L)*SECANG(L)
          FPRED7(3,L)=SECANG(L)*TR
          FPRED7(4,L)=SECANG(L)*TJUNKS
          FPRED7(5,L)=TR
          FPRED7(6,L)=TJUNKS
          FPRED7(7,L)=SECANG(L)*TRZ
          FPRED7(8,L)=SQRT(SECANG(L))
!
!
! -----
! Ozone
! -----
          OJUNKA=SECANG(L)*A_O
!
! ozone predictors for FCOW = set4
          OPRED4(1,L)=OJUNKA
          OPRED4(2,L)=SQRT( OJUNKA )
          OPRED4(3,L)=OJUNKA*DT
!
! ozone predictors for FWO sun bfsw = set5
          OPRED5(1,L)=OJUNKA
!
! ozone predictors for FWO sun mfmw = set6
          OPRED6(1,L)=OJUNKA
!
! ozone predictors for FWO sun mfbw = set7
          OPRED7(1,L)=OJUNKA
!
!
! -----
! Water
! -----
          WJUNKA=SECANG(L)*A_W
          WJUNKR=SQRT( WJUNKA )
          WJUNKS=WJUNKA*WJUNKA
          WJUNKZ=WJUNKA*A_W/AZ_W
          WJUNK4=SQRT( WJUNKR )
!
! water predictors for FCOW = set4
          WPRED4( 1,L)=WJUNKA
          WPRED4( 2,L)=A_W
          WPRED4( 3,L)=WJUNKR
          WPRED4( 4,L)=WJUNKA*DT
          WPRED4( 5,L)=WJUNKS
          WPRED4( 6,L)=WJUNKR*DT
          WPRED4( 7,L)=WJUNK4
          WPRED4( 8,L)=WJUNKZ
          WPRED4( 9,L)=WJUNKA*SECANG(L)
          WPRED4(10,L)=WJUNKS*WJUNKA
          WPRED4(11,L)=WJUNKA*AZ_C*SECANG(L)
          WPRED4(12,L)=WJUNKZ/WJUNKR
          WPRED4(13,L)=WJUNKA*DT*SECANG(L)
!
! Water predictors for FWO sun bfsw = set5
          WPRED5( 1,L)=WJUNKA
          WPRED5( 2,L)=WJUNKA*WJUNKR
          WPRED5( 3,L)=WJUNKA*DT
!
! Water predictors for FWO sun mfmw = set6
          WPRED6( 1,L)=WJUNKA
          WPRED6( 2,L)=WJUNKA*WJUNKR
          WPRED6( 3,L)=WJUNKA*DT
          WPRED6( 4,L)=WJUNKS
          WPRED6( 5,L)=WJUNKA*WJUNKR*DT
          WPRED6( 6,L)=WJUNKA*WJUNKS
          WPRED6( 7,L)=WJUNKA*SECANG(L)
!
! Water predictors for FWO sun mfbw = set7
          WPRED7( 1,L)=WJUNKA
          WPRED7( 2,L)=WJUNKA*WJUNKR
          WPRED7( 3,L)=WJUNKA*DT
          WPRED7( 4,L)=WJUNKS
          WPRED7( 5,L)=WJUNKA*WJUNKR*DT
          WPRED7( 6,L)=WJUNKA*WJUNKS
          WPRED7( 7,L)=WJUNKA*SECANG(L)
          WPRED7( 8,L)=WJUNKZ
          WPRED7( 9,L)=WJUNKZ*WJUNKR
          WPRED7(10,L)=WJUNKA*WJUNK4
          WPRED7(11,L)=WJUNKA*WJUNKZ
          WPRED7(12,L)=WJUNKA*A_W
          WPRED7(13,L)=WJUNKS/WJUNK4
!
! ---------------
! Water continuum (for FWO, FOW, FMW, FCOW)
! ---------------
          CONPRD(1,L)=WJUNKA/TJUNKS
          CONPRD(2,L)=CONPRD(1,L)*A_W/TJUNKS
          CONPRD(3,L)=WJUNKA/TR
          CONPRD(4,L)=CONPRD(3,L)*A_W
          CONPRD(5,L)=CONPRD(1,L)*A_W
          CONPRD(6,L)=CONPRD(1,L)/TJUNKS
          CONPRD(7,L)=WJUNKA
!
!
! ---------------
! Carbon monoxide for FCOW = set4
! ---------------
          CJUNKA=SECANG(L)*A_C
          CJUNKR=SQRT( CJUNKA )
          CJUNKS=CJUNKA*CJUNKA
          CJUNKZ=CJUNKA*A_C/AZ_C
!
          CPRED4(1,L)=CJUNKA
          CPRED4(2,L)=CJUNKR
          CPRED4(3,L)=CJUNKA*DT
          CPRED4(4,L)=CJUNKS
          CPRED4(5,L)=CJUNKZ
          CPRED4(6,L)=CJUNKR*DT
          CPRED4(7,L)=SQRT( CJUNKR )
          CPRED4(8,L)=CJUNKZ/CJUNKR
          CPRED4(9,L)=A_C
          CPRED4(10,L)=CJUNKA*SECANG(L)
          CPRED4(11,L)=CJUNKR*SECANG(L)
!
!
! ---------------
! trace gas perturbation coefs
! ---------------
! The first 4 trace predictors are used by all trace gases
          TRCPRD(1,L)=SECANG(L)
          TRCPRD(2,L)=TR
          TRCPRD(3,L)=SECANG(L)*TR
          TRCPRD(4,L)=SECANG(L)*TJUNKS
! The last 3 trace predictors are only used by N2O
          TRCPRD(5,L)=SECANG(L)*SECANG(L)
          TRCPRD(6,L)=1.0
          TRCPRD(7,L)=SQRT( SECANG(L) )
!
       ENDDO
! End loop over layers
!
       RETURN
       END SUBROUTINE SUNPAR
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! TUNMLT
!
!F77====================================================================
!ROUTINE NAME:
! TUNMLT
!ABSTRACT:
! Apply transmittance tuning multipliers to the coefficients
!CALL PROTOCOL
! note: call is identical to RDCOEF
! TUNMLT ( AORB, IOUN, NCHAN, INDCHN, SETCHN,
! $ NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7,
! $ CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
! $ COEF1, COEF2, COEF3, COEF4, COEF5, COEF6, COEF7,
! $ NCHCO2, INDCO2, CLICO2, COFCO2,
! $ NCHN2O, INDN2O, CLIN2O, COFN2O,
! $ NCHSO2, INDSO2, CLISO2, COFSO2,
! $ NCHHNO, INDHNO, CLIHNO, COFHNO,
! $ NCHH2O, INDH2O, CLIH2O, WAZOP, WAVGOP, COFH2O,
! $ FREQ, COEFF, NCHNTE, CLISTN, COEFN, FX)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! CHAR*1 AORB specify A or B none
! INT arr CLICO2 CO2 channel list none
! INT arr CLIHNO HNO3 channel list none
! INT arr CLIH2O H2O channel list none
! INT arr CLIN2O N2O channel list none
! INT arr CLISO2 SO2 channel list none
! INT arr CLIST1 set1 channel list none
! INT arr CLIST2 set2 channel list none
! INT arr CLIST3 set3 channel list none
! INT arr CLIST4 set4 channel list none
! INT arr CLIST5 set5 channel list none
! INT arr CLIST6 set6 channel list none
! INT arr CLIST7 set7 channel list none
! INT arr CLISTN non-LTE channel list none
! REAL arr FREQ channel freqs cm-1
! REAL arr FX fixed gases adjustment none
! INT arr INDCHN indices of channels none
! INT arr INDCO2 CO2 pert channel indices none
! INT arr INDH2O OPTRAN H2O channel indices none
! INT arr INDHNO HNO3 pert channel indices none
! INT arr INDN2O N2O pert channel indices none
! INT arr INDSO2 SO2 pert channel indices none
! INTEGER IOUN I/O unit number none
! INTEGER NCHAN number of channels none
! INTEGER NCHN1 set1 number of channels none
! INTEGER NCHN2 set2 number of channels none
! INTEGER NCHN3 set3 number of channels none
! INTEGER NCHN4 set4 number of channels none
! INTEGER NCHN5 set5 number of channels none
! INTEGER NCHN6 set6 number of channels none
! INTEGER NCHN7 set7 number of channels none
! INTEGER NCHCO2 number of CO2 pert chans none
! INTEGER NCHH2O number of OPTRAN H2O chans none
! INTEGER NCHHNO number of HNO3 pert chans none
! INTEGER NCHN2O number of N2O pert chans none
! INTEGER NCHNTE non-LTE number of channels none
! INTEGER NCHSO2 number of SO2 pert chans none
! INT arr SETCHN set# (1-7) chan belongs to none (integer, 1 - 7)
! REAL arr WAZOP OPTRAN water grid kiloMoles/cm^2
! REAL arr WAVGOP OPTRAN water pred averges various
!INPUT/OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL arr COEF1 set1 fast trans coefs various
! REAL arr COEF2 set2 fast trans coefs various
! REAL arr COEF3 set3 fast trans coefs various
! REAL arr COEF4 set4 fast trans coefs various
! REAL arr COEF5 set5 fast trans coefs various
! REAL arr COEF6 set6 fast trans coefs various
! REAL arr COEF7 set7 fast trans coefs various
! REAL arr COEFF thermal "F" factor coefs various
! REAL arr COEFN non-LTE coefficients various
! REAL arr COFCO2 CO2 perturbation coefs various
! REAL arr COFH2O OPTRAN H2O coefs various
! REAL arr COFHNO HNO3 perturbation coefs various
! REAL arr COFN2O N2O perturbation coefs various
! REAL arr COFSO2 SO2 perturbation coefs various
!OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
! unit IOUN : input text file of tuning adjustments.
!COMMON BLOCKS
! none
!DESCRIPTION:
! August 2009 version of SARTA v1.08 code by L.Strow/S.Hannon.
!
! The routine reads a text file of tuning multipliers and
! applies them to the fast transmittace coefficients.
!
! The tuning multiplier file must consist of MXCHAN lines of data
! sorted (in ascending order) by channel ID and containing the
! following 9 columns:
! 1 2 3 4 5 6 7 8 9
! ID RJUNK XF XH2OL XH2OC XO3 XCO XCH4 XNTE
! where
! ID = integer, channel ID number
! RJUNK = real, value is ignored, eg perhaps channel freq
! XF,XH2OL,XH2OC,XO3,XCO,XCH4,XNTE = real, tuning multipler
! for fixed, water lines, water continuum, ozone, carbon
! monoxide, methane, and non-LTE. The value should
! be of the order of one.
! Comment lines may be included anywhere in the tuning multiplier
! by using a "!" or "%" as the first character on the line.
!
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 06 Feb 2003 Scott Hannon Created
! 11 Nov 2005 Scott Hannon Replace CO2pert tuning with non-LTE
! 08 Dec 2005 Scott Hannon Add SO2, HNO3, & N2O to arguments so
! tunmlt call is identical to rdcoef
! 09 May 2008 Scott Hannon Call changed to match 09 May 2008 rdcoef;
! add code for AORB filename
! 16 May 2008 Scott Hannon Add NJUNK to prevent adjusting 7th nonLTE
! coef
! 02 Nov 2008 Scott Hannon Bug fix for NJUNK
! 03 Aug 2009 Scott Hannon Add CLI<gas> for CO2,HNO,H2O,N2O,SO2 to
! call to match rdcoef; simplify OPTRAN
! adjustment using CLIH2O
!END====================================================================
! =================================================================
       SUBROUTINE TUNMLT ( AORB, IOUN, NCHAN, INDCHN, SETCHN, &
           NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7, &
          CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
           COEF1, COEF2, COEF3, COEF4, COEF5, COEF6, COEF7, &
          NCHCO2, INDCO2, CLICO2, COFCO2, &
          NCHN2O, INDN2O, CLIN2O, COFN2O, &
          NCHSO2, INDSO2, CLISO2, COFSO2, &
          NCHHNO, INDHNO, CLIHNO, COFHNO, &
          NCHH2O, INDH2O, CLIH2O, WAZOP, WAVGOP, COFH2O, &
            FREQ, COEFF, NCHNTE, CLISTN, COEFN, FX)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! INPUT
       CHARACTER*1 AORB
       INTEGER IOUN
       INTEGER NCHAN
       INTEGER INDCHN(MXCHAN)
       INTEGER SETCHN(MXCHAN)
       INTEGER NCHN1
       INTEGER NCHN2
       INTEGER NCHN3
       INTEGER NCHN4
       INTEGER NCHN5
       INTEGER NCHN6
       INTEGER NCHN7
       INTEGER CLIST1(MXCHN1)
       INTEGER CLIST2(MXCHN2)
       INTEGER CLIST3(MXCHN3)
       INTEGER CLIST4(MXCHN4)
       INTEGER CLIST5(MXCHN5)
       INTEGER CLIST6(MXCHN6)
       INTEGER CLIST7(MXCHN7)
       INTEGER NCHCO2
       INTEGER INDCO2(MXCHAN)
       INTEGER CLICO2(MXCHNC)
       INTEGER NCHN2O
       INTEGER INDN2O(MXCHAN)
       INTEGER CLIN2O(MXCHNN)
       INTEGER NCHSO2
       INTEGER INDSO2(MXCHAN)
       INTEGER CLISO2(MXCHNS)
       INTEGER NCHHNO
       INTEGER INDHNO(MXCHAN)
       INTEGER CLIHNO(MXCHNH)
       INTEGER NCHH2O
       INTEGER INDH2O(MXCHAN)
       INTEGER CLIH2O(MXCHNW)
       REAL WAZOP(MXOWLY)
       REAL WAVGOP(NOWAVG,MXOWLY)
       REAL FREQ(MXCHAN)
       REAL COEFF(NFCOEF,MXCHAN)
       INTEGER NCHNTE
       INTEGER CLISTN(MXCNTE)
       REAL FX(MAXLAY)
!
! INPUT/OUTPUT
       REAL COEF1(N1COEF,MAXLAY,MXCHN1)
       REAL COEF2(N2COEF,MAXLAY,MXCHN2)
       REAL COEF3(N3COEF,MAXLAY,MXCHN3)
       REAL COEF4(N4COEF,MAXLAY,MXCHN4)
       REAL COEF5(N5COEF,MAXLAY,MXCHN5)
       REAL COEF6(N6COEF,MAXLAY,MXCHN6)
       REAL COEF7(N7COEF,MAXLAY,MXCHN7)
       REAL COFCO2( NCO2,MAXLAY,MXCHNC)
       REAL COFN2O( NN2O,MAXLAY,MXCHNN)
       REAL COFSO2( NSO2,MAXLAY,MXCHNS)
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH)
       REAL COFH2O( NH2O,MXOWLY,MXCHNW)
       REAL COEFN(NNCOEF,MXCNTE)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       CHARACTER*80 FNTMLT
       CHARACTER*80 CLINE
       REAL RJUNK
       REAL XMULTS(7,MXCHAN)
       INTEGER I
       INTEGER ICHAN
       INTEGER IERR
       INTEGER J
       INTEGER L
       INTEGER NJUNK
       LOGICAL USEMLT(7)
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! Tuning multiplier filename
       IF (AORB .EQ. 'A') THEN
          FNTMLT=FATMLT
       ELSE
          FNTMLT=FBTMLT
       ENDIF
! -------------------------------
! Open the tuning multiplier file
! -------------------------------
       OPEN(UNIT=IOUN,FILE=FNTMLT,FORM='FORMATTED',STATUS='OLD', &
          IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(6,1020) IERR, FNTMLT
 1020 FORMAT('Error ',I5,' opening file:',/,A80)
          STOP
       ENDIF
!
! Initialize USEMLT
       DO J=1,7
          USEMLT(J)=.FALSE.
       ENDDO
! Initialize the channel counter
       I=0
! -------------------------
! Read the tuning mult file
! -------------------------
! Read a line of text from the file
 10 READ(IOUN,9000,END=910) CLINE
 9000 FORMAT(A80)
!
! Determine if the text line is data or a comment
       IF (CLINE(1:1) .NE. '!' .AND. CLINE(1:1) .NE. '%') THEN
!
! It's data, so increment the channel counter
          I=I+1
!
! Read the data from the text line
          READ(CLINE,*) ICHAN, RJUNK, (XMULTS(J,I),J=1,7)
!
! Check that ICHAN agrees with I
          IF (ICHAN .NE. I) THEN
             WRITE(6,1040) I, MXCHAN
 1040 FORMAT('Error reading tuning multipler file:',/, &
             'Expected channel ID ',I4,' but file has ID ',I4)
             STOP
          ENDIF
!
! Update USEMLT
          DO J=1,7
             IF (XMULTS(J,I) .LT. 0.9999 .OR. &
                 XMULTS(J,I) .GT. 1.0001) USEMLT(J)=.TRUE.
          ENDDO
!
       ENDIF
!
       GOTO 10
 910 CLOSE(IOUN)
!cc
! write(6,*) 'usemlt(fixed,H2Oline,H2Ocon,O3,CO,CH4,NTE)=',
! $ (USEMLT(J),J=1,7)
!cc
! ------------
! Adjust Fixed
! ------------
       IF (USEMLT(1)) THEN
! Set 1
          DO I=1,NCHN1
             ICHAN=CLIST1(I)
             DO L=1,MAXLAY
                DO J=8,15
                   COEF1(J,L,I)=XMULTS(1,ICHAN)*COEF1(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 2
          DO I=1,NCHN2
             ICHAN=CLIST2(I)
             DO L=1,MAXLAY
                DO J=8,15
                   COEF2(J,L,I)=XMULTS(1,ICHAN)*COEF2(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 3
          DO I=1,NCHN3
             ICHAN=CLIST3(I)
             DO L=1,MAXLAY
                DO J=8,15
                   COEF3(J,L,I)=XMULTS(1,ICHAN)*COEF3(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 4
          DO I=1,NCHN4
             ICHAN=CLIST4(I)
             DO L=1,MAXLAY
                DO J=8,18
                   COEF4(J,L,I)=XMULTS(1,ICHAN)*COEF4(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 5
          DO I=1,NCHN5
             ICHAN=CLIST5(I)
             DO L=1,MAXLAY
                DO J=8,18
                   COEF5(J,L,I)=XMULTS(1,ICHAN)*COEF5(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 6
          DO I=1,NCHN6
             ICHAN=CLIST6(I)
             DO L=1,MAXLAY
                DO J=8,15
                   COEF6(J,L,I)=XMULTS(1,ICHAN)*COEF6(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 7
          DO I=1,NCHN7
             ICHAN=CLIST7(I)
             DO L=1,MAXLAY
                DO J=8,15
                   COEF7(J,L,I)=XMULTS(1,ICHAN)*COEF7(J,L,I)
                ENDDO
             ENDDO
          ENDDO
       ENDIF ! end of Fixed
!
! ---------------
! Adjust H2O line
! ---------------
       IF (USEMLT(2)) THEN
! Set 1
          DO I=1,NCHN1
             ICHAN=CLIST1(I)
             DO L=1,MAXLAY
                DO J=16,26
                   COEF1(J,L,I)=XMULTS(2,ICHAN)*COEF1(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 2
          DO I=1,NCHN2
             ICHAN=CLIST2(I)
             DO L=1,MAXLAY
                DO J=26,36
                   COEF2(J,L,I)=XMULTS(2,ICHAN)*COEF2(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 3
          DO I=1,NCHN3
             ICHAN=CLIST3(I)
             DO L=1,MAXLAY
                DO J=25,35
                   COEF3(J,L,I)=XMULTS(2,ICHAN)*COEF3(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 4
          DO I=1,NCHN4
             ICHAN=CLIST4(I)
             DO L=1,MAXLAY
                DO J=33,45
                   COEF4(J,L,I)=XMULTS(2,ICHAN)*COEF4(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 5
          DO I=1,NCHN5
             ICHAN=CLIST5(I)
             DO L=1,MAXLAY
                DO J=19,21
                   COEF5(J,L,I)=XMULTS(2,ICHAN)*COEF5(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 6
          DO I=1,NCHN6
             ICHAN=CLIST6(I)
             DO L=1,MAXLAY
                DO J=16,22
                   COEF6(J,L,I)=XMULTS(2,ICHAN)*COEF6(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 7
          DO I=1,NCHN7
             ICHAN=CLIST7(I)
             DO L=1,MAXLAY
                DO J=16,28
                   COEF7(J,L,I)=XMULTS(2,ICHAN)*COEF7(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! OPTRAN
          DO I=1,NCHH2O
             ICHAN=CLIH2O(I)
             DO L=1,MXOWLY
                DO J=1,NH2O
                   COFH2O(J,L,I)=XMULTS(2,ICHAN)*COFH2O(J,L,I)
                ENDDO
             ENDDO
          ENDDO
!cc Old OPTRAN code
! DO ICHAN=1,MXCHAN
! IF (INDH2O(ICHAN) .GT. 0) THEN
! I=INDH2O(ICHAN)
! DO L=1,MXOWLY
! DO J=1,NH2O
! COFH2O(J,L,I)=XMULTS(2,ICHAN)*COFH2O(J,L,I)
! ENDDO
! ENDDO
! ENDIF
! ENDDO
!cc
       ENDIF ! end of H2O line
! --------------
! Adjust H2O con
! --------------
       IF (USEMLT(3)) THEN
! Set 1
          DO I=1,NCHN1
             ICHAN=CLIST1(I)
             DO L=1,MAXLAY
                DO J=1,7
                   COEF1(J,L,I)=XMULTS(3,ICHAN)*COEF1(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 2
          DO I=1,NCHN2
             ICHAN=CLIST2(I)
             DO L=1,MAXLAY
                DO J=1,7
                   COEF2(J,L,I)=XMULTS(3,ICHAN)*COEF2(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 3
          DO I=1,NCHN3
             ICHAN=CLIST3(I)
             DO L=1,MAXLAY
                DO J=1,7
                   COEF3(J,L,I)=XMULTS(3,ICHAN)*COEF3(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 4
          DO I=1,NCHN4
             ICHAN=CLIST4(I)
             DO L=1,MAXLAY
                DO J=1,7
                   COEF4(J,L,I)=XMULTS(3,ICHAN)*COEF4(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 5
          DO I=1,NCHN5
             ICHAN=CLIST5(I)
             DO L=1,MAXLAY
                DO J=1,7
                   COEF5(J,L,I)=XMULTS(3,ICHAN)*COEF5(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 6
          DO I=1,NCHN6
             ICHAN=CLIST6(I)
             DO L=1,MAXLAY
                DO J=1,7
                   COEF6(J,L,I)=XMULTS(3,ICHAN)*COEF6(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 7
          DO I=1,NCHN7
             ICHAN=CLIST7(I)
             DO L=1,MAXLAY
                DO J=1,7
                   COEF7(J,L,I)=XMULTS(3,ICHAN)*COEF7(J,L,I)
                ENDDO
             ENDDO
          ENDDO
       ENDIF ! end of H2O con
! ------------
! Adjust ozone
! ------------
       IF (USEMLT(4)) THEN
! Set 1
          DO I=1,NCHN1
             ICHAN=CLIST1(I)
             DO L=1,MAXLAY
                DO J=27,31
                   COEF1(J,L,I)=XMULTS(4,ICHAN)*COEF1(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 2
          DO I=1,NCHN2
             ICHAN=CLIST2(I)
             DO L=1,MAXLAY
                DO J=16,25
                   COEF2(J,L,I)=XMULTS(4,ICHAN)*COEF2(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 3: none
! Set 4
          DO I=1,NCHN4
             ICHAN=CLIST4(I)
             DO L=1,MAXLAY
                DO J=30,32
                   COEF4(J,L,I)=XMULTS(4,ICHAN)*COEF4(J,L,I)
                ENDDO
             ENDDO
          ENDDO
! Set 5
          J=22
          DO I=1,NCHN5
             ICHAN=CLIST5(I)
             DO L=1,MAXLAY
                COEF5(J,L,I)=XMULTS(4,ICHAN)*COEF5(J,L,I)
             ENDDO
          ENDDO
! Set 6
          J=23
          DO I=1,NCHN6
             ICHAN=CLIST6(I)
             DO L=1,MAXLAY
                COEF6(J,L,I)=XMULTS(4,ICHAN)*COEF6(J,L,I)
             ENDDO
          ENDDO
! Set 7
          J=29
          DO I=1,NCHN7
             ICHAN=CLIST7(I)
             DO L=1,MAXLAY
                COEF7(J,L,I)=XMULTS(4,ICHAN)*COEF7(J,L,I)
             ENDDO
          ENDDO
       ENDIF ! end of O3
! ---------
! Adjust CO
! ---------
       IF (USEMLT(5)) THEN
! Set 4
          DO I=1,NCHN4
             ICHAN=CLIST4(I)
             DO L=1,MAXLAY
                DO J=19,29
                   COEF4(J,L,I)=XMULTS(5,ICHAN)*COEF4(J,L,I)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
! --------------
! Adjust methane
! --------------
       IF (USEMLT(6)) THEN
! Set 3
          DO I=1,NCHN3
             ICHAN=CLIST3(I)
             DO L=1,MAXLAY
                DO J=16,24
                   COEF3(J,L,I)=XMULTS(6,ICHAN)*COEF3(J,L,I)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
! --------------
! Adjust non-LTE
! --------------
       IF (USEMLT(7)) THEN
! wrong NJUNK=MIN0(NCHNTE,6)
! wrong DO I=1,NJUNK ! do not adjust 7th coef
          NJUNK=MIN0(NNCOEF,6)
          DO I=1,NCHNTE
             ICHAN=CLISTN(I)
! wrong DO J=1,NNCOEF
             DO J=1,NJUNK
                COEFN(J,I)=XMULTS(7,ICHAN)*COEFN(J,I)
             ENDDO
          ENDDO
       ENDIF
!cc
!C ---------------
!C Adjust CO2 pert
!C ---------------
! IF (USEMLT(7)) THEN
! DO ICHAN=1,MXCHAN
! IF (INDCO2(ICHAN) .GT. 0) THEN
! I=INDCO2(ICHAN)
! DO L=1,MAXLAY
! DO J=1,NCO2
!C wrong COFCO2(J,L,I)=XMULTS(2,ICHAN)*COFCO2(J,L,I)
! COFCO2(J,L,I)=XMULTS(7,ICHAN)*COFCO2(J,L,I)
! ENDDO
! ENDDO
! ENDIF
! ENDDO
! ENDIF
!cc
!
       RETURN
       END SUBROUTINE TUNMLT
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!----------------------------------------------------------------------
! Utility functions & subroutines
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

       LOGICAL FUNCTION STR2BO(BUF)
!
! Translate "true"/"false" string to boolean
!
! include 'incFTC.f'
!
       CHARACTER*(*) BUF
!
       INTEGER I
       INTEGER J
       CHARACTER*1 XCHAR
! -------------------
!
       XCHAR=BUF(1:1)
       IF (XCHAR .EQ. 'T' .OR. XCHAR .EQ. 't') THEN
          STR2BO=.TRUE.
       ELSEIF (XCHAR .EQ. 'F' .OR. XCHAR .EQ. 'f') THEN
          STR2BO=.FALSE.
       ELSE
          I=LEN(BUF)
          J=MIN(I,5)
          WRITE(IOERR,1010) BUF(1:J)
 1010 FORMAT('Error, unexpected string boolean=',A5)
          STOP
       ENDIF
!
       RETURN
       END FUNCTION STR2BO
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!
!----------------------------------------------------------------------
!
       SUBROUTINE UPCASE(BUF)
!
! Convert a string to upper case
!
       CHARACTER*(*) BUF
!
       INTEGER I
       INTEGER IC
       INTEGER J
! -----------------
!
       DO I=1,LEN(BUF)
          IC=ICHAR( BUF(I:I) )
          IF (IC .GE. ICHAR('a') .AND. IC .LE. ICHAR('z')) THEN
             J=IC + (ICHAR('A') - ICHAR('a'))
             BUF(I:I)=CHAR(J)
          ENDIF
       ENDDO
!
       RETURN
       END SUBROUTINE UPCASE
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!
!----------------------------------------------------------------------
!
       INTEGER FUNCTION LENNB(BUF)
!
! Find the index of the last non-blank char in a string
!
       CHARACTER*(*) BUF
!
       INTEGER I
! -----------------
!
       I=LEN(BUF)
 10 IF (I .GT. 0) THEN
          IF (BUF(I:I) .EQ. ' ') THEN
             I=I - 1
             GOTO 10
          ENDIF
       ENDIF
       LENNB=I
!
       RETURN
       END FUNCTION LENNB
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!C----------------------------------------------------------------------
!
       INTEGER FUNCTION LENDSH(BUF)
!
! Find the index of the last dash '/' char in a string, or
! zero if no dash is found.
!
       CHARACTER*(*) BUF
!
       INTEGER I
! -----------------
!
       I=LEN(BUF)
 10 IF (I .GT. 0) THEN
          IF (BUF(I:I) .EQ. '/') THEN
             I=I - 1
             GOTO 10
          ENDIF
       ENDIF
       LENDSH=I
!
       RETURN
       END FUNCTION LENDSH
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!
!----------------------------------------------------------------------
!
       SUBROUTINE N2BITS(NUMBER,LFLAGS)
!
! Converts an integer number into a logical(32) array
!
! Input:
       INTEGER*4 NUMBER
!
! Ouput:
       LOGICAL LFLAGS(32)
!
! Local:
       INTEGER N
       INTEGER I2NM1
       INTEGER IWORK
! -----------------
!
       IWORK=NUMBER
! Note: 1st bit which is for +- sign
       IF (IWORK .LT. 0) THEN
          IWORK=ABS(IWORK)
          LFLAGS(32)=.TRUE.
       ELSE
          LFLAGS(32)=.FALSE.
       ENDIF
       DO N=31,1,-1
          I2NM1=2**(N-1)
          IF (IWORK .GE. I2NM1) THEN
             LFLAGS(N)=.TRUE.
             IWORK=IWORK - I2NM1
          ELSE
             LFLAGS(N)=.FALSE.
          ENDIF
       ENDDO
!
       RETURN
       END SUBROUTINE N2BITS
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!
!-----------------------------------------------------------------------
!
       SUBROUTINE BITS2N(NUMBER,LFLAGS)
!
! Converts a logical(32) array into an integer number
!
! Input:
       LOGICAL LFLAGS(32)
!
! Ouput:
       INTEGER*4 NUMBER
!
! Local:
       INTEGER N
       INTEGER I2NM1
       INTEGER IWORK
! ---------------
!
       IWORK=0
! Note: ignore 1st bit which is for +- sign
       DO N=1,31
          I2NM1=2**(N-1)
          IF (LFLAGS(N)) THEN
             IWORK=IWORK + I2NM1
          ENDIF
       ENDDO
       IF (LFLAGS(32)) IWORK=-IWORK
!
       NUMBER=IWORK
!
       RETURN
       END SUBROUTINE BITS2N
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!
!-----------------------------------------------------------------------
!

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!
!----------------------------------------------------------------------
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore Country (UMBC)
!
! AIRS
!
! VACONV
!
!F77====================================================================
!ROUTINE NAME:
! VACONV (real function)
!ABSTRACT:
! Function to convert the AIRS satellite viewing angle into the
! local path angle.
!CALL PROTOCOL
! VACONV( SVA, SALT, ALT )
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL SVA Satellite viewing angle degrees
! REAL SALT Satellite altitude kilometers
! REAL ALT Average layer altitude meters
!OUTPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! REAL fun VACONV local path angle radians
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S):
! USEFAST
!ROUTINES CALLED:
! none
!FILES ACCESSED:
! none
!COMMON BLOCKS
! none
!DESCRIPTION:
! March 1998 version of the 100 layer AIRS Fast Transmittance
! Code by L.Strow/S.Hannon.
!
! ===================================================================
! Function to convert the AIRS satellite viewing angle into a local
! path angle. The local path angle generally varies slightly with
! altitude due to the curvature of the Earth and its atmosphere.
! The effect is largest at the maximum satellite viewing angle,
! and goes to zero as the viewing angle approaches 0 degrees.
!
! For AIRS, the maximum difference in the local path angle secant
! between the bottom and top layers is around 3 percent.
!
! Currently this function only considers the geometry of the
! situation, and no refractive effects are included.
!
! The layers of the atmosphere may be considered as concentric
! rings with some average altitude. A ray traced thru these rings
! at any viewing angle other than nadir will have a slightly
! different angle (relative to the outward radial at the point
! of intersection) in each ring.
!
! If the Earth is treated as a perfect sphere of radius RE (hard
! coded into this routine), then the local angle may be calculated
! using trigonometry if we know:
! The satellite viewing angle
! The satellite's altitude above the Earth's surface
! The layer's altitude above the Earth's surface.
!
! The solution uses the law of sines and sin(180 - x) = sin(x)
! ===================================================================
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! No refractive effects have been included.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- ----------------------------------------
! 10 Apr 1995 Scott Hannon Created
! 1 Apr 1997 Scott Hannon Fix error in Earth radius (was 3.67E+3)
! 27 Feb 1998 Scott Hannon Simplified; made SALT an input var
!END====================================================================
! =================================================================
       REAL FUNCTION VACONV( LAT, SVA, SALT, ALT )
! =================================================================
!
! Viewing Angle CONVersion
!
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
       REAL LAT
       REAL SVA
       REAL SALT
       REAL ALT
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       REAL CONV
       REAL RE
       REAL RS
       REAL RA
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
! ------------------
! Assign some values
! ------------------
! CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02
!
! RE = radius of the Earth (in km)
!       RE=6.37E+03

! More precise parameter for earth radius. (Toshi Matsui Oct 2013.)
    call earth_radius_single( LAT*const_degrad, RE )


!
! RA = radius of the point to calc the angle at (in km)
       RA=RE + (ALT/1000.0)
!
! RS = radius of the satellite orbit (in km)
       RS=RE + SALT
!
! -----------------
! Do the conversion
! -----------------
!
       VACONV=ASIN( (RS/RA) * SIN(CONV*SVA) )
!
       RETURN
       END FUNCTION VACONV
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

       FUNCTION C4TOR4(C4IN) RESULT(R4)
!
!      Convert four characters into a REAL*4 number
!
!      Input:
       CHARACTER(len=1) :: C4IN(4)
       REAL :: R4
!
!      Local:

!       UNION
!          MAP
             !CHARACTER*4 C4
!          ENDMAP
!          MAP
!             REAL R4
!          ENDMAP
!       ENDUNION
!      Note: the Absoft compiler requires the union be done on
!      local variables; it crashes if tried on the input argument
!
!       C4=C4IN
!       C4TOR4=R4
!
!
       RETURN
       END FUNCTION C4TOR4



!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County [UMBC]
!
! AIRS
!
! WRTRTP
!
!F77====================================================================
!ROUTINE NAME: WRTRTP
!ABSTRACT:
! Write a profile to a previously openned RTP file
!CALL PROTOCOL:
! WRTRTP(IP, IOPCO, NCHAN, RAD, PROF)
!INPUT PARAMETERS:
! type name purpose units
! -------- ------ -------------------------- ---------------------
! INTEGER IP profile count so far none
! INTEGER IOPCO input RTP file I/O number none
! INTEGER NCHAN # of channels none
! REAL arr RAD radiance W/m^2/cm^-1/sterad
! STRUCT PROF RTP profile structure (see attributes)
!OUTPUT PARAMETERS:
! none
!INPUT/OUTPUT PARAMETERS: none
!RETURN VALUES: none
!PARENT(S): KLAYERS
!ROUTINES CALLED: none
!FILES ACCESSED:
! Output RTP file with I/O number IOPCO
! unit IOERR: error message
!COMMON BLOCKS: none
!DESCRIPTION:
! Writes a single profile to a previously openned RTP file.
!ALGORITHM REFERENCES: see DESCRIPTION
!KNOWN BUGS AND LIMITATIONS:
! none
!ROUTINE HISTORY:
! Date Programmer Comments
!------------ ----------------- ----------------------------------------
! 14 Feb 2001 Scott Hannon created based on version for klayers
! 23 Oct 2008 Scott Hannon Minor update for rtpV201
!END====================================================================
! =================================================================
       SUBROUTINE WRTRTP(IP, IOPCO, NCHAN, RAD, PROF)
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
! include 'rtpdefs.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
! none
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! Input parameters:
       INTEGER IP
       INTEGER IOPCO
       INTEGER NCHAN
       REAL RAD(MXCHAN)
!
! Profile data structure
! RECORD /RTPPROF/ PROF
       type( RTPPROF ) PROF
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER ISTAT
       INTEGER rtpwrite
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE begins below
!***********************************************************************
!***********************************************************************
! -------------------------
! Load up the new PROF data
! -------------------------
! Loop over the channels
       DO I=1,NCHAN
! Convert from Watts/m^2/cm^-1 to milliWatts/m^2/cm^-1
          PROF%rcalc(I)=RAD(I)*1000.0
       ENDDO
!
! -------------------------
! Write the current profile
! -------------------------
! ISTAT=rtpwrite(IOPCO, PROF)
!
       IF (ISTAT .EQ. -1) THEN
          WRITE(IOERR,1010) IP
 1010 FORMAT('ERROR! unable to write PROF data for prof ',I5)
       ENDIF
!
       RETURN
       END SUBROUTINE WRTRTP
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
! This version with rtp input profiles and command-line arguments
!=======================================================================
!=======================================================================
!
! University of Maryland Baltimore County (UMBC)
!
! AIRS
!
! SARTA version with trace gases
!
!F77====================================================================
!ROUTINE NAME:
! SARTA
!ABSTRACT:
! Program to quickly compute simulated AIRS radiances.
!CALL PROTOCOL
! none (main program)
!INPUT PARAMETERS:
! none
!OUTPUT PARAMETERS:
! none
!INPUT/OUTPUT PARAMETERS:
! none
!RETURN VALUES:
! none
!PARENT(S)
! none
!ROUTINES CALLED:
! CALOWP : calc OPTRAN water predictors
! CALPAR : calculate a profile's predictors
! CALRAD : calc radiance
! CALT1 : calc effective layer trans for set1 (FWO)
! CALT2 : calc effective layer trans for set2 (FOW)
! CALT3 : calc effective layer trans for set3 (FMW)
! CALT4 : calc effective layer trans for set4 (FCOW)
! CALT5 : calc effective layer trans for set5 (FWO bfsw)
! CALT6 : calc effective layer trans for set6 (FWO mfmw)
! CALT7 : calc effective layer trans for set7 (FWO mfbw)
! FAKETZ : calc a "fake" (rough approx) surface-to-space trans
! RDCOEF : read the fast transmittance coefficients
! RDPROF : read profile data
! RDSUN : read the solar radiance datafile
! SUNPAR : calc a profile's predictors for sets4-7 (sun channels)
! CALNTE : calc radiance contribution for non-LTE
!FILES ACCESSED:
! incFTC.f : include file of parameter statements accessed during
! compilation only.
! unit IOUN: used by routines RDCOEF and RDPROF.
! unit 6: USEFAST text messages to the screen
! unit 5: USEFAST user input instructions, etc
! unit 10: USEFAST output radiance, text file(s)
!COMMON BLOCKS
! COMLEV : layer boundary pressure levels
!DESCRIPTION:
! May 2001 version of the SARTA_RTP (Stand-Alone Rapid
! Transmittance Algorith with RTP I/O) by
! L.L.Strow, S.Hannon, and H.Mottler
!
! Computes radiances for the layers profiles contained in the
! input RTP file. This is the main program, and consists
! primarily of calls to the external routines to do most of the
! computing.
!ALGORITHM REFERENCES:
! none
!KNOWN BUGS AND LIMITATIONS:
! This program is only intended as a demo of the fast model.
!ROUTINE HISTORY:
! Date Programmer Comments
! ----------- -------------- -------------------------------------------
! 01 Dec 1994 Scott Hannon Created
! 10 Apr 1995 Scott Hannon New header comments; added ALT; new
! external function VACONV; SECANG may
! vary with layer
! 03 Jul 1995 Scott Hannon Add parameter DZ/RDZ to RDPROF call
! 03 Feb 1997 Scott Hannon Re-written for FWO+FOW+FMW+FCOW
! 12 Sep 1997 Scott Hannon Re-written for 7 sets and reflected sun
! and downwelling thermal
! 30 Sep 1997 Scott Hannon Added variable CO2
! 27 Feb 1998 Scott Hannon Added OPTRAN water
! 26 Aug 1998 Scott Hannon Added LBOT to calls to CALPAR, CALOWP,
! and SUNPAR; rename TBOT to TSURF; calc
! fractional bottom layer temperature and
! put it in TEMP(LBOT)
! 15 Oct 1999 Scott Hannon Add ANGMAX and re-arrange angle conv
! 31 Mar 2000 Scott Hannon Redid calpar for FIXMUL and added getbot
! 15 Mar 2001 Scott Hannon Major re-write for RTP
! 03 May 2001 Scott Hannon Add COMLEV; add PLEV to getbot call
! 13 Sep 2001 Scott Hannon Changes to check of FCHAN vs FREQ
! 01 Nov 2002 Scott Hannon Added SATZEN & SALT to RDRTP call, and
! if valid use SATZEN rather than SATANG
! 03 Jan 2003 Scott Hannon Delete SUNSEC, add XZ & SUNFDG & code
! to fudge large sun angles (previously
! sunang>80 were treated as no sun).
! 24 Jul 2003 Scott Hannon Fix error in TEMP(LBOT) calc for
! bottom fractional layer; add PLAY
! 06 Feb 2004 Scott Hannon Add call to TUNMLT; add call to MEAN_T
! and associated prep code; add PTYPE
! to OPNRTP call; add LRHOT to RDINFO,
! OPNRTP, & SETEMS calls.
! 20 Dec 2004 Scott Hanonn Add NLAY to getbot.f call; add PTYPE
! to rdrtp_so2.f call
! 18 May 2005 Scott Hannon Add HNO3 based on SO2 code
! 28 Jun 2005 Scott Hannon "trace" version for CO2,SO2,HNO3,N2O
! 13 Oct 2005 Scott Hannon Add non-LTE
! 08 Dec 2005 Scott Hannon Update tunmlt call for non-LTE tuning
! 02 May 2007 Scott Hannon Replace hardcoded default SALT value
! with XSALT from incFTC.
! 23 Jan 2008 Scott Hannon Add LCO2PM to allow CO2 profile in ppmv;
! add LCO2,LN2O,LSO2,LHNO3 switches
! 31 Jan 2008 Scott Hannon Fix bug so LCO2,LN2O,LSO2,LHNO3 are LOGICAL
! as intended instead of REAL
! 05 May 2008 Scott Hannon Set SUNCOS=0 when DOSUN=false
! 09 May 2008 Scott Hannon Add freq adjustment code including A & B
! variables, two calls to rdcoef & tunmlt,
! new call to sumcof, and add yoff to rdrtp.
! 14 May 2008 Scott Hannon Add default freq FRQDEF and WGTA,WGTB; add
! CO2TOP to calpar.f & calnte.f calls;
! add CO2PPM to calpar.f call; move no
! prof CO2MLT calc to calpar.f; add AORB
! to RDSUN
! 24 Oct 2008 Scott Hannon Update for rtpV201
! 04 Aug 2009 Scott Hannon Add call to SPECCAL; add DFCAL; add
! CLI<gas> for CO2,HNO,H2O,N2O,SO2 to
! RDCOEF & TUNMLT calls; add CLI<*> to
! SUMCOF call; adjust head.vchan
! 06 Aug 2009 Scott Hannon Remove NCHAN and INDCHAN from SPECCAL call
! 11 Oct 2011 Scott Hannon Replace call to speccal with rddyab; add
! call for setyof; modify arguments of
! rdrtp and sumcof; change HEAD%vchan
! 21 Oct 2011 Scott Hannon Bug fix: SUMCOF call uses YOFCHN (was YOFF)
!END====================================================================
! =================================================================
 SUBROUTINE SARTA_DRIVER( proc )
 use module_mpi
! =================================================================
!-----------------------------------------------------------------------
! IMPLICIT NONE
!-----------------------------------------------------------------------
 IMPLICIT NONE
 character(len=*),intent(in) :: proc ! init or calc
 SAVE  !all parameters are saved.

!-----------------------------------------------------------------------
! INCLUDE FILES
!-----------------------------------------------------------------------
! include 'incFTC.f'
! include 'rtpdefs.f'
!-----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
       !REAL VACONV
       !REAL SACONV
!-----------------------------------------------------------------------
! ARGUMENTS
!-----------------------------------------------------------------------
! none (main program)
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
!
       INTEGER IOUN ! I/O unit number
!
! for RDINFO
       CHARACTER*80 FIN ! input RTP filename
       CHARACTER*80 FOUT ! output RTP filename
       LOGICAL LRHOT ! force refl therm rho=(1-emis)/pi?
!       INTEGER NWANTP ! number of wanted profiles (-1=all)
!       INTEGER LISTP(MAXPRO) ! list of wanted profiles
!
! for OPNRTP
       INTEGER PTYPE ! profile type
       INTEGER NCHAN ! # of selected channels
       REAL FCHAN(MXCHAN) ! chan center frequency
       INTEGER LSTCHN(MXCHAN) ! list of selected channels
       INTEGER INDCHN(MXCHAN) ! array indices for all channels
       INTEGER IH2O ! index of H2O in gamnt
       INTEGER IO3 ! index of O3 in gamnt
       INTEGER ICO ! index of CO in gamnt
       INTEGER ICH4 ! index of CH4 in gamnt
       INTEGER ICO2 ! index of CO2 in gamnt
       INTEGER ISO2 ! index of SO2 in gamnt
       INTEGER IHNO3 ! index of HNO3 in gamnt
       INTEGER IN2O ! index of N2O in gamnt
       INTEGER IOPCO ! output RTP unit
       LOGICAL LCO2PM ! CO2 profile in ppmv?
!
! for RDCOEF ! Info for selected channels only
       CHARACTER*1 AORB ! A or B
       INTEGER SETCHN(MXCHAN) ! set # for each channel
       INTEGER NCHN1 ! # of set1 channels
       INTEGER NCHN2 ! # of set2 channels
       INTEGER NCHN3 ! # of set3 channels
       INTEGER NCHN4 ! # of set4 channels
       INTEGER NCHN5 ! # of set5 channels
       INTEGER NCHN6 ! # of set6 channels
       INTEGER NCHN7 ! # of set7 channels
       INTEGER CLIST1(MXCHN1) ! list of set1 channels
       INTEGER CLIST2(MXCHN2) ! list of set2 channels
       INTEGER CLIST3(MXCHN3) ! list of set3 channels
       INTEGER CLIST4(MXCHN4) ! list of set4 channels
       INTEGER CLIST5(MXCHN5) ! list of set5 channels
       INTEGER CLIST6(MXCHN6) ! list of set6 channels
       INTEGER CLIST7(MXCHN7) ! list of set7 channels
       INTEGER CLICO2(MXCHNC) ! list of CO2 channels
       INTEGER CLIHNO(MXCHNH) ! list of HNO3 channels
       INTEGER CLIH2O(MXCHNW) ! list of H2O channels
       INTEGER CLIN2O(MXCHNN) ! list of N2O channels
       INTEGER CLISO2(MXCHNS) ! list of SO2 channels
       INTEGER LABOVE(MXCHAN) ! chan downwelling thermal layer above
       REAL FREQ(MXCHAN) ! chan center frequency
       REAL COEF1(N1COEF,MAXLAY,MXCHN1) ! coefs for set1 chans
       REAL COEF2(N2COEF,MAXLAY,MXCHN2) ! coefs for set2 chans
       REAL COEF3(N3COEF,MAXLAY,MXCHN3) ! coefs for set3 chans
       REAL COEF4(N4COEF,MAXLAY,MXCHN4) ! coefs for set4 chans
       REAL COEF5(N5COEF,MAXLAY,MXCHN5) ! coefs for set5 chans
       REAL COEF6(N6COEF,MAXLAY,MXCHN6) ! coefs for set6 chans
       REAL COEF7(N7COEF,MAXLAY,MXCHN7) ! coefs for set7 chans
       REAL COEFF(NFCOEF,MXCHAN) ! coefs for chan "F" factor
       INTEGER NCHCO2 ! number of CO2 chans
       INTEGER NCHN2O ! number of N2O chans
       INTEGER NCHSO2 ! number of SO2 chans
       INTEGER NCHHNO ! number of HNO3 chans
       INTEGER INDCO2(MXCHAN) ! chan indices for CO2 pert
       INTEGER INDN2O(MXCHAN) ! chan indices for N2O pert
       INTEGER INDSO2(MXCHAN) ! chan indices for SO2 pert
       INTEGER INDHNO(MXCHAN) ! chan indices for HNO3 pert
       REAL COFCO2( NCO2,MAXLAY,MXCHNC) ! coefs for CO2 pert
       REAL COFN2O( NN2O,MAXLAY,MXCHNN) ! coefs for N2O pert
       REAL COFSO2( NSO2,MAXLAY,MXCHNS) ! coefs for SO2 pert
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH) ! coefs for HNO3 pert
       INTEGER NCHH2O ! number of OPTRAN H2O chans
       INTEGER INDH2O(MXCHAN) ! chan indices for OPTRAN H2O
       REAL WAZOP(MXOWLY) ! OPTRAN water l-to-s amounts
       REAL WAVGOP(NOWAVG,MXOWLY) ! OPTRAN raw predictor averages
       REAL COFH2O( NH2O,MXOWLY,MXCHNW) ! coefs for OPTRAN H2O
       REAL FX(MAXLAY) ! fixed gases adjustment
       INTEGER NCHNTE ! number of non-LTE channels
       INTEGER CLISTN(MXCNTE) ! non-LTE channel list
       REAL COEFN(NNCOEF,MXCNTE) ! non-LTE coefficients
!
! database "A"
       REAL FREQA(MXCHAN) ! chan center frequency
       REAL HSUNA(MXCHAN) ! solar irradiance
       REAL COEFFA(NFCOEF,MXCHAN) ! coefs for chan "F" factor
       REAL COEFNA(NNCOEF,MXCNTE) ! non-LTE coefficients
       REAL COEF1A(N1COEF,MAXLAY,MXCHN1) ! coefs for set1 chans
       REAL COEF2A(N2COEF,MAXLAY,MXCHN2) ! coefs for set2 chans
       REAL COEF3A(N3COEF,MAXLAY,MXCHN3) ! coefs for set3 chans
       REAL COEF4A(N4COEF,MAXLAY,MXCHN4) ! coefs for set4 chans
       REAL COEF5A(N5COEF,MAXLAY,MXCHN5) ! coefs for set5 chans
       REAL COEF6A(N6COEF,MAXLAY,MXCHN6) ! coefs for set6 chans
       REAL COEF7A(N7COEF,MAXLAY,MXCHN7) ! coefs for set7 chans
       REAL COFCOA( NCO2,MAXLAY,MXCHNC) ! coefs for CO2 pert
       REAL COFN2A( NN2O,MAXLAY,MXCHNN) ! coefs for N2O pert
       REAL COFSOA( NSO2,MAXLAY,MXCHNS) ! coefs for SO2 pert
       REAL COFHNA( NHNO3,MAXLAY,MXCHNH) ! coefs for HNO3 pert
       REAL COFH2A( NH2O,MXOWLY,MXCHNW) ! coefs for OPTRAN H2O
!
! database "B"
       REAL FREQB(MXCHAN) ! chan center frequency
       REAL HSUNB(MXCHAN) ! solar irradiance
       REAL COEFFB(NFCOEF,MXCHAN) ! coefs for chan "F" factor
       REAL COEFNB(NNCOEF,MXCNTE) ! non-LTE coefficients
       REAL COEF1B(N1COEF,MAXLAY,MXCHN1) ! coefs for set1 chans
       REAL COEF2B(N2COEF,MAXLAY,MXCHN2) ! coefs for set2 chans
       REAL COEF3B(N3COEF,MAXLAY,MXCHN3) ! coefs for set3 chans
       REAL COEF4B(N4COEF,MAXLAY,MXCHN4) ! coefs for set4 chans
       REAL COEF5B(N5COEF,MAXLAY,MXCHN5) ! coefs for set5 chans
       REAL COEF6B(N6COEF,MAXLAY,MXCHN6) ! coefs for set6 chans
       REAL COEF7B(N7COEF,MAXLAY,MXCHN7) ! coefs for set7 chans
       REAL COFCOB( NCO2,MAXLAY,MXCHNC) ! coefs for CO2 pert
       REAL COFN2B( NN2O,MAXLAY,MXCHNN) ! coefs for N2O pert
       REAL COFSOB( NSO2,MAXLAY,MXCHNS) ! coefs for SO2 pert
       REAL COFHNB( NHNO3,MAXLAY,MXCHNH) ! coefs for HNO3 pert
       REAL COFH2B( NH2O,MXOWLY,MXCHNW) ! coefs for OPTRAN H2O
!
! for RDDYAB
       DOUBLE PRECISION ETIME(MAXABW)
       REAL FREQAB(MXCHAN)
       REAL WIDTH(MXCHAN)
       INTEGER INDMOD(MXCHAN)
       INTEGER IABWGT(MXCHAN,MAXABW)
       REAL DELYAB(MXCHAN,3)
!
! for rtpopen
       INTEGER rtpopen
       CHARACTER*1 MODE
!
! for FAKETZ
       INTEGER NFAKE ! # of channels to "fake"
       INTEGER INDFAK(MXCHAN) ! indices of channels to fake
!
! for RDPROF; reference profile
       CHARACTER*40 RPNAM ! ref prof name/ID
       REAL RALT(MAXLAY) ! ref prof layer altitude
       REAL RDZ(MAXLAY) ! ref prof layer thickness
       REAL RPRES(MAXLAY) ! ref prof layer average pressure
       REAL RTEMP(MAXLAY) ! ref prof layer average temperature
       REAL RFAMNT(MAXLAY) ! ref prof layer "fixed" (CO2) amount
       REAL RWAMNT(MAXLAY) ! ref prof layer water (H2O) amount
       REAL ROAMNT(MAXLAY) ! ref prof layer ozone (O3) amount
       REAL RCAMNT(MAXLAY) ! ref prof layer carbon monoxide (CO) amount
       REAL RMAMNT(MAXLAY) ! ref prof layer methane (CH4) amount
       REAL RSAMNT(MAXLAY) ! ref prof layer sulfer dioxide (SO2) amount
       REAL RHAMNT(MAXLAY) ! ref prof layer nitric acid (HNO3) amount
       REAL RNAMNT(MAXLAY) ! ref prof layer nitrous oxide (N2O) amount
!
! for RDRTP; profile to calculate
       INTEGER NLAY ! number of layers in profile
       REAL LAT ! prof latitude
       REAL LON ! prof longitude
       REAL ALT(MAXLAY) ! prof layer altitudes
       REAL TEMP(MAXLAY) ! prof layer average temperature
       REAL WAMNT(MAXLAY) ! prof layer water (H2O) amount
       REAL OAMNT(MAXLAY) ! prof layer ozone (O3) amount
       REAL CAMNT(MAXLAY) ! prof layer carbon monoxide (CO) amount
       REAL MAMNT(MAXLAY) ! prof layer methane (CH4) amount
       REAL FAMNT(MAXLAY) ! prof layer CO2 amount
       REAL SAMNT(MAXLAY) ! prof layer SO2 amount
       REAL HAMNT(MAXLAY) ! prof layer HNO3 amount
       REAL NAMNT(MAXLAY) ! prof layer N2O amount
!
! for surface
       INTEGER LBOT ! bottom layer index number
       INTEGER NEMIS ! # of emis pts
       REAL PSURF ! surface pressure
       REAL BLMULT ! bottom layer fractional multiplier
       REAL FEMIS(MXEMIS) ! emis freq pts
       REAL XEMIS(MXEMIS) ! emis pts
       REAL XRHO(MXEMIS) ! reflec pts
!
! for MEAN_T
       REAL TPSEUD(MAXLAY)
!
! for CALPAR
       LOGICAL LCO2 ! CO2 profile switch
       LOGICAL LN2O ! N2O profile switch
       LOGICAL LSO2 ! SO2 profile switch
       LOGICAL LHNO3 ! HNO3 profile switch
       REAL SECANG(MAXLAY) ! local path angle secant
       REAL FIXMUL(MAXLAY) ! "fixed" amount multiplier (~1)
       REAL CONPRD( N1CON,MAXLAY) ! water continuum predictors
       REAL FPRED1( N1FIX,MAXLAY) ! set1 "fixed" predictors
       REAL FPRED2( N2FIX,MAXLAY) ! set2 "fixed" predictors
       REAL FPRED3( N3FIX,MAXLAY) ! set3 "fixed" predictors
       REAL FPRED4( N4FIX,MAXLAY) ! set4 "fixed" predictors
       REAL FPRED5( N5FIX,MAXLAY) ! set5 "fixed" predictors
       REAL FPRED6( N6FIX,MAXLAY) ! set6 "fixed" predictors
       REAL FPRED7( N7FIX,MAXLAY) ! set7 "fixed" predictors
       REAL WPRED1( N1H2O,MAXLAY) ! set1 water predictors
       REAL WPRED2( N2H2O,MAXLAY) ! set2 water predictors
       REAL WPRED3( N3H2O,MAXLAY) ! set3 water predictors
       REAL WPRED4( N4H2O,MAXLAY) ! set4 water predictors
       REAL WPRED5( N5H2O,MAXLAY) ! set5 water predictors
       REAL WPRED6( N6H2O,MAXLAY) ! set6 water predictors
       REAL WPRED7( N7H2O,MAXLAY) ! set7 water predictors
       REAL OPRED1( N1O3,MAXLAY) ! set1 ozone predictors
       REAL OPRED2( N2O3,MAXLAY) ! set2 ozone predictors
       REAL OPRED4( N4O3,MAXLAY) ! set4 ozone predictors
       REAL OPRED5( N5O3,MAXLAY) ! set5 ozone predictors
       REAL OPRED6( N6O3,MAXLAY) ! set6 ozone predictors
       REAL OPRED7( N7O3,MAXLAY) ! set7 ozone predictors
       REAL MPRED3( N3CH4,MAXLAY) ! set3 methane predictors
       REAL CPRED4( N4CO,MAXLAY) ! set4 carbon monoxide predictors
       REAL TRCPRD(NTRACE,MAXLAY) ! trace gas pert perdictors
       REAL CO2MLT(MAXLAY) ! CO2 perturbation multiplier
       REAL SO2MLT(MAXLAY) ! SO2 perturbation multiplier
       REAL HNOMLT(MAXLAY) ! HNO3 perturbation multiplier
       REAL N2OMLT(MAXLAY) ! N2O perturbation multiplier
       REAL CO2TOP ! top layers CO2 mixing ratio
!
! for CALOWP
       REAL WAANG(MAXLAY)
       INTEGER LOPMIN
       INTEGER LOPMAX
       REAL H2OPRD( NH2O,MXOWLY)
       LOGICAL LOPUSE(MXOWLY)
       INTEGER LOPLOW(MAXLAY)
       REAL DAOP(MAXLAY)
!
! for CALT
       REAL TAU(MAXLAY,MXCHAN) ! chan layer effective trans
       REAL TAUZ(MXCHAN) ! chan surface-to-space trans
       REAL WAOP(MXOWLY) ! OPTRAN abs coef scaling factor
       REAL XZ ! optical depth multiplier for TAUZ
       LOGICAL LTAU ! Calc all layer transmittances?
!
! for CALRAD
       REAL TSURF ! surface temperature
       REAL EMIS(MXCHAN) ! chan surface emissivity
       REAL RHOSUN(MXCHAN) ! chan reflectivity for sun
       REAL RHOTHR(MXCHAN) ! chan reflectivity for downwelling thermal
       REAL RAD(MXCHAN) ! chan radiance
       REAL BT(MXCHAN) ! chan brightness temperature
!
! for RDSUN
       REAL HSUN(MXCHAN) ! sun radiance (direct from sun)
!
! Other variables for the sun
       REAL SUNANG ! solar zenith angle (at 0 altitude)
       REAL SZALAY ! solar zenith angle in some layer
       REAL SUNCOS ! cosine of sun zenith angle
       REAL SCOS1 ! cosine of sun zenith angle at layer1
       REAL SUNFDG ! fudge factor for large solar angles
       REAL SECSUN(MAXLAY) ! secant of effective sun local path angle
       REAL DISTES ! distance of Earth from the sun
       REAL TAUZSN(MXCHAN) ! chan eff sun angle surface-to-space trans
       LOGICAL DOSUN ! do sun calc?
!
! for satellite viewing angle
       REAL SATANG ! input satellite scan angle (degrees)
       REAL SATZEN ! input satellite zenith angle (degrees)
       REAL SALT ! input satellite altitude (kilometers)
       REAL SVA ! satellite viewing angle (degrees)
!
! for RDRTP
       INTEGER IPROF ! profile loop counter
!       LOGICAL LWANT ! do you want this profile?
       REAL YOFF ! overall yoffset
       REAL YMODP1(18) ! yoffset of modules+1
       REAL DOPPPM ! Doppler shift in ppm
       REAL OPHASE ! orbit phase
!
! for SETYOF
       REAL YOFCHN(MXCHAN) ! channel total Yoffset
! used locally only
       INTEGER I ! loop counter
       INTEGER L ! loop counter
       INTEGER rtpclose ! for call to RTP close interface routine
       REAL EVA ! (Earth) local viewing angle
       REAL CONV ! degrees to radians conversion factor
       REAL ANGMAX ! maximum allowed viewing angle
       REAL FRQDEF ! Default freq
       REAL RJUNK1 ! junk/work
       REAL RJUNK2 ! another junk/work
       REAL CO2PPM ! Profile mean dry air CO2 mixing ratio
       REAL PLAY(MAXLAY) ! layer mean pressure
       REAL YOFPRE ! previous YOFF
!
! Profile data structure
       INTEGER ISTAT
! RECORD /RTPPROF/ PROF ! profile
! RECORD /RTPHEAD/ HEAD ! header data
! RECORD /RTPATTR/ HATT(MAXNATTR) ! header attributes
! RECORD /RTPATTR/ PATT(MAXNATTR) ! profile attributes
       type ( RTPPROF ) :: PROF ! profile

       type ( RTPHEAD ) :: HEAD ! header data
       type ( RTPPROF ) :: PATT(MAXNATTR) ! profile attributes


! G-SDSU plugging
 integer :: i_dom, j_dom   ! domain loop indices

!
! Boundary pressure levels
!-----------------------------------------------------------------------
! SAVE STATEMENTS
!-----------------------------------------------------------------------
! none
!
!***********************************************************************
!***********************************************************************
! EXECUTABLE CODE
!***********************************************************************
!***********************************************************************

 select case(trim(proc))
 case('init')

#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG SARTA_DRIVER: Initialize SARTA module.',myrank, masterproc

!
! CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02
       YOFPRE=-999.0 ! dummy
!
! --------------------------
! Assign the I/O unit number
! --------------------------
       IOUN=11
!
! ---------
! Calc PLAY
! ---------

! Mean layer pressure (KLAYERS definition) (Don't need this for G-SDSU)
       DO L=1,MAXLAY
          PLAY(L) = ( PLEV(L+1) - PLEV(L) )/LOG( PLEV(L+1)/PLEV(L) )
       ENDDO

! -----------------------------
! Read in the reference profile 
!
! This routine reads the reference profile data from a simple text file.
! After skipping over the header, the routine reads in the profile variables
! (altitude, presure, temperature, and CO2, H2O, O3, CO, and CH4 amounts).
! These are read in the order they occur in the file (that is, 1=lowest
! altitude to 100=highest altitude), but they are reversed in the arrays
! (ie element 1=highest altitude and 100=lowest altitude).
! Reference parameters are used for CALPAR and SUNPAR.
! -----------------------------

       CALL RDPROF(IOUN, FNPREF, RPNAM, RALT, RDZ, RPRES, RTEMP, &
                   RFAMNT, RWAMNT, ROAMNT, RCAMNT, RMAMNT, RSAMNT,RHAMNT,RNAMNT)


!
!  initialize emissivity table (tunable)
!
    call sarta_init_emissivity

!
! ---------------------
! Get run-time command-line info: 'FIN', 'FOUT', 'LRHOT', 'LISTP'  
! (Don't need this routine for G-SDSU. Toshi Masui Oct 2013)
!
! to read get command-line processing info, e.g., input / output name
! ---------------------
!       CALL RDINFO(FIN, FOUT, LRHOT, NWANTP, LISTP)


!
! ---------------------------
! Open & check input RTP file 
!
! This routine opens the input RTP file and performs some error checking
! to make sure the file contains the required info.
! ---------------------------
 !print*,'call OPNRTP'

       CALL OPNRTP(PTYPE, NCHAN, FCHAN, LSTCHN, INDCHN, &
          IH2O, IO3, ICO, ICH4, ICO2, ISO2, IHNO3, IN2O, &
          HEAD, PATT, LCO2PM)
!print*,'ICO2',ICO2

!
! ------------------------
! Read the coef data files
!
! The seven main binary data files containing the AIRS fast transmittance
! coefficients are opened and read one channel at a time.  The channel
! list LSTCHN determines what channels to read/use.  The seven sets of
! coefs are each stored in their own arrays.  Next, the peturbation
! coefficients for the trace gases CO2,N2O,SO2,& HNO3 are read in from
! their corresponding data files.  Next, the OPTRAN water coefficients
! (and OPTRAN water grid and raw predictor average values) are read in
! from the OPTRAN data file.  Next, the downwelling thermal "F factor"
! coefficients read in from their file.  Next, the "fixed gases"
! adjustment factor "fx" is read in from its file.  Finally, the non-LTE
! coefficients are read in from their file.
! ------------------------
! Read the "A" coef database
 !print*,'call RDCOEF A'

       AORB='A'
       CALL RDCOEF ( AORB, IOUN, NCHAN, INDCHN, SETCHN, &
           NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7, &
          CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
          COEF1A, COEF2A, COEF3A, COEF4A, COEF5A, COEF6A, COEF7A, &
          NCHCO2, INDCO2, CLICO2, COFCOA, &
          NCHN2O, INDN2O, CLIN2O, COFN2A, &
          NCHSO2, INDSO2, CLISO2, COFSOA, &
          NCHHNO, INDHNO, CLIHNO, COFHNA, &
          NCHH2O, INDH2O, CLIH2O, WAZOP, WAVGOP, COFH2A, &
           FREQA, COEFFA, NCHNTE, CLISTN, COEFNA, FX)


!
! Get and apply multipler tuning to "A" coefficients
 !print*,'call TUNMLT'
       CALL TUNMLT(AORB, IOUN, NCHAN, INDCHN, SETCHN, &
           NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7, &
          CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
          COEF1A, COEF2A, COEF3A, COEF4A, COEF5A, COEF6A, COEF7A, &
          NCHCO2, INDCO2, CLICO2, COFCOA, &
          NCHN2O, INDN2O, CLIN2O, COFN2A, &
          NCHSO2, INDSO2, CLISO2, COFSOA, &
          NCHHNO, INDHNO, CLIHNO, COFHNA, &
          NCHH2O, INDH2O, CLIH2O, WAZOP, WAVGOP, COFH2A, &
           FREQA, COEFFA, NCHNTE, CLISTN, COEFNA, FX)
! Read the "B" coef database
 !print*,'call RDCOEF B'

       AORB='B'
       CALL RDCOEF (AORB, IOUN, NCHAN, INDCHN, SETCHN, &
           NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7, &
          CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
          COEF1B, COEF2B, COEF3B, COEF4B, COEF5B, COEF6B, COEF7B, &
          NCHCO2, INDCO2, CLICO2, COFCOB, &
          NCHN2O, INDN2O, CLIN2O, COFN2B, &
          NCHSO2, INDSO2, CLISO2, COFSOB, &
          NCHHNO, INDHNO, CLIHNO, COFHNB, &
          NCHH2O, INDH2O, CLIH2O, WAZOP, WAVGOP, COFH2B, &
           FREQB, COEFFB, NCHNTE, CLISTN, COEFNB, FX)
 !print*,'call TUNMLT'

!
! Get and apply multipler tuning to "B" coefficients
       CALL TUNMLT( AORB, IOUN, NCHAN, INDCHN, SETCHN, &
           NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7, &
          CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
          COEF1B, COEF2B, COEF3B, COEF4B, COEF5B, COEF6B, COEF7B, &
          NCHCO2, INDCO2, CLICO2, COFCOB, &
          NCHN2O, INDN2O, CLIN2O, COFN2B, &
          NCHSO2, INDSO2, CLISO2, COFSOB, &
          NCHHNO, INDHNO, CLIHNO, COFHNB, &
          NCHH2O, INDH2O, CLIH2O, WAZOP, WAVGOP, COFH2B, &
           FREQB, COEFFB, NCHNTE, CLISTN, COEFNB, FX)
!
! Calc OPTRAN absorption coefficient scaling factor WAOP
       WAOP(1)=WAZOP(1)
       DO L=2,MXOWLY
          WAOP(L)=WAZOP(L) - WAZOP(L-1)
       ENDDO
!
! -------------------------
! Read the A/B weight table
! -------------------------

 !print*,'call RDDYAB'
       CALL RDDYAB(IOUN, ETIME, FREQAB, WIDTH, INDMOD, IABWGT, &
          DELYAB)
!
! --------------------------
! Read in the solar radiance
! --------------------------
 !print*, 'call RDSUN'

       AORB='A'
       CALL RDSUN(AORB, IOUN, INDCHN, HSUNA)
       AORB='B'
       CALL RDSUN(AORB, IOUN, INDCHN, HSUNB)
!
       DISTES=1.496E+11 ! distance Earth to Sun 
!
! --------------------
! Check FREQ and FCHAN
! --------------------
! Note: FREQ comes the coef data, while FCHAN comes
! from the input RTP file read by OPNRTP. It is possible
! that FCHAN is "nodata", so we check the first element.
!
       IF (FCHAN(1) .GT. 640 .AND. FCHAN(1) .LT. 2670) THEN
          DO I=1,NCHAN
! Check old channel freqs before updating
             FRQDEF=FREQAB(LSTCHN(I))
             RJUNK1=ABS(FRQDEF - FCHAN(I))
             RJUNK2=0.05*FRQDEF/1200.0 ! ~5% of a channel fullwidth
             IF (RJUNK1 .GT. RJUNK2) THEN
                WRITE(IOINFO,1010) I, LSTCHN(I), FRQDEF, FCHAN(I)
 1010 FORMAT('Warning &! index=',I4,', chan ID=',I4, &
                ', fastmodel default freq=',F8.3,', RTP freq=',F8.3)
             ENDIF
             HEAD%vchan(I)=FRQDEF
          ENDDO
       ELSE
! Assign freqs without checking old freqs
          DO I=1,NCHAN
             HEAD%vchan(I)=FREQAB(LSTCHN(I))
          ENDDO
       ENDIF



!
! ------------------------
! Open the output RTP file
! ------------------------
       MODE='c'

! ISTAT=rtpopen(FOUT, MODE, HEAD, HATT, PATT, IOPCO)
!cc
! print *, 'read open status = ', ISTAT
!cc
!
! -----------------------------------------------
! All channels from sets 1, 2, and 3 are to use a
! fake effective sun angle layer-to-space trans
! -----------------------------------------------
       NFAKE=0
!
       DO I=1,NCHN1
          NFAKE=NFAKE + 1
          INDFAK(NFAKE)=INDCHN( CLIST1(I) )
       ENDDO
!
       DO I=1,NCHN2
          NFAKE=NFAKE + 1
          INDFAK(NFAKE)=INDCHN( CLIST2(I) )
       ENDDO
!
       DO I=1,NCHN3
          NFAKE=NFAKE + 1
          INDFAK(NFAKE)=INDCHN( CLIST3(I) )
       ENDDO
!

!
! ###########################  CALCLUATION  #############################################
!

 case('calc')

#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG SARTA_DRIVER: Calculate radiance.'

!-------------------------------------------------------------
!                  Beggining of domain loop
!-------------------------------------------------------------

 DOMAIN_LOOP: do i_dom = myi_start , myi_end ; do j_dom = myj_start , myj_end 

!
! --------------
! Read input RTP (--> input is now G-SDSU atmosphere-land parameters)
! --------------
 !print*, 'call RDRTP',ICO2


       CALL RDRTP( atmos(i_dom,j_dom,:), atmos_stag(i_dom,j_dom,:), surface(i_dom,j_dom), &
                   emiss_sarta(i_dom,j_dom,:), HEAD, PATT(1),  &
          IH2O, IO3, ICO, ICH4, ICO2, ISO2, IHNO3, IN2O, PTYPE, &
          RALT, LCO2PM, NLAY, NEMIS, LAT, LON, SATANG, SATZEN, &
          SALT, SUNANG, PSURF, TSURF, CO2PPM, FEMIS, XEMIS, XRHO, &
          TEMP, WAMNT, OAMNT, CAMNT, MAMNT, FAMNT, SAMNT, HAMNT, NAMNT, &
           ALT, PROF, ISTAT, YOFF, YMODP1, DOPPPM, OPHASE )



!

!(Don't need this routine for G-SDSU. Toshi Masui Oct 2013)
!       IF (ISTAT .EQ. -1) GOTO 9999 ! reached End Of File  
!       IF (.NOT. LWANT) THEN     
!          IPROF=IPROF+ 1
!          GOTO 10
!       ENDIF

!
! ------------------------
! Set the channel Yoffsets (toshii: this may go to initialization)
! ------------------------
 !print*, 'call SETYOF'
       CALL SETYOF(NCHAN, LSTCHN, PROF%rtime, YOFF, YMODP1, &
          DOPPPM, ETIME, FREQAB, WIDTH, INDMOD, IABWGT, &
          DELYAB, YOFCHN)

! -------------------
! Do freq adjustments (toshii: this may go to initialization)
! -------------------
       IF (YOFF .NE. YOFPRE) THEN
!print*,'call SUMCOF'
          CALL SUMCOF( NCHAN, NCHNTE, NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, &
          NCHN6, NCHN7, NCHCO2, NCHN2O, NCHSO2, NCHHNO, NCHH2O, &
         LSTCHN, CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
         CLISTN, CLICO2, CLIN2O, CLISO2, CLIHNO, CLIH2O, YOFCHN, &
          FREQA, HSUNA, COEFFA, COEFNA, COEF1A, COEF2A, COEF3A, COEF4A, &
         COEF5A, COEF6A, COEF7A, COFCOA, COFN2A, COFSOA, COFHNA, COFH2A, &
          FREQB, HSUNB, COEFFB, COEFNB, COEF1B, COEF2B, COEF3B, COEF4B, &
         COEF5B, COEF6B, COEF7B, COFCOB, COFN2B, COFSOB, COFHNB, COFH2B, &
           FREQ, HSUN, COEFF, COEFN, COEF1, COEF2, COEF3, COEF4, &
          COEF5, COEF6, COEF7, COFCO2, COFN2O, COFSO2, COFHNO, COFH2O)
!
! Update YOFPRE for the newly summed coef data
          YOFPRE=YOFF
       ENDIF


!
! -------------------------------------
! Determine bottom layer, CO2, & angles
! -------------------------------------
!print*,'call GETBOT'
       CALL GETBOT(NLAY, PLEV, PSURF, LBOT, BLMULT) 

!LBOT is used 
!

! Calc the fractional bottom layer air temperature
!cc
! TEMP(LBOT)=TEMP(LBOT-1) + BLMULT*( TEMP(LBOT) - TEMP(LBOT-1) )
! Above line commented out & replaced by Scott Hannon, 24 July 2003.
! Mistakenly treats T at the center of the layer above as T at the
! bottom of the layer above.
!cc

! CO2 profile switch
       IF (ICO2 .LT. 1) THEN
          LCO2=.FALSE.
       ELSE
          LCO2=.TRUE.
       ENDIF
! N2O profile switch
       IF (IN2O .LT. 1) THEN
          LN2O=.FALSE.
       ELSE
          LN2O=.TRUE.
       ENDIF
! SO2 profile switch
       IF (ISO2 .LT. 1) THEN
          LSO2=.FALSE.
       ELSE
          LSO2=.TRUE.
       ENDIF
! HNO3 profile switch
       IF (IHNO3 .LT. 1) THEN
          LHNO3=.FALSE.
       ELSE
          LHNO3=.TRUE.
       ENDIF


!
!  fractional bottom layer mean temperature 
!  (Don't need this for G-SDSU; Toshi Matsui Oct 2013)
!
!       IF (PTYPE .EQ. AIRSLAY) THEN
!         ! Copy pseudo level temperatures to another array
!          DO I=1,LBOT
!             TPSEUD(I)=TEMP(I)
!          ENDDO
          ! Convert temperatures
!          CALL MEAN_T(LBOT, PLEV, PSURF, TPSEUD, TEMP)
!       ELSE
         ! Calc mean pressure for bottom fractional layer
!          RJUNK1 = ( PSURF - PLEV(LBOT) )/LOG( PSURF/PLEV(LBOT) )
         ! Do interpolation for fractional bottom layer mean temperature
          ! assuming T is in linear in log(P)
!          RJUNK2=( TEMP(LBOT) - TEMP(LBOT-1) )/ &
!             LOG( PLAY(LBOT)/PLAY(LBOT-1) ) ! slope
!          TEMP(LBOT)=RJUNK2*LOG( RJUNK1/PLAY(LBOT-1) ) + TEMP(LBOT - 1)  
!       ENDIF


!
! Check satellite elevation
!
       IF (SALT .GT. 0) THEN
          ! Warn and use default if invalid
          IF (SALT .LT. XSALT-150 .OR. SALT .GT. XSALT+150) THEN
             WRITE(IOINFO,1020) IPROF, SALT, XSALT
 1020 FORMAT('Warning &! Profile',I5, &
                ': replacing invalid input satellite altitude ', &
                1PE11.4,' with default ',1PE11.4,' km')
             SALT=XSALT
          ENDIF
       ELSE
          SALT=XSALT
       ENDIF

!
! Convert SATZEN or SATANG to viewing angle
!
       IF (SATZEN .GE. 0 .AND. SATZEN .LT. 63) THEN
          ! Convert zenith angle at surface to view angle at satellite
          SVA=SACONV( SATZEN, SALT*1000 )/CONV
       ELSE
          ! Check if scan angle is valid
          IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
             ! View angle should be within a few degrees of scan angle
             SVA=ABS( SATANG )
          ELSE
             WRITE(IOERR,1030) IPROF, SATZEN, SATANG
 1030 FORMAT('Error &! Profile',I5, &
                ': invalid angles for SATZEN ',1PE11.4, &
                ' and SATANG ',E11.4)
             STOP
          ENDIF
       ENDIF

       ANGMAX=53.E0 ! max satellite view angle (49.5 scan + 3.5 spacecraft)
       IF (SVA .GT. ANGMAX) THEN
          ! Truncate angle if too big
          WRITE(IOINFO,1040) IPROF, SVA
 1040 FORMAT('Warning &! Profile',I5,': truncating view angle ', &
             1PE11.4,' to 53 degrees')
          SVA=ANGMAX
       ENDIF

       ! Convert from satellite to earth viewing angle (in radians)

       DO L=1,LBOT
          EVA=VACONV(LAT, SVA, SALT, ALT(L))
          SECANG(L)=1.0E+0/COS(EVA)
       ENDDO

!
! Calc total sun angle secant 
!
       DOSUN=.FALSE.
       IF (SUNANG .GE. 0.0 .AND. SUNANG .LT. 89.9) DOSUN=.TRUE.
       IF (DOSUN) THEN
          SUNCOS=COS(CONV*SUNANG)     !cosine of solar zenith angle
          SZALAY=SACONV(SUNANG,ALT(1))
          SCOS1=COS(SZALAY)
          RJUNK2=SECANG(LBOT) + 1.0/SUNCOS ! Total secant

          ! Calc non-unity fudge factor if total secant > 9
          IF (RJUNK2 .GT. 9.0) THEN
          ! fudge factor = true_total_secant/calc_total_secant
             SUNFDG=RJUNK2/9.0
          ! truncated solar angle to use to calc SECSUN
             RJUNK1=ACOS( 1.0/(9.0 - SECANG(LBOT)) )/CONV
          ELSE
             SUNFDG=1.0
             RJUNK1=SUNANG
          ENDIF

          DO L=1,LBOT
             SZALAY=SACONV(RJUNK1,ALT(L))
             SECSUN(L)=SECANG(L) + 1.0E+0/COS(SZALAY)
          ENDDO
       ENDIF


!
! OK Final set up of profile is done. From here  conversion into 
! optical properties
!

 !print*,'Compute IR transmittances '
!
! -----------------------------------
! Calculate the fast trans predictors
! -----------------------------------
!
       CALL CALPAR (LBOT, &
          RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RSAMNT,RHAMNT,RNAMNT, &
           TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT, &
          RPRES,SECANG, LAT, FX, RDZ, &
           LCO2, LN2O, LSO2, LHNO3,LCO2PM,CO2PPM,CO2TOP,FIXMUL,CONPRD, &
         FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7, &
         WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7, &
         OPRED1,OPRED2, OPRED4,OPRED5,OPRED6,OPRED7, &
         MPRED3,CPRED4,TRCPRD,CO2MLT,SO2MLT,HNOMLT,N2OMLT )
!
! -----------------------------------
! Calculate the OPTRAN H2O predictors
! -----------------------------------
       CALL CALOWP ( LBOT, WAMNT, RPRES, TEMP, SECANG, WAZOP, WAVGOP, &
          WAANG, LOPMIN, LOPMAX, LOPUSE, H2OPRD, LOPLOW, DAOP )
!
! ----------------------------------
! Calculate the layer transmittances
! ----------------------------------
! Calculate TAU for set 1 thru 7
!
       CALL CALT1( INDCHN, LBOT, BLMULT, NCHN1, CLIST1, COEF1, &
           FIXMUL, CONPRD, FPRED1, WPRED1, OPRED1, TRCPRD, &
           INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
           INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, &
           INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX, LOPLOW, &
           LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ)
!
       CALL CALT2( INDCHN, LBOT, BLMULT, NCHN2, CLIST2, COEF2, &
          FIXMUL, CONPRD, FPRED2, OPRED2, WPRED2, TRCPRD, &
          INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
          INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, TAU, TAUZ)
!
       CALL CALT3( INDCHN, LBOT, BLMULT, NCHN3, CLIST3, COEF3, &
           FIXMUL, CONPRD, FPRED3, MPRED3, WPRED3, TRCPRD, &
           INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT, &
           INDN2O, COFN2O, N2OMLT, &
           INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE, &
             WAOP, DAOP, WAANG, TAU, TAUZ)
!
       LTAU=.TRUE.
       XZ=1.0
!
       CALL CALT4( LTAU, INDCHN, LBOT, BLMULT, NCHN4, CLIST4, &
           COEF4, FIXMUL, CONPRD, FPRED4, CPRED4, OPRED4, WPRED4, &
          TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
              XZ, TAU, TAUZ )
!
       CALL CALT5( LTAU, INDCHN, LBOT, BLMULT, NCHN5, CLIST5, &
           COEF5, FIXMUL, CONPRD, FPRED5, WPRED5, OPRED5, &
          TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
              XZ, TAU, TAUZ )
!
       CALL CALT6( LTAU, INDCHN, LBOT, BLMULT, NCHN6, CLIST6, &
           COEF6, FIXMUL, CONPRD, FPRED6, WPRED6, OPRED6, TRCPRD, &
          INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
          INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZ )
!
       CALL CALT7( LTAU, INDCHN, LBOT, BLMULT, NCHN7, CLIST7, &
           COEF7, FIXMUL, CONPRD, FPRED7, WPRED7, OPRED7, &
          TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
              XZ, TAU, TAUZ )

!
       DOSUN_IF: IF (DOSUN) THEN

!print*,'Compute NIR transmittance'

! ---------------------------------------------
! Calculate the fast trans predictors *for sun*
! ---------------------------------------------
!
          CALL SUNPAR ( LBOT, &
             RTEMP, RWAMNT, ROAMNT, RCAMNT, &
              TEMP, WAMNT, OAMNT, CAMNT, &
             RPRES, SECSUN, CONPRD, &
             FPRED4, FPRED5, FPRED6, FPRED7, &
             WPRED4, WPRED5, WPRED6, WPRED7, &
             OPRED4, OPRED5, OPRED6, OPRED7, &
             CPRED4, TRCPRD )
! --------------------------------------------
! Calculate the layer transmittances *for sun*
! --------------------------------------------
!
! Calc fake TAUZSN for sets 1, 2, and 3
          RJUNK1=SUNFDG*SECSUN(LBOT)
          CALL FAKETZ( NFAKE, INDFAK, TAUZ, SECANG(LBOT), &
             RJUNK1, TAUZSN)
!
! Calculate TAUZSN for sets 4 thru 7
!
          LTAU=.FALSE.
          XZ=SUNFDG
!
          CALL CALT4( LTAU, INDCHN, LBOT, BLMULT, NCHN4, CLIST4, &
             COEF4, FIXMUL, CONPRD, FPRED4, CPRED4, OPRED4, WPRED4, &
             TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
             XZ, TAU, TAUZSN )
!
          CALL CALT5( LTAU, INDCHN, LBOT, BLMULT, NCHN5, CLIST5, &
             COEF5, FIXMUL, CONPRD, FPRED5, WPRED5, OPRED5, &
             TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
             XZ, TAU, TAUZSN )
!
          CALL CALT6( LTAU, INDCHN, LBOT, BLMULT, NCHN6, CLIST6, &
             COEF6, FIXMUL, CONPRD, FPRED6, WPRED6, OPRED6, &
             TRCPRD, INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
             INDN2O, COFN2O, N2OMLT, XZ, TAU, TAUZSN )
!
          CALL CALT7( LTAU, INDCHN, LBOT, BLMULT, NCHN7, CLIST7, &
             COEF7, FIXMUL, CONPRD, FPRED7, WPRED7, OPRED7, &
             TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
             XZ, TAU, TAUZSN )
!
       ELSE
! No sun; set the sun surface-to-space trans to zero
          SUNCOS=0.0
          DO I=1,NCHAN
             TAUZSN(I)=0.0
          ENDDO
       ENDIF DOSUN_IF

!
! ---------------------------------------------------
! Set the emissivity & reflectivity for every channel
! ---------------------------------------------------
   !print*,'call SETEMS'
       CALL SETEMS( NCHAN, NEMIS, FREQ, FEMIS, XEMIS, &
          XRHO, LRHOT, EMIS, RHOSUN, RHOTHR)
!
! ----------------------
! Calculate the radiance
! ----------------------
   !print*,'call CALRAD'
       CALL CALRAD ( NCHAN, FREQ, TAU, TEMP, TSURF, EMIS, LBOT, &
          SUNCOS, RHOSUN, DISTES, HSUN, TAUZSN, &
          SECANG(LBOT), RHOTHR, LABOVE, COEFF, TAUZ, RAD, BT)

!
! -----------------
! Calculate non-LTE
! -----------------
       IF (DOSUN) THEN
          !print*,'call CALNTE'
          CALL CALNTE ( INDCHN, FREQ, TEMP, SUNCOS, SCOS1, SECANG(1), &
             NCHNTE, CLISTN, COEFN, CO2TOP, RAD , BT)
       ENDIF


!
! OUTPUT to G-SDSU 
!
 if( i_dom == myi_start .and. j_dom == myj_start ) then
     sarta_chan_id  (1:NCHAN) = LSTCHN(1:NCHAN)  ! AIRS channel ID [-]
     sarta_chan_freq(1:NCHAN) = FREQ(1:NCHAN)    ! AIRS channel frequenc [cm-1]

! For testing....
!     do I=1,NCHAN, 100
!        write(*,FMT='(I4,3x,F8.3,3x,F6.2)') LSTCHN(I), FREQ(I), BT(I)
!     enddo

!     do I=1,NCHAN
!        write(111,FMT='(I4,3x,F8.3,3x,F6.2)') LSTCHN(I), FREQ(I), BT(I)
!     enddo

 endif

 if    ( airs_unit == 'tb' ) then
    sarta_out (i_dom,j_dom,1:NCHAN) = BT (1:NCHAN)  ! AIRS brightness temperature [K]
 elseif(  airs_unit == 'rad' ) then
    sarta_out (i_dom,j_dom,1:NCHAN) = RAD(1:NCHAN)  ! AIRS radiance [Watts/m^2/cm^-1]
 endif




 enddo ; enddo DOMAIN_LOOP

!-------------------------------------------------------------
!                  End of domain loop
!-------------------------------------------------------------

 end select

#if MPI == 2
 if(masterproc) &
#endif
 print*,''


 END SUBROUTINE SARTA_DRIVER

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine sarta_init_emissivity
 implicit none
 integer,parameter :: io = 100
 character(len=200)   :: io_file
 character( len=200 ) :: dummy

!
! read emissivity data 
!
 io_file = trim(sarta_dir_data)//'Data_AIRS_apr08/sarta_emis_table.txt'

#if MPI == 2
 if(masterproc) &
#endif
 print*,'MSG sarta_init_emissivity: reading',trim(io_file)


 call open_oldfile(io,trim(io_file))
 read(io,*) dummy                              ! read string
 read(io,*) dummy,dummy, nemis_ap              ! number of emissivity anchor point 
 read(io,*) dummy,dummy, efreq_ap(1:nemis_ap)  ! frequency of emissivity anchor point
 read(io,*) dummy,dummy, emis_ap (1:nemis_ap)  ! emissivity at anchor point
 read(io,*) dummy,dummy, rho_ap  (1:nemis_ap)  ! albedo at anchor point
 close(io)

 return
 end subroutine sarta_init_emissivity

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
! not used any more.....
 subroutine read_sarta_param(HEAD,PATT)
 implicit none
 type( RTPHEAD ),intent(inout) :: HEAD ! sarta header data
 type( RTPPROF ),intent(inout) :: PATT(MAXNATTR) ! profile attributes
 integer,parameter :: io = 100
 character(len=200)   :: io_file
 character( len=200 ) :: dummy
 integer :: z


!
! read profile data
!
 io_file = trim(sarta_dir_data)//'Data_AIRS_apr08/sarta_p.txt'
 print*,'reading',trim(io_file)
 call open_oldfile(io,trim(io_file))

 read(io,*)  dummy            !read string
 read(io,*)  dummy,dummy, PATT(1)%plat  !used for gases gravity correction
 read(io,*)  dummy,dummy, PATT(1)%plon  !not used
 read(io,*)  dummy,dummy, PATT(1)%ptime !not used
 read(io,*)  dummy,dummy, PATT(1)%stemp !important
 read(io,*)  dummy,dummy, PATT(1)%salti !not used currently...
 read(io,*)  dummy,dummy, PATT(1)%spres !important
 read(io,*)  dummy,dummy, PATT(1)%nemis !important (done)
 read(io,*)  dummy,dummy, PATT(1)%efreq(1:PATT(1)%nemis)  !important (done)
 read(io,*)  dummy,dummy, PATT(1)%emis (1:PATT(1)%nemis)  !important (done)
 read(io,*)  dummy,dummy, PATT(1)%rho  (1:PATT(1)%nemis)  !important (done)
 read(io,*)  dummy,dummy, PATT(1)%nlevs
 read(io,*)  dummy,dummy, PATT(1)%plevs(1:PATT(1)%nlevs)
 read(io,*)  dummy,dummy, PATT(1)%ptemp(1:PATT(1)%nlevs)
 read(io,*)  dummy,dummy, PATT(1)%gamnt(1:PATT(1)%nlevs,1)
 read(io,*)  dummy,dummy, PATT(1)%gamnt(1:PATT(1)%nlevs,2)
 read(io,*)  dummy,dummy, PATT(1)%gamnt(1:PATT(1)%nlevs,3)
 read(io,*)  dummy,dummy, PATT(1)%gamnt(1:PATT(1)%nlevs,4)
 read(io,*)  dummy,dummy, PATT(1)%gamnt(1:PATT(1)%nlevs,5)
 read(io,*)  dummy,dummy, PATT(1)%gamnt(1:PATT(1)%nlevs,6)
 read(io,*)  dummy,dummy, PATT(1)%gamnt(1:PATT(1)%nlevs,7)
 read(io,*)  dummy,dummy, PATT(1)%gamnt(1:PATT(1)%nlevs,8)
 read(io,*)  dummy,dummy, PATT(1)%gxover(1:ngas)  !not used
 read(io,*)  dummy,dummy, PATT(1)%scanang !scan angle
 read(io,*)  dummy,dummy, PATT(1)%satzen !satellite zenigt angle -9999
 read(io,*)  dummy,dummy, PATT(1)%zobs  !satellight height Convert SATZEN or SATANG to viewing angle

 print*,'done reading'
 close(io)

!
! derive altitude from hydrostatic equation.
!Z2 - Z1 = RT/g * ln(p1/p2) ..... Eq(b)
!replace p1 and p2 by 1000 hPa and 500 hPa respectively
!therefore ln(2) = 0.69
!R = universal gas constant for dry air = 2.87 * 10^2 kg^1 K^1
!g = 9.81 ms^2
!T = mean temperature through layer (K)
!Z1, Z2 = heights of isobaric surfaces p1,p2 (in metres)
!by substitution, and allowing for the fact that T in the original equation is Kelvin, we have, for the 500/1000 hPa layer....
!mean temperature (degC) = (Thickness/20.3) - 273 ..... Eq(c)
!
 PATT(1)%palts(PATT(1)%nlevs) = 0.

 do z = PATT(1)%nlevs-1, 1, -1

   PATT(1)%palts(z)  = 2.87e2*(0.5*PATT(1)%ptemp(z+1)+0.5*PATT(1)%ptemp(z))/9.81&
                       *log(PATT(1)%plevs(z+1)/PATT(1)%plevs(z)) + PATT(1)%palts(z+1)
!   print*,z,PATT(1)%palts(z)
 enddo


 return
 end subroutine read_sarta_param

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 end module module_sarta

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU
