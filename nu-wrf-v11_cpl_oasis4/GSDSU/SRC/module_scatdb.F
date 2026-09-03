!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_scatdb
 use module_floating_point
 use module_const
 implicit none

!--------------------------------------------------------------------------------------------------
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
!  This module is designed to wire G. Liu's non-sherical single-scattering database (SCATDB) and the 
!  microwave optics module in the G-SDSU. Essentially, the G-SDSU treat particle as sphere, and 
!  SCATDB's LUT must be slightly arranged for interpolation. 
!  One caution is that SCATDB's particle size is quite limitted. So, if particle, size goes beyond the 
!  available size ranges, it return logical output success = .false. 
!
!  For more details, please the note in scatdb.c
!
! History:
! 01/2011  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
! Liu, Guosheng, 2008: A Database of Microwave Single-Scattering Properties for Nonspherical Ice Particles.
!    Bull. Amer. Meteor. Soc., 89, 1563-1570. 
!
!-----------------------------------------------------------------------------------------------------

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
  private   ! encapsulate all variables and subourtines (non accessible)

!
! public (accessible) subroutines
!
  public :: scatdb_prep, & ! create non-spherical microwave single-scattering properties LUTs
            scatdb_interp  ! interpolate single-scattering propertes for a givne sphere radius

!
!
! SCATDB V2
!  anchor points frequencies are: 
!  3, 5, 9, 10, 13.4, 15, 19, 24.1, 35.6, 50, 60, 70, 80, 85.5, 90, 94, 118, 150, 166, 183, 220, 340 GHz. 
!
!  SHAPE-ID            SHAPE#          SHAPE                         
!   hexl                 0             long hexagonal column l/d=4   
!   hexs                 1             short hexagonal col   l/d=2   
!   hexb                 2             block hex col         l/d=1   
!   hexf                 3             thick hex plate       l/d=0.2 
!   hexp                 4             thin  hex plate      l/d=0.05 
!   ros3                 5             3-bullet rosette              
!   ros4                 6             4-bullet rosette              
!   ros5                 7             5-bullet rosette              
!   ros6                 8             6-bullet rosette              
!   sstr                 9             sector-like snowflake         
!   sden                10             dendrite snowflake            
!
!
 integer :: i,d


 type lut_parameter
   real :: d_mm    ! sphere diameter [mm]
   real :: dmax    ! maximum dimension [mm]
   real :: Cext    ! total extincion cross section [m2]
   real :: Csca    ! scattering cross section [m2]
   real :: Cbsca   ! backscattering cross section [m2]
   real :: asym    ! asymetry parameter [-]
 end type lut_parameter
 type ( lut_parameter ) :: lut(20)


 integer,parameter :: nmax(0:10) = (/7, 7, 7, 7, 7,20,20,20,20,20,20/) 

! Anchor points of maximum dimension (μm) of ice particles used for DDA computations.
 real :: dmax_bin(0:10,1:21)
 data ((dmax_bin(i,d),i= 0,0),d=1,20) / &  ! Long columns
   121.e0,  242.e0,  484.e0,  967.e0, 2417.e0, 3626.e0, 4835.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0/
 data ((dmax_bin(i,d),i= 1,1),d=1,20) / &  ! Short columns
    83.e0,  165.e0,  330.e0,  661.e0, 1652.e0, 2478.e0, 3304.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0/
 data ((dmax_bin(i,d),i= 2,2),d=1,20) / &  ! Block columns
    66.e0,  132.e0,  263.e0,  527.e0, 1316.e0, 1974.e0, 2632.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0/
 data ((dmax_bin(i,d),i= 3,3),d=1,20) / &  ! Thick plates
    81.e0,  162.e0,  325.e0,  649.e0, 1623.e0, 2434.e0, 3246.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0/
 data ((dmax_bin(i,d),i= 4,4),d=1,20) / &  ! Thin plates
   127.e0,  253.e0,  506.e0, 1012.e0, 2529.e0, 3794.e0, 5059.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0,&
  -999.e0, -999.e0, -999.e0, -999.e0, -999.e0, -999.e0/
 data ((dmax_bin(i,d),i= 5,5),d=1,20) / &  ! 3-bullet rosette 
    50.e0,  100.e0,  200.e0,  300.e0,  400.e0,  500.e0,  750.e0,&
  1000.e0, 1500.e0, 2000.e0, 2500.e0, 3000.e0, 3500.e0, 4000.e0,&
  5000.e0, 6000.e0, 7000.e0, 8000.e0, 9000.e0,10000.e0/
 data ((dmax_bin(i,d),i= 6,6),d=1,20) / &  ! 4-bullet rosette 
    50.e0,  100.e0,  200.e0,  300.e0,  400.e0,  500.e0,  750.e0,& 
  1000.e0, 1500.e0, 2000.e0, 2500.e0, 3000.e0, 3500.e0, 4000.e0,& 
  5000.e0, 6000.e0, 7000.e0, 8000.e0, 9000.e0,10000.e0/
 data ((dmax_bin(i,d),i= 7,7),d=1,20) / &  ! 5-bullet rosette 
    50.e0,  100.e0,  200.e0,  300.e0,  400.e0,  500.e0,  750.e0,& 
  1000.e0, 1500.e0, 2000.e0, 2500.e0, 3000.e0, 3500.e0, 4000.e0,& 
  5000.e0, 6000.e0, 7000.e0, 8000.e0, 9000.e0,10000.e0/
 data ((dmax_bin(i,d),i= 8,8),d=1,20) / &  ! 6-bullet rosette 
    50.e0,  100.e0,  200.e0,  300.e0,  400.e0,  500.e0,  750.e0,& 
  1000.e0, 1500.e0, 2000.e0, 2500.e0, 3000.e0, 3500.e0, 4000.e0,& 
  5000.e0, 6000.e0, 7000.e0, 8000.e0, 9000.e0,10000.e0/
 data ((dmax_bin(i,d),i= 9,9),d=1,20) / &  ! Sector snowflake
    50.e0,  100.e0,  200.e0,  300.e0,  400.e0,  500.e0,  750.e0,& 
  1000.e0, 1500.e0, 2000.e0, 2500.e0, 3000.e0, 3500.e0, 4000.e0,& 
  5000.e0, 6000.e0, 7000.e0, 8000.e0, 9000.e0,10000.e0/
 data ((dmax_bin(i,d),i=10,10),d=1,20) / &  ! Dendrite snowflake
    75.e0,  101.e0,  201.e0,  277.e0,  369.e0,  526.e0,  621.e0,&
   759.e0,  854.e0,  976.e0, 1073.e0, 1430.e0, 2076.e0, 3113.e0,&
  3558.e0, 4151.e0, 4981.e0, 6227.e0, 8302.e0,12454.e0/

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine scatdb_prep(freq,temp,type_shape )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine open SCATDB database, and put in LUT as a function of frequency, temperature and 
! crystal type.
!
! History:
! 01/2011  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 real,intent(in) :: freq             ! freqeuncy [GHz]
 real,intent(in) :: temp             ! air temperature [K]
 integer,intent(in) :: type_shape    ! 1~10 shape type

 integer :: n !loop index
 integer :: is_loaded = 0
 integer ::  nshp,iret
 real :: cabss,csca,cbsc,g,re,ph(37)
 real :: t,f,dmax
 real,parameter :: pi = 3.141592653589793

!
! intialize
!
 lut(:)%d_mm = -999.e0   ! sphere diameter (ice density = 0.917g/cm3) [mm]
 lut(:)%dmax = -999.e0   ! maximum dimension [mm]
 lut(:)%Cext = -999.e0   ! total extincion cross section [m2]
 lut(:)%Csca = -999.e0   ! scattering cross section [m2]
 lut(:)%Cbsca= -999.e0   ! backscattering cross section [m2]
 lut(:)%asym = -999.e0   ! asymetry params [-]  
 
!
! filter input first
!
 if( freq < 3. .or. freq > 340. ) then
    print*, 'MSG scatdb_prep: freq=',freq,'It must be within 3~340 (GHz).'
    stop
 endif
 f    = freq  ! frequency [GHz]

 t    = temp  ![K]

 if( type_shape < 0 .or. type_shape > 10 ) then
    print*, 'MSG scatdb_prep: type_shape=',type_shape,'It must be within 0~10.'
    stop
 endif
 nshp = type_shape


  do n = 1, nmax(nshp) !dmax bin loop

     dmax = dmax_bin(nshp,n) !maximum dimesion bin [micron]

     !
     ! derive SCATDB parameters
     !
     call scatdb(f,t,nshp,dmax,cabss,csca,cbsc,g,ph,re,iret,is_loaded)


     !
     ! Conver following parameters and units compatible to G-SDSU.
     !  cabs - absorption cross section for particle m^2  real/float                   
     !  csca - scattering cross section for particle m^2  real/float                   
     !  cbsc - backscattering cross section for particle m^2 real/float         
     !  g - asymmetry parameter for particle            real/float  
     !

     lut(n)%d_mm = 2.e0*re*1.e-3  ! solid-sphere diameter [mm]
     lut(n)%dmax = dmax*1.e-3     ! maximum dimension [mm]
     lut(n)%Cext = (cabss+csca)   ! total extincion cross section [m2]
     lut(n)%Csca = (      csca)   ! scattering cross section [m2]
     lut(n)%Cbsca= (      cbsc)/(4.e0*3.14159265)   ! backscattering cross section [m2] / 4pi radiance
     lut(n)%asym = g              ! asymetry parameter [-]


!
! Extinciton efficiency = Extinciton cross section / (pai * r*r )
!
  enddo


 end subroutine scatdb_prep

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine scatdb_interp(type_shape, d_mm, cross, qext, qsca, asym, qbsca, success)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine interplote/extrapolate prepared LUTs from SCATDB database.
!
! History:
! 03/2013  Toshi Matsui@NASA GSFC : corrected cross section and effieincy term.
! 01/2011  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: type_shape    ! 1~10 shape type
 real,intent(in)  :: d_mm    ! sphere diameter [mm]
 real,intent(in)  :: cross    ! cross section of fully particle [m2] (= pai*r*r)
 real,intent(out) :: qext    ! total extincion efficiency [-]
 real,intent(out) :: qsca    ! scattering efficiency [-]
 real,intent(out) :: asym    ! asymetry parameter [-]
 real,intent(out) :: qbsca   ! backscattering efficiency [-] / 4pi radiance
 logical,intent(out) :: success ! status of interpolation
 integer :: nshp  !type 
 integer :: n !loop inex
 real :: wgt1, wgt2  !weight for interpolation
 type ( lut_parameter ) :: dlut

!
! filtering input
!
 if( type_shape < 0 .or. type_shape > 10 ) then
    print*, 'MSG scatdb_prep: type_shape=',type_shape,'It must be within 0~10.'
    stop
 endif
 nshp = type_shape


!
! initialize output
!
 qext = -999. ; qsca = -999. ; qbsca= -999. ; asym = -999. ; success = .false.

 if( d_mm < lut(1)%d_mm ) then !use Mie routine
   return
 elseif( d_mm >= lut( nmax(nshp) )%d_mm  ) then 
   return 
 else  !interpolate

  !
  ! interpolate medium particle
  !
   size_loop: do n = 1, nmax(nshp)-1 !dmax bin loop

      if( d_mm >= lut(n)%d_mm .and. d_mm < lut(n+1)%d_mm ) then

          wgt2 = (d_mm - lut(n)%d_mm)/( lut(n+1)%d_mm - lut(n)%d_mm )  ![-]
          wgt1 = 1.- wgt2       
 
          !
          ! weiting interplation and convert cross section to efficiency
          !
          qext  =  (wgt1 * lut(n)%Cext  + wgt2 * lut(n+1)%Cext )/cross ![-]
          qsca  =  (wgt1 * lut(n)%Csca  + wgt2 * lut(n+1)%Csca )/cross ![-]
          qbsca =  (wgt1 * lut(n)%Cbsca + wgt2 * lut(n+1)%Cbsca)/cross * 4.*const_pi ![-]
          asym  =  (wgt1 * lut(n)%asym  + wgt2 * lut(n+1)%asym )      ![-]

          !interpolation results
          success = .true.

          exit size_loop

      endif

   enddo  size_loop

 endif


 end subroutine scatdb_interp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_scatdb

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
