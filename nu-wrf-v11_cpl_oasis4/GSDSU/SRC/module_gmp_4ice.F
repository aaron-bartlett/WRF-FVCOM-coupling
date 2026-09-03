! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 module module_gmp_4ice
 use module_const
 implicit none


 private  !everything is unless declare it as public. 


!
! public subroutines
!
 public  ::               &   
    gmp13_size_rho_diag , &  ! diagnostic size and density mapping for GMP13 scheme. 
    gmp13_bin_corretion      ! diagnostic terminal velocity for GMP13 scheme.     


!
! define module global parameters and constants in below. 
!
 real,parameter :: cpi=4.*atan(1.)

 real,parameter :: roqr=1.   ! rain density [g/cm3]
 real,parameter :: roqs=.05  ! snow density [g/cm3]
 real,parameter :: roqg=.3e0 ! graupel density [g/cm3]
 real,parameter :: roqh=0.9  ! hail density [g/cm3]
 real,parameter :: tnw=.08   ! intercept of rain [1/cm4]
 real,parameter :: tns=.10   ! intercept of snow [1/cm4]
 real,parameter :: tng=.04   ! intercept of graupel [1/cm4]

! real,parameter :: cmin=1.e-37 !minimum condensate [g/m3]
 real,parameter :: cmin=q_min_condensate !minimum condensate [g/m3]

!
! snow size mapping
!

 real,parameter :: tslopes=0.18420681 !0.08317766   ! increase tns by 8.0 from 0 to -25C !sg10
 real,parameter :: tslopeg=0.05756866   ! increase tns by 7.5 from 0 to -35C !sg9
 real,parameter :: dsnomin=0.0110 !minimum snow diameter (cm) !sg10
 real,parameter :: dgrpmin=0.0135 !minimum graupel diameter (cm) !sg9
 real,parameter :: draimax=0.0500**4. *roqr *cpi
 real,parameter :: dsnomin4=dsnomin**4.*0.900*cpi  !density at 100 microns
 real,parameter :: dgrpmin4=dgrpmin**4.*roqg *cpi
 real,parameter :: xs=0.97      !coeff [-]
 real,parameter :: sno11=1.25   !sg10
 real,parameter :: sno00=0.00
 real,parameter :: dsno11=4.00   !sg10
 real,parameter :: dsno00=1.00
 real,parameter :: sexp11=1.2
 real,parameter :: sexp00=0.6
 real,parameter :: stt=-25.    ! temperature threshold [degC]
 real,parameter :: stexp=0.50 !sg10
 real,parameter :: slim=0.80
 real,parameter :: sbase=0.240 !0.01200 !sg10

!
! graupel size mapping
!
 real,parameter :: xg=0.98
 real,parameter :: grp11=0.70
 real,parameter :: grp00=0.60
 real,parameter :: dgrp11=4.30
 real,parameter :: dgrp00=2.60
 real,parameter :: gexp11=0.50
 real,parameter :: gexp00=0.50
 real,parameter :: gtt=-20.   !temperature threshold [degC]
 real,parameter :: gtexp=0.30
 real,parameter :: glim=0.90
 real,parameter :: gbase=0.0095
 real,parameter :: roqg2=0.5
 real,parameter :: qrog2=2.0*1.e-6 !g/m**3


!
! terminal velocity constatn
!
!c***   define the coefficients used in terminal velocity
 real,parameter :: ag=330.22    !for bulk density of 0.3
 real,parameter :: bg=.36
 real,parameter :: ag2=544.83   !for bulk density of 0.5
 real,parameter :: bg2=.54

 real,parameter :: as=151.01
 real,parameter :: bs=.24
 real,parameter :: aw=2115.
 real,parameter :: bw=.8

 real,parameter :: GRVT=980.                   !hail
 real,parameter :: CD1=6.E-1                   !hail
 real,parameter :: CD2=4.*GRVT/(3.*CD1)        !hail
 real,parameter :: ah=sqrt(cd2*roqh)           !hail
 real,parameter :: bh=0.5                      !hail


! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 contains

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine gmp13_size_rho_diag (spc, tair, q, re, rhoe) 
 implicit none

!
! subroutine for intercept of particle size distributions. 
!
 character(len=*) :: spc   ! hydrometeor species
 real,intent(in)  :: tair  ! air temperature [K]
 real,intent(in)  :: q     ! mass concentration [g/m3]

 real,intent(out)  :: re   ! effective radius [micron]
 real,intent(out)  :: rhoe ! effective density [g/cm3]

 real :: tairc, taird, ftnsT, sno1, dsno1, sexp1
 real :: xx, ftnsg, ftnsQ, fexp, qsg, tnsmax
 real :: ftngT, dgrp1, grp1, gexp1, ftngQ, tngmax
 real :: intercept  !intercept  [1/cm4]  
 real :: lambda  ! lambda [1/cm]



 if (q < cmin) then
   re = 0. 
   if(trim(spc) == 'snow'   ) rhoe = 0.9   ![g/cm3]
   if(trim(spc) == 'graupel') rhoe = roqg  ![g/cm3]
   return
 endif

 !
 ! set air temperature
 !
  tairc = tair - 273.1   !air temperature [degC]
  taird=min(0.,max(-35.,tairc)+0.0)


 select case(trim(spc))
 case('snow')

  taird = min(0.,max(-25.,tairc)+0.0) ! regulate air temp [degC]
  ftnsT = exp(-1.*tslopes*taird)      ! snow intercept scaling factor in temp component [-]
  sno1  = sno11
  dsno1 = dsno11
  sexp1 = sexp11
  if(taird > stt) then !temperature dependence of coeffs
     sno1  = sno00-(sno00-sno11)    *(taird/stt)**stexp
     dsno1 = dsno00-(dsno00-dsno11)*(taird/stt)**stexp
     sexp1 = sexp00-(sexp00-sexp11)*(taird/stt)**stexp
  endif

  xx     = xs-xs*min(slim,max(0.,(q-sno1)/dsno1)**sexp1)  ! [-]
  ftnsT  = ftnsT**xx   !update ftnsT by scaling of xx  [-]
  fexp   = xx
  ftnsQ  = (q/sbase)**fexp      ! mixing ratio component [-]
  ftnsg  = ftnsT*ftnsQ             ! total component of scaling
  tnsmax = (q*1.e-6)/dsnomin4  ! maximum intercept [1/cm4] 
  if(ftnsg*tns > tnsmax) ftnsg = tnsmax/tns  !scaling parameter for snow intercept [-]

  intercept = tns * ftnsg  ! scaled intercept [1/cm4]  

!
! output effective snow density from Heymsfield et al. [2004]
! This relation is equivalent to following relationship. 
!  rhoe = 142.36 * re^(-0.95)
!  rhoe[g/cm3] and re [micron]
! 
  rhoe = min( 0.9 , 0.004168*(ftnsg*tns/ (q*1.e-6)  )**0.3115 ) ! [g/cm3]
 ! rhoe = min( 0.9 , 0.004168*(ftnsg*tns/ max(tiny(1.0),q*1.e-6)  )**0.3115 ) ! [g/cm3]

!
! output effective radius 
!
  lambda = ( intercept * cpi * rhoe / (q*1.e-6) ) ** 0.25e0  ! lambda [1/cm]
  !lambda = ( intercept * cpi * rhoe / max(tiny(1.0),q*1.e-6)  ) ** 0.25e0  ! lambda [1/cm]

 re = 3.e0 / 2.e0 / lambda * 1.e+4
!  re = 3.e0 / 2.e0 / max(tiny(1.0),lambda) * 1.e+4   !effective radius [micron]
!  re = min(huge(1.0),3.e0 / 2.e0 / max(tiny(1.0),lambda) * 1.e+4)   !effective radius [micron]


 case('graupel')

  ftngT = exp(-1.*tslopeg*taird)   ! graupel intercept scaling factor in temp component [-]
  grp1  = grp11
  dgrp1 = dgrp11
  gexp1 = gexp11
  if(taird > gtt)then
     grp1  = grp00 -(grp00 -grp11 )*(taird/gtt)**gtexp
     dgrp1 = dgrp00-(dgrp00-dgrp11)*(taird/gtt)**gtexp
     gexp1 = gexp00-(gexp00-gexp11)*(taird/gtt)**gtexp
  endif

  xx     = xg-xg*min(glim,max(0.0,(q-grp1)/dgrp1)**gexp1)
  ftngT  = ftngT**xx  !update ftnsT by scaling of xx  [-]
  fexp   = xx
  ftngQ  = 1.0
  ftngQ  = (q/gbase)**fexp
  ftnsg  = ftngT*ftngQ        ! total scaling
  tngmax = (q*1.e-6)/dgrpmin4
  if(ftnsg*tng > tngmax) ftnsg = tngmax/tng

  intercept = ftnsg*tng ! scaled intercept [1/cm4]  

  rhoe = roqg  !graupel effective density [g/cm3]

!
! output effective radius 
!
  lambda = ( intercept * cpi * rhoe / (q*1.e-6) ) ** 0.25e0  ! lambda [1/cm]
  !lambda = ( intercept * cpi * rhoe / max(tiny(1.0),q*1.e-6) ) ** 0.25e0  ! lambda [1/cm]

 re = 3.e0 / 2.e0 / lambda * 1.e+4
!  re = min(huge(1.0),3.e0 / 2.e0 / max(tiny(1.0),lambda) * 1.e+4)   !effective radius [micron]

 case default

   stop 'MSG gmp13_size_rho_diag: there is no such spc'

 end select 


 return
 end subroutine gmp13_size_rho_diag

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 subroutine gmp13_bin_corretion( tair, rh, q_cloud, q_rain, rho_dair, re )
 implicit none

!
! subroutine for intercept of rain particle size distributions due to bin correction.
!
 real,intent(in)  :: tair     ! air temperature [K]
 real,intent(in)  :: rh       ! relative humidity [%]
 real,intent(in)  :: q_cloud  ! mass concentration of clouds [g/m3]
 real,intent(in)  :: q_rain   ! mass concentration of rain [g/m3]
 real,intent(in)  :: rho_dair ! dry air density    [g/m3]
 real,intent(inout)  :: re      ! effective radius [micron]

 real :: bin_factor
 real :: ftnw 
 real :: ftnwmin
 real :: lambda
 real :: qr ! mass mixing ratio [g/g]
 real :: zrc

 qr = q_rain / rho_dair       ! mass mixing ratio [g/g] 

 if (qr < cmin .or. tair < 273.1 .or. rh >= 100. .or. q_cloud > cmin ) then
   ! do nothing
   return
 endif

 bin_factor=min( 1.30,  0.11*(1000.*qr)**(-1.27) + 0.98 )  !ranges from 0.98 to 1.30

 ftnw=1./bin_factor**3.35     !scaling factor ranges from 1.07 ~ 0.4 

!
! updating lambda (q_rain*1.e-6 is unit [g/cm3])
!
 zrc=(cpi*roqr*tnw)**0.25

 lambda  =  zrc  / ( ( q_rain*1.e-6 )**0.25 ) * ftnw**0.25  ! corrected slope (lambda) [1/cm]

 re = 3.e0 / 2.e0 / lambda * 1.e+4   !effective radius [micron]

 return
 end subroutine gmp13_bin_corretion

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

 end module module_gmp_4ice

! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP
! GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP GMP

