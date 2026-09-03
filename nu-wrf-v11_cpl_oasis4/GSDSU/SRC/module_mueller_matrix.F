!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 module module_mueller_matrix
 use module_const
 use module_tmatrix
 implicit none

!-----------------------------------------------------------------------------------------------------------
!
! Method (Calling sequence) : original
!
! driver_mueller( main driver of mueller matrix )
!     !-- get_mueller_inp1  ( get input parameters from mueller.inp1, and open output file check.out )
!     |-- zero ( zero out all output arrays )
!     |-- select_tmatrix( select the tmatrix for the particle type, open mueller scattering matrix file )
!     |     |-- band( get band labels  )
!     |-- init_orientdist_weight( derive user-input of particle orienting angle  )
!     |     |-- shotyp ( simple harmonic oscillator distribution )
!     |     |-- gtype  ( gaussian distribution )
!     |     |-- ftype  ( fisher distribution  )
!     |     |-- glquad  ( gauss-legendre quadrature weights  )
!     |     |-- dgquad  ( double gauss-legendre quadrasture weights  )
!     |     |-- lbquad  ( lobatto-quadrature weights  )
!     |     |-- qdwts   ( user-specified weights  )
!     |     |-- gauss   ( gaussian  weights  )
!     |-- zero_ass_afscat( zero out output parameters )
!     |-- init_tmatrix ( initialize constant and read in t-matrix )
!     |-- finddsd ( derive gamma DSD weight = N(D)*factor )
!     |-- zero_mueller ( zero out mueller matrix )
!     |-- thetaint ( integration of theta --> particle tilting )
!     |     |-- zero_output ( zero out output parameters )
!     |     |-- phiint ( integrate over phi--> particle rotation  )
!     |           |- scfld ( calculate scattered field for current particle orientation  )
!     |           |   |- rotate ( rotate particle orientation vector for +Z-incident vector  )
!     |           |   |- genang ( generate lab-body frame transformation angle )
!     |           |   |- addprc (calculate differential scattering cross section  )
!     |           |        |- genlgp (  )
!     |           |- getprob ( get orientation probability density function )
!     |           |- mmcalc ( integration of phi component of mueller matrix )
!     |               |- mmat ( calculate mueller matrix from scattering-amplitude matrix )
!     |-- crossec ( compute the scattering, extinction and absorption cross section )
!     |-- dsdint ( integrate over drop size distribution using weight(i) )
!     |-- store ( storing accumulated mueller matrix array )
!     |-- tormp ( integrate mueller matrix array over different hydrometeor species. )
!            |-- rmp ( calculate radar measurable parameters )
!                 |- conjugate ( conjugate covariance and mueller matrix )
!                 |- transmtx ( compute the 2x2 transmission matrix  )
!                 |- emat ( calculate extinction matrix from forward amplitude matrix F  )
!
!
! History;
!
! 9/2016   Toshi Matsui@NASA GSFC : Complete cleaning up, by removing old subroutines. 
! 3/2016   Toshi Matsui@NASA GSFC : re-orgniazed code before adapting G-SDSU.
! 1/2016   Toshi Matsui@NASA GSFC : modified the entire code into F90 format.
!  1993    John Beaver, John Hubbert, Larry Carey @CSU ; additional development
!  1990    Vivekanandan@CSU       : initial development.
!
!
! References:
!
!  Vivekanandan, J., W. M. Adams, and V. N. Bringi, 1991: Rigorous Approach to Polarimetric Radar Modeling 
!    of Hydrometeor Orientation Distributions. J. Appl. Meteor., 30, 1053-1063. 
!
!-----------------------------------------------------------------------------------------------------------

!
! Encapsulation control 
!
  private   ! encapsulate all module parameters and subourtines.


!
! public (accessible) subroutine
!
  public :: driver_mueller         , & ! drive mueller matrix
            derive_output_polradar     ! derive polarimetric radar measurables

!
! floating point
!
 integer,parameter :: fps = SELECTED_REAL_KIND( 6, 37)  !single precision precision kind
 integer,parameter :: fpd = SELECTED_REAL_KIND(15,307)  !double precision precision kind

!
! control debuggig  (for develoopment --> true )
!
 logical,parameter :: debug_mueller = .false.

!
! io device number
!
 integer,parameter :: io_check_out    =  9
 integer,parameter :: io_scatmat      = 20
 integer,parameter :: io_radarobs     = 18

!
! new dimension index for ensemble parameters 
!
! integer,parameter,public :: max_stdev = 10
 integer,parameter,public :: max_ens = 10

 integer,parameter,public :: max_elv = 91  !consistent to LUT 46


! file name
 character(len=200) :: io_file_append


! this is for T-matrix data
!mtxcom
 integer :: nrank, nranki
 complex :: tmat(12,20,40,40)   ! toshii this must be store for all species x 14. --> add 12 for size bin

!fnccom
 real(fpd) :: pnmllg(41)

!cmvcom
 integer :: nm, kmv
 real(fpd) :: cmi(20), cmv, cm2, twm, prodm

!thtcom
 real(fpd) :: theta,sinth,costh

!uvccom
 real(fpd) :: angint, anginp, utheta, uphi, rtsfct 
 integer :: ip  !1:     2: orthogonal polarization

!angcom
 real(fpd) :: tang(21), pang(21), cosb(21), sinb(21), scv, sch, exv, exh
 complex  :: acans(21,2,2)   !scattered field (2nd dimension 1:theta polarized scattering field  2:phi polarized scatterinig field)

!main01
 real(fpd) :: abv, abh, conk, lambda, schpt, scvpt, exhpt, exvpt, phnorm

!main01a
 real(fpd) :: s(21,4,4), ss(21,4,4) !, exmpt(4,4)  !ss is 4x4 scattering (mueller) matrix
 complex :: fscat(2,2), fscatt(2,2)   !2x2 forward scattering matrix

!main01b
 real(fpd) :: scvp,schp,exvp,exhp

!main02
 integer,parameter :: ntheta = 20 ! order of integration over theta (this is also consistent for Langevin Unimodal)
 integer,parameter :: nphi   = 20 ! order of integration over phi   (this is also consistent for Langevin Unimodal)
! integer :: ntheta, nphi
 real(fpd) :: tasc(20), twt(20), pasc(20), pwt(20)

!main03
 integer :: distyp  ! distribution type
 real(fpd) :: thetam, pnorm, mean, sigma, tcl, tcu

!main04
 integer :: nang, nuang, check, iang
 real(fpd) :: oldtht, oldphi

!main05
 character(len=20) :: fname

!main06
 real(fpd) :: eastrt, eastop, eainc, elvang, kappa, theta0, phi0
 integer :: nelang  !number of elevation angle

! const
 real(fpd),parameter :: degrad = 0.017453292519943
 real(fpd),parameter :: raddeg = 57.295779513
 real(fpd),parameter :: pi     = 3.14159265358980


!infile
 integer,parameter  :: npart = 1!# of particle types
 character(len=80) :: pathname

!infile2
 integer :: pathlen

!params
 real(fpd) ::  mugrp(max_elv,21,4,4), ass(max_ens,21,4,4)
 complex :: afscat(max_ens,2,2),fscatp(max_elv,2,2)

!dsd
 real(fpd) :: n0(180), d0(180), mgam(180), weight
 
!dsd1
 real(fpd) :: dmin, dmax, dstep

!dsd3
 real(fpd) :: dpar1, dpar2, dpar3

!trans1 (fixed...)
 integer,parameter :: transflag  =    1    ! compute transmission matrix and out in file (1=yes,0=no)
 real(fpd),parameter :: dist_trans = 1000. ! distance used to compute the transmission matrix (in meters)

!trans2
 character(len=25) :: medium, orient, orientmd(10), mediummd(10)
 character(len=2) ::  mdfile

!getphi
 integer,parameter :: phiflag = 0   ! 0 is ramdom phi distributions; 1 is fixed only for langevin distributions
 real(fpd),parameter :: phiang = 90. ! fixed phi angle (not used usually)
 
!ang
 real(fpd) :: eldeg
 integer,parameter :: nmu=4

!nang = nmu * 2 = 8  (1 is 0deg, each deg increment 45 deg)
!backscatteirng angle is nmu 5 -->  (5-1)*45 = 180
! 
 integer,parameter :: nback = 5

!int
 real(fpd) :: qmu(nmu),qwt(nmu),qmue(2*nmu),qwte(2*nmu)

!
! outfile
!
 real(fpd) :: dpar1md(10),dpar2md(10),dpar3md(10),dminmd(10),dmaxmd(10),dstepmd(10)

 real(fpd) :: wx(0:1000)

! toshii -- special for new modification
 integer :: ndset_global
 integer :: nds_global


!
! Atmosphere parameters at staggered level, dimension(is:ie,js:je,ks-1:ke) 
!
 public :: polradar_parameter
 type polradar_parameter
    real(fpd) :: zhh      ! reflectivity [dBZ] 
    real(fpd) :: zdr      ! differential reflectivity [dB]
    real(fpd) :: kdp      ! specific diff phase [deg/km] 
    real(fpd) :: rhohv    ! copolar (hv) correlation coefficient [-]
    real(fpd) :: ldr      ! linear depol ratio [dB]
    real(fpd) :: zdp      ! difference reflectivity [dB]
    real(fpd) :: delhv    ! phaseshift upon backscatter [deg]
    real(fpd) :: datt     ! differential attenuation [dB/km]
    real(fpd) :: atth     ! specific hor attenuation [dB/km]
    real(fpd) :: attv     ! specific ver attenuation [dB/km]
    real(fpd) :: degpolh  ! deg. polarization h
    real(fpd) :: degpolv  ! deg. polarization v
!    complex   :: rhoh     ! deg. coherence h
!    complex   :: rhov     ! deg. coherence v
 end type polradar_parameter

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 contains

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine driver_mueller ( axis_ratio, dcreal, dcimag, diam_sph, dwave, &       
                             ityp_dist, mean_orient, num_stdev, stdev_orient , &
                             tmat_out, conk, sigma_temp ,mxelev_in, elev_in , &
                             mueller_mat_out, fscat_mat_out  )
 implicit none
! original driver of mueller matrix
!
!     ... .................................................................
!     CALCULATE MUELLER MATRIX (AVERAGED STOKES MATRIX) FOR A DISTRIBUTION
!     OF SCATTERERS WITH A PRESCRIBED ORIENTATION DISTRIBUTION.
!
!     ALSO CALCULATES THE RADAR OBSERVABLES FROM THE MUELLER MATRICES
!     FOR DIFFERENT TYPES (SPECIES) OF PARTICLES. THE MUELLER MATRICES
!     OBTAINED ARE INTEGRATED OVER THE SPECIFIED DSD FOR EACH TYPE OF
!     PARTICLE. THE REQUIRED ELEVATION (LOOK) ANGLES, DSD CAN BE SPECIFIED
!     IN THE PROGRAM AS ALSO THE PARTICLE TYPE AND THE DSD USED FOR EACH TYPE.
!     RAVI: 12/30/91
!
!     Note:  jb=John Beaver; jh=John Hubbert; ldc=Larry Carey
!     ... .................................................................
!     *** CAVEAT: MAX PARMATERS THIS PROGRAM CAN HANDLE DUE
!     ***         TO MEMORY LIMITATIONS ARE : NDSD = 180
!     ***                                     NELANG = 10
!     RAVI: 12/30/91
!     In this version of the mueller code the NDSD loop has been
!     removed -- to save memory and allow the program to be run on
!     all sun machines. Previously, only Everest had enough memory
!     to run the program.                                    jb 3/9/95
!----------------------------------------------------------------------------

 real(fpd),intent(in) :: axis_ratio    ! particle axis ratio [-]
 real(fpd),intent(in) :: dcreal,dcimag ! real component of effective permittivity [-]
 real(fpd),intent(in) :: diam_sph      ! spherical diameter [cm]
 real(fpd),intent(in) :: dwave         ! wavelength [cm] 

 integer,intent(in) :: ityp_dist  ! particle orienting angle distribution type
 real(fpd),intent(in) :: mean_orient! mean orienting angle [deg]
 integer,intent(in) :: num_stdev  ! number of orienting angle standard deviations
 real(fpd),intent(in) :: stdev_orient(:)  ! standard deviation / theta_m of orienting angle [deg]

 complex,intent(in)   :: tmat_out(20,40,40)   ! tmatrix [-]
 real(fpd),intent(in) :: conk, sigma_temp     ! 
 integer,intent(in) :: mxelev_in              ! maximum dimension of elevation angle 
 real(fpd),intent(in) :: elev_in( mxelev_in ) ! elevation angles [deg]

 real(fpd),intent(out) :: mueller_mat_out (max_ens,max_elv,4,4) ! 4x4 mueller matrix (backscattering component only)
 complex,intent(out)   :: fscat_mat_out   (max_ens,max_elv,2,2) ! 2x2 forward scattering matrix 

 integer :: n_htyp     !loop index for particle type
 integer :: nelev      !loop index for elevation angle
 integer :: nds ! # of input t-matrix
 integer :: ndset
 integer :: norient, nstdev

 if(debug_mueller) print*,'MSG driver_mueller: start mueller matrix routine'

!
! initialize output
!
 mueller_mat_out = 0. !4x4 mueller matrix 
 fscat_mat_out   = cmplx(0.,0.) !2x2 forward scattering matrix 

 n_htyp = 1  ! always single hydrometeor type


!
! loop for standard deviation (or theta_m)
!
 loop_stdev: do nstdev = 1, num_stdev
   if(debug_mueller) print*,'nstdev=',nstdev,num_stdev

  !
  ! Distribiton types (0. Gaussian, 1: Harmonic Oscillator, 2. 2sigma-truncated Gaussian, 3 Langevin Unimodal)
  ! 
  distyp = ityp_dist   
  dpar1  = mean_orient           ! mean orienting angle [deg]
  dpar2  = stdev_orient(nstdev)  ! standard deviation of orienting angle  [deg]
  dpar3  = 00.00               ! mean phi (only for dtyp1=3)

  !
  ! initialize output parameters
  !
  mugrp(:,:,:,:) = 0.     !mueller matrix for diffent elevation angles
  fscatp(:,:,:)  = cmplx(0.,0.)

  !
  ! transfer some parameters (orientation distribution type) and creating output name 
  !
   nds = 1
   if(debug_mueller) print*,'*** the # of t-matrices to be processed are ',nds



  !
  ! Calculate weight of particle orientation distributions for a given osccilation type/function.
  !
  if(debug_mueller) print*,'init_orientdist_weight'
  call init_orientdist_weight

  !
  ! initialie ass and afscat and other parameters
  !
  ass(:,:,:,:)=0.
  afscat(:,:,:)=cmplx(0.0)  !not used anymore
  ndset = 1
  ndset_global = ndset
  iang=1

  !
  ! Transfer tmatrix parameters into global parameters
  !
  if(debug_mueller) print*,'init_tmatrix'
  call init_tmatrix ( axis_ratio, dcreal, dcimag, diam_sph, dwave, &
                     tmat_out, conk, sigma_temp )

  !
  ! elevation loop
  !
  loop_elevation: do nelev = 1, mxelev_in  

    elvang = elev_in( nelev ) * degrad  ! elevation angle [rad] 
    if(debug_mueller) print*,'elevation angle [deg] = ', elev_in( nelev )

    !
    ! initialize mueller and scattering matrix and ext/scat/abs cross sections
    !
    scvpt = 0.0 ; schpt = 0.0 ; exvpt = 0.0 ; exhpt = 0.0
    ss(:,:,:) = 0.0  ! 4x4 mueller matrix [m2] 
    fscatt(:,:) = cmplx(0.,0.)  !2x2 forward scatttering matrix 

    !
    ! integration of theta and phi in mueller matrix (this is MAJOR subroutine)
    !
    if(debug_mueller) print*,'thetaint'
    call thetaint

    !
    ! compute the scattering, extinction and absorption cross sections (short)
    !
    if(debug_mueller) print*,'crossec'
    call crossec

    !
    ! parameters transfered for tormp
    !
    mugrp(nelev,1:nang,1:4,1:4) = ss(1:nang,1:4,1:4)
    fscatp(nelev,1:2,1:2) = fscatt(1:2,1:2)

    !
    ! output bundled mueller (backscatter component only) and scattering matrices
    !
    mueller_mat_out(nstdev,nelev,1:4,1:4) =  ss(nback,1:4,1:4)  !4x4 mueller matrix 
    fscat_mat_out  (nstdev,nelev,1:2,1:2) = fscatt   (1:2,1:2)  !2x2 forward scattering matrix 

  enddo loop_elevation


!
! computer polarimetric radar observables --> this will be placed in module_driver_polarris.F
! total_opt subroutine
!
 if(debug_mueller) then
    print*,'tormp'
    call tormp( mxelev_in, elev_in, nang )
 endif


 enddo loop_stdev



 return
 end subroutine driver_mueller

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine band(lambda,fband)
 implicit none

!
! get radar band 'label' based on lambda
!

 real(fpd) :: lambda
 character(len=2) :: fband

!
 lambda = lambda*.01
 if(lambda .ge. .3 .and. lambda .lt. 1)then
    fband='pb'
 elseif(lambda .ge. .15 .and. lambda .lt. .3 )then
    fband='lb'
 elseif(lambda .ge. .075 .and. lambda .lt. .15 )then
    fband='sb'
 elseif(lambda .ge. .0375 .and. lambda .lt. .075 )then
    fband='cb'
 elseif(lambda .ge. .024 .and. lambda .lt. .0375 )then
    fband='xb'
 elseif(lambda .ge. .0167 .and. lambda .lt. .024 )then
    fband='ku'
 elseif(lambda .ge. .0110 .and. lambda .lt. .0167 )then
    fband='kb'
 elseif(lambda .ge. .0075 .and. lambda .lt. .0110 )then
    fband='ka'
 endif
 if(debug_mueller) print*,'MSG band: band name=',trim(fband)

 return
 end subroutine band

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine init_orientdist_weight 
 implicit none

!
! initialize weight of specific orienting-angle distributiions function. 
!
 integer :: i
 character(len=1) :: qdtype
 character(len=12) :: quadnam
 real(fpd) :: a,b,lower,upper,lowerp,upperp

!
! derive input of various orienting angle distribution functions
!
 if (distyp .eq. 0) then
   ! RANDOM DISTRIBUTION
   lower = -1.0
   upper = 1.0
   pnorm = 1.0
 else if (distyp .eq. 1) then
   ! SIMPLE HARMONIC OSCILLATOR:
   call shotyp (mean,thetam,pnorm,lower,upper,check,tcl,tcu)  
 else if (distyp .eq. 2) then
   ! QUASI-GAUSSIAN: (lower,upper is derived from user.in)
   call gtype (mean,sigma,pnorm,lower,upper,check,tcl,tcu)
 else if (distyp .eq. 3) then
   ! FISHER:
   call ftype (kappa,theta0,phi0,lower,upper,lowerp, pnorm,upperp,check)
   ! PNORM IS SET TO ONE FOR THIS DISTRIBUTION
   ! The fisher dis. is nearly one, .998 depending on mean values and variance
   ! used. Better to use phnorm.
   !JH      pnorm = 1.0
 endif

!
! specify quadrature type --> currently always latatto-quad 
! (toshii--> this will change output significantly)
!
 qdtype='l'

!
! SET THE NAME OF THE QUADRATURE ACCORDING TO SELECTION
!
 if(qdtype.eq.'g') then
   quadnam='gauss-legend'
 elseif(qdtype.eq.'d') then
   quadnam='double-gauss'
 elseif(qdtype.eq.'l') then
   quadnam='lobatto-quad'
 else
   quadnam='user specify'
 endif

!
! MAKE THE DESIRED QUADRATURE ABSCISSAS AND WEIGHTS
!
 if (qdtype(1:1) .eq. 'g') then
   call glquad (nmu, qmu, qwt) 
 elseif (qdtype(1:1) .eq. 'd') then
   call dgquad (nmu, qmu, qwt)
 elseif (qdtype(1:1) .eq. 'l') then
   call lbquad (nmu, qmu, qwt)
 else
   do i = 1, nmu
      stop 'MSG : You must specify qmu(i) '
   enddo
   call qdwts (nmu, qmu, qwt)
 endif


!
! specify mueller matrix angle bin 
!
 nang = 2*nmu
 nuang = nang

!
! get weights for main program integration
!
 if (distyp .ne. 3) then
    ! NOT LANGEVIN:
    lowerp=0.
    upperp=2*pi
! over ride nphi for LANGEVIN
!     nphi=20
!     ntheta=20
 endif

 a=acos(upper)*raddeg
 b=acos(lower)*raddeg

 if( debug_mueller ) then
   write(6,*) 'n lo up th',ntheta,a,b
   write(6,*) '       phi',nphi,lowerp*raddeg,upperp*raddeg
 endif


!
! get weight and point of gauss 
!
 call gauss(twt,tasc,ntheta,lower,upper)
 call gauss(pwt,pasc,nphi,lowerp,upperp)


! 900  format (1x,f7.2)
! 901  format (1x,i2)
! 903  format (a20)
! 904  format (1x,a20)
! 905  format (1x,'elevation angle start, stop, increment = ', &
!              2(f7.2,', '),f7.2)
! 906  format (1x,'theta distribution type = ',i1,/,' (0 = random;', &
!            /,'  1 = simple harmonic oscillator;',/,'  2 = gaussian)')
! 907  format (1x,'cos(theta) ranges from ',f7.2,' to ',f7.2)
! 908  format (1x,'integration orders:  theta -- ',i2,'; phi -- ',i2)
 !909  format (/,1x,'random distribution in theta:')

 return
 end subroutine init_orientdist_weight

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine shotyp (mean,thetam,pnorm,lower,upper,check,tcl,tcu)
 implicit none

!  GET INPUT ON SIMPLE HARMONIC OSCILLATOR DISTRIBUTION

 integer :: check
 real(fpd) :: mean, thetam, pnorm, lower, upper, tcl, tcu
    
 real(fpd) :: twt(20), tasc(20), pwt(20), pasc(20), theta, phi
 real(fpd) ::  prob, psum, tsum
 integer :: pcount, tcount
 
!
!
!     GET MEAN, MAXIMUM AMPLITUDE OF THE OSCILLATION
!     WRITE (6,'(\A)') 'ENTER MEAN THETA OF OSCILLATION:'
!     TYPE *, 'ENTER MEAN THETA OF OSCILLATION:'
! *** READ (4,*) MEAN
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91
      mean = dpar1*degrad
!     WRITE (6,900) MEAN*RADDEG
!
!     TYPE *, 'ENTER AMPLITUDE OF OSCILLATION (DEGREES): '
! *** READ (4,*) THETAM
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91
      thetam = dpar2*degrad
!     WRITE (6,900) THETAM*RADDEG

!     NORMALIZE THE DISTRIBUTION OVER GIVEN THETA RANGE
      tcl = mean + thetam
      if (tcl .gt. pi) then
         tcl = 2*pi  - tcl
         lower = -1.00
      else
         tcl =  pi
         lower = cos((mean + thetam))
      endif

      tcu = mean - thetam
      if (tcu .lt. 0.0) then
         tcu = -tcu
         upper = 1.00
      else
         tcu = 0.00
         upper = cos((mean - thetam))
      endif

      call gauss (twt, tasc, 20, lower, upper)
      call gauss (pwt, pasc, 20, 0.0d0, 2*pi )

      psum = 0.0
      do pcount = 1,20
         if(phiflag .eq. 0)then
            phi = pasc(pcount)
         else
            phi = phiang*degrad
         endif
         tsum = 0.0
         do tcount = 1,20
            theta = 2.0*atan(sqrt(1.0-tasc(tcount)**2)/ (1.0 + tasc(tcount)))

            prob = shop(theta,mean,thetam,phi,1.00d0)
!           ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE
            if (theta .lt. tcu) then
               prob = prob + shop(-theta,mean,thetam,(phi+pi),1.00d0)
            else if (theta .gt. tcl) then
               prob = prob + shop((2*pi -theta),mean,thetam,(phi+pi),1.00d0)
            endif
            tsum = tsum + prob*twt(tcount)
         enddo
         psum = psum + tsum*pwt(pcount)
      enddo
      pnorm = psum

!     MAKE SURE NORMALIZATION HAS BEEN DONE PROPERLY
      psum = 0.0
      do pcount = 1,20
         if(phiflag .eq. 0)then
            phi = pasc(pcount)
         else
            phi = phiang*degrad
         endif
         tsum = 0.0
         do tcount = 1,20
            theta = 2.0*atan(sqrt(1.0-tasc(tcount)**2)/ (1.0 + tasc(tcount)))

            prob = shop(theta,mean,thetam,phi,pnorm)
!           ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE
            if (theta .lt. tcu) then
               prob = prob + &
                      shop(-theta,mean,thetam,(phi+pi),pnorm)
            else if (theta .gt. tcl) then
               prob = prob + &
                      shop((2*pi -theta),mean,thetam,(phi+pi),pnorm)
            endif
            tsum = tsum + prob*twt(tcount)
         enddo
         psum = psum + tsum*pwt(pcount)
      enddo

!     OUTPUT TO DIAGNOSTIC FILE
if(debug_mueller) then
      write(io_check_out,901)
      write(io_check_out,902) mean*raddeg, thetam*raddeg
      write(io_check_out,903) psum
endif

! 900  format (1x,f7.2)
 901  format (/,1x,'simple harmonic oscillator distribution in theta:')
 902  format (1x,'mean theta = ',f7.2,'; amplitude of oscillation = ', &
              f7.2)
 903  format (1x,'normalization of p(theta,phi) yields ',f6.4)

 return
 end subroutine shotyp 

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine gtype (mean,sigma,pnorm,lower,upper,check,tcl,tcu)
 implicit none

!
! GET USER INPUT ON GAUSSIAN DISTRIBUTION FUNCTION
!

 integer :: check, pcount, tcount
 real(fpd) :: lower, upper, mean, pnorm, prob, sigma
 real(fpd) :: twt(20), tasc(20), pwt(20), pasc(20), &
              theta, phi
 real(fpd) :: tcl, tcu, psum, tsum

!     GET MEAN, SIGMA OF THE OSCILLATION
!     WRITE (6,*) 'ENTER MEAN THETA OF OSCILLATION:'
! *** READ (4,*) MEAN
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI:12/26/91
 mean = dpar1*degrad
!     WRITE (6,900) MEAN*RADDEG
!     PRINT*, 'MEAN = ',MEAN
!     TYPE *, 'ENTER SIGMA OF GAUSSIAN (DEGREES): '
! *** READ (4,*) SIGMA
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI:12/26/91
 sigma = dpar2*degrad


!     PRINT*, 'SIGMA = ',SIGMA*RADDEG
!     WRITE (6,900) SIGMA

!     NORMALIZE THE DISTRIBUTION OVER GIVEN THETA RANGE
 tcl = mean + 2.0*sigma

 if (tcl .gt. pi/2.) then
         tcl = pi/2.
!jh      LOWER = -1.00
         lower = 0.00
 else
!jh      TCL = PI
         lower = cos((mean + 2.0*sigma))
 endif

 tcu = mean - 2.0*sigma
 if (tcu .lt. 0.0) then
         tcu = 0.0
         upper = 1.00
 else
!        TCU = 0.00
         upper = cos((mean - 2.0*sigma))
 endif

      call gauss (twt, tasc, 20, lower, upper)
      call gauss(pwt, pasc, 20, 0.0d0, 2*pi )
!     NORMALIZE THE DISTRIBUTION OVER GIVEN THETA RANGE
      psum = 0.0
      do 200 pcount = 1,20
         if(phiflag .eq. 0)then
            phi = pasc(pcount)
         else
            phi = phiang*degrad
         endif
         tsum = 0.0
         do 100 tcount = 1,20
            theta = 2.0*atan(sqrt(1.0-tasc(tcount)**2)/ &
                    (1.0 + tasc(tcount)))

            prob = gp(theta,mean,sigma,phi,1.00d0)

!           ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE
!jh         IF (THETA .LT. TCU) THEN
!jh            PROB = PROB +
!jh  >                GP(-THETA,MEAN,SIGMA,(PHI+PI),1.00)
!jh         ELSE IF (THETA .GT. TCL) THEN
!jh            PROB = PROB +
!jh  >                GP((2*PI -THETA),MEAN,SIGMA,(PHI+PI),1.00)
!jh         ENDIF
            tsum = tsum + prob*twt(tcount)
 100     continue
         psum = psum + tsum*pwt(pcount)
 200  continue
      pnorm = psum

!     MAKE SURE NORMALIZATION HAS BEEN DONE PROPERLY
      psum = 0.0
      do 400 pcount = 1,20
         if(phiflag .eq. 0)then
            phi = pasc(pcount)
         else
            phi = phiang*degrad
         endif

         tsum = 0.0
         do 300 tcount = 1,20
            theta = 2.0*atan(sqrt(1.0-tasc(tcount)**2)/ &
                    (1.0 + tasc(tcount)))

            prob = gp(theta,mean,sigma,phi,pnorm)
!           ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE
!jh         IF (THETA .LT. TCU) THEN
!jh            PROB = PROB +
!jh  >                GP(-THETA,MEAN,SIGMA,(PHI+PI),PNORM)
!jh         ELSE IF (THETA .GT. TCL) THEN
!jh            PROB = PROB +
!jh  >                GP((2*PI -THETA),MEAN,SIGMA,(PHI+PI),PNORM)
!jh         ENDIF
            tsum = tsum + prob*twt(tcount)
 300     continue
         psum = psum + tsum*pwt(pcount)
 400  continue

if( debug_mueller ) then
      write(6,901)
      write(6,902) mean*raddeg, sigma*raddeg
      write(6,903) pnorm
      write(6,904) psum
!     OUTPUT TO DIAGNOSTIC FILE
      write(io_check_out,901)
      write(io_check_out,902) mean*raddeg, sigma*raddeg
      write(io_check_out,903) pnorm
      write(io_check_out,904) psum
endif

      return
! 900  format (1x,f7.2)
 901  format (/,1x,'gaussian distribution in theta:')
 902  format (1x,'mean theta = ',f7.2,'; sigma = ',f7.2)
 903  format (1x,'area under surface = ',f6.3)
 904  format (1x,'normalization of p(theta,phi) yields ',f6.4)

 return
 end subroutine gtype

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine ftype (kappa,theta0,phi0,lower,upper,lowerp, &
                        pnorm,upperp,check)
 implicit none

!   this subroutine defines the integration limits for theta and
!   phi.
!
!  This function will always integrate to 1 and therefore PNORM
!  will be one. Thus this subroutine serves as a check for the
!  program, ie if the integration does not yield one then there
!  is a problem in the software.
!                                 john hubbert 20 oct 1993
!
 integer :: check, pcount, tcount
 real(fpd) :: lower, upper, kappa, theta0, phi0, &
      pnorm, prob, upperp, lowerp
 real(fpd) :: twt(20), tasc(20), pwt(20), pasc(20), &
      theta, phi
 real(fpd) :: psum, tsum, phicorect, summ, tmp

!
!
!     GET KAPPA, THETA0, AND PHI0 OF THE OSCILLATION
!     WRITE (6,*) 'ENTER KAPPA:'
! *** READ (4,*) KAPPA
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91
      kappa = dpar1
!     WRITE (6,*) 'ENTER THETA0 OF DISTRIBUTION:'
! *** READ (4,*) THETA0
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91
      theta0 = dpar2*degrad
!     WRITE (6,*) 'ENTER PHI0 OF DISTRIBUTION:'
! *** READ (4,*) PHI0
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91
      phi0 = dpar3*degrad


! For narrow distributions, integrating over the entire
! sphere can cause integration problems unless the integration
! intervals are reduced. Thus theta is integrated over 3
! "standard deviations" and phi will be intgrated over
! the equivalent of 3*sigma.
! This integration routine and the so defined limits were
! tested up to kappa=7000 (sigma is about 1 degree) for
! any mean orientation angles. The integration always yielded
! .9995 and  more typically higher.

!  calculate equivalent sigma

      sigma = 1.08*2*(1/(sqrt(2*kappa)))
      if(debug_mueller) print *,'equivalent sigma  for fisher dis. is:',sigma*raddeg
      if(debug_mueller) write(io_check_out,111)sigma*raddeg
 111   format('equivalent sigma  for fisher dis. is:',f7.3)
!
! DETERMINE LIMITS
!   theta:

      if(theta0 + 3*sigma.gt. pi) then
         upper= pi
      else
         upper = (theta0 + 3*sigma)
      endif

      if(theta0 - 3*sigma.lt. 0.0) then
         lower  = 0.0
      else
         lower = theta0 - 3*sigma
      endif

! transform integration limits by u=cos(theta) and
! switch integration limits

      tmp=upper
      upper=cos(lower)
      lower=cos(tmp)

! phi:
!       correct phi limits as a function of theta

      phicorect = atan(tan(3*sigma)/abs(sin(theta0)))
      lowerp = phi0 - phicorect
      upperp = phi0 + phicorect
      if(upperp - lowerp .gt. 2*pi) then
         lowerp = 0.0
         upperp = 2*pi
      endif
      if (theta0 - 3*sigma .lt. 0.0) then
         lowerp = 0.0
         upperp = 2*pi
      endif

!  get weights and points

      call gauss (twt, tasc, 20, lower, upper)


      call gauss (pwt, pasc, 20, lowerp, upperp)

! integrate

      psum = 0.0
      do pcount = 1,20
         phi = pasc(pcount)
         tsum = 0.0
         do tcount = 1,20
            theta = 2.0*atan(sqrt(1.0-tasc(tcount)**2)/ &
                    (1.0 + tasc(tcount)))

            prob = fp(kappa, theta, phi, theta0, phi0, 1.00d0)
            tsum = tsum + prob*twt(tcount)
         enddo
         psum = psum + tsum*pwt(pcount)
      enddo

      pnorm = psum
      if(debug_mueller) write(6,*) 'area under fisher curve:',pnorm

!     MAKE SURE NORMALIZATION HAS BEEN DONE PROPERLY
!
       summ=0.0

       psum = 0.0
       do pcount = 1,20
          phi = pasc(pcount)
          tsum = 0.0
          do tcount = 1,20
             theta = 2.0*atan(sqrt(1.0-tasc(tcount)**2)/ &
                    (1.0 + tasc(tcount)))

             prob = fp(kappa, theta, phi, theta0, phi0, pnorm)
             tsum = tsum + prob*twt(tcount)
          enddo
          psum = psum + tsum*pwt(pcount)
       enddo
       summ=summ+psum

       if(debug_mueller) write(6,*) 'normalized fisher integrates to: ',summ
       psum=summ

!     OUTPUT TO DIAGNOSTIC FILE
if(debug_mueller) then
      write(io_check_out,901)
      write(io_check_out,902) kappa, theta0*raddeg, phi0*raddeg
      write(io_check_out,903) pnorm
      write(io_check_out,904) psum
endif

      return
! 900  format (1x,f7.2)
 901  format (/,1x,'langevin unimodal distribution:')
 902  format (1x,'kappa = ',f4.1,'; theta0 = ',f7.2, &
              '; phi0 = ',f7.2)
 903  format (1x,'area under surface = ',d12.5)
 904  format (1x,'normalization of p(theta,phi) yields ',f6.4)

 return
 end subroutine ftype

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine glquad (num, asbcis, wts)
 implicit none
!      GENERATES THE ABSCISSAS AND WEIGHTS FOR AN EVEN 2*NUM POINT
!      GAUSS-LEGENDRE QUADRATURE.  ONLY THE NUM POSITIVE POINTS ARE RETURNED.

 integer,intent(in) ::  num
 real(fpd),intent(out) :: asbcis(num), wts(num)

 integer ::  n, i, j, k
 real(fpd) :: z, z1, p1, p2, p3, pp
 real(fpd),parameter :: eps=3.0e-14
!
      n=2*num
      do 12 i=1,num
        z=cos(3.14159265400*(i-.2500)/(n+.500))
        k=0
 1       continue
          p1=1.00
          p2=0.00
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.00*j-1.00)*z*p2-(j-1.00)*p3)/j
 11        continue
          pp=n*(z*p1-p2)/(z*z-1.00)
          z1=z
          z=z1-p1/pp
          k=k+1
        if (abs(z-z1).gt.eps .and. k.lt.10) go to 1
        asbcis(num+1-i)=z
        wts(num+1-i)=2.00/((1.00-z*z)*pp*pp)
 12    continue
!
 return
 end subroutine glquad

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 


 subroutine dgquad (num, asbcis, wts)
 implicit none

!        GENERATES THE ABSCISSAS AND WEIGHTS FOR AN NUM POINT
!      GAUSS-LEGENDRE QUADRATURE BETWEEN 0 AND 1.

 integer,intent(in)::  num
 real(fpd) ::  asbcis(num), wts(num)

 integer ::  n, i, j, m, k
 real(fpd) ::   z, z1, p1, p2, p3, pp
 real(fpd), parameter :: eps=3.00-14

 n=num
 m=(n+1)/2

 do i=1,m
    z=cos(3.14159265400*(i-.2500)/(n+.500))
    k=0
 1       continue
    p1=1.00
    p2=0.00
    do j=1,n
       p3=p2
       p2=p1
       p1=((2.00*j-1.00)*z*p2-(j-1.00)*p3)/j
    enddo

    pp=n*(z*p1-p2)/(z*z-1.00)
    z1=z
    z=z1-p1/pp
    k=k+1
    if(abs(z-z1).gt.eps .and. k.lt.10) go to 1
    asbcis(i)=.500-.500*z
    asbcis(n+1-i)=.500+.500*z
    wts(i)=1.00/((1.00-z*z)*pp*pp)
    wts(n+1-i)=wts(i)
 enddo

 return
 end subroutine dgquad


!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine lbquad (num, asbcis, wts)
 implicit none
!        GENERATES THE ABSCISSAS AND WEIGHTS FOR AN EVEN 2*NUM POINT
!      LOBATTO QUADRATURE.  ONLY THE NUM POSITIVE POINTS ARE RETURNED.

 integer,intent(in) :: num
 real(fpd),intent(out) :: asbcis(num), wts(num)
 integer ::  n, n1, i, j, k
 real(fpd) :: z, z1, p1, p2, p3, pp, ppp, ci
 real(fpd), parameter :: eps=3.00-14

  n=2*num
  n1=n-1
  ci = 0.50
  if (mod(n,2) .eq. 1)  ci = 1.00

  do i=1,num-1
     z=sin(3.14159265400*(i-ci)/(n-.5))
     k=0

 1   continue
     p1=1.00
     p2=0.00
     do j=1,n1
        p3=p2
        p2=p1
        p1=((2.00*j-1.00)*z*p2-(j-1.00)*p3)/j
     enddo
     pp=n1*(z*p1-p2)/(z*z-1.00)
     ppp=(2.00*z*pp-n1*(n1+1)*p1)/(1.00-z*z)
     z1=z
     z=z1-pp/ppp
     k=k+1
     if(abs(z-z1).gt.eps .and. k.lt.10) go to 1
     asbcis(i)=z

     wts(i)=2.00/(n*(n-1)*p1*p1)
  enddo

  asbcis(num)=1.00
  wts(num)=2.00/(n*(n-1))

 return
 end subroutine lbquad

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine qdwts (num, asbcis, wts)
 implicit none

 integer,intent(in) :: num
 real(fpd),intent(in) ::  asbcis(num)
 real(fpd),intent(out) ::  wts(num)

 real(fpd) :: x(100), alphas(100)
 integer :: i

!
 do i = 1, num
    x(i) = asbcis(i)**2
 enddo

 do i = 1, num
    alphas(i) = 1.00/(2*i-1)
 enddo

 call vander (x, wts, alphas, num)
!
 return
 end subroutine qdwts 

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine vander(x,w,q,n)
 implicit none

 integer,intent(in) :: n

!     IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
 integer, parameter :: nmax=100
 real(fpd),parameter :: zero=0.000,one=1.000
 real(fpd) :: x(n),w(n),q(n),c(nmax)
 real(fpd) :: s,b,t,xx
 integer :: i,j,k,k1

 if(n.eq.1)then
    w(1)=q(1)
 else

    do i=1,n
       c(i)=zero
    enddo
    c(n)=-x(1)
    do i=2,n
       xx=-x(i)
       do j=n+1-i,n-1
          c(j)=c(j)+xx*c(j+1)
       enddo
       c(n)=c(n)+xx
    enddo

    do i=1,n
       xx=x(i)
       t=one
       b=one
       s=q(n)
       k=n
       do j=2,n
          k1=k-1
          b=c(k)+xx*b
          s=s+q(k1)*b
          t=xx*t+b
          k=k1
       enddo
       w(i)=s/t
    enddo
 endif

 return
 end subroutine vander

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine gauss (wt,asc,n,aa,bb)
 implicit none

 integer,intent(in) :: n
 real(fpd),intent(in) :: aa,bb
 real(fpd),intent(out) :: wt(n),asc(n)

 integer :: inn, k, ndiv2, nm1, nm2, nnp1, np1  
 real(fpd) :: appfct, b, bfroot, bisq, con1, con2, der1p, der2p,dn
 real(fpd) :: pm1, pm2, prod, ratio, x, xi, aux, p

 real(fpd),parameter :: const = .148678816357, &
                        tol=1.0e-07, &
                        c1=.125 , &
                        c2=-.0807291666 , &
                        c3=.2460286458 , &
                        c4=-1.824438767

 if (n .eq. 1) then
         asc(1) = 0.5773502692
         wt(1) = 1.0
 else
         dn = dble(n)
         ndiv2 = n/2
         np1 = n + 1
         nnp1 = n*np1
         appfct = 1./sqrt((n + 0.5)**2 + const)
         con1 = 0.5*(bb - aa)
         con2 = 0.5*(bb + aa)

         do k = 1,ndiv2
            b = (k - 0.25)*pi
            bisq = 1.0/(b*b)
            bfroot = b*(1.0 + bisq*(c1 + bisq*(c2 + bisq*(c3 &
                     + c4*bisq))))
            x = cos(appfct*bfroot)
!           REPEAT
 100           pm2 = 1.0
               pm1 = x
               do inn = 2,n
                  p = ((2*inn - 1)*x*pm1 - (inn - 1)*pm2)/inn
                  pm2 = pm1
                  pm1 = p
               enddo
               pm1 = pm2
               aux = 1.0/(1.0 - x*x)
               der1p = dn*(pm1 - x*p)*aux
               der2p = (2.0*x*der1p - nnp1*p)*aux
               ratio = p/der1p
               xi = x - ratio*(1.0 + ratio*der2p/(2.0*der1p))
               if ((abs(xi - x) - tol) .gt. 0.0) then
                  x = xi
                  go to 100
               endif
!           UNTIL (ABS(XI - X) - TOL) <= 0.0
            asc(k) = -x
            wt(k) = 2.0*(1.0 - x*x)/(dn*pm1)**2
            asc(np1-k) = -asc(k)
            wt(np1-k) = wt(k)
         enddo

         if (mod(n,2) .ne. 0) then
            asc(ndiv2 + 1) = 0.0
            nm1 = n-1
            nm2 = n-2
            prod = dn
            do k = 1,nm2,2
               prod = float(nm1 - k)/float(n - k)*prod
            enddo
            wt(ndiv2 + 1) = 2.0/prod**2
         endif

         do k = 1,n
            asc(k) = con1*asc(k) + con2
            wt(k) = con1*wt(k)
         enddo

 endif
 return
 end subroutine gauss

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine init_tmatrix( axis_ratio, dcreal, dcimag, diam_sph, dwave, &
                          tmat_out, conk_in, sigma_in )
 implicit none

 real(fpd),intent(in) :: axis_ratio    ! particle axis ratio [-]
 real(fpd),intent(in) :: dcreal,dcimag ! real component of effective permittivity [-]
 real(fpd),intent(in) :: diam_sph      ! spherical diameter [cm]
 real(fpd),intent(in) :: dwave         ! wavelength [cm] 
 complex,intent(in)   :: tmat_out(20,40,40)   ! maybe dynamic allocation
 real(fpd),intent(in) :: conk_in, sigma_in ! 

! Objective: Initialize tmat() array.

!     DATA INITIALIZATION
! *** THERE ARE SOME CHANGES INTRODUCED HERE       *** RAVI: 12/30/91
!
 real(fpd) :: treal(40,40), timag(40,40)
 real(fpd) :: aoverb, deq, epi, epr,mode, sigma_temp
 integer :: i,j,nr2,mctr

 nrank = input_nrank  ! from module_tmatrix
 nm    = input_nm     ! from module_tmatrix
 conk  = conk_in      ! constant K
 lambda= dwave*0.01   ! wavelength [m]
 deq   = diam_sph     ! spherical/equivalent diameter [cm]
 aoverb= axis_ratio   ! axis ratio
 sigma_temp = sigma_in !
 epr = dcreal
 epi = dcimag
 
!
! read t-matrix
!
 tmat(1, 1:20,1:40,1:40 ) = tmat_out( 1:20,1:40,1:40 )

!
! this must be initialized...
!
 do i = 1,nm
    cmi(i) = float(i-1)
 enddo

 nr2    = 2*nrank
 nranki = nrank+1


 return
 end subroutine init_tmatrix

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine finddsd(ndset,n_htyp,nds)
 implicit none
!
! INTEGRATION OVER THE REQUIRED DSD IS PERFORMED BY THIS ROUTINE
!                                                 RAVI: 12/24/91
!
 integer,intent(in) :: ndset,n_htyp,nds

 integer :: nds1
 real(fpd) :: wtint, curdia, fndset, tempnd, fnds
 real(fpd) ::  tempnds

!
! common /dsd/ n0, d0, mgam, weight
! common /dsd1/ dmin, dmax, dstep

!
!
! *** INITIALIZE THE VARIABLE 'WTINT'
      wtint=0.
!
! *** FIRST CHECK IF NDSET IS EVEN, TO ASSIGN WEIGHTS 4(EVEN) OR 2(ODD).
      fndset = float(ndset)/2.0
      tempnd = abs(fndset-int(fndset))
      if(tempnd.lt.0.01) then
        wtint=4.0
      else
        wtint=2.0
      endif
!
      if(ndset.eq.1) then
        wtint=1
      endif
!
! *** CHECK IF NDS IS ODD WHICH IS REQUIRED FOR SIMPSON'S RULE TO WORK.
! *** IF NOT ADD A TRAPEZOIDAL STEP AT THE END
      fnds = float(nds)/2.0
      tempnds = abs(fnds-int(fnds))
      if(tempnds.gt.0.01) then
! ***   IF NDS IS ODD, THEN ADD WTINT=1 FOR THE LAST VALUE
        if(ndset.eq.nds) then
          wtint=1
        endif
      else
! ***   THIS IS FOR NDS EVEN ( ADD A TRAPEZOIDAL STEP )
        nds1 = nds-1
          if(ndset.eq.nds1) then
            wtint=2.50
          endif
          if(ndset.eq.nds) then
            wtint=1.50
          endif
      endif

!
! *** PRINT ERROR MESSAGE IF WTINT=0.                    RAVI: 12/24/91
      if ( wtint .eq. 0. ) then
        if(debug_mueller) print*,' ****   wtint = 0  *****'
      endif
!
!
! *** GAMMA DISTRIBUTION IS IMPLEMENTED HERE!
! *** N(D) = N0*D^(M)*EXP(-(3.67+M)*D/D0)
! *** N0 IS IN CM^(-1-M) M^(-3)
! *** D0 IS IN CMS. MGAM (IS M) IS UNITLESS
!
! *** COMPUTE THE CURRENT PARTICLE DIAMETER BEING PROCESSED!
      curdia = dmin + (ndset-1)*dstep
!
! *** THE VARIABLE 'WEIGHT(I)' INCLUDES THE N(D) FORM FOR THE GAMMA DSD.
! *** THE FACTOR H/3 (DSTEP/3.) FOR INTEGRATION USING SIMPSON'S RULE,
! *** IS ALSO INCLUDED IN THE VARIABLE 'WEIGHT(I)'.
!
! *** WEIGHTS TO ALL N(D),(N0,D0,M-TRIPLETS), IS COMPUTED FOR THE
! *** CURRENT PARTICLE DIAMETER BEING PROCESSED!
!
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      Gamma DSD
       weight=exp((-1.0)*(3.67+mgam(n_htyp))*curdia/d0(n_htyp))
       weight=wtint*n0(n_htyp)*(curdia**mgam(n_htyp))*weight*dstep/3.0

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!        WEIGHT=EXP((-1.0)*D0(II)*CURDIA)
!        WEIGHT=WTINT*N0(II)*WEIGHT*DSTEP/3.0
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  return
 end subroutine finddsd

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine thetaint
 implicit none

! OUTER LOOP:  INTEGRATE OVER THETA

 integer :: ith,is,ir,ic
 real(fpd) :: prob


!
! theta loop ntheta = 20
!
 do ith = 1,ntheta    
!   CALCULATE THETA OF CURRENT PARTICLE ORIENTATION
    angint = 2.0*atan(sqrt(1.0-tasc(ith)**2)/ (1.0 + tasc(ith)))
    oldtht = angint
    if(debug_mueller) write(io_check_out,901)
    if(debug_mueller) write(io_check_out,902)

!        initialize variables
    call zero_output(nang)

!        do inner phi loop (integrating phi in mueller matrix)
    call phiint(prob,ith)

!        multiply by the theta weights determined by type of orienting angle distributions
    do is = 1,nang
       do ir = 1,4
          do ic = 1,4
             ss(is,ir,ic) = ss(is,ir,ic) + twt(ith)*s(is,ir,ic)
          enddo
       enddo
    enddo
!
!  do forward amplitudes scattering 2X2 matrix (twt is ith gaussian weight of theta angle
!
    do ir = 1,2
       do ic = 1,2
          fscatt(ir,ic) = fscatt(ir,ic) + twt(ith)*fscat(ir,ic)
       enddo
    enddo
    schpt = schpt + prob*twt(ith)*schp    !scattering cross section (H)
    exhpt = exhpt + prob*twt(ith)*exhp    !extinction cross section (H)
    scvpt = scvpt + prob*twt(ith)*scvp    !scattering cross section (V)
    exvpt = exvpt + prob*twt(ith)*exvp    !extinction cross section (V)

 enddo

 901  format (/,1x,'theta',5x,'phi',4x,'original',2x,'original',2x, &
             'rotated',2x,'rotated')
 902  format ('counter',2x,'counter',3x,'theta',6x,'phi',6x, &
              'theta',5x,'phi',/)

 return
 end subroutine thetaint

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine phiint(prob,ith)
 implicit none
!     INNER LOOP:  INTEGRATE OVER PHI

 integer,intent(in) :: ith
 real(fpd),intent(out) :: prob

 integer :: iphi

 do iphi = 1,nphi
!  CALCULATE PHI OF CURRENT PARTICLE ORIENTATION; RECOVER
!  UNROTATED THETA OF PARTICLE ORIENTATION
    if (phiflag .eq. 0 .or. distyp .eq. 3)then
        anginp = pasc(iphi)
    else
        anginp = phiang*degrad
    endif

    oldphi = anginp
    angint = oldtht
!   CALCULATE SCATTERED FIELD FOR CURRENT THETA, PHI ORIENTATION
!   OF PARTICLE
    call scfld(ith,iphi)
!
    call getprob(iphi,prob)

!   ACCUMULATE TOTAL SCATTERING IN MUELLER-MATRIX ACCUMULATOR
!   itegration of phi is done in mmcalc

    call mmcalc(iphi,schp, scvp, exhp, exvp, prob)
 enddo

 return
 end subroutine phiint

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine scfld (tcount, pcount)
 implicit none

 integer,intent(in) :: tcount, pcount

 integer :: i, iz, im, ifct
 real(fpd) ::  quanm
 real(fpd) :: tempth, tempph, wvnmbr

!     ADJUST ORIENTATION VECTOR FOR +Z-INCIDENT WAVE
 tempth = angint
 tempph = anginp
 call rotate(elvang,angint,anginp)
 if(debug_mueller) write(io_check_out,900) tcount, pcount, tempth, tempph, angint, anginp

!     ZERO SCATTERED-FIELD ARRAY
 do i=1,nuang
    do iz = 1,2
       acans(i,1,iz) = 0.0
       acans(i,2,iz) = 0.0
    enddo
 enddo
 sch = 0.0

!     SET POLARIZATION FOR FIRST LOOP THROUGH
 utheta = 0.0
 uphi = 1.0
 ip = 1

!     GENERATE TRANSFORMATION ANGLES BETWEEN LAB AND BODY FRAMES
 call genang(nuang,tang,pang,cosb,sinb)

!     REPEAT
!        ONE CALL TO ADDPRC FOR EACH MODE:

 300     do im = 1,nm
            cmv = cmi(im)
            kmv = int(cmv)  !toshii add int for type change
            cm2 = cmv**2
            prodm = 1.0
            if (kmv .ne. 0) then
               quanm = cmv
               do ifct = 1,kmv
                  quanm = quanm + 1.0
                  prodm = quanm*prodm/2.0
               enddo
            endif
            twm = 2.0*cmv

            call addprc
         enddo

         if (ip .eq. 1) then
!           SET UP FOR NEXT RUN, ORTHOGONAL POLARIZATION
            scv = 0.0

            utheta = 1.0
            uphi = 0.0
            ip = 2
            go to 300
         endif
!     UNTIL IP = 2

!     CALCULATE WAVENUMBER AND NORMALIZING FACTOR
      wvnmbr = 2.0*pi/lambda
      do i=1,nuang
         do iz = 1,2
            acans(i,1,iz) = (4.0/wvnmbr)*acans(i,1,iz)
            acans(i,2,iz) = (4.0/wvnmbr)*acans(i,2,iz)
         enddo
      enddo

!JH   PRINT SCATTERING 2X2 MATRICIES HERE

 900  format (2x,i2,7x,i2,6x,f6.2,4x,f6.2,3x,f6.2,3x,f6.2)

 return
 end subroutine scfld

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine rotate (elvang, thetap, phip)
 implicit none
!     IN THIS PROGRAM, THE INCIDENT WAVE TRAVELS IN THE PLUS-Z DIRECTION.
!     THIS ROUTINE, ASSUMING A MINUS-X-DIRECTED INCIDENT WAVE AND AN
!     ELEVATION ANGLE, ROTATES THE PARTICLE ORIENTATION VECTOR SO THAT
!     THE PLUS-Z-DIRECTED WAVE "SEES" THE SAME ORIENTATION AS THE USER-
!     SPECIFIED WAVE.

 real(fpd),intent(in) :: elvang
 real(fpd),intent(inout) :: thetap, phip

 real(fpd) :: magxy
 real(fpd) :: elvrad, therad, phirad, x, y, z, xprime, yprime, zprime


!     CONVERT TO RADIANS
!     ELVRAD = DEGRAD*ELVANG
!     THERAD = DEGRAD*THETAP
!     PHIRAD = DEGRAD*PHIP
!
!     all angles are now in radians: oct 25, 1993 j. hubbert
!
 elvrad = elvang
 therad = thetap
 phirad = phip

!     CALCULATE CARTESIAN COMPONENTS OF UNROTATED ORIENTATION VECTOR
 x = sin(therad)*cos(phirad)
 y = sin(therad)*sin(phirad)
 z = cos(therad)

!     MULTIPLY ORIENTATION VECTOR BY ROTATION MATRIX
 xprime = sin(elvrad)*x + cos(elvrad)*z
 yprime = y
 zprime = -cos(elvrad)*x + sin(elvrad)*z
 if (abs(xprime) .lt. 1.0e-03) xprime = 0.0
 if (abs(yprime) .lt. 1.0e-03) yprime = 0.0
 if (abs(zprime) .lt. 1.0e-03) zprime = 0.0

!     CONVERT ROTATED CARTESIAN COMPONENTS INTO ANGLES THETA AND PHI
!     THETA:
 magxy = sqrt(xprime*xprime + yprime*yprime)
 if (abs(zprime) .gt. 0.0) then
     thetap = atan(magxy/zprime)
     if (thetap .lt. 0.0) thetap = thetap + pi
     if ((thetap .eq. 0.0) .and. (zprime .lt. 0.0)) thetap = pi
 elseif (abs(zprime) .eq. 0.0) then
     if (abs(magxy) .gt. 0.0) then
         thetap = pi/2.0
     else
         thetap = 0.0
     endif
 endif

!     PHI:
 if (yprime .eq. 0.0) then
     if (xprime .lt. 0.0) then
         phip = pi
     else
         phip = 0.0
     endif
 else
     phip = 2.0*atan(yprime/(xprime + magxy))
     if (phip .lt. 0.0) phip = phip + 2.0*pi
 endif

 return
 end  subroutine rotate

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine genang (n,tsp,psp,cosbta,sinbta)
 implicit none
!     A ROUTINE TO TRANSFORM SCATTERING ANGLES FROM THE LAB FRAME TO THE
!     BODY FRAME  AND TO GENERATE THE TRANSFORMATION MATRIX FOR THE
!     POLARIZATION VECTORS

 integer,intent(in) :: n
 real(fpd),dimension(n) :: tsp,psp,cosbta,sinbta

 integer :: i

 real(fpd) :: chi, alpha, sinchi, coschi, sinalf, cosalf
 real(fpd) :: thetas, sints, costs, cosps, xprime, yprime, &
              zprime, magxy
 real(fpd) :: phisp, sinpsp, cospsp

!     CALCULATE CHI AND ALPHA (WANG'S THESIS, EQNS. 72 AND 73)
 chi = anginp - pi
 alpha = -angint

!     CALCULATE SINES AND COSINES OF CHI AND ALPHA
 sinchi = sin(chi)
 coschi = cos(chi)
 sinalf = sin(alpha)
 cosalf = cos(alpha)

!        write(80,1000)chi,alpha
! 1000   format(2(e15.6))

!     CALCULATE THETA (SCATTERED) PRIME -- "TSP" -- AND PHI
!     (SCATTERED) PRIME -- "PSP" -- FOR ALL LAB-FRAME THETA
!     (SCATTERED) -- "THETAS" -- ANGLES
 do i = 1,n
!    THETA (SCATTERED) AND ITS SINE AND COSINE:
    if(i.le.(n/2)) then
!U            QMUE(I)=QMU(I)
       qmue(i)=qmu(nmu+1-i)
    else
!U             QMUE(I)=-QMU(I-NMU)
       qmue(i)=-qmu(n+1-i)
    endif

    thetas = acos(qmue(i))
    sints = sin(thetas)
    costs = cos(thetas)

!        PHI (SCATTERED) AND ITS COSINE (ACTUALLY ONLY NEED COSINE OF PHI
!        (SCATTERED)).  THIS WAS ADDED TO IMPOSE THE REQUIREMENT THAT
!        PHI (SCATTERED) BE EQUAL TO ZERO EXCEPT AT THETA (SCATTERED)
!        EQUAL TO PI, WHERE PHI (SCATTERED) IS REQUIRED TO EQUAL PI.
!        NOTE THAT THIS CHANGE MAKES NO DIFFERENCE IN THIS LOOP, SINCE
!        THE DIFFERENT COS(PS) IS MULTIPLIED BY SIN(PI), WHICH IS ZERO
!        ANYWAY.  REAL DIFFERENCE COMES IN THE SINE BETA/ COSINE BETA
!        LOOP AT END OF THIS SUBROUTINE.
    if (i .lt. n) then
!           ALL OTHER ANGLES OTHER THAN BACKSCATTER:
        cosps = 1.00
    else if (i .eq. n) then
!           BACKSCATTER:  PHI (SCATTERED) = PI
       cosps = -1.00
    endif

!        X', Y' AND Z' (WANG, EQNS. 84 MODIFIED FOR PHI (SCATTERED))
    xprime = sints*cosalf*coschi*cosps - costs*sinalf
    yprime = -sints*sinchi*cosps
    zprime = sints*sinalf*coschi*cosps + costs*cosalf

!        THETA (SCATTERED) PRIME (WANG, EQN. 85A)
    magxy = sqrt(xprime*xprime + yprime*yprime)
    if (abs(zprime) .gt. 0.0) then
        tsp(i) = atan(magxy/zprime)
        if (tsp(i) .lt. 0.0) tsp(i) = tsp(i) + pi
        if ((tsp(i) .eq. 0.0) .and. (zprime .lt. 0.0)) tsp(i) = pi
    elseif (abs(zprime) .eq. 0.0) then
        if (abs(magxy) .gt. 0.0) then
            tsp(i) = pi/2.0
        else
            tsp(i) = 0.0
        endif
    endif

!        PHI (SCATTERED) PRIME (WANG, EQN. 85B -- NOTE ERROR IN
!        EQUATION:  SHOULD BE ATAN(Y'/X'))
    if (abs(xprime) .gt. 0.0) then
        if (abs(xprime + magxy) .lt. 1.0e-20) then
            psp(i) = pi
        else
            psp(i) = 2.0*atan(yprime/(xprime + magxy))
        endif
        if (psp(i) .lt. 0.0) psp(i) = psp(i) + 2.0*pi
    elseif (abs(xprime) .eq. 0.0) then
        if (yprime .lt. 0.0) then
            psp(i) = 1.5*pi
        elseif (yprime .eq. 0.0) then
            psp(i) = 0.0
        elseif (yprime .gt. 0.0) then
            psp(i) = 0.5*pi
        endif
    endif

 enddo

!     CALCULATE THE ELEMENTS OF THE TRANSFORMATION
!     MATRIX BETWEEN LAB AND BODY FRAMES FOR THE
!     DIRECTION OF POLARIZATION
 do i = 1,n
!        PHI (SCATTERED) PRIME AND SINE AND COSINE
    phisp = psp(i)
    sinpsp = sin(phisp)
    cospsp = cos(phisp)

!        THETA (SCATTERED) AND SINE AND COSINE; COSINE OF PHI (SCATTERED)
    thetas = acos(qmue(i))
    sints = sin(thetas)
    costs = cos(thetas)
!        SEE COMMENTS IN "DO 100" LOOP ABOUT PHI (SCATTERED)
    if (i .lt. n) then
!           ALL OTHER ANGLES OTHER THAN BACKSCATTER:
       cosps = 1.00
    elseif (i .eq. n) then
!           BACKSCATTER:  PHI (SCATTERED) = PI
       cosps = -1.00
    endif

!        COSINE BETA AND SINE BETA (WANG, EQNS. 93 -- NOTE COUPLE OF
!        ERRORS IN 93B)
    cosbta(i) = -sinpsp*cosalf*sinchi*cosps + cospsp*coschi*cosps
    sinbta(i) = sinpsp*cosalf*coschi*costs*cosps &
               + sinalf*sints*sinpsp + sinchi*costs*cospsp*cosps
 enddo

 return
 end subroutine genang

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine addprc
 implicit none

!     A ROUTINE TO OBTAIN THE SCATTERED FIELD COEFFICIENTS AND CALCULATE
!     THE DIFFERENTIAL SCATTERING CROSS SECTION.
 real(fpd) :: cmxnrm(40)  !this is used in one subroutine....
 complex :: ab1(40),ab2(40),ci,c1,c2,s1,s2,cim
 complex :: fg1(40),fg2(40)

 integer :: i,ibfct,iefct,irow,iu,j,kup,lfct,n,n1,np,nr2
 
 real(fpd) :: aa,bb,cc,ckrow,cn,cosphi,dd,em,fctki,fprod,p1,p2,phi,sinphi,uu1,uu2

      nr2 = 2*nrank
      ci = (0.0,1.0)
!     TRANSFORM THE PROBLEM FROM THE LAB FRAME TO THE BODY FRAME
!     TANG(1) = INCIDENT ANGLE IN THE BODY SYSTEM
!     GENERATE THE LEGENDRE FUNCTIONS FOR THE INCIDENT ANGLE.
      theta=tang(1)
      costh = cos(theta)
      sinth = sin(theta)
      call genlgp

!     TRANSFORM THE INCIDENT FIELD POLARIZATION VECTOR FROM THE LAB FRAME TO
!     THE BODY FRAME:
!     UU1 = COMPONENT IN THETA DIRECTION; UU2 = COMPONENT IN PHI DIRECTION
      uu1 = cosb(1)*utheta + sinb(1)*uphi
      uu2 = -sinb(1)*utheta + cosb(1)*uphi
!        write(80,1000)cosb(1),sinb(1)
! 1000   format(2(e15.6))

!     GENERATE THE INCIDENT FIELD COEFFICIENTS -- AB1 = THETA POLARIZATION
!     AND AB2 = PHI POLARIZATION.
      cn = 0.0
      do n = 1,nrank
         np = n + nrank
         cn = cn + 1.0
         n1 = n + 1
         c1 = ci**n
         c2 = ci**n1
         p1 = cn*costh*pnmllg(n1) - (cn+cmv)*pnmllg(n)
         p2 = cmv*pnmllg(n1)
         ab1(n) = -c1*p2*uu1
         ab1(np) = c2*p1*uu1
         ab2(n) = c1*p1*uu2
         ab2(np) = -c2*p2*uu2
      enddo

!     THE SCATTERED FIELD COEFFICIENTS = THE TRANSITION MATRIX TIMES THE
!     INCIDENT FIELD COEFFICIENTS. (f = Ta, *** refer routine INIT ***)
!                                                        RAVI: 1/9/92
      do i = 1,nr2
         s1 = 0.0
         s2 = 0.0
         do j = 1,nr2
            s1 = s1 + tmat(ndset_global,kmv+1,j,i)*ab1(j)
            s2 = s2 + tmat(ndset_global,kmv+1,j,i)*ab2(j)
         enddo
         fg1(i) = s1
         fg2(i) = s2
      enddo

!     CALCULATE THE NORMALIZATION FACTOR
      if (kmv .eq. 0) then
         em = 1.0
      else
         em = 2.0
      endif
      do 450 irow = 1,nrank
      ckrow = irow
      if(kmv.gt.0) go to 426
      fctki = 1.0
      go to 440

  426 if(irow.ge.kmv) go to 430
      cmxnrm(irow) = 1.0
      go to 450
  430 ibfct = irow-kmv+1
      iefct = irow+kmv
      fprod = ibfct
      fctki = 1.0
      do 432 lfct = ibfct,iefct
      fctki = fctki*fprod
      fprod = fprod+1.0
  432 continue
  440 cmxnrm(irow) = 4.0*ckrow*(ckrow+1.0)*fctki/(em*(2.0*ckrow+1.0))
!      WRITE(6,"(' CMXNRM', I5,D15.7)") IROW,CMXNRM(IROW)
  450 continue

!     EVALUATE THE SCATTERED FIELD AT EACH SCATTERING ANGLE.
      do iu = 1,nuang
!        GENERATE THE LEGENDRE MULTIPLIERS.
         theta = tang(iu)
         sinth = sin(theta)
         costh = cos(theta)
         if (tang(iu) .eq. pi) then
            theta = 0.0
            costh = -1.0
         endif
         call genlgp
         phi = cmv*pang(iu)
         sinphi = sin(phi)
         cosphi = cos(phi)
         s1 = 0.0
         s2 = 0.0
         cn = 0.0
         do n = 1,nrank
            np = n + nrank
            if ((iu .eq. 1) .and. (ip .eq. 1)) then
               sch = sch + (cabs(fg1(n))**2 + cabs(fg1(np))**2 &
                    + cabs(fg2(n))**2 + cabs(fg2(np))**2)/cmxnrm(n)
            endif

            if ((iu .eq. 1) .and. (ip .eq. 2)) then
               scv = scv + (cabs(fg1(n))**2 + cabs(fg1(np))**2 &
                    + cabs(fg2(n))**2 + cabs(fg2(np))**2)/cmxnrm(n)
            endif
            n1 = n + 1
            cn = cn + 1.0
            p1 = cn*costh*pnmllg(n1) - (cn + cmv)*pnmllg(n)
            p2 = cmv*pnmllg(n1)
            aa = sinphi*p1
            bb = cosphi*p1
            cc = sinphi*p2
            dd = cosphi*p2
            cim = (-ci)**n1

!           SOLVE FOR THE THETA-POLARIZED SCATTERED FIELD
!           IN THE BODY FRAME
            s1 = s1 + cim*(fg1(n)*dd + ci*fg1(np)*bb - fg2(n)*cc &
                 - ci*fg2(np)*aa)/cmxnrm(n)

!           SOLVE FOR THE PHI-POLARIZED SCATTERING FIELD
!           IN THE BODY FRAME
            s2 = s2 - cim*(fg1(n)*aa + ci*fg1(np)*cc + fg2(n)*bb &
                 + ci*fg2(np)*dd)/cmxnrm(n)
        enddo

         acans(iu,1,ip) = acans(iu,1,ip) + s1
         acans(iu,2,ip) = acans(iu,2,ip) + s2
      enddo

!     TRANSFORM THE SCATTERED FIELD (ACANS) BACK TO THE LAB FRAME
 if (kmv .eq. (nm-1)) then
     do kup = 1,nuang
            s1 = acans(kup,1,ip)
            s2 = acans(kup,2,ip)
            acans(kup,1,ip) = cosb(kup)*s1 - sinb(kup)*s2
            acans(kup,2,ip) = sinb(kup)*s1 + cosb(kup)*s2
     enddo
     if (ip .eq. 1) exh = aimag(acans(1,2,1))
     if (ip .eq. 2) exv = aimag(acans(1,1,2))
 endif

 return
 end subroutine addprc

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine genlgp
 implicit none

 real(fpd) :: cnm, cnmm, cnmul, dtwm, pla, plb, plc
 integer :: ibeg, ilg, ilgr


   dtwm = twm+1.0
   if(theta.ne.0.0) go to 16
   if(kmv.eq.1) go to 12

   do ilg = 1,nranki
      pnmllg(ilg) = 0.0
   enddo

      go to 88
   12 pnmllg(1) = 0.0
      pnmllg(2) = 1.0
      pla = 1.0
      go to 48
   16 if(kmv.gt.0) go to 40
!     THE SPECIAL CASE WHEN M = 0.
      pla = 1.0/sinth
      plb = costh*pla
      pnmllg(1) = pla
      pnmllg(2) = plb
      ibeg = 3
      go to 60
!     GENERAL CASE FOR M NOT EQUAL TO 0.
   40 do 44 ilg = 1,kmv
      pnmllg(ilg) = 0.0
   44 continue
      pla = prodm*sinth**(kmv-1)
      pnmllg(kmv+1) = pla
   48 plb = dtwm*costh*pla
      pnmllg(kmv+2) = plb
      ibeg = kmv+3
!     DO RECURSION FORMULA FOR ALL REMAINING LEGENDRE POLYNOMIALS.
   60 cnmul = ibeg+ibeg-3  !toshii is this OK ?? integer --> double
      cnm = 2.0
      cnmm = dtwm

   do ilgr = ibeg,nranki
      plc = (cnmul*costh*plb-cnmm*pla)/cnm
      pnmllg(ilgr) = plc
      pla = plb
      plb = plc
      cnmul = cnmul+2.0
      cnm = cnm+1.0
      cnmm = cnmm+1.0
   enddo
   88 return

 end subroutine genlgp

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine getprob(iphi,prob)
 implicit none

 integer,intent(in) :: iphi
 real(fpd),intent(out) :: prob


 if(iphi==1) then !toshii just in case.
 endif

!           ACCUMULATE PHI-AVERAGED CONTRIBUTION TO CURRENT THETA TERM:
!           GET ORIENTATION PROBABILITY DENSITY FUNCTION, TAKING CARE
!           TO USE THE ORIGINAL (NOT ROTATED) THETA AND PHI.
 if (distyp .eq. 0) then
!    P(THETA) = 0.5; P(PHI) = 1.0/2*PI (BOTH RANDOM)
     prob = 0.5/(2*pi)
 elseif (distyp .eq. 1) then
!    SIMPLE HARMONIC OSCILLATOR MODEL
     prob = shop(oldtht, mean, thetam, oldphi, pnorm)
!    ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE

     if (oldtht .lt. tcu) then
         prob = prob + shop(-oldtht,mean,thetam,(oldphi+pi),pnorm)
     elseif (oldtht .gt. tcl) then
         prob = prob + shop((2.*pi -oldtht),mean,thetam,(oldphi+pi),pnorm)
     endif
 elseif (distyp .eq. 2) then
!    2-SIGMA-TRUNCATED GAUSSIAN MODEL

     prob = gp(oldtht, mean, sigma, oldphi, pnorm)

!              ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE
!jh            IF (OLDTHT .LT. TCU) THEN
!jh               PROB = PROB +
!jh  >            GP(-OLDTHT,MEAN,SIGMA,(OLDPHI+PI),PNORM)
!jh            ELSE IF (OLDTHT .GT. TCL) THEN
!jh               PROB = PROB +
!jh  >            GP((2*PI -OLDTHT),MEAN,SIGMA,(OLDPHI+PI),PNORM)
!jh            ENDIF
 elseif (distyp .eq. 3) then
!              LANGEVIN UNIMODAL DISTRIBUTION
     prob = fp(kappa, oldtht, oldphi, theta0, phi0, pnorm)
 endif

! 903  format (1x,'(orientation probability density function = ',f8.5, &
!              ';',/,2x,'normalization factor = ',d12.5,')',/)

 return
 end subroutine getprob

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 function shop(theta, thbar, thmax, phi, norm)
 implicit none
!     SIMPLE HARMONIC OSCILLATOR ORIENTATION PROBABILITY DENSITY FUNCTION
!     THETA, THBAR, THETA-MAX (THMAX) AND PHI ARE ASSUMED GIVEN IN DEGREES
!     "NORM" IS CALCULATED IN SUBROUTINE "USERIN" AND
!     NORMALIZES THE FUNCTION TO UNITY OVER THE RANGE
!     OF THETA AND PHI SPECIFIED.
!     NOTE THAT P(THETA,PHI) IS ASSUMED SEPARABLE IN THETA AND PHI.

 real(fpd) :: shop
 real(fpd),intent(in) :: theta, thbar, thmax, phi, norm
 
 real(fpd) :: const, ratio
 real(fpd) :: ptheta, pphi

 const = 1.0/(2.0*pi*pi*thmax)
 ratio = (theta - thbar)/thmax
 ptheta = const/sqrt(1.0 - ratio*ratio)
 if(phiflag .eq. 0)then
    pphi = 1.0/(2*pi)
 else
    pphi = 1.
 endif
 shop = ptheta*pphi/norm

 return
 end function shop

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 function gp (theta, thbar, sigma, phi, norm)
 implicit none
!     GAUSSIAN ORIENTATION PROBABILITY DENSITY FUNCTION
!     THETA, THBAR, SIGMA AND PHI ARE ASSUMED GIVEN IN DEGREES
!     "NORM" IS CALCULATED IN SUBROUTINE "USERIN" AND
!     NORMALIZES THE FUNCTION TO UNITY OVER THE RANGE
!     OF THETA AND PHI SPECIFIED.
!     NOTE THAT P(THETA,PHI) IS ASSUMED SEPARABLE IN THETA AND PHI.

 real(fpd) :: gp
 real(fpd),intent(in) :: theta, thbar, sigma, phi, norm

 real(fpd) :: const, ratio, ptheta, pphi

 if( sigma == 0. ) print*,'MSG gp: warning input sigma is zero'
 const = 1.0/(sqrt(2.0*pi)*sigma)
 ratio = (theta - thbar)/(sqrt(2.0)*sigma)
 ptheta = const*exp(-ratio*ratio)
 if(phiflag .eq. 0)then
    pphi = 1.0/(2*pi)
 else
    pphi = 1.
 endif

 gp = ptheta*pphi/norm

 return
 end function gp

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 function fp(kappa, theta, phi, theta0, phi0, norm)
 implicit none
!     LANGEVIN ORIENTATION PROBABILITY DISTRIBUTION
!     THETA, PHI, THETA0, AND PHI0 ARE ASSUMED GIVEN IN DEGREES
!     "NORM" IS CALCULATED IN SUBROUTINE "USERIN" AND
!     NORMALIZES THE FUNCTION TO UNITY.

 real(fpd) :: fp
 real(fpd),intent(in) :: kappa, theta, phi, theta0, phi0, norm

 real(fpd) :: c, term1, term2

 c=kappa/(2*pi*(1.0- exp(-2*kappa)))
 term1 = cos(theta0)*cos(theta)
 term2 = sin(theta0)*sin(theta)* cos((phi - phi0))
 fp = c*exp(kappa*(term1 + term2 -1.))/(norm)

 return
 end function fp

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine mmcalc(iphi,schp,scvp,exhp, exvp,prob)
 implicit none

!     ACCUMULATE CURRENT CONTRIBUTION TO MUELLER-MATRIX ACCUMULATOR

 integer,intent(in)    :: iphi
 real(fpd),intent(in)  :: prob
 real(fpd),intent(out) :: schp, scvp, exhp, exvp


 integer j, row, col

 complex  :: f(2,2)

 real(fpd) :: mmatrx(4,4)

!
! loop over entire scattering angle loop (1: forward scattering)
! nang = nmu * 2 = 8  (1 is 0deg, each deg increment 45 deg)
! backscatteirng angle is nmu 5 -->  (5-1)*45 = 180

 do j = 1,nang

!        PREPARE TO PASS SCATTERING-AMPLITUDE MATRIX f TO "MMAT"
    f(1,1) = acans(j,1,1)
    f(1,2) = acans(j,1,2)
    f(2,1) = acans(j,2,1)
    f(2,2) = acans(j,2,2)

!        CALL "MMAT"; GET MUELLER MATRIX BACK
!
!     CALCULATE (I,Q,U,V) MUELLER MATRIX FROM SCATTERING-AMPLITUDE MATRIX F.
!     THE SCATTERING-AMPLITUDE MATRIX IS DEFINED IN TERMS OF THE E-FIELD
!     COMPONENTS PARALLEL AND PERPENDICULAR TO THE PLANE OF SCATTERING.
!
    call mmat(f,mmatrx)

!        WEIGHT SUMMATION TERMS APPROPRIATELY -- TOTAL OF 16
!        (NOTE:  NORMALIZATION OF PHI INTEGRATION IS TAKEN CARE
!        OF BY THE VARIABLE "PNORM" IN THE THETA INTEGRATION.)
    do row = 1,4
       do col = 1,4
          s(j,row,col) = s(j,row,col) + prob*pwt(iphi)*mmatrx(row,col) ! pwt(iphi) is phi weight... ramdom --> does not matter
       enddo
    enddo

 enddo


!****call to create rayleigh extinction matrix(jb)
!     call emat(fr,rexm,lambda)

! integrate the forward scattering amplitudes
! minus sign converts from tmatrix notation to an optical FSA one (jb)
! SHH
  f(1,1) = acans(1,2,1)
! SHV
  f(1,2) = acans(1,2,2)
! SVH
  f(2,1) = acans(1,1,1)
! SVV
  f(2,2) = acans(1,1,2)

  do row = 1,2
     do col = 1,2
        fscat(row,col) = fscat(row,col) + prob*pwt(iphi)*f(row,col)
     enddo
  enddo

  schp = schp + pwt(iphi)*sch
  exhp = exhp + pwt(iphi)*exh
  scvp = scvp + pwt(iphi)*scv
  exvp = exvp + pwt(iphi)*exv


 return
 end subroutine mmcalc

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine mmat (f, m )
 implicit none

!     CALCULATE (I,Q,U,V) MUELLER MATRIX FROM SCATTERING-AMPLITUDE MATRIX F.
!     THE SCATTERING-AMPLITUDE MATRIX IS DEFINED IN TERMS OF THE E-FIELD
!     COMPONENTS PARALLEL AND PERPENDICULAR TO THE PLANE OF SCATTERING.

!     REFERENCE:
!        LIGHT SCATTERING BY SMALL PARTICLES
!        H. C. VAN DE HULST
!        DOVER PUBLICATIONS, INC. (NEW YORK), 1981;
!        PP. 34, 43-44, 46.
!        (NOTE VAN DE HULST'S NOTATION FOR THE M AND A MATRICES; THE
!        COMPONENTS OF THE F MATRIX WILL BE COPIED INTO THE APPROPRIATE
!        A-MATRIX COMPONENTS FIRST TO AVOID CONFUSION.)
!     W. M. ADAMS     CSU -- FT. COLLINS, CO     31 JULY 1988
!     INPUT TO SUBROUTINE:
!        F  --  COMPLEX 2-BY-2 SCATTERING-AMPLITUDE MATRIX
!        LAMBDA  --  WAVELENGTH IN MM
!     RETURNED BY SUBROUTINE:
!        M  --  REAL 4-BY-4 MUELLER MATRIX
!     LOCAL VARIABLES (SEE REFERENCE, PAGE 44 FOR NOTATION
!     USED IN LOCAL VARIABLES):
!        A1   --  ELEMENT OF TRANSFORMATION MATRIX A = F(2,1)
!        A2   --  ELEMENT OF TRANSFORMATION MATRIX A = F(1,2)
!        A3   --  ELEMENT OF TRANSFORMATION MATRIX A = -F(1,1)
!        A4   --  ELEMENT OF TRANSFORMATION MATRIX A = -F(2,2)
!        M1   --  SQUARED MAGNITUDE OF A1
!        M2   --  SQUARED MAGNITUDE OF A2
!        M3   --  SQUARED MAGNITUDE OF A3
!        M4   --  SQUARED MAGNITUDE OF A4
!        S21  --  (A2*CONJ(A1) + A1*CONJ(A2))/2
!        S23  --  (A2*CONJ(A3) + A3*CONJ(A2))/2
!        S24  --  (A2*CONJ(A4) + A4*CONJ(A2))/2
!        S31  --  (A3*CONJ(A1) + A1*CONJ(A3))/2
!        S34  --  (A3*CONJ(A4) + A4*CONJ(A3))/2
!        S41  --  (A4*CONJ(A1) + A1*CONJ(A4))/2
!        D21  --  (I/2)*(A2*CONJ(A1) - A1*CONJ(A2))
!        D23  --  (I/2)*(A2*CONJ(A3) - A3*CONJ(A2))
!        D24  --  (I/2)*(A2*CONJ(A4) - A4*CONJ(A2))
!        D31  --  (I/2)*(A3*CONJ(A1) - A1*CONJ(A3))
!        D34  --  (I/2)*(A3*CONJ(A4) - A4*CONJ(A3))
!        D41  --  (I/2)*(A4*CONJ(A1) - A1*CONJ(A4))
!     VERY IMPORTANT NOTES:
!     ->  VAN DE HULST USES L (PARALLEL) AND R (PERPENDICULAR) NOTATION:
!         VECTOR R CROSS VECTOR L = DIRECTION OF PROPAGATION.
!     ->  THIS PROGRAM USES UNIT VECTORS A(THETA) AND A(PHI); A(THETA)
!         CROSS A(PHI) = DIRECTION OF PROPAGATION.
!     ->  VECTOR R AND A(PHI) POINT IN OPPOSITE DIRECTIONS, NECESSITATING* *
!         A SIGN CHANGE (WHEN CONVERTING TO VAN DE HULST'S NOTATION) EACH
!         TIME THE PERPENDICULAR POLARIZATION IS INVOLVED.  SEE COPY INTO
!         VAN DE HULST'S "A" MATRIX, BELOW.
!****************************************************************
!****************************************************************
!***************************************************************
! changed by john hubbert sept. 22, 1993
!**************************************************************
!  This program, as written by Adams, uses van de Hulst's
!  coordinate system convention for mueller matrix constuction.
!  Incidence is along the positive z direction with x cross y
!  giving the z+ direction. However, this program uses a spherical
!  coordinate system for the scattering computations. Thus:
!
!      phi_inc = phi_scat = y (rotation from x in x-y plane)
!      theta_inc = x ; theta_scat = -x
!      Adams used van de Hulst coordinates:
!
!      parl_inc = theta_inc ; parl_scat = theta_scat
!      perp_inc = -phi_inc  ; perp_scat = -phi_scat
!
!   In order to apply the results of this program to polarimetric
!   studies, the following is introduced so that a consistantcy
!   of definitions is maintained. Note that before this modification,
!   V x H did not give propagation direction.
!
!     vertical_inc = theta_inc ; vertical_scat = -theta_scat
!     horiz_inc    = phi_inc ;      horiz_scat = -phi_scat
!
!    Thus since vertical_scat = -theta_scat and
!    horiz_scat = -phi_scat , a minus sign will
!    be introduced where the theta_scat and horiz_scat = -phi_scat
!    term appears. The scattering
!    covariances will thus be given in an H/V basis for the optic
!    or FSA convention.
!    To switch to the radar or BSA convention, negative signs
!    will need to be introduced on the S_hh and S_hv back scatter
!    amplitudes.
!
!****************************************************************
!***************************************************************
!***************************************************************

 complex :: f(2,2)

 real(fpd) ::  m(4,4)

! real(fpd) :: lambda

 complex :: a1, a2, a3, a4, shh,shv,svh,svv
 real(fpd) :: m1, m2, m3, m4


!     SET TRANSFORMATION MATRIX VALUES

!     F(1,1):  PARALLEL COMPONENT, PERPENDICULAR INCIDENT POLARIZATION
!                (VAN DE HULST'S A3 -- ONE SIGN CHANGE)
!jh      A3 = -F(1,1)

!     F(1,2):  PARALLEL COMPONENT, PARALLEL INCIDENT POLARIZATION
!                (VAN DE HULST'S A2 -- NO SIGN CHANGES)
!jh      A2 = F(1,2)

!     F(2,1):  PERPENDICULAR COMPONENT, PERPENDICULAR INCIDENT POLARIZATION
!                (VAN DE HULST'S A1 -- TWO SIGN CHANGES WHICH CANCEL)
!jh      A1 = F(2,1)

!     F(2,2):  PERPENDICULAR COMPONENT, PARALLEL INCIDENT POLARIZATION
!                (VAN DE HULST'S A4 -- ONE SIGN CHANGE)
!jh      A4 = -F(2,2)

!**************************************************************
!*     changed by john hubbert sept. 22, 1993
!*************************************************************
!      NEW SET TRANSFORMATION MATRIX VALUES (sept. 22, 1993)

!     F(1,1):  PARALLEL COMPONENT, PERPENDICULAR INCIDENT POLARIZATION
!              (Svh amplitude) one sign change

      a3 = -f(1,1)

!     F(1,2):  PARALLEL COMPONENT, PARALLEL INCIDENT POLARIZATION
!               (Svv amplitude) one sign change

      a2 = -f(1,2)

!     F(2,1):  PERPENDICULAR COMPONENT, PERPENDICULAR INCIDENT POLARIZATION
!               (Shh amplitude) one sign change
      a1 = -f(2,1)

!     F(2,2):  PERPENDICULAR COMPONENT, PARALLEL INCIDENT POLARIZATION
!               (Shv amplitude) one sign change
      a4 = -f(2,2)

!    Therefore we have:

       shh=a1
       shv=a4
       svh=a3
       svv=a2

!    where shh, shv, svh, and shh are members of the FSA (optical)
!    backscatter matrix.

       m1=cabs(shh)**2
       m2=cabs(svv)**2
       m3=cabs(shv)**2
       m4=cabs(svh)**2

       m(1,1) = (m1 + m2 + m3 + m4)/2
       m(1,2) = (m1 - m2 - m3 + m4)/2
       m(2,1) = (m1 - m2 + m3 - m4)/2
       m(2,2) = (m1 + m2 - m3 - m4)/2

       m(1,3) = dble (shh*conjg(shv) + svh*conjg(svv))
       m(1,4) = aimag(shh*conjg(shv) + svh*conjg(svv))
       m(2,3) = dble (shh*conjg(shv) - svh*conjg(svv))
       m(2,4) = aimag(shh*conjg(shv) - svh*conjg(svv))

       m(3,1) = dble (shh*conjg(svh) + shv*conjg(svv))
       m(3,2) = dble (shh*conjg(svh) - shv*conjg(svv))
       m(3,3) = dble (shh*conjg(svv) + shv*conjg(svh))
       m(3,4) = aimag(shh*conjg(svv) + svh*conjg(shv))

       m(4,1) = aimag(svh*conjg(shh) + svv*conjg(shv))
       m(4,2) = aimag(svh*conjg(shh) + shv*conjg(svv))
       m(4,3) = aimag(svv*conjg(shh) + svh*conjg(shv))
       m(4,4) = dble (svv*conjg(shh) - shv*conjg(svh))
!
!******************************************************************
!    end change
!    begin old code
!*****************************************************************
!     CALCULATE CONSTANTS TO BE USED IN MUELLER MATRIX CALCULATION --
!     -> M TERMS:

      m1 = cabs(a1)*cabs(a1)
      m2 = cabs(a2)*cabs(a2)
      m3 = cabs(a3)*cabs(a3)
      m4 = cabs(a4)*cabs(a4)
!     -> S TERMS:
!jh      S21 = REAL(A2)*REAL(A1) + AIMAG(A2)*AIMAG(A1)
!jh      S23 = REAL(A2)*REAL(A3) + AIMAG(A2)*AIMAG(A3)
!jh      S24 = REAL(A2)*REAL(A4) + AIMAG(A2)*AIMAG(A4)
!jh      S31 = REAL(A3)*REAL(A1) + AIMAG(A3)*AIMAG(A1)
!jh      S34 = REAL(A3)*REAL(A4) + AIMAG(A3)*AIMAG(A4)
!jh      S41 = REAL(A4)*REAL(A1) + AIMAG(A4)*AIMAG(A1)
!     -> D TERMS:
!jh      D21 = REAL(A2)*AIMAG(A1) - AIMAG(A2)*REAL(A1)
!jh      D23 = REAL(A2)*AIMAG(A3) - AIMAG(A2)*REAL(A3)
!jh      D24 = REAL(A2)*AIMAG(A4) - AIMAG(A2)*REAL(A4)
!jh      D31 = REAL(A3)*AIMAG(A1) - AIMAG(A3)*REAL(A1)
!jh      D34 = REAL(A3)*AIMAG(A4) - AIMAG(A3)*REAL(A4)
!jh      D41 = REAL(A4)*AIMAG(A1) - AIMAG(A4)*REAL(A1)

!     CALCULATE MUELLER MATRIX PARAMETERS --
!     -> ROW 1:
!jh      M(1,1) = (M2 + M3 + M4 + M1)/2.0
!jh      M(1,2) = (M2 - M3 + M4 - M1)/2.0
!jh      M(1,3) = S23 + S41
!jh      M(1,4) = -D23 - D41
!     -> ROW 2:
!jh      M(2,1) = (M2 + M3 - M4 - M1)/2.0
!jh      M(2,2) = (M2 - M3 - M4 + M1)/2.0
!jh      M(2,3) = S23 - S41
!jh      M(2,4) = -D23 + D41
!     -> ROW 3:
!jh      M(3,1) = S24 + S31
!jh      M(3,2) = S24 - S31
!jh      M(3,3) = S21 + S34
!jh      M(3,4) = -D21 + D34
!     -> ROW 4:
!jh      M(4,1) = D24 + D31
!jh      M(4,2) = D24 - D31
!jh      M(4,3) = D21 + D34
!jh      M(4,4) = S21 - S34

!     MULTIPLY MUELLER-MATRIX TERMS BY APPROPRIATE SCALING FACTOR
 return
 end  subroutine mmat

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine crossec
 implicit none

!
!     SCATTERING, EXTINCTION AND ABSORPTION CROSS SECTIONS
!     here units are in
!
 real(fpd) :: wvnmbr, abi, abq


 wvnmbr = 2.*pi/lambda
 schpt = schpt*16.0*pi/(wvnmbr*wvnmbr)  !scattering cross section (H) 
 exhpt = exhpt*16.0*pi/(wvnmbr*wvnmbr)  !extinction cross section (H)
 scvpt = scvpt*16.0*pi/(wvnmbr*wvnmbr)  !scattering cross section (V)
 exvpt = exvpt*16.0*pi/(wvnmbr*wvnmbr)  !extinction cross section (V)

 abh = exhpt - schpt  !absorption cross section (H)
 abv = exvpt - scvpt  !absorption cross section (V)
 if (abh .lt. 1.0e-08) abh = 0.0
 if (abv .lt. 1.0e-08) abv = 0.0
 abi = abv + abh
 abq = abv - abh

 if(debug_mueller) write(io_check_out,8001) abi,abq
 8001 format('abi=',d15.7,1x,'abq=',d15.7)

 return
 end subroutine crossec

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine tormp (mxelev_in, elev_in, nang)
 implicit none

! *** THIS ROUTINE DOES THE REQUIRES ADDITION OF THE MUELLER
! *** AND EXTINCTION MATRICES CORRESPONDING TO THE DIFFERENT
! *** PARTICLE GROUPS (TYPES) BEFORE CALLING ROUTINE RMP
! ***                                                    RAVI: 12/24/91
!
! *** THESE ARE THE ARRAYS THAT CONTAIN THE FINAL MUELLER AND
! *** EXTINCTION MATRICES.

 integer,intent(in) :: mxelev_in              ! maximum dimension of elevation angle 
 real(fpd),intent(in) :: elev_in( mxelev_in ) ! elevation angles [deg]
 integer :: nang
 integer :: m, k,ir,ic

 type ( polradar_parameter ) :: output ! polarimeterica radar measureable


! eldeg = eastrt

 loop_elev_angle: do m = 1,mxelev_in   !loop for elevation angle

!    eldeg = elev_in(m) * degrad  ! output writing purpose [rad]
!    write(io_check_out,fmt ='(" elevation angle deg = ",f7.2)') eldeg*raddeg
!    write(io_check_out,*)

    !
    ! particle size - type integrated output is transfered into ss and fscat for final output
    !
    do k=1,nang
       do ir=1,4
          do ic=1,4
             ss(k,ir,ic) = mugrp(m,k,ir,ic)  !4x4 scattering (mueller) matrix [m2]
             if (k .eq. 1) then
                if(ir.lt.3 .and. ic.lt.3) then
                   fscat(ir,ic)= fscatp(m,ir,ic)
                endif
             endif
          enddo
       enddo
    enddo

   ! toshi
   ! 1bin width is 45 deg
   ! 1st bin = (1-1) * 45 =   0 deg (forwad scatterg)
   ! 5th bin = (5-1) * 45 = 180 deg (backward scattering) 
   !

   ! HERE THE INTEGRATED MUELLER MATRIX IS FED TO THE PROGRAM RMP
   ! NMU + 1 = 5, IS THE BACKSCATTER ANGLE (180 DEGS)
   !                                                RAVI: 12/30/91
   ! CALCULATE RADAR MEASURABLE PARAMETERS
   ! CALL RMP(NMU+1)
   !

   !
   ! input only backscattering componet of mueller matrix 
   !
    call derive_output_polradar( ss(5,1:4,1:4), fscat(1:2,1:2), lambda, output  )

 enddo loop_elev_angle

 return
 end subroutine tormp

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine derive_output_polradar( bss_in, fscat_in, lambda_in, output  )
 implicit none
!
! CALCULATE RADAR MEASURABLE PARAMETERS
! FROM MUELLER MATRIX (AT BACKSCATTER)
!
 real(fpd),intent(inout) :: bss_in(4,4)  ! backscattering component of 4x4 mueller matrix [m2] (actually m2/m3)
 complex,intent(in)   :: fscat_in(2,2)   ! 2x2 forward scattering matrix [m] (re, im)          (actually m/m3)
 real(fpd),intent(in) :: lambda_in       ! wavelength [m]
 type ( polradar_parameter ),intent(out) :: output ! polarimeterica radar measureable


 integer :: n_htyp,i,j,k    ! loop
 complex :: covarm(3,3),rhoh,rhov
 complex :: tm(2,2)
 complex :: fscat_temp(2,2)
 real(fpd) :: exmpt(4,4)  !4x4 extinction matrix K
 real(fpd) :: hnorm
 real(fpd) :: atth, attv, datt, delhv, denom, dpha, ldr, magksq, mrho, nfac, num, &
              zdp, zdpt, zdr, zhh, degpolh, degpolv


!
!
!
 fscat_temp = fscat_in

 if(debug_mueller) print *,'in derive_output_polradar'
!
!
! AS IN RADAR CASE MAGINITUDE OF |K|^2 IS TAKEN TO BE .92
!
 magksq = .92    !toshii --> must be user input later

!
!     1.e6 is the factor to convert from cm units to mm
!     1.e18 is the factor to convert from m units to mm
!
 nfac =1.d18*4.*pi*(lambda_in**4)/((pi**5)*(magksq))

!***********************************************************8
!     changed by john hubbert sept 22, 1993
!     calculate the covariance matrix << IN BSA >>
!***********************************************************

! mag(Shh)**2
 covarm(1,1) =cmplx( nfac*(bss_in(1,1) + bss_in(1,2) + &
                     bss_in(2,1) +bss_in(2,2))/2 , 0.0)

! mag(Svv)**2
 covarm(3,3) = cmplx( nfac*(bss_in(1,1) - bss_in(1,2) - &
                      bss_in(2,1) + bss_in(2,2))/2, 0.0)

! 2*mag(Shv)**2
 covarm(2,2) = cmplx( nfac*(bss_in(1,1) - bss_in(2,2)), &
                            0.0)

! sqrt(2)Shh*conjg(Shv)
 covarm(1,2) = nfac*sqrt(2.)*(cmplx(bss_in(1,3) + bss_in(2,3), &
                              +bss_in(1,4) + bss_in(2,4)))/2

! Shh*conjg(Svv)
 covarm(1,3) = nfac*(cmplx(-bss_in(3,3) - bss_in(4,4), &
                      bss_in(4,3) - bss_in(3,4)))/2

 covarm(2,1) = conjg(covarm(1,2))

! sqrt(2)Shv*conjg(Svv)
 covarm(2,3) = nfac*sqrt(2.)*(cmplx(bss_in(1,3) - bss_in(2,3), &
                              bss_in(1,4) - bss_in(2,4)))/2

 covarm(3,1) = conjg(covarm(1,3))

 covarm(3,2) = conjg(covarm(2,3))

!  THE COVARIANCE MATRIX IS NOW CONJUGATED. THE SCATTERING PROGRAM
!  USES THE +e^-jwt and polarization definitions demand e^jwt. This
!  MEANS THAT THE COVARIANCE MATRIX SHOULD BE CONJUGATED. THE RESULT
!  IS THAT -delta IS CHANGED TO +delta.
!  -----john hubbert 7 june, 1994
!  subroutine to conjugate the covariance matrix and the mueller matrix
!  to obtain an +jwt convention

! call conjugate(covarm,ss_in,nang)
 do k=1,3
    do j=1,3
       covarm(k,j) = conjg(covarm(k,j))
    enddo
 enddo

 do i=1,4
    do j=1,4
       if(j .eq. 4) bss_in(i,j) = -(bss_in(i,j))
       if(i .eq. 4) bss_in(i,j) = -(bss_in(i,j))
! there used be a whole mueller matrix 
!       if(j .eq. 4) ss_in(1,i,j)    = -(ss_in(1,i,j))
!       if(i .eq. 4) ss_in(1,i,j)    = -(ss_in(1,i,j))
    enddo
 enddo

!
!  ABSOLUTE REFLECTIVITY IN LOG(MM^6/M^3)
!jh       ZHH = 0.5*NFAC*(SS(NANG,1,1) - SS(NANG,1,2) - SS(NANG,2,1) +
!jh     >           SS(NANG,2,2))

 zhh = cabs(covarm(1,1))
 if (zhh .gt. 0.0) then
     zhh = 10.0*dlog10(zhh)
 else
     zhh = undefined_dble  !undefined value
 endif
 output%zhh = zhh 

!
! if zhh is undefined, all undefined
!
 if( output%zhh == undefined_dble ) then
    output%zhh     = undefined_dble
    output%zdr     = undefined_dble
    output%kdp     = undefined_dble
    output%rhohv   = undefined_dble
    output%ldr     = undefined_dble
    output%zdp     = undefined_dble
    output%delhv   = undefined_dble
    output%datt    = undefined_dble
    output%atth    = undefined_dble
    output%attv    = undefined_dble
    output%degpolh = undefined_dble
    output%degpolv = undefined_dble
    return
 endif 


!
!     DIFFERENTIAL REFLECTIVITY
!jh   NUM = SS(NANG,1,1) - SS(NANG,1,2) - SS(NANG,2,1) + SS(NANG,2,2)
!jh   DENOM = SS(NANG,1,1) + SS(NANG,1,2) + SS(NANG,2,1) + SS(NANG,2,2)

!
!  num = bss_in(1,1) - bss_in(1,2) - bss_in(2,1) + bss_in(2,2)
!  demon = bss_in(1,1) + bss_in(1,2) + bss_in(2,1) + bss_in(2,2)
!
 num = cabs(covarm(1,1))
 denom = cabs(covarm(3,3))
 if( num/denom > 0. ) then 
   zdr = 10.0*dlog10(num/denom)
 else
   zdr = undefined_dble
 endif
 output%zdr = zdr

!
!     DIFFERENCE REFLECTIVITY
! *** THIS IS TO ACCOUNT FOR ZDPT=0. AT AN ELEVATION ANGLE = 90DEG
! ***                                                   RAVI  12/23/91.
!     ldc 10/20/95:  Zdp's are not correct (too high). I believe problem
!     came about when jh changed method of calculating observables from
!     the covariance terms.  The NFAC*0.5 factor is included above in
!     the covarm(1,1) and covarm(3,3) terms and does not need to be
!     repeated as is below.
!     old code:
!
!      ZDPT = NUM - DENOM
!      IF(ZDPT .GT. 0.) THEN
!         ZDP = 10.*DLOG10(.5*NFAC*ZDPT)       !what is ZDP?
!      ELSEIF(ZDPT .LT. 0.) THEN
!         ZDP = -10.*DLOG10(.5*NFAC*(-ZDPT))
!      ELSE
!         ZDP = 0.
!      ENDIF
!
!     ldc 10/20/95:  Corrected code for Zdp
 zdpt = num - denom
 if(zdpt .gt. 0.) then
    zdp = 10.*dlog10(zdpt)
 elseif(zdpt .lt. 0.) then
    zdp = -10.*dlog10(-1.*(zdpt))
 else
    zdp = 0.
 endif
 output%zdp = zdp


!
!     LINEAR DEPOLARIZATION RATIO
!jh   NUM = SS(NANG,1,1) - SS(NANG,1,2) + SS(NANG,2,1) - SS(NANG,2,2)
!jh   DENOM = SS(NANG,1,1) - SS(NANG,1,2) - SS(NANG,2,1) + SS(NANG,2,2)
!
 denom = cabs(covarm(1,1))
 num = cabs(covarm(2,2)/2)
 ldr = num/denom
 if (ldr .gt. 0.0) then
     ldr = 10.0*dlog10(ldr)
 else
     ldr = undefined_dble
 endif
 output%ldr = ldr

!
!     RHO
!jh   MAG1 = SS(NANG,1,1) - SS(NANG,1,2) - SS(NANG,2,1) + SS(NANG,2,2)
!jh   MAG2 = SS(NANG,1,1) + SS(NANG,1,2) + SS(NANG,2,1) + SS(NANG,2,2)
!jh   DENOM = SQRT(MAG1)*SQRT(MAG2)
!     -1 IS INTRODUCED IN RHORE BECAUSE OF S(3,3), S(4,4) ARE -VE IN
!     SCATTERING MATRIX. S(3,3) AND S(4,4) SHOULD BE +VE .
!jh   RHORE =-1.*((SS(NANG,3,3) + SS(NANG,4,4)))/DENOM
!jh   RHOIM = -(SS(NANG,4,3) - SS(NANG,3,4))/DENOM
!jh   RHOHV = CMPLX(RHORE,RHOIM)
!jh   DELHV = ATAN(AIMAG(RHOHV)/REAL(RHOHV))*RADDEG
!jh   MRHO = CABS(RHOHV)
!
! modified to avoid NaN output. 
 if( cabs(covarm(3,1))  == 0. ) then
   mrho = 0.
 else
   mrho = cabs(covarm(3,1))/ sqrt(cabs(covarm(1,1))*cabs(covarm(3,3)))
 endif
 if( (mrho > const_largest_real) .or. &  !check if infinity
     (mrho < -1*const_largest_real)     ) then
        mrho = 0.
  endif

 output%rhohv = mrho


!
!      delhv= atan2( aimag(covarm(3,1)), dble(covarm(3,1) ))*RADDEG
! toshi bug fixxed below
!
 delhv = atan2( dble(aimag(covarm(3,1))), dble( covarm(3,1) ) )*raddeg
 output%delhv = delhv

!
! also calculate the degree of polarization for vertical and
! horizontal incidence
! for horizontal incidence
!
 num = 4*(dble(covarm(1,1)*covarm(2,2))/2 - &
          dble(covarm(1,2)*covarm(2,1))/2)
 denom =dble((covarm(1,1)+covarm(2,2)/2))**2
 degpolh = sqrt( 1. - num/denom)
 output%degpolh = degpolh

! for vertical incidence

 num = 4*(dble(covarm(3,3)*covarm(2,2)/2) - &
          dble(covarm(2,3)*covarm(3,2))/2)
 denom = dble((covarm(3,3)+covarm(2,2)/2))**2
 degpolv = sqrt( 1. - num/denom)
 output%degpolv = degpolv

!  new degree of pol. based on cov. matrix
!  calculate degree of coherence
!  (toshii tends to be overflow)
! rhoh= covarm(1,2)/  (sqrt(covarm(1,1)*covarm(2,2))) 
! rhov= covarm(2,3)/  (sqrt(covarm(3,3)*covarm(2,2))) 
! output%rhoh = rhoh
! output%rhov = rhov

!
!     GET EXTINCTION MATRIX  
!
 if(transflag .eq. 1)then
    call transmtx(fscat_temp,tm,lambda_in)
 endif

!
!     conjugate is to put the foward scattering amplitudes in
!     a +jwt time convention, to be consistent with other outputs
!
 do j=1,2
    do k=1,2
       fscat_temp(j,k)=conjg(fscat_temp(j,k))
    enddo
 enddo

!
! CALCULATE EXTINCTION MATRIX FROM FORWARD AMPLITUDE MATRIX 'F'
! RETURNS 4-BY-4 REAL EXTINCTION MATRIX IN THE ARRAY 'EX' [m2]
!
 call emat(fscat_temp, exmpt, lambda_in)

!
!     since the extinction matrix is given in m^2/m^3, the diff
!     phase needs to be mulitiplied by 1000 to obtain deg/km
!     ldc 10/20/95:  After consulting with jb about the negative Kdp's
!                    in rain, he instructed me to make the following change
!                    (Reference e-mail from 10/20/95)
!     old code for Kdp
!         DPHA = RADDEG*EXMPT(3,4)*1000.0
!     ldc 10/20/95: corrected code for Kdp
! toshi
!
! raddeg = 180/pi [radince to degree]
! exmpt(4,3) = Im(Mhh - Mvv)  [m2/m3]
!

 dpha = raddeg*exmpt(4,3)*1000.0
 output%kdp = dpha     ![deg/km]


!     old code for DATT
!         DATT = -.4343*EXMPT(1,2)
!     ldc 10/24/95: corrected code for Differential Attenuation
!                   as instructed by jb.
!     ldc 10/26/95: multiply by 1000 to get db/km (according to jb)
!
!         DATT = 1000.0*0.4343*EXMPT(1,2)
!     ldc 03/07/96:  values for DATT (and ATTH, ATTV below) seem to
!                    be low by an 1 order of magnitude (factor of 10)
!                    compared to Bringi et al. (1990,JAOT) for gamma DSD's.
!                    I've explored literature and it seems to me
!                    that it should be 10/ln(10)*EXMPT(1,2)*1000
!                    assuming EXMPT is in m**-1 (where 10/ln(10)=4.343).
!                    NOTE for me:  LOG=ln
!
!    ldc 03/13/96:   Since lambda_in is converted to m and f is in meters
!                    (see sub EMAT), then EXMPT is in m**2.  Therefore,
!                    the correct multiplictive constant is indeed given
!                    by 4.343E+03 = [10/1n(10)]*1000 in order to get
!                    dB/km for DATT, ATTH, and ATTV for EXMPT in m**2.
!                    See Doviak and Zrnic (1993), Oguchi (1983), or
!                    Battan (1973, watch units).  Notice that all past
!                    Bringi and collaborator papers have 0.4343 as the
!                    constant...but this is for EXMPT in cm**2.

!
!          DATT = (10.0/LOG(10.0))*1000.0*EXMPT(1,2)
!     ldc 04/05/96:  multiplicative factor should be 20/log(10) (see below)
!

 datt = (20.0/log(10.0))*1000.0*exmpt(1,2)

!
!     ldc 03/04/96:  added code for Ah: specific horizontal attenuation
!                                   Av: specific vertical attenuation
!                    [dB/km]
!
!          ATTH = 0.4343*(0.5)*(EXMPT(1,1) + EXMPT(1,2))*1000.0
!
!          ATTV = 0.4343*(0.5)*(EXMPT(1,1) - EXMPT(1,2))*1000.0
!     ldc 03/07/96:  fix code for Ah and Av for order of magnitude
!                    problem discussed above.
!          ATTH = (10.0/LOG(10.0))*(0.5)*(EXMPT(1,1) +
!     +            EXMPT(1,2))*1000.0
!
!          ATTV = (10.0/LOG(10.0))*(0.5)*(EXMPT(1,1) -
!     +            EXMPT(1,2))*1000.0
!     ldc 04/05/96:  After comparison with Bringi et al (1990),
!                    I came to the conclusion that the results
!                    for ATTH,ATTV, and DATT are too low.  After
!                    consultation with jb and bringi, determined
!                    that it should be 20/log(10)*1000 because
!                    we are dealing with E fields such that P ~ E**2
!                    so 10 log P = 20 log E.  Tadaa....

 atth = (20.0/log(10.0))*(0.5)*(exmpt(1,1) + exmpt(1,2))*1000.0
 output%atth = atth

 attv = (20.0/log(10.0))*(0.5)*(exmpt(1,1) - exmpt(1,2))*1000.0
 output%attv = attv

!
! output for checking
!
 if(debug_mueller) then
  print*,"zh:",zhh
  print*,"zdr:",zdr
  print*,"dpha:",dpha
  print*,"zdp:",zdp
  print*,"mrho:",mrho
  print*,"delhv:",delhv
  print*,"atth:",atth
  print*,"datt:",datt
  print*,"ldr:",ldr
 endif


!         write(io_scatmat,*)"4x4 extinction matrix units are m^2"
!         write(io_scatmat,'(x,4d14.6,x,d14.6,x,a)') (exmpt(1,j),j=1,4)
!         write(io_scatmat,'(x,4d14.6,d14.6,x,a)') (exmpt(2,j),j=1,4)
!         write(io_scatmat,'(x,4d14.6,x,f7.2,x,f5.2,x,a,x,a)') (exmpt(3,j),j=1,4)
!         write(io_scatmat,'(x,4d14.6,x,d14.6,x,a)') (exmpt(4,j),j=1,4)

!         write(io_scatmat,*)"4x4 forward scattering matrix units are m^2"
!         write(io_scatmat,'(x,4d14.6,x,d14.6,x,a)') (ss_in(1,1,j),j=1,4)
!         write(io_scatmat,'(x,4d14.6,x,d14.6,x,a)') (ss_in(1,2,j),j=1,4)
!         write(io_scatmat,'(x,4d14.6,x,d14.6,x,a)') (ss_in(1,3,j),j=1,4)
!         write(io_scatmat,'(x,4d14.6,x,d14.6,x,a)') (ss_in(1,4,j),j=1,4)

!         write(io_scatmat,*)"4x4 backscatter matrix units are m^2"
!         write(io_scatmat,'(x,4d14.6,2x,d14.6,x,a)') (bss_in(1,j),j=1,4)
!         write(io_scatmat,'(x,4d14.6,2x,d14.6,x,a)') (bss_in(2,j),j=1,4)
!         write(io_scatmat,'(x,4d14.6,x,d14.6,x,a)') (bss_in(3,j),j=1,4)
!         write(io_scatmat,'(x,4d14.6,x,d14.6,x,a)') (bss_in(4,j),j=1,4)

!         write(io_scatmat,*) "2x2 fowrward scattering matrix is (m) (re, im):"
!         do j=1,2
!            write(io_scatmat,1011)(fscat_temp(j,k), k=1,2)
!         enddo

!         write(io_scatmat,*) "2x2 transmission  matrix is (re, im):"
!         do j=1,2
!            write(io_scatmat,1011)(tm(j,k), k=1,2)
!         enddo

! 1011    format(4e13.5)
!        output covariance matrix

!         hnorm=cabs(covarm(1,1))
!         do j=1,3
!            do k=1,3
!               covarm(j,k)=covarm(j,k)/hnorm
!            enddo
!         enddo

!         write(io_scatmat,*)"3x3 covariance matrix and results (re, im):"
!         do j=1,3
!            write(io_scatmat,1010) (covarm(j,k), k=1,3)
!         enddo

! 1010    format(6e13.5)
!        compute and write magnitude and phase of the covariance matrix
!         write(io_scatmat,*)"3x3 covariance matrix and results (mag, phase):"

!         do j=1,3
! debbug toshi (add dble before aimag)
!            write(io_scatmat,1010)(cabs(covarm(j,k)), atan2(dble(aimag(covarm(j,k))), &
!                           dble(covarm(j,k)) ),  k=1,3)
!         enddo
! 1014    format(2e13.5)

!     ldc 3/5/96:  output radar observables to separate file w/ no text

!      write(io_radarobs,1020) zhh,zdr,dpha,zdp,mrho,delhv,atth,datt,ldr
! 1020  format( f7.2,1x,f7.3,1x,f9.5,1x,f7.2,1x,f6.4,1x,f7.2,1x, f10.8,1x,f11.9,1x,f6.2)

!      write(io_scatmat,*)
!      write(io_scatmat,*)"  radar parameters"
!      write(io_scatmat,fmt ='(d22.13,8x,"reflectivity dbz")')zhh
!      write(io_scatmat,fmt ='(d22.13,8x,"differential reflectivity db")')zdr
!      write(io_scatmat,fmt ='(d22.13,8x,"linear depol ratio db")')ldr
!      write(io_scatmat,fmt ='(d22.13,8x,"specific diff phase deg/km")')dpha
!      write(io_scatmat,fmt ='(d22.13,8x,"difference reflectivity db")')zdp
!      write(io_scatmat,fmt ='(d22.13,8x,"hv correlation")')mrho
!      write(io_scatmat,fmt ='(d22.13,8x,"phaseshift upon backscatter deg")') delhv
!      write(io_scatmat,fmt ='(d22.13,8x,"differential attenuation db/km")') datt
!      write(io_scatmat,fmt ='(d22.13,8x,"specific hor attenuation db/km")') atth
!      write(io_scatmat,fmt ='(d22.13,8x,"specific ver attenuation db/km")') attv
!      write(io_scatmat,*)
!      write(io_scatmat,*)
!      write(io_scatmat,fmt ='(d22.13,8x,"extinction cross-section h")')exhpt
!      write(io_scatmat,fmt ='(d22.13,8x,"extinction cross-section v")')exvpt
!      write(io_scatmat,fmt ='(d22.13,8x,"scattering cross-section h")')schpt
!      write(io_scatmat,fmt ='(d22.13,8x,"scattering cross-section v")')scvpt
!      write(io_scatmat,*)
!      write(io_scatmat,fmt ='("  deg. polarization h",d22.13,8x)')degpolh
!      write(io_scatmat,fmt ='("  deg. polarization v",d22.13,8x)')degpolv
!      write(io_scatmat,fmt ='("  deg. coherence h",2d22.13,8x)')rhoh
!      write(io_scatmat,fmt ='("  deg. coherence v",2d22.13,8x)')rhov
!      write(io_scatmat,*)

!        Output 2x2 transmission matrix information to file
!         write(io_scatmat,1007)npart  ! # of particle types

!
! distributions type name (output purpose) --> remove later..
!
! if (distyp .eq. 0)then
!     orientmd(n_htyp) = 'random'
! elseif(distyp .eq. 1)then
!     orientmd(n_htyp) = 'simple harmonic'
! elseif(distyp .eq. 2)then
!     orientmd(n_htyp) = 'quasi gaussian'
! elseif(distyp .eq. 3)then
!     orientmd(n_htyp) = 'langevin'
! endif

!
! distributions type name (output purpose) --> remove later..
!
! if (distyp .eq. 0)then
!     orientmd(1) = 'random'
! elseif(distyp .eq. 1)then
!     orientmd(1) = 'simple harmonic'
! elseif(distyp .eq. 2)then
!     orientmd(1) = 'quasi gaussian'
! elseif(distyp .eq. 3)then
!     orientmd(1) = 'langevin'
! endif
!
! global parameters for output purpose
!
! medium = ''
! mdfile = ''
! dpar1md(n_htyp)  = dpar1
! dpar2md(n_htyp)  = dpar2
! dpar3md(n_htyp)  = dpar3
! dminmd(n_htyp)   = dmin
! dmaxmd(n_htyp)   = dmax
! dstepmd(n_htyp)  = dstep
! mediummd(n_htyp) = medium

!         do n_htyp=1,npart !# of particle types
!            write(io_scatmat,1008) n_htyp,mediummd(n_htyp)
!            write(io_scatmat,*)'the orientation distribution is ',orientmd(n_htyp)
!            if(distyp .eq. 3)then
!               write(io_scatmat,*)'with a kappa of ', dpar1md(n_htyp)
!               write(io_scatmat,*)'mean theta is ',dpar2md(n_htyp)
!            else
!               write(io_scatmat,*)'with a mean orienting angle (theta) of ', &
!                        dpar1md(n_htyp)
!               write(io_scatmat,*)'and a standard deviation of ',dpar2md(n_htyp)
!            endif
!            if(distyp .eq. 3)then
!               write(io_scatmat,*)'mean phi is ',dpar3md(n_htyp)
!            else
!               if(phiflag .eq. 1)then
!                  write(io_scatmat,*)'phi angle is fixed at ',phiang
!               else
!                  write(io_scatmat,*)'phi angle is random'
!               endif
!            endif
!            write(io_scatmat,*)'elevation angle is ',eldeg*raddeg
!            write(io_scatmat,*)'gamma model used for size dist.'
!            write(io_scatmat,1000)dminmd(n_htyp),dmaxmd(n_htyp),dstepmd(n_htyp)
!            write(io_scatmat,999)n0(n_htyp), d0(n_htyp), mgam(n_htyp)
!            write(io_scatmat,*)
!         enddo



! 1000    format('dmin = ',f5.3,'(cm)',2x,'dmax = ',f5.3,'(cm)' &
!      ,2x,'dstep =',f5.3,'(cm)')
! 1007    format(i1,1x,'particle types in this resolution volume')
! 1008    format('particle',i2,1x,'is ',a)
! 999     format('n0 = ',f16.1,'(/cm/m^3)',2x,'d0 = ',f6.3,'(cm)' &
!      ,2x,'m = ',f4.1)
! 900  format (1x,'zhh = ',f6.1,' db')
! 901  format (1x,'zdr = ',f6.1,' db')
! 902  format (1x,'ldr = ',f6.1,' db')
! 903  format (1x,'rho (h,v) = ',2(f6.2,2x))
! 904  format (/,20x,'null',9x,'spread    depth',/,15x,'(mag)   (phase)')
! 905  format (1x,'co-pol:     ',f8.3,2x,f7.2,2(3x,f6.3))
! 906  format (1x,'cross-pol:  ',f8.3,2x,f7.2,2(3x,f6.3))
 return
 end subroutine derive_output_polradar

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine emat(f, ex, lambda_in)
 implicit none
!
!     CALCULATE EXTINCTION MATRIX FROM FORWARD AMPLITUDE MATRIX 'F'
!
!     RETURNS 4-BY-4 REAL EXTINCTION MATRIX IN THE ARRAY 'EX'
!

 complex,intent(inout) :: f(2,2)    ! 2x2 foward scattering amplitude [m]
 real(fpd),intent(in) :: lambda_in  ! wavelength [m]
 real(fpd),intent(out) :: ex(4,4)   ! 4x4 extinction matrix [m2]

 complex :: mhh,mhv,mvh,mvv,i
 data i/(0.,1.)/

!     i=cmplx(0.,1.)


!     ldc  03/13/96:  lambda_in is in meters (m).  f was already in m. -->if you consider per volume --> [m/m3]
!                     As a result, mii is in m**2 (and therefore EX(i,j)
!                     or EXMPT(i,j) also in m**2). --> if yo consder per volume --. [m2/m3]
! lambda_in is 2*pi/k component in Viveck paper. toshi
!

 f(1,1)=lambda_in*f(1,1)
 f(1,2)=lambda_in*f(1,2)
 f(2,1)=lambda_in*f(2,1)
 f(2,2)=lambda_in*f(2,2)

 mhh=i*f(1,1)
 mhv=i*f(1,2)
 mvh=i*f(2,1)
 mvv=i*f(2,2)

 ex(1,1) = dble(mhh + mvv)
 ex(2,2) = ex(1,1)
 ex(3,3) = ex(1,1)
 ex(4,4) = ex(1,1)

 ex(1,2) = dble(mhh-mvv)
 ex(2,1) = ex(1,2)

 ex(3,4) = aimag(mvv - mhh)
 ex(4,3) = -ex(3,4)

 ex(1,3) = dble(mvh + mhv)
 ex(3,1) = ex(1,3)

 ex(1,4) = aimag(mhv - mvh)
 ex(4,1) = -ex(1,4)

 ex(2,3) = dble(mhv-mvh)
 ex(3,2) = -ex(2,3)

 ex(2,4) = aimag(mhv+mvh)
 ex(4,2) = -ex(2,4)

 return
 end subroutine emat

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine transmtx(fscat_in,tm_out,lambda_in)
 implicit none

!     A subroutine to compute the 2x2 transmission matrix using
!     using Oguchi's formulations found in 'Electromagnetic Wave
!     Propagation and Scattering in Rain and Other Hydrometeors'
!     Proceedings of IEEE, vol71, No. 9, Sept 1983'
!     Addition by John Beaver November 8, 1994.
!
! toshii--> This routine has many warning of precision conversion from double complex to single complex. 
!       m11=-j*ko-j*((2*pi)/ko)*shh
!          1
! Warning: Possible change of value in conversion from COMPLEX(8) to COMPLEX(4) at (1) [-Wconversion]
!
!
 complex,intent(in) :: fscat_in(2,2)
 real(fpd),intent(in) :: lambda_in
 complex,intent(out) :: tm_out(2,2)


 real(fpd) :: ko
 complex :: m11,m22,m12,m21,eigen1,eigen2,j
 complex :: g,phi0,arg, shh, svv, shv, svh

! common /trans1/ dist_trans,transflag

 ko=(2.*pi)/lambda_in
 j=(0,1.)
!     conjugate is to convert from -jwt to +jwt  Oguchi's
!     equations are in the +jwt
 shh = conjg(fscat_in(1,1))
 svv = conjg(fscat_in(2,2))
 shv = conjg(fscat_in(1,2))
 svh = conjg(fscat_in(2,1))
 m11=-j*ko-j*((2*pi)/ko)*shh
 m22=-j*ko-j*((2*pi)/ko)*svv
 m12=-j*((2*pi)/ko)*svh
 m21=m12

!
! compute the eigenvalues
!
 eigen1=((m11+m22)/2.)-csqrt((((m11-m22)/2.)**2)+m12*m21)
 eigen2=((m11+m22)/2.)+csqrt((((m11-m22)/2.)**2)+m12*m21)

!
!     compute phi and g
!
 arg = (2.*m12)/(m11-m22)   !toshii (this function often got invalid floating process with debbugging mode). Got to fix... 
 phi0=.5*(.5*j*clog((j+arg)/(j-arg)))


!
!     subtract 90  degrees if phi0 is positive(tilt angle is greater than 45)
!     For example, the above formula for phi0 will give an angle of 40
!     for a -50 degree tilt(oguchi notation) due to the nature of the log
!     function.
!
 if(dble(phi0) .gt. 0) phi0= phi0-90.*degrad
 g = cexp((eigen1 - eigen2)*cmplx(dist_trans))

!
!     compUTE THE TRANSMISSION MATrix
!
 tm_out(1,1) = (cos(phi0)**2)*(cexp(eigen2*cmplx(dist_trans))*(g + &
                              (csin(phi0)/ccos(phi0))**2))
 tm_out(2,2) = (cos(phi0)**2)*(cexp(eigen2*cmplx(dist_trans))*(1 + &
                            g*(csin(phi0)/ccos(phi0))**2))
 tm_out(1,2) = -(cos(phi0)**2)*(cexp(eigen2*cmplx(dist_trans))*(1 - g)* &
                               (csin(phi0)/ccos(phi0)))
 tm_out(2,1) = tm_out(1,2)

!     tm_out(1,1) = Thh
!     tm_out(2,2) = Tvv
!     tm_out(1,2) = Thv
!     tm_out(2,1) = Tvh
!     transmission matrix is output in a +jwt, FSA convention
!     using a (h,v) basis. Output file will be hh,hv,vh,vv

 return
 end subroutine transmtx

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 subroutine zero_output(nang)
 implicit none

 integer,intent(in) :: nang

 integer :: iz,ir,ic

!        ZERO AVERAGE MUELLER-MATRIX ACCUMULATOR

  do  iz = 1,nang
      do  ir = 1,4
          do  ic = 1,4
              s(iz,ir,ic) = 0.0
          enddo 
      enddo 
  enddo 

  fscat(1,1)=cmplx(0.)
  fscat(1,2)=cmplx(0.)
  fscat(2,1)=cmplx(0.)
  fscat(2,2)=cmplx(0.)

  scvp = 0.0
  schp = 0.0
  exvp = 0.0
  exhp = 0.0

 return
 end subroutine zero_output

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 

 end module module_mueller_matrix

!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
!MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX MUELLER MATRIX 
