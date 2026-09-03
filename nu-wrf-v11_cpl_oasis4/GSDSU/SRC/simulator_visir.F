#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine simulator_visir(indg,inda,indt,indp,imthd,nda,na1u,am1u,na0,am0, &
                    nfi,fi,nln,thk,omg,nlgn1,g,nang,ang,phsf,epsp,epsu,galb,fsol, &
                    nplk1,cplk,bgnd,ntau,utau,scr,sci,flxd,flxu,ai_sdsu)
 use module_mpi
 use module_opt_visir, only: kna1u,kna1,kna0,kndm,knfi,klgn1,knln,kntau,&
     kplk1,knang,knln1,kndm2,pi,rad
 implicit none
!----------------------------------------------------------------------------------------------------
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
! Solve the radiative transfer in atmosphere system based on Nakajima and Tanaka (1986).
!  -use ftrn21-
! complex refractive index for ocean surface
! useful for microwave region
! by teruyuki nakajima
! 94.4.4 observation of ftrn21: 
! insufficeient correction of solar single scattering from ocean surface.
! grnd0 should be corrected in the future.
! this version of grndo calculates r(i,j) and sr(i,j0) by integrating
! the reflected intensity at ocean surface over the mu-interval
! (mu(i-1/2), mu(i+1/2)).  after solving the transfer with this
! reflection and source matrices of ocean surface, we subtracted
! the singly scattered intensity with sr(i,j0).  after that,
! we add the true single reflected intensity.
! however, this trick will have some error if mu(i)<>mu(j0),
! since the strongly peaked reflected intensity in the direction
! of -mu(j0).  this is distributed among the light along the 
! quadrature stream -mu(i).  in more accurate approximation,
! we should treat the peaked reflected intensity as a reflected
! direct radiation not distributed amond the quadrature streams.
!
! 94.9.2  the symmetric operation for r(ocean) in 'oceanrf'
! can introduce a difference of order about 0.3 % from a version
! without such symmetric operation.
!
!
! Methods: Calling sequence 
!      simulator_visir (solve the radiative transfer in atmosphere system)
!         |-- chkrt (check-in variables for -rtrn1- and set some variables)
!         |      qgausn (compute weights and abscissae for ordinary gaussian quadrature)
!         |-- lgndf3 (Legendre expansion ,same as lgndf2 but generating g-moments)
!         |     qgausn (compute weights and abscissae for ordinary gaussian quadrature)
!         |-- ftrn21 (solve the radiative transfer in atmosphere system for each fourier component)
!         |     |-- plgnd (normalized associated legendre polynomials)
!         |     |-- phas2 ()
!         |     |-- homog2 (solve the transfer in a homogeneous scattering and emitting medium)
!         |     |     |-- cpcon (machine constants of computer, overflow, underflow error)
!         |     |     |-- getqm (solve    xy = q zeig**2 inverse(q))
!         |     |     |-- csfn (get c sl and sl1 functions)
!         |     |     |-- tnvss2 (inversion of real matrix. sweep out, complete positioning)
!         |     |-- grndl3 (lambert surface reflection)
!         |     |-- grndo1 (ocean surface reflection  )
!         |     |     |-- ocnr11 (reflection matrix of ocean surface)
!         |     |     |-- ocnr31 (reflection matrix of ocean surface)
!         |     |-- trn1 (solve the radiative transfer in the multi-sublayer system by the adding method.)
!         |     |-- plgnd (normalized associated legendre polynomials)
!         |     |-- cingr (get integral constants: alfa and beta)
!         |     |-- adisc (intensity at a user defined depth)
!         |     |-- phas2 ()
!         |     |-- ainteg ()
!         |     |     |-- exint (integral dt exp(-(thk-t)/am1u-t/am0))
!         |     |     |-- csint (integral dt exp(-(dpt-t)/am1u) c(t, thk))
!         |     |-- cspl1 (getting coefficients of natural cubic spline fitting)
!         |-- convu (check convergence of intensity in the direction of mu)
!         |--intcr1 (sigle and second scattering correction assume delta-m method)
!
!
! History:
!  09/2014  Toshi Matsui@NASA GSFC : slight modification to add rt_failed flag. 
!  03/2009  Toshi Matsui@NASA GSFC : RSTAR6bF90 was modulated to SDSU framework
!  ~11/2008  Miho Sekiguchi        : F90 version of RSTAR 6b (See origional LONG history below)
!
!--- original history
! 89.11. 7 created
! 90.11.22 if(indp <= 0 .and. (inda >= 1 .and. imthd >= 1)) indpl=1 ->
!                  eq
!          change codes between -do 4 - and -do 8 -;
!          if(abs(cs1) > 1.0) cs1=sign(1.0,cs1); change loop-13,14
!          adding loop-47.
! 93. 4. 5 nsf meaning changed, set nsf=0 after loop-16
!     5. 4 naming rtrn2 from higurashi's rtrn with use of ftrn21
!           and corrected single scattering for ocean surface reflection
!    10.19 ok comparing with hasumi's program.
! 94. 9.27 debug if(am1u(iu) < 0) then a3= ...
! 95. 6. 2 gerated from rtrn2 with scr and sci
!           useful for microwave region
! 98.1.7
!   u10 > 0.01 for rstr and pstr.  u10 < 0.01 will have significant
!   error in radiances due to angular integration errors for surface
!   reflection.  recommend that u10 is set to be 0.01 for flat surface
!   case
! 08. 9.28 modified Fortran 90 free-style form
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
!----------------------------------------------------------------------------------------------------


!--- input
! indg       i      -1: no ground surface
!                    0: lambert surface
!                    1: ocean surface initialization
!                    2: ocean surface with no initialization !Do not use for SDSU since u10 is grid dependent (toshi)
!                    when indg>0 and imthd>0 then single scattering correction
!                      for ocean surface reflection
! inda       i       0: flux only.
!                    1: intensity using -am1u-.
!                    2: intensity at gaussian quadrature points.
!                       -na1u- and -am1u- are set automatically as
!                       2*nda and (-amua, +amua).
! indt       i       1: intensity/flux at user defined depth -utau-.
!                    2: intensity/flux at sublayer interfaces.
!                       -ntau- and -utau- are set automaticaaly as
!                       -nln1- and -dpt-.
! indp       i      -1: give -g- and use -g- for intensity interpolation
!                       from gaussian points to user points.
!                    0: give -g- and construct phase function for
!                       intensity interpolation.
!                    1: give phase function -phsf-.
! imthd      i      -1: nt,  0: dms-method  for intensity/flux
!                    1: ms,  2:tms,  3:ims-method for intensity.
!                    when indg>0 and imthd>0 then single scattering correction
!                      for ocean surface reflection
! nda        i       nunber of nadir-quadrature angles in the
!                    hemisphere.
! na1u       i       number of emergent nadir angles in sphere.
! am1u    r(kna1u)   consine of the emergent nadir angles.
!                      - for upward, + for downward.
! na0        i       number of solar incidences.
! am0     r(na0)     consine of solar zenith angles  >  0.
! nfi        i       number of azimuthal angles.
! fi     r(knfi)     azimuthal angles in degrees.
! nln        i       number of atmospheric sublayers.
! thk     r(knln)    optical thickness of sublayers from top to bottom.
! omg     r(knln)    single scattering albedo.
! nlgn1   i(knln)    maximum order of moments + 1.
!                      give when -indp-  <=  0.
!                      give when imthd=3 regardless of -indp-.
!                      otherwise, a value  <=  2*nda+1 is given
!                      automatically by the routine.
! g   r(klgn1,knln)  legendre moments of phase function.
!                      give when indp  <=  0.
! nang       i       number of scattering angles for phase functions.
!                      give when indp=1.
!                      give when indp=0 .and. (inda >= 1) .and.
!                                             (imthd >= 1).
! ang     r(knang)   scattering angles in degrees.
!                      give when indp=1.
!                      give when indp=0 .and. (inda >= 1) .and.
!                                             (imthd >= 1).
!                      angles should have enough resolution for legendre
!                      expansion (indp=1) and interpolation by cubic
!                      polynomial.
! phsf    r(knang,   phase function. give when indp=1.
! epsp       r       truncation criterion of phase function moments.
! epsu       r       convergence criterion of intensity (inda > 0)
! galb       r       ground albedo if indg=0
!                    u10 (m/sec)   if indg>0; u10>0.01
! fsol       r       solar irradiance at the system top.
! nplk1      i       number of order to approximate plank + 1.
!                      if 0 then no thermal.
! cplk r(kplk1,knln) plank function
!                      = sum(k1=1,nplk1) cplk(k1,l) * tau**(k1-1)
!                      tau is the optical depth measured from the top
!                      of the sublayer (same unit as fsol).
! bgnd       r       (1-galb)*b when indg=0
!                    b          when indg>0
!                    where b is plank function
!                      (same unit as fsol).
! ntau       i       number of user defined optical depths.
! utau    r(kntau)   optical depths where the field is calculted.
!                      top to bottom.
! scr       r    relative refractive index of the media
!               about 1.33 for atmosphere to ocean incidence,
!               and 1/1.33 for ocean to atmosphere incidence.
! sci       r    relative refractive index for imaginary part
!               m = cr + i*ci

!--- output
! indg               if 1 then 2
! na1u, am1u         same as 2*nda, (-auma, +amua) when inda=2.
! nlgn1, g           normalized by g(1), and cut off by the criterion
!                     -epsp-.
! ntau, utau         same as nln+1 and -dpt- when indt=1.
! phsf               normalized phase function (1 over unit sphere).
!                    reconstructed when (inda >= 1) .and. (imthd >= 1)
!                    .and. (indp=0).
! flxd    r(kna0,    downward flux at -utau-.
!           kntau)
! flxu               upward   flux at -utau-.
! ai      r(kna1u,   i(mu1u(i), mu0(j), fi(k), utau(it))
!           kna0,    intensity at -utau- (watt/m2/micron).
!           knfi,
!           kntau)
!--- correspondence between variables and parameters
! kna1u       na1u
! kna1        na1
! kna0        na0
! kndm        nda
! knfi        nfi
! knln        nln
! knln1       knln+1
! kntau       ntau
! klgn1       mxlgn1 = max(nlgn1)
!              this should be very large when ims-methd is used
!              so as to approximate p**2, or when indp=1 is assigned
!              to reconstruct the original phase function from -g-.
! klgt1       mxlgt1 = same as mxlgn1 but for that of truncated
!              variables
! kplk1       nplk1
!--

! input
  integer,intent(inout):: indg   !! surface flag (0:Lambert, 1: Ocean)
  integer,intent(in):: inda      !! radiance flag (0: flux only, 1: flux & rad)
  integer,intent(in):: indt      !! flag of optical depth (1:utau, 2:dpt)
  integer,intent(in):: indp      !! flag of phase function
  integer,intent(in):: imthd     !! flux/intensity calculation modes
  integer,intent(in):: nda       !! number of stream
  integer,intent(in):: na1u      !! number of zenith angles of emerging rad
  real(8),intent(in):: am1u(kna1u)   !! cosine of the emergent nadir angles
  integer,intent(in):: na0           !! number of solar zenith angles
  real(8),intent(in):: am0(kna0)     !! cosine of solar zenith angles
  integer,intent(in):: nfi           !! azimuthal angles
  real(8),intent(in):: fi(knfi)      !! zenith angles (degree)

  integer,intent(in):: nln           !! number of sublayers
  real(8),intent(in):: thk(knln)     !! opthical thickness
  real(8),intent(in):: omg(knln)     !! single scattering albedo
  integer,intent(out):: nlgn1(knln) !! maximum order of moments + 1 (indp<=0;imthd=3)
  real(8),intent(out):: g(klgn1,knln) !! Legendre moments of phase function (indp<=0)
  integer,intent(in):: nang          !! number of scaattering angles
  real(8),intent(in):: ang(knang)    !! scattering angles [degree]
  real(8),intent(inout):: phsf(knang,knln)  !! phase function (indp=1)
  real(8),intent(in):: epsp      !! truncation criterion of phase func. moments
  real(8),intent(in):: epsu      !! convergence criterion of intensity (inda>0)
  real(8),intent(in):: galb      !! ground albedo (indg=0)/ U10 [m/sec] (indg=1)
  real(8),intent(in):: fsol      !! solar irradiance at the system top
  integer,intent(in):: nplk1     !! number of order to approximate plank + 1
  real(8),intent(in):: cplk(kplk1,knln)  !! coefficient of planck function 
  real(8),intent(in):: bgnd      !! planck function at ground 
  integer,intent(in):: ntau      !! number of user defined optical depth
  real(8),intent(in):: utau(kntau)   !! optical depths where the field is calc
  real(8),intent(in):: scr        !! relative reflactive index (real part)
  real(8),intent(in):: sci        !! relative reflactive index (imaginary part)

! output
  real(8),intent(out):: flxd(kna0,kntau)     !! downward flux
  real(8),intent(out):: flxu(kna0,kntau)     !! upward flux

  real(8) :: ai(kna1u,kna0,knfi,kntau) !! intensity
  real(8),intent(out):: ai_sdsu(kntau) !! intensity


! parameters in this routine
  integer,parameter:: nchk1=3, nchk2=3
  integer,parameter:: klgt1=2*kndm

! initialization
  integer,save:: init0=1
  integer:: init
  integer:: nln1
  integer:: nda2,nda21,mxlgt1,mxlgn1
  integer:: indpl

! for chkrt
  real(8):: amua(kndm)
  real(8):: wa(kndm)
  real(8):: dpt(knln1)
  integer:: mxlgn2

! for work
  integer:: l,k2
  real(8):: cs(knang),yy(knang)
  real(8):: gg0,gg(klgn1),sum
  integer:: ichk,ichkt
  real(8):: ff(knln),thkt(knln)
  real(8):: dtau,dtaut
  real(8):: plgd,searf1,pint4

  integer:: mmax1
  integer:: iichk0(kna0),jm0(kna0)
  integer:: iichk1(kna1u),im1u(kna1u)

! for ftrn
  integer:: m
  real(8):: fm
  integer:: nb0
  real(8):: bm0(kna0)
  integer:: nb1u
  real(8):: bm1u(kna1u)
  integer:: nlgt1(knln)
  real(8):: gt(kndm2,knln)
  real(8):: dptt(knln1)
  real(8):: omgt(knln)
  real(8):: cplkt(kplk1,knln)
  real(8):: utaut(kntau)
  real(8):: aif(kna1u,kna0,kntau)

! for intcr1
  real(8):: cor(kntau)
  real(8):: sgl2(kntau)
  real(8):: phsb(knln)

! work  
  integer:: m1,jb,ib,k,i,j,it
  real(8):: eu0(kna0),eu1(kna1u)
  real(8):: cosm(knfi)
  real(8):: dai,eai

  integer:: is,iu,iz,l9,k1
  real(8):: cs1,ang1
  real(8):: phs(knln),phst(knln)

  real(8):: a3,am11,fis
  logical :: rt_failed  !logic if radiative tarnsfer is failed, --> true (default is fault)

!--exec
  if(init0 > 0) then
     init0=0
     if(indg > 0) indg=1
  endif

! check and set of variables
  call chkrt(inda,indt,indp,imthd,nda,na1u,am1u,na0,am0,nfi,fi,&
       nln,thk,omg,nlgn1,g,nang,ang,phsf,epsp,epsu,galb,fsol, &
       nplk1,cplk,bgnd,ntau,utau,amua,wa,dpt,mxlgn2)
  nln1=nln+1

! fourier expansion parameters
  nda2=2*nda
  nda21=nda2+1
  mxlgt1=nda2
  if(mxlgt1 > klgt1) then
     stop 'MSG simulator_visir: -klgt1- is too small'
  endif

  if(inda <= 0) then

!! flux
     if(imthd < 0) then
        mxlgn1=min(nda2,mxlgn2)
     else
        mxlgn1=min(nda21,mxlgn2)
     endif

  else

!! intensity
     if(indp <= 0) then
        mxlgn1=mxlgn2
     else
        if(imthd <= -1) then
           mxlgn1=nda2
        else
           if(imthd <= 2) then
              mxlgn1=nda21
           else
              mxlgn1=mxlgn2
           endif
        endif
     endif
  endif

  if(mxlgn1 > klgn1) then
     stop 'MSG simulator_visir: -klgn1- is too small'
  endif

! indicater for re-calculating phase function using -g-
  indpl=0
  if(indp == 0 .and. (inda >= 1 .and. imthd >= 1)) indpl=1

! clear variables
  if(inda >= 1) ai(1:na1u,1:na0,1:nfi,1:ntau)=0.d0

! set cos(scattering angle)
  if((indp == 1) .or. (indpl >= 1)) cs(1:nang)=cos(ang(1:nang)*rad)
! loop for sublayers
  do l=1,nln
!! Legendre expansion of the phase functions
!!  we truncate the series when abs(g) becomes smaller than -epsp-
!!   for successive three times.
     if(indp == 1) then
        yy(1:nang)=phsf(1:nang,l)
        call lgndf3(mxlgn1,nang,cs,yy,gg)
        gg0=gg(1)
        ichk=0
        do k1=1,mxlgn1
           g(k1,l)=gg(k1)/gg0
           if(abs(g(k1,l)) <= epsp) then
              ichk=ichk+1
              if(ichk >= nchk1) goto 6
           endif
        enddo
        k1=mxlgn1
    6   nlgn1(l)=k1

! cs is from 1 to -1 for lgndf3, so that gg0<0 (90.11.22)
        gg0=abs(gg0)
        phsf(1:nang,l)=phsf(1:nang,l)/gg0/4.d0/pi

!! construction of phase function from -g- when indp=0 and imthd >= 1
     else
        gg0=g(1,l)
        ichk=0
        k2=nlgn1(l)
        do k1=1,k2


           g(k1,l)=g(k1,l)/gg0
           if(abs(g(k1,l)) <= epsp) then
              ichk=ichk+1
              if(ichk >= nchk1) goto 46
           endif
        enddo
        k1=mxlgn1
  46    nlgn1(l)=k1
        if(indpl >= 1) then
           do i=1,nang
              cs1=cs(i)
              init=1
              sum=0
              do k1=1,nlgn1(l)
                 sum=sum+dble(2*k1-1)*g(k1,l)*plgd(init,cs1)
              enddo
              phsf(i,l)=sum/4.d0/pi
           enddo
        endif
     endif

!! delta-m moments (ichkt=1 means we do truncation at least one time)
     ichkt=0
     nlgt1(l)=min(mxlgt1,nlgn1(l))
     if(imthd >= 0 .and. nlgn1(l) > mxlgt1) then
        ff(l)=g(mxlgt1+1,l)
        ichkt=1
     else
        ff(l)=0.d0
     endif
        
     gt(1:nlgt1(l),l)=(g(1:nlgt1(l),l)-ff(l))/(1.d0-ff(l))
     omgt(l)=(1.d0-ff(l))*omg(l)/(1.d0-ff(l)*omg(l))
     thkt(l)=(1.d0-ff(l)*omg(l))*thk(l)


  enddo

! reset -mxlgn1- and -mxlgt1-
  mxlgn1=1
  mxlgt1=1
  do l=1,nln
     mxlgn1=max(mxlgn1,nlgn1(l))
     mxlgt1=max(mxlgt1,nlgt1(l))
  enddo

! set truncated user defined depths
  dptt(1)=dpt(1)
  do l=1,nln
     dptt(l+1)=dptt(l)+thkt(l)
  enddo

  do it=1,ntau
     do l=1,nln
        if((utau(it)-dpt(l))*(utau(it)-dpt(l+1)) <= 0.d0) then
           dtau=utau(it)-dpt(l)
           dtaut=(1.d0-ff(l)*omg(l))*dtau
           utaut(it)=dptt(l)+dtaut
        endif
     enddo
     if(utaut(it) > dptt(nln+1)) utaut(it)=dptt(nln+1)
  enddo

! scaling coefficients for thermal emission
  if(nplk1 > 0) then
     do k1=1,nplk1
        cplkt(k1,1:nln)=cplk(k1,1:nln)/(1.d0-omg(1:nln)*ff(1:nln))**(k1-1)
     enddo
  endif

! set -mmax1-
  mmax1=1
  if(inda >= 1) mmax1=mxlgt1

! initialize angle index for checking convergence
  nb0=na0

  iichk0(1:na0)=0
  do j=1,na0
        jm0(j)=j
  enddo

  bm0(1:na0)=am0(1:na0)
  if(inda > 0) then
     nb1u=na1u
     iichk1(1:na1u)=0
     do i=1,na1u; im1u(i)=i; enddo
     bm1u(1:na1u)=am1u(1:na1u)
  endif

! fourier sum
  init=1
  do m1=1,mmax1
     m=m1-1
     fm=pi
     if(m == 0) fm=2.d0*pi

     call ftrn21(indg,init,m,inda,indt,imthd,nb1u,bm1u,nb0,bm0, &
          nda,amua,wa,nln,nlgt1,gt,dptt,omgt,nplk1,cplkt,galb, &
          bgnd,fsol,ntau,utaut,scr,sci,flxd,flxu,aif,rt_failed)
     if( rt_failed ) then
         flxd = 0. ; flxu = 0.
         ai_sdsu = 0.  !zero out, when radiative trnasfer is failed
         return
     endif


     if(inda > 0) then
        eu0(1:nb0)=0.d0
        eu1(1:nb1u)=0.d0
        cosm(1:nfi)=cos(dble(m)*fi(1:nfi)*rad)/fm
        
        do jb=1,nb0
           j=jm0(jb)
           if(j==0) stop 'MSG simulator_visir: jm0 became 0 --> check code'
           do ib=1,nb1u
              i=im1u(ib)
              do k=1,nfi
                 do it=1,ntau
                    dai=aif(ib,jb,it)*cosm(k)
                    ai(i,j,k,it)=ai(i,j,k,it)+dai
                    if(abs(ai(i,j,k,it)) > 0) then
                       eai=abs(dai/ai(i,j,k,it))
                    else
                       if(abs(dai) <= 0) then
                          eai=0
                       else
                          eai=100
                       endif
                    endif
                    eu0(jb)=max(eu0(jb),eai)
                    eu1(ib)=max(eu1(ib),eai)
                 enddo
              enddo
           enddo
        enddo

! check convergence for am0
        call convu(nb0,kna0,bm0,jm0,eu0,iichk0,epsu,nchk2)

        if(nb0 <= 0) goto 30
        if(inda /= 2) then
           call convu(nb1u,kna1u,bm1u,im1u,eu1,iichk1,epsu,nchk2)
           if(nb1u <= 0) goto 30
        endif

     endif
  enddo

!toshii
 if(indg==1) indg=2  !Do not use for SDSU since u10 is grid dependent (toshi)
                      !also this cause segmentation fault

30 if(imthd <= 0 .or. inda <= 0 .or. ichkt <= 0) return

!--- intensity correction.
  do is=1,na0; do iu=1,na1u; do iz =1,nfi

!! cos(scattering angle)
     cs1=am0(is)*am1u(iu)+sqrt((1.d0-am0(is)**2)*(1.d0-am1u(iu)**2)) &
          *cos(fi(iz)*rad)
     if(abs(cs1) > 1.d0) cs1=sign(1.d0,cs1)
     ang1=acos(cs1)/rad
     
     do l=1,nln
        l9=l
        init=1

!! interpolation of original phase function from given phase function.
        if(indp >= 0) then
           phs(l)=pint4(init,ang1,nang,l9,ang,phsf)

!! interpolation of original phase function from -g-.
        else
           sum=0
           do k1=1,nlgn1(l)
              sum=sum+dble(2*k1-1)*g(k1,l)*plgd(init,cs1)
           enddo
           phs(l)=sum/4.d0/pi
        endif

!! interpolation of truncated phase function from -gt-.
        init=1
        sum=0
        do k1=1,nlgt1(l)
           sum=sum+dble(2*k1-1)*gt(k1,l)*plgd(init,cs1)
        enddo
        phst(l)=sum/4.d0/pi
     enddo
     call intcr1(imthd,am0(is),am1u(iu),cs1,ntau,utau,utaut, &
          nln,thk,thkt,omg,omgt,phs,phst,ff,mxlgn1,nlgn1, &
          nlgt1,g,epsp,nchk1,cor,sgl2,phsb)

     do it=1,ntau
! corrected single scattering from ocean surface
        if(indg <= 0) then
           a3=0.d0
        else
           am11=abs(am1u(iu))
           fis=fi(iz)*rad
           if(am1u(iu) < 0) then
              a3=fsol *searf1(am11,am0(is),fis,galb,scr,sci) &
                   *exp(-dptt(nln1)/am0(is)-(dptt(nln1)-utaut(it))/am11)
           else
              a3=0.d0
           endif
        endif
        ai(iu,is,iz,it)=ai(iu,is,iz,it)+cor(it)*fsol+a3
     enddo

! angular loop end
  enddo; enddo; enddo
! normal end

!
! toshii (I needed to put new output variables to stablize code for SDSU...)
!
do it=1,ntau
  ai_sdsu(it) = ai(1,1,1,it)
enddo


  return
end subroutine simulator_visir 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine cspl1(n,x,y,a,b,c,d)
! getting coefficients of natural cubic spline fitting
! use a linear interpolation outside the meaningful range of x
! x-increasing order
!--- history
! 88. 1. 4  created
! 08.10.29 modified fortran 90 version.
!--- input
! n      i     nbr of data
! x    r(n)    independent variable data
!              x-increasing order
! y    r(n)      dependent variable data
!--- output
! a    r(n)    lambda -> a   where  y=a+x*(b+x*(c+d*x))
! b    r(n)    d      -> b   i-th for range (x(i-1), x(i))
! c    r(n)    m         c   (a,b,c,d for i=1 are same as those for i=2)
! d    r(n)    m         d
!--- notes
! ref-1   p. f. davis and philip rabinowitz (1984)
!         methods pf numerical integration, second edition
!         academic press, inc., pp612.
!--
  implicit none
  integer,intent(in)::n
  real(8),intent(in)::x(n),y(n)
  real(8),intent(out)::a(n),b(n),c(n),d(n)
  integer::i
  real(8)::h1,h2,s1,s2,q1,q2,u1,u2,p1,p2,c1,c2,x1,x2,x3,x4,z1,z2,z3,z4
!--
  if(n > 4) then
! n > 4
     a(1)=0.d0; a(n)=1.d0; b(1)=0.d0; b(n)=0.d0
     h2=x(2)-x(1); s2=(y(2)-y(1))/h2
     do i=2,n-1
        h1=h2; h2=x(i+1)-x(i)
        s1=s2; s2=(y(i+1)-y(i))/h2
        a(i)=1.d0/(1.d0+h1/h2)
        b(i)=6.d0*(s2-s1)/(h2+h1)
     enddo

     q2=0.d0; u2=0.d0
     do i=1,n
        q1=q2; u1=u2
        p2=(1.d0-a(i))*q1+2.d0
        q2=-a(i)/p2; u2=(b(i)-(1.0-a(i))*u1)/p2
        a(i)=q2; b(i)=u2
     enddo
     c(n)=b(n)
     do i=n-1,1,-1
        c(i)=a(i)*c(i+1)+b(i)
     enddo

     x2=x(1); c2=c(1)/6.d0
     do i=2,n
        h2=x(i)-x(i-1); x1=x2
        x2=x(i); c1=c2; c2=c(i)/6.d0
        p1=y(i-1)/h2-c1*h2; p2=y(i  )/h2-c2*h2
        a(i)=     ( c1*x2**3-c2*x1**3)/h2+p1*x2-p2*x1
        b(i)=3.d0*(-c1*x2**2+c2*x1**2)/h2-p1   +p2
        c(i)=3.d0*( c1*x2   -c2*x1   )/h2
        d(i)=     (-c1      +c2      )/h2
     enddo

  elseif(n==1) then
! n=1
     a(1)=y(1); b(1)=0.d0; c(1)=0.d0; d(1)=0.d0
     return

  elseif(n==2) then
! n=2
     x1=x(1); x2=x(2); z1=y(1)/(x1-x2); z2=y(2)/(x2-x1)
     a(2)=-x2*z1-x1*z2; b(2)=z1+z2; c(2)=0.d0; d(2)=0.d0

  elseif(n==3) then
! n=3
     x1=x(1); x2=x(2); x3=x(3)
     z1=y(1)/(x1-x2)/(x1-x3); z2=y(2)/(x2-x3)/(x2-x1)
     z3=y(3)/(x3-x1)/(x3-x2)
     a(2)=x2*x3*z1+x3*x1*z2+x1*x2*z3
     b(2)=-(x2+x3)*z1-(x3+x1)*z2-(x1+x2)*z3
     c(2)=z1+z2+z3; d(2)=0.d0
     a(3)=a(2); b(3)=b(2); c(3)=c(2); d(3)=d(2)

  elseif(n==4) then
! n=4
     x1=x(1); x2=x(2); x3=x(3); x4=x(4)
     z1=y(1)/(x1-x2)/(x1-x3)/(x1-x4)
     z2=y(2)/(x2-x3)/(x2-x4)/(x2-x1)
     z3=y(3)/(x3-x4)/(x3-x1)/(x3-x2)
     z4=y(4)/(x4-x1)/(x4-x2)/(x4-x3)
     a(2:4)=-x2*x3*x4*z1-x3*x4*x1*z2-x4*x1*x2*z3-x1*x2*x3*z4
     b(2:4)=(x2*x3+x3*x4+x4*x2)*z1+(x3*x4+x4*x1+x1*x3)*z2 &
          +(x4*x1+x1*x2+x2*x4)*z3+(x1*x2+x2*x3+x3*x1)*z4
     c(2:4)=-(x2+x3+x4)*z1-(x3+x4+x1)*z2-(x4+x1+x2)*z3-(x1+x2+x3)*z4
     d(2:4)=z1+z2+z3+z4
  endif

  a(1)=a(2); b(1)=b(2); c(1)=c(2); d(1)=d(2)
  return
end subroutine cspl1

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function cspli(x1,n,x,a,b,c,d)
! get interpolated value using cubic spline (pair with cspl1)
!--- history
! 88. 6. 6  registered by t. nakajima
! 08.10.29 modified fortran 90 version.
!--- input
! x1      r     interpolation point
! n       i     no. of intervals + 1
! x     r(n)    division points of the intervals.
! a, b, c, d
!       r(n)    y=a+x*(b+x*(c+x*d)))
!--- output
! cspli  rf     interpolated value
!--
  implicit none
  integer,intent(in)::n
  real(8),intent(in)::x1,x(n),a(n),b(n),c(n),d(n)
  integer::j
  real(8)::cspli
!--exec
  do j=1,n
     if(x1 <= x(j)) then
        cspli=a(j)+x1*(b(j)+x1*(c(j)+x1*d(j)))
        return
     endif
  enddo
  cspli=a(n)+x1*(b(n)+x1*(c(n)+x1*d(n)))
  return
end function cspli

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function erfc(x)
! complementary error function=(2/sqrt(pi) integral(x,inf) exp(-t**2)dt
! 92.12.22: created by Nakajima
!--- input
! x       r      independent variable 0 to inf
!--- output
! erfc   rf      complementary error function
!--
  implicit none

! input & output
  real(8),intent(in)::x
  real(8):: erfc

! work
  real(8),parameter:: c1 = 7.05230784d-2
  real(8),parameter:: c2 = 4.22820123d-2
  real(8),parameter:: c3 = 9.2705272d-3
  real(8),parameter:: c4 = 1.520143d-4
  real(8),parameter:: c5 = 2.765672d-4
  real(8),parameter:: c6 = 4.30638d-5

  real(8)::v
!--
  v=abs(x)
  if(v <= 7.d0) then
     erfc=1.d0/(1.d0+v*(c1+v*(c2+v*(c3+v*(c4+v*(c5+c6*v))))))**16
  else
     erfc=0.d0
  endif
  if(x < 0.d0) erfc=2.d0-erfc
  return
end function erfc

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function expfn(x)
! exponential function with over/under flow setting.
!--- history
! 89. 5. 4   created by T. Nakajima.
! 90. 1.17   updated with cpcon
! 03. 9. 9   changed fpr free source term
!--- input
! x        R         independent variable.
!--- output
! expfn    F         exp(x).
!                     if x.le. vmn then exp(x) is reset as   0.
!                     if x.ge. vmx then exp(x) is reset as exp(vmx).
!--- parameters
! system set the -vmn- and -vmx- by the function r1mach.
!--
  implicit none
  integer,save::init
  real(8),save::vmn,vmx,expmn,expmx
  real(8),intent(in)::x
  real(8)::ccp(3),expfn
!--exec
  init=1

! set vmn      r         enderflow limit.
!     vmx      r         overflow limt.
  if(init > 0) then
     init=0
     call cpcon(ccp)
     vmn=ccp(2)*0.8*2.3; vmx=ccp(3)*0.8*2.3
     expmn=0.0d0; expmx=exp(vmx)
  endif

  if(x <= vmn) then
     expfn=expmn
  else
     if(x >= vmx) then
        expfn=expmx
     else
        expfn=exp(x)
     endif
  endif
  return
end function expfn

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine cpcon(c)
! machine constants of computer
!--- output
! c     r(3)      (1)  minimum positive x for  1+x      .ne. 1
!                      average as a result of complex arithmetic
!                      operations.
!                 (2)  minimum exponent y for  10.0**y  .ne. 0
!                 (3)  maximum exponent z for  10.0**z  is max. value
!                  if init=1 (data statement) then  set as z=y
!                  if init=2 then this routine gets actual value of z.
!                  - see note -
!--- history
! 90. 1.20  created
!     6.27  change the algorithm to get x and y taking into account the
!           high accuracy co-processor and graceful underflow.
! 92. 3.21  n=1000 from n=200
!     7. 9  bug in x-definition
! 03. 9. 9   changed fpr free source term
!--- note
! this program will generate -underflow error- and -overflow error-
!  messages.  on some computer -overflow error- message may be
!  fatal error.  in that case, please set init = 1 in the data
!  satement for suppressing the procedure of getting c(3).
!--
  use module_opt_visir, only: pi
  implicit none
  real(8),intent(out)::c(3)
  real(8),save::x,y,z
  integer,save::init=1
  real(8)::y1,y2,y3,z1,z2
  integer::n,m0,m,k,i
  !character ch*80
!--exec
! resolution of computation
  if(init <= 0) then
     c(1)=x; c(2)=y; c(3)=z
     return
  endif

!! test sum(k=1,m) cos((2k-1)*pi/(2m+1)) = 0.5
!! simple check, x = x + e, is not valide when the computer
!!  use a high accurate co-processor.
  n=500; m0=10; x=0.d0
  do m=1,m0
     y=0.d0
     do k=1,m
        y=y+cos((2.d0*k-1)*pi/(2.d0*m+1))
     enddo
     y=abs(2.d0*y-1)
     x=x+y
  enddo
  x=x/m0
  c(1)=x

! exponent for minimum positive value
! this procedure will generate -underflow error message-
  y2=1.d0; n=1000
  do i=1,n
     y1=y2; y3=y1/10.d0
!! for graceful underflow
!! even y2 becomes 0 as output, y2 is not 0 inside
!! computer when graceful underflow is applied.
!! so we replace the value of y2 by output.
!         ch='0'
!         write(ch,107) y3
!  107    format(1p,e12.5)
!         y2=0
!         read(ch,*) y2
     y2=y3
     if(abs(10.d0*y2/y1-1.d0) > 5.d0*x) goto 104
  enddo
  i=n+1
104 y=1-i
  c(2)=y

! exponent for maximum positive value
! this procedure will generate -overflow message-
  if(init <= 1) then
     z=-y
  else
     z2=1
     do i=1,n
        z1=z2; z2=z1*10.d0
        if(abs(z2/z1/10.d0-1.d0) > 5.d0*x) goto 106
     enddo
     i=n+1
106  z=i-1
  endif
  c(3)=z

  init=0
  return
end subroutine cpcon

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function plgd(init,x)
! Legendre polynomials
!--- history
! 87.11.12 created
! 90. 1.16 save statement
!--- input
! init   i    if 1 then l=0
!             if 0 then l=l+1
! x      r    (-1,1)
!--- out
! init=0
  implicit none
  integer,intent(inout)::init
  real(8),intent(in)::x
  integer,save::l
  real(8),save::pl,pl1,pl2
  real(8)::plgd
!--
  if(init > 0) then
     init=0; l=-1
  endif
  l=l+1
  if(l == 0) then
     pl=1
  elseif(l == 1) then
     pl1=pl; pl=x
  else
     pl2=pl1; pl1=pl
     pl=(dble(2*l-1)*x*pl1-dble(l-1)*pl2)/dble(l)
  endif
  plgd=pl
  return
end function plgd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine plgnd(m1,mmx1,nx,nx0,x,pl)
! normalized associated legendre polynomials
!------ history
! 1987.03.04
! 08.10.24 modified fortran 90 version.
!------ input variables
! variable  type    interpretation
! m1         i      order of the fourier series + 1
! mmx1       i      max order of m  + 1
! nx         i      nbr of x
! nx0        i      declared nbr of nx
! x       r(nx)     xi, i=1,nx
!------ output variables
! pl      r(nx0,    normalized associated legendre polynomials
!           mmx1)
  use module_opt_visir, only:pi
  implicit none
  integer,intent(in)::m1,mmx1,nx,nx0
  real(8),intent(in)::x(nx)
  real(8),intent(out)::pl(nx0,mmx1)
  real(8)::eta,epsi,expcc,c0,c1,c2
  integer::m,l1,l,k,j,i
!--exec
! normalized associated legendre polynomials.
!  cc=p(m,l,x)*sqrt((l-m)|/(l+m)|), where x= cos(mu).
  m=m1-1
  do l1=m1,mmx1
     l=l1-1
     if(l < (m+1)) then
        if(m <= 0) then
           pl(1:nx,l1)=1.d0
        else
           k=2; eta=1.d0
           do j=1,m
              epsi=1.d0-1.d0/dble(k)
              eta=eta*epsi
              k=k+2
           enddo
           eta=sqrt(eta)
           do i=1,nx
              if(x(i) > 1.d0) then
                 pl(i,l1)=0.d0
              else
                 expcc=(dble(m)*log(  max(1.d-15,abs(1.d0-x(i)**2)) ))/2.d0  !add max to avoid log(zero) toshii
                 if(expcc <= -100.d0) then
                    pl(i,l1)=0.d0
                 else
                    pl(i,l1)=eta*exp(expcc)
                 endif
              endif
           enddo
        endif
     elseif(l == m+1) then
        c0=sqrt(dble(2*m+1))
        pl(1:nx,l1)=c0*pl(1:nx,l1-1)*x(1:nx)
        c0=sqrt(dble((l-m)*(l+m)))
     else
        c1=sqrt(dble((l-m)*(l+m)))
        c2=dble(2*l-1)
        pl(1:nx,l1)=(c2*x(1:nx)*pl(1:nx,l1-1)-c0*pl(1:nx,l1-2))/c1
        c0=c1
     endif
  enddo
  return
end subroutine plgnd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine  qgausn( gwt, gmu, m )
!  compute weights and abscissae for ordinary gaussian quadrature
!   (no weight function inside integral) on the interval (0,1)
!--- history
! 90. 1.17  registered
! 08. 9.25 modified Fortran 90 free-style form
!--- input
! m        i       order of quadrature rule
!--- output
! gmu    r(m)      array of abscissae (0, 1)
! gwt    r(m)      array of weights   sum=1
!--- notes
! reference:  Davis,p.j. and p. Rabinowitz, methods of numerical
!             integration,academic press, new york, 1975, pp. 87
! method:     compute the abscissae as roots of the legendre
!             polynomial p-sub-n using a cubically convergent
!             refinement of newton's method.  compute the
!             weights from eq. 2.7.3.8 of davis/rabinowitz.
!             accuracy:  at least 13 significant digits.
!--- internal variables
! pm2,pm1,p : 3 successive legendre polynomials
! ppr       : derivative of legendre polynomial
! p2pri     : 2nd derivative of legendre polynomial
! tol       : convergence criterion
! x,xi      : successive iterates in cubically-
!             convergent version of newton's method
!            ( seeking roots of legendre polynomial )
!--
  use module_opt_visir, only:pi
  implicit none

! input & output
  integer,intent(inout)::m
  real(8),intent(out)::gmu(m),gwt(m)

! work
  integer::lim,np1,k,nn
  real(8)::cona, t, en, nnp1, p, pm1, pm2, ppr, p2pri, prod, tmp, x, xi
  real(8)::tol=1.0d-13
!--exec
  if ( m <= 1 )  then
     m = 1
     gmu( 1 ) = 0.5d0
     gwt( 1 ) = 1.0d0
     return
  end if

  en   = m
  np1  = m + 1
  nnp1 = m * np1
  cona = dble( m-1 ) / ( 8.d0 * m**3 )
!+---------------------------------------------------------------------+
!!         initial guess for k-th root of legendre polynomial,         |
!!         from davis/rabinowitz  eq. (2.7.3.3a)                       |
!+---------------------------------------------------------------------+
  lim  = m / 2
  do k = 1, lim
     t = dble( 4*k - 1 ) * pi / dble( 4*m + 2 )
     x = cos ( t + cona / tan( t ) )

!+---------------------------------------------------------------------+
!!             recursion relation for legendre polynomials             |
!!       initialize legendre polynomials: (pm2) p-sub-0, (pm1) p-sub-1 |
!+---------------------------------------------------------------------+
10   pm2 = 1.d0
     pm1 = x
     do nn = 2, m
!            p   = ( ( 2*nn - 1 ) * x * pm1 - ( nn-1 ) * pm2 ) / nn
        p   =  2* x * pm1 - pm2 - ( x * pm1 - pm2 ) / nn
        pm2 = pm1
        pm1 = p
     enddo

     tmp   = 1.d0 / ( 1.d0 - x**2 )
     ppr   = en * ( pm2 - x * p ) * tmp
     p2pri = ( 2.d0 * x * ppr - nnp1 * p ) * tmp
     xi    = x - ( p / ppr ) * ( 1.d0 + ( p / ppr ) &
          * p2pri / ( 2.d0 * ppr ) )

     if ( dabs(xi-x) > tol ) then
!!          check for convergence
        x = xi
        go to 10
     end if

!       ** iteration finished--calc. weights, abscissae for (-1,1)
     gmu( k ) = - x
     gwt( k ) = 2.d0 / ( tmp * ( en * pm2 )**2 )
     gmu( np1 - k ) = - gmu( k )
     gwt( np1 - k ) =   gwt( k )
  enddo

  if ( mod( m,2 ) /= 0 )  then
!!       set middle abscissa and weight for rules of odd order
     gmu( lim + 1 ) = 0.d0
     prod = 1.d0
     do k = 3, m, 2
        prod = prod * k / ( k-1 )
     enddo
     gwt( lim + 1 ) = 2.d0 / prod**2
  end if

  do k = 1, m
!!       convert from (-1,1) to (0,1)
     gmu( k ) = 0.5d0 * gmu( k ) + 0.5d0
     gwt( k ) = 0.5d0 * gwt( k )
  enddo
  return
end subroutine qgausn

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine symtrx(rowij,m,root,eigv,ni,wk,ier)
! solves eigenfunction problem for symmetric matrix
!--- history
! 89.12. 4 modified with cncpu
! 90. 1. 6 cncpu is replaced by prec.
! 08.10.29 modified fortran 90 version.
!--- input
! m       i      order of original symmetric matrix
! ni      i      initial dimension of -root-, -eigv- and -wk-
! rowij  r(*)    symmetric storage mode of order m*(m+1)/2
!--- output
! eigv   r(ni,m) eigenvectors of original symmetric matrix
! ier     i      index for root(j) failed to converge (j=ier-128)
! root   r(m)    eigenvalues of original symmetric matrix
! rowij          storage of householder reduction elements
! wk             work area
!--      
  use module_opt_visir,only: krowij
  implicit none

! input & output
  integer,intent(in):: m
  integer,intent(in):: ni
  real(8),intent(inout):: rowij(krowij)
  real(8),intent(out):: root(ni)
  real(8),intent(out):: eigv(ni,ni)
  real(8),intent(out):: wk(ni)
  integer,intent(out):: ier

! work
  real(8)::ccp(3),rdelp,h,scale,f,g,hh,b,sum1,c,s,p,r
  integer::mp1,mm,mbeg,i,l,ii,mk,k,jk1,ik,jk,jp1,ia,nn1,j,n,l1,nn1pl
!--exec
!+++ add epscp
  call cpcon(ccp)
  rdelp=ccp(1)*10.d0
!+++
  ier = 0; mp1 = m + 1
  mm = (m*mp1)/2 - 1
  mbeg = mm + 1- m

!+---------------------------------------------------------------------+
!|          loop-100 reduce -rowij- (symmetric storage mode) to a      |
!|          symmetric tridiagonal form by householder method           |
!|                      cf. wilkinson, j.h., 1968,                     |
!|              the algebraic eigenvalue problem, pp 290-293.          |
!|          loop-30&40 and 50 form element of a*u and element p        |
!+---------------------------------------------------------------------+
! (loop-100)
  do ii=1,m
     i = mp1 - ii
     l = i - 1
     h = 0.d0; scale = 0.d0
     mk = mm
     do k=1,l
        scale = scale + abs(rowij(mk))
        mk = mk - 1
     enddo

     if (l < 1 .or. scale <= 0.0) then
! in l case   scale row (algol tol then not needed)
        wk(i) = 0.0
     else
        mk = mm
        do k = 1,l
           rowij(mk) = rowij(mk)/scale
           h = h + rowij(mk)*rowij(mk)
           mk = mk - 1
        enddo
        wk(i) = scale*scale*h
        f = rowij(mm)
        g = - sign(sqrt(h),f)
        wk(i) = scale*g
        h = h - f*g
        rowij(mm) = f - g
        if (l > 1) then
           f = 0.d0; jk1 = 1
! (loop-50)
           do j=1,l
              g = 0.d0
              ik = mbeg + 1; jk = jk1
! (loop-30)               
              do k=1,j
                 g = g + rowij(jk)*rowij(ik)
                 jk = jk + 1; ik = ik + 1
              enddo
              jp1 = j + 1
              if (l >= jp1) then
                 jk = jk + j - 1
! (loop-40)
                 do k=jp1,l
                    g = g + rowij(jk)*rowij(ik)
                    jk = jk + k; ik = ik + 1
                 enddo
              endif
              wk(j) = g/h
              f = f + wk(j)*rowij(mbeg+j)
              jk1 = jk1 + j
           enddo
           hh = f/(h+h)

           jk = 1
           do j=1,l
              f = rowij(mbeg+j)
              g = wk(j) - hh*f
              wk(j) = g
              do k=1,j
                 rowij(jk) = rowij(jk) - f*wk(k) - g*rowij(mbeg+k)
                 jk = jk + 1
              enddo
           enddo
        endif

        do k=1,l
           rowij(mbeg+k) = scale*rowij(mbeg+k)
        enddo
     end if
     root(i) = rowij(mbeg+i)
     rowij(mbeg+i) = h*scale*scale
     mbeg = mbeg - i + 1
     mm = mm - i
  enddo

!+---------------------------------------------------------------------+
!|          loop-210 compute eigenvalues and eigenvectors              |
!|          setup work area location eigv to the identity matrix       |
!|          loop-140 for finding small sub-diagonal element            |
!|          loop-160 for convergence of eigenvalue j (max. 30 times)   |
!|          loop-190 for ql transformation and loop-180 form vectors   |
!+---------------------------------------------------------------------+
  do i=1,m-1
     wk(i) = wk(i+1)
  enddo
  wk(m) = 0.d0; b = 0.d0; f = 0.d0
  eigv(1:m,1:m) = 0.d0
  do i=1,m
     eigv(i,i) = 1.d0
  enddo

! (loop-210)
  do l=1,m
     h = rdelp*(abs(root(l))+abs(wk(l)))
     if (b < h) b = h
! (loop-140)
     do n=l,m
        k = n
        if (abs(wk(k)) <= b) exit
     enddo
     n = k

     if (n/=l)  then
        do j=1,30
           if (j==31) then
              ier = 128 + l
              return
           end if

           l1 = l + 1
           g = root(l); p = (root(l1)-g)/(wk(l)+wk(l))
           r = abs(p)
           if (rdelp*abs(p).lt.1.0) r = sqrt(p*p+1.0)
           root(l) = wk(l)/(p+sign(r,p))
           h = g - root(l)
           root(l1:m) = root(l1:m) - h
           f = f + h

           p = root(n); c = 1.d0; s = 0.d0
           nn1 = n - 1; nn1pl = nn1 + l
           if (l <= nn1) then
              do ii=l,nn1
                 i = nn1pl - ii
                 g = c*wk(i)
                 h = c*p
                 if (abs(p) < abs(wk(i))) then
                    c = p/wk(i)
                    r = sqrt(c*c+1.d0)
                    wk(i+1) = s*wk(i)*r
                    s = 1.d0/r
                    c = c*s
                 else
                    c = wk(i)/p
                    r = sqrt(c*c+1.d0)
                    wk(i+1) = s*p*r
                    s = c/r
                    c = 1.d0/r
                 end if
                 p = c*root(i) - s*g
                 root(i+1) = h + s*(c*g+s*root(i))
                 if (ni >= m) then
                    do k=1,m
                       h = eigv(k,i+1)
                       eigv(k,i+1) = s*eigv(k,i) + c*h
                       eigv(k,i) = c*eigv(k,i) - s*h
                    enddo
                 end if
              enddo
           end if
           wk(l) = s*p; root(l) = c*p
           if (abs(wk(l)) <= b) exit
        enddo
     endif
     root(l) = root(l) + f
  enddo

!+---------------------------------------------------------------------+
!|          back transform eigenvectors of the original matrix from    |
!|          eigenvectors 1 to m of the symmetric tridiagonal matrix    |
!+---------------------------------------------------------------------+
  do i=2,m
     l = i - 1
     ia = (i*l)/2
     if (abs(rowij(ia+i)) > 0.0) then
        do j=1,m
           sum1=sum(rowij(ia+1:ia+l)*eigv(1:l,j))/rowij(ia+i)
           do k=1,l
              eigv(k,j) = eigv(k,j) - sum1*rowij(ia+k)
           enddo
        enddo
     end if
  enddo

  return
end subroutine symtrx


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function phintg(nang,ang,ph)
!---input
! nang    i      nbr of scattering angle
! ang  r(nang)   independent variables (scat. ang.)
! ph   r(nang)   dependent variables
!---output
! phintg  fr     integral
!--
  use module_opt_visir, only:knang,pi,rad,pstd
  implicit none

! input & output
  integer,intent(in):: nang  
  real(8),intent(in):: ang(knang)
  real(8),intent(in):: ph(knang)
  real(8):: phintg

! work
  integer:: ia
  real(8):: dx
!--
  phintg=0.d0
  do ia=nang,2,-1
     dx=cos(ang(ia-1)*rad)-cos(ang(ia)*rad)
     phintg=phintg+dx*(ph(ia-1)+ph(ia))*0.5d0
  enddo
  phintg=phintg*2.d0*pi
  return
end function phintg

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 function rawb(rh1,r,rho,naw,aw,rmmd)
! growth of wet aerosols
! raw is r(aw) in eq.(2) in shettle and fenn.
!--- history
! 94.12.26 created by i. lensky
! 95. 9.13 modified by t. nakajima
! 96. 5. 5 drop caw, number of iterations 10 -> 6
! 08. 9.19 modified Fortran 90 free-style form
!--- input
! rh1    r           relative humidity 0 <= rh1 < 1
! r      r           dry aerosol radius (cm)
! rho    r           paricle density relative to water
! naw    i           number of aw (if 0 then no growth)
! aw     r(kaw)      water activity
! rmmd   r(kaw)      Hanel's water uptake data
!--- output
! raw    r           r(aw) mode radius of wet aerosols (cm)
!--
  use module_opt_visir, only: kaw
  implicit none

! input & output
  real(8),intent(in):: rh1       !! relative humidity
  real(8),intent(in):: r         !! dry aerosol radius [cm]
  real(8),intent(in):: rho       !! particle density relative to water
  integer,intent(in):: naw       !! number of aw
  real(8),intent(in):: aw(kaw)   !! water activity
  real(8),intent(in):: rmmd(kaw) !! Hanel's water uptake data
  real(8):: rawb                 !! mode radius of wet aerosols [cm]

! work
  integer:: it,iaw
  real(8):: aw1,caw,r1
!--exec
  if (naw < 0 .or. rh1 < aw(1)) then
     rawb = r
     return
  endif
! iteration
  aw1 = rh1
  do it=1,6
! eq.(5)  of shettle & fenn
     if(aw1 <= aw(1)) then
        iaw=1
     else if(aw1 >= aw(naw)) then
        iaw=naw
     else
        do iaw=1,naw-1
           if(aw1 >= aw(iaw) .and. aw1 < aw(iaw+1)) exit
        enddo
     endif
     caw=log(rmmd(iaw+1)/rmmd(iaw))/log((1-aw(iaw+1))/(1-aw(iaw)))
     r1=rmmd(iaw)*((1-aw1)/(1-aw(iaw)))**caw
     rawb=r*(1.d0+rho*r1)**(1.d0/3.d0)
! aw1 is aw, eq.(4) in shettle & fenn. rawb in cm.
     aw1=rh1*exp(-0.001056d0/(1.d4*rawb))
  enddo
  return
 end function rawb

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine frnlr(cr,ci,wi,wr,r1,r2)
! fresnel reflection coefficients
!--- history
! 88. 6.16  created
!--- input
! cr      r      real part of the complex refraction index
! ci      r      imaginary part m = cr + i*ci
! wi      r      incident angle in degrees
!--- output
! wr      r      refraction angle in degrees
!                if total reflection, then wr=999
! r1      r      reflectivity for polarization
! r2      r      reflectivity for polarization
!--- reference
! k.-n. liou
!--
  use module_opt_visir, only: pi,rad
  implicit none

! input 
  real(8),intent(in):: cr
  real(8),intent(in):: ci
  real(8),intent(in):: wi

! output
  real(8),intent(out):: wr
  real(8),intent(out):: r1
  real(8),intent(out):: r2

! work
  real(8):: ci1,c,s
  real(8):: en1,en2,en3,u,v,u2
!--  
  ci1=-ci
  c=cos(wi*rad); s=sin(wi*rad)
  wr=s/cr
  if(wr.gt.1.d0) then
     wr=999.d0
  else
     wr=sin(wr)/rad
  endif

  en1=cr**2-ci1**2; en2=2.d0*cr*ci1; en3=en1-s**2
  u=sqrt(abs(en3+sqrt(en3**2+en2**2))/2.d0)
  v=sqrt(abs(u**2-en3))
  u2=u**2+v**2

  r1=1.d0/(u*c/(u2+c**2)+0.5d0)-1.d0
  r2=1.d0/((en1*u+en2*v)*c/((en1**2+en2**2)*c**2+u2)+0.5d0)-1.d0
  return
end subroutine frnlr

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine add(it,na0,n1,n2,n3,r1d,r1u,t1d,t1u,s1d,s1u,r2d,t2d,t2u, &
     s2d,s2u,tu,rd,sd,su)
! adding two layers 1 and 2.
!--- history
! 86.10.15  check ok.
! 89.10.30  re-edit.
! 90. 2. 2  eleminate knd00
! 08.10.27 modified fortran 90 version.
!--- input
! it         i      indicator for calculation of tu and sd.
! na0        i      no. of solar directions.
! n1,n2,n3   i       ----------------------------     mu1(i), i=1,n1
! r1d, r1u r(kndm,   r1d, r1u, t1d, t1u, s1d, s1u    (layer-1)
! t1d, t1u   kndm)   ----------------------------     mu2(i), i=1,n2
! r2d, r2u           r2d, r2u, t2d, t2u, s2d, s2u    (layer-2)
! t2d, t2u           ----------------------------     mu3(i),i=1,n3
! s1d, s1u r(kndm,   suffix u = upgoing,   d = downgoing incidences.
! s2d, s2u   kna0)          r = reflection,t = transmission matrices.
!                           s = source matrix.
!--- output
! rd       r(kndm,   -----------------    mu1(i), i=1,n1
! tu         kndm)   rd, tu, sd, su       (layer 1+2)
! sd       r(kndm,   -----------------    mu3(i), i=1,n3
! su         kna0)
!--- parameter
! kna0       i       number of solar zenith angles.
! kndm       i       number of quadrature points.
!--- areas for this routine

  use module_opt_visir,only:kna0,kndm
  implicit none

! input & output  
  integer,intent(in):: it
  integer,intent(in):: na0
  integer,intent(in):: n1
  integer,intent(in):: n2
  integer,intent(in):: n3
  real(8),intent(in):: r1d(kndm,kndm)
  real(8),intent(in):: r1u(kndm,kndm)
  real(8),intent(in):: t1d(kndm,kndm)
  real(8),intent(in):: t1u(kndm,kndm)
  real(8),intent(in):: s1d(kndm,kna0)
  real(8),intent(in):: s1u(kndm,kna0)
  real(8),intent(in):: r2d(kndm,kndm)
  real(8),intent(in):: t2d(kndm,kndm)
  real(8),intent(in):: t2u(kndm,kndm)
  real(8),intent(in):: s2d(kndm,kna0)
  real(8),intent(in):: s2u(kndm,kna0)

  real(8),intent(out):: rd(kndm,kndm)
  real(8),intent(out):: tu(kndm,kndm)
  real(8),intent(out):: sd(kndm,kna0)
  real(8),intent(out):: su(kndm,kna0)

!--- working areas
  real(8)::aa(kndm,kna0),bb(kndm,kndm),cc(kndm,kndm),dd(kndm,kndm)
!--exec
  aa(1:n2,1:na0)=matmul(r2d(1:n2,1:n2),s1d(1:n2,1:na0))
  su(1:n2,1:na0)=aa(1:n2,1:na0)+s2u(1:n2,1:na0)
  cc(1:n2,1:n2)=matmul(r2d(1:n2,1:n2),r1u(1:n2,1:n2))

  call multi(n2,cc,bb)
!  if(err /= '') then
!     err='error in multi of add'; return
!  endif

  sd(1:n2,1:na0)=matmul(bb(1:n2,1:n2),su(1:n2,1:na0))
  su(1:n1,1:na0)=matmul(t1u(1:n1,1:n2),sd(1:n2,1:na0))
  su(1:n1,1:na0)=su(1:n1,1:na0)+s1u(1:n1,1:na0)
  cc(1:n1,1:n2)=matmul(t1u(1:n1,1:n2),bb(1:n2,1:n2))
  dd(1:n2,1:n1)=matmul(r2d(1:n2,1:n2),t1d(1:n2,1:n1))
  bb(1:n1,1:n1)=matmul(cc(1:n1,1:n2),dd(1:n2,1:n1))
  rd(1:n1,1:n1)=r1d(1:n1,1:n1)+bb(1:n1,1:n1)
  if(it <= 0) return

  tu(1:n1,1:n3)=matmul(cc(1:n1,1:n2),t2u(1:n2,1:n3))
  aa(1:n2,1:na0)=matmul(r1u(1:n2,1:n2),sd(1:n2,1:na0))
  aa(1:n2,1:na0)=aa(1:n2,1:na0)+s1d(1:n2,1:na0)
  sd(1:n3,1:na0)=matmul(t2d(1:n3,1:n2),aa(1:n2,1:na0))
  sd(1:n3,1:na0)=sd(1:n3,1:na0)+s2d(1:n3,1:na0)
  return
end subroutine add

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine adisc(m,l,nda,na0,am0,wmp,tc,t1,zee,qe,qie,vpe,vme,alfa,beta, &
     nplk1,dpe,dme,udn,uup)
! intensity at a user defined depth.
!--- history
! 87. 3. 9
! 89.11. 2 modified
! 08.10.27 modified fortran 90 version.
!--- input
! m        i          fourier order
! l        i          layer number.
! nda      i          number of quadrature points.
! na0      i          number of solar zenith angles.
! am0    r(kna0)      cos(solar zenith angles).
! wmp    r(kndm)      sqrt(w*m)
! tc       r          optical thickness of the layer.
! t1       r          opitcal depth of interpolation measured from top.
! zee    r(kndm,knln)  eigenvalues
! qe    r(kndm,kndm,knln)
! qie   r(kndm,kndm,knln)
! vpe   r(kndm,kna0,knln)
! vme   r(kndm,kna0,knln)
! dpe   r(kndm,kplk1,knln)
! dme   r(kndm,kplk1,knln)
! alfa  r(kndm,kna0)  integral constant alfa.
! beta  r(kndm,kna0)  integral constatn beta.
!--- output
! uup   r(kndm,kna0)  upward intensity (unscale)
! udn   r(kndm,kna0)  dnward intensity (unscale)
!--
  use module_opt_visir, only:knln,kna0,kndm,kplk1
  implicit none
 
! input
  integer,intent(in):: m
  integer,intent(in):: l
  integer,intent(in):: nda
  integer,intent(in):: na0
  real(8),intent(in):: am0(na0)
  real(8),intent(in):: wmp(kndm)
  real(8),intent(in):: tc
  real(8),intent(in):: t1
  real(8),intent(in):: zee(kndm,knln)
  real(8),intent(in):: qe(kndm,kndm,knln)
  real(8),intent(in):: qie(kndm,kndm,knln)
  real(8),intent(in):: vpe(kndm,kna0,knln)
  real(8),intent(in):: vme(kndm,kna0,knln)
  integer,intent(in):: nplk1
  real(8),intent(in):: dpe(kndm,kplk1,knln)
  real(8),intent(in):: dme(kndm,kplk1,knln)
  real(8),intent(in):: alfa(kndm,kna0)
  real(8),intent(in):: beta(kndm,kna0)

! output
  real(8),intent(out):: uup(kndm,kna0)
  real(8),intent(out):: udn(kndm,kna0)
! work areas
  integer::i,j
  real(8)::c(kndm),sl(kndm),sl1(kndm),ex1,e0(kna0),ap(kndm),am(kndm)
  real(8):: bp(kndm),bm(kndm),dp,dm,taun,sum1,sum2,expfn
!--exec
! a+-, b+-
  do i=1,nda
     call csfn(tc,t1,zee(i,l),c(i),sl(i),sl1(i))
  enddo
  do j=1,na0
     ex1=-t1/am0(j)
     e0(j)=expfn(ex1)
  enddo

! u+ = s1d,  u- = s1u
  do i=1,nda
     ap(1:nda)=qe(i,1:nda,l)*c  (1:nda) - qie(1:nda,i,l)*sl(1:nda)
     am(1:nda)=qe(i,1:nda,l)*c  (1:nda) + qie(1:nda,i,l)*sl(1:nda)
     bp(1:nda)=qe(i,1:nda,l)*sl1(1:nda) - qie(1:nda,i,l)*c (1:nda)
     bm(1:nda)=qe(i,1:nda,l)*sl1(1:nda) + qie(1:nda,i,l)*c (1:nda)
     dp=0.d0; dm=0.d0
     if(m==0 .and. nplk1 > 0) then
        taun=1.d0
        do j=1,nplk1
           dp=dp+dpe(i,j,l)*taun; dm=dm+dme(i,j,l)*taun
           taun=taun*t1
        enddo
     endif
     do j=1,na0
        sum1=sum(ap(1:nda)*alfa(1:nda,j)+bp(1:nda)*beta(1:nda,j))
        sum2=sum(am(1:nda)*alfa(1:nda,j)+bm(1:nda)*beta(1:nda,j))
        udn(i,j)=(sum1+vpe(i,j,l)*e0(j)+dp)/wmp(i)
        uup(i,j)=(sum2+vme(i,j,l)*e0(j)+dm)/wmp(i)
     enddo
  enddo
  return
end subroutine adisc

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine ainteg(m,l,nda,na0,i,fsol,am1u,am0,w,t,t1,t2,zee,qe,qie,vpe,vme, &
     pt10,pr10,pt1,pr1,alfa,beta,nplk1,dpe,dme,cplk,ai)
! 08.10.27 modified fortran 90 version.
!--- input
! m       i         fourier order
! l       i         layer number
! nda
! na0
! fsol
! i       i         stream number for am1(i)
! am1u    r         user defined direction
! am0   r(kna0)
! w       r         single scattering albedo
! t       r         optical depth for interpolation
! t1      r         optical depth at the top of the sublayer
! t2      r         optical depth at the bottom of the sublayer.
!--
  use module_opt_visir, only:kndm,kna1u,kna0,knln,kplk1,kna1,pi
  implicit none

! input & output
  integer,intent(in):: m
  integer,intent(in):: l
  integer,intent(in):: nda
  integer,intent(in):: na0
  integer,intent(in):: i
  real(8),intent(in):: fsol
  real(8),intent(in):: am1u
  real(8),intent(in):: am0(kna0)
  real(8),intent(in):: w
  real(8),intent(in):: t
  real(8),intent(in):: t1
  real(8),intent(in):: t2
  real(8),intent(in):: zee(kndm,knln)
  real(8),intent(in):: qe(kndm,kndm,knln)
  real(8),intent(in):: qie(kndm,kndm,knln)
  real(8),intent(in):: vpe(kndm,kna0,knln)
  real(8),intent(in):: vme(kndm,kna0,knln)
  real(8),intent(in):: pt10(kna1,kna0)
  real(8),intent(in):: pr10(kna1,kna0)
  real(8),intent(in):: pt1(kna1,kndm)
  real(8),intent(in):: pr1(kna1,kndm)
  real(8),intent(in):: alfa(kndm,kna0)
  real(8),intent(in):: beta(kndm,kna0)

  integer,intent(in):: nplk1
  real(8),intent(in):: dpe(kndm,kplk1,knln)
  real(8),intent(in):: dme(kndm,kplk1,knln)
  real(8),intent(in):: cplk(kplk1,knln)
  real(8),intent(out):: ai(kna0)

! working areas
  integer:: j
  real(8):: tt1,tt2,tt,am1,sum1,sum2,h1(kndm),h2(kndm),c(kndm)
  real(8):: sl(kndm),sl1(kndm),ei(kna0),pki(kplk1),pdu,trns0,hu,pvu
  real(8):: expfn
!--exec
  tt1=0.d0; tt2=t2-t1
  tt =t -t1; am1=abs(am1u)
  do j=1,nda
     sum1=sum((pt1(i,1:nda)+pr1(i,1:nda))*qe (1:nda,j,l))
     sum2=sum((pt1(i,1:nda)-pr1(i,1:nda))*qie(j,1:nda,l))
     h1(j)=w*sum1; h2(j)=w*sum2
  enddo

! downward intensity
  if(am1u > 0.d0) then
     do j=1,na0
        call exint(am1u,am0(j),tt,tt1,tt,ei(j))
     enddo
     do j=1,nda
        call csint(am1u,zee(j,l),tt2,tt,tt1,tt,c(j),sl(j),sl1(j))
     enddo
     pdu=0.d0
     if(m==0 .and. nplk1 > 0) then
        call pkint(am1u,tt,tt1,tt,nplk1,pki)
        do j=1,nplk1
           sum1=sum(pt1(i,1:nda)*dpe(1:nda,j,l)+pr1(i,1:nda)*dme(1:nda,j,l))
           pdu=pdu+(w*sum1+2.d0*pi*(1.d0-w)*cplk(j,l))*pki(j)
        enddo
     endif
     do j=1,na0
        trns0=-t1/am0(j)
        trns0=expfn(trns0)*fsol
        pvu=pt10(i,j)*trns0
        hu=sum((h1(1:nda)*c(1:nda)-h2(1:nda)*sl(1:nda))*alfa(1:nda,j) &
             +(h1(1:nda)*sl1(1:nda)-h2(1:nda)*c(1:nda))*beta(1:nda,j))
        pvu=pvu+sum(pt1(i,1:nda)*vpe(1:nda,j,l)+pr1(i,1:nda)*vme(1:nda,j,l))
        pvu=w*pvu*ei (j)
        ai(j)=(hu+pvu+pdu)/am1
     enddo
  else

! upward intensity
     do j=1,na0
        call exint(am1u,am0(j),tt,tt,tt2,ei(j))
     enddo
     do j=1,nda
        call csint(am1u,zee(j,l),tt2,tt,tt,tt2,c(j),sl(j),sl1(j))
     enddo
     pdu=0.d0
     if(m==0 .and. nplk1 > 0) then
        call pkint(am1u,tt,tt,tt2,nplk1,pki)
        do j=1,nplk1
           sum1=sum(pr1(i,1:nda)*dpe(1:nda,j,l)+pt1(i,1:nda)*dme(1:nda,j,l))
           pdu=pdu+(w*sum1+2*pi*(1-w)*cplk(j,l))*pki(j)
        enddo
     endif
     do j=1,na0
        trns0=-t1/am0(j)
        trns0=expfn(trns0)*fsol
        pvu=pr10(i,j)*trns0
        hu=sum((h1(1:nda)*c(1:nda)+h2(1:nda)*sl(1:nda))*alfa(1:nda,j) &
             +(h1(1:nda)*sl1(1:nda)+h2(1:nda)*c(1:nda))*beta(1:nda,j))
        pvu=pvu+sum(pr1(i,1:nda)*vpe(1:nda,j,l)+pt1(i,1:nda)*vme(1:nda,j,l))
        pvu=w*pvu*ei (j)
        ai(j)=(hu+pvu+pdu)/am1
     enddo
  endif
  return
end subroutine ainteg

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine chkrt(inda,indt,indp,imthd,nda,na1u,am1u,na0,am0,nfi,fi,&
     nln,thk,omg,nlgn1,g,nang,ang,phsf,epsp,epsu,galb,fsol,nplk1,cplk,&
     bgnd,ntau,utau,amua,wa,dpt,mxlgn2)
! check-in variables for -rtrn1- and set some variables.
!--- history
! 90. 1.20   created
!    11.22   add: else mxlgn2=2*nda+1
!    12. 1   add loop-19.
! 93. 5. 4   delete the condition galb<=1 because galb is used as wind
!            velocity when indg>0.
! input/output see main routine.
! 08. 9.29 modified Fortran 90 free-style form
!--
  use module_opt_visir, only: kna1u,kna0,kndm,knfi,klgn1,knln,kntau,kplk1,knang, &
       knln1,pi
  use module_simulator
  implicit none

! input   
  integer,intent(in):: inda      !! radiance flag (0: flux only, 1: flux & rad)
  integer,intent(in):: indt      !! flag of optical depth (1:utau, 2:dpt)
  integer,intent(in):: indp      !! flag of phase function
  integer,intent(in):: imthd     !! flux/intensity calculation modes
  integer,intent(in):: nda       !! number of stream
  integer,intent(inout):: na1u   !! number of zenith angles of emerging rad
  real(8),intent(inout):: am1u(kna1u) !! cosine of the emergent nadir angles
  integer,intent(in):: na0           !! number of solar zenith angles
  real(8),intent(in):: am0(kna0)     !! cosine of solar zenith angles
  integer,intent(in):: nfi           !! azimuthal angles
  real(8),intent(in):: fi(knfi)      !! zenith angles (degree)

  integer,intent(in):: nln           !! number of sublayers
  real(8),intent(in):: thk(knln)     !! opthical thickness
  real(8),intent(in):: omg(knln)     !! single scattering albedo
  integer,intent(in):: nlgn1(knln) !! max order of moments+1 (indp<=0;imthd=3)
  real(8),intent(in):: g(klgn1,knln) !! Legendre moments of phase function (indp<=0)
  integer,intent(in):: nang          !! number of scaattering angles
  real(8),intent(in):: ang(knang)    !! scattering angles [degree]
  real(8),intent(in):: phsf(knang,knln)  !! phase function (indp=1)
  real(8),intent(in):: epsp      !! truncation criterion of phase func. moments
  real(8),intent(in):: epsu      !! convergence criterion of intensity (inda>0)
  real(8),intent(in):: galb      !! ground albedo (indg=0)/ U10 [m/sec] (indg=1)
  real(8),intent(in):: fsol      !! solar irradiance at the system top
  integer,intent(in):: nplk1     !! number of order to approximate plank + 1
  real(8),intent(inout):: cplk(kplk1,knln)  !! coefficient of planck function 
  real(8),intent(inout):: bgnd   !! planck function at ground 
  integer,intent(inout):: ntau   !! number of user defined optical depth
  real(8),intent(inout):: utau(kntau) !! optical depths where the field is calc

  real(8),intent(out):: amua(kndm)   !! 
  real(8),intent(out):: wa(kndm)     !! 
  real(8),intent(out):: dpt(knln1)   !! 
  integer,intent(out):: mxlgn2   !! 
  real(8) :: p !phase function
! areas for this routine
  real(8):: ccp(3)
  real(8):: eps
  integer:: i,l,it,n1
!--exec
  call cpcon(ccp)
  eps=ccp(1)*10.d0

! epsp, epsu, fsol
  if(epsp < 0.d0) then
     stop 'MSG chkrt: illegal value of epsp'
  endif
  if(epsu < 0.d0) then
     stop 'MSG chkrt: illegal value of epsu'
  endif
  if(fsol < 0.d0) then
     stop 'MSG chkrt: illegal value of fsol'
  endif

! inda
  if(inda < 0 .or. inda > 2) then
     stop 'MSG chkrt: illegal value of inda'
  endif

! indt
  if(indt < 0 .or. indt > 2) then
     stop 'MSG chkrt: illegal value of indt'
  endif

! indp
  if(indp < -1 .or. indp > 1) then
     stop 'MSG chkrt: illegal value of indp'
  endif

! imthd
  if(imthd > 3) then
     stop 'MSG chkrt: illegal value of imthd'
  endif

! nda
  if(nda <= 0 .or. nda > kndm) then
     stop 'MSG chkrt: illegal value of nda'
  endif
 
! set quadrature
!! shifted gaussian quadrature on (0, 1).
  call qgausn(wa,amua,nda)
!! reordering.
  amua(1:nda)=amua(nda:1:-1)
  wa(1:nda)=wa(nda:1:-1)

! na0, am0
  if(na0 <= 0 .or. na0 > kna0) then
     stop 'MSG chkrt: illegal value of na0'
  endif

  do i=1,na0
     if(am0(i) <= 0.d0 .or. am0(i) > 1.d0) then
        stop 'MSG chkrt: illegal value of am0'
     endif
  enddo

! inda, am1u, fi
  if(inda > 0) then
     if(nfi <= 0 .or. nfi > knfi) then
        stop 'MSG chkrt: illegal value of nfi'
     endif
     if(inda == 2) na1u=2*nda
     if(na1u <= 0 .or. na1u > kna1u) then
        stop 'MSG chkrt: illegal value of na1u'
     endif
     if(inda == 2) then
        am1u(1:nda)=-amua(1:nda)
        am1u(na1u:na1u-nda+1:-1)=amua(1:nda)
     endif
  endif

! nln
  if(nln <= 0 .or. nln > knln) then
     print*,'nln=',nln,'knln=',knln
     stop 'MSG chkrt: illegal value of nln'
  endif

  dpt(1)=0.d0
  do l=1,nln
     dpt(l+1)=dpt(l)+thk(l)
  enddo

! ntau
  if(indt == 2) ntau=nln+1
  if(ntau <= 0 .or. ntau > kntau) then
     stop 'MSG chkrt: illegal value of ntau'
  endif

! utau
  if(indt == 2) then
     utau(1:ntau)=dpt(1:ntau)
  endif

! nlgn1,  mxlgn2=max(nlgn1)
  if(indp <= 0 .or. (inda > 0 .and. imthd == 3)) then
     mxlgn2=1
     do l=1,nln
        n1=nlgn1(l)
        if(n1 <= 0 .or. n1 > klgn1) then
           stop 'MSG chkrt: illegal value of nlgn1'
        endif
        mxlgn2=max(mxlgn2,n1)
     enddo
  else
     mxlgn2=2*nda+1
  endif


! nang
  if(indp > 0 .or. (imthd >= 1 .and. inda >= 1)) then
     if(nang <= 3 .or. nang > knang) then
        print*,'nang=',nang,'knang=',knang
        stop 'MSG chkrt: you should set at least four angles for this condition'
     endif
     do i=2,nang
        if(ang(i) <= ang(i-1)) then
           stop 'MSG chkrt: your should set scattering angle from 0 to 180 degrees'
        endif
     enddo
  endif

! phsf
  do i=1,nang ; do l=1,nln
     p=phsf(i,l)
     if( p==0. .or. p>0. .or. p<0. .or. p>=0. .or. p<=0. ) then
     else
       print*,'MSG chkrt: find NaN in phsf, nang=',i,' lyr=',l,' phsf=',phsf(i,l)
     endif
  enddo ; enddo 

! check order of dpt and utau.
  do l=1,nln
     if(dpt(l+1) < dpt(l)) then
        stop 'MSG chkrt: dpt should be set from top to bottom'
     endif
  enddo

  if(ntau >= 2) then
     do it=2,ntau
        if(utau(it) < utau(it-1)) then
           stop 'MSG chkrt: utau should be set from top to bottom'
        endif
     enddo
  endif
  do it=1,ntau
     if(utau(it) > dpt(nln+1)) then
        if(abs(utau(it)-dpt(nln+1)) <= eps) then
           utau(it)=dpt(nln+1)
        else
           stop 'MSG chkrt: utau is out of bounds'
        endif
     endif
  enddo

! reset cplk and bgnd, galb
  if(galb < 0.d0) then
     stop 'MSG chkrt: illegal value of galb'
  endif

  if(nplk1 > kplk1) then
     stop 'MSG chkrt: too large -nplk1-, change -kplk1-'
  else
     if(nplk1 <= 0) then
        bgnd=0.d0
        cplk(1:nplk1,1:nln)=0.d0
     endif
  endif
  return
end subroutine chkrt

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine cingr(nda,na0,l,am0,tau,rdn,rup,vpe,vme,c1e,c2e,nplk1,dpe,dme,alfa,beta)
! get alfa and beta (integral constants).
!--- input
! nda        i       no. of zenith-quadrature angles in the hemisphere
! na0        i       number of the solar incidences.
! l
! am0
! tau
! rup   r(kndm,kna0,knsb1) upwelling internal intensities.
! rdn   r(kndm,kna0,knsb1) downgoing internal intensities.
! vpe
! vme
! c1e
! c2e
! nplk1      i       number of order to approximate plank + 1.
! dpe
! dme
!--- output
! alfa
! beta
!--
  use module_opt_visir, only:kna0,kndm,knln,kplk1,knln1,knlnm,knlnm1,pi,rad
  implicit none

! input 
  integer,intent(in):: nda
  integer,intent(in):: na0
  integer,intent(in):: l
  real(8),intent(in):: am0(kna0)
  real(8),intent(in):: tau
  real(8),intent(in):: rup(kndm,kna0,knlnm1)
  real(8),intent(in):: rdn(kndm,kna0,knlnm1)
  real(8),intent(in):: vpe(kndm,kna0,knln)
  real(8),intent(in):: vme(kndm,kna0,knln)
  real(8),intent(in):: c1e(kndm,kndm,knln)
  real(8),intent(in):: c2e(kndm,kndm,knln)
  integer,intent(in):: nplk1
  real(8),intent(in):: dpe(kndm,kplk1,knln)
  real(8),intent(in):: dme(kndm,kplk1,knln)

! output
  real(8),intent(out):: alfa(kndm,kna0)
  real(8),intent(out):: beta(kndm,kna0)

! work areas
  integer::l1,i,j,k
  real(8)::ex1,sum1,sum2,sum3,buf1(kndm),buf2(kndm),expfn
!--exec
  l1=l+1
  do j=1,na0
     ex1=-tau/am0(j)
     ex1=expfn(ex1)
     do i=1,nda
        sum1=rdn(i,j,l)-vpe(i,j,l)
        sum2=rup(i,j,l1)-vme(i,j,l)*ex1
        if(nplk1 > 0) then
           sum1=sum1-dpe(i,1,l); sum3=0.d0
           do k=1,nplk1
              sum3=sum3+dme(i,k,l)*tau**(k-1)
           enddo
           sum2=sum2-sum3
        endif
        buf1(i)=sum2+sum1
        buf2(i)=sum2-sum1
     enddo
     do i=1,nda
        alfa(i,j)=sum(c1e(i,1:nda,l)*buf1(1:nda))
        beta(i,j)=sum(c2e(i,1:nda,l)*buf2(1:nda))
     enddo
  enddo
  return
end subroutine cingr

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine convu(n,kn,am,jj,er,ic,epsu,nchk)
! check convergence of intensity in the direction of mu
!--- history
! 90. 1.20 created
! 08.10.20 devided n into size and flag, modified fortran 90 version.
!--- input
! n      i        number of streams to be checked
! am    r(n)      direction of stream
! jj    i(n)      stream number in the original order
! er    r(n)      maximum error for the stream regardless of other
!                  angles and layers
! ic    i(n)      number of consective series with -er- less than -epsu-
! epsu    r       convergence criterion
! nchk    i       max. number for the value of -ic- by which the
!                  routine confirms convergence
!--- output
! n, am, jj, ic   updated after dropping convergent streams
!--

!
! Note bug report  input parametere jj must be declared as integer. 
! Then, yet for some reasons, jj values behaves crazy within do loop
! that I don't quite understand yet. 
! I temporally fixed the issue by adding jj_convu. 
!
  implicit none
  save
! input & output
  integer,intent(inout):: n
  integer,intent(in):: kn
  real(8),intent(inout):: am(kn)
  integer,intent(inout):: jj(kn)
  real(8),intent(inout):: er(kn)
  integer,intent(inout):: ic(kn)  !debug from EK
  real(8),intent(in):: epsu
  real(8),intent(in):: nchk

! work
  integer::j,j1
  integer :: jj_convu(kn)

!--exec
  jj_convu = jj
  j1=0
  do j=1,n
     if(abs(abs(am(j))-1.d0) <= 0.d0) then
        ic(j)=999
        er(j)=0.d0
     endif

     if(er(j) < epsu) then
        ic(j)=ic(j)+1
        if(ic(j) >= nchk) cycle
     else
        ic(j)=0
     endif
     j1=j1+1
     ic(j1)=ic(j)
     jj_convu(j1)=jj_convu(j)
     am(j1)=am(j)
     er(j1)=er(j)
  enddo

  jj = jj_convu

  n=j1
  return
end subroutine convu

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine csfn(tc,t,ze,c,sl,sl1)
! get c sl and sl1 functions.
!--- history
! 89.11. 2  created
! 08.10.27 modified fortran 90 version.
!---
! 2 c = e(tc-t) + e(t),   2s = e(tc-t) - e(t)
!   e(t) = exp(-ze*t)
! sl = s l,   sl1 = s l**-1
!--
  implicit none
  real(8),intent(in)::tc,t,ze
  real(8),intent(out)::c,sl,sl1
  integer,save::init=1
  real(8),save::eps
  real(8)::ccp(3),s,exp1,exp2,exp3,exp4,expfn
!--exec
  if(init==1) then
     init=0
     call cpcon(ccp)
     eps=ccp(1)*10.d0
  endif
  exp1=-ze*(tc-t); exp2=-ze*t
  exp3=expfn(exp1); exp4=expfn(exp2)
  c  =(exp3+exp4)/2.d0; s  =(exp3-exp4)/2.d0
  sl =ze*s
  if(abs(exp1) <= eps .and. abs(exp2) <= eps) then
     sl1 =t-tc/2.d0
  else
     sl1 =s/ze
  endif
  return
end subroutine csfn

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine csint(am1u,ze,thk,dpt,t1,t2,c,sl,sl1)
! integral(from t1 to t2) dt exp(-(dpt-t)/am1u) c(t, thk)
!   where c(t, thk) = (exp(-z*(thk-t)+exp(-z*t))/2
! we also define the integral for the functions: s*z, s/z.
!--- history
! 90. 1.13 created
! 08.10.29 modified fortran 90 version.
!--- input
! am1u     r      userdefined emergent mu.
!                 if  >  0 then down,  if  <  0 then upward.
! ze       r      eigenvalue (l).
! thk      r      optical thickness of layer.
! dpt      r      optical depth for interpolation.
! t1       r      lower limit of integration.
! t2       r      upper limit of integration.
!--- output
! c        r      integration of c.
! sl       r      integration of s*l.
! sl1      r      integration of s/l.
!--
  implicit none

! input 
  real(8),intent(in):: am1u
  real(8),intent(in):: ze
  real(8),intent(in):: thk
  real(8),intent(in):: dpt
  real(8),intent(in):: t1
  real(8),intent(in):: t2

! output
  real(8),intent(out):: c
  real(8),intent(out):: sl
  real(8),intent(out):: sl1

! work
  integer,save:: init=1
  real(8),save:: eps
  real(8):: ccp(3)
  real(8):: ex,x,ex1,ex2,ex3

! for function
  real(8):: exx,expfn
!--exec
  if(init == 1) then
     init=0
     call cpcon(ccp)
     eps=ccp(1)*100.d0
  endif
! c
  ex=-dpt/am1u-ze*thk
  x=1.d0/am1u+ze
  ex1=exx(ex,t1,t2,x)
  ex=-dpt/am1u
  x=1.d0/am1u-ze
  ex2=exx(ex,t1,t2,x)
  c =(ex1+ex2)/2.d0
  if(abs(ze*t1) > eps .or. abs(ze*t2) > eps) then
     sl1=(ex1-ex2)/2.d0/ze
  else
     ex=(t1-dpt)/am1u
     ex1=expfn(ex)
     ex=(t2-dpt)/am1u
     ex2=expfn(ex)
     ex=-dpt/am1u
     x=1.d0/am1u
     ex3=exx(ex,t1,t2,x)
     sl1=am1u*(t2*ex2-t1*ex1)-(am1u+thk/2.d0)*ex3
  endif
  sl =ze**2*sl1
  return
end subroutine csint

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine exint(am1u,am0,dpt,t1,t2,c)
! integral(from t1 to t2) dt exp(-(thk-t)/am1u-t/am0)
!--- history
! 90. 1.13 created
! 08.10.29 modified fortran 90 version.
!--- input
! am1u     r      userdefined emergent mu.
!                 if  >  0 then down,  if  <  0 then upward.
! am0      r      solar direction.
! dpt      r      optical depth for interpolation.
! t1       r      lower limit of integration.
! t2       r      upper limit of integration.
!--- output
! c        r      integration.
!--
  implicit none

! input
  real(8),intent(in):: am1u
  real(8),intent(in):: am0
  real(8),intent(in):: dpt
  real(8),intent(in):: t1
  real(8),intent(in):: t2

! output
  real(8),intent(out):: c

! work
  real(8):: ex,x,exx
!--exec  
  ex=-dpt/am1u
  x=1/am1u-1/am0
  c=exx(ex,t1,t2,x)
  return
end subroutine exint

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function exx(ex,t1,t2,x)
! exp(ex)*integ(from t1 to t2) dt exp(xt)
!--- history
! 89.11. 3 modified
! 08.10.29 modified fortran 90 version.
!--
  implicit none

! input & output
  real(8),intent(in):: ex
  real(8),intent(in):: t1
  real(8),intent(in):: t2
  real(8),intent(in):: x
  real(8):: exx

! work
  integer,save:: init=1
  real(8),save:: eps
  real(8),save:: exmn 
  real(8):: ccp(3)

  real(8):: ex1,ex2,ex3,t22,expfn
!--
  if(init == 1) then
     init=0
     call cpcon(ccp)
     eps=ccp(1)*100.d0
     exmn=ccp(2)*0.8d0*2.3d0
  endif
  ex1=ex+x*t1
  ex2=expfn(ex1)
  t22=t2-t1
  if(abs(t22*x) <= eps) then
     exx=ex2*(t22+x*t22**2/2.d0)
  else
     ex3=ex+x*t2
     if(ex3 <= exmn) then
        exx=-ex2/x
     else
        exx=(exp(ex3)-ex2)/x
     endif
  endif
  return
end function exx

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine ftrn21(indg,init,m,inda,indt,imthd,na1u,am1u,na0,am0, &
     nda,amua,wa,nln,nlgn1,g,dpt,omg,nplk1,cplk,galb,bgnd,fsol, &
     ntau,utau,scr,sci,flxd,flxu,ai,rt_failed)
! solve the radiative transfer in atmosphere system for each fourier
! component.
!  -dom and adding method-
! by Teruyuki Nakajima
! since this system does not include a ocean surface,
! knlt=knln in trn1.
!--- history
! 89.11. 2 created from strn7 with thermal emission
! 90.12. 1 call aint(m,l,  ->  call aint(m,lt,
! 92.12.   add ocean surface
! 93. 1.   change upward adding
!            use cspl1 and cspli
!     3.10 bug ( xx,yy used k1 )
!     4. 5 nsf = 1 for initialization (terry)
!     5. 4 indg instead of nsf and losw
!          subtracted single scattering if indg>0
! 94. 5. 7 bug for interpolating aib(j)
! 95. 5.25 replace grndl by grndl3
! 95. 6. 2 with grndo1 and with scr and sci
! 96. 3.17 add ert by Takashi Nakajima to debug thermal. 
! 02.05.14 aint --> ainteg
! 08.10.24 modified fortran 90 version.
!--- input
! indg       i       0: lambert surface
!                    1: ocean surface initialized
!                    2: ocean surface with no-initialization
!                    when indg>0 and imthd>0 then single scattering correction
!                      for ocean surface reflection
! init       i       1: initialize the part (dependent amua, utau)
!                    0:  bypass the m-independent part.
! m          i       forier order.
! inda       i       0: flux only.
!                    1: intensity using am1u.
!                    2: na1u and am1u are supposed to be
!                       2*nda and (-amua, +amua).
! indt       i       0: set user defined depth.
!                    1: same as above.
!                    2: set ntau and utau as nln1 and dpt.
! imthd      i      -1: nt,  0: dms-method  for intensity/flux
!                    1: ms,  2:tms,  3:ims-method for intensity.
!                    when indg>0 and imthd>0 then single scattering correction
!                      for ocean surface reflection
! na1u       i       number of emergent zenith angles in the sphere.
!                      na1u=2*nda when inda=2.
! am1u    r(kna1u)   consine of the emergent zenith angles.
!                      + for downward, - for upward.
!                      am1u = (-amua, +amua) when inda=2.
! na0        i       number of the solar incidences.
! am0     r(na0)     consine of solar zenith angles  > 0.
! nda        i       no. of zenith-quadrature angles in the hemisphere.
! amua    r(kndm)    cosine of the zenith-quadrature angles. 1 to 0.
! wa      r(kndm)    corresponding quadrature weights.
! nln        i       number of atmospheric sublayers.
! nlgn1   i(knln)    maximum order of the legendre sesies of the phase
!                      function + 1
! g     r(klgn1,     legendre moments of phase function.
!           knln)      g(1,l)=1
! dpt     r(knln1)   optical depth at the interfaces between sublayers.
!                      top to bottom (top = 0 for normal application).
! omg     r(knln)    single scattering albedo.
! nplk1      i       number of order to approximate plank + 1.
!                      if 0 then no thermal.
! cplk    r(kplk1    plank function (b) =
!         ,knln)       sum(ib=1,nplk1) cplk(ib,l) tau**(ib-1).
!                      tau is optical depth measured from
!                      the top of the subsurface.
! galb       r       ground albedo if indg=0
!                    u10 (m/sec)   if indg>0
! bgnd       r       (1-galb)*b when indg=0
!                    b          when indg>0
!                    where b is plank function
!                      (same unit as fsol).
! fsol       r       solar irradiance (w/m2/micron).
! ntau       i       number of user defined optical depths.
!                      ntau=nln+1 when indt=2.
! utau    r(kntau)   optical depths where the field is calculted.
!                      top to bottom.
!                      utau=dpt when indt=2.
!--- output
! indg               if 1 then 2
! init               0
! flxd  r(kna0,kntau) downward flux at utau.
! flxu               same as flxd but for upward flux.
! ai      r(kna1u,   i(mu1u(i), mu0(j), l)
!           kna0,    intensity at utau.
!           kntau)   subtracted single scattering if indg>0
!--- correspondence between variables and parameters
! kna1u       na1u
! kna1        na1
! kna0        na0
! kndm        nda
! knln        nln
! knln1       knln+1
! kntau       ntau
! klgn1       nlgn1
! kplk1       nplk1
!--- notes for local variables
! na1        i       number of streams for am1.
! am1     r(kna1)    abs(am1u).
! iiam1   i(kna1u)   direction number of am1 for each am1u.
!                      + for downward, - for upward.
! iitau   i(kntau)   sublayer number for user defined depths.
!--
  use module_opt_visir, only: kna1u,kna0,kndm,kna1,knln,kntau,kplk1,knln1,klgn1, &
       knltm,knlnm,knlnm1,kndm2,pi,rad
  implicit none

! input
  integer,intent(in):: indg      !! surface flag (0:Lambert, 1 or 2: Ocean)
  integer,intent(inout):: init      !! initialization flag 
  integer,intent(in):: m         !! stream number
  integer,intent(in):: inda      !! radiance flag (0: flux only, 1: flux & rad)
  integer,intent(in):: indt      !! flag of optical depth (1:utau, 2:dpt)
  integer,intent(in):: imthd     !! flux/intensity calculation modes

  integer,intent(in):: na1u      !! number of zenith angles of emerging rad
  real(8),intent(in):: am1u(kna1u)   !! cosine of the emergent nadir angles
  integer,intent(in):: na0           !! number of solar zenith angles
  real(8),intent(in):: am0(kna0)     !! cosine of solar zenith angles
  integer,intent(in):: nda           !! number of stream
  real(8),intent(in):: amua(kndm)   !! cosine of the zenith-quadrature angles. 1 to 0.
  real(8),intent(in):: wa(kndm)     !! corresponding quadrature weights.

  integer,intent(in):: nln           !! number of sublayers
  integer,intent(in):: nlgn1(knln) !! maximum order of moments + 1 (indp<=0;imthd=3)
  real(8),intent(in):: g(kndm2,knln) !! Legendre moments of phase function (indp<=0)
  real(8),intent(in):: dpt(knln1)   !! 
  real(8),intent(in):: omg(knln)     !! single scattering albedo
  integer,intent(in):: nplk1     !! number of order to approximate plank + 1
  real(8),intent(in):: cplk(kplk1,knln)  !! coefficient of planck function 

  real(8),intent(in):: galb      !! ground albedo (indg=0)/ U10 [m/sec] (indg=1)
  real(8),intent(in):: bgnd      !! planck function at ground 
  real(8),intent(in):: fsol      !! solar irradiance at the system top
  integer,intent(in):: ntau      !! number of user defined optical depth
  real(8),intent(in):: utau(kntau)   !! optical depths where the field is calc
  real(8),intent(in):: scr        !! relative reflactive index (real part)
  real(8),intent(in):: sci        !! relative reflactive index (imaginary part)

! output
  real(8),intent(out):: flxd(kna0,kntau)     !! downward flux
  real(8),intent(out):: flxu(kna0,kntau)     !! upward flux
  real(8),intent(out):: ai(kna1u,kna0,kntau) !! intensity
  logical,intent(inout) :: rt_failed  !logic if radiative tarnsfer is failed, --> true (default is fault)

! for cpcon
  real(8):: ccp(3)
  real(8),save:: eps
  real(8),save:: wmp(kndm),wmm(kndm)

  integer:: it,l
  real(8):: odp,eps1
  integer,save:: iitau(kntau)

  integer:: iu,j,k
  integer:: na1,iiam1(kna1u)
  real(8):: am1(kna1),x

  integer:: m1,mxlgn1
! for plgnd
  real(8):: pl0(kna0,kndm2),pla(kndm,kndm2)

  integer:: lb,nal1,ndd(knlnm1)
  integer:: iup(knlnm),idn(knlnm)
  real(8):: t1,t2,w,cpl(kplk1),gbuf(klgn1)
! for phas2
  real(8):: pt(kndm,kndm),pr(kndm,kndm),pt0(kndm,kna0),pr0(kndm,kna0)

! for homog2
  real(8):: r(kndm,kndm),t(kndm,kndm),er(kndm,kna0),et(kndm,kna0),zeig(kndm)
  real(8):: q(kndm,kndm),qi(kndm,kndm),c11(kndm,kndm),c22(kndm,kndm)
  real(8):: vp(kndm,kna0),vm(kndm,kna0),dp(kndm,kplk1),dm(kndm,kplk1)

  real(8):: re(kndm,kndm,knltm),  te(kndm,kndm,knltm)
  real(8):: ser(kndm,kna0,knlnm),set(kndm,kna0,knlnm)
  real(8):: qe(kndm,kndm,knln),qie(kndm,kndm,knln)
  real(8):: c1e(kndm,kndm,knln),c2e(kndm,kndm,knln),zee(kndm,knln)
  real(8):: vpe(kndm,kna0,knln),vme(kndm,kna0,knln)
  real(8):: dpe(kndm,kplk1,knln),dme(kndm,kplk1,knln)

  integer:: nlt,i,nln1,nlt1
  real(8):: ert(kndm,kna0)
  real(8):: rup(kndm,kna0,knlnm1),rdn(kndm,kna0,knlnm1),rup1

  real(8):: ex1,trns0

  real(8):: pl1(kna1,kndm2)
  integer:: lt
  real(8):: dpth,tau
! for cingr  
  real(8):: alfa(kndm,kna0),beta(kndm,kna0)

  integer:: l1
  real(8):: dtau
! for adisc
  real(8):: uup(kndm,kna0),udn(kndm,kna0)
  integer:: lu

  real(8):: pr1(kna1,kndm),pt1(kna1,kndm),pr10(kna1,kna0),pt10(kna1,kna0)
! for ainteg
  real(8):: aii(kna1u,kna0,knln1),aib(kna0)

  integer:: k1
  real(8):: xx(kndm),yy(kndm),a(kndm),b(kndm),c(kndm),d(kndm)

! for function
  real(8):: expfn,cspli
  real(8):: x1
!--exec
  rt_failed = .false. !default

  if(init > 0) then
     init=0
     call cpcon(ccp)
     eps=ccp(1)*10.d0

     wmp(1:nda)=sqrt(wa(1:nda)*amua(1:nda))
     wmm(1:nda)=sqrt(wa(1:nda)/amua(1:nda))

! set sublayer where user-defined depth resides -iitau-.
     do it=1,ntau
        odp=utau(it)
        eps1=odp*eps

        do l=1,nln
!print*,l,nln,(odp-dpt(l+1)),odp,dpt(l+1),it
           if((odp-dpt(l))*(odp-dpt(l+1)) <= 0.d0) exit
           if(l==nln) then
              print*,'MSG ftrn21: out of bounds in layer setting or optical properties are NaN'
              print*,'odp=',odp,'  dpt(l)=',dpt(l),'  dpt(l+1)=',dpt(l+1)
              rt_failed = .true. ; return !add skip comments here. 
              !stop
           endif
        enddo

        if(abs(odp-dpt(l)) <= eps1) then
           iitau(it)=-l
        else
           if(abs(odp-dpt(l+1)) <= eps1) then
              iitau(it)=-(l+1)
           else
              iitau(it)=l
           endif
        endif
     enddo
  endif

! set am1, iiam1
  if(inda > 0) then
     na1=1
     am1(1)=abs(am1u(1))
     iiam1(1)= 1
     if(na1u >= 2) then
        do iu=2,na1u
           x=abs(am1u(iu))
           do j=1,na1
              if(abs(x-am1(j)) <= eps) then
                 iiam1(iu)= j; exit
              endif
           enddo
           na1=na1+1
           if(na1 > kna1) then
              stop 'MSG ftrn21: setting error of am1u'
           endif
           am1(na1)=x; iiam1(iu)= na1
        enddo
     endif
  endif

! check maximum nlgn1
  mxlgn1=maxval(nlgn1(1:nln))
  if(mxlgn1 <= 0) then
     stop 'MSG ftrn21: nlgn1 are all zero'
  endif
  m1=m+1

! set legendre polynomials.
  call plgnd(m1,mxlgn1,na0,kna0,am0 ,pl0)
  call plgnd(m1,mxlgn1,nda,kndm,amua,pla)

! scaling for symmetlicity
  do j=1,mxlgn1
     pla(1:nda,j)=wmm(1:nda)*pla(1:nda,j)
  enddo

! solve the eigenvalue problem of atmospheric sublayers.
  do l=1,nln
     lb=l
     nal1=nlgn1(l);  ndd(l)=nda
     iup(l)=l;  idn(l)=l
     t1=dpt(l); t2=dpt(l+1)
     w=omg(l)
     if(nplk1 > 0) then
        cpl(1:nplk1)=2.d0*pi*(1.d0-w)*cplk(1:nplk1,l)
     endif
!! scattering media
!! pt, pr
     gbuf(1:nal1)=g(1:nal1,lb)
     call phas2(m1,nal1,nda,nda,kndm,kndm,kndm,gbuf,pla,pla,pt,pr)

!! pt0, pr0
     call phas2(m1,nal1,nda,na0,kndm,kndm,kna0,gbuf,pla,pl0,pt0,pr0)

!! eigenvalue problem
     call homog2(m,t1,t2,w,nda,amua,wmm,na0,am0,pr,pt,pr0,pt0,fsol, &
          nplk1,cpl,r,t,er,et,zeig,q,qi,c11,c22,vp,vm,dp,dm,rt_failed)
 
     if( rt_failed ) then
        flxd=0. ; flxu=0. ; ai = 0. 
     endif
 
     re(1:nda,1:nda,lb)=r(1:nda,1:nda)
     te(1:nda,1:nda,lb)=t(1:nda,1:nda)
     ser(1:nda,1:na0,lb)=er(1:nda,1:na0)
     set(1:nda,1:na0,lb)=et(1:nda,1:na0)

! store extra for internal field at arbitrary depth.
     if(indt /= 2 .or. inda /= 2) then

        qe(1:nda,1:nda,lb)=q(1:nda,1:nda)
        qie(1:nda,1:nda,lb)=qi(1:nda,1:nda)
        c1e(1:nda,1:nda,lb)=c11(1:nda,1:nda)
        c2e(1:nda,1:nda,lb)=c22(1:nda,1:nda)
        vpe(1:nda,1:na0,lb)=vp(1:nda,1:na0)
        vme(1:nda,1:na0,lb)=vm(1:nda,1:na0)
        zee(1:nda,lb)=zeig(1:nda)
        if(nplk1 > 0) then
           dpe(1:nda,1:nplk1,lb)=dp(1:nda,1:nplk1)
           dme(1:nda,1:nplk1,lb)=dm(1:nda,1:nplk1)
        endif
     endif
  enddo

! Lambert surface
  nln1=nln+1
  if(indg < 0) then
     nlt=nln
  else
     nlt=nln1
     ndd(nlt)=nda
     iup(nlt)=nlt; idn(nlt)=nlt
     t1=dpt(nlt)
     if(indg <= 0)then  !land (lambertian)
        call grndl3(fsol,galb,bgnd,t1,m,nda,amua,wa,na0,am0,r,t,er,et)
     else  !ocean (indg is 1 or 2)
        call grndo1(indg,fsol,galb,bgnd,t1,m,nda,wmp,na0,am0,scr,sci, &
             r,t,er,ert,et,amua,wa)
     endif
     re(1:nda,1:nda,nlt)=r(1:nda,1:nda)
     te(1:nda,1:nda,nlt)=t(1:nda,1:nda)
     ser(1:nda,1:na0,nlt)=er(1:nda,1:na0)
     set(1:nda,1:na0,nlt)=et(1:nda,1:na0)
  endif
  nlt1=nlt+1
  ndd(nlt1)=nda

! Adding of the sublayers.
  call trn1(nlt,ndd,na0,iup,idn,re,te,ser,set,rup,rdn)

  if(indt == 2 .and. inda /= 1) then
     if(inda == 2) then
        do l=1,nln1
           do j=1,na0
              if(imthd <= 0 .or. indg <= 0) then
                 ai(1:nda,j,l)=rup(1:nda,j,l)/wmp(1:nda)
              else
                 do i=1,nda
! don't change rup itself, because it is used for flux calculation later
! 97.3.17 debug by takashi nakajima
!            rup1=rup(i,j,l)-ser(i,j,nln1)
!     &                     *exp(-(dpt(nln1)-dpt(l))/amua(i))
                    rup1=rup(i,j,l)+(ert(i,j)-ser(i,j,nln1)) &
                         *exp(-(dpt(nln1)-dpt(l))/amua(i))
                    ai(i,j,l)=rup1/wmp(i)
                 enddo
              endif
              ai(na1u:na1u+1-nda:-1,j,l)=rdn(1:nda,j,l)/wmp(1:nda)
           enddo
        enddo
     endif
     if(m == 0) then
        do l=1,nln1
           do j=1,na0
              ex1=-utau(l)/am0(j)
              trns0=expfn(ex1)*fsol
!! for scaled intensity
              flxu(j,l)=sum(wmp(1:nda)*rup(1:nda,j,l))
              flxd(j,l)=sum(wmp(1:nda)*rdn(1:nda,j,l))+trns0*am0(j)
           enddo
        enddo
     endif
     return
  endif

!--- interpolation of the field.
  if(inda == 1) call plgnd(m1,mxlgn1,na1,kna1,am1 ,pl1)
  do l=1,nln
     lt=l
     lb=lt+1
     dpth=dpt(lt)
     tau=dpt(lb)-dpth
     t1=dpt(lt); t2=dpt(lb)
     w=omg(l)

! integral constants:  alfa, beta
     call cingr(nda,na0,lt,am0,tau,rdn,rup,vpe,vme,c1e,c2e,nplk1, &
          dpe,dme,alfa,beta)
     
! fluxes or intensity when inda=2
     if(m == 0 .or. inda == 2) then
        do it=1,ntau
           l1=iabs(iitau(it))
           if(l1 == nln1) l1=nln
           if(l1 /= l) cycle
           dtau=utau(it)-dpth
           if(iitau(it) > 0) then
!! in the sublayer
              call adisc(m,l1,nda,na0,am0,wmp,tau,dtau,zee,qe,qie, &
                   vpe,vme,alfa,beta,nplk1,dpe,dme,udn,uup)
           else
!! just interface
              lu=iabs(iitau(it))
              do j=1,na0
                 udn(1:nda,j)=rdn(1:nda,j,lu)/wmp(1:nda)
                 uup(1:nda,j)=rup(1:nda,j,lu)/wmp(1:nda)
              enddo
           endif

!! intensity when inda=2
           if(inda == 2) then
              do j=1,na0
                 if(imthd <= 0 .or. indg <= 0) then
                    ai(1:nda,j,it)=uup(1:nda,j)
                 else
! don't change uup itself, because it is used for flux calculation later
! 97.3.17 debug by takashi nakajima
!   56         ai(i,j,it)=uup(i,j)-ser(i,j,nln1)/wmp(i)
!     &                  *exp(-(dpt(nln1)-utau(it))/amua(i))
                    ai(1:nda,j,it)=uup(1:nda,j) &
                         +(ert(1:nda,j)-ser(1:nda,j,nln1)/wmp(1:nda)) &
                         *exp(-(dpt(nln1)-utau(it))/amua(1:nda))
                 endif
                 ai(na1u:na1u+1-nda:-1,j,it)=udn(1:nda,j)
              enddo
           endif

!! flux
           if(m == 0) then
              do j=1,na0
                 ex1=-utau(it)/am0(j)
                 trns0=expfn(ex1)*fsol
                 flxu(j,it)=sum(amua(1:nda)*wa(1:nda)*uup(1:nda,j))
                 flxd(j,it)=sum(amua(1:nda)*wa(1:nda)*udn(1:nda,j))+trns0*am0(j)
              enddo
           endif
        enddo
     endif
     
!--- intensities in the user defined directions (without contribution
!     from the interfaces).
! phase functions for angular interpolation
     if(inda == 1) then
        nal1=nlgn1(l)
        gbuf(1:nal1)=g(1:nal1,lt)

! pt10, pr10
        call phas2(m1,nal1,na1,na0,kna1,kna1,kna0,gbuf,pl1,pl0,pt10,pr10)
       
! pt1, pr1
        call phas2(m1,nal1,na1,nda,kna1,kna1,kndm,gbuf,pl1,pla,pt1,pr1)
       
        do iu=1,na1u
! at interface.
           i=iiam1(iu)
           if(am1u(iu) > 0.d0) then
!! downward increment at bottom of sublayer (stamnes integration)
              call ainteg(m,lt,nda,na0,i,fsol,am1u(iu),am0,w,t2,t1,t2, &
                   zee,qe,qie,vpe,vme,pt10,pr10,pt1,pr1,alfa,beta, &
                   nplk1,dpe,dme,cplk,aib)
              aii(iu,1:na0,lb)=aib(1:na0)
           else
!! upward increment at top of sublayer (stamnes integration)
              call ainteg(m,lt,nda,na0,i,fsol,am1u(iu),am0,w,t1,t1,t2, &
                   zee,qe,qie,vpe,vme,pt10,pr10,pt1,pr1,alfa,beta, &
                   nplk1,dpe,dme,cplk,aib)
              aii(iu,1:na0,l)=aib(1:na0)
           endif

! at user defined depth
           do it=1,ntau
              l1=iitau(it)
!!  we do not calculate the field at this stage if the user defined
!!     depth is exactly on the interfaces (l1 < 0).
              if(l1 == l) then
                 call ainteg(m,lt,nda,na0,i,fsol,am1u(iu),am0,w,utau(it), &
                      t1,t2,zee,qe,qie,vpe,vme,pt10,pr10,pt1,pr1,alfa, &
                      beta,nplk1,dpe,dme,cplk,aib)
                 ai(iu,1:na0,it)=aib(1:na0)
              endif
           enddo
        enddo
     endif
  enddo
  if(inda <= 0 .or. inda == 2) return

!--- adding intensities (inda=1)
  do iu=1,na1u
     i=iiam1(iu)
     if(am1u(iu) > 0.d0) then
!! downward
        aib(1:na0)=0.d0
        do l=1,nln1
           if(l > 1) then
              ex1=-(dpt(l)-dpt(l-1))/am1(i)
              ex1=expfn(ex1)
              aib(1:na0)=aib(1:na0)*ex1+aii(iu,1:na0,l)
           endif
           do it=1,ntau
              lu=iitau(it)
              if(iabs(lu) == l) then
                 if(lu < 0) then
                    ai(iu,1:na0,it)=aib(1:na0)
                 else
                    dtau=utau(it)-dpt(l)
                    ex1=-dtau/am1(i)
                    ex1=expfn(ex1)
                    ai(iu,1:na0,it)=ai(iu,1:na0,it)+aib(1:na0)*ex1
                 endif
              endif
           enddo
        enddo

!! upward adding
     else
        do j=1,na0
! 93.5.4
           if(indg < 0) then
              aib(j)=0.d0
           else
              if(indg == 0) then
                 aib(j)=rup(1,j,nln1)/wmp(1)
              else
! interpolation of reflected intensities subtracted single scattering
                 if(imthd <= 0) then
                    do k=1,nda
                       k1=nda-k+1
! 94.5.7 bug  yy(k1)=rup(k,j,nln1)
! corrected
                       yy(k1)=rup(k,j,nln1)/wmp(k)
                       xx(k1)=amua(k)
                    enddo
                 else
                    do k=1,nda
                       k1=nda-k+1
! 94.5.7 bug  yy(k1)=rup(k,j,nln1)-ser(k,j,nln1)
! corrected
! 97.3.17 debug by takashi nakajima
!              yy(k1)=(rup(k,j,nln1)-ser(k,j,nln1))/wmp(k)
                       yy(k1)=(rup(k,j,nln1)-ser(k,j,nln1))/wmp(k) &
                            + ert(k,j)/wmp(k)
                       xx(k1)=amua(k)
                    enddo
                 endif
                 call cspl1(nda,xx,yy,a,b,c,d)
                 x1=am1(i)
! 94.5.7 bug  aib(j)=cspli(x1,nda,xx,a,b,c,d)/wmp(1)
! corrected
                 aib(j)=cspli(x1,nda,xx,a,b,c,d)
              endif
           endif
! 93.5.4 end
        enddo
        
        do l=nln1,1,-1
           if(l <= nln) then
              ex1=-(dpt(l+1)-dpt(l))/am1(i)
              ex1=expfn(ex1)
              aib(1:na0)=aib(1:na0)*ex1+aii(iu,1:na0,l)
           endif
           do it=1,ntau
              lu=iitau(it)
              if(lu < 0) then
                 if(iabs(lu) == l) ai(iu,1:na0,it)=aib(1:na0)
              else
                 if(lu == l-1) then
                    dtau=dpt(l)-utau(it)
                    ex1=-dtau/am1(i)
                    ex1=expfn(ex1)
                    ai(iu,1:na0,it)=ai(iu,1:na0,it)+aib(1:na0)*ex1
                 endif
              endif
           enddo
        enddo
     endif
  enddo
  return
end subroutine ftrn21

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine getqm(iw0,n,x,y,zeig,q,qi,xi,imn,rt_failed )
! solve    xy = q zeig**2 inverse(q)
! root decomposition method
!--- history
! 89. 8. 4 created
!--- input
! iw0       i        if 1 then renomalization (for m=0 and w0=1)
! n         i        order of matrices
! x     r(kndm,kndm) symmetric matrix
! y     r(kndm,kndm) symmetric matrix
!--- output
! zeig    r(kndm)    sqrt(eigenvalue)
! q     r(kndm,kndm) rotation matix
! qi    r(kndm,kndm) inverse of q
! xi    r(kndm,kndm) inverse of x
! imn       i        location of minimum eigenvalue
!--- prpc-parameter
! kndm      i        declared size of matrices
!--      
  use module_opt_visir,only:kndm,krowij
  implicit none

! input & output
  integer,intent(in):: iw0
  integer,intent(in):: n
  real(8),intent(in):: x(kndm,kndm)
  real(8),intent(in):: y(kndm,kndm)

  real(8),intent(out):: zeig(kndm)
  real(8),intent(out):: q(kndm,kndm)
  real(8),intent(out):: qi(kndm,kndm)
  real(8),intent(out):: xi(kndm,kndm)
  integer,intent(out):: imn
  logical,intent(inout) :: rt_failed !logic if radiative tarnsfer is failed, --> true (default is fault)

! working area
  integer:: k,i,j,ierr
  real(8):: v(kndm,kndm),sqx(kndm,kndm),sqxi(kndm,kndm)
  real(8):: rowij(krowij),wk(kndm),zmn
!--exec
 rt_failed = .false.  !default

! root decomposition of x (use v and zeing for u and xeig).
  k=0
  do i=1,n
     do j=1,i
        k=k+1
        rowij(k)=x(i,j)
     enddo
  enddo

  call symtrx(rowij,n,zeig,v,kndm,wk,ierr)

  if(ierr.gt.128) then
     print*, 'MSG getqm: error in decomposition of x in getqm'
     rt_failed = .true.  ; return
  endif

  do i=1,n
     if(zeig(i) <= 0.d0) then
        write(*,*) i,zeig(i)
        print*,'MSG getqm: non-positive eigenvalue of x'
        rt_failed = .true.  ; return
     endif
     zeig(i)=sqrt(zeig(i))
  enddo

  do i=1,n
     do j=1,n
        sqx (i,j)=sum(v(i,1:n)*zeig(1:n)*   v(j,1:n))
        sqxi(i,j)=sum(v(i,1:n)/zeig(1:n)*   v(j,1:n))
        xi  (i,j)=sum(v(i,1:n)/zeig(1:n)**2*v(j,1:n))
     enddo
  enddo

  v(1:n,1:n)=matmul(sqx(1:n,1:n),y(1:n,1:n))
  q(1:n,1:n)=matmul(v(1:n,1:n),sqx(1:n,1:n))

! root decomposition of z (use q for z).
  k=0
  do i=1,n
     do j=1,i
        k=k+1
        rowij(k)=q(i,j)
     enddo
  enddo

  call symtrx(rowij,n,zeig,v,kndm,wk,ierr)
  if(ierr > 128) then
     print*, 'MSG getqm: error in decomposition of z in getqm'
     rt_failed = .true.  ; return
  endif

! check minimum eigenvalue
  imn=1; zmn=zeig(1)
  if(n >= 2) then
     do j=2,n
        if(zeig(j) < zmn) then
           imn=j; zmn=zeig(j)
        endif
     enddo
  endif

! renormalization
  if(iw0==1) zeig(imn)=0.d0
  do i=1,n
     if(zeig(i) < 0.d0) then
        print*, 'MSG getqm: non-positive eigenvalue of z'
        rt_failed = .true.  ; return
     endif
     zeig(i)=sqrt(zeig(i))
  enddo

! q-matrices
  do i=1,n
     do j=1,n
        q (i,j)=sum(sqx(i,1:n)*v(1:n,j))
        qi(i,j)=sum(v(1:n,i)*sqxi(1:n,j))
     enddo
  enddo
  return
end subroutine getqm

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine grndl3(fsol,galb,bgnd,dpt,m,n1,am,w,n0,am0,r,t,er,et)
! lambert surface
!--- history
! 95. 5.26 generated from grndl introducing am and w
!            and r=2*galb... -> sumwm*galb...
! 08.10.26 modified fortran 90 version.
!--- input
! fsol       r           solar irradiance at the system top.
! galb       r           flux reflectivity of the surface.
! bgnd       r           thermal intensity from the surface=(1-r)*plank.
! dpt        r           optical depth at the surface.
! m          i           fourier order.
! n1         i           number of quadrature points.
! am      r(kndm)        mu(i)
! w       r(kndm)        w (i)
! n0         i           number of solar angles.  (== 1 for SDSU)
! am0     r(kna0)        cos(solar zenith angle)  (kna0 == 1 for SDSU)
!--- output
! r      r(kndm,kndm)    scaled refelction matrix.
! t      r(kndm,kndm)    scaled transmission matrix.      (always zero for surface)
! er     r(kndm,kna0)    scaled source matrix for upward radiance.    
! et     r(kndm,kna0)    scaled source matrix for downward radiance.  (always zero for surface) 
!--
  use module_opt_visir, only:kndm,kna0,pi
  implicit none
  integer,intent(in)::m,n1,n0
  real(8),intent(in)::fsol,galb,bgnd,dpt,am(kndm),w(kndm),am0(kna0)
  real(8),intent(out)::r(kndm,kndm),t(kndm,kndm),er(kndm,kna0),et(kndm,kna0)

  real(8)::sumwm,trns0,x
  integer::j
!--exec
  t(1:n1,1:n1)=0.d0  ! scaled transmission matrix is always zero at surface

  sumwm=sum(am(1:n1)*w(1:n1)) 
  
  if(galb <= 0.d0 .or. m > 0) then !in case of black surface or forier order /= 0
     r(1:n1,1:n1)=0.d0
  else
     do j=1,n1  !quadrature (= stream #) integration
        r (1:n1,j)=galb*sqrt(am(1:n1)*w(1:n1)*am(j)*w(j))/sumwm  ! galb is rescaled here (for 0th forier order)
     enddo
  endif

  et(1:n1,1:n0)=0.d0    ! scaled source matrix for downward radiance is alwasy zero at surface

  if(m > 0) then
     er(1:n1,1:n0)=0.d0  ! scaled source matrix for upward is zero for forier order /= 0
  else
     do j=1,n0  !n0 == 1 for SDSU
        trns0=exp(-dpt/am0(j))*fsol             ! atmos transmission  (fsol=0 for IR. So trns0=0 for IR) 
        x=galb*am0(j)*trns0/sumwm+2.d0*pi*bgnd  ! galb is here too
        er(1:n1,j)=sqrt(am(1:n1)*w(1:n1))*x     ! scaled source matrix for upward radiance
     enddo
  endif

  return
 end subroutine grndl3

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine grndo1(indg,fsol,u10,bgnd,dpt,m,nda,wmp,na0,am0,cr,ci, &
                   r,t,er,ert,et,amua,wa)

! water surface
!--- history
! 89.11. 2
! 93. 2.26  bug bg
!     3. 2  chenged calculation of emissivity
!     3.24  sub.ocnrf3 changed parameter ( add wa )
!     4. 5  put save
!     4. 5 nsf = 1 for initialization (terry)
!     5. 4 nsf -> indg; n1, n0 change name as nda, na0
!          if(nsf.gt.0) -> if(indg.eq.1)
! 95. 6. 2 generated from grndo with cr and ci
! 96. 1.12  ems(kna0)
! 97. 3.17 add ert by takashi nakajima to thermal. 
! 98. 2. 4  replace ems(kna0)->ems(kndm)
! 08.10.24 modified fortran 90 version.
! 09.05.18 Bug fixed (Reported by Hiro Masunaga)
! 09.09,04 Create LUT approach (matsui)
!--- input
! indg       i           1 then initialization of ocean surface matrix
!                        otherwise use the calculated matrices.
! fsol       r           solar irradiance at the system top.
! u10        r           wind velocity in m/s at 10m height
! bgnd       r           plank intensity from the surface
!                        emitted thermal is calculated as (1-r)*bgnd
!                        where r is unidirectional reflectivity.
!                        note the definition of bgnd is different
!                        from that for grndl.
! dpt        r           optical depth at the surface.
! m          i           fourier order.
! nda        i           number of quadrature points.
! wmp     r(kndm)        sqrt(wm)
! na0        i           number of solar angles.
! am0     r(kna0)        cos(solar zenith angle)  (kna0 == 1 for SDSU)
!                        needs all the m-fourier components
! cr         r           relative refractive index of the media
!                        about 1.33 for atmosphere to ocean incidence,
!                        and 1/1.33 for ocean to atmosphere incidence.
! ci         r           relative refractive index for imaginary part
!                        m = cr + i*ci
!--- output
! r      r(kndm,kndm)    scaled refelction matrix.
! t      r(kndm,kndm)    scaled transmission matrix.
! er     r(kndm,kna0)    scaled source matrix for upward radiance.
! et     r(kndm,kna0)    scaled source matrix for downward radiance.
!--
  use module_opt_visir, only: kndm,kna0,kndm2,pi,iw_global

  implicit none

! input
  integer,intent(in):: indg      !! surface flag (0:Lambert, 1: Ocean)
  real(8),intent(in):: fsol      !! solar irradiance at the system top
  real(8),intent(in):: u10       !! U10 [m/sec]
  real(8),intent(in):: bgnd      !! planck function at ground 
  real(8),intent(in):: dpt       !! 
  integer,intent(in):: m         !! Fourier order
  integer,intent(in):: nda           !! number of stream
  real(8),intent(in):: wmp(kndm)
  integer,intent(in):: na0           !! number of solar zenith angles (na0  == 1 for SDSU)
  real(8),intent(in):: am0(kna0)     !! cosine of solar zenith angles (kna0 == 1 for SDSU)
  real(8),intent(in):: cr           !! relative reflactive index (real part)
  real(8),intent(in):: ci           !! relative reflactive index (imaginary part)
  real(8),intent(in):: amua(kndm)   !! cosine of the zenith-quadrature angles. 1 to 0.
  real(8),intent(in):: wa(kndm)     !! corresponding quadrature weights.

! output
  real(8),intent(out):: r(kndm,kndm)  !!scaled refelction matrix
  real(8),intent(out):: t(kndm,kndm)  !! scaled transmission matrix (zero at surface)
  real(8),intent(out):: er(kndm,kna0)  !! scaled source matrix for upward radiance
  real(8),intent(out):: ert(kndm,kna0) !!
  real(8),intent(out):: et(kndm,kna0)  !! scaled source matrix for downward radiance (zero at surface)

!  real(8),save:: orf(kndm,kndm,kndm2)
!  real(8),save:: srr(kndm,kna0,kndm2)
  real(8),save,allocatable :: orf(:,:,:)
  real(8),save,allocatable :: srr(:,:,:)

  integer:: i,j,m1 !looping indice
  real(8):: bg     !background emission * 2pi
  real(8):: x,trns0,sr(kndm)
  real(8):: ems(kndm) !water emissivity (1-reflection)
  integer :: ierr  !allocation stat
  logical,parameter :: lut_grnd = .true.  ! true allow to use orf_srr LUT (5x faster)

!
! allocate orf & srr
!
  if( .not. allocated(orf) ) then
      allocate( orf(1:kndm,1:kndm,1:kndm2) , stat=ierr )
      if (ierr /= 0) stop 'MSG grndo1: allocation error -> Terminate program. '
  endif
  if( .not. allocated(srr) ) then
      allocate( srr(1:kndm,1:kna0,1:kndm2) , stat=ierr )
      if (ierr /= 0) stop 'MSG grndo1: allocation error -> Terminate program. '
  endif

!--exec
  m1=m+1


 if(lut_grnd .eqv. .true.) then

!
! get r and sr
!
     call interp_sr_r(m,nda,am0(1),u10,r,sr)

     orf(1:nda,1:nda,m1)=r(1:nda,1:nda)  !reflection matrix used for scaled reflection matrix and emissivity 
     srr(1:nda,1,m1)    =sr(1:nda)       !source matrix 

 else

!
! initialize orf and srr for indg==1
!
   if(indg == 1) then

     call ocnr11(m,nda,amua,wa,cr,ci,u10,r)
     orf(1:nda,1:nda,m1)=r(1:nda,1:nda)  !reflection matrix used for scaled reflection matrix and emissivity 

     do j=1,na0  !na0==1 for SDSU
        call ocnr31(m,nda,amua,wa,am0(j),cr,ci,u10,sr)  
        srr(1:nda,j,m1)=sr(1:nda)
     enddo
 
   endif

 endif

!
! scaled source matrix for downward radiance (zero at surface)
!
  t(1:nda,1:nda)=0.d0  

!
! scaled refelction matrix
!
  do j=1,nda !stream #
     r(1:nda,j)=wmp(1:nda)*orf(1:nda,j,m1)/wmp(j)  !scaled refelction matrix
  enddo

!
! emissivity and background emission
!
  if(m == 0)then !for 0th Fourier order (assuming Lambertian ??)
     bg=2.d0*pi*bgnd !background emission
     do i=1,nda !stream #
        ems(i)=1.d0-sum(orf(i,1:nda,1))  !emissivity of water
     enddo
  else
     bg=0.d0 !background emission is zero for high-order Fourier (assuming Lambertian??)
  endif

!
! scaled source matrix for downward radiance (zero at surface)
!
  et(1:nda,1:nda)=0.d0

!
! scaled source matrix for upward radiance
!
  do j=1,na0  !na0==1 for SDSU
     x=-dpt/am0(j)
     trns0=exp(x)*fsol  !tran0 is zero for IR, because fsol=0.
     do i=1,nda !stream #
        ert(i,j)=wmp(i)*(ems(i)*bg)
        er(i,j) =wmp(i)*(srr(i,j,m1)*trns0+ems(i)*bg)
     enddo
  enddo

  return
 end subroutine grndo1

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine interp_sr_r(m,nda,am0,u10,r,sr) 
 use module_opt_visir
 implicit none

! input
 integer,intent(in):: m        !! Fourier order
 integer,intent(in):: nda      !! number of stream
 real(8),intent(in):: am0      !! cosine of solar zenith angles (kna0 == 1 for SDSU)
 real(8),intent(in):: u10      !! U10 [m/sec]

! output
 real(8),intent(out):: r(kndm,kndm)  !!scaled refelction matrix
 real(8),intent(out):: sr(kndm)      !! scaled transmission matrix (zero at surface)

 integer :: m1
 real(8) :: wgt1_u10, wgt2_u10
 real(8) :: u10_filter  !Filetered U10 [m/sec]
 real(8) :: du10        !increment of u10 [m/s]
 real(8) :: wgt1_ia, wgt2_ia
 integer :: iu, iu_low, iu_hig
 integer :: ia_low, ia_hig

 m1 = m+1

!
! u10 filter
!
  du10 = 25.d0 / REAL(mxpts_u10)
  u10_filter = min( REAL(mxpts_u10) * du10 ,  max( REAL(1) * du10 ,  u10 ) )  

!
! find iu index
!
  iu_low = 0 ; iu_hig = 0
  iu_low = max( 1, min( mxpts_u10-1, INT( u10_filter / du10)) )
  iu_hig = iu_low+1
  wgt1_u10 = ( REAL(iu_hig)*du10 - u10_filter ) / du10
  wgt2_u10 = 1.d0 - wgt1_u10

!
! find r
!
  r(1:nda,1:nda) =  & !  orf_lut(1:nda,1:nda,iu,m1,iw)=r(1:nda,1:nda)  ! ORF LUT
    orf_lut(1:nda,1:nda,iu_low,m1,iw_global) * wgt1_u10 + &
    orf_lut(1:nda,1:nda,iu_hig,m1,iw_global) * wgt2_u10 


!
! find am0 index
!
  ia_low = min(99, INT(am0*100))
  ia_hig = ia_low + 1
  wgt1_ia = REAL(ia_hig)*1.d-2 - am0
  wgt2_ia = 1. - wgt1_ia
 
!
! find sr
!
  sr(1:nda) =  &
              (srr_lut(1:nda,ia_low,iu_low,m1,iw_global) * wgt1_u10 + &
               srr_lut(1:nda,ia_low,iu_hig,m1,iw_global) * wgt2_u10 ) * wgt1_ia &
             + &
              (srr_lut(1:nda,ia_hig,iu_low,m1,iw_global) * wgt1_u10 + &
               srr_lut(1:nda,ia_hig,iu_hig,m1,iw_global) * wgt2_u10 ) * wgt2_ia


 return
 end subroutine interp_sr_r

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function hf(tau,am1,am2,am3)
! geometrical factor for the scondery scattering eq.(24) of nt.
!    hf=integ(0,tau)dt*integ(0,t)dt1*exp(t*(1/mu1-1/mu2)
!             + t1*(1/mu2-1/mu3))*exp(-tau/am1)/am1/am2
!--- reference
! nt:  T. Nakajima and M. Tanaka, 1988, jqsrt, 40, 51-69
!--- history
! 88. 9.22  created by T. Nakajima
! 89. 5. 4  use expfn
! 08.10.24  modified fortran 90 version.
!--- input
! tau      r         optical thickness of the layer.
! am1      r         cos(zenith angle-1).
! am2      r         cos(zenith angle-2).
! am3      r         cos(zenith angle-3).
!--- output
! hf       f         geometrical factor.
!--
  implicit none

! input & output
  real(8),intent(in):: tau
  real(8),intent(in):: am1
  real(8),intent(in):: am2
  real(8),intent(in):: am3
  real(8):: hf

! work
  integer,save:: init = 1
  real(8),save:: eps
  real(8):: x1,x2,x3,ex1,ex2,ex3
  real(8):: ai1,ai2

! cpcon
  real(8)::ccp(3)
  real(8):: expfn
!--exec
! set eps: if abs(1/am1 - 1/am0)*tau  <=  eps then
!                      the routine sets almucantar condition-ialm.
  if(init > 0) then
     init=0
     call cpcon(ccp)
     eps=ccp(1)*30.d0
  endif

  x1=1.d0/am1-1.d0/am2; x2=1.d0/am2-1.d0/am3;  x3=1.d0/am1-1.d0/am3
  ex1=-tau/am1; ex2=-tau/am2; ex3=-tau/am3
  ex1=expfn(ex1); ex2=expfn(ex2); ex3=expfn(ex3)

  if(abs(x2*tau) > eps) then
! x2 <> 0
!! i1
     if(abs(x3*tau) <= eps) then
        ai1=ex1*(tau+x3*tau*tau/2.d0)
     else
        ai1=(ex3-ex1)/x3
     endif
!! i2
     if(abs(x1*tau) <= eps) then
        ai2=ex1*(tau+x1*tau*tau/2.d0)
     else
        ai2=(ex2-ex1)/x1
     endif
     hf=(ai1-ai2)/am1/am2/x2
     return
  endif

! x2 =  0
  if(abs(x1*tau) <= eps) then
     hf=tau**2*(0.5d0-x1*tau/3.d0)*ex1/am1/am2
  else
     hf=((tau-1.d0/x1)*ex2+ex1/x1)/am1/am2/x1
  endif
  return
end function hf

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine homog2(m,t1,t2,omg,n1,am1,wmm,n2,am2,pr,pt,pr0,pt0,fsol,nplk1,cplk, &
     r,t,er,et,zeig,q,qi,c11,c22,vp,vm,dp,dm, rt_failed)
! solve the transfer in a homogeneous scattering and emitting
!  medium by the discrete ordinate method.
!--- history
! 89.10.31 created from homog1 including thermal radiation.
! 95.11.22 eps=ccp(1)*10 -> eps=exp(log(ccp(1))*0.8)
! 08.10.26 modified fortran 90 version.
!--- input
! t1         r       optical depth at the layer top.
! t2         r       optical depth at the layer bottom.
! omg        r       single scattering albedo.
! n1         i       no. of the quadrature streams.
! am1     r(kndm)    mu (i), i=1,n1.
! wmm     r(kndm)    sqrt(w1/m1)
! n2         i       no. of the solar directions.
! am2     r(kna0)    mu0(i), i=1,n2.
! pr      r(kndm,    scaled p+-(i,j) i,j=1,n1
!           kndm)
! pt                 scaled p++(i,j)
! pr0     r(kndm,    scaled p0+-(i,j)  i=1,n1; j=1,n2.
! pt0       kndm)    scaled p0++(i,j)  i=1,n1; j=1,n2.
! fsol       r       solar irradiance at the top of the system.
! nplk1      i       max order of plank function expansion by tau + 1.
!                      if nplk1=0 then no thermal.
! cplk    r(nplk1)   2*pi*(1-w)*b(n)
!--- output
! r       r(kndm,    reflection   matrix   rij, i,j=1,n1.
! t         kndm)    transmission matrix   tij.
! er      r(kndm,    upgoing   source matrix eu(i,j), i=1,n1; j=1,n2.
! et        kna0)    downgoing source matrix ed(i,j).
! zeig    r(kndm)    root of the eigenvalues of z.
! q                  q-matrix
! qi                 inverse of q
! c11                c1 = inverse of 2a-
! c22                c2 = inverse of 2b-
! vp                 vs+
! vm                 vs-
! dp   r(kndm,kplk1) thermal emission intensity exapnsion (dnward)
! dm   r(kndm,kplk1) thermal emission intensity exapnsion (upward)
!--- parameter
! kndm       i       declared size for nda
! kna0       i       declared size for na0
! kplk1      i       declared size for nplk1
!--
  use module_opt_visir,only:kndm,kna0,kplk1,kndm2,krowij,pi,rad
  implicit none

! input
  integer,intent(in):: m
  real(8),intent(in):: t1
  real(8),intent(in):: t2
  real(8),intent(in):: omg 
  integer,intent(in):: n1 
  real(8),intent(in):: am1(kndm) 
  real(8),intent(in):: wmm(kndm) 
  integer,intent(in):: n2 
  real(8),intent(in):: am2(kna0)
  real(8),intent(in):: pr(kndm,kndm)
  real(8),intent(in):: pt(kndm,kndm)
  real(8),intent(in):: pr0(kndm,kna0)
  real(8),intent(in):: pt0(kndm,kna0)
  real(8),intent(in):: fsol
  integer,intent(in):: nplk1
  real(8),intent(in):: cplk(kplk1)

! output
  real(8),intent(out)::r(kndm,kndm)
  real(8),intent(out)::t(kndm,kndm)
  real(8),intent(out)::er(kndm,kna0)
  real(8),intent(out)::et(kndm,kna0)
  real(8),intent(out)::zeig(kndm)
  real(8),intent(out)::q(kndm,kndm)
  real(8),intent(out)::qi(kndm,kndm)
  real(8),intent(out)::c11(kndm,kndm)
  real(8),intent(out)::c22(kndm,kndm)
  real(8),intent(out)::vp(kndm,kna0)
  real(8),intent(out)::vm(kndm,kna0)
  real(8),intent(out)::dp(kndm,kplk1)
  real(8),intent(out)::dm(kndm,kplk1)
  logical,intent(inout) :: rt_failed !logic if radiative tarnsfer is failed, --> true (default is fault)

!--- working areas
  integer::iw0,ipk,imn,i,j,k,iw(kndm2)
  real(8)::ccp(3),eps,tau,x(kndm,kndm),y(kndm,kndm),xi(kndm,kndm), &
       cplk1,sum1,sp(kndm,kna0),sm(kndm,kna0),g1(kndm,kna0), &
       gam(kndm,kna0),e0(kna0),trns0,sum2,sum3,sum4,c(kndm), &
       sl(kndm),sl1(kndm),ap(kndm,kndm),am(kndm,kndm),bp(kndm,kndm), &
       bm(kndm,kndm),dt,dp0(kndm),dp1(kndm),dm0(kndm),dm1(kndm),vp0,vm1
  real(8)::expfn
!--exec
  rt_failed = .false.
! precision
  call cpcon(ccp)
  eps=exp(log(ccp(1))*0.8d0)

  tau=t2-t1
! x, y matrices
  do i=1,n1
     x(i,1:n1)=-omg*(pt(i,1:n1)-pr(i,1:n1))
     y(i,1:n1)=-omg*(pt(i,1:n1)+pr(i,1:n1))
     x(i,i)=1.d0/am1(i)+x(i,i)
     y(i,i)=1.d0/am1(i)+y(i,i)
  enddo

! decomposition of xy
  if(m==0 .and. 1.d0-omg <= eps) then
     iw0=1
  else
     iw0=0
  endif
  if(m==0 .and. iw0==0 .and. nplk1 > 0) then
     ipk=1
  else
     ipk=0
  endif

  call getqm(iw0,n1,x,y,zeig,q,qi,xi,imn,rt_failed)

  if( rt_failed ) return


! thermal source
!! c
  if(ipk==0) then
     if(nplk1 > 0) then
        dp(1:n1,1:nplk1)=0.d0
        dm(1:n1,1:nplk1)=0.d0
     endif
  else
     do i=1,n1
        sum1=sum(q(1:n1,i)*wmm(1:n1))
        do j=nplk1,1,-1
           if(j+2 > nplk1) then
              cplk1=0.d0
           else
              cplk1=dp(i,j+2)
           endif
           dp(i,j)=(dble((j+1)*j)*cplk1+sum1*cplk(j))/zeig(i)**2
        enddo
     enddo
!! d
     dm(1:n1,1:nplk1)=matmul(q(1:n1,1:n1),dp(1:n1,1:nplk1))
     do j=1,nplk1
        sum1=0.d0
        do i=1,n1
           if(j+1 <= nplk1) sum1=sum(qi(1:n1,i)*dp(1:n1,j+1))
           dp(i,j)=dm(i,j)-dble(j)*sum1
           dm(i,j)=dm(i,j)+dble(j)*sum1
        enddo
!            write(6,'(i1,4e15.7)') j,dp(1:n1,j),dm(1:n1,j)
     enddo
  endif

! sigma+ - (for single scattering)
  do i=1,n1
     sp(i,1:n2)=omg*(pt0(i,1:n2)+pr0(i,1:n2))
     sm(i,1:n2)=omg*(pt0(i,1:n2)-pr0(i,1:n2))
  enddo

! lower g
  do i=1,n1
     do j=1,n2
        g1(i,j)=-sum(x(i,1:n1)*sp(1:n1,j))-sm(i,j)/am2(j)
     enddo
  enddo

! gamma
  gam(1:n1,1:n2)=matmul(qi(1:n1,1:n1),g1(1:n1,1:n2))
  do i=1,n1
     gam(i,1:n2)=gam(i,1:n2)/(1.d0/am2(1:n2)**2-zeig(i)**2)
  enddo

! vs+ and -
  do j=1,n2
     trns0=-t1/am2(j)
     trns0=expfn(trns0)*fsol
     do i=1,n1
        sum1=sum(q (i,1:n1)*gam(1:n1,j))
        sum2=sum(qi(1:n1,i)*gam(1:n1,j)/am2(j)+xi(i,1:n1)*sm(1:n1,j))
        vp(i,j)=(sum1+sum2)/2.d0*trns0
        vm(i,j)=(sum1-sum2)/2.d0*trns0
     enddo
  enddo

! e0
  do j=1,n2
     e0(j)=expfn(-tau/am2(j))
  enddo

! base function  c(tau) and s(tau).
  do i=1,n1
     call csfn(tau,tau,zeig(i),c(i),sl(i),sl1(i))
  enddo
      
! a+-, b+-
  do i=1,n1
     do j=1,n1
        sum1=q (i,j)*c  (j); sum2=qi(j,i)*sl (j)
        sum3=q (i,j)*sl1(j); sum4=qi(j,i)*c  (j)
        ap(i,j)=sum1-sum2; am(i,j)=sum1+sum2
        bp(i,j)=sum3-sum4; bm(i,j)=sum3+sum4
     enddo
  enddo

! c11 and c22 -> their inversion.
  c11(1:n1,1:n1)=2.d0*am(1:n1,1:n1)
  c22(1:n1,1:n1)=2.d0*bm(1:n1,1:n1)

  call tnvss2(n1,c11,dt,0.d0,kndm,iw)
!  if(err /= '') then
!     err='error to get -c11- (homog2)'; return
!  endif
 
  call tnvss2(n1,c22,dt,0.d0,kndm,iw)
!  if(err /= '') then
!     err='error to get -c22- (homog2)'; return
!  endif

! r, t matrices
  do i=1,n1
     do j=1,n1
        sum1=sum(ap(i,1:n1)*c11(1:n1,j))
        sum2=sum(bp(i,1:n1)*c22(1:n1,j))
        r(i,j)=sum1+sum2
        t(i,j)=sum1-sum2
     enddo
  enddo

! er, et matrices
  dp1(1:n1)=0.d0; dm1(1:n1)=0.d0
  if(ipk==1) then
     do j=1,nplk1
        dp1(1:n1)=dp1(1:n1)+dp(1:n1,j)*tau**(j-1)
        dm1(1:n1)=dm1(1:n1)+dm(1:n1,j)*tau**(j-1)
     enddo
     dp0(1:n1)=dp(1:n1,1)
     dm0(1:n1)=dm(1:n1,1)
  else
     dp0(1:n1)=0.d0; dm0(1:n1)=0.d0
  endif

  do j=1,n2
     do i=1,n1
        sum1=0.d0; sum2=0.d0
        do k=1,n1
           vp0=vp(k,j)      +dp0(k)
           vm1=vm(k,j)*e0(j)+dm1(k)
           sum1=sum1+r(i,k)*vp0+t(i,k)*vm1
           sum2=sum2+t(i,k)*vp0+r(i,k)*vm1
        enddo
        er(i,j)=vm(i,j)      +dm0(i)-sum1
        et(i,j)=vp(i,j)*e0(j)+dp1(i)-sum2
     enddo
  enddo
  return
end subroutine homog2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine intcr1(imthd,am0,am1u,cs1,ntau,utau,utaut,nln,thk,thkt, &
     omg,omgt,phs,phst,ff,mxlgn1,nlgn1,nlgt1,g,epsp,nchk1,cor, &
     sgl2,phsb)
! sigle and second scattering correction
! assume delta-m method
!--- history
! 90. 1.28 created
! 08.10.20 modified fortran 90 version.
!--- input
! imthd      i         1: ms, 2: tms, 3: ims
! am0        r         cos(solar zenith angle)
! am1u       r         cos(emergent nadir angle)
!                      < o: upward,   > 0: downward
! cs1        r         cos(scattering angle)
! ntau       i         number of user defined angle
! utau     r(ntau)     user defined optical depth
! utaut    r(ntau)     truncated user defined depth
! nln        i         number of layer
! thk      r(nln)      optical thickness of sublayer
! thkt     r(nln)      truncated thickness
! omg      r(nln)      single scattering alebedo
! omgt     r(nln)      truncated single scattering albedo
! phs      r(nln)      phase function  (integral = 1 over unit sphere)
! phst     r(nln)      truncated phase function (1 over unit sphere)
! ff       r(nln)      truncation fraction
! klgn1      i         first argument size of g
! mxlgn1     i         max(nlgn1) (only for imthd=3)
! nlgn1    i(nln)      max order of legender series +1 in each sublayer
!                        (only for imthd=3)
! nlgt1    i(nln)      same as -nlgn1- but for truncation
! g      r(klgn1,nln)  phase function moment (only for imthd=3)
! epsp       r         convergence criterion of legendre seris of
!                        phase function ** 2
! nchk1      i         number of consective convergences before
!                        final decision
!--- output
! cor      r(ntau)     ai= ai+cor*fsol is corrected intensity
!                        at each user
!                       defined depth (utau)
!--- work
! sgl2     r(ntau)
! phsb     r(nln)
!--- area for this routine
  use module_opt_visir, only: kntau,knln,klgn1,pi,rad
  implicit none

! input
  integer,intent(in):: imthd
  real(8),intent(in):: am0
  real(8),intent(in):: am1u
  real(8),intent(in):: cs1
  integer,intent(in):: ntau
  real(8),intent(in):: utau(kntau)
  real(8),intent(in):: utaut(kntau)
  integer,intent(in):: nln
  real(8),intent(in):: thk(knln)
  real(8),intent(in):: thkt(knln)
  real(8),intent(in):: omg(knln)
  real(8),intent(in):: omgt(knln)
  real(8),intent(in):: phs(knln)
  real(8),intent(in):: phst(knln)
  real(8),intent(in):: ff(knln)
  integer,intent(in):: mxlgn1
  integer,intent(in):: nlgn1(knln)
  integer,intent(in):: nlgt1(knln)
  real(8),intent(in):: g(klgn1,knln)
  real(8),intent(in):: epsp
  integer,intent(in):: nchk1

! output
  real(8),intent(out):: cor(kntau)
  real(8),intent(out):: sgl2(kntau)
  real(8),intent(out):: phsb(knln)

! work
  integer:: l,it,k1
  real(8):: sgl1
  real(8):: utaus,eh,sh,shh,phspk,dpt2,dpt1
  real(8):: wh,fh,am3,phspk2,gpk,tau,gp,gpk2, sgl3
  integer:: init,ichk

! function
  real(8):: sglr,plgd,hf
!--exec
  cor(1:ntau)=0.d0
  if(imthd <= 0 .or. imthd >= 4) return

!--- ms-method in reference-nt
  if(imthd == 1) then
! sgl2 = u1wave of eq.(14)
     phsb(1:nln)=(1.d0-ff(1:nln))*phst(1:nln)
     do it=1,ntau
        sgl2(it)=sglr(am1u,am0,knln,thk,omg,phsb,utau(it))
     enddo

! correction by ms-method.  see eq.(14)
     do it=1,ntau
!! sgl1 = u1 of eq.(14)
        sgl1=sglr(am1u,am0,knln,thk,omg,phs,utau(it))
        cor(it)=-sgl2(it)+sgl1
     enddo
     return
  endif

!--- TMS and IMS methods
!  sgl2 = u1* of eq.(15)
  do it=1,ntau
     sgl2(it)=sglr(am1u,am0,knln,thkt,omgt,phst,utaut(it))
  enddo
  phsb(1:nln)=phs(1:nln)/(1.d0-ff(1:nln))

! correction of intensity by tms-method.  see eq.(15)
  do it=1,ntau
!!  sgl1 = u1wave* of eq.(15)
     sgl1=sglr(am1u,am0,knln,thkt,omgt,phsb,utaut(it))
     cor(it)=-sgl2(it)+sgl1
  enddo

  if(imthd == 2 .or. am1u <= 0.0) return

!--- secondary scattering correction for IMS-method
  do it=1,ntau
     utaus=utau(it)
! getting mean optical constants above the user defined level-utaus
     eh=0.d0; sh=0.d0
     shh=0.d0; phspk=0.d0
     dpt2=0.d0
     do l=1,nln
        dpt1=dpt2
        dpt2=dpt1+thk(l)
        if(utaus <= dpt1) exit
        if(utaus < dpt2) then
           tau=utaus-dpt1
        else
           tau=thk(l)
        endif
        eh= eh     +       tau
        sh= sh     +omg(l)*tau
        shh=shh    +omg(l)*tau*ff(l)
        phspk=phspk+omg(l)*tau*(phs(l)-(1-ff(l))*phst(l))
     enddo
     if(abs(eh) <= 0.d0) cycle

! wh: mean single scattering albedo
     wh=sh/eh

! fh: mean truncated fraction
     fh=shh/sh
     if(fh <= 0.d0) cycle

! phspk: mean trunceted peak of phase function
     phspk=phspk/shh

! am3: variable appearing in eq.(23)
     am3=am0/(1.d0-fh*wh)

! Legendre sum for eq.(23)
!  we truncate the series when abs(gpk**2) becomes smaller than -epsp-
!   for successive three times.
     phspk2=0.d0
     init=1
     ichk=0
     do k1=1,mxlgn1
!! phspk: 2phat-phat2 of eq.(23)
!!  mean value for the layer above the user-defined level
        gpk=0.d0; dpt2=0.d0
        do l=1,nln
           dpt1=dpt2
           dpt2=dpt1+thk(l)
           if(k1 <= nlgn1(l) .and. omg(l) > 0.d0) then
              if(utaus <= dpt1) exit
              if(utaus < dpt2) then
                 tau=utaus-dpt1
              else
                 tau=thk(l)
              endif
              if(k1 <= nlgt1(l)) then
                 gp=ff(l)
              else
                 gp=g(k1,l)
              endif
              gpk=gpk+gp*omg(l)*tau
           endif
        enddo
        gpk=gpk/shh
        gpk2=gpk**2
        if(gpk2 <= epsp) then
           ichk=ichk+1
           if(ichk >= nchk1) exit
        endif
        phspk2=phspk2+dble(2*k1-1)*gpk2*plgd(init,cs1)
     enddo
     phspk2=phspk2/4.d0/pi
!! correction by uu = uu + uhat in ims-method
!!   sgl3: uhat of eq.(23)
     sgl3=(fh*wh)**2/(1.d0-fh*wh)*(2.d0*phspk-phspk2) &
          * hf(utaus, am1u, am3, am3)
     cor(it)=cor(it)-sgl3
  enddo
  return
end subroutine intcr1

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine lgndf3(lmax1,n,x,y,g)
! Legendre expansion (same as lgndf2 but generating g-moments).
!--- history
! 90. 1.20 created from lgndf2, use expfn.
!          direction of integration from x(1) to x(n).
! 91. 2.16 strip ng from save statement.
! 92. 4. 3 kill the statement of gw=gw/2 after qgausn
!     6.22 add gw/2 again because the above change is mistake
! 08.10.20 modified fortran 90 version.
!--- input
! lmax1      i      maximum order + 1.
! n          i      number of data on (-1, 1). .ge. 4.
! x        r(na)    independent variables on (-1, 1)
! y        r(na)    y(x).
!--- output
! g    r(lmax1)     y = sum(l1=1,lmax1) (2*l1-1)*g(l1)*pl(l1)
!                      where pl(l1) is (l1-1)th order legendre
!                      polynomial.
!--
  use module_opt_visir, only: pi,rad
  implicit none
  integer,intent(in)::lmax1,n
  real(8),intent(in)::x(n),y(n)
  real(8),intent(out)::g(lmax1)
! working areas.
  integer,save::init=1
  integer,parameter::ng=5
  real(8),save::gw(ng),gx(ng)
  integer::i,j,i1,i2,i3,i4,isign,l1
  real(8)::x1,x2,x3,x4,alp1,alp2,alp3,alp4,xx,ww,pp,pl,pl1,pl2
  real(8)::expfn
!--exec
! shifted gaussian quadrature.
  if(init >= 1) then
     init=0
     call qgausn(gw,gx,ng)
     gw(1:ng)=gw(1:ng)/2.0d0
  endif
! clear
  g(1:lmax1)=0.d0
! loop for angle
  do i=1,n-1
! cubic interpolation
     if(i <= 2) then
        i1=1; i4=4
     else
        if(i <= n-2) then
           i1=i-1; i4=i+2
        else
           i1=n-3; i4=n
        endif
     endif

     i2=i1+1; i3=i2+1
     x1=x(i1); x2=x(i2); x3=x(i3); x4=x(i4)
     if( (y(i1) <= 0.) .or. (y(i2) <= 0.) .or. (y(i3) <= 0.) .or.  &
          (y(i4) <= 0.) ) then
        isign=-1; alp1=y(i1); alp2=y(i2); alp3=y(i3); alp4=y(i4)
     else
        isign=1;  alp1=log(y(i1)); alp2=log(y(i2))
        alp3=log(y(i3)); alp4=log(y(i4))
     endif

! loop for gaussian integration
     do j=1,ng
!! interpolated value of y
        xx=x(i)+gx(j)*(x(i+1)-x(i))
        ww=gw(j)*(x(i+1)-x(i))
        pp =  (xx-x2)*(xx-x3)*(xx-x4)/(x1-x2)/(x1-x3)/(x1-x4)*alp1 &
             +(xx-x1)*(xx-x3)*(xx-x4)/(x2-x1)/(x2-x3)/(x2-x4)*alp2 &
             +(xx-x1)*(xx-x2)*(xx-x4)/(x3-x1)/(x3-x2)/(x3-x4)*alp3 &
             +(xx-x1)*(xx-x2)*(xx-x3)/(x4-x1)/(x4-x2)/(x4-x3)*alp4
        if(isign == 1) pp=expfn(pp)

! legendre sum
        pl1=0.d0; pl=1.d0
        g(1)=g(1)+pp*ww
        if(lmax1 >= 2) then
           do l1=2,lmax1
              pl2=pl1; pl1=pl
              pl=(dble(2*l1-3)*xx*pl1-dble(l1-2)*pl2)/dble(l1-1)
              g(l1)=g(l1)+pp*pl*ww
           enddo
        endif
     enddo
  enddo
  return
end subroutine lgndf3

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine multi(nd,cd,cc)
! calculation of multiple reflection between two layers.
!--- history
! 86.10.15  check ok.
! 89.10.30  re-edit.
! 94. 5. 7  err*64 -> err*(*)
! 08.10.27 modified fortran 90 version.
!--- input
! nd          i     no. of streams at the interface between two layers.
! cd      r(kndm,   r = r1 * r2.
!           kndm)
!--- output
! cc      r(kndm,   cc = ( 1 - cd )**-1 = 1 + cd + cd**2 + cd**3 + ...
!           kndm)
!--- parameter
! kndm        i     number of quadrature points.
!--- areas for this routine
  use module_opt_visir, only: kndm,kndm2
  implicit none

! input & output
  integer,intent(in)::nd
  real(8),intent(in)::cd(kndm,kndm)
  real(8),intent(out)::cc(kndm,kndm)

! working area
  integer::iw(kndm2),i,j
  real(8)::eps,dt
!---exec
  do i=1,nd
     cc(i,1:nd)= -cd(i,1:nd)
     cc(i,i)=1.d0 - cd(i,i)
  enddo

! inversion of -cc-.
  eps=0.d0
  call tnvss2(nd,cc,dt,eps,kndm,iw)
!  if(err /= '') err='error in multi'
  return
end subroutine multi

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine ocnr11(m,nda,amua,wa,cr,ci,u10,r)
! reflection matrix of ocean surface
!--- history
! 92. 9. 1 created by hasumi
!    12.23 modified by nakajima
! 95. 6. 2 generated from ocnrf1
! 08.10.27 modified fortran 90 version.
!--- input
! m      i     fourier order
! kndm   i     declared dimension of r and amua
! nda    i     used dimension of r and amua
! amua  r(nda) quadrature points in hemisphere
!              decreasing order (zenith to horizon, or, 1 -> 0)
! wa    r(nda) quadrature weights
! em     r     relative refractive index
!              about 1.33 from atmosphere to ocean
!              about 1/1.33 from ocean to atmosphere
! u10    r     wind velocity at 10 m above the surface
! cr       r    relative refractive index of the media
!               about 1.33 for atmosphere to ocean incidence,
!               and 1/1.33 for ocean to atmosphere incidence.
! ci       r    relative refractive index for imaginary part
!               m = cr + i*ci
!--- output
! r    r(knda,nda)  reflection matrix for m-th fourier order
!                    ur = r * ui
!--
  use module_opt_visir, only: kndm,pi
  implicit none

! input & output
  integer,intent(in):: m
  integer,intent(in):: nda
  real(8),intent(in):: amua(kndm)
  real(8),intent(in):: wa(kndm)
  real(8),intent(in):: cr
  real(8),intent(in):: ci
  real(8),intent(in):: u10
  real(8),intent(out):: r(kndm,kndm)

! local parameter
  integer,parameter:: kn=30
  integer,parameter:: kndm1=kndm+1

! local variables
  integer,save:: init = 1
  integer,save:: n
  real(8),save:: gx(kn),gw(kn)

  real(8):: sigma,dispa,cosds
  real(8):: dfi2,fi1(kn),fi2(kn),cosm1(kn),cosm2(kn)

  real(8):: x(kndm1)
  integer:: i,ns1,ns(2)
  real(8):: b,x1,x2,xb(5,2)
  integer:: j,ieq,is,ii,k
  real(8):: rij,ami,dx,ame,w,rr(kndm,kndm),rrr

! for function
  real(8):: searf1
!--exec
  if(init > 0) then
     n=kn
     init=0
     call qgausn(gw, gx, n )
  endif

! parameters for integration
  sigma=sqrt(0.00534d0*u10)
  dispa=atan(sigma)
  cosds=cos(dispa)

  do k=1,n
     dfi2=pi-dispa
     fi1(k)=dispa+dfi2*gx(k)
     fi2(k)=dispa*gx(k)
     cosm1(k)=2.d0*cos(fi1(k)*dble(m))*dfi2 *gw(k)
     cosm2(k)=2.d0*cos(fi2(k)*dble(m))*dispa*gw(k)
  enddo

  x(1)=1.d0
  if(nda >= 2) then
     x(2:nda)=(amua(1:nda-1)+amua(2:nda))/2.d0
  endif
  x(nda+1)=0.d0

  do i = 1, nda
!! x1=cos(theta(i)-dispa), x2=cos(theta(i)+dispa)
!!   where amua(i)=cos(theta(i))
     b=sqrt((1.d0-amua(i)**2)*(1.d0-cosds**2))
     x1=amua(i)*cosds+b
     x2=amua(i)*cosds-b
!! seting mu-baoudary for mu-inegration
     ns1=1
     xb(1,1)=x(i)
     if(x1 < x(i)) then
        ns1=ns1+1
        xb(ns1,1)=x1
     endif
     ns1=ns1+1
     xb(ns1,1)=amua(i)
     if(x2 > x(i+1)) then
        ns1=ns1+1
        xb(ns1,1)=x2
     endif
     ns1=ns1+1
     xb(ns1,1)=x(i+1)
     ns(1)=ns1-1
     ns(2)=1
     xb(1,2)=x(i)
     xb(2,2)=x(i+1)

     do j = 1, nda
        if(i ==j) then
           ieq=1
        else
           ieq=2
        endif
!! mu-integration
        rij=0.d0
        ami = amua(j)
        do is=1,ns(ieq)
           dx=xb(is,ieq)-xb(is+1,ieq)
           do ii = 1, n
              ame=xb(is+1,ieq)+dx*gx(ii)
              w=dx*gw(ii)*ame
              do k = 1, n
                 rij=rij+(cosm1(k)*searf1(ame,ami,fi1(k),u10,cr,ci) &
                      +cosm2(k)*searf1(ame,ami,fi2(k),u10,cr,ci))*w
              enddo
           enddo
        enddo

        rr(i,j)=rij/wa(i)
     enddo
  enddo

! symmetric operation
  do i=1,nda
     do j=1,i
        rrr=(rr(i,j)+rr(j,i))/2.d0
        r(i,j)=rrr/amua(i)*wa(j)
        r(j,i)=rrr/amua(j)*wa(i)
     enddo
  enddo
  return
end subroutine ocnr11

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine ocnr31(m,nda,amua,wa,am0,cr,ci,u10,sr)
! reflection matrix of ocean surface
!--- history
! 92. 9. 1 created by hasumi
!    12.23 modified by nakajima
! 93. 3.22 /wa(i) debugged by takashi
!     3.29 ami -> am0
!--- input
! m      i     fourier order
! kndm   i     declared dimension of r and amua
! nda    i     used dimension of r and amua
! amua  r(nda) quadrature points in hemisphere
!              decreasing order (zenith to horizon, or, 1 -> 0)
! wa    r(nda) quadrature weights
! am0    r     cos (solar zenith angle)
! cr       r    relative refractive index of the media
!               about 1.33 for atmosphere to ocean incidence,
!               and 1/1.33 for ocean to atmosphere incidence.
! ci       r    relative refractive index for imaginary part
!               m = cr + i*ci
! u10    r     wind velocity at 10 m above the surface
!--- output
! sr    sr(knda)  reflection source matrix for m-th fourier order
!--
  use module_opt_visir, only: kndm,pi
  implicit none

! input & output
  integer,intent(in):: m
  integer,intent(in):: nda
  real(8),intent(in):: amua(kndm)
  real(8),intent(in):: wa(kndm)
  real(8),intent(in):: am0
  real(8),intent(in):: cr
  real(8),intent(in):: ci
  real(8),intent(in):: u10
  real(8),intent(out):: sr(kndm)

! local parameter
  integer,parameter:: kn=30
  integer,parameter:: kndm1=kndm+1

! local variables
  integer,save:: init = 1
  integer,save:: n
  real(8),save:: gx(kn),gw(kn)

  real(8):: sigma,dispa,cosds
  real(8):: dfi2,fi1(kn),fi2(kn),cosm1(kn),cosm2(kn)

  real(8):: x(kndm1)
  integer:: i,ns1,ns(2)
  real(8):: b,x1,x2,xb(5,2)
  integer:: j,ieq,is,ii,k
  real(8):: rij,dx,ame,w

! for function
  real(8):: searf1
!--exec
  if(init > 0) then
     n=kn
     init=0
     call qgausn(gw, gx, n )
  endif

! parameters for integration
  sigma=sqrt(0.00534d0*u10)
  dispa=atan(sigma)
  cosds=cos(dispa)

  do k=1,n
     dfi2=pi-dispa
     fi1(k)=dispa+dfi2*gx(k)
     fi2(k)=dispa*gx(k)
     cosm1(k)=2.d0*cos(fi1(k)*dble(m))*dfi2 *gw(k)
     cosm2(k)=2.d0*cos(fi2(k)*dble(m))*dispa*gw(k)
  enddo

  x(1)=1.d0
  if(nda >= 2) then
     x(2:nda)=(amua(1:nda-1)+amua(2:nda))/2.d0
  endif
  x(nda+1)=0.d0

  do i = 1, nda
!! x1=cos(theta(i)-dispa), x2=cos(theta(i)+dispa)
!!   where amua(i)=cos(theta(i))
     b=sqrt((1.d0-amua(i)**2)*(1.d0-cosds**2))
     x1=amua(i)*cosds+b
     x2=amua(i)*cosds-b
!! seting mu-baoudary for mu-inegration
     ns1=1
     xb(1,1)=x(i)
     if(x1 < x(i)) then
        ns1=ns1+1
        xb(ns1,1)=x1
     endif
     ns1=ns1+1
     xb(ns1,1)=amua(i)
     if(x2 > x(i+1)) then
        ns1=ns1+1
        xb(ns1,1)=x2
     endif
     ns1=ns1+1
     xb(ns1,1)=x(i+1)
     ns(1)=ns1-1
     ns(2)=1
     xb(1,2)=x(i)
     xb(2,2)=x(i+1)

     if(amua(i) == am0) then
        ieq=1
     else
        ieq=2
     endif

!! mu-integration
     rij=0.d0
     do is=1,ns(ieq)
        dx=xb(is,ieq)-xb(is+1,ieq)
        do ii = 1, n
           ame=xb(is+1,ieq)+dx*gx(ii)
           w=dx*gw(ii)*ame
           do k = 1, n
              rij=rij+(cosm1(k)*searf1(ame,am0,fi1(k),u10,cr,ci) &
                   +cosm2(k)*searf1(ame,am0,fi2(k),u10,cr,ci))*w
           enddo
        enddo
     enddo
     sr(i)=rij/amua(i)/wa(i)
  enddo
  return
end subroutine ocnr31

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine phas2(m1,mmax1,n1,n2,knp,kn1,kn2,g,pl1,pl2,pt,pr)
!--- history
! 90. 1.20  generated from phase (same as phase but using the moments g)
! 08.10.24  modified fortran 90 version.
! 08.11. 3  Fixed some dimensions.
!--- input
! m1         i       m + 1    fourier order + 1.
! mmax1      i       mmax+1   max fourier order + 1.
! n1         i       no. of emergent zenith angles.
! n2         i       no. of incident zenith angles.
! knp        i       size of n1 for array pt and pr.
! kn1        i       size of n1 for array pl1.
! kn2        i       size of n2 for array pl2.
! g     r(mmax1)     legendre moments of phase function.
!                      g(1)=1
! pl1   r(kn1,mmax1) legendre polynomials for emergent direction.
!                    (i,m1) = (emergent directions, order+1).
! pl2   r(kn2,mmax1) legendre polynomials for incident direction.
!--- output
! pt    r(knp,n2)    phase matrix for transmission.
! pr    r(knp,n2)    phase matrix for reflection.
!--
  use module_opt_visir,only: kndm2,pi
  implicit none
  integer,intent(in)::m1,mmax1,n1,n2,knp,kn1,kn2
  real(8),intent(in)::g(kndm2),pl1(kn1,kndm2),pl2(kn2,kndm2)
  real(8),intent(out)::pt(knp,kn2),pr(knp,kn2)
  integer::i,j,k1
  real(8)::sign,c4
!--      
  pt(1:n1,1:n2)=0.d0; pr(1:n1,1:n2)=0.d0
  if(m1 <= mmax1) then
     do j=1,n2
        do i=1,n1
           sign=-1.d0
           do k1=m1,mmax1
              sign=-sign
              c4=dble(2*k1-1)/2.d0*g(k1)*pl1(i,k1)*pl2(j,k1)
              pt(i,j)=pt(i,j)+c4
              pr(i,j)=pr(i,j)+sign*c4
           enddo
        enddo
     enddo
  endif
  return
end subroutine phas2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function pint4(init,ang1,na,l,ang,p)
! interpolation of the phase function.
!--- history
! 89.11. 8 created from pint3. change interpolation x-> ang
! 90. 1.23 debug
! 08.10.20 modified fortran 90 version.
!--- input
! init       i      1 then search ang1-interval else not search.
! ang1       r      scattering angle in degree for interpolation.
! na         i      no. of scattering angles.
! ang      r(na)    scattering angles in degree.
! p     r(kna,l)    phase function
! l          i      layer number
!--- output variables
! init       i      0
! pint4      r      interpolated value of p at ang1.
!--- variables for the routine.
  use module_opt_visir, only: knang, knln, pi, rad
  implicit none

! input & output  
  integer,intent(inout):: init
  real(8),intent(in):: ang1
  integer,intent(in):: na
  integer,intent(in):: l
  real(8),intent(in):: ang(knang)
  real(8),intent(in):: p(knang,knln)
  real(8):: pint4

! work
  integer,save:: i1,i2,i3
  real(8):: alp1,alp2,alp3
  real(8):: x1,x2,x3,xx
  integer:: i,isign
  real(8):: pp,expfn
!--exec
! search near angles  
  if(init >= 1) then
     init=0
     do i=1,na-1
        if((ang1-ang(i))*(ang1-ang(i+1)) <= 0.d0) then
           if(i<=1) then
              i1=1; i3=3
           elseif(i<na-1) then
              i1=i-1; i3=i+1
           else
              i1=na-2; i3=na
           endif
           i2=i1+1
           exit
        endif
     enddo
     if(i==na) then
        i1=na-2; i2=na-1; i3=na
     endif
  endif
  xx=ang1
  x1=ang(i1)
  x2=ang(i2)
  x3=ang(i3)
  alp1=p(i1,l)
  alp2=p(i2,l)
  alp3=p(i3,l)
  isign=-1
  if(alp1 > 0.d0 .and. alp2 > 0.d0 .and. alp3 > 0.d0) then
     isign=1
     alp1=log(alp1)
     alp2=log(alp2)
     alp3=log(alp3)
  endif
  pp=(xx-x2)*(xx-x3)/(x1-x2)/(x1-x3)*alp1 &
       +(xx-x1)*(xx-x3)/(x2-x1)/(x2-x3)*alp2 &
       +(xx-x1)*(xx-x2)/(x3-x1)/(x3-x2)*alp3
  if(isign >= 1) pp=expfn(pp)
  pint4=pp
  return
end function pint4

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine pkint(am1u,dpt,t1,t2,nplk1,c)
! integral(from t1 to t2) dt exp(-(dpt-t)/am1u) t**n
!--- history
! 90. 1.13 created
! 08.10.31 modified fortran 90 version.
!--- input
! am1u     r      userdefined emergent mu.
!                 if .gt. 0 then down,  if .lt. 0 then upward.
! dpt      r      optical depth for interpolation.
! t1       r      lower limit of integration.
! t2       r      upper limit of integration.
! nplk1    i      .ge. 1
!--- output
! c     r(nplk1)  integrations
!--
  implicit none

! input & output
  real(8),intent(in):: am1u
  real(8),intent(in):: dpt
  real(8),intent(in):: t1
  real(8),intent(in):: t2
  integer,intent(in)::nplk1
  real(8),intent(out)::c(nplk1)

! work
  real(8)::ex,x,ex1,ex2,expfn,exx
  integer::n1
!--exec  
  ex=-dpt/am1u; x=1/am1u
  c(1)=exx(ex,t1,t2,x)
  if(nplk1 >= 2) then
     ex=(t1-dpt)/am1u; ex1=expfn(ex)
     ex=(t2-dpt)/am1u; ex2=expfn(ex)
     do n1=2,nplk1
        ex1=t1*ex1; ex2=t2*ex2
        c(n1)=am1u*(ex2-ex1-dble(n1-1)*c(n1-1))
     enddo
  endif
  return
end subroutine pkint

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine quada(nda,amua,wa)
! construction of quadrature in the atmosphere.
!------ history
! 1986.10.2  use qgausn from Stamnes and Tsay.
! 2008.10.6  Rewrite fortran 90 free form.
!--- input
! nda        i      no. of quadrature streams in the hemisphere of atmos
!--- output
! amua    r(knda)   mua(i), i=1,n1    mua(1) > mua(2) > ...
! wa      r(knda)   corresponding weights.
!--- variables for the routine
  use module_opt_visir, only:kndm
  implicit none
  integer,intent(in)::nda
  real(8),intent(out)::amua(kndm),wa(kndm)
  integer::i,i1
  real(8)::x
!--exec
! shifted gaussian quadrature on (0, 1).
  call qgausn(wa,amua,nda)
! reordering.
  amua(1:nda)=amua(nda:1:-1)
  wa(1:nda)=wa(nda:1:-1)
  return
end subroutine quada


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function searf1(amue,amui,phi,u10,cr,ci)

! bidirectional reflectivity of ocean surface
! shadowing factor for single scattering
! so that energy conservation is not satisfied
!--- history
! 92. 9. 1 created by Hasumi
!    12.23 modified by Nakajima
! 95. 6. 2 generated from searef with cr and ci
! 00. 2. 3 debuug parameter for frnlr.
!--- input
! amue     r    cosine of zenith angle of emergent ray .gt. 0
! amui     r    cosine of zenith angle of incident ray .gt. 0
! phi      r    azimurthal angle  (radian)
! u10      r    wind velocity in m/s at 10m height
! cr       r    relative refractive index of the media
!               about 1.33 for atmosphere to ocean incidence,
!               and 1/1.33 for ocean to atmosphere incidence.
! ci       r    relative refractive index for imaginary part
!               m = cr + i*ci
!--- output
! searef  rf    bidirectional reflection function
!--
  use module_opt_visir, only: pi,rad
  implicit none

! input & output
  real(8),intent(in):: amue !cos (emergent ray)
  real(8),intent(in):: amui !cos (solar zenith angle)
  real(8),intent(in):: phi
  real(8),intent(in):: u10
  real(8),intent(in):: cr
  real(8),intent(in):: ci
  real(8):: searf1

! work
  real(8):: sqrtpi
  real(8):: amue1,alpha,cosw,amun
  real(8):: sigma,vie,ve,fe,vii,vi,fi,g
  real(8):: r,p

! function
  real(8):: erfc

! for frnlr
  real(8):: r1,r2,wi,wr
!--exec

  if( amue == 0.d0 .or. amui == 0.d0 ) then
    searf1 = 0.d0
    return
  endif


  sqrtpi=sqrt(pi)

  amue1=max(amue,1.0d-7)
  alpha=-amue1*amui+sqrt((1.d0-amue1**2)*(1.d0-amui**2))*cos(phi)
  cosw=sqrt((1-alpha)/2.d0)
  amun=(amue1+amui)/2.d0/cosw

! shadowing factor
  sigma=sqrt(5.34d-3*u10)
  vie =sigma*sqrt(1.d0-amue1**2)/amue1
  if (vie > 0.15d0) then
     ve =1.d0/vie
     fe=(exp(- ve**2)/sqrtpi/ve-erfc(ve))/2.d0
  else
     fe=0.d0
  endif

  vii =sigma*sqrt(1.d0-amui**2)/amui

  if (vii > 0.15d0) then
     vi =1.d0/vii
     fi=(exp(- vi**2)/sqrtpi/vi-erfc(vi))/2.d0
  else
     fi=0.d0
  endif
  g=1.d0/(1.d0+fe+fi)

! 00.2.3 debug (T.Y.Nakajima. By report of Norman G. Loeb (NASA LaRC))
  wi=acos(cosw)/rad
!  call frnlr(cr,ci,cosw,coswt,r1,r2)
  call frnlr(cr,ci,wi,wr,r1,r2)

  r=(r1+r2)/2.d0
  p=1.d0/pi/sigma**2/amun**3*exp(-(1.d0-amun**2)/(sigma*amun)**2)
  searf1=1.d0/4.d0/amue1/amun*g*p*r
  return
 end function searf1

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

function sglr (am1, am0, nl, tc, w, p, ut)
! singly scattered intensity (multi-layer) by eq. (12) of nt.
!--- reference
! nt:  T. Nakajima and M. Tanaka, 1988, jqsrt, 40, 51-69
!--- history
! 89. 2.22   created by T. Nakajima
!              from aids by changing the meaning of t, w, p,
!              and ut.  ut is a level inside the multi-layer system.
!     5. 4   patch for exponential overflow.
!     6.20   reform style.
!    12.20   with cpex(3)
! 90.12. 1   put if(ialm... out of do-loop (for optimization)
! followings are definitions of symbols in the routine, with the format:
! name    type       content
! where type = s (subroutine), f (function), i (integer), r (real).
! for the type = i and r, the size of array is shown by ( ).
!--- input
! am1      r         cos(emergent zenith angle) with the sign of
!                    transmitted( .gt.0),  reflected( .lt.0).
! am0      r         cos(incident zenith angle) .gt. 0.
! nl       i         no. of sublayers.
! tc     r(nl)       optical thickness of sublayer.
! w      r(nl)       single scattering albedo of sublayer.
! p      r(nl)       phase function of sublayer.
!                    integration over unit sphere should be 1.
! ut       r         user defined optical depth for the intensity
!                    within the multi-layered system.
!--- output
! sglr     f         singly scattered intensity.
!--
  implicit none
  integer,intent(in)::nl
  real(8),intent(in)::am1,am0,tc(nl),w(nl),p(nl),ut
  integer,save::init=1
  real(8),save::eps
  integer::l
  real(8)::ccp(3),x,t1,t2,e1,e2,expfn,sglr
  
  logical ialm,ist
!--
  if(init > 0) then
     init=0
     call cpcon(ccp)
     eps=ccp(1)*30.d0
  endif

! set eps: if abs(1/am1 - 1/am0) .le. eps then
!                      the routine sets almucantar condition-ialm.
  sglr=0.d0; x=1.d0/am1-1.d0/am0
  if(abs(x) <= eps) then
     ialm=.true.
  else
     ialm=.false.
  endif
  t2=0.d0; ist=.true.

! x=0 (transmission)
  if(ialm) then
     do l=1,nl
        t1=t2; t2=t1+tc(l)
        if(ut <= t1) then
           e1=-ut/am1
           sglr=sglr*expfn(e1)/abs(am1)
           return
        endif
        if(ut < t2) t2=ut
        if(ist) then
           ist=.false.; e2=t1
        endif
        e1=e2; e2=t2
        sglr=sglr+w(l)*p(l)*(e2-e1)
     enddo
  else
! x<>0
     do l=1,nl
        t1=t2; t2=t1+tc(l)
!! reflection
        if(am1 < 0) then
           if(ut >= t2) cycle
           if(ut > t1) t1=ut
!! transmission
        else
           if(ut <= t1) then
              e1=-ut/am1
              sglr=sglr*expfn(e1)/abs(am1)
              return
           endif
           if(ut < t2) t2=ut
        endif

        if(ist) then
           ist=.false.
           e2=t1*x
           e2=expfn(e2)/x
        endif
        e1=e2
        e2=t2*x
        e2=expfn(e2)/x
        sglr=sglr+w(l)*p(l)*(e2-e1)
     enddo
  endif
  e1=-ut/am1
  sglr=sglr*expfn(e1)/abs(am1)
  return
end function sglr

subroutine    tnvss2(n,a,dt,e,nn,iw)
! inversion of real matrix. sweep out, complete positioning.
!--- history
! 71. 4.30 created by sakata masato
! 89.11.10 added err
! 90. 1. 6 err*(*)
! 08.10.27 modified fortran 90 version.
!--- input
! n       i        dimension of the matrix
! a     r(nn,n)    matrix
! nn      i        size of first argument
! e       r        convergence criterion (0 is ok)
!--- output
! dt      r        determination of the matrix
! iw    i(2*n)     working area
!--
  use module_opt_visir, only: kndm
  implicit none

! input & output
  integer,intent(in)::n,nn
  real(8),intent(in)::e
  real(8),intent(inout)::a(kndm,kndm)
  integer,intent(out)::iw(2*n)
  real(8),intent(out)::dt

  real(8)::eps,piv,work,aa,az
  integer::ipiv,jpiv,k,kk,ij,i,j,nnn
!--exec      
  if(n <= 0) then
     stop 'MSG tnvss2: error in tinvss: n.le.0'
  elseif(n == 1) then
     dt=a(1,1); k=1
     if(dt == 0.d0) then
        stop 'MSG tinvss: ill conditioned matrix with rank 0'
     endif
     a(1,1)=1.d0/a(1,1)
     return
  else
     if(n > nn) then
        stop 'MSG tnvss2: error in tinvss: n.gt.nn'
     endif
     eps=0.d0; dt=1.d0
     do k=1,n
        piv=0.d0
        do i=k,n
           do j=k,n
              if(abs(a(i,j)) <= abs(piv)) cycle
              ipiv=i; jpiv=j; piv=a(i,j)
           enddo
        enddo
        dt=dt*piv
        if(abs(piv) <= eps) then
           dt=0.d0; nnn=k-1
           print*, 'MSG tinvss: ill conditioned matrix with rank ',nnn
           stop
        endif
        if(k==1)   eps=abs(piv)*e
        if(ipiv/=k) then
           dt=-dt
           do j=1,n
              work=a(ipiv,j)
              a(ipiv,j)=a(k,j)
              a(k,j)=work
           enddo
        endif
        if(jpiv/=k) then
           dt=-dt
           do i=1,n
              work=a(i,jpiv)
              a(i,jpiv)=a(i,k)
              a(i,k)=work
           enddo
        endif
        iw(2*k-1)=ipiv; aa=1.d0/piv
        iw(2*k)=jpiv
        a(k,1:n)=a(k,1:n)*aa
        do i=1,n
           if(i==k) cycle
           az=a(i,k)
           if(az==0.d0) cycle
           a(i,1:n)=a(i,1:n)-a(k,1:n)*az
           a(i,k)=-aa*az
        enddo
        a(k,k)=aa
     enddo
     do kk=2,n
        k=n+1-kk; ij=iw(2*k)
        if(ij/=k) then
           do j=1,n
              work=a(ij,j)
              a(ij,j)=a(k,j)
              a(k,j)=work
           enddo
        endif
           ij=iw(2*k-1)
        if(ij==k) cycle
        do i=1,n
           work=a(i,ij)
           a(i,ij)=a(i,k)
           a(i,k)=work
        enddo
     enddo
     return
  endif
end subroutine tnvss2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine trn1(nsb,ndd,na0,iup,idn,re,te,ser,set,rup,rdn)
! solve the radiative transfer in the multi-sublayer system
!  by the adding method.
!--- history
! 87. 3. 4 for intensity calculation.
! 89.10.30 re-edit.
! 96.01.19   eq22(s2u,st,n2,na0,kndm,kndm,kndm,kndm)
!          ->eq22(s2u,st,n2,na0,kndm,kna0,kndm,kna0)
! 08.10.27 modified fortran 90 version.
!--- input
! nsb       i              number of sublayers.
! ndd    r(nsb1)           number of quadrature points at interfaces.
! iup    i(nsb)            element number of upwelling operators.
! idn    i(nsb)            element number of downgoing operators.
! na0       i              number of solar zenith angles.
! re    r(kndm,kndm,knlt)  reflection matrices of sublayers.
! te    r(kndm,kndm,knlt)  transmission matrices of sublayers.
! ser   r(kndm,kna0,knsb)  upwelling source matrices.
! set   r(kndm,kna0,knsb)  downgoing source matrices.
!--- output
! rup   r(kndm,kna0,knsb1) upwelling internal intensities.
! rdn   r(kndm,kna0,knsb1) downgoing internal intensities.
!--- parameter
! kna0      i              number of solar zenith angles.
! kndm      i              number of quadrature points.
! knsb      i              number of sublayers.
! knlt      i              total number of elementary operators
!                           taking porality into account .ge. knsb.
!--- areas for this routine
  use module_opt_visir,only:kna0,kndm,knlnm,knlnm1
  implicit none
  integer,intent(in)::nsb
  integer,intent(in)::ndd(knlnm1)
  integer,intent(in)::iup(knlnm)
  integer,intent(in)::idn(knlnm)
  integer,intent(in)::na0

  real(8),intent(in)::re(kndm,kndm,knlnm)
  real(8),intent(in)::te(kndm,kndm,knlnm)
  real(8),intent(in)::ser(kndm,kna0,knlnm)
  real(8),intent(in)::set(kndm,kna0,knlnm)

  real(8),intent(out)::rup(kndm,kna0,knlnm1)
  real(8),intent(out)::rdn(kndm,kna0,knlnm1)

!--- working areas
  real(8)::r1d(kndm,kndm),r1u(kndm,kndm),t1d(kndm,kndm),t1u(kndm,kndm), &
       s1d(kndm,kna0),s1u(kndm,kna0),r2d(kndm,kndm), &
       t2d(kndm,kndm),t2u(kndm,kndm),s2d(kndm,kna0),s2u(kndm,kna0),&
       tu(kndm,kndm),rd(kndm,kndm),sd(kndm,kna0),su(kndm,kna0), &
       rl(kndm,kndm,knlnm),sl(kndm,kna0,knlnm),rt(kndm,kndm),st(kndm,kna0)
  integer::nsb1,id,n1,n2,n3,l1,iu,l,i,n
!--exec
!--- upward adding
  nsb1=nsb+1; id=idn(nsb)
  n2=ndd(nsb); n3=ndd(nsb1)
  rl(1:n2,1:n2,nsb)=re(1:n2,1:n2,id)
  sl(1:n2,1:na0,nsb)=ser(1:n2,1:na0,nsb)
  if(nsb >= 2) then
     do l=nsb-1,1,-1
        l1=l+1; n1=ndd(l); iu=iup(l); id=idn(l); n2=ndd(l1)
        r1d(1:n1,1:n1) = re(1:n1,1:n1,id)
        r1u(1:n2,1:n2) = re(1:n2,1:n2,iu)
        t1d(1:n2,1:n1) = te(1:n2,1:n1,id)
        t1u(1:n1,1:n2) = te(1:n1,1:n2,iu)
        s1d(1:n2,1:na0)=set(1:n2,1:na0,l)
        s1u(1:n1,1:na0)=ser(1:n1,1:na0,l)
        r2d(1:n2,1:n2) = rl(1:n2,1:n2,l1)
        s2u(1:n2,1:na0)= sl(1:n2,1:na0,l1)
        call add(0,na0,n1,n2,n3,r1d,r1u,t1d,t1u,s1d,s1u,r2d,t2d,t2u, &
             s2d,s2u,tu,rd,sd,su)
        !if(err /= '') then
        !   err='error in add for upward adding (trn1)'; return
        !endif
        rl(1:n1,1:n1,l) =rd(1:n1,1:n1)
        sl(1:n1,1:na0,l)=su(1:n1,1:na0)
     enddo
  endif

! field at the top of the system
  n1=ndd(1)
  rup(1:n1,1:na0,1)=sl(1:n1,1:na0,1)
  rdn(1:n1,1:na0,1)=0.d0
  
!--- downward adding
  iu=iup(1); n1=ndd(1); n2=ndd(2)
  rt(1:n2,1:n2) = re(1:n2,1:n2,iu)
  st(1:n2,1:na0)=set(1:n2,1:na0,1)
  if(nsb >= 2) then
     do l=2,nsb
        n2=ndd(l)
! internal field
        s2u(1:n2,1:na0)=sl(1:n2,1:na0,l)
        r2d(1:n2,1:n2) =rl(1:n2,1:n2,l)
        sd(1:n2,1:na0)=matmul(rt(1:n2,1:n2),s2u(1:n2,1:na0))
        sd(1:n2,1:na0)=sd(1:n2,1:na0)+st(1:n2,1:na0)
        rd(1:n2,1:n2)=matmul(rt(1:n2,1:n2),r2d(1:n2,1:n2))
        call multi(n2,rd,r1d)
!           if(l==2) write(*,*) n2,rd(1,1),r1d(1,1)
        !if(err /= '') then
        !   err='error in multi for internal field (trn1)'; return
        !endif
        s1d(1:n2,1:na0)=matmul(r1d(1:n2,1:n2),sd(1:n2,1:na0))
        s1u(1:n2,1:na0)=matmul(r2d(1:n2,1:n2),s1d(1:n2,1:na0))
        s1u(1:n2,1:na0)=s1u(1:n2,1:na0)+s2u(1:n2,1:na0)
        rdn(1:n2,1:na0,l)=s1d(1:n2,1:na0)
        rup(1:n2,1:na0,l)=s1u(1:n2,1:na0)
        
! upside down direction for application of the routine add.
        iu=iup(l); id=idn(l); n3=ndd(l+1)
        r1d(1:n3,1:n3)=re(1:n3,1:n3,iu)
        r1u(1:n2,1:n2)=re(1:n2,1:n2,id)
        t1d(1:n2,1:n3)=te(1:n2,1:n3,iu)
        t1u(1:n3,1:n2)=te(1:n3,1:n2,id)
        r2d(1:n2,1:n2)=rt(1:n2,1:n2)
        s1d(1:n2,1:na0)=ser(1:n2,1:na0,l)
        s1u(1:n3,1:na0)=set(1:n3,1:na0,l)
        s2u(1:n2,1:na0)=st(1:n2,1:na0)
        call add(0,na0,n3,n2,n1,r1d,r1u,t1d,t1u,s1d,s1u,r2d,t2d,t2u, &
             s2d,s2u,tu,rt,sd,st)
        !if(err /= '') then
        !   err='error in add for downward adding (trn1)'; return
        !endif
     enddo
  endif

! field at the bottom of the system.
  rdn(1:n3,1:na0,nsb1)=st(1:n3,1:na0)
  rup(1:n3,1:na0,nsb1)=0.d0
!--- downward adding end
  return
end subroutine trn1

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
