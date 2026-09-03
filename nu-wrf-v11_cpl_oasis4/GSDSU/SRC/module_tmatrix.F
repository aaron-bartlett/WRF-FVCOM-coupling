!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 module module_tmatrix

!--------------------------------------------------------------------------------------------
!
!
! Method (Calling sequence) : original
!
! driver_tmatrix ( driver of Tmatrix )
!     !-- init_tmatrix ( read input parameters )
!     |     |-- ray ( get dielectric constant of ice or water  )
!     |-- quad ( Numerical integration using the Gauss-Legendre method. )
!     |     |-- gener (   )
!     |          |-- genbkr ( )
!     |                |-- dvdd ()
!     |-- prcssm ( Calculate T(1,1)= B(inverse)*A and store in B )
!     |-- addprc ( obtain the scattered field coefficients and the differential scattering cross section  )
!           |-- genlgp ( generate Legendre polynomials )
!
!-----------------------------------------------------------------------------------------------

!
! Encapsulation control 
!
  private   ! encapsulate all module parameters and subourtines.

!
! public (accessible) subroutine
!
  public :: driver_tmatrix ! drive tmatrix  


 logical,parameter :: debug_tmatrix = .false.

!
! default value used for CSU HID algorithm database (don't change below....)
!
 integer,parameter :: input_nuang= 2  ! Number of scattering angles (1-180): default=2
 integer,parameter,public :: input_nrank= 17 ! N-value (matrix order): 17
 integer,parameter :: input_nsect= 2  ! Number of sections in the body
 integer,parameter :: input_ib= 8     ! Symmetry code : 8-minor symetry about theta=90deg, 9-for general shaped body
 integer,parameter,public :: input_nm= 7     ! Number of m values 
 real(8) :: input_anginc      = 90.   ! Angle of incidence of incidence wave: 90

 contains

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine driver_tmatrix( axis_ratio, dcreal, dcimag, diam_sph, dwave , &       !input
            tmat_out, nrank_out, nm_out, conk_out, aovrb_out, sigma_out )  !output

!     Electromagnetic scattering from a lossy dielectric imbedded within
!     another lossy dielectric.
!     The exp(-jwt) time convention is used.
!     Perform the numerical integration and fill the A,B,C,D matrices
!     for the outer surface and X,Y matrices for the inner surface.

!     Unit 8 outputs the size & RCS
!     Unit 7 outputs T-matrix
!     Unit 6 outputs input data, RCS and intermediate parameters

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

 implicit double precision (a-h)
 implicit double precision (o-z)
 implicit integer (i-n)

 complex :: tmat_out(20,40,40)   ! output tmatrix B 
!
! keep old common block
!
 common dtr,rtd,cpi
 common /bdycom/ dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib
 common /cmvcom/ nm,kmv,cmi(20),cmv,cm2,twm,prodm
 common /variam/qem,b89
 common /mtxcom/nrank,nranki,a(80,80,2),b(80,80,2),cmxnrm(80)
 common /uvccom/anginc,acans(181,2,2),uang(181),rtsfct,dltang,nuang
 common /rowcol/irow,irow1,crow,crowm,crow1,icol,icol1,ccol,ccolm, &
        ccol1,crij,crssij
 common /brin/ml,ijk,height
 common /endpnt/epps(4),nsect
 common /scatt/summ1,summ2

!
! local parameter
!
 double precision temp(16,8,4)
 double precision clrmtx(25600),clrtot(724),result(16,8)
 double precision buff(5)
 equivalence (a(1,1,1),clrmtx(1)),(acans(1,1,1),clrtot(1))



!
! old CSU tmatrix text output 
!
 if(debug_tmatrix) then
    open(unit=6,file='tmatrix_output_6')
    open(unit=7,file='tmatrix_output_7')
    open(unit=8,file='tmatrix_output_8')
    open(unit=50,file='tmatrix_output_mp.dat')
    open(unit=10,file='tmatrix_output_scatamp.dat')
    write(10,2001)
 endif


!
! initialization constant and other paramters
!
 anginc=0.
 ml = 1  ! number of sizes computed using T-matrix
 cpi = 4.0*atan(1.0)
 dtr = cpi/180.0
 rtd = 180.0/cpi

 if(debug_tmatrix) write(7,*)ml

!
! initialize tmatrix run
!
  call init_tmatrix(axis_ratio, dcreal, dcimag, diam_sph, dwave)

  if(debug_tmatrix) iflag=0  

!
! Set the number of scattering angles (nuang= 181 maximum and
! nuang-1 should be divisible into 180 by a whole number)
!
  nuang= input_nuang
  dltang= float(180/(nuang-1))
  uang(1)= anginc

  do i= 2,nuang
     uang(i)= uang(i-1)+dltang
  enddo 

!
!     Clear the accumulating answer register (used in addprc)
!
 do j=1,724
    clrtot(j)= 0.0
 enddo

 summ1= 0.0
 summ2= 0.0
 rtsfct= 8.0/conk

!
! Set multiplier b89 dependent on ib value (symmetry indicator)
!
 b89 = 1.0
 if (ib.eq.8) b89= 2.0
 bdyfct=1.0

!
! Set up a loop for each m value
!
 m_loop: do im= 1,nm

   !
   ! Set m dependent variables
   !
    cmv= cmi(im)
    kmv= cmv
    cm2= cmv**2
    prodm= 1.0

    if (kmv.gt.0) then
      em= 2.0
      quanm = cmv
      do ifct = 1,kmv
         quanm= quanm+1.0
         prodm= quanm*prodm/2.0
      enddo
      qem= -2.0/em
      twm= 2.0*cmv
    else
      em= 1.0
      qem= -2.0/em
      twm= 2.0*cmv
    endif
    
    !
    ! Initialize all matrix areas to zero
    !
    do i=1,25600
       clrmtx(i)= 0.0
    enddo

    !
    ! Set up a loop for each row of the matrices
    !
    crow=0.0
    crowm = cmv
    do 600 irow= 1,nrank
       irow1 = irow+nrank
       crow=crow+1.0
       crowm=crowm+1.0
       crow1=crow+1.0

       !
       ! Set up a loop for each column of the matrices
       !
       ccol=0.0
       ccolm = cmv

       loop_icol: do icol = 1,nrank
          icol1 = icol+nrank
          ccol=ccol+1.0
          ccolm=ccolm+1.0
          ccol1=ccol+1.0

          ! Calculate matrices A,B associated with the outer surface.
          ! Following notation by Peterson and Strom, q1(out,re) is stored
          ! in A, q1(re,re) is stored in B.
          ! All matrices are transposed in the following code.
          ! Perform integration using a sequence of 31
          ! point extended Gauss-type quadrature formulae.
          ! result(16,k) contains the values of the integrals.
          ! There are 16 integrations to
          ! be performed for each looping through irow and icol.
          ! These correspond to 4 sub-matrix elements for each of the
          ! two matrices A,B and assiciated real and imaginary parts.

          nsect1=nsect-1
          do j=1,nsect1
             js=j
             ith=0
             call quad(epps(j),epps(j+1),k,result,npts,ith,js)
             do i=1,16
                temp(i,k,j)=result(i,k)
             enddo
          enddo

          do j=1,nsect1
              a(icol,irow1,1)= temp(1,k,j) + a(icol,irow1,1)
              a(icol,irow1,2)= temp(2,k,j) + a(icol,irow1,2)
              b(icol,irow1,1)= temp(3,k,j) + b(icol,irow1,1)
              b(icol,irow1,2)= temp(4,k,j) + b(icol,irow1,2)
              a(icol1,irow,1)= temp(5,k,j) + a(icol1,irow,1)
              a(icol1,irow,2)= temp(6,k,j) + a(icol1,irow,2)
              b(icol1,irow,1)= temp(7,k,j) + b(icol1,irow,1)
              b(icol1,irow,2)= temp(8,k,j) + b(icol1,irow,2)
              a(icol1,irow1,1)= temp(9,k,j) + a(icol1,irow1,1)
              a(icol1,irow1,2)= temp(10,k,j) + a(icol1,irow1,2)
              b(icol1,irow1,1)= temp(11,k,j) + b(icol1,irow1,1)
              b(icol1,irow1,2)= temp(12,k,j) + b(icol1,irow1,2)
              a(icol,irow,1)= temp(13,k,j) + a(icol,irow,1)
              a(icol,irow,2)= temp(14,k,j) + a(icol,irow,2)
              b(icol,irow,1)= temp(15,k,j) + b(icol,irow,1)
              b(icol,irow,2)= temp(16,k,j) + b(icol,irow,2)
          enddo

          if(debug_tmatrix) then
            if (iflag.eq.0) write(6,101) npts
            iflag=iflag+1
          endif
       enddo loop_icol

       !
       ! Calculate the normalization factor (used in addprc)
       !
       ckrow = irow

       if (kmv.gt.0) go to 426

       fctki=1.0
       go to 440
  426     if (irow.ge.kmv) go to 430
          cmxnrm(irow)=1.0
          go to 600
  430     ibfct= irow-kmv+1
          iefct= irow+kmv
          fprod= ibfct
          fctki= 1.0
          do lfct= ibfct,iefct
            fctki= fctki*fprod
            fprod=fprod+1.0
          enddo
  440     cmxnrm(irow)= 4.0*ckrow*(ckrow+1.0)*fctki/(em*(2.0*ckrow+1.0))
  600   continue


    nrankx2 = 2*nrank !nrankx2

    !
    ! Process computed matrices
    ! Calculate T(1,1)= B(inverse)*A and store in B
    !
    call prcssm(a,b,nrank,nranki)

    !
    ! tmatrix output for mueller matrix
    !  b(icom2,icom3,1)  is real component
    !  b(icom2,icom3,2)  is imagenary component
    !
    do icom2= 1,nrankx2 ; do icom3=1,nrankx2
       tmat_out(im,icom2,icom3) = cmplx( b(icom2,icom3,1) ,b(icom2,icom3,2) )  !output
    enddo ; enddo  


    if(debug_tmatrix) then
       write(7,11) cmi(im)   !cmi is double float of im
       do icom1= 1,2
          do icom2= 1,nrankx2  
             if(debug_tmatrix) write(7,9) (b(icom2,icom3,icom1),icom3=1,nrankx2)
          enddo
       enddo
    endif

    !
    ! Subroutine to obtain the scattered field coefficients and calculate
    ! the differential scattering cross section in the azimuthal plane.
    !
    if(debug_tmatrix) call addprc

 enddo m_loop

! anginc=anginc+10.
! close(20)



 if(debug_tmatrix) then
    write(6,*) 'end of input'
    close(6)
    close(7)
    close(8)
    close(50)
 endif

!
! output for mueller matrix
!
 nrank_out = nrank
 nm_out    = nm
 conk_out  = conk
 aovrb_out = aovrb 
 sigma_out = sigma


!
! output format statement
!
 9     format(2x,8d15.7)
 11    format(2x,f10.5)
 101   format('number of gauss points used= ',i5,/)
 2001  format('#ang',2x,'de',7x,'re(forw)',5x,'im(forw)', 5x,'re(back)',5x,'im(back)')

 return
 end subroutine driver_tmatrix

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine init_tmatrix (axis_ratio, dcreal, dcimag, diam_sph, dwave)
!  subroutine to read input data for the tmatrix program

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

 implicit double precision (a-h)
 implicit double precision (o-z)
 implicit integer (i-n)

 common dtr,rtd,cpi
 common /gauss/kgauss
 common /mtxcom/nrank,nranki,a(80,80,2),b(80,80,2),cmxnrm(80)
 common /cmvcom/ nm,kmv,cmi(20),cmv,cm2,twm,prodm
 common /bdycom/ dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib
 common /uvccom/ anginc,acans(181,2,2),uang(181),rtsfct,dtlang,nuang
 common /brin/ml,ijk,height
 common /endpnt/epps(4),nsect
 common /vivek2/dpart

 double precision epdeg(4),epr,epi
 integer :: ph
 complex rhoia


!     Read necessary input data.
!     Card1:  nm = number of m values,nrank = n value(matrix order),
!     nsect = number of sections in the body,ib = symmetry code  ib = 8
!     for mirror symmetry about theta= 90 degrees,ib = 9 for general
!     shaped body, anginc = angle of incidence of the incident wave.

 nrank  = input_nrank
 nsect  = input_nsect
 ib     = input_ib
 nm     = input_nm
 anginc = input_anginc

 nranki = nrank+1

!     Read in parameters for the outer surface.
!     Card 2:  conk= ka of body, aovrb = a/b ratio of spheroids or any
!     needed constant, sigma = 4/k (used with scattered field addprc) or
!     any needed constant, dcnr = real part of relative dielectric constant,
!     dcni= imaginary part of relative dielectric constant.
!     dcnr and dcni are both positive since exp(-jwt) time convention is
!     assumed.

 dcnr=dcreal
 dcni=dcimag

!
! diameter of hydrometeor
!
 dpart = diam_sph  !spherical diameter [cm]

!
! ambient temperature (not used any more) (C) 
!     T=20C for MCTEX rain
!     T=15C to compare to Bringi 90
!     T=0C for wet-growth hail
!     Temperature (tempt, deg C)
 tempt=20.0

!
! Wavelength (alamd, cm)
!     CSU-CHILL, S-band: 11.0 cm
!     C-pol, C-band: 5.5 - 5.625 GHZ (5.3 - 5.5 cm)
 alamd=dwave  !wavelength in 1cm

 sigma=2.*alamd/(cpi*100.)


! RAIN a/b (axis_ratio)
!      Investigate new D vs a/b relationship developed by
!      Keenan et al. for use w/ MCTEX data where they fit
!      a cubic equation to the Pruppacher and Pitter (1971),
!      Green (1975), and Jameson (1983) where Jameson's
!      results were significantly less oblate (not equilibrium
!      drop shape but measured by camera in Florida storms).
!      a/b=0.9934+0.0822D-1.8726D^2+1.4685D^3
!      LDC: 1/15/97
!     *** CHANGE THIS:  PICK one aovrb relationship! ***
!      aovrb=.9934+.0822*dpart-1.8726*dpart**2+
!     +      1.4685*dpart**3

 aovrb = axis_ratio  ! axis ratio is now from subroutine input

 conk= cpi*dpart/(alamd*(aovrb**(1./3.)))

!
! text output
!
 if(debug_tmatrix) then
      write(6,"(i15,  5x,'symmetry code')")ib
      write(6,"(i15,  5x,'number of sections')")nsect
      write(6,"(i15,  5x,'number of sizes')")ml
      write(6,"(i15,  5x,'nrank')")nrank
      write(6,"(i15,  5x,'nmode')")nm
      write(6,"(d15.7,5x,'conkb')")conk
      write(6,"(d15.7,5x,'wavelength in cm')")alamd
      write(6,"(d15.7,5x,'temperature in c')")tempt
      write(6,"(d15.7,5x,'equivalent particle diameter in cm')")dpart
      write(6,"(d15.7,5x,'a over b ratio')")aovrb
      write(6,"(d15.7,5x,'sigma')")sigma
      write(6,"(d15.7,5x,'real part of dielectric const')")dcnr
      write(6,"(d15.7,5x,'imag part of dielectric const')")dcni
 endif

!     Unit 7 eventually outputs the T-matrix for each size
!     ldc 10/17/95 moved this write to main program
!      write(7,*)ml
  if(debug_tmatrix) then
      write(7,"(i15,  5x,'nrank')")nrank
      write(7,"(i15,  5x,'nmode')")nm
      write(7,"(d15.7,5x,'conkb')")conk
      write(7,"(d15.7,5x,'wavelength in cm')")alamd
      write(7,"(d15.7,5x,'equivalent particle diameter in cm')")dpart
      write(7,"(d15.7,5x,'aoverb')")aovrb
      write(7,"(d15.7,5x,'sigma')")sigma
      write(7,"(d15.7,5x,'real part of dielectric const')")dcnr
      write(7,"(d15.7,5x,'imag part of dielectric const')")dcni
  endif

  if(debug_tmatrix) then
!     Unit 8 eventually outputs the field amplitudes
      write(8,"(d15.7,5x,'conk')")conk
      write(8,"(d15.7,5x,'equivalent particle diameter in cm')")dpart
      write(8,"(d15.7,5x,'aoverb')")aovrb
      write(8,"(d15.7,5x,'sigma')")sigma
      write(8,"(d15.7,5x,'real part of dielectric const')")dcnr
      write(8,"(d15.7,5x,'imag part of dielectric const')")dcni
  endif

!     Card 4:  cmi = m values (only m=1 is required for anginc=0,but
!     m=0,1,2, ... is required for general incidence angle).

  if (nm.eq.1) then
     cmi(1)= 1.0
  else
     do 7 i= 1,nm
 7   cmi(i)= float(i-1)
  endif

  kgauss=5
  call calenp
  do i=1,nsect
     epdeg(i)=rtd*epps(i)
  enddo

  if(debug_tmatrix) write(6,148) (epdeg(i),i=1,nsect)
148  format('end points ',8d12.4,/(1x,23x,8d12.4))

 return
 end subroutine init_tmatrix

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine gener(xxx,t,ith,js)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      common dtr,rtd,cpi
      common /thtcom/theta,sinth,costh
      common /mtxcom/nrank,nranki,a(80,80,2),b(80,80,2),cmxnrm(80)
      common /bdycom/ dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib
      common /fnccom/pnmllg(81),bsslsp(81,31,3),cneumn(81,31,3), &
       bslkpr(81,31,3),bslkpi(81,31,3),cneumr(81,31,3),cneumi(81,31,3)
      common /rowcol/irow,irow1,crow,crowm,crow1,icol,icol1,ccol,ccolm, &
       ccol1,crij,crssij
      common /cmvcom/ nm,kmv,cmi(20),cmv,cm2,twm,prodm
      common /variam/qem,b89

      double precision rhankl(80,2),rbskpr(80,2),rbessl(80),hanklr(81), &
        hankli(81),rhskpr(80,2)
      double precision t(4,2,2)
      xxxr=xxx*rtd
      do 1101 i=1,4
        do 1101 j=1,2
          do 1101 k=1,2
            t(i,j,k)= 0.0
1101  continue
      sqreal= crootr(dcnr,dcni)
      sqimag= crooti(dcnr,dcni)
      qsreal= cddvdr(1.0d0  ,0.0d0  ,sqreal,sqimag)
      qsimag= cddvdi(1.0d0  ,0.0d0  ,sqreal,sqimag)
      theta= xxx
      costh= cos(theta)
      sinth= sin(theta)
      srmsin= sinth
      tempp= cpi-theta
      if (abs(tempp).lt.1.0e-8) sinth=0.0

!     Generate the Legendre polynomials.
      call genlgp
!     Evaluate kr and its derivative as a function of theta.
  348 call genkr
!     Generate all necessary Bessel and Neumann functions and their ratios.
      if((irow.eq.1).and.(icol.eq.1)) call genbsl(ith,js)

!     Calculate k*r for outer surface.
      ckprr = sqreal*ckr
      ckpri = sqimag*ckr
      iswt=1
      if ((irow.eq.1).and.(icol.eq.1)) call genbkr(ckprr,ckpri,iswt, &
       ith,js)
      do 347 k=1,nranki
!v      write(6,'(a3,i4,4d15.7,/)') 'new',k,bslkpr(k,ith,js),
!v     1bslkpi(k,ith,js),cneumr(k,ith,js),cneumi(k,ith,js)
      hanklr(k)=bslkpr(k,ith,js)-cneumi(k,ith,js)
      hankli(k)=bslkpi(k,ith,js)+cneumr(k,ith,js)
347    continue
      do 350 k = 1,nrank
      rbessl(k) = bsslsp(k,ith,js)/bsslsp(k+1,ith,js)
      rhankl(k,1)=cddvdr(bsslsp(k,ith,js),cneumn(k,ith,js), &
      bsslsp(k+1,ith,js),cneumn(k+1,ith,js))
      rhankl(k,2) = cddvdi(bsslsp(k,ith,js),cneumn(k,ith,js), &
      bsslsp(k+1,ith,js),cneumn(k+1,ith,js))
      rbskpr(k,1) = cddvdr(bslkpr(k,ith,js),bslkpi(k,ith,js), &
      bslkpr(k+1,ith,js),bslkpi(k+1,ith,js))
      rbskpr(k,2) = cddvdi(bslkpr(k,ith,js),bslkpi(k,ith,js), &
      bslkpr(k+1,ith,js),bslkpi(k+1,ith,js))
      bkr = cdmpyr(sqreal,sqimag,rbskpr(k,1),rbskpr(k,2))
      bki = cdmpyi(sqreal,sqimag,rbskpr(k,1),rbskpr(k,2))
      rbskpr(k,1) = bkr
      rbskpr(k,2) = bki
      tempnr=bslkpr(k,ith,js)-cneumi(k,ith,js)
      tempni=bslkpi(k,ith,js)+cneumr(k,ith,js)
      tempdr=bslkpr(k+1,ith,js)-cneumi(k+1,ith,js)
      tempdi=bslkpi(k+1,ith,js)+cneumr(k+1,ith,js)
      rhskpr(k,1)=cddvdr(tempnr,tempni,tempdr,tempdi)
      rhskpr(k,2)=cddvdi(tempnr,tempni,tempdr,tempdi)
      hkr=cdmpyr(sqreal,sqimag,rhskpr(k,1),rhskpr(k,2))
      hki=cdmpyi(sqreal,sqimag,rhskpr(k,1),rhskpr(k,2))
      rhskpr(k,1)=hkr
      rhskpr(k,2)=hki
  350 continue
      br = rbessl(irow)
      hr = rhankl(irow,1)
      hi = rhankl(irow,2)

!     Calculate frequently used variable combinations for use in a,b
!     matrices.
      crij = crow+ccol
      crssij = crow*ccol
      cmcrco = cm2-qem*crssij*costh**2
      pnr0c0 = pnmllg(irow)*pnmllg(icol)
      pnr0c1 = pnmllg(irow)*pnmllg(icol+1)
      pnr1c0 = pnmllg(irow+1)*pnmllg(icol)
      pnr1c1 = pnmllg(irow+1)*pnmllg(icol+1)
      b1a = crow*costh*pnr1c1-crowm*pnr0c1
      b1b = ccol*costh*pnr1c1-ccolm*pnr1c0
      bkr = rbskpr(icol,1)
      bki = rbskpr(icol,2)
      hkr=rhskpr(icol,1)
      hki=rhskpr(icol,2)
      hbkmlr=cdmpyr(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js), &
      bslkpr(icol+1,ith,js),bslkpi(icol+1,ith,js))
      hbkmli=cdmpyi(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js), &
      bslkpr(icol+1,ith,js),bslkpi(icol+1,ith,js))
      bbkmlr = bsslsp(irow+1,ith,js)*bslkpr(icol+1,ith,js)
      bbkmli = bsslsp(irow+1,ith,js)*bslkpi(icol+1,ith,js)
      hhkmlr=cdmpyr(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js), &
      hanklr(icol+1),hankli(icol+1))
      hhkmli=cdmpyi(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js), &
      hanklr(icol+1),hankli(icol+1))
      bhkmlr=bsslsp(irow+1,ith,js)*hanklr(icol+1)
      bhkmli=bsslsp(irow+1,ith,js)*hankli(icol+1)
      hepsr = cdmpyr(qsreal,qsimag,hbkmlr,hbkmli)
      hepsi = cdmpyi(qsreal,qsimag,hbkmlr,hbkmli)
      bepsr = cdmpyr(qsreal,qsimag,bbkmlr,bbkmli)
      bepsi = cdmpyi(qsreal,qsimag,bbkmlr,bbkmli)
      hhepsr=cdmpyr(qsreal,qsimag,hhkmlr,hhkmli)
      hhepsi=cdmpyi(qsreal,qsimag,hhkmlr,hhkmli)
      bbepsr=cdmpyr(qsreal,qsimag,bhkmlr,bhkmli)
      bbepsi=cdmpyi(qsreal,qsimag,bhkmlr,bhkmli)
       ftheta=theta*rtd
!      write(6,902) js,ith,irow,icol,ftheta,bkpr2(irow+1,ith,js),
!     1 bkpi2(irow+1,ith,js),bsslsp(irow+1,ith,js)
      if(ib.eq.9) go to 380
!     if ib = 8 (mirror symmetry body), i=l=0 if irow and icol are both
!     odd or both even, j=k=0 if irow and icol are odd,even or even,odd.
      if((irow+icol).eq.((irow+icol)/2)*2) go to 392
!     test for m=0 (if m=0 the i and l submatrices are zero).
  380 if(kmv.eq.0) go to 390
!     fill out elements for equivalent i-submatrix position.
      b1 = b1a+b1b
      htbkr = cdmpyr(hr,hi,bkr,bki)
      htbki = cdmpyi(hr,hi,bkr,bki)
      btbkr = br*bkr
      btbki = br*bki
      tempp=(crow*crow1*bkr+ccol*ccol1*hr-crssij*(crij+2.0  )/ckr &
      )*dckr*sinth
      sumar=tempp*pnr1c1
      sumr=(ckr*(1.0  +htbkr)-ccol*hr-crow*bkr+crssij/ckr)*b1*ckr+sumar
      sumai = pnr1c1*(crow*crow1*bki+ccol*ccol1*hi)*dckr*sinth
      sumi = (ckr*htbki-ccol*hi-crow*bki)*b1*ckr+sumai
      t(1,1,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,hbkmlr,hbkmli)
      t(1,1,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,hbkmlr,hbkmli)
      sumbr = pnr1c1*(crow*crow1*bkr+ccol*ccol1*br-crssij*(crij+2.0d0)/ckr)*dckr*sinth
      sumr=(ckr*(1.0  +btbkr)-ccol*br-crow*bkr+crssij/ckr)*b1*ckr+sumbr
      sumbi = pnr1c1*crow*crow1*bki*dckr*sinth
      sumi = (ckr*btbki-crow*bki)*b1*ckr+sumbi
      t(1,2,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,bbkmlr,bbkmli)
      t(1,2,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,bbkmlr,bbkmli)
      bthkr=br*hkr
      bthki=br*hki

!     Fill out elements for equivalent l-submatrix position.
      sumr=(ckr*(dcnr+htbkr)-ccol*hr-crow*bkr+crssij/ckr)*b1*ckr+sumar
      sumi=(ckr*(dcni+htbki)-ccol*hi-crow*bki)*b1*ckr+sumai
      t(2,1,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,hepsr,hepsi)
      t(2,1,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,hepsr,hepsi)
      sumr = (ckr*(dcnr+btbkr)-ccol*br-crow*bkr+crssij/ckr)*b1*ckr+sumbr
      sumi = (ckr*(dcni+btbki)-crow*bki)*b1*ckr+sumbi
      t(2,2,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,bepsr,bepsi)
      t(2,2,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,bepsr,bepsi)
390   if (ib.eq.8) go to 400

!     Fill out elements for eqiivalent j-submatrix position.
  392 a12=cmcrco*pnr1c1+qem*(crow*ccolm*costh*pnr1c0+ccol*crowm*costh*pnr0c1-crowm*ccolm*pnr0c0)
      b1a = ccol*ccol1*b1a
      b1b = crow*crow1*b1b
      b1 = (b1a-b1b)*sinth
      dd=-qem*dckr
      cr = cdmpyr(dcnr,dcni,hr,hi)
      ci = cdmpyi(dcnr,dcni,hr,hi)
      sumr=(ckr*(bkr-cr)+dcnr*crow-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*dd
      sumi=(ckr*(bki-ci)+dcni*crow)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,1,1)=b89*srmsin*cdmpyr(sumr,sumi,hepsr,hepsi)
      t(3,1,2)=b89*srmsin*cdmpyi(sumr,sumi,hepsr,hepsi)
      cr = br*dcnr
      ci = br*dcni
      sumr=(ckr*(bkr-cr)+dcnr*crow-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*dd
      sumi=(ckr*(bki-ci)+dcni*crow)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,2,1)=b89*srmsin*cdmpyr(sumr,sumi,bepsr,bepsi)
      t(3,2,2)=b89*srmsin*cdmpyi(sumr,sumi,bepsr,bepsi)

!     Fill out elements for equivalent k-submatrix.
      sumr = (ckr*(bkr-hr)+crow-ccol)*a12*ckr+b1*dd
      sumi = (ckr*(bki-hi))*a12*ckr
      t(4,1,1)=b89*srmsin*cdmpyr(sumr,sumi,hbkmlr,hbkmli)
      t(4,1,2)=b89*srmsin*cdmpyi(sumr,sumi,hbkmlr,hbkmli)
      sumr = (ckr*(bkr-br)+crow-ccol)*a12*ckr+b1*dd
      sumi = bki*a12*ckr**2
      t(4,2,1)=b89*srmsin*cdmpyr(sumr,sumi,bbkmlr,bbkmli)
      t(4,2,2)=b89*srmsin*cdmpyi(sumr,sumi,bbkmlr,bbkmli)

400   return

 end subroutine gener

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine quad(a,b,k,result,npts,ith,js)
!     Numerical integration using the Gauss-Legendre method.

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

!     double precision p,d1,d2,d3,d4,d5,d6,d7,d8
      common/gauss/kgauss
      double precision p(381),funct(16,127),fzero(4,2,2),acum(16), &
      test(16),t(4,2,2),result(16,8)
      double precision d1(54),d2(54),d3(54),d4(54),d5(54),d6(54), &
      d7(54),d8(3)
      equivalence (p(1),d1(1)),(p(55),d2(1)),(p(109),d3(1)),(p(163),d4 &
      (1)),(p(217),d5(1)),(p(271),d6(1)),(p(325),d7(1)),(p(379),d8(1))

      data d1/ &
       7.74596669241483d-01, 5.55555555555557d-01, 8.88888888888889d-01, &
       2.68488089868333d-01, 9.60491268708019d-01, 1.04656226026467d-01, &
       4.34243749346802d-01, 4.01397414775962d-01, 4.50916538658474d-01, &
       1.34415255243784d-01, 5.16032829970798d-02, 2.00628529376989d-01, &
       9.93831963212756d-01, 1.70017196299402d-02, 8.88459232872258d-01, &
       9.29271953151245d-02, 6.21102946737228d-01, 1.71511909136392d-01, &
       2.23386686428967d-01, 2.19156858401588d-01, 2.25510499798206d-01, &
       6.72077542959908d-02, 2.58075980961766d-02, 1.00314278611795d-01, &
       8.43456573932111d-03, 4.64628932617579d-02, 8.57559200499902d-02, &
       1.09578421055925d-01, 9.99098124967666d-01, 2.54478079156187d-03, &
       9.81531149553739d-01, 1.64460498543878d-02, 9.29654857429739d-01, &
       3.59571033071293d-02, 8.36725938168868d-01, 5.69795094941234d-02, &
       7.02496206491528d-01, 7.68796204990037d-02, 5.31319743644374d-01, &
       9.36271099812647d-02, 3.31135393257977d-01, 1.05669893580235d-01, &
       1.12488943133187d-01, 1.11956873020953d-01, 1.12755256720769d-01, &
       3.36038771482077d-02, 1.29038001003512d-02, 5.01571393058995d-02, &
       4.21763044155885d-03, 2.32314466399103d-02, 4.28779600250078d-02, &
       5.47892105279628d-02, 1.26515655623007d-03, 8.22300795723591d-03/

      data d2/ &
       1.79785515681282d-02, 2.84897547458336d-02, 3.84398102494556d-02, &
       4.68135549906281d-02, 5.28349467901166d-02, 5.59784365104763d-02, &
       9.99872888120358d-01, 3.63221481845531d-04, 9.97206259372224d-01, &
       2.57904979468569d-03, 9.88684757547428d-01, 6.11550682211726d-03, &
       9.72182874748583d-01, 1.04982469096213d-02, 9.46342858373402d-01, &
       1.54067504665595d-02, 9.10371156957005d-01, 2.05942339159128d-02, &
       8.63907938193691d-01, 2.58696793272147d-02, 8.06940531950218d-01, &
       3.10735511116880d-02, 7.39756044352696d-01, 3.60644327807826d-02, &
       6.62909660024781d-01, 4.07155101169443d-02, 5.77195710052045d-01, &
       4.49145316536321d-02, 4.83618026945841d-01, 4.85643304066732d-02, &
       3.83359324198731d-01, 5.15832539520484d-02, 2.77749822021825d-01, &
       5.39054993352661d-02, 1.68235251552208d-01, 5.54814043565595d-02, &
       5.63443130465928d-02, 5.62776998312542d-02, 5.63776283603847d-02, &
       1.68019385741038d-02, 6.45190005017574d-03, 2.50785696529497d-02, &
       2.10881524572663d-03, 1.16157233199551d-02, 2.14389800125039d-02, &
       2.73946052639814d-02, 6.32607319362634d-04, 4.11150397865470d-03, &
       8.98927578406411d-03, 1.42448773729168d-02, 1.92199051247278d-02, &
       2.34067774953141d-02, 2.64174733950583d-02, 2.79892182552381d-02/

      data d3/ &
       1.80739564445388d-04, 1.28952408261042d-03, 3.05775341017553d-03, &
       5.24912345480885d-03, 7.70337523327974d-03, 1.02971169579564d-02, &
       1.29348396636074d-02, 1.55367755558440d-02, 1.80322163903913d-02, &
       2.03577550584721d-02, 2.24572658268161d-02, 2.42821652033366d-02, &
       2.57916269760242d-02, 2.69527496676331d-02, 2.77407021782797d-02, &
       2.81388499156271d-02, 9.99982430354891d-01, 5.05360952078625d-05, &
       9.99598799671912d-01, 3.77746646326985d-04, 9.98316635318407d-01, &
       9.38369848542380d-04, 9.95724104698407d-01, 1.68114286542147d-03, &
       9.91495721178104d-01, 2.56876494379402d-03, 9.85371499598521d-01, &
       3.57289278351730d-03, 9.77141514639705d-01, 4.67105037211432d-03, &
       9.66637851558417d-01, 5.84344987583563d-03, 9.53730006425761d-01, &
       7.07248999543356d-03, 9.38320397779592d-01, 8.34283875396818d-03, &
       9.20340025470011d-01, 9.64117772970252d-03, 8.99744899776941d-01, &
       1.09557333878379d-02, 8.76513414484705d-01, 1.22758305600827d-02, &
       8.50644494768350d-01, 1.35915710097655d-02, 8.22156254364980d-01, &
       1.48936416648152d-02, 7.91084933799848d-01, 1.61732187295777d-02, &
       7.57483966380512d-01, 1.74219301594641d-02, 7.21423085370098d-01, &
       1.86318482561388d-02, 6.82987431091078d-01, 1.97954950480975d-02/

      data d4/ &
       6.42276642509760d-01, 2.09058514458120d-02, 5.99403930242243d-01, &
       2.19563663053178d-02, 5.54495132631931d-01, 2.29409642293877d-02, &
       5.07687757533716d-01, 2.38540521060385d-02, 4.59130011989833d-01, &
       2.46905247444876d-02, 4.08979821229888d-01, 2.54457699654648d-02, &
       3.57403837831532d-01, 2.61156733767061d-02, 3.04576441556714d-01, &
       2.66966229274503d-02, 2.50678730303482d-01, 2.71855132296248d-02, &
       1.95897502711100d-01, 2.75797495664819d-02, 1.40424233152560d-01, &
       2.78772514766137d-02, 8.44540400837110d-02, 2.80764557938172d-02, &
       2.81846489497457d-02, 2.81763190330167d-02, 2.81888141801924d-02, &
       8.40096928705192d-03, 3.22595002508787d-03, 1.25392848264749d-02, &
       1.05440762286332d-03, 5.80786165997757d-03, 1.07194900062519d-02, &
       1.36973026319907d-02, 3.16303660822264d-04, 2.05575198932735d-03, &
       4.49463789203206d-03, 7.12243868645840d-03, 9.60995256236391d-03, &
       1.17033887476570d-02, 1.32087366975291d-02, 1.39946091276191d-02, &
       9.03727346587510d-05, 6.44762041305726d-04, 1.52887670508776d-03, &
       2.62456172740443d-03, 3.85168761663987d-03, 5.14855847897819d-03, &
       6.46741983180368d-03, 7.76838777792199d-03, 9.01610819519566d-03, &
       1.01788775292361d-02, 1.12286329134080d-02, 1.21410826016683d-02/

      data d5/ &
       1.28958134880121d-02, 1.34763748338165d-02, 1.38703510891399d-02, &
       1.40694249578135d-02, 2.51578703842806d-05, 1.88873264506505d-04, &
       4.69184924247851d-04, 8.40571432710723d-04, 1.28438247189701d-03, &
       1.78644639175865d-03, 2.33552518605716d-03, 2.92172493791781d-03, &
       3.53624499771678d-03, 4.17141937698409d-03, 4.82058886485126d-03, &
       5.47786669391895d-03, 6.13791528004137d-03, 6.79578550488277d-03, &
       7.44682083240758d-03, 8.08660936478883d-03, 8.71096507973207d-03, &
       9.31592412806942d-03, 9.89774752404876d-03, 1.04529257229060d-02, &
       1.09781831526589d-02, 1.14704821146939d-02, 1.19270260530193d-02, &
       1.23452623722438d-02, 1.27228849827324d-02, 1.30578366883530d-02, &
       1.33483114637252d-02, 1.35927566148124d-02, 1.37898747832410d-02, &
       1.39386257383068d-02, 1.40382278969086d-02, 1.40881595165083d-02, &
       9.99997596379750d-01, 6.93793643241083d-06, 9.99943996207055d-01, &
       5.32752936697805d-05, 9.99760490924434d-01, 1.35754910949228d-04, &
       9.99380338025023d-01, 2.49212400482998d-04, 9.98745614468096d-01, &
       3.89745284473282d-04, 9.97805354495956d-01, 5.54295314930373d-04, &
       9.96514145914890d-01, 7.40282804244503d-04, 9.94831502800622d-01, &
       9.45361516858527d-04, 9.92721344282788d-01, 1.16748411742996d-03/

      data d6/ &
       9.90151370400771d-01, 1.40490799565515d-03, 9.87092527954033d-01, &
       1.65611272815445d-03, 9.83518657578632d-01, 1.91971297101387d-03, &
       9.79406281670862d-01, 2.19440692536384d-03, 9.74734459752401d-01, &
       2.47895822665757d-03, 9.69484659502459d-01, 2.77219576459345d-03, &
       9.63640621569812d-01, 3.07301843470258d-03, 9.57188216109859d-01, &
       3.38039799108691d-03, 9.50115297521293d-01, 3.69337791702565d-03, &
       9.42411565191083d-01, 4.01106872407503d-03, 9.34068436157727d-01, &
       4.33264096809299d-03, 9.25078932907077d-01, 4.65731729975685d-03, &
       9.15437587155765d-01, 4.98436456476553d-03, 9.05140358813263d-01, &
       5.31308660518706d-03, 8.94184568335557d-01, 5.64281810138445d-03, &
       8.82568840247341d-01, 5.97291956550816d-03, 8.70293055548114d-01, &
       6.30277344908575d-03, 8.57358310886234d-01, 6.63178124290190d-03, &
       8.43766882672707d-01, 6.95936140939044d-03, 8.29522194637402d-01, &
       7.28494798055382d-03, 8.14628787655138d-01, 7.60798966571904d-03, &
       7.99092290960843d-01, 7.92794933429486d-03, 7.82919394118284d-01, &
       8.24430376303287d-03, 7.66117819303759d-01, 8.55654356130769d-03, &
       7.48696293616938d-01, 8.86417320948252d-03, 7.30664521242183d-01, &
       9.16671116356077d-03, 7.12033155362253d-01, 9.46368999383007d-03/

      data d7/ &
       6.92813769779114d-01, 9.75465653631741d-03, 6.73018830230419d-01, &
       1.00391720440569d-02, 6.52661665410019d-01, 1.03168123309476d-02, &
       6.31756437711193d-01, 1.05871679048852d-02, 6.10318113715188d-01, &
       1.08498440893373d-02, 5.88362434447664d-01, 1.11044611340069d-02, &
       5.65905885423653d-01, 1.13506543159806d-02, 5.42965666498311d-01, &
       1.15880740330440d-02, 5.19559661537457d-01, 1.18163858908302d-02, &
       4.95706407918762d-01, 1.20352707852796d-02, 4.71425065871658d-01, &
       1.22444249816120d-02, 4.46735387662029d-01, 1.24435601907140d-02, &
       4.21657686626164d-01, 1.26324036435421d-02, 3.96212806057616d-01, &
       1.28106981638774d-02, 3.70422087950079d-01, 1.29782022395374d-02, &
       3.44307341599437d-01, 1.31346900919602d-02, 3.17890812068477d-01, &
       1.32799517439305d-02, 2.91195148518247d-01, 1.34137930851101d-02, &
       2.64243372410927d-01, 1.35360359349562d-02, 2.37058845589829d-01, &
       1.36465181025713d-02, 2.09665238243181d-01, 1.37450934430019d-02, &
       1.82086496759252d-01, 1.38316319095064d-02, 1.54346811481378d-01, &
       1.39060196013255d-02, 1.26470584372302d-01, 1.39681588065169d-02, &
       9.84823965981194d-02, 1.40179680394566d-02, 7.04069760428552d-02, &
       1.40553820726499d-02, 4.22691647653637d-02, 1.40803519625536d-02/

      data d8/ &
       1.40938864107825d-02, 1.40928450691604d-02, 1.40944070900962d-02/

      if (a.eq.b) go to 107
      sum=(b+a)/2.0
      diff=(b-a)/2.0
!     One point formula.
!     Set up variable combinations for use in evaluation of integrands
      ith=ith+1
      call gener(sum,t,ith,js)
      do 1010 ii=1,4
        do 1020 jj=1,2
          do 1030 kk=1,2
!           jj=1,2 corresponds to matrices a and b . ii=1,4
!           corresponds to filling out eq. i,l,j,k positions or sub-matrices
!           kk= 1,2 corresponds to real and imaginary parts.
            fzero(ii,jj,kk)= t(ii,jj,kk)
            l= (ii-1)*4+(jj-1)*2+kk
            result(l,1)= 2.0*t(ii,jj,kk)*diff
1030      continue
1020    continue
1010  continue
      i=0
      iold=0
      inew=1
      k=2
      do 1040 n=1,16
        acum(n)=0.0
1040  continue
      go to 103
101   continue
      if (k.eq.kgauss) go to 105
      k=k+1
      do 1050 n=1,16
        acum(n)=0.0
1050  continue

!     Contribution from function values already computed.
      do 102 j=1,iold
        i=i+1
        do 1060 n=1,16
          acum(n)=acum(n)+p(i)*funct(n,j)
1060    continue
102   continue

!     Contribution from new values.
103   continue
      iold=iold+inew
      do 104 j=inew,iold
        i=i+1
        x=p(i)*diff
        temp1=sum+x
        ith=ith+1
        call gener(temp1,t,ith,js)
        do 1070 ii=1,4
          do 1080 jj=1,2
            do 1090 kk=1,2
              l=(ii-1)*4+(jj-1)*2+kk
!             l goes from 1 to 16.
              funct(l,j)=t(ii,jj,kk)
1090        continue
1080      continue
1070    continue
        temp2=sum-x
        ith=ith+1
        call gener(temp2,t,ith,js)
        do 1071 ii=1,4
          do 1081 jj=1,2
            do 1091 kk=1,2
              l=(ii-1)*4+(jj-1)*2+kk
              funct(l,j)=funct(l,j)+t(ii,jj,kk)
1091        continue
1081      continue
1071    continue
        i=i+1
        do 1100 n=1,16
          acum(n)=acum(n)+p(i)*funct(n,j)
1100    continue
104   continue
      inew=iold+1
      i=i+1
      do 1200 ii=1,4
        do 1300 jj=1,2
          do 1400 kk=1,2
            n=(ii-1)*4+(jj-1)*2+kk
            result(n,k)=(acum(n)+p(i)*fzero(ii,jj,kk))*diff
1400      continue
1300    continue
1200  continue
      go to 101
105   continue

!     Normal termination
106   continue
      npts= inew+iold
      return

!     Trivial case
 107  continue
      k=2
      do 1600 m=1,16
        do 1700 n=1,2
          result(m,n)=0.0
1700    continue
1600  continue
      npts=0
      return

 end subroutine quad

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine genlgp

!     subroutine to generate Legendre polynomials.
!     The index on the function is incremented by one.

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      common dtr,rtd,cpi
      common/mtxcom/nrank,nranki,a(80,80,2),b(80,80,2),cmxnrm(80)
      common/fnccom/pnmllg(81),bsslsp(81,31,3),cneumn(81,31,3), &
      bslkpr(81,31,3),bslkpi(81,31,3),cneumr(81,31,3),cneumi(81,31,3)
      common /cmvcom/ nm,kmv,cmi(20),cmv,cm2,twm,prodm
      common /thtcom/ theta,sinth,costh

      dtwm=twm+1.0
!     This is special case when theta equals cpiand m=0.
!     When theta equals cpi all integrands are 0 and any values can be
!     put in pnmllg(41). Here we have put them equal to 0.
      if((sinth.eq.0.0  ).and.(kmv.eq.0)) go to 6
!     At this point theta lies strictly between 0 and cpi.
      if (theta) 16,4,16
    4 if (kmv-1) 6,12,6
    6 do 8 ilg = 1,nranki
      pnmllg(ilg)= 0.0
    8 continue
      go to 88
12    pnmllg(1)= 0.0
      pnmllg(2)= 1.0
      pla= 1.0
      go to 48
   16 if (kmv) 20,20,40

!     The special case when m= 0.
20    pla=1.0/sinth
      plb = costh*pla
      pnmllg(1) = pla
      pnmllg(2) = plb
      ibeg = 3
      go to 60

!     General case for m not equal to 0.
   40 do 44 ilg = 1,kmv
      pnmllg(ilg)=0.0
   44 continue
      if((sinth.eq.0.0  ).and.(kmv.eq.1)) go to 1001
      pla = prodm*sinth**(kmv-1)
      go to 1002
1001  pla=0.0
1002  continue
      pnmllg(kmv+1) = pla
   48 plb = dtwm*costh*pla
      pnmllg(kmv+2) = plb
      ibeg = kmv+3

!     Do recursion formula for all remaining Legendre polynomials.
   60 cnmul = ibeg+ibeg-3
      cnm= 2.0
      cnmm = dtwm
      do 80 ilgr = ibeg,nranki
      plc = (cnmul*costh*plb-cnmm*pla)/cnm
      pnmllg(ilgr) = plc
      pla = plb
      plb = plc
      cnmul= cnmul+2.0
      cnm= cnm+1.0
      cnmm= cnmm+1.0
   80 continue
   88 return

 end subroutine genlgp

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 double precision function cdabx(a,b)

!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      if(a) 4,22,4
    4 if(b) 8,30,8
!   8 e = amax1(a,b)
!     f = amin1(a,b)
    8 e = max(a,b)
      f = min(a,b)
      g = f/e
      cdabx= abs(e)*sqrt(1.0  +g*g)
      return
   22 if(b) 28,26,28
26    cdabx= 0.0
      return
28    cdabx= abs(b)
      return
30    cdabx= abs(a)
      return

 end function cdabx

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 double precision function cdmpyr(a,b,c,d)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      cdmpyr = a*c-b*d
      return
 end function cdmpyr

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 double precision function cdmpyi(a,b,c,d)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      cdmpyi = b*c+a*d
      return
 end function cdmpyi

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 double precision function cddvdr(a,b,c,d)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      e = c*c+d*d
      f=a*c+b*d
      cddvdr=f/e
      return

 end function cddvdr

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 double precision function cddvdi(a,b,c,d)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      e = c*c+d*d
      f=b*c-a*d
      cddvdi=f/e
      return
 end function cddvdi

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine genbsl(ith,js)

!     Generate Bessel and Neumann functions for real arguments.
!     The index on the function is incremented by one.

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      common dtr,rtd,cpi
      common/mtxcom/nrank,nranki,a(80,80,2),b(80,80,2),cmxnrm(80)
      common /thtcom/ theta,sinth,costh
      common /bdycom/ dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib
      common/fnccom/pnmllg(81),bsslsp(81,31,3),cneumn(81,31,3), &
      bslkpr(81,31,3),bslkpi(81,31,3),cneumr(81,31,3),cneumi(81,31,3)

!     Set up a loop to get 2 successive Bessel functions
      nval=nrank-1
      pckr=ckr
      do 40 i=1,4
      call bessel(nval,pckr,answr,ierror)
      if(ierror) 20,20,32
20    ansa=answr
      nval=nval+1
      call bessel(nval,pckr,answr,ierror)
      if(ierror) 24,24,28
24    ansb=answr
      go to 60
28    nval=nval-1
32    nval=nval+nrank
40    continue

!     Program unable to generate Bessel functions
      if(debug_tmatrix) write(6,1001)
1001  format('unable to generate bessel functions')
!     Set up for proper recursion of the bessel functons
60    if(nval-nrank)100,100,64
64    iend=nval-nrank
      conn=2*(nval-1)+1.0
      do 72 ip=1,iend
      ansc=conn*ansa/pckr-ansb
      conn=conn-2.0
      ansb=ansa
      ansa=ansc
72     continue

!     Program is ready to recurse downward into Bessel function
100   bsslsp(nranki,ith,js)=ansb
      bsslsp(nranki-1,ith,js)=ansa
      conn=    (float(nrank+nrank-1))
      ie=nranki-2
      je=ie
      do 120 jb=1,je
      ansc=conn*ansa/pckr-ansb
      bsslsp(ie,ith,js)=ansc
      ansb=ansa
      ansa=ansc
      ie=ie-1
      conn=conn-2.0
120   continue

!     Generate Neumann functions
      cskrx= cos(pckr)/pckr
      snkrx= sin(pckr)/pckr
      ckr2=pckr**2
      cmuln=3.0
      snsa=-cskrx
      snsb=-cskrx/pckr-snkrx
      cneumn(1,ith,js)=snsa
      cneumn(2,ith,js)=snsb
      do 280 i=3,nranki
      snsc=cmuln*snsb/pckr-snsa
      cneumn(i,ith,js)=snsc
      snsa=snsb
      snsb=snsc
      cmuln=cmuln+2.0
280   continue

!     Perform Wronskian test on orders 0 and 1 and orders nrank-1 and nrank
      quanbt=abs(ckr2*(bsslsp(2,ith,js)*cneumn(1,ith,js)- &
      bsslsp(1,ith,js)*cneumn(2,ith,js))-1.0)
      quannt=abs(ckr2*(bsslsp(nranki,ith,js)*cneumn(nrank,ith,js)- &
      bsslsp(nrank,ith,js)*cneumn(nranki,ith,js))-1.0  )
      if (quanbt-1.0e-10) 360,352,352
  352 thtprt = rtd*theta
      if(debug_tmatrix) write(6,356) thtprt,pckrr,quanbt,quannt
356   format('theta= ',f9.4,' kr= ',f10.4,' bessel test= ',d12.5, &
      ' neumann test= ',d12.5)
      go to 362
  360 if (quannt-1.0e-10) 362,352,352
  362 return

 end subroutine genbsl

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine bessel(norder,argmnt,answr,ierror)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      ierror=0
      n=norder
      x=argmnt
      cn=n
      sum=1.0
      apr=1.0
      topr=-0.5*x*x
      ci=1.0
      cni=float(2*n) +3.0
      cni=float(2*n) +3.0
      do 60 i=1,100
      acr=topr*apr/(ci*cni)
      sum=sum+acr
      if( abs(acr/sum)-1.0e-20) 100,100,40
40    apr=acr
      ci=ci+1.0
      cni=cni+2.0
60    continue
      ierror=1
      go to 200
!     the series has converged
100   prod=    (float(2*n))+1.0
      fact=1.0d0
      if(n) 160,160,120
120    do 140 ifct=1,n
      fact=fact*x/prod
      prod=prod-2.0
140   continue
160   answr=fact*sum
200   return

 end subroutine bessel

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 double precision function crootr(a,b)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      dmag=(a*a+b*b)**0.25
      angle=0.5*atan2(b,a)
      crootr=dmag*cos(angle)
      return

 end function crootr

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 double precision function crooti(a,b)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      dmag=(a*a+b*b)**0.25
      angle=0.5*atan2(b,a)
      crooti=dmag*sin(angle)
      return
 end function crooti

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine prcssm(a,b,nr,nri)

!     subroutine to solve the equation t = (a-inverse)*b  (All matrices
!     are transposed) using Gauss-Jordan elimination.

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      double precision a(80,80,2),b(80,80,2)
      double precision aijmax(2),arat(2)
      equivalence (l,fl),(k,fk)
      n = 2*nr

!     Start reduction of the a matrix.
      do 80 i = 1,n
!     Search for the maximum element in the ith row of the a-matrix.
      aijmax(1) = a(i,1,1)
      aijmax(2) = a(i,1,2)
      jmax = 1
      do 10 j = 2,n
      if(cdabx(a(i,j,1),a(i,j,2)).le.cdabx(aijmax(1),aijmax(2))) goto 10
      aijmax(1) = a(i,j,1)
      aijmax(2) = a(i,j,2)
      jmax = j
   10 continue

!     If aijmax is zero ( as it will be for any row (or column) where the
!     index m is .gt. the index n, i.e., the legendre functions force those
!     matrix elements to zero),then the matrix is singular so solve the
!     reduced matrix (order = 2*(nrank-m)).

      if (cdabx(aijmax(1),aijmax(2)).gt.0.0  ) go to 20
      jmax = i
      go to 75

!     Normalize the ith row by aijmax (jmax element of the ith row).
   20 do 30 j = 1,n
      t1 = a(i,j,1)
      t2 = a(i,j,2)
      a(i,j,1) = cddvdr(t1,t2,aijmax(1),aijmax(2))
      a(i,j,2) = cddvdi(t1,t2,aijmax(1),aijmax(2))
!     normalize the ith row of b.
      t1 = b(i,j,1)
      t2 = b(i,j,2)
      b(i,j,1) = cddvdr(t1,t2,aijmax(1),aijmax(2))
      b(i,j,2) = cddvdi(t1,t2,aijmax(1),aijmax(2))
   30 continue

!     Use row transformations to get zeros above and below the jmax
!     element of the ith row of a.  Apply same row transformations
!     to the b matrix.
      do 70 k = 1,n
      if(k.eq.i) go to 70
      arat(1) = -a(k,jmax,1)
      arat(2) = -a(k,jmax,2)
      do 50 j = 1,n
      if (cdabx(a(i,j,1),a(i,j,2)).le.0.0  ) go to 50
      a(k,j,1) = cdmpyr(arat(1),arat(2),a(i,j,1),a(i,j,2))+a(k,j,1)
      a(k,j,2) = cdmpyi(arat(1),arat(2),a(i,j,1),a(i,j,2))+a(k,j,2)
   50 continue
      a(k,jmax,1)=0.0
      a(k,jmax,2)=0.0
      do 60 j=1,n
      if (cdabx(b(i,j,1),b(i,j,2)).le.0.0  ) go to 60
      b(k,j,1) = cdmpyr(arat(1),arat(2),b(i,j,1),b(i,j,2))+b(k,j,1)
      b(k,j,2) = cdmpyi(arat(1),arat(2),b(i,j,1),b(i,j,2))+b(k,j,2)
   60 continue
   70 continue

!     Store row counter (i) in top element of jmax column.  Thus,
!     the top row of a will contain the location of the pivot
!     (unity) element of each column (after reduction).
   75 l = i
!     Store the integer i in the top row of a.
      a(1,jmax,1) = fl
   80 continue

!     The reduction of a is complete.  Perform row interchanges
!     as indicated in the first row of a.
      do 120 i = 1,n
      k=i
!     Put the integer value in a into k.
   90 fk = a(1,k,1)
      if (k-i) 90,120,100

!     If k(1,i) is less than i, then that row has already been
!     involved in an interchange, and we use k(1,k) until we get
!     a value of k greater than i (corresponding to a row stored
!     below the ith row).
  100 do 110 j=1,n
      arat(1) = b(i,j,1)
      arat(2) = b(i,j,2)
      b(i,j,1) = b(k,j,1)
      b(i,j,2) = b(k,j,2)
      b(k,j,1) = arat(1)
      b(k,j,2) = arat(2)
  110 continue
  120 continue

      return
 end subroutine prcssm

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine addprc

!     subroutine to obtain the scattered field coefficients and calculate
!     the differential scattering cross section in the azimuthal plane.

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      complex ci,cim
      double precision fff
      common dtr,rtd,cpi
      common/mtxcom/nrank,nranki,a(80,80,2),tmat(80,80,2),cmxnrm(80)
      common/fnccom/pnmllg(81),bsslsp(81,31,3),cneumn(81,31,3), &
      bslkpr(81,31,3),bslkpi(81,31,3),cneumr(81,31,3),cneumi(81,31,3)
      common /bdycom/ dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib
      common /cmvcom/ nm,kmv,cmi(20),cmv,cm2,twm,prodm
      common /thtcom/ theta,sinth,costh
      common /uvccom/anginc,acans(181,2,2),uang(181),rtsfct,dltang,nuang
      common/scatt/summ1,summ2
!      common/vivek2/twoa
      common /vivek2/dpart
      double precision zxold(181),zyold(181),ab1(80,2),ab2(80,2), &
      fg1(80,2),fg2(180,2),fgans(181,2,2)
      logical test
      data test/.true./
      ci = (0.0,1.0)

!     Generate the Legendre functions for the incident angle.
      if(anginc) 15,5,15
5     costh=1.0
10     sinth=0.0
      theta=0.0
      go to 30
   15 if(anginc-180.0  ) 25,20,25
20    costh=-1.0
      go to 10
   25 theta = dtr*anginc
      sinth=sin(theta)
      costh=cos(theta)
   30 call genlgp
!     Generate the incident field coefficients -- ab1 = theta polarization
!     and ab2 = phi polarization.

      cn=0.0
      do 35 n=1,nrank
      np = n+nrank
      cn=cn+1.0
      n1 = n+1
      ci1r=dble(ci**n)
      ci1i=aimag(ci**n)
      ci2r=dble(ci**n1)
      ci2i=aimag(ci**n1)
      p1 = cn*costh*pnmllg(n1)-(cn+cmv)*pnmllg(n)
      p2 = cmv*pnmllg(n1)
      ab1(n,1) = -ci1r*p2
      ab1(n,2) = -ci1i*p2
      ab1(np,1) = ci2r*p1
      ab1(np,2) = ci2i*p1
      ab2(n,1) = ci1r*p1
      ab2(n,2) = ci1i*p1
      ab2(np,1) = -ci2r*p2
      ab2(np,2) = -ci2i*p2
35    continue

!     The scattered field coefficients = the row vector of incident field
!     coefficients times the t-transposed matrix.

      nr2 = 2*nrank
      do 45 j = 1,nr2
      s1r=0.0
      s1i=0.0
      s2r=0.0
      s2i=0.0
      do 40 i = 1,nr2
      s1r = s1r+cdmpyr(ab1(i,1),ab1(i,2),tmat(i,j,1),tmat(i,j,2))
      s1i = s1i+cdmpyi(ab1(i,1),ab1(i,2),tmat(i,j,1),tmat(i,j,2))
      s2r = s2r+cdmpyr(ab2(i,1),ab2(i,2),tmat(i,j,1),tmat(i,j,2))
      s2i = s2i+cdmpyi(ab2(i,1),ab2(i,2),tmat(i,j,1),tmat(i,j,2))
   40 continue
      fg1(j,1) = s1r
      fg1(j,2) = s1i
      fg2(j,1) = s2r
      fg2(j,2) = s2i
   45 continue

!     Calculate scattering cossections normalized for parallel and
!     perpendicular polarizations

      sum1=0.0
      sum2=0.0
      do 1001 i=1,nrank
        ii=i+nrank
        temp1=fg1(i,1)**2+fg1(i,2)**2+fg1(ii,1)**2+fg1(ii,2)**2
        temp2=fg2(i,1)**2+fg2(i,2)**2+fg2(ii,1)**2+fg2(ii,2)**2
        sum1=sum1+temp1/cmxnrm(i)
        sum2=sum2+temp2/cmxnrm(i)
1001  continue

!     Normalize scattering crossections
      sum1=(rtsfct*2.0/conk)*sum1
      sum2=(rtsfct*2.0/conk)*sum2

!     Normalize wrt equivalent spherical diameter
      cnorm=aovrb**(-2./3.)
      sum1=sum1*cnorm
      sum2=sum2*cnorm

!     Accumulate results for each m value
      summ1=sum1+summ1
      summ2=sum2+summ2


!     Evaluate the scattered field at each scattering angle.
      do 170 iu = 1,nuang
 if(debug_tmatrix) print*,uang(iu),'here1'
!     Generate the Legendre multipliers
      if(uang(iu)) 95,85,95
85    costh=1.0
90    sinth=0.0
      theta=0.0
      go to 110
   95 if(uang(iu)-180.0  ) 105,100,105
100   costh=-1.0
      go to 90
  105 theta = dtr*uang(iu)
      sinth=sin(theta)
      costh=cos(theta)
  110 call genlgp
      fgans(iu,1,1)=0.0
      fgans(iu,1,2)=0.0
      fgans(iu,2,1)=0.0
      fgans(iu,2,2)=0.0
      cn=0.0
      do 160 n = 1,nrank
      np = n+nrank
      n1 = n+1
      cn=cn+1.0
      p1 = cn*costh*pnmllg(n1)-(cn+cmv)*pnmllg(n)
      p2 = cmv*pnmllg(n1)
      cim = (-ci)**n1
      cir=dble(cim)
      cii=aimag(cim)
      f1r = fg1(n,1)*p2
      f1i = fg1(n,2)*p2
      g1r = -fg1(np,2)*p1
      g1i = fg1(np,1)*p1
      fgans(iu,1,1)= fgans(iu,1,1)+cdmpyr(cir,cii,f1r+g1r,f1i+g1i)/cmxnrm(n)
      fgans(iu,1,2)= fgans(iu,1,2)+cdmpyi(cir,cii,f1r+g1r,f1i+g1i)/cmxnrm(n)
      f2r = fg2(n,1)*p1
      f2i = fg2(n,2)*p1
      g2r = -fg2(np,2)*p2
      g2i = fg2(np,1)*p2
      fgans(iu,2,1)= fgans(iu,2,1)-cdmpyr(cir,cii,f2r+g2r,f2i+g2i)/cmxnrm(n)
      fgans(iu,2,2)= fgans(iu,2,2)-cdmpyi(cir,cii,f2r+g2r,f2i+g2i)/cmxnrm(n)
  160 continue

!     The normalized diff scat cross sect is given by ((8/ka)*fgans)**2
!     Scale fgans to calculate diff scat cross sect (rtsfct = 8/ka)

      fgans(iu,1,1) = rtsfct*fgans(iu,1,1)
      fgans(iu,1,2) = rtsfct*fgans(iu,1,2)
      fgans(iu,2,1) = rtsfct*fgans(iu,2,1)
      fgans(iu,2,2) = rtsfct*fgans(iu,2,2)
  170 continue

      if(debug_tmatrix) write(6,175) kmv,anginc
  175 format(10x,'accumulations for m= ',i3,' angle of inc= ',f6.2)

!     Accumulate the results for each m value.
      do 172 iup = 1,nuang
      acans(iup,1,1) = acans(iup,1,1)+fgans(iup,1,1)
      acans(iup,1,2) = acans(iup,1,2)+fgans(iup,1,2)
      acans(iup,2,1) = acans(iup,2,1)+fgans(iup,2,1)
      acans(iup,2,2) = acans(iup,2,2)+fgans(iup,2,2)

      if(debug_tmatrix) write(6,202) uang(iup),acans(iup,1,1),acans(iup,1,2)
 202  format('parallel      ',f7.3,2x,2d15.7)
      if(debug_tmatrix) write(6,203) uang(iup),acans(iup,2,1),acans(iup,2,2)
 203  format('perpendicular ',f7.3,2x,2d15.7)
  172 continue

!     Calculate the extinction crossections
      extpp=acans(1,1,2)*rtsfct/4.0
      extper=acans(1,2,2)*rtsfct/4.0

!     Normalize wrt equivelent spherical diameter
      extpp=extpp*cnorm
      extper=extper*cnorm

!     Calculate forward and backward amplitude in far zone
!     sigma equals 4.0/k
      forrp=sigma*acans(1,1,1)/rtsfct
      forip=sigma*acans(1,1,2)/rtsfct
      forpe=sigma*acans(1,2,1)/rtsfct
      foripe=sigma*acans(1,2,2)/rtsfct
      borrp=sigma*acans(nuang,1,1)/rtsfct
      borip=sigma*acans(nuang,1,2)/rtsfct
      borrpe=sigma*acans(nuang,2,1)/rtsfct
      boripe=sigma*acans(nuang,2,2)/rtsfct

!     Calculate normalized radar crossections for both polarizations
      xhor=acans(nuang,1,1)**2+acans(nuang,1,2)**2
      yver=acans(nuang,2,1)**2+acans(nuang,2,2)**2

!     Normalize wrt equivalent spherical diameter
      xhor=xhor*cnorm
      yver=yver*cnorm

 if(debug_tmatrix) write(50,*)nrank,kmv,xhor,yver

!     Write the scattering results

if(debug_tmatrix) then
      write(6,2001)
      write(6,2002) summ1,extpp,xhor,forrp,forip,borrp,borip
      write(6,2003) summ2,extper,yver,forpe,foripe,borrpe,boripe
      write(6,*)
endif


!      if(kmv .eq. 4)then
!      fff=4./(cpi*dpart**2)
!      write(10,*)kmv
!      write(10,2030)anginc, dpart,forrp,forip, borrp, borip
!      write(10,2030)anginc, dpart,forpe,foripe,borrpe,boripe
if(debug_tmatrix)  then
      write(8,2001)
      write(8,2002) summ1,extpp,xhor,forrp,forip,borrp,borip
      write(8,2003) summ2,extper,yver,forpe,foripe,borrpe,boripe
      write(8,*)
endif


!      write(8,2020) forrp,forip,borrp,borip
!      write(8,2020) forpe,foripe,borrpe,boripe

2001  format(18x,'scat',9x,'ext',10x,'radar',8x,'re(forw)',5x,'im(forw)' &
       ,5x,'re(back)',5x,'im(back)')
2002  format('parallel      ',7(d12.5,1x))
2003  format('perpendicular ',7(d12.5,1x))
2030  format(f4.0,f5.2,8(d12.5,1x))
2020  format(4d15.7)

      return
 end subroutine addprc

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine printm(p,n,nd)
!     subroutine to print a matrix

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      double precision p(nd,nd,2)
      do 20 k=1,2
        do 10 i=1,n
          if(debug_tmatrix) write(6,100) i,(p(i,j,k),j=1,n)
   10   continue
   20 continue
  100 format (/,' rowi3',2x,1p8d15.6,/,(10x,8d15.6))

      return
 end subroutine printm

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine genkr
!     calculate ckr and dckr as a function of theta for a oblate spheroid

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      common dtr,rtd,cpi
      common /bdycom/ dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib
      common /thtcom/ theta,sinth,costh

      bovra=1.0/aovrb
      qb = 1.000/sqrt((bovra*costh)**2+sinth**2)
      ckr = conk*qb
      dckr = conk*costh*sinth*(bovra**2-1.000)*qb**3
!       thet1=theta*rtd
 10     format(2x,3(d15.7,4x))

      return

 end subroutine genkr

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine calenp

!toshii: i don't understand this subroutine... why epps dimension is just 2...?

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

 implicit double precision (a-h)
 implicit double precision (o-z)
 implicit integer (i-n)

!     *****

 common dtr,rtd,cpi
 common/endpnt/epps(4),nsect
 common/bdycom/dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib

 epps(1)=0.
 epps(2)=cpi/2.

 return
 end  subroutine calenp

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine epslon(alam,temp,reps,aimeps)

!     This subroutine calculates the complex dielectric function for
!     water at precipitation and cloud radar wavelengths.
!     The equations are taken from Ray, P. S., 1972: Broadband complex
!     refractive indices of ice and water.  Applied Optics, 11(8), 1836-1844.

!     Note Temperature:  strictly applicable for water : -20C to 50C


!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

      sigma=12.5664d08
      cpi=3.141592
      einf=5.27137+.0216474*temp-.00131198*(temp**2)
      alpha=-16.8129/(temp+273.)+.0609265
      brin=2513.98d0/(temp+273.000)
      als=.00033836*exp(brin)
      t=temp-25.0
      eps=78.54*(1.-4.579e-03*t+1.19e-05*(t**2)-2.8e-08*(t**3))
!     set up common parameters.
      slam=(als/alam)**(1.0-alpha)
      sinal=sin(alpha*cpi/2.000)
      cosal=cos(alpha*cpi/2.000)
!     calculate re(epslon)
      reps=(eps-einf)*(1.0+slam*sinal)
      x=1.0+2.0*slam*sinal+slam**2
      reps=reps/x
      reps=reps+einf
!     calculate im(epslon)
      aimeps=(eps-einf)*slam*cosal
      aimeps=(aimeps/x)+sigma*alam/18.8496d10
      return

 end subroutine epslon

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine spice(alam,temp,dcnr,dcni,dcnr2,dcni2)
!      subroutine to compute the complex relative dielectric const of ice
!      solid ice

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

       call epslon(alam,temp,dcnr,dcni)
       if(debug_tmatrix) write(6,2) dcnr,dcni
 2     format(1x,'dielectric const of water= ',d15.8,' + j ',d15.8)
       emx= reps
       emy= -aimeps

!     dielectric constant of air/ice mixture with rho=600 kg/m^-3
!     hardwire from debye theory where solid ice was obtained from
!     Ray (1982) with lambda=5.5cm and T=-10C  (not strong function of
!     either however).
!      dcnr2=2.135
!      dcni2=0.0003210


!      Assign the inclusion epsilon.
!      Value taken from Tiuri, Trans IEEE, GE-20, page 52
!      Compute the real part of epsilon.
       dcnr2=1.508**3
!      Compute the imaginary part of epsilon.
       dcni2=(0.34*0.001)/((1.0-0.417)**2.0)

       if(debug_tmatrix) write(6,30) dcnr2,dcni2
 30    format(2x,'dielectric const of ice= ',2(d15.7,4x))
       return

 end subroutine spice

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine genbkr(ar,ai,iswt,ith,js)
!     4th january 1989

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

       common/mtxcom/nrank,n,a(80,80,2),b(80,80,2),cmxnrm(80)
       common/fnccom/pnmllg(81),bsslsp(81,31,3),cneumn(81,31,3), &
       bslkpr(81,31,3),bslkpi(81,31,3),cneumr(81,31,3),cneumi(81,31,3)
       double precision sjr(2000),sji(2000),syr(100),syi(100),&
!       shr(1), shi(1)
!toshii it should be 2
       shr(2), shi(2)
    
!      logical lt,lf
!      write(6,'(a6,2d15.7,/)') 'newarg',ar,ai
       nak=2
!     write (5,104)
104   format ('spherical bessel functions for complex argument')
      zero=0.0d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      iz=0
      dr=ar*ar-ai*ai
        di=two*ar*ai
       cc=two
       eps=1.0d-16
       wunr=one
       wuni=zero
        call dvdd(wunr,wuni,dr,di,t1,t2)
!         write(6,*) 'aft csp dv1'
        srarg=sqrt(ar*ar+ai*ai)
        if(srarg.gt.0.5d0)go to 29
        np=n+1
        call mltd(ar,ai,ar,ai,zr,zi)
        zr=zr/two
        zi=zi/two
        fdnm=three
        hdn=one
        hdnm=one
        hdni=zero
        do 14 i=1,np
         nn=i-1
        en=nn
!       calculate..
        if(nn-1)2,6,3
6       fnr=ar/three
        fni=ai/three
        go to 5
2       fnr=one
        fni=zero
        go to 5
3       call mltd(fnr,fni,ar,ai,fnr,fni)
        fdnm=fdnm+two
        fnr=fnr/fdnm
        fni=fni/fdnm
5       ansr=one
        ansi=zero
        pansr=one
        pansi=zero
        trm=-one
        tim=zero
        ab=one
        ba=three
7       gnu=ab*(two*en+ba)
        zrs=-zr/gnu
        zis=-zi/gnu
        call mltd(trm,tim,zrs,zis,trm,tim)
        ansr=ansr-trm
        ansi=ansi-tim
        if(ansr.eq.zero)go to 15
        if (ansi.eq.zero)go to 16
      if (abs((pansr-ansr)/ansr).le.eps.and.abs((pansi-ansi)/ansi) &
          .le.eps)go to 8
        go to 17
15      if (abs((pansi-ansi)/ansi).le.eps)go to 8
        go to 17
16      if (abs((pansr-ansr)/ansr).le.eps)go to 8
17      pansr=ansr
        pansi=ansi
        ab=ab+one
        ba=ba+two
        go  to 7
8       call mltd(fnr,fni,ansr,ansi,sjr(i),sji(i))
!        calculate...
        if(nn-1)4,10,9
4       gdr=-one
        gdi=zero
        call dvdd(gdr,gdi,ar,ai,hr,hi)
!        write(6,*) 'aft csp dv2'
        go to 11
10      hdr=ar
        hdi=ai
9       call mltd(hdr,hdi,ar,ai,hdr,hdi)
        hdnm=hdnm*hdn
        hdn=hdn+two
        call dvdd(hdnm,hdni,hdr,hdi,hr,hi)
!          write(6,*) 'aft csp dv3'
        hr=-hr
        hi=-hi
11      alsr=one
        alsi=zero
        palsr=one
        palsi=zero
        trn=-one
        tin=zero
        ac=one
        ca=one
12      hnu=ac*(ca-two*en)
        xrs=-zr/hnu
        xis=-zi/hnu
        call mltd(trn,tin,xrs,xis,trn,tin)
        alsr=alsr-trn
        alsi=alsi-tin
        if(alsr.eq.zero)go to 18
        if(alsi.eq.zero)go to 19
      if(abs((palsr-alsr)/alsr).le.eps.and.abs((palsi-alsi)/alsi) &
            .le.eps)go to 13
        go to 20
18      if(abs((palsi-alsi)/alsi).le.eps)go to 13
       go to 20
19       if (abs((palsr-alsr)/alsr).le.eps)go to 13
20       palsr=alsr
         palsi=alsi
         ac=ac+one
         ca=ca+two
         go to 12
13       call mltd(hr,hi,alsr,alsi,syr(i),syi(i))
!
!
!
!
!
         if(nak.eq.2)go to 14
         if(nak.eq.5)go to 50
         yrr=syr(i)
         yii=syi(i)
         if(ai.lt.zero)go to 51
         syr(i)=sjr(i)-yii
         syi(i)=sji(i)+yrr
         goto 14
51      syr(i)=sjr(i)+yii
        syi(i)=sji(i)-yrr
        go to 14
50        if(ai.lt.zero)go to 48
          shr(i)=sjr(i)-syi(i)
          shi(i)=sji(i)+syr(i)
          go to 14
48        shr(i)=sjr(i)+syi(i)
        shi(i)=sji(i)-syr(i)
14      continue
        do 545 ka=1,n
        bslkpr(ka,ith,js)=sjr(ka)
        bslkpi(ka,ith,js)=sji(ka)
        cneumr(ka,ith,js)=syr(ka)
        cneumi(ka,ith,js)=syi(ka)
545     continue
        return
29      dsn=sin(ar)
        dcs=cos(ar)
        exyl=exp(ai)
        exys=exp(-ai)
        xsn=ar*dsn
        xco=ar*dcs
        ysn=ai*dsn
        yco=ai*dcs
        zxy=ar*ar+ai*ai
        tzxy=two*zxy
        sjzrl=(xsn+yco)/tzxy
        sjzrs=(xsn-yco)/tzxy
        sjzil=(xco-ysn)/tzxy
        sjzis=(-xco-ysn)/tzxy
        syzrl=exyl*(-sjzil)
        syzrs=exys*sjzis
        syzil=exyl*sjzrl
        syzis=exys*(-sjzrs)
        sjzrl=exyl*sjzrl
        sjzrs=exys*sjzrs
        sjzil=exyl*sjzil
        sjzis=exys*sjzis
        sjr(1)=sjzrl+sjzrs
        sji(1)=sjzil+sjzis
       sjr(2)=((ar*sjzrl+ai*sjzil)/zxy+syzrl)+ &
           ((ar*sjzrs+ai*sjzis)/zxy+syzrs)
       sji(2)=((-ai*sjzrl+ar*sjzil)/zxy+syzil)+ &
            ((-ai*sjzrs+ar*sjzis)/zxy+syzis)
       nho=0
       if(abs(ai).lt.5.0d0)go to 43
!      .....
       nho=1
       yex=exp(-abs(ai))
       anur=yex*dsn
       anui=yex*dcs
       if(ai.ge.zero)anui=-anui
       call  dvdd(anur,anui,ar,ai,hrz,hiz)
!       write(6,*) 'aft csp dv4'
       call  mltd(ar,ai,ar,ai,zsr,zsi)
       call  dvdd(anur,anui,zsr,zsi,hrw,hiw)
!          write(6,*) 'aft csp dv5'
       if(ai)38,39,39
38      anur=-anur
        go to 40
39       anui=-anui
40      call  dvdd(anui,anur,ar,ai,hoa,hob)
!          write(6,*) 'aft csp dv6'
! ....
! ....
! ....
         if(nak.lt.5)go to 54
         shr(1)=hrz
         shi(1)=hiz
         shr(2)=hrw-hoa
         shi(2)=hiw-hob
       go to 55
54     if (nak.eq.2)go to 56
       syr(1)=hrz
       syi(1)=hiz
        syr(2)=hrw-hoa
        syi(2)=hiw-hob
        go to 36
56      hrw=hrw-hoa
        hiw=hiw-hob
        syr(1)=-sji(1)+hiz
        syi(1)=sjr(1)-hrz
        syr(2)=-sji(2)+hiw
        syi(2)=sjr(2)-hrw
        go to 57
55      syr(1)=-sji(1)+shi(1)
        syi(1)=sjr(1)-shr(1)
        syr(2)=-sji(2)+shi(2)
        syi(2)=sjr(2)-shr(2)
57      if(ai.ge.zero)go to 36
        syr(1)=-syr(1)
        syi(1)=-syi(1)
        syr(2)=-syr(2)
        syi(2)=-syi(2)
       go to 36
! ....
43      syr(1)=syzrl+syzrs
        syi(1)=syzil+syzis
        syr(2)=((ar*syzrl+ai*syzil)/zxy-sjzrl)+ &
            ((ar*syzrs+ai*syzis)/zxy-sjzrs)
      syi(2)=((-ai*syzrl+ar*syzil)/zxy-sjzil)+ &
          ((-ai*syzrs+ar*syzis)/zxy-sjzis)
!     ...
!     ...
!     ...
42      if(nak.eq.2)go to 36
        if(nak.eq.5)go to 52
        yrz=syr(1)
        yiz=syi(1)
        yrw=syr(2)
        yiw=syi(2)
        if(ai.lt.zero)go to 53
        syr(1)=sjr(1)-yiz
        syi(1)=sji(1)+yrz
        syr(2)=sjr(2)-yiw
        syi(2)=sji(2)+yrw
        go to 36
53      syr(1)=sjr(1)+yiz
        syi(1)=sji(1)-yrz
        syr(2)=sjr(2)+yiw
        syi(2)=sji(2)-yrw
        go to 36
52      if(ai.lt.zero)go to 41
        shr(1)=sjr(1)-syi(1)
        shi(1)=sji(1)+syr(1)
        shr(2)=sjr(2)-syi(2)
        shi(2)=sji(2)+syr(2)
        go to 36
41      shr(1)=sjr(1)+syi(1)
        shi(1)=sji(1)-syr(1)
        shr(2)=sjr(2)+syi(2)
        shi(2)=sji(2)-syr(2)
36      if(n.gt.1)goto 111
        do 551 kb=1,n
        bslkpr(kb,ith,js)=sjr(kb)
        bslkpi(kb,ith,js)=sji(kb)
        cneumr(kb,ith,js)=syr(kb)
        cneumi(kb,ith,js)=syi(kb)
551     continue
        return
!       .....
111     m=n+1
!       .....
        nn=srarg+30
        if ((n+30).gt.nn)nn=n+30
        gdr=sjr(2)
        gdi=sji(2)
30      sjr(nn)=zero
        sji(nn)=zero
        sjr(nn-1)=1.0d-20
        sji(nn-1)=1.0d-20
        nm=nn-2
        do 31 k=2,nm
        kk=nn-k
        call dvdd(sjr(kk+1),sji(kk+1),ar,ai,sjr(kk),sji(kk))
!             write (6,*) 'aft csp dv7'
!       call errset(72,lt,lf,lf,lf,)
!       call errsns
        sjr(kk)=(cc*kk+one)*sjr(kk)-sjr(kk+2)
!       call errsns(num,,,)
!       if(num.eq.72)go to 24
!       call errsns
        sji(kk)=(cc*kk+one)*sji(kk)-sji(kk+2)
!     call errset(72,lf,lf,lf,lt,)
!     call errsns(num,,,)
!     if(num.eq.72)go to 24
31    continue
      call dvdd(gdr,gdi,sjr(2),sji(2),rar,rai)
!          write (6,*) 'aft csp dv8'
! ...
! ...
      if(rar.ne.zero.and.rai.ne.zero)go to 67
      if(abs(sjr(2)).lt.abs(sji(2)))go to 68
      if(rar.ne.zero)go to 69
      if(gdr.eq.zero.and.sji(2).eq.zero)go to 69
      go to 24
69      if(rai.ne.zero)go to 67
        if(gdi.eq.zero.and.sji(2).eq.zero)go to 67
        go to 24
68      if(rar.ne.zero)go to 70
        if(gdi.eq.zero.and.sjr(2).eq.zero)go to 70
        go to 24
70     if(rai.ne.zero)go to  67
         if(gdr.eq.zero.and.sjr(2).eq.zero)go to 67
         go to 24
67       do 32 k=3,m
         tr=sjr(k)
         ti=sji(k)
32       call mltd(tr,ti,rar,rai,sjr(k),sji(k))
         sjr(2)=gdr
         sji(2)=gdi
! ..
         if(nho.eq.1)go to 44
!
!
!
         if(nak.eq.3)go to 66
22       do 23 k=3,m
       call dvdd(syr(k-1),syi(k-1),ar,ai,syr(k),syi(k))
!          write (6,*) 'aft csp dv9'
       syr(k)=(cc*k-three)*syr(k)-syr(k-2)
       syi(k)=(cc*k-three)*syi(k)-syi(k-2)
        if(nak.eq.2)go to 23
         if(ai.lt.zero)go to 45
47       shr(k)=sjr(k)-syi(k)
         shi(k)=sji(k)+syr(k)
        go to 23
45      shr(k)=sjr(k)+syi(k)
        shi(k)=sji(k)-syr(k)
23      continue
        do 555 kc=1,n
        bslkpr(kc,ith,js)=sjr(kc)
        bslkpi(kc,ith,js)=sji(kc)
        cneumr(kc,ith,js)=syr(kc)
        cneumi(kc,ith,js)=syi(kc)
555     continue
        return
66      do 60 k=3,m
        call dvdd(yrw,yiw,ar,ai,yrt,yit)
!            write(6,*) 'aft csp dv10'
        yrt=(cc*k-three)*yrt-yrz
         yit=(cc*k-three)*yit-yiz
         if(ai.lt.zero)go to 58
        syr(k)=sjr(k)-yit
        syi(k)=sji(k)-yrt
        go to 59
58      syr(k)=sjr(k)+yit
        syi(k)=sji(k)-yrt
59     yrz=yrw
       yiz=yiw
       yrw=yrt
       yiw=yit
60     continue
        do 565 kd=1,n
        bslkpr(kd,ith,js)=sjr(kd)
        bslkpi(kd,ith,js)=sji(kd)
        cneumr(kd,ith,js)=syr(kd)
        cneumi(kd,ith,js)=syi(kd)
565     continue
       return
!
!
!
44     if(nak.ne.5)go to 61
       do 46 k=3,m
       call dvdd(shr(k-1),shi(k-1),ar,ai,shr(k),shi(k))
!           write (6,*) 'aft csp dv11'
       shr(k)=(cc*k-three)*shr(k)-shr(k-2)
        shi(k)=(cc*k-three)*shi(k)-shi(k-2)
        syr(k)=-sji(k)+shi(k)
        syi(k)=sjr(k)-shr(k)
        if(ai.ge.zero)go to 46
        syr(k)=-syr(k)
        syi(k)=-syi(k)
46        continue
        do 575 ke=1,n
        bslkpr(ke,ith,js)=sjr(ke)
        bslkpi(ke,ith,js)=sji(ke)
        cneumr(ke,ith,js)=syr(ke)
        cneumi(ke,ith,js)=syi(ke)
575     continue
          return
61      if(nak.eq.3)go to 62
        do 63 k=3 ,m
        call dvdd(hrw,hiw,ar,ai,hrt,hit)
!           write(6,*) 'aft csp dv12'
        hrt=(cc*k-three)*hrt-hrz
        hit=(cc*k-three)*hit-hiz
      syr(k)=-sji(k)+hit
      syi(k)=sjr(k)-hrt
       if(ai.ge.zero)go to 64
        syr(k)=-syr(k)
        syi(k)=-syi(k)
64      hrz=hrw
        hiz=hiw
        hrw=hrt
        hiw=hit
63      continue
        do 585 kf=1,n
        bslkpr(kf,ith,js)=sjr(kf)
        bslkpi(kf,ith,js)=sji(kf)
        cneumr(kf,ith,js)=syr(kf)
        cneumi(kf,ith,js)=syi(kf)
585     continue
        return
62      do 65 k=3,m
        call dvdd(syr(k-1),syi(k-1),ar,ai,syr(k),syi(k))
!           write (6,*) 'aft csp dv13'
        syr(k)=(cc*k-three)*syr(k)-syr(k-2)
        syi(k)=(cc*k-three)*syi(k)-syi(k-2)
65      continue
        do 595 kh=1,n
        bslkpr(kh,ith,js)=sjr(kh)
        bslkpi(kh,ith,js)=sji(kh)
        cneumr(kh,ith,js)=syr(kh)
        cneumi(kh,ith,js)=syi(kh)
595     continue
        return
24      nn=nn-1
        if(debug_tmatrix) write(6,26)nn
26      format (1x,'nn reduced to',i6)
        iz=iz+1
        if(iz.le.25)goto 211
        do 591 kg=1,n
        bslkpr(kg,ith,js)=sjr(kg)
        bslkpi(kg,ith,js)=sji(kg)
        cneumr(kg,ith,js)=syr(kg)
        cneumi(kg,ith,js)=syi(kg)
591     continue
        return
211     go to 30
 end subroutine genbkr

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine dvdd(xa,ya,xb,yb,xc,yc)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****

!
!       implicit real*8 (a-h,o-z)
!      logical lt,lf
!      data lt/.true./,lf/.false./
       zero=0.0d0
       if(xb.ne.zero.or.yb.ne.zero)go to 3
       if(debug_tmatrix) write(6,100)
 100   format('both real and imaginary parts of denominator are zero')
       return
!        write(6,500)
!        format('error 1 in dvdd')
!      call errset(73,lt,lf,lf,lf,)
!      call errset(74,lt,lf,lf,lf,)
3      denom=xb*xb+yb*yb
!         if(denom.eq.zero) return
       if(denom.eq.zero)go to 1
       xx=(xa/denom)*xb+(ya/denom)*yb
       if(xx.eq.zero)go to 1
       yc=ya*(xb/denom)-(xa/denom)*yb
       if(yc.eq.zero)go to 1
       xc=xx
        return
!            write(6,501)
!        format('error 3 in dvdd')
!       call errset(73,lf,lf,lf,lt,)
!       call errset(74,lf,lf,lf,lt,)
 1        if(debug_tmatrix)  write(6,*) 'overflow,denom=0'
       if(abs(xb).lt.abs(yb))go to 2
    8  dc=yb/xb
        ac=xa/xb
        bc=ya/xb
!       call errset(74,lt,lf,lf,lf,)
        denom=1.0d0+dc*dc
!
       xc=(ac+bc*dc)/denom

       yc=(bc-ac*dc)/denom
!      call errset(74,lf,lf,lf,lt,)
       return
2     ad=xa/yb
       cd=xb/yb
       bd=ya/yb
!      call errset(74,lt,lf,lf,lf,)
       denom=1.0d0+cd*cd
!
       xc=(bd+ad*cd)/denom
       yc=(-ad+bd*cd)/denom
!      call errset(74,lf,lf,lf,lt,)
       return
 end subroutine dvdd

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine mltd(xa,ya,xb,yb,xc,yc)

!     *****
!     Declare all a-h, o-z variables as implicity double precision real
!     Carey 9/15/95

      implicit double precision (a-h)
      implicit double precision (o-z)
      implicit integer (i-n)

!     *****


!      implicit real*8 (a-h,o-z)
       xx=xa*xb-ya*yb
       yc=xa*yb+ya*xb
       xc=xx
   return
 end subroutine mltd

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

!---------------------------------------------------------------
!      Determines the hailstone axis ratio (a/b) using a look-up
!      table based on hailstone diameter from Knight (1986).
!      Note that average a/b is used and 95% confidence level ranges
!      are not utilized.
!      ldc 11/03/95

!      dpart is in cm here

 subroutine hailab(dpart,aovrb)

       double precision dpart, aovrb

       if(dpart.lt.0.1)then
         aovrb=0.99
       elseif(dpart.le.0.5)then
         aovrb=0.86
       elseif(dpart.le.1.0)then
         aovrb=0.85
       elseif(dpart.le.1.5)then
         aovrb=0.85
       elseif(dpart.le.2.0)then
         aovrb=0.81
       elseif(dpart.le.2.5)then
         aovrb=0.78
       elseif(dpart.le.3.0)then
         aovrb=0.74
       elseif(dpart.le.3.5)then
         aovrb=0.76
       elseif(dpart.le.4.0)then
         aovrb=0.84
       elseif(dpart.le.4.5)then
         aovrb=0.79
       elseif(dpart.le.5.0)then
         aovrb=0.76
       elseif(dpart.le.5.5)then
         aovrb=0.71
       elseif(dpart.le.6.0)then
         aovrb=0.70
       else
         aovrb=0.80
       endif

!      note that the aovrb for d<1mm and d>60mm are made up
!      with an attempt to use known info.  Pea size hail
!      should be nearly spherical and according to Knight
!      (1986) extremely large hail tend back toward a
!      spherical shape so 0.8 was a compromise.  Most likely
!      calculations should be 1mm < D < 60 mm anyway.

       return
 end subroutine hailab

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine ray(ph,lm,t,er,ei)

!     Carey 6/97
!
!     This program calculates the complex dielectric function for ice
!     and water at precipitation and cloud radar wavelengths.
!     The equations are taken from Ray, P. S., 1972: Broadband complex
!     refractive indices of ice and water.  Applied Optics, 11(8), 1836-1844.

!     Input:  1. water vs. ice
!             2. Wavelength
!             3. Temperature:  strictly applicable for ice : -20C to 0C
!                                                      water: -20C to 50C

!     Output: Dielectric function (real,imaginary)

      real pi
      parameter (pi=3.1415927)

      integer ph
      double precision t,lm
      double precision einf,alph,lams,sigm,epss
      double precision er,ei

!     set functions for water or ice respectively

      if(ph.eq.0)then
!       water

        epss=78.54*(1.0-4.579e-3*(t-25.0)+1.19e-5*(t-25.0)**2- &
             2.8e-8*(t-25.0)**3)

        sigm=12.5664e8

        lams=0.00033836*exp(2513.98/(t+273))

        alph=-16.8129/(t+273)+0.0609265

        einf=5.27137+0.0216474*t-0.00131198*t**2

      else
!       ice

        epss=203.168+2.5*t+0.15*t**2

        sigm=1.26*exp(-12500.0/((t+273.0)*1.9869))

        lams=9.990288e-4*exp(13200.0/((t+273.0)*1.9869))

        alph=0.288+0.0052*t+0.00023*t**2

        einf=3.168

      endif

!     calculate the dielectric function eps = er - i*ei

!     real

      er=einf+((epss-einf)*(1+(lams/lm)**(1-alph)*sin(alph*pi/2))) &
        / (1+2*(lams/lm)**(1-alph)*sin(alph*pi/2)+(lams/lm)** &
        (2*(1-alph)))

!     imaginary

      ei=((epss-einf)*(lams/lm)**(1-alph)*cos(alph*pi/2))/ &
         (1+2*(lams/lm)**(1-alph)*sin(alph*pi/2)+(lams/lm)** &
         (2*(1-alph)))+(sigm*lm)/18.8496e10

      if(ph.eq.0)then
        if(debug_tmatrix) write(6,*) 'liquid water'
      else
        if(debug_tmatrix) write(6,*) 'ice'
      endif
if(debug_tmatrix) then
      write(6,*) 't = ',t,' c'
      write(6,*) 'lambda = ',lm,' cm'
      write(6,*)
      write(6,*) 'dielectric function:'
      write(6,*) 'real component = ',er
      write(6,*) 'imaginary component =',ei
endif

      return

 end subroutine ray

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 subroutine debye(rho,eprsol,episol,epr,epi)

!     Debye Theory is used to determine the relative permittivity/dielectric
!     of an air/ice mixture given the density of this mixture.  See Battan
!     (1973) p. 40.
!     LDC 3/21/98


      complex m2, m2sol, msol, kovrrho, rho, rhosol
      double precision nsol,ksol,n,k,epr,epi,eprsol,episol

!     calculates refractive index and of an air/ice mixture

!     Initialization

!     Densitiy of solid ice
      rhosol=(0.917,0)

!     Complex index of refraction for water and solid ice (rho=0.917 g/cm^3)
!     is obtained from subroutine ray which gives relative permittivity.
!     Must convert first.
!     approximate since n**2>>k**2 (n**2-k**2 ~= n**2) for ice
      nsol=sqrt(eprsol)
      ksol=episol/(2*nsol)
      msol=cmplx(nsol,ksol)
      m2sol=msol**2

!     Determin ice refractive index and dielectric using Debye theory (1929)
!     which postulataes that K/rho is constant for an ice/air mixture
!     where rho is density and  K is (m**2-1)/(m**2+2) Battan (1973) p. 40.

      kovrrho=((m2sol-(1,0))/(m2sol+(2,0)))/(rhosol)

      m2=((2,0)*kovrrho*rho+(1,0))/((1,0)-kovrrho*rho)
      write(*,*)
      write(*,*) 'ice/air mixture:'
      write(*,*) 'rho = ',rho,'   m= ',sqrt(m2)
      write(*,*)
      write(*,*) 'refractive index:'
      n=real(sqrt(m2))
      write(*,*) 'n = ',n
      k=imag(sqrt(m2))
      write(*,*) 'k = ',k
      write(*,*)
      write(*,*) 'relative permittivity/dielectric:'
      epr=n**2-k**2
      write(*,*) 'er = ',epr
      epi=2*n*k
      write(*,*) 'ei = ',epi

      return
 end subroutine debye

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 

 end module module_tmatrix

!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
!TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX TMATRIX 
