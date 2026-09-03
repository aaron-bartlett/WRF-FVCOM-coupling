!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_openssp_ext
 use module_const
!---------------------------------------------------------------------------------------------------
! Comments:  
!   This module read and interpolate OpenSSP-EXT single-scattering database. 
!   OpenSSP-EXT data is composed of 
!     1)  averages of Kuo's DDSCAT-computed single scattering from snow aggregate ensembles for smaller particles, and
!     2) extrapolated value of single scattering from DDSCAT-fitted Liao's T-matrix table for larger particles.
!
! History:
! 5/2020  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
! K.-S. Kuo, W.S. Olson, B.T. Johnson, M. Grecu, L. Tian, T.L. Clune, B.H. van Aartsen, A.J. Heymsfield, L. Liao, R. Meneghini.
!     The Microwave Radiative Properties of Falling Snow Derived from Nonspherical Ice Particle Models. Part I: An Extensive 
!     Database of Simulated Pristine Crystals and Aggregate Particles, and Their Scattering Properties. Journal of Applied 
!     Meteorology and Climatology, March 2016, Vol. 55, No. 3
!
! Liao, L. and R. Meneghini, 2013: Examination of Effective Dielectric Constants of Nonspherical 
!     Mixed-Phase Hydrometeors. J. Appl. Meteor. Climatol., 52, 197-212. 
!     doi: http://dx.doi.org/10.1175/JAMC-D-11-0244.1 
!
!----------------------------------------------------------------------------------------------------
 implicit none 

 private   !all subroutines and parameters are encapsulated.

 public :: read_openssp_ext, interp_openssp_ext

 save

!
! currently OpenSSP has 13 frequency anchor point
!
 integer,parameter :: anchor_point = 13
 character(len=*),parameter :: freq_char(anchor_point) = &
    (/ '011GHz', '014GHz', '019GHz', '024GHz', '036GHz', '089GHz',&
       '094GHz', '150GHz', '166GHz', '176GHz', '180GHz', '186GHz', '190GHz' /)
 real,parameter :: freq_real(anchor_point) = &
    (/ 10.65e0, 13.6e0, 18.7e0, 23.8e0, 35.5e0, 89.e0, &
       94.e0  , 150.e0,165.5e0,176.e0 ,180.e0, 186.e0, 190.e0 /)


! real,parameter :: size_max_ext = 5100.e0 ! maximum bin size[um]
! real,parameter :: dsize    =  50.e0 ! bin width [um]
 integer,parameter  :: nbin_max_ext = 102  


 real :: size_bin_ext(nbin_max_ext)      ! bin size [um] (Solid-sphere radius)
 real :: rhoe_ddscat_ext(nbin_max_ext)   ! ellipoidal (true) rhoe [g/cm3] 

 integer,parameter :: iQ_abs=1, iQ_sca=2, ig=3, iQ_bk=4
 real :: lut(nbin_max_ext,anchor_point,4)  !LUT of OpenSSP scattering (last 4 means abs, sca, g, bk)
 real :: flut(nbin_max_ext,4)              !frequenly-interpolated LUT


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_openssp_ext(masterproc, dir)
!---------------------------------------------------------------------------------------------------
! Comments:  
!  Read OpenSSP-EXT database in binary format.
!
! History:
! 5/2020  Toshi Matsui@NASA GSFC ; Initial
!           
! References: 
!----------------------------------------------------------------------------------------------------
 logical,intent(in) :: masterproc    !if true, print statement
 character(len=*),intent(in) :: dir  !input directory

 integer :: f, irec
 integer,parameter :: io = 981
 character(len=300) :: outfile

!
! dump output for G-SDSU 
!

 do f = 1, anchor_point

  !
  ! toshii- uncertainties of openssp can be diagnosed by choosing below choice
  !
  outfile = trim(dir)//'AB_ddscat_ext_'//trim(freq_char(f))//'.dat'  !average value
!  outfile = trim(dir)//'HB_ddscat_ext_'//trim(freq_char(f))//'.dat'  !upper bound
!  outfile = trim(dir)//'LB_ddscat_ext_'//trim(freq_char(f))//'.dat'  !lower bound

  open(io,file=trim(outfile),access='direct',status='old',recl=nbin_max_ext*4)
  if(masterproc) print*,'MSG read_openssp_ext: reading OpenSSP-EXT datafiles --> ',trim(outfile)

!
! read binary file
!
  irec = 1
  call read_bin( size_bin_ext    (1:nbin_max_ext)     , io, irec )
  call read_bin( rhoe_ddscat_ext (1:nbin_max_ext)     , io, irec )
  call read_bin(   lut(1:nbin_max_ext,f,iQ_abs)       , io, irec )
  call read_bin(   lut(1:nbin_max_ext,f,iQ_sca)       , io, irec )
  call read_bin(   lut(1:nbin_max_ext,f,ig    )       , io, irec )
  call read_bin(   lut(1:nbin_max_ext,f,iQ_bk )       , io, irec )

!print*, lut(1:nbin_max_ext,f,iQ_sca)

  close(io)

 enddo

 return
 end subroutine read_openssp_ext

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_bin( var1d, io, irec)
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
 real,dimension(:), intent(out)  :: var1d  ! assumed length
 integer,intent(in)    :: io   ! io index
 integer,intent(inout) :: irec ! binary record #
 integer :: i
 integer :: bnd(1)        ! bound of array

 bnd = UBOUND(var1d)
 read(io,rec=irec) ( var1d(i),i=1,bnd(1) )
 irec=irec+1
 end subroutine read_bin

!DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT 
!DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT DDSCAT 

 subroutine interp_openssp_ext(freq, area, d_mm, qext, qsca, asym, qbsca, success)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine interplote prepared OpenSSP LUTs for a given frequency, solid-ice sphere diameter.
!
! History:
! 08/2020  Toshi Matsui@NASA GSFC ; Add conversion of scattering efficiency for fluffy X-section. 
! 06/2020  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 real,intent(in)  :: freq    ! frequency [GHz]
 real,intent(in)  :: area    ! cross-section of fully particle [m2]
 real,intent(in)  :: d_mm    ! solid ice sphere diameter [mm]
 real,intent(out) :: qext    ! total extincion efficiency [-]
 real,intent(out) :: qsca    ! scattering efficiency [-]
 real,intent(out) :: asym    ! asymetry parameter [-]
 real,intent(out) :: qbsca   ! backscattering efficiency [-] / 4pi radiance
 logical,intent(out) :: success ! status of interpolation

 integer :: f,n   !loop indice
 real :: wgt1, wgt2  !interpolatin weight
 real :: Q_abs    ! Absorption efficiency = C_abs / (pi* aeff^2)
 real :: Q_sca    ! Scattering efficiency = C_sca / (pi* aeff^2)
 real :: g        ! Asymetry parameters (1) = cos(theta) , (2) = cos2(theta)
 real :: Q_bk     ! backscattering efficiency ( C_back (area per str) / (pi* aeff^2) )
 real :: r_micron ! radius in micron [micron]
 real :: area_solid !cross section of solid sphere [m2]

 r_micron = 0.5* d_mm * 1.e3

 area_solid = const_pi * (d_mm*0.5e-3)*(d_mm*0.5e-3) !cross section of solid sphere [m2]

!
! initialize output
!
 qext = undefined ; qsca = undefined ; qbsca= undefined ; asym = undefined ; success = .false.

!
!if negative or too large particle, failed -> return
!
 if( r_micron <= 0. .or. r_micron > size_bin_ext(nbin_max_ext-1) ) return


!
! frequency interpolation first
!
 flut = undefined
 if( freq <= freq_real(1) ) then
     flut(1:nbin_max_ext,1:4) = lut(1:nbin_max_ext,1,1:4) 
 elseif( freq > freq_real(anchor_point) ) then
     flut(1:nbin_max_ext,1:4) = lut(1:nbin_max_ext,anchor_point,1:4)
 else
     freq_loop: do f = 1, anchor_point-1    
        if( freq_real(f) < freq .and. freq <= freq_real(f+1) ) then
            wgt1 = ( freq_real(f+1) - freq ) / ( freq_real(f+1) - freq_real(f) )
            wgt2 = 1.-wgt1
            if(wgt2<0. .or. wgt2>1.) stop 'MSG interp_openssp_ext: strange wgt2 in freq_loop'
            flut(1:nbin_max_ext,1:4) = wgt1 * lut(1:nbin_max_ext,f,1:4) + wgt2 * lut(1:nbin_max_ext,f+1,1:4)
            exit freq_loop
        endif
     enddo freq_loop

 endif

!
! particle size interpolation
!
 if( r_micron <= size_bin_ext(2) ) then  !avoid 1st bin due to undefined value

   wgt1 = ( size_bin_ext(2) - r_micron ) / size_bin_ext(2)
   wgt2 = 1.-wgt1
   if(wgt2<0. .or. wgt2>1.) stop 'MSG interp_openssp_ext: strange wgt2 in size '

   Q_abs =  wgt2* flut(2,iQ_abs)
   Q_sca =  wgt2* flut(2,iQ_sca)
   g     =  wgt2* flut(2,ig    )
   Q_bk  =  wgt2* flut(2,iQ_bk )
   success = .true.

 else

     size_loop: do n = 2, nbin_max_ext-2  !avoid 1st and final bin due to undefined value
        if( size_bin_ext(n) < r_micron .and. r_micron <= size_bin_ext(n+1) ) then
            wgt1 = ( size_bin_ext(n+1) - r_micron ) / (  size_bin_ext(n+1) -  size_bin_ext(n) )
            wgt2 = 1.-wgt1
            if(wgt2<0. .or. wgt2>1.) stop 'MSG interp_openssp_ext: strange wgt2 in size_loop'
            Q_abs =  wgt1* flut(n,iQ_abs) + wgt2* flut(n+1,iQ_abs)
            Q_sca =  wgt1* flut(n,iQ_sca) + wgt2* flut(n+1,iQ_sca)
            g     =  wgt1* flut(n,ig    ) + wgt2* flut(n+1,ig    )
            Q_bk  =  wgt1* flut(n,iQ_bk ) + wgt2* flut(n+1,iQ_bk )
            success = .true.
            exit size_loop
        endif
     enddo size_loop
 
 endif

!
! convert scattering efficienty calculated by solid sphere to scattering efficienty by actual fluffy sphere
!
 if( success ) then
      qext = (Q_abs + Q_sca)*area_solid / area
      qsca = Q_sca *area_solid / area
      qbsca= Q_bk  *area_solid / area * 4.*const_pi
      asym = g 
 endif

 return
 end subroutine interp_openssp_ext

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_openssp_ext

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
