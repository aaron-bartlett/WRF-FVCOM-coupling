!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 MODULE module_tmatrix_database
 IMPLICIT NONE
!--------------------------------------------------------------------------------
! MODULE read_sca_para_snow_table.f90 --> module_tmatrix
! 
! Read single scattering parameters of randomly oriented oblate particle 
!  from T-matrilx calculation.
!
!  To access scattering data from snow table, use following SUBROUTINE:
!
!     SUBROUTINE read_table_snow(filename,qb,qs,qe,g,snow_density,diameter,&
!                                temperature,wavelength)
!
!        Input:
!           filename:         name of snow scattering look-up table 
!        Output:
!           wavelength:      (mm)
!           snow_density:    snow density (g/cm^3)
!           diameter:        melted particle diameter (mm)
!           aspect_ratio:    particle aspect ratio (oblate<1,sphere=1,prolate>1)
!           temperature:     temperature (C)
!           nRho:            number of mass densities 
!           nsizebin:        number of particle size
!           qb:              backscattering xsection (mm^2)
!           qs:              scattering xsection (mm^2)
!           qe:              extinction xsection (mm^2)
!           g:               asymmetry factor (dimensionless)
!
!
!
!  History:
!  01/2013  Toshi Matsui@NASA GSFC: Modified for various paramters setting and structure. 
!  03/2012  Liang Liao@NASA GSFC  : Init
!
! Reference
! Liao, L. and R. Meneghini, 2013: Examination of Effective Dielectric Constants of Nonspherical 
!      Mixed-Phase Hydrometeors. J. Appl. Meteor. Climatol., 52, 197-212. 
!      doi: http://dx.doi.org/10.1175/JAMC-D-11-0244.1 
!
!-------------------------------------------------------------------------------- 

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
  private   ! encapsulate all variables and subourtines (non accessible)

!
! public (accessible) subroutines
!
  public :: tmatrix_prep, & ! create oblate microwave single-scattering properties LUTs
            tmatrix_interp  ! interpolate single-scattering propertes for a givne sphere radius



!
! T-matrix LUT parameters
!

 INTEGER, PARAMETER :: nSizeBin=400, &  !number of melted diameters
                       nRho=10           !number of snow mass densities

 type single_scatter
   real :: aspect_ratio              !aspect ratio of ellipsoid
   real :: rhoe_ellip(1:nRho)        ! Snow Density [g.cm3], based on circumscribing ellipsoid
   real :: aeff(1:nSizeBin) ! solid-ice radius [micron]
   real :: wave       ! wavelength [um]
   real :: freq       ! frequency [GHz]
   real :: Q_ext(1:nRho,1:nSizeBin)      ! Extinction cross section [mm2]
   real :: Q_abs(1:nRho,1:nSizeBin)      ! Absorption cross section [mm2]
   real :: Q_sca(1:nRho,1:nSizeBin)      ! Scattering cross section [mm2]
   real :: g    (1:nRho,1:nSizeBin)      ! Asymetry parameters (1) = cos(theta) , (2) = cos2(theta)
   real :: Q_bk (1:nRho,1:nSizeBin)      ! backscattering cross section / (4*pi) [mm2]
 end type single_scatter
 type (single_scatter) :: sample

 integer,parameter :: Tanchor_point = 13
 character(len=*),parameter :: Tfreq_char(Tanchor_point) = &
    (/ '10.65G', '13.6G ', '18.7G ', '23.8G ', '35.5G ', '89G   ', &
       '94G   ', '150G  ', '165.5G', '176G  ', '180G  ', '186G  ', '190G  ' /)

 type (single_scatter) :: Tsample(Tanchor_point)


 REAL, DIMENSION(1:nSizeBin) :: diameter
 REAL, DIMENSION(1:nRho,1:nSizeBin) :: qb_snow,qs_snow,qe_snow,g_snow
 REAL, DIMENSION(1:nRho) :: snow_density
 REAL :: wavelength,t_snow,aspect_ratio

!
! other parameters
!
 character(len=200),dimension(:),allocatable :: inp_dir  !input directory
 integer :: nmax_dir
 real,parameter :: undefined = -999.e0
 character(len=3),parameter :: undefined_cha = 'NAN'
 integer,parameter :: undefined_int = -999
 integer,parameter  :: data_len = 7500

 real,parameter :: pi = 3.141592653589793
         
 CONTAINS

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine tmatrix_prep(i_can_say,dir)
 implicit none
 logical,intent(in) :: i_can_say
 character(len=*),intent(in) :: dir   !input data directory
 integer :: f
 character(len=5) :: aspect_ratio
 character(len=300) :: io_file 


 if(i_can_say) print*,'MSG tmatrix_prep, Preparing microwave T-matrix database'

!
! currently apect ratio is always 0.5 (will be changed future)
!
 aspect_ratio = '0.5'


!
! looping over frequency anchor points
!
 do f = 1, Tanchor_point 

  !
  ! input file and directory
  !
  io_file = trim(dir)//'Tmatrix/scaPara_spheroid_aspectRatio-'//trim(aspect_ratio)//'_snow_table_'//trim(Tfreq_char(f))

  !
  ! reading file and store in module parameters
  !
  call driver_read_file(f,io_file)

 enddo

 end subroutine tmatrix_prep

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine driver_read_file(f,io_file)
 IMPLICIT NONE
!--------------------------------------------------------------------------------
! PROGRAM mainProgram_read_sca_para_table.f90
! 
!  To access scattering results from snow table, use following SUBROUTINE:
!
!     SUBROUTINE read_table_rain(filename,qb,qs,qe,g,temperature,diameter,wavelength)
!
!        Input:
!           filename:       name of rain table 
!        Output:
!           wavelength:     (mm)
!           temperature:     rain temperature (C)
!           diameter:        particle diameter (mm)
!           nTemp:           number of temperatures
!           nsizebin:        number of particle size
!           qb:              backscattering xsection (mm^2)
!           qs:              scattering xsection (mm^2)
!           qe:              extinction xsection (mm^2)
!           g:               asymmetry factor (dimensionless)
!
!-------------------------------------------------------------------------------- 
 integer,intent(in) :: f ! frequency index
 character(len=*),intent(in) :: io_file
 REAL, DIMENSION(1:nSizeBin) :: diameter
 REAL, DIMENSION(1:nRho,1:nSizeBin) :: qb,qs,qe,g
 REAL, DIMENSION(1:nRho) :: snow_density
 REAL :: wavelength,t_snow,aspect_ratio
 INTEGER :: i,j,n,m
 real :: pia2
 real :: r

!
! -- Read T-matrix Single Scattering table for Snow
!   
 call read_table_snow(io_file,qb,qs,qe,g,snow_density, &
                      diameter,t_snow,aspect_ratio,wavelength)

!
! convert parameters to the global parameters
!
 Tsample(f)%wave = wavelength * 1.e3  ![micron]
 Tsample(f)%aspect_ratio = aspect_ratio  !aspect ratio [-]
 Tsample(f)%freq = 2.997925e2/wavelength

!
! Note that original diameter is melted (water) shere diameter 
! So convert into solid-ice diameter
!
 Tsample(f)%aeff(1:nSizeBin) = (1.0/0.917)**(1./3.) * 0.5* diameter(1:nSizeBin) * 1.e3 !solid-ice radius [micron]

 Tsample(f)%rhoe_ellip(1:nRho) = snow_density(1:nRho) ! Snow Density [g.cm3] on ellipsoid

!Density ahncor points are 
! 0.05    0.1000000      0.2000000      0.3000000      0.4000000    
!  0.5000000      0.6000000      0.7000000      0.8000000      0.9160000 
!

 do m = 1,nSizeBin ; do n = 1,nRho
!    r = (1.0/0.917)**(1./3.) * 0.5* diameter(m)   !solid-ice sphere radius [mm]
!    pia2 = pi * r * r
! now stay as cross section [mm2]
    Tsample(f)%Q_ext(n,m) = qe(n,m)             
    Tsample(f)%Q_abs(n,m) = (qe(n,m)-qs(n,m)) 
    Tsample(f)%Q_sca(n,m) = qs(n,m)           
    Tsample(f)%g    (n,m) = g(n,m)
    Tsample(f)%Q_bk (n,m) = qb(n,m)           / (4.*pi)  ! per steradian
 enddo ; enddo


 return
 end subroutine driver_read_file

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

!-----------------------------------------------------------------------------------
! SUBROUTINE read_table_snow()
! 
! Read single scattering parameters from look-up tables of spheroidal snow particles
!
! Names of variables:
!   Input:
!      filename:   name of table
!   Output:
!      wavelength:      (mm)
!      snow_density:    snow density (g/cm^3)
!      diameter:        melted particle diameter (mm)
!      aspect_ratio:    particle aspect ratio (oblate<1,sphere=1,prolate>1)
!      temperature:     temperature (C)
!      nRho:            number of mass densities       
!      nsizebin:        number of particle size
! 
!      qb:              backscattering xsection (mm^2)
!      qs:              scattering xsection (mm^2)
!      qe:              extinction xsection (mm^2)
!      g:               asymmetry factor (dimensionless)
! 
!-----------------------------------------------------------------------------------

   SUBROUTINE read_table_snow(filename,qb,qs,qe,g,snow_density,diameter,temperature, &
                              aspect_ratio,wavelength)
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: filename
      REAL, INTENT(OUT), DIMENSION(1:nRho,1:nSizeBin) :: qb,qs,qe,g
      REAL, INTENT(OUT), DIMENSION(1:nRho) :: snow_density
      REAL, INTENT(OUT), DIMENSION(1:nSizeBin) :: diameter
      REAL, INTENT(OUT) :: wavelength,temperature,aspect_ratio
      CHARACTER(LEN=60) :: v
      INTEGER :: n, i, nRho_table, nSizeBin_table
      
      OPEN(UNIT=10, FILE=filename, STATUS='OLD')
      
      READ(10,'(1x,a17,f8.4)') v,wavelength         !wavelength in mm
      READ(10,'(1x,a26, i3)') v,nRho_table         !number of temperatures
      IF (nRho_table /= nRho) THEN
         PRINT*, 'Number of snow densities are unmatched with snow-table'
	 STOP
      ENDIF
      READ(10,'(1x,a26, i4)') v,nSizeBin_table      !number of size bins from table
      IF (nSizeBin_table /= nSizeBin) THEN
         PRINT*, 'Number of particle sizes are unmatched with rain-table'
	 PRINT*, 'nSizeBin, nSizeBin_table =', nSizeBin, nSizeBin_table
	 STOP
      ENDIF
      READ(10, '(1x,a17,f5.3)') v,temperature
      READ(10, '(1x,a39,f4.2)'),v,aspect_ratio
      READ(10, '(1x,a22,10(f5.3,3x))'),v,snow_density
      READ(10, '(1x,a14,400(f6.3,3x),/,/,/)') v,diameter
      
      DO n=1, nRho
         READ(10,'(1x,400(e17.10,2x))') (qb(n,i),i=1,nSizeBin)
      ENDDO 
             
      READ(10,'(1x,a2,/,/)') v
      DO n=1, nRho
         READ(10,'(1x,400(e17.10,2x))') (qs(n,i),i=1,nSizeBin)
      ENDDO

      READ(10,'(1x,a2,/,/)') v
      DO n=1, nRho
         READ(10,'(1x,400(e17.10,2x))') (qe(n,i),i=1,nSizeBin)
      ENDDO

      READ(10,'(1x,a2,/,/)') v
      DO n=1, nRho
         READ(10,'(1x,400(e17.10,2x))') (g(n,i),i=1,nSizeBin)
      ENDDO
      
      CLOSE(UNIT=10)

   END SUBROUTINE read_table_snow

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_table(dir)
 IMPLICIT NONE
!--------------------------------------------------------------------------------
! PROGRAM mainProgram_read_sca_para_table.f90
! 
!  To access scattering results from snow table, use following SUBROUTINE:
!
!     SUBROUTINE read_table_rain(filename,qb,qs,qe,g,temperature,diameter,wavelength)
!
!        Input:
!           filename:       name of rain table 
!        Output:
!           wavelength:     (mm)
!           temperature:     rain temperature (C)
!           diameter:        particle diameter (mm)
!           nTemp:           number of temperatures
!           nsizebin:        number of particle size
!           qb:              backscattering xsection (mm^2)
!           qs:              scattering xsection (mm^2)
!           qe:              extinction xsection (mm^2)
!           g:               asymmetry factor (dimensionless)
!
!-------------------------------------------------------------------------------- 
 character(len=*),intent(in) :: dir

 integer :: i,j,n,m, f
 character(len=100) :: outfile
 integer,parameter :: io = 400

!
! input file directory and file name. 
!
 outfile = './../OUTPUT/'//trim(dir)//'.txt'

 print*,'MSG write_table: output file -->',trim(outfile)
 open(io,file=trim(outfile),status='replace')
 write(io,*)&
 '          aeff    rhoe(ellip)         Q_abs          Q_sca              g           Q_bk'

 do n = 1,nRho ; do m = 1,nSizeBin 

 write(io,100) sample%aeff(m), sample%rhoe_ellip(n), &
            sample%Q_abs(n,m), sample%Q_sca(n,m), sample%g(n,m), sample%Q_bk(n,m)
100 format(6E15.6)
 enddo ; enddo

 close(io)

 end subroutine write_table

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine write_table_grads(dir)
 implicit none
!--------------------------------------------------------------------------------
! PROGRAM mainProgram_read_sca_para_table.f90
! 
!  To access scattering results from snow table, use following SUBROUTINE:
!
!     SUBROUTINE read_table_rain(filename,qb,qs,qe,g,temperature,diameter,wavelength)
!
!        Input:
!           filename:       name of rain table 
!        Output:
!           wavelength:     (mm)
!           temperature:     rain temperature (C)
!           diameter:        particle diameter (mm)
!           nTemp:           number of temperatures
!           nsizebin:        number of particle size
!           qb:              backscattering xsection (mm^2)
!           qs:              scattering xsection (mm^2)
!           qe:              extinction xsection (mm^2)
!           g:               asymmetry factor (dimensionless)
!
!-------------------------------------------------------------------------------- 
 character(len=*),intent(in) :: dir
 integer :: i,j,n,m,f,irec
 character(len=100) :: outfile
 integer,parameter  :: io = 500
 character(len=3)   :: rrr  !density character

!
! input file directory and file name. 
!
 outfile = './../OUTPUT/'//trim(dir)//'.bin'

 print*,'MSG write_table_grads: output file -->',trim(outfile)
 open(io,file=trim(outfile),access='direct',status='replace',recl=data_len*4)

 irec = 1
 call dump( sample%aeff( 1:nSizeBin ) , io, irec ) 

 do n = 1,nRho
    call dump( sample%Q_abs(n, 1:nSizeBin ) , io, irec )
    call dump( sample%Q_sca(n, 1:nSizeBin ) , io, irec ) 
    call dump( sample%g    (n, 1:nSizeBin ) , io, irec ) 
    call dump( sample%Q_bk (n, 1:nSizeBin ) , io, irec ) 
 enddo 

 close(io)

!
! write grads control file
!
 outfile = './../OUTPUT/'//trim(dir)//'.ctl'

 print*,'MSG write_table_grads: GrADS control file -->',trim(outfile)
 open(io,file=trim(outfile),status='replace')

 outfile = trim(dir)//'.bin'

 write(io,*) 'dset  '//trim(outfile)
 write(io,*) 'title  '//trim(dir)
 write(io,*) 'undef  ',undefined
 write(io,*) 'xdef   ',data_len,' linear 1. 1.'
 write(io,*) 'ydef   1   linear   1.   1. '
 write(io,*) 'zdef   1   linear   1.   1. '
 write(io,*) 'tdef    1 linear 01:00Z06aug2006 1hr'
 write(io,*) 'vars   ', nRho*4+1
 write(io,*) 'aeff           0  0  Solid-sphere radius [micron]'

 do n = 1,nRho
    write(rrr,"(I3.3)") INT(sample%rhoe_ellip(n)*100.)
    write(io,*) 'Qabs'//rrr,'          0  0  Absorption efficiency [-]'
    write(io,*) 'Qsca'//rrr,'          0  0  Scattering efficiency [-]'
    write(io,*) '   g'//rrr,'          0  0  Asymetry parameter    [-]'
    write(io,*) '  Qb'//rrr,'          0  0  Backscatter efficiency[-]'
 enddo

 write(io,"(a)") 'endvars '
 close(io)

 end subroutine write_table_grads

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dump( var1d, io, irec)
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

 real,dimension(:), intent(in)  :: var1d  ! assumed length
 integer,intent(in)    :: io   ! io index
 integer,intent(inout) :: irec ! binary record #
 integer :: i
 real :: temp(1:data_len)
 integer :: bnd(1)        ! bound of array
 
 bnd = UBOUND(var1d)

 !
 ! init
 !
 temp = undefined

 !
 ! transfer
 !
 do i = 1, bnd(1)
    temp(i) = var1d(i)
 enddo

 write(io,rec=irec) (temp(i),i=1,data_len)
 irec=irec+1

 end subroutine dump

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine tmatrix_interp(i_can_say, proc, freq, rhoe, rad_solid, area, &
                           qext, qsca, asym, qbsca, success )
 implicit none
 logical,intent(in) :: i_can_say
 character(len=*),intent(in) :: proc   
 real,intent(in),optional :: freq  !frequency [GHz]
 real,intent(in),optional :: rhoe  !effective density [g/cm3]
 real,intent(in),optional :: rad_solid   ! input solid-sphere ice radius [micron] (for LUT interpolation)
 real,intent(in),optional :: area     ! area of fluffy particle [mm2] 

 real,intent(out), optional :: qext  ! total extincion efficiency [-]
 real,intent(out), optional :: qsca  ! scattering efficiency [-]
 real,intent(out), optional :: asym  ! asymetry parameter [-]
 real,intent(out), optional :: qbsca ! backscattering efficiency / (4*pi) [-]
 logical,intent(out),optional :: success ! status of interpolation


 real :: qext_x  ! total extincion cross section [mm2] 
 real :: qsca_x  ! scattering cross section [mm2] 
 real :: asym_x  ! asymetry parameter [-]
 real :: qbsca_x ! backscattering cross section  / (4*pi) [mm2]
 integer :: f, ns, nr
 real :: wgt1, wgt2, dgt1, dgt2  !weiting function 
 

 select case(trim(proc))

 case ('prep')  ! ---------------------  interpolate for frequency -------------------------------


   !
   ! check input parameters
   !
   if( .not. present(freq) ) &
    stop 'MSG tmatrix_interp: rhoe or rs not present in interp case'

   if(freq < 10.e0 .or. freq > 190.) then
    if(i_can_say) print*,'MSG tmatrix_interp: Warning --> frequency is out of range (10 ~ 190GHz)',freq 
   endif

   !
   ! interpolate LUT for a given frequency
   !
   if( freq <= Tsample(1)%freq  ) then

     sample%rhoe_ellip = Tsample(1)%rhoe_ellip
     sample%aeff   = Tsample(1)%aeff
     sample%Q_ext  = Tsample(1)%Q_ext
     sample%Q_abs  = Tsample(1)%Q_abs
     sample%Q_sca  = Tsample(1)%Q_sca
     sample%g      = Tsample(1)%g    
     sample%Q_bk   = Tsample(1)%Q_bk 

   elseif(  freq >= Tsample(Tanchor_point)%freq ) then

     sample%rhoe_ellip = Tsample(Tanchor_point)%rhoe_ellip
     sample%aeff   = Tsample(Tanchor_point)%aeff
     sample%Q_ext  = Tsample(Tanchor_point)%Q_ext
     sample%Q_abs  = Tsample(Tanchor_point)%Q_abs
     sample%Q_sca  = Tsample(Tanchor_point)%Q_sca
     sample%g      = Tsample(Tanchor_point)%g   
     sample%Q_bk   = Tsample(Tanchor_point)%Q_bk     
 
   else  ! interpolation

     wgt1 = undefined ; wgt2 = undefined
     freq_loop: do f = 1, Tanchor_point-1
        if( freq >= Tsample(f)%freq .and. freq <= Tsample(f+1)%freq ) then
            wgt1 = (Tsample(f+1)%freq - freq  ) / ( Tsample(f+1)%freq -  Tsample(f)%freq )
            wgt2 = 1.e0 - wgt1

            sample%rhoe_ellip = Tsample(f)%rhoe_ellip  !no interppolation for rhoe
            sample%aeff       = Tsample(f)%aeff        !no interpolation for aeff
            sample%Q_ext  = wgt1*Tsample(f)%Q_ext + wgt2*Tsample(f+1)%Q_ext  ! (1:nRho,1:nSizeBin)
            sample%Q_abs  = wgt1*Tsample(f)%Q_abs + wgt2*Tsample(f+1)%Q_abs
            sample%Q_sca  = wgt1*Tsample(f)%Q_sca + wgt2*Tsample(f+1)%Q_sca
            sample%g      = wgt1*Tsample(f)%g     + wgt2*Tsample(f+1)%g
            sample%Q_bk   = wgt1*Tsample(f)%Q_bk  + wgt2*Tsample(f+1)%Q_bk
            exit freq_loop
        endif
     enddo freq_loop

     if( wgt1 == undefined ) then
         stop 'MSG tmatrix_interp: cannnot find weight function '
     endif


   endif


 case ('interp') ! ----------------------interpolate for density and radius---------------------------

  !
  ! intialize output
  !
  qext=undefined; qsca=undefined; asym=undefined; qbsca=undefined; success=.false.

  if( .not. present(rhoe) .or. .not. present(rad_solid) ) &
     stop 'MSG tmatrix_interp: rhoe or rs not present in interp case'

   !
   ! return if bad input ( success --> .false.)
   !
   if( rad_solid <  sample%aeff(1) .or. rad_solid > sample%aeff(nSizeBin) )  return
   if( rhoe < sample%rhoe_ellip(1) .or. rhoe > sample%rhoe_ellip(nRho) ) return
    
   !
   ! interpolate from freq-interpolated LUT
   ! 
   wgt1 = undefined ; wgt2 = undefined ; dgt1 = undefined ; dgt2 = undefined
   ns_loop: do ns = 1, nSizeBin-1
      if( rad_solid >= sample%aeff(ns) .and. rad_solid <= sample%aeff(ns+1) ) then
          wgt1 = ( sample%aeff(ns+1) - rad_solid ) / ( sample%aeff(ns+1) - sample%aeff(ns) ) 
          wgt2 = 1.e0 - wgt1
          
          do nr = 1, nRho-1
             if( rhoe >= sample%rhoe_ellip(nr) .and. rhoe <= sample%rhoe_ellip(nr+1) ) then
                dgt1 = ( sample%rhoe_ellip(nr+1) - rhoe ) / ( sample%rhoe_ellip(nr+1) - sample%rhoe_ellip(nr) )
                dgt2 = 1.e0 - dgt1 

                qext_x =  dgt1 * ( wgt1*sample%Q_ext(nr  ,ns  )  + wgt2*sample%Q_ext(nr  ,ns+1) )    &
                        + dgt2 * ( wgt1*sample%Q_ext(nr+1,ns  )  + wgt2*sample%Q_ext(nr+1,ns+1) )

                qsca_x =  dgt1 * ( wgt1*sample%Q_sca(nr  ,ns  )  + wgt2*sample%Q_sca(nr  ,ns+1) )    &
                        + dgt2 * ( wgt1*sample%Q_sca(nr+1,ns  )  + wgt2*sample%Q_sca(nr+1,ns+1) )

                asym_x =  dgt1 * ( wgt1*sample%g    (nr  ,ns  )  + wgt2*sample%g    (nr  ,ns+1) )    &
                        + dgt2 * ( wgt1*sample%g    (nr+1,ns  )  + wgt2*sample%g    (nr+1,ns+1) )

                qbsca_x =  dgt1 * ( wgt1*sample%Q_bk (nr  ,ns  )  + wgt2*sample%Q_bk (nr  ,ns+1) )    &
                         + dgt2 * ( wgt1*sample%Q_bk (nr+1,ns  )  + wgt2*sample%Q_bk (nr+1,ns+1) )

                !
                ! convert X-section to efficiency for a given fluffy radius in mm
                !
                qext =  qext_x / area  !exctinction efficiency [-]
                qsca =  qsca_x / area  !scattering efficiency [-]
                asym =  asym_x         !asymetry paraemter [-]  (asym is not crossection. So no change.)
                qbsca=  qbsca_x / area *4.*pi !backscattering efficiency [-]

                ! interpolation results
                success=.true.

                exit ns_loop 
             endif
          enddo
      endif
   enddo ns_loop



 case default
  stop 'MSG tmatrix_interp: no such proc case'
 end select 


 return
 end subroutine tmatrix_interp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 END MODULE module_tmatrix_database

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

