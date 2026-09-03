!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_icescat_yang_visir
 use module_floating_point
 use module_mpi
 use netcdf
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
!  This module is designed to process Ping Yang UV-VIS-IR non-sherical ice single-scattering database. 
!  vis-ir optics module in the G-SDSU.
!
! History:
! 07/2015  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
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
  public :: prep_icescat_yang_visir , &  ! prepare (read and archive) uv-vis-ir single-scattering properties.
            interp_wave_icescat_yang_visir , &
            interp_size_icescat_yang_visir


! ==========================================================================
! ============= species registration =======================================
!       Add new species consitently between name_species and id_species
! ==========================================================================
 integer,parameter :: nSpecies = 8

 character(len=22),dimension(nSpecies), parameter :: name_species(nSpecies) = (/&! name of species
    'droxtals              ', &
    'solid_columns         ', &
    'plates                ', &
    'hollow_bullet_rosettes', &
    'solid_bullet_rosettes ', &
    'aggregate_columns     ', &
    'small_aggregate_plates', &
    'large_aggregate_plates'  &
    /)


! ==========================================================================
! ========== end of species registration ===================================
! ==========================================================================

!
! dimension bounds
!
 integer,parameter :: nWaveLen = 396 
 integer,parameter :: nRad = 189 
 integer,parameter,public :: nScatAng = 74 

!
! microphysics and single scattering parameters and anchor points from original database
! Note that this is already size-integrated via gamma PSD function. 
!
 real,allocatable :: phase_angles(:) !scattering angle [degrees]
 real,allocatable :: wavelengths(:)  !wavelength [microns]

 real,allocatable :: radius(:,:)       !concete spherical ice radius [microns]
 real,allocatable :: max_diam(:,:)     !maximum diameter [microns2]
 real,allocatable :: rhoe(:,:)         !effective density [g/cm3]
 real,allocatable :: extinction_efficiency   (:,:,:)  ![-]
 real,allocatable :: single_scattering_albedo(:,:,:)  ![-]
 real,allocatable :: asymmetry_parameter     (:,:,:) ![-]
 real,allocatable :: p11_phase_function(:,:,:,:) ![-]
 real,allocatable :: p21_phase_function(:,:,:,:) ![-]
 real,allocatable :: p22_phase_function(:,:,:,:) ![-]
 real,allocatable :: p33_phase_function(:,:,:,:) ![-]
 real,allocatable :: p43_phase_function(:,:,:,:) ![-]
 real,allocatable :: p44_phase_function(:,:,:,:) ![-]

!
! interpolarted single scattering parameters and anchor points from original database
!
 real,allocatable :: cext_icescat (:,:)  ! size-integrated extinction coeff per unit mass [cm2/g]
 real,allocatable :: cabst_icescat(:,:)  ! size-integrated absorption coeff per unit mass [cm2/g]
 real,allocatable :: ph_icescat   (:,:,:)! size-integrated phase function [-] 
 real,allocatable :: wavel_icescat(:)    ! wavelgnth
 integer,save     :: nwav_icescat        ! total number of wavelength
 integer,save     :: nang_rstar

!
! linear depolarization ratio (depol)
!
! depol = [ P11(180deg)-P22(180deg) ] / [ P11(180deg)+P22(180deg) ]
! For spherical particle  P11(180deg) = P22(180deg), so that depol = 0.0
!

!
! wave interpolated single scattering properties
! 
 real :: rad (1:nRad) ! concrete spherical radius [um]
 real :: qext(1:nRad) ! extinction efficiencey [-]
 real :: w0  (1:nRad) ! single scattering albeod (0-1) [-]
 real :: p11 (1:nRad,1:nScatAng) ! P11 phase function [-]
 real :: p22 (1:nRad,1:nScatAng) ! P22 phase function [-]


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine prep_icescat_yang_visir(dir, ss_opt )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine prepare to open Ping Yang Ice Scattering database, and stored in LUT as a function of frequency, temperature and 
! crystal type.
!
! History:
! 07/2015  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: dir  !input database directory
 integer,intent(in) :: ss_opt  !single scattering table options 
 integer :: ierr     !allocation stat

 character(len=300) :: io_file  !input file
 integer :: id


 if(masterproc) print*,'MSG prep_icescat_yang_visir: Prepare Ping Yang VIS-IR IceScat database '

!
! allocate single scattering parameters
!
 if( .not. allocated(phase_angles) ) then

 allocate( &
  phase_angles(nScatAng) ,&  ![degrees]
  wavelengths(nWaveLen)  ,& ![microns]
  radius(nRad,nSpecies) ,& ![microns]
  max_diam(nRad,nSpecies) ,& ![microns]
  rhoe(nRad,nSpecies) ,& ! 
  extinction_efficiency   (nWaveLen, nRad, nSpecies) ,&
  single_scattering_albedo(nWaveLen, nRad, nSpecies) ,&
  asymmetry_parameter     (nWaveLen, nRad, nSpecies) ,&
  p11_phase_function      (nWaveLen, nRad, nScatAng, nSpecies) ,&
  p22_phase_function      (nWaveLen, nRad, nScatAng, nSpecies) ,&  !needed for depolarization ratio
 stat=ierr )

 if(ierr /= 0) call stop_sdsu(myrank,'MSG prep_icescat_yang_visir: cannot allocate the variable')
 endif

 !
 ! read database
 !
 do id = 1, nSpecies

    io_file=trim(dir)//'icescat_yang_'//trim(name_species(id))//'_Rough050.nc'

    call read_icescat_visir( trim(io_file), id )

 enddo

 !
 ! interpolate and store in memory
 !
 ! call store_icescat_visir( mxwavel, wavel(1:mxwavel), nang, ang )

 
 return
 end subroutine prep_icescat_yang_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_icescat_visir(io_icescat,id_spec)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine open Ping Yang Ice Scattering database, and stored in LUT as a function of frequency, temperature and 
! crystal type.
!
! History:
! 07/2015  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: io_icescat   !io_file name
 integer,intent(in) :: id_spec

 integer :: i,j,k,n  !loop
 integer :: ncid     !# of variables to inquire  
 integer :: varid    !# of variable id
 logical :: missing_para
 logical :: fexist

!
! file name
!

 if(masterproc) print("(a)"), 'MSG read_icescat_visir: reading...'//trim(io_icescat)
 inquire(file=trim(io_icescat),exist=fexist)
 if( .not. fexist ) then 
!     print*,'MSG read_icescat_visir; Cannot Find Database ',trim(io_icescat)
     call stop_sdsu(myrank,'MSG read_icescat_visir; Cannot Find Database '//trim(io_icescat))
 endif

!
! open netCDF file
!
 call check( nf90_open(trim(io_icescat), nf90_nowrite, ncid) )

!
! phase function angles (nScatAng) [degree]
!
 call check( nf90_inq_varid(ncid,'phase_angles',varid ),'phase_angles', missing_para) 
 call check( nf90_get_var(ncid, varid, phase_angles ))

!
! wave length (nWaveLen) [micron]
!
 call check( nf90_inq_varid(ncid,'wavelengths',varid ),'wavelengths', missing_para) 
 call check( nf90_get_var(ncid, varid, wavelengths ))

!
! concrete spherical ice radius (nRad) [micron]
!
 call check( nf90_inq_varid(ncid,'radius',varid ),'radius', missing_para)
 call check( nf90_get_var(ncid, varid, radius(1:nRad,id_spec) ))

!
!  maximum diameter (nRad) [micron2]
!
 call check( nf90_inq_varid(ncid,'max_diam',varid ),'max_diam', missing_para)
 call check( nf90_get_var(ncid, varid, max_diam(1:nRad,id_spec) ))

!
!  effectie diameter (nRad) [micron2]
!
 call check( nf90_inq_varid(ncid,'rhoe',varid ),'rhoe', missing_para)
 call check( nf90_get_var(ncid, varid, rhoe(1:nRad,id_spec) ))


!
! extinction_efficiency (nWaveLen, nRad) [-] 
!
 call check( nf90_inq_varid(ncid,'extinction_efficiency',varid ),'extinction_efficiency', missing_para)
 call check( nf90_get_var(ncid, varid, extinction_efficiency(1:nWaveLen,1:nRad,id_spec) ))

!
! single_scattering_albedo (nWaveLen, nRad) [-] 
!
 call check( nf90_inq_varid(ncid,'single_scattering_albedo',varid ),'single_scattering_albedo', missing_para)
 call check( nf90_get_var(ncid, varid, single_scattering_albedo(1:nWaveLen,1:nRad,id_spec) ))

!
! asymmetry_parameter (nWaveLen, nRad) [-] 
!
 call check( nf90_inq_varid(ncid,'asymmetry_parameter',varid ),'asymmetry_parameter', missing_para)
 call check( nf90_get_var(ncid, varid, asymmetry_parameter(1:nWaveLen,1:nRad,id_spec) ))

!
! p11_phase_function (nWaveLen, nRad, nScatAng) [-] 
!
 call check( nf90_inq_varid(ncid,'p11_phase_function',varid ),'p11_phase_function', missing_para)
 call check( nf90_get_var(ncid, varid, p11_phase_function(1:nWaveLen,1:nRad,1:nScatAng,id_spec) ))

!
! p22_phase_function (nWaveLen, nRad, nScatAng) [-] 
!
 call check( nf90_inq_varid(ncid,'p22_phase_function',varid ),'p22_phase_function', missing_para)
 call check( nf90_get_var(ncid, varid, p22_phase_function(1:nWaveLen,1:nRad,1:nScatAng,id_spec) ))


 return
 end subroutine read_icescat_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine interp_wave_icescat_yang_visir(wave_in,species,success)
 implicit none
 real,intent(in) :: wave_in ! wavelength [um] 
 character(len=*),intent(in) :: species !
 logical,intent(inout) :: success !success criteria of interpolation

 integer :: n , nw
 integer :: id
 real :: wgt1_wave, wgt2_wave

 success = .false.

!
! find matched species name and id
!
 id = -999
 do n = 1, nSpecies
    if( trim(species) == trim( name_species(n) ) ) then
       id = n
       exit
    endif
 enddo

 if( id == -999 ) then
    print*, 'MSG interp_wave_icescat_yang_visir: Warning. species name does not match up with database',trim(species)
    return 
 endif

!
! check input wavelength bounds
!
 if( wave_in < wavelengths(1) ) then
   print*, 'MSG interp_wave_icescat_yang_visir: Warning. input wavelength is out of bounds....'
  print*, 'min_wave[um]=',wavelengths(1),'max_wave[um]=',wavelengths(nWaveLen),'wave_in',wave_in
   print*, 'proceed anyway'
   qext(1:nRad) =  extinction_efficiency(1  , 1:nRad, id) 
   w0(1:nRad)   =  single_scattering_albedo(1  , 1:nRad, id) 
   p11(1:nRad,1:nScatAng)  = p11_phase_function(1  ,1:nRad, 1:nScatAng, id) 
   p22(1:nRad,1:nScatAng)  = p22_phase_function(1  ,1:nRad, 1:nScatAng, id) 

 elseif( wave_in > wavelengths(nWaveLen) ) then
   print*, 'MSG interp_wave_icescat_yang_visir: Warning. input wavelength is out of bounds....'
   print*, 'min_wave[um]=',wavelengths(1),'max_wave[um]=',wavelengths(nWaveLen),'wave_in',wave_in
   print*, 'proceed anyway'
   qext(1:nRad) =  extinction_efficiency(nWaveLen  , 1:nRad, id)      
   w0(1:nRad)   =  single_scattering_albedo(nWaveLen  , 1:nRad, id)      
   p11(1:nRad,1:nScatAng)  = p11_phase_function(nWaveLen  ,1:nRad, 1:nScatAng, id)      
   p22(1:nRad,1:nScatAng)  = p22_phase_function(nWaveLen  ,1:nRad, 1:nScatAng, id)      

 else

 wave_loop:  do nw = 1, nWaveLen-1  

     wave_if: if( wave_in >= wavelengths(nw) .and. wave_in <= wavelengths(nw+1) ) then  ! interpolate

        ! weighting points of wavelength
        wgt1_wave = ( wavelengths(nw+1) - wave_in)  / ( wavelengths(nw+1) - wavelengths(nw) )
        wgt2_wave = 1. - wgt1_wave

        qext(1:nRad) = wgt1_wave * extinction_efficiency(nw  , 1:nRad, id) + &
                       wgt2_wave * extinction_efficiency(nw+1, 1:nRad, id)

        w0(1:nRad)   = wgt1_wave * single_scattering_albedo(nw  , 1:nRad, id) + &
                       wgt2_wave * single_scattering_albedo(nw+1, 1:nRad, id)

        p11(1:nRad,1:nScatAng)  = wgt1_wave * p11_phase_function(nw  ,1:nRad, 1:nScatAng, id) + &
                                  wgt2_wave * p11_phase_function(nw+1,1:nRad, 1:nScatAng, id)

        p22(1:nRad,1:nScatAng)  = wgt1_wave * p22_phase_function(nw  ,1:nRad, 1:nScatAng, id) + &
                                  wgt2_wave * p22_phase_function(nw+1,1:nRad, 1:nScatAng, id)

    endif wave_if

 enddo wave_loop

 endif

 rad(1:nRad) = radius(1:nRad, id)
 success = .true.


 return
 end subroutine interp_wave_icescat_yang_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine interp_size_icescat_yang_visir(rad_in, qext_interp, w0_interp, p11_interp, p22_interp )
 implicit none
 real,intent(in) :: rad_in ! radius [um] 
 real,intent(out) :: qext_interp  !interpolated extinction efficiency [-]
 real,intent(out) :: w0_interp    !interpolated single scattering albedo [-]
 real,intent(out) :: p11_interp(nScatAng) !P11 phase function [-]
 real,intent(out) :: p22_interp(nScatAng) !P22 phase function [-]

 integer :: nd
 real :: wgt1_rad, wgt2_rad

 if(rad_in < rad(1) ) then  !small size

   qext_interp = qext(1) 
   w0_interp   = w0(1)   
   p11_interp(1:nScatAng) = p11(1,1:nScatAng)  
   p22_interp(1:nScatAng) = p22(1,1:nScatAng)  

 elseif( rad_in > rad(nRad) ) then ! large size

   qext_interp = qext(nRad) 
   w0_interp   = w0(nRad)   
   p11_interp(1:nScatAng) = p11(nRad,1:nScatAng)       
   p22_interp(1:nScatAng) = p22(nRad,1:nScatAng) 

 else

    nrad_loop: do nd = 1, nRad - 1

      nrad_if: if( rad_in >= rad(nd) .and. &
               rad_in <= rad(nd+1) ) then

        !weighting points of Deff
        wgt1_rad = ( rad(nd+1) - rad_in     ) / &
                   ( rad(nd+1) - rad(nd) )
        wgt2_rad = 1.-wgt1_rad

        qext_interp = wgt1_rad * qext(nd) + wgt2_rad * qext(nd+1)
        w0_interp   = wgt1_rad * w0  (nd)  + wgt2_rad * w0  (nd+1)
        p11_interp(1:nScatAng) = wgt1_rad * p11(nd,1:nScatAng) + wgt2_rad * p11(nd+1,1:nScatAng)
        p22_interp(1:nScatAng) = wgt1_rad * p22(nd,1:nScatAng) + wgt2_rad * p22(nd+1,1:nScatAng)

      endif nrad_if

    enddo nrad_loop

 endif

 return
 end subroutine interp_size_icescat_yang_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine check(status,para,missing_para)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Check status of reading process of netCDF. 
!
! History:
! 05/2014  Toshi Matsui@NASA GSFC ; included optional statement for skiping
! 05/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!---------------------------------------------------------------------------------------------------
  integer, intent (in) :: status
  character(len=*),intent(in),optional :: para  !optional parameter for missing, but skippable parameters.
  logical,intent(out),optional :: missing_para  !optional parameter for missing, but skippable parameters.

  if( present(missing_para) )  missing_para = .false.

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    if( present(para) ) then
          if(masterproc) print*,'Missing ',trim(para)
          if( present(missing_para) ) missing_para = .true.
          return
    endif
    stop "Stopped missing parameters!"
  end if
 end subroutine check

 end module module_icescat_yang_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
