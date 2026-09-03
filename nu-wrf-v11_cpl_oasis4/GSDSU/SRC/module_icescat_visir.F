!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_icescat_visir
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
!  This module is designed to process Bryan Baum UV-VIS-IR non-sherical ice single-scattering database. 
!  vis-ir optics module in the G-SDSU.
!
! History:
! 07/2015  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
! Baum, B.A., P. Yang, A.J. Heymsfield, A.R. Bansemer, B.H. Cole, A. Merrelli, C.G. Schmitt, 
!    and C. Wang, 2014: Ice cloud single-scattering property models with the full phase matrix 
!    at wavelengths from 0.2 to 100 µm. Journal of Quantitative Spectroscopy & Radiative Transfer
!    , 146, 123-139, DOI: 10.1016/j.jqsrt.2014.02.029.
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
  public :: prep_icescat_visir, &  ! prepare (read and archive) uv-vis-ir single-scattering properties.
            interp_icescat_visir   ! interpolate single scattering tables

!
! dimension bounds
!
 integer :: nDeff = 23 ;
 integer :: nScatAng = 498 ;
 integer :: nWaveLen = 445 ;

!
! microphysics and single scattering parameters and anchor points from original database
! Note that this is already size-integrated via gamma PSD function. 
!
 real,allocatable :: phase_angles(:)  ![degrees]
 real,allocatable :: wavelengths(:)  ![microns]
 real,allocatable :: effective_diameter(:) ![microns]
 real,allocatable :: ice_water_content(:,:)  ! [gm-3]
 real,allocatable :: total_area(:, :)  ! [microns2]
 real,allocatable :: single_scattering_albedo(:, :) 
 real,allocatable :: extinction_efficiency(: , :) 
 real,allocatable :: extinction_coefficient_over_iwc(:, :) ; ![ m2/g ]
 real,allocatable :: asymmetry_parameter(:, :) 
 real,allocatable :: p11_phase_function(:,:,:) 
 real,allocatable :: p21_phase_function(:,:,:) 
 real,allocatable :: p22_phase_function(:,:,:)
 real,allocatable :: p33_phase_function(:,:,:) 
 real,allocatable :: p43_phase_function(:,:,:) 
 real,allocatable :: p44_phase_function(:,:,:) 

!
! interpolarted single scattering parameters and anchor points from original database
!
 real,allocatable :: cext_icescat (:,:)  ! size-integrated extinction coeff per unit mass [cm2/g]
 real,allocatable :: cabst_icescat(:,:)  ! size-integrated absorption coeff per unit mass [cm2/g]
 real,allocatable :: ph_icescat   (:,:,:)! size-integrated phase function [-] 
 real,allocatable :: wavel_icescat(:)    ! wavelgnth
 integer,save          :: nwav_icescat        ! total number of wavelength
 integer,save          :: nang_rstar



!
! linear depolarization ratio (depol)
!
! depol = [ P11(180deg)-P22(180deg) ] / [ P11(180deg)+P22(180deg) ]
! For spherical particle  P11(180deg) = P22(180deg), so that depol = 0.0
!

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine prep_icescat_visir(dir,ss_opt, mxwavel, wavel, nang, ang )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine prepare to open Bryan Baum Ice Scattering database, and stored in LUT as a function of frequency, temperature and 
! crystal type.
!
! History:
! 07/2015  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!
! Baum, B.A., P. Yang, A.J. Heymsfield, A.R. Bansemer, B.H. Cole, A. Merrelli, C.G. Schmitt, 
!    and C. Wang, 2014: Ice cloud single-scattering property models with the full phase matrix 
!    at wavelengths from 0.2 to 100 µm. Journal of Quantitative Spectroscopy & Radiative Transfer
!    , 146, 123-139, DOI: 10.1016/j.jqsrt.2014.02.029.
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: dir  !input database directory
 integer,intent(in) :: ss_opt  !single scattering table options 
 integer,intent(in) :: mxwavel !total channel of visir simulator
 real,intent(in)    :: wavel(mxwavel) ! channel wavelength [micron]
 integer,intent(in) ::nang
 real,intent(in)    :: ang(1:nang)

 character(len=300) :: io_file  !input file


 if(masterproc) print*,'MSG prep_icescat_visir: Prepare Bryan Baum VIS-IR IceScat database '

 !
 ! read database
 !
 select case(ss_opt)
 case(1)  

   io_file=trim(dir)//'AggregateSolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc'
   call read_icescat_visir( trim(io_file) )

 case(2)

   io_file=trim(dir)//'GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc'
   call read_icescat_visir( trim(io_file) )

 case(3)

   io_file=trim(dir)//'SolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc'
   call read_icescat_visir( trim(io_file) )

 case default

   stop 'strange options. it must be 0< ss_opt_visir  < 4'

 end select


 !
 ! interpolate and store in memory
 !
 call store_icescat_visir( mxwavel, wavel(1:mxwavel), nang, ang )

 
 return
 end subroutine prep_icescat_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_icescat_visir(io_icescat)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine open Bryan Baum Ice Scattering database, and stored in LUT as a function of frequency, temperature and 
! crystal type.
!
! History:
! 07/2015  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!
! Baum, B.A., P. Yang, A.J. Heymsfield, A.R. Bansemer, B.H. Cole, A. Merrelli, C.G. Schmitt, 
!    and C. Wang, 2014: Ice cloud single-scattering property models with the full phase matrix 
!    at wavelengths from 0.2 to 100 µm. Journal of Quantitative Spectroscopy & Radiative Transfer
!    , 146, 123-139, DOI: 10.1016/j.jqsrt.2014.02.029.
!
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: io_icescat   !io_file name

 integer :: i,j,k,n  !loop
 integer :: ncid     !# of variables to inquire  
 integer :: varid    !# of variable id
 integer :: ierr     !allocation stat
 logical :: missing_para
 logical :: fexist

!
! file name
!

 if(masterproc) print*, 'MSG read_icescat_visir: reading...'//trim(io_icescat)
 inquire(file=trim(io_icescat),exist=fexist)
 if( .not. fexist ) then 
!     print*,'MSG read_icescat_visir; Cannot Find Database ',trim(io_icescat)
     call stop_sdsu(myrank,'MSG read_icescat_visir; Cannot Find Database '//trim(io_icescat))
 endif

!
! allocate single scattering parameters
!
 if( .not. allocated(phase_angles) ) then

 allocate( &
  phase_angles(nScatAng) ,&  ![degrees]
  wavelengths(nWaveLen)  ,& ![microns]
  effective_diameter(nDeff) ,& ![microns]
  ice_water_content(nWaveLen, nDeff) ,&  !ice water content  [gm-3]
  total_area       (nWaveLen, nDeff) ,& ! total projected area [microns^2]
  extinction_efficiency(nWaveLen, nDeff)  ,&
  single_scattering_albedo(nWaveLen, nDeff) ,&
  extinction_coefficient_over_iwc(nWaveLen, nDeff) ,& ![ m^2 g-1 ]
  asymmetry_parameter(nWaveLen, nDeff)  ,&
  p11_phase_function(nWaveLen, nDeff, nScatAng),&
!  p21_phase_function(nWaveLen, nDeff, nScatAng) ,&
  p22_phase_function(nWaveLen, nDeff, nScatAng) ,&  !needed for depolarization ratio
!  p33_phase_function(nWaveLen, nDeff, nScatAng) ,&
!  p43_phase_function(nWaveLen, nDeff, nScatAng) ,&
!  p44_phase_function(nWaveLen, nDeff, nScatAng) ,&
 stat=ierr )

 if(ierr /= 0) call stop_sdsu(myrank,'MSG read_icescat_visir: cannot allocate the variable')
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
! effective diameter (nDeff) [micron]
!
 call check( nf90_inq_varid(ncid,'effective_diameter',varid ),'effective_diameter', missing_para)
 call check( nf90_get_var(ncid, varid, effective_diameter ))

!
! ice water content  (nWaveLen, nDeff) [g/m3]
!
 call check( nf90_inq_varid(ncid,'ice_water_content',varid ),'ice_water_content', missing_para)
 call check( nf90_get_var(ncid, varid, ice_water_content ))

!
!  total projected area (nWaveLen, nDeff) [micron2]
!
 call check( nf90_inq_varid(ncid,'total_area',varid ),'total_area', missing_para)
 call check( nf90_get_var(ncid, varid, total_area ))

!
! extinction_efficiency (nWaveLen, nDeff) [-] 
!
 call check( nf90_inq_varid(ncid,'extinction_efficiency',varid ),'extinction_efficiency', missing_para)
 call check( nf90_get_var(ncid, varid, extinction_efficiency ))

!
!   extinction_coefficient_over_iwc(nWaveLen, nDeff) ,& ![ m^2 g-1 ]
!  (nWaveLen, nDeff) [-] 
!
 call check( nf90_inq_varid(ncid,'extinction_coefficient_over_iwc',varid))
 call check( nf90_get_var(ncid, varid, extinction_coefficient_over_iwc   ))


!
! single_scattering_albedo (nWaveLen, nDeff) [-] 
!
 call check( nf90_inq_varid(ncid,'single_scattering_albedo',varid ),'single_scattering_albedo', missing_para)
 call check( nf90_get_var(ncid, varid, single_scattering_albedo ))

!
! asymmetry_parameter (nWaveLen, nDeff) [-] 
!
 call check( nf90_inq_varid(ncid,'asymmetry_parameter',varid ),'asymmetry_parameter', missing_para)
 call check( nf90_get_var(ncid, varid, asymmetry_parameter ))

!
! p11_phase_function (nWaveLen, nDeff, nScatAng) [-] 
!
 call check( nf90_inq_varid(ncid,'p11_phase_function',varid ),'p11_phase_function', missing_para)
 call check( nf90_get_var(ncid, varid, p11_phase_function ))

!
! p22_phase_function (nWaveLen, nDeff, nScatAng) [-] 
!
 call check( nf90_inq_varid(ncid,'p22_phase_function',varid ),'p22_phase_function', missing_para)
 call check( nf90_get_var(ncid, varid, p22_phase_function ))


 return
 end subroutine read_icescat_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine  store_icescat_visir( mxwavel, wavel, nang, ang )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine interpolate Bryan Baum's Ice Scattering database, compatible to R-star
! (G-SDSU visible-IR simulator) units, and stored it in memory. 
!
! History:
! 07/2015  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: mxwavel !total channel of visir simulator
 real,intent(in)    :: wavel(mxwavel) ! channel wavelength [micron]
 integer,intent(in) :: nang
 real,intent(in)    :: ang(1:nang)  !angle for phase function [deg]
 integer :: ierr
 integer :: iw, nw, nd
 real :: wgt1, wgt2
 integer :: iang, nsang

!
! allocate single scattering parameters
!
 if( .not. allocated(ph_icescat) ) then

   allocate( &
    cext_icescat (nWaveLen, nDeff ),&
    cabst_icescat(nWaveLen, nDeff ),&
    ph_icescat   (nWaveLen, nDeff,nang), &
   stat=ierr )
   if(ierr /= 0) call stop_sdsu(myrank,'MSG store_icescat_visir: cannot allocate the variable')


 endif

 nang_rstar = nang


!
!
! change unit and phase angle consistent to R-star
!

  cext_icescat(1:nWaveLen,1:nDeff) = &  !total extinction cross sectoin per unit mass [cm2/g]
             extinction_coefficient_over_iwc(1:nWaveLen,1:nDeff) * 1.e+4 

  cabst_icescat(1:nWaveLen,1:nDeff) =  & !total absorption cross section per unit mass [cm2/g]
             (1.- single_scattering_albedo(1:nWaveLen,1:nDeff)) * cext_icescat(1:nWaveLen,1:nDeff)

!
! interpolate phase function from icescat phase angle to rstart phase angle
!
  ph_icescat(1:nWaveLen, 1:nDeff,    1) = p11_phase_function(1:nWaveLen, 1:nDeff, 1)       !0deg forward
  ph_icescat(1:nWaveLen, 1:nDeff, nang) = p11_phase_function(1:nWaveLen, 1:nDeff, nScatAng)!180deg backscatter

  ang_rstar: do iang = 2, nang-1  !r-star angle loop

     ang_icescat:  do nsang = 1, nScatAng - 1  !icescat database angle loop    

       if( ang(iang) >= phase_angles(nsang) .and. &
           ang(iang) <  phase_angles(nsang+1) ) then

            wgt1 = ( phase_angles(nsang+1) - ang(iang) ) / &
                   ( phase_angles(nsang+1) - phase_angles(nsang) ) 
            wgt2 = 1. - wgt1

            ph_icescat(1:nWaveLen, 1:nDeff, iang) = &  ! phase function  [-]
                    wgt1*p11_phase_function(1:nWaveLen,1:nDeff, nsang  ) &
                  + wgt2*p11_phase_function(1:nWaveLen,1:nDeff, nsang+1)
            exit ang_icescat
       endif

     enddo ang_icescat

  enddo ang_rstar


 return
 end subroutine  store_icescat_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine interp_icescat_visir(q, re, wlc, cextp, cabsp, php, success )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
! This routine calculate single-scattering properties for visible-IR simulator in G-SDSU by 
! interpolating Bryan Baum's Ice Scattering database.
!
! History:
! 07/2015  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 real(8),intent(in):: q         ! mass mixing ratio of particle [g/cm3]
 real(8),intent(in):: re        ! particle effective radius [micron]
 real(8),intent(in) :: wlc      ! wavelength [cm]
 real(8),intent(out):: cextp    ! extinction cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: cabsp    ! absorption cross section [cm-1] = [cm2/cm3]
 real(8),intent(out):: php(:)   ! volume scattering phase function [cm-1 str-1]
 logical, intent(out):: success ! interpolation sucess

 real :: q_real  ![g/cm3]
 real :: de_real ![micron]
 real :: wave
 integer :: nw, nw_int, nd
 real :: wgt1_deff, wgt2_deff
 real :: wgt1_wave, wgt2_wave

!
! initialize outpout
!
 success = .false.
 cextp = 0.e0 ; cabsp = 0.e0 ; php = 0.e0

!
! convert units and precision
! 
 q_real = real(q)         ! ice mixing ratio [g/cm3]
 de_real = real(re) * 2.  ! effective diameter [um]
 wave = real(wlc) * 1.e+4 ! wavelength [um]

!
! check effective diameter is within database anchor points.
! 
 if( de_real < effective_diameter(1) .or. de_real > effective_diameter(nDeff) ) return

!
! check wavelength is within database anchor points
!
 if( wave < wavelengths(1) .or. wave > wavelengths(nWaveLen) ) return


!
! interpolate single scattering properties for a given effective diameter
!
 ndeff_loop: do nd = 1, nDeff - 1

  ndeff_if: if( de_real >= effective_diameter(nd) .and. &
      de_real <= effective_diameter(nd+1) ) then

        !weighting points of Deff
        wgt1_deff = ( effective_diameter(nd+1) - de_real ) / &
               ( effective_diameter(nd+1) - effective_diameter(nd) ) 
        wgt2_deff = 1.-wgt1_deff

        wave_loop:  do nw = 1, nWaveLen-1  !~445

           wave_if: if( wave >= wavelengths(nw) .and. wave <= wavelengths(nw+1) ) then  ! interpolate

             ! weighting points of wavelength
              wgt1_wave = ( wavelengths(nw+1) - wave)  / ( wavelengths(nw+1) - wavelengths(nw) )
              wgt2_wave = 1. - wgt1_wave


             ! units conversion [cm2/cm3]=   [cm2/g] * [g/cm3]
              cextp = &
               (  ( wgt1_deff*cext_icescat(nw  ,nd) + wgt2_deff*cext_icescat(nw  ,nd+1) ) * wgt1_wave &
                 +( wgt1_deff*cext_icescat(nw+1,nd) + wgt2_deff*cext_icescat(nw+1,nd+1) ) * wgt2_wave ) * q_real

              cabsp = &
               (  ( wgt1_deff*cabst_icescat(nw  ,nd) + wgt2_deff*cabst_icescat(nw  ,nd+1) ) * wgt1_wave &
                 +( wgt1_deff*cabst_icescat(nw+1,nd) + wgt2_deff*cabst_icescat(nw+1,nd+1) ) * wgt2_wave ) * q_real

             !phase function weight by scattering cross section  [cm-1 str-1]
              php(1:nang_rstar) = (cextp-cabsp) * (  &
              ( wgt1_deff*ph_icescat(nw  ,nd,1:nang_rstar) + wgt2_deff*ph_icescat(nw  ,nd+1,1:nang_rstar) ) * wgt1_wave&
            + ( wgt1_deff*ph_icescat(nw+1,nd,1:nang_rstar) + wgt2_deff*ph_icescat(nw+1,nd+1,1:nang_rstar) ) * wgt2_wave&
               ) / ( 4.*3.14159265)

!              print*,'MSG icescat:  ', wgt1_deff, wgt2_deff, de_real, wgt1_wave, wgt2_wave, wave
               success = .true.
               exit ndeff_loop
            endif wave_if

        enddo wave_loop

  endif ndeff_if

 enddo ndeff_loop


 return
 end subroutine interp_icescat_visir

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

 end module module_icescat_visir

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

