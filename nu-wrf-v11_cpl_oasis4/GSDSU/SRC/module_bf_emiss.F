#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_bf_emiss
 use module_floating_point
 use module_const
 use netcdf
 use module_mpi
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
!  This is an IO interface module to handle Baseline Fit (BF) Global Infrared 
!  Land Surface Emissivity from CIMSS, University of Wisconsin.
!
! Database:  (http://cimss.ssec.wisc.edu/iremis/)
!  This global database of infrared land surface emissivity is derived using input from the 
!  Moderate Resolution Imaging Spectroradiometer (MODIS) operational land surface emissivity 
!  product (MOD11). The baseline fit method (Seemann et al., 2007), based on a conceptual model 
!  developed from laboratory measurements of surface emissivity, is applied to fill in the 
!  spectral gaps between the six emissivity wavelengths available in MOD11. The six available 
!  MOD11 wavelengths span only three spectral regions (3.8-4, 8.6, and 11-12 microns), while 
!  the retrievals of atmospheric temperature and moisture from satellite IR sounder radiances 
!  require surface emissivity at higher spectral resolution. Emissivity in the baseline fit 
!  database presented here is available globally at ten wavelengths (3.6, 4.3, 5.0, 5.8, 7.6, 
!  8.3, 9.3, 10.8, 12.1, and 14.3 microns) with 0.05 degree spatial resolution. The wavelengths 
!  in the database were chosen as hinge points to capture as much of the shape of the higher 
!  resolution emissivity spectra as possible between 3.6 and 14.3 microns.
!
!
! History:
!  07/2014  Toshi Matsui@NASA GSFC ; Initial
!
!
! References:
!  Seemann, S.W., E. E. Borbas, R. O. Knuteson, G. R. Stephenson, H.-L. Huang, 2008: Development 
!    of a Global Infrared Land Surface Emissivity Database for Application to Clear Sky Sounding 
!    Retrievals from Multi-spectral Satellite Radiance Measurements.J. Appl. Meteor. Climatol., 
!    Vol. 47, 108-123.
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
 public  :: read_bf_emiss, &      ! read emissivity files and allocate sub-domain value in memory
            interp_bf_emiss, &    ! interpolate emissiivty for a given lat/lon and wavelength 
            interp_bf_emiss_all   ! interpolate all-channel emisssivity for a given lat/lon. 



!
! module global parameters
!

 integer,parameter,public :: anchor_bs_emiss = 10 !anchor wavelength number
 real,parameter,public :: wave_bs_emiss(anchor_bs_emiss) = &  ! anchor wavelengths [micron]
        (/3.6, 4.3, 5.0, 5.8, 7.6, 8.3, 9.3, 10.8, 12.1, 14.3/)
 character(len=4),parameter,public :: wave_bs_emiss_cha(anchor_bs_emiss) = &  ! anchor wavelengths 
        (/' 3.6',' 4.3',' 5.0',' 5.8',' 7.6',' 8.3',' 9.3','10.8','12.1','14.3'/)


!
! original database parameters
!
 integer,parameter :: imax_bs_emiss=7200, jmax_bs_emiss=3600  !database domain index
 integer :: is_data, ie_data, js_data, je_data                !subsetted database domain index
 real,parameter :: dgrid_bs_emiss = 0.05  !grid size of database [deg]

 real,allocatable :: bf_emiss(:,:,:)  !baseline-fitted emissivity value [-] (j,i,wavelength)
                                      !j --> lat, i-->lon 
                                      !this sub-domain emissivity is stored in memory. 


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine read_bf_emiss(io_emiss,lat_model,lon_model)
 implicit none
!-----------------------------------------------------------------------------------------------------
! Comments: 
!  This subroutine read land surface emissivity over subset domain, corresponding to model domain. 
!  Land Surface Emissivity from CIMSS, University of Wisconsin.
!
! BF EMISSIVITY DATA FORMAT:
!
! The global emissivity data are available in NETCDF format, with one file per 
! month. 
!
! Each file contains monthly emissivity data at 10 inflection points on a 
! 0.05 x 0.05 degree grid (lon from -180deg ~ 180deg).  The wavelengths saved in the files are 
! 3.6, 4.3, 5.0, 5.8, 7.6, 8.3, 9.3, 10.8, 12.1, and 14.3 microns. 
! Since these 10 wavelengths were chosen as inflection points to capture as much of the shape 
! of the higher resolution spectra as possible, these values can be interpolated between inflection 
! points to arrive at any emissivity between 3.6 and 14.3 microns.
!
! The monthly files also contain a flag for each grid point which gives 
! information about how the emissivity data was derived, as explained below:
!
! flag = 0  	for ocean/coastlines where there was no MYD11 data, and thus no 
!		baseline fit-derived emissivity.
!
! flag = 1  	for regions where the baseline fit method was used to compute 
!		emissivity.
!
! flag = 2  	where MYD11 data was not available at this lat/lon, and the 
!		emissivity average over the two adjacent months (before and after) 		
!		for the corresponding lat/lon and wavelength was used.
!
! flag = 3  	where no MYD11 data was available in either adjacent month for the 
!		corresponding lat/lon, and emissivity was filled in with the 
!		2003 annual average for this lat/lon and wavelength.
!
! flag = 4  	where no MYD11 data was available all year for this lat/lon (this 
!		occurs in some Antarctic regions only due to MYD11 algorithm 			
!		constraints).  Emissivity was filled in with the average of the 
!		annual average over all emissivity at a given wavelength with 			
!		latitude < -80S.
!
! History:
!  07/2014  Toshi Matsui@NASA GSFC ; Initial
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: io_emiss                 !io_file name
 real,intent(in),dimension(:,:) :: lat_model, lon_model  !lat and lon of model in assumed dimension

 integer :: i,j,w  !looping
 integer :: ncid                !# of variables to inquire  
 integer :: varid               !# of variable id
 integer :: bnd(2)  !upper memory bound for 3D array

 real ::  min_lat,max_lat,min_lon,max_lon  ! model domain bondary lat/lon [deg]  
 real,allocatable :: tmp2d(:,:)
 character(len=2) :: char_wave
 character(len=10) :: para_char


!
! latitude and longitude of data
!
 bnd = UBOUND(lat_model)  ! find upper bound from assumped 2D array

!
! if domain croos -180deg (or 180deg) longitude  --> stop 
!
 if( ( lon_model(1,     1) > 0.e0 .and. lon_model(bnd(1),     1) < 0.e0 ) .or. &
     ( lon_model(1,bnd(2)) > 0.e0 .and. lon_model(bnd(1),bnd(2)) < 0.e0 ) ) then
    call stop_sdsu(myrank,'MSG read_bf_emiss: Warning Ooops domain cross -180deg Not supported in this code.' )
 endif

!
! estimate minimum and maximum domain of model
!
 min_lat = minval( lat_model ) ! minimum latitude [deg]  
 max_lat = maxval( lat_model ) ! maximum latitude [deg]
 min_lon = minval( lon_model ) ! minimum longitude [deg]  
 max_lon = maxval( lon_model ) ! maximum longitude [deg]


!
! estimate starting and ending i j index to read subset of database.
!
 js_data = INT(  (  89.975 - max_lat) / dgrid_bs_emiss    ) + 1
 je_data = INT(  (  89.975 - min_lat) / dgrid_bs_emiss    ) + 1
 is_data = INT( -(-179.975 - min_lon) / dgrid_bs_emiss    ) + 1
 ie_data = INT( -(-179.975 - max_lon) / dgrid_bs_emiss    ) + 1

!
! extend domain with extra five grid in case. 
!
 js_data = min( max( js_data - 5 , 1 ) , jmax_bs_emiss )   
 je_data = min( max( je_data + 5 , 1 ) , jmax_bs_emiss )    
 is_data = min( max( is_data - 5 , 1 ) , imax_bs_emiss )   
 ie_data = min( max( ie_data + 5 , 1 ) , imax_bs_emiss )


!
! define IO directory and file
!
 if(masterproc) print("(a)"),'MSG read_bf_emiss: input ->',trim(io_emiss)


!
! allocate global parameters
!
 if( .not. allocated(tmp2d) ) then
    allocate( tmp2d   (js_data:je_data,is_data:ie_data) , &
              bf_emiss(js_data:je_data,is_data:ie_data,anchor_bs_emiss)  )
 endif

!
! open netCDF file
!
 call check( nf90_open(trim(io_emiss), nf90_nowrite, ncid) )

!
! read emissivity data 
!
 do w = 1,  10 !anchor_bs_emiss

    if( w >= 10 )  then
       write(char_wave,"(I2.2)") w
       para_char = 'emis'//char_wave
    else
       write(char_wave,"(I2.1)") w
       para_char = 'emis'//char_wave(2:2)
    endif

    !
    ! read value of subsetted domain
    !
    call check( nf90_inq_varid(ncid, trim(para_char) , varid ) ) !IR emissivity [-]
    call check( nf90_get_var(ncid, varid,   tmp2d(js_data:je_data,is_data:ie_data), &
                start=(/js_data,is_data/), count=(/je_data-js_data+1,ie_data-is_data+1/)) )

    do j = js_data,je_data ; do i = is_data,ie_data 
       if( tmp2d(j,i) == -9999.0) then     
          bf_emiss(j,i,w) = undefined
       else
          bf_emiss(j,i,w) = tmp2d(j,i) * 0.001  !emissivity of sub-domain [-]
       endif
    enddo ; enddo

 enddo

!
! close netCDF file
!
  call check( nf90_close(ncid) )

 end subroutine read_bf_emiss

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

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine interp_bf_emiss(wave,lat_in,lon_in, emis_ir)
 implicit none
!-----------------------------------------------------------------------------------------------------
! Comments: 
!  This subroutine spatial interpolate land-surface IR emissivity database for a given lat-lon point
!  for visir simulator. 
!  Currently, this employ 4pixel weighted interpolation. So, emissivity database is 0.05deg (~5km)
!  , meaining CRM grid size should not be greater than 10km or so. 
!
! History:
!  07/2014  Toshi Matsui@NASA GSFC ; Initial
!-----------------------------------------------------------------------------------------------------
 real,intent(in) :: wave   !wavelength [micron]
 real,intent(in) :: lat_in,lon_in !latitude and longitude [deg]
 real,intent(out) :: emis_ir  !interpolated emissivity [-]

 integer :: i,j, w
 integer :: i_data, j_data
 real :: i_real, j_real
 real ::  dist, wgt_sum, wgt1, wgt2
 real :: bf_sum(anchor_bs_emiss)
 real :: bf_emiss_interp(anchor_bs_emiss)

!
! position in floating point. 
!
 j_real =   (  89.975 - lat_in) / dgrid_bs_emiss   + 1.e0
 i_real =  -(-179.975 - lon_in) / dgrid_bs_emiss   + 1.e0

!
! spatial interpolation from nearest 4 point. 
!
 wgt_sum = 0.e0 ; bf_sum = 0.e0
 do i = int(i_real), int(i_real)+1  ; do j = int(j_real), int(j_real)+1
    if( bf_emiss(j,i,1 ) /= undefined ) then
       dist = sqrt( (j_real-real(j))**2 + (i_real-real(i))**2 )
       wgt_sum = wgt_sum + dist
       bf_sum(:)  = bf_sum(:) +  bf_emiss(j,i,:) * dist
    endif
 enddo ; enddo

 if( wgt_sum > 0.e0 ) then
   bf_emiss_interp(:) = bf_sum(:) / wgt_sum 
 else
   !
   ! water pixels surrouned in database  --> so it later treated as water in module_opt_visir.F
   !
   bf_emiss_interp(:) = undefined 
   emis_ir = undefined
   return
 endif


!
! interpolate for wavelength
!
 if    ( wave <= wave_bs_emiss(1) ) then
    emis_ir = bf_emiss_interp(1)   ! below 3.6micron case 
 elseif( wave >= wave_bs_emiss(anchor_bs_emiss) ) then
    emis_ir = bf_emiss_interp(anchor_bs_emiss)  !abvoe 14.3micron case
 else !interpolate

  do w = 1, anchor_bs_emiss-1
     if( wave >= wave_bs_emiss(w) .and. wave <= wave_bs_emiss(w+1) ) then
         wgt1 = (wave_bs_emiss(w+1) - wave) / (wave_bs_emiss(w+1)-wave_bs_emiss(w))
         wgt2 = 1.e0 - wgt1 
         emis_ir = wgt1 * bf_emiss_interp(w) + wgt2 * bf_emiss_interp(w+1)
     endif
  enddo 

 endif

 end subroutine interp_bf_emiss

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine interp_bf_emiss_all(lat_in,lon_in, bf_emiss_interp)
 implicit none
!-----------------------------------------------------------------------------------------------------
! Comments: 
!  This subroutine spatial interpolate land-surface IR emissivity database for a given lat-lon point
!  for sarta simulator (use all anchor-wavelength values). 
!  Currently, this employ 4pixel weighted interpolation. So, emissivity database is 0.05deg (~5km)
!  , meaining CRM grid size should not be greater than 10km or so.
!
! History:
!  07/2014  Toshi Matsui@NASA GSFC ; Initial
!-----------------------------------------------------------------------------------------------------
 real,intent(in) :: lat_in,lon_in !latitude and longitude [deg]
 real,intent(out) :: bf_emiss_interp(anchor_bs_emiss) !spatially interpoltaed spectrum emissivity [-]

 integer :: i,j, w
 integer :: i_data, j_data
 real :: i_real, j_real
 real ::  dist, wgt_sum, wgt1, wgt2
 real :: bf_sum(anchor_bs_emiss)


!
! position in floating point. 
!
 j_real =   (  89.975 - lat_in) / dgrid_bs_emiss   + 1.e0
 i_real =  -(-179.975 - lon_in) / dgrid_bs_emiss   + 1.e0

!
! spatial interpolation from nearest 4 point. 
!
 wgt_sum = 0.e0 ; bf_sum = 0.e0
 do i = int(i_real), int(i_real)+1  ; do j = int(j_real), int(j_real)+1
    if( bf_emiss(j,i,1 ) /= undefined ) then
       dist = sqrt( (j_real-real(j))**2 + (i_real-real(i))**2 )
       wgt_sum = wgt_sum + dist
       bf_sum(:)  = bf_sum(:) +  bf_emiss(j,i,:) * dist
    endif
 enddo ; enddo

 if( wgt_sum > 0.e0 ) then
   bf_emiss_interp(:) = bf_sum(:) / wgt_sum 
 else
   !
   ! water pixels surrouned in database  --> so it later treated as water in module_opt_visir.F
   !
   bf_emiss_interp(:) = undefined 
   return
 endif


 end subroutine interp_bf_emiss_all

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end  module module_bf_emiss

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
