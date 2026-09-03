!###############################################################################!
!  CCATT-BRAMS/MCGA-CPTEC/WRF emission model    CPTEC/INPE                      !
!  Version 1.0.0: 12/nov/2010                                                   !
!  Coded by Saulo Freitas and Karla Longo                                       !
!  Contact: gmai@cptec.inpe.br - http://meioambiente.cptec.inpe.br              !
!###############################################################################!

module gocart_backgr
!---------------------------------------------------------------------------
integer, parameter :: maxnspecies = 200 , nspecies = 6 ! NUWRF...Added PS

integer, parameter ::         &
 DMS		  =1  &
,EROD		  =2  &
,H2O2		  =3  & 
,OH		  =4  &
,NO3		  =5  &
,PS               =6  ! NUWRF...Added PS

!---------------------------------------------------------------------------
character(LEN=25),dimension(nspecies),parameter :: spc_name= &
! '1234567890123456789012345'
! EMK NUWRF...Added spaces to appease gfortran
(/                          & 
  'DMS                      '&
, 'EROD                     '&
, 'H2O2                     '&
, 'OH                       '&
, 'NO3                      '&
, 'PS                       '&  ! NUWRF...Added PS
/)

integer :: nvlev_go ! NUWRF

! NUWRF...nlevels_netcdf now set in subroutine init_gocart_bg
integer,dimension(nspecies) :: nlevels_netcdf
!(/12 &  !12 months
! ,3  &  ! 3 layers
! ,nvlev_go &  ! 55 vertical levels
! ,nvlev_go &
! ,nvlev_go &
! ,nvlev_go &
!/)

real, allocatable, dimension(:) :: lev ! NUWRF...lev now allocatable


! NUWRF...netcdf_spc_name now set in subroutine init_gocart_bg
character(LEN=25),dimension(nspecies) :: netcdf_spc_name
! '1234567890123456789012345'
!(/                          & 
!  'DMSO '&
!, 'EROD '&
!, 'H2O2 '&
!, 'OH   '&
!, 'NO3  '&
!, 'PS   '&
!/)

!REAL,PARAMETER,DIMENSION(nspecies) :: convert_to_kg_per_day=(/&
!    1.e6/365.  ,   & ! OC  => Gg/year per grid box = 1.e6/365. kg/day
!    1.e6/365.  ,   & ! BC  => Gg/year per grid box = 1.e6/365. kg/day
!    1./365.    ,   & ! SO2 => kg/year per grid box = 1./365.   kg/day
!    1.e6/365.  ,   & ! BC  => Gg/year per grid box = 1.e6/365. kg/day
!    1./365.        & ! SO2 => kg/year per grid box = 1./365.   kg/day
!/)
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------

  type gocart_bg_vars   
     real, pointer, dimension(:,:,:)  :: src
!-----------

  end type gocart_bg_vars

  type (gocart_bg_vars), allocatable :: gocart_bg_g(:)

contains
  !---------------------------------------------------------------
  ! NUWRF...New subroutine for setting netcdf_spc_name, lev, and
  ! nlevels_netcdf depending on whether 'old' or 'new' gocart background
  ! data are read in.
  subroutine init_gocart_bg(n1,n2,gocart_bg_data_type)
     implicit none
     integer :: n1,n2
     character(len=3) :: gocart_bg_data_type    

     if(trim(gocart_bg_data_type) == "new") then
        ! NUWRF...Name of PS could be overwritten later if non-2006 data is
        ! used.
        nvlev_go = 72
        netcdf_spc_name = &
          (/              & 
             'DMSO '      &
           , 'EROD '      &
           , 'H2O2 '      &
           , 'OH   '      &
           , 'NO3  '      &
           , 'PS   '      &
         /)
    else
        nvlev_go = 55
        netcdf_spc_name = &
          (/              & 
             'DMSO '      &
           , 'EROD '      &
           , 'h2o2 '      &
           , 'oh   '      &
           , 'no3  '      &
           , 'ps   '      &
         /)
    endif

    nlevels_netcdf = &
     (/12 &  !12 months
      ,3  &  ! 3 layers
      ,nvlev_go &  ! 55 vertical levels
      ,nvlev_go &
      ,nvlev_go &
      ,nvlev_go &
     /)

    call mem_gocart_bg(n1,n2,nvlev_go)

! EMK NUWRF: Make sure lev isn't cleared out for inner grids.
    if (.not. allocated(lev)) then
       allocate(lev(nvlev_go))
       lev = 0.0
    end if

  end subroutine init_gocart_bg

  !---------------------------------------------------------------

  subroutine alloc_gocart_bg(gocart_bg,n1,n2,n3)

    implicit none

    type (gocart_bg_vars),dimension(nspecies)  :: gocart_bg
    integer,intent(in) :: n1,n2,n3
    integer ispc, max_nlevels_netcdf
   
!    max_nlevels_netcdf = max(nlevels_netcdf)
   
    do ispc=1,nspecies
     allocate (gocart_bg(ispc)%src(n1,n2,n3))
     gocart_bg(ispc)%src = 0.
    enddo

    return
  end subroutine alloc_gocart_bg

  !---------------------------------------------------------------

  subroutine nullify_gocart_bg(gocart_bg)

    implicit none

    type (gocart_bg_vars),dimension(nspecies)  :: gocart_bg
    integer ispc

    do ispc=1,nspecies
       if (associated(gocart_bg(ispc)%src))    nullify (gocart_bg(ispc)%src)
    enddo

    return
  end subroutine nullify_gocart_bg

!---------------------------------------------------------------
!---------------------------------------------------------------
  ! NUWRF...Added mem_gocart_bg to module.  Added n3 argument.
  subroutine mem_gocart_bg(n1,n2,n3)
!    use gocart_backgr ! NUWRF...Not needed
    implicit none
    integer i
    integer, intent(in) :: n1,n2,n3 ! NUWRF...Added n3

    if(.not. allocated(gocart_bg_g)) allocate(gocart_bg_g(nspecies))
    !do i=1,nspecies
    ! if(associated(gocart_bg_g(i)%src)) deallocate(gocart_bg_g(i)%src)
    !enddo

    call nullify_gocart_bg(gocart_bg_g(:))      
    call alloc_gocart_bg  (gocart_bg_g(:),n1,n2,n3) ! NUWRF changed last dimension
  end subroutine mem_gocart_bg

end module gocart_backgr

!---------------------------------------------------------------
!---------------------------------------------------------------

! NUWRF...Added gocart_bg_data_type argument.  
subroutine read_gocart_bg(iyear,imon,iday,ng,ngrids,n1,n2,n3,rlat,rlon,rland,deltax,deltay&
                            ,xt,yt,xm,ym,plat,plon,gocart_bg_data_type)
use grid_dims_out, only: grid_type, gocart_bg_data_dir
use gocart_backgr
use mem_grid, only :grid_g
use netcdf
implicit none
!include 'netcdf.inc'
integer, parameter ::  nlon = 288, nlat=181, nmonths=12
integer, parameter ::  ndlon = 288, ndlat=181 !ndlon = 1440, ndlat=720
integer, intent (in) :: iyear,imon,iday,ng,n1,n2,n3,ngrids
real, intent (in), dimension(n1,n2)  :: rlat,rlon,rland
real, intent (in) :: deltax,deltay
real,intent(in) ::  xt(n1), yt(n2),xm(n1), ym(n2)
integer :: ispc,im,i,k,iread,nc,ilaea,jlaea,kk,k1,k2,ii,i1,i2,j,imonx
integer :: ic,jc,j1,j2,var_id,ncid
character*240 prefix,suffix,filename(nspecies)
character*10 dummy
character*3 gocart_bg_data_type ! NUWRF
character*2 cmon(12)
real longgocart_bg(nlon),latgocart_bg(nlat)
real longgocart_dust(ndlon),latgocart_dust(ndlat)
!NUWRF EMK...Fixed resolution of EROD data
!real, parameter:: ilatn=1.,ilonn=1.25, idlonn = 0.25,idlatn = 0.25
real, parameter:: ilatn=1.,ilonn=1.25, idlonn = 1.25,idlatn = 1.
 data (cmon(i),i=1,12) /'01','02','03','04','05','06','07','08','09','10',  &
                       '11','12'/

!real lat,lon,src_dummy(nlon,nlat,maxval(nlevels_netcdf))
real lat,lon
real         src_dust (ndlon,ndlat,3,1)
real         src_dust_dummy (ndlon,ndlat,nvlev_go) ! NUWRF edit
! real rrlat,rrlon,dlon1,dlon2,dlat1,dlat2,TX(maxval(nlevels_netcdf),nspecies),plat,plon,area
! NUWRF...Changed first dim for TX
real rrlat,rrlon,dlon1,dlon2,dlat1,dlat2,TX(nvlev_go,nspecies),plat,plon,area 

real, parameter ::                    &
        pi180    = 3.1415927 / 180.   &
    ,   r_earth  = 6370000.

real, allocatable :: src_dummy(:,:,:)
real,allocatable ,save :: RAWsrc(:,:,:,:)
real,allocatable ,save :: RAWdustsrc(:,:,:,:)

integer :: a_i, a_len, a_inds(1), a_err

! NUWRF...New variables for calculating 3D pressure
real          p0GMI
real          hyamGMI   (nvlev_go)
real          hybmGMI   (nvlev_go)
real          surfPress (nlon,nlat)
real          press3c   (nlon,nlat,nvlev_go)

character(len=4) :: cyear ! EMK NUWRF
integer :: ierr ! EMK NUWRF
logical :: found_curyear ! EMK NUWRF
character(len=240) :: test_filename ! EMK NUWRF

    print *, "IN READ_GOCART_BG...", nvlev_go, gocart_bg_data_type
!--- for gocart_bg there is not monthly variation, 
imonx=1  

!--- lat e lon gocart_bg (corner mais ao sul e mais a oeste)		 		      
do k=1,nlon;  longgocart_bg(k)=-180. + (k-1)*ilonn; enddo
do i=1,nlat;  latgocart_bg(i)= -89.75 + (i-1)*ilatn; enddo
!do k=1,ndlon;  longgocart_dust(k)=-179.875 + (k-1)*idlonn; enddo
! NUWRF EMK...Fixed lat/lon of EROD
!do i=1,ndlat;  latgocart_dust(i)= -89.875 + (i-1)*idlatn; enddo
!do k=1,ndlon;  longgocart_dust(k)=-179.875 + (k-1)*idlonn; enddo
do i=1,ndlat;  latgocart_dust(i)=   -90 + (i-1)*idlatn; enddo
do k=1,ndlon;  longgocart_dust(k)= -180 + (k-1)*idlonn; enddo

if( .not. allocated (RAWsrc)) then
!  allocate(  RAWsrc(nlon,nlat,maxval(nlevels_netcdf),nspecies) )
   allocate(  RAWsrc(nlon,nlat,nvlev_go,nspecies) ) ! NUWRF...Changed third dim
   allocate(  RAWdustsrc(ndlon,ndlat,nvlev_go,nspecies) ) ! NUWRF...Changed third dim
   RAWsrc= 0.0
ENDIF

if(ng == 1) then  ! need to read just one time

  !- name of dataset
  suffix='.nc'

  ! EMK NUWRF
  ! If processing "new" data, we need to figure out if we are using standard
  ! 2006 data or more recent year/month specific data. Unfortunately the
  ! "more recent" data is split into two files. It is convenient to answer the
  ! question here before we start reading variables.
  found_curyear = .false. ! First guess
  if (trim(gocart_bg_data_type) == 'new' .or. &
      trim(gocart_bg_data_type) == 'NEW') then
     found_curyear = .true. ! Second guess, tested below
     write(cyear,'(I4.4)') iyear
     do i = 1, 2
        if (i == 1) then
           prefix= 'gocart_bg_new/geos5_met_1MAVG_'//cyear//cmon(IMON)
        else
           prefix= 'gocart_bg_new/gmi_merra_oxidants_'//cyear//cmon(IMON) &
                //'_1.25x1'
        end if
        test_filename=trim(gocart_bg_data_dir)//'/'//trim(prefix)//suffix
        ierr = nf90_open(trim(test_filename), NF90_NOWRITE, ncid)
        if (ierr .eq. NF90_NOERR) then ! File exists
           call check(nf90_close(ncid))
        else ! File not found. Will use 2006 data
           found_curyear=.false.
           print*,'WARNING, cannot open file ',trim(test_filename)
           print*, trim(nf90_strerror(ierr))
           print*,'Will switch to 2006 data.'
           exit ! Jump out of loop
        end if
     end do
  end if
  
  do ispc=1,nspecies

    nc=len_trim(spc_name(ispc))

    !print*,'nc=',nc,spc_name(ispc)
    if(trim( spc_name(ispc) ) == 'DMS') prefix= 'dms_data/dms_1x1.25'
    if(trim( spc_name(ispc) ) == 'EROD') prefix= 'erod/GAO_source_3cl'

    ! EMK NUWRF...Changed logic to use special subdirectory for 72-level data,
    ! which may be older 2006 data or newer year/month specific data.
    if(trim( spc_name(ispc) ) == 'H2O2' .or.&
       trim( spc_name(ispc) ) == 'OH'   .or.&
       trim( spc_name(ispc) ) == 'NO3'  .or.&
       trim( spc_name(ispc) ) == 'PS' ) then
       if (trim(gocart_bg_data_type) == 'new' .or. &
           trim(gocart_bg_data_type) == 'NEW') then
          if (found_curyear) then
             if ( trim (spc_name(ispc) ) == 'PS' ) then
                prefix= 'gocart_bg_new/geos5_met_1MAVG_'//cyear//cmon(IMON)
             else
                prefix= 'gocart_bg_new/gmi_merra_oxidants_'//cyear//cmon(IMON) &
                  //'_1.25x1'
             end if
          else
             prefix= 'gocart_bg_new/gmi_2006'//cmon(IMON)
          end if
       else
          prefix= 'gocart_bg/gmi_2006'//cmon(IMON)
       end if
    end if

    filename(ispc)=trim(gocart_bg_data_dir)//'/'//trim(prefix)//suffix
    
    print *,'opening   ', trim(filename(ispc)),' - specie= ',trim(spc_name(ispc))
    
! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
    call check( nf90_open(TRIM(filename(ispc)), NF90_NOWRITE, ncid) )
    
!     print*,'1',ncid,trim(netcdf_spc_name(ISPC))
!pause
! Get the varid of the data variable, based on its name.
! EMK NUWRF..."New" year/month specific data has lowercase ps, while "new"
! 2006 data has uppercase.
    if (found_curyear .and. trim(spc_name(ispc)) .eq. 'PS') then
       call check( nf90_inq_varid(ncid, 'ps', var_id) )
    else
       call check( nf90_inq_varid(ncid, trim(netcdf_spc_name(ispc)), var_id) )
    end if

    if(trim( spc_name(ispc) ) == 'EROD') then
      call check( nf90_get_var(ncid, var_id, src_dust  ) )
      src_dust_dummy(1:ndlon,1:ndlat,1:3) = src_dust(:,:,1:3,1)
      do i=1,ndlon;do j=1,ndlat
	    RAWdustsrc(i,j,1:nlevels_netcdf(ispc),ispc)=src_dust_dummy(i,j,1:nlevels_netcdf(ispc))
      enddo;enddo      
! NUWRF
    elseif(trim( spc_name(ispc) ) == 'PS') then
      ! Get the surface pressure
!       call check( nf_inq_varid(ncid, 'PS', var_id) )
       call check( nf90_get_var(ncid, var_id, surfPress ) )
      ! Get p0
       call check( nf90_inq_varid(ncid, 'p0', var_id) )
       call check( nf90_get_var(ncid, var_id, p0GMI ) )
      ! Get hyam
       call check( nf90_inq_varid(ncid, 'hyam', var_id) )
       call check( nf90_get_var(ncid, var_id, hyamGMI ) )
      ! Get hybm
       call check( nf90_inq_varid(ncid, 'hybm', var_id) )
       call check( nf90_get_var(ncid, var_id, hybmGMI ) )
      ! Compute the atmospheric pressure
      call CalcPress3dCenter(p0GMI, hyamGMI, hybmGMI, surfPress, press3c, nlon, nlat, nvlev_go)
      RAWsrc(1:nlon,1:nlat,1:nlevels_netcdf(ispc),ispc)= &
           press3c(1:nlon,1:nlat,1:nlevels_netcdf(ispc))
! NUWRF END
    else 
      allocate(src_dummy(nlon,nlat,nlevels_netcdf(ispc)))
      call check( nf90_get_var(ncid, var_id, src_dummy  ) )      
      do i=1,nlon;do j=1,nlat
	    RAWsrc(i,j,1:nlevels_netcdf(ispc),ispc)=src_dummy(i,j,1:nlevels_netcdf(ispc))
      enddo;enddo
!      print*,'EMK: maxval ',ispc,' ',trim( spc_name(ispc)), &
!           maxval(src_dummy(:,:,1:nlevels_netcdf(ispc)))
      deallocate(src_dummy)
    endif

    !-special section for reading 'lev' 
! NUWRF...Discriminate between "old" and "new" GOCART background fields
 if(trim(gocart_bg_data_type) == "old") then
    if(trim( spc_name(ispc) ) == 'H2O2') then
       call check( nf90_inq_varid(ncid, 'lev', var_id) )
       !  call check( nf90_get_var(ncid, var_id,lev ) )

      ! read the values for lev one layer at a time
      !do a_i = 1,nvlev_go      ! nlevels_netcdf(ispc)
       do a_i = 1,maxval(nlevels_netcdf)
          a_inds(1) = a_i
          call check(nf90_get_var(ncid, var_id, lev(a_i), a_inds ))
          !print*, 'NetCDF Error message: ', nf90_strerror( a_err )
          !print*, 'ind = ', a_i, 'lev(ind) = ', lev(a_i)
       end do
    endif
 endif

! Close the file, freeing all resources.
      call check( nf90_close(ncid) )
  enddo ! nspecies

! convert to mass/area (m^2) , only EROD
  do j=1,nlat
    area=cos(latgocart_dust(j)*pi180)* (r_earth**2) * idlatn * idlonn * pi180**2
 
    do i=1,nlon   
      do ispc=erod,erod
           RAWdustsrc(i,j,1:nlevels_netcdf(ispc),ispc)= RAWdustsrc(i,j,1:nlevels_netcdf(ispc),ispc)/area
      enddo
   enddo
  enddo 

endif ! ng==1

!--- performs the interpolation to model grid box
do i=1,n1
  do j=1,n2

    call get_index1(i,j,nlon,nlat,n1,n2,rlat,rlon,longgocart_bg,latgocart_bg &
                 ,ilatn, ilonn ,dlat1,dlat2,dlon1,dlon2,i1,j1,i2,j2,ic,jc)

!   print *,'1', i,j,n1,n2,ic,jc,nlat,nlon,maxval(nlevels_netcdf),ilatn,ilonn,imon,nmonths,nspecies

    ! EMK NUWRF...Added Lat/Lon of GOCART background fields
    call interpol2_gocart_bg(i,j,n1,n2,rlon,rlat,ic,jc,nlat,nlon,maxval(nlevels_netcdf)&
                   ,ilatn,ilonn ,imon,nmonths,nspecies,RAWsrc,tx(:,:) &
                   ,longgocart_bg,latgocart_bg)
!   call interpol2_gocart_bg(i,j,n1,n2,rlon,rlat,ic,jc,nlat,nlon,nvlev_go                   &
!                  ,ilatn,ilonn ,imon,nmonths,nspecies,RAWsrc,tx(:,:))
! 
! TX is the value interpolated to the model grid box.
!-obs: convert to kg / day
    do ispc = 1, nspecies 
     if(ispc .ne. erod) then
     gocart_bg_g(ispc)%src(i,j,1:nlevels_netcdf(ispc))=TX(1:nlevels_netcdf(ispc),ispc) 
     endif

    enddo

  enddo
enddo 

! print *,'interpolate dust to grid'
!--- performs the interpolation to model grid box
do i=1,n1
  do j=1,n2

    call get_index1(i,j,ndlon,ndlat,n1,n2,rlat,rlon,longgocart_dust,latgocart_dust &
                 ,idlatn, idlonn ,dlat1,dlat2,dlon1,dlon2,i1,j1,i2,j2,ic,jc)

! print *,'grid point',ic,jc
    ! NUWRF EMK...Added Lat/Lon of GOCART background fields
    call interpol2_gocart_bg(i,j,n1,n2,rlon,rlat,ic,jc,ndlat,ndlon,maxval(nlevels_netcdf) &
                   ,idlatn,idlonn ,imon,nmonths,nspecies,RAWdustsrc,tx(:,:) &
                   ,longgocart_bg,latgocart_bg)
! TX is the value interpolated to the model grid box.
!-obs: convert to kg / day
    do ispc = erod, erod 
     gocart_bg_g(ispc)%src(i,j,1:nlevels_netcdf(ispc))=TX(1:nlevels_netcdf(ispc),ispc) 
    enddo
  enddo
enddo 

!print*,' MAX=',maxval(gocart_bg_g(EROD)%src(:,:,1)),maxval(gocart_bg_g(EROD)%src(:,:,1)),maxval(gocart_bg_g(SO2)%src(:,:,1))

if(ng==ngrids) then
   deallocate (RAWsrc)
   deallocate (RAWdustsrc)
endif


end subroutine read_gocart_bg
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
subroutine get_gocart_bg_indentity(spc_name,ident)
!use chem1_list
use gocart_backgr, only :  gocart_bg_nspecies=>nspecies&
                           ,gocart_bg_spc_name=>spc_name!&
!, gocart_bg_OC     => OC	  &
!, gocart_bg_BC    => BC	  &    
!, gocart_bg_SO2    => SO2	    


implicit none
integer isp
character (len=*), intent(in)  :: spc_name
integer          , intent(out) :: ident

do isp = 1,gocart_bg_nspecies
  ident=-1
  if(spc_name == gocart_bg_spc_name(isp)) then
      print*,'==>gocart_bg found for ',spc_name
      ident=isp
      return
   endif
enddo

!print*,'chem1-list specie ',trim(spc_name), ' does not match if any one of gocart_bg'
!stop 444
end subroutine get_gocart_bg_indentity
!---------------------------------------------------------------
subroutine interpol2_gocart_bg(i,j,n1,n2,rlon,rlat,ic,jc,nlat,nlon,nlev,ilatn,ilonn&
	                ,imon,nmonths,nspecies,RAWsrc,tx,longGOCART,latGOCART)
use grid_dims_out, only: grid_type
implicit none
integer n1,n2,ic,jc,nlat,nlon,i,j,imon,nmonths,nspecies,ispc,nlev
real, dimension(n1,n2) :: rlat,rlon
real, dimension(nlon,nlat,nlev,nspecies) :: RAWsrc
real ilatn,ilonn,tx(nlev,nspecies),delta
real longGOCART(nlon),latGOCART(nlat) ! NUWRF EMK
!-local var
real dlonr,dlatr,usdum
integer qi1,qi2,qj1,qj2,ncount,ii,jj
real :: x1,x2,y1,y2,denom,difflon1,difflon2,difflat1,difflat2 ! NUWRF EMK
integer :: i1,i2,j1,j2 ! NUWRF EMK

if(grid_type == 'fim') then

 ncount = 0
 dlonr = 0.
 dlatr = 0.
 do ii=1,n1-1
     ncount = ncount+1
     dlonr= dlonr + rlon(ii+1,1)-rlon(ii,1)
     dlatr= dlatr + rlat(ii+1,1)-rlat(ii,1)
 enddo
 dlonr = 0.5*dlonr/(float(ncount) + 1.E-10)
 dlatr = 0.5*dlatr/(float(ncount) + 1.E-10)
elseif(grid_type == 'mercator') then
 dlonr=0.5*(rlon(2,j)-rlon(1,j))
 dlatr=0.5*(rlat(i,n2)-rlat(i,1))/float(n2-1)
else    	     
 delta = .01*(int(100.*rlon(n1,j))-int(100.*rlon(1,j)))
 if (delta .gt. rlon(2,j)-rlon(1,j)) then            
   dlonr=0.5*(rlon(n1,j)-rlon(1,j))/float(n1-1)
 else
   dlonr=180./float(n1-1)
 endif
 dlatr=0.5*(rlat(i,n2)-rlat(i,1))/float(n2-1)
endif

qi1=int(dlonr/ilonn+0.5)
qi2=int(dlonr/ilonn+0.5)
qj1=int(dlatr/ilatn+0.5)
qj2=int(dlatr/ilatn+0.5)
 
ncount = 0
TX(:,:)  = 0.

do jj = min(max(1,jc-qj1),nlat),min(nlat,jc+qj2)
   do ii = min(max(1,ic-qi1),nlon),min(nlon,ic+qi2)   
   
     	    ncount = ncount + 1
     	    TX(1:nlev,:) = TX(1:nlev,:) + RAWsrc(ii,jj,1:nlev,:)  

   enddo
enddo
TX(1:nlev,:) = TX(1:nlev,:) / (float(ncount) + 1.E-10) ! interpolated rate

! EMK NUWRF...If WRF is at finer resolution than GOCART background field,
! use linear or bilinear interpolation.
if (ncount .lt. 2) then

   TX(:,:) = 0.

   ! Special handling of longitudes to handle interpolation near 
   ! International Date Line
   if (longGOCART(nlon) .le. rlon(i,j)) then
      i1 = nlon
      i2 = 1
   else
      do ii = 1,nlon-1
         if (longGOCART(ii  ) .le. rlon(i,j) .and. &
             longGOCART(ii+1) .gt. rlon(i,j)) then
            i1 = ii
            i2 = ii+1
            exit
         end if
      end do
   end if
   x1 = longGOCART(i1)
   x2 = longGOCART(i2)
   difflon1 = rlon(i,j) - x1
   if (difflon1 .lt. -180.) then
      difflon1 = difflon1 + 360.
   else if (difflon1 .gt. 180.) then
      difflon1 = difflon1 - 360.
   end if
   difflon2 = x2 - rlon(i,j)
   if (difflon2 .gt. 180) then
      difflon2 = difflon2 - 360.
   else if (difflon2 .lt. -180) then
      difflon2 = difflon2 + 360.
   end if
   
   ! Use linear interpolation near poles
   if (latGOCART(1)    .gt. rlat(i,j) .or. &
       latGOCART(nlat) .le. rlat(i,j)) then
      denom = (ilonn)
      if (latGOCART(1) .gt. rlat(i,j)) then
         TX(1:nlev,:) = &
              (RAWsrc(i1,1,1:nlev,:)*difflon2) + &
              (RAWsrc(i2,1,1:nlev,:)*difflon1)
      end if
      if (latGOCART(nlat) .le. rlat(i,j)) then
         TX(1:nlev,:) = &
              (RAWsrc(i1,nlat,1:nlev,:)*difflon2) + &
              (RAWsrc(i2,nlat,1:nlev,:)*difflon1)
      end if
      TX(1:nlev,:) = TX(1:nlev,:) / denom
   else ! Use bilinear interpolation away from poles
      denom = (ilonn)*(ilatn)
      do jj = 1,nlat-1
         if (latGOCART(jj  ) .le. rlat(i,j) .and. &
             latGOCART(jj+1) .gt. rlat(i,j)) then
            y1 = latGOCART(jj  )
            y2 = latGOCART(jj+1)

            TX(1:nlev,:) = &
                 (RAWsrc(i1,jj  ,1:nlev,:)*(difflon2)*(y2 - rlat(i,j))) + &
                 (RAWsrc(i2,jj  ,1:nlev,:)*(difflon1)*(y2 - rlat(i,j))) + &
                 (RAWsrc(i1,jj+1,1:nlev,:)*(difflon2)*(rlat(i,j) - y1)) + &
                 (RAWsrc(i2,jj+1,1:nlev,:)*(difflon1)*(rlat(i,j) - y1))
            TX(1:nlev,:) = TX(1:nlev,:) / denom

         end if
      end do
   end if
end if
end subroutine interpol2_gocart_bg

! NUWRF...Added subroutine CalcPress3dCenter to calculate pressure of GOCART
! background data.
!-----------------------------------------------------------------------------
!BOP
!
! !ROUTUNE: CalcPress3dCenter
!
! !INTERFACE
!
      subroutine CalcPress3dCenter (pt, hyam, hybm, psf, press3c, i2, j2, k2)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i2
      integer, intent(in) :: j2
      integer, intent(in) :: k2
      real   , intent(in) :: pt
      real   , intent(in) :: hyam(1:k2)
      real   , intent(in) :: hybm(1:k2)
      real   , intent(in) :: psf(1:i2, 1:j2) ! surface pressure field
!
! !OUTPUT PARAMETERS:
      real   , intent(out) :: press3c(1:i2, 1:j2, 1:k2+1) ! atmospheric pressure at the 
                                                          ! center of each grid box
!
! !DESCRIPTION:
! Calculates the atmospheric pressure at the center of each
! grid box assuming a hybrid sigma/pressure coordinate system.
!
!   pressure = (hyam * pt) + (hybm * psf)
!
! !LOCAL VARIABLES:
      integer :: il, ij, ik
!EOP
!-----------------------------------------------------------------------------
!BOC

      do ik = 1, k2
        do ij = 1, j2
          do il = 1, i2
              press3c(il,ij,ik) = (hyam(ik) * pt) + (hybm(ik) * psf(il,ij))
          end do
        end do
              write(300,*) ik,press3c(i2,j2,ik)
      end do

      return

      end subroutine CalcPress3dCenter

! NUWRF...Added subroutine CalcPress3dEdge to calculate pressures at edges of
! model layers in GOCART background fields.
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !ROUTUNE: CalcPress3dEdge
!
! !INTERFACE
!
      subroutine CalcPress3dEdge (pt, hyai, hybi, psf, press3e, i2, j2, k2)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i2
      integer, intent(in) :: j2
      integer, intent(in) :: k2
      real   , intent(in) :: pt
      real   , intent(in) :: hyai(1:k2+1)
      real   , intent(in) :: hybi(1:k2+1)
      real   , intent(in) :: psf(1:i2, 1:j2) ! surface pressure field
!
! !OUTPUT PARAMETERS:
      real   , intent(out) :: press3e(1:i2, 1:j2, 1:k2+1) ! atmospheric pressure at the 
                                                          ! edge of each grid box
!
! !DESCRIPTION:
! Calculates the atmospheric pressure at the edge of each
! grid box assuming a hybrid sigma/pressure coordinate system.
!
!   pressure = (hyai * pt) + (hybi * psf)
!
! !LOCAL VARIABLES:
      integer :: il, ij, ik
!EOP
!-----------------------------------------------------------------------------
!BOC

      do ik = 1, k2
        do ij = 1, j2
          do il = 1, i2
              press3e(il,ij,ik) = (hyai(ik) * pt) + (hybi(ik) * psf(il,ij))
          end do
        end do
      end do

      return

      end subroutine CalcPress3dEdge
!EOC
!-----------------------------------------------------------------------------
