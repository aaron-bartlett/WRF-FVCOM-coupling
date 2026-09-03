!------------------------------------------------------------------------------
! NASA/GSFC, Computational and Information Science and Technology Office,
! Code 606
!------------------------------------------------------------------------------
!
! MODULE: qfed_emissions
!
! AUTHOR:
! Eric Kemp, NASA CISTO/SSAI
!
! DESCRIPTION:
! Contains code to read in and interpolate NASA QFED (Quick Fire Emissions 
! Dataset) data. Based heavily on gfedv3_emissions module.
!
! REVISION:
! 14 May 2015 - Initial version
!------------------------------------------------------------------------------

module qfed_emissions

   ! Import modules
   use AeM_emission_factors, nspecies=>AeM_nspecies, N2_nitrogenio=>N2

   ! Change defaults
   implicit none

   ! Public structures
   type qfed_vars
      real, pointer, dimension(:,:,:) :: src
   end type qfed_vars
   public :: qfed_vars

   type (qfed_vars), allocatable :: qfed_g(:)
   public :: qfed_g

   ! Public routines
   public :: alloc_qfed
   public :: nullify_qfed

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  alloc_qfed
   !
   ! DESCRIPTION:  Allocates member pointers of qfed_vars structure array.
   !
   !---------------------------------------------------------------------------

   subroutine alloc_qfed(qfed,n1,n2,n3,nspecies)

      implicit none

      ! Arguments
      type(qfed_vars),dimension(nspecies),intent(inout) :: qfed
      integer,intent(in) :: n1,n2,n3

      ! Local variables
      integer :: ispc,nspecies

      do ispc=1,nspecies
         allocate(qfed(ispc)%src(n1,n2,n3))
      end do

      return
   end subroutine alloc_qfed

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  nullify_qfed
   !
   ! DESCRIPTION: Nullifies member pointers in qfed_vars structure array.
   !
   !---------------------------------------------------------------------------

   subroutine nullify_qfed(qfed,nspecies)

      implicit none

      ! Arguments
      type(qfed_vars),dimension(nspecies),intent(inout) :: qfed
      integer, intent(in) :: nspecies

      ! Local variables
      integer :: ispc

      do ispc=1,nspecies
         if (associated(qfed(ispc)%src)) nullify(qfed(ispc)%src)
      end do

      return
   end subroutine nullify_qfed

end module qfed_emissions

! FIXME: Put into qfed_emissions module
!------------------------------------------------------------------------------
!
! ROUTINE:  mem_qfed
!
! DESCRIPTION: Driver level allocator for qfed_g array.
!
!------------------------------------------------------------------------------

subroutine mem_qfed(n1,n2,n3)

   ! Import modules
   use qfed_emissions

   ! Change defaults
   implicit none

   ! Arguments
   integer,intent(in) :: n1,n2,n3

   if (.not. allocated(qfed_g)) allocate(qfed_g(nspecies))
   call nullify_qfed(qfed_g,nspecies)
   call alloc_qfed(qfed_g,n1,n2,n3,nspecies)

   return
end subroutine mem_qfed

! FIXME: Put in qfed_emissions module.
!------------------------------------------------------------------------------
!
! ROUTINE:  read_qfed
!
! DESCRIPTION: Reads and interpolates QFED emissions. Borrows heavily from 
! GFEDV3 processing code, which in turn follows GFEDV2.
!
!------------------------------------------------------------------------------

subroutine read_qfed(iyear,imon,iday,qfed_names,ng,ngrids,n1,n2,n3, &
     rlat,rlon,rland,deltax,deltay,xt,yt,plat,plon)

   ! Import modules
   use grid_dims_out, only: maxfiles,grid_type,qfed_data_dir, &
        use_bbem_plumerise
   use qfed_emissions
   use bbbem_plumerise

   ! Change defaults
   implicit none

   ! Arguments
   integer, intent(in) :: iyear,imon,iday,ng,ngrids,n1,n2,n3
   character(len=*),intent(in) :: qfed_names ! "comma" delimited species names
   real, intent(in), dimension(n1,n2) :: rlat,rlon,rland
   real, intent(in) :: deltax,deltay
   real, intent(inout) :: xt(n1),yt(n1),plat,plon

   ! Local variables
   integer :: numConsideredSpecies,iConsideredSpecies
   character(len=20),allocatable :: consideredSpecies(:)
   real,allocatable,save,dimension(:,:) :: carb,veg
   real,allocatable,save,dimension(:,:,:) :: drybbc
   real,allocatable,save,dimension(:) :: longqfed,latqfed
   integer,save :: nlat,nlon
   real, save :: typical_fire_size(0:3)
   integer,allocatable,save,dimension(:,:) :: veg_type
   real, parameter :: UNIT_CONVERSION = 1000.*86400.   ! kg/s/m-2 to g/m-2/day
   integer :: ident, ident2
   integer :: iveg_qg,ispc,icount
   character(len=240) :: qfedFilename
   integer :: i1,i2,j1,j2,ic,jc
   integer :: i,j
   real :: TX(nspecies)
   real :: ilatn,ilonn,iveg_ag
   real :: dlat1,dlat2,dlon1,dlon2

   integer :: ierr
   integer :: start(3),count(3)
   integer :: ncid,varid,irank,dimids(1)
   include "netcdf.inc"

   ! Process list of QFED species to process. Uses same subroutines here as
   ! for GFED3
   call getNumConsideredSpecies(qfed_names,numConsideredSpecies)
   allocate(consideredSpecies(numConsideredSpecies))
   call getConsideredSpecies(qfed_names,numConsideredSpecies,consideredSpecies)

   ! Allocate memory for plumerise arrays
   if (use_bbem_plumerise == 1) call mem_bbbem_plume(n1,n2,n3)

   !--
   if (ng == 1) then

      print*,'*** NASA QFED Emissions ***'

      if (.not. allocated(carb)) then

         ! Get lat/lon information. Use data from first QFED species file
         write(qfedFileName, &
              '(A,A,I4.4,A,I2.2,A,A,A,I3.3,A,I4.4,I2.2,I2.2,A)') &
              trim(qfed_data_dir),'/Y',iyear,'/M',imon,'/qfed2.emis_', &
              trim(consideredSpecies(1)),'.',5,'.',iyear,imon,iday,'.nc4'

         ! Note: longqfed and latqfed are allocated in these subroutine calls,
         ! based on dimensions read in from netCDF file.
         call get1drealdims(trim(qfedFileName),'lon',nlon)
         allocate(longqfed(nlon))
         call read1dreal(trim(qfedFileName),'lon',nlon,longqfed)
         call get1drealdims(trim(qfedFileName),'lat',nlat)
         allocate(latqfed(nlat))
         call read1dreal(trim(qfedFileName),'lat',nlat,latqfed)

         allocate(carb(nlon,nlat))
         allocate(drybbc(nspecies,nlon,nlat))
         allocate(veg(nlon,nlat))
         allocate(veg_type(nlon,nlat))

      end if

      carb = 0.0
      drybbc = 0.0
      veg = 0.0   
      veg_type = 0

      do iConsideredSpecies = 1, numConsideredSpecies

         ! Match QFED species to AeM list
         call get_qfed_identity(trim(consideredSpecies(iConsideredSpecies)), &
              ident,ident2)

         if (ident < 0) cycle ! No match

         ! Read daily biomass burning emissions
         write(qfedFileName, &
              '(A,A,I4.4,A,I2.2,A,A,A,I3.3,A,I4.4,I2.2,I2.2,A)') &
              trim(qfed_data_dir),'/Y',iyear,'/M',imon,'/qfed2.emis_', &
              trim(consideredSpecies(iConsideredSpecies)),'.',5,'.', &
              iyear,imon,iday,'.nc4'

         print*,'Reading from ',trim(qfedFileName)
         call read2dreal(qfedFileName,'biomass',nlon,nlat,carb)

         ! Copy to drybbc array
         if (ident > 0 .and. ident2 > 0) then
            ! 50/50 split, change to g/m2/day
            drybbc(ident,:,:) = 0.5*carb*UNIT_CONVERSION  
            drybbc(ident2,:,:) = 0.5*carb*UNIT_CONVERSION 
         else if (ident > 0) then
            drybbc(ident,:,:) = carb*UNIT_CONVERSION ! Change to g/m2/day
         end if
      end do

      ! Set veg_type
      where(nint(veg(:,:)) == 1) veg_type(:,:)=3
      where(nint(veg(:,:)) == 2) veg_type(:,:)=1
      where(nint(veg(:,:)) == 3) veg_type(:,:)=2

      !- for plumerise - typical fire size
      typical_fire_size(0)=0.
      typical_fire_size(1)=20.* 1.e4 ! 20ha for tropical forest
      typical_fire_size(2)=20.* 1.e4 ! 20ha for extra tropical
      typical_fire_size(3)=5. * 1.e4 ! 5ha for savana/past ...

      call yrevert(nlon,nlat,1,veg_type)

   end if ! (ng == 1)

   ! Make sure this is calculated for every grid
   ilonn = longqfed(2) - longqfed(1)
   ilatn = latqfed(2) - latqfed(1)

   ! Interpolate to model grid box
   do i=1,n1
      do j=1,n2

         call get_index1(i,j,nlon,nlat,n1,n2,rlat,rlon,longqfed,latqfed &
              ,ilatn, ilonn ,dlat1,dlat2,dlon1,dlon2,i1,j1,i2,j2,ic,jc)
         if(ic == 0 .or. jc == 0. .or. ic > nlon .or. jc > nlat) cycle
         ! Use same interpolation routine as for GFEDV3.
         call interpol_gfedv3(i,j,n1,n2,rlon,rlat,ic,jc,nlat,nlon, &
              ilatn,ilonn,nspecies,drybbc,veg_type,tx)

         ! TX is the value interpolated to the model grid box.
         ! Convert from g/m^2/day to kg/m^2/day
         do ispc=1,nspecies
            qfed_g(ispc)%src(i,j,1)=TX(ispc) * 1.e-3 
         enddo

         !------- plumerise section - only for areas with burning
         if(use_bbem_plumerise /= 1 .or. veg_type(ic,jc) == 0 &
              .or. qfed_g(1)%src(i,j,1) < 1.e-6) cycle

         iveg_ag = veg_type(ic,jc)

         do ispc=1,nspecies
            bbbem_plume_g(ispc,iveg_ag)%src(i,j,1) = &
                 qfed_g(ispc)%src(i,j,1)*flaming(iveg_ag)
         end do

         !-fires properties : accumulated size
         bbbem_plume_fire_prop_g(qarea_agreg,iveg_ag)%fire_prop(i,j) = &
              typical_fire_size(veg_type(ic,jc))

         !-fires properties : fire number
         bbbem_plume_fire_prop_g(qfires_agreg,iveg_ag)%fire_prop(i,j) = 1.     
         !------- end of plumerise section

      end do
   end do

   !------- plumerise section 2
   if(use_bbem_plumerise == 1) then
      !- produce some statistical quantities for the plumerise model
      !
      do i=1,n1
         do j=1,n2
            do iveg_ag=1,nveg_agreg
               ! calculate the mean size of aggregate fire

               if( bbbem_plume_fire_prop_g(qfires_agreg,iveg_ag) % &
                    fire_prop(i,j) .ge. 1.) then 

                  bbbem_plume_fire_prop_g(qarea_agreg ,iveg_ag) % &
                       fire_prop(i,j) = & ! mean  size =
                       bbbem_plume_fire_prop_g(qarea_agreg ,iveg_ag) % &
                       fire_prop(i,j) / & ! total size / nfires
                       bbbem_plume_fire_prop_g(qfires_agreg,iveg_ag) % &
                       fire_prop(i,j)
               else   
                  bbbem_plume_fire_prop_g(qarea_agreg ,iveg_ag) % &
                       fire_prop(i,j) = 0.
               endif

               !- convert from mass consumed during the flaming phase to 
               !  fraction of the total mass consumed
               do ispc=1,nspecies
                  if(qfed_g(ispc)%src(i,j,1) > 0.) &    
                                !- fraction consumed at the flaming phase:
                       bbbem_plume_g(ispc,iveg_ag)%src(i,j,1) = &
                       bbbem_plume_g(ispc,iveg_ag)%src(i,j,1) &
                       / qfed_g(ispc)%src(i,j,1)
               enddo
            enddo
         enddo
      enddo

      !- calculate the mean value for bbbem_plume_g (to save memory)
      do i=1,n1
         do j=1,n2
            do iveg_ag=1,nveg_agreg
               icount = 0
               bbbem_plume_mean_g(iveg_ag)%src(i,j,1)=0.

               do ispc=1,nspecies
                  if(qfed_g(ispc)%src(i,j,1) > 0.) then
                     icount = icount + 1
                     bbbem_plume_mean_g(iveg_ag)%src(i,j,1) = &
                          bbbem_plume_mean_g(  iveg_ag)%src(i,j,1)+ &
                          bbbem_plume_g(ispc,iveg_ag)%src(i,j,1)
                  endif
               enddo
               bbbem_plume_mean_g(iveg_ag)%src(i,j,1) = &
                    bbbem_plume_mean_g(iveg_ag)%src(i,j,1) / &
                    (float(icount)+1.e-10)
            enddo
         enddo
      enddo

   endif
   !------- end of plumerise section 2

   if(grid_type == 'rams' .or. grid_type == 'polar') then 
      ! only for 'rams' until rland is also defined for others grids 
      do ispc=1,nspecies 
         call apply_land_restriction(n1,n2,rland,qfed_g(ispc)%src(:,:,1))
      enddo
   endif

   ! Clean up
   if(ng==ngrids) deallocate (carb,drybbc,veg,veg_type,longqfed,latqfed)
   deallocate(consideredSpecies)

   return
end subroutine read_qfed

!FIXME: Put this in qfed_emissions module.
!----------------------------------------------------------------------------
!
! ROUTINE:  get1drealdims
!
! DESCRIPTION:  Reads and returns dimension length for requested 1-D array
! from QFED netCDF file.
!
!----------------------------------------------------------------------------

subroutine get1drealdims(filename,varname,nx)
   implicit none

   ! Arguments
   character(len=*),intent(in) :: filename
   character(len=*),intent(in) :: varname
   integer,intent(inout) :: nx

   ! Local variables
   integer :: ierr,irank
   integer :: ncid,varid,dimids(1)
   include 'netcdf.inc'

   ierr = nf_open(trim(filename),NF_NOWRITE,ncid)
   if (ierr .ne. NF_NOERR) then
      print*, 'ERROR from nf_open while trying to open '//trim(filename)
      print*, trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_inq_varid(ncid,trim(varname),varid)
   if (ierr .ne. NF_NOERR) then
      print*, 'ERROR from nf_inq_varid while searching for '//trim(varname)
      print*, trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_inq_varndims(ncid,varid, irank)
   if (ierr .ne. NF_NOERR) then
      print*, 'ERROR from nf_inq_varndims while finding rank of '// &
           trim(varname)
      print*, trim(nf_strerror(ierr))
      stop 1     
   end if
   if (irank .ne. 1) then
      print*,'ERROR, rank of '//trim(varname)//' is not 1, but is ',irank
      stop 1
   end if

   ierr = nf_inq_vardimid(ncid,varid,dimids)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_inq_vardimid while finding dim ids for '// &
           trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_inq_dimlen(ncid,dimids,nx)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_inq_dimlen while finding length of '//trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_close(ncid)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_close while trying to close '//trim(filename)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

end subroutine get1drealdims

!FIXME: Put this in qfed_emissions module.
!----------------------------------------------------------------------------
!
! ROUTINE:  read1dreal
!
! DESCRIPTION:  Reads and returns 1-d real array from QFED netCDF file.
! Will check to ensure assumed dimension matches that in file.
!
!----------------------------------------------------------------------------

subroutine read1dreal(filename,varname,nx,array1d)
   implicit none

   ! Arguments
   character(len=*),intent(in) :: filename
   character(len=*),intent(in) :: varname
   integer,intent(in) :: nx
   real,intent(inout) :: array1d(nx)

   ! Local variables
   integer :: ierr,irank
   integer :: ncid,varid,dimids(1)
   integer :: tmp_nx
   integer :: start(1),count(1)
   include 'netcdf.inc'

   ierr = nf_open(trim(filename),NF_NOWRITE,ncid)
   if (ierr .ne. NF_NOERR) then
      print*, 'ERROR from nf_open while trying to open '//trim(filename)
      print*, trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_inq_varid(ncid,trim(varname),varid)
   if (ierr .ne. NF_NOERR) then
      print*, 'ERROR from nf_inq_varid while searching for '//trim(varname)
      print*, trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_inq_varndims(ncid,varid, irank)
   if (ierr .ne. NF_NOERR) then
      print*, 'ERROR from nf_inq_varndims while finding rank of '// &
           trim(varname)
      print*, trim(nf_strerror(ierr))
      stop 1     
   end if
   if (irank .ne. 1) then
      print*,'ERROR, rank of '//trim(varname)//' is not 1, but is ',irank
      stop 1
   end if

   ierr = nf_inq_vardimid(ncid,varid,dimids)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_inq_vardimid while finding dim ids for '// &
           trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_inq_dimlen(ncid,dimids,tmp_nx)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_inq_dimlen while finding length of '//trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

   if (nx .ne. tmp_nx) then
      print*,'ERROR, dimension mismatch!'
      print*,'Expected ',nx,' found ',tmp_nx
      stop 1
   end if

   count(1) = nx
   start(1) = 1

   ierr = nf_get_vara_real(ncid,varid,start,count,array1d)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_get_vara_real while trying to read '//trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_close(ncid)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_close while trying to close '//trim(filename)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

   return
end subroutine read1dreal

!FIXME: Put this in qfed_emissions module.
!----------------------------------------------------------------------------
!
! ROUTINE:  read2dreal
!
! DESCRIPTION:  Reads and returns 2-d real array from QFED netCDF file.
! Will check to ensure assumed dimensions matches that in file.
!
!----------------------------------------------------------------------------

subroutine read2dreal(filename,varname,nx,ny,array2d)
   implicit none

   ! Arguments
   character(len=*),intent(in) :: filename
   character(len=*),intent(in) :: varname
   integer,intent(in) :: nx,ny
   real,intent(inout) :: array2d(nx,ny)

   ! Local variables
   integer :: ntime
   integer :: ierr,irank
   real, allocatable :: tmp_3d(:,:,:)
   integer :: tmp_nx,tmp_ny
   integer :: ncid,varid,dimids(3)
   integer :: start(3),count(3)
   include 'netcdf.inc'

   ierr = nf_open(trim(filename),NF_NOWRITE,ncid)
   if (ierr .ne. NF_NOERR) then
      print*, 'ERROR from nf_open while trying to open '//trim(filename)
      print*, trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_inq_varid(ncid,trim(varname),varid)
   if (ierr .ne. NF_NOERR) then
      print*, 'ERROR from nf_inq_varid while searching for '//trim(varname)
      print*, trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_inq_varndims(ncid,varid, irank)
   if (ierr .ne. NF_NOERR) then
      print*, 'ERROR from nf_inq_varndims while finding rank of '// &
           trim(varname)
      print*, trim(nf_strerror(ierr))
      stop 1     
   end if
   if (irank .ne. 3) then
      print*,'ERROR, rank of '//trim(varname)//' is not 3, but is ',irank
      stop 1
   end if

   ierr = nf_inq_vardimid(ncid,varid,dimids)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_inq_vardimid while finding dim ids for '// &
           trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

   ierr = nf_inq_dimlen(ncid,dimids(1),tmp_nx)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_inq_dimlen while finding dimensions of '// &
           trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if
   ierr = nf_inq_dimlen(ncid,dimids(2),tmp_ny)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_inq_dimlen while finding dimensions of '// &
           trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if
   if (nx .ne. tmp_nx .or. ny .ne. tmp_ny) then
      print*,'ERROR, dimension mismatch for '//trim(varname)
      print*,'Expected ',nx,' ',ny
      print*,'Found ',nx,' ',ny
      stop 1
   end if

   ierr = nf_inq_dimlen(ncid,dimids(3),ntime)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_inq_dimlen while finding dimensions of '// &
           trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if
   if (ntime .ne. 1) then
      print*,'ERROR, expected single time level in '//trim(filename)
      stop 1
   end if

   allocate(tmp_3d(nx,ny,ntime))
   count(1) = nx ; count(2) = ny ; count(3) = 1
   start(:) = 1

   ierr = nf_get_vara_real(ncid,varid,start,count,tmp_3d)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_get_vara_real while trying to read '//trim(varname)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

   array2d(:,:) = tmp_3d(:,:,1)
   deallocate(tmp_3d)

   ierr = nf_close(ncid)
   if (ierr .ne. NF_NOERR) then
      print*,'ERROR from nf_close while trying to close '//trim(filename)
      print*,trim(nf_strerror(ierr))
      stop 1
   end if

   return
end subroutine read2dreal

!FIXME: Put this in qfed_emissions module.
!----------------------------------------------------------------------------
!
! ROUTINE:  get_qfed_identity
!
! DESCRIPTION: Matches ID of QFED variable with that from AeM emissions list.
! It is possible that two IDs will be returned if a QFED species is to be
! partitioned.
!
!----------------------------------------------------------------------------

subroutine get_qfed_identity(spc_name,ident,ident2)
   use qfed_emissions, only : qfed_nspecies=>nspecies &
        ,qfed_spc_name=>AeM_spc_name
   implicit none

   ! Arguments
   character(len=*),intent(in) :: spc_name
   integer,intent(out) :: ident, ident2

   ! Local variables
   integer :: isp

   ! Match species with that specified in AeM list
   ident = -1
   ident2 = -1

   ! First, look for exact match.
   do isp = 1,qfed_nspecies
      if (spc_name == qfed_spc_name(isp)) then
         print*,'Matched with AeM species ',trim(spc_name)
         ident = isp
         return
      end if
   end do

   ! Next, look for equivalents.
   if (trim(spc_name) == 'acet') then
      ident = special_identity(trim(spc_name),'Acetone',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'ald2') then
      ident = special_identity(trim(spc_name),'Acetald',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'alk4') then
      ident = special_identity(trim(spc_name),'n_butane',qfed_spc_name, &
           qfed_nspecies)
      ident2 = special_identity(trim(spc_name),'i-butane',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'bc') then
      ident = special_identity(trim(spc_name),'BC',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'c2h6') then
      ident = special_identity(trim(spc_name),'C2H6',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'c3h6') then
      ident = special_identity(trim(spc_name),'C3H6',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'c3h8') then
      ident = special_identity(trim(spc_name),'C3H8',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'ch2o') then
      ident = special_identity(trim(spc_name),'Formaldehyde',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'ch4') then
      ident = special_identity(trim(spc_name),'CH4',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'co') then
      ident = special_identity(trim(spc_name),'CO',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'co2') then
      ident = special_identity(trim(spc_name),'CO2',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'mek') then
      ident = special_identity(trim(spc_name),'2_Butanone',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'nh3') then
      ident = special_identity(trim(spc_name),'NH3',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'no') then
      ident = special_identity(trim(spc_name),'NOx',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'oc') then
      ident = special_identity(trim(spc_name),'OC',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'pm25') then
      ident = special_identity(trim(spc_name),'BBURN2',qfed_spc_name, &
           qfed_nspecies)
   else if (trim(spc_name) == 'so2') then
      ident = special_identity(trim(spc_name),'SO2',qfed_spc_name, &
           qfed_nspecies)
   end if
   !  ! It's possible no match exists.
   !  if (ident < 0) then
   !     print*,'Warning, cannot process QFED species ',trim(spc_name)
   !  end if

   return
contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE: special_identity
   !
   ! DESCRIPTION: Internal function for special comparisons with AeM species 
   ! list. Assumes a match will be found with AeM. A fatal error will occur
   ! if this is not the case -- useful for detecting future changes to the    
   ! AeM species list.
   !
   !---------------------------------------------------------------------------

   function special_identity(orig_spc_name,spc_name,qfed_spc_name, &
        qfed_nspecies) result(ident)

      implicit none

      ! Arguments
      character(len=*),intent(in) :: orig_spc_name
      character(len=*),intent(in) :: spc_name
      integer,intent(in) :: qfed_nspecies
      character(len=*),intent(in) :: qfed_spc_name(qfed_nspecies)

      ! Return variable
      integer :: ident

      ! Local variables
      integer :: ii

      ident= -1
      do ii = 1,qfed_nspecies
         if (trim(spc_name) == trim(qfed_spc_name(ii))) then
            print*,'Matched QFED species ',trim(orig_spc_name), &
                 ' with AeM species ',trim(qfed_spc_name(ii))
            ident=ii
            return               
         end if
      end do
      print*,'INTERNAL ERROR, cannot match QFED species ', &
           trim(orig_spc_name),' with AeM species!'
      stop 1
   end function special_identity

end subroutine get_qfed_identity

