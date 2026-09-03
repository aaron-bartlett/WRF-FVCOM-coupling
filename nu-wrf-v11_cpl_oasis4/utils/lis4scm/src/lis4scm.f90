!------------------------------------------------------------------------------
! NASA/GSFC, Computational and Information Science and Technology Office, 
! Code 606.
!------------------------------------------------------------------------------
!
! PROGRAM:  lis4scm
!
! AUTHOR:
! Eric Kemp, NASA CISTO/SSAI
!
! DESCRIPTION:
! Program for preparing a LIS restart file, a LDT preprocessing file for
! REAL, and a LDT parameters file for LIS for use in NU-WRF Single Column
! Mode.  LIS and LDT must have been run with a 2x2 unstaggered grid, and the
! data from teh southwest grid box are copied to the other points, resulting
! in horizontally homogeneous conditions.
!
!------------------------------------------------------------------------------

program lis4scm

   ! Imports
   use FileUnit_mod, only: PATH_MAX
   use netcdf

   ! Change defaults
   implicit none

   ! Local variables
   character(len=PATH_MAX) :: lis_restart4lis_filename     ! Input to WRF-LIS
   character(len=PATH_MAX) :: ldt_preproc4real_filename    ! Input to IDEAL
   character(len=PATH_MAX) :: ldt_parameters4lis_filename  ! Input to WRF-LIS

   ! Set filenames via namelist.lis4scm
   call read_namelist(lis_restart4lis_filename,ldt_preproc4real_filename, &
        ldt_parameters4lis_filename)

   ! Update the data files.
   call update_lis_restart4lis_file(lis_restart4lis_filename)
   call update_ldt_file(ldt_preproc4real_filename)
   call update_ldt_file(ldt_parameters4lis_filename)

contains
   
   !---------------------------------------------------------------------------

   subroutine read_namelist(lis_restart4lis_filename, &
        ldt_preproc4real_filename, &
        ldt_parameters4lis_filename)

      ! Imports
      use FileUnit_mod, only: PATH_MAX, select_file_unit, close_file_unit

      ! Change defaults
      implicit none

      ! Arguments
      character(len=PATH_MAX), intent(inout) :: lis_restart4lis_filename
      character(len=PATH_MAX), intent(inout) :: ldt_preproc4real_filename
      character(len=PATH_MAX), intent(inout) :: ldt_parameters4lis_filename

      ! Local variables
      integer :: file_unit

      ! Local namelists
      namelist /lis4scm/ lis_restart4lis_filename, ldt_preproc4real_filename, &
           ldt_parameters4lis_filename

      ! Local constants
      character(len=*), parameter :: FILENAME = 'namelist.lis4scm'

      ! Initialize namelist variables
      lis_restart4lis_filename = 'NULL'
      ldt_preproc4real_filename = 'NULL'
      ldt_parameters4lis_filename = 'NULL'

      ! Read from namelist file
      file_unit = select_file_unit()
      open(unit=file_unit,file=trim(FILENAME),delim='APOSTROPHE',status='OLD')
      read(unit=file_unit,nml=lis4scm)
      call close_file_unit(file_unit)
      
      ! Sanity checks
      if (trim(lis_restart4lis_filename) == 'NULL') then
         print*,'FATAL, lis_restart4lis_filename not read from namelist!'
         stop 1
      end if
      if (trim(ldt_preproc4real_filename) == 'NULL') then
         print*,'FATAL, ldt_preproc4real_filename not read from namelist!'
         stop 1
      end if
      if (trim(ldt_parameters4lis_filename) == 'NULL') then
         print*,'FATAL, ldt_parameters4lis_filename not read from namelist!'
         stop 1
      end if

   end subroutine read_namelist

   !---------------------------------------------------------------------------

   subroutine handle_netcdf_error(status)

      ! Change defaults
      implicit none

      ! Arguments
      integer,intent(in) :: status

      ! Check return status
      if (status .ne. NF90_NOERR) then
         print*,"FATAL, error returned from netCDF library!"
         print*,trim(nf90_strerror(status))
         stop 1
      end if

   end subroutine handle_netcdf_error

   !---------------------------------------------------------------------------

   subroutine update_lis_restart4lis_file(lis_restart4lis_filename)

      ! Change defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: lis_restart4lis_filename

      ! Local variables
      integer :: ncid
      integer :: nvars
      integer :: status
      integer :: ntiles, ntiles_id
      integer, allocatable :: varids(:)
      integer :: xtype
      integer :: ndims
      real, allocatable :: var1d(:),var2d(:,:)
      integer, allocatable :: dims(:),dimids(:)
      integer, allocatable :: start(:),count(:)
      character(len=NF90_MAX_NAME) :: varname
      real :: tmpvar
      logical :: found
      integer :: nv,i,j

      ! Open file and and check number of tiles
      status = nf90_open(trim(lis_restart4lis_filename),NF90_NOWRITE,ncid)
      call handle_netcdf_error(status)
      status = nf90_inq_dimid(ncid,"ntiles",ntiles_id)
      call handle_netcdf_error(status)
      status = nf90_inquire_dimension(ncid,ntiles_id,len=ntiles)
      call handle_netcdf_error(status)
      if (ntiles .ne. 4) then
         print*,'FATAL, ntiles is not 4, suggesting LIS is not a 2x2 grid!'
         print*,'Cannot process for NU-WRF Single Column Model'
         stop 1
      end if

      ! Close and reopen the LIS restart file, now in write mode
      status = nf90_close(ncid)
      call handle_netcdf_error(status)
      status = nf90_open(trim(lis_restart4lis_filename),NF90_WRITE,ncid)
      call handle_netcdf_error(status)

      ! Get the number of variables in this file and start looping through 
      ! them.
      status = nf90_inquire(ncid,nVariables=nvars)
      call handle_netcdf_error(status)
      allocate(varids(nvars))
      status = nf90_inq_varids(ncid,nvars,varids)
      call handle_netcdf_error(status)
      do nv = 1, nvars

         ! Check type of present variable
         status = nf90_inquire_variable(ncid,varids(nv),xtype=xtype)
         call handle_netcdf_error(status)
         if (xtype .ne. NF90_FLOAT) then
            print*,'FATAL, LIS4SCM currently only supports floats from LIS!'
            stop 1
         end if

         ! Get name of present variable
         status = nf90_inquire_variable(ncid,varids(nv),name=varname)
         call handle_netcdf_error(status)
         !print*,nv,trim(varname)
         
         ! Get dimensions
         status = nf90_inquire_variable(ncid,varids(nv),ndims=ndims)
         call handle_netcdf_error(status)
         if (ndims .ne. 1 .and. ndims .ne. 2) then
            print*,'FATAL, LIS4SCM only supports 1D and 2D variables from LIS!'
            stop 1
         end if
         allocate(dimids(ndims))
         status = nf90_inquire_variable(ncid,varids(nv),dimids=dimids)
         call handle_netcdf_error(status)

         ! See if ntiles is one of the dimensions of the array.  If not, move
         ! on to the next one.
         found = .false.
         do i = 1,ndims
            if (dimids(i) .eq. ntiles_id) then
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            deallocate(dimids)
            cycle 
         end if

         ! Get the variable dimensions
         allocate(dims(ndims))
         status = nf90_inquire_dimension(ncid,dimids(1),len=dims(1))
         call handle_netcdf_error(status)
         if (ndims .gt. 1) then
            status = nf90_inquire_dimension(ncid,dimids(2),len=dims(2))
            call handle_netcdf_error(status)
         end if

         ! Now read the array, copy the southwest tile data across the domain, 
         ! and write back to file.
         if (ndims .eq. 1) then
            allocate(var1d(dims(1)))
            status = nf90_get_var(ncid,varids(nv),values=var1d)
            call handle_netcdf_error(status)
            tmpvar = var1d(1)
            var1d(:) = tmpvar
            status = nf90_put_var(ncid,varids(nv),var1d)
            call handle_netcdf_error(status)
            deallocate(var1d)

         else if (ndims .eq. 2) then
            allocate(var2d(dims(1),dims(2)))
            allocate(start(ndims))
            allocate(count(ndims))
            start(:) = 1
            count(:) = dims(:)
            status = nf90_get_var(ncid,varids(nv),var2d,start=start, &
                 count=count)
            call handle_netcdf_error(status)
            do j = 1,dims(2)
               tmpvar = var2d(1,j)
               var2d(:,j) = tmpvar
            end do
            status = nf90_put_var(ncid,varids(nv),var2d)
            call handle_netcdf_error(status)
            deallocate(count)
            deallocate(start)
            deallocate(var2d)
         end if

         ! Clean up before next variable
         deallocate(dims)
         deallocate(dimids)

      end do

      ! Close the restart file
      status = nf90_close(ncid)
      call handle_netcdf_error(status)

      ! Clean up
      deallocate(varids)

   end subroutine update_lis_restart4lis_file

   !---------------------------------------------------------------------------

   ! Note:  Used for LDT preproc4real files and for LDT parameters4lis file.
   subroutine update_ldt_file(ldt_filename)

      ! Change defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: ldt_filename

      ! Local variables
      integer :: ncid
      integer :: nvars
      integer :: status
      integer :: nc, nc_dimid
      integer :: nr, nr_dimid
      integer :: ncb, ncb_dimid
      integer :: nrb, nrb_dimid
      integer, allocatable :: varids(:)
      integer :: xtype
      integer :: ndims
      real, allocatable :: var2d(:,:),var3d(:,:,:)
      integer, allocatable :: dims(:),dimids(:)
      integer, allocatable :: start(:),count(:)
      character(len=NF90_MAX_NAME) :: varname
      real :: tmpvar
      integer :: nv,k

      ! Open file and and check number of rows and columns
      status = nf90_open(trim(ldt_filename),NF90_NOWRITE,ncid)
      call handle_netcdf_error(status)
      status = nf90_inq_dimid(ncid,"north_south",nr_dimid)
      call handle_netcdf_error(status)
      status = nf90_inquire_dimension(ncid,nr_dimid,len=nr)
      call handle_netcdf_error(status)
      if (nr .ne. 2) then
         print*,'FATAL, number of rows is not 2, LIS is not using 2x2 grid!'
         print*,'Cannot process for NU-WRF Single Column Model'
         stop 1
      end if
      status = nf90_inq_dimid(ncid,"east_west",nc_dimid)
      call handle_netcdf_error(status)
      status = nf90_inquire_dimension(ncid,nc_dimid,len=nc)
      call handle_netcdf_error(status)
      if (nc .ne. 2) then
         print*,'FATAL, number of columns is not 2, LIS is not using 2x2 grid!'
         print*,'Cannot process for NU-WRF Single Column Model'
         stop 1
      end if

      ! See if this file also has "border" dimensions (for halos).  These
      ! exist in the LDT parameter file for LIS.  If these are absent,
      ! just record that fact.
      status = nf90_inq_dimid(ncid,"north_south_b",nrb_dimid)
      if (status .eq. NF90_EBADDIM) then
         nrb_dimid = -9999
         nrb = -9999
      else
         call handle_netcdf_error(status)
         status = nf90_inquire_dimension(ncid,nrb_dimid,len=nrb)
         call handle_netcdf_error(status)
         if (nrb .ne. 6) then
            print*,'FATAL, number of border rows is not 6'
            print*,' LIS is not using 2x2 grid!'
            print*,'Cannot process for NU-WRF Single Column Model'
            stop 1
         end if
      end if
      status = nf90_inq_dimid(ncid,"east_west_b",ncb_dimid)
      if (status .eq. NF90_EBADDIM) then
         ncb_dimid = -9999
         ncb = -9999
      else
         call handle_netcdf_error(status)
         status = nf90_inquire_dimension(ncid,ncb_dimid,len=ncb)
         call handle_netcdf_error(status)
         if (ncb .ne. 6) then
            print*,'FATAL, number of border columns is not 6'
            print*,' LIS is not using 2x2 grid!'
            print*,'Cannot process for NU-WRF Single Column Model'
            stop 1
         end if
      end if
      
      ! Close and reopen the LDT file, now in write mode
      status = nf90_close(ncid)
      call handle_netcdf_error(status)
      status = nf90_open(trim(ldt_filename),NF90_WRITE,ncid)
      call handle_netcdf_error(status)

      ! Get the number of variables in this file and start looping through 
      ! them.
      status = nf90_inquire(ncid,nVariables=nvars)
      call handle_netcdf_error(status)
      allocate(varids(nvars))
      status = nf90_inq_varids(ncid,nvars,varids)
      call handle_netcdf_error(status)
      do nv = 1, nvars

         ! Check type of present variable
         status = nf90_inquire_variable(ncid,varids(nv),xtype=xtype)
         call handle_netcdf_error(status)
         if (xtype .ne. NF90_FLOAT) then
            print*,'FATAL, LIS4SCM currently only supports floats from LDT!'
            stop 1
         end if

         ! Get name of present variable
         status = nf90_inquire_variable(ncid,varids(nv),name=varname)
         call handle_netcdf_error(status)
         !print*,nv,trim(varname)
         
         ! Get rank of variable
         status = nf90_inquire_variable(ncid,varids(nv),ndims=ndims)
         call handle_netcdf_error(status)
         if (ndims .ne. 1 .and. ndims .ne. 2 .and. ndims .ne. 3) then
            print*,&
                 'FATAL, LIS4SCM does not support ',ndims,&
                 ' variables from LDT!'
            stop 1
         end if

         ! Skip the 1-D array if it is 'time'.  Otherwise, this may be
         ! the wrong type of file -- report error and stop.
         if (ndims .eq. 1) then
            if (trim(varname) .eq. 'time') cycle
            print*,&
                 'FATAL, LIS4SCM does not support 1-D variables from LDT!'
            print*,&
                 '(Except for time variable)'
            stop 1
         end if

         ! Get the dimensions
         allocate(dimids(ndims))
         status = nf90_inquire_variable(ncid,varids(nv),dimids=dimids)
         call handle_netcdf_error(status)

         ! Check that the variable dimensions are in fact the column and
         ! row dimensions.  Better safe than sorry.
         if (dimids(1) .ne. nc_dimid .and. dimids(1) .ne. ncb_dimid) then
            print*,'FATAL, first dimension of ',trim(varname), &
                 ' is not the column dimension!'
            stop 1
         end if
         if (dimids(2) .ne. nr_dimid .and. dimids(2) .ne. nrb_dimid) then
            print*,'FATAL, second dimension of ',trim(varname), &
                 ' is not the row dimension!'
            stop 1
         end if

         ! Get the variable dimensions
         allocate(dims(ndims))
         status = nf90_inquire_dimension(ncid,dimids(1),len=dims(1))
         call handle_netcdf_error(status)
         status = nf90_inquire_dimension(ncid,dimids(2),len=dims(2))
         call handle_netcdf_error(status)
         if (ndims .eq. 3) then
            status = nf90_inquire_dimension(ncid,dimids(3),len=dims(3))
            call handle_netcdf_error(status)
         end if

         ! Now read the array, copy the southwest tile data across the domain, 
         ! and write back to file.
         allocate(start(ndims))
         allocate(count(ndims))
         start(:) = 1
         count(:) = dims(:)

         if (ndims .eq. 2) then            
            allocate(var2d(dims(1),dims(2)))
            status = nf90_get_var(ncid,varids(nv),var2d,start=start, &
                 count=count)
            call handle_netcdf_error(status)
            ! Special handling of 2D "border" fields (including halos)
            if (dims(1) .eq. ncb .and. dims(2) .eq. nrb) then
               tmpvar = var2d(3,3)
            else
               tmpvar = var2d(1,1)
            end if
            var2d(:,:) = tmpvar
            status = nf90_put_var(ncid,varids(nv),var2d)
            call handle_netcdf_error(status)
            deallocate(var2d)

         else if (ndims .eq. 3) then
            allocate(var3d(dims(1),dims(2),dims(3)))
            status = nf90_get_var(ncid,varids(nv),var3d,start=start, &
                 count=count)
            call handle_netcdf_error(status)
            do k = 1, dims(3)
               tmpvar = var3d(1,1,k)
               var3d(:,:,k) = tmpvar
            end do
            status = nf90_put_var(ncid,varids(nv),var3d)
            call handle_netcdf_error(status)
            deallocate(var3d)
         end if

         ! Clean up before next variable
         deallocate(count)
         deallocate(start)
         deallocate(dims)
         deallocate(dimids)

      end do

      ! Close the LDT file
      status = nf90_close(ncid)
      call handle_netcdf_error(status)

      ! Clean up
      deallocate(varids)

   end subroutine update_ldt_file

end program lis4scm
