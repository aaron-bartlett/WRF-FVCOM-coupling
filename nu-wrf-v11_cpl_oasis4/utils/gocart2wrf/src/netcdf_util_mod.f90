!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration and Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  netcdf_util_mod
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Utility subroutines for calling the netCDF library.
!
!------------------------------------------------------------------------------

module netcdf_util_mod
 
   ! Reset defaults
   implicit none
   private

   ! Include netCDF constants and routines.  Keep private.
   include "netcdf.inc"   

   ! Public subroutines
   public :: handle_netcdf_error
   public :: open_netcdf_readfile
   public :: open_netcdf_writefile
   public :: close_netcdf_file
   public :: close_netcdf_file_if_open
   public :: set_netcdf_define_mode
   public :: unset_netcdf_define_mode
   public :: read_netcdf_dim_id
   public :: read_netcdf_dim_len
   public :: define_netcdf_real_array_4d
   public :: write_netcdf_integer_attribute
   public :: write_netcdf_text_attribute
   public :: write_netcdf_real_array_4d
   public :: read_netcdf_var_id
   public :: read_netcdf_real_array_4d
   public :: read_netcdf_real_array_3d
   public :: read_netcdf_real_array_2d
   public :: read_netcdf_real_array_1d
   public :: read_netcdf_double_array_1d
   public :: read_netcdf_character_array_1d
   public :: read_netcdf_dimension
   public :: write_netcdf_dimension
   public :: open_netcdf_newfile
   public :: check_define_netcdf_real_array_4d
   public :: get_real_array_4d_netcdf4_settings
   public :: define_real_array_4d_netcdf4_settings
   public :: read_netcdf_real_variable
contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  handle_netcdf_error
   !
   ! DESCRIPTION:  Checks return status of netCDF library call.  If error
   ! indicated, print description of error and stop program.
   !
   !---------------------------------------------------------------------------

   subroutine handle_netcdf_error(status)

      ! Arguments
      integer,intent(in) :: status

      ! Check return status
      if (status .ne. NF_NOERR) then
         print*,"ERROR returned from netCDF library!!!"
         print*,trim(nf_strerror(status))
         stop
      end if

   end subroutine handle_netcdf_error

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  open_netcdf_readfile
   !
   ! DESCRIPTION:  Opens a netCDF file in read mode, returns file ID.
   !
   !---------------------------------------------------------------------------

   function open_netcdf_readfile(filename) result (ncid)
 
      ! Arguments
      character(len=*),intent(in) :: filename

      ! Return value
      integer :: ncid

      ! Local variables
      integer :: status

      ! Open the file
      status = nf_open(trim(filename),NF_NOWRITE,ncid)
      call handle_netcdf_error(status)

      return
   end function open_netcdf_readfile

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  open_netcdf_readfile
   !
   ! DESCRIPTION:  Opens a netCDF file in write mode, returns file ID.
   !
   !---------------------------------------------------------------------------

   function open_netcdf_writefile(filename) result (ncid)
 
      ! Arguments
      character(len=*),intent(in) :: filename

      ! Return value
      integer :: ncid

      ! Local variables
      integer :: status

      ! Open the file
      status = nf_open(trim(filename),NF_WRITE,ncid)
      call handle_netcdf_error(status)

      return
   end function open_netcdf_writefile

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  close_netcdf_file
   !
   ! DESCRIPTION:  Closes a netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine close_netcdf_file(ncid)
      
      ! Arguments
      integer,intent(in) :: ncid

      ! Local variables
      integer :: status

      ! Close the file
      status = nf_close(ncid)
      call handle_netcdf_error(status)
      
   end subroutine close_netcdf_file

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  close_netcdf_file
   !
   ! DESCRIPTION:  Closes a netCDF file.  This version ignores an error code
   ! associated with the file not being open--basically, nothing will happen.
   !
   !---------------------------------------------------------------------------

   subroutine close_netcdf_file_if_open(ncid)

      ! Arguments
      integer,intent(in) :: ncid

      ! Local variables
      integer :: status

      ! Close the file.  Ignore status message indicating file was never
      ! open.
      status = nf_close(ncid)
      if (status .ne. NF_NOERR .and. status .ne. NF_EBADID) then
         call handle_netcdf_error(status)
      end if

   end subroutine close_netcdf_file_if_open

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  set_netcdf_define_mode
   !
   ! DESCRIPTION:  Puts an open netCDF file in define mode.
   !
   !---------------------------------------------------------------------------

   subroutine set_netcdf_define_mode(ncid)
      
      ! Arguments
      integer, intent(in) :: ncid

      ! Local variables
      integer :: status

      status = nf_redef(ncid)
      call handle_netcdf_error(status)

   end subroutine set_netcdf_define_mode

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  unset_netcdf_define_mode
   !
   ! DESCRIPTION:  Turns off define mode in an open netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine unset_netcdf_define_mode(ncid)
      
      ! Arguments
      integer, intent(in) :: ncid

      ! Local variables
      integer :: status

      status = nf_enddef(ncid)
      call handle_netcdf_error(status)

   end subroutine unset_netcdf_define_mode

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_dim_id
   !
   ! DESCRIPTION:  Reads in and returns ID of requested dimension from netCDF
   ! file.
   !
   !---------------------------------------------------------------------------

   function read_netcdf_dim_id(ncid,dimension_name) result (dim_id)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: dimension_name

      ! Return variable
      integer :: dim_id

      ! Local variables
      integer :: status

      status = nf_inq_dimid(ncid,trim(dimension_name),dim_id)
      call handle_netcdf_error(status)

      return
   end function read_netcdf_dim_id

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_dim_len
   !
   ! DESCRIPTION:  Reads in and returns length of requested dimension from 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   function read_netcdf_dim_len(ncid,dim_id) result (dim_len)

      ! Arguments
      integer,intent(in) :: ncid
      integer,intent(in) :: dim_id

      ! Return variable
      integer :: dim_len

      ! Local variables
      integer :: status

      status = nf_inq_dimlen(ncid,dim_id,dim_len)
      call handle_netcdf_error(status)

      return
   end function read_netcdf_dim_len

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  define_netcdf_real_array_4d
   !
   ! DESCRIPTION:  Defines a 4D real array in a netCDF file.  Values are
   ! written separately.
   !
   !---------------------------------------------------------------------------

   function define_netcdf_real_array_4d(ncid,variable_name, &
        dim1_name,dim2_name,dim3_name,dim4_name) result (var_id)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: variable_name
      character(len=*),intent(in) :: dim1_name
      character(len=*),intent(in) :: dim2_name
      character(len=*),intent(in) :: dim3_name
      character(len=*),intent(in) :: dim4_name

      ! Return variable
      integer :: var_id

      ! Local variables
      integer :: dim_ids(4)
      integer :: status

      ! Get dimension IDs
      dim_ids(1) = read_netcdf_dim_id(ncid,trim(dim1_name))
      dim_ids(2) = read_netcdf_dim_id(ncid,trim(dim2_name))
      dim_ids(3) = read_netcdf_dim_id(ncid,trim(dim3_name))
      dim_ids(4) = read_netcdf_dim_id(ncid,trim(dim4_name))

      ! Now define the variable
      status = nf_def_var(ncid,trim(variable_name),NF_REAL,4,dim_ids,var_id)
      if (status .ne. NF_NOERR .and. status .ne. NF_EINDEFINE) then
         call handle_netcdf_error(status)
      end if

      return
   end function define_netcdf_real_array_4d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  write_netcdf_integer_attribute
   !
   ! DESCRIPTION:  Writes an integer attribute to a netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine write_netcdf_integer_attribute(ncid,var_id,attr_name, &
      length, ivals) 
      
      ! Arguments
      integer,intent(in) :: ncid
      integer,intent(in) :: var_id
      character(len=*),intent(in) :: attr_name
      integer,intent(in) :: length
      integer,intent(in) :: ivals(length)

      ! Local variables
      integer :: status

      ! Enter define mode, add the attribute, and then close define mode.
      status = nf_put_att_int(ncid,var_id,trim(attr_name),NF_INT,length,ivals)
      call handle_netcdf_error(status)

   end subroutine write_netcdf_integer_attribute

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  write_netcdf_text_attribute
   !
   ! DESCRIPTION:  Writes a text attribute to a netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine write_netcdf_text_attribute(ncid,var_id,attr_name,length,text) 
      
      ! Arguments
      integer,intent(in) :: ncid
      integer,intent(in) :: var_id
      character(len=*),intent(in) :: attr_name
      integer,intent(in) :: length
      character(len=*),intent(in) :: text

      ! Local variables
      integer :: status
      
      ! Add the attribute
      status = nf_put_att_text(ncid,var_id,trim(attr_name),length, &
           text(1:length))
      call handle_netcdf_error(status)

   end subroutine write_netcdf_text_attribute

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  write_netcdf_real_array_4d
   !
   ! DESCRIPTION:  Writes a 4D real array to a netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine write_netcdf_real_array_4d(ncid,var_id,dim1,dim2,dim3,dim4,array)

      ! Arguments
      integer,intent(in) :: ncid
      integer,intent(in) :: var_id
      integer,intent(in) :: dim1
      integer,intent(in) :: dim2
      integer,intent(in) :: dim3
      integer,intent(in) :: dim4
      real,intent(in) :: array(dim1,dim2,dim3,dim4)

      ! Local variables
      integer :: start(4)
      integer :: count(4)
      integer :: status

      ! Write the array
      start = (/1,1,1,1/)
      count = (/dim1,dim2,dim3,dim4/)
      status = nf_put_vara_real(ncid,var_id,start,count,array)
      call handle_netcdf_error(status)

   end subroutine write_netcdf_real_array_4d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_var_id
   !
   ! DESCRIPTION:  Retrieves ID for a netCDF variable.
   !
   !---------------------------------------------------------------------------

   function read_netcdf_var_id(ncid,varname) result(var_id)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: varname

      ! Return variable
      integer :: var_id

      ! Local variables
      integer :: status

      ! Get netCDF ID of array
      status = nf_inq_varid(ncid,trim(varname),var_id)
      call handle_netcdf_error(status)

      return
   end function read_netcdf_var_id

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_real_array_4d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 4D real array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine read_netcdf_real_array_4d(ncid,varname,dim1,dim2,dim3,dim4, &
        array, startarg, countarg)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: varname
      integer,intent(in) :: dim1
      integer,intent(in) :: dim2
      integer,intent(in) :: dim3
      integer,intent(in) :: dim4
      real, dimension(dim1,dim2,dim3,dim4), intent(out) :: array
      integer,intent(in),optional :: startarg(4)
      integer,intent(in),optional :: countarg(4)

      ! Local variables
      integer :: start(4)
      integer :: count(4)
      integer :: varid
      integer :: status

      ! Get netCDF ID of array
      status = nf_inq_varid(ncid,trim(varname),varid)
      call handle_netcdf_error(status)

      ! Now read the array
      if (present(startarg)) then
         start = startarg
      else
         start = (/ 1, 1, 1, 1 /)
      end if
      if (present(countarg)) then
         count = countarg
      else
         count = (/ dim1, dim2, dim3, dim4 /)
      end if
      status = nf_get_vara_real(ncid, varid, start, count, array)
      call handle_netcdf_error(status)

   end subroutine read_netcdf_real_array_4d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_real_array_3d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 3D real array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine read_netcdf_real_array_3d(ncid,varname,dim1,dim2,dim3,array, &
        startarg, countarg)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: varname
      integer,intent(in) :: dim1
      integer,intent(in) :: dim2
      integer,intent(in) :: dim3
      real, dimension(dim1,dim2,dim3), intent(out) :: array
      integer,intent(in),optional :: startarg(3)
      integer,intent(in),optional :: countarg(3)

      ! Local variables
      integer :: start(3)
      integer :: count(3)
      integer :: varid
      integer :: status

      ! Get netCDF ID of array
      status = nf_inq_varid(ncid,trim(varname),varid)
      call handle_netcdf_error(status)

      ! Now read the array
      if (present(startarg)) then
         start = startarg(:)
      else
         start = (/ 1, 1, 1 /)         
      end if
      if (present(countarg)) then
         count = countarg(:)
      else
         count = (/ dim1, dim2, dim3 /)
      end if
      status = nf_get_vara_real(ncid, varid, start, count, array)
      call handle_netcdf_error(status)

   end subroutine read_netcdf_real_array_3d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_real_array_2d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 2D real array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine read_netcdf_real_array_2d(ncid,varname,dim1,dim2,array)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: varname
      integer,intent(in) :: dim1
      integer,intent(in) :: dim2
      real, dimension(dim1,dim2), intent(out) :: array

      ! Local variables
      integer :: start(2)
      integer :: count(2)
      integer :: varid
      integer :: status

      ! Get netCDF ID of array
      status = nf_inq_varid(ncid,trim(varname),varid)
      call handle_netcdf_error(status)

      ! Now read the array
      start = (/ 1, 1 /)
      count = (/ dim1, dim2 /)
      status = nf_get_vara_real(ncid, varid, start, count, array)
      call handle_netcdf_error(status)

   end subroutine read_netcdf_real_array_2d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_real_array_1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D real array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine read_netcdf_real_array_1d(ncid,varname,dim1,array)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: varname
      integer,intent(in) :: dim1
      real, dimension(dim1), intent(out) :: array

      ! Local variables
      integer :: start(1)
      integer :: count(1)
      integer :: varid
      integer :: status

      ! Get netCDF ID of array
      status = nf_inq_varid(ncid,trim(varname),varid)
      call handle_netcdf_error(status)

      ! Now read the array
      start = (/ 1 /)
      count = (/ dim1 /)
      status = nf_get_vara_real(ncid, varid, start, count, array)
      call handle_netcdf_error(status)

   end subroutine read_netcdf_real_array_1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_double_array_1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D double array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine read_netcdf_double_array_1d(ncid,varname,dim1,array)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: varname
      integer,intent(in) :: dim1
      double precision, dimension(dim1), intent(out) :: array

      ! Local variables
      integer :: start(1)
      integer :: count(1)
      integer :: varid
      integer :: status

      ! Get netCDF ID of array
      status = nf_inq_varid(ncid,trim(varname),varid)
      call handle_netcdf_error(status)

      ! Now read the array
      start = (/ 1 /)
      count = (/ dim1 /)
      status = nf_get_vara_double(ncid, varid, start, count, array)
      call handle_netcdf_error(status)

   end subroutine read_netcdf_double_array_1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_character_array_1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D character array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine read_netcdf_character_array_1d(ncid,varname,charlen,dim1,array)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: varname
      integer,intent(in) :: charlen ! Number of characters
      integer,intent(in) :: dim1
      character(len=*), dimension(dim1), intent(out) :: array

      ! Local variables
      integer :: start(2) ! charlen and dim1 
      integer :: count(2) ! charlen and dim1
      integer :: varid
      integer :: status

      ! Get netCDF ID of array
      status = nf_inq_varid(ncid,trim(varname),varid)
      call handle_netcdf_error(status)

      ! Now read the array
      start = (/ 1, 1 /)
      count = (/ charlen, dim1 /)
      status = nf_get_vara_text(ncid, varid, start, count, array)
      call handle_netcdf_error(status)

   end subroutine read_netcdf_character_array_1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_dimension
   !
   ! DESCRIPTION:  Retrieves length of a netCDF dimension.
   !
   !---------------------------------------------------------------------------

   function read_netcdf_dimension(ncid,dim_name) result (dim_len)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: dim_name

      ! Return variable
      integer :: dim_len

      ! Local variable
      integer :: dim_id

      ! Get dimension ID
      dim_id = read_netcdf_dim_id(ncid,trim(dim_name))

      ! Get the dimension length
      dim_len = read_netcdf_dim_len(ncid,dim_id)

      return
   end function read_netcdf_dimension

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  write_netcdf_dimension
   !
   ! DESCRIPTION:  Writes a dimension in a netCDF file.
   !
   !---------------------------------------------------------------------------

   function write_netcdf_dimension(nc_id,dim_name,dim) result (dim_id)

      ! Arguments
      integer,intent(in) :: nc_id
      character(len=*),intent(in) :: dim_name
      integer,intent(in) :: dim

      ! Return variable
      integer :: dim_id

      ! Local variables
      integer :: status

      ! Define the dimension in the file
      status = nf_def_dim(nc_id,trim(dim_name),dim,dim_id)
      call handle_netcdf_error(status)

      return
   end function write_netcdf_dimension

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  open_netcdf_newfile
   !
   ! DESCRIPTION:  Opens a new netCDF file.
   !
   !---------------------------------------------------------------------------

   function open_netcdf_newfile(filename) result (nc_id)

      ! Arguments
      character(len=*),intent(in) :: filename

      ! Return variable
      integer :: nc_id

      ! Local variables
      integer :: status

      ! Open new netCDF file
      status = nf_create(trim(filename),NF_CLOBBER,nc_id)
      call handle_netcdf_error(status)

      return
   end function open_netcdf_newfile


   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  check_define_netcdf_real_array_4d
   !
   ! DESCRIPTION:  Checks if a 4D real array with specified dimension has
   ! already been defined in a netCDF file.  If not, defines it.
   ! Values are written separately.
   !
   !---------------------------------------------------------------------------

   function check_define_netcdf_real_array_4d(ncid,variable_name, &
        dim1_name,dim2_name,dim3_name,dim4_name) result (var_id)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: variable_name
      character(len=*),intent(in) :: dim1_name
      character(len=*),intent(in) :: dim2_name
      character(len=*),intent(in) :: dim3_name
      character(len=*),intent(in) :: dim4_name

      ! Return variable
      integer :: var_id

      ! Local variables
      integer :: dim_ids(4), dim_ids2(4)
      integer :: xtype
      integer :: ndims
      integer :: status
      integer :: i,j
      logical :: found

      ! See if variable with same name already exists in netCDF file.
      status = nf_inq_varid(ncid,trim(variable_name),var_id)
      if (status .ne. NF_NOERR .and. status .ne. NF_ENOTVAR) then
         call handle_netcdf_error(status)
      end if

      ! If there is no name conflict, define as a new array and return.
      if (status .eq. NF_ENOTVAR) then
         var_id = define_netcdf_real_array_4d(ncid,variable_name, &
              dim1_name,dim2_name,dim3_name,dim4_name) 
         return
      end if

      ! There is a name conflict.  See if the data types conflict.
      status = nf_inq_vartype(ncid,var_id,xtype)
      call handle_netcdf_error(status)
      if (xtype .ne. NF_FLOAT) then
         print*,'ERROR, variable ',trim(variable_name),' already exists'
         print*,'in netCDF file and has a different data type!'
         print*,'Aborting...'
         stop
      end if

      ! See if the dimensions conflict
      status = nf_inq_varndims(ncid,var_id,ndims)
      call handle_netcdf_error(status)
      if (ndims .ne. 4) then
         print*,'ERROR, variable ',trim(variable_name),' already exists'
         print*,'in netCDF file and has a different rank!'
         print*,'Found ',ndims,' dimensions, expected 4!'
         stop
      end if
      status = nf_inq_vardimid(ncid,var_id,dim_ids2)
      call handle_netcdf_error(status)
      dim_ids(1) = read_netcdf_dim_id(ncid,trim(dim1_name))
      dim_ids(2) = read_netcdf_dim_id(ncid,trim(dim2_name))
      dim_ids(3) = read_netcdf_dim_id(ncid,trim(dim3_name))
      dim_ids(4) = read_netcdf_dim_id(ncid,trim(dim4_name))
      do i = 1,ndims
         found=.false.
         do j = 1,ndims
            if (dim_ids(i) .eq. dim_ids2(j)) then
               found=.true.
               exit
            end if
         end do
         if (.not. found) then
            print*,'ERROR, variable ',trim(variable_name),' already exists'
            print*,'in netCDF file and has different dimensions!'
            stop
         end if
      end do

      ! Although the variable already exists in the netCDF file, it has
      ! the same data type and dimensions.  We'll just return the var_id.
      return
   end function check_define_netcdf_real_array_4d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  check_real_array_4d_netcdf4_settings
   !
   ! DESCRIPTION:  Checks 4D real array in opened netCDF file for netCDF4
   ! specific information (chunking and compression settings).
   !
   !---------------------------------------------------------------------------

   subroutine get_real_array_4d_netcdf4_settings(ncid,varname, &
        is_netcdf4,storage,chunksizes,shuffle,deflate,deflate_level)
      implicit none

      ! Arguments
      integer, intent(in) :: ncid
      character(len=*),intent(in) :: varname
      logical,intent(out) :: is_netcdf4
      integer,intent(out) :: storage
      integer,intent(out) :: chunksizes(4)
      integer,intent(out) :: shuffle
      integer,intent(out) :: deflate
      integer,intent(out) :: deflate_level

      ! Local variables
      integer :: var_id
      integer :: status
      integer :: iformat

      ! Initialize output
      is_netcdf4 = .false.
      storage = 0
      chunksizes = 0
      shuffle = 0
      deflate = 0
      deflate_level = 0

      ! Use default settings if file is not netCDF4
      status = nf_inq_format(ncid,iformat)
      call handle_netcdf_error(status)
      if ( ( iformat .eq. NF_FORMAT_CLASSIC) .or. &
           ( iformat .eq. NF_FORMAT_64BIT  ) ) return

      ! At this point we've established the file is netCDF4.
      is_netcdf4 = .true.

      ! Get netCDF variable ID of requested variable.
      status = nf_inq_varid(ncid,trim(varname),var_id)
      call handle_netcdf_error(status)

      ! Get chunking settings
      status = nf_inq_var_chunking(ncid,var_id,storage,chunksizes)
      call handle_netcdf_error(status)

      ! Get compression settings
      status = nf_inq_var_deflate(ncid,var_id,shuffle,deflate,deflate_level)
      call handle_netcdf_error(status)
      
      ! The end
      return
   end subroutine get_real_array_4d_netcdf4_settings

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  define_real_array_4d_netcdf4_settings
   !
   ! DESCRIPTION:  Defines netCDF4 chunking and compression settings for new
   ! 4D real array.
   !
   !---------------------------------------------------------------------------

   subroutine define_real_array_4d_netcdf4_settings(ncid,var_id, &
        storage,chunksizes,shuffle,deflate,deflate_level)

      implicit none

      ! Arguments
      integer,intent(in) :: ncid
      integer,intent(in) :: var_id
      integer,intent(in) :: storage
      integer,intent(in) :: chunksizes(4)
      integer,intent(in) :: shuffle
      integer,intent(in) :: deflate
      integer,intent(in) :: deflate_level

      ! Local variables
      integer :: status

      status = nf_def_var_chunking(ncid,var_id,storage,chunksizes)
      ! It's possible chunking has already been set and can't be redefined.
      ! Treat this scenario as benign.
      if (status .eq. NF_ELATEDEF) then
!         print*,'Warning, chunking already set for variable in file.'
      else
         call handle_netcdf_error(status)
      end if

      status = nf_def_var_deflate(ncid,var_id,shuffle,deflate,deflate_level)
      ! It's possible deflate has already been set and can't be redefined.
      ! Treat this scenario as benign.
      if (status .eq. NF_ELATEDEF) then
!         print*,'Warning, deflate already set for variable in file.'
      else
         call handle_netcdf_error(status)
      end if

      return
   end subroutine define_real_array_4d_netcdf4_settings

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  read_netcdf_real_variable
   !
   ! DESCRIPTION:  Reads in and returns real variable from netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine read_netcdf_real_variable(ncid,varname,var)

      ! Arguments
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: varname
      real, intent(out) :: var

      ! Local variables
      integer :: varid
      integer :: status

      ! Get netCDF ID of variable
      status = nf_inq_varid(ncid,trim(varname),varid)
      call handle_netcdf_error(status)

      ! Now read the variable
      status = nf_get_var_real(ncid, varid, var)
      call handle_netcdf_error(status)

   end subroutine read_netcdf_real_variable

end module netcdf_util_mod
