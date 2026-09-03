#ifdef WRF_LIS
!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  NetcdfUtils_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Utility subroutines for calling the netCDF library.
!
!------------------------------------------------------------------------------

module NetcdfUtils_mod

   ! Import modules
   use module_wrf_error

   ! Change defaults
   implicit none
   private

   ! Include netCDF constants and routines.  Keep private.
   include 'netcdf.inc'

   ! Public subroutines
   public :: closeNetcdfFile
   public :: openNetcdfReadfile
   public :: openNetcdfWritefile
   public :: readNetcdfDoubleArray1d
   public :: readNetcdfGlobalIntegerAttribute
   public :: readNetcdfGlobalTextAttribute
   public :: readNetcdfIntegerArray1d
   public :: readNetcdfRealArray1d
   public :: readNetcdfRealArray2d
   public :: checkNetcdfRealArray2d
   public :: readNetcdfRealArray3d
   public :: readNetcdfRealArray4d
   public :: writeNetcdfRealArray3dUpdate
   public :: writeNetcdfRealArray4dUpdate
   public :: readNetcdfRealGlobalAttribute
   public :: readNetcdfTextGlobalAttribute
   public :: readNetcdfIntegerGlobalAttribute
   public :: findVarID

   ! Local constants
   integer,parameter :: MAX_MESSAGE_LEN = 256

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  openNetcdfReadfile
   !
   ! DESCRIPTION:  Opens a netCDF file in read mode, returns file ID.
   !
   !---------------------------------------------------------------------------

   function openNetcdfReadfile(filename) result (ncID)

      ! Arguments
      character(len=*),intent(in) :: filename

      ! Return value
      integer :: ncID

      ! Local variables
      integer :: status
      character(len=MAX_MESSAGE_LEN) :: message

      ! Open the file
      status = nf_open(trim(filename),NF_NOWRITE,ncID)
      if (status .ne. NF_NOERR) then
         write(message,*)'Problem opening file ',trim(filename), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if 

      return
   end function openNetcdfReadfile

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  openNetcdfWritefile
   !
   ! DESCRIPTION:  Opens an existing netCDF file in write mode, returns file 
   ! ID.
   !
   !---------------------------------------------------------------------------

   function openNetcdfWritefile(filename) result (ncID)

      ! Arguments
      character(len=*),intent(in) :: filename

      ! Return value
      integer :: ncID

      ! Local variables
      integer :: status
      character(len=MAX_MESSAGE_LEN) :: message

      ! Open the file
      status = nf_open(trim(filename),NF_WRITE,ncID)
      if (status .ne. NF_NOERR) then
         write(message,*)'Problem opening file ',trim(filename), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if 

      return
   end function openNetcdfWritefile

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  closeNetcdfFile
   !
   ! DESCRIPTION:  Closes a netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine closeNetcdfFile(ncID)

      ! Arguments
      integer,intent(in) :: ncID

      ! Local variables
      integer :: status
      character(len=MAX_MESSAGE_LEN) :: message

      ! Close the file
      status = nf_close(ncID)
      if (status .ne. NF_NOERR) then
         write(message,*)'Problem closing netCDF file, ', &
              trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      return
   end subroutine closeNetcdfFile
 
!
!---------------------------------------------------------------------------
   !
   ! ROUTINE:  checkVarId
   !
   ! DESCRIPTION:  Check the ID of a particular netCDF variable.  Private 
   ! method.
   !
   !---------------------------------------------------------------------------

   function checkVarId(ncID,variableName) result (CheckStatus)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer :: varID
      ! Return variable
      integer :: CheckStatus
      ! Local variables
      integer :: status

      ! Get ID of netCDF variable
      status = nf_inq_varid(ncID,trim(variableName),varID)
      if (status .ne. NF_NOERR) then
         !write(message,*)'Problem finding variable ',trim(variableName), &
         !     ', ',trim(nf_strerror(status))
              CheckStatus = -1
         else
              CheckStatus = 1
      end if         
      return
    end function checkVarId
!
     
   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  findVarId
   !
   ! DESCRIPTION:  Finds the ID of a particular netCDF variable.  Private 
   ! method.
   !
   !---------------------------------------------------------------------------

   function findVarId(ncID,variableName,rank,dataType) result (varID)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(in) :: rank
      integer,intent(in) :: dataType

      ! Return variable
      integer :: varID

      ! Local variables
      integer :: status
      integer :: tmp_rank, tmp_dataType
      character(len=MAX_MESSAGE_LEN) :: message

      ! Get ID of netCDF variable
      status = nf_inq_varid(ncID,trim(variableName),varID)
      if (status .ne. NF_NOERR) then
         write(message,*)'Problem finding variable ',trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Sanity check the rank
      status = nf_inq_varndims(ncID, varID, tmp_rank)
      if (status .ne. NF_NOERR) then
         write(message,*)'Problem finding rank of variable ', &
              trim(variableName),', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
      if (tmp_rank .ne. rank) then
         write(message,*)'Rank mismatch for variable ',trim(variableName), &
              ', expected ',rank,', found ',tmp_rank
         call wrf_error_fatal(trim(message))
      end if

      ! Sanity check the data type
      status = nf_inq_vartype(ncID, varID, tmp_dataType)
      if (status .ne. NF_NOERR) then
         write(message,*)'Problem finding type of variable ', &
              trim(variableName),', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
      if (tmp_dataType .ne. dataType) then
         write(message,*)'Type mismatch for variable ', &
              trim(variableName),', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      return
   end function findVarId

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNetcdfIntegerArray1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D integer array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readNetcdfIntegerArray1d(ncID,variableName,dim1,array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1
      integer, allocatable, intent(out) :: array(:)

      ! Local variables
      integer :: start(1)
      integer :: count(1)
      integer :: varID
      integer :: dimIDs(1)
      integer :: dimensionSizes(1)
      integer :: rank, dataType
      integer :: status
      integer :: i
      character(len=MAX_MESSAGE_LEN) :: message

      ! Set basic requirements for variable
      rank = 1
      dataType = NF_INT

      ! Get ID of variable
      varID = findVarID(ncID,variableName,rank,dataType)
      ! Retrieve the dimension IDs
      status = nf_inq_vardimid(ncID,varID,dimIDs)
      if (status .ne. NF_NOERR) then
         write(message,*)'Problem retrieving dimension IDs for variable ', &
              trim(variableName),', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Get dimensions
      do i = 1, rank
         status = nf_inq_dimlen(ncID,dimIDs(i),dimensionSizes(i))
         if (status .ne. NF_NOERR) then
            write(message,*)'Problem retrieving dimension length from netCDF',&
                 ', ',trim(nf_strerror(status))
            call wrf_error_fatal(trim(message))
         end if
      end do

      ! Allocate the memory
      dim1 = dimensionSizes(1)
      allocate(array(dim1),stat=status)
      if (status .ne. 0) then
         write(message,*)'Allocation failure!'
         call wrf_error_fatal(trim(message))
      end if

      ! Read the data
      start(1:1) = 1
      count(1:1) = dimensionSizes(1:1)
      status = nf_get_vara_int(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*)'Problem reading netCDF variable ', &
              trim(variableName),', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
         
      return
   end subroutine readNetcdfIntegerArray1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNetcdfDoubleArray1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D double array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readNetcdfDoubleArray1d(ncID,variableName,dim1,array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1
      double precision, allocatable, intent(out) :: array(:)

      ! Local variables
      integer :: start(1)
      integer :: count(1)
      integer :: varID
      integer :: dimIDs(1)
      integer :: dimensionSizes(1)
      integer :: rank, dataType
      integer :: status
      integer :: i
      character(len=MAX_MESSAGE_LEN) :: message

      ! Set basic requirements for variable
      rank = 1
      dataType = NF_DOUBLE

      ! Get ID of variable
      varID = findVarID(ncID,variableName,rank,dataType)

      ! Retrieve the dimension IDs
      status = nf_inq_vardimid(ncID,varID,dimIDs)
      if (status .ne. NF_NOERR) then
         write(message,*)'Problem retrieving dimension IDs for variable ', &
              trim(variableName),', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Get dimensions
      do i = 1, rank
         status = nf_inq_dimlen(ncID,dimIDs(i),dimensionSizes(i))
         if (status .ne. NF_NOERR) then
            write(message,*) &
                 'Problem retrieving dimension length from netCDF ', &
                 ', ',trim(nf_strerror(status))
            call wrf_error_fatal(trim(message))
         end if
      end do

      ! Allocate the memory
      dim1 = dimensionSizes(1)
      allocate(array(dim1),stat=status)
      if (status .ne. 0) then
         write(message,*)'Allocation failure!'
         call wrf_error_fatal(trim(message))
      end if

      ! Read the data
      start(1:1) = 1
      count(1:1) = dimensionSizes(1:1)
      status = nf_get_vara_double(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving variable ', trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
         
      return
   end subroutine readNetcdfDoubleArray1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  checkNetcdfRealArray2d
   !
   ! DESCRIPTION:  Reads in and checks if entire requested 2D real array 
   ! exists in the netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine checkNetcdfRealArray2d(ncID,variableName,CheckStatus)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: CheckStatus

      ! Get ID of variable
      CheckStatus = checkVarID(ncID,variableName)
      return
   end subroutine checkNetcdfRealArray2d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNetcdfRealArray1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D real array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readNetcdfRealArray1d(ncID,variableName,dim1,array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1
      real, allocatable, intent(out) :: array(:)

      ! Local variables
      integer :: start(1)
      integer :: count(1)
      integer :: varID
      integer :: dimIDs(1)
      integer :: dimensionSizes(1)
      integer :: rank, dataType
      integer :: status
      integer :: i
      character(len=MAX_MESSAGE_LEN) :: message

      ! Set basic requirements for variable
      rank = 1
      dataType = NF_FLOAT

      ! Get ID of variable
      varID = findVarID(ncID,variableName,rank,dataType)

      ! Retrieve the dimension IDs
      status = nf_inq_vardimid(ncID,varID,dimIDs)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving dimension IDs for variable ', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Get dimensions
      do i = 1, rank
         status = nf_inq_dimlen(ncID,dimIDs(i),dimensionSizes(i))
         if (status .ne. NF_NOERR) then
            write(message,*) &
                 'Problem retrieving dimension length from netCDF', &
                 ', ',trim(nf_strerror(status))
            call wrf_error_fatal(trim(message))
         end if
      end do

      ! Allocate the memory
      dim1 = dimensionSizes(1)
      allocate(array(dim1),stat=status)
      if (status .ne. 0) then
         write(message,*)'Allocation failure!'
         call wrf_error_fatal(trim(message))
      end if

      ! Read the data
      start(1) = 1
      count(1) = dimensionSizes(1)
      status = nf_get_vara_real(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem reading variable', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
         
      return
   end subroutine readNetcdfRealArray1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNetcdfRealArray2d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 2D real array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readNetcdfRealArray2d(ncID,variableName,dim1,dim2,array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1,dim2
      real, allocatable, intent(out) :: array(:,:)

      ! Local variables
      integer :: start(2)
      integer :: count(2)
      integer :: varID
      integer :: dimIDs(2)
      integer :: dimensionSizes(2)
      integer :: rank, dataType
      integer :: status
      integer :: i
      character(len=MAX_MESSAGE_LEN) :: message

      ! Set basic requirements for variable
      rank = 2
      dataType = NF_FLOAT

      ! Get ID of variable
      varID = findVarID(ncID,variableName,rank,dataType)

      ! Retrieve the dimension IDs
      status = nf_inq_vardimid(ncID,varID,dimIDs)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving dimension IDs for variable ', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Get dimensions
      do i = 1, rank
         status = nf_inq_dimlen(ncID,dimIDs(i),dimensionSizes(i))
         if (status .ne. NF_NOERR) then
            write(message,*) &
                 'Problem retrieving dimension length from netCDF', &
                 ', ',trim(nf_strerror(status))
            call wrf_error_fatal(trim(message))
         end if
      end do

      ! Allocate the memory
      dim1 = dimensionSizes(1)
      dim2 = dimensionSizes(2)
      allocate(array(dim1,dim2),stat=status)
      if (status .ne. 0) then
         write(message,*)'Allocation failure!'
         call wrf_error_fatal(trim(message))
      end if

      ! Read the data
      start(1:2) = 1
      count(1:2) = dimensionSizes(1:2)
      status = nf_get_vara_real(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem reading variable', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
         
      return
   end subroutine readNetcdfRealArray2d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNetcdfRealArray3d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 3D real array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readNetcdfRealArray3d(ncID,variableName,dim1,dim2,dim3,array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1,dim2,dim3
      real, allocatable, intent(out) :: array(:,:,:)

      ! Local variables
      integer :: start(3)
      integer :: count(3)
      integer :: varID
      integer :: dimIDs(3)
      integer :: dimensionSizes(3)
      integer :: rank, dataType
      integer :: status
      integer :: i
      character(len=MAX_MESSAGE_LEN) :: message

      ! Set basic requirements for variable
      rank = 3
      dataType = NF_FLOAT

      ! Get ID of variable
      varID = findVarID(ncID,variableName,rank,dataType)

      ! Retrieve the dimension IDs
      status = nf_inq_vardimid(ncID,varID,dimIDs)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving dimension IDs for variable', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Get dimensions
      do i = 1, rank
         status = nf_inq_dimlen(ncID,dimIDs(i),dimensionSizes(i))
         if (status .ne. NF_NOERR) then
            write(message,*) &
                 'Problem retrieving dimension length from netCDF', &
                 ', ',trim(nf_strerror(status))
            call wrf_error_fatal(trim(message))
         end if
      end do

      ! Allocate the memory
      dim1 = dimensionSizes(1)
      dim2 = dimensionSizes(2)
      dim3 = dimensionSizes(3)
      allocate(array(dim1,dim2,dim3),stat=status)
      if (status .ne. 0) then
         write(message,*)'Allocation failure!'
         call wrf_error_fatal(trim(message))
      end if

      ! Read the data
      start(1:3) = 1
      count(1:3) = dimensionSizes(1:3)
      status = nf_get_vara_real(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem reading variable', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
         
      return
   end subroutine readNetcdfRealArray3d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNetcdfRealArray4d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 4D real array in 
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readNetcdfRealArray4d(ncID,variableName,dim1,dim2,dim3,dim4, &
        array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1,dim2,dim3,dim4
      real, allocatable, intent(out) :: array(:,:,:,:)

      ! Local variables
      integer :: start(4)
      integer :: count(4)
      integer :: varID
      integer :: dimIDs(4)
      integer :: dimensionSizes(4)
      integer :: rank, dataType
      integer :: status
      integer :: i
      character(len=MAX_MESSAGE_LEN) :: message

      ! Set basic requirements for variable
      rank = 4
      dataType = NF_FLOAT

      ! Get ID of variable
      varID = findVarID(ncID,variableName,rank,dataType)

      ! Retrieve the dimension IDs
      status = nf_inq_vardimid(ncID,varID,dimIDs)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving dimension IDs for variable', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Get dimensions
      do i = 1, rank
         status = nf_inq_dimlen(ncID,dimIDs(i),dimensionSizes(i))
         if (status .ne. NF_NOERR) then
            write(message,*) &
                 'Problem retrieving dimension length from netCDF', &
                 ', ',trim(nf_strerror(status))
            call wrf_error_fatal(trim(message))
         end if
      end do

      ! Allocate the memory
      dim1 = dimensionSizes(1)
      dim2 = dimensionSizes(2)
      dim3 = dimensionSizes(3)
      dim4 = dimensionSizes(4)
      allocate(array(dim1,dim2,dim3,dim4),stat=status)
      if (status .ne. 0) then
         write(message,*)'Allocation failure!'
         call wrf_error_fatal(trim(message))
      end if

      ! Read the data
      start(1:4) = 1
      count(1:4) = dimensionSizes(1:4)
      status = nf_get_vara_real(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem reading variable ', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
         
      return
   end subroutine readNetcdfRealArray4d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNetcdfGlobalTextAttribute
   !
   ! DESCRIPTION:  Reads in and returns global text attribute from
   ! netCDF file.
   !
   !---------------------------------------------------------------------------
      
   subroutine readNetcdfGlobalTextAttribute(ncID,name,value)
      
      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: name
      character(len=*),intent(out) :: value

      ! Local variables
      integer :: status
      character(len=MAX_MESSAGE_LEN) :: message

      status = nf_get_att_text(ncID,NF_GLOBAL,trim(name),value)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem fetching global attribute ', &
              trim(name), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if 

      return      
   end subroutine readNetcdfGlobalTextAttribute

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readNetcdfGlobalIntegerAttribute
   !
   ! DESCRIPTION:  Reads in and returns global integer attribute from
   ! netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readNetcdfGlobalIntegerAttribute(ncID,name,value)
      
      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: name
      integer,intent(out) :: value

      ! Local variables
      integer :: status
      character(len=MAX_MESSAGE_LEN) :: message

      status = nf_get_att_int1(ncID,NF_GLOBAL,trim(name),value)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem fetching global attribute ', &
              trim(name), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if 

      return      
   end subroutine readNetcdfGlobalIntegerAttribute

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  writeNetcdfRealArray3dUpdate
   !
   ! DESCRIPTION:  Writes updated values of a 3D real array already written
   ! to a netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine writeNetcdfRealArray3dUpdate(ncID,variableName,dim1,dim2,dim3, &
        array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(in) :: dim1,dim2,dim3
      real, intent(in) :: array(dim1,dim2,dim3)

      ! Local variables
      integer :: start(3)
      integer :: count(3)
      integer :: varID
      integer :: dimIDs(3)
      integer :: dimensionSizes(3)
      integer :: rank, dataType
      integer :: status
      integer :: i
      character(len=MAX_MESSAGE_LEN) :: message

      ! Set basic requirements for variable
      rank = 3
      dataType = NF_FLOAT

      ! Get ID of variable
      varID = findVarID(ncID,variableName,rank,dataType)

      ! Retrieve the dimension IDs
      status = nf_inq_vardimid(ncID,varID,dimIDs)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving dimension IDs for variable ', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Get dimensions
      do i = 1, rank
         status = nf_inq_dimlen(ncID,dimIDs(i),dimensionSizes(i))
         if (status .ne. NF_NOERR) then
            write(message,*) &
                 'Problem retrieving dimension length from netCDF', &
                 ', ',trim(nf_strerror(status))
            call wrf_error_fatal(trim(message))
         end if
      end do

      ! Make sure dimensions match
      if (dim1 .ne. dimensionSizes(1) .or. &
          dim2 .ne. dimensionSizes(2) .or. &
          dim3 .ne. dimensionSizes(3)) then
         write(message,*) &
              'Dimension mismatch for variable ', &
              trim(variableName),', expected ',dim1,dim2,dim3, &
              ', found ',dimensionSizes(1:3), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Read the data
      start(1:3) = 1
      count(1:3) = dimensionSizes(1:3)
      status = nf_put_vara_real(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem updating variable ', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
         
      return
   end subroutine writeNetcdfRealArray3dUpdate

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  writeNetcdfRealArray4dUpdate
   !
   ! DESCRIPTION:  Writes updated values of a 4D real array already written
   ! to a netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine writeNetcdfRealArray4dUpdate(ncID,variableName,dim1,dim2,dim3, &
        dim4,array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: variableName
      integer,intent(in) :: dim1,dim2,dim3,dim4
      real, intent(in) :: array(dim1,dim2,dim3,dim4)

      ! Local variables
      integer :: start(4)
      integer :: count(4)
      integer :: varID
      integer :: dimIDs(4)
      integer :: dimensionSizes(4)
      integer :: rank, dataType
      integer :: status
      integer :: i
      character(len=MAX_MESSAGE_LEN) :: message

      ! Set basic requirements for variable
      rank = 4
      dataType = NF_FLOAT

      ! Get ID of variable
      varID = findVarID(ncID,variableName,rank,dataType)

      ! Retrieve the dimension IDs
      status = nf_inq_vardimid(ncID,varID,dimIDs)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving dimension IDs for variable ', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Get dimensions
      do i = 1, rank
         status = nf_inq_dimlen(ncID,dimIDs(i),dimensionSizes(i))
         if (status .ne. NF_NOERR) then
            write(message,*) &
                 'Problem retrieving dimension length from netCDF ', &
            ', ',trim(nf_strerror(status))
            call wrf_error_fatal(trim(message))
         end if
      end do

      ! Make sure dimensions match
      if (dim1 .ne. dimensionSizes(1) .or. &
          dim2 .ne. dimensionSizes(2) .or. &
          dim3 .ne. dimensionSizes(3) .or. &
          dim4 .ne. dimensionSizes(4)) then
         write(message,*) &
              'Dimension mismatch for variable ', &
              trim(variableName),', expected ',dim1,dim2,dim3,dim4, &
              ', found ',dimensionSizes(1:4), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      ! Read the data
      start(1:4) = 1
      count(1:4) = dimensionSizes(1:4)
      status = nf_put_vara_real(ncID,varID,start,count,array)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem updating variable ', &
              trim(variableName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
         
      return
   end subroutine writeNetcdfRealArray4dUpdate

   ! Public
   subroutine readNetcdfRealGlobalAttribute(ncID,attributeName,dim1,array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: attributeName
      integer,intent(out) :: dim1
      real, allocatable, intent(out) :: array(:)

      ! Local variables
      integer :: status
      integer :: attributeType
      integer :: attributeLen
      character(len=MAX_MESSAGE_LEN) :: message

      ! Confirm attribute is present and is correct type.
      status = nf_inq_atttype(ncID,NF_GLOBAL,trim(attributeName), &
           attributeType)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving type for attribute ', &
              trim(attributeName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
      if (attributeType .ne. NF_FLOAT) then
         write(message,*) &
              'Invalid data type for attribute ', &
              trim(attributeName), &
              ', expected reals'
         call wrf_error_fatal(trim(message))
      end if

      ! Get length of attribute
      status = nf_inq_attlen(ncID, NF_GLOBAL, trim(attributeName), &
           attributeLen)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving attribute length for ', &
              trim(attributeName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      dim1 = attributeLen
      allocate(array(attributeLen),stat=status)
      if (status .ne. 0) then
         write(message,*)'Allocation failure!'
         call wrf_error_fatal(trim(message))
      end if

      ! Get attribute
      status = nf_get_att_real(ncID, NF_GLOBAL, trim(attributeName), array)
      if (status .ne. 0) then
         write(message,*)'Error reading global attribute ', &
              trim(attributeName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
      
   end subroutine readNetcdfRealGlobalAttribute

   ! Public
   subroutine readNetcdfTextGlobalAttribute(ncID,attributeName,stringMaxLen,&
        string)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: attributeName
      integer,intent(in) :: stringMaxLen
      character(len=stringMaxLen),intent(out) :: string

      ! Local variables
      integer :: status
      integer :: attributeType
      integer :: attributeLen
      character(len=MAX_MESSAGE_LEN) :: message

      ! Confirm attribute is present and is correct type.
      status = nf_inq_atttype(ncID,NF_GLOBAL,trim(attributeName), &
           attributeType)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving type for attribute ', &
              trim(attributeName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
      if (attributeType .ne. NF_CHAR) then
         write(message,*) &
              'Invalid data type for attribute ', &
              trim(attributeName), &
              ', expected characters'
         call wrf_error_fatal(trim(message))
      end if

      ! Get length of attribute
      status = nf_inq_attlen(ncID, NF_GLOBAL, trim(attributeName), &
           attributeLen)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving attribute length for ', &
              trim(attributeName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
      if (attributeLen .gt. stringMaxLen) then
         write(message,*) &
              'Attribute length for ', &
              trim(attributeName), &
              ', is too long, exceeds ',stringMaxLen
         call wrf_error_fatal(trim(message))
      end if

      ! Get attribute
      string = ''
      status = nf_get_att_text(ncID, NF_GLOBAL, trim(attributeName), string)
      if (status .ne. 0) then
         write(message,*)'Error reading global attribute ', &
              trim(attributeName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

   end subroutine readNetcdfTextGlobalAttribute

   ! Public
   subroutine readNetcdfIntegerGlobalAttribute(ncID,attributeName,dim1,array)

      ! Arguments
      integer,intent(in) :: ncID
      character(len=*),intent(in) :: attributeName
      integer,intent(out) :: dim1
      integer, allocatable, intent(out) :: array(:)

      ! Local variables
      integer :: status
      integer :: attributeType
      integer :: attributeLen
      character(len=MAX_MESSAGE_LEN) :: message

      ! Confirm attribute is present and is correct type.
      status = nf_inq_atttype(ncID,NF_GLOBAL,trim(attributeName), &
           attributeType)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving type for attribute ', &
              trim(attributeName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
      if (attributeType .ne. NF_INT) then
         write(message,*) &
              'Invalid data type for attribute ', &
              trim(attributeName), &
              ', expected reals'
         call wrf_error_fatal(trim(message))
      end if

      ! Get length of attribute
      status = nf_inq_attlen(ncID, NF_GLOBAL, trim(attributeName), &
           attributeLen)
      if (status .ne. NF_NOERR) then
         write(message,*) &
              'Problem retrieving attribute length for ', &
              trim(attributeName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if

      dim1 = attributeLen
      allocate(array(attributeLen),stat=status)
      if (status .ne. 0) then
         write(message,*)'Allocation failure!'
         call wrf_error_fatal(trim(message))
      end if

      ! Get attribute
      status = nf_get_att_int(ncID, NF_GLOBAL, trim(attributeName), array)
      if (status .ne. 0) then
         write(message,*)'Error reading global attribute ', &
              trim(attributeName), &
              ', ',trim(nf_strerror(status))
         call wrf_error_fatal(trim(message))
      end if
      
   end subroutine readNetcdfIntegerGlobalAttribute

end module NetcdfUtils_mod
#else
module NetcdfUtils_mod
end module NetcdfUtils_mod
#endif

