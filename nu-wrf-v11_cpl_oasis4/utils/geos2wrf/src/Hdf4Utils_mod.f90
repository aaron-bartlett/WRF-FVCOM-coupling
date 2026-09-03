!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  Hdf4Utils_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! Utility subroutines for calling the HDF4 library.
!
!------------------------------------------------------------------------------

module Hdf4Utils_mod

   ! Change defaults
   implicit none
   private

   ! Include HDF4 constants and routines.  Keep private.
   include 'hdf.f90'
   include 'mffunc.f90'

   ! Public subroutines
   public :: openHdf4Readfile
   public :: closeHdf4file
   public :: readHdf4IntegerArray1d
   public :: readHdf4DoubleArray1d
   public :: readHdf4RealArray3d
   public :: readHdf4RealArray4d

   ! Useful constants (found in HDF4 C header files, but not in Fortran)
   integer, parameter :: MAX_VAR_DIMS = 32
   integer, parameter :: MAX_SDS_NAME_LEN = 64

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  openHdf4Readfile
   !
   ! DESCRIPTION:  Opens a HDF4 file in read mode, returns file ID.
   !
   !---------------------------------------------------------------------------

   function openHdf4Readfile(filename) result (sdID)
 
      ! Arguments
      character(len=*),intent(in) :: filename

      ! Return value
      integer :: sdID

      ! Open the file
      sdID = sfstart(trim(filename),DFACC_READ)
      if (sdID == FAIL) then
         print*,'ERROR opening HDF4 file ',trim(filename)
         stop 1
      end if

      return
   end function openHdf4Readfile

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  closeHdf4File
   !
   ! DESCRIPTION:  Closes a HDF4 file.
   !
   !---------------------------------------------------------------------------

   subroutine closeHdf4File(sdID)
      
      ! Arguments
      integer,intent(in) :: sdID

      ! Local variables
      integer :: status

      ! Close the file
      status = sfend(sdID)
      if (status == FAIL) then
         print*,'ERROR closing HDF4 file!'
         stop 1
      end if
      
   end subroutine closeHdf4File

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  findSdsId
   !
   ! DESCRIPTION:  Finds the ID of a particular SDS array.  Private method.
   !
   !---------------------------------------------------------------------------

   function findSdsId(sdID,variableName,rank,dataType) result (sdsID)

      ! Arguments
      integer,intent(in) :: sdID
      character(len=*),intent(in) :: variableName
      integer,intent(in) :: rank
      integer,intent(in) :: dataType

      ! Return variable
      integer :: sdsID

      ! Local variables
      integer :: tmp_sdsID
      integer :: status
      character(len=MAX_SDS_NAME_LEN) :: sdsName
      integer,dimension(MAX_VAR_DIMS) :: dimensionSizes
      integer :: numberDatasets, numberFileAttributes
      integer :: tmp_rank, tmp_dataType, numberAttributes
      integer :: i

      ! Get contents of file
      status = sffinfo(sdID, numberDatasets, numberFileAttributes)
      if (status == FAIL) then
         print*,'ERROR getting contents from file!'
         stop 1
      end if

      sdsID = -1

      ! Search for variable
      do i = 0, numberDatasets - 1

         ! First, get the ID of the current variable
         tmp_sdsID = sfselect(sdID,i)
         if (tmp_sdsID == FAIL) then
            print*,'ERROR accessing dataset in file!'
            stop 1
         end if

         ! Next, get basic information about the array.
         status = sfginfo(tmp_sdsID, sdsName, tmp_rank, dimensionSizes, &
              tmp_dataType, numberAttributes)
         if (status == FAIL) then
            print*,'ERROR retrieving dataset information from file!'
            stop 1
         end if

         ! Check the name of the SDS array
         if (trim(sdsName) .ne. trim(variableName)) cycle

         ! Check the rank of the SDS array
         if (tmp_rank .ne. rank) cycle

         ! Check the datatype of the SDS array
         if (tmp_dataType .ne. dataType) cycle

         ! If we're here, we found the array we're looking for.
         sdsID = tmp_sdsID
         exit

      end do

      ! Sanity check
      if (sdsID .eq. -1) then
         print*,'ERROR, cannot find array in HDF4 file!'
         print*,'Looked for name ',trim(variableName)
         print*,'With rank ',rank
         print*,'And type ',dataType
         stop 1
      end if

      return
   end function findSdsId

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readHdf4IntegerArray1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D integer array in 
   ! HDF4 file.
   !
   !---------------------------------------------------------------------------

   subroutine readHdf4IntegerArray1d(sdID,variableName,dim1,array, &
        subset,istart1,iend1)

      use CheckSubsetBounds_mod

      ! Arguments
      integer,intent(in) :: sdID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1
      integer, allocatable, intent(out) :: array(:)
      logical,intent(in),optional :: subset
      integer,intent(in),optional :: istart1,iend1

      ! Local variables
      integer :: start(1)
      integer :: stride(1)
      integer :: edges(1)
      integer :: sdsID
      integer :: status
      character(len=MAX_SDS_NAME_LEN) :: sdsName
      integer,dimension(MAX_VAR_DIMS) :: dimensionSizes
      integer :: rank, dataType, numberAttributes

      ! Sanity check optional arguments
      if (present(subset)) then
         if (subset) then
            if (.not. present(istart1) .or. &
                .not. present(iend1)) then
               print*,'ERROR, subset is set to true, but missing index bounds!'
               stop 1
            end if
         end if
      end if

      ! Set basic requirements of SDS array
      rank = 1
      dataType = DFNT_INT32
      
      ! Get ID of variable
      sdsID = findSdsId(sdID,variableName,rank,dataType)

      ! Retrieve the dimension sizes
      status = sfginfo(sdsID, sdsName, rank, dimensionSizes, dataType, &
           numberAttributes)
      if (status == FAIL) then
         print*,'ERROR retrieving dimensions from HDF4 file for variable  ', &
              trim(variableName)
         stop 1
      end if

      ! Subset if requested
      start(1:1) = 0
      if (present(subset)) then
         if (subset) then
            call checkSubsetBounds1d(istart1,iend1)
            dimensionSizes(1) = iend1 - istart1 + 1
            start(1) = istart1-1
         end if
      end if

      ! Allocate memory for the SDS array
      dim1 = dimensionSizes(1)
      allocate(array(dim1))

      ! Read the data
      edges(1:1) = dimensionSizes(1:1)
      stride(1:1) = 1
      status = sfrdata(sdsID, start, stride, edges, array)
      if (status == FAIL) then
         print*,'ERROR reading array ',trim(variableName), ' from HDF4 file!'
         stop 1
      end if
         
      return
   end subroutine readHdf4IntegerArray1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readHdf4DoubleArray1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D double array in 
   ! HDF4 file.
   !
   !---------------------------------------------------------------------------

   subroutine readHdf4DoubleArray1d(sdID,variableName,dim1,array, &
        subset,istart1,iend1)

      use CheckSubsetBounds_mod

      ! Arguments
      integer,intent(in) :: sdID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1
      double precision, allocatable, intent(out) :: array(:)
      logical,intent(in),optional :: subset
      integer,intent(in),optional :: istart1,iend1

      ! Local variables
      integer :: start(1)
      integer :: stride(1)
      integer :: edges(1)
      integer :: sdsID
      integer :: status
      character(len=MAX_SDS_NAME_LEN) :: sdsName
      integer,dimension(MAX_VAR_DIMS) :: dimensionSizes
      integer :: rank, dataType, numberAttributes

      ! Sanity check optional arguments
      if (present(subset)) then
         if (subset) then
            if (.not. present(istart1) .or. &
                .not. present(iend1)) then
               print*,'ERROR, subset is set to true, but missing index bounds!'
               stop 1
            end if
         end if
      end if

      ! Set basic requirements of SDS array
      rank = 1
      dataType = DFNT_FLOAT64
      
      ! Get ID of variable
      sdsID = findSdsId(sdID,variableName,rank,dataType)

      ! Retrieve the dimension sizes
      status = sfginfo(sdsID, sdsName, rank, dimensionSizes, dataType, &
           numberAttributes)
      if (status == FAIL) then
         print*,'ERROR retrieving dimensions from HDF4 file for variable  ', &
              trim(variableName)
         stop 1
      end if

      ! Subset if requested
      start(1:1) = 0
      if (present(subset)) then
         if (subset) then
            call checkSubsetBounds1d(istart1,iend1)
            dimensionSizes(1) = iend1 - istart1 + 1
            start(1) = istart1 - 1
         end if
      end if

      ! Allocate memory for the SDS array
      dim1 = dimensionSizes(1)
      allocate(array(dim1))

      ! Read the data

      edges(1:1) = dimensionSizes(1:1)
      stride(1:1) = 1
      status = sfrdata(sdsID, start, stride, edges, array)
      if (status == FAIL) then
         print*,'ERROR reading array ',trim(variableName), ' from HDF4 file!'
         stop 1
      end if
         
      return
   end subroutine readHdf4DoubleArray1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readHdf4RealArray3d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 3D real array in 
   ! HDF4 file.
   !
   !---------------------------------------------------------------------------

   subroutine readHdf4RealArray3d(sdID,variableName,dim1,dim2,dim3,array, &
        subset,istart1,istart2,istart3,iend1,iend2,iend3)

      use CheckSubsetBounds_mod

      ! Arguments
      integer,intent(in) :: sdID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1,dim2,dim3
      real, allocatable, intent(out) :: array(:,:,:)
      logical,intent(in),optional :: subset
      integer,intent(in),optional :: istart1,istart2,istart3
      integer,intent(in),optional :: iend1,iend2,iend3

      ! Local variables
      integer :: start(3)
      integer :: stride(3)
      integer :: edges(3)
      integer :: sdsID
      integer :: status
      character(len=MAX_SDS_NAME_LEN) :: sdsName
      integer,dimension(MAX_VAR_DIMS) :: dimensionSizes
      integer :: rank, dataType, numberAttributes

      ! Sanity check optional arguments
      if (present(subset)) then
         if (subset) then
            if (.not. present(istart1) .or. &
                .not. present(istart2) .or. &
                .not. present(istart3) .or. &
                .not. present(iend1) .or. &
                .not. present(iend2) .or. &
                .not. present(iend3)) then
               print*,'ERROR, subset is set to true, but missing index bounds!'
               stop 1
            end if
         end if
      end if

      ! Set basic requirements of SDS array
      rank = 3
      dataType = DFNT_FLOAT32
      
      ! Get ID of variable
      sdsID = findSdsId(sdID,variableName,rank,dataType)

      ! Retrieve the dimension sizes
      status = sfginfo(sdsID, sdsName, rank, dimensionSizes, dataType, &
           numberAttributes)
      if (status == FAIL) then
         print*,'ERROR retrieving dimensions from HDF4 file for variable  ', &
              trim(variableName)
         stop 1
      end if

      ! Subset if requested
      start(1:3) = 0
      if (present(subset)) then
         if (subset) then
            call checkSubsetBounds3d(istart1,iend1,istart2,iend2,istart3,iend3)
            dimensionSizes(1) = iend1 - istart1 + 1
            dimensionSizes(2) = iend2 - istart2 + 1
            dimensionSizes(3) = iend3 - istart3 + 1
            start(1) = istart1 - 1
            start(2) = istart2 - 1
            start(3) = istart3 - 1
         end if
      end if

      ! Allocate memory for the SDS array
      dim1 = dimensionSizes(1)
      dim2 = dimensionSizes(2)
      dim3 = dimensionSizes(3)
      allocate(array(dim1,dim2,dim3))

      ! Read the data
      edges(1:3) = dimensionSizes(1:3)
      stride(1:3) = 1
      status = sfrdata(sdsID, start, stride, edges, array)
      if (status == FAIL) then
         print*,'ERROR reading array ',trim(variableName), ' from HDF4 file!'
         stop 1
      end if
         
      return
   end subroutine readHdf4RealArray3d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readHdf4RealArray4d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 4D real array in 
   ! HDF4 file.
   !
   !---------------------------------------------------------------------------

   subroutine readHdf4RealArray4d(sdID,variableName,dim1,dim2,dim3,dim4, &
        array, &
        subset,istart1,istart2,istart3,istart4,iend1,iend2,iend3,iend4)

      use CheckSubsetBounds_mod

      ! Arguments
      integer,intent(in) :: sdID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1,dim2,dim3,dim4
      real, allocatable, intent(out) :: array(:,:,:,:)
      logical,intent(in),optional :: subset
      integer,intent(in),optional :: istart1,istart2,istart3,istart4
      integer,intent(in),optional :: iend1,iend2,iend3,iend4

      ! Local variables
      integer :: start(4)
      integer :: stride(4)
      integer :: edges(4)
      integer :: sdsID
      integer :: status
      character(len=MAX_SDS_NAME_LEN) :: sdsName
      integer,dimension(MAX_VAR_DIMS) :: dimensionSizes
      integer :: rank, dataType, numberAttributes

      ! Sanity check optional arguments
      if (present(subset)) then
         if (subset) then
            if (.not. present(istart1) .or. &
                .not. present(istart2) .or. &
                .not. present(istart3) .or. &
                .not. present(istart4) .or. &
                .not. present(iend1) .or. &
                .not. present(iend2) .or. &
                .not. present(iend3) .or. &
                .not. present(iend4)) then
               print*,'ERROR, subset is set to true, but missing index bounds!'
               stop 1
            end if
         end if
      end if

      ! Set basic requirements of SDS array
      rank = 4
      dataType = DFNT_FLOAT32
      
      ! Get ID of variable
      sdsID = findSdsId(sdID,variableName,rank,dataType)

      ! Retrieve the dimension sizes
      status = sfginfo(sdsID, sdsName, rank, dimensionSizes, dataType, &
           numberAttributes)
      if (status == FAIL) then
         print*,'ERROR retrieving dimensions from HDF4 file for variable  ', &
              trim(variableName)
         stop 1
      end if

      ! Subset if requested
      start(1:4) = 0
      if ( present(subset)) then
         if (subset) then
            call checkSubsetBounds4d(istart1,iend1,istart2,iend2, &
                 istart3,iend3, &
                 istart4,iend4)
            dimensionSizes(1) = iend1 - istart1 + 1
            dimensionSizes(2) = iend2 - istart2 + 1
            dimensionSizes(3) = iend3 - istart3 + 1
            dimensionSizes(4) = iend4 - istart4 + 1
            start(1) = istart1
            start(2) = istart2
            start(3) = istart3
            start(4) = istart4
         end if
      end if

      ! Allocate memory for the SDS array
      dim1 = dimensionSizes(1)
      dim2 = dimensionSizes(2)
      dim3 = dimensionSizes(3)
      dim4 = dimensionSizes(4)
      allocate(array(dim1,dim2,dim3,dim4))

      ! Read the data
      edges(1:4) = dimensionSizes(1:4)
      stride(1:4) = 1
      status = sfrdata(sdsID, start, stride, edges, array)
      if (status == FAIL) then
         print*,'ERROR reading array ',trim(variableName), ' from HDF4 file!'
         stop 1
      end if
         
      return
   end subroutine readHdf4RealArray4d

end module Hdf4Utils_mod
