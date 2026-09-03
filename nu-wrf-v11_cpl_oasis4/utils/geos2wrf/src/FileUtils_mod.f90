!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  File4Utils_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! Utility subroutines for calling the HDF4 or netCDF libraries.
!
!------------------------------------------------------------------------------

module FileUtils_mod

   ! Import lower-level wrapper libraries
   use HDF4Utils_mod
   use NetcdfUtils_mod

   ! Change defaults
   implicit none
   private

   ! Public subroutines
   public :: openReadfile
   public :: closeFile
   public :: readIntegerArray1d
   public :: readDoubleArray1d
   public :: readRealArray3d
   public :: readRealArray4d

   ! Public constants
   integer,parameter,public :: HDF4_FORMAT    = 1
   integer,parameter,public :: NETCDF_FORMAT  = 2
   integer,parameter,public :: HDFEOS2_FORMAT = 4

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  printError
   !
   ! DESCRIPTION:  Prints error message about file format selection.
   !
   !---------------------------------------------------------------------------

   subroutine printError(fileFormat)
      
      ! Arguments
      integer,intent(in) :: fileFormat

      print*,'ERROR, invalid file format selected!'
      print*,'For HDF4    : ',HDF4_FORMAT
      print*,'For NETCDF  : ',NETCDF_FORMAT
      print*,'For HDFEOS2 : ',HDFEOS2_FORMAT
      print*,'Selected    : ',fileFormat
      
      return
   end subroutine printError

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  openReadFile
   !
   ! DESCRIPTION:  Opens a HDF4 or netCDF file in read mode, returns file ID.
   !
   !---------------------------------------------------------------------------

   function openReadFile(fileFormat,filename) result (fileID)

      ! Arguments
      integer, intent(in) :: fileFormat
      character(len=*),intent(in) :: filename

      ! Return value
      integer :: fileID

      ! Open the file
      select case (fileFormat)
      case (HDF4_FORMAT, HDFEOS2_FORMAT)
         fileID = openHdf4ReadFile(filename)
      case (NETCDF_FORMAT)
         fileID = openNetcdfReadFile(filename)
      case default
         call printError(fileFormat)
         stop 1
      end select

      return
   end function openReadFile

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  closeFile
   !
   ! DESCRIPTION:  Closes a HDF4 or netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine closeFile(fileFormat,fileID)

      ! Arguments
      integer,intent(in) :: fileFormat
      integer,intent(in) :: fileID

      ! Close the file
      select case (fileFormat)
      case (HDF4_FORMAT, HDFEOS2_FORMAT)
         call closeHdf4File(fileID)
      case (NETCDF_FORMAT)
         call closeNetcdfFile(fileID)
      case default
         call printError(fileFormat)
         stop 1
      end select

      return
   end subroutine closeFile

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readIntegerArray1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D integer array in 
   ! HDF4 or netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readIntegerArray1d(fileFormat,fileID,variableName,dim1,array, &
        subset,istart1,iend1)

      ! Arguments
      integer,intent(in) :: fileFormat
      integer,intent(in) :: fileID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1
      integer, allocatable, intent(out) :: array(:)
      logical,intent(in),optional :: subset
      integer,intent(in),optional :: istart1,iend1

      ! Read array
      select case (fileFormat)
      case (HDF4_FORMAT, HDFEOS2_FORMAT)
         call readHdf4IntegerArray1d(fileID,variableName,dim1,array, &
              subset,istart1,iend1)
      case (NETCDF_FORMAT)
         call readNetcdfIntegerArray1d(fileID,variableName,dim1,array, &
              subset,istart1,iend1)
      case default
         call printError(fileFormat)
         stop 1
      end select

      return
   end subroutine readIntegerArray1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readDoubleArray1d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 1D double array in 
   ! HDF4 or netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readDoubleArray1d(fileFormat,fileID,variableName,dim1,array, &
        subset,istart1,iend1)

      ! Arguments
      integer,intent(in) :: fileFormat
      integer,intent(in) :: fileID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1
      double precision, allocatable, intent(out) :: array(:)
      logical,intent(in),optional :: subset
      integer,intent(in),optional :: istart1,iend1

      ! Read array
      select case (fileFormat)
      case (HDF4_FORMAT, HDFEOS2_FORMAT)
         call readHdf4DoubleArray1d(fileID,variableName,dim1,array, &
              subset,istart1,iend1)
      case (NETCDF_FORMAT)
         call readNetcdfDoubleArray1d(fileID,variableName,dim1,array, &
              subset,istart1,iend1)
      case default
         call printError(fileFormat)
         stop 1
      end select

      return
   end subroutine readDoubleArray1d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readRealArray3d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 3D array in 
   ! HDF4 or netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readRealArray3d(fileFormat,fileID,variableName,dim1,dim2,dim3, &
        array, &
        subset,istart1,istart2,istart3,iend1,iend2,iend3)

      ! Arguments
      integer,intent(in) :: fileFormat
      integer,intent(in) :: fileID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1,dim2,dim3
      real, allocatable, intent(out) :: array(:,:,:)
      logical,intent(in),optional :: subset
      integer,intent(in),optional :: istart1,istart2,istart3
      integer,intent(in),optional :: iend1,iend2,iend3

      ! Read array
      select case (fileFormat)
      case (HDF4_FORMAT, HDFEOS2_FORMAT)
         call readHdf4RealArray3d(fileID,variableName,dim1,dim2,dim3,array,&
              subset,istart1,istart2,istart3,iend1,iend2,iend3)
      case (NETCDF_FORMAT)
         call readNetcdfRealArray3d(fileID,variableName,dim1,dim2,dim3, &
              array, &
              subset,istart1,istart2,istart3,iend1,iend2,iend3)
      case default
         call printError(fileFormat)
         stop 1
      end select

      return
   end subroutine readRealArray3d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  readRealArray4d
   !
   ! DESCRIPTION:  Reads in and returns entire requested 4D array in 
   ! HDF4 or netCDF file.
   !
   !---------------------------------------------------------------------------

   subroutine readRealArray4d(fileFormat,fileID,variableName,dim1,dim2,dim3, &
        dim4,array, &
        subset,istart1,istart2,istart3,istart4,iend1,iend2,iend3,iend4)

      ! Arguments
      integer,intent(in) :: fileFormat
      integer,intent(in) :: fileID
      character(len=*),intent(in) :: variableName
      integer,intent(out) :: dim1,dim2,dim3,dim4
      real, allocatable, intent(out) :: array(:,:,:,:)
      logical,intent(in),optional :: subset
      integer,intent(in),optional :: istart1,istart2,istart3,istart4
      integer,intent(in),optional :: iend1,iend2,iend3,iend4

      ! Read array
      select case (fileFormat)
      case (HDF4_FORMAT, HDFEOS2_FORMAT)
         call readHdf4RealArray4d(fileID,variableName,dim1,dim2,dim3,dim4, &
              array, &
              subset,istart1,istart2,istart3,istart4, &
              iend1,iend2,iend3,iend4)
      case (NETCDF_FORMAT)
         call readNetcdfRealArray4d(fileID,variableName,dim1,dim2,dim3,dim4, &
              array, &
              subset,istart1,istart2,istart3,istart4, &
                 iend1,iend2,iend3,iend4)
      case default
         call printError(fileFormat)
         stop 1
      end select

      return
   end subroutine readRealArray4d

end module FileUtils_mod
