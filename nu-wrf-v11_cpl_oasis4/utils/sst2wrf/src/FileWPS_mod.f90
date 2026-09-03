!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  FileWPS_mod
!
! AUTHOR:
! Eric Kemp, NGIS 
!
! DESCRIPTION:
! Defines FileWPS data type and associated routines, for writing data in
! WPS intermediate format.
!
! REVISION HISTORY:
! 07 Apr 2010 - Initial version
!
!------------------------------------------------------------------------------

module FileWPS_mod

  ! Use external modules
  use FieldWPS_mod, only: FieldWPS, &
       checkFieldWPS, &
       IPROJ_LATLON, IPROJ_MERCATOR, IPROJ_LAMBERT_CONFORMAL, &
       IPROJ_GAUSSIAN, IPROJ_POLAR_STEREOGRAPHIC

  ! Force explicit variable declarations
  implicit none
 
  ! Force explicit public declarations
  private

  ! FileWPS data type
  public :: FileWPS
  type FileWPS
     integer :: fileUnit
     character(len=150) :: fileName
     logical :: isOpen
     character(len=24) :: hdate
  end type FileWPS

  ! Public methods
  public :: createFileWPS
  public :: destroyFileWPS
  public :: writeFileWPS

  ! Internal parameters
  integer,parameter :: MAX_FILENAME = 132
  integer,parameter :: MAX_PREFIX = MAX_FILENAME - 14
  integer,parameter :: FILE_WPS_MISSING = -9999

contains

  !----------------------------------------------------------------------------
  ! ROUTINE:
  ! createFileWPS
  !
  ! DESCRIPTION:
  ! Public "constructor method" for FileWPS data type.
  !----------------------------------------------------------------------------

  function createFileWPS(fileUnit,outputDirectory,prefix, &
       year,month,day,hour) result (this)

    ! Force explicit variable declarations
    implicit none

    ! Arguments
    integer,intent(in) :: fileUnit
    character(len=*),intent(in) :: outputDirectory
    character(len=*),intent(in) :: prefix
    integer,intent(in) :: year
    integer,intent(in) :: month
    integer,intent(in) :: day
    integer,intent(in) :: hour

    ! Return variable
    type(FileWPS) :: this

    ! Local variables
    integer :: status
    integer :: length

    ! Sanity check the prefix length
    if (len(prefix) > MAX_PREFIX) then
       print*,'ERROR, prefix must be less than ',MAX_PREFIX,' characters!'
       print*,'Current length is ',len(prefix)
       stop
    end if

    ! Build full file name
    write(this%filename,'(A,A,I4.4,A,I2.2,A,I2.2,A,I2.2)') &
         trim(prefix),':',year,'-',month,'-',day,'_',hour
    length = len(trim(outputDirectory)) + len("/") + len(trim(this%filename))
    if (length > MAX_FILENAME) then
       print*,'ERROR, output file name is too long, must not exceed ', &
            MAX_FILENAME, 'characters!'
       stop
    end if
    
    this%filename = trim(outputDirectory) // "/" // trim(this%filename)
    
    ! Set up data structure and open file
    print*,'opening wps file: ',trim(this%fileName)
    open(unit=fileUnit,file=trim(this%fileName),iostat=status, &
         form="unformatted",convert="big_endian")
    if (status /= 0) then
       print*,'ERROR opening ',trim(this%fileName)
       stop
    end if

    this%fileUnit = fileUnit
    this%isOpen = .true.

    write(this%hdate,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A)') &
         year,":",month,":",day,"_",hour,":00:00"

    return
  end function createFileWPS

  !----------------------------------------------------------------------------
  ! ROUTINE:
  ! destroyFileWPS
  !
  ! DESCRIPTION:
  ! Public "destructor method" for FileWPS data type.
  !----------------------------------------------------------------------------

  subroutine destroyFileWPS(this)

    ! Force explicit variable declarations
    implicit none

    ! Arguments
    type(FileWPS),intent(inout) :: this

    ! Local variables
    integer :: status

    ! Clean up
    if (this%isOpen) then
       close(this%fileUnit,iostat=status)
       if (status /= 0) then
          print*,'ERROR closing file ',trim(this%fileName)
          stop
       end if
    end if
    this%isOpen = .false.
    this%fileUnit = FILE_WPS_MISSING
    this%fileName = "NULL"
    this%hdate = "NULL"

    return
  end subroutine destroyFileWPS

  !----------------------------------------------------------------------------
  ! ROUTINE:
  ! writeFileWPS
  !
  ! DESCRIPTION:
  ! Writes data from FieldWPS to file associated with FileWPS.
  !----------------------------------------------------------------------------

  subroutine writeFileWPS(this,field)

    ! Force explicit variable declarations
    implicit none

    ! Arguments
    type(FileWPS),intent(in) :: this
    type(FieldWPS),intent(in) :: field

    ! Local variables
    integer :: status

    ! Sanity check the field
    call checkFieldWPS(field)

    ! Make sure file is open
    if (.not. this%isOpen) then
       print*,'ERROR, file ',trim(this%fileName),' is not open for writing!'
       stop
    end if

    ! Make sure current field is for current file
    if (trim(this%hdate) /= trim(field%hdate)) then
       print*,'ERROR, field date does not match file date!'
       print*,'field%hdate = ',trim(field%hdate)
       print*,'file%hdate = ',trim(this%hdate)
       stop
    end if

    print*,'Writing field ',trim(field%field),' to file ',trim(this%fileName)
    write(unit=this%fileUnit,iostat=status) field%version
    if (status /= 0) then
       print*,'ERROR writing to WPS intermediate format file!'
       print*,'status = ',status
       stop
    endif

    ! Write metadata
    write(unit=this%fileUnit,iostat=status) field%hdate,field%xfcst, &
         field%map_source,field%field,field%units,field%desc,field%xlvl, &
         field%nx, field%ny,field%iproj
    if (status /= 0) then
       print*,'ERROR writing to WPS intermediate format file!'
       stop
    endif

    if (field%iproj == IPROJ_LATLON) then
     ! Cylindrical equidistant (lat/lon)
       write(unit=this%fileUnit,iostat=status) field%startloc, &
            field%startlat, field%startlon, field%deltalat, field%deltalon, &
            field%earth_radius
    else if (field%iproj == IPROJ_MERCATOR) then
       ! Mercator
       write(unit=this%fileUnit,iostat=status) field%startloc, &
            field%startlat, field%startlon, field%dx, field%dy, &
            field%truelat1, field%earth_radius
    else if (field%iproj == IPROJ_LAMBERT_CONFORMAL) then
       ! Lambert conformal
       write(unit=this%fileUnit,iostat=status) field%startloc, &
            field%startlat, field%startlon, field%dx, field%dy, field%xlonc, &
            field%truelat1, field%truelat2, field%earth_radius
    else if (field%iproj == IPROJ_GAUSSIAN) then
       ! Gaussian
       write(unit=this%fileUnit,iostat=status) field%startloc, &
            field%startlat, field%startlon, field%nlats, field%deltalon, &
            field%earth_radius
    else if (field%iproj == IPROJ_POLAR_STEREOGRAPHIC) then
       ! Polar stereographic
       write(unit=this%fileUnit,iostat=status) field%startloc, &
            field%startlat, field%startlon, field%xlonc, field%truelat1, &
            field%earth_radius
    else
       print*,'ERROR, invalid iproj value!'
       print*,'field%iproj = ',field%iproj
       stop
    end if
    if (status /= 0) then
       print*,'ERROR writing to WPS intermediate format file!'
       stop
    endif

    ! Write wind rotation flag
    write(unit=this%fileUnit,iostat=status) field%is_wind_grid_rel
    if (status /= 0) then
       print*,'ERROR writing to WPS intermediate format file!'
       stop
    endif

    ! Write 2-d array of data
    write(unit=this%fileUnit,iostat=status) field%slab
    if (status /= 0) then
       print*,'ERROR writing to WPS intermediate format file!'
       stop
    endif

    return
  end subroutine writeFileWPS

end module FileWPS_mod
