!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
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
! 17 Apr 2012 - Added option to open WPS file without overwriting.  Added
!               new subroutine to read 2D slab from file and handle 
!               end-of-file.
! 16 Aug 2016 - Now uses FileUtil_mod for PATH_MAX parameter, specifying
!               maximum length of file path.  Changed ERROR to FATAL in
!               error messages where the program halts.
!------------------------------------------------------------------------------

module FileWPS_mod

  ! Use external modules
  use FieldWPS_mod, only: FieldWPS, &
       checkFieldWPS, destroyFieldWPS, &
       IPROJ_LATLON, IPROJ_MERCATOR, IPROJ_LAMBERT_CONFORMAL, &
       IPROJ_GAUSSIAN, IPROJ_POLAR_STEREOGRAPHIC
  use FileUnit_mod, only: PATH_MAX
  
  ! Force explicit variable declarations
  implicit none
 
  ! Force explicit public declarations
  private

  ! FileWPS data type
  public :: FileWPS
  type FileWPS
     integer :: fileUnit
     character(len=PATH_MAX) :: fileName
     logical :: isOpen
     character(len=24) :: hdate
     logical :: isOld
  end type FileWPS

  ! Public methods
  public :: createFileWPS
  public :: destroyFileWPS
  public :: writeFileWPS
  public :: readFileWPS

  ! Internal parameters
  integer,parameter :: MAX_PREFIX = PATH_MAX - 14
  integer,parameter :: FILE_WPS_MISSING = -9999

contains

  !----------------------------------------------------------------------------
  !
  ! ROUTINE:  createFileWPS
  !
  ! DESCRIPTION:  Public "constructor method" for FileWPS data type.
  !
  !----------------------------------------------------------------------------

  function createFileWPS(fileUnit,outputDirectory,prefix, &
       year,month,day,hour,preserve) result (this)

    ! Arguments
    integer,intent(in) :: fileUnit
    character(len=*),intent(in) :: outputDirectory
    character(len=*),intent(in) :: prefix
    integer,intent(in) :: year
    integer,intent(in) :: month
    integer,intent(in) :: day
    integer,intent(in) :: hour
    logical,optional,intent(in) :: preserve

    ! Return variable
    type(FileWPS) :: this

    ! Local variables
    logical :: oldFile
    integer :: status
    integer :: length

    ! See if the open statement should have status=old
    oldFile=.false.
    if (present(preserve)) then
       oldFile=preserve
    end if
    this%isOld = oldFile

    ! Sanity check the prefix length
    if (len(prefix) > MAX_PREFIX) then
       print*,'FATAL, prefix must be less than ',MAX_PREFIX,' characters!'
       print*,'Current length is ',len(prefix)
       stop 1
    end if

    ! Build full file name
    write(this%filename,'(A,A,I4.4,A,I2.2,A,I2.2,A,I2.2)') &
         trim(prefix),':',year,'-',month,'-',day,'_',hour
    length = len(trim(outputDirectory)) + len(trim(this%filename))
    if (length > PATH_MAX) then
       print*,'FATAL, output file name is too long, must not exceed ', &
            PATH_MAX, 'characters!'
       stop 1
    end if
    
    this%filename = trim(outputDirectory) // trim(this%filename)
    
    ! Set up data structure and open file
    if (oldFile) then
       print*,'opening old WPS file: ',trim(this%fileName)
       open(unit=fileUnit,file=trim(this%fileName),iostat=status, &
            form="unformatted",convert="big_endian",status='old')
       rewind(fileUnit)
    else
       print*,'opening new WPS file: ',trim(this%fileName)
       open(unit=fileUnit,file=trim(this%fileName),iostat=status, &
            form="unformatted",convert="big_endian")
    end if
    if (status /= 0) then
       print*,'FATAL, cannot open ',trim(this%fileName)
       stop 1
    end if

    this%fileUnit = fileUnit
    this%isOpen = .true.

    write(this%hdate,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A)') &
         year,":",month,":",day,"_",hour,":00:00"

    return
  end function createFileWPS

  !----------------------------------------------------------------------------
  !
  ! ROUTINE:  destroyFileWPS
  !
  ! DESCRIPTION:  Public "destructor method" for FileWPS data type.
  !
  !----------------------------------------------------------------------------

  subroutine destroyFileWPS(this)

    ! Arguments
    type(FileWPS),intent(inout) :: this

    ! Local variables
    integer :: status

    ! Clean up
    if (this%isOpen) then
       close(this%fileUnit,iostat=status)
       if (status /= 0) then
          print*,'FATAL, cannot close ',trim(this%fileName)
          stop 1
       end if
    end if
    this%isOpen = .false.
    this%isOld = .false.
    this%fileUnit = FILE_WPS_MISSING
    this%fileName = "NULL"
    this%hdate = "NULL"

    return
  end subroutine destroyFileWPS

  !----------------------------------------------------------------------------
  !
  ! ROUTINE:  writeFileWPS
  !
  ! DESCRIPTION:  Public method for writeing data from FieldWPS to file 
  ! associated with FileWPS.
  !
  !----------------------------------------------------------------------------

  subroutine writeFileWPS(this,field)

    ! Arguments
    type(FileWPS),intent(in) :: this
    type(FieldWPS),intent(in) :: field

    ! Local variables
    integer :: status

    ! Sanity check the field
    call checkFieldWPS(field)

    ! Make sure file is open
    if (.not. this%isOpen) then
       print*,'FATAL, file ',trim(this%fileName),' is not open for writing!'
       stop 1
    end if

    ! Make sure current field is for current file
    if (trim(this%hdate) /= trim(field%hdate)) then
       print*,'FATAL, field date does not match file date!'
       print*,'field%hdate = ',trim(field%hdate)
       print*,'file%hdate = ',trim(this%hdate)
       stop 1
    end if

    !print*,'Writing field ',trim(field%field),' to file ',trim(this%fileName)
    write(unit=this%fileUnit,iostat=status) field%version
    if (status /= 0) then
       print*,'FATAL, cannot write to WPS intermediate format file!'
       print*,'status = ',status
       stop 1
    endif

    ! Write metadata
    write(unit=this%fileUnit,iostat=status) field%hdate,field%xfcst, &
         field%map_source,field%field,field%units,field%desc,field%xlvl, &
         field%nx, field%ny,field%iproj
    if (status /= 0) then
       print*,'FATAL, cannot write to WPS intermediate format file!'
       stop 1
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
       print*,'FATAL, invalid iproj value!'
       print*,'field%iproj = ',field%iproj
       stop 1
    end if
    if (status /= 0) then
       print*,'FATAL, cannot write to WPS intermediate format file!'
       stop 1
    endif

    ! Write wind rotation flag
    write(unit=this%fileUnit,iostat=status) field%is_wind_grid_rel
    if (status /= 0) then
       print*,'FATAL, cannot write to WPS intermediate format file!'
       stop 1
    endif

    ! Write 2-d array of data
    write(unit=this%fileUnit,iostat=status) field%slab
    if (status /= 0) then
       print*,'FATAL, cannot write to WPS intermediate format file!'
       stop 1
    endif

    return
  end subroutine writeFileWPS

  !----------------------------------------------------------------------------
  !
  ! ROUTINE:  readFileWPS
  !
  ! DESCRIPTION:  Public method for reading data file associated with FileWPS 
  ! to FieldWPS.
  !
  !----------------------------------------------------------------------------

  subroutine readFileWPS(this,field,endOfFile) 

     ! Arguments
     type(FileWPS), intent(in) :: this
     type(FieldWPS), intent(out) :: field
     logical, intent(out) :: endOfFile

     ! Local variables
     integer :: status

     ! Initialize field
     call destroyFieldWPS(field)

     ! Make sure file is open for reading
     if (.not. this%isOpen .or. .not. this%isOld) then
        print*,'FATAL, file ',trim(this%fileName),' is not open for reading!'
        stop 1
     end if

     ! Read version
     endOfFile = .false.
     read(unit=this%fileUnit,iostat=status,end=100) field%version
     if (status .ne. 0) then                                                 
        print*,'FATAL, cannot read WPS intermediate file!'
        stop 1
     end if
     
     ! Read header
     read(unit=this%fileUnit,iostat=status) field%hdate,field%xfcst, &
          field%map_source,field%field,field%units,field%desc,field%xlvl, &
          field%nx,field%ny,field%iproj
      if (status .ne. 0) then
         print*,'FATAL, cannot read WPS intermediate file!'
         stop 1
      end if

      ! Read map project data
      if (field%iproj == IPROJ_LATLON) then
         ! Cylindrical equidistant (lat/lon)    
         read(unit=this%fileUnit,iostat=status) field%startloc, &
              field%startlat, field%startlon, field%deltalat, field%deltalon, &
              field%earth_radius
      else if (field%iproj == IPROJ_MERCATOR) then
         ! Mercator                           
         read(unit=this%fileUnit,iostat=status) field%startloc, &
              field%startlat, field%startlon, field%dx, field%dy, &
              field%truelat1, field%earth_radius
      else if (field%iproj == IPROJ_LAMBERT_CONFORMAL) then
         ! Lambert conformal                   
         read(unit=this%fileUnit,iostat=status) field%startloc, &
              field%startlat, field%startlon, field%dx, field%dy, &
              field%xlonc, field%truelat1, field%truelat2, field%earth_radius
      else if (field%iproj == IPROJ_GAUSSIAN) then
         ! Gaussian
         read(unit=this%fileUnit,iostat=status)  field%startloc, &
              field%startlat, field%startlon, field%nlats, field%deltalon, &
              field%earth_radius
      else if (field%iproj == IPROJ_POLAR_STEREOGRAPHIC) then
         ! Polar stereographic                      
         read(unit=this%fileUnit,iostat=status) field%startloc, &
              field%startlat, field%startlon, field%xlonc, field%truelat1, &
              field%earth_radius
      else
         print*,'FATAL, invalid iproj value!'
         print*,'field%iproj = ',field%iproj
         stop 1
      end if
      if (status .ne. 0) then
         print*,'FATAL, cannot read WPS intermediate file!'
         stop 1
      end if

      ! Read wind rotation flag
      read(unit=this%fileUnit,iostat=status) field%is_wind_grid_rel
      if (status .ne. 0) then
         print*,'FATAL, cannot read WPS intermediate file!'
         stop 1
      end if

      ! Read 2-d array of data
      allocate(field%slab(field%nx,field%ny))
      read(unit=this%fileUnit,iostat=status) field%slab
      if (status .ne. 0) then
         print*,'FATAL, cannot read WPS intermediate file!'
         stop 1
      end if

      ! Sanity check the field
      call checkFieldWPS(field)

      return

      ! Handle end of file
      100 continue
      endOfFile = .true.
      return

   end subroutine readFileWPS
end module FileWPS_mod
