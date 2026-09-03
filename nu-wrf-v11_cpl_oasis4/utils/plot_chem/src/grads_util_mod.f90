!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE: grads_util_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! Contains utility routines for reading GrADS data.
!
!------------------------------------------------------------------------------

module grads_util_mod

   ! Defaults
   implicit none
   private

   ! Public data type
   public :: gradsctl
   type gradsctl
      private
      character(len=132) :: dset
      real :: undef
      character(len=132) :: title
      integer :: xdef_dim
      integer :: ydef_dim
      integer :: zdef_dim
      integer :: tdef_dim
      integer :: vars_dim
      character(len=132),allocatable :: vars_name(:)
   end type gradsctl

   ! Public routines
   public :: read_grads_binary
   public :: create_gradsctl
   public :: get_dims_gradsctl
   public :: get_binfilename_gradsctl
   public :: get_varname_gradsctl
   public :: get_mapfilename_gradsctl

contains

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  read_grads_binary
   !
   ! DESCRIPTION:  Reads in data from GrADS binary file.
   !
   !--------------------------------------------------------------------------
   
   subroutine read_grads_binary(filename,nx,ny,nz,nt,nvar,grads_data)

      character(len=*), intent(in) :: filename
      integer,intent(in) :: nx,ny,nz,nt,nvar
      real,intent(inout) :: grads_data(nx,ny,nz,nt,nvar)

      integer :: iunit=20
      integer :: i,j,k,m,n
      integer :: nrec
      integer :: output_byte_size

      ! FIXME:  Make this portable to systems other than Discover
      inquire(iolength=output_byte_size) 4 ! 4-byte = 32-bit

      open(iunit,file=trim(filename),form="unformatted",access="direct", &
           status="old",recl=output_byte_size*nx*ny)

      nrec = 0
      do n = 1,nvar
         do m = 1,nt
            do k = 1,nz
               nrec = nrec+1
               read(iunit,rec=nrec) ((grads_data(i,j,k,m,n),i=1,nx),j=1,ny)

!                PRINT*, '*****************************'
!                PRINT*, 'k,m,n = ',k,m,n
!                PRINT*, 'max:', maxval(grads_data(:,:,k,m,n))
!                PRINT*, 'min:', minval(grads_data(:,:,k,m,n)) 
!                print*, 'true min:',minval(grads_data(:,:,k,m,n))
!                PRINT*, '*****************************'
            end do
         end do
      end do

      close(iunit)
      return
   end subroutine read_grads_binary

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  create_gradsctl
   !
   ! DESCRIPTION:  Reads in data from GrADS control file.
   !
   ! WARNING:  This is not a robust subroutine.  It is intended to read
   ! the control file from prep_chem_sources and extract specific information.
   !
   !--------------------------------------------------------------------------

   function create_gradsctl(filename) result(this)

      ! Arguments
      character(len=*),intent(in) :: filename

      ! Return variable
      type(gradsctl) :: this

      ! Local variables
      character(len=132) :: ctmp
      integer :: iunit
      integer :: i,ii

      iunit = 10
      open(unit=iunit,file=trim(filename),status='old')

      read(unit=iunit,fmt='(A)') ctmp
      if (ctmp(1:4) .ne. "dset") then
         print*,'ERROR reading control file ',trim(filename)
         print*,'Read ',ctmp
         stop
      end if
      this%dset = trim(ctmp(5:))

      read(unit=iunit,fmt='(A)') ctmp
      if (ctmp(1:7) .eq. 'OPTIONS') then
         read(unit=iunit,fmt='(A)') ctmp
      end if

      if (ctmp(1:5) .ne. "undef") then
         print*,'ERROR reading control file ',trim(filename)
         print*,'Read ',ctmp
         stop
      end if
      read(ctmp,'(6x,F5.5)') this%undef

      read(unit=iunit,fmt='(A)') ctmp
      if (ctmp(1:5) .ne. "title") then
         print*,'ERROR reading control file ',trim(filename)
         print*,'Read ',ctmp
         stop
      end if
      this%title = trim(ctmp(6:))

      read(unit=iunit,fmt='(A)') ctmp
      if (ctmp(1:4) .ne. "xdef") then
         print*,'Read ',ctmp
         print*,'ERROR reading control file ',trim(filename)
         stop
      end if
      read(ctmp,'(5x,i4)') this%xdef_dim

      read(unit=iunit,fmt='(A)') ctmp
      if (ctmp(1:4) .ne. "ydef") then
         print*,'ERROR reading control file ',trim(filename)
         print*,'Read ',ctmp
         stop
      end if
      read(ctmp,'(5x,i4)') this%ydef_dim

      read(unit=iunit,fmt='(A)') ctmp
      if (ctmp(1:4) .ne. "zdef") then
         print*,'ERROR reading control file ',trim(filename)
         print*,'Read ',ctmp
         stop
      end if
      read(ctmp,'(5x,i4)') this%zdef_dim

      read(unit=iunit,fmt='(A)') ctmp
      if (ctmp(1:4) .ne. "tdef") then
         print*,'ERROR reading control file ',trim(filename)
         print*,'Read ',ctmp
         stop
      end if
      read(ctmp,'(5x,i4)') this%tdef_dim

      read(unit=iunit,fmt='(A)') ctmp
      if (ctmp(1:4) .ne. "vars") then
         print*,'ERROR reading control file ',trim(filename)
         print*,'Read ',ctmp
         stop
      end if
      read(ctmp,'(5x,i4)') this%vars_dim

      allocate(this%vars_name(this%vars_dim))

      do i = 1, this%vars_dim
         read(unit=iunit,fmt='(A)') ctmp
         print*,'Read ',ctmp
         ctmp = adjustl(ctmp)
         this%vars_name(i) = trim(ctmp)
      end do

      close(unit=iunit)

      return
   end function create_gradsctl

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  get_dims_gradsctl
   !
   ! DESCRIPTION:  Retrieve dimensions from gradsctl structure.
   !
   !--------------------------------------------------------------------------

   subroutine get_dims_gradsctl(this,nx,ny,nz,nt,nvar)

      ! Arguments
      type(gradsctl),intent(in) :: this
      integer,intent(out) :: nx
      integer,intent(out) :: ny
      integer,intent(out) :: nz
      integer,intent(out) :: nt
      integer,intent(out) :: nvar

      ! Copy the data
      nx = this%xdef_dim
      ny = this%ydef_dim
      nz = this%zdef_dim
      nt = this%tdef_dim
      nvar = this%vars_dim

      return
   end subroutine get_dims_gradsctl

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  get_binfilename_gradsctl
   !
   ! DESCRIPTION:  Retrieve GrADS binary filename from gradsctl structure.
   !
   !--------------------------------------------------------------------------

   function get_binfilename_gradsctl(this) result(filename)

      ! Arguments
      type(gradsctl), intent(in) :: this

      ! Return variable
      character(len=132) :: filename

      ! Local variables
      integer :: i

      filename = adjustl(trim(this%dset))
      
      i = 1
      if (filename(1:1) == '^') then
         i = 2
         filename = trim(filename(i:))
      end if      

      return
   end function get_binfilename_gradsctl

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  get_varname_gradsctl
   !
   ! DESCRIPTION:  Retrieve GrADS variable name from gradsctl structure.
   !
   !--------------------------------------------------------------------------

   function get_varname_gradsctl(this,i) result (varname)

      ! Arguments
      type(gradsctl),intent(in) :: this 
      integer,intent(in) :: i

      ! Return variable
      character(len=132) :: varname

      ! Sanity check
      if (i < 1 .or. i > this%vars_dim) then
         print*,'ERROR, invalid variable number!'
         print*,'var number = ',i
         print*,'Total number of variables = ',this%vars_dim
         stop
      end if

      varname = trim(this%vars_name(i))
      
      return
   end function get_varname_gradsctl

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  get_mapfilename_gradsctl
   !
   ! DESCRIPTION:  Retrieve map filename from gradsctl structure.
   !
   !--------------------------------------------------------------------------

   function get_mapfilename_gradsctl(this) result(filename)

      ! Arguments
      type(gradsctl), intent(in) :: this

      ! Return variable
      character(len=132) :: filename

      ! Local variables
      integer :: length
      integer :: index
      integer :: i

      filename = adjustl(trim(this%dset))
      
      i = 1
      if (filename(1:1) == '^') then
         i = 2
         filename = trim(filename(i:))
      end if      

      length = len_trim(filename)
      index = length - 2
      filename(index:) = 'map'

      return
   end function get_mapfilename_gradsctl

end module grads_util_mod
