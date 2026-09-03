!------------------------------------------------------------------------------
! NASA/GSFC, Code 606
!------------------------------------------------------------------------------
!
! MODULE:  nml
!
! DESCRIPTION:
! Data and methods to read fortran namelists.
!
! NU-WRF preprocessors read different namelists each with different data
! types  (real, integer, logical, string) and formats (scalars, arrays).
! This module provides the mechanism to instantiate a namelist object (nml)
! and retrieve its associated data using a single method (read). Upon reading,
! the read method returns the field value in the desired data type/format.
!
!------------------------------------------------------------------------------
module nml_mod

   implicit none 
   private 

   public :: nml    
   public :: nml_new

   type nml
      integer :: file_unit
   contains
      procedure :: read_string
      procedure :: read_double
      procedure :: read_float
      procedure :: read_integer
      procedure :: read_logical
      procedure :: read_string_vector
      procedure :: read_double_vector
      procedure :: read_float_vector
      procedure :: read_integer_vector
      procedure :: read_logical_vector
      generic :: read => read_string, read_double, &
           read_float , read_integer, read_logical, &
           read_string_vector, read_double_vector, &
           read_float_vector, read_integer_vector, &
           read_logical_vector
      procedure :: print_string
      procedure :: print_double
      procedure :: print_float
      procedure :: print_integer
      procedure :: print_logical
      procedure :: print_string_vector
      procedure :: print_double_vector
      procedure :: print_float_vector
      procedure :: print_integer_vector
      procedure :: print_logical_vector
      generic :: print => print_string , print_double, &
           print_float, print_integer, print_logical, &
           print_string_vector, print_double_vector, &
           print_float_vector, print_integer_vector, &
           print_logical_vector
      procedure :: set_verbose
   end type nml

   logical :: VERBOSE = .false.
   integer, parameter :: MAX_NML_LEN = 5000
   ! Missing values
   integer, parameter :: MISSING_INT  = -9999
   real,    parameter :: MISSING_REAL = -9999.
   character(len=8), parameter :: MISSING_STR = 'MISSING'
   
contains

   function nml_new(file_name) result (this)
      use fileunit_mod, only: select_file_unit
      ! Return variable
      type(nml) :: this
      character(len=*), intent(in) :: file_name
      
      this%file_unit = select_file_unit()
      open(unit=this%file_unit, file=file_name, delim='APOSTROPHE')
      
   end function nml_new
   
   subroutine set_verbose(this, option)
      class (nml), intent(inout) :: this
      logical :: option
      VERBOSE = option
   end subroutine set_verbose

   ! =============================================================
   !
   ! nml parameter reading functions
   !
   ! =============================================================
   
   ! This is the basic nml reading subroutine
   ! All interfaces use this to read the parameter, then it
   ! is converted to the correct type
   subroutine read_internal(file_unit,group,name,value,comment)
      integer, intent(IN)             :: file_unit
      character(len=*), intent(IN)    :: group, name 
      character(len=*), intent(INOUT) :: value 
      character(len=*), intent(INOUT), optional :: comment   
      integer :: iostat, l, ltype 
      character(len=1000) :: line, name1, value1, comment1 
      logical :: ingroup

      ingroup = .FALSE. 

      rewind(file_unit)
      do l = 1, MAX_NML_LEN
         read(file_unit,"(a1000)", iostat=iostat) line 
         if (iostat /= 0) exit 
         call parse_line(line, ltype, name1, value1, comment1 )

         ! Check if the parameter has been found
         if (ingroup .and. ltype == 3) then 
            if (trim(name1) == trim(name)) then 
               value   = trim(value1)
               comment = trim(comment1)
               exit 
            end if
         end if

         ! Open and close group as necessary
         if (ltype == 1 .and. trim(name1) == trim(group)) ingroup = .TRUE. 
         if (ltype == 2)                                  ingroup = .FALSE. 
         if (l == MAX_NML_LEN) &
              write(*,*) "ERROR: max nml length reached."
      end do

   end subroutine read_internal

   subroutine read_string(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      character(len=*), intent(INOUT) :: value 
      character(len=*), intent(IN)    :: group, name 
      character(len=*), intent(INOUT), optional :: comment 
      logical, optional :: init  
      logical :: init_var 
      character(len=256) :: value_str 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value = MISSING_STR

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)

      if (value_str /= "") value = trim(value_str)

      if (VERBOSE) &
           call this%print_string(name,value,comment)  

   end subroutine read_string

   subroutine read_double(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      double precision, intent(INOUT) :: value 
      character(len=*), intent(IN)    ::  group, name 
      character(len=*), intent(INOUT), optional :: comment   
      logical, optional :: init  
      logical :: init_var 
      character(len=256) :: value_str 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value = MISSING_REAL

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)
      if (value_str /= "") value = string_to_double(value_str)

      if (VERBOSE) &
           call this%print_string(name,value_str,comment)  

   end subroutine read_double

   subroutine read_float(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      real(4), intent(INOUT) :: value 
      character(len=*), intent(IN)    ::  group, name 
      character(len=*), intent(INOUT), optional :: comment  
      logical, optional :: init  
      logical :: init_var  
      character(len=256) :: value_str 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value = MISSING_REAL

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)
      if (value_str /= "") value = real(string_to_double(value_str))

      if (VERBOSE) &
           call this%print_string(name,value_str,comment)  

   end subroutine read_float

   subroutine read_integer(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      integer, intent(INOUT) :: value 
      character(len=*), intent(IN)    ::  group, name 
      character(len=*), intent(INOUT), optional :: comment 
      logical, optional :: init  
      logical :: init_var   
      character(len=256) :: value_str 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value = MISSING_INT

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)

      if (value_str /= "") then 
         value = nint(string_to_double(value_str))
      end if

      if (VERBOSE) &
           call this%print_string(name,value_str,comment)  

   end subroutine read_integer

   subroutine read_logical(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      logical, intent(INOUT) :: value 
      character(len=*), intent(IN)    ::  group, name 
      character(len=*), intent(INOUT), optional :: comment 
      logical, optional :: init  
      logical :: init_var   
      character(len=256) :: value_str 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value = .FALSE. 

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)
      if (value_str /= "") value = string_to_logical(value_str)

      if (VERBOSE) &
           call this%print_string(name,value_str,comment)  

   end subroutine read_logical

   !! Vectors 

   subroutine read_string_vector(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      character(len=*), intent(INOUT) :: value(:) 
      character(len=*), intent(IN)    ::  group, name 
      character(len=*), intent(INOUT), optional :: comment 
      logical, optional :: init  
      logical :: init_var  
      character(len=256) :: value_str 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value(:) = "" 

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)

      if (value_str /= "") call string_to_vector(value_str,value)

      if (VERBOSE) &
           call this%print_string(name,value_str,comment)  

   end subroutine read_string_vector

   subroutine read_double_vector(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      double precision, intent(INOUT) :: value(:) 
      character(len=*), intent(IN)    ::  group, name 
      character(len=*), intent(INOUT), optional :: comment   
      character(len=256) :: value_str, value_str_vec(size(value))
      logical, optional :: init  
      logical :: init_var  
      integer :: q 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value(:) = 0.d0

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)

      if (value_str /= "") then
         call string_to_vector(value_str,value_str_vec)
         do q = 1, size(value)
            if (trim(value_str_vec(q)) /= "") then 
               value(q) = string_to_double(trim(adjustl(value_str_vec(q))))
            end if

         end do
      end if

      if (VERBOSE) &
           call this%print_string(name,value_str,comment)  

   end subroutine read_double_vector

   subroutine read_float_vector(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      real(4), intent(INOUT) :: value(:) 
      character(len=*), intent(IN)    ::  group, name 
      character(len=*), intent(INOUT), optional :: comment   
      character(len=256) :: value_str, value_str_vec(size(value)) 
      logical, optional :: init  
      logical :: init_var
      integer :: q 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value(:) = 0.0

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)

      if (value_str /= "") then
         call string_to_vector(value_str,value_str_vec)
         do q = 1, size(value)
            if (trim(value_str_vec(q)) /= "") then 
               value(q) = &
                    real(string_to_double(trim(adjustl(value_str_vec(q)))))
            end if

         end do
      end if

      if (VERBOSE) &
           call this%print_string(name,value_str,comment)  

   end subroutine read_float_vector

   subroutine read_integer_vector(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      integer, intent(INOUT) :: value(:) 
      character(len=*), intent(IN)    ::  group, name 
      character(len=*), intent(INOUT), optional :: comment   
      character(len=256) :: value_str, value_str_vec(size(value))
      logical, optional :: init  
      logical :: init_var 
      integer :: q 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value(:) = 0

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)

      if (value_str /= "") then
         call string_to_vector(value_str,value_str_vec)
         do q = 1, size(value)
            if (trim(value_str_vec(q)) /= "") then 
               value(q) = &
                    nint(string_to_double(trim(adjustl(value_str_vec(q)))))
            end if

         end do
      end if

      if (VERBOSE) &
           call this%print_string(name,value_str,comment)  

   end subroutine read_integer_vector

   subroutine read_logical_vector(this,group,name,value,comment,init)
      class (nml), intent(inout) :: this
      logical, intent(INOUT) :: value(:) 
      character(len=*), intent(IN)    ::  group, name 
      character(len=*), intent(INOUT), optional :: comment   
      character(len=256) :: value_str, value_str_vec(size(value))
      logical, optional :: init  
      logical :: init_var 
      integer :: q 

      init_var = .FALSE. 
      if (present(init)) init_var = init 
      if (init_var) value(:) = .FALSE.

      ! First find parameter value as a string 
      value_str = ""
      call read_internal(this%file_unit,group,name,value_str,comment)

      if (value_str /= "") then
         call string_to_vector(value_str,value_str_vec)
         do q = 1, size(value)
            if (trim(value_str_vec(q)) /= "") then 
               value(q) = string_to_logical(trim(adjustl(value_str_vec(q))))
            end if

         end do
      end if

      if (VERBOSE) &
           call this%print_string(name,value_str,comment)  

   end subroutine read_logical_vector

   ! =============================================================
   !
   ! nml line printing functions
   !
   ! =============================================================
   
   ! This is the basic routine for printing a parameter to a formatted line
   ! All other interfaces use this routine after converting to a string.
   subroutine print_string(this, name,value,comment,io,no_quotes)
      class (nml), intent(inout) :: this
      character(len=*), intent(in) :: name, value 
      character(len=*), intent(in), optional :: comment 
      integer, intent(in), optional :: io 
      integer :: io_val 
      logical, intent(in), optional :: no_quotes
      logical :: no_quotes_loc
      character(len=1000) :: line
      character(len=500)  :: comment1 
      character(len=len(value))  :: val_repr

      io_val = 6 
      if (present(io)) io_val = io 
      val_repr = value

      no_quotes_loc = .false.
      if (present(no_quotes)) no_quotes_loc = no_quotes
      if (.not.no_quotes_loc) val_repr = '"'//trim(value)//'"'

      comment1 = "" 
      if (present(comment)) comment1 = "   "//trim(comment)

      write(line,"(a)") &
           "    "//trim(name)//" = "//trim(val_repr)//trim(comment1)
      write(io_val,*) trim(line)

   end subroutine print_string

   subroutine print_double(this, name,value,comment,io)
      class (nml), intent(inout) :: this
      double precision :: value
      character(len=*) :: name 
      character(len=*), optional :: comment
      integer, optional :: io 
      character(len=500) :: value_str  

      write(value_str,*) value 
      if (VERBOSE) &
           call this%print_string(name,value_str,comment,io,no_quotes=.true.)

   end subroutine print_double

   subroutine print_float(this, name,value,comment,io)
      class (nml), intent(inout) :: this
      real(4) :: value
      character(len=*) :: name 
      character(len=*), optional :: comment
      integer, optional :: io 
      character(len=500) :: value_str  

      write(value_str,*) value 
      if (VERBOSE) &
           call this%print_string(name,value_str,comment,io)

   end subroutine print_float

   subroutine print_integer(this, name,value,comment,io)
      class (nml), intent(inout) :: this
      integer :: value
      character(len=*) :: name 
      character(len=*), optional :: comment
      integer, optional :: io 
      character(len=500) :: value_str  

      write(value_str,*) value 
      if (VERBOSE) &
           call this%print_string(name,value_str,comment,io,no_quotes=.true.)

   end subroutine print_integer

   subroutine print_logical(this, name,value,comment,io)
      class (nml), intent(inout) :: this
      logical :: value
      character(len=*) :: name 
      character(len=*), optional :: comment
      integer, optional :: io 
      character(len=500) :: value_str  

      value_str = "F"
      if (value) value_str = "T" 
      if (VERBOSE) &
           call this%print_string(name,value_str,comment,io,no_quotes=.true.)

   end subroutine print_logical

   !! Vectors

   subroutine print_string_vector(this, name,value,comment,io)
      class (nml), intent(inout) :: this
      character(len=*) :: value(:)
      character(len=*) :: name 
      character(len=*), optional :: comment
      integer, optional :: io 
      character(len=500) :: value_str  
      integer :: q 

      value_str = '"'//trim(value(1))//'"'
      do q = 2, size(value)
         write(value_str,*) trim(value_str)//" "//'"'//trim(value(q))//'"'
      end do

      if (VERBOSE) &
           call this%print_string(name,value_str,comment,io,no_quotes=.true.)

   end subroutine print_string_vector

   subroutine print_double_vector(this, name,value,comment,io)
      class (nml), intent(inout) :: this
      double precision :: value(:)
      character(len=*) :: name 
      character(len=*), optional :: comment
      integer, optional :: io 
      character(len=500) :: value_str  
      integer :: q 

      value_str = ""
      do q = 1, size(value)
         write(value_str,"(a,g12.3)") trim(value_str)//" ",value(q)
      end do

      if (VERBOSE) &
           call this%print_string(name,value_str,comment,io,no_quotes=.true.)

   end subroutine print_double_vector

   subroutine print_float_vector(this, name,value,comment,io)
      class (nml), intent(inout) :: this
      real(4) :: value(:)
      character(len=*) :: name 
      character(len=*), optional :: comment
      integer, optional :: io 
      character(len=500) :: value_str  
      integer :: q 

      value_str = ""
      do q = 1, size(value)
         write(value_str,"(a,g12.3)") trim(value_str)//" ",value(q)
      end do

      if (VERBOSE) &
           call this%print_string(name,value_str,comment,io,no_quotes=.true.)
   end subroutine print_float_vector

   subroutine print_integer_vector(this, name,value,comment,io)
      class (nml), intent(inout) :: this
      integer :: value(:)
      character(len=*) :: name 
      character(len=*), optional :: comment
      integer, optional :: io 
      character(len=500) :: value_str  
      integer :: q 

      value_str = ""
      do q = 1, size(value)
         write(value_str,"(a,i12)") trim(value_str)//" ",value(q)
      end do

      if (VERBOSE) &
           call this%print_string(name,value_str,comment,io,no_quotes=.true.)

   end subroutine print_integer_vector

   subroutine print_logical_vector(this, name,value,comment,io)
      class (nml), intent(inout) :: this
      logical :: value(:)
      character(len=*) :: name 
      character(len=*), optional :: comment
      integer, optional :: io 
      character(len=500) :: value_str  
      integer :: q 

      value_str = ""
      do q = 1, size(value)
         if (value(q)) then 
            write(value_str,"(a,a1)") trim(value_str)//" ","T"
         else
            write(value_str,"(a,a1)") trim(value_str)//" ","F"
         end if
      end do

      if (VERBOSE) &
           call this%print_string(name,value_str,comment,io,no_quotes=.true.)
   end subroutine print_logical_vector


   ! =============================================================
   !
   ! Type conversion functions
   !
   ! =============================================================

   function string_to_double(string) result(value)
      character(len=*), intent(IN) :: string 
      double precision :: value 

      character(len=256) :: tmpstr 
      integer :: stat, n
      double precision :: x 

      tmpstr = trim(adjustl(string))
      n      = len_trim(tmpstr)

      read(tmpstr(1:n),*,IOSTAT=stat) x

      value = 0
      if (stat .eq. 0) then 
         value = x 
      else
         n = len_trim(tmpstr)-1
         READ(tmpstr(1:n),*,IOSTAT=stat) x
         if (stat .ne. 0) then 
            write(*,*) "nml:: ","Error converting string to number!"
            write(*,*) "|",trim(tmpstr),"|",n,stat,x
         else
            value = x 
         end if
      end if

   end function string_to_double

   function string_to_logical(string) result(value)
      character(len=*), intent(IN) :: string 
      logical :: value 

      character(len=256) :: tmpstr 
      integer :: stat, n
      double precision :: x 

      tmpstr = trim(adjustl(string))

      select case(trim(tmpstr))
      case("T","True","TRUE","true",".TRUE.", ".true.")
         value = .TRUE. 
      case("F","False","FALSE","false",".FALSE.", ".false.")
         value = .FALSE. 
      case DEFAULT
         write(*,*) "nml:: Error reading logical parameter."
         stop 
      end select

   end function string_to_logical

   subroutine string_to_vector(string,value)
      character(len=*), intent(IN) :: string 
      character(len=*) :: value(:)
      character(len=256) :: tmpvec(size(value))
      character(len=256) :: tmpstr, fmt 
      integer :: stat, n, q, q1, q2, j 

      tmpstr = trim(adjustl(string))
      n      = len_trim(tmpstr)+2

      tmpvec(:) = "" 

      q1 = 1 
      do q = 1, size(tmpvec)
         q2 = index(tmpstr(q1:n)," ") + q1
         if (q2 .gt. q1 .and. q2 .le. n) then 
            tmpvec(q) = tmpstr(q1:q2-1)
            q1 = q2

            ! Make sure gaps of more than one space are properly handled
            do j = 1, 1000
               if (tmpstr(q1:q1) == " ") q1 = q1+1
               if (q1 .ge. n) exit 
            end do

            ! Remove quotes around string if they exist 
            call remove_quotes_comma(tmpvec(q))

         end if
      end do

      value = tmpvec 

   end subroutine string_to_vector

   ! =============================================================
   !
   ! Helper functions
   !
   ! =============================================================
   
   ! This function parses a namelist line into pieces, determining
   ! whether it is a blank line (-2), comment (-1), group name (1)
   ! end-of-group (2), or a parameter line (3)
   subroutine parse_line(line,linetype,name,value,comment)
      character(len=*), intent(IN)    :: line
      character(len=*), intent(INOUT) :: name, value, comment 
      integer :: linetype

      character(len=1000) :: line1 
      integer :: q, q1, q2

      name     = ""
      value    = ""
      comment  = "" 

      line1 = trim(adjustl(line))

      if (trim(line1) == "") then         ! Blank line
         linetype = -2 

      else if (line1(1:1) == "!") then    ! Comment line 
         linetype = -1 
         comment = trim(line1)

      else if (line1(1:1) == "&") then    ! Group name 
         linetype = 1 
         q = len_trim(line1)
         name     = line1(2:q)

      else if (line1(1:1) == "/") then    ! End of group 
         linetype = 2 

      else   ! Line must contain parameter to read
         linetype = 3

         q = index(line1,"=")
         if (q == 0) then 
            write(*,*) "nml:: Error reading namelist file."
            write(*,*) "No '=' found on parameter line."
            stop 
         end if

         name = trim(adjustl(line1(1:q-1)))

         q1 = index(line1,"!")
         q2 = len_trim(line1)

         if (q1 > 0) then 
            comment = trim(adjustl(line1(q1:q2)))
            value   = trim(adjustl(line1(q+1:q1-1)))
         else
            value   = trim(adjustl(line1(q+1:q2)))
         end if

         ! Remove quotes around string, and final line comma, if they exist
         call remove_quotes_comma(value)

      end if

   end subroutine parse_line

   subroutine remove_quotes_comma(string)
      character(len=*), intent(INOUT) :: string 
      integer :: i, n 

      ! Eliminate quotes
      n = len_trim(string)
      do i = 1,n 
         if (string(i:i) == '"' .or. string(i:i) == "'") string(i:i) = " "
      end do
      string = trim(adjustl(string))

      ! Remove final comma too
      n = len_trim(string)
      if (n > 0) then 
         if (string(n:n) == ",") string(n:n) = " "
         string = trim(adjustl(string))
      end if

   end subroutine remove_quotes_comma

end module nml_mod

