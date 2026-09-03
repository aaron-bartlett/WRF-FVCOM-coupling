!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration and Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE:  sorted_datetimes_mod
!
! AUTHOR:
! Eric Kemp, NASA SIVO/Northrop Grumman
!
! DESCRIPTION:
! Defines a data type and attendant methods for saving datetimes (in order)
! as well as their association with wrfinput and wrfbdy files.
!------------------------------------------------------------------------------

module sorted_datetimes_mod

   ! Reset defaults
   implicit none
   private

   ! Public parameters
   integer, parameter :: MAX_DATE_LEN = 19
   public :: MAX_DATE_LEN

   ! Public data type
   type sorted_datetimes
      integer :: num_datetimes
      integer :: num_wrf_domains
      character(len=MAX_DATE_LEN),allocatable :: datetimes(:)
      logical,allocatable :: for_wrfbdy(:)
      logical,allocatable :: for_wrfinput(:,:)
   end type sorted_datetimes
   public :: sorted_datetimes

   ! Public methods
   public :: create_sorted_datetimes
   public :: destroy_sorted_datetimes
   public :: add_to_sorted_datetimes
   public :: get_sorted_datetimes
   public :: compare_datetimes

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  destroy_sorted_datetimes
   !
   ! DESCRIPTION:  Public destructor method for sorted_datetimes data type.
   !
   !---------------------------------------------------------------------------

   subroutine destroy_sorted_datetimes(this)

      ! Arguments
      type(sorted_datetimes), intent(inout) :: this

      ! Destroy
      if (allocated(this%for_wrfinput)) deallocate(this%for_wrfinput)
      if (allocated(this%for_wrfbdy)) deallocate(this%for_wrfbdy)
      if (allocated(this%datetimes)) deallocate(this%datetimes)
      this%num_wrf_domains = 0
      this%num_datetimes = 0

   end subroutine destroy_sorted_datetimes

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  create_sorted_datetimes
   !
   ! DESCRIPTION:  Public constructor method for sorted_datetimes data type.
   !
   !---------------------------------------------------------------------------

   function create_sorted_datetimes(num_wrf_domains) result (this)

      ! Arguments
      integer,intent(in) :: num_wrf_domains

      ! Return variable
      type(sorted_datetimes) :: this

      ! Initialize
      call destroy_sorted_datetimes(this)

      ! Sanity check
      if (num_wrf_domains < 1) then
         print*,'ERROR, number of wrf domains must be positive!'
         print*,'num_wrf_domains = ',num_wrf_domains
         stop
      end if
      this%num_wrf_domains = num_wrf_domains
   end function create_sorted_datetimes
   
   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  add_to_sorted_datetimes
   !
   ! DESCRIPTION:  Public method for adding a datetime to the sorted_datetimes
   ! data type, plus information on whether a wrfbdy or wrfinput file is
   ! associated with the datetime.  This is a generic handler--the data type
   ! may be empty, the datetime may already be in the data type, or the 
   ! data type arrays may need to be reallocated to add the new datetime.
   !
   !---------------------------------------------------------------------------

   subroutine add_to_sorted_datetimes(this,datetime,for_wrfbdy,for_wrfinput, &
        wrf_domain_id)

      ! Arguments
      type(sorted_datetimes),intent(inout) :: this
      character(len=MAX_DATE_LEN),intent(in) :: datetime
      logical,intent(in) :: for_wrfbdy
      logical,intent(in) :: for_wrfinput
      integer,intent(in) :: wrf_domain_id

      ! Local variables
      logical :: found_match

      ! Sanity checks
      if (this%num_wrf_domains < 1) then
         print*,'ERROR, number of wrf domains must be positive'
         print*,'Is sorted_datetimes data type initialized?'
         print*,'this%num_wrf_domains = ',this%num_wrf_domains
         stop
      end if
      if ( .not. for_wrfbdy .and. .not. for_wrfinput) then
         print*,'ERROR, datetime is neither for wrfbdy nor wrfinput!'
         print*,'datetime = ',datetime
         stop
      end if
      if (.not. for_wrfinput .and. wrf_domain_id .ne. 1) then
         print*,'ERROR, wrf_domain_id must be 1 if datetime is for wrfbdy!'
         print*,'datetime = ',datetime
         print*,'wrf_domain_id = ',wrf_domain_id
         stop
      end if
      if (for_wrfinput) then
         if (wrf_domain_id < 1 .or. wrf_domain_id > this%num_wrf_domains) then
            print*,'ERROR, invalid wrf domain ID used!'
            print*,'Must be from 1 to ',this%num_wrf_domains
            print*,'wrf_domain_id = ',wrf_domain_id
            stop
         end if
      end if

      ! Easy case:  No sorted datetimes exist yet.  Just allocate the
      ! array for single value, assign, and return.
      if (.not. allocated(this%datetimes)) then
         call add_first_sorted_datetimes(this,datetime,for_wrfbdy, &
              for_wrfinput, wrf_domain_id)
         return
      end if
      
      ! Harder case:  Need to see if datetime already exists in sorted list.
      ! If so, update for_wrfbdy or for_wrfinput
      found_match = update_usage_sorted_datetimes(this,datetime,for_wrfbdy, &
        for_wrfinput, wrf_domain_id) 
      if (found_match) return

      ! Hard case:  Need to add datetime to existing list.  
      call insert_to_sorted_datetimes(this,datetime,for_wrfbdy, &
           for_wrfinput, wrf_domain_id)

   end subroutine add_to_sorted_datetimes

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  add_first_sorted_datetimes
   !
   ! DESCRIPTION:  Private method for adding the very first datetime to the
   ! data type.  Assumes data type is empty (no arrays allocated).
   !
   !---------------------------------------------------------------------------

   subroutine add_first_sorted_datetimes(this,datetime,for_wrfbdy, &
        for_wrfinput, wrf_domain_id)

      ! Arguments
      type(sorted_datetimes),intent(inout) :: this
      character(len=MAX_DATE_LEN),intent(in) :: datetime
      logical,intent(in) :: for_wrfbdy
      logical,intent(in) :: for_wrfinput
      integer,intent(in) :: wrf_domain_id

      ! Put the datetime into the currently empty data type.
      ! Calling routine is required to sanity check.
      this%num_datetimes = 1
      allocate(this%datetimes(1))
      allocate(this%for_wrfbdy(1))
      allocate(this%for_wrfinput(this%num_wrf_domains,1))
      this%datetimes(1) = datetime
      this%for_wrfbdy(:) = .false.
      if (for_wrfbdy) then
         this%for_wrfbdy(1) = for_wrfbdy
      end if
      this%for_wrfinput(:,:) = .false.
      if (for_wrfinput) then
         this%for_wrfinput(wrf_domain_id,1) = for_wrfinput
      end if
      
   end subroutine add_first_sorted_datetimes

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  update_usage_sorted_datetimes
   !
   ! DESCRIPTION:  Private method for determining if datetime is already in
   ! data type.  If so, update the usage (wrfbdy or wrfinput).  Return
   ! a logical indicating whether a match was found.
   !
   !---------------------------------------------------------------------------

   function update_usage_sorted_datetimes(this,datetime,for_wrfbdy, &
        for_wrfinput, wrf_domain_id) result(found_match)

      ! Arguments
      type(sorted_datetimes),intent(inout) :: this
      character(len=MAX_DATE_LEN),intent(in) :: datetime
      logical,intent(in) :: for_wrfbdy
      logical,intent(in) :: for_wrfinput
      integer,intent(in) :: wrf_domain_id

      ! Return variable
      logical :: found_match

      ! Local variables
      integer :: ii

      ! Loop through the data type to see if current datetime is already 
      ! present.  Sanity check is the responsibility of the calling routine.
      found_match = .false.
      do ii = 1, this%num_datetimes
         if (this%datetimes(ii) .eq. datetime) then
            found_match = .true.
            if (for_wrfbdy) then
               this%for_wrfbdy(ii) = for_wrfbdy
            end if
            if (for_wrfinput) then
               this%for_wrfinput(wrf_domain_id,ii) = for_wrfinput
            end if
            exit
         end if
      end do

   end function update_usage_sorted_datetimes

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  insert_to_sorted_datetimes
   !
   ! DESCRIPTION:  Private method for putting new datetime into data type.
   ! Requires reallocating the arrays in the data type to fit all the
   ! information.
   !
   !---------------------------------------------------------------------------

   subroutine insert_to_sorted_datetimes(this,datetime,for_wrfbdy, &
        for_wrfinput, wrf_domain_id)

      ! Arguments
      type(sorted_datetimes),intent(inout) :: this
      character(len=MAX_DATE_LEN),intent(in) :: datetime
      logical,intent(in) :: for_wrfbdy
      logical,intent(in) :: for_wrfinput
      integer,intent(in) :: wrf_domain_id

      ! Local variables
      character(len=MAX_DATE_LEN),allocatable :: tmp_datetimes(:)
      logical,allocatable :: tmp_for_wrfbdy(:)
      logical,allocatable :: tmp_for_wrfinput(:,:)
      integer :: order
      integer :: new_index
      integer :: ii,jj

      do ii = 1, this%num_datetimes
         order = compare_datetimes(this%datetimes(ii),datetime)
!         if (compare_datetimes(datetime,this%datetimes(ii)) == -1) then
         if (order == -1) then
            new_index = ii
            exit
!         else if (compare_datetimes(datetime,this%datetimes(ii)) == 0) then
         else if (order == 0) then
            print*,'ERROR, datetime is already in data type!'
            print*,'datetime = ',datetime
            stop
         else
            new_index = ii+1
         end if
      end do

      ! Need some temporary space
      allocate(tmp_datetimes(this%num_datetimes+1))
      allocate(tmp_for_wrfbdy(this%num_datetimes+1))
      allocate(tmp_for_wrfinput(this%num_wrf_domains,this%num_datetimes+1))

      ! Copy all datetime info before new datetime
      do jj = 1, new_index-1
         tmp_datetimes(jj) = this%datetimes(jj)
         tmp_for_wrfbdy(jj) = this%for_wrfbdy(jj)
         do ii = 1, this%num_wrf_domains
            tmp_for_wrfinput(ii,jj) = this%for_wrfinput(ii,jj)
         end do
      end do

      ! Copy new datetime
      tmp_datetimes(new_index) = datetime
      tmp_for_wrfbdy(new_index) = for_wrfbdy
      tmp_for_wrfinput(:,new_index) = .false.
      tmp_for_wrfinput(wrf_domain_id,new_index) = for_wrfinput

      ! Copy datetime info after new datetime
      do jj = new_index,this%num_datetimes
         tmp_datetimes(jj+1) = this%datetimes(jj)
         tmp_for_wrfbdy(jj+1) = this%for_wrfbdy(jj)
         do ii = 1, this%num_wrf_domains
            tmp_for_wrfinput(ii,jj+1) = this%for_wrfinput(ii,jj)            
         end do
      end do

      ! Now reallocate memory in data type, and copy data back in
      deallocate(this%datetimes)
      deallocate(this%for_wrfbdy)
      deallocate(this%for_wrfinput)
      this%num_datetimes = this%num_datetimes + 1
      allocate(this%datetimes(this%num_datetimes))
      allocate(this%for_wrfbdy(this%num_datetimes))
      allocate(this%for_wrfinput(this%num_wrf_domains,this%num_datetimes))
      this%datetimes(:) = tmp_datetimes(:)
      this%for_wrfbdy(:) = tmp_for_wrfbdy(:)
      this%for_wrfinput(:,:) = tmp_for_wrfinput(:,:)
      deallocate(tmp_datetimes)
      deallocate(tmp_for_wrfbdy)
      deallocate(tmp_for_wrfinput)

   end subroutine insert_to_sorted_datetimes

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  compare_datetimes
   !
   ! DESCRIPTION: Function for determining order of two datetimes.
   !
   !---------------------------------------------------------------------------

   function compare_datetimes(datetime1,datetime2) result(diff)

      ! Arguments
      character(len=MAX_DATE_LEN),intent(in) :: datetime1
      character(len=MAX_DATE_LEN),intent(in) :: datetime2

      ! Return variable
      integer :: diff !  1 = datetime1 before datetime2
                      !  0 = datetime1 .eq. datetime2
                      ! -1 = datetime1 after datetime2
                            
      ! Local variables
      integer :: year1, year2
      integer :: month1, month2
      integer :: day1, day2
      integer :: hour1, hour2
      integer :: minute1, minute2
      integer :: second1, second2

      diff = 0
      if (datetime1 == datetime2) return

      diff = 1
      read(datetime1,"(I4.4,X,I2.2,X,I2.2,x,I2.2,x,I2.2,X,I2.2)") &
           year1,month1,day1,hour1,minute1,second1 
      read(datetime2,"(I4.4,X,I2.2,X,I2.2,x,I2.2,x,I2.2,X,I2.2)") &
           year2,month2,day2,hour2,minute2,second2 

      if (year1 > year2) then
         diff = -1
      else if (year1 < year2) then
         diff =  1
      else if (month1 > month2) then
         diff = -1
      else if (month1 < month2) then
         diff =  1
      else if (day1 > day2) then
         diff = -1
      else if (day1 < day2) then
         diff =  1
      else if (hour1 > hour2) then
         diff = -1
      else if (hour1 < hour2) then
         diff =  1
      else if (minute1 > minute2) then
         diff = -1
      else if (minute1 < minute2) then
         diff =  1
      else if (second1 > second2) then
         diff = -1
      else if (second1 < second2) then
         diff =  1
      end if

      return
   end function compare_datetimes

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  get_sorted_datetimes
   !
   ! DESCRIPTION: Public method for copying contents of data type into
   ! arrays.
   !
   !---------------------------------------------------------------------------

   subroutine get_sorted_datetimes(this,num_datetimes,num_wrf_domains, &
        datetimes, for_wrfbdy, for_wrfinput)

      ! Arguments
      type(sorted_datetimes),intent(in) :: this
      integer,intent(out) :: num_datetimes
      integer,intent(out) :: num_wrf_domains
      character(len=MAX_DATE_LEN),allocatable,intent(out) :: datetimes(:)
      logical,allocatable,intent(out) :: for_wrfbdy(:)
      logical,allocatable,intent(out) :: for_wrfinput(:,:)

      ! Sanity checks
      if (allocated(datetimes)) then
         print*,'ERROR, datetimes is already allocated!'
         stop
      end if
      if (allocated(for_wrfbdy)) then
         print*,'ERROR, for_wrfbdy is already allocated!'
         stop
      end if
      if (allocated(for_wrfinput)) then
         print*,'ERROR, for_wrfinput is already allocated!'
         stop
      end if
      if (this%num_wrf_domains == 0) then
         print*,'ERROR, number of wrf domains is zero.'
         print*,'Is data type initialized?'
         stop
      end if

      ! Copy the contents of the current data type
      num_datetimes = this%num_datetimes
      num_wrf_domains = this%num_wrf_domains
      allocate(datetimes(num_datetimes))
      allocate(for_wrfbdy(num_datetimes))
      allocate(for_wrfinput(num_wrf_domains,num_datetimes))
      datetimes(:) = this%datetimes(:)
      for_wrfbdy(:) = this%for_wrfbdy(:)
      for_wrfinput(:,:) = this%for_wrfinput(:,:)

   end subroutine get_sorted_datetimes
end module sorted_datetimes_mod
