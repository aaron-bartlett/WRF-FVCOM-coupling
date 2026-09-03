!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE: ncarg_util_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! Contains utility routines for interacting with NCAR Graphics.
!
!------------------------------------------------------------------------------

module ncarg_util_mod

   ! Defaults
   implicit none
   private

   ! Public data type for storing NCAR Graphics map projection data
   public :: map
   type map
      private
      integer :: grid_id
      integer :: map_proj
      real :: pole_latitude
      real :: pole_longitude
      real :: rotation
      real :: rlat_start
      real :: rlon_start
      real :: rlat_end
      real :: rlon_end
   end type map
   
   ! Public routines
   public :: create_map
   public :: start_ncarg, end_ncarg, plot_map, plot_contours,plot_title
   
contains

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  create_map
   !
   ! DESCRIPTION:  Creates data structure storing input map projection data
   ! for NCAR Graphics.  Reads ASCII file from prep_chem_sources.
   !
   !--------------------------------------------------------------------------

   function create_map(filename) result (this)

      ! Arguments
      character(len=*),intent(in) :: filename
      
      ! Return variable
      type(map) :: this

      ! Temporary variables
      integer :: grid
      integer :: map_proj
      real :: stdlat1,stdlat2,stdlon
      real :: centlat,centlon
      real :: rlats,rlons,rlate,rlone
      character(len=132) :: ckeys(9)
      character(len=5) :: ctmp_grid
      character(len=11) :: ctmp_projection
      character(len=19) :: ctmp_std
      character(len=16) :: ctmp_centlat
      character(len=17) :: ctmp_centlon
      character(len=18) :: ctmp_startlat
      character(len=19) :: ctmp_startlon
      character(len=16) :: ctmp_endlat
      character(len=17) :: ctmp_endlon
      integer :: iunit
      character(len=132) :: ctmp
      integer :: itmp
      real :: rtmp
      integer :: i

      ! Read from file
      iunit = 10
      open(unit=iunit,file=trim(filename),status='old')

      ckeys(1:2) = (/ 'grid:      ', &
                      'projection:' /)
      i = 1

      read(unit=iunit,fmt='(a5,1x,i3)') ctmp,itmp      
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      grid = itmp

      read(unit=iunit,fmt='(a11,1x,i3)') ctmp,itmp     
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      map_proj = itmp
         
      ckeys(:) = (/ 'standard_latitude1:', &
                    'standard_latitude2:', &
                    'standard_longitude:', &
                    'center_latitude:   ', &
                    'center_longitude:  ', &
                    'starting_latitude: ', &
                    'starting_longitude:', &
                    'ending_latitude:   ', &
                    'ending_longitude:  ' /)
      i = 1

      read(unit=iunit,fmt='(a19,1x,f15.3)') ctmp,rtmp   
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      stdlat1 = rtmp

      read(unit=iunit,fmt='(a19,1x,f15.3)') ctmp,rtmp   
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      stdlat2 = rtmp

      read(unit=iunit,fmt='(a19,1x,f15.3)') ctmp,rtmp   
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      stdlon = rtmp

      read(unit=iunit,fmt='(a16,1x,f15.3)') ctmp,rtmp      
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      centlat = rtmp

      read(unit=iunit,fmt='(a17,1x,f15.3)') ctmp,rtmp      
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      centlon = rtmp

      read(unit=iunit,fmt='(a18,1x,f15.3)') ctmp,rtmp      
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      rlats = rtmp

      read(unit=iunit,fmt='(a19,1x,f15.3)') ctmp,rtmp      
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      rlons = rtmp

      read(unit=iunit,fmt='(a16,1x,f15.3)') ctmp,rtmp      
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      rlate = rtmp

      read(unit=iunit,fmt='(a17,1x,f15.3)') ctmp,rtmp      
      if (trim(ctmp) .ne. trim(ckeys(i))) then
         print*,'ERROR reading map projection data!'
         stop
      end if
      i = i + 1
      rlone = rtmp

      close(unit=iunit)

      ! Populate the map structure
      this%grid_id = grid
      this%map_proj = map_proj
      this%rlat_start =  rlats
      this%rlon_start =  rlons
      this%rlat_end =    rlate
      this%rlon_end =    rlone
      this%pole_longitude = stdlon
      if (map_proj == 1) then ! Polar stereographic
         this%rotation = 0.
         this%pole_latitude = sign(90.,centlat)
      else if (map_proj == 3) then ! Lambert conformal
         this%rotation = stdlat1
         this%pole_latitude = stdlat2
      else if (map_proj == 9) then ! Mercator
         this%rotation = 0.
         this%pole_latitude = 0.
      else
         print*,'ERROR, invalid map projection!'
         print*,'map_proj = ',map_proj
         stop
      end if

      return
   end function create_map

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  start_ncarg
   !
   ! DESCRIPTION:  Initializes NCAR Graphics.
   !
   !--------------------------------------------------------------------------

   subroutine start_ncarg()
      
      ! Open GKS for output
      call opngks()

      ! Set some colors
      call gscr(1, 0, 1.00, 1.00, 1.00) ! White background
      call gscr(1, 1, 0.00, 0.00, 0.00) ! Black foreground

      return
   end subroutine start_ncarg

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  end_ncarg
   !
   ! DESCRIPTION:  Closes NCAR Graphics.
   !
   !--------------------------------------------------------------------------
   
   subroutine end_ncarg()

      ! Close GKS after all output
      call clsgks()

      return
   end subroutine end_ncarg

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  plot_map
   !
   ! DESCRIPTION:  Plots a simple map with NCAR Graphics.
   !
   !--------------------------------------------------------------------------

   subroutine plot_map(this,nx,ny)
      
      ! Arguments
      type(map), intent(in) :: this
      integer,intent(in) :: nx,ny

      ! Local variables
      integer :: jlts,jgrid,iusout,idot,ier
      real :: xa,xb,ya,yb,xxa,xxy,yya,yyb
      integer :: ltype
     
      ! Create map plot
      jlts =      2
      jgrid =    10
      iusout =    4
      idot =      0
      call supmap(this%map_proj,this%pole_latitude,this%pole_longitude, &
           this%rotation, &
           this%rlat_start,this%rlon_start,this%rlat_end,this%rlon_end, &
           jlts,jgrid,iusout,idot,ier)

      if (ier .ne. 0) then
         print*,'ERROR returned from supmap!'
         print*,'ier = ',ier
         stop
      end if
      
      ! Set line width
      call setusv('LW',1000)

      ! Relate the user coordinate system to the fractional coordinate system
      call getset(xa,xb,ya,yb,xxa,xxy,yya,yyb,ltype)
      call    set(xa,xb,ya,yb,1.,real(nx),1.,real(ny),ltype)

      return
   end subroutine plot_map

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  plot_contours
   !
   ! DESCRIPTION:  Plots a simple contour plot with NCAR Graphics.
   !
   !--------------------------------------------------------------------------

   subroutine plot_contours(nx,ny,data2d)

      ! Arguments
      integer,intent(in) :: nx,ny
      real,intent(in) :: data2d(nx,ny)

      call cpcnrc(data2d,nx,nx,ny,0.0,0.0,0,-1,0,0)

      return
   end subroutine plot_contours

   !--------------------------------------------------------------------------
   !
   ! ROUTINE:  plot_title
   !
   ! DESCRIPTION:  Plots a simple title with NCAR Graphics.
   !
   !--------------------------------------------------------------------------
   
   subroutine plot_title(title1,title2)
      
      ! Arguments
      character(len=*),intent(in) :: title1
      character(len=*),intent(in),optional :: title2

      ! Make sure the function code separator is set to something other
      ! than ":", which is found in the GrADS control file.
      call pcsetc('FC',"!")

      ! Put label at the top of the plot
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.05,trim(title1),0.01,0.,0.)
      if (present(title2)) then
         call plchhq(0.5,0.02,trim(title2),0.01,0.,0.)         
      end if

   end subroutine plot_title
end module ncarg_util_mod
