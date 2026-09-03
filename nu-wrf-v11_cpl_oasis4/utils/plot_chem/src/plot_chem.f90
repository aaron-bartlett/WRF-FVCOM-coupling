!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  plot_chem
!
! AUTHOR:
! Eric Kemp, NASA SSSO/Northrop Grumman
!
! DESCRIPTION:
! Simple program for plotting remapped emissions data from prep_chem_sources.
! Reads in GrADS binary file, GrADS control file, and a special ASCII "map"
! look-up table.  Plots are generated using NCAR Graphics.
!
! Users should create a symbolic link called 'grads.ctl' to point to the
! GrADS control file.  The program will use the data in this file to determine
! the name of the GrADS binary file and the "map" file.
!
!------------------------------------------------------------------------------

program plot_chem

   ! Modules
   use ncarg_util_mod
   use grads_util_mod

   ! Set default for declarations
   implicit none

   ! Local variables
   type(map) :: map_obj
   type(gradsctl) :: gradsctl_obj
   character(len=132) :: grads_binary_file
   character(len=132) :: grads_var_name
   character(len=132) :: mapfile
   real,allocatable :: data5d(:,:,:,:,:)
   character(len=80) :: level_title
   integer :: nx,ny,nz,nt,nvar
   integer :: i,j,k,m,n

   ! Get dimensions and names from GrADS control file.
   gradsctl_obj = create_gradsctl('grads.ctl')
   call get_dims_gradsctl(gradsctl_obj,nx,ny,nz,nt,nvar)
   grads_binary_file = get_binfilename_gradsctl(gradsctl_obj)

   ! Read gridded data from GrADS binary file.
   allocate(data5d(nx,ny,nz,nt,nvar))
   call read_grads_binary(trim(grads_binary_file),nx,ny,nz,nt,nvar,data5d)

   ! Create map object
   mapfile = get_mapfilename_gradsctl(gradsctl_obj)
   map_obj = create_map(trim(mapfile))

   ! Open NCAR Graphics
   call start_ncarg()

   ! Plot the data
   do n = 1,nvar
      do m = 1,nt
         do k = 1,nz

            ! Plot the map
            call plot_map(map_obj,nx,ny)

            ! Overlay contours
            call plot_contours(nx,ny,data5d(1,1,k,m,n))

            ! Plot title
            grads_var_name = get_varname_gradsctl(gradsctl_obj,n)
            write(level_title,'(A,I3,A,I3)') 'k = ',k,' of ',nz
            call plot_title(trim(grads_var_name),title2=trim(level_title))

            ! Finally, advance the frame
            call frame()
            
         end do
      end do
   end do

   ! Close NCAR Graphics
   call end_ncarg()
   
   ! The end
   deallocate(data5d)

   stop
end program plot_chem
