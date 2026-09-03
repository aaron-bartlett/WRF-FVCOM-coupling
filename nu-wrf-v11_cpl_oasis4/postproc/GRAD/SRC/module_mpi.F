#include <define.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_mpi
 use module_floating_point
 use netcdf
#if MPI > 0 
  use mpi
#endif
 implicit none

!--------------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
! Note that you MUST NOT use sdsu_fps for this procedure, since module_mpi is the first module to 
! be compiled. 
!
! Comments: 
!  If MPI is defined = 1 or = 2 in define.h, this module deals with MPI decomposition parameters. 
!  Currently MPI handles 1) file decomposition, and 2) domain decomposition.  
!  If you have a large number of file to simulate, option 1 ( mpi==1) gain best performance. 
!  If you have a large domain to simulate, option 2 (mpi==2) gain best performaince.    
!
! History:
! 01/2011  Toshi Matsui@NASA GSFC : Add defined-only gatherv routines. 
! 10/2009  Toshi Matsui@NASA GSFC : Add domain decomposition routines.
! 06/2008  Toshi Matsui@NASA GSFC : Add file decomposition routines. 
! 05/2008  Toshi Matsui@NASA GSFC : Initial
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 public   ! all variables and subourtines are accessible in module_simulator.

!
! These parameters are needed, even if not using MPI (this is for memory)
!
 logical :: masterproc   ! if true -> master processor, if not -> slave processors
 integer :: myrank       ! myrank=0 -> master processors, myrank>0 slave processors
 integer :: numproc_tot  ! # of total processors
 integer :: myn_start, myn_end            ! starting/ending file # for each processor
 integer :: myi_start, myi_end, myi_size  ! starting/ending j index for each processor (eastward direction)
 integer :: myj_start, myj_end, myj_size  ! starting/ending j index for each processor (northward direction)
 integer :: myk_start, myk_end, myk_size  ! starting/ending j index for each processor (vertical direction)


 


!
! memory bounds
!
 type memory_bounds
   integer :: is,ie
   integer :: js,je
   integer :: totnum !total tile number
 end type memory_bounds
 type (memory_bounds),allocatable,public :: memory(:) , buffer(:), domain(:), overlap(:,:)

!
! For OMP/MIC tile (this is for thread)
!
 type (memory_bounds),allocatable,public :: omp_tile(:) 

 
 integer,allocatable,public :: iproc_dom(:,:)  !MPI processor-rank ID for each domains. 

 integer, public :: ibs,ibe,jbs,jbe !extra buffer domain bounds


#if MPI > 0 

!
! MPI module interface
!

!
! Gather memory-size calculated value within domain-size allocated array into master processor
! Use : MPI_GATHERV with option 'TO_MASTER'  (for output)
! Use : MPI_ALLGATHERV with option "TO_ALL"  (for broadcast)
!
  interface mpi_sdsu_communicate
      module procedure mpi_sdsu_communicate_2d
      module procedure mpi_sdsu_communicate_3d
      module procedure mpi_sdsu_communicate_4d
      module procedure mpi_sdsu_communicate_2d_dble
      module procedure mpi_sdsu_communicate_3d_dble
      module procedure mpi_sdsu_communicate_4d_dble
  end interface

  interface mpi_sdsu_communicate_buf
      module procedure mpi_sdsu_communicate_buf_2d
      module procedure mpi_sdsu_communicate_buf_3d
      module procedure mpi_sdsu_communicate_buf_4d
  end interface


  interface mpi_sdsu_communicate_defined
      module procedure mpi_sdsu_communicate_defined_2d
      module procedure mpi_sdsu_communicate_defined_3d
      module procedure mpi_sdsu_communicate_defined_4d
  end interface

!
! Collect memory-size calculated value within memory-size allocated array into new domain-size array in master processor
! Use : MPI_GATHERVI
!
  interface mpi_sdsu_collect_tile  
      module procedure mpi_sdsu_collect_tile_2d
      module procedure mpi_sdsu_collect_tile_3d
  end interface

!
  interface mpi_sdsu_reduce
     module procedure mpi_sdsu_reduce_real_0d
     module procedure mpi_sdsu_reduce_real_1d
     module procedure mpi_sdsu_reduce_real_2d
     module procedure mpi_sdsu_reduce_real_3d
  end interface

  interface mpi_sdsu_allreduce
     module procedure mpi_sdsu_allreduce_real_0d
  end interface

  integer,private :: ierr,rc   !index for MPI error statistics

  integer,private ::   i_start,   i_end,   i_size  !starting/ending j index for domain (eastward direction)
  integer,private ::   j_start,   j_end,   j_size  !starting/ending j index for domain (northward direction)
  integer,private ::   k_start,   k_end,   k_size  !starting/ending j index for domain (vertical direction)

  integer,allocatable,private :: myistr(:), myilen(:)  ! starting index and length in lon direction
                                               ! for each processors
  integer,allocatable,private :: myjstr(:), myjlen(:)  ! starting index and length in lat direction
                                               ! for each processors
  integer,allocatable,private :: mykstr(:), myklen(:)  ! starting index and length in vertical direction
                                               ! for each processors

#endif 


 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

#if MPI > 0 

 subroutine mpi_sdsu_init
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Initialize basic properties of MPI parameter. 
!
! History:
! 06/2008  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------

!
! initialize MPI
!
   call MPI_INIT(ierr)
   if (ierr /= MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_init: Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if

!
! get rank for each processor
!
   call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
   if (ierr /= MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_init: Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if

   if ( myrank == 0 ) then
      masterproc = .true.
   else
      masterproc = .false.
   endif

!
! get total number of processors (numproc_tot)
!
   call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc_tot, ierr)
   if(masterproc) print*,'MSG mpi_sdsu_init: total number of processor is ',numproc_tot
   if (ierr /= MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_init: Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if


   return

 end subroutine mpi_sdsu_init

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_file(nmax,exit_proc)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Compute file decomposition. 
!
! History:
! 06/2008  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: nmax  !maximum input file number
 logical,intent(out) :: exit_proc  !exit processing 
 integer :: dfile  !file # increment


!
! define file increment for each processor (this must be >= 1)
!
   dfile = max(1 ,  nint( real(nmax)/real(numproc_tot) ) )

!
! Derive Starting/Ending file index for each processor rank
!
   myn_start =  max( 1, dfile*(myrank  )+1 )
   myn_end   =  min( nmax, dfile*(myrank+1)   )
   exit_proc = .false.
   if( myn_start > nmax ) then
     myn_end = myn_start
     print*,'MSG mpi_sdsu_file: Wasting processors at myrank=',myrank, 'exit the process'
     print*,'myn_end=',myn_end, myn_start, 'dfile=',dfile, nmax, numproc_tot
     exit_proc = .true.
   endif

   write(*,'(A28,I4,A16,I4,A15,I4)') &
    'MSG mpi_sdsu(file): myrank =', myrank,'  Start file # =', myn_start, '   End file # =', myn_end

   return
 end subroutine mpi_sdsu_file

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_tile(imin,imax,jmin,jmax,kmin,kmax)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Compute domain decomposition in i-j-direction.  
!  This algorithm intends to maximize i-size of tile shape for better CPU cash use. 
!
! History:
! 10/2009  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: imin, imax
 integer,intent(in) :: jmin, jmax 
 integer,intent(in) :: kmin, kmax

 integer :: t,j        ! looping
 integer :: di ,dj     ! j-index # increment
 integer :: tot_grid, tot_grid_i, tot_grid_j
 integer :: resid      ! residual processors
 integer :: tile_size  ! tile size (total grid #)
 logical :: cannot_find_size
 integer :: irank
 integer :: num_i, num_j
 integer :: it,jt

!
! Total domain parameters
!
 i_start=imin ; i_end=imax ; j_start=jmin ; j_end=jmax ; k_start=kmin ; k_end=kmax
 i_size=i_end-i_start+1    ; j_size=j_end-j_start+1    ; k_size=k_end-k_start+1

!
! allocate
!
 if(.not. allocated(myistr)) allocate(myistr(0:numproc_tot-1))
 if(.not. allocated(myilen)) allocate(myilen(0:numproc_tot-1))
 if(.not. allocated(myjstr)) allocate(myjstr(0:numproc_tot-1))
 if(.not. allocated(myjlen)) allocate(myjlen(0:numproc_tot-1))
 if(.not. allocated(mykstr)) allocate(mykstr(0:numproc_tot-1))
 if(.not. allocated(myklen)) allocate(myklen(0:numproc_tot-1))
 if(.not. allocated(memory)) allocate(memory(0:numproc_tot-1))
 if(.not. allocated(memory)) allocate(iproc_dom(imin:imax,jmin:jmax))
 if(.not. allocated(domain)) allocate(domain(0:numproc_tot-1))

 tot_grid_i = (imax-imin+1)
 tot_grid_j = (jmax-jmin+1) 
 tot_grid = tot_grid_i * tot_grid_j  !total grid number

 resid = MOD( tot_grid,  numproc_tot )



 if(resid == 0 ) then  ! no residual case ---------------------------------------------------

 tile_size = tot_grid / numproc_tot  !tile grid (total)

 cannot_find_size = .true.

!
! this loop changes tile size from (j=1,i=tile_size) -> (j=tile_size,i=1)
! to find the best tile shape
!
   tile_loop: do t = 1, tile_size

    if( MOD(tile_size, t) == 0 ) then
       dj =  t
       di =  tile_size / t

       if( MOD( tot_grid_i , di ) == 0 .and. MOD( tot_grid_j , dj ) == 0 ) then
             num_i = tot_grid_i / di
             num_j = tot_grid_j / dj

             irank = 0 
             do jt = 1, num_j  
                do it = 1, num_i
                   myistr(irank) = di*(it-1) +1
                   myilen(irank) = di
                   myjstr(irank) = dj*(jt-1) +1
                   myjlen(irank) = dj
                   irank = irank + 1 
                enddo
             enddo

           !
           ! Domain decomposition (asign each start i,j index here)
           !
           myi_start =  myistr(myrank)
           myi_end   =  myistr(myrank) + myilen(myrank) - 1
           myi_size  = myi_end - myi_start + 1
           myj_start =  myjstr(myrank) 
           myj_end   =  myjstr(myrank) + myjlen(myrank) - 1
           myj_size  = myj_end - myj_start + 1
           cannot_find_size = .false.
           exit tile_loop 
       endif

    endif

   enddo tile_loop

    if(cannot_find_size) then
      call mpi_wait_for_master
      if(masterproc) then
        print*,'MSG mpi_sdsu_tile: cannot find appropriate tile size'
        print*,'Abort MPI --> Change number of processors '
        call example_sdsu_tile(imin,imax,jmin,jmax,kmin,kmax)
      endif
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr) 
    endif

 else    !residual processor ------------------------------------------------------------------

   call mpi_wait_for_master
   if(masterproc) then
      print*,'MSG mpi_sdsu_tile: cannot find appropriate tile size'
      print*,' --> Change the number of CPUs and nodes'
      call example_sdsu_tile(imin,imax,jmin,jmax,kmin,kmax)
   endif

   call mpi_sdsu_abort('MSG mpi_sdsu_tile: cannot find appropriate tile size', myrank)
 
 endif

!
! k index is identical between the processors
!
 myk_start =  kmin 
 myk_end   =  kmax 
 myk_size  =  kmax - kmin + 1

 mykstr(0:numproc_tot-1)=kmin
 myklen(0:numproc_tot-1)=kmax - kmin + 1

!
! memory domain
!
 do irank = 0, numproc_tot-1

    memory(irank)%is  =  myistr(irank)
    memory(irank)%ie  =  myistr(irank) + myilen(irank) - 1
    memory(irank)%js  =  myjstr(irank)
    memory(irank)%je  =  myjstr(irank) + myjlen(irank) - 1
    memory(irank)%totnum = myilen(irank) *  myjlen(irank)

    domain(irank)%is  = imin
    domain(irank)%ie  = imax
    domain(irank)%js  = jmin
    domain(irank)%je  = jmax
    domain(irank)%totnum = (imax-imin+1)*(jmax-jmin+1)

 enddo

 return
 end subroutine mpi_sdsu_tile

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine example_sdsu_tile(imin,imax,jmin,jmax,kmin,kmax)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Compute domain decomposition in i-j-direction.  
!  This algorithm intends to maximize i-size of tile shape for better CPU cash use. 
!
! History:
! 10/2009  Toshi Matsui@NASA GSFC : Initial 
!
! References:
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: imin, imax
 integer,intent(in) :: jmin, jmax
 integer,intent(in) :: kmin, kmax
 integer :: numproc_tot  ! # of total processors

 integer :: n,nt,t,j        ! looping
 integer :: di ,dj     ! j-index # increment
 integer :: tot_grid, tot_grid_i, tot_grid_j
 integer :: resid      ! residual processors
 integer :: tile_size  ! tile size (total grid #)
 logical :: cannot_find_size
 integer :: irank
 integer :: num_i, num_j
 integer :: it,jt
 integer :: node_number(1:12)
 integer :: total_cpu

 tot_grid_i = (imax-imin+1)
 tot_grid_j = (jmax-jmin+1)
 tot_grid = tot_grid_i * tot_grid_j  !total grid number

 total_cpu = 1024

 write(*,*), ''
 write(*,*), ' See Example below for DISCOVER/PLEIADES cluster'
 write(*,*), ''
 write(*,FMT='(A107)'), '       ncpus=1 ncpus=2 ncpus=3 ncpus=4 ncpus=5 ncpus=6 ncpus=7 ncpus=8 ncpus=9 ncpus=10 ncpus=11 ncpus=12'
 do n = 1, total_cpu

    !
    ! TEST one
    !
    numproc_tot = n
    resid = MOD( tot_grid,  numproc_tot )

    if(resid == 0) then

       tile_size = tot_grid / numproc_tot  !tile grid (total)
       cannot_find_size = .true.

       !
       ! this loop changes tile size from (j=1,i=tile_size) -> (j=tile_size,i=1)
       ! to find the best tile shape
       !
       tile_loop: do t = 1, tile_size

       if( MOD(tile_size, t) == 0 ) then
         dj =  t
         di =  tile_size / t
         if( MOD( tot_grid_i , di ) == 0 .and. MOD( tot_grid_j , dj ) == 0 ) then
             num_i = tot_grid_i / di
             num_j = tot_grid_j / dj
             irank = 0
!             print*,'find out total number=',t,'di',di,'dj',dj,num_i,num_j
             cannot_find_size = .false.
             exit tile_loop
         endif
       endif

       enddo tile_loop

      !
      ! write right node number for a given processor number. 
      !
      if( .not. cannot_find_size ) then
        do nt = 1, 12
           resid = MOD( numproc_tot, nt )
           if(resid /= 0)  then
              node_number(nt) = 0
           else
              node_number(nt) = numproc_tot/nt
           endif
        enddo
        write(*,FMT='(A8,12I8)'),'node # =', node_number(:)
      endif

    endif

 enddo

 return
 end subroutine example_sdsu_tile

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_tile_omp(num_tiles, imin,imax,jmin,jmax,kmin,kmax)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Compute domain decomposition in i-j-direction.  
!  This algorithm intends to maximize i-size of tile shape for better CPU cash use. 
!
! History:
! 10/2009  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: num_tiles   !number of tiles for OMP/MIC
 integer,intent(in) :: imin, imax  !memory domain
 integer,intent(in) :: jmin, jmax  !memory domain
 integer,intent(in) :: kmin, kmax  !memory domain


 integer :: t,j        ! looping
 integer :: di ,dj     ! j-index # increment
 integer :: tot_grid, tot_grid_i, tot_grid_j
 integer :: resid      ! residual processors
 integer :: tile_size  ! tile size (total grid #)
 logical :: cannot_find_size
 integer :: irank
 integer :: num_i, num_j
 integer :: it,jt
 integer :: itile

!
! comments
!
 if(masterproc) then
   print*,'MSG mpi_sdsu_tile_omp: OMP tile decomposition'
   print*,'total memory domain size is i x j =',imax-imin+1,'x',jmax-jmin+1
 endif

!
! Total domain parameters
!
 i_start=imin ; i_end=imax ; j_start=jmin ; j_end=jmax ; k_start=kmin ; k_end=kmax
 i_size=i_end-i_start+1    ; j_size=j_end-j_start+1    ; k_size=k_end-k_start+1

!
! allocate
!
 if(.not. allocated(omp_tile)) allocate(omp_tile(num_tiles))


 tot_grid_i = (imax-imin+1)
 tot_grid_j = (jmax-jmin+1) 
 tot_grid = tot_grid_i * tot_grid_j  !total grid number of memory domain

 resid = MOD( tot_grid,  num_tiles )

 print*,'MSG mpi_sdsu_tile_omp: omp decomp', tot_grid_i, tot_grid_j, tot_grid, resid, numproc_tot


 if(resid == 0 ) then  ! no residual case ---------------------------------------------------

 tile_size = tot_grid / num_tiles  ! tile grid (total)

 cannot_find_size = .true.

!
! this loop changes tile size from (j=1,i=tile_size) -> (j=tile_size,i=1)
! to find the best tile shape
!
   tile_loop: do t = 1, tile_size

    if( MOD(tile_size, t) == 0 ) then
       dj =  t
       di =  tile_size / t

       if( MOD( tot_grid_i , di ) == 0 .and. MOD( tot_grid_j , dj ) == 0 ) then
             num_i = tot_grid_i / di
             num_j = tot_grid_j / dj
             if(masterproc) print*,'MSG mpi_sdsu_tile_omp: OMP tile size =',di,'x',dj
             if(masterproc) print*,''

             itile = 1   !tile index start from 1
             do jt = 1, num_j  
                do it = 1, num_i
                   omp_tile(itile)%is = di*(it-1) + imin
                   omp_tile(itile)%ie = di*it     + imin - 1
                   omp_tile(itile)%js = dj*(jt-1) + jmin
                   omp_tile(itile)%je = dj*jt     + jmin - 1
                   omp_tile(itile)%totnum = di*dj
                   if( myrank == 10 ) print*,itile, omp_tile(itile)
                   itile = itile + 1 
                enddo
             enddo

           cannot_find_size = .false.
           exit tile_loop 
       endif

    endif

   enddo tile_loop

    if(cannot_find_size) then
      call mpi_wait_for_master
      if(masterproc) then
        print*,'MSG mpi_sdsu_tile_omp: cannot find appropriate omp tile size'
        print*,'Abort MPI --> Change num_tiles_omp in Configuration file '
      endif
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr) 
    endif

 else    !residual processor ------------------------------------------------------------------

      call mpi_wait_for_master
      if(masterproc) then
        print*,'MSG mpi_sdsu_tile_omp: cannot find appropriate omp tile size'
        print*,'Abort MPI --> Change num_tiles_omp in Configuration file '
      endif
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)

 endif


 end subroutine mpi_sdsu_tile_omp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_allreduce_real_0d(method, real_0d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program reduce (sum) 0dimensional real array from all slave processors to mater processor.
!   Typically used for statistical purpose.
!
! History:
! 07/2010 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: method     !MPI colletable method
 real,intent(inout)  :: real_0d   ! input real parameter
 real :: tmp

!
! send parameter, and reduce by various method 
!
 select case(trim(method))
 case('SUM')

   call MPI_ALLREDUCE(real_0d, tmp, 1,  mpi_real,   &
                  mpi_sum, mpi_comm_world, ierr  )

 case('MIN')

   call MPI_ALLREDUCE(real_0d, tmp, 1,  mpi_real,   &
                  mpi_min, mpi_comm_world, ierr  )

 case('MAX')

   call MPI_ALLREDUCE(real_0d, tmp, 1,  mpi_real,   &
                  mpi_max, mpi_comm_world, ierr  )

 case default

     print *,'MSG mpi_sdsu_allreduce_real_0d: Error MPI_REDUCE. no such method',myrank
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)

 end select

 if (ierr .ne. MPI_SUCCESS) then
     print *,'MSG mpi_sdsu_reduce_real_0d: Error MPI_REDUCE. Terminating. myrank=',myrank
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
 endif

 real_0d = tmp

 return
 end subroutine mpi_sdsu_allreduce_real_0d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_reduce_real_0d(method, real_0d ) 
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program reduce (sum) 0dimensional real array from all slave processors to mater processor.
!   Typically used for statistical purpose.
!
! History:
! 07/2010 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,intent(inout)  :: real_0d   ! input real parameter
 real :: tmp

!
! send parameter, and reduce by summing all of them
!
 call MPI_REDUCE(real_0d, tmp, 1,  mpi_real,   &
                 mpi_sum, 0, mpi_comm_world, ierr  )
 if (ierr .ne. MPI_SUCCESS) then
     print *,'MSG mpi_sdsu_reduce_real_0d: Error MPI_REDUCE. Terminating. myrank=',myrank
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
 endif
 real_0d = tmp

 return 
 end subroutine mpi_sdsu_reduce_real_0d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_reduce_real_1d(method, real_1d ) 
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program reduce (sum) 1dimensional real array from all slave processors to mater processor.
!   Typically used for statistical purpose.
!
! History:
! 07/2010 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,dimension(:),intent(inout)  :: real_1d   ! input real 1d

 integer :: size_vector   ! elevental size of vector
 integer :: n, i1         ! looping indice
 integer :: bnd(1)        ! bound of array
 real,allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector

 !
 ! dimension bound
 ! 
 bnd = UBOUND(real_1d)

 size_vector = bnd(1)   !size of vector

!
! allocate tempoeral sending and receiving vector 
! 
 allocate( vecsnd(1:size_vector))
 allocate( vecrcv(1:size_vector))

!
! initializing sending vector
! 
 n = 1 
 do i1 = 1, bnd(1)
    vecsnd(n) = real_1d(i1)  !sending vector
    n=n+1
 enddo
 
!
! send parameter, and reduce by summing all of them
!
 call MPI_REDUCE(vecsnd, vecrcv, size_vector,  mpi_real,   &
                 mpi_sum, 0, mpi_comm_world, ierr  )
 if (ierr .ne. MPI_SUCCESS) then
     print *,'MSG mpi_sdsu_reduce_real_1d: Error MPI_REDUCE. Terminating. myrank=',myrank
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
 endif

!
! outputing receiving vector into real_1d 
!
 n = 1
 do i1 = 1, bnd(1)
    real_1d(i1) = vecrcv(n) !output sum vector
    n=n+1
 enddo 

 deallocate(vecsnd,vecrcv)

 return 
 end subroutine mpi_sdsu_reduce_real_1d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_reduce_real_2d(method, real_2d ) 
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program reduce (sum) 2dimensional real array from all slave processors to mater processor.
!   Typically used for statistical purpose.
!
! History:
! 07/2010 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,dimension(:,:),intent(inout)  :: real_2d   ! input real 2d

 integer :: size_vector   ! elevental size of vector
 integer :: n, i1,i2         ! looping indice
 integer :: bnd(2)        ! bound of array
 real,allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector

 !
 ! dimension bound
 ! 
 bnd = UBOUND(real_2d)

 size_vector = bnd(1)*bnd(2)   !size of vector

!
! allocate tempoeral sending and receiving vector 
! 
 allocate( vecsnd(1:size_vector))
 allocate( vecrcv(1:size_vector))

!
! initializing sending vector
! 
 n = 1 
 do i2 = 1, bnd(2) ; do i1 = 1, bnd(1)
    vecsnd(n) = real_2d(i1,i2)  ! sending vector
    n=n+1
 enddo ; enddo
 
!
! send parameter, and reduce by summing all of them
!
 call MPI_REDUCE(vecsnd, vecrcv, size_vector,  mpi_real,   &
                 mpi_sum, 0, mpi_comm_world, ierr  )
 if (ierr .ne. MPI_SUCCESS) then
     print *,'MSG mpi_sdsu_reduce_real_2d: Error MPI_REDUCE. Terminating. myrank=',myrank
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
 endif

!
! outputing receiving vector into real_2d 
!
 n = 1
 do i2 = 1, bnd(2) ; do i1 = 1, bnd(1)
    real_2d(i1,i2) = vecrcv(n) !output sum vector
    n=n+1
 enddo ; enddo

 deallocate(vecsnd,vecrcv)

 return 
 end subroutine mpi_sdsu_reduce_real_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_reduce_real_3d( method, real_3d ) 
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program reduce (sum) 3dimensional real array from all slave processors to mater processor.
!   Typically used for statistical purpose.
!
! History:
! 07/2010 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,dimension(:,:,:),intent(inout)  :: real_3d   ! input real 3d

 integer :: size_vector   ! elevental size of vector
 integer :: n, i1,i2,i3         ! looping indice
 integer :: bnd(3)        ! bound of array
 real,allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector

 !
 ! dimension bound
 ! 
 bnd = UBOUND(real_3d)

 size_vector = bnd(1)*bnd(2)*bnd(3)   !size of vector

!
! allocate tempoeral sending and receiving vector 
! 
 allocate( vecsnd(1:size_vector))
 allocate( vecrcv(1:size_vector))

!
! initializing sending vector
! 
 n = 1 
 do i3 = 1, bnd(3) ; do i2 = 1, bnd(2) ; do i1 = 1, bnd(1)
    vecsnd(n) = real_3d(i1,i2,i3)  ! sending vector
    n=n+1
 enddo ; enddo ; enddo
 
!
! send parameter, and reduce by summing all of them
!
 call MPI_REDUCE(vecsnd, vecrcv, size_vector,  mpi_real,   &
                 mpi_sum, 0, mpi_comm_world, ierr  )
 if (ierr .ne. MPI_SUCCESS) then
     print *,'MSG mpi_sdsu_reduce_real_3d: Error MPI_REDUCE. Terminating. myrank=',myrank
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
 endif

!
! outputing receiving vector into real_3d 
!
 n = 1
 do i3 = 1, bnd(3) ; do i2 = 1, bnd(2) ; do i1 = 1, bnd(1)
    real_3d(i1,i2,i3) = vecrcv(n) !output sum vector
    n=n+1
 enddo ; enddo ; enddo

 deallocate(vecsnd,vecrcv)

 return 
 end subroutine mpi_sdsu_reduce_real_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_collect_tile_2d( var2d , out2d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program distribute tile 2d array from slave processors to mater processor.
!   Typically used for output purporse. 
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 real,dimension(myi_start:myi_end, myj_start:myj_end), intent(in)  :: var2d   !tile 2d input 
 real,dimension(i_start:i_end,j_start:j_end), intent(out) :: out2d   !domain 2d output

 real,allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2       !1~2nd memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 

!
! allocate disp totallength point
!
 if( .not. allocated(totlen)) allocate(totlen(0:numproc_tot-1))
 if( .not. allocated(disp  )) allocate(disp  (0:numproc_tot-1))

!
! allocate sending/receiving vectors 
!
 size_vecsnd = SIZE( var2d )  !sending vector
 if( .not. allocated(vecsnd))  allocate(  vecsnd(0: size_vecsnd-1 )  )

 size_vecrcv = SIZE( out2d )  !receiving vector (here)
 if( .not. allocated(vecrcv))  allocate(  vecrcv(0: size_vecrcv-1 )  )

! print*,'size_vecsnd=', size_vecsnd, 'size_vecrcv=', size_vecrcv

!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 if( masterproc ) then
     do irank = 0, numproc_tot-1
        totlen(irank) = myilen(irank) * myjlen(irank)  !total length
     enddo
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif

!
! Vectorize var2d
!
 n = 0
 do i2 = myj_start, myj_end ; do i1 = myi_start, myi_end
    vecsnd(n) = var2d(i1,i2)  !sending vector
    n=n+1
 enddo ; enddo


!
! send vector to master node
!
 
!  print*,'mpi module vecsnd',vecsnd
!  print*,'myrank=',myrank,'mpi module',size_vecsnd, mpi_real

  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                   vecrcv, totlen, disp, mpi_real, 0, mpi_comm_world, ierr)

!  print*,'Finished gatherv  myrank=',myrank,'mpi module',totlen,disp
  

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_collect_tile_2d: Error MPI_GATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

!
! Distribute vector to each node
!
!   call MPI_ALLGATHERV

!
! Vector -> output 2d Array
!

 if(masterproc) then
    do irank = 0, numproc_tot-1
       is = myistr(irank)
       ie = myistr(irank) + myilen(irank) - 1
       js = myjstr(irank)
       je = myjstr(irank) + myjlen(irank) - 1
       n=disp(irank)
       do i2 = js, je ; do i1 = is, ie
          out2d(i1,i2) = vecrcv(n)  !sending vector (here)
          n=n+1
       enddo ; enddo

    enddo
 endif
 

 return
 end subroutine mpi_sdsu_collect_tile_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_collect_tile_3d( kmax, var3d , out3d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program distribute tile 3d array from slave processors to mater processor.
!   Typically used for output purporse. 
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: kmax   !maximum vertical profile
 real,dimension(myi_start:myi_end, myj_start:myj_end,1:kmax),intent(in)::var3d !tile 3d input
 real,dimension(i_start:i_end,j_start:j_end,1:kmax), intent(out) :: out3d   !domain 3d output

 real,allocatable,dimension(:),save :: vecsnd, vecrcv ! temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp    ! vector length and disp points
 integer :: bnd(3)      ! upper memory bound for 3D array
 integer :: size_vecsnd, size_vecrcv  ! size of sending/receiving vector
 integer :: i1,i2,i3    ! 1~3th memory loop
 integer :: n           ! vector count
 integer :: irank       ! myrank loop
 integer :: is,ie,js,je,ks,ke ! starting/ending i j indices 

 if( .not. allocated(totlen)) allocate(totlen(0:numproc_tot-1))
 if( .not. allocated(disp  )) allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(out3d)  !(here)

!
! allocate sending/receiving vector
!
 size_vecsnd = SIZE( var3d )  !sending vector
 if( .not. allocated(vecsnd)) allocate(  vecsnd(0: size_vecsnd-1 )  ) 

 size_vecrcv = SIZE( out3d )  !receiving vector (here)
 if( .not. allocated(vecrcv)) allocate(  vecrcv(0: size_vecrcv-1 )  )


!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!

 if( masterproc ) then  
     do irank = 0, numproc_tot-1  
        totlen(irank) = myilen(irank) * myjlen(irank) * bnd(3)  !total length
!         totlen(irank) = myilen(irank) * myjlen(irank) * myklen(irank)
     enddo
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif

!
! Vectorize var3d
!
 n = 0
 do i3 = 1, kmax ; do i2 = myj_start, myj_end ; do i1 = myi_start, myi_end
    vecsnd(n) = var3d(i1,i2,i3)    ! sending vector
    n=n+1
 enddo ; enddo ; enddo


!
! send vector to master node
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                   vecrcv, totlen, disp, mpi_real, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_collect_tile_3d: Error MPI_GATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif


!
! Distribute vector to each node
!
!   call MPI_ALLGATHERV

!
! Vector -> 3d Array
!

 if(masterproc) then

    do irank = 0, numproc_tot-1

       is = myistr(irank) 
       ie = myistr(irank) + myilen(irank) - 1
       js = myjstr(irank) 
       je = myjstr(irank) + myjlen(irank) - 1
       ks = 1 
       ke = kmax
       n=disp(irank)
       do i3 = ks, ke ; do i2 = js, je ; do i1 = is, ie 
          out3d(i1,i2,i3) = vecrcv(n)  ! sending vector (here)
          n=n+1
       enddo ; enddo ; enddo

    enddo

 endif


 return
 end subroutine mpi_sdsu_collect_tile_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_buf_2d( ibs,ibe,jbs,jbe, var2d  )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 2d array to mater (slave) processor via MPI library.
!   Note that about var2d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)         :: ibs,ibe,jbs,jbe  !extended buffer bounds
 real,dimension(ibs:ibe,jbs:jbe), intent(inout)  :: var2d   !assumed-shape 2d array

 real,allocatable,dimension(:) :: vecsnd, vecrcv   !temporal sending & receiving vector
 integer :: bnd(2)  !upper memory bound for 3D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2   !1~2nd memory loop
 integer :: n           !vector count
 integer :: i,j
 integer :: irank_src,irank_dest
 integer status(MPI_STATUS_SIZE) 
 integer ,parameter :: TAG=1
!
! allocate sending and receiving vector using maximum possible size
!
 allocate( vecsnd( 1: maxval( memory(:)%totnum ) ) ) 
 allocate( vecrcv( 1: maxval( memory(:)%totnum ) ) )                


 do irank_src = 0, numproc_tot-1
 do irank_dest = 0, numproc_tot-1


    if( overlap(irank_src,irank_dest)%totnum > 0 .and. & ! enough grid overlapped
        irank_src /= irank_dest ) then ! src and dest proc are different


        if( myrank == irank_src ) then  !source irank is equal rank

            size_vecsnd = overlap(irank_src,irank_dest)%totnum

            n = 0
            do j = overlap(irank_src,irank_dest)%js, overlap(irank_src,irank_dest)%je
            do i = overlap(irank_src,irank_dest)%is, overlap(irank_src,irank_dest)%ie
               n = n + 1
               vecsnd(n) = var2d(i,j)
            enddo ; enddo

            call MPI_SEND( vecsnd(1:size_vecsnd), size_vecsnd, MPI_REAL, irank_dest, &
                            TAG, MPI_COMM_WORLD, ierr ) 

            if (ierr .ne. MPI_SUCCESS) then
               print *,'MSG mpi_sdsu_communicate_buf_2d: Error MPI_SEND. Terminating. myrank=',myrank
               call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
            endif

        endif

       if( myrank == irank_dest ) then

            size_vecrcv = overlap(irank_src,irank_dest)%totnum 

            call MPI_RECV( vecrcv(1:size_vecrcv), size_vecrcv, MPI_REAL, irank_src, &
                          TAG, MPI_COMM_WORLD, status, ierr )   

            if (ierr .ne. MPI_SUCCESS) then
               print *,'MSG mpi_sdsu_communicate_buf_2d: Error MPI_RECV. Terminating. myrank=',myrank
               call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
            endif

            n = 0
            do j = overlap(irank_src,irank_dest)%js, overlap(irank_src,irank_dest)%je
            do i = overlap(irank_src,irank_dest)%is, overlap(irank_src,irank_dest)%ie
               n = n + 1
               var2d(i,j) = vecrcv(n)
            enddo ; enddo

       endif

    endif


 enddo
 enddo


 deallocate( vecsnd, vecrcv )                


 return
 end subroutine mpi_sdsu_communicate_buf_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_buf_3d( ibs,ibe,jbs,jbe, var3d  )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 2d array to mater (slave) processor via MPI library.
!   Note that about var2d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)         :: ibs,ibe,jbs,jbe  !extended buffer bounds
 real,dimension(ibs:,jbs:,:), intent(inout)  :: var3d   !assumed-shape 3d array

 real,allocatable,dimension(:) :: vecsnd, vecrcv   !temporal sending & receiving vector
 integer :: bnd(3)  !upper memory bound for 3D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2,i3   !1~2nd memory loop
 integer :: n           !vector count
 integer :: i,j
 integer :: irank_src,irank_dest
 integer status(MPI_STATUS_SIZE) 
 integer ,parameter :: TAG=1

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(var3d)

!
! allocate sending and receiving vector using maximum possible size
!
 allocate( vecsnd( 1: bnd(3)*maxval( memory(:)%totnum ) ) ) 
 allocate( vecrcv( 1: bnd(3)*maxval( memory(:)%totnum ) ) )                


 do irank_src = 0, numproc_tot-1
 do irank_dest = 0, numproc_tot-1


    if( overlap(irank_src,irank_dest)%totnum > 0 .and. & ! enough grid overlapped
        irank_src /= irank_dest ) then ! src and dest proc are different


        if( myrank == irank_src ) then  !source irank is equal rank

            size_vecsnd = bnd(3)*overlap(irank_src,irank_dest)%totnum

            n = 0
            do i3= 1,bnd(3)
            do j = overlap(irank_src,irank_dest)%js, overlap(irank_src,irank_dest)%je
            do i = overlap(irank_src,irank_dest)%is, overlap(irank_src,irank_dest)%ie
               n = n + 1
               vecsnd(n) = var3d(i,j,i3)
            enddo ; enddo ; enddo

            call MPI_SEND( vecsnd(1:size_vecsnd), size_vecsnd, MPI_REAL, irank_dest, &
                            TAG, MPI_COMM_WORLD, ierr ) 

            if (ierr .ne. MPI_SUCCESS) then
               print *,'MSG mpi_sdsu_communicate_buf_3d: Error MPI_SEND. Terminating. myrank=',myrank
               call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
            endif

        endif

       if( myrank == irank_dest ) then

            size_vecrcv = bnd(3)*overlap(irank_src,irank_dest)%totnum 

            call MPI_RECV( vecrcv(1:size_vecrcv), size_vecrcv, MPI_REAL, irank_src, &
                          TAG, MPI_COMM_WORLD, status, ierr )   
            if (ierr .ne. MPI_SUCCESS) then
               print *,'MSG mpi_sdsu_communicate_buf_3d: Error MPI_RECV. Terminating. myrank=',myrank
               call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
            endif

            n = 0
            do i3= 1,bnd(3)
            do j = overlap(irank_src,irank_dest)%js, overlap(irank_src,irank_dest)%je
            do i = overlap(irank_src,irank_dest)%is, overlap(irank_src,irank_dest)%ie
               n = n + 1
               var3d(i,j,i3) = vecrcv(n)
            enddo ; enddo ; enddo

       endif

    endif


 enddo
 enddo


 deallocate( vecsnd, vecrcv )                


 return
 end subroutine mpi_sdsu_communicate_buf_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_buf_4d( ibs,ibe,jbs,jbe, var4d  )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 2d array to mater (slave) processor via MPI library.
!   Note that about var2d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in)         :: ibs,ibe,jbs,jbe  !extended buffer bounds
 real,dimension(ibs:,jbs:,:,:), intent(inout)  :: var4d   !assumed-shape 4d array

 real,allocatable,dimension(:) :: vecsnd, vecrcv   !temporal sending & receiving vector
 integer :: bnd(4)  !upper memory bound for 4D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2,i3,i4   !1~2nd memory loop
 integer :: n           !vector count
 integer :: i,j
 integer :: irank_src,irank_dest
 integer status(MPI_STATUS_SIZE) 
 integer ,parameter :: TAG=1

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(var4d)

!
! allocate sending and receiving vector using maximum possible size
!
 allocate( vecsnd( 1: bnd(4)*bnd(3)*maxval( memory(:)%totnum ) ) ) 
 allocate( vecrcv( 1: bnd(4)*bnd(3)*maxval( memory(:)%totnum ) ) )                


 do irank_src = 0, numproc_tot-1
 do irank_dest = 0, numproc_tot-1


    if( overlap(irank_src,irank_dest)%totnum > 0 .and. & ! enough grid overlapped
        irank_src /= irank_dest ) then ! src and dest proc are different


        if( myrank == irank_src ) then  !source irank is equal rank

            size_vecsnd = bnd(4)*bnd(3)*overlap(irank_src,irank_dest)%totnum 

            n = 0
            do i4= 1,bnd(4)
            do i3= 1,bnd(3)
            do j = overlap(irank_src,irank_dest)%js, overlap(irank_src,irank_dest)%je
            do i = overlap(irank_src,irank_dest)%is, overlap(irank_src,irank_dest)%ie
               n = n + 1
               vecsnd(n) = var4d(i,j,i3,i4)
            enddo ; enddo ; enddo ; enddo

            call MPI_SEND( vecsnd(1:size_vecsnd), size_vecsnd, MPI_REAL, irank_dest, &
                            TAG, MPI_COMM_WORLD, ierr ) 

            if (ierr .ne. MPI_SUCCESS) then
               print *,'MSG mpi_sdsu_communicate_buf_4d: Error MPI_SEND. Terminating. myrank=',myrank
               call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
            endif

        endif

       if( myrank == irank_dest ) then

            size_vecrcv = bnd(4)*bnd(3)*overlap(irank_src,irank_dest)%totnum 

            call MPI_RECV( vecrcv(1:size_vecrcv), size_vecrcv, MPI_REAL, irank_src, &
                          TAG, MPI_COMM_WORLD, status, ierr )   

            if (ierr .ne. MPI_SUCCESS) then
               print *,'MSG mpi_sdsu_communicate_buf_4d: Error MPI_SEND. Terminating. myrank=',myrank
               call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
            endif

            n = 0
            do i4= 1,bnd(4)
            do i3= 1,bnd(3)
            do j = overlap(irank_src,irank_dest)%js, overlap(irank_src,irank_dest)%je
            do i = overlap(irank_src,irank_dest)%is, overlap(irank_src,irank_dest)%ie
               n = n + 1
               var4d(i,j,i3,i4) = vecrcv(n)
            enddo ; enddo ; enddo ; enddo

       endif

    endif


 enddo
 enddo


 deallocate( vecsnd, vecrcv )                

 end subroutine mpi_sdsu_communicate_buf_4d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine buffer_overlapp( m, b, o )
 implicit none
 type(memory_bounds),intent(in) :: m  !memory grid bounds
 type(memory_bounds),intent(in) :: b  !buffer grid bounds
 type(memory_bounds),intent(out) :: o  !buffer grid bounds
 
 integer :: im,jm
 integer :: cnt_overlap, cnt_max
 integer :: totlen
 integer,allocatable,dimension(:) :: ivec,jvec

 totlen = (m%je - m%js +1) * (m%ie - m%is +1) 
 allocate( ivec(totlen ),jvec(totlen) )
 !
 ! First determine memory and buffer grids overlapps or not.
 !
 cnt_overlap = 0     
 do jm = m%js,m%je
   do im = m%is,m%ie

       if( (jm >= b%js .and. jm <= b%je) .and. (im >= b%is .and. im <= b%ie) ) then
            cnt_overlap = cnt_overlap + 1
            ivec( cnt_overlap ) = im
            jvec( cnt_overlap ) = jm 
       endif

   enddo
 enddo

 cnt_max = cnt_overlap

 if( cnt_max > 0 ) then

!
! overlapped domain
!
   o%is = minval( ivec( 1:cnt_max ) )
   o%ie = maxval( ivec( 1:cnt_max ) )
   o%js = minval( jvec( 1:cnt_max ) )
   o%je = maxval( jvec( 1:cnt_max ) )
   o%totnum = cnt_max

 else

   o%is = -999
   o%ie = -999
   o%js = -999
   o%je = -999
   o%totnum = 0

 endif

 deallocate( ivec, jvec )

 return
 end subroutine buffer_overlapp

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_defined_2d( method, undefined, var2d  )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 2d array to mater (slave) processor via MPI library.
!   Note that about var2d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,intent(in)                     :: undefined     ! undefined value 
 real,dimension(:,:), intent(inout)  :: var2d   !assumed-shape 2d array

 real,allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer, allocatable,dimension(:) :: i1_vecsnd, i2_vecsnd
 integer, allocatable,dimension(:) :: i1_vecrcv, i2_vecrcv 
 integer :: bnd(2)  !upper memory bound for 3D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2 !1~2nd memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 
 integer :: cnt

 allocate(totlen(0:numproc_tot-1))
 allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(var2d)

!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 cnt = 0
 do i2 = 1,bnd(2) ; do i1 = 1,bnd(1)
    if( var2d(i1,i2) /= undefined ) cnt = cnt + 1
 enddo ; enddo

!
! distribute defined-only vector size to other memory. 
!
 call MPI_ALLGATHER(   cnt, 1, mpi_integer,                  &
                    totlen, 1, mpi_integer, mpi_comm_world, ierr)


!
! determine memory bound
! 
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif


!
! total vector size for tiles
!
 size_vecsnd = cnt  !sending vector (tile)
 allocate(    vecsnd(0: size_vecsnd-1 ) ,&
           i1_vecsnd(0: size_vecsnd-1 ) ,&
           i2_vecsnd(0: size_vecsnd-1 ) ,&
           stat=ierr )

! size_vecrcv = SIZE( var2d( 1:bnd(1) , 1:bnd(2)  ) )  !receiving vector (domain)
 size_vecrcv = sum( totlen )
 allocate(     vecrcv(0: size_vecrcv-1 ) ,&
            i1_vecrcv(0: size_vecrcv-1 ) ,&
            i2_vecrcv(0: size_vecrcv-1 ) ,&
            stat=ierr )


!
! Vectorize var2d
!
 n = 0
 do i2 = 1,bnd(2) ; do i1 = 1,bnd(1)
    if( var2d(i1,i2) /= undefined ) then
         vecsnd(n) = var2d(i1,i2)  !sending vector
      i1_vecsnd(n) = i1
      i2_vecsnd(n) = i2
      n=n+1
    endif
 enddo ; enddo


 if( trim(method) == 'TO_MASTER' ) then
!
! send vector to the master thread
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                   vecrcv, totlen, disp, mpi_real, 0, mpi_comm_world, ierr)

  call MPI_GATHERV(i1_vecsnd, size_vecsnd, mpi_integer,                  &
                   i1_vecrcv, totlen, disp, mpi_integer, 0, mpi_comm_world, ierr)

  call MPI_GATHERV(i2_vecsnd, size_vecsnd, mpi_integer,                  &
                   i2_vecrcv, totlen, disp, mpi_integer, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_defined_2d: Error MPI_GATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 elseif( trim(method) == 'TO_ALL' ) then
!
! send vector to all threads
!
  call MPI_ALLGATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                      vecrcv, totlen, disp, mpi_real, mpi_comm_world, ierr)

  call MPI_ALLGATHERV(i1_vecsnd, size_vecsnd, mpi_integer,                  &
                      i1_vecrcv, totlen, disp, mpi_integer, mpi_comm_world, ierr)

  call MPI_ALLGATHERV(i2_vecsnd, size_vecsnd, mpi_integer,                  &
                      i2_vecrcv, totlen, disp, mpi_integer, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_defined_2d: Error MPI_ALLGATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 endif

!
! Vector -> 2d Array
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then

     do n = 0 , size_vecrcv-1
        i1 = i1_vecrcv(n) 
        i2 = i2_vecrcv(n)
        var2d(i1,i2) =  vecrcv(n)
     enddo

 endif

!
! Deallocate all
! 
 deallocate( vecsnd, vecrcv, totlen, disp, &
             i1_vecsnd, i2_vecsnd, i1_vecrcv, i2_vecrcv   )


 return
 end subroutine mpi_sdsu_communicate_defined_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_defined_3d( method, undefined, var3d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 3d array to mater (slave) processor via MPI library.
!   Note that about var3d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,intent(in)                     :: undefined     ! undefined value 
 real,dimension(:,:,:), intent(inout)  :: var3d   !assumed-shape 3d array

 real,allocatable,dimension(:),save :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer, allocatable,dimension(:) :: i1_vecsnd, i2_vecsnd, i3_vecsnd
 integer, allocatable,dimension(:) :: i1_vecrcv, i2_vecrcv, i3_vecrcv
 integer :: bnd(3)  !upper memory bound for 3D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2,i3 !1~3th memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 
 integer :: cnt

 allocate(totlen(0:numproc_tot-1))
 allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(var3d)

!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 cnt = 0
 do i3 = 1, bnd(3) ;  do i2 = 1,bnd(2) ; do i1 = 1,bnd(1)
    if( var3d(i1,i2,i3) /= undefined ) cnt = cnt + 1
 enddo ; enddo ; enddo

!
! distribute defined-only vector size to other memory. 
!
 call MPI_ALLGATHER(   cnt, 1, mpi_integer,                  &
                    totlen, 1, mpi_integer, mpi_comm_world, ierr)

!
! determine memory bound
! 
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif

!
! total vector size for tiles
!
 size_vecsnd = cnt  !sending vector (tile)
 allocate(    vecsnd(0: size_vecsnd-1 ) ,&
           i1_vecsnd(0: size_vecsnd-1 ) ,&
           i2_vecsnd(0: size_vecsnd-1 ) ,&
           i3_vecsnd(0: size_vecsnd-1 ) ,&
           stat=ierr )

 size_vecrcv = sum( totlen )
 allocate(     vecrcv(0: size_vecrcv-1 ) ,&
            i1_vecrcv(0: size_vecrcv-1 ) ,&
            i2_vecrcv(0: size_vecrcv-1 ) ,&
            i3_vecrcv(0: size_vecrcv-1 ) ,&
            stat=ierr )

!
! Vectorize var3d
!
 n = 0
 do i3 = 1,bnd(3) ; do i2 = 1,bnd(2) ; do i1 = 1,bnd(1)
    if( var3d(i1,i2,i3) /= undefined ) then
         vecsnd(n) = var3d(i1,i2,i3)  !sending vector
      i1_vecsnd(n) = i1
      i2_vecsnd(n) = i2
      i3_vecsnd(n) = i3
      n=n+1
    endif
 enddo ; enddo ; enddo


 if( trim(method) == 'TO_MASTER' ) then
!
! send vector to the master thread
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                   vecrcv, totlen, disp, mpi_real, 0, mpi_comm_world, ierr)

  call MPI_GATHERV(i1_vecsnd, size_vecsnd, mpi_integer,                  &
                   i1_vecrcv, totlen, disp, mpi_integer, 0, mpi_comm_world, ierr)

  call MPI_GATHERV(i2_vecsnd, size_vecsnd, mpi_integer,                  &
                   i2_vecrcv, totlen, disp, mpi_integer, 0, mpi_comm_world, ierr)

  call MPI_GATHERV(i3_vecsnd, size_vecsnd, mpi_integer,                  &
                   i3_vecrcv, totlen, disp, mpi_integer, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_defined_3d: Error MPI_GATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 elseif( trim(method) == 'TO_ALL' ) then
!
! send vector to all threads
!
  call MPI_ALLGATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                      vecrcv, totlen, disp, mpi_real, mpi_comm_world, ierr)

  call MPI_ALLGATHERV(i1_vecsnd, size_vecsnd, mpi_integer,                  &
                      i1_vecrcv, totlen, disp, mpi_integer, mpi_comm_world, ierr)

  call MPI_ALLGATHERV(i2_vecsnd, size_vecsnd, mpi_integer,                  &
                      i2_vecrcv, totlen, disp, mpi_integer, mpi_comm_world, ierr)

  call MPI_ALLGATHERV(i3_vecsnd, size_vecsnd , mpi_integer,                  &
                      i3_vecrcv, totlen, disp, mpi_integer, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_defined_3d: Error MPI_ALLGATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 endif


!
! Vector -> 3d Array
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then

     do n = 0 , size_vecrcv-1
        i1 = i1_vecrcv(n)
        i2 = i2_vecrcv(n)
        i3 = i3_vecrcv(n)
        var3d(i1,i2,i3) = vecrcv(n)
     enddo

 endif

!
! Deallocate all
! 
 deallocate( vecsnd, vecrcv, totlen, disp, &
             i1_vecsnd, i2_vecsnd, i3_vecsnd, i1_vecrcv, i2_vecrcv, i3_vecrcv  )


 return
 end subroutine mpi_sdsu_communicate_defined_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_defined_4d( method, undefined, var4d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 4d array to mater (slave) processor via MPI library.
!   Note that about var4d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!   Makesure 3rd and 4th dimensions are identical between var4d and out4d
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,intent(in)                     :: undefined     ! undefined value 
 real,dimension(:,:,:,:), intent(inout)  :: var4d   !assumed-shape 4d array

 real,allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer, allocatable,dimension(:) :: i1_vecsnd, i2_vecsnd, i3_vecsnd, i4_vecsnd
 integer, allocatable,dimension(:) :: i1_vecrcv, i2_vecrcv, i3_vecrcv, i4_vecrcv
 integer :: bnd(4)  !upper memory bound for 4D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2,i3,i4 !1~4th memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 
 integer :: cnt

 allocate(totlen(0:numproc_tot-1))
 allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 4D array
!
 bnd = UBOUND(var4d)

!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 cnt = 0
 do i4 = 1, bnd(4) ;  do i3 = 1, bnd(3) ;  do i2 = 1,bnd(2) ; do i1 = 1,bnd(1)
    if( var4d(i1,i2,i3,i4) /= undefined ) cnt = cnt + 1
 enddo ; enddo ; enddo ; enddo

!
! distribute defined-only vector size to other memory. 
!
 call MPI_ALLGATHER(   cnt, 1, mpi_integer,                  &
                    totlen, 1, mpi_integer, mpi_comm_world, ierr)

!
! determine memory bound
! 
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif


!
! total vector size for tiles
!
 size_vecsnd = cnt  !sending vector (tile)
 allocate(    vecsnd(0: size_vecsnd-1 ) ,&
           i1_vecsnd(0: size_vecsnd-1 ) ,&
           i2_vecsnd(0: size_vecsnd-1 ) ,&
           i3_vecsnd(0: size_vecsnd-1 ) ,&
           i4_vecsnd(0: size_vecsnd-1 ) ,&
           stat=ierr )

 size_vecrcv = sum( totlen )
 allocate(     vecrcv(0: size_vecrcv-1 ) ,&
            i1_vecrcv(0: size_vecrcv-1 ) ,&
            i2_vecrcv(0: size_vecrcv-1 ) ,&
            i3_vecrcv(0: size_vecrcv-1 ) ,&
            i4_vecrcv(0: size_vecrcv-1 ) ,&
            stat=ierr )

!
! Vectorize var3d
!
 n = 0
 do i4 = 1,bnd(4) ;  do i3 = 1,bnd(3) ; do i2 = 1,bnd(2) ; do i1 = 1,bnd(1)
    if( var4d(i1,i2,i3,i4) /= undefined ) then
         vecsnd(n) = var4d(i1,i2,i3,i4)  !sending vector
      i1_vecsnd(n) = i1
      i2_vecsnd(n) = i2
      i3_vecsnd(n) = i3
      i4_vecsnd(n) = i4
      n=n+1
    endif
 enddo ; enddo ; enddo ; enddo

 if( trim(method) == 'TO_MASTER' ) then
!
! send vector to the master thread
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                   vecrcv, totlen, disp, mpi_real, 0, mpi_comm_world, ierr)

  call MPI_GATHERV(i1_vecsnd, size_vecsnd, mpi_integer,                  &
                   i1_vecrcv, totlen, disp, mpi_integer, 0, mpi_comm_world, ierr)

  call MPI_GATHERV(i2_vecsnd, size_vecsnd, mpi_integer,                  &
                   i2_vecrcv, totlen, disp, mpi_integer, 0, mpi_comm_world, ierr)

  call MPI_GATHERV(i3_vecsnd, size_vecsnd, mpi_integer,                  &
                   i3_vecrcv, totlen, disp, mpi_integer, 0, mpi_comm_world, ierr)

  call MPI_GATHERV(i4_vecsnd, size_vecsnd, mpi_integer,                  &
                   i4_vecrcv, totlen, disp, mpi_integer, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_defined_4d: Error MPI_GATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 elseif( trim(method) == 'TO_ALL' ) then
!
! send vector to all threads
!
  call MPI_ALLGATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                      vecrcv, totlen, disp, mpi_real, mpi_comm_world, ierr)

  call MPI_ALLGATHERV(i1_vecsnd, size_vecsnd, mpi_integer,                  &
                      i1_vecrcv, totlen, disp, mpi_integer, mpi_comm_world, ierr)

  call MPI_ALLGATHERV(i2_vecsnd, size_vecsnd, mpi_integer,                  &
                      i2_vecrcv, totlen, disp, mpi_integer, mpi_comm_world, ierr)

  call MPI_ALLGATHERV(i3_vecsnd, size_vecsnd , mpi_integer,                  &
                      i3_vecrcv, totlen, disp, mpi_integer, mpi_comm_world, ierr)

  call MPI_ALLGATHERV(i4_vecsnd, size_vecsnd , mpi_integer,                  &
                      i4_vecrcv, totlen, disp, mpi_integer, mpi_comm_world, ierr)


  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_defined_4d: Error MPI_ALLGATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 endif



!
! Vector -> 3d Array
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then

     do n = 0 , size_vecrcv-1
        i1 = i1_vecrcv(n)
        i2 = i2_vecrcv(n)
        i3 = i3_vecrcv(n)
        i4 = i4_vecrcv(n)
        var4d(i1,i2,i3,i4) = vecrcv(n)
     enddo

 endif

!
! Deallocate all
! 
 deallocate( vecsnd, vecrcv, totlen, disp, &
             i1_vecsnd, i2_vecsnd, i3_vecsnd, i4_vecsnd, i1_vecrcv, i2_vecrcv, i3_vecrcv, i4_vecrcv  )


 return
 end subroutine mpi_sdsu_communicate_defined_4d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_2d( method, var2d  )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 2d array to mater (slave) processor via MPI library.
!   Note that about var2d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,dimension(:,:), intent(inout)  :: var2d   !assumed-shape 2d array

 real,allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer :: bnd(2)  !upper memory bound for 3D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2 !1~2nd memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 

 allocate(totlen(0:numproc_tot-1))
 allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 3D array
!
    bnd = UBOUND(var2d)

!
! total vector size for tiles
!
 size_vecsnd = SIZE( var2d( myi_start:myi_end , myj_start:myj_end ) )  !sending vector
 allocate(  vecsnd(0: size_vecsnd-1 )  )

 size_vecrcv = SIZE( var2d( 1:bnd(1) , 1:bnd(2)  ) )  !receiving vector
 allocate(  vecrcv(0: size_vecrcv-1 )  )

!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then
     do irank = 0, numproc_tot-1
        totlen(irank) = myilen(irank) * myjlen(irank)  !total length
     enddo
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif

!
! Vectorize var2d
!
 n = 0
 do i2 = myj_start, myj_end ; do i1 = myi_start, myi_end
    vecsnd(n) = var2d(i1,i2)  !sending vector
    n=n+1
 enddo ; enddo

 if( trim(method) == 'TO_MASTER' ) then
!
! send vector to the master thread
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                   vecrcv, totlen, disp, mpi_real, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_2d: Error MPI_GATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 elseif( trim(method) == 'TO_ALL' ) then
!
! send vector to all threads
!
  call MPI_ALLGATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                      vecrcv, totlen, disp, mpi_real, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_2d: Error MPI_ALLGATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 endif

!
! Vector -> 2d Array
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then

    do irank = 0, numproc_tot-1

       is = myistr(irank)
       ie = myistr(irank) + myilen(irank) - 1
       js = myjstr(irank)
       je = myjstr(irank) + myjlen(irank) - 1
       n=disp(irank)
       do i2 = js, je ; do i1 = is, ie
          var2d(i1,i2) = vecrcv(n)  !sending vector
          n=n+1
       enddo ; enddo

    enddo

 endif
 

 deallocate( vecsnd, vecrcv, totlen, disp)

 return
 end subroutine mpi_sdsu_communicate_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_3d( method, var3d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 3d array to mater (slave) processor via MPI library.
!   Note that about var3d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!   Makesure 3rd dimension is identical between var3d and out3d
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,dimension(:,:,:), intent(inout)  :: var3d   !assumed-shape 3d array

 real,allocatable,dimension(:),save :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer :: bnd(3)  !upper memory bound for 3D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2,i3 !1~3th memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 

 allocate(totlen(0:numproc_tot-1))
 allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(var3d)

!
! total vector size for tiles
!
 size_vecsnd = SIZE( var3d( myi_start:myi_end , myj_start:myj_end , 1:bnd(3) ) )  !sending vector
 allocate(  vecsnd(0: size_vecsnd-1 )  ) 

 size_vecrcv = SIZE( var3d( 1:bnd(1) , 1:bnd(2) , 1:bnd(3) ) )  !receiving vector
 allocate(  vecrcv(0: size_vecrcv-1 )  )


!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or.  &
     (trim(method) == 'TO_ALL')   ) then
     do irank = 0, numproc_tot-1  
        totlen(irank) = myilen(irank) * myjlen(irank) * bnd(3)  !total length
     enddo
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif

!
! Vectorize var3d
!
 n = 0
 do i3 = 1, bnd(3) ; do i2 = myj_start, myj_end ; do i1 = myi_start, myi_end
    vecsnd(n) = var3d(i1,i2,i3)  !sending vector
    n=n+1
 enddo ; enddo ; enddo


 if( trim(method) == 'TO_MASTER' ) then

!
! send vector to the master thread
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                   vecrcv, totlen, disp, mpi_real, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_3d: Error MPI_GATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 elseif( trim(method) == 'TO_ALL' ) then

!
! send vector to all threads
!
  call MPI_ALLGATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                      vecrcv, totlen, disp, mpi_real, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_3d: Error MPI_ALLGATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 endif

!
! Vector -> 3d Array
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then

    do irank = 0, numproc_tot-1

       is = myistr(irank) 
       ie = myistr(irank) + myilen(irank) - 1
       js = myjstr(irank) 
       je = myjstr(irank) + myjlen(irank) - 1
       n=disp(irank)
       do i3 = 1, bnd(3) ; do i2 = js, je ; do i1 = is, ie 
          var3d(i1,i2,i3) = vecrcv(n)  !sending vector
          n=n+1
       enddo ; enddo ; enddo

    enddo

 endif

 deallocate( vecsnd, vecrcv, totlen, disp)

 return
 end subroutine mpi_sdsu_communicate_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_4d( method, var4d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 4d array to mater (slave) processor via MPI library.
!   Note that about var4d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!   Makesure 3rd and 4th dimensions are identical between var4d and out4d
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real,dimension(:,:,:,:), intent(inout)  :: var4d   !assumed-shape 4d array

 real,allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer :: bnd(4)  !upper memory bound for 4D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2,i3,i4 !1~4th memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 

 allocate(totlen(0:numproc_tot-1))
 allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 4D array
!
 bnd = UBOUND(var4d)

!
! total vector size for tiles
!
 size_vecsnd = SIZE( var4d( myi_start:myi_end , myj_start:myj_end , 1:bnd(3) , 1:bnd(4) ) )  !sending vector
 allocate(  vecsnd(0: size_vecsnd-1 )  ) 

 size_vecrcv = SIZE( var4d( 1:bnd(1) , 1:bnd(2) , 1:bnd(3) , 1:bnd(4) ) )  !receiving vector
 allocate(  vecrcv(0: size_vecrcv-1 )  )

!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or.  &
     (trim(method) == 'TO_ALL')   ) then
     do irank = 0, numproc_tot-1  
        totlen(irank) = myilen(irank) * myjlen(irank) * bnd(3) * bnd(4)  !total length
     enddo
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif

!
! Vectorize var4d
!
 n = 0
 do i4 = 1, bnd(4) ; do i3 = 1, bnd(3) ; do i2 = myj_start, myj_end ; do i1 = myi_start, myi_end
    vecsnd(n) = var4d(i1,i2,i3,i4)  !sending vector
    n=n+1
 enddo ; enddo ; enddo ; enddo


 if( trim(method) == 'TO_MASTER' ) then
!
! send vector to the master thread
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                   vecrcv, totlen, disp, mpi_real, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_4d: Error MPI_GATHERV in subroutine xy_gather. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 elseif( trim(method) == 'TO_ALL' ) then

!
! send vector to all threads
!
  call MPI_ALLGATHERV(vecsnd, size_vecsnd, mpi_real,                  &
                      vecrcv, totlen, disp, mpi_real, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_2d: Error MPI_ALLGATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 endif

!
! Vector -> 4d Array
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then

    do irank = 0, numproc_tot-1

       is = myistr(irank) 
       ie = myistr(irank) + myilen(irank) - 1
       js = myjstr(irank) 
       je = myjstr(irank) + myjlen(irank) - 1
       n=disp(irank)
       do i4 = 1, bnd(4) ; do i3 = 1, bnd(3) ; do i2 = js, je ; do i1 = is, ie 
          var4d(i1,i2,i3,i4) = vecrcv(n)  !sending vector
          n=n+1
       enddo ; enddo ; enddo ; enddo

    enddo

 endif

 deallocate( vecsnd, vecrcv, totlen, disp)

 return
 end subroutine mpi_sdsu_communicate_4d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_2d_dble( method, var2d  )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 2d array to mater (slave) processor via MPI library.
!   Note that about var2d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real(sdsu_fpd),dimension(:,:), intent(inout)  :: var2d   !assumed-shape 2d array

 real(sdsu_fpd),allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer :: bnd(2)  !upper memory bound for 3D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2 !1~2nd memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 

 allocate(totlen(0:numproc_tot-1))
 allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 3D array
!
    bnd = UBOUND(var2d)

!
! total vector size for tiles
!
 size_vecsnd = SIZE( var2d( myi_start:myi_end , myj_start:myj_end ) )  !sending vector
 allocate(  vecsnd(0: size_vecsnd-1 )  )

 size_vecrcv = SIZE( var2d( 1:bnd(1) , 1:bnd(2)  ) )  !receiving vector
 allocate(  vecrcv(0: size_vecrcv-1 )  )

!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then
     do irank = 0, numproc_tot-1
        totlen(irank) = myilen(irank) * myjlen(irank)  !total length
     enddo
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif

!
! Vectorize var2d
!
 n = 0
 do i2 = myj_start, myj_end ; do i1 = myi_start, myi_end
    vecsnd(n) = var2d(i1,i2)  !sending vector
    n=n+1
 enddo ; enddo

 if( trim(method) == 'TO_MASTER' ) then
!
! send vector to the master thread
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real8,                  &
                   vecrcv, totlen, disp, mpi_real8, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_2d_dble: Error MPI_GATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 elseif( trim(method) == 'TO_ALL' ) then
!
! send vector to all threads
!
  call MPI_ALLGATHERV(vecsnd, size_vecsnd, mpi_real8,                  &
                      vecrcv, totlen, disp, mpi_real8, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_2d_dble: Error MPI_ALLGATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 endif

!
! Vector -> 2d Array
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then

    do irank = 0, numproc_tot-1

       is = myistr(irank)
       ie = myistr(irank) + myilen(irank) - 1
       js = myjstr(irank)
       je = myjstr(irank) + myjlen(irank) - 1
       n=disp(irank)
       do i2 = js, je ; do i1 = is, ie
          var2d(i1,i2) = vecrcv(n)  !sending vector
          n=n+1
       enddo ; enddo

    enddo

 endif
 

 deallocate( vecsnd, vecrcv, totlen, disp)

 return
 end subroutine mpi_sdsu_communicate_2d_dble

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_3d_dble( method, var3d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 3d array to mater (slave) processor via MPI library.
!   Note that about var3d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!   Makesure 3rd dimension is identical between var3d and out3d
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real(sdsu_fpd),dimension(:,:,:), intent(inout)  :: var3d   !assumed-shape 3d array

 real(sdsu_fpd),allocatable,dimension(:),save :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer :: bnd(3)  !upper memory bound for 3D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2,i3 !1~3th memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 

 allocate(totlen(0:numproc_tot-1))
 allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 3D array
!
 bnd = UBOUND(var3d)

!
! total vector size for tiles
!
 size_vecsnd = SIZE( var3d( myi_start:myi_end , myj_start:myj_end , 1:bnd(3) ) )  !sending vector
 allocate(  vecsnd(0: size_vecsnd-1 )  ) 

 size_vecrcv = SIZE( var3d( 1:bnd(1) , 1:bnd(2) , 1:bnd(3) ) )  !receiving vector
 allocate(  vecrcv(0: size_vecrcv-1 )  )


!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or.  &
     (trim(method) == 'TO_ALL')   ) then
     do irank = 0, numproc_tot-1  
        totlen(irank) = myilen(irank) * myjlen(irank) * bnd(3)  !total length
     enddo
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif

!
! Vectorize var3d
!
 n = 0
 do i3 = 1, bnd(3) ; do i2 = myj_start, myj_end ; do i1 = myi_start, myi_end
    vecsnd(n) = var3d(i1,i2,i3)  !sending vector
    n=n+1
 enddo ; enddo ; enddo


 if( trim(method) == 'TO_MASTER' ) then

!
! send vector to the master thread
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real8,                  &
                   vecrcv, totlen, disp, mpi_real8, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_3d_dble: Error MPI_GATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 elseif( trim(method) == 'TO_ALL' ) then

!
! send vector to all threads
!
  call MPI_ALLGATHERV(vecsnd, size_vecsnd, mpi_real8,                  &
                      vecrcv, totlen, disp, mpi_real8, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_3d_dble: Error MPI_ALLGATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 endif

!
! Vector -> 3d Array
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then

    do irank = 0, numproc_tot-1

       is = myistr(irank) 
       ie = myistr(irank) + myilen(irank) - 1
       js = myjstr(irank) 
       je = myjstr(irank) + myjlen(irank) - 1
       n=disp(irank)
       do i3 = 1, bnd(3) ; do i2 = js, je ; do i1 = is, ie 
          var3d(i1,i2,i3) = vecrcv(n)  !sending vector
          n=n+1
       enddo ; enddo ; enddo

    enddo

 endif

 deallocate( vecsnd, vecrcv, totlen, disp)

 return
 end subroutine mpi_sdsu_communicate_3d_dble

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_communicate_4d_dble( method, var4d )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!   This program send/distribute 4d array to mater (slave) processor via MPI library.
!   Note that about var4d 1st memory must be i direction, 2nd memory must be j direction.
!   For input var, each dimension is assumed to start from 1.
!   Makesure 3rd and 4th dimensions are identical between var4d and out4d
!
! History:
! 10/2009 Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in)         :: method  !MPI colletable method
 real(sdsu_fpd),dimension(:,:,:,:), intent(inout)  :: var4d   !assumed-shape 4d array

 real(sdsu_fpd),allocatable,dimension(:) :: vecsnd, vecrcv !temporal sending & receiving vector
 integer, allocatable,dimension(:) :: totlen, disp !vector length and disp points
 integer :: bnd(4)  !upper memory bound for 4D array
 integer :: size_vecsnd, size_vecrcv  !size of sending/receiving vector
 integer :: i1,i2,i3,i4 !1~4th memory loop
 integer :: n           !vector count
 integer :: irank       !myrank loop
 integer :: is,ie,js,je !starting/ending i j indices 

 allocate(totlen(0:numproc_tot-1))
 allocate(disp(0:numproc_tot-1))

!
! find upper bound from assumped 4D array
!
 bnd = UBOUND(var4d)

!
! total vector size for tiles
!
 size_vecsnd = SIZE( var4d( myi_start:myi_end , myj_start:myj_end , 1:bnd(3) , 1:bnd(4) ) )  !sending vector
 allocate(  vecsnd(0: size_vecsnd-1 )  ) 

 size_vecrcv = SIZE( var4d( 1:bnd(1) , 1:bnd(2) , 1:bnd(3) , 1:bnd(4) ) )  !receiving vector
 allocate(  vecrcv(0: size_vecrcv-1 )  )

!
! determin memory bound of vector for each processor (memory index of vector start from 0)
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or.  &
     (trim(method) == 'TO_ALL')   ) then
     do irank = 0, numproc_tot-1  
        totlen(irank) = myilen(irank) * myjlen(irank) * bnd(3) * bnd(4)  !total length
     enddo
     disp(0) = 0
     do irank = 1, numproc_tot-1
        disp(irank) = disp(irank-1)+totlen(irank-1)  !disp points
     enddo
 endif

!
! Vectorize var4d
!
 n = 0
 do i4 = 1, bnd(4) ; do i3 = 1, bnd(3) ; do i2 = myj_start, myj_end ; do i1 = myi_start, myi_end
    vecsnd(n) = var4d(i1,i2,i3,i4)  !sending vector
    n=n+1
 enddo ; enddo ; enddo ; enddo


 if( trim(method) == 'TO_MASTER' ) then
!
! send vector to the master thread
!
  call MPI_GATHERV(vecsnd, size_vecsnd, mpi_real8,                  &
                   vecrcv, totlen, disp, mpi_real8, 0, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_4d_dble: Error MPI_GATHERV in subroutine xy_gather. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 elseif( trim(method) == 'TO_ALL' ) then

!
! send vector to all threads
!
  call MPI_ALLGATHERV(vecsnd, size_vecsnd, mpi_real8,                  &
                      vecrcv, totlen, disp, mpi_real8, mpi_comm_world, ierr)

  if (ierr .ne. MPI_SUCCESS) then
      print *,'MSG mpi_sdsu_communicate_2d_dble: Error MPI_ALLGATHERV. Terminating. myrank=',myrank
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  endif

 endif

!
! Vector -> 4d Array
!
 if( (trim(method) == 'TO_MASTER' .and. masterproc) .or. &
     (trim(method) == 'TO_ALL')   ) then

    do irank = 0, numproc_tot-1

       is = myistr(irank) 
       ie = myistr(irank) + myilen(irank) - 1
       js = myjstr(irank) 
       je = myjstr(irank) + myjlen(irank) - 1
       n=disp(irank)
       do i4 = 1, bnd(4) ; do i3 = 1, bnd(3) ; do i2 = js, je ; do i1 = is, ie 
          var4d(i1,i2,i3,i4) = vecrcv(n)  !sending vector
          n=n+1
       enddo ; enddo ; enddo ; enddo

    enddo

 endif

 deallocate( vecsnd, vecrcv, totlen, disp)

 return
 end subroutine mpi_sdsu_communicate_4d_dble

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_wait_for_master
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Let slave processors wait for master processor.
!
! History:
! 09/2009  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=1) ::  inmsg

   call MPI_BCAST(inmsg, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

 return
 end subroutine mpi_wait_for_master

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_abort(msg,irank)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Abort MPI process
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------
 character(len=*),intent(in) :: msg
 integer,intent(in) :: irank   

 integer :: ierror

 print*,'myrank=',irank,' ',msg
 call MPI_FINALIZE(ierror)
 stop 'Abort SDSU'

 end subroutine mpi_sdsu_abort

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine mpi_sdsu_finale
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: 
!  Finalize MPI process
!
! History:
! 06/2008  Toshi Matsui@NASA GSFC : Initial 
!
! References:
!  S. Vetter, Y. Aoyama, J. Nakano, 1999, RS/600 SP: Practical MPI Programming, 
!     IBM Redbooks, 1999, ISBN: 0738413658, 238p. 
!-----------------------------------------------------------------------------------------------------

 if(masterproc) print*,'MSG mpi_sdsu_finale: finalize MPI process'
  call MPI_FINALIZE(ierr)
 return
 end subroutine mpi_sdsu_finale

#endif 

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine stop_sdsu(irank,error_msg)
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!  Terminate simulator. 
! 
! History:
! 07/2014   Toshi Matsui@NASA GSFC : moved from module_simulator to module_mpi
! 11/2011  Toshi Matsui@NASA GSFC : Initial. 
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 integer,intent(in) :: irank
 character(len=*),intent(in) :: error_msg


#if MPI > 0
 !
 ! For MPI
 !
 call mpi_sdsu_abort(error_msg,irank)

#else
 !
 ! For non-MPI
 !
 print*,trim(error_msg)
 print*,'Abort SDSU'
 stop

#endif 

 end subroutine stop_sdsu

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine check(status,para,missing_para)
 implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Check status of reading process of netCDF. 
!
! History:
! 05/2014  Toshi Matsui@NASA GSFC ; included optional statement for skiping
! 05/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!---------------------------------------------------------------------------------------------------
  integer, intent (in) :: status
  character(len=*),intent(in),optional :: para  !optional parameter for missing, but skippable parameters.
  logical,intent(out),optional :: missing_para  !optional parameter for missing, but skippable parameters.

  if( present(missing_para) )  missing_para = .false.

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    if( present(para) ) then
          if(masterproc) print*,'Missing ',trim(para)
          if( present(missing_para) ) missing_para = .true.
          return
    endif
    call stop_sdsu(myrank, "Stopped missing parameters!"  )
  end if
 end subroutine check

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU

 subroutine check_io(status,bad_io)
  implicit none
!---------------------------------------------------------------------------------------------------
! Comments:  
! Check status of reading process of netCDF. 
!
! History:
! 05/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!---------------------------------------------------------------------------------------------------
  integer, intent (in) :: status
  logical, intent (out) :: bad_io
  bad_io = .false.
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    bad_io = .true.
  end if
 end subroutine check_io

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_mpi

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
