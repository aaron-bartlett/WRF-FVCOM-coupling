 module find_thread_number
 implicit none

 public   ! all variables and subourtines are accessible in module_simulator.

 contains

 subroutine mpi_sdsu_tile(mode,imin,imax,jmin,jmax,kmin,kmax)
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
 integer,intent(in) :: mode  ! 1 ~ 4
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
 integer,parameter :: max_cpu = 28
 integer :: node_number(1:max_cpu)
 integer :: total_cpu
 character(len=2) :: cc
 character(len=4) :: nnnn
 character(len=5) :: ttttt


 tot_grid_i = (imax-imin+1)
 tot_grid_j = (jmax-jmin+1)
 tot_grid = tot_grid_i * tot_grid_j  !total grid number


 if    (mode == 1) then ! debug
    total_cpu = 32
 elseif(mode == 2) then ! general_small
    total_cpu = 16
 elseif(mode == 3) then ! general
    total_cpu = 1024
 elseif(mode == 4) then ! general_hig
    total_cpu = 3072
 else
  stop 'no such mode'
 endif

 write(*,*), ''
 write(*,*), ' Find appropriate node and processor number for DISCOVER/PLEIADES cluster'
 write(*,*), ''
! write(*,FMT='(A107)'), '       ncpus=1 ncpus=2 ncpus=3 ncpus=4 ncpus=5 ncpus=6 ncpus=7 ncpus=8 ncpus=9 ncpus=10 ncpus=11 ncpus=12'

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
        do nt = 1, max_cpu
           resid = MOD( numproc_tot, nt )
           if(resid /= 0)  then
              node_number(nt) = 0
           else
              node_number(nt) = numproc_tot/nt
              write(cc,"(I2)") nt
              write(nnnn,"(I4)") node_number(nt)
              write(ttttt,"(I5)") numproc_tot
              write(*,*) cc//'CPU per node X '//nnnn//'nodes = '//ttttt//' total CPUs' 
           endif
        enddo
!        write(*,FMT='(A8,12I8)'),'node # =', node_number(:)
        

      endif

    endif

 enddo

 return
 end subroutine mpi_sdsu_tile
 end module find_thread_number

 program main
 use find_thread_number
 implicit none
 integer :: mode, imin,imax,jmin,jmax,kmin,kmax

 print*,'put qsub mode 1~4'
 print*,'1-debug  2-general_small 3-general 4-general_hig'
 read(*,*) mode
 imin=1 ; jmin=1 ; kmin=1 ; kmax=60

 print*,'Please put domain size'
 print*,'imax = '
 read(*,*) imax
 print*,'jmax = '
 read(*,*) jmax

 call  mpi_sdsu_tile(mode, imin,imax,jmin,jmax,kmin,kmax)

 end program main

