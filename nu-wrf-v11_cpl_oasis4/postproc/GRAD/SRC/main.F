#include <define.h>

 program main
 use module_radiation_driver
 use module_ra_goddard_2014
 use module_ra_goddard_2011
 use module_runtime_param
 use module_mpi

!
! Main program that drives entire radiation testbed program
! Toshi Matsui @ NASA GSFC
!
 implicit none
 character(len=100),dimension(1000) :: file_list ! input file names
 integer :: n , nmax, nproc
 integer :: inp
! real(sdsu_fpd) :: ts, te  ! starting & ending wall-clock time
 logical :: exit_proc  !if true exit entire rest of the processes. 
 real :: ts_rad, te_rad

!
!  read run-time parameter
!
  call read_runtime_param

!
! read input file name list
!
   print*,'MSG main reading '//'inpfile_names'
   open(unit=99,file='inpfile_names',status='old') !Open a namelist file for CRM input
   n = 0
   do !forever 
     n = n + 1
     read(unit=99,FMT='(A100)',END=999) file_list(n)  !Read CRM input file name
   enddo !forever
999 nmax = n - 1
   close(99)

!
! Define MPI decomposition parametrs 
!
  exit_proc = .false. !initialize as false
#if MPI == 1

 call mpi_sdsu_init  ! initialize MPI
 call mpi_sdsu_file(nmax,exit_proc)  ! soring input files for each processor
  myi_start = ids_wrf  ; myi_end   = ide_wrf
  myj_start = jds_wrf  ; myj_end   = jde_wrf
  myk_start = 1        ; myk_end   = kmax

#elif MPI == 2
 
 call mpi_sdsu_init  ! initialize MPI

 call mpi_sdsu_tile(ids_wrf,ide_wrf,jds_wrf,jde_wrf,1,kmax) ! sorting domain into tiles for each processor

 myn_start = 1  ; myn_end   = nmax

!
! OMP decomposition in case of MPI
!
 call mpi_sdsu_tile_omp( num_tiles_omp, myi_start, myi_end, myj_start, myj_end, myk_start, myk_end ) 


#else

  ! Default values of MPI parameters for single CPU run 
  masterproc = .true.
  myrank = 0
  myn_start = 1        ; myn_end   = nmax
  myi_start = ids_wrf  ; myi_end   = ide_wrf
  myj_start = jds_wrf  ; myj_end   = jde_wrf
  myk_start = 1        ; myk_end   = kmax

#endif


#if MPI > 0
  call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
#endif

!
! allocate all paramters
!
  call allocate_all


!
! initialize goddard radiation 2014 (LUT of single-scattering properties)
!
 if ( lw_physics==goddardlwscheme_2014 .or. sw_physics== goddardswscheme_2014 ) &
 call init_goddardrad_2014(.true.,  &
                           ids, ide, jds, jde, kds, kde,                &
                           ims, ime, jms, jme, kms, kme,                &
                           its, ite, jts, jte, kts, kte                 )



!------------- BEGIN FILE LOOP ---------------------------------

  file_loop: do n = myn_start, myn_end

    if( exit_proc ) cycle  !this ensure to cycle proces, when processor is not assigned for task in MPI=1 options


!  do n = 1, nmax
    CRM_file_list = trim(file_list(n))

    !
    ! read WRF output file (for input of radiation scheme)
    !
    call read_netcdf

    !
    ! call radiation driver
    !
    if(masterproc) call cpu_time(ts_rad)

    call radiation_driver

    if(masterproc) call cpu_time(te_rad)
    if(masterproc) print*,'MSG main: total computing time of radiation_driver is ',te_rad-ts_rad
    !
    ! write out simulated radiation flux and heating rate.
    !
    call write_out


  enddo file_loop


!-------------- END OF FILE LOOP -------------------------------------

#if MPI > 0
!
! terminate mpi process
!
 call mpi_sdsu_finale  !terminate mpi process
#endif


end program main

