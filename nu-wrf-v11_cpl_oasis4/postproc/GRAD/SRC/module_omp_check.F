#include <define.h>

 module module_omp_check

 contains

!******************************************************************************
! FILE: omp_getEnvInfo.f
! DESCRIPTION:
!   OpenMP Example - Get Environment Information - Fortran Version
!   The master thread queries and prints selected environment information.
! AUTHOR: Blaise Barney  7/06
! LAST REVISED: 07/12/06
!******************************************************************************
 
 subroutine getinfo_omp

 integer :: nthreads, tid, omp_get_num_threads, &
        omp_get_thread_num, omp_get_num_procs, omp_get_max_threads, &
        omp_in_parallel, omp_get_dynamic, omp_get_nested, &
        procs, maxt

!     These are for AIX compilations
!     INTEGER INPAR, DYNAMIC, NESTED
!     These are for non-AIX compilations
 logical inpar, dynamic, nested

#ifdef OMP

!     Start parallel region
!$OMP PARALLEL PRIVATE(nthreads, tid) 

!     Obtain thread number
   tid = omp_get_thread_num()

!     Only master thread does this
      if (tid .eq. 0) then

        print *, 'thread',tid,'getting environment information'

!     Get environment information
        procs = omp_get_num_procs()
        nthreads = omp_get_num_threads()
        maxt = omp_get_max_threads()
        inpar = omp_in_parallel()
        dynamic = omp_get_dynamic()
        nested = omp_get_nested()

!     Print environment information

        print *, 'number of processors = ', procs
        print *, 'number of threads = ', nthreads
        print *, 'max threads = ', maxt
        print *, 'in parallel? = ', inpar
        print *, 'dynamic threads enabled? = ', dynamic
        print *, 'nested parallelism supported? = ', nested

      end if

!     Done
!$OMP END PARALLEL

#endif

 end subroutine getinfo_omp

 end module module_omp_check
