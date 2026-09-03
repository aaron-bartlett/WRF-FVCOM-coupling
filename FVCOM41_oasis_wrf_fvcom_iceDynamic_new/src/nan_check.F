!--ch201902 added for nan check during simulation
SUBROUTINE nan_check(vname)

  !==============================================================================|
  USE LIMS
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_PAR  
  IMPLICIT NONE
  CHARACTER(LEN=*),intent(in)  :: vname
  INTEGER :: IERR
  !==============================================================================|

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: nan_check"

#  if defined (MULTIPROCESSOR)
   if (PAR) then
      CALL MPI_BARRIER(il_commlocal,IERR)
   end if
#  endif

   if (trim(vname)=='EL') then
      if(maxval(EL1)>30.0.or.minval(EL1)<-30.0) then
         print *,"maxval(EL1)",maxval(EL1)
         print *,"minval(EL1)",minval(EL1)
         print *,"error,check --> on,free surface elevation"
         CHECK_index(1)=1.0
      end if
   elseif (trim(vname)=='UV') then
      if(maxval(u)>30.0.or.minval(u)<-30.0.or.&
         maxval(v)>30.0.or.minval(v)<-30.0) then
         print *,"maxval(u)",maxval(u)
         print *,"minval(u)",minval(u)
         print *,"maxval(v)",maxval(v)
         print *,"minval(v)",minval(v)
         print *,"error,check --> on,u or v"
         CHECK_index(1)=1.0
      end if
   elseif (trim(vname)=='TEMP') then
      if(maxval(T)>50.0.or.minval(T)<-10.0) then
         print *,"maxval(T)",maxval(T)
         print *,"minval(T)",minval(T)
         print *,"error,check --> on,temperature"
         CHECK_index(1)=1.0
      end if
   elseif (trim(vname)=='SALINITY') then
      if(maxval(S)>50.0.or.minval(S)<-2.0) then
         print *,"maxval(S)",maxval(S)
         print *,"minval(S)",minval(S)
         print *,"error,check --> on,salinity"
         CHECK_index(1)=1.0
      end if
   end if

#  if defined (MULTIPROCESSOR)
   if (PAR) then
      CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,check_index,check_index_global)
      CALL MPI_BCAST(CHECK_index_global,MGL,MPI_F,0,il_commlocal,IERR)
      CHECK_ERROR=int(sum(CHECK_index_global))
   end if
#  endif

   if (trim(vname)=='check') then
    if (check_error>0) then
      print *,"CHECK_ERROR:", CHECK_ERROR
      print *,"error,call pstop"
      CALL FATAL_ERROR("find extreme value,stop before NaN")
      call pstop
    end if
   end if


   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: nan_check"
   
   return
   END SUBROUTINE nan_check
!==============================================================================|
