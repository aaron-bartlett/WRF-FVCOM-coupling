MODULE module_casaco2_fluxes
   USE module_configure
   USE module_state_description

   IMPLICIT NONE

   ! 09/12/2013 This module contains parameterizations to calculate CASA CO2 fluxes:
   ! Author: Jossy P. Jacob / Tao Zhining 
   !

CONTAINS

   SUBROUTINE add_casaco2_fluxes (  ids,ide, jds,jde, kds,kde,             &     ! Domain dimensions
        ims,ime, jms,jme, kms,kme,             &     ! Memory dimensions 
        its,ite, jts,jte, kts,kte,             &     ! Tile dimensions
        dtstep, dz8w, config_flags, rho_phy,   &
        ktau, stepcasa,                        &  !add by Z. Tao
        chem, emis_ant )


      TYPE(grid_config_rec_type), INTENT(IN   )    :: config_flags

      INTEGER, INTENT(IN   ) ::   ids,ide, jds,jde, kds,kde,        &
           ims,ime, jms,jme, kms,kme,        &
           its,ite, jts,jte, kts,kte

      REAL, INTENT(IN   )    ::   dtstep ! ktau, stepcasa !add by Z. Tao
      INTEGER, INTENT(IN   )    ::   ktau, stepcasa !added by JJACOB

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN ) ::   rho_phy, dz8w
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), INTENT(INOUT ) ::      chem
      REAL, DIMENSION( ims:ime, kms:config_flags%kemit, jms:jme, num_emis_ant ), INTENT(INOUT ) :: emis_ant !modified by Z. Tao

      INTEGER :: i,j,k
      REAL    :: conv_rho
      character *256 :: message
      real :: chem_old ! EMK
!      print *, 'num_emis_ant, size(emis_ant), p_e_casaco2, p_et_casaco2 are '
!      print *, num_emis_ant, size(emis_ant), p_e_casaco2, p_et_casaco2
!       print *, 'dtstep,stepcasa,ktau ',dtstep, stepcasa, ktau
      ! For both CASACO2 options
      DO j=jts,jte
        DO i=its,ite
         !  3D anthropogenic fluxes
          DO k=kts,min(config_flags%kemit,kte )
            chem_old = chem(i,k,j,p_casaco2) ! EMK
            !write (message,'(A,3E12.4)') 'In module_casaco2_fluxes2', &
            ! emis_ant(i,k,j,p_e_casaco2),emis_ant(i,k,j,p_et_casaco2), chem(i,k,j,p_casaco2)
            !call wrf_message(message)
            !write (message,'(A,3E12.4)') 'In module_casaco2_fluxes2 p_e_casaco2, p_et_casaco2, size(emis_ant,4)', &
            ! p_e_casaco2, p_et_casaco2, size(emis_ant,4) 
            !call wrf_message(message)
            conv_rho=8.0461e-6/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)  ! 8.0461e-6=molar_mass(air)/3600, [g/mol/s]

! add by Z. Tao  (Tendency is per hour).
            if(MOD(ktau,stepcasa) .ne. 0) then
              emis_ant(i,k,j,p_e_casaco2) = emis_ant(i,k,j,p_e_casaco2) + (dtstep/3600.) * emis_ant(i,k,j,p_et_casaco2)
            end if

            chem(i,k,j,p_casaco2)= chem(i,k,j,p_casaco2) + conv_rho* emis_ant(i,k,j,p_e_casaco2)
            !write (message,'(A,I3.3,1X,I3.3,1X,I3.3,1X,3E12.4)') &
            !       'In module_casaco2_fluxes2,i,k,j,emis_ant,chem_old,chem ', &
            !        i,k,j,emis_ant(i,k,j,p_e_casaco2),chem_old,chem(i,k,j,p_casaco2)
            !call wrf_message(message)
            !write (message,'(A,I3.3,1X,I3.3,1X,I3.3,1X,4E12.4)') &
            !      'In module_casaco2_fluxes2,i,k,j,rho_phy,dtstep,dz8w,conv_rho ', &
            !       i,k,j,rho_phy(i,k,j),dtstep,dz8w(i,k,j),conv_rho
            !call wrf_message(message)
          ENDDO
        ENDDO
      ENDDO

!write (message,'(A,2I3.3,2E12.4)') 'In module_casaco2_fluxes2,its,jts', &
!     its,jts,emis_ant(its,1,jts,p_e_casaco2), chem(its,1,jts,p_casaco2)
!call wrf_message(message)

END SUBROUTINE add_casaco2_fluxes

END MODULE module_casaco2_fluxes
