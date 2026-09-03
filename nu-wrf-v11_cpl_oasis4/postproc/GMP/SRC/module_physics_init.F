!WRF:MODEL_LAYER:INITIALIZATION
!

!  This MODULE holds the routines which are used to perform model start-up operations
!  for the individual domains.  This is the stage after inputting wrfinput and before
!  calling 'integrate'.

!  This MODULE CONTAINS the following routines:


MODULE module_physics_init

 use module_runtime_param
 USE module_mp_goddard_unified, only : goddard_microphysics_init

 type grid_config_rec_type
  integer :: mp_physics
 end type grid_config_rec_type

CONTAINS

!=====================================================================
   SUBROUTINE mp_init(config_flags, allowed_to_read, dt, z_at_q,  &
                      ids, ide, jds, jde, kds, kde,               &
                      ims, ime, jms, jme, kms, kme,               &
                      its, ite, jts, jte, kts, kte                )
!------------------------------------------------------------------
!   USE module_mp_wsm3
!------------------------------------------------------------------
   IMPLICIT NONE
!------------------------------------------------------------------
! Arguments
 TYPE (grid_config_rec_type) :: config_flags
 LOGICAL , INTENT(IN)        :: allowed_to_read
 REAL,INTENT(IN)             :: dt ! time step 
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN):: z_at_q  !height at q level (m)

   INTEGER , INTENT(IN)        :: ids, ide, jds, jde, kds, kde,   &
                                  ims, ime, jms, jme, kms, kme,   &
                                  its, ite, jts, jte, kts, kte


   mp_select: SELECT CASE(config_flags%mp_physics)

     CASE (NUWRF4ICESCHEME)

     CASE(GMPUNIFIED)
         call goddard_microphysics_init (allowed_to_read, dt, &
                       ids, ide, jds, jde, kds, kde,      &
                       ims, ime, jms, jme, kms, kme,      &
                       its, ite, jts, jte, kts, kte )


     CASE DEFAULT

   END SELECT mp_select


 END SUBROUTINE mp_init

END MODULE module_physics_init
