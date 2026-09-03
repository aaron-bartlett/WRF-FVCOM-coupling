#include <define.h>

!WRF:MEDIATION_LAYER:PHYSICS
! *** add new modules of schemes here
!
 MODULE module_microphysics_driver
 CONTAINS

 SUBROUTINE microphysics_driver( itimestep &
!                          ids,ide, jds,jde, kds,kde                &
!                         ,ims,ime, jms,jme, kms,kme                &
!                         ,its,ite, jts,jte, kts,kte                &
                                ) 

 use module_runtime_param
 use module_wrf_offline
 use module_mpi

!#ifdef DM_PARALLEL
!  USE module_dm, ONLY : &
!                 local_communicator, mytask,  wrf_dm_min_real, wrf_dm_max_real
!#endif

! Model Layer
!   USE module_model_constants
!   USE module_wrf_error
!   USE module_configure, only: grid_config_rec_type

!EMK NUWRF...For MPI safe CSI diagnostic
!#ifdef DM_PARALLEL
!   use module_comm_dm, only: halo_em_phys_mp_0_sub,halo_em_phys_mp_1_sub, &
!        halo_em_phys_mp_2_sub,halo_em_phys_mp_3_sub,halo_em_phys_mp_4_sub
!   use module_dm, only: ntasks_x,ntasks_y,local_communicator,mytask,ntasks
!#endif
!EMK NUWRF

! *** add new modules of schemes here
 USE module_mp_gsfcgce_4ice_nuwrf 
 USE module_mp_goddard_unified, only : goddard_microphysics_wrf

 
!   USE module_mp_gsfcgce_4ice_nuwrf, only: gsfcgce_4ice_nuwrf

! For checking model timestep is history time (for radar reflectivity)
!   USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm
!   USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep
!   use module_domain, only : domain ! EMK NUWRF...for MPI-safe CSI

   ! EMK NUWRF...MPI-safe CSI
!   use module_diag_csi, only: cal_sepa_part1, cal_sepa_part2, &
!        cal_sepa_part3, var_sepa

!----------------------------------------------------------------------
   ! This driver calls subroutines for the microphys.

!----------------------------------------------------------------------
   IMPLICIT NONE
!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Definitions
!-----------
! Rho_d      dry density (kg/m^3)
! Theta_m    moist potential temperature (K)
! Qv         water vapor    mixing ratio (kg/kg)
! Qc         cloud water    mixing ratio (kg/kg)
! Qr         rain water     mixing ratio (kg/kg)
! Qi         cloud ice      mixing ratio (kg/kg)
! Qs         snow           mixing ratio (kg/kg)
! Qg         graupel        mixing ratio (kg/kg)
! Qh         hail           mixing ratio (kg/kg)
! Qndrop     droplet number mixing ratio (#/kg)
! Qni        cloud ice number concentration (#/kg)
! Qns        snow      number concentration (#/kg)
! Qnr        rain      number concentration (#/kg)
! Qng        graupel   number concentration (#/kg)
! Qnh        hail      number concentration (#/kg)

! Qzr        rain             reflectivity (m6/kg)
! Qzi        ice              reflectivity (m6/kg)
! Qzs        snow             reflectivity (m6/kg)
! Qzg        graupel          reflectivity (m6/kg)
! Qzh        hail             reflectivity (m6/kg)

! Qvolg        graupel   particle volume (m3/kg)

!
!----------------------------------------------------------------------
!-- th        potential temperature    (K)
!-- moist_new     updated moisture array   (kg/kg)
!-- moist_old     Old moisture array       (kg/kg)
!-- rho           density of air           (kg/m^3)
!-- pi_phy        exner function           (dimensionless)
!-- p             pressure                 (Pa)
!-- RAINNC        grid scale precipitation (mm)
!-- RAINNCV       one time step grid scale precipitation (mm/step)
!-- SNOWNC        grid scale snow and ice (mm)
!-- SNOWNCV       one time step grid scale snow and ice (mm/step)
!-- GRAUPELNC     grid scale graupel (mm)
!-- GRAUPELNCV    one time step grid scale graupel (mm/step)
!-- HAILNC        grid scale hail (mm)
!-- HAILNCV       one time step grid scale hail (mm/step)
!-- SR            one time step mass ratio of snow to total precip
!-- z             Height above sea level   (m)
!-- dt            Time step              (s)
!-- G             acceleration due to gravity  (m/s^2)
!-- CP            heat capacity at constant pressure for dry air (J/kg/K)
!-- R_d           gas constant for dry air (J/kg/K)
!-- R_v           gas constant for water vapor (J/kg/K)
!-- XLS           latent heat of sublimation   (J/kg)
!-- XLV           latent heat of vaporization  (J/kg)
!-- XLF           latent heat of melting       (J/kg)
!-- rhowater      water density                      (kg/m^3)
!-- rhosnow       snow density               (kg/m^3)
!-- F_ICE_PHY     Fraction of ice.
!-- F_RAIN_PHY    Fraction of rain.
!-- F_RIMEF_PHY   Mass ratio of rimed ice (rime factor)
!-- t8w           temperature at layer interfaces
!-- cldfra, cldfra_old, current, previous cloud fraction
!-- exch_h        vertical diffusivity (m2/s)
!-- qlsink        Fractional cloud water sink (/s)
!-- precr         rain precipitation rate at all levels (kg/m2/s)
!-- preci         ice precipitation rate at all levels (kg/m2/s)
!-- precs         snow precipitation rate at all levels (kg/m2/s)
!-- precg         graupel precipitation rate at all levels (kg/m2/s)                             &
!-- P_QV          species index for water vapor
!-- P_QC          species index for cloud water
!-- P_QR          species index for rain water
!-- P_QI          species index for cloud ice
!-- P_QS          species index for snow
!-- P_QG          species index for graupel
!-- P_QH          species index for hail
!-- P_QNDROP      species index for cloud drop mixing ratio
!-- P_QNR         species index for rain number concentration,
!-- P_QNI         species index for cloud ice number concentration
!-- P_QNS         species index for snow number concentration,
!-- P_QNG         species index for graupel number concentration,
!-- P_QNH         species index for hail number concentration,
!-- P_QZR         species index for rain    reflectivity
!-- P_QZI         species index for ice     reflectivity
!-- P_QZS         species index for snow    reflectivity
!-- P_QZG         species index for graupel reflectivity
!-- P_QZH         species index for hail    reflectivity
!-- P_QVOLG       species index for graupel particle volume,
!-- id            grid id number
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!-- diagflag      Logical to tell us when to produce diagnostics for history or restart
!
!======================================================================

 integer,intent(in) :: itimestep
 integer :: i,j,k


!  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(INOUT) ::    &
!     qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr,qh_curr




     micro_select: SELECT CASE(mp_physics)

        CASE (GMPUNIFIED)
print*, 'MSG microphysics_driver: drive goddard unified microphysics scheme'
!             IF ( PRESENT( QV_CURR ) .AND. PRESENT ( QC_CURR ) .AND.  &
!                  PRESENT( QR_CURR ) .AND. PRESENT ( QI_CURR ) .AND.  &
!                  PRESENT( QS_CURR ) .AND. PRESENT ( QG_CURR ) .AND.  &
!                  PRESENT( QH_CURR )                           .AND.  &
!                  PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
!                  PRESENT( HAIL    ) .AND. PRESENT ( ICE2    ) .AND.  &
!                  PRESENT( Z       ) .AND. PRESENT ( W       )  ) THEN
               CALL goddard_microphysics_wrf(                       &
                  TH=th                                             &
                 ,QV=qv_curr                                        &
                 ,QC=qc_curr                                        &
                 ,QR=qr_curr                                        &
                 ,QI=qi_curr                                        &
                 ,QS=qs_curr                                        &
                 ,QG=qg_curr                                        &
                 ,QH=qh_curr                                        &
                 ,RHO=rho, PII=pi_phy, P=p, DT_IN=dt, Z=z           &
                 ,HT=ht, DZ8W=dz8w, W=w                             &
                 ,ITIMESTEP=itimestep, XLAND=xland, DX=dx           &
                 ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
                 ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
                 ,RAINNC=rainnc, RAINNCV=rainncv                    &
                 ,SNOWNC=snownc, SNOWNCV=snowncv ,SR=sr             &
                 ,GRAUPELNC=graupelnc ,GRAUPELNCV=graupelncv        &
                 ,HAILNC=hailnc, HAILNCV=hailncv                    &
                 ,REFL_10CM=refl_10cm                               &  ! added for radar reflectivity
                 ,diagflag=diagflag                                 &  ! added for radar reflectivity
                 ,do_radar_ref=do_radar_ref                         &  ! added for radar reflectivity
                 ,PHYSC=physc, PHYSE=physe, PHYSD=physd             &
                 ,PHYSS=physs, PHYSM=physm, PHYSF=physf             &
                 ,ACPHYSC=acphysc, ACPHYSE=acphyse, ACPHYSD=acphysd  &
                 ,ACPHYSS=acphyss, ACPHYSM=acphysm, ACPHYSF=acphysf  &
                 ,REFC=refc, REFR=refr, REFI=refi                   &
                 ,REFS=refs, REFG=refg, REFH=refh                   &
                 ,PRECR3D=precr3d,PRECI3D=preci3d,PRECS3D=precs3d   &
                 ,PRECG3D=precg3d,PRECH3D=prech3d                   &
#if ( WRF_CHEM == 1)
                 ,AERO=aero                                         &
                 ,ICN_DIAG=icn_diag, NC_DIAG=nc_diag, GID=ID        &
                 ,CHEM_OPT=chem_opt                                 &
                 ,GSFCGCE_GOCART_COUPLING=gsfcgce_gocart_coupling   &
#endif
                  )

!               call calc_comdbz(comdbz,refl_10cm,diagflag,do_radar_ref &
!                    ,ids,ide,jds,jde,kds,kde &
!                    ,ims,ime,jms,jme,kms,kme &
!                    ,its,ite,jts,jte,kts,kte)
               do j=jts,jte
                  do k=kts,kte
                     do i=its,ite
                        ! Snapshot total latent heating rate [K/s]
                        phys_tot(i,k,j) = physc(i,k,j) + physe(i,k,j) + &
                             physd(i,k,j) + physs(i,k,j) + physm(i,k,j) + &
                             physf(i,k,j)
                        ! Accumulated total latent heating [K]
                        acphys_tot(i,k,j) = acphysc(i,k,j) + acphyse(i,k,j) + &
                             acphysd(i,k,j) + acphyss(i,k,j) + acphysm(i,k,j) + &
                             acphysf(i,k,j)

                     enddo
                  enddo
               enddo

!             ELSE
!                CALL wrf_error_fatal ( 'arguments not present for calling NUWRF4ICESFCGCE' )
!             ENDIF


        CASE (NUWRF4ICESCHEME)

print*, 'MSG microphysics_driver: drive gsfcgce_4ice_nuwrf scheme'
!             IF ( PRESENT( QV_CURR ) .AND. PRESENT ( QC_CURR ) .AND.  &
!                  PRESENT( QR_CURR ) .AND. PRESENT ( QI_CURR ) .AND.  &
!                  PRESENT( QS_CURR ) .AND. PRESENT ( QG_CURR ) .AND.  &
!                  PRESENT( QH_CURR )                           .AND.  &
!                  PRESENT( RAINNC  ) .AND. PRESENT ( RAINNCV ) .AND.  &
!                  PRESENT( HAIL    ) .AND. PRESENT ( ICE2    ) .AND.  &
!                  PRESENT( Z       ) .AND. PRESENT ( W       )  ) THEN
               CALL gsfcgce_4ice_nuwrf(                             &
                  TH=th                                             &
                 ,QV=qv_curr                                        &
                 ,QL=qc_curr                                        &
                 ,QR=qr_curr                                        &
                 ,QI=qi_curr                                        &
                 ,QS=qs_curr                                        &
                 ,QG=qg_curr                                        &
                 ,QH=qh_curr                                        &
                 ,RHO=rho, PII=pi_phy, P=p, DT_IN=dt, Z=z           &
                 ,HT=ht, DZ8W=dz8w, GRAV=g, W=w                     &
                 ,RHOWATER=rhowater, RHOSNOW=rhosnow                &
                 ,ITIMESTEP=itimestep, XLAND=xland, DX=dx           &
                 ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
                 ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
                 ,RAINNC=rainnc, RAINNCV=rainncv                    &
                 ,SNOWNC=snownc, SNOWNCV=snowncv ,SR=sr             &
                 ,GRAUPELNC=graupelnc ,GRAUPELNCV=graupelncv        &
                 ,HAILNC=hailnc, HAILNCV=hailncv                    &
                 ,REFL_10CM=refl_10cm                               &  ! added for radar reflectivity
                 ,diagflag=diagflag                                 &  ! added for radar reflectivity
                 ,do_radar_ref=do_radar_ref                         &  ! added for radar reflectivity
!NUWRF JJS 20110525 vvvvv
                 ,PHYSC=physc, PHYSE=physe, PHYSD=physd             &
                 ,PHYSS=physs, PHYSM=physm, PHYSF=physf             &
                 ,ACPHYSC=acphysc, ACPHYSE=acphyse, ACPHYSD=acphysd  &
                 ,ACPHYSS=acphyss, ACPHYSM=acphysm, ACPHYSF=acphysf  &
                 ,REFC=refc, REFR=refr, REFI=refi                   &
                 ,REFS=refs, REFG=refg, REFH=refh                   &
                 ,PRECR3D=precr3d,PRECI3D=preci3d,PRECS3D=precs3d   &
                 ,PRECG3D=precg3d,PRECH3D=prech3d                   &
#if ( WRF_CHEM == 1)
                 ,AERO=aero                                         &
                 ,ICN_DIAG=icn_diag, NC_DIAG=nc_diag, GID=ID        &
                 ,CHEM_OPT=chem_opt                                 &
                 ,GSFCGCE_GOCART_COUPLING=gsfcgce_gocart_coupling   &
#endif
!NUWRF JJS 20110525 ^^^^^
                                                                    )
               !EMK NUWRF
!               call calc_comdbz(comdbz,refl_10cm,diagflag,do_radar_ref &
!                    ,ids,ide,jds,jde,kds,kde &
!                    ,ims,ime,jms,jme,kms,kme &
!                    ,its,ite,jts,jte,kts,kte)
               !EMK NUWRF

               do j=jts,jte
                  do k=kts,kte
                     do i=its,ite
                        ! Snapshot total latent heating rate [K/s]
                        phys_tot(i,k,j) = physc(i,k,j) + physe(i,k,j) + &
                             physd(i,k,j) + physs(i,k,j) + physm(i,k,j) + &
                             physf(i,k,j)
                        ! Accumulated total latent heating [K]
                        acphys_tot(i,k,j) = acphysc(i,k,j) + acphyse(i,k,j) + &
                             acphysd(i,k,j) + acphyss(i,k,j) + acphysm(i,k,j) + &
                             acphysf(i,k,j)

                     end do
                  end do
               end do

!             ELSE
!                CALL wrf_error_fatal ( 'arguments not present for calling NUWRF4ICESFCGCE' )
!             ENDIF


      CASE DEFAULT

!         WRITE( wrf_err_message , * ) 'The microphysics option does not exist: mp_physics = ', mp_physics
!         CALL wrf_error_fatal ( wrf_err_message )

      END SELECT micro_select

!   CALL wrf_debug ( 200 , 'microphysics_driver: returning from' )

   RETURN

 END SUBROUTINE microphysics_driver

END  MODULE module_microphysics_driver

