!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: noahmp401_setwrfexport.F90
!
! !DESCRIPTION:  
!  Defines the export states from Noah to WRF in coupled mode
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 17 Nov 2008; Sujay Kumar, Modified for the ESMF coupled version
! 11 Dec 2015; Eric Kemp, updated for Noah 3.6.
! 05 Feb 2020; Carlos Cruz, Adapted for NoahMP 4.0.1.
! 
! !INTERFACE:
subroutine noahmp401_setwrfexport(n)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_logMod
  use LISWRFGridCompMod, only : LISWRF_export
  use noahmp401_lsmMod
  
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n ! nest index
!EOP
  integer               :: t ! tile index
  real, allocatable     :: a_array(:)
  real                  :: q2
  
  allocate(a_array(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  ! trad -> avgsurft_t -> TSK
  ! average surface temperature
  a_array = noahmp401_struc(n)%noahmp401%trad
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%avgsurft_t, a_array)
  
  ! hfx -> qh_t -> HFX
  ! total sensible heat to atmosphere
  a_array = noahmp401_struc(n)%noahmp401%hfx
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%qh_t, a_array)

  ! qfx -> eta_kinematic_t -> QFX
  ! upward moisture flux
  a_array = noahmp401_struc(n)%noahmp401%qfx
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%eta_kinematic_t, a_array)

  ! lh -> qle_t -> LH
  ! total latent heat flux
  a_array = noahmp401_struc(n)%noahmp401%lh
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%qle_t, a_array)
     
  ! grdflx -> qg_t -> GRDFLX
  ! ground/snow heat flux to soil
  a_array = noahmp401_struc(n)%noahmp401%grdflx
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%qg_t, a_array)

  ! albedo -> albedo_t -> ALBEDO
  ! surface albedo
  a_array = noahmp401_struc(n)%noahmp401%albedo
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%albedo_t, a_array)

  ! z0 -> znt_t -> ZNT
  ! roughness length 
  a_array = noahmp401_struc(n)%noahmp401%z0
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%znt_t, a_array)

  ! qsfc -> q1_t -> QSFC
  ! mixing ratio   ! calculate based on qair (spec. hum) as: mr = qair / (1-qair)
  a_array = noahmp401_struc(n)%noahmp401%qair/(1.0-noahmp401_struc(n)%noahmp401%qair)
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%q1_t, a_array)

  ! ch -> chs2_t -> CHS2
  ! Sensible heat exchange coefficient
  a_array = noahmp401_struc(n)%noahmp401%ch
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%chs2_t, a_array)
     
  ! cm -> cqs2_t -> CQS2
  ! momemtum drag coefficient
  a_array = noahmp401_struc(n)%noahmp401%cm
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%cqs2_t, a_array)
     
  ! snowc -> snocvr_t -> SNOWC
  ! snow cover fraction
  a_array = noahmp401_struc(n)%noahmp401%snowc
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%snocvr_t, a_array)

  ! sneqv -> snow_t -> SNOW
  ! snow water equivalent
  a_array = noahmp401_struc(n)%noahmp401%sneqv
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%snow_t, a_array)

  ! snowh -> snowh_t -> SNOWH
  ! snow height
  a_array = noahmp401_struc(n)%noahmp401%snowh
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%snowh_t, a_array)

  ! smcmax -> lispor_t -> LISPOROSITY
  ! porosity, saturated value of soil moisture (volumetric)
  a_array = 0.01 ! need to get smcmax from land model
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%lispor_t, a_array)

  ! canwat -> cmc_t -> CANWAT 
  ! intercepted liquid water (mm)
  a_array = noahmp401_struc(n)%noahmp401%canwat
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%cmc_t, a_array)

  ! emiss -> emissi -> EMISS
  ! surface emissivity
  a_array = noahmp401_struc(n)%noahmp401%emiss
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%emiss_t, a_array)


  ! fveg -> fveg_t -> FVEGXY
  ! vegetation fraction
  a_array = noahmp401_struc(n)%noahmp401%fveg
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%fveg_t, a_array)

  ! t2mb -> t2mb_t -> T2MBXY
  ! 2m air temperature - bare ground
  a_array = noahmp401_struc(n)%noahmp401%t2mb
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%t2mb_t, a_array)

  ! t2mv -> t2mv_t -> T2MVXY
  ! 2m air temperature - over vegetation
  a_array = noahmp401_struc(n)%noahmp401%t2mv
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%t2mv_t, a_array)

  ! q2b -> q2mb_t -> Q2MBXY
  ! 2m air spec. hum. - bare ground
  a_array = noahmp401_struc(n)%noahmp401%q2mb
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%q2mb_t, a_array)

  ! q2v -> q2mv_t -> Q2MVXY
  ! 2m air spec. hum. - over vegetation
  a_array = noahmp401_struc(n)%noahmp401%q2mv
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%q2mv_t, a_array)

  ! smc -> smc_t -> SMOIS
  ! volumetric soil moisture, ice + liquid
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%smc1_t, a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp401_struc(n)%noahmp401(t)%smc(2)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%smc2_t, a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp401_struc(n)%noahmp401(t)%smc(3)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%smc3_t, a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp401_struc(n)%noahmp401(t)%smc(4)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%smc4_t, a_array)
     
  ! tslb -> stc_t -> TSLB
  ! soil layer temperature
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = NOAHMP401_struc(n)%noahmp401(t)%tslb(1) !soil_temp(1)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%stc1_t, a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = NOAHMP401_struc(n)%noahmp401(t)%tslb(2) !soil_temp(2)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%stc2_t, a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = NOAHMP401_struc(n)%noahmp401(t)%tslb(3) !soil_temp(3)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%stc3_t, a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = NOAHMP401_struc(n)%noahmp401(t)%tslb(4) !soil_temp(4)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%stc4_t, a_array)
  
  ! sh2o -> sh2o_t -> SH2O
  ! volumetric liquid soil moisture 
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp401_struc(n)%noahmp401(t)%sh2o(1)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%sh2o1_t,a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp401_struc(n)%noahmp401(t)%sh2o(2)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%sh2o2_t,a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp401_struc(n)%noahmp401(t)%sh2o(3)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%sh2o3_t,a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp401_struc(n)%noahmp401(t)%sh2o(4)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%sh2o4_t,a_array)
  
  deallocate(a_array)

end subroutine noahmp401_setwrfexport
 
