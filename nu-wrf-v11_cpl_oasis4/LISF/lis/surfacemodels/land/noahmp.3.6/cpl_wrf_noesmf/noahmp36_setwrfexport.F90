!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: noahmp36_setwrfexport.F90
!
! !DESCRIPTION:  
!  Defines the export states from Noah to WRF in coupled mode
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 17 Nov 2008; Sujay Kumar, Modified for the ESMF coupled version
! 11 Dec 2015; Eric Kemp, updated for Noah 3.6.
! 08 Nov 2018; Carlos Cruz, Adapted for NoahMP 3.6.
! 
! !INTERFACE:
subroutine noahmp36_setwrfexport(n)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_logMod
  use LISWRFGridCompMod, only : LISWRF_export
  use noahmp36_lsmMod
  
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
  a_array = noahmp36_struc(n)%noahmp36%trad
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%avgsurft_t, a_array)
  
  ! fsh -> qh_t -> HFX
  ! total sensible heat to atmosphere
  a_array = noahmp36_struc(n)%noahmp36%fsh
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%qh_t, a_array)

  ! qfx -> eta_kinematic_t -> QFX
  ! upward moisture flux
  a_array = noahmp36_struc(n)%noahmp36%qfx
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%eta_kinematic_t, a_array)

  ! lh -> qle_t -> LH
  ! total latent heat flux
  a_array = noahmp36_struc(n)%noahmp36%lh
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%qle_t, a_array)
     
  ! ssoil -> qg_t -> GRDFLX
  ! ground/snow heat flux to soil
  a_array = noahmp36_struc(n)%noahmp36%ssoil
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%qg_t, a_array)

  ! albedo -> albedo_t -> ALBEDO
  ! surface albedo
  a_array = noahmp36_struc(n)%noahmp36%albedo
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%albedo_t, a_array)

  ! z0 -> znt_t -> ZNT
  ! roughness length 
  a_array = noahmp36_struc(n)%noahmp36%z0
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%znt_t, a_array)

  ! qsfc -> q1_t -> QSFC
  ! mixing ratio
  a_array = noahmp36_struc(n)%noahmp36%qsfc
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%q1_t, a_array)

  ! ch -> chs2_t -> CHS2
  ! Sensible heat exchange coefficient
  a_array = noahmp36_struc(n)%noahmp36%ch
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%chs2_t, a_array)
     
  ! cm -> cqs2_t -> CQS2
  ! momemtum drag coefficient
  a_array = noahmp36_struc(n)%noahmp36%cm
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%cqs2_t, a_array)
     
  ! fsno -> snocvr_t -> SNOWC
  ! snow cover fraction
  a_array = noahmp36_struc(n)%noahmp36%fsno
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%snocvr_t, a_array)

  ! sneqv -> snow_t -> SNOW
  ! snow water equivalent
  a_array = noahmp36_struc(n)%noahmp36%sneqv
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%snow_t, a_array)

  ! snowh -> snowh_t -> SNOWH
  ! snow height
  a_array = noahmp36_struc(n)%noahmp36%snowh
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%snowh_t, a_array)

  ! smcmax -> lispor_t -> LISPOROSITY
  ! porosity, saturated value of soil moisture (volumetric)
  a_array = noahmp36_struc(n)%noahmp36%smcmax
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%lispor_t, a_array)

  ! canwat=canliq+canice -> cmc_t -> CANWAT 
  ! intercepted liquid water (mm)
  a_array = noahmp36_struc(n)%noahmp36%canliq + noahmp36_struc(n)%noahmp36%canice
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%cmc_t, a_array)

  ! emissi -> emissi -> EMISS
  ! surface emissivity
  a_array = noahmp36_struc(n)%noahmp36%emissi
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%emiss_t, a_array)

  ! fveg -> fveg_t -> FVEGXY
  ! vegetation fraction
  a_array = noahmp36_struc(n)%noahmp36%fveg
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%fveg_t, a_array)

  ! t2mb -> t2mb_t -> T2MBXY
  ! 2m air temperature - bare ground
  a_array = noahmp36_struc(n)%noahmp36%t2mb
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%t2mb_t, a_array)

  ! t2mv -> t2mv_t -> T2MVXY
  ! 2m air temperature - over vegetation
  a_array = noahmp36_struc(n)%noahmp36%t2mv
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%t2mv_t, a_array)

  ! q2b -> q2mb_t -> Q2MBXY
  ! 2m air spec. hum. - bare ground
  ! convert to mixing ratio (module_sf_noahmpdrv)
  a_array = noahmp36_struc(n)%noahmp36%q2b/(1.0-noahmp36_struc(n)%noahmp36%q2b)
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%q2mb_t, a_array)

  ! q2v -> q2mv_t -> Q2MVXY
  ! 2m air spec. hum. - over vegetation
  ! convert to mixing ratio (module_sf_noahmpdrv)
  a_array = noahmp36_struc(n)%noahmp36%q2v/(1.0-noahmp36_struc(n)%noahmp36%q2v)
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%q2mv_t, a_array)

  ! smc -> smc_t -> SMOIS
  ! volumetric soil moisture, ice + liquid
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp36_struc(n)%noahmp36(t)%smc(1)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%smc1_t, a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp36_struc(n)%noahmp36(t)%smc(2)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%smc2_t, a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp36_struc(n)%noahmp36(t)%smc(3)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%smc3_t, a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp36_struc(n)%noahmp36(t)%smc(4)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%smc4_t, a_array)
     
  ! sstc -> stc_t -> TSLB
  ! soil layer temperature
  block
    integer :: i0 ! first index of soil layer temp
    i0 = NOAHMP36_struc(n)%nsnow+1 
    do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
       a_array(t) = NOAHMP36_struc(n)%noahmp36(t)%sstc(i0) !soil_temp(1)
    end do
    call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%stc1_t, a_array)
    do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
       a_array(t) = NOAHMP36_struc(n)%noahmp36(t)%sstc(i0+1) !soil_temp(2)
    end do
    call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%stc2_t, a_array)
    do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
       a_array(t) = NOAHMP36_struc(n)%noahmp36(t)%sstc(i0+2) !soil_temp(3)
    end do
    call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%stc3_t, a_array)
    do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
       a_array(t) = NOAHMP36_struc(n)%noahmp36(t)%sstc(i0+3) !soil_temp(4)
    end do
    call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%stc4_t, a_array)
  end block
  
  ! sh2o -> sh2o_t -> SH2O
  ! volumetric liquid soil moisture 
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp36_struc(n)%noahmp36(t)%sh2o(1)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%sh2o1_t,a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp36_struc(n)%noahmp36(t)%sh2o(2)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%sh2o2_t,a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp36_struc(n)%noahmp36(t)%sh2o(3)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%sh2o3_t,a_array)
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     a_array(t) = noahmp36_struc(n)%noahmp36(t)%sh2o(4)
  end do
  call LIS_patch2tile(n,LIS_rc%lsm_index, LISWRF_export(n)%sh2o4_t,a_array)
  
  deallocate(a_array)

end subroutine noahmp36_setwrfexport
 
