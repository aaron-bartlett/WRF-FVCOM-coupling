!###############################################################################!
!  CCATT-BRAMS/MCGA-CPTEC/WRF emission model    CPTEC/INPE                      !
!  Version 1.0.0: 12/nov/2010                                                   !
!  Coded by Saulo Freitas and Karla Longo                                       !
!  Contact: gmai@cptec.inpe.br - http://meioambiente.cptec.inpe.br              !
!###############################################################################!

subroutine convert_AeM_to_wrf_reac(isp,iespc,ident,source_type,spc_name)
!use chem1_list
use grid_dims_out, only:  use_bbem, use_gfedv2, use_gfedv3, use_gfedv4, use_fwbawb, use_qfed ! NUWRF
use emiss_vars_emissions
use bbbem_emissions , only:  bbbem_g                     
use gfedv2_emissions, only:  gfedv2_g                     
use gfedv3_emissions, only:  gfedv3_g ! NUWRF
use gfedv4_emissions, only:  gfedv4_g ! NUWRF
use qfed_emissions, only: qfed_g ! NUWRF
use fwbawb_emissions, only : fwbawb_g

use AeM_emission_factors 


!implicit none  
integer, intent(in) :: isp
integer, intent(inout) :: iespc,ident
character (len=*) source_type
character (len=*), intent(in)  :: spc_name  !kml 


!--     AeM    |    RACM
!--     Emissions for RACM weighted by Agregation factor Agg (only organics)
!--     Agregation factors for RACM can be found in Stockwell et al 1997, 
!-- 	J. Geos Res. 102, 25847-25879


!- NOx	=>   NO
!      change MP 201108 : NO=NOX/2,   NO2 = NOX/2 in molecules (40% and 60% in mass)
if(spc_name == 'NO') then
   ident = NOx
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.4* bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.4*gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.4*gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.4*gfedv4_g(ident)%src(:,:,1) ! NUWRF
   ! Special QFED rule: Split NOx into 90% NO and 10% NO2
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.9*qfed_g(ident)%src(:,:,1) ! NUWRF
   
   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1 ) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
      
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + 0.4*fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1


   print*,'==> converted from AeM - found for ',spc_name

   return
endif
if(spc_name == 'NO2') then
   ident = NOx
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.6* bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.6*gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.6*gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.6*gfedv4_g(ident)%src(:,:,1) ! NUWRF
   ! Special QFED rule: Split NOx into 90% NO and 10% NO2
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.1*qfed_g(ident)%src(:,:,1) ! NUWRF
   
   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1 ) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
      
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + 0.6*fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1


   print*,'==> converted from AeM - found for ',spc_name

   return
endif


!----------------------------------------------------------------------------
!- C2H6 =>	ETH (Agg=1)
if(spc_name == 'ETH') then
   ident = C2H6
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!- C2H4	=> ETE (Agg=1)
if(spc_name == 'ETE') then
   ident = C2H4
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1 ) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!- Xylene   =>     XYL (Agg=1)
if(spc_name == 'XYL') then
   ident = xylenes
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1 ) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
! - CSL (m=108)(Agg=1)
! Phenol (m=94)(Agg=1)
!RACM)
if(spc_name == 'CSL') then
   ident = Phenol
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   = bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1 ) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!Formaldehyde	      =>  HCHO (Agg=1)
if(spc_name == 'HCHO') then
   ident = Formaldehyde
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!HFo	=>       ORA1 (Agg=1)
if(spc_name == 'ORA1') then
   ident = HFo
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif



!----------------------------------------------------------------------------
!Isoprene  =>	ISO (Agg=1)
if(spc_name == 'ISO') then
   ident = Isoprene
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1 ) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!---------------------------------------------------------------------------
!-ORA2 (m=60, Agg=1)
! HAc (m=60, Agg=1)
! Propanoic (m=74, Agg=1)
! HAc	     =>   ORA2
! Propanoic  =>   ORA2
if(spc_name == 'ORA2') then
   ident = HAc
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)

   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(HAc	     )%src(:,:,1)  + &
                                                            bbbem_g(Propanoic)%src(:,:,1)     
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  gfedv2_g(HAc      )%src(:,:,1) + &
   							    gfedv2_g(Propanoic)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  gfedv3_g(HAc      )%src(:,:,1) + &
   							    gfedv3_g(Propanoic)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  gfedv4_g(HAc      )%src(:,:,1) + &
   							    gfedv4_g(Propanoic)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  qfed_g(HAc      )%src(:,:,1) + &
   							    qfed_g(Propanoic)%src(:,:,1) ! NUWRF
							    
   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
                                                                                            fwbawb_g(HAc      )%src(:,:,1) + &
   							                                   fwbawb_g(Propanoic)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!-TOL (m=92)
! benzene (m=82, Agg=0.29)
! toluene (m=92, Agg=1)
! ethylbenzene (m=106, Agg=1)
! styrene (m=104, Agg=1)
! PAH (m=164, Agg=1)
!benzene	  =>      TOL
!toluene	  =>      TOL
!ethylbenzene	  =>      TOL
!styrene	  =>      TOL
!PAH	          =>      TOL
!    change MP 201108    Toluene only => TOL
if(spc_name == 'TOL') then
   ident = toluene
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)

   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.29*bbbem_g(benzene	)%src(:,:,1) + &
                                                           bbbem_g(toluene	)%src(:,:,1) + &
                                                           bbbem_g(ethylbenzene)%src(:,:,1) + &
                                                           bbbem_g(styrene	)%src(:,:,1) + &
                                                           bbbem_g(PAH 	)%src(:,:,1)    
      
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.29*gfedv2_g(benzene	 )%src(:,:,1) + &
   							    gfedv2_g(toluene	 )%src(:,:,1) + &
   							    gfedv2_g(ethylbenzene)%src(:,:,1) + &
   							    gfedv2_g(styrene	 )%src(:,:,1) + &
   							    gfedv2_g(PAH 	 )%src(:,:,1)	
							    
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.29*gfedv3_g(benzene	 )%src(:,:,1) + &
   							    gfedv3_g(toluene	 )%src(:,:,1) + &
   							    gfedv3_g(ethylbenzene)%src(:,:,1) + &
   							    gfedv3_g(styrene	 )%src(:,:,1) + &
   							    gfedv3_g(PAH 	 )%src(:,:,1) ! NUWRF
	
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.29*gfedv4_g(benzene	 )%src(:,:,1) + &
   							    gfedv4_g(toluene	 )%src(:,:,1) + &
   							    gfedv4_g(ethylbenzene)%src(:,:,1) + &
   							    gfedv4_g(styrene	 )%src(:,:,1) + &
   							    gfedv4_g(PAH 	 )%src(:,:,1) ! NUWRF
	
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.29*qfed_g(benzene	 )%src(:,:,1) + &
   							    qfed_g(toluene	 )%src(:,:,1) + &
   							    qfed_g(ethylbenzene)%src(:,:,1) + &
   							    qfed_g(styrene	 )%src(:,:,1) + &
   							    qfed_g(PAH 	 )%src(:,:,1) ! NUWRF
						    
   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1 ) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
                                                          0.29*fwbawb_g(benzene	)%src(:,:,1) + &
   							  fwbawb_g(toluene	)%src(:,:,1) + &
   							  fwbawb_g(ethylbenzene )%src(:,:,1) + &
   							  fwbawb_g(styrene	)%src(:,:,1) + &
   							  fwbawb_g(PAH  	)%src(:,:,1)   

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!API        =>	0.5*terpenes (Agg=1)
if(spc_name == 'API') then
   ident = TERPENES
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.5* bbbem_g(TERPENES)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.5*gfedv2_g(TERPENES)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.5*gfedv3_g(TERPENES)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.5*gfedv4_g(TERPENES)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.5*qfed_g(TERPENES)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                     + 0.5* fwbawb_g(TERPENES)%src(:,:,1)
   
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!LIM        =>	0.2*terpenes (Agg=1)
if(spc_name == 'LIM') then
   ident = TERPENES
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.2* bbbem_g(TERPENES)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.2*gfedv2_g(TERPENES)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.2*gfedv3_g(TERPENES)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.2*gfedv4_g(TERPENES)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  0.2*qfed_g(TERPENES)%src(:,:,1) ! NUWRF
							    
   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1 ) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                     + 0.2* fwbawb_g(TERPENES)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!-OLT (m=42)  (Agg=1)
! C3H6 (m=42)
! 1_butene (m=56)
! 1_pentene (m=70)
! 2_Me_Butene (m=70)
! 4_me_1_pentene (m=84)
! 2_me_1_pentene (m=84)
! 1_hexene (m=84)
! octenes (m=112)
!-C3H6	          =>	OLT
!-1_butene	  =>	OLT
!-1_pentene	  =>	OLT
!-2_Me_Butene	  =>	OLT
!-4_me_1_pentene  =>	OLT
!-2_me_1_pentene  =>	OLT
!-1_hexene	  =>    OLT
!-octenes	  =>    OLT
!     change MP 201108 C3H6 only => OLT
if(spc_name == 'OLT') then
   ident = C3H6
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = bbbem_g(C3H6	      )%src(:,:,1) + &
 							bbbem_g(butene_1      )%src(:,:,1) + &
							bbbem_g(pentene_1     )%src(:,:,1) + &
  							bbbem_g(Butene_2_Me   )%src(:,:,1) + &
  							bbbem_g(pentene_4_me_1)%src(:,:,1) + &
  							bbbem_g(pentene_2_me_1)%src(:,:,1) + &
 							bbbem_g(hexene_1      )%src(:,:,1) + & 
							bbbem_g(octenes       )%src(:,:,1)    	       

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv2_g(C3H6	       )%src(:,:,1) + &
  							gfedv2_g(butene_1      )%src(:,:,1) + &
  							gfedv2_g(pentene_1     )%src(:,:,1) + &
   							gfedv2_g(Butene_2_Me   )%src(:,:,1) + &
   							gfedv2_g(pentene_4_me_1)%src(:,:,1) + &
  							gfedv2_g(pentene_2_me_1)%src(:,:,1) + &
 							gfedv2_g(hexene_1      )%src(:,:,1) + & 
							gfedv2_g(octenes       )%src(:,:,1)    	       

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv3_g(C3H6	       )%src(:,:,1) + &
  							gfedv3_g(butene_1      )%src(:,:,1) + &
  							gfedv3_g(pentene_1     )%src(:,:,1) + &
   							gfedv3_g(Butene_2_Me   )%src(:,:,1) + &
   							gfedv3_g(pentene_4_me_1)%src(:,:,1) + &
  							gfedv3_g(pentene_2_me_1)%src(:,:,1) + &
 							gfedv3_g(hexene_1      )%src(:,:,1) + & 
							gfedv3_g(octenes       )%src(:,:,1) ! NUWRF
   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv4_g(C3H6	       )%src(:,:,1) + &
  							gfedv4_g(butene_1      )%src(:,:,1) + &
  							gfedv4_g(pentene_1     )%src(:,:,1) + &
   							gfedv4_g(Butene_2_Me   )%src(:,:,1) + &
   							gfedv4_g(pentene_4_me_1)%src(:,:,1) + &
  							gfedv4_g(pentene_2_me_1)%src(:,:,1) + &
 							gfedv4_g(hexene_1      )%src(:,:,1) + & 
							gfedv4_g(octenes       )%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1) = qfed_g(C3H6	       )%src(:,:,1) + &
  							qfed_g(butene_1      )%src(:,:,1) + &
  							qfed_g(pentene_1     )%src(:,:,1) + &
   							qfed_g(Butene_2_Me   )%src(:,:,1) + &
   							qfed_g(pentene_4_me_1)%src(:,:,1) + &
  							qfed_g(pentene_2_me_1)%src(:,:,1) + &
 							qfed_g(hexene_1      )%src(:,:,1) + & 
							qfed_g(octenes       )%src(:,:,1) ! NUWRF


   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
                                                          fwbawb_g(C3H6 	 )%src(:,:,1) + &
  						          fwbawb_g(butene_1	 )%src(:,:,1) + &
  						          fwbawb_g(pentene_1	 )%src(:,:,1) + &
    						          fwbawb_g(Butene_2_Me   )%src(:,:,1) + &
   						          fwbawb_g(pentene_4_me_1)%src(:,:,1) + &
  						          fwbawb_g(pentene_2_me_1)%src(:,:,1) + &
 						          fwbawb_g(hexene_1	 )%src(:,:,1) + & 
						          fwbawb_g(octenes	 )%src(:,:,1)	 

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!-KET (m=72)
!-Acetone (m=58)	  =>KET (Agg=0.33)
!-2_Butanone (m=72)       =>KET (Agg=1.61)
!-2_3_Butanedione (m=86)  =>KET (Agg=1.61)
!-Pentanones (m=86)	  =>KET (Agg=1.61)
!-Hexanones (m=100)	  =>KET (Agg=1.61)
!-Heptanones (m=114)	  =>KET (Agg=1.61)
!-Octanones (m=128	  =>KET (Agg=1.61)
if(spc_name == 'KET') then
   ident = Acetone
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.33*bbbem_g(Acetone        )%src(:,:,1) + &
  							1.61*bbbem_g(Butanone_2     )%src(:,:,1) + &
  							1.61*bbbem_g(Butanedione_2_3)%src(:,:,1) + &
  							1.61*bbbem_g(Pentanones     )%src(:,:,1) + &
   							1.61*bbbem_g(Hexanones      )%src(:,:,1) + &
   							1.61*bbbem_g(Heptanones     )%src(:,:,1) + &
  							1.61*bbbem_g(Octanones      )%src(:,:,1)		

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.33*gfedv2_g(Acetone	)%src(:,:,1) + &
  							1.61*gfedv2_g(Butanone_2	)%src(:,:,1) + &
  							1.61*gfedv2_g(Butanedione_2_3)%src(:,:,1) + &
  							1.61*gfedv2_g(Pentanones	)%src(:,:,1) + &
   							1.61*gfedv2_g(Hexanones	)%src(:,:,1) + &
   							1.61*gfedv2_g(Heptanones	)%src(:,:,1) + &
  							1.61*gfedv2_g(Octanones	)%src(:,:,1)		

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.33*gfedv3_g(Acetone	)%src(:,:,1) + &
  							1.61*gfedv3_g(Butanone_2	)%src(:,:,1) + &
  							1.61*gfedv3_g(Butanedione_2_3)%src(:,:,1) + &
  							1.61*gfedv3_g(Pentanones	)%src(:,:,1) + &
   							1.61*gfedv3_g(Hexanones	)%src(:,:,1) + &
   							1.61*gfedv3_g(Heptanones	)%src(:,:,1) + &
  							1.61*gfedv3_g(Octanones	)%src(:,:,1) ! NUWRF
   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.33*gfedv4_g(Acetone	)%src(:,:,1) + &
  							1.61*gfedv4_g(Butanone_2	)%src(:,:,1) + &
  							1.61*gfedv4_g(Butanedione_2_3)%src(:,:,1) + &
  							1.61*gfedv4_g(Pentanones	)%src(:,:,1) + &
   							1.61*gfedv4_g(Hexanones	)%src(:,:,1) + &
   							1.61*gfedv4_g(Heptanones	)%src(:,:,1) + &
  							1.61*gfedv4_g(Octanones	)%src(:,:,1) ! NUWRF
   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.33*qfed_g(Acetone	)%src(:,:,1) + &
  							1.61*qfed_g(Butanone_2	)%src(:,:,1) + &
  							1.61*qfed_g(Butanedione_2_3)%src(:,:,1) + &
  							1.61*qfed_g(Pentanones	)%src(:,:,1) + &
   							1.61*qfed_g(Hexanones	)%src(:,:,1) + &
   							1.61*qfed_g(Heptanones	)%src(:,:,1) + &
  							1.61*qfed_g(Octanones	)%src(:,:,1) ! NUWRF


   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
                                                          0.33*fwbawb_g(Acetone	  )%src(:,:,1) + &
  							  1.61*fwbawb_g(Butanone_2	  )%src(:,:,1) + &
  							  1.61*fwbawb_g(Butanedione_2_3)%src(:,:,1) + &
  							  1.61*fwbawb_g(Pentanones	  )%src(:,:,1) + &
   							  1.61*fwbawb_g(Hexanones	  )%src(:,:,1) + &
  							  1.61*fwbawb_g(Heptanones	  )%src(:,:,1) + &
  							  1.61*fwbawb_g(Octanones	  )%src(:,:,1)

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif
!----------------------------------------------------------------------------
!-DIEN (m=54) (Agg=1)
!butadiene (m=54)	       => DIEN
!pentadienes (m=68)	       => DIEN
!cyclopentadiene (m=66)       => DIEN
!hexadienes (m=82)	       => DIEN
if(spc_name == 'DIEN') then
   ident = butadiene
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = bbbem_g(butadiene      )%src(:,:,1) + &
  							bbbem_g(pentadienes    )%src(:,:,1) + &
  							bbbem_g(cyclopentadiene)%src(:,:,1) + &
  							bbbem_g(hexadienes     )%src(:,:,1) 		


   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv2_g(butadiene	)%src(:,:,1) + &
  							gfedv2_g(pentadienes	)%src(:,:,1) + &
  							gfedv2_g(cyclopentadiene)%src(:,:,1) + &
  							gfedv2_g(hexadienes	)%src(:,:,1)		

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv3_g(butadiene	)%src(:,:,1) + &
  							gfedv3_g(pentadienes	)%src(:,:,1) + &
  							gfedv3_g(cyclopentadiene)%src(:,:,1) + &
  							gfedv3_g(hexadienes	)%src(:,:,1) ! NUWRF
   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv4_g(butadiene	)%src(:,:,1) + &
  							gfedv4_g(pentadienes	)%src(:,:,1) + &
  							gfedv4_g(cyclopentadiene)%src(:,:,1) + &
  							gfedv4_g(hexadienes	)%src(:,:,1) ! NUWRF
   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = qfed_g(butadiene	)%src(:,:,1) + &
  							qfed_g(pentadienes	)%src(:,:,1) + &
  							qfed_g(cyclopentadiene)%src(:,:,1) + &
  							qfed_g(hexadienes	)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
                                                          fwbawb_g(butadiene	  )%src(:,:,1) + &
  							  fwbawb_g(pentadienes	  )%src(:,:,1) + &
  							  fwbawb_g(cyclopentadiene)%src(:,:,1) + &
  							  fwbawb_g(hexadienes	  )%src(:,:,1)


   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!-HC8 (m=114) 
!heptane (m=100)	 => HC8 (Agg=0.94)
!Butanols (m=74)	 => HC8 (Agg=0.91) Not in the list!
! Calculated using (16) of Stockwell 1997
! k_i from Atkinson and Arey, 2003
! Chemical Reviews, Vol. 103, No. 12
!cyclopentanol (m=86)	 => HC8 (Agg=1.01) Not in the list!
! Calculated using (16) of Stockwell 1997
! k_i from Atkinson and Arey, 2003
! Chemical Reviews, Vol. 103, No. 12

if(spc_name == 'HC8') then
   ident = heptane
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.94*bbbem_g(heptane  )%src(:,:,1) + &
  							0.91*bbbem_g(Butanols  )%src(:,:,1) + &
  							1.01*bbbem_g(cyclopentanol)%src(:,:,1)	   

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.94*gfedv2_g(heptane      )%src(:,:,1) + &
  							0.91*gfedv2_g(Butanols     )%src(:,:,1) + &
  							1.01*gfedv2_g(cyclopentanol)%src(:,:,1)	    

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.94*gfedv3_g(heptane      )%src(:,:,1) + &
  							0.91*gfedv3_g(Butanols     )%src(:,:,1) + &
  							1.01*gfedv3_g(cyclopentanol)%src(:,:,1)	! NUWRF

   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.94*gfedv4_g(heptane      )%src(:,:,1) + &
  							0.91*gfedv4_g(Butanols     )%src(:,:,1) + &
  							1.01*gfedv4_g(cyclopentanol)%src(:,:,1)	! NUWRF

   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.94*qfed_g(heptane      )%src(:,:,1) + &
  							0.91*qfed_g(Butanols     )%src(:,:,1) + &
  							1.01*qfed_g(cyclopentanol)%src(:,:,1)	! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
 							  0.94*fwbawb_g(heptane      )%src(:,:,1) + &
  							  0.91*fwbawb_g(Butanols     )%src(:,:,1) + &
  							  1.01*fwbawb_g(cyclopentanol)%src(:,:,1)	    


   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif
!----------------------------------------------------------------------------
!-OLI (m=68)
!tr_2_butene (m=56) (butene_tr_2)  => OLI  (Agg=1) (internal alkene)
!cis_2_butene (m=56)	 (butene_cis_2) => OLI (Agg=1) (internal alkene)
!-2_pentene (m=70)	                => OLI (Agg=1) (internal alkene)
!cyclopentene (m=68)	                => OLI (Agg=0.5) (mix alkene, cyclic)
!i-butene (butene_i)(m=56)	        => OLT or OLI (Agg=0.5) (terminal (not internal) alkene)
if(spc_name == 'OLI') then
   ident = butene_tr_2
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = bbbem_g(butene_tr_2 )%src(:,:,1) + &
  							bbbem_g(butene_cis_2)%src(:,:,1) + &
  							bbbem_g(pentene_2   )%src(:,:,1) + &
  							0.5*bbbem_g(cyclopentene)%src(:,:,1) + &
  							0.5*bbbem_g(butene_i    )%src(:,:,1)	  

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv2_g(butene_tr_2 )%src(:,:,1) + &
  							gfedv2_g(butene_cis_2)%src(:,:,1) + &
  							gfedv2_g(pentene_2     )%src(:,:,1) + &
  							0.5*gfedv2_g(cyclopentene)%src(:,:,1) + &
  							0.5*gfedv2_g(butene_i    )%src(:,:,1)	  

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv3_g(butene_tr_2 )%src(:,:,1) + &
  							gfedv3_g(butene_cis_2)%src(:,:,1) + &
  							gfedv3_g(pentene_2     )%src(:,:,1) + &
  							0.5*gfedv3_g(cyclopentene)%src(:,:,1) + &
  							0.5*gfedv3_g(butene_i    )%src(:,:,1) ! NUWRF

   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv4_g(butene_tr_2 )%src(:,:,1) + &
  							gfedv4_g(butene_cis_2)%src(:,:,1) + &
  							gfedv4_g(pentene_2     )%src(:,:,1) + &
  							0.5*gfedv4_g(cyclopentene)%src(:,:,1) + &
  							0.5*gfedv4_g(butene_i    )%src(:,:,1) ! NUWRF

   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = qfed_g(butene_tr_2 )%src(:,:,1) + &
  							qfed_g(butene_cis_2)%src(:,:,1) + &
  							qfed_g(pentene_2     )%src(:,:,1) + &
  							0.5*qfed_g(cyclopentene)%src(:,:,1) + &
  							0.5*qfed_g(butene_i    )%src(:,:,1) ! NUWRF


   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1 ) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
  							  fwbawb_g(pentene_2   )%src(:,:,1) + &
                                                          fwbawb_g(butene_tr_2 )%src(:,:,1) + &
  							  fwbawb_g(butene_cis_2)%src(:,:,1) + &
  							  0.5*fwbawb_g(cyclopentene)%src(:,:,1) + &
  							  0.5*fwbawb_g(butene_i    )%src(:,:,1)	  

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!-ALD (m=44) (Agg=1)
!Acetald (m=44)	                 =>   ALD (Agg=1)
!Hydroxyacetaldehyde (m=60)	=>   ALD (Agg=1)  
!Propanal (m=58)	   	=>   ALD (Agg=1)
!Butanals (m=72)	   	=>   ALD (Agg=1)
!Hexanals (m=86)             	=>   ALD (Agg=1)
!Heptanals (m=100)	   	=>   ALD (Agg=1)
!Benzaldehyde (m=106)	   	=>   ALD (Agg=1)
if(spc_name == 'ALD') then
   ident = Acetald
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = bbbem_g(Acetald 	   )%src(:,:,1) + &
 							bbbem_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							bbbem_g(Propanal	   )%src(:,:,1) + &
  							bbbem_g(Butanals	   )%src(:,:,1) + &
  							bbbem_g(Hexanals	   )%src(:,:,1) + &
  							bbbem_g(Heptanals	   )%src(:,:,1) + &
  							bbbem_g(Benzaldehyde	   )%src(:,:,1) 	 

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv2_g(Acetald	    )%src(:,:,1) + &
  							gfedv2_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							gfedv2_g(Propanal	    )%src(:,:,1) + &
  							gfedv2_g(Butanals	    )%src(:,:,1) + &
  							gfedv2_g(Hexanals	    )%src(:,:,1) + &
 							gfedv2_g(Heptanals	    )%src(:,:,1) + &
  							gfedv2_g(Benzaldehyde	    )%src(:,:,1)	  

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv3_g(Acetald	    )%src(:,:,1) + &
  							gfedv3_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							gfedv3_g(Propanal	    )%src(:,:,1) + &
  							gfedv3_g(Butanals	    )%src(:,:,1) + &
  							gfedv3_g(Hexanals	    )%src(:,:,1) + &
 							gfedv3_g(Heptanals	    )%src(:,:,1) + &
  							gfedv3_g(Benzaldehyde	    )%src(:,:,1) ! NUWRF

   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv4_g(Acetald	    )%src(:,:,1) + &
  							gfedv4_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							gfedv4_g(Propanal	    )%src(:,:,1) + &
  							gfedv4_g(Butanals	    )%src(:,:,1) + &
  							gfedv4_g(Hexanals	    )%src(:,:,1) + &
 							gfedv4_g(Heptanals	    )%src(:,:,1) + &
  							gfedv4_g(Benzaldehyde	    )%src(:,:,1) ! NUWRF

   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = qfed_g(Acetald	    )%src(:,:,1) + &
  							qfed_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							qfed_g(Propanal	    )%src(:,:,1) + &
  							qfed_g(Butanals	    )%src(:,:,1) + &
  							qfed_g(Hexanals	    )%src(:,:,1) + &
 							qfed_g(Heptanals	    )%src(:,:,1) + &
  							qfed_g(Benzaldehyde	    )%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
							  fwbawb_g(Acetald	      )%src(:,:,1) + &
               						  fwbawb_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							  fwbawb_g(Propanal	      )%src(:,:,1) + &
  							  fwbawb_g(Butanals	      )%src(:,:,1) + &
  							  fwbawb_g(Hexanals	      )%src(:,:,1) + &
 							  fwbawb_g(Heptanals	      )%src(:,:,1) + &
  							  fwbawb_g(Benzaldehyde       )%src(:,:,1)


   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!-HC3 (m=44)
!2_Me-butane (m=72)	 => HC3 (Agg=1.11) (2,2-dimethylbutane)
!C3H8 (m=44)		 => HC3 (Agg=0.57) (Propane)
!n_butane (m=58)	 => HC3 (Agg=1.11)
!i-butane (m=58)	 => HC3 (Agg=1.11)
!Methanol (m=32)	 => HC3 (Agg=0.49)
!Ethanol (m=46) 	 => HC3 (Agg=1.37)
!C2H2 (m=26)		 => HC3 (Agg=0.41) (Acetylene)
if(spc_name == 'HC3') then
   ident = butane_2_Me
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = 1.11*bbbem_g(butane_2_Me)%src(:,:,1) + &
  							0.57*bbbem_g(C3H8       )%src(:,:,1) + &
  							1.11*bbbem_g(butane_n   )%src(:,:,1) + &
  							1.11*bbbem_g(butane_i   )%src(:,:,1) + &
  							0.49*bbbem_g(Methanol   )%src(:,:,1) + &
  							1.37*bbbem_g(Ethanol    )%src(:,:,1) + &
  							0.41*bbbem_g(C2H2       )%src(:,:,1) 	 

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = 1.11*gfedv2_g(butane_2_Me)%src(:,:,1) + &
  							0.57*gfedv2_g(C3H8	    )%src(:,:,1) + &
  							1.11*gfedv2_g(butane_n   )%src(:,:,1) + &
  							1.11*gfedv2_g(butane_i   )%src(:,:,1) + &
  							0.49*gfedv2_g(Methanol   )%src(:,:,1) + &
 							1.37*gfedv2_g(Ethanol    )%src(:,:,1) + &
  							0.41*gfedv2_g(C2H2	    )%src(:,:,1)	  

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = 1.11*gfedv3_g(butane_2_Me)%src(:,:,1) + &
  							0.57*gfedv3_g(C3H8	    )%src(:,:,1) + &
  							1.11*gfedv3_g(butane_n   )%src(:,:,1) + &
  							1.11*gfedv3_g(butane_i   )%src(:,:,1) + &
  							0.49*gfedv3_g(Methanol   )%src(:,:,1) + &
 							1.37*gfedv3_g(Ethanol    )%src(:,:,1) + &
  							0.41*gfedv3_g(C2H2	    )%src(:,:,1) ! NUWRF

   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = 1.11*gfedv4_g(butane_2_Me)%src(:,:,1) + &
  							0.57*gfedv4_g(C3H8	    )%src(:,:,1) + &
  							1.11*gfedv4_g(butane_n   )%src(:,:,1) + &
  							1.11*gfedv4_g(butane_i   )%src(:,:,1) + &
  							0.49*gfedv4_g(Methanol   )%src(:,:,1) + &
 							1.37*gfedv4_g(Ethanol    )%src(:,:,1) + &
  							0.41*gfedv4_g(C2H2	    )%src(:,:,1) ! NUWRF

   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = 1.11*qfed_g(butane_2_Me)%src(:,:,1) + &
  							0.57*qfed_g(C3H8	    )%src(:,:,1) + &
  							1.11*qfed_g(butane_n   )%src(:,:,1) + &
  							1.11*qfed_g(butane_i   )%src(:,:,1) + &
  							0.49*qfed_g(Methanol   )%src(:,:,1) + &
 							1.37*qfed_g(Ethanol    )%src(:,:,1) + &
  							0.41*qfed_g(C2H2	    )%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
                                                          1.11*fwbawb_g(butane_2_Me)%src(:,:,1) + &
			       				  0.57*fwbawb_g(C3H8       )%src(:,:,1) + &
			        			  1.11*fwbawb_g(butane_n   )%src(:,:,1) + &
				        		  1.11*fwbawb_g(butane_i   )%src(:,:,1) + &
					        	  0.49*fwbawb_g(Methanol   )%src(:,:,1) + &
							  1.37*fwbawb_g(Ethanol    )%src(:,:,1) + &
							  0.41*fwbawb_g(C2H2       )%src(:,:,1) 
!                                                        
!                                                                                                               

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!-HC5 (m=72)
!n_pentane (m=72)	=> HC5 (Agg=0.97)
!n_hexane (m=86)	=> HC5 (Agg=0.97)
!isohexanes (m=86)	=> HC5 (Agg=0.97)
!1_propanol (m=60)	=> HC5 (Agg=1.07)
!2-propanol (m=60)	=> HC5 (Agg=1.07)

if(spc_name == 'HC5') then
   ident = pentane_n
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
!   print*,'EMK: iespc,bburn,emiss_spc_name = ',iespc,bburn,emiss_spc_name(iespc,bburn)
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.97*bbbem_g(pentane_n )%src(:,:,1) + &
  							0.97*bbbem_g(hexane_n  )%src(:,:,1) + &
  							0.97*bbbem_g(isohexanes)%src(:,:,1) + &
  							1.07*bbbem_g(Propanol_1)%src(:,:,1) + &
  							1.07*bbbem_g(propanol_2)%src(:,:,1) 	

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.97*gfedv2_g(pentane_n )%src(:,:,1) + &
  							0.97*gfedv2_g(hexane_n  )%src(:,:,1) + &
  							0.97*gfedv2_g(isohexanes)%src(:,:,1) + &
  							1.07*gfedv2_g(Propanol_1)%src(:,:,1) + &
  							1.07*gfedv2_g(propanol_2)%src(:,:,1)

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.97*gfedv3_g(pentane_n )%src(:,:,1) + &
  							0.97*gfedv3_g(hexane_n  )%src(:,:,1) + &
  							0.97*gfedv3_g(isohexanes)%src(:,:,1) + &
  							1.07*gfedv3_g(Propanol_1)%src(:,:,1) + &
  							1.07*gfedv3_g(propanol_2)%src(:,:,1) ! NUWRF

   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.97*gfedv4_g(pentane_n )%src(:,:,1) + &
  							0.97*gfedv4_g(hexane_n  )%src(:,:,1) + &
  							0.97*gfedv4_g(isohexanes)%src(:,:,1) + &
  							1.07*gfedv4_g(Propanol_1)%src(:,:,1) + &
  							1.07*gfedv4_g(propanol_2)%src(:,:,1) ! NUWRF

   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.97*qfed_g(pentane_n )%src(:,:,1) + &
  							0.97*qfed_g(hexane_n  )%src(:,:,1) + &
  							0.97*qfed_g(isohexanes)%src(:,:,1) + &
  							1.07*qfed_g(Propanol_1)%src(:,:,1) + &
  							1.07*qfed_g(propanol_2)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1
      
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
							  0.97*fwbawb_g(pentane_n )%src(:,:,1) + &
  							  0.97*fwbawb_g(hexane_n  )%src(:,:,1) + &
  							  0.97*fwbawb_g(isohexanes)%src(:,:,1) + &
  							  1.07*fwbawb_g(Propanol_1)%src(:,:,1) + &
  							  1.07*fwbawb_g(propanol_2)%src(:,:,1)

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name
   return
endif

end subroutine convert_AeM_to_wrf_reac
