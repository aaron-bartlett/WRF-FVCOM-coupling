!###############################################################################!
!  CCATT-BRAMS/MCGA-CPTEC/WRF emission model    CPTEC/INPE                      !
!  Version 1.0.0: 12/nov/2010                                                   !
!  Coded by Saulo Freitas and Karla Longo                                       !
!  Contact: gmai@cptec.inpe.br - http://meioambiente.cptec.inpe.br              !
!###############################################################################!

subroutine convert_AeM_to_CB07(isp,iespc,ident,source_type,spc_name)
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


implicit none  
integer, intent(in) :: isp
integer, intent(inout) :: iespc,ident
character (len=*) source_type
character (len=*), intent(in)  :: spc_name  !kml 


!--     AeM    |    CB2002

!- NOx	=>   NO
if(spc_name == 'NO') then
   ident = NOx
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)

   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name

   return
endif

!----------------------------------------------------------------------------
!- Phenol =>	CRES
if(spc_name == 'CRES') then
   ident = Phenol
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif


!----------------------------------------------------------------------------
!- C2H4 =>	ETH
if(spc_name == 'ETH') then
   ident = C2H4
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif


!----------------------------------------------------------------------------
!- Formaldehyde	=> FORM
if(spc_name == 'FORM') then
   ident = Formaldehyde
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)

   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!- Xylenes   =>     XYL
if(spc_name == 'XYL') then
   ident = xylenes
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)

   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!- Methanol   =>     MEOH
if(spc_name == 'MEOH') then
   ident = Methanol
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)

   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!- Ethanol   =>     ETOH
if(spc_name == 'ETOH') then
   ident = Ethanol
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)

   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!-----------------------------------------------------------------------------
!Isoprene  =>	ISOP
if(spc_name == 'ISOP') then
   ident = Isoprene
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(ident)%src(:,:,1) ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
! styrene       => MGLY (methylglyoxal and other aromatics)
! PAH           => MGLY
!if(spc_name == 'MGLY') then
!   ident = MGLY
!   iespc=iespc+1
!   emiss_spc_name(iespc,bburn)           = spc_name
!
!   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(styrene   )%src(:,:,1) + &
!                                                            bbbem_g(PAH   )%src(:,:,1)     
!      
!   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  gfedv2_g(styrene  )%src(:,:,1) + &
!   							    gfedv2_g(PAH   )%src(:,:,1)
!							    
!   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  gfedv3_g(styrene  )%src(:,:,1) + &
!   							    gfedv3_g(PAH   )%src(:,:,1) ! NUWRF
!   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  qfed_g(styrene  )%src(:,:,1) + &
!   							    qfed_g(PAH   )%src(:,:,1) ! NUWRF
!
!   !- this emission is included in the "antro" category 
!   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
!                                                                                            fwbawb_g(styrene  )%src(:,:,1) + &
!   							                                    fwbawb_g(PAH   )%src(:,:,1)
!   found_emiss_spc(iespc,bburn)    = 1
!   print*,'==> converted from AeM - found for ',spc_name
!   return
!endif

!----------------------------------------------------------------------------
!toluene	  =>      TOL
!ethylbenzene	  =>      TOL + PAR
if(spc_name == 'TOL') then
   ident = toluene
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name

   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =      bbbem_g(toluene	)%src(:,:,1) + &
                                                          0.868*bbbem_g(ethylbenzene)%src(:,:,1)
      
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =      gfedv2_g(toluene	 )%src(:,:,1) + &
   				                          0.868*gfedv2_g(ethylbenzene)%src(:,:,1)
   
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =      gfedv3_g(toluene	 )%src(:,:,1) + &
        0.868*gfedv3_g(ethylbenzene)%src(:,:,1) ! NUWRF
   
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =      gfedv4_g(toluene	 )%src(:,:,1) + &
        0.868*gfedv4_g(ethylbenzene)%src(:,:,1) ! NUWRF
   
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   =      qfed_g(toluene	 )%src(:,:,1) + &
        0.868*qfed_g(ethylbenzene)%src(:,:,1) ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) =    &
                                                          emiss_g(iespc)%src_antro(:,:,1)    + &
   							  fwbawb_g(toluene	)%src(:,:,1) + &
   				                    0.868*fwbawb_g(ethylbenzene )%src(:,:,1)
   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!-C3H6	          =>	OLE +  PAR
!-1_butene	  =>	OLE +2*PAR
!-1_pentene	  =>	OLE +3*PAR
!-2_me_Butene     =>    OLE +3*PAR
!-4_me_1_pentene  =>    OLE +4*PAR
!-2_me_1_pentene  =>    OLE +4*PAR
!-1_hexene	  =>    OLE +4*PAR
!-octenes	  =>    OLE +6*PAR
!-cis_2_butene	  =>    OLE +2*PAR 
!-2_pentene       =>    OLE +3*PAR
!-i_buten         =>    OLE +2*PAR
if(spc_name == 'OLE') then
   ident = C3H6
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) =0.666*bbbem_g(C3H6	      )%src(:,:,1) + &
  						       0.500*bbbem_g(butene_1      )%src(:,:,1) + &
  						       0.400*bbbem_g(pentene_1     )%src(:,:,1) + &
   						       0.400*bbbem_g(Butene_2_Me   )%src(:,:,1) + &
   						       0.333*bbbem_g(pentene_4_me_1)%src(:,:,1) + &
  						       0.333*bbbem_g(pentene_2_me_1)%src(:,:,1) + &
 						       0.333*bbbem_g(hexene_1      )%src(:,:,1) + & 
						       0.250*bbbem_g(octenes       )%src(:,:,1) + &   	       
  						       0.500*bbbem_g(butene_cis_2)%src(:,:,1) + &
  						       0.400*bbbem_g(pentene_2   )%src(:,:,1) + &
  						       0.500*bbbem_g(butene_i    )%src(:,:,1)	  

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) =0.666*gfedv2_g(C3H6	       )%src(:,:,1) + &
  						       0.500*gfedv2_g(butene_1      )%src(:,:,1) + &
  						       0.400*gfedv2_g(pentene_1     )%src(:,:,1) + &
   						       0.400*gfedv2_g(Butene_2_Me   )%src(:,:,1) + &
   						       0.333*gfedv2_g(pentene_4_me_1)%src(:,:,1) + &
  						       0.333*gfedv2_g(pentene_2_me_1)%src(:,:,1) + &
 						       0.333*gfedv2_g(hexene_1      )%src(:,:,1) + & 
						       0.250*gfedv2_g(octenes       )%src(:,:,1) + &   	       
  						       0.500*gfedv2_g(butene_cis_2)%src(:,:,1) + &
  						       0.400*gfedv2_g(pentene_2   )%src(:,:,1) + &
  						       0.500*gfedv2_g(butene_i    )%src(:,:,1)	  

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) =0.666*gfedv3_g(C3H6	       )%src(:,:,1) + &
  						       0.500*gfedv3_g(butene_1      )%src(:,:,1) + &
  						       0.400*gfedv3_g(pentene_1     )%src(:,:,1) + &
   						       0.400*gfedv3_g(Butene_2_Me   )%src(:,:,1) + &
   						       0.333*gfedv3_g(pentene_4_me_1)%src(:,:,1) + &
  						       0.333*gfedv3_g(pentene_2_me_1)%src(:,:,1) + &
 						       0.333*gfedv3_g(hexene_1      )%src(:,:,1) + & 
						       0.250*gfedv3_g(octenes       )%src(:,:,1) + &   	       
  						       0.500*gfedv3_g(butene_cis_2)%src(:,:,1) + &
  						       0.400*gfedv3_g(pentene_2   )%src(:,:,1) + &
  						       0.500*gfedv3_g(butene_i    )%src(:,:,1)	  ! NUWRF
   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) =0.666*gfedv4_g(C3H6	       )%src(:,:,1) + &
  						       0.500*gfedv4_g(butene_1      )%src(:,:,1) + &
  						       0.400*gfedv4_g(pentene_1     )%src(:,:,1) + &
   						       0.400*gfedv4_g(Butene_2_Me   )%src(:,:,1) + &
   						       0.333*gfedv4_g(pentene_4_me_1)%src(:,:,1) + &
  						       0.333*gfedv4_g(pentene_2_me_1)%src(:,:,1) + &
 						       0.333*gfedv4_g(hexene_1      )%src(:,:,1) + & 
						       0.250*gfedv4_g(octenes       )%src(:,:,1) + &   	       
  						       0.500*gfedv4_g(butene_cis_2)%src(:,:,1) + &
  						       0.400*gfedv4_g(pentene_2   )%src(:,:,1) + &
  						       0.500*gfedv4_g(butene_i    )%src(:,:,1)	  ! NUWRF
   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) =0.666*qfed_g(C3H6	       )%src(:,:,1) + &
  						       0.500*qfed_g(butene_1      )%src(:,:,1) + &
  						       0.400*qfed_g(pentene_1     )%src(:,:,1) + &
   						       0.400*qfed_g(Butene_2_Me   )%src(:,:,1) + &
   						       0.333*qfed_g(pentene_4_me_1)%src(:,:,1) + &
  						       0.333*qfed_g(pentene_2_me_1)%src(:,:,1) + &
 						       0.333*qfed_g(hexene_1      )%src(:,:,1) + & 
						       0.250*qfed_g(octenes       )%src(:,:,1) + &   	       
  						       0.500*qfed_g(butene_cis_2)%src(:,:,1) + &
  						       0.400*qfed_g(pentene_2   )%src(:,:,1) + &
  						       0.500*qfed_g(butene_i    )%src(:,:,1) ! NUWRF

   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) =     &
                                                          emiss_g(iespc)%src_antro(:,:,1)     + &
                                                    0.666*fwbawb_g(C3H6 	 )%src(:,:,1) + &
  						    0.500*fwbawb_g(butene_1	 )%src(:,:,1) + &
  						    0.400*fwbawb_g(pentene_1	 )%src(:,:,1) + &
    						    0.400*fwbawb_g(Butene_2_Me   )%src(:,:,1) + &
   						    0.333*fwbawb_g(pentene_4_me_1)%src(:,:,1) + &
  						    0.333*fwbawb_g(pentene_2_me_1)%src(:,:,1) + &
 						    0.333*fwbawb_g(hexene_1	 )%src(:,:,1) + & 
						    0.250*fwbawb_g(octenes	 )%src(:,:,1) + &	 
  					            0.500*fwbawb_g(butene_cis_2)%src(:,:,1) + &
  					            0.400*fwbawb_g(pentene_2   )%src(:,:,1) + &
  					            0.500*fwbawb_g(butene_i    )%src(:,:,1)	  

   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!butadiene	       => DIEN
!pentadienes	       => DIEN
!cyclopentadiene       => DIEN
!hexadienes	       => DIEN
!if(spc_name == 'DIEN') then
!   ident = DIEN
!   iespc=iespc+1
!   emiss_spc_name(iespc,bburn)           = spc_name
!   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = bbbem_g(butadiene      )%src(:,:,1) + &
!  							bbbem_g(pentadienes    )%src(:,:,1) + &
!  							bbbem_g(cyclopentadiene)%src(:,:,1) + &
!  							bbbem_g(hexadienes     )%src(:,:,1) 		
!
!   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv2_g(butadiene	)%src(:,:,1) + &
!  							gfedv2_g(pentadienes	)%src(:,:,1) + &
!  							gfedv2_g(cyclopentadiene)%src(:,:,1) + &
!  							gfedv2_g(hexadienes	)%src(:,:,1)		
!
!   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv3_g(butadiene	)%src(:,:,1) + &
!  							gfedv3_g(pentadienes	)%src(:,:,1) + &
!  							gfedv3_g(cyclopentadiene)%src(:,:,1) + &
!  							gfedv3_g(hexadienes	)%src(:,:,1) ! NUWRF
!   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = qfed_g(butadiene	)%src(:,:,1) + &
!  						       qfed_g(pentadienes	)%src(:,:,1) + &
!  						       qfed_g(cyclopentadiene)%src(:,:,1) + &
!  						       qfed_g(hexadienes	)%src(:,:,1) ! NUWRF
!
!   
   !- this emission is included in the "antro" category 
!   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) =      &
!                                                         emiss_g(iespc)%src_antro(:,:,1)      + &
!                                                         fwbawb_g(butadiene	  )%src(:,:,1) + &
!  							  fwbawb_g(pentadienes	  )%src(:,:,1) + &
!  							  fwbawb_g(cyclopentadiene)%src(:,:,1) + &
!  							  fwbawb_g(hexadienes	  )%src(:,:,1)
!
!   found_emiss_spc(iespc,bburn)    = 1
!   print*,'==> converted from AeM - found for ',spc_name
!   return
!endif
!
!----------------------------------------------------------------------------
!Acetald	        =>   ALD2
!Propanal	   	=>   ALD2 +  PAR
!Butanals	   	=>   ALD2 +2*PAR
!Hexanals	   	=>   ALD2 +4*PAR
!Heptanals	   	=>   ALD2 +5*PAR
!tr_2_butene	        => 2*ALD2 

if(spc_name == 'ALD2') then
   ident = Acetald
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = bbbem_g(Acetald 	   )%src(:,:,1) + &
  						  0.759*bbbem_g(Propanal	   )%src(:,:,1) + &
  						  0.611*bbbem_g(Butanals	   )%src(:,:,1) + &
  					          0.523*bbbem_g(Hexanals	   )%src(:,:,1) + &
  					          0.458*bbbem_g(Heptanals	   )%src(:,:,1) + &
  					        	bbbem_g(butene_tr_2	   )%src(:,:,1) 	 

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv2_g(Acetald	    )%src(:,:,1) + &
  						  0.759*gfedv2_g(Propanal	    )%src(:,:,1) + &
  						  0.611*gfedv2_g(Butanals	    )%src(:,:,1) + &
  					          0.523*gfedv2_g(Hexanals	    )%src(:,:,1) + &
 						  0.458*gfedv2_g(Heptanals	    )%src(:,:,1) + &
  						        gfedv2_g(butene_tr_2	    )%src(:,:,1)	  

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv3_g(Acetald	    )%src(:,:,1) + &
  						  0.759*gfedv3_g(Propanal	    )%src(:,:,1) + &
  						  0.611*gfedv3_g(Butanals	    )%src(:,:,1) + &
  					          0.523*gfedv3_g(Hexanals	    )%src(:,:,1) + &
 						  0.458*gfedv3_g(Heptanals	    )%src(:,:,1) + &
  						        gfedv3_g(butene_tr_2	    )%src(:,:,1) ! NUWRF
   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv4_g(Acetald	    )%src(:,:,1) + &
  						  0.759*gfedv4_g(Propanal	    )%src(:,:,1) + &
  						  0.611*gfedv4_g(Butanals	    )%src(:,:,1) + &
  					          0.523*gfedv4_g(Hexanals	    )%src(:,:,1) + &
 						  0.458*gfedv4_g(Heptanals	    )%src(:,:,1) + &
  						        gfedv4_g(butene_tr_2	    )%src(:,:,1) ! NUWRF
   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = qfed_g(Acetald	    )%src(:,:,1) + &
  						  0.759*qfed_g(Propanal	    )%src(:,:,1) + &
  						  0.611*qfed_g(Butanals	    )%src(:,:,1) + &
  					          0.523*qfed_g(Hexanals	    )%src(:,:,1) + &
 						  0.458*qfed_g(Heptanals	    )%src(:,:,1) + &
  						        qfed_g(butene_tr_2	    )%src(:,:,1) ! NUWRF

   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) =          &
                                                          emiss_g(iespc)%src_antro(:,:,1)        +   &
							  fwbawb_g(Acetald	      )%src(:,:,1) + &
  						    0.759*fwbawb_g(Propanal	      )%src(:,:,1) + &
  					            0.611*fwbawb_g(Butanals	      )%src(:,:,1) + &
  					            0.523*fwbawb_g(Hexanals	      )%src(:,:,1) + &
 				                    0.458*fwbawb_g(Heptanals	      )%src(:,:,1) + &
  							  fwbawb_g(butene_tr_2       )%src(:,:,1)


   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif

!----------------------------------------------------------------------------
!
!ethylbenzene	 => TOL + PAR
!C3H8		 => 3*PAR
!n_butane	 => 4*PAR
!i-butane	 => 4*PAR
!n_pentane       => 5*PAR
!n_hexane	 => 6*PAR
!isohexanes	 => 6*PAR
!heptane         => 7*PAR
!C2H6            => 2*PAR
!Propanal   	 =>   ALD2 +  PAR
!Butanals	 =>   ALD2 +2*PAR
!Hexanals	 =>   ALD2 +4*PAR
!Heptanals	 =>   ALD2 +5*PAR
!-C3H6	          =>	OLE +  PAR
!-1_butene	  =>	OLE +2*PAR
!-1_pentene	  =>	OLE +3*PAR
!-2_me_Butene     =>    OLE +3*PAR
!-4_me_1_pentene  =>    OLE +4*PAR
!-2_me_1_pentene  =>    OLE +4*PAR
!-1_hexene	  =>    OLE +4*PAR
!-octenes	  =>    OLE +6*PAR
!-cis_2_butene	  =>    OLE +2*PAR 
!-2_pentene       =>    OLE +3*PAR
!-i_buten         =>    OLE +2*PAR

if(spc_name == 'PAR') then
   ident = C2H6
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) =0.132*bbbem_g(ethylbenzene)%src(:,:,1) + &
  						       0.955*bbbem_g(C3H8       )%src(:,:,1) + &
  						       0.965*bbbem_g(butane_n   )%src(:,:,1) + &
  						       0.965*bbbem_g(butane_i   )%src(:,:,1) + &
                                                       0.972*bbbem_g(pentane_n )%src(:,:,1) + &
  						       0.976*bbbem_g(hexane_n  )%src(:,:,1) + &
  						       0.976*bbbem_g(isohexanes)%src(:,:,1) + &
  						       0.980*bbbem_g(heptane   )%src(:,:,1) + &
  						       0.933*bbbem_g(C2H6      )%src(:,:,1) + &
  						       0.241*bbbem_g(Propanal	   )%src(:,:,1) + &
  						       0.389*bbbem_g(Butanals	   )%src(:,:,1) + &
  						       0.476*bbbem_g(Hexanals	   )%src(:,:,1) + &
  						       0.542*bbbem_g(Heptanals	   )%src(:,:,1) + &
                                                       0.333*bbbem_g(C3H6	      )%src(:,:,1) + &
  						       0.500*bbbem_g(butene_1      )%src(:,:,1) + &
  						       0.600*bbbem_g(pentene_1     )%src(:,:,1) + &
   						       0.600*bbbem_g(Butene_2_Me   )%src(:,:,1) + &
   						       0.666*bbbem_g(pentene_4_me_1)%src(:,:,1) + &
  						       0.666*bbbem_g(pentene_2_me_1)%src(:,:,1) + &
 						       0.666*bbbem_g(hexene_1      )%src(:,:,1) + & 
						       0.750*bbbem_g(octenes       )%src(:,:,1) + &   	       
  						       0.500*bbbem_g(butene_cis_2)%src(:,:,1) + &
  						       0.600*bbbem_g(pentene_2   )%src(:,:,1) + &
  						       0.500*bbbem_g(butene_i    )%src(:,:,1)	  

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.132*gfedv2_g(ethylbenzene)%src(:,:,1) + &
  							0.955*gfedv2_g(C3H8	    )%src(:,:,1) + &
  							0.965*gfedv2_g(butane_n   )%src(:,:,1) + &
					        	0.965*gfedv2_g(butane_i   )%src(:,:,1) + &
                                                        0.972*gfedv2_g(pentane_n )%src(:,:,1) + &
  							0.976*gfedv2_g(hexane_n  )%src(:,:,1) + &
  							0.976*gfedv2_g(isohexanes)%src(:,:,1) + &
  							0.980*gfedv2_g(heptane   )%src(:,:,1) + &
  							0.933*gfedv2_g(C2H6      )%src(:,:,1) + &
  							0.241*gfedv2_g(Propanal	    )%src(:,:,1) + &
  							0.389*gfedv2_g(Butanals	    )%src(:,:,1) + &
  							0.476*gfedv2_g(Hexanals	    )%src(:,:,1) + &
 							0.542*gfedv2_g(Heptanals	    )%src(:,:,1) + &
                                                        0.333*gfedv2_g(C3H6	       )%src(:,:,1) + &
  							0.500*gfedv2_g(butene_1      )%src(:,:,1) + &
  							0.600*gfedv2_g(pentene_1     )%src(:,:,1) + &
   							0.600*gfedv2_g(Butene_2_Me   )%src(:,:,1) + &
   							0.666*gfedv2_g(pentene_4_me_1)%src(:,:,1) + &
  							0.666*gfedv2_g(pentene_2_me_1)%src(:,:,1) + &
 							0.666*gfedv2_g(hexene_1      )%src(:,:,1) + & 
							0.750*gfedv2_g(octenes       )%src(:,:,1) + &   	       
  							0.500*gfedv2_g(butene_cis_2)%src(:,:,1) + &
  							0.600*gfedv2_g(pentene_2   )%src(:,:,1) + &
  							0.500*gfedv2_g(butene_i    )%src(:,:,1)	  

   if(use_gfedv3== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.132*gfedv3_g(ethylbenzene)%src(:,:,1) + &
  							0.955*gfedv3_g(C3H8	    )%src(:,:,1) + &
  							0.965*gfedv3_g(butane_n   )%src(:,:,1) + &
					        	0.965*gfedv3_g(butane_i   )%src(:,:,1) + &
                                                        0.972*gfedv3_g(pentane_n )%src(:,:,1) + &
  							0.976*gfedv3_g(hexane_n  )%src(:,:,1) + &
  							0.976*gfedv3_g(isohexanes)%src(:,:,1) + &
  							0.980*gfedv3_g(heptane   )%src(:,:,1) + &
  							0.933*gfedv3_g(C2H6      )%src(:,:,1) + &
  							0.241*gfedv3_g(Propanal	    )%src(:,:,1) + &
  							0.389*gfedv3_g(Butanals	    )%src(:,:,1) + &
  							0.476*gfedv3_g(Hexanals	    )%src(:,:,1) + &
 							0.542*gfedv3_g(Heptanals	    )%src(:,:,1) + &
                                                        0.333*gfedv3_g(C3H6	       )%src(:,:,1) + &
  							0.500*gfedv3_g(butene_1      )%src(:,:,1) + &
  							0.600*gfedv3_g(pentene_1     )%src(:,:,1) + &
   							0.600*gfedv3_g(Butene_2_Me   )%src(:,:,1) + &
   							0.666*gfedv3_g(pentene_4_me_1)%src(:,:,1) + &
  							0.666*gfedv3_g(pentene_2_me_1)%src(:,:,1) + &
 							0.666*gfedv3_g(hexene_1      )%src(:,:,1) + & 
							0.750*gfedv3_g(octenes       )%src(:,:,1) + &   	       
  							0.500*gfedv3_g(butene_cis_2)%src(:,:,1) + &
  							0.600*gfedv3_g(pentene_2   )%src(:,:,1) + &
  							0.500*gfedv3_g(butene_i    )%src(:,:,1)	  ! NUWRF
   if(use_gfedv4== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.132*gfedv4_g(ethylbenzene)%src(:,:,1) + &
  							0.955*gfedv4_g(C3H8	    )%src(:,:,1) + &
  							0.965*gfedv4_g(butane_n   )%src(:,:,1) + &
					        	0.965*gfedv4_g(butane_i   )%src(:,:,1) + &
                                                        0.972*gfedv4_g(pentane_n )%src(:,:,1) + &
  							0.976*gfedv4_g(hexane_n  )%src(:,:,1) + &
  							0.976*gfedv4_g(isohexanes)%src(:,:,1) + &
  							0.980*gfedv4_g(heptane   )%src(:,:,1) + &
  							0.933*gfedv4_g(C2H6      )%src(:,:,1) + &
  							0.241*gfedv4_g(Propanal	    )%src(:,:,1) + &
  							0.389*gfedv4_g(Butanals	    )%src(:,:,1) + &
  							0.476*gfedv4_g(Hexanals	    )%src(:,:,1) + &
 							0.542*gfedv4_g(Heptanals	    )%src(:,:,1) + &
                                                        0.333*gfedv4_g(C3H6	       )%src(:,:,1) + &
  							0.500*gfedv4_g(butene_1      )%src(:,:,1) + &
  							0.600*gfedv4_g(pentene_1     )%src(:,:,1) + &
   							0.600*gfedv4_g(Butene_2_Me   )%src(:,:,1) + &
   							0.666*gfedv4_g(pentene_4_me_1)%src(:,:,1) + &
  							0.666*gfedv4_g(pentene_2_me_1)%src(:,:,1) + &
 							0.666*gfedv4_g(hexene_1      )%src(:,:,1) + & 
							0.750*gfedv4_g(octenes       )%src(:,:,1) + &   	       
  							0.500*gfedv4_g(butene_cis_2)%src(:,:,1) + &
  							0.600*gfedv4_g(pentene_2   )%src(:,:,1) + &
  							0.500*gfedv4_g(butene_i    )%src(:,:,1)	  ! NUWRF
   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = 0.132*qfed_g(ethylbenzene)%src(:,:,1) + &
  							0.955*qfed_g(C3H8	    )%src(:,:,1) + &
  							0.965*qfed_g(butane_n   )%src(:,:,1) + &
					        	0.965*qfed_g(butane_i   )%src(:,:,1) + &
                                                        0.972*qfed_g(pentane_n )%src(:,:,1) + &
  							0.976*qfed_g(hexane_n  )%src(:,:,1) + &
  							0.976*qfed_g(isohexanes)%src(:,:,1) + &
  							0.980*qfed_g(heptane   )%src(:,:,1) + &
  							0.933*qfed_g(C2H6      )%src(:,:,1) + &
  							0.241*qfed_g(Propanal	    )%src(:,:,1) + &
  							0.389*qfed_g(Butanals	    )%src(:,:,1) + &
  							0.476*qfed_g(Hexanals	    )%src(:,:,1) + &
 							0.542*qfed_g(Heptanals	    )%src(:,:,1) + &
                                                        0.333*qfed_g(C3H6	       )%src(:,:,1) + &
  							0.500*qfed_g(butene_1      )%src(:,:,1) + &
  							0.600*qfed_g(pentene_1     )%src(:,:,1) + &
   							0.600*qfed_g(Butene_2_Me   )%src(:,:,1) + &
   							0.666*qfed_g(pentene_4_me_1)%src(:,:,1) + &
  							0.666*qfed_g(pentene_2_me_1)%src(:,:,1) + &
 							0.666*qfed_g(hexene_1      )%src(:,:,1) + & 
							0.750*qfed_g(octenes       )%src(:,:,1) + &   	       
  							0.500*qfed_g(butene_cis_2)%src(:,:,1) + &
  							0.600*qfed_g(pentene_2   )%src(:,:,1) + &
  							0.500*qfed_g(butene_i    )%src(:,:,1) ! NUWRF

   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1)  = &
                                                          emiss_g(iespc)%src_antro(:,:,1)  + &
                                                          0.132*fwbawb_g(ethylbenzene)%src(:,:,1) + &
							  0.955*fwbawb_g(C3H8       )%src(:,:,1) + &
							  0.965*fwbawb_g(butane_n   )%src(:,:,1) + &
							  0.965*fwbawb_g(butane_i   )%src(:,:,1) + &
                                                          0.972*fwbawb_g(pentane_n )%src(:,:,1) + &
  							  0.976*fwbawb_g(hexane_n  )%src(:,:,1) + &
  							  0.976*fwbawb_g(isohexanes)%src(:,:,1) + &
  							  0.980*fwbawb_g(heptane   )%src(:,:,1) + &
  							  0.933*fwbawb_g(C2H6      )%src(:,:,1) + &
  							  0.241*fwbawb_g(Propanal	      )%src(:,:,1) + &
  							  0.389*fwbawb_g(Butanals	      )%src(:,:,1) + &
  							  0.476*fwbawb_g(Hexanals	      )%src(:,:,1) + &
 							  0.542*fwbawb_g(Heptanals	      )%src(:,:,1) + &
                                                          0.333*fwbawb_g(C3H6 	 )%src(:,:,1) + &
  						          0.500*fwbawb_g(butene_1	 )%src(:,:,1) + &
  						          0.600*fwbawb_g(pentene_1	 )%src(:,:,1) + &
    						          0.600*fwbawb_g(Butene_2_Me   )%src(:,:,1) + &
   						          0.666*fwbawb_g(pentene_4_me_1)%src(:,:,1) + &
  						          0.666*fwbawb_g(pentene_2_me_1)%src(:,:,1) + &
 						          0.666*fwbawb_g(hexene_1	 )%src(:,:,1) + & 
						          0.750*fwbawb_g(octenes	 )%src(:,:,1) + &	 
  							  0.500*fwbawb_g(butene_cis_2)%src(:,:,1) + &
  							  0.600*fwbawb_g(pentene_2   )%src(:,:,1) + &
  							  0.500*fwbawb_g(butene_i    )%src(:,:,1)	  

   found_emiss_spc(iespc,bburn)    = 1
   print*,'==> converted from AeM - found for ',spc_name
   return
endif


end subroutine convert_AeM_to_CB07
