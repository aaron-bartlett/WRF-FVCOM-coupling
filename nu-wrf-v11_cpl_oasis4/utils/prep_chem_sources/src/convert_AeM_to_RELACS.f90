!###############################################################################!
!  CCATT-BRAMS/MCGA-CPTEC/WRF emission model    CPTEC/INPE                      !
!  Version 1.0.0: 12/nov/2010                                                   !
!  Coded by Saulo Freitas and Karla Longo                                       !
!  Contact: gmai@cptec.inpe.br - http://meioambiente.cptec.inpe.br              !
!###############################################################################!

subroutine convert_AeM_to_relacs(isp,iespc,ident,source_type,spc_name_dummy)
use chem1_list!, only : alke, bio,ora2,aro,ket,alka,ald
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
character (len=*), intent(in)  :: spc_name_dummy  !kml 

!print*,alke, bio,ora2,aro,ket,alka,ald, "x"; pause 1111


!--     AeM    |    RACM 

!- NOx	=>   NO
!      change MP 201108 : NO=NOX/2,   NO2 = NOX/2 in molecules (40% and 60% in mass)
if(spc_name_dummy == 'NO') then
   ident = NOx
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.4* bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.4*gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.4*gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.4*gfedv4_g(ident)%src(:,:,1) ! NUWRF
   ! Special QFED rule: Split NOx into 90% NO and 10% NO2
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.9*qfed_g(ident)%src(:,:,1) ! NUWRF
   
   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
      
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + 0.4*fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1


   print*,'==> converted from AeM - found for ',spc_name_dummy

   return
endif
if(spc_name_dummy == 'NO2') then
   ident = NOx
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.6* bbbem_g(ident)%src(:,:,1)
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.6*gfedv2_g(ident)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.6*gfedv3_g(ident)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.6*gfedv4_g(ident)%src(:,:,1) ! NUWRF
   ! Special QFED rule: Split NOx into 90% NO and 10% NO2
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 0.1*qfed_g(ident)%src(:,:,1) ! NUWRF
   
   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
      
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + 0.6*fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1


   print*,'==> converted from AeM - found for ',spc_name_dummy

   return
endif


!----------------------------------------------------------------------------
!- C2H6 =>	ETH
if(spc_name_dummy == 'ETH') then
   ident = C2H6
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(ident)%src(:,:,1)

   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(ident)%src(:,:,1)

   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(ident)%src(:,:,1) ! NUWRF

   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(ident)%src(:,:,1) ! NUWRF

   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   =qfed_g(ident)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(ident)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name_dummy
   return
endif

!----------------------------------------------------------------------------
!- ALKE  (m=33)  
!- C2H4	(m=28)				=> ALKE
!tr_2_butene (m=56) (butene_tr_2)   	=> ALKE
!cis_2_butene (m=56)     (butene_cis_2) => ALKE
!-2_pentene (m=70)                      => ALKE
!cyclopentene (m=68)                    => ALKE
!i-butene (butene_i)(m=56)              => ALKE
!butadiene (m=54)              		=> ALKE
!pentadienes (m=68)            		=> ALKE
!cyclopentadiene (m=66)       		=> ALKE
!hexadienes (m=82)             		=> ALKE
! C3H6 (m=42)   			=> ALKE
! 1_butene (m=56)   			=> ALKE
! 1_pentene (m=70)  			=> ALKE
! 2_Me_Butene (m=70)			=> ALKE
! 4_me_1_pentene (m=84) 		=> ALKE
! 2_me_1_pentene (m=84)			=> ALKE
! 1_hexene (m=84)			=> ALKE
! octenes (m=112)			=> ALKE



if(spc_name_dummy == 'ALKE') then
   ident = ALKE   
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  33/28*bbbem_g(C2H4)%src(:,:,1)+ &
                                                        33/56*bbbem_g(butene_tr_2 )%src(:,:,1) + &
                                                        33/56*bbbem_g(butene_cis_2)%src(:,:,1) + &
                                                        33/70*bbbem_g(pentene_2   )%src(:,:,1) + &
                                                        33/68*bbbem_g(cyclopentene)%src(:,:,1) + &
                                                        33/56*bbbem_g(butene_i    )%src(:,:,1) + &
                                                        33/54*bbbem_g(butadiene      )%src(:,:,1) + &
                                                        33/68*bbbem_g(pentadienes    )%src(:,:,1) + &
                                                        33/66*bbbem_g(cyclopentadiene)%src(:,:,1) + &
                                                        33/82*bbbem_g(hexadienes     )%src(:,:,1) + &
							33/42*bbbem_g(C3H6	      )%src(:,:,1) + &
 							33/56*bbbem_g(butene_1      )%src(:,:,1) + &
							33/70*bbbem_g(pentene_1     )%src(:,:,1) + &
  							33/70*bbbem_g(Butene_2_Me   )%src(:,:,1) + &
  							33/84*bbbem_g(pentene_4_me_1)%src(:,:,1) + &
  							33/84*bbbem_g(pentene_2_me_1)%src(:,:,1) + &
 							33/84*bbbem_g(hexene_1      )%src(:,:,1) + & 
							33/112*bbbem_g(octenes       )%src(:,:,1)   

   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 33/28*gfedv2_g(C2H4)%src(:,:,1)+ &
                                                        33/56*gfedv2_g(butene_tr_2 )%src(:,:,1) + &
                                                        33/56*gfedv2_g(butene_cis_2)%src(:,:,1) + &
                                                        33/70*gfedv2_g(pentene_2     )%src(:,:,1) + &
                                                        33/68*gfedv2_g(cyclopentene)%src(:,:,1) + &
                                                        33/56*gfedv2_g(butene_i    )%src(:,:,1) + &
							33/54*gfedv2_g(butadiene	)%src(:,:,1) + &
  							33/68*gfedv2_g(pentadienes	)%src(:,:,1) + &
  							33/66*gfedv2_g(cyclopentadiene)%src(:,:,1) + &
  							33/82*gfedv2_g(hexadienes	)%src(:,:,1) + &
							33/42*gfedv2_g(C3H6	       )%src(:,:,1) + &
  							33/56*gfedv2_g(butene_1      )%src(:,:,1) + &
  							33/70*gfedv2_g(pentene_1     )%src(:,:,1) + &
   							33/70*gfedv2_g(Butene_2_Me   )%src(:,:,1) + &
   							33/84*gfedv2_g(pentene_4_me_1)%src(:,:,1) + &
  							33/84*gfedv2_g(pentene_2_me_1)%src(:,:,1) + &
 							33/84*gfedv2_g(hexene_1      )%src(:,:,1) + & 
							33/112*gfedv2_g(octenes       )%src(:,:,1) 
							    	       

   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 33/28*gfedv3_g(C2H4)%src(:,:,1)+ &
                                                        33/56*gfedv3_g(butene_tr_2 )%src(:,:,1) + &
                                                        33/56*gfedv3_g(butene_cis_2)%src(:,:,1) + &
                                                        33/70*gfedv3_g(pentene_2     )%src(:,:,1) + &
                                                        33/68*gfedv3_g(cyclopentene)%src(:,:,1) + &
                                                        33/56*gfedv3_g(butene_i    )%src(:,:,1) + &
							33/54*gfedv3_g(butadiene	)%src(:,:,1) + &
  							33/68*gfedv3_g(pentadienes	)%src(:,:,1) + &
  							33/66*gfedv3_g(cyclopentadiene)%src(:,:,1) + &
  							33/82*gfedv3_g(hexadienes	)%src(:,:,1) + &
							33/42*gfedv3_g(C3H6	       )%src(:,:,1) + &
  							33/56*gfedv3_g(butene_1      )%src(:,:,1) + &
  							33/70*gfedv3_g(pentene_1     )%src(:,:,1) + &
   							33/70*gfedv3_g(Butene_2_Me   )%src(:,:,1) + &
   							33/84*gfedv3_g(pentene_4_me_1)%src(:,:,1) + &
  							33/84*gfedv3_g(pentene_2_me_1)%src(:,:,1) + &
 							33/84*gfedv3_g(hexene_1      )%src(:,:,1) + & 
							33/112*gfedv3_g(octenes       )%src(:,:,1) ! NUWRF

   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 33/28*gfedv4_g(C2H4)%src(:,:,1)+ &
                                                        33/56*gfedv4_g(butene_tr_2 )%src(:,:,1) + &
                                                        33/56*gfedv4_g(butene_cis_2)%src(:,:,1) + &
                                                        33/70*gfedv4_g(pentene_2     )%src(:,:,1) + &
                                                        33/68*gfedv4_g(cyclopentene)%src(:,:,1) + &
                                                        33/56*gfedv4_g(butene_i    )%src(:,:,1) + &
							33/54*gfedv4_g(butadiene	)%src(:,:,1) + &
  							33/68*gfedv4_g(pentadienes	)%src(:,:,1) + &
  							33/66*gfedv4_g(cyclopentadiene)%src(:,:,1) + &
  							33/82*gfedv4_g(hexadienes	)%src(:,:,1) + &
							33/42*gfedv4_g(C3H6	       )%src(:,:,1) + &
  							33/56*gfedv4_g(butene_1      )%src(:,:,1) + &
  							33/70*gfedv4_g(pentene_1     )%src(:,:,1) + &
   							33/70*gfedv4_g(Butene_2_Me   )%src(:,:,1) + &
   							33/84*gfedv4_g(pentene_4_me_1)%src(:,:,1) + &
  							33/84*gfedv4_g(pentene_2_me_1)%src(:,:,1) + &
 							33/84*gfedv4_g(hexene_1      )%src(:,:,1) + & 
							33/112*gfedv4_g(octenes       )%src(:,:,1) ! NUWRF

   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 33/28*qfed_g(C2H4)%src(:,:,1)+ &
                                                        33/56*qfed_g(butene_tr_2 )%src(:,:,1) + &
                                                        33/56*qfed_g(butene_cis_2)%src(:,:,1) + &
                                                        33/70*qfed_g(pentene_2     )%src(:,:,1) + &
                                                        33/68*qfed_g(cyclopentene)%src(:,:,1) + &
                                                        33/56*qfed_g(butene_i    )%src(:,:,1) + &
							33/54*qfed_g(butadiene	)%src(:,:,1) + &
  							33/68*qfed_g(pentadienes	)%src(:,:,1) + &
  							33/66*qfed_g(cyclopentadiene)%src(:,:,1) + &
  							33/82*qfed_g(hexadienes	)%src(:,:,1) + &
							33/42*qfed_g(C3H6	       )%src(:,:,1) + &
  							33/56*qfed_g(butene_1      )%src(:,:,1) + &
  							33/70*qfed_g(pentene_1     )%src(:,:,1) + &
   							33/70*qfed_g(Butene_2_Me   )%src(:,:,1) + &
   							33/84*qfed_g(pentene_4_me_1)%src(:,:,1) + &
  							33/84*qfed_g(pentene_2_me_1)%src(:,:,1) + &
 							33/84*qfed_g(hexene_1      )%src(:,:,1) + & 
							33/112*qfed_g(octenes       )%src(:,:,1) ! NUWRF
							    	       

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) +&
                                                        33/28*fwbawb_g(C2H4)%src(:,:,1)+ &
                                                        33/56*fwbawb_g(butene_tr_2 )%src(:,:,1) + &
                                                        33/56*fwbawb_g(butene_cis_2)%src(:,:,1) + &
                                                        33/70*fwbawb_g(pentene_2     )%src(:,:,1) + &
                                                        33/68*fwbawb_g(cyclopentene)%src(:,:,1) + &
                                                        33/56*fwbawb_g(butene_i    )%src(:,:,1) + &
							33/54*fwbawb_g(butadiene	)%src(:,:,1) + &
  							33/68*fwbawb_g(pentadienes	)%src(:,:,1) + &
  							33/66*fwbawb_g(cyclopentadiene)%src(:,:,1) + &
  							33/82*fwbawb_g(hexadienes	)%src(:,:,1) + &
							33/42*fwbawb_g(C3H6	       )%src(:,:,1) + &
  							33/56*fwbawb_g(butene_1      )%src(:,:,1) + &
  							33/70*fwbawb_g(pentene_1     )%src(:,:,1) + &
   							33/70*fwbawb_g(Butene_2_Me   )%src(:,:,1) + &
   							33/84*fwbawb_g(pentene_4_me_1)%src(:,:,1) + &
  							33/84*fwbawb_g(pentene_2_me_1)%src(:,:,1) + &
 							33/84*fwbawb_g(hexene_1      )%src(:,:,1) + & 
							33/112*fwbawb_g(octenes       )%src(:,:,1)+ &
                                                        68/56*fwbawb_g(pentene_2   )%src(:,:,1) 
   
							  

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name_dummy
   return
endif



!----------------------------------------------------------------------------
!Formaldehyde	      =>  HCHO
if(spc_name_dummy == 'HCHO') then
   ident = Formaldehyde
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
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

   print*,'==> converted from AeM - found for ',spc_name_dummy
   return
endif

!----------------------------------------------------------------------------
!HFo	=>       ORA1
if(spc_name_dummy == 'ORA1') then
   ident = HFo
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
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
! BIO (m=68)
!Isoprene  (m=68)		=>	BIO
!API (m=136) 0.5*terpenes	=>	BIO
!LIM (m=136) 0.2*terpenes	=>	BIO

if(spc_name_dummy == 'BIO') then
   ident = BIO 
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(Isoprene)%src(:,:,1) + 68/136*0.7*bbbem_g(TERPENES)%src(:,:,1) 
                                                            
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv2_g(Isoprene)%src(:,:,1) + 68/136*0.7*gfedv2_g(TERPENES)%src(:,:,1)

   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv3_g(Isoprene)%src(:,:,1) + 68/136*0.7*gfedv3_g(TERPENES)%src(:,:,1) ! NUWRF

   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   = gfedv4_g(Isoprene)%src(:,:,1) + 68/136*0.7*gfedv4_g(TERPENES)%src(:,:,1) ! NUWRF

   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   = qfed_g(Isoprene)%src(:,:,1) + 68/136*0.7*qfed_g(TERPENES)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1) &
                                                                                          + fwbawb_g(Isoprene)%src(:,:,1) + &
											  68/136*0.7*fwbawb_g(TERPENES)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name_dummy
   return
endif

!---------------------------------------------------------------------------
!-ORA2 (m=60)
! HAc (m=60)
! Propanoic (m=74)
! HAc	     =>   ORA2
! Propanoic  =>   ORA2
if(spc_name_dummy == 'ORA2') then
   ident = ORA2
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy

   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  bbbem_g(HAc	     )%src(:,:,1)  + &
                                                           60/74* bbbem_g(Propanoic)%src(:,:,1)     
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  gfedv2_g(HAc      )%src(:,:,1) + &
   							    60/74*gfedv2_g(Propanoic)%src(:,:,1)
   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  gfedv3_g(HAc      )%src(:,:,1) + &
   							    60/74*gfedv3_g(Propanoic)%src(:,:,1) ! NUWRF
   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  gfedv4_g(HAc      )%src(:,:,1) + &
   							    60/74*gfedv4_g(Propanoic)%src(:,:,1) ! NUWRF
   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  qfed_g(HAc      )%src(:,:,1) + &
   							    60/74*qfed_g(Propanoic)%src(:,:,1) ! NUWRF
							    
   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF

   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
                                                                                            fwbawb_g(HAc      )%src(:,:,1) + &
   							                                   60/74*fwbawb_g(Propanoic)%src(:,:,1)
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name_dummy
   return
endif

!----------------------------------------------------------------------------
!-ARO (m=98)
!-TOL (m=92)			=> ARO
! benzene (m=82)		=> ARO
! toluene (m=92)		=> ARO
! ethylbenzene (m=106)		=> ARO
! styrene (m=104)		=> ARO
! PAH (m=164)			=> ARO
! XYL 	(m=106)			=> ARO
! Phenol (m=94)	(CSL,(m=108))	=> ARO

if(spc_name_dummy == 'ARO') then
   ident = ARO
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy

   if(use_bbem   == 1) emiss_g(iespc)%src_bburn(:,:,1)   = 98/82*bbbem_g(benzene	)%src(:,:,1) + &
                                                           98/92*bbbem_g(toluene	)%src(:,:,1) + &
                                                           98/106*bbbem_g(ethylbenzene)%src(:,:,1) + &
                                                           98/104*bbbem_g(styrene	)%src(:,:,1) + &
                                                           98/164*bbbem_g(PAH 	)%src(:,:,1)  + &
							   98/106*bbbem_g(xylenes)%src(:,:,1) + &
							   98/94*bbbem_g(Phenol)%src(:,:,1) 
      
   if(use_gfedv2 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  98/82*gfedv2_g(benzene	)%src(:,:,1) + &
                                                           98/92*gfedv2_g(toluene	)%src(:,:,1) + &
                                                           98/106*gfedv2_g(ethylbenzene)%src(:,:,1) + &
                                                           98/104*gfedv2_g(styrene	)%src(:,:,1) + &
                                                           98/164*gfedv2_g(PAH 	)%src(:,:,1)  + &
							   98/106*gfedv2_g(xylenes)%src(:,:,1) + &
							   98/94*gfedv2_g(Phenol)%src(:,:,1) 

   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  98/82*gfedv3_g(benzene	)%src(:,:,1) + &
                                                           98/92*gfedv3_g(toluene	)%src(:,:,1) + &
                                                           98/106*gfedv3_g(ethylbenzene)%src(:,:,1) + &
                                                           98/104*gfedv3_g(styrene	)%src(:,:,1) + &
                                                           98/164*gfedv3_g(PAH 	)%src(:,:,1)  + &
							   98/106*gfedv3_g(xylenes)%src(:,:,1) + &
							   98/94*gfedv3_g(Phenol)%src(:,:,1) ! NUWRF

   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  98/82*gfedv4_g(benzene	)%src(:,:,1) + &
                                                           98/92*gfedv4_g(toluene	)%src(:,:,1) + &
                                                           98/106*gfedv4_g(ethylbenzene)%src(:,:,1) + &
                                                           98/104*gfedv4_g(styrene	)%src(:,:,1) + &
                                                           98/164*gfedv4_g(PAH 	)%src(:,:,1)  + &
							   98/106*gfedv4_g(xylenes)%src(:,:,1) + &
							   98/94*gfedv4_g(Phenol)%src(:,:,1) ! NUWRF

   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1)   =  98/82*qfed_g(benzene	)%src(:,:,1) + &
                                                           98/92*qfed_g(toluene	)%src(:,:,1) + &
                                                           98/106*qfed_g(ethylbenzene)%src(:,:,1) + &
                                                           98/104*qfed_g(styrene	)%src(:,:,1) + &
                                                           98/164*qfed_g(PAH 	)%src(:,:,1)  + &
							   98/106*qfed_g(xylenes)%src(:,:,1) + &
							   98/94*qfed_g(Phenol)%src(:,:,1) ! NUWRF

							   
   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
                                                          98/82*fwbawb_g(benzene	)%src(:,:,1) + &
                                                           98/92*fwbawb_g(toluene	)%src(:,:,1) + &
                                                           98/106*fwbawb_g(ethylbenzene)%src(:,:,1) + &
                                                           98/104*fwbawb_g(styrene	)%src(:,:,1) + &
                                                           98/164*fwbawb_g(PAH 	)%src(:,:,1)  + &
							   98/106*fwbawb_g(xylenes)%src(:,:,1) + &
							   98/94*fwbawb_g(Phenol)%src(:,:,1) 

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name_dummy
   return
endif



!----------------------------------------------------------------------------
!-KET (m=72)
!-Acetone (m=58)	  =>KET
!-2_Butanone (m=72)       =>KET
!-2_3_Butanedione (m=86)  =>KET
!-Pentanones (m=86)	  =>KET
!-Hexanones (m=100)	  =>KET
!-Heptanones (m=114)	  =>KET
!-Octanones (m=128	  =>KET
if(spc_name_dummy == 'KET') then
   ident = KET
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = 72/58*bbbem_g(Acetone        )%src(:,:,1) + &
  							bbbem_g(Butanone_2     )%src(:,:,1) + &
  							72/86*bbbem_g(Butanedione_2_3)%src(:,:,1) + &
  							72/86*bbbem_g(Pentanones     )%src(:,:,1) + &
   							72/100*bbbem_g(Hexanones      )%src(:,:,1) + &
   							72/114*bbbem_g(Heptanones     )%src(:,:,1) + &
  							72/128*bbbem_g(Octanones      )%src(:,:,1)		

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = 72/58*gfedv2_g(Acetone	)%src(:,:,1) + &
  							gfedv2_g(Butanone_2	)%src(:,:,1) + &
  							72/86*gfedv2_g(Butanedione_2_3)%src(:,:,1) + &
  							72/86*gfedv2_g(Pentanones	)%src(:,:,1) + &
   							72/100*gfedv2_g(Hexanones	)%src(:,:,1) + &
   							72/114*gfedv2_g(Heptanones	)%src(:,:,1) + &
  							72/128*gfedv2_g(Octanones	)%src(:,:,1)		

   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1) = 72/58*gfedv3_g(Acetone	)%src(:,:,1) + &
  							gfedv3_g(Butanone_2	)%src(:,:,1) + &
  							72/86*gfedv3_g(Butanedione_2_3)%src(:,:,1) + &
  							72/86*gfedv3_g(Pentanones	)%src(:,:,1) + &
   							72/100*gfedv3_g(Hexanones	)%src(:,:,1) + &
   							72/114*gfedv3_g(Heptanones	)%src(:,:,1) + &
  							72/128*gfedv3_g(Octanones	)%src(:,:,1) ! NUWRF

   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1) = 72/58*gfedv4_g(Acetone	)%src(:,:,1) + &
  							gfedv4_g(Butanone_2	)%src(:,:,1) + &
  							72/86*gfedv4_g(Butanedione_2_3)%src(:,:,1) + &
  							72/86*gfedv4_g(Pentanones	)%src(:,:,1) + &
   							72/100*gfedv4_g(Hexanones	)%src(:,:,1) + &
   							72/114*gfedv4_g(Heptanones	)%src(:,:,1) + &
  							72/128*gfedv4_g(Octanones	)%src(:,:,1) ! NUWRF

   if(use_qfed== 1) emiss_g(iespc)%src_bburn(:,:,1) = 72/58*qfed_g(Acetone	)%src(:,:,1) + &
  							qfed_g(Butanone_2	)%src(:,:,1) + &
  							72/86*qfed_g(Butanedione_2_3)%src(:,:,1) + &
  							72/86*qfed_g(Pentanones	)%src(:,:,1) + &
   							72/100*qfed_g(Hexanones	)%src(:,:,1) + &
   							72/114*qfed_g(Heptanones	)%src(:,:,1) + &
  							72/128*qfed_g(Octanones	)%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
                                                          72/58*fwbawb_g(Acetone	  )%src(:,:,1) + &
  							  fwbawb_g(Butanone_2	  )%src(:,:,1) + &
  							  72/86*fwbawb_g(Butanedione_2_3)%src(:,:,1) + &
  							  72/86*fwbawb_g(Pentanones	  )%src(:,:,1) + &
   							  72/100*fwbawb_g(Hexanones	  )%src(:,:,1) + &
  							  72/114*fwbawb_g(Heptanones	  )%src(:,:,1) + &
  							  72/128*fwbawb_g(Octanones	  )%src(:,:,1)

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name_dummy
   return
endif

!----------------------------------------------------------------------------
!- ALKA (m=62)
!heptane (m=100) (HC8 (m=114))  	=> ALKA
!Butanols (m=74) (HC8 (m=114))  	=> ALKA	
!cyclopentanol (m=86) (HC8 (m=114)) 	=> ALKA
!2_Me-butane (m=72) HC3 (m=44) 		=> ALKA
!C3H8 (m=44)	    HC3 (m=44) 		=> ALKA
!n_butane (m=58)    HC3 (m=44) 		=> ALKA
!i-butane (m=58)    HC3 (m=44) 		=> ALKA
!Methanol (m=32)    HC3 (m=44) 		=> ALKA
!Ethanol (m=46)     HC3 (m=44) 		=> ALKA
!C2H2 (m=26)        HC3 (m=44) 		=> ALKA
!n_pentane (m=72)   HC5 (m=72) 		=> ALKA
!n_hexane (m=86)    HC5 (m=72) 		=> ALKA
!isohexanes (m=86)  HC5 (m=72) 		=> ALKA
!1_propanol (m=60)  HC5 (m=72)	  	=> ALKA
!2-propanol (m=60)  HC5 (m=72) 		=> ALKA


if(spc_name_dummy == 'ALKA') then
   ident = ALKA
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = 62/100*bbbem_g(heptane         )%src(:,:,1) + &
  							62/74*bbbem_g(Butanols     )%src(:,:,1) + &
  							62/86*bbbem_g(cyclopentanol)%src(:,:,1) + &
							62/72*bbbem_g(butane_2_Me)%src(:,:,1) + &
  							62/44*bbbem_g(C3H8       )%src(:,:,1) + &
  							62/58*bbbem_g(butane_n   )%src(:,:,1) + &
  							62/58*bbbem_g(butane_i   )%src(:,:,1) + &
  							62/32*bbbem_g(Methanol   )%src(:,:,1) + &
  							62/46*bbbem_g(Ethanol    )%src(:,:,1) + &
  							62/26*bbbem_g(C2H2       )%src(:,:,1) + &
							62/72*bbbem_g(pentane_n )%src(:,:,1) + &
  							62/86*bbbem_g(hexane_n  )%src(:,:,1) + &
  							62/86*bbbem_g(isohexanes)%src(:,:,1) + &
  							62/60*bbbem_g(Propanol_1)%src(:,:,1) + &
  							62/60*bbbem_g(propanol_2)%src(:,:,1) 	
	 	   

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = 62/100*gfedv2_g(heptane         )%src(:,:,1) + &
  							62/74*gfedv2_g(Butanols     )%src(:,:,1) + &
  							62/86*gfedv2_g(cyclopentanol)%src(:,:,1) + &
							62/72*gfedv2_g(butane_2_Me)%src(:,:,1) + &
  							62/44*gfedv2_g(C3H8       )%src(:,:,1) + &
  							62/58*gfedv2_g(butane_n   )%src(:,:,1) + &
  							62/58*gfedv2_g(butane_i   )%src(:,:,1) + &
  							62/32*gfedv2_g(Methanol   )%src(:,:,1) + &
  							62/46*gfedv2_g(Ethanol    )%src(:,:,1) + &
  							62/26*gfedv2_g(C2H2       )%src(:,:,1) + &
							62/72*gfedv2_g(pentane_n )%src(:,:,1) + &
  							62/86*gfedv2_g(hexane_n  )%src(:,:,1) + &
  							62/86*gfedv2_g(isohexanes)%src(:,:,1) + &
  							62/60*gfedv2_g(Propanol_1)%src(:,:,1) + &
  							62/60*gfedv2_g(propanol_2)%src(:,:,1) 	
							

   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1) = 62/100*gfedv3_g(heptane         )%src(:,:,1) + &
  							62/74*gfedv3_g(Butanols     )%src(:,:,1) + &
  							62/86*gfedv3_g(cyclopentanol)%src(:,:,1) + &
							62/72*gfedv3_g(butane_2_Me)%src(:,:,1) + &
  							62/44*gfedv3_g(C3H8       )%src(:,:,1) + &
  							62/58*gfedv3_g(butane_n   )%src(:,:,1) + &
  							62/58*gfedv3_g(butane_i   )%src(:,:,1) + &
  							62/32*gfedv3_g(Methanol   )%src(:,:,1) + &
  							62/46*gfedv3_g(Ethanol    )%src(:,:,1) + &
  							62/26*gfedv3_g(C2H2       )%src(:,:,1) + &
							62/72*gfedv3_g(pentane_n )%src(:,:,1) + &
  							62/86*gfedv3_g(hexane_n  )%src(:,:,1) + &
  							62/86*gfedv3_g(isohexanes)%src(:,:,1) + &
  							62/60*gfedv3_g(Propanol_1)%src(:,:,1) + &
  							62/60*gfedv3_g(propanol_2)%src(:,:,1) ! NUWRF

   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1) = 62/100*gfedv4_g(heptane         )%src(:,:,1) + &
  							62/74*gfedv4_g(Butanols     )%src(:,:,1) + &
  							62/86*gfedv4_g(cyclopentanol)%src(:,:,1) + &
							62/72*gfedv4_g(butane_2_Me)%src(:,:,1) + &
  							62/44*gfedv4_g(C3H8       )%src(:,:,1) + &
  							62/58*gfedv4_g(butane_n   )%src(:,:,1) + &
  							62/58*gfedv4_g(butane_i   )%src(:,:,1) + &
  							62/32*gfedv4_g(Methanol   )%src(:,:,1) + &
  							62/46*gfedv4_g(Ethanol    )%src(:,:,1) + &
  							62/26*gfedv4_g(C2H2       )%src(:,:,1) + &
							62/72*gfedv4_g(pentane_n )%src(:,:,1) + &
  							62/86*gfedv4_g(hexane_n  )%src(:,:,1) + &
  							62/86*gfedv4_g(isohexanes)%src(:,:,1) + &
  							62/60*gfedv4_g(Propanol_1)%src(:,:,1) + &
  							62/60*gfedv4_g(propanol_2)%src(:,:,1) ! NUWRF

   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1) = 62/100*qfed_g(heptane         )%src(:,:,1) + &
  							62/74*qfed_g(Butanols     )%src(:,:,1) + &
  							62/86*qfed_g(cyclopentanol)%src(:,:,1) + &
							62/72*qfed_g(butane_2_Me)%src(:,:,1) + &
  							62/44*qfed_g(C3H8       )%src(:,:,1) + &
  							62/58*qfed_g(butane_n   )%src(:,:,1) + &
  							62/58*qfed_g(butane_i   )%src(:,:,1) + &
  							62/32*qfed_g(Methanol   )%src(:,:,1) + &
  							62/46*qfed_g(Ethanol    )%src(:,:,1) + &
  							62/26*qfed_g(C2H2       )%src(:,:,1) + &
							62/72*qfed_g(pentane_n )%src(:,:,1) + &
  							62/86*qfed_g(hexane_n  )%src(:,:,1) + &
  							62/86*qfed_g(isohexanes)%src(:,:,1) + &
  							62/60*qfed_g(Propanol_1)%src(:,:,1) + &
  							62/60*qfed_g(propanol_2)%src(:,:,1) ! NUWRF
							
   if(use_bbem  == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
  							62/100*fwbawb_g(heptane         )%src(:,:,1) + &
  							62/74*fwbawb_g(Butanols     )%src(:,:,1) + &
  							62/86*fwbawb_g(cyclopentanol)%src(:,:,1) + &
							62/72*fwbawb_g(butane_2_Me)%src(:,:,1) + &
  							62/44*fwbawb_g(C3H8       )%src(:,:,1) + &
  							62/58*fwbawb_g(butane_n   )%src(:,:,1) + &
  							62/58*fwbawb_g(butane_i   )%src(:,:,1) + &
  							62/32*fwbawb_g(Methanol   )%src(:,:,1) + &
  							62/46*fwbawb_g(Ethanol    )%src(:,:,1) + &
  							62/26*fwbawb_g(C2H2       )%src(:,:,1) + &
							62/72*fwbawb_g(pentane_n )%src(:,:,1) + &
  							62/86*fwbawb_g(hexane_n  )%src(:,:,1) + &
  							62/86*fwbawb_g(isohexanes)%src(:,:,1) + &
  							62/60*fwbawb_g(Propanol_1)%src(:,:,1) + &
  							62/60*fwbawb_g(propanol_2)%src(:,:,1) 	
							

   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name_dummy
   return
endif

!----------------------------------------------------------------------------
!-ALD (m=44)
!Acetald (m=44)	                 =>   ALD
!Hydroxyacetaldehyde (m=60)	=>   ALD
!Propanal (m=58)	   	=>   ALD
!Butanals (m=72)	   	=>   ALD
!Hexanals (m=86)             	=>   ALD
!Heptanals (m=100)	   	=>   ALD
!Benzaldehyde (m=106)	   	=>   ALD
if(spc_name_dummy == 'ALD') then
   ident = ALD
   iespc=iespc+1
   emiss_spc_name(iespc,bburn)           = spc_name_dummy
   if(use_bbem	== 1) emiss_g(iespc)%src_bburn(:,:,1) = bbbem_g(Acetald 	   )%src(:,:,1) + &
 							44/60*bbbem_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							44/58*bbbem_g(Propanal	   )%src(:,:,1) + &
  							44/72*bbbem_g(Butanals	   )%src(:,:,1) + &
  							44/86*bbbem_g(Hexanals	   )%src(:,:,1) + &
  							44/100*bbbem_g(Heptanals	   )%src(:,:,1) + &
  							44/106*bbbem_g(Benzaldehyde	   )%src(:,:,1) 	 

   if(use_gfedv2== 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv2_g(Acetald	    )%src(:,:,1) + &
  							44/60*gfedv2_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							44/58*gfedv2_g(Propanal	    )%src(:,:,1) + &
  							44/72*gfedv2_g(Butanals	    )%src(:,:,1) + &
  							44/86*gfedv2_g(Hexanals	    )%src(:,:,1) + &
 							44/100*gfedv2_g(Heptanals	    )%src(:,:,1) + &
  							44/106*gfedv2_g(Benzaldehyde	    )%src(:,:,1)	  

   if(use_gfedv3 == 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv3_g(Acetald	    )%src(:,:,1) + &
  							44/60*gfedv3_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							44/58*gfedv3_g(Propanal	    )%src(:,:,1) + &
  							44/72*gfedv3_g(Butanals	    )%src(:,:,1) + &
  							44/86*gfedv3_g(Hexanals	    )%src(:,:,1) + &
 							44/100*gfedv3_g(Heptanals	    )%src(:,:,1) + &
  							44/106*gfedv3_g(Benzaldehyde	    )%src(:,:,1) ! NUWRF

   if(use_gfedv4 == 1) emiss_g(iespc)%src_bburn(:,:,1) = gfedv4_g(Acetald	    )%src(:,:,1) + &
  							44/60*gfedv4_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							44/58*gfedv4_g(Propanal	    )%src(:,:,1) + &
  							44/72*gfedv4_g(Butanals	    )%src(:,:,1) + &
  							44/86*gfedv4_g(Hexanals	    )%src(:,:,1) + &
 							44/100*gfedv4_g(Heptanals	    )%src(:,:,1) + &
  							44/106*gfedv4_g(Benzaldehyde	    )%src(:,:,1) ! NUWRF

   if(use_qfed == 1) emiss_g(iespc)%src_bburn(:,:,1) = qfed_g(Acetald	    )%src(:,:,1) + &
  							44/60*qfed_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							44/58*qfed_g(Propanal	    )%src(:,:,1) + &
  							44/72*qfed_g(Butanals	    )%src(:,:,1) + &
  							44/86*qfed_g(Hexanals	    )%src(:,:,1) + &
 							44/100*qfed_g(Heptanals	    )%src(:,:,1) + &
  							44/106*qfed_g(Benzaldehyde	    )%src(:,:,1) ! NUWRF

   if(use_bbem   == 1 .or. use_gfedv2 == 1 .or. use_gfedv3 == 1 .or. use_gfedv4 == 1 .or. use_qfed == 1) found_emiss_spc(iespc,bburn) = 1 ! NUWRF
   
   !- this emission is included in the "antro" category 
   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') emiss_g(iespc)%src_antro(:,:,1) = emiss_g(iespc)%src_antro(:,:,1)+ &
							  fwbawb_g(Acetald	      )%src(:,:,1) + &
               						  44/60*fwbawb_g(Hydroxyacetaldehyde)%src(:,:,1) + &
  							  44/58*fwbawb_g(Propanal	      )%src(:,:,1) + &
  							  44/72*fwbawb_g(Butanals	      )%src(:,:,1) + &
  							  44/86*fwbawb_g(Hexanals	      )%src(:,:,1) + &
 							  44/100*fwbawb_g(Heptanals	      )%src(:,:,1) + &
  							  44/106*fwbawb_g(Benzaldehyde       )%src(:,:,1)


   if(use_fwbawb == 1 .and. trim(source_type) == 'antro') found_emiss_spc(iespc,antro)    = 1

   print*,'==> converted from AeM - found for ',spc_name_dummy
   return
endif


end subroutine convert_AeM_to_relacs
