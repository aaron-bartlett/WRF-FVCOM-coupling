!###############################################################################!
!  CCATT-BRAMS/MCGA-CPTEC/WRF emission model    CPTEC/INPE                      !
!  Version 1.0.0: 12/nov/2010                                                   !
!  Coded by Saulo Freitas and Karla Longo                                       !
!  Contact: gmai@cptec.inpe.br - http://meioambiente.cptec.inpe.br              !
!###############################################################################!
!  Created by Rafael Stockler: 19/oct/2010
!  Version: 1.0.0
!  Rotinas para WRF
!---------------------------------------------------------------------------

#if RADM_WRF_FIM

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!---------------- WRF or FIM Pos-Processing        -----------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine write_cptec_to_wrf_fim(file_in,iday,imon,iyear)
use chem1_list
use grid_dims_out, only: grid_type  

implicit none
integer, parameter :: nsrc=4!,nspecies=79, 
real, allocatable :: var2d(:,:,:,:),plume_mean_fct(:,:,:),plume_firesize(:,:,:),volc_prop(:,:,:), volc_degass_prop(:,:,:)
character(len=*) file_in
character*240 file_out_prefix

integer :: iunit,ispc,isrc,iveg_agreg,i,iday,imon,iyear
character(len=20) read_spc_name,read_src_name,section,read_units,date
character(len=20) read_mean,read_veg_name,c_dummy
integer read_ident_chem_mec,read_ident_src,dummy,read_ident_veg,ii,iii
integer read_ident_aer,read_aer_mode, imode
real ash_size_distr_dummy(10)
integer nxp,nyp
real dep_glon(2), dep_glat(2)
integer     :: &
  antro   = 01, & ! anthropogenic sources
  bburn   = 02, & ! biomass burning sources 
  bioge   = 03, & ! biogenic sources ! must be equal to "nsrc"
  geoge   = 04    ! geoge sources ! must be equal to "nsrc"


iunit=2


if(trim(chemical_mechanism)/='RADM_WRF_FIM')then
 print*,'==> WRF or FIM output only prepared for RADM chemical mechanism'
 stop 'compile again using CHEM=RADM_WRF_FIM'
endif

! input file name

open(iunit,status='old',form='formatted',file=file_in) 

read(iunit,*)  nxp,(dep_glon(i),i=1,2)
read(iunit,*)  nyp,(dep_glat(i),i=1,2)
read(iunit,*)  date

print*,'-----------------------------------------------------'
print*,'-----------------------------------------------------'
print*,'WRF/FIM pos-proc=',nxp,nyp,trim(file_in)


allocate(var2d(nxp,nyp,nspecies,nsrc));var2d=0.
allocate(plume_mean_fct(nxp,nyp,4));plume_mean_fct=0.
allocate(plume_firesize(nxp,nyp,4));plume_firesize=0.
allocate(volc_prop(nxp,nyp,3));volc_prop=0.
allocate(volc_degass_prop(nxp,nyp,3));volc_degass_prop=0.
ii=0; iii=0
   do i=1,nspecies*nsrc*5
       read(iunit,*,end=100) section

         !- emission section  --------------------------------
       if(trim(section) == 'emission' 		&
         .OR. trim(section) == 'chemistry' ) then
	 
         read(iunit,*) , read_spc_name &
		       , read_ident_chem_mec &
		       , read_src_name       &
		       , read_ident_src      &
		       , read_units
         
	 ispc = read_ident_chem_mec
         isrc = read_ident_src
         print*,read_spc_name,read_src_name,ispc,isrc
         if(trim(read_spc_name)=='ASH' ) then
	   read(iunit,*) c_dummy
	   read(iunit,*) ash_size_distr_dummy(:)
	   print*,"ash_size_distr=",ash_size_distr_dummy
	 endif
	 CALL vfirec(iunit,var2d(1,1,ispc,isrc),nxp*nyp,'LIN')
	 
	 where( var2d(:,:,ispc,isrc) <0.)  var2d(:,:,ispc,isrc)=0.
         !print*,var2d
       
       elseif (trim(section) == 'aerosol') then
       
         READ(iunit,*)   read_spc_name       &
               , read_ident_aer      &
               , read_aer_mode       &
               , read_src_name       &
               , read_ident_src      &
               , read_units
          ispc = read_ident_aer
          imode= read_aer_mode
          isrc = read_ident_src         
	 
         print*,read_spc_name,read_src_name,ispc,isrc
         if(trim(read_spc_name)=='V_ASH1' ) then
	   read(iunit,*) c_dummy
	   read(iunit,*) ash_size_distr_dummy(:)
	   !print*,"ash_size_distr_dummy(:)=",ash_size_distr_dummy
	 endif
         
	 CALL vfirec(iunit,var2d(1,1,ispc,isrc),nxp*nyp,'LIN')
	 
	 where( var2d(:,:,ispc,isrc) <0.)  var2d(:,:,ispc,isrc)=0.
         !print*,var2d      
	
       !- plume section  --------------------------------
       elseif(trim(section) == 'plume') then

         read(iunit,*) , read_mean &
		       , dummy &
		       , read_veg_name       &
		       , read_ident_veg      &
		       , read_units
         print*,'reading ',read_mean,' for ',read_veg_name

	 iveg_agreg = read_ident_veg       
         if(trim(read_mean) == 'mean_fct' ) then
	   call vfirec(iunit,plume_mean_fct(:,:,iveg_agreg),nxp*nyp,'LIN')
	   where( plume_mean_fct(:,:,iveg_agreg) <0.) plume_mean_fct(:,:,iveg_agreg)=0.
           !print*,plume_mean_fct

	 elseif(trim(read_mean) == 'firesize' ) then
	   call vfirec(iunit,plume_firesize(:,:,iveg_agreg),nxp*nyp,'LIN')
 	   where( plume_firesize(:,:,iveg_agreg) <0.) plume_firesize(:,:,iveg_agreg)=0.
           !print*,plume_firesize
        endif

      !- volcanoes eruption section  --------------------------------
       elseif(trim(section) == 'volcanic-eruption') then
        ii=ii+1
	 read(iunit,*) , read_mean &
		       , read_units 
         print*,'reading ',read_mean,' for volcanoes eruption'
	 call vfirec(iunit,volc_prop(:,:,ii),nxp*nyp,'LIN')
 	 where( volc_prop(:,:,ii) <0.) volc_prop(:,:,ii)=0.
	 print*,'maxval=',maxval(volc_prop(:,:,ii))
      
      !- volcanoes  degassing section  --------------------------------
       elseif(trim(section) == 'volcanic-degassing') then
        iii=iii+1
	 read(iunit,*) , read_mean &
		       , read_units 
         print*,'reading ',read_mean,' for volcanoes degassing'
	 call vfirec(iunit,volc_degass_prop(:,:,iii),nxp*nyp,'LIN')
 	 where( volc_prop(:,:,iii) <0.) volc_degass_prop(:,:,iii)=0.
	 print*,'maxval=',maxval(volc_degass_prop(:,:,iii))	
       endif

    enddo
100 continue

!- split bburn emissions into flaming/smoldering parts
 call emis_flam_smold_WRF(nxp,nyp,nspecies,nsrc,var2d,plume_mean_fct,plume_firesize)


if(trim(grid_type)=='fim') then
    call write_to_FIM(nxp,nyp,nspecies,nsrc,var2d,plume_mean_fct&
           ,plume_firesize,volc_prop,file_in,iday,imon,iyear,volc_degass_prop&
	   ,ash_size_distr_dummy)

else
!- convert to WRF Format 
 call write_to_WRFCHEM(nxp,nyp,nspecies,nsrc,var2d,plume_mean_fct&
           ,plume_firesize,volc_prop,file_in,iday,imon,iyear,volc_degass_prop&
	   ,ash_size_distr_dummy)

endif


close(iunit)
end subroutine write_cptec_to_wrf_fim
!-------------------------------------------------------------------------
  subroutine emis_flam_smold_WRF(n2,n3,nspecies,nsrc,var2d,plume_mean_fct,plume_firesize)
    implicit none
    integer,intent(IN) :: n2,n3,nspecies,nsrc
    real,dimension(n2,n3) :: smold_frac 
    real,dimension(n2,n3,4) :: plume_mean_fct,plume_firesize
    real,dimension(n2,n3,nspecies,nsrc) :: var2d 
    
    integer     :: &
      antro   = 01, & ! anthropogenic sources
      bburn   = 02, & ! biomass burning sources 
      bioge   = 03    ! biogenic sources 
    integer, parameter :: nveg_agreg      = 4
    integer, parameter :: tropical_forest = 1
    integer, parameter :: boreal_forest   = 2
    integer, parameter :: savannah        = 3
    integer, parameter :: grassland       = 4

    integer iv,ispc,i,j
    
    !-----  
    !- calcula a emissao smoldering e fatores para obtencao da fracao
    !- flaming em funcao da emissao smoldering


      smold_frac(:,:) = 1.- ( plume_mean_fct(:,:,tropical_forest) + &
  			      plume_mean_fct(:,:,boreal_forest  ) + &
  			      plume_mean_fct(:,:,savannah	) + &
  			      plume_mean_fct(:,:,grassland	)   )	 
       do ispc = 1,nspecies

  	   
  	!- convert from 'total' emisson to 'smoldering'
	var2d(:,:,ispc,bburn) = smold_frac(:,:) * var2d(:,:,ispc,bburn) 
	 !do i=1,n2; do j=1,n3
	 ! if(smold_frac(i,j) < 1. .and. var2d(i,j,ispc,bburn) > 0.)  print*,smold_frac(i,j) ,var2d(i,j,ispc,bburn)
	 !enddo;enddo
  						     
       enddo

       !- convert from flaming fraction to relationship with phase smoldering emission
       do iv = 1, nveg_agreg
         plume_mean_fct(:,:,iv) = plume_mean_fct(:,:,iv)/(1.e-8+smold_frac(:,:))
       enddo
       
  end subroutine  emis_flam_smold_WRF    
!---------------------------------------------------------
subroutine write_to_FIM(n2,n3,nspecies,nsrc,var2d,plume_mean_fct &
                             ,plume_firesize,volc_prop,file_in,iday,imon,iyear&
			     ,volc_degass_prop,ash_size_distr)

                            
    use grid_dims_out
    use chem1_list, only : spc_name, weight, BBURN2,BBURN3,OC,BC,URBAN2,URBAN3
    use gocart_backgr, only :  gocart_bg_nspecies=>nspecies&
                              ,gocart_bg_spc_name=>spc_name&
			      ,nlevels_netcdf,gocart_bg_g,lev
    implicit none
    integer,intent(IN) :: n2,n3,nspecies,nsrc
    real,dimension(n2,n3) :: smold_frac
    real,dimension(n2,n3,4) :: plume_mean_fct,plume_firesize
    real,dimension(n2,n3,3) :: volc_prop,volc_degass_prop
    real,dimension(n2,n3,nspecies,nsrc) :: var2d
    real, dimension(n2,n3) :: dummy
    character(*) ::file_in
    character(len=240) :: fileAB,fileBB,fileBG
    integer     :: &
      antro   = 01, & ! anthropogenic sources
      bburn   = 02, & ! biomass burning sources
      bioge   = 03 ,&   ! biogenic sources
      geoge   = 04    ! biogenic sources
    integer, parameter :: nveg_agreg      = 4
    integer, parameter :: tropical_forest = 1
    integer, parameter :: boreal_forest   = 2
    integer, parameter :: savannah        = 3
    integer, parameter :: grassland       = 4

    integer, parameter :: nspecies_wrf = 29
    real :: ash_size_distr(10)
    integer ispc,isrc,ispc_wrf,iv,iveg,ifound,iday,imon,iyear
    character(len=20) ename
    integer itime,i,j
    real fx
!from     /home/poluicao/WRFV2-CHEM/chem/module_input_chem_data.F
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspecies_wrf) :: spc_name_wrf=(/ &
 'SO2     ' &!
,'NO2     ' &!
,'NO      ' &!
,'ALD     ' &!
,'HCHO    ' &!
,'ORA2    ' &!
,'NH3     '  &! XXXXX
,'HC3     ' &!
,'HC5     ' &!
,'HC8     ' &!
,'ETH     ' &!
,'CO      ' &!
,'OL2     ' &! XXXXX
,'OLT     ' &!
,'OLI     ' &!
,'TOL     ' &!
,'XYL     ' &!
,'KET     ' &!
,'CSL     ' &!
,'ISO     ' &!
,'BBURN2  ' &!PM2.5
,'BBURN3  ' &!PM10.
,'URBAN2  ' &!PM2.5
,'URBAN3  ' &!PM10.
,'OC      ' &!
,'BC      ' &
,'DMS     ' &
,'SO4     ' &! SULF changed SO4
,'ASH     ' &

/)
!
!,'PM25I ' &!
!,'PM25J ' &!
!,'SO4I	' &!
!,'SO4J	' &!
!,'NO3I	' &!
!,'NO3J	' &!
!,'ORGI	' &!
!,'ORGJ	' &!
!,'ECI	' &!
!,'ECJ	' &!
!,'PM10	'


!antro + bioge section
ename='XXXXXXXX'
itime = 0
!fileAB=file_in(1:len_trim(file_in)-4)//'-ab.bin'

!open (91,file=fileAB,form='unformatted')
!     write(91)nspecies_wrf
!     write(91)ename
!     write(91)itime
iv=0
do ispc_wrf=1,nspecies_wrf
  !print*,'spc_name_wrf(ispc_wrf)',spc_name_wrf(ispc_wrf)
  ifound=0
  do ispc = 1,nspecies
     !if(spc_alloc(src,ispc) == 0) cycle

     !print*,'YY    ',spc_name(ispc),spc_name_wrf(ispc_wrf)

     if(trim(spc_name_wrf(ispc_wrf)) == trim(spc_name(ispc))) then
     	 if(trim(spc_name_wrf(ispc_wrf)) == 'BBURN2' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'BBURN3') then
	    ifound=1
	    cycle  
	 endif
         
	 print*,'FOUND ',spc_name_wrf(ispc_wrf),spc_name(ispc)
         !-for WRF/FIM, 'biogenic' emissions for CO is not included
	 !- 
	 if(trim(spc_name(ispc)) == 'CO') then
	     dummy(:,:) = var2d(:,:,ispc,antro)
	 else
	     dummy(:,:) = var2d(:,:,ispc,antro)+var2d(:,:,ispc,bioge)
	 endif

	 if(trim(spc_name_wrf(ispc_wrf)) == 'URBAN2' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'URBAN3' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'OC'   .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'BC'  ) then
	    fx= 1.e+9  / 86400.	 !  ug/(m^2 sec)
	 else
	    fx= 1.e+3 * 1.e6 / weight(ispc) / 24. ! emission unit: mol/(km^2 hour)
	 endif

         fileAB=file_in(1:len_trim(file_in)-7)//'-'//trim(spc_name_wrf(ispc_wrf))//'-ab.bin'
	 open (91,file=fileAB,form='unformatted')
         write(91) dummy * fx
         iv=iv+1
	 ifound=1
	 close(91)
	 exit
     endif
  enddo
  !if(ifound == 0) then
   ! print*,'NOT FOUND, zero emission',spc_name_wrf(ispc_wrf)
  !   dummy(:,:) =0.
  !   write(91) dummy
  !   iv=iv+1
  !   ifound=1
  !endif

enddo
!close (91)
!print*,'Anthro + biogenic emission for wrf: ',fileAB

!bburn + plumerise section
!fileBB=file_in(1:len_trim(file_in)-4)//'-bb.bin'
!open (91,file=fileBB,form='unformatted')
!     write(91)nspecies_wrf + 8
!     write(91)ename
!     write(91)itime
iv=0
do ispc_wrf=1,nspecies_wrf

   ifound=0
   do ispc = 1,nspecies
     if(trim(spc_name_wrf(ispc_wrf)) == trim(spc_name(ispc))) then
     	 if(trim(spc_name_wrf(ispc_wrf)) == 'URBAN2' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'URBAN3') then
	    ifound=1
	    cycle
	 endif
        ! print*,' specie=',trim(spc_name(ispc))
         dummy(:,:) = var2d(:,:,ispc,bburn)


	 !-PM25 for wrf:

	 if(trim(spc_name_wrf(ispc_wrf)) == 'BBURN2' ) then

	 dummy(:,:) = max(0.,var2d(:,:,BBURN2,bburn)  -  &
	                     var2d(:,:,OC  ,bburn)  -  &
	                     var2d(:,:,BC  ,bburn)     )

	 endif

	 if(trim(spc_name_wrf(ispc_wrf)) == 'BBURN2' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'BBURN3' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'OC'   .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'BC'  ) then
	    fx= 1.e+9  / 86400.	 !  ug/(m^2 sec)
	 else
	    fx= 1.e+3 * 1.e6 / weight(ispc) / 24.! emission unit: mol/(km^2 hour)
	 endif

         print*,' specie=',trim(spc_name(ispc)),maxval(dummy*fx),maxloc(dummy*fx)!,,minval(dummy*fx)

         fileBB=file_in(1:len_trim(file_in)-7)//'-'//trim(spc_name_wrf(ispc_wrf))//'-bb.bin'
         open (91,file=fileBB,form='unformatted')
	 write(91) dummy*fx
	 iv=iv+1
	 ifound=1
        exit
     endif

  enddo


!  if(ifound == 0) then
!     !print*,'not found, zero emission ',spc_name_wrf(ispc_wrf)
!     dummy(:,:) =0.
!     write(91) dummy
!     iv=iv+1
!     ifound=1
!  endif

enddo
fileBB=file_in(1:len_trim(file_in)-7)//'-plume.bin'
open (91,file=fileBB,form='unformatted')
do iveg=1,nveg_agreg
   write(91) plume_mean_fct(:,:,iveg)
	 iv=iv+1
enddo
do iveg=1,nveg_agreg
   write(91) plume_firesize(:,:,iveg)
	 iv=iv+1
enddo
close (91)
!print*,'BB: total 2d fields=',iv
!print*,' Biomass burning emission + plumerise data for wrf: ',fileBB

!-------------------------- volcanoes -------------------------
print*,'=> Volcanoes Eruption Section ---------------------------'

if(use_volcanoes == 1) then
   ename='volc-eruption'
   fileBB=file_in(1:len_trim(file_in)-4)//'-volc.bin'
   open (91,file=fileBB,form='unformatted')
     write(91) 5 !fields
     write(91)ename
     write(91)begin_eruption
     write(91) ash_size_distr(:)
     
     iv=0
     do ispc_wrf=1,nspecies_wrf
      do ispc = 1,nspecies
        if(trim(spc_name_wrf(ispc_wrf)) == trim(spc_name(ispc))) then
         if(trim(spc_name_wrf(ispc_wrf)) == 'ASH' .or. &
            trim(spc_name_wrf(ispc_wrf)) == 'SO2'      ) then
	 dummy=0.
	 print*,' specie=',trim(spc_name(ispc))
	 
         do j=1,n3 ; do i=1,n2
          if( volc_prop(i,j,2) > 0.001) &! time duration
	  
          dummy(i,j) = var2d(i,j,ispc,geoge)*1.e+9  / volc_prop(i,j,3) !  ug/(m^2 sec)
         
         enddo;enddo
	 write(91) dummy
       endif;endif
     enddo;enddo
     ! inj height
     write(91) volc_prop(:,:,1)
     ! elvation
     write(91) volc_prop(:,:,2)
     ! duration
     write(91) volc_prop(:,:,3)

   close (91)
print*,'volcanic SO2 and ASH emission, height, elevation and duration  data for FIM: ',fileBB
print*,'volcanic begin time=',begin_eruption
print*,'Volc ASH size distribution (fraction)'
print*, ash_size_distr(:)

endif
!---------------------------------------------------------------------------

if(use_degass_volcanoes == 1) then
   ename='volc-degass'
   begin_eruption= "NOT-DEFINED"  ! dummy
   fileBB=file_in(1:len_trim(file_in)-4)//'-volc.bin'
   open (91,file=fileBB,form='unformatted')
     write(91) 3 !fields
     write(91)ename
     write(91)begin_eruption
     iv=0
     do ispc_wrf=1,nspecies_wrf
      do ispc = 1,nspecies
        if(trim(spc_name_wrf(ispc_wrf)) == trim(spc_name(ispc))) then
         if(trim(spc_name_wrf(ispc_wrf)) == 'SO2') then
	 dummy=0.
	 print*,' specie=',trim(spc_name(ispc))
	 do j=1,n3 ; do i=1,n2
          if( volc_degass_prop(i,j,2) > 0.001) &! time duration
	  
          dummy(i,j) = var2d(i,j,ispc,geoge) / volc_degass_prop(i,j,2) !  kg/(m^2 sec)
         enddo;enddo
	 write(91) dummy
       endif;endif
     enddo;enddo
     ! inj height
     write(91) volc_degass_prop(:,:,1)
     ! elevation
     write(91) volc_degass_prop(:,:,2)
   close (91)
print*,'volcanic degassing SO2 emission, height and elevation data for wrf: ',fileBB
endif
!---------------------------------------------------------------------------

return  !<<<<<<<<<<
!----------- GOCART background data --------------
if(use_gocart_bg == 1) then

fileBG=file_in(1:len_trim(file_in)-4)//'-gocartBG.bin'
open (91,file=fileBG,form='unformatted')

 do ispc=1,gocart_bg_nspecies
   print*,'writing spc=',gocart_bg_spc_name(ispc),n2,n3,nlevels_netcdf(ispc)

   !-DMS
   if(trim(gocart_bg_spc_name(ispc))=='DMS') then
       print*,'DMS for month=',imon
       print*,'gocart_bg_spc_name',gocart_bg_spc_name(ispc)
       write(91) gocart_bg_g(ispc)%src(1:n2,1:n3,imon)
       cycle
   endif

   !special section of writing 'lev'
   if(trim(gocart_bg_spc_name(ispc))=='H2O2') then
       write(91) lev(1:nlevels_netcdf(ispc))
   endif


   write(91) gocart_bg_g(ispc)%src(1:n2,1:n3,1:nlevels_netcdf(ispc))

 enddo
close (91)

endif

!    use gocart_backgr, only :  gocart_bg_nspecies=>nspecies&
!                              ,gocart_bg_spc_name=>spc_name&
!			      ,nlevels_netcdf,gocart_bg_g,lev!

end subroutine write_to_FIM
!---------------------------------------------------------
subroutine write_to_WRFCHEM(n2,n3,nspecies,nsrc,var2d,plume_mean_fct &
                             ,plume_firesize,volc_prop,file_in,iday,imon,iyear&
			     ,volc_degass_prop,ash_size_distr)
    use grid_dims_out
    use chem1_list, only : spc_name, weight, BBURN2,BBURN3,OC,BC, URBAN2,URBAN3
    use gocart_backgr, only :  gocart_bg_nspecies=>nspecies&
                              ,gocart_bg_spc_name=>spc_name&
			      ,nlevels_netcdf,gocart_bg_g,lev
    implicit none
    integer,intent(IN) :: n2,n3,nspecies,nsrc
    real,dimension(n2,n3) :: smold_frac
    real,dimension(n2,n3,4) :: plume_mean_fct,plume_firesize
    real,dimension(n2,n3,3) :: volc_prop,volc_degass_prop
    real,dimension(n2,n3,nspecies,nsrc) :: var2d
    real, dimension(n2,n3) :: dummy
    character(*) ::file_in
    character(len=240) :: fileAB,fileBB,fileBG
    integer     :: &
      antro   = 01, & ! anthropogenic sources
      bburn   = 02, & ! biomass burning sources
      bioge   = 03 ,&   ! biogenic sources
      geoge   = 04    ! biogenic sources
    integer, parameter :: nveg_agreg      = 4
    integer, parameter :: tropical_forest = 1
    integer, parameter :: boreal_forest   = 2
    integer, parameter :: savannah        = 3
    integer, parameter :: grassland       = 4

    integer, parameter :: nspecies_wrf = 29 

    integer ispc,isrc,ispc_wrf,iv,iveg,ifound,iday,imon,iyear,i,j
    character(len=20) ename
    real ash_size_distr(10)
    integer itime
    real fx
!from     /home/poluicao/WRFV2-CHEM/chem/module_input_chem_data.F
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspecies_wrf) :: spc_name_wrf=(/ &
 'SO2     ' &!
,'NO2     ' &!
,'NO      ' &!
,'ALD     ' &!
,'HCHO    ' &!
,'ORA2    ' &!
,'NH3     ' &! XXXXX
,'HC3     ' &!
,'HC5     ' &!
,'HC8     ' &!
,'ETH     ' &!
,'CO      ' &!
,'OL2     ' &! XXXXX ! NUWRF...Not read by convert_emiss 3.4; see below.
,'OLT     ' &!
,'OLI     ' &!
,'TOL     ' &!
,'XYL     ' &!
,'KET     ' &!
,'CSL     ' &!
,'ISO     ' &!
,'BBURN2  ' &!
,'BBURN3  ' &!)
,'URBAN2  '&
,'URBAN3  '&
,'OC      ' &!
,'BC      ' &
,'DMS     ' &
,'SO4     ' &! SO4
,'ASH     ' &

/)

  ! EMK NUWRF...Extra emissions for WRF-Chem not included in community code.
  integer,parameter :: nspecies_wrf_extra = 7
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspecies_wrf_extra) :: &
       spc_name_wrf_extra=(/ &
 'ETE     ' &
,'CH4     ' &
,'CO2     ' &
,'API     ' &
,'LIM     ' &
,'ORA1    ' &
,'DM      ' &
/)

!
!,'PM25I ' &!
!,'PM25J ' &!
!,'SO4I	' &!
!,'SO4J	' &!
!,'NO3I	' &!
!,'NO3J	' &!
!,'ORGI	' &!
!,'ORGJ	' &!
!,'ECI	' &!
!,'ECJ	' &!
!,'PM10	'

!antro + bioge section-------------------------------------
print*,'=> anthro + bioge section ---------------------------'

ename='XXXXXXXX'
itime = 0
fileAB=file_in(1:len_trim(file_in)-4)//'-ab.bin'

open (91,file=fileAB,form='unformatted')
     write(91)nspecies_wrf
     write(91)ename
     write(91)itime
iv=0
do ispc_wrf=1,nspecies_wrf
  print*,'spc_name_wrf(ispc_wrf)',spc_name_wrf(ispc_wrf)
  ifound=0
  do ispc = 1,nspecies
     !if(spc_alloc(src,ispc) == 0) cycle

     !print*,'YY    ',spc_name(ispc),spc_name_wrf(ispc_wrf)

     if(trim(spc_name_wrf(ispc_wrf)) == trim(spc_name(ispc))) then

     	 if(trim(spc_name_wrf(ispc_wrf)) == 'BBURN2' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'BBURN3') then
	    ifound=1
	    cycle
	 endif
	 
         print*,'FOUND ',spc_name_wrf(ispc_wrf),spc_name(ispc)
	 
         !-for WRF/FIM, 'biogenic' emissions for CO is not included
	 !- 
	 if(trim(spc_name(ispc)) == 'CO') then
	     dummy(:,:) = var2d(:,:,ispc,antro)
	 else
	     dummy(:,:) = var2d(:,:,ispc,antro)+var2d(:,:,ispc,bioge)
	 endif
	    
	 if(trim(spc_name_wrf(ispc_wrf)) == 'URBAN2' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'URBAN3' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'OC'   .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'BC'  ) then
	    fx= 1.e+9  / 86400.	 !  ug/(m^2 sec)
	 else
	    fx= 1.e+3 * 1.e6 / weight(ispc) / 24. ! emission unit: mol/(km^2 hour)
	 endif

	 print*, 'max/min values= ', maxval(dummy), minval(dummy),iv+1

         write(91) dummy * fx
         iv=iv+1
	 ifound=1
	 exit
     endif
  enddo
  if(ifound == 0) then
    print*,'NOT FOUND, zero emission',spc_name_wrf(ispc_wrf),iv+1

!     ! NU-WRF...Changed convert_emiss to intentionally read in OL2, so
      ! the if-block below is unnecessary.
!     if (trim(spc_name_wrf(ispc_wrf)) == 'OL2') then
!        ifound=1
!        cycle
!     end if

     dummy(:,:) =0.
     write(91) dummy
     iv=iv+1
     ifound=1
  endif

enddo
close (91)
print*,'Anthro + biogenic emission for wrf: ',fileAB


!--------------------------       bburn + plumerise section
print*,'=> bburn + plumerise section ---------------------------'

fileBB=file_in(1:len_trim(file_in)-4)//'-bb.bin'
open (91,file=fileBB,form='unformatted')
! NUWRF EMK...Add extra species
!     write(91)nspecies_wrf + 8
     write(91)nspecies_wrf + 8 + nspecies_wrf_extra
     write(91)ename
     write(91)itime
iv=0
do ispc_wrf=1,nspecies_wrf

   ifound=0
   do ispc = 1,nspecies
     if(trim(spc_name_wrf(ispc_wrf)) == trim(spc_name(ispc))) then
     	 if(trim(spc_name_wrf(ispc_wrf)) == 'URBAN2' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'URBAN3') then
	    ifound=1
	    cycle
	 endif
         print*,' specie=',trim(spc_name(ispc))
         dummy(:,:) = var2d(:,:,ispc,bburn)


	 !-PM25 for wrf:

	 if(trim(spc_name_wrf(ispc_wrf)) == 'BBURN2' ) then

	 dummy(:,:) = max(0.,var2d(:,:,BBURN2,bburn)  -  &
	                     var2d(:,:,OC  ,bburn)  -  &
	                     var2d(:,:,BC  ,bburn)     )

	 endif

	 if(trim(spc_name_wrf(ispc_wrf)) == 'BBURN2' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'BBURN3' .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'OC'   .or. &
	    trim(spc_name_wrf(ispc_wrf)) == 'BC'  ) then
	    fx= 1.e+9  / 86400.	 !  ug/(m^2 sec)
	 else
	    fx= 1.e+3 * 1.e6 / weight(ispc) / 24.! emission unit: mol/(km^2 hour)
	 endif
	 
	 print*, 'BB = ', maxval(dummy), minval(dummy)

	 write(91) dummy*fx
	 iv=iv+1
	 ifound=1
	 !print*,'spc=',iv,trim(spc_name(ispc))
        exit
     endif

  enddo


  if(ifound == 0) then
     print*,'not found, zero emission ',spc_name_wrf(ispc_wrf)
     dummy(:,:) =0.
     write(91) dummy
     iv=iv+1
     ifound=1
  endif

enddo
do iveg=1,nveg_agreg
   write(91) plume_mean_fct(:,:,iveg)
	 iv=iv+1
	 !print*,'iveg plume_fct', iveg,maxval(plume_mean_fct(:,:,iveg)),iv
enddo
do iveg=1,nveg_agreg
   write(91) plume_firesize(:,:,iveg)
	 iv=iv+1
	 !print*,'iveg plume_firesize', iveg,maxval(plume_firesize(:,:,iveg)),iv
enddo

! NUWRF EMK...Add extra GFEDV3 species
do ispc_wrf=1,nspecies_wrf_extra
   ifound=0
   do ispc = 1,nspecies
      if(trim(spc_name_wrf_extra(ispc_wrf)) == trim(spc_name(ispc))) then
         print*,' species=',trim(spc_name(ispc))
         dummy(:,:) = var2d(:,:,ispc,bburn)
         ! NUWRF EMK...Change DM units from kg/m^2/day to kg/km^2/hour. We don't
         ! have a molecular weight for this category.
         if (trim(spc_name_wrf_extra(ispc_wrf)) == 'DM') then
            fx = 1.e6 / 24.
         else
            fx= 1.e+3 * 1.e6 / weight(ispc) / 24.! emission unit: mol/(km^2 hour)
         end if
	 print*, 'BB = ', maxval(dummy), minval(dummy)
         write(91) dummy*fx
         iv=iv+1
         ifound=1
         exit
      endif
   enddo
   if(ifound == 0) then
      print*,'not found, zero emission ',spc_name_wrf_extra(ispc_wrf)
      dummy(:,:) =0.
      write(91) dummy
      iv=iv+1
      ifound=1
   endif
enddo
! NUWRF EMK END...Add extra GFEDV3 species

close (91)
!print*,'BB: total 2d fields=',iv
print*,' Biomass burning emission + plumerise data for wrf: ',fileBB

!-------------------------- volcanoes -------------------------
print*,'=> ------------------------------------------------------'
print*,'=> Volcanoes Eruption Section ---------------------------'

if(use_volcanoes == 1) then
   ename='volc-eruption'
   fileBB=file_in(1:len_trim(file_in)-4)//'-volc.bin'
   open (91,file=fileBB,form='unformatted')
     write(91) 5 !fields
     write(91)ename
     write(91) trim(begin_eruption)
     write(91) ash_size_distr(:)
     print*,'ash_size_distr=',ash_size_distr(:),trim(begin_eruption)

     iv=0
     do ispc_wrf=1,nspecies_wrf
      do ispc = 1,nspecies
        if(trim(spc_name_wrf(ispc_wrf)) == trim(spc_name(ispc))) then
         if(trim(spc_name_wrf(ispc_wrf)) == 'ASH' .or. &
            trim(spc_name_wrf(ispc_wrf)) == 'SO2'      ) then
	 dummy=0.
	 print*,' specie=',trim(spc_name(ispc))
	 
         do j=1,n3 ; do i=1,n2
          if( volc_prop(i,j,3) > 0.001) then ! time duration

	   if( trim(spc_name_wrf(ispc_wrf)) == 'SO2') then
            fx= 1.e+3 * 1.e6 / weight(ispc) *3600.
            dummy(i,j) = var2d(i,j,ispc,geoge)*fx  / volc_prop(i,j,3) !  mol/(km^2 hour)

           else

            dummy(i,j) = var2d(i,j,ispc,geoge)*1.e+9  / volc_prop(i,j,3) !  ug/(m^2 sec)
           endif
          endif

         enddo;enddo
	 write(91) dummy
	 print*,'max = ',maxval(dummy), trim(spc_name_wrf(ispc_wrf))
       endif;endif
     enddo;enddo
     ! inj height
     write(91) volc_prop(:,:,1)
     ! elvation
     write(91) volc_prop(:,:,2)
     ! duration
     write(91) volc_prop(:,:,3)

   close (91)
print*,'volcanic SO2 and ASH emission, height,elevation and duration   data for WRF: ',fileBB
print*,'volcanic begin time=',begin_eruption
print*,'Volc ASH size distribution (fraction)'
print*, ash_size_distr(:)

endif
!---------------------------------------------------------------------------

if(use_degass_volcanoes == 1) then
   ename='volc-degass'
   begin_eruption= "NOT-DEFINED"  ! dummy
   fileBB=file_in(1:len_trim(file_in)-4)//'-volc.bin'
   open (91,file=fileBB,form='unformatted')
     write(91) 3 !fields
     write(91)ename
     write(91)begin_eruption
     iv=0
     do ispc_wrf=1,nspecies_wrf
      do ispc = 1,nspecies
        if(trim(spc_name_wrf(ispc_wrf)) == trim(spc_name(ispc))) then
         if(trim(spc_name_wrf(ispc_wrf)) == 'SO2') then
	 dummy=0.
	 print*,' specie=',trim(spc_name(ispc))
	 do j=1,n3 ; do i=1,n2
          if( volc_degass_prop(i,j,2) > 0.001) &! time duration
	  
          dummy(i,j) = var2d(i,j,ispc,geoge) / volc_degass_prop(i,j,2) !  kg/(m^2 sec)
         enddo;enddo
	 write(91) dummy
       endif;endif
     enddo;enddo
     ! inj height
     write(91) volc_degass_prop(:,:,1)
     ! elevation
     write(91) volc_degass_prop(:,:,2)
   close (91)
print*,'volcanic degassing SO2 emission, height and elevation data for wrf: ',fileBB
endif
!---------------------------------------------------------------------------

!print*,' The files *.ctl, *.gra and *.vfm are only auxiliary files'
!print*,' and, actually, are not used by WRF model.'


!----------- GOCART background data --------------
if(use_gocart_bg == 1) then

print*,'=> ------------------------------------------------------'
print*,'=> Gocart Background  Section ---------------------------'


fileBG=file_in(1:len_trim(file_in)-4)//'-gocartBG.bin'
open (91,file=fileBG,form='unformatted')
 do ispc=1,gocart_bg_nspecies
   print*,'writing spc = ', trim(gocart_bg_spc_name(ispc)),n2,n3,nlevels_netcdf(ispc)

   !-DMS
   if(trim(gocart_bg_spc_name(ispc))=='DMS') then
       print*,'DMS for month = ',imon
       print*,'gocart_bg_spc_name ', trim(gocart_bg_spc_name(ispc))
       write(91) gocart_bg_g(ispc)%src(1:n2,1:n3,imon)
       cycle
   endif

   !special section of writing 'lev'
   if(trim(gocart_bg_spc_name(ispc))=='H2O2') then
       write(91) lev(1:nlevels_netcdf(ispc))
   endif


   write(91) gocart_bg_g(ispc)%src(1:n2,1:n3,1:nlevels_netcdf(ispc))

 enddo
close (91)

endif

!    use gocart_backgr, only :  gocart_bg_nspecies=>nspecies&
!                              ,gocart_bg_spc_name=>spc_name&
!			      ,nlevels_netcdf,gocart_bg_g,lev!

end subroutine write_to_WRFCHEM
!---------------------------------------------------------

#endif

