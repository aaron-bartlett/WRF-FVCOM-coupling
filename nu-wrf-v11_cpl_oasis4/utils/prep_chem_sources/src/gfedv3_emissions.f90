!#############################################################################
!
!  CCATT-BRAMS/MCGA-CPTEC/WRF emission model    CPTEC/INPE                   
!
!  Version 1.0.0: 12/nov/2010                                                 
!
!  Coded by Saulo Freitas and Karla Longo                                     
!
!  Contact: gmai@cptec.inpe.br - http://meioambiente.cptec.inpe.br            
!
!#############################################################################


module gfedv3_emissions    

   use AeM_emission_factors,  nspecies=>AeM_nspecies, N2_nitrogenio=>N2

   type gfedv3_vars   
      real, pointer, dimension(:,:,:)  :: src
   end type gfedv3_vars

   type (gfedv3_vars), allocatable :: gfedv3_g(:)

contains

   subroutine alloc_gfedv3(gfedv3,n1,n2,n3,nspecies)

      implicit none

      type (gfedv3_vars),dimension(nspecies)  :: gfedv3
      integer,intent(in) :: n1,n2,n3
      integer ispc,nspecies

      do ispc=1,nspecies
         allocate (gfedv3(ispc)%src(n1,n2,n3))
      enddo

      return
   end subroutine alloc_gfedv3

   subroutine nullify_gfedv3(gfedv3,nspecies)

      implicit none

      type (gfedv3_vars),dimension(nspecies)  :: gfedv3
      integer ispc,nspecies

      do ispc=1,nspecies
         if (associated(gfedv3(ispc)%src))    nullify (gfedv3(ispc)%src)
      enddo

      return
   end subroutine nullify_gfedv3

end module gfedv3_emissions

subroutine mem_gfedv3(n1,n2,n3)
   use gfedv3_emissions
   implicit none
   integer i
   integer, intent(in) :: n1,n2,n3

   if(.not. allocated(gfedv3_g)) allocate(gfedv3_g(nspecies))

   call nullify_gfedv3(gfedv3_g(:),nspecies)      
   call alloc_gfedv3  (gfedv3_g(:),n1,n2,n3,nspecies) 
end subroutine mem_gfedv3

subroutine read_gfedv3(iyear,imon,iday,gfedv3_suffix,ng,ngrids,n1,n2,n3, &
     rlat,rlon,rland,deltax,deltay,xt,yt,plat,plon)

   use grid_dims_out, only: maxfiles,grid_type,gfedv3_data_dir, &
        use_bbem_plumerise
   use gfedv3_emissions
   use bbbem_plumerise

   implicit none

   integer, parameter ::  nlon = 720, nlat=360, nmonths=12 
   integer, intent (in) :: iyear,imon,iday,ng,n1,n2,n3,ngrids
   real, intent (in), dimension(n1,n2)  :: rlat,rlon,rland
   real, intent (in) :: deltax,deltay
   integer :: ispc,im,i,k,iread,nc,ilaea,jlaea,kk,k1,k2,ii,i1,i2,j,igbox &
        ,jgbox,j1,j2,ic,jc,iveg_ag, icount

   character cjday*3,cyear*4,cmon*2,cday*2
   character*240 filename, dayFractionFile
   character(len=*) :: gfedv3_suffix   ! "comma" delimitted species names
   logical there
   integer jday,jday_avail
   integer iveg, veg_type(nlon,nlat)
   integer  itimes(maxfiles),nfiles
   real longgfedv3(nlon),latgfedv3(nlat)
   real, parameter :: ilatn=.5 , ilonn=.5  
   real rrlat,rrlon,dlon1,dlon2,dlat1,dlat2,TX(nspecies),rdummy
   real dx,dy,  dxm ,dym, xt(n1), yt(n2),plat,plon  
   real, dimension(12) :: mondays
   data mondays/31,28,31,30,31,30,31,31,30,31,30,31/
   real, allocatable, save, dimension(:,:) :: carb,veg,carb_mean, &
        fractionOfEmissions
   real, allocatable, save, dimension(:,:,:) :: drybbc
   real, save :: typical_fire_size(0:3)
   character(len=20), allocatable :: consideredSpecies(:)
   integer :: numConsideredSpecies,iConsideredSpecies,ident,ident2

   call getNumConsideredSpecies(gfedv3_suffix,numConsideredSpecies)
   allocate(consideredSpecies(numConsideredSpecies))
   call getConsideredSpecies(gfedv3_suffix,numConsideredSpecies, &
        consideredSpecies)

   if(.not. allocated(carb)) allocate(carb(nlon,nlat), &
        drybbc(nspecies,nlon,nlat),veg(nlon,nlat),carb_mean(nlon,nlat), &
        fractionOfEmissions(nlon,nlat))

   !- allocate memory for plumerise arrays
   if(use_bbem_plumerise == 1) call mem_bbbem_plume(n1,n2,n3)

   !-lat e lon gfedv3 (corner mais ao sul e mais a oeste) 
   do i=1,nlon
      longgfedv3(i)=-179.75 + (i-1)*ilonn
   enddo
   do j=1,nlat
      latgfedv3(j)= -89.75 + (j-1)*ilatn
   enddo

   !- avoid negative values for emission factor
   where(emission_factor <0.) emission_factor=0. 

   if(ng == 1) then  

      print*,'==================   GFEDv3  ==================================='
      !- get vegetation dataset 
      veg = 0.0 ! not used for gfedv3; should clean

      !- determine if the day fraction data is available
      call dataset_avail_dayFraction(iyear,imon,iday,dayFractionFile)
      !
      dayFractionFile=trim(gfedv3_data_dir)//'/'//trim(dayFractionFile)
      !print*,'day fraction file= ',trim(dayFractionFile)
      inquire(file=dayFractionFile,exist=there)
      if(.not.there) then	
         !- calculate the time average over the existent dataset for the same 
         !  8-days period
         print*, 'FILE= ',trim(dayFractionFile),' does not exist, will stop'
         stop "stopping...."
      else
         call readDailyFractionGFEDv3(dayFractionFile,fractionOfEmissions, &
              nlon,nlat)
      endif

      drybbc = 0.0
      do iConsideredSpecies = 1, numConsideredSpecies
         !- determine the nearest monthly data available
         call dataset_avail_gfedv3(iyear,imon, &
              consideredSpecies(iConsideredSpecies),filename)
         !
         !- reading Carbon emission dataset 
         filename=trim(gfedv3_data_dir)//'/'//trim(filename)
         !print*,'file= ',trim(filename)
         inquire(file=filename,exist=there)
         if(.not.there) then	
            !- calculate the time average over the existent dataset for the 
            ! same 8-days period
            print*, 'FILE= ',trim(filename),' does not exist, will stop'
            stop "stopping...."
         else
            !print *,'GFEDv3 source opening   ', trim(filename)
            open(30, FILE=filename, STATUS='old')
            do j=1,360
               read(30,'(720e16.7)') carb(:,j)	
            enddo
            close(30)
         endif

         !convert from carbon flux to dry-biomass-consumed per day
         !  drybbc=carb/(0.45*8.)  ! unit = g[DM]/m^2/day         
         ! Convert from g/m^2/month to g/m^2/day using daily fraction of
         ! monthly emissions.
         call get_gfedv3_identity(consideredSpecies(iConsideredSpecies), &
              ident,ident2)
         ! Several GFEDV3 species are divided equally into two AeM species.
         if (ident2 > 0 .and. ident > 0) then
            drybbc(ident,:,:)=0.5*fractionOfEmissions*carb ! unit = g[DM]/m^2/day
            drybbc(ident2,:,:)=0.5*fractionOfEmissions*carb ! unit = g[DM]/m^2/day                                 
         ! EMK BUG FIX: Make sure valid index returned
         else if (ident > 0) then
            drybbc(ident,:,:)=fractionOfEmissions*carb ! unit = g[DM]/m^2/day
         end if
      enddo

      ! compatibilize a legenda de vegetacao com a adoptada no veg_agreg
      !gfedv3: 1=sav, 2=trop for, 3=extra trop forest
      ! agreg: 1 =trop forest, 2 =extra trop for, 3 =sav, 4 = past, 5=deserto,
      ! 0 = agua/gelo
      ! veg dominante por grid de 1x1 grau (data from gfedv3)
      where(nint(veg(:,:)) == 1) veg_type(:,:)=3
      where(nint(veg(:,:)) == 2) veg_type(:,:)=1
      where(nint(veg(:,:)) == 3) veg_type(:,:)=2

      !- for plumerise - typical fire size
      typical_fire_size(0)=0.
      typical_fire_size(1)=20.* 1.e4 ! 20ha for tropical forest
      typical_fire_size(2)=20.* 1.e4 ! 20ha for extra tropical
      typical_fire_size(3)=5. * 1.e4 ! 5ha for savana/past ...

      !-revert y-dir
      call yrevert(nlon,nlat,1,veg_type)
      do iConsideredSpecies = 1, numConsideredSpecies
         call get_gfedv3_identity(consideredSpecies(iConsideredSpecies), &
              ident,ident2)
         ! EMK NUWRF...Some GFED species are split into two AeM species
         if (ident2 > 0 .and. ident > 0) then
            call yrevert(nlon,nlat,1,drybbc(ident,:,:))
            call yrevert(nlon,nlat,1,drybbc(ident2,:,:))
         ! EMK NUWRF...Bug fix
         else if (ident > 0) then
            call yrevert(nlon,nlat,1,drybbc(ident,:,:))
         end if
      enddo

   endif !(ng==1) 


   !--- performs the interpolation to model grid box

   do i=1,n1
      do j=1,n2

         call get_index1(i,j,nlon,nlat,n1,n2,rlat,rlon,longgfedv3,latgfedv3 &
              ,ilatn, ilonn ,dlat1,dlat2,dlon1,dlon2,i1,j1,i2,j2,ic,jc)
         !print*,'ic jc=',ic,jc,rlat(i,j),rlon(i,j)
         if(ic == 0 .or. jc == 0. .or. ic > 720 .or. jc > 360) cycle
         ! RAWsrc contains CO data 
         call interpol_gfedv3(i,j,n1,n2,rlon,rlat,ic,jc,nlat,nlon, &
              ilatn,ilonn,nspecies,drybbc,veg_type,tx)
         ! 
         ! TX is the value interpolated to the model grid box.
         ! Convert from g/m^2/day to kg/m^2/day
         do ispc=1,nspecies
            gfedv3_g(ispc)%src(i,j,1)=TX(ispc) * 1.e-3 
         enddo

         !------- plumerise section - only for areas with burning
         if(use_bbem_plumerise /= 1 .or. veg_type(ic,jc) == 0 &
              .or. gfedv3_g(1)%src(i,j,1) < 1.e-6) cycle

         iveg_ag = veg_type(ic,jc) 

         do ispc=1,nspecies
            bbbem_plume_g(ispc,iveg_ag)%src(i,j,1) = &
                 gfedv3_g(ispc)%src(i,j,1)*flaming(iveg_ag)
         enddo

         !-fires properties : accumulated size
         bbbem_plume_fire_prop_g(qarea_agreg,iveg_ag)%fire_prop(i,j) = &
              typical_fire_size(veg_type(ic,jc))

         !-fires properties : fire number
         bbbem_plume_fire_prop_g(qfires_agreg,iveg_ag)%fire_prop(i,j) = 1.     
         !------- end of plumerise section

      enddo
   enddo

   !------- plumerise section 2
   if(use_bbem_plumerise == 1) then
      !- produce some statistical quantities for the plumerise model
      !
      do i=1,n1
         do j=1,n2
            do iveg_ag=1,nveg_agreg
               ! calculate the mean size of aggregate fire

               if( bbbem_plume_fire_prop_g(qfires_agreg,iveg_ag) % &
                    fire_prop(i,j) .ge. 1.) then 

                  bbbem_plume_fire_prop_g(qarea_agreg ,iveg_ag) % &
                       fire_prop(i,j) = & ! mean  size =
                       bbbem_plume_fire_prop_g(qarea_agreg ,iveg_ag) % &
                       fire_prop(i,j) / & ! total size / nfires
                       bbbem_plume_fire_prop_g(qfires_agreg,iveg_ag) % &
                       fire_prop(i,j)
               else   
                  bbbem_plume_fire_prop_g(qarea_agreg ,iveg_ag) % &
                       fire_prop(i,j) = 0.
               endif
               !- convert from mass consumed during the flaming phase to 
               !  fraction of the total mass consumed


               do ispc=1,nspecies
                  if(gfedv3_g(ispc)%src(i,j,1) > 0.) &	  
                                !- fraction consumed at the flaming phase:
                       bbbem_plume_g(ispc,iveg_ag)%src(i,j,1) = &
                       bbbem_plume_g(ispc,iveg_ag)%src(i,j,1) &
                       / gfedv3_g(ispc)%src(i,j,1)
               enddo
            enddo
         enddo
      enddo

      !- calculate the mean value for bbbem_plume_g (to save memory)
      do i=1,n1
         do j=1,n2
            do iveg_ag=1,nveg_agreg
               icount = 0
               bbbem_plume_mean_g(iveg_ag)%src(i,j,1)=0.

               do ispc=1,nspecies
                  if(gfedv3_g(ispc)%src(i,j,1) > 0.) then
                     icount = icount + 1
                     bbbem_plume_mean_g(iveg_ag)%src(i,j,1) = &
                          bbbem_plume_mean_g(  iveg_ag)%src(i,j,1)+ &
                          bbbem_plume_g(ispc,iveg_ag)%src(i,j,1)
                  endif
               enddo
               bbbem_plume_mean_g(iveg_ag)%src(i,j,1) = &
                    bbbem_plume_mean_g(iveg_ag)%src(i,j,1) / &
                    (float(icount)+1.e-10)
            enddo
         enddo
      enddo

   endif
   !------- end of plumerise section 2

   !print*,'gfedv3_g(CO)%src(i,j,1)=',maxval(gfedv3_g(CO)%src(:,:,1))
   !pause '1 '
   !return
   if(grid_type == 'rams' .or. grid_type == 'polar') then 
      ! only for 'rams' until rland is also defined for others grids 
      do ispc=1,nspecies 
         call apply_land_restriction(n1,n2,rland,gfedv3_g(ispc)%src(:,:,1))
      enddo
   endif
   !- deallocate memory that we do not need anymore    
   if(ng==ngrids) deallocate (carb,drybbc,veg,carb_mean,fractionOfEmissions)

end subroutine read_gfedv3

subroutine get_gfedv3_identity(spc_name,ident,ident2)
   !use chem1_list
   use gfedv3_emissions, only : gfedv3_nspecies=>nspecies& 
        !don't use AeM_nspecies
        ,gfedv3_spc_name=>AeM_spc_name
   implicit none
   integer isp
   character (len=*), intent(in)  :: spc_name
   integer                        :: ident
   integer                        :: ident2
   integer :: ii

   do isp = 1,gfedv3_nspecies
      ident=-1
      ident2=-1
      if(spc_name == gfedv3_spc_name(isp)) then
         print*,'Matched GFEDV3 species with AeM species ',trim(spc_name)
         ident=isp
         return
      endif
   end do
   ! EMK NUWRF Bug fix...Not all GFEDV3 names exactly match the AeM list.
   ! So some special double checks are needed.
   ! We'll have extensive checking to better detect breakage if the AeM
   ! list changes in the future.
   ! TODO: Handle C4, C, DM, Higher_Alkanes, Higher_Alkenes, 

   if (trim(spc_name) == 'C2H4O') then
      ident = special_identity('Acetald',gfedv3_spc_name, &
           gfedv3_nspecies)         
   else if (trim(spc_name) == 'C2H5OH') then
      ident = special_identity('Ethanol',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'C2H6S') then
      ident = special_identity('DMS',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'C3H6O') then
      ident  = special_identity('Propanal',gfedv3_spc_name, &
           gfedv3_nspecies)
      ident2 = special_identity('Acetone',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'C5H8') then
      ident  = special_identity('cyclopentene',gfedv3_spc_name, &
           gfedv3_nspecies)
      ident2 = special_identity('Isoprene',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'CH2O') then
      ident = special_identity('Formaldehyde',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'CH3OH') then
      ident = special_identity('Methanol',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'NMHC') then
      ident = special_identity('NHMC',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'PM2p5') then
      ident = special_identity('BBURN2',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'Terpenes') then
      ident = special_identity('terpenes',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'Toluene_lump') then
      ident = special_identity('toluene',gfedv3_spc_name, &
           gfedv3_nspecies)
   else if (trim(spc_name) == 'TPM') then
      ident = special_identity('BBURN3',gfedv3_spc_name, &
           gfedv3_nspecies)
   end if
   if (ident < 0) then
      print*,'Warning, cannot process GFEDV3 species ',trim(spc_name)
   end if

   
contains
   ! EMK NUWRF
   ! Internal function for special comparisons with AeM species list.
   ! Assumes a match will be found with AeM. A fatal error will occur
   ! if this is not the case -- useful for detecting future changes to the
   ! AeM species list.
   function special_identity(spc_name,gfedv3_spc_name, &
        gfedv3_nspecies) result(ident)
      implicit none
      character(len=*),intent(in) :: spc_name
      integer,intent(in) :: gfedv3_nspecies
      character(len=*),intent(in) :: gfedv3_spc_name(gfedv3_nspecies)
      integer :: ident
      integer :: ii
      ident=-1
      do ii = 1,gfedv3_nspecies
         if (trim(spc_name) == trim(gfedv3_spc_name(ii))) then
            print*,'Matched GFEDV3 species with AeM species ',trim(spc_name)
            ident=ii
            return               
         end if
      end do
      print*,'INTERNAL ERROR, cannot match GFEDV3 species with AeM species', &
           trim(spc_name)
      stop 
   end function special_identity

end subroutine get_gfedv3_identity

subroutine interpol_gfedv3(i,j,n1,n2,rlon,rlat,ic,jc,nlat,nlon,ilatn,ilonn &
     ,nspecies,drybbc,veg_type &
     ,tx)
   use AeM_emission_factors, N2_nitrogenio=>N2
   use grid_dims_out, only: grid_type
   implicit none
   integer n1,n2,ic,jc,nlat,nlon,i,j,nspecies,ispc
   real, dimension(n1,n2) :: rlat,rlon
   real ilatn,ilonn,tx(nspecies)
   real, dimension(nspecies,nlon,nlat) :: drybbc
   integer  veg_type(nlon,nlat)

   !-local var
   real dlonr,dlatr,usdum
   integer qi1,qi2,qj1,qj2,ncount,ii,jj
   integer iveg , veg_type_model(n1,n2)

   !if(ic.ge.0 .and. jc .ge. 0) then 

   if(grid_type == 'fim') then

      ncount = 0
      dlonr = 0.
      dlatr = 0.
      do ii=1,n1-1
         ncount = ncount+1
         dlonr= dlonr + rlon(ii+1,1)-rlon(ii,1)
         dlatr= dlatr + rlat(ii+1,1)-rlat(ii,1)
      enddo
      dlonr = 0.5*dlonr/(float(ncount) + 1.E-10)
      dlatr = 0.5*dlatr/(float(ncount) + 1.E-10)
   else    	     

      dlonr=0.5*(rlon(n1,j)-rlon(1,j))/float(n1-1)
      dlatr=0.5*(rlat(i,n2)-rlat(i,1))/float(n2-1)
   endif
   qi1=int(dlonr/ilonn+0.5)
   qi2=int(dlonr/ilonn+0.5)
   qj1=int(dlatr/ilatn+0.5)
   qj2=int(dlatr/ilatn+0.5)

   ncount = 0
   TX(:)  = 0.

   do jj =max(1,jc-qj1),min(nlat,jc+qj2)
      do ii = max(1,ic-qi1),min(nlon,ic+qi2)

         !     if(veg_type(ii,jj) > 0)  then 

         ncount = ncount + 1

         ! use emission factors to convert from  others species ...
         !       TX(:) = TX(:) + drybbc(ii,jj) * & 
         ! Em factor is in g/kg (1.e-3 => kg/kg]
         !               emission_factor(veg_type(ii,jj),: )*1.e-3
         ! NUWRF EMK TODO...Add option here to adjust emissions factors using
         ! DM. But for now, just use the emissions calculated by the GFEDV3
         ! team.
         TX(:) = TX(:) + drybbc(:,ii,jj)
         !      endif 
      enddo
   enddo
   TX(:) = TX(:) / (float(ncount) + 1.E-10) ! interpolated rate
end subroutine interpol_gfedv3

subroutine getNumConsideredSpecies(gfedv3_suffix, numConsideredSpecies)
   character(len=*), intent(in) :: gfedv3_suffix
   integer :: numConsideredSpecies

   character(len=len_trim(gfedv3_suffix)) :: local_gfedv3_suffix
   integer :: localPosition,localLength

   numConsideredSpecies = 0
   local_gfedv3_suffix = trim(gfedv3_suffix)
   do
      localPosition = scan(trim(local_gfedv3_suffix),",")
      localLength = len_trim(local_gfedv3_suffix)

      if (localPosition > 0) then
         numConsideredSpecies = numConsideredSpecies + 1
         local_gfedv3_suffix = &
              local_gfedv3_suffix(localPosition+1:localLength) 
         localPostion = 0
      else
         exit
      endif
   enddo

   numConsideredSpecies = numConsideredSpecies+1

end subroutine getNumConsideredSpecies
!---------------------------------------------------------------

subroutine getConsideredSpecies(gfedv3_suffix, &
     numConsideredSpecies,consideredSpecies)
   character(len=*), intent(in) :: gfedv3_suffix
   integer :: numConsideredSpecies
   character(len=*) :: consideredSpecies(numConsideredSpecies)

   character(len=len_trim(gfedv3_suffix)) :: local_gfedv3_suffix
   integer :: localPosition,localLength

   local_gfedv3_suffix = trim(gfedv3_suffix)
   do iconsideredSpecies = 1, numConsideredSpecies - 1

      localPosition = scan(trim(local_gfedv3_suffix),",")
      localLength = len_trim(local_gfedv3_suffix)
      if (localPosition > 0) then
         consideredSpecies(iconsideredSpecies) = &
              local_gfedv3_suffix(1:localPosition-1)
         local_gfedv3_suffix = &
              local_gfedv3_suffix(localPosition+1:localLength) 
         localPosition = 0
         print *,'DEBUG: ',trim(consideredSpecies(iconsideredSpecies))
      endif

   enddo

   localPosition = scan(trim(gfedv3_suffix),",",back=.true.)
   consideredSpecies(numConsideredSpecies) = &
        gfedv3_suffix(localPosition+1:len_trim(gfedv3_suffix))

end subroutine getConsideredSpecies

subroutine dataset_avail_gfedv3(iyear,imon,gfedv3_suffix,filename)
   implicit none
   integer :: iyear,imon
   character (len=*)  :: filename,gfedv3_suffix
   character (len=4)  :: cyear
   character (len=2)  :: cmon

   !print*,'date required =',iyear,imon

   write(cyear,'(i4.4)') iyear
   write(cmon,'(i2.2)') imon

   filename='GFED3.1_'//trim(cyear)//trim(cmon)//'_'//trim(gfedv3_suffix)// &
        '.txt'
end subroutine dataset_avail_gfedv3

subroutine dataset_avail_dayFraction(iyear,imon,iday,dayFractionFile)
   implicit none
   integer :: iyear,imon,iday
   character (len=*)  :: dayFractionFile
   character (len=4)  :: cyear
   character (len=2)  :: cmon,cday

   !print*,'date required =',iyear,imon,iday

   write(cyear,'(i4.4)') iyear
   write(cmon,'(i2.2)') imon
   write(cday,'(i2.2)') iday

   dayFractionFile='fraction_emissions'//'_'//trim(cyear)//trim(cmon)// &
        trim(cday)//'.nc'
end subroutine dataset_avail_dayFraction

subroutine readDailyFractionGFEDv3(fractionFile,fractionOfEmissions,nlon,nlat)
   implicit none
   include "netcdf.inc"
   character(len=*), intent(in) :: fractionFile
   integer, intent(in) :: nlon,nlat
   real, intent(out) :: fractionOfEmissions(nlon,nlat)
   ! local
   integer             :: cnt (2), strt(2)
   integer             :: ncid
   integer             :: ierr, varid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Open the netCDF file
   ierr = Nf_Open (trim(fractionFile), NF_NOWRITE, ncid)

   if (ierr /= NF_NOERR) then
      print *, 'Nf_Open, cannot open:  '// Trim(fractionFile)
      stop "Code stopped from Nf_Open"
   end if

   strt(:) = (/1,1/)
   cnt (:) = (/nlon,nlat/)

   ierr = Nf_Inq_Varid (ncid, 'Fraction_of_Emissions', varid)

   if (ierr /= NF_NOERR) then
      print *, 'Nf_Inq_Varid  Fraction_of_Emissions' // ', ' // &
           Nf_Strerror (ierr)
      stop "Code stopped from Nf_Inq_Varid"
   end if

   ierr = Nf_Get_Vara_Real(ncid, varid, strt, cnt, fractionOfEmissions)

   if (ierr /= NF_NOERR) then
      print *, 'Nf_Get_Vara_Real ' // Nf_Strerror (ierr)
      stop "Code stopped from Nf_Get_Vara_Real"
   end if

   ierr = Nf_Close (ncid)

end subroutine readDailyFractionGFEDv3
