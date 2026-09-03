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

module gfedv4_emissions    
   use AeM_emission_factors,  nspecies=>AeM_nspecies, N2_nitrogenio=>N2
   implicit none
   private
   public mem_gfedv4
   public read_gfedv4
   public interpol_gfedv4
   public get_gfedv4_identity
   public gfedv4_g

   type gfedv4_vars   
      real, pointer, dimension(:,:,:)  :: src
   end type gfedv4_vars

   type (gfedv4_vars), allocatable :: gfedv4_g(:)
   character(len=64), parameter :: Iam = " gfedv4_emissions: "

contains

   subroutine mem_gfedv4(n1,n2,n3)
      integer, intent(in) :: n1,n2,n3
      integer ispc

      if(.not. allocated(gfedv4_g)) allocate(gfedv4_g(nspecies))
      do ispc=1,nspecies
         if (associated(gfedv4_g(ispc)%src)) nullify(gfedv4_g(ispc)%src)
      enddo
      do ispc=1,nspecies
         allocate (gfedv4_g(ispc)%src(n1,n2,n3))
      enddo

   end subroutine mem_gfedv4

   subroutine read_gfedv4(iyear,imon,iday,gfedv4_suffix,source_option, &
        ng,ngrids,n1,n2,n3,rlat,rlon,rland)
      ! Read GFEDv4 data.
      ! Emission factors array (species, source) are read froma text file  
      ! Daily emission fraction for 1997-2016 are read from HDF5 files.
      !    Read DM emissions and fractional contribution of each source      
      ! On output, gfedv4_g is updated and used by prep_chem_sources. 
      use grid_dims_out, only: maxfiles, grid_type, gfedv4_data_dir, &
           use_bbem_plumerise
      use bbbem_plumerise, only: bbbem_plume_g, bbbem_plume_fire_prop_g, nveg_agreg, &
           qfires_agreg, qarea_agreg, bbbem_plume_mean_g, flaming
      use ReadH5dataset

      integer, intent (in) :: iyear,imon,iday,ng,n1,n2,n3,ngrids
      real, intent (in), dimension(n1,n2)  :: rlat,rlon,rland
      character(len=*) :: gfedv4_suffix   ! "comma" delimitted species names
      integer, intent(in) :: source_option

      integer, parameter :: nlon = 1440, nlat=720
      integer, parameter :: n_species = 41, n_sources=6
      character(len=4), parameter :: source_name(n_sources) = &
           [character(len=4) :: "SAVA","BORF","TEMF","DEFO","PEAT","AGRI"]
      real, parameter :: ilatn=.25 , ilonn=.25
      real, parameter :: e_lon = -179.875, s_lat = -89.875

      ! Emission factors for individual sources comprise a 41x6 matrix
      real :: ef_sources(n_species, n_sources), ef
      real :: DM_part(nlon,nlat,n_sources)
      real :: DM_emiss(nlon,nlat), frac_contrib(nlon,nlat)
      real :: drybbc(nspecies,nlon,nlat)
      real :: tmp_array(nlon,nlat)
      real, save :: typical_fire_size(0:3)
      real :: lons(nlon), lats(nlat)
      real :: dlon1,dlon2,dlat1,dlat2,TX(nspecies)
      logical there
      integer :: i, j, n, i1, i2, j1, j2
      integer :: iveg_ag, ispc, ic, jc, icount
      integer :: iveg, veg_type(nlon,nlat)
      integer :: numSpecies, id1, id2
      integer :: num_days_month
      ! names of the different gas and aerosol species
      character(len=20), allocatable :: Sources(:)
      character(len=20), allocatable :: Species(:)
      character(len=2) :: cmon
      character(len=:), allocatable :: cday
      character(len=256) :: ef_filename, hdf5File, text

      call get_num_species(gfedv4_suffix, numSpecies)
      allocate(Species(numSpecies))
      call get_species(gfedv4_suffix, numSpecies, Species)

      if(use_bbem_plumerise == 1) call mem_bbbem_plume(n1,n2,n3)

      do i=1,nlon
         lons(i)= e_lon + (i-1)*ilonn
      enddo
      do j=1,nlat
         lats(j)= s_lat + (j-1)*ilatn
      enddo

      if (ng == 1) then  

         print*,'==================   GFEDv4  ==================================='

         ef_filename = trim(gfedv4_data_dir)//'/GFED4_Emission_Factors.txt'
         inquire(file=ef_filename, exist=there)
         If (.not.there) then	
            print*, trim(Iam)//trim(ef_filename)//' does not exist.'
            stop "stopping...."
         end if
         hdf5File = get_hdf5_filename(gfedv4_data_dir, iyear)
         inquire(file=trim(hdf5File), exist=there)
         if (.not.there) then	
            print*, trim(Iam)//trim(hdf5File)//' does not exist.'
            stop "stopping...."
         end if
         
         write(6,*)trim(Iam)//' Read emission factors file'
         call ef_read(ef_filename, ef_sources, n_species, n_sources)

         write(6,*)trim(Iam)//' Get HDF5 file handle'
         hdf5File = get_hdf5_filename(gfedv4_data_dir, iyear)
         
         frac_contrib = 0.0
         DM_part = 0.0
         DM_emiss = 0.0
         drybbc = 0.0

         write(cmon,'(i2.2)') imon
         text = "/emissions/"//cmon//"/DM"
         write(6,'(1x,a)') trim(Iam)//' Read '//trim(text)
         call H5ReadDataset(hdf5File, trim(text), DM_emiss)
         
         ! Biomass DM emissions per month (in kg DM/m2/month)
         !DM_emiss = hdf5_ptr

         if (source_option == 7) then
            ! If source_option option equals 7 (TOTL) then we have to read emissions
            ! for all biome sources
            do i = 1, n_sources
               text = "/emissions/"//cmon//"/partitioning/DM_"//source_name(i)
               write(6,'(1x,a)') trim(Iam)//' Read '//trim(text)
               ! DM emissions by fire source (fraction)
               call H5ReadDataset(hdf5File, trim(text), DM_part(:,:,i))
            end do
         else 
            ! just get the emissions from the specified source_option
            text = "/emissions/"//cmon//"/partitioning/DM_"//source_name(source_option)
            write(6,'(1x,a)') trim(Iam)//' Read '//trim(text)
            ! DM emissions by fire source (fraction)
            call H5ReadDataset(hdf5File, trim(text), DM_part(:,:,1))
         end if
         
         if (iyear > 2002) then
            if (iday < 10) then
               Allocate(character(len=1) :: cday)
               write(cday,'(i1.1)') iday
            else
               allocate(character(len=2) :: cday)
               write(cday,'(i2.2)') iday
            end if
            text = "/emissions/"//cmon//"/daily_fraction/day_"//cday
            deallocate(cday)
            write(6,'(1x,a)') trim(Iam)//' Read '//trim(text)
            call H5ReadDataset(hdf5File, trim(text), frac_contrib)
         end if
         
         num_days_month = getDaysInMonth_(iyear, imon)
         ! Calculate emissions as the product of monthly biomass burning dry matter
         ! (DM) emissions (kg DM per m2 per month), the fraction the specific source
         ! contributes to this (frac_contrib, unitless) and the emission factor
         ! (g species per kg DM burned).
         do n = 1, numSpecies
            tmp_array = 0.0
            write(6,*)trim(Iam)//' Calculate emissions for ', trim(Species(n))
            
            ! Several GFEDV4 species are divided equally into two AeM species.
            call get_gfedv4_identity(Species(n), id1, id2)
            
            if (id2 > 0 .and. id1 > 0) then
               
               ! Note that for years <= 2002 we do not have frac_contrib, so we scale
               ! total emissions by the number of days in the month.
               if (source_option == 7) then
                  do i = 1, n_sources
                     ef = ef_sources(n, i)
                     if (iyear > 2002) then
                        tmp_array = tmp_array + &
                             DM_part(:,:,i) * DM_emiss * ef * frac_contrib
                     else
                        tmp_array = tmp_array + &
                             DM_part(:,:,i) * DM_emiss * ef / float(num_days_month)
                     end if
                  end do
               else
                  ef = ef_sources(n, source_option)
                  if (iyear > 2002) then
                     tmp_array = DM_part(:,:,1) * DM_emiss * ef * frac_contrib
                  else
                     tmp_array = DM_part(:,:,1) * DM_emiss * ef / float(num_days_month)
                  end if
               end if
               drybbc(id1,:,:) = 0.5 * tmp_array ! g[DM]/m^2/day
               drybbc(id2,:,:) = 0.5 * tmp_array ! g[DM]/m^2/day

            else if (id1 > 0) then

               if (source_option == 7) then
                  do i = 1, n_sources
                     ef = ef_sources(n, i)
                     ! For years <= 2002 we do not have frac_contrib, so we scale
                     ! total emissions by the number of days in the month.
                     if (iyear > 2002) then
                        tmp_array = tmp_array + &
                             DM_part(:,:,i) * DM_emiss * ef * frac_contrib
                     else
                        tmp_array = tmp_array + &
                             DM_part(:,:,i) * DM_emiss * ef / float(num_days_month)
                     end if
                  end do
               else
                  ef = ef_sources(n, source_option)
                  if (iyear > 2002) then
                     tmp_array = DM_part(:,:,1) * DM_emiss * ef * frac_contrib
                  else
                     tmp_array = DM_part(:,:,1) * DM_emiss * ef / float(num_days_month)
                  end if
               end if
               drybbc(id1,:,:) = tmp_array ! g[DM]/m^2/day

            end if
         enddo
         
         !- for plumerise - typical fire size
         typical_fire_size(0)=0.
         typical_fire_size(1)=20.* 1.e4 ! 20ha for tropical forest
         typical_fire_size(2)=20.* 1.e4 ! 20ha for extra tropical
         typical_fire_size(3)=5. * 1.e4 ! 5ha for savana/past ...

         !-revert y-dir
         call yrevert(nlon,nlat,1,veg_type)
         do i = 1, numSpecies
            call get_gfedv4_identity(Species(i), id1, id2)
            ! EMK NUWRF...Some GFED species are split into two AeM species
            if (id2 > 0 .and. id1 > 0) then
               call yrevert(nlon,nlat,1,drybbc(id1,:,:))
               call yrevert(nlon,nlat,1,drybbc(id2,:,:))
            else if (id1 > 0) then
               call yrevert(nlon,nlat,1,drybbc(id1,:,:))
            end if
         enddo
      endif !(ng==1) 


      !--- performs the interpolation to model grid box

      do i=1,n1
         do j=1,n2

            call get_index1(i,j,nlon,nlat,n1,n2,rlat,rlon,lons,lats &
                 ,ilatn,ilonn,dlat1,dlat2,dlon1,dlon2,i1,j1,i2,j2,ic,jc)
            !print*,'ic jc=',ic,jc,rlat(i,j),rlon(i,j)
            if(ic == 0 .or. jc == 0. .or. ic > nlon .or. jc > nlat) cycle
            ! RAWsrc contains CO data 
            call interpol_gfedv4(i,j,n1,n2,rlon,rlat,ic,jc,nlat,nlon, &
                 ilatn,ilonn,nspecies,drybbc,veg_type,tx)
            ! 
            ! TX is the value interpolated to the model grid box.
            ! Convert from g/m^2/day to kg/m^2/day
            do ispc=1,nspecies
               gfedv4_g(ispc)%src(i,j,1) = TX(ispc)*1.e-3
            enddo

            !------- plumerise section - only for areas with burning
            if(use_bbem_plumerise /= 1 .or. veg_type(ic,jc) == 0 &
                 .or. gfedv4_g(1)%src(i,j,1) < 1.e-6) cycle

            iveg_ag = veg_type(ic,jc) 

            do ispc=1,nspecies
               bbbem_plume_g(ispc,iveg_ag)%src(i,j,1) = &
                    gfedv4_g(ispc)%src(i,j,1)*flaming(iveg_ag)
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
                     if(gfedv4_g(ispc)%src(i,j,1) > 0.) &	  
                                !- fraction consumed at the flaming phase:
                          bbbem_plume_g(ispc,iveg_ag)%src(i,j,1) = &
                          bbbem_plume_g(ispc,iveg_ag)%src(i,j,1) &
                          / gfedv4_g(ispc)%src(i,j,1)
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
                     if(gfedv4_g(ispc)%src(i,j,1) > 0.) then
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

      if(grid_type == 'rams' .or. grid_type == 'polar') then 
         ! only for 'rams' until rland is also defined for others grids 
         do ispc=1,nspecies 
            call apply_land_restriction(n1,n2,rland,gfedv4_g(ispc)%src(:,:,1))
         enddo
      endif
      
   contains

      function getDaysInMonth_(year, month)
         integer, intent(in) :: year, month
         integer :: getDaysInMonth_

         integer :: extraDay
         integer :: daysInMonth(12)
         ! Leap year determination
         extraDay = 0
         if (mod(year,4).eq.0) then
            extraDay = 1
            if (mod(year,100).eq.0) then
               extraDay = 0
               if (mod(year,400).eq.0) then
                  extraDay = 1
               endif
            endif
         endif

         ! Set number of days in each month
         daysInMonth(1)  = 31
         daysInMonth(2)  = 28 + extraDay
         daysInMonth(3)  = 31
         daysInMonth(4)  = 30
         daysInMonth(5)  = 31
         daysInMonth(6)  = 30
         daysInMonth(7)  = 31
         daysInMonth(8)  = 31
         daysInMonth(9)  = 30
         daysInMonth(10) = 31
         daysInMonth(11) = 30
         daysInMonth(12) = 31

         getDaysInMonth_ = daysInMonth(month)

      end function getDaysInMonth_

   end subroutine read_gfedv4

   subroutine get_gfedv4_identity(spc_name,ident,ident2)
      use AeM_emission_factors, only : gfedv4_spc_name=>AeM_spc_name
      use AeM_emission_factors,  gfedv4_nspecies=>AeM_nspecies
      integer isp
      character (len=*), intent(in)  :: spc_name
      integer                        :: ident
      integer                        :: ident2
      integer :: ii

      do isp = 1,gfedv4_nspecies
         ident=-1
         ident2=-1
         if(spc_name == gfedv4_spc_name(isp)) then
            print*,'Matched GFEDV4 species with AeM species ',trim(spc_name)
            ident=isp
            return
         endif
      end do
      ! EMK NUWRF Bug fix...Not all GFEDV4 names exactly match the AeM list.
      ! So some special double checks are needed.
      ! We'll have extensive checking to better detect breakage if the AeM
      ! list changes in the future.
      ! ZTao modified AeM -> gfedv4 species mapping
      ! TODO: Handle C, DM, CH3COCHO

      if (trim(spc_name) == 'C2H4O') then
         ident = special_identity('Acetald',gfedv4_spc_name, &
              gfedv4_nspecies)         
      else if (trim(spc_name) == 'C2H5OH') then
         ident = special_identity('Ethanol',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'C2H6S') then
         ident = special_identity('DMS',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'C3H6O') then
         ident  = special_identity('Propanal',gfedv4_spc_name, &
              gfedv4_nspecies)
         ident2 = special_identity('Acetone',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'C5H8') then
         ident = special_identity('Isoprene',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'CH2O') then
         ident = special_identity('Formaldehyde',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'CH3OH') then
         ident = special_identity('Methanol',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'NMHC') then
         ident = special_identity('NHMC',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'PM2.5') then
         ident = special_identity('BBURN2',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'C10H16') then
         ident = special_identity('terpenes',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'C8H10') then
         ident  = special_identity('xylenes',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'C6H6') then
         ident = special_identity('benzene',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'C7H8') then
         ident = special_identity('toluene',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'Toluene_lump') then
         ident = special_identity('toluene',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'TPM') then
         ident = special_identity('BBURN3',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'TPC') then
         ident = special_identity('TC',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'HCOOH') then
         ident = special_identity('HFo',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'CH3COOH') then
         ident = special_identity('HAc',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'MEK') then
         ident = special_identity('2_Butanone',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'HOCH2CHO') then
         ident = special_identity('Hydroxyacetaldehyde',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'Higher_Alkenes') then
         ident = special_identity('1_hexene',gfedv4_spc_name, &
              gfedv4_nspecies)
         ident2 = special_identity('octenes',gfedv4_spc_name, &
              gfedv4_nspecies)
      else if (trim(spc_name) == 'Higher_Alkanes') then
         ident = special_identity('n_hexane',gfedv4_spc_name, &
              gfedv4_nspecies)
         ident2 = special_identity('heptane',gfedv4_spc_name, &
              gfedv4_nspecies)
      end if
      if (ident < 0) then
         print*,'Warning, cannot process GFEDV4 species ',trim(spc_name)
      end if

   contains

      ! EMK NUWRF
      ! Internal function for special comparisons with AeM species list.
      ! Assumes a match will be found with AeM. A fatal error will occur
      ! if this is not the case -- useful for detecting future changes to the
      ! AeM species list.
      function special_identity(spc_name,gfedv4_spc_name, gfedv4_nspecies) &
           result(ident)
         character(len=*),intent(in) :: spc_name
         integer,intent(in) :: gfedv4_nspecies
         character(len=*),intent(in) :: gfedv4_spc_name(gfedv4_nspecies)
         integer :: ident
         integer :: ii
         ident=-1
         do ii = 1,gfedv4_nspecies
            if (trim(spc_name) == trim(gfedv4_spc_name(ii))) then
               print*,'Matched GFEDV4 species with AeM species ',trim(spc_name)
               ident=ii
               return               
            end if
         end do
         print*,'INTERNAL ERROR, cannot match GFEDV4 species with AeM species', &
              trim(spc_name)
         stop 
      end function special_identity

   end subroutine get_gfedv4_identity

   subroutine interpol_gfedv4(i,j,n1,n2,rlon,rlat,ic,jc,nlat,nlon,ilatn,ilonn &
        ,nspecies,drybbc,veg_type &
        ,tx)
      use AeM_emission_factors, N2_nitrogenio=>N2
      use grid_dims_out, only: grid_type
      integer n1,n2,ic,jc,nlat,nlon,i,j,nspecies,ispc
      real, dimension(n1,n2) :: rlat,rlon
      real ilatn,ilonn,tx(nspecies)
      real, dimension(nspecies,nlon,nlat) :: drybbc
      integer  veg_type(nlon,nlat)

      !-local var
      real dlonr,dlatr,usdum
      integer qi1,qi2,qj1,qj2,ncount,ii,jj
      integer iveg , veg_type_model(n1,n2)

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
            ! DM. But for now, just use the emissions calculated by the GFEDV4
            ! team.
            TX(:) = TX(:) + drybbc(:,ii,jj)
            !      endif 
         enddo
      enddo
      TX(:) = TX(:) / (float(ncount) + 1.E-10) ! interpolated rate
      
   end subroutine interpol_gfedv4

   subroutine get_num_species(gfedv4_suffix, numSpecies)
      character(len=*), intent(in) :: gfedv4_suffix
      integer :: numSpecies

      character(len=len_trim(gfedv4_suffix)) :: local_gfedv4_suffix
      integer :: localPosition, localLength

      numSpecies = 0
      local_gfedv4_suffix = trim(gfedv4_suffix)
      do
         localPosition = scan(trim(local_gfedv4_suffix),",")
         localLength = len_trim(local_gfedv4_suffix)

         if (localPosition > 0) then
            numSpecies = numSpecies + 1
            local_gfedv4_suffix = &
                 local_gfedv4_suffix(localPosition+1:localLength) 
            localPosition = 0
         else
            exit
         endif
      enddo

      numSpecies = numSpecies+1

   end subroutine get_num_species

   subroutine get_species(gfedv4_suffix,  numSpecies,Species)
      character(len=*), intent(in) :: gfedv4_suffix
      integer :: numSpecies
      character(len=*) :: Species(numSpecies)

      character(len=len_trim(gfedv4_suffix)) :: local_gfedv4_suffix
      integer :: localPosition, localLength, iSpecies

      local_gfedv4_suffix = trim(gfedv4_suffix)
      do iSpecies = 1, numSpecies - 1

         localPosition = scan(trim(local_gfedv4_suffix),",")
         localLength = len_trim(local_gfedv4_suffix)
         if (localPosition > 0) then
            Species(iSpecies) = &
                 local_gfedv4_suffix(1:localPosition-1)
            local_gfedv4_suffix = &
                 local_gfedv4_suffix(localPosition+1:localLength) 
            localPosition = 0
            print *,'DEBUG: ',trim(Species(iSpecies))
         endif

      enddo

      localPosition = scan(trim(gfedv4_suffix),",",back=.true.)
      Species(numSpecies) = &
           gfedv4_suffix(localPosition+1:len_trim(gfedv4_suffix))

   end subroutine get_species

   function get_hdf5_filename(data_dir, iyear) result(filename)
      ! Given the data path and year, construct hdf5 filename 
      character (len=*), intent(in)  :: data_dir
      integer, intent(in) :: iyear
      character (len=256)  :: filename
      character (len=4)  :: cyear      
      character (len=:), allocatable :: str1, str2

      write(cyear,'(i4.4)') iyear
      str1 = trim(data_dir)
      str2 = trim('/GFED4.1s_'//cyear//'.hdf5')      
      filename = trim(str1)//trim(str2)

   end function get_hdf5_filename

   subroutine ef_read(ef_filename, efs, n_sp, n_so)
      ! Read emission factors for individual species by biome (fire source)
      ! from text file.
      character (len=*), intent(in)  :: ef_filename
      integer, intent(in) :: n_sp, n_so
      real, intent(inout) :: efs(n_sp, n_so) ! 41 species, 6 sources
      character(len=128) :: buffer, label
      character(LEN=80)  :: token
      integer :: pos, i, j
      character(len=1)  :: delimiter
      ! strtok is in utils/shared
      external strtok
      character(len=255) :: strtok

      efs = 0.0
      !print *,trim(ef_filename),shape(efs)
      delimiter = " "
      open(unit=10, file=trim(ef_filename), form='formatted')
      ! Loop over entries in file.
      i = 1
      do 
         read (10, '(a)', end=99) buffer
         ! First of all skip all the lines starting with a '#'
         if (buffer(1:1) == '#') then
            cycle
         else
            pos = index(buffer, ' ')
            j = 0
            ! First column is a string
            token = strtok(buffer, delimiter)
            do while (token .ne. char(0))
               ! Extract next 6 tokens: those are the emission factors
               if (j>0) read(token(1:len_trim(token)), '(f12.0)') efs(i, j)
               token = strtok(char(0), delimiter)
               j = j + 1
            end do
            i = i + 1
         end if
      end do
99    close(10)
      print *,'Read emission factors:'
      do i=1,n_sp
         write(*,'(20f10.2)') efs(i,:)
      end do

   end subroutine ef_read

end module gfedv4_emissions
