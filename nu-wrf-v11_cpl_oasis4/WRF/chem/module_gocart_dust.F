MODULE GOCART_DUST
  

  USE module_data_gocart_dust

CONTAINS
  subroutine gocart_dust_driver(ktau,dt,config_flags,julday,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,dustin,           &
         ivgtyp,isltyp,vegfra,xland,xlat,xlong,gsw,dx,g,emis_dust,        &
#ifdef WRF_LIS
         lisPorosity, & ! NUWRF
#endif
         ust, ust_t, & ! NUWRF, DEAD scheme
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_configure
  USE module_state_description
  USE module_model_constants, ONLY: mwdry
  IMPLICIT NONE
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: julday, ktau,                     &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL,&
         INTENT(INOUT ) ::                                                 &
         emis_dust
  REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) ,      &
      INTENT(INOUT) ::                               smois
   REAL,  DIMENSION( ims:ime , jms:jme, 3 )                   ,               &
          INTENT(IN   ) ::    erod
   REAL,  DIMENSION( ims:ime , jms:jme, 5 )                   ,               &
          INTENT(INout   ) ::    dustin
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &  !NU-WRF DEAD scheme
          INTENT(INout   ) ::   ust_t                                         !registered by uoc 
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     ust,                  &
                                                     gsw,                  &
                                                  vegfra,                  &
                                                     xland,                &
                                                     xlat,                 &
                                                     xlong
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                      t_phy,               &
                                                     dz8w,p8w,             &
                                              u_phy,v_phy,rho_phy
 
  REAL, INTENT(IN   ) :: dt,dx,g
#ifdef WRF_LIS
  ! NU-WRF from LIS
  real, intent(in) :: lisPorosity(ims:ime,jms:jme)
#endif
!
! local variables
!
  integer :: nmx,i,j,k,ndt,imx,jmx,lmx
  integer,dimension (1,1) :: ilwi
  real*8, DIMENSION (1,1,3,1) :: erodin
  real*8, DIMENSION (5) :: tc,bems
  real*8, DIMENSION (5) :: tc_old,dtc  !dkim6, NU-WRF
  integer*8             :: dustEmissLayer,dl   !dkim6, NU-WRF
  real*8, dimension (1,1) :: my_smois,my_porosity !dkim6, NU-WRF, DEAD scheme
  real,   dimension (1,1) :: my_ust               !dkim6, NU-WRF, DEAD scheme
  real                    :: my_ust_t             !dkim6, NU-WRF, DEAD scheme
  real*8, dimension (1,1) :: w10m,gwet,airden,airmas
  real*8, dimension (1) :: dxy
  real*8  conver,converi
  real dttt
  integer ibegc,jbegc,iendc,jendc
  jbegc=max(jts,jds+3)
  jendc=min(jte,jde-4)
  ibegc=max(its,ids+3)
  iendc=min(ite,ide-4)
! conver=1.e-9*mwdry
! converi=1.e9/mwdry
  conver=1.e-9
  converi=1.e9
!
! number of dust bins
!
  imx=1
  jmx=1
  lmx=1
  nmx=5 
  k=kts
  do j=jbegc,jendc
  do i=ibegc,iendc
!
! no dust over water!!!
!
     if(xland(i,j).lt.1.5)then
     ilwi(1,1)=1
     if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
        tc(:)=1.e-16*conver
     else
        tc(1)=chem(i,kts,j,p_dust_1)*conver
        tc(2)=chem(i,kts,j,p_dust_2)*conver
        tc(3)=chem(i,kts,j,p_dust_3)*conver
        tc(4)=chem(i,kts,j,p_dust_4)*conver
        tc(5)=chem(i,kts,j,p_dust_5)*conver
        ! dkim6, 20141115, distirbute to dust Emission layers
        tc_old(1)=tc(1)
        tc_old(2)=tc(2)
        tc_old(3)=tc(3)
        tc_old(4)=tc(4)
        tc_old(5)=tc(5)
     endif
     w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
     airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g
!
! don't trust the u10,v10 values, is model layers are very thin near surface
!
     if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
     erodin(1,1,1,1)=erod(i,j,1)!/dx/dx
     erodin(1,1,2,1)=erod(i,j,2)!/dx/dx
     erodin(1,1,3,1)=erod(i,j,3)!/dx/dx
!
!  volumetric soil moisture over porosity
!
! NU-WRF:  If running LIS, use the porosity values from that software
     if (config_flags%sf_surface_physics .eq. NUWRFLISSCHEME) then
#ifdef WRF_LIS
        gwet(1,1)=smois(i,1,j)/lisPorosity(i,j)
#endif
     else
        gwet(1,1)=smois(i,1,j)/porosity(isltyp(i,j))
     end if

     !--dkim6, for DEAD scheme
     if (config_flags%sf_surface_physics .eq. NUWRFLISSCHEME) then
#ifdef WRF_LIS
        my_porosity(1,1)=lisPorosity(i,j)
#endif
     else
        my_porosity(1,1)=porosity(isltyp(i,j))
     end if
     my_smois(1,1)=smois(i,1,j)
     my_ust(1,1)=ust(i,j)

     gwet(1,1)=smois(i,1,j)/porosity(isltyp(i,j))
     ndt=ifix(dt)
     airden(1,1)=rho_phy(i,kts,j)
     dxy(1)=dx*dx
!    if(erod(i,j,1).gt.0.)write(0,*)'er1=',p_dust_1,num_chem,erod(i,j,1),tc(2)
!    if(erod(i,j,1).gt.0.)write(0,*)'er1=',dt,dxy(1),u10(i,j),w10m(1,1)
!    erodin(1,1,1,1)=  0.149748762553862
!    erodin(1,1,2,1)=  7.487438878070708E-002
!    erodin(1,1,3,1)=  7.487438878070708E-002
!    ilwi(1,1)=           1
!    dxy(1)=   54585850453.7552
!    w10m(1,1)=   10.6305338763678
!    gwet(1,1)=  9.136307984590530E-002
!    airden(1,1)=   1.16423276395132
!    airmas(1,1)=   8114017750938.79
!    tc (1) =  1.000000000000000D-030
!    tc (2) =  1.000000000000000d-030
!    tc (3) =  1.000000000000000d-030
!    tc(4)  =  1.000000000000000d-030
!    tc(5)  =  1.000000000000000d-030
!    dttt=3600.
  !dkim6, gocart_dustemiss_scheme (1:GOCART, 2:DEAD) 
  if (config_flags%gocart_dustemiss_scheme == 1 ) then
    call source_du( imx,jmx,lmx,nmx, dt, tc, &
                     erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,config_flags%start_month,g, &
                     config_flags%gocart_dustemiss_gwetthresh, &
                     config_flags%gocart_dustemiss_resfactor)
  else if (config_flags%gocart_dustemiss_scheme == 2 ) then 
    call source_du_dead( imx,jmx,lmx,nmx, dt, tc, &
                     erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,config_flags%start_month,g, &
                     my_ust, my_ust_t, my_smois, my_porosity, &
                     config_flags%gocart_dustemiss_gwetthresh, &
                     config_flags%gocart_dustemiss_resfactor)
  else
    write(wrf_err_message,*) 'ERROR, gocart_dustemiss_scheme = ', &
          config_flags%gocart_dustemiss_scheme, ', should be 1:GOCART or 2:DEAD'
    call wrf_error_fatal(wrf_err_message)
  endif

!    write(0,*)tc(1)
!    write(0,*)tc(2)
!    write(0,*)tc(3)
!    write(0,*)tc(4)
!    write(0,*)tc(5)
!    if(erod(i,j,1).gt.0.)write(0,*)'er2=',i,j,erod(i,j,1),tc(2)
     if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
      dustin(i,j,1:5)=tc(1:5)*converi
     else
!dkim6, comment out
!     chem(i,kts,j,p_dust_1)=tc(1)*converi
!     chem(i,kts,j,p_dust_2)=tc(2)*converi
!     chem(i,kts,j,p_dust_3)=tc(3)*converi
!     chem(i,kts,j,p_dust_4)=tc(4)*converi
!     chem(i,kts,j,p_dust_5)=tc(5)*converi
      !dkim6, distribute dust to upper layers 
      dtc(:)=tc(:)-tc_old(:)
      dustEmissLayer=config_flags%gocart_dustemiss_layer !currently 1, 5, or 8
      if (dustEmissLayer >= kte) then ! NUWRF EMK...Earlier was kte-1
          write(wrf_err_message,*) 'ERROR, gocart_dustemiss_layer = ', &
          dustEmissLayer,', should be less than ',kte,'!'
          call wrf_error_fatal(wrf_err_message)
      end if
      do dl = kts, dustEmissLayer ! NUWRF EMK...Earlier was dustEmissLayer-1 
       chem(i,dl,j,p_dust_1)= chem(i,dl,j,p_dust_1)+dtc(1)*converi/dustEmissLayer !update chem here, using chem(i,k,j,m), not tc(m)
       chem(i,dl,j,p_dust_2)= chem(i,dl,j,p_dust_2)+dtc(2)*converi/dustEmissLayer
       chem(i,dl,j,p_dust_3)= chem(i,dl,j,p_dust_3)+dtc(3)*converi/dustEmissLayer
       chem(i,dl,j,p_dust_4)= chem(i,dl,j,p_dust_4)+dtc(4)*converi/dustEmissLayer
       chem(i,dl,j,p_dust_5)= chem(i,dl,j,p_dust_5)+dtc(5)*converi/dustEmissLayer
      enddo
     endif
! for output diagnostics

! NU-WRF:  Output this, but make sure optional variable is available.
     if (present(emis_dust)) then
        emis_dust(i,1,j,p_edust1)=bems(1)
        emis_dust(i,1,j,p_edust2)=bems(2)
        emis_dust(i,1,j,p_edust3)=bems(3)
        emis_dust(i,1,j,p_edust4)=bems(4)
        emis_dust(i,1,j,p_edust5)=bems(5)
        ust_t(i,j)=my_ust_t
     end if
     endif
  enddo
  enddo
!

end subroutine gocart_dust_driver

  
  SUBROUTINE source_du( imx,jmx,lmx,nmx, dt1, tc, &
                     erod, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,month,g0, &
                     gocart_dustemiss_gwetthresh,gocart_dustemiss_resfactor)
                     ! NU-WRF:  Added dust emission tuning factors

! ****************************************************************************
! *  Evaluate the source of each dust particles size classes  (kg/m3)        
! *  by soil emission.
! *  Input:
! *         EROD      Fraction of erodible grid cell                (-)
! *                   for 1: Sand, 2: Silt, 3: Clay
! *         DUSTDEN   Dust density                                  (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         AIRVOL    Volume occupy by each grid boxes              (m3)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *         u_tresh   Threshold velocity for particule uplifting    (m/s)
! *         CH_dust   Constant to fudge the total emission of dust  (s2/m2)
! *      
! *  Output:
! *         DSRC      Source of each dust type           (kg/timestep/cell) 
! *
! *  Working:
! *         SRC       Potential source                   (kg/m/timestep/cell)
! *
! ****************************************************************************

! USE module_data_gocart
! USE module_data_gocart_dust

  use module_wrf_error, only:  wrf_err_message
  

  INTEGER, INTENT(IN)    :: nmx,imx,jmx,lmx
  REAL*8,    INTENT(IN)    :: erod(imx,jmx,ndcls,ndsrc)
  INTEGER, INTENT(IN)    :: ilwi(imx,jmx),month

  REAL*8,    INTENT(IN)    :: w10m(imx,jmx), gwet(imx,jmx)
  REAL*8,    INTENT(IN)    :: dxy(jmx)
  REAL*8,    INTENT(IN)    :: airden(imx,jmx,lmx), airmas(imx,jmx,lmx)
  REAL*8,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(OUT)   :: bems(imx,jmx,nmx) 

  ! NU-WRF:  Added tunable dust emission factors
  real, intent(in) :: gocart_dustemiss_gwetthresh
  real, intent(in) :: gocart_dustemiss_resfactor

  REAL*8    :: den(nmx), diam(nmx)
  REAL*8    :: tsrc, u_ts0, cw, u_ts, dsrc, srce
  REAL, intent(in)    :: g0
  REAL    :: rhoa, g,dt1
  INTEGER :: i, j, n, m, k


  REAL*8                  :: tcmw(nmx), ar(nmx), tcvv(nmx)
  REAL*8                  :: ar_wetdep(nmx), kc(nmx)
  CHARACTER(LEN=20)     :: tcname(nmx), tcunits(nmx)
  LOGICAL               :: aerosol(nmx)


! REAL*8 :: tc1(imx,jmx,lmx,nmx)
! REAL*8, TARGET :: tcms(imx,jmx,lmx,nmx) ! tracer mass (kg; kgS for sulfur case)
! REAL*8, TARGET :: tcgm(imx,jmx,lmx,nmx) ! g/m3

  !-----------------------------------------------------------------------  
  ! sea salt specific
  !-----------------------------------------------------------------------  
! REAL*8, DIMENSION(nmx) :: ssaltden, ssaltreff, ra, rb
! REAL*8 :: ch_ss(nmx,12)

  !-----------------------------------------------------------------------  
  ! emissions (input)
  !-----------------------------------------------------------------------  
! REAL*8 :: e_an(imx,jmx,2,nmx), e_bb(imx,jmx,nmx), &
!         e_ac(imx,jmx,lmx,nmx)

  !-----------------------------------------------------------------------  
  ! diagnostics (budget)
  !-----------------------------------------------------------------------
!  ! tendencies per time step and process
!  REAL, TARGET :: bems(imx,jmx,nmx), bdry(imx,jmx,nmx), bstl(imx,jmx,nmx)
!  REAL, TARGET :: bwet(imx,jmx,nmx), bcnv(imx,jmx,nmx)
!
!  ! integrated tendencies per process
!  REAL, TARGET :: tems(imx,jmx,nmx), tstl(imx,jmx,nmx)
!  REAL, TARGET :: tdry(imx,jmx,nmx), twet(imx,jmx,nmx), tcnv(imx,jmx,nmx)

  ! global mass balance per time step 
  REAL*8 :: tmas0(nmx), tmas1(nmx)
  REAL*8 :: dtems(nmx), dttrp(nmx), dtdif(nmx), dtcnv(nmx)
  REAL*8 :: dtwet(nmx), dtdry(nmx), dtstl(nmx)
  REAL*8 :: dtems2(nmx), dttrp2(nmx), dtdif2(nmx), dtcnv2(nmx)
  REAL*8 :: dtwet2(nmx), dtdry2(nmx), dtstl2(nmx)




  ! executable statemenst

  ! NU-WRF:  Sanity check tuning factors
  if (gocart_dustemiss_gwetthresh < 0 .or. &
      gocart_dustemiss_gwetthresh > 1) then
     write(wrf_err_message,*) 'ERROR, gocart_dustemiss_gwetthresh = ', &
          gocart_dustemiss_gwetthresh,', should be from 0 to 1!'
     call wrf_error_fatal(wrf_err_message)
  end if
  if (gocart_dustemiss_resfactor < 0) then
     write(wrf_err_message,*) 'ERROR, gocart_dustemiss_resfactor = ', &
          gocart_dustemiss_resfactor, ', should be > 0'
     call wrf_error_fatal(wrf_err_message)
  end if

  DO n = 1, nmx
     ! Threshold velocity as a function of the dust density and the diameter
     ! from Bagnold (1941)
     den(n) = den_dust(n)*1.0D-3
     diam(n) = 2.0*reff_dust(n)*1.0D2
     g = g0*1.0E2
     ! Pointer to the 3 classes considered in the source data files
     m = ipoint(n)
     tsrc = 0.0
     DO k = 1, ndsrc
        ! No flux if wet soil 
        DO i = 1,imx
           DO j = 1,jmx
              rhoa = airden(i,j,1)*1.0D-3
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 
!             write(0,*)u_ts0,den(n),diam(n),rhoa,g
              ! Fraction of emerged surfaces (subtract lakes, coastal ocean,..)
!              cw = 1.0 - water(i,j)
              
              ! Case of surface dry enough to erode
              ! NU-WRF:  User-defined threshold from namelist
              IF (gwet(i,j) < gocart_dustemiss_gwetthresh) THEN
!              IF (gwet(i,j) < 0.5) THEN  
                 u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet(i,j)))))
              ELSE
                 ! Case of wet surface, no erosion
                 u_ts = 100.0
              END IF
              srce = frac_s(n)*erod(i,j,m,k)*dxy(j)  ! (m2)
              IF (ilwi(i,j) == 1 ) THEN
                 dsrc = ch_dust(n,month)*srce*w10m(i,j)**2 &
                      * (w10m(i,j) - u_ts)*dt1  ! (kg)
!                write(0,*)ch_dust(n,month),srce,w10m(i,j),u_ts,gwet(i,j)
              ELSE 
                 dsrc = 0.0
              END IF
!              dsrc = cw*ch_dust(k)*srce*w10m(i,j)**2 &
!                   * (w10m(i,j) - u_ts)*dt1  ! (kg)
!              dsrc = cw*ch_dust(n,dt(1)%mn)*srce*w10m(i,j)**2 &
!                   * (w10m(i,j) - u_ts)*dt1  ! (kg)
              IF (dsrc < 0.0) dsrc = 0.0
              
              ! Update dust mixing ratio at first model level.
              ! NU-WRF:  User defined factor
!              tc(i,j,1,n) = tc(i,j,1,n) + .2*dsrc / airmas(i,j,1)
!              bems(i,j,n) = .2*dsrc
              tc(i,j,1,n) = tc(i,j,1,n) + &
                   gocart_dustemiss_resfactor*dsrc / airmas(i,j,1)
              bems(i,j,n) = gocart_dustemiss_resfactor*dsrc

           END DO
        END DO
     END DO
  END DO
  
END SUBROUTINE source_du

  SUBROUTINE source_du_dead( imx,jmx,lmx,nmx, dt1, tc, &
                     erod, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,month,g0, &
                     my_ust, my_ust_t, my_smois, my_porosity, &
                     gocart_dustemiss_gwetthresh,gocart_dustemiss_resfactor)
! ****************************************************************************
! *  Evaluate the source of each dust particles size classes  (kg/m3)        
! *  by soil emission.
! *      
! *  Output:
! *         DSRC      Source of each dust type           (kg/timestep/cell) 
! *
! *  Modified based on GEOS, DEAD scheme
! *  20180724, first crack, Dongchul Kim
! ****************************************************************************

  use module_wrf_error, only:  wrf_err_message

  INTEGER, INTENT(IN)    :: nmx,imx,jmx,lmx
  REAL*8,    INTENT(IN)    :: erod(imx,jmx,ndcls,ndsrc)
  INTEGER, INTENT(IN)    :: ilwi(imx,jmx),month

  REAL*8,    INTENT(IN)    :: w10m(imx,jmx), gwet(imx,jmx)
  REAL*8,    INTENT(IN)    :: dxy(jmx)
  REAL*8,    INTENT(IN)    :: airden(imx,jmx,lmx), airmas(imx,jmx,lmx)
  REAL*8,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(OUT)   :: bems(imx,jmx,nmx)
  REAL :: my_ust_t

  ! NU-WRF:  Added tunable dust emission factors
  real, intent(in) :: gocart_dustemiss_gwetthresh
  real, intent(in) :: my_ust(imx,jmx)
  real*8, intent(in) :: my_smois(imx,jmx), my_porosity(imx,jmx)

  REAL*8    :: den(nmx), diam(nmx)
  REAL*8    :: tsrc, u_ts0, cw, u_ts, dsrc, srce
  REAL, intent(in)    :: g0
  REAL    :: rhoa, g,dt1
  INTEGER :: i, j, n, m, k
  ! frac_s_dead, soil fraction from D'Almeida,(1987), used in GEOS
  REAL, DIMENSION (5), PARAMETER :: frac_s_dead(5)=(/0.112,0.232,0.296,0.277,0.064/)

   !  Variables and parameters specific to Zender source implementation
   real            ::  fd               ! drag partitioning eff. factor
                     ! dust effective diameter for monomodal soil
   real, parameter ::  soil_diameter = 75.e-6
                     ! soil grain particle and water density [kg m-3]
   real, parameter ::  soil_density  = 2650.
   real, parameter ::  water_density = 1000.
                     ! scheme defendency [Zender eq 17, Nowottnick 2010 Talbe 1]
                     !                   [tested for Arizonda dust storm case]
   real, parameter ::  scheme_factor = 20.0E4 
   real, parameter ::  zoms = 33.e-6     ! smooth roughness length [m]
   real, parameter ::  zom = 100.e-6     ! roughness length [m]
   real, parameter ::  mclay = 0.2       ! mass fraction of clay
   real            ::  fw               ! water eff. factor
   real            ::  wt               ! Threshold water content
   real            ::  u_thresh_drag
   real            ::  u_thresh_drag_water
   real            ::  soil_bulk_density
   real            ::  gvsmois
   real            ::  ustars           ! Modified fric. vel for Owen Effect
   real            ::  k_z              ! Log profile for nonsaltating
   real            ::  wt10m            ! Threshold 10m wind speed
   real, parameter ::  cs=2.61          ! Constant for horiz flux
   real            ::  rat              ! Ratio of threshold to friction vel
   real            ::  horiz_flux       ! Horizontal Mass Flux
   real            ::  vert_flux        ! Vertical Mass Flux
   real            ::  alpha            ! Vertical Flux Conversion Factor [m-1]

  ! executable statemenst

  !print *, 'I am in DEAD scheme ',nmx,ndsrc,imx,jmx

  ! NU-WRF:  Sanity check tuning factors
  if (gocart_dustemiss_gwetthresh < 0 .or. &
      gocart_dustemiss_gwetthresh > 1) then
     write(wrf_err_message,*) 'ERROR, gocart_dustemiss_gwetthresh = ', &
          gocart_dustemiss_gwetthresh,', should be from 0 to 1!'
     call wrf_error_fatal(wrf_err_message)
  end if
  if (gocart_dustemiss_resfactor < 0) then
     write(wrf_err_message,*) 'ERROR, gocart_dustemiss_resfactor = ', &
          gocart_dustemiss_resfactor, ', should be > 0'
     call wrf_error_fatal(wrf_err_message)
  end if

  !  Calculate drag partitioning efficiency to represent sink of atmos momentum 
  !  to nonerodible roughness elements
  !  Assumes constant smooth roughness length zoms and constant roughness
  !  length for momentum transfer zom
  !  Zender 2003 eq. 3
     fd = ( 1.0 - &
            ( log( zom/zoms ) / &
              log( 0.35 * ( (0.1/zoms)**0.8 ) ) &
            ) &
          ) ** (-1.)

  !  Volumetric to gravimetric soil moisture conversion (Darmenova 2009 eq. 4)
     soil_bulk_density=soil_density*(1.-my_porosity(imx,jmx))
     gvsmois=my_smois(imx,jmx)*water_density/soil_bulk_density

  !  Calculate the threshold water content following Fecan et al. [1999]
  !  Assumes a globally constant mass fraction of clay in soil, mclay
  !  Zender 2003 eq. 5 (implicit a = 1)
     wt = 0.17 * mclay + 0.14 * mclay**2
  !
  !  Compute parameter for scaling horizontal to vertical mass flux,
  !  alpha = [1/m], Zender 2003 eq. 11
     alpha = 100. * exp( (13.4*mclay - 6.)*log(10.) )
  !
  !  Calculate constant k for nonsaltating profile
  !  Gillette et al. [1998] eq. 3
     k_z = 0.4 / log(10./zom)   ! k_z ~ 0.035, where zom=100 (um), wt10m=ts0*fd*fw*28.6, always w10m<wt10m


  DO n = 1, nmx
     ! Threshold velocity as a function of the dust density and the diameter
     ! from Bagnold (1941)
     ! assuming soil particle is uniform single size 
     den(n) = soil_density*1.0D-3          ![g m-3]
     diam(n) = soil_diameter*1.0D2         ![cm]
     g = g0*1.0E2                          ![cm s-2]
     ! Pointer to the 3 classes considered in the source data files
     m = ipoint(n)
     tsrc = 0.0
     DO k = 1, ndsrc
        !  note: imx=1, jmx=1
        DO i = 1,imx
           DO j = 1,jmx
              rhoa = airden(i,j,1)*1.0D-3  ![g/m3]
              !  "1.0D-2", u_ts0 is calculated in cgs, outcome is in [m/s]
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0)

              !  Apply the drag partitioning correction
              !  Equivalent to Marticorena [1997] eq. 4, where our fd = 1/(their)feff
              u_thresh_drag = u_ts0*fd

              !  Modify the threshold depending on soil moisture as in Fecan 1999
              !  Zender 2003 eq. 6, Darmenova 2009 eq. 4
              !if (gwet(i,j) <= wt) then
              if (gvsmois <= wt) then
                fw = 1.0
              else
                !
                !fw = sqrt( 1. + 1.21 * (100. * (gwet(i,j)-wt) )**0.68)   !GEOS
                fw = sqrt( 1. + 1.21 * (100. * (gvsmois-wt) )**0.68)      !NUWRF
              endif
              u_thresh_drag_water = u_thresh_drag*fw

              !  output final threshold friction velocity (ust_t)
              my_ust_t = u_thresh_drag_water

              !  Modify friction velocity for Owen Effect
              !  Assumption of stable atmospheric profile to go from saltation
              !  wind speed to equivalent threshold at z = 10m
              !  Gillette et al. [1998] eq. 3
              wt10m = u_thresh_drag_water/k_z
              if (w10m(i,j) >= wt10m) then
                ustars = my_ust(i,j) + 0.003*((w10m(i,j)-wt10m)**2)
              else
                ustars = my_ust(i,j)
              endif

              !  Calculate the horizontal mass flux of dust [kg m-1 s-1]
              !  Marticorena et al. 1997 eq. 5 
              rat = u_thresh_drag_water / ustars
              if (rat < 1.0) then
                horiz_flux = cs * airden(i,j,1) * ustars**3 /g0 * &
                              (1 - rat**2) * (1+rat)
                !print*, 'gvsmois ',gvsmois,my_smois(i,j),wt,gwet(i,j)
                !print*, 'ust_t ',u_ts0,my_ust_t,fd,fw
                !print*, 'ustars ',ustars,my_ust(i,j),w10m(i,j),wt10m,k_z
                !print*, 'horiz_flux ',horiz_flux,ustars,rat,g0,cs,airden(i,j,1)
              else
                horiz_flux = 0.0
              endif

              !  Calculate the vertical mass flux of dust and scale to source [kg m-2 s-1]
              !  alpha = 0.0479
              vert_flux = alpha * horiz_flux

              !  This is optional...
              IF (gwet(i,j) >= gocart_dustemiss_gwetthresh) THEN
                vert_flux = 0.
              ENDIF

              srce = frac_s_dead(n)*erod(i,j,m,k)*dxy(j)  ! (m2)
              IF (ilwi(i,j) == 1 ) THEN
                 !dsrc = ch_dust(n,month)*srce*w10m(i,j)**2 &           [GOCART]
                 !     * (w10m(i,j) - u_ts)*dt1  ! (kg)
                 dsrc = ch_dust(n,month)*srce*vert_flux*dt1 ! (kg)      [DEAD]
                 dsrc = dsrc*scheme_factor                  ! (kg)      [DEAD]
                 !print*, 'dsrc ',ch_dust(n,month),srce,vert_flux,dt1
              ELSE
                 dsrc = 0.0
              END IF
              IF (dsrc < 0.0) dsrc = 0.0

              ! Update dust mixing ratio at first model level.
              tc(i,j,1,n) = tc(i,j,1,n) + &
                   gocart_dustemiss_resfactor*dsrc / airmas(i,j,1)
              bems(i,j,n) = gocart_dustemiss_resfactor*dsrc
              !print*, 'bems ',dsrc,gocart_dustemiss_resfactor,bems(i,j,n)

           END DO
        END DO
     END DO
  END DO

END SUBROUTINE source_du_dead

END MODULE GOCART_DUST
