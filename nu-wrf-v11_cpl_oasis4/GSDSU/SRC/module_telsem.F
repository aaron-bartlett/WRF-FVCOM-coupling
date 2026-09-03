!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

 module module_telsem

!MODULE mod_mwatlas
!   This module has been written by
!   Catherine Prigent
!   Elodie Joumouille
!   Filipe Aires
!
!   September 2009
!
  IMPLICIT none

  TYPE atlas_emis_mw
    ! number of lines in the atlas
    INTEGER :: ndat
    ! number of channels in the atlas
    INTEGER :: nchan=7
    ! name of the atlas (including version number)
    CHARACTER(len=22) :: name
    ! month of the atlas
    INTEGER :: month                  
    ! resolution of the atlas (equal-area)
    REAL :: dlat=0.25
    ! number of cells per lat band
    INTEGER, POINTER :: ncells(:)
    ! the first cellnumber of lat band
    INTEGER, POINTER :: firstcell(:)
    ! limits of the spatial domain (flaged if global)
    REAL :: lat1, lat2, lon1, lon2 
    ! Emissivities
    REAL, POINTER :: emis(:,:)!emis(ndat,nchan)
    ! Correlations
    REAL, POINTER :: correl(:,:,:)!correl(10,nchan,nchan)
    ! Emissivity uncertainties (std)
    REAL, POINTER :: emis_err(:,:)!emis_err(ndat,nchan)
    ! Surface classe
    INTEGER, POINTER :: class(:)
    ! cellnumber of each of the pixels in the atlas
    INTEGER, POINTER :: cellnum(:)
    ! "Correspondance" vector indicating that for the ith element, the j so that EMIS(j,...) is the emissivity of cellnumber i.
    INTEGER :: correspondance(660066) 
  END TYPE atlas_emis_mw

    type neibor_pix
     REAL :: lat, lon ![deg]
     INTEGER :: cellnum
     REAL :: dist ![km] 
     REAL :: wgt 
    end type neibor_pix

 logical,parameter :: vervose_telsem = .false.

!=======================================================
!ROUTINES ==============================================
!=======================================================
CONTAINS

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE test_inputs(month,lat,lon,theta,freq)
    !======to test if inputs are correct
    INTEGER, INTENT (IN) :: month
    REAL, INTENT (IN)    :: lat, lon, theta, freq
    IF (month<1 .OR. month>12) THEN
       PRINT*,'Pb - month=', month
       STOP
    END IF
    IF (lat<-90 .OR. lat>90) THEN
       PRINT*,'Pb - latitude=',lat
       STOP
    END IF
    IF (lon<0 .OR. lon>360) THEN
       PRINT*,'Pb - longitude=',lon
       STOP
    END IF
    IF (theta<0 .OR. theta>60) THEN
       PRINT*,'Pb - angle=', theta
       STOP
    END IF
    IF (freq<15 .OR. freq>100) THEN
       PRINT*,'Pb - frequency=', freq
       STOP
    END IF
  END SUBROUTINE test_inputs

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE rttov_readmw_atlas(dir,month,atlas,error_status,lat1,lat2,lon1,lon2)
    !======Read a monthly atlas
    INTEGER, INTENT (IN)                :: month
    CHARACTER(len=*), INTENT (IN)       :: dir
    REAL,OPTIONAL, INTENT (IN)          :: lat1, lat2, lon1, lon2 
    TYPE(atlas_emis_mw), INTENT (INOUT) :: atlas
    INTEGER, INTENT (INOUT)             :: error_status !TYPE RTTOV existe?
    !----TRANSITORY VARIABLES
    INTEGER   :: ipos
    INTEGER   :: j,i
    INTEGER   :: iiin=21    ! unit for input
    INTEGER   :: ios
    CHARACTER :: month2*2
    INTEGER   :: cellnum
    REAL      :: ssmi(2*7)
    REAL      :: lat,lon
    INTEGER   :: cur_class
    INTEGER   :: take ! flag to take or not the pixel atlas for location constrains
    !initialisation
    error_status=0
    !
    WRITE(month2,'(I2.2)') month
    atlas%month=month
    atlas%lat1=-777
    atlas%lat2=-777
    atlas%lon1=-777
    atlas%lon2=-777
    IF (PRESENT(lat1)) atlas%lat1=lat1
    IF (PRESENT(lat2)) atlas%lat2=lat2
    IF (PRESENT(lon1)) atlas%lon1=lon1
    IF (PRESENT(lon2)) atlas%lon2=lon2
    !----ALLOCATION SPECIFIC TO SSMI ATLAS
    ALLOCATE(atlas%ncells(int(180./atlas%dlat)))
    ALLOCATE(atlas%firstcell(int(180./atlas%dlat)) )
    CALL equare(atlas%dlat,atlas%ncells,atlas%firstcell)
    atlas%nchan=7
    atlas%name='ssmi_mean_emis_climato'
    !----INITIALISATIONS
    DO j=1,660066
       atlas%correspondance(j)=-777
    END DO
    !----DEFINING NUMBER OF DATA
    if( vervose_telsem ) WRITE(0,*) 'Reading number of data in atlas...'
    OPEN(iiin,file=trim(dir)//atlas%name//'_'//month2//'_cov_interpol',status='old',form='formatted',iostat=ios)

     !print*, trim(dir)//atlas%name//'_'//month2//'_cov_interpol'

    IF( ios /= 0 ) THEN
       WRITE(0,*) 'Error opening the input file ios= ',ios
       !STOP
       error_status=2
    ENDIF
    READ(iiin,*,IOSTAT=ios) j
    atlas%ndat=j
    if( vervose_telsem ) WRITE(0,*) 'Nb data=',atlas%ndat

    !----ALLOCATE VARIABLES
    ALLOCATE(atlas%emis(atlas%ndat,atlas%nchan)) 
    ALLOCATE(atlas%emis_err(atlas%ndat,atlas%nchan)) 
    ALLOCATE(atlas%class(atlas%ndat)) 
    ALLOCATE(atlas%cellnum(atlas%ndat)) 

    ipos=0
    DO j=1,atlas%ndat
       READ(iiin,*) cellnum, (ssmi(i),i=1,2*atlas%nchan),cur_class

       take=1
       IF (PRESENT(lat1)) THEN
          CALL get_coordinates(cellnum,atlas,lat,lon)
          IF ((lat<lat1).OR.(lat>lat2).OR.(lon<lon1).OR.(lon>lon2)) THEN
             take=0
          END IF
       END IF

       IF ((cur_class > 0).AND.(ipos<atlas%ndat).AND.(take==1)) THEN
          ipos=ipos+1
          DO i=1,atlas%nchan
             atlas%emis(ipos,i)=ssmi(i)
             atlas%emis_err(ipos,i)=sqrt(ssmi(atlas%nchan+i))
          END DO
          atlas%cellnum(ipos)=cellnum 
          atlas%class(ipos)=cur_class 
          atlas%correspondance(cellnum)=ipos 
       END IF          
    END DO
    atlas%ndat=ipos;
    CLOSE(iiin)

    !Correlation of uncertainties
    ALLOCATE(atlas%correl(10,atlas%nchan,atlas%nchan)) 
    if( vervose_telsem )  WRITE(0,*) 'reading classes...'
    OPEN(iiin,file=trim(dir)//'correlations', status='old')
    DO i=1,10
       READ(iiin,*) 
       DO j=1,7
          READ(iiin,'(7F5.2)') atlas%correl(i,j,1),atlas%correl(i,j,2), &
               & atlas%correl(i,j,3),atlas%correl(i,j,4),atlas%correl(i,j,5),&
               & atlas%correl(i,j,6),atlas%correl(i,j,7)
       END DO
    END DO
    CLOSE(iiin)

  END SUBROUTINE rttov_readmw_atlas

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE rttov_closemw_atlas(atlas)
    !======free a monthly atlas
    TYPE(atlas_emis_mw), INTENT (INOUT) :: atlas
    DEALLOCATE(atlas%emis) 
    DEALLOCATE(atlas%emis_err) 
    DEALLOCATE(atlas%class) 
    DEALLOCATE(atlas%cellnum) 
    DEALLOCATE(atlas%ncells) 
    DEALLOCATE(atlas%firstcell)
    DEALLOCATE(atlas%correl) 
  END SUBROUTINE rttov_closemw_atlas

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE equare(DLAT,NCELLS,FIRSTCELL)
    !======computes the number of cells in a lattitude band
    !======and the first cellnumber of a lattitude band
    IMPLICIT none
    ! EQUAL-AREA COMPUTATIONS                                           
    ! I  TOTCEL             TOTAL NUMBER OF E.A. BOXES                  
    !  NCELLS(720)        NUMBER OF E.A. BOXES PER LAT ZONE           
    !  TOCELL(1440,720)   BOX NUMBER OF E.A. BOX (1-659600)         
    REAL, INTENT (IN)      :: dlat
    INTEGER,INTENT (INOUT) :: NCELLS(:)
    INTEGER,INTENT (INOUT) :: FIRSTCELL(:)
    INTEGER :: maxlat,maxlon
    INTEGER :: TOTCEL
    INTEGER, ALLOCATABLE   :: TOCELL(:,:)
    INTEGER :: maxlt2,lat,icellr,lat1,lat2,numcel,numcls,lon
    DOUBLE PRECISION    :: rcells,rearth,pi,aezon,hezon,aecell,xlatb
    DOUBLE PRECISION    :: rlatb,rlate,xlate,htb,hte,htzone,rcelat,azone
    !DOUBLE PRECISION    :: dlongr,acell,asq,twopi,halfpi,rcellr
    INTEGER :: I

    maxlat=FLOOR(180./DLAT)
    maxlon=FLOOR(360./DLAT)

    ALLOCATE(TOCELL(maxlon,maxlat))

    !COMMON /EQUCOM/TOCEll,NCELLS
    REARTH = 6371.2d0                                                   
    PI     = 2.0d0 * ASIN(1.d0)   
    !HALFPI=PI/2.d0
    !TWOPI=2.d0*PI                                                       
    RCELAT=(DLAT*PI)/180.d0
    TOTCEL=0d0                                                          
    ! CALCULATE HEIGHT AND AREA OF EQUATORIAL ZONE                          
    HEZON=REARTH*SIN(RCELAT)                                          
    AEZON=2.d0*PI*REARTH*HEZON                                          
    ! CALCULATE AREA OF EQUATORIAL CELL                                    
    AECELL=(AEZON*DLAT)/360.d0                                          
    ! print*,'aecell',aecell
    ! COMPUTE LONGITUDE ZONES FOR EACH LATITUDE ZONE                       
    MAXLT2=MAXLAT/2                                                   
    DO LAT=1,MAXLT2                                               
       XLATB=(LAT-1)*DLAT                                                
       XLATE=XLATB+DLAT                                                  
       RLATB=(2.0d0*PI*XLATB)/360.0d0
       RLATE=(2.0d0*PI*XLATE)/360.0d0
       !CALCULATE HEIGHTS OF LATB,LATE,ZONE                                  
       HTB=REARTH*SIN(RLATB)                                             
       HTE=REARTH*SIN(RLATE)                                             
       HTZONE=HTE-HTB                                                    
       AZONE=2.d0*PI*REARTH*HTZONE                                         
       !CALCULATE NUMBER OF CELLS 
       RCELLS=AZONE/AECELL                                               
       ICELLR=FLOOR(RCELLS+.50d0)                                                
       !AUGMENT TOTAL # GRID CELLS (BOTH HEMISPHERES)                     
       TOTCEL=TOTCEL+2*ICELLR                                            
       !RCELLR=ICELLR                                                     
       !DLONGR=360.0d0/RCELLR                                                
       !ACELL=AZONE/RCELLR                                                
       !ASQ=AZONE/MAXLON                                                  
       !CREATE TABLE OF LONGITUDES
       LAT1=LAT+MAXLT2                                                   
       LAT2=MAXLT2+1-LAT                                                 
       NCELLS(LAT1)=ICELLR                                               
       NCELLS(LAT2)=ICELLR                                               
    END DO
    NUMCEL = 0                                                     
    ! THROUGH EACH LAT ZONE                                        
    DO LAT=1,MAXLAT                                              
       NUMCLS = NCELLS(LAT)                                           
       ! if (lat.eq.9) print*,lat,numcls,360./numcls
       ! LOOP THROUGH EACH LON FOR THIS LAT ZONE                           
       DO LON=1,NUMCLS                                           
          NUMCEL = NUMCEL + 1                                         
          ! FILL TOCELL ARRAY WITH STARTING CELL NUMBER
          TOCELL(LON,LAT) = NUMCEL
       END DO
    END DO
    DEALLOCATE(TOCELL)

    ! search for the first cellnumber in each lat band
    FIRSTCELL(1)=1
    DO I=2,maxlat
       FIRSTCELL(I)=FIRSTCELL(I-1)+NCELLS(I-1)
    END DO
  END SUBROUTINE equare

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM


  !------------------------------------------------------------------
  FUNCTION calc_cellnum(lat,lon,atlas)
    !======computes the cellnumber given the lat and long
    !======using the ncells included in the atlas
    IMPLICIT none
    REAL, INTENT (IN) :: lat, lon
    type(atlas_emis_mw), INTENT (IN) :: atlas
    INTEGER :: calc_cellnum
    INTEGER :: ilat,ilon, i
    ! search for the cellnum in the older file
    calc_cellnum=0

!    ilat=int((lat+90.)/atlas%dlat)+1
!    ilon=int(lon/(360./atlas%ncells(ilat)))+1
! bug.. int must be nint
    ilat=nint((lat+90.)/atlas%dlat)+1
    ilon=nint(lon/(360./atlas%ncells(ilat)))+1


    do i=1,ilat-1
       calc_cellnum=atlas%ncells(i)+calc_cellnum
    end do
    calc_cellnum=calc_cellnum+ilon
  END FUNCTION calc_cellnum

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE get_coordinates(cellnum,atlas,lat,lon)
    !======computes lat and lon given the cellnumber
    IMPLICIT none
    INTEGER, INTENT (IN) :: cellnum
    type(atlas_emis_mw), INTENT (IN) :: atlas
    REAL, INTENT (OUT) :: lat, lon
    !INTEGER :: ilat,ilon
    INTEGER :: i
    REAL :: res_lat ! latitude resolution
    INTEGER :: index_lat_max,index_lat,index_lon

    res_lat = atlas%dlat 

    index_lat_max  = int(180/res_lat)
 
    IF (cellnum >= atlas%firstcell(index_lat_max)) THEN
       index_lat = index_lat_max
       lat =(index_lat - 0.5)*res_lat - 90
       index_lon = cellnum - atlas%firstcell(index_lat_max) 
       lon = (index_lon - 0.5)*(360.0/atlas%ncells(index_lat))
    ELSE
       DO i=1,index_lat_max-1 
          IF ( (cellnum>=atlas%firstcell(i)) .AND. (cellnum<atlas%firstcell(i+1)) ) THEN
	     index_lat = i
	     lat = (index_lat - 0.5)*res_lat- 90
	     !	cout << i << "  " << *lat <<endl;
	     index_lon = cellnum - atlas%firstcell(i)
	     lon = (index_lon - 0.5)*(360.0/atlas%ncells(index_lat))
          END IF
       END DO
    END IF
  END SUBROUTINE get_coordinates

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE calc_cellnum_mult(lat,lon,resol,atlas,cell_num_mult,nb_cell)
    !======computes the cellnumbers given the lat, long and resol
    !======using the ncells and firstcells included in the atlas
    IMPLICIT none
    REAL, INTENT (IN) :: lat, lon
    ! lat(-90 90)  lon(-180 180) of the pixel center
    ! resol = spatial resolution of the new grid
    type(atlas_emis_mw), INTENT (IN) :: atlas
    INTEGER, POINTER :: cell_num_mult(:)
    INTEGER, INTENT (OUT) :: nb_cell
    REAL, INTENT(IN)    :: resol
    INTEGER :: ilat,ilon, i
    INTEGER :: calc_cellnum
    !--------new variables
    INTEGER :: nbcel,cell(200) 
    !                 200 is here the maximum number of cells inside the "big" pixel
    INTEGER :: i2lon,i3lon,i2lat,nbreslat,nbreslon,i4lon
    REAL    :: maxlat
    ! search for the cellnum in the older file
    calc_cellnum=0
    ilat=INT((lat+90.)/atlas%dlat)+1
    ilon=INT(lon/(360./atlas%ncells(ilat)))+1
    DO i=1,ilat-1
       calc_cellnum=atlas%ncells(i)+calc_cellnum
    END DO
    calc_cellnum=calc_cellnum+ilon

    ! search for the nbcel pixels that are in the box of resolution 
    ! 'resol' around lat and lon.         
    maxlat=180./atlas%dlat
    !maxlon=360./atlas%dlat

    nbcel=1
    cell(nbcel)=calc_cellnum
    nbreslat=INT(resol/atlas%dlat/2.)
    IF (nbreslat>=1) THEN

       DO i2lat=ilat-nbreslat,ilat+nbreslat
          IF ( (i2lat<1).OR.(i2lat>INT(maxlat)) ) THEN
             !PRINT*,'Congratulation, you reached the pole!' 
          ELSE
             IF (ABS((i2lat-.5)*atlas%dlat-90-lat)<=resol/2) THEN
                i2lon=INT(lon/(360./atlas%ncells(i2lat)))+1
                nbreslon=INT(resol/(360./(1.*atlas%ncells(i2lat)))/2.)
                DO i3lon=i2lon-nbreslon,i2lon+nbreslon
                   IF (amod(ABS((i3lon-.5)*(360./atlas%ncells(i2lat))-lon),360.)<=resol/2.) THEN 
                      nbcel=nbcel+1
                      i4lon=i3lon
                      IF (i3lon<1.) i4lon=atlas%ncells(i2lat)+i3lon
                      IF (i3lon>=atlas%ncells(i2lat)) i4lon=i3lon-atlas%ncells(i2lat)
                      cell(nbcel)=atlas%firstcell(i2lat)+ i4lon-1
                      IF (cell(nbcel)==cell(1)) nbcel=nbcel-1
                   END IF
                END DO
             END IF
          END IF
       END DO

    END IF
    nb_cell=nbcel
    ALLOCATE(cell_num_mult(int(nb_cell)))
    DO i=1,nb_cell
       cell_num_mult(i)=cell(i)
    END DO
  END SUBROUTINE calc_cellnum_mult

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE emis_interp_ind_sing(lat,lon,theta,freq,atlas,ev,eh,stdv,stdh,verb)
    !======interpolates emissivities 
    !          IND: individual cellnumber atlas
    !          SING: singular channel
    Implicit none     
    REAL, INTENT (IN)    :: lat, lon, theta, freq
    TYPE(atlas_emis_mw), INTENT (IN) :: atlas
    REAL, INTENT (OUT)   :: ev, eh   !resultat de l'interpolation
    REAL, OPTIONAL, INTENT (OUT)   :: stdv, stdh
    INTEGER, INTENT (IN) :: verb
    INTEGER :: ipos
    INTEGER :: cellnum
    REAL    :: ev_a(3),eh_a(3)      !emissivity in the atlas
    REAL    :: stdv_a(3),stdh_a(3)  !std in the atlas

    !Initialisations
    ev=0
    eh=0
    IF (PRESENT(stdv)) stdv=0
    IF (PRESENT(stdh)) stdh=0

    cellnum=calc_cellnum(lat,lon,atlas)
    ipos=atlas%correspondance(cellnum)
    IF (ipos>0) THEN
       ev_a(1)=atlas%emis(ipos,1)
       eh_a(1)=atlas%emis(ipos,2)
       ev_a(2)=atlas%emis(ipos,4)
       eh_a(2)=atlas%emis(ipos,5)
       ev_a(3)=atlas%emis(ipos,6)
       eh_a(3)=atlas%emis(ipos,7)
       CALL emis_interp(lat,lon,theta,freq,atlas%class(ipos),ev_a,eh_a,ev,eh)
       IF (PRESENT(stdv)) THEN
          stdv_a(1)=atlas%emis_err(ipos,1)
          stdv_a(2)=atlas%emis_err(ipos,4)
          stdv_a(3)=atlas%emis_err(ipos,6)
          CALL interp_freq2(stdv_a(1),stdv_a(2),stdv_a(3),freq,stdv)
       END IF
       IF (PRESENT(stdh)) THEN
          stdh_a(1)=atlas%emis_err(ipos,2)
          stdh_a(2)=atlas%emis_err(ipos,5)
          stdh_a(3)=atlas%emis_err(ipos,7)
          CALL interp_freq2(stdh_a(1),stdh_a(2),stdh_a(3),freq,stdh)
       END IF
    END IF
    IF (verb==1) THEN
       WRITE(0,*) 'Cellnum=',cellnum,' lat=',lat,' lon=',lon,' ipos=',ipos
    END IF
  END SUBROUTINE emis_interp_ind_sing

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE emis_interp_ind_mult(lat,lon,theta,freq,n_chan,atlas,ev,eh,std,verb)
    !======interpolates emissivities 
    !          IND: individual cellnumber atlas
    !          MULT: multiple channels
    IMPLICIT none     
    INTEGER, INTENT (IN) :: n_chan
    REAL, INTENT (IN)    :: lat, lon, theta, freq(:)
    TYPE(atlas_emis_mw), INTENT (IN) :: atlas
    REAL, INTENT (OUT)   :: ev(:), eh(:)   !resultat de l'interpolation
    REAL, OPTIONAL, INTENT (OUT)   :: std(:,:)
    INTEGER, INTENT (IN) :: verb

    INTEGER :: ipos, i, j, i2, j2
    INTEGER :: cellnum
    REAL    :: ev_a(3),eh_a(3)      !emissivity in the atlas
    REAL    :: stdv_a(3),stdh_a(3)  !std in the atlas
    REAL    :: std2
    REAL    :: FIM(2*n_chan,3*2)    !Frequency Interpolation Matrix
    REAL    :: trans_std(3*2,2*n_chan)        !transpose Frequency Interpolation Matrix
    REAL    :: a,b,c                !frequency interpolation coefficients
    REAL    :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas

    REAL    :: new_FIM(3*2,2*n_chan)          !transpose Frequency Interpolation Matrix
    REAL    :: correlation(2*n_chan,2*n_chan) !correlation

    !initialisations    
    DO i=1,n_chan
       ev(i)=0
       eh(i)=0
    END DO
    IF (PRESENT(std)) THEN
       DO i=1,2*n_chan
          DO j=1,2*n_chan
             std(i,j)=0
          END DO
       END DO
    END IF

    !localization
    cellnum=calc_cellnum(lat,lon,atlas)
    ipos=atlas%correspondance(cellnum)

    !compute emissivities
    IF (ipos>0) THEN
       ev_a(1)=atlas%emis(ipos,1)
       eh_a(1)=atlas%emis(ipos,2)
       ev_a(2)=atlas%emis(ipos,4)
       eh_a(2)=atlas%emis(ipos,5)
       ev_a(3)=atlas%emis(ipos,6)
       eh_a(3)=atlas%emis(ipos,7)
       stdv_a(1)=atlas%emis_err(ipos,1)  !don't matter to put
       stdh_a(1)=atlas%emis_err(ipos,2)  !real data here, we 
       stdv_a(2)=atlas%emis_err(ipos,4)  !just need a, b, c  
       stdh_a(2)=atlas%emis_err(ipos,5)  !don't matter to put
       stdv_a(3)=atlas%emis_err(ipos,6)  !real data here, we 
       stdh_a(3)=atlas%emis_err(ipos,7)  !just need a, b, c  

       DO i=1,n_chan
          !verifier que j'ai ici la bonne interpretation des canaux !!!!
          CALL emis_interp(lat,lon,theta,freq(i),atlas%class(ipos),ev_a,eh_a,ev(i),eh(i))
       END DO

       IF (PRESENT(std)) THEN !compute uncertainties on emissivities
          !compute covariance matrix on the atlas frequencies
          DO i=1,6
             IF (i>2) THEN
                i2=i+1
             ELSE
                i2=i
             END IF
             DO j=1,6
                IF (j>2) THEN
                   j2=j+1
                ELSE
                   j2=j
                END IF
                cov(i,j)=atlas%correl(atlas%class(ipos),i2,j2)* &
                    & SQRT(atlas%emis_err(ipos,i2)*atlas%emis_err(ipos,j2))
             END DO
          END DO
          DO i=1,6
             DO j=1,6
                std(i,j)=cov(i,j)/SQRT(cov(i,i)*cov(j,j))
             END DO
          END DO
          !compute the Frequency Linear Matrix
          DO i=1,n_chan
             CALL interp_freq2(stdv_a(1),stdv_a(2),stdv_a(3),freq(i),std2,a,b,c)
             FIM(i,1)=a
             FIM(i,2)=b
             FIM(i,3)=c
          END DO
          DO i=1,n_chan
             DO j=1,3
                FIM(i,3+j)=0
                FIM(n_chan+i,j)=0
                FIM(n_chan+i,j+3)=FIM(i,j)
             END DO
          END DO
          new_FIM=TRANSPOSE(FIM)
          trans_std=MATMUL(cov,new_FIM)
          !std=MATMUL(FIM,MATMUL(cov,new_FIM)) 
          std=MATMUL(FIM,trans_std) 
       END IF
       !verbose ?
       IF (verb==1) THEN
          WRITE(0,*) 'Cellnum=',cellnum,' lat=',lat,' lon=',lon,' ipos=',ipos
          WRITE(0,*)    'emis_SSMI_V(1)=',ev_a(1),'/',stdv_a(1)
          WRITE(0,*)    'emis_SSMI_V(2)=',ev_a(2),'/',stdv_a(2)
          WRITE(0,*)    'emis_SSMI_V(3)=',ev_a(3),'/',stdv_a(3)
          WRITE(0,*)    'emis_SSMI_H(1)=',eh_a(1),'/',stdh_a(1)
          WRITE(0,*)    'emis_SSMI_H(2)=',eh_a(2),'/',stdh_a(2)
          WRITE(0,*)    'emis_SSMI_H(3)=',eh_a(3),'/',stdh_a(3)
          DO i=1,n_chan
             WRITE(0,*) 'IND(',freq(i),') EV=', ev(i),' EH=', eh(i)
          END DO
          WRITE(0,*) 'Matrice de covariance 6x6'
          DO i=1,6
             WRITE(0,'(I3,a,6F8.5)') i,': ',(cov(i,j),j=1,6)
          END DO
          WRITE(0,*) 'Matrice de correl'
          DO i=1,6
             WRITE(0,'(I3,a,6F8.5)') i,': ',(std(i,j),j=1,6)
          END DO
          !the Frequency Linear Matrix
          DO i=1,n_chan
             IF (PRESENT(std)) THEN
                WRITE(0,'(3F8.5)') FIM(i,1),FIM(i,2),FIM(i,3)
             END IF
          END DO
          WRITE(0,*) 'FIM'
          DO i=1,2*n_chan
             WRITE(0,'(I3,a,6F8.5)') i,': ',(FIM(i,j),j=1,6)
          END DO
          WRITE(0,*) 'STD'
          DO i=1,2*n_chan
             WRITE(0,'(I3,a,10F8.5)') i,': ',(std(i,j),j=1,2*n_chan)
          END DO
          WRITE(0,*) 'Matrice de correl'
          DO i=1,2*n_chan
             DO j=1,2*n_chan
                correlation(i,j)=std(i,j)/SQRT(std(i,i)*std(j,j))
             END DO
             WRITE(0,'(I3,a,10F8.5)') i,': ',(correlation(i,j),j=1,2*n_chan)
          END DO
       END IF  !verb
    END IF   !ipos
  END SUBROUTINE emis_interp_ind_mult

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE emis_interp_int_sing(lat,lon,resol,theta,freq,atlas,ev,eh,stdv,stdh,verb)
    !======interpolates emissivities 
    !          INT: integrate atlas, i.e. multiple cellnumber atlas
    !          SING: singular channel
    IMPLICIT none     
    REAL, INTENT (IN)                :: lat, lon, theta, freq
    TYPE(atlas_emis_mw), INTENT (IN) :: atlas
    REAL, INTENT (OUT)               :: ev, eh
    REAL, OPTIONAL, INTENT (OUT)     :: stdv,stdh
    INTEGER, INTENT (IN)             :: verb
    REAL, INTENT (IN)                :: resol

    INTEGER :: ipos
    REAL    :: ev_a(3),eh_a(3)     !emissivities in the atlas
    REAL    :: stdv_a(3),stdh_a(3) !std emissivities in the atlas
    INTEGER :: ii
    INTEGER :: nb_cell
    INTEGER, POINTER :: cellnum_mult(:) => NULL()
    REAL    :: ev_mean, eh_mean
    REAL    :: stdv_mean, stdh_mean
    INTEGER :: inumb

    !initialisations
    ev=0
    eh=0
    IF (PRESENT(stdv)) stdv=0
    IF (PRESENT(stdh)) stdh=0
    !computes the list of cells that need to be integrated in the atlas
    CALL calc_cellnum_mult(lat,lon,resol,atlas,cellnum_mult,nb_cell)
    ev_mean=0
    eh_mean=0
    stdv_mean=0
    stdh_mean=0
    inumb=0
    DO ii=1,nb_cell
       ipos=atlas%correspondance(cellnum_mult(ii))
       IF (ipos>0) THEN
          inumb=inumb+1
          ev_a(1)=atlas%emis(ipos,1)
          eh_a(1)=atlas%emis(ipos,2)
          ev_a(2)=atlas%emis(ipos,4)
          eh_a(2)=atlas%emis(ipos,5)
          ev_a(3)=atlas%emis(ipos,6)
          eh_a(3)=atlas%emis(ipos,7)
          ev=0
          eh=0
          CALL emis_interp(lat,lon,theta,freq,atlas%class(ipos),ev_a,eh_a,ev,eh)
          !WRITE(0,*) 'Cellnum(',ii,')=',cellnum_mult(ii),' ',atlas%correspondance(cellnum_mult(ii)),' ',ev,' ',eh
          ev_mean=ev_mean+ev
          eh_mean=eh_mean+eh
          IF (PRESENT(stdv)) THEN
             stdv_a(1)=atlas%emis_err(ipos,1)
             stdv_a(2)=atlas%emis_err(ipos,4)
             stdv_a(3)=atlas%emis_err(ipos,6)
             CALL interp_freq2(stdv_a(1),stdv_a(2),stdv_a(3),freq,stdv)
             stdv_mean=stdv_mean+stdv
          END IF
          IF (PRESENT(stdh)) THEN
             stdh_a(1)=atlas%emis_err(ipos,2)
             stdh_a(2)=atlas%emis_err(ipos,5)
             stdh_a(3)=atlas%emis_err(ipos,7)
             CALL interp_freq2(stdh_a(1),stdh_a(2),stdh_a(3),freq,stdh)
             stdh_mean=stdh_mean+stdh
          END IF
       END IF
    END DO
    IF (inumb>0) THEN
       ev=ev_mean/inumb
       eh=eh_mean/inumb
       IF (PRESENT(stdv)) THEN
          stdv=stdv_mean/inumb
       ENDIF
       IF (PRESENT(stdh)) THEN
          stdh=stdh_mean/inumb
       ENDIF
    END IF
    DEALLOCATE(cellnum_mult) 
  END SUBROUTINE emis_interp_int_sing

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  !------------------------------------------------------------------
  SUBROUTINE emis_interp_int_mult(lat,lon,resol,theta,freq,n_chan,atlas,ev,eh,std,verb)
    !======interpolates emissivities 
    !          INT: integrate atlas, i.e. multiple cellnumber atlas
    !          MULT: multiple channel
    IMPLICIT none     
    INTEGER, INTENT (IN)             :: n_chan
    REAL, INTENT (IN)                :: lat, lon, theta, freq(:)
    TYPE(atlas_emis_mw), INTENT (IN) :: atlas
    REAL, INTENT (OUT)               :: ev(:), eh(:)
    REAL, OPTIONAL, INTENT (OUT)     :: std(:,:)
    INTEGER, INTENT (IN)             :: verb
    REAL, INTENT(IN)                 :: resol
    INTEGER :: ipos
    REAL    :: ev_a(3),eh_a(3)     !emissivities in the atlas
    REAL    :: stdv_a(3) !std emissivities in the atlas
    INTEGER :: ii, i2, j2
    INTEGER :: nb_cell
    INTEGER, POINTER :: cellnum_mult(:) => NULL()
    REAL    :: ev_mean(10), eh_mean(10)
    REAL    :: std_mean(10,10)
    REAL    :: std2(10,10)
    REAL    :: ev2, eh2
    INTEGER :: inumb, i,j
  
    REAL    :: FIM(2*n_chan,3*2)    !Frequency Interpolation Matrix
    REAL    :: a,b,c                !frequency interpolation coefficients
    REAL    :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas
    REAL    :: new_FIM(3*2,2*n_chan)          !transpose Frequency Interpolation Matrix
    REAL    :: trans_std(3*2,2*n_chan)        !transpose Frequency Interpolation Matrix
    !REAL    :: correlation(2*n_chan,2*n_chan) !correlation
  
    !initialisations    
    ev(:)=0
    eh(:)=0
    ev_mean(:)=0
    eh_mean(:)=0
    IF (PRESENT(std)) THEN
       std_mean(:,:)=0
       std(:,:)=0
    END IF

    !computes the list of cells that need to be integrated in the atlas
    CALL calc_cellnum_mult(lat,lon,resol,atlas,cellnum_mult,nb_cell)
  
    inumb=0
    DO ii=1,nb_cell
       ipos=atlas%correspondance(cellnum_mult(ii))
       IF (ipos>0) THEN
          inumb=inumb+1
          ev_a(1)=atlas%emis(ipos,1)
          eh_a(1)=atlas%emis(ipos,2)
          ev_a(2)=atlas%emis(ipos,4)
          eh_a(2)=atlas%emis(ipos,5)
          ev_a(3)=atlas%emis(ipos,6)
          eh_a(3)=atlas%emis(ipos,7)
          stdv_a(1)=atlas%emis_err(ipos,1)  !don't matter to put
          !stdh_a(1)=atlas%emis_err(ipos,2)  !real data here, we 
          stdv_a(2)=atlas%emis_err(ipos,4)  !just need a, b, c  
          !stdh_a(2)=atlas%emis_err(ipos,5)  !don't matter to put
          stdv_a(3)=atlas%emis_err(ipos,6)  !real data here, we 
          !stdh_a(3)=atlas%emis_err(ipos,7)  !just need a, b, c  
          DO i=1,n_chan
             !initialisations    
             ev2=0
             eh2=0
             CALL emis_interp(lat,lon,theta,freq(i),atlas%class(ipos),ev_a,eh_a,ev2,eh2)
             !WRITE(0,*) 'Cellnum(',ii,')=',cellnum_mult(ii),' ',atlas%correspondance(cellnum_mult(ii)),' ',ev2,' ',eh2
             ev_mean(i)=ev_mean(i)+ev2
             eh_mean(i)=eh_mean(i)+eh2
          END DO

          IF (PRESENT(std)) THEN !compute uncertainties on emissivities
             !compute covariance matrix on the atlas frequencies
             DO i=1,6
                IF (i>2) THEN
                   i2=i+1
                ELSE
                   i2=i
                END IF
                DO j=1,6
                   IF (j>2) THEN
                      j2=j+1
                   ELSE
                      j2=j
                   END IF
                   cov(i,j)=atlas%correl(atlas%class(ipos),i2,j2)* &
                       & SQRT(atlas%emis_err(ipos,i2)*atlas%emis_err(ipos,j2))
                END DO
             END DO
             DO i=1,6
                DO j=1,6
                   std2(i,j)=cov(i,j)/SQRT(cov(i,i)*cov(j,j))
                END DO
             END DO
             !compute the Frequency Linear Matrix
             DO i=1,n_chan
                CALL interp_freq2(stdv_a(1),stdv_a(2),stdv_a(3),freq(i),ev2,a,b,c) !ev2=dummy variable
                FIM(i,1)=a
                FIM(i,2)=b
                FIM(i,3)=c
             END DO
             DO i=1,n_chan
                DO j=1,3
                   FIM(i,3+j)=0
                   FIM(n_chan+i,j)=0
                   FIM(n_chan+i,j+3)=FIM(i,j)
                END DO
             END DO
             new_FIM=TRANSPOSE(FIM)
             trans_std=MATMUL(cov,new_FIM)
             !std2=MATMUL(FIM,MATMUL(cov,new_FIM)) 
             std2=MATMUL(FIM,trans_std) 

             Do i=1,2*n_chan
                DO j=1,2*n_chan
                   std(i,j)=std_mean(i,j)+std2(i,j)
                END DO
             END DO
          END IF !std ou pas
       END IF !ipos
    END DO  !nb_cell

    IF (inumb>0) THEN
       DO i=1,n_chan
          ev(i)=ev_mean(i)/inumb
          eh(i)=eh_mean(i)/inumb
          !WRITE(0,*) 'INT: EV=',ev,' EH=',eh
       END DO
       IF (PRESENT(std)) THEN
          DO i=1,2*n_chan
             DO j=1,2*n_chan
                std_mean(i,j)=std_mean(i,j)/inumb
             END DO
          END DO
       END IF
    END IF

    DEALLOCATE(cellnum_mult) 
  END SUBROUTINE emis_interp_int_mult

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM
 
  !------------------------------------------------------------------
  SUBROUTINE interp_freq2(emiss19,emiss37,emiss85,f,emiss,an,bn,cn)  
    !======linear interpolation of emissivity given the freq
    !======and the atlas values on that cellnum
    IMPLICIT none
    REAL, INTENT (IN) :: emiss19,emiss37,emiss85,f
    REAL, INTENT (OUT) :: emiss
    REAL, OPTIONAL, INTENT (OUT) :: an, bn, cn
    REAL :: a, b, c
    IF (f<=19.35) THEN
       a=1
       b=0
       c=0
       emiss = emiss19
    ELSE IF ((19.35<f).AND.(f<=37.)) THEN
       a=(37.-f  )/(37.-19.35)
       b=(f-19.35)/(37.-19.35)
       c=0
       emiss = a*emiss19+b*emiss37
    ELSE IF ((f>37.).AND.(f<85.5)) THEN
       a=0
       b=(85.5-f )/(85.5-37)
       c=(f-37   )/(85.5-37)
       emiss = b*emiss37+c*emiss85
    ELSE IF(85.5<=f) THEN
       a=0
       b=0
       c=1
       emiss = emiss85
    END IF
    IF (PRESENT(an)) an=a
    IF (PRESENT(bn)) bn=b
    IF (PRESENT(cn)) cn=c
  END SUBROUTINE interp_freq2

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

 !------------------------------------------------------------------
   SUBROUTINE emis_interp(lat,lon,theta,freq,classe,ev,eh,emiss_interp_v,emiss_interp_h)
    !======interpolation of emissivities for angle, freq
    IMPLICIT none     
    INTEGER, INTENT (IN) :: classe
    REAL, INTENT (IN)    :: lat, lon, theta, freq, ev(3), eh(3)
    REAL, INTENT (OUT)   :: emiss_interp_v, emiss_interp_h
    REAL    :: e0, theta0, theta53 , emiss_scal_v(3), emiss_scal_h(3)
    REAL    :: S1_v, S1_h, S2_v, S2_h, S_v, S_h, a0, a1, a2, a3, b0, b1
    REAL    :: b2, b3, em53_v, em53_h, emtheta_v, emtheta_h
    REAL    :: a0_k0(3,10),a0_k1(3,10),a0_k2(3,10)
    REAL    :: a0_eveh(3,10),a1_eveh(3,10),a2_eveh(3,10),a3_eveh(3,10)
    REAL    :: b0_eveh(3,10),b1_eveh(3,10),b2_eveh(3,10),b3_eveh(3,10)
    INTEGER :: j
    !COMMON /EMISSIVITE/emiss_interp_v,emiss_interp_h
    data a0_k0/0.11509,0.091535,0.34796,0.10525,0.16627,0.24434, &
         & 0.29217,0.23809,0.28954,0.17516,0.19459,0.28697, &
         & 0.10521,0.12126,0.30278,0.18212,0.19625,0.14551, &
         & -0.19202,0.5411,0.03739,0.10292,0.5486,-0.058937, &
         & -0.022672,0.44492,-0.058448,-0.33894,-0.17621,0.14742/
    data a0_k1/0.61168,0.59095,0.7918,0.60271,0.69213,0.62218, &
         &  0.32728,0.34334,0.37062,0.51217,0.4491,0.50101, &
         & 0.48913,0.41932,0.29734,0.64474,0.30637,0.031107, &
         & 1.0405,0.17538,1.3215,0.61819,0.31298,1.7218, &
         & 0.87761,0.47583,1.2583,1.0959,0.92842,0.51033/
    data a0_k2/0.26726,0.32033,-0.14778,0.28547,0.13592,0.13193, &
         & 0.37178,0.41813,0.33875,0.30203,0.35479,0.20189, &
         & 0.40663,0.47493,0.40668,0.14811,0.52382,0.86634, &
         & 0.14286,0.27164,-0.37947,0.2737,0.12001,-0.67315, &
         & 0.13492,0.065463,-0.19316,0.24905,0.25475,0.34637/
    data a0_eveh/0.9592599869E+00,0.9565299749E+00,0.9511899948E+00, &
         & 0.9560700059E+00,0.9541199803E+00,0.9483199716E+00, &
         & 0.9461100101E+00,0.9439799786E+00,0.9387800097E+00, &
         & 0.9317600131E+00,0.9289000034E+00,0.9236800075E+00, &
         & 0.9208700061E+00,0.9190599918E+00,0.9105200171E+00, &
         & 0.9162799716E+00,0.8937299848E+00,0.8014699817E+00, &
         & 0.9570500255E+00,0.9213600159E+00,0.7893999815E+00, &
         & 0.9639400244E+00,0.9530599713E+00,0.8850200176E+00, &
         & 0.9685299993E+00,0.9622600079E+00,0.9118800163E+00, &
         & 0.8997200131E+00,0.9012699723E+00,0.9107499719E+00/
    data a1_eveh/0.3627802414E-07,-0.7778328204E-08,0.4396108011E-07, &
         & 0.2503205394E-06,0.1996262995E-06,0.2929977541E-06, &
         & 0.4190530660E-06,0.3655744649E-06,0.3519195673E-06, &
         & 0.5574374313E-06,0.5273076340E-06,0.5376484182E-06, &
         & 0.1026844529E-05,0.9679998811E-06,0.8616486866E-06, &
         & 0.3180800832E-06,0.2886778532E-06,0.2310362675E-06, &
         & -0.1118036366E-06,-0.1502856577E-06,0.4842232926E-07, &
         & -0.8410978580E-08,-0.3478669441E-07,0.2209441590E-06, &
         & 0.2485776633E-06,0.1800235907E-06,0.2510202251E-06, &
         & 0.2687000915E-06,0.1740325644E-06,0.3562134339E-06/
    data a2_eveh/0.3067140824E-05,0.2520012231E-05,0.4831396382E-05, &
         & 0.8213598448E-05,0.7378375358E-05,0.1022081960E-04, &
         & 0.1225889173E-04,0.1165553113E-04,0.1188659007E-04, &
         & 0.1693615741E-04,0.1648317448E-04,0.1715818144E-04, &
         & 0.2744720041E-04,0.2642072104E-04,0.2671847506E-04, &
         & 0.1349592094E-04,0.1261523357E-04,0.5447756394E-05, &
         & 0.2064244654E-05,0.1919016057E-06,0.5940860319E-06, &
         & 0.5334760772E-05,0.4130339221E-05,0.4104662821E-05, &
         & 0.6530796327E-05,0.5727014013E-05,0.7451782039E-05, &
         & 0.1071246970E-04,0.9539280654E-05,0.1034286015E-04/
    data a3_eveh/-0.2004991551E-07,-0.6895366056E-07, &
         & -0.2047409282E-06, &
         & -0.7322448425E-07,-0.1273002681E-06,-0.2729916844E-06, &
         & -0.9421125213E-07,-0.1683332300E-06,-0.2726891637E-06, &
         & -0.1317753799E-06,-0.2107972250E-06,-0.3556060904E-06, &
         & -0.1889465580E-06,-0.2757958271E-06,-0.4909850304E-06, &
         & 0.7339644004E-08,-0.4058669560E-06,-0.4146343997E-06, &
         & 0.6170279931E-07,-0.1998567996E-06,-0.4713119139E-07, &
         & -0.1361754887E-07,-0.1765622955E-06,-0.2348146637E-06, &
         & -0.3901189061E-07,-0.1305666189E-06,-0.1533838798E-06, &
         & -0.2679148992E-07,-0.4441960044E-07,-0.1815613899E-06/
    data b0_eveh/0.9592599869E+00,0.9565299749E+00,0.9511899948E+00, &
         & 0.9560700059E+00,0.9541199803E+00,0.9483199716E+00, &
         & 0.9461100101E+00,0.9439799786E+00,0.9387800097E+00, &
         & 0.9317600131E+00,0.9289000034E+00,0.9236800075E+00, &
         & 0.9208700061E+00,0.9190599918E+00,0.9105200171E+00, &
         & 0.9162799716E+00,0.8937299848E+00,0.8014699817E+00, &
         & 0.9570500255E+00,0.9213600159E+00,0.7893999815E+00, &
         & 0.9639400244E+00,0.9530599713E+00,0.8850200176E+00, &
         & 0.9685299993E+00,0.9622600079E+00,0.9118800163E+00, &
         & 0.8997200131E+00,0.9012699723E+00,0.9107499719E+00/
    data b1_eveh/0.3626608347E-07,-0.7786279177E-08,0.4393379172E-07, &
         & 0.2502746099E-06,0.1995944388E-06,0.2929554341E-06, &
         & 0.4189516289E-06,0.3655020180E-06,0.3518483140E-06, &
         & 0.5572838404E-06,0.5271903092E-06,0.5375342766E-06, &
         & 0.1026605219E-05,0.9677979733E-06,0.8614680951E-06, &
         & 0.3179358714E-06,0.2884899004E-06,0.2308632219E-06, &
         & -0.1118781370E-06,-0.1503948681E-06,0.4834672396E-07, &
         & -0.8455684153E-08,-0.3485171618E-07,0.2208606134E-06, &
         & 0.2485595019E-06,0.1799959364E-06,0.2509846695E-06, &
         & 0.2686167306E-06,0.1739760478E-06,0.3561317214E-06/
    data b2_eveh/0.3065537157E-05,0.2518960400E-05,0.4829731552E-05, &
         & 0.8209894986E-05,0.7375769655E-05,0.1021809931E-04, &
         & 0.1225203869E-04,0.1165053800E-04,0.1188218721E-04, &
         & 0.1692612022E-04,0.1647546378E-04,0.1715117833E-04, &
         & 0.2743142431E-04,0.2640772436E-04,0.2670711910E-04, &
         & 0.1348545720E-04,0.1260529825E-04,0.5439695997E-05, &
         & 0.2058213340E-05,0.1860650656E-06,0.5898303925E-06, &
         & 0.5330772183E-05,0.4126528893E-05,0.4100859314E-05, &
         & 0.6528573977E-05,0.5725009032E-05,0.7449450095E-05, &
         & 0.1070590315E-04,0.9534271157E-05,0.1033751869E-04/
    data b3_eveh/-0.1370247134E-06,-0.1436897747E-06, &
         & -0.2954870411E-06, &
         & -0.3118435643E-06,-0.2916583242E-06,-0.4311032171E-06, &
         & -0.5048401022E-06,-0.4662823869E-06,-0.5206445053E-06, &
         & -0.7210980471E-06,-0.6662896794E-06,-0.7548637200E-06, &
         & -0.1110204039E-05,-0.1030801400E-05,-0.1140921199E-05, &
         & -0.6330818110E-06,-0.9186441048E-06,-0.7947813856E-06, &
         & -0.3242539890E-06,-0.5027602583E-06,-0.2777987334E-06, &
         & -0.2747250676E-06,-0.3811997260E-06,-0.4102405455E-06, &
         & -0.1994112324E-06,-0.2555484855E-06,-0.2842682534E-06, &
         & -0.4413041665E-06,-0.3717419474E-06,-0.4975536854E-06/
    ! Interpolation en angle
    DO j = 1,3
       ! Calcul par regression multilineaire de la valeur e0 en theta=0
       e0 = a0_k0(j,classe)+a0_k1(j,classe)*ev(j)+a0_k2(j,classe)*eh(j)
       ! Lecture des coefficients des polynomes ev et eh
       a0 = a0_eveh(j,classe)
       a1 = a1_eveh(j,classe)
       a2 = a2_eveh(j,classe)
       a3 = a3_eveh(j,classe)
       b0 = b0_eveh(j,classe)
       b1 = b1_eveh(j,classe)
       b2 = b2_eveh(j,classe)
       b3 = b3_eveh(j,classe)
       theta0 = 0.
       theta53 = 53.
       ! Polarisation verticale
       S1_v = ((theta-theta53)/(theta0-theta53)) * ((e0-a0)/a0)
       em53_v = a3*(theta53**3) + a2*(theta53**2) + a1*theta53 + a0
       S2_v =((theta-theta0)/(theta53-theta0))*((ev(j)-em53_v)/em53_v)
       S_v = 1 + S1_v + S2_v
       emtheta_v = a3*(theta**3) + a2*(theta**2) + a1*theta + a0
       emiss_scal_v(j) = S_v * emtheta_v     
       ! Polarisation horizontale
       S1_h = ((theta-theta53)/(theta0-theta53)) * ((e0-b0)/b0)
       em53_h = b3*(theta53**3) + b2*(theta53**2) + b1*theta53 + b0
       S2_h =((theta-theta0)/(theta53-theta0))*((eh(j)-em53_h)/em53_h)
       S_h = 1 + S1_h + S2_h
       emtheta_h = b3*(theta**3) + b2*(theta**2) + b1*theta + b0
       emiss_scal_h(j) = S_h * emtheta_h     
    END DO
    ! Interpolation en frequence
    !emiss_interp_v=interp_freq(emiss_scal_v(1),emiss_scal_v(2),emiss_scal_v(3),freq)
    !emiss_interp_h=interp_freq(emiss_scal_h(1),emiss_scal_h(2),emiss_scal_h(3),freq)  
    CALL interp_freq2(emiss_scal_v(1),emiss_scal_v(2),emiss_scal_v(3),freq,emiss_interp_v)
    CALL interp_freq2(emiss_scal_h(1),emiss_scal_h(2),emiss_scal_h(3),freq,emiss_interp_h)  
    ! Cas ev<eh: on fait la moyenne entre les deux
    IF (emiss_interp_v < emiss_interp_h) THEN
       emiss_interp_v = (emiss_interp_v + emiss_interp_h)/2.
       emiss_interp_h =  emiss_interp_v
    END IF
  END SUBROUTINE emis_interp

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

  SUBROUTINE emis_interp_ind_sing_bilinear(lat,lon,theta,freq,atlas,ev,eh,stdv,stdh,verb)
    !======interpolates emissivities 
    !          IND: individual cellnumber atlas
    !         SING: singular channel
    Implicit none
    REAL, INTENT (IN)    :: lat, lon, theta, freq
    TYPE(atlas_emis_mw), INTENT (IN) :: atlas
    REAL, INTENT (OUT)   :: ev, eh   !resultat de l'interpolation
    REAL, OPTIONAL, INTENT (OUT)   :: stdv, stdh
    INTEGER, INTENT (IN) :: verb
    INTEGER :: ipos
    INTEGER :: cellnum
    REAL    :: ev_a(3),eh_a(3)      !emissivity in the atlas
    REAL    :: stdv_a(3),stdh_a(3)  !std in the atlas
    type(neibor_pix) :: npix(9)


    !Initialisations
    ev=0
    eh=0
    IF (PRESENT(stdv)) stdv=0
    IF (PRESENT(stdh)) stdh=0
    call calc_cellnum_npix(lat,lon,atlas,npix)
    ipos = 1

    IF (ipos>0) THEN

       call bilinear_interpolation( npix, atlas, atlas%emis(:,1), ev_a(1) )
       call bilinear_interpolation( npix, atlas, atlas%emis(:,2), eh_a(1) )
       call bilinear_interpolation( npix, atlas, atlas%emis(:,4), ev_a(2) )
       call bilinear_interpolation( npix, atlas, atlas%emis(:,5), eh_a(2) )
       call bilinear_interpolation( npix, atlas, atlas%emis(:,6), ev_a(3) )
       call bilinear_interpolation( npix, atlas, atlas%emis(:,7), eh_a(3) )

       CALL emis_interp(lat,lon,theta,freq,atlas%class(ipos),ev_a,eh_a,ev,eh)
       IF (PRESENT(stdv)) THEN

          call bilinear_interpolation( npix, atlas, atlas%emis_err(:,1), stdv_a(1) )
          call bilinear_interpolation( npix, atlas, atlas%emis_err(:,4), stdv_a(2) )
          call bilinear_interpolation( npix, atlas, atlas%emis_err(:,6), stdv_a(3) )

          CALL interp_freq2(stdv_a(1),stdv_a(2),stdv_a(3),freq,stdv)
       END IF
       IF (PRESENT(stdh)) THEN
          call bilinear_interpolation( npix, atlas, atlas%emis_err(:,2), stdh_a(1) )
          call bilinear_interpolation( npix, atlas, atlas%emis_err(:,5), stdh_a(2) )
          call bilinear_interpolation( npix, atlas, atlas%emis_err(:,7), stdh_a(3) )

          CALL interp_freq2(stdh_a(1),stdh_a(2),stdh_a(3),freq,stdh)
       END IF
    END IF
    IF (verb==1) THEN
       WRITE(0,*) 'Cellnum=',cellnum,' lat=',lat,' lon=',lon,' ipos=',ipos
    END IF
  END SUBROUTINE emis_interp_ind_sing_bilinear

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

 subroutine bilinear_interpolation( npix, atlas, var_in , var_out )
 implicit none
 type(neibor_pix),intent(inout) :: npix(9)
 TYPE(atlas_emis_mw), INTENT (IN) :: atlas
 real,intent(in) :: var_in(:)
 real,intent(out) :: var_out
 integer :: n
 integer :: ipos(9)

 do n = 1, 9
    ipos(n)=atlas%correspondance(npix(n)%cellnum)
 enddo

 var_out = sum( npix(1:9)%wgt *  var_in(ipos(1:9)) ) / sum( npix(1:9)%wgt )

 if ( isnan( var_out ) ) then
    print*,'var_out',var_out,  npix(1:9)%wgt 
 endif

 return 
 end subroutine bilinear_interpolation

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

 subroutine calc_cellnum_npix(lat_f,lon_f,atlas, npix)
 !======computes the cellnumber given the lat and long
 !======using the ncells included in the atlas
 IMPLICIT none
 REAL, INTENT (IN) :: lat_f, lon_f  ! lat lon point of focus
 type(atlas_emis_mw), INTENT (IN) :: atlas
 type(neibor_pix),intent(inout) :: npix(9) !neibor pixel
 INTEGER :: cellsum
 INTEGER :: i,n
 INTEGER :: ilat, ilon, ilon_m
 INTEGER :: n_loc(1)
 integer :: ilat_c,ilon_c     !center-cell lat/lon index
 real    :: lat_c, lon_c      !center-cell lat/lon [deg]

 !
 ! Atlas center-cell index from focused lat, lon point (nearest of focus point)
 !
 ilat_c=nint((lat_f+90.)/atlas%dlat)+1
 ilon_c=nint(lon_f/(360./atlas%ncells(ilat_c)))+1

 !
 ! Atlas center-cell latitude/longitude 
 !
 lat_c = real(ilat_c-1)*atlas%dlat-90.  ![deg]
 lon_c = real(ilon_c-1)*(360./real(atlas%ncells(ilat_c))) ![deg]

 !
 ! search 9 points surrounding pixcels
 !
 n = 0
 do ilat = ilat_c-1,ilat_c+1

 !
 ! middle cell locatioin (ilon_m)
 !
 ilon_m=nint(lon_c/(360./atlas%ncells(ilat)))+1  !using ilat and  center-cell longitude 

    do ilon = ilon_m-1,ilon_m+1
       n = n + 1

       !
       ! surrounding pixcel lat/lon
       !
       npix(n)%lat = real(ilat-1)*atlas%dlat-90.
       npix(n)%lon = real(ilon-1)*(360./real(atlas%ncells(ilat)))

       !
       !derive arc distance between surrdouning Atlas pixcel and focused point. 
       !
       call get_distance_latlon( lat_f,lon_f, npix(n)%lat, npix(n)%lon, npix(n)%dist )

       cellsum=0
       do i=1,ilat-1
          cellsum=atlas%ncells(i)+cellsum
      end do
      npix(n)%cellnum=cellsum+ilon
    enddo
 enddo

!
! remove 5 large pixel (just use 4 neibor pixcel) 
!

 do n = 1, 5
    n_loc(1:1) = maxloc( npix(:)%dist )
    npix(n_loc)%dist = -1.e9
 enddo

!
! assigne weight for 4 surrounding pixcels
!
 do n = 1, 9
    if( npix(n)%dist == -1.e9 ) then
      npix(n)%wgt = 0
    else
      npix(n)%wgt = 1./max( npix(n)%dist, 1.e-6 )
    endif
 enddo


  end subroutine calc_cellnum_npix

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

 subroutine get_distance_latlon( lat,lon, lat2, lon2, dist )
 implicit none
 real,intent(in) :: lat, lon, lat2, lon2 ![deg]
 real,intent(out) :: dist ![km] 
 real :: v(3), v2(3)
 real,parameter :: pi=3.14159265358979323846e0  ! pai 
 real,parameter :: degrad=pi/180.d0             ! conversion from degree to radiance

 call latlon2cartesian_s(lat*degrad,lon*degrad,v)
 call latlon2cartesian_s(lat2*degrad,lon2*degrad,v2)
 call get_arc_distance_s( v, v2, dist)

 return
 end subroutine get_distance_latlon

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

 subroutine latlon2cartesian_s(lat,lon,v)
 implicit none
!-----------------------------------------------------------------------------------
! Comments: Using latitude and longitude, find out vector point on reference ellipsoid
!           in ECEF cartesian coordinate. 
!
!           Note that (x,y,z) = (rs, 0, 0) => lat= 0deg, lon= 0deg
!                     (x,y,z) = ( 0,rs, 0) => lat= 0deg, lon=90deg
!                     (x,y,z) = ( 0, 0,rs) => lat=90deg, lon= 0deg
!
! History:
! 12/2010  Toshi Matsui@NASA GSFC : Initial. 
!
! References:
!-----------------------------------------------------------------------------------
 real,intent(in) :: lat, lon !earth's latitude and longitude [rad] 
 real,intent(out) :: v(3) !cartesian cordinate position [km]
 real :: re !distance from earth's core to ellipsoidal surface  [km]

!
! Get distance from the earth's core to the ellipoidal surface. 
!
 call earth_radius_s(lat, re)

!
! Derive vector x,y,z
!
 v(1) = re*cos(lat)*cos(lon)  ! x [km]

 v(2) = re*cos(lat)*sin(lon)  ! y [km]

 v(3) = re*sin(lat)           ! z [km]

 return
 end subroutine latlon2cartesian_s

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

 subroutine earth_radius_s( lat, re )
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Compute earth's radius as spherical oblate (for a given latitude). 
!           Namely, from the earth's core of earth to reference ellipsoid based ono WGS 1984. 
! History:
! 10/2010  Toshi Matsui@NASA GSFC : Initial
!
! References:
!-----------------------------------------------------------------------------------------------------
 real,intent(in) :: lat ! latitude [rad]
 real,intent(out) :: re ! distance between earth's core and ellipsoid [km]

 real :: a,b,c,d
 real,parameter :: req=6378.1370d0     ! earth's equatorial radius [km]
 real,parameter :: rpl=6356.7523142d0  ! earth's polar radius [km]

 a = req*req*cos(lat)*req*req*cos(lat)
 b = rpl*rpl*sin(lat)*rpl*rpl*sin(lat)
 c = req*cos(lat)*req*cos(lat)
 d = rpl*sin(lat)*rpl*sin(lat)

 re = sqrt( (a+b)/(c+d) )   ! distance between earth's core and ellipsoid [km]

 return
 end subroutine earth_radius_s

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

 subroutine get_arc_distance_s( v1, v2, arc_dist)
 implicit none
!---------------------------------------------------------------------
! Comments:  Estimate arc distance between two vector points on earth's surface. 
!
! History:
! 11/2011  Toshi Matsui@NASA GSFC         : initial
!
! References:
!----------------------------------------------------------------------
 real,intent(in)  :: v1(3),v2(3)
 real,intent(out) :: arc_dist
 real :: a,b,c  ![km]
 real :: alpha ![rad]
 real,parameter :: v_ec(3) = (/0.d0,0.d0,0.d0/) ! ECEF position of earth's core. [km]

 call get_distance_s( v1 ,   v2, a)
 call get_distance_s( v1 , v_ec, b)
 call get_distance_s( v2 , v_ec, c)

 alpha = acos( ( b*b + c*c - a*a ) / (2.*b*c) )  ![rad]

 arc_dist = alpha * 0.5*(b+c)

 return
 end subroutine get_arc_distance_s

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

 subroutine get_distance_s( v1, v2, dist)
 implicit none
 real,intent(in) :: v1(3)  !3D vector 1
 real,intent(in) :: v2(3)  !3D vector 2
 real,intent(out) :: dist

!
! distance [unit is same as v1 and v2 
!
 dist = sqrt( (v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 + (v1(3)-v2(3))**2 )

 return
 end subroutine get_distance_s

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

!END MODULE mod_mwatlas
end  module module_telsem

!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM 
!TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM TELSEM

