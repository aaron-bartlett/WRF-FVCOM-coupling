!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

module module_openssp_cmb
implicit none


private   !all subroutines are encapsulated.

public :: readTablesLiang2, getsnowp2_sdsu


save

!module microwFreq
   real, allocatable :: mfreq(:) ! stores the passive microwave frequencies
   integer :: nmfreqR            ! # of passive microwave frequencies
                                 ! these variables are set in readTables but have to be consistent
                                 ! with their definition in mie2
!end module microwFreq

!module Tables2                              ! look up table module
  ! this module defines the structure of the look up tables
  ! to see how these tables are used have a look at subroutines 
  ! integratecvHB and integratestHB in retTablesInt.f90 
  integer :: nbins, nbinG, nbinS2, nbinH, &                       ! # of entries (rows) in the lookup tables
       nmu                                  ! # of mus in the look up tables
  real :: zmin,&                            ! min reflectivity in the tables (dBZ)
       dzbin                                ! reflectivity bin size (dB) 
                                            ! this is used to quickly locate an entry in the 
                                            ! lookup table as i=(z-zmin)/dzbin

  integer, allocatable :: i0max(:)          ! the last rows in the look up tables contain data only
  integer, allocatable :: i0BBmax(:) 
                                            ! for rain and melting snow
                                            ! so (zsnow-zmin)/dzbin cannot be greater than i0max
  integer :: nmfreq
  integer, allocatable :: ibranch(:)
  real :: d013Smax, d013BBmax, d013max
  real :: d013Smin, d013BBmin, d013min
 
  integer :: i0HTS(400), i0HT(400), i0HTBB(400)

  real,allocatable :: z13Table(:,:), att13Table(:,:), &
       pwc13Table(:,:)  ! rain reflectivity (dBZ), attenuation (dB/km), and
                        ! precipitation water conter(g/m3) at Ku-band 
                       

  real,allocatable :: z13TableS(:,:), att13TableS(:,:), &
       pwc13TableS(:,:) ! snow reflectivity (dBZ), attenuation (dB/km), and
                        ! precipitation water conter(g/m3) at Ku-band 

  real,allocatable :: z13TableBB(:,:), att13TableBB(:,:), &
       pwc13TableBB(:,:)! BB reflectivity (dBZ), attenuation (dB/km), and
                        ! precipitation water conter(g/m3) at Ku-band 

  ! the next three statements have same meaning as the previous but for Ka-band
  
  real,allocatable :: z35Table(:,:), att35Table(:,:), &   
       pwc35Table(:,:)
  
  real,allocatable :: z35TableS(:,:), att35TableS(:,:), &
      pwc35TableS(:,:)

  real,allocatable :: z35TableBB(:,:), att35TableBB(:,:), &
       pwc35TableBB(:,:)

  real,allocatable :: kextTableBB(:,:,:), salbTableBB(:,:,:), &
       asymTableBB(:,:,:)  ! extinction, scattering albedo and asymmetry fact tables for the BB

  real,allocatable :: kextTableS(:,:,:), salbTableS(:,:,:), &
       asymTableS(:,:,:)   ! extinction, scattering albedo and asymmetry fact tables for snow
!------------------------------graupel tables---------------------!
  real,allocatable :: kextTableG(:,:,:), salbTableG(:,:,:), &
       asymTableG(:,:,:)   ! extinction, scattering albedo and asymmetry fact tables for snow

  real,allocatable :: pr13TableG(:,:), d013TableG(:,:)
  real,allocatable :: pr35TableG(:,:), d035TableG(:,:)
  real,allocatable :: z13TableG(:,:), att13TableG(:,:), &
       pwc13TableG(:,:)

  real,allocatable :: z35TableG(:,:), att35TableG(:,:), &
       pwc35TableG(:,:)

!------------------------------snow exp tables---------------------!
  real,allocatable :: kextTableS2(:,:,:), salbTableS2(:,:,:), &
       asymTableS2(:,:,:)   ! extinction, scattering albedo and asymmetry fact tables for snow

  real,allocatable :: pr13TableS2(:,:), d013TableS2(:,:)
  real,allocatable :: pr35TableS2(:,:), d035TableS2(:,:)
  real,allocatable :: z13TableS2(:,:), att13TableS2(:,:), &
       pwc13TableS2(:,:)

  real,allocatable :: z35TableS2(:,:), att35TableS2(:,:), &
       pwc35TableS2(:,:)

!-----------------------------hail-------------------------------------!

  real,allocatable :: kextTableH(:,:,:), salbTableH(:,:,:), &
       asymTableH(:,:,:)   ! extinction, scattering albedo and asymmetry fact tables for snow

  real,allocatable :: pr13TableH(:,:), d013TableH(:,:)
  real,allocatable :: pr35TableH(:,:), d035TableH(:,:)
  real,allocatable :: z13TableH(:,:), att13TableH(:,:), &
       pwc13TableH(:,:)

  real,allocatable :: z35TableH(:,:), att35TableH(:,:), &
       pwc35TableH(:,:)
!-------------------------------------------------------------------------!

  real,allocatable :: kextTable(:,:,:), salbTable(:,:,:), &
       asymTable(:,:,:)    ! extinction, scattering albedo and asymmetry fact tables for rain

  real, allocatable :: pwc35min(:), dpwc35(:),pwc35minS(:), &
       dpwc35S(:),pwc35minBB(:), dpwc35BB(:)  ! these variables end up not being used at all

  real, allocatable :: pr13Table(:,:),pr13TableBB(:,:),pr13TableS(:,:), & !precipitation rates
       d013Table(:,:),d013TableBB(:,:),d013TableS(:,:), &                 !d0
       logN13Table(:,:),logN13TableBB(:,:),logN13TableS(:,:)              ! log(Nw*f(mu)/D0^mu)

!begin  WSO 8/8/13
  real, allocatable :: mu_tab(:)
!end    WSO 8/8/13

!end module Tables2

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

!!SJM 12/19/2013
!!DSD Working Group Tables - rain only
subroutine readTablesDSDWG(extnmu,extnmfreq)
  !! SJM This routine should be used after readTablesLiang2, since the structures are assumed to be already allocated
  ! this reads in the look up tables and stores them in the structures define in
  ! module Tables2
  ! to see how these tables are used have a look at subroutines
  ! integratecvHB and integratestHB in retTablesInt.f90
  !use Tables2
  !use microwFreq
  implicit none
  integer :: i, j, imu, extnmu,extnmfreq
  real :: mu
  real :: pwc13, pwc13s, pwc13bb
  character*4 mus
  nmu=extnmu
  nmfreq=extnmfreq
  nmfreqR=extnmfreq
  nbins=289                     ! the number of rows in the lookup tables
  zmin=-12                        ! minimum reflectivity in the lookup tables
  dzbin=0.25                    ! reflectivity bin size in the lookup tables


  kextTable=0
  salbTable=0
  asymTable=0

!  print*, nmfreq, extnmu, nmu

  do imu=1,extnmu

     write(mus,'(A2,I2.2)') '0p',imu*6+11
     open(10,file='TablesN/DmSm_aSy'//mus//'_sphere_Dmax3p0Dm_rain.dat')
     print*, 'TablesN/DmSm_aSy'//mus//'_sphere_Dmax3p0Dm_rain.dat'

     do i=1,nbins       
        read(10,*) &
             pwc13table(i,imu), z13table(i,imu), att13Table(i,imu), &
             z35Table(i,imu), att35Table(i,imu), pr13Table(i,imu), d013Table(i,imu), &
             kextTable(i,1,imu), salbTable(i,1,imu), asymTable(i,1,imu), &
             kextTable(i,2,imu), salbTable(i,2,imu), asymTable(i,2,imu), &
             kextTable(i,3,imu), salbTable(i,3,imu), asymTable(i,3,imu), &
             kextTable(i,4,imu), salbTable(i,4,imu), asymTable(i,4,imu), &
             kextTable(i,5,imu), salbTable(i,5,imu), asymTable(i,5,imu)
        !write(*,*) z13Table(i,imu), pr13Table(i,imu), d013Table(i,imu)
     enddo

     close(10)

     !To do: Add in 166 and 183 GHz tables. Need to compute values from Mie code, since they aren't in Christopher's tables.

  enddo

end subroutine readTablesDSDWG

!!Johnson/Kuo ice tables
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


!!end SJM 12/19/2013

subroutine readTables2(extnmu,extnmfreq)
  ! this reads in the look up tables and stores them in the structures define in
  ! module Tables2
  ! to see how these tables are used have a look at subroutines 
  ! integratecvHB and integratestHB in retTablesInt.f90 
  !use Tables2
  !use microwFreq
  implicit none
  integer :: i, j, imu, extnmu,extnmfreq
  real :: mu
  real :: pwc13, pwc13s, pwc13bb
  character*5 mus
  nmu=extnmu
  nmfreq=extnmfreq
  nmfreqR=extnmfreq
  nbins=240                     ! the number of rows in the lookup tables
  zmin=-5                        ! minimum reflectivity in the lookup tables
  dzbin=0.25                    ! reflectivity bin size in the lookup tables
  allocate(i0max(nmu))
  allocate(ibranch(nmu))
  allocate(mfreq(nmfreq))

  allocate( z13Table(nbins,nmu), att13Table(nbins,nmu), &
       pwc13Table(nbins,nmu))
  allocate( z13TableS(nbins,nmu), att13TableS(nbins,nmu), &
       pwc13TableS(nbins,nmu))
  allocate( z13TableBB(nbins,nmu), att13TableBB(nbins,nmu), &
       pwc13TableBB(nbins,nmu))
  
  allocate( z13TableG(nbins,nmu), att13TableG(nbins,nmu), &
       pwc13TableG(nbins,nmu))

  allocate( z35TableG(nbins,nmu), att35TableG(nbins,nmu), &
       pwc35TableG(nbins,nmu))

  allocate( z13TableH(nbins,nmu), att13TableH(nbins,nmu), &
       pwc13TableH(nbins,nmu))
  
  allocate( z35TableH(nbins,nmu), att35TableH(nbins,nmu), &
       pwc35TableH(nbins,nmu))
  
  allocate( z35Table(nbins,nmu), att35Table(nbins,nmu), &
       pwc35Table(nbins,nmu))
  allocate( z35TableS(nbins,nmu), att35TableS(nbins,nmu), &
       pwc35TableS(nbins,nmu))
  allocate( z35TableBB(nbins,nmu), att35TableBB(nbins,nmu), &
       pwc35TableBB(nbins,nmu))

  allocate( kextTable(nbins,nmfreq,nmu), salbTable(nbins,nmfreq,nmu), &
       asymTable(nbins,nmfreq,nmu))
  
  allocate( kextTableG(nbins,nmfreq,nmu), salbTableG(nbins,nmfreq,nmu), &
       asymTableG(nbins,nmfreq,nmu))

  allocate( kextTableH(nbins,nmfreq,nmu), salbTableH(nbins,nmfreq,nmu), &
       asymTableH(nbins,nmfreq,nmu))

  allocate( kextTableBB(nbins,nmfreq,nmu), salbTableBB(nbins,nmfreq,nmu), &
       asymTableBB(nbins,nmfreq,nmu))
  allocate( kextTableS(nbins,nmfreq,nmu), salbTableS(nbins,nmfreq,nmu), &
       asymTableS(nbins,nmfreq,nmu))

  allocate(pr13Table(nbins,nmu),pr13TableBB(nbins,nmu),pr13TableS(nbins,nmu), &
       d013Table(nbins,nmu),d013TableBB(nbins,nmu),d013TableS(nbins,nmu), &
       logN13Table(nbins,nmu),logN13TableBB(nbins,nmu),&
       logN13TableS(nbins,nmu),pr13TableG(nbins,nmu),d013TableG(nbins,nmu),&
       pr35TableG(nbins,nmu),d035TableG(nbins,nmu),&
       d013TableH(nbins,nmu),&
       pr35TableH(nbins,nmu),d035TableH(nbins,nmu))

  allocate(pwc35min(nmu), dpwc35(nmu),pwc35minS(nmu), &
       dpwc35S(nmu),pwc35minBB(nmu), dpwc35BB(nmu))

!begin  WSO 8/8/13
  allocate (mu_tab(nmu))
!end    WSO 8/8/13

  
  20 format((F4.1))
  do imu=1,extnmu
     mu=imu-3.
!begin  WSO 8/8/13
     mu_tab(imu) = mu
     write(*, '("imu: ", i5, "  mu: ", f10.4, "  mu_tab: ", f10.4)') imu, mu, mu_tab(imu)
!end    WSO 8/8/13
     write(mus,20) mu
     if(mus(1:1).ne.'-') mus(1:1)='0'
     
     open(10,file='TablesN/tables.13-35GHz.mu.'//mus)
     !print*, 'TablesN/tables.13-35GHz.mu.'//mus, imu, extnmu
     read(10,*) dpwc35(imu), pwc35min(imu)
     read(10,*) dpwc35BB(imu), pwc35minBB(imu)
     read(10,*) dpwc35S(imu), pwc35minS(imu)
     ibranch(imu)=1
     do i=1,nbins
        read(10,*) pwc13Table(i,imu), pwc13TableBB(i,imu),&
             pwc13TableS(i,imu), &
             z13Table(i,imu), z13TableBB(i,imu), z13TableS(i,imu), &
             att13Table(i,imu), att13TableBB(i,imu), att13TableS(i,imu), &
             z35Table(i,imu), z35TableBB(i,imu), z35TableS(i,imu), &
             att35Table(i,imu), att35TableBB(i,imu), att35TableS(i,imu), &
             pr13Table(i,imu),pr13TableBB(i,imu),pr13TableS(i,imu), &
             d013Table(i,imu),d013TableBB(i,imu),d013TableS(i,imu), &
             logN13Table(i,imu),logN13TableBB(i,imu),logN13TableS(i,imu)
        if( z13Table(i,imu)- z35Table(i,imu)<z13Table(ibranch(imu),imu)- z35Table(ibranch(imu),imu)) then
           ibranch(imu)=i
        endif
        
     enddo
     !print*, z13TableS(1:5,imu)
     !print*, z35TableS(1:5,imu)
     read(10,*) i0max(imu)
 
     close(10)
     mfreq(7)=186.
     mfreq(8)=190.
     open(10,file='TablesN/tables.microw.mu.'//mus)
     read(10,*) nmfreq
     read(10,*) (mfreq(i),i=1,nmfreq-2)       ! reads the microwave frequencies
     do i=1,nbins
        read(10,*) pwc13, pwc13bb,pwc13s, &
             (kextTable(i,j,imu),j=1,nmfreq),(salbTable(i,j,imu),j=1,nmfreq),&
             (asymTable(i,j,imu),j=1,nmfreq),&
             (kextTableBB(i,j,imu),j=1,nmfreq),(salbTableBB(i,j,imu),j=1,nmfreq),&
             (asymTableBB(i,j,imu),j=1,nmfreq),&
             (kextTableS(i,j,imu),j=1,nmfreq),(salbTableS(i,j,imu),j=1,nmfreq),&
             (asymTableS(i,j,imu),j=1,nmfreq)
     enddo
  enddo
  close(10)
  do imu=1,-extnmu
     open(10,file='RadarOnly/TablesN/scatTables.Radar.LLiang')
     do i=1,-nbins
        read(10,*) pwc13Table(i,imu), &
             z13Table(i,imu),  &
             att13Table(i,imu), &
             z35Table(i,imu),  &
             att35Table(i,imu),  &
             pr13Table(i,imu), &
             d013Table(i,imu)
     enddo
     close(10)
  enddo
  !stop
  i0max=i0max-1
  !write(*,*) i0max
  ! stop
10 format(30(F10.6))

end subroutine readTables2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine readTablesLiang(extnmu,extnmfreq)
  ! this reads in the look up tables and stores them in the structures define in
  ! module Tables2
  ! to see how these tables are used have a look at subroutines 
  ! integratecvHB and integratestHB in retTablesInt.f90 
  !use Tables2
  !use microwFreq
  implicit none
  integer :: i, j, imu, extnmu,extnmfreq
  real :: mu
  real :: pwc13, pwc13s, pwc13bb
  character*5 mus
  nmu=extnmu
  nmfreq=extnmfreq
  nmfreqR=extnmfreq
  nbins=700                     ! the number of rows in the lookup tables
  zmin=-10                        ! minimum reflectivity in the lookup tables
  dzbin=0.25                    ! reflectivity bin size in the lookup tables
  allocate(i0max(nmu))
  allocate(i0BBmax(nmu))
  allocate(ibranch(nmu))
  allocate(mfreq(nmfreq))

  allocate( z13Table(nbins,nmu), att13Table(nbins,nmu), &
       pwc13Table(nbins,nmu))
  allocate( z13TableS(nbins,nmu), att13TableS(nbins,nmu), &
       pwc13TableS(nbins,nmu))
  allocate( z13TableBB(nbins,nmu), att13TableBB(nbins,nmu), &
       pwc13TableBB(nbins,nmu))

   allocate( z35Table(nbins,nmu), att35Table(nbins,nmu), &
       pwc35Table(nbins,nmu))
  allocate( z35TableS(nbins,nmu), att35TableS(nbins,nmu), &
       pwc35TableS(nbins,nmu))
  allocate( z35TableBB(nbins,nmu), att35TableBB(nbins,nmu), &
       pwc35TableBB(nbins,nmu))

  allocate( kextTable(nbins,nmfreq,nmu), salbTable(nbins,nmfreq,nmu), &
       asymTable(nbins,nmfreq,nmu))
  allocate( kextTableBB(nbins,nmfreq,nmu), salbTableBB(nbins,nmfreq,nmu), &
       asymTableBB(nbins,nmfreq,nmu))
  allocate( kextTableS(nbins,nmfreq,nmu), salbTableS(nbins,nmfreq,nmu), &
       asymTableS(nbins,nmfreq,nmu))

  allocate(pr13Table(nbins,nmu),pr13TableBB(nbins,nmu),pr13TableS(nbins,nmu), &
       d013Table(nbins,nmu),d013TableBB(nbins,nmu),d013TableS(nbins,nmu), &
       logN13Table(nbins,nmu),logN13TableBB(nbins,nmu),logN13TableS(nbins,nmu))

  allocate(pwc35min(nmu), dpwc35(nmu),pwc35minS(nmu), &
       dpwc35S(nmu),pwc35minBB(nmu), dpwc35BB(nmu))

!begin  WSO 8/8/13
  allocate (mu_tab(nmu))
!end    WSO 8/8/13
  
  20 format((F4.1))


  do imu=1,extnmu
     mu=imu-3.
!begin  WSO 8/8/13
     mu_tab(imu) = mu
     write(*, '("imu: ", i5, "  mu: ", f10.4, "  mu_tab: ", f10.4)') imu, mu, mu_tab(imu)
!end    WSO 8/8/13
     write(mus,20) mu
     if(mus(1:1).ne.'-') mus(1:1)='0'
    
     open(10,file='nw_d0_dsd/liang.table.rho01mu2.mu4')
   
     ibranch(imu)=1
     do i=1,nbins
        read(10,*) pwc13Table(i,imu), pwc13TableBB(i,imu),&
             pwc13TableS(i,imu), &
             z13Table(i,imu), z13TableBB(i,imu), z13TableS(i,imu), &
             att13Table(i,imu), att13TableBB(i,imu), att13TableS(i,imu), &
             z35Table(i,imu), z35TableBB(i,imu), z35TableS(i,imu), &
             att35Table(i,imu), att35TableBB(i,imu), att35TableS(i,imu), &
             pr13Table(i,imu),pr13TableBB(i,imu),pr13TableS(i,imu), &
             d013Table(i,imu),d013TableBB(i,imu),d013TableS(i,imu)
     
        if( z13Table(i,imu)- z35Table(i,imu)<z13Table(ibranch(imu),imu)- z35Table(ibranch(imu),imu)) then
           ibranch(imu)=i
        endif
     enddo
    
     close(10)
    

  enddo
!  write(*,*) ibranch
!  stop
!  write(*,*) att35Table(:,1)
!  write(*,*) att35TableS(:,1)
!  write(*,*) att35TableBB(:,1)
!  stop
  close(10)
  
  !stop
  
  !write(*,*) i0max
  ! stop
10 format(30(F10.6))

end subroutine readTablesLiang

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

subroutine readTablesLiang2(masterproc,dir, extnmu,extnmfreq)
  ! this reads in the look up tables and stores them in the structures define in
  ! module Tables2
  ! to see how these tables are used have a look at subroutines 
  ! integratecvHB and integratestHB in retTablesInt.f90 
  !use Tables2
  !use microwFreq
  implicit none
  logical,intent(in) :: masterproc
  character(len=*),intent(in) :: dir  !directory
  integer :: i, j, imu, extnmu,extnmfreq, k
  real :: mu
  real :: pwc13, pwc13s, pwc13bb, dum1, dum2, dum3
  character*5 mus
  integer :: idum
  real :: dumfreq
  nmu=extnmu
  nmfreq=extnmfreq
  nmfreqR=extnmfreq
  nbins=289                     ! the number of rows in the lookup tables
  zmin=-12                        ! minimum reflectivity in the lookup tables
  dzbin=0.25                    ! reflectivity bin size in the lookup tables
  if(allocated(i0max)) return

  if(masterproc) print*,'MSG readTablesLiang2: reading CMB OpenSSP table'

  allocate(i0max(nmu))
  allocate(i0BBmax(nmu))
  allocate(ibranch(nmu))
  allocate(mfreq(nmfreq))
  print*, nmfreq
  !stop
  allocate( z13Table(nbins,nmu), att13Table(nbins,nmu), &
       pwc13Table(nbins,nmu))
  allocate( z13TableS(nbins,nmu), att13TableS(nbins,nmu), &
       pwc13TableS(nbins,nmu))
  allocate( z13TableS2(nbins,nmu), att13TableS2(nbins,nmu), &
       pwc13TableS2(nbins,nmu))

  allocate( z13TableBB(nbins,nmu), att13TableBB(nbins,nmu), &
       pwc13TableBB(nbins,nmu))
  
  allocate( z13TableG(nbins,nmu), att13TableG(nbins,nmu), &
       pwc13TableG(nbins,nmu))
  allocate( z35TableG(nbins,nmu), att35TableG(nbins,nmu), &
       pwc35TableG(nbins,nmu))

  allocate( z13TableH(nbins,nmu), att13TableH(nbins,nmu), &
       pwc13TableH(nbins,nmu))
  allocate( z35TableH(nbins,nmu), att35TableH(nbins,nmu), &
       pwc35TableH(nbins,nmu))

  allocate( z35Table(nbins,nmu), att35Table(nbins,nmu), &
       pwc35Table(nbins,nmu))
  allocate( z35TableS(nbins,nmu), att35TableS(nbins,nmu), &
       pwc35TableS(nbins,nmu))
  allocate( z35TableS2(nbins,nmu), att35TableS2(nbins,nmu), &
       pwc35TableS2(nbins,nmu))
  allocate( z35TableBB(nbins,nmu), att35TableBB(nbins,nmu), &
       pwc35TableBB(nbins,nmu))

  allocate( kextTableG(nbins,nmfreq,nmu), salbTableG(nbins,nmfreq,nmu), &
       asymTableG(nbins,nmfreq,nmu))
  allocate( kextTableH(nbins,nmfreq,nmu), salbTableH(nbins,nmfreq,nmu), &
       asymTableH(nbins,nmfreq,nmu))

  allocate( kextTable(nbins,nmfreq,nmu), salbTable(nbins,nmfreq,nmu), &
       asymTable(nbins,nmfreq,nmu))
  allocate( kextTableBB(nbins,nmfreq,nmu), salbTableBB(nbins,nmfreq,nmu), &
       asymTableBB(nbins,nmfreq,nmu))
  allocate( kextTableS(nbins,nmfreq,nmu), salbTableS(nbins,nmfreq,nmu), &
       asymTableS(nbins,nmfreq,nmu))
  allocate( kextTableS2(nbins,nmfreq,nmu), salbTableS2(nbins,nmfreq,nmu), &
       asymTableS2(nbins,nmfreq,nmu))

  allocate(pr13Table(nbins,nmu),pr13TableBB(nbins,nmu),pr13TableS(nbins,nmu), &
       d013Table(nbins,nmu),d013TableBB(nbins,nmu),d013TableS(nbins,nmu), &
       logN13Table(nbins,nmu),logN13TableBB(nbins,nmu),&
       logN13TableS(nbins,nmu),&
       pr13TableG(nbins,nmu),d013TableG(nbins,nmu),&
       pr35TableG(nbins,nmu),d035TableG(nbins,nmu),&
       pr13TableH(nbins,nmu),d013TableH(nbins,nmu),&
       pr35TableH(nbins,nmu),d035TableH(nbins,nmu),&
       pr13TableS2(nbins,nmu),d013TableS2(nbins,nmu),&
       pr35TableS2(nbins,nmu),d035TableS2(nbins,nmu))

  allocate(pwc35min(nmu), dpwc35(nmu),pwc35minS(nmu), &
       dpwc35S(nmu),pwc35minBB(nmu), dpwc35BB(nmu))

!begin  WSO 8/8/13
  allocate (mu_tab(nmu))
!end    WSO 8/8/13

  20 format((F4.1))

  i0max=nbins
  i0BBmax=nbins

  kextTableS=0
  kextTable=0
  kextTableBB=0
  salbTableS=0
  asymTableS=0
  salbTable=0
  asymTable=0
  salbTableBB=0
  asymTableBB=0

!  print*, nmfreq, extnmu, nmu

  mfreq(7)=186.
  mfreq(8)=190.
  do imu=1,extnmu
     mu=imu-3.
!begin  WSO 8/8/13
     mu_tab(imu) = mu
 if(masterproc)     write(*, '("imu: ", i5, "  mu: ", f10.4, "  mu_tab: ", f10.4)') imu, mu, mu_tab(imu)
!end    WSO 8/8/13
 if(masterproc)     write(mus,20) mu
     
     if(mus(1:1).ne.'-') then 
        mus(1:1)='0'
        open(10,file=trim(dir)//'tables.microw.mu.'//mus)
        read(10,*) nmfreq
        read(10,*) (mfreq(i),i=1,nmfreq)  
        nmfreq=extnmfreq
        mfreq(6)=165.5
        close(10)
     endif
     !print*, imu
!begin  WSO 12/30/13 insert Dm-indexed tables
!     open(10,file='nw_d0_dsd/snowgTables.rho01.mu2.2')

     open(10,file=trim(dir)//'snowTableRho0.3')
     read(10,*) nbinG
     read(10,*) dumfreq
     do i=1,nbinG
        read(10,*) idum, z13TableG(i,imu), att13TableG(i,imu), &
             d013TableG(i,imu), pwc13TableG(i,imu), pr13TableG(i,imu) 
        pwc13TableG(i,imu)=log10(pwc13TableG(i,imu))
     enddo
     read(10,*) dumfreq
     do i=1,nbinG
        read(10,*) idum, z35TableG(i,imu), att35TableG(i,imu), &
             d035TableG(i,imu), pwc35TableG(i,imu), pr35TableG(i,imu) 
     enddo
     do j=1,7
        read(10,*) dumfreq
        !print*, dumfreq
        do i=1,nbinG
           read(10,*) idum, kextTableG(i,j,imu), salbTableG(i,j,imu), &
                asymTableG(i,j,imu), dum1, dum2
           if(kextTableG(i,j,imu)>0) then
              salbTableG(i,j,imu)=salbTableG(i,j,imu)/kextTableG(i,j,imu)
           endif
        enddo
     enddo
     salbTableG(:,8,imu)=salbTableG(:,7,imu)
     kextTableG(:,8,imu)=kextTableG(:,7,imu)
     asymTableG(:,8,imu)=asymTableG(:,7,imu)
     close(10)

     
     open(10,file=trim(dir)//'snowTableRho0.8')
     read(10,*) nbinH
     read(10,*) dumfreq
     do i=1,nbinH
        read(10,*) idum, z13TableH(i,imu), att13TableH(i,imu), &
             d013TableH(i,imu), pwc13TableH(i,imu), pr13TableH(i,imu) 
        pwc13TableH(i,imu)=log10(pwc13TableH(i,imu))
     enddo
     read(10,*) dumfreq
     do i=1,nbinH
        read(10,*) idum, z35TableH(i,imu), att35TableH(i,imu), &
             d035TableH(i,imu), pwc35TableH(i,imu), pr35TableH(i,imu) 
     enddo
     do j=1,7
        read(10,*) dumfreq
        !print*, dumfreq
        do i=1,nbinH
           read(10,*) idum, kextTableH(i,j,imu), salbTableH(i,j,imu), &
                asymTableH(i,j,imu), dum1, dum2
           if(kextTableH(i,j,imu)>0) then
              salbTableH(i,j,imu)=salbTableH(i,j,imu)/kextTableH(i,j,imu)
           endif
        enddo
     enddo
     salbTableH(:,8,imu)=salbTableH(:,7,imu)
     kextTableH(:,8,imu)=kextTableH(:,7,imu)
     asymTableH(:,8,imu)=asymTableH(:,7,imu)
     close(10)

!begin WSO 12/12/16 replace spherical snow tables with nonspherical snow tables
!     open(10,file='nw_dm_dsd/snowTableRho0.1')
     open(10,file=trim(dir)//'table.dda.266.mu2.dm.Dec1.2016')
!end   WSO 12/12/16
     read(10,*) nbinS2
     read(10,*) dumfreq
     do i=1,nbinS2
        read(10,*) idum, z13TableS2(i,imu), att13TableS2(i,imu), &
             d013TableS2(i,imu), pwc13TableS2(i,imu), pr13TableS2(i,imu) 
        pwc13TableS2(i,imu)=log10(pwc13TableS2(i,imu))
     enddo
     read(10,*) dumfreq
     do i=1,nbinS2
        read(10,*) idum, z35TableS2(i,imu), att35TableS2(i,imu), &
             d035TableS2(i,imu), pwc35TableS2(i,imu), pr35TableS2(i,imu) 
     enddo
     do j=1,7
        read(10,*) dumfreq
        print*, dumfreq
        do i=1,nbinS2
           read(10,*) idum, kextTableS2(i,j,imu), salbTableS2(i,j,imu), &
                asymTableS2(i,j,imu), dum1, dum2
           if(kextTableS2(i,j,imu)>0) then
              salbTableS2(i,j,imu)=salbTableS2(i,j,imu)/kextTableS2(i,j,imu)
           endif
        enddo
     enddo
     salbTableS2(:,8,imu)=salbTableS2(:,7,imu)
     kextTableS2(:,8,imu)=kextTableS2(:,7,imu)
     asymTableS2(:,8,imu)=asymTableS2(:,7,imu)
     close(10)


     open(10,file=trim(dir)//'snowgTables.Dm.rho01.mu2.2')!end    WSO 12/30/13      
     
     do i=1,nbins
        read(10,*) &
             z13TableS(i,imu), att13TableS(i,imu),  &
             pwc13TableS(i,imu), &
             z35TableS(i,imu), att35TableS(i,imu),  &
             kextTableS(i,4,imu), salbTableS(i,4,imu), &
             asymTableS(i,4,imu),  &      
             kextTableS(i,5,imu), salbTableS(i,5,imu), &
             asymTableS(i,5,imu),  &      
             pr13TableS(i,imu), d013TableS(i,imu), &
             kextTableS(i,2,imu), salbTableS(i,2,imu), &
             asymTableS(i,2,imu),  &      
             kextTableS(i,3,imu), salbTableS(i,3,imu), &
             asymTableS(i,3,imu), &
             kextTableS(i,6,imu), salbTableS(i,6,imu), &  
             asymTableS(i,6,imu), & 
             kextTableS(i,7,imu), salbTableS(i,7,imu), &  
             asymTableS(i,7,imu), &
             kextTableS(i,1,imu), salbTableS(i,1,imu), &  
             asymTableS(i,1,imu)
        !print*, kextTableS(i,5:7,imu)
        salbTableS(i,1:7,imu)=salbTableS(i,1:7,imu)/kextTableS(i,1:7,imu)
        salbTables(i,8,imu)=salbTables(i,7,imu)
        kextTables(i,8,imu)=kextTables(i,7,imu)
        asymTables(i,8,imu)=asymTables(i,7,imu)

        if(i>2) then
           if( att13TableS(i-1,imu)>0 .and. att13TableS(i,imu)<-90) then
              i0max(imu)=i
           endif
        endif
     enddo
     !stop
     close(10)
     
!begin  WSO 12/30/13 insert Dm-indexed tables
     open(10,file=trim(dir)//'raingTables.Dm.rho01.mu2.2') !Sept 17, 2015 MG
     !open(10,file='nw_dm_dsd/raingTables.Dm.rho01.u.mu3')
!end    WSO 12/30/13
     
     
     do i=1,nbins
        read(10,*) &
             z13Table(i,imu), att13Table(i,imu),  &
             pwc13Table(i,imu), &
             z35Table(i,imu), att35Table(i,imu),  &
             kextTable(i,4,imu), salbTable(i,4,imu), &
             asymTable(i,4,imu),  &      
             kextTable(i,5,imu), salbTable(i,5,imu), &
             asymTable(i,5,imu),  &      
             pr13Table(i,imu), d013Table(i,imu), &
             kextTable(i,2,imu), salbTable(i,2,imu), &
             asymTable(i,2,imu),  &      
             kextTable(i,3,imu), salbTable(i,3,imu), &
             asymTable(i,3,imu), &
             kextTable(i,6,imu), salbTable(i,6,imu), &
             asymTable(i,6,imu), &
             kextTable(i,7,imu), salbTable(i,7,imu), &
             asymTable(i,7,imu), &
             kextTable(i,1,imu), salbTable(i,1,imu), &  
             asymTable(i,1,imu)

        salbTable(i,1:7,imu)=salbTable(i,1:7,imu)/kextTable(i,1:7,imu)
        salbTable(i,8,imu)=salbTable(i,7,imu)
        kextTable(i,8,imu)=kextTable(i,7,imu)
        asymTable(i,8,imu)=asymTable(i,7,imu)
        !write(*,*) z13Table(i,imu), pr13Table(i,imu), d013Table(i,imu) 
     enddo
     
     close(10)
     
!begin  WSO 12/30/13 insert Dm-indexed tables
!     open(10,file='nw_d0_dsd/bbgTables.rho01.mu2.2')
     open(10,file=trim(dir)//'bbgTables.Dm.rho02.mu2.2')
!end    WSO 12/30/13

     
     
     do i=1,nbins
        read(10,*) &
             z13Tablebb(i,imu), att13Tablebb(i,imu),  &
             pwc13Tablebb(i,imu), &
             z35Tablebb(i,imu), att35Tablebb(i,imu),  &
             kextTablebb(i,4,imu), salbTablebb(i,4,imu), &
             asymTablebb(i,4,imu),  &      
             kextTablebb(i,5,imu), salbTablebb(i,5,imu), &
             asymTablebb(i,5,imu),  &      
             pr13Tablebb(i,imu), d013TableBB(i,imu), &
             kextTablebb(i,2,imu), salbTablebb(i,2,imu), &
             asymTablebb(i,2,imu),  &      
             kextTablebb(i,3,imu), salbTablebb(i,3,imu), &
             asymTablebb(i,3,imu), &
             kextTablebb(i,6,imu), salbTablebb(i,6,imu), &
             asymTablebb(i,6,imu), &
             kextTablebb(i,7,imu), salbTablebb(i,7,imu), &
             asymTablebb(i,7,imu), &
             kextTablebb(i,1,imu), salbTablebb(i,1,imu), &  
             asymTablebb(i,1,imu)
        !print*, i, d013TableBB(i,1)
        salbTablebb(i,1:7,imu)=salbTablebb(i,1:7,imu)/kextTablebb(i,1:7,imu)
        salbTablebb(i,8,imu)=salbTablebb(i,7,imu)
        kextTablebb(i,8,imu)=kextTablebb(i,7,imu)
        asymTablebb(i,8,imu)=asymTablebb(i,7,imu)
        if(i>2) then
           if( att13TableBB(i-1,imu)>0 .and. att13TableBB(i,imu)<-90) then
              i0BBmax(imu)=i
           endif
        endif
        
        
     enddo
     do i=1,-nbins
       if(masterproc) print*, z13Tablebb(i,imu)-z35Tablebb(i,imu), z13Tablebb(i,imu), &
            pr13Tablebb(i,imu)
    enddo
    !stop
    close(10)
     
     
  enddo
  do i=1,nbins
     !write(*,*) kextTableS(i,5:8,1)
  enddo
  do k=7,-8
     kextTables(:,k,:)=kextTables(:,6,:)
     kextTable(:,k,:)=kextTable(:,6,:)
     kextTableBB(:,k,:)=kextTableBB(:,6,:)
     salbTables(:,k,:)=salbTables(:,6,:)
     salbTable(:,k,:)=salbTable(:,6,:)
     salbTableBB(:,k,:)=salbTableBB(:,6,:)
     asymTableS(:,k,:)=asymTableS(:,6,:)
     asymTable(:,k,:)=asymTable(:,6,:)
     asymTableBB(:,k,:)=asymTableBB(:,6,:)
  enddo
  if(masterproc)  write(*,*) mfreq
!  write(*,*) mfreq
!  stop
!  write(*,*) ibranch
!  stop
!  write(*,*) att35Table(:,1)
!  write(*,*) att35TableS(:,1)
!  write(*,*) att35TableBB(:,1)
!  stop
!  close(10)
  
  !stop
  
  !write(*,*) i0max
  !stop
  i0max=i0max-1
  !print*, i0max
  !stop
  i0BBmax=i0BBmax-1
  !print*, i0BBmax
10 format(30(F10.6))


end subroutine readTablesLiang2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


subroutine makehashTables()
 !use Tables2
 implicit none
 integer :: i,i0
 real    :: d0i
 d013Smax=maxval(d013TableS)
 d013BBmax=maxval(d013TableBB)
 d013max=maxval(d013Table)
 d013Smin=minval(d013TableS(1:i0max(1),1))
 d013BBmin=minval(d013TableBB(1:i0BBmax(1),1))
 d013min=minval(d013Table)
 do i=1,400
    d0i= d013Smin+(d013Smax-d013Smin)/400.*i
    call bisection2(d013TableS(:,1),i0max(1),d0i, i0)
    i0HTS(i)=i0
    d0i= d013BBmin+(d013BBmax-d013BBmin)/400.*i
    call bisection2(d013TableBB(:,1),i0BBmax(1),d0i, i0)
    !print*, d0i, d013TableBB(1:i0,1)
    i0HTBB(i)=i0
    d0i= d013min+(d013max-d013min)/400.*i
    call bisection2(d013Table(:,1),nbins,d0i, i0)
    i0HT(i)=i0
    !print*, i, i0HT(i), i0HTBB(i), i0HTS(i)
 enddo
 !print*, d013TableBB(1,:),  d013BBmin,   d013BBmax
 !print*, i0HTS
 !stop
end subroutine makehashTables

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


subroutine bisection2(a,n,x, ifind)
integer :: n, i, ifind
real    :: a(n), x
integer :: i1, i2, imid
i1=1
i2=n
if(x>=a(n)) then
   ifind=n
   return
endif
if(x<=a(1)) then
   ifind=1
   return
endif
do while (i2-i1>1)
   imid=(i1+i2)/2
   if(a(imid)>x) then
      i2=imid
   else
      if (a(imid)<x) then
         i1=imid
      else
         ifind=imid
         return
      endif
   endif
end do

ifind=i1

end subroutine bisection2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine getsnowp2(n0w,iwc,kext,salb,asym,nfreq)

  !use tables2
 real :: n0w  !log10(n0w) [log10(1/m4)]
 real :: iwc  !ice water content [g/m3]
 real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
 integer :: nfreq
 real  :: dpia13,dpia35,z35,z13,z94,dm
 integer :: imu, i0

 imu=3
 call bisection2(pwc13TableS2(:,imu),nbinS2,log10(iwc)-n0w, i0)

!5 mus (2,2,2,2,2) 8 frequencies (10,19,22,36,89,166,183.3+3,183.3+7)

 kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w     !10.**n0w = 10.**log10(nw) = nw 
 salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
 asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)

 dpia13=att13TableS2(i0,imu)*10.**n0w
 dpia35=att35TableS2(i0,imu)*10.**n0w
 z35=z35TableS2(i0,imu)+10*n0w
 !z94=z94TableS2(i0,imu)+10*n0w
 z13=zmin+i0*dzbin+10.*n0w
 dm=d013TableS2(i0,imu)

 return
 end subroutine getsnowp2

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine getsnowp2_sdsu(n0w, iwc, freq, kext_out, salb_out, asym_out)

 real :: n0w  !log10(n0w) [log10(1/m4)]
 real :: iwc  !ice water content [g/m3]
 real :: freq !frequency [GHz]
 real, intent(out) :: kext_out,salb_out, asym_out
 
 integer,parameter :: imu = 3 ! (-> mu=0)
 integer,parameter :: nfreq = 8
 real :: kext(nfreq),salb(nfreq), asym(nfreq)
 integer :: n

 integer :: i0
 real :: wgt1, wgt2


 if( n0w == 0. .or. iwc == 0. ) then
    kext_out=0. ; salb_out=0. ; asym_out=0.
    return
 endif

 call bisection2( pwc13TableS2(:,imu), nbinS2, log10(iwc)-n0w, i0 )

print*, 'i0=',i0, 'imu',imu


!5 mus (2,2,2,2,2) 8 frequencies (10,19,22,36,89,166,183.3+3,183.3+7)

 kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w     !10.**n0w = 10.**log10(nw) = nw 
 salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
 asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)

!print*,'freq=',freq, 'mfreq=',mfreq(:)

 if( freq < mfreq(1) ) then
    kext_out = kext(1) ; salb_out = salb(1) ; asym_out = asym(1)
 elseif( freq >= mfreq(nfreq) ) then
    kext_out = kext(nfreq) ; salb_out = salb(nfreq) ; asym_out = asym(nfreq)
 else

    freq_loop: do n = 1, nfreq-1
       if( mfreq(n) <= freq .and. freq <= mfreq(n+1) ) then
          wgt1 = ( mfreq(n+1) - freq ) / ( mfreq(n+1) - mfreq(n) )
          wgt2 = 1.-wgt1 
          if(wgt2<0. .or. wgt2>1.) stop 'MSG getsnowp2_sdsu: strange wgt2'
          kext_out = wgt1 * kext(n) + wgt2 * kext(n+1) 
          salb_out = wgt1 * salb(n) + wgt2 * salb(n+1)
          asym_out = wgt1 * asym(n) + wgt2 * asym(n+1)
!print*, wgt1, wgt2

          exit freq_loop
       endif
    enddo freq_loop
   
 endif

 return
 end subroutine getsnowp2_sdsu

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_openssp_cmb

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
