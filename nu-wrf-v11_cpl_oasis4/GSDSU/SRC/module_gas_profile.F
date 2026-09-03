!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_gas_profile
 use module_floating_point
 use module_const
 use module_simulator
 implicit none

!
! initialize (reading) gaseous profiles.
!

 private  ! all are encapsulate.

 public :: init_gas_profile  , &! initialize gaseous profiles of climatology
           interp_gas_profile   ! interpolate gas concentration 

 save     ! all parameters are saved.

!
!
!
 integer :: nmax_profile    ! number of standard atmospheric condition
 integer :: kmax_gas        ! max number of sounding toshii
 integer :: nmax_majgas     ! number of major gas
 integer :: nmax_mingas     ! number of minor gas
 integer :: nmax_totgas     ! total number of gaseous species

! for mlatm
 character(len=20),allocatable :: id (:)  ! gas id
 integer,allocatable :: idms(:,:)
 real(sdsu_fps),allocatable :: wmol(:,:)    ! Molecular weight [g/mol] 
 real(sdsu_fps),allocatable :: rams(:,:)  
 real(sdsu_fps),allocatable :: alt(:)       ! altitude [km]
 real(sdsu_fps),allocatable :: pres(:,:)    ! air pressure
 real(sdsu_fps),allocatable :: temp(:,:)    ! air temperature
 real(sdsu_fps),allocatable :: dnsty(:,:)   ! air density  
 real(sdsu_fps),allocatable :: gas_major(:,:,:)  ! major gas conc [ppmv]
 real(sdsu_fps),allocatable :: gas_minor(:,:)    ! minor gas conc [ppmv]

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine init_gas_profile
 implicit none


!--- work area
 integer,parameter ::io=13
 logical :: logic_opened
 integer :: i,j,k,l,m, ns
 integer :: ierr
 real(sdsu_fps) :: airm  ! air molecular weight [ ] 
 character(len=1) :: dummy ! dummy character argument
 character(len=200) :: input_dir_file     ! file name

!
! initilization start
!

 input_dir_file = trim(sdsu_dir_data)//'MLATMD'
 inquire(file=trim(input_dir_file),opened=logic_opened)

 if( .not. logic_opened)  then
    if(masterproc .and. verbose_SDSU) print*,'MSG mlatm: open/read MLATMD'

    call open_oldfile(io, trim(input_dir_file) )
    
    !
    ! read dimension indices 
    !
    read(io,'(a1)') dummy 
    read(io,*) nmax_majgas, nmax_mingas, nmax_profile, kmax_gas  ! dimension indices
    nmax_totgas = nmax_majgas + nmax_mingas

    !
    ! allocate all parameters
    !
    if( .not. allocated(alt) ) &
      allocate( &
            id   (nmax_totgas   ),&   ! gas id
            idms (nmax_totgas,10),&
            wmol (nmax_totgas,10),&   !Molecular weight [g/mol]
            rams (nmax_totgas,10),&
            alt  (kmax_gas),&         ! altitude [km]
            pres (kmax_gas,            nmax_profile),& ! 
            temp (kmax_gas,            nmax_profile),& ! 
            dnsty(kmax_gas,            nmax_profile),& ! 
            gas_major (kmax_gas,nmax_majgas,nmax_profile),& ! major gas conc [ppmv]
            gas_minor (kmax_gas,nmax_mingas),&  ! minor gas conc []
            stat=ierr )
    if (ierr /= 0) stop 'MSG init_gas_profile: allocation error > Terminate program. '

    !
    ! read molecular information (not used)
    !
    read(io,'(a1)') dummy 
    read(io,*) airm  ! air molecular weight [ ]
    read(io,'(a1)') dummy 

    !
    ! gas information (not used)
    !
    idms(:,:)=0
    rams(:,:)=0.0d0
    do i = 1, nmax_totgas
       read(io,*) id(i)
       read(io,*) ns
       read(io,*) idms(i,1:ns)  !??
       read(io,*) wmol(i,1:ns)  ! Molecular weight [g/mol]
       read(io,*) rams(i,1:ns)  !??
    enddo

    !
    ! read altitude
    !
    read(io,'(a1)') dummy 
    read(io,*) alt(1:kmax_gas)

    !
    ! read pressure
    !
    do j = 1, nmax_profile
       read(io,'(a1)') dummy 
       read(io,*) pres(1:kmax_gas,j)
    enddo

    !
    ! read temperature
    !
    do j = 1, nmax_profile
       read(io,'(a1)') dummy 
       read(io,*) temp(1:kmax_gas,j)
    enddo

    !
    ! read molecular profiles
    !
    do k = 1, nmax_majgas
       do j=1, nmax_profile
          read(io,'(a1)') dummy 
          read(io,*) gas_major(1:kmax_gas,k,j)  ! [PPMV]
       enddo
    enddo

    do j=1, nmax_profile
       read(io,'(a1)') dummy 
       read(io,*) dnsty(1:kmax_gas,j)  ! [AIR MOLECULES / cm3]
    enddo

    !
    ! read minor gase profiles
    !
    do  k=1,nmax_mingas
       read(io,'(a1)') dummy 
       read(io,*) gas_minor(1:kmax_gas,k)
    enddo

 endif !logic_opened

 return
 end subroutine init_gas_profile

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine interp_gas_profile( month, latitude, p_in, gas_id, unit_id, gas_conc )
 implicit none
 integer,intent(in) :: month               ! month (1 ~12)
 real(sdsu_fps),intent(in)    :: latitude            ! latitude [deg] (-90 ~ 90)
 real(sdsu_fps),dimension(:),intent(in)  :: p_in     ! input pressure profile [mb]
 character(len=*),intent(in)   :: gas_id   ! character gas id
 character(len=*),intent(in)   :: unit_id  ! unit id 'ppmv' or 'kmol/cm3'

 real(sdsu_fps),dimension(:),intent(out) :: gas_conc ! output gas concentration 

 integer :: ip
 integer :: k, k_gas, n, n_gas
 integer :: n_gas_major, n_gas_minor
 integer :: kmax
 logical :: find_gas
 character(len=5) :: gas_type
 real(sdsu_fps) :: wgt1, wgt2
 real(sdsu_fps) :: p_wgt, t_wgt, gas_wgt
 real(sdsu_fps),parameter :: avo = 6.02214199e+23  ! Avogadro number



!
! filter input
!
 if( month < 1 .or. month > 12 ) then
    print*, 'MSG interp_gas_profile: bad month value =',month ; stop
 elseif( latitude < -90.0 .or. latitude > 90.0 ) then
    print*, 'MSG interp_gas_profile: bad latitude value =',latitude ; stop
 endif


!
! identify gas species
!
 gas_loop: do n = 1, nmax_totgas
    find_gas = .false.
    if( trim(id(n)) == trim(gas_id) ) then
      find_gas = .true.
      n_gas = n
!      print*,'Gas ID = ', trim(id(n)) 
      exit gas_loop
    endif
 enddo gas_loop

 if( .not. find_gas ) then
    print*,'MSG interp_gas_profile: no such gas id =', trim(gas_id)
    stop
 endif


!
! major or minor gas
!
 if( n_gas <= nmax_majgas ) then
    gas_type = 'major'
    n_gas_major = n_gas
 else
    gas_type = 'minor'
    n_gas_minor = n_gas - nmax_majgas
 endif



!
! Find climatology index: ip
!
 select case(gas_type)
 case('major')

  !
  ! For major gas (region and season dependent)
  !  1: tropical
  !  2: mid-latitude summer    3: mid-latitude winter
  !  4: high-lat summer        5: high-lat winter

 if( latitude >= 60. ) then
   ! NH high latitude  
   if( month >= 5 .and. month <= 10 ) then
     ! summer
     ip = 4
   else
     ip = 5
   endif
 elseif( latitude > 30. .and. latitude < 60.  ) then
   ! NH mid latitude
   if( month >= 5 .and. month <= 10 ) then
     ! summer
     ip = 2 
   else
     ip = 3
   endif
 elseif( abs(latitude) <= 30. ) then  
   ! tropics
   ip = 1
 elseif( latitude > -60. .and. latitude < -30.  ) then
   ! SH mid latitude
   if( month >= 5 .and. month <= 10 ) then
     ! winter
     ip = 3 
   else
     ip = 2 
   endif
 elseif( latitude <= -60. ) then
   if( month >= 5 .and. month <= 10 ) then
     ! winter
     ip = 5 
   else
     ip = 4 
   endif
 endif 

 case('minor')

  !
  ! For minor gas (use US stanradr profile; indepdent from season and region)
  !  6: us standard
  !
  ip = 6

 end select


!
! interpolate pressure, temperature, and gaseous concentrations from input pressure. 
!
 kmax = size(p_in) ! vertical dimension index
 do k = 1, kmax

    if ( p_in(k) < pres(kmax_gas, ip) ) then ! TOA
        p_wgt  = pres(kmax_gas, ip)   ! weighted pressure [mb]
        t_wgt  = temp(kmax_gas, ip)   ! weighted temperature [K]
        select case(gas_type)
        case('major')
                 gas_wgt = gas_major (kmax_gas,n_gas_major,ip)   ! major gas conc [ppmv]
        case('minor')
                 gas_wgt = gas_minor (kmax_gas,n_gas_minor)      ! minor gas conc [ppmv]
        end select

    elseif( p_in(k) > pres(1,        ip) ) then ! BOA
        p_wgt  = pres(       1, ip)   ! weighted pressure [mb]
        t_wgt  = temp(       1, ip)   ! weighted temperature [K]
        select case(gas_type)
        case('major')
                 gas_wgt = gas_major (1,n_gas_major,ip)   ! major gas conc [ppmv]
        case('minor')
                 gas_wgt = gas_minor (1,n_gas_minor)      ! minor gas conc [ppmv]
        end select

    else  ! interpolate
       do k_gas = 1, kmax_gas-1
          if( p_in(k) <=  pres(k_gas, ip) .and. p_in(k) > pres(k_gas+1, ip) ) then
              wgt1 = (p_in(k) - pres(k_gas+1, ip) ) / ( pres(k_gas, ip) - pres(k_gas+1, ip) )
              wgt2 = 1.0 - wgt1 
              if( wgt2 < 0.0 .or. wgt2 > 1.0 ) stop 'MSG interp_gas_profile: strange wgt2'
              p_wgt = wgt1* pres(k_gas, ip) + wgt2 * pres(k_gas+1, ip)  !weighted pressure [mb]
              t_wgt = wgt1* temp(k_gas, ip) + wgt2 * temp(k_gas+1, ip)  !weighted temperature [K]

              select case(gas_type)
              case('major')
                 gas_wgt =  wgt1 * gas_major (k_gas  ,n_gas_major,ip) &     ! major gas conc [ppmv]
                          + wgt2 * gas_major (k_gas+1,n_gas_major,ip)
              case('minor')   
                 gas_wgt =  wgt1 * gas_minor (k_gas  ,n_gas_minor) &     ! minor gas conc [ppmv]
                          + wgt2 * gas_minor (k_gas+1,n_gas_minor)
              end select

          endif

       enddo

    endif


    !
    ! output unit options....
    !
    select case(trim(unit_id)) 
    case('ppmv') 
      gas_conc(k) = gas_wgt    ! interpolated gaseous concentration [ppmv]
    case('kmol/cm3') 
      gas_conc(k) = gas_wgt*1.e-6 * (p_wgt * 0.01 )/(8.315*t_wgt) * 1.e-9 * 1.e4  ! [kilo mole/cm3]
    case default
      print*,'MSG interp_gas_profile: no such unit',trim(unit_id)
      stop
    end select

 enddo


 return
 end subroutine interp_gas_profile

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_gas_profile

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
