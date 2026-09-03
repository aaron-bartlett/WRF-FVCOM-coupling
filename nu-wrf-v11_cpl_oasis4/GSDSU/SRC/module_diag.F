#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_diag
 use module_floating_point
 use module_const
 use module_simulator
 use module_mpi
 implicit none

!
! initialize (reading) gaseous profiles.
!

 save     ! all parameters are saved.


 interface dmain_mean_profile
    module procedure dmain_mean_profile_gen
 end interface dmain_mean_profile


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine diag_domain_mask
 implicit none
! create domain mask for analysis and simulator
 integer :: i,j

!
! create domain mask
!
 surface%diag_mask = 0
 
 do j = myj_start, myj_end ; do i = myi_start, myi_end  
    
    select case(trim(mask_method))
    case('GRID')
     
      if( is_diag <= i .and. i <= ie_diag ) then 
        if(  js_diag <= j .and. j <= je_diag ) then

             surface(i,j)%diag_mask = 1

        endif
      endif
    
    case('LATLON')

      if( minlat_diag <= surface(i,j)%lat .and. surface(i,j)%lat <= maxlat_diag ) then
        if( minlon_diag <= surface(i,j)%lon .and. surface(i,j)%lon <= maxlon_diag ) then

             surface(i,j)%diag_mask = 1

        endif
      endif
    case('NONE')

             surface(i,j)%diag_mask = 1

    case default
      print*,'mask_method=',mask_method
      stop 'MSG: diag_domain_mask: no such option'
    end select


 enddo ; enddo

 return
 end subroutine diag_domain_mask

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine dmain_mean_profile_gen(proc, kmax_in, q3d, n3d, w, hgt_stag, dhgt_crm, &
                sctot, lh, cldnuc, rhodair, mask, iclass ) 
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive domain-mean vertical profile
!           
! History:
! 04/2019  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 character*(*) :: proc  ! process
 integer,intent(in) :: kmax_in
 type ( particle_gen ), dimension(:,:,:), optional, intent(in) :: q3d, n3d ! parameter assumed 3d array
 real, dimension(:,:,:), optional, intent(in) :: w    ! vertical velocity [m/s]
 real, dimension(:,:,:), optional, intent(in) :: hgt_stag ! staggered height [km]
 real, dimension(:,:,:), optional, intent(in) :: dhgt_crm ! thickness of layer [km]

 real, dimension(:,:,:), optional, intent(in) :: sctot !total space charge density [C/m3]  NSSL microphysics only..
 real, dimension(:,:,:), optional, intent(in) :: lh !latent heat [K/day] .
 real, dimension(:,:,:), optional, intent(in) :: cldnuc !cloud nucleation rate [kg/kg/s] P3 only..

 real, dimension(:,:,:), optional, intent(in) :: rhodair !air density [g/m3]
 integer, dimension(:,:), optional, intent(in) :: mask    !0-no analysis 1-analysi
 integer, dimension(:,:), optional, intent(in) :: iclass  !0-undefined, 1-low, 2-mid, 3-high, 4-convective, 5-deep stratiform

 integer :: bnd(3)  !upper memory bound for 2D array

 character(len=10) :: case_name

 integer :: i,j,k,kk,np,nclass !loop indice
 integer :: kmax, mxlyr_crm
 integer :: irec
 integer,parameter :: io = 111

 real :: tot_conden  !total condensate [g/kg]
 integer,parameter :: nmax = 5     !total cloud category
 integer,parameter :: np_max = 16  !total physical parameter

! integer,parameter :: mxlyr_inst = 20*4   ! every 250m for 20km
! real,parameter :: dH = 0.25  !thickness of layer

 integer,parameter :: mxlyr_inst = 20   ! every 1km for 20km
 real,parameter :: dH = 1.0  !thickness of layer

 real :: hgt_stag_inst(0:mxlyr_inst)  !staggered height level [km]

 real,allocatable :: bundle_crm(:,:)  !bundled parameters at CRM level
 real :: bundle_inst(np_max,mxlyr_inst) ! bundled parameters at instrument level

 real,dimension(0:nmax,1:np_max,1:mxlyr_inst),save :: domain_sum, domain_mean
 real,dimension(nmax,mxlyr_inst),save :: domain_cnt

 real,dimension(0:nmax,1:np_max),save :: domain_column_sum, domain_column_mean
 real,dimension(nmax),save :: domain_column_cnt

 real :: hgt_stag_crm(0:kmax_in)
 real :: tot_cnt
 real :: w_inst
 real :: hgt_crm
 
 logical,parameter :: new_method = .true.  ! new method better for coarese vertical resolution

 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

  if(masterproc) print*,'MSG dmain_mean_profile_gen: initialize parameters '

  domain_cnt  = 0.
  domain_sum  = 0.
  domain_mean = 0.

  domain_column_cnt  = 0.
  domain_column_sum  = 0.
  domain_column_mean = 0.

 !
 ! initialize constant height
 !
  do k = 0, mxlyr_inst
    hgt_stag_inst(k) = real(k)*dH       !20x4 layers (20km for every 250m)
  enddo

 case('add') !----------------------------------------------------------------------------------

  if(masterproc) print*, 'MSG dmain_mean_profile_gen: add'


  if(present(q3d) ) then

   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(q3d)
   mxlyr_crm = kmax_in

   if( .not. allocated(bundle_crm) ) then
      allocate ( bundle_crm(np_max, mxlyr_crm) )
   endif

   do j = 1,bnd(2) ; do i = 1,bnd(1)

      !skip analysis if outside of domain or no defined cloud class
      if( mask(i,j) == 0 .or. iclass(i,j) == 0) cycle

      np = 1
      bundle_crm(np,1:mxlyr_crm) = w    (i,j,1:mxlyr_crm)          ; np = np+1
      bundle_crm(np,1:mxlyr_crm) = sctot(i,j,1:mxlyr_crm)          ; np = np+1
      bundle_crm(np,1:mxlyr_crm) = lh   (i,j,1:mxlyr_crm)/(24. * 3600.)  ; np = np+1   ![K/s]
      bundle_crm(np,1:mxlyr_crm) = cldnuc(i,j,1:mxlyr_crm)          ; np = np+1

      bundle_crm(np,1:mxlyr_crm) = q3d(i,j,1:mxlyr_crm)%cloud / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [g/kg] 
      bundle_crm(np,1:mxlyr_crm) = q3d(i,j,1:mxlyr_crm)%ice   / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [g/kg] 
      bundle_crm(np,1:mxlyr_crm) = q3d(i,j,1:mxlyr_crm)%rain  / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [g/kg] 
      bundle_crm(np,1:mxlyr_crm) = q3d(i,j,1:mxlyr_crm)%snow  / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [g/kg] 
      bundle_crm(np,1:mxlyr_crm) = q3d(i,j,1:mxlyr_crm)%graupel/rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [g/kg] 
      bundle_crm(np,1:mxlyr_crm) = q3d(i,j,1:mxlyr_crm)%hail  / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [g/kg] 
      bundle_crm(np,1:mxlyr_crm) = n3d(i,j,1:mxlyr_crm)%cloud / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [#/kg] 
      bundle_crm(np,1:mxlyr_crm) = n3d(i,j,1:mxlyr_crm)%ice   / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [#/kg] 
      bundle_crm(np,1:mxlyr_crm) = n3d(i,j,1:mxlyr_crm)%rain  / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [#/kg] 
      bundle_crm(np,1:mxlyr_crm) = n3d(i,j,1:mxlyr_crm)%snow  / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [#/kg] 
      bundle_crm(np,1:mxlyr_crm) = n3d(i,j,1:mxlyr_crm)%graupel/rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [#/kg] 
      bundle_crm(np,1:mxlyr_crm) = n3d(i,j,1:mxlyr_crm)%hail  / rhodair(i,j,1:mxlyr_crm) *1.e3 ; np = np+1 ! [#/kg] 


     nclass = iclass(i,j) !column class 1-low, 2-mid, 3-high, 4-convective, 5-deep stratiform

     ! column integrated statistics
     do k = 1, mxlyr_crm
        tot_conden = sum( bundle_crm(5:10,k) )  ! [g/kg]
        w_inst = bundle_crm(1,k)   !vertical velocity [m/s]
!        if( tot_conden > 0.1 .and. .not. isnan(tot_conden) ) then !ACPCin-cloud definition
        if( tot_conden > 0.1 .and. .not. isnan(tot_conden) .and. w_inst >= 1.0 ) then !ACPC with vertical velocity
            domain_column_cnt(nclass) = domain_column_cnt(nclass) + 1.   !count for different cloud class

            domain_column_sum(nclass,1:np_max) = domain_column_sum(nclass,1:np_max) &
                                            + bundle_crm(1:np_max,k)  !different cloud
            domain_column_sum(0,1:np_max) = domain_column_sum(0,1:np_max) &
                                            + bundle_crm(1:np_max,k)  !all cloud
        endif

     enddo



   if( new_method ) then  !========  new method better for coarese vertical resolution

     do k = 0, mxlyr_crm
        hgt_stag_crm(k) = hgt_stag(i,j,k+1)
     enddo

     do k = 1,mxlyr_crm

        hgt_crm = 0.5*( hgt_stag_crm(k-1) + hgt_stag_crm(k) ) !km

        kk = int(hgt_crm) + 1   !toshii temporal


        if(kk < 1 .or. kk > 20) cycle

        tot_conden = sum( bundle_crm(5:10,k) )  ! [g/kg]
        w_inst = bundle_crm(1,k)   !vertical velocity [m/s]
!        if( tot_conden > 0.1 .and. .not. isnan(tot_conden) ) then !ACPCin-cloud definition
        if( tot_conden > 0.1 .and. .not. isnan(tot_conden) .and. w_inst >= 1.0 ) then !ACPC with vertical velocity

!          do kk = 0, mxlyr_inst-1
!             if( hgt_stag_inst(kk) < hgt_crm .and. hgt_crm <= hgt_stag_inst(kk+1) ) then
!                 print*, 'kk=',kk, hgt_stag_inst(kk), hgt_crm, hgt_stag_inst(kk+1)
 
                domain_cnt(nclass,kk) = domain_cnt(nclass,kk) + 1.   !count for different cloud class

                domain_sum(nclass,1:np_max,kk) = domain_sum(nclass,1:np_max,kk) + bundle_crm(1:np_max,k)  !different cloud

                domain_sum(     0,1:np_max,kk) = domain_sum(     0,1:np_max,kk) + bundle_crm(1:np_max,k)  !all cloud
 
!             endif 
!          enddo 

        endif

     enddo



   else ! =========  old method (better for finer vertical resolution)
     !
     ! refine vertical layer
     !
     do k = 0, mxlyr_crm
        hgt_stag_crm(k) = hgt_stag(i,j,k+1)
     enddo
     call refine_range_bin( mxlyr_crm, hgt_stag_crm(0:mxlyr_crm), dhgt_crm(i,j,1:mxlyr_crm), &
                            mxlyr_inst, hgt_stag_inst, &
                            np_max, bundle_crm(1:np_max,1:mxlyr_crm), bundle_inst(1:np_max,1:mxlyr_inst) )


     ! profile statistics
     do k = 1, mxlyr_inst

        tot_conden = sum( bundle_inst(5:10,k) )  ! [g/kg]
        w_inst = bundle_inst(1,k)   !vertical velocity [m/s]


        if( tot_conden > 0.1 .and. .not. isnan(tot_conden) ) then !ACPCin-cloud definition
! ACPCin-cloud definition with vertical velocity 
!        if( tot_conden > 0.1 .and. .not. isnan(tot_conden) .and. w_inst >= 1.0 ) then !ACPC with vertical velocity

            domain_cnt(nclass,k) = domain_cnt(nclass,k) + 1.   !count for different cloud class

            domain_sum(nclass,1:np_max,k) = domain_sum(nclass,1:np_max,k) + bundle_inst(1:np_max,k)  !different cloud

            domain_sum(     0,1:np_max,k) = domain_sum(     0,1:np_max,k) + bundle_inst(1:np_max,k)  !all cloud

        endif

     enddo 

   endif !new_method


   enddo ; enddo



  endif


 case('finish') !-------------------------------------------------------------------------------

  if(masterproc) print*,'MSG dmain_mean_profile_gen: finalize parameters '

#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', domain_cnt        )
  call mpi_sdsu_reduce( 'SUM', domain_sum      )

  call mpi_sdsu_reduce( 'SUM', domain_column_cnt        )
  call mpi_sdsu_reduce( 'SUM', domain_column_sum      )

#endif

  MPI_MASTERPROC: if(masterproc) then

!
! calcualte domain mean column integrated values
!
  do nclass = 1, nmax

    if( domain_column_cnt(nclass) > 0. ) then
        domain_column_mean(nclass,1:np_max) = domain_column_sum(nclass,1:np_max) / domain_column_cnt(nclass)
    else
        domain_column_mean(nclass,1:np_max) = undefined
    endif

  enddo

    tot_cnt = sum(domain_column_cnt(1:nmax) )
    if( tot_cnt > 0. ) then
        domain_column_mean(0,1:np_max) = domain_column_sum(0,1:np_max) / tot_cnt
    else
        domain_column_mean(0,1:np_max) = undefined
    endif

!
! calcualte domain mean profile
!
  do k = 1, mxlyr_inst ; do nclass = 1, nmax

    if( domain_cnt(nclass,k) > 0. ) then
        domain_mean(nclass,1:np_max,k) = domain_sum(nclass,1:np_max,k)   / domain_cnt(nclass,k)
    else
        domain_mean(nclass,1:np_max,k) = undefined
    endif

  enddo ; enddo

  do k = 1, mxlyr_inst 

    tot_cnt = sum(domain_cnt(1:nmax,k) ) 
    if( tot_cnt > 0. ) then
        domain_mean(0,1:np_max,k) = domain_sum(0,1:np_max,k)   / tot_cnt
    else
        domain_mean(0,1:np_max,k) = undefined
    endif

  enddo 


!
! bundle output
!

 file_loop: do nclass = 0, nmax   

   select case(nclass)
   case(0) ; case_name = 'tot'
   case(1) ; case_name = 'low'
   case(2) ; case_name = 'mid'
   case(3) ; case_name = 'hig'
   case(4) ; case_name = 'conv'
   case(5) ; case_name = 'dstr'
   end select

!
! write out domain profile
!
 sdsu_io_file = trim(sdsu_dir_output)//'hydrometeor_profile_'//trim(case_name)//trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=mxlyr_inst*4)
 print("(a)"),'MSG dmain_mean_profile_gen: see output ->',trim(sdsu_io_file)

 irec = 1
 do np = 1, np_max !parameer loop
    write(io,rec=irec) ( domain_mean(nclass, np, k) , k=1,mxlyr_inst) ; irec = irec + 1
 enddo
 close(io)


 if( write_grads_ctl  ) then
    sdsu_io_file = trim(sdsu_dir_output)//'hydrometeor_profile_'//trim(case_name)//trim(output_suffix)//'.ctl'
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = 'hydrometeor_profile_'//trim(case_name)//trim(output_suffix)//'.grads'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  domain mean  hydrometeor profiles '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   1 linear 1. ',gridsize
    write(io,*) 'ydef   1 linear 1. ',gridsize
    write(io,*) 'zdef   ',mxlyr_inst  ,' linear 0. ',dH
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',irec-1

    write(io,*) 'w    ',mxlyr_inst,'  0  vertical velocity [m/sec]'
    write(io,*) 'sctot',mxlyr_inst,'  0  total space charge density [C/m3]'
    write(io,*) 'LH'   ,mxlyr_inst,'  0  latent heat  [K/s]'
    write(io,*) 'cldnuc',mxlyr_inst,'  0  cloud nucleation rate [kg/kg/s]'
    write(io,*) 'qc   ',mxlyr_inst,'  0  cloud mass conc. [g/kg]'
    write(io,*) 'qi   ',mxlyr_inst,'  0  ice mass conc. [g/kg]'
    write(io,*) 'qr   ',mxlyr_inst,'  0  rain mass conc. [g/kg]'
    write(io,*) 'qs   ',mxlyr_inst,'  0  snow mass conc. [g/kg]'
    write(io,*) 'qg   ',mxlyr_inst,'  0  graupel mass conc. [g/kg]'
    write(io,*) 'qh   ',mxlyr_inst,'  0  hail mass conc. [g/kg]'
    write(io,*) 'nc   ',mxlyr_inst,'  0  cloud number conc. [#/m3]'
    write(io,*) 'ni   ',mxlyr_inst,'  0  ice number conc. [#/m3]'
    write(io,*) 'nr   ',mxlyr_inst,'  0  rain number conc. [#/m3]'
    write(io,*) 'ns   ',mxlyr_inst,'  0  snow number conc. [#/m3]'
    write(io,*) 'ng   ',mxlyr_inst,'  0  graupel number conc. [#/m3]'
    write(io,*) 'nh   ',mxlyr_inst,'  0  hail number conc. [#/m3]'

    write(io,*) 'endvars '
    close(io)
 endif

!
! write out domain column integrated
!
 sdsu_io_file = trim(sdsu_dir_output)//'hydrometeor_column_'//trim(case_name)//trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=4)
 print("(a)"),'MSG dmain_mean_profile_gen: see output ->',trim(sdsu_io_file)

 irec = 1
 do np = 1, np_max !parameer loop
    write(io,rec=irec) domain_column_mean(nclass, np)  ; irec = irec + 1
 enddo

 close(io)


 if( write_grads_ctl  ) then
    sdsu_io_file = trim(sdsu_dir_output)//'hydrometeor_column_'//trim(case_name)//trim(output_suffix)//'.ctl'
    open(io,file=trim(sdsu_io_file))
    sdsu_io_file = 'hydrometeor_column_'//trim(case_name)//trim(output_suffix)//'.grads'

    write(io,"(a)") 'dset  '//trim(sdsu_io_file)
    write(io,*) 'title  domain mean  hydrometeor profiles '
    write(io,*) 'undef  ',undefined
    write(io,*) 'xdef   1 linear 1. ',gridsize
    write(io,*) 'ydef   1 linear 1. ',gridsize
    write(io,*) 'zdef   1 linear 0. 1 '
    write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
    write(io,*) 'vars ',irec-1

    write(io,*) 'w        0  0  vertical velocity [m/sec]'
    write(io,*) 'sctot    0  0  total space charge density [C/m3]'
    write(io,*) 'LH       0  0  latent heat  [K/s]'
    write(io,*) 'cldnuc   0  0  cloud nucleation rate [kg/kg/s]'
    write(io,*) 'qc       0  0  cloud mass conc. [g/kg]'
    write(io,*) 'qi       0  0  ice mass conc. [g/kg]'
    write(io,*) 'qr       0  0  rain mass conc. [g/kg]'
    write(io,*) 'qs       0  0  snow mass conc. [g/kg]'
    write(io,*) 'qg       0  0  graupel mass conc. [g/kg]'
    write(io,*) 'qh       0  0  hail mass conc. [g/kg]'
    write(io,*) 'nc       0  0  cloud number conc. [#/m3]'
    write(io,*) 'ni       0  0  ice number conc. [#/m3]'
    write(io,*) 'nr       0  0  rain number conc. [#/m3]'
    write(io,*) 'ns       0  0  snow number conc. [#/m3]'
    write(io,*) 'ng       0  0  graupel number conc. [#/m3]'
    write(io,*) 'nh       0  0  hail number conc. [#/m3]'

    write(io,*) 'endvars '
    close(io)
 endif

 enddo file_loop

 endif MPI_MASTERPROC


 case default ; call stop_sdsu(myrank,'MSG dmain_mean_profile_gen: There is no such proc')
 end select

 return
 end subroutine dmain_mean_profile_gen

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


 subroutine spatial_mean_2d(proc,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive spatial mean map of 2D parameters.
!           
! History:
! 11/2016  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 character*(*) :: proc ! process
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 2d array

 integer :: bnd(2)  !upper memory bound for 2D array

!
! statisticsl parameters
!
 type( surface_parameter ), allocatable, dimension(:,:),save :: spatial_sum, spatial_mean
 real(sdsu_fps), allocatable, dimension(:,:),save :: spatial_cnt
 real(sdsu_fps), allocatable, dimension(:),save :: spatial_lat, spatial_lon

 integer,save :: ibin_max
 integer :: i,j, ibin, ierr, irec
 integer,parameter :: io = 149

 integer,save :: imax, jmax
 integer :: ii, jj

 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

   if(masterproc) print*, 'MSG spatial_mean_2d: init'

   imax = int( ( maxlon - minlon ) / res_latlon ) + 1
   jmax = int( ( maxlat - minlat ) / res_latlon ) + 1

   if( imax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_2d: maxlon is smaller than minlon')
   if( jmax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_2d: maxlat is smaller than minlat')

  
   if( allocated(spatial_sum) ) deallocate(spatial_sum, spatial_cnt, spatial_mean, spatial_lat, spatial_lon)
 
   allocate( spatial_sum(imax,jmax)  , &
             spatial_cnt(imax,jmax)  , &
             spatial_mean(imax,jmax) , &
             spatial_lat(jmax)  , &
             spatial_lon(imax)  , &
             stat=ierr )

   ! initialized
   do j = 1, jmax
      spatial_lat(j) = minlat + res_latlon * real(j-1) ![deg]
   enddo
   do i = 1, imax
      spatial_lon(i) = minlon + res_latlon * real(i-1) ![deg]
   enddo

  !
  ! zero out other parameters
  !
  spatial_cnt = 0.

  spatial_sum%rain_rate      = 0. ; spatial_sum%h2o_col       = 0.
  spatial_sum%path_cloud_liq = 0. ; spatial_sum%path_rain_liq = 0.
  spatial_sum%path_cloud_ice = 0. ; spatial_sum%path_rain_ice = 0.
  spatial_sum%toa_lw_up      = 0. ; spatial_sum%toa_lw_dn     = 0.
  spatial_sum%toa_sw_up      = 0. ; spatial_sum%toa_sw_dn     = 0.
  spatial_sum%boa_lw_up      = 0. ; spatial_sum%boa_lw_dn     = 0.
  spatial_sum%boa_sw_up      = 0. ; spatial_sum%boa_sw_dn     = 0.

  spatial_mean%rain_rate      = undefined ; spatial_mean%h2o_col       = undefined
  spatial_mean%path_cloud_liq = undefined ; spatial_mean%path_rain_liq = undefined
  spatial_mean%path_cloud_ice = undefined ; spatial_mean%path_rain_ice = undefined
  spatial_mean%toa_lw_up      = undefined ; spatial_mean%toa_lw_dn     = undefined
  spatial_mean%toa_sw_up      = undefined ; spatial_mean%toa_sw_dn     = undefined
  spatial_mean%boa_lw_up      = undefined ; spatial_mean%boa_lw_dn     = undefined
  spatial_mean%boa_sw_up      = undefined ; spatial_mean%boa_sw_dn     = undefined


 case('add') !-------------------------------------------------------------------------------------

  if(masterproc) print*, 'MSG spatial_mean_2d: add'

  if(present(var2d) ) then

   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var2d)

   do j = 1,bnd(2) ; do i = 1,bnd(1)

      !
      ! loop for latlon_grid domain. 
      !
      do jj = 1, jmax
        within_grid: if( var2d(i,j)%lat >  (spatial_lat(jj)-0.5*res_latlon) .and. &
                         var2d(i,j)%lat <= (spatial_lat(jj)+0.5*res_latlon) ) then

!        if(masterproc) print*, 'MSG spatial_mean_2d: found lat ', (spatial_lat(jj)-0.5*res_latlon), &
!         var2d(i,j)%lat, (spatial_lat(jj)+0.5*res_latlon)

        do ii = 1, imax
          if( var2d(i,j)%lon >  (spatial_lon(ii)-0.5*res_latlon) .and. &
              var2d(i,j)%lon <= (spatial_lon(ii)+0.5*res_latlon) ) then

!             if(masterproc) print*, 'MSG spatial_mean_2d: found lon ',(spatial_lon(ii)-0.5*res_latlon),&
!                var2d(i,j)%lon, (spatial_lon(ii)+0.5*res_latlon)
 
             spatial_cnt(ii,jj) = spatial_cnt(ii,jj) + 1.e0  !update total count

             spatial_sum(ii,jj)%rain_rate      = spatial_sum(ii,jj)%rain_rate + var2d(i,j)%rain_rate
             spatial_sum(ii,jj)%h2o_col        = spatial_sum(ii,jj)%h2o_col   + var2d(i,j)%h2o_col
             spatial_sum(ii,jj)%path_cloud_liq = spatial_sum(ii,jj)%path_cloud_liq + var2d(i,j)%path_cloud_liq
             spatial_sum(ii,jj)%path_rain_liq  = spatial_sum(ii,jj)%path_rain_liq  + var2d(i,j)%path_rain_liq
             spatial_sum(ii,jj)%path_cloud_ice = spatial_sum(ii,jj)%path_cloud_ice + var2d(i,j)%path_cloud_ice
             spatial_sum(ii,jj)%path_rain_ice  = spatial_sum(ii,jj)%path_rain_ice  + var2d(i,j)%path_rain_ice
             spatial_sum(ii,jj)%toa_lw_up      = spatial_sum(ii,jj)%toa_lw_up + var2d(i,j)%toa_lw_up
             spatial_sum(ii,jj)%toa_lw_dn      = spatial_sum(ii,jj)%toa_lw_dn + var2d(i,j)%toa_lw_dn
             spatial_sum(ii,jj)%toa_sw_up      = spatial_sum(ii,jj)%toa_sw_up + var2d(i,j)%toa_sw_up
             spatial_sum(ii,jj)%toa_sw_dn      = spatial_sum(ii,jj)%toa_sw_dn + var2d(i,j)%toa_sw_dn
             spatial_sum(ii,jj)%boa_lw_up      = spatial_sum(ii,jj)%boa_lw_up + var2d(i,j)%boa_lw_up
             spatial_sum(ii,jj)%boa_lw_dn      = spatial_sum(ii,jj)%boa_lw_dn + var2d(i,j)%boa_lw_dn
             spatial_sum(ii,jj)%boa_sw_up      = spatial_sum(ii,jj)%boa_sw_up + var2d(i,j)%boa_sw_up
             spatial_sum(ii,jj)%boa_sw_dn      = spatial_sum(ii,jj)%boa_sw_dn + var2d(i,j)%boa_sw_dn

          endif
        enddo
         
        endif within_grid
      enddo


   enddo ; enddo


  endif

 case('finish') !---------------------------------------------------------------------------------------

#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', spatial_cnt  )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%h2o_col   )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_cloud_liq )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_rain_liq )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_cloud_ice )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_rain_ice )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_lw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_sw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_sw_dn )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_lw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_sw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_sw_dn )

  MPI_MASTERPROC: if(masterproc) then
#endif

  do jj = 1, jmax ; do ii = 1, imax

    if( spatial_cnt(ii,jj) > 0. ) then

     spatial_mean(ii,jj)%rain_rate      = spatial_sum(ii,jj)%rain_rate / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%h2o_col        = spatial_sum(ii,jj)%h2o_col   / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_cloud_liq = spatial_sum(ii,jj)%path_cloud_liq / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_rain_liq  = spatial_sum(ii,jj)%path_rain_liq  / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_cloud_ice = spatial_sum(ii,jj)%path_cloud_ice / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_rain_ice  = spatial_sum(ii,jj)%path_rain_ice  / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_lw_up      = spatial_sum(ii,jj)%toa_lw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_lw_dn      = spatial_sum(ii,jj)%toa_lw_dn / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_sw_up      = spatial_sum(ii,jj)%toa_sw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_sw_dn      = spatial_sum(ii,jj)%toa_sw_dn / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_lw_up      = spatial_sum(ii,jj)%boa_lw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_lw_dn      = spatial_sum(ii,jj)%boa_lw_dn / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_sw_up      = spatial_sum(ii,jj)%boa_sw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_sw_dn      = spatial_sum(ii,jj)%boa_sw_dn / spatial_cnt(ii,jj)
    
    endif

  enddo ; enddo


!
! write out lat-lon mean value
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM2D'//trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=imax*jmax*4)
 print("(a)"),'MSG zonal_mean2d: see output ->',trim(sdsu_io_file)

 irec = 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%rain_rate      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%h2o_col        , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_cloud_liq , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_rain_liq  , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_cloud_ice , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_rain_ice  , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_lw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_lw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_sw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_sw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_lw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_lw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_sw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_sw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1

 close(io)


!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM2D'//trim(output_suffix)//'.ctl'
 print("(a)"),'MSG spatial_mean2d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'spatial_mean_CRM2D'//trim(output_suffix)//'.grads'
 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   spatial mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    ',imax,' linear ',minlon,' ',res_latlon
 write(io,*) 'ydef    ',jmax,' linear ',minlat,' ',res_latlon
 write(io,*) 'zdef    1 linear 0. 0.1'
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars ',irec-1
 write(io,*) 'rain    0  0  rain rate [mm/hr]'
 write(io,*) 'h2ocol  0  0  column water vapro [kg/m2]'
 write(io,*) 'clwp    0  0  cloud liquid water path [kg/m2]'
 write(io,*) 'lrwp    0  0  liquid rain water path [kg/m2]'
 write(io,*) 'ciwp    0  0  cloud ice water path [kg/m2]'
 write(io,*) 'irwp    0  0  ice rain water path [kg/m2]'
 write(io,*) 'lwut   0 0  TOA LW upwelling   flux [W/m2]'
 write(io,*) 'lwdt   0 0  TOA LW downwelling flux [W/m2]'
 write(io,*) 'swut   0 0  TOA SW upwelling flux   [W/m2]'
 write(io,*) 'swdt   0 0  TOA SW downwelling flux [W/m2]'
 write(io,*) 'lwus   0 0  Surface LW upwelling   flux [W/m2]'
 write(io,*) 'lwds   0 0  Surface LW downwelling flux [W/m2]'
 write(io,*) 'swus   0 0  Surface SW upwelling   flux [W/m2]'
 write(io,*) 'swds   0 0  Surface SW downwelling flux [W/m2]'
 write(io,*) 'endvars '
 close(io)


#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG spatial_mean_2d: There is no such proc')
 end select


 return
 end subroutine spatial_mean_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine spatial_mean_daily_2d(proc,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive daily spatial mean map of 2D parameters.
!           
! History:
! 06/2017  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 character*(*) :: proc ! process
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 2d array

 integer :: bnd(2)  !upper memory bound for 2D array

!
! statisticsl parameters
!
 type( surface_parameter ), allocatable, dimension(:,:),save :: spatial_sum, spatial_mean
 real(sdsu_fps), allocatable, dimension(:,:),save :: spatial_cnt
 real(sdsu_fps), allocatable, dimension(:),save :: spatial_lat, spatial_lon

 integer,save :: ibin_max
 integer :: i,j, ibin, ierr, irec
 integer,parameter :: io = 149

 integer,save :: imax, jmax
 integer :: ii, jj

 character :: yyyy*4, mm*2, dd*2

 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

   if(masterproc) print*, 'MSG spatial_mean_daily_2d: init'

   imax = int( ( maxlon - minlon ) / res_latlon ) + 1
   jmax = int( ( maxlat - minlat ) / res_latlon ) + 1

   if( imax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_daily_2d: maxlon is smaller than minlon')
   if( jmax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_daily_2d: maxlat is smaller than minlat')

  
   if( allocated(spatial_sum) ) deallocate(spatial_sum, spatial_cnt, spatial_mean, spatial_lat, spatial_lon)
 
   allocate( spatial_sum(imax,jmax)  , &
             spatial_cnt(imax,jmax)  , &
             spatial_mean(imax,jmax) , &
             spatial_lat(jmax)  , &
             spatial_lon(imax)  , &
             stat=ierr )

   ! initialized
   do j = 1, jmax
      spatial_lat(j) = minlat + res_latlon * real(j-1) ![deg]
   enddo
   do i = 1, imax
      spatial_lon(i) = minlon + res_latlon * real(i-1) ![deg]
   enddo

  !
  ! zero out other parameters
  !
  spatial_cnt = 0.

  spatial_sum%rain_rate      = 0. ; spatial_sum%h2o_col       = 0.
  spatial_sum%path_cloud_liq = 0. ; spatial_sum%path_rain_liq = 0.
  spatial_sum%path_cloud_ice = 0. ; spatial_sum%path_rain_ice = 0.
  spatial_sum%toa_lw_up      = 0. ; spatial_sum%toa_lw_dn     = 0.
  spatial_sum%toa_sw_up      = 0. ; spatial_sum%toa_sw_dn     = 0.
  spatial_sum%boa_lw_up      = 0. ; spatial_sum%boa_lw_dn     = 0.
  spatial_sum%boa_sw_up      = 0. ; spatial_sum%boa_sw_dn     = 0.

  spatial_mean%rain_rate      = undefined ; spatial_mean%h2o_col       = undefined
  spatial_mean%path_cloud_liq = undefined ; spatial_mean%path_rain_liq = undefined
  spatial_mean%path_cloud_ice = undefined ; spatial_mean%path_rain_ice = undefined
  spatial_mean%toa_lw_up      = undefined ; spatial_mean%toa_lw_dn     = undefined
  spatial_mean%toa_sw_up      = undefined ; spatial_mean%toa_sw_dn     = undefined
  spatial_mean%boa_lw_up      = undefined ; spatial_mean%boa_lw_dn     = undefined
  spatial_mean%boa_sw_up      = undefined ; spatial_mean%boa_sw_dn     = undefined


 case('add') !-------------------------------------------------------------------------------------

  if(masterproc) print*, 'MSG spatial_mean_daily_2d: add'

  if(present(var2d) ) then

   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var2d)

   do j = 1,bnd(2) ; do i = 1,bnd(1)

      !
      ! loop for latlon_grid domain. 
      !
      do jj = 1, jmax
        within_grid: if( var2d(i,j)%lat >  (spatial_lat(jj)-0.5*res_latlon) .and. &
                         var2d(i,j)%lat <= (spatial_lat(jj)+0.5*res_latlon) ) then

!        if(masterproc) print*, 'MSG spatial_mean_daily_2d: found lat ', (spatial_lat(jj)-0.5*res_latlon), &
!         var2d(i,j)%lat, (spatial_lat(jj)+0.5*res_latlon)

        do ii = 1, imax
          if( var2d(i,j)%lon >  (spatial_lon(ii)-0.5*res_latlon) .and. &
              var2d(i,j)%lon <= (spatial_lon(ii)+0.5*res_latlon) ) then

!             if(masterproc) print*, 'MSG spatial_mean_daily_2d: found lon ',(spatial_lon(ii)-0.5*res_latlon),&
!                var2d(i,j)%lon, (spatial_lon(ii)+0.5*res_latlon)
 
             spatial_cnt(ii,jj) = spatial_cnt(ii,jj) + 1.e0  !update total count

             spatial_sum(ii,jj)%rain_rate      = spatial_sum(ii,jj)%rain_rate + var2d(i,j)%rain_rate
             spatial_sum(ii,jj)%h2o_col        = spatial_sum(ii,jj)%h2o_col   + var2d(i,j)%h2o_col
             spatial_sum(ii,jj)%path_cloud_liq = spatial_sum(ii,jj)%path_cloud_liq + var2d(i,j)%path_cloud_liq
             spatial_sum(ii,jj)%path_rain_liq  = spatial_sum(ii,jj)%path_rain_liq  + var2d(i,j)%path_rain_liq
             spatial_sum(ii,jj)%path_cloud_ice = spatial_sum(ii,jj)%path_cloud_ice + var2d(i,j)%path_cloud_ice
             spatial_sum(ii,jj)%path_rain_ice  = spatial_sum(ii,jj)%path_rain_ice  + var2d(i,j)%path_rain_ice
             spatial_sum(ii,jj)%toa_lw_up      = spatial_sum(ii,jj)%toa_lw_up + var2d(i,j)%toa_lw_up
             spatial_sum(ii,jj)%toa_lw_dn      = spatial_sum(ii,jj)%toa_lw_dn + var2d(i,j)%toa_lw_dn
             spatial_sum(ii,jj)%toa_sw_up      = spatial_sum(ii,jj)%toa_sw_up + var2d(i,j)%toa_sw_up
             spatial_sum(ii,jj)%toa_sw_dn      = spatial_sum(ii,jj)%toa_sw_dn + var2d(i,j)%toa_sw_dn
             spatial_sum(ii,jj)%boa_lw_up      = spatial_sum(ii,jj)%boa_lw_up + var2d(i,j)%boa_lw_up
             spatial_sum(ii,jj)%boa_lw_dn      = spatial_sum(ii,jj)%boa_lw_dn + var2d(i,j)%boa_lw_dn
             spatial_sum(ii,jj)%boa_sw_up      = spatial_sum(ii,jj)%boa_sw_up + var2d(i,j)%boa_sw_up
             spatial_sum(ii,jj)%boa_sw_dn      = spatial_sum(ii,jj)%boa_sw_dn + var2d(i,j)%boa_sw_dn

          endif
        enddo
         
        endif within_grid
      enddo


   enddo ; enddo


  endif

 case('finish') !---------------------------------------------------------------------------------------

#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', spatial_cnt  )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%h2o_col   )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_cloud_liq )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_rain_liq )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_cloud_ice )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_rain_ice )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_lw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_sw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_sw_dn )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_lw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_sw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_sw_dn )

  MPI_MASTERPROC: if(masterproc) then
#endif

  do jj = 1, jmax ; do ii = 1, imax

    if( spatial_cnt(ii,jj) > 0. ) then

     spatial_mean(ii,jj)%rain_rate      = spatial_sum(ii,jj)%rain_rate / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%h2o_col        = spatial_sum(ii,jj)%h2o_col   / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_cloud_liq = spatial_sum(ii,jj)%path_cloud_liq / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_rain_liq  = spatial_sum(ii,jj)%path_rain_liq  / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_cloud_ice = spatial_sum(ii,jj)%path_cloud_ice / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_rain_ice  = spatial_sum(ii,jj)%path_rain_ice  / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_lw_up      = spatial_sum(ii,jj)%toa_lw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_lw_dn      = spatial_sum(ii,jj)%toa_lw_dn / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_sw_up      = spatial_sum(ii,jj)%toa_sw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_sw_dn      = spatial_sum(ii,jj)%toa_sw_dn / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_lw_up      = spatial_sum(ii,jj)%boa_lw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_lw_dn      = spatial_sum(ii,jj)%boa_lw_dn / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_sw_up      = spatial_sum(ii,jj)%boa_sw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_sw_dn      = spatial_sum(ii,jj)%boa_sw_dn / spatial_cnt(ii,jj)
    
    endif

  enddo ; enddo


!
! write out cloud frequency diurnal profile
!
 write(yyyy,"(I4.4)") int( sdsu_yyyy )
 write(mm  ,"(I2.2)") int( sdsu_mm   )
 write(dd  ,"(I2.2)") int( sdsu_dd   )

!
! write out cloud frequency zonal profile
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM2D_'//yyyy//'_'//mm//'_'//dd&
                 //trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=imax*jmax*4)
 print("(a)"),'MSG spatial_mean_daily_2d: see output ->',trim(sdsu_io_file)

 irec = 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%rain_rate      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%h2o_col        , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_cloud_liq , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_rain_liq  , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_cloud_ice , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_rain_ice  , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_lw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_lw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_sw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_sw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_lw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_lw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_sw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_sw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1

 close(io)


!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM2D_'//yyyy//'_'//mm//'_'//dd&
                //trim(output_suffix)//'.ctl'
 print("(a)"),'MSG spatial_mean_daily_2d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'spatial_mean_CRM2D_'//yyyy//'_'//mm//'_'//dd&
                 //trim(output_suffix)//'.grads'
 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   spatial mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    ',imax,' linear ',minlon,' ',res_latlon
 write(io,*) 'ydef    ',jmax,' linear ',minlat,' ',res_latlon
 write(io,*) 'zdef    1 linear 0. 0.1'
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars ',irec-1
 write(io,*) 'rain    0  0  rain rate [mm/hr]'
 write(io,*) 'h2ocol  0  0  column water vapro [kg/m2]'
 write(io,*) 'clwp    0  0  cloud liquid water path [kg/m2]'
 write(io,*) 'lrwp    0  0  liquid rain water path [kg/m2]'
 write(io,*) 'ciwp    0  0  cloud ice water path [kg/m2]'
 write(io,*) 'irwp    0  0  ice rain water path [kg/m2]'
 write(io,*) 'lwut   0 0  TOA LW upwelling   flux [W/m2]'
 write(io,*) 'lwdt   0 0  TOA LW downwelling flux [W/m2]'
 write(io,*) 'swut   0 0  TOA SW upwelling flux   [W/m2]'
 write(io,*) 'swdt   0 0  TOA SW downwelling flux [W/m2]'
 write(io,*) 'lwus   0 0  Surface LW upwelling   flux [W/m2]'
 write(io,*) 'lwds   0 0  Surface LW downwelling flux [W/m2]'
 write(io,*) 'swus   0 0  Surface SW upwelling   flux [W/m2]'
 write(io,*) 'swds   0 0  Surface SW downwelling flux [W/m2]'
 write(io,*) 'endvars '
 close(io)

#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG spatial_mean_daily_2d: There is no such proc')
 end select

 return
 end subroutine spatial_mean_daily_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine spatial_mean_hourly_2d(proc,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive hourly spatial mean map of 2D parameters.
!           
! History:
! 06/2017  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 character*(*) :: proc ! process
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 2d array

 integer :: bnd(2)  !upper memory bound for 2D array

!
! statisticsl parameters
!
 type( surface_parameter ), allocatable, dimension(:,:),save :: spatial_sum, spatial_mean
 real(sdsu_fps), allocatable, dimension(:,:),save :: spatial_cnt
 real(sdsu_fps), allocatable, dimension(:),save :: spatial_lat, spatial_lon

 integer,save :: ibin_max
 integer :: i,j, ibin, ierr, irec
 integer,parameter :: io = 149

 integer,save :: imax, jmax
 integer :: ii, jj

 character :: yyyy*4, mm*2, dd*2, hh*2

 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

   if(masterproc) print*, 'MSG spatial_mean_hourly_2d: init'

   imax = int( ( maxlon - minlon ) / res_latlon ) + 1
   jmax = int( ( maxlat - minlat ) / res_latlon ) + 1

   if( imax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_hourly_2d: maxlon is smaller than minlon')
   if( jmax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_hourly_2d: maxlat is smaller than minlat')

  
   if( allocated(spatial_sum) ) deallocate(spatial_sum, spatial_cnt, spatial_mean, spatial_lat, spatial_lon)
 
   allocate( spatial_sum(imax,jmax)  , &
             spatial_cnt(imax,jmax)  , &
             spatial_mean(imax,jmax) , &
             spatial_lat(jmax)  , &
             spatial_lon(imax)  , &
             stat=ierr )

   ! initialized
   do j = 1, jmax
      spatial_lat(j) = minlat + res_latlon * real(j-1) ![deg]
   enddo
   do i = 1, imax
      spatial_lon(i) = minlon + res_latlon * real(i-1) ![deg]
   enddo

  !
  ! zero out other parameters
  !
  spatial_cnt = 0.

  spatial_sum%rain_rate      = 0. ; spatial_sum%h2o_col       = 0.
  spatial_sum%path_cloud_liq = 0. ; spatial_sum%path_rain_liq = 0.
  spatial_sum%path_cloud_ice = 0. ; spatial_sum%path_rain_ice = 0.
  spatial_sum%toa_lw_up      = 0. ; spatial_sum%toa_lw_dn     = 0.
  spatial_sum%toa_sw_up      = 0. ; spatial_sum%toa_sw_dn     = 0.
  spatial_sum%boa_lw_up      = 0. ; spatial_sum%boa_lw_dn     = 0.
  spatial_sum%boa_sw_up      = 0. ; spatial_sum%boa_sw_dn     = 0.

  spatial_mean%rain_rate      = undefined ; spatial_mean%h2o_col       = undefined
  spatial_mean%path_cloud_liq = undefined ; spatial_mean%path_rain_liq = undefined
  spatial_mean%path_cloud_ice = undefined ; spatial_mean%path_rain_ice = undefined
  spatial_mean%toa_lw_up      = undefined ; spatial_mean%toa_lw_dn     = undefined
  spatial_mean%toa_sw_up      = undefined ; spatial_mean%toa_sw_dn     = undefined
  spatial_mean%boa_lw_up      = undefined ; spatial_mean%boa_lw_dn     = undefined
  spatial_mean%boa_sw_up      = undefined ; spatial_mean%boa_sw_dn     = undefined


 case('add') !-------------------------------------------------------------------------------------

  if(masterproc) print*, 'MSG spatial_mean_hourly_2d: add'

  if(present(var2d) ) then

   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var2d)

   do j = 1,bnd(2) ; do i = 1,bnd(1)

      !
      ! loop for latlon_grid domain. 
      !
      do jj = 1, jmax
        within_grid: if( var2d(i,j)%lat >  (spatial_lat(jj)-0.5*res_latlon) .and. &
                         var2d(i,j)%lat <= (spatial_lat(jj)+0.5*res_latlon) ) then

!        if(masterproc) print*, 'MSG spatial_mean_hourly_2d: found lat ', (spatial_lat(jj)-0.5*res_latlon), &
!         var2d(i,j)%lat, (spatial_lat(jj)+0.5*res_latlon)

        do ii = 1, imax
          if( var2d(i,j)%lon >  (spatial_lon(ii)-0.5*res_latlon) .and. &
              var2d(i,j)%lon <= (spatial_lon(ii)+0.5*res_latlon) ) then

!             if(masterproc) print*, 'MSG spatial_mean_hourly_2d: found lon ',(spatial_lon(ii)-0.5*res_latlon),&
!                var2d(i,j)%lon, (spatial_lon(ii)+0.5*res_latlon)
 
             spatial_cnt(ii,jj) = spatial_cnt(ii,jj) + 1.e0  !update total count

             spatial_sum(ii,jj)%rain_rate      = spatial_sum(ii,jj)%rain_rate + var2d(i,j)%rain_rate
             spatial_sum(ii,jj)%h2o_col        = spatial_sum(ii,jj)%h2o_col   + var2d(i,j)%h2o_col
             spatial_sum(ii,jj)%path_cloud_liq = spatial_sum(ii,jj)%path_cloud_liq + var2d(i,j)%path_cloud_liq
             spatial_sum(ii,jj)%path_rain_liq  = spatial_sum(ii,jj)%path_rain_liq  + var2d(i,j)%path_rain_liq
             spatial_sum(ii,jj)%path_cloud_ice = spatial_sum(ii,jj)%path_cloud_ice + var2d(i,j)%path_cloud_ice
             spatial_sum(ii,jj)%path_rain_ice  = spatial_sum(ii,jj)%path_rain_ice  + var2d(i,j)%path_rain_ice
             spatial_sum(ii,jj)%toa_lw_up      = spatial_sum(ii,jj)%toa_lw_up + var2d(i,j)%toa_lw_up
             spatial_sum(ii,jj)%toa_lw_dn      = spatial_sum(ii,jj)%toa_lw_dn + var2d(i,j)%toa_lw_dn
             spatial_sum(ii,jj)%toa_sw_up      = spatial_sum(ii,jj)%toa_sw_up + var2d(i,j)%toa_sw_up
             spatial_sum(ii,jj)%toa_sw_dn      = spatial_sum(ii,jj)%toa_sw_dn + var2d(i,j)%toa_sw_dn
             spatial_sum(ii,jj)%boa_lw_up      = spatial_sum(ii,jj)%boa_lw_up + var2d(i,j)%boa_lw_up
             spatial_sum(ii,jj)%boa_lw_dn      = spatial_sum(ii,jj)%boa_lw_dn + var2d(i,j)%boa_lw_dn
             spatial_sum(ii,jj)%boa_sw_up      = spatial_sum(ii,jj)%boa_sw_up + var2d(i,j)%boa_sw_up
             spatial_sum(ii,jj)%boa_sw_dn      = spatial_sum(ii,jj)%boa_sw_dn + var2d(i,j)%boa_sw_dn

          endif
        enddo
         
        endif within_grid
      enddo


   enddo ; enddo


  endif

 case('finish') !---------------------------------------------------------------------------------------

#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', spatial_cnt  )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%h2o_col   )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_cloud_liq )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_rain_liq )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_cloud_ice )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%path_rain_ice )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_lw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_sw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%toa_sw_dn )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_lw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_sw_up )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%boa_sw_dn )

  MPI_MASTERPROC: if(masterproc) then
#endif

  do jj = 1, jmax ; do ii = 1, imax

    if( spatial_cnt(ii,jj) > 0. ) then

     spatial_mean(ii,jj)%rain_rate      = spatial_sum(ii,jj)%rain_rate / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%h2o_col        = spatial_sum(ii,jj)%h2o_col   / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_cloud_liq = spatial_sum(ii,jj)%path_cloud_liq / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_rain_liq  = spatial_sum(ii,jj)%path_rain_liq  / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_cloud_ice = spatial_sum(ii,jj)%path_cloud_ice / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%path_rain_ice  = spatial_sum(ii,jj)%path_rain_ice  / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_lw_up      = spatial_sum(ii,jj)%toa_lw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_lw_dn      = spatial_sum(ii,jj)%toa_lw_dn / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_sw_up      = spatial_sum(ii,jj)%toa_sw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%toa_sw_dn      = spatial_sum(ii,jj)%toa_sw_dn / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_lw_up      = spatial_sum(ii,jj)%boa_lw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_lw_dn      = spatial_sum(ii,jj)%boa_lw_dn / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_sw_up      = spatial_sum(ii,jj)%boa_sw_up / spatial_cnt(ii,jj)
     spatial_mean(ii,jj)%boa_sw_dn      = spatial_sum(ii,jj)%boa_sw_dn / spatial_cnt(ii,jj)
    
    endif

  enddo ; enddo


!
! write out cloud frequency diurnal profile
!
 write(yyyy,"(I4.4)") int( sdsu_yyyy )
 write(mm  ,"(I2.2)") int( sdsu_mm   )
 write(dd  ,"(I2.2)") int( sdsu_dd   )
 write(hh  ,"(I2.2)") int( sdsu_hh   )

!
! write out binary data
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM2D_'//yyyy//'_'//mm//'_'//dd//'_'//hh&
                 //trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=imax*jmax*4)
 print("(a)"),'MSG spatial_mean_hourly_2d: see output ->',trim(sdsu_io_file)

 irec = 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%rain_rate      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%h2o_col        , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_cloud_liq , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_rain_liq  , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_cloud_ice , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%path_rain_ice  , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_lw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_lw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_sw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%toa_sw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_lw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_lw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_sw_up      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1
 write(io,rec=irec) ((spatial_mean(ii,jj)%boa_sw_dn      , ii=1,imax ),jj=1,jmax) ; irec = irec + 1

 close(io)


!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM2D_'//yyyy//'_'//mm//'_'//dd//'_'//hh &
                //trim(output_suffix)//'.ctl'
 print("(a)"),'MSG spatial_mean_hourly_2d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'spatial_mean_CRM2D_'//yyyy//'_'//mm//'_'//dd//'_'//hh&
                 //trim(output_suffix)//'.grads'
 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   spatial mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    ',imax,' linear ',minlon,' ',res_latlon
 write(io,*) 'ydef    ',jmax,' linear ',minlat,' ',res_latlon
 write(io,*) 'zdef    1 linear 0. 0.1'
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars ',irec-1
 write(io,*) 'rain    0  0  rain rate [mm/hr]'
 write(io,*) 'h2ocol  0  0  column water vapro [kg/m2]'
 write(io,*) 'clwp    0  0  cloud liquid water path [kg/m2]'
 write(io,*) 'lrwp    0  0  liquid rain water path [kg/m2]'
 write(io,*) 'ciwp    0  0  cloud ice water path [kg/m2]'
 write(io,*) 'irwp    0  0  ice rain water path [kg/m2]'
 write(io,*) 'lwut   0 0  TOA LW upwelling   flux [W/m2]'
 write(io,*) 'lwdt   0 0  TOA LW downwelling flux [W/m2]'
 write(io,*) 'swut   0 0  TOA SW upwelling flux   [W/m2]'
 write(io,*) 'swdt   0 0  TOA SW downwelling flux [W/m2]'
 write(io,*) 'lwus   0 0  Surface LW upwelling   flux [W/m2]'
 write(io,*) 'lwds   0 0  Surface LW downwelling flux [W/m2]'
 write(io,*) 'swus   0 0  Surface SW upwelling   flux [W/m2]'
 write(io,*) 'swds   0 0  Surface SW downwelling flux [W/m2]'
 write(io,*) 'endvars '
 close(io)

#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG spatial_mean_hourly_2d: There is no such proc')
 end select

 return
 end subroutine spatial_mean_hourly_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine spatial_mean_3d(proc,var3d,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive spatial mean map of 2D parameters.
!           
! History:
! 11/2016  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 character*(*) :: proc ! process
 type ( atmos_parameter ), dimension(:,:,:), optional :: var3d ! input surface parameter assumed 3d array
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 3d array

 integer :: bnd(3)  !upper memory bound for 2D array

!
! statisticsl parameters
!
 type( atmos_parameter ), allocatable, dimension(:,:,:),save :: spatial_sum, spatial_mean
 real(sdsu_fps), allocatable, dimension(:,:,:),save :: spatial_cnt
 real(sdsu_fps), allocatable, dimension(:),save :: spatial_lat, spatial_lon

 integer,save :: ibin_max
 integer :: i,j,k, ibin, ierr, irec
 integer,parameter :: io = 149


 integer,parameter :: play = 20, plev = 21
 real(sdsu_fps),parameter :: pressure_lay(play) =    (/1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,500.,&
                                                     450.,400.,350.,300.,250.,200.,150.,100.,50./)
 real(sdsu_fps),parameter :: pressure_lev(plev) = (/1025., 975.,925.,875.,825.,775.,725.,675.,625.,575.,525.,475.,&
                                                        425.,375.,325.,275.,225.,175.,125.,75.,25./)

 integer,save :: imax, jmax, kmax
 integer :: ii, jj, kk

 real(sdsu_fps) :: pr


 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

   if(masterproc) print*, 'MSG spatial_mean_3d: init'

   imax = int( ( maxlon - minlon ) / res_latlon ) + 1
   jmax = int( ( maxlat - minlat ) / res_latlon ) + 1
   kmax = play

   if( imax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_3d: maxlon is smaller than minlon')
   if( jmax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_3d: maxlat is smaller than minlat')

   if( allocated(spatial_sum) ) deallocate(spatial_sum, spatial_cnt, spatial_mean, spatial_lat, spatial_lon)

   allocate( spatial_sum(imax,jmax,kmax)  , &
             spatial_cnt(imax,jmax,kmax)  , &
             spatial_mean(imax,jmax,kmax) , &
             spatial_lat(jmax)  , &
             spatial_lon(imax)  , &
             stat=ierr )

   ! initialized
   do j = 1, jmax
      spatial_lat(j) = minlat + res_latlon * real(j-1) ![deg]
   enddo
   do i = 1, imax
      spatial_lon(i) = minlon + res_latlon * real(i-1) ![deg]
   enddo

  !
  ! zero out other parameters
  !
  spatial_cnt = 0.

  spatial_sum%rain_rate      = 0. ; spatial_sum%t_air       = 0.
  spatial_sum%sh             = 0. ; spatial_sum%u           = 0.
  spatial_sum%v              = 0. ; spatial_sum%w           = 0.
  spatial_sum%lw_heat        = 0. ; spatial_sum%sw_heat     = 0.
  spatial_sum%tot_cloud      = 0.

  spatial_mean%rain_rate      = undefined ; spatial_mean%t_air       = undefined
  spatial_mean%sh             = undefined ; spatial_mean%u           = undefined
  spatial_mean%v              = undefined ; spatial_mean%w           = undefined
  spatial_mean%lw_heat        = undefined ; spatial_mean%sw_heat     = undefined
  spatial_mean%tot_cloud      = undefined 

 case('add') !-------------------------------------------------------------------------------------

  if(masterproc) print*, 'MSG spatial_mean_3d: add'

  if(present(var3d) ) then

   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var3d)

   do j = 1,bnd(2) ; do i = 1,bnd(1)

      !
      ! loop for latlon_grid domain. 
      !
      do jj = 1, jmax
        within_grid: if( var2d(i,j)%lat >  (spatial_lat(jj)-0.5*res_latlon) .and. &
                         var2d(i,j)%lat <= (spatial_lat(jj)+0.5*res_latlon) ) then

!        if(masterproc) print*, 'MSG spatial_mean_2d: found lat ', (spatial_lat(jj)-0.5*res_latlon), &
!         var2d(i,j)%lat, (spatial_lat(jj)+0.5*res_latlon)

        do ii = 1, imax
          if( var2d(i,j)%lon >  (spatial_lon(ii)-0.5*res_latlon) .and. &
              var2d(i,j)%lon <= (spatial_lon(ii)+0.5*res_latlon) ) then

!             if(masterproc) print*, 'MSG spatial_mean_2d: found lon ',(spatial_lon(ii)-0.5*res_latlon),&
!                var2d(i,j)%lon, (spatial_lon(ii)+0.5*res_latlon)

           do kk = 1, kmax

               do k = 1, bnd(3)
                  if( var3d(i,j,k)%press <= pressure_lev(kk) .and. &
                      var3d(i,j,k)%press >  pressure_lev(kk+1)  ) then

                    spatial_cnt(ii,jj,kk) = spatial_cnt(ii,jj,kk) + 1.e0  !update total count
                    spatial_sum(ii,jj,kk)%rain_rate  = spatial_sum(ii,jj,kk)%rain_rate + var3d(i,j,k)%rain_rate
                    spatial_sum(ii,jj,kk)%t_air      = spatial_sum(ii,jj,kk)%t_air     + var3d(i,j,k)%t_air
                    spatial_sum(ii,jj,kk)%sh         = spatial_sum(ii,jj,kk)%sh        + var3d(i,j,k)%sh
                    spatial_sum(ii,jj,kk)%u          = spatial_sum(ii,jj,kk)%u         + var3d(i,j,k)%u
                    spatial_sum(ii,jj,kk)%v          = spatial_sum(ii,jj,kk)%v         + var3d(i,j,k)%v
                    spatial_sum(ii,jj,kk)%w          = spatial_sum(ii,jj,kk)%w         + var3d(i,j,k)%w
                    spatial_sum(ii,jj,kk)%lw_heat    = spatial_sum(ii,jj,kk)%lw_heat   + var3d(i,j,k)%lw_heat
                    spatial_sum(ii,jj,kk)%sw_heat    = spatial_sum(ii,jj,kk)%sw_heat   + var3d(i,j,k)%sw_heat
                    spatial_sum(ii,jj,kk)%tot_cloud  = spatial_sum(ii,jj,kk)%tot_cloud + var3d(i,j,k)%tot_cloud

                  endif
               enddo

           enddo
           endif
        enddo

        endif within_grid
      enddo


   enddo ; enddo


  endif

 case('finish') !---------------------------------------------------------------------------------------

#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', spatial_cnt  )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%t_air )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%sh   )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%u )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%v )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%w )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%lw_heat )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%sw_heat )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%tot_cloud )

  MPI_MASTERPROC: if(masterproc) then
#endif

  do kk = 1, kmax ; do jj = 1, jmax ; do ii = 1, imax 

    if( spatial_cnt(ii,jj,kk) > 0. ) then

     spatial_mean(ii,jj,kk)%rain_rate = spatial_sum(ii,jj,kk)%rain_rate / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%t_air     = spatial_sum(ii,jj,kk)%t_air / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%sh        = spatial_sum(ii,jj,kk)%sh   / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%u         = spatial_sum(ii,jj,kk)%u    / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%v         = spatial_sum(ii,jj,kk)%v    / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%w         = spatial_sum(ii,jj,kk)%w    / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%lw_heat   = spatial_sum(ii,jj,kk)%lw_heat / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%sw_heat   = spatial_sum(ii,jj,kk)%sw_heat / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%tot_cloud = spatial_sum(ii,jj,kk)%tot_cloud / spatial_cnt(ii,jj,kk)

    endif

  enddo ; enddo ; enddo


!
! write out cloud frequency zonal profile
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM3D'//trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=imax*jmax*kmax*4)
 print("(a)"),'MSG spatial_mean3d: see output ->',trim(sdsu_io_file)

 irec = 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%rain_rate , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%t_air     , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%sh        , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%u         , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%v         , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%w         , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%lw_heat   , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%sw_heat   , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%tot_cloud , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1

 close(io)


!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM3D'//trim(output_suffix)//'.ctl'
 print("(a)"),'MSG spatial_mean3d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'spatial_mean_CRM3D'//trim(output_suffix)//'.grads'
 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   spatial mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    ',imax,' linear ',minlon,' ',res_latlon
 write(io,*) 'ydef    ',jmax,' linear ',minlat,' ',res_latlon
 write(io,*) 'zdef    ',play,' levels ', pressure_lay(1:play)
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars ',irec-1
 write(io,*) 'rain   ',play,'  0  rain rate [mm/hr]'
 write(io,*) 'tair   ',play,'  0  air temperature [degK]'
 write(io,*) 'sh   ',play ,'  0  specific humidity [g/g]'
 write(io,*) 'u    ',play,'  0  u wind [m/s]'
 write(io,*) 'v    ',play,'  0  v wind [m/s]'
 write(io,*) 'w    ',play,'  0  w wind [m/s]'
 write(io,*) 'lwheat',play,'  0  LW heating rate [K/day]'
 write(io,*) 'swheat',play,'  0  SW heating rate [K/day]'
 write(io,*) 'totcloud',play,'  0  total cloud + ice [g/m3]'

 write(io,*) 'endvars '
 close(io)


#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG spatial_mean_3d: There is no such proc')
 end select


 return
 end subroutine spatial_mean_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 


 subroutine spatial_mean_daily_3d(proc,var3d,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive spatial mean daily map of 3D parameters.
!           
! History:
! 11/2016  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 character*(*) :: proc ! process
 type ( atmos_parameter ), dimension(:,:,:), optional :: var3d ! input surface parameter assumed 3d array
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 3d array

 integer :: bnd(3)  !upper memory bound for 2D array

!
! statisticsl parameters
!
 type( atmos_parameter ), allocatable, dimension(:,:,:),save :: spatial_sum, spatial_mean
 real(sdsu_fps), allocatable, dimension(:,:,:),save :: spatial_cnt
 real(sdsu_fps), allocatable, dimension(:),save :: spatial_lat, spatial_lon

 integer,save :: ibin_max
 integer :: i,j,k, ibin, ierr, irec
 integer,parameter :: io = 149


 integer,parameter :: play = 20, plev = 21
 real(sdsu_fps),parameter :: pressure_lay(play) =    (/1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,500.,&
                                                     450.,400.,350.,300.,250.,200.,150.,100.,50./)
 real(sdsu_fps),parameter :: pressure_lev(plev) = (/1025., 975.,925.,875.,825.,775.,725.,675.,625.,575.,525.,475.,&
                                                        425.,375.,325.,275.,225.,175.,125.,75.,25./)

 integer,save :: imax, jmax, kmax
 integer :: ii, jj, kk

 real(sdsu_fps) :: pr

 character :: yyyy*4, mm*2, dd*2


 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

   if(masterproc) print*, 'MSG spatial_mean_daily_3d: init'

   imax = int( ( maxlon - minlon ) / res_latlon ) + 1
   jmax = int( ( maxlat - minlat ) / res_latlon ) + 1
   kmax = play

   if( imax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_daily_3d: maxlon is smaller than minlon')
   if( jmax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_daily_3d: maxlat is smaller than minlat')

   if( allocated(spatial_sum) ) deallocate(spatial_sum, spatial_cnt, spatial_mean, spatial_lat, spatial_lon)

   allocate( spatial_sum(imax,jmax,kmax)  , &
             spatial_cnt(imax,jmax,kmax)  , &
             spatial_mean(imax,jmax,kmax) , &
             spatial_lat(jmax)  , &
             spatial_lon(imax)  , &
             stat=ierr )

   ! initialized
   do j = 1, jmax
      spatial_lat(j) = minlat + res_latlon * real(j-1) ![deg]
   enddo
   do i = 1, imax
      spatial_lon(i) = minlon + res_latlon * real(i-1) ![deg]
   enddo

  !
  ! zero out other parameters
  !
  spatial_cnt = 0.

  spatial_sum%rain_rate      = 0. ; spatial_sum%t_air       = 0.
  spatial_sum%sh             = 0. ; spatial_sum%u           = 0.
  spatial_sum%v              = 0. ; spatial_sum%w           = 0.
  spatial_sum%lw_heat        = 0. ; spatial_sum%sw_heat     = 0.
  spatial_sum%tot_cloud      = 0.

  spatial_mean%rain_rate      = undefined ; spatial_mean%t_air       = undefined
  spatial_mean%sh             = undefined ; spatial_mean%u           = undefined
  spatial_mean%v              = undefined ; spatial_mean%w           = undefined
  spatial_mean%lw_heat        = undefined ; spatial_mean%sw_heat     = undefined
  spatial_mean%tot_cloud      = undefined

 case('add') !-------------------------------------------------------------------------------------

  if(masterproc) print*, 'MSG spatial_mean_daily_3d: add'

  if(present(var3d) ) then

   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var3d)

   do j = 1,bnd(2) ; do i = 1,bnd(1)

      !
      ! loop for latlon_grid domain. 
      !
      do jj = 1, jmax
        within_grid: if( var2d(i,j)%lat >  (spatial_lat(jj)-0.5*res_latlon) .and. &
                         var2d(i,j)%lat <= (spatial_lat(jj)+0.5*res_latlon) ) then

!        if(masterproc) print*, 'MSG spatial_mean_2d: found lat ', (spatial_lat(jj)-0.5*res_latlon), &
!         var2d(i,j)%lat, (spatial_lat(jj)+0.5*res_latlon)

        do ii = 1, imax
          if( var2d(i,j)%lon >  (spatial_lon(ii)-0.5*res_latlon) .and. &
              var2d(i,j)%lon <= (spatial_lon(ii)+0.5*res_latlon) ) then

!             if(masterproc) print*, 'MSG spatial_mean_2d: found lon ',(spatial_lon(ii)-0.5*res_latlon),&
!                var2d(i,j)%lon, (spatial_lon(ii)+0.5*res_latlon)

           do kk = 1, kmax

               do k = 1, bnd(3)
                  if( var3d(i,j,k)%press <= pressure_lev(kk) .and. &
                      var3d(i,j,k)%press >  pressure_lev(kk+1)  ) then

                    spatial_cnt(ii,jj,kk) = spatial_cnt(ii,jj,kk) + 1.e0  !update total count
                    spatial_sum(ii,jj,kk)%rain_rate  = spatial_sum(ii,jj,kk)%rain_rate + var3d(i,j,k)%rain_rate
                    spatial_sum(ii,jj,kk)%t_air      = spatial_sum(ii,jj,kk)%t_air     + var3d(i,j,k)%t_air
                    spatial_sum(ii,jj,kk)%sh         = spatial_sum(ii,jj,kk)%sh        + var3d(i,j,k)%sh
                    spatial_sum(ii,jj,kk)%u          = spatial_sum(ii,jj,kk)%u         + var3d(i,j,k)%u
                    spatial_sum(ii,jj,kk)%v          = spatial_sum(ii,jj,kk)%v         + var3d(i,j,k)%v
                    spatial_sum(ii,jj,kk)%w          = spatial_sum(ii,jj,kk)%w         + var3d(i,j,k)%w
                    spatial_sum(ii,jj,kk)%lw_heat    = spatial_sum(ii,jj,kk)%lw_heat   + var3d(i,j,k)%lw_heat
                    spatial_sum(ii,jj,kk)%sw_heat    = spatial_sum(ii,jj,kk)%sw_heat   + var3d(i,j,k)%sw_heat
                    spatial_sum(ii,jj,kk)%tot_cloud  = spatial_sum(ii,jj,kk)%tot_cloud + var3d(i,j,k)%tot_cloud

                  endif
               enddo

           enddo
           endif
        enddo

        endif within_grid
      enddo


   enddo ; enddo


  endif

 case('finish') !---------------------------------------------------------------------------------------

#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', spatial_cnt  )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%t_air )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%sh   )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%u )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%v )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%w )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%lw_heat )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%sw_heat )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%tot_cloud )

  MPI_MASTERPROC: if(masterproc) then
#endif

  do kk = 1, kmax ; do jj = 1, jmax ; do ii = 1, imax 

    if( spatial_cnt(ii,jj,kk) > 0. ) then

     spatial_mean(ii,jj,kk)%rain_rate = spatial_sum(ii,jj,kk)%rain_rate / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%t_air     = spatial_sum(ii,jj,kk)%t_air / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%sh        = spatial_sum(ii,jj,kk)%sh   / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%u         = spatial_sum(ii,jj,kk)%u    / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%v         = spatial_sum(ii,jj,kk)%v    / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%w         = spatial_sum(ii,jj,kk)%w    / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%lw_heat   = spatial_sum(ii,jj,kk)%lw_heat / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%sw_heat   = spatial_sum(ii,jj,kk)%sw_heat / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%tot_cloud = spatial_sum(ii,jj,kk)%tot_cloud / spatial_cnt(ii,jj,kk)

    endif

  enddo ; enddo ; enddo


!
! write out cloud frequency diurnal profile
!
 write(yyyy,"(I4.4)") int( sdsu_yyyy )
 write(mm  ,"(I2.2)") int( sdsu_mm   )
 write(dd  ,"(I2.2)") int( sdsu_dd   )


!
! write out cloud frequency zonal profile
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM3D_'//yyyy//'_'//mm//'_'//dd&
                //trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=imax*jmax*kmax*4)
 print("(a)"),'MSG spatial_mean3d: see output ->',trim(sdsu_io_file)

 irec = 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%rain_rate , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%t_air     , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%sh        , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%u         , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%v         , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%w         , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%lw_heat   , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%sw_heat   , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%tot_cloud , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1

 close(io)


!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM3D_'//yyyy//'_'//mm//'_'//dd&
                //trim(output_suffix)//'.ctl'
 print("(a)"),'MSG spatial_mean3d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'spatial_mean_CRM3D_'//yyyy//'_'//mm//'_'//dd&
                //trim(output_suffix)//'.grads'
 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   spatial mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    ',imax,' linear ',minlon,' ',res_latlon
 write(io,*) 'ydef    ',jmax,' linear ',minlat,' ',res_latlon
 write(io,*) 'zdef    ',play,' levels ', pressure_lay(1:play)
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars ',irec-1
 write(io,*) 'rain   ',play,'  0  rain rate [mm/hr]'
 write(io,*) 'tair   ',play,'  0  air temperature [degK]'
 write(io,*) 'sh   ',play ,'  0  specific humidity [g/g]'
 write(io,*) 'u    ',play,'  0  u wind [m/s]'
 write(io,*) 'v    ',play,'  0  v wind [m/s]'
 write(io,*) 'w    ',play,'  0  w wind [m/s]'
 write(io,*) 'lwheat',play,'  0  LW heating rate [K/day]'
 write(io,*) 'swheat',play,'  0  SW heating rate [K/day]'
 write(io,*) 'totcloud',play,'  0  total cloud + ice [g/m3]'

 write(io,*) 'endvars '
 close(io)


#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG spatial_mean_daily_3d: There is no such proc')
 end select


 return
 end subroutine spatial_mean_daily_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine spatial_mean_hourly_3d(proc,var3d,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive spatial mean hourly map of 3D parameters.
!           
! History:
! 06/2017  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 character*(*) :: proc ! process
 type ( atmos_parameter ), dimension(:,:,:), optional :: var3d ! input surface parameter assumed 3d array
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 3d array

 integer :: bnd(3)  !upper memory bound for 2D array

!
! statisticsl parameters
!
 type( atmos_parameter ), allocatable, dimension(:,:,:),save :: spatial_sum, spatial_mean
 real(sdsu_fps), allocatable, dimension(:,:,:),save :: spatial_cnt
 real(sdsu_fps), allocatable, dimension(:),save :: spatial_lat, spatial_lon

 integer,save :: ibin_max
 integer :: i,j,k, ibin, ierr, irec
 integer,parameter :: io = 149


 integer,parameter :: play = 20, plev = 21
 real(sdsu_fps),parameter :: pressure_lay(play) =    (/1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,500.,&
                                                     450.,400.,350.,300.,250.,200.,150.,100.,50./)
 real(sdsu_fps),parameter :: pressure_lev(plev) = (/1025., 975.,925.,875.,825.,775.,725.,675.,625.,575.,525.,475.,&
                                                        425.,375.,325.,275.,225.,175.,125.,75.,25./)

 integer,save :: imax, jmax, kmax
 integer :: ii, jj, kk

 real(sdsu_fps) :: pr

 character :: yyyy*4, mm*2, dd*2, hh*2


 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

   if(masterproc) print*, 'MSG spatial_mean_hourly_3d: init'

   imax = int( ( maxlon - minlon ) / res_latlon ) + 1
   jmax = int( ( maxlat - minlat ) / res_latlon ) + 1
   kmax = play

   if( imax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_hourly_3d: maxlon is smaller than minlon')
   if( jmax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_hourly_3d: maxlat is smaller than minlat')

   if( allocated(spatial_sum) ) deallocate(spatial_sum, spatial_cnt, spatial_mean, spatial_lat, spatial_lon)

   allocate( spatial_sum(imax,jmax,kmax)  , &
             spatial_cnt(imax,jmax,kmax)  , &
             spatial_mean(imax,jmax,kmax) , &
             spatial_lat(jmax)  , &
             spatial_lon(imax)  , &
             stat=ierr )

   ! initialized
   do j = 1, jmax
      spatial_lat(j) = minlat + res_latlon * real(j-1) ![deg]
   enddo
   do i = 1, imax
      spatial_lon(i) = minlon + res_latlon * real(i-1) ![deg]
   enddo

  !
  ! zero out other parameters
  !
  spatial_cnt = 0.

  spatial_sum%rain_rate      = 0. ; spatial_sum%t_air       = 0.
  spatial_sum%sh             = 0. ; spatial_sum%u           = 0.
  spatial_sum%v              = 0. ; spatial_sum%w           = 0.
  spatial_sum%lw_heat        = 0. ; spatial_sum%sw_heat     = 0.
  spatial_sum%tot_cloud      = 0.

  spatial_mean%rain_rate      = undefined ; spatial_mean%t_air       = undefined
  spatial_mean%sh             = undefined ; spatial_mean%u           = undefined
  spatial_mean%v              = undefined ; spatial_mean%w           = undefined
  spatial_mean%lw_heat        = undefined ; spatial_mean%sw_heat     = undefined
  spatial_mean%tot_cloud      = undefined

 case('add') !-------------------------------------------------------------------------------------

  if(masterproc) print*, 'MSG spatial_mean_hourly_3d: add'

  if(present(var3d) ) then

   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var3d)

   do j = 1,bnd(2) ; do i = 1,bnd(1)

      !
      ! loop for latlon_grid domain. 
      !
      do jj = 1, jmax
        within_grid: if( var2d(i,j)%lat >  (spatial_lat(jj)-0.5*res_latlon) .and. &
                         var2d(i,j)%lat <= (spatial_lat(jj)+0.5*res_latlon) ) then

!        if(masterproc) print*, 'MSG spatial_mean_2d: found lat ', (spatial_lat(jj)-0.5*res_latlon), &
!         var2d(i,j)%lat, (spatial_lat(jj)+0.5*res_latlon)

        do ii = 1, imax
          if( var2d(i,j)%lon >  (spatial_lon(ii)-0.5*res_latlon) .and. &
              var2d(i,j)%lon <= (spatial_lon(ii)+0.5*res_latlon) ) then

!             if(masterproc) print*, 'MSG spatial_mean_2d: found lon ',(spatial_lon(ii)-0.5*res_latlon),&
!                var2d(i,j)%lon, (spatial_lon(ii)+0.5*res_latlon)

           do kk = 1, kmax

               do k = 1, bnd(3)
                  if( var3d(i,j,k)%press <= pressure_lev(kk) .and. &
                      var3d(i,j,k)%press >  pressure_lev(kk+1)  ) then

                    spatial_cnt(ii,jj,kk) = spatial_cnt(ii,jj,kk) + 1.e0  !update total count
                    spatial_sum(ii,jj,kk)%rain_rate  = spatial_sum(ii,jj,kk)%rain_rate + var3d(i,j,k)%rain_rate
                    spatial_sum(ii,jj,kk)%t_air      = spatial_sum(ii,jj,kk)%t_air     + var3d(i,j,k)%t_air
                    spatial_sum(ii,jj,kk)%sh         = spatial_sum(ii,jj,kk)%sh        + var3d(i,j,k)%sh
                    spatial_sum(ii,jj,kk)%u          = spatial_sum(ii,jj,kk)%u         + var3d(i,j,k)%u
                    spatial_sum(ii,jj,kk)%v          = spatial_sum(ii,jj,kk)%v         + var3d(i,j,k)%v
                    spatial_sum(ii,jj,kk)%w          = spatial_sum(ii,jj,kk)%w         + var3d(i,j,k)%w
                    spatial_sum(ii,jj,kk)%lw_heat    = spatial_sum(ii,jj,kk)%lw_heat   + var3d(i,j,k)%lw_heat
                    spatial_sum(ii,jj,kk)%sw_heat    = spatial_sum(ii,jj,kk)%sw_heat   + var3d(i,j,k)%sw_heat
                    spatial_sum(ii,jj,kk)%tot_cloud  = spatial_sum(ii,jj,kk)%tot_cloud + var3d(i,j,k)%tot_cloud

                  endif
               enddo

           enddo
           endif
        enddo

        endif within_grid
      enddo


   enddo ; enddo


  endif

 case('finish') !---------------------------------------------------------------------------------------

#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', spatial_cnt  )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%t_air )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%sh   )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%u )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%v )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%w )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%lw_heat )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%sw_heat )
  call mpi_sdsu_reduce( 'SUM', spatial_sum%tot_cloud )

  MPI_MASTERPROC: if(masterproc) then
#endif

  do kk = 1, kmax ; do jj = 1, jmax ; do ii = 1, imax 

    if( spatial_cnt(ii,jj,kk) > 0. ) then

     spatial_mean(ii,jj,kk)%rain_rate = spatial_sum(ii,jj,kk)%rain_rate / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%t_air     = spatial_sum(ii,jj,kk)%t_air / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%sh        = spatial_sum(ii,jj,kk)%sh   / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%u         = spatial_sum(ii,jj,kk)%u    / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%v         = spatial_sum(ii,jj,kk)%v    / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%w         = spatial_sum(ii,jj,kk)%w    / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%lw_heat   = spatial_sum(ii,jj,kk)%lw_heat / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%sw_heat   = spatial_sum(ii,jj,kk)%sw_heat / spatial_cnt(ii,jj,kk)
     spatial_mean(ii,jj,kk)%tot_cloud = spatial_sum(ii,jj,kk)%tot_cloud / spatial_cnt(ii,jj,kk)

    endif

  enddo ; enddo ; enddo


!
! write out cloud frequency diurnal profile
!
 write(yyyy,"(I4.4)") int( sdsu_yyyy )
 write(mm  ,"(I2.2)") int( sdsu_mm   )
 write(dd  ,"(I2.2)") int( sdsu_dd   )
 write(hh  ,"(I2.2)") int( sdsu_hh   )

!
! write out cloud frequency zonal profile
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM3D_'//yyyy//'_'//mm//'_'//dd//'_'//hh&
                //trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=imax*jmax*kmax*4)
 print("(a)"),'MSG spatial_mean3d: see output ->',trim(sdsu_io_file)

 irec = 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%rain_rate , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%t_air     , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%sh        , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%u         , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%v         , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%w         , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%lw_heat   , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%sw_heat   , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1
 write(io,rec=irec) ((( spatial_mean(ii,jj,kk)%tot_cloud , ii=1,imax ),jj=1,jmax),kk=1,kmax) ; irec = irec + 1

 close(io)


!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'spatial_mean_CRM3D_'//yyyy//'_'//mm//'_'//dd//'_'//hh&
                //trim(output_suffix)//'.ctl'
 print("(a)"),'MSG spatial_mean3d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'spatial_mean_CRM3D_'//yyyy//'_'//mm//'_'//dd//'_'//hh&
                //trim(output_suffix)//'.grads'
 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   spatial mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    ',imax,' linear ',minlon,' ',res_latlon
 write(io,*) 'ydef    ',jmax,' linear ',minlat,' ',res_latlon
 write(io,*) 'zdef    ',play,' levels ', pressure_lay(1:play)
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars ',irec-1
 write(io,*) 'rain   ',play,'  0  rain rate [mm/hr]'
 write(io,*) 'tair   ',play,'  0  air temperature [degK]'
 write(io,*) 'sh   ',play ,'  0  specific humidity [g/g]'
 write(io,*) 'u    ',play,'  0  u wind [m/s]'
 write(io,*) 'v    ',play,'  0  v wind [m/s]'
 write(io,*) 'w    ',play,'  0  w wind [m/s]'
 write(io,*) 'lwheat',play,'  0  LW heating rate [K/day]'
 write(io,*) 'swheat',play,'  0  SW heating rate [K/day]'
 write(io,*) 'totcloud',play,'  0  total cloud + ice [g/m3]'

 write(io,*) 'endvars '
 close(io)


#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG spatial_mean_hourly_3d: There is no such proc')
 end select


 return
 end subroutine spatial_mean_hourly_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine zonal_mean_2d(proc,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive zonal mean composite of 2D parameters.
!           
! History:
! 06/2016  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------

 character*(*) :: proc ! process
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 2d array

 integer :: bnd(2)  !upper memory bound for 2D array

 integer,save :: ibin_max

!
! statisticsl parameters
!
 type( surface_parameter ), allocatable, dimension(:),save :: zonal_sum,zonal_mean
 real(sdsu_fps), allocatable, dimension(:),save :: zonal_lat, zonal_cnt
 integer,save :: jmax

 integer :: i,j, ibin, ierr, irec
 integer,parameter :: io = 151

 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

 if(masterproc) print*, 'MSG zonal_mean_2d: init'

!
! maximum bin of zone
!
 jmax  =int( (maxlat - minlat) / res_latlon ) + 1

 if( jmax < 1 ) call stop_sdsu(myrank,'MSG spatial_mean_2d: maxlat is smaller than minlat')
 
!
! allocate
!
 if( .not. allocated(zonal_sum) ) then
     allocate( zonal_sum (jmax), &
               zonal_mean(jmax), &
               zonal_cnt (jmax), &
               zonal_lat (jmax), &
               stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG zonal_mean_2d: allocation error -> Terminate program.')
 endif

!
! assign latitude bin
!
  do ibin = 1, ibin_max 
     zonal_lat(ibin) = minlat + res_latlon * real(ibin) - 0.5 * res_latlon
  enddo



!
! zero out other parameters
!
 zonal_cnt = 0.

 zonal_sum%rain_rate      = 0. ; zonal_sum%h2o_col       = 0.  
 zonal_sum%path_cloud_liq = 0. ; zonal_sum%path_rain_liq = 0.
 zonal_sum%path_cloud_ice = 0. ; zonal_sum%path_rain_ice = 0.
 zonal_sum%toa_lw_up      = 0. ; zonal_sum%toa_lw_dn     = 0.
 zonal_sum%toa_sw_up      = 0. ; zonal_sum%toa_sw_dn     = 0.
 zonal_sum%boa_lw_up      = 0. ; zonal_sum%boa_lw_dn     = 0.
 zonal_sum%boa_sw_up      = 0. ; zonal_sum%boa_sw_dn     = 0.

 zonal_mean%rain_rate      = undefined ; zonal_mean%h2o_col       = undefined
 zonal_mean%path_cloud_liq = undefined ; zonal_mean%path_rain_liq = undefined
 zonal_mean%path_cloud_ice = undefined ; zonal_mean%path_rain_ice = undefined
 zonal_mean%toa_lw_up      = undefined ; zonal_mean%toa_lw_dn     = undefined
 zonal_mean%toa_sw_up      = undefined ; zonal_mean%toa_sw_dn     = undefined
 zonal_mean%boa_lw_up      = undefined ; zonal_mean%boa_lw_dn     = undefined
 zonal_mean%boa_sw_up      = undefined ; zonal_mean%boa_sw_dn     = undefined


 case('add') !-------------------------------------------------------------------------------------


 if(present(var2d) ) then


   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var2d)


   do j = 1,bnd(2) ; do i = 1,bnd(1)

      do ibin = 1, ibin_max
         if( var2d(i,j)%lat >  zonal_lat(ibin) -  0.5 * res_latlon  .and. &
             var2d(i,j)%lat <= zonal_lat(ibin) +  0.5 * res_latlon  ) then

             zonal_cnt(ibin) = zonal_cnt(ibin) + 1.e0  !update total count

             zonal_sum(ibin)%rain_rate      = zonal_sum(ibin)%rain_rate + var2d(i,j)%rain_rate
             zonal_sum(ibin)%h2o_col        = zonal_sum(ibin)%h2o_col   + var2d(i,j)%h2o_col
             zonal_sum(ibin)%path_cloud_liq = zonal_sum(ibin)%path_cloud_liq + var2d(i,j)%path_cloud_liq
             zonal_sum(ibin)%path_rain_liq  = zonal_sum(ibin)%path_rain_liq  + var2d(i,j)%path_rain_liq
             zonal_sum(ibin)%path_cloud_ice = zonal_sum(ibin)%path_cloud_ice + var2d(i,j)%path_cloud_ice
             zonal_sum(ibin)%path_rain_ice  = zonal_sum(ibin)%path_rain_ice  + var2d(i,j)%path_rain_ice
             zonal_sum(ibin)%toa_lw_up      = zonal_sum(ibin)%toa_lw_up + var2d(i,j)%toa_lw_up
             zonal_sum(ibin)%toa_lw_dn      = zonal_sum(ibin)%toa_lw_dn + var2d(i,j)%toa_lw_dn
             zonal_sum(ibin)%toa_sw_up      = zonal_sum(ibin)%toa_sw_up + var2d(i,j)%toa_sw_up
             zonal_sum(ibin)%toa_sw_dn      = zonal_sum(ibin)%toa_sw_dn + var2d(i,j)%toa_sw_dn
             zonal_sum(ibin)%boa_lw_up      = zonal_sum(ibin)%boa_lw_up + var2d(i,j)%boa_lw_up
             zonal_sum(ibin)%boa_lw_dn      = zonal_sum(ibin)%boa_lw_dn + var2d(i,j)%boa_lw_dn
             zonal_sum(ibin)%boa_sw_up      = zonal_sum(ibin)%boa_sw_up + var2d(i,j)%boa_sw_up
             zonal_sum(ibin)%boa_sw_dn      = zonal_sum(ibin)%boa_sw_dn + var2d(i,j)%boa_sw_dn

         endif
      enddo

   enddo ; enddo



 endif

 case('finish') !---------------------------------------------------------------------------------------


#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', zonal_cnt  )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%h2o_col   )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%path_cloud_liq )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%path_rain_liq )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%path_cloud_ice )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%path_rain_ice )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%toa_lw_up )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%toa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%toa_sw_up )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%toa_sw_dn )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%boa_lw_up )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%boa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%boa_sw_up )
  call mpi_sdsu_reduce( 'SUM', zonal_sum%boa_sw_dn )

  MPI_MASTERPROC: if(masterproc) then
#endif

!
! derive zonal mean values
! 
 do ibin = 1, ibin_max

    if( zonal_cnt(ibin) > 0. ) then

     zonal_mean(ibin)%rain_rate      = zonal_sum(ibin)%rain_rate / zonal_cnt(ibin)
     zonal_mean(ibin)%h2o_col        = zonal_sum(ibin)%h2o_col   / zonal_cnt(ibin)
     zonal_mean(ibin)%path_cloud_liq = zonal_sum(ibin)%path_cloud_liq / zonal_cnt(ibin)
     zonal_mean(ibin)%path_rain_liq  = zonal_sum(ibin)%path_rain_liq  / zonal_cnt(ibin)
     zonal_mean(ibin)%path_cloud_ice = zonal_sum(ibin)%path_cloud_ice / zonal_cnt(ibin)
     zonal_mean(ibin)%path_rain_ice  = zonal_sum(ibin)%path_rain_ice  / zonal_cnt(ibin)
     zonal_mean(ibin)%toa_lw_up      = zonal_sum(ibin)%toa_lw_up / zonal_cnt(ibin)
     zonal_mean(ibin)%toa_lw_dn      = zonal_sum(ibin)%toa_lw_dn / zonal_cnt(ibin)
     zonal_mean(ibin)%toa_sw_up      = zonal_sum(ibin)%toa_sw_up / zonal_cnt(ibin)
     zonal_mean(ibin)%toa_sw_dn      = zonal_sum(ibin)%toa_sw_dn / zonal_cnt(ibin)
     zonal_mean(ibin)%boa_lw_up      = zonal_sum(ibin)%boa_lw_up / zonal_cnt(ibin)
     zonal_mean(ibin)%boa_lw_dn      = zonal_sum(ibin)%boa_lw_dn / zonal_cnt(ibin)
     zonal_mean(ibin)%boa_sw_up      = zonal_sum(ibin)%boa_sw_up / zonal_cnt(ibin)
     zonal_mean(ibin)%boa_sw_dn      = zonal_sum(ibin)%boa_sw_dn / zonal_cnt(ibin)

    endif

 enddo


!
! write out cloud frequency zonal profile
!

 sdsu_io_file = trim(sdsu_dir_output)//'zonal_mean_CRM2D'//trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=ibin_max*4)
 print("(a)"),'MSG zonal_mean2d: see output ->',trim(sdsu_io_file)

 irec = 1
 write(io,rec=irec) ( zonal_mean(ibin)%rain_rate      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%h2o_col        , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%path_cloud_liq , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%path_rain_liq  , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%path_cloud_ice , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%path_rain_ice  , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%toa_lw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%toa_lw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%toa_sw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%toa_sw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%boa_lw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%boa_lw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%boa_sw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( zonal_mean(ibin)%boa_sw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1

 close(io)


!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'zonal_mean_CRM2D'//trim(output_suffix)//'.ctl'
 print("(a)"),'MSG zonal_mean2d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'zonal_mean_CRM2D'//trim(output_suffix)//'.grads'
 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   zonal mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    1 linear 1. 1'
 write(io,*) 'ydef    ',ibin_max,' linear ',minlat,' ',res_latlon
 write(io,*) 'zdef    1 linear 0. 0.1'
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars ',irec-1
 write(io,*) 'rain    0  0  rain rate [mm/hr]'

 write(io,*) 'h2ocol  0  0  column water vapro [kg/m2]'
 write(io,*) 'clwp    0  0  cloud liquid water path [kg/m2]'
 write(io,*) 'lrwp    0  0  liquid rain water path [kg/m2]'
 write(io,*) 'ciwp    0  0  cloud ice water path [kg/m2]'
 write(io,*) 'irwp    0  0  ice rain water path [kg/m2]'

 write(io,*) 'lwut   0 0  TOA LW upwelling   flux [W/m2]'
 write(io,*) 'lwdt   0 0  TOA LW downwelling flux [W/m2]'
 write(io,*) 'swut   0 0  TOA SW upwelling flux   [W/m2]'
 write(io,*) 'swdt   0 0  TOA SW downwelling flux [W/m2]'
 write(io,*) 'lwus   0 0  Surface LW upwelling   flux [W/m2]'
 write(io,*) 'lwds   0 0  Surface LW downwelling flux [W/m2]'
 write(io,*) 'swus   0 0  Surface SW upwelling   flux [W/m2]'
 write(io,*) 'swds   0 0  Surface SW downwelling flux [W/m2]'

 write(io,*) 'endvars '
 close(io)


#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG zonal_mean_2d: There is no such proc')
 end select


 return
 end subroutine zonal_mean_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine diurnal_mean_2d(proc,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive diurnal composite of 2D parameters.
!           
! History:
! 06/2016  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 character*(*) :: proc ! process
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 2d array

 integer :: bnd(2)  !upper memory bound for 2D array


!
! statisticsl parameters
!
 type( surface_parameter ), allocatable, dimension(:),save :: diurnal_sum,diurnal_mean
 real(sdsu_fps), allocatable, dimension(:),save :: diurnal_lst, diurnal_cnt
 integer,save :: ibin_max
 real(sdsu_fps) :: hour_loc

 integer :: i,j, ibin, ierr, irec
 integer,parameter :: io = 154

 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------
 if(masterproc) print*,'MSG diurnal_mean_2d: init'

!
! maximum bin of zone
!
  ibin_max  = 24

!
! allocate
!
 if( .not. allocated(diurnal_sum) ) then
     allocate( diurnal_sum (ibin_max), &
               diurnal_mean(ibin_max), &
               diurnal_cnt (ibin_max), &
               diurnal_lst (ibin_max), &
               stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG diurnal_mean_2d: allocation error -> Terminate program.')
 endif

!
! assign latitude bin
!
  do ibin = 1, ibin_max
     diurnal_lst(ibin) = real(ibin) -1.   !0~23LST (starting -0.5~0.5LST and ending 22.5~23.5LST
  enddo



!
! zero out other parameters
!
 diurnal_cnt = 0.

 diurnal_sum%rain_rate      = 0. ; diurnal_sum%h2o_col       = 0.
 diurnal_sum%path_cloud_liq = 0. ; diurnal_sum%path_rain_liq = 0.
 diurnal_sum%path_cloud_ice = 0. ; diurnal_sum%path_rain_ice = 0.
 diurnal_sum%toa_lw_up      = 0. ; diurnal_sum%toa_lw_dn     = 0.
 diurnal_sum%toa_sw_up      = 0. ; diurnal_sum%toa_sw_dn     = 0.
 diurnal_sum%boa_lw_up      = 0. ; diurnal_sum%boa_lw_dn     = 0.
 diurnal_sum%boa_sw_up      = 0. ; diurnal_sum%boa_sw_dn     = 0.

 diurnal_mean%rain_rate      = undefined ; diurnal_mean%h2o_col       = undefined
 diurnal_mean%path_cloud_liq = undefined ; diurnal_mean%path_rain_liq = undefined
 diurnal_mean%path_cloud_ice = undefined ; diurnal_mean%path_rain_ice = undefined
 diurnal_mean%toa_lw_up      = undefined ; diurnal_mean%toa_lw_dn     = undefined
 diurnal_mean%toa_sw_up      = undefined ; diurnal_mean%toa_sw_dn     = undefined
 diurnal_mean%boa_lw_up      = undefined ; diurnal_mean%boa_lw_dn     = undefined
 diurnal_mean%boa_sw_up      = undefined ; diurnal_mean%boa_sw_dn     = undefined


 case('add') !-------------------------------------------------------------------------------------

 if(masterproc) print*,'MSG diurnal_mean_2d: add'

 if(present(var2d) ) then



   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var2d)


   do j = 1,bnd(2) ; do i = 1,bnd(1)

    land_if: if( var2d(i,j)%iland == 1 ) then !over land

    latlon_if: if( var2d(i,j)%lat >= minlat .and. var2d(i,j)%lat <= maxlat .and.&
                   var2d(i,j)%lon >= minlon .and. var2d(i,j)%lon <= maxlon ) then


         hour_loc = sdsu_gmt + var2d(i,j)%lon/15.
         if( hour_loc < -0.5 ) hour_loc = hour_loc + 24.
         if( hour_loc > 23.5 ) hour_loc = hour_loc -24.


      do ibin = 1, ibin_max
         if( hour_loc >  diurnal_lst(ibin) -  0.5  .and. &
             hour_loc <= diurnal_lst(ibin) +  0.5  ) then

             diurnal_cnt(ibin) = diurnal_cnt(ibin) + 1.e0  !update total count

             diurnal_sum(ibin)%rain_rate      = diurnal_sum(ibin)%rain_rate + var2d(i,j)%rain_rate
             diurnal_sum(ibin)%h2o_col        = diurnal_sum(ibin)%h2o_col   + var2d(i,j)%h2o_col
             diurnal_sum(ibin)%path_cloud_liq = diurnal_sum(ibin)%path_cloud_liq + var2d(i,j)%path_cloud_liq
             diurnal_sum(ibin)%path_rain_liq  = diurnal_sum(ibin)%path_rain_liq  + var2d(i,j)%path_rain_liq
             diurnal_sum(ibin)%path_cloud_ice = diurnal_sum(ibin)%path_cloud_ice + var2d(i,j)%path_cloud_ice
             diurnal_sum(ibin)%path_rain_ice  = diurnal_sum(ibin)%path_rain_ice  + var2d(i,j)%path_rain_ice
             diurnal_sum(ibin)%toa_lw_up      = diurnal_sum(ibin)%toa_lw_up + var2d(i,j)%toa_lw_up
             diurnal_sum(ibin)%toa_lw_dn      = diurnal_sum(ibin)%toa_lw_dn + var2d(i,j)%toa_lw_dn
             diurnal_sum(ibin)%toa_sw_up      = diurnal_sum(ibin)%toa_sw_up + var2d(i,j)%toa_sw_up
             diurnal_sum(ibin)%toa_sw_dn      = diurnal_sum(ibin)%toa_sw_dn + var2d(i,j)%toa_sw_dn
             diurnal_sum(ibin)%boa_lw_up      = diurnal_sum(ibin)%boa_lw_up + var2d(i,j)%boa_lw_up
             diurnal_sum(ibin)%boa_lw_dn      = diurnal_sum(ibin)%boa_lw_dn + var2d(i,j)%boa_lw_dn
             diurnal_sum(ibin)%boa_sw_up      = diurnal_sum(ibin)%boa_sw_up + var2d(i,j)%boa_sw_up
             diurnal_sum(ibin)%boa_sw_dn      = diurnal_sum(ibin)%boa_sw_dn + var2d(i,j)%boa_sw_dn

         endif
      enddo

    endif latlon_if
    endif land_if

   enddo ; enddo

 endif


 case('finish') !---------------------------------------------------------------------------------------



#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', diurnal_cnt  )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%h2o_col   )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%path_cloud_liq )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%path_rain_liq )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%path_cloud_ice )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%path_rain_ice )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%toa_lw_up )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%toa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%toa_sw_up )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%toa_sw_dn )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%boa_lw_up )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%boa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%boa_sw_up )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%boa_sw_dn )

  MPI_MASTERPROC: if(masterproc) then
#endif

!
! derive diurnal mean values
! 
 do ibin = 1, ibin_max

    if( diurnal_cnt(ibin) > 0. ) then

     diurnal_mean(ibin)%rain_rate      = diurnal_sum(ibin)%rain_rate / diurnal_cnt(ibin)
     diurnal_mean(ibin)%h2o_col        = diurnal_sum(ibin)%h2o_col   / diurnal_cnt(ibin)
     diurnal_mean(ibin)%path_cloud_liq = diurnal_sum(ibin)%path_cloud_liq / diurnal_cnt(ibin)
     diurnal_mean(ibin)%path_rain_liq  = diurnal_sum(ibin)%path_rain_liq  / diurnal_cnt(ibin)
     diurnal_mean(ibin)%path_cloud_ice = diurnal_sum(ibin)%path_cloud_ice / diurnal_cnt(ibin)
     diurnal_mean(ibin)%path_rain_ice  = diurnal_sum(ibin)%path_rain_ice  / diurnal_cnt(ibin)
     diurnal_mean(ibin)%toa_lw_up      = diurnal_sum(ibin)%toa_lw_up / diurnal_cnt(ibin)
     diurnal_mean(ibin)%toa_lw_dn      = diurnal_sum(ibin)%toa_lw_dn / diurnal_cnt(ibin)
     diurnal_mean(ibin)%toa_sw_up      = diurnal_sum(ibin)%toa_sw_up / diurnal_cnt(ibin)
     diurnal_mean(ibin)%toa_sw_dn      = diurnal_sum(ibin)%toa_sw_dn / diurnal_cnt(ibin)
     diurnal_mean(ibin)%boa_lw_up      = diurnal_sum(ibin)%boa_lw_up / diurnal_cnt(ibin)
     diurnal_mean(ibin)%boa_lw_dn      = diurnal_sum(ibin)%boa_lw_dn / diurnal_cnt(ibin)
     diurnal_mean(ibin)%boa_sw_up      = diurnal_sum(ibin)%boa_sw_up / diurnal_cnt(ibin)
     diurnal_mean(ibin)%boa_sw_dn      = diurnal_sum(ibin)%boa_sw_dn / diurnal_cnt(ibin)

    endif

 enddo


!
! write out cloud frequency diurnal profile
!


 sdsu_io_file = trim(sdsu_dir_output)//'diurnal_mean_CRM2D'//trim(output_suffix)//'.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=ibin_max*4)
 print("(a)"),'MSG diurnal_mean2d: see output ->',trim(sdsu_io_file)

 irec = 1
 write(io,rec=irec) ( diurnal_mean(ibin)%rain_rate      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%h2o_col        , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%path_cloud_liq , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%path_rain_liq  , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%path_cloud_ice , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%path_rain_ice  , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%toa_lw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%toa_lw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%toa_sw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%toa_sw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%boa_lw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%boa_lw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%boa_sw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%boa_sw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1

 close(io)


!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'diurnal_mean_CRM2D'//trim(output_suffix)//'.ctl'
 print("(a)"),'MSG diurnal_mean2d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'diurnal_mean_CRM2D'//trim(output_suffix)//'.grads'
 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   diurnal mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    1 linear 1. 1'
 write(io,*) 'ydef    ',ibin_max,' linear 0   1. '
 write(io,*) 'zdef    1 linear 0. 0.1'
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars ',irec-1
 write(io,*) 'rain    0  0  rain rate [mm/hr]'

 write(io,*) 'h2ocol  0  0  column water vaper [kg/m2]'
 write(io,*) 'clwp    0  0  cloud liquid water path [kg/m2]'
 write(io,*) 'lrwp    0  0  liquid rain water path [kg/m2]'
 write(io,*) 'ciwp    0  0  cloud ice water path [kg/m2]'
 write(io,*) 'irwp    0  0  ice rain water path [kg/m2]'

 write(io,*) 'lwut   0 0  TOA LW upwelling   flux [W/m2]'
 write(io,*) 'lwdt   0 0  TOA LW downwelling flux [W/m2]'
 write(io,*) 'swut   0 0  TOA SW upwelling flux   [W/m2]'
 write(io,*) 'swdt   0 0  TOA SW downwelling flux [W/m2]'
 write(io,*) 'lwus   0 0  Surface LW upwelling   flux [W/m2]'
 write(io,*) 'lwds   0 0  Surface LW downwelling flux [W/m2]'
 write(io,*) 'swus   0 0  Surface SW upwelling   flux [W/m2]'
 write(io,*) 'swds   0 0  Surface SW downwelling flux [W/m2]'

 write(io,*) 'endvars '
 close(io)


#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG diurnal_mean_2d: There is no such proc')
 end select


 return
 end subroutine diurnal_mean_2d


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine diurnal_daily_2d(proc,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive diurnal composite of 2D parameters.
!           
! History:
! 06/2016  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 character*(*) :: proc ! process
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 2d array

 integer :: bnd(2)  !upper memory bound for 2D array

 character :: yyyy*4, mm*2, dd*2

!
! statisticsl parameters
!
 type( surface_parameter ), allocatable, dimension(:),save :: diurnal_sum,diurnal_mean
 real(sdsu_fps), allocatable, dimension(:),save :: diurnal_lst, diurnal_cnt
 integer,save :: ibin_max
 real(sdsu_fps) :: hour_loc

 integer :: i,j, ibin, ierr, irec
 integer,parameter :: io = 153

 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------

 if( masterproc ) print*,'MSG diurnal_daily_2d: init'
!
! maximum bin of zone
!
  ibin_max  = 24

!
! allocate
!
 if( .not. allocated(diurnal_sum) ) then
     allocate( diurnal_sum (ibin_max), &
               diurnal_mean(ibin_max), &
               diurnal_cnt (ibin_max), &
               diurnal_lst (ibin_max), &
               stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG diurnal_daily_2d: allocation error -> Terminate program.')
 endif

!
! assign latitude bin
!
  do ibin = 1, ibin_max
     diurnal_lst(ibin) = real(ibin) -1.   !0~23LST (starting -0.5~0.5LST and ending 22.5~23.5LST
  enddo


!
! zero out other parameters
!
 diurnal_cnt = 0.

 diurnal_sum%rain_rate      = 0. ; diurnal_sum%h2o_col       = 0.
 diurnal_sum%path_cloud_liq = 0. ; diurnal_sum%path_rain_liq = 0.
 diurnal_sum%path_cloud_ice = 0. ; diurnal_sum%path_rain_ice = 0.
 diurnal_sum%toa_lw_up      = 0. ; diurnal_sum%toa_lw_dn     = 0.
 diurnal_sum%toa_sw_up      = 0. ; diurnal_sum%toa_sw_dn     = 0.
 diurnal_sum%boa_lw_up      = 0. ; diurnal_sum%boa_lw_dn     = 0.
 diurnal_sum%boa_sw_up      = 0. ; diurnal_sum%boa_sw_dn     = 0.

 diurnal_mean%rain_rate      = undefined ; diurnal_mean%h2o_col       = undefined
 diurnal_mean%path_cloud_liq = undefined ; diurnal_mean%path_rain_liq = undefined
 diurnal_mean%path_cloud_ice = undefined ; diurnal_mean%path_rain_ice = undefined
 diurnal_mean%toa_lw_up      = undefined ; diurnal_mean%toa_lw_dn     = undefined
 diurnal_mean%toa_sw_up      = undefined ; diurnal_mean%toa_sw_dn     = undefined
 diurnal_mean%boa_lw_up      = undefined ; diurnal_mean%boa_lw_dn     = undefined
 diurnal_mean%boa_sw_up      = undefined ; diurnal_mean%boa_sw_dn     = undefined


 case('add') !-------------------------------------------------------------------------------------


 if(present(var2d) ) then

  if( masterproc ) print*,'MSG diurnal_daily_2d: add'


   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var2d)


   do j = 1,bnd(2) ; do i = 1,bnd(1)

    land_if: if( var2d(i,j)%iland == 1 ) then !over land

    latlon_if: if( var2d(i,j)%lat >= minlat .and. var2d(i,j)%lat <= maxlat .and.&
                   var2d(i,j)%lon >= minlon .and. var2d(i,j)%lon <= maxlon ) then

     !local solar time
     hour_loc = sdsu_gmt + var2d(i,j)%lon/15.
     if( hour_loc < -0.5 ) hour_loc = hour_loc + 24.
     if( hour_loc > 23.5 ) hour_loc = hour_loc -24.


      do ibin = 1, ibin_max
         if( hour_loc >  diurnal_lst(ibin) -  0.5  .and. &
             hour_loc <= diurnal_lst(ibin) +  0.5  ) then

             diurnal_cnt(ibin) = diurnal_cnt(ibin) + 1.e0  !update total count

             diurnal_sum(ibin)%rain_rate      = diurnal_sum(ibin)%rain_rate + var2d(i,j)%rain_rate
             diurnal_sum(ibin)%h2o_col        = diurnal_sum(ibin)%h2o_col   + var2d(i,j)%h2o_col
             diurnal_sum(ibin)%path_cloud_liq = diurnal_sum(ibin)%path_cloud_liq + var2d(i,j)%path_cloud_liq
             diurnal_sum(ibin)%path_rain_liq  = diurnal_sum(ibin)%path_rain_liq  + var2d(i,j)%path_rain_liq
             diurnal_sum(ibin)%path_cloud_ice = diurnal_sum(ibin)%path_cloud_ice + var2d(i,j)%path_cloud_ice
             diurnal_sum(ibin)%path_rain_ice  = diurnal_sum(ibin)%path_rain_ice  + var2d(i,j)%path_rain_ice
             diurnal_sum(ibin)%toa_lw_up      = diurnal_sum(ibin)%toa_lw_up + var2d(i,j)%toa_lw_up
             diurnal_sum(ibin)%toa_lw_dn      = diurnal_sum(ibin)%toa_lw_dn + var2d(i,j)%toa_lw_dn
             diurnal_sum(ibin)%toa_sw_up      = diurnal_sum(ibin)%toa_sw_up + var2d(i,j)%toa_sw_up
             diurnal_sum(ibin)%toa_sw_dn      = diurnal_sum(ibin)%toa_sw_dn + var2d(i,j)%toa_sw_dn
             diurnal_sum(ibin)%boa_lw_up      = diurnal_sum(ibin)%boa_lw_up + var2d(i,j)%boa_lw_up
             diurnal_sum(ibin)%boa_lw_dn      = diurnal_sum(ibin)%boa_lw_dn + var2d(i,j)%boa_lw_dn
             diurnal_sum(ibin)%boa_sw_up      = diurnal_sum(ibin)%boa_sw_up + var2d(i,j)%boa_sw_up
             diurnal_sum(ibin)%boa_sw_dn      = diurnal_sum(ibin)%boa_sw_dn + var2d(i,j)%boa_sw_dn

         endif
      enddo

    endif latlon_if
    endif land_if

   enddo ; enddo

 endif


 case('finish') !---------------------------------------------------------------------------------------

  if( masterproc ) print*,'MSG diurnal_daily_2d: finish'


#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', diurnal_cnt  )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%h2o_col   )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%path_cloud_liq )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%path_rain_liq )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%path_cloud_ice )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%path_rain_ice )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%toa_lw_up )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%toa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%toa_sw_up )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%toa_sw_dn )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%boa_lw_up )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%boa_lw_dn )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%boa_sw_up )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%boa_sw_dn )

  MPI_MASTERPROC: if(masterproc) then
#endif

!
! derive diurnal mean values
! 
 do ibin = 1, ibin_max

    if( diurnal_cnt(ibin) > 0. ) then

     diurnal_mean(ibin)%rain_rate      = diurnal_sum(ibin)%rain_rate / diurnal_cnt(ibin)
     diurnal_mean(ibin)%h2o_col        = diurnal_sum(ibin)%h2o_col   / diurnal_cnt(ibin)
     diurnal_mean(ibin)%path_cloud_liq = diurnal_sum(ibin)%path_cloud_liq / diurnal_cnt(ibin)
     diurnal_mean(ibin)%path_rain_liq  = diurnal_sum(ibin)%path_rain_liq  / diurnal_cnt(ibin)
     diurnal_mean(ibin)%path_cloud_ice = diurnal_sum(ibin)%path_cloud_ice / diurnal_cnt(ibin)
     diurnal_mean(ibin)%path_rain_ice  = diurnal_sum(ibin)%path_rain_ice  / diurnal_cnt(ibin)
     diurnal_mean(ibin)%toa_lw_up      = diurnal_sum(ibin)%toa_lw_up / diurnal_cnt(ibin)
     diurnal_mean(ibin)%toa_lw_dn      = diurnal_sum(ibin)%toa_lw_dn / diurnal_cnt(ibin)
     diurnal_mean(ibin)%toa_sw_up      = diurnal_sum(ibin)%toa_sw_up / diurnal_cnt(ibin)
     diurnal_mean(ibin)%toa_sw_dn      = diurnal_sum(ibin)%toa_sw_dn / diurnal_cnt(ibin)
     diurnal_mean(ibin)%boa_lw_up      = diurnal_sum(ibin)%boa_lw_up / diurnal_cnt(ibin)
     diurnal_mean(ibin)%boa_lw_dn      = diurnal_sum(ibin)%boa_lw_dn / diurnal_cnt(ibin)
     diurnal_mean(ibin)%boa_sw_up      = diurnal_sum(ibin)%boa_sw_up / diurnal_cnt(ibin)
     diurnal_mean(ibin)%boa_sw_dn      = diurnal_sum(ibin)%boa_sw_dn / diurnal_cnt(ibin)

    endif

 enddo

!
! write out cloud frequency diurnal profile
!
 write(yyyy,"(I4.4)") int( sdsu_yyyy )
 write(mm  ,"(I2.2)") int( sdsu_mm   )
 write(dd  ,"(I2.2)") int( sdsu_dd   )

 sdsu_io_file = trim(sdsu_dir_output)//'diurnal_mean_CRM2D_'//yyyy//'_'//mm//'_'//dd&
                //trim(output_suffix)//'.grads'
 print("(a)"),'MSG diurnal_mean2d: see output ->',trim(sdsu_io_file)

 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=ibin_max*4)

 irec = 1
 write(io,rec=irec) ( diurnal_mean(ibin)%rain_rate      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%h2o_col        , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%path_cloud_liq , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%path_rain_liq  , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%path_cloud_ice , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%path_rain_ice  , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%toa_lw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%toa_lw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%toa_sw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%toa_sw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%boa_lw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%boa_lw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%boa_sw_up      , ibin = 1, ibin_max ) ; irec = irec + 1
 write(io,rec=irec) ( diurnal_mean(ibin)%boa_sw_dn      , ibin = 1, ibin_max ) ; irec = irec + 1

 close(io)


!
! write grads control file
!


 sdsu_io_file = trim(sdsu_dir_output)//'diurnal_mean_CRM2D_'//yyyy//'_'//mm//'_'//dd&
                //trim(output_suffix)//'.ctl'

 print("(a)"),'MSG diurnal_mean2d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'diurnal_mean_CRM2D_'//yyyy//'_'//mm//'_'//dd&
                //trim(output_suffix)//'.grads'

 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   diurnal mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    1 linear 1. 1'
 write(io,*) 'ydef    ',ibin_max,' linear 0   1. '
 write(io,*) 'zdef    1 linear 0. 0.1'
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars ',irec-1
 write(io,*) 'rain    0  0  rain rate [mm/hr]'

 write(io,*) 'h2ocol  0  0  column water vaper [kg/m2]'
 write(io,*) 'clwp    0  0  cloud liquid water path [kg/m2]'
 write(io,*) 'lrwp    0  0  liquid rain water path [kg/m2]'
 write(io,*) 'ciwp    0  0  cloud ice water path [kg/m2]'
 write(io,*) 'irwp    0  0  ice rain water path [kg/m2]'

 write(io,*) 'lwut   0 0  TOA LW upwelling   flux [W/m2]'
 write(io,*) 'lwdt   0 0  TOA LW downwelling flux [W/m2]'
 write(io,*) 'swut   0 0  TOA SW upwelling flux   [W/m2]'
 write(io,*) 'swdt   0 0  TOA SW downwelling flux [W/m2]'
 write(io,*) 'lwus   0 0  Surface LW upwelling   flux [W/m2]'
 write(io,*) 'lwds   0 0  Surface LW downwelling flux [W/m2]'
 write(io,*) 'swus   0 0  Surface SW upwelling   flux [W/m2]'
 write(io,*) 'swds   0 0  Surface SW downwelling flux [W/m2]'

 write(io,*) 'endvars '
 close(io)


#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG diurnal_daily_2d: There is no such proc')
 end select


 return
 end subroutine diurnal_daily_2d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine diurnal_mean_3d(proc,var3d,press_stag,var2d)
 implicit none
!--------------------------------------------------------------------------------------------------
! Comments: Derive diurnal composite of 3D parameters.
!           
! History:
! 06/2016  Toshi Matsui@NASA GSFC ; Initial.
!
! References:
!-----------------------------------------------------------------------------------------------------
 character*(*) :: proc ! process
 type ( atmos_parameter ), dimension(:,:,:), optional :: var3d ! input surface parameter assumed 3d array
 real, dimension(:,:,:), optional :: press_stag ! input pressure sraggered array
 type ( surface_parameter ), dimension(:,:), optional :: var2d ! input surface parameter assumed 3d array

 integer :: bnd(3)  !upper memory bound for 2D array

!
! statisticsl parameters
!
 type( atmos_parameter ), allocatable, dimension(:,:),save :: diurnal_sum,diurnal_mean
 real(sdsu_fps), allocatable, save :: diurnal_lst(:), diurnal_cnt(:,:)
 integer,save :: ibin_max
 real(sdsu_fps) :: hour_loc

 integer :: i,j,k,p, ibin, ierr, irec
 integer,parameter :: io = 154

 integer,parameter :: plev = 20
 real(sdsu_fps),parameter :: pressure_levels(20) = (/1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,500.,&
                                                     450.,400.,350.,300.,250.,200.,150.,100.,50./)

 real(sdsu_fps) :: pr

 select case(trim(proc))
 case('init') !---------------------------------------------------------------------------------
 if(masterproc) print*,'MSG diurnal_mean_3d: initialize parameters'

!
! maximum bin of zone
!
  ibin_max  = 24

!
! allocate
!
 if( .not. allocated(diurnal_sum) ) then
     allocate( diurnal_sum (ibin_max,plev), &
               diurnal_mean(ibin_max,plev), &
               diurnal_cnt (ibin_max,plev), &
               diurnal_lst (ibin_max), &
               stat=ierr )
      if (ierr /= 0) call stop_sdsu(myrank,'MSG diurnal_mean_3d: allocation error -> Terminate program.')
 endif

!
! assign latitude bin
!
  do ibin = 1, ibin_max
     diurnal_lst(ibin) = real(ibin) -1.   !0~23LST (starting -0.5~0.5LST and ending 22.5~23.5LST
  enddo


!
! zero out other parameters
!
 diurnal_cnt = 0.

 diurnal_sum%rain_rate      = 0. ; diurnal_sum%sh            = 0.
 diurnal_sum%w              = 0. ; diurnal_sum%latent_heat   = 0.
 diurnal_sum%lw_heat        = 0. ; diurnal_sum%sw_heat       = 0.

 diurnal_mean%rain_rate      = undefined ; diurnal_mean%sh            = undefined
 diurnal_mean%w              = undefined ; diurnal_mean%latent_heat   = undefined
 diurnal_mean%lw_heat        = undefined ; diurnal_mean%sw_heat       = undefined


 case('add') !-------------------------------------------------------------------------------------


 if(present(var3d) ) then

  if(masterproc) print*,'MSG diurnal_mean_3d: add parameters'

   !
   ! find upper bound from assumped 3D array
   !
   bnd = UBOUND(var3d)


   do j = 1,bnd(2) ; do i = 1,bnd(1)


    land_if: if( var2d(i,j)%iland == 1 ) then !over land

    latlon_if: if( var2d(i,j)%lat >= minlat .and. var2d(i,j)%lat <= maxlat .and.&
                   var2d(i,j)%lon >= minlon .and. var2d(i,j)%lon <= maxlon ) then

      !
      ! estimate local solar hour
      !
      hour_loc = sdsu_gmt + var2d(i,j)%lon/15.

      if( hour_loc < -0.5 ) hour_loc = hour_loc + 24.
      if( hour_loc > 23.5 ) hour_loc = hour_loc -24.

      do ibin = 1, ibin_max

         if( hour_loc >  diurnal_lst(ibin) -  0.5  .and. &
             hour_loc <= diurnal_lst(ibin) +  0.5  ) then

            !
            ! interpolation level
            !
            do p = 1, plev
               pr = pressure_levels(p)
               do k = 1, bnd(3)
                  if( pr <= press_stag(i,j,k) .and. pr > press_stag(i,j,k+1) ) then
                     diurnal_cnt(ibin,p)           = diurnal_cnt(ibin,p) + 1.
                     diurnal_sum(ibin,p)%rain_rate = diurnal_sum(ibin,p)%rain_rate + var3d(i,j,k)%rain_rate
                     diurnal_sum(ibin,p)%sh        = diurnal_sum(ibin,p)%sh        + var3d(i,j,k)%sh
                     diurnal_sum(ibin,p)%w         = diurnal_sum(ibin,p)%w         + var3d(i,j,k)%w
                     diurnal_sum(ibin,p)%latent_heat = diurnal_sum(ibin,p)%latent_heat + var3d(i,j,k)%latent_heat
                     diurnal_sum(ibin,p)%lw_heat   = diurnal_sum(ibin,p)%lw_heat   + var3d(i,j,k)%lw_heat
                     diurnal_sum(ibin,p)%sw_heat   = diurnal_sum(ibin,p)%sw_heat   + var3d(i,j,k)%sw_heat
                  endif
               enddo

            enddo

         endif
      enddo

    endif latlon_if
    endif land_if

   enddo ; enddo

 endif

!call mpi_wait_for_master !this process hold slave processors to wait for master proc. 
!stop

 case('finish') !---------------------------------------------------------------------------------------


#if MPI > 0 
!
! MPI reduce to add all together
!
  call mpi_sdsu_reduce( 'SUM', diurnal_cnt  )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%rain_rate )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%sh        )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%w         )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%latent_heat)
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%lw_heat   )
  call mpi_sdsu_reduce( 'SUM', diurnal_sum%sw_heat   )

  MPI_MASTERPROC: if(masterproc) then
#endif

!
! derive diurnal mean values
! 
 do p = 1, plev ; do ibin = 1, ibin_max

    if( diurnal_cnt(ibin,p) > 0. ) then

     diurnal_mean(ibin,p)%rain_rate      = diurnal_sum(ibin,p)%rain_rate / diurnal_cnt(ibin,p)
     diurnal_mean(ibin,p)%sh             = diurnal_sum(ibin,p)%sh        / diurnal_cnt(ibin,p)
     diurnal_mean(ibin,p)%w              = diurnal_sum(ibin,p)%w         / diurnal_cnt(ibin,p)
     diurnal_mean(ibin,p)%latent_heat    = diurnal_sum(ibin,p)%latent_heat / diurnal_cnt(ibin,p)
     diurnal_mean(ibin,p)%lw_heat        = diurnal_sum(ibin,p)%lw_heat   / diurnal_cnt(ibin,p)
     diurnal_mean(ibin,p)%sw_heat        = diurnal_sum(ibin,p)%sw_heat   / diurnal_cnt(ibin,p)

    endif

 enddo ; enddo


!
! write out cloud frequency diurnal profile
!

 sdsu_io_file = trim(sdsu_dir_output)//'diurnal_mean_CRM3D.grads'
 open(io,file=trim(sdsu_io_file),access='DIRECT',status='replace',recl=ibin_max*plev*4)
 print("(a)"),'MSG diurnal_mean3d: see output ->',trim(sdsu_io_file)

 irec = 1
 write(io,rec=irec) ( ( diurnal_mean(ibin,p)%rain_rate   , ibin = 1, ibin_max ), p =1,plev) ; irec = irec + 1
 write(io,rec=irec) ( ( diurnal_mean(ibin,p)%sh          , ibin = 1, ibin_max ), p =1,plev) ; irec = irec + 1
 write(io,rec=irec) ( ( diurnal_mean(ibin,p)%w           , ibin = 1, ibin_max ), p =1,plev) ; irec = irec + 1
 write(io,rec=irec) ( ( diurnal_mean(ibin,p)%latent_heat , ibin = 1, ibin_max ), p =1,plev) ; irec = irec + 1
 write(io,rec=irec) ( ( diurnal_mean(ibin,p)%lw_heat     , ibin = 1, ibin_max ), p =1,plev) ; irec = irec + 1
 write(io,rec=irec) ( ( diurnal_mean(ibin,p)%sw_heat     , ibin = 1, ibin_max ), p =1,plev) ; irec = irec + 1

 close(io)


!
! write grads control file
!
 sdsu_io_file = trim(sdsu_dir_output)//'diurnal_mean_CRM3D'//trim(output_suffix)//'.ctl'
 print("(a)"),'MSG diurnal_mean2d: see output.ctl ->',trim(sdsu_io_file)
 open(io,file=trim(sdsu_io_file))

 sdsu_io_file = 'diurnal_mean_CRM3D'//trim(output_suffix)//'.grads'
 write(io,*) 'dset  ',trim(sdsu_io_file)
 write(io,*) 'title   diurnal mean statistics '
 write(io,*) 'undef   -999.             '
 write(io,*) 'xdef    1 linear 1. 1'
 write(io,*) 'ydef    ',ibin_max,' linear 0   1. '
 write(io,*) 'zdef    ',plev,' levels ', pressure_levels(1:plev)
 write(io,*) 'tdef    1 linear ',grads_time,' 1hr'
 write(io,*) 'vars   ',irec-1
 write(io,*) 'rain   ',plev,'  0  rain rate [mm/hr]'
 write(io,*) 'sh     ',plev,'  0  specific humidity [kg/kg]'
 write(io,*) 'w      ',plev,'  0  vertical velocity [m/s]'
 write(io,*) 'lhr    ',plev,'  0  latent heat release [K/day]'
 write(io,*) 'lwh    ',plev,'  0  long wave heating [K/day]'
 write(io,*) 'swh    ',plev,'  0  short wave heating [K/day]'
 write(io,*) 'endvars'
 close(io)


#if MPI > 0 
  endif MPI_MASTERPROC
#endif

 case default ; call stop_sdsu(myrank,'MSG diurnal_mean_3d: There is no such proc')
 end select


 return
 end subroutine diurnal_mean_3d

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_diag


