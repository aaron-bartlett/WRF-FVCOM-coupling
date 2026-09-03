program read_CO2_data2

   ! This program read the CO2 concentration data in binary format
   ! Convert it to netcdf files with timestamp
   ! Notice that LON0 =-180+dx/2 and LAT0=-90 are hardcoded in read_data2 
   ! dx and dy (grid space for dataset) is computed in this code, from the Num of longitude and latitude 
   ! points which are read from the header of the datafile. 
   !
   ! Author: Jossy P. Jacob (ASTG/CISTO/GSFC/NASA/SSAI)   
   ! 
   !Usage: compile: make -f Makefile_CO2_conc
   ! run: ./Read_CO2_conc.x filename indir sdt sodt sndays
   !
   ! filename: Name of the file (e624.co2.b20100302.e20100307)
   ! indir : name of the directory from rundir where the output data need to be created, and indir/binary exists
   ! dt : the time interval of input data in hours 
   ! odt : the time interval of output data in hours
   ! sndays : number of days of data in one output file  
   !
   !  ./Read_CO2_conc.x e624.co2.b20100302.e20100307 conc 1 6 5
   !
   use netcdf
   implicit none
   INCLUDE 'netcdf.inc'

   ! Input Arguments
   integer           :: IARGC
!   integer, parameter :: NLON = 288, NLAT = 181  
   character*256                :: file1,indir,infile1,outfile,ofile,outdir
   character*30                 :: sdate
   character*2                  :: sdt, sodt,sndays
   
   integer                      :: vardim, status
   integer                      :: ntime,k,i
   integer                      :: TVarId, SVarId
   character*100                :: VarUnit1, VarUnit2
   real                         :: maxdiff, dx, dy
   integer                      :: year, dt,bdate,edate,ndays,it,odt,ntime2
   integer                      :: NLON, NLAT,NLEV
 !  real ::  LAT(NLAT), LON(NLON)
   real, allocatable  			:: LEV(:),LAT(:),LON(:)
   real (kind=4), allocatable :: VAR1(:,:,:,:)
   real (kind=4), allocatable :: VAR2(:,:,:,:)
   real (kind=4), allocatable :: ps(:,:,:)
   real (kind=4), allocatable :: ps2(:,:,:)
   integer, allocatable :: TIME(:)
   integer, allocatable :: TIME2(:)
   character *19, allocatable :: TIMES(:)
   character *19, allocatable :: TIMES2(:)

   ! ********************************************************************************************
   call GETARG (1, file1)
   print *,trim(file1)
   call GETARG (2, indir)
   print *,trim(indir)
   ! call GETARG (2, syear)
   ! read(syear,*)year
   ! print *, year
    call GETARG (3, sdt)
    read(sdt,*)dt
    call GETARG (4, sodt)
    read(sodt,*)odt 
    call GETARG (5, sndays)
    read(sndays,*)ndays
    !!read(sdx,'(f5.2)')dx
    !print *, sdx, dx 
    !call GETARG (6, sdy)
    !read(sdy,'(f5.2)')dy 
    !print *, sdy, dy 

   read(file1(11:14),*) year
   read(file1(11:18),*) bdate
   read(file1(21:28),*) edate
!   dt=1 ! hrly files are input
!   odt=6 ! 6 hrly files for output
   !ndays = 5 ! 5days of data in 1 file
   print *,'input data interval=',dt, 'output data interval=',odt,  'ndays= ',ndays
   sdate = file1(10:28)
   !infile1 = 'CO2Conc/'//trim(file1)
   !infile1 = 'CASA/conc/'//trim(file1)
   infile1 = trim(indir)//'/binary/'//trim(file1)
   ofile = 'CO2conc_'//file1(11:28)
   outfile = trim(indir)//'/'//trim(ofile)//'.nc'
   ! outdir = 'CO2c/'
   outdir = indir
   !  print *, 'file1: ',trim(infile1)
   print *, 'OUT: ',trim(outfile)
   !
   ! For 2D Variables:
   !
   !	ntime=4
   call get_timestamp_ndays(dt,ndays,file1,time,times,ntime)
   call read_data2(infile1,VAR1,ps,LON,LAT,NLON,NLAT,dt,year,dx,dy,NLEV,LEV,time,ntime)
   print *, 'read 2D CO2'
   ! Extract the data for output
   ntime2=24/odt* ndays
   allocate (VAR2(1:NLON, 1:NLAT, 1: NLEV, 1:ntime2))
   allocate (ps2(1:NLON, 1:NLAT,  1:ntime2))
   allocate (TIME2(1:ntime2))
   allocate (TIMES2(1:ntime2))
   it = 0
!   do i = 1,ntime,2
!by Z. Tao: skip first 5 hour data because data starts at 1 UTC
   do i = odt,ntime,odt
      it = it+1
      VAR2(:,:,:,it)= VAR1(:,:,:,i)
      ps2(:,:,it)=ps(:,:,i)
      TIME2(it)=TIME(i)
      TIMES2(it) = TIMES(i)
   enddo
   ! Write the outfile
   call write_nc3(outdir,VAR2,ps2,LON,LAT,time2,times2,NLON,NLAT,NLEV,LEV,ntime2,sdate,ndays,odt)
contains
   !-----------------------
   subroutine check(status)
      !-----------------------
      use netcdf

      integer, intent (in) :: status

      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if

   end subroutine check
   !---------------------
   subroutine read_data2(infile,VAR,ps,LON,LAT,NLON,NLAT,dt,year,dx,dy,NLEV,LEV,time,ntime)
      character*256,   intent(in)      ::infile 
      integer, intent(in)             :: year, dt
      integer, intent(out)            :: NLON,NLAT,NLEV
      real, intent(out)               :: dx, dy
      integer, intent(in)            :: ntime
!      real, intent(out):: LON(NLON)
!      real, intent(out):: LAT(NLAT)
      integer  :: time(ntime)

      real (kind=4), allocatable :: temp(:,:,:)
      real (kind=4), allocatable :: temp2(:,:)
      character*12                :: varname
      integer                     :: ncid,  i,j,m, n, itime,nt,ic
      integer                     :: ieof,recpos
      real 			:: LON0, LAT0
      real (kind=4)               :: head(337)
      real (kind=4),allocatable		:: ps(:,:,:) 
      real (kind=4), allocatable 		:: var(:,:,:,:) 
      real, allocatable  			:: LEV(:),LAT(:),LON(:)
            !    print *, trim(infile)
      open(unit=11, file=trim(infile),access='sequential',convert='big_endian',status='old', form='unformatted')
      read(11) head
      close(11)
      NLON=head(9)
      NLAT=head(10)
      NLEV=head(11)
      ic=head(40)

      allocate (ps(1:NLON,1:NLAT,1:ntime))
      allocate (temp2(1:NLON,1:NLAT)) 
      allocate (LEV(1:nlev))
      allocate (LON(1:NLON))
      allocate (LAT(1:NLAT))
      allocate (VAR(1:NLON,1:NLAT,1:NLEV,1:ntime))
      allocate (temp(1:NLON,1:NLAT,1:NLEV))
      !allocate (time(1:ntime))
      dx = 360.0/NLON
      dy = 180/(NLAT-1)
      LON0 = -180.+dx/2.0
      LAT0 = -90.
      print *, 'dx = ', dx, 'dy = ',dy , 'NLON=' , NLON, 'NLAT=',NLAT,'NLEV=',NLEV
      print *, 'LON0=',LON0,'LAT0=', LAT0

      do i=1,NLEV
         LEV(i)=i
      enddo
      do i=1,NLON
         LON(i)= LON0+ dx * float(i-1)
      enddo
      do j=1,NLAT
         LAT(j)= LAT0+ dy * float(j-1)
      enddo 
      print *, 'LON: ',LON
      print *, 'LAT: ',LAT
      print *, 'NTIME=',NTIME
      open(unit=11, file=trim(infile),access='sequential',convert='big_endian',status='old', form='unformatted')
      ieof =0
      print *, infile,ieof
      do nt = 1,ntime
         read(11) head,temp2
         read(unit=11,iostat=ieof) temp
         if (ieof == 0) then 
            VAR(:,:,:,nt)= temp
            ps(:,:,nt)=temp2
         else
            !                print *, 'File reading complete, ntime=',nt
            return
         endif
      enddo ! end of ntime
      close(11)
      print *, 'Finished reading data'
   end subroutine read_data2
   !---------------------

   subroutine get_timestamp_ndays(dt,ndays,file1,time,stime,ntime)
      integer                     :: ntime, dt,ndays
      integer                     :: nt , itime,nd
      integer                     :: year, mon, day, day1, dd,hour,hr, isleapyear
      integer, dimension(12)      :: mdd
      integer, intent(out), allocatable        :: time(:) 
      character *19, intent(out), allocatable        :: stime(:) 
      character *256 :: file1
      ! stime format 1999-05-06_12:00:00  
      character :: syear*4, smon*2 , sday*2, shr*2 
      read(file1(11:14),*) year
      read(file1(15:16),*) mon
      read(file1(17:18),*) day
      print *, year, mon, day
      print *,(year/4)*4
!      if (((year/4)*4) == year) then 
      if(mod(year,400) == 0 .or. (mod(year,4) == 0 .and. (mod(year,100) .ne. 0))) then		!by ZTao
         isleapyear = 1
         mdd = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
         print *, year, ' leapyear'
      else 
         isleapyear = 0
         mdd = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
         print *, year, ' not leapyear'
      endif
      ntime = 24/dt*ndays ! for 5 days 
      allocate (time(ntime))
      allocate (stime(ntime))
      ! Get timestamp for 5days frim bdate to edate 20100302 to 20100307 
      itime=0
      dd=day
      do nd = 1, ndays
         do hr = 1, 23,dt		!first record of input file is at 1 UTC	by Z. Tao
            itime=itime+1 
            hour = hr
            time(itime) = year * 10**6 + mon *10**4 + dd * 10**2 + hour
            write(syear,'(I4)') year
            write(smon,'(I2.2)') mon
            write(sday,'(I2.2)') dd
            write(shr,'(I2.2)') hour
            stime(itime) = syear//'-'//smon//'-'//sday//'_'//shr//':00:00'
            !time yearmmddhh 2010030100
            !stime  2010-03-01_00:00:00
         enddo ! end of hour
         dd= dd+1 
         if (dd > mdd(mon)) then 
            dd =1 
            mon = mon + 1
! by Z.Tao	    
	    if(mon > 12) then
	      mon = 1
	      year = year + 1
	    end if
         endif

!by Z. Tao for every 24th record that is 0 UTC of next day
	itime = itime + 1
        time(itime) = year * 10**6 + mon *10**4 + dd * 10**2
        write(syear,'(I4)') year
        write(smon,'(I2.2)') mon
        write(sday,'(I2.2)') dd
        stime(itime) = syear//'-'//smon//'-'//sday//'_00:00:00'
      enddo ! end ndays
      !              print *,'file1 = ',trim(file1)
                    print *,'ntime = ',ntime
                   print *,'time = ',time 
   end subroutine get_timestamp_ndays

   !..........................................................................................

   !............................................................................................
   subroutine write_nc3(outdir,TVAR,pvar,LON,LAT,JTIME,STIME,NLON,NLAT,nlev,Lev,ntime,sdate,ndays,odt)
      !............................................................................................
      use netcdf
      implicit none
      !    INCLUDE 'netcdf.inc'
      integer, intent(in)          :: NLON,NLAT,ntime,nlev,ndays
      character*256                :: outdir,fname
      character *30                :: sdate
      real, intent(in)             :: TVAR(nlon,nlat,nlev,ntime)
      real, intent(in)             :: pvar(nlon,nlat,ntime)
      real, intent(in)             :: LON(nlon)
      real, intent(in)             :: LAT(nlat)
      real, intent(in)             :: LEV(nlev)
      integer, intent(in)          :: JTIME(ntime)
      character *19, intent(in)          :: STIME(ntime)
      real              :: pvar_out(nlon,nlat)
      real              :: tvar_out(nlon,nlat,nlev)

      integer                     :: ncid, status
      integer                     :: TVarId,LONvarid,LATvarid,TIMEvarid,LevVarid,PVarID,dim2id,TIME2varid
      integer                     :: m,n , nt, nd,odt
      integer                     :: LatDimID,LonDimID,LevdimID,TimeDimID
      integer , dimension (3) :: count1 , start1 
      integer , dimension (4) :: count2 , start2 
      character (len = *), parameter :: LAT_NAME = "lat"
      character (len = *), parameter :: LON_NAME = "lon"
      character (len = *), parameter :: LEV_NAME = "lev"
      character (len = *), parameter :: TIME_NAME = "time"
      character (len = *), parameter :: STIME_NAME = "Times"
      character (len =*), parameter :: TVAR_NAME = "CO2c" 
      character (len =*), parameter :: PVAR_NAME = "PS" 

      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: TVAR_UNITS = "ppmv"
      character (len = *), parameter :: PVAR_UNITS = "mb"
      character (len = *), parameter :: LAT_UNITS = "degrees_north"
      character (len = *), parameter :: LON_UNITS = "degrees_east"
      character (len = *), parameter :: LEV_UNITS = "Levels"
      character (len = *), parameter :: TIME_UNITS = "YYYYMMDDHH"
      character (len = *), parameter :: TIME_UNITS2 = "YYYY-MM-DD_HH:MI:SS"

      real :: TVAR1(nlon,nlat,nlev,1)
      real :: PVAR1(nlon,nlat,1)
      character *19 ::stime1 
      ! Write each data in one data file 6hourly 
      do nt = 1,ntime
         fname = trim(outdir)//'/CASACO2.'//stime(nt)//'.nc'  
         TVAR1(:,:,:,1) = TVAR(:,:,:,nt)
         PVAR1(:,:,1) = pvar(:,:,nt)
         stime1 = stime(nt)
         !	print *, fname, stime1

         !	print *,'size of TVAR=',size(TVAR)
         !	print *, 'Starting to write the data'
         ! Open the file to write: 
         call check(nf90_create(trim(fname),nf90_clobber,ncid))

         !	print *, 'Opened file to write the data'
         ! Define the dimensions: 
         call check(nf90_def_dim(ncid, LAT_NAME,NLAT,LatDimID))
         call check(nf90_def_dim(ncid, LON_NAME,NLON,LonDimID))
         call check(nf90_def_dim(ncid, LEV_NAME,NLEV,LevDimID))
         call check(nf90_def_dim(ncid, TIME_NAME,1,TimeDimID))
         call check(nf90_def_dim(ncid, 'date_len',19,dim2id))

         !	print *, 'Defined Dimensions'

         !Define the coordinate variables 
         call check(nf90_def_var(ncid, LAT_NAME,NF90_REAL, LatDimID,Latvarid))
         call check(nf90_def_var(ncid, LON_NAME,NF90_REAL,LonDimID, Lonvarid))
         call check(nf90_def_var(ncid, LEV_NAME,NF90_INT,LevDimID,Levvarid))
         ! call check(nf90_def_var(ncid, TIME_NAME,NF90_INT,TimeDimID,Timevarid))
         call check(nf90_def_var(ncid, STIME_NAME,NF90_CHAR,(/ dim2id,TimeDimID /),Time2varid))
         call check(nf90_def_var(ncid,TVAR_NAME,NF90_REAL,(/ LonDimId, LatDimID, LevDimID, TimeDimID /), TVarId))
         call check(nf90_def_var(ncid,PVAR_NAME,NF90_REAL,(/ LonDimId, LatDimID, TimeDimID /), PVarId))
         !	print *, 'Defined Coordinate variables'

         call check(nf90_put_att(ncid, nf90_global, "Data_Interval", sdate))
         ! Assign units attributes to coordinate variables.
         call check( nf90_put_att(ncid, latvarid, UNITS, LAT_UNITS) )
         call check( nf90_put_att(ncid, lonvarid, UNITS, LON_UNITS) )
         call check( nf90_put_att(ncid, levvarid, UNITS, LEV_UNITS) )
         call check( nf90_put_att(ncid, Time2varid, UNITS, TIME_UNITS2) )

         !    call check( nf90_put_att(ncid, Timevarid, UNITS, TIME_UNITS) )
         ! Assign units attributes to the netCDF variables.
         call check( nf90_put_att(ncid, TVarid, UNITS, TVAR_UNITS) )
         call check( nf90_put_att(ncid, PVarid, UNITS, PVAR_UNITS) )
         ! End define mode.
         call check( nf90_enddef(ncid) )

         ! Write the coordinate variable data. This will put the latitudes
         ! and longitudes of our data grid into the netCDF file.
         call check( nf90_put_var(ncid, LATvarid, LAT) )
         !	print *, 'Wrote LAT'
         call check( nf90_put_var(ncid, LONvarid, LON) )
         call check( nf90_put_var(ncid, LEVvarid, LEV) )
         !	print *, 'Wrote LON'
         !  call check( nf90_put_var(ncid, TIMEvarid, JTIME) )
         call check( nf90_put_var(ncid, TIME2varid, STIME1) )
         !	print *, 'Wrote TIME'
         ! count1 = (/ NLON, NLAT, 1 /)
         ! start1 = (/ 1, 1, 1 /)
         ! count2 = (/ NLON, NLAT, nlev, 1 /)
         ! start2 = (/ 1, 1, 1, 1 /)
         !  
         !    ! Write the pretend data. This will write our surface pressure and
         !    ! surface temperature data. The arrays only hold one timestep worth
         !    ! of data. We will just rewrite the same data for each timestep. In
         !    ! a real :: application, the data would change between timesteps.
         !    do nt = 1, ntime
         !      PVAR_out = PVAR(:,:,nt)
         !      TVAR_out = TVAR(:,:,:,nt)
         !      start2(4) = nt
         !      start1(3) = nt
         !       call check( nf90_put_var(ncid, TVarid, TVAR_out, start = start2, &
         !                                count = count2) )
         !       call check( nf90_put_var(ncid, PVarid, PVAR_out, start = start1, &
         !                                count = count1) )
         !    end do
         call check( nf90_put_var(ncid, TVarid, TVAR1))
         call check( nf90_put_var(ncid, PVarid, PVAR1))

         ! Write the variables
         !    call check(nf90_put_var(ncid,Tvarid, TVAR))
         !	print *, 'Wrote TVAR'
         call check(nf90_close(ncid))
      enddo ! end ntime	
   end subroutine write_nc3

END PROGRAM read_CO2_data2
