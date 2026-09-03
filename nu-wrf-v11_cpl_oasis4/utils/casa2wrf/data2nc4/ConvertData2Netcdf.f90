program ConvertData2Netcdf
!=============================================================
! This program read the datasets for CASACO2 pre-processor and  
! Convert it to netcdf files. 
!=============================================================
! Author: Jossy P. Jacob (ASTG/CISTO/GSFC/NASA/SSAI)   
! April 2015
!=============================================================
   use netcdf
   implicit none
   INCLUDE 'netcdf.inc'

   ! Input Arguments
   integer           :: IARGC
   character*256                :: file1,indir,outdir,infile1,infile2,outfile,ofile,casadir
   character*30                 :: sdate, VarName
   character                    :: syear*4, smon*2, sday*2, shr*2
   character*10                  :: sdate1

   integer                      :: vardim, status
   integer                      :: ntime, NLON,NLAT,k,i,j,ntime1,ntime2,it, ndays,nmon
   integer                      :: TVarId, SVarId
   character*100                :: VarUnit1, VarUnit2
   real                         :: maxdiff,dt1
   integer                      :: year, mon, day, hour,  dt,ntime3,enddate,startdate

   real, allocatable ::  LAT(:), LON(:)
   real, allocatable :: VAR1(:,:,:), FossilFuel(:,:,:),OceanCO2(:,:,:)
   real, allocatable ::  LAT1(:), LON1(:)
   real, allocatable :: VAR2(:,:,:), WildFire(:,:,:), BioFuel(:,:,:), RespCO2(:,:,:), NPP(:,:,:)


   integer, allocatable :: TIME(:),TIME2(:)
   integer :: mdd(12), istart,iend
   character(len=19), allocatable :: STIME(:),TIMES2(:)

   ! ********************************************************************************************
  ! call GETARG (1, file1)
   !print *,trim(file1)
   !call GETARG (2, indir)
   !print *,trim(indir)
   ! call GETARG (2, syear)
   ! read(syear,*)year
   ! print *, year
   ! call GETARG (3, sdt)
   ! read(sdt,*)dt
  ! read(file1(29:32),*) year
  ! read(file1(15:15),*) dt
  ! read(file1(15:15),*) dt1
  ! read(file1(21:23),*) NLON
  ! read(file1(25:27),*) NLAT
  ! print *,year,dt,NLON, NLAT
  ! sdate = file1(15:15)
   !!casadir = 'CASA/flux/'
   !casadir = indir 
!==========Fossil Fuel Data =======================
! D_float.2010 (direct access binary)
! FossilFuel_2010.dat (binary)
!==================================================
   indir = 'Binary/'
   outdir = 'Netcdf_data/'
   infile1 = trim(indir)//'FossilFuel_2010.dat'
  ! NLON = 288
  ! NLAT = 181
  ! ntime = 1 
   year = 2010
   VarName = 'FossilFuel' 
   !infile1 = trim(casadir)//'/binary/'//trim(file1)
   !
   ! For 2D Variables:
   !
    print *, 'Read datafile: ',infile1
 !  allocate(VAR1(1:NLON, 1:NLAT, 1:ntime)) 
   call read_data(infile1,VAR1,NLON,NLAT, LON, LAT,ntime)
   allocate(FossilFuel(1:NLON, 1:NLAT, 1:ntime)) 
   allocate(TIME(1:ntime))
   allocate(stime(1:ntime)) 
   FossilFuel = VAR1
   outfile = trim(outdir)//'FFUE_2010_year.nc'
   print *, 'OUT: ',trim(outfile)
   time = 2010010100
   dt1 =  0
   sdate = 'Annual Mean'
   write(syear,'(I4)') year
   stime(1) = syear//'-01-01_00:00:00'
   call write_nc3(outfile,FossilFuel,LON,LAT,time,stime,NLON,NLAT,ntime,sdate,dt1, VarName)
   print *, 'Write datafile: ',outfile
   deallocate (VAR1, Time, stime)
   deallocate (LAT, LON)
!
!=================OCEAN CO2=======================
!     Read the Ocean CO2 emissions (taka* files)
!   infile2 ='Taka02.1x1.25.data' 
!=================================================
   infile2 = trim(indir)//'OceanCO2_2010.dat' 
!   NLON = 288
!   NLAT = 181
!   ntime = 12  
   VarName = 'OceanCO2' 
   print *, 'Read datafile: ',infile2
  ! allocate(VAR1(1:NLON, 1:NLAT, 1:ntime)) 
   call read_data(infile2,VAR1,NLON,NLAT, LON, LAT,ntime)
   allocate(OceanCO2(1:NLON, 1:NLAT, 1:ntime)) 
   allocate(TIME(1:ntime))
   allocate(stime(1:ntime)) 
   OceanCO2 = VAR1
   outfile = trim(outdir)//'OCO2_2010_mon0.nc'
   print *, 'OUT: ',trim(outfile)
  
   do i = 1,ntime
      mon = i
      day = 1
      hour = 0
      time(i) = year * 10**6 + mon *10**4 + day * 10**2 + hour
      write(syear,'(I4)') year
      write(smon,'(I2.2)') mon
      write(sday,'(I2.2)') day
      write(shr,'(I2.2)') hour
      stime(i) = syear//'-'//smon//'-'//sday//'_'//shr//':00:00'
   enddo
   dt1 =  1
   sdate = 'Monthly mean'
   call write_nc3(outfile,OceanCO2,LON,LAT,time,stime,NLON,NLAT,ntime,sdate,dt1,VarName)
   print *, 'Write datafile: ',outfile
   deallocate (VAR1, Time, stime)
   deallocate (LAT, LON)
!=================NPP DATA=======================
! NPP2010.txt 
!NPP.dat (binary)
!================================================
   infile1 =trim(indir)//'NPP_2010.dat' 
   VarName = 'NPP' 
   !NLON = 720
   !NLAT = 360
   !ntime = 12
   print *, 'Read datafile: ',infile1
   call read_data(infile1,VAR2,NLON,NLAT, LON1, LAT1,ntime)
   allocate(NPP(1:NLON, 1:NLAT, 1:ntime)) 
   allocate(TIME(1:ntime))
   allocate(stime(1:ntime)) 
   outfile = trim(outdir)//'NPP0_2010_mon0.nc'
   NPP = VAR2
   print *, 'NPP ', NLON, NLAT, ntime
   print *, size(NPP, 1),  size(NPP, 2),  size(NPP, 3) 
   do i = 1,ntime
      mon = i
      day = 1
      hour = 0
      time(i) = year * 10**6 + mon *10**4 + day * 10**2 + hour
      write(syear,'(I4)') year
      write(smon,'(I2.2)') mon
      write(sday,'(I2.2)') day
      write(shr,'(I2.2)') hour
      stime(i) = syear//'-'//smon//'-'//sday//'_'//shr//':00:00'
   enddo
   dt1 =  1
   sdate = 'Monthly mean'
   call write_nc3(outfile,NPP,LON1,LAT1,time,stime,NLON,NLAT,ntime,sdate,dt1,VarName)
   print *, 'Wrote datafile: ',outfile
   deallocate (VAR2, Time, stime)
   deallocate (LAT1, LON1)

!=================Bio Fuel Emissions data========
! FUE2010.txt 
! BFUEL_2010.dat (binary)
!================================================
   infile1 = trim(indir)//'BFUEL_2010.dat' 
   VarName = 'BioFuel' 
   !NLON = 720
   !NLAT = 360
   !ntime = 12

   print *, 'Read datafile: ',infile1
   call read_data(infile1,VAR2,NLON,NLAT, LON1, LAT1,ntime)
   allocate(BioFuel(1:NLON, 1:NLAT, 1:ntime)) 
   allocate(TIME(1:ntime))
   allocate(stime(1:ntime)) 
   outfile = trim(outdir)//'BFUE_2010_mon0.nc'
   BioFuel = VAR2
   do i = 1,ntime
      mon = i
      day = 1
      hour = 0
      time(i) = year * 10**6 + mon *10**4 + day * 10**2 + hour
      write(syear,'(I4)') year
      write(smon,'(I2.2)') mon
      write(sday,'(I2.2)') day
      write(shr,'(I2.2)') hour
      stime(i) = syear//'-'//smon//'-'//sday//'_'//shr//':00:00'
   enddo
   dt1 =  1
   sdate = 'Monthly Climatology'
   call write_nc3(outfile,BioFuel,LON1,LAT1,time,stime,NLON,NLAT,ntime,sdate,dt1,VarName)
   print *, 'Wrote datafile: ',outfile
   deallocate (VAR2, Time, stime)
   deallocate (LAT1, LON1)


!=================Resp Emissions data========
! resp2010.txt 
! RESP_2010.dat (binary)
!================================================
   infile1 = trim(indir)//'RESP_2010.dat' 
   VarName = 'RespCO2' 
   !NLON = 720
   !NLAT = 360
   !ntime = 12

   print *, 'Read datafile: ',infile1
   call read_data(infile1,VAR2,NLON,NLAT, LON1, LAT1,ntime)
   allocate(RespCO2(1:NLON, 1:NLAT, 1:ntime)) 
   allocate(TIME(1:ntime))
   allocate(stime(1:ntime)) 
   outfile = trim(outdir)//'RESP_2010_mon0.nc'
   RespCO2 = VAR2
   do i = 1,ntime
      mon = i
      day = 1
      hour = 0
      time(i) = year * 10**6 + mon *10**4 + day * 10**2 + hour
      write(syear,'(I4)') year
      write(smon,'(I2.2)') mon
      write(sday,'(I2.2)') day
      write(shr,'(I2.2)') hour
      stime(i) = syear//'-'//smon//'-'//sday//'_'//shr//':00:00'
   enddo
   dt1 =  1
   sdate = 'Monthly Climatology'
   call write_nc3(outfile,RespCO2,LON1,LAT1,time,stime,NLON,NLAT,ntime,sdate,dt1,VarName)
   print *, 'Wrote datafile: ',outfile
   deallocate (VAR2, Time, stime)
   deallocate (LAT1, LON1)


!=================Wild Fire Emissions daily data========
! FIRE2010.txt 
! FIRE2010.dat (binary)
!================================================
   infile1 =trim(indir)//'FIRE_2010.dat' 
   VarName = 'WildFire' 
   !NLON = 720
   !NLAT = 360
   !ntime = 365

   print *, 'Read datafile: ',infile1
   call read_data(infile1,VAR2,NLON,NLAT, LON1, LAT1,ntime)
   allocate(WildFire(1:NLON, 1:NLAT, 1:ntime)) 
   !allocate(TIME(1:ntime))
   !allocate(stime(1:ntime)) 
   outfile = trim(outdir)//'FIRE_2010_daily_01.nc'
   WildFire = VAR2
   dt = 1
   year = 2010
   call get_timestamp(dt,time,year,ntime,stime,mdd)
   print *, 'Size of time, stime are:',size(time), size(stime)
   sdate = 'Daily Emissions'
   if (ntime >= 365) then 
      nmon = 12
      istart = 1
      ntime1 = 0
      do mon = 1, nmon
          write(smon,'(I2.2)') mon
          outfile = trim(outdir)//'FIRE_2010_daily_'//smon//'.nc'
          print *, 'outfile name: ',outfile, 'mdd(mon)',mdd(mon)
          allocate(VAR1(1:NLON, 1:NLAT,1:mdd(mon)))
          istart = istart + ntime1 
          iend = istart + mdd(mon) - 1
          VAR1 = WildFire(:,:,istart:iend);
          ntime1 = iend - istart +1
          print *, 'mysize',ntime1, size(VAR1), 'istart/iend= ',istart, iend
          call write_nc3(outfile,VAR1,LON1,LAT1,time(istart:iend),stime(istart:iend),NLON,NLAT,ntime1,sdate,dt1,VarName)
          print *, 'Wrote datafile: ',outfile
          deallocate(VAR1)
      enddo
   endif

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
   subroutine read_data(infile,VAR,NLON,NLAT, LON, LAT,ntime)
      character*256,   intent(in)      ::infile 
      integer, intent(out)             :: NLON,NLAT
      integer, intent(out)            :: ntime
      real, intent(out), allocatable  :: VAR(:,:,:)
      real, intent(out), allocatable :: LON(:)
      real, intent(out), allocatable :: LAT(:) 
    
      real (kind =4), allocatable :: temp(:,:)
      character*12                :: varname
      integer                     :: ncid,  i,j,m, n, itime,nt
      integer                     :: ieof,recpos
      real 			:: dx, dy, LON0, LAT0

      !NLON = 288
      !NLAT = 181
      ieof =0
      open(unit=11, file=trim(infile), form='unformatted',convert='big_endian')
      read(11,iostat=ieof) NLON, NLAT, ntime
      print *,'myfile', infile
      print *, 'LON: ',NLON
      print *, 'LAT: ',NLAT
      print *, 'NTIME=',ntime
      dx = 360.0/(NLON)
      dy = 180./(NLAT-1)
      LAT0= -90. 
      !LON0 = -180.
      if (NLON == 288 ) then 
         LON0 = 0.
         !LON0 = LON0+ dx/2
         dy = 180./(NLAT-1)
         LAT0 = -90.0 
      elseif (NLON == 720) then 
         LON0 = -180.
         LON0 = LON0 + dx/2
         dy = 180./(NLAT)
         LAT0=-89.75
      endif
      !call get_timestamp(dt,time,year,ntime,stime)
      allocate (LON(1:NLON))
      allocate (LAT(1:NLAT))
      allocate (VAR(1:NLON,1:NLAT,1:ntime))
      allocate (temp(1:NLON,1:NLAT))
      do i=1,NLON
         LON(i)= LON0+ dx * float(i-1)
      enddo
      do j=1,NLAT
         LAT(j)= LAT0+ dy * float(j-1)
      enddo
      print *, 'LON: ',LON
      print *, 'LAT: ',LAT
      !print *, 'NTIME=',ntime
     ! open(unit=11, file=trim(infile), form='unformatted',convert='big_endian')
      !ntime = 24/dt*365 ! 1 year data in 1 fie
     ! recpos = 0 
      ieof =0
      do nt = 1, ntime
        recpos = recpos + 1 
        read(unit=11,iostat=ieof) temp
        if (ieof == 0) then 
           VAR(:,:,nt)= temp
        !print *, 'size(tmp)=',size(temp), 'ntime =', nt
        !print *, temp(100,:)
        !print *, 'max, min for the variable: ',maxval(temp), minval(temp)
        else
           print *, 'File reading complete for ',infile, ',ntime = ', nt 
           exit
        endif
      enddo
   end subroutine read_data
!
   !---------------------
   subroutine read_data_DirectAccess(infile,VAR,NLON,NLAT, LON, LAT,ntime)
      character*256,   intent(in)      ::infile 
      integer, intent(in)             :: NLON,NLAT
      integer, intent(in)            :: ntime
      real, intent(out), allocatable  :: VAR(:,:,:)
      real, intent(out), allocatable :: LON(:)
      real, intent(out), allocatable :: LAT(:) 
    
      real (kind =4), allocatable :: temp(:,:)
      character*12                :: varname
      integer                     :: ncid,  i,j,m, n, itime,nt
      integer                     :: ieof,recpos
      real 			:: dx, dy, LON0, LAT0

      !NLON = 288
      !NLAT = 181
      dx = 360.0/(NLON)
      dy = 180./(NLAT-1)
      !LON0 = -180.
      LON0 = 0.
      LON0 = LON0+ dx/2
      LAT0= -90. 

      !call get_timestamp(dt,time,year,ntime,stime)
      allocate (LON(1:NLON))
      allocate (LAT(1:NLAT))
      allocate (VAR(1:NLON,1:NLAT,1:ntime))
      allocate (temp(1:NLON,1:NLAT))
      do i=1,NLON
         LON(i)= LON0+ dx * float(i-1)
      enddo
      do j=1,NLAT
         LAT(j)= LAT0+ dy * float(j-1)
      enddo
      print *, 'LON: ',LON
      print *, 'LAT: ',LAT
      print *, 'NTIME=',ntime
      open(unit=11, file=trim(infile), form='unformatted',convert='big_endian',access='direct',recl=4*NLON*NLAT)
      !ntime = 24/dt*365 ! 1 year data in 1 fie
      recpos = 0 
      ieof =0
      do nt = 1, ntime
        recpos = recpos + 1 
        read(unit=11,iostat=ieof,rec=recpos) temp
        !if (ieof == 0) then 
           VAR(:,:,nt)= temp
        print *, 'size(tmp)=',size(temp), 'ntime =', nt
        print *, temp(100,:)
        print *, 'max, min for the variable: ',maxval(temp), minval(temp)
        !else
        !   print *, 'File reading complete for ',infile, ',ntime = ', nt 
        !   exit
        !endif
      enddo
   end subroutine read_data_DirectAccess

   !---------------------
 
   !---------------------
   subroutine get_timestamp(dt,time,year,ntime,stime,mdd)
      integer                     :: ntime, dt
      integer                     :: nt , itime
      integer 			:: year, mon, day, hour,hr, isleapyear
      integer, dimension(12)      :: mdd
      integer, intent(out), allocatable        :: time(:)
      character *19, intent(out), allocatable  :: stime(:)
      ! stime format 1999-05-06_12:00:00  
      character :: syear*4, smon*2 , sday*2, shr*2 
      print *,(year/4)*4
!      if (((year/4)*4) == year) then 
      if(mod(year,400) == 0 .or. (mod(year,4) == 0 .and. (mod(year,100) .ne. 0))) then          
         isleapyear = 1
      else 
         isleapyear = 0
      endif
      if (isleapyear == 1) then
         mdd = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
         print *, year, ' leapyear'
         !ntime = 24/dt * 366
      elseif (isleapyear == 0) then 
         mdd = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
         !ntime = 24/dt * 365
         print *, year, ' not leapyear'
         print *, 'ntime', ntime
      endif
      allocate (time(ntime))
      allocate (stime(ntime))
      ! Get timestamp for 1 year
      itime=0
      do mon = 1,12 
         do day = 1, mdd(mon) 
           ! do hr = 1, 24,dt
               hour = 0
               itime=itime+1 
               !time yearmmddhh
               time(itime) = year * 10**6 + mon *10**4 + day * 10**2 + hour
               write(syear,'(I4)') year
               write(smon,'(I2.2)') mon
               write(sday,'(I2.2)') day
               write(shr,'(I2.2)') hour
               stime(itime) = syear//'-'//smon//'-'//sday//'_'//shr//':00:00'
           ! enddo
         enddo
      enddo
      !              print *,'ntime = ',ntime
       !print *,'time -new = ',time
       !print *,'stime-new = ',stime
   end subroutine get_timestamp

   !............................................................................................
   subroutine write_nc3(fname_out,TVAR,LON,LAT,TIME,STIME,NLON,NLAT,ntime,sdate,dt,VarName)
      !............................................................................................
      use netcdf
      implicit none
      !    INCLUDE 'netcdf.inc'
      integer, intent(in)          :: NLON,NLAT,ntime
      character*256                :: fname_out
      character *30                :: sdate, VarName
      real, intent(in)             :: TVAR(nlon,nlat,ntime)
      real, intent(in)             :: LON(nlon)
      real, intent(in)             :: LAT(nlat)
      real, intent(in)             :: dt
      integer, intent(in)          :: TIME(ntime)
      character(len=19), intent(in)          :: STIME(ntime)

      integer                     :: ncid, status
      integer                     :: TVarId,LONvarid,LATvarid,TIMEvarid,TIME2varid
      integer                     :: m,n   
      integer                     :: LatDimID,LonDimID,TimeDimID,dim1id,dim2id
      character (len = *), parameter :: LAT_NAME = "lat"
      character (len = *), parameter :: LON_NAME = "lon"
      character (len = *), parameter :: TIME_NAME = "time"
      character (len = *), parameter :: STIME_NAME = "Times"
      character (len =30) :: TVAR_NAME =  "VAR"

      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: TVAR_UNITS = "ppmv"
      character (len = *), parameter :: LAT_UNITS = "degrees_north"
      character (len = *), parameter :: LON_UNITS = "degrees_east"
      !character (len = *), parameter :: TIME_UNITS = "YYYYMMDDHH"
      character (len = *), parameter :: TIME_UNITS2 = "YYYY-MM-DD_HH:MI:SS"
      TVAR_NAME = VarName
      print *, 'Starting to write the data'
      !  print *, 'STIME', STIME
      ! Open the file to write: 
      call check(nf90_create(trim(fname_out),nf90_clobber,ncid))

      !	print *, 'Opened file to write the data'
      ! Define the dimensions: 
      call check(nf90_def_dim(ncid, LAT_NAME,NLAT,LatDimID))
      call check(nf90_def_dim(ncid, LON_NAME,NLON,LonDimID))
      call check(nf90_def_dim(ncid, TIME_NAME,NTIME,TimeDimID))
      call check(nf90_def_dim(ncid, 'dateDim',19,dim2id))
      print *, 'Defined Dimensions'

      !Define the coordinate variables 
      call check(nf90_def_var(ncid, LAT_NAME,NF90_REAL, LatDimID,Latvarid))
      call check(nf90_def_var(ncid, LON_NAME,NF90_REAL,LonDimID, Lonvarid))
      !call check(nf90_def_var(ncid, TIME_NAME,NF90_INT,TimeDimID,Timevarid))
      call check(nf90_def_var(ncid, STIME_NAME,NF90_CHAR,(/ dim2id,TimeDimID /),Time2varid))
      call check(nf90_def_var(ncid,TVAR_NAME,NF90_REAL,(/ LonDimId, LatDimID, TimeDimID /), TVarId))

      call check(nf90_put_att(ncid, nf90_global, "DataInterval_hours", dt))
      call check(nf90_put_att(ncid, nf90_global, "Data Type", sdate))
      ! Assign units attributes to coordinate variables.
      call check( nf90_put_att(ncid, latvarid, UNITS, LAT_UNITS) )
      call check( nf90_put_att(ncid, lonvarid, UNITS, LON_UNITS) )
      !call check( nf90_put_att(ncid, Timevarid, UNITS, TIME_UNITS) )
      print *, 'Wrote TIME UNITS'
      call check( nf90_put_att(ncid, Time2varid, UNITS, TIME_UNITS2) )
      call check( nf90_put_att(ncid, TVarid, UNITS, TVAR_UNITS) )
      !	print *, 'Defined Attributes'
      ! End define mode.
      call check( nf90_enddef(ncid) )
      ! Write the coordinate variable data. This will put the latitudes
      ! and longitudes of our data grid into the netCDF file.
      call check( nf90_put_var(ncid, LATvarid, LAT) )
      !	print *, 'Wrote LAT'
      call check( nf90_put_var(ncid, LONvarid, LON) )
      !	print *, 'Wrote LON'
      !call check( nf90_put_var(ncid, TIMEvarid, TIME) )
      print *, 'Wrote TIME'
      call check( nf90_put_var(ncid, TIME2varid, STIME) )
      print *, 'Wrote STIME'

      ! Write the variables
      call check(nf90_put_var(ncid,Tvarid, TVAR))
      print *, 'Wrote TVAR'
      call check(nf90_close(ncid))
   end subroutine write_nc3

END PROGRAM ConvertData2Netcdf
