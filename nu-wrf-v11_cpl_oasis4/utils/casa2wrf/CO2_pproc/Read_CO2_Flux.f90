program read_CO2

   ! This program read the CO2 data in binary format
   ! Convert it to netcdf files with timestamp
   !
   !Author: Jossy P. Jacob (ASTG/CISTO/GSFC/NASA/SSAI)   
   !
   ! Compile: make -f Makefile_CO2_Flux 
   !Run:./Read_CO2_Flux.x filename indir 
   !Example:./Read_CO2_Flux.x flux_e624_co2_3hrly_288x181_2010 flux/
   ! filename: name of the flux file (yearly)
   ! indir : name of the directory where the output will be written as well as flux/binary/filename exists.   
   ! 
   use netcdf
   implicit none
   INCLUDE 'netcdf.inc'

   ! Input Arguments
   integer           :: IARGC
   character*256                :: file1,indir,infile1,outfile,ofile,casadir
   character*30                 :: sdate
   character*10                  :: sdate1

   integer                      :: vardim, status
   integer                      :: ntime, NLON,NLAT,k,i,ntime1,ntime2,it
   integer                      :: TVarId, SVarId
   character*100                :: VarUnit1, VarUnit2
   real                         :: maxdiff,dt1
   integer                      :: year, dt,ntime3,enddate,startdate

   real, allocatable ::  LAT(:), LON(:)
   real, allocatable :: VAR1(:,:,:) 
   real, allocatable :: VAR2(:,:,:)

   integer, allocatable :: TIME(:),TIME2(:)

   character(len=19), allocatable :: STIME(:),TIMES2(:)

   ! ********************************************************************************************
   call GETARG (1, file1)
   print *,trim(file1)
   call GETARG (2, indir)
   print *,trim(indir)
   ! call GETARG (2, syear)
   ! read(syear,*)year
   ! print *, year
   ! call GETARG (3, sdt)
   ! read(sdt,*)dt
   read(file1(29:32),*) year
   read(file1(15:15),*) dt
   read(file1(15:15),*) dt1
   read(file1(21:23),*) NLON
   read(file1(25:27),*) NLAT
   print *,year,dt,NLON, NLAT
   sdate = file1(15:15)
   !casadir = 'CASA/flux/'
   casadir = indir 
   infile1 = trim(casadir)//'/binary/'//trim(file1)
   ofile = 'CO2flux_'//file1(29:32)
   outfile =trim(casadir)//trim(ofile)//'.nc'
   print *, 'file1: ',trim(infile1)
   print *, 'OUT: ',trim(outfile)
   !
   ! For 2D Variables:
   !
   call read_data(infile1,VAR1,LON,LAT,NLON,NLAT,dt,year,time,stime,ntime)
   print *, 'read 2D CO2'
   call write_nc3(outfile,VAR1,LON,LAT,time,stime,NLON,NLAT,ntime,sdate,dt1)
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
   subroutine read_data(infile,VAR,LON,LAT,NLON,NLAT,dt,year,time,stime,ntime)
      character*256,   intent(in)      ::infile 
      integer, intent(in)             :: year, dt,NLON,NLAT
      integer, intent(out)            :: ntime
      real, intent(out), allocatable  :: VAR(:,:,:)
      real, intent(out), allocatable :: LON(:)
      real, intent(out), allocatable :: LAT(:) 
      integer, intent(out), allocatable  :: time(:)
      character(len=19), intent(out), allocatable  :: stime(:)

      real, allocatable :: temp(:,:), temp_prev(:,:)
      character*12                :: varname
      integer                     :: ncid,  i,j,m, n, itime,nt
      integer                     :: ieof,recpos
      real 			:: dx, dy, LON0, LAT0

      !NLON = 288
      !NLAT = 181
      dx = 360.0/NLON
      dy = 180./(NLAT-1)
      LON0 = -180.+ dx/2
      LAT0= -90.
      call get_timestamp(dt,time,year,ntime,stime)
      allocate (LON(1:NLON))
      allocate (LAT(1:NLAT))
      allocate (VAR(1:NLON,1:NLAT,1:ntime))
      allocate (temp(1:NLON,1:NLAT))
      allocate (temp_prev(1:NLON,1:NLAT))		!by Z. Tao to record flux at previous time stemp
      do i=1,NLON
         LON(i)= LON0+ dx * float(i-1)
      enddo
      do j=1,NLAT
         LAT(j)= LAT0+ dy * float(j-1)
      enddo
      !print *, 'LON: ',LON
      !print *, 'LAT: ',LAT
      print *, 'NTIME=',NTIME
      open(unit=11, file=trim(infile), form='unformatted',convert='big_endian',access='direct',recl=NLON*NLAT)
      !ntime = 24/dt*365 ! 1 year data in 1 fie
!      recpos = 0 
      ieof =0

! by Z. Tao: flux at every time step is actually at the middle of each time interval, e.g., flux for hr 0-3 is at hr 1.5, and so on.
! Therefore we take the following procedure to allocate flux to each hour:
! hour = 0 of Jan. 1, flux equals to 1st record
! for other hours at 0, 3, 6, 9, 12, 15, 18, 21, flux = 0.5 (previous time flux + current time flux)
 
      recpos = 1 			!by Z. Tao for very first record
      read(unit=11,iostat=ieof,rec=recpos) temp
      if (ieof == 0) then 
         VAR(:,:,1)= temp
	 temp_prev = temp
      else
         print *, 'File reading error, ntime=1'
         return
      endif

      do nt = 2, ntime			!by Z. Tao, starting from 2nd time step
         recpos = recpos + 1 
         read(unit=11,iostat=ieof,rec=recpos) temp
         if (ieof == 0) then 
            VAR(:,:,nt)= 0.5 * (temp + temp_prev)
	    temp_prev = temp
         else
            print *, 'File reading complete, ntime=',nt
            return
         endif
         !print *,'nhours, min,max are ',nt,minval(VAR(:,:,1),mask=VAR(:,:,1)>-999),maxval(VAR(:,:i,1),mask=VAR(:,:,1)<999) 
      enddo

   end subroutine read_data
   !---------------------
   subroutine get_timestamp(dt,time,year,ntime,stime)
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
      if(mod(year,400) == 0 .or. (mod(year,4) == 0 .and. (mod(year,100) .ne. 0))) then          !by ZTao
         isleapyear = 1
      else 
         isleapyear = 0
      endif
      if (isleapyear == 1) then
         mdd = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
         print *, year, ' leapyear'
         ntime = 24/dt * 366
      elseif (isleapyear == 0) then 
         mdd = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
         ntime = 24/dt * 365
         print *, year, ' not leapyear'
         print *, 'ntime', ntime
      endif
      allocate (time(ntime))
      allocate (stime(ntime))
      ! Get timestamp for 1 year
      itime=0
      do mon = 1,12 
         do day = 1, mdd(mon) 
            do hr = 1, 24,dt
               hour = (hr-1)
               itime=itime+1 
               !time yearmmddhh
               time(itime) = year * 10**6 + mon *10**4 + day * 10**2 + hour
               write(syear,'(I4)') year
               write(smon,'(I2.2)') mon
               write(sday,'(I2.2)') day
               write(shr,'(I2.2)') hour
               stime(itime) = syear//'-'//smon//'-'//sday//'_'//shr//':00:00'
            enddo
         enddo
      enddo
      !              print *,'ntime = ',ntime
       print *,'time -new = ',time
       print *,'stime-new = ',stime
   end subroutine get_timestamp


   !............................................................................................
   subroutine write_nc3(fname_out,TVAR,LON,LAT,TIME,STIME,NLON,NLAT,ntime,sdate,dt)
      !............................................................................................
      use netcdf
      implicit none
      !    INCLUDE 'netcdf.inc'
      integer, intent(in)          :: NLON,NLAT,ntime
      character*256                :: fname_out
      character *30                :: sdate
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
      character (len =*), parameter :: TVAR_NAME = "CO2" 

      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: TVAR_UNITS = "ppmv"
      character (len = *), parameter :: LAT_UNITS = "degrees_north"
      character (len = *), parameter :: LON_UNITS = "degrees_east"
      !character (len = *), parameter :: TIME_UNITS = "YYYYMMDDHH"
      character (len = *), parameter :: TIME_UNITS2 = "YYYY-MM-DD_HH:MI:SS"
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

END PROGRAM read_CO2
