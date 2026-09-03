	subroutine read_rss_oisst_v4(file_name,sst_data,error_data,mask_data,file_exists)
!	subroutine read_rss_mw_sst(file_name,sst_data,error_data,mask_data,file_exists)

!	This routine reads version-4 RSS MW-only, daily SST files made from any of the available 
!     microwave radiometers capable of measuring SST: TMI, AMSR-E, AMSR2, WindSat, GMI
!	You must UNZIP FILES before reading them
!
!	INPUT
! 	file_name  with path in form mw.fusion.yyyy.doy.v04.0
!	where yyyy = year
!		 doy = day of year
!
!	OUTPUT  
!	sst_data  (a 1440 x 720 real*4 array of Sea Surface Temperature)
!	error_data(a 1440 x 720 real*4 array of the interpolation error estimate)
!	mask_data (a 1440 x 720 integer*1 array of data masking information)
!			mask_data 
!				bit 0 = 1 if land is present			
!				bit 1 = 1 if sea ice is present
!				bit 2 = 1 if IR satellite data were utilized for this pixel (will not happen in this data set)
!				bit 3 = 1 if MW data are utilized for this pixel (if no MW data are available, only previous day data exist)	
!				bit 4 = 1 if the pixel is near land and data is suspect (bad data)	
!				bit 5 = not used
!				bit 6 = not used
!				bit 7 = not used
!
!     file_exists  = 1 if file is read and the data returned,  = 0 no file
!
!	xcell = grid cell values between 1 and 1440
!     ycell = grid cell values between 1 and  720
!     dx = 360./1440.
!     dy = 180./720.
!	Center of grid cell Longitude  = dx*xcell-dx/2.    (degrees east)
!	Center of grid cell Latitude   = dy*ycell-(90+dy/2.)  (-90 to 90)
!
!	Please read the data description on www.remss.com
!	To contact RSS support:
!	http://www.remss.com/support
!
!	updated 07/2014 d.smith


	integer*4,parameter		:: xdim = 1440
	integer*4,parameter		:: ydim = 720
	character(len=150)		:: file_name
	real*4,dimension(xdim,ydim)	:: sst_data,error_data
	integer*1,dimension(xdim,ydim):: mask_data
	integer(4)				:: file_exists    

	character(len=1),dimension(1440,720):: buffer
	logical lexist
	
	real(4),parameter			:: scale_sst  =  0.15
	real(4),parameter			:: offset_sst = -3.0
	real*4,parameter			:: scale_error  = 0.005
	real*4,parameter			:: offset_error = 0.0

	 
!	check to see if file exists -- if not return a -1 in file_exists
	file_exists=0
	inquire(file=file_name,exist=lexist)
!	if(.not. lexist) return
!        file_exists=1
        if(.not. lexist) then
                file_exists = -1
                return
        endif

!	open the file and read in character data
	write(*,*) 'reading sst file: ', file_name
	
!	open(3,file=file_name,status='OLD',RECL=1036800,access='DIRECT',form='UNFORMATTED')
!	read(3,rec=1) buffer
!	sst_data = real(ichar(buffer))	!convert character data to real
!	read(3,rec=2) buffer
!	error_data = real(ichar(buffer))	! convert character data to real
!	read(3,rec=3) buffer
!	mask_data = ichar(buffer)		!convert character data to integer1
!	close(3)

!	if you have trouble reading the file, try this alternate set of read code
! Original, not standard Fortran
!	open(3,file=file_name,status='OLD',form='binary')
	open(3,file=file_name,status='OLD',form='unformatted',access='stream') ! NUWRF

	read(3) buffer
	sst_data = real(ichar(buffer))	!convert character data to real
	read(3) buffer
	error_data = real(ichar(buffer))	!convert character data to real
	read(3) buffer
	mask_data=ichar(buffer)			!convert character data to integer1
	close(3)

!	scale the data 
	where(sst_data<=250)
!		scale SSTs using scaling and offset parameters
		sst_data = sst_data * scale_sst + offset_sst
	endwhere
		
	where(error_data<=250)
!		scale Errors using scaling and offset parameters
		error_data = error_data * scale_error + offset_error
	endwhere
!
!	section to set land/ice/bad data to flag values in SST and error data arrays using the mask information
	do i=1,xdim
	do j=1,ydim

!		add mask to SST and ERROR data
		if(btest(mask_data(i,j),0)) then
			sst_data(i,j)=255	!land
			error_data(i,j)=255
		endif
		if(btest(mask_data(i,j),1)) then
			sst_data(i,j)=252 !sea ice
			error_data(i,j)=252
		endif
		if(btest(mask_data(i,j),4)) then
			sst_data(i,j)=254	!bad data
			error_data(i,j)=254
		endif
		
	enddo !j
	enddo !i


	return
	end


