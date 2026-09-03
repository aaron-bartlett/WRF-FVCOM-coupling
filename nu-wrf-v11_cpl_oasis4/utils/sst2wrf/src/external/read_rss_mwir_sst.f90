	subroutine read_rss_mwir_sst(file_name,sst_data,error_data,mask_data,file_exists) 
	
!	This routine reads version-4.0 RSS MW-IR Optimal Interpolation SST fusion daily files
!	You must UNZIP FILES before reading them
!
!	INPUT
! 	file_name  with full path in form: mw_ir.fusion.yyyy.doy.v04.0
!       yyyy= year
!	  doy	= day of year
!
!	OUTPUT  
!	sst_data   (a 4096 x 2048 real*4 array of SST data)
!	error_data (a 4096 x 2048 real*4 array of the interpolation error estimate)
!	mask_data  (a 4096 x 2048 integer*1 array of data masking information)
!		   mask_data bits:
!			bit 0 = 1 if land is present			
!			bit 1 = 1 if sea ice is present
!			bit 2 = 1 if IR satellite data were utilized for this pixel (if no IR data are available, only MW data exist)
!			bit 3 = 1 if MW data are utilized for this pixel (if no MW data are available, only IR data exist)	
!			bit 4 = 1 if the pixel is near land and data is suspect (bad data)	
!			bit 5 = 1 if the pixel is in an unclassified sea ice region near land.  This could be sea ice or not.
!			bit 6 = not used
!			bit 7 = not used
!
!	file_exists  = 1 if file was read and data returned,  = 0 if no file found
!
!	xcell = grid cell values between 1 and 4096
!	ycell = grid cell values between 1 and 2048
!	dx = 360./4096.  ! ~9km lat/lon grid 
!	dy = 180./2048.
!	Center of grid cell Longitude = is dx*xcell-dx/2.       degrees east
!	Center of grid cell Latitude  = is dy*ycell-(90+dy/2.)  -90 to 90
!
! In this sample program the mask values are used to set data in the SST and Error arrays to:
!      252 = sea ice present
!      254 = missing data, near land
!      255 = land in this pixel
!      2252 = possible sea ice (EMK NUWRF)
!
!	Please read the data description on www.remss.com/measurements/sea-surface-temperature
!	To contact RSS support:
!	http://www.remss.com/support
!
!	4/2010 c.gentemann
!    updated 8/2010 d.smith
!    updated 1/2014 d.smith
!    updated to version-4.0 11/2014 d.smith


	integer*4,parameter		:: xdim = 4096
	integer*4,parameter		:: ydim = 2048
	character(len=150)		:: file_name
	real*4,dimension(xdim,ydim)	:: sst_data, error_data
	integer*1,dimension(xdim,ydim):: mask_data
	integer*4 				:: file_exists   

	character*1				:: buffer(xdim,ydim)
	logical lexist

	real*4,parameter			:: scale_sst    = 0.15
	real*4,parameter			:: offset_sst   = -3.0
	real*4,parameter			:: scale_error  = 0.005
	real*4,parameter			:: offset_error = 0.0

	 
!	check to see if file exists -- if not return a -1 in file_exists
	file_exists = 0
	inquire(file=trim(file_name),exist=lexist)
!	if(.not. lexist) return
!	file_exists = 1
        
        if(.not. lexist) then
                file_exists = -1
                return
        endif

!	open the file and read in character data
	write(*,*) 'reading sst file: ', trim(file_name)

! Original, which is not standard Fortran
!	open(3,file=trim(file_name),status='OLD',form='binary')
	open(3,file=trim(file_name),status='OLD',form='unformatted',access='stream') ! NUWRF

	read(3) buffer
	sst_data = real(ichar(buffer))  !convert character data to real
	read(3) buffer
	error_data = real(ichar(buffer))  !convert character data to real
	read(3) buffer
	mask_data=ichar(buffer)		!convert character data to integer1
	close(3)


!     if you have trouble reading the file, try this set of code instead
!	open(3,file=trim(file_name),status='OLD',RECL=8388608,access='DIRECT',form='UNFORMATTED')
!	read(3,rec=1) buffer
!	sst_data = real(ichar(buffer))
!	read(3,rec=2) buffer
!     error_data=real(ichar(buffer))
!	read(3,rec=3) buffer
!	mask_data=ichar(buffer)
!	close(3)
!
!	scale the data
	where(sst_data<=250)
!		scale SSTs using scaling and offset parameters
		sst_data = sst_data * scale_sst + offset_sst
	endwhere
      where(error_data<=250)
!		scale Errors using scaling and offset parameters
		error_data = error_data * scale_error + offset_error
	endwhere

!     section set land/ice/bad data to flag values using the mask information
	do i=1,xdim
	do j=1,ydim
!				
!		add mask to data
		if(btest(mask_data(i,j),0)) then
			sst_data(i,j)=255  !land
			error_data(i,j)=255
		endif
		if(btest(mask_data(i,j),1)) then
			sst_data(i,j)=252   !sea ice
			error_data(i,j)=252
		endif
		if(btest(mask_data(i,j),4)) then
			sst_data(i,j)=254	  !bad data 
			error_data(i,j)=254
		endif
		
		!users are advised to make their own judgement whether to flag these data (comment out to use)
                ! EMK NUWRF...Preserve SST value.  Create "possible sea ice"
                ! flag.
		if(btest(mask_data(i,j),5)) then
!			sst_data(i,j)=252	  !no good ice mask to know, but could be sea ice 
!			error_data(i,j)=252
			error_data(i,j)=2252

		endif		
		
	enddo  !j
	enddo	 !i

	return
	end


