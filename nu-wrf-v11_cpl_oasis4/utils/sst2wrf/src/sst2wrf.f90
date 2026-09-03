!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  SST2WRF
!
! AUTHOR:
! Eric Kemp, NGIS
! Incorporating read routines from Remote Sensing Systems.
!
! DESCRIPTION:
! Reads binary file containing sea surface temperature data from Remote
! Sensing Systems, and outputs in WPS Intermediate Format for use by
! WRF.
!
! REVISION HISTORY:
! 05 Apr 2010 - Initial version
! 15 Apr 2010 - Changed missing/bad SST values to 0 K (will be rejected by 
!               REAL).
!------------------------------------------------------------------------------

program sst2wrf

  ! Use modules
  use nml_mod
  use FieldSST_RSS_mod, only:  FieldSST_RSS, &
       createFieldSST_RSS, destroyFieldSST_RSS, convertToKelvin, setMissingSST
  use OutputHours_mod, only:  OutputHours, &
       createOutputHours,destroyOutputHours
  use FieldWPS_mod, only: FieldWPS, &
       createFieldWPS, destroyFieldWPS, &
       IPROJ_LATLON, SWCORNER
  use FileWPS_mod, only: FileWPS, &
       createFileWPS, destroyFileWPS, writeFileWPS

  ! Force explicit variable declarations
  implicit none

  ! Variable declarations
  type(FieldSST_RSS) :: sstData
  type(OutputHours) :: allHours
  type(FileWPS) :: wpsFileData
  type(FieldWPS) :: wpsFieldData
  character(len=46) :: description
  real :: earthRadius
  integer :: i

  type(nml) :: nml_s2w
  character(len=20) :: instrument
  integer :: year
  integer :: dayOfYear
  character(len=5) :: version
  character(len=256) :: inputDirectory
  character(len=256) :: outputDirectory
  character(len=80) :: prefixWPS
  integer :: numFakeHours
  integer,allocatable :: fakeHours(:)

  nml_s2w = nml_new('namelist.sst2wrf')

  call nml_s2w%read("input","instrument",instrument)
  call nml_s2w%read("input","year",year)
  call nml_s2w%read("input","dayOfYear",dayOfYear)
  call nml_s2w%read("input","version",version)
  call nml_s2w%read("input","inputDirectory",inputDirectory)
  call nml_s2w%read("output","outputDirectory",outputDirectory)
  call nml_s2w%read("output","prefixWPS",prefixWPS)
  call nml_s2w%read("fakeoutput","numFakeHours",numFakeHours)
  allocate(fakeHours(numFakeHours))
  call nml_s2w%read("fakeoutput","fakeHours",fakeHours)
  earthRadius = 111.111111*180./(4.*atan(1.))
  

  ! Fetch SST field
  sstData = createFieldSST_RSS(trim(inputDirectory), &
       trim(instrument), year, &
       dayOfYear, trim(version))
  call convertToKelvin(sstData)
  call setMissingSST(sstData)

  ! Create list of hours to output
  allHours = createOutputHours(numFakeHours, fakeHours, sstData)
  
  ! Loop through all hours
  do i = 1,allHours%numHours
     
     ! Open WPS file
     wpsFileData = createFileWPS(fileUnit=10, &
          outputDirectory=trim(outputDirectory), &
          prefix=trim(prefixWPS), &
          year=sstData%year, &
          month=sstData%month, &
          day=sstData%day, &
          hour=allHours%hours(i))
     
     ! Build WPS SST Field 
     description="Sea Surface Temperature (" // &
          trim(instrument) // &
          " " // trim(version) // ")"
     wpsFieldData = createFieldWPS( &
          year=sstData%year, &
          month=sstData%month, &
          day=sstData%day, &
          hour=allHours%hours(i), &
          xfcst=0., &
          map_source="Remote Sensing Systems", &
          field="SST", &
          units="K", &
          desc=trim(description), &
          xlvl=200100., &
          nx=sstData%nx, &
          ny=sstData%ny, &
          iproj=IPROJ_LATLON, &
          startloc=SWCORNER, &
          startlat=sstData%swLat, &
          startlon=sstData%swLon, &
          earth_radius=earthRadius, &
          is_wind_grid_rel=.false. , &
          slab=sstData%sstData, &
          deltalat=sstData%dLat, &
          deltalon=sstData%dLon)
     
     ! Write the SST field
     call writeFileWPS(wpsFileData,wpsFieldData)
     call destroyFieldWPS(wpsFieldData)
     
     ! Build WPS SST Field 
     description="SST Mask (0=SST,1=land,2=seaice,3=missing)"
     wpsFieldData = createFieldWPS( &
          year=sstData%year, &
          month=sstData%month, &
          day=sstData%day, &
          hour=allHours%hours(i), &
          xfcst=0., &
          map_source="Remote Sensing Systems", &
          field="SSTMask", &
          units=" ", &
          desc=trim(description), &
          xlvl=200100., &
          nx=sstData%nx, &
          ny=sstData%ny, &
          iproj=IPROJ_LATLON, &
          startloc=SWCORNER, &
          startlat=sstData%swLat, &
          startlon=sstData%swLon, &
          earth_radius=earthRadius, &
          is_wind_grid_rel=.false. , &
          slab=sstData%sstMask, &
          deltalat=sstData%dLat, &
          deltalon=sstData%dLon)
  
     ! Write the SST mask field
     call writeFileWPS(wpsFileData,wpsFieldData)
     call destroyFieldWPS(wpsFieldData)
  end do

  ! Clean up
  deallocate(fakeHours)
  call destroyFileWPS(wpsFileData)
  call destroyOutputHours(allHours)
  call destroyFieldSST_RSS(sstData)

end program sst2wrf
