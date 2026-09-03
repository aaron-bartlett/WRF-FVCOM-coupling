; builds a forcing file for the SCM from met_em* input.
; 1.  Stuffs in data
; 2.  Modifies global attributes as appropriate
; 3.  Writes in initial profile for soil and atmosphere 
; 4.  Writes a suggested scm namelist
;
; NOTES: 
; 1.  Required environment variables:
;     - the ensemble-center date/time in "SCM_CENTER" following the WRF
;       file naming convention, e.g. setenv SCM_CENTER 2000-06-10_00:00:00
;     - output file name in FORCE_FILE_NAME
;     - path of directory containing met_em* or other input files in METPATH
;     - length of the simulation (seconds) in SIMLENGTH
;     - random seeds in RANDSEEDa and RANDSEEDb, for reproducability 
; 2.  If the ensemble member is specified as 0,
; no perturbation will occur - you will get a single file for the center date
; 3.  Because it's a pain to put environment variables on the NCL command line
; (even possible?), this uses getenv instead.
; 4.  Use a random seed if you want repeatable results
; 5.  Make sure that the number of levels in the forcing file is equal to the
; input data.  Note that metgrid files put the surface diagnostics in the
; first layer.  We don't want to include those.
; 6.  For now we assume that a 3x3 grid is already set up, and we just take the center grid point.  The input stencil should be:
;      v   v
;    u M u M u
;      v X v   <--- forcing is valid at the X
;    u M u M u
;      v   v
; 7.  Built on NARR data, which has the sfc values in the first level. These
; are stripped out for the forcing but not for input_sounding.  Is this
; general to all or at least most met_em* files?
; 8.  Right now this does not do any extrapolation, which will be a problem
; in the case that the lowest pressure-level data (e.g. 1000 mb) higher than
; the lowest WRF level. **
; 9.  Rotates the winds from the input grid to zonal and meridional.
; 10. No consideration of the map scale factors at this point.  Make sure
; that your input met_em* are centered on the SCM location and that you are
; not introducing much grid stretching with the projection.  **
;
; ** to be handled in future versions.
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./build_forcing_functions.ncl"
load "./rotwind.ncl"

begin

; any fixed parameters?
maxEnsembleSize = 1000

;--------------------------------------------------------
; END user modifications
;--------------------------------------------------------

; check for required environment variables

fName = getenv("FORCE_FILE_NAME")
if ( ismissing(fName) ) then
  print("Set environment variable FORCE_FILE_NAME")
  exit
end if
print("Filling forcing file "+fName)

metPath = getenv("METPATH")
if ( ismissing(metPath) ) then
  print("Set environment variable METPATH")
  exit
end if
print("With input from directory "+metPath)

ensMem = stringtointeger(getenv("ENSEMBLE_MEMBER"))
perturb = True
if ( ismissing(ensMem) .or. ensMem .lt. 1 ) then
  print("No perturbation")
  perturb = False
end if

dateStr = getenv("CENTER_DATE")
if ( ismissing(dateStr) ) then
  print("Set environment variable CENTER_DATE")
  exit
  ; NOTE: could use this one to default to a random center
end if

simLen = stringtodouble(getenv("SIMLENGTH"))
if ( ismissing(simLen) ) then
  print("Set environment variable SIMLENGTH")
  exit
end if

xll = stringtointeger(getenv("XLL"))
if ( ismissing(xll) ) then
  print("Set environment variable XLL")
  exit
end if

yll = stringtointeger(getenv("YLL"))
if ( ismissing(yll) ) then
  print("Set environment variable YLL")
  exit
end if

rSeed1 = stringtointeger(getenv("RANDSEEDa"))
rSeed2 = stringtointeger(getenv("RANDSEEDb"))

if ( .not. ismissing(rSeed1) .and. .not. ismissing(rSeed2) ) then
  random_setallseed(rSeed1,rSeed2)
else
  print("Using default NCL random seed")
end if

; open the forcing file
forceFile = addfile(fName,"rw")

; get date string length
timeDims = getfilevardimsizes(forceFile,"Times")
dateStrLen = timeDims(1)

; get times for all files in directory
fList = systemfunc("ls "+metPath+"/met_em.d01*")
nFiles = dimsizes(fList)
fileDates = new((/nFiles,dateStrLen/),"character")
do iFile = 0, nFiles - 1
  inF = addfile(fList(iFile),"r")
  fTime = inF->Times(0,:)
  fileDates(iFile,:) = fTime
  delete(fTime)
  delete(inF)
end do 

; convert to integers
fileYrList = stringtointeger(chartostring(fileDates(:,0:3)))
fileMoList = stringtointeger(chartostring(fileDates(:,5:6)))
fileDyList = stringtointeger(chartostring(fileDates(:,8:9)))
fileHrList = stringtointeger(chartostring(fileDates(:,11:12)))
fileMmList = stringtointeger(chartostring(fileDates(:,14:15)))
fileSsList = stringtointeger(chartostring(fileDates(:,17:18)))
fileJul = greg2jul(fileYrList,fileMoList,fileDyList,fileHrList)

; find start file for center
dateChar = stringtochar(dateStr)
centerYrStart = stringtointeger(chartostring(dateChar(0:3)))
centerMoStart = stringtointeger(chartostring(dateChar(5:6)))
centerDyStart = stringtointeger(chartostring(dateChar(8:9)))
centerHrStart = stringtointeger(chartostring(dateChar(11:12)))
centerMmStart = stringtointeger(chartostring(dateChar(14:15)))
centerSsStart = stringtointeger(chartostring(dateChar(17:18)))

startJul = greg2jul(centerYrStart,centerMoStart,centerDyStart,centerHrStart)
startJul = startJul + centerMmStart/1440 + centerSsStart/86400

; find all files within range
; NOTE: need last time to compute tendency only - is not written out
endJul = startJul + simLen/86400
centerFileInds = ind(fileJul .ge. startJul .and. fileJul .le. endJul)
centerFiles = fList(centerFileInds)
centerJul   = fileJul(centerFileInds)
nCenterFiles  = dimsizes(centerFileInds)

; make sure beginning and ending are in the list
if ( startJul .ne. centerJul(0) ) then
  print("Could not find a file to start simulation at "+dateStr)
  exit
end if
if ( endJul .ne. centerJul(nCenterFiles-1) ) then
  print("Could not find a file to end simulation of "+simLen+" s")
  exit
end if

nForceTimes = nCenterFiles-1
vdims = getfilevardimsizes(forceFile,"Z_FORCE")
nForceLevels = vdims(1)
z_force = new((/nForceTimes+1,nForceLevels/),"float")
z_force  = forcing_z(centerFiles,z_force,xll,yll)
zCenterInit = z_force(0,:)

inF = addfile(centerFiles(0),"r")
map_proj = inF@MAP_PROJ
std_lon = inF@STAND_LON
true_lat1 = inF@TRUELAT1
true_lat2 = inF@TRUELAT2

; rotate winds?
rotate = False
if ( map_proj .lt. 4 ) then
  rotate = True
end if


; need lat/long for rotating winds and to write namelist
latTmp = inF->XLAT_M(0,:,:)
lat = 0.25*(latTmp(yll,xll) + latTmp(yll+1,xll) + latTmp(yll,xll+1) + latTmp(yll+1,xll+1))
delete(latTmp)


lonTmp = inF->XLONG_M(0,:,:)
lon = 0.25*(lonTmp(yll,xll) + lonTmp(yll+1,xll) + lonTmp(yll,xll+1) + lonTmp(yll+1,xll+1))
delete(lonTmp)

;-------------------------------------------------------------
; Now deal with initial conditions.
; These are directly from the NARR so the elevation is different.
;-------------------------------------------------------------

; first get surface
f = addfile(centerFiles(0),"r")
sfcvars = get_surface_vars(f,xll,yll)

if ( rotate ) then 
   utmp = gridwind_to_truewind("U",sfcvars@u10,sfcvars@v10,lat,lon,map_proj,std_lon,true_lat1,true_lat2)
   vtmp = gridwind_to_truewind("V",sfcvars@u10,sfcvars@v10,lat,lon,map_proj,std_lon,true_lat1,true_lat2)
   sfcvars@u10 = (/utmp/)
   sfcvars@v10 = (/vtmp/)
   delete(utmp)
   delete(vtmp)
end if
 
outdat = new((/1,6/),"float")
outdat(0,:) = (/sfcvars@terrain,sfcvars@u10,sfcvars@v10,sfcvars@t2,sfcvars@qv2, \\
           sfcvars@psurf/)


mopt = True
mopt@fout = "surface_init.txt"
write_matrix(outdat,"6f12.4",mopt)
delete(outdat)

; now profile
profvars = get_profile_vars(f,zCenterInit,xll,yll)

if ( rotate ) then 
   utmp = gridwind_to_truewind("U",profvars@u,profvars@v,lat,lon,map_proj,std_lon,true_lat1,true_lat2)
   vtmp = gridwind_to_truewind("V",profvars@u,profvars@v,lat,lon,map_proj,std_lon,true_lat1,true_lat2)
   profvars@u = (/utmp/)
   profvars@v = (/vtmp/)
   delete(utmp)
   delete(vtmp)
end if
 

nz = dimsizes(profvars@z)
outdat = new((/nz,5/),"float")
outdat(:,0) = profvars@z
outdat(:,1) = profvars@u
outdat(:,2) = profvars@v
outdat(:,3) = profvars@t
outdat(:,4) = profvars@qv


mopt@fout = "profile_init.txt"
write_matrix(outdat,"5f12.6",mopt)
delete(outdat)


end

