;------------------------------------------------------------------------------
; NASA/GSFC, Computational and Information Science and Technology Office, 
; Code 606
;------------------------------------------------------------------------------
;
; PROGRAM: ke_spectra_arw.ncl
;
; AUTHOR:
; Eric Kemp, NASA CISTO/SSAI
;
; DESCRIPTION:
; Calculates and plots mean horizontal kinetic energy spectra from WRF-ARW
; 3D output. Loosely based on ke_spectra_dd.ncl, which performs analogous 
; work on global lat/lon data using spherical harmonics, and on approach
; described by Skamarock (2004).
;
; REFERENCES:
; Errico, R. M., 1985: Spectra computed from a limited area grid. Mon. Wea.
;   Rev., 113, 1554-1562. [See section 2a for detrending description.]
; Koshyk, J. N., and K. Hamilton, 2001: The horizontal kinetic energy spectrum
;   and spectral budget simulated by a high-resolution troposphere-
;   stratosphere-mesosphere GCM. J. Atmos. Sci., 58, 329-348. [See appendix
;   for plane geometry equations.]
; Skamarock, W. C., 2004: Evaluating mesoscale NWP models using kinetic energy
;   spectra. Mon. Wea. Rev., 132, 3019-3032.
; http://www.ncl.ucar.edu/Document/Functions/Built-in/cfftf.shtml [Summarizes
;   function for forward complex discrete Fourier transform.]
;
; REVISION:
; 22 Jul 2015 - First version.
; 03 Sep 2015 - Bug fix for copying data to v_subset (Eric Kemp, Weile Wang)
;
;------------------------------------------------------------------------------

script_name="ke_spectra_arw.ncl"

;------------------------------------------------------------------------------
; Load libraries.
;------------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;------------------------------------------------------------------------------
;
; MAIN PROGRAM DRIVER
;
;------------------------------------------------------------------------------

begin

  ;----------------------------------------------------------------------------
  ;
  ; Check command line arguments
  ;
  ;----------------------------------------------------------------------------

  if (isvar("fnamWRF")) then
      print(script_name+": fnamWRF="+fnamWRF)
  else
      print(script_name+": fnamWRF is a required variable on the command line")
      print(script_name+": Please specify the filename:  fnamWRF")
      exit
  end if

  if (isvar("diri")) then            ; input dir on command line?
      diri = diri + "/"              ; append /
  else
      if (.not.ismissing(getenv("DIRI"))) then   ; environment variable set?
          diri = getenv("DIRI") +"/" ; use environment variable .. if present
      else
          diri = ""                ; default
      end if
  end if

  if (isvar("diro")) then            ; output dir on command line?
      diro = diro + "/"              ; append /
  else
      if (.not.ismissing(getenv("DIRO"))) then   ; environment variable set?
          diro = getenv("DIRO") +"/" ; use environment variable .. if present
      else
          diro = "./"                ; default
      end if
  end if

  if (isvar("pfmt")) then            ; plot format specified on command line?
      pltType = pfmt                 ; command line
  else
      if (.not.ismissing(getenv("PFMT"))) then   ; environment variable set?
          pltType = getenv("PFMT")   ; use environment variable .. if present
      else
          pltType = "x11"            ; default, others: eps, png, pdf 
      end if
  end if

  if (.not.isvar("mlvl")) then       ; model level on command line?
      if (.not.ismissing(getenv("KE_MLEVEL"))) then
          ilvl = getenv("KE_MLEVEL")     ; use environment variable .. if present
          mlvl = stringtointeger( ilvl )
     else
          mlvl =  21                ; default  
      end if
  end if
  pltLevel = "Eta Level "+mlvl

  if (isvar("ptit")) then            ; is plot title specified on the command line?
      pltTitle = ptit                ; command line override if present
  else
      pltTitle = "WRF"
  end if

  if (isvar("pcolor")) then
      pltColor = pcolor
  else
      pltColor = "Turquoise"
  end if

  if (isvar("pnam")) then            ; plot name specified on command line?
      pltName = diro+pnam            ; command line
  else
      if (.not.ismissing(getenv("PNAME"))) then  ; environment variable set?
          pltName = getenv("PNAME")  ; use environment variable .. if present
      else
          pltName = diro+script_name
      end if
  end if
  print ("out image="+pltName+"."+pltType)

;------------------------------------------------------------------------------
; Read WRF file
;------------------------------------------------------------------------------

  print ("diri="+diri+" fnamWRF="+fnamWRF)
  fpsWRF = addfile( diri+fnamWRF+".nc" , "r")

  DX = fpsWRF@DX*0.001 ; in km
  DY = fpsWRF@DY*0.001 ; in km
 
  ; Just grab first time
  U = fpsWRF->U(0,:,:,:)
  V = fpsWRF->V(0,:,:,:)
  MAPFAC_M = fpsWRF->MAPFAC_M(0,:,:)
  Time = fpsWRF->Times(0,:)

  dims_a = dimsizes(MAPFAC_M)
  ny = dims_a(0)
  nx = dims_a(1)

  dims_u = dimsizes(U)
  nz = dims_u(0)

  ; Interpolate winds from staggered C-grid points to A-grid
  u_a = new((/nz,ny,nx/),"float")
  v_a = new((/nz,ny,nx/),"float")
  u_a(:,:,:) = 0.5 * ( U(:,  :,      0:(nx-1)) + U(:,  :,  1:nx) )
  v_a(:,:,:) = 0.5 * ( V(:, 0:(ny-1), :      ) + V(:, 1:ny, :  ) )
  delete(U)
  delete(V)
  delete(MAPFAC_M)

  ; Take spatial subset of winds. Skip 15 points closest to lateral boundaries
  ; following Skamarock (2004). 
  istart = 15
  iend = nx - 1 - 15
  jstart = 15
  jend = ny - 1 - 15

  if (mlvl .lt. 1 .or. mlvl .gt. nz) then
    print("ERROR, invalid model level selected!")
    print("mlvl must be positive and less than " + nz)
    exit()
  end if
  k = mlvl

  u_subset = new((/ny-30,nx-30/), "float")
  v_subset = new((/ny-30,nx-30/), "float")
  u_subset(:,:) = u_a(k, jstart:jend, istart:iend )
  v_subset(:,:) = v_a(k, jstart:jend, istart:iend )
  delete(u_a)
  delete(v_a)

  ; We only need half of the wavenumbers, corresponding to wavelengths >= 
  ; 2*dx. These also correspond to the positive frequencies processed by the 
  ; complex Fourier transform (see on-line page for NCL function cfftf 
  ; referenced above).
  trunc_num = (nx - 30 - 1) / 2

  ; Calculate spectra for each nominal east-west row, and average them across
  ; nominal north-south dimension of domain.
  u_row = new((/nx-30/), "float")
  v_row = new((/nx-30/), "float")
  spectrum_row =  new((/3,trunc_num/), "float")
  mean_spectrum = new((/3,trunc_num/), "float")

  mean_spectrum(:,:) = 0.0

  do j = 0,ny-30-1
    u_row(:) = u_subset(j,:)
    v_row(:) = v_subset(j,:)

    ; Detrend data along row. Follows Errico (1985) eqns 1 and 2, but is only
    ; performed in one dimension.
    s_u_row = (u_row(nx-30-1) - u_row(0)) / ((nx-30) - 1)
    s_v_row = (v_row(nx-30-1) - v_row(0)) / ((nx-30) - 1)
    do i = 0, nx-30-1
      u_row(i) = u_row(i) - 0.5*(2*i - (nx - 30) - 1)*s_u_row
      v_row(i) = v_row(i) - 0.5*(2*i - (nx - 30) - 1)*s_v_row
    end do

    ; TODO - Rotate the winds

    ; Get complex Fourier coefficients.
    ; First coef dimension: 0 = real coefficients, 1 = imaginary
    coef_u = cfftf(u_row, 0.0, 0)
    coef_v = cfftf(v_row, 0.0, 0)

    ; Normalize the coefficients. See web page description for NCL
    ; function cfftf referenced above.
    dims_u = dimsizes(u_row)
    dims_v = dimsizes(v_row)
    coef_u = coef_u / dims_u(0)
    coef_v = coef_v / dims_v(0)

    ; First spectrum_row dimension: 0 = u, 1 = v, 2 = total
    spectrum_row(0,0) = 0.0
    spectrum_row(1,0) = 0.0
    do in = 1, trunc_num-1
      spectrum_row(0,in) = coef_u(0,in)^2 + coef_u(1,in)^2
      spectrum_row(1,in) = coef_v(0,in)^2 + coef_v(1,in)^2

      ; Koshyk and Hamilton (2001), eqn A1, but broken down by u and v
      spectrum_row(0,in) = spectrum_row(0,in) / 4.0
      spectrum_row(1,in) = spectrum_row(1,in) / 4.0
    end do
    spectrum_row(2,:) = spectrum_row(0,:) + spectrum_row(1,:)

    mean_spectrum(:,:) = mean_spectrum(:,:) + spectrum_row(:,:)
  end do
  ; Normalize by number of rows
  mean_spectrum(:,:) = mean_spectrum(:,:) / (ny-30)

  ; Set axes
  x = ispan(0,trunc_num-1,1)*1.0 ; Nominal east-west domain wavenumbers
  x2 = (nx-1-30)*DX / (ispan(1,trunc_num-1,1)) ; Nominal east-west wavelengths
  mean_spectrum!1 = "x"
  x!0        = "x"

  ;----------------------------------------------------------------------------
  ; Plot
  ;----------------------------------------------------------------------------

  wks    = gsn_open_wks(pltType,pltName)

  lthick = 2.0

  ; create custom color map
  colors = (/"white","black","OrangeRed4","Red3",     \
             "red1","Red3","Red4","Blue",             \
             "OrangeRed2","OrangeRed3","OrangeRed4","DarkOrange",\
             "red1","red2","red3","red4", "OrangeRed1",          \
             "DarkTurquoise","Turquoise","Cyan","blue",          \
             "gold","greenyellow","GoldenRod","DarkGoldenRod1",  \
             "lightgoldenrod1","lightyellow","lightgoldenrodyellow","khaki1"/)
  gsn_define_colormap(wks, colors)

  res          = True
  res@xyXStyle = "Log"
  res@xyYStyle = "Log"
  res@tmXBMinorPerMajor = 8
  res@tmXTMinorPerMajor = 8
  res@tmYRMinorPerMajor = 8
  res@tmYLMinorPerMajor = 8

  res@tiMainString      = pltTitle+"    "+pltLevel
  res@tiMainString      = pltTitle+"    "+pltLevel+"~C~Valid "+Time+" UTC"

  ; res@gsnCenterString   = pltLevel
  res@tiXAxisString     = "Wavenumber"
  res@tiYAxisString     = "Mean KE per unit mass (m^2/s^2)"
  res@gsnFrame          = False
  res@xyLineThicknessF  = lthick

  yrangeMin = 1.e-10
  yrangeMax = 1.e2

  res@trXMinF = 1.e0
  res@trXMaxF = 1.e4
  res@trXReverse = False ; reverse X-axis
  res@trYMinF = yrangeMin
  res@trYMaxF = yrangeMax

  res@tmXMajorGrid                 = False   ; Turn on vertical lines
  res@tmXMajorGridThicknessF       = 0.5
  res@tmXMajorGridLineDashPattern  = 1
  res@tmXMinorGrid                 = False   ; Turn on vertical lines
  res@tmXMinorGridThicknessF       = 0.5
  res@tmXMinorGridLineDashPattern  = 2

  resT          = True
  resT@xyXStyle = "Log"
  resT@xyYStyle = "Log"
  resT@tmXBMinorPerMajor = 8
  resT@tmXTMinorPerMajor = 8
  resT@tmYRMinorPerMajor = 8
  resT@tmYLMinorPerMajor = 8
  resT@tiXAxisString     = "Wavelength (km)"
  resT@gsnFrame          = False
  resT@xyLineThicknessF  = lthick
  resT@trXMinF = (nx-1-30)*DX/10000.
  resT@trXMaxF = (nx-1-30)*DX
  resT@trXReverse = True ; reverse X-axis
  resT@trYMinF = yrangeMin
  resT@trYMaxF = yrangeMax
  resT@tmXMajorGrid                 = False   ; Turn on vertical lines
  resT@tmXMajorGridThicknessF       = 0.5
  resT@tmXMajorGridLineDashPattern  = 1
  resT@tmXMinorGrid                 = False   ; Turn on vertical lines
  resT@tmXMinorGridThicknessF       = 0.5
  resT@tmXMinorGridLineDashPattern  = 2

  system("echo [Notice]: plotting")

  res@xyDashPattern = 0
  res@xyLineColors  = (/pltColor,"blue"/)
  plot    = gsn_csm_x2y (wks,x(1:),x2,mean_spectrum(2,1:),res,resT) 

;   ; U component
;   polyres                  = True
;   polyres@gsLineThicknessF = lthick
;   polyres@gsLineColor = "blue"
;   gsn_polyline(wks,plot,x(1:),mean_spectrum(0,1:),polyres)

;   ; V component
;   polyres                  = True
;   polyres@gsLineThicknessF = lthick
;   polyres@gsLineColor = "red1"
;   gsn_polyline(wks,plot,x(1:),mean_spectrum(1,1:),polyres)

  ; Total
  polyres                  = True
  polyres@gsLineThicknessF = lthick
  polyres@gsLineColor = "DarkTurquoise"
  gsn_polyline(wks,plot,x(1:),mean_spectrum(2,1:),polyres)

  system("echo [Notice]: plotting")

  ; -3 line
  xref    = (/1., 200. /)
  c = (0.01^(1./3.))
  spc_ref = (c*xref)^(-3)
  polyres                  = True
  polyres@gsLineThicknessF = lthick
  polyres@gsLineDashPattern = 2
  polyres@gsLineColor      = "black"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline

  ; -5/3 line
  xref    = (/1. , 200./)
  c = (100.)^(3./5.)
  spc_ref = (c*xref)^(-5./3.)
  polyres                  = True
  polyres@gsLineThicknessF = lthick
  polyres@gsLineDashPattern = 1
  polyres@gsLineColor      = "black"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline

  ; WRF 600km nudging line
  xWRFx   = 13 ; New nudging wave number for DSCALE
  xref    = (/xWRFx,xWRFx/)
  spc_ref = (/yrangeMin,yrangeMax/)
  polyres                   = True
  polyres@gsLineThicknessF  = 4
  polyres@gsLineDashPattern = 1
  polyres@gsLineColor       = "Green4"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline

  ; WRF 2000km line
  xWRFx    = 4 ; Original WRF nudging for DSCALE
  xref    = (/xWRFx,xWRFx/)
  spc_ref = (/yrangeMin,yrangeMax/)
  polyres                   = True
  polyres@gsLineThicknessF  = 4
  polyres@gsLineDashPattern = 1
  polyres@gsLineColor       = "DarkOrange2"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline

  system("echo [Finished]: plotting")

  frame(wks)                                     ; now advance

end