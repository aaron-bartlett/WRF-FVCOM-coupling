!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_thermo
 implicit none

!---------------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
! Comments:  
!   This module compute thermodynamics indices.
! 
! History: 
!  03/2018  Toshi Matsui @NASA GSFC: Initial
!
! References: 
!-----------------------------------------------------------------------
 private  !everything is unless declare it as public. 

!
! public subroutines
!
 public  :: compute_cape

 save


 real,parameter :: no_data_merra = 1.e+15


!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine compute_cape(no_data, kmax, q1d, h1d, dh1d, t1d, rh1d, p1d , CAPE, NCAPE, LCL)
 implicit none


 real,intent(in) :: no_data
 integer,intent(in) :: kmax  
 real,dimension(kmax),intent(in) :: q1d  ! water vapor mixing ratio [g/g]
 real,dimension(kmax),intent(in) :: h1d  ! height [m]
 real,dimension(kmax),intent(in) :: dh1d ! thickness [m]
 real,dimension(kmax),intent(in) :: t1d  ! air temperature [K]
 real,dimension(kmax),intent(in) :: rh1d ! relative humidity [-] (0~1)
 real,dimension(kmax),intent(in) :: p1d  ! air pressure [mb]

 real,intent(out) :: CAPE   ! convective available potential energy [J/kg]
 real,intent(out) :: NCAPE  ! Normalized CAPE (column mean buoyant force) [J/kg/m] or [m/s2]
 real,intent(out) :: LCL    ! lifting condensation level [m]

 real,dimension(kmax) :: bf1d  !buoyant force [m/s2]

 integer :: k
 real,dimension(kmax) :: tv1d, due_dT1d, qs1d
 real :: es, e, Td
 real :: Zlcl ! lifting condensate model [km]
 real :: T_surf ! temperature at surface [K]
 real :: T_parcel ! parcel temp [K]
 integer :: k_surf, k_lcl
 real,dimension(kmax) :: t1d_parcel, tv1d_parcel
 real :: q_surf  !surface H2o mixing ratio [g/g]
 real :: mlap1,mlap2,mlap  ! saturated-adiabatic lapse rate [K/km]
 real :: totwork  !total buoyancy [J/kg]
 real :: tothgt   !total depth [m]

!
! initialize outpout data
!
 CAPE = no_data ; NCAPE = no_data
 LCL  = no_data 

 t1d_parcel(:) = 10.e0  !initialization

!
! Derive extra 1D parameter at pressure level
!

 do k = 1, kmax
        es = 0.611e0 * exp( 5423.e0 * (1.e0/273.e0 - 1.e0 / t1d(k)) ) !saturated vapor pressure [kPa]
        e  = rh1d(k) * es ! vapour pressure [kPa]
        Td = 1.e0  /  (  1.e0/273.e0  - 0.0001844e0 * log( e / 0.611e0) ) !due-point temperature [K]
        due_dT1d(k) = (t1d(k) - Td) ! due-point depression [K]
        tv1d(k) =  t1d(k) * ( 1.e0 + 0.61e0*q1d(k) )  !virtual temp [K]

        qs1d(k) = 0.622e0 * es / (p1d(k)*1.0e-1 - es)  !saturated h2o mixing rato [g/g]
 enddo

!
! let rise surface parcel to estimate ( parcel_tv1d(k) profile through dry/moist adiabatice process ) 
!

 !
 ! estimate height of lifting condensate level
 !
 k_surf = 1 !initialize
 k = k_surf
 Zlcl = MAX(0.01e0, 0.125e0 * due_dT1d(k) ) + h1d(k)*1e-3  !LCL geopotential height [km] 
 T_surf = t1d(k)   ! surface air temperature [K]
 q_surf =  q1d(k_surf)   ! surface H2O mixing ratio [g/g]

 !
 ! Dry Adiabatice Ascend (SURF-> LCL) estimate k level of lifting condensate level, and also lift parcel until Zlcl
 ! bug found here
 k_lcl = -999 ; LCL = no_data 
   do k = k_surf+1, kmax-1

    T_parcel = -9.75e0 *  ( h1d(k) - h1d(k_surf) ) * 1.e-3 + T_surf !adiabatic parcel temperature [K]
    t1d_parcel(k) = T_parcel
    tv1d_parcel(k) =  T_parcel * (1.e0 + 0.61e0*q_surf)  !adiabatice parcel virtual temp [K]

    if( Zlcl*1.e+3 >= 0.5e0*(h1d(k-1)+h1d(k)) .and. Zlcl*1.e+3 <  0.5e0*(h1d(k)+h1d(k+1))  ) then
        k_lcl = k  !LCL k level
        LCL = h1d(k) ! LCL [m]
        exit
    endif
   enddo

 if( k_lcl == no_data ) then
    CAPE=0. ;  NCAPE=0. ;  LCL=0.  
    return 
 endif


 !
 ! Moist Adiabatic Ascend ( LCL -> LOF ) 
 !

 do k = k_lcl+1, kmax

    call moist_lapse( qs1d(k), t1d(k), mlap1 )

    call moist_lapse( qs1d(k-1), t1d(k-1), mlap2 )


    mlap = 0.5e0 * (mlap1 + mlap2)  ! mean saturated-adiabatic lapse rate [K/km]


    t1d_parcel(k) = MAX( 10.e0 , 1.e-3 * ( h1d(k)-h1d(k-1) ) * mlap + t1d_parcel(k-1) )

!print*,k, qs1d(k), t1d(k), mlap1, mlap2, t1d_parcel(k)

!    print*,'t1d_parcel(k)',k,t1d_parcel(k), 1.e-3 * ( h1d(k)-h1d(k-1) ) * mlap + t1d_parcel(k-1)

 enddo


!
! compute buyancy profile
!
 do k = k_surf, kmax

if( t1d_parcel(k) == 0. ) print*,'t1d_parcel(k) is zero'

    es = 0.611e0 * exp( 5423.e0 * (1.e0/273.e0 - 1.e0 / t1d_parcel(k)) ) !saturated vapor pressure [kPa]
    qs1d(k) = 0.622e0 * es / (p1d(k)*1.0e-1 - es)                        !saturated h2o mixing rato for parcel [g/g]
    tv1d_parcel(k) = t1d_parcel(k) * (1.e0 + 0.61*qs1d(k))
    bf1d(k) = max(0.e0, 9.8 * (tv1d_parcel(k) - tv1d(k)) / tv1d(k))     !buoyant force [m/s2] OUTPUT (postive only)
 enddo

!
! compute CAPE (integrate only positive value)
!
 totwork=0.e0  ; tothgt = 0.e0

 do k = k_surf+1, kmax
    if( bf1d(k) > 0. ) then
        totwork = totwork + bf1d(k) * dh1d(k)  ! [m2/s2] --> [J/kg]
        tothgt  = tothgt  + dh1d(k)            ! [m]
    endif
!    totwork = totwork + 0.5e0*(bf1d(k)+bf1d(k-1)) * (h1d(k) - h1d(k-1))   ! [m2/s2] --> [J/kg]
 enddo

!
! output
!
 CAPE = totwork  ! total work -> CAPE [J/kg] OUTPUT

 if( tothgt > 0. ) then
   NCAPE = CAPE / tothgt ! Normalized CAPE (column mean buoyant force) [J/kg/m] or [m/s2]
 else
   NCAPE = 0.
 endif

 return
 end subroutine compute_cape

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine moist_lapse( r_h2o, T, moist_lapse_rate )
 implicit none
 real,intent(in) :: r_h2o !water vapor mixing ratio [g/g] 
 real,intent(in) :: T     !temperature
 real,intent(out) :: moist_lapse_rate  !saturated-adiabatic lapse rate [K/km]

 real :: Cp  ! specific heat  [J / kg / K]
 real,parameter :: g     = -9.807e0  ! [m/s2]  --> [ J / m / kg ]
 real,parameter :: Rd    =  287.053  ! gas constnt for dry air [J/K/kg]
 real,parameter :: Lv    =  2.501e+6 ! latent heat of vaporiz at 0degC [ J / kg ]
 real,parameter :: eta   =  0.622e0
 real :: a,b
 real :: g_Cp

 Cp = 1004.67e0 * (1.e0 + 0.84 * r_h2o )

 g_Cp = 1.e3 * (g / Cp)

 a  = ( 1.e0 + ( r_h2o * (Lv/Rd) ) / T )

 b  = ( 1.e0 + ( (Lv/Cp) * (Lv/Rd) * r_h2o * eta ) / (T*T) )

 moist_lapse_rate = g_CP * a / b  ![K/km]

 return
 end subroutine moist_lapse

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_thermo

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

