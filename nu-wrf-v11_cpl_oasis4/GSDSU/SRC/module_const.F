#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
 module module_const
 use module_floating_point
 implicit none

!--------------------------------------------------------------------------------------------------
! Comments:                             
!         This module manages many G-SDSU universal constant parameters.
!
!  History:
!  03/2012  Toshi Matsui@NASA GSFC : Init 
!
! References:
!-----------------------------------------------------------------------------------------------------

 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 public   ! ALL variables and subourtines are accessible.

!
! SDSU Constants 
!
 integer,parameter :: undefined_i2 = -999
 real(sdsu_fps), parameter :: undefined = -999.e0
 real(sdsu_fps), parameter :: q_min_condensate = 1.0e-5  ! siginificant minimum condensate mixing ratio [g/m3]
 real(sdsu_fps), parameter :: q_min_aerosol    = 1.0e-8  ! siginificant minimum aerosol mixing ratio [g/m3]
 real(sdsu_fps), parameter :: const_g = 9.807e0 ! graviational acceleration on earth at sea level [m/s2]
 real(sdsu_fps), parameter :: const_Rd = 287.053 ! dry-air gas constant [J/K/kg] 
 real(sdsu_fps), parameter :: const_Lv_Rv = 5423.e0 ! Clausius-Clapeyron parameter
 real(sdsu_fps), parameter :: const_Rd_Rv = 0.622e0 ! epsilon 
 real(sdsu_fps), parameter :: const_Kel2Cel  = 273.15e0 !conversion constant from Kelvin to Celcius
 real(sdsu_fps), parameter :: const_pi = 3.14159e0 ! pai 
 real(sdsu_fps), parameter :: const_degrad=const_pi/180.e0  ! conversion factor for degrees to radian
 real(sdsu_fps), parameter :: const_smallest_real = tiny(const_degrad) ! smallest floating point
 real(sdsu_fps), parameter :: const_largest_real = huge(const_degrad)  ! largest floating point
 real(sdsu_fps), parameter :: const_ghz_mm = 2.997925e2

 real(sdsu_fpd), parameter :: const_dpi = 3.141592653589793d0 ! pai 
 real(sdsu_fpd), parameter :: const_smallest_dble = tiny(const_dpi)  ! smallest double precision point
 real(sdsu_fpd), parameter :: const_largest_dble  = huge(const_dpi)  ! largest doubple precision point

 real(sdsu_fpd), parameter :: undefined_dble = dble( undefined )
 real(sdsu_fps), parameter :: tb_or_rad = 3.3  !threshold seperate visible and infrared wavelength [micron]

 end module module_const

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gamma_reff(x,ga_out)

!---------------------------------------------------------------------------------------------------
! Comments:                             
!   compute the gamma function a(x) for single precision floating point. 
!       input :  x  --- argument of a(x)
!                       ( x is not equal to 0,-1,-2,... )
!       output:  ga --- a(x)
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Adapted to SDSU
!           
! References: 
!----------------------------------------------------------------------------------------------------
 implicit double precision (a-h,o-z)
 dimension g(26)
 data g/1.0d0,0.5772156649015329d0, &
       -0.6558780715202538d0, -0.420026350340952d-1, &
        0.1665386113822915d0,-.421977345555443d-1, &
        -.96219715278770d-2, .72189432466630d-2, &
        -.11651675918591d-2, -.2152416741149d-3, &
        .1280502823882d-3, -.201348547807d-4, &
        -.12504934821d-5, .11330272320d-5, &
        -.2056338417d-6, .61160950d-8, &
         .50020075d-8, -.11812746d-8, &
        .1043427d-9, .77823d-11, &
        -.36968d-11, .51d-12, &
        -.206d-13, -.54d-14, .14d-14, .1d-15/
 real(4) :: x, ga_out
 pi=3.141592653589793d0
 if (x.eq.int(x)) then
     if (x.gt.0.0d0) then
         ga=1.0d0
         m1=x-1
        do k=2,m1
           ga=ga*k
        enddo
     else
        ga=1.0d+300
     endif
  else

     if (dabs(dble(x)).gt.1.0d0) then
         z=dabs(dble(x))
         m=int(z)
         r=1.0d0
        do k=1,m
           r=r*(z-k)
        enddo
        z=z-m
     else
        z=dble(x)
     endif
     gr=g(26)
     do k=25,1,-1
        gr=gr*z+g(k)
     enddo
     ga=1.0d0/(gr*z)
     if (dabs(dble(x)).gt.1.0d0) then
         ga=ga*r
         if (x.lt.0.0d0) ga=-pi/(x*ga*dsin(pi*x))
     endif
  endif

  ga_out = real(ga)

  return
 end subroutine gamma_reff

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gamma_reff_r8(x,ga_out)
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   compute the gamma function a(x) for single precision floating point. 
!       input :  x  --- argument of a(x)
!                       ( x is not equal to 0,-1,-2,... )
!       output:  ga --- a(x)
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Adapted to SDSU
!           
! References: 
!----------------------------------------------------------------------------------------------------
 implicit double precision (a-h,o-z)
 dimension g(26)
 data g/1.0d0,0.5772156649015329d0, &
        -0.6558780715202538d0, -0.420026350340952d-1, &
        0.1665386113822915d0,-.421977345555443d-1, &
        -.96219715278770d-2, .72189432466630d-2, &
        -.11651675918591d-2, -.2152416741149d-3, &
        .1280502823882d-3, -.201348547807d-4, &
        -.12504934821d-5, .11330272320d-5, &
        -.2056338417d-6, .61160950d-8, &
        .50020075d-8, -.11812746d-8, &
        .1043427d-9, .77823d-11, &
      -.36968d-11, .51d-12, &
      -.206d-13, -.54d-14, .14d-14, .1d-15/
 real(8) ::  x, ga_out
 pi=3.141592653589793d0
 if (x.eq.int(x)) then
    if (x.gt.0.0d0) then
       ga=1.0d0
       m1=x-1
      do k=2,m1
         ga=ga*k
      enddo
    else
      ga=1.0d+300
    endif
 else

    if (dabs(x).gt.1.0d0) then
        z=dabs(x)
        m=int(z)
        r=1.0d0
        do k=1,m
           r=r*(z-k)
        enddo
        z=z-m
    else
        z=x
    endif
    gr=g(26)
    do k=25,1,-1
       gr=gr*z+g(k)
    enddo
    ga=1.0d0/(gr*z)
    if (dabs(x).gt.1.0d0) then
        ga=ga*r
        if (x.lt.0.0d0) ga=-pi/(x*ga*dsin(pi*x))
    endif
 endif

 ga_out = ga

 return
 end subroutine gamma_reff_r8

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 function gammln_reff(xx)
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   compute the gamma function 
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Adapted to SDSU
!           
! References: 
!----------------------------------------------------------------------------------------------------
 real*8 :: cof(6)
 data cof, stp/76.18009173d0, -86.50532033d0, 24.01409822d0,  &
      -1.231739516d0, .120858003d-2, -.536382d-5, 2.50662827465d0/
 data half, one, fpf/0.5d0,  1.0d0,  5.5d0/

 x=xx-one
 tmp=x+fpf
 tmp=(x+half)*log(tmp)-tmp
 ser=one
 do j=1,6
    x=x+one
    ser=ser+cof(j)/x
 enddo
 gammln_reff=tmp+log(stp*ser)
 return
 end function gammln_reff

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gamma_function(gnu,d,dn, fgamma)

!---------------------------------------------------------------------------------------------------
! Comments:                             
!   compute the gamma function for single-precision floating pint
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Adapted to SDSU
!           
! References: 
!----------------------------------------------------------------------------------------------------

 implicit none
 real(4) ,intent(in)  :: gnu    ! PSD shape parameter for generalized gamma distribution
 real(4) ,intent(in)  :: d      ! diameter [m] 
 real(4) ,intent(in)  :: dn     ! characteristic diameter [m]
 real(4) ,intent(out) :: fgamma ! gamma function [1/m]

 real(4) :: gfac

 call gamma_reff(gnu,gfac)

 fgamma = (1.e0/gfac) * ( (d/dn)**(gnu-1.e0) ) * (1.e0/dn) * exp(-d/dn)  ![1/m]

 return
 end subroutine gamma_function

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine gamma_function_r8(gnu,d,dn, fgamma)
!---------------------------------------------------------------------------------------------------
! Comments:                             
!   compute the gamma function for double-precision floating pint
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Adapted to SDSU
!           
! References: 
!----------------------------------------------------------------------------------------------------
 implicit none
 real(8)   ,intent(in) :: gnu    ! PSD shape parameter for generalized gamma distribution
 real(8)   ,intent(in) :: d      ! diameter [m] 
 real(8)   ,intent(in) :: dn     ! characteristic diameter [m]
 real(8)   ,intent(out):: fgamma ! gamma function [1/m]
 real(8) :: gfac

 call gamma_reff_r8(gnu,gfac)

 fgamma = (1.d0/gfac) * ( (d/dn)**(gnu-1.d0) ) * (1.d0/dn) * exp(-d/dn)  ![1/m]

 return
 end subroutine gamma_function_r8

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

