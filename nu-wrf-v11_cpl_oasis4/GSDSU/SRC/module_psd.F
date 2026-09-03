#include <define_CPP.h>
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 module module_psd
 use module_floating_point

 implicit none

!----------------------------------------------------------------------------------------------------
!              = Goddard Satellite Data Simulator Unit =
!
! NASA GSFC makes no representations about the suitability of software for any purpose. 
! It is provided as is without express or implied warranty. Neither NASA GSFC (the US 
! government) nor Principal Developers (their organizations) shall be liable for any 
! damages suffered by the user of this software. In addition, please do not distribute 
! the software to third party.
!
! Comments:
!          Module to handle some PSD functions. 
!
! History:
! 04/2012  Toshi Matsui@NASA GSFC ; Initial.
!
! References: 
!----------------------------------------------------------------------------------------------------


 save     ! all module parameters will be saved

!
! Encapsulation control 
!
 private   ! ALL variables and subourtines are non-accessible, except....

 public :: get_params_bisect 
 public :: gamma_parameters


!
! gamma parameters
!
 type gamma_parameters
   real :: mu     ! mu        [-]
   real :: lambda ! slope     [1/m]
   real :: n0     ! intercept [1/m3]
   real :: nw     ! intercept [1/m3]
   real :: d0     ! the median volume drop diameter (i.e. D0=(mu+3.67)/lambda)  [mm]
 end type gamma_parameters


 real, parameter :: mu_min=0., mu_max=80., tolerance=0.01
 real,allocatable,dimension(:) :: xb1,xb2

 contains

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine get_params_bisect(m0,m3,m6,output)
 implicit none
 real,intent(in) :: m0,m3,m6 !0th 3rd 6th moments
 type(gamma_parameters), intent(out) :: output
 real   :: f1,f2
 real   :: root,smu
 integer:: imu,i,j,testv
 real*8 :: ag6, ag3, ag0, ag4
! real, external :: gamma
 real :: x1,x2
 real :: f_mu
 integer :: nroot
 real*8 :: sup_n0
 real,allocatable,dimension(:) :: mroot
 real*8 :: imsiaaa,imsibbb,imsiccc,imsiddd,imsieee
 real :: gfac3
 real*8 :: term_dble

!      write(*,*)  'this program can solves equation f(x) = 0'
!      write(*,*)  'please enter two values left and right such that '
!      write(*,*)  'f(left) and f(right) have opposite signs.'
!      write(*,*)
!      write(*,*)  'left and right please --> '
!      read(*,*)   left, right              ! read in left and right
!      call zbrac(x1,x2)
 call zbrak(m0,m3,m6,nroot)

 if(nroot.eq.0) then
    allocate(mroot(1))
    f1=funct(mu_min,m0,m3,m6)
    f2=funct(mu_max,m0,m3,m6)
    if(f2.le.f1) then
       mroot(1)=mu_max
    else
       mroot(1)=mu_min
    endif
!    write(44,fmt='(5i10,40f10.2)') xp,yp,zp,spp,nroot,(mroot(i),i=1,1)
 else
    allocate(mroot(nroot))
    do i=1,nroot
        mroot(i)=rtbis(xb1(i),xb2(i),m0,m3,m6)
    enddo
!    write(44,fmt='(5i10,40f10.2)') xp,yp,zp,spp,nroot,(mroot(i),i=1,nroot)
 endif

 deallocate(xb1,xb2)

 output%mu=mroot(1)

 if( output%mu >= 10.) then
     output%mu = -999.
     return
 endif
 

 smu=mroot(1)
 imu=nint(smu)
 ag6=ngamma(dble(imu+6+1))
 ag3=ngamma(dble(imu+3+1))
 ag0=ngamma(dble(imu+0+1))
 ag4=ngamma(dble(imu+0+4))


 print*,'mu',smu,imu,ag0,ag3,ag6

 !lambda(xp,yp,zp,spp)=((m0(xp,yp,zp,spp)*real(ag6))/(m6(xp,yp,zp,spp)*real(ag0)))**(1./6.)

 output%lambda=(m0/m6)**(1./6.)*(real(ag6/ag0))**(1./6.) !unit(1/m)
 output%lambda=(m0/m6)**(1./6.)*(real(ag6/ag0))**(1./6.) !unit(1/mm)


 !lambda(xp,yp,zp,spp)=(m3(xp,yp,zp,spp)/m6(xp,yp,zp,spp))**(1./3.)*(real(ag6)/real(ag3))**(1./3.)
 !imsiaaa = (lambda(xp,yp,zp,spp)**(imu+3+1))/ag3

 imsiaaa = output%lambda  ! slope (1/m)
if( output%lambda == 0.) print*, 'lambda is zero', (m0/m6)**(1./6.),(real(ag6)/real(ag0))**(1./6.)

! imsibbb = imsiaaa**(imu+4.)/ag3
! imsiccc = m3
! sup_n0=1./(10.**(3*imu+3))
! imsiddd = imsibbb * imsiccc * sup_n0  ! n0
! imsieee = imsiddd*ag0/(imsiaaa*0.001) ! nt
! output%n0=imsieee ! unit=1/m3 mm(mu+1)

! term_dble = dble((output%lambda**(output%mu+3.+1.)))

 output%n0 = m3 * real( dble(output%lambda**(output%mu+3.+1.))    /ag3)    ! [1/m3]


print*, 'n0=',output%n0, m3 , (output%lambda**(output%mu+3.+1.)) , real(ag3)
 call gamma_reff(real(imu)+3.e0+1.e0,gfac3)
print*,'gfac=',gfac3,imu, output%mu


 output%d0=(imu+3.67)/output%lambda*1000. !median diameter   [mm]

      !write(88,*) imsiaaa,imsiddd,imsieee,d0(xp,yp,zp,spp)


 f_mu = 6. * ( (3.67 + output%mu)**(output%mu+4.) )  / ( (3.67**4.) * ag4 )

 output%nw = output%n0 * output%d0 / f_mu

 deallocate(mroot)

 end subroutine get_params_bisect

!rDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 subroutine zbrak(m0,m3,m6,nb)
 implicit none
 real,intent(in) :: m0,m3,m6 !0th 3rd 6th moments
 real :: x1,x2,dx
 integer :: nb, i
 integer, parameter :: n=100
 real,dimension(0:n) :: f,x
 logical,dimension(1:n) :: mask

 x1=mu_min; x2=mu_max
 dx=(x2-x1)/n
 do i=0,n
     x(i)=x1+dx*i
     f(i)=funct(x(i),m0,m3,m6)
 enddo

 mask=f(1:n)*f(0:n-1) <= 0.0
 nb=count(mask)
      !if(associated(xb1)) deallocate(xb1)
      !if(associated(xb2)) deallocate(xb2)
 if(nb.eq.0) then
    allocate(xb1(1),xb2(1))
    !f1=funct(mu_min)
    !f2=funct(mu_max)
    !if(f2.le.f1) then
    xb1(1)=-999.
    !else
    xb2(1)=-999.
        !endif
    else
      allocate(xb1(nb),xb2(nb))
      xb1(1:nb)=pack(x(0:n-1),mask)
      xb2(1:nb)=pack(x(1:n),mask)
    endif
 end subroutine

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 function  funct(x,m0,m3,m6) result(value)
 implicit  none
 real, intent(in) :: x
 real, intent(in) :: m0, m3, m6  !0th 3rd 6th moments 
 real             :: value
 real             :: term1,term2
 if(x==-1.) stop 'MSG funct in module_psd; input x must not be -1'
 term1=(((x+6.)*(x+5.)*(x+4.))/((x+3.)*(x+2.)*(x+1.)))**3.  !x must be non -1
 term2=m0*m6/m3**2.
 value = term1-term2
 !print *, "values :", first, second, first-second,x
 end function  funct

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 function rtbis(x11,x22,m0,m3,m6) result(value)
 implicit none
 real,intent(in) :: x11,x22
 real,intent(in) :: m0, m3, m6  !0th 3rd 6th moments 
 integer, parameter :: jmax=40
 real :: value,fmid,f,dx,xmid !,x11,x22
 integer :: j

 fmid=funct(x22,m0,m3,m6)
 f=funct(x11,m0,m3,m6)
 if(f.lt.0) then
    value=x11
    dx=x22-x11
 else
    value=x22
    dx=x11-x22
 endif

 do j=1,jmax
     dx=dx/2
     xmid=value+dx
     fmid=funct(xmid,m0,m3,m6)
     if(fmid.le.0.) value=xmid
     if(abs(dx).lt.tolerance .or. fmid.eq.0.) return
 enddo

 end function rtbis

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

  function nGamma(xx) result(value)
!*******************************************
!*           FUNCTION  GAMMA(X)            *
!* --------------------------------------- *
!* Returns the value of Gamma(x) in double *
!* precision as EXP(LN(GAMMA(X))) for X>0. *
!*******************************************

  real*8 :: value
  real*8 xx,cof(6),stp,half,one,fpf,x,tmp,ser
  DATA cof,stp /76.18009173d0,-86.50532033d0,24.01409822d0,  &
       -1.231739516d0,0.120858003d-2,-0.536382d-5,2.50662827465d0/
  DATA half,one,fpf /0.5d0,1.0d0,5.5d0/
 integer :: j
  x=xx-one
  tmp=x+fpf
  tmp=(x+half)*DLOG(tmp)-tmp
  ser=one
  do j=1,6
    x=x+one
    ser=ser+cof(j)/x
  end do
  !Gamma = DEXP(tmp+DLOG(stp*ser))
  value = DEXP(tmp+DLOG(stp*ser))
end function nGamma

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

 end module module_psd

!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 
!SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU SDSU 

