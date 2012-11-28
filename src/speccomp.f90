!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================! 

! To calculate the comptonization.
!
! all the formula can be found in paper: 
!     Coppi & Blandford 1990 MNRAS 245,453 (CB)
!
subroutine compn(nrd)
  use const
  use specvars
  implicit none
  integer(kind=8) nrd
  real(kind=8), external:: distribution,elect_dist
  real(kind=8) omig_in,omig_out,gam,theta,temp_gamm,temp_omig,            &
              taoes,numb_dens_omig_in(nleg_omig)
  integer(kind=8) inu,iomig,igam
      
  theta=kb*tecd(nrd)/mec2
  taoes=(hcd(nrd)*rhocd(nrd)*Medd/c/Rg)/mue/mp*sigmat
      
  do iomig=1,nleg_omig
     numb_dens_omig_in(iomig)=fnu_unred(iomig)/(h*nu_unred(iomig))
  end do
           
  do 50 inu=1,nnu
    omig_out=nu_unred(inu)*h/mec2
    temp_gamm=0.0d0
      
    do 40 igam=nleg_gamm,1,-1
      gam=10.0d0**(xleg_gamm(igam))
      numb_elec(igam)=elect_dist(gam,theta) 
!      if((gam-1.0).lt.omig_out)goto 40  
          
      temp_omig=0.0d0
      do 30 iomig=1,nleg_omig
      
      omig_in=10.0**(xleg_omig(iomig))
      
      if(omig_in.gt.omig_out)goto 30
      
      temp_omig=temp_omig+wleg_omig(iomig)                          &
               *distribution(omig_out,iomig,igam)                   &
               *rate(iomig,igam)**2                                 &
               *omig_in*numb_dens_omig_in(iomig)
30    continue
      temp_gamm=temp_gamm+wleg_gamm(igam)*numb_elec(igam)           &
                         *gam*taoes*temp_omig
40    continue

      flux_comp(inu)=temp_gamm*dlog(10.d0)**2*h*nu_unred(inu)
!     &                        *(1.0+taoes)
50    continue
end subroutine compn


!\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\           
subroutine rate_numb()
      use specvars
      implicit none
            
      real(kind=8), external:: ratecal, omig_mean,omig2_mean
      real(kind=8) omi,omi2
      
      real(kind=8) omig_in,gam,temp
      integer(kind=8) iomig,igam 
      
      write(*,*)"call rate_numb(nrd)"     
      do igam=1,nleg_gamm
      write(*,*)"igam=",igam
      gam=10.0d0**(xleg_gamm(igam))
          
      do iomig=1,nleg_omig
      omig_in=10.0d0**(xleg_omig(iomig))
      rate(iomig,igam)=ratecal(omig_in,gam)
      omi=omig_mean(omig_in,gam)
      omi2=omig2_mean(omig_in,gam)
      temp=dsqrt(3.0*(omi2-omi**2)) 
      omigm(iomig,igam)=omi
      domi(iomig,igam)=min(temp,omi)    
      end do
      end do
                 
end subroutine rate_numb

      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! To cal the electron population
function elect_dist(gam,theta)
      use const
      implicit none
      real*8 elect_dist,gam,theta
      real*8 beta_gam,bessk
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      if(theta.lt.1.0d-2)then
      elect_dist=gam**2*beta_gam*dexp(-(gam-1.0d0)/theta)             &
           /dsqrt(PI/2.0d0)/theta**1.5d0    
      else
      elect_dist=gam**2*beta_gam*dexp(-gam/theta)                     &
                /theta/bessk(2,1.0/theta)
      end if
      return
end function elect_dist
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To cal the scatted photon distribution
! see equation 2.7 in CB
function distribution(omig_out,iomig,igam)
      use specvars
      implicit none    
      real(kind=8) distribution,omig_out
      integer(kind=8) iomig,igam
      real(kind=8) delta_omig,omi,disp
      
      omi=omigm(iomig,igam)
      disp=domi(iomig,igam)
!      write(*,*)disp,dabs(omig_out-omi)      
      if(disp.gt.dabs(omig_out-omi))then
      distribution=0.5d0/disp
      else
      distribution=0.0d0
      end if     
end function distribution
! above is to cal the scattederd photon
! distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
!################################################
! To cal the second momentum of the mean scattered
! photon energy.
!///////////////////////
! function omig2_mean 
! see equation 2.16 in CB 
!//////////////////////    
function omig2_mean(omiga,gam)
      use specvars
      implicit none    
      real(kind=8) omig2_mean,omiga,gam
      real(kind=8) beta_gam,omig2_alpha,ratecal
      integer(kind=8) i
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      omig2_mean=0.0d0
      do i=1,nleg/2
      omig2_mean=omig2_mean                                               &
     +wleg(i)*omig2_alpha(omiga,gam,xleg(i))*(1.0-beta_gam*xleg(i))/2.0   &
     +wleg(i)*omig2_alpha(omiga,gam,xleg(nleg-i+1))                       &
     *(1.0-beta_gam*xleg(nleg+1-i))/2.0
      end do
!      omig2_mean=omig2_mean/ratecal(omiga,gam)
      return
end function omig2_mean
!///////////////////////
! function omig2_alpha 
! see equation 2.16 in CB 
!//////////////////////   
function omig2_alpha(omiga,gam,mu)
      use specvars
      implicit none    
      real(kind=8) omiga,gam,omig2_alpha,mu
      real(kind=8) beta_gam
      real(kind=8), external:: fun_omig2
      
      integer(kind=8) i

      beta_gam=dsqrt(1.0-1.0/gam/gam)
      omig2_alpha=0.0d0
      do i=1,nleg/2
      omig2_alpha=omig2_alpha                               &
      +wleg(i)*fun_omig2(xleg(i),omiga,gam,mu)              &
      +wleg(i)*fun_omig2(xleg(nleg-i+1),omiga,gam,mu)       
      end do
      
      omig2_alpha=omig2_alpha*3.0d0/8.0
      return      
end function omig2_alpha
!///////////////////////
! function fun_omig2
! see equation 2.16 in CB 
!//////////////////////                         
function fun_omig2(t,omiga,gam,mu)
      implicit none
      real*8 fun_omig2,t,omiga,gam,mu
      real*8 xx,xr,beta_gam,omig2averg
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      xx=gam*omiga*(1.0-beta_gam*mu)
      xr=1.0+(1.0-t)*xx
      
      omig2averg=gam**2*omiga**2*                                 &
      (gam**2*(1.0-beta_gam*mu+beta_gam*t*(mu-beta_gam))**2       &     
      +0.5*beta_gam**2*(1-t**2)*(1-mu**2))/xr**2               
      fun_omig2=omig2averg/xr**2*(xr+1.0/xr-1.0+t**2)
      return      
end function fun_omig2
! above is to cal the second moment of the mean 
! scattered
!################################################

!**********************************************
! To cal the mean scattered photon energy.
!///////////////////////
! function omig_mean 
! see equation 2.9 in CB 
!//////////////////////    
function omig_mean(omiga,gam)
      use specvars
      implicit none     
      real*8 omig_mean,omiga,gam

      real*8 beta_gam,omig_alpha,ratecal
      integer i
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      omig_mean=0.0d0
      do i=1,nleg/2
      omig_mean=omig_mean                                                  &
      +wleg(i)*omig_alpha(omiga,gam,xleg(i))*(1.0-beta_gam*xleg(i))/2.0    &
      +wleg(i)*omig_alpha(omiga,gam,xleg(nleg-i+1))                        &
      *(1.0-beta_gam*xleg(nleg+1-i))/2.0
      end do
!      omig_mean=omig_mean/ratecal(omiga,gam)
      return
end function omig_mean
!///////////////////////
! function omig_alpha 
! see equation 2.9 in CB 
!//////////////////////   
function omig_alpha(omiga,gam,mu)
      use specvars
      implicit none     
      real*8 omiga,gam,omig_alpha,mu
      real*8 beta_gam
      real*8 fun_omig
      
      integer i

      beta_gam=dsqrt(1.0-1.0/gam/gam)
      omig_alpha=0.0d0
      do i=1,nleg/2
      omig_alpha=omig_alpha                                     &
         +wleg(i)*fun_omig(xleg(i),omiga,gam,mu)                &
         +wleg(i)*fun_omig(xleg(nleg-i+1),omiga,gam,mu)
      end do
      
      omig_alpha=omig_alpha*3.0d0/8.0
      return      
end function omig_alpha
!///////////////////////
! function fun_omig 
! see equation 2.10 in CB
!//////////////////////                         
function fun_omig(t,omiga,gam,mu)
      implicit none
      real*8 fun_omig,t,omiga,gam,mu
      real*8 xx,xr,beta_gam,omigaverg
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      xx=gam*omiga*(1.0-beta_gam*mu)
      xr=1.0+(1.0-t)*xx
      omigaverg=gam**2*omiga*(1.0-beta_gam*mu+beta_gam*t        &
               *(mu-beta_gam))/xr
      fun_omig=omigaverg/xr**2*(xr+1.0/xr-1.0+t**2)
      return      
end function fun_omig
! above is to cal the mean scattered photon energy.
!*********************************************************
      
!==========================================
! To cal the angle-averaged scattering rate 
!///////////////////////
! function ratecal 
! see equation (2.3) in CB
! the unit is c*sigmat
!//////////////////////       
function ratecal(omiga,gam)
      use specvars
      implicit none   
      real*8 ratecal,omiga,gam
!      integer nleg
!      parameter(nleg=400)
!      real*8 xleg(nleg),wleg(nleg)
!      common/intleg/xleg,wleg,nleg
      real*8 beta_gam,section
      integer i
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      ratecal=0.0d0
      do i=1,nleg/2
      ratecal=ratecal                                                   &
       +wleg(i)*section(omiga,gam,xleg(i))*(1.0-beta_gam*xleg(i))/2.0   &
       +wleg(i)*section(omiga,gam,xleg(nleg-i+1))                       &
       *(1.0-beta_gam*xleg(nleg+1-i))/2.0
      end do
      return
end function ratecal
!///////////////////////
! function section  
! see equation (2.9) in CB
!//////////////////////   
function section(omiga,gam,mu)
      use specvars
      implicit none      
      real*8 omiga,gam,section,mu
      real*8 beta_gam
      real*8 x,funsec
!      integer nleg
!      parameter(nleg=400)
!      real*8 xleg(nleg),wleg(nleg)
!      common/intleg/xleg,wleg,nleg
      
      integer i,n
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      x=omiga*gam*(1.0-beta_gam*mu)

      section=0.0d0
      do i=1,nleg/2
      section=section                               &
      +wleg(i)*funsec(xleg(i),x)                    &
      +wleg(i)*funsec(xleg(nleg-i+1),x)             
      end do
      section=section*3.0d0/8.0
      return      
end function section
!///////////////////////
! function funsec 
! see equation (2.11) in CB 
!//////////////////////                         
function funsec(t,xx)
      implicit none
      real*8 funsec,t,xx
      real*8 xr
      xr=1.0+(1.0-t)*xx
      funsec=1.0/xr**2*(xr+1.0/xr-1.0+t**2)
      return      
end function funsec
! above is to cal the angle-averaged scattering rate
!====================================================
