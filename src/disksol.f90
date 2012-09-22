!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

SUBROUTINE SOLVE
  use const
  use diskvars
  implicit none 
! for rk       
  real(kind=8) k1(2),k2(2),k3(2),k4(2),k5(2),k6(2),x,                        &
       y1,y2,x0,y10,y20,yerr1,yerr2, errmax,htemp,                           &                               
       a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,                               &
       b51,b52,b53,b54,b61,b62,b63,b64,b65,                                  &
       c1,c2,c3,c4,c5,c6,dc1,dc2,dc3,dc4,dc5,dc6
  parameter(a2=0.2d0,a3=0.3d0,a4=0.6d0,a5=1.0d0,a6=0.875d0,                  &
       b21=0.2d0,b31=0.075d0,b32=0.225d0, b41=0.3d0,b42=-0.9d0, b43=1.2d0,   &
       b51=-11.0d0/54.0d0,b52=2.5d0,b53=-70.0d0/27.0d0,                      &
       b54=35.0d0/27.0d0,b61=1631.0d0/55296.0d0,                             &
       b62=175.0d0/512.0d0,b63=575.0d0/13824.0d0,                            &
       b64=44275.0d0/110592.0d0,b65=253.0d0/4096.0d0,                        &
       c1=37.0d0/378.0d0,c2=0.0d0,c3=250.0d0/621.0d0,                        &
       c4=125.0d0/594.0d0,c5=0.0d0,c6=512.0d0/1771.0d0,                      &
       dc1=c1-2825.0d0/27648.0d0,dc2=0.0d0,                                  &
       dc3=c3-18575.0d0/48384.0d0,dc4=c4-13525.0d0/55296.0d0,                &
       dc5=-277.0d0/14336.0d0,dc6=c6-0.25d0)     
  real(kind=8) safety,pgrow,pshrnk,errcon,eacc
  parameter(safety=0.9d0,pgrow=-0.20d0,pshrnk=-0.25d0,errcon=1.89d-4,        &
       eacc=1.0d-4)  
!===============================================================================
  real(kind=8) hstep, rdp, sigp, wtp, df(2), df0(2), vr, cs, vr2cs, factor
  integer(kind=8) ird, irdi, irde, ird0, i, ns
  logical istrs, flag

  istrs=.false.
  flag=.false.

20  irdi=1
  irde=nd
  hstep=-drstep_min
  rdp=Rd(irdi)
  irdi=irdi+1
  
  write(*,*)"Start!"
  ellin=0.5d0*(ellin_min+ellin_max) 
  write(*,*)ellin_min, ellin_max, ellin
  do ird=irdi, irde, 1
    if(rdp.le.rin)exit
    ird0=ird-1
   
! rk step 1
40    rdp=Rd(ird0)
      sigp=Sigd(ird0)
      wtp=Wtd(ird0)
      
      x0=dlog(rdp)      
      y10=dlog(wtp)
      y20=dlog(sigp)     
       
      call deriva(df,rdp,wtp,sigp,ird0) 
           
      do i=1,2
      k1(i)=hstep*df(i)
      end do  
! rk step 2        
      x=x0+a2*hstep
      y1=y10+b21*k1(1)
      y2=y20+b21*k1(2)
      
      rdp=dexp(x)
      wtp=dexp(y1)
      sigp=dexp(y2)
      
      call deriva(df,rdp,wtp,sigp,ird0) 
      do i=1,2
      k2(i)=hstep*df(i)
      end do      
! rk step 3
      x=x0+a3*hstep
      y1=y10+b31*k1(1)+b32*k2(1)
      y2=y20+b31*k1(2)+b32*k2(2)
      
      rdp=dexp(x)
      wtp=dexp(y1)
      sigp=dexp(y2)
      call deriva(df,rdp,wtp,sigp,ird0)  
      do i=1,2
      k3(i)=hstep*df(i)
      end do 
                   
! rk step 4
      x=x0+a4*hstep
      y1=y10+b41*k1(1)+b42*k2(1)+b43*k3(1)
      y2=y20+b41*k1(2)+b42*k2(2)+b43*k3(2)
      
      rdp=dexp(x)
      wtp=dexp(y1)
      sigp=dexp(y2)
      call deriva(df,rdp,wtp,sigp,ird0) 
      do i=1,2
      k4(i)=hstep*df(i)
      end do 
! rk step 5

      x=x0+a5*hstep
      y1=y10+b51*k1(1)+b52*k2(1)+b53*k3(1)+b54*k4(1)
      y2=y20+b51*k1(2)+b52*k2(2)+b53*k3(2)+b54*k4(2)
      
      rdp=dexp(x)
      wtp=dexp(y1)
      sigp=dexp(y2)
      call deriva(df,rdp,wtp,sigp,ird0) 
      do i=1,2
      k5(i)=hstep*df(i)
      end do 
! rk step 6

      x=x0+a6*hstep
      y1=y10+b61*k1(1)+b62*k2(1)+b63*k3(1)+b64*k4(1)+b65*k5(1)
      y2=y20+b61*k1(2)+b62*k2(2)+b63*k3(2)+b64*k4(2)+b65*k5(2)
      
      rdp=dexp(x)
      wtp=dexp(y1)
      sigp=dexp(y2)
      
      call deriva(df,rdp,wtp,sigp,ird0) 
      do i=1,2
      k6(i)=hstep*df(i)
      end do 

      y1=c1*k1(1)+c3*k3(1)+c4*k4(1)+c6*k6(1)
      y2=c1*k1(2)+c3*k3(2)+c4*k4(2)+c6*k6(2)
      
      df0(1)=y1/hstep
      df0(2)=y2/hstep
           
!      write(90,*)y1,y2,y3      
      yerr1=dc1*k1(1)+dc3*k3(1)+dc4*k4(1)+dc5*k5(1)+dc6*k6(1)                      
      yerr2=dc1*k1(2)+dc3*k3(2)+dc4*k4(2)+dc5*k5(2)+dc6*k6(2)
      
      errmax=0.0d0
      errmax=max(dabs(yerr1),dabs(yerr2))
      
      errmax=errmax/eacc
      if(errmax.gt.1.0d0)then
      htemp=safety*hstep*(errmax**pshrnk)/10.0
      htemp=sign(max(dabs(htemp),0.1*dabs(hstep)),hstep)
      if(dabs(htemp).lt.drstep_min)then
         write(*,*)'   hstep underflow'
         htemp=-drstep_min
!         nstep=ird0
         goto 60
         end if
      hstep=htemp
      goto 40 
      else
         if(errmax.gt.errcon)then
         htemp=safety*hstep*(errmax**pgrow)
         else
         htemp=5.0*hstep
         end if
         if(dabs(htemp).gt.drstep_max)then
         htemp=-drstep_max 
         end if
         
         rdp=rd(ird0)*dexp(hstep)
         rd(ird)=rdp     
         Wtd(ird)=Wtd(ird0)*dexp(y1)
         Sigd(ird)=Sigd(ird0)*dexp(y2)
         dfd(ird)=df0(1)
         call VRCAL(rdp, wtp, sigp, factor)
         vr=mdot/(2.0*PI*rdp*Sigd(ird))
         cs=dsqrt(Wtd(ird)/Sigd(ird))
         vr2cs=vr/cs/factor
!         if(mod(ird, 100).eq.0)write(*,*)ird, rdp
         if(flag)write(12,*)Rd(ird), Wtd(ird), Sigd(ird), vr2cs
         
         if(vr2cs.gt.0.95d0.and.vr2cs.lt.1.05d0)then
           ns=ird0
           istrs=.true.
           call crossonic(ird0,df0,-1.0d-1)
           rdp=rd(ird0)*dexp(-1.0d-1)
           write(*,*)'Transonic!', rdp, vr2cs
         endif
         
         hstep=htemp
         if(rdp*dexp(hstep).lt.rin)then
         hstep=dlog(rin/rdp)
         end if
      endif
  enddo
  if(istrs.and.(ellin_max-ellin_min).gt.1.0d-6)then
    if(dabs((dfd(ns+1)-dfd(ns))/(rd(ns+1)-rd(ns))).gt.10.0d0)then
      ellin_max=ellin
      istrs=.false.
      write(*,*)dabs((dfd(ns+1)-dfd(ns))/(rd(ns+1)-rd(ns)))
      goto 20
    end if
  else
    write(*,*)"No transonic"
    ellin_min=ellin
    istrs=.false.
    goto 20
  endif
  nt=ird-1
  write(*,*)"End!", ird, rdp,dabs((dfd(ns+1)-dfd(ns))/(rd(ns+1)-rd(ns)))
  call output
60 stop  
END SUBROUTINE SOLVE

SUBROUTINE DERIVA(df,rdp,wtp,sigp,ird0)
  use const
  use diskvars 
  implicit none
  real(kind=8) df(2),rdp,sigp,wtp
  integer(kind=8) ird0
  
  real(kind=8), external::gettemp
  real(kind=8) a11, a12, a21, a22, c1, c2, wgp, wrp, Tp, fun, dfun, heigp,          &
      OmgK, beta, gam1, gam3, aleff, Qrad, kappa, rhop, ell, ellk, dOmgk, temp1, temp2
  
  OmgK=1.0d0/dsqrt(rdp) / (rdp-2.0d0)
  ellk=rdp*rdp*OmgK
  dOmgK=-0.5d0/rdp-1.0d0/(rdp-2.0d0)
  heigp=dsqrt(wtp/sigp)/OmgK
  rhop=0.5d0*sigp/heigp
   
  Tp=gettemp(rdp, wtp, sigp)
  wgp=sigp * kb*Tp/muave/mp /c/c
  wrp=max(wtp-wgp, 0.0d0)
  beta=wgp/wtp
  beta=min(beta, 1.0d0)
  aleff=alpha*beta**mu
  gam1=beta + (gam-1.0d0)*(4.0d0-3.0d0*beta)**2.0d0 / (beta + 12.0d0*(gam-1.0d0)*(1.0d0-beta))
  gam3=1.0d0 + (gam-1.0d0)*(4.0d0-3.0d0*beta) / (beta + 12.0d0*(gam-1.0d0)*(1.0d0-beta))
  kappa=0.40d0 + 0.64d23*(rhop*Medd/c/Rg/Rg)/Tp**(3.5d0)
  Qrad=8.0d0*c*(1.0d0-beta)* (wtp/sigp/heigp)/kappa  * Rg/Medd
  ell=ellin + 2.0d0*PI*aleff*wtp*rdp*rdp/mdot
!  write(*,*)wgp, wrp, beta, gam1, gam3, aleff
!  write(*,*)kappa, rhop*Medd/c/Rg/Rg, Qrad
!=================================================  
  a11=wtp/sigp
  a12=-(mdot/(2.0d0*PI*rdp*sigp))**2.0d0
  a21= (gam1+1.0d0)/2.0d0/(gam3-1.0d0) * mdot*wtp/sigp                                &
      -(2.0d0*PI*aleff*rdp*wtp)**2.0d0/mdot * (mu*(1.0d0+beta)/2.0d0/(4.0d0-3.0d0*beta) + (1.0d0-mu))
  a22= (3.0d0*gam1-1.0d0)/2.0d0/(gam3-1.0d0) * mdot*wtp/sigp                                &
      +(2.0d0*PI*aleff*rdp*wtp)**2.0d0 * (9.0d0*mu*(1.0d0-beta)/2.0d0/(4.0d0-3.0d0*beta)) 
  a22=-a22
  
  c1=(ell*ell-ellk*ellk)/rdp**3.0d0  - wtp/sigp * dOmgK + mdot*mdot/(4.0d0*PI*PI*rdp**3.0d0 * sigp*sigp)
  
  temp1=-mdot*(ell-ellin)*ellin/rdp**3.0d0 + 2.0*PI*rdp*Qrad   
  temp2= (gam1-1.0d0)/(gam3-1.0d0)*mdot*wtp/sigp                                         & 
       + (2.0d0*PI*aleff*rdp*wtp)**2.0d0/mdot * mu *(1.0d0-beta)/(4.0d0-3.0d0*beta)
  c2=temp1 + temp2*dOmgK
  
  df(1)=(a22*c1-a12*c2)/(a11*a22-a12*a21) * rdp
  df(2)=(a11*c2-a21*c1)/(a11*a22-a12*a21) * rdp

END SUBROUTINE DERIVA

FUNCTION GETTEMP(rdp, wtp, sigp)
  use const
  use diskvars
  implicit none
  real(kind=8) gettemp, rdp, wtp, sigp
  real(kind=8) fun, Tp, dfun, OmgK, heigp
  integer(kind=8) nloop

  OmgK=1.0d0/dsqrt(rdp) / (rdp-2.0d0)
  heigp=dsqrt(wtp/sigp)/OmgK 
!=================================================
! solve for wg and wr  
  fun=1.0d10
  Tp= ( wtp * 3.0d0 / 2.0d0 /heigp/ab * (Medd*c/Rg/Rg) )**0.25d0
  nloop=0
  do while (dabs(fun).gt.eps)
  fun= sigp * kb*Tp/muave/mp /c/c + 2.0d0*heigp/3.0d0 * ab * Tp**4.0d0 / (Medd*c/Rg/Rg) - wtp
  dfun = sigp*kb/muave/mp  /c/c + 8.0d0*heigp/3.0d0 * ab * Tp**3.0d0 / (Medd*c/Rg/Rg)
  Tp=Tp - fun/dfun
  nloop=nloop+1
  if(nloop.gt.10000)then
!    write(*,*)"nloop=10000!", Tp, fun
    exit
  end if
  enddo 
  gettemp=Tp  
END FUNCTION GETTEMP


SUBROUTINE VRCAL(rdp, wtp, sigp, factor)
  use const
  use diskvars
  implicit none
  real(kind=8) rdp, wtp, sigp, factor
  
  real(kind=8), external::gettemp
  real(kind=8) OmgK, heigp, fun, dfun, Tp, wgp, beta, aleff, gam1, gam3,   &
         b1, b2, b3, b4
  
  OmgK=1.0d0/dsqrt(rdp) / (rdp-2.0d0)
  heigp=dsqrt(wtp/sigp)/OmgK
   
  Tp=gettemp(rdp, wtp, sigp)
  wgp=sigp * kb*Tp/muave/mp /c/c
  beta=wgp/wtp
  beta=min(beta, 1.0d0)
  aleff=alpha*beta**mu
  gam1=beta + (gam-1.0d0)*(4.0d0-3.0d0*beta)**2.0d0 / (beta + 12.0d0*(gam-1.0d0)*(1.0d0-beta))
  gam3=1.0d0 + (gam-1.0d0)*(4.0d0-3.0d0*beta) / (beta + 12.0d0*(gam-1.0d0)*(1.0d0-beta))

  b1=(3.0d0*gam1-1.0d0)/2.0d0/(gam3-1.0d0)
  b2=aleff*aleff*9.0d0*mu*(1.0d0-beta)/2.0d0/(4.0d0-3.0d0*beta)
  b3=(gam1+1.0d0)/2.0d0/(gam3-1.0d0)
  b4=aleff*aleff*(mu*(1.0d0+beta)/2.0d0/(4.0d0-3.0d0*beta) + (1.0d0-mu))
  if(b2.lt.1.0d-10)then
  factor=b3/(b1+b4)
  else
  factor=(-(b1+b4)+dsqrt((b1+b4)**2.0d0+4.0d0*b2*b3))/2.0d0/b2
  endif
  factor=1.0d0/dsqrt(factor)
!  write(*,*)factor, b1, b2, aleff
!================================================= 
END SUBROUTINE VRCAL

SUBROUTINE OUTPUT
  use const
  use diskvars
  implicit none
  integer(kind=8) ird
  
  real(kind=8), external::gettemp
  real(kind=8) rdp, wtp, sigp
  real(kind=8) OmgK, heigp, fun, dfun, Tp, wgp, beta, aleff, gam1, gam3,   &
         Qrad, kappa, rhop, teffp, omgp
  
  do ird=1, nt
  rdp=Rd(ird)
  sigp=Sigd(ird)
  wtp=Wtd(ird)
  
  OmgK=1.0d0/dsqrt(rdp) / (rdp-2.0d0)
  heigp=dsqrt(wtp/sigp)/OmgK
  rhop=0.5d0*sigp/heigp
  Tp=gettemp(rdp, wtp, sigp)
  wgp=sigp * kb*Tp/muave/mp /c/c
  beta=wgp/wtp
  beta=min(beta, 1.0d0)
  aleff=alpha*beta**mu
  gam1=beta + (gam-1.0d0)*(4.0d0-3.0d0*beta)**2.0d0 / (beta + 12.0d0*(gam-1.0d0)*(1.0d0-beta))
  gam3=1.0d0 + (gam-1.0d0)*(4.0d0-3.0d0*beta) / (beta + 12.0d0*(gam-1.0d0)*(1.0d0-beta))
  kappa=0.40d0 + 0.64d23*(rhop*Medd/c/Rg/Rg)/Tp**(3.5d0)
  Qrad=8.0d0*c*(1.0d0-beta)* (wtp/sigp/heigp * c*c/Rg )/kappa
  teffp=(0.5d0*Qrad/sigmab)**0.25d0
  omgp=(2.0d0*PI*aleff*wtp/mdot + ellin)/rdp/rdp
  
  Td(ird)=Tp
  Heigd(ird)=heigp
  Wgd(ird)=sigp * kb*Tp/muave/mp /c/c
  Wrd(ird)=2.0d0*heigp/3.0d0 * ab * Tp**4.0d0 / (Medd*c/Rg/Rg)
  Teff(ird)=teffp
  Omgd(ird)=omgp
  write(14, *) rdp, sigp, Tp, teffp, omgp/Omgk, heigp, beta, aleff
  enddo
END SUBROUTINE OUTPUT


subroutine crossonic(ird0,df,hstep)
  use const
  use diskvars
  implicit none
  integer(kind=8) ird0
  real(kind=8) df(2),hstep
  integer(kind=8) j,ird
  real(kind=8) y1,y2,y3,rdp
      
  y1=df(1)*hstep
  y2=df(2)*hstep

  ird=ird0+1
      
  rdp=rd(ird0)*dexp(hstep)
  Rd(ird)=rdp     
  Wtd(ird)=Wtd(ird0)*dexp(y1)      
  Sigd(ird)=Sigd(ird0)*dexp(y2)    
end subroutine crossonic











