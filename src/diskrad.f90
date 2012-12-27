!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

!******************************************************************************!
! cal the spectrum from the cold disk
!
!******************************************************************************!
SUBROUTINE DISKSPEC
  use const
  use diskvars
  implicit none
  real(kind=8) rdp, sigp, heigp, Tp, Teffp,  Ts, &
               drlog, fd(nnu_max)
  integer(kind=8) i, j
  
  open(unit=10, file='../data/sed.dat')

  diskfnu=0.0d0
  do i=1, nt-1
    rdp=Rd(i)
    if(rdp.le.1.0d3.and.rdp.ge.6.0d0)then
      drlog=dlog(Rd(i-1))-dlog(Rd(i))
      Teffp=Teffd(i)
      Tp=Td(i)
      sigp=Sigd(i)
      heigp=Heigd(i)
      call diskrad3(tp, teffp, sigp, heigp, Ts, fd)
      write(98,*)rdp, Ts, Tp, Teffp
      do j=1, nnu
        diskfnu(j)=diskfnu(j) + fd(j) * 2.0d0*PI*rdp*rdp*drlog
      end do
    endif
  enddo
  do i=1,nnu
    write(10, *)nu(i), 2.0d0*(diskfnu(i)*Rg*Rg+eps)*nu(i)
  end do
  close(10)
END SUBROUTINE DISKSPEC

!*****************************************************************************!
! cal the local flux given temperature and density
! input: tep--temperature, teffp--effective temperature
!        sigp--surface density,  heigp--height
! output: fd
!*****************************************************************************!
SUBROUTINE DISKRAD(tep, teffp, sigp, heigp, Tsur, fd)
  use const
  use diskvars
  implicit none
  real(kind=8) tep, teffp, sigp, heigp, Tsur, fd(nnu_max)
  real(kind=8), external::diskrad_fun
  real(kind=8) Ts, Tsold, Tsd, fun, dfun
  integer(kind=8) nloop
  
  Tsold=Tep
  fun=1.0d10
  nloop=0
  do while (dabs(fun).gt.1.0d-1)
  nloop=nloop+1
  Ts=Tsold
  fun=diskrad_fun(ts, teffp, sigp, heigp, fd)
  Tsd=Ts*1.05
  dfun=diskrad_fun(tsd, teffp, sigp, heigp, fd)
  dfun=(dfun-fun)/(Tsd-Ts)
  Tsold=Ts-0.5*fun/dfun
  Tsold=max(Tsold, Teffp)
  if(nloop.gt.100000)then
    write(*,*)"nloop.gt.10000", fun
    exit
  endif
  enddo
  Tsur=Ts
END SUBROUTINE DISKRAD

SUBROUTINE DISKRAD2(tep, teffp, sigp, heigp, Tsur, fd)
  use const
  use diskvars
  implicit none
  real(kind=8) tep, teffp, sigp, heigp, Tsur, fd(nnu_max)
  real(kind=8) nup, x, rhop, Qrad, kappaff, kappaes, taunu, taueff,   &
               temp1, Ts, Ts1, Ts2, Bnu, Inu, ftot, Fc, Comp, Xie, fff, fth
  integer(kind=8) j
  
  rhop=sigp/heigp * (Medd/c/Rg/Rg)
  Qrad=sigmab*Teffp**4.0d0
  kappaes=0.4d0
  temp1=2.0d0*h/c/c

  Ts1=teffp
  Ts2=tep*100.0
1 Ts=dsqrt(Ts1*Ts2)

  Xie=4.0*kb*Ts/mec2
  Fc=0.0d0
  do j=1, nnu
    nup=nu(j)
    x=h*nup/kb/Ts
    if(x.ge.1.0d-6)then
      Bnu=temp1*nup**3.0d0 / (dexp(x)-1.0d0)
      kappaff=1.5d25*rhop/Ts**(3.5d0) / x**3.0d0 *(1.0d0-dexp(-x))
    else
      Bnu=temp1*nup**3.0d0 / (x+0.5d0*x*x)
      kappaff=1.5d25*rhop/Ts**(3.5d0) / x**3.0d0 *(x-0.5d0*x*x) 
    end if
    taueff=dsqrt((kappaes+kappaff)*kappaff)* sigp * (Medd/c/Rg)
    taunu=(kappaff+kappaes)/(1.0d0+taueff)* sigp * (Medd/c/Rg)
    
    fth=0.0d0
    if(x.lt.1.0d0)then
    if(Xie.lt.1.0d-6)then
      fth=dexp(-dlog(1.0/x)/(taunu*taunu*(Xie + Xie*Xie)))
    else
      fth=dexp(-dlog(1.0/x)/(taunu*taunu*dlog(1.0d0+Xie+Xie*Xie)))
    endif
    endif
    fth=0.0d0

    fff=2.0d0* (1.0d0-dexp(-2.0*taueff)) / (1.0d0+dsqrt(1.0d0+kappaes/kappaff)) 
    Inu= fff*(1.0d0-fth)* PI*Bnu
    Fc=Fc+ fth * fff*PI*Bnu * wleg_nu(j)
    fd(j)= Inu
   end do

   Fc=Fc/h * 3.0d0*kb*Ts * dlog(10.0d0) 

   Comp=Fc/sigmab/Ts**4.0d0
   ftot=0.0d0
   do j=1, nnu
    nup=nu(j)
    x=h*nup/kb/Ts
    if(x.ge.1.0d-6)then
      Bnu=temp1*nup**3.0d0 / (dexp(x)-1.0d0)
    else
      Bnu=temp1*nup**3.0d0 / (x+0.5d0*x*x)
    end if
     fd(j)= fd(j) + PI * Bnu * Comp
     ftot=ftot + fd(j) * nup * dlog(10.0d0) * wleg_nu(j)
   enddo

   if(abs(Ts2-Ts1).lt.1.0d-6)then
     write(99, *)Ts, ftot/Qrad, Teffp, Tep
     Tsur=Ts
     return
   endif

   if(ftot/Qrad.lt.1.0)then
    Ts1=Ts
   else
    Ts2=Ts
   endif
   goto 1

END SUBROUTINE DISKRAD2

SUBROUTINE DISKRAD3(tep, teffp, sigp, heigp, Tsur, fd)
  use const
  use diskvars
  implicit none
  real(kind=8) tep, teffp, sigp, heigp, Tsur, fd(nnu_max)
  real(kind=8) nup, x, Bnu, temp1
  integer(kind=8) j
  Tsur=Teffp
  temp1=2.0d0*h/c/c
  do j=1, nnu
    nup=nu(j)
    x=h*nup/kb/Teffp
    if(x.ge.1.0d-6)then
      Bnu=temp1*nup**3.0d0 / (dexp(x)-1.0d0)
    else
      Bnu=temp1*nup**3.0d0 / (x+0.5d0*x*x)
    end if
     fd(j)= PI * Bnu
  enddo 
END SUBROUTINE DISKRAD3

FUNCTION DISKRAD_FUN(ts, teffp, sigp, heigp, fd)
  use const
  use diskvars
  implicit none
  real(kind=8) diskrad_fun, ts, sigp, teffp, heigp, fd(nnu_max)
  real(kind=8) Qrad, kappaes, kappaff, temp1, taueff,  x, Bnu, nup, rhop, Inu, &
               taunu, Xie, fth, fff, Comp, Fc, ftot
  integer(kind=8) j

  rhop=sigp/heigp * (Medd/c/Rg/Rg)
  Qrad=sigmab*Teffp**4.0d0
  kappaes=0.4d0
  temp1=2.0d0*h/c/c

  Xie=4.0*kb*Ts/mec2
  Fc=0.0d0
  do j=1, nnu
    nup=nu(j)
    x=h*nup/kb/Ts
    if(x.ge.1.0d-6)then
      Bnu=temp1*nup**3.0d0 / (dexp(x)-1.0d0)
      kappaff=1.5d25*rhop/Ts**(3.5d0) / x**3.0d0 *(1.0d0-dexp(-x))
    else
      Bnu=temp1*nup**3.0d0 / (x+0.5d0*x*x)
      kappaff=1.5d25*rhop/Ts**(3.5d0) / x**3.0d0 *(x-0.5d0*x*x) 
    end if
    taueff=dsqrt((kappaes+kappaff)*kappaff)* sigp * (Medd/c/Rg)
    taunu=(kappaff+kappaes)/(1.0d0+taueff)* sigp * (Medd/c/Rg)
    
    fth=0.0d0
    if(x.lt.1.0d0)then
    if(Xie.lt.1.0d-6)then
      fth=dexp(-dlog(1.0/x)/(taunu*taunu*(Xie + Xie*Xie)))
    else
      fth=dexp(-dlog(1.0/x)/(taunu*taunu*dlog(1.0d0+Xie+Xie*Xie)))
    endif
    endif

    fff=2.0d0* (1.0d0-dexp(-2.0*taueff)) / (1.0d0+dsqrt(1.0d0+kappaes/kappaff)) 
    Inu= fff*(1.0d0-fth)* PI*Bnu
    Fc=Fc+ fth * fff*PI*Bnu * wleg_nu(j)
    fd(j)= Inu
   end do

   Fc=Fc/h * 3.0d0*kb*Ts * dlog(10.0d0) 

   Comp=Fc/sigmab/Ts**4.0d0
   ftot=0.0d0
   do j=1, nnu
    nup=nu(j)
    x=h*nup/kb/Ts
    if(x.ge.1.0d-6)then
      Bnu=temp1*nup**3.0d0 / (dexp(x)-1.0d0)
    else
      Bnu=temp1*nup**3.0d0 / (x+0.5d0*x*x)
    end if
     fd(j)= fd(j) + PI * Bnu * Comp
     ftot=ftot + fd(j) * nup * dlog(10.0d0) * wleg_nu(j)
   enddo

   diskrad_fun=ftot/Qrad -1.0d0

END FUNCTION DISKRAD_FUN

!**********************************************************************!
! cal the correction factor
!
!**********************************************************************!
SUBROUTINE CORR_FACTOR
  use const
  use diskvars
  implicit none
  real(kind=8) Lbol, Lopt, nuopt, derr
  integer(kind=8) i, np
  nuopt=c/(5100.0*1.0d-8)

  Lbol=0.0d0
  do i=1, nnu
    Lbol=Lbol+diskfnu(i)*nu(i)*dlog(10.0d0)*wleg_nu(i)
  enddo
  
  call hunt(nu, nnu, nuopt, np)

  call polint(nu(np), diskfnu(np), 2, nuopt, Lopt, derr)

  write(*,*)"correction", Lbol/(Lopt*nuopt) , Lbol* 2.0 *Rg*Rg, Lopt*nuopt* 2.0 *Rg*Rg,dlog10(nuopt)
END SUBROUTINE CORR_FACTOR