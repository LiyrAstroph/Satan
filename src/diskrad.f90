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
  real(kind=8) rdp, sigp, rhop, heigp, Tp, Teffp, nup, Bnu, x, kappaes, kappaff,    &
               temp1, Inu, drlog
  integer(kind=8) i, j
  
  open(unit=10, file='../data/sed.dat')

  kappaes=0.4d0
  temp1=2.0d0*h/c/c
  diskfnu=0.0d0
  do i=1, nt-1
    rdp=Rd(i)
    if(rdp.le.1.0d3.and.rdp.ge.6.0d0)then
      drlog=dlog(Rd(i))-dlog(Rd(i+1))
      Teffp=Teffd(i)
      Tp=Td(i)
      sigp=Sigd(i)
      heigp=Heigd(i)
      rhop=sigp/heigp * (Medd/c/Rg/Rg)
      do j=1, nnu
        nup=nu(j)
        x=h*nup/kb/Teffp
        if(x.ge.1.0d-6)then
        Bnu=temp1*nup**3.0d0 / (dexp(x)-1.0d0)
        kappaff=1.5d25*rhop/Teffp**(3.5d0) / x**3.0d0 *(1.0d0-dexp(-x))
        else
        Bnu=temp1*nup**3.0d0 / (x+0.5d0*x*x)
        kappaff=1.5d25*rhop/Teffp**(3.5d0) / x**3.0d0 *(x-0.5d0*x*x) 
        end if
        Inu=2.0d0/(1.0d0+dsqrt(1.0d0+kappaes/kappaff)) * Bnu
        diskfnu(j)=diskfnu(j) + Inu * 2.0d0*PI*rdp*rdp*drlog
      end do
    endif
  enddo
  do i=1,nnu
    write(10, *)nu(i), (2.0d0*PI*diskfnu(i)*Rg*Rg+eps)*nu(i)
  end do
  close(10)
END SUBROUTINE DISKSPEC

!*****************************************************************************!
! cal the local flux given temperature and density
! input: tep--temperature, teffp--effective temperature
!        sigp--surface density,  heigp--height
! output: fd
!*****************************************************************************!
SUBROUTINE DISKRAD(tep, teffp, sigp, heigp, fd)
  use const
  use diskvars
  implicit none
  real(kind=8) tep, teffp, sigp, heigp, fd(nnu_max)
  real(kind=8) rhop, nup, x, Bnu, Inu, kappaff, kappaes, temp1
  integer(kind=8) j

  rhop=sigp/heigp * (Medd/c/Rg/Rg)
  kappaes=0.4d0
  temp1=2.0d0*h/c/c
  do j=1, nnu
    nup=nu(j)
    x=h*nup/kb/Teffp
    if(x.ge.1.0d-6)then
      Bnu=temp1*nup**3.0d0 / (dexp(x)-1.0d0)
      kappaff=1.5d25*rhop/Teffp**(3.5d0) / x**3.0d0 *(1.0d0-dexp(-x))
    else
      Bnu=temp1*nup**3.0d0 / (x+0.5d0*x*x)
      kappaff=1.5d25*rhop/Teffp**(3.5d0) / x**3.0d0 *(x-0.5d0*x*x) 
    end if
    Inu=2.0d0/(1.0d0+dsqrt(1.0d0+kappaes/kappaff)) * Bnu
    fd(j)= Inu * PI
   end do
END SUBROUTINE DISKRAD