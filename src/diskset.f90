!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

SUBROUTINE INIT
  use const
  use diskvars
  implicit none
  open(unit=10, file="init.dat", status="unknown")
  
  read(*,slim)
  write(10,slim)
  read(*,mesh)
  write(10,mesh)
  
  Ledd=mbh * 4.0d0*PI*G*Msun*mp*c/sigmaT
  Medd=Ledd/0.1d0/c/c
  Rg=mbh * G*Msun/c/c
  
  muave=1.0d0/mui + 1.0d0/mue
  muave=1.0d0/muave
  
  write(10,"(a, e10.2)")"Ledd=	", Ledd
  write(10,"(a, e10.2)")"Medd=	", Medd
  write(10,"(a, e10.2)")"Rg=	", Rg
  write(10,"(a, e10.2)")"muave=	", muave
  close(10)
END SUBROUTINE INIT

SUBROUTINE BOUNDARY
  use const
  use diskvars
  implicit none
  
  integer(kind=8):: ir
  real(kind=8):: OmgK, r, ell, hd, W, Wg, Wr, T, Sigp, aT4, Wl
  real(kind=8):: S1, S2, S, F1, F2, F, temp1, heigp
  integer(kind=8):: i, n
  
  ir=1
  r=rout
  OmgK=1.0d0/dsqrt(r)/(r-2.0d0)
  Omgd(ir)=OmgK/dsqrt(5.0d0)
  Rd(ir)=r
  ell=r*r*Omgd(ir)
  T=1.37d8*(mdot/16.0d0/mbh)**0.25d0 * (10.0d0/r)**(5.0d0/8.0d0) 
  aT4= ab*T**4.0d0/3.0d0 /(Medd*c/Rg/Rg)
  Wl=mdot*(ell-ellin)/(2.0*PI*alpha*r*r)
  
  write(*,*)T
!========================================================================
! initial guess
  Wg=Wl
  Sigp=Wg*muave*mp/(kb*T) * c*c
  heigp=dsqrt(Wg/sigp)/OmgK
  Wr=2.0d0*heigp * aT4
  W=Wr+Wg
!==========================================================================  
  Td(ir)=T
  Sigd(ir)=Sigp
  Wgd(ir)=Wg
  Wtd(ir)=W
  Heigd(ir)=Heigp
  Wrd(ir)=Wr
  write(*,*)Heigd(ir), Wgd(ir), Wrd(ir), Sigp
END SUBROUTINE BOUNDARY
