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
  
  if(ellin_min.gt.ellin_max)then
    write(*,*)"Bad initial ellin!"
    stop
  end if

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
  
  integer(kind=8):: ird
  real(kind=8):: OmgK, rdp, ell, hd, W, Wg, Wr, T, Sigp, aT4
  real(kind=8):: S1, S2, S, F1, F2, F, temp1, heigp
  integer(kind=8):: i, n
  
  ird=1
  rdp=rout
  OmgK=1.0d0/dsqrt(rdp)/(rdp-2.0d0)
! solution for S-S disk  
  Sigp=1.4d5 * alpha**(-0.8d0) * mbh**(0.2d0) * (10.0d0*mdot)**(0.7d0)          &
      * (rdp/2.0d0)**(-0.75d0) /(Medd/c/Rg)
  T=6.9d7 * alpha**(-0.2d0) * mbh**(-0.2d0) * (10.0d0*mdot)**0.3d0              &
      * (rdp/2.0d0)**(-0.75d0) 
  aT4= ab*T**4.0d0/3.0d0 /(Medd*c/Rg/Rg)
  Wg=Sigp*kb*T/muave/mp /c/c
  heigp=dsqrt(Wg/sigp)/OmgK
  Wr=2.0d0*heigp * aT4
  W=Wr+Wg
  ell=ellin + 2.0*PI*alpha*rdp*rdp *W
  
  write(*,*)T
!========================================================================
! initial guess
!  Wg=Wl
!  Sigp=Wg*muave*mp/(kb*T) * c*c
!  heigp=dsqrt(Wg/sigp)/OmgK
!  Wr=2.0d0*heigp * aT4
!  W=Wr+Wg
!==========================================================================  
  Rd(ird)=rdp
  Omgd(ird)=ell/rdp/rdp
  Td(ird)=T
  Sigd(ird)=Sigp
  Wgd(ird)=Wg
  Wtd(ird)=W
  Heigd(ird)=Heigp
  Wrd(ird)=Wr
  write(*,*)Heigd(ird), Wgd(ird), Wrd(ird), Sigp
END SUBROUTINE BOUNDARY
