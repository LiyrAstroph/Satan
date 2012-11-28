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
  real(kind=8) nulog_min, nulog_max

  open(unit=10, file="../data/init.dat", status="unknown")
  
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
  
  nulog_min=dlog10(numin)
  nulog_max=dlog10(numax)
  call gauleg(nulog_min, nulog_max, xleg_nu, wleg_nu, nnu)
  nu=10.0d0**(xleg_nu)
  
END SUBROUTINE INIT

SUBROUTINE BOUNDARY
  use const
  use diskvars
  implicit none
  
  integer(kind=8):: ird
  real(kind=8):: OmgK, rdp, ell, Wt, Wg, Wr, T, Sigp, aT4
  real(kind=8):: heigp, factor, vr2cs, vr, cs
  
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
  Wt=Wr+Wg
  ell=ellin + 2.0*PI*alpha*rdp*rdp *Wt/mdot
  
  write(*,*)T, sigp, ell/OmgK/rdp/rdp

!  ell=OmgK*rdp*rdp*0.0001
!  Wt=mdot*(ell-ellin)/(2.0d0*PI*alpha*rdp*rdp)
!  T=1.0d3
!  aT4= ab*T**4.0d0/3.0d0 /(Medd*c/Rg/Rg)
!  Sigp=Wt * muave * mp / kb /T * c*c
!  heigp=dsqrt(Wt/sigp)/OmgK
!  Wr=2.0d0*heigp * aT4
!  Wg=Wt
!  Wt=Wr+Wg
!  write(*,*)T, sigp, ell/OmgK/rdp/rdp

  call VRCAL(rdp, wt, sigp, factor)
  vr=mdot/(2.0*PI*rdp*Sigp)
  cs=dsqrt(Wt/Sigp)
  vr2cs=vr/cs/factor
  if(vr2cs.gt.1.0d0)then
    write(*,*)"vr2cs .gt.1.0", vr2cs
    stop
  end if
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
  Wtd(ird)=Wt
  Heigd(ird)=Heigp
  Wrd(ird)=Wr
  write(*,*)Heigd(ird), Wgd(ird), Wrd(ird), Sigp
END SUBROUTINE BOUNDARY
