!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.10.07
!==============================================================================!


SUBROUTINE CORSOLVE
  use const
  use diskvars
  implicit none
  real(kind=8), external::corcomp, Qie
  real(kind=8) rdp, tep, teffp, sigp, heigp, tecp, ticp, hcp, rhocp,        &
               OmgK, Pcp, ne, ni, y, qcorp
  real(kind=8) fd(nnu_max)
  integer(kind=8) ird, np, np2, nwt, ndw
  
  open(unit=20, file="../data/corona_eng.dat")
  open(unit=22, file="../data/corona.dat")
  open(unit=40, file="../data/data_for_spec.dat")
  write(40,*)mbh, rin, epsilon
   
  np=0
  do ird=1, nt
    rdp=Rd(ird)
    if(rdp.lt.1.0d3)then
      np=ird
      exit
    endif
  end do

  do ird=1, nt
    rdp=Rd(ird)
    if(rdp.le.6.0d0)then
      np2=ird-1
      exit
    endif
  end do
  nwt=np2-np+1
  ndw=nwt/120
  write(*,*)Rd(np), np2, nwt, ndw

  tecp=1.0d9
  ticp=3.0d9
  do ird=np, np2, ndw
    rdp=Rd(ird)
    tep=Td(ird)
    teffp=Teffd(ird)
    sigp=Sigd(ird)
    heigp=Heigd(ird)
    OmgK=1.0d0/dsqrt(rdp-2.0d0)/rdp
    qcorp=Qcord(ird)*Medd*c*c/Rg/Rg
    Pcp=Pbased(ird)
    call diskrad(tep, teffp, sigp, heigp, fd)
    call corona(rdp, pcp, qcorp, fd, tecp, ticp) 
!    tecp=1.0d10/rdp
!    ticp=3.0d11/rdp
    Pcp=pcp
    if(tecp.gt.1.0d6)then
    rhocp=Pcp*mp/(kb*(tecp/mue + ticp/mui)) *c*c
    hcp=dsqrt(Pcp/rhocp)/OmgK
    else
    rhocp=1.0d-20
    hcp=Rg
    endif
    ne=(rhocp*Medd/c/Rg/Rg)/mp/mue
    ni=ne*mue/mui
    y=ne*sigmat*hcp*Rg * 4.0d0*kb*tecp/mec2
!               1     2     3    4     5     6       7                 8       
    write(22,*)rdp, tecp, ticp, hcp, rhocp, pcp, 0.5*sigp/heigp, 0.5*Wtd(ird)/heigp
    write(20,*)rdp, rhocp, hcp, y, corcomp(tecp, rhocp, hcp, fd),                  &
    Qie(tecp, ticp, ne, ni)*2.0d0*hcp*Rg, 2.0d0*sigmab*teffp**4.0d0,               &
    Qcord(ird)*Medd*c*c/Rg/Rg, ne*sigmat*hcp*Rg
    
    write(40, *)rdp, teffp, sigp, heigp, tecp, ticp, hcp, rhocp
  enddo
  close(20)
  close(22)
  close(40)
END SUBROUTINE CORSOLVE


FUNCTION CORCOMP(tecp, rhocp, hcp, fd)
  use const
  use diskvars
  implicit none
  real(kind=8) tecp, rhocp, hcp, corcomp, fd(nnu_max)
  real(kind=8) etad(nnu_max), fcool(nnu_max), etat
  integer(kind=8) inu

  call etaComp(tecp, rhocp, hcp, etad)
  
  fcool=fd*(etad-1.0d0)

  corcomp=0.0d0
  etat=0.0d0
  do inu=1, nnu
    corcomp=corcomp+ fcool(inu) * nu(inu) * wleg_nu(inu)
    etat=etat + etad(inu)
  end do
  corcomp=corcomp * 2.0d0 * dlog(10.0d0)
  return
END FUNCTION CORCOMP

!***********************************************************************!
! cal the enhancement factor for Comptonization
!   
!***********************************************************************!
subroutine etaComp(tecp, rhocp, hcp, etad)
  use const
  use diskvars
  implicit none
  real(kind=8) etad(nnu_max), tecp, rhocp, hcp
  real(kind=8), external::gammp
  real(kind=8) A, thetae, taues, s, eta, etamax, jm, temp_nu, ne,    &
               rhop, heigp
  integer(kind=8) inu

  rhop=rhocp*(Medd/c/Rg/Rg)
  heigp=hcp*Rg

  ne=rhop/mp/mue 
  thetae=kb*tecp/mec2
  A=1.0+4.0d0*thetae+16.0d0*thetae**2.0d0
  taues=ne*sigmat*heigp
  s=taues+taues**2.0d0
  do inu=1, nnu
  temp_nu=h*nu(inu)/kb/tecp  
  etamax=3.0d0/temp_nu
  if(etamax.ge.1.0d0)then
  jm=dlog(etamax)/dlog(A)
  jm=min(jm, 100.0d0)
  if(s*A.gt.1.0d-6)then
    eta=dexp(s*(A-1.0d0))*(1.0d0-gammp(jm+1.0d0,A*s))+etamax*gammp(jm+1.0d0,s)
  else
    eta=1.0d0 + s*(A-1.0d0)
  endif
  eta=max(1.0d0,eta)
  else
  jm=0.0d0
  eta=1.0d0      
  end if
  etad(inu)=eta
  enddo
  return               
end subroutine etaComp

!***********************************************************************!
! function radheat
! note the unit in cgs
!***********************************************************************!      
FUNCTION Qie(te,ti,ne,ni)
  use const
  implicit none
  real(kind=8) ti,te,ne,ni,Qie
  real(kind=8) bessk,bessk0,bessk1
  real(kind=8) thetai,thetae,temp
  
  Qie=0.0d0  
  thetai=kb*ti/mpc2
  thetae=kb*te/mec2
  if((thetae.ge.1.0d-2).and.(thetai.ge.1.0d-2))then
  Qie=5.61d-32*ne*ni*(ti-te)/(bessk(2,1.0/thetae)*bessk(2,1.0/thetai))         &
          *((2.0*(thetae+thetai)**2+1.0)/(thetae+thetai)                       &
          *bessk1((thetae+thetai)/thetae/thetai)                               &
          +2.0*bessk0((thetae+thetai)/thetae/thetai))                          
  return
  end if
      
  temp=thetae*thetai/(thetae+thetai) 
           
  if((thetae.ge.1.0d-2).and.(thetai.lt.1.0d-2))then
  Qie =5.61d-32*ne*ni*(ti-te)/bessk(2,1.0/thetae)*dexp(-1.0/thetae)            &
      *dsqrt(thetae/(thetae+thetai))                                           &
      *((2.0*(thetae+thetai)**2+1)/(thetae+thetai)                             &
      *(1.0+3.0/8.0*temp-15.0/128.0*temp**2+15.0*21.0/6.0/8.0**3*temp**3)      &
      +2.0*(1.0-1.0/8.0*temp+9.0/128.0*temp**2-9.0*25.0/6.0/8.0**3*temp**3))   &
      /(1.0+15.0/8.0*thetai+15.0*7.0/128.0*thetai**2                           &
      -15.0*7.0*9.0/6.0/8.0**3*thetai**3)          
  return
  end if 
      
  if((thetae.lt.1.0d-2).and.(thetai.ge.1.0d-2))then
  Qie=5.61d-32*ne*ni*(ti-te)/bessk(2,1.0/thetai)*dexp(-1.0/thetai)             &
      *dsqrt(thetai/(thetae+thetai))                                           &
      *((2.0*(thetae+thetai)**2+1)/(thetae+thetai)                             &
      *(1.0+3.0/8.0*temp-15.0/128.0*temp**2+15.0*21.0/6.0/8.0**3*temp**3)      &
      +2.0*(1.0-1.0/8.0*temp+9.0/128.0*temp**2-9.0*25.0/6.0/8.0**3*temp**3))   &
      /(1.0+15.0/8.0*thetae+15.0*7.0/128.0*thetae**2                           &
      -15.0*7.0*9.0/6.0/8.0**3*thetae**3)           
  return
  end if          
      
  if((thetae.lt.1.0d-2).and.(thetai.lt.1.0d-2))then  
  Qie=5.61d-32*ne*ni*(ti-te)*dsqrt(2.0/pi)/dsqrt(thetae+thetai)                &
      *((2.0*(thetae+thetai)**2+1)/(thetae+thetai)                             &
      *(1.0+3.0/8.0*temp-15.0/128.0*temp**2+15.0*21.0/6.0/8.0**3*temp**3)      &
      +2.0*(1.0-1.0/8.0*temp+9.0/128.0*temp**2                                 &
           -9.0*25.0/6.0/8.0**3*temp**3))                                      &
      /(1.0+15.0/8.0*thetai+15.0*7.0/128.0*thetai**2                           &
          -15.0*7.0*9.0/6.0/8.0**3*thetai**3)                                  &
      /(1.0+15.0/8.0*thetae+15.0*7.0/128.0*thetae**2                           &
          -15.0*7.0*9.0/6.0/8.0**3*thetae**3)             
  return
  end if    
      
end function Qie
