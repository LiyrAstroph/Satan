!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

include "const.f90"
include "specvars.f90"

PROGRAM MAIN
  implicit none
  call specset
  call speccal
  call phindex
END PROGRAM MAIN

SUBROUTINE SPECSET
  use const
  use specvars
  implicit none
  real(kind=8) mbh, Ledd, rdin, rdout
  integer(kind=8) i

  nu_min=1.0d13
  nu_max=1.0d21

  omig_in_min=dlog10(h*nu_min/mec2)
  omig_in_max=dlog10(h*nu_max/mec2)
  
  gam_min=dlog10(1.001d0)
  gam_max=dlog10(50.0d0)

  call gauleg(omig_in_min,omig_in_max,xleg_omig,wleg_omig,nleg_omig)    
  call gauleg(gam_min,gam_max,xleg_gamm,wleg_gamm,nleg_gamm)
  call gauleg(-1.0d0,1.0d0,xleg,wleg,nleg)

  do i=1,nleg_omig
    nu(i)=10.0d0**(xleg_omig(i))*mec2/h
    nulnu(i)=0.0d0
    fnu_unred(i)=0.0d0
    flux_comp(i)=0.0d0
  end do

  open(unit=40, file="../data/data_for_spec.dat", status="old")
  read(40, *)mbh, rdin, epsilon
  do i=1, nnr 
    read(40,*, end=10)Rd(i), teffd(i), sigdd(i), heigdd(i), tecd(i), ticd(i),  &
        hcd(i), rhocd(i)
  enddo
10 nr=i-1
  
  Ledd=mbh * 4.0d0*PI*G*Msun*mp*c/sigmaT
  Medd=Ledd/0.1d0/c/c
  Rg=mbh * G*Msun/c/c
  rdout=1.0d3
  rdin=6.0d0
  if((Rd(1).gt.rdout).or.(Rd(nr).lt.rdin))then
    stop "rdin and rdout"
  endif

  darea(1)=pi*((rd(1)+rdout)**2/4.0-(rd(1)+rd(2))**2/4.0)*Rg**2      
  do i=2,nr-1
    darea(i)=pi*((rd(i-1)+rd(i))**2/4.0-(rd(i)+rd(i+1))**2/4.0)*Rg**2
  end do
  darea(nr)=pi*((rd(nr)+rd(nr-1))**2/4.0-(rdin+rd(nr))**2/4.0)*Rg**2

  close(40)
END SUBROUTINE SPECSET

SUBROUTINE SPECCAL
  use const
  use specvars
  implicit none
  real(kind=8) tep, teffp, sigp, heigp, tecp, rhocp, hcp
  real(kind=8) fluxr(nnu), fbr(nnu), fin, fcom
  integer(kind=8) inu, ird
  open(unit=50, file="../data/spectrum.dat")

  call rate_numb()

  do inu=1, nnu
    nu_unred(inu)=nu(inu)
  enddo

  do ird=1, nr
    write(*,*)ird, Rd(ird)
    teffp=teffd(ird)
    sigp=sigdd(ird)
    heigp=heigdd(ird)
    tep=teffp

    tecp=Tecd(ird)
    rhocp=rhocd(ird)
    hcp=hcd(ird)
    call diskrad(tep, teffp, sigp, heigp)
    call corbrem(tecp, rhocp, hcp, fbr)
    
    fnu_unred=fnu_unred+fbr
    do inu=1,nnu
      fluxr(inu)=fnu_unred(inu)
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*fnu_unred(inu)*darea(ird)*2.0d0
    enddo
    
    fin=0.0d0
    do inu=1, nnu-1
      fin=fin + fnu_unred(inu) * nu(inu) * (dlog(nu(inu+1))-dlog(nu(inu)))
    end do 
! first comptonization
    call compn(ird)
    do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)         &
                *2.0*(1.0d0-epsilon)
      write(24,*)nu(inu), fluxr(inu)*nu(inu), fnu_unred(inu)*nu(inu), flux_comp(inu)*nu(inu)
    enddo
    fcom=0.0d0
    do inu=1, nnu
      fcom=fcom + flux_comp(inu) * nu(inu) * (dlog(nu(inu+1))-dlog(nu(inu)))
    end do 
    write(*,*)fcom/fin

! second comptonization
    do inu=1,nnu
      fnu_unred(inu)=flux_comp(inu)
    enddo
    call compn(ird)
    do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)         &
                *2.0*(1.0d0-epsilon)
    enddo

! third comptonization
    do inu=1,nnu
      fnu_unred(inu)=flux_comp(inu)
    enddo
    call compn(ird)
    do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)         &
                *2.0*(1.0d0-epsilon)
    enddo

! forth comptonization
    do inu=1,nnu
      fnu_unred(inu)=flux_comp(inu)
    enddo
    call compn(ird)
    do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)         &
                *2.0*(1.0d0-epsilon)
    enddo
    write(24,*)
  enddo

  do inu=1,nnu
    write(50, *)nu(inu), nulnu(inu)
  end do

  close(50)
END SUBROUTINE SPECCAL

SUBROUTINE DISKRAD(tep, teffp, sigp, heigp)
  use const
  use specvars
  implicit none
  real(kind=8) tep, teffp, sigp, heigp
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
    fnu_unred(j)= Inu * PI
   end do
END SUBROUTINE DISKRAD

SUBROUTINE CORBREM(tecp, rhocp, hcp, fbr)
  use const
  use specvars
  implicit none
  real(kind=8) tecp, rhocp, hcp, fbr(nnu)
  real(kind=8) ftheta, qbr, ne, gaunt, chibr, thetae, temp_nu, temp_gaun
  integer(kind=8) inu
  
  ne=(rhocp*Medd/c/Rg/Rg)/mue/mp 
  thetae=kb*Tecp/mec2
  temp_gaun=h/kb/tecp

  if(thetae.lt.1.0d0)then
    ftheta=4.0*dsqrt(2.0*thetae/pi**3)*(1.0+1.781*thetae**(1.34d0))                 &
          +1.73*thetae**(1.5d0)*(1.0+1.1*thetae+thetae**2-1.25*thetae**(2.5d0))
  else
    ftheta=4.5d0*thetae/pi*(dlog(0.48+1.123*thetae)+1.5)                            &
     +2.30*thetae*(dlog(1.123*thetae)+1.28)
  end if
  qbr=1.48d-22*ne**2*ftheta

  do inu=1, nnu
    temp_nu=h*nu(inu)/kb/tecp
    if(temp_nu.lt.1.0d0)then
      gaunt=temp_gaun*dsqrt(3.0d0)/pi*dlog(4.0d0/1.781d0/temp_nu)
    else
      gaunt=temp_gaun*dsqrt(3.0/pi/temp_nu)
    endif
    chibr=qbr*gaunt*dexp(-temp_nu)
    fbr(inu)=chibr*hcp*Rg
  enddo
END SUBROUTINE CORBREM
