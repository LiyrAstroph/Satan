!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

SUBROUTINE PHINDEX
  use const
  use specvars
  implicit none
  real(kind=8) nu1, nu2, nufit(nnu), nulnufit(nnu)
  real(kind=8) sig, mwt, a, b, siga, sigb, chi2, q
  integer(kind=8) inu, ip
  
  nu1=2.0d3*ev2Hz
  nu2=10.0d3*ev2Hz
  
  ip=0
  do inu=1, nnu
    if(nu(inu).ge.nu1.and.nu(inu).le.nu2)then
      ip=ip+1
      nufit(ip)=log10(nu(inu))
      nulnufit(ip)=log10(nulnu(inu))
    endif
  enddo
  
  call fit(nufit, nulnufit, ip, sig,mwt,a,b,siga,sigb,chi2,q)
  write(*,*)a, -b+2
END SUBROUTINE PHINDEX
