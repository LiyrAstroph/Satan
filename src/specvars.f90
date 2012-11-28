!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

MODULE SPECVARS
  integer(kind=8), parameter::nnu=200, nnr=200
  integer(kind=8) nr
  real(kind=8) nu(nnu), nu_unred(nnu), fnu_unred(nnu), flux_comp(nnu), nulnu(nnu)
  real(kind=8) Rd(nnr), teffd(nnr), sigdd(nnr), heigdd(nnr), tecd(nnr), ticd(nnr),    &
               rhocd(nnr), hcd(nnr), darea(nnr)

  integer(kind=8), parameter::nleg=100
  real(kind=8) xleg(nleg), wleg(nleg)

  integer(kind=8), parameter::nleg_gamm=200, nleg_omig=nnu

  real(kind=8):: rate(nleg_omig,nleg_gamm),   numb_elec(nleg_gamm),                    &
      domi(nleg_omig,nleg_gamm),omigm(nleg_omig,nleg_gamm)

  real(kind=8):: xleg_gamm(nleg_gamm),wleg_gamm(nleg_gamm)
  real(kind=8):: xleg_omig(nleg_omig),wleg_omig(nleg_omig)

  real(kind=8) omig_in_min,omig_in_max,nu_min,nu_max
  real(kind=8) gam_min,gam_max

  real(kind=8) Medd, Rg, epsilon
END MODULE SPECVARS