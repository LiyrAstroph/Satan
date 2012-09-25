!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

MODULE DISKVARS
  implicit none
  real(kind=8):: Ledd, Medd, Rg
  
  real(kind=8):: alpha, mu, mbh, mdot, gam, betam, epsilon, ar, ellin_min,         &
      ellin_max
  namelist /slim/ alpha, mu, mbh, mdot, gam, betam, epsilon, ar, ellin_min,        &
      ellin_max
  real(kind=8):: rin, rout, drstep_min, drstep_max, numin, numax
  integer(kind=8):: nnu
  namelist /mesh/ rin, rout, drstep_min, drstep_max, nnu, numin, numax
  
! slim disk
  real(kind=8)::ellin, muave 
  integer(kind=8), parameter::nd=2000000
  real(kind=8):: Wtd(nd), Wgd(nd), Wrd(nd), Sigd(nd), Td(nd), Omgd(nd),Heigd(nd),  &
                 Rd(nd), Teff(nd), dfd(nd), dWtd(nd), dSigd(nd)
  integer(kind=8):: nt             
  real(kind=8), parameter::eps=1.0d-20   

! corona
  real(kind=8), parameter::Ximin=0.05d0              
END MODULE DISKVARS
