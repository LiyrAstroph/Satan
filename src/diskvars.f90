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
                 Rd(nd), Teffd(nd), dWtd(nd), dSigd(nd), Qcord(nd), Qvisd(nd), Pbased(nd)
  integer(kind=8):: nt   

  integer(kind=8), parameter:: nnu_max=500
  real(kind=8):: nu(nnu_max), diskfnu(nnu_max)

  real(kind=8), parameter::eps=1.0d-10   
  
! corona
  real(kind=8), parameter::Ximin=0.5d0 

  real(kind=8) xleg_nu(nnu_max), wleg_nu(nnu_max)             
END MODULE DISKVARS
