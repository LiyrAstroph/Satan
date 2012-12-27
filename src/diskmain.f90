!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

include "const.f90"
include "diskvars.f90"

PROGRAM MAIN
  use diskvars
  implicit none
  call init
  call boundary
  call solve
  call diskspec
  call corr_factor
  if(betam.gt.0.0d0)call corsolve
END PROGRAM MAIN

