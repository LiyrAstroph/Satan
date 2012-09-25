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
  implicit none
  call init
  call boundary
  call solve
END PROGRAM MAIN

