!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

MODULE CONST
  implicit none
  real(kind=8), parameter:: G=6.672d-8, c=2.9979d10, mp=1.672d-24, me=9.1095d-28,  &
          kb=1.38066d-16,                                                          &
          sigmab=5.6703d-5, ab=7.56566d-15, sigmaT=6.6524d-25, Msun=1.9891d33,     &
          h=6.6262d-27, PI=2.0d0*dacos(0.0d0)
  real(kind=8), parameter:: mpc2=mp*c**2,mec2=me*c**2    
  real(kind=8), parameter:: mui=1.23, mue=1.14
END MODULE CONST
