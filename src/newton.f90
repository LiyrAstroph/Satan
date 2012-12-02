!==============================================================================!
!*************************** Satan   ******************************************!
!
!                         By Yan-Rong LI
!                       liyropt@gmail.com
!                     2012.09.19-2012.09.23
!==============================================================================!

module corvars
  use diskvars
  real(kind=8) r, pc, qcor, fdnu(nnu_max)
end module corvars

subroutine corona(rdp, pcp, qcorp, fd, tecp, ticp)
  use const
  use diskvars
  use corvars
  implicit none
  real(kind=8) rdp, pcp, qcorp, fd(nnu_max), tecp, ticp
  integer(kind=8) ntrial, n
  real(kind=8) x(2), tolx, tolf
  
  r=rdp
  pc=pcp
  qcor=qcorp
  fdnu=fd

  n=2
  ntrial=100000
  x(1)=tecp
  x(2)=ticp

  tolx=0.01*x(1)
  tolf=1.0d-10
  call mnewt(ntrial, x, n, tolx, tolf)
  tecp=x(1)
  ticp=x(2)
  write(*,'(3e15.7)')rdp, x(1), x(2)
end subroutine corona

subroutine usrfun(x, n, NP, fvec, fjac)
  implicit none
  integer(kind=8) n, NP
  real(kind=8) x(n), fvec(NP), fjac(NP, NP)

  call funcv(n, x, fvec, NP)
  call fdjac(n, x, fvec, np, fjac)
!  write(30,*)fjac(1, 1), fjac(1, 2), fjac(2,1), fjac(2,2)
!  write(32,*)x(1), x(2), fvec(1), fvec(2)
end 

subroutine funcv(n, x, f, NP)
  use const
  use diskvars
  use corvars
  implicit none
  integer(kind=8) n, NP
  real(kind=8) x(n), f(NP)
  real(kind=8), external::Qie, corcomp
  real(kind=8) tecp, ticp, OmgK, rhocp, hcp, ne, ni, Qiep,        &
           Qcp, pcp

  tecp=x(1)
  ticp=x(2)
  pcp=Pc
  rhocp=Pcp*mp/(kb*(tecp/mue + ticp/mui)) * c*c  
  OmgK=1.0d0/dsqrt(r-2.0d0)/r
  hcp=dsqrt(Pcp/rhocp)/OmgK

  ne=(rhocp*Medd/c/Rg/Rg)/mp/mue
  ni=ne*mue/mui
  Qiep=Qie(tecp, ticp, ne, ni)
  Qiep=Qiep*2.0d0*hcp*Rg
  Qcp=corcomp(tecp, rhocp, hcp, fdnu)

!  write(30,*)"func", Qiep, qcor, qcp, ne
  f(1)=Qiep-(0.5*qcor)
  f(2)=Qcp-(Qiep+0.5*qcor)

!  write(*,*)r, tecp, ticp, hcp, qiep, qcp, qcor
  return
end 

SUBROUTINE mnewt(ntrial,x,n,tolx,tolf)
      INTEGER(kind=8) n,ntrial,NP
      REAL(kind=8) tolf,tolx,x(n)
      PARAMETER (NP=15)
!U    USES lubksb,ludcmp,usrfun
      INTEGER(kind=8) i,k,indx(NP)
      REAL(kind=8) d,errf,errx,fjac(NP,NP),fvec(NP),p(NP)
      do 14  k=1,ntrial
        call usrfun(x,n,NP,fvec,fjac)
        errf=0.
        do 11 i=1,n
          errf=errf+dabs(fvec(i))
11      continue
        if(errf.le.tolf)return
        do 12 i=1,n
          p(i)=-fvec(i)
12      continue
        call ludcmp(fjac,n,NP,indx,d)
        call lubksb(fjac,n,NP,indx,p)
        errx=0.
        do 13 i=1,n
          errx=errx+dabs(p(i))
          x(i)=x(i)+1.0d-3*p(i)
          x(i)=max(x(i), 1.0d6)
!          write(*,'(I6, e15.7)')k, x(i)
13      continue
        if(errx.le.tolx)return
14    continue
      write(*,*)"trial times", ntrial
      return
END

SUBROUTINE fdjac(n,x,fvec,np,df)
      INTEGER(kind=8) n,np,NMAX
      REAL(kind=8) df(np,np),fvec(n),x(n),EPS
      PARAMETER (NMAX=40,EPS=1.d-4)
!U    USES funcv
      INTEGER(kind=8) i,j
      REAL(kind=8) h,temp,f(n)
      do 12 j=1,n
        temp=x(j)
        h=EPS*dabs(temp)
        if(h.eq.0.)h=EPS
        x(j)=temp+h
        h=x(j)-temp
        call funcv(n,x,f, NMAX)
        x(j)=temp
        do 11 i=1,n
          df(i,j)=(f(i)-fvec(i))/h
11      continue
12    continue
      return
END

SUBROUTINE ludcmp(a,n,np,indx,d)
  INTEGER(kind=8) n,np,indx(n),NMAX
  REAL(kind=8) d,a(np,np),TINY
  PARAMETER (NMAX=500,TINY=1.0d-20)
  INTEGER(kind=8) i,imax,j,k
  REAL(kind=8) aamax,dum,sum,vv(NMAX)
  imax=0
  d=1.
  do 12 i=1,n
    aamax=0.
    do 11 j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11  continue
    if (aamax.eq.0.) then
    write(*,*)'singular matrix in ludcmp'
    stop
    endif
    vv(i)=1./aamax
12  continue
    do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19  continue
  return
END

SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER(kind=8) n,np,indx(n)
      REAL(kind=8) a(np,np),b(n)
      INTEGER(kind=8) i,ii,j,ll
      REAL(kind=8) sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
END
