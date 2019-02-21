!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! SUBROUTINE INTEGRATE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integrate(z_w,q,z_iso,iqu,iql,L,M,N)

      implicit none
      integer L, M, N
      real*8 z_w(N+1,M,L)
      real*8 z_iso(M,L)
      real*8 q(N,M,L)
      real*8 iqu(M,L)
      real*8 iql(M,L)
!f2py intent(out) iqu
!f2py intent(out) iql
!f2py intent(hide) :: L
!f2py intent(hide) :: M
!f2py intent(hide) :: N
      integer i, j, k
      real*8 dz(N)
      real*8 dzp

      do i=1,L
        do j=1,M
          do k=1,N
            dz(k)=z_w(k+1,j,i)-z_w(k,j,i)
          enddo
          iqu(j,i)=0.0
          iql(j,i)=0.0
          do k=1,N
            if (z_w(k,j,i).gt.z_iso(j, i)) then
              iqu(j,i) = iqu(j,i) + dz(k)*q(k,j,i)
            else if (z_w(k+1,j,i).gt.(z_iso(j,i))) then
              dzp = z_w(k+1,j,i) - z_iso(j,i)
              iqu(j,i) = iqu(j,i) + dzp*q(k,j,i)
              dzp = z_iso(j,i) - z_w(k,j,i)
              iql(j,i) = iql(j,i) + dzp*q(k,j,i)
            else
              iql(j,i) = iql(j,i) + dz(k)*q(k,j,i)
            endif
          enddo
        enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! SUBROUTINE SURFACE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine surface(z, q, q0, z_iso, L, M, N)
! Assume q is sorted

      implicit none
      integer L, M, N
      real*8 z(N,M,L)
      real*8 q(N,M,L)
      real*8 q0(M,L)
      real*8 z_iso(M,L)
!f2py intent(out) z_iso
!f2py intent(hide) :: L
!f2py intent(hide) :: M
!f2py intent(hide) :: N
      integer i, j, k
      real*8 dz, dq, dq0

      do i=1,L
        do j=1,M
          z_iso(j,i)=1.0d20 ! default value - isoline not in profile
          do k=1,N-1
            if ( (q(k,j,i).lt.q0(j,i).and.q(k+1,j,i).gt.q0(j,i)).or.
     &           (q(k,j,i).gt.q0(j,i).and.q(k+1,j,i).lt.q0(j,i)) ) then
              dz = z(k+1,j,i) - z(k,j,i)
              dq = q(k+1,j,i) - q(k,j,i)
              dq0 = q0(j,i) - q(k,j,i)
              z_iso(j,i) = z(k,j,i) + dz*dq0/dq
            endif
          enddo
        enddo
      enddo

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! SUBROUTINE ZSLICE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zslice (z,f3d,depths,f2d,vinterp,L,M,N)

      implicit none
      integer L, M, N, vinterp
      real*8 depths(M,L)
      real*8 f3d(N,M,L), z(N,M,L)
      real*8 f2d(M,L)
!f2py intent(out) f2d
!f2py intent(hide) :: L
!f2py intent(hide) :: M
!f2py intent(hide) :: N
      integer i, j, k
      real*8 zk(N), fk(N)
      real*8 wk(N), dfdz(N)

!  Linear Interpolation.
      if (vinterp.eq.0) then
        do i=1,L
          do j=1,M
            do k=1,N
              fk(k)=f3d(k,j,i)
              zk(k)=z(k,j,i)
            enddo
            if ((zk(1).le.depths(j,i)).and.(depths(j,i).le.zk(N))) then
              call lintrp (N,zk,fk,1,depths(j,i),f2d(j,i))
            else
              f2d(j,i)=1.0d20
            endif
          enddo
        enddo
!  Cubic spline interpolation.
      else if (vinterp.eq.1) then
        do i=1,L
          do j=1,M
            do k=1,N
              fk(k)=f3d(k,j,i)
              zk(k)=z(k,j,i)
            enddo
            call spline (zk,fk,N,1.0d30,1.0d30,wk)
            if ((zk(1).le.depths(j,i)).and.(depths(j,i).le.zk(N))) then
              call splint (zk,fk,wk,N,depths(j,i),f2d(j,i),dfdz)
            else
              f2d(j,i)=1.0d20
            endif
          enddo
        enddo
      endif

      return
      end




      subroutine lintrp (n,x,y,ni,xi,yi)
!
!=======================================================================
!  Copyright (c) 1996 Rutgers University                             ===
!=======================================================================
!                                                                    ===
!  Given arrays X and Y of length N, which tabulate a function,      ===
!  Y = F(X),  with the Xs  in ascending order, and given array       ===
!  XI of lenght NI, this routine returns a linear interpolated       ===
!  array YI.                                                         ===
!                                                                    ===
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define local variable.
!-----------------------------------------------------------------------
!
      implicit none
      integer i, ii, j, n, ni
      real*8 d1, d2
      real*8 x(n), y(n), xi(ni), yi(ni)
!f2py intent(out) yi
!f2py intent(hide) n
!f2py intent(hide) ni
!
!-----------------------------------------------------------------------
!  Begin executable code.
!-----------------------------------------------------------------------
!
      do 30 j=1,ni
        if (xi(j).le.x(1)) then
          ii=1
          yi(j)=y(1)
        elseif (xi(j).ge.x(n))then
          yi(j)=y(n)
        else
          do 10 i=1,n-1
            if ((x(i).lt.xi(j)).and.(xi(j).le.x(i+1))) then
              ii=i
              goto 20
            endif
 10       continue
 20       d1=xi(j)-x(ii)
          d2=x(ii+1)-xi(j)
          yi(j)=(d1*y(ii+1)+d2*y(ii))/(d1+d2)
        endif
 30   continue
      return
      end


      subroutine spline (x,y,n,yp1,ypn,y2)
!
!=======================================================================
!  Copyright (c) 1996 Rutgers University                             ===
!=======================================================================
!                                                                    ===
!  Given X, Y of length N containing a tabulated function,  Y=f(X),  ===
!  with the Xs  in ascending order,  and given values  Yp1 and  Ypn  ===
!  for the first derivative of the interpolating function at points  ===
!  1 and N, respectively this routine returns an array Y2 of length  ===
!  N  which contains the  second  derivatives of the  interpolating  ===
!  function at the tabulated points X.  If Yp1 and/or Ypn are equal  ===
!  to  1.0E+30  or larger,  the routine  is  signalled  to  set the  ===
!  corresponding boundary condition for a natural spline, with zero  ===
!  second derivative on that boundary.                               ===
!                                                                    ===
!  Reference :                                                       ===
!                                                                    ===
!  Press, W.H, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,   ===
!        1986: Numerical Recipes, the art of scientific computing.   ===
!        Cambridge University Press.                                 ===
!                                                                    ===
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define global data.
!-----------------------------------------------------------------------
!
      real*8 cm1,cm3,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c20,c25,c50,
     &     c90,c100,c180,c200,c255,c300,c360,c366,c500,c1000,c5000,
     &     c10000,c1em9,c1em10,c1em12,c1em20,c1ep30,p006,p009,p035,
     &     p015,p012,p08,p06,p5,p25,p75,p98,r3,r10,r20,r33,r35,r40,
     &     r50,r80,r100,r200,r250,r400,r1000
      real*8 day2sec,deg2rad,grav,cm2m,m2cm,m2km,pi,rad2deg,re,root2,
     &     sec2day,spval0,spval1,spval2,spvgeo
!
      parameter (cm1=-1.0,cm3=-3.0,c0=0.0,c1=1.0,c2=2.0,c3=3.0,c4=4.0,
     &           c5=5.0,c6=6.0,c7=7.0,c8=8.0,c9=9.0,c10=10,c11=11.0,
     &           c20=20.0,c25=25.0,c50=50.0,c90=90.0,c100=100.0,
     &           c180=180.0,c200=200.0,c255=255.0,c300=300.0,c360=360.0,
     &           c366=366.0,c500=500,c1000=1000.0,c5000=5000.0,
     &           c10000=10000.0,c1em9=1.0e-9,c1em10=1.0e-10,
     &           c1em12=1.0e-12,c1em20=1.0e-20,c1ep30=1.0e+30,
     &           p006=0.006,p009=0.009,p012=0.012,p015=0.015,p035=0.035,
     &           p06=0.06,p08=0.08,p5=0.5,p25=0.25,p75=0.75,p98=0.98,
     &           r3=c1/c3,r10=0.1,r20=0.05,r33=c1/33.0,r35=c1/35.0,
     &           r40=0.025,r50=0.02,r80=0.0125,r100=0.01,r200=0.005,
     &           r250=0.004,r400=0.0025,r1000=0.001)
      parameter (day2sec=86400.0,cm2m=r100,grav=9.8,m2cm=c100,
     &           m2km=r1000,pi=3.14159265358979323846,re=637131500.0,
     &           root2=1.41421356237309504880,sec2day=c1/86400.0,
     &           spvgeo=999.0,spval0=0.99e+35,spval1=1.0e+35,
     &           spval2=0.99e+30)
      parameter (deg2rad=pi/c180,rad2deg=c180/pi)

!
!-----------------------------------------------------------------------
!  Define local data.  Change NMAX as desired to be the largest
!  anticipated value of N.
!-----------------------------------------------------------------------
!
      integer i, k, n, nmax
      parameter (nmax=10000)
      real*8 p, qn, sig, un, ypn, yp1
      real*8 x(n), y(n), y2(n), u(nmax)
!f2py intent(out) :: y2
!f2py intent(hide) :: n
!
!-----------------------------------------------------------------------
!  Begin excutable code.
!-----------------------------------------------------------------------
!
      if (n.gt.nmax) then
        print 10, n,nmax
 10     format(/' SPLINE: underdimensioned array, N, NMAX = ',2i5)
        call crash ('SPLINE',1)
      endif
!
!  The lower boundary condition is set either to be "natural" or else
!  to have a specified first derivative.
!
      if (yp1.gt.spval2) then
        y2(1)=c0
        u(1)=c0
      else
        y2(1)=-p5
        u(1)=(c3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
!
!  This is the decomposition loop of the tridiagonal algorithm. Y2 and
!  U are used for temporary storage of the decomposition factors.
!
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+c2
        y2(i)=(sig-c1)/p
        u(i)=(c6*((y(i+1)-y(i))/(x(i+1)-x(i))-
     &           (y(i)-y(i-1))/(x(i)-x(i-1)))/
     &           (x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
!
!  The upper boundary condition is set either to be "natural" or else
!  to have a specified first derivative.
!
      if (ypn.gt.spval2) then
        qn=c0
        un=c0
      else
        qn=p5
        un=(c3/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+c1)
!
!  This is the back-substitution loop of the tridiagonal algorithm.
!
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      end

      subroutine splint (x,y,y2,n,xx,yy,dydx)
!
!=======================================================================
!  Copyright (c) 1996 Rutgers University                             ===
!=======================================================================
!                                                                    ===
!  Given arrays X and Y of length N, which tabulate a function,      ===
!  Y=f(X), with the Xs  in ascending order, and given the array      ===
!  Y2 which contains the second derivative of the interpolating      ===
!  function  at the  tabulated points X as computed  by routine      ===
!  SPLINE, and given a value  XX, this routine returns a cubic-      ===
!  spline interpolated value YY.                                     ===
!                                                                    ===
!  Reference :                                                       ===
!                                                                    ===
!  Press, W.H, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,   ===
!         1986: Numerical Recipes, the art of scientific computing.  ===
!         Cambridge University Press.                                ===
!                                                                    ===
!  Modified by H.G. Arango (1989) to output the first derivative     ===
!  DYDX at a given value XX.                                         ===
!                                                                    ===
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define global data.
!-----------------------------------------------------------------------
!
      implicit none
      real*8 cm1,cm3,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c20,c25,c50,
     &     c90,c100,c180,c200,c255,c300,c360,c366,c500,c1000,c5000,
     &     c10000,c1em9,c1em10,c1em12,c1em20,c1ep30,p006,p009,p035,
     &     p015,p012,p08,p06,p5,p25,p75,p98,r3,r10,r20,r33,r35,r40,
     &     r50,r80,r100,r200,r250,r400,r1000
      real*8 day2sec,deg2rad,grav,cm2m,m2cm,m2km,pi,rad2deg,re,root2,
     &     sec2day,spval0,spval1,spval2,spvgeo
!
      parameter (cm1=-1.0,cm3=-3.0,c0=0.0,c1=1.0,c2=2.0,c3=3.0,c4=4.0,
     &           c5=5.0,c6=6.0,c7=7.0,c8=8.0,c9=9.0,c10=10,c11=11.0,
     &           c20=20.0,c25=25.0,c50=50.0,c90=90.0,c100=100.0,
     &           c180=180.0,c200=200.0,c255=255.0,c300=300.0,c360=360.0,
     &           c366=366.0,c500=500,c1000=1000.0,c5000=5000.0,
     &           c10000=10000.0,c1em9=1.0e-9,c1em10=1.0e-10,
     &           c1em12=1.0e-12,c1em20=1.0e-20,c1ep30=1.0e+30,
     &           p006=0.006,p009=0.009,p012=0.012,p015=0.015,p035=0.035,
     &           p06=0.06,p08=0.08,p5=0.5,p25=0.25,p75=0.75,p98=0.98,
     &           r3=c1/c3,r10=0.1,r20=0.05,r33=c1/33.0,r35=c1/35.0,
     &           r40=0.025,r50=0.02,r80=0.0125,r100=0.01,r200=0.005,
     &           r250=0.004,r400=0.0025,r1000=0.001)
      parameter (day2sec=86400.0,cm2m=r100,grav=9.8,m2cm=c100,
     &           m2km=r1000,pi=3.14159265358979323846,re=637131500.0,
     &           root2=1.41421356237309504880,sec2day=c1/86400.0,
     &           spvgeo=999.0,spval0=0.99e+35,spval1=1.0e+35,
     &           spval2=0.99e+30)
      parameter (deg2rad=pi/c180,rad2deg=c180/pi)

!
!-----------------------------------------------------------------------
!  Define local data.
!-----------------------------------------------------------------------
!
      integer k, khi, klo, n
      real*8 a, b, c, d, dydx, e, f, h, xx, yy
      real*8 x(n), y(n), y2(n)
!f2py intent(out) :: yy
!f2py intent(hide) :: n
!
!-----------------------------------------------------------------------
!  Begin executable code.
!-----------------------------------------------------------------------
!
!  Found the right place of XX in the table by means of bisection.
!
      klo=1
      khi=n
  10  if ((khi-klo).gt.1) then
        k=(khi+klo)/2
        if(x(k).gt.xx) then
          khi=k
        else
          klo=k
        endif
        goto 10
      endif
!
!  KLO and KHI now bracket the input value XX.
!
      h=x(khi)-x(klo)
      if (h.eq.c0) then
        print *, ' SPLINT: bad X input, they must be distinct.'
        call crash ('SPLINT',1)
      endif
!
!  Evaluate cubic spline polynomial.
!
      a=(x(khi)-xx)/h
      b=(xx-x(klo))/h
      c=(a*a*a-a)*(h*h)/c6
      d=(b*b*b-b)*(h*h)/c6
      e=(c3*(a*a)-c1)*h/c6
      f=(c3*(b*b)-c1)*h/c6
      yy=a*y(klo)+b*y(khi)+c*y2(klo)+d*y2(khi)
      dydx=(y(khi)-y(klo))/h-e*y2(klo)+f*y2(khi)
      return
      end

      subroutine crash (string,ierr)
      integer ierr
      character*(*) string

      if (string(1:4).eq.'DONE') then
        print 10, string
  10    format(a)
        ierr=0
      else
        print 20, string
  20    format(/,' Execution abnormally terminated in module: ',a,/)
        ierr=1
      endif
      return
      end




