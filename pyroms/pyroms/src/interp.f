!=======================================================================
!=======================================================================
!
      SUBROUTINE xhslice (f2d,f3d,z,depth,mask,im,jm,km,vintrp,spval)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine extract a horizontal slice from a 3D field at the      !
!  requested depth via interpolation.                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     depth    Depth to interpolate.                                   !
!     f3d      3D field.                                               !
!     z        Depths of 3D field.                                     !
!     im       Size of inner dimension.                                !
!     jm       Size of outter dimension.                               !
!     km       Number of vertical levels.                              !
!     vintrp   Interpolation scheme:                                   !
!                vintrp = 0  linear                                    !
!                vintrp = 1  cubic splines                             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     f2d      Interpolated field.                                     !
!                                                                      !
!  Calls:                                                              !
!                                                                      !
!     linterp                                                          !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define global variables.
!-----------------------------------------------------------------------
!
      implicit none
      integer NH, NK, NMSK, NV, NX
      parameter (NH=6000000,NK=100,NMSK=100000,NV=7000000,NX=10000)
!     parameter (NH=1400000,NK=100,NMSK=100000,NV=31000000,NX=10000)
      integer cntplt, secplt
      parameter (cntplt=1, secplt=2)
!
      real*8 cm1, cm3, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,&
     &     c20, c25, c50, c90, c100, c180, c200, c255, c300, c360, c366,&
     &     c500, c1000, c5000, c10000, c1em9, c1em10, c1em12, c1em20,   &
     &     c1ep30, p006, p009, p035, p015, p012, p08, p06, p5, p25, p75,&
     &     p98, r3, r10, r20, r33, r35, r40, r50, r80, r100, r200, r250,&
     &     r400, r1000
      real*8 day2sec, deg2rad, grav, cm2m, cm2mm, m2cm, m2km, m2knot,   &
     &     m2mm, m2smallknot, pi, rad2deg, re, root2, sec2day, spval0,  &
     &     spval1, spval2, spvgeo
!
      parameter (cm1=-1.0, cm3=-3.0, c0=0.0, c1=1.0, c2=2.0, c3=3.0,    &
     &           c4=4.0, c5=5.0, c6=6.0, c7=7.0, c8=8.0, c9=9.0, c10=10,&
     &           c11=11.0, c20=20.0, c25=25.0, c50=50.0, c90=90.0,      &
     &           c100=100.0, c180=180.0, c200=200.0, c255=255.0,        &
     &           c300=300.0, c360=360.0, c366=366.0, c500=500,          &
     &           c1000=1000.0, c5000=5000.0, c10000=10000.0,            &
     &           c1em9=1.0e-9, c1em10=1.0e-10, c1em12=1.0e-12,          &
     &           c1em20=1.0e-20, c1ep30=1.0e+30, p006=0.006, p009=0.009,&
     &           p012=0.012, p015=0.015, p035=0.035, p06=0.06, p08=0.08,&
     &           p5=0.5, p25=0.25, p75=0.75, p98=0.98, r3=c1/c3,        &
     &           r10=0.1, r20=0.05, r33=c1/33.0, r35=c1/35.0, r40=0.025,&
     &           r50=0.02, r80=0.0125, r100=0.01, r200=0.005,           &
     &           r250=0.004, r400=0.0025, r1000=0.001)
      parameter (day2sec=86400.0, cm2m=r100, cm2mm=r10, grav=9.801,     &
     &           m2cm=c100, m2km=r1000, m2knot=1.9459, m2mm=c1000,      &
     &           m2smallknot=m2knot*c5, pi=3.14159265358979323846,      &
     &           re=637131500.0, root2=1.41421356237309504880,          &
     &           sec2day=c1/86400.0, spvgeo=999.0, spval0=0.99e+35,     &
     &           spval1=1.0e+35, spval2=0.99e+30)
      parameter (deg2rad=pi/c180, rad2deg=c180/pi)
!
!-----------------------------------------------------------------------
!  Define local variables.
!-----------------------------------------------------------------------
!
      integer i, im, j, jm, k, km, vintrp
      real*8 der1, derkm, frstd, spval
      real*8 f3d(km,jm,im), z(km,jm,im)
      real*8 depth(jm,im), mask(jm,im)
      real*8 f2d(jm,im)
!f2py intent(out) f2d
!f2py intent(hide) im
!f2py intent(hide) jm
      real*8 fk(NK), zk(NK), wk(NK)
      parameter (der1=c1ep30,derkm=c1ep30)
!
!=======================================================================
!  Begin executable code.
!=======================================================================
!
!  Linear Interpolation.
!
      IF (vintrp.eq.0) THEN
        DO j=1,jm
          DO i=1,im
            DO k=1,km
              fk(k)=f3d(k,j,i)
              zk(k)=z(k,j,i)
            END DO
            IF (((zk(1).le.depth(j,i)).and.(depth(j,i).le.zk(km)))      &
     &          .and. mask(j,i).eq.1) THEN
              CALL lintrp (km,zk,fk,1,depth(j,i),f2d(j,i))
            ELSE
              f2d(j,i)=spval
            END IF
          END DO
        END DO
!
!  Cubic spline interpolation.
!
      ELSE IF (vintrp.eq.1) THEN
        DO j=1,jm
          DO i=1,im
            DO k=1,km
              fk(k)=f3d(k,j,i)
              zk(k)=z(k,j,i)
            END DO
            CALL spline (zk,fk,km,der1,derkm,wk)
            IF (((zk(1).le.depth(j,i)).and.(depth(j,i).le.zk(km)))      &
     &          .and. mask(j,i).eq.1) THEN
              CALL splint (zk,fk,wk,km,depth(j,i),f2d(j,i),frstd)
            ELSE
              f2d(j,i)=spval
            END IF
          END DO
        END DO
      END IF
      RETURN
      END


!=======================================================================
!=======================================================================
!
      SUBROUTINE lintrp (n,x,y,ni,xi,yi)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Given arrays X and Y of length N, which tabulate a function,        !
!  Y = F(X),  with the Xs  in ascending order, and given array         !
!  XI of lenght NI, this routine returns a linear interpolated         !
!  array YI.                                                           !
!                                                                      !
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
      DO j=1,ni
        IF (xi(j).le.x(1)) THEN
          ii=1
          yi(j)=y(1)
        ELSE IF (xi(j).ge.x(n))then
          yi(j)=y(n)
        ELSE
          DO i=1,n-1
            IF ((x(i).lt.xi(j)).and.(xi(j).le.x(i+1))) THEN
              ii=i
              GO TO 10
            END IF
          END DO
 10       d1=xi(j)-x(ii)
          d2=x(ii+1)-xi(j)
          yi(j)=(d1*y(ii+1)+d2*y(ii))/(d1+d2)
        END IF
      END DO

      RETURN
      END


!=======================================================================
!=======================================================================
!
      SUBROUTINE spline (x,y,n,yp1,ypn,y2)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Given X, Y of length N containing a tabulated function,  Y=f(X),    !
!  with the Xs  in ascending order,  and given values  Yp1 and  Ypn    !
!  for the first derivative of the interpolating function at points    !
!  1 and N, respectively this routine returns an array Y2 of length    !
!  N  which contains the  second  derivatives of the  interpolating    !
!  function at the tabulated points X.  If Yp1 and/or Ypn are equal    !
!  to  1.0E+30  or larger,  the routine  is  signalled  to  set the    !
!  corresponding boundary condition for a natural spline, with zero    !
!  second derivative on that boundary.                                 !
!                                                                      !
!  Reference :                                                         !
!                                                                      !
!  Press, W.H, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,     !
!        1986: Numerical Recipes, the art of scientific computing.     !
!        Cambridge University Press.                                   !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define global data.
!-----------------------------------------------------------------------
!
      implicit none
      integer cntplt, secplt
      parameter (cntplt=1, secplt=2)
!
      real*8 cm1, cm3, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,&
     &     c20, c25, c50, c90, c100, c180, c200, c255, c300, c360, c366,&
     &     c500, c1000, c5000, c10000, c1em9, c1em10, c1em12, c1em20,   &
     &     c1ep30, p006, p009, p035, p015, p012, p08, p06, p5, p25, p75,&
     &     p98, r3, r10, r20, r33, r35, r40, r50, r80, r100, r200, r250,&
     &     r400, r1000
      real*8 day2sec, deg2rad, grav, cm2m, cm2mm, m2cm, m2km, m2knot,   &
     &     m2mm, m2smallknot, pi, rad2deg, re, root2, sec2day, spval0,  &
     &     spval1, spval2, spvgeo
!
      parameter (cm1=-1.0, cm3=-3.0, c0=0.0, c1=1.0, c2=2.0, c3=3.0,    &
     &           c4=4.0, c5=5.0, c6=6.0, c7=7.0, c8=8.0, c9=9.0, c10=10,&
     &           c11=11.0, c20=20.0, c25=25.0, c50=50.0, c90=90.0,      &
     &           c100=100.0, c180=180.0, c200=200.0, c255=255.0,        &
     &           c300=300.0, c360=360.0, c366=366.0, c500=500,          &
     &           c1000=1000.0, c5000=5000.0, c10000=10000.0,            &
     &           c1em9=1.0e-9, c1em10=1.0e-10, c1em12=1.0e-12,          &
     &           c1em20=1.0e-20, c1ep30=1.0e+30, p006=0.006, p009=0.009,&
     &           p012=0.012, p015=0.015, p035=0.035, p06=0.06, p08=0.08,&
     &           p5=0.5, p25=0.25, p75=0.75, p98=0.98, r3=c1/c3,        &
     &           r10=0.1, r20=0.05, r33=c1/33.0, r35=c1/35.0, r40=0.025,&
     &           r50=0.02, r80=0.0125, r100=0.01, r200=0.005,           &
     &           r250=0.004, r400=0.0025, r1000=0.001)
      parameter (day2sec=86400.0, cm2m=r100, cm2mm=r10, grav=9.801,     &
     &           m2cm=c100, m2km=r1000, m2knot=1.9459, m2mm=c1000,      &
     &           m2smallknot=m2knot*c5, pi=3.14159265358979323846,      &
     &           re=637131500.0, root2=1.41421356237309504880,          &
     &           sec2day=c1/86400.0, spvgeo=999.0, spval0=0.99e+35,     &
     &           spval1=1.0e+35, spval2=0.99e+30)
      parameter (deg2rad=pi/c180, rad2deg=c180/pi)
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
      IF (n.gt.nmax) THEN
        PRINT 10, n,nmax
 10     FORMAT (/' SPLINE: underdimensioned array, N, NMAX = ',2i5)
        !CALL crash ('SPLINE',1)
      END IF
!
!  The lower boundary condition is set either to be "natural" or else
!  to have a specified first derivative.
!
      IF (yp1.gt.spval2) THEN
        y2(1)=c0
        u(1)=c0
      ELSE
        y2(1)=-p5
        u(1)=(c3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      END IF
!
!  This is the decomposition loop of the tridiagonal algorithm. Y2 and
!  U are used for temporary storage of the decomposition factors.
!
      DO i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+c2
        y2(i)=(sig-c1)/p
        u(i)=(c6*((y(i+1)-y(i))/(x(i+1)-x(i))-                          &
     &           (y(i)-y(i-1))/(x(i)-x(i-1)))/                          &
     &           (x(i+1)-x(i-1))-sig*u(i-1))/p
      END DO
!
!  The upper boundary condition is set either to be "natural" or else
!  to have a specified first derivative.
!
      IF (ypn.gt.spval2) THEN
        qn=c0
        un=c0
      ELSE
        qn=p5
        un=(c3/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      END IF
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+c1)
!
!  This is the back-substitution loop of the tridiagonal algorithm.
!
      DO k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      END DO
      RETURN
      END



!=======================================================================
!=======================================================================
!
      SUBROUTINE splint (x,y,y2,n,xx,yy,dydx)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Given arrays X and Y of length N, which tabulate a function,        !
!  Y=f(X), with the Xs  in ascending order, and given the array        !
!  Y2 which contains the second derivative of the interpolating        !
!  function  at the  tabulated points X as computed  by routine        !
!  SPLINE, and given a value  XX, this routine returns a cubic-        !
!  spline interpolated value YY.                                       !
!                                                                      !
!  Reference :                                                         !
!                                                                      !
!  Press, W.H, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,     !
!         1986: Numerical Recipes, the art of scientific computing.    !
!         Cambridge University Press.                                  !
!                                                                      !
!  Modified by H.G. Arango (1989) to output the first derivative       !
!  DYDX at a given value XX.                                           !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define global data.
!-----------------------------------------------------------------------
!
      implicit none
      integer cntplt, secplt
      parameter (cntplt=1, secplt=2)
!
      real*8 cm1, cm3, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,&
     &     c20, c25, c50, c90, c100, c180, c200, c255, c300, c360, c366,&
     &     c500, c1000, c5000, c10000, c1em9, c1em10, c1em12, c1em20,   &
     &     c1ep30, p006, p009, p035, p015, p012, p08, p06, p5, p25, p75,&
     &     p98, r3, r10, r20, r33, r35, r40, r50, r80, r100, r200, r250,&
     &     r400, r1000
      real*8 day2sec, deg2rad, grav, cm2m, cm2mm, m2cm, m2km, m2knot,   &
     &     m2mm, m2smallknot, pi, rad2deg, re, root2, sec2day, spval0,  &
     &     spval1, spval2, spvgeo
!
      parameter (cm1=-1.0, cm3=-3.0, c0=0.0, c1=1.0, c2=2.0, c3=3.0,    &
     &           c4=4.0, c5=5.0, c6=6.0, c7=7.0, c8=8.0, c9=9.0, c10=10,&
     &           c11=11.0, c20=20.0, c25=25.0, c50=50.0, c90=90.0,      &
     &           c100=100.0, c180=180.0, c200=200.0, c255=255.0,        &
     &           c300=300.0, c360=360.0, c366=366.0, c500=500,          &
     &           c1000=1000.0, c5000=5000.0, c10000=10000.0,            &
     &           c1em9=1.0e-9, c1em10=1.0e-10, c1em12=1.0e-12,          &
     &           c1em20=1.0e-20, c1ep30=1.0e+30, p006=0.006, p009=0.009,&
     &           p012=0.012, p015=0.015, p035=0.035, p06=0.06, p08=0.08,&
     &           p5=0.5, p25=0.25, p75=0.75, p98=0.98, r3=c1/c3,        &
     &           r10=0.1, r20=0.05, r33=c1/33.0, r35=c1/35.0, r40=0.025,&
     &           r50=0.02, r80=0.0125, r100=0.01, r200=0.005,           &
     &           r250=0.004, r400=0.0025, r1000=0.001)
      parameter (day2sec=86400.0, cm2m=r100, cm2mm=r10, grav=9.801,     &
     &           m2cm=c100, m2km=r1000, m2knot=1.9459, m2mm=c1000,      &
     &           m2smallknot=m2knot*c5, pi=3.14159265358979323846,      &
     &           re=637131500.0, root2=1.41421356237309504880,          &
     &           sec2day=c1/86400.0, spvgeo=999.0, spval0=0.99e+35,     &
     &           spval1=1.0e+35, spval2=0.99e+30)
      parameter (deg2rad=pi/c180, rad2deg=c180/pi)
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
  10  IF ((khi-klo).gt.1) THEN
        k=(khi+klo)/2
        IF (x(k).gt.xx) THEN
          khi=k
        ELSE
          klo=k
        END IF
        GO TO 10
      END IF
!
!  KLO and KHI now bracket the input value XX.
!
      h=x(khi)-x(klo)
      IF (h.eq.c0) THEN
        PRINT *, ' SPLINT: bad X input, they must be distinct.'
        !CALL crash ('SPLINT',1)
      END IF
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
      RETURN
      END



!=======================================================================
!=======================================================================
!
      SUBROUTINE get_bottom (bottom,var,mask,im,jm,km,spval)
!
!
!=============================================== Frederic Castruccio ===

      integer i, im, j, jm, k, km
      real*8 spval, bot
      real*8 var(km,jm,im)
      real*8 mask(jm,im), bottom(jm,im)
!f2py intent(out) bottom
!f2py intent(hide) :: jm
!f2py intent(hide) :: im

      DO j=1,jm
        DO i=1,im
          k=1
          bot=spval
          IF (mask(j,i).eq.1) THEN
            DO WHILE ((bot.eq.spval).and.(k.le.km))
              IF (var(k,j,i).ne.spval) THEN
                bot=k
              ELSE
                k=k+1
              END IF
            END DO
          END IF
          ! python index start at 0...
          bottom(j,i)=bot-1
        END DO
      END DO

      RETURN
      END


!=======================================================================
!=======================================================================
!
      SUBROUTINE get_surface (surface,var,mask,im,jm,km,spval)
!
!
!=============================================== Frederic Castruccio ===

      integer i, im, j, jm, k, km
      real*8 spval, surf
      real*8 var(km,jm,im)
      real*8 mask(jm,im), surface(jm,im)
!f2py intent(out) surface
!f2py intent(hide) :: jm
!f2py intent(hide) :: im

      DO j=1,jm
        DO i=1,im
          k=km
          surf=spval
          IF (mask(j,i).eq.1) THEN
            DO WHILE ((surf.eq.spval).and.(k.ge.1))
              IF (var(k,j,i).ne.spval) THEN
                surf=k
              ELSE
                k=k-1
              END IF
            END DO
          END IF
          ! python index start at 0...
          surface(j,i)=surf-1
        END DO
      END DO

      RETURN
      END


