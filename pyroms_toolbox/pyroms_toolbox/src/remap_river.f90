      subroutine remap_river(runoff, x, y, litpt, litx, lity,           &
     &                 adj_runoff, im, jm, nbpt, nblitpt)

!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!
!     input arrays
!
!-----------------------------------------------------------------------

      real*8, dimension(nbpt), intent(in) :: runoff
      real*8, dimension(nbpt), intent(in) :: x, y

      integer, dimension(nblitpt, 2), intent(in) :: litpt
      real*8, dimension(im, jm), intent(in) :: litx, lity

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      real*8, dimension(im, jm), intent(out) :: adj_runoff

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer :: im, jm, nbpt, nblitpt, n, m
      integer :: ipt, jpt, ilitpt, jlitpt, dmin_idx, iclose, jclose

      real*8, dimension(nblitpt) :: d
      real*8 :: alph
      integer :: flag = 1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      adj_runoff = 0.0

      ! loop over the runoff point
      do n=1,nbpt
        ! point is not in the littoral band
        ! compute littoral point distance
        do m=1,nblitpt
          ilitpt = litpt(m,1)
          jlitpt = litpt(m,2)
!         d(m) = sqrt( (x(n) - litx(ilitpt,jlitpt)) *                   &
!    &                 (x(n) - litx(ilitpt,jlitpt)) *                   &
!    &               + (y(n) - lity(ilitpt,jlitpt)) *                   &
!    &                 (y(n) - lity(ilitpt,jlitpt)) )
          call geodesic_dist(x(n), y(n),                                &
     &             litx(ilitpt,jlitpt), lity(ilitpt,jlitpt),            &
     &             flag, d(m), alph)
        enddo
        dmin_idx = minloc(d,1)
        iclose = litpt(dmin_idx,1)
        jclose = litpt(dmin_idx,2)

        adj_runoff(iclose,jclose) = adj_runoff(iclose,jclose) +       &
     &            runoff(n)
!       if (runoff(n) > 5000.0) then
!         print *, "found something to move", iclose, jclose
!         print *, "after in move_river_t", runoff(n),                &
!    &         adj_runoff(iclose, jclose)
!       endif
      enddo

!-----------------------------------------------------------------------

      end subroutine remap_river

!-----------------------------------------------------------------------
      subroutine geodesic_dist (lon1,lat1,lon2,lat2,flag,dist,alpha)
!
!=======================================================================
!                                                                     ==
!  Inverse, non-iterative solutions for distance and geodesic azimuth ==
!  between two points on the ellipsoid (The Earth) from the equations,==
!  second order in spheroidal flatttening, given by:                  ==
!                                                                     ==
!  Sodano , E.M., and T. Robinson, 1963: Direct and inverse solutions ==
!   of geodesics, Army Map Service Technical Report No. 7, AD 657591. ==
!                                                                     ==
!  On input:    Longitude is positive to the east and negative to the ==
!               west.  Latitude is positive to the north and negative ==
!               to the south.                                         ==
!                                                                     ==
!     LON1    Longitude point 1 (decimal degrees, real*8)             ==
!     LAT1    Latitude  point 1 (decimal degrees, real*8)             ==
!     LON2    Longitude point 2 (decimal degrees, real*8)             ==
!     LAT2    Latitude  point 2 (decimal degrees, real*8)             ==
!     FLAG    flag for distance units on output (integer)             ==
!                                                                     ==
!  On output:                                                         ==
!                                                                     ==
!     GALPHA  Geodesic azimuth from point 1 to point 2 clockwise from ==
!             North (decimal degrees, real*8)                         ==
!     GDIST   Geodesic distance between point 1 and point 2 (real*8)  ==
!                                                                     ==
!             Units of distance                                       ==
!                                                                     ==
!               Flag    Units                                         ==
!              ------  -------                                        ==
!                 1     Meters                                        ==
!                 2     Nautical Miles                                ==
!                 3     Feet                                          ==
!                 4     Kilometers                                    ==
!                 5     Statute Mile                                  ==
!                                                                     ==
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define local data
!-----------------------------------------------------------------------
!
      logical first
      integer flag
      real*8 alpha, dist, lat1, lat2, lon1, lon2, r_lat1, r_lat2, &
     &        delta, l, beta1, beta2, a, b, c, ct, st, t, m, sob, &
     &        lambda, cott, adist
      real*8 c0, c1, c4, c90, c180, c360, deg2rad, pi, rad2deg, smin
      parameter (c0=0.d0, c1=1.0d0, c4=4.0d0, c90=90.0d0, &
     &        c180=180.0d0, c360=360.0d0  )
!#if ELLIPSOID
!      real*8 f, smaj, q, w, x, y, z, p5, p0625, p125, p25, c5
!      parameter (c5=5.0d0, p5=0.5d0, p25=0.25d0, p125=0.125d0, &
!     &        p0625=0.0625d0  )
!      save    f, smaj
!#endif  /* ELLIPSOID */
      save deg2rad, pi, rad2deg, smin
      data first /.true./
!
!=======================================================================
!  Begin executable code.
!=======================================================================
!
!  Define parameters on first pass (SMIN: Ellipsoid semi-minor axis in
!  meters; SMAJ: Ellipsoid semi-major axis in meters; F: spheroidal
!  flattening).
!
      if (first) then
        smin = 6356750.52d0
!#if ELLIPSOID
!        smaj = 6378135.0d0
!        f = c1-(smin/smaj)
!#endif  /* ELLIPSOID */
        pi = c4*atan(c1)
        deg2rad = pi/c180
        rad2deg = c180/pi
        first = .false.
      endif
!
!  Determine proper longitudinal shift.
!
      delta = lon2-lon1
      l = abs(delta)
      if (l.ge.c180) l = c360-abs(lon1-lon2)
!
!  Convert Decimal degrees to radians.
!
      r_lat1 = lat1*deg2rad
      r_lat2 = lat2*deg2rad
      l = l*deg2rad
!
!  Calculate S/Bo subformulas.
!
!#if ELLIPSOID
!      beta1 = atan(tan(r_lat1)*(c1-f))
!      beta2 = atan(tan(r_lat2)*(c1-f))
!      a = sin(beta1)*sin(beta2)
!      b = cos(beta1)*cos(beta2)
!      ct = a+b*cos(l)
!      st = sqrt(((sin(l)*cos(beta2))**2)+(((sin(beta2)*cos(beta1)) &
!     &           -(sin(beta1)*cos(beta2)*cos(l)))**2))
!      t = asin(st)
!      c = (b*sin(l))/st
!      m = c1-(c*c)
!#else
      beta1 = r_lat1
      beta2 = r_lat2
      a = sin(beta1)*sin(beta2)
      b = cos(beta1)*cos(beta2)
      ct = a+b*cos(l)
      st = sqrt(((sin(l)*cos(beta2))**2)+(((sin(beta2)*cos(beta1)) &
     &           -(sin(beta1)*cos(beta2)*cos(l)))**2))
      t = asin(st)
      c = (b*sin(l))/st
      m = c1-(c*c)
!#endif  /* ELLIPSOID */
!
!  Calculate S/Bo term.
!
!#if ELLIPSOID
!      q = f+(f*f)
!      z = f*f*p5
!      x = f*f*p0625
!      y = f*f*p125
!      w = f*f*p25
!!
!      sob = ((c1+q)*t)+(a*((q*st)-(z*(t*t)*(c1/sin(t))))) &
!     &    +(m*(((-q*p5)*t)-((q*p5)*st*ct)+(z*(t*t)*(c1/tan(t))))) &
!     &    +((a**2)*(-z*st*ct)) &
!     &    +((m**2)*((x*t)+(x*st*ct)-(z*(t*t)*(c1/tan(t))) &
!     &    -(y*st*(ct*ct*ct)))) &
!     &    +((a*m)*((z*(t*t)*(c1/sin(t)))+(z*st*(ct*ct))))
!#else
      sob = t
!#endif  /* ELLIPSOID */
!
!  Compute geodesic azimuth from point 1 to point 2 clockwise from
!  North, alpha.
!
!#if ELLIPSOID
!      lambda = q*t+a*(-z*st-f*f*t*t/sin(t)) &
!     &          +m*(-c5*w*t+w*st*cos(t)+f*f*t*t/tan(t))
!      lambda = c*lambda+l
!#else
      lambda = l
!#endif  /* ELLIPSOID */
      if (lambda.eq.c0) then
        if (lat1.lt.lat2) alpha = c0
        if (lat1.gt.lat2) alpha = c180
        goto 1
      endif

      cott = (sin(beta2)*cos(beta1)- &
     &      cos(lambda)*sin(beta1)*cos(beta2))/ &
     &     (sin(lambda)*cos(beta2))
      if (cott.eq.c0) then
        alpha = c90
      else
        alpha = atan(c1/cott)*rad2deg
      endif

!       Compute heading from point#1 to point#2 clockwise from north

        if (delta .gt. 0.d0) then
            if (cott .gt. 0.d0) then
                alpha = alpha                 ! first quadrant
            else if (cott .lt. 0.d0) then
                alpha = 180.d0 + alpha        ! second quadrant
            end if
        end if
        if (delta .lt. 0.d0) then
            if (cott .lt. 0.d0) then
                alpha = 180.d0 - alpha        ! third quadrant
            else if (cott .gt. 0.d0) then
                alpha = 360.d0 - alpha        ! fourth quadrant
            end if
        end if

!       Calculate distance from point 1 to point 2

   1    adist = sob * smin

!       Check flag for proper output units

        if (flag .eq. 1) dist = adist                 ! meters
        if (flag .eq. 2) dist = adist * 5.396d-4      ! nautical miles
        if (flag .eq. 3) dist = adist * 3.281d0       ! feet
        if (flag .eq. 4) dist = adist * 1.d-3         ! kilometers
        if (flag .eq. 5) dist = adist * 6.214d-4      ! statute mile

        return
        end

