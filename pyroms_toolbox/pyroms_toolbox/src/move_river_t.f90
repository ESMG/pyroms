      subroutine move_river_t(runoff, pt, litpt, maskl, x, y, dx, dy,
     &                 adj_runoff, im, jm, nbpt, nblitpt)

!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!
!     input arrays
!
!-----------------------------------------------------------------------

      integer, dimension(nbpt, 2), intent(in) :: pt

      integer, dimension(nblitpt, 2), intent(in) :: litpt

      real*8, dimension(im, jm), intent(in) :: runoff, maskl
      real*8, dimension(im, jm), intent(in) :: x, y, dx, dy

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

      real*8 :: adj

      real*8, dimension(nblitpt) :: d

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      adj_runoff = runoff

      ! loop over the runoff point
      do n=1,nbpt
        ipt = pt(n,1)
        jpt = pt(n,2)
        ! point is not in the littoral band
        if (maskl(ipt, jpt).eq.0) then
          ! compute littoral point distance
          do m=1,nblitpt
            ilitpt = litpt(m,1)
            jlitpt = litpt(m,2)
            d(m) = sqrt( (x(ipt,jpt) - x(ilitpt,jlitpt)) * 
     &                   (x(ipt,jpt) - x(ilitpt,jlitpt))
     &                 + (y(ipt,jpt) - y(ilitpt,jlitpt)) * 
     &                   (y(ipt,jpt) - y(ilitpt,jlitpt)) )
          enddo
          dmin_idx = minloc(d,1)
          iclose = litpt(dmin_idx,1)
          jclose = litpt(dmin_idx,2)
!          ! we are conservative
!          adj = (runoff(ipt,jpt) * dx(ipt,jpt) * dy(ipt,jpt)) /
!     &            (dx(iclose,jclose) * dy(iclose,jclose))
!          adj_runoff(iclose,jclose) = adj_runoff(iclose,jclose) + adj

          ! we are not conservative
          adj_runoff(iclose,jclose) = runoff(ipt,jpt)
          adj_runoff(ipt, jpt) = 0
        endif
      enddo

!-----------------------------------------------------------------------

      end subroutine move_river_t

!-----------------------------------------------------------------------
