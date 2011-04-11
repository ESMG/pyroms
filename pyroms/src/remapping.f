      subroutine flood(zslice, wet, dry, x, y, dmax,
     &                 flooded_zslice, im, jm, nbwet, nbdry)

!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!
!     input arrays
!
!-----------------------------------------------------------------------

      integer, dimension(nbwet, 2), intent(in) :: wet

      integer, dimension(nbdry, 2), intent(in) :: dry

      real*8, dimension(im, jm), intent(in) :: zslice, x, y

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      real*8, dimension(im, jm), intent(out) :: flooded_zslice

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer :: im, jm, nbwet, nbdry, n, m, dmax
      integer :: idry, jdry, iwet, jwet, dmin_idx, iclose, jclose

      real*8, dimension(nbwet) :: d

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      flooded_zslice = zslice

      ! loop over the dry point that need to be flooded
      do n=1,nbdry
        idry = dry(n,1)
        jdry = dry(n,2) 
        ! compute wet point distance
        do m=1,nbwet
          iwet = wet(m,1)
          jwet = wet(m,2)
          d(m) = sqrt( (x(idry,jdry) - x(iwet,jwet)) * 
     &                 (x(idry,jdry) - x(iwet,jwet))
     &               + (y(idry,jdry) - y(iwet,jwet)) * 
     &                 (y(idry,jdry) - y(iwet,jwet)) )
        enddo
        dmin_idx = minloc(d,1)
        if (dmax .eq. 0) then
          iclose = wet(dmin_idx,1)
          jclose = wet(dmin_idx,2)
          flooded_zslice(idry,jdry) = zslice(iclose,jclose)
        else
          if (d(dmin_idx) < dmax) then
            iclose = wet(dmin_idx,1)
            jclose = wet(dmin_idx,2)
            flooded_zslice(idry,jdry) = zslice(iclose,jclose)
          endif
        endif

      enddo

!-----------------------------------------------------------------------

      end subroutine flood

!-----------------------------------------------------------------------
