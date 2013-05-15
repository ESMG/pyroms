      subroutine flood(zslice, dry, wet_mask, delta,
     &                 flooded_zslice, im, jm, nbdry)

!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!
!     input arrays
!
!-----------------------------------------------------------------------

      integer :: delta

      integer, dimension(nbdry, 2), intent(in) :: dry

      real*8, dimension(jm, im), intent(in) :: zslice, wet_mask

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      real*8, dimension(jm, im), intent(out) :: flooded_zslice

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer :: im, jm, nbdry, n
      integer :: imin, imax, jmin, jmax, dxy
      integer :: idry, jdry

      real*8 :: wet_sum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      flooded_zslice = zslice

      ! loop over the dry point that need to be flooded
      do n=1,nbdry
        idry = dry(n,2)
        jdry = dry(n,1)
        ! compute wet point distance
        dxy = 0
        wet_sum = 0
        do while (wet_sum .eq. 0)
          dxy =  dxy + delta
          imin = MAX(idry-dxy,1)
          imax = MIN(idry+dxy,im)
          jmin = MAX(jdry-dxy,1)
          jmax = MIN(jdry+dxy,jm)
          wet_sum = SUM(wet_mask(jmin:jmax,imin:imax))
        enddo
        flooded_zslice(jdry,idry) = SUM(zslice(jmin:jmax,imin:imax)
     &                                  * wet_mask(jmin:jmax,imin:imax))
     &                                  / wet_sum
      enddo

!-----------------------------------------------------------------------

      end subroutine flood

!-----------------------------------------------------------------------
