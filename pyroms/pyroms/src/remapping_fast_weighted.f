      subroutine flood(zslice, dry, wet_mask, delta,
     &                 flooded_zslice, im, jm, nbdry)

!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!
!     input arrays
!
!-----------------------------------------------------------------------

      integer, intent(in) :: delta

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
      integer :: is, ie, js, je

      real*8 :: wet_sum

      real*8, dimension(jm,im) :: binom_coef, mycoef

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
        binom_coef = 0.
        mycoef = 0.
        binom_coef(1:2*dxy+1,1:2*dxy+1) = binc(2*dxy)
        is = dxy - (idry-imin) + 1
        ie = is + (imax-imin)
        js = dxy - (jdry-jmin) + 1
        je = js + (jmax-jmin)
        mycoef(1:jmax-jmin+1,1:imax-imin+1) = binom_coef(js:je,is:ie)
        flooded_zslice(jdry,idry) = 
     &       SUM(zslice(jmin:jmax,imin:imax) 
     &         * wet_mask(jmin:jmax,imin:imax) 
     &         *  mycoef(1:jmax-jmin+1,1:imax-imin+1)) 
     &     / SUM(wet_mask(jmin:jmax,imin:imax) 
     &         * mycoef(1:jmax-jmin+1,1:imax-imin+1))
      enddo


!-----------------------------------------------------------------------

      contains

          function bin(n, k)

              implicit none
              integer, intent (in) :: n
              integer, intent (in) :: k
              integer :: i
              real*8 :: bin

              bin = 1

              do i=1,k
                  bin = bin * (n - i + 1) / i
              enddo

              return

          end function bin



          function binc(n)

              implicit none
              integer, intent (in) :: n
              integer :: i,j
              real*8, dimension(:,:), allocatable :: binc

              allocate(binc(n+1,n+1))

              do j=0,n
                  do i=0,n
                      binc(i+1,j+1) = bin(n,j) * bin(n,i)
                  enddo
              enddo

              return

          end function binc


!-----------------------------------------------------------------------

      end subroutine flood

!-----------------------------------------------------------------------

