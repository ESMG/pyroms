
SUBROUTINE cslf(varin,spval,validmin,validmax,nx,ny,varout)

  IMPLICIT NONE

  REAL(8),DIMENSION(nx,ny),INTENT(IN)  :: varin
  REAL(8),DIMENSION(nx,ny),INTENT(OUT) :: varout
  INTEGER,INTENT(IN)                   :: nx,ny
  REAL(8),INTENT(IN)                   :: spval,validmin,validmax

  REAL(8),DIMENSION(nx,ny)             :: ztmp
  REAL(8)                              :: zval, minvarin, maxvarin
  INTEGER,DIMENSION(nx,ny)             :: zmask, zmask_save
  INTEGER,DIMENSION(0:nx+1,0:ny+1)         :: zmask_halo2pt
  REAL,DIMENSION(0:nx+1,0:ny+1)            :: ztmp_halo2pt
  INTEGER                              :: jj, ji    ! loop index
  INTEGER                              :: jjp1, jjm1, jip1, jim1 ! neighbors
  INTEGER                              :: ctot, c1, c2, c3, c4, c5, c6, c7, c8, ct
  INTEGER                              :: nt, nmax=3000 ! customizable

  ztmp(:,:)  = varin(:,:)
  minvarin   = MINVAL(varin(:,:))
  maxvarin   = MAXVAL(varin(:,:))

  ! Define mask
  zmask(:,:) = 1
  WHERE( varin == spval    ) zmask = 0
  WHERE( varin <= validmin ) zmask = 0
  WHERE( varin >= validmax ) zmask = 0

  ! save for smoother use
  zmask_save = zmask

  ! Define arrays with 1 point halo...
  zmask_halo2pt(:,:) = 0
  ztmp_halo2pt(:,:)  = spval
  ! and setting values
  ztmp_halo2pt(1:nx,1:ny)  = varin(:,:)
  zmask_halo2pt(1:nx,1:ny) = zmask(:,:)

  ! Here would be the right place for an MPI call
  ! to fill halo with values from neighbors (if any)

  ! c1 -- c2 -- c3
  ! |
  ! c4          c5
  ! |
  ! c6 -- c7 -- c8


  nt = 0

  DO WHILE ( ANY(zmask(:,:) == 0 ) .AND. ( nt < nmax ) )

     DO jj=1,ny
       DO ji=1,nx

          ! compute these indices once
          jjm1 = jj-1 ; jjp1 = jj+1
          jim1 = ji-1 ; jip1 = ji+1

          IF (zmask(ji,jj) == 0 ) THEN

             ! reads along the fast axis
             c6 = 1 * zmask_halo2pt(jim1,jjm1)
             c7 = 2 * zmask_halo2pt(ji  ,jjm1)
             c8 = 1 * zmask_halo2pt(jip1,jjm1)

             c4 = 2 * zmask_halo2pt(jim1,jj  )
             c5 = 2 * zmask_halo2pt(jip1,jj  )

             c1 = 1 * zmask_halo2pt(jim1,jjp1)
             c2 = 2 * zmask_halo2pt(ji  ,jjp1)
             c3 = 1 * zmask_halo2pt(jip1,jjp1)

             ctot = c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8

             IF (ctot >= 3 ) THEN
                ! compute the new value for this point
                    zval  = ( c6 * ztmp_halo2pt(jim1,jjm1) + &
               &              c7 * ztmp_halo2pt(ji  ,jjm1) + &
               &              c8 * ztmp_halo2pt(jip1,jjm1) + &
               &              c4 * ztmp_halo2pt(jim1,jj  ) + &
               &              c5 * ztmp_halo2pt(jip1,jj  ) + &
               &              c1 * ztmp_halo2pt(jim1,jjp1) + &
               &              c2 * ztmp_halo2pt(ji  ,jjp1) + &
               &              c3 * ztmp_halo2pt(jip1,jjp1) ) / &
               &              ( ctot )

                ! update value in field array
                ztmp(ji,jj) = zval
                ! set the mask to sea
                zmask(ji,jj) = 1


             ENDIF

          ENDIF
       ENDDO
     ENDDO
     DO jj=1,ny
       DO ji=1,nx
          ztmp_halo2pt(ji,jj) = ztmp(ji,jj)
          zmask_halo2pt(ji,jj) = zmask(ji,jj)
       ENDDO
     ENDDO

  nt = nt + 1

  ENDDO

  IF (nt >= nmax ) THEN
     PRINT *, 'WARNING: flood did not converge after number of iteration = ', nmax
     PRINT *, 'you should increase nmax in the flooding routine'
  ENDIF

  ztmp(:,:) = ztmp_halo2pt(1:nx,1:ny)

!RD: work on this is in progress
  !CALL smooth(ztmp,zmask_save,nx,ny,varout)

  ! bound the values with min/max of the input
  WHERE( ztmp < minvarin ) ztmp = minvarin
  WHERE( ztmp > maxvarin ) ztmp = maxvarin

  varout = ztmp

END SUBROUTINE

SUBROUTINE smooth(varin,mask,nx,ny,varout)

  IMPLICIT NONE

  REAL(8),DIMENSION(nx,ny),INTENT(in)   :: varin
  INTEGER,DIMENSION(nx,ny),INTENT(in)   :: mask
  REAL(8),DIMENSION(nx,ny),INTENT(out)  :: varout
  INTEGER,INTENT(in)                    :: nx,ny

  REAL(8), DIMENSION(nx,ny)             :: ztmp, ztmp2
  REAL(8)                               :: zval
  INTEGER                               :: n=0, npass=1000
  INTEGER                               :: ji, jj
  REAL(8)                               :: c1, c2, ctot

  c1 = 1.
  c2 = 2.

  ztmp = varin

  !PRINT *, 'Entering smoother'
  !PRINT *, MINVAL(mask) , MAXVAL(mask)

  n = 0
  DO WHILE ( n < npass )

     DO jj=2,ny-1
       DO ji=2,nx-1

       !   ztmp2(ji,jj) = c1 * ztmp(ji-1,jj-1) + &
       !  &               c1 * ztmp(ji,  jj-1) + &
       !  &              c1 * ztmp(ji+1,jj-1) + &
       !  &              c1 * ztmp(ji-1,jj  ) + &
       !  &              c1 * ztmp(ji  ,jj  ) + &
       !  &              c1 * ztmp(ji+1,jj  ) + &
       !  &              c1 * ztmp(ji-1,jj+1) + &
       !  &              c1 * ztmp(ji  ,jj+1) + &
       !  &              c1 * ztmp(ji+1,jj+1)

          ! direct neighbor have a superior weight / nonmasked points also
          zval  = c1 * ztmp(ji-1,jj-1) * ( 1 + mask(ji-1,jj-1) ) + &
         &        c2 * ztmp(ji,  jj-1) * ( 1 + mask(ji,  jj-1) ) + &
         &        c1 * ztmp(ji+1,jj-1) * ( 1 + mask(ji+1,jj-1) ) + &
         &        c2 * ztmp(ji-1,jj  ) * ( 1 + mask(ji-1,jj  ) ) + &
         &        c2 * ztmp(ji  ,jj  ) * ( 1 + mask(ji  ,jj  ) ) + &
         &        c2 * ztmp(ji+1,jj  ) * ( 1 + mask(ji+1,jj  ) ) + &
         &        c1 * ztmp(ji-1,jj+1) * ( 1 + mask(ji-1,jj+1) ) + &
         &        c2 * ztmp(ji  ,jj+1) * ( 1 + mask(ji  ,jj+1) ) + &
         &        c1 * ztmp(ji+1,jj+1) * ( 1 + mask(ji+1,jj+1) )

          ctot  = c1 * ( 1 + mask(ji-1,jj-1) ) + &
         &        c2 * ( 1 + mask(ji,  jj-1) ) + &
         &        c1 * ( 1 + mask(ji+1,jj-1) ) + &
         &        c2 * ( 1 + mask(ji-1,jj  ) ) + &
         &        c2 * ( 1 + mask(ji  ,jj  ) ) + &
         &        c2 * ( 1 + mask(ji+1,jj  ) ) + &
         &        c1 * ( 1 + mask(ji-1,jj+1) ) + &
         &        c2 * ( 1 + mask(ji  ,jj+1) ) + &
         &        c1 * ( 1 + mask(ji+1,jj+1) )

          ztmp2(ji,jj) = zval / ctot

       ENDDO
     ENDDO

  WHERE( mask == 1 ) ztmp2 = varin
  ztmp(2:nx-1,2:ny-1) = ztmp2(2:nx-1,2:ny-1)

  n = n + 1
  !PRINT *, 'smoother iter = ', n
  ENDDO

  varout = ztmp


END SUBROUTINE
