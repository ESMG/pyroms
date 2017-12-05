C File avgerage.f90
      subroutine avg3d(avg,dataset,incavg,counter,spval,dim1,dim2,dim3)

      IMPLICIT NONE

      !input variables
      INTEGER::dim1,dim2,dim3
      REAL*8,DIMENSION(dim1,dim2,dim3),intent(in)::dataset
      REAL*8,DIMENSION(dim1,dim2,dim3),intent(in)::incavg

      !output variables
      REAL*8,DIMENSION(dim1,dim2,dim3),intent(out)::avg

      !local variables
      REAL*8::counter
      REAL*8::part1,part2,spval
      INTEGER::i,j,k

      do i=1,dim1
         do j=1,dim2
            do k=1,dim3
               if (dataset(i,j,k).lt.spval) then
                 !This is the code for performing an incremental
                 !average broken into two parts. It follows the formula:
                 !avg(I+1) = avg(I) + ((data(I+1)-avg(I))/I+1)
                 part1 = dataset(i,j,k)-incavg(i,j,k)
                 part2 = part1 / counter
                 avg(i,j,k) = incavg(i,j,k)+part2
               else
                 !This keeps the missing value in the return array
                 avg(i,j,k) = dataset(i,j,k)
               endif
            enddo
         enddo
      enddo

      return

      end

      subroutine avg2d(avg,dataset,incavg,counter,spval,dim1,dim2)

      IMPLICIT NONE

      !input variables
      INTEGER::dim1,dim2
      REAL*8,DIMENSION(dim1,dim2),intent(in)::dataset
      REAL*8,DIMENSION(dim1,dim2),intent(in)::incavg

      !output variables
      REAL*8,DIMENSION(dim1,dim2),intent(out)::avg

      !local variables
      REAL*8::counter
      REAL*8::part1,part2,spval
      INTEGER::i,j

      do i=1,dim1
         do j=1,dim2
               if (dataset(i,j).lt.spval) then
                 !This is the code for performing an incremental
                 !average broken into two parts. It follows the formula:
                 !avg(I+1) = avg(I) + ((data(I+1)-avg(I))/I+1)
                 part1 = dataset(i,j)-incavg(i,j)
                 part2 = part1 / counter
                 avg(i,j) = incavg(i,j)+part2
               else
                 !This keeps the missing value in the return array
                 avg(i,j) = dataset(i,j)
               endif
         enddo
      enddo

      return

      end
