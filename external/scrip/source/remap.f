!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this routine performs a remapping based on addresses and weights
!     computed in a setup phase
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap.f,v 1.5 2000/04/19 21:56:25 pwjones Exp $
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!***********************************************************************

      module remap_mod

!-----------------------------------------------------------------------
!
!     this module contains the routines for performing the actual
!     remappings
!
!-----------------------------------------------------------------------

      use kinds_mod    ! defines common data types
      use constants    ! defines common constants

      implicit none

!***********************************************************************

      contains

!***********************************************************************

      subroutine remap(dst_array, map_wts, dst_add, src_add, 
     &                 src_array, src_grad1, src_grad2, src_grad3)

!-----------------------------------------------------------------------
!
!     performs the remapping based on weights computed elsewhere
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input arrays
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), dimension(:), intent(in) ::
     &     dst_add,     ! destination address for each link
     &     src_add      ! source      address for each link

      real (kind=dbl_kind), dimension(:,:), intent(in) ::
     &     map_wts      ! remapping weights for each link

      real (kind=dbl_kind), dimension(:), intent(in) ::
     &     src_array    ! array with source field to be remapped

      real (kind=dbl_kind), dimension(:), intent(in), optional ::
     &     src_grad1    ! gradient arrays on source grid necessary for
     &,    src_grad2    ! higher-order remappings
     &,    src_grad3

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), dimension(:), intent(inout) ::
     &     dst_array    ! array for remapped field on destination grid

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n, iorder

!-----------------------------------------------------------------------
!
!     check the order of the interpolation
!
!-----------------------------------------------------------------------

      if (present(src_grad1)) then
        iorder = 2
      else
        iorder = 1
      endif

!-----------------------------------------------------------------------
!
!     first order remapping
!
!-----------------------------------------------------------------------

      dst_array = zero

      select case (iorder)
      case(1)

        do n=1,size(dst_add)
          dst_array(dst_add(n)) = dst_array(dst_add(n)) + 
     &                            src_array(src_add(n))*map_wts(1,n)
        end do

!-----------------------------------------------------------------------
!
!     second order remapping
!
!-----------------------------------------------------------------------

      case(2)

        if (size(map_wts,DIM=1) == 3) then
          do n=1,size(dst_add)
            dst_array(dst_add(n)) = dst_array(dst_add(n)) +
     &                              src_array(src_add(n))*map_wts(1,n) +
     &                              src_grad1(src_add(n))*map_wts(2,n) +
     &                              src_grad2(src_add(n))*map_wts(3,n)
          end do
        else if (size(map_wts,DIM=1) == 4) then
          do n=1,size(dst_add)
            dst_array(dst_add(n)) = dst_array(dst_add(n)) +
     &                              src_array(src_add(n))*map_wts(1,n) +
     &                              src_grad1(src_add(n))*map_wts(2,n) +
     &                              src_grad2(src_add(n))*map_wts(3,n) +
     &                              src_grad3(src_add(n))*map_wts(4,n)
          end do
        endif

      end select

!-----------------------------------------------------------------------

      end subroutine remap

!***********************************************************************

      end module remap_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
