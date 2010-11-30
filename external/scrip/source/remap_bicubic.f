!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains necessary routines for performing an 
!     bicubic interpolation.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_bicubic.f,v 1.5 2001/08/22 18:20:41 pwjones Exp $
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

      module remap_bicubic

!-----------------------------------------------------------------------

      use kinds_mod     ! defines common data types
      use constants     ! defines common constants
      use grids         ! module containing grid info
      use remap_vars    ! module containing remap info

      implicit none

!-----------------------------------------------------------------------

      integer (kind=int_kind), parameter ::
     &    max_iter = 100   ! max iteration count for i,j iteration

      real (kind=dbl_kind), parameter ::
     &     converge = 1.e-10_dbl_kind  ! convergence criterion

!***********************************************************************

      contains

!***********************************************************************

      subroutine remap_bicub

!-----------------------------------------------------------------------
!
!     this routine computes the weights for a bicubic interpolation.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n,icount,
     &     dst_add,        ! destination address
     &     iter,           ! iteration counter
     &     nmap            ! index of current map being computed

      integer (kind=int_kind), dimension(4) :: 
     &     src_add         ! address for the four source points

      real (kind=dbl_kind), dimension(4)  ::
     &     src_lats,       ! latitudes  of four bilinear corners
     &     src_lons        ! longitudes of four bilinear corners

      real (kind=dbl_kind), dimension(4,4)  ::
     &     wgts            ! bicubic weights for four corners

      real (kind=dbl_kind) ::
     &     plat, plon,       ! lat/lon coords of destination point
     &     iguess, jguess,   ! current guess for bilinear coordinate
     &     thguess, phguess, ! current guess for lat/lon coordinate
     &     deli, delj,       ! corrections to i,j
     &     dth1, dth2, dth3, ! some latitude  differences
     &     dph1, dph2, dph3, ! some longitude differences
     &     dthp, dphp,       ! difference between point and sw corner
     &     mat1, mat2, mat3, mat4, ! matrix elements
     &     determinant,      ! matrix determinant
     &     sum_wgts,         ! sum of weights for normalization
     &     w1,w2,w3,w4,w5,w6,w7,w8, ! 16 bicubic weight functions
     &     w9,w10,w11,w12,w13,w14,w15,w16

!-----------------------------------------------------------------------
!
!     compute mappings from grid1 to grid2
!
!-----------------------------------------------------------------------

      nmap = 1
      if (grid1_rank /= 2) then
        stop 'Can not do bicubic interpolation when grid_rank /= 2'
      endif

      !***
      !*** loop over destination grid 
      !***

      grid_loop1: do dst_add = 1, grid2_size

        if (.not. grid2_mask(dst_add)) cycle grid_loop1

        plat = grid2_center_lat(dst_add)
        plon = grid2_center_lon(dst_add)

!-----------------------------------------------------------------------
!
!       find nearest square of grid points on source grid
!
!-----------------------------------------------------------------------

        call grid_search_bicub(src_add, src_lats, src_lons, 
     &                         plat, plon, grid1_dims,
     &                         grid1_center_lat, grid1_center_lon,
     &                         grid1_bound_box, bin_addr1, bin_addr2)

        !***
        !*** check to see if points are land points
        !***

        if (src_add(1) > 0) then
          do n=1,4
            if (.not. grid1_mask(src_add(n))) src_add(1) = 0
          end do
        endif

!-----------------------------------------------------------------------
!
!       if point found, find local i,j coordinates for weights
!
!-----------------------------------------------------------------------

        if (src_add(1) > 0) then

          grid2_frac(dst_add) = one

          !***
          !*** iterate to find i,j for bicubic approximation
          !***

          dth1 = src_lats(2) - src_lats(1)
          dth2 = src_lats(4) - src_lats(1)
          dth3 = src_lats(3) - src_lats(2) - dth2

          dph1 = src_lons(2) - src_lons(1)
          dph2 = src_lons(4) - src_lons(1)
          dph3 = src_lons(3) - src_lons(2)

          if (dph1 >  three*pih) dph1 = dph1 - pi2
          if (dph2 >  three*pih) dph2 = dph2 - pi2
          if (dph3 >  three*pih) dph3 = dph3 - pi2
          if (dph1 < -three*pih) dph1 = dph1 + pi2
          if (dph2 < -three*pih) dph2 = dph2 + pi2
          if (dph3 < -three*pih) dph3 = dph3 + pi2

          dph3 = dph3 - dph2

          iguess = half
          jguess = half

          iter_loop1: do iter=1,max_iter

            dthp = plat - src_lats(1) - dth1*iguess -
     &                    dth2*jguess - dth3*iguess*jguess
            dphp = plon - src_lons(1)

            if (dphp >  three*pih) dphp = dphp - pi2
            if (dphp < -three*pih) dphp = dphp + pi2

            dphp = dphp - dph1*iguess - dph2*jguess - 
     &                    dph3*iguess*jguess

            mat1 = dth1 + dth3*jguess
            mat2 = dth2 + dth3*iguess
            mat3 = dph1 + dph3*jguess
            mat4 = dph2 + dph3*iguess

            determinant = mat1*mat4 - mat2*mat3

            deli = (dthp*mat4 - mat2*dphp)/determinant
            delj = (mat1*dphp - dthp*mat3)/determinant

            if (abs(deli) < converge .and. 
     &          abs(delj) < converge) exit iter_loop1

            iguess = iguess + deli
            jguess = jguess + delj

          end do iter_loop1

          if (iter <= max_iter) then

!-----------------------------------------------------------------------
!
!           successfully found i,j - compute weights
!
!-----------------------------------------------------------------------

            wgts(1,1) = (one - jguess**2*(three-two*jguess))*
     &                  (one - iguess**2*(three-two*iguess))
            wgts(1,2) = (one - jguess**2*(three-two*jguess))*
     &                         iguess**2*(three-two*iguess)
            wgts(1,3) =        jguess**2*(three-two*jguess)*
     &                         iguess**2*(three-two*iguess)
            wgts(1,4) =        jguess**2*(three-two*jguess)*
     &                  (one - iguess**2*(three-two*iguess))
            wgts(2,1) = (one - jguess**2*(three-two*jguess))*
     &                         iguess*(iguess-one)**2
            wgts(2,2) = (one - jguess**2*(three-two*jguess))*
     &                         iguess**2*(iguess-one)
            wgts(2,3) =        jguess**2*(three-two*jguess)*
     &                         iguess**2*(iguess-one)
            wgts(2,4) =        jguess**2*(three-two*jguess)*
     &                         iguess*(iguess-one)**2
            wgts(3,1) =        jguess*(jguess-one)**2*
     &                  (one - iguess**2*(three-two*iguess))
            wgts(3,2) =        jguess*(jguess-one)**2*
     &                         iguess**2*(three-two*iguess)
            wgts(3,3) =        jguess**2*(jguess-one)*
     &                         iguess**2*(three-two*iguess)
            wgts(3,4) =        jguess**2*(jguess-one)*
     &                  (one - iguess**2*(three-two*iguess))
            wgts(4,1) =        iguess*(iguess-one)**2*
     &                         jguess*(jguess-one)**2
            wgts(4,2) =        iguess**2*(iguess-one)*
     &                         jguess*(jguess-one)**2
            wgts(4,3) =        iguess**2*(iguess-one)*
     &                         jguess**2*(jguess-one)
            wgts(4,4) =        iguess*(iguess-one)**2*
     &                         jguess**2*(jguess-one)

            call store_link_bicub(dst_add, src_add, wgts, nmap)

          else
            stop 'Iteration for i,j exceed max iteration count'
          endif

!-----------------------------------------------------------------------
!
!       search for bilinear failed - use a distance-weighted
!       average instead (this is typically near the pole)
!
!-----------------------------------------------------------------------

        else if (src_add(1) < 0) then

          src_add = abs(src_add)

          icount = 0
          do n=1,4
            if (grid1_mask(src_add(n))) then
              icount = icount + 1
            else
              src_lats(n) = zero
            endif
          end do

          if (icount > 0) then
            !*** renormalize weights

            sum_wgts = sum(src_lats)
            wgts(1,1) = src_lats(1)/sum_wgts
            wgts(1,2) = src_lats(2)/sum_wgts
            wgts(1,3) = src_lats(3)/sum_wgts
            wgts(1,4) = src_lats(4)/sum_wgts
            wgts(2:4,:) = zero

            grid2_frac(dst_add) = one
            call store_link_bicub(dst_add, src_add, wgts, nmap)
          endif

        endif
      end do grid_loop1

!-----------------------------------------------------------------------
!
!     compute mappings from grid2 to grid1 if necessary
!
!-----------------------------------------------------------------------

      if (num_maps > 1) then

      nmap = 2
      if (grid2_rank /= 2) then
        stop 'Can not do bicubic interpolation when grid_rank /= 2'
      endif

      !***
      !*** loop over destination grid 
      !***

      grid_loop2: do dst_add = 1, grid1_size

        if (.not. grid1_mask(dst_add)) cycle grid_loop2

        plat = grid1_center_lat(dst_add)
        plon = grid1_center_lon(dst_add)

        !***
        !*** find nearest square of grid points on source grid
        !***

        call grid_search_bicub(src_add, src_lats, src_lons, 
     &                         plat, plon, grid2_dims,
     &                         grid2_center_lat, grid2_center_lon,
     &                         grid2_bound_box, bin_addr2, bin_addr1)

        !***
        !*** check to see if points are land points
        !***

        if (src_add(1) > 0) then
          do n=1,4
            if (.not. grid2_mask(src_add(n))) src_add(1) = 0
          end do
        endif

        !***
        !*** if point found, find i,j coordinates for weights
        !***

        if (src_add(1) > 0) then

          grid1_frac(dst_add) = one

          !***
          !*** iterate to find i,j for bilinear approximation
          !***

          dth1 = src_lats(2) - src_lats(1)
          dth2 = src_lats(4) - src_lats(1)
          dth3 = src_lats(3) - src_lats(2) - dth2

          dph1 = src_lons(2) - src_lons(1)
          dph2 = src_lons(4) - src_lons(1)
          dph3 = src_lons(3) - src_lons(2)

          if (dph1 >  pi) dph1 = dph1 - pi2
          if (dph2 >  pi) dph2 = dph2 - pi2
          if (dph3 >  pi) dph3 = dph3 - pi2
          if (dph1 < -pi) dph1 = dph1 + pi2
          if (dph2 < -pi) dph2 = dph2 + pi2
          if (dph3 < -pi) dph3 = dph3 + pi2

          dph3 = dph3 - dph2

          iguess = zero
          jguess = zero

          iter_loop2: do iter=1,max_iter

            dthp = plat - src_lats(1) - dth1*iguess -
     &                    dth2*jguess - dth3*iguess*jguess
            dphp = plon - src_lons(1)

            if (dphp >  pi) dphp = dphp - pi2
            if (dphp < -pi) dphp = dphp + pi2

            dphp = dphp - dph1*iguess - dph2*jguess - 
     &                    dph3*iguess*jguess

            mat1 = dth1 + dth3*jguess
            mat2 = dth2 + dth3*iguess
            mat3 = dph1 + dph3*jguess
            mat4 = dph2 + dph3*iguess

            determinant = mat1*mat4 - mat2*mat3

            deli = (dthp*mat4 - mat2*dphp)/determinant
            delj = (mat1*dphp - dthp*mat3)/determinant

            if (abs(deli) < converge .and. 
     &          abs(delj) < converge) exit iter_loop2

            iguess = iguess + deli
            jguess = jguess + delj

          end do iter_loop2

          if (iter <= max_iter) then

            !***
            !*** successfully found i,j - compute weights
            !***

            wgts(1,1) = (one - jguess**2*(three-two*jguess))*
     &                  (one - iguess**2*(three-two*iguess))
            wgts(1,2) = (one - jguess**2*(three-two*jguess))*
     &                         iguess**2*(three-two*iguess)
            wgts(1,3) =        jguess**2*(three-two*jguess)*
     &                         iguess**2*(three-two*iguess)
            wgts(1,4) =        jguess**2*(three-two*jguess)*
     &                  (one - iguess**2*(three-two*iguess))
            wgts(2,1) = (one - jguess**2*(three-two*jguess))*
     &                         iguess*(iguess-one)**2
            wgts(2,2) = (one - jguess**2*(three-two*jguess))*
     &                         iguess**2*(iguess-one)
            wgts(2,3) =        jguess**2*(three-two*jguess)*
     &                         iguess**2*(iguess-one)
            wgts(2,4) =        jguess**2*(three-two*jguess)*
     &                         iguess*(iguess-one)**2
            wgts(3,1) =        jguess*(jguess-one)**2*
     &                  (one - iguess**2*(three-two*iguess))
            wgts(3,2) =        jguess*(jguess-one)**2*
     &                         iguess**2*(three-two*iguess)
            wgts(3,3) =        jguess**2*(jguess-one)*
     &                         iguess**2*(three-two*iguess)
            wgts(3,4) =        jguess**2*(jguess-one)*
     &                  (one - iguess**2*(three-two*iguess))
            wgts(4,1) =        iguess*(iguess-one)**2*
     &                         jguess*(jguess-one)**2
            wgts(4,2) =        iguess**2*(iguess-one)*
     &                         jguess*(jguess-one)**2
            wgts(4,3) =        iguess**2*(iguess-one)*
     &                         jguess**2*(jguess-one)
            wgts(4,4) =        iguess*(iguess-one)**2*
     &                         jguess**2*(jguess-one)

            call store_link_bicub(dst_add, src_add, wgts, nmap)

          else
            stop 'Iteration for i,j exceed max iteration count'
          endif

        !***
        !*** search for bilinear failed - us a distance-weighted
        !*** average instead
        !***

        else if (src_add(1) < 0) then

          src_add = abs(src_add)

          icount = 0
          do n=1,4
            if (grid2_mask(src_add(n))) then
              icount = icount + 1
            else
              src_lats(n) = zero
            endif
          end do

          if (icount > 0) then
            !*** renormalize weights

            sum_wgts = sum(src_lats)
            wgts(1,1) = src_lats(1)/sum_wgts
            wgts(1,2) = src_lats(2)/sum_wgts
            wgts(1,3) = src_lats(3)/sum_wgts
            wgts(1,4) = src_lats(4)/sum_wgts
            wgts(2:4,:) = zero

            grid1_frac(dst_add) = one
            call store_link_bicub(dst_add, src_add, wgts, nmap)
          endif

        endif
      end do grid_loop2

      endif ! nmap=2

!-----------------------------------------------------------------------

      end subroutine remap_bicub

!***********************************************************************

      subroutine grid_search_bicub(src_add, src_lats, src_lons, 
     &                             plat, plon, src_grid_dims,
     &                             src_center_lat, src_center_lon,
     &                             src_bound_box,
     &                             src_bin_add, dst_bin_add)

!-----------------------------------------------------------------------
!
!     this routine finds the location of the search point plat, plon
!     in the source grid and returns the corners needed for a bicubic
!     interpolation.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), dimension(4), intent(out) ::
     &        src_add  ! address of each corner point enclosing P

      real (kind=dbl_kind), dimension(4), intent(out) ::
     &        src_lats, ! latitudes  of the four corner points
     &        src_lons  ! longitudes of the four corner points

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::
     &        plat,   ! latitude  of the search point
     &        plon    ! longitude of the search point

      integer (kind=int_kind), dimension(2), intent(in) ::
     &        src_grid_dims  ! size of each src grid dimension

      real (kind=dbl_kind), dimension(:), intent(in) ::
     &        src_center_lat, ! latitude  of each src grid center 
     &        src_center_lon  ! longitude of each src grid center

      real (kind=dbl_kind), dimension(:,:), intent(in) ::
     &        src_bound_box   ! bounding box for src grid search

      integer (kind=int_kind), dimension(:,:), intent(in) ::
     &        src_bin_add,    ! search bins for restricting
     &        dst_bin_add     ! searches

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n, next_n, srch_add,   ! dummy indices
     &    nx, ny,            ! dimensions of src grid
     &    min_add, max_add,  ! addresses for restricting search
     &    i, j, jp1, ip1, n_add, e_add, ne_add  ! addresses

      real (kind=dbl_kind) ::  ! vectors for cross-product check
     &      vec1_lat, vec1_lon,
     &      vec2_lat, vec2_lon, cross_product, cross_product_last,
     &      coslat_dst, sinlat_dst, coslon_dst, sinlon_dst,
     &      dist_min, distance ! for computing dist-weighted avg

!-----------------------------------------------------------------------
!
!     restrict search first using search bins. 
!
!-----------------------------------------------------------------------

      src_add = 0

      min_add = size(src_center_lat)
      max_add = 1
      do n=1,num_srch_bins
        if (plat >= bin_lats(1,n) .and. plat <= bin_lats(2,n) .and.
     &      plon >= bin_lons(1,n) .and. plon <= bin_lons(2,n)) then
          min_add = min(min_add, src_bin_add(1,n))
          max_add = max(max_add, src_bin_add(2,n))
        endif
      end do
 
!-----------------------------------------------------------------------
!
!     now perform a more detailed search 
!
!-----------------------------------------------------------------------

      nx = src_grid_dims(1)
      ny = src_grid_dims(2)

      srch_loop: do srch_add = min_add,max_add

        if (plat <= src_bound_box(2,srch_add) .and. 
     &      plat >= src_bound_box(1,srch_add) .and.
     &      plon <= src_bound_box(4,srch_add) .and. 
     &      plon >= src_bound_box(3,srch_add)) then

          !***
          !*** we are within bounding box so get really serious
          !***

          !*** find N,S and NE points to this grid point

          j = (srch_add - 1)/nx +1
          i = srch_add - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
            ip1 = 1
          endif

          if (j < ny) then
            jp1 = j+1
          else
            jp1 = 1
          endif

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

          !***
          !*** find N,S and NE lat/lon coords and check bounding box
          !***

          src_lats(1) = src_center_lat(srch_add)
          src_lats(2) = src_center_lat(e_add)
          src_lats(3) = src_center_lat(ne_add)
          src_lats(4) = src_center_lat(n_add)

          src_lons(1) = src_center_lon(srch_add)
          src_lons(2) = src_center_lon(e_add)
          src_lons(3) = src_center_lon(ne_add)
          src_lons(4) = src_center_lon(n_add)

          !***
          !*** for consistency, we must make sure all lons are in
          !*** same 2pi interval
          !***

          vec1_lon = src_lons(1) - plon
          if (vec1_lon > pi) then
            src_lons(1) = src_lons(1) - pi2
          else if (vec1_lon < -pi) then
            src_lons(1) = src_lons(1) + pi2
          endif
          do n=2,4
            vec1_lon = src_lons(n) - src_lons(1)
            if (vec1_lon > pi) then
              src_lons(n) = src_lons(n) - pi2
            else if (vec1_lon < -pi) then
              src_lons(n) = src_lons(n) + pi2
            endif
          end do

          corner_loop: do n=1,4
            next_n = MOD(n,4) + 1

            !***
            !*** here we take the cross product of the vector making 
            !*** up each box side with the vector formed by the vertex
            !*** and search point.  if all the cross products are 
            !*** same sign, the point is contained in the box.
            !***

            vec1_lat = src_lats(next_n) - src_lats(n)
            vec1_lon = src_lons(next_n) - src_lons(n)
            vec2_lat = plat - src_lats(n)
            vec2_lon = plon - src_lons(n)

            !***
            !*** check for 0,2pi crossings
            !***

            if (vec1_lon >  three*pih) then
              vec1_lon = vec1_lon - pi2
            else if (vec1_lon < -three*pih) then
              vec1_lon = vec1_lon + pi2
            endif
            if (vec2_lon >  three*pih) then
              vec2_lon = vec2_lon - pi2
            else if (vec2_lon < -three*pih) then
              vec2_lon = vec2_lon + pi2
            endif

            cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***

            if (n==1) cross_product_last = cross_product
            if (cross_product*cross_product_last < zero) then
              exit corner_loop
            else
              cross_product_last = cross_product
            endif

          end do corner_loop

          !***
          !*** if cross products all positive, we found the location
          !***

          if (n > 4) then
            src_add(1) = srch_add
            src_add(2) = e_add
            src_add(3) = ne_add
            src_add(4) = n_add

            return
          endif

          !***
          !*** otherwise move on to next cell
          !***

        endif !bounding box check
      end do srch_loop

      !***
      !*** if no cell found, point is likely either in a box that
      !*** straddles either pole or is outside the grid.  fall back
      !*** to a distance-weighted average of the four closest
      !*** points.  go ahead and compute weights here, but store
      !*** in src_lats and return -add to prevent the parent
      !*** routine from computing bilinear weights
      !***

      coslat_dst = cos(plat)
      sinlat_dst = sin(plat)
      coslon_dst = cos(plon)
      sinlon_dst = sin(plon)

      dist_min = bignum
      src_lats = bignum
      do srch_add = min_add,max_add
        distance = acos(coslat_dst*cos(src_center_lat(srch_add))*
     &                 (coslon_dst*cos(src_center_lon(srch_add)) +
     &                  sinlon_dst*sin(src_center_lon(srch_add)))+
     &                  sinlat_dst*sin(src_center_lat(srch_add)))

        if (distance < dist_min) then
          sort_loop: do n=1,4
            if (distance < src_lats(n)) then
              do i=4,n+1,-1
                src_add (i) = src_add (i-1)
                src_lats(i) = src_lats(i-1)
              end do
              src_add (n) = -srch_add
              src_lats(n) = distance
              dist_min = src_lats(4)
              exit sort_loop
            endif
          end do sort_loop
        endif
      end do

      src_lons = one/(src_lats + tiny)
      distance = sum(src_lons)
      src_lats = src_lons/distance

!-----------------------------------------------------------------------

      end subroutine grid_search_bicub 

!***********************************************************************

      subroutine store_link_bicub(dst_add, src_add, weights, nmap)

!-----------------------------------------------------------------------
!
!     this routine stores the address and weight for four links 
!     associated with one destination point in the appropriate address 
!     and weight arrays and resizes those arrays if necessary.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(in) ::
     &        dst_add,  ! address on destination grid
     &        nmap      ! identifies which direction for mapping

      integer (kind=int_kind), dimension(4), intent(in) ::
     &        src_add   ! addresses on source grid

      real (kind=dbl_kind), dimension(4,4), intent(in) ::
     &        weights ! array of remapping weights for these links

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n, ! dummy index
     &       num_links_old          ! placeholder for old link number

!-----------------------------------------------------------------------
!
!     increment number of links and check to see if remap arrays need
!     to be increased to accomodate the new link.  then store the
!     link.
!
!-----------------------------------------------------------------------

      select case (nmap)
      case(1)

        num_links_old  = num_links_map1
        num_links_map1 = num_links_old + 4

        if (num_links_map1 > max_links_map1) 
     &     call resize_remap_vars(1,resize_increment)

        do n=1,4
          grid1_add_map1(num_links_old+n) = src_add(n)
          grid2_add_map1(num_links_old+n) = dst_add
          wts_map1    (:,num_links_old+n) = weights(:,n)
        end do

      case(2)

        num_links_old  = num_links_map2
        num_links_map2 = num_links_old + 4

        if (num_links_map2 > max_links_map2) 
     &     call resize_remap_vars(2,resize_increment)

        do n=1,4
          grid1_add_map2(num_links_old+n) = dst_add
          grid2_add_map2(num_links_old+n) = src_add(n)
          wts_map2    (:,num_links_old+n) = weights(:,n)
        end do

      end select

!-----------------------------------------------------------------------

      end subroutine store_link_bicub

!***********************************************************************

      end module remap_bicubic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
