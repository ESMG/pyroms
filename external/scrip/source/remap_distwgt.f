!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains necessary routines for performing an 
!     interpolation using a distance-weighted average of n nearest
!     neighbors.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_distwgt.f,v 1.3 2000/04/19 21:56:26 pwjones Exp $
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

      module remap_distance_weight

!-----------------------------------------------------------------------

      use kinds_mod     ! defines common data types
      use constants     ! defines common constants
      use grids         ! module containing grid info
      use remap_vars    ! module containing remap info

      implicit none

!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), parameter :: 
     &     num_neighbors=4  ! num nearest neighbors to interpolate from

      real (kind=dbl_kind), dimension(:), allocatable, save ::
     &     coslat, sinlat, ! cosine, sine of grid lats (for distance)
     &     coslon, sinlon, ! cosine, sine of grid lons (for distance)
     &     wgtstmp         ! an array to hold the link weight

!***********************************************************************

      contains

!***********************************************************************

      subroutine remap_distwgt

!-----------------------------------------------------------------------
!
!     this routine computes the inverse-distance weights for a
!     nearest-neighbor interpolation.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      logical (kind=log_kind), dimension(num_neighbors) ::
     &     nbr_mask        ! mask at nearest neighbors

      integer (kind=int_kind) :: n,
     &     dst_add,        ! destination address
     &     nmap            ! index of current map being computed

      integer (kind=int_kind), dimension(num_neighbors) ::
     &     nbr_add         ! source address at nearest neighbors

      real (kind=dbl_kind), dimension(num_neighbors) ::
     &     nbr_dist        ! angular distance four nearest neighbors

      real (kind=dbl_kind) ::
     &     coslat_dst,     ! cos(lat) of destination grid point
     &     coslon_dst,     ! cos(lon) of destination grid point
     &     sinlat_dst,     ! sin(lat) of destination grid point
     &     sinlon_dst,     ! sin(lon) of destination grid point
     &     dist_tot        ! sum of neighbor distances (for normalizing)

!-----------------------------------------------------------------------
!
!     compute mappings from grid1 to grid2
!
!-----------------------------------------------------------------------

      nmap = 1

      !***
      !*** allocate wgtstmp to be consistent with store_link interface
      !***

      allocate (wgtstmp(num_wts))

      !***
      !*** compute cos, sin of lat/lon on source grid for distance
      !*** calculations
      !***

      allocate (coslat(grid1_size), coslon(grid1_size),
     &          sinlat(grid1_size), sinlon(grid1_size))

      coslat = cos(grid1_center_lat)
      coslon = cos(grid1_center_lon)
      sinlat = sin(grid1_center_lat)
      sinlon = sin(grid1_center_lon)

      !***
      !*** loop over destination grid 
      !***

      grid_loop1: do dst_add = 1, grid2_size

        if (.not. grid2_mask(dst_add)) cycle grid_loop1

        coslat_dst = cos(grid2_center_lat(dst_add))
        coslon_dst = cos(grid2_center_lon(dst_add))
        sinlat_dst = sin(grid2_center_lat(dst_add))
        sinlon_dst = sin(grid2_center_lon(dst_add))

        !***
        !*** find nearest grid points on source grid and
        !*** distances to each point
        !***

        call grid_search_nbr(nbr_add, nbr_dist, 
     &                       grid2_center_lat(dst_add),
     &                       grid2_center_lon(dst_add),
     &                       coslat_dst, coslon_dst, 
     &                       sinlat_dst, sinlon_dst,
     &                       bin_addr1, bin_addr2)

        !***
        !*** compute weights based on inverse distance
        !*** if mask is false, eliminate those points
        !***

        dist_tot = zero
        do n=1,num_neighbors
          if (grid1_mask(nbr_add(n))) then
            nbr_dist(n) = one/nbr_dist(n)
            dist_tot = dist_tot + nbr_dist(n)
            nbr_mask(n) = .true.
          else
            nbr_mask(n) = .false.
          endif
        end do

        !***
        !*** normalize weights and store the link
        !***

        do n=1,num_neighbors
          if (nbr_mask(n)) then
            wgtstmp(1) = nbr_dist(n)/dist_tot
            call store_link_nbr(nbr_add(n), dst_add, wgtstmp, nmap)
            grid2_frac(dst_add) = one
          endif
        end do

      end do grid_loop1

      deallocate (coslat, coslon, sinlat, sinlon)

!-----------------------------------------------------------------------
!
!     compute mappings from grid2 to grid1 if necessary
!
!-----------------------------------------------------------------------

      if (num_maps > 1) then

      nmap = 2

      !***
      !*** compute cos, sin of lat/lon on source grid for distance
      !*** calculations
      !***

      allocate (coslat(grid2_size), coslon(grid2_size),
     &          sinlat(grid2_size), sinlon(grid2_size))

      coslat = cos(grid2_center_lat)
      coslon = cos(grid2_center_lon)
      sinlat = sin(grid2_center_lat)
      sinlon = sin(grid2_center_lon)

      !***
      !*** loop over destination grid 
      !***

      grid_loop2: do dst_add = 1, grid1_size

        if (.not. grid1_mask(dst_add)) cycle grid_loop2

        coslat_dst = cos(grid1_center_lat(dst_add))
        coslon_dst = cos(grid1_center_lon(dst_add))
        sinlat_dst = sin(grid1_center_lat(dst_add))
        sinlon_dst = sin(grid1_center_lon(dst_add))

        !***
        !*** find four nearest grid points on source grid and
        !*** distances to each point
        !***

        call grid_search_nbr(nbr_add, nbr_dist,
     &                       grid1_center_lat(dst_add),
     &                       grid1_center_lon(dst_add),
     &                       coslat_dst, coslon_dst, 
     &                       sinlat_dst, sinlon_dst,
     &                       bin_addr2, bin_addr1)

        !***
        !*** compute weights based on inverse distance
        !*** if mask is false, eliminate those points
        !***

        dist_tot = zero
        do n=1,num_neighbors
          if (grid2_mask(nbr_add(n))) then
            nbr_dist(n) = one/nbr_dist(n)
            dist_tot = dist_tot + nbr_dist(n)
            nbr_mask(n) = .true.
          else
            nbr_mask(n) = .false.
          endif
        end do

        !***
        !*** normalize weights and store the link
        !***

        do n=1,num_neighbors
          if (nbr_mask(n)) then
            wgtstmp(1) = nbr_dist(n)/dist_tot
            call store_link_nbr(dst_add, nbr_add(n), wgtstmp, nmap)
            grid1_frac(dst_add) = one
          endif
        end do

      end do grid_loop2

      deallocate (coslat, coslon, sinlat, sinlon)

      endif

      deallocate(wgtstmp)

!-----------------------------------------------------------------------

      end subroutine remap_distwgt

!***********************************************************************

      subroutine grid_search_nbr(nbr_add, nbr_dist, plat, plon, 
     &               coslat_dst, coslon_dst, sinlat_dst, sinlon_dst,
     &               src_bin_add, dst_bin_add, check_mask)

!-----------------------------------------------------------------------
!
!     this routine finds the closest num_neighbor points to a search 
!     point and computes a distance to each of the neighbors.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), dimension(num_neighbors), intent(out) ::
     &        nbr_add  ! address of each of the closest points

      real (kind=dbl_kind), dimension(num_neighbors), intent(out) ::
     &        nbr_dist ! distance to each of the closest points

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), dimension(:,:), intent(in) ::
     &        src_bin_add, ! search bins for restricting search
     &        dst_bin_add   

      logical (kind=log_kind), intent(in), optional ::
     &        check_mask

      real (kind=dbl_kind), intent(in) ::
     &        plat,         ! latitude  of the search point
     &        plon,         ! longitude of the search point
     &        coslat_dst,   ! cos(lat)  of the search point
     &        coslon_dst,   ! cos(lon)  of the search point
     &        sinlat_dst,   ! sin(lat)  of the search point
     &        sinlon_dst    ! sin(lon)  of the search point

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n, nmax, nadd, nchk, ! dummy indices
     &        min_add, max_add, nm1, np1, i, j, ip1, im1, jp1, jm1

      real (kind=dbl_kind) ::
     &        distance      ! angular distance

!-----------------------------------------------------------------------
!
!     loop over source grid and find nearest neighbors
!
!-----------------------------------------------------------------------

      !***
      !*** restrict the search using search bins
      !*** expand the bins to catch neighbors
      !***

      select case (restrict_type)
      case('latitude')

        do n=1,num_srch_bins
          if (plat >= bin_lats(1,n) .and. plat <= bin_lats(2,n)) then
            min_add = src_bin_add(1,n)
            max_add = src_bin_add(2,n)

            nm1 = max(n-1,1)
            np1 = min(n+1,num_srch_bins)

            min_add = min(min_add,src_bin_add(1,nm1))
            max_add = max(max_add,src_bin_add(2,nm1))
            min_add = min(min_add,src_bin_add(1,np1))
            max_add = max(max_add,src_bin_add(2,np1))
          endif
        end do

      case('latlon')

        n = 0
        nmax = nint(sqrt(real(num_srch_bins)))
        do j=1,nmax
        jp1 = min(j+1,nmax)
        jm1 = max(j-1,1)
        do i=1,nmax
          ip1 = min(i+1,nmax)
          im1 = max(i-1,1)

          n = n+1
          if (plat >= bin_lats(1,n) .and. plat <= bin_lats(2,n) .and.
     &        plon >= bin_lons(1,n) .and. plon <= bin_lons(3,n)) then
            min_add = src_bin_add(1,n)
            max_add = src_bin_add(2,n)

            nm1 = (jm1-1)*nmax + im1
            np1 = (jp1-1)*nmax + ip1
            nm1 = max(nm1,1)
            np1 = min(np1,num_srch_bins)

            min_add = min(min_add,src_bin_add(1,nm1))
            max_add = max(max_add,src_bin_add(2,nm1))
            min_add = min(min_add,src_bin_add(1,np1))
            max_add = max(max_add,src_bin_add(2,np1))
          endif
        end do
        end do

      end select

      !***
      !*** initialize distance and address arrays
      !***

      nbr_add = 0
      nbr_dist = bignum

      do nadd=min_add,max_add

        !***
        !*** find distance to this point
        !***
!       jw adding this mask test to return only valid unmasked
!       neighbors as an option, used right now only in the 
!       modified bilinear remap
        if (present(check_mask)) then
          if (check_mask .AND. (.not. grid1_mask(nadd))) cycle
        endif

        distance = acos(sinlat_dst*sinlat(nadd) +
     &                  coslat_dst*coslat(nadd)*
     &                 (coslon_dst*coslon(nadd) +
     &                  sinlon_dst*sinlon(nadd)) )

        !***
        !*** store the address and distance if this is one of the
        !*** smallest four so far
        !***

        check_loop: do nchk=1,num_neighbors
          if (distance .lt. nbr_dist(nchk)) then
            do n=num_neighbors,nchk+1,-1
              nbr_add(n) = nbr_add(n-1)
              nbr_dist(n) = nbr_dist(n-1)
            end do
            nbr_add(nchk) = nadd
            nbr_dist(nchk) = distance
            exit check_loop
          endif
        end do check_loop

      end do

!-----------------------------------------------------------------------

      end subroutine grid_search_nbr 

!***********************************************************************

      subroutine store_link_nbr(add1, add2, weights, nmap)

!-----------------------------------------------------------------------
!
!     this routine stores the address and weight for this link in
!     the appropriate address and weight arrays and resizes those
!     arrays if necessary.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(in) ::
     &        add1,  ! address on grid1
     &        add2,  ! address on grid2
     &        nmap   ! identifies which direction for mapping

      real (kind=dbl_kind), dimension(:), intent(in) ::
     &        weights ! array of remapping weights for this link

!-----------------------------------------------------------------------
!
!     increment number of links and check to see if remap arrays need
!     to be increased to accomodate the new link.  then store the
!     link.
!
!-----------------------------------------------------------------------

      select case (nmap)
      case(1)

        num_links_map1  = num_links_map1 + 1

        if (num_links_map1 > max_links_map1) 
     &     call resize_remap_vars(1,resize_increment)

        grid1_add_map1(num_links_map1) = add1
        grid2_add_map1(num_links_map1) = add2
        wts_map1    (:,num_links_map1) = weights

      case(2)

        num_links_map2  = num_links_map2 + 1

        if (num_links_map2 > max_links_map2) 
     &     call resize_remap_vars(2,resize_increment)

        grid1_add_map2(num_links_map2) = add1
        grid2_add_map2(num_links_map2) = add2
        wts_map2    (:,num_links_map2) = weights

      end select

!-----------------------------------------------------------------------

      end subroutine store_link_nbr

!***********************************************************************

      end module remap_distance_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
