!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains necessary routines for computing addresses
!     and weights for a conservative interpolation  between any two 
!     grids on a sphere.  the weights are computed by performing line 
!     integrals around all overlap regions of the two grids.  see 
!     Dukowicz and Kodis, SIAM J. Sci. Stat. Comput. 8, 305 (1987) and
!     Jones, P.W. Monthly Weather Review (submitted).
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_conserv.f,v 1.10 2001/08/21 21:05:13 pwjones Exp $
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

      module remap_conservative

!-----------------------------------------------------------------------

      use kinds_mod    ! defines common data types
      use constants    ! defines common constants
      use timers       ! module for timing
      use grids        ! module containing grid information
      use remap_vars   ! module containing remap information

      implicit none

!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), save :: 
     &        num_srch_cells ! num cells in restricted search arrays

      integer (kind=int_kind), dimension(:), allocatable, save :: 
     &        srch_add       ! global address of cells in srch arrays

      real (kind=dbl_kind), parameter :: 
     &     north_thresh = 1.45_dbl_kind, ! threshold for coord transf.
     &     south_thresh =-2.00_dbl_kind  ! threshold for coord transf.

      real (kind=dbl_kind), dimension(:,:), allocatable, save ::
     &     srch_corner_lat,  ! lat of each corner of srch cells
     &     srch_corner_lon   ! lon of each corner of srch cells

!***********************************************************************

      contains

!***********************************************************************

      subroutine remap_conserv

!-----------------------------------------------------------------------
!
!     this routine traces the perimeters of every grid cell on each
!     grid checking for intersections with the other grid and computing
!     line integrals for each subsegment.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), parameter :: 
     &        max_subseg = 10000 ! max number of subsegments per segment
                                 ! to prevent infinite loop

      integer (kind=int_kind) :: 
     &        grid1_add,  ! current linear address for grid1 cell
     &        grid2_add,  ! current linear address for grid2 cell
     &        min_add,    ! addresses for restricting search of
     &        max_add,    !   destination grid
     &        n, nwgt,    ! generic counters
     &        corner,     ! corner of cell that segment starts from
     &        next_corn,  ! corner of cell that segment ends on
     &        num_subseg  ! number of subsegments 

      logical (kind=log_kind) :: 
     &        lcoinc,  ! flag for coincident segments
     &        lrevers, ! flag for reversing direction of segment
     &        lbegin   ! flag for first integration of a segment

      logical (kind=log_kind), dimension(:), allocatable ::
     &        srch_mask   ! mask for restricting searches

      real (kind=dbl_kind) ::
     &     intrsct_lat, intrsct_lon,       ! lat/lon of next intersect
     &     beglat, endlat, beglon, endlon, ! endpoints of current seg.
     &     norm_factor                     ! factor for normalizing wts

      real (kind=dbl_kind), dimension(:), allocatable ::
     &       grid2_centroid_lat, grid2_centroid_lon, ! centroid coords
     &       grid1_centroid_lat, grid1_centroid_lon  ! on each grid

      real (kind=dbl_kind), dimension(2) :: begseg ! begin lat/lon for
                                                   ! full segment

      real (kind=dbl_kind), dimension(6) :: weights ! local wgt array

!-----------------------------------------------------------------------
!
!     initialize centroid arrays
!
!-----------------------------------------------------------------------

      allocate( grid1_centroid_lat(grid1_size),
     &          grid1_centroid_lon(grid1_size),
     &          grid2_centroid_lat(grid2_size),
     &          grid2_centroid_lon(grid2_size))

      grid1_centroid_lat = zero
      grid1_centroid_lon = zero
      grid2_centroid_lat = zero
      grid2_centroid_lon = zero

!-----------------------------------------------------------------------
!
!     integrate around each cell on grid1
!
!-----------------------------------------------------------------------

      allocate(srch_mask(grid2_size))

      print *,'grid1 sweep '
      do grid1_add = 1,grid1_size

        !***
        !*** restrict searches first using search bins
        !***

        call timer_start(1)
        min_add = grid2_size
        max_add = 1
        do n=1,num_srch_bins
          if (grid1_add >= bin_addr1(1,n) .and.
     &        grid1_add <= bin_addr1(2,n)) then
            min_add = min(min_add, bin_addr2(1,n))
            max_add = max(max_add, bin_addr2(2,n))
          endif
        end do

        !***
        !*** further restrict searches using bounding boxes
        !***

        num_srch_cells = 0
        do grid2_add = min_add,max_add
          srch_mask(grid2_add) = (grid2_bound_box(1,grid2_add) <= 
     &                            grid1_bound_box(2,grid1_add)) .and.
     &                           (grid2_bound_box(2,grid2_add) >= 
     &                            grid1_bound_box(1,grid1_add)) .and.
     &                           (grid2_bound_box(3,grid2_add) <= 
     &                            grid1_bound_box(4,grid1_add)) .and.
     &                           (grid2_bound_box(4,grid2_add) >= 
     &                            grid1_bound_box(3,grid1_add))

          if (srch_mask(grid2_add)) num_srch_cells = num_srch_cells+1
        end do

        !***
        !*** create search arrays
        !***

        allocate(srch_add(num_srch_cells),
     &           srch_corner_lat(grid2_corners,num_srch_cells),
     &           srch_corner_lon(grid2_corners,num_srch_cells))

        n = 0
        gather1: do grid2_add = min_add,max_add
          if (srch_mask(grid2_add)) then
            n = n+1
            srch_add(n) = grid2_add
            srch_corner_lat(:,n) = grid2_corner_lat(:,grid2_add)
            srch_corner_lon(:,n) = grid2_corner_lon(:,grid2_add)
          endif
        end do gather1
        call timer_stop(1)

        !***
        !*** integrate around this cell
        !***

        do corner = 1,grid1_corners
          next_corn = mod(corner,grid1_corners) + 1

          !***
          !*** define endpoints of the current segment
          !***

          beglat = grid1_corner_lat(corner,grid1_add)
          beglon = grid1_corner_lon(corner,grid1_add)
          endlat = grid1_corner_lat(next_corn,grid1_add)
          endlon = grid1_corner_lon(next_corn,grid1_add)
          lrevers = .false.

          !***
          !*** to ensure exact path taken during both
          !*** sweeps, always integrate segments in the same 
          !*** direction (SW to NE).
          !***

          if ((endlat < beglat) .or.
     &        (endlat == beglat .and. endlon < beglon)) then 
            beglat = grid1_corner_lat(next_corn,grid1_add)
            beglon = grid1_corner_lon(next_corn,grid1_add)
            endlat = grid1_corner_lat(corner,grid1_add)
            endlon = grid1_corner_lon(corner,grid1_add)
            lrevers = .true.
          endif

          begseg(1) = beglat
          begseg(2) = beglon
          lbegin = .true.
          num_subseg = 0

          !***
          !*** if this is a constant-longitude segment, skip the rest 
          !*** since the line integral contribution will be zero.
          !***

          if (endlon /= beglon) then

          !***
          !*** integrate along this segment, detecting intersections 
          !*** and computing the line integral for each sub-segment
          !***

          do while (beglat /= endlat .or. beglon /= endlon)

            !***
            !*** prevent infinite loops if integration gets stuck
            !*** near cell or threshold boundary
            !***

            num_subseg = num_subseg + 1
            if (num_subseg > max_subseg) then
              stop 'integration stalled: num_subseg exceeded limit'
            endif

            !***
            !*** find next intersection of this segment with a grid
            !*** line on grid 2.
            !***

            call timer_start(2)
            call intersection(grid2_add,intrsct_lat,intrsct_lon,lcoinc,
     &                        beglat, beglon, endlat, endlon, begseg, 
     &                        lbegin, lrevers)
            call timer_stop(2)
            lbegin = .false.

            !***
            !*** compute line integral for this subsegment.
            !***

            call timer_start(3)
            if (grid2_add /= 0) then
              call line_integral(weights, num_wts,
     &                         beglon, intrsct_lon, beglat, intrsct_lat,
     &                         grid1_center_lat(grid1_add), 
     &                         grid1_center_lon(grid1_add),
     &                         grid2_center_lat(grid2_add), 
     &                         grid2_center_lon(grid2_add))
            else
              call line_integral(weights, num_wts,
     &                         beglon, intrsct_lon, beglat, intrsct_lat,
     &                         grid1_center_lat(grid1_add), 
     &                         grid1_center_lon(grid1_add),
     &                         grid1_center_lat(grid1_add), 
     &                         grid1_center_lon(grid1_add))
            endif
            call timer_stop(3)

            !***
            !*** if integrating in reverse order, change
            !*** sign of weights
            !***

            if (lrevers) then
              weights = -weights
            endif

            !***
            !*** store the appropriate addresses and weights. 
            !*** also add contributions to cell areas and centroids.
            !***

            !if (grid1_add == 119247) then
            !  print *,grid1_add,grid2_add,corner,weights(1)
            !  print *,grid1_corner_lat(:,grid1_add)
            !  print *,grid1_corner_lon(:,grid1_add)
            !  print *,grid2_corner_lat(:,grid2_add)
            !  print *,grid2_corner_lon(:,grid2_add)
            !  print *,beglat,beglon,intrsct_lat,intrsct_lon
            !endif

            if (grid2_add /= 0) then
              if (grid1_mask(grid1_add)) then
                call timer_start(4)
                call store_link_cnsrv(grid1_add, grid2_add, weights)
                call timer_stop(4)
                grid1_frac(grid1_add) = grid1_frac(grid1_add) + 
     &                                  weights(1)
                grid2_frac(grid2_add) = grid2_frac(grid2_add) + 
     &                                  weights(num_wts+1)
              endif

            endif

            grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1)
            grid1_centroid_lat(grid1_add) = 
     &      grid1_centroid_lat(grid1_add) + weights(2)
            grid1_centroid_lon(grid1_add) = 
     &      grid1_centroid_lon(grid1_add) + weights(3)

            !***
            !*** reset beglat and beglon for next subsegment.
            !***

            beglat = intrsct_lat
            beglon = intrsct_lon
          end do

          endif

          !***
          !*** end of segment
          !***

        end do

        !***
        !*** finished with this cell: deallocate search array and
        !*** start on next cell

        deallocate(srch_add, srch_corner_lat, srch_corner_lon)

      end do

      deallocate(srch_mask)

!-----------------------------------------------------------------------
!
!     integrate around each cell on grid2
!
!-----------------------------------------------------------------------

      allocate(srch_mask(grid1_size))

      print *,'grid2 sweep '
      do grid2_add = 1,grid2_size

        !***
        !*** restrict searches first using search bins
        !***

        call timer_start(5)
        min_add = grid1_size
        max_add = 1
        do n=1,num_srch_bins
          if (grid2_add >= bin_addr2(1,n) .and.
     &        grid2_add <= bin_addr2(2,n)) then
            min_add = min(min_add, bin_addr1(1,n))
            max_add = max(max_add, bin_addr1(2,n))
          endif
        end do

        !***
        !*** further restrict searches using bounding boxes
        !***

        num_srch_cells = 0
        do grid1_add = min_add, max_add
          srch_mask(grid1_add) = (grid1_bound_box(1,grid1_add) <= 
     &                            grid2_bound_box(2,grid2_add)) .and.
     &                           (grid1_bound_box(2,grid1_add) >= 
     &                            grid2_bound_box(1,grid2_add)) .and.
     &                           (grid1_bound_box(3,grid1_add) <= 
     &                            grid2_bound_box(4,grid2_add)) .and.
     &                           (grid1_bound_box(4,grid1_add) >= 
     &                            grid2_bound_box(3,grid2_add))

          if (srch_mask(grid1_add)) num_srch_cells = num_srch_cells+1
        end do

        allocate(srch_add(num_srch_cells),
     &           srch_corner_lat(grid1_corners,num_srch_cells),
     &           srch_corner_lon(grid1_corners,num_srch_cells))

        n = 0
        gather2: do grid1_add = min_add,max_add
          if (srch_mask(grid1_add)) then
            n = n+1
            srch_add(n) = grid1_add
            srch_corner_lat(:,n) = grid1_corner_lat(:,grid1_add)
            srch_corner_lon(:,n) = grid1_corner_lon(:,grid1_add)
          endif
        end do gather2
        call timer_stop(5)

        !***
        !*** integrate around this cell
        !***

        do corner = 1,grid2_corners
          next_corn = mod(corner,grid2_corners) + 1

          beglat = grid2_corner_lat(corner,grid2_add)
          beglon = grid2_corner_lon(corner,grid2_add)
          endlat = grid2_corner_lat(next_corn,grid2_add)
          endlon = grid2_corner_lon(next_corn,grid2_add)
          lrevers = .false.

          !***
          !*** to ensure exact path taken during both
          !*** sweeps, always integrate in the same direction
          !***

          if ((endlat < beglat) .or.
     &        (endlat == beglat .and. endlon < beglon)) then 
            beglat = grid2_corner_lat(next_corn,grid2_add)
            beglon = grid2_corner_lon(next_corn,grid2_add)
            endlat = grid2_corner_lat(corner,grid2_add)
            endlon = grid2_corner_lon(corner,grid2_add)
            lrevers = .true.
          endif

          begseg(1) = beglat
          begseg(2) = beglon
          lbegin = .true.

          !***
          !*** if this is a constant-longitude segment, skip the rest 
          !*** since the line integral contribution will be zero.
          !***

          if (endlon /= beglon) then
          num_subseg = 0

          !***
          !*** integrate along this segment, detecting intersections 
          !*** and computing the line integral for each sub-segment
          !***

          do while (beglat /= endlat .or. beglon /= endlon)

            !***
            !*** prevent infinite loops if integration gets stuck
            !*** near cell or threshold boundary
            !***

            num_subseg = num_subseg + 1
            if (num_subseg > max_subseg) then
              stop 'integration stalled: num_subseg exceeded limit'
            endif

            !***
            !*** find next intersection of this segment with a line 
            !*** on grid 2.
            !***

            call timer_start(6)
            call intersection(grid1_add,intrsct_lat,intrsct_lon,lcoinc,
     &                        beglat, beglon, endlat, endlon, begseg,
     &                        lbegin, lrevers)
            call timer_stop(6)
            lbegin = .false.

            !***
            !*** compute line integral for this subsegment.
            !***

            call timer_start(7)
            if (grid1_add /= 0) then
              call line_integral(weights, num_wts,
     &                         beglon, intrsct_lon, beglat, intrsct_lat,
     &                         grid1_center_lat(grid1_add), 
     &                         grid1_center_lon(grid1_add),
     &                         grid2_center_lat(grid2_add), 
     &                         grid2_center_lon(grid2_add))
            else
              call line_integral(weights, num_wts,
     &                         beglon, intrsct_lon, beglat, intrsct_lat,
     &                         grid2_center_lat(grid2_add), 
     &                         grid2_center_lon(grid2_add),
     &                         grid2_center_lat(grid2_add), 
     &                         grid2_center_lon(grid2_add))
            endif
            call timer_stop(7)

            if (lrevers) then
              weights = -weights
            endif

            !***
            !*** store the appropriate addresses and weights. 
            !*** also add contributions to cell areas and centroids.
            !*** if there is a coincidence, do not store weights
            !*** because they have been captured in the previous loop.
            !*** the grid1 mask is the master mask
            !***

            !if (grid1_add == 119247) then
            !  print *,grid1_add,grid2_add,corner,weights(1)
            !  print *,grid1_corner_lat(:,grid1_add)
            !  print *,grid1_corner_lon(:,grid1_add)
            !  print *,grid2_corner_lat(:,grid2_add)
            !  print *,grid2_corner_lon(:,grid2_add)
            !  print *,beglat,beglon,intrsct_lat,intrsct_lon
            !endif

            if (.not. lcoinc .and. grid1_add /= 0) then
              if (grid1_mask(grid1_add)) then
                call timer_start(8)
                call store_link_cnsrv(grid1_add, grid2_add, weights)
                call timer_stop(8)
                grid1_frac(grid1_add) = grid1_frac(grid1_add) + 
     &                                  weights(1)
                grid2_frac(grid2_add) = grid2_frac(grid2_add) + 
     &                                  weights(num_wts+1)
              endif

            endif

            grid2_area(grid2_add) = grid2_area(grid2_add) + 
     &                                      weights(num_wts+1)
            grid2_centroid_lat(grid2_add) = 
     &      grid2_centroid_lat(grid2_add) + weights(num_wts+2)
            grid2_centroid_lon(grid2_add) = 
     &      grid2_centroid_lon(grid2_add) + weights(num_wts+3)

            !***
            !*** reset beglat and beglon for next subsegment.
            !***

            beglat = intrsct_lat
            beglon = intrsct_lon
          end do

          endif

          !***
          !*** end of segment
          !***

        end do

        !***
        !*** finished with this cell: deallocate search array and
        !*** start on next cell

        deallocate(srch_add, srch_corner_lat, srch_corner_lon)

      end do

      deallocate(srch_mask)

!-----------------------------------------------------------------------
!
!     correct for situations where N/S pole not explicitly included in
!     grid (i.e. as a grid corner point). if pole is missing from only
!     one grid, need to correct only the area and centroid of that 
!     grid.  if missing from both, do complete weight calculation.
!
!-----------------------------------------------------------------------

      !*** North Pole
      weights(1) =  pi2
      weights(2) =  pi*pi
      weights(3) =  zero
      weights(4) =  pi2
      weights(5) =  pi*pi
      weights(6) =  zero

      grid1_add = 0
      pole_loop1: do n=1,grid1_size
        if (grid1_area(n) < -three*pih .and.
     &      grid1_center_lat(n) > zero) then
          grid1_add = n
          exit pole_loop1
        endif
      end do pole_loop1

      grid2_add = 0
      pole_loop2: do n=1,grid2_size
        if (grid2_area(n) < -three*pih .and.
     &      grid2_center_lat(n) > zero) then
          grid2_add = n
          exit pole_loop2
        endif
      end do pole_loop2

      if (grid1_add /=0) then
        grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1)
        grid1_centroid_lat(grid1_add) = 
     &  grid1_centroid_lat(grid1_add) + weights(2)
        grid1_centroid_lon(grid1_add) =
     &  grid1_centroid_lon(grid1_add) + weights(3)
      endif

      if (grid2_add /=0) then
        grid2_area(grid2_add) = grid2_area(grid2_add) + 
     &                                  weights(num_wts+1)
        grid2_centroid_lat(grid2_add) = 
     &  grid2_centroid_lat(grid2_add) + weights(num_wts+2)
        grid2_centroid_lon(grid2_add) =
     &  grid2_centroid_lon(grid2_add) + weights(num_wts+3)
      endif

      if (grid1_add /= 0 .and. grid2_add /=0) then
        call store_link_cnsrv(grid1_add, grid2_add, weights)

        grid1_frac(grid1_add) = grid1_frac(grid1_add) + 
     &                          weights(1)
        grid2_frac(grid2_add) = grid2_frac(grid2_add) + 
     &                          weights(num_wts+1)
      endif

      !*** South Pole
      weights(1) =  pi2
      weights(2) = -pi*pi
      weights(3) =  zero
      weights(4) =  pi2
      weights(5) = -pi*pi
      weights(6) =  zero

      grid1_add = 0
      pole_loop3: do n=1,grid1_size
        if (grid1_area(n) < -three*pih .and.
     &      grid1_center_lat(n) < zero) then
          grid1_add = n
          exit pole_loop3
        endif
      end do pole_loop3

      grid2_add = 0
      pole_loop4: do n=1,grid2_size
        if (grid2_area(n) < -three*pih .and.
     &      grid2_center_lat(n) < zero) then
          grid2_add = n
          exit pole_loop4
        endif
      end do pole_loop4

      if (grid1_add /=0) then
        grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1)
        grid1_centroid_lat(grid1_add) = 
     &  grid1_centroid_lat(grid1_add) + weights(2)
        grid1_centroid_lon(grid1_add) =
     &  grid1_centroid_lon(grid1_add) + weights(3)
      endif

      if (grid2_add /=0) then
        grid2_area(grid2_add) = grid2_area(grid2_add) + 
     &                                  weights(num_wts+1)
        grid2_centroid_lat(grid2_add) = 
     &  grid2_centroid_lat(grid2_add) + weights(num_wts+2)
        grid2_centroid_lon(grid2_add) =
     &  grid2_centroid_lon(grid2_add) + weights(num_wts+3)
      endif

      if (grid1_add /= 0 .and. grid2_add /=0) then
        call store_link_cnsrv(grid1_add, grid2_add, weights)

        grid1_frac(grid1_add) = grid1_frac(grid1_add) + 
     &                          weights(1)
        grid2_frac(grid2_add) = grid2_frac(grid2_add) + 
     &                          weights(num_wts+1)
      endif

!-----------------------------------------------------------------------
!
!     finish centroid computation
!
!-----------------------------------------------------------------------

      where (grid1_area /= zero)
        grid1_centroid_lat = grid1_centroid_lat/grid1_area
        grid1_centroid_lon = grid1_centroid_lon/grid1_area
      end where

      where (grid2_area /= zero)
        grid2_centroid_lat = grid2_centroid_lat/grid2_area
        grid2_centroid_lon = grid2_centroid_lon/grid2_area
      end where

!-----------------------------------------------------------------------
!
!     include centroids in weights and normalize using destination
!     area if requested
!
!-----------------------------------------------------------------------

      do n=1,num_links_map1
        grid1_add = grid1_add_map1(n)
        grid2_add = grid2_add_map1(n)
        do nwgt=1,num_wts
          weights(        nwgt) = wts_map1(nwgt,n)
          if (num_maps > 1) then
            weights(num_wts+nwgt) = wts_map2(nwgt,n)
          endif
        end do

        select case(norm_opt)
        case (norm_opt_dstarea)
          if (grid2_area(grid2_add) /= zero) then
            if (luse_grid2_area) then
              norm_factor = one/grid2_area_in(grid2_add)
            else
              norm_factor = one/grid2_area(grid2_add)
            endif
          else
            norm_factor = zero
          endif
        case (norm_opt_frcarea)
          if (grid2_frac(grid2_add) /= zero) then
            if (luse_grid2_area) then
              norm_factor = grid2_area(grid2_add)/
     &                     (grid2_frac(grid2_add)*
     &                      grid2_area_in(grid2_add))
            else
              norm_factor = one/grid2_frac(grid2_add)
            endif
          else
            norm_factor = zero
          endif
        case (norm_opt_none)
          norm_factor = one
        end select

        wts_map1(1,n) =  weights(1)*norm_factor
        wts_map1(2,n) = (weights(2) - weights(1)*
     &                              grid1_centroid_lat(grid1_add))*
     &                              norm_factor
        wts_map1(3,n) = (weights(3) - weights(1)*
     &                              grid1_centroid_lon(grid1_add))*
     &                              norm_factor

        if (num_maps > 1) then
          select case(norm_opt)
          case (norm_opt_dstarea)
            if (grid1_area(grid1_add) /= zero) then
              if (luse_grid1_area) then
                norm_factor = one/grid1_area_in(grid1_add)
              else
                norm_factor = one/grid1_area(grid1_add)
              endif
            else
              norm_factor = zero
            endif
          case (norm_opt_frcarea)
            if (grid1_frac(grid1_add) /= zero) then
              if (luse_grid1_area) then
                norm_factor = grid1_area(grid1_add)/
     &                       (grid1_frac(grid1_add)*
     &                        grid1_area_in(grid1_add))
              else
                norm_factor = one/grid1_frac(grid1_add)
              endif
            else
              norm_factor = zero
            endif
          case (norm_opt_none)
            norm_factor = one
          end select

          wts_map2(1,n) =  weights(num_wts+1)*norm_factor
          wts_map2(2,n) = (weights(num_wts+2) - weights(num_wts+1)*
     &                                grid2_centroid_lat(grid2_add))*
     &                                norm_factor
          wts_map2(3,n) = (weights(num_wts+3) - weights(num_wts+1)*
     &                                grid2_centroid_lon(grid2_add))*
     &                                norm_factor
        endif

      end do

      print *, 'Total number of links = ',num_links_map1

      where (grid1_area /= zero) grid1_frac = grid1_frac/grid1_area
      where (grid2_area /= zero) grid2_frac = grid2_frac/grid2_area

!-----------------------------------------------------------------------
!
!     perform some error checking on final weights
!
!-----------------------------------------------------------------------

      grid2_centroid_lat = zero
      grid2_centroid_lon = zero

      do n=1,grid1_size
        if (grid1_area(n) < -.01) then
          print *,'Grid 1 area error: ',n,grid1_area(n)
        endif
        if (grid1_centroid_lat(n) < -pih-.01 .or.
     &      grid1_centroid_lat(n) >  pih+.01) then
          print *,'Grid 1 centroid lat error: ',n,grid1_centroid_lat(n)
        endif
        grid1_centroid_lat(n) = zero
        grid1_centroid_lon(n) = zero
      end do

      do n=1,grid2_size
        if (grid2_area(n) < -.01) then
          print *,'Grid 2 area error: ',n,grid2_area(n)
        endif
        if (grid2_centroid_lat(n) < -pih-.01 .or.
     &      grid2_centroid_lat(n) >  pih+.01) then
          print *,'Grid 2 centroid lat error: ',n,grid2_centroid_lat(n)
        endif
        grid2_centroid_lat(n) = zero
        grid2_centroid_lon(n) = zero
      end do

      do n=1,num_links_map1
        grid1_add = grid1_add_map1(n)
        grid2_add = grid2_add_map1(n)
        
        if (wts_map1(1,n) < -.01) then
          print *,'Map 1 weight < 0 ',grid1_add,grid2_add,wts_map1(1,n)
        endif
        if (norm_opt /= norm_opt_none .and. wts_map1(1,n) > 1.01) then
          print *,'Map 1 weight > 1 ',grid1_add,grid2_add,wts_map1(1,n)
        endif
        grid2_centroid_lat(grid2_add) = 
     &  grid2_centroid_lat(grid2_add) + wts_map1(1,n)

        if (num_maps > 1) then
          if (wts_map2(1,n) < -.01) then
            print *,'Map 2 weight < 0 ',grid1_add,grid2_add,
     &                                  wts_map2(1,n)
          endif
          if (norm_opt /= norm_opt_none .and. wts_map2(1,n) > 1.01) then
            print *,'Map 2 weight < 0 ',grid1_add,grid2_add,
     &                                  wts_map2(1,n)
          endif
          grid1_centroid_lat(grid1_add) = 
     &    grid1_centroid_lat(grid1_add) + wts_map2(1,n)
        endif
      end do

      do n=1,grid2_size
        select case(norm_opt)
        case (norm_opt_dstarea)
          norm_factor = grid2_frac(grid2_add)
        case (norm_opt_frcarea)
          norm_factor = one
        case (norm_opt_none)
          if (luse_grid2_area) then
            norm_factor = grid2_area_in(grid2_add)
          else
            norm_factor = grid2_area(grid2_add)
          endif
        end select
        if (abs(grid2_centroid_lat(grid2_add)-norm_factor) > .01) then
          print *,'Error: sum of wts for map1 ',grid2_add,
     &            grid2_centroid_lat(grid2_add),norm_factor
        endif
      end do

      if (num_maps > 1) then
        do n=1,grid1_size
          select case(norm_opt)
          case (norm_opt_dstarea)
            norm_factor = grid1_frac(grid1_add)
          case (norm_opt_frcarea)
            norm_factor = one
          case (norm_opt_none)
            if (luse_grid1_area) then
              norm_factor = grid1_area_in(grid1_add)
            else
              norm_factor = grid1_area(grid1_add)
            endif
          end select
          if (abs(grid1_centroid_lat(grid1_add)-norm_factor) > .01) then
            print *,'Error: sum of wts for map2 ',grid1_add,
     &              grid1_centroid_lat(grid1_add),norm_factor
          endif
        end do
      endif

      deallocate(grid1_centroid_lat,
     &           grid1_centroid_lon,
     &           grid2_centroid_lat,
     &           grid2_centroid_lon)


!-----------------------------------------------------------------------

      end subroutine remap_conserv

!***********************************************************************

      subroutine intersection(location,intrsct_lat,intrsct_lon,lcoinc,
     &                        beglat, beglon, endlat, endlon, begseg,
     &                        lbegin, lrevers)

!-----------------------------------------------------------------------
!
!     this routine finds the next intersection of a destination grid 
!     line with the line segment given by beglon, endlon, etc.
!     a coincidence flag is returned if the segment is entirely 
!     coincident with an ocean grid line.  the cells in which to search
!     for an intersection must have already been restricted in the
!     calling routine.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      logical (kind=log_kind), intent(in) ::
     &     lbegin, ! flag for first integration along this segment
     &     lrevers ! flag whether segment integrated in reverse

      real (kind=dbl_kind), intent(in) :: 
     &     beglat, beglon,  ! beginning lat/lon endpoints for segment
     &     endlat, endlon   ! ending    lat/lon endpoints for segment

      real (kind=dbl_kind), dimension(2), intent(inout) :: 
     &     begseg ! begin lat/lon of full segment

!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(out) ::
     &        location  ! address in destination array containing this
                        ! segment

      logical (kind=log_kind), intent(out) ::
     &        lcoinc    ! flag segments which are entirely coincident
                        ! with a grid line

      real (kind=dbl_kind), intent(out) ::
     &     intrsct_lat, intrsct_lon ! lat/lon coords of next intersect.

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n, next_n, cell, srch_corners, pole_loc

      integer (kind=int_kind), save :: 
     &     last_loc  ! save location when crossing threshold

      logical (kind=log_kind) :: 
     &     loutside  ! flags points outside grid

      logical (kind=log_kind), save :: 
     &     lthresh = .false.  ! flags segments crossing threshold bndy

      real (kind=dbl_kind) ::
     &     lon1, lon2,       ! local longitude variables for segment
     &     lat1, lat2,       ! local latitude  variables for segment
     &     grdlon1, grdlon2, ! local longitude variables for grid cell
     &     grdlat1, grdlat2, ! local latitude  variables for grid cell
     &     vec1_lat, vec1_lon, ! vectors and cross products used
     &     vec2_lat, vec2_lon, ! during grid search
     &     cross_product, 
     &     eps, offset,        ! small offset away from intersect
     &     s1, s2, determ,     ! variables used for linear solve to
     &     mat1, mat2, mat3, mat4, rhs1, rhs2  ! find intersection

      real (kind=dbl_kind), save ::
     &     intrsct_lat_off, intrsct_lon_off ! lat/lon coords offset 
                                            ! for next search

!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      location = 0
      lcoinc = .false.
      intrsct_lat = endlat
      intrsct_lon = endlon

      if (num_srch_cells == 0) return

      if (beglat > north_thresh .or. beglat < south_thresh) then

        if (lthresh) location = last_loc
        call pole_intersection(location,
     &               intrsct_lat,intrsct_lon,lcoinc,lthresh,
     &               beglat, beglon, endlat, endlon, begseg, lrevers)
        if (lthresh) then
          last_loc = location
          intrsct_lat_off = intrsct_lat
          intrsct_lon_off = intrsct_lon
        endif
        return

      endif

      loutside = .false.
      if (lbegin) then
        lat1 = beglat
        lon1 = beglon
      else
        lat1 = intrsct_lat_off
        lon1 = intrsct_lon_off
      endif
      lat2 = endlat
      lon2 = endlon
      if ((lon2-lon1) > three*pih) then
        lon2 = lon2 - pi2
      else if ((lon2-lon1) < -three*pih) then
        lon2 = lon2 + pi2
      endif
      s1 = zero

!-----------------------------------------------------------------------
!
!     search for location of this segment in ocean grid using cross
!     product method to determine whether a point is enclosed by a cell
!
!-----------------------------------------------------------------------

      call timer_start(12)
      srch_corners = size(srch_corner_lat,DIM=1)
      srch_loop: do

        !***
        !*** if last segment crossed threshold, use that location
        !***

        if (lthresh) then
          do cell=1,num_srch_cells
            if (srch_add(cell) == last_loc) then
              location = last_loc
              eps = tiny
              exit srch_loop
            endif
          end do
        endif

        !***
        !*** otherwise normal search algorithm
        !***

        cell_loop: do cell=1,num_srch_cells
          corner_loop: do n=1,srch_corners
            next_n = MOD(n,srch_corners) + 1

            !***
            !*** here we take the cross product of the vector making 
            !*** up each cell side with the vector formed by the vertex
            !*** and search point.  if all the cross products are 
            !*** positive, the point is contained in the cell.
            !***

            vec1_lat = srch_corner_lat(next_n,cell) - 
     &                 srch_corner_lat(n     ,cell)
            vec1_lon = srch_corner_lon(next_n,cell) - 
     &                 srch_corner_lon(n     ,cell)
            vec2_lat = lat1 - srch_corner_lat(n,cell)
            vec2_lon = lon1 - srch_corner_lon(n,cell)

            !***
            !*** if endpoint coincident with vertex, offset
            !*** the endpoint
            !***

            if (vec2_lat == 0 .and. vec2_lon == 0) then
              lat1 = lat1 + 1.d-10*(lat2-lat1)
              lon1 = lon1 + 1.d-10*(lon2-lon1)
              vec2_lat = lat1 - srch_corner_lat(n,cell)
              vec2_lon = lon1 - srch_corner_lon(n,cell)
            endif

            !***
            !*** check for 0,2pi crossings
            !***

            if (vec1_lon >  pi) then
              vec1_lon = vec1_lon - pi2
            else if (vec1_lon < -pi) then
              vec1_lon = vec1_lon + pi2
            endif
            if (vec2_lon >  pi) then
              vec2_lon = vec2_lon - pi2
            else if (vec2_lon < -pi) then
              vec2_lon = vec2_lon + pi2
            endif

            cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat

            !***
            !*** if the cross product for a side is zero, the point 
            !***   lies exactly on the side or the side is degenerate
            !***   (zero length).  if degenerate, set the cross 
            !***   product to a positive number.  otherwise perform 
            !***   another cross product between the side and the 
            !***   segment itself. 
            !*** if this cross product is also zero, the line is 
            !***   coincident with the cell boundary - perform the 
            !***   dot product and only choose the cell if the dot 
            !***   product is positive (parallel vs anti-parallel).
            !***

            if (cross_product == zero) then
              if (vec1_lat /= zero .or. vec1_lon /= zero) then
                vec2_lat = lat2 - lat1
                vec2_lon = lon2 - lon1

                if (vec2_lon >  pi) then
                  vec2_lon = vec2_lon - pi2
                else if (vec2_lon < -pi) then
                  vec2_lon = vec2_lon + pi2
                endif

                cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat
              else
                cross_product = one
              endif

              if (cross_product == zero) then
                lcoinc = .true.
                cross_product = vec1_lon*vec2_lon + vec1_lat*vec2_lat
                if (lrevers) cross_product = -cross_product
              endif
            endif

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***

            if (cross_product < zero) exit corner_loop

          end do corner_loop

          !***
          !*** if cross products all positive, we found the location
          !***

          if (n > srch_corners) then
            location = srch_add(cell)

            !***
            !*** if the beginning of this segment was outside the
            !*** grid, invert the segment so the intersection found
            !*** will be the first intersection with the grid
            !***

            if (loutside) then
              lat2 = beglat
              lon2 = beglon
              location = 0
              eps  = -tiny
            else
              eps  = tiny
            endif

            exit srch_loop
          endif

          !***
          !*** otherwise move on to next cell
          !***

        end do cell_loop

        !***
        !*** if still no cell found, the point lies outside the grid.
        !***   take some baby steps along the segment to see if any
        !***   part of the segment lies inside the grid.  
        !***

        loutside = .true.
        s1 = s1 + 0.001_dbl_kind
        lat1 = beglat + s1*(endlat - beglat)
        lon1 = beglon + s1*(lon2   - beglon)

        !***
        !*** reached the end of the segment and still outside the grid
        !*** return no intersection
        !***

        if (s1 >= one) return

      end do srch_loop
      call timer_stop(12)

!-----------------------------------------------------------------------
!
!     now that a cell is found, search for the next intersection.
!     loop over sides of the cell to find intersection with side
!     must check all sides for coincidences or intersections
!
!-----------------------------------------------------------------------

      call timer_start(13)
      intrsct_loop: do n=1,srch_corners
        next_n = mod(n,srch_corners) + 1

        grdlon1 = srch_corner_lon(n     ,cell)
        grdlon2 = srch_corner_lon(next_n,cell)
        grdlat1 = srch_corner_lat(n     ,cell)
        grdlat2 = srch_corner_lat(next_n,cell)

        !***
        !*** set up linear system to solve for intersection
        !***

        mat1 = lat2 - lat1
        mat2 = grdlat1 - grdlat2
        mat3 = lon2 - lon1
        mat4 = grdlon1 - grdlon2
        rhs1 = grdlat1 - lat1
        rhs2 = grdlon1 - lon1

        if (mat3 >  pi) then
          mat3 = mat3 - pi2
        else if (mat3 < -pi) then
          mat3 = mat3 + pi2
        endif
        if (mat4 >  pi) then
          mat4 = mat4 - pi2
        else if (mat4 < -pi) then
          mat4 = mat4 + pi2
        endif
        if (rhs2 >  pi) then
          rhs2 = rhs2 - pi2
        else if (rhs2 < -pi) then
          rhs2 = rhs2 + pi2
        endif

        determ = mat1*mat4 - mat2*mat3

        !***
        !*** if the determinant is zero, the segments are either 
        !***   parallel or coincident.  coincidences were detected 
        !***   above so do nothing.
        !*** if the determinant is non-zero, solve for the linear 
        !***   parameters s for the intersection point on each line 
        !***   segment.
        !*** if 0<s1,s2<1 then the segment intersects with this side.
        !***   return the point of intersection (adding a small
        !***   number so the intersection is off the grid line).
        !***

        if (abs(determ) > 1.e-30) then

          s1 = (rhs1*mat4 - mat2*rhs2)/determ
          s2 = (mat1*rhs2 - rhs1*mat3)/determ

          if (s2 >= zero .and. s2 <= one .and.
     &        s1 >  zero. and. s1 <= one) then

            !***
            !*** recompute intersection based on full segment
            !*** so intersections are consistent for both sweeps
            !***

            if (.not. loutside) then
              mat1 = lat2 - begseg(1)
              mat3 = lon2 - begseg(2)
              rhs1 = grdlat1 - begseg(1)
              rhs2 = grdlon1 - begseg(2)
            else
              mat1 = begseg(1) - endlat
              mat3 = begseg(2) - endlon
              rhs1 = grdlat1 - endlat
              rhs2 = grdlon1 - endlon
            endif

            if (mat3 >  pi) then
              mat3 = mat3 - pi2
            else if (mat3 < -pi) then
              mat3 = mat3 + pi2
            endif
            if (rhs2 >  pi) then
              rhs2 = rhs2 - pi2
            else if (rhs2 < -pi) then
              rhs2 = rhs2 + pi2
            endif

            determ = mat1*mat4 - mat2*mat3

            !***
            !*** sometimes due to roundoff, the previous 
            !*** determinant is non-zero, but the lines
            !*** are actually coincident.  if this is the
            !*** case, skip the rest.
            !***

            if (determ /= zero) then
              s1 = (rhs1*mat4 - mat2*rhs2)/determ
              s2 = (mat1*rhs2 - rhs1*mat3)/determ

              offset = s1 + eps/determ
              if (offset > one) offset = one

              if (.not. loutside) then
                intrsct_lat = begseg(1) + mat1*s1
                intrsct_lon = begseg(2) + mat3*s1
                intrsct_lat_off = begseg(1) + mat1*offset
                intrsct_lon_off = begseg(2) + mat3*offset
              else
                intrsct_lat = endlat + mat1*s1
                intrsct_lon = endlon + mat3*s1
                intrsct_lat_off = endlat + mat1*offset
                intrsct_lon_off = endlon + mat3*offset
              endif
              exit intrsct_loop
            endif

          endif
        endif

        !***
        !*** no intersection this side, move on to next side
        !***

      end do intrsct_loop
      call timer_stop(13)

!-----------------------------------------------------------------------
!
!     if the segment crosses a pole threshold, reset the intersection
!     to be the threshold latitude.  only check if this was not a
!     threshold segment since sometimes coordinate transform can end
!     up on other side of threshold again.
!
!-----------------------------------------------------------------------

      if (lthresh) then
        if (intrsct_lat < north_thresh .or. intrsct_lat > south_thresh)
     &      lthresh = .false.
      else if (lat1 > zero .and. intrsct_lat > north_thresh) then
        intrsct_lat = north_thresh + tiny
        intrsct_lat_off = north_thresh + eps*mat1
        s1 = (intrsct_lat - begseg(1))/mat1
        intrsct_lon     = begseg(2) + s1*mat3
        intrsct_lon_off = begseg(2) + (s1+eps)*mat3
        last_loc = location
        lthresh = .true.
      else if (lat1 < zero .and. intrsct_lat < south_thresh) then
        intrsct_lat = south_thresh - tiny
        intrsct_lat_off = south_thresh + eps*mat1
        s1 = (intrsct_lat - begseg(1))/mat1
        intrsct_lon     = begseg(2) + s1*mat3
        intrsct_lon_off = begseg(2) + (s1+eps)*mat3
        last_loc = location
        lthresh = .true.
      endif

!-----------------------------------------------------------------------

      end subroutine intersection

!***********************************************************************

      subroutine pole_intersection(location,
     &                 intrsct_lat,intrsct_lon,lcoinc,lthresh,
     &                 beglat, beglon, endlat, endlon, begseg, lrevers)

!-----------------------------------------------------------------------
!
!     this routine is identical to the intersection routine except
!     that a coordinate transformation (using a Lambert azimuthal
!     equivalent projection) is performed to treat polar cells more
!     accurately.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) :: 
     &     beglat, beglon,  ! beginning lat/lon endpoints for segment
     &     endlat, endlon   ! ending    lat/lon endpoints for segment

      real (kind=dbl_kind), dimension(2), intent(inout) :: 
     &     begseg ! begin lat/lon of full segment

      logical (kind=log_kind), intent(in) ::
     &        lrevers   ! flag true if segment integrated in reverse

!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(inout) ::
     &        location  ! address in destination array containing this
                        ! segment -- also may contain last location on
                        ! entry

      logical (kind=log_kind), intent(out) ::
     &        lcoinc    ! flag segment coincident with grid line

      logical (kind=log_kind), intent(inout) ::
     &        lthresh   ! flag segment crossing threshold boundary

      real (kind=dbl_kind), intent(out) ::
     &     intrsct_lat, intrsct_lon ! lat/lon coords of next intersect.

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n, next_n, cell, srch_corners, pole_loc

      logical (kind=log_kind) :: loutside ! flags points outside grid

      real (kind=dbl_kind) :: pi4, rns, ! north/south conversion
     &     x1, x2,       ! local x variables for segment
     &     y1, y2,       ! local y variables for segment
     &     begx, begy,   ! beginning x,y variables for segment
     &     endx, endy,   ! beginning x,y variables for segment
     &     begsegx, begsegy,   ! beginning x,y variables for segment
     &     grdx1, grdx2, ! local x variables for grid cell
     &     grdy1, grdy2, ! local y variables for grid cell
     &     vec1_y, vec1_x, ! vectors and cross products used
     &     vec2_y, vec2_x, ! during grid search
     &     cross_product, eps, ! eps=small offset away from intersect
     &     s1, s2, determ,     ! variables used for linear solve to
     &     mat1, mat2, mat3, mat4, rhs1, rhs2  ! find intersection

      real (kind=dbl_kind), dimension(:,:), allocatable ::
     &     srch_corner_x,  ! x of each corner of srch cells
     &     srch_corner_y   ! y of each corner of srch cells

      !***
      !*** save last intersection to avoid roundoff during coord
      !*** transformation
      !***

      logical (kind=log_kind), save :: luse_last = .false.

      real (kind=dbl_kind), save :: 
     &     intrsct_x, intrsct_y  ! x,y for intersection

      !***
      !*** variables necessary if segment manages to hit pole
      !***

      integer (kind=int_kind), save :: 
     &     avoid_pole_count = 0  ! count attempts to avoid pole

      real (kind=dbl_kind), save :: 
     &     avoid_pole_offset = tiny  ! endpoint offset to avoid pole

!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      if (.not. lthresh) location = 0
      lcoinc = .false.
      intrsct_lat = endlat
      intrsct_lon = endlon

      loutside = .false.
      s1 = zero

!-----------------------------------------------------------------------
!
!     convert coordinates
!
!-----------------------------------------------------------------------

      allocate(srch_corner_x(size(srch_corner_lat,DIM=1),
     &                       size(srch_corner_lat,DIM=2)),
     &         srch_corner_y(size(srch_corner_lat,DIM=1),
     &                       size(srch_corner_lat,DIM=2)))

      if (beglat > zero) then
        pi4 = quart*pi
        rns = one
      else
        pi4 = -quart*pi
        rns = -one
      endif

      if (luse_last) then
        x1 = intrsct_x
        y1 = intrsct_y
      else
        x1 = rns*two*sin(pi4 - half*beglat)*cos(beglon)
        y1 =     two*sin(pi4 - half*beglat)*sin(beglon)
        luse_last = .true.
      endif
      x2 = rns*two*sin(pi4 - half*endlat)*cos(endlon)
      y2 =     two*sin(pi4 - half*endlat)*sin(endlon)
      srch_corner_x = rns*two*sin(pi4 - half*srch_corner_lat)*
     &                        cos(srch_corner_lon)
      srch_corner_y =     two*sin(pi4 - half*srch_corner_lat)*
     &                        sin(srch_corner_lon)

      begx = x1
      begy = y1
      endx = x2
      endy = y2
      begsegx = rns*two*sin(pi4 - half*begseg(1))*cos(begseg(2))
      begsegy =     two*sin(pi4 - half*begseg(1))*sin(begseg(2))
      intrsct_x = endx
      intrsct_y = endy

!-----------------------------------------------------------------------
!
!     search for location of this segment in ocean grid using cross
!     product method to determine whether a point is enclosed by a cell
!
!-----------------------------------------------------------------------

      call timer_start(12)
      srch_corners = size(srch_corner_lat,DIM=1)
      srch_loop: do

        !***
        !*** if last segment crossed threshold, use that location
        !***

        if (lthresh) then
          do cell=1,num_srch_cells
            if (srch_add(cell) == location) then
              eps = tiny
              exit srch_loop
            endif
          end do
        endif

        !***
        !*** otherwise normal search algorithm
        !***

        cell_loop: do cell=1,num_srch_cells
          corner_loop: do n=1,srch_corners
            next_n = MOD(n,srch_corners) + 1

            !***
            !*** here we take the cross product of the vector making 
            !*** up each cell side with the vector formed by the vertex
            !*** and search point.  if all the cross products are 
            !*** positive, the point is contained in the cell.
            !***

            vec1_x = srch_corner_x(next_n,cell) - 
     &               srch_corner_x(n     ,cell)
            vec1_y = srch_corner_y(next_n,cell) - 
     &               srch_corner_y(n     ,cell)
            vec2_x = x1 - srch_corner_x(n,cell)
            vec2_y = y1 - srch_corner_y(n,cell)

            !***
            !*** if endpoint coincident with vertex, offset
            !*** the endpoint
            !***

            if (vec2_x == 0 .and. vec2_y == 0) then
              x1 = x1 + 1.d-10*(x2-x1)
              y1 = y1 + 1.d-10*(y2-y1)
              vec2_x = x1 - srch_corner_x(n,cell)
              vec2_y = y1 - srch_corner_y(n,cell)
            endif

            cross_product = vec1_x*vec2_y - vec2_x*vec1_y

            !***
            !*** if the cross product for a side is zero, the point 
            !***   lies exactly on the side or the length of a side
            !***   is zero.  if the length is zero set det > 0.
            !***   otherwise, perform another cross 
            !***   product between the side and the segment itself. 
            !*** if this cross product is also zero, the line is 
            !***   coincident with the cell boundary - perform the 
            !***   dot product and only choose the cell if the dot 
            !***   product is positive (parallel vs anti-parallel).
            !***

            if (cross_product == zero) then
              if (vec1_x /= zero .or. vec1_y /= 0) then
                vec2_x = x2 - x1
                vec2_y = y2 - y1
                cross_product = vec1_x*vec2_y - vec2_x*vec1_y
              else
                cross_product = one
              endif

              if (cross_product == zero) then
                lcoinc = .true.
                cross_product = vec1_x*vec2_x + vec1_y*vec2_y
                if (lrevers) cross_product = -cross_product
              endif
            endif

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***

            if (cross_product < zero) exit corner_loop

          end do corner_loop

          !***
          !*** if cross products all positive, we found the location
          !***

          if (n > srch_corners) then
            location = srch_add(cell)

            !***
            !*** if the beginning of this segment was outside the
            !*** grid, invert the segment so the intersection found
            !*** will be the first intersection with the grid
            !***

            if (loutside) then
              x2 = begx
              y2 = begy
              location = 0
              eps  = -tiny
            else
              eps  = tiny
            endif

            exit srch_loop
          endif

          !***
          !*** otherwise move on to next cell
          !***

        end do cell_loop

        !***
        !*** if no cell found, the point lies outside the grid.
        !***   take some baby steps along the segment to see if any
        !***   part of the segment lies inside the grid.  
        !***

        loutside = .true.
        s1 = s1 + 0.001_dbl_kind
        x1 = begx + s1*(x2 - begx)
        y1 = begy + s1*(y2 - begy)

        !***
        !*** reached the end of the segment and still outside the grid
        !*** return no intersection
        !***

        if (s1 >= one) then
          deallocate(srch_corner_x, srch_corner_y)
          luse_last = .false.
          return
        endif

      end do srch_loop
      call timer_stop(12)

!-----------------------------------------------------------------------
!
!     now that a cell is found, search for the next intersection.
!     loop over sides of the cell to find intersection with side
!     must check all sides for coincidences or intersections
!
!-----------------------------------------------------------------------

      call timer_start(13)
      intrsct_loop: do n=1,srch_corners
        next_n = mod(n,srch_corners) + 1

        grdy1 = srch_corner_y(n     ,cell)
        grdy2 = srch_corner_y(next_n,cell)
        grdx1 = srch_corner_x(n     ,cell)
        grdx2 = srch_corner_x(next_n,cell)

        !***
        !*** set up linear system to solve for intersection
        !***

        mat1 = x2 - x1
        mat2 = grdx1 - grdx2
        mat3 = y2 - y1
        mat4 = grdy1 - grdy2
        rhs1 = grdx1 - x1
        rhs2 = grdy1 - y1

        determ = mat1*mat4 - mat2*mat3

        !***
        !*** if the determinant is zero, the segments are either 
        !***   parallel or coincident or one segment has zero length.  
        !***   coincidences were detected above so do nothing.
        !*** if the determinant is non-zero, solve for the linear 
        !***   parameters s for the intersection point on each line 
        !***   segment.
        !*** if 0<s1,s2<1 then the segment intersects with this side.
        !***   return the point of intersection (adding a small
        !***   number so the intersection is off the grid line).
        !***

        if (abs(determ) > 1.e-30) then

          s1 = (rhs1*mat4 - mat2*rhs2)/determ
          s2 = (mat1*rhs2 - rhs1*mat3)/determ

          if (s2 >= zero .and. s2 <= one .and.
     &        s1 >  zero. and. s1 <= one) then

            !***
            !*** recompute intersection using entire segment
            !*** for consistency between sweeps
            !***

            if (.not. loutside) then
              mat1 = x2 - begsegx
              mat3 = y2 - begsegy
              rhs1 = grdx1 - begsegx
              rhs2 = grdy1 - begsegy
            else 
              mat1 = x2 - endx
              mat3 = y2 - endy
              rhs1 = grdx1 - endx
              rhs2 = grdy1 - endy
            endif

            determ = mat1*mat4 - mat2*mat3

            !***
            !*** sometimes due to roundoff, the previous 
            !*** determinant is non-zero, but the lines
            !*** are actually coincident.  if this is the
            !*** case, skip the rest.
            !***

            if (determ /= zero) then
              s1 = (rhs1*mat4 - mat2*rhs2)/determ
              s2 = (mat1*rhs2 - rhs1*mat3)/determ

              if (.not. loutside) then
                intrsct_x = begsegx + s1*mat1
                intrsct_y = begsegy + s1*mat3
              else 
                intrsct_x = endx + s1*mat1
                intrsct_y = endy + s1*mat3
              endif

              !***
              !*** convert back to lat/lon coordinates
              !***

              intrsct_lon = rns*atan2(intrsct_y,intrsct_x)
              if (intrsct_lon < zero) 
     &          intrsct_lon = intrsct_lon + pi2

              if (abs(intrsct_x) > 1.d-10) then
                intrsct_lat = (pi4 - 
     &            asin(rns*half*intrsct_x/cos(intrsct_lon)))*two
              else if (abs(intrsct_y) > 1.d-10) then
                intrsct_lat = (pi4 - 
     &            asin(half*intrsct_y/sin(intrsct_lon)))*two
              else
                intrsct_lat = two*pi4
              endif

              !***
              !*** add offset in transformed space for next pass.
              !***

              if (s1 - eps/determ < one) then
                intrsct_x = intrsct_x - mat1*(eps/determ)
                intrsct_y = intrsct_y - mat3*(eps/determ)
              else
                if (.not. loutside) then
                  intrsct_x = endx
                  intrsct_y = endy
                  intrsct_lat = endlat
                  intrsct_lon = endlon
                else 
                  intrsct_x = begsegx
                  intrsct_y = begsegy
                  intrsct_lat = begseg(1)
                  intrsct_lon = begseg(2)
                endif
              endif

              exit intrsct_loop
            endif
          endif
        endif

        !***
        !*** no intersection this side, move on to next side
        !***

      end do intrsct_loop
      call timer_stop(13)

      deallocate(srch_corner_x, srch_corner_y)

!-----------------------------------------------------------------------
!
!     if segment manages to cross over pole, shift the beginning 
!     endpoint in order to avoid hitting pole directly
!     (it is ok for endpoint to be pole point)
!
!-----------------------------------------------------------------------

      if (abs(intrsct_x) < 1.e-10 .and. abs(intrsct_y) < 1.e-10 .and.
     &    (endx /= zero .and. endy /=0)) then
        if (avoid_pole_count > 2) then
           avoid_pole_count = 0
           avoid_pole_offset = 10.*avoid_pole_offset
        endif

        cross_product = begsegx*(endy-begsegy) - begsegy*(endx-begsegx)
        intrsct_lat = begseg(1)
        if (cross_product*intrsct_lat > zero) then
          intrsct_lon = beglon    + avoid_pole_offset
          begseg(2)   = begseg(2) + avoid_pole_offset
        else
          intrsct_lon = beglon    - avoid_pole_offset
          begseg(2)   = begseg(2) - avoid_pole_offset
        endif

        avoid_pole_count = avoid_pole_count + 1
        luse_last = .false.
      else
        avoid_pole_count = 0
        avoid_pole_offset = tiny
      endif

!-----------------------------------------------------------------------
!
!     if the segment crosses a pole threshold, reset the intersection
!     to be the threshold latitude and do not reuse x,y intersect
!     on next entry.  only check if did not cross threshold last
!     time - sometimes the coordinate transformation can place a
!     segment on the other side of the threshold again
!
!-----------------------------------------------------------------------

      if (lthresh) then
        if (intrsct_lat > north_thresh .or. intrsct_lat < south_thresh)
     &    lthresh = .false.
      else if (beglat > zero .and. intrsct_lat < north_thresh) then
        mat4 = endlat - begseg(1)
        mat3 = endlon - begseg(2)
        if (mat3 >  pi) mat3 = mat3 - pi2
        if (mat3 < -pi) mat3 = mat3 + pi2
        intrsct_lat = north_thresh - tiny
        s1 = (north_thresh - begseg(1))/mat4
        intrsct_lon = begseg(2) + s1*mat3
        luse_last = .false.
        lthresh = .true.
      else if (beglat < zero .and. intrsct_lat > south_thresh) then
        mat4 = endlat - begseg(1)
        mat3 = endlon - begseg(2)
        if (mat3 >  pi) mat3 = mat3 - pi2
        if (mat3 < -pi) mat3 = mat3 + pi2
        intrsct_lat = south_thresh + tiny
        s1 = (south_thresh - begseg(1))/mat4
        intrsct_lon = begseg(2) + s1*mat3
        luse_last = .false.
        lthresh = .true.
      endif

      !***
      !*** if reached end of segment, do not use x,y intersect 
      !*** on next entry
      !***

      if (intrsct_lat == endlat .and. intrsct_lon == endlon) then
        luse_last = .false.
      endif

!-----------------------------------------------------------------------

      end subroutine pole_intersection

!***********************************************************************

      subroutine line_integral(weights, num_wts, 
     &                       in_phi1, in_phi2, theta1, theta2,
     &                       grid1_lat, grid1_lon, grid2_lat, grid2_lon)

!-----------------------------------------------------------------------
!
!     this routine computes the line integral of the flux function 
!     that results in the interpolation weights.  the line is defined
!     by the input lat/lon of the endpoints.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in):
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(in) ::
     &        num_wts  ! number of weights to compute

      real (kind=dbl_kind), intent(in) :: 
     &     in_phi1, in_phi2,     ! longitude endpoints for the segment
     &     theta1, theta2,       ! latitude  endpoints for the segment
     &     grid1_lat, grid1_lon, ! reference coordinates for each
     &     grid2_lat, grid2_lon  ! grid (to ensure correct 0,2pi interv.

!-----------------------------------------------------------------------
!
!     intent(out):
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), dimension(2*num_wts), intent(out) ::
     &     weights   ! line integral contribution to weights

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind) :: dphi, sinth1, sinth2, costh1, costh2, fac,
     &                        phi1, phi2, phidiff1, phidiff2, sinint
      real (kind=dbl_kind) :: f1, f2, fint

!-----------------------------------------------------------------------
!
!     weights for the general case based on a trapezoidal approx to
!     the integrals.
!
!-----------------------------------------------------------------------

      sinth1 = SIN(theta1)
      sinth2 = SIN(theta2)
      costh1 = COS(theta1)
      costh2 = COS(theta2)

      dphi = in_phi1 - in_phi2
      if (dphi >  pi) then
        dphi = dphi - pi2
      else if (dphi < -pi) then
        dphi = dphi + pi2
      endif
      dphi = half*dphi

!-----------------------------------------------------------------------
!
!     the first weight is the area overlap integral. the second and
!     fourth are second-order latitude gradient weights.
!
!-----------------------------------------------------------------------

      weights(        1) = dphi*(sinth1 + sinth2)
      weights(num_wts+1) = dphi*(sinth1 + sinth2)
      weights(        2) = dphi*(costh1 + costh2 + (theta1*sinth1 +
     &                                              theta2*sinth2))
      weights(num_wts+2) = dphi*(costh1 + costh2 + (theta1*sinth1 +
     &                                              theta2*sinth2))

!-----------------------------------------------------------------------
!
!     the third and fifth weights are for the second-order phi gradient
!     component.  must be careful of longitude range.
!
!-----------------------------------------------------------------------

      f1 = half*(costh1*sinth1 + theta1)
      f2 = half*(costh2*sinth2 + theta2)

      phi1 = in_phi1 - grid1_lon
      if (phi1 >  pi) then
        phi1 = phi1 - pi2
      else if (phi1 < -pi) then
        phi1 = phi1 + pi2
      endif

      phi2 = in_phi2 - grid1_lon
      if (phi2 >  pi) then
        phi2 = phi2 - pi2
      else if (phi2 < -pi) then
        phi2 = phi2 + pi2
      endif

      if ((phi2-phi1) <  pi .and. (phi2-phi1) > -pi) then
        weights(3) = dphi*(phi1*f1 + phi2*f2)
      else
        if (phi1 > zero) then
          fac = pi
        else
          fac = -pi
        endif
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi)
        weights(3) = half*phi1*(phi1-fac)*f1 -
     &               half*phi2*(phi2+fac)*f2 +
     &               half*fac*(phi1+phi2)*fint
      endif

      phi1 = in_phi1 - grid2_lon
      if (phi1 >  pi) then
        phi1 = phi1 - pi2
      else if (phi1 < -pi) then
        phi1 = phi1 + pi2
      endif

      phi2 = in_phi2 - grid2_lon
      if (phi2 >  pi) then
        phi2 = phi2 - pi2
      else if (phi2 < -pi) then
        phi2 = phi2 + pi2
      endif

      if ((phi2-phi1) <  pi .and. (phi2-phi1) > -pi) then
        weights(num_wts+3) = dphi*(phi1*f1 + phi2*f2)
      else
        if (phi1 > zero) then
          fac = pi
        else
          fac = -pi
        endif
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi)
        weights(num_wts+3) = half*phi1*(phi1-fac)*f1 -
     &                       half*phi2*(phi2+fac)*f2 +
     &                       half*fac*(phi1+phi2)*fint
      endif

!-----------------------------------------------------------------------

      end subroutine line_integral

!***********************************************************************

      subroutine store_link_cnsrv(add1, add2, weights)

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
     &        add2   ! address on grid2

      real (kind=dbl_kind), dimension(:), intent(in) ::
     &        weights ! array of remapping weights for this link

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: nlink, min_link, max_link ! link index

      integer (kind=int_kind), dimension(:,:), allocatable, save ::
     &        link_add1,  ! min,max link add to restrict search
     &        link_add2   ! min,max link add to restrict search

      logical (kind=log_kind), save :: first_call = .true.

!-----------------------------------------------------------------------
!
!     if all weights are zero, do not bother storing the link
!
!-----------------------------------------------------------------------

      if (all(weights == zero)) return

!-----------------------------------------------------------------------
!
!     restrict the range of links to search for existing links
!
!-----------------------------------------------------------------------

      if (first_call) then
        allocate(link_add1(2,grid1_size), link_add2(2,grid2_size))
        link_add1 = 0
        link_add2 = 0
        first_call = .false.
        min_link = 1
        max_link = 0
      else
        min_link = min(link_add1(1,add1),link_add2(1,add2))
        max_link = max(link_add1(2,add1),link_add2(2,add2))
        if (min_link == 0) then
          min_link = 1
          max_link = 0
        endif
      endif

!-----------------------------------------------------------------------
!
!     if the link already exists, add the weight to the current weight
!     arrays
!
!-----------------------------------------------------------------------

      do nlink=min_link,max_link
        if (add1 == grid1_add_map1(nlink)) then
        if (add2 == grid2_add_map1(nlink)) then

          wts_map1(:,nlink) = wts_map1(:,nlink) + weights(1:num_wts)
          if (num_maps == 2) then
            wts_map2(:,nlink) = wts_map2(:,nlink) + 
     &                                  weights(num_wts+1:2*num_wts)
          endif
          return

        endif
        endif
      end do

!-----------------------------------------------------------------------
!
!     if the link does not yet exist, increment number of links and 
!     check to see if remap arrays need to be increased to accomodate 
!     the new link.  then store the link.
!
!-----------------------------------------------------------------------

      num_links_map1  = num_links_map1 + 1
      write(*,*) 'jw adding a new link, number #',num_links_map1
      if (num_links_map1 > max_links_map1) 
     &   call resize_remap_vars(1,resize_increment)

      grid1_add_map1(num_links_map1) = add1
      grid2_add_map1(num_links_map1) = add2
      wts_map1    (:,num_links_map1) = weights(1:num_wts)

      if (num_maps > 1) then
        num_links_map2  = num_links_map2 + 1
        if (num_links_map2 > max_links_map2) 
     &     call resize_remap_vars(2,resize_increment)

        grid1_add_map2(num_links_map2) = add1
        grid2_add_map2(num_links_map2) = add2
        wts_map2    (:,num_links_map2) = weights(num_wts+1:2*num_wts)
      endif

      if (link_add1(1,add1) == 0) link_add1(1,add1) = num_links_map1
      if (link_add2(1,add2) == 0) link_add2(1,add2) = num_links_map1
      link_add1(2,add1) = num_links_map1
      link_add2(2,add2) = num_links_map1

!-----------------------------------------------------------------------

      end subroutine store_link_cnsrv

!***********************************************************************

      end module remap_conservative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
