!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module reads in and initializes two grids for remapping.
!     NOTE: grid1 must be the master grid -- the grid that determines
!           which cells participate (e.g. land mask) and the fractional
!           area of grid2 cells that participate in the remapping.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: grids.f,v 1.6 2001/08/21 21:06:41 pwjones Exp $
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

      module grids

!-----------------------------------------------------------------------

      use kinds_mod    ! defines data types
      use constants    ! common constants
      use iounits      ! I/O unit manager
      use netcdf_mod   ! netCDF stuff

      implicit none

!-----------------------------------------------------------------------
!
!     variables that describe each grid
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), save ::
     &             grid1_size, grid2_size, ! total points on each grid
     &             grid1_rank, grid2_rank, ! rank of each grid
     &             grid1_corners, grid2_corners ! number of corners
                                                ! for each grid cell

      integer (kind=int_kind), dimension(:), allocatable, save ::
     &             grid1_dims, grid2_dims  ! size of each grid dimension

      character(char_len), save :: 
     &             grid1_name, grid2_name  ! name for each grid

      character (char_len), save :: 
     &             grid1_units, ! units for grid coords (degs/radians)
     &             grid2_units  ! units for grid coords

      real (kind=dbl_kind), parameter ::
     &      deg2rad = pi/180.   ! conversion for deg to rads

!-----------------------------------------------------------------------
!
!     grid coordinates and masks
!
!-----------------------------------------------------------------------

      logical (kind=log_kind), dimension(:), allocatable, save ::
     &             grid1_mask,        ! flag which cells participate
     &             grid2_mask         ! flag which cells participate

      logical (kind=log_kind), dimension(2), save ::
     &             grid1_periodic,    ! grid periodicity in x and y dirs
     &             grid2_periodic     ! grid periodicity in x and y dirs

      real (kind=dbl_kind), dimension(:), allocatable, save ::
     &             grid1_center_lat,  ! lat/lon coordinates for
     &             grid1_center_lon,  ! each grid center in radians
     &             grid2_center_lat, 
     &             grid2_center_lon,
     &             grid1_area,        ! tot area of each grid1 cell
     &             grid2_area,        ! tot area of each grid2 cell
     &             grid1_area_in,     ! area of grid1 cell from file
     &             grid2_area_in,     ! area of grid2 cell from file
     &             grid1_frac,        ! fractional area of grid cells
     &             grid2_frac         ! participating in remapping

      real (kind=dbl_kind), dimension(:,:), allocatable, save ::
     &             grid1_corner_lat,  ! lat/lon coordinates for
     &             grid1_corner_lon,  ! each grid corner in radians
     &             grid2_corner_lat, 
     &             grid2_corner_lon

      logical (kind=log_kind), save ::
     &             luse_grid_centers ! use centers for bounding boxes
     &,            luse_grid1_area   ! use area from grid file
     &,            luse_grid2_area   ! use area from grid file

      real (kind=dbl_kind), dimension(:,:), allocatable, save ::
     &             grid1_bound_box,  ! lat/lon bounding box for use
     &             grid2_bound_box   ! in restricting grid searches

!-----------------------------------------------------------------------
!
!     bins for restricting searches
!
!-----------------------------------------------------------------------

      character (char_len), save ::
     &        restrict_type  ! type of bins to use

      integer (kind=int_kind), save ::
     &        num_srch_bins  ! num of bins for restricted srch

      integer (kind=int_kind), dimension(:,:), allocatable, save ::
     &        bin_addr1, ! min,max adds for grid1 cells in this lat bin
     &        bin_addr2  ! min,max adds for grid2 cells in this lat bin

      real(kind=dbl_kind), dimension(:,:), allocatable, save ::
     &        bin_lats   ! min,max latitude for each search bin
     &,       bin_lons   ! min,max longitude for each search bin

!***********************************************************************

      contains

!***********************************************************************

      subroutine grid_init(grid1_file, grid2_file)

!-----------------------------------------------------------------------
!
!     this routine reads grid info from grid files and makes any
!     necessary changes (e.g. for 0,2pi longitude range)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(char_len), intent(in) :: 
     &             grid1_file, grid2_file  ! grid data files

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: 
     &  n      ! loop counter
     &, nele   ! element loop counter
     &, iunit  ! unit number for opening files
     &, i,j    ! logical 2d addresses
     &, ip1,jp1
     &, n_add, e_add, ne_add
     &, nx, ny

      integer (kind=int_kind) :: 
     &         ncstat,           ! netCDF status variable
     &         nc_grid1_id,       ! netCDF grid file id
     &         nc_grid2_id,       ! netCDF grid file id
     &         nc_grid1size_id,   ! netCDF grid size dim id
     &         nc_grid2size_id,   ! netCDF grid size dim id
     &         nc_grid1corn_id,   ! netCDF grid corner dim id
     &         nc_grid2corn_id,   ! netCDF grid corner dim id
     &         nc_grid1rank_id,   ! netCDF grid rank dim id
     &         nc_grid2rank_id,   ! netCDF grid rank dim id
     &         nc_grid1area_id,   ! netCDF grid rank dim id
     &         nc_grid2area_id,   ! netCDF grid rank dim id
     &         nc_grid1dims_id,   ! netCDF grid dimension size id
     &         nc_grid2dims_id,   ! netCDF grid dimension size id
     &         nc_grd1imask_id,   ! netCDF grid imask var id
     &         nc_grd2imask_id,   ! netCDF grid imask var id
     &         nc_grd1crnrlat_id, ! netCDF grid corner lat var id
     &         nc_grd2crnrlat_id, ! netCDF grid corner lat var id
     &         nc_grd1crnrlon_id, ! netCDF grid corner lon var id
     &         nc_grd2crnrlon_id, ! netCDF grid corner lon var id
     &         nc_grd1cntrlat_id, ! netCDF grid center lat var id
     &         nc_grd2cntrlat_id, ! netCDF grid center lat var id
     &         nc_grd1cntrlon_id, ! netCDF grid center lon var id
     &         nc_grd2cntrlon_id  ! netCDF grid center lon var id

      integer (kind=int_kind), dimension(:), allocatable :: 
     &                            imask ! integer mask read from file

      real (kind=dbl_kind) :: 
     &  dlat,dlon           ! lat/lon intervals for search bins

      real (kind=dbl_kind), dimension(4) ::
     &  tmp_lats, tmp_lons  ! temps for computing bounding boxes

!-----------------------------------------------------------------------
!
!     open grid files and read grid size/name data
!
!-----------------------------------------------------------------------

      ncstat = nf_open(grid1_file, NF_NOWRITE, nc_grid1_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_open(grid2_file, NF_NOWRITE, nc_grid2_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_grid1_id, 'grid_size', nc_grid1size_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_grid1_id, nc_grid1size_id, grid1_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_grid2_id, 'grid_size', nc_grid2size_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_grid2_id, nc_grid2size_id, grid2_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_grid1_id, 'grid_rank', nc_grid1rank_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_grid1_id, nc_grid1rank_id, grid1_rank)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_grid2_id, 'grid_rank', nc_grid2rank_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_grid2_id, nc_grid2rank_id, grid2_rank)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_grid1_id,'grid_corners',nc_grid1corn_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_grid1_id,nc_grid1corn_id,grid1_corners)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_grid2_id,'grid_corners',nc_grid2corn_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_grid2_id,nc_grid2corn_id,grid2_corners)
      call netcdf_error_handler(ncstat)

      allocate( grid1_dims(grid1_rank),
     &          grid2_dims(grid2_rank))

      ncstat = nf_get_att_text(nc_grid1_id, nf_global, 'title',
     &                         grid1_name)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_att_text(nc_grid2_id, nf_global, 'title',
     &                         grid2_name)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     allocate grid coordinates/masks and read data
!
!-----------------------------------------------------------------------

      allocate( grid1_mask      (grid1_size),
     &          grid2_mask      (grid2_size),
     &          grid1_center_lat(grid1_size), 
     &          grid1_center_lon(grid1_size),
     &          grid2_center_lat(grid2_size), 
     &          grid2_center_lon(grid2_size),
     &          grid1_area      (grid1_size),
     &          grid2_area      (grid2_size),
     &          grid1_frac      (grid1_size),
     &          grid2_frac      (grid2_size),
     &          grid1_corner_lat(grid1_corners, grid1_size),
     &          grid1_corner_lon(grid1_corners, grid1_size),
     &          grid2_corner_lat(grid2_corners, grid2_size),
     &          grid2_corner_lon(grid2_corners, grid2_size),
     &          grid1_bound_box (4            , grid1_size),
     &          grid2_bound_box (4            , grid2_size))

      allocate(imask(grid1_size))

      ncstat = nf_inq_varid(nc_grid1_id, 'grid_dims', nc_grid1dims_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid1_id, 'grid_imask', nc_grd1imask_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid1_id, 'grid_center_lat', 
     &                                   nc_grd1cntrlat_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid1_id, 'grid_center_lon', 
     &                                   nc_grd1cntrlon_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid1_id, 'grid_corner_lat', 
     &                                   nc_grd1crnrlat_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid1_id, 'grid_corner_lon', 
     &                                   nc_grd1crnrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_grid1_id, nc_grid1dims_id, grid1_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_grid1_id, nc_grd1imask_id, imask)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_grid1_id, nc_grd1cntrlat_id, 
     &                                       grid1_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_grid1_id, nc_grd1cntrlon_id, 
     &                                       grid1_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_grid1_id, nc_grd1crnrlat_id, 
     &                                       grid1_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_grid1_id, nc_grd1crnrlon_id, 
     &                                       grid1_corner_lon)
      call netcdf_error_handler(ncstat)

      if (luse_grid1_area) then
        allocate (grid1_area_in(grid1_size))
        ncstat = nf_inq_varid(nc_grid1_id, 'grid_area', nc_grid1area_id)
        call netcdf_error_handler(ncstat)
        ncstat = nf_get_var_double(nc_grid1_id, nc_grid1area_id, 
     &                                          grid1_area_in)
        call netcdf_error_handler(ncstat)
      endif

      grid1_area = zero
      grid1_frac = zero

!-----------------------------------------------------------------------
!
!     initialize logical mask and convert lat/lon units if required
!
!-----------------------------------------------------------------------

      where (imask == 1)
        grid1_mask = .true.
      elsewhere
        grid1_mask = .false.
      endwhere
      deallocate(imask)

      grid1_units = ' '
      ncstat = nf_get_att_text(nc_grid1_id, nc_grd1cntrlat_id, 'units',
     &                         grid1_units)
      call netcdf_error_handler(ncstat)

      select case (grid1_units(1:7))
      case ('degrees')

        grid1_center_lat = grid1_center_lat*deg2rad
        grid1_center_lon = grid1_center_lon*deg2rad

      case ('radians')

        !*** no conversion necessary

      case default

        print *,'unknown units supplied for grid1 center lat/lon: '
        print *,'proceeding assuming radians'

      end select

      grid1_units = ' '
      ncstat = nf_get_att_text(nc_grid1_id, nc_grd1crnrlat_id, 'units',
     &                         grid1_units)
      call netcdf_error_handler(ncstat)

      select case (grid1_units(1:7))
      case ('degrees')

        grid1_corner_lat = grid1_corner_lat*deg2rad
        grid1_corner_lon = grid1_corner_lon*deg2rad

      case ('radians')

        !*** no conversion necessary

      case default

        print *,'unknown units supplied for grid1 corner lat/lon: '
        print *,'proceeding assuming radians'

      end select

      ncstat = nf_close(nc_grid1_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     read data for grid 2
!
!-----------------------------------------------------------------------

      allocate(imask(grid2_size))

      ncstat = nf_inq_varid(nc_grid2_id, 'grid_dims', nc_grid2dims_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid2_id, 'grid_imask', nc_grd2imask_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid2_id, 'grid_center_lat', 
     &                                   nc_grd2cntrlat_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid2_id, 'grid_center_lon', 
     &                                   nc_grd2cntrlon_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid2_id, 'grid_corner_lat', 
     &                                   nc_grd2crnrlat_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_varid(nc_grid2_id, 'grid_corner_lon', 
     &                                   nc_grd2crnrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_grid2_id, nc_grid2dims_id, grid2_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_grid2_id, nc_grd2imask_id, imask)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_grid2_id, nc_grd2cntrlat_id, 
     &                                       grid2_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_grid2_id, nc_grd2cntrlon_id, 
     &                                       grid2_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_grid2_id, nc_grd2crnrlat_id, 
     &                                       grid2_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_grid2_id, nc_grd2crnrlon_id, 
     &                                       grid2_corner_lon)
      call netcdf_error_handler(ncstat)

      if (luse_grid2_area) then
        allocate (grid2_area_in(grid2_size))
        ncstat = nf_inq_varid(nc_grid2_id, 'grid_area', nc_grid2area_id)
        call netcdf_error_handler(ncstat)
        ncstat = nf_get_var_double(nc_grid2_id, nc_grid2area_id, 
     &                                          grid2_area_in)
        call netcdf_error_handler(ncstat)
      endif

      grid2_area = zero
      grid2_frac = zero

!-----------------------------------------------------------------------
!
!     initialize logical mask and convert lat/lon units if required
!
!-----------------------------------------------------------------------

      where (imask == 1)
        grid2_mask = .true.
      elsewhere
        grid2_mask = .false.
      endwhere
      deallocate(imask)

      grid2_units = ' '
      ncstat = nf_get_att_text(nc_grid2_id, nc_grd2cntrlat_id, 'units',
     &                         grid2_units)
      call netcdf_error_handler(ncstat)

      select case (grid2_units(1:7))
      case ('degrees')

        grid2_center_lat = grid2_center_lat*deg2rad
        grid2_center_lon = grid2_center_lon*deg2rad

      case ('radians')

        !*** no conversion necessary

      case default

        print *,'unknown units supplied for grid2 center lat/lon: '
        print *,'proceeding assuming radians'

      end select

      grid2_units = ' '
      ncstat = nf_get_att_text(nc_grid2_id, nc_grd2crnrlat_id, 'units',
     &                         grid2_units)
      call netcdf_error_handler(ncstat)

      select case (grid2_units(1:7))
      case ('degrees')

        grid2_corner_lat = grid2_corner_lat*deg2rad
        grid2_corner_lon = grid2_corner_lon*deg2rad

      case ('radians')

        !*** no conversion necessary

      case default

        print *,'no units supplied for grid2 corner lat/lon: '
        print *,'proceeding assuming radians'

      end select

      ncstat = nf_close(nc_grid2_id)
      call netcdf_error_handler(ncstat)


!-----------------------------------------------------------------------
!
!     convert longitudes to 0,2pi interval
!
!-----------------------------------------------------------------------

      where (grid1_center_lon .gt. pi2)  grid1_center_lon =
     &                                   grid1_center_lon - pi2
      where (grid1_center_lon .lt. zero) grid1_center_lon =
     &                                   grid1_center_lon + pi2
      where (grid2_center_lon .gt. pi2)  grid2_center_lon =
     &                                   grid2_center_lon - pi2
      where (grid2_center_lon .lt. zero) grid2_center_lon =
     &                                   grid2_center_lon + pi2
      where (grid1_corner_lon .gt. pi2)  grid1_corner_lon =
     &                                   grid1_corner_lon - pi2
      where (grid1_corner_lon .lt. zero) grid1_corner_lon =
     &                                   grid1_corner_lon + pi2
      where (grid2_corner_lon .gt. pi2)  grid2_corner_lon =
     &                                   grid2_corner_lon - pi2
      where (grid2_corner_lon .lt. zero) grid2_corner_lon =
     &                                   grid2_corner_lon + pi2

!-----------------------------------------------------------------------
!
!     make sure input latitude range is within the machine values
!     for +/- pi/2 
!
!-----------------------------------------------------------------------

      where (grid1_center_lat >  pih) grid1_center_lat =  pih
      where (grid1_corner_lat >  pih) grid1_corner_lat =  pih
      where (grid1_center_lat < -pih) grid1_center_lat = -pih
      where (grid1_corner_lat < -pih) grid1_corner_lat = -pih

      where (grid2_center_lat >  pih) grid2_center_lat =  pih
      where (grid2_corner_lat >  pih) grid2_corner_lat =  pih
      where (grid2_center_lat < -pih) grid2_center_lat = -pih
      where (grid2_corner_lat < -pih) grid2_corner_lat = -pih

!-----------------------------------------------------------------------
!
!     compute bounding boxes for restricting future grid searches
!
!-----------------------------------------------------------------------

      if (.not. luse_grid_centers) then
        grid1_bound_box(1,:) = minval(grid1_corner_lat, DIM=1)
        grid1_bound_box(2,:) = maxval(grid1_corner_lat, DIM=1)
        grid1_bound_box(3,:) = minval(grid1_corner_lon, DIM=1)
        grid1_bound_box(4,:) = maxval(grid1_corner_lon, DIM=1)

        grid2_bound_box(1,:) = minval(grid2_corner_lat, DIM=1)
        grid2_bound_box(2,:) = maxval(grid2_corner_lat, DIM=1)
        grid2_bound_box(3,:) = minval(grid2_corner_lon, DIM=1)
        grid2_bound_box(4,:) = maxval(grid2_corner_lon, DIM=1)

      else

        nx = grid1_dims(1)
        ny = grid1_dims(2)

        do n=1,grid1_size

          !*** find N,S and NE points to this grid point

          j = (n - 1)/nx +1
          i = n - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
            !*** assume cyclic
            ip1 = 1
            !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid1_center_lat(e_add) - 
     &              grid1_center_lat(n   )) > pih) then
              ip1 = i
            endif
          endif

          if (j < ny) then
            jp1 = j+1
          else
            !*** assume cyclic
            jp1 = 1
            !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid1_center_lat(n_add) - 
     &              grid1_center_lat(n   )) > pih) then
              jp1 = j
            endif
          endif

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

          !*** find N,S and NE lat/lon coords and check bounding box

          tmp_lats(1) = grid1_center_lat(n)
          tmp_lats(2) = grid1_center_lat(e_add)
          tmp_lats(3) = grid1_center_lat(ne_add)
          tmp_lats(4) = grid1_center_lat(n_add)

          tmp_lons(1) = grid1_center_lon(n)
          tmp_lons(2) = grid1_center_lon(e_add)
          tmp_lons(3) = grid1_center_lon(ne_add)
          tmp_lons(4) = grid1_center_lon(n_add)

          grid1_bound_box(1,n) = minval(tmp_lats)
          grid1_bound_box(2,n) = maxval(tmp_lats)
          grid1_bound_box(3,n) = minval(tmp_lons)
          grid1_bound_box(4,n) = maxval(tmp_lons)
        end do

        nx = grid2_dims(1)
        ny = grid2_dims(2)

        do n=1,grid2_size

          !*** find N,S and NE points to this grid point

          j = (n - 1)/nx +1
          i = n - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
            !*** assume cyclic
            ip1 = 1
            !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid2_center_lat(e_add) - 
     &              grid2_center_lat(n   )) > pih) then
              ip1 = i
            endif
          endif

          if (j < ny) then
            jp1 = j+1
          else
            !*** assume cyclic
            jp1 = 1
            !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid2_center_lat(n_add) - 
     &              grid2_center_lat(n   )) > pih) then
              jp1 = j
            endif
          endif

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

          !*** find N,S and NE lat/lon coords and check bounding box

          tmp_lats(1) = grid2_center_lat(n)
          tmp_lats(2) = grid2_center_lat(e_add)
          tmp_lats(3) = grid2_center_lat(ne_add)
          tmp_lats(4) = grid2_center_lat(n_add)

          tmp_lons(1) = grid2_center_lon(n)
          tmp_lons(2) = grid2_center_lon(e_add)
          tmp_lons(3) = grid2_center_lon(ne_add)
          tmp_lons(4) = grid2_center_lon(n_add)

          grid2_bound_box(1,n) = minval(tmp_lats)
          grid2_bound_box(2,n) = maxval(tmp_lats)
          grid2_bound_box(3,n) = minval(tmp_lons)
          grid2_bound_box(4,n) = maxval(tmp_lons)
        end do

      endif

      where (abs(grid1_bound_box(4,:) - grid1_bound_box(3,:)) > pi)
        grid1_bound_box(3,:) = zero
        grid1_bound_box(4,:) = pi2
      end where

      where (abs(grid2_bound_box(4,:) - grid2_bound_box(3,:)) > pi)
        grid2_bound_box(3,:) = zero
        grid2_bound_box(4,:) = pi2
      end where

      !***
      !*** try to check for cells that overlap poles
      !***

      where (grid1_center_lat > grid1_bound_box(2,:))
     &  grid1_bound_box(2,:) = pih

      where (grid1_center_lat < grid1_bound_box(1,:))
     &  grid1_bound_box(1,:) = -pih

      where (grid2_center_lat > grid2_bound_box(2,:))
     &  grid2_bound_box(2,:) = pih

      where (grid2_center_lat < grid2_bound_box(1,:))
     &  grid2_bound_box(1,:) = -pih

!-----------------------------------------------------------------------
!
!     set up and assign address ranges to search bins in order to 
!     further restrict later searches
!
!-----------------------------------------------------------------------

      select case (restrict_type)

      case ('latitude')
        write(stdout,*) 'Using latitude bins to restrict search.'

        allocate(bin_addr1(2,num_srch_bins))
        allocate(bin_addr2(2,num_srch_bins))
        allocate(bin_lats (2,num_srch_bins))
        allocate(bin_lons (2,num_srch_bins))

        dlat = pi/num_srch_bins

        do n=1,num_srch_bins
          bin_lats(1,n) = (n-1)*dlat - pih
          bin_lats(2,n) =     n*dlat - pih
          bin_lons(1,n) = zero
          bin_lons(2,n) = pi2
          bin_addr1(1,n) = grid1_size + 1
          bin_addr1(2,n) = 0
          bin_addr2(1,n) = grid2_size + 1
          bin_addr2(2,n) = 0
        end do

        do nele=1,grid1_size
          do n=1,num_srch_bins
            if (grid1_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid1_bound_box(2,nele) >= bin_lats(1,n)) then
              bin_addr1(1,n) = min(nele,bin_addr1(1,n))
              bin_addr1(2,n) = max(nele,bin_addr1(2,n))
            endif
          end do
        end do

        do nele=1,grid2_size
          do n=1,num_srch_bins
            if (grid2_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid2_bound_box(2,nele) >= bin_lats(1,n)) then
              bin_addr2(1,n) = min(nele,bin_addr2(1,n))
              bin_addr2(2,n) = max(nele,bin_addr2(2,n))
            endif
          end do
        end do

      case ('latlon')
        write(stdout,*) 'Using lat/lon boxes to restrict search.'

        dlat = pi /num_srch_bins
        dlon = pi2/num_srch_bins

        allocate(bin_addr1(2,num_srch_bins*num_srch_bins))
        allocate(bin_addr2(2,num_srch_bins*num_srch_bins))
        allocate(bin_lats (2,num_srch_bins*num_srch_bins))
        allocate(bin_lons (2,num_srch_bins*num_srch_bins))

        n = 0
        do j=1,num_srch_bins
        do i=1,num_srch_bins
          n = n + 1

          bin_lats(1,n) = (j-1)*dlat - pih
          bin_lats(2,n) =     j*dlat - pih
          bin_lons(1,n) = (i-1)*dlon
          bin_lons(2,n) =     i*dlon
          bin_addr1(1,n) = grid1_size + 1
          bin_addr1(2,n) = 0
          bin_addr2(1,n) = grid2_size + 1
          bin_addr2(2,n) = 0
        end do
        end do

        num_srch_bins = num_srch_bins**2

        do nele=1,grid1_size
          do n=1,num_srch_bins
            if (grid1_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid1_bound_box(2,nele) >= bin_lats(1,n) .and.
     &          grid1_bound_box(3,nele) <= bin_lons(2,n) .and.
     &          grid1_bound_box(4,nele) >= bin_lons(1,n)) then
              bin_addr1(1,n) = min(nele,bin_addr1(1,n))
              bin_addr1(2,n) = max(nele,bin_addr1(2,n))
            endif
          end do
        end do

        do nele=1,grid2_size
          do n=1,num_srch_bins
            if (grid2_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid2_bound_box(2,nele) >= bin_lats(1,n) .and.
     &          grid2_bound_box(3,nele) <= bin_lons(2,n) .and.
     &          grid2_bound_box(4,nele) >= bin_lons(1,n)) then
              bin_addr2(1,n) = min(nele,bin_addr2(1,n))
              bin_addr2(2,n) = max(nele,bin_addr2(2,n))
            endif
          end do
        end do

      case default
        stop 'unknown search restriction method'
      end select

!-----------------------------------------------------------------------

      end subroutine grid_init

!***********************************************************************

      end module grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

