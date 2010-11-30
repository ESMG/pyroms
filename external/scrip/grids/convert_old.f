!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This file converts a POP grid.dat file to a remapping grid file
!     in netCDF format.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: convert_old.f,v 1.2 2000/04/19 22:05:57 pwjones Exp $
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
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
!***********************************************************************

      program convertPOPT

!-----------------------------------------------------------------------
!
!     This file converts a POP grid.dat file to a remapping grid file.
!
!-----------------------------------------------------------------------

      use kinds_mod
      use constants
      use iounits
      use netcdf_mod

      implicit none

!-----------------------------------------------------------------------
!
!     variables that describe the grid
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) ::
     &             grid_size, grid_rank, grid_corners

      integer (kind=int_kind), dimension(2) ::
     &             grid_dims   ! size of each dimension

      character(char_len) :: 
     &             grid_name

      character(char_len), parameter :: 
     &             grid_file_in  = 'remap_grid_WWice.dat',
     &             grid_file_out = 'remap_grid_WWice.nc'

!-----------------------------------------------------------------------
!
!     grid coordinates and masks
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), dimension(:), allocatable ::
     &             grid_imask

      real (kind=dbl_kind), dimension(:), allocatable ::
     &             grid_center_lat,  ! lat/lon coordinates for
     &             grid_center_lon   ! each grid center in radians

      real (kind=dbl_kind), dimension(:,:), allocatable ::
     &             grid_corner_lat,  ! lat/lon coordinates for
     &             grid_corner_lon   ! each grid corner in radians

!-----------------------------------------------------------------------
!
!     other local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: i, j, n, iunit, ocn_add, im1, jm1

      integer (kind=int_kind) ::
     &        ncstat,            ! general netCDF status variable
     &        nc_grid_id,        ! netCDF grid dataset id
     &        nc_gridsize_id,    ! netCDF grid size dim id
     &        nc_gridcorn_id,    ! netCDF grid corner dim id
     &        nc_gridrank_id,    ! netCDF grid rank dim id
     &        nc_griddims_id,    ! netCDF grid dimensions id
     &        nc_grdcntrlat_id,  ! netCDF grid center lat id
     &        nc_grdcntrlon_id,  ! netCDF grid center lon id
     &        nc_grdimask_id,    ! netCDF grid mask id
     &        nc_grdcrnrlat_id,  ! netCDF grid corner lat id
     &        nc_grdcrnrlon_id   ! netCDF grid corner lon id

      integer (kind=int_kind), dimension(2) ::
     &        nc_dims2_id        ! netCDF dim id array for 2-d arrays

      real (kind=dbl_kind) :: tmplon

!-----------------------------------------------------------------------
!
!     read in grid info
!     lat/lon info is on velocity points which correspond
!     to the NE corner (in logical space) of the grid cell.
!
!-----------------------------------------------------------------------

      call get_unit(iunit)
      open(unit=iunit, file=grid_file_in, status='old', 
     &     form='unformatted')

      read(iunit) grid_name
      read(iunit) grid_size, grid_corners, grid_rank, grid_dims

      allocate( grid_center_lat(grid_size),
     &		grid_center_lon(grid_size),
     &		grid_imask     (grid_size),
     &		grid_corner_lat(grid_corners, grid_size),
     &		grid_corner_lon(grid_corners, grid_size) )

      read(iunit) grid_center_lat
      read(iunit) grid_center_lon
      read(iunit) grid_corner_lat
      read(iunit) grid_corner_lon
      read(iunit) grid_imask
      call release_unit(iunit)

!-----------------------------------------------------------------------
!
!     set up attributes for netCDF file
!
!-----------------------------------------------------------------------

      !***
      !*** create netCDF dataset for this grid
      !***

      ncstat = nf_create (grid_file_out, NF_CLOBBER,
     &                    nc_grid_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_grid_id, NF_GLOBAL, 'title',
     &                          len_trim(grid_name), grid_name)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid size dimension
      !***

      ncstat = nf_def_dim (nc_grid_id, 'grid_size', grid_size, 
     &                     nc_gridsize_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid rank dimension
      !***

      ncstat = nf_def_dim (nc_grid_id, 'grid_rank', grid_rank, 
     &                     nc_gridrank_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid corner dimension
      !***

      ncstat = nf_def_dim (nc_grid_id, 'grid_corners', grid_corners, 
     &                     nc_gridcorn_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid dim size array
      !***

      ncstat = nf_def_var (nc_grid_id, 'grid_dims', NF_INT,
     &                     1, nc_gridrank_id, nc_griddims_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid center latitude array
      !***

      ncstat = nf_def_var (nc_grid_id, 'grid_center_lat', NF_DOUBLE,
     &                     1, nc_gridsize_id, nc_grdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcntrlat_id, 'units',
     &                          7, 'radians')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid center longitude array
      !***

      ncstat = nf_def_var (nc_grid_id, 'grid_center_lon', NF_DOUBLE,
     &                     1, nc_gridsize_id, nc_grdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcntrlon_id, 'units',
     &                          7, 'radians')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid mask
      !***

      ncstat = nf_def_var (nc_grid_id, 'grid_imask', NF_INT,
     &                     1, nc_gridsize_id, nc_grdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdimask_id, 'units',
     &                          8, 'unitless')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid corner latitude array
      !***

      nc_dims2_id(1) = nc_gridcorn_id
      nc_dims2_id(2) = nc_gridsize_id

      ncstat = nf_def_var (nc_grid_id, 'grid_corner_lat', NF_DOUBLE,
     &                     2, nc_dims2_id, nc_grdcrnrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcrnrlat_id, 'units',
     &                          7, 'radians')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid corner longitude array
      !***

      ncstat = nf_def_var (nc_grid_id, 'grid_corner_lon', NF_DOUBLE,
     &                     2, nc_dims2_id, nc_grdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcrnrlon_id, 'units',
     &                          7, 'radians')
      call netcdf_error_handler(ncstat)

      !***
      !*** end definition stage
      !***

      ncstat = nf_enddef(nc_grid_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     write grid data
!
!-----------------------------------------------------------------------

      ncstat = nf_put_var_int(nc_grid_id, nc_griddims_id, grid_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_grid_id, nc_grdcntrlat_id, 
     &                           grid_center_lat)
      ncstat = nf_put_var_int(nc_grid_id, nc_grdimask_id, grid_imask)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_grid_id, nc_grdcntrlat_id, 
     &                           grid_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_grid_id, nc_grdcntrlon_id, 
     &                           grid_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_grid_id, nc_grdcrnrlat_id, 
     &                           grid_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_grid_id, nc_grdcrnrlon_id, 
     &                           grid_corner_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_close(nc_grid_id)

!***********************************************************************

      end program convertPOPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
