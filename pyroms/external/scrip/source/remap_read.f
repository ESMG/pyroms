!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This routine reads remapping information from files written
!     by remap_setup.  If remapping in both directions are required,
!     two input files must be specified.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_read.f,v 1.6 2000/04/19 21:56:26 pwjones Exp $
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

      module remap_read

!-----------------------------------------------------------------------
!
!     contains routines for reading a remap file
!
!-----------------------------------------------------------------------

      use kinds_mod     ! defines common data types
      use constants     ! defines useful constants
      use grids         ! includes all grid information
      use netcdf_mod    ! module with netcdf vars and utilities
      use remap_vars    ! module for all required remapping variables

      implicit none

!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     various netCDF ids for files variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), private :: ! netCDF ids
     &         ncstat, nc_file_id,
     &         nc_srcgrdsize_id, nc_dstgrdsize_id,
     &         nc_srcgrdcorn_id, nc_dstgrdcorn_id,
     &         nc_srcgrdrank_id, nc_dstgrdrank_id,
     &         nc_srcgrddims_id, nc_dstgrddims_id,
     &         nc_numlinks_id, nc_numwgts_id, 
     &         nc_srcgrdimask_id, nc_dstgrdimask_id,
     &         nc_srcgrdcntrlat_id, nc_srcgrdcntrlon_id,
     &         nc_srcgrdcrnrlat_id, nc_srcgrdcrnrlon_id,
     &         nc_srcgrdarea_id, nc_srcgrdfrac_id,
     &         nc_dstgrdcntrlat_id, nc_dstgrdcntrlon_id,
     &         nc_dstgrdcrnrlat_id, nc_dstgrdcrnrlon_id,
     &         nc_dstgrdarea_id, nc_dstgrdfrac_id,
     &         nc_srcgrdadd_id, nc_dstgrdadd_id, nc_rmpmatrix_id

!***********************************************************************

      contains

!***********************************************************************

      subroutine read_remap(map_name, interp_file)

!-----------------------------------------------------------------------
!
!     this driver routine reads some global attributes and then
!     calls a specific read routine based on file conventions
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(char_len), intent(in) ::
     &  interp_file        ! filename for remap data

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      character(char_len), intent(out) ::
     &  map_name            ! name for mapping

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character(char_len) :: 
     &   map_method       ! character string for map_type
     &,  normalize_opt    ! character string for normalization option
     &,  convention       ! character string for output convention

!-----------------------------------------------------------------------
!
!     open file and read some global information
!
!-----------------------------------------------------------------------

      ncstat = nf_open(interp_file, NF_NOWRITE, nc_file_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** map name
      !***
      map_name = ' '
      ncstat = nf_get_att_text(nc_file_id, NF_GLOBAL, 'title',
     &                         map_name)
      call netcdf_error_handler(ncstat)

      print *,'Reading remapping:',trim(map_name)
      print *,'From file:',trim(interp_file)

      !***
      !*** normalization option
      !***
      normalize_opt = ' '
      ncstat = nf_get_att_text(nc_file_id, NF_GLOBAL, 'normalization',
     &                         normalize_opt)
      call netcdf_error_handler(ncstat)

      select case(normalize_opt)
      case ('none')
        norm_opt = norm_opt_none
      case ('fracarea')
        norm_opt = norm_opt_frcarea
      case ('destarea')
        norm_opt = norm_opt_dstarea
      case default
        print *,'normalize_opt = ',normalize_opt
        stop 'Invalid normalization option'
      end select

      !***
      !*** map method
      !***
      map_method = ' '
      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'map_method',
     &                          map_method)
      call netcdf_error_handler(ncstat)

      select case(map_method)
      case('Conservative remapping')
        map_type = map_type_conserv
      case('Bilinear remapping')
        map_type = map_type_bilinear
      case('Distance weighted avg of nearest neighbors')
        map_type = map_type_distwgt
      case('Bicubic remapping')
        map_type = map_type_bicubic
      case default
        print *,'map_type = ',map_method
        stop 'Invalid Map Type'
      end select

      !***
      !*** file convention
      !***
      convention = ' '
      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'conventions',
     &                          convention)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     call appropriate read routine based on output convention
!
!-----------------------------------------------------------------------

      select case(convention)
      case ('SCRIP')
        call read_remap_scrip
      case ('NCAR-CSM')
        call read_remap_csm
      case default
        print *,'convention = ',convention
        stop 'unknown output file convention'
      end select

!-----------------------------------------------------------------------

      end subroutine read_remap

!***********************************************************************

      subroutine read_remap_scrip

!-----------------------------------------------------------------------
!
!     the routine reads a netCDF file to extract remapping info
!     in SCRIP format
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (char_len) ::
     &  grid1_name           ! grid name for source grid
     &, grid2_name           ! grid name for dest   grid

      integer (kind=int_kind) ::  
     &  n                    ! dummy index

      integer (kind=int_kind), dimension(:), allocatable ::
     &  grid1_mask_int,      ! integer masks to determine
     &  grid2_mask_int       ! cells that participate in map

!-----------------------------------------------------------------------
!
!     read some additional global attributes
!
!-----------------------------------------------------------------------

      !***
      !*** source and destination grid names
      !***

      grid1_name = ' '
      grid2_name = ' '
      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'source_grid',
     &                          grid1_name)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'dest_grid',
     &                          grid2_name)
      call netcdf_error_handler(ncstat)

      print *,' '
      print *,'Remapping between:',trim(grid1_name)
      print *,'and ',trim(grid2_name)
      print *,' '

!-----------------------------------------------------------------------
!
!     read dimension information
!
!-----------------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'src_grid_size', 
     &                      nc_srcgrdsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdsize_id, grid1_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'dst_grid_size', 
     &                      nc_dstgrdsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdsize_id, grid2_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'src_grid_corners', 
     &                      nc_srcgrdcorn_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdcorn_id, 
     &                       grid1_corners)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'dst_grid_corners', 
     &                      nc_dstgrdcorn_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdcorn_id, 
     &                       grid2_corners)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'src_grid_rank', 
     &                      nc_srcgrdrank_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdrank_id, 
     &                       grid1_rank)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'dst_grid_rank', 
     &                      nc_dstgrdrank_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdrank_id, 
     &                       grid2_rank)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'num_links', 
     &                      nc_numlinks_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numlinks_id, 
     &                       num_links_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'num_wgts', 
     &                      nc_numwgts_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numwgts_id, num_wts)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     allocate arrays
!
!-----------------------------------------------------------------------

      allocate( grid1_dims      (grid1_rank),
     &          grid1_center_lat(grid1_size), 
     &          grid1_center_lon(grid1_size),
     &          grid1_area      (grid1_size),
     &          grid1_frac      (grid1_size),
     &          grid1_mask      (grid1_size),
     &          grid1_mask_int  (grid1_size),
     &          grid1_corner_lat(grid1_corners, grid1_size),
     &          grid1_corner_lon(grid1_corners, grid1_size) )

      allocate( grid2_dims      (grid2_rank),
     &          grid2_center_lat(grid2_size), 
     &          grid2_center_lon(grid2_size),
     &          grid2_area      (grid2_size),
     &          grid2_frac      (grid2_size),
     &          grid2_mask      (grid2_size),
     &          grid2_mask_int  (grid2_size),
     &          grid2_corner_lat(grid2_corners, grid2_size),
     &          grid2_corner_lon(grid2_corners, grid2_size) )

      allocate( grid1_add_map1(num_links_map1),
     &          grid2_add_map1(num_links_map1),
     &          wts_map1(num_wts,num_links_map1) )

!-----------------------------------------------------------------------
!
!     get variable ids
!
!-----------------------------------------------------------------------

      ncstat = nf_inq_varid(nc_file_id, 'src_grid_dims', 
     &                      nc_srcgrddims_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'src_grid_imask', 
     &                      nc_srcgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'src_grid_center_lat', 
     &                                   nc_srcgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'src_grid_center_lon', 
     &                                   nc_srcgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'src_grid_corner_lat', 
     &                                   nc_srcgrdcrnrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'src_grid_corner_lon', 
     &                                   nc_srcgrdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'src_grid_area', 
     &                                   nc_srcgrdarea_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'src_grid_frac', 
     &                                   nc_srcgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_dims', 
     &                      nc_dstgrddims_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_imask', 
     &                      nc_dstgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_center_lat', 
     &                                   nc_dstgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_center_lon', 
     &                                   nc_dstgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_corner_lat', 
     &                                   nc_dstgrdcrnrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_corner_lon', 
     &                                   nc_dstgrdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_area', 
     &                                   nc_dstgrdarea_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_frac', 
     &                                   nc_dstgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'src_address', 
     &                                   nc_srcgrdadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_address', 
     &                                   nc_dstgrdadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'remap_matrix', 
     &                                   nc_rmpmatrix_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     read all variables
!
!-----------------------------------------------------------------------

      ncstat = nf_get_var_int(nc_file_id, nc_srcgrddims_id, 
     &                        grid1_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_srcgrdimask_id, 
     &                        grid1_mask_int)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdcntrlat_id, 
     &                                       grid1_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdcntrlon_id, 
     &                                       grid1_center_lon)
      call netcdf_error_handler(ncstat)

      grid1_units = ' '
      ncstat = nf_get_att_text(nc_file_id, nc_srcgrdcntrlat_id, 'units',
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

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdcrnrlat_id, 
     &                                       grid1_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdcrnrlon_id, 
     &                                       grid1_corner_lon)
      call netcdf_error_handler(ncstat)

      grid1_units = ' '
      ncstat = nf_get_att_text(nc_file_id, nc_srcgrdcrnrlat_id, 'units',
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

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdarea_id, 
     &                                       grid1_area)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdfrac_id, 
     &                                       grid1_frac)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_dstgrddims_id, 
     &                        grid2_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_dstgrdimask_id, 
     &                        grid2_mask_int)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdcntrlat_id, 
     &                                       grid2_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdcntrlon_id, 
     &                                       grid2_center_lon)
      call netcdf_error_handler(ncstat)

      grid2_units = ' '
      ncstat = nf_get_att_text(nc_file_id, nc_dstgrdcntrlat_id, 'units',
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

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdcrnrlat_id, 
     &                                       grid2_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdcrnrlon_id, 
     &                                       grid2_corner_lon)
      call netcdf_error_handler(ncstat)

      grid2_units = ' '
      ncstat = nf_get_att_text(nc_file_id, nc_dstgrdcrnrlat_id, 'units',
     &                         grid2_units)
      call netcdf_error_handler(ncstat)

      select case (grid2_units(1:7))
      case ('degrees')
        grid2_corner_lat = grid2_corner_lat*deg2rad
        grid2_corner_lon = grid2_corner_lon*deg2rad
      case ('radians')
        !*** no conversion necessary
      case default
        print *,'unknown units supplied for grid2 corner lat/lon: '
        print *,'proceeding assuming radians'
      end select

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdarea_id, 
     &                                       grid2_area)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdfrac_id, 
     &                                       grid2_frac)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_srcgrdadd_id, 
     &                        grid1_add_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_dstgrdadd_id, 
     &                        grid2_add_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                                       wts_map1)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     initialize logical mask 
!
!-----------------------------------------------------------------------

      where (grid1_mask_int == 1)
        grid1_mask = .true.
      elsewhere
        grid1_mask = .false.
      endwhere
      where (grid2_mask_int == 1)
        grid2_mask = .true.
      elsewhere
        grid2_mask = .false.
      endwhere
      deallocate(grid1_mask_int, grid2_mask_int)

!-----------------------------------------------------------------------
!
!     close input file
!
!-----------------------------------------------------------------------

      ncstat = nf_close(nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------

      end subroutine read_remap_scrip

!***********************************************************************

      subroutine read_remap_csm

!-----------------------------------------------------------------------
!
!     the routine reads a netCDF file to extract remapping info
!     in NCAR-CSM format
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (char_len) ::
     &  grid1_name           ! grid name for source grid
     &, grid2_name           ! grid name for dest   grid

      integer (kind=int_kind) ::
     &  nc_numwgts1_id    ! extra netCDF id for num_wgts > 1 
     &, nc_rmpmatrix2_id  ! extra netCDF id for high-order remap matrix

      real (kind=dbl_kind), dimension(:),allocatable ::
     &  wts1              ! CSM wants single array for 1st-order wts

      real (kind=dbl_kind), dimension(:,:),allocatable ::
     &  wts2              ! write remaining weights in different array

      integer (kind=int_kind) ::  
     &  n                    ! dummy index

      integer (kind=int_kind), dimension(:), allocatable ::
     &  grid1_mask_int,      ! integer masks to determine
     &  grid2_mask_int       ! cells that participate in map

!-----------------------------------------------------------------------
!
!     read some additional global attributes
!
!-----------------------------------------------------------------------

      !***
      !*** source and destination grid names
      !***

      grid1_name = ' '
      grid2_name = ' '
      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'domain_a',
     &                          grid1_name)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'domain_b',
     &                          grid2_name)
      call netcdf_error_handler(ncstat)

      print *,' '
      print *,'Remapping between:',trim(grid1_name)
      print *,'and ',trim(grid2_name)
      print *,' '

!-----------------------------------------------------------------------
!
!     read dimension information
!
!-----------------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'n_a', nc_srcgrdsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdsize_id, grid1_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'n_b', nc_dstgrdsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdsize_id, grid2_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'nv_a', nc_srcgrdcorn_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdcorn_id, 
     &                       grid1_corners)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'nv_b', nc_dstgrdcorn_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdcorn_id, 
     &                       grid2_corners)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'src_grid_rank', 
     &                      nc_srcgrdrank_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdrank_id, 
     &                       grid1_rank)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'dst_grid_rank', 
     &                      nc_dstgrdrank_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdrank_id, 
     &                       grid2_rank)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'n_s', 
     &                      nc_numlinks_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numlinks_id, 
     &                       num_links_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'num_wgts', 
     &                      nc_numwgts_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numwgts_id, num_wts)
      call netcdf_error_handler(ncstat)

      if (num_wts > 1) then
        ncstat = nf_inq_dimid(nc_file_id, 'num_wgts1', 
     &                        nc_numwgts1_id)
        call netcdf_error_handler(ncstat)
      endif

!-----------------------------------------------------------------------
!
!     allocate arrays
!
!-----------------------------------------------------------------------

      allocate( grid1_dims      (grid1_rank),
     &          grid1_center_lat(grid1_size), 
     &          grid1_center_lon(grid1_size),
     &          grid1_area      (grid1_size),
     &          grid1_frac      (grid1_size),
     &          grid1_mask      (grid1_size),
     &          grid1_mask_int  (grid1_size),
     &          grid1_corner_lat(grid1_corners, grid1_size),
     &          grid1_corner_lon(grid1_corners, grid1_size) )

      allocate( grid2_dims      (grid2_rank),
     &          grid2_center_lat(grid2_size), 
     &          grid2_center_lon(grid2_size),
     &          grid2_area      (grid2_size),
     &          grid2_frac      (grid2_size),
     &          grid2_mask      (grid2_size),
     &          grid2_mask_int  (grid2_size),
     &          grid2_corner_lat(grid2_corners, grid2_size),
     &          grid2_corner_lon(grid2_corners, grid2_size) )

      allocate( grid1_add_map1(num_links_map1),
     &          grid2_add_map1(num_links_map1),
     &          wts_map1(num_wts,num_links_map1),
     &          wts1(num_links_map1),
     &          wts2(num_wts-1,num_links_map1) )

!-----------------------------------------------------------------------
!
!     get variable ids
!
!-----------------------------------------------------------------------

      ncstat = nf_inq_varid(nc_file_id, 'src_grid_dims', 
     &                      nc_srcgrddims_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'mask_a', 
     &                      nc_srcgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'yc_a', nc_srcgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'xc_a', nc_srcgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'yv_a', nc_srcgrdcrnrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'xv_a', nc_srcgrdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'area_a', nc_srcgrdarea_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'frac_a', nc_srcgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_dims', 
     &                      nc_dstgrddims_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'mask_b', 
     &                      nc_dstgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'yc_b', nc_dstgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'xc_b', nc_dstgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'yv_b', nc_dstgrdcrnrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'xv_b', nc_dstgrdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'area_b', nc_dstgrdarea_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'frac_b', nc_dstgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'col', nc_srcgrdadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'row', nc_dstgrdadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'S', nc_rmpmatrix_id)
      call netcdf_error_handler(ncstat)

      if (num_wts > 1) then
        ncstat = nf_inq_varid(nc_file_id, 'S2', nc_rmpmatrix2_id)
        call netcdf_error_handler(ncstat)
      endif

!-----------------------------------------------------------------------
!
!     read all variables
!
!-----------------------------------------------------------------------

      ncstat = nf_get_var_int(nc_file_id, nc_srcgrddims_id, 
     &                        grid1_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_srcgrdimask_id, 
     &                        grid1_mask_int)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdcntrlat_id, 
     &                                       grid1_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdcntrlon_id, 
     &                                       grid1_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_att_text(nc_file_id, nc_srcgrdcntrlat_id, 'units',
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

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdcrnrlat_id, 
     &                                       grid1_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdcrnrlon_id, 
     &                                       grid1_corner_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_att_text(nc_file_id, nc_srcgrdcrnrlat_id, 'units',
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

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdarea_id, 
     &                                       grid1_area)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_srcgrdfrac_id, 
     &                                       grid1_frac)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_dstgrddims_id, 
     &                        grid2_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_dstgrdimask_id, 
     &                        grid2_mask_int)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdcntrlat_id, 
     &                                       grid2_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdcntrlon_id, 
     &                                       grid2_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_att_text(nc_file_id, nc_dstgrdcntrlat_id, 'units',
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

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdcrnrlat_id, 
     &                                       grid2_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdcrnrlon_id, 
     &                                       grid2_corner_lon)
      call netcdf_error_handler(ncstat)


      ncstat = nf_get_att_text(nc_file_id, nc_dstgrdcrnrlat_id, 'units',
     &                         grid2_units)
      call netcdf_error_handler(ncstat)

      select case (grid2_units(1:7))
      case ('degrees')
        grid2_corner_lat = grid2_corner_lat*deg2rad
        grid2_corner_lon = grid2_corner_lon*deg2rad
      case ('radians')
        !*** no conversion necessary
      case default
        print *,'unknown units supplied for grid2 corner lat/lon: '
        print *,'proceeding assuming radians'
      end select

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdarea_id, 
     &                                       grid2_area)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdfrac_id, 
     &                                       grid2_frac)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_srcgrdadd_id, 
     &                        grid1_add_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_dstgrdadd_id, 
     &                        grid2_add_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                                       wts1)
      wts_map1(1,:) = wts1
      deallocate(wts1)

      if (num_wts > 1) then
        ncstat = nf_get_var_double(nc_file_id, nc_rmpmatrix2_id, 
     &                                         wts2)
        wts_map1(2:,:) = wts2
        deallocate(wts2)
      endif
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     initialize logical mask 
!
!-----------------------------------------------------------------------

      where (grid1_mask_int == 1)
        grid1_mask = .true.
      elsewhere
        grid1_mask = .false.
      endwhere
      where (grid2_mask_int == 1)
        grid2_mask = .true.
      elsewhere
        grid2_mask = .false.
      endwhere
      deallocate(grid1_mask_int, grid2_mask_int)

!-----------------------------------------------------------------------
!
!     close input file
!
!-----------------------------------------------------------------------

      ncstat = nf_close(nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------

      end subroutine read_remap_csm

!***********************************************************************

      end module remap_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
