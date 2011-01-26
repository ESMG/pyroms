!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module contains routines for writing the remapping data to 
!     a file.  Before writing the data for each mapping, the links are 
!     sorted by destination grid address.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_write.f,v 1.7 2001/08/21 21:06:42 pwjones Exp $
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

      module remap_write

!-----------------------------------------------------------------------

      use kinds_mod     ! defines common data types
      use constants     ! defines common scalar constants
      use grids         ! module containing grid information
      use remap_vars    ! module containing remap information
      use netcdf_mod    ! module with netCDF stuff

      implicit none

!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

      character(char_len), private :: 
     &   map_method       ! character string for map_type
     &,  normalize_opt    ! character string for normalization option
     &,  history          ! character string for history information
     &,  convention       ! character string for output convention

      character(8), private :: 
     &   cdate            ! character date string

      integer (kind=int_kind), dimension(:), allocatable, private ::
     &   src_mask_int     ! integer masks to determine
     &,  dst_mask_int     ! cells that participate in map

!-----------------------------------------------------------------------
!
!     various netCDF identifiers used by output routines
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), private ::
     &   ncstat               ! error flag for netCDF calls 
     &,  nc_file_id           ! id for netCDF file
     &,  nc_srcgrdsize_id     ! id for source grid size
     &,  nc_dstgrdsize_id     ! id for destination grid size
     &,  nc_srcgrdcorn_id     ! id for number of source grid corners
     &,  nc_dstgrdcorn_id     ! id for number of dest grid corners
     &,  nc_srcgrdrank_id     ! id for source grid rank
     &,  nc_dstgrdrank_id     ! id for dest grid rank
     &,  nc_numlinks_id       ! id for number of links in mapping
     &,  nc_numwgts_id        ! id for number of weights for mapping
     &,  nc_srcgrddims_id     ! id for source grid dimensions
     &,  nc_dstgrddims_id     ! id for dest grid dimensions
     &,  nc_srcgrdcntrlat_id  ! id for source grid center latitude
     &,  nc_dstgrdcntrlat_id  ! id for dest grid center latitude
     &,  nc_srcgrdcntrlon_id  ! id for source grid center longitude
     &,  nc_dstgrdcntrlon_id  ! id for dest grid center longitude
     &,  nc_srcgrdimask_id    ! id for source grid mask
     &,  nc_dstgrdimask_id    ! id for dest grid mask
     &,  nc_srcgrdcrnrlat_id  ! id for latitude of source grid corners
     &,  nc_srcgrdcrnrlon_id  ! id for longitude of source grid corners
     &,  nc_dstgrdcrnrlat_id  ! id for latitude of dest grid corners
     &,  nc_dstgrdcrnrlon_id  ! id for longitude of dest grid corners
     &,  nc_srcgrdarea_id     ! id for area of source grid cells
     &,  nc_dstgrdarea_id     ! id for area of dest grid cells
     &,  nc_srcgrdfrac_id     ! id for area fraction on source grid
     &,  nc_dstgrdfrac_id     ! id for area fraction on dest grid
     &,  nc_srcadd_id         ! id for map source address
     &,  nc_dstadd_id         ! id for map destination address
     &,  nc_rmpmatrix_id      ! id for remapping matrix

      integer (kind=int_kind), dimension(2), private ::
     &   nc_dims2_id  ! netCDF ids for 2d array dims

!***********************************************************************

      contains

!***********************************************************************

      subroutine write_remap(map1_name, map2_name, 
     &                       interp_file1, interp_file2, output_opt)

!-----------------------------------------------------------------------
!
!     calls correct output routine based on output format choice
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(char_len), intent(in) ::
     &            map1_name,    ! name for mapping grid1 to grid2
     &            map2_name,    ! name for mapping grid2 to grid1
     &            interp_file1, ! filename for map1 remap data
     &            interp_file2, ! filename for map2 remap data
     &            output_opt    ! option for output conventions

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     define some common variables to be used in all routines
!
!-----------------------------------------------------------------------

      select case(norm_opt)
      case (norm_opt_none)
        normalize_opt = 'none'
      case (norm_opt_frcarea)
        normalize_opt = 'fracarea'
      case (norm_opt_dstarea)
        normalize_opt = 'destarea'
      end select

      select case(map_type)
      case(map_type_conserv)
        map_method = 'Conservative remapping'
      case(map_type_bilinear)
        map_method = 'Bilinear remapping'
      case(map_type_distwgt)
        map_method = 'Distance weighted avg of nearest neighbors'
      case(map_type_bicubic)
        map_method = 'Bicubic remapping'
      case default
        stop 'Invalid Map Type'
      end select

      call date_and_time(date=cdate)
      write (history,1000) cdate(5:6),cdate(7:8),cdate(1:4)
 1000 format('Created: ',a2,'-',a2,'-',a4)

!-----------------------------------------------------------------------
!
!     sort address and weight arrays
!
!-----------------------------------------------------------------------

      call sort_add(grid2_add_map1, grid1_add_map1, wts_map1)
      if (num_maps > 1) then
        call sort_add(grid1_add_map2, grid2_add_map2, wts_map2)
      endif

!-----------------------------------------------------------------------
!
!     call appropriate output routine
!
!-----------------------------------------------------------------------

      select case(output_opt)
      case ('scrip')
        call write_remap_scrip(map1_name, interp_file1, 1)
      case ('ncar-csm')
        call write_remap_csm  (map1_name, interp_file1, 1)
      case default
        stop 'unknown output file convention'
      end select

!-----------------------------------------------------------------------
!
!     call appropriate output routine for second mapping if required
!
!-----------------------------------------------------------------------

      if (num_maps > 1) then
        select case(output_opt)
        case ('scrip')
          call write_remap_scrip(map2_name, interp_file2, 2)
        case ('ncar-csm')
          call write_remap_csm  (map2_name, interp_file2, 2)
        case default
          stop 'unknown output file convention'
        end select
      endif

!-----------------------------------------------------------------------

      end subroutine write_remap

!***********************************************************************

      subroutine write_remap_scrip(map_name, interp_file, direction)

!-----------------------------------------------------------------------
!
!     writes remap data to a netCDF file using SCRIP conventions
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(char_len), intent(in) ::
     &            map_name     ! name for mapping 
     &,           interp_file  ! filename for remap data

      integer (kind=int_kind), intent(in) ::
     &  direction              ! direction of map (1=grid1 to grid2
                               !                   2=grid2 to grid1)

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character(char_len) ::
     &  grid1_ctmp        ! character temp for grid1 names
     &, grid2_ctmp        ! character temp for grid2 names

      integer (kind=int_kind) ::
     &  itmp1             ! integer temp
     &, itmp2             ! integer temp
     &, itmp3             ! integer temp
     &, itmp4             ! integer temp

!-----------------------------------------------------------------------
!
!     create netCDF file for mapping and define some global attributes
!
!-----------------------------------------------------------------------

      ncstat = nf_create (interp_file, NF_CLOBBER, nc_file_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** map name
      !***
      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, 'title',
     &                          len_trim(map_name), map_name)
      call netcdf_error_handler(ncstat)

      !***
      !*** normalization option
      !***
      ncstat = nf_put_att_text(nc_file_id, NF_GLOBAL, 'normalization',
     &                         len_trim(normalize_opt), normalize_opt)
      call netcdf_error_handler(ncstat)

      !***
      !*** map method
      !***
      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, 'map_method',
     &                          len_trim(map_method), map_method)
      call netcdf_error_handler(ncstat)

      !***
      !*** history
      !***
      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, 'history',
     &                          len_trim(history), history)
      call netcdf_error_handler(ncstat)

      !***
      !*** file convention
      !***
      convention = 'SCRIP'
      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, 'conventions',
     &                          len_trim(convention), convention)
      call netcdf_error_handler(ncstat)

      !***
      !*** source and destination grid names
      !***

      if (direction == 1) then
        grid1_ctmp = 'source_grid'
        grid2_ctmp = 'dest_grid'
      else
        grid1_ctmp = 'dest_grid'
        grid2_ctmp = 'source_grid'
      endif

      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, trim(grid1_ctmp),
     &                          len_trim(grid1_name), grid1_name)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, trim(grid2_ctmp),
     &                          len_trim(grid2_name), grid2_name)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     prepare netCDF dimension info
!
!-----------------------------------------------------------------------

      !***
      !*** define grid size dimensions
      !***

      if (direction == 1) then
        itmp1 = grid1_size
        itmp2 = grid2_size
      else
        itmp1 = grid2_size
        itmp2 = grid1_size
      endif

      ncstat = nf_def_dim (nc_file_id, 'src_grid_size', itmp1, 
     &                     nc_srcgrdsize_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'dst_grid_size', itmp2, 
     &                     nc_dstgrdsize_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid corner dimension
      !***

      if (direction == 1) then
        itmp1 = grid1_corners
        itmp2 = grid2_corners
      else
        itmp1 = grid2_corners
        itmp2 = grid1_corners
      endif

      ncstat = nf_def_dim (nc_file_id, 'src_grid_corners', 
     &                     itmp1, nc_srcgrdcorn_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'dst_grid_corners', 
     &                     itmp2, nc_dstgrdcorn_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid rank dimension
      !***

      if (direction == 1) then
        itmp1 = grid1_rank
        itmp2 = grid2_rank
      else
        itmp1 = grid2_rank
        itmp2 = grid1_rank
      endif

      ncstat = nf_def_dim (nc_file_id, 'src_grid_rank', 
     &                     itmp1, nc_srcgrdrank_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'dst_grid_rank', 
     &                     itmp2, nc_dstgrdrank_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define map size dimensions
      !***

      if (direction == 1) then
        itmp1 = num_links_map1
      else
        itmp1 = num_links_map2
      endif

      ncstat = nf_def_dim (nc_file_id, 'num_links', 
     &                     itmp1, nc_numlinks_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'num_wgts', 
     &                     num_wts, nc_numwgts_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid dimensions
      !***

      ncstat = nf_def_var (nc_file_id, 'src_grid_dims', NF_INT,
     &                     1, nc_srcgrdrank_id, nc_srcgrddims_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'dst_grid_dims', NF_INT,
     &                     1, nc_dstgrdrank_id, nc_dstgrddims_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     define all arrays for netCDF descriptors
!
!-----------------------------------------------------------------------

      !***
      !*** define grid center latitude array
      !***

      ncstat = nf_def_var (nc_file_id, 'src_grid_center_lat', 
     &                     NF_DOUBLE, 1, nc_srcgrdsize_id, 
     &                     nc_srcgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'dst_grid_center_lat', 
     &                     NF_DOUBLE, 1, nc_dstgrdsize_id, 
     &                     nc_dstgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid center longitude array
      !***

      ncstat = nf_def_var (nc_file_id, 'src_grid_center_lon', 
     &                     NF_DOUBLE, 1, nc_srcgrdsize_id, 
     &                     nc_srcgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'dst_grid_center_lon', 
     &                     NF_DOUBLE, 1, nc_dstgrdsize_id, 
     &                     nc_dstgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid corner lat/lon arrays
      !***

      nc_dims2_id(1) = nc_srcgrdcorn_id
      nc_dims2_id(2) = nc_srcgrdsize_id

      ncstat = nf_def_var (nc_file_id, 'src_grid_corner_lat', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_srcgrdcrnrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'src_grid_corner_lon', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_srcgrdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      nc_dims2_id(1) = nc_dstgrdcorn_id
      nc_dims2_id(2) = nc_dstgrdsize_id

      ncstat = nf_def_var (nc_file_id, 'dst_grid_corner_lat', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_dstgrdcrnrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'dst_grid_corner_lon', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_dstgrdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define units for all coordinate arrays
      !***

      if (direction == 1) then
        grid1_ctmp = grid1_units
        grid2_ctmp = grid2_units
      else
        grid1_ctmp = grid2_units
        grid2_ctmp = grid1_units
      endif

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdcntrlat_id, 
     &                          'units', 7, grid1_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdcntrlat_id, 
     &                          'units', 7, grid2_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdcntrlon_id, 
     &                          'units', 7, grid1_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdcntrlon_id, 
     &                          'units', 7, grid2_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdcrnrlat_id, 
     &                          'units', 7, grid1_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdcrnrlon_id, 
     &                          'units', 7, grid1_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdcrnrlat_id, 
     &                          'units', 7, grid2_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdcrnrlon_id, 
     &                          'units', 7, grid2_ctmp)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid mask
      !***

      ncstat = nf_def_var (nc_file_id, 'src_grid_imask', NF_INT,
     &                     1, nc_srcgrdsize_id, nc_srcgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdimask_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'dst_grid_imask', NF_INT,
     &                     1, nc_dstgrdsize_id, nc_dstgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdimask_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid area arrays
      !***

      ncstat = nf_def_var (nc_file_id, 'src_grid_area', 
     &                     NF_DOUBLE, 1, nc_srcgrdsize_id, 
     &                     nc_srcgrdarea_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdarea_id, 
     &                          'units', 14, 'square radians')
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'dst_grid_area', 
     &                     NF_DOUBLE, 1, nc_dstgrdsize_id, 
     &                     nc_dstgrdarea_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdarea_id, 
     &                          'units', 14, 'square radians')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid fraction arrays
      !***

      ncstat = nf_def_var (nc_file_id, 'src_grid_frac', 
     &                     NF_DOUBLE, 1, nc_srcgrdsize_id, 
     &                     nc_srcgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdfrac_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'dst_grid_frac', 
     &                     NF_DOUBLE, 1, nc_dstgrdsize_id, 
     &                     nc_dstgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdfrac_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      !***
      !*** define mapping arrays
      !***

      ncstat = nf_def_var (nc_file_id, 'src_address', 
     &                     NF_INT, 1, nc_numlinks_id, 
     &                     nc_srcadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'dst_address', 
     &                     NF_INT, 1, nc_numlinks_id, 
     &                     nc_dstadd_id)
      call netcdf_error_handler(ncstat)

      nc_dims2_id(1) = nc_numwgts_id
      nc_dims2_id(2) = nc_numlinks_id

      ncstat = nf_def_var (nc_file_id, 'remap_matrix', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_rmpmatrix_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** end definition stage
      !***

      ncstat = nf_enddef(nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     compute integer masks
!
!-----------------------------------------------------------------------

      if (direction == 1) then
        allocate (src_mask_int(grid1_size),
     &            dst_mask_int(grid2_size))

        where (grid2_mask)
          dst_mask_int = 1
        elsewhere
          dst_mask_int = 0
        endwhere

        where (grid1_mask)
          src_mask_int = 1
        elsewhere
          src_mask_int = 0
        endwhere
      else
        allocate (src_mask_int(grid2_size),
     &            dst_mask_int(grid1_size))

        where (grid1_mask)
          dst_mask_int = 1
        elsewhere
          dst_mask_int = 0
        endwhere

        where (grid2_mask)
          src_mask_int = 1
        elsewhere
          src_mask_int = 0
        endwhere
      endif

!-----------------------------------------------------------------------
!
!     change units of lat/lon coordinates if input units different
!     from radians
!
!-----------------------------------------------------------------------

      if (grid1_units(1:7) == 'degrees' .and. direction == 1) then
        grid1_center_lat = grid1_center_lat/deg2rad
        grid1_center_lon = grid1_center_lon/deg2rad
        grid1_corner_lat = grid1_corner_lat/deg2rad
        grid1_corner_lon = grid1_corner_lon/deg2rad
      endif

      if (grid2_units(1:7) == 'degrees' .and. direction == 1) then
        grid2_center_lat = grid2_center_lat/deg2rad
        grid2_center_lon = grid2_center_lon/deg2rad
        grid2_corner_lat = grid2_corner_lat/deg2rad
        grid2_corner_lon = grid2_corner_lon/deg2rad
      endif

!-----------------------------------------------------------------------
!
!     write mapping data
!
!-----------------------------------------------------------------------

      if (direction == 1) then
        itmp1 = nc_srcgrddims_id
        itmp2 = nc_dstgrddims_id
      else
        itmp2 = nc_srcgrddims_id
        itmp1 = nc_dstgrddims_id
      endif

      ncstat = nf_put_var_int(nc_file_id, itmp1, grid1_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_int(nc_file_id, itmp2, grid2_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_int(nc_file_id, nc_srcgrdimask_id, 
     &                        src_mask_int)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_int(nc_file_id, nc_dstgrdimask_id,
     &                        dst_mask_int)
      call netcdf_error_handler(ncstat)

      deallocate(src_mask_int, dst_mask_int)

      if (direction == 1) then
        itmp1 = nc_srcgrdcntrlat_id
        itmp2 = nc_srcgrdcntrlon_id
        itmp3 = nc_srcgrdcrnrlat_id
        itmp4 = nc_srcgrdcrnrlon_id
      else
        itmp1 = nc_dstgrdcntrlat_id
        itmp2 = nc_dstgrdcntrlon_id
        itmp3 = nc_dstgrdcrnrlat_id
        itmp4 = nc_dstgrdcrnrlon_id
      endif

      ncstat = nf_put_var_double(nc_file_id, itmp1, grid1_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp2, grid1_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp3, grid1_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp4, grid1_corner_lon)
      call netcdf_error_handler(ncstat)

      if (direction == 1) then
        itmp1 = nc_dstgrdcntrlat_id
        itmp2 = nc_dstgrdcntrlon_id
        itmp3 = nc_dstgrdcrnrlat_id
        itmp4 = nc_dstgrdcrnrlon_id
      else
        itmp1 = nc_srcgrdcntrlat_id
        itmp2 = nc_srcgrdcntrlon_id
        itmp3 = nc_srcgrdcrnrlat_id
        itmp4 = nc_srcgrdcrnrlon_id
      endif

      ncstat = nf_put_var_double(nc_file_id, itmp1, grid2_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp2, grid2_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp3, grid2_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp4, grid2_corner_lon)
      call netcdf_error_handler(ncstat)

      if (direction == 1) then
        itmp1 = nc_srcgrdarea_id
        itmp2 = nc_srcgrdfrac_id
        itmp3 = nc_dstgrdarea_id
        itmp4 = nc_dstgrdfrac_id
      else
        itmp1 = nc_dstgrdarea_id
        itmp2 = nc_dstgrdfrac_id
        itmp3 = nc_srcgrdarea_id
        itmp4 = nc_srcgrdfrac_id
      endif

      if (luse_grid1_area) then
        ncstat = nf_put_var_double(nc_file_id, itmp1, grid1_area_in)
      else
        ncstat = nf_put_var_double(nc_file_id, itmp1, grid1_area)
      endif
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp2, grid1_frac)
      call netcdf_error_handler(ncstat)

      if (luse_grid2_area) then
        ncstat = nf_put_var_double(nc_file_id, itmp3, grid2_area_in)
      else
        ncstat = nf_put_var_double(nc_file_id, itmp3, grid2_area)
      endif
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp4, grid2_frac)
      call netcdf_error_handler(ncstat)

      if (direction == 1) then
        ncstat = nf_put_var_int(nc_file_id, nc_srcadd_id, 
     &                          grid1_add_map1)
        call netcdf_error_handler(ncstat)

        ncstat = nf_put_var_int(nc_file_id, nc_dstadd_id, 
     &                          grid2_add_map1)
        call netcdf_error_handler(ncstat)

        ncstat = nf_put_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                             wts_map1)
        call netcdf_error_handler(ncstat)
      else
        ncstat = nf_put_var_int(nc_file_id, nc_srcadd_id, 
     &                          grid2_add_map2)
        call netcdf_error_handler(ncstat)

        ncstat = nf_put_var_int(nc_file_id, nc_dstadd_id, 
     &                          grid1_add_map2)
        call netcdf_error_handler(ncstat)

        ncstat = nf_put_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                             wts_map2)
        call netcdf_error_handler(ncstat)
      endif

      ncstat = nf_close(nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------

      end subroutine write_remap_scrip

!***********************************************************************

      subroutine write_remap_csm(map_name, interp_file, direction)

!-----------------------------------------------------------------------
!
!     writes remap data to a netCDF file using NCAR-CSM conventions
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(char_len), intent(in) ::
     &            map_name     ! name for mapping 
     &,           interp_file  ! filename for remap data

      integer (kind=int_kind), intent(in) ::
     &  direction              ! direction of map (1=grid1 to grid2
                               !                   2=grid2 to grid1)

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character(char_len) ::
     &  grid1_ctmp        ! character temp for grid1 names
     &, grid2_ctmp        ! character temp for grid2 names

      integer (kind=int_kind) ::
     &  itmp1             ! integer temp
     &, itmp2             ! integer temp
     &, itmp3             ! integer temp
     &, itmp4             ! integer temp
     &, nc_numwgts1_id    ! extra netCDF id for additional weights
     &, nc_src_isize_id   ! extra netCDF id for ni_a
     &, nc_src_jsize_id   ! extra netCDF id for nj_a
     &, nc_dst_isize_id   ! extra netCDF id for ni_b
     &, nc_dst_jsize_id   ! extra netCDF id for nj_b
     &, nc_rmpmatrix2_id  ! extra netCDF id for high-order remap matrix

      real (kind=dbl_kind), dimension(:),allocatable ::
     &  wts1              ! CSM wants single array for 1st-order wts

      real (kind=dbl_kind), dimension(:,:),allocatable ::
     &  wts2              ! write remaining weights in different array

!-----------------------------------------------------------------------
!
!     create netCDF file for mapping and define some global attributes
!
!-----------------------------------------------------------------------

      ncstat = nf_create (interp_file, NF_CLOBBER, nc_file_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** map name
      !***
      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, 'title',
     &                          len_trim(map_name), map_name)
      call netcdf_error_handler(ncstat)

      !***
      !*** normalization option
      !***
      ncstat = nf_put_att_text(nc_file_id, NF_GLOBAL, 'normalization',
     &                         len_trim(normalize_opt), normalize_opt)
      call netcdf_error_handler(ncstat)

      !***
      !*** map method
      !***
      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, 'map_method',
     &                          len_trim(map_method), map_method)
      call netcdf_error_handler(ncstat)

      !***
      !*** history
      !***
      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, 'history',
     &                          len_trim(history), history)
      call netcdf_error_handler(ncstat)

      !***
      !*** file convention
      !***
      convention = 'NCAR-CSM'
      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, 'conventions',
     &                          len_trim(convention), convention)
      call netcdf_error_handler(ncstat)

      !***
      !*** source and destination grid names
      !***

      if (direction == 1) then
        grid1_ctmp = 'domain_a'
        grid2_ctmp = 'domain_b'
      else
        grid1_ctmp = 'domain_b'
        grid2_ctmp = 'domain_a'
      endif

      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, trim(grid1_ctmp),
     &                          len_trim(grid1_name), grid1_name)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, NF_GLOBAL, trim(grid2_ctmp),
     &                          len_trim(grid2_name), grid2_name)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     prepare netCDF dimension info
!
!-----------------------------------------------------------------------

      !***
      !*** define grid size dimensions
      !***

      if (direction == 1) then
        itmp1 = grid1_size
        itmp2 = grid2_size
      else
        itmp1 = grid2_size
        itmp2 = grid1_size
      endif

      ncstat = nf_def_dim (nc_file_id, 'n_a', itmp1, nc_srcgrdsize_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'n_b', itmp2, nc_dstgrdsize_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid corner dimension
      !***

      if (direction == 1) then
        itmp1 = grid1_corners
        itmp2 = grid2_corners
      else
        itmp1 = grid2_corners
        itmp2 = grid1_corners
      endif

      ncstat = nf_def_dim (nc_file_id, 'nv_a', itmp1, nc_srcgrdcorn_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'nv_b', itmp2, nc_dstgrdcorn_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid rank dimension
      !***

      if (direction == 1) then
        itmp1 = grid1_rank
        itmp2 = grid2_rank
      else
        itmp1 = grid2_rank
        itmp2 = grid1_rank
      endif

      ncstat = nf_def_dim (nc_file_id, 'src_grid_rank', 
     &                     itmp1, nc_srcgrdrank_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'dst_grid_rank', 
     &                     itmp2, nc_dstgrdrank_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define first two dims as if 2-d cartesian domain
      !***

      if (direction == 1) then
        itmp1 = grid1_dims(1)
        if (grid1_rank > 1) then
          itmp2 = grid1_dims(2)
        else
          itmp2 = 0
        endif
        itmp3 = grid2_dims(1)
        if (grid2_rank > 1) then
          itmp4 = grid2_dims(2)
        else
          itmp4 = 0
        endif
      else
        itmp1 = grid2_dims(1)
        if (grid2_rank > 1) then
          itmp2 = grid2_dims(2)
        else
          itmp2 = 0
        endif
        itmp3 = grid1_dims(1)
        if (grid1_rank > 1) then
          itmp4 = grid1_dims(2)
        else
          itmp4 = 0
        endif
      endif

      ncstat = nf_def_dim (nc_file_id, 'ni_a', itmp1, nc_src_isize_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'nj_a', itmp2, nc_src_jsize_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'ni_b', itmp3, nc_dst_isize_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'nj_b', itmp4, nc_dst_jsize_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define map size dimensions
      !***

      if (direction == 1) then
        itmp1 = num_links_map1
      else
        itmp1 = num_links_map2
      endif

      ncstat = nf_def_dim (nc_file_id, 'n_s', itmp1, nc_numlinks_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_dim (nc_file_id, 'num_wgts', 
     &                     num_wts, nc_numwgts_id)
      call netcdf_error_handler(ncstat)

      if (num_wts > 1) then
        ncstat = nf_def_dim (nc_file_id, 'num_wgts1', 
     &                       num_wts-1, nc_numwgts1_id)
        call netcdf_error_handler(ncstat)
      endif

      !***
      !*** define grid dimensions
      !***

      ncstat = nf_def_var (nc_file_id, 'src_grid_dims', NF_INT,
     &                     1, nc_srcgrdrank_id, nc_srcgrddims_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'dst_grid_dims', NF_INT,
     &                     1, nc_dstgrdrank_id, nc_dstgrddims_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     define all arrays for netCDF descriptors
!
!-----------------------------------------------------------------------

      !***
      !*** define grid center latitude array
      !***

      ncstat = nf_def_var (nc_file_id, 'yc_a',
     &                     NF_DOUBLE, 1, nc_srcgrdsize_id, 
     &                     nc_srcgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'yc_b', 
     &                     NF_DOUBLE, 1, nc_dstgrdsize_id, 
     &                     nc_dstgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid center longitude array
      !***

      ncstat = nf_def_var (nc_file_id, 'xc_a', 
     &                     NF_DOUBLE, 1, nc_srcgrdsize_id, 
     &                     nc_srcgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'xc_b', 
     &                     NF_DOUBLE, 1, nc_dstgrdsize_id, 
     &                     nc_dstgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid corner lat/lon arrays
      !***

      nc_dims2_id(1) = nc_srcgrdcorn_id
      nc_dims2_id(2) = nc_srcgrdsize_id

      ncstat = nf_def_var (nc_file_id, 'yv_a', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_srcgrdcrnrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'xv_a', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_srcgrdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      nc_dims2_id(1) = nc_dstgrdcorn_id
      nc_dims2_id(2) = nc_dstgrdsize_id

      ncstat = nf_def_var (nc_file_id, 'yv_b', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_dstgrdcrnrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'xv_b', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_dstgrdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** CSM wants all in degrees
      !***

      grid1_units = 'degrees'
      grid2_units = 'degrees'

      if (direction == 1) then
        grid1_ctmp = grid1_units
        grid2_ctmp = grid2_units
      else
        grid1_ctmp = grid2_units
        grid2_ctmp = grid1_units
      endif

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdcntrlat_id, 
     &                          'units', 7, grid1_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdcntrlat_id, 
     &                          'units', 7, grid2_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdcntrlon_id, 
     &                          'units', 7, grid1_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdcntrlon_id, 
     &                          'units', 7, grid2_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdcrnrlat_id, 
     &                          'units', 7, grid1_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdcrnrlon_id, 
     &                          'units', 7, grid1_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdcrnrlat_id, 
     &                          'units', 7, grid2_ctmp)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdcrnrlon_id, 
     &                          'units', 7, grid2_ctmp)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid mask
      !***

      ncstat = nf_def_var (nc_file_id, 'mask_a', NF_INT,
     &                     1, nc_srcgrdsize_id, nc_srcgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdimask_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'mask_b', NF_INT,
     &                     1, nc_dstgrdsize_id, nc_dstgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdimask_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid area arrays
      !***

      ncstat = nf_def_var (nc_file_id, 'area_a', 
     &                     NF_DOUBLE, 1, nc_srcgrdsize_id, 
     &                     nc_srcgrdarea_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdarea_id, 
     &                          'units', 14, 'square radians')
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'area_b', 
     &                     NF_DOUBLE, 1, nc_dstgrdsize_id, 
     &                     nc_dstgrdarea_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdarea_id, 
     &                          'units', 14, 'square radians')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid fraction arrays
      !***

      ncstat = nf_def_var (nc_file_id, 'frac_a', 
     &                     NF_DOUBLE, 1, nc_srcgrdsize_id, 
     &                     nc_srcgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_srcgrdfrac_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'frac_b', 
     &                     NF_DOUBLE, 1, nc_dstgrdsize_id, 
     &                     nc_dstgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_file_id, nc_dstgrdfrac_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      !***
      !*** define mapping arrays
      !***

      ncstat = nf_def_var (nc_file_id, 'col', 
     &                     NF_INT, 1, nc_numlinks_id, 
     &                     nc_srcadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'row', 
     &                     NF_INT, 1, nc_numlinks_id, 
     &                     nc_dstadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_file_id, 'S', 
     &                     NF_DOUBLE, 1, nc_numlinks_id, 
     &                     nc_rmpmatrix_id)
      call netcdf_error_handler(ncstat)

      if (num_wts > 1) then
        nc_dims2_id(1) = nc_numwgts1_id
        nc_dims2_id(2) = nc_numlinks_id

        ncstat = nf_def_var (nc_file_id, 'S2', 
     &                     NF_DOUBLE, 2, nc_dims2_id, 
     &                     nc_rmpmatrix2_id)
        call netcdf_error_handler(ncstat)
      endif

      !***
      !*** end definition stage
      !***

      ncstat = nf_enddef(nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     compute integer masks
!
!-----------------------------------------------------------------------

      if (direction == 1) then
        allocate (src_mask_int(grid1_size),
     &            dst_mask_int(grid2_size))

        where (grid2_mask)
          dst_mask_int = 1
        elsewhere
          dst_mask_int = 0
        endwhere

        where (grid1_mask)
          src_mask_int = 1
        elsewhere
          src_mask_int = 0
        endwhere
      else
        allocate (src_mask_int(grid2_size),
     &            dst_mask_int(grid1_size))

        where (grid1_mask)
          dst_mask_int = 1
        elsewhere
          dst_mask_int = 0
        endwhere

        where (grid2_mask)
          src_mask_int = 1
        elsewhere
          src_mask_int = 0
        endwhere
      endif

!-----------------------------------------------------------------------
!
!     change units of lat/lon coordinates if input units different
!     from radians. if this is the second mapping, the conversion has
!     alread been done.
!
!-----------------------------------------------------------------------

      if (grid1_units(1:7) == 'degrees' .and. direction == 1) then
        grid1_center_lat = grid1_center_lat/deg2rad
        grid1_center_lon = grid1_center_lon/deg2rad
        grid1_corner_lat = grid1_corner_lat/deg2rad
        grid1_corner_lon = grid1_corner_lon/deg2rad
      endif

      if (grid2_units(1:7) == 'degrees' .and. direction == 1) then
        grid2_center_lat = grid2_center_lat/deg2rad
        grid2_center_lon = grid2_center_lon/deg2rad
        grid2_corner_lat = grid2_corner_lat/deg2rad
        grid2_corner_lon = grid2_corner_lon/deg2rad
      endif

!-----------------------------------------------------------------------
!
!     write mapping data
!
!-----------------------------------------------------------------------

      if (direction == 1) then
        itmp1 = nc_srcgrddims_id
        itmp2 = nc_dstgrddims_id
      else
        itmp2 = nc_srcgrddims_id
        itmp1 = nc_dstgrddims_id
      endif

      ncstat = nf_put_var_int(nc_file_id, itmp1, grid1_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_int(nc_file_id, itmp2, grid2_dims)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_int(nc_file_id, nc_srcgrdimask_id, 
     &                        src_mask_int)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_int(nc_file_id, nc_dstgrdimask_id,
     &                        dst_mask_int)
      call netcdf_error_handler(ncstat)

      deallocate(src_mask_int, dst_mask_int)

      if (direction == 1) then
        itmp1 = nc_srcgrdcntrlat_id
        itmp2 = nc_srcgrdcntrlon_id
        itmp3 = nc_srcgrdcrnrlat_id
        itmp4 = nc_srcgrdcrnrlon_id
      else
        itmp1 = nc_dstgrdcntrlat_id
        itmp2 = nc_dstgrdcntrlon_id
        itmp3 = nc_dstgrdcrnrlat_id
        itmp4 = nc_dstgrdcrnrlon_id
      endif

      ncstat = nf_put_var_double(nc_file_id, itmp1, grid1_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp2, grid1_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp3, grid1_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp4, grid1_corner_lon)
      call netcdf_error_handler(ncstat)

      if (direction == 1) then
        itmp1 = nc_dstgrdcntrlat_id
        itmp2 = nc_dstgrdcntrlon_id
        itmp3 = nc_dstgrdcrnrlat_id
        itmp4 = nc_dstgrdcrnrlon_id
      else
        itmp1 = nc_srcgrdcntrlat_id
        itmp2 = nc_srcgrdcntrlon_id
        itmp3 = nc_srcgrdcrnrlat_id
        itmp4 = nc_srcgrdcrnrlon_id
      endif

      ncstat = nf_put_var_double(nc_file_id, itmp1, grid2_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp2, grid2_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp3, grid2_corner_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp4, grid2_corner_lon)
      call netcdf_error_handler(ncstat)

      if (direction == 1) then
        itmp1 = nc_srcgrdarea_id
        itmp2 = nc_srcgrdfrac_id
        itmp3 = nc_dstgrdarea_id
        itmp4 = nc_dstgrdfrac_id
      else
        itmp1 = nc_dstgrdarea_id
        itmp2 = nc_dstgrdfrac_id
        itmp3 = nc_srcgrdarea_id
        itmp4 = nc_srcgrdfrac_id
      endif

      if (luse_grid1_area) then
        ncstat = nf_put_var_double(nc_file_id, itmp1, grid1_area_in)
      else
        ncstat = nf_put_var_double(nc_file_id, itmp1, grid1_area)
      endif
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp2, grid1_frac)
      call netcdf_error_handler(ncstat)

      if (luse_grid2_area) then
        ncstat = nf_put_var_double(nc_file_id, itmp3, grid2_area)
      else
        ncstat = nf_put_var_double(nc_file_id, itmp3, grid2_area)
      endif
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_file_id, itmp4, grid2_frac)
      call netcdf_error_handler(ncstat)

      if (direction == 1) then
        ncstat = nf_put_var_int(nc_file_id, nc_srcadd_id, 
     &                          grid1_add_map1)
        call netcdf_error_handler(ncstat)

        ncstat = nf_put_var_int(nc_file_id, nc_dstadd_id, 
     &                          grid2_add_map1)
        call netcdf_error_handler(ncstat)

        if (num_wts == 1) then
          ncstat = nf_put_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                               wts_map1)
          call netcdf_error_handler(ncstat)
        else
          allocate(wts1(num_links_map1),wts2(num_wts-1,num_links_map1))

          wts1 = wts_map1(1,:)
          wts2 = wts_map1(2:,:)

          ncstat = nf_put_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                               wts1)
          call netcdf_error_handler(ncstat)
          ncstat = nf_put_var_double(nc_file_id, nc_rmpmatrix2_id, 
     &                               wts2)
          call netcdf_error_handler(ncstat)
          deallocate(wts1,wts2)
        endif
      else
        ncstat = nf_put_var_int(nc_file_id, nc_srcadd_id, 
     &                          grid2_add_map2)
        call netcdf_error_handler(ncstat)

        ncstat = nf_put_var_int(nc_file_id, nc_dstadd_id, 
     &                          grid1_add_map2)
        call netcdf_error_handler(ncstat)

        if (num_wts == 1) then
          ncstat = nf_put_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                               wts_map2)
          call netcdf_error_handler(ncstat)
        else
          allocate(wts1(num_links_map2),wts2(num_wts-1,num_links_map2))

          wts1 = wts_map2(1,:)
          wts2 = wts_map2(2:,:)

          ncstat = nf_put_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                               wts1)
          call netcdf_error_handler(ncstat)
          ncstat = nf_put_var_double(nc_file_id, nc_rmpmatrix2_id, 
     &                               wts2)
          call netcdf_error_handler(ncstat)
          deallocate(wts1,wts2)
        endif
      endif

      ncstat = nf_close(nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------

      end subroutine write_remap_csm

!***********************************************************************

      subroutine sort_add(add1, add2, weights)

!-----------------------------------------------------------------------
!
!     this routine sorts address and weight arrays based on the
!     destination address with the source address as a secondary
!     sorting criterion.  the method is a standard heap sort.
!
!-----------------------------------------------------------------------

      use kinds_mod     ! defines common data types
      use constants     ! defines common scalar constants

      implicit none

!-----------------------------------------------------------------------
!
!     Input and Output arrays
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(inout), dimension(:) ::
     &        add1,       ! destination address array (num_links)
     &        add2        ! source      address array

      real (kind=dbl_kind), intent(inout), dimension(:,:) ::
     &        weights     ! remapping weights (num_wts, num_links)

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) ::
     &          num_links,          ! num of links for this mapping
     &          num_wts,            ! num of weights for this mapping
     &          add1_tmp, add2_tmp, ! temp for addresses during swap
     &          nwgt,
     &          lvl, final_lvl,     ! level indexes for heap sort levels
     &          chk_lvl1, chk_lvl2, max_lvl

      real (kind=dbl_kind), dimension(SIZE(weights,DIM=1)) ::
     &          wgttmp              ! temp for holding wts during swap

!-----------------------------------------------------------------------
!
!     determine total number of links to sort and number of weights
!
!-----------------------------------------------------------------------

      num_links = SIZE(add1)
      num_wts   = SIZE(weights, DIM=1)

!-----------------------------------------------------------------------
!
!     start at the lowest level (N/2) of the tree and sift lower 
!     values to the bottom of the tree, promoting the larger numbers
!
!-----------------------------------------------------------------------

      do lvl=num_links/2,1,-1

        final_lvl = lvl
        add1_tmp = add1(lvl)
        add2_tmp = add2(lvl)
        wgttmp(:) = weights(:,lvl)

        !***
        !*** loop until proper level is found for this link, or reach
        !*** bottom
        !***

        sift_loop1: do

          !***
          !*** find the largest of the two daughters
          !***

          chk_lvl1 = 2*final_lvl
          chk_lvl2 = 2*final_lvl+1
          if (chk_lvl1 .EQ. num_links) chk_lvl2 = chk_lvl1

          if ((add1(chk_lvl1) >  add1(chk_lvl2)) .OR.
     &       ((add1(chk_lvl1) == add1(chk_lvl2)) .AND.
     &        (add2(chk_lvl1) >  add2(chk_lvl2)))) then
            max_lvl = chk_lvl1
          else 
            max_lvl = chk_lvl2
          endif

          !***
          !*** if the parent is greater than both daughters,
          !*** the correct level has been found
          !***

          if ((add1_tmp .GT. add1(max_lvl)) .OR.
     &       ((add1_tmp .EQ. add1(max_lvl)) .AND.
     &        (add2_tmp .GT. add2(max_lvl)))) then
            add1(final_lvl) = add1_tmp
            add2(final_lvl) = add2_tmp
            weights(:,final_lvl) = wgttmp(:)
            exit sift_loop1

          !***
          !*** otherwise, promote the largest daughter and push
          !*** down one level in the tree.  if haven't reached
          !*** the end of the tree, repeat the process.  otherwise
          !*** store last values and exit the loop
          !***

          else 
            add1(final_lvl) = add1(max_lvl)
            add2(final_lvl) = add2(max_lvl)
            weights(:,final_lvl) = weights(:,max_lvl)

            final_lvl = max_lvl
            if (2*final_lvl > num_links) then
              add1(final_lvl) = add1_tmp
              add2(final_lvl) = add2_tmp
              weights(:,final_lvl) = wgttmp(:)
              exit sift_loop1
            endif
          endif
        end do sift_loop1
      end do

!-----------------------------------------------------------------------
!
!     now that the heap has been sorted, strip off the top (largest)
!     value and promote the values below
!
!-----------------------------------------------------------------------

      do lvl=num_links,3,-1

        !***
        !*** move the top value and insert it into the correct place
        !***

        add1_tmp = add1(lvl)
        add1(lvl) = add1(1)

        add2_tmp = add2(lvl)
        add2(lvl) = add2(1)

        wgttmp(:) = weights(:,lvl)
        weights(:,lvl) = weights(:,1)

        !***
        !*** as above this loop sifts the tmp values down until proper 
        !*** level is reached
        !***

        final_lvl = 1

        sift_loop2: do

          !***
          !*** find the largest of the two daughters
          !***

          chk_lvl1 = 2*final_lvl
          chk_lvl2 = 2*final_lvl+1
          if (chk_lvl2 >= lvl) chk_lvl2 = chk_lvl1

          if ((add1(chk_lvl1) >  add1(chk_lvl2)) .OR.
     &       ((add1(chk_lvl1) == add1(chk_lvl2)) .AND.
     &        (add2(chk_lvl1) >  add2(chk_lvl2)))) then
            max_lvl = chk_lvl1
          else 
            max_lvl = chk_lvl2
          endif

          !***
          !*** if the parent is greater than both daughters,
          !*** the correct level has been found
          !***

          if ((add1_tmp >  add1(max_lvl)) .OR.
     &       ((add1_tmp == add1(max_lvl)) .AND.
     &        (add2_tmp >  add2(max_lvl)))) then
            add1(final_lvl) = add1_tmp
            add2(final_lvl) = add2_tmp
            weights(:,final_lvl) = wgttmp(:)
            exit sift_loop2

          !***
          !*** otherwise, promote the largest daughter and push
          !*** down one level in the tree.  if haven't reached
          !*** the end of the tree, repeat the process.  otherwise
          !*** store last values and exit the loop
          !***

          else 
            add1(final_lvl) = add1(max_lvl)
            add2(final_lvl) = add2(max_lvl)
            weights(:,final_lvl) = weights(:,max_lvl)

            final_lvl = max_lvl
            if (2*final_lvl >= lvl) then
              add1(final_lvl) = add1_tmp
              add2(final_lvl) = add2_tmp
              weights(:,final_lvl) = wgttmp(:)
              exit sift_loop2
            endif
          endif
        end do sift_loop2
      end do

      !***
      !*** swap the last two entries
      !***


      add1_tmp = add1(2)
      add1(2)  = add1(1)
      add1(1)  = add1_tmp

      add2_tmp = add2(2)
      add2(2)  = add2(1)
      add2(1)  = add2_tmp

      wgttmp (:)   = weights(:,2)
      weights(:,2) = weights(:,1)
      weights(:,1) = wgttmp (:)

!-----------------------------------------------------------------------

      end subroutine sort_add

!***********************************************************************

      end module remap_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
