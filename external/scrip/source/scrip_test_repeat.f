!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this program is a short driver that tests the remappings using
!     a simple analytic field.  the results are written in netCDF
!     format.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: scrip_test_repeat.f,v 1.3 2000/04/19 21:56:26 pwjones Exp $
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

      program remap_test_repeat

!-----------------------------------------------------------------------

      use kinds_mod    ! defines common data types
      use constants    ! defines common constants
      use iounits      ! I/O unit manager
      use netcdf_mod   ! netcdf I/O stuff
      use grids        ! module containing grid info
      use remap_vars   ! module containing remapping info

      implicit none

!-----------------------------------------------------------------------
!
!     interface for remap routine
!
!-----------------------------------------------------------------------

      interface
        subroutine remap(dst_array, map_wts, dst_add, src_add, 
     &                   src_array, src_grad1, src_grad2, src_grad3)

        use kinds_mod
        use constants

        implicit none

        integer (kind=int_kind), dimension(:), intent(in) ::
     &       dst_add, src_add 

        real (kind=dbl_kind), dimension(:,:), intent(in) :: map_wts

        real (kind=dbl_kind), dimension(:), intent(in) :: src_array

        real (kind=dbl_kind), dimension(:), intent(in), optional ::
     &     src_grad1, src_grad2, src_grad3

        real (kind=dbl_kind), dimension(:), intent(inout) :: dst_array

        end subroutine remap
      end interface

!-----------------------------------------------------------------------
!
!     input namelist variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) ::
     &           field_choice, ! choice of field to be interpolated
     &           num_repeats   ! number of times to repeat remappings

      character (char_len) :: 
     &           interp_file1, ! filename containing remap data (map1)
     &           interp_file2, ! filename containing remap data (map2)
     &           output_file   ! filename containing output test data

      namelist /remap_inputs/ field_choice, num_repeats,
     &                        interp_file1, interp_file2, output_file

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (char_len) :: 
     &           map_name1,    ! name for mapping from grid1 to grid2
     &           map_name2     ! name for mapping from grid2 to grid1

      integer (kind=int_kind) ::    ! netCDF ids for files and arrays
     &        n, ncstat, nc_outfile_id, 
     &        nc_srcgrdcntrlat_id, nc_srcgrdcntrlon_id,
     &        nc_dstgrdcntrlat_id, nc_dstgrdcntrlon_id,
     &        nc_srcgrdrank_id, nc_dstgrdrank_id,
     &        nc_srcgrdimask_id, nc_dstgrdimask_id,
     &        nc_srcgrdarea_id, nc_dstgrdarea_id,
     &        nc_srcgrdfrac_id, nc_dstgrdfrac_id,
     &        nc_srcarray_id, nc_srcgradlat_id, nc_srcgradlon_id,
     &        nc_dstarray1_id, nc_dstarray2_id,
     &        nc_dsterror1_id, nc_dsterror2_id

      integer (kind=int_kind), dimension(:), allocatable ::
     &        nc_grid1size_id, nc_grid2size_id

!-----------------------------------------------------------------------

      character (char_len) :: 
     &          dim_name    ! netCDF dimension name

      integer (kind=int_kind) :: i,j,n,
     &    iunit  ! unit number for namelist file

      integer (kind=int_kind), dimension(:), allocatable ::
     &    grid1_imask, grid2_imask

      real (kind=dbl_kind) ::
     &    length            ! length scale for cosine hill test field

      real (kind=dbl_kind), dimension(:), allocatable ::
     &    grid1_array, 
     &    grid1_tmp, 
     &    grad1_lat, 
     &    grad1_lon, 
     &    grid2_array, 
     &    grid2_err,
     &    grid2_tmp,
     &    grad2_lat, 
     &    grad2_lon  

!-----------------------------------------------------------------------
!
!     read namelist for file and mapping info
!
!-----------------------------------------------------------------------

      call get_unit(iunit)
      open(iunit, file='repeat_test_in', status='old', form='formatted')
      read(iunit, nml=remap_inputs)
      call release_unit(iunit)
      write(*,nml=remap_inputs)

!-----------------------------------------------------------------------
!
!     read remapping data
!
!-----------------------------------------------------------------------

      num_maps = 2

      call read_remap(map_name1, map_name2, interp_file1, interp_file2)

!-----------------------------------------------------------------------
!
!     allocate arrays
!
!-----------------------------------------------------------------------

      allocate (grid1_array    (grid1_size), 
     &          grid1_tmp      (grid1_size),
     &          grad1_lat      (grid1_size), 
     &          grad1_lon      (grid1_size), 
     &          grid1_imask    (grid1_size),
     &          grid2_array    (grid2_size), 
     &          grid2_err      (grid2_size),
     &          grid2_tmp      (grid2_size),
     &          grad2_lat      (grid2_size), 
     &          grad2_lon      (grid2_size), 
     &          grid2_imask    (grid2_size))

      where (grid1_mask)
        grid1_imask = 1
      elsewhere
        grid1_imask = 0
      endwhere
      where (grid2_mask)
        grid2_imask = 1
      elsewhere
        grid2_imask = 0
      endwhere

!-----------------------------------------------------------------------
!
!     setup a NetCDF file for output
!
!-----------------------------------------------------------------------

      !***
      !*** create netCDF dataset 
      !***

      ncstat = nf_create (output_file, NF_CLOBBER, nc_outfile_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, NF_GLOBAL, 'title',
     &                          len_trim(map_name1), map_name1)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid size dimensions
      !***

      allocate( nc_grid1size_id(grid1_rank),
     &          nc_grid2size_id(grid2_rank))

      do n=1,grid1_rank
        write(dim_name,1000) 'grid1_dim',n
        ncstat = nf_def_dim (nc_outfile_id, dim_name, 
     &                       grid1_dims(n), nc_grid1size_id(n))
        call netcdf_error_handler(ncstat)
      end do

      do n=1,grid2_rank
        write(dim_name,1000) 'grid2_dim',n
        ncstat = nf_def_dim (nc_outfile_id, dim_name, 
     &                       grid2_dims(n), nc_grid2size_id(n))
        call netcdf_error_handler(ncstat)
      end do
 1000 format(a9,i1)

      !***
      !*** define grid center latitude array
      !***

      ncstat = nf_def_var (nc_outfile_id, 'src_grid_center_lat', 
     &                     NF_DOUBLE, grid1_rank, nc_grid1size_id, 
     &                     nc_srcgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_srcgrdcntrlat_id, 
     &                          'units', 7, 'radians')
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'dst_grid_center_lat', 
     &                     NF_DOUBLE, grid2_rank, nc_grid2size_id, 
     &                     nc_dstgrdcntrlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_dstgrdcntrlat_id, 
     &                          'units', 7, 'radians')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid center longitude array
      !***

      ncstat = nf_def_var (nc_outfile_id, 'src_grid_center_lon', 
     &                     NF_DOUBLE, grid1_rank, nc_grid1size_id, 
     &                     nc_srcgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_srcgrdcntrlon_id, 
     &                          'units', 7, 'radians')
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'dst_grid_center_lon', 
     &                     NF_DOUBLE, grid2_rank, nc_grid2size_id, 
     &                     nc_dstgrdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_dstgrdcntrlon_id, 
     &                          'units', 7, 'radians')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid mask
      !***

      ncstat = nf_def_var (nc_outfile_id, 'src_grid_imask', NF_INT,
     &               grid1_rank, nc_grid1size_id, nc_srcgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_srcgrdimask_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'dst_grid_imask', NF_INT,
     &               grid2_rank, nc_grid2size_id, nc_dstgrdimask_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_outfile_id, nc_dstgrdimask_id, 
     &                          'units', 8, 'unitless')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid area arrays
      !***

      ncstat = nf_def_var (nc_outfile_id, 'src_grid_area', 
     &                     NF_DOUBLE, grid1_rank, nc_grid1size_id, 
     &                     nc_srcgrdarea_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'dst_grid_area', 
     &                     NF_DOUBLE, grid2_rank, nc_grid2size_id, 
     &                     nc_dstgrdarea_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid fraction arrays
      !***

      ncstat = nf_def_var (nc_outfile_id, 'src_grid_frac', 
     &                     NF_DOUBLE, grid1_rank, nc_grid1size_id, 
     &                     nc_srcgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'dst_grid_frac', 
     &                     NF_DOUBLE, grid2_rank, nc_grid2size_id, 
     &                     nc_dstgrdfrac_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define source array
      !***

      ncstat = nf_def_var (nc_outfile_id, 'src_array', 
     &                     NF_DOUBLE, grid1_rank, nc_grid1size_id, 
     &                     nc_srcarray_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define gradient arrays
      !***

      ncstat = nf_def_var (nc_outfile_id, 'src_grad_lat', 
     &                     NF_DOUBLE, grid1_rank, nc_grid1size_id, 
     &                     nc_srcgradlat_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'src_grad_lon', 
     &                     NF_DOUBLE, grid1_rank, nc_grid1size_id, 
     &                     nc_srcgradlon_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define destination arrays
      !***

      ncstat = nf_def_var (nc_outfile_id, 'dst_array1', 
     &                     NF_DOUBLE, grid2_rank, nc_grid2size_id, 
     &                     nc_dstarray1_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'dst_array2', 
     &                     NF_DOUBLE, grid2_rank, nc_grid2size_id, 
     &                     nc_dstarray2_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define error arrays
      !***

      ncstat = nf_def_var (nc_outfile_id, 'dst_error1', 
     &                     NF_DOUBLE, grid2_rank, nc_grid2size_id, 
     &                     nc_dsterror1_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_def_var (nc_outfile_id, 'dst_error2', 
     &                     NF_DOUBLE, grid2_rank, nc_grid2size_id, 
     &                     nc_dsterror2_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** end definition stage
      !***

      ncstat = nf_enddef(nc_outfile_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     write some grid info
!
!-----------------------------------------------------------------------

      !***
      !*** write grid center latitude array
      !***

      ncstat = nf_put_var_double(nc_outfile_id, nc_srcgrdcntrlat_id,
     &                           grid1_center_lat)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_dstgrdcntrlat_id,
     &                           grid2_center_lat)
      call netcdf_error_handler(ncstat)

      !***
      !*** write grid center longitude array
      !***

      ncstat = nf_put_var_double(nc_outfile_id, nc_srcgrdcntrlon_id,
     &                           grid1_center_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_dstgrdcntrlon_id,
     &                           grid2_center_lon)
      call netcdf_error_handler(ncstat)

      !***
      !*** write grid mask
      !***

      ncstat = nf_put_var_int(nc_outfile_id, nc_srcgrdimask_id,
     &                        grid1_imask)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_int(nc_outfile_id, nc_dstgrdimask_id,
     &                        grid2_imask)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid area arrays
      !***

      ncstat = nf_put_var_double(nc_outfile_id, nc_srcgrdarea_id,
     &                           grid1_area)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_dstgrdarea_id,
     &                           grid2_area)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid fraction arrays
      !***

      ncstat = nf_put_var_double(nc_outfile_id, nc_srcgrdfrac_id,
     &                           grid1_frac)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_dstgrdfrac_id,
     &                           grid2_frac)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     set up fields for test cases based on user choice
!
!-----------------------------------------------------------------------

      select case (field_choice)
      case(1)  !*** cosine hill at lon=pi and lat=0

        length = 0.1*pi2

        grid1_array = cos(grid1_center_lat)*cos(grid1_center_lon)
        grid2_array = cos(grid2_center_lat)*cos(grid2_center_lon)

        grid1_tmp = acos(-grid1_array)/length
        grid2_tmp = acos(-grid2_array)/length

        where (grid1_tmp <= one)
          grad1_lat   = (pi/length)*sin(pi*grid1_tmp)*
     &                  sin(grid1_center_lat)*cos(grid1_center_lon)/
     &                  sqrt(one-grid1_array**2)
          grad1_lon   = (pi/length)*sin(pi*grid1_tmp)*
     &                  sin(grid1_center_lon)/
     &                  sqrt(one-grid1_array**2)
          grid1_array = two + cos(pi*grid1_tmp)
        elsewhere
          grid1_array = one
          grad1_lat   = zero
          grad1_lon   = zero
        endwhere
        
        where (grid2_tmp <= one)
          grad2_lat   = (pi/length)*sin(pi*grid2_tmp)*
     &                  sin(grid2_center_lat)*cos(grid2_center_lon)/
     &                  sqrt(one-grid2_array**2)
          grad2_lon   = (pi/length)*sin(pi*grid2_tmp)*
     &                  sin(grid2_center_lon)/
     &                  sqrt(one-grid2_array**2)
          grid2_array = two + cos(pi*grid2_tmp)
        elsewhere
          grid2_array = one
          grad2_lat   = zero
          grad2_lon   = zero
        endwhere
        
        where (.not. grid1_mask)
          grid1_array = zero
          grad1_lat   = zero
          grad1_lon   = zero
        end where

        where (.not. grid2_mask)
          grid2_array = zero
          grad2_lat   = zero
          grad2_lon   = zero
        end where

      case(2)  !*** pseudo-spherical harmonic l=2,m=2

        where (grid1_mask)
          grid1_array = two + cos(grid1_center_lat)**2*
     &                    cos(two*grid1_center_lon)
          grad1_lat   = -sin(two*grid1_center_lat)*
     &                   cos(two*grid1_center_lon)
          grad1_lon   = -two*cos(grid1_center_lat)*
     &                   sin(two*grid1_center_lon)
        elsewhere
          grid1_array = zero
          grad1_lat   = zero
          grad1_lon   = zero
        end where

        where (grid2_mask)
          grid2_array = two + cos(grid2_center_lat)**2*
     &                        cos(two*grid2_center_lon)
          grad2_lat   = -sin(two*grid2_center_lat)*
     &                   cos(two*grid2_center_lon)
          grad2_lon   = -two*cos(grid2_center_lat)*
     &                   sin(two*grid2_center_lon)
        elsewhere
          grid2_array = zero
          grad2_lat   = zero
          grad2_lon   = zero
        end where

      case(3)  !*** pseudo-spherical harmonic l=32, m=16

        where (grid1_mask)
          grid1_array = two + sin(two*grid1_center_lat)**16*
     &                        cos(16.*grid1_center_lon)
          grad1_lat   = 32.*sin(two*grid1_center_lat)**15*
     &                      cos(two*grid1_center_lat)*
     &                      cos(16.*grid1_center_lon)
          grad1_lon   = -32.*sin(two*grid1_center_lat)**15*
     &                       sin(grid1_center_lat)*
     &                   sin(16.*grid1_center_lon)
        elsewhere
          grid1_array = zero
          grad1_lat   = zero
          grad1_lon   = zero
        end where

        where (grid2_mask)
          grid2_array = two + sin(two*grid2_center_lat)**16*
     &                        cos(16.*grid2_center_lon)
          grad2_lat   = 32.*sin(two*grid2_center_lat)**15*
     &                      cos(two*grid2_center_lat)*
     &                      cos(16.*grid2_center_lon)
          grad2_lon   = -32.*sin(two*grid2_center_lat)**15*
     &                       sin(grid2_center_lat)*
     &                   sin(16.*grid2_center_lon)
        elsewhere
          grid2_array = zero
          grad2_lat   = zero
          grad2_lon   = zero
        end where

      case default

        stop 'Bad choice for field to interpolate'

      end select

!-----------------------------------------------------------------------
!
!     test repeated first-order maps between grid1 and grid2
!
!-----------------------------------------------------------------------

      call remap(grid2_tmp, wts_map1, grid2_add_map1, grid1_add_map1,
     &           grid1_array)
      do n=1,num_repeats
        call remap(grid1_tmp, wts_map2, grid1_add_map2, grid2_add_map2,
     &             grid2_tmp)
        call remap(grid2_tmp, wts_map1, grid2_add_map1, grid1_add_map1,
     &             grid1_tmp)
      end do

      where (grid2_frac > .999)
        grid2_err = (grid2_tmp - grid2_array)/grid2_array
      elsewhere 
        grid2_err = zero
      end where

      print *,'First order mapping from grid1 to grid2:'
      print *,'----------------------------------------'
      print *,'Grid1 min,max: ',minval(grid1_array),maxval(grid1_array)
      print *,'Grid2 min,max: ',minval(grid2_tmp  ),maxval(grid2_tmp  )
      print *,' Err2 min,max: ',minval(grid2_err),maxval(grid2_err)
      print *,' Err2    mean: ',sum(abs(grid2_err))/
     &                          count(grid2_frac > .001)

      !***
      !*** Conservation Test
      !***

      print *,'Conservation:'
      print *,'Grid1 Integral = ',sum(grid1_array*grid1_area*grid1_frac)
      print *,'Grid2 Integral = ',sum(grid2_tmp  *grid2_area*grid2_frac)

!-----------------------------------------------------------------------
!
!     write results to NetCDF file
!
!-----------------------------------------------------------------------

      ncstat = nf_put_var_double(nc_outfile_id, nc_srcarray_id,
     &                           grid1_array)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_dstarray1_id,
     &                           grid2_tmp  )
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_dsterror1_id,
     &                           grid2_err)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     test repeated second-order mapppings between grid1 and grid2
!
!-----------------------------------------------------------------------

      if (num_wts > 1) then

      call remap(grid2_tmp  , wts_map1, grid2_add_map1, grid1_add_map1,
     &           grid1_array, src_grad1=grad1_lat,
     &                        src_grad2=grad1_lon) 

      do n=1,num_repeats
        call remap(grid1_tmp, wts_map2, grid1_add_map2, grid2_add_map2,
     &             grid2_tmp, src_grad1=grad2_lat,
     &                        src_grad2=grad2_lon) 

        call remap(grid2_tmp, wts_map1, grid2_add_map1, grid1_add_map1,
     &             grid1_tmp, src_grad1=grad1_lat,
     &                        src_grad2=grad1_lon) 

      end do

      where (grid2_frac > .999)
        grid2_err = (grid2_tmp - grid2_array)/grid2_array
      elsewhere 
        grid2_err = zero
      end where

      print *,'Second order mapping from grid1 to grid2:'
      print *,'-----------------------------------------'
      print *,'Grid1 min,max: ',minval(grid1_array),maxval(grid1_array)
      print *,'Grid2 min,max: ',minval(grid2_tmp  ),maxval(grid2_tmp  )
      print *,' Err2 min,max: ',minval(grid2_err),maxval(grid2_err)
      print *,' Err2    mean: ',sum(abs(grid2_err))/
     &                          count(grid2_frac > .001)

      !***
      !*** Conservation Test
      !***

      print *,'Conservation:'
      print *,'Grid1 Integral = ',sum(grid1_array*grid1_area*grid1_frac)
      print *,'Grid2 Integral = ',sum(grid2_tmp  *grid2_area*grid2_frac)

!-----------------------------------------------------------------------
!
!     write results to NetCDF file
!
!-----------------------------------------------------------------------

      ncstat = nf_put_var_double(nc_outfile_id, nc_srcgradlon_id,
     &                           grad1_lon)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_dstarray2_id,
     &                           grid2_tmp  )
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_var_double(nc_outfile_id, nc_dsterror2_id,
     &                           grid2_err)
      call netcdf_error_handler(ncstat)

      endif

!-----------------------------------------------------------------------
!
!     close netCDF file
!
!-----------------------------------------------------------------------

      ncstat = nf_close(nc_outfile_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------

      end program remap_test_repeat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
