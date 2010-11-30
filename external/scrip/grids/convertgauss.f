!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This program creates a remapping grid file for Gaussian lat/lon
!     grids (for spectral transform codes).
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: convertgauss.f,v 1.3 2000/04/19 22:05:57 pwjones Exp $
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

      program convert_gauss

!-----------------------------------------------------------------------
!
!     This file creates a remapping grid file for a Gaussian grid
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
!     T42: nx=128 ny=64
!     T62: nx=192 ny=94
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), parameter ::
     &             nx = 192, ny = 94,
     &             grid_size = nx*ny,
     &             grid_rank = 2,
     &             grid_corners = 4

      character(char_len), parameter :: 
     &             grid_name = 'T62 Gaussian Grid',
     &             grid_file_out = 'remap_grid_T62.nc'

      integer (kind=int_kind), dimension(grid_rank) ::
     &             grid_dims

!-----------------------------------------------------------------------
!
!     grid coordinates and masks
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), dimension(grid_size) ::
     &             grid_imask

      real (kind=dbl_kind), dimension(grid_size) ::
     &             grid_center_lat,  ! lat/lon coordinates for
     &             grid_center_lon   ! each grid center in degrees

      real (kind=dbl_kind), dimension(grid_corners,grid_size) ::
     &             grid_corner_lat,  ! lat/lon coordinates for
     &             grid_corner_lon   ! each grid corner in degrees

!-----------------------------------------------------------------------
!
!     other local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: i, j, iunit, atm_add

      integer (kind=int_kind) ::
     &        ncstat,            ! general netCDF status variable
     &        nc_grid_id,        ! netCDF grid dataset id
     &        nc_gridsize_id,    ! netCDF grid size dim id
     &        nc_gridcorn_id,    ! netCDF grid corner dim id
     &        nc_gridrank_id,    ! netCDF grid rank dim id
     &        nc_griddims_id,    ! netCDF grid dimension size id
     &        nc_grdcntrlat_id,  ! netCDF grid center lat id
     &        nc_grdcntrlon_id,  ! netCDF grid center lon id
     &        nc_grdimask_id,    ! netCDF grid mask id
     &        nc_grdcrnrlat_id,  ! netCDF grid corner lat id
     &        nc_grdcrnrlon_id   ! netCDF grid corner lon id

      integer (kind=int_kind), dimension(2) ::
     &        nc_dims2_id        ! netCDF dim id array for 2-d arrays

      real (kind=dbl_kind) :: dlon, minlon, maxlon, centerlon,
     &                              minlat, maxlat, centerlat

      real (kind=dbl_kind), dimension(ny) :: gauss_root, gauss_wgt

!-----------------------------------------------------------------------
!
!     compute longitudes of cell centers and corners.  set up alon
!     array for search routine.
!
!-----------------------------------------------------------------------

      grid_dims(1) = nx
      grid_dims(2) = ny

      dlon = 360./nx

      do i=1,nx

        centerlon = (i-1)*dlon
        minlon = centerlon - half*dlon
        maxlon = centerlon + half*dlon

        do j=1,ny
          atm_add = (j-1)*nx + i

          grid_center_lon(atm_add  ) = centerlon
          grid_corner_lon(1,atm_add) = minlon
          grid_corner_lon(2,atm_add) = maxlon
          grid_corner_lon(3,atm_add) = maxlon
          grid_corner_lon(4,atm_add) = minlon
        end do

      end do

!-----------------------------------------------------------------------
!
!     compute Gaussian latitudes and store in gauss_wgt.
!
!-----------------------------------------------------------------------

      call gquad(ny, gauss_root, gauss_wgt)
      do j=1,ny
        gauss_wgt(j) = pih - gauss_root(ny+1-j)
      end do

!-----------------------------------------------------------------------
!
!     compute latitudes at cell centers and corners.  set up alat 
!     array for search routine.
!
!-----------------------------------------------------------------------

      do j=1,ny
        centerlat = gauss_wgt(j)

        if (j .eq. 1) then
          minlat = -pih
        else
          minlat = ATAN((COS(gauss_wgt(j-1)) - 
     &                   COS(gauss_wgt(j  )))/
     &                  (SIN(gauss_wgt(j  )) - 
     &                   SIN(gauss_wgt(j-1))))
        endif

        if (j .eq. ny) then
          maxlat = pih
        else
          maxlat = ATAN((COS(gauss_wgt(j  )) - 
     &                   COS(gauss_wgt(j+1)))/
     &                  (SIN(gauss_wgt(j+1)) - 
     &                   SIN(gauss_wgt(j  ))))
        endif

        do i=1,nx
          atm_add = (j-1)*nx + i
          grid_center_lat(atm_add  ) = centerlat*360./pi2
          grid_corner_lat(1,atm_add) = minlat*360./pi2
          grid_corner_lat(2,atm_add) = minlat*360./pi2
          grid_corner_lat(3,atm_add) = maxlat*360./pi2
          grid_corner_lat(4,atm_add) = maxlat*360./pi2
        end do

      end do

!-----------------------------------------------------------------------
!
!     define mask
!
!-----------------------------------------------------------------------

      grid_imask = 1

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
      !*** define grid corner dimension
      !***

      ncstat = nf_def_dim (nc_grid_id, 'grid_corners', grid_corners, 
     &                     nc_gridcorn_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid rank dimension
      !***

      ncstat = nf_def_dim (nc_grid_id, 'grid_rank', grid_rank, 
     &                     nc_gridrank_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid dimension size array
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
     &                          7, 'degrees')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid center longitude array
      !***

      ncstat = nf_def_var (nc_grid_id, 'grid_center_lon', NF_DOUBLE,
     &                     1, nc_gridsize_id, nc_grdcntrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcntrlon_id, 'units',
     &                          7, 'degrees')
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
     &                          7, 'degrees')
      call netcdf_error_handler(ncstat)

      !***
      !*** define grid corner longitude array
      !***

      ncstat = nf_def_var (nc_grid_id, 'grid_corner_lon', NF_DOUBLE,
     &                     2, nc_dims2_id, nc_grdcrnrlon_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_put_att_text (nc_grid_id, nc_grdcrnrlon_id, 'units',
     &                          7, 'degrees')
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
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------

      end program convert_gauss

!***********************************************************************

      subroutine gquad(l,root,w)

!-----------------------------------------------------------------------
!
!     This subroutine finds the l roots (in theta) and gaussian weights 
!     associated with the legendre polynomial of degree l > 1.
!
!-----------------------------------------------------------------------

      use kinds_mod
      use constants

      implicit none

!-----------------------------------------------------------------------
!
!     intent(in)
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: l

!-----------------------------------------------------------------------
!
!     intent(out)
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), dimension(l), intent(out) ::
     &      root, w

!-----------------------------------------------------------------------
!
!     local
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: l1, l2, l22, l3, k, i, j

      real (kind=dbl_kind) :: 
     &     del,co,p1,p2,p3,t1,t2,slope,s,c,pp1,pp2,p00

!-----------------------------------------------------------------------
!
!     Define useful constants.
!
!-----------------------------------------------------------------------

      del= pi/float(4*l)
      l1 = l+1
      co = float(2*l+3)/float(l1**2)
      p2 = 1.0
      t2 = -del
      l2 = l/2
      k = 1
      p00 = one/sqrt(two)

!-----------------------------------------------------------------------
!
!     Start search for each root by looking for crossing point.
!
!-----------------------------------------------------------------------

      do i=1,l2
   10    t1 = t2
         t2 = t1+del
         p1 = p2
         s = sin(t2)
         c = cos(t2)
         pp1 = 1.0
         p3 = p00
         do j=1,l1
            pp2 = pp1
            pp1 = p3
            p3 = 2.0*sqrt((float(j**2)-0.250)/float(j**2))*c*pp1-
     &           sqrt(float((2*j+1)*(j-1)*(j-1))/
     &           float((2*j-3)*j*j))*pp2
         end do
         p2 = pp1
         if ((k*p2).gt.0) goto 10

!-----------------------------------------------------------------------
!
!        Now converge using Newton-Raphson.
!
!-----------------------------------------------------------------------

         k = -k
   20    continue
            slope = (t2-t1)/(p2-p1)
            t1 = t2
            t2 = t2-slope*p2
            p1 = p2
            s = sin(t2)
            c = cos(t2)
            pp1 = 1.0
            p3 = p00
            do j=1,l1
               pp2 = pp1
               pp1 = p3
               p3 = 2.0*sqrt((float(j**2)-0.250)/float(j**2))*c*pp1-
     &              sqrt(float((2*j+1)*(j-1)*(j-1))/
     &              float((2*j-3)*j*j))*pp2
            end do
            p2 = pp1
         if (abs(p2).gt.1.e-10) goto 20
         root(i) = t2
         w(i) = co*(sin(t2)/p3)**2
      end do

!-----------------------------------------------------------------------
!
!     If l is odd, take care of odd point.
!
!-----------------------------------------------------------------------

      l22 = 2*l2
      if (l22 .ne. l) then
         l2 = l2+1
         t2 = pi/2.0
         root(l2) = t2
         s = sin(t2)
         c = cos(t2)
         pp1 = 1.0
         p3 = p00
         do j=1,l1
            pp2 = pp1
            pp1 = p3
            p3 = 2.0*sqrt((float(j**2)-0.250)/float(j**2))*c*pp1-
     &           sqrt(float((2*j+1)*(j-1)*(j-1))/
     &           float((2*j-3)*j*j))*pp2
         end do
         p2 = pp1
         w(l2) = co/p3**2
      endif

!-----------------------------------------------------------------------
!
!     Use symmetry to compute remaining roots and weights.
!
!-----------------------------------------------------------------------

      l3 = l2+1
      do i=l3,l
         root(i) = pi-root(l-i+1)
         w(i) = w(l-i+1)
      end do

!-----------------------------------------------------------------------

      end subroutine gquad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
