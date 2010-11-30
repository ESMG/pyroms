!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This routine is the driver for computing the addresses and weights 
!     for interpolating between two grids on a sphere.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: scrip.f,v 1.6 2001/08/21 21:06:44 pwjones Exp $
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

      program scrip

!-----------------------------------------------------------------------

      use kinds_mod                  ! module defining data types
      use constants                  ! module for common constants
      use iounits                    ! I/O unit manager
      use timers                     ! CPU timers
      use grids                      ! module with grid information
      use remap_vars                 ! common remapping variables
      use remap_conservative         ! routines for conservative remap
      use remap_distance_weight      ! routines for dist-weight remap
      use remap_bilinear             ! routines for bilinear interp
      use remap_bicubic              ! routines for bicubic  interp
      use remap_write                ! routines for remap output

      implicit none

!-----------------------------------------------------------------------
!
!     input namelist variables
!
!-----------------------------------------------------------------------

      character (char_len) :: 
     &           grid1_file,   ! filename of grid file containing grid1
     &           grid2_file,   ! filename of grid file containing grid2
     &           interp_file1, ! filename for output remap data (map1)
     &           interp_file2, ! filename for output remap data (map2)
     &           map1_name,    ! name for mapping from grid1 to grid2
     &           map2_name,    ! name for mapping from grid2 to grid1
     &           map_method,   ! choice for mapping method
     &           normalize_opt,! option for normalizing weights
     &           output_opt    ! option for output conventions

      integer (kind=int_kind) ::
     &           nmap          ! number of mappings to compute (1 or 2)

      namelist /remap_inputs/ grid1_file, grid2_file, 
     &                        interp_file1, interp_file2,
     &                        map1_name, map2_name, num_maps,
     &                        luse_grid1_area, luse_grid2_area,
     &                        map_method, normalize_opt, output_opt,
     &                        restrict_type, num_srch_bins,
     &                        grid1_periodic, grid2_periodic

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n,     ! dummy counter
     &                           iunit  ! unit number for namelist file

!-----------------------------------------------------------------------
!
!     initialize timers
!
!-----------------------------------------------------------------------

      call timers_init
      do n=1,max_timers
        call timer_clear(n)
      end do

!-----------------------------------------------------------------------
!
!     read input namelist
!
!-----------------------------------------------------------------------

      grid1_file    = 'unknown'
      grid2_file    = 'unknown'
      interp_file1  = 'unknown'
      interp_file2  = 'unknown'
      map1_name     = 'unknown'
      map2_name     = 'unknown'
      luse_grid1_area   = .false.
      luse_grid2_area   = .false.
      grid1_periodic(1) = .true.
      grid1_periodic(2) = .true.
      grid2_periodic(1) = .true.
      grid2_periodic(2) = .true.
      num_maps      = 2
      map_type      = 1
      normalize_opt = 'fracarea'
      output_opt    = 'scrip'
      restrict_type = 'latitude'
      num_srch_bins = 900

      call get_unit(iunit)
      open(iunit, file='scrip_in', status='old', form='formatted')
      read(iunit, nml=remap_inputs)
      call release_unit(iunit)

      select case(map_method)
      case ('conservative')
        map_type = map_type_conserv
        luse_grid_centers = .false.
      case ('bilinear')
        map_type = map_type_bilinear
        luse_grid_centers = .true.
      case ('bicubic')
        map_type = map_type_bicubic
        luse_grid_centers = .true.
      case ('distwgt')
        map_type = map_type_distwgt
        luse_grid_centers = .true.
      case default
        stop 'unknown mapping method'
      end select

      select case(normalize_opt(1:4))
      case ('none')
        norm_opt = norm_opt_none
      case ('frac')
        norm_opt = norm_opt_frcarea
      case ('dest')
        norm_opt = norm_opt_dstarea
      case default
        stop 'unknown normalization option'
      end select

!-----------------------------------------------------------------------
!
!     initialize grid information for both grids
!
!-----------------------------------------------------------------------

      call grid_init(grid1_file, grid2_file)

      write(stdout, *) ' Computing remappings between: ',grid1_name
      write(stdout, *) '                          and  ',grid2_name

!-----------------------------------------------------------------------
!
!     initialize some remapping variables.
!
!-----------------------------------------------------------------------

      call init_remap_vars

!-----------------------------------------------------------------------
!
!     call appropriate interpolation setup routine based on type of
!     remapping requested.
!
!-----------------------------------------------------------------------

      select case(map_type)
      case(map_type_conserv)
        call remap_conserv
      case(map_type_bilinear)
        call remap_bilin
      case(map_type_distwgt)
        call remap_distwgt
      case(map_type_bicubic)
        call remap_bicub
      case default
        stop 'Invalid Map Type'
      end select

!-----------------------------------------------------------------------
!
!     reduce size of remapping arrays and then write remapping info
!     to a file.
!
!-----------------------------------------------------------------------

      if (num_links_map1 /= max_links_map1) then
        call resize_remap_vars(1, num_links_map1-max_links_map1)
      endif
      if ((num_maps > 1) .and. (num_links_map2 /= max_links_map2)) then
        call resize_remap_vars(2, num_links_map2-max_links_map2)
      endif

      call write_remap(map1_name, map2_name, 
     &                 interp_file1, interp_file2, output_opt)

!-----------------------------------------------------------------------

      end program scrip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
