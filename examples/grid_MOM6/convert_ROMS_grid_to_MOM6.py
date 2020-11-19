#!/bin/env python
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 :

import sys
import numpy
import netCDF4
import os.path
import datetime
import subprocess
import warnings

import Spherical

# This script converts a ROMS horizontal grid into a MOM6 horizontal
# grid.

# Both ROMS and MOM6 horizontal grids use an Arakawa C-grid, with four
# types of points:
#   rho: the centers of the cells
#   psi: the corners of the cells, located diagonally between the
#        'rho' points
#   u:   the u-velocity points, located between 'rho' points in the
#        east/west direction
#   v:   the v-velocity points, located between 'rho' points in the
#        north/south direction

# The main differences between the two grids are:
#  * the outermost points of the ROMS grid are the 'rho' points, while
#    the outermost points of the MOM6 grid are the 'psi' points (both
#    with interspersed 'u' and 'v' points); and
#  * the MOM6 grid interleaves all four types of points into a single
#    "supergrid", while ROMS stores them as separate grids.

# The ROMS grid looks like this, with an extra layer of 'rho' points
# around the outside:
# (see https://www.myroms.org/wiki/Numerical_Solution_Technique)
#
#       p - p - p - p - p
#    3  | + | + | + | + |     p = rho (center) points
#       p - p - p - p - p     + = psi (corner) points
#    2  | + | + | + | + |     - = u points
#       p - p - p - p - p     | = v points
#    1  | + | + | + | + |
#       p - p - p - p - p
#
#         1   2   3   4

# The MOM6 grid has 'psi' points on the outside, not 'rho':
#
#    3    + | + | + | +       p = rho (center) points
#         - p - p - p -       + = psi (corner) points
#    2    + | + | + | +       - = u points
#         - p - p - p -       | = v points
#    1    + | + | + | +
#
#         1   2   3   4

def read_ROMS_grid(roms_grid_filename):
    """Load the ROMS grid from a NetCDF file."""

    subgrids = ['psi', 'rho', 'u', 'v'] # also available: 'vert'
    fields = ['lon', 'lat', 'x', 'y', 'mask']

    roms_grid = dict()
    with netCDF4.Dataset(roms_grid_filename) as roms_ds:
        # extract fields named based on which grid they are on
        for subgrid in subgrids:
            roms_grid[subgrid] = dict()
            for field in fields:
                var_name = field + '_' + subgrid
                roms_grid[subgrid][field] = roms_ds.variables[var_name][:]
                if (field == 'x') or (field == 'y'):
                    units = roms_ds.variables[var_name].units.lower()
                    assert units.startswith('meter')
                elif (field == 'lat') or (field == 'lon'):
                    units = roms_ds.variables[var_name].units.lower()
                    assert units.startswith('degree')

        # extract fields that don't follow the above naming pattern
        roms_grid['rho']['h'] = roms_ds.variables['h'][:] # on the rho grid, but not called "h_rho"

        roms_grid['metadata'] = dict()

        spherical = roms_ds.variables['spherical'][:]
        if (spherical == 0) or (spherical == b'F') or (spherical == b'f'):
            roms_grid['metadata']['is_spherical'] = False
        elif (spherical == 1) or (spherical == b'T') or (spherical == b't'):
            roms_grid['metadata']['is_spherical'] = True
        else:
            warnings.warn('Unrecognized value for spherical in ROMS grid: %s' % str(spherical))

    return roms_grid

def trim_ROMS_grid(old_grid):
    """Remove extraneous points on the outside of the ROMS grid."""

    trim_subgrid = dict()
    # remove the outer:   ( rows,  cols)
    trim_subgrid['psi'] = (False, False) # Cell corners (leave alone)
    trim_subgrid['rho'] = ( True,  True) # Cell centers (remove outer row and column)
    trim_subgrid[ 'u' ] = ( True, False) # U-points (remove outer row)
    trim_subgrid[ 'v' ] = (False,  True) # V-points (remove outer column)

    new_grid = dict()
    for subgrid in old_grid.keys():
        if subgrid == 'metadata':
            new_grid[subgrid] = dict(old_grid[subgrid])
            continue
        new_grid[subgrid] = dict()
        trim_rows,trim_cols = trim_subgrid[subgrid]
        for field in old_grid[subgrid].keys():
            if trim_rows and trim_cols:
                new_grid[subgrid][field] = old_grid[subgrid][field][1:-1,1:-1]
            elif trim_rows:
                new_grid[subgrid][field] = old_grid[subgrid][field][1:-1, :  ]
            elif trim_cols:
                new_grid[subgrid][field] = old_grid[subgrid][field][ :  ,1:-1]
            else:
                new_grid[subgrid][field] = old_grid[subgrid][field][ :  , :  ]

    return new_grid

def get_git_repo_version_info():
    """Describe the current version of this script as known by Git."""
    repo_name = 'ESMG/PyCNAL'
    git_command = ['git', 'describe', '--all', '--long', '--dirty', '--abbrev=10']
    description =  subprocess.check_output(git_command, universal_newlines=True).rstrip()
    return repo_name + ': ' + description

def get_history_entry(argv):
    """Construct an entry for the global 'history' attribute of a NetCDF file,
    which is a date and the command used."""
    today = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d')
    command = ' '.join(argv)
    return today + ': ' + command

def setup_MOM6_grid(argv):
    tile_str = 'tile1'
    atmos_tile = 'atmos_mosaic_' + tile_str
    land_tile  =  'land_mosaic_' + tile_str
    ocean_tile = 'ocean_mosaic_' + tile_str

    mom6_grid = dict()

    mom6_grid['filenames'] = dict()

    # Always use './' as the directory for files, since FMS always runs from the
    # "INPUT" directory
    mom6_grid['filenames']['directory']            = './'

    mom6_grid['filenames']['supergrid']            = 'ocean_hgrid.nc'
    mom6_grid['filenames']['topography']           = 'ocean_topog.nc'
    mom6_grid['filenames']['mosaic']               = 'ocean_mosaic.nc'
    mom6_grid['filenames']['land_mask']            = 'land_mask.nc'
    mom6_grid['filenames']['ocean_mask']           = 'ocean_mask.nc'
    mom6_grid['filenames']['atmos_land_exchange']  = '%sX%s.nc' % (atmos_tile,  land_tile)
    mom6_grid['filenames']['atmos_ocean_exchange'] = '%sX%s.nc' % (atmos_tile, ocean_tile)
    mom6_grid['filenames']['land_ocean_exchange']  = '%sX%s.nc' % ( land_tile, ocean_tile)

    mom6_grid['netcdf_info'] = dict()
    mom6_grid['netcdf_info']['tile_str']      = tile_str
    mom6_grid['netcdf_info']['string_length'] = 255
    mom6_grid['netcdf_info']['grid_version']  = '0.2' # taken from make_solo_mosaic
    #mom6_grid['netcdf_info']['code_version']  = '$Name: tikal $' ### for testing
    mom6_grid['netcdf_info']['code_version']  = get_git_repo_version_info() ### for production
    mom6_grid['netcdf_info']['history_entry'] = get_history_entry(argv)

    mom6_grid['supergrid'] = dict()

    mom6_grid['cell_grid'] = dict()

    return mom6_grid

def convert_ROMS_to_MOM6(mom6_grid, roms_grid):
    """Convert the ROMS grid data into a skeleton MOM6 grid, mainly by
    merging the four sets of point locations from the ROMS grid
    into a single supergrid for MOM6."""

    ny, nx = roms_grid['rho']['lon'].shape # trimmed
    mom6_grid['cell_grid']['nx'] = nx
    mom6_grid['cell_grid']['ny'] = ny

    # Double the size of the *trimmed* ROMS grid to merge the four
    # sets of points.
    nx *= 2
    ny *= 2

    mom6_grid['supergrid']['nx'] = nx
    mom6_grid['supergrid']['ny'] = ny

    if roms_grid['metadata']['is_spherical']:
        copy_fields = ['lon', 'lat']
    else:
        copy_fields = ['x', 'y']

    # Copy points from ROMS grid
    for field in copy_fields:
        mom6_grid['supergrid'][field] = numpy.zeros((ny+1,nx+1))
        mom6_grid['supergrid'][field][ ::2, ::2] = roms_grid['psi'][field] # outer
        mom6_grid['supergrid'][field][1::2,1::2] = roms_grid['rho'][field] # inner
        mom6_grid['supergrid'][field][1::2, ::2] = roms_grid[ 'u' ][field] # between e/w
        mom6_grid['supergrid'][field][ ::2,1::2] = roms_grid[ 'v' ][field] # between n/s

    # 0 = land, 1 = water, but sometimes some huge number indicates
    # "missing" values, which we'll assume to be water
    mask = roms_grid['rho']['mask']
    mask[mask != 0] = 1
    mom6_grid['cell_grid']['depth'] = roms_grid['rho']['h'] * mask
    mom6_grid['cell_grid']['ocean_mask'] = mask
    mom6_grid['cell_grid']['land_mask'] = numpy.logical_not(mask)

    return mom6_grid

def _fill_in_MOM6_supergrid_metrics_spherical(mom6_grid):
    """Fill in missing MOM6 supergrid metrics by computing best guess
    values based on latitude and longitude coordinates."""

    lat = mom6_grid['supergrid']['lat']
    lon = mom6_grid['supergrid']['lon']

    # Approximate edge lengths as great arcs
    R = 6370.e3 # Radius of sphere
    mom6_grid['supergrid']['dx'][:,:] = R * Spherical.angle_through_center( (lat[ :,1:],lon[ :,1:]), (lat[:  ,:-1],lon[:  ,:-1]) )
    mom6_grid['supergrid']['dy'][:,:] = R * Spherical.angle_through_center( (lat[1:, :],lon[1:, :]), (lat[:-1,:  ],lon[:-1,:  ]) )

    # Approximate angles using centered differences in interior, and side differences on left/right edges
    # TODO: Why do something different at the edges when we have extra ROMS points available?
    # Because we're using a big enough footprint to need to.
    cos_lat = numpy.cos(numpy.radians(lat))
#   mom6_grid['supergrid']['angle'][:,1:-1] = numpy.arctan( (lat[:,2:] - lat[:,:-2]) / ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
#   mom6_grid['supergrid']['angle'][:, 0  ] = numpy.arctan( (lat[:, 1] - lat[:, 0 ]) / ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
#   mom6_grid['supergrid']['angle'][:,-1  ] = numpy.arctan( (lat[:,-1] - lat[:,-2 ]) / ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )
    # Compute it twice to recover from dateline problems, if any
    angle = numpy.zeros(lat.shape)
    angle2 = numpy.zeros(lat.shape)
#   angle[:,1:-1] = numpy.arctan2( (lat[:,2:] - lat[:,:-2]) , ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
#   angle[:, 0  ] = numpy.arctan2( (lat[:, 1] - lat[:, 0 ]) , ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
#   angle[:,-1  ] = numpy.arctan2( (lat[:,-1] - lat[:,-2 ]) , ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )
    lon = numpy.where(lon < 0., lon+360, lon)
    angle2[:,1:-1] = numpy.arctan2( (lat[:,2:] - lat[:,:-2]) , ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
    angle2[:, 0  ] = numpy.arctan2( (lat[:, 1] - lat[:, 0 ]) , ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
    angle2[:,-1  ] = numpy.arctan2( (lat[:,-1] - lat[:,-2 ]) , ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )
    mom6_grid['supergrid']['angle'][:,:] = numpy.maximum(angle, angle2)

    # Approximate cell areas as that of spherical polygon
    mom6_grid['supergrid']['area'][:,:] = R * R * Spherical.quad_area(lat, lon)

    return mom6_grid

def _fill_in_MOM6_supergrid_metrics_cartesian(mom6_grid):
    """Fill in missing MOM6 supergrid metrics by computing best guess
    values based on x and y coordinates."""

    x = mom6_grid['supergrid']['x']
    y = mom6_grid['supergrid']['y']

    # Compute edge lengths
    mom6_grid['supergrid']['dx'][:,:] = numpy.sqrt( (x[:,1:] - x[:,:-1])**2 + (y[:,1:] - y[:,:-1])**2 )
    mom6_grid['supergrid']['dy'][:,:] = numpy.sqrt( (x[1:,:] - x[:-1,:])**2 + (y[1:,:] - y[:-1,:])**2 )

    # Compute angles using centered differences in interior, and side differences on left/right edges
    # TODO: Why do something different at the edges when we have extra ROMS points available?
    mom6_grid['supergrid']['angle'][:,1:-1] = numpy.arctan2( (y[:,2:] - y[:,:-2]), (x[:,2:] - x[:,:-2]) )
    mom6_grid['supergrid']['angle'][:, 0  ] = numpy.arctan2( (y[:, 1] - y[:, 0 ]), (x[:, 1] - x[:, 0 ]) )
    mom6_grid['supergrid']['angle'][:,-1  ] = numpy.arctan2( (y[:,-1] - y[:,-2 ]), (x[:,-1] - x[:,-2 ]) )

    # Compute cell areas
    mom6_grid['supergrid']['area'][:,:] = mom6_grid['supergrid']['dx'][:-1, :] * mom6_grid['supergrid']['dy'][:, :-1]

    return mom6_grid

def _calculate_MOM6_cell_grid_area(mom6_grid):
    """Compute the area for the MOM6 cells (not the sub-cells of the
    supergrid)."""

    # Combine areas of smaller supergrid cells into areas of cell-grid cells
    a00 = mom6_grid['supergrid']['area'][0::2,0::2]
    a01 = mom6_grid['supergrid']['area'][0::2,1::2]
    a10 = mom6_grid['supergrid']['area'][1::2,0::2]
    a11 = mom6_grid['supergrid']['area'][1::2,1::2]
    mom6_grid['cell_grid']['area'] = a00 + a01 + a10 + a11

    return mom6_grid

def approximate_MOM6_grid_metrics(mom6_grid):
    """Fill in missing MOM6 supergrid metrics by computing best guess
    values."""

    nx = mom6_grid['supergrid']['nx']
    ny = mom6_grid['supergrid']['ny']

    # Declare shapes
    mom6_grid['supergrid']['dx']    = numpy.zeros((ny+1,nx  ))
    mom6_grid['supergrid']['dy']    = numpy.zeros((ny,  nx+1))
    mom6_grid['supergrid']['angle'] = numpy.zeros((ny+1,nx+1))
    mom6_grid['supergrid']['area']  = numpy.zeros((ny,  nx  ))

    if 'lat' in mom6_grid['supergrid']:
        mom6_grid = _fill_in_MOM6_supergrid_metrics_spherical(mom6_grid)
    else:
        mom6_grid = _fill_in_MOM6_supergrid_metrics_cartesian(mom6_grid)

    mom6_grid = _calculate_MOM6_cell_grid_area(mom6_grid)

    return mom6_grid

def _add_global_attributes(mom6_grid, netcdf_dataset):
    netcdf_dataset.grid_version = mom6_grid['netcdf_info']['grid_version']
    netcdf_dataset.code_version = mom6_grid['netcdf_info']['code_version']
    netcdf_dataset.history      = mom6_grid['netcdf_info']['history_entry']

def write_MOM6_supergrid_file(mom6_grid):
    """Save the MOM6 supergrid data into its own file."""

    nx = mom6_grid['supergrid']['nx']
    ny = mom6_grid['supergrid']['ny']
    string_len = len(mom6_grid['netcdf_info']['tile_str'])

    with netCDF4.Dataset(mom6_grid['filenames']['supergrid'], 'w', format='NETCDF3_CLASSIC') as hgrid_ds:
        # Dimensions
        hgrid_ds.createDimension('nx',  nx)
        hgrid_ds.createDimension('nxp', nx+1)
        hgrid_ds.createDimension('ny',  ny)
        hgrid_ds.createDimension('nyp', ny+1)
        hgrid_ds.createDimension('string', string_len)

        # Variables & Values
        hx = hgrid_ds.createVariable('x', 'f4', ('nyp','nxp',))
        hy = hgrid_ds.createVariable('y', 'f4', ('nyp','nxp',))

        if 'lon' in mom6_grid['supergrid']:
            hx.units = 'degrees_east'
            hx[:] = mom6_grid['supergrid']['lon']

            hy.units = 'degrees_north'
            hy[:] = mom6_grid['supergrid']['lat']
        else:
            hx.units = 'meters'
            hx[:] = mom6_grid['supergrid']['x']

            hy.units = 'meters'
            hy[:] = mom6_grid['supergrid']['y']

        hdx = hgrid_ds.createVariable('dx', 'f4', ('nyp','nx',))
        hdx.units = 'meters'
        hdx[:] = mom6_grid['supergrid']['dx']

        hdy = hgrid_ds.createVariable('dy', 'f4', ('ny','nxp',))
        hdy.units = 'meters'
        hdy[:] = mom6_grid['supergrid']['dy']

        harea = hgrid_ds.createVariable('area', 'f4', ('ny','nx',))
        harea.units = 'meters^2'
        harea[:] = mom6_grid['supergrid']['area']

        hangle = hgrid_ds.createVariable('angle_dx', 'f4', ('nyp','nxp',))
        hangle.units = 'radians'
        hangle[:] = mom6_grid['supergrid']['angle']

        htile = hgrid_ds.createVariable('tile', 'c', ('string',))
        htile[:] = mom6_grid['netcdf_info']['tile_str']

        # Global attributes
        _add_global_attributes(mom6_grid, hgrid_ds)

def write_MOM6_topography_file(mom6_grid):
    """Save the MOM6 ocean topography field in a separate file."""

    nx = mom6_grid['cell_grid']['nx']
    ny = mom6_grid['cell_grid']['ny']
    with netCDF4.Dataset(mom6_grid['filenames']['topography'], 'w', format='NETCDF3_CLASSIC') as topog_ds:
        # Dimensions
        topog_ds.createDimension('nx', nx)
        topog_ds.createDimension('ny', ny)
        topog_ds.createDimension('ntiles', 1)

        # Variables & Values
        hdepth = topog_ds.createVariable('depth', 'f4', ('ny','nx',))
        hdepth.units = 'm'
        hdepth[:] = mom6_grid['cell_grid']['depth']

        # Global attributes
        _add_global_attributes(mom6_grid, topog_ds)

def write_MOM6_solo_mosaic_file(mom6_grid):
    """Write the "solo mosaic" file, which describes to the FMS infrastructure
     where to find the grid file(s).  Based on tools in version 5 of MOM
     (http://www.mom-ocean.org/)."""

    # NOTE: This function is very basic, since we're skipping the
    # finding of "contact regions" between the tiles that the real
    # make_solo_mosaic tool performs.  It's not needed right now,
    # since we only have one (regional) tile, but I think this feature
    # will be needed if we ever use a tripolar grid.

    with netCDF4.Dataset(mom6_grid['filenames']['mosaic'], 'w', format='NETCDF3_CLASSIC') as mosaic_ds:
        # Dimenisons
        mosaic_ds.createDimension('ntiles', 1)
        mosaic_ds.createDimension('string', mom6_grid['netcdf_info']['string_length'])

        # Variables & Values
        hmosaic = mosaic_ds.createVariable('mosaic', 'c', ('string',))
        hmosaic.standard_name = 'grid_mosaic_spec'
        hmosaic.children = 'gridtiles'
        hmosaic.contact_regions = 'contacts'
        hmosaic.grid_descriptor = ''
        dirname,filename = os.path.split(mom6_grid['filenames']['mosaic'])
        filename,ext = os.path.splitext(filename)
        hmosaic[:len(filename)] = filename

        hgridlocation = mosaic_ds.createVariable('gridlocation', 'c', ('string',))
        hgridlocation.standard_name = 'grid_file_location'
        this_dir = mom6_grid['filenames']['directory']
        hgridlocation[:len(this_dir)] = this_dir

        hgridfiles = mosaic_ds.createVariable('gridfiles', 'c', ('ntiles', 'string',))
        hgridfiles[0, :len(mom6_grid['filenames']['supergrid'])] = mom6_grid['filenames']['supergrid']

        hgridtiles = mosaic_ds.createVariable('gridtiles', 'c', ('ntiles', 'string',))
        hgridtiles[0, :len(mom6_grid['netcdf_info']['tile_str'])] = mom6_grid['netcdf_info']['tile_str']

        # Global attributes
        _add_global_attributes(mom6_grid, mosaic_ds)

def write_MOM6_land_mask_file(mom6_grid):
    """Write the land mask file.  Based on 'make_quick_mosaic' tool in version
    5 of MOM (http://www.mom-ocean.org/)."""

    nx = mom6_grid['cell_grid']['nx']
    ny = mom6_grid['cell_grid']['ny']
    with netCDF4.Dataset(mom6_grid['filenames']['land_mask'], 'w', format='NETCDF3_CLASSIC') as land_mask_ds:
        # Dimenisons (of grid cells, not supergrid)
        land_mask_ds.createDimension('nx', nx)
        land_mask_ds.createDimension('ny', ny)

        # Variables & Values
        hmask = land_mask_ds.createVariable('mask', 'd', ('ny', 'nx'))
        hmask.standard_name = 'land fraction at T-cell centers'
        hmask.units = 'none'
        hmask[:] = mom6_grid['cell_grid']['land_mask']

        # Global attributes
        _add_global_attributes(mom6_grid, land_mask_ds)

def write_MOM6_ocean_mask_file(mom6_grid):
    """Write the ocean mask file.  Based on 'make_quick_mosaic' tool in version
    5 of MOM (http://www.mom-ocean.org/)."""

    nx = mom6_grid['cell_grid']['nx']
    ny = mom6_grid['cell_grid']['ny']
    with netCDF4.Dataset(mom6_grid['filenames']['ocean_mask'], 'w', format='NETCDF3_CLASSIC') as ocean_mask_ds:
        # Dimenisons (of grid cells, not supergrid)
        ocean_mask_ds.createDimension('nx', nx)
        ocean_mask_ds.createDimension('ny', ny)

        # Variables & Values
        hmask = ocean_mask_ds.createVariable('mask', 'd', ('ny', 'nx'))
        hmask.standard_name = 'ocean fraction at T-cell centers'
        hmask.units = 'none'
        hmask[:] = mom6_grid['cell_grid']['ocean_mask']

        # Global attributes
        _add_global_attributes(mom6_grid, ocean_mask_ds)

def write_MOM6_exchange_grid_file(mom6_grid, name1, name2):
    """Write one of the three exchange grid files (depending on name1 and name2).
    Based on 'make_quick_mosaic' tool in version 5 of MOM (http://www.mom-ocean.org/)."""

    # calculate the exchange grid

    mask = mom6_grid['cell_grid'][name2 + '_mask']

    tile_cells_j, tile_cells_i = numpy.where(mask == 1)
    tile_cells = numpy.column_stack((tile_cells_i, tile_cells_j)) + 1 # +1 converts from Python indices to Fortran
    xgrid_area = mom6_grid['cell_grid']['area'][mask == 1]
    ncells = len(xgrid_area)
    tile_dist = numpy.zeros((ncells,2))

    # write out exchange grid file

    filename_key = '%s_%s_exchange' % (name1, name2)
    filename = mom6_grid['filenames'][filename_key]

    with netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC') as exchange_ds:
        exchange_ds.createDimension('string', mom6_grid['netcdf_info']['string_length'])
        exchange_ds.createDimension('ncells', ncells)
        exchange_ds.createDimension('two', 2)

        contact_str = '{0}_mosaic:{2}::{1}_mosaic:{2}'.format(name1, name2, mom6_grid['netcdf_info']['tile_str'])
        hcontact = exchange_ds.createVariable('contact', 'c', ('string',))
        hcontact.standard_name = 'grid_contact_spec'
        hcontact.contact_type = 'exchange'
        hcontact.parent1_cell = 'tile1_cell'
        hcontact.parent2_cell = 'tile2_cell'
        hcontact.xgrid_area_field = 'xgrid_area'
        hcontact.distant_to_parent1_centroid = 'tile1_distance'
        hcontact.distant_to_parent2_centroid = 'tile2_distance'
        hcontact[0:len(contact_str)] = contact_str

        htile1_cell = exchange_ds.createVariable('tile1_cell', 'i', ('ncells', 'two'))
        htile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
        htile1_cell[:] = tile_cells

        htile2_cell = exchange_ds.createVariable('tile2_cell', 'i', ('ncells', 'two'))
        htile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
        htile2_cell[:] = tile_cells

        hxgrid_area = exchange_ds.createVariable('xgrid_area', 'd', ('ncells',))
        hxgrid_area.standard_name = 'exchange_grid_area'
        hxgrid_area.units = 'm2'
        hxgrid_area[:] = xgrid_area

        htile1_dist = exchange_ds.createVariable('tile1_distance', 'd', ('ncells', 'two'))
        htile1_dist.standard_name = 'distance_from_parent1_cell_centroid'
        htile1_dist[:] = tile_dist

        htile2_dist = exchange_ds.createVariable('tile2_distance', 'd', ('ncells', 'two'))
        htile2_dist.standard_name = 'distance_from_parent2_cell_centroid'
        htile2_dist[:] = tile_dist

        # Global attributes
        _add_global_attributes(mom6_grid, exchange_ds)

def write_MOM6_coupler_mosaic_file(mom6_grid):
    """Write the compler mosaic file, which references all the rest.
    Based on 'make_quick_mosaic' tool in version 5 of MOM (http://www.mom-ocean.org/)."""

    def add_string_var_1d(fid, var_name, standard_name, value):
        """Creates and stores values for a 1-dimensional string variable."""
        id_var = fid.createVariable(var_name, 'c', ('string',))
        id_var.standard_name = standard_name
        id_var[0:len(value)] = value

    def add_string_var_2d(fid, var_name, dim_name, standard_name, value):
        """Creates and stores values for a 2-dimensional string variable, for
        which the first dimension has length 1."""
        id_var = fid.createVariable(var_name, 'c', (dim_name, 'string'))
        id_var.standard_name = standard_name
        id_var[0, 0:len(value)] = value

    with netCDF4.Dataset('mosaic.nc', 'w', format='NETCDF3_CLASSIC') as mosaic_ds:
        mosaic_ds.createDimension('string', mom6_grid['netcdf_info']['string_length'])
        mosaic_ds.createDimension('nfile_aXo', 1)
        mosaic_ds.createDimension('nfile_aXl', 1)
        mosaic_ds.createDimension('nfile_lXo', 1)

        # same mosaic file for all three -- just like when "make_solo_mosaic" is used

        add_string_var_1d(mosaic_ds, 'atm_mosaic_dir',  'directory_storing_atmosphere_mosaic', mom6_grid['filenames']['directory'])
        add_string_var_1d(mosaic_ds, 'atm_mosaic_file', 'atmosphere_mosaic_file_name',         mom6_grid['filenames']['mosaic'])
        add_string_var_1d(mosaic_ds, 'atm_mosaic',      'atmosphere_mosaic_name',              'atmos_mosaic')

        add_string_var_1d(mosaic_ds, 'lnd_mosaic_dir',  'directory_storing_land_mosaic',       mom6_grid['filenames']['directory'])
        add_string_var_1d(mosaic_ds, 'lnd_mosaic_file', 'land_mosaic_file_name',               mom6_grid['filenames']['mosaic'])
        add_string_var_1d(mosaic_ds, 'lnd_mosaic',      'land_mosaic_name',                    'land_mosaic')

        add_string_var_1d(mosaic_ds, 'ocn_mosaic_dir',  'directory_storing_ocean_mosaic',      mom6_grid['filenames']['directory'])
        add_string_var_1d(mosaic_ds, 'ocn_mosaic_file', 'ocean_mosaic_file_name',              mom6_grid['filenames']['mosaic'])
        add_string_var_1d(mosaic_ds, 'ocn_mosaic',      'ocean_mosaic_name',                   'ocean_mosaic')

        add_string_var_1d(mosaic_ds, 'ocn_topog_dir',   'directory_storing_ocean_topog',       mom6_grid['filenames']['directory'])
        add_string_var_1d(mosaic_ds, 'ocn_topog_file',  'ocean_topog_file_name',               mom6_grid['filenames']['topography'])

        add_string_var_2d(mosaic_ds, 'aXo_file', 'nfile_aXo', 'atmXocn_exchange_grid_file', mom6_grid['filenames']['atmos_ocean_exchange'])
        add_string_var_2d(mosaic_ds, 'aXl_file', 'nfile_aXl', 'atmXlnd_exchange_grid_file', mom6_grid['filenames']['atmos_land_exchange'])
        add_string_var_2d(mosaic_ds, 'lXo_file', 'nfile_lXo', 'lndXocn_exchange_grid_file', mom6_grid['filenames']['land_ocean_exchange'])

        # Global attributes
        _add_global_attributes(mom6_grid, mosaic_ds)

def main(argv):
    """Take a ROMS grid file and output a set of files to represent the MOM6 grid."""

    mom6_grid = setup_MOM6_grid(argv)

    if len(argv) == 2:
        roms_grid_filename = argv[1]
    else:
        print('Usage: %s RGRID' % argv[0])
        print('')
        print('Converts the ROMS horizontal grid stored in the NetCDF file RGRID into')
        print('a collection of NetCDF files representing the MOM6 horizontal grid:')
        print(' * supergrid file ("%s")' % mom6_grid['filenames']['supergrid'])
        print(' * topography file ("%s")' % mom6_grid['filenames']['topography'])
        print(' * land and ocean mask files ("%s", "%s")' % (mom6_grid['filenames']['land_mask'], mom6_grid['filenames']['ocean_mask']))
        print(' * coupler mosaic file ("%s")' % mom6_grid['filenames']['mosaic'])
        print(' * coupler exchange files:')
        print('    - ' + mom6_grid['filenames']['atmos_land_exchange'])
        print('    - ' + mom6_grid['filenames']['atmos_ocean_exchange'])
        print('    - ' + mom6_grid['filenames']['land_ocean_exchange'])
        print('')
        print('Files are placed in the current diretory.')
        exit(1)

    roms_grid = read_ROMS_grid(roms_grid_filename)
    roms_grid = trim_ROMS_grid(roms_grid)
    mom6_grid = convert_ROMS_to_MOM6(mom6_grid, roms_grid)
    mom6_grid = approximate_MOM6_grid_metrics(mom6_grid)

    write_MOM6_supergrid_file(mom6_grid)
    write_MOM6_topography_file(mom6_grid)
    write_MOM6_solo_mosaic_file(mom6_grid)
    write_MOM6_land_mask_file(mom6_grid)
    write_MOM6_ocean_mask_file(mom6_grid)
    write_MOM6_exchange_grid_file(mom6_grid, 'atmos',  'land')
    write_MOM6_exchange_grid_file(mom6_grid, 'atmos', 'ocean')
    write_MOM6_exchange_grid_file(mom6_grid,  'land', 'ocean')
    write_MOM6_coupler_mosaic_file(mom6_grid)

if __name__ == "__main__":
    main(sys.argv)
