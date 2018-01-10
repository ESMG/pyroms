import numpy as np
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import pyroms


def make_remap_grid_file(grid, Cpos='rho', irange=None, jrange=None):
    '''
    make_remap_grid_file(grid)

    generate grid file to be used with scrip to compute
    the weights for remapping.
    '''

    # get grid
    if type(grid).__name__ == 'ROMS_Grid':
        grd = grid
    else:
        grd = pyroms.grid.get_ROMS_grid(grid)

    Mp, Lp = grd.vgrid.h.shape
    if irange is None:
        irange = (0,Lp-1)
    else:
        assert irange[1]-irange[0] > 0, \
               'irange must span a positive range'

    if jrange is None:
        jrange = (0,Mp-1)
    else:
        assert jrange[1]-jrange[0] > 0, \
               'jrange must span a positive range'

    #create remap file
    remap_filename = 'remap_grid_' + grd.name + '_' + Cpos + '.nc'
    nc = netCDF.Dataset(remap_filename, 'w', format='NETCDF3_CLASSIC')
    nc.Description = 'remap grid file on' + Cpos + 'points'
    nc.Author = 'pyroms.remapping.make_remap_grid_file'
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = grd.name

    if Cpos == 'rho':
        if jrange != (0,Mp-1) or irange != (0,Lp-1):
            lon_corner = grd.hgrid.lon_vert[jrange[0]:jrange[1]+1, \
                  irange[0]:irange[1]+1]
            lat_corner = grd.hgrid.lat_vert[jrange[0]:jrange[1]+1, \
                  irange[0]:irange[1]+1]
            grid_center_lon = grd.hgrid.lon_rho[jrange[0]:jrange[1], \
                  irange[0]:irange[1]].flatten()
            grid_center_lat = grd.hgrid.lat_rho[jrange[0]:jrange[1], \
                  irange[0]:irange[1]].flatten()
            grid_imask = grd.hgrid.mask_rho[jrange[0]:jrange[1], \
                  irange[0]:irange[1]].flatten()
            Lp = irange[1] - irange[0]
            Mp = jrange[1] - jrange[0]
        else:
            lon_corner = grd.hgrid.lon_vert
            lat_corner = grd.hgrid.lat_vert
            grid_center_lon = grd.hgrid.lon_rho.flatten()
            grid_center_lat = grd.hgrid.lat_rho.flatten()
            grid_imask = grd.hgrid.mask_rho.flatten()
            Mp, Lp = grd.hgrid.mask_rho.shape
    elif Cpos == 'u':
        if jrange != (0,Mp-1) or irange != (0,Lp-1):
            lon_corner = 0.5 * \
                  (grd.hgrid.lon_vert[jrange[0]:jrange[1]+1,irange[0]:irange[1]] + \
                   grd.hgrid.lon_vert[jrange[0]:jrange[1]+1,1+irange[0]:irange[1]+1])
            lat_corner = 0.5 * \
                  (grd.hgrid.lat_vert[jrange[0]:jrange[1]+1,irange[0]:irange[1]] + \
                   grd.hgrid.lat_vert[jrange[0]:jrange[1]+1,1+irange[0]:irange[1]+1])
            grid_center_lon = grd.hgrid.lon_u[jrange[0]:jrange[1], \
                   irange[0]:irange[1]-1].flatten()
            grid_center_lat = grd.hgrid.lat_u[jrange[0]:jrange[1], \
                   irange[0]:irange[1]-1].flatten()
            grid_imask = grd.hgrid.mask_u[jrange[0]:jrange[1], \
                   irange[0]:irange[1]-1].flatten()
            Lp = irange[1] - irange[0] - 1
            Mp = jrange[1] - jrange[0]
        else:
            lon_corner = 0.5 * (grd.hgrid.lon_vert[:,:-1] + \
                                grd.hgrid.lon_vert[:,1:])
            lat_corner = 0.5 * (grd.hgrid.lat_vert[:,:-1] + \
                                grd.hgrid.lat_vert[:,1:])
            grid_center_lon = grd.hgrid.lon_u.flatten()
            grid_center_lat = grd.hgrid.lat_u.flatten()
            grid_imask = grd.hgrid.mask_u.flatten()
            Mp, Lp = grd.hgrid.mask_u.shape
    elif Cpos == 'v':
        if jrange != (0,Mp-1) or irange != (0,Lp-1):
            lon_corner = 0.5 * \
                  (grd.hgrid.lon_vert[jrange[0]:jrange[1],irange[0]:irange[1]+1] + \
                   grd.hgrid.lon_vert[1+jrange[0]:jrange[1]+1,irange[0]:irange[1]+1])
            lat_corner = 0.5 * \
                  (grd.hgrid.lat_vert[jrange[0]:jrange[1],irange[0]:irange[1]+1] + \
                   grd.hgrid.lat_vert[1+jrange[0]:jrange[1]+1,irange[0]:irange[1]+1])
            grid_center_lon = grd.hgrid.lon_v[jrange[0]:jrange[1]-1, \
                                irange[0]:irange[1]].flatten()
            grid_center_lat = grd.hgrid.lat_v[jrange[0]:jrange[1]-1, \
                                irange[0]:irange[1]].flatten()
            grid_imask = grd.hgrid.mask_v[jrange[0]:jrange[1]-1, \
                                irange[0]:irange[1]].flatten()
            Lp = irange[1] - irange[0]
            Mp = jrange[1] - jrange[0] - 1
        else:
            lon_corner = 0.5 * (grd.hgrid.lon_vert[:-1,:] + \
                                grd.hgrid.lon_vert[1:,:])
            lat_corner = 0.5 * (grd.hgrid.lat_vert[:-1,:] + \
                                grd.hgrid.lat_vert[1:,:])
            grid_center_lon = grd.hgrid.lon_v.flatten()
            grid_center_lat = grd.hgrid.lat_v.flatten()
            grid_imask = grd.hgrid.mask_v.flatten()
            Mp, Lp = grd.hgrid.mask_v.shape
    else:
        raise ValueError('Cpos must be rho, u or v')

    grid_size = Lp * Mp
    print('grid shape', Mp, Lp)

    grid_corner_lon = np.zeros((grid_size, 4))
    grid_corner_lat = np.zeros((grid_size, 4))
    k = 0
    for j in range(Mp):
        for i in range(Lp):
            grid_corner_lon[k,0] = lon_corner[j,i]
            grid_corner_lat[k,0] = lat_corner[j,i]
            grid_corner_lon[k,1] = lon_corner[j,i+1]
            grid_corner_lat[k,1] = lat_corner[j,i+1]
            grid_corner_lon[k,2] = lon_corner[j+1,i+1]
            grid_corner_lat[k,2] = lat_corner[j+1,i+1]
            grid_corner_lon[k,3] = lon_corner[j+1,i]
            grid_corner_lat[k,3] = lat_corner[j+1,i]
            k = k + 1

    #Write netcdf file
    nc.createDimension('grid_size', grid_size)
    nc.createDimension('grid_corners', 4)
    nc.createDimension('grid_rank', 2)

    nc.createVariable('grid_dims', 'i4', ('grid_rank'))
    nc.variables['grid_dims'].long_name = 'grid size along x and y axis'
    nc.variables['grid_dims'].units = 'None'
    nc.variables['grid_dims'][:] = [(Lp, Mp)]

    nc.createVariable('grid_center_lon', 'f8', ('grid_size'))
    nc.variables['grid_center_lon'].long_name = 'longitude of cell center'
    nc.variables['grid_center_lon'].units = 'degrees'
    nc.variables['grid_center_lon'][:] = grid_center_lon

    nc.createVariable('grid_center_lat', 'f8', ('grid_size'))
    nc.variables['grid_center_lat'].long_name = 'latitude of cell center'
    nc.variables['grid_center_lat'].units = 'degrees'
    nc.variables['grid_center_lat'][:] = grid_center_lat

    nc.createVariable('grid_imask', 'i4', ('grid_size'))
    nc.variables['grid_imask'].long_name = 'mask'
    nc.variables['grid_imask'].units = 'None'
    nc.variables['grid_imask'][:] = grid_imask

    nc.createVariable('grid_corner_lon', 'f8', ('grid_size', 'grid_corners'))
    nc.variables['grid_corner_lon'].long_name = 'longitude of cell corner'
    nc.variables['grid_corner_lon'].units = 'degrees'
    nc.variables['grid_corner_lon'][:] = grid_corner_lon

    nc.createVariable('grid_corner_lat', 'f8', ('grid_size', 'grid_corners'))
    nc.variables['grid_corner_lat'].long_name = 'latitude of cell corner'
    nc.variables['grid_corner_lat'].units = 'degrees'
    nc.variables['grid_corner_lat'][:] = grid_corner_lat

    nc.close()
