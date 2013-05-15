import numpy as np
from mpl_toolkits.basemap import pyproj
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import pyroms


def make_remap_grid_file(Bgrd, Bpos='t'):

    #create remap file
    remap_filename = 'remap_grid_' + Bgrd.name + '_' + Bpos + '.nc'
    nc = netCDF.Dataset(remap_filename, 'w', format='NETCDF3_CLASSIC')
    nc.Description = 'remap grid file for POP'
    nc.Author = 'pyroms_toolbox.BGrid_POP.make_remap_grid_file'
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = Bgrd.name

    lon_corner = Bgrd.lon_t_vert
    lat_corner = Bgrd.lat_t_vert
    grid_center_lon = Bgrd.lon_t.flatten()
    grid_center_lat = Bgrd.lat_t.flatten()
    Mp, Lp = Bgrd.lon_t.shape
    if Bpos == 't':
        lon_corner = Bgrd.lon_t_vert
        lat_corner = Bgrd.lat_t_vert
        grid_center_lon = Bgrd.lon_t.flatten()
        grid_center_lat = Bgrd.lat_t.flatten()
        grid_imask = Bgrd.mask_t[0,:].flatten()
    elif Bpos == 'uv':
        lon_corner = Bgrd.lon_u_vert
        lat_corner = Bgrd.lat_u_vert
        grid_center_lon = Bgrd.lon_u.flatten()
        grid_center_lat = Bgrd.lat_u.flatten()
        grid_imask = Bgrd.mask_u[0,:].flatten()

    grid_size = Lp * Mp

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

