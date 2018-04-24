import numpy as np
from datetime import datetime
import netCDF4 as netCDF

import pyroms
import pyroms_toolbox


##  load 2-dimentional interannual discharge data 
##  from 1948-2007. See Dai and Trenberth (2002) and Dai et al. (2009)
print('Load interannual discharge data')
nc_data = netCDF.Dataset('/archive/u1/uaf/kate/CORE2/runoff.daitren.iaf.10FEB2011.nc', 'r')
runoff = nc_data.variables['runoff'][:]
lon = nc_data.variables['xc'][:]
lat = nc_data.variables['yc'][:]
lon_corner = nc_data.variables['xv'][:]
lat_corner = nc_data.variables['yv'][:]
mask = nc_data.variables['mask'][:]

##  create data remap file for scrip
print('Create remap grid file for Dai and Trenberth runoff')
remap_filename = 'remap_grid_daitren.nc'
nc = netCDF.Dataset(remap_filename, 'w', format='NETCDF3_CLASSIC')
nc.Description = 'remap grid file for Dai and Trenberth runoff data'
nc.Author = 'build_runoff'
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'Dai and Trenberth runoff'

grid_center_lon = lon.flatten()
grid_center_lat = lat.flatten()
Mp, Lp = lon.shape
grid_imask = mask.flatten()
grid_size = Lp * Mp
grid_corner_lon = np.zeros((grid_size, 4))
grid_corner_lat = np.zeros((grid_size, 4))
k = 0
for j in range(Mp):
    for i in range(Lp):
        grid_corner_lon[k,0] = lon_corner[0,j,i]
        grid_corner_lat[k,0] = lat_corner[0,j,i]
        grid_corner_lon[k,1] = lon_corner[1,j,i]
        grid_corner_lat[k,1] = lat_corner[1,j,i]
        grid_corner_lon[k,2] = lon_corner[2,j,i]
        grid_corner_lat[k,2] = lat_corner[2,j,i]
        grid_corner_lon[k,3] = lon_corner[3,j,i]
        grid_corner_lat[k,3] = lat_corner[3,j,i]
        k = k + 1

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


##  create SCS remap file for scrip
print('Create remap grid file for SCS grid')
dstgrd = pyroms.grid.get_ROMS_grid('SCS')
dstgrd.hgrid.mask_rho = np.ones(dstgrd.hgrid.mask_rho.shape)
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')


## compute remap weights
print('compute remap weights using scrip')
# input namelist variables for conservative remapping at rho points
grid1_file = 'remap_grid_daitren.nc'
grid2_file = 'remap_grid_SCS_rho.nc'
interp_file1 = 'remap_weights_daitren_to_SCS_conservative_nomask.nc'
interp_file2 = 'remap_weights_SCS_to_daitren_conservative_nomask.nc'
map1_name = 'daitren to SCS conservative Mapping'
map2_name = 'SCS to daitren conservative Mapping'
num_maps = 1
map_method = 'conservative'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)
