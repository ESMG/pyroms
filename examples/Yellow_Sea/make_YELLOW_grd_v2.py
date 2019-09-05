import numpy as np
import netCDF4 as netCDF
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime
from scipy.signal import medfilt2d
import pyroms
import pyroms_toolbox
from bathy_smoother import *


# get horizontal grid from v1
grd_v1 = pyroms.grid.get_ROMS_grid('YELLOW')
hgrd = grd_v1.hgrid

# generate the bathy
data = netCDF.Dataset('/Volumes/R1/Data/Bathymetries/SmithSandwell_30s/topo30.swap_YELLOW.grd', 'r')
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
lon_corner = 0.5 * (lon[1:] + lon[:-1])
lat_corner = 0.5 * (lat[1:] + lat[:-1])
lon = lon[1:-1]
lat = lat[1:-1]
z = data.variables['z'][1:-1,1:-1]
z = np.array(z, dtype='float32')

# median filter
z = medfilt2d(z, 5)

# interpolate bathymetry
lon, lat = np.meshgrid(lon, lat)
lon_corner, lat_corner = np.meshgrid(lon_corner, lat_corner)
mask = np.ones(lon.shape)

# create SRTM remap file for scrip
remap_filename = 'remap_grid_SRTM.nc'
nc = netCDF.Dataset(remap_filename, 'w', format='NETCDF3_64BIT')
nc.Description = 'remap grid file for SRTM 30s bathymetry'
nc.Author = 'make_YELLOW_grd_v2.py'
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'SRTM 30s bathymetry'

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
        grid_corner_lon[k,0] = lon_corner[j,i]
        grid_corner_lat[k,0] = lat_corner[j,i]
        grid_corner_lon[k,1] = lon_corner[j,i+1]
        grid_corner_lat[k,1] = lat_corner[j,i+1]
        grid_corner_lon[k,2] = lon_corner[j+1,i+1]
        grid_corner_lat[k,2] = lat_corner[j+1,i+1]
        grid_corner_lon[k,3] = lon_corner[j+1,i]
        grid_corner_lat[k,3] = lat_corner[j+1,i]
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

# create coral remap file for scrip
dstgrd = pyroms.grid.get_ROMS_grid('YELLOW')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_SRTM.nc'
grid2_file = 'remap_grid_YELLOW_rho.nc'
interp_file1 = 'remap_weights_SRTM_to_YELLOW_bilinear.nc'
interp_file2 = 'remap_weights_YELLOW_to_SRTM_bilinear.nc'
map1_name = 'SRTM to YELLOW bilinear Mapping'
map2_name = 'YELLOW to SRTM bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

# remap bathymetry using scrip
h = pyroms.remapping.remap(z, 'remap_weights_SRTM_to_YELLOW_bilinear.nc', \
                           spval=1e37)
h = -h
hmin = 5
h = np.where(h < hmin, hmin, h)

# save raw bathymetry
hraw = h.copy()

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())
hsmooth = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, 0.3)
RoughMat = bathy_tools.RoughnessMatrix(hsmooth, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# insure that depth is always deeper than hmin
h = np.where(h < hmin, hmin, h)

# set depth to hmin where masked
idx = np.where(hgrd.mask_rho == 0)
h[idx] = hmin

# vertical grd
theta_b = 0.4
theta_s = 5.0
Tcline = 5
N = 30
vgrd = pyroms.vgrid.s_coordinate(hsmooth, theta_b, theta_s, Tcline, N, hraw=hraw)

# ROMS grid
grd_name = 'YELLOW'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)


# write grid to netcdf file
pyroms.grid.write_ROMS_grid(grd, 'YELLOW_grd_v2.nc')
