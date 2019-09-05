import numpy as np
import netCDF4
import sys
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
from scipy.signal import medfilt2d
from scipy.interpolate import interp2d
#from scipy.interpolate import griddata

import pyroms
import pyroms_toolbox
from bathy_smoother import *
import creep

# For converting a WRF grid to a ROMS grid.
# WRF is more restrictive in the grids it can handle, plus
# has a nice grid generation tool. This way, you can have
# grids that match for both models.
#
# Only supports lcc grids so far - I'd welcome a pull request
# for more options.
#
fname = "geo_em.d01.nc"
bathyfile = '/import/AKWATERS/kshedstrom/bathy/ARDEMv2.0.nc'
f2name = "Bering_WRF_grid.nc"
ncid = netCDF4.Dataset(fname, "r")

rlat = ncid.variables['XLAT_M'][0,:,:]
rlon = ncid.variables['XLONG_M'][0,:,:]
latv = ncid.variables['XLAT_V'][0,1:-1,:]
lonv = ncid.variables['XLONG_V'][0,1:-1,:]
latu = ncid.variables['XLAT_U'][0,:,1:-1]
lonu = ncid.variables['XLONG_U'][0,:,1:-1]
f = ncid.variables['F'][:]
cosang = ncid.variables['COSALPHA'][:]
sinang = ncid.variables['SINALPHA'][:]
mask = ncid.variables['LANDMASK'][:]
f = ncid.variables['F'][:]
corner_lats = ncid.corner_lats
corner_lons = ncid.corner_lons

map_proj = ncid.MAP_PROJ
if (map_proj == 1):
   my_proj = 'lcc'
   lat_1 = ncid.TRUELAT1
   lat_2 = ncid.TRUELAT2
   lon_0 = ncid.STAND_LON
   llcrnrlon = corner_lons[12]
   llcrnrlat = corner_lats[12]
   urcrnrlon = corner_lons[14]
   urcrnrlat = corner_lats[14]

ncid.close()

map = Basemap(projection=my_proj, lat_1=lat_1, lat_2=lat_2, lon_0=lon_0, \
llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, \
resolution='h')

# define the 4 corners of the grid
# first point is the top left corner then counter clock wise rotation
lon0=corner_lons[13] ; lat0=corner_lats[13]
lon1=corner_lons[12] ; lat1=corner_lats[12]
lon2=corner_lons[15] ; lat2=corner_lats[15]
lon3=corner_lons[14] ; lat3=corner_lats[14]

#generate the new grid
lonp=np.array([lon0, lon1, lon2, lon3])
latp=np.array([lat0, lat1, lat2, lat3])

# shift data so lons go from 0 to 360 instead of -180 to 180.
lonp = np.where(lonp < 0, lonp+360, lonp)

beta = np.array([1, 1, 1, 1])

Mp, Lp  = rlon.shape
hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mp+1,Lp+1), proj=map)

lonv, latv = map(hgrd.x_vert, hgrd.y_vert, inverse=True)
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)
hgrd.lon_rho = np.where(hgrd.lon_rho < 0, hgrd.lon_rho+360, hgrd.lon_rho)

hgrd.mask_rho = (1.0 - mask)

ncid = netCDF4.Dataset(bathyfile, "r")

lons = ncid.variables['lon'][:]
lats = ncid.variables['lat'][:]
topo = ncid.variables['z'][:]
ncid.close()

# depth positive
topo = -topo

# fix minimum depth
hmin = 5
topo = np.where(topo < hmin, hmin, topo)

# interpolate new bathymetry
lon, lat = meshgrid(lons, lats)
h = griddata(lon.flat,lat.flat,topo.flat,hgrd.lon_rho,hgrd.lat_rho)
#h = griddata((lon,lat),topo.flat,(hgrd.lon_rho,hgrd.lat_rho),method='nearest')
#interpf = interp2d(lons,lats,topo)
#lon_rho = reshape(hgrd.lon_rho,Lp*Mp)
#lat_rho = reshape(hgrd.lat_rho,Lp*Mp)
#h = interpf(lon_rho,lat_rho)

# insure that depth is always deeper than hmin
h = np.where(h < hmin, hmin, h)

# check bathymetry roughness
hgrd.mask_rho = reshape(hgrd.mask_rho, (Mp,Lp))
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()
#h = creep.cslf(h, nan, -200., 200.)
h = np.where(isnan(h), 5500.0, h)
hraw = h.copy()

# smooth the raw bathy using the direct iterative method from Martinho and Batteen
# (2006)
rx0_max = 0.35
h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

# check bathymetry roughness again
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()
h = pyroms_toolbox.shapiro_filter.shapiro2(h, 2)

hgrd.h = h
# define vertical grd
Vtrans = 2
theta_s = 7.0
theta_b = 0.1
Tcline = 250
N = 50
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)
vgrd.h = h   # what the hell??

#ROMS grid
grd_name = 'Bering_WRF'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

#write grid to netcdf file
pyroms.grid.write_ROMS_grid(grd, filename=f2name)

