# This file is designed to be cut and pasted into an ipython --pylab
# session. Otherwise, you'll need to "import numpy as np" then
# convert "array" to "np.array".
import os
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4

import pyroms
import pyroms_toolbox
from bathy_smoother import *


#Grid dimension
Lp = 150
Mp = 80

## define the 4 corners of the grid
## first point is the top left corner then counter clock wise rotation
#lon0=-79.59 ; lat0=36.95
#lon1=-75.84 ; lat1=33.75
#lon2=-67.70 ; lat2=39.89
#lon3=-71.43 ; lat3=42.87

# define the 4 corners of the grid
# first point is the top left corner then counter clock wise rotation
lon0=98. ; lat0=25.
lon1=98. ; lat1=-10.
lon2=125. ; lat2=-10.
lon3=125. ; lat3=25.

#define map projection (here mercator)
lon_min = min(lon0,lon1,lon2,lon3)
lon_max = max(lon0,lon1,lon2,lon3)
lon_0 = (lon_min + lon_max) / 2.

lat_min = min(lat0,lat1,lat2,lat3)
lat_max = max(lat0,lat1,lat2,lat3)
lat_0 = (lat_min + lat_max) / 2.

map = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min, \
         urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0, \
         resolution='i')

#generate the new grid
lonp=array([lon0, lon1, lon2, lon3])
latp=array([lat0, lat1, lat2, lat3])

beta = array([1, 1, 1, 1])

hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mp+3,Lp+3), proj=map)

# if you want to us ethe graphical interface
#map.drawcoastlines()
#xp, yp = map(lonp, latp)
#bry = pyroms.hgrid.BoundaryInteractor(xp, yp, beta, shp=(Mp+3,Lp+3), proj=map)
#hgrd = bry.grd


lonv, latv = list(map(hgrd.x_vert, hgrd.y_vert, inverse=True))
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)

# generate the mask
#for verts in map.coastsegs:
#    hgrd.mask_polygon(verts)
# alternate version from johan.navarro.padron

for xx,yy in map.coastpolygons:
    xa = np.array(xx, np.float32)
    ya = np.array(yy,np.float32)
    vv = np.zeros((xa.shape[0],2))
    vv[:, 0] = xa
    vv[:, 1] = ya

# Edit the mask for change
pyroms.grid.edit_mask_mesh(hgrd, proj=map)

# generate the bathy
# read in topo data (on a regular lat/lon grid)
# this topo come with basemap so you should have it on your laptop.
# just update datadir with the appropriate path
# you can get this data from matplolib svn with
# svn co https://matplotlib.svn.sourceforge.net/svnroot/matplotlib/trunk/htdocs/screenshots/data/"
datadir = '/home/frederic/python/basemap-0.99.4/examples/'
topo = np.loadtxt(os.path.join(datadir, 'etopo20data.gz'))
lons = np.loadtxt(os.path.join(datadir, 'etopo20lons.gz'))
lats = np.loadtxt(os.path.join(datadir, 'etopo20lats.gz'))

# shift data so lons go from -180 to 180 instead of 20 to 380.
topo,lons = shiftgrid(180.,topo,lons,start=False)

# keep only the US
topo = topo[270:525, 30:400]
lons = lons[30:400]
lats = lats[270:525]

# depth positive
topo = -topo

# fix minimum depth
hmin = 5
topo = pyroms_toolbox.change(topo, '<', hmin, hmin)

# interpolate new bathymetry
lon, lat = meshgrid(lons, lats)
h = griddata(lon.flat,lat.flat,topo.flat,hgrd.lon_rho,hgrd.lat_rho)

# insure that depth is always deeper than hmin
h = pyroms_toolbox.change(h, '<', hmin, hmin)

# check bathymetry roughness
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
rx0_max = 0.35
h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

# check bathymetry roughness again
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

hgrd.h = h

# define vertical grd
hc = 5.0
theta_b = 0.4
theta_s = 5.0
Tcline = 5
N = 36
vgrd = pyroms.vgrid.s_coordinate(h, hc, theta_b, theta_s, Tcline, N)

#ROMS grid
grd_name = 'test'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

#write grid to netcdf file
pyroms.grid.write_ROMS_grid(grd, filename='grd.nc')






