import os
from pyroms import _iso
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy.interpolate import griddata
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4

import pyroms
from bathy_smoother import *

# Grid dimension
Lm = 140
Mm = 120

lon0=117.5 ; lat0 = 41.
lon1=117.5 ; lat1 = 34.5
lon2 = 127. ; lat2 = 34.5
lon3 = 127. ; lat3 = 41.
map = Basemap(projection='lcc', lat_0=35., lat_1=30., lat_2=40, lon_0 =123, \
              width=2000000, height=2000000, resolution='i')

lonp = np.array([lon0, lon1, lon2, lon3])
latp = np.array([lat0, lat1, lat2, lat3])
beta = np.array([1, 1, 1, 1])

#generate the new grid
# Do this if you aren't going to move the grid corners interactively.
hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)
# Do this if you are going to use the Boundary Interactor
#map.drawcoastlines()
#xp, yp = map(lonp, latp)
#bry = pyroms.hgrid.BoundaryInteractor(xp, yp, beta, shp=(Mm+3,Lm+3), proj=map)
#hgrd=bry.grd

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
    hgrd.mask_polygon(vv,mask_value=0)

# Edit the land mask interactively.
#pyroms.grid.edit_mask_mesh(hgrd, proj=map)
#edit_mask_mesh_ij is a faster version using imshow... but no map projection.
coast = pyroms.utility.get_coast_from_map(map)
pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)


#### Use the following to interpolate from etopo2 bathymetry.
# generate the bathy
# read in topo data (on a regular lat/lon grid)
# this topo come with basemap so you should have it on your laptop.
# just update datadir with the appropriate path
# you can get this data from matplolib svn with
# svn co https://matplotlib.svn.sourceforge.net/svnroot/matplotlib/trunk/htdocs/screenshots/data/"

datadir = 'data/'
topo = np.loadtxt(os.path.join(datadir, 'etopo20data.gz'))
lons = np.loadtxt(os.path.join(datadir, 'etopo20lons.gz'))
lats = np.loadtxt(os.path.join(datadir, 'etopo20lats.gz'))

# depth positive
topo = -topo

# fix minimum depth
hmin = 5
topo = np.where(topo < hmin, hmin, topo)

# interpolate new bathymetry
lon, lat = np.meshgrid(lons, lats)
h = griddata((lon.flat,lat.flat),topo.flat,(hgrd.lon_rho,hgrd.lat_rho), method='linear')

# insure that depth is always deeper than hmin
h = np.where(h < hmin, hmin, h)

# set depth to hmin where masked
idx = np.where(hgrd.mask_rho == 0)
h[idx] = hmin

# save raw bathymetry
hraw = h.copy()

# check bathymetry roughness
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
rx0_max = 0.35
h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

# check bathymetry roughness again
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# vertical coordinate
theta_b = 2
theta_s = 7.0
Tcline = 50
N = 30
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)

# ROMS grid
grd_name = 'YELLOW'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# write grid to netcdf file
pyroms.grid.write_ROMS_grid(grd, filename='YELLOW_grd_v1.nc')
