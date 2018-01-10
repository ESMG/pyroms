import numpy as np
from mpl_toolkits.basemap import pyproj
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF



class BGrid_GFDL(object):
    """
    Arakawa B-Grid for GFDL CM2.1
    """


    def __init__(self, lon_t, lat_t, lon_uv, lat_uv, \
                       mask_t, mask_uv, h, z_t, z_t_edges, \
                       z_uv, z_uv_edges, f, name, xrange, yrange):

        self.name = name

        self.xrange = xrange
        self.yrange = yrange

        self.lon_t = lon_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_t = lat_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lon_uv = lon_uv[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_uv = lat_uv[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.mask_t = mask_t[:,yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.mask_uv = mask_uv[:,yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.h = h[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.z_t = z_t
        self.z_t_edges = z_t_edges
        self.z_uv = z_uv
        self.z_uv_edges = z_uv_edges
        self.f = f[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_t_vert = lon_uv[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1]
        self.lat_t_vert = lat_uv[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1]

        self.lon_uv_vert = lon_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2]
        self.lat_uv_vert = lat_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2]

        self._calculate_grid_angle()


    def _calculate_grid_angle(self):
        geod = pyproj.Geod(ellps='WGS84')
        az_forward, az_back, dx = geod.inv(self.lon_t_vert[:,:-1], self.lat_t_vert[:,:-1], \
                                           self.lon_t_vert[:,1:], self.lat_t_vert[:,1:])

        angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
        self.angle = (90 - angle) * np.pi/180.
