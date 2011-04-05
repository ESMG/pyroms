import numpy as np
from mpl_toolkits.basemap import pyproj
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import pyroms



class BGrid_POP(object):
    """
    Arakawa B-Grid for POP
    """


    def __init__(self, lon_t, lat_t, lon_uv, lat_uv, \
                       mask_t, z_t, angle, name):

        self.name = name

        self.xrange = xrange
        self.yrange = yrange

        self.lon_t = lon_t
        self.lat_t = lat_t
        self.lon_uv = lon_uv
        self.lat_uv = lat_uv
        self.mask_t = mask_t
#        self.mask_uv = mask_uv
#        self.h = h
        self.z_t = z_t
#        self.z_t_edges = z_t_edges
#        self.z_uv = z_uv
#        self.z_uv_edges = z_uv_edges
#        self.f = f

        self.lon_t_vert = lon_uv
        self.lat_t_vert = lat_uv

        self.lon_uv_vert = lon_t
        self.lat_uv_vert = lat_t

        self.angle = angle


    def _calculate_grid_angle(self):
        geod = pyproj.Geod(ellps='WGS84')
        az_forward, az_back, dx = geod.inv(self.lon_t_vert[:,:-1], self.lat_t_vert[:,:-1], \
                                           self.lon_t_vert[:,1:], self.lat_t_vert[:,1:])

        angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
        self.angle = (90 - angle) * np.pi/180.
