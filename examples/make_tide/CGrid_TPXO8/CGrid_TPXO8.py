import numpy as np
from mpl_toolkits.basemap import pyproj
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
  import pyroms

class CGrid_TPXO8(object):

    # CGrid object for TPXO8 v1

    def __init__(self, lon_t, lat_t, lon_u, lat_u, lon_v, lat_v, \
                 mask_t, mask_u, mask_v, z_t, z_u, z_v, missing_value, name, xrange, yrange):

        self.name = name

        self.missing_value = -9999.

        self.xrange = xrange
        self.yrange = yrange

        self.z_t = z_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.z_u = z_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.z_v = z_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_t = lon_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_t = lat_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_u = lon_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_u = lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lon_v = lon_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_v = lat_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_t_vert = 0.5 * (lon_t[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lon_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_t_vert = 0.5 * (lat_t[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lat_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])

        self.lon_u_vert = 0.5 * (lon_u[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lon_u[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_u_vert = 0.5 * (lat_u[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lat_u[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lon_v_vert = 0.5 * (lon_v[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lon_v[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_v_vert = 0.5 * (lat_v[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lat_v[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])

        self.mask_t = mask_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.mask_u = mask_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.mask_v = mask_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        ones = np.ones(self.z_t.shape)
        a1 = lat_u[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] - \
        lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        a2 = lon_u[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] - \
        lon_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        a3 = 0.5*(lat_u[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] + \
        lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1])
        a2 = np.where(a2 > 180*ones, a2 - 360*ones, a2)
        a2 = np.where(a2 < -180*ones, a2 + 360*ones, a2)
        a2 = a2 * np.cos(np.pi/180.*a3)
        self.angle = np.arctan2(a1, a2)
