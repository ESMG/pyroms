import numpy as np
from mpl_toolkits.basemap import pyproj
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import pyroms

class BGrid_SODA(object):
    """
    BGrid object for SODA
    """

    def __init__(self, lon_t, lat_t, lon_uv, lat_uv, mask_t, mask_uv, depth, depth_bnds, h, name, xrange, yrange):

        self.name = name

        self.xrange = xrange
        self.yrange = yrange

        self.h = h[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        lon_t, lat_t = np.meshgrid(lon_t, lat_t)
        lon_uv, lat_uv = np.meshgrid(lon_uv, lat_uv)

        self.lon_t = lon_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_t = lat_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_uv = lon_uv[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_uv = lat_uv[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_t_vert = 0.5 * (lon_t[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lon_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_t_vert = 0.5 * (lat_t[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lat_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])

        self.lon_uv_vert = 0.5 * (lon_uv[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lon_uv[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_uv_vert = 0.5 * (lat_uv[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lat_uv[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])

        self.mask_t = mask_t[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.mask_uv = mask_uv[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.z_t = np.tile(depth,(self.mask_t.shape[2],self.mask_t.shape[1],1)).T

        self.z_t_bnds = np.tile(depth_bnds,(self.mask_t.shape[2],self.mask_t.shape[1],1)).T
