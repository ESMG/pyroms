import numpy as np
from mpl_toolkits.basemap import pyproj
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
  import pyroms

class SODA_Grid(object):
    """
    Grid object for SODA
    """

    def __init__(self, lon, lat, mask, depth, depth_bnds, h, name, xrange, yrange):

        self.name = name

        self.xrange = xrange
        self.yrange = yrange

        self.h = h[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        lon, lat = np.meshgrid(lon, lat)

        self.lon = lon[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat = lat[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_vert = 0.5 * (lon[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lon[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_vert = 0.5 * (lat[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lat[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])

        self.mask = mask[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.z = np.tile(depth,(self.mask.shape[2],self.mask.shape[1],1)).T

        self.z_bnds = np.tile(depth_bnds,(self.mask.shape[2],self.mask.shape[1],1)).T
