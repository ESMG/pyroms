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
    BGrid object for POP
    """

    def __init__(self, lon_t, lat_t, lon_u, lat_u, angle, h_t, h_u, z_t, z_w, name, xrange, yrange):

        self.name = name

        self.xrange = xrange
        self.yrange = yrange

        self.lon_t = lon_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_t = lat_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_u = lon_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_u = lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.angle = angle[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.h_t = h_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1] / 100.

        self.h_u = h_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1] / 100.

        self.z_t = np.tile(z_t,(self.h_t.shape[1],self.h_t.shape[0],1)).T / 100.

        self.z_w = np.tile(z_w,(self.h_t.shape[1],self.h_t.shape[0],1)).T / 100.

        self.lon_t_vert = lon_u[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1]
        self.lat_t_vert = lat_u[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1]

        self.lon_u_vert = lon_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2]
        self.lat_u_vert = lat_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2]

        # compute the mask at t point from h_t
        self.mask_t = np.ones(self.z_t.shape)
        for n in range(self.z_t.shape[0]):
            depth = self.z_w[n,0,0]
            rtol=1e-6
            midx = np.where(np.abs(self.h_t - depth) <= rtol * np.abs(depth))
            self.mask_t[n:,midx[0],midx[1]] = 0.

        # compute the mask at u point from h_u
        self.mask_u = np.ones(self.z_t.shape)
        for n in range(self.z_t.shape[0]):
            depth = self.z_w[n,0,0]
            rtol=1e-6
            midx = np.where(np.abs(self.h_u - depth) <= rtol * np.abs(depth))
            self.mask_u[n:,midx[0],midx[1]] = 0.

