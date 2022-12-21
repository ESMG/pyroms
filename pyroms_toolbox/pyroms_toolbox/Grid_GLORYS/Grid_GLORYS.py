import numpy as np
from mpl_toolkits.basemap import pyproj
from datetime import datetime

class Grid_GLORYS(object):
    """
    A-Grid object for GLORYS
    """

    def __init__(self, lon_t, lat_t, mask_t, depth, depth_bnds, h, name, irange, jrange):

        self.name = name

        self.irange = irange
        self.jrange = jrange

        self.h = h[jrange[0]:jrange[1]+1, irange[0]:irange[1]+1]

        self.lon_t = lon_t[irange[0]:irange[1]+1]
        self.lat_t = lat_t[jrange[0]:jrange[1]+1]

        self.lon_t_vert = 0.5 * (lon_t[irange[0]-1:irange[1]+1] + lon_t[irange[0]:irange[1]+2])
        self.lat_t_vert = 0.5 * (lat_t[jrange[0]-1:jrange[1]+1] + lat_t[jrange[0]:jrange[1]+2])

        self.mask_t = mask_t[:, jrange[0]:jrange[1]+1, irange[0]:irange[1]+1]

        self.z_t = np.tile(depth,(self.mask_t.shape[2],self.mask_t.shape[1],1)).T

        self.z_t_bnds = np.tile(depth_bnds,(self.mask_t.shape[2],self.mask_t.shape[1],1)).T

        self.angle = np.zeros(self.h.shape)
