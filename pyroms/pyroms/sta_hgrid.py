# encoding: utf-8
'''Tools for creating and working with Line (Station) Grids'''
__docformat__ = "restructuredtext en"

import os
import sys
import ctypes
import pickle
from warnings import warn
from copy import deepcopy

import numpy as np

import pyroms
from pyroms.vgrid import *


class Sta_CGrid(object):
    """
    Stations Grid

    EXAMPLES:
    --------

    >>> x = arange(8)
    >>> y = arange(8)*2-1
    >>> grd = pyroms.grid.StaGrid(x, y)
    >>> print grd.x_rho
    [4.5 4.5 4.5 4.5 4.5 4.5 4.5]
    """

    def __init__(self, x_rho, y_rho, angle_rho=None):

        assert np.ndim(x_rho)==1 and np.ndim(y_rho)==1 and \
            np.shape(x_rho)==np.shape(y_rho), \
            'x and y must be 2D arrays of the same size.'

        if np.any(np.isnan(x_rho)) or np.any(np.isnan(y_rho)):
            x_rho = np.ma.masked_where( (isnan(x_rho)) | (isnan(y_rho)) , x_rho)
            y_rho = np.ma.masked_where( (isnan(x_rho)) | (isnan(y_rho)) , y_rho)

        self.x_rho = x_rho
        self.y_rho = y_rho

        self.spherical = 'F'

        if angle_rho is None:
            self.angle_rho = np.zeros(len(self.y_rho))
        else:
            self.angle_rho = angle_rho


    x = property(lambda self: self.x_rho, None, None, 'Return x_rho')
    y = property(lambda self: self.y_rho, None, None, 'Return y_rho')


class Sta_CGrid_geo(Sta_CGrid):
    """
    Curvilinear Arakawa C-grid defined in geographic coordinates

    For a geographic grid, a projection may be specified, or The default
    projection for will be defined by the matplotlib.toolkits.Basemap
    projection:

    proj = Basemap(projection='merc', resolution=None, lat_ts=0.0)
    """

    def __init__(self, lon_rho, lat_rho, proj, \
                    angle_rho=None):

        x, y = proj(lon_rho, lat_rho)
        self.lon_rho = lon_rho
        self.lat_rho = lat_rho
        self.proj = proj

        #calculate cartesian position
        self.x_rho, self.y_rho = proj(lon_rho, lat_rho)
        self.x_rho, self.y_rho = proj(lon_rho, lat_rho)

        if angle_rho is None:
            self.angle_rho = np.zeros(len(self.x_rho))
        else:
            self.angle_rho = angle_rho

        self.spherical = 'T'

