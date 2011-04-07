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
                       mask_t, mask_uv, z_t, angle, name):

        self.name = name

        self.lon_t = lon_t
        self.lat_t = lat_t
        self.lon_uv = lon_uv
        self.lat_uv = lat_uv
        self.mask_t = mask_t
        self.mask_uv = mask_uv
#        self.h = h
        self.z_t = z_t
#        self.z_t_edges = z_t_edges
#        self.z_uv = z_uv
#        self.z_uv_edges = z_uv_edges
#        self.f = f

        self._calculate_t_vert()
        self._calculate_uv_vert()

        self.angle = angle

    def _calculate_t_vert(self):
        Mm, Lm = self.lon_t.shape
        lon = np.zeros((Mm+1,Lm+1))
        lat = np.zeros((Mm+1,Lm+1))

        lon[1:, 1:] = self.lon_uv[:,:]
        lat[1:, 1:] = self.lat_uv[:,:]

        #South edge
        lon[0,0:-1] = self.lon_t[0,:] - ( self.lon_uv[0,:] - self.lon_t[0,:] )
        lon[0,-1] = self.lon_t[0,-1] - ( self.lon_uv[0,-2] - self.lon_t[0,-1] )
        lat[0,0:-1] = self.lat_t[0,:] - ( self.lat_uv[0,:] - self.lat_t[0,:] )
        lat[0,-1] = self.lat_t[0,-1] - ( self.lat_uv[0,-2] - self.lat_t[0,-1] )

        #West edge
        lon[0:-1,0] = self.lon_t[:,0] - ( self.lon_uv[:,0] - self.lon_t[:,0] )
        lon[-1,0] = self.lon_t[-1,0] - ( self.lon_uv[-2,0] - self.lon_t[-1,0] )
        lat[0:-1,0] = self.lat_t[:,0] - ( self.lat_uv[:,0] - self.lat_t[:,0] )
        lat[-1,0] = self.lat_t[-1,0] - ( self.lat_uv[-2,0] - self.lat_t[-1,0] )

        self.lon_t_vert = lon
	self.lat_t_vert = lat


    def _calculate_uv_vert(self):
        Mm, Lm = self.lon_uv.shape
        lon = np.zeros((Mm+1,Lm+1))
        lat = np.zeros((Mm+1,Lm+1))

        lon[:-1, :-1] = self.lon_t[:,:]
        lat[:-1, :-1] = self.lat_t[:,:]

        #North edge
        lon[-1,0:-2] = self.lon_uv[-1,:-1] - ( self.lon_t[-1,1:] - self.lon_uv[-1,:-1] )
        lon[-1,-2:] = self.lon_uv[-1,-2:] - ( self.lon_t[-1,-2:] - self.lon_uv[-1,-2:] )
        lat[-1,0:-2] = self.lat_uv[-1,:-1] - ( self.lat_t[-1,1:] - self.lat_uv[-1,:-1] )
        lat[-1,-2:] = self.lat_uv[-1,-2:] - ( self.lat_t[-1,-2:] - self.lat_uv[-1,-2:] )

        #East edge
        lon[0:-2,-1] = self.lon_uv[:-1,-1] - ( self.lon_t[1:,-1] - self.lon_uv[:-1,-1] )
        lon[-2,-1] = self.lon_uv[-2:-1,-1] - ( self.lon_t[-2:-1,-1] - self.lon_uv[-2:-1,-1] )
        lat[0:-2,-1] = self.lat_uv[:-1,-1] - ( self.lat_t[1:,-1] - self.lat_uv[:-1,-1] )
        lat[-2,-1] = self.lat_uv[-2:-1,-1] - ( self.lat_t[-2:-1,-1] - self.lat_uv[-2:-1,-1] )

        self.lon_uv_vert = lon
        self.lat_uv_vert = lat


    def _calculate_grid_angle(self):
        geod = pyproj.Geod(ellps='WGS84')
        az_forward, az_back, dx = geod.inv(self.lon_t[:,:-1], self.lat_t[:,:-1], \
                                           self.lon_t[:,1:], self.lat_t[:,1:])

        angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
        self.angle = (90 - angle) * np.pi/180.
        

    def _calculate_grid_angle(self):
        geod = pyproj.Geod(ellps='WGS84')
        az_forward, az_back, dx = geod.inv(self.lon_t_vert[:,:-1], self.lat_t_vert[:,:-1], \
                                           self.lon_t_vert[:,1:], self.lat_t_vert[:,1:])

        angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
        self.angle = (90 - angle) * np.pi/180.
