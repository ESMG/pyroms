import numpy as np
from mpl_toolkits.basemap import pyproj


def get_cell_area(lon, lat):

    if np.size(lon.shape) == 1:
        lon, lat  = np.meshgrid(lon, lat)


    #compute position of the vertices of the cell
    lonv = np.zeros((lon.shape[0]+1, lon.shape[1]+1))
    latv = np.zeros((lat.shape[0]+1, lat.shape[1]+1))

    lonv[1:-1,1:-1] = 0.25 * (lon[1:,1:] + lon[:-1,:-1] + lon[1:,:-1] + lon[:-1,1:]) 
    lonv[0,:] = lonv[1,:]
    lonv[-1,:] = lonv[-2,:]
    lonv[:,0] = lonv[:,1] - (lonv[:,2] - lonv[:,1])
    lonv[:,-1] = lonv[:,-2] + (lonv[:,-2] - lonv[:,-3])

    latv[1:-1,1:-1] = 0.25 * (lat[1:,1:] + lat[:-1,:-1] + lat[1:,:-1] + lat[:-1,1:])
    latv[:,0] = latv[:,1]
    latv[:,-1] = latv[:,-2]
    latv[0,:] = latv[1,:] - (latv[2,:] - latv[1,:])
    latv[-1,:] = latv[-2,:] + (latv[-2,:] - latv[-3,:])

    #get dx and dy based on great circle distances
    ellipse='WGS84'
    geod = pyproj.Geod(ellps=ellipse)
    az_forward, az_back, dx = geod.inv(lonv[:,1:], latv[:,1:], \
                                       lonv[:,:-1], latv[:,:-1])
    dx = 0.5*(dx[1:,:]+dx[:-1,:])
    az_forward, az_back, dy = geod.inv(lonv[1:,:], latv[1:,:], \
                                       lonv[:-1,:], latv[:-1,:])
    dy = 0.5*(dy[:,1:]+dy[:,:-1])

    area = dx * dy

    return dx, dy, area
