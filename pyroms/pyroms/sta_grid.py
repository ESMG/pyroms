# encoding: utf-8

import sys
import os
import numpy as np
from mpl_toolkits.basemap import Basemap
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF

import pyroms
from pyroms.sta_hgrid import *
from pyroms.vgrid import *
from pyroms.grid import *
from pyroms import io

class Stations_Grid(object):
    """
    grd = Stations_Grid(sta_hgrid, vgrid)

    Stations Grid object combining horizontal and vertical grid
    """

    def __init__(self, name, hgrid=Sta_CGrid, vgrid=s_coordinate):
        self.name = name
        self.hgrid = hgrid
        self.vgrid = vgrid


def get_Stations_hgrid(gridid, sta_file):
    """
    hgrid = get_Stations_hgrid(gridid, sta_file)

    Load Stations horizontal grid object
    """

    gridinfo = ROMS_gridinfo(gridid)
    grdfile = gridinfo.grdfile

    nc = io.Dataset(sta_file)

    #Check for cartesian or geographical grid
    spherical = nc.variables['spherical'][:]

    #Get horizontal grid 
    if ((spherical == 0) or (spherical == 'F')):
        #cartesian grid
        print('Load cartesian grid from file')
        if 'x_rho' in list(nc.variables.keys()) and 'y_rho' in list(nc.variables.keys()):
            x_rho = nc.variables['x_rho'][:]
            y_rho = nc.variables['y_rho'][:]
            try: angle = nc.variables['angle'][:]
            except: angle = np.zeros(x_rho.shape)
        else:
            raise ValueError('NetCDF file must contain x_rho and y_rho \
                     and possibly angle for a cartesian grid')

        x_rho = nc.variables['x_rho'][:]
        y_rho = nc.variables['y_rho'][:]

        if 'angle' in list(nc.variables.keys()):
            angle = nc.variables['angle'][:]
        else:
            angle = None

        #Get cartesian grid
        hgrd = Sta_CGrid(x_rho, y_rho, angle_rho=angle)

    else:
        #geographical grid
        print('Load geographical grid from file')
        proj = Basemap(projection='merc', resolution=None, lat_0=0, lon_0=0)
        if 'lon_rho' in list(nc.variables.keys()) and 'lat_rho' in list(nc.variables.keys()):
            lon_rho = nc.variables['lon_rho'][:]
            lat_rho = nc.variables['lat_rho'][:]
        else:
            raise ValueError('NetCDF file must contain lon_rho and lat_rho \
                  for a geographical grid')

        lon_rho = nc.variables['lon_rho'][:]
        lat_rho = nc.variables['lat_rho'][:]

        if 'angle' in list(nc.variables.keys()):
            angle = nc.variables['angle'][:]
        else:
            angle = None

        #Get geographical grid   
        hgrd = Sta_CGrid_geo(lon_rho, lat_rho, proj, \
                         angle_rho=angle)

    return hgrd



def get_Stations_grid(gridid, sta_file, zeta=None):
    """
    grd = get_Stations_grid(gridid, sta_file, zeta=None)

    Load Stations grid object.

    gridid is a string with the name of the grid in it.  sta_file
    is the name of a stations file to read.

    grd.vgrid is a s_coordinate or
    a z_coordinate object, depending on gridid.grdtype.
    grd.vgrid.z_r and grd.vgrid.z_w (grd.vgrid.z for a 
    z_coordinate object) can be indexed in order to retreive the 
    actual depths. The free surface time series zeta can be provided 
    as an optional argument. Note that the values of zeta are not 
    calculated until z is indexed, so a netCDF variable for zeta may 
    be passed, even if the file is large, as only the values that 
    are required will be retrieved from the file.
    """

    gridinfo = ROMS_gridinfo(gridid, grid_file=sta_file, hist_file=sta_file)
    name = gridinfo.name

    hgrd = get_Stations_hgrid(gridid, sta_file)
    vgrid = get_ROMS_vgrid(gridid, zeta=zeta)

    #Get Stations grid
    return Stations_Grid(name, hgrd, vgrid)


def write_Stations_grid(grd, filename='roms_grd.nc'):
    """
    write_Stations_grid(grd, filename)

    Write Stations_CGrid class on a NetCDF file.
    """

    Sm = grd.hgrid.x_rho.shape[0]

    
    # Write Stations grid to file
    nc = netCDF.Dataset(filename, 'w', format='NETCDF3_64BIT')
    nc.Description = 'Stations grid'
    nc.Author = 'pyroms.grid.write_grd'
    nc.Created = datetime.now().isoformat()
    nc.type = 'Stations grid file'

    nc.createDimension('station', Sm)

    if hasattr(grd.vgrid, 's_rho') is True and grd.vgrid.s_rho is not None:
        N, = grd.vgrid.s_rho.shape
        nc.createDimension('s_rho', N)
        nc.createDimension('s_w', N+1)

    def write_nc_var(var, name, dimensions, long_name=None, units=None):
        nc.createVariable(name, 'f8', dimensions)
        if long_name is not None:
            nc.variables[name].long_name = long_name
        if units is not None:
            nc.variables[name].units = units
        nc.variables[name][:] = var
        print(' ... wrote ', name)

    if hasattr(grd.vgrid, 's_rho') is True and grd.vgrid.s_rho is not None:
        write_nc_var(grd.vgrid.theta_s, 'theta_s', (), 'S-coordinate surface control parameter')
        write_nc_var(grd.vgrid.theta_b, 'theta_b', (), 'S-coordinate bottom control parameter')
        write_nc_var(grd.vgrid.Tcline, 'Tcline', (), 'S-coordinate surface/bottom layer width', 'meter')
        write_nc_var(grd.vgrid.hc, 'hc', (), 'S-coordinate parameter, critical depth', 'meter')
        write_nc_var(grd.vgrid.s_rho, 's_rho', ('s_rho'), 'S-coordinate at RHO-points')
        write_nc_var(grd.vgrid.s_w, 's_w', ('s_w'), 'S-coordinate at W-points')
        write_nc_var(grd.vgrid.Cs_r, 'Cs_r', ('s_rho'), 'S-coordinate stretching curves at RHO-points')
        write_nc_var(grd.vgrid.Cs_w, 'Cs_w', ('s_w'), 'S-coordinate stretching curves at W-points')

    write_nc_var(grd.vgrid.h, 'h', ('station'), 'bathymetry at RHO-points', 'meter')
    write_nc_var(grd.hgrid.x_rho, 'x_rho', ('station'), 'x location of RHO-points', 'meter')
    write_nc_var(grd.hgrid.y_rho, 'y_rho', ('station'), 'y location of RHO-points', 'meter')

    if hasattr(grd.hgrid, 'lon_rho'):
        write_nc_var(grd.hgrid.lon_rho, 'lon_rho', ('station'), 'longitude of RHO-points', 'degree_east')
        write_nc_var(grd.hgrid.lat_rho, 'lat_rho', ('station'), 'latitude of RHO-points', 'degree_north')

    nc.createVariable('spherical', 'c')
    nc.variables['spherical'].long_name = 'Grid type logical switch'
    nc.variables['spherical'][:] = grd.hgrid.spherical
    print(' ... wrote ', 'spherical')

    write_nc_var(grd.hgrid.angle_rho, 'angle', ('station'), 'angle between XI-axis and EAST', 'radians')

    nc.close()
