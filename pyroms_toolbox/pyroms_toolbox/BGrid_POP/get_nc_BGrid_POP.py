import numpy as np
import pyroms
from pyroms_toolbox.BGrid_POP import BGrid_POP


def get_nc_BGrid_POP(grdfile, name='POP_NEP', \
                         xrange=(170,270), yrange=(240, 350)):
    """
    grd = get_nc_BGrid_POP(grdfile)

    Load Bgrid object for POP from netCDF file
    """

    nc = pyroms.io.Dataset(grdfile)

    lon_t = nc.variables['TLONG'][:]
    lat_t = nc.variables['TLAT'][:]

    lon_u = nc.variables['ULONG'][:]
    lat_u = nc.variables['ULAT'][:]

    angle = nc.variables['ANGLET'][:]

    h_t = nc.variables['HT'][:]
    h_u = nc.variables['HU'][:]

    z_t = nc.variables['z_t'][:]
    
    z_w_top = nc.variables['z_w_top'][:]
    z_w_bot = nc.variables['z_w_bot'][:]
    z_w = np.zeros(z_t.size + 1)
    z_w[:-1] = z_w_top
    z_w[-1] = z_w_bot[-1]


    return BGrid_POP(lon_t, lat_t, lon_u, lat_u, angle, h_t, h_u, z_t, z_w, \
                        name, xrange, yrange)
