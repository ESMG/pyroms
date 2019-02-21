import numpy as np
import pyroms
from .CGrid_TPXO8 import CGrid_TPXO8


def get_nc_CGrid_TPXO8(grdfile, name='TPXO8', \
                       xrange=(4350, 4550), yrange=(6600, 6900), missing_value=-9999):
    """
    grd = get_nc_CGrid_TPXO8(grdfile)

    Load Cgrid object for TPXO8 from netCDF file
    """

    nc = pyroms.io.Dataset(grdfile)

    lonh = nc.variables['lon_z'][:]
    lath = nc.variables['lat_z'][:]
    lonu = nc.variables['lon_u'][:]
    latu = nc.variables['lat_u'][:]
    lonv = nc.variables['lon_v'][:]
    latv = nc.variables['lat_v'][:]

    zh = nc.variables['hz'][:]
    zu = nc.variables['hu'][:]
    zv = nc.variables['hv'][:]
    nc.close()

    # land mask
    h_msk = zh!=0
    u_msk = zu!=0
    v_msk = zv!=0

    # longitude from -180 to 180
    lonh[lonh>180] = lonh[lonh>180]-360
    lonu[lonu>180] = lonu[lonu>180]-360
    lonv[lonv>180] = lonv[lonv>180]-360

    lathh, lonhh = np.meshgrid(lath, lonh)
    latuu, lonuu = np.meshgrid(latu, lonu)
    latvv, lonvv = np.meshgrid(latv, lonv)

    # generate tpxo8 grid
    # xrange = [4400, 4600]
    # yrange = [6600, 6900]
    return CGrid_TPXO8(lonhh, lathh, lonuu, latuu, lonvv, latvv, \
                                   h_msk, u_msk, v_msk, \
                                   zh, zu, zv, missing_value, 'TPXO8', xrange, yrange)

