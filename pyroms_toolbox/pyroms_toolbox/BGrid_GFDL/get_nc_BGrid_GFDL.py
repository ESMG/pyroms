import numpy as np
import pyroms
from pyroms_toolbox.BGrid_GFDL import BGrid_GFDL


def get_nc_BGrid_GFDL(grdfile, name='GFDL_CM2.1_North_Pacific', \
#                             xrange=(80,189), yrange=(96,198)):
                             xrange=(60,175), yrange=(120, 190)):
    """
    Bgrd = get_nc_BGrid_GFDL(grdfile)

    Load B-Grid grid object for GFDL CM2.1 from netCDF grid file
    """

    nc = pyroms.io.Dataset(grdfile)

    lon_t = nc.variables['geolon_t'][:]
    lat_t = nc.variables['geolat_t'][:]
    lon_uv = nc.variables['geolon_c'][:]
    lat_uv = nc.variables['geolat_c'][:]

    h = nc.variables['ht'][:]

    f = nc.variables['coriolis_param'][:]

    kmt = nc.variables['kmt'][:]
    z_t = nc.variables['st_ocean'][:]
    z_t_edges = nc.variables['st_edges_ocean'][:]

    kmu = nc.variables['kmu'][:]
    z_uv = nc.variables['sw_ocean'][:]
    z_uv_edges = nc.variables['sw_edges_ocean'][:]

    # compute mask at t-point
    M_t, L_t = kmt.shape
    N_t = z_t.shape[0]
    mask_t = np.zeros((N_t, M_t, L_t))
    for j in range(M_t):
        for i in range(L_t):
            try:
                mask_t[0:kmt[j,i], j,i] = 1
            except:
                mask_t[:, j,i] = 0

    # compute mask at uv-point
    M_uv, L_uv = kmu.shape
    N_uv = z_uv.shape[0]
    mask_uv = np.zeros((N_uv, M_uv, L_uv))
    for j in range(M_uv):
        for i in range(L_uv):
            try:
                mask_uv[0:kmu[j,i], j,i] = 1
            except:
                mask_uv[:, j,i] = 0

    return BGrid_GFDL(lon_t, lat_t, lon_uv, lat_uv, \
                       mask_t, mask_uv, h, z_t, z_t_edges, \
                       z_uv, z_uv_edges, f, \
                       name, xrange=xrange, yrange=yrange)
