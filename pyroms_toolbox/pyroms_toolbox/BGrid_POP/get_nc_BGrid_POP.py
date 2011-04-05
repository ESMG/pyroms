import numpy as np
import pyroms
from pyroms_toolbox.BGrid_POP import BGrid_POP


def get_nc_BGrid_POP(grdfile, name='POP_CHUKCHI'):
    """
    Bgrd = get_nc_BGrid_POP(grdfile)

    Load B-Grid grid object for POP from netCDF grid file
    """

    nc = pyroms.io.Dataset(grdfile)

    lon_t = nc.variables['clon'][:]
    lat_t = nc.variables['clat'][:]
    lon_uv = nc.variables['ulon'][:]
    lat_uv = nc.variables['ulat'][:]
    angle = nc.variables['angle'][:]
    kmt = nc.variables['kmt'][:]
    dz = nc.variables['dz'][:]

#    h = nc.variables['ht'][:]

#    f = nc.variables['coriolis_param'][:]

#    z_t = nc.variables['st_ocean'][:]
#    z_t_edges = nc.variables['st_edges_ocean'][:]

#    kmu = nc.variables['kmu'][:]
#    z_uv = nc.variables['sw_ocean'][:]
#    z_uv_edges = nc.variables['sw_edges_ocean'][:]

    # compute mask at t-point
    M_t, L_t = kmt.shape
    N_t = dz.shape[0]
    dz = dz * 0.01
    z_t = dz
    for k in range(1:N_t):
        z_t(k) = z_t(k-1) + dz(k)

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

    return BGrid_POP(lon_t, lat_t, lon_uv, lat_uv, \
                       mask_t, z_t, angle, name)
#    return BGrid_POP(lon_t, lat_t, lon_uv, lat_uv, \
#                       mask_t, mask_uv, h, z_t, z_t_edges, \
#                       z_uv, z_uv_edges, f, \
#                       name)
