import numpy as np
import pyroms
from pyroms_toolbox.SODA import SODA_Grid


def get_nc_SODA_Grid(grdfile, name='SODA_2.0.2_CORAL', \
                         xrange=(185,340), yrange=(100, 210)):
    """
    grd = get_nc_SODA_Grid(grdfile)

    Load grid object for SODA 2.0.2 from netCDF file
    """

    nc = pyroms.io.Dataset(grdfile)

    lon = nc.variables['LON'][:]
    lat = nc.variables['LAT'][:]

    depth = nc.variables['DEPTH'][:]
    dep = nc.variables['DEPTH_bnds'][:]
    depth_bnds = np.zeros(depth.shape[0]+1)
    depth_bnds[:-1] = dep[:,0]
    depth_bnds[-1] = dep[-1,1]

    temp = nc.variables['TEMP']
    mask = np.array(~temp[0,:].mask, dtype='int')

    bottom = pyroms.utility.get_bottom(temp[0,::-1,:,:], mask[0,:], spval=temp.missing_value)
    nlev = mask.shape[0]
    bottom = (nlev-1) - bottom
    h = np.zeros(mask[0,:].shape)
    for i in range(mask[0,:].shape[1]):
        for j in range(mask[0,:].shape[0]):
            if mask[0,j,i] == 1:
                h[j,i] = depth_bnds[bottom[j,i]]

    return SODA_Grid(lon, lat, mask, depth, depth_bnds, h, \
                        name, xrange, yrange)
