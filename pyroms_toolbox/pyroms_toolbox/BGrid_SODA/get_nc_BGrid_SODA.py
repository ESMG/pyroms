import numpy as np
import pyroms
from pyroms_toolbox.BGrid_SODA import BGrid_SODA


def get_nc_BGrid_SODA(grdfile, name='SODA_2.1.6_CORAL', area='regional', \
                         xrange=(185,340), yrange=(100, 210), ystart=245):
    """
    grd = get_nc_BGrid_SODA(grdfile)

    Load Bgrid object for SODA 2.1.6 from netCDF file
    """

    nc = pyroms.io.Dataset(grdfile)

    lon_t = nc.variables['LON'][:]
    lat_t = nc.variables['LAT'][:]

    # All the data have been interpolated at the t-point
    # lon_t = lon_uv, lat_t = lat_uv
    #lon_uv = 0.5 * (lon_t[1:] + lon_t[:-1])
    #lat_uv = 0.5 * (lat_t[1:] + lat_t[:-1])
    lon_uv = lon_t
    lat_uv = lat_t

    depth = nc.variables['DEPTH'][:]
    dep = nc.variables['DEPTH_bnds'][:]
    depth_bnds = np.zeros(depth.shape[0]+1)
    depth_bnds[:-1] = dep[:,0]
    depth_bnds[-1] = dep[-1,1]

    nc_mask_t = nc.variables['MASK_T']
    mask_t = np.array(~nc_mask_t[:].mask, dtype='int')

    nc_mask_uv = nc.variables['MASK_UV']
    mask_uv = np.array(~nc_mask_uv[:].mask, dtype='int')

    bottom = pyroms.utility.get_bottom(nc_mask_t[::-1,:,:], mask_t[0,:], spval=nc_mask_t.missing_value)
    nlev = mask_t.shape[0]
    bottom = (nlev-1) - bottom
    h = np.zeros(mask_t[0,:].shape)
    for i in range(mask_t[0,:].shape[1]):
        for j in range(mask_t[0,:].shape[0]):
            if mask_t[0,j,i] == 1:
                h[j,i] = depth_bnds[int(bottom[j,i])]

    if area == 'global':
        #add one row in the north and the south
        lon_t = lon_t[np.r_[0,:len(lon_t),-1]]
        lon_t[0] = lon_t[1] - (lon_t[2]-lon_t[1])
        lon_t[-1] = lon_t[-2] + (lon_t[-2]-lon_t[-3])
        lat_t = lat_t[np.r_[0,0,:len(lat_t),-1,-1]]
        lat_t[0] = -85
        lat_t[1] = -80
        lat_t[-2] = 90
        lat_t[-1] = 91
        lon_uv = lon_t
        lat_uv = lat_t
        mask_t = mask_t[:,np.r_[0,0,:np.size(mask_t,1),-1,-1],:]
        mask_t = mask_t[:,:,np.r_[0,:np.size(mask_t,2),-1]]
        mask_t[:,:,0] = mask_t[:,:,-2]
        mask_t[:,:,-1] = mask_t[:,:,1]
        mask_uv = mask_uv[:,np.r_[0,0,:np.size(mask_uv,1),-1,-1],:]
        mask_uv = mask_uv[:,:,np.r_[0,:np.size(mask_uv,2),-1]]
        mask_uv[:,:,0] = mask_uv[:,:,-2]
        mask_uv[:,:,-1] = mask_uv[:,:,1]
        h = h[np.r_[0,0,:np.size(h,0),-1,-1]]
        h = h[:,np.r_[0,:np.size(h,1),-1]]
        h[:,0] = h[:,-2]
        h[:,-1] = h[:,1]
        m,l = h.shape
        xrange=(1,l-2)
        yrange=(1,m-2)

    if area == 'npolar':
        #add one row in the north and the south
        lon_t = lon_t[np.r_[0,:len(lon_t),-1]]
        lon_t[0] = lon_t[1] - (lon_t[2]-lon_t[1])
        lon_t[-1] = lon_t[-2] + (lon_t[-2]-lon_t[-3])
        lat_t = lat_t[np.r_[0,0,:len(lat_t),-1,-1]]
        lat_t[0] = -85
        lat_t[1] = -80
        lat_t[-2] = 90
        lat_t[-1] = 91
        lon_uv = lon_t
        lat_uv = lat_t
        mask_t = mask_t[:,np.r_[0,0,:np.size(mask_t,1),-1,-1],:]
        mask_t = mask_t[:,:,np.r_[0,:np.size(mask_t,2),-1]]
        mask_t[:,:,0] = mask_t[:,:,-2]
        mask_t[:,:,-1] = mask_t[:,:,1]
        mask_uv = mask_uv[:,np.r_[0,0,:np.size(mask_uv,1),-1,-1],:]
        mask_uv = mask_uv[:,:,np.r_[0,:np.size(mask_uv,2),-1]]
        mask_uv[:,:,0] = mask_uv[:,:,-2]
        mask_uv[:,:,-1] = mask_uv[:,:,1]
        h = h[np.r_[0,0,:np.size(h,0),-1,-1]]
        h = h[:,np.r_[0,:np.size(h,1),-1]]
        h[:,0] = h[:,-2]
        h[:,-1] = h[:,1]
        m,l = h.shape
        xrange=(1,l-2)
        yrange=(ystart+2,m-2)

    return BGrid_SODA(lon_t, lat_t, lon_uv, lat_uv, mask_t, mask_uv, depth, depth_bnds, h, \
                        name, xrange, yrange)

