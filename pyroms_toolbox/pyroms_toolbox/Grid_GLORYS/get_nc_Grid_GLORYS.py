import numpy as np
import pyroms
from pyroms_toolbox.Grid_GLORYS import Grid_GLORYS
import pdb


def get_nc_Grid_GLORYS(grdfile, name='GLORYS_NWGOA', area='regional', \
                         irange=(270,460), jrange=(150, 328), ystart=245):
    """
    grd = get_nc_Grid_GLORYS(grdfile)

    Load A-grid object for GLORYS from netCDF file
    """

    nc = pyroms.io.Dataset(grdfile)

    lon_t = nc.variables['longitude'][:]
    lat_t = nc.variables['latitude'][:]
#   lat_u = nc.variables['gphiu'][:]
#   lon_u = nc.variables['glamu'][:]
#   lat_v = nc.variables['gphiv'][:]
#   lon_v = nc.variables['glamv'][:]

    depth = nc.variables['depth'][:]
#   depth_w = nc.variables['gdepw_0'][:]
    depth_bnds = np.zeros(depth.shape[0]+1)
#   depth_bnds[:-1] = depth_w[:]
    depth_bnds[-1] = 6000.

    nc_mask_t = nc.variables['zos']
#    mask_t = np.array(~nc_mask_t[:].mask, dtype='int')
    mask_t = np.where(np.array(nc_mask_t[:], dtype='int')==nc_mask_t._FillValue, 0, 1)
#   pdb.set_trace()
    nc_mask_t = nc.variables['thetao']

    bottom = pyroms.utility.get_bottom(nc_mask_t[0,::-1,:,:], mask_t[0,:], spval=nc_mask_t._FillValue)
    nlev = mask_t.shape[0]
    bottom = (nlev-1) - bottom
    h = np.zeros(mask_t[0,:].shape)
    for i in range(mask_t[0,:].shape[1]):
        for j in range(mask_t[0,:].shape[0]):
            if mask_t[0,j,i] == 1:
                h[j,i] = depth_bnds[int(bottom[j,i])]

    if area == 'global':
        #add rows in the north and the south, east and west
        lon_t = lon_t[:,np.r_[0,:np.size(lon_t,1),-1]]
        lon_t[:,0] = lon_t[:,1] - (lon_t[:,2]-lon_t[:,1])
        lon_t[:,-1] = lon_t[:,-2] + (lon_t[:,-2]-lon_t[:,-3])
        lon_t = lon_t[np.r_[0,0,:np.size(lon_t,0),-1,-1]]

        lat_t = lat_t[np.r_[0,0,:np.size(lat_t,0),-1,-1]]
        lat_t[-1,:] = -80
        lat_t[0,:] = -85
        lat_t[-2,:] = lat_t[-3,:]
        lat_t[-1,:] = lat_t[-4,:]
        lat_t = lat_t[:,np.r_[0,:np.size(lat_t,1),-1]]

        mask_t = mask_t[:,np.r_[0,0,:np.size(mask_t,1),-1,-1],:]
        mask_t = mask_t[:,:,np.r_[0,:np.size(mask_t,2),-1]]
        mask_t[:,:,0] = mask_t[:,:,-2]
        mask_t[:,:,-1] = mask_t[:,:,1]
        h = h[np.r_[0,0,:np.size(h,0),-1,-1]]
        h = h[:,np.r_[0,:np.size(h,1),-1]]
        h[:,0] = h[:,-2]
        h[:,-1] = h[:,1]
        m,l = h.shape
        irange=(1,l-2)
        jrange=(1,m-2)

    if area == 'npolar':
        #add rows in the north and the south, east and west
        lon_t = lon_t[:,np.r_[0,:np.size(lon_t,1),-1]]
        lon_t[:,0] = lon_t[:,1] - (lon_t[:,2]-lon_t[:,1])
        lon_t[:,-1] = lon_t[:,-2] + (lon_t[:,-2]-lon_t[:,-3])
        lon_t = lon_t[np.r_[0,0,:np.size(lon_t,0),-1,-1]]

        lat_t = lat_t[np.r_[0,0,:np.size(lat_t,0),-1,-1]]
        lat_t[-1,:] = -80
        lat_t[0,:] = -85
        lat_t[-2,:] = lat_t[-3,:]
        lat_t[-1,:] = lat_t[-4,:]
        lat_t = lat_t[:,np.r_[0,:np.size(lat_t,1),-1]]

        mask_t = mask_t[:,np.r_[0,0,:np.size(mask_t,1),-1,-1],:]
        mask_t = mask_t[:,:,np.r_[0,:np.size(mask_t,2),-1]]
        mask_t[:,:,0] = mask_t[:,:,-2]
        mask_t[:,:,-1] = mask_t[:,:,1]
        h = h[np.r_[0,0,:np.size(h,0),-1,-1]]
        h = h[:,np.r_[0,:np.size(h,1),-1]]
        h[:,0] = h[:,-2]
        h[:,-1] = h[:,1]
        m,l = h.shape
        irange=(1,l-2)
        jrange=(ystart+2,m-2)

    return Grid_GLORYS(lon_t, lat_t, mask_t, depth, depth_bnds, h, \
                        name, irange, jrange)

