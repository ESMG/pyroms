import numpy as np
import pyroms
from pyroms_toolbox.CGrid_GLORYS import CGrid_GLORYS


def get_nc_CGrid_GLORYS(grdfile, name='GLORYS_CORAL', area='regional', \
                         xrange=(185,340), yrange=(100, 210), ystart=245):
    """
    grd = get_nc_CGrid_GLORYS(grdfile)

    Load Cgrid object for GLORYS from netCDF file
    """

    nc = pyroms.io.Dataset(grdfile)

    lon_t = nc.variables['nav_lon'][:]
    lat_t = nc.variables['nav_lat'][:]
    lat_u = nc.variables['gphiu'][:]
    lon_u = nc.variables['glamu'][:]
    lat_v = nc.variables['gphiv'][:]
    lon_v = nc.variables['glamv'][:]

    depth = nc.variables['gdept_0'][:]
    depth_w = nc.variables['gdepw_0'][:]
    depth_bnds = np.zeros(depth.shape[0]+1)
    depth_bnds[:-1] = depth_w[:]
    depth_bnds[-1] = 6000.

    nc_mask_t = nc.variables['tmask']
#    mask_t = np.array(~nc_mask_t[:].mask, dtype='int')
    mask_t = np.array(nc_mask_t[:], dtype='int')

    nc_mask_u = nc.variables['umask']
#    mask_u = np.array(~nc_mask_u[:].mask, dtype='int')
    mask_u = np.array(nc_mask_u[:], dtype='int')

    nc_mask_v = nc.variables['vmask']
#    mask_v = np.array(~nc_mask_v[:].mask, dtype='int')
    mask_v = np.array(nc_mask_v[:], dtype='int')

    bottom = pyroms.utility.get_bottom(nc_mask_t[::-1,:,:], mask_t[0,:], spval=nc_mask_t.missing_value)
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

        lon_u = lon_u[:,np.r_[0,:np.size(lon_u,1),-1]]
        lon_u[:,0] = lon_u[:,1] - (lon_u[:,2]-lon_u[:,1])
        lon_u[:,-1] = lon_u[:,-2] + (lon_u[:,-2]-lon_u[:,-3])
        lon_u = lon_u[np.r_[0,0,:np.size(lon_u,0),-1,-1]]

        lat_u = lat_u[np.r_[0,0,:np.size(lat_u,0),-1,-1]]
        lat_u[-1,:] = -80
        lat_u[0,:] = -85
        lat_u[-2,:] = lat_u[-3,:]
        lat_u[-1,:] = lat_u[-4,:]
        lat_u = lat_u[:,np.r_[0,:np.size(lat_u,1),-1]]

        lon_v = lon_v[:,np.r_[0,:np.size(lon_v,1),-1]]
        lon_v[:,0] = lon_v[:,1] - (lon_v[:,2]-lon_v[:,1])
        lon_v[:,-1] = lon_v[:,-2] + (lon_v[:,-2]-lon_v[:,-3])
        lon_v = lon_v[np.r_[0,0,:np.size(lon_v,0),-1,-1]]

        lat_v = lat_v[np.r_[0,0,:np.size(lat_v,0),-1,-1]]
        lat_v[-1,:] = -80
        lat_v[0,:] = -85
        lat_v[-2,:] = lat_v[-3,:]
        lat_v[-1,:] = lat_v[-4,:]
        lat_v = lat_v[:,np.r_[0,:np.size(lat_v,1),-1]]

        mask_t = mask_t[:,np.r_[0,0,:np.size(mask_t,1),-1,-1],:]
        mask_t = mask_t[:,:,np.r_[0,:np.size(mask_t,2),-1]]
        mask_t[:,:,0] = mask_t[:,:,-2]
        mask_t[:,:,-1] = mask_t[:,:,1]
        mask_u = mask_u[:,np.r_[0,0,:np.size(mask_u,1),-1,-1],:]
        mask_u = mask_u[:,:,np.r_[0,:np.size(mask_u,2),-1]]
        mask_u[:,:,0] = mask_u[:,:,-2]
        mask_u[:,:,-1] = mask_u[:,:,1]
        mask_v = mask_v[:,np.r_[0,0,:np.size(mask_v,1),-1,-1],:]
        mask_v = mask_v[:,:,np.r_[0,:np.size(mask_v,2),-1]]
        mask_v[:,:,0] = mask_v[:,:,-2]
        mask_v[:,:,-1] = mask_v[:,:,1]
        h = h[np.r_[0,0,:np.size(h,0),-1,-1]]
        h = h[:,np.r_[0,:np.size(h,1),-1]]
        h[:,0] = h[:,-2]
        h[:,-1] = h[:,1]
        m,l = h.shape
        xrange=(1,l-2)
        yrange=(1,m-2)

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

        lon_u = lon_u[:,np.r_[0,:np.size(lon_u,1),-1]]
        lon_u[:,0] = lon_u[:,1] - (lon_u[:,2]-lon_u[:,1])
        lon_u[:,-1] = lon_u[:,-2] + (lon_u[:,-2]-lon_u[:,-3])
        lon_u = lon_u[np.r_[0,0,:np.size(lon_u,0),-1,-1]]

        lat_u = lat_u[np.r_[0,0,:np.size(lat_u,0),-1,-1]]
        lat_u[-1,:] = -80
        lat_u[0,:] = -85
        lat_u[-2,:] = lat_u[-3,:]
        lat_u[-1,:] = lat_u[-4,:]
        lat_u = lat_u[:,np.r_[0,:np.size(lat_u,1),-1]]

        lon_v = lon_v[:,np.r_[0,:np.size(lon_v,1),-1]]
        lon_v[:,0] = lon_v[:,1] - (lon_v[:,2]-lon_v[:,1])
        lon_v[:,-1] = lon_v[:,-2] + (lon_v[:,-2]-lon_v[:,-3])
        lon_v = lon_v[np.r_[0,0,:np.size(lon_v,0),-1,-1]]

        lat_v = lat_v[np.r_[0,0,:np.size(lat_v,0),-1,-1]]
        lat_v[-1,:] = -80
        lat_v[0,:] = -85
        lat_v[-2,:] = lat_v[-3,:]
        lat_v[-1,:] = lat_v[-4,:]
        lat_v = lat_v[:,np.r_[0,:np.size(lat_v,1),-1]]

        mask_t = mask_t[:,np.r_[0,0,:np.size(mask_t,1),-1,-1],:]
        mask_t = mask_t[:,:,np.r_[0,:np.size(mask_t,2),-1]]
        mask_t[:,:,0] = mask_t[:,:,-2]
        mask_t[:,:,-1] = mask_t[:,:,1]
        mask_u = mask_u[:,np.r_[0,0,:np.size(mask_u,1),-1,-1],:]
        mask_u = mask_u[:,:,np.r_[0,:np.size(mask_u,2),-1]]
        mask_u[:,:,0] = mask_u[:,:,-2]
        mask_u[:,:,-1] = mask_u[:,:,1]
        mask_v = mask_v[:,np.r_[0,0,:np.size(mask_v,1),-1,-1],:]
        mask_v = mask_v[:,:,np.r_[0,:np.size(mask_v,2),-1]]
        mask_v[:,:,0] = mask_v[:,:,-2]
        mask_v[:,:,-1] = mask_v[:,:,1]
        h = h[np.r_[0,0,:np.size(h,0),-1,-1]]
        h = h[:,np.r_[0,:np.size(h,1),-1]]
        h[:,0] = h[:,-2]
        h[:,-1] = h[:,1]
        m,l = h.shape
        xrange=(1,l-2)
        yrange=(ystart+2,m-2)

    return CGrid_GLORYS(lon_t, lat_t, lon_u, lat_u, lon_v, lat_v, mask_t, mask_u, mask_v, depth, depth_bnds, h, \
                        name, xrange, yrange)

