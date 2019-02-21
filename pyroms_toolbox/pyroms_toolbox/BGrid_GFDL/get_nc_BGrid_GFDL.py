import numpy as np
import pyroms
from pyroms_toolbox.BGrid_GFDL import BGrid_GFDL


def get_nc_BGrid_GFDL(grdfile, name='GFDL_CM2.1_North_Pacific', area='regional', \
                      xrange=(60,175), yrange=(120, 190), ystart=235):
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
                mask_t[0:int(kmt[j,i]), j,i] = 1
            except:
                mask_t[:, j,i] = 0

    # compute mask at uv-point
    M_uv, L_uv = kmu.shape
    N_uv = z_uv.shape[0]
    mask_uv = np.zeros((N_uv, M_uv, L_uv))
    for j in range(M_uv):
        for i in range(L_uv):
            try:
                mask_uv[0:int(kmu[j,i]), j,i] = 1
            except:
                mask_uv[:, j,i] = 0

    if area == 'npolar':
        #add two rows in the north and the south
        lon_t = lon_t[np.r_[0,0,:np.size(lon_t,0),-1,-1]]
        lon_t = lon_t[:,np.r_[0,:np.size(lon_t,1),-1]]
        lon_t[:,0] = lon_t[:,1] - (lon_t[:,2]-lon_t[:,1])
        lon_t[:,-1] = lon_t[:,-2] + (lon_t[:,-2]-lon_t[:,-3])
        lat_t = lat_t[np.r_[0,0,:np.size(lat_t,0),-1,-1]]
        lat_t = lat_t[:,np.r_[0,:np.size(lat_t,1),-1]]
        lat_t[0,:] = -85
        lat_t[1,:] = -80
        lat_t[-2,:] = 90
        lat_t[-1,:] = 91
        lon_uv = lon_uv[np.r_[0,0,:np.size(lon_uv,0),-1,-1]]
        lon_uv = lon_uv[:,np.r_[0,:np.size(lon_uv,1),-1]]
        lon_uv[:,0] = lon_uv[:,1] - (lon_uv[:,2]-lon_t[:,1])
        lon_uv[:,-1] = lon_uv[:,-2] + (lon_uv[:,-2]-lon_uv[:,-3])
        lat_uv = lat_uv[np.r_[0,0,:np.size(lat_uv,0),-1,-1]]
        lat_uv = lat_uv[:,np.r_[0,:np.size(lat_uv,1),-1]]
        lat_uv[0,:] = -85
        lat_uv[1,:] = -80
        lat_uv[-2,:] = 90
        lat_uv[-1,:] = 91
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
        f = f[np.r_[0,0,:np.size(f,0),-1,-1]]
        f = f[:,np.r_[0,:np.size(f,1),-1]]
        f[:,0] = f[:,-2]
        f[:,-1] = f[:,1]
        m,l = h.shape
        xrange=(1,l-2)
        yrange=(ystart+2,m-2)

    if area == 'tripole':
        #add two rows in the north and the south
        fold1 = L_t//2
        lon_t = lon_t[np.r_[0,0,:np.size(lon_t,0),-1,-1]]
        lon_t[-2,:fold1] = lon_t[-3,L_t:fold1-1:-1]
        lon_t[-2,L_t:fold1-1:-1] = lon_t[-3,:fold1]
        lon_t[-1,:fold1] = lon_t[-4,L_t:fold1-1:-1]
        lon_t[-1,L_t:fold1-1:-1] = lon_t[-4,:fold1]

        lon_t = lon_t[:,np.r_[0,:np.size(lon_t,1),-1]]
        lon_t[:,0] = lon_t[:,1] - (lon_t[:,2]-lon_t[:,1])
        lon_t[:,-1] = lon_t[:,-2] + (lon_t[:,-2]-lon_t[:,-3])
        lat_t = lat_t[np.r_[0,0,:np.size(lat_t,0),-1,-1]]
        lat_t = lat_t[:,np.r_[0,:np.size(lat_t,1),-1]]
        lat_t[0,:] = -85
        lat_t[1,:] = -80
        lat_t[-2,:] = lat_t[-3,:]
        lat_t[-1,:] = lat_t[-4,:]
        lon_uv = lon_uv[np.r_[0,0,:np.size(lon_uv,0),-1,-1]]

        lon_uv[-2,:fold1] = lon_uv[-4,L_t:fold1-1:-1]
        lon_uv[-2,L_t:fold1-1:-1] = lon_uv[-4,:fold1]
        lon_uv[-1,:fold1] = lon_uv[-5,L_t:fold1-1:-1]
        lon_uv[-1,L_t:fold1-1:-1] = lon_uv[-5,:fold1]

        lon_uv = lon_uv[:,np.r_[0,:np.size(lon_uv,1),-1]]
        lon_uv[:,0] = lon_uv[:,1] - (lon_uv[:,2]-lon_t[:,1])
        lon_uv[:,-1] = lon_uv[:,-2] + (lon_uv[:,-2]-lon_uv[:,-3])
        lat_uv = lat_uv[np.r_[0,0,:np.size(lat_uv,0),-1,-1]]
        lat_uv = lat_uv[:,np.r_[0,:np.size(lat_uv,1),-1]]
        lat_uv[0,:] = -85
        lat_uv[1,:] = -80
        lat_uv[-2,:] = lat_uv[-3,:]
        lat_uv[-1,:] = lat_uv[-4,:]
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
        f = f[np.r_[0,0,:np.size(f,0),-1,-1]]
        f = f[:,np.r_[0,:np.size(f,1),-1]]
        f[:,0] = f[:,-2]
        f[:,-1] = f[:,1]
        m,l = h.shape
        xrange=(1,l-2)
        yrange=(ystart+2,m-2)

    return BGrid_GFDL(lon_t, lat_t, lon_uv, lat_uv, \
                       mask_t, mask_uv, h, z_t, z_t_edges, \
                       z_uv, z_uv_edges, f, \
                       name, xrange=xrange, yrange=yrange)
