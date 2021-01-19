# encoding: utf-8

import numpy as np
from . import _iso
from . import _obs_interp

import pyroms


def zslice(var, depth, grd, Cpos='rho', vert=False, mode='linear'):
    """
    zslice, lon, lat = zslice(var, depth, grd)

    optional switch:
      - Cpos='rho', 'u', 'v' or 'w'  specify the C-grid position where
                                     the variable rely
      - vert=True/False              If True, return the position of
                                     the verticies
      - mode='linear' or 'spline'    specify the type of interpolation

    return a constant-z slice at depth depth from 3D variable var
    lon and lat contain the C-grid position of the slice for plotting.
    If vert=True, lon and lat contain contain the position of the
    verticies (to be used with pcolor)
    """

    if mode=='linear':
        imode=0
    elif mode=='spline':
        imode=1
    else:
        imode=0
        raise Warning('%s not supported, defaulting to linear' % mode)


    # compute the depth on Arakawa-C grid position

    if Cpos == 'u':
        # average z_r at Arakawa-C u points
        z = 0.5 * (grd.vgrid.z_r[0,:,:,:-1] + grd.vgrid.z_r[0,:,:,1:])
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
                y = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:,:-1] + grd.hgrid.x_vert[:,1:])
                y = 0.5 * (grd.hgrid.y_vert[:,:-1] + grd.hgrid.y_vert[:,1:])
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_u[:]
                y = grd.hgrid.lat_u[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_u[:]
                y = grd.hgrid.y_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos == 'v':
        # average z_r at Arakawa-C v points
        z = 0.5 * (grd.vgrid.z_r[0,:,:-1,:] + grd.vgrid.z_r[0,:,1:,:])
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
                y = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:-1,:] + grd.hgrid.x_vert[1:,:])
                y = 0.5 * (grd.hgrid.y_vert[:-1,:] + grd.hgrid.y_vert[1:,:])
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_v[:]
                y = grd.hgrid.lat_v[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_v[:]
                y = grd.hgrid.y_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos == 'w':
        z = grd.vgrid.z_w[0,:]
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[:]
                y = grd.hgrid.lat_vert[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_vert[:]
                y = grd.hgrid.y_vert[:]
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    elif Cpos == 'rho':
        # for temp, salt, rho
        z = grd.vgrid.z_r[0,:]
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[:]
                y = grd.hgrid.lat_vert[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_vert[:]
                y = grd.hgrid.y_vert[:]
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)

    assert len(z.shape) == 3, 'z must be 3D'
    assert len(var.shape) == 3, 'var must be 3D'
    assert z.shape == var.shape, 'data and prop must be the same size'

    depth = -abs(depth)
    depth = depth * np.ones(z.shape[1:])

    zslice = _iso.zslice(z, var, depth, imode)

    # mask land
    zslice = np.ma.masked_where(mask == 0, zslice)
    # mask region with shalower depth than requisted depth
    zslice = np.ma.masked_where(zslice == 1e20, zslice)

    return zslice, x, y


def sslice(var, sindex, grd, Cpos='rho', vert=False):
    """
    sslice, lon, lat = sslice(var, sindex, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'       specify the C-grid position where
                                     the variable rely
      - vert=True/False              If True, return the position of
                                     the verticies
      - mode='linear' or 'spline'    specify the type of interpolation

    return a constant-s slice at index sindex from 3D variable var
    lon and lat contain the C-grid position of the slice for plotting.
    If vert=True, lon and lat contain contain the position of the
    verticies (to be used with pcolor)
    """

    # compute the depth on Arakawa-C grid position

    if Cpos == 'u':
        # average z_r at Arakawa-C u points
        z = 0.5 * (grd.vgrid.z_r[0,:,:,:-1] + grd.vgrid.z_r[0,:,:,1:])
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
                y = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:,:-1] + grd.hgrid.x_vert[:,1:])
                y = 0.5 * (grd.hgrid.y_vert[:,:-1] + grd.hgrid.y_vert[:,1:])
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_u[:]
                y = grd.hgrid.lat_u[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_u[:]
                y = grd.hgrid.y_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos == 'v':
        # average z_r at Arakawa-C v points
        z = 0.5 * (grd.vgrid.z_r[0,:,:-1,:] + grd.vgrid.z_r[0,:,1:,:])
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
                y = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:-1,:] + grd.hgrid.x_vert[1:,:])
                y = 0.5 * (grd.hgrid.y_vert[:-1,:] + grd.hgrid.y_vert[1:,:])
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_v[:]
                y = grd.hgrid.lat_v[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_v[:]
                y = grd.hgrid.y_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos == 'w':
        z = grd.vgrid.z_w[0,:]
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[:]
                y = grd.hgrid.lat_vert[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_vert[:]
                y = grd.hgrid.y_vert[:]
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    elif Cpos == 'rho':
        # for temp, salt, rho
        z = grd.vgrid.z_r[0,:]
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[:]
                y = grd.hgrid.lat_vert[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_vert[:]
                y = grd.hgrid.y_vert[:]
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)

    assert len(var.shape) == 3, 'var must be 3D'

    sslice = var[sindex,:,:]

    # mask land
    sslice = np.ma.masked_where(mask == 0, sslice)

    return sslice, x, y



def islice(var, iindex, grd, Cpos='rho', vert=False):
    """
    islice, z, lon, lat = islice(var, iindex, grd)

    optional switch:
      - Cpos='rho', 'u', 'v' or 'w'  specify the C-grid position where
                                     the variable rely
      - vert=True/False              If True, return the position of
                                     the verticies


    return a constant-i slice at index iindex from 3D variable var
    lon, lat and z contain the C-grid position of the slice for plotting.
    If vert=True, lon, lat and z contain contain the position of the
    verticies (to be used with pcolor)
    """

    # compute the depth on Arakawa-C grid position

    if Cpos == 'u':
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = np.concatenate((z[:,0:1,:], z, z[:,-1:,:]), 1)
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[:,1:-1]
                y = grd.hgrid.lat_vert[:,1:-1]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_vert[:,1:-1]
                y = grd.hgrid.y_vert[:,1:-1]
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_u[:]
                y = grd.hgrid.lat_u[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_u[:]
                y = grd.hgrid.y_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos == 'v':
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho
                y = grd.hgrid.lat_rho
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho
                y = grd.hgrid.y_rho
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_v[:]
                y = grd.hgrid.lat_v[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_v[:]
                y = grd.hgrid.y_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos == 'w':
        # for w, AKt, ...
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:-1,:,:] + z[1:,:,:])
            z = np.concatenate((np.array(grd.vgrid.z_w[0,0,:,:], ndmin=3), \
                                z, \
                                np.array(grd.vgrid.z_w[0,-1,:,:], ndmin=3)), 0)
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-1:]), 2)
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
                y = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:-1,:] + grd.hgrid.x_vert[1:,:])
                y = 0.5 * (grd.hgrid.y_vert[:-1,:] + grd.hgrid.y_vert[1:,:])
        else:
            z = grd.vgrid.z_r[0,:]
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]


    elif Cpos == 'rho':
        # for temp, salt, rho, ...
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = np.concatenate((z[:,0:1,:], z, z[:,-1:,:]), 1)
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
                y = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:,:-1] + grd.hgrid.y_vert[:,1:])
                y = 0.5 * (grd.hgrid.y_vert[:,:-1] + grd.hgrid.y_vert[:,1:])
        else:
            z = grd.vgrid.z_r[0,:]
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)

    # get constant-i slice
    vari = var[:,:,iindex]
    zi = z[:,:,iindex]
    xi = np.tile(x[:,iindex], (zi.shape[0], 1))
    yi = np.tile(y[:,iindex], (zi.shape[0], 1))

    # land/sea mask
    maski = np.tile(mask[:,iindex], (vari.shape[0], 1))
    vari = np.ma.masked_where(maski[:,:] == 0, vari[:,:])

    return vari, zi, xi, yi


def jslice(var, jindex, grd, Cpos='rho', vert=False):
    """
    jslice, z, lon, lat = jslice(var, jindex, grd)

    optional switch:
      - Cpos='rho', 'u', 'v' or 'w'  specify the C-grid position where
                                     the variable rely
      - vert=True/False              If True, return the position of
                                     the verticies


    return a constant-j slice at index jindex from 3D variable var
    lon, lat and z contain the C-grid position of the slice for plotting.
    If vert=True, lon, lat and z contain contain the position of the
    verticies (to be used with pcolor)
    """

    # compute the depth on Arakawa-C grid position

    if Cpos == 'u':
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho
                y = grd.hgrid.lat_rho
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho
                y = grd.hgrid.y_rho
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_u[:]
                y = grd.hgrid.lat_u[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_u[:]
                y = grd.hgrid.y_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos == 'v':
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-1:]), 2)
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[1:-1,:]
                y = grd.hgrid.lat_vert[1:-1,:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_vert[1:-1,:]
                y = grd.hgrid.y_vert[1:-1,:]
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_v[:]
                y = grd.hgrid.lat_v[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_v[:]
                y = grd.hgrid.y_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos == 'w':
        # for w, AKt, ...
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:-1,:,:] + z[1:,:,:])
            z = np.concatenate((np.array(grd.vgrid.z_w[0,0,:,:], ndmin=3), \
                                z, \
                                np.array(grd.vgrid.z_w[0,-1,:,:], ndmin=3)), 0)
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-1:]), 2)
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
                y = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:-1,:] + grd.hgrid.x_vert[1:,:])
                y = 0.5 * (grd.hgrid.y_vert[:-1,:] + grd.hgrid.y_vert[1:,:])
        else:
            z = grd.vgrid.z_r[0,:]
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    elif Cpos == 'rho':
        # for temp, salt, rho, ...
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-1:]), 2)
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
                y = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:-1,:] + grd.hgrid.x_vert[1:,:])
                y = 0.5 * (grd.hgrid.y_vert[:-1,:] + grd.hgrid.y_vert[1:,:])
        else:
            z = grd.vgrid.z_r[0,:]
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)

    # get constant-j slice
    varj = var[:,jindex,:]
    zj = z[:,jindex,:]
    xj = np.tile(x[jindex,:], (zj.shape[0], 1))
    yj = np.tile(y[jindex,:], (zj.shape[0], 1))

    # land/sea mask
    maskj = np.tile(mask[jindex,:], (varj.shape[0], 1))
    varj = np.ma.masked_where(maskj[:,:] == 0, varj[:,:])

    return varj, zj, xj, yj



def isoslice(var,prop,isoval, grd, Cpos='rho', masking=True, vert=False):
    """
    isoslice, lon, lat = isoslice(variable,property, isoval, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'       specify the C-grid position where
                                     the variable rely
      - masking=True                 mask the output if True
      - vert=True/False              If True, return the position of
                                     the verticies
      - mode='linear' or 'spline'    specify the type of interpolation


    result is a projection of variable at property == isoval in the first
    nonsingleton dimension.  In the case when there is more than one zero
    crossing, the results are averaged.
    lon, and lat contain the C-grid position of the slice for plotting.
    If vert=True, lon and lat and z contain contain the position of the
    verticies (to be used with pcolor)

    EXAMPLE:
    s_at_m5  = isoslice(s,z,-5);        # s at z == -5
    h_at_s30 = isoslice(z,s,30);        # z at s == 30
    """
    if (len(var.squeeze().shape)<=2):
        raise ValueError('variable must have at least two dimensions')
    if not prop.shape == var.shape:
        raise ValueError('dimension of var and prop must be identical')

    # compute the depth on Arakawa-C grid position

    if Cpos == 'u':
        # average z_r at Arakawa-C u points
        z = 0.5 * (grd.vgrid.z_r[0,:,:,:-1] + grd.vgrid.z_r[0,:,:,1:])
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
                y = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:,:-1] + grd.hgrid.x_vert[:,1:])
                y = 0.5 * (grd.hgrid.y_vert[:,:-1] + grd.hgrid.y_vert[:,1:])
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_u[:]
                y = grd.hgrid.lat_u[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_u[:]
                y = grd.hgrid.y_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos == 'v':
        # average z_r at Arakawa-C v points
        z = 0.5 * (grd.vgrid.z_r[0,:,:-1,:] + grd.vgrid.z_r[0,:,1:,:])
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
                y = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:-1,:] + grd.hgrid.x_vert[1:,:])
                y = 0.5 * (grd.hgrid.y_vert[:-1,:] + grd.hgrid.y_vert[1:,:])
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_v[:]
                y = grd.hgrid.lat_v[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_v[:]
                y = grd.hgrid.y_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos == 'w':
        z = grd.vgrid.z_w[0,:]
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[:]
                y = grd.hgrid.lat_vert[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_vert[:]
                y = grd.hgrid.y_vert[:]
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    elif Cpos == 'rho':
        # for temp, salt, rho
        z = grd.vgrid.z_r[0,:]
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[:]
                y = grd.hgrid.lat_vert[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_vert[:]
                y = grd.hgrid.y_vert[:]
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)

    prop = prop-isoval
    sz = np.shape(var)
    var = var.reshape(sz[0],-1)
    prop = prop.reshape(sz[0],-1)
    #find zero-crossings (zc == 1)
    zc =  np.where( (prop[:-1,:] * prop[1:,:])<0 ,1., 0.)
    varl = var[:-1,:] * zc
    varh = var[1:,:] * zc
    propl = prop[:-1,:] * zc
    proph = prop[1:,:] * zc
    isoslice = varl - propl * (varh - varl) / (proph - propl)
    isoslice = np.where(zc==1., isoslice, 0.)
    szc = zc.sum(axis=0)
    szc = np.where(szc==0., 1, szc)
    isoslice = isoslice.sum(axis=0)/szc
    if masking:
        isoslice = np.ma.masked_where(zc.sum(axis=0)==0, isoslice)
        if all(isoslice.mask):
            raise Warning('property==%f out of range (%f, %f)' % \
                           (isoval, (prop+isoval).min(), (prop+isoval).max()))
    isoslice = isoslice.reshape(sz[1:])

    # mask land
    isoslice = np.ma.masked_where(mask == 0, isoslice)

    return isoslice, x, y



def transect(var, istart, iend, jstart, jend, grd, Cpos='rho', vert=False, \
            spval=1e37):
    """
    transect, z, lon, lat = transect(var, istart, iend, jstart, jend, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'       specify the C-grid position where
                                     the variable rely
      - vert=True/False              If True, return the position of
                                     the verticies
      - spval                        special value
      - rtol                         tolerance parameter


    return a vertical transect between the points P1=(istart, jstart)
    and P2=(iend, jend) from 3D variable var
    lon, lat and z contain the C-grid position of the section for plotting.
    If vert=True, lon, lat and z contain contain the position of the
    verticies (to be used with pcolor)
    """

    # compute the depth on Arakawa-C grid position and get grid information

    if Cpos == 'u':
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = np.concatenate((z[:,0:1,:], z, z[:,-1:,:]), 1)
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
                y = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:,:-1] + grd.hgrid.x_vert[:,1:])
                y = 0.5 * (grd.hgrid.y_vert[:,:-1] + grd.hgrid.y_vert[:,1:])
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_u[:]
                y = grd.hgrid.lat_u[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_u[:]
                y = grd.hgrid.y_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos == 'v':
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-1:]), 2)
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
                y = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:-1,:] + grd.hgrid.x_vert[1:,:])
                y = 0.5 * (grd.hgrid.y_vert[:-1,:] + grd.hgrid.y_vert[1:,:])
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,-1:,:] + z[:,1:,:])
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_v[:]
                y = grd.hgrid.lat_v[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_v[:]
                y = grd.hgrid.y_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos == 'rho':
        # for temp, salt, rho
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-1:]), 2)
            z = np.concatenate((z[:,0:1,:], z, z[:,-1:,:]), 1)
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[:]
                y = grd.hgrid.lat_vert[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.lon_vert[:]
                y = grd.hgrid.lat_vert[:]
        else:
            z = grd.vgrid.z_r[0,:]
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)


    # Find the nearest point between P1 (imin,jmin) and P2 (imax, jmax)
    # -----------------------------------------------------------------
    # Initialization
    i0=istart; j0=jstart; i1=iend;  j1=jend
    istart = float(istart); iend = float(iend)
    jstart = float(jstart); jend = float(jend)

    # Compute equation:  j = aj i + bj
    if istart != iend:
        aj = (jend - jstart ) / (iend - istart)
        bj = jstart - aj * istart
    else:
        aj=10000.
        bj=0.

    # Compute equation:  i = ai j + bi
    if jstart != jend:
        ai = (iend - istart ) / ( jend - jstart )
        bi = istart - ai * jstart
    else:
        ai=10000.
        bi=0.

    # Compute the integer pathway:
    # Chose the strait line with the smallest slope
    if (abs(aj) <=  1 ):
        # Here, the best line is y(x)
        print('Here, the best line is y(x)')
        # If i1 < i0 swap points and remember it has been swapped
        if (i1 <  i0 ):
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j

        # compute the nearest j point on the line crossing at i
        n=0
        near = np.zeros(((i1-i0+1),4))
        for i in range(i0,i1+1):
            jj = aj*i + bj
            near[n,0] = i
            near[n,1] = jj
            near[n,2] = np.floor(jj)
            near[n,3] = np.ceil(jj)
            n = n + 1

        if vert == False:
            nearp = np.zeros(((i1-i0+1),4))
            nearp = near
        else:
            # compute the nearest j vert point on the line crossing at i
            n=0
            nearp = np.zeros(((i1-i0+2),4))
            for i in range(i0,i1+2):
                jj = aj*(i-0.5) + bj
                nearp[n,0] = i
                nearp[n,1] = jj
                nearp[n,2] = np.floor(jj)
                nearp[n,3] = np.ceil(jj)
                n = n + 1

    else:
        # Here, the best line is x(y)
        print('Here, the best line is x(y)')
        # If j1 < j0 swap points and remember it has been swapped
        if (j1 <  j0 ):
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j

        # compute the nearest i point on the line crossing at j
        n=0
        near = np.zeros(((j1-j0+1),4))
        for j in range(j0,j1+1):
            ii = ai*j + bi
            near[n,0] = j
            near[n,1] = ii
            near[n,2] = np.floor(ii)
            near[n,3] = np.ceil(ii)
            n = n + 1

        if vert == False:
            nearp = np.zeros(((j1-j0+1),4))
            nearp = near
        else:
            # compute the nearest i vert point on the line crossing at j
            n=0
            nearp = np.zeros(((j1-j0+2),4))
            for j in range(j0,j1+2):
                ii = ai*(j-0.5) + bi
                nearp[n,0] = j
                nearp[n,1] = ii
                nearp[n,2] = np.floor(ii)
                nearp[n,3] = np.ceil(ii)
                n = n + 1


    # Now interpolate between the nearest point through the section
    # -------------------------------------------------------------
    # Initialize output variables
    nlev = z.shape[0]

    transect = np.zeros((grd.vgrid.N, near.shape[0]))
    zs = np.zeros((nlev, nearp.shape[0]))
    xs = np.zeros((nlev, nearp.shape[0]))
    ys = np.zeros((nlev, nearp.shape[0]))

    # mask variable
    for k in range(var.shape[0]):
        var[k,:,:] = np.ma.masked_where(mask == 0, var[k,:,:])

    for n in range(near.shape[0]):
        if (abs(aj) <=  1 ):
            # check if our position match a grid cell
            if (near[n,2] == near[n,3]):
                transect[:,n] = var[:, near[n,2], near[n,0]]
            else:
                if mask[near[n,3], near[n,0]] == 0 or mask[near[n,2], near[n,0]] == 0:
                    transect[:,n] = spval
                else:
                    transect[:,n] = (near[n,1] - near[n,2]) * var[:, near[n,3], near[n,0]] + \
                                   (near[n,3] - near[n,1]) * var[:, near[n,2], near[n,0]]
        else:
            # check if our position match a grid cell
            if (near[n,2] == near[n,3]):
                transect[:,n] = var[:, near[n,0], near[n,2]]
            else:
                if mask[near[n,0], near[n,3]] == 0 or mask[near[n,0], near[n,2]] == 0:
                    transect[:,n] = spval
                else:
                    transect[:,n] = (near[n,1] - near[n,2]) * var[:, near[n,0], near[n,3]] + \
                                   (near[n,3] - near[n,1]) * var[:, near[n,0], near[n,2]]

    for n in range(nearp.shape[0]):
        if (abs(aj) <=  1 ):
            # check if our position match a grid cell
            if (nearp[n,2] == nearp[n,3]):
                zs[:,n] = z[:, nearp[n,2], nearp[n,0]]
                xs[:,n] = x[nearp[n,2], nearp[n,0]]
                ys[:,n] = y[nearp[n,2], nearp[n,0]]
            else:
                zs[:,n] = (nearp[n,1] - nearp[n,2]) * z[:, nearp[n,3], nearp[n,0]] + \
                          (nearp[n,3] - nearp[n,1]) * z[:, nearp[n,2], nearp[n,0]]
                xs[:,n] = (nearp[n,1] - nearp[n,2]) * x[nearp[n,3], nearp[n,0]] + \
                            (nearp[n,3] - nearp[n,1]) * x[nearp[n,2], nearp[n,0]]
                ys[:,n] = (nearp[n,1] - nearp[n,2]) * y[nearp[n,3], nearp[n,0]] + \
                            (nearp[n,3] - nearp[n,1]) * y[nearp[n,2], nearp[n,0]]
        else:
            # check if our position match a grid cell
            if (nearp[n,2] == nearp[n,3]):
                zs[:,n] = z[:, nearp[n,0], nearp[n,2]]
                xs[:,n] = x[nearp[n,0], nearp[n,2]]
                ys[:,n] = y[nearp[n,0], nearp[n,2]]
            else:
                zs[:,n] = (nearp[n,1] - nearp[n,2]) * z[:, nearp[n,0], nearp[n,3]] + \
                          (nearp[n,3] - nearp[n,1]) * z[:, nearp[n,0], nearp[n,2]]
                xs[:,n] = (nearp[n,1] - nearp[n,2]) * x[nearp[n,0], nearp[n,3]] + \
                            (nearp[n,3] - nearp[n,1]) * x[nearp[n,0], nearp[n,2]]
                ys[:,n] = (nearp[n,1] - nearp[n,2]) * y[nearp[n,0], nearp[n,3]] + \
                            (nearp[n,3] - nearp[n,1]) * y[nearp[n,0], nearp[n,2]]

    # mask transect
    transect = np.ma.masked_values(transect, spval)


    return transect, zs, xs, ys




def lonslice(var, longitude, grd, Cpos='rho', vert=False, spval=1e37):
    """
    lonslice, z, lon, lat = lonslice(var, longitude, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'       specify the C-grid position where
                                     the variable rely
      - vert=True/False              If True, return the position of
                                     the verticies
      - spval                        special value
      - rtol                         tolerance parameter


    return a longitudinal slice along longitude=longitude from 3D variable var
    lon, lat and z contain the C-grid position of the section for plotting.
    If vert=True, lon, lat and z contain contain the position of the
    verticies (to be used with pcolor)
    Returns a longitudinal slice of the grid
    """

    if Cpos == 'u':
        lon = grd.hgrid.lon_u
        lat = grd.hgrid.lat_u
    elif Cpos == 'v':
        lon = grd.hgrid.lon_v
        lat = grd.hgrid.lat_v
    elif Cpos == 'rho':
        lon = grd.hgrid.lon_rho
        lat = grd.hgrid.lat_rho
    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)

    edge = np.concatenate((lon[1,1:-1], \
                           lon[1:-1,-2], \
                           lon[-2,-2:0:-1], \
                           lon[-2:0:-1,1]))
    idx =  np.concatenate((list(range(1,lon[0,:].shape[0]-1)), \
                           list(range(1,lon[:,-1].shape[0]-1)), \
                           list(range(1,lon[-1,::-1].shape[0]-1))[::-1], \
                           list(range(1,lon[::-1,0].shape[0]-1))[::-1]))

    d = np.zeros(edge.shape)
    for i in range (edge.shape[0]):
        d[i] = edge[i] - longitude
        d[i] = d[i] / abs(d[i])
    d = np.diff(d)

    pt_idx = np.where(d != 0)[0]

    Mp, Lp = lon.shape

    if len(pt_idx) != 2:
        raise ValueError('this function only works for simple quadrangle')

    # determine if latitude line is crossing a i or j edge
    side = np.zeros(2)
    if pt_idx[0] < Lp: side[0] = 1
    if pt_idx[0] >= Lp and pt_idx[0] < Lp+Mp: side[0] = 2
    if pt_idx[0] >= Lp+Mp and pt_idx[0] < Lp+Mp+Lp: side[0] = 3
    if pt_idx[0] >= Lp+Mp+Lp: side[0] = 4
    if pt_idx[1] < Lp: side[1] = 1
    if pt_idx[1] >= Lp and pt_idx[1] < Lp+Mp: side[1] = 2
    if pt_idx[1] >= Lp+Mp and pt_idx[1] < Lp+Mp+Lp: side[1] = 3
    if pt_idx[1] >= Lp+Mp+Lp: side[1] = 4


    if side[0] == 1 and side[1] == 2:
        lonslice, z, lon, lat = pyroms.tools.transect(var, \
                                idx[pt_idx[0]], Lp-2, \
                                1, idx[pt_idx[1]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 1 and side[1] == 3:
        lonslice, z, lon, lat = pyroms.tools.transect(var, \
                                idx[pt_idx[0]], idx[pt_idx[1]], \
                                1, Mp-2, \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 1 and side[1] == 4:
        lonslice, z, lon, lat = pyroms.tools.transect(var, \
                                idx[pt_idx[0]], 1, \
                                1, idx[pt_idx[1]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 2 and side[1] == 3:
        lonslice, z, lon, lat = pyroms.tools.transect(var, \
                                Lp-2, idx[pt_idx[1]], \
                                idx[pt_idx[0]], Mp-2, \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 2 and side[1] == 4:
        lonslice, z, lon, lat = pyroms.tools.transect(var, \
                                Lp-2, 1, \
                                idx[pt_idx[0]], idx[pt_idx[0]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 3 and side[1] == 4:
        lonslice, z, lon, lat = pyroms.tools.transect(var, \
                                idx[pt_idx[0]], 1, \
                                Mp-2, idx[pt_idx[1]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)


    return lonslice, z, lon, lat



def latslice(var, latitude, grd, Cpos='rho', vert=False, spval=1e37):
    """
    latslice, z, lon, lat = latslice(var, latitude, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'       specify the C-grid position where
                                     the variable rely
      - vert=True/False              If True, return the position of
                                     the verticies
      - spval                        special value
      - rtol                         tolerance parameter


    return a latitudinal slice along latitude=latitude from 3D variable var
    lon, lat and z contain the C-grid position of the section for plotting.
    If vert=True, lon, lat and z contain contain the position of the
    verticies (to be used with pcolor)
    Returns a longitudinal slice of the grid
    """

    if Cpos == 'u':
        lon = grd.hgrid.lon_u
        lat = grd.hgrid.lat_u
    elif Cpos == 'v':
        lon = grd.hgrid.lon_v
        lat = grd.hgrid.lat_v
    elif Cpos == 'rho':
        lon = grd.hgrid.lon_rho
        lat = grd.hgrid.lat_rho
    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)

    edge = np.concatenate((lat[1,1:-1], \
                           lat[1:-1,-2], \
                           lat[-2,-2:0:-1], \
                           lat[-2:0:-1,1]))
    idx =  np.concatenate((list(range(1,lat[0,:].shape[0]-1)), \
                           list(range(1,lat[:,-1].shape[0]-1)), \
                           list(range(1,lat[-1,::-1].shape[0]-1))[::-1], \
                           list(range(1,lat[::-1,0].shape[0]-1))[::-1]))

    d = np.zeros(edge.shape)
    for i in range (edge.shape[0]):
        d[i] = edge[i] - latitude
        d[i] = d[i] / abs(d[i])
    d = np.diff(d)

    pt_idx = np.where(d != 0)[0]

    Mp, Lp = lon.shape

    if len(pt_idx) != 2:
        raise ValueError('this function only works for simple quadrangle')

    # determine if latitude line is crossing a i or j edge
    side = np.zeros(2)
    if pt_idx[0] < Lp: side[0] = 1
    if pt_idx[0] >= Lp and pt_idx[0] < Lp+Mp: side[0] = 2
    if pt_idx[0] >= Lp+Mp and pt_idx[0] < Lp+Mp+Lp: side[0] = 3
    if pt_idx[0] >= Lp+Mp+Lp: side[0] = 4
    if pt_idx[1] < Lp: side[1] = 1
    if pt_idx[1] >= Lp and pt_idx[1] < Lp+Mp: side[1] = 2
    if pt_idx[1] >= Lp+Mp and pt_idx[1] < Lp+Mp+Lp: side[1] = 3
    if pt_idx[1] >= Lp+Mp+Lp: side[1] = 4


    if side[0] == 1 and side[1] == 2:
        latslice, z, lon, lat = pyroms.tools.section(var, \
                                idx[pt_idx[0]], Lp-2, \
                                1, idx[pt_idx[1]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 1 and side[1] == 3:
        latslice, z, lon, lat = pyroms.tools.section(var, \
                                idx[pt_idx[0]], idx[pt_idx[1]], \
                                1, Mp-2, \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 1 and side[1] == 4:
        latslice, z, lon, lat = pyroms.tools.section(var, \
                                idx[pt_idx[0]], 1, \
                                1, idx[pt_idx[1]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 2 and side[1] == 3:
        latslice, z, lon, lat = pyroms.tools.section(var, \
                                Lp-2, idx[pt_idx[1]], \
                                idx[pt_idx[0]], Mp-2, \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 2 and side[1] == 4:
        latslice, z, lon, lat = pyroms.tools.section(var, \
                                Lp-2, 1, \
                                idx[pt_idx[0]], idx[pt_idx[0]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 3 and side[1] == 4:
        latslice, z, lon, lat = pyroms.tools.section(var, \
                                idx[pt_idx[0]], 1, \
                                Mp-2, idx[pt_idx[1]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)

    return latslice, z, lon, lat


def zlayer(var, grd, h1=None, h2=None, Cpos='rho', vert=False):
    """
    extract a z layer between depth h1 and h2
    """

    # compute the depth on Arakawa-C grid position

    if Cpos == 'u':
        # average z_r at Arakawa-C u points
        z = 0.5 * (grd.vgrid.z_w[0,:,:,:-1] + grd.vgrid.z_w[0,:,:,1:])
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
                y = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:,:-1] + grd.hgrid.x_vert[:,1:])
                y = 0.5 * (grd.hgrid.y_vert[:,:-1] + grd.hgrid.y_vert[:,1:])
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_u[:]
                y = grd.hgrid.lat_u[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_u[:]
                y = grd.hgrid.y_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos == 'v':
        # average z_r at Arakawa-C v points
        z = 0.5 * (grd.vgrid.z_w[0,:,:-1,:] + grd.vgrid.z_w[0,:,1:,:])
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
                y = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
            elif grd.hgrid.spherical == 'F':
                x = 0.5 * (grd.hgrid.x_vert[:-1,:] + grd.hgrid.x_vert[1:,:])
                y = 0.5 * (grd.hgrid.y_vert[:-1,:] + grd.hgrid.y_vert[1:,:])
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_v[:]
                y = grd.hgrid.lat_v[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_v[:]
                y = grd.hgrid.y_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos == 'rho':
        # for temp, salt, rho
        z = grd.vgrid.z_w[0,:]
        if vert == True:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_vert[:]
                y = grd.hgrid.lat_vert[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_vert[:]
                y = grd.hgrid.y_vert[:]
        else:
            if grd.hgrid.spherical == 'T':
                x = grd.hgrid.lon_rho[:]
                y = grd.hgrid.lat_rho[:]
            elif grd.hgrid.spherical == 'F':
                x = grd.hgrid.x_rho[:]
                y = grd.hgrid.y_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)



    #set var to zero where tra is masked for the sum
    mask = np.tile(mask, (var.shape[0],1,1))
    var = np.where(mask == 1, var, 0)

    zlayer = np.zeros((var.shape[1], var.shape[2]))
    dz = z[1:,:,:] - z[:-1,:,:]

    for jj in range(var.shape[1]):
        for ji in range(var.shape[2]):
            ratio = np.ones(var.shape[0])
            if h1 is not None:
                if np.abs(h1) < np.abs(z[0,jj,ji]):
                    idx1 = np.where(z[:,jj,ji] > -np.abs(h1))
                    if np.any(idx1):
                        idx1 = idx1[0][0]
                        r = (-np.abs(h1) - z[idx1-1,jj,ji]) / \
                                (z[idx1,jj,ji] -z[idx1-1,jj,ji])
                        ratio[idx1:] = 0.
                        ratio[idx1-1] = r
                else:
                    ratio[:] = 0.
            if h2 is not None:
                idx2 = np.where(z[:,jj,ji] < -np.abs(h2))
                if np.any(idx2):
                    idx2 = idx2[0][-1]
                    r = (z[idx2+1,jj,ji] - (-np.abs(h2))) / \
                                 (z[idx2+1,jj,ji] -z[idx2,jj,ji])
                    ratio[:idx2] = 0.
                    ratio[idx2] = r
            zlayer[jj,ji] = np.sum(var[:,jj,ji] * dz[:,jj,ji] * ratio[:])/ \
                           np.sum(dz[:,jj,ji] * ratio[:])

    return zlayer, x, y


def section_transport(u, v, grd, istart, iend, jstart, jend):
    """
    transpu, transpv = section_transport(u, v, grd, istart, iend, jstart, jend)

    compute the transport through the section defined between
    the point P1 (istart,jstart) and P2 (iend, jend).
    P1 and P2 are Arakawa-C psi points.
    The transpot is positive right handside of the section.
    """


    # Find the nearest point between P1 (imin,jmin) and P2 (imax, jmax)
    # -----------------------------------------------------------------
    # Initialization
    i0=istart; j0=jstart; i1=iend;  j1=jend
    istart = float(istart); iend = float(iend)
    jstart = float(jstart); jend = float(jend)

    # Compute equation:  j = aj i + bj
    if istart != iend:
        aj = (jend - jstart ) / (iend - istart)
        bj = jstart - aj * istart
    else:
        aj=10000.
        bj=0.

    # Compute equation:  i = ai j + bi
    if jstart != jend:
        ai = (iend - istart ) / ( jend - jstart )
        bi = istart - ai * jstart
    else:
        ai=10000.
        bi=0.

    # Compute the integer pathway:
    # Chose the strait line with the smallest slope
    if (abs(aj) <=  1 ):
        # Here, the best line is y(x)
        print('Here, the best line is y(x)')
        # If i1 < i0 swap points and remember it has been swapped
        if i1 <  i0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if j1 >= j0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 1; jst = 0
            norm_u = -1; norm_v = -1

        near = []
        # compute the nearest j point on the line crossing at i
        for i in range(i0,i1+1):
            j = aj*i + bj
            near.append(i + round(j)*1j)

    else:
        # Here, the best line is x(y)
        print('Here, the best line is x(y)')
        # If j1 < j0 swap points and remember it has been swapped
        if j1 <  j0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if i1 >= i0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 0; jst = 1
            norm_u = 1; norm_v = 1

        near = []
        # compute the nearest i point on the line crossing at j
        for j in range(j0,j1+1):
            i = ai*j + bi
            near.append(round(i) + j*1j)


    # Look for intermediate points to be added
    # -------------------------------------------------------------

    inear = np.copy(near)

    n = len(near)
    nn=1

    for k in range(1,n):
        # distance between 2 neighbour points
        d = abs(inear[k] - inear[k-1])

        if ( d > 1 ):
            # intermediate points required if d>1
            neari = interm_pt(inear, k, ai, bi, aj, bj)
            near.insert(nn,neari)
            nn=nn+1

        nn=nn+1


    # Now extract the transport through a section
    # -------------------------------------------

    #get metrics
    dx = grd.hgrid.dx
    dy = grd.hgrid.dy
    z_w = grd.vgrid.z_w[0,:]
    # average z_w at Arakawa-C u points
    zu = 0.5 * (z_w[:,:,:-1] + z_w[:,:,1:])
    dzu = zu[1:,:,:] - zu[:-1,:,:]
    # average z_w at Arakawa-C v points
    zv = 0.5 * (z_w[:,:-1,:] + z_w[:,1:,:])
    dzv = zv[1:,:,:] - zv[:-1,:,:]

    #set u and v to zero where u and v are masked for the sum
    for k in range(u.shape[0]):
        u[k,:] = np.where(grd.hgrid.mask_u == 1, u[k,:], 0)
        v[k,:] = np.where(grd.hgrid.mask_v == 1, v[k,:], 0)

    n = len(near)
    transpu = 0
    transpv = 0

    for l in range(0,n-1):
        ii = int(np.real(near[l])); jj = int(np.imag(near[l]))
        for k in range(0, dzu.shape[0]):
            if np.real(near[l]) == np.real(near[l+1]):
                trans = u[k, jj+jst, ii] * dy[jj+jst, ii] * \
                        dzu[k, jj+jst, ii] * norm_u * norm
                transpu = transpu + trans

            elif np.imag(near[l]) == np.imag(near[l+1]):
                trans = v[k, jj, ii+ist] * dx[jj, ii+ist] * \
                        dzv[k, jj, ii+ist] * norm_v * norm
                transpv = transpv + trans


    return transpu, transpv


def section_transport_z(u, v, grd, istart, iend, jstart, jend, h1=None, h2=None):
    """
    transpu, transpv = section_transport(u, v, grd, istart, iend, jstart, jend, h1, h2)

    compute the transport between depth h1 and h2 (if specified) through the section
    defined between the point P1 (istart,jstart) and P2 (iend, jend).
    P1 and P2 are Arakawa-C psi points.
    The transpot is positive right handside of the section.
    """


    # Find the nearest point between P1 (imin,jmin) and P2 (imax, jmax)
    # -----------------------------------------------------------------
    # Initialization
    i0=istart; j0=jstart; i1=iend;  j1=jend
    istart = float(istart); iend = float(iend)
    jstart = float(jstart); jend = float(jend)

    # Compute equation:  j = aj i + bj
    if istart != iend:
        aj = (jend - jstart ) / (iend - istart)
        bj = jstart - aj * istart
    else:
        aj=10000.
        bj=0.

    # Compute equation:  i = ai j + bi
    if jstart != jend:
        ai = (iend - istart ) / ( jend - jstart )
        bi = istart - ai * jstart
    else:
        ai=10000.
        bi=0.

    # Compute the integer pathway:
    # Chose the strait line with the smallest slope
    if (abs(aj) <=  1 ):
        # Here, the best line is y(x)
        print('Here, the best line is y(x)')
        # If i1 < i0 swap points and remember it has been swapped
        if i1 <  i0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if j1 >= j0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 1; jst = 0
            norm_u = -1; norm_v = -1

        near = []
        # compute the nearest j point on the line crossing at i
        for i in range(i0,i1+1):
            j = aj*i + bj
            near.append(i + round(j)*1j)

    else:
        # Here, the best line is x(y)
        print('Here, the best line is x(y)')
        # If j1 < j0 swap points and remember it has been swapped
        if j1 <  j0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if i1 >= i0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 0; jst = 1
            norm_u = 1; norm_v = 1

        near = []
        # compute the nearest i point on the line crossing at j
        for j in range(j0,j1+1):
            i = ai*j + bi
            near.append(round(i) + j*1j)


    # Look for intermediate points to be added
    # -------------------------------------------------------------

    inear = np.copy(near)

    n = len(near)
    nn=1

    for k in range(1,n):
        # distance between 2 neighbour points
        d = abs(inear[k] - inear[k-1])

        if ( d > 1 ):
            # intermediate points required if d>1
            neari = interm_pt(inear, k, ai, bi, aj, bj)
            near.insert(nn,neari)
            nn=nn+1

        nn=nn+1


    # Now extract the transport through a section
    # -------------------------------------------

    #get metrics
    dx = grd.hgrid.dx
    dy = grd.hgrid.dy
    z_w = grd.vgrid.z_w[0,:]
    # average z_w at Arakawa-C u points
    zu = 0.5 * (z_w[:,:,:-1] + z_w[:,:,1:])
    dzu = zu[1:,:,:] - zu[:-1,:,:]
    # average z_w at Arakawa-C v points
    zv = 0.5 * (z_w[:,:-1,:] + z_w[:,1:,:])
    dzv = zv[1:,:,:] - zv[:-1,:,:]

    #set u and v to zero where u and v are masked for the sum
    for k in range(u.shape[0]):
        u[k,:] = np.where(grd.hgrid.mask_u == 1, u[k,:], 0)
        v[k,:] = np.where(grd.hgrid.mask_v == 1, v[k,:], 0)

    n = len(near)
    transpu = 0
    transpv = 0

    for l in range(0,n-1):

        ii = int(np.real(near[l])); jj = int(np.imag(near[l]))

        if np.real(near[l]) == np.real(near[l+1]):
            vel = u[:, jj+jst, ii]
            ds = dy[jj+jst, ii]
            z = zu[:, jj+jst, ii]
            dz = dzu[:, jj+jst, ii]
            norm_su = norm_u * norm
            norm_sv = 0.
        elif np.imag(near[l]) == np.imag(near[l+1]):
            vel = v[:, jj, ii+ist]
            ds = dx[jj, ii+ist]
            z = zv[:, jj, ii+ist]
            dz = dzv[:, jj, ii+ist]
            norm_su = 0.
            norm_sv = norm_v * norm

        ratio = np.ones(dz.shape)
        if h1 is not None:
            if np.abs(h1) < np.abs(z[0]):
                idx1 = np.where(z > -np.abs(h1))
                if np.any(idx1):
                    idx1 = idx1[0][0]
                    r = (-np.abs(h1) - z[idx1-1]) / (z[idx1] -z[idx1-1])
                    ratio[idx1:] = 0.
                    ratio[idx1-1] = r
            else:
                ratio[:] = 0.
        if h2 is not None:
            idx2 = np.where(z < -np.abs(h2))
            if np.any(idx2):
                idx2 = idx2[0][-1]
                r = (z[idx2+1] - (-np.abs(h2))) / (z[idx2+1] -z[idx2])
                ratio[:idx2] = 0.
                ratio[idx2] = r

        for k in range(0, dzu.shape[0]):
            transu = vel[k] * ds * dz[k] * norm_su * ratio[k]
            transv = vel[k] * ds * dz[k] * norm_sv * ratio[k]
            transpu = transpu + transu
            transpv = transpv + transv

    return transpu, transpv


def section_tracer_transport_z(u, v, tracer, grd, istart, iend, jstart, jend, h1=None, h2=None):
    """
    transpu, transpv = section_transport(u, v, tracer, grd, istart, iend, jstart, jend, h1, h2)

    compute the tracer transport between depth h1 and h2 (if specified) through the section
    defined between the point P1 (istart,jstart) and P2 (iend, jend).
    P1 and P2 are Arakawa-C psi points.
    The transpot is positive right handside of the section.
    """


    # Find the nearest point between P1 (imin,jmin) and P2 (imax, jmax)
    # -----------------------------------------------------------------
    # Initialization
    i0=istart; j0=jstart; i1=iend;  j1=jend
    istart = float(istart); iend = float(iend)
    jstart = float(jstart); jend = float(jend)

    # Compute equation:  j = aj i + bj
    if istart != iend:
        aj = (jend - jstart ) / (iend - istart)
        bj = jstart - aj * istart
    else:
        aj=10000.
        bj=0.

    # Compute equation:  i = ai j + bi
    if jstart != jend:
        ai = (iend - istart ) / ( jend - jstart )
        bi = istart - ai * jstart
    else:
        ai=10000.
        bi=0.

    # Compute the integer pathway:
    # Chose the strait line with the smallest slope
    if (abs(aj) <=  1 ):
        # Here, the best line is y(x)
        print('Here, the best line is y(x)')
        # If i1 < i0 swap points and remember it has been swapped
        if i1 <  i0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if j1 >= j0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 1; jst = 0
            norm_u = -1; norm_v = -1

        near = []
        # compute the nearest j point on the line crossing at i
        for i in range(i0,i1+1):
            j = aj*i + bj
            near.append(i + round(j)*1j)

    else:
        # Here, the best line is x(y)
        print('Here, the best line is x(y)')
        # If j1 < j0 swap points and remember it has been swapped
        if j1 <  j0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if i1 >= i0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 0; jst = 1
            norm_u = 1; norm_v = 1

        near = []
        # compute the nearest i point on the line crossing at j
        for j in range(j0,j1+1):
            i = ai*j + bi
            near.append(round(i) + j*1j)


    # Look for intermediate points to be added
    # -------------------------------------------------------------

    inear = np.copy(near)

    n = len(near)
    nn=1

    for k in range(1,n):
        # distance between 2 neighbour points
        d = abs(inear[k] - inear[k-1])

        if ( d > 1 ):
            # intermediate points required if d>1
            neari = interm_pt(inear, k, ai, bi, aj, bj)
            near.insert(nn,neari)
            nn=nn+1

        nn=nn+1


    # Now extract the transport through a section
    # -------------------------------------------

    #get metrics
    dx = grd.hgrid.dx
    dy = grd.hgrid.dy
    z_w = grd.vgrid.z_w[0,:]
    # average z_w at Arakawa-C u points
    zu = 0.5 * (z_w[:,:,:-1] + z_w[:,:,1:])
    dzu = zu[1:,:,:] - zu[:-1,:,:]
    # average z_w at Arakawa-C v points
    zv = 0.5 * (z_w[:,:-1,:] + z_w[:,1:,:])
    dzv = zv[1:,:,:] - zv[:-1,:,:]

    #tracer value at u and v position
    trau = 0.5 * (tracer[:,:,1:] + tracer[:,:,:-1])
    trav = 0.5 * (tracer[:,1:,:] + tracer[:,:-1,:])

    #set u and v to zero where u and v are masked for the sum
    for k in range(u.shape[0]):
        u[k,:] = np.where(grd.hgrid.mask_u == 1, u[k,:], 0)
        v[k,:] = np.where(grd.hgrid.mask_v == 1, v[k,:], 0)
        trau[k,:] = np.where(grd.hgrid.mask_u == 1, trau[k,:], 0)
        trav[k,:] = np.where(grd.hgrid.mask_v == 1, trav[k,:], 0)

    n = len(near)
    transpu = 0
    transpv = 0

    for l in range(0,n-1):

        ii = int(np.real(near[l])); jj = int(np.imag(near[l]))

        if np.real(near[l]) == np.real(near[l+1]):
            vel = u[:, jj+jst, ii]
            tra = trau[:, jj+jst, ii]
            ds = dy[jj+jst, ii]
            z = zu[:, jj+jst, ii]
            dz = dzu[:, jj+jst, ii]
            norm_su = norm_u * norm
            norm_sv = 0.
        elif np.imag(near[l]) == np.imag(near[l+1]):
            vel = v[:, jj, ii+ist]
            tra = trav[:, jj, ii+ist]
            ds = dx[jj, ii+ist]
            z = zv[:, jj, ii+ist]
            dz = dzv[:, jj, ii+ist]
            norm_su = 0.
            norm_sv = norm_v * norm

        ratio = np.ones(dz.shape)
        if h1 is not None:
            if np.abs(h1) < np.abs(z[0]):
                idx1 = np.where(z > -np.abs(h1))
                if np.any(idx1):
                    idx1 = idx1[0][0]
                    r = (-np.abs(h1) - z[idx1-1]) / (z[idx1] -z[idx1-1])
                    ratio[idx1:] = 0.
                    ratio[idx1-1] = r
            else:
                ratio[:] = 0.
        if h2 is not None:
            idx2 = np.where(z < -np.abs(h2))
            if np.any(idx2):
                idx2 = idx2[0][-1]
                r = (z[idx2+1] - (-np.abs(h2))) / (z[idx2+1] -z[idx2])
                ratio[:idx2] = 0.
                ratio[idx2] = r

        for k in range(0, dzu.shape[0]):
            transu = vel[k] * tra[k] * ds * dz[k] * norm_su * ratio[k]
            transv = vel[k] * tra[k] * ds * dz[k] * norm_sv * ratio[k]
            transpu = transpu + transu
            transpv = transpv + transv

    return transpu, transpv



def interm_pt(pnear, pk, pai, pbi, paj, pbj):
    ### FIND THE BEST INTERMEDIATE POINT ON A PATHWAY
    # -----------------------------
    # pnear   : vector of the position of the nearest point
    # pk      : current working index
    # pai, pbi: slope and original ordinate of x(y)
    # paj, pbj: slope and original ordinate of y(x)
    # pneari  : vector holding the position of intermediate point
    # -----------------------------

    # 1 - Compute intermediate point

    # Determine whether we use y(x) or x(y):
    if (abs(paj) <= 1):
        # y(x)
        # possible intermediate point
        ylptmp1 = pnear[pk-1] + 1
        ylptmp2 = pnear[pk-1] + (paj/abs(paj))*1j
        # M is the candidate point:
        zxm = np.real(ylptmp1)
        zym = np.imag(ylptmp1)
        za0 = paj
        zb0 = pbj
        #
        za1 = -1./za0
        zb1 = zym-za1*zxm
        # P is the projection of M in the strait line
        zxp = -(zb1-zb0)/(za1-za0)
        zyp = za0*zxp+zb0
        # zd1 is the distance MP
        zd1 = (zxm-zxp) * (zxm-zxp) + (zym-zyp) * (zym-zyp)
        #
        # M is the candidate point:
        zxm = np.real(ylptmp2)
        zym = np.imag(ylptmp2)
        za1 = -1./za0
        zb1 = zym-za1*zxm
        # P is the projection of M in the strait line
        zxp = -(zb1-zb0)/(za1-za0)
        zyp = za0*zxp+zb0
        # zd1 is the distance MP
        zd2 = (zxm-zxp) * (zxm-zxp) + (zym-zyp) * (zym-zyp)
        #
        # choose the smallest (zd1,zd2)
        if (zd2 <= zd1):
            pneari = ylptmp2
        else:
            pneari = ylptmp1
        #
    else:
        # x(y)
        ylptmp1 = pnear[pk-1] + (pai/abs(pai))
        ylptmp2 = pnear[pk-1] + 1*1j
        # M is the candidate point:
        zxm = np.real(ylptmp1)
        zym = np.imag(ylptmp1)
        za0 = pai
        zb0 = pbi
        #
        za1 = -1./za0
        zb1 = zxm-za1*zym
        # P is the projection of M in the strait line
        zyp = -(zb1-zb0)/(za1-za0)
        zxp = za0*zyp+zb0
        # zd1 is the distance MP
        zd1 = (zxm-zxp) * (zxm-zxp) + (zym-zyp) * (zym-zyp)
        #
        # M is the candidate point:
        zxm = np.real(ylptmp2)
        zym = np.imag(ylptmp2)
        za1 = -1./za0
        zb1 = zxm-za1*zym
        # P is the projection of M in the strait line
        zyp = -(zb1-zb0)/(za1-za0)
        zxp = za0*zyp+zb0
        # zd2 is the distance MP
        zd2 = (zxm-zxp) * (zxm-zxp) + (zym-zyp) * (zym-zyp)
        #
        # choose the smallest (zd1,zd2)
        if (zd2 <= zd1):
            pneari = ylptmp2
        else:
            pneari = ylptmp1

    return pneari


def hindices(lon, lat, grd, Cpos='rho', rectangular=0, spval=1e37):
    """
    """

    if type(grd).__name__ == 'ROMS_Grid':
        spherical = grd.hgrid.spherical
        if Cpos == 'u':
            long = grd.hgrid.lon_u
            latg = grd.hgrid.lat_u
            angle = 0.5 * (grd.hgrid.angle_rho[:,1:] + grd.hgrid.angle_rho[:,:-1])
        elif Cpos == 'v':
            long = grd.hgrid.lon_v
            latg = grd.hgrid.lat_v
            angle = 0.5 * (grd.hgrid.angle_rho[1:,:] + grd.hgrid.angle_rho[:-1,:])
        elif Cpos == 'rho':
            long = grd.hgrid.lon_rho
            latg = grd.hgrid.lat_rho
            angle = grd.hgrid.angle_rho
        else:
            raise Warning('%s bad position. Valid Arakawa-C are \
                               rho, u or v.' % Cpos)

    if type(grd).__name__ == 'CGrid_geo':
        spherical = grd.spherical
        if Cpos == 'u':
            long = grd.lon_u
            latg = grd.lat_u
            angle = 0.5 * (grd.angle_rho[:,1:] + grd.angle_rho[:,:-1])
        elif Cpos == 'v':
            long = grd.lon_v
            latg = grd.lat_v
            angle = 0.5 * (grd.angle_rho[1:,:] + grd.angle_rho[:-1,:])
        elif Cpos == 'rho':
            long = grd.lon_rho
            latg = grd.lat_rho
            angle = grd.angle_rho
        else:
            raise Warning('%s bad position. Valid Arakawa-C are \
                               rho, u or v.' % Cpos)


    lon = np.matrix(lon)
    lat = np.matrix(lat)

    ipos, jpos = _obs_interp.hindices(spherical, angle.T, long.T, latg.T, \
                                      lon, lat, spval, rectangular)

    # python indexing start with zero...
    ipos = np.squeeze(ipos) - 1
    jpos = np.squeeze(jpos) - 1

    ipos = np.ma.masked_where(ipos == spval, ipos)
    jpos = np.ma.masked_where(jpos == spval, jpos)

    return ipos, jpos


def obs_interp2d(Finp, lon, lat, grd, Cpos='rho', rectangular=0, spval=1e37):
    """
    """

    Iout, Jout = hindices(lon, lat, grd, Cpos=Cpos, rectangular=rectangular, spval=spval)

    # fortran indexing start with one...
    Iout = Iout + 1
    Jout = Jout + 1

    Fout = _obs_interp.linterp2d(Finp.T,Iout,Jout)

    return Fout
