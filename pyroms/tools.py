# encoding: utf-8

import numpy as np
import _iso

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
        raise Warning, '%s not supported, defaulting to linear' % mode
        

    # compute the depth on Arakawa-C grid position

    if Cpos is 'u':
        # average z_r at Arakawa-C u points
        z = 0.5 * (grd.vgrid.z_r[0,:,:,:-1] + grd.vgrid.z_r[0,:,:,1:])
        if vert == True: 
            lon = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
            lat = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
        else:
            lon = grd.hgrid.lon_u[:]
            lat = grd.hgrid.lat_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos is 'v':
        # average z_r at Arakawa-C v points
        z = 0.5 * (grd.vgrid.z_r[0,:,:-1,:] + grd.vgrid.z_r[0,:,1:,:])
        if vert == True: 
            lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
            lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        else:
            lon = grd.hgrid.lon_v[:]
            lat = grd.hgrid.lat_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos is 'w':
        # for temp, salt, rho
        z = grd.vgrid.z_w[0,:]
        if vert == True:
            lon = grd.hgrid.lon_vert[:]
            lat = grd.hgrid.lat_vert[:]
        else:
            lon = grd.hgrid.lon_rho[:]
            lat = grd.hgrid.lat_rho[:]
        mask = grd.hgrid.mask_rho[:]

    elif Cpos is 'rho':
        # for temp, salt, rho
        z = grd.vgrid.z_r[0,:]
        if vert == True: 
            lon = grd.hgrid.lon_vert[:]
            lat = grd.hgrid.lat_vert[:]
        else:
            lon = grd.hgrid.lon_rho[:]
            lat = grd.hgrid.lat_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning, '%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos

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

    return zslice, lon, lat


def sslice(var, sindex, grd, Cpos='rho', vert=False):
    """ 
    sslice, lon, lat = sslice(var, sindex, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'	     specify the C-grid position where 
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

    if Cpos is 'u':
        # average z_r at Arakawa-C u points
        z = 0.5 * (grd.vgrid.z_r[0,:,:,:-1] + grd.vgrid.z_r[0,:,:,1:])
        if vert == True: 
            lon = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
            lat = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
        else:
            lon = grd.hgrid.lon_u[:]
            lat = grd.hgrid.lat_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos is 'v':
        # average z_r at Arakawa-C v points
        z = 0.5 * (grd.vgrid.z_r[0,:,:-1,:] + grd.vgrid.z_r[0,:,1:,:])
        if vert == True: 
            lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
            lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        else:
            lon = grd.hgrid.lon_v[:]
            lat = grd.hgrid.lat_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos is 'rho':
        # for temp, salt, rho, w
        z = grd.vgrid.z_r[0,:]
        if vert == True: 
            lon = grd.hgrid.lon_vert[:]
            lat = grd.hgrid.lat_vert[:]
        else:
            lon = grd.hgrid.lon_rho[:]
            lat = grd.hgrid.lat_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning, '%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos

    assert len(var.shape) == 3, 'var must be 3D'

    sslice = var[sindex,:,:]

    # mask land
    sslice = np.ma.masked_where(mask == 0, sslice)
    
    return sslice, lon, lat



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

    if Cpos is 'u':
        # average z_r at Arakawa-C u points
        if vert == True: 
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = np.concatenate((z[:,0:1,:], z, z[:,-2:-1,:]), 1)
            lon = grd.hgrid.lon_vert[:,1:-1]
            lat = grd.hgrid.lat_vert[:,1:-1]
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            lon = grd.hgrid.lon_u[:]
            lat = grd.hgrid.lat_u[:]         
        mask = grd.hgrid.mask_u[:]

    elif Cpos is 'v':
        # average z_r at Arakawa-C v points
        if vert == True: 
            z = grd.vgrid.z_w[0,:]
            lon = grd.hgrid.lon_rho
            lat = grd.hgrid.lat_rho
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            lon = grd.hgrid.lon_v[:]
            lat = grd.hgrid.lat_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos is 'w':
        # for w, AKt, ...
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:-1,:,:] + z[1:,:,:])
            z = np.concatenate((np.array(grd.vgrid.z_w[0,0,:,:], ndmin=3), \
                                z, \
                                np.array(grd.vgrid.z_w[0,-1,:,:], ndmin=3)), 0)
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-2:-1]), 2)
            lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
            lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        else:
            z = grd.vgrid.z_w[0,:]
            lon = grd.hgrid.lon_rho[:]
            lat = grd.hgrid.lat_rho[:]
        mask = grd.hgrid.mask_rho[:]


    elif Cpos is 'rho':
        # for temp, salt, rho, ...
        if vert == True: 
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = np.concatenate((z[:,0:1,:], z, z[:,-2:-1,:]), 1)
            lon = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
            lat = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
        else:
            z = grd.vgrid.z_r[0,:]
            lon = grd.hgrid.lon_rho[:]
            lat = grd.hgrid.lat_rho[:]      
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning, '%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos

    # get constant-i slice
    vari = var[:,:,iindex]
    zi = z[:,:,iindex]
    loni = np.tile(lon[:,iindex], (zi.shape[0], 1))
    lati = np.tile(lat[:,iindex], (zi.shape[0], 1))

    # land/sea mask
    maski = np.tile(mask[:,iindex], (vari.shape[0], 1))
    vari = np.ma.masked_where(maski[:,:] == 0, vari[:,:])

    return vari, zi, loni, lati
    

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

    if Cpos is 'u':
        # average z_r at Arakawa-C u points
        if vert == True: 
            z = grd.vgrid.z_w[0,:]
            lon = grd.hgrid.lon_rho
            lat = grd.hgrid.lat_rho
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            lon = grd.hgrid.lon_u[:]
            lat = grd.hgrid.lat_u[:]         
        mask = grd.hgrid.mask_u[:]

    elif Cpos is 'v':
        # average z_r at Arakawa-C v points
        if vert == True: 
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-2:-1]), 2)
            lon = grd.hgrid.lon_vert[1:-1,:]
            lat = grd.hgrid.lat_vert[1:-1,:]
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            lon = grd.hgrid.lon_v[:]
            lat = grd.hgrid.lat_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos is 'w':
        # for w, AKt, ...
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:-1,:,:] + z[1:,:,:])
            z = np.concatenate((np.array(grd.vgrid.z_w[0,0,:,:], ndmin=3), \
                                z, \
                                np.array(grd.vgrid.z_w[0,-1,:,:], ndmin=3)), 0)
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-2:-1]), 2)
            lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
            lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        else:
            z = grd.vgrid.z_w[0,:]
            lon = grd.hgrid.lon_rho[:]
            lat = grd.hgrid.lat_rho[:]
        mask = grd.hgrid.mask_rho[:]

    elif Cpos is 'rho':
        # for temp, salt, rho, ...
        if vert == True: 
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-2:-1]), 2)
            lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
            lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        else:
            z = grd.vgrid.z_r[0,:]
            lon = grd.hgrid.lon_rho[:]
            lat = grd.hgrid.lat_rho[:]      
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning, '%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos

    # get constant-j slice
    varj = var[:,jindex,:]
    zj = z[:,jindex,:]
    lonj = np.tile(lon[jindex,:], (zj.shape[0], 1))
    latj = np.tile(lat[jindex,:], (zj.shape[0], 1))

    # land/sea mask
    maskj = np.tile(mask[jindex,:], (varj.shape[0], 1))
    varj = np.ma.masked_where(maskj[:,:] == 0, varj[:,:])

    return varj, zj, lonj, latj



def isoslice(var,prop,isoval, grd, Cpos='rho', masking=True, vert=False):
    """
    isoslice, lon, lat = isoslice(variable,property, isoval, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'	     specify the C-grid position where 
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
        raise ValueError, 'variable must have at least two dimensions'
    if not prop.shape == var.shape:
        raise ValueError, 'dimension of var and prop must be identical'

    # compute the depth on Arakawa-C grid position
    if Cpos is 'u':
        # average z_r at Arakawa-C u points
        z = 0.5 * (grd.vgrid.z_r[0,:,:,:-1] + grd.vgrid.z_r[0,:,:,1:])
        if vert == True: 
            lon = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
            lat = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
        else:
            lon = grd.hgrid.lon_u[:]
            lat = grd.hgrid.lat_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos is 'v':
        # average z_r at Arakawa-C v points
        z = 0.5 * (grd.vgrid.z_r[0,:,:-1,:] + grd.vgrid.z_r[0,:,1:,:])
        if vert == True: 
            lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
            lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        else:
            lon = grd.hgrid.lon_v[:]
            lat = grd.hgrid.lat_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos is 'rho':
        # for temp, salt, rho, w
        z = grd.vgrid.z_r[0,:]
        if vert == True: 
            lon = grd.hgrid.lon_vert[:]
            lat = grd.hgrid.lat_vert[:]
        else:
            lon = grd.hgrid.lon_rho[:]
            lat = grd.hgrid.lat_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning, '%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos
   
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
            raise Warning, 'property==%f out of range (%f, %f)' % \
                           (isoval, (prop+isoval).min(), (prop+isoval).max())
    isoslice = isoslice.reshape(sz[1:])

    # mask land
    isoslice = np.ma.masked_where(mask == 0, isoslice)

    return isoslice, lon, lat



def transect(var, istart, iend, jstart, jend, grd, Cpos='rho', vert=False, \
            spval=1e37):
    """
    transect, z, lon, lat = transect(var, istart, iend, jstart, jend, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'	     specify the C-grid position where 
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

    if Cpos is 'u':
        # average z_r and z_w at Arakawa-C u points
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = np.concatenate((z[:,0:1,:], z, z[:,-2:-1,:]), 1)
            lon = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
            lat = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            lon = grd.hgrid.lon_u[:]
            lat = grd.hgrid.lat_u[:]
        mask = grd.hgrid.mask_u[:]

    elif Cpos is 'v':
        # average z_r and z_w at Arakawa-C v points
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-2:-1]), 2)
            lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
            lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        else:
            z = grd.vgrid.z_r[0,:]
            z = 0.5 * (z[:,-1:,:] + z[:,1:,:])
            lon = grd.hgrid.lon_v[:]
            lat = grd.hgrid.lat_v[:]
        mask = grd.hgrid.mask_v[:]

    elif Cpos is 'rho':
        # for temp, salt, rho
        if vert == True:
            z = grd.vgrid.z_w[0,:]
            z = 0.5 * (z[:,:,:-1] + z[:,:,1:])
            z = 0.5 * (z[:,:-1,:] + z[:,1:,:])
            z = np.concatenate((z[:,:,0:1], z, z[:,:,-2:-1]), 2)
            z = np.concatenate((z[:,0:1,:], z, z[:,-2:-1,:]), 1)
            lon = grd.hgrid.lon_vert[:]
            lat = grd.hgrid.lat_vert[:]
        else:
            z = grd.vgrid.z_r[0,:]
            lon = grd.hgrid.lon_rho[:]
            lat = grd.hgrid.lat_rho[:]
        mask = grd.hgrid.mask_rho[:]

    else:
        raise Warning, '%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos


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
        print 'Here, the best line is y(x)'
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
        print 'Here, the best line is x(y)'
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
    lons = np.zeros((nlev, nearp.shape[0]))
    lats = np.zeros((nlev, nearp.shape[0]))

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
                lons[:,n] = lon[nearp[n,2], nearp[n,0]]
                lats[:,n] = lat[nearp[n,2], nearp[n,0]]
            else:
                zs[:,n] = (nearp[n,1] - nearp[n,2]) * z[:, nearp[n,3], nearp[n,0]] + \
                          (nearp[n,3] - nearp[n,1]) * z[:, nearp[n,2], nearp[n,0]]
                lons[:,n] = (nearp[n,1] - nearp[n,2]) * lon[nearp[n,3], nearp[n,0]] + \
                            (nearp[n,3] - nearp[n,1]) * lon[nearp[n,2], nearp[n,0]]
                lats[:,n] = (nearp[n,1] - nearp[n,2]) * lat[nearp[n,3], nearp[n,0]] + \
                            (nearp[n,3] - nearp[n,1]) * lat[nearp[n,2], nearp[n,0]]
        else:
            # check if our position match a grid cell
            if (nearp[n,2] == nearp[n,3]):
                zs[:,n] = z[:, nearp[n,0], nearp[n,2]]
                lons[:,n] = lon[nearp[n,0], nearp[n,2]]
                lats[:,n] = lat[nearp[n,0], nearp[n,2]]
            else:
                zs[:,n] = (nearp[n,1] - nearp[n,2]) * z[:, nearp[n,0], nearp[n,3]] + \
                          (nearp[n,3] - nearp[n,1]) * z[:, nearp[n,0], nearp[n,2]]
                lons[:,n] = (nearp[n,1] - nearp[n,2]) * lon[nearp[n,0], nearp[n,3]] + \
                            (nearp[n,3] - nearp[n,1]) * lon[nearp[n,0], nearp[n,2]]
                lats[:,n] = (nearp[n,1] - nearp[n,2]) * lat[nearp[n,0], nearp[n,3]] + \
                            (nearp[n,3] - nearp[n,1]) * lat[nearp[n,0], nearp[n,2]]

    # mask transect
    transect = np.ma.masked_values(transect, spval)


    return transect, zs, lons, lats
            

            

def lonslice(var, longitude, grd, Cpos='rho', vert=False, spval=1e37):
    """
    lonslice, z, lon, lat = lonslice(var, longitude, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'	     specify the C-grid position where 
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
        raise Warning, '%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos

    edge = np.concatenate((lon[1,1:-1], \
                           lon[1:-1,-2], \
                           lon[-2,-2:0:-1], \
                           lon[-2:0:-1,1]))
    idx =  np.concatenate((range(1,lon[0,:].shape[0]-1), \
                           range(1,lon[:,-1].shape[0]-1), \
                           range(1,lon[-1,::-1].shape[0]-1)[::-1], \
                           range(1,lon[::-1,0].shape[0]-1)[::-1]))

    d = np.zeros(edge.shape)
    for i in range (edge.shape[0]):
        d[i] = edge[i] - longitude
        d[i] = d[i] / abs(d[i])
    d = np.diff(d)

    pt_idx = np.where(d != 0)[0]

    Mp, Lp = lon.shape

    if len(pt_idx) != 2:
        raise ValueError, 'this function only works for simple quadrangle'

    # determine is latitude ligne is crossing a i or j edge
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
        lonslice, z, lon, lat = pyroms.tools.section(var, \
                                idx[pt_idx[0]], Lp-2, \
                                1, idx[pt_idx[1]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 1 and side[1] == 3:
        lonslice, z, lon, lat = pyroms.tools.section(var, \
                                idx[pt_idx[0]], idx[pt_idx[1]], \
                                1, Mp-2, \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 1 and side[1] == 4:
        lonslice, z, lon, lat = pyroms.tools.section(var, \
                                idx[pt_idx[0]], 1, \
                                1, idx[pt_idx[1]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 2 and side[1] == 3:
        lonslice, z, lon, lat = pyroms.tools.section(var, \
                                Lp-2, idx[pt_idx[1]], \
                                idx[pt_idx[0]], Mp-2, \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 2 and side[1] == 4:
        lonslice, z, lon, lat = pyroms.tools.section(var, \
                                Lp-2, 1, \
                                idx[pt_idx[0]], idx[pt_idx[0]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)
    elif side[0] == 3 and side[1] == 4:
        lonslice, z, lon, lat = pyroms.tools.section(var, \
                                idx[pt_idx[0]], 1, \
                                Mp-2, idx[pt_idx[1]], \
                                grd, Cpos=Cpos, vert=vert, spval=spval)


    return lonslice, z, lon, lat



def latslice(var, latitude, grd, Cpos='rho', vert=False, spval=1e37):
    """
    latslice, z, lon, lat = latslice(var, latitude, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'	     specify the C-grid position where 
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
        raise Warning, '%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos

    edge = np.concatenate((lat[1,1:-1], \
                           lat[1:-1,-2], \
                           lat[-2,-2:0:-1], \
                           lat[-2:0:-1,1]))
    idx =  np.concatenate((range(1,lat[0,:].shape[0]-1), \
                           range(1,lat[:,-1].shape[0]-1), \
                           range(1,lat[-1,::-1].shape[0]-1)[::-1], \
                           range(1,lat[::-1,0].shape[0]-1)[::-1]))

    d = np.zeros(edge.shape)
    for i in range (edge.shape[0]):
        d[i] = edge[i] - latitude
        d[i] = d[i] / abs(d[i])
    d = np.diff(d)

    pt_idx = np.where(d != 0)[0]

    Mp, Lp = lon.shape

    if len(pt_idx) != 2:
        raise ValueError, 'this function only works for simple quadrangle'

    # determine is latitude ligne is crossing a i or j edge
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



def section_transport(u, v, istart, iend, jstart, jend, grd):
    """
    transpu, transpv = section_transport(u, v, istart, iend, jstart, jend, grd)

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
        print 'Here, the best line is y(x)'
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
        print 'Here, the best line is x(y)'
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

    inear = copy(near)

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
    # average z_w at Arakawa-C u points  				           
    zu = 0.5 * (grd.vgrid.z_w[:,:,:-1] + grd.vgrid.z_w[:,:,1:])    			           
    dzu = zu[1:,:,:] - zu[:-1,:,:]
    # average z_w at Arakawa-C v points
    zv = 0.5 * (grd.vgrid.z_w[:,:-1,:] + grd.vgrid.z_w[:,1:,:])
    dzv = zv[1:,:,:] - zv[:-1,:,:]

    #set u and v to zero where u and v are masked for the sum
    for k in range(u.shape[0]):
        u[k,:] = np.where(grd.hgrid.mask_u == 1, u[k,:], 0)
        v[k,:] = np.where(grd.hgrid.mask_v == 1, v[k,:], 0)

    n = len(near)
    transpu = 0
    transpv = 0

    for l in range(0,n-1):
        ii = int(real(near[l])); jj = int(imag(near[l]))
        for k in range(0, dzu.shape[0]):
            if real(near[l]) == real(near[l+1]):
                trans = u[k, jj+jst, ii] * dy[jj+jst, ii] * \
                        dzu[k, jj+jst, ii] * norm_u * norm
                transpu = transpu + trans

            elif imag(near[l]) == imag(near[l+1]):
                trans = v[k, jj, ii+ist] * dx[jj, ii+ist] * \
                        dzv[k, jj, ii+ist] * norm_v * norm
                transpv = transpv + trans


    return transpu, transpv


def interm_pt(pnear, pk, pai, pbi, paj, pbj):
    ### FIND THE BEST INTERMEDIATE POINT ON A PATHWAY
    #		-----------------------------
    #	pnear	: vector of the position of the nearest point
    #	pk	: current working index
    #	pai, pbi: slope and original ordinate of x(y)
    #	paj, pbj: slope and original ordinate of y(x)
    #	pneari	: vector holding the position of intermediate point
    #		-----------------------------
    
    # 1 - Compute intermediate point
    
    # Determine whether we use y(x) or x(y):
    if (abs(paj) <= 1):
        # y(x)
	# possible intermediate point
	ylptmp1 = pnear[pk-1] + 1
	ylptmp2 = pnear[pk-1] + (paj/abs(paj))*1j
	# M is the candidate point:
	zxm = real(ylptmp1)
	zym = imag(ylptmp1)
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
	zxm = real(ylptmp2)
	zym = imag(ylptmp2)
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
	zxm = real(ylptmp1)
	zym = imag(ylptmp1)
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
	zxm = real(ylptmp2)
	zym = imag(ylptmp2)
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
