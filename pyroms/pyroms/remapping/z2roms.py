# encoding: utf-8

import numpy as np
import _interp
import _remapping

import pyroms

def z2roms(varz, grdz, grd, Cpos='rho', irange=None, jrange=None, \
           spval=1e37, flood=True, dmax=0, cdepth=0, kk=0, \
           mode='linear'):
    """
    var = z2roms(var, grdz, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'       specify the C-grid position where
                                     the variable rely
      - irange                       specify grid sub-sample for i direction
      - jrange                       specify grid sub-sample for j direction
      - spval=1e37                   define spval value
      - dmax=0                       if dmax>0, maximum horizontal
                                     flooding distance
      - cdepth=0                     critical depth for flooding
                                     if depth<cdepth => no flooding
      - kk
      - mode='linear' or 'spline'    specify the type of interpolation

    Interpolate the variable from z vertical grid grdz to ROMS grid grd
    """

    varz = varz.copy()

    assert len(varz.shape) == 3, 'var must be 3D'

    if mode=='linear':
        imode=0
    elif mode=='spline':
        imode=1
    else:
        raise Warning('%s not supported, defaulting to linear' % mode)

    if Cpos is 'rho':
        z = grdz.vgrid.z[:]
        depth = grd.vgrid.z_r[0,:]
        mask = grd.hgrid.mask_rho
    elif Cpos is 'u':
        z = 0.5 * (grdz.vgrid.z[:,:,:-1] + grdz.vgrid.z[:,:,1:])
        depth = 0.5 * (grd.vgrid.z_r[0,:,:,:-1] + grd.vgrid.z_r[0,:,:,1:])
        mask = grd.hgrid.mask_u
    elif Cpos is 'v':
        z = 0.5 * (grdz.vgrid.z[:,:-1,:] + grdz.vgrid.z[:,1:,:])
        depth = 0.5 * (grd.vgrid.z_r[0,:,:-1,:] + grd.vgrid.z_r[0,:,1:,:])
        mask = grd.hgrid.mask_v
    elif Cpos is 'w':
        z = grdz.vgrid.z[:]
        depth = grd.vgrid.z_w[0,:]
        mask = grd.hgrid.mask_rho
    else:
        raise Warning('%s bad position. Use depth at Arakawa-C \
                             rho points instead.' % Cpos)

    nlev, Mm, Lm = varz.shape
    Nm = depth.shape[0]

    if irange is None:
        irange = (0,Lm)
    else:
        assert varz.shape[2] == irange[1]-irange[0], \
               'var shape and irange must agree'

    if jrange is None:
        jrange = (0,Mm)
    else:
        assert varz.shape[1] == jrange[1]-jrange[0], \
               'var shape and jrange must agree'

    # flood varz if requested
    if flood is True:
        varz = pyroms.remapping.flood(varz, grdz, Cpos=Cpos, \
                 irange=irange, jrange=jrange, spval=spval, \
                 dmax=dmax, cdepth=cdepth, kk=kk)

    varz = np.concatenate((varz[0:1,:,:], varz, varz[-1:,:,:]), 0)
    z = np.concatenate((-9999*np.ones((1,z.shape[1], z.shape[2])), \
           z, \
           100*np.ones((1,z.shape[1], z.shape[2]))), 0)

    var = np.ma.zeros((Nm, Mm, Lm))

    for k in range(Nm):
        var[k,:,:] = _interp.xhslice(varz, \
                      z[:,jrange[0]:jrange[1], irange[0]:irange[1]], \
                      depth[k,jrange[0]:jrange[1], irange[0]:irange[1]], \
                      mask[jrange[0]:jrange[1], irange[0]:irange[1]], \
                      imode, spval)
        #mask
        var = np.ma.masked_values(var, spval, rtol=1e-5)
        #var[k,:,:] = np.ma.masked_where(mask == 0, var[k,:,:])

    return var
