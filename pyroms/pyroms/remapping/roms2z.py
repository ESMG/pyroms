# encoding: utf-8

import numpy as np
import _interp

def roms2z(var, grd, grdz, Cpos='rho', irange=None, jrange=None, \
           spval=1e37, mode='linear'):
    """
    varz = roms2z(var, grd, grdz)

    optional switch:
      - Cpos='rho', 'u' 'v' or 'w'   specify the C-grid position where
                                     the variable rely
      - irange                       specify grid sub-sample for i direction
      - jrange                       specify grid sub-sample for j direction
      - spval=1e37                   define spval value
      - mode='linear' or 'spline'    specify the type of interpolation

    Interpolate the variable from ROMS grid grd to z vertical grid grdz
    """

    var = var.copy()

    assert len(var.shape) == 3, 'var must be 3D'

    if mode=='linear':
        imode=0
    elif mode=='spline':
        imode=1
    else:
        imode=0
        raise Warning('%s not supported, defaulting to linear' % mode)


    if Cpos is 'rho':
        z = grd.vgrid.z_r[0,:]
        depth = grdz.vgrid.z
        mask = grd.hgrid.mask_rho
    elif Cpos is 'u':
        z = 0.5 * (grd.vgrid.z_r[0,:,:,:-1] + grd.vgrid.z_r[0,:,:,1:])
        depth = 0.5 * (grdz.vgrid.z[:,:,:-1] + grdz.vgrid.z[:,:,1:])
        mask = grd.hgrid.mask_u
    elif Cpos is 'v':
        z = 0.5 * (grd.vgrid.z_r[0,:,:-1,:] + grd.vgrid.z_r[0,:,1:,:])
        depth = 0.5 * (grdz.vgrid.z[:,:-1,:] + grdz.vgrid.z[:,1:,:])
        mask = grd.hgrid.mask_v
    elif Cpos is 'w':
        z = grd.vgrid.z_w[0,:]
        depth = grdz.vgrid.z
        mask = grd.hgrid.mask_rho
    else:
        raise Warning('%s unknown position. Cpos must be rho, u, v or w.' % Cpos)

    Nm, Mm, Lm = var.shape
    nlev = grdz.vgrid.N

    var = np.concatenate((var, var[-2:-1,:,:]), 0)
    z = np.concatenate((z, 100*np.ones((1,z.shape[1], z.shape[2]))), 0)

    if irange is None:
        irange = (0,Lm)
    else:
        assert var.shape[2] == irange[1]-irange[0], \
               'var shape and irange must agree'

    if jrange is None:
        jrange = (0,Mm)
    else:
        assert var.shape[1] == jrange[1]-jrange[0], \
               'var shape and jrange must agree'

    varz = np.zeros((nlev, jrange[1]-jrange[0], irange[1]-irange[0]))

    for k in range(nlev):
        varz[k,:,:] = _interp.xhslice(var, \
                        z[:,jrange[0]:jrange[1], irange[0]:irange[1]], \
                        depth[k,jrange[0]:jrange[1], irange[0]:irange[1]], \
                        mask[jrange[0]:jrange[1], irange[0]:irange[1]], \
                        imode, spval)

    #mask
    idx = np.where(abs((varz-spval)/spval)<=1e-5)
    varz[idx] = spval
    #varz = np.ma.masked_values(varz, spval, rtol=1e-5)

    return varz
