# encoding: utf-8

import numpy as np
from .. import _interp
import pdb

def sta2z(var, grd, grdz, Cpos='rho', srange=None, \
           spval=1e37, mode='linear'):
    """
    varz = roms2z(var, grd, grdz)

    optional switch:
      - Cpos                         specify vertical grid type
      - srange                       specify grid sub-sample of stations
      - spval=1e37                   define spval value
      - mode='linear' or 'spline'    specify the type of interpolation

    Interpolate the variable from stations grid grd to z vertical grid grdz
    """

    var = var.copy()

    assert len(var.shape) == 2, 'var must be 2D'

    if mode=='linear':
        imode=0
    elif mode=='spline':
        imode=1
    else:
        imode=0
        raise Warning('%s not supported, defaulting to linear' % mode)

    if Cpos == 'rho':
        z = grd.vgrid.z_r[0,:]
        depth = grdz.vgrid.z
    elif Cpos == 'w':
        z = grd.vgrid.z_w[0,:]
        depth = grdz.vgrid.z
    else:
        raise Warning('%s unknown position. Cpos must be rho or w.' % Cpos)

    var = var.T
    Nm, Sm = var.shape
    nlev = grdz.vgrid.N

    # putting in a fake i dimension
    var = np.dstack(var).T
    z = np.dstack(z).T
    depth = np.dstack(depth).T
    mask = np.ones((Sm, 1))

    # copy surface value to high in the sky
    var = np.concatenate((var, var[-2:-1,:,:]), 0)
    z = np.concatenate((z, 100*np.ones((1,z.shape[1], z.shape[2]))), 0)
#    print 'nlev', nlev, 'var.shape', var.shape, srange
#    print 'z.shape', z.shape
#    print 'mask.shape', mask.shape
#    print 'depth.shape', depth.shape

    if srange is None:
        srange = (0,Sm)
    else:
        assert var.shape[1] == srange[1]-srange[0], \
               'var shape and srange must agree'

    varz = np.zeros((nlev, srange[1]-srange[0], 1))

    for k in range(nlev):
        varz[k,:] = _interp.xhslice(var, \
                        z[:, srange[0]:srange[1], :], \
                        depth[k, srange[0]:srange[1], :], \
                        mask[srange[0]:srange[1], :], \
                        imode, spval)

#    pdb.set_trace()
    #mask
    idx = np.where(abs((varz-spval)/spval)<=1e-5)
    varz[idx] = spval
    #varz = np.ma.masked_values(varz, spval, rtol=1e-5)

    return varz
