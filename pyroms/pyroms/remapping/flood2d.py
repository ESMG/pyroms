# encoding: utf-8

import numpy as np
from .. import _remapping

def flood2d(varz, grdz, Cpos='rho', irange=None, jrange=None, \
          spval=1e37, dmax=0, cdepth=0, kk=0):
    """
    var = flood(var, grdz)

    optional switch:
      - Cpos='rho', 'u' or 'v'       specify the C-grid position where
                                     the variable rely
      - irange                       specify grid sub-sample for i direction
      - jrange                       specify grid sub-sample for j direction
      - spval=1e37                   define spval value
      - dmax=0                       if dmax>0, maximum horizontal
                                     flooding distance
    Flood varz on gridz
    """

    varz = varz.copy()
    varz = np.array(varz)

    assert len(varz.shape) == 2, 'var must be 2D'

    # replace spval by nan
    idx = np.where(abs((varz-spval)/spval)<=1e-5)
    varz[idx] = np.nan

    if Cpos == 'rho':
        x = grdz.hgrid.lon_rho
        y = grdz.hgrid.lat_rho
        z = grdz.vgrid.z[:]
        mask = grdz.hgrid.mask_rho
    elif Cpos == 'u':
        x = grdz.hgrid.lon_u
        y = grdz.hgrid.lat_u
        z = 0.5 * (grdz.vgrid.z[:,:,:-1] + grdz.vgrid.z[:,:,1:])
        mask = grdz.hgrid.mask_u
    elif Cpos == 'v':
        x = grdz.hgrid.lon_v
        y = grdz.hgrid.lat_v
        z = 0.5 * (grdz.vgrid.z[:,:-1,:] + grdz.vgrid.z[:,1:,:])
        mask = grdz.hgrid.mask_v
    elif Cpos == 'w':
        x = grdz.hgrid.lon_rho
        y = grdz.hgrid.lat_rho
        z = grdz.vgrid.z[:]
        mask = grdz.hgrid.mask_rho
    else:
        raise Warning('%s bad position. Use depth at Arakawa-C rho points instead.' % Cpos)

    Mm, Lm = varz.shape

    if irange is None:
        irange = (0,Lm)
    else:
        assert varz.shape[1] == irange[1]-irange[0], \
               'var shape and irange must agree'

    if jrange is None:
        jrange = (0,Mm)
    else:
        assert varz.shape[0] == jrange[1]-jrange[0], \
               'var shape and jrange must agree'

    x = x[jrange[0]:jrange[1], irange[0]:irange[1]]
    y = y[jrange[0]:jrange[1], irange[0]:irange[1]]
    z = z[:,jrange[0]:jrange[1], irange[0]:irange[1]]
    mask = mask[jrange[0]:jrange[1], irange[0]:irange[1]]

    # Finding nearest values in horizontal
    # critical deph => no change if depth is less than specified value
    msk = mask.copy()
#    c1 = np.array(msk, dtype=bool)
    c2 = np.isnan(varz[:,:]) == 1
#    c = c1 & c2
#    idxnan = np.where(c == True)
    idxnan = np.where(c2 == True)
    idx = np.where(c2 == False)
    if list(idx[0]):
#       print "inside test", len(idx[0]), len(idxnan[0])
        wet = np.zeros((len(idx[0]),2))
        dry = np.zeros((len(idxnan[0]),2))
        wet[:,0] = idx[0]+1
        wet[:,1] = idx[1]+1
        dry[:,0] = idxnan[0]+1
        dry[:,1] = idxnan[1]+1

        varz = _remapping.flood(varz, wet, dry, x, y, dmax)

    return varz
