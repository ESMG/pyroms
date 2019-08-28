# encoding: utf-8

import numpy as np
from pyroms import _remapping

import pyroms

def flood(varz, Cgrd, Cpos='t', irange=None, jrange=None, \
          spval=-9.99e+33, dmax=0, cdepth=0, kk=0):
    """
    var = flood(var, Cgrd)

    optional switch:
      - Cpos='t', 'u', 'v'           specify the grid position where
                                     the variable resides
      - irange                       specify grid sub-sample for i direction
      - jrange                       specify grid sub-sample for j direction
      - spval=1e35                   define spval value
      - dmax=0                       if dmax>0, maximum horizontal
                                     flooding distance
      - cdepth=0                     critical depth for flooding
                                     if depth<cdepth => no flooding
      - kk

    Flood varz on Cgrd
    """

    varz = varz.copy()
    varz = np.array(varz)

    assert len(varz.shape) == 3, 'var must be 3D'

    # replace spval by nan
    idx = np.where(abs((varz-spval)/spval)<=1e-5)
    varz[idx] = np.nan

    x = Cgrd.lon_t
    y = Cgrd.lat_t
    h = Cgrd.h
    if Cpos is 't':
        mask = Cgrd.mask_t[0,:,:]
    elif Cpos is 'u':
        mask = Cgrd.mask_u[0,:,:]
    elif Cpos is 'v':
        mask = Cgrd.mask_v[0,:,:]

    nlev, Mm, Lm = varz.shape

    if irange is None:
        irange = (0,Lm)
    else:
        assert varz.shape[2] == irange[1]-irange[0], \
               'var shape and irange must agreed'

    if jrange is None:
        jrange = (0,Mm)
    else:
        assert varz.shape[1] == jrange[1]-jrange[0], \
               'var shape and jrange must agreed'

    x = x[jrange[0]:jrange[1], irange[0]:irange[1]]
    y = y[jrange[0]:jrange[1], irange[0]:irange[1]]
    h = h[jrange[0]:jrange[1], irange[0]:irange[1]]
    mask = mask[jrange[0]:jrange[1], irange[0]:irange[1]]

    # Finding nearest values in horizontal
    # critical deph => no change if depth is less than specified value
    cdepth = abs(cdepth)
    if cdepth != 0:
        idx = np.where(h >= cdepth)
        msk = np.zeros(mask.shape)
        msk[idx] = 1
    else:
        msk = mask.copy()
    for k in range(nlev-1,0,-1):
        c1 = np.array(msk, dtype=bool)
        c2 = np.isnan(varz[k,:,:]) == 1
        if kk == 0:
            c3 = np.ones(mask.shape).astype(bool)
        else:
            c3 = np.isnan(varz[min(k-kk,0),:,:]) == 0
        c = c1 & c2 & c3
        idxnan = np.where(c == True)
        idx = np.where(c2 == False)
        if list(idx[0]):
            wet = np.zeros((len(idx[0]),2))
            dry = np.zeros((len(idxnan[0]),2))
            wet[:,0] = idx[0]+1
            wet[:,1] = idx[1]+1
            dry[:,0] = idxnan[0]+1
            dry[:,1] = idxnan[1]+1

            varz[k,:] = _remapping.flood(varz[k,:], wet, dry, x, y, dmax)

    # drop the deepest values down
    idx = np.where(np.isnan(varz) == 1)
    varz[idx] = spval
    bottom = pyroms.utility.get_bottom(varz[::-1,:,:], mask, spval=spval)
    bottom = (nlev-1) - bottom
    for i in range(Lm):
        for j in range(Mm):
            if mask[j,i] == 1:
                varz[bottom[j,i]:,j,i] = varz[bottom[j,i],j,i]

    return varz
