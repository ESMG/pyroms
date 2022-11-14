# encoding: utf-8

import numpy as np
#from pyroms import _remapping

import creep

def flood2d(varz, Cgrd, irange=None, jrange=None, \
          spval=-9.99e+33):
    """
    var = flood(var, Cgrd)

    optional switch:
      - irange                       specify grid sub-sample for i direction
      - jrange                       specify grid sub-sample for j direction
      - spval=1e35                   define spval value

    Flood varz on Cgrd
    """

    varz = varz.copy()
    varz = np.array(varz)

    assert len(varz.shape) == 2, 'var must be 2D'

    # replace spval by nan
    idx = np.where(abs((varz-spval)/spval)<=1e-5)
    varz[idx] = np.nan

    x = Cgrd.lon_t
    y = Cgrd.lat_t
    h = Cgrd.h
    mask = Cgrd.mask_t[0,:,:]

    Mm, Lm = varz.shape

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

    x = x[irange[0]:irange[1]]
    y = y[jrange[0]:jrange[1]]
    h = h[jrange[0]:jrange[1], irange[0]:irange[1]]
    mask = mask[jrange[0]:jrange[1], irange[0]:irange[1]]

    # Finding nearest values in horizontal
    # critical deph => no change if depth is less than specified value
    msk = mask.copy()

    # replace spval by nan
    c1 = np.array(msk, dtype=bool)
    idx = np.where(np.logical_not(c1))
    varz[idx] = spval

    c2 = np.isnan(varz[:,:]) == 1
    c3 = np.ones(mask.shape).astype(bool)
    c = c1 & c3
    idxnan = np.where(c == True)
    idx = np.where(c2 == False)
    if list(idx[0]):
        wet = np.zeros((len(idx[0]),2))
        dry = np.zeros((len(idxnan[0]),2))
        wet[:,0] = idx[0]+1
        wet[:,1] = idx[1]+1
        dry[:,0] = idxnan[0]+1
        dry[:,1] = idxnan[1]+1

        varz[:] = creep.cslf(varz[:],spval,-200.,200.)

    print("after flood", varz.shape, varz[-1,1])
    return varz
