# encoding: utf-8

import numpy as np
from pyroms import _remapping_fast
import xarray as xr

import pyroms

import creep

def flood_fast(varz, grd, pos='t', irange=None, jrange=None, \
          spval=1.2676506e+30, dxy=5, cdepth=0, kk=0):
    """
    var = flood_fast(var, Bgrd)

    optional switch:
      - Bpos='t'                     specify the grid position where
                                     the variable rely
      - irange                       specify grid sub-sample for i direction
      - jrange                       specify grid sub-sample for j direction
      - spval=1e35                   define spval value
      - dmax=0                       if dmax>0, maximum horizontal
                                     flooding distance
      - cdepth=0                     critical depth for flooding
                                     if depth<cdepth => no flooding
      - kk

    Flood varz on Bgrd
    """

    varz = varz.copy()
    varz = np.array(varz)

    assert len(varz.shape) == 3, 'var must be 3D'

    # replace spval by nan
    idx = np.where(abs((varz-spval)/spval)<=1e-5)
    varz[idx] = np.nan

    x = grd.lon_t
    y = grd.lat_t
    h = grd.h
    if pos == 't':
        mask = grd.mask_t[0,:,:]

    nlev, Mm, Lm = varz.shape

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

    x = x[jrange[0]:jrange[1], irange[0]:irange[1]]
    y = y[jrange[0]:jrange[1], irange[0]:irange[1]]
    h = h[jrange[0]:jrange[1], irange[0]:irange[1]]
    mask = mask[jrange[0]:jrange[1], irange[0]:irange[1]]

    # Finding nearest values in horizontal
    # critical depth => no change if depth is less than specified value
    cdepth = abs(cdepth)
    if cdepth != 0:
        idx = np.where(h >= cdepth)
        msk = np.zeros(mask.shape)
        msk[idx] = 1
    else:
        msk = mask.copy()

    for k in range(nlev-1,-1,-1):
        c1 = np.array(msk, dtype=bool)
        c2 = np.isnan(varz[k,:,:]) == 1
        if kk == 0:
            c3 = np.ones(mask.shape).astype(bool)
        else:
            c3 = np.isnan(varz[max(k-kk,0),:,:]) == 0
        c = c1 & c2 & c3
        idxnan = np.where(c == True)
        idx = np.where(c2 == False)
        if list(idx[0]):
            wet_mask = np.zeros(mask.shape)
            wet_mask[idx] = 1
            dry = np.zeros((len(idxnan[0]),2))
            dry[:,0] = idxnan[0]+1
            dry[:,1] = idxnan[1]+1

            idx = np.where(np.isnan(varz[k]) == 1)
            varz[k][idx] = spval
#           varz[k,:] = _remapping_fast.flood(varz[k,:], dry, wet_mask, dxy)
            varz[k,:] = creep.cslf(varz[k,:],spval,-200.,200.)
            print(k, varz[k,:].min() , varz[k,:].max())

    # Debugging output
#   flooded = xr.DataArray(varz)
#   flooded.to_netcdf("flooded_before_bot.nc")

    # drop the deepest values down
    idx = np.where(np.isnan(varz) == 1)
    varz[idx] = spval
    bottom = pyroms.utility.get_bottom(varz[::-1,:,:], mask, spval=spval)
    bottom = (nlev-1) - bottom
    for i in range(Lm):
        for j in range(Mm):
            if mask[j,i] == 1:
                varz[int(bottom[j,i]):,j,i] = varz[int(bottom[j,i]),j,i]

    return varz
