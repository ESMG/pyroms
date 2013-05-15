import numpy as np
import pyroms

def compute_moc(v, grd, mask_subregion=None):


    assert len(v.shape) == 3, 'v must be 3D'

    M, L = v[0,:].shape

    zw = np.array([0, 1.023907, 2.10319, 3.251309, 4.485053, 5.825238, 7.297443, 8.932686, \
                   10.7679, 12.84599, 15.21527, 17.92792, 21.03757, 24.59599, 28.64965, \
                   33.23697, 38.3871, 44.12101, 50.45447, 57.40257, 64.9846, 73.2287, \
                   82.17556, 91.88141, 102.4202, 113.8852, 126.3909, 140.074, 155.095, \
                   171.6402, 189.9228, 210.1845, 232.697, 257.7629, 285.7158, 316.9199, \
                   351.768, 390.6786, 434.0905, 482.4563, 536.2332, 595.8721, 661.8052, \
                   734.4321, 814.1057, 901.118, 995.6885, 1097.954, 1207.963, 1325.672, \
                   1450.95, 1583.582, 1723.28, 1869.693, 2022.425, 2181.044, 2345.101, \
                   2514.137, 2687.699, 2865.347, 3046.659, 3231.24, 3418.723, 3608.769, \
                   3801.072, 3995.354, 4191.367, 4388.89, 4587.726, 4787.702, 4988.667, \
                   5190.488, 5393.049, 5596.249, 5800, 6004.2285])

    zw = -zw[::-1]

    dz = np.diff(zw)
    dz = np.tile(dz, (L,M,1)).T

    zr = 0.5 * (zw[1:] + zw[:-1])

    nzr = len(zr)
    zcoord = pyroms.vgrid.z_coordinate(grd.vgrid.h, zr, nzr)
    grdz = pyroms.grid.ROMS_Grid(grd.name+'_Z', grd.hgrid, zcoord)

    vz = pyroms.remapping.roms2z(v, grd, grdz, Cpos='v', spval=0)

    #sub-region
    if mask_subregion is not None:
        if len(mask_subregion.shape) == 3:
            nb_subregion = mask_subregion.shape[0]
        elif len(mask_subregion.shape) == 2:
            mask_subregion = mask_subregion[np.newaxis, :]
            nb_subregion = mask_subregion.shape[0]
    else:
        mask_subregion = np.ones((1,M,L))
        nb_subregion = mask_subregion.shape[0]

    z = np.tile(zr, (nb_subregion,1))

    lat = np.zeros((nb_subregion,M))
    for r in range(nb_subregion):
        latv = grd.hgrid.lat_v * mask_subregion[r]
        lat[r] = latv.sum(axis=1) / mask_subregion[r].sum(axis=1)

    #mask_subregion = np.tile(mask_subregion, (1,nzr,1,1))
    mask_bassin = np.zeros((nb_subregion,nzr,M,L))
    for r in range(nb_subregion):
        mask_bassin[r,:] = np.tile(mask_bassin[r,:], (nzr,1,1))

    #dy at V-pos
    dy = grd.hgrid.dy
    dy = 0.5 * (dy[1:,:] + dy[:-1,:])
    dy = np.tile(dy,(nzr,1,1))

    #initialize moc
    vzonavg = np.zeros((nb_subregion,nzr,M))
    moc = np.zeros((nb_subregion,nzr,M))

    for r in range(nb_subregion):
        vzr = -vz * dy * dz * mask_bassin[r,:]
        vzonavg[r] = vzr.sum(axis=2) #integrate zonally
        moc[r] = vzonavg[r].copy()
        for k in range(1,nzr):
            moc[r,k,:] = moc[r,k,:] + moc[r,k-1,:]

    moc = moc / 1e6  #Sverdrup

    return lat, z, moc


