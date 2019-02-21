import numpy as np

def get_coast_line(grd, Cpos='rho'):
    '''
    coast = get_coast_line(grd)

    return the coastline from the grid object grid 
    '''

    if Cpos is 'rho':
        lon = grd.hgrid.lon_vert
        lat = grd.hgrid.lat_vert
        mask = grd.hgrid.mask_rho
    elif Cpos is 'u':
        lon = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
        lat = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
        mask = grd.hgrid.mask_u
    elif Cpos is 'v':
        lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
        lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        mask = grd.hgrid.mask_v
    elif Cpos is 'psi':
        lon = grd.hgrid.lon_rho
        lat = grd.hgrid.lat_rho
        mask = grd.hgrid.mask_psi
    else:
        raise Warning('%s bad position. Valid Arakawa-C are \
                           rho, u or v.' % Cpos)

    jidx, iidx = np.where(mask == 0)

    coast = []

    for i in range(iidx.shape[0]):
        if jidx[i] != mask.shape[0]-1:
            if mask[jidx[i], iidx[i]] != mask[jidx[i]+1, iidx[i]]:
                lonc = ([lon[jidx[i]+1,iidx[i]], lon[jidx[i]+1,iidx[i]+1]])
                latc = ([lat[jidx[i]+1,iidx[i]], lat[jidx[i]+1,iidx[i]+1]])
                seg = list(zip(lonc,latc))
                coast.append(seg)

        if jidx[i] != 0:
            if mask[jidx[i], iidx[i]] != mask[jidx[i]-1, iidx[i]]:
                lonc = ([lon[jidx[i],iidx[i]], lon[jidx[i],iidx[i]+1]])
                latc = ([lat[jidx[i],iidx[i]], lat[jidx[i],iidx[i]+1]])
                seg = list(zip(lonc,latc))
                coast.append(seg)

        if iidx[i] != mask.shape[1]-1:
            if mask[jidx[i], iidx[i]] != mask[jidx[i], iidx[i]+1]:
                lonc = ([lon[jidx[i],iidx[i]+1], lon[jidx[i]+1,iidx[i]+1]])
                latc = ([lat[jidx[i],iidx[i]+1], lat[jidx[i]+1,iidx[i]+1]])
                seg = list(zip(lonc,latc))
                coast.append(seg)

        if iidx[i] != 0:
            if mask[jidx[i], iidx[i]] != mask[jidx[i], iidx[i]-1]:
                lonc = ([lon[jidx[i],iidx[i]], lon[jidx[i]+1,iidx[i]]])
                latc = ([lat[jidx[i],iidx[i]], lat[jidx[i]+1,iidx[i]]])
                seg = list(zip(lonc,latc))
                coast.append(seg)

    return np.array(coast)
