#encoding: utf-8

import sys
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import time
from datetime import datetime
#from matplotlib.nxutils import pnpoly
from scipy import interpolate

import pyroms
import _interp


def get_lonlat(iindex, jindex, grd, Cpos='rho'):
    """
    lon, lat = get_lonlat(iindex, jindex, grd)

    return the longitude (degree east) and latitude (degree north)
    for grid point (iindex, jindex)
    """

    if Cpos is 'u':
        lon = grd.hgrid.lon_u[:,:]
        lat = grd.hgrid.lat_u[:,:]
    elif Cpos is 'v':
        lon = grd.hgrid.lon_v[:,:]
        lat = grd.hgrid.lat_v[:,:]
    elif Cpos is 'rho':
        lon = grd.hgrid.lon_rho[:,:]
        lat = grd.hgrid.lat_rho[:,:]
    elif Cpos is 'psi':
        lon = grd.hgrid.lon_psi[:,:]
        lat = grd.hgrid.lat_psi[:,:]
    else:
        raise Warning('%s bad position. Cpos must be rho, psi, u or v.' % Cpos)

    return lon[jindex, iindex], lat[jindex, iindex]


def get_ij(longitude, latitude, grd, Cpos='rho'):
    """
    i, j = get_ij(longitude, latitude, grd)

    return the index of the closest point on the grid from the
    point (longitude,latitude) in degree
    """

    if Cpos is 'u':
        lon = grd.hgrid.lon_u[:,:]
        lat = grd.hgrid.lat_u[:,:]
    elif Cpos is 'v':
        lon = grd.hgrid.lon_v[:,:]
        lat = grd.hgrid.lat_v[:,:]
    elif Cpos is 'rho':
        lon = grd.hgrid.lon_rho[:,:]
        lat = grd.hgrid.lat_rho[:,:]
    elif Cpos is 'psi':
        lon = grd.hgrid.lon_psi[:,:]
        lat = grd.hgrid.lat_psi[:,:]
    else:
        raise Warning('%s bad position. Cpos must be rho, psi, u or v.' % Cpos)

    lon = lon[:,:] - longitude
    lat = lat[:,:] - latitude

    diff = (lon * lon) + (lat * lat)

    jindex, iindex = np.where(diff==diff.min())

    return iindex[0], jindex[0]



def find_nearestgridpoints(longitude, latitude, grd, Cpos='rho'):


    if type(grd).__name__ == 'ROMS_Grid':

        if Cpos is 'u':
            lon = grd.hgrid.lon_u[:,:]
            lat = grd.hgrid.lat_u[:,:]
        elif Cpos is 'v':
            lon = grd.hgrid.lon_v[:,:]
            lat = grd.hgrid.lat_v[:,:]
        elif Cpos is 'rho':
            lon = grd.hgrid.lon_rho[:,:]
            lat = grd.hgrid.lat_rho[:,:]
        elif Cpos is 'vert':
            lon = grd.hgrid.lon_vert[:,:]
            lat = grd.hgrid.lat_vert[:,:]
        else:
            raise Warning('%s bad position. Cpos must be rho, u or v.' % Cpos)


    if type(grd).__name__ == 'CGrid_geo':

        if Cpos is 'u':
            lon = grd.lon_u[:,:]
            lat = grd.lat_u[:,:]
        elif Cpos is 'v':
            lon = grd.lon_v[:,:]
            lat = grd.lat_v[:,:]
        elif Cpos is 'rho':
            lon = grd.lon_rho[:,:]
            lat = grd.lat_rho[:,:]
        elif Cpos is 'vert':
            lon = grd.lon_vert[:,:]
            lat = grd.lat_vert[:,:]
        else:
            raise Warning('%s bad position. Cpos must be rho, u or v.' % Cpos)


    dlon = lon[:,:] - longitude
    dlat = lat[:,:] - latitude

    diff = (dlon * dlon) + (dlat * dlat)

    jidx, iidx = np.where(diff==diff.min())

    iidx = iidx[0] # take element 1 in case min dist is not unique
    jidx = jidx[0]

    try:

        iindex = [iidx, iidx+1, iidx+1, iidx]
        jindex = [jidx, jidx, jidx+1, jidx+1]
        xp = lon[jindex, iindex]
        yp = lat[jindex, iindex]
        verts = []
        for n in range(4):
            verts.append([xp[n], yp[n]])
        #inside = pnpoly(longitude, latitude, verts)
        inside = mpl.path.Path(verts).contains_point([longitude, latitude])

        if inside == 0:
            iindex = [iidx, iidx+1, iidx+1, iidx]
            jindex = [jidx-1, jidx-1, jidx, jidx]
            xp = lon[jindex, iindex]
            yp = lat[jindex, iindex]
            verts = []
            for n in range(4):
                verts.append([xp[n], yp[n]])
            #inside = pnpoly(longitude, latitude, verts)
            inside = mpl.path.Path(verts).contains_point([longitude, latitude])

            if inside == 0:
                iindex = [iidx-1, iidx, iidx, iidx-1]
                jindex = [jidx-1, jidx-1, jidx, jidx]
                xp = lon[jindex, iindex]
                yp = lat[jindex, iindex]
                verts = []
                for n in range(4):
                    verts.append([xp[n], yp[n]])
                #inside = pnpoly(longitude, latitude, verts)
                inside = mpl.path.Path(verts).contains_point([longitude, latitude])

                if inside == 0:
                    iindex = [iidx-1, iidx, iidx, iidx-1]
                    jindex = [jidx, jidx, jidx+1, jidx+1]
                    xp = lon[jindex, iindex]
                    yp = lat[jindex, iindex]
                    verts = []
                    for n in range(4):
                        verts.append([xp[n], yp[n]])
                    #inside = pnpoly(longitude, latitude, verts)
                    inside = mpl.path.Path(verts).contains_point([longitude, latitude])

                    if inside == 0:
                        raise ValueError('well where is it then?')

        iindex = iindex[:2]
        jindex = jindex[1:3]

    except:
        #print 'point (%f, %f) is not in the grid' %(longitude, latitude)
        iindex = []
        jindex = []

    return iindex, jindex



def get_coast_from_map(map):

    coast = []

    kk=len(map.coastsegs)
    for k in range(kk):
        ll = len(map.coastsegs[k])
        for l in range(ll):
            c = list(map(map.coastsegs[k][l][0], map.coastsegs[k][l][1], inverse=True))
            coast.append(c)
        coast.append((np.nan, np.nan))

    return np.array(coast)


def ijcoast(coast, grd):

    i, j = pyroms.tools.hindices(coast[:,0], coast[:,1], grd, spval=np.nan)
    ijcoast = np.concatenate([[i],[j]], axis=0).T

    return ijcoast


def ijcoast_old(coast, grd):


    if type(grd).__name__ == 'ROMS_Grid':
        lon = grd.hgrid.lon_vert
        lat = grd.hgrid.lat_vert
    if type(grd).__name__ == 'CGrid_geo':
        lon = grd.lon_vert
        lat = grd.lat_vert

    ijcoast = []

    for k in range(coast.shape[0]):
        if np.isnan(coast[k,0]):
            ijcoast.append([np.nan, np.nan])
        else:
            iindex, jindex = find_nearestgridpoints( \
                            coast[k,0], coast[k,1], grd, Cpos='vert')
            if iindex:
                i, j = np.meshgrid(iindex, jindex)
                x = lon[j,i]
                y = lat[j,i]
                funct_i = interpolate.interp2d(x.flatten(), y.flatten(), i.flatten())
                funct_j = interpolate.interp2d(x.flatten(), y.flatten(), j.flatten())
                i_coast = funct_i(coast[k,0], coast[k,1])[0]
                j_coast = funct_j(coast[k,0], coast[k,1])[0]
                ijcoast.append([i_coast, j_coast])

    return np.array(ijcoast)



def get_grid_proj(grd, type='merc', resolution='h', **kwargs):
    """
    map = get_grid_proj(grd)

    optional arguments:
      - type           set projection type (default is merc)
      - resolution     set resolution parameter (default is high)

    return a Basemap object that can be use for plotting
    """

    lon_min = grd.hgrid.lon_vert.min()
    lon_max = grd.hgrid.lon_vert.max()
    lon_0 = (lon_min + lon_max) / 2.

    lat_min = grd.hgrid.lat_vert.min()
    lat_max = grd.hgrid.lat_vert.max()
    lat_0 = (lat_min + lat_max) / 2.

    x_min = grd.hgrid.x_vert.min()
    x_max = grd.hgrid.x_vert.max()
    width = x_max - x_min

    y_max = grd.hgrid.y_vert.max()
    y_min = grd.hgrid.y_vert.min()
    height = y_max - y_min

    lat_1 = lat_min
    lat_2 = lat_max

    if type == 'lcc' or type == 'stere':
        map = Basemap(projection=type, width=width, height=height, \
                      lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0, \
                      resolution=resolution, **kwargs)
    else:
        map = Basemap(projection=type, llcrnrlon=lon_min, llcrnrlat=lat_min, \
                      urcrnrlon=lon_max, urcrnrlat=lat_max, \
                      lat_0=lat_0, lon_0=lon_0, \
                      resolution=resolution, **kwargs)

    return map


def get_nc_var(varname, filename):
    """
    var = roms_nc_var(varname, filename)

    a simple wraper for netCDF4
    """

    data = pyroms.io.Dataset(filename)
    var = data.variables[varname]

    return var


def roms_varlist(option):
    """
    varlist = roms_varlist(option)

    Return ROMS varlist.
    """

    if option == 'physics':
        varlist = (['temp','salt','u','v','ubar','vbar','zeta'])
    elif option == 'physics2d':
        varlist = (['ubar','vbar','zeta'])
    elif option == 'physics3d':
        varlist = (['temp','salt','u','v'])
    elif option == 'mixing3d':
        varlist = (['AKv','AKt','AKs'])
    elif option == 's-param':
        varlist = (['theta_s','theta_b','Tcline','hc'])
    elif option == 's-coord':
        varlist = (['s_rho','s_w','Cs_r','Cs_w'])
    elif option == 'coord':
        varlist = (['lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v'])
    elif option == 'grid':
        varlist = (['h','f','pm','pn','angle','lon_rho','lat_rho', \
          'lon_u','lat_u','lon_v','lat_v','lon_psi','lat_psi', \
          'mask_rho','mask_u','mask_v','mask_psi'])
    elif option == 'hgrid':
        varlist = (['f','dx','dy','angle_rho','lon_rho','lat_rho', \
          'lon_u','lat_u','lon_v','lat_v','lon_psi','lat_psi', \
          'mask_rho','mask_u','mask_v','mask_psi'])
    elif option == 'vgrid':
        varlist = (['h','s_rho','s_w','Cs_r','Cs_w', \
          'theta_s','theta_b','Tcline','hc'])
    else:
        raise Warning('Unknow varlist id')

    return varlist


def get_bottom(varz, mask, spval=1e37):

    assert len(varz.shape) == 3, 'var must be 3D'

    N = varz.shape[0]
    jj = varz.shape[1]
    ii = varz.shape[2]

    bottom = spval * np.ones((jj, ii))

    bottom[:,:] = _interp.get_bottom(varz,mask,spval)

    return bottom


def get_surface(varz, mask, spval=1e37):

    assert len(varz.shape) == 3, 'var must be 3D'

    N = varz.shape[0]
    jj = varz.shape[1]
    ii = varz.shape[2]

    surface = spval * np.ones((jj, ii))

    surface[:,:] = _interp.get_surface(varz,mask,spval)

    return surface


def move2grid(varin, init_grid, final_grid):
    '''
    tempu = move2grid(temp, 'rho', 'u')

    Move var from init_grid to final_grid.
    '''

    ndim = len(varin.shape)

    if ndim == 2:

        if (init_grid == 'rho' and final_grid == 'u'):
            varout = 0.5 * (varin[:,1:] + varin[:,:-1])
        elif (init_grid == 'rho' and final_grid == 'v'):
            varout = 0.5 * (varin[1:,:] + varin[:-1,:])
        elif (init_grid == 'rho' and final_grid == 'psi'):
            varout = 0.25 * (varin[1:,1:] + varin[:-1,:-1] + \
                             varin[1:,:-1] + varin[:-1,1:])
        elif (init_grid == 'u' and final_grid == 'psi'):
            varout = 0.5 * (varin[1:,:] + varin[:-1,:])
        elif (init_grid == 'v' and final_grid == 'psi'):
            varout = 0.5 * (varin[:,1:] + varin[:,:-1])
        else:
            raise ValueError('Undefined combination for init_grid and final_grid')

    elif ndim == 3:

        if (init_grid == 'rho' and final_grid == 'u'):
            varout = 0.5 * (varin[:,:,1:] + varin[:,:,:-1])
        elif (init_grid == 'rho' and final_grid == 'v'):
            varout = 0.5 * (varin[:,1:,:] + varin[:,:-1,:])
        elif (init_grid == 'rho' and final_grid == 'psi'):
            varout = 0.25 * (varin[:,1:,1:] + varin[:,:-1,:-1] + \
                             varin[:,1:,:-1] + varin[:,:-1,1:])
        elif (init_grid == 'u' and final_grid == 'psi'):
            varout = 0.5 * (varin[:,1:,:] + varin[:,:-1,:])
        elif (init_grid == 'v' and final_grid == 'psi'):
            varout = 0.5 * (varin[:,:,1:] + varin[:,:,:-1])
        else:
            raise ValueError('Undefined combination for init_grid and final_grid')

    else:
        raise ValueError('varin must be 2D or 3D')

    return varout


def get_date_tag(roms_time, ref=(2006, 0o1, 0o1), format="%d %b %Y at %H:%M:%S"):
    '''
    tag = get_date_tag(roms_time)

    return date tag for roms_time (in second since initialisation).
    default reference time is January 1st 2006.
    '''

    ref = time.mktime(datetime(ref[0], ref[1], ref[2]).timetuple())
    timestamp = ref + roms_time
    tag = datetime.fromtimestamp(timestamp).strftime(format)

    return tag


def apply_mask_change(file, grd):
    '''
    Apply mask change saved by edit_mesh_mask in the mask_change.txt file
    '''

    mask_changes = open(file,'r')
    lines = mask_changes.readlines()
    mask_changes.close()

    for line in lines:
        s = line.split()
        i = int(s[0])
        j = int(s[1])
        mask = float(s[2])
        grd.hgrid.mask_rho[j,i] = mask
