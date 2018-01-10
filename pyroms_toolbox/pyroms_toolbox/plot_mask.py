import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pyroms
import pyroms_toolbox


def plot_mask(gridid, Cpos='rho', proj=None, **kwargs):


    # get grid
    if type(gridid).__name__ == 'ROMS_Grid':
        grd = gridid
    else:
        grd = pyroms.grid.get_ROMS_grid(gridid)

    Cpos = str(Cpos)
    print(Cpos)

    # get grid information
    if Cpos == 'rho':
        lon = grd.hgrid.lon_vert
        lat = grd.hgrid.lat_vert
        mask = grd.hgrid.mask_rho

    elif Cpos == 'u':
        lon = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
        lat = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
        mask = grd.hgrid.mask_u

    elif Cpos == 'v':
        lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
        lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        mask = grd.hgrid.mask_v

    else:
        raise Warning('Cpos must be rho, u or v')

    # defined color map
    land_color = kwargs.pop('land_color', (0.6, 1.0, 0.6))
    sea_color = kwargs.pop('sea_color', (0.6, 0.6, 1.0))

    cm = plt.matplotlib.colors.ListedColormap([land_color, sea_color],
                                             name='land/sea')



    if proj is None:
        plt.pcolor(lon, lat, mask, cmap=cm, vmin=0, vmax=1, \
                   edgecolor='k', **kwargs)
        pyroms_toolbox.plot_coast_line(grd)
    else:
        x, y = proj(lon, lat)
        Basemap.pcolor(proj, x, y, mask, cmap=cm, vmin=0, vmax=1, \
                       edgecolor='k', **kwargs)
        pyroms_toolbox.plot_coast_line(grd, proj=proj)

        lon_min = lon.min()
        lon_max = lon.max()
        lat_min = lat.min()
        lat_max = lat.max()

        proj.drawmeridians(np.arange(lon_min,lon_max,(lon_max-lon_min)/5.001), \
                          labels=[0,0,0,1], fmt='%.1f')
        proj.drawparallels(np.arange(lat_min,lat_max,(lat_max-lat_min)/5.001), \
                          labels=[1,0,0,0], fmt='%.1f')

