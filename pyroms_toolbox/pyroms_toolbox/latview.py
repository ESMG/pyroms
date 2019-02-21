import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import cm, colors
from mpl_toolkits.basemap import Basemap
import pyroms


def latview(var, tindex, latitude, gridid, filename=None, \
            cmin=None, cmax=None, clev=None, fill=False, \
            contour=False, d=4, lonrange=None, hrange=None,\
            fts=None, title=None, map=False, \
            pal=None, clb=True, outfile=None):
    """
    latview(var, tindex, latitude, gridid, {optional switch})

    optional switch:
      - filename         if defined, load the variable from file
      - cmin             set color minimum limit
      - cmax             set color maximum limit
      - clev             set the number of color step
      - fill             use contourf instead of pcolor
      - contour          overlay contour (request fill=True)
      - d                contour density (default d=4)
      - lonrange         longitude range
      - hrange           h range
      - fts              set font size (default: 12)
      - title            add title to the plot
      - map              if True, draw a map showing islice location
      - pal              set color map (default: cm.jet)
      - clb              add colorbar (defaul: True)
      - outfile          if defined, write figure to file

    plot a constante-latitudinal slice of variable var. If filename
    is provided, var must be a string and the variable will be load
    from the file.
    grid can be a grid object or a gridid. In the later case, the grid
    object correponding to the provided gridid will be loaded.
    """

    # get grid
    if type(gridid).__name__ == 'ROMS_Grid':
        grd = gridid
    else:
        grd = pyroms.grid.get_ROMS_grid(gridid)

    # get variable
    if filename == None:
        var = var
    else:
        data = pyroms.io.Dataset(filename)

        var = data.variables[var]

    Np, Mp, Lp = grd.vgrid.z_r[0,:].shape

    if tindex is not -1:
        assert len(var.shape) == 4, 'var must be 4D (time plus space).'
        K, N, M, L = var.shape
    else:
        assert len(var.shape) == 3, 'var must be 3D (no time dependency).'
        N, M, L = var.shape

    # determine where on the C-grid these variable lies
    if N == Np and M == Mp and L == Lp:
        Cpos='rho'
        lon = grd.hgrid.lon_vert
        lat = grd.hgrid.lat_vert
        mask = grd.hgrid.mask_rho

    if N == Np and M == Mp and L == Lp-1:
        Cpos='u'
        lon = 0.5 * (grd.hgrid.lon_vert[:,:-1] + grd.hgrid.lon_vert[:,1:])
        lat = 0.5 * (grd.hgrid.lat_vert[:,:-1] + grd.hgrid.lat_vert[:,1:])
        mask = grd.hgrid.mask_u

    if N == Np and M == Mp-1 and L == Lp:
        Cpos='v'
        lon = 0.5 * (grd.hgrid.lon_vert[:-1,:] + grd.hgrid.lon_vert[1:,:])
        lat = 0.5 * (grd.hgrid.lat_vert[:-1,:] + grd.hgrid.lat_vert[1:,:])
        mask = grd.hgrid.mask_v

    # get constante-lat slice
    if tindex == -1:
        var = var[:,:,:]
    else:
        var = var[tindex,:,:,:]

    if fill == True:
        latslice, zs, lons, lats, = pyroms.tools.latslice(var, latitude, grd, Cpos)
    else:
        latslice, zs, lons, lats, = pyroms.tools.latslice(var, latitude, grd, Cpos, vert=True)


    # plot
    if cmin is None:
        cmin = latslice.min()
    else:
        cmin = float(cmin)

    if cmax is None:
        cmax = latslice.max()
    else:
        cmax = float(cmax)

    if clev is None:
        clev = 100.
    else:
        clev = float(clev)

    dc = (cmax - cmin)/clev ; vc = np.arange(cmin,cmax+dc,dc)

    if pal is None:
        pal = cm.jet
    else:
        pal = pal

    if fts is None:
        fts = 12
    else:
        fts = fts

    #pal.set_over('w', 1.0)
    #pal.set_under('w', 1.0)
    #pal.set_bad('w', 1.0)

    pal_norm = colors.BoundaryNorm(vc,ncolors=256, clip = False)

    # clear figure
    #plt.clf()

    if map is True:
        # set axes for the main plot in order to keep space for the map
        if fts < 12:
            ax=None
        else:
            ax = plt.axes([0.15, 0.08, 0.8, 0.65])
    else:
        if fts < 12:
            ax=None
        else:
            ax=plt.axes([0.15, 0.1, 0.8, 0.8])


    if fill is True:
        cf = plt.contourf(lons, zs, latslice, vc, cmap = pal, norm = pal_norm, axes=ax)
    else:
        cf = plt.pcolor(lons, zs, latslice, cmap = pal, norm = pal_norm, axes=ax)

    if clb is True:
        clb = plt.colorbar(cf, fraction=0.075,format='%.2f')
        for t in clb.ax.get_yticklabels():
            t.set_fontsize(fts)

    if contour is True:
        if fill is not True:
            raise Warning('Please run again with fill=True for overlay contour.')
        else:
            plt.contour(lons, zs, latslice, vc[::d], colors='k', linewidths=0.5, linestyles='solid', axes=ax)


    if lonrange is not None:
        plt.xlim(lonrange)

    if hrange is not None:
        plt.ylim(hrange)

    if title is not None:
        if map is True:
            # move the title on the right
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            xt = xmin - (xmax-xmin)/9.
            yt = ymax + (ymax-ymin)/7.
            plt.text(xt, yt, title, fontsize=fts+4)
        else:
            plt.title(title, fontsize=fts+4)

    plt.xlabel('Latitude', fontsize=fts)
    plt.ylabel('Depth', fontsize=fts)

    if map is True:
        # draw a map with constant-i slice location
        ax_map = plt.axes([0.4, 0.76, 0.2, 0.23])
        varm = np.ma.masked_where(mask[:,:] == 0, var[var.shape[0]-1,:,:])
        lon_min = lon.min()
        lon_max = lon.max()
        lon_0 = (lon_min + lon_max) / 2.
        lat_min = lat.min()
        lat_max = lat.max()
        lat_0 = (lat_min + lat_max) / 2.
        map = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min, \
                 urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0, \
                 resolution='i', area_thresh=10.)
        x, y = list(map(lon,lat))
        if lonrange is None:
            xs, ys = list(map(lons[0,:],lats[0,:]))
        else:
            c1 = lats[0,:] >= lonrange[0]
            c2 = lats[0,:] <= lonrange[1]
            c = c1 & c2
            idx = np.where(c == True)
            xs, ys = list(map(lons[0,idx[0]],lats[0,idx[0]]))
        # fill land and draw coastlines
        map.drawcoastlines()
        map.fillcontinents(color='grey')
        #map.drawmapboundary()
        Basemap.pcolor(map, x, y, varm, axes=ax_map)
        Basemap.plot(map, xs, ys, 'k-', linewidth=3, axes=ax_map)


    if outfile is not None:
        if outfile.find('.png') != -1 or outfile.find('.svg') != -1 or outfile.find('.eps') != -1:
            print('Write figure to file', outfile)
            plt.savefig(outfile, dpi=200, facecolor='w', edgecolor='w', orientation='portrait')
        else:
            print('Unrecognized file extension. Please use .png, .svg or .eps file extension.')


    return
