from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pyroms


def quiver(uvar, vvar, tindex, depth, gridid, \
           filename=None, proj=None, d=2, uscale=None, max_length=None, \
           xkey=0.9, ykey=0.1, ukey=1, outfile=None):
    """
    quiver(uvar, vvar, tindex, depth, gridid)

    optional switch:
      - filename         if defined, load the variable from file
      - proj             Basemap object returned by sview, zview, ...
      - d                arrow density parameter
      - uscale           data units per arrow length unit parameter
      - max_length       if defined, set maximum arrow length displayed
      - xkey             x location of the key
      - ykey             y location of the key
      - ukey             length of the key
      - outfile          if defined, write figure to file


    overlay a 2-D field of arrows for velocity (uvar, vvar) above an 
    existing horizontal 2D plot. If filename is provided, uvar and vvar 
    must be strings and the variables will be load from the file. 
    grid can be a grid object or a gridid. In the later case, the grid 
    object correponding to the provided gridid will be loaded. 
    For projection, use proj=map, map being the Basemap object returned 
    by sview, zview, ...

    Note: if quiver is called before any other part of the plot has been
    created, you must create an axis which covers the region to be plotted.
    to do this, you can call axis([Longitude_min,Longitude_max,Latitude_min,Latitude_max]
    where Longitude_min, etc, are replaced with the appropriate longitudes and latitudes.
    """


    # get grid
    if type(gridid).__name__ == 'ROMS_Grid':
        grd = gridid
    else:
        grd = pyroms.grid.get_ROMS_grid(gridid)
    lon = grd.hgrid.lon_rho
    lat = grd.hgrid.lat_rho
    mask = grd.hgrid.mask_rho    

    # get u and v
    if filename == None:

        if tindex is not -1:
            assert len(uvar.shape) == 4, 'uvar must be 4D (time plus space).'
            assert len(vvar.shape) == 4, 'vvar must be 4D (time plus space).'
        else:
            assert len(uvar.shape) == 3, 'uvar must be 3D (no time dependency).'
            assert len(vvar.shape) == 3, 'vvar must be 3D (no time dependency).'

        if tindex == -1:
            u = uvar[:,:,:]
            v = vvar[:,:,:]
        else:
            u = uvar[tindex,:,:,:]
            v = vvar[tindex,:,:,:]
    else:
        data = pyroms.io.Dataset(filename)
        u = data.variables[uvar][tindex,:,:,:]
        v = data.variables[vvar][tindex,:,:,:]

    # get u and v slice at requested depth
    zsliceu, lonu, latu = pyroms.tools.zslice(u, depth, grd, Cpos='u')
    zslicev, lonv, latv = pyroms.tools.zslice(v, depth, grd, Cpos='v')
 
    # average field at rho point position
    zsliceu = 0.5 * (zsliceu[:,:-1] + zsliceu[:,1:])    
    zsliceu = zsliceu[:,r_[0,:size(zsliceu,1),-1]]
    zsliceu = ma.masked_where(mask == 0, zsliceu)
    zsliceu = ma.masked_where(zsliceu >= 1000, zsliceu)
    zslicev = 0.5 * (zslicev[:-1,:] + zslicev[1:,:])       
    zslicev = zslicev[r_[0,:size(zslicev,0),-1],:]
    zslicev = ma.masked_where(mask == 0, zslicev)
    zslicev = ma.masked_where(zslicev >= 1000, zslicev)

    U = zsliceu + 1j * zslicev

    # rotate velocity vector according to grid angle
    U = U * exp(1j * grd.hgrid.angle_rho)

    #limit arrow length to max length if requested
    if max_length is not None:
        idx = where(abs(U) > max_length)
        U[idx] = U[idx] * max_length / abs(U[idx])

    # plot
    if proj is not None:
        x, y = proj(lon,lat)
    else:
        range = plt.axis()

    if uscale is None:
        if proj is not None: 
            qv = Basemap.quiver(proj, x[::d,::d], y[::d,::d], \
                                real(U[::d,::d]), imag(U[::d,::d]), \
                                linewidths=0.01)
        else: 
            qv = plt.quiver(lon[::d,::d], lat[::d,::d], \
                            real(U[::d,::d]), imag(U[::d,::d]), \
                            linewidths=0.01)
    else:
        if proj is not None:
            qv = Basemap.quiver(proj, x[::d,::d], y[::d,::d], \
                                real(U[::d,::d]), imag(U[::d,::d]), \
                                scale=uscale, linewidths=0.01)
        else:
            qv = plt.quiver(lon[::d,::d], lat[::d,::d], \
                            real(U[::d,::d]), imag(U[::d,::d]), \
                            scale=uscale, linewidths=0.01)

    if proj is None:
        plt.axis(range)

    plt.quiverkey(qv, xkey, ykey, ukey, str(ukey) + ' ms$^{-1}$')

    if outfile is not None:
        if outfile.find('.png') != -1 or outfile.find('.svg') != -1 or \
           outfile.find('.eps') != -1:
            print('Write figure to file', outfile)
            plt.savefig(outfile, dpi=100, facecolor='w', edgecolor='w', \
                        orientation='portrait')
        else:
            print('Unrecognized file extension. Please use .png, .svg or .eps file extension.')

    return
