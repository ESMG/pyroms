from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pyroms


def quiver(uvar, vvar, tindex, depth, gridid, \
           filename=None, proj=None, d=2, uscale=None, \
           xkey=0.9, ykey=0.1, ukey=1, outfile=None):
    """
    quiver(uvar, vvar, tindex, depth, gridid)

    optional switch:
      - filename         if defined, load the variable from file
      - proj             Basemap object returned by sview, zview, ...
      - d                arrow density parameter
      - uscale           data units per arrow length unit parameter
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
    zslicev = 0.5 * (zslicev[:-1,:] + zslicev[1:,:])       

    # correct dimension with zeros at edges
    zsliceu = concatenate((zeros((zsliceu.shape[0],1)), zsliceu, \
                           zeros((zsliceu.shape[0],1))),1)
    zsliceu = ma.masked_where(mask == 0, zsliceu)
    zsliceu = ma.masked_where(zsliceu >= 1000, zsliceu)
    zslicev = concatenate((zeros((1,zslicev.shape[1])), zslicev, \
                           zeros((1,zslicev.shape[1]))),0)
    zslicev = ma.masked_where(mask == 0, zslicev)
    zslicev = ma.masked_where(zslicev >= 1000, zslicev)

    # rotate velocity vector according to grid angle
    rotzsliceu = zsliceu * cos(grd.hgrid.angle_rho) - \
                 zslicev * sin(grd.hgrid.angle_rho)
    rotzslicev = zsliceu * sin(grd.hgrid.angle_rho) + \
                 zslicev * cos(grd.hgrid.angle_rho)

    # plot
    if proj is not None:
        x, y = proj(lon,lat)
    else:
        range = plt.axis()

    if uscale is None:
        if proj is not None: 
            qv = Basemap.quiver(proj, x[::d,::d], y[::d,::d], \
                                rotzsliceu[::d,::d], rotzslicev[::d,::d], \
                                linewidths=0.01)
        else: 
            qv = plt.quiver(lon[::d,::d], lat[::d,::d], \
                            rotzsliceu[::d,::d], rotzslicev[::d,::d], \
                            linewidths=0.01)
    else:
        if proj is not None:
            qv = Basemap.quiver(proj, x[::d,::d], y[::d,::d], \
                                      rotzsliceu[::d,::d], rotzslicev[::d,::d], \
                                     scale=uscale, linewidths=0.01)
        else:
            qv = plt.quiver(lon[::d,::d], lat[::d,::d], \
                            rotzsliceu[::d,::d], rotzslicev[::d,::d], \
                            scale=uscale, linewidths=0.01)

    if proj is None:
        plt.axis(range)

    plt.quiverkey(qv, xkey, ykey, ukey, str(ukey) + ' ms$^{-1}$')

    if outfile is not None:
        if outfile.find('.png') != -1 or outfile.find('.svg') != -1 or \
           outfile.find('.eps') != -1:
            print 'Write figure to file', outfile
            plt.savefig(outfile, dpi=100, facecolor='w', edgecolor='w', \
                        orientation='portrait')
        else:
            print 'Unrecognized file extension. Please use .png, .svg or .eps file extension.'

    return
