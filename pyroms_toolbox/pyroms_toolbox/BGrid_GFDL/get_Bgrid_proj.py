from mpl_toolkits.basemap import Basemap

def get_Bgrid_proj(Bgrd, type='merc', resolution='h'):
    """
    map = get_Bgrid_proj(Bgrd)

    optional arguments:
      - type           set projection type (default is merc)
      - resolution     set resolution parameter (default is high)

    return a Basemap object that can be use for plotting
    """

    lon_min = Bgrd.lon_t_vert.min()
    lon_max = Bgrd.lon_t_vert.max()
    lon_0 = (lon_min + lon_max) / 2.

    lat_min = Bgrd.lat_t_vert.min()
    lat_max = Bgrd.lat_t_vert.max()
    lat_0 = (lat_min + lat_max) / 2.

    map = Basemap(projection=type, llcrnrlon=lon_min, llcrnrlat=lat_min, \
         urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0, \
         resolution=resolution)

    return map

