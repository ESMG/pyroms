import numpy as np
from mpl_toolkits.basemap import pyproj
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import pyroms



class BGrid_GFDL(object):
    """
    Arakawa B-Grid for GFDL CM2.1
    """


    def __init__(self, lon_t, lat_t, lon_uv, lat_uv, \
                       mask_t, mask_uv, h, z_t, z_t_edges, \
                       z_uv, z_uv_edges, f):

        self.name = 'GFDL_CM2.1'
        self.lon_t = lon_t
        self.lat_t = lat_t
        self.lon_uv = lon_uv
        self.lat_uv = lat_uv
        self.mask_t = mask_t
        self.mask_uv = mask_uv
        self.h = h
        self.z_t = z_t
        self.z_t_edges = z_t_edges
        self.z_uv = z_uv
        self.z_uv_edges = z_uv_edges
        self.f = f

        self._calculate_t_vert()
        self._calculate_uv_vert()
        self._calculate_grid_angle()


    def _calculate_t_vert(self):
        Mm, Lm = self.lon_t.shape
        lon = np.zeros((Mm+1,Lm+1))
        lat = np.zeros((Mm+1,Lm+1))

        lon[1:, 1:] = self.lon_uv[:,:]
        lat[1:, 1:] = self.lat_uv[:,:]

        #South edge
        lon[0,0:-1] = self.lon_t[0,:] - ( self.lon_uv[0,:] - self.lon_t[0,:] )
        lon[0,-1] = self.lon_t[0,-1] - ( self.lon_uv[0,-2] - self.lon_t[0,-1] )
        lat[0,0:-1] = self.lat_t[0,:] - ( self.lat_uv[0,:] - self.lat_t[0,:] )
        lat[0,-1] = self.lat_t[0,-1] - ( self.lat_uv[0,-2] - self.lat_t[0,-1] )

        #West edge
        lon[0:-1,0] = self.lon_t[:,0] - ( self.lon_uv[:,0] - self.lon_t[:,0] )
        lon[-1,0] = self.lon_t[-1,0] - ( self.lon_uv[-2,0] - self.lon_t[-1,0] )
        lat[0:-1,0] = self.lat_t[:,0] - ( self.lat_uv[:,0] - self.lat_t[:,0] )
        lat[-1,0] = self.lat_t[-1,0] - ( self.lat_uv[-2,0] - self.lat_t[-1,0] )

        self.lon_t_vert = lon
        self.lat_t_vert = lat


    def _calculate_uv_vert(self):
        Mm, Lm = self.lon_uv.shape
        lon = np.zeros((Mm+1,Lm+1))
        lat = np.zeros((Mm+1,Lm+1))

        lon[:-1, :-1] = self.lon_t[:,:]
        lat[:-1, :-1] = self.lat_t[:,:]

        #North edge
        lon[-1,0:-2] = self.lon_uv[-1,:-1] - ( self.lon_t[-1,1:] - self.lon_uv[-1,:-1] )
        lon[-1,-2:] = self.lon_uv[-1,-2:] - ( self.lon_t[-1,-2:] - self.lon_uv[-1,-2:] )
        lat[-1,0:-2] = self.lat_uv[-1,:-1] - ( self.lat_t[-1,1:] - self.lat_uv[-1,:-1] )
        lat[-1,-2:] = self.lat_uv[-1,-2:] - ( self.lat_t[-1,-2:] - self.lat_uv[-1,-2:] )

        #East edge
        lon[0:-2,-1] = self.lon_uv[:-1,-1] - ( self.lon_t[1:,-1] - self.lon_uv[:-1,-1] )
        lon[-2,-1] = self.lon_uv[-2:-1,-1] - ( self.lon_t[-2:-1,-1] - self.lon_uv[-2:-1,-1] )
        lat[0:-2,-1] = self.lat_uv[:-1,-1] - ( self.lat_t[1:,-1] - self.lat_uv[:-1,-1] )
        lat[-2,-1] = self.lat_uv[-2:-1,-1] - ( self.lat_t[-2:-1,-1] - self.lat_uv[-2:-1,-1] )

        self.lon_uv_vert = lon
        self.lat_uv_vert = lat


    def _calculate_grid_angle(self):
        geod = pyproj.Geod(ellps='WGS84')
        az_forward, az_back, dx = geod.inv(self.lon_t[:,:-1], self.lat_t[:,:-1], \
                                           self.lon_t[:,1:], self.lat_t[:,1:])

        angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
        self.angle = (90 - angle) * np.pi/180.



def get_nc_BGrid_GFDL(grdfile):
    """
    Bgrd = get_nc_BGrid_GFDL(grdfile)

    Load B-Grid grid object for GFDL CM2.1 from netCDF grid file
    """

    nc = pyroms.io.Dataset(grdfile)

    lon_t = nc.variables['geolon_t'][:]
    lat_t = nc.variables['geolat_t'][:]
    lon_uv = nc.variables['geolon_c'][:]
    lat_uv = nc.variables['geolat_c'][:]

    h = nc.variables['ht'][:]

    f = nc.variables['coriolis_param'][:]

    kmt = nc.variables['kmt'][:]
    z_t = nc.variables['st_ocean'][:]
    z_t_edges = nc.variables['st_edges_ocean'][:]

    kmu = nc.variables['kmu'][:]
    z_uv = nc.variables['sw_ocean'][:]
    z_uv_edges = nc.variables['sw_edges_ocean'][:]

    # compute mask at t-point
    M_t, L_t = kmt.shape
    N_t = z_t.shape[0]
    mask_t = np.zeros((N_t, M_t, L_t))
    for j in range(M_t):
        for i in range(L_t):
            try:
                mask_t[0:kmt[j,i], j,i] = 1
            except:
                mask_t[:, j,i] = 0

    # compute mask at uv-point
    M_uv, L_uv = kmu.shape
    N_uv = z_uv.shape[0]
    mask_uv = np.zeros((N_uv, M_uv, L_uv))
    for j in range(M_uv):
        for i in range(L_uv):
            try:
                mask_uv[0:kmt[j,i], j,i] = 1
            except:
                mask_uv[:, j,i] = 0

    return BGrid_GFDL(lon_t, lat_t, lon_uv, lat_uv, \
                       mask_t, mask_uv, h, z_t, z_t_edges, \
                       z_uv, z_uv_edges, f)



def make_remap_BGrid_GFDL_file(Bgrd, Bpos='t'):

    #create remap file
    remap_filename = 'remap_grid_' + Bgrd.name + '_' + Bpos + '.nc'
    nc = netCDF.Dataset(remap_filename, 'w', format='NETCDF3_CLASSIC')
    nc.Description = 'remap grid file on' + Bpos + 'points'
    nc.Author = 'pyroms.remapping.make_remap_grid_file'
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = Bgrd.name

    if Bpos == 't':
        lon_corner = Bgrd.lon_t_vert
        lat_corner = Bgrd.lat_t_vert
        grid_center_lon = Bgrd.lon_t.flatten()
        grid_center_lat = Bgrd.lat_t.flatten()
        grid_imask = Bgrd.mask_t[0,:].flatten()
        Mp, Lp = Bgrd.lon_t.shape
    elif Bpos == 'uv':
        lon_corner = Bgrd.lon_uv_vert
        lat_corner = Bgrd.lat_uv_vert
        grid_center_lon = Bgrd.lon_uv.flatten()
        grid_center_lat = Bgrd.lat_uv.flatten()
        grid_imask = Bgrd.mask_uv[0,:].flatten()
        Mp, Lp = Bgrd.lon_uv.shape
    else:
        raise ValueError('Bpos must be t or uv')

    grid_size = Lp * Mp

    grid_corner_lon = np.zeros((grid_size, 4))
    grid_corner_lat = np.zeros((grid_size, 4))
    k = 0
    for j in range(Mp):
        for i in range(Lp):
            grid_corner_lon[k,0] = lon_corner[j,i]
            grid_corner_lat[k,0] = lat_corner[j,i]
            grid_corner_lon[k,1] = lon_corner[j,i+1]
            grid_corner_lat[k,1] = lat_corner[j,i+1]
            grid_corner_lon[k,2] = lon_corner[j+1,i+1]
            grid_corner_lat[k,2] = lat_corner[j+1,i+1]
            grid_corner_lon[k,3] = lon_corner[j+1,i]
            grid_corner_lat[k,3] = lat_corner[j+1,i]
            k = k + 1

    #Write netcdf file
    nc.createDimension('grid_size', grid_size)
    nc.createDimension('grid_corners', 4)
    nc.createDimension('grid_rank', 2)

    nc.createVariable('grid_dims', 'i4', ('grid_rank'))
    nc.variables['grid_dims'].long_name = 'grid size along x and y axis'
    nc.variables['grid_dims'].units = 'None'
    nc.variables['grid_dims'][:] = [(Lp, Mp)]

    nc.createVariable('grid_center_lon', 'f8', ('grid_size'))
    nc.variables['grid_center_lon'].long_name = 'longitude of cell center'
    nc.variables['grid_center_lon'].units = 'degrees'
    nc.variables['grid_center_lon'][:] = grid_center_lon

    nc.createVariable('grid_center_lat', 'f8', ('grid_size'))
    nc.variables['grid_center_lat'].long_name = 'latitude of cell center'
    nc.variables['grid_center_lat'].units = 'degrees'
    nc.variables['grid_center_lat'][:] = grid_center_lat

    nc.createVariable('grid_imask', 'i4', ('grid_size'))
    nc.variables['grid_imask'].long_name = 'mask'
    nc.variables['grid_imask'].units = 'None'
    nc.variables['grid_imask'][:] = grid_imask

    nc.createVariable('grid_corner_lon', 'f8', ('grid_size', 'grid_corners'))
    nc.variables['grid_corner_lon'].long_name = 'longitude of cell corner'
    nc.variables['grid_corner_lon'].units = 'degrees'
    nc.variables['grid_corner_lon'][:] = grid_corner_lon

    nc.createVariable('grid_corner_lat', 'f8', ('grid_size', 'grid_corners'))
    nc.variables['grid_corner_lat'].long_name = 'latitude of cell corner'
    nc.variables['grid_corner_lat'].units = 'degrees'
    nc.variables['grid_corner_lat'][:] = grid_corner_lat

    nc.close()

