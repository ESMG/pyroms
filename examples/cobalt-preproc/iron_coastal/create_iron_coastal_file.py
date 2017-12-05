import pyroms
import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt
from scipy.ndimage import morphology as morph
import scipy.interpolate as scipyint
import mpl_toolkits.basemap as bmap

class iron_coastal():

    def __init__(self,domain):
        self.grd = pyroms.grid.get_ROMS_grid(domain)
        self.fileout = 'iron_fecoast_' + domain + '.nc'
        self.ny, self.nx = self.grd.hgrid.lon_rho.shape
        self.spval = 1.0e+15
        return None

    def __call__(self,mom_grid_file):
        self.maskcoast = self.create_mask_coast_MOM(mom_grid_file)
        if self.grd.name == 'NWGOA3':
            self.ana_iron_ccs1()
        else:
            print('domain not supported') ; pass
        self.write_nc_file(self)
        return None

    def create_mask_coast_MOM(self,mom_grid_file):
        fidmask = nc.Dataset(mom_grid_file,'r')
        kmt = fidmask.variables['kmt'][:].squeeze()
        lon_mom = fidmask.variables['geolon_t'][:].squeeze()
        lat_mom = fidmask.variables['geolat_t'][:].squeeze()
        fidmask.close()

        if self.grd.name == 'NWGOA3':
            lon_mom = lon_mom + 360.

        mask = np.zeros(kmt.shape)
        mask[np.where(kmt.mask)] = 1
        #plt.figure() ; plt.pcolormesh(mask) ; plt.colorbar() 

        mask_reduced = morph.binary_erosion(mask, structure=None, iterations=1, mask=None, output=None, border_value=0, origin=0, brute_force=False)
        mask_extended = morph.binary_dilation(mask, structure=None, iterations=1, mask=None, output=None, border_value=0, origin=0, brute_force=False)
        mask_coast = mask_extended - mask_reduced
        #plt.figure() ; plt.pcolormesh(mask_coast) ; plt.colorbar()

        # interpolate mask to target grid
        nymom, nxmom = lon_mom.shape

        sizetot_ini  = nymom * nxmom
        position_ini = np.zeros((sizetot_ini ,2))

        a = lon_mom.flat 
        b = lat_mom.flat

        for k in np.arange(sizetot_ini):
                position_ini[k,0] = a[k]
                position_ini[k,1] = b[k]

        sizetot_final  = self.grd.hgrid.lon_rho.shape[0]*self.grd.hgrid.lon_rho.shape[1]
        position_final = np.zeros((sizetot_final ,2))

        c = self.grd.hgrid.lon_rho.flat
        d = self.grd.hgrid.lat_rho.flat

        for k in np.arange(sizetot_final):
                position_final[k,0] = c[k]
                position_final[k,1] = d[k]

        tmp = scipyint.griddata(position_ini,mask_coast.flat,position_final,method='nearest')
        mask_coast_interp = np.reshape(tmp,self.grd.hgrid.lon_rho.shape)

        #plt.figure() ; plt.pcolormesh(mask_coast_interp) ; plt.colorbar() ; plt.show()
        return mask_coast_interp
        
    def ana_iron_ccs1(self):
        h = self.grd.vgrid.h
        one_fe_coast = h.copy()
        # set to constant value everywhere
        one_fe_coast = 1.0e-11 # mol Fe m kg-1 s-1
        # mask the array
        one_fe_coast = one_fe_coast * self.maskcoast * self.grd.hgrid.mask_rho
        one_fe_coast[np.where( self.grd.hgrid.mask_rho == 0)] = self.spval
        #plt.figure() ; plt.pcolormesh(one_fe_coast) ; plt.colorbar() ; plt.show()
        #
        self.fe_coast = np.zeros((12,self.ny,self.nx))
        for kt in np.arange(12):
            self.fe_coast[kt,:,:] = one_fe_coast[:,:]
        return None

    def write_nc_file(self,fileout):
        fid = nc.Dataset(self.fileout, 'w', format='NETCDF3_CLASSIC')
        fid.description = 'Iron coastal source file (raphael@esm.rutgers.edu)'
        # dimensions
        fid.createDimension('lat', self.ny)
        fid.createDimension('lon', self.nx)
        fid.createDimension('fecoast_time', None)
        # variables
        latitudes  = fid.createVariable('lat', 'f4', ('lat','lon',))
        longitudes = fid.createVariable('lon', 'f4', ('lat','lon',))
        times      = fid.createVariable('fecoast_time', 'f4', ('fecoast_time',))
        times.units = "days since 1900-01-01 00:00"
        times.cycle_length = 365.25

        variable   = fid.createVariable('fecoast', 'f8', ('fecoast_time','lat','lon',),fill_value=self.spval)
        variable.coordinates="lon lat fecoast_time"
        variable.missing_value = 1e+15
        variable.time ="fecoast_time"
        variable.units = "mol Fe m kg-1 s-1"
        variable.long_name = "iron coastal source"

        # data
        latitudes[:,:]    = self.grd.hgrid.lat_rho
        longitudes[:,:]   = self.grd.hgrid.lon_rho
        times[:]        = [15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5]
        variable[:,:,:] = self.fe_coast
        # close
        fid.close()
        return None

#-------------------------------------------------------------------------------
# iron coastal source for NWGOA3

coastal_NWGOA3 = iron_coastal('NWGOA3')
coastal_NWGOA3('/archive/u1/uaf/kate/COBALT/GFDL_CM2.1_grid.nc')



