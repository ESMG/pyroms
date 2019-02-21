import pyroms
import netCDF4 as nc
import numpy as npy

class iron_sediment():

    def __init__(self,domain):
        self.grd = pyroms.grid.get_ROMS_grid(domain)
        self.fileout = 'iron_sediment_' + domain + '.nc'
        self.ny, self.nx = self.grd.hgrid.lon_rho.shape
        return None

    def __call__(self):
        if self.grd.name == 'NWGOA3':
            self.ana_iron_ccs1()
        else:
            print('domain not supported') ; pass
        self.write_nc_file(self)
        return None


    def ana_iron_ccs1(self):
        h = self.grd.vgrid.h
        iron_sed = h.copy()
        # set to constant value everywhere
        # jk moore et al., global BGC cycles, 2004
        # suggest a value of 5 micromol/m2/day for california shelf
        iron_sed[:,:] = 0e-6 / 86400. # mol/m2/s
#        iron_sed[:,:] = 5e-6 / 86400. # mol/m2/s
        # set to zero deeper than cutoff depth
        iron_sed[npy.where( h > 1100. ) ] = 0.
        # mask the array
        iron_sed = iron_sed * self.grd.hgrid.mask_rho
        #
        self.iron_flx = npy.zeros((12,self.ny,self.nx))
        print(self.iron_flx.shape)
        for kt in npy.arange(12):
            self.iron_flx[kt,:,:] = iron_sed[:,:]
        return None

    def write_nc_file(self,fileout):
        fid = nc.Dataset(self.fileout, 'w', format='NETCDF3_CLASSIC')
        fid.description = 'Iron flux from sediment file (raphael@esm.rutgers.edu)'
        # dimensions
        fid.createDimension('lat', self.ny)
        fid.createDimension('lon', self.nx)
        fid.createDimension('ironsed_time', None)
        # variables
        latitudes  = fid.createVariable('lat', 'f4', ('lat','lon',))
        longitudes = fid.createVariable('lon', 'f4', ('lat','lon',))
        times      = fid.createVariable('ironsed_time', 'f4', ('ironsed_time',))
        times.units = "days since 1900-01-01 00:00"
        times.cycle_length = 365.25

        variable   = fid.createVariable('ironsed', 'f8', ('ironsed_time','lat','lon',))
        variable.coordinates="lon lat ironsed"
        variable.missing_value = 1e+15
        variable.time ="ironsed_time"
        variable.units = "mol/m2/s"
        variable.long_name = "iron flux from sediments"

        # data
        latitudes[:,:]    = self.grd.hgrid.lat_rho
        longitudes[:,:]   = self.grd.hgrid.lon_rho
        times[:]        = [15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5]
        variable[:,:,:] = self.iron_flx
        # close
        fid.close()
        return None


#-------------------------------------------------------------------------------
# iron sediment for NWGOA3

sediment_NWGOA3 = iron_sediment('NWGOA3')
sediment_NWGOA3()

