#import matplotlib
#matplotlib.use('Agg')

#import numpy as np
import netCDF4
from datetime import datetime
#import pyroms
#import pyroms_toolbox
#import sys



# get HYCOM Northeast Pacific data from 2007 to 2011

invarname = 'angle'
outvarname = 'angle'

#read grid and variable attributes from the first file
url='http://viz.clusters.umaine.edu:8080/thredds/dodsC/pacific/1991-2008'
fp = netCDF4.Dataset(url)
lon = fp.variables['lon_rho'][880:1040,790:970]
lat = fp.variables['lat_rho'][880:1040,790:970]
mask = fp.variables['mask_rho'][880:1040,790:970]
lonu = fp.variables['lon_u'][880:1040,790:970-1]
latu = fp.variables['lat_u'][880:1040,790:970-1]
masku = fp.variables['mask_u'][880:1040,790:970-1]
lonv = fp.variables['lon_v'][880:1040-1,790:970]
latv = fp.variables['lat_v'][880:1040-1,790:970]
maskv = fp.variables['mask_v'][880:1040-1,790:970]
lonpsi = fp.variables['lon_psi'][880:1040-1,790:970-1]
latpsi = fp.variables['lat_psi'][880:1040-1,790:970-1]

h = fp.variables['h'][880:1040,790:970]
spherical = fp.variables['spherical']
angle = fp.variables['angle'][880:1040,790:970]
angle = 0.
units = fp.variables[invarname].units
long_name = fp.variables[invarname].long_name
fp.close()


year = 1991
day = 1


#create netCDF file
outfile = 'Pacific_NWGOA_grid.nc'
nc = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'U Maine Pacific with Cosine'

#create dimensions
Mp, Lp = lon.shape
N = 30
nc.createDimension('xi_rho', Lp)
nc.createDimension('eta_rho', Mp)
nc.createDimension('xi_u', Lp-1)
nc.createDimension('eta_u', Mp)
nc.createDimension('xi_v', Lp)
nc.createDimension('eta_v', Mp-1)
nc.createDimension('xi_psi', Lp-1)
nc.createDimension('eta_psi', Mp-1)
nc.createDimension('s_rho', N)
nc.createDimension('ocean_time', None)

#create variables        
nc.createVariable('lon_rho', 'f', ('eta_rho', 'xi_rho'))
nc.variables['lon_rho'].long_name = 'longitude'
nc.variables['lon_rho'].units = 'degrees_east'
nc.variables['lon_rho'][:] = lon

nc.createVariable('lat_rho', 'f', ('eta_rho', 'xi_rho'))
nc.variables['lat_rho'].long_name = 'latitude'
nc.variables['lat_rho'].units = 'degrees_north'
nc.variables['lat_rho'][:] = lat

nc.createVariable('mask_rho', 'f', ('eta_rho', 'xi_rho'))
nc.variables['mask_rho'].long_name = 'land mask'
nc.variables['mask_rho'][:] = mask

nc.createVariable('lon_u', 'f', ('eta_u', 'xi_u'))
nc.variables['lon_u'].long_name = 'longitude'
nc.variables['lon_u'].units = 'degrees_east'
nc.variables['lon_u'][:] = lonu

nc.createVariable('lat_u', 'f', ('eta_u', 'xi_u'))
nc.variables['lat_u'].long_name = 'latitude'
nc.variables['lat_u'].units = 'degrees_north'
nc.variables['lat_u'][:] = latu

nc.createVariable('mask_u', 'f', ('eta_u', 'xi_u'))
nc.variables['mask_u'].long_name = 'land mask'
nc.variables['mask_u'][:] = masku

nc.createVariable('lon_v', 'f', ('eta_v', 'xi_v'))
nc.variables['lon_v'].long_name = 'longitude'
nc.variables['lon_v'].units = 'degrees_east'
nc.variables['lon_v'][:] = lonv

nc.createVariable('lat_v', 'f', ('eta_v', 'xi_v'))
nc.variables['lat_v'].long_name = 'latitude'
nc.variables['lat_v'].units = 'degrees_north'
nc.variables['lat_v'][:] = latv

nc.createVariable('mask_v', 'f', ('eta_v', 'xi_v'))
nc.variables['mask_v'].long_name = 'land mask'
nc.variables['mask_v'][:] = maskv

nc.createVariable('lon_psi', 'f', ('eta_psi', 'xi_psi'))
nc.variables['lon_psi'].long_name = 'longitude'
nc.variables['lon_psi'].units = 'degrees_east'
nc.variables['lon_psi'][:] = lonpsi

nc.createVariable('lat_psi', 'f', ('eta_psi', 'xi_psi'))
nc.variables['lat_psi'].long_name = 'latitude'
nc.variables['lat_psi'].units = 'degrees_north'
nc.variables['lat_psi'][:] = latpsi

nc.createVariable('h', 'f', ('eta_rho', 'xi_rho'))
nc.variables['h'].long_name = 'depth'
nc.variables['h'].units = 'meter'
nc.variables['h'][:] = h

nc.createVariable('spherical', 'd')
nc.variables['spherical'].long_name = 'Grid type logical switch'
nc.variables['spherical'].option_T = 'spherical'
nc.variables['spherical'][:] = 1

nc.createVariable(outvarname, 'f', ('eta_rho', 'xi_rho'))
nc.variables[outvarname].long_name = long_name
nc.variables[outvarname].units = units
nc.variables[outvarname][:] = angle
        
nc.close()

