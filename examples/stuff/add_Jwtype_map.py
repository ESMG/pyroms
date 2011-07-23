import numpy as np
import netCDF4
import sys

ncfile = sys.argv[1]
nc = netCDF4.Dataset(ncfile, 'a', format='NETCDF3_CLASSIC')

mask = nc.variables['mask_rho'][:]
h = nc.variables['h'][:]
spval = 1e30

#wtype = 1 (Jerlov water type I)
wtype = np.ones(h.shape)

#wtype = 2 (Jerlov water type IA)
idx = np.where(h <= 1500)
wtype[idx] = 2

#wtype = 3 (Jerlov water type IB)
idx = np.where(h <= 750)
wtype[idx] = 3

#wtype = 4 (Jerlov water type II)
idx = np.where(h <= 250)
wtype[idx] = 4

#wtype = 5 (Jerlov water type III)
idx = np.where(h <= 100)
wtype[idx] = 5

#mask
idx = np.where(mask == 0)
wtype[idx] = spval

nc.createVariable('wtype_grid', 'f8', ('eta_rho', 'xi_rho'), fill_value=spval)
nc.variables['wtype_grid'].long_name = 'spatially variable Jerlov water type.'
nc.variables['wtype_grid'].units = 'N/A'
nc.variables['wtype_grid'][:] = wtype

nc.close()
