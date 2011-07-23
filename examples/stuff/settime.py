import numpy as np
import netCDF4
import sys

#ncfile = sys.argv[1]
ncfile = 'SODA_2.1.6_19841230-19850104_ic_CHUKCHI.nc'
nc = netCDF4.Dataset(ncfile, 'a', format='NETCDF3_CLASSIC')

time = nc.variables['ocean_time'][:]
time = 31052.

nc.variables['ocean_time'][:] = time

nc.close()
