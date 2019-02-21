import numpy as np
import numpy.ma as ma
import netCDF4
#import sys

#ncfile = sys.argv[1]
infile = "goa_discharge_narr_bc_1980.nc"
nc = netCDF4.Dataset(infile, 'r', format='NETCDF3_CLASSIC')
maskfile = "../lat_lon.nc"
mc = netCDF4.Dataset(maskfile, 'r', format='NETCDF3_64BIT')
outfile = "runoff.nc"
out = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')

slow = nc.variables['discharge_cmd'][:]
#hours = nc.variables['time'][:]
mask = mc.variables['coast_cells'][0,:,:]
#print hours.shape
#print hours[:]

ntime, numy, numx = slow.shape

out.createDimension('x', numx)
out.createDimension('y', numy)
out.createDimension('time', ntime)

out.createVariable('time', 'f8', ('time'))
out.variables['time'].units = 'days since 1900-01-01'
out.createVariable('runoff', 'f4', ('time', 'y', 'x'), fill_value=-9999.)

for it in np.arange(ntime):
    time = it
    runoff = slow[it,:,:]
    runoff = ma.masked_where((mask < -100.), runoff)
    out.variables['time'][it] = time
    out.variables['runoff'][it,:,:] = runoff
    if it == 180:
        raw_180 = runoff[280:600,160:460]

print('Sum 1', np.sum(raw_180))

nc.close()
mc.close()
out.close()
