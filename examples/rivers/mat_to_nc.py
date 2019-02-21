import scipy.io
import numpy as np
import netCDF4
import sys

file = sys.argv[1]
mat = scipy.io.loadmat(file)

outfile = "coast_cells.nc"
out = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')

coast_cells = mat['coast_cells']
nine = -9999.0*np.ones(coast_cells.shape)
coast_cells = np.where(np.isnan(coast_cells),nine,coast_cells)

numy, numx = coast_cells.shape
print(numy, numx)
ntime = 1

out.createDimension('x', numx)
out.createDimension('y', numy)
out.createDimension('time', ntime)

out.createVariable('coast_cells', 'f4', ('time', 'y', 'x'), fill_value=-9999.)
 
out.variables['coast_cells'][0,:,:] = coast_cells

out.close()
