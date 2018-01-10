import scipy.io
import numpy as np
import netCDF4
import sys

#file = sys.argv[1]
file = 'grid_ll.mat'
mat = scipy.io.loadmat(file)

outfile = "lat_lon.nc"
out = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')

lat = mat['lat_dd']
lon = mat['lon_dd']
lon = lon + 360.
lat2 = np.zeros(lat.shape)
lon2 = np.zeros(lat.shape)

numy, numx = lat.shape
print(numy, numx)

for j in range(numy):
     lat2[j,:] = lat[numy-1-j,:]
     lon2[j,:] = lon[numy-1-j,:]

out.createDimension('x', numx)
out.createDimension('y', numy)

out.createVariable('lat', 'f4', ('y', 'x'))
out.createVariable('lon', 'f4', ('y', 'x'))
 
out.variables['lat'][:] = lat2
out.variables['lon'][:] = lon2

out.close()
