import numpy as np
import netCDF4 as netCDF
from datetime import datetime

import pyroms
import pyroms_toolbox


# load 2-dimentional interannual discharge data 
# from 1948-2007. See Dai and Trenberth (2002) and Dai et al. (2009)
print('Load interannual discharge data')

runoff_file = 'NWGOA_runoff.nc'
nc = netCDF.Dataset(runoff_file, 'a', format='NETCDF3_64BIT')
runoff_raw = nc.variables['Runoff_raw'][:]
runoff = nc.variables['Runoff'][:]

raw_180 = runoff_raw[180,:,:]
runoff_180 = runoff[180,:,:]

print('Sum 1', np.sum(raw_180))
print('Sum 2', np.sum(runoff_180))

susitna = runoff[:,349,590]
copper = runoff[:,158,688]
copper2 = runoff[:,157,688]

print('copper', copper[180])
print('copper2', copper2[180])
print('susitna', susitna[180])

runoff[:,350,591] += susitna/3.0
runoff[:,348,590] += susitna/3.0
runoff[:,349,590] = susitna/3.0

runoff[:,157,688] += copper/3.0
runoff[:,159,686] += copper/3.0
runoff[:,158,688] = copper/3.0

runoff_180 = runoff[180,:,:]
print('Sum 3', np.sum(runoff_180))

print('copper', copper[180])
print('copper2', copper2[180])
print('susitna', susitna[180])

nc.variables['Runoff'][:] = runoff

# close netcdf file
nc.close()
