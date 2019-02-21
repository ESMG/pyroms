import re
import numpy as np
import netCDF4
import sys
import pdb

outfile = sys.argv[1]

# Read the river temperatures
f = open('KenaiRiverTemps.dat', 'r')
# Eat first two lines
f.readline()
f.readline()

# These are for the ROMS sources file
ttime = []
temp = []
salt = []

#pdb.set_trace()

for line in f:
    nul, a, b, c = re.split('\s+', line)
    ttime.append([a]) 
    temp.append(float(b)) 
    salt.append(0.0) 

print(temp)

# create file with all the objects
out = netCDF4.Dataset(outfile, 'a', format='NETCDF3_64BIT')

#out.createDimension('river_tracer_time', len(ttime))

#times = out.createVariable('river_tracer_time', 'f8', ('river_tracer_time'))
#times.units = 'day'
#times.cycle_length = 365.25
#times.long_name = 'river tracer time'

#temp = out.createVariable('river_temp', 'f8', ('river_tracer_time'))
#temp.long_name = 'river runoff potential temperature'
#temp.units = 'Celsius'
#temp.time = 'river_tracer_time'

#salt = out.createVariable('river_salt', 'f8', ('river_tracer_time'))
#salt.long_name = 'river runoff salinity'
#salt.time = 'river_tracer_time'

out.variables['river_tracer_time'][:] = ttime
out.variables['river_temp'][:] = temp
out.variables['river_salt'][:] = salt

out.close()
