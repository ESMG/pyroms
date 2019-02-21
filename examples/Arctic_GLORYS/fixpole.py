import numpy as np
import netCDF4
import sys

ncfile = sys.argv[1]
nc = netCDF4.Dataset(ncfile, 'a', format='NETCDF3_64BIT')

t = nc.variables['temp'][0,:,759:762,277:281]
s = nc.variables['salt'][0,:,759:762,277:281]
h = nc.variables['hice'][0,759:762,277:281]
a = nc.variables['aice'][0,759:762,277:281]

print('temp shape', t.shape)
t1 = (t[:,0,1] + t[:,1,0] + t[:,2,2] + t[:,2,1] + t[:,0,2] + t[:,1,3])/6.
s1 = (s[:,0,1] + s[:,1,0] + s[:,2,2] + s[:,2,1] + s[:,0,2] + s[:,1,3])/6.
h1 = (h[0,1] + h[1,0] + h[2,2] + h[2,1] + h[0,2] + h[1,3])/6.
a1 = (a[0,1] + a[1,0] + a[2,2] + a[2,1] + a[0,2] + a[1,3])/6.
print('temp 0', t[:,1,1])
print('temp 0a', t[:,1,2])
print('temp 1', t[:,1,0])
print('temp 2', t1)
print('salt 0', s[:,1,1])
print('salt 1', s[:,1,0])
print('salt 2', s1)

print('hice 0', h[1,1])
print('hice 0a', h[1,2])
print('hice 1', h[1,0])
print('hice 2', h1)
print('aice 0', a[1,1])
print('aice 1', a[1,0])
print('aice 2', a1)

nc.variables['temp'][0,:,760,278] = t1
nc.variables['salt'][0,:,760,278] = s1
nc.variables['temp'][0,:,760,279] = t1
nc.variables['salt'][0,:,760,279] = s1

nc.variables['hice'][0,760,278] = h1
nc.variables['aice'][0,760,278] = a1
nc.variables['hice'][0,760,279] = h1
nc.variables['aice'][0,760,279] = a1

nc.close()
