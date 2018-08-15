#/usr/bin/env python

# A Python Implementation of Adding an Additional Variable to a NetCDF File,
# and adding values to that file.

import numpy as np
import netCDF4 as nc
import time

# Have the netCDF file listed below in same directory as this script
root = nc.Dataset('HYCOM_GLBa0.08_2011_181_ic_SVALBARD.nc', 'a')

temp = root.variables['temp'][:]
for j in range(63,111):
   for i in range(933,959):
       if temp[0,0,j,i] < 100:
           temp[0,:,j,i] = 0.0

salt = root.variables['salt'][:]
for j in range(165,169):
   for i in range(925,931):
       if salt[0,0,j,i] < 100:
           salt[0,:,j,i] = salt[0,:,165,925]
for j in range(63,149):
   for i in range(942,994):
       if salt[0,0,j,i] < 100:
           salt[0,:,j,i] = salt[0,:,120,949]

root.variables['temp'][:] = temp
root.variables['salt'][:] = salt

ubar = root.variables['ubar'][:]
u = root.variables['u'][:]

u[0,:,189,946] = 0
ubar[0,189,946] = 0
u[0,:,190,947] = 0
ubar[0,190,947] = 0

for j in range(216,566):
   for i in range(706,898):
       if np.fabs(ubar[0,j,i]) > 10:
           if np.fabs(ubar[0,j,i+1]) > 10:
               ubar[0,j,i] = 2/3.0*ubar[0,j,i-1]+1/3.0*ubar[0,j,i+2]
               ubar[0,j,i+1] = 1/3.0*ubar[0,j,i-1]+2/3.0*ubar[0,j,i+2]
               u[0,:,j,i] = 2/3.0*u[0,:,j,i-1]+1/3.0*u[0,:,j,i+2]
               u[0,:,j,i+1] = 1/3.0*u[0,:,j,i-1]+2/3.0*u[0,:,j,i+2]
           else:
               ubar[0,j,i] = 0.5*(ubar[0,j,i-1]+ubar[0,j,i+1])
               u[0,:,j,i] = 0.5*(u[0,:,j,i-1]+u[0,:,j,i+1])
         

vbar = root.variables['vbar'][:]
v = root.variables['v'][:]

v[0,:,189,946] = 0
vbar[0,189,946] = 0
v[0,:,188,947] = 0
vbar[0,188,947] = 0
v[0,:,189,947] = 0
vbar[0,189,947] = 0
v[0,:,190,948] = 0
vbar[0,190,948] = 0

for j in range(216,565):
   for i in range(708,900):
       if np.fabs(vbar[0,j,i]) > 10:
           if np.fabs(vbar[0,j,i+1]) > 10:
               vbar[0,j,i] = 2/3.0*vbar[0,j,i-1]+1/3.0*vbar[0,j,i+2]
               vbar[0,j,i+1] = 1/3.0*vbar[0,j,i-1]+2/3.0*vbar[0,j,i+2]
               v[0,:,j,i] = 2/3.0*v[0,:,j,i-1]+1/3.0*v[0,:,j,i+2]
               v[0,:,j,i+1] = 1/3.0*v[0,:,j,i-1]+2/3.0*v[0,:,j,i+2]
           else:
               vbar[0,j,i] = 0.5*(vbar[0,j,i-1]+vbar[0,j,i+1])
               v[0,:,j,i] = 0.5*(v[0,:,j,i-1]+v[0,:,j,i+1])
     
root.variables['ubar'][:] = ubar
root.variables['u'][:] = u
root.variables['vbar'][:] = vbar
root.variables['v'][:] = v

root.close()

