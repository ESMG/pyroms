#/usr/bin/env python

# A Python Implementation of Adding an Additional Variable to a NetCDF File,
# and adding values to that file.

import numpy as np
import netCDF4 as nc
import time

# Have the netCDF file listed below in same directory as this script
root = nc.Dataset('HYCOM_bdry_svalbard_2012.nc', 'a')

ubar = root.variables['ubar_north'][:]
u = root.variables['u_north'][:]

i=707
ubar[:,i] = 2/3.0*ubar[:,i-1]+1/3.0*ubar[:,i+2]
ubar[:,i+1] = 1/3.0*ubar[:,i-1]+2/3.0*ubar[:,i+2]
u[:,:,i] = 2/3.0*u[:,:,i-1]+1/3.0*u[:,:,i+2]
u[:,:,i+1] = 1/3.0*u[:,:,i-1]+2/3.0*u[:,:,i+2]

root.variables['ubar_north'][:] = ubar
root.variables['u_north'][:] = u

vbar = root.variables['vbar_north'][:]
v = root.variables['v_north'][:]

i=707
vbar[:,i] = 3/4.0*vbar[:,i-1]+1/4.0*vbar[:,i+3]
vbar[:,i+1] = 1/2.0*vbar[:,i-1]+1/2.0*vbar[:,i+3]
vbar[:,i+2] = 1/4.0*vbar[:,i-1]+3/4.0*vbar[:,i+3]
v[:,:,i] = 3/4.0*v[:,:,i-1]+1/4.0*v[:,:,i+3]
v[:,:,i+1] = 1/2.0*v[:,:,i-1]+1/2.0*v[:,:,i+3]
v[:,:,i+2] = 1/4.0*v[:,:,i-1]+3/4.0*v[:,:,i+3]

root.variables['vbar_north'][:] = vbar
root.variables['v_north'][:] = v

root.close()

