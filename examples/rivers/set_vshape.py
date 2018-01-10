import re
import numpy as np
import netCDF4
import sys
import pdb

outfile = sys.argv[1]

# Set the vertical distribution of the river transport.

out = netCDF4.Dataset(outfile, 'a', format='NETCDF3_64BIT')
N = len(out.dimensions['s_rho'])
Nr = len(out.dimensions['river'])

vshape = np.zeros([N, Nr])
for k in range(N):
    vshape[k,:] = k

area = sum(vshape[:,0])
vshape = (1.0/area)*vshape
print(vshape[:,0])

vshape2 = np.zeros([N])
for k in range(N):
    vshape2[k] = 1

area = sum(vshape2[:])
vshape2 = (1.0/area)*vshape2
print(vshape2)

# Copper River
#for k in range(5490,5497):
#    vshape[:,k] = vshape2

# Susitna River
#for k in range(3279,3284):
#    vshape[:,k] = vshape2

# Others
#for k in range(3371,3375):
#    vshape[:,k] = vshape2
#vshape[:,3255] = vshape2

out.variables['river_Vshape'][:] = vshape

out.close()
