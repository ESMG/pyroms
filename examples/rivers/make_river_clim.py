import numpy as np
import netCDF4 as netCDF
from datetime import datetime

import pyroms
import pyroms_toolbox


# load 2-dimentional discharge data
print('Load discharge data')
nc_data = netCDF.Dataset('CI_runoff.nc', 'r')
nc_rivers = netCDF.Dataset('Cook_Inlet_rivers.nc', 'a')
data = nc_data.variables['Runoff'][:]
time = nc_data.variables['runoff_time'][:]
sign = nc_rivers.variables['river_sign'][:]
xi = nc_rivers.variables['river_Xposition'][:]
eta = nc_rivers.variables['river_Eposition'][:]
dir = nc_rivers.variables['river_direction'][:]

# load NWGOA grid object
grd = pyroms.grid.get_ROMS_grid('COOK_INLET_LYON')


# define some variables
nt = data.shape[0]
Nr = sign.shape[0]
Mp, Lp = grd.hgrid.mask_rho.shape
runoff = np.zeros((Nr))
count = np.zeros(grd.hgrid.mask_rho.shape, dtype=np.int32)

# from a Python forum - create an array of lists
filler = np.frompyfunc(lambda x: list(), 1, 1)
rivers = np.empty((Mp, Lp), dtype=np.object)
filler(rivers, rivers)

for k in range(Nr):
    if (sign[k]==1):
        count[eta[k],xi[k]] += 1
        rivers[eta[k],xi[k]].append(k)
    elif (sign[k]==-1 and dir[k]==0):
        count[eta[k],xi[k]-1] += 1
        rivers[eta[k],xi[k]-1].append(k)
    elif (sign[k]==-1 and dir[k]==1):
        count[eta[k]-1,xi[k]] += 1
        rivers[eta[k]-1,xi[k]].append(k)


nct=0
for t in range(nt):
    print('Remapping runoff for time %f' %time[t])
    for j in range(Mp):
        for i in range(Lp):
            for n in range(count[j,i]):
                frac = 1.0/count[j,i]
                k = rivers[j,i][n]
                runoff[k] = frac*data[t,j,i]


    if t==180:
        sum180 = np.sum(runoff)
    runoff = runoff*sign/86400.

    # write data in destination file
    nc_rivers.variables['river_transport'][nct] = runoff
    nc_rivers.variables['river_time'][nct] = time[nct]
    nct = nct + 1


# close netcdf file
nc_rivers.close()

print('sum 4', sum180)
