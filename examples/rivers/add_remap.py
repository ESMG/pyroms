import re
import numpy as np
import netCDF4
import sys
import pdb
from datetime import datetime
import pyroms

def add_to_lists(pairs, i, j, sign, dir):
    x1, y1 = pairs[0]

    for it in range(1,len(pairs)):
        x2, y2 = pairs[it]

        if x2 > x1:
        # negative v-velocity
            i.append(x1)
            j.append(y1)
            sign.append(-1)
            dir.append(1)
        elif x1 > x2:
        # positive v-velocity
            i.append(x2)
            j.append(y1)
            sign.append(1)
            dir.append(1)
        elif y2 > y1:
        # positive u-velocity
            i.append(x1)
            j.append(y1)
            sign.append(1)
            dir.append(0)
        elif y1 > y2:
        # negative u-velocity
            i.append(x1)
            j.append(y2)
            sign.append(-1)
            dir.append(0)
        x1 = x2
        y1 = y2

#outfile = sys.argv[1]
outfile = 'remap_grid_CI_rivers.nc'

# We need to parse the output of the maskedge program for two
# different purposes:
#  1. Create the rivers file for ROMS, at least the locations part.
#  2. Create a scrip grid file for the river locations.
# This routine will only do #2.

# Read the landmask boundaries
f = open('maskedge.out', 'r')
pairs = []
# Eat first line so we don't trigger the add_to_lists routine
f.readline()

# These are for the ROMS sources file
i = []
j = []
sign = []
dir = []

#pdb.set_trace()

for line in f:
    a, b, c = re.split('\s+', line)
    if a=='-10':
        # wrap up object
        add_to_lists(pairs, i, j, sign, dir)
    elif (a=='-1' or a=='-3'):
        # wrap up object
        add_to_lists(pairs, i, j, sign, dir)
        # start new object
        pairs = []
    else:
        pairs.append([int(a),int(b)])

# set up grid coords
grd = pyroms.grid.get_ROMS_grid('COOK_INLET_LYON')
Mp, Lp = grd.vgrid.h.shape
# number of rivers
Nr = len(i)

grid_center_lon = np.zeros(Nr)
grid_center_lat = np.zeros(Nr)
grid_imask = np.ones(Nr)
grid_corner_lon = np.zeros((Nr, 4))
grid_corner_lat = np.zeros((Nr, 4))

for it in range(Nr):
    if (dir[it] == 0):
        grid_center_lon[it] = grd.hgrid.lon_u[j[it],i[it]-1]
        grid_center_lat[it] = grd.hgrid.lat_u[j[it],i[it]-1]
        grid_corner_lon[it,0] = grd.hgrid.lon_rho[j[it],i[it]-1]
        grid_corner_lat[it,0] = grd.hgrid.lat_rho[j[it],i[it]-1]
        grid_corner_lon[it,1] = grd.hgrid.lon_psi[j[it]-1,i[it]-1]
        grid_corner_lat[it,1] = grd.hgrid.lat_psi[j[it]-1,i[it]-1]
        grid_corner_lon[it,2] = grd.hgrid.lon_rho[j[it],i[it]]
        grid_corner_lat[it,2] = grd.hgrid.lat_rho[j[it],i[it]]
        grid_corner_lon[it,3] = grd.hgrid.lon_psi[j[it],i[it]-1]
        grid_corner_lat[it,3] = grd.hgrid.lat_psi[j[it],i[it]-1]
    else:
        grid_center_lon[it] = grd.hgrid.lon_v[j[it]-1,i[it]]
        grid_center_lat[it] = grd.hgrid.lat_v[j[it]-1,i[it]]
        grid_corner_lon[it,0] = grd.hgrid.lon_psi[j[it]-1,i[it]-1]
        grid_corner_lat[it,0] = grd.hgrid.lat_psi[j[it]-1,i[it]-1]
        grid_corner_lon[it,1] = grd.hgrid.lon_rho[j[it]-1,i[it]]
        grid_corner_lat[it,1] = grd.hgrid.lat_rho[j[it]-1,i[it]]
        grid_corner_lon[it,2] = grd.hgrid.lon_psi[j[it]-1,i[it]]
        grid_corner_lat[it,2] = grd.hgrid.lat_psi[j[it]-1,i[it]]
        grid_corner_lon[it,3] = grd.hgrid.lon_rho[j[it],i[it]]
        grid_corner_lat[it,3] = grd.hgrid.lat_rho[j[it],i[it]]


# create file with all the objects
nc = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')
nc.type = 'ROMS remap file'
nc.title = grd.name
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.source = 'David Hill and Jordan Beamer'

nc.createDimension('grid_size', len(i))
nc.createDimension('grid_corners', 4)
nc.createDimension('grid_rank', 1)

nc.createVariable('grid_dims', 'i4', ('grid_rank'))
nc.variables['grid_dims'].long_name = 'grid size along river axis'
nc.variables['grid_dims'].units = 'None'
nc.variables['grid_dims'][:] = [(len(i))]

nc.createVariable('grid_center_lon', 'f8', ('grid_size'))
nc.variables['grid_center_lon'].long_name = 'longitude of cell center'
nc.variables['grid_center_lon'].units = 'degrees'
nc.variables['grid_center_lon'][:] = grid_center_lon

nc.createVariable('grid_center_lat', 'f8', ('grid_size'))
nc.variables['grid_center_lat'].long_name = 'latitude of cell center'
nc.variables['grid_center_lat'].units = 'degrees'
nc.variables['grid_center_lat'][:] = grid_center_lat

nc.createVariable('grid_imask', 'i4', ('grid_size'))
nc.variables['grid_imask'].long_name = 'mask'
nc.variables['grid_imask'].units = 'None'
nc.variables['grid_imask'][:] = grid_imask

nc.createVariable('grid_corner_lon', 'f8', ('grid_size', 'grid_corners'))
nc.variables['grid_corner_lon'].long_name = 'longitude of cell corner'
nc.variables['grid_corner_lon'].units = 'degrees'
nc.variables['grid_corner_lon'][:] = grid_corner_lon

nc.createVariable('grid_corner_lat', 'f8', ('grid_size', 'grid_corners'))
nc.variables['grid_corner_lat'].long_name = 'latitude of cell corner'
nc.variables['grid_corner_lat'].units = 'degrees'
nc.variables['grid_corner_lat'][:] = grid_corner_lat

nc.close()

