import numpy as np
import netCDF4 as netCDF
from datetime import datetime

import pyroms
import pyroms_toolbox


# load 2-dimentional interannual discharge data 
# from Hill and Beamer.
print('Load interannual discharge data')
nc_data = netCDF.Dataset('runoff.nc', 'r')
time = nc_data.variables['time'][:]
data = nc_data.variables['runoff'][:]

## time: cyclic year (365.25 days)
#time = np.array([15.21875, 45.65625, 76.09375, 106.53125, 136.96875, 167.40625, \
#    197.84375, 228.28125, 258.71875, 289.15625, 319.59375, 350.03125])


# load CI grid object
grd = pyroms.grid.get_ROMS_grid('COOK_INLET_LYON')


# define some variables
wts_file = 'remap_weights_runoff_to_CI_conservative_nomask.nc'
nt = data.shape[0]
Mp, Lp = grd.hgrid.mask_rho.shape
spval = -1e30
runoff_raw = np.zeros((Mp,Lp))
runoff = np.zeros((Mp,Lp))
rspread = 6

# create runoff file
#runoff_file = 'runoff_CI_daitren_inter_annual_2002-2004.nc'
runoff_file = 'CI_runoff.nc'
nc = netCDF.Dataset(runoff_file, 'w', format='NETCDF3_64BIT')
nc.Description = 'Hill & Beamer monthly climatology river discharge'
nc.Author = 'make_runoff_clim.py'
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'Hill & Beamer river discharge'

# create dimensions and variables
nc.createDimension('xi_rho', np.size(grd.hgrid.mask_rho,1))
nc.createDimension('eta_rho', np.size(grd.hgrid.mask_rho,0))
nc.createDimension('runoff_time', (365))

nc.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
nc.variables['lon_rho'].long_name = 'longitude of RHO-points'
nc.variables['lon_rho'].units = 'degree_east'
nc.variables['lon_rho'].field = 'lon_rho, scalar'
nc.variables['lon_rho'][:] = grd.hgrid.lon_rho

nc.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
nc.variables['lat_rho'].long_name = 'latitude of RHO-points'
nc.variables['lat_rho'].units = 'degree_north'
nc.variables['lat_rho'].field = 'lat_rho, scalar'
nc.variables['lat_rho'][:] = grd.hgrid.lat_rho

nc.createVariable('runoff_time', 'f8', ('runoff_time'))
nc.variables['runoff_time'].long_name = 'time'
nc.variables['runoff_time'].units = 'days since 1900-01-01 00:00:00'
nc.variables['runoff_time'].cycle_length = 365.25

nc.createVariable('Runoff_raw', 'f8', ('runoff_time', 'eta_rho', 'xi_rho'))
nc.variables['Runoff_raw'].long_name = 'Hill_Beamer River Runoff raw'
nc.variables['Runoff_raw'].missing_value = str(spval)
nc.variables['Runoff_raw'].units = 'kg/s/m^2'

nc.createVariable('Runoff', 'f8', ('runoff_time', 'eta_rho', 'xi_rho'))
nc.variables['Runoff'].long_name = 'Hill_Beamer River Runoff'
nc.variables['Runoff'].missing_value = str(spval)
nc.variables['Runoff'].units = 'kg/s/m^2'


# get littoral (here just 1 cell wide, no diagonals)
width = 1
idx = []
idy = []
maskl = grd.hgrid.mask_rho.copy()
for w in range(width):
    lit = pyroms_toolbox.get_littoral2(maskl)
    idx.extend(lit[0])
    idy.extend(lit[1])
    maskl[lit] = 0

littoral_idx = (np.array(idx), np.array(idy))
maskl = np.zeros(grd.hgrid.mask_rho.shape)
maskl[littoral_idx] = 1

mask_idx = np.where(grd.hgrid.mask_rho == 0)

nct=0
# Do January first
#for t in range(nt):
for t in range(nt-243,nt):
    flow = np.sum(data[t,280:600,160:460])
    print(nct+1, 'Remapping runoff for time %f' %time[t])
#    print 'Remapping runoff for time %f' %time[nct]
    # conservative horizontal interpolation using scrip
    runoff_raw = pyroms.remapping.remap(data[t,:,:], wts_file, \
                                           spval=spval)
    # Scale runoff to match incoming in Cook Inlet
    nflow = np.sum(runoff_raw)
    runoff_raw = runoff_raw*flow/nflow
    idx = np.where(runoff_raw != 0)
    runoff = pyroms_toolbox.move_runoff(runoff_raw, \
                  np.array(idx).T + 1, np.array(littoral_idx).T + 1, maskl, \
                  grd.hgrid.x_rho, grd.hgrid.y_rho, grd.hgrid.dx, grd.hgrid.dy)

    # write data in destination file
    nc.variables['Runoff'][nct] = runoff
    nc.variables['Runoff_raw'][nct] = runoff_raw
    nc.variables['runoff_time'][nct] = nct+1
# HACK    nc.variables['runoff_time'][nct] = time[nct]

    if t==180:
        print('Sum 2', np.sum(runoff_raw))
        print('Sum 3', np.sum(runoff))
    if nct==180:
        print('Sum 2 new 180', np.sum(runoff_raw))
        print('Sum 3 new 180', np.sum(runoff))
    nct = nct + 1

# Get rest of year
for t in range(nt-243):
    flow = np.sum(data[t,280:600,160:460])
    print(nct+1, 'Remapping runoff for time %f' %time[t])

    # conservative horizontal interpolation using scrip
    runoff_raw = pyroms.remapping.remap(data[t,:,:], wts_file, \
                                           spval=spval)
    # Scale runoff to match incoming in Cook Inlet
    nflow = np.sum(runoff_raw)
    runoff_raw = runoff_raw*flow/nflow
    idx = np.where(runoff_raw != 0)
    runoff = pyroms_toolbox.move_runoff(runoff_raw, \
                  np.array(idx).T + 1, np.array(littoral_idx).T + 1, maskl, \
                  grd.hgrid.x_rho, grd.hgrid.y_rho, grd.hgrid.dx, grd.hgrid.dy)
    # write data in destination file
    nc.variables['Runoff'][nct] = runoff
    nc.variables['Runoff_raw'][nct] = runoff_raw
    nc.variables['runoff_time'][nct] = nct+1
# HACK    nc.variables['runoff_time'][nct] = time[nct]
    nct = nct + 1

    if t==180:
        print('Sum 2', np.sum(runoff_raw))
        print('Sum 3', np.sum(runoff))

# close netcdf file
nc.close()
