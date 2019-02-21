import numpy as np
import netCDF4 as netCDF
from datetime import datetime

import pyroms
import pyroms_toolbox


# load 2-dimentional interannual discharge data 
# from 1948-2007. See Dai and Trenberth (2002) and Dai et al. (2009)
print('Load interannual discharge data')
nc_data = netCDF.Dataset('/archive/u1/uaf/kate/CORE2/runoff.daitren.clim.10FEB2011.nc', 'r')
data = nc_data.variables['runoff'][:]

# time: cyclic year (365.25 days)
time = np.array([15.21875, 45.65625, 76.09375, 106.53125, 136.96875, 167.40625, \
    197.84375, 228.28125, 258.71875, 289.15625, 319.59375, 350.03125])


# load ARCTIC2 grid object
grd = pyroms.grid.get_ROMS_grid('ARCTIC2')


# define some variables
wts_file = 'remap_weights_daitren_to_ARCTIC2_conservative_nomask.nc'
nt = data.shape[0]
Mp, Lp = grd.hgrid.mask_rho.shape
spval = -1e30
runoff_raw = np.zeros((Mp,Lp))
runoff = np.zeros((Mp,Lp))
rspread = 6

# create runoff file
#runoff_file = 'runoff_ARCTIC2_daitren_inter_annual_2002-2004.nc'
runoff_file = 'runoff_ARCTIC2_daitren_clim.nc'
nc = netCDF.Dataset(runoff_file, 'w', format='NETCDF3_64BIT')
nc.Description = 'Dai & Trenberth monthly climatology river discharge'
nc.Author = 'make_ARCTIC2_runoff.py'
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'Dai & Trenberth river discharge'

# creat dimensions and variables
nc.createDimension('xi_rho', np.size(grd.hgrid.mask_rho,1))
nc.createDimension('eta_rho', np.size(grd.hgrid.mask_rho,0))
nc.createDimension('runoff_time', (12))

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
nc.variables['runoff_time'][:] = time

nc.createVariable('Runoff_raw', 'f8', ('runoff_time', 'eta_rho', 'xi_rho'))
nc.variables['Runoff_raw'].long_name = 'Dai_Trenberth River Runoff raw'
nc.variables['Runoff_raw'].missing_value = str(spval)
nc.variables['Runoff_raw'].units = 'kg/s/m^2'

nc.createVariable('Runoff', 'f8', ('runoff_time', 'eta_rho', 'xi_rho'))
nc.variables['Runoff'].long_name = 'Dai_Trenberth River Runoff'
nc.variables['Runoff'].missing_value = str(spval)
nc.variables['Runoff'].units = 'kg/s/m^2'


# get littoral (here 4 cells width)
width = 4
idx = []
idy = []
maskl = grd.hgrid.mask_rho.copy()
for w in range(width):
    lit = pyroms_toolbox.get_littoral(maskl)
    idx.extend(lit[0])
    idy.extend(lit[1])
    maskl[lit] = 0

littoral_idx = (np.array(idx), np.array(idy))
maskl = np.zeros(grd.hgrid.mask_rho.shape)
maskl[littoral_idx] = 1

mask_idx = np.where(grd.hgrid.mask_rho == 0)

nct=0
for t in range(nt):
    print('Remapping runoff for time %f' %time[nct])
    # conservative horizontal interpolation using scrip
    runoff_raw = pyroms.remapping.remap(data[t,:,:], wts_file, \
                                           spval=spval)
    idx = np.where(runoff_raw != 0)
    runoff = pyroms_toolbox.move_runoff(runoff_raw, \
                  np.array(idx).T + 1, np.array(littoral_idx).T + 1, maskl, \
                  grd.hgrid.x_rho, grd.hgrid.y_rho, grd.hgrid.dx, grd.hgrid.dy)

    # spread the runoff within the littoral band
    runoff_spread = np.zeros((Mp,Lp))
    idx = np.where(runoff != 0)
    for p in range(np.size(idx,1)):
        j = list(range(max(0,idx[0][p]-rspread), min(Mp-1,idx[0][p]+rspread+1)))
        i = list(range(max(0,idx[1][p]-rspread), min(Lp-1,idx[1][p]+rspread+1)))
        ji = np.meshgrid(j,i)
        sidx = np.where(maskl[ji] == 1)
        nbpt = np.size(sidx) / 2
        rpt = runoff[idx[0][p],idx[1][p]] * grd.hgrid.dx[idx[0][p],idx[1][p]] * grd.hgrid.dy[idx[0][p],idx[1][p]]
        rpt = rpt / nbpt
        for pa in range(nbpt):
            pai = sidx[0][pa] + ji[1].min()
            paj = sidx[1][pa] + ji[0].min()
            runoff_spread[paj, pai] = runoff_spread[paj, pai] + \
                           rpt / (grd.hgrid.dx[paj, pai] * grd.hgrid.dy[paj, pai])


    # spval
    runoff_spread[mask_idx] = spval

    # write data in destination file
    nc.variables['Runoff'][nct] = runoff_spread
    nc.variables['Runoff_raw'][nct] = runoff_raw

    nct = nct + 1

# close netcdf file
nc.close()
