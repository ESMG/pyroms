import numpy as np
import os
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import matplotlib.pyplot as plt
import time
from datetime import datetime
from matplotlib.dates import date2num, num2date

import pyroms
import pyroms_toolbox
from pyroms import _remapping
import xarray as xr
import xesmf
from regrid_Hill import regrid_Hill

class nctime(object):
    pass

def remap_river(src_file, src_varname, cells_file, dst_grd, dst_dir='./', irange=None, jrange=None):

    cf = xr.open_dataset(cells_file)
    cell_mask = cf['coast_cells'][0,:,:]
    cell_mask = cell_mask.fillna(0)
    cell_mask = cell_mask.data

    dx = dst_grd.hgrid.dx
    dy = dst_grd.hgrid.dy
    area_dst = dx * dy
    mask_rho = dst_grd.hgrid.mask_rho

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    cdf = netCDF.Dataset(src_file)
    time = cdf.variables['time'][:]
    time = time + 25567
    print('time (days) =', time[0])

    # create runoff file
    dst_file = src_file.rsplit('/')[-1]
    dst_file = dst_dir + dst_file[:-3] + '_' + src_varname + '_Hill_' + dst_grd.name + '.nc'
    print('\nCreating file', dst_file)
    if os.path.exists(dst_file) is True:
        os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)

    # open River file
    nc = netCDF.Dataset(dst_file, 'a')

    #load var
    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]

    #get missing value
    spval = src_var._FillValue

    dst_varname = 'Runoff_raw'
    dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
    long_name = 'river discharge'
    units = 'meter^3/sec'

    # create variable in file
    print('Creating variable', dst_varname)
    nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    nc.variables[dst_varname].long_name = long_name
    nc.variables[dst_varname].units = 'm^3/day'

    nc.createVariable('Runoff', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'))
    nc.variables['Runoff'].long_name = 'Hill River Runoff'
    nc.variables['Runoff'].units = 'kg/m^2/sec'

    # get littoral (here just 1 cell wide, no diagonals)
    width = 1
    idx = []
    idy = []
    maskl = dst_grd.hgrid.mask_rho.copy()
    for w in range(width):
        lit = pyroms_toolbox.get_littoral2(maskl)
        idx.extend(lit[0])
        idy.extend(lit[1])
        maskl[lit] = 0

    littoral_idx = (np.array(idx), np.array(idy))
    maskl = np.zeros(mask_rho.shape)
    maskl[littoral_idx] = 1

#   flooded = xr.DataArray(maskl)
#   flooded.to_netcdf("maskl.nc")

    src_var_all = cdf.variables[src_varname][:]

    ntimes = len(time)
    for it in range(ntimes):
        src_var = src_var_all[it,:,:] * cell_mask
        net_flow = np.sum(src_var)
        src_var = xr.DataArray(src_var)
        src_var = src_var.fillna(0)
        src_var = src_var.data
#       if (it == 0):
#           coastal = xr.DataArray(src_var)
#           coastal.to_netcdf("coastal.nc")

        # determine variable dimension
        #ndim = len(src_var.dimensions)
        ndim = len(src_var.shape)
        if irange is not None and jrange is not None:
            if ndim == 3:
                src_var = src_var[:,jrange[0]:jrange[1],irange[0]:irange[1]]
            else:
                src_var = src_var[jrange[0]:jrange[1],irange[0]:irange[1]]

        Mp, Lp = mask_rho.shape

        # remapping
        print('remapping', dst_varname, 'from Hill rivers to NGOA')
        print('time =', time[it])

        # horizontal interpolation using xesmf
        print('horizontal interpolation using xesmf')

        runoff_raw = regrid_Hill(src_var, method='nearest_d2s', irange=irange, jrange=jrange)
        runoff_raw[0,:] = 0
        runoff_raw[:,-1] = 0
        flows = runoff_raw
        net_flow2 = np.sum(flows)
        print("flow before", net_flow)
        print("regridded flow", net_flow2)

        # write data in destination file
        nc.variables['ocean_time'][it] = time[it]
        nc.variables[dst_varname][it] = runoff_raw
#       if (it == 0):
#           coastal = xr.DataArray(runoff_raw)
#           coastal.to_netcdf("runoff_raw.nc")

        idx = np.where(runoff_raw[:,:] != 0)
#       breakpoint()
        runoff = pyroms_toolbox.move_river_t(runoff_raw, \
                  np.array(idx).T + 1, np.array(littoral_idx).T + 1, maskl, \
                  dst_grd.hgrid.lon_rho, dst_grd.hgrid.lat_rho, dx, dy)
        net_flow3 = np.sum(runoff)
        print("final flow", net_flow3)
        runoff = runoff * 1000.0 / area_dst  # get the units right
        nc.variables['Runoff'][it] = runoff

    # close destination file
    nc.close()

