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
from regrid_GLORYS import regrid_GLORYS

class nctime(object):
    pass

def remap_clm(src_file, src_varname, src_grd, dst_grd, cdepth=0, kk=0, dst_dir='./', irange=None, jrange=None):

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname][0]
    time = cdf.variables['time'][0]
    time = time / 24 + 18262

    # create IC file
    dst_file = src_file.rsplit('/')[-1]
    dst_file = dst_dir + dst_file[:-3] + '_' + src_varname + '_clim_' + dst_grd.name + '.nc'
    print('\nCreating file', dst_file)
    if os.path.exists(dst_file) is True:
        os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)

    # open IC file
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')

    #load var
    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]

    #get missing value
    spval = src_var._FillValue
    src_var = src_var[0]

    # determine variable dimension
    #ndim = len(src_var.dimensions)
    ndim = len(src_var.shape)
    if irange is not None and jrange is not None:
        if ndim == 3:
            src_var = src_var[:,jrange[0]:jrange[1],irange[0]:irange[1]]
        else:
            src_var = src_var[jrange[0]:jrange[1],irange[0]:irange[1]]


    Cpos = 'rho'
    z = src_grd.z_t
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    if src_varname == 'zos':
        dst_varname = 'zeta'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'free-surface'
        units = 'meter'
        field = 'free-surface, scalar, series'
        vartime = 'ocean_time'
    elif src_varname == 'thetao':
        dst_varname = 'temp'
        dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        long_name = 'potential temperature'
        units = 'Celsius'
        field = 'temperature, scalar, series'
        vartime = 'ocean_time'
    elif src_varname == 'so':
        dst_varname = 'salt'
        dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        long_name = 'salinity'
        units = 'PSU'
        field = 'salinity, scalar, series'
        vartime = 'ocean_time'
    else:
        raise ValueError('Undefined src_varname')


    if ndim == 3:
        # build intermediate zgrid
        zlevel = -z[::-1,0,0]
        nzlevel = len(zlevel)
        dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
        dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)


    # create variable in file
    print('Creating variable', dst_varname)
    nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    nc.variables[dst_varname].long_name = long_name
    nc.variables[dst_varname].units = units
    nc.variables[dst_varname].field = field
    nc.variables[dst_varname].time = vartime


    # remapping
    print('remapping', dst_varname, 'from', src_grd.name, \
              'to', dst_grd.name)
    print('time =', time)

    if ndim == 3:
        # flood the grid
        print('flood the grid, spval = ', spval)
        src_varz = pyroms_toolbox.Grid_GLORYS.flood_fast(src_var, src_grd, spval=spval, \
                                cdepth=cdepth, kk=kk)
    else:
        src_varz = pyroms_toolbox.Grid_GLORYS.flood2d(src_var, src_grd, spval=spval)

    # horizontal interpolation using xesmf
    print('horizontal interpolation using xesmf')

    dst_varz = regrid_GLORYS(src_varz, method='nearest_s2d', irange=irange, jrange=jrange)

    if ndim == 3:
        # vertical interpolation from standard z level to sigma
        print('vertical interpolation from standard z level to sigma')
        dst_var = pyroms.remapping.z2roms(dst_varz[::-1,:,:], dst_grdz, \
                          dst_grd, Cpos=Cpos, spval=spval, flood=False)

    # land mask
        idxu = np.where(dst_grd.hgrid.mask_rho == 0)
        for n in range(dst_grd.vgrid.N):
            dst_var[n,idxu[0], idxu[1]] = spval
    else:
        dst_var = dst_varz

    # write data in destination file
    print('write data in destination file')
    nc.variables['ocean_time'][0] = time
    nc.variables[dst_varname][0] = dst_var

    # close destination file
    nc.close()

    if src_varname == 'zos':
        return dst_varz
