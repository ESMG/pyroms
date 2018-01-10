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
import _remapping

class nctime(object):
    pass

def remap_clm(src_file, src_varname, src_grd, dst_grd, dmax=0, cdepth=0, kk=0, dst_dir='./'):
    
    ystart=240

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    # time reference "days since 1900-01-01 00:00:00"
    ref = datetime(1900, 1, 1, 0, 0, 0)
    ref = date2num(ref)
# For IC
#    tag = src_file.rsplit('/')[-1].rsplit('_')[-1].rsplit('-')[0]
#    year = int(tag[:4])
#    month = int(tag[4:6])
#    day = int(tag[6:])
#    time = datetime(year, month, day, 0, 0, 0)
# For CLM
    year = int(src_file.rsplit('/')[-1].rsplit('_')[-2])
    month = int(src_file.rsplit('/')[-1].rsplit('_')[-1][0:2])
    day = np.array([15,14,15,15,15,15,15,15,15,15,15,15])
    hour = np.array([12,0,12,0,12,0,12,12,0,12,0,12])
    if year%4 == 0:
        hour = np.array([12,12,12,0,12,0,12,12,0,12,0,12])
    time = datetime(year, month, day[month-1], hour[month-1], 0, 0)
    time = date2num(time)
    time = time - ref
#    time = time + 2.5 # 5-day average

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

    # determine variable dimension
    ndim = len(src_var.dimensions)

    # global grid
    if ndim == 3:
        src_var = src_var[:]
        src_var = src_var[:,np.r_[ystart:np.size(src_var,1),-1],:]
    elif ndim == 2:
        src_var = src_var[:]
        src_var = src_var[np.r_[ystart:np.size(src_var,0),-1],:]

    if src_varname == 'ssh':
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_SODA_2.1.6_to_ARCTIC2_bilinear_t_to_rho.nc'
        dst_varname = 'zeta'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'free-surface'
        units = 'meter'
        field = 'free-surface, scalar, series'
        vartime = 'ocean_time'
    elif src_varname == 'temp':
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_SODA_2.1.6_to_ARCTIC2_bilinear_t_to_rho.nc'
        dst_varname = 'temp'
        dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        long_name = 'potential temperature'
        units = 'Celsius'
        field = 'temperature, scalar, series'
        vartime = 'ocean_time'
    elif src_varname == 'salt':
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_SODA_2.1.6_to_ARCTIC2_bilinear_t_to_rho.nc'
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
        print('flood the grid')
        src_varz = pyroms_toolbox.BGrid_SODA.flood(src_var, src_grd, Bpos=Bpos, spval=spval, \
                                dmax=dmax, cdepth=cdepth, kk=kk)
    else:
        src_varz = src_var

    # horizontal interpolation using scrip weights
    print('horizontal interpolation using scrip weights')
    dst_varz = pyroms.remapping.remap(src_varz, wts_file, \
                                          spval=spval)

    if ndim == 3:
        # vertical interpolation from standard z level to sigma
        print('vertical interpolation from standard z level to sigma')
        dst_var = pyroms.remapping.z2roms(dst_varz[::-1,:,:], dst_grdz, \
                          dst_grd, Cpos=Cpos, spval=spval, flood=False)
    else:
        dst_var = dst_varz

    # write data in destination file
    print('write data in destination file')
    nc.variables['ocean_time'][0] = time
    nc.variables[dst_varname][0] = dst_var

    # close destination file
    nc.close()

    if src_varname == 'SSH':
        return dst_varz
