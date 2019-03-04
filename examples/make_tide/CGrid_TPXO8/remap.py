import numpy as np
import os
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import matplotlib.pyplot as plt
import time
import datetime as dt
from matplotlib.dates import date2num, num2date

import pyroms
import pyroms_toolbox
from pyroms import _remapping

class nctime(object):
    pass

def remap(src_varname, src_file, src_grd, dst_grd, dst_file, dmax=0, cdepth=0, kk=0, dst_dir='./'):

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'

    # create tide file
    print('\nCreating tide file', dst_file)
    if os.path.exists(dst_file) is True:
        os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)

    # open tide file
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')

    # load var
    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]

    # get missing value
    spval = src_var.missing_value

    if src_varname == 'h_re':
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_t_to_rho.nc'
        dst_varname = 'h_re'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'tidal constituent amplitude (real)'
        units = 'meter'
        field = 'free-surface, scalar, series'
    if src_varname == 'h_im':
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_t_to_rho.nc'
        dst_varname = 'h_im'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'tidal constituent amplitude (imaginary)'
        units = 'meter'
        field = 'free-surface, scalar, series'

    if src_varname == 'u_re':
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_u_to_rho.nc'
        dst_varname = 'h_re'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'tidal constituent x-velocity (real)'
        units = 'meter per second'
        field = 'x-velocity, scalar, series'
    if src_varname == 'u_im':
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_u_to_rho.nc'
        dst_varname = 'u_im'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'tidal constituent x-velocity (imaginary)'
        units = 'meter per second'
        field = 'x-velocity, scalar, series'

    if src_varname == 'v_re':
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_v_to_rho.nc'
        dst_varname = 'v_re'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'tidal constituent y-velocity (real)'
        units = 'meter per second'
        field = 'y-velocity, scalar, series'
    if src_varname == 'v_im':
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_v_to_rho.nc'
        dst_varname = 'v_im'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'tidal constituent y-velocity (imaginary)'
        units = 'meter per second'
        field = 'y-velocity, scalar, series'

    else:
        raise ValueError('Undefined src_varname')

    # create variable in file
    print('Creating variable', dst_varname)
    nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    # nc.createVariable(dst_varname, 'f8', dimensions)
    nc.variables[dst_varname].long_name = long_name
    nc.variables[dst_varname].units = units
    nc.variables[dst_varname].field = field

    # remapping
    print('remapping', dst_varname, 'from', src_grd.name, \
              'to', dst_grd.name)

    # horizontal interpolation using scrip weights
    print('horizontal interpolation using scrip weights')
    dst_var = pyroms.remapping.remap(src_var, wts_file, spval=spval)

    # write data in destination file
    print('write data in destination file\n')
    nc.variables[dst_varname][:] = dst_var

    # close file
    nc.close()
    cdf.close()

