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

def remap_bdry_bio_woa(argdict, src_grd, dst_grd, dmax=0, cdepth=0, kk=0, dst_dir='./'):

    # NWGOA3 grid sub-sample
    xrange=src_grd.xrange; yrange=src_grd.yrange
    src_varname = argdict['tracer']
    tracer = src_varname
    src_file = argdict['file']
    units = argdict['units']
    longname = argdict['longname']

    kt = argdict['frame']

    if src_varname == 'sio4':
       src_varname = 'si'

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'

    # create boundary file
    dst_file = tracer + '.nc'
    dst_file = dst_dir + dst_grd.name + '_bdry_bio_' + dst_file
    print('Creating boundary file', dst_file)
    if os.path.exists(dst_file) is True:
        os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_bdry_file(dst_file, dst_grd, nctime)

    # open boundary file
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')

    #load var
    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]

    # correct time to some classic value
    days_in_month = np.array([31,28.25,31,30,31,30,31,31,30,31,30,31])
    time = days_in_month[:kt].sum() + days_in_month[kt] / 2.

    #get missing value
    spval = src_var._FillValue
    spval2 = -1.0e+10

    # determine variable dimension
    ndim = len(src_var.dimensions) - 1

    # NWGOA3 grid sub-sample
    if ndim == 3:
        src_var = src_var[kt,:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    elif ndim == 2:
        src_var = src_var[kt,yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

    if tracer == 'no3':
       unit_conversion = 1. / 1e6 / 1.035
    elif tracer == 'po4':
       unit_conversion = 1. / 1e6 / 1.035
    elif tracer == 'o2':
       unit_conversion = 1. / 1035 / 22391.6 * 1000.0
    elif tracer == 'sio4':
       unit_conversion = 1. / 1e6 / 1.035

    src_var = src_var * unit_conversion

    Bpos = 't'
    Cpos = 'rho'
    z = src_grd.z_t
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    wts_file = 'remap_weights_ESM2M_to_NWGOA3_bilinear_t_to_rho.nc'
    dst_varname = tracer
    dst_varname_north = tracer + '_north'
    dimensions_north = ('ocean_time', 's_rho', 'xi_rho')
    long_name_north = longname + ' north boundary condition'
    field_north = tracer + '_north, scalar, series'
    dst_varname_south = tracer + '_south'
    dimensions_south = ('ocean_time', 's_rho', 'xi_rho')
    long_name_south = longname + ' south boundary condition'
    field_south = tracer + '_south, scalar, series'
    dst_varname_east = tracer + '_east'
    dimensions_east = ('ocean_time', 's_rho', 'eta_rho')
    long_name_east = longname + ' east boundary condition'
    field_east = tracer + '_east, scalar, series'
    dst_varname_west = tracer + '_west'
    dimensions_west = ('ocean_time', 's_rho', 'eta_rho')
    long_name_west = longname + ' west boundary condition'
    field_west = tracer + '_west, scalar, series'
    units = units


    if ndim == 3:
        # build intermediate zgrid
        zlevel = -z[::-1]
        nzlevel = len(zlevel)
        dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
        dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)


    # create variable in boudary file
    print('Creating variable', dst_varname_north)
    nc.createVariable(dst_varname_north, 'f8', dimensions_north, fill_value=spval2)
    nc.variables[dst_varname_north].long_name = long_name_north
    nc.variables[dst_varname_north].units = units
    nc.variables[dst_varname_north].field = field_north

    print('Creating variable', dst_varname_south)
    nc.createVariable(dst_varname_south, 'f8', dimensions_south, fill_value=spval2)
    nc.variables[dst_varname_south].long_name = long_name_south
    nc.variables[dst_varname_south].units = units
    nc.variables[dst_varname_south].field = field_south

    print('Creating variable', dst_varname_west)
    nc.createVariable(dst_varname_west, 'f8', dimensions_west, fill_value=spval2)
    nc.variables[dst_varname_west].long_name = long_name_west
    nc.variables[dst_varname_west].units = units
    nc.variables[dst_varname_west].field = field_west

    print('Creating variable', dst_varname_east)
    nc.createVariable(dst_varname_east, 'f8', dimensions_east, fill_value=spval2)
    nc.variables[dst_varname_east].long_name = long_name_east
    nc.variables[dst_varname_east].units = units
    nc.variables[dst_varname_east].field = field_east


    # remapping
    print('remapping', dst_varname, 'from', src_grd.name, \
              'to', dst_grd.name)

    if ndim == 3:
        # flood the grid
        print('flood the grid')
        src_varz = pyroms_toolbox.BGrid_GFDL.flood(src_var, src_grd, Bpos=Bpos, spval=spval, \
                                dmax=dmax, cdepth=cdepth, kk=kk)
    else:
        src_varz = src_var

    # horizontal interpolation using scrip weights
    print('horizontal interpolation using scrip weights')
    dst_varz = pyroms.remapping.remap(src_varz, wts_file, spval=spval)

    if ndim == 3:
        # vertical interpolation from standard z level to sigma
        print('vertical interpolation from standard z level to sigma')
        dst_var_north = pyroms.remapping.z2roms(dst_varz[::-1, Mp-1:Mp, 0:Lp], \
                          dst_grdz, dst_grd, Cpos=Cpos, spval=spval, \
                          flood=False, irange=(0,Lp), jrange=(Mp-1,Mp))
        dst_var_south = pyroms.remapping.z2roms(dst_varz[::-1, 0:1, :], \
                          dst_grdz, dst_grd, Cpos=Cpos, spval=spval, \
                          flood=False, irange=(0,Lp), jrange=(0,1))
        dst_var_east = pyroms.remapping.z2roms(dst_varz[::-1, :, Lp-1:Lp], \
                          dst_grdz, dst_grd, Cpos=Cpos, spval=spval, \
                          flood=False, irange=(Lp-1,Lp), jrange=(0,Mp))
        dst_var_west = pyroms.remapping.z2roms(dst_varz[::-1, :, 0:1], \
                          dst_grdz, dst_grd, Cpos=Cpos, spval=spval, \
                          flood=False, irange=(0,1), jrange=(0,Mp))
    else:
        dst_var_north = dst_varz[-1, :]
        dst_var_south = dst_varz[0, :]
        dst_var_east = dst_varz[:, -1]
        dst_var_west = dst_varz[:, 0]

    dst_var_north[np.where(dst_var_north == spval)] = spval2
    dst_var_south[np.where(dst_var_south == spval)] = spval2
    dst_var_east[np.where(dst_var_east == spval)]   = spval2
    dst_var_west[np.where(dst_var_west == spval)]   = spval2

    # write data in destination file
    print('write data in destination file\n')
    nc.variables['ocean_time'][0] = time
    nc.variables['ocean_time'].cycle_length = 365.25
    nc.variables[dst_varname_north][0] = np.squeeze(dst_var_north)
    nc.variables[dst_varname_south][0] = np.squeeze(dst_var_south)
    nc.variables[dst_varname_east][0] = np.squeeze(dst_var_east)
    nc.variables[dst_varname_west][0] = np.squeeze(dst_var_west)

    # close file
    nc.close()
    cdf.close()

    if src_varname == 'eta':
        return dst_varz
