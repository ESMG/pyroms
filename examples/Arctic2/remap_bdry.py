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

def remap_bdry(src_file, src_varname, src_grd, dst_grd, dmax=0, cdepth=0, kk=0, dst_dir='./'):
    
    ystart=240

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    # time reference "days since 1900-01-01 00:00:00"
    ref = datetime(1900, 1, 1, 0, 0, 0)
    ref = date2num(ref)
    tag = src_file.rsplit('/')[-1].rsplit('_')[-1].rsplit('-')[0]
    year = int(tag[:4])
    month = int(tag[4:6])
    day = int(tag[6:])
    time = datetime(year, month, day, 0, 0, 0)
    time = date2num(time)
    time = time - ref
    time = time + 2.5 # 5-day average

    # create boundary file
    dst_file = src_file.rsplit('/')[-1]
    dst_file = dst_dir + dst_file[:-4] + '_' + src_varname + '_bdry_' + dst_grd.name + '.nc'
    print('\nCreating boundary file', dst_file)
    if os.path.exists(dst_file) is True:
        os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)

    # open boundary file
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
        dst_varname_north = 'zeta_north'
        dimensions_north = ('ocean_time', 'xi_rho')
        long_name_north = 'free-surface north boundary condition'
        field_north = 'zeta_north, scalar, series'
        dst_varname_south = 'zeta_south'
        dimensions_south = ('ocean_time', 'xi_rho')
        long_name_south = 'free-surface south boundary condition'
        field_south = 'zeta_south, scalar, series'
        dst_varname_east = 'zeta_east'
        dimensions_east = ('ocean_time', 'eta_rho')
        long_name_east = 'free-surface east boundary condition'
        field_east = 'zeta_east, scalar, series'
        dst_varname_west = 'zeta_west'
        dimensions_west = ('ocean_time', 'eta_rho')
        long_name_west = 'free-surface west boundary condition'
        field_west = 'zeta_west, scalar, series'
        units = 'meter'
    elif src_varname == 'temp':
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_SODA_2.1.6_to_ARCTIC2_bilinear_t_to_rho.nc'
        dst_varname = 'temperature'
        dst_varname_north = 'temp_north'
        dimensions_north = ('ocean_time', 's_rho', 'xi_rho')
        long_name_north = 'potential temperature north boundary condition'
        field_north = 'temp_north, scalar, series'
        dst_varname_south = 'temp_south'
        dimensions_south = ('ocean_time', 's_rho', 'xi_rho')
        long_name_south = 'potential temperature south boundary condition'
        field_south = 'temp_south, scalar, series'
        dst_varname_east = 'temp_east'
        dimensions_east = ('ocean_time', 's_rho', 'eta_rho')
        long_name_east = 'potential temperature east boundary condition'
        field_east = 'temp_east, scalar, series'
        dst_varname_west = 'temp_west'
        dimensions_west = ('ocean_time', 's_rho', 'eta_rho')
        long_name_west = 'potential temperature west boundary condition'
        field_west = 'temp_west, scalar, series'
        units = 'Celsius'
    elif src_varname == 'salt':
        Bpos = 't'
        Cpos = 'rho'
        z = src_grd.z_t
        Mp, Lp = dst_grd.hgrid.mask_rho.shape
        wts_file = 'remap_weights_SODA_2.1.6_to_ARCTIC2_bilinear_t_to_rho.nc'
        dst_varname = 'salinity'
        dst_varname_north = 'salt_north'
        dimensions_north = ('ocean_time', 's_rho', 'xi_rho')
        long_name_north = 'salinity north boundary condition'
        field_north = 'salt_north, scalar, series'
        dst_varname_south = 'salt_south'
        dimensions_south = ('ocean_time', 's_rho', 'xi_rho')
        long_name_south = 'salinity south boundary condition'
        field_south = 'salt_south, scalar, series'
        dst_varname_east = 'salt_east'
        dimensions_east = ('ocean_time', 's_rho', 'eta_rho')
        long_name_east = 'salinity east boundary condition'
        field_east = 'salt_east, scalar, series'
        dst_varname_west = 'salt_west'
        dimensions_west = ('ocean_time', 's_rho', 'eta_rho')
        long_name_west = 'salinity west boundary condition'
        field_west = 'salt_west, scalar, series'
        units = 'PSU'
    else:
        raise ValueError('Undefined src_varname')


    if ndim == 3:
        # build intermediate zgrid
        zlevel = -z[::-1,0,0]
        nzlevel = len(zlevel)
        dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
        dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)


    # create variable in boudary file
    print('Creating variable', dst_varname_north)
    nc.createVariable(dst_varname_north, 'f8', dimensions_north, fill_value=spval)
    nc.variables[dst_varname_north].long_name = long_name_north
    nc.variables[dst_varname_north].units = units
    nc.variables[dst_varname_north].field = field_north

    print('Creating variable', dst_varname_south)
    nc.createVariable(dst_varname_south, 'f8', dimensions_south, fill_value=spval)
    nc.variables[dst_varname_south].long_name = long_name_south
    nc.variables[dst_varname_south].units = units
    nc.variables[dst_varname_south].field = field_south

    print('Creating variable', dst_varname_east)
    nc.createVariable(dst_varname_east, 'f8', dimensions_east, fill_value=spval)
    nc.variables[dst_varname_east].long_name = long_name_east
    nc.variables[dst_varname_east].units = units
    nc.variables[dst_varname_east].field = field_east

    print('Creating variable', dst_varname_west)
    nc.createVariable(dst_varname_west, 'f8', dimensions_west, fill_value=spval)
    nc.variables[dst_varname_west].long_name = long_name_west
    nc.variables[dst_varname_west].units = units
    nc.variables[dst_varname_west].field = field_west


    # remapping
    print('remapping', dst_varname, 'from', src_grd.name, \
              'to', dst_grd.name)
    print('time =', time)


    if ndim == 3:
        # flood the grid
        print('flood the grid')
        print(src_var.shape)
        print(src_grd.z_t.shape)
        src_varz = pyroms_toolbox.BGrid_SODA.flood(src_var, src_grd, Bpos=Bpos, spval=spval, \
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

    # write data in destination file
    print('write data in destination file')
    nc.variables['ocean_time'][0] = time
    nc.variables[dst_varname_north][0] = np.squeeze(dst_var_north)
    nc.variables[dst_varname_south][0] = np.squeeze(dst_var_south)
    nc.variables[dst_varname_east][0] = np.squeeze(dst_var_east)
    nc.variables[dst_varname_west][0] = np.squeeze(dst_var_west)

    # close file
    nc.close()
    cdf.close()

    if src_varname == 'ssh':
        return dst_varz
