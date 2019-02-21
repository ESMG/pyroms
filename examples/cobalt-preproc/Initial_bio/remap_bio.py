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

def remap_bio(argdict, src_grd, dst_grd, dmax=0, cdepth=0, kk=0, dst_dir='./'):

    # NWGOA3 grid sub-sample
    xrange=src_grd.xrange; yrange=src_grd.yrange

    src_varname = argdict['tracer']
    tracer = src_varname
    src_file = argdict['file']
    units = argdict['units']
    longname = argdict['longname']
    nframe = argdict['nframe']

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'

    # create clim file
    dst_file = tracer + '.nc'
    dst_file = dst_dir + dst_grd.name + '_ic_bio_' + dst_file
    print('Creating clim file', dst_file)
    if os.path.exists(dst_file) is True:
        os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)

    # open clim file
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')

    #load var
    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]

    time = cdf.variables['time'][nframe]

    # to be in sync with physics, add +0.5 day
    #time = time + 0.5
    # time will be given by physics anyway

    #get missing value
    spval = src_var._FillValue

    # determine variable dimension
    ndim = len(src_var.dimensions) - 1

    # NWGOA3 grid sub-sample
    if ndim == 3:
        src_var = src_var[nframe,:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    elif ndim == 2:
        src_var = src_var[nframe,yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]


    Bpos = 't'
    Cpos = 'rho'
    z = src_grd.z_t
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    wts_file = 'remap_weights_ESM2M_to_NWGOA3_bilinear_t_to_rho.nc'
    dst_varname = tracer
    dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
    long_name = longname
    field = tracer + ', scalar, series'
    units = units

    if ndim == 3:
        # build intermediate zgrid
        zlevel = -z[::-1]
        nzlevel = len(zlevel)
        dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
        dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)


    # create variable in file
    print('Creating variable', dst_varname)
    nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    nc.variables[dst_varname].long_name = long_name
    nc.variables[dst_varname].units = units
    nc.variables[dst_varname].field = field
    #nc.variables[dst_varname_north]._FillValue = spval


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
        dst_var = pyroms.remapping.z2roms(dst_varz[::-1,:,:], dst_grdz, \
                          dst_grd, Cpos=Cpos, spval=spval, flood=False)
    else:
        dst_var = dst_varz

    # write data in destination file
    print('write data in destination file\n')
    nc.variables['ocean_time'][0] = time
    nc.variables[dst_varname][0] = dst_var

    # close file
    nc.close()
    cdf.close()

    if src_varname == 'eta':
        return dst_varz
