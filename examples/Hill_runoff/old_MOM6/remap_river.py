import numpy as np
import os
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import matplotlib.pyplot as plt
import time
import datetime
from matplotlib.dates import date2num, num2date

import pyroms
import pyroms_toolbox
from pyroms import _remapping
import xarray as xr

class nctime(object):
    pass

def remap_river(src_file, src_varname, dst_grd, dst_mask, dst_dir='./'):

    ocn_lon = xr.open_dataset(dst_grd).x[1::2,1::2]  # Cell-center longitudes (cell centers)
    ocn_lat = xr.open_dataset(dst_grd).y[1::2,1::2]  # Cell-center latitudes (cell centers)
    lon_b = xr.open_dataset(dst_grd).x[::2,::2]      # Corner longitudes
    lat_b = xr.open_dataset(dst_grd).y[::2,::2]      # Corner latitudes
    coords = xr.Dataset({"lon": ocn_lon, "lat": ocn_lat})
    coords2 = xr.Dataset({"lon_b": lon_b, "lat_b": lat_b })
    coords = coords.rename({'nyp': 'ny', 'nxp': 'nx'})
    coords = coords.merge(coords2)

    dx = xr.open_dataset(dst_grd).dx[1::2,::2] + xr.open_dataset(dst_grd).dx[1::2,1::2]
    dy = xr.open_dataset(dst_grd).dy[::2,1::2] + xr.open_dataset(dst_grd).dy[1::2,1::2]
    area_dst = dx.data * dy.data
    mask_rho = xr.open_dataset(dst_mask).mask

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    time0 = datetime.date(1900, 1, 1)

    # create runoff file
    dst_file = src_file.rsplit('/')[-1]
    dst_file = dst_dir + dst_file[:-3] + '_' + src_varname + '_Hill_NGOA.nc'
    print('\nCreating file', dst_file)
    if os.path.exists(dst_file) is True:
        os.remove(dst_file)

    # open River file
    nc = netCDF.Dataset(dst_file, 'w', format='NETCDF4')
    nc.Author = 'pyroms_toolbox.nc_create_roms_file'
    nc.Created = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'MOM6 runoff file'

    Mp, Lp = (lon_b.shape)
    nc.createDimension('IQ', Lp)
    nc.createDimension('JQ', Mp)
    nc.createDimension('i', Lp-1)
    nc.createDimension('j', Mp-1)
    nc.createDimension('time', None)

    #load var
    cf = xr.open_dataset(src_file)
    years = cf.year.data
    months = cf.month.data
    days = cf.day.data
    src_lon = cf.lon.data + 360.0
    src_lat = cf.lat.data

    src_var_all = cf.q.data

    #get missing value
    spval = 1.e30

#   dst_varname = 'Runoff_raw'
#   dimensions = ('time', 'j', 'i')
#   long_name = 'river discharge'
#   units = 'meter^3/sec'

    # create variable in file
    nc.createVariable('time', 'f8', ('time'))
    nc.variables['time'].units = 'days since 1900-01-01'
    nc.variables['time'].calendar = 'gregorian'
#   print('Creating variable', dst_varname)
#   nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
#   nc.variables[dst_varname].long_name = long_name
#   nc.variables[dst_varname].units = 'm^3/day'

    nc.createVariable('Runoff', 'f8', ('time', 'j', 'i'), fill_value=spval)
    nc.variables['Runoff'].long_name = 'Hill River Runoff'
    nc.variables['Runoff'].units = 'kg/m^2/sec'

    # get littoral (here just 1 cell wide, with diagonals)
    width = 1
    idx = []
    idy = []
    maskl = mask_rho.copy()
    for w in range(width):
        lit = pyroms_toolbox.get_littoral(maskl)
        idx.extend(lit[0])
        idy.extend(lit[1])
        maskl[lit] = 0

    littoral_idx = (np.array(idx), np.array(idy))

    ntimes = len(years)
    for it in range(ntimes):
        src_var = src_var_all[it,:]/86400.0
        net_flow = np.sum(src_var)
        src_var = xr.DataArray(src_var)
        src_var = src_var.fillna(0)
        src_var = src_var.data
        print(src_var.shape)

        time = datetime.date(years[it], months[it], days[it])
        print('time =', time)

        # horizontal interpolation
        print('horizontal interpolation using brute force')


        # write data in destination file
        print('write data in destination file')
        nc.variables['time'][it] = (time - time0).days

        runoff = pyroms_toolbox.remap_river(src_var, src_lon, src_lat,  \
                  np.array(littoral_idx).T + 1, \
                  ocn_lon, ocn_lat)
        net_flow2 = np.sum(runoff)
        print("flow before", net_flow)
        print("regridded flow", net_flow2)

        # Runoff_raw for debugging
#       nc.variables[dst_varname][it] = runoff
## zero out contribution from outside the domain
        runoff[0,:] = 0
        runoff[:,-1] = 0
        net_flow3 = np.sum(runoff)
        print("final flow", net_flow3)
        runoff = runoff * 1000.0 / area_dst  # convert units
        nc.variables['Runoff'][it] = runoff

    # close destination file
    nc.close()

