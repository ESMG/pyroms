import matplotlib
matplotlib.use('Agg')

import numpy as np
import netCDF4
from datetime import datetime
import pyroms
import pyroms_toolbox
import sys



def create_HYCOM_file(name, time, lon, lat, z, var):

    print('Write with file %s' %name)

    #create netCDF file
    nc = netCDF4.Dataset(name, 'w', format='NETCDF3_64BIT')
    nc.Author = sys._getframe().f_code.co_name
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'HYCOM + NCODA Global 1/12 Analysis (GLBa0.08)'

    #create dimensions
    Mp, Lp = lon.shape
    N = len(z)
    nc.createDimension('lon', Lp)
    nc.createDimension('lat', Mp)
    nc.createDimension('z', N)
    nc.createDimension('ocean_time', None)

    #create variables        
    nc.createVariable('lon', 'f', ('lat', 'lon'))
    nc.variables['lon'].long_name = 'longitude'
    nc.variables['lon'].units = 'degrees_east'
    nc.variables['lon'][:] = lon

    nc.createVariable('lat', 'f', ('lat', 'lon'))
    nc.variables['lat'].long_name = 'latitude'
    nc.variables['lat'].units = 'degrees_north'
    nc.variables['lat'][:] = lat

    nc.createVariable('z', 'f', ('z'))
    nc.variables['z'].long_name = 'depth'
    nc.variables['z'].units = 'meter'
    nc.variables['z'][:] = z

    nc.createVariable('ocean_time', 'f', ('ocean_time'))
    nc.variables['ocean_time'].units = 'days since 1900-01-01 00:00:00'
    nc.variables['ocean_time'].calendar = 'LEAP'
    nc.variables['ocean_time'][0] = time

    nc.createVariable(outvarname, 'f', ('ocean_time', 'z', 'lat', 'lon'), fill_value=spval)
    nc.variables[outvarname].long_name = long_name
    nc.variables[outvarname].units = units
    nc.variables[outvarname].coordinates = 'lon lat'
    nc.variables[outvarname][0] = var

    nc.close()

    print('Done with file %s' %name)




# get HYCOM Northeast Pacific data from 2007 to 2009

year = 2014
retry='True'

invarname = 'temperature'
outvarname = 'temp'

#read grid and variable attributes from the first file
url='http://tds.hycom.org/thredds/dodsC/datasets/GLBa0.08/expt_90.6/2009/temp/archv.2009_001_00_3zt.nc'
dataset = netCDF4.Dataset(url)
lon = dataset.variables['Longitude'][2100:,550:4040]
lat = dataset.variables['Latitude'][2100:,550:4040]
z = dataset.variables['Depth'][:]
#spval = dataset.variables[invarname]._FillValue
units = dataset.variables[invarname].units
long_name = dataset.variables[invarname].long_name
dataset.close()


retry_day = []

# loop over daily files
if year%4 == 0:
    daysinyear = 366
else:
    daysinyear = 365
daysinyear = 94
for day in range(1,daysinyear+1):
#for day in range(95,daysinyear+1):
    print('Processing file for %s, day %03d, year %04d' %(invarname, day, year))
    url='http://tds.hycom.org/thredds/dodsC/datasets/GLBa0.08/expt_91.0/2014/temp/archv.%04d_%03d_00_3zt.nc' %(year,day)
#    url='http://tds.hycom.org/thredds/dodsC/datasets/GLBa0.08/expt_91.1/2014/temp/archv.%04d_%03d_00_3zt.nc' %(year,day)
    #get data from server
    try:
        dataset = netCDF4.Dataset(url)
        var = dataset.variables[invarname][0,:,2100:,550:4040]
        spval = var.get_fill_value()
        dataset.close()
        print('Got %s from server...' %invarname)
    except:
        print('No file on the server... We skip this day.')
        retry_day.append(day)
        continue

    #create netCDF file
    outfile = 'data/HYCOM_GLBa0.08_%s_%04d_%03d.nc' %(outvarname,year,day)
    jday = pyroms_toolbox.date2jday(datetime(year, 1, 1)) + day - 1
    create_HYCOM_file(outfile, jday, lon, lat, z, var)


if retry == 'True':
    if len(retry_day) != 0:
        print("Some file have not been downloded... Let's try again")
    while len(retry_day) != 0:
        for day in retry_day:
            print('Retry file for %s, day %03d, year %04d' %(invarname, day, year))
            url='http://tds.hycom.org/thredds/dodsC/datasets/GLBa0.08/expt_91.0/2014/temp/archv.%04d_%03d_00_3zt.nc' %(year,day)
#            url='http://tds.hycom.org/thredds/dodsC/datasets/GLBa0.08/expt_91.1/2014/temp/archv.%04d_%03d_00_3zt.nc' %(year,day)
            #get data from server
            try:
                dataset = netCDF4.Dataset(url)
                var = dataset.variables[invarname][0,:,2100:,550:4040]
                spval = var.get_fill_value()
                dataset.close()
                print('Got %s from server...' %invarname)
            except:
                print('No file on the server... We skip this day.')
                continue

            #create netCDF file
            outfile = 'data/HYCOM_GLBa0.08_%s_%04d_%03d.nc' %(outvarname,year,day)
            jday = pyroms_toolbox.date2jday(datetime(year, 1, 1)) + day - 1
            create_HYCOM_file(outfile, jday, lon, lat, z, var)

            retry_day.remove(day)


