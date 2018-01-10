import matplotlib
matplotlib.use('Agg')

import numpy as np
import netCDF4
from datetime import datetime
import pyroms
import pyroms_toolbox
import sys



def create_HYCOM_file(name, time, lon, lat, var):

    #create netCDF file
    nc = netCDF4.Dataset(name, 'w', format='NETCDF3_64BIT')
    nc.Author = sys._getframe().f_code.co_name
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'HYCOM + NCODA Global 1/12 Analysis (GLBa0.08)'

    #create dimensions
    Mp, Lp = lon.shape
    nc.createDimension('lon', Lp)
    nc.createDimension('lat', Mp)
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

    nc.createVariable('ocean_time', 'f', ('ocean_time'))
    nc.variables['ocean_time'].units = 'days since 1900-01-01 00:00:00'
    nc.variables['ocean_time'].calendar = 'LEAP'
    nc.variables['ocean_time'][0] = time

    nc.createVariable(outvarname, 'f', ('ocean_time', 'lat', 'lon'), fill_value=spval)
    nc.variables[outvarname].long_name = long_name
    nc.variables[outvarname].units = units
    nc.variables[outvarname].coordinates = 'lon lat'
    nc.variables[outvarname][0] = var

    nc.close()

    print('Done with file %s' %name)




# get HYCOM Northeast Pacific data from 2007 to 2011

year = 2015
retry='True'

invarname = 'ssh'
outvarname = 'ssh'

#read grid and variable attributes from the first file
url='http://tds.hycom.org/thredds/dodsC/datasets/GLBa0.08/expt_91.1/2015/2d/archv.2015_001_00_2d.nc'
dataset = netCDF4.Dataset(url)
lon = dataset.variables['Longitude'][1500-9:1800,600:940]
lat = dataset.variables['Latitude'][1500-9:1800,600:940]
#spval = dataset.variables[invarname]._FillValue
units = dataset.variables[invarname].units
long_name = dataset.variables[invarname].long_name
dataset.close()


retry_day = []

# loop over daily files
if year%4 == 0:
    daysinyear = 366
else:
#    daysinyear = 365
    daysinyear = 32
for day in range(1,daysinyear+1):
    print('Processing file for day %03d, year %04d' %(day, year))
    url='http://tds.hycom.org/thredds/dodsC/datasets/GLBa0.08/expt_91.1/2015/2d/archv.%04d_%03d_00_2d.nc' %(year,day)
    #get data from server
    try:
        dataset = netCDF4.Dataset(url)
        var = dataset.variables[invarname][0,1500-9:1800,600:940]
        spval = var.get_fill_value()
        dataset.close()
    except:
        print('No file on the server... We skip this day.')
        retry_day.append(day)
        continue

    #create netCDF file
    outfile = 'data/HYCOM_GLBa0.08_%s_%04d_%03d.nc' %(outvarname,year,day)
    jday = pyroms_toolbox.date2jday(datetime(year, 1, 1)) + day - 1
    create_HYCOM_file(outfile, jday, lon, lat, var)


if retry == 'True':
    if len(retry_day) != 0:
        print("Some files have not been downloded... Let's try again")
    while len(retry_day) != 0:
        for day in retry_day:
            print('Retry file for day %03d, year %04d' %(day, year))
            url='http://tds.hycom.org/thredds/dodsC/datasets/GLBa0.08/expt_91.1/2015/2d/archv.%04d_%03d_00_2d.nc' %(year,day)
            #get data from server
            try:
                dataset = netCDF4.Dataset(url)
                var = dataset.variables[invarname][0,1500-9:1800,600:940]
                spval = var.get_fill_value()
                dataset.close()
            except:
                print('No file on the server... We skip this day.')
            continue

            #create netCDF file
            outfile = 'data/HYCOM_GLBa0.08_%s_%04d_%03d.nc' %(outvarname,year,day)
            jday = pyroms_toolbox.date2jday(datetime(year, 1, 1)) + day - 1
            create_HYCOM_file(outfile, jday, lon, lat, var)

            retry_day.remove(day)


