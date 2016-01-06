import matplotlib
matplotlib.use('Agg')

import numpy as np
import netCDF4
from datetime import datetime
import pyroms
import pyroms_toolbox
import sys



# get HYCOM Northeast Pacific data from 2007 to 2011

invarname = 'temperature'
outvarname = 'temp'

#read grid and variable attributes from the first file
url='http://tds.hycom.org/thredds/dodsC/datasets/GLBa0.08/expt_91.1/2015/temp/archv.2015_001_00_3zt.nc'
dataset = netCDF4.Dataset(url)
lon = dataset.variables['Longitude'][1500-10:1800+1,600-1:940+1]
lat = dataset.variables['Latitude'][1500-10:1800+1,600-1:940+1]
z = dataset.variables['Depth'][:]
#spval = dataset.variables[invarname]._FillValue
var = dataset.variables[invarname][0,:,1500-10:1800+1,600-1:940+1]
spval = var.get_fill_value()
units = dataset.variables[invarname].units
long_name = dataset.variables[invarname].long_name
dataset.close()


year = 2007
day = 1


#create netCDF file
outfile = 'HYCOM_GLBa0.08_PALAU_grid.nc'
nc = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')
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
jday = pyroms_toolbox.date2jday(datetime(year, 1, 1)) + day - 1
nc.variables['ocean_time'][0] = jday

nc.createVariable(outvarname, 'f', ('ocean_time', 'z', 'lat', 'lon'), fill_value=spval)
nc.variables[outvarname].long_name = long_name
nc.variables[outvarname].units = units
nc.variables[outvarname][0] = var
        
nc.close()





