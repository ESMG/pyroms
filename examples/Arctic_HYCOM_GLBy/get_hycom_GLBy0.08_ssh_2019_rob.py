import matplotlib
matplotlib.use('Agg')

import numpy as np
import netCDF4
from netCDF4 import num2date, date2num
from datetime import datetime, timedelta
import pyroms
import pyroms_toolbox
import os, sys, time

def create_HYCOM_file(name, time, lon, lat, vard, spval):

    print('Write with file %s' %name)

    #create netCDF file
    nc = netCDF4.Dataset(name, 'w', format='NETCDF3_64BIT')
    nc.Author = sys._getframe().f_code.co_name
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'HYCOM + NCODA Global 1/12 Analysis (GLBy0.08) expt_93.0'

    #create dimensions
    #Mp, Lp = lon.shape
    Lp = len(lon)
    Mp = len(lat)
    nc.createDimension('lon', Lp)
    nc.createDimension('lat', Mp)
    nc.createDimension('ocean_time', None)

    #create variables
    #nc.createVariable('lon', 'f', ('lat', 'lon'))
    nc.createVariable('lon', 'f', ('lon'))
    nc.variables['lon'].long_name = 'longitude'
    nc.variables['lon'].units = 'degrees_east'
    nc.variables['lon'][:] = lon

    #nc.createVariable('lat', 'f', ('lat', 'lon'))
    nc.createVariable('lat', 'f', ('lat'))
    nc.variables['lat'].long_name = 'latitude'
    nc.variables['lat'].units = 'degrees_north'
    nc.variables['lat'][:] = lat

    nc.createVariable('ocean_time', 'f', ('ocean_time'))
    nc.variables['ocean_time'].units = 'days since 1900-01-01 00:00:00'
    nc.variables['ocean_time'].calendar = 'gregorian'
    nc.variables['ocean_time'][0] = time

    nc.createVariable(outvarname, 'f', ('ocean_time', 'lat', 'lon'), fill_value=spval)
    nc.variables[outvarname].long_name = long_name
    nc.variables[outvarname].units = units
    nc.variables[outvarname].coordinates = 'lon lat'
    nc.variables[outvarname][0] = vard

    nc.close()

    print('Done with file %s' % name)




# get HYCOM Northeast Pacific data from 2019

year = 2019

invarname = 'surf_el'
outvarname = 'ssh'

#read grid and variable attributes from thredds server
url='http://tds.hycom.org/thredds/dodsC/datasets/GLBy0.08/expt_93.0/data/hindcasts/2019/hycom_glby_930_2019120412_t000_ssh.nc'
dataset = netCDF4.Dataset(url)
lon = dataset.variables['lon'][:]
# Lat 45.0+N
lat = dataset.variables['lat'][:]
lat_coord = np.where(lat>=45.0)
latStart = lat_coord[0][0]
latEnd = lat_coord[0][-1]+1
lat = dataset.variables['lat'][latStart:latEnd]
mt = dataset.variables['time'][:]
mtStart = mt[0]
mtEnd = mt[-1]
nowNum = mtStart
timeUnits = dataset.variables['time'].units
# Convert times and units to dates we can use
mtDates = num2date(mt,timeUnits)
units = dataset.variables[invarname].units
long_name = dataset.variables[invarname].long_name
dataset.close()

# loop over daily files
if year%4 == 0:
    daysinyear = 366
else:
    daysinyear = 365

# Use python dates to specify start and stop time
startTime = datetime(year, 1, 1, 0, 0)
iTime = startTime
endTime = datetime(year, 2, 5, 0, 0)
#endTime = datetime(year, 12, 31, 0, 0)
missingDates = [
]
# Convert start and end Time(s) to indicies from thredds server
#startNum = date2num(startTime,timeUnits)
#endNum = date2num(endTime,timeUnits)
# we assume numeric times within the dataset on the server are
# consecutive and that none are missing.
#if startNum < mtStart:
#    print("Requested start: %s (index=%d)" % (startTime.strftime("%Y-%m-%d %H:%M:%S"),startNum))
#    print("Dataset start: %s (index=%d)" % (mtDates[0].strftime("%Y-%m-%d %H:%M:%S"),mtStart))
#    print("Dataset reported time units: %s" % (timeUnits))
#    sys.exit("Start date is outside available time index of dataset.")
#if endNum > mtEnd:
#    print("Requested end: %s (index=%d)" % (endTime.strftime("%Y-%m-%d %H:%M:%S"),endNum))
#    print("Dataset end: %s (index=%d)" % (mtDates[-1].strftime("%Y-%m-%d %H:%M:%S"),mtEnd))
#    print("Dataset reported time units: %s" % (timeUnits))
#    sys.exit("End date is outside available time index of dataset.")

#sys.exit()
#import pdb; pdb.set_trace()

# Setup the correct indicies to request from the server
# 0 for thredds is the first index
#istartNum = int(startNum) - int(mtStart)
#iendNum = int(endNum) - int(mtStart)

# Thredds server limits amount of data that can be recalled at one time

#for day in range(1,daysinyear+1):
#for datasetIndex in range(istartNum,iendNum+1):
while iTime <= endTime:

    #nowNum = datasetIndex + int(mtStart)

    # create netCDF output filename
    #day = nowNum - int(date2num(datetime(year, 1, 1, 0, 0),timeUnits)) + 1
    day = iTime.timetuple().tm_yday
    outfile = 'data/HYCOM_GLBy0.08_%s_%04d_%03d.nc' % (outvarname,year,day)

    # sometimes years have specific days missing
    if iTime in missingDates:
        iTime = iTime + timedelta(days=1)
        continue

    # skip files that exist
    if os.path.isfile(outfile):
        iTime = iTime + timedelta(days=1)
        continue

    print('Processing file for %s, %s' % (invarname, iTime.strftime("%Y-%m-%d")))

    #get data from server
    url='http://tds.hycom.org/thredds/dodsC/datasets/GLBy0.08/expt_93.0/data/hindcasts/%04d/hycom_glby_930_%s12_t000_ssh.nc' % (year, iTime.strftime("%Y%m%d"))
    writeFile = True
    for layerChunk in range(0,1):

        retries = 0
        success = False

        while success == False:
            try:
                dataset = netCDF4.Dataset(url)
                # datasetIndex,lat,lon
                var = dataset.variables[invarname][:,latStart:latEnd,:]
                #print(var.shape)

                # Server sometimes delivers blocks back with all zeros.
                # If server is returning zeros, try to give it a [90*number_of_failures second] break.
                vmin = var.min()
                vmax = var.max()
                if vmin == 0.0 and vmax == 0.0:
                    print("Field arrived as all zeros, retrying...")
                    success = False
                    retries = retries + 1
                    time.sleep(90*retries)
                    continue
                if layerChunk == 0:
                  vard = var
                  spval = var.get_fill_value()
                else:
                  vard = np.concatenate((vard,var),axis=0)
                #print(vard.shape)
                dataset.close()
                print('Got %s(chunk %d) from server...' % (invarname,layerChunk+1))
                success = True
            except:
                # On any server failure, wait 120 seconds.
                err = sys.exc_info()[1]
                # If it is a missing data file, automatically skip it
                if err.strerror == 'NetCDF: file not found':
                    #iTime = iTime + timedelta(days=1)
                    writeFile = False
                    success = True
                else:
                    print("Unexpected error:",err.strerror)
                    print('Server failed to deliver data, retrying...')
                    retries = retries + 1
                    time.sleep(120)
                    continue
                #import pdb; pdb.set_trace()

    if writeFile:
        #create netCDF file
        jday = pyroms_toolbox.date2jday(datetime(year, 1, 1)) + day - 1
        create_HYCOM_file(outfile, jday, lon, lat, vard, spval)
    else:
        print("File not found on server %s" % (url))

    iTime = iTime + timedelta(days=1)
