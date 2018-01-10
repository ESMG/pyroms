import matplotlib
matplotlib.use('Agg')

import numpy as np
import netCDF4
from datetime import datetime
import pyroms
import pyroms_toolbox
import sys

'''
Usage: python get_pacific_3d.py [year] [field]
    For instance: python get_pacific_3d.py 1991 sio4
'''

year = int(sys.argv[1])
field = sys.argv[2]

invarname = field
outvarname = field

def create_HYCOM_file(name):
    global nc
    print('Creating file %s' %name)

    #create netCDF file
    nc = netCDF4.Dataset(name, 'w', format='NETCDF3_64BIT')
    nc.Author = sys._getframe().f_code.co_name
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'ROMS + CoSiNE Pacific (U Maine)'

    #create dimensions
    Mp, Lp = zeta0.shape
    N = 20
    theta_s = 5
    theta_b = 0
    Tcline = 50
    nc.createDimension('xi_rho', Lp)
    nc.createDimension('eta_rho', Mp)
    nc.createDimension('s_rho', N)
    nc.createDimension('ocean_time', None)

    #create variables        
#    nc.createVariable('lon_rho', 'f', ('eta_rho', 'xi_rho'))
#    nc.variables['lon_rho'].long_name = 'longitude'
#    nc.variables['lon_rho'].units = 'degrees_east'
#    nc.variables['lon_rho'][:] = lon
#
#    nc.createVariable('lat_rho', 'f', ('eta_rho', 'xi_rho'))
#    nc.variables['lat_rho'].long_name = 'latitude'
#    nc.variables['lat_rho'].units = 'degrees_north'
#    nc.variables['lat_rho'][:] = lat
#
#    nc.createVariable('s_rho', 'f', ('s_rho'))
#    nc.variables['s_rho'].standard_name = 'ocean_s_coordinate'
#    nc.variables['s_rho'].formula_terms = 's: s_rho eta: zeta depth: h a: theta_s b: theta_b depth_c: Tcline'
#    nc.variables['s_rho'][:] = s_rho

    nc.createVariable('theta_s', 'f')
    nc.variables['theta_s'][:] = theta_s

    nc.createVariable('theta_b', 'f')
    nc.variables['theta_b'][:] = theta_b

    nc.createVariable('Tcline', 'f')
    nc.variables['Tcline'][:] = Tcline

    nc.createVariable('ocean_time', 'f', ('ocean_time'))
    nc.variables['ocean_time'].units = 'days since 1900-01-01 00:00:00'
    nc.variables['ocean_time'].calendar = 'LEAP'
    nc.variables['ocean_time'].long_name = 'ocean time'
#    nc.variables['ocean_time'][:] = time

    nc.createVariable(outvarname, 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
    nc.variables[outvarname].long_name = long_name
    nc.variables[outvarname].units = units
    nc.variables[outvarname].field = field
    nc.variables[outvarname].time = 'ocean_time'
    nc.variables[outvarname].coordinates = 'ocean_time s_rho lon_rho lat_rho'


    print('Done with header for file %s' %name)



# get Pacific data from 1991 to 2008

#year = 1991
retry='True'
spval = 1.e30
rec_start = (year-1948)*120
rec_end = rec_start + 120

#read grid and variable attributes from the first file
#url='http://viz.clusters.umaine.edu:8080/thredds/dodsC/pacific/1991-2008'
url='http://viz.clusters.umaine.edu:8080/thredds/dodsC/pacific/roms50-avg'
dataset = netCDF4.Dataset(url)
zeta0 = dataset.variables['zeta'][0,170:215,195:265]
#lon = dataset.variables['lon_rho'][170:215,195:265]
#lat = dataset.variables['lat_rho'][170:215,195:265]
#h = dataset.variables['h'][170:215,195:265]
#theta_s = dataset.variables['theta_s'][:]
#theta_b = dataset.variables['theta_b'][:]
#Tcline = dataset.variables['Tcline'][:]
#s_rho = dataset.variables['s_rho'][:]
time = dataset.variables['scrum_time'][:]
#spval = dataset.variables[invarname]._FillValue
units = dataset.variables[invarname].units
field = dataset.variables[invarname].field
long_name = dataset.variables[invarname].long_name
#dataset.close()


retry_day = []

# loop over records, 73 hours apart
outfile = 'data/Pacific_%s_%04d.nc' %(outvarname,year)
create_HYCOM_file(outfile)
day_out = 0
for day in range(rec_start,rec_end):
    print('Processing file for %s, day %d, year %04d' %(invarname, day_out*3, year))
    #get data from server
    try:
        var = dataset.variables[invarname][day,:,170:215,195:265]
#        spval = var.get_fill_value()
#        dataset.close()
        print('Got %s from server...' %invarname)
    except:
        print('No file on the server... We skip this day.')
        retry_day.append((day,day_out))
        continue

    #create netCDF file
    nc.variables[outvarname][day_out] = var
    jday = pyroms_toolbox.date2jday(datetime(year, 1, 1)) + day_out*73/24.
    nc.variables['ocean_time'][day_out] = jday
    day_out += 1


if retry == 'True':
    if len(retry_day) != 0:
        print("Some file have not been downloded... Let's try again")
    while len(retry_day) != 0:
        for (day,day_out) in retry_day:
            print('Retry file for %s, day %03d, year %04d' %(invarname, day_out, year))
            #get data from server
            try:
                var = dataset.variables[invarname][day,:,170:215,195:265]
#                spval = var.get_fill_value()
                print('Got %s from server...' %invarname)
            except:
                print('No file on the server... We skip this day.')
                continue

            #create netCDF file
            jday = pyroms_toolbox.date2jday(datetime(year, 1, 1)) + day_out*73/24.
            nc.variables[outvarname][day_out] = var
            nc.variables['ocean_time'][day_out] = jday

            retry_day.remove((day,day_out))


dataset.close()
nc.close()
