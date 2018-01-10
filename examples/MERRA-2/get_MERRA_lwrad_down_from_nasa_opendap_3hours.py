import matplotlib
matplotlib.use('Agg')
import numpy as np
import netCDF4
from datetime import datetime
import pyroms
import pyroms_toolbox

import http.cookiejar
import netrc
import urllib.request, urllib.error, urllib.parse
import re
import pydap.lib
from pydap.exceptions import ClientError
#import logging
import sys

#logging.info("Starting logger")

#log = logging.getLogger(__name__)
#log.setLevel(logging.DEBUG)

# Set the debug level for urllib2.
debuglevel=1

def install_basic_client(uri='', user='', passwd='', use_netrc=True):
    # Create special opener with support for Cookies
    cj = http.cookiejar.CookieJar()
    # Create the password manager and load with the credentials using
    pwMgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    # Get passwords from the .netrc file nless use_netrc is False
    if use_netrc:
        logins = netrc.netrc()
        accounts = logins.hosts # a dist of hosts and tuples
        for host, info in accounts.items():
            login, account, password = info
#            log.debug('Host: %s; login: %s; account: %s; password: %s' % (host, login, account, password))
            pwMgr.add_password(None, host, login, password)
    if uri and user and passwd:
        pwMgr.add_password(None, uri, user, passwd)
    opener = urllib.request.build_opener(urllib.request.HTTPBasicAuthHandler(pwMgr), urllib.request.HTTPCookieProcessor(cj))
    opener.addheaders = [('User-agent', pydap.lib.USER_AGENT)]
    urllib.request.install_opener(opener)
    def new_request(url):
        if url[-1] is '&': url = url[0:-1]
#        log.debug('Opening %s (install_basic_client)' % url)
        r = urllib.request.urlopen(url)
        resp = r.headers.dict
        resp['status'] = str(r.code)
        data = r.read()
        # When an error is returned, we parse the error message from the
        # server and return it in a ``ClientError`` exception.
        if resp.get("content-description") == "dods_error":
            m = re.search('code = (?P<code>\d+);\s*message = "(?P<msg>.*)"', data, re.DOTALL | re.MULTILINE)
            msg = 'Server error %(code)s: "%(msg)s"' % m.groupdict()
            raise ClientError(msg)
        return resp, data
    from pydap.util import http
    http.request = new_request
# END BASIC AUTH MODULE CODE

install_basic_client()
from pydap.client import open_url

year = int(sys.argv[1])

invarname = 'LWGAB'

outvarname = 'lwrad_down'
outtimename = 'lrf_time'

server = 'http://goldsmr4.sci.gsfc.nasa.gov'

leap = year%4
if leap == 0:
    daysinmonth = ([31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    daysinyear = 366
else:
    daysinmonth = ([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    daysinyear = 365

#if (year <= 1992):
#    file_tag = 'MERRA101'
#elif ((year >= 1993) & (year <= 2000)):
#    file_tag = 'MERRA201'
#elif ((year >= 2001) & (year <= 2009)):
#    file_tag = 'MERRA301'
#elif (year >= 2010):
#    file_tag = 'MERRA300'
file_tag = 'MERRA2_400'

#read grid and variable attributes from the first file
year_tag = '%04d' %year
month_tag = '01'
day_tag = '01'
date_tag = year_tag + month_tag + day_tag
url = server + '/opendap/MERRA2/M2T1NXRAD.5.12.4/' + year_tag + '/' + month_tag + '/' + \
      file_tag + '.tavg1_2d_rad_Nx.' + date_tag + '.nc4'
dataset = open_url(url)
lon = dataset['lon'][:]
#shift data between 0 and 360 deg.
gidx = np.where(np.abs(lon) < 1.0e-10)[0][0]
lon = lon + 180
lat = dataset['lat'][:]
spval = dataset[invarname].missing_value
units = dataset[invarname].units
long_name = dataset[invarname].long_name

#get data from NASA opendap
for month in range(12):
    nday = 0

    #create ROMS forcing file
    month_tag = '%02d' %(month+1)
    outfile = 'Forcings/MERRA_' + outvarname + '_3hours_' + year_tag + '_' + month_tag + '.nc'
    nc = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')
    nc.Author = sys._getframe().f_code.co_name
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'MERRA-2 dataset. Modern Era Retrospective-analysis'

    nc.createDimension('lon', np.size(lon))
    nc.createDimension('lat', np.size(lat))
    nc.createDimension(outtimename, None)

    nc.createVariable('lon', 'f8', ('lon'))
    nc.variables['lon'].long_name = 'longitude'
    nc.variables['lon'].units = 'degrees_east'
    nc.variables['lon'][:] = lon

    nc.createVariable('lat', 'f8', ('lat'))
    nc.variables['lat'].long_name = 'latitude'
    nc.variables['lat'].units = 'degrees_north'
    nc.variables['lat'][:] = lat

    nc.createVariable(outtimename, 'f8', (outtimename))
    nc.variables[outtimename].units = 'days since 1900-01-01 00:00:00'
    nc.variables[outtimename].calendar = 'LEAP'
    dstart = pyroms_toolbox.date2jday(datetime(year, month+1, 1, 1, 30))
    roms_time = np.arange(dstart, dstart+daysinmonth[month], 3./24)
    nc.variables[outtimename][:] = roms_time

    nc.createVariable(outvarname, 'f', (outtimename, 'lat', 'lon'), fill_value=spval)
    nc.variables[outvarname].missing_value = spval
    nc.variables[outvarname].long_name = long_name
    nc.variables[outvarname].units = units
    nc.variables[outvarname].coordinates = 'lon lat'

    for day in range(daysinmonth[month]):
#    if year == 2010:
#      if ((month+1 >= 6) & (month+1 <= 8)):
#        file_tag = 'MERRA301'
#      else:
#        file_tag = 'MERRA300'
        day_tag = '%02d' %(day+1)
        date_tag = year_tag + month_tag + day_tag
        url = server + '/opendap/MERRA2/M2T1NXRAD.5.12.4/' + year_tag + '/' + month_tag + '/' + \
              file_tag + '.tavg1_2d_rad_Nx.' + date_tag + '.nc4'
        dataset = open_url(url)
        var = dataset[invarname][:]
        #shift data between 0 and 360 deg.
        svar = np.zeros(var.shape)
        svar[:,:,:len(lon)-gidx] = var[:,:,gidx:]
        svar[:,:,len(lon)-gidx:] = var[:,:,:gidx]
        var_3hr = (svar[::3] + svar[1::3] + svar[2::3]) / 3.
        nc.variables[outvarname][nday*8:(nday+1)*8,:,:] = var_3hr
        nday = nday + 1
#       dataset.close()

    nc.close()
