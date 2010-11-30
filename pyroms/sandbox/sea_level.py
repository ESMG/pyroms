
# sample url
#http://tidesandcurrents.noaa.gov/data_listing.shtml?bdate=20080323&edate=20080423&datum=6&unit=0&shift=g&stn=8762075&type=Tide%20Data&format=View+Data&listing=1

import urllib2
from datetime import datetime

import numpy as np

class sea_level(object):
    """docstring for sea_level"""
    
    data_dict =  {'bdate'   : '20080323',
                  'edate'   : '20080423',
                  'datum'   : '6',
                  'unit'    : '0',
                  'shift'   : 'g',
                  'stn'     : '8762075',
                  'type'    : 'Tide%20Data',
                  'format'  : 'View+Data',
                  'listing' : '1'}
    
    root = 'http://tidesandcurrents.noaa.gov/data_listing.shtml?'
    
    def __init__(self, station_id, start_date, end_date=None):
        if end_date is None:
            end_date = datetime.now()
        self.data_dict['stn'] = str(station_id)
        
        url = self.root + '&'.join([ '='.join(keyval) for (keyval) in 
                                    self.data_dict.iteritems()])
        lines = urllib2.urlopen(url).readlines()
        
        pr = False
        date = []
        ssh = []
        for line in lines:
            data = line.split()
            if len(data) == 0:
                continue
            if data[0] == self.data_dict['stn']:
                yr = int(data[1][:4])
                mo = int(data[1][4:6])
                day = int(data[1][6:])
                hr = int(data[2][:2])
                mn = int(data[2][3:])
                date.append(datetime(yr, mo, day, hr, mn))
                ssh.append(float(data[3]))
            
        self.ssh = np.asarray(ssh)
        self.date = np.asarray(date)
        



sl = sea_level('8762075', 'foo')
print sl.date
print sl.ssh
        