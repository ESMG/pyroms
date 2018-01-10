import matplotlib
matplotlib.use('Agg')
import subprocess
from multiprocessing import Pool
import os
import sys
import string

import pyroms
from station_bound import *
import subprocess
import pdb

irange = None
jrange = None

def do_file(file):
    print('file is: ' + file)
    var_list = ['u', 'v', 'temp', 'salt', 'zeta']
#    pdb.set_trace()
    dst_var = station_bound(var_list, file,\
                     src_grd, dst_grd)

# Change src_filename to your directory for the file's containing variable data
data_dir = '/archive/u1/uaf/kate/NGOA/run05/'
lst = subprocess.getoutput('ls ' + data_dir + 'nwgoa_sta.nc')
lst_file = lst.split()

src_grd = pyroms.sta_grid.get_Stations_grid('NWGOA3', lst_file[0])
dst_grd = pyroms.grid.get_ROMS_grid('GLACIER_BAY')
# Outfile is a parameter to allow you to place these created remap files in a different
# directory than the one that is default which is where the file came from.
#processes = 1
##processes = 2
#p = Pool(processes)
#results = p.map(do_file, lst_file)

for file in lst_file:
    do_file(file)
