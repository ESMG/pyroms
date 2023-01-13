import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
import numpy as np

import pyroms
import pyroms_toolbox

from remap_river import remap_river

lst_year = sys.argv[1:]

data_dir = '/import/AKWATERS/kshedstrom/hydroflow/new_2019/'
dst_dir='./'

lst_file = []

for year in lst_year:
    year = np.str(year)
    lst = subprocess.getoutput('ls ' + data_dir + 'goa_dischargex_*' + year + '.nc')
    lst = lst.split()
    lst_file = lst_file + lst

print('Build CLM file from the following file list:')
print(lst_file)
print(' ')

dst_grd = pyroms.grid.get_ROMS_grid('NWGOA')

for file in lst_file:
# remap
    remap_river(file, 'q', dst_grd, dst_dir=dst_dir)
