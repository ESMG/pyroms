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
dst_grd = '/import/c1/AKWATERS/kate/ESMG/ESMG-configs/NGOA/INPUT/ocean_hgrid.nc'
dst_mask = '/import/c1/AKWATERS/kate/ESMG/ESMG-configs/NGOA/INPUT/ocean_mask.nc'
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


for file in lst_file:
# remap
    remap_river(file, 'q', dst_grd, dst_mask, dst_dir=dst_dir)
