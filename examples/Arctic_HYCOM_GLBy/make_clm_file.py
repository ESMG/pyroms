import subprocess
import os
import sys
import subprocess
import numpy as np
from datetime import datetime
import matplotlib
matplotlib.use('Agg')

import pyroms
import pyroms_toolbox

from remap_clm import remap_clm
from remap_clm_uv import remap_clm_uv

lst_year = sys.argv[1:]

data_dir = '/import/AKWATERS/kshedstrom/HYCOM/Svalbard/Monthly_avg/'
dst_dir='./clm/'

lst_file = []

for year in lst_year:
    year = np.str(year)
    lst = subprocess.getoutput('ls ' + data_dir + '*GLBy*' + year + '*')
#   lst = subprocess.getoutput('ls ' + data_dir + '*GLBy*' + year + '*12.nc')
    lst = lst.split()
    lst_file = lst_file + lst

print('Build CLM file from the following file list:')
print(lst_file)
print(' ')

src_grd_file = data_dir + '../data/HYCOM_GLBy0.08_2018_345.nc'
src_grd = pyroms_toolbox.Grid_HYCOM.get_nc_Grid_HYCOM2(src_grd_file)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC4')

for file in lst_file:
# remap
    zeta = remap_clm(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
    dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC4', zeta=zeta)
    remap_clm(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
    remap_clm(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
    remap_clm_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

# merge file
    clim_file = dst_dir + file.rsplit('/')[-1][:-3] + '_clim_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_ssh_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_temp_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_salt_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
