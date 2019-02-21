import subprocess
import os
import sys
import subprocess
import numpy as np

import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

year = int(sys.argv[1])
lst_year = [year]

data_dir = '/Volumes/R1/Data/SODA_2.1.6/'
dst_dir='./'

lst_file = []

for year in lst_year:
    year = np.str(year)
    lst = subprocess.getoutput('ls ' + data_dir + 'SODA_2.1.6_' + year + '*')
    lst = lst.split()
    lst_file = lst_file + lst

print('Build OBC file from the following file list:')
print(lst_file)
print(' ')

src_grd_file = data_dir + 'SODA_grid.cdf'
src_grd = pyroms_toolbox.BGrid_SODA.get_nc_BGrid_SODA('/Volumes/R1/DATA/SODA_2.1.6/SODA_grid.cdf', name='SODA_2.1.6_YELLOW', xrange=(225, 275), yrange=(190, 240))
dst_grd = pyroms.grid.get_ROMS_grid('YELLOW')

for file in lst_file:
    zeta = remap_bdry(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
    dst_grd = pyroms.grid.get_ROMS_grid('YELLOW', zeta=zeta)
    remap_bdry(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
    remap_bdry(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
    remap_bdry_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

    # merge file
    bdry_file = dst_dir + file.rsplit('/')[-1][:-4] + '_bdry_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_ssh_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, bdry_file) 
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_temp_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_salt_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_u_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_v_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.check_call(command)
    os.remove(out_file)
