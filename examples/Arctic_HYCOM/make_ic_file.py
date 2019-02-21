import subprocess
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')

import pyroms
import pyroms_toolbox

from remap import remap
from remap_uv import remap_uv


file = '/archive/u1/uaf/kate/HYCOM/Svalbard/data/HYCOM_GLBa0.08_2009_001.nc'
dst_dir='./'

print('Build IC file from the following file:')
print(file)
print(' ')

src_grd = pyroms_toolbox.Grid_HYCOM.get_nc_Grid_HYCOM('/archive/u1/uaf/kate/HYCOM/Svalbard/HYCOM_GLBa0.08_North_grid2.nc')
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2')

# remap
zeta = remap(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2', zeta=zeta)
remap(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
remap(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
remap_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

# merge file
ic_file = dst_dir + file.rsplit('/')[-1][:-3] + '_ic_' + dst_grd.name + '.nc'

out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_ssh_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-O', out_file, ic_file) 
print(command)
subprocess.check_call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_temp_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
print(command)
subprocess.check_call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_salt_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
print(command)
subprocess.check_call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
print(command)
subprocess.check_call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
print(command)
subprocess.check_call(command)
os.remove(out_file)
