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


file = '/archive/u1/uaf/kate/HYCOM/SCS/data/HYCOM_GLBa0.08_2015_001.nc'
data_dir = '/archive/u1/uaf/kate/HYCOM/SCS/data/'
dst_dir='./'

print('Build IC file from the following file:')
print(file)
print(' ')

src_grd_file = data_dir + '../HYCOM_GLBa0.08_PALAU_grid.nc'
src_grd = pyroms_toolbox.Grid_HYCOM.get_nc_Grid_HYCOM(src_grd_file)
dst_grd = pyroms.grid.get_ROMS_grid('PALAU1')

# remap
zeta = remap(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
dst_grd = pyroms.grid.get_ROMS_grid('PALAU1', zeta=zeta)
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
