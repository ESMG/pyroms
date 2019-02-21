import subprocess
import os
import subprocess
import _iso
import numpy as np

import pyroms
import pyroms_toolbox

from remap import remap
from remap_uv import remap_uv


file = '/center/w/kate/SODA/SODA_2.1.6_20071230-20080104.cdf'
dst_dir='./'

print('Build IC file from the following file:')
print(file)
print(' ')

src_grd = pyroms_toolbox.BGrid_SODA.get_nc_BGrid_SODA('/center/w/SODA/SODA_grid.cdf', name='SODA_2.1.6_YELLOW', xrange=(225, 275), yrange=(190, 240))
dst_grd = pyroms.grid.get_ROMS_grid('YELLOW')

# remap
zeta = remap(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
dst_grd = pyroms.grid.get_ROMS_grid('YELLOW', zeta=zeta)
remap(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
remap(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
remap_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

# merge file
ic_file = dst_dir + file.rsplit('/')[-1][:-4] + '_ic_' + dst_grd.name + '.nc'

out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_ssh_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-O', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_temp_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_salt_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_u_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_v_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.check_call(command)
os.remove(out_file)
