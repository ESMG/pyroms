import matplotlib
matplotlib.use('Agg')
import subprocess
import os
import subprocess
import numpy as np

import pyroms
import pyroms_toolbox

from remap import remap
from remap_uv import remap_uv

#file = '/Volumes/R2/Data/SODA/SODA_2.1.6/SODA_2.1.6_20030714-20030719.cdf'
file = '/nfs/P1/Data/SODA/SODA_2.1.6/SODA_2.1.6_20031231-20040105.cdf'
dst_dir='./'

print('Build IC file from the following file:')
print(file)
print(' ')

src_grd = pyroms_toolbox.BGrid_SODA.get_nc_BGrid_SODA('/nfs/P1/Data/SODA/SODA_2.1.6/SODA_grid.cdf', name='SODA_2.1.6', area='npolar', ystart=240)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2')

# remap
zeta = remap(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2', zeta=zeta)
remap(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
remap(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
remap_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

# merge file
ic_file = dst_dir + file.rsplit('/')[-1][:-4] + '_ic_' + dst_grd.name + '.nc'

out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_ssh_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-O', out_file, ic_file) 
#subprocess.check_call(command)
subprocess.call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_temp_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
#subprocess.check_call(command)
subprocess.call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_salt_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
#subprocess.check_call(command)
subprocess.call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_u_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
#subprocess.check_call(command)
subprocess.call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_v_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
#subprocess.check_call(command)
subprocess.call(command)
os.remove(out_file)
