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

file = '/archive/u1/uaf/kate/GLORYS/GLORYS2V3_1dAV_19980101_19980102.nc'

dst_dir='./'

print('Build IC file from the following file:')
print(file)
print(' ')

src_grd = pyroms_toolbox.CGrid_GLORYS.get_nc_CGrid_GLORYS('/archive/u1/uaf/kate/GLORYS/GL2V1_mesh_mask_new.nc', name='GLORYS', area='npolar', ystart=690)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2')

# remap
zeta = remap(file, 'sossheig', src_grd, dst_grd, dst_dir=dst_dir)
remap(file, 'iicethic', src_grd, dst_grd, dst_dir=dst_dir)
remap(file, 'ileadfra', src_grd, dst_grd, dst_dir=dst_dir)
remap(file, 'iicetemp', src_grd, dst_grd, dst_dir=dst_dir)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2', zeta=zeta)
remap(file, 'votemper', src_grd, dst_grd, dst_dir=dst_dir)
remap(file, 'vosaline', src_grd, dst_grd, dst_dir=dst_dir)
#remap_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

# merge file
ic_file = dst_dir + file.rsplit('/')[-1][:-4] + '_ic_' + dst_grd.name + '.nc'

out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_sossheig_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-O', out_file, ic_file) 
subprocess.call(command)
os.remove(out_file)

out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_iicethic_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.call(command)
os.remove(out_file)

out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_ileadfra_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.call(command)
os.remove(out_file)

out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_iicetemp_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
subprocess.call(command)
os.remove(out_file)

out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_votemper_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file) 
#subprocess.check_call(command)
subprocess.call(command)
os.remove(out_file)
out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_vosaline_ic_' + dst_grd.name + '.nc'
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
