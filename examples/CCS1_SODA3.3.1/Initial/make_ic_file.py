import subprocess
import os
import sys
import subprocess
import numpy as np
import pyroms
import pyroms_toolbox

from remap import remap
from remap_uv import remap_uv

#date to remap
data_dir = '/import/archive/u1/uaf/AKWATERS/kshedstrom/SODA/'
dst_dir='./'

tag='2010_12_26'
filein=data_dir + 'soda3.3.1_5dy_ocean_reg_' + tag + '.nc'

# load grids
src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir + 'SODA3_0.5deg_grid.nc', name='SODA3.3.1', xrange=(400, 500), yrange=(180, 280))
dst_grd = pyroms.grid.get_ROMS_grid('CCS_10')

print('\nBuild IC file from %s' %filein)

zeta_dst_file = dst_dir + dst_grd.name + '_ic_zeta_' + tag + '_' + src_grd.name + '.nc'
temp_dst_file = dst_dir + dst_grd.name + '_ic_temp_' + tag + '_' + src_grd.name + '.nc'
salt_dst_file = dst_dir + dst_grd.name + '_ic_salt_' + tag + '_' + src_grd.name + '.nc'
u_dst_file    = dst_dir + dst_grd.name + '_ic_u_'    + tag + '_' + src_grd.name + '.nc'
v_dst_file    = dst_dir + dst_grd.name + '_ic_v_'    + tag + '_' + src_grd.name + '.nc'

# remap ssh
zeta = remap('ssh', filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

# reload grid with zeta (more accurate)
dst_grd = pyroms.grid.get_ROMS_grid('CCS_10', zeta=zeta)

# regrid temp, salt and velocities
remap('temp',filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
remap('salt',filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
remap_uv(filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

# merge file
ic_file = dst_dir + dst_grd.name + '_ic_' + tag + '_' + src_grd.name + '.nc'

command1 = 'mv '      + zeta_dst_file + ' '    + ic_file
command2 = 'ncks -A ' + temp_dst_file + ' -o ' + ic_file
command3 = 'ncks -A ' + salt_dst_file + ' -o ' + ic_file
command4 = 'ncks -A ' + u_dst_file    + ' -o ' + ic_file
command5 = 'ncks -A ' + v_dst_file    + ' -o ' + ic_file

subprocess.call(command1,shell=True)
subprocess.call(command2,shell=True)
subprocess.call(command3,shell=True)
subprocess.call(command4,shell=True)
subprocess.call(command5,shell=True)

# clean up
os.remove(temp_dst_file)
os.remove(salt_dst_file)
os.remove(u_dst_file)
os.remove(v_dst_file)
