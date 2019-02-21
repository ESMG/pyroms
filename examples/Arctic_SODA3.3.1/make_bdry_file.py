import matplotlib
matplotlib.use('Agg')
import subprocess
import os
import sys
import commands
import numpy as np

#increase the maximum number of open files allowed
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE, (3000,-1))

import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

my_year=int(sys.argv[-1])

data_dir = '/archive/u1/uaf/AKWATERS/kshedstrom/SODA/'
data_dir_year = data_dir + str(my_year) + '/'
dst_dir='bdry/'

filelst = subprocess.check_output(['ls', data_dir_year]).replace('/n',' ').split()

#src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir + 'grid/SODA3_0.5deg_grid.nc', name='SODA3.3.1', xrange=(285, 500), yrange=(180, 300) )
#src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir + 'grid/SODA3_0.5deg_grid.nc', name='SODA3.3.1', area=npolar, ystart=236)
src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir + 'SODA3_0.5deg_grid.nc', name='SODA3.3.1', area='npolar')
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC4')

for filein in filelst:
    tag=filein.replace('soda3.3.1_5dy_ocean_reg_','').replace('.nc','')
    print '\nBuild OBC file for time %s' %filein
    zeta_dst_file = dst_dir + dst_grd.name + '_bdry_zeta_' + tag + '_' + src_grd.name + '.nc'
    temp_dst_file = dst_dir + dst_grd.name + '_bdry_temp_' + tag + '_' + src_grd.name + '.nc'
    salt_dst_file = dst_dir + dst_grd.name + '_bdry_salt_' + tag + '_' + src_grd.name + '.nc'
    u_dst_file    = dst_dir + dst_grd.name + '_bdry_u_'    + tag + '_' + src_grd.name + '.nc'
    v_dst_file    = dst_dir + dst_grd.name + '_bdry_v_'    + tag + '_' + src_grd.name + '.nc'

    # remap ssh
    zeta = remap_bdry('ssh', data_dir_year + filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

    # reload grid with zeta (more accurate)
    dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC4', zeta=zeta)

    # regrid temp, salt and velocities
    remap_bdry('temp',data_dir_year + filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
    remap_bdry('salt',data_dir_year + filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
    remap_bdry_uv(data_dir_year + filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

    # merge file
    bdry_file = dst_dir + dst_grd.name + '_bdry_' + tag + '_' + src_grd.name + '.nc'

    command1 = 'mv '      + zeta_dst_file + ' '    + bdry_file
    command2 = 'ncks -A ' + temp_dst_file + ' -o ' + bdry_file
    command3 = 'ncks -A ' + salt_dst_file + ' -o ' + bdry_file
    command4 = 'ncks -A ' + u_dst_file    + ' -o ' + bdry_file
    command5 = 'ncks -A ' + v_dst_file    + ' -o ' + bdry_file

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
