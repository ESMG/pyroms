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

from remap import remap
from remap_uv import remap_uv

my_year=int(sys.argv[-1])

data_dir = '/archive/u1/uaf/AKWATERS/kshedstrom/SODA/'
data_dir_year = data_dir + '/monthly/*' + str(my_year) + '*'
dst_dir='clm/'

filelst = commands.getoutput('ls ' + data_dir_year)
filelst = filelst.split()
#filelst = subprocess.check_output(['ls', data_dir_year]).replace('/n',' ').split()

#src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir + 'grid/SODA3_0.5deg_grid.nc', name='SODA3.3.1', xrange=(285, 500), yrange=(180, 300) )
#src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir + 'grid/SODA3_0.5deg_grid.nc', name='SODA3.3.1', area=npolar, ystart=236)
src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir + 'SODA3_0.5deg_grid.nc', name='SODA3.3.1', area='npolar')
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC4')
print src_grd.name
print dst_grd.name

for filein in filelst:
    tag=filein.replace(data_dir+'/monthly/','').replace('soda3.3.1_monthly_ocean_reg_','').replace('.nc','')
    print '\nBuild OBC file for time %s' %tag
    zeta_dst_file = dst_dir + dst_grd.name + '_clm_zeta_' + tag + '_' + src_grd.name + '.nc'
    temp_dst_file = dst_dir + dst_grd.name + '_clm_temp_' + tag + '_' + src_grd.name + '.nc'
    salt_dst_file = dst_dir + dst_grd.name + '_clm_salt_' + tag + '_' + src_grd.name + '.nc'
    u_dst_file    = dst_dir + dst_grd.name + '_clm_u_'    + tag + '_' + src_grd.name + '.nc'
    v_dst_file    = dst_dir + dst_grd.name + '_clm_v_'    + tag + '_' + src_grd.name + '.nc'

    # remap ssh
    zeta = remap('ssh', filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

    # reload grid with zeta (more accurate)
    dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC4', zeta=zeta)

    # regrid temp, salt and velocities
    remap('temp', filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
    remap('salt', filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
    remap_uv(filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

    # merge file
    clm_file = dst_dir + dst_grd.name + '_clm_' + tag + '_' + src_grd.name + '.nc'

    command1 = 'mv '      + zeta_dst_file + ' '    + clm_file
    command2 = 'ncks -A ' + temp_dst_file + ' -o ' + clm_file
    command3 = 'ncks -A ' + salt_dst_file + ' -o ' + clm_file
    command4 = 'ncks -A ' + u_dst_file    + ' -o ' + clm_file
    command5 = 'ncks -A ' + v_dst_file    + ' -o ' + clm_file

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
