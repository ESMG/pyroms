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

#my_year=int(sys.argv[-1])
my_year = '2010'

data_dir = '/archive/AKWATERS/kshedstrom/SODA/'
data_dir_year = data_dir + '/' + str(my_year) + '/soda3.3.1_5dy_ocean_reg_2010_12_26.nc'
dst_dir='./'

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
    tag=filein.replace(data_dir+'/'+my_year+'/','').replace('soda3.3.1_5dy_ocean_reg_','').replace('.nc','')
    print '\nBuild OBC file for time %s' %tag
    zeta_dst_file = dst_dir + dst_grd.name + '_IC_zeta_' + tag + '_' + src_grd.name + '.nc'
    temp_dst_file = dst_dir + dst_grd.name + '_IC_temp_' + tag + '_' + src_grd.name + '.nc'
    salt_dst_file = dst_dir + dst_grd.name + '_IC_salt_' + tag + '_' + src_grd.name + '.nc'
    u_dst_file    = dst_dir + dst_grd.name + '_IC_u_'    + tag + '_' + src_grd.name + '.nc'
    v_dst_file    = dst_dir + dst_grd.name + '_IC_v_'    + tag + '_' + src_grd.name + '.nc'

    # remap ssh
    zeta = remap('ssh', filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

    # reload grid with zeta (more accurate)
    dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC4', zeta=zeta)

    # regrid temp, salt and velocities
    remap('temp', filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
    remap('salt', filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
    remap_uv(filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

    # merge file
    IC_file = dst_dir + dst_grd.name + '_IC_' + tag + '_' + src_grd.name + '.nc'

    command1 = 'mv '      + zeta_dst_file + ' '    + IC_file
    command2 = 'ncks -A ' + temp_dst_file + ' -o ' + IC_file
    command3 = 'ncks -A ' + salt_dst_file + ' -o ' + IC_file
    command4 = 'ncks -A ' + u_dst_file    + ' -o ' + IC_file
    command5 = 'ncks -A ' + v_dst_file    + ' -o ' + IC_file

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
