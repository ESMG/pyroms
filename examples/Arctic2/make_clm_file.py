import matplotlib
matplotlib.use('Agg')
import subprocess
import os
import sys
import commands
import numpy as np
from datetime import datetime

import pyroms
import pyroms_toolbox

from remap_clm import remap_clm
from remap_clm_uv import remap_clm_uv

#increase the maximum number of open files allowed
import resource
resource.setrlimit(resource.RLIMIT_NOFILE, (3000,-1))

start = datetime.now()

#year = int(sys.argv[1])
#lst_year = [year]
lst_year = sys.argv[1:]

data_dir = '/nfs/P1/Data/SODA/SODA_2.1.6/Monthly_avg/'
#data_dir = '/Volumes/R2/Data/SODA/SODA_2.1.6/Monthly_avg/'
dst_dir='./SODA/'

lst_file = []

for year in lst_year:
    year = np.str(year)
#    lst = commands.getoutput('ls ' + data_dir + 'SODA_2.1.6_' + year + '_0*')
    lst = commands.getoutput('ls ' + data_dir + 'SODA_2.1.6_' + year + '*')
    lst = lst.split()
    lst_file = lst_file + lst

print 'Build CLM file from the following file list:'
print lst_file
print ' '

src_grd_file = data_dir + '../SODA_grid.cdf'
src_grd = pyroms_toolbox.BGrid_SODA.get_nc_BGrid_SODA(src_grd_file, name='SODA_2.1.6_ARCTIC2', area='npolar', ystart=240)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2')

for file in lst_file:
    zeta = remap_clm(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
    dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2', zeta=zeta)
    remap_clm(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
    remap_clm(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
    remap_clm_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

    # merge file
    clim_file = dst_dir + file.rsplit('/')[-1][:-3] + '_clim_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_ssh_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, clim_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_temp_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_salt_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)
