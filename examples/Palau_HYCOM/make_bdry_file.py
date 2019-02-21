import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
from multiprocessing import Pool
from functools import partial
import numpy as np
#import pdb

#increase the maximum number of open files allowed
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE, (3000,-1))

import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

data_dir = '/archive/u1/uaf/kate/HYCOM/SCS/data/'
dst_dir='./bdry/'

def do_file(file, src_grd, dst_grd):
    zeta = remap_bdry(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
    dst_grd = pyroms.grid.get_ROMS_grid('PALAU1', zeta=zeta)
    remap_bdry(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
    remap_bdry(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
#    pdb.set_trace()
    remap_bdry_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

    # merge file
    bdry_file = dst_dir + file.rsplit('/')[-1][:-3] + '_bdry_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_ssh_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_temp_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_salt_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)

year = int(sys.argv[1])
#lst_year = sys.argv[1:]
lst_year = [year]
lst_file = []

for year in lst_year:
    year = np.str(year)
    command = 'ls ' + data_dir + 'HYCOM_GLBa0.08_' + year + '*'
    lst = subprocess.check_output(command, shell=True)
    lst = lst.split()
    lst_file = lst_file + lst

print('Build OBC file from the following file list:')
print(lst_file)
print(' ')

src_grd_file = data_dir + '../HYCOM_GLBa0.08_PALAU_grid.nc'
src_grd = pyroms_toolbox.Grid_HYCOM.get_nc_Grid_HYCOM(src_grd_file)
dst_grd = pyroms.grid.get_ROMS_grid('PALAU1')

processes = 4
p = Pool(processes)
# Trick to pass more than one arg
partial_do_file = partial(do_file, src_grd=src_grd, dst_grd=dst_grd)
results = p.map(partial_do_file, lst_file)
