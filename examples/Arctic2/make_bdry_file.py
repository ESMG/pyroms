import matplotlib
matplotlib.use('Agg')
import subprocess
from multiprocessing import Pool
from functools import partial
import os
import sys
import numpy as np

import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

def do_file(file, src_grd, dst_grd):
    dst_dir='./SODA/'

    zeta = remap_bdry(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
    dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2', zeta=zeta)
    remap_bdry(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
    remap_bdry(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
    remap_bdry_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

    # merge file
    bdry_file = dst_dir + file.rsplit('/')[-1][:-4] + '_bdry_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_ssh_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, bdry_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_temp_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_salt_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_u_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-4] + '_v_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
#    subprocess.check_call(command)
    subprocess.call(command)
    os.remove(out_file)

year = int(sys.argv[1])
lst_year = [year]

data_dir = '/nfs/P1/Data/SODA/SODA_2.1.6/'

lst_file = []

for year in lst_year:
    year = np.str(year)
#    command = 'ls ' + data_dir + 'SODA_2.1.6_' + year + '0[1-6]*'
#    command = 'ls ' + data_dir + 'SODA_2.1.6_' + year + '1*'
#    command = 'ls ' + data_dir + 'SODA_2.1.6_' + year + '0628*'
    command = 'ls ' + data_dir + 'SODA_2.1.6_' + year + '*'
    lst = subprocess.check_output(command, shell=True)
    lst = lst.split()
    lst_file = lst_file + lst

print('Build OBC file from the following file list:')
print(lst_file)
print(' ')

src_grd_file = data_dir + 'SODA_grid.cdf'
src_grd = pyroms_toolbox.BGrid_SODA.get_nc_BGrid_SODA(src_grd_file, name='SODA_2.1.6_ARCTIC2', area='npolar', ystart=240)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC2')

processes = 4
p = Pool(processes)
# Trick to pass more than one arg
partial_do_file = partial(do_file, src_grd=src_grd, dst_grd=dst_grd)
results = p.map(partial_do_file, lst_file)
