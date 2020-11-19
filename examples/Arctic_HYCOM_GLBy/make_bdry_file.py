import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
import numpy as np
from multiprocessing import Pool
#import pdb

#increase the maximum number of open files allowed
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE, (3000,-1))

import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

lst_year = sys.argv[1:]

data_dir = '/import/AKWATERS/kshedstrom/HYCOM/Svalbard/data/'
dst_dir='./bdry/'

lst_file = []

for year in lst_year:
    lst = subprocess.getoutput('ls ' + data_dir + 'HYCOM_GLBy0.08_' + year + '_[01]*')
#   lst = subprocess.getoutput('ls ' + data_dir + 'HYCOM_GLBy0.08_' + year + '*')
#    lst = subprocess.getoutput('ls ' + data_dir + 'HYCOM_GLBy0.08_' + year + '_0*')
#    lst = subprocess.getoutput('ls ' + data_dir + 'HYCOM_GLBy0.08_' + year + '_0[4-9]*')
    lst = lst.split()
    lst_file = lst_file + lst

print('Build OBC file from the following file list:')
print(lst_file)
print(' ')

src_grd_file = data_dir + 'HYCOM_GLBy0.08_2018_345.nc'
src_grd = pyroms_toolbox.Grid_HYCOM.get_nc_Grid_HYCOM2(src_grd_file)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC4')

def do_file(file):
    zeta = remap_bdry(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
    dst_grd2 = pyroms.grid.get_ROMS_grid('ARCTIC4', zeta=zeta)
    remap_bdry(file, 'temp', src_grd, dst_grd2, dst_dir=dst_dir)
    remap_bdry(file, 'salt', src_grd, dst_grd2, dst_dir=dst_dir)
#   pdb.set_trace()
    remap_bdry_uv(file, src_grd, dst_grd2, dst_dir=dst_dir)

    # merge file
    bdry_file = dst_dir + file.rsplit('/')[-1][:-3] + '_bdry_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_ssh_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_temp_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_salt_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)

processes = 4
p = Pool(processes)
results = p.map(do_file, lst_file)
