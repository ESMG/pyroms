import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
import numpy as np
import pdb


import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

lst_year = sys.argv[1:]

data_dir = '/import/AKWATERS/kshedstrom/glorys/'
dst_dir='./bdry/'

lst_file = []

for year in lst_year:
    lst = subprocess.getoutput('ls ' + data_dir + 'GLORYS_REANALYSIS_' + year + '-01-*')
#   lst = subprocess.getoutput('ls ' + data_dir + 'GLORYS_REANALYSIS_' + year + '*')
#    lst = subprocess.getoutput('ls ' + data_dir + 'GLORYS_REANALYSIS_' + year + '-0*')
#    lst = subprocess.getoutput('ls ' + data_dir + 'GLORYS_REANALYSIS_' + year + '-0[4-9]*')
    lst = lst.split()
    lst_file = lst_file + lst

print('Build OBC file from the following file list:')
print(lst_file)
print(' ')

irange = (230, 460)
jrange = (150, 328)

src_grd_file = data_dir + 'GLORYS_REANALYSIS_2006-09-30.nc'
src_grd = pyroms_toolbox.Grid_GLORYS.get_nc_Grid_GLORYS(src_grd_file, irange=irange, jrange=jrange)
dst_grd = pyroms.grid.get_ROMS_grid('NWGOA')

def do_file(file):
    pdb.set_trace()
    zeta = remap_bdry(file, 'zos', src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)
    dst_grd2 = pyroms.grid.get_ROMS_grid('NWGOA', zeta=zeta)
    remap_bdry_uv(file, src_grd, dst_grd2, dst_dir=dst_dir, irange=irange, jrange=jrange)
    remap_bdry(file, 'thetao', src_grd, dst_grd2, dst_dir=dst_dir, irange=irange, jrange=jrange)
    remap_bdry(file, 'so', src_grd, dst_grd2, dst_dir=dst_dir, irange=irange, jrange=jrange)

    # merge file
    bdry_file = dst_dir + file.rsplit('/')[-1][:-3] + '_bdry_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_zos_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_thetao_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_so_bdry_' + dst_grd.name + '.nc'
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

#processes = 4
#p = Pool(processes)
results = do_file(lst_file[1])
