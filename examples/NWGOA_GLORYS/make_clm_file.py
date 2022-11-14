import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
import numpy as np
from multiprocessing import Pool

import pyroms
import pyroms_toolbox

from remap_clm import remap_clm
from remap_clm_uv import remap_clm_uv

lst_year = sys.argv[1:]

data_dir = '/import/AKWATERS/kshedstrom/glorys/monthly/'
dst_dir='./clm/'

lst_file = []

for year in lst_year:
    year = np.str(year)
    lst = subprocess.getoutput('ls ' + data_dir + 'GLORYS_REANALYSIS_*' + year + '*.nc')
#   lst = subprocess.getoutput('ls ' + data_dir + 'GLORYS_REANALYSIS_*' + year + '-12*.nc')
    lst = lst.split()
    lst_file = lst_file + lst

print('Build CLM file from the following file list:')
print(lst_file)
print(' ')

irange = (230, 460)
jrange = (150, 328)

src_grd_file = data_dir + 'GLORYS_REANALYSIS_1993-01-avg.nc'
src_grd = pyroms_toolbox.Grid_GLORYS.get_nc_Grid_GLORYS(src_grd_file, irange=irange, jrange=jrange)

def do_file(file):
# remap
    dst_grd = pyroms.grid.get_ROMS_grid('NWGOA')
    zeta = remap_clm(file, 'zos', src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)
    dst_grd = pyroms.grid.get_ROMS_grid('NWGOA', zeta=zeta)
    remap_clm(file, 'thetao', src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)
    remap_clm(file, 'so', src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)
    remap_clm_uv(file, src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)

# merge file
    clim_file = dst_dir + file.rsplit('/')[-1][:-3] + '_clim_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_zos_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_thetao_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_so_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)

processes = 4
p = Pool(processes)
results = p.map(do_file, lst_file)
