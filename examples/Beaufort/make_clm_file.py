import matplotlib
matplotlib.use('Agg')
import subprocess
from multiprocessing import Pool
from functools import partial
import os
import sys
import string

import pyroms
import pyroms_toolbox
import subprocess

irange=(370,580)
jrange=(460,580)

months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

src_varname = ['aice', 'hice', 'snow_thick', 'uice_eastward', 'vice_northward']
part_filename = '/archive/u1/uaf/kate/Arctic2/run46/averages/arctic2_avg_1998-'
src_grd = pyroms.grid.get_ROMS_grid('ARCTIC2')
dst_grd = pyroms.grid.get_ROMS_grid('BEAUFORT2')

def do_file(month):
    global src_varname
    wts_file = './remap_weights_ARCTIC2_to_BEAUFORT2_bilinear_*'
    src_filename = part_filename + month + '*.nc'
    lcopy = list(src_varname)
    dst_var = pyroms_toolbox.remapping(lcopy, src_filename, \
                     wts_file, src_grd, dst_grd, rotate_uv=True, \
                     irange=irange, jrange=jrange, \
                     uvar='uice_eastward', vvar='vice_northward', \
                     rotate_part=True)

#processes = 1
processes = 4
p = Pool(processes)
results = p.map(do_file, months)
