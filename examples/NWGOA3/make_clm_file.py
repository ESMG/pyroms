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

irange = None
jrange = None

def do_file(file):
    print('file is: ' + file)
    start = string.find(file,'_')
    end = string.find(file,'_',start+1)
    var = file[start+1:end]
    print('var is: ' + var)
    if var == 'uv':
        dst_var = pyroms_toolbox.remapping(['u', 'v'], file,\
                     wts_file,src_grd,dst_grd,rotate_uv=True,\
                     uvar='u', vvar='v')
    else:
        dst_var = pyroms_toolbox.remapping([var], file,\
                     wts_file,src_grd,dst_grd)

lst_file = []

# Change src_filename to your directory for the file's containing variable data
data_dir = '/archive/u1/uaf/kate/COSINE/months/'
lst = subprocess.getoutput('ls ' + data_dir + 'Pac*_2008.nc')
#lst = commands.getoutput('ls ' + data_dir + 'Pacific_z*2_1987.nc')
#lst = commands.getoutput('ls ' + data_dir + 'Pacific_[tuvz]*_2008.nc')
lst = lst.split()
lst_file = lst_file + lst

wts_file = "./remap_weights_PP_COSINE_to_NWGOA3_bilinear_*"
src_grd = pyroms.grid.get_ROMS_grid('PP_COSINE')
dst_grd = pyroms.grid.get_ROMS_grid('NWGOA3')
# Outfile is a parameter to allow you to place these created remap files in a different
# directory than the one that is default which is where the file came from.
#processes = 1
processes = 4
p = Pool(processes)
results = p.map(do_file, lst_file)

