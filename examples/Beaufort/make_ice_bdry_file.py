import matplotlib
matplotlib.use('Agg')
import subprocess
from multiprocessing import Pool

import pyroms
import pyroms_toolbox

src_varname = ['aice','hice','tisrf','snow_thick', \
               'ti','uice', 'vice']
src_sigma = ['sig11', 'sig22', 'sig12']

irange=(420,580)
jrange=(470,570)
#irange = None
#jrange = None

#months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
months = ['02', '03', '04', '05']

# Change src_filename to your directory for the file's containing variable data
#src_filename = '/archive/u1/uaf/kate/Arctic2/run24/averages2/arctic2_avg2_1992-12-3*.nc'
part_filename = '/archive/u1/uaf/kate/Arctic2/run24/averages2/arctic2_avg2_1993-'

wts_file = "./remap_weights_ARCTIC2_to_BEAUFORT_bilinear_*"
src_grd = pyroms.grid.get_ROMS_grid('ARCTIC2')
dst_grd = pyroms.grid.get_ROMS_grid('BEAUFORT')
# Outfile is a parameter to allow you to place these created remap files in a different
# directory than the one that is default which is where the file came from.

def do_file(month):
    src_filename = part_filename + month + '*.nc'
    lcopy = list(src_varname)
    dst_var = pyroms_toolbox.remapping_bound(lcopy, src_filename,\
                     wts_file,src_grd,dst_grd,rotate_uv=True,\
                     irange=irange,jrange=jrange, \
                     uvar='uice', vvar='vice')
#                     uvar='uice_eastward', vvar='vice_northward', rotate_part=True)
    dst_var = pyroms_toolbox.remapping_bound_sig(src_sigma, src_filename,\
                     wts_file,src_grd,dst_grd,rotate_sig=True,\
                     irange=irange,jrange=jrange)

processes = 4
p = Pool(processes)
results = p.map(do_file, months)
