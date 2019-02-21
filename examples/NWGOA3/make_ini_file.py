import matplotlib
matplotlib.use('Agg')

import pyroms
import pyroms_toolbox

#src_varname = ['u', 'v']
src_varname = ['zeta', 'temp', 'salt', 'u', 'v', 'bac', 'c1', 'c2', 'c3', \
               'chl1', 'chl2', 'chl3', 'cldoc', 'csdoc', 'ddc', 'ddca', \
               'ddn', 'ddsi', 'ldoc', 'ldon', 'nh4', 'no3', 'ox', 'po4', \
               's1', 's2', 's3', 'sdoc', 'sdon', 'sio4', 'talk', 'tco2', \
               'zz1', 'zz2', 'zzc1', 'zzc2']
print('Number of variables', len(src_varname))
irange = None
jrange = None

# Need to pick up ubar/vbar and uice/vice later

# Change src_filename to your directory for the file's containing variable data
src_filename = '/archive/u1/uaf/kate/COSINE/data/day_000.nc'
wts_file = "./remap_weights_PP_COSINE_to_NWGOA3_bilinear_*"
src_grd = pyroms.grid.get_ROMS_grid('PP_COSINE')
dst_grd = pyroms.grid.get_ROMS_grid('NWGOA3')
# Outfile is a parameter to allow you to place these created remap files in a different
# directory than the one that is default which is where the file came from.
dst_var = pyroms_toolbox.remapping(src_varname, src_filename,\
                wts_file,src_grd,dst_grd,rotate_uv=True)
