import matplotlib
matplotlib.use('Agg')
import pyroms
import pyroms_toolbox

src_varname = ['zeta','temp','salt','u_eastward','v_northward']

irange=(370,580)
jrange=(460,580)
#irange = None
#jrange = None

# Change src_filename to your directory for the file's containing variable data
src_filename = '/archive/u1/uaf/kate/Arctic2/run46/averages/arctic2_avg_1999-*.nc'

wts_file = "./remap_weights_ARCTIC2_to_BEAUFORT2_bilinear_*"
src_grd = pyroms.grid.get_ROMS_grid('ARCTIC2')
dst_grd = pyroms.grid.get_ROMS_grid('BEAUFORT2')
# Outfile is a parameter to allow you to place these created remap files in a different
# directory than the one that is default which is where the file came from.
dst_var = pyroms_toolbox.remapping_bound(src_varname, src_filename,\
                     wts_file,src_grd,dst_grd,rotate_uv=True,\
                     irange=irange,jrange=jrange, \
                     uvar='u_eastward', vvar='v_northward', rotate_part=True)
