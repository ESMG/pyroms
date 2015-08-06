import pyroms
import pyroms_toolbox

src_varname = ['zeta','aice','hice','tisrf','snow_thick', \
               'ageice','ti','t0mk','s0mk', 'temp','salt','u_eastward','v_northward']
irange=(420,580)
jrange=(470,570)
#irange = None
#jrange = None

# Need to pick up ubar/vbar and uice/vice later

# Change src_filename to your directory for the file's containing variable data
src_filename = '/archive/u1/uaf/kate/Arctic2/run24/averages/arctic2_avg_1992-12-31T00:00:00.nc'
#src_filename = '/archive/u1/uaf/kate/NEP5/run42/weeks_2000/nep5_04[89].nc', '/archive/u1/uaf/kate/NEP5/run42/weeks_2000/nep5_05?.nc'
wts_file = "./remap_weights_ARCTIC2_to_BEAUFORT_bilinear_*"
src_grd = pyroms.grid.get_ROMS_grid('ARCTIC2')
dst_grd = pyroms.grid.get_ROMS_grid('BEAUFORT')
# Outfile is a parameter to allow you to place these created remap files in a different
# directory than the one that is default which is where the file came from.
dst_var = pyroms_toolbox.remapping(src_varname, src_filename, wts_file, \
                src_grd, dst_grd, rotate_uv=True, irange=irange, jrange=jrange, \
                uvar='u_eastward', vvar='v_northward', rotate_part=True)
#dst_var = pyroms_toolbox.remapping(['uice','vice'], src_filename, \
#                wts_file, src_grd, dst_grd, rotate_uv=True, \
#                irange=irange, jrange=jrange, uvar='uice', vvar='vice')
dst_var = pyroms_toolbox.remapping(['uice_eastward','vice_northward'], src_filename, \
                wts_file, src_grd, dst_grd, rotate_uv=True, \
                irange=irange, jrange=jrange, \
                uvar='uice_eastward', vvar='vice_northward', rotate_part=True)
dst_var = pyroms_toolbox.remapping_tensor(['sig11','sig22','sig12'], src_filename, \
                wts_file, src_grd, dst_grd, rotate_sig=True, \
                irange=irange, jrange=jrange, shapiro=False)
