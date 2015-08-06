import pyroms
import pyroms_toolbox

src_varname = ['zeta','uice','vice','aice','hice','tisrf','snow_thick', \
               'sfwat','ageice','ti','sig11','sig12','sig22','t0mk', \
               's0mk','temp','salt','u','v']

# Change src_filename to your directory for the file's containing variable data
src_filename = '/archive/u1/uaf/kate/NEP5/run42/weeks_2000/nep5_001.nc'
#src_filename = '/archive/u1/uaf/kate/NEP5/run42/weeks_2000/nep5_04[89].nc', '/archive/u1/uaf/kate/NEP5/run42/weeks_2000/nep5_05?.nc'
wts_file = "./remap_weights_NEP5_to_BERING_bilinear_*"
src_grd = pyroms.grid.get_ROMS_grid('NEP5')
dst_grd = pyroms.grid.get_ROMS_grid('BERING')
# Outfile is a parameter to allow you to place these created remap files in a different
# directory than the one that is default which is where the file came from.
dst_var = pyroms_toolbox.remapping(src_varname, src_filename,\
                                   wts_file,src_grd,dst_grd,rotate_uv=False,\
                                   irange=(20,210),jrange=(340,610))
