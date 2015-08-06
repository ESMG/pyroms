import pyroms
import pyroms_toolbox

src_varname = ['zeta','temp','salt','u','v','ubar','vbar','uice', \
               'vice','aice','hice','tisrf','snow_thick','sfwat', \
               'ageice','ti','sig11','sig12','sig22','t0mk','s0mk']

# Change src_filename to your directory for the file's containing variable data
src_filename = '/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_00?.nc'
#src_filename = ['/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_001.nc']
#src_filename = [ '/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_046.nc', \
#        '/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_047.nc', \
#        '/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_048.nc', \
#        '/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_049.nc', \
#        '/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_050.nc', \
#        '/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_051.nc', \
#        '/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_052.nc', \
#        '/archive/u1/uaf/kate/NEP5/run42/weeks_2001/nep5_053.nc']

wts_file = "./remap_weights_NEP5_to_BERING_bilinear_*"
src_grd = pyroms.grid.get_ROMS_grid('NEP5')
dst_grd = pyroms.grid.get_ROMS_grid('BERING')
# Outfile is a parameter to allow you to place these created remap files in a different
# directory than the one that is default which is where the file came from.
dst_var = pyroms_toolbox.remapping_bound(src_varname, src_filename,\
                                   wts_file,src_grd,dst_grd,rotate_uv=False,\
                                   irange=(20,210),jrange=(340,610))
