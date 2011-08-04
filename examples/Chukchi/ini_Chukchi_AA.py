import pyroms
import pyroms_toolbox

src_varname = ['zeta','aice','hice','tisrf','snow_thick','sfwat', \
               'ageice','ti','sig11','sig12','sig22','t0mk','s0mk', \
               'temp','salt','u','v']
#irange=(55,116)
#jrange=(230,283)
irange = None
jrange = None

# Need to pick up ubar/vbar and uice/vice later

# Change src_filename to your directory for the file's containing variable data
src_filename = '/archive/u1/uaf/kate/ARC_NATL/run01/weeks_1985/arc_natl_1985_001.nc'
#src_filename = '/archive/u1/uaf/kate/NEP5/run42/weeks_2000/nep5_04[89].nc', '/archive/u1/uaf/kate/NEP5/run42/weeks_2000/nep5_05?.nc'
wts_file = "./remap_weights_ARC_NATL_to_CHUKCHI_bilinear_*"
src_grd = pyroms.grid.get_ROMS_grid('ARC_NATL')
dst_grd = pyroms.grid.get_ROMS_grid('CHUKCHI')
# Outfile is a parameter to allow you to place these created remap files in a different
# directory than the one that is default which is where the file came from.
dst_var = pyroms_toolbox.remapping(src_varname, src_filename,\
                wts_file,src_grd,dst_grd,rotate_uv=True,\
                trange=(0,0),irange=irange,jrange=jrange)
dst_var = pyroms_toolbox.remapping(['uice','vice'], src_filename,\
                wts_file,src_grd,dst_grd,rotate_uv=True,\
                trange=(0,0),irange=irange, \
		jrange=jrange,uvar='uice',vvar='vice')
