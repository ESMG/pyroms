import pyroms
import pyroms_toolbox

# load the grid
srcgrd = pyroms_toolbox.BGrid_SODA.get_nc_BGrid_SODA('/Volumes/r2/data/SODA/SODA_2.1.6/SODA_grid.cdf', name='SODA_2.1.6_ARCTIC2', area='npolar', ystart=240)
dstgrd = pyroms.grid.get_ROMS_grid('ARCTIC2')

# make remap grid file for scrip
pyroms_toolbox.BGrid_SODA.make_remap_grid_file(srcgrd, Bpos='t')
pyroms_toolbox.BGrid_SODA.make_remap_grid_file(srcgrd, Bpos='uv')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_SODA_2.1.6_ARCTIC2_t.nc'
grid2_file = 'remap_grid_ARCTIC2_rho.nc'
interp_file1 = 'remap_weights_SODA_2.1.6_to_ARCTIC2_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_ARCTIC2_to_SODA_2.1.6_bilinear_rho_to_t.nc'
map1_name = 'SODA_2.1.6 to ARCTIC2 Bilinear Mapping'
map2_name = 'ARCTIC2 to SODA_2.1.6 Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.true.', grid2_periodic='.true.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_SODA_2.1.6_ARCTIC2_uv.nc'
grid2_file = 'remap_grid_ARCTIC2_rho.nc'
interp_file1 = 'remap_weights_SODA_2.1.6_to_ARCTIC2_bilinear_uv_to_rho.nc'
interp_file2 = 'remap_weights_ARCTIC2_to_SODA_2.1.6_bilinear_rho_to_uv.nc'
map1_name = 'SODA_2.1.6 to ARCTIC2 Bilinear Mapping'
map2_name = 'ARCTIC2 to SODA_2.1.6 Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.true.', grid2_periodic='.true.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_SODA_2.1.6_ARCTIC2_t.nc'
grid2_file = 'remap_grid_ARCTIC2_u.nc'
interp_file1 = 'remap_weights_SODA_2.1.6_to_ARCTIC2_bilinear_t_to_u.nc'
interp_file2 = 'remap_weights_ARCTIC2_to_SODA_2.1.6_bilinear_u_to_t.nc'
map1_name = 'SODA_2.1.6 to ARCTIC2 Bilinear Mapping'
map2_name = 'ARCTIC2 to SODA_2.1.6 Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.true.', grid2_periodic='.true.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_SODA_2.1.6_ARCTIC2_t.nc'
grid2_file = 'remap_grid_ARCTIC2_v.nc'
interp_file1 = 'remap_weights_SODA_2.1.6_to_ARCTIC2_bilinear_t_to_v.nc'
interp_file2 = 'remap_weights_ARCTIC2_to_SODA_2.1.6_bilinear_v_to_t.nc'
map1_name = 'SODA_2.1.6 to ARCTIC2 Bilinear Mapping'
map2_name = 'ARCTIC2 to SODA_2.1.6 Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.true.', grid2_periodic='.true.')

