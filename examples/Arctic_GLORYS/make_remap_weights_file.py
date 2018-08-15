import matplotlib
matplotlib.use('Agg')
import pyroms
import pyroms_toolbox
#import CGrid_GLORYS

# load the grid
#srcgrd = CGrid_GLORYS.get_nc_CGrid_GLORYS('/Volumes/P1/Data/GLORYS/data/GL2V1_mesh_mask_new.nc', name='GLORYS_ARCTIC2', area='npolar', ystart=690)
#srcgrd = pyroms_toolbox.CGrid_GLORYS.get_nc_CGrid_GLORYS('/Volumes/P1/Data/GLORYS/data/GL2V1_mesh_mask_new.nc', name='GLORYS_ARCTIC2', area='npolar', ystart=690)
srcgrd = pyroms_toolbox.CGrid_GLORYS.get_nc_CGrid_GLORYS('/archive/u1/uaf/kate/GLORYS/GL2V1_mesh_mask_new.nc', name='GLORYS_ARCTIC2', area='npolar', ystart=690)
dstgrd = pyroms.grid.get_ROMS_grid('ARCTIC2')

# make remap grid file for scrip
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='t')
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='u')
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='v')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_GLORYS_ARCTIC2_t.nc'
grid2_file = 'remap_grid_ARCTIC2_rho.nc'
interp_file1 = 'remap_weights_GLORYS_to_ARCTIC2_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_ARCTIC2_to_GLORYS_bilinear_rho_to_t.nc'
map1_name = 'GLORYS to ARCTIC2 Bilinear Mapping'
map2_name = 'ARCTIC2 to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.true.', grid2_periodic='.true.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_GLORYS_ARCTIC2_u.nc'
grid2_file = 'remap_grid_ARCTIC2_rho.nc'
interp_file1 = 'remap_weights_GLORYS_to_ARCTIC2_bilinear_u_to_rho.nc'
interp_file2 = 'remap_weights_ARCTIC2_to_GLORYS_bilinear_rho_to_u.nc'
map1_name = 'GLORYS to ARCTIC2 Bilinear Mapping'
map2_name = 'ARCTIC2 to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.true.', grid2_periodic='.true.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_GLORYS_ARCTIC2_v.nc'
grid2_file = 'remap_grid_ARCTIC2_rho.nc'
interp_file1 = 'remap_weights_GLORYS_to_ARCTIC2_bilinear_v_to_rho.nc'
interp_file2 = 'remap_weights_ARCTIC2_to_GLORYS_bilinear_rho_to_v.nc'
map1_name = 'GLORYS to ARCTIC2 Bilinear Mapping'
map2_name = 'ARCTIC2 to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.true.', grid2_periodic='.true.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_GLORYS_ARCTIC2_t.nc'
grid2_file = 'remap_grid_ARCTIC2_u.nc'
interp_file1 = 'remap_weights_GLORYS_to_ARCTIC2_bilinear_t_to_u.nc'
interp_file2 = 'remap_weights_ARCTIC2_to_GLORYS_bilinear_u_to_t.nc'
map1_name = 'GLORYS to ARCTIC2 Bilinear Mapping'
map2_name = 'ARCTIC2 to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.true.', grid2_periodic='.true.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_GLORYS_ARCTIC2_t.nc'
grid2_file = 'remap_grid_ARCTIC2_v.nc'
interp_file1 = 'remap_weights_GLORYS_to_ARCTIC2_bilinear_t_to_v.nc'
interp_file2 = 'remap_weights_ARCTIC2_to_GLORYS_bilinear_v_to_t.nc'
map1_name = 'GLORYS to ARCTIC2 Bilinear Mapping'
map2_name = 'ARCTIC2 to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.true.', grid2_periodic='.true.')

