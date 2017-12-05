import pyroms
import pyroms_toolbox

# load the grid
srcgrd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL('/center/d/kate/COBALT/GFDL_CM2.1_grid.nc', name='ESM2M_NWGOA3')
dstgrd = pyroms.grid.get_ROMS_grid('NWGOA3')
dstgrd.hgrid.lon_rho = dstgrd.hgrid.lon_rho - 360.
dstgrd.hgrid.lon_u = dstgrd.hgrid.lon_u - 360.
dstgrd.hgrid.lon_v = dstgrd.hgrid.lon_v - 360.
dstgrd.hgrid.lon_psi = dstgrd.hgrid.lon_psi - 360.
dstgrd.hgrid.lon_vert = dstgrd.hgrid.lon_vert - 360.

# make remap grid file for scrip
pyroms_toolbox.BGrid_GFDL.make_remap_grid_file(srcgrd, Bpos='t')
pyroms_toolbox.BGrid_GFDL.make_remap_grid_file(srcgrd, Bpos='uv')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_ESM2M_NWGOA3_t.nc'
grid2_file = 'remap_grid_NWGOA3_rho.nc'
interp_file1 = 'remap_weights_ESM2M_to_NWGOA3_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_NWGOA3_to_ESM2M_bilinear_rho_to_t.nc'
map1_name = 'ESM2M to NWGOA3 Bilinear Mapping'
map2_name = 'NWGOA3 to ESM2M Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_ESM2M_NWGOA3_uv.nc'
grid2_file = 'remap_grid_NWGOA3_rho.nc'
interp_file1 = 'remap_weights_ESM2M_to_NWGOA3_bilinear_uv_to_rho.nc'
interp_file2 = 'remap_weights_NWGOA3_to_ESM2M_bilinear_rho_to_uv.nc'
map1_name = 'ESM2M to NWGOA3 Bilinear Mapping'
map2_name = 'NWGOA3 to ESM2M Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_ESM2M_NWGOA3_t.nc'
grid2_file = 'remap_grid_NWGOA3_u.nc'
interp_file1 = 'remap_weights_ESM2M_to_NWGOA3_bilinear_t_to_u.nc'
interp_file2 = 'remap_weights_NWGOA3_to_ESM2M_bilinear_u_to_t.nc'
map1_name = 'ESM2M to NWGOA3 Bilinear Mapping'
map2_name = 'NWGOA3 to ESM2M Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_ESM2M_NWGOA3_t.nc'
grid2_file = 'remap_grid_NWGOA3_v.nc'
interp_file1 = 'remap_weights_ESM2M_to_NWGOA3_bilinear_t_to_v.nc'
interp_file2 = 'remap_weights_NWGOA3_to_ESM2M_bilinear_v_to_t.nc'
map1_name = 'ESM2M to NWGOA3 Bilinear Mapping'
map2_name = 'NWGOA3 to ESM2M Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

