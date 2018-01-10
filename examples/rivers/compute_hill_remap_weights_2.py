import numpy as np
from datetime import datetime
import netCDF4 as netCDF
import pdb

import pyroms
import pyroms_toolbox


##  create CI remap file for scrip
#print 'Create remap grid file for CI grid'
#dstgrd = pyroms.grid.get_ROMS_grid('CI')
#dstgrd.hgrid.mask_rho = np.ones(dstgrd.hgrid.mask_rho.shape)
#pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
#
#
### compute remap weights
print('compute remap weights using scrip')
# input namelist variables for conservative remapping at rho points
grid1_file = '../version1/remap_grid_runoff.nc'
grid2_file = 'remap_grid_CI_rho.nc'
interp_file1 = 'remap_weights_runoff_to_CI_conservative_nomask.nc'
interp_file2 = 'remap_weights_CI_to_runoff_conservative_nomask.nc'
map1_name = 'runoff to CI conservative Mapping'
map2_name = 'CI to runoff conservative Mapping'
num_maps = 1
map_method = 'conservative'

#pdb.set_trace()
pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)
