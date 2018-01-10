import pyroms

srcgrd = pyroms.grid.get_ROMS_grid('NEP5')
dstgrd = pyroms.grid.get_ROMS_grid('BERING')

pyroms.remapping.make_remap_grid_file(srcgrd)
pyroms.remapping.make_remap_grid_file(srcgrd,Cpos='u')
pyroms.remapping.make_remap_grid_file(srcgrd,Cpos='v')

pyroms.remapping.make_remap_grid_file(dstgrd)
pyroms.remapping.make_remap_grid_file(dstgrd,Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd,Cpos='v')

type = ['rho','u','v']

for typ in type:
        grid1_file = 'remap_grid_NEP5_'+str(typ)+'.nc'
        grid2_file = 'remap_grid_BERING_'+str(typ)+'.nc'
        interp_file1 = 'remap_weights_NEP5_to_BERING_bilinear_'+str(typ)+'_to_'+str(typ)+'.nc'
        interp_file2 = 'remap_weights_BERING_to_NEP5_bilinear_'+str(typ)+'_to_'+str(typ)+'.nc'
        map1_name = 'NEP5 to BERING Bilinear Mapping'
        map2_name = 'BERING to NEP5 Bilinear Mapping'
        num_maps = 1
        map_method = 'bilinear'
            
        print("Making "+str(interp_file1)+"...")
            
        pyroms.remapping.compute_remap_weights(grid1_file,grid2_file,\
                         interp_file1,interp_file2,map1_name,\
                         map2_name,num_maps,map_method)
