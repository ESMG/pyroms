import pyroms

# Part of NEP grid containing the Bering
irange=(55,110)
jrange=(230,283)
irange=None
jrange=None

srcgrd = pyroms.grid.get_ROMS_grid('ARC_NATL')
dstgrd = pyroms.grid.get_ROMS_grid('CHUKCHI')

pyroms.remapping.make_remap_grid_file(srcgrd,irange=irange,jrange=jrange)
pyroms.remapping.make_remap_grid_file(srcgrd,Cpos='u',irange=irange,jrange=jrange)
pyroms.remapping.make_remap_grid_file(srcgrd,Cpos='v',irange=irange,jrange=jrange)

pyroms.remapping.make_remap_grid_file(dstgrd)
pyroms.remapping.make_remap_grid_file(dstgrd,Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd,Cpos='v')

type = ['rho','u','v']

for typ in type:
    for tip in type:
        grid1_file = 'remap_grid_ARC_NATL_'+str(typ)+'.nc'
        grid2_file = 'remap_grid_CHUKCHI_'+str(tip)+'.nc'
        interp_file1 = 'remap_weights_ARC_NATL_to_CHUKCHI_bilinear_'+str(typ)+'_to_'+str(tip)+'.nc'
        interp_file2 = 'remap_weights_CHUKCHI_to_ARC_NATL_bilinear_'+str(tip)+'_to_'+str(typ)+'.nc'
        map1_name = 'ARC_NATL to CHUKCHI Bilinear Mapping'
        map2_name = 'CHUKCHI to ARC_NATL Bilinear Mapping'
        num_maps = 1
        map_method = 'bilinear'
            
        print "Making "+str(interp_file1)+"..."
            
        pyroms.remapping.compute_remap_weights(grid1_file,grid2_file,\
                         interp_file1,interp_file2,map1_name,\
                         map2_name,num_maps,map_method)
