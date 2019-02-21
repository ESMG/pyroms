import pyroms

# Part of Arctic2 grid containing the Beaufort
irange=(370,580)
jrange=(460,580)
#irange=None
#jrange=None

srcgrd = pyroms.grid.get_ROMS_grid('ARCTIC2')
dstgrd = pyroms.grid.get_ROMS_grid('BEAUFORT2')

pyroms.remapping.make_remap_grid_file(srcgrd,irange=irange,jrange=jrange)
pyroms.remapping.make_remap_grid_file(srcgrd,Cpos='u',irange=irange,jrange=jrange)
pyroms.remapping.make_remap_grid_file(srcgrd,Cpos='v',irange=irange,jrange=jrange)

pyroms.remapping.make_remap_grid_file(dstgrd)
pyroms.remapping.make_remap_grid_file(dstgrd,Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd,Cpos='v')

type = ['rho','u','v']

for typ in type:
    for tip in type:
        grid1_file = 'remap_grid_ARCTIC2_'+str(typ)+'.nc'
        grid2_file = 'remap_grid_BEAUFORT2_'+str(tip)+'.nc'
        interp_file1 = 'remap_weights_ARCTIC2_to_BEAUFORT2_bilinear_'+str(typ)+'_to_'+str(tip)+'.nc'
        interp_file2 = 'remap_weights_BEAUFORT2_to_ARCTIC2_bilinear_'+str(tip)+'_to_'+str(typ)+'.nc'
        map1_name = 'ARCTIC2 to BEAUFORT Bilinear Mapping'
        map2_name = 'BEAUFORT to ARCTIC2 Bilinear Mapping'
        num_maps = 1
        map_method = 'bilinear'
            
        print("Making "+str(interp_file1)+"...")
            
        pyroms.remapping.compute_remap_weights(grid1_file,grid2_file,\
                         interp_file1,interp_file2,map1_name,\
                         map2_name,num_maps,map_method)
