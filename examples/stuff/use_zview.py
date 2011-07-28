import pyroms
import pyroms_toolbox

grd = pyroms.grid.get_ROMS_grid('NEP5')

# In order to make an horizontal plot of a 2-D spacial field (like
# zeta) you need to use pyroms_toolbox.twoDview. To plot an horizontal
# plot of a 3-D spacial field (like temp) you can use zview (for a
# cte-z slice) or sview to plot sigma layer. If you plot from a netcdf
# file, you always have a time dependency, but if you past an array to
# the function without time dependency, then put -1 for the time index.

# For the grip projection, you can use pyroms.utility.get_grid_proj to
# get a Basemap object with all the projection parameters from a grid
# object. You can specify the type of projection using the optional
# type argument, mercator being the default projection.
# Use map = pyroms.utility.get_grid_proj(grd, type='lcc') for lambert
# conformal conic or
# map = pyroms.utility.get_grid_proj(grd, type='stere') for
# stereographic.
# For now I have almost always use the default mercator projection so
# you may have to modify a low bit the function as the parameter are
# not the same for all the projection type but this should be
# straightforward. You can also change the projection type in zview
# and sview using the optional argument proj. Again I have always use
# the default mercator projection so you may have to tune the function
# if you use a different projection. 

# plot temp at 5 metre depth for the first time idx in nep5_avg_00001.nc
pyroms_toolbox.zview('temp', 0, 5, grd, filename='nep5_avg_00001.nc')

# plot vertical section along i=20
# you can past an array instead of reading from a file
temp = pyroms.utility.get_nc_var('temp', 'nep5_avg_00001.nc')
pyroms_toolbox.iview(temp[:], 0, 20, grd)
