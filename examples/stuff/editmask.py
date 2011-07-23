import pyroms
import pyroms_toolbox
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4

grd = pyroms.grid.get_ROMS_grid('CHUKCHI')

# dir(grd), dir(grd.hgrid)
#
# Resolution values are 'c', 'l', 'i', 'h', and 'f'.
#m = Basemap(projection='npstere', boundinglat=60, lon_0=205, resolution='h')
# S_Africa
#m = Basemap(projection='lcc',lat_1=-35, lat_2=-10, lon_0=25, lat_0=-20, width=7000000, height=7000000, resolution='h')
# Benguela
#m = Basemap(projection='lcc',lat_1=-35, lat_2=-10, lon_0=10, lat_0=-13, width=4000000, height=5000000, resolution='h')
#m = Basemap(projection='lcc',lat_1=30, lat_2=40, lon_0=-120, lat_0=35, width=2000000, height=5000000, resolution='h')

coast = pyroms.utility.get_coast_from_map(m)
pyroms.grid.edit_mask_mesh_ij(grd.hgrid, coast=coast)


#pyroms.grid.edit_mask_mesh(grd.hgrid, proj=m)

pyroms.grid.write_ROMS_grid(grd, filename='grid_py.nc')
# ncks -v mask_rho,mask_u,mask_v,mask_psi grid_py.nc grid_Chukchi.nc


