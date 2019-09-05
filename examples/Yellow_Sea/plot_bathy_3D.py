import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
import matplotlib.collections as collections
from mpl_toolkits.mplot3d.art3d import Line3DCollection, line_collection_2d_to_3d

import pyroms
import pyroms_toolbox

grd = pyroms.grid.get_ROMS_grid('YELLOW')


lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
h = grd.vgrid.h
mask = grd.hgrid.mask_rho
idx = np.where(mask == 0)
h[idx] = 0


fig = plt.figure(figsize=(10, 5))
ax = fig.gca(projection='3d', azim=-50, elev=70)
pal = cm.Spectral
pal.set_over('#666666', 1.0)
vc = np.arange(-80,-5+0.1,0.1)
pal_norm = colors.BoundaryNorm(vc,ncolors=256, clip = False)
surf = ax.plot_surface(lon, lat, -h, rstride=1, cstride=1, cmap=pal, norm=pal_norm, linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.9, aspect=40)

coast = pyroms_toolbox.get_coast_line(grd)
coast = np.array(coast)
coast = collections.LineCollection(coast)
line_collection_2d_to_3d(coast, zs=-1, zdir='z')
coast.set_color('k')
ax.add_collection3d(coast, zs=-1, zdir='z')

outfile='YELLOW_bathy_3D.png'
plt.savefig(outfile, dpi=300, facecolor='w', edgecolor='w', \
                 orientation='portrait')

