import numpy as np
import netCDF4
import os
import sys
import subprocess
#import Ngl
import pyroms
from pyroms_toolbox import jday2date
import pyroms_toolbox
import matplotlib.pyplot as plt

#year = int(sys.argv[1])
#lst_year = [year]

lst_file = []

#for year in lst_year:
#    year = np.str(year)
#lst = subprocess.getoutput('ls averages/*74??.nc')
lst = subprocess.getoutput('ls months/*.nc')
lst = lst.split()
lst_file = lst_file + lst

grd = pyroms.grid.get_ROMS_grid('ARCTIC2')

#clat = grd.hgrid.lat_rho
#clon = grd.hgrid.lon_rho

istart = 356
jstart = 373
iend = 387
jend = 411

for file in lst_file:
    print("Plotting "+file)
    nc = netCDF4.Dataset(file, "r")
    temp = nc.variables["temp"][0,:,:,:]
    time = nc.variables["ocean_time"][0]
    myday = jday2date(time/86400.)
#    date_tag = myday.strftime('%d %B %Y')
    date_tag = myday.strftime('%Y_%m_%d')
    print(date_tag)
    plotout = date_tag + '.png'
#    plot = Ngl.contour_map(wks, aice, res)
#    Ngl.text_ndc(wks, date_tag, 0.85, 0.84, txres)
    plt.clf()
    pyroms_toolbox.transectview(temp, -1, istart, iend, jstart, \
          jend, 'ARCTIC2', outfile=plotout, cmin=-1.8, cmax=10.25)
    nc.close()
#    Ngl.frame(wks)

#Ngl.end()
