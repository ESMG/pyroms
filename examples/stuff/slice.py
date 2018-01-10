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
#lst = commands.getoutput('ls averages/*74??.nc')
lst = subprocess.getoutput('ls months/*.nc')
lst = lst.split()
lst_file = lst_file + lst

grd = pyroms.grid.get_ROMS_grid('ARCTIC2')

#clat = grd.hgrid.lat_rho
#clon = grd.hgrid.lon_rho

#********************************************
#wks  = Ngl.open_wks("ncgm", "aice")            # open workstation
#cmap = ["white", "black", \
#	  "(/0., 0.,.85/)", "(/ 0., 0.,1./)", "(/.30,.30,1./)", \
#	  "(/.42,.42,1./)", "(/.55,.55,1./)", \
#	  "(/.64,.64,1./)", "(/.72,.72,1./)", "(/.80,.80,1./)", \
#	  "(/.88,.88,1./)", "(/.9, .9, .9/)", "burlywood"]

#rlist = Ngl.Resources()
#rlist.wkColorMap = cmap
#Ngl.set_values(wks, rlist)

#ncolors = len(cmap) - 3
#print ncolors
#;  gsn_define_colormap(wks, cmap)
#
#zmin = 0.0
#zmax = 1.0
#cnLevels = np.zeros(ncolors-1)
#for i in range(ncolors-1):
#    cnLevels[i] = zmin + (zmax-zmin)*(i+1)/ncolors
#spcng = (zmax-zmin)/ncolors

#
#; Names is used when drawing the labelbar
#;  names = new(dimsizes(cnLevels), string)
#;  do i=0,dimsizes(cnLevels)-1
#;    names(i) = sprintf("%5.1f", doubletofloat(cnLevels(i)))
#;  end do

#i = NhlNewColor(wks,0.8,0.8,0.8)                #; add gray to colormap 
                                          
#res                        = Ngl.Resources()    #; plot mods desired
#res.nglMaximize            = True
#res.sfXArray               = clon
#res.sfYArray               = clat
#
#res.cnFillOn               = True               #; color fill  
#res.cnLinesOn              = True              #; no contour lines
#res.cnLineLabelsOn         = False             # ; no contour labels
#res.cnFillDrawOrder        = "PreDraw"         # ; put continents on top
#
##res.gsnSpreadColors        = True              # ; use total colormap
##res.gsnSpreadColorEnd      = -3
## res.cnInfoLabelOn          = False              ; no contour info label
#res.cnLevelSelectionMode = "ManualLevels"
#res.cnMinLevelValF       = cnLevels[0]
#res.cnMaxLevelValF       = cnLevels[ncolors-2]
#res.cnLevelSpacingF      = spcng

## Chukchi plot parameters
#res.mpProjection        = "Stereographic"
#res.mpCenterLatF = 90
#res.mpCenterLonF = 205 
#
#res.mpLimitMode         = "Corners"
#res.mpLeftCornerLatF    = 60.0
#res.mpLeftCornerLonF    = 180.0
#res.mpRightCornerLatF   = 74. 
#res.mpRightCornerLonF   = 250.0
#res.mpDataBaseVersion = "MediumRes"
#
#res.mpFillOn            = True
#res.mpFillColors        = ["background","transparent","burlywood","transparent"]
#
#res.nglFrame     = False                    # we can attach some text.
#res.tiMainString = "ROMS Arctic Simulation"
#res.tiMainOffsetYF    = 0.04
#res.tiMainFontHeightF = 0.02
#res.nglSpreadColorEnd   = -2
#res.lbOrientation            = "Horizontal"
#res.lbTitleString    = "Sea Ice Concentration"
#res.lbTitleFontHeightF        = 0.012
#res.lbLabelFontHeightF        = 0.015
#res.pmLabelBarOrthogonalPosF = +0.02
#res.pmLabelBarHeightF        = 0.1
#res.pmLabelBarWidthF         = 0.6
#
#txres             = Ngl.Resources()          # Text resources desired
#txres.txFontHeightF = 0.015

#txres.txFontColor = "OrangeRed"

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
