import numpy as np
import netCDF4
import sys
import subprocess
import Ngl
import pyroms
from pyroms_toolbox import jday2date

#year = int(sys.argv[1])
#lst_year = [year]

#for year in lst_year:
#    year = np.str(year)
lst = subprocess.check_output('ls months/*.nc', shell=True)
lst_file = lst.split()

grd = pyroms.grid.get_ROMS_grid('ARCTIC2')

clat = grd.hgrid.lat_rho
clon = grd.hgrid.lon_rho

#********************************************
wks  = Ngl.open_wks("ncgm", "aice")            # open workstation
#cmap = ["white", "black", \
#	  "(/.9, .9, .9/)", "(/1.,.88,.88/)", "(/1.,.80,.80/)", \
#	  "(/1.,.72,.72/)", "(/1.,.64,.64/)", "(/1.,.55,.55/)", \
#	  "(/1.,.42,.42/)", "(/1.,.30,.30/)", \
#	  "(/1., 0., 0./)", "(/.85,0., 0./)"]
cmap = ["white", "black", \
          "(/0., 0.,.85/)", "(/ 0., 0.,1./)", "(/.30,.30,1./)", \
          "(/.42,.42,1./)", "(/.55,.55,1./)", \
          "(/.64,.64,1./)", "(/.72,.72,1./)", "(/.80,.80,1./)", \
          "(/.88,.88,1./)", "(/.9, .9, .9/)", "burlywood"]

rlist = Ngl.Resources()
rlist.wkColorMap = cmap
Ngl.set_values(wks, rlist)

ncolors = len(cmap) - 3
#Levels             = ( 0.5, 1., 2., 3., 5., 8., 10., 12., 15. )
zmin = 0.0
zmax = 1.0
cnLevels = np.zeros(ncolors-1)
for i in range(ncolors-1):
    cnLevels[i] = zmin + (zmax-zmin)*(i+1)/ncolors
spcng = (zmax-zmin)/ncolors

#
#i = NhlNewColor(wks,0.8,0.8,0.8)                #; add gray to colormap 
                                          
res                        = Ngl.Resources()    #; plot mods desired
res.nglMaximize            = True
res.sfXArray               = clon
res.sfYArray               = clat

res.cnFillOn               = True               #; color fill  
res.cnLinesOn              = True              #; no contour lines
res.cnLineLabelsOn         = False             # ; no contour labels
res.cnFillDrawOrder        = "PreDraw"         # ; put continents on top

#res.gsnSpreadColors        = True              # ; use total colormap
#res.gsnSpreadColorEnd      = -3
# res.cnInfoLabelOn          = False              ; no contour info label
#res.cnLevelSelectionMode = "ExplicitLevels"
#res.cnLevels             = Levels
res.cnLevelSelectionMode = "ManualLevels"
res.cnMinLevelValF       = cnLevels[0]
res.cnMaxLevelValF       = cnLevels[ncolors-2]
res.cnLevelSpacingF      = spcng

# Bering plot parameters
#res.mpProjection        = "LambertConformal"
#res.mpLambertParallel1F = 40
#res.mpLambertParallel2F = 60
#res.mpLambertMeridianF  = 180
 
#res.mpLeftCornerLatF    = 50.0
#res.mpLeftCornerLonF    = 165.0
#res.mpRightCornerLatF   = 65.
#res.mpRightCornerLonF   = 215.0


# Arctic2
res.mpProjection        = "Stereographic"
res.mpCenterLatF = 90
res.mpCenterLonF = -180

res.mpLimitMode         = "Corners"
res.mpLeftCornerLatF    = 40.0
res.mpLeftCornerLonF    = -210.0
res.mpRightCornerLatF   = 50.
res.mpRightCornerLonF   = -50.0

res.mpDataBaseVersion = "MediumRes"

res.mpFillOn            = True
res.mpFillColors        = ["background","transparent","burlywood","transparent"]

res.nglFrame     = False                    # we can attach some text.
res.tiMainString = "ROMS Arctic Simulation"
res.tiMainOffsetYF    = 0.04
res.tiMainFontHeightF = 0.02
res.nglSpreadColorEnd   = -2
res.lbTitleString    = "Sea Ice Concentration"
res.lbTitleFontHeightF        = 0.012
res.lbLabelFontHeightF        = 0.015
res.pmLabelBarOrthogonalPosF = +0.02
#res.lbOrientation            = "Horizontal"
#res.pmLabelBarHeightF        = 0.1
#res.pmLabelBarWidthF         = 0.6
res.pmLabelBarHeightF        = 0.6
res.pmLabelBarWidthF         = 0.1


txres             = Ngl.Resources()          # Text resources desired
txres.txFontHeightF = 0.015

for file in lst_file:
    print("Plotting "+file)
    nc = netCDF4.Dataset(file, "r")
    aice = nc.variables["aice"][0,:,:]
    time = nc.variables["ocean_time"][0]
    myday = jday2date(time/86400.)
    date_tag = myday.strftime('%d %B %Y')
    print(date_tag)
    plot = Ngl.contour_map(wks, aice, res)
    Ngl.text_ndc(wks, date_tag, 0.85, 0.94, txres)
    nc.close()
    Ngl.frame(wks)

Ngl.end()
