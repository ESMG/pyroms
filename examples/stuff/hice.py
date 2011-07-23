import numpy, os
import netCDF4
import Ngl

file1 = netCDF4.Dataset("grid_zhang.nc", "r")
file2 = netCDF4.Dataset("ice_1985.nc", "r")

clat = file1.variables["clat"][:]
clon = file1.variables["clon"][:]

hice = file2.variables["hice"][12,:,:]

#********************************************
wks  = Ngl.open_wks("ncgm", "uice")            # open workstation
cmap = numpy.array([[1., 1., 1.], [0., 0., 0.], [.85,0., 0.], [1., 0., 0.], \
          [1.,.13,.13], [1.,.25,.25], [1.,.38,.38], [1.,.48,.48], \
          [1.,.58,.58], [1.,.67,.67], [1.,.76,.76], [1.,.85,.85], \
          [.9, .9, .9], [.75,1., 1.], [.5, 1., 1.], [.25,1., 1.], \
          [0., .85,1.], [0. ,.76,1.], [0., .67,1.], [0. ,.58,1.], \
          [0., .38,1.], [0., 0., 1.], [0., 0.,.85]], 'f')

rlist = Ngl.Resources()
rlist.wkColorMap = cmap
Ngl.set_values(wks, rlist)

#;  ncolors = dimsizes(cmap(:,0)) - 2
#;  gsn_define_colormap(wks, cmap)
#
#;  zmin = min(uice)
#;  zmax = max(uice)
#;  cnLevels = new((/ncolors-1/), "double")
#;  do i=0,ncolors-2
#;    cnLevels(i) = zmin + (zmax-zmin)*(i+1)/ncolors
#;  end do
#;  spcng = (zmax-zmin)/ncolors
#
#; Names is used when drawing the labelbar
#;  names = new(dimsizes(cnLevels), string)
#;  do i=0,dimsizes(cnLevels)-1
#;    names(i) = sprintf("%5.1f", doubletofloat(cnLevels(i)))
#;  end do

#i = NhlNewColor(wks,0.8,0.8,0.8)                #; add gray to colormap 
                                          
res                        = Ngl.Resources()    #; plot mods desired
res.sfXArray               = clon
res.sfYArray               = clat

res.cnFillOn               = True               #; color fill  
res.cnLinesOn              = True              #; no contour lines
res.cnLineLabelsOn         = False             # ; no contour labels
res.cnFillDrawOrder        = "PreDraw"         # ; put continents on top

#res.gsnSpreadColors        = True              # ; use total colormap
#res.gsnSpreadColorEnd      = -3
# res@cnInfoLabelOn          = False              ; no contour info label
# res@cnLevelSelectionMode = "ManualLevels"
# res@cnMinLevelValF       = cnLevels(0)
# res@cnMaxLevelValF       = cnLevels(ncolors-2)
# res@cnLevelSpacingF      = spcng

# Chukchi plot parameters
res.mpProjection        = "Stereographic"
res.mpCenterLatF = 90
res.mpCenterLonF = 205
 
res.mpLimitMode         = "Corners"
res.mpLeftCornerLatF    = 60.0
res.mpLeftCornerLonF    = 180.0
res.mpRightCornerLatF   = 74.
res.mpRightCornerLonF   = 250.0

res.mpLandFillColor        = "Tan1"
res.mpOceanFillColor       = "SkyBlue"
res.mpInlandWaterFillColor = "SkyBlue"

plot = Ngl.contour_map(wks, hice, res)

res.cnFillOn                      = True
res.cnLinesOn                     = False
res.cnLineLabelsOn                = False
res.cnFillMode                    = "CellFill"
res.cnCellFillEdgeColor           = "Black"
res.cnMonoFillColor               = True
res.cnFillColor                   = "Transparent"
res.cnCellFillMissingValEdgeColor = "Red"

plot = Ngl.contour_map(wks, hice, res)
#plot = Ngl.contour_map(wks, hice[:-2,:-2], res)

Ngl.end()
