import numpy as np
import sys
import netCDF4
import pyroms
import pyroms_toolbox
import Ngl


def plot_Robinson_pyngl(var,lon,lat,wks_name='plot', clim=None, cspace=None, cmap=None):

    #
    #  Select a colormap and open a workstation.
    #
    rlist            = Ngl.Resources()

    wks_type = "png"
    wks = Ngl.open_wks(wks_type,wks_name,rlist)

    if cmap is None:
        mycmap = 'BlAqGrYeOrReVi200'
        Ngl.define_colormap(wks,mycmap)
    else:
        try:
            mycmap = '/Users/frederic/python/cmap/' + cmap    
            mycmap = np.loadtxt(mycmap)
        except:
            mycmap = cmap
        try:
            Ngl.define_colormap(wks,mycmap)
        except:
            raise Warning('Unknown colormap')

    #
    #  The next set of resources will apply to the contour plot and the labelbar.
    #
    resources = Ngl.Resources()

    resources.sfXArray            = lon
    resources.sfYArray            = lat

    resources.gsnMaximize         = True                # use full page


    resources.cnFillOn            = True
    resources.cnFillMode          = "RasterFill"
    resources.cnMaxLevelCount     = 255 
    resources.cnLinesOn           = False
    resources.cnLineLabelsOn      = False
    if clim is not None:
        resources.cnLevelSelectionMode = "ManualLevels"          # set manual contour levels
        resources.cnMinLevelValF       = clim[0]                 # set min contour level
        resources.cnMaxLevelValF       = clim[1]                 # set max contour level
        if cspace is None:
            LevelSpacing = (clim[1] - clim[0]) / 100.
            resources.cnLevelSpacingF      =  LevelSpacing       # set contour spacing
        else:
            resources.cnLevelSpacingF      =  cspace             # set contour spacing

    resources.lbOrientation      = "Horizontal"     # Default is vertical.
    resources.lbBoxLinesOn       = False
    resources.lbLabelFontHeightF = 0.01            # label font height
    resources.lbBoxMinorExtentF  = 0.15

    #
    # The contour plot is not very interesting, so don't draw it.
    # 
    resources.nglDraw  = False
    resources.nglFrame = False

    contour = Ngl.contour(wks,var,resources)

    #
    # Retrieve the actual lat/lon end points of the scalar array so
    # we know where to overlay on map.
    #
    xs = Ngl.get_float(contour.sffield,"sfXCActualStartF")
    xe = Ngl.get_float(contour.sffield,"sfXCActualEndF")
    ys = Ngl.get_float(contour.sffield,"sfYCActualStartF")
    ye = Ngl.get_float(contour.sffield,"sfYCActualEndF")

    resources.nglDraw           = True        # Turn these resources back on.
    resources.nglFrame          = True

    resources.mpProjection      = "Robinson"
    resources.mpCenterLonF = 270
    resources.mpCenterLatF =   0

    resources.mpGeophysicalLineThicknessF = 2

    map = Ngl.contour_map(wks,var,resources)

    Ngl.end()
