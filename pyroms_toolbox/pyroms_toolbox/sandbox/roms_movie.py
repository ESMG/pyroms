# encoding: utf-8

import os
import sys
import subprocess
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import octant 



def make_movie(filelst, varname, cmin, cmax, view, lev=0, istart=None, iend=None, grd=None, \
                proj='merc', imode='on', title='', clean=False):
    """
       make a movie using a 2D horizontal slice of variable varname from nc files in filelst
    """


    # get variable
    data = netCDF4.MFDataset(filelst)
    var = data.variables[varname]

    # get grid
    if grd is None:
        grd = octant.roms.roms_grid.get_roms_grd(filelst[0])
    else:
        grd = grd

    # determine where on the C-grid these variable lies
    if var.dimensions[2].find('_rho') != -1:
        Cpos='rho'
        mask = grd.mask_rho[:]

    if var.dimensions[2].find('_u') != -1:
        Cpos='u'
        mask = grd.mask_u[:]

    if var.dimensions[2].find('_v') != -1:
        Cpos='v'
        mask = grd.mask_v[:]

    # get time 
    time = data.variables['ocean_time'][:]


    if istart is not None:
        istart = istart
    else:
        istart = 0

    if iend is not None:
        iend = iend
    else:
        iend = time.shape[0]


    if cmin is None:
        cmin = var.min()
    else:
        cmin = float(cmin)

    if cmax is None:
        cmax = var.max()
    else:
        cmax = float(cmax)

    if imode is 'off':
        print('Turn interactive mode off')
        plt.ioff()

    for tindex in range(istart, iend, 1):
        if view is 'sview':
            octant.roms.roms_plot.varsview(var, grd, tindex, lev, Cpos, cmin=cmin, cmax=cmax, proj=proj, \
                                   title=title, outfile='plot.png')
            
        elif view is 'zview':
            octant.roms.roms_plot.varzview(var, grd, tindex, lev, Cpos, cmin=cmin, cmax=cmax, proj=proj, \
                                   title=title, outfile='plot.png')

        elif view is 'view2D':
            octant.roms.roms_plot.varview2D(var, grd, tindex, Cpos, cmin=cmin, cmax=cmax, proj=proj, \
                                    title=title, outfile='plot.png')

        else:
            print('Option not available. view must be set to sview, zview or view2D') 


        outfile = str('%05d' % tindex) + '.png'
        os.rename('plot.png',outfile)        


    command = ('mencoder',
               'mf://0*.png',
               '-mf',
               'type=png:w=800:h=600:fps=25',
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec=mpeg4',
               '-oac',
               'copy',
               '-o',
               'output.avi')

    print("\n\nabout to execute:\n%s\n\n" % ' '.join(command))
    subprocess.check_call(command)

    print("\n\n The movie was written to 'output.avi'")

    if imode is 'off':
        print('Turn interactive mode on again')
        plt.ion()

    if clean is True:
        for tindex in range(istart, iend, 1):
            os.remove(str('%04d' % tindex) + '.png')


    return



def make_big_movie(filelst, varname, cmin, cmax, Cpos, view, lev=0, grd=None, \
                    proj='merc', imode='on', title='', clean=False):
    """
       make a movie using a 2D horizontal slice of variable varname from nc files in filelst
    """

    # get grid
    if grd is None:
        grd = octant.roms.roms_grid.get_roms_grd(filelst[0])
    else:
        grd = grd

    # get mask
    if Cpos is 'rho':
        mask = grd.mask_rho[:]

    if Cpos is 'u':
        mask = grd.mask_u[:]

    if Cpos is 'v':
        mask = grd.mask_v[:]

    nfile = len(filelst)

    if imode is 'off':
        print('Turn interactive mode off')
        plt.ioff()

    counter = 0

    for ifile in range(nfile):

        # get variable
        data = netCDF4.Dataset(filelst[ifile], 'r')
        var = data.variables[varname]

        # get time 
        time = data.variables['ocean_time'][:]

        for tindex in range(time.shape[0]):
            if view is 'sview':
                octant.roms.roms_plot.varsview(var, grd, tindex, lev, Cpos, cmin=cmin, cmax=cmax, proj=proj, \
                                       title=title, outfile='plot.png')
            
            elif view is 'zview':
                octant.roms.roms_plot.varzview(var, grd, tindex, lev, Cpos, cmin=cmin, cmax=cmax, proj=proj, \
                                       title=title, outfile='plot.png')

            elif view is 'view2D':
                octant.roms.roms_plot.varview2D(var, grd, tindex, Cpos, cmin=cmin, cmax=cmax, proj=proj, \
                                        title=title, outfile='plot.png')

            else:
                print('Option not available. view must be set to sview, zview or view2D') 

            Tindex = counter + tindex

            outfile = str('%05d' % Tindex) + '.png'
            os.rename('plot.png',outfile)        

            
        counter = counter + time.shape[0]


    command = ('mencoder',
               'mf://0*.png',
               '-mf',
               'type=png:w=800:h=600:fps=25',
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec=mpeg4',
               '-oac',
               'copy',
               '-o',
               'output.avi')

    if imode is 'off':
        print('Turn interactive mode on again')
        plt.ion()

    print("\n\nabout to execute:\n%s\n\n" % ' '.join(command))
    subprocess.check_call(command)
 
    print("\n\n The movie was written to 'output.avi'")

    if clean is True:
        for tindex in range(counter):
            os.remove(str('%05d' % tindex) + '.png')


    return
