# encoding: utf-8

import os
import numpy as np
import glob
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF

import pyroms
import pyroms_toolbox
import _remapping

import matplotlib.pyplot as plt

import datetime

def remapping_tensor(varname, srcfile, wts_files, srcgrd, dstgrd, \
              rotate_sig=False, trange=None, irange=None, jrange=None, \
              dstdir='./', shapiro=False):
    '''
    A remapping function to go from a ROMS grid to another ROMS grid.
    This is for 2D tensors: internal ice stress, hard-coding for sig11, sig22,
    sig12.
    '''

    # get input and output grid
    if type(srcgrd).__name__ == 'ROMS_Grid':
        srcgrd = srcgrd
    else:
        srcgrd = pyroms.grid.get_ROMS_grid(srcgrd)
    if type(dstgrd).__name__ == 'ROMS_Grid':
        dstgrd = dstgrd
    else:
        dstgrd = pyroms.grid.get_ROMS_grid(dstgrd)

    # varname argument
    if type(varname).__name__ == 'list':
        nvar = len(varname)
    elif type(varname).__name__ == 'str':
        varname = [varname]
        nvar = len(varname)
    else:
        raise ValueError('varname must be a str or a list of str')

    # srcfile argument
    if type(srcfile).__name__ == 'list':
        nfile = len(srcfile)
    elif type(srcfile).__name__ == 'str':
        srcfile = sorted(glob.glob(srcfile))
        nfile = len(srcfile)
    else:
        raise ValueError('src_srcfile must be a str or a list of str')

    # get wts_file
    if type(wts_files).__name__ == 'str':
        wts_files = sorted(glob.glob(wts_files))
 
    # loop over the srcfile
    for nf in range(nfile):
        print('Working with file', srcfile[nf], '...')

        # get time 
        ocean_time = pyroms.utility.get_nc_var('ocean_time', srcfile[nf])
        ntime = len(ocean_time[:])

        # trange argument
        if trange is None:
            trange = list(range(ntime))

        # create destination file
        dstfile = dstdir + os.path.basename(srcfile[nf])[:-3] + '_' + dstgrd.name + '.nc'
        if os.path.exists(dstfile) is False:
            print('Creating destination file', dstfile)
            pyroms_toolbox.nc_create_roms_file(dstfile, dstgrd, ocean_time)

        # open destination file
        nc = netCDF.Dataset(dstfile, 'a', format='NETCDF3_64BIT')

        nctidx = 0
        # loop over time
        for nt in trange:

            nc.variables['ocean_time'][nctidx] = ocean_time[nt]

            # loop over variable
            for nv in range(nvar):
                print(' ')
                print('remapping', varname[nv], 'from', srcgrd.name, \
                      'to', dstgrd.name)
                print('time =', ocean_time[nt])   

                # get source data
                src_var = pyroms.utility.get_nc_var(varname[nv], srcfile[nf])

                # get spval
                try:
                    spval = src_var._FillValue
                except:
                    raise Warning('Did not find a _FillValue attribute.') 

                # irange
                if irange is None:
                    iirange = (0,src_var.shape[-1])
                else:
                    iirange = irange

                # jrange
                if jrange is None:
                    jjrange = (0,src_var.shape[-2])
                else:
                    jjrange = jrange

                # determine where on the C-grid these variable lies
                if src_var.dimensions[2].find('_rho') != -1:
                    Cpos='rho'
                else:
                    print("Sigma should be on rho points")

                print('Arakawa C-grid position is', Cpos)

                # create variable in _destination file
                if nt == trange[0]:
                    print('Creating variable', varname[nv])
                    nc.createVariable(varname[nv], 'f8', src_var.dimensions, fill_value=spval)
                    nc.variables[varname[nv]].long_name = src_var.long_name
                    try:
                        nc.variables[varname[nv]].units = src_var.units
                    except:
                        print(varname[nv]+' has no units')
                    nc.variables[varname[nv]].time = src_var.time
                    nc.variables[varname[nv]].coordinates = \
                        src_var.coordinates
                    nc.variables[varname[nv]].field = src_var.field
#                    nc.variables[varname[nv]]._FillValue = spval

                # get the right remap weights file
                for s in range(len(wts_files)):
                    if wts_files[s].__contains__(Cpos+'_to_'+Cpos+'.nc'):
                        wts_file = wts_files[s]
                        break
                    else:
                        if s == len(wts_files) - 1:
                            raise ValueError('Did not find the appropriate remap weights file')


                # write data in destination file
#                print 'write data in destination file'
#                nc.variables[varname[nv]][nctidx] = dst_var

            # rotate the velocity field if requested
#            print datetime.datetime.now()
            print(' ') 
            print('remapping and rotating sigma from', srcgrd.name, \
                  'to', dstgrd.name)

            # get source data
            src_11 = pyroms.utility.get_nc_var(varname[0], srcfile[nf])
            # get spval
            try:
                spval = src_11._FillValue
            except:
                raise Warning('Did not find a _FillValue attribute.') 

            src_11 = src_11[nt,jjrange[0]:jjrange[1],iirange[0]:iirange[1]]

            src_22 = pyroms.utility.get_nc_var(varname[1], srcfile[nf])
            src_22 = src_22[nt,jjrange[0]:jjrange[1],iirange[0]:iirange[1]]

            src_12 = pyroms.utility.get_nc_var(varname[2], srcfile[nf])
            src_12 = src_12[nt,jjrange[0]:jjrange[1],iirange[0]:iirange[1]]

            print("before", src_11[-1,30], src_12[-1,30], src_22[-1,30])
            if shapiro:
                src_11 = pyroms_toolbox.shapiro_filter.shapiro2(src_11,2)
                src_22 = pyroms_toolbox.shapiro_filter.shapiro2(src_22,2)
                src_12 = pyroms_toolbox.shapiro_filter.shapiro2(src_12,2)
            print("after", src_11[-1,30], src_12[-1,30], src_22[-1,30])

            # horizontal interpolation using scrip weights
            print('horizontal interpolation using scrip weights')
            dst_11 = pyroms.remapping.remap(src_11, wts_file, \
                                              spval=spval)
            dst_22 = pyroms.remapping.remap(src_22, wts_file, \
                                              spval=spval)
            dst_12 = pyroms.remapping.remap(src_12, wts_file, \
                                              spval=spval)
            print("after remapping", dst_11[-1,30], dst_12[-1,30], dst_22[-1,30])

            if rotate_sig is True:
                # rotate stress tensor
                src_ang = srcgrd.hgrid.angle_rho[jjrange[0]:jjrange[1],iirange[0]:iirange[1]]
                src_angle = pyroms.remapping.remap(src_ang, wts_file)
                dst_angle = dstgrd.hgrid.angle_rho
                angle = dst_angle - src_angle
                cos_ang = np.cos(angle)
                sin_ang = np.sin(angle)
                Lp = cos_ang.shape[-1]
                Mp = cos_ang.shape[-2]
                print("Lp, Mp", Lp, Mp)

                for j in range(Mp):
                    for i in range(Lp):
                        Qrot = [[cos_ang[j,i], sin_ang[j,i]],
                               [-sin_ang[j,i], cos_ang[j,i]]]
                        QrotT = [[cos_ang[j,i], -sin_ang[j,i]],
                                 [sin_ang[j,i],  cos_ang[j,i]]]
#                        Qrot = [[cos_ang[j,i], -sin_ang[j,i]],
#                               [sin_ang[j,i], cos_ang[j,i]]]
#                        QrotT = [[cos_ang[j,i], sin_ang[j,i]],
#                                 [-sin_ang[j,i],  cos_ang[j,i]]]
                        sig = [[dst_11[j,i], dst_12[j,i]],
                               [dst_12[j,i], dst_22[j,i]]]
                        sig_rot = np.dot(np.dot(Qrot, sig), QrotT)
                        dst_11[j,i] = sig_rot[0,0]
                        dst_12[j,i] = sig_rot[0,1]
                        dst_22[j,i] = sig_rot[1,1]
                print("after rotating", dst_11[-1,30], dst_12[-1,30], dst_22[-1,30])


            # spval
            idx = np.where(dstgrd.hgrid.mask_rho == 0)
            dst_11[idx[0], idx[1]] = spval
            dst_12[idx[0], idx[1]] = spval
            dst_22[idx[0], idx[1]] = spval

            # write data in destination file
            print('write data in destination file')
            nc.variables['sig11'][nctidx] = dst_11
            nc.variables['sig12'][nctidx] = dst_12
            nc.variables['sig22'][nctidx] = dst_22

        nctidx = nctidx + 1
        nc.sync()
 
    # close destination file
    nc.close()

    return
