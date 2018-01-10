# encoding: utf-8

import os
import numpy as np
import glob
import re
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF

import pyroms
import pyroms_toolbox
import _remapping

import matplotlib.pyplot as plt

import datetime

def remapping_bound_sig(varname, srcfile, wts_files, srcgrd, dst_grd, \
              rotate_sig=False, trange=None, irange=None, jrange=None, \
              dstdir='./' ,zlevel=None, dmax=0, cdepth=0, kk=0):
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
    if type(dst_grd).__name__ == 'ROMS_Grid':
        dst_grd = dst_grd
    else:
        dst_grd = pyroms.grid.get_ROMS_grid(dst_grd)

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

    sides = ['_west','_east','_north','_south']
    long = {'_west':'Western', '_east':'Eastern', \
            '_north':'Northern', '_south':'Southern'}
    dimexcl = {'_west':'xi', '_east':'xi', \
            '_north':'eta', '_south':'eta'}

    nctidx = 0
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
        if nctidx == 0:
            dstfile = dstdir + os.path.basename(srcfile[nf])[:-3] + '_' \
                   + dst_grd.name + '_bdry.nc'
            if os.path.exists(dstfile) is False:
                print('Creating destination file', dstfile)
                pyroms_toolbox.nc_create_roms_file(dstfile, dst_grd, \
                    ocean_time, lgrid=False)

            # open destination file
            nc = netCDF.Dataset(dstfile, 'a', format='NETCDF3_64BIT')

        # loop over time
        for nt in trange:

            nc.variables['ocean_time'][nctidx] = ocean_time[nt]

            # loop over variable
            for nv in range(nvar):
                print(' ')
                print('remapping', varname[nv], 'from', srcgrd.name, \
                      'to', dst_grd.name)
                print('time =', ocean_time[nt])
                Mp, Lp = dst_grd.hgrid.mask_rho.shape

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
                if nctidx == 0:
                    for sid in sides:
                       varn = varname[nv]+str(sid)
                       print('Creating variable', varn)
                       dimens = [i for i in src_var.dimensions]
                       for dim in dimens:
                           if re.match(dimexcl[sid],dim):
                               dimens.remove(dim)
                       nc.createVariable(varn, 'f8', dimens, \
                           fill_value=spval)
                       nc.variables[varn].long_name = varname[nv] + \
                            ' ' + long[sid] + ' boundary condition'
                       try:
                           nc.variables[varn].units = src_var.units
                       except:
                           print(varn+' has no units')
                       nc.variables[varn].time = src_var.time
                       nc.variables[varn].coordinates = \
                           str(dimens.reverse())
                       nc.variables[varn].field = src_var.field

                # get the right remap weights file
                for s in range(len(wts_files)):
                    if wts_files[s].__contains__(Cpos+'_to_'+Cpos+'.nc'):
                        wts_file = wts_files[s]
                        break
                    else:
                        if s == len(wts_files) - 1:
                            raise ValueError('Did not find the appropriate remap weights file')

#                print datetime.datetime.now()
                # horizontal interpolation using scrip weights
#                print 'horizontal interpolation using scrip weights'
            if not rotate_sig:
                dst_var = pyroms.remapping.remap(tmp_src_var, wts_file, \
                                                  spval=spval)

                dst_var_north = dst_var[-1, :]
                dst_var_south = dst_var[0, :]
                dst_var_east = dst_var[:, -1]
                dst_var_west = dst_var[:, 0]

                # write data in destination file
                print('write data in destination file')
                sid = '_west'
                varn = varname[nv]+str(sid)
                nc.variables[varn][nctidx] = np.squeeze(dst_var_west)

                sid = '_east'
                varn = varname[nv]+str(sid)
                nc.variables[varn][nctidx] = np.squeeze(dst_var_east)

                sid = '_north'
                varn = varname[nv]+str(sid)
                nc.variables[varn][nctidx] = np.squeeze(dst_var_north)

                sid = '_south'
                varn = varname[nv]+str(sid)
                nc.variables[varn][nctidx] = np.squeeze(dst_var_south)

            # rotate the velocity field if requested
            if rotate_sig:
                print(' ')
                print('remapping and rotating sigma from', srcgrd.name, \
                      'to', dst_grd.name)

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


                # horizontal interpolation using scrip weights
                print('horizontal interpolation using scrip weights')
                dst_11 = pyroms.remapping.remap(src_11, wts_file, \
                                                  spval=spval)
                dst_22 = pyroms.remapping.remap(src_22, wts_file, \
                                                  spval=spval)
                dst_12 = pyroms.remapping.remap(src_12, wts_file, \
                                                  spval=spval)
                Mp, Lp = dst_grd.hgrid.mask_rho.shape

                dst_11_north = dst_11[Mp-1, 0:Lp]
                dst_22_north = dst_22[Mp-1, 0:Lp]
                dst_12_north = dst_12[Mp-1, 0:Lp]

                dst_11_south = dst_11[0, 0:Lp]
                dst_22_south = dst_22[0, 0:Lp]
                dst_12_south = dst_12[0, 0:Lp]

                dst_11_east = dst_11[0:Mp, Lp-1]
                dst_22_east = dst_22[0:Mp, Lp-1]
                dst_12_east = dst_12[0:Mp, Lp-1]

                dst_11_west = dst_11[0:Mp, 0]
                dst_22_west = dst_22[0:Mp, 0]
                dst_12_west = dst_12[0:Mp, 0]

                # rotate stress tensor
                for s in range(len(wts_files)):
                    if wts_files[s].__contains__('rho_to_rho.nc'):
                        wts_file = wts_files[s]
                src_ang = srcgrd.hgrid.angle_rho[jjrange[0]:jjrange[1],iirange[0]:iirange[1]]
                src_angle = pyroms.remapping.remap(src_ang, wts_file)
                dst_angle = dst_grd.hgrid.angle_rho
                angle = dst_angle - src_angle
                cos_ang = np.cos(angle)
                sin_ang = np.sin(angle)
                Lp = cos_ang.shape[-1]
                Mp = cos_ang.shape[-2]
                print("Lp, Mp", Lp, Mp)

                if rotate_sig:
                    # North
                    for i in range(Lp):
                        Qrot = [[cos_ang[Mp-1,i], sin_ang[Mp-1,i]],
                               [-sin_ang[Mp-1,i], cos_ang[Mp-1,i]]]
                        QrotT = [[cos_ang[Mp-1,i], -sin_ang[Mp-1,i]],
                                 [sin_ang[Mp-1,i],  cos_ang[Mp-1,i]]]
                        sig = [[dst_11_north[i], dst_12_north[i]],
                               [dst_12_north[i], dst_22_north[i]]]
                        sig_rot = np.dot(np.dot(Qrot, sig), QrotT)
                        dst_11_north[i] = sig_rot[0,0]
                        dst_12_north[i] = sig_rot[0,1]
                        dst_22_north[i] = sig_rot[1,1]

                    # South
                    for i in range(Lp):
                        Qrot = [[cos_ang[0,i], sin_ang[0,i]],
                               [-sin_ang[0,i], cos_ang[0,i]]]
                        QrotT = [[cos_ang[0,i], -sin_ang[0,i]],
                                 [sin_ang[0,i],  cos_ang[0,i]]]
                        sig = [[dst_11_south[i], dst_12_south[i]],
                               [dst_12_south[i], dst_22_south[i]]]
                        sig_rot = np.dot(np.dot(Qrot, sig), QrotT)
                        dst_11_south[i] = sig_rot[0,0]
                        dst_12_south[i] = sig_rot[0,1]
                        dst_22_south[i] = sig_rot[1,1]

                    # East
                    for j in range(Mp):
                        Qrot = [[cos_ang[j,Lp-1], sin_ang[j,Lp-1]],
                               [-sin_ang[j,Lp-1], cos_ang[j,Lp-1]]]
                        QrotT = [[cos_ang[j,Lp-1], -sin_ang[j,Lp-1]],
                                 [sin_ang[j,Lp-1],  cos_ang[j,Lp-1]]]
                        sig = [[dst_11_east[j], dst_12_east[j]],
                               [dst_12_east[j], dst_22_east[j]]]
                        sig_rot = np.dot(np.dot(Qrot, sig), QrotT)
                        dst_11_east[j] = sig_rot[0,0]
                        dst_12_east[j] = sig_rot[0,1]
                        dst_22_east[j] = sig_rot[1,1]

                    # West
                    for j in range(Mp):
                        Qrot = [[cos_ang[j,0], sin_ang[j,0]],
                               [-sin_ang[j,0], cos_ang[j,0]]]
                        QrotT = [[cos_ang[j,0], -sin_ang[j,0]],
                                 [sin_ang[j,0],  cos_ang[j,0]]]
                        sig = [[dst_11_west[j], dst_12_west[j]],
                               [dst_12_west[j], dst_22_west[j]]]
                        sig_rot = np.dot(np.dot(Qrot, sig), QrotT)
                        dst_11_west[j] = sig_rot[0,0]
                        dst_12_west[j] = sig_rot[0,1]
                        dst_22_west[j] = sig_rot[1,1]


                # spval
                idx_north = np.where(dst_grd.hgrid.mask_rho[-1,:] == 0)
                idx_south = np.where(dst_grd.hgrid.mask_rho[0,:] == 0)
                idx_east = np.where(dst_grd.hgrid.mask_rho[:,-1] == 0)
                idx_west = np.where(dst_grd.hgrid.mask_rho[:,0] == 0)

                dst_11_north[idx_north[0]] = spval
                dst_22_north[idx_north[0]] = spval
                dst_12_north[idx_north[0]] = spval
                dst_11_south[idx_south[0]] = spval
                dst_22_south[idx_south[0]] = spval
                dst_12_south[idx_south[0]] = spval
                dst_11_east[idx_east[0]] = spval
                dst_22_east[idx_east[0]] = spval
                dst_12_east[idx_east[0]] = spval
                dst_11_west[idx_west[0]] = spval
                dst_22_west[idx_west[0]] = spval
                dst_12_west[idx_west[0]] = spval

                # write data in destination file
                print('write data in destination file')
                sid = '_west'
                varn = 'sig11'+str(sid)
                nc.variables[varn][nctidx] = dst_11_west
                varn = 'sig22'+str(sid)
                nc.variables[varn][nctidx] = dst_22_west
                varn = 'sig12'+str(sid)
                nc.variables[varn][nctidx] = dst_12_west

                sid = '_north'
                varn = 'sig11'+str(sid)
                nc.variables[varn][nctidx] = dst_11_north
                varn = 'sig22'+str(sid)
                nc.variables[varn][nctidx] = dst_22_north
                varn = 'sig12'+str(sid)
                nc.variables[varn][nctidx] = dst_12_north

                sid = '_east'
                varn = 'sig11'+str(sid)
                nc.variables[varn][nctidx] = dst_11_east
                varn = 'sig22'+str(sid)
                nc.variables[varn][nctidx] = dst_22_east
                varn = 'sig12'+str(sid)
                nc.variables[varn][nctidx] = dst_12_east

                sid = '_south'
                varn = 'sig11'+str(sid)
                nc.variables[varn][nctidx] = dst_11_south
                varn = 'sig22'+str(sid)
                nc.variables[varn][nctidx] = dst_22_south
                varn = 'sig12'+str(sid)
                nc.variables[varn][nctidx] = dst_12_south

            nctidx = nctidx + 1
            nc.sync()
        # close files here? how?

    # close destination file
    nc.close()

    return
